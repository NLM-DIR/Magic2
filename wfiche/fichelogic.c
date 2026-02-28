#include "../wac/ac.h"
#include "../wfiche/biolog.h"
#include "../wfiche/gtitle.h"
#include "table.h"
#include "pick.h"
#include "query.h"
#define vONLINEMAXURLSIZE 2048
#define vSTRMAXNAME 2048 

#define _00    "\0\0" 
#define _0     "\0" 

#define TABLEBORDER "2"

#define HTML5 1

static void ficheMRNAExpressionProfileParagraphContent (vTXT blkp, GMP *gmp, BOOL isMrna, BOOL isBold) ;

char genomeRelease[1024] ; 
static char *dbNam="" ; 
static char *ficheMarkupContent[]= {
	/* gone _sMarkup */"",					"",
	/* gone _PN */  "<b>",					"</b><br/>\n",  /* part names markup */
	/* gone _PN2 */  "<font color=#00007f><b>",	       	"</b></font><br/>\n",  /* part names markup */
	/* gone _PB */  "<hr></hr>",					"<hr></hr>",               /* part body markup */
	/* gone _SN */   "\n\n<font color=#00007f><b>",		"</b></font><br/>\n",  /* section name markup */
	/* gone _SN2 */  "\n\n<font color=#00007f><b>",		"</b></font><br/>\n",  /* section name markup */
	/* gone _SB */  "",						"<br/><br/>",    /* section body markup */
	/* _TB */  "<table width=\"98%%\" border="TABLEBORDER">","</table>\n",          /* table  markup */
	/* _TR */  "<tr>",					"</tr>\n",            /* table row markup */
	/* _TD */  "<td VALIGN=TOP>",				"</td> ",      /* table cell markup */
	/* gone _SE */  "<font face='courier' color=#707070><pre>",	"</pre></font>",       /* sequence markup */
	/* gone _BR */  "\n<BR/>",					"", /* sequence markup */
	/* gone _LI */  "<LI>    ",					"", 
	/* gone _UL */  "<UL>",					"</UL>\n", 
    "",""} ; 

static int ficheOperonDistance=450 ; 

/* ATTENTION: mail form Donna
If you have applications that point to LocusLink, e.g. 
https://www.ncbi.nlm.nih.gov/LocusLink/LocRpt.cgi?l=6606

please be ready to replace it with the link to Gene:
https://www.ncbi.nlm.nih.gov/corehtml/query.fcgi?db=gene&cmd=Retrieve&dopt=Graphics&list_uids=6606

*/
/*
  struct tagSPECIESINFO {  char * speciesName,* speciesGenomicName;	
   int		long5P,long3P;}
*/
SPECIESINFO SpcI[SPECIESCNT] =
{ 
	/* Human */
	{"Homo sapiens"_0"human"_00
	 , "the human genome"
	 , 1094, 3395 /* measured 2007_09_16 tablemaker -f 5p3pUTRlength */
	},
	/* Worm */
	{"Caenorhabditis elegans"_00
	 , "C.elegans"
	 , 650, 695  /* measured 2007_09_16 tablemaker -f 5p3pUTRlength */
	} ,
	/* ARA */
	{"Arabidopsis thaliana"_00
	 , "A.thaliana"
	 , 520, 701  /* measured 2007_10_15 tablemaker -f 5p3pUTRlength */
	} ,
        /* DROSO */
	{"Drosophila melanogaster"_0"oikopleura dioica"_00
	 , "D.melanogaster"
	 , 600, 700 /* guessed */
	} ,
        /* COLI */
	{"Escherichia coli"_00
	 , "E.doli"
	 , 100, 100 /* guessed */
	} ,
	/* MOUSE */
	{"Mus musculus"_00
	 , "M.musculus"
	 , 832, 3400   /* measured 2007_09_16 tablemaker -f 5p3pUTRlength */
	} ,
	/* RAT */
	{"Rattus norvegicus"_00
	 , "R.norvegicus"
	 , 832, 3400   /* not measured 2007_09_16 tablemaker -f 5p3pUTRlength */
	} 
} ; 

typedef struct { int icl, isNm, isTiling ; KEY est, clone, tissue, mrna, hinv_lib ; int ali, nerr ; AC_OBJ oMrna ; } CLALI ;

static BOOL ficheNewGeneIntronsStatement (vTXT blkp, GMP *gmp, BOOL decorate, BOOL showTitle) ;
static int ficheTGOperonParagraphContent (vTXT blkp, GMP *gmp, int operonDistance, int isTitles, int isStrands, int previous, BOOL details) ;
static int ficheNewGeneMappingSentence (vTXT blkp, GMP *gmp, int prefix) ;
static void ficheGmpYkImageStatement (vTXT blkp, GMP *gmp) ;
static void ficheGmpPatternStatement (vTXT blkp, GMP *gmp) ;
static void ficheGmpExpressionPatternStatement (vTXT blkp, GMP *gmp) ;
static int ficheNewGeneAliasStatement (vTXT blkp, GMP *gmp) ;
static int ficheNewGeneComplexLocusStatement (vTXT blkp, GMP *gmp) ;
static void ficheNewGeneCoversValue (vTXT blkp, GMP *gmp, int type) ;
static int ficheNewGeneExpressionTissueStatement (vTXT blkp, GMP *gmp, BOOL decorate) ;
static int ficheNewGeneAltVariantStatement (vTXT blkp, GMP *gmp, BOOL decorate) ;
static int ficheNewGeneProteinStatement (vTXT blkp, GMP *gmp, BOOL decorate) ;
static int ficheNewGenePfamPsortStatement (vTXT blkp, GMP *gmp, BOOL isCloud, BOOL fromSummary) ;
static int ficheNewGeneNonGoodVariantStatement (vTXT blkp, GMP *gmp) ;
static int ficheNewGeneNmdStatement (vTXT blkp, GMP *gmp, BOOL decorate) ;
static int ficheNewGeneKozakStatement (vTXT blkp, GMP *gmp) ;
static int ficheNewGeneUorfStatement (vTXT blkp, GMP *gmp) ;
static int ficheNewGenePhosphositeStatement (vTXT blkp, GMP *gmp) ;
static int ficheNewGeneLocus_PhenotypeStatement (vTXT blkp, GMP *gmp) ;
static int ficheNewGeneFunctionStatement (vTXT blkp, GMP *gmp) ;
static int ficheNewGeneAltFeatureStatement (vTXT blkp, GMP *gmp, BOOL decorate, BOOL fromSummary) ;
static void markupLinkPubmed (vTXT blkp, GMP *gmp, const char *id, const char *txt) ;
static void ficheNewGenePleaseQuote  (vTXT blkp, GMP *gmp) ;
static int ficheNewGeneAnnotationOfVariantsTable (vTXT blkp, GMP *gmp, BOOL justSequences) ;
static void ficheNewGeneIntronsParagraph (vTXT blkp, GMP *gmp) ;
static void ficheNewGenePolyAParagraph (vTXT blkp, GMP *gmp) ;
static void ficheNewTissuesParagraph (vTXT blkp, GMP *gmp, Array tissues
				      , const char *title, const char *subtitle) ;
static void ficheNewTissuesSentence (vTXT blkp, GMP *gmp, Array tissues, const char *subtitle) ;
static void ficheNewMrna5PrimeParagraph (vTXT blkp, GMP *gmp) ;
static void ficheNewMrna3PrimeParagraph (vTXT blkp, GMP *gmp) ;
static void ficheNewMrnaPeptideParagraph (vTXT blkp, GMP *gmp) ;
static void ficheNewMrnaDnaParagraph (vTXT blkp, GMP *gmp) ;
void ficheNewMrnaDecoratedPeptide (vTXT blkp, GMP *gmp, BOOL hasStop) ;
void ficheNewMrnaDnaDecoratedSequence (vTXT blkp, GMP *gmp, int x1, int x2) ;
static void ficheNewGeneAceKogSubSection (vTXT blkp, GMP *gmp) ;
static AC_KEYSET ficheNewGeneDiseaseKs (AC_OBJ gene, AC_HANDLE h) ;

/* -===================================- /
* -=  The ID of Creature              =- /
* -===================================- */	

int ficheDBCreatureID (AC_DB db0)
{
  AC_ITER lCloneMain = 0 ;
  AC_OBJ oClo ; 
  const char *tmpNom ; 
  char *cp ;
  int  iSpc=HUMAN ; 
  AC_DB db = db0 ;
    
  if (!db)
    db =  ac_open_db ("local", 0) ;
  lCloneMain = ac_query_iter (db, TRUE, "find clone strategy", 0, 0) ; 
  if ((oClo = ac_next_obj (lCloneMain)))
    {
      tmpNom = ac_tag_text  (oClo, "Species", "toto") ; 
      
      if  (!strcasecmp  (tmpNom, "human")) iSpc=HUMAN ; 
      else if  (strstr  (tmpNom, "Mouse")) iSpc=MOUSE ; 
      else if  (strstr  (tmpNom, "Rat")) iSpc=RAT ; 
      else if  (!strcasecmp  (tmpNom, "coli")) iSpc=COLI ; 
      else if  (!strcasecmp  (tmpNom, "Worm")) iSpc=WORM ; 
      else if  (strstr  (tmpNom, "Droso")) iSpc=DROSO ; 
      else if  (strstr  (tmpNom, "Arabidopsis")) iSpc=ARA ;
      
      tmpNom = ac_tag_text  (oClo, "MainTitle", "") ; 
      if (!tmpNom)
	tmpNom = ac_tag_text  (oClo, "Genome_release", 0) ;
      if  (dbNam && *dbNam) ac_free (dbNam) ; 
      if  (tmpNom && tmpNom[0]) dbNam = strnew (tmpNom, 0) ; else dbNam="" ; 
      if  ( (cp = strstr (dbNam, ", with"))) *cp = 0 ;
      ac_free  (oClo) ;
    } 
  ac_free (lCloneMain) ;

  if (!db0) ac_db_close (db) ; 
  return iSpc ; 
}

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  Service non-fichelogic functions first 
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


/* -===================================- /
* -=  Print arbitrable Table          =- /
* -===================================- */	

#define	ooo  "\\\\\\"
#define	TbbDim		50
#define	TBB(_v_r,_v_c)		 array (Tbb,  (1+  (_v_r))*  (TbbDim)+  ( (_v_c)+1), int)

static BOOL isVoyelle (char c1)
{
  int c = ace_lower (c1) ;
  switch (c)
    {
    case 'a':    case 'e':    case 'i':    case 'o':    case 'u':
      return TRUE ;
    }
  return FALSE ;
}


static char *isOne (int nn)
{
  static char buf[12] ;

  if (nn == 0) return "zero" ;
  if (nn == 1) return "one" ;
  sprintf (buf, "%d", nn) ;
  return buf ;
}

/******************************************************************************/
/******************************************************************************/
/*** utilities (former MACROS ***/

static void markupStart (vTXT blkp, GMP *gmp, int type)
{ 
  if (gmp->markup)
    {
      vtxtPrintf (blkp, ficheMarkupContent [2 * type]) ;
    }
}

static void markupEnd (vTXT blkp, GMP *gmp, int type)
{
  if (gmp->markup)
    {
      vtxtPrintf (blkp, ficheMarkupContent [2 * type + 1]) ;
    }
}

static void markupText (vTXT blkp, GMP *gmp, int type, char *txt)
{
 if (gmp->markup)
    {
      markupStart (blkp, gmp, type) ;
      if (txt && *txt) vtxtPrint (blkp, txt) ;
      markupEnd (blkp, gmp, type) ;
    }
}

/*******************************************************/
/* 2025_03 Tranfer the printing of the table to ac_table_display */
static int fichePrintSquareTable (GMP *gmp, char style, Array Tbb, vTXT blkp, vTXT bfr, int doEmpty, int colColor, int sR, int eR, int colNum, ... )
{
  AC_HANDLE h = ac_new_handle () ;
  va_list marker ;
  vtxtMarkup (blkp) ;
  AC_TABLE table = 0 ;
  Array titles = arrayHandleCreate (16, const char *, h) ;
  int iC, nCol = -1 ;

  iC = colNum ; 
  va_start  (marker, colNum) ; 
  while (iC != -1)
    {
      nCol++ ;
      /* clean up the table, replace ooo cellends by zero */
      char *txt = vtxtPtr  (bfr) + TBB  (0, iC) ;
      char *ptr = strstr  (txt, ooo) ;
      if (ptr) *ptr=0 ;
      /* insert in table */
      if  (*txt)
	array (titles, nCol, const char *) = txt ;
      
      iC= va_arg  (marker, int) ; 
    }
  va_end  ( marker ) ; 
  table = ac_empty_table (1, arrayMax (titles), h) ; 
  va_start  (marker, colNum) ;
  nCol = -1 ;
  iC = colNum ;
  while (iC != -1)
    {
      nCol++ ;
      for  (int jj = 0, iR = sR == 0 ? 1 : sR ; iR < eR ; iR++, jj++) /* jump title line */
	{
	  if  (TBB (iR, iC))
	    {
	      /* clean up the table, replace ooo cellends by zero */
	      char *txt = vtxtPtr  (bfr) + TBB  (iR, iC) ;
	      char *ptr = strstr  (txt, ooo) ;
	      if (ptr) *ptr=0 ;
	      /* insert in table */
	      if  (*txt)
		ac_table_insert_text (table, jj, nCol, txt) ;
	    }
	}
      iC= va_arg  (marker, int) ; 
    }
  va_end  ( marker ) ; 
  vtxtPrint (blkp, "\n<br>\n") ;
  int nnn = ac_table_display (blkp, table, arrayp (titles, 0, const char *), 0, 0, 0, 0, 0, 0) ;
  ac_free (h) ;
  return nnn ;
} /* fichePrintSquareTable2 */


/* -===================================- /
* -=  separates the fiche into title, protein title and coment sections
* -===================================- */	

void ficheSectionizeFiche (char *ficheComments, char ** titleFromFiche,  char **commentSection)
{
  char *ptr ; 
  
  if (titleFromFiche)*titleFromFiche=0 ; 
   if (commentSection)*commentSection=0 ; 
  
  if (!ficheComments)return ; 
  
  if (titleFromFiche)*titleFromFiche=ficheComments ; 
  if ((ptr=strstr (ficheComments, "\n\n")))
    {
      *ptr=0 ; 
      if (commentSection)*commentSection=ptr+2 ; 
    }
}

/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  STATEMENTS FUNCTIONS
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

#define	ficheListOperons(v_dict,v_blk,v_gmp,v_tg, v_dist, v_titl, v_strand)	 (ficheTGListOperonsDirect ((v_dict),(v_blk),(v_gmp),(v_tg), (v_dist), 0, (v_titl), (v_strand))+ficheTGListOperonsDirect ((v_dict),(v_blk),(v_gmp), (v_tg), - (v_dist), 0, (v_titl), (v_strand)))


/* -===================================- /
 * -=  TG Best Name                    =- /
 * -===================================- */	

void ficheTGBestNameStatement (vTXT blkp, GMP *gmp, AC_OBJ tg, int isLocalLink)
{
  AC_OBJ oTmp ; 
  AC_HANDLE h = ac_new_handle () ;
  char	linkBuf[vONLINEMAXURLSIZE] ;
  const char *tmpTxt ; 
  
  if (gmp->Spc == HUMAN && (oTmp = ac_tag_obj (tg, "Gene", h)) && 
      ac_tag_printable (oTmp, "LocusLink", 0)) 
    {
      tmpTxt = ac_tag_printable (oTmp, "GeneId", 0) ; 
      if (!tmpTxt) tmpTxt = ac_tag_printable (oTmp, "LocusId", 0) ; 
      if (!tmpTxt) tmpTxt = ac_name (oTmp) ; 
      sprintf (linkBuf, ENTREZ_LINK, tmpTxt) ; 

      if (! (tmpTxt = ac_tag_printable (oTmp, "LocusLink", 0)))
	tmpTxt = ac_name (oTmp) ; 
      if (isLocalLink)
	gmpObjLink (blkp, gmp, tg, tmpTxt) ; 
      else gmpURL (blkp, gmp, linkBuf, tmpTxt) ; 
    }
  else if ( gmp->Spc==HUMAN && (oTmp = ac_tag_obj (tg, "Matching_genefinder_gene", h)))
    {
      tmpTxt = ac_tag_printable (oTmp, "GeneId", 0) ; 
      if (!tmpTxt) tmpTxt = ac_tag_printable (oTmp, "LocusId", 0) ; 
      if (!tmpTxt) tmpTxt = ac_name (oTmp) ; 
      sprintf (linkBuf, ENTREZ_LINK, tmpTxt) ; 
      if (isLocalLink)
	gmpObjLink (blkp, gmp, tg, ac_name (oTmp)) ; 
      else gmpURL (blkp, gmp, linkBuf, ac_name (oTmp)) ; 
    }
  else 
    {
      gmpObjLink (blkp, gmp, tg, 0) ; 
    }
  ac_free (h) ;

  return ; 
}
    
/* -===================================- /
* -=  TG Best Name                    =- /
* -===================================- */	

static int ficheTGTitleStatement (vTXT blkp, GMP *gmp, AC_OBJ tg)
{
  int ok = 0  ; 
  char *cp ;
  AC_OBJ gene = ac_tag_obj (tg, "Gene", 0) ;
  GMP *gmp2 = gmpCreate (gmp->db, gene, 0, 0, 0, 0, gmp->style, 'g') ;

  if (gmp2)
    {
      gtGeneTitle (blkp, gmp2, FALSE) ;
      if ((cp = vtxtPtr (blkp)))
	{      
	  cp = cp + strlen(cp) - 1 ;
	  if (*cp == '.') *cp = 0 ; /* clean up terminal . */
	}
      ok = 1 ;
      gmpDestroy (gmp2) ; /* will kill gene */
    }
  else
    ac_free (gene) ;
  return ok  ; 
} 

/* -===================================- /
* -=  COUNT ATGCs                     =- /
* -===================================- */	

static void ficheCountATGCs (vTXT blkp, char *seq, int leng)
{
  int       ic, cntA, cntT, cntG, cntC, cntAll ; 
  char    *let ; 
  
  for (ic = 0, cntA = cntT = cntG = cntC = 0, let = seq ; *let && ic < leng ; let++, ic++)
    {
      if (*let == 'a')cntA++ ; 
      else if (*let == 't')cntT++ ; 
      else if (*let == 'g')cntG++ ; 
      else if (*let == 'c')cntC++ ; 
    }
  cntAll = cntA+cntT+cntC+cntG ; 
  vtxtPrintf (blkp, 
	     "%d%% A, %d%% T, %d%% G, %d%% C"
	     , (cntA*100)/cntAll
	     , (cntT*100)/cntAll
	     , (cntG*100)/cntAll
	     , (cntC*100)/cntAll
	     ) ; 
}

/* -===================================- /
* -=  MRNA Completeness               =- /
* -===================================- */	

static void ficheNewMrnaCompletenessStatement (vTXT blkp, GMP *gmp)
{
  int isComplete=0 ; 
  int isCompleteCDS=0 ; 
  AC_HANDLE h  = ac_new_handle () ;
  int ng = 0 ;

  if (gmp->mrna)
    {
      ng = ac_keyset_count (ac_objquery_keyset (gmp->mrna, ">Product ; (Best_product && Good_product) || very_good_product", h)) ;
    }

  if (ac_has_tag (gmp->mrna, "Complete") || (ac_has_tag (gmp->mrna, "Found5p") && ac_has_tag (gmp->mrna, "Found3p") ))
    {
       if (ac_keyset_count (ac_objquery_keyset (gmp->mrna, ">product ; mRNA_5p_complete && best_product && !at_position_1", h))) ;
       else
	 isComplete=1 ; 
    }
  if (ac_has_tag (gmp->product, "Complete") ||
      (ac_has_tag (gmp->product, "NH2_Complete") && ac_has_tag (gmp->product, "at_position_1")  && ac_has_tag (gmp->product, "COOH_Complete")))
    isCompleteCDS=1 ; 
  
  if (ng)
    {
      if (isComplete && isCompleteCDS)
	vtxtPrintf (blkp, "complete mRNA") ; 
      else if (isCompleteCDS)
	vtxtPrintf (blkp, "complete CDS mRNA") ; 
      
      else if (ac_has_tag (gmp->product, "NH2_Complete") && ac_has_tag (gmp->product, "at_position_1") && !ac_has_tag (gmp->product, "COOH_Complete") )
	vtxtPrintf (blkp, "partial mRNA") ;  /* "partial mRNA, 3' incomplete" */
      
      else if ((!ac_has_tag (gmp->product, "NH2_Complete") || !ac_has_tag (gmp->product, "at_position_1")) && 
	       ac_has_tag (gmp->product, "COOH_Complete") )
	vtxtPrintf (blkp, "mRNA") ; 
      else 
	vtxtPrintf (blkp, "partial mRNA") ; /* "partial mRNA, internal fragment" */
    }
  else
    {
      if (isComplete && isCompleteCDS)
	vtxtPrintf (blkp, "apparently non-coding mRNA") ; 
      else 
	vtxtPrint (blkp, "partial mRNA") ;
    }
  ac_free (h) ;
  return ;
} /* ficheNewMrnaCompletenessStatement */

/* -===================================- /
* -=  MRNA Alternatives               =- /
* -===================================- */	


/* -===================================- /
* -=  TG locus description            =- /
* -===================================- */	

/* -===================================- /
* -=  READS transplition              =- /
* -===================================- */	

static int ficheDoREADListTransplitionStatement (vTXT blkp, GMP *gmp, AC_OBJ gene, char *prefix, AC_KEYSET kRead, int cntRead, int getSeq, char *separ, BOOL *isSl2sp)
{
  int		ir, whichSL, vectorClip, isAny, length, ixSL ; 
  int		lenSL[128], seqLen, frqSL[128], indSL[128] ; 
  vTXT    blkSL[128] ; 
  AC_OBJ oRead = 0, oSage ;
  AC_ITER iter ;
  AC_TABLE gTanspliced_to, gFrom ; 
  AC_HANDLE h = ac_new_handle () ;
  BOOL ok ;

  memset (lenSL, 0, sizeof (lenSL)) ; 
  memset (blkSL, 0, sizeof (blkSL)) ; 
  memset (frqSL, 0, sizeof (frqSL)) ; 
    
  iter = ac_keyset_iter (kRead, TRUE, 0) ;
  isAny = 0 ;
  while (ac_free (oRead), (oRead = ac_next_obj (iter)))
    {
      if (! (gTanspliced_to = ac_tag_table (oRead, "Transpliced_to", h)))
	continue ; 
      vectorClip = ac_tag_int (oRead, "Vector_Clipping", 0) ; 
      
      ok = FALSE ;
      gFrom = ac_tag_table (oRead, "From_gene", h) ;
      for (ir = 0 ; !ok && gFrom && ir < gFrom->rows ; ir++)
	if (!strcmp (ac_name(gene), ac_table_printable (gFrom, ir, 0, "")))
	  {
	    if (ac_table_int (gFrom, ir, 2, 1000) < 3)
	      ok = TRUE ;
	  }
      if (!ok)
	continue ;
      for (ir = 0 ; ir < gTanspliced_to->rows ; ir++)
	{
	  whichSL = 0 ; 
	  if (!sscanf (ac_table_printable (gTanspliced_to, ir, 0, ""), "SL%d", &whichSL) || !whichSL)
	    continue  ; 
	  if (whichSL > 1)
	    *isSl2sp = TRUE  ; 
	  length = ac_table_int (gTanspliced_to, ir, 1, 0) - vectorClip + 1 ; 
	  if (length>lenSL[whichSL])lenSL[whichSL] = length ; 
	  if (!blkSL[whichSL])
	    {
	      blkSL[whichSL] = vtxtCreate () ;
	      if (gmp->markup) vtxtMarkup (blkSL[whichSL]) ;
	    }
	  if (frqSL[whichSL] < 12)
	    {
	      frqSL[whichSL]++ ; 
	      if (vtxtPtr ( (blkSL[whichSL])))vtxtPrintf (blkSL[whichSL], ", ") ; 
	      vtxtPrintf ( (blkSL[whichSL]), "%s", gtYkName (ac_name (oRead))) ; 
	      /* gmpObjLink (& (blkSL[whichSL]), oRead, 0) ;  */
	    }
	  else if (frqSL[whichSL]++ == 12) 
	    vtxtPrintf ( (blkSL[whichSL]), "...") ;	  
	}
    }
  ac_free (oRead) ;
  ac_free (iter) ;

  /* sorting*/
  {
    int lSL, tmp ; 
    for (whichSL = 0 ; whichSL< 128 ; whichSL++)indSL[whichSL]=whichSL ;  /* indexing */
    for (whichSL = 0 ;  whichSL  < 128 ;  whichSL++)
      {
	for (lSL = 0 ;  lSL < 128 ;  lSL++)
	  {
	    if (frqSL[indSL[whichSL]] < frqSL[indSL[lSL]])
	      {
		tmp=indSL[whichSL] ; 
		indSL[whichSL]=indSL[lSL] ; 
		indSL[lSL]=tmp ; 
	      }
	  }
      }
  }
  
  isAny = 0 ; 
  for (ixSL = 128 - 1 ;  ixSL >= 1 ;  ixSL --)
    {
      whichSL=indSL[ixSL] ; 
      if (lenSL[whichSL] > 0)
	{
	  vtxtPrintf (blkp, "%s%sSL%d", prefix, (isAny>0 ? separ : " "), whichSL) ; 
	  prefix = "" ; 
	  if (getSeq) 
	    {
	      const char *sequenceSL ;
	      oSage = ac_get_obj (gmp->db, "Sage"
				  , messprintf ("SL%d", whichSL)
				  , 0) ;
	      if (oSage && (sequenceSL = ac_tag_printable (oSage, "Motif", 0)))
		{
		  seqLen = strlen (sequenceSL) ; 
		  if (lenSL[whichSL]>seqLen)
		    lenSL[whichSL] = seqLen ;
		  if (!strncmp (vtxtPtr(blkSL[whichSL]), "XS_",3) && lenSL[whichSL]< 10)
		    lenSL[whichSL] = 10 ;
		  vtxtPrintf (blkp, " (%s visible in"
			      , sequenceSL+seqLen-lenSL[whichSL]
			      ) ; 
		  if (!blkSL[whichSL])
		    {
		      blkSL[whichSL] = vtxtCreate () ;
		      if (gmp->markup) vtxtMarkup (blkSL[whichSL]) ;
		    }
		  if (vtxtPtr (blkSL[whichSL]))
		    vtxtPrintf (blkp, " %s", vtxtPtr (blkSL[whichSL])) ; 
		  vtxtPrintf (blkp, ") ") ; 
		}
	      ac_free (oSage) ; 
	    }
	  isAny++ ; 
	}
    }
  for (whichSL = 0 ;  whichSL  < 128 ;  whichSL++)
    {
      vtxtDestroy (blkSL[whichSL]) ; 
    }
  
  ac_free (h) ;
  return isAny ; 
}

static int ficheREADListTransplitionStatement (vTXT blkp, GMP *gmp, AC_OBJ gene, char *prefix, AC_KEYSET kRead, int cntRead, int getSeq, char *separ, BOOL *ip)
{
  int ok = 0  ; 
  char *ptr  ; 

  vTXT buf = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (buf) ;

  ficheDoREADListTransplitionStatement (buf, gmp, gene, prefix, kRead, cntRead, getSeq, separ, ip) ; 
  
  if ((ptr = vtxtPtr (buf)))
    {
      vtxtPrintf (blkp, ptr) ;
      ok = 1  ; 
    }
  vtxtDestroy (buf) ; 
  return ok  ; 
}


/* -===================================- /
* -=  PSORT Aliases                   =- /
* -===================================- */	

static const char *fichePSORTNameAlias (const char *nm)
{
  struct {char *find, *repl ; } pSortAliasStruct[]={
  {"cytoplasmic", "in the cytoplasm"}, 
  {"mitochondrial", "in the mitochondria"}, 
  {"nuclear", "in the nucleus"}, 
  {"endoplasmic_reticulum", "in the endoplasmic reticulum"}, 
  {"endoplasmic reticulum", "in the endoplasmic reticulum"}, 
  {"vesicles_of_secretory_system", "in vesicles of the secretory system"}, 
  {"vesicles of secretory_system", "in vesicles of the secretory system"}, 
  {"Golgi", "in the Golgi"}, 
  {"Golgiapparatus", "in the Golgi apparatus"},
  {"Membrane", "in a membrane"}, 
  {"peroxisomal", "in peroxisomes"}, 
  {"cytoskeletal", "in the cytoskeleton"}, 
  {"plasma_membrane", "in the plasma membrane"}, 
  {"plasma membrane", "in the plasma membrane"}, 
  {"vacuolar", "in vacuoles"}, 
  {"extracellular, including cell wall", "extracellular or secreted"}, 
  {"Secreted", "secreted"}, 
  {0, 0}
} ; 
  int iFnd ; 
  
  for (iFnd=0 ;  pSortAliasStruct[iFnd].find ; iFnd++)
    if (!strcasecmp (nm,  pSortAliasStruct[iFnd].find))
      break ; 
  
  return pSortAliasStruct[iFnd].find ? pSortAliasStruct[iFnd].repl : nm ; 
}



/*_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  PARAGRAPHS OF FICHES
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


/* -===================================- /
* -=  MRNA Title                      =- /
* -===================================- */	

int ficheMrnaTitleParagraph (vTXT blkp, GMP *gmp)
{
  char *ptr ;
  vTXT buf = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (buf) ;

  vtxtHr (blkp, 0, 0) ;
  if (gmp->pg)
    {
      if (ac_has_tag (gmp->pg, "CDS"))
	gtPredictedMrnaTitle (buf, gmp) ; 
      else
	gtPredictedTrnaTitle (buf, gmp, TRUE, 0, 0) ;
    }
  else if (gmp->mrna)
    gtMrnaTitle (buf, gmp) ; 
  
  if ((ptr = vtxtPtr (buf))) /* contains a dot. nasty when bold */
    {	
      if (!strstr(ptr, messprintf("(%s)", ac_name (gmp->gene))) &&
	  !strstr(ptr, messprintf(">%s<", ac_name (gmp->gene))))
	vtxtPrintf (buf, " (%s)", ac_name (gmp->gene)) ;
      vtxtDot (buf) ;
      ptr = vtxtPtr (buf) ; /* could be reallocated */
      vtxtBreak (blkp) ;
      vtxtBold (blkp, ptr) ;
      vtxtPrintf (blkp,"\n") ; /* avoid a second dot */
    }
  vtxtHr (blkp, 0, 0) ;
     
  vtxtDestroy (buf) ;
  return 1 ;
}

/* -===================================- /
* -=  MRNA Overview                   =- /
* -===================================- */	
/*****************************************************************************/

static BOOL ficheNewMrnaOverviewCompletenessAndSupport (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  int n, nbp = 0 ;
  BOOL ok = FALSE ;
  AC_TABLE oDNA ;
  AC_KEYSET clones ;

  oDNA = ac_tag_table (gmp->mrna, "DNA", h) ;
  
  if (oDNA)
    nbp = ac_table_int (oDNA, 0, 1, 0) ;
  vtxtDot (blkp) ;
  
  /* completeness */
  vtxtPrintf (blkp,"This ") ;
  ficheNewMrnaCompletenessStatement (blkp, gmp) ;

  /* length */
  vtxtPrintf (blkp," is ") ;

  if (blkp->markUp)
    vtxtPrintf (blkp,"<a %s>%d bp</a>"
		, "href=\"javascript:openAnchor(0,'mRNA_sequences')\"" 
		, nbp) ;
  else
    vtxtPrintf (blkp,"%d bp", nbp) ;

  vtxtPrintf (blkp," long") ;
  
  /* clone support */
  clones = ac_objquery_keyset (gmp->mrna, ">cdna_clone", h) ;
  n = ac_keyset_count (clones) ;
  if (n) 
    {
      if (blkp->markUp)
	vtxtPrintf (blkp,". It is reconstructed from %s <a href=\"javascript:openAnchor(0,'mRNA_expression')\">cDNA clone%s</a>", isOne (n), _multi (n)) ;
      else
	vtxtPrintf (blkp,". It is reconstructed from %s cDNA clone%s", isOne (n), _multi (n)) ;
    
      /* tissues */
      {
	ficheNewGeneExpressionTissue (blkp, gmp, clones, 5, 8, gmp->view, ac_keyset_count (clones)) ;
      }
    }
  ac_free (h) ;

  return ok ;
} /* ficheNewMrnaOverviewCompletenessAndSupport */

/*********************************************************************/
/*********************************************************************/

typedef struct cloSup_struct { int nam, delta, nv, type, lib ; } CLOSUP ;

static int cloSupOrder (const void *va, const void *vb)
{
  const CLOSUP *up = (const CLOSUP*) va ;
  const CLOSUP *vp = (const CLOSUP*) vb ;
  int n ;
  
  if (up->lib != vp->lib) 
    return up->lib - vp->lib ; 
  n = up->delta + up->nv - vp->delta - vp->nv ;
  if (n) return n ;
  return up->nam - vp->nam ;
}

/*******************/

static int ficheNewMrnaCloneSupport (vTXT blkp, GMP *gmp)
{
  int nnLib[12] ;
  int nn = 0, nn1, ii, ir, lib ;
  int nMgc_r = 0, nMgc = 0, nDkfz = 0, nFlj = 0, nKiaa = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  Array reads =  arrayHandleCreate (32, CLOSUP, h) ;
  AC_OBJ Est = 0, Clone ;
  AC_TABLE gSup = ac_tag_table (gmp->product, "Covered_by", h) ;
  CLOSUP *up ;
  DICT *dict = dictHandleCreate (32, h) ;
  char *libNam [7] ;
  const char *cp, *prefix ;
  AC_OBJ Seq ;
  char	linkBuf[vONLINEMAXURLSIZE] ; 

  memset (nnLib, 0, sizeof (nnLib)) ;

  if (gmp->markup)
    {
      libNam [0] = "" ;
      libNam [1] = "NM" ;
      libNam [2] = "<a href=\"http://www.kazusa.or.jp/huge/clone.req\">KIAA</a>" ;
      libNam [3] = "<a href=\"http://www.rzpd.de\">DKFZ</a>" ;
      libNam [4] = "<a href=\"http://www.nbrc.nite.go.jp/e/hflcdna2-e.html\">FLJ</a>" ;
      libNam [5] = "<a href=\"https://mgc.nci.nih.gov\">MGC</a>" ;
      libNam [6] = "<a href=\"https://mgc.nci.nih.gov\">MGC'</a>" ;
    }
  else
    {
      libNam [0] = "" ;
      libNam [1] = "NM" ;
      libNam [2] = "KIAA" ;
      libNam [3] = "DKFZ" ;
      libNam [4] = "FLJ" ;
      libNam [5] = "MGC" ;
      libNam [6] = "MGC'" ;
    }
  dictAdd (dict, "toto", &ii) ; /* avoid zero */
  for (ir = nn = 0 ; gSup && ir < gSup->rows ; ir++)
    {
      cp = ac_table_printable (gSup, ir, 0, 0) ;
      if (!cp)
	continue ;
      dictAdd (dict, cp, &ii) ;
      up = arrayp (reads, nn++, CLOSUP) ;
      up->nam = ii ;
      up->delta = ac_table_int (gSup, ir, 2, -1) ;
      up->nv = ac_table_int (gSup, ir, 4, -1) ;
      cp = ac_table_printable (gSup, ir, 5, 0) ; 
      if (cp)
	{
	  dictAdd (dict, cp, &ii) ;
	  up->type = ii ;
	}
      cp = ac_table_printable (gSup, ir, 0, 0) ;
      if (!strncmp (cp, "NM_", 3))
	{ up->lib = 1 ; }
      else
	if (gmp->Spc == HUMAN && (!nMgc_r || !nMgc || !nDkfz || !nFlj || !nKiaa))
	  {
	    Est = ac_table_obj (gSup, ir, 0, h) ;
	    Clone = Est ? ac_tag_obj (Est, "cDNA_clone", h) : 0 ;
	    if (Clone && ac_has_tag (Clone, "KIAA"))
	      { nKiaa++ ; up->lib = 2 ; }
	    if (Clone && ac_has_tag (Clone, "DKFZ"))
	      { nDkfz++ ; up->lib = 3 ; }
	    if (Clone && ac_has_tag (Clone, "FLJ"))
	      {  nFlj++ ; up->lib = 4 ; }
	    if (Clone && ac_has_tag (Clone, "MGC_r")) 
	      { nMgc_r++ ; up->lib = 5 ; }
	    else if (Clone && ac_has_tag (Clone, "MGC"))
	      { nMgc++ ; up->lib = 6 ; }
	    ac_free (Est) ; ac_free (Clone) ;
	  }
    }
  arraySort (reads, cloSupOrder) ;
  nn = 0 ; 
  prefix = ". It is exactly encoded by " ;

  if (gmp->Spc == HUMAN )
    for (nn = 0, lib = 1 ; lib < (nMgc_r > 0 ? 6 : 7) ; lib++)
      {
	nnLib [lib] = 0 ;
	for (ii = 0 ; nnLib[lib] == 0  && ii < arrayMax (reads) ; ii++)
	  { 
	    up = arrp (reads, ii, CLOSUP) ;
	    if (up->lib == lib && up->delta == 0 && up->nv == 0)
	      {
		if (nn++)
		  vtxtPrint (blkp, ", ") ;
		else
		  vtxtPrint (blkp, prefix) ; 
		nnLib[lib]++ ;

		Seq = ac_get_obj (gmp->db, "Sequence", dictName (dict, up->nam), 0) ;
		sprintf (linkBuf, GENBANK_LINK, dictName (dict, up->nam)) ;
		gmpURL (blkp, gmp, linkBuf, dictName (dict, up->nam)) ;
		ac_free (Seq) ;

		if (lib > 1)
		  vtxtPrintf (blkp, " (available from %s)", libNam[up->lib]) ; 
		up->lib = 99 ; /* consume this clone */
	      }
	  }
      }
  /* nothing from the main libs */
  for (ii = 0 ; nn < 5 && ii < arrayMax (reads) ; ii++)
    { 
      up = arrp (reads, ii, CLOSUP) ;
      if (! up->lib && up->nam && up->delta == 0 && up->nv == 0)
	{
	  if (nn++)
	    vtxtPrint (blkp, ", ") ;
	  else
	    vtxtPrint (blkp, prefix) ; 
	  up->lib = 99 ; /* consume this clone */
	  
	  if (gmp->Spc != WORM) /* not worm */
	    {
	      Seq = ac_get_obj (gmp->db, "Sequence", dictName (dict, up->nam), 0) ;
	      sprintf (linkBuf, GENBANK_LINK, dictName (dict, up->nam)) ;
	      gmpURL (blkp, gmp, linkBuf, dictName (dict, up->nam)) ;
	      ac_free (Seq) ;
	    }
	  else
	    vtxtPrintf (blkp, " %s", dictName (dict, up->nam)) ;
	  if (nn >= 5)
	    break ;
	}
    }

  if (nn) 
    prefix = ", and with some amino-acid variations by " ;
  else
    prefix = ". It is encoded with some amino-acid variations by " ;
  nn1 = nn ;
  nn = 0 ;
  if (gmp->Spc == HUMAN )
    {
      for (lib = 1 ; nn < 5 && lib < (nnLib[5] ? 6 : 7) ; lib++)
	for (ii = 0 ; nnLib[lib] == 0 && ii < arrayMax (reads) ; ii++)
	  { 
	    up = arrp (reads, ii, CLOSUP) ;
	    if (lib == up->lib)
	      {
		if (nn++)
		  vtxtPrint (blkp, ", ") ;
		else
		  vtxtPrint (blkp, prefix) ; 
		nnLib[lib]++ ;

		Seq = ac_get_obj (gmp->db, "Sequence", dictName (dict, up->nam), 0) ;
		sprintf (linkBuf, GENBANK_LINK, dictName (dict, up->nam)) ;
		gmpURL (blkp, gmp, linkBuf, dictName (dict, up->nam)) ;
		ac_free (Seq) ;

		if (up->type)
		  vtxtPrintf (blkp, " [%s]", dictName (dict, up->type)) ;
		if (lib > 1)
		  vtxtPrintf (blkp, " (available from %s)", libNam[up->lib]) ; 
		up->lib = 99 ; /* consume this clone */
	      }
	  }
    }
  /* nothing from the main libs */
  for (ii = 0 ; nn < 5 && ii < arrayMax (reads) ; ii++)
    { 
      up = arrp (reads, ii, CLOSUP) ;
      if (! up->lib && up->nam)
	{
	  if (nn++)
	    vtxtPrint (blkp, ", ") ;
	  else
	    vtxtPrint (blkp, prefix) ; 
	  up->lib = 99 ; /* consume this clone */

	  if (gmp->Spc != WORM)
	    {
	      Seq = ac_get_obj (gmp->db, "Sequence", dictName (dict, up->nam), 0) ;
	      sprintf (linkBuf, GENBANK_LINK, dictName (dict, up->nam)) ;
	      gmpURL (blkp, gmp, linkBuf, dictName (dict, up->nam)) ;
	      ac_free (Seq) ;
	    }
	  else
	    vtxtPrintf (blkp, " %s", dictName (dict, up->nam)) ;
	  if (up->type)
	    vtxtPrintf (blkp, " [%s]", dictName (dict, up->type)) ;
	  if (nn >= 5)
	    break ;
	}
    }
  nn1 += nn ;
  nn = arrayMax(reads) - nn1 ;
  if (nn && nn1)
    {
      vtxtPrintf (blkp, " and %d other sequence%s", nn, _multi(nn)) ;
    }
  
  ac_free (h) ;

  return nn ;
} /* ficheNewMrnaCloneSupport */

/*********************************************************************/
/*********************************************************************/

static int ficheMRNAOverviewParagraphContent (vTXT blkp, GMP *gmp)
{
  AC_TABLE gExpasy, gIntMap ;
  AC_TABLE gTranspliced_to ; 
  AC_TABLE gSplicing ; 
  AC_TABLE gPfam/*, * gLocalization*/ ; 
  const char *prefix ;
  const char *txtAccession ;
  char *ptr, familyName[vSTRMAXNAME], *extName, linkBuf[vONLINEMAXURLSIZE] ; 
  int x1, x2, dx=0, nPfam = 0  ; 
  int nValid3p = 0, nValid5p = 0 ;
  int  ir, Length_3prime_UTR, Length_5prime_UTR, nexons = 0, level=0 ; 
  BOOL isSl2s = FALSE,  hasPsort = FALSE,  hasPfam = FALSE ;
  AC_HANDLE h = ac_new_handle () ;
  const char * UTR ;
  const char * lExon ; 
  const char * lExons ;
  extName = strrchr (ac_name (gmp->mrna), '.') ;  if (!extName) extName = "" ;  else extName ++ ; 
  
  UTR = blkp->markUp ? "<a href=\"javascript:openAnchor('fmol','mRNAs')\">UTR</a>" : "UTR" ;
  lExon = blkp->markUp ? "<a href=\"javascript:openAnchor(0,'tg_introns')\">exon</a>" : "exon" ;
  lExons = blkp->markUp ? "<a href=\"javascript:openAnchor(0,'tg_introns')\">exons</a>" : "exons" ;

  if (! ac_has_tag (gmp->product, "Best_product"))
    {
      vtxtPrintf (blkp, "We annotate here a putative secondary product of this mRNA") ;
      goto laba ;
    }
  if (gmp->tg &&
      ficheNewMrnaOverviewCompletenessAndSupport (blkp, gmp))
    vtxtDot (blkp) ;
  
  if (gmp->mrna && gmp->tg)
    {
      vtxtDot (blkp) ;
      prefix = ". It" ; 
      if ((gTranspliced_to = ac_tag_table (gmp->mrna, "Transpliced_to", h)))
	{ /* TRANSPLICED INFO */
	  AC_KEYSET kRead ;
	  
	  kRead = ac_objquery_keyset (gmp->mrna
				   , messprintf ("Follow cDNA_clone  ;  Follow Read  ;  From_gene == \"%s\" " 
						 , ac_name (gmp->tg)) 
				   , 0) ;
	  if (kRead && ac_keyset_count (kRead))
	    {
	      ficheREADListTransplitionStatement 
		(blkp, gmp, gmp->gene, "It is transpliced to", kRead, 0, 0, " or ", &isSl2s) ; 
	    }
	  ac_free (kRead) ; 
	}
      else 
	{
	  if (ac_has_tag (gmp->mrna, "Found5p"))
	    {
	      if (gmp->Spc == WORM)
		{ 
		  vtxtDot (blkp) ;
		  vtxtPrintf (blkp, " It does not appear to be transpliced to a leader") ; 
		}
	    }
	  else if (gmp->product && !ac_has_tag (gmp->product, "Complete") && !ac_has_tag (gmp->product, "NH2_Complete") )
	    {
	      vtxtDot (blkp) ;
	      vtxtPrintf (blkp, " It may be incomplete at the 5' end") ; /* N terminus */
	      prefix=" and" ; 
	    }
	}
      if (gmp->product && 
	  !ac_has_tag (gmp->product, "Complete") &&
	  !ac_has_tag (gmp->product, "COOH_Complete") )
	{
	  if (*prefix == '.')
	    vtxtDot (blkp) ;
	  vtxtPrintf (blkp, " It is incomplete at the 3' end") ; /*  C terminus */
	  prefix=", " ; 
	}
    }
  prefix=". " ; 
  if (gmp->mrna)
    {
      nValid3p = ac_tag_count (gmp->mrna, "Valid3p") ;
      nValid5p = ac_tag_count (gmp->mrna, "Valid5p") ;
      if (nValid5p)
	{
	  if (ac_keyset_count (ac_objquery_keyset (gmp->mrna, ">product ; mRNA_5p_complete && best_product && !at_position_1", h))) 
	    nValid5p = 0 ;
	}
      if (nValid5p > 1 ||  nValid3p > 1)
	{
	  vtxtDot (blkp) ;
	  vtxtPrintf (blkp, "It is a composite of similar variants differing only in the length of their %ss: ", UTR) ;
	  
	  if (nValid5p > 1)
	    vtxtPrintf (blkp, "there are %d validated transcription start sites", nValid5p) ;
	  if (nValid3p > 1)
	    vtxtPrintf (blkp, "%s %d alternative polyadenylation sites"
		       , nValid5p > 1 ? " and" : "there are"
		       , nValid3p
			) ;
	}
    }

  if (gmp->mrna)
    {
      if ((Length_5prime_UTR = ac_tag_int (gmp->mrna, "Length_5prime_UTR", 0)) >= SpcI[gmp->Spc].long5P)
	{
	  vtxtDot (blkp) ;
	  if ((Length_3prime_UTR = ac_tag_int (gmp->mrna, "Length_3prime_UTR", 0)) >= SpcI[gmp->Spc].long3P)
	    vtxtPrintf (blkp, " Both the 3\' and 5\' %s are very long", UTR) ;
	  else
	    vtxtPrintf (blkp, " The 5\' %s is very long", UTR) ;
	}
      else if ((Length_3prime_UTR = ac_tag_int (gmp->mrna, "Length_3prime_UTR", 0)) >= SpcI[gmp->Spc].long3P)
	{
	  vtxtDot (blkp) ;
	  vtxtPrintf (blkp, " The %s3\' UTR is very long", nValid3p > 1 ? "longest " : "") ;
	}
    }
 

  if (gmp->mrna && ac_has_tag (gmp->mrna, "NMD"))
    {  
      vtxtDot (blkp) ;
      vtxtPrint (blkp, "This mRNA could be a target of nonsense mediated mRNA decay") ;
    }

 
  if (1)
    { 
      gSplicing = ac_tag_table (gmp->mrna, "Splicing", h) ; 
      
      for (ir = nexons = 0 ; gSplicing && ir < gSplicing->rows ; ir++)
	if (strstr (ac_table_printable (gSplicing, ir, 4, ""), "xon"))
	  nexons++ ;
      gIntMap = gmp->mrna ? ac_tag_table (gmp->mrna, "IntMap", h) : ac_tag_table (gmp->pg, "IntMap", h) ;
      if (gIntMap && gIntMap->rows && gIntMap->cols >= 3)
	{
	  x1 = ac_table_int (gIntMap, 0, 1, 0) ;
	  x2 = ac_table_int (gIntMap, 0, 2, 0) ;
	  dx = x2 - x1 ;
	  if (dx < 0) dx = - dx ;
	  dx++ ; 
	}
     vtxtDot (blkp) ;
      vtxtPrintf (blkp, " The") ;

      if (gmp->tg)
	level = 4 ;
      else if (gmp->pg && ac_has_tag (gmp->pg, "CDS"))
	level = gtPredictedMrnaSupport (gmp->pg, 0) ;
      if (level < 4)
	vtxtBold (blkp, " predicted") ;
      vtxtPrintf (blkp, " %s"
		  , gmp->tg ? "premessenger" : "CDS") ;

      if (nexons <= 1)
	vtxtPrintf (blkp, " has a single %s", lExon) ;
      else
	vtxtPrintf (blkp, " has %d %s", nexons, lExons) ;
      if (gmp->pg)
	{
	  AC_TABLE gDNA = ac_tag_table (gmp->mrna, "DNA", h) ;
  
	  if (gDNA)
	    {
	      int nbp = ac_table_int (gDNA, 0, 1, 0) ;
	      vtxtDot (blkp) ;
	      vtxtPrintf (blkp, " It is %d bp long", nbp) ;
	    }
	  vtxtPrintf (blkp, " and covers ") ;
	  ficheNewGeneCoversValue (blkp, gmp, 4) ;
	}
      else
	{
	  vtxtPrintf (blkp, " and covers ") ;
	  ficheNewGeneCoversValue (blkp, gmp, 3) ;
	}
      vtxtPrintf (blkp, " on the %s%sgenome"
		  , dbNam ? dbNam : "" , dbNam ? " " : ""
		  ) ;
    }

laba:
  if ( gmp->product)
    {
      AC_TABLE tbl = ac_tag_table (gmp->product, "Peptide", h) ;
      int len = tbl ? ac_table_int (tbl, 0, 1, 0) : 0 ; 
      BOOL isGood = gmp->product && ac_has_tag (gmp->product, "Good_product")  ;
      int score =  gmp->product ? ac_tag_int (gmp->product, "Quality", 0) : -99 ;
      char *nScore [] = { "poor", "good", "very good"} ;
      int nice = 0 ;
      int codingLn =  gmp->product ? ac_tag_int (gmp->product, "Coding_length", 0) : 0 ;
      
      /* 2021_01_28: some proteins (CCND2.a) are annotated with a negative kozak */
      codingLn /= 3 ;
      if (codingLn && codingLn < len)
	len = codingLn ;
      if (isGood)
	{
	  nice = 1 ;
	  if (ac_has_tag (gmp->product, "Very_good_product"))
	    nice = 2 ;
	}
      
      if (gmp->tg &&
	  !ac_has_tag (gmp->product, "Complete") &&
	  !ac_has_tag (gmp->product, "COOH_Complete"))
	prefix= ". The predicted partial protein" ;
      else 
	{
	  if (nice)
	    prefix = ". The predicted protein" ; 
	  else
	    prefix = ". It appears to be non coding, the best  predicted protein" ;
	}  
      vtxtPrint (blkp, prefix) ;
      if (ac_has_tag (gmp->product, "Coding_gap"))
	vtxtPrintf (blkp," has a sequence gap") ;
      else
	{
	  char old = gmp->view ;
	  gmp->view = 'm' ;
	  if (ficheNewGeneKozakStatement (blkp, gmp))
	    vtxtPrint (blkp, ",") ;
	  gmp->view = old ;
	  
	  vtxtPrintf (blkp," %s ", isGood ? "has" : "would have") ;
	  if (blkp->markUp)
	    vtxtPrintf (blkp,"<a %s>%d aa</a>"
			, "href=\"javascript:openAnchor(0,'mRNA_sequences')\"" 
			, len) ;
	  else
	    vtxtPrintf (blkp,"%d aa", len) ;
	    
	  /*
	    ,acTag (gmp->gene, "Use_AM")  
	    ? "derived for the cDNA consensus" 
	    : "derived from the genome sequence"
	  */
	  
	  if ((gExpasy = ac_tag_table (gmp->product, "Expasy", h)))
	    vtxtPrintf (blkp
			," (%.1f kDa, pI %.1f)" 
			, ac_table_float (gExpasy, 0, 0, 0)
			, ac_table_float (gExpasy, 0, 1, 0)
			) ; 
	  if (! ac_has_tag (gmp->product, "Complete") && nice > 0)
	    vtxtPrintf (blkp, " and a %s coding score (%d)", nScore[nice], score) ;
	}
      
      prefix = ". It" ; 
      
      /* tout buggue 
	 if (gmp->tg &&
	 ac_has_tag (gmp->product, "First_Kozak") &&
	 ac_has_tag (gmp->product, "First_atg") &&
	 (ac_has_tag (gmp->product, "Complete") || ac_has_tag (gmp->product, "NH2_Complete")) &&
	 (atg = ac_tag_int (gmp->product, "First_ATG", "")) > 1
	 )
	 {
	 AC_TABLE gNtg = ac_tag_table (gmp->product, "First_Kozak", h) ;
	 char *codon = ac_table_printable (gNtg, 0, 1, 0) ;
	 
	 vtxtPrintf (blkp," if we use the first putative initiator codon%s%s following the Kozak consensus (%d aa less if we use the first AUG)"
	 , *codon ? " " : "", codon, atg - 1) ;
	 prefix = ". It" ; 
	 }
	 else
	 prefix = ", it" ;
      */ 
      ac_free (tbl) ;
    
  
      if ( 
	  (
	   ac_has_tag (gmp->product, "very_good_product") ||
	   (ac_has_tag (gmp->product, "Best_product") && ac_has_tag (gmp->product, "good_product"))
	   ) &&
	  ac_has_tag (gmp->product, "Covered_by") &&
	  ! ac_has_tag (gmp->product, "Coding_gap"))
	{
	  ficheNewMrnaCloneSupport (blkp, gmp) ;
	  vtxtBreak (blkp) ;
	  if (strlen(prefix) > 2 && *(prefix+1)==' ') prefix = prefix+2 ;
	}
      hasPsort = hasPfam = FALSE ;
      if (gmp->product && (gPfam = ac_tag_table (gmp->product, "Pfam", h)))
	{
	  int     cntType, jr ; 
	  AC_OBJ oPfam ; 
	  const char *lastName = "" ;
	  
	  vtxtPrintf (blkp, "%s contains", prefix) ;
	  prefix = "" ; 
	  for (ir = 0 ; ir < gPfam->rows ; ir++)
	    {
	      oPfam = ac_table_obj (gPfam, ir, 0, h) ; 
	      
	      if (strcmp (lastName, ac_name(oPfam)))
		{
		  lastName = ac_name (oPfam) ;
		  for (cntType = 1, jr = 1 
			 ; 	 (jr+ir) <  (gPfam->rows) && 
			 !strcmp (lastName, ac_table_printable (gPfam, ir + jr, 0, ""))
			 ; jr++ )
		    {
		      cntType++ ; 
		    }
		}
	      else continue ; 
	      
	      nPfam++ ;
	      strcpy (familyName, ac_tag_printable (oPfam, "Definition", "")) ; 
	      if (!*familyName)
		{
		  /* try to get it from a synonymous pfam */
		  AC_KEYSET ks = ac_objquery_keyset (oPfam, ">accession; >quoted_in; Definition", h) ;
		  AC_ITER iter = ac_keyset_iter (ks, TRUE, h) ;
		  AC_OBJ oPfam2 = ac_next_obj (iter) ;
		  
		  strcpy (familyName, ac_tag_printable (oPfam2, "Definition", "")) ;
		  ac_free (oPfam2) ;
		}
	      if (cntType == 1)
		vtxtPrintf (blkp, "%s one ", _comma (ir)) ; 
	      else 
		vtxtPrintf (blkp, "%s %d ", _comma (ir), cntType) ; 
	      
	      /*            gmpObjLink (blkp, gmp, oPfam, familyName) ; */
	      
	      if (gmp->markup &&
		  (txtAccession = ac_tag_printable (oPfam, "Accession", 0)))
		{
		  char *cq = strstr (txtAccession , ".") ; /* drop pfam version number */
		  if (cq) *cq = 0 ;
		  sprintf (linkBuf, PFAM_LINK, txtAccession) ; 
		  gmpURL (blkp, gmp, linkBuf, familyName) ; 
		}
	      else if (*familyName)
		vtxtPrint (blkp, familyName) ; 
	      else 
		vtxtPrint (blkp, "Pfam") ; 
	      
	      
	      if (!strstr (familyName, "domain") && !strstr (familyName, "motif"))
		vtxtPrintf (blkp, " domain%s", _multi (cntType)) ; 
	    } 
	  hasPfam = TRUE ;      
	}
      if (hasPfam) vtxtPrint (blkp, " [Pfam]") ;
      if (ac_has_tag (gmp->product, "Psort_domain"))
	{
	  int ii, jj ;
	  AC_TABLE gPsort ;
	  char **cpp, *psNN [] =
	    { 
	      "N_terminal_signal_domain", "N-terminal peptide signal",

	      "Transmembrane_domain", "transmembrane domain", 
	      "N_myristoylation_domain", "N-myristoylation domain",
	      "prenylation_domain", "prenylation domain",
	      "Golgi_transport_domain", "Golgi transport domain", 
	      
	      "peroxisomal_domain", "peroxisomal domain", 
	      "2nd_peroximal_domain", "second peroximal domain", 
	      "vacuolar_domain", "vacuolar domain", 
	      
	      "ER_retention_domain", "ER_retention domain", 
	      "Coiled_coil_region", "coiled coil stretch", 
	      "Leucine_zipper_domain", "leucine zipper domain",
	      "actin_binding_1_domain", "actin binding domain", 
	      "actin_binding_2_domain", "actin binding domain", 
	      "ATP_binding_domain", "ATP binding domain",   
	      
	      "RNA_binding_domain", "RNA binding domain",                           
	      0, 0
	    } ;
	  
	  prefix = nPfam ? "It also" : "It" ;
	  if (nPfam) jj = 1 ; /* implies a simple comma */
	  else jj = 0 ;

	  for (ii = 0, cpp = psNN ; *cpp ; ii++, cpp += 2)
	    if ((gPsort = ac_tag_table (gmp->product, *cpp, h)))
	      {
		if (!jj++)
		  {
		    vtxtDot (blkp) ;
		    vtxtPrintf (blkp, "%s contains", prefix) ;
		  }
		else
		  vtxtPrintf (blkp, ", ") ;
		if (gPsort->rows > 1)
		  vtxtPrintf (blkp, " %s%d", ii==0 ? " at least " : "", gPsort->rows) ;
	    else
	      vtxtPrintf (blkp, " a%s", isVoyelle (**(cpp+1)) ? "n" : "") ;
		vtxtPrintf (blkp, " %s", *(cpp+1)) ;
		if (gPsort->rows > 1)
		  vtxtPrintf (blkp, "s") ;
		hasPsort = TRUE ;
	      }
	  if (hasPsort) vtxtPrintf (blkp, " [Psort2]") ;
	}
      
      if (!hasPfam && gmp->kantor && ac_has_tag (gmp->kantor, "Pfam_date"))
	{
	  vtxtDot (blkp) ;
	  if (!hasPsort &&
	      ac_has_tag (gmp->kantor, "Psort_date"))
	    vtxtPrintf (blkp, "%s contains no protein domain or characteristic Psort motif", prefix) ; 
	  else
	    vtxtPrintf (blkp, "%s contains no Pfam protein domain", prefix) ; 
	}
      
      if ((prefix = ac_tag_printable (gmp->product, "Psort_title", 0)))
	{
	  const char * lLoc = "href=\"javascript:openAnchor(0,'Psort')\"" ;

	  vtxtDot (blkp) ;
	  if (strstr(fichePSORTNameAlias (prefix),"ecreted"))
	    vtxtPrintf (blkp
			, "It is predicted to be <a %s>%s</a> [Psort2]"
			, lLoc
			, fichePSORTNameAlias (prefix)
			) ; 
	  else
	    vtxtPrintf (blkp
			, "It is predicted to <a %s>localize</a> %s [Psort2]"
			, lLoc
			, fichePSORTNameAlias (prefix)
			) ; 
	}
      
      /* https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&name=Arabidopsis+thaliana&lvl=0&srchmode=1 */
    if ((ptr = ficheTAXClosestAncestorStatement (gmp)))
      {
	char buf2[1000] ;
	char *cp ;
	
	cp = ptr ;
	while (*cp == ' ') cp++ ;
	sprintf (buf2
		 , "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&lvl=2&srchmode=1&name=%s"
		 , cp) ;
	vtxtPrintf (blkp, ". ") ;
	
	if (gmp->markup && gmp->mrna)
	  vtxtPrintf (blkp, "<a href=\"javascript:openAnchor(0,'BlastP')\">BlastP results</a>", ac_name (gmp->mrna)) ; 
	else
	  vtxtPrintf (blkp, "BlastP results") ;
	vtxtPrint (blkp, " identify related proteins (threshold .001) in ") ;
	/*fprintf (stderr, "ERROR \"%s\"", vtxtPtr (bfr)) ; fflush (0) ; */
	cp = ptr ;
	if (!strcasecmp (ptr, "teleostomi"))
	  cp = "vertebrates" ;
	if (1 ||  /* never link from the summary */
	    strstr (ptr, " and ") || !strcasecmp (ptr, "sushi")
	    ) /* composite name of AceView invention */
	  vtxtPrint (blkp, cp)  ;
	else
	  gmpURL (blkp, gmp, buf2, cp) ;
	messfree (ptr) ;
      }
    }

  /* uorf */
  {
    AC_KEYSET ks = ac_objquery_keyset (gmp->mrna, "COUNT {>product; good_product;} > 0 && COUNT {>product; uorf_candidate && !good_product && !best_product;} > 0", h) ;
    if (ks && ac_keyset_count (ks) > 0) 
      {
	vtxtDot (blkp) ;
	vtxtPrintf (blkp, "Efficacy of translation of the protein is likely impaired, in a controlled"
		    " fashion, by the presence of a shorter translated product  (") ;
	vtxtPrint (blkp, "<a href=\"javascript:openAnchor ('fmol','Proteins')\">uORF</a>") ;
	vtxtPrintf (blkp, ") initiating at an AUG upstream of the main open reading frame") ;
      }
  }

  ac_free (h) ;

  return 1 ; 
} /* ficheMRNAOverviewParagraphContent */

/***************************************************************************************/
/***************************************************************************************/
/* tye: 1=Disease, 2=Pathways, 3=Go_b, 4=Process, , 5=Pfam, Psort=6, 7=protein interactions 8=gene interactio */
typedef struct itrsStruct { KEY gene ; AC_OBJ Gene ; int nn, type ; double cc ; } ITRC ;
#define FDEBUG FALSE
#define  Ndisease      1000
#define  Nbiocarta     4000
#define  Nkegg         4000

#define  Nexpression_tissue       3000
#define  Nexpression_pathology    3000
#define  Nexpression_external     3000

#define  Ngob_ace      2000
#define  Ngob_pfam     1000
#define  Ngob_iea      0000

#define  Ngom_ace      2000
#define  Ngom_pfam     1000
#define  Ngom_iea      0000

#define  Ngoc_ace      2000
#define  Ngoc_psort    1000
#define  Ngoc_pfam     1000
#define  Ngoc_iea      0000

#define  Npfam         1000
#define  Npsort        1000
#define  NPinteract1    6000  /* gene->with_protein */
#define  NGinteract1    3000  /* gene->with_gene */
#define  Ninteract2    1000 /* do not contribute to friends */

#define Dtype  0x01
#define Wtype  0x02
#define Gtype  0x04
#define Ltype  0x08
#define Mtype  0x10
#define IPtype  0x20
#define IGtype  0x40
#define Ttype  0x80
#define Etype  0x100

static int itrcOrder (const void *a, const void *b)
{
  const ITRC *aa = (const ITRC *)a , *bb = (const ITRC *)b ;
  int nn = aa->nn - bb->nn ;
  double cc = aa->cc - bb->cc ;
  if (cc > 0) return -1 ; /* large number of interactions first */
  if (cc < 0) return  1 ; /* large number of interactions first */
  if (nn) return -nn ; /* large number of interactions first */
  return aa->gene < bb->gene ? -1 : 1 ;
}

static int itrcGeneOrder (const void *a, const void *b)
{
  const ITRC *aa = (const ITRC *)a , *bb = (const ITRC *)b ;
  double cc = aa->cc - bb->cc ;
  int nn = aa->gene - bb->nn ;
  if (aa->gene < bb->gene) return -1 ;
  if (aa->gene > bb->gene) return 1 ;
  if (cc > 0) return -1 ; /* large number of interactions first */
  if (cc < 0) return  1 ; /* large number of interactions first */
  if (nn) return -nn ; /* large number of interactions first */
  return aa->type - bb->type ;
}

/***************/

static void ficheNewGeneDiseasePathwaysBioProcessTableDisease (AC_TABLE tbl1, GMP *gmp, int justOmim, vTXT bfrPap, int *nGwithPapp, int *nnGwithPapp, KEYSET ksGPap, BOOL *hasVotep, Array aa, AC_HANDLE haa)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *ccp ;
  char linkBuf[2000] ;
  char *ccr, *ccg ;
  /*   const char *cgid ; */
  int ir, jr, jr1, row, ng ;
  AC_OBJ obj, obj1, mesh1 ;
  AC_TABLE tbl ;
  AC_ITER iter, iter1 ;
  AC_KEYSET ks, ksOmim = 0, ksKegg, ksGad, ksPaper, ksGene, ksGeneU ;
  vTXT bfr = vtxtHandleCreate (h) ;  
  vTXT bfr1 = vtxtHandleCreate (h) ;  
  enum { COL_TYPE=0, COL_MESH, COL_EVIDENCE, COL_SOURCE, COL_GENES, COL_VOTE, COL_LAST} ; 

  if (gmp->markup) vtxtMarkup (bfr) ;
  if (gmp->markup) vtxtMarkup (bfr1) ;

  ccg = ac_protect (ac_name(gmp->gene), h) ;
  ks = ficheNewGeneDiseaseKs (gmp->gene, h) ;
  iter = ac_keyset_iter (ks, TRUE, h) ;
 
  obj = 0 ; 
  while (ac_free (obj), (obj = ac_iter_obj (iter)))
    {
      ccr = ac_protect (ac_name(obj), h) ;
      /* baddies */
      if (strstr (ccr, "Genetic Predisposition to Disease"))
	continue ;
      
      if (gmp->Spc == HUMAN)
	ksOmim = ac_objquery_keyset (gmp->gene, messprintf ("{>Extern ; OMIM_disease;} $| {>Disease; >OMIM_alias}; gene = %s ;  COUNT {{>OMIM_title; IS zero} $| {>OMIM_alias} ; {!Hidden_alias_of} $| {>Hidden_alias_of} ; {! alias_of} $| {>alias_of}; {IS *} $| {>meshkey;>parent;>mesh}; {IS *} $| {>meshkey;>parent;>mesh} ; {IS *} $| {>meshkey;>parent;>mesh} ; IS %s} > 0", ccg, ccr) , h) ; 
      else if (gmp->Spc == MOUSE)
	ksOmim = ac_objquery_keyset (gmp->gene, messprintf ("{>Extern ; MGI;} ; gene = %s && Disease = %s ", ccg, ccr) , h) ;
      else if (gmp->Spc == RAT)
	ksOmim = ac_objquery_keyset (gmp->gene, messprintf ("{>Extern ; RATGI;} ; gene = %s && Disease = %s ", ccg, ccr) , h) ;
      /* display the OMIM supported didease first */
      if (ksOmim && ac_keyset_count (ksOmim))
	{ if (justOmim>0) continue ; }
      else
	{ if (justOmim==0) continue ; }

      ksKegg = ac_objquery_keyset (gmp->gene, messprintf (" >Extern ; KEGG_disease ; COUNT {>pathway ; {!Hidden_alias_of} $| {>Hidden_alias_of} ; {! alias_of} $| {>alias_of}; {IS *} $| {>meshkey;>parent;>mesh} ; {IS *} $| {>meshkey;>parent;>mesh} ; {IS *} $| {>meshkey;>parent;>mesh}; IS %s} > 0", ccr) , h) ;
      
     if (ksKegg && ac_keyset_count (ksKegg))
	{ if (justOmim>1) continue ; }
      else
	{ if (justOmim==1) continue ; }
       
      ksGad = ac_objquery_keyset (gmp->gene, messprintf (">Extern ; GAD && NOT AntiGad ; COUNT {>GAD_TITLE NOT AntiGad ; {!Hidden_alias_of} $| {>Hidden_alias_of} ; {! alias_of} $| {>alias_of}; {IS *} $| {>meshkey;>parent;>mesh}; {IS *} $| {>meshkey;>parent;>mesh} ; {IS *} $| {>meshkey;>parent;>mesh} ; IS %s} > 0", ccr) , h) ;
 
      if (ksGad && ac_keyset_count (ksGad))
	{ if (justOmim>2) continue ; }
      else
	{ if (justOmim==2) continue ; }
      
      if (justOmim==0 && gmp->Spc == MOUSE)
	{ /* very special way to store the biblio */
	  const char *err_message = 0 ;
	  AC_TABLE tbl2 = ac_bql_table (gmp->db, hprintf(h, "select p from mgi in @active:1, d in mgi->disease where d like %s, p in d[1]", ccr), ksOmim, 0, &err_message, h) ;
	  ksPaper = ac_table_keyset (gmp->db, tbl2, 0, h) ;
	}
      else
	ksPaper = ac_objquery_keyset (obj, messprintf ("{IS *} $| {>alias_of} ; {IS *} $| {>meshkey;>child;>mesh};  {IS *} $| {>Extern NOT antigad} $| {>OMIM_title; IS zero} $| {>GAD_title NOT AntiGad} $| {>KEGG_pathway KEGG_disease}  ;  >Reference  IS pm* && COUNT gene < 3 ; gene = %s || COUNT {>extern NOT AntiGad;gene=%s} > 0", ccg, ccg), h) ;
      
      row = tbl1->rows ;
      vtxtClear (bfr) ;       vtxtClear (bfr1) ;
      
      if (!row)
	ac_table_insert_text (tbl1, row, COL_TYPE, "Disease") ;
      ac_table_insert_text (tbl1, row, COL_LAST,  "1") ;
      ccp = ac_tag_printable (obj, "Meshkey", 0) ;
      if (ccp) 
	{
	  sprintf (linkBuf, MESH_LINK, ac_name (obj)) ;
	  gmpURL (bfr, gmp, linkBuf, ac_name (obj)) ;
	}
      else
	vtxtPrint (bfr, ac_name (obj)) ;
      
      ksGeneU = ac_objquery_keyset (obj, "{IS *} $| {>hidden_alias} $| {>hidden_alias_of} $| {>alias} $| {>alias_of} ; {IS *} $| {>reference ; COUNT Gene < 3} $| {>Extern NOT antigad} $| {>OMIM_title; IS zero} $| {>GAD_title} $| {>KEGG_pathway KEGG_disease}  ; >Gene", h) ;
      ng = ksGeneU ? ac_keyset_count (ksGeneU) : 0 ;
      if (ng)
	{
	  gmpObjLink (bfr1, gmp, obj
		      , messprintf ("%s gene%s", isOne (ng), _multi(ng))) ;
	}

      iter1 = ac_objquery_iter (gmp->gene, messprintf("{>disease} $| {>Mesh} $| {>Reference ;  mesh && COUNT gene < 3 ; >Mesh meshkey } $| {>Extern  NOT AntiGad;>Disease} $| {>Extern ; KEGG_disease ; >pathway } $| {>Extern ; GAD && NOT AntiGad ; >GAD_TITLE NOT AntiGad} $| {>Extern ; OMIM_disease; >OMIM_title; IS zero} ; ! IS %s ;  COUNT {{!Hidden_alias_of} $| {>Hidden_alias_of} ; {meshkey || ! alias_of} $| {>alias_of} ;  {IS *} SETOR {>meshkey;>parent;>mesh}; {IS *} $| {>meshkey;>parent;>mesh}; {IS *} $| {>meshkey;>parent;>mesh} ; IS %s} > 0", ccr, ccr) , h) ;
      mesh1 = 0 ; jr = jr1 = 0 ;
      while (ac_free (mesh1), iter1 && (mesh1 = ac_iter_obj (iter1)))
	{
	  vtxtPrint (bfr, jr++ ? "; " : " (including: ") ;
	  ccp = ac_tag_printable (mesh1, "Meshkey", 0) ;
	  if (ccp) 
	    {
	      sprintf (linkBuf, MESH_LINK, ac_name (mesh1)) ;
	      gmpURL (bfr, gmp, linkBuf, ac_name (mesh1)) ;
	      
	      {
		ksGene = ac_objquery_keyset (mesh1, "{IS *} $| {>hidden_alias} $| {>hidden_alias_of} $| {>alias} $| {>alias_of} ; {IS *} $| {>reference ; COUNT Gene < 3} $| {>Extern NOT AntiGad} $| {>OMIM_title; IS zero} $| {>GAD_title NOT AntiGad} $| {>KEGG_pathway KEGG_disease}  ; >Gene", h) ;
		ng = ksGene ? ac_keyset_count (ksGene) : 0 ;
		if (ng)
		  {
		    if (0)
		      {
			if (! ksGeneU)
			  ksGeneU = ac_new_keyset (gmp->db, h) ;
			ac_keyset_or (ksGeneU, ksGene) ;
		      }
		    if (!jr1++) vtxtPrint (bfr1, " (and also: ") ;
		    else vtxtPrint (bfr1, "; ") ;
		    gmpObjLink (bfr1, gmp, mesh1
				, messprintf ("%s gene%s", isOne (ng), _multi(ng))) ;
		  }

		if (ksGene)
		  {
		    AC_TABLE tblU = ac_keyset_table (ksGene, 0, -1, 0, h) ;
		    ITRC *itrc ;
		    int iU, nItrc1 = arrayMax (aa), nItrc2 ;
		    KEY kU, k0 = ac_obj_key (gmp->gene) ;
		    double cc = 0 ;
		    iU = tblU->rows ;
		    if (iU < 2) iU = 2 ;
		    if (iU < 10000) cc = - log(((double)iU)/10000.0) ;
		    cc *= Ndisease ;
		    
		    for (iU = 0; iU < tblU->rows ; iU++)
		      {		
			kU = ac_table_key (tblU, iU, 0, 0) ;
			for (nItrc2 = 0 ; nItrc2 < nItrc1 ; nItrc2++)
			  {
			    itrc = arrp (aa, nItrc2, ITRC) ;
			    if (itrc->gene == kU)
			      { itrc->nn++ ; itrc->cc += cc ; itrc->type |= Dtype ; break ; }
			  }
			if (kU != k0 && nItrc2 >= nItrc1)
			  {
			    itrc = arrayp (aa, nItrc1++, ITRC) ;
			    itrc->gene = kU ;
			    itrc->Gene = ac_table_obj (tblU, iU, 0, haa) ; 
			    itrc->nn = 1 ;
			    itrc->cc = cc ;
			    itrc->type |= Dtype ;
			  }
		      }
		    ac_free (tblU) ;
		  }


		ac_free (ksGene) ;
	      }
	    }
	  else
	    vtxtPrint (bfr, ac_name (mesh1)) ;
	}
      if (jr)  vtxtPrint (bfr, ").") ; 
      if (jr1) vtxtPrint (bfr1, ").") ; 
      if (ksGeneU && ac_keyset_count (ksGeneU))
	{
	  AC_TABLE tblU = ac_keyset_table (ksGeneU, 0, -1, 0, h) ;
	  ITRC *itrc ;
	  int iU, nItrc1 = arrayMax (aa), nItrc2 ;
	  KEY kU, k0 = ac_obj_key (gmp->gene) ;
	  double cc = 0 ;
	  iU = tblU->rows ;
	  if (iU < 2) iU = 2 ;
	  if (iU < 10000) cc = - log(((double)iU)/10000.0) ;
	  cc *= Ndisease ;

	  for (iU = 0; iU < tblU->rows ; iU++)
	    {		
	      kU = ac_table_key (tblU, iU, 0, 0) ;
	      for (nItrc2 = 0 ; nItrc2 < nItrc1 ; nItrc2++)
		{
		  itrc = arrp (aa, nItrc2, ITRC) ;
		  if (itrc->gene == kU)
		    { itrc->nn++ ; itrc->cc += cc ; itrc->type |= Dtype ; break ; }
		}
	      if (kU != k0 && nItrc2 >= nItrc1)
		{
		  itrc = arrayp (aa, nItrc1++, ITRC) ;
		  itrc->gene = kU ;
		  itrc->Gene = ac_table_obj (tblU, iU, 0, haa) ; 
		  itrc->nn = 1 ;
		  itrc->cc = cc ;
		  itrc->type |= Dtype ;
		}
	    }
	  ac_free (tblU) ;
	  ac_free (ksGeneU) ;
	}
      if (ksGad && ac_keyset_count (ksGad))
	{
	  AC_TABLE ksGads = ac_keyset_table (ksGad, 0, -1, 0, h) ;
	  AC_TABLE ksGadTbl ;
	  AC_OBJ myGad ;
	  int iGad, jGad, jjGad = 0, j2Gad = 0 ;
	  for (iGad = 0 ; ksGads && iGad < ksGads->rows ; iGad++)
	    {
	      myGad = ac_table_obj (ksGads, iGad, 0, h) ;
              if (ac_has_tag (myGad, "AntiGad"))
		continue ;
	      ksGadTbl = myGad ? ac_tag_table (myGad, "GAD_alias", h) : 0 ;
	      for (jGad = 0 ; ksGadTbl && jGad < 10 && jGad < ksGadTbl->rows ; jGad++)
		{
		  if (!jjGad++) 
		    { 
		      vtxtBreak (bfr) ; 
		      vtxtPrint (bfr, "<span class='explain2'>") ;
		      if ((ksPaper && ac_keyset_count (ksPaper))  || (ksOmim && ac_keyset_count (ksOmim))) 
			{ vtxtPrintf (bfr, "e.g. from GAD ") ; vtxtBreak (bfr) ; }
		    }
		  else 
		    {
		      if (j2Gad) 
			vtxtComma (bfr) ;
		      else
			vtxtBreak (bfr) ;
		    }
		  j2Gad++ ;
		  vtxtPrintf (bfr, "%s", ac_table_printable (ksGadTbl, jGad, 0, "")) ;
		}
	      ksGadTbl = myGad ? ac_tag_table (myGad, "Properties", h) : 0 ;
	      for (jGad = 0 ; ksGadTbl && jGad < 1 && jGad < ksGadTbl->rows ; jGad++)
		{
		  if (!jjGad++) 
		    { 
		      vtxtBreak (bfr) ; 
		      vtxtPrint (bfr, "<span class='explain2'>") ;
		      if ((ksPaper && ac_keyset_count (ksPaper))  || (ksOmim && ac_keyset_count (ksOmim))) 
			{ vtxtPrintf (bfr, "e.g. from GAD ") ; vtxtBreak (bfr) ; }
		    }
		  else { if (!j2Gad) vtxtBreak (bfr) ; else vtxtPrint (bfr, ": ") ; j2Gad = 0 ; }
		  vtxtPrintf (bfr, "%s", ac_table_printable (ksGadTbl, jGad, 0, "")) ;
		}
	    }
	  if (jjGad)
	    vtxtPrint (bfr, "</span>") ;
	}
      ac_table_insert_text (tbl1, row, COL_MESH, vtxtPtr (bfr)) ;
      if (vtxtPtr (bfr1))
	ac_table_insert_text (tbl1, row, COL_GENES, vtxtPtr (bfr1)) ;
    
      vtxtClear (bfr) ; jr = 0 ;
      jr = 0 ; 
      if (gmp->Spc == HUMAN && ksOmim && ac_keyset_count (ksOmim))
	{
	  tbl = ac_keyset_table (ksOmim, 0, -1, 0, h) ;
	  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	    {
	      if (jr++) vtxtBreak (bfr) ;
	      obj1 = ac_table_obj (tbl, ir, 0, h) ;
	      sprintf (linkBuf, OMIM_LINK, ac_name (obj1) + 5) ; 
	      gmpURL (bfr, gmp, linkBuf, ac_name (obj1)) ;
	    }
	}
      if (ksKegg && ac_keyset_count (ksKegg))
	{
	  int ok = 1 ;
	  if (jr++)  vtxtBreak (bfr1) ;
	  tbl = ac_keyset_table (ksKegg, 0, -1, 0, h) ;
	  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	    {  
	      const char *cgid = ac_tag_printable (gmp->gene, "GeneId", "") ;
	      if (jr++) vtxtBreak (bfr) ;
	      obj1 = ac_table_obj (tbl, ir, 0, h) ;
	      if (gmp->Spc == HUMAN)
		sprintf (linkBuf, KEGG_LINK_HUMAN, ac_name (obj1) + 5, cgid) ; 
	      else if (gmp->Spc == MOUSE)
		sprintf (linkBuf, KEGG_LINK_MOUSE, ac_name (obj1) + 5, cgid) ; 
	      else if (gmp->Spc == RAT)
		sprintf (linkBuf, KEGG_LINK_RAT, ac_name (obj1) + 5, cgid) ; 
	      else if (gmp->Spc == ARA)
		sprintf (linkBuf, KEGG_LINK_ARA, ac_name (obj1) + 5, cgid) ; 
	      else
		{ ok = 0; vtxtPrintf (bfr, "KEGG") ; }
	      if (ok) 
		gmpURL (bfr, gmp, linkBuf, ac_name (obj1)) ;
	    }
	}
      /* vanished 2016
	if (ksGad && ac_keyset_count (ksGad))
	{
	  if (jr++) vtxtBreak (bfr) ;
	  tbl = ac_keyset_table (ksGad, 0, 1, 0, h) ;
	  obj1 = ac_table_obj (tbl, 0, 0, h) ;
	  cgid =  ac_tag_printable (obj1, "GeneId", "") ;
	  sprintf (linkBuf, GAD_LINK, cgid) ;
	  gmpURL (bfr, gmp, linkBuf, "GAD") ;
	}
      */
      if (ksPaper && ac_keyset_count (ksPaper))
	{
	  int ir1 ;
	  (*nGwithPapp)++  ;
	  vtxtClear (bfr1) ; 
	  vtxtPrint (bfr1, PUBMED_LINK) ; /* no %s included */
	  tbl = ac_keyset_table (ksPaper, 0, -1, 0, h) ;
	  for (ir = tbl->rows - 1, ir1 = 0 ; ir >= 0 ; ir--)
	    {
	      const char *cp = ac_table_printable (tbl, ir, 0, 0) ;
	      KEY kPap = ac_table_key (tbl, ir, 0, 0) ;
	      if (cp && !strncasecmp (cp, "pmp", 2))
		{
		  if (ir1++) vtxtPrint (bfr1, ",") ;
		  vtxtPrint (bfr1, cp+2) ; 
		  if (keySetInsert (ksGPap, kPap))
		    {
		      if (keySetMax (ksGPap))
			vtxtPrint (bfrPap, ",") ;
		      vtxtPrint (bfrPap, cp + 2) ;
		      (*nnGwithPapp)++ ;
		    }
		}
	    }
	  if (ir1)
	    {
	      if (jr++) vtxtBreak (bfr) ;
	      gmpURL (bfr, gmp, vtxtPtr (bfr1)
		      , messprintf ("%d article%s", tbl->rows, _multi(tbl->rows))) ;
	    }
	} 
      if (jr)
	ac_table_insert_text (tbl1, row, COL_EVIDENCE, vtxtPtr (bfr)) ;
      
      vtxtClear (bfr) ;
      
      jr = 0 ;
      if (gmp->Spc == HUMAN && ksOmim && ac_keyset_count (ksOmim))
	{
	  if (jr++) vtxtBreak (bfr) ;
	  vtxtPrint (bfr, "OMIM") ;
	}
      
      if (gmp->Spc == MOUSE && ksOmim && ac_keyset_count (ksOmim))
	{
	  if (jr++) vtxtBreak (bfr) ;
	  tbl = ac_keyset_table (ksOmim, 0, 1, 0, h) ;
	  obj1 = ac_table_obj (tbl, 0, 0, h) ; 
	  ccp = ac_name (obj1) ;
	  if (!ccp || strncmp (ccp, "MGI_", 4))
	      continue ;
	  sprintf (linkBuf, MGI_ALLELE_LINK, ccp + 4) ;
	  gmpURL (bfr, gmp, linkBuf, ac_name(obj1)) ;
	}

      if (ksKegg && ac_keyset_count (ksKegg))
	{
	  if (jr++) vtxtBreak (bfr) ;
	  tbl = ac_keyset_table (ksKegg, 0, 1, 0, h) ;
	  obj1 = ac_table_obj (tbl, 0, 0, h) ; 
	  {
	    const char *cgid = ac_tag_printable (gmp->gene, "GeneId", "") ;
	    if (gmp->Spc == HUMAN)
	      sprintf (linkBuf, KEGG_LINK_HUMAN, ac_name (obj1) + 5, cgid) ; 
	    else if (gmp->Spc == MOUSE)
	      sprintf (linkBuf, KEGG_LINK_MOUSE, ac_name (obj1) + 5, cgid) ; 
	    else if (gmp->Spc == RAT)
	      sprintf (linkBuf, KEGG_LINK_RAT, ac_name (obj1) + 5, cgid) ; 
	    else if (gmp->Spc == ARA)
	      sprintf (linkBuf, KEGG_LINK_ARA, ac_name (obj1) + 5, cgid) ; 
	  }
	  gmpURL (bfr, gmp, linkBuf, "KEGG") ;
	}

      /* vanished 2016
      if (ac_keyset_count (ksGad))
	{
	  if (jr++) vtxtBreak (bfr) ;
	  tbl = ac_keyset_table (ksGad, 0, 1, 0, h) ;
	  obj1 = ac_table_obj (tbl, 0, 0, h) ;
	  cgid =  ac_tag_printable (obj1, "GeneId", "") ;
	  sprintf (linkBuf, GAD_LINK, cgid) ;
	  gmpURL (bfr, gmp, linkBuf, "GAD") ;
	}
      */
      
      ks = ac_objquery_keyset (gmp->gene, messprintf (" >Reference ; COUNT gene < 3 ; COUNT {>Mesh ; {!Hidden_alias_of} $| {>Hidden_alias_of} ; {meshkey || ! alias_of} $| {>alias_of}; {IS *} $| {>meshkey;>parent;>mesh}; {IS *} $| {>meshkey;>parent;>mesh} ; {IS *} $| {>meshkey;>parent;>mesh}; IS %s} > 0", ccr) , h) ;
      if (ac_keyset_count (ks))
	{
	  if (jr++) vtxtBreak (bfr) ;
	  vtxtPrintf (bfr, "PubMed") ;
	}
      
      if (vtxtPtr (bfr))
	ac_table_insert_text (tbl1, row, COL_SOURCE, vtxtPtr (bfr)) ;

      if (0 && row < tbl1->rows && ! (ksOmim && ac_keyset_count (ksOmim)))
	{
	  vtxtClear (bfr) ;
	  if (hasVotep) *hasVotep = TRUE ;
	  gmpURL (bfr, gmp, messprintf ("mailto:mieg@ncbi.nlm.nih.gov?subject=vote_gene%s", ac_name(gmp->gene)), "<font color='red'><i>Vote</i></font>") ;
	  ac_table_insert_text (tbl1, row, COL_VOTE, vtxtPtr (bfr)) ;
	}
    }
  ac_free (obj) ;
  ac_free (h) ;
  return ;
} /* ficheNewGeneDiseasePathwaysBioProcessTableDisease */

/***************/
/* worm case, export Locus_phenotype, RNAi etc */
static void ficheNewGeneDiseasePathwaysBioProcessTablePhenotype (AC_TABLE tbl1, GMP *gmp, int justOmim, vTXT bfrPap, int *nGwithPapp, int *nnGwithPapp, KEYSET ksGPap, BOOL *hasVotep, Array aa, AC_HANDLE haa)
{
  AC_HANDLE h = ac_new_handle () ;
  int row, ng ;
  AC_OBJ obj ;
  AC_TABLE tbl ;
  AC_KEYSET ksGeneU ;
  vTXT bfr = vtxtHandleCreate (h) ;  
  vTXT bfr1 = vtxtHandleCreate (h) ;  
  enum { COL_TYPE=0, COL_MESH, COL_EVIDENCE, COL_SOURCE, COL_GENES, COL_VOTE, COL_LAST} ; 

  if (gmp->markup) vtxtMarkup (bfr) ;
  if (gmp->markup) vtxtMarkup (bfr1) ;

  tbl = ac_tag_table (gmp->gene, "Locus_description", h) ;
  for (row = 0 ; tbl && row < tbl->rows ; row++)
    {
      vtxtClear (bfr) ;       vtxtClear (bfr1) ;
      if (!row)
	ac_table_insert_text (tbl1, row, COL_TYPE, "Phenotype") ;
      ac_table_insert_text (tbl1, row, COL_LAST,  "1") ;
      obj = ac_table_obj (tbl, row, 0, h) ;
      vtxtPrint (bfr, ac_name (obj)) ;
      ac_table_insert_text (tbl1, row, COL_MESH, vtxtPtr (bfr)) ;

      ksGeneU = ac_objquery_keyset (obj, ">Gene", h) ;
      ng = ksGeneU ? ac_keyset_count (ksGeneU) : 0 ;
      if (ng)
	{
	  gmpObjLink (bfr1, gmp, obj
		      , messprintf ("%s gene%s", isOne (ng), _multi(ng))) ;
	}
      if (vtxtPtr (bfr1))
	ac_table_insert_text (tbl1, row, COL_GENES, vtxtPtr (bfr1)) ;
      ac_table_insert_text (tbl1, row, COL_EVIDENCE, "Manual") ;
      ac_table_insert_text (tbl1, row, COL_SOURCE, "AceView") ;
      if (ksGeneU && ac_keyset_count (ksGeneU))
	{
	  AC_TABLE tblU = ac_keyset_table (ksGeneU, 0, -1, 0, h) ;
	  ITRC *itrc ;
	  int iU, nItrc1 = arrayMax (aa), nItrc2 ;
	  KEY kU, k0 = ac_obj_key (gmp->gene) ;
	  double cc = 0 ;
	  iU = tblU->rows ;
	  if (iU < 2) iU = 2 ;
	  if (iU < 10000) cc = - log(((double)iU)/10000.0) ;
	  cc *= Ndisease ;

	  for (iU = 0; iU < tblU->rows ; iU++)
	    {		
	      kU = ac_table_key (tblU, iU, 0, 0) ;
	      for (nItrc2 = 0 ; nItrc2 < nItrc1 ; nItrc2++)
		{
		  itrc = arrp (aa, nItrc2, ITRC) ;
		  if (itrc->gene == kU)
		    { itrc->nn++ ; itrc->cc += cc ; itrc->type |= Dtype ; break ; }
		}
	      if (kU != k0 && nItrc2 >= nItrc1)
		{
		  itrc = arrayp (aa, nItrc1++, ITRC) ;
		  itrc->gene = kU ;
		  itrc->Gene = ac_table_obj (tblU, iU, 0, haa) ; 
		  itrc->nn = 1 ;
		  itrc->cc = cc ;
		  itrc->type |= Dtype ;
		}
	    }
	  ac_free (tblU) ;
	  ac_free (ksGeneU) ;
	}

    }
  ac_free (h) ;
  return ;
} /* ficheNewGeneDiseasePathwaysBioProcessTablePhenotype */

/***************************************************************************************/

static void ficheNewGeneDiseaseKeggPathwaysTableNew (AC_TABLE tbl1, GMP *gmp, Array aa)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *ccp ;
  char linkBuf[2000], *ccr ;
  int ir, jr, iType, row, pass ;
  AC_OBJ obj, obj1, gene ;
  AC_TABLE tbl ;
  AC_ITER iter = 0 ;
  AC_KEYSET ksKegg ;
  vTXT bfr = vtxtHandleCreate (h) ;  
  vTXT bfr1 = vtxtHandleCreate (h) ;  
  enum { COL_TYPE=0, COL_MESH, COL_EVIDENCE, COL_SOURCE, COL_GENES, COL_VOTE, COL_LAST} ; 

  if (gmp->markup) vtxtMarkup (bfr) ;   
  if (gmp->markup) vtxtMarkup (bfr1) ;   

  /* KEGG + BIOCARTA */
  iType = 0 ;
  for (pass = 0 ; pass < 2 ; ac_free (iter), pass++)
    {
      switch (pass)
	{
	case 0: /* KEGG */
	  iter = ac_objquery_iter (gmp->gene, ">Extern ; KEGG && ! KEGG_disease; COUNT {>pathway;  ; {!Hidden_alias_of} $| {>Hidden_alias_of} ; {meshkey || ! alias_of} $| {>alias_of} ; meshkey} = 0", h) ;
	  break ;
	case 1: /* BioCarta */
	   iter = ac_objquery_iter (gmp->gene, ">Extern ; BIOCARTA", h) ;
	   break ;
	}
      
      obj = 0 ; vtxtClear (bfr) ; 
      while (ac_free (obj), (obj = ac_iter_obj (iter)))
	{

	  row = tbl1->rows ;
	  if (!iType++)
	    ac_table_insert_text (tbl1, row, COL_TYPE, "<a href='#' name='tg_Pathway'></a>Pathway") ;
	  ac_table_insert_text (tbl1, row, COL_LAST,  "2") ;

	  if (pass==1)
	    {	      
	      vtxtClear (bfr) ;	 
	      
	      sprintf (linkBuf, BIOCARTA_LINK, ac_name (obj)) ;
	      gmpURL (bfr, gmp, linkBuf, ac_name (obj)) ; 
	      ac_table_insert_text (tbl1, row, COL_EVIDENCE, vtxtPtr (bfr)) ;
	      vtxtClear (bfr) ;	 
	      gmpURL (bfr, gmp, linkBuf, "BioCarta") ;
	      ac_table_insert_text (tbl1, row, COL_SOURCE, vtxtPtr (bfr)) ;

	      ksKegg = ac_objquery_keyset (obj, ">geneid ; > gene", h) ;
	      vtxtClear (bfr) ;	 
	      if (1) /* link inside acedb */
		gmpObjLink (bfr, gmp, obj,  messprintf ("%d genes.", ac_keyset_count (ksKegg))) ;
	      else /* link back to biocarta */
		gmpURL (bfr, gmp, linkBuf,  messprintf ("%d genes.", ac_keyset_count (ksKegg))) ;
	      ac_table_insert_text (tbl1, row, COL_GENES, vtxtPtr (bfr)) ;

	      ccp = ac_tag_printable (obj, "Biocarta_title", 0) ;
	      if (ccp) 
		{
		  vtxtClear (bfr) ;
		  gmpURL (bfr, gmp, linkBuf,  ccp) ;
		  ac_table_insert_text (tbl1, row, COL_MESH, vtxtPtr (bfr)) ;
		}

	      tbl = ac_tag_table (obj, "Gene", h) ;
	      if (tbl && tbl->rows)
		{
		  ITRC *itrc ;
		  int nItrc1 = arrayMax (aa), nItrc2 ;
		  KEY kU, k0 = ac_obj_key (gmp->gene) ;
		  double cc = 0 ;
		  ir = tbl->rows ;
		  if (ir < 2) ir = 2 ;
		  if (ir < 10000) cc = - log(((double)ir)/10000.0) ;
		  cc *= Nbiocarta ;
	      
		  for (ir = jr = 0 ; tbl && ir < tbl->rows ; ir++)
		    {
		      gene = ac_table_obj (tbl, ir, 0, h) ; 
		      kU = ac_table_key (tbl, ir, 0, 0) ; 
		      for (nItrc2 = 0 ; nItrc2 < nItrc1 ; nItrc2++)
			{
			  itrc = arrp (aa, nItrc2, ITRC) ;
			  if (itrc->gene == kU)
			    { itrc->nn ++ ; itrc->cc += cc ; itrc->type |= Wtype ; break ; }
			}
		      if (kU != k0 && nItrc2 >= nItrc1)
			{
			  itrc = arrayp (aa, nItrc1++, ITRC) ;
			  itrc->gene = kU ;
			  itrc->Gene = 0 ;
			  itrc->nn = 1 ;
			  itrc->cc = cc ;
			  itrc->type |= Wtype ;
			}
		    }
		}
	      continue ;
	    }

	  ccr = ac_protect (ac_tag_printable (obj, "Pathway", ""), h) ;
	  ksKegg = ac_objquery_keyset (gmp->gene, messprintf (" >Extern ; KEGG  && ! KEGG_disease ; COUNT {>pathway ; {!Hidden_alias_of} $| {>Hidden_alias_of} ; IS %s} > 0", ccr) , h) ;
	  
	  vtxtClear (bfr) ;
	  {
	    const char *cgid = ac_tag_printable (gmp->gene, "GeneId", "") ;
	    if (gmp->Spc == HUMAN)
	      sprintf (linkBuf, KEGG_LINK_HUMAN, ac_name (obj) + 5, cgid) ; 
	    else if (gmp->Spc == MOUSE)
	      sprintf (linkBuf, KEGG_LINK_MOUSE, ac_name (obj) + 5, cgid) ; 
	    else if (gmp->Spc == RAT)
	      sprintf (linkBuf, KEGG_LINK_RAT, ac_name (obj) + 5, cgid) ; 
	    else if (gmp->Spc == ARA)
	      sprintf (linkBuf, KEGG_LINK_ARA, ac_name (obj) + 5, cgid) ; 
	  }
	  ccp = ac_tag_printable (obj, "Pathway", ac_name (obj)) ;
	  gmpURL (bfr, gmp, linkBuf, ccp) ;
	  ac_table_insert_text (tbl1, row, COL_MESH, vtxtPtr (bfr)) ;
	  
	  vtxtClear (bfr) ; jr = 0 ;
	  if (ksKegg && ac_keyset_count (ksKegg))
	    {
	      if (jr++)  vtxtBreak (bfr) ;
	      tbl = ac_keyset_table (ksKegg, 0, -1, 0, h) ;
	      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
		{
		  const char *cgid = ac_tag_printable (gmp->gene, "GeneId", "") ;
		  obj1 = ac_table_obj (tbl, ir, 0, h) ;
		  if (gmp->Spc == HUMAN)
		    sprintf (linkBuf, KEGG_LINK_HUMAN, ac_name (obj1) + 5, cgid) ; 
		  else if (gmp->Spc == MOUSE)
		    sprintf (linkBuf, KEGG_LINK_MOUSE, ac_name (obj1) + 5, cgid) ; 
		  else if (gmp->Spc == RAT)
		    sprintf (linkBuf, KEGG_LINK_RAT, ac_name (obj1) + 5, cgid) ; 
		  else if (gmp->Spc == ARA)
		    sprintf (linkBuf, KEGG_LINK_ARA, ac_name (obj1) + 5, cgid) ; 
		  gmpURL (bfr, gmp, linkBuf, ac_name (obj1)) ;
		}
	    }
	  ac_table_insert_text (tbl1, row, COL_EVIDENCE, vtxtPtr (bfr)) ;
	  
	  vtxtClear (bfr) ;
	  gmpURL (bfr, gmp, linkBuf, "KEGG") ;
	  ac_table_insert_text (tbl1, row, COL_SOURCE, vtxtPtr (bfr)) ;
	  
	  tbl = ac_tag_table (obj, "Gene", h) ;
	  vtxtClear (bfr) ;
	  vtxtClear (bfr1) ;
	  if (tbl && tbl->rows)
	    {
	      BOOL show = TRUE ;
	      ITRC *itrc ;
	      int nItrc1 = arrayMax (aa), nItrc2 ;
	      KEY kU, k0 = ac_obj_key (gmp->gene) ;
	      double cc = 0 ;
	      ir = tbl->rows ;
	      if (ir < 2) ir = 2 ;
	      if (ir < 10000) cc = - log(((double)ir)/10000.0) ;
	      cc *= Nkegg ;
	      
	      vtxtPrintf (bfr1, "%s gene%s", isOne (tbl->rows) , _multi (tbl->rows)) ;
	      for (ir = jr = 0 ; tbl && ir < tbl->rows ; ir++)
		{
		  gene = ac_table_obj (tbl, ir, 0, h) ; 
		  kU = ac_table_key (tbl, ir, 0, 0) ; 
		  for (nItrc2 = 0 ; nItrc2 < nItrc1 ; nItrc2++)
		    {
		      itrc = arrp (aa, nItrc2, ITRC) ;
		      if (itrc->gene == kU)
			{ itrc->nn ++ ; itrc->cc += cc ; itrc->type |= Wtype ; break ; }
		    }
		  if (kU != k0 && nItrc2 >= nItrc1)
		    {
		      itrc = arrayp (aa, nItrc1++, ITRC) ;
		      itrc->gene = kU ;
		      itrc->Gene = 0 ;
		      itrc->nn = 1 ;
		      itrc->cc = cc ;
		      itrc->type |= Wtype ;
		    }
		  
		  if (show &&  tbl->rows < 3)
		    {
		      vtxtPrint (bfr1, jr++ ? ", " : ": ") ;
		      gmpObjLink (bfr1, gmp, gene, ac_name (gene)) ;
		      
		      if (jr > 4 &&  tbl->rows > 8)
			{ vtxtPrint (bfr1, "...") ; show = FALSE ; }
		    }
		}
	    }
	  vtxtDot (bfr1) ;
	  gmpObjLink (bfr, gmp, obj, vtxtPtr (bfr1)) ;
	  ac_table_insert_text (tbl1, row, COL_GENES, vtxtPtr (bfr)) ;
	}
      ac_free (obj) ;
    }

  ac_free (h) ;
  return ;
} /* ficheNewGeneDiseaseKeggPathwaysTableNew */

/***************************************************************************************/

static int ficheNewGeneProcessFunctionLocalizationTableNew (AC_TABLE tbl1, GMP *gmp, vTXT bfrPap, int *nGwithPapp, int *nnGwithPapp, KEYSET ksGPap, int isB, Array aa)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *ccp, *ccp2, *ptr  ;
  char linkBuf[2000], *ccr ;
  const char *goClass, *goAce, *goFollow, *goQuery, *goIea, *goType ;
  int Nitrc = 0, mytype = 0, ir, ir2, jr2, iType, row, pass, row0 = tbl1->rows ;
  AC_OBJ obj, gene ;
  AC_KEYSET ks ;
  AC_TABLE tbl = 0, tbl2 ;
  vTXT bfr = vtxtHandleCreate (h) ;  
  vTXT bfr1 = vtxtHandleCreate (h) ;  
  vTXT bfrSrc = vtxtHandleCreate (h) ;  
  vTXT tBfr = vtxtHandleCreate (h) ;  
  enum { COL_TYPE=0, COL_MESH, COL_EVIDENCE, COL_SOURCE, COL_GENES, COL_VOTE, COL_LAST} ; 

  if (gmp->markup) vtxtMarkup (bfr) ;   
  if (gmp->markup) vtxtMarkup (bfr1) ;   
  goClass = goAce = goFollow = goQuery = goIea = goType = 0 ;

  switch (isB)
    {
    case 0:  /* Process */
      goClass = "Go_b" ;
      goFollow = ">gene_ace" ;
      goAce = "Go_b_ace" ; 
      goQuery = "{>Go_b_pfam} SETMINUS {>Go_b_ace}" ;
      goIea = "{>Go_b_iea} SETMINUS { {>Go_b_pfam} SETOR {>Go_b_ace}}" ;
      goType = "<a href='#' name='tg_Process'></a>Process" ;
      Nitrc = 9999999 ;
      mytype = Gtype ;
      break ;
    case 1: /* function */
      goClass = "Go_m" ; 
      goFollow = "{>gene_ace}" ;
      goAce = "Go_m_ace" ; 
      goQuery = "{>Go_m_pfam} SETMINUS {>Go_m_ace}" ;
      goIea = "{>Go_m_iea} SETMINUS { {>Go_m_pfam} SETOR {>Go_m_ace}}" ;
      goType =  "Function" ;
      Nitrc = 9999999 ;
      mytype = Gtype ;
      break ;
    case 2: /* localization */
      goClass = "Go_c" ;
      goAce = "Go_c_ace" ;
      goFollow = "{CLASS go_c ; >gene_ace} SETOR { CLASS go_c ; >gene_psort}" ;
      goQuery = "{>Go_c_pfam} SETMINUS {>Go_c_ace}" ;
      goIea = "{>Go_c_iea} SETMINUS { {>Go_c_pfam} SETOR {>Go_c_ace}}" ;
      goType =  "Localization" ;
      Nitrc = 9999999 ;
      mytype = Ltype ;
      break ;
    }
  if (gmp->markup) vtxtMarkup (bfr) ;

  /* Go_b_pfam */
  iType = 0 ; vtxtPrint (tBfr, "toto ") ;

  for (pass = 0 ; pass < 4 ; pass++)
    {
      switch (pass)
	{
	case 0:
	  tbl = ac_tag_table (gmp->gene, goAce, h) ;
	  if (isB== 0) Nitrc = Ngob_ace ;
	  if (isB== 1) Nitrc = Ngom_ace ;
	  if (isB== 2) Nitrc = Ngoc_ace ;
	  break ;
	case 1: 
	  if (isB == 0) continue ;
	  else if (isB == 1)
	    {
	      continue ;
	      /* Locus_description is exported in the function  ...Phenotype  below
	       * i leave this here as a template if we need to add here something else
	       */
	      Nitrc = Ngom_ace ;
	      tbl = ac_tag_table (gmp->gene, "Locus_description", h) ;
	    }
	  else if (isB == 2) 
	    {
	      Nitrc = Ngoc_psort ;
	      tbl = ac_tag_table (gmp->gene, "go_c_psort", h) ;
	    }
	  break ;
	case 2: 
	  ks = ac_objquery_keyset (gmp->gene, goQuery, h) ;
	  tbl = ac_keyset_table (ks, 0, -1 , 1, h) ;
	  if (isB== 0) Nitrc = Ngob_pfam ;
	  if (isB== 1) Nitrc = Ngom_pfam ;
	  if (isB== 2) Nitrc = Ngoc_pfam ;
	  break ;
	case 3: 
	  ks = ac_objquery_keyset (gmp->gene, goIea, h) ;
	  tbl = ac_keyset_table (ks, 0, -1 , 1, h) ;
	  if (isB== 0) Nitrc = Ngob_iea ;
	  if (isB== 1) Nitrc = Ngom_iea ;
	  if (isB== 2) Nitrc = Ngoc_iea ;
	  break ;
	}
 
      for (ir=0; tbl && ir < tbl->rows && tbl->cols >= 1; ir++)
	{
	  if (! (ptr = ac_table_printable (tbl, ir, 0, 0)) ||
	      strstr (ptr, "unknown") ||
	      pickMatch (vtxtPtr (tBfr), messprintf("*%s*", gtLowerCleanUp(ptr))))
	    continue ;
	  vtxtPrint (tBfr, ptr) ; /* avoid duplications */

	  row = tbl1->rows ;
	  if (!iType++)
	    ac_table_insert_text (tbl1, row, COL_TYPE,  goType) ;
	  ac_table_insert_text (tbl1, row, COL_LAST,  goAce) ;
	  ccr = gmp->Spc == WORM ? gtLowerCleanUp(ptr) : gtCleanUp(ptr) ;
	  obj =  ac_table_obj (tbl, ir, 0, h) ;

	  vtxtClear (bfr) ;
	  ccp = ac_tag_printable (obj, "GO", 0) ;
	  if (ccp)
	    {
	      char *ccp2, *ccp3 ;
	      ccp2 = strstr(ccp, ",") ;
	      if (ccp2)
		{
		  ccp3 = messprintf ("%s", ccp) ;
		  ccp2 = strstr(ccp3, ",") ;
		  *ccp2 = 0 ;
		  sprintf (linkBuf, GO_LINK, ccp3) ; 
		}
	      else
		sprintf (linkBuf, GO_LINK, ccp) ; 
	      gmpURL (bfr, gmp, linkBuf, ccr) ;
	    }
	  else 
	    {
	      if (strstr (ccr, "extracellular"))
		ccr = "secreted or extracellular" ;
	      vtxtPrint (bfr, ccr) ;
	    }

	  ac_table_insert_text (tbl1, row, COL_MESH, vtxtPtr (bfr)) ;
	  
	  switch (pass)
	    {
	    case 1: 
	      ac_table_insert_text (tbl1, row, COL_EVIDENCE, isB == 2 ? "Psort" : "Manual") ;
	      ac_table_insert_text (tbl1, row, COL_SOURCE,  "AceView") ;
	      break ;
	    case 2: 
	      ac_table_insert_text (tbl1, row, COL_EVIDENCE, "Pfam") ;
	      ac_table_insert_text (tbl1, row, COL_SOURCE,  "AceView") ;
	      break ;
	    case 3:
	      ac_table_insert_text (tbl1, row, COL_SOURCE
				    , "<a href='http://www.geneontology.org/GO.evidence.shtml'>GOA/IEA</a>") ;
	      break ;
	    case 0:  
	      vtxtClear (bfrSrc) ;
	      ac_table_insert_text (tbl1, row, COL_SOURCE,  "<a href='http://www.geneontology.org/GO.evidence.shtml'>GOA/IEA</a>GOA</a>") ;
	      if (gmp->Spc == WORM)
		vtxtPrint (bfrSrc, "Manual/TAS") ;
	      else
		{
		  int jr1 ;
		  int ir1 = 0 ;

		  vtxtClear (bfr1) ; vtxtPrint (bfr1, PUBMED_LINK) ; /* no %s included */
		  const char *ccp11, *ccp10 = ac_table_printable (tbl, ir, 0, 0) ;
		  vtxtClear (bfr) ;
		  /* subloop, because a single go_c_ace may have several types of support: ABCA2 */
		  for (jr1 = ir1 ; jr1 < tbl->rows ; jr1++)
		    {
		      ccp11 = ac_table_printable (tbl, jr1, 0, 0) ;
		      if (!ccp10 || !ccp11 || strcmp (ccp10, ccp11))
			continue ;
		      ccp2 = ac_table_printable (tbl, jr1, 1, 0) ;
		      if (ccp2 &&
			  (!vtxtPtr (bfrSrc) || !strstr (vtxtPtr (bfrSrc), ccp2))
			  )
			{
			  if (!vtxtPtr (bfrSrc))
			    vtxtPrint (bfrSrc, "<a href='http://www.geneontology.org/GO.evidence.shtml'>GOA/") ;
			  else
			    vtxtComma (bfrSrc) ;
			  vtxtPrint (bfrSrc, ccp2) ;
			}
		      ccp2 = ac_table_printable (tbl, jr1, 2, 0) ;
		      if (ccp2)
			{ /* loop on pm broken on | */
			  char *cq, *cp = strnew (ccp2, h) ;
			  AC_OBJ Pap ;
			  KEY kPap ;
			  (*nGwithPapp)++ ;
			  while (cp && *cp)
			    {
			      cq = strstr (cp, "|") ; 
			      if (cq) *cq = 0 ;
			      if (ir1++) vtxtPrint (bfr1, ",") ; 
			      vtxtPrint (bfr1, cp) ;
			      if (!strncasecmp(cp, "pm", 2))
				Pap = ac_get_obj (gmp->db, "Paper", cp, h) ;
			      else
				Pap = ac_get_obj (gmp->db, "Paper", messprintf("pm%s",cp), h) ;
			      if (Pap)
				{
				  kPap = ac_obj_key (Pap) ;
				  if (keySetInsert (ksGPap, kPap))
				    {
				      if (keySetMax (ksGPap)>1)
				      vtxtPrint (bfrPap, ",") ;
				      if (!strncasecmp(cp, "pm", 2)) cp+=2 ;
				      vtxtPrint (bfrPap, cp) ;
				      (*nnGwithPapp)++ ;
				    }
				  ac_free (Pap) ;
				}
			      if (cq) cp = cq + 1 ; 
			      else cp = 0 ;
			    }
			}
		    }
		  if (ir1)
		    {
		      vtxtBreak (bfr) ;
		      gmpURL (bfr, gmp, vtxtPtr (bfr1)
			      , messprintf ("%d article%s", ir1, _multi(ir1))) ;
		    }

		  if (vtxtPtr (bfr))
		    ac_table_insert_text (tbl1, row, COL_EVIDENCE, vtxtPtr (bfr)) ;
		}
	      if (vtxtPtr (bfrSrc))
		{
		  vtxtPrint (bfrSrc, "</a>") ;
		  ac_table_insert_text (tbl1, row, COL_SOURCE,  vtxtPtr (bfrSrc)) ;
		}
	      break ;
	    }
	  vtxtClear (bfr) ; 
	  ks = ac_objquery_keyset (obj, goFollow, h) ;
	  
	  tbl2 = ks ? ac_keyset_table (ks, 0, -1, 0, h) : 0 ;
	  vtxtClear (bfr) ; 
	  if (tbl2 && tbl2->rows)
	    {
	      BOOL show = TRUE ;
	      ITRC *itrc ;
	      int nItrc1 = aa ? arrayMax (aa) : 0, nItrc2 ;
	      KEY kU, k0 = ac_obj_key (gmp->gene) ;
	      char *cp = messprintf ("%s gene%s", isOne (tbl2->rows) , _multi (tbl2->rows)) ; 
	      double cc = 0 ;
	      
	      ir2 = tbl2->rows ;
	      if (ir2 < 2) ir2 = 2 ;
	      if (ir2 < 10000) cc = - log(((double)ir2)/10000.0) ;
	      cc *= Nitrc ;
	      
	      gmpFakeObjLink (bfr, gmp, goClass, obj, cp) ;
	      for (ir2 = jr2 = 0 ; ir2 < tbl2->rows ; ir2++)
		{
		  gene = ac_table_obj (tbl2, ir2, 0, h) ;
		  if (aa) 
		    {
		      kU = ac_table_key (tbl2, ir2, 0, 0) ; 
		      for (nItrc2 = 0 ; nItrc2 < nItrc1 ; nItrc2++)
			{
			  itrc = arrp (aa, nItrc2, ITRC) ;
			  if (itrc->gene == kU)
			    { itrc->nn++ ; itrc->cc += cc ; itrc->type |= mytype ; break ; }
			}
		      if (kU != k0 && nItrc2 >= nItrc1)
			{
			  itrc = arrayp (aa, nItrc1++, ITRC) ;
			  itrc->gene = kU ;
			  itrc->Gene = 0 ;
			  itrc->nn = 1 ;
			  itrc->cc = cc ;
			  itrc->type |= mytype ;
			}
		    }
	      
		  if (show &&  tbl2->rows < 3)
		    {
		      vtxtPrint (bfr, jr2++ ? ", " : ": ") ;
		      gmpObjLink (bfr, gmp, gene, ac_name (gene)) ;
		    
		      if (jr2 > 4 &&  tbl2->rows > 8)
			{ vtxtPrint (bfr, "...") ; show = FALSE ; }
		    }
		}
	    }
	  if (pass < 2 && isB <2 && vtxtPtr (bfr)) vtxtPrintf (bfr, " *") ;
	  ac_table_insert_text (tbl1, row, COL_GENES, vtxtPtr (bfr)) ;
	}
    }

  if (tbl1->rows - row0 >= 1 && isB == 2 &&
      (
       tbl1->rows - row0 > 1 ||
       ac_keyset_count (ac_objquery_keyset (gmp->gene, ">product; good_product && best_product ; > Kantor", h)) > 1
       )
      )
    {
      row = tbl1->rows ; 
      ac_table_insert_text (tbl1, row, COL_LAST, "goAce") ;
      ac_table_insert_text (tbl1, row, COL_MESH, "Different localizations may apply to different protein isoforms.") ;
    }
  ac_free (h) ;
  return tbl1->rows - row0 ;
} /* ficheNewGeneProcessFunctionLocalizationTableNew */

/***************************************************************************************/

static int ficheNewGeneDiseasePathwaysBioProcessTableNew (vTXT blkp, GMP *gmp, BOOL *hasVotep, Array bb)
{
  AC_HANDLE h = ac_new_handle () ;
  Array aa = arrayHandleCreate (256, ITRC, h) ;
  AC_TABLE tbl ; 
  vTXT bfr = vtxtHandleCreate (h) ; 
  vTXT bfrPap = vtxtHandleCreate (h) ; 
  int nGwithPap = 0, nnGwithPap = 0 ; 
  KEYSET ksGPap = keySetHandleCreate (h) ;
  int row = 0, cols[7] = {1,2,3,4, 6,5, 0} ;
  enum { COL_TYPE=0, COL_MESH, COL_EVIDENCE, COL_SOURCE, COL_GENES, COL_VOTE,  COL_LAST} ; 
  const char *colNames[]={ "Type", 
			   "Description", 
			   "Evidence",
			   "Source",
			   "All candidate genes for this disease",
			   "Vote",
			   0} ; 

  /* contruct the table */
  tbl = ac_empty_table (80, 6, h) ; /* tbl->rows is a hint, not a hard limit */

  /* format the table and add the http links */
  nGwithPap = 0, nnGwithPap = 0 ;  
  vtxtClear (bfrPap) ;
  vtxtPrint (bfrPap, PUBMED_LINK) ; /* no %s included */ 

  ficheNewGeneDiseasePathwaysBioProcessTableDisease (tbl, gmp, 0 /* OMIM */, bfrPap, &nGwithPap, &nnGwithPap, ksGPap, hasVotep, aa, h) ;
  ficheNewGeneDiseasePathwaysBioProcessTableDisease (tbl, gmp, 1 /* KEGG */, bfrPap, &nGwithPap, &nnGwithPap, ksGPap, hasVotep, aa, h) ;
  ficheNewGeneDiseasePathwaysBioProcessTableDisease (tbl, gmp, 2 /* GAD */, bfrPap, &nGwithPap, &nnGwithPap, ksGPap, hasVotep, aa, h) ;
  ficheNewGeneDiseasePathwaysBioProcessTableDisease (tbl, gmp, 3, bfrPap, &nGwithPap, &nnGwithPap, ksGPap, hasVotep, aa, h) ;
  ficheNewGeneDiseasePathwaysBioProcessTablePhenotype (tbl, gmp, 4 /* WORM */, bfrPap, &nGwithPap, &nnGwithPap, ksGPap, hasVotep, aa, h) ;

  if (nGwithPap > 1)
    {
      int row = tbl->rows ;
      int ir = nnGwithPap ;
      ac_table_insert_text (tbl, row, COL_MESH, "Cumulated literature") ;
      
      vtxtClear (bfr) ; 
      gmpURL (bfr, gmp, vtxtPtr (bfrPap)
	      , messprintf ("%d article%s", ir, _multi(ir))) ;
      ac_table_insert_text (tbl, row, COL_EVIDENCE, vtxtPtr (bfr)) ;
    }
  
  if (aa && arrayMax (aa))
    {
      int nn, row, ir, ir2 ;
      ITRC *itrc, *itrc2 ;
      int nItrc1 = arrayMax (aa), nItrc2 ;
      double cc, oldcc ;

      arraySort (aa, itrcOrder) ;  
      if (bb)
	for (nItrc2 = 0, ir2 = arrayMax (bb) ; nItrc2 < nItrc1 ; ir2++, nItrc2++)
	  {
	    itrc = arrp (aa, nItrc2, ITRC) ;
	    itrc2 = arrayp (bb, ir2, ITRC) ;
	    *itrc2 = *itrc ; 
	    itrc2->Gene = 0 ;
	  }
	  
      itrc = arrp (aa, 0, ITRC) ;
      nn = itrc->nn ;
      if (nn > 1)
	{
	  row = tbl->rows ;
	  ac_table_insert_text (tbl, row, COL_MESH
			    , "<font color='red'>Genes with most related pattern of diseases</font>") ;
      
	  cc = oldcc = itrc->cc ;
	  for (nItrc2 = ir = 0 ; nItrc2 < nItrc1 ; nItrc2++)
	    {
	      itrc = arrp (aa, nItrc2, ITRC) ;
	      if (2*itrc->cc >= cc) ir++ ;
	    }
	  
	  vtxtClear (bfr) ; 
	  for (nItrc2 = ir2 = 0 ; nItrc2 < nItrc1 ; nItrc2++)
	    {
	      itrc = arrp (aa, nItrc2, ITRC) ;
	      if (2*itrc->cc >= cc || FDEBUG)
		{
		  if (itrc->cc < oldcc)
		    vtxtPrint (bfr, ir2++ ? " | " : "") ;
		  else
		    vtxtPrint (bfr, ir2++ ? ", " : "") ;
		  oldcc = itrc->cc ;
		  gmpObjLink (bfr, gmp, itrc->Gene, ac_name (itrc->Gene)) ;
	      
		  if (ir2 > 8 &&  ir > 12 && ! FDEBUG)
		    { vtxtPrint (bfr, "...") ; break ; }
		}
	    }
	  ac_table_insert_text (tbl, row, COL_GENES, vtxtPtr (bfr)) ;
	}
    }
  
  /* export */
  if ((row = tbl->rows))
    ac_table_display (blkp 
		      , tbl, colNames
		      , cols, 7
		      , 0, 0, 0
		      , 0
		      ) ;

  ac_free (h) ; 
  return row ;
} /* ficheNewGeneDiseasePathwaysBioProcessTableNew */

/***************************************************************************************/
/***************************************************************************************/
/* the same as below product, but this time as a text */

/* the same as above + product, but this time as a table */
#ifdef JUNK
static int ficheNewGeneDiseasePathwaysBioProcessTableJUNK (vTXT blkp, GMP *gmp)
{
  AC_TABLE oTmp, oTmp2, gCogInfo, gProducts ;
  AC_OBJ oOmim, oCog, oGad ;
  int ir, ir1, ir2, jr, jt, jt2, nProducts ;
  int jLocalization = 0 ;
  const char *ccp, *ptr, *goType = "" ;
  char linkBuf[2000] ;
  int  maxRows, maxCols ; 
  Array Tbb = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  char *colNames[]={ "Type", 
		     "Description linking to related genes", "Source",
		     0} ; 
  enum {COL_TYPE=0, COL_VALUE, COL_ORIGIN, COL_LAST} ; 
   
  vTXT tBfr, bfr ; 

  oTmp = ac_tag_table (gmp->gene, "Product", h) ;
  if (!oTmp) return 0 ;
  tBfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (tBfr) ;
  bfr = vtxtCreate () ;
  if (gmp->markup) vtxtMarkup (bfr) ;

  nProducts = oTmp ? oTmp->rows : 0 ;
  maxRows = nProducts + 200 ;
   
  /* initialize the column titles */
  jt = 1 ; jt2 = 0 ;
  Tbb = arrayCreate ((TbbDim+1)* (maxRows+2), int) ;

  vtxtPrintf (tBfr, "toto") ; /* avoid zero */
  if (gmp->markup)
    for (maxCols=0 ; colNames[maxCols] ; maxCols++)
      TBB (0, maxCols) = vtxtPrintf (tBfr, "%s"ooo, colNames[maxCols]) ; 
 

  /* table */
  /* phenotype */
  ir1 = 0 ;
  if ((gmp->Spc == WORM  || gmp->Spc == ARA) &&
      (oTmp = ac_tag_table (gmp->gene, "Locus_Description", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
      {
	if (!ir1++)
	  TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Phenotype"ooo) ;
	ptr = ac_table_printable (oTmp, ir, 0, 0) ;
	if (ptr && *ptr &&
	!pickMatch (vtxtPtr (tBfr), messprintf("*%s*", gtLowerCleanUp (ptr))) &&
	    !strstr(ptr,"ssential gene")) /* essential gene is reported separately */
	  {
	    if (0)
	      TBB (jt++, COL_VALUE) = vtxtPrintf (tBfr, "%s"ooo 
						  , gtLowerCleanUp(ptr)) ;
	    else
	      { /* this works but the link must be activated in aceviewmain */
		TBB (jt++, COL_VALUE) = 
		  gmpObjLink (tBfr, gmp, ac_table_obj (oTmp, ir, 0, h)
			      , gmp->Spc == WORM ? gtLowerCleanUp(ptr) : gtCleanUp(ptr)) ;
		vtxtPrintf (tBfr, ooo) ;
	      }
	  }
      }
  /* OMIM */
  ir1 = 0 ;
  if ((oTmp = ac_tag_table (gmp->gene, "Extern", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
      {
	ccp = ac_table_printable (oTmp, ir, 0, "toto") ;
	if (strncmp (ccp, "OMIM_", 5))
	  continue ;
	oOmim = ac_table_obj (oTmp, ir, 0, h) ;
	if (!ac_has_tag (oOmim, "OMIM_disease"))
	  continue ;
	if (!ir1++)
	  TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Disease"ooo) ;
	sprintf (linkBuf, OMIM_LINK, ac_name (oOmim) + 5) ; 
	TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "OMIM") ;
	vtxtPrintf (tBfr, ooo) ;
	if (
	    (
	     ((ptr = ac_tag_printable (oOmim, "OMIM_title", 0))) &&
	     !pickMatch (vtxtPtr (tBfr), messprintf("*%s*",ptr))
	     ) 
	    )
	  TBB (jt++, COL_VALUE) = gmpURL (tBfr, gmp, linkBuf, ptr) ;
	vtxtPrintf (tBfr, ooo) ;
	ac_free (oOmim) ;
      }

  /* EXTERN */
  ir2 = 0 ;
  if ((oTmp = ac_tag_table (gmp->gene, "Extern", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
      {
	ccp = ac_table_printable (oTmp, ir, 0, "toto") ;
	if (!ccp || strncmp (ccp, "GAD_", 4))
	  continue ;
	ccp = ac_tag_printable (gmp->gene, "GeneId", 0) ;
	if (!ccp)
	  continue ;
	oGad = ac_table_obj (oTmp, ir, 0, h) ; 
	if (ac_has_tag (oGad, "AntiGad"))
	  continue ;
	if (!ac_has_tag (oGad, "Disease_class"))
	  continue ;
	oTmp2 = oGad ? ac_tag_table (oGad, "GAD_title", h) : 0 ;
	for (jr = 0 ; oTmp2 && jr < 1 && jr < oTmp2->rows ; jr++)
	  {
	    ptr = ac_table_printable (oTmp2, jr, 0, 0) ; 
	    if (!ptr || pickMatch (vtxtPtr (tBfr), messprintf("*%s*",ptr)))
	      continue ;
	    if (!ir1++)
	      TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Disease"ooo) ;
	    if (!ir2++)
	      {
		TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "GAD") ;
		vtxtPrintf (tBfr, ooo) ;
	      }
	    else
	      vtxtBreak (tBfr) ;
	    /* sprintf (linkBuf, GAD_LINK, ccp) ;  gone 2016 */
	    TBB (jt, COL_VALUE) = gmpURL (tBfr, gmp, linkBuf, ptr) ;
	  }
	if (ir2) { jt++ ; vtxtPrintf (tBfr, ooo) ; }

	ac_free (oTmp2) ;
	ac_free (oGad) ;
      }

  /* KEGG */
  ir1 = 0 ;
  if ((oTmp = ac_tag_table (gmp->gene, "Extern", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
      {
	ccp = ac_table_printable (oTmp, ir, 0, "toto") ;
	if (strncmp (ccp, "KEGG_", 5))
	  continue ;
	oOmim = ac_table_obj (oTmp, ir, 0, h) ;
	if (!ir1++)
	  TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Pathway"ooo) ;
	{
	  const char *cgid = ac_tag_printable (gmp->gene, "GeneId", "") ;
	  if (gmp->Spc == HUMAN)
	    sprintf (linkBuf, KEGG_LINK_HUMAN, ac_name (oOmim) + 5, cgid) ; 
	  else if (gmp->Spc == MOUSE)
	    sprintf (linkBuf, KEGG_LINK_MOUSE, ac_name (oOmim) + 5, cgid) ; 
	  else if (gmp->Spc == RAT)
	    sprintf (linkBuf, KEGG_LINK_RAT, ac_name (oOmim) + 5, cgid) ; 
	  else if (gmp->Spc == ARA)
	    sprintf (linkBuf, KEGG_LINK_ARA, ac_name (oOmim) + 5, cgid) ; 
	}
	TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "KEGG") ;
	vtxtPrintf (tBfr, ooo) ;
	if (
	    (
	     ((ptr = ac_tag_printable (oOmim, "Pathway", 0))) &&
	     !pickMatch (vtxtPtr (tBfr), messprintf("*%s*",ptr))
	     ) 
	    )
	  TBB (jt++, COL_VALUE) = gmpURL (tBfr, gmp, linkBuf, ptr) ;
	vtxtPrintf (tBfr, ooo) ;
	ac_free (oOmim) ;
      }
  
  /* motif */
  if ((oTmp = ac_tag_table (gmp->gene, "Pfam", h)))
    for (ir=jt2=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
      if ((ptr = ac_table_printable (oTmp, ir, 0, 0)))
	{
	  AC_OBJ oFam = ac_table_obj (oTmp, ir, 0, h) ;
	  const char *txtAccession = ac_tag_printable (oFam, "Accession", 0) ;

	  if (txtAccession)
	    {
	      char *cp = 0, *cq = strstr (txtAccession , ".") ; /* drop pfam version number */
	      if (cq) 
		{
		  cp = strnew (txtAccession, 0) ;
		  cq = strstr (cp, ".") ; /* drop pfam version number in non const char* */
		  *cq = 0 ;
		  txtAccession = cp ;
		}
	      sprintf (linkBuf, "http://pfam.xfam.org/family/%s", txtAccession) ; 
	      TBB (jt, COL_ORIGIN) = gmpURL (tBfr, gmp, linkBuf, "Pfam"ooo) ; 
	      messfree (cp) ;
	    }
	  else
	    TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "PFam"ooo) ;
	  if (!jt2++)
	    TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Motif"ooo) ;
	  
	  ptr = ac_tag_printable (oFam, "Definition", ptr) ;
	  TBB (jt, COL_VALUE) = 
	    gmpObjLink (tBfr, gmp, ac_table_obj (oTmp, ir, 0, h)
			, gtCleanUp(ptr)) ;
	  jt++ ;
	  vtxtPrintf (tBfr, ooo) ;
	}
 
  if ((oTmp = ac_tag_table (gmp->gene, "Psort_domain", h)))
    for (ir=0;ir < oTmp->rows ;ir++)
      {
	AC_OBJ oProduct = ac_table_obj (oTmp, ir, 0, h) ;
	if ((oTmp2 = ac_tag_table (oProduct, "Psort_domain", h))) 
	  for (jr=0;jr < oTmp2->rows ; jr++)
	    if ((ptr = ac_table_printable (oTmp2, jr, 0, 0)))
	      {
		if (strcasecmp (ptr, "Nuclear_localization_domain"))
		  continue ;
		if (strcasecmp (ptr, "Dileucine_domain"))
		  continue ;
		/* http://psort.org */
		sprintf (linkBuf, "http://psort.nibb.ac.jp") ;
		TBB (jt, COL_ORIGIN) = gmpURL (tBfr, gmp, linkBuf, "Psort"ooo) ;
		if (!jt2++)
		  TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Motif"ooo) ;
		
		TBB (jt, COL_VALUE) = 
		  gmpObjLink (tBfr, gmp, ac_table_obj (oTmp2, jr, 0, h)
			      , gtCleanUp(ptr)) ;
		jt++ ;
		vtxtPrintf (tBfr, ooo) ;
	      }
      }
  
  /* function */
  jt2 = 0 ;
  if ((oTmp = ac_tag_table (gmp->gene, "Descriptor", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
      if ((ptr = ac_table_printable (oTmp, ir, 0, 0)) &&
	  ! strstr (ptr, "unknown") &&
	  !pickMatch (vtxtPtr (tBfr), messprintf("*%s*",gtLowerCleanUp(ptr))))
	{
	   if (ir != -1)
	      {
		if (gmp->Spc == HUMAN)
		  TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Entrez"ooo) ;
		else
		  TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Manual"ooo) ;
	      }
	   if (!jt2++)
	    TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Function"ooo) ;
	   goType = ac_table_printable (oTmp, ir, 1, "") ;
	   TBB (jt++, COL_VALUE) = 
	    gmpObjLink (tBfr, gmp, ac_table_obj (oTmp, ir, 0, h)
			, gmp->Spc == WORM ? gtLowerCleanUp(ptr) : messprintf ("[GOA%s%s] %s", goType[0] ? " " : "", goType,  gtCleanUp(ptr))) ;
	  vtxtPrintf (tBfr, ooo) ;
	}
 
  if ((oTmp = ac_tag_table (gmp->gene, "Pathway", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
      if ((ptr = ac_table_printable (oTmp, ir, 0, 0)) &&
	  ! strstr (ptr, "unknown") &&
	  !pickMatch (vtxtPtr (tBfr), messprintf("*%s*",gtLowerCleanUp(ptr))))
	{
	   if (ir != -1)
	      {
		if (gmp->Spc == HUMAN)
		  TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Entrez"ooo) ;
		else
		  TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Manual"ooo) ;
	      }
	   if (!jt2++)
	    TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Function"ooo) ;
	  TBB (jt++, COL_VALUE) = 
	    gmpObjLink (tBfr, gmp, ac_table_obj (oTmp, ir, 0, h)
			, gmp->Spc == WORM ? gtLowerCleanUp(ptr) : gtCleanUp(ptr)) ;
	  vtxtPrintf (tBfr, ooo) ;
	}
 
  if ((oTmp = ac_tag_table (gmp->gene, "Go_b_ace", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
      if ((ptr = ac_table_printable (oTmp, ir, 0, 0)) &&
	  ! strstr (ptr, "unknown") &&
	  !pickMatch (vtxtPtr (tBfr), messprintf("*%s*",gtLowerCleanUp(ptr))))
	{
	  if (ir != -1)
	      {
		if (gmp->Spc == HUMAN)
		  TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Entrez"ooo) ;
		else
		  TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Manual"ooo) ;
	      }
	  if (!jt2++)
	    TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Function"ooo) ;
	  goType = ac_table_printable (oTmp, ir, 1, "") ;
	  TBB (jt++, COL_VALUE) = gmpFakeObjLink (tBfr, gmp, "Go_b", ac_table_obj (oTmp, ir, 0, h)
					      , gmp->Spc == WORM ? gtLowerCleanUp(ptr) : messprintf ("[GOA%s%s] %s", goType[0] ? " " : "", goType,  gtCleanUp(ptr))) ;
	  vtxtPrintf (tBfr, ooo) ;
	}

  if ((oTmp = ac_tag_table (gmp->gene, "Go_b_pfam", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
      if ((ptr = ac_table_printable (oTmp, ir, 0, 0)) &&
	  ! strstr (ptr, "unknown") &&
	  !pickMatch (vtxtPtr (tBfr), messprintf("*%s*",gtLowerCleanUp(ptr)) ))
	{
	  if (ir != -1)
	    TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Pfam"ooo) ;
	  if (!jt2++)
	    TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Function"ooo) ;
	  TBB (jt++, COL_VALUE) = gmpFakeObjLink (tBfr, gmp, "Go_b", ac_table_obj (oTmp, ir, 0, h)
					      , gmp->Spc == WORM ? gtLowerCleanUp(ptr) : gtCleanUp(ptr)) ;
	  vtxtPrintf (tBfr, ooo) ;
	}

  /* GO_locuslink */
  if ((oTmp = ac_tag_table (gmp->gene, "GO_locuslink", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
      if (((ptr = ac_table_printable (oTmp, ir, 0, 0)) &&
	   ! strstr (ptr, "unknown") &&
	   !pickMatch (vtxtPtr (tBfr), messprintf("*%s*",gtLowerCleanUp(ptr)))))
	{
	  if (!jt2++)
	    TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Function"ooo) ;
	  TBB (jt, COL_VALUE) = vtxtPrintf (tBfr, "%s"ooo, 
					      gtCleanUp(ptr)) ;
	  if ((ptr = ac_table_printable (oTmp, ir, 1, 0)))
	    TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "%s"ooo, ptr) ;
	  jt++ ;
	}

   /* COG */
  if (0 && (oTmp = ac_tag_table (gmp->gene, "COG", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1; ir++)
      {
	oCog = ac_table_obj (oTmp, ir, 0, h) ;
	sprintf (linkBuf, COG_LINK, ac_name (oCog)) ; 
	TBB (jt, COL_ORIGIN) = gmpURL (tBfr, gmp, linkBuf, ac_name (oCog)) ; 
	vtxtPrintf (tBfr, ooo) ;
	if ((gCogInfo = ac_tag_table (oCog, "Info", h)))
	  for (jr=0;jr < gCogInfo->rows && gCogInfo->cols >= 2;jr++)
	    if ((ptr = ac_table_printable (gCogInfo, jr, 1, 0)) &&
		!pickMatch (vtxtPtr (tBfr), messprintf("*%s*",gtLowerCleanUp(ptr)))  &&
		!strstr(ptr,"rediction_only"))
	      {
		if (!jt2++)
		  TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Function"ooo) ;
		TBB (jt++, COL_VALUE) = vtxtPrintf (tBfr, "%s"ooo, 
						    gtCleanUp(ptr)) ;
	      }
	ac_free (oCog) ;
      }
 
  if ((oTmp = ac_tag_table (gmp->gene, "Go_m_ace", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1; ir++)
      {
	if ((ptr = ac_table_printable (oTmp, ir, 0, 0)) &&
	    ! strstr (ptr, "unknown") &&
	    !pickMatch (vtxtPtr (tBfr), messprintf("*%s*",gtLowerCleanUp(ptr))))
	  {
	    goType = ac_table_printable (oTmp, ir, 1, "") ;
	    if (ir != -1)
	      {
		if (gmp->Spc == HUMAN)
		  TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Entrez"ooo) ;
		else
		  TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Manual"ooo) ;
	      }
	    if (!jt2++)
	      TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Function"ooo) ;
	    TBB (jt++, COL_VALUE) = gmpFakeObjLink (tBfr, gmp, "Go_m", ac_table_obj (oTmp, ir, 0, h)
						, gmp->Spc == WORM ? gtLowerCleanUp(ptr) : messprintf ("[GOA%s%s] %s", goType[0] ? " " : "", goType,  gtCleanUp(ptr))) ;
	    vtxtPrintf (tBfr, ooo) ;
	  }
      }

  if ((oTmp = ac_tag_table (gmp->gene, "Go_m_pfam", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1; ir++)
      {
	if ((ptr = ac_table_printable (oTmp, ir, 0, 0)) &&
	    ! strstr (ptr, "unknown") &&
	    !pickMatch (vtxtPtr (tBfr), messprintf("*%s*",gtLowerCleanUp(ptr))) )
	  {
	    if (ir != -1)
	      TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Pfam"ooo) ;
	    if (!jt2++)
	      TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Function"ooo) ;
	    TBB (jt++, COL_VALUE) = gmpFakeObjLink (tBfr, gmp, "Go_m", ac_table_obj (oTmp, ir, 0, h)
						, gmp->Spc == WORM ? gtLowerCleanUp(ptr) : gtCleanUp(ptr)) ;
	    vtxtPrintf (tBfr, ooo) ;
	  }
      }

  jt2 = 0 ;
  jLocalization = 0 ;
  if ((oTmp = ac_tag_table (gmp->gene, "Go_c_ace", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1; ir++)
      {
	if ((ptr = ac_table_printable (oTmp, ir, 0, 0)) &&
	    !pickMatch (vtxtPtr (tBfr), messprintf("*%s*",gtLowerCleanUp(ptr))))
	  {
	    jLocalization++ ;
	    if (ir != -1)
	      {
		if (gmp->Spc == HUMAN)
		  TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Entrez"ooo) ;
		else
		  TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Manual"ooo) ;
	      }
	    if (!jt2++)
	      TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Localization"ooo) ;
	    goType = ac_table_printable (oTmp, ir, 1, "") ;
	    TBB (jt++, COL_VALUE) = gmpFakeObjLink (tBfr, gmp, "Go_c", ac_table_obj (oTmp, ir, 0, h)
						    , gmp->Spc == WORM ? gtLowerCleanUp(ptr) : messprintf ("[GOA%s%s] %s", goType[0] ? " " : "", goType,  gtCleanUp(ptr))) ;
	    vtxtPrintf (tBfr, ooo) ;
	  }
      }
  
  if ((oTmp = ac_tag_table (gmp->gene, "Go_c_pfam", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1; ir++)
      {
	if ((ptr = ac_table_printable (oTmp, ir, 0, 0)) &&
	    ! strstr (ptr, "unknown") &&
	    !pickMatch (vtxtPtr (tBfr), messprintf("*%s*",gtLowerCleanUp(ptr))))
	  {
	    jLocalization++ ;
	    if (ir != -1)
	      TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Pfam"ooo) ;
	    if (!jt2++)
	      TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Localization"ooo) ;
	    TBB (jt++, COL_VALUE) = gmpFakeObjLink (tBfr, gmp, "Go_c", ac_table_obj (oTmp, ir, 0, h)
						, gmp->Spc == WORM ? gtLowerCleanUp(ptr) : gtCleanUp(ptr)) ;
	    vtxtPrintf (tBfr, ooo) ;
	  }
      }

 
  if ((oTmp = ac_tag_table (gmp->gene, "Go_c_psort", h)))
    for (ir=0;ir < oTmp->rows && oTmp->cols >= 1; ir++)
      {
	if ((ptr = ac_table_printable (oTmp, ir, 0, 0)) &&
	    ! strstr (ptr, "unknown") &&
	    !pickMatch (vtxtPtr (tBfr), messprintf("*%s*",gtLowerCleanUp(ptr))))
	  {
	    jLocalization++ ;
	    if (ir != -1)
	      TBB (jt, COL_ORIGIN) = vtxtPrintf (tBfr, "Psort"ooo) ;
	    if (!jt2++)
	      TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Localization"ooo) ;
	    TBB (jt++, COL_VALUE) = gmpFakeObjLink (tBfr, gmp, "Go_c", ac_table_obj (oTmp, ir, 0, h)
						    , gmp->Spc == WORM ? gtLowerCleanUp(ptr) : gtCleanUp(ptr)) ;
	    vtxtPrintf (tBfr, ooo) ;
	  }
      }
  
  /* product */
  jt2 = 0 ;
  if (gmp->view == 'g' && 0 && gmp->Spc == WORM && /* this table is too repetitive */
      (gProducts = ac_tag_table (gmp->gene, "product", h)))
    {
      DICT *dict = 0 ;
      AC_OBJ obj, oProduct, oMrna ;

      if (gProducts->rows > 1)
	dict = dictCreate (0) ;
      for (ir = 0 ; ir < gProducts->rows ; ir++)
	{
	  oProduct = ac_table_obj (gProducts, ir, 0, h) ;
	  if ((obj = ac_tag_obj (oProduct, "mRNA", h)))
	    {
	      ptr = ac_name (oProduct) ;
	      if (dict && !dictAdd (dict, ptr, 0))
		continue ;
	      if (!jt2++)
		TBB (jt, COL_TYPE) = vtxtPrintf (tBfr, "Protein"ooo) ;
	      oMrna = ac_tag_obj (oProduct, "mRNA", h) ;
	      if (oMrna)
		{
		  TBB (jt, COL_ORIGIN) = gmpObjLink (tBfr, gmp, oMrna, ptr) ;
		  vtxtPrintf (tBfr, ooo) ;
		}
	      if ((ptr = gtProductMolecularName (bfr, gmp, oProduct)))
		TBB (jt, COL_VALUE) = vtxtPrintf (tBfr, "%s"ooo, ptr) ;
	      jt++ ;
	      ac_free (obj) ;
	    }
	  ac_free (oProduct) ;
	}
      ac_free (gProducts) ;
      dictDestroy (dict) ;
    }

  vtxtBreak (blkp) ;
  if (jt < 2)
    {
      if (0)
	vtxtPrintf (blkp, "Its function is still unknown") ;
    }
  else 
    {      
      if (gmp->view == 'm')
	vtxtPrintf (blkp, "The following table provides links to other genes in the same species with the same annotation:") ;
      else
	vtxtPrintf (blkp, "The following table provides links to other genes in the same species with similar annotations:") ;
    }

  vtxtEmptyLine (blkp, 1) ;
  if (1)fichePrintSquareTable (gmp, gmp->style, Tbb, blkp, tBfr, 0, 0, 0, maxRows, 
			 COL_TYPE, COL_VALUE, COL_ORIGIN, -1) ; 
  arrayDestroy (Tbb) ;
  vtxtDestroy (tBfr) ; 
  vtxtDestroy (bfr) ; 

  ac_free (h) ;

  return jt > 1 ? 1 : 0 ;
} /* ficheNewGeneDiseasePathwaysBioProcessTable */
#endif
/***************************************************************************************/
/***************************************************************************************/

static int ficheGeneExpressionLevelStatement (vTXT blkp, GMP *gmp, int previous)
{
  BOOL first = TRUE ;
  int ncl, level = -2 ;
  float magic = 0, xClo = -1 ; 
  AC_HANDLE h = ac_new_handle () ;
  const char *ptr, *cp[] =
  { "expressed at low level"
    , "moderately expressed"
    , "well expressed"
    , "expressed at high level"
    , "expressed at very high level"
  } ;
  
  if (first)
    {
      int ncl = 0, ntg = 0 ;
       
      if (gmp->Spc == HUMAN)
	{
	  /* tg = non cloud genes, clones = all clones in non cloud genes
	   * BUILD           30     31/33       34       35
	   * #tg          66830     89459   284760   100926 
	   * #clones    1722358   2605832  3630566  3606568
           * ratio           25        29       13       36
	   */
	  /* 34  */ ntg = 45757 ; ncl = 3237203 ;
	  /* 35g */ ntg = 52954 ; ncl = 4018074 ; /* each read counts as one clone */
	  /* 36a */ ntg = 57926 ; ncl = 3867490 ; /* read pairs used */
	  /* 36a */ ntg = 36268 ; ncl = 3925078 ; /* tg with introns, clones+buried clones */
	}       
      else if (gmp->Spc == WORM)
	{
	  /* BUILD WS130  WS140  140-aug05  dec-2006  feb-2007 
	   *      11547   15372    16704    15732   15892 
	   *     145451  173904   179041   176069   177331
	   *        query find gene balise && transcribed_gene (not just a model)
	   *        query {find est is_buried_under ; > cdna_clone ; ! from_gene} SETOR {find cdna_clone from_gene}
	   */
	  ntg = 15892 ; ncl = 177331 ;
	  ntg = 16420 ; ncl = 216139 ; /* 2008_10_18 */
	}
      else if (gmp->Spc == ARA)
	{
	  /* BUILD ara07   ara08   araSep08
	   *      27600    29064    27660 query find gene balise && transcribed_gene (not just a model)
	   *     977682  1100576
	   *        query find est is_buried_under ; > cdna_clone ; ! from_gene
	   *  SETOR query find cdna_clone from_gene
	   */
	  ntg = 29064 ; ncl = 1100576 ; /* 2008_02_19 */
	  ntg = 27660 ; ncl =  977682 ; /* 2008_09_29 */
	}
      else if (gmp->Spc == MOUSE)
	{
	  /* Clone per spliced gene */
	  ntg = 32249 ; ncl = 3332028 ;
	}
      else if (gmp->Spc == RAT)
	{
	  /* BUILD araSep08
	   *       27202        query find gene balise && transcribed_gene (not just a model)
	   *      682260
	   *        query find est is_buried_under ; > cdna_clone ; ! from_gene
	   *  SETOR query find cdna_clone from_gene
	   */
	  ntg = 27202 ; ncl = 682260 ; /* 2008_09_29 */
	}
      else if (gmp->db)
	{
	  ncl = ac_keyset_count (ac_dbquery_keyset (gmp->db, "Find cdna_clone from_gene", h)) ;
	  ntg = ac_keyset_count (ac_dbquery_keyset (gmp->db, "Find tg", h)) ;
	}
      
      if (ntg > 1000) magic = (ncl + 1.0)/(ntg + 1.0) ;
      first = FALSE ; 
      if (0) printf ("magic = %g ncl=%d ntg=%d\n", magic, ncl, ntg) ;
    }
  
  if (magic && gmp->gene)
    {
      AC_KEYSET oTmp = 0 ;

      if (gmp->gene && ac_has_tag (gmp->gene, "Has_cDNA_clone"))
	oTmp = ac_objquery_keyset (gmp->gene, ">Has_cDNA_clone", h) ;
      else
	oTmp =  ac_objquery_keyset (gmp->gene, ">transcribed_gene ; >read ; {IS *} SETOR {>buries;>buried_est;}; !Ref_seq ; > cdna_clone", h) ;

      ncl = oTmp ? ac_keyset_count (oTmp) : 0 ;
      xClo = ncl/(float) magic ;
      if (ncl == 0 && gmp->tg && ac_has_tag (gmp->tg, "cDNA_clone")) 
	{
	  oTmp =  ac_objquery_keyset (gmp->gene, ">transcribed_gene ; >read ; {IS *} SETOR {>buries;>buried_est;}; Ref_seq ; > cdna_clone", h) ;
	  ncl = oTmp ? ac_keyset_count (oTmp) : 0 ;
	  if (ncl)
	    level = -1 ;
	}
      else
	{
	  if (ncl > 0.0) level = 0 ;
	  if (ncl > 0.2 * magic) level = 1 ;
	  if (ncl > 0.4 * magic) level = 2 ;
	  if (ncl > 1.4 * magic) level = 3 ;
	  if (ncl > 4.0 * magic) level = 4 ;
	}
    }

  if (level == -1)
    {
      vtxtPrintf (blkp, "This gene is predicted and so far supported only by %sRefSeq", ncl == 1 ? "a " : "") ;
    }
  else
    {
      if (level >= 0)
	{
	  if (0 && gmp->Spc == WORM)
	    vtxtPrintf (blkp, "According to the Worm Transcriptome Project") ;
	  else
	    {
	      char *species ;
	      
	      switch (gmp->Spc)
		{
		case WORM: species = "worm" ; break ;
		default: species = "human" ; break ;
		}
	      vtxtPrintf (blkp, "According to ") ;
	      if (gmp->style == 'r')
		{
		  vtxtPrintf (blkp
			      , "<a href=\"https://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/av.cgi?db=%s&q=%s\">AceView</a>"
			      , species, ac_name(gmp->gene)) ;
		}
	      else
		vtxtPrintf (blkp, "AceView") ;
	    }
	  if (blkp->markUp)
	    vtxtPrintf (blkp, "%s is <a href=\"javascript:openAnchor('fexp','tg_expression')\">%s</a>, ", FALSE && previous ? ", it" : ", this gene", cp[level]) ;
	  else
	    vtxtPrintf (blkp, "%s is %s, ", FALSE && previous ? ", it" : ", this gene", cp[level]) ;
	  if (xClo > 1.4)
	    { 
	      vtxtPrint (blkp, messprintf ("%2.1f times", xClo)) ;
	    }
	  else if (xClo > .4)
	    {
	      vtxtPrint (blkp, messprintf ("%2.1f times", xClo)) ;
	    }
	  else if (xClo > 0)
	    {
	      vtxtPrint (blkp, messprintf ("only %2.1f%%", (int) 100*xClo)) ;
	      vtxtPrint (blkp, " of") ;
	    }
	  vtxtPrintf (blkp, " the average gene in this release") ;
	}
      
      if ((ptr = ac_tag_printable (gmp->gene, "Expression_title", 0)))
	{
	  if (level++ >= 0) vtxtPrintf (blkp, ", %s", ptr) ;
	  else { vtxtDot (blkp) ; vtxtPrintf (blkp, "This gene is expressed %s", ptr) ; }
	}
      if (gmp->Spc == WORM)
	{ 
	  ficheMRNAExpressionProfileParagraphContent (blkp, gmp, gmp->view == 'm', FALSE) ;
	}

      if (gmp->Spc == WORM)
	{ 
	  const char *ccp ;
	  int ir, n2, number ;
	  DICT *dict = gtGeneAliases (gmp, FALSE) ;
	  char	linkBuf[vONLINEMAXURLSIZE] ;

	  if (dict && (n2 = dictMax (dict)) > 1)
	    {
	      for (ir = 1  ; ir <= n2  ; ir++)
		{
		  ccp = dictName (dict, ir) ;
		  
		  if (!strncmp(ccp, "YK", 2) && (sscanf(ccp+2,"%d",&number) == 1)) 
		    {
		      vtxtDot (blkp) ; vtxtPrintf (blkp, "See the in situ hybridization pattern in ") ; 
		      sprintf (linkBuf, "http://nematode.lab.nig.ac.jp/db2/ShowGeneInfo.php?celk=CELK%05d",number) ;
		      gmpURL (blkp, gmp, linkBuf, "Kohara NextDB")  ;
		      break ;
		    }
		}
	    }
	  dictDestroy (dict) ;
	}
      
    }
  ac_free (h) ;
  return level >= -1 ? 1 : 0 ;
}

static BOOL ficheGmpHighlyMutableStatement (vTXT blkp, GMP *gmp, BOOL details)
{ 
  AC_TABLE oo = ac_tag_table (gmp->gene, "Highly_mutable", 0) ;
  const char *cp = ac_tag_printable (gmp->gene, "Highly_mutable", 0) ;
  if (cp)
    {
      if (details) vtxtBreak (blkp) ;
      vtxtPrintf (blkp, "mRNAs from this gene (or the genome) appear ") ;
      vtxtBold (blkp, "unstable") ;
      if (details)
	vtxtPrintf (blkp, " %s", cp) ;
      ac_free (oo) ;
      return TRUE ;
    }
  ac_free (oo) ;
  return FALSE ;
} /* ficheGmpHighlyMutableStatement */

static BOOL ficheGmpG26Statement (vTXT blkp, GMP *gmp, BOOL details)
{ 
  BOOL ok = FALSE ;
  int x, y1 ;
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl3, tbl2 = 0, tbl = ac_tag_table (gmp->gene, "Anti_G26", h) ;

  if (tbl && (x = ac_table_int (tbl, 0, 0, 0)) > 4)
    {
      tbl2 =  ac_tag_table (gmp->gene, "Exonic_RNA_antisense", h) ;   /* mieg 2011_10_02, this tag is G26 specific */
      y1 = tbl2 ? tbl2->rows : 0 ;
      vtxtBreak (blkp) ;
      vtxtPrintf (blkp, "This gene appears to be the target of ") ;
      
      if (! details)
	{
	  if (blkp->markUp)
	    vtxtPrintf (blkp, "<a href=\"javascript:openAnchor ('ffunc','tg_regulation')\">'26G' endo si-RNA</a>") ;
	  else
	    vtxtPrintf (blkp, "'26G' endo si-RNA") ;
	}
      else
	vtxtBold (blkp, "'26G' endo si-RNA") ;

      vtxtPrintf (blkp, "[Han et al. 2009, see the <a href='https://www.aceview.org/av.cgi?db=worm&c=extern&q=g26'>top target genes</a>]:  %d '26G' RNA sequences, from %d distinct positions, are antisense to the mRNA. Expression profile: %d tags were sequenced from embryos, %d from oocytes, %d from sperm and %d from mixed stages N2; %s 26G %s observed in glp-4(bn2) germ line defective mutant", x, y1, ac_table_int (tbl, 0, 6, 0), ac_table_int (tbl, 0, 5, 0), ac_table_int (tbl, 0, 4, 0), ac_table_int (tbl, 0, 1, 0), isOne (ac_table_int (tbl, 0, 3, 0)), _waser(ac_table_int (tbl, 0, 3, 0))) ;

      if (details)
	{
	  vtxtBreak (blkp) ;
	  vtxtPrintf (blkp, "[communicated by John Kim, 2009] The 26G endogenous small interfering RNAs (siRNAs) are 26 nt in length and possess a 5' guanine nucleotide (Ruby et al 2006).  They are all expressed in the antisense orientation to their target mRNAs and appear to degrade their targets via the canonical RNA interference mechanism.  Biogenesis of the 26G endosiRNAs requires the transcription of their target mRNAs and depends on the activity of the RRF-3 RNA-dependent RNA polymerase, as well as the ERI-1 exonuclease and the DCR-1 ribonuclease.  Two largely non-overlapping subclasses of 26G endosiRNAs are generated in the developing male and hermaphrodite germlines. Class I 26G endosiRNAs target genes are expressed during spermatogenesis, whereas class II 26G RNAs are maternally inherited and silence gene expression during oogenesis and subsequent zygotic development (Han et al 2009)" ) ; 
	  
	  tbl3 = ac_tag_table (gmp->gene, "Exonic_RNA_antisense", h) ;  /* mieg 2011_10_02, this tag is G26 specific */
	  if (tbl3 && tbl3->rows)
	    {
	      AC_TABLE tbl2 ;
	      vTXT bfr = vtxtHandleCreate (h) ;
	      const char *ccp, *ccq ;
	      int ir, jr ;
	      AC_OBJ Seq ;
	      int maxCols = 6 ;
	      Array Tbb =  arrayHandleCreate (12*100, int, h) ;
	      char *secColNames[] = {"26G endo si-RNA", "Total count", "Mixed stage N2", "eri-1", "glp-4(bn2)", "sperm", "oocyte", "embryo", 0} ;
	      enum {COL_KIM, COL1, COL2, COL3, COL4, COL5, COL6, COL7, COL_LAST} ;
	      
	      /* initialize the column titles */
	      vtxtPrintf (bfr, "toto"ooo) ;
	      for (maxCols=0 ; secColNames[maxCols] ; maxCols++)
		TBB(0, maxCols) = vtxtPrintf (bfr, "%s"ooo, secColNames[maxCols]) ; 
	      
	      for (ir = 0 ; ir < tbl3->rows ; ir++)
		{
		  ccp = ac_table_printable (tbl3, ir, 0, 0) ;
		  ccq = strstr (ccp, "_g") ;
		  if (! ccq) continue ;
		  TBB (ir+1, 0) = vtxtPrintf (bfr, "%s"ooo, ccq+1) ;
		  if ((Seq = ac_table_obj (tbl3, ir, 0, h)))
		    {
		      tbl2 = ac_tag_table (Seq, "Expression_profile", h) ;
		      for (jr = 0 ; tbl2 && jr < 7 ; jr++)
			TBB(ir+1, jr+1) = vtxtPrintf (bfr, "%s"ooo, ac_table_printable (tbl2, 0, jr, "-")) ;
		      ac_free (tbl2) ;
		      ac_free (Seq) ;
		    }
		}
	      fichePrintSquareTable (gmp, gmp->style, Tbb, blkp, bfr, 0, 0, 0, ir+1,
				     COL_KIM, COL1, COL2, COL3, COL4, COL5, COL6, COL7, -1) ;
	    }
	}
      
      ok = TRUE ;
    }
  ac_free (h) ;
  return ok ;
} /* ficheGmpG26Statement */

static BOOL ficheGmpRnaEditingStatement (vTXT blkp, GMP *gmp, BOOL details)
{ 
  AC_TABLE oo = ac_tag_table (gmp->gene, "RNA_editing", 0) ;
  const char *cp = ac_tag_printable (gmp->gene, "RNA_editing", 0) ;
  if (cp)
    {
      vtxtBreak (blkp) ;
      vtxtPrintf (blkp, "Conspicuous ") ;
      vtxtBold (blkp, "RNA editing") ;
      vtxtPrintf (blkp, " was observed in this gene by inspection of the cDNA sequencing traces") ;
      if (details)
	vtxtPrintf (blkp, ". Our notes mention %s", cp) ;
      ac_free (oo) ;
      return TRUE ;
    }
  ac_free (oo) ;
  return FALSE ;
} /* ficheGmpRnaEditingStatement */

static BOOL ficheGmpSelenocysteineStatement (vTXT blkp, GMP *gmp, BOOL details)
{ 
  AC_TABLE oo = ac_tag_table (gmp->gene, "Selenocysteine", 0) ;
  const char *cp = ac_tag_printable (gmp->gene, "Selenocysteine", 0) ;
  if (cp) 
    {
      vtxtBreak (blkp) ;
      vtxtPrintf (blkp, "Translating ") ;
      vtxtBold (blkp, "one stop codon") ;
      vtxtPrintf (blkp, " into an amino acid, standard or modified, would significantly increase the ORF") ;
      if (details)
	vtxtPrintf (blkp, ": %s", cp) ;
      ac_free (oo) ;
      return TRUE ;
    }
  ac_free (oo) ;
  return FALSE ;
} /* ficheGmpSelenocysteineStatement */

static BOOL ficheGmpTranslational_frameshiftStatement (vTXT blkp, GMP *gmp, BOOL details)
{ 
  AC_TABLE oo = ac_tag_table (gmp->gene, "Translational_frameshift", 0) ;
  const char *cp = ac_tag_printable (gmp->gene, "Translational_frameshift", 0) ;
  if (cp) 
    { 
      vtxtDot (blkp) ;
      if (details)
	vtxtPrintf (blkp, "Use of translational frameshift would significantly increase the ORF: %s.<br/>\n", cp) ;
      else
	vtxtPrintf (blkp, "Use of translational frameshift would significantly increase the ORF") ;
      ac_free (oo) ;
      return TRUE ;
    }
  ac_free (oo) ;
  return FALSE ;
} /* ficheGmpTranslational_frameshiftStatement */

/*************/

static char *ficheChromName (const char *mapNam)
{
  static char buf[1024] ;
  char *cp = 0 ;
  
  if (mapNam)
    {
      strncpy (buf, mapNam, 1023) ;
      cp = strstr (buf, "|") ;
      if (cp) *cp = 0 ;
      cp = buf ;
      if (!strncasecmp (buf, "chromosome_", 11))
	cp += 11 ;
    }
  return cp && *cp ? cp : 0 ; ;
}  /* ficheChromName */

/* -===================================- /
 * -=  GMP Map Description            =- /
 * -===================================- */	
/* -===================================- /
* -=  Position on Chromosome          =- /
* -===================================- */	
/* -===================================- /
* -=  GENE Phenotype and Descriptor   =- /
* -===================================- */


static void ficheGmpYkImageStatement (vTXT blkp, GMP *gmp)
{
  /* hum-5: gene->tg->cdna_clone->image */
  int ir, nn = 0, n2, number ;
  const char *cp ;
  char	linkBuf[vONLINEMAXURLSIZE] ;
  DICT *dict = 0 ;

  dict = gtGeneAliases (gmp, FALSE) ;

  if (dict && (n2 = dictMax (dict)) > 1)
    {
      for (ir = 1  ; !nn && ir <= n2  ; ir++)
	{
	  cp = dictName (dict, ir) ;

	  if (gmp->Spc == WORM && !strncmp(cp, "YK", 2) && (sscanf(cp+2,"%d",&number) == 1)) 
	    { 
	      vtxtBreak (blkp) ;
	      vtxtBold (blkp, "In situ") ;
	      vtxtBold (blkp, " hybridisation pictures to all stages of development are available from Kohara ") ;

  sprintf (linkBuf, "http://nematode.lab.nig.ac.jp/db2/ShowGeneInfo.php?celk=CELK%05d",number) ;
	      gmpURL (blkp, gmp, linkBuf, "NextDB")  ;
              if (gmp->markup)
		vtxtPrintf (blkp
			    , " (may be slow, but please be patient, it is worth it !)"
			    ". Note that this extraordinary resource is still" 
			    " unpublished: Users of the data presented on the NextDB web pages should" 
			    " not publish the information without ") ;
	      if (gmp->markup)
		vtxtPrintf (blkp, "<a href=\"mailto:ykohara@lab.nig.ac.jp, tshini@genes.nig.ac.jp\")>Kohara's</a>") ; 
	      else
		vtxtPrintf (blkp, "Kohara's (ykohara@lab.nig.ac.jp, tshini@genes.nig.ac.jp") ;
	      vtxtPrintf (blkp
			  , " permission and appropriate acknowledgment." 
			  " They also appreciate any suggestions and comments."
			  ) ;
	      nn++ ;
	    }
	}
    }
  dictDestroy (dict) ;
} /* ficheGmpYkImageStatement */

static void ficheGmpPatternStatement (vTXT blkp, GMP *gmp)
{
  /* ace-1 */
  int ir, nn = 0 ;
  const char *cp ;
  AC_TABLE gPat = ac_tag_table (gmp->gene, "Pattern", 0) ;
  
  if (gPat)
    {
      for (ir = 0  ; ir < gPat->rows ; ir++)
	{
	  cp = ac_table_printable (gPat, ir, 0, 0) ;
	  if (cp)
	    {
	      if (!nn++)
		{
		  vtxtBreak (blkp) ;
		  vtxtBold (blkp, "Pattern:") ;
		  vtxtPrintf (blkp, " %s"
			      , gtCleanUp (cp)) ;
		}
	      else
		vtxtPrintf (blkp, ". %s", gtCleanUp (cp)) ;
	    }
	}
    }
  ac_free (gPat) ;
} /* ficheGmpPatternStatement */

/* -===================================- /
 * -=  EXPRESIION PROFILE              =- /
 * -===================================- */	

static void ficheMRNAExpressionProfileParagraphContent (vTXT blkp, GMP *gmp, BOOL isMrna, BOOL isBold)
{
  int ic, jj ; 
  AC_TABLE gExpression_profile = 0 ;
  int profNum ; 
  AC_HANDLE h = 0 ;

  if (gmp->Spc != WORM ) return ;
  h = ac_new_handle () ;
  if (isMrna && gmp->mrna && gmp->tg)
    gExpression_profile = ac_tag_table (gmp->mrna, "Expression_profile", h) ;
  else if (!isMrna && gmp->tg)
    gExpression_profile = ac_tag_table (gmp->tg, "Expression_profile", h) ;
  else if (!isMrna && gmp->pg && gmp->gene)
    gExpression_profile = ac_tag_table (gmp->gene, "Has_expression_profile", h) ;
    
  if (gExpression_profile && gExpression_profile->cols >= 3)
    {
      vtxtDot (blkp) ;  
      vtxtPrintf (blkp, "The ") ;
      if (isBold)
	vtxtBold (blkp, "expression profile") ;
      else
	vtxtPrintf (blkp, "expression profile") ;
      if (!gmp->markup || !isMrna)
	vtxtPrintf (blkp, " for the gene, derived from the proportion of animals at each stage in each Kohara library is: ") ; 
      else 
	vtxtPrintf (blkp, " for mRNA%s%s is: "
		    , gmp->variant ? " variant " : ""
		    , gmp->variant) ; 
      
      for (ic = jj = profNum = 0 ; ic < gExpression_profile->cols && ic < 7 ; ic++)
	{
	  profNum += ac_table_int (gExpression_profile, 0, ic, 0) ;
	  switch (ic)
	    {
	    case 0:
	      if (profNum)
		vtxtPrintf (blkp, "%s%s %d%%", _comma (jj++), "embryos", profNum) ;
	      profNum = 0 ;
	      break ;
	    case 2:
	      if (profNum)
		vtxtPrintf (blkp, "%s%s %d%%", _comma (jj++), "L1 or L2 larvae", profNum) ;
	      profNum = 0 ;
	      break ;
	    case 6:
	      if (ac_table_int (gExpression_profile, 0, 6, 0))
		vtxtPrintf (blkp, "%s%s %d%%", _comma (jj++), "L3 to adult (including dauer)", profNum) ;
	      else
		vtxtPrintf (blkp, "%s%s %d%%", _comma (jj++), "L3 to adult", profNum) ;
	      profNum = 0 ;
	      break ;
	    }
	}
    }
  ac_free (h) ;
  return ;
}

/* -===================================- /
* -=  EXPRESSION PATTERN              =- /
* -===================================- */	

static void ficheGmpExpressionPatternStatement (vTXT blkp, GMP *gmp)
{
  int ir, jr ;
  AC_TABLE oPat ;
  AC_OBJ oP ;
  char linkBuf[vONLINEMAXURLSIZE] ; 
  AC_HANDLE h = ac_new_handle () ;

  if ((oPat = ac_tag_table (gmp->gene, "Expr_pattern", h)) &&
      ac_has_tag (gmp->gene, "Expr_pattern"))
    {
      for (ir = jr = 0 ; ir < oPat->rows ; ir++)
	{
	  oP = ac_table_obj (oPat, ir, 0, h) ;
	  if (! ac_has_tag (oP, "Not_done"))
	    continue ;
	  if (!jr++)
	    { 
	      vtxtBreak (blkp) ;
	      vtxtPrintf (blkp, "For a detailed expression pattern description, see Wormbase ") ; 
	    }
	  sprintf (linkBuf, "http://www.wormbase.org/db/gene/expression?name=%s;class=Expr_pattern"
		   , ac_name (oP)) ;
	  if (ir >= 1 ) vtxtPrintf (blkp, ", ") ;
	  gmpURL (blkp, gmp, linkBuf, ac_name (oP)) ; 
	}
    }
  ac_free (h) ;  
}
 

/* -===================================- /
* -=  CONCEPTUAL TRNASLATION          =- /
* -===================================- */	
void ficheMRNAConceptualTranslationParagraphContent (vTXT blkp, GMP *gmp)
{
  int lenGap = ac_tag_int (gmp->mrna, "Gap_length", 0) ; 
  int lenC = 0, lenCoding = 0, isComplete = 0, firstMet = 0, nn = 0 ; 
  AC_TABLE gMW = 0 ; 
  AC_HANDLE h = 0 ;
  
  if (!gmp->product)
    return ;
  if (!gmp->peptide)
    gmp->peptide = ac_obj_peptide (gmp->product, gmp->h) ; 
  if (gmp->peptide)
    lenCoding = strlen (gmp->peptide) ;
  lenC = gmp->product ? ac_tag_int (gmp->product, "Coding_Length", 0) : 0 ;
  lenC /= 3 ;
  if (lenC > 0 && lenC < lenCoding)
    lenCoding = lenC ;
   h = ac_new_handle () ;
   if ( lenGap == 0)
    {
      vtxtPrintf (blkp, "The ") ; 
      if ( ac_has_tag (gmp->product, "at_position_1") &&
	   (
	   ac_has_tag (gmp->product, "Complete") ||
	   ( ac_has_tag (gmp->product, "NH2_Complete") && ac_has_tag (gmp->product, "COOH_Complete"))
	   )
	   )
	{
	  vtxtPrintf (blkp, "complete") ; isComplete=1 ; 
	}
      else vtxtPrintf (blkp, "partial") ; 
      
      if (!ac_has_tag (gmp->product, "at_position_1") ||
	  (!ac_has_tag (gmp->product, "Complete") && !ac_has_tag (gmp->product, "COOH_Complete"))
	  )
	vtxtPrintf (blkp, " ORF") ; 
      else {if (isComplete)vtxtPrintf (blkp, " protein") ; else vtxtPrintf (blkp, " CDS") ; }
      vtxtPrintf (blkp, " encoded between the first") ; 
      
      if (ac_has_tag (gmp->product, "Met") && ac_has_tag (gmp->product, "at_position_1"))
	{
	  int x1, ir ;
	  char buf[4] ;
	  AC_OBJ oProd ;
	  AC_TABLE oTmp = ac_tag_table (gmp->mrna, "Product", 0) ;
	  
	  for (ir = 0 ; oTmp && ir < oTmp->rows; ir++)
	    {
	      oProd = ac_table_obj (oTmp, ir, 0, h) ;
	      if (! ac_has_tag (oProd, "Best_product"))
		continue ;
	      if (!strcmp (ac_name (oProd), ac_name (gmp->product))&&
		  (x1 = ac_table_int (oTmp, ir, 1, 0)) &&
		  ( gmp->dna || (gmp->dna = ac_obj_dna (gmp->mrna, gmp->h))))
		{
		  strncpy(buf, gmp->dna + x1 - 1, 3) ;
		  buf[3] = 0 ;
		  if (*buf && ace_lower (*(buf+1)) == 't' && ace_lower (*(buf+2)) == 'g')
		    vtxtPrintf (blkp, " Met (%s)", buf) ;
		  else
		    vtxtPrintf (blkp, " Met") ;
		}
	      else
		vtxtPrintf (blkp, " Met") ;
	    }
	  nn++ ;
	  ac_free (oTmp) ;
	}
      else 
	vtxtPrintf (blkp, " in frame amino acid") ; 
      
      vtxtPrintf (blkp, " and the %s codon contains "
		 , ac_has_tag (gmp->product, "Complete") || ac_has_tag (gmp->product, "COOH_Complete") 
		 ? "stop":"last") ;
      {
	char buf1[256], buf2[256] ;
	
	sprintf (buf1, "PEP_Product::") ;
	sprintf (buf2, "%d", lenCoding) ;
	gmpFakeObjLink (blkp, gmp, buf1, gmp->product, buf2) ;
      }
      vtxtPrint (blkp, " residues") ;
      
      if (!nn && gmp->product)
	{
	  int openL = ac_tag_int (gmp->product, "Open_length", 0) ;
	  firstMet = ac_tag_int (gmp->product, "First_Met", 0) + ac_tag_int (gmp->product, "First_ATG", 0) ;
	  if (3*firstMet < openL)
	    vtxtPrintf (blkp
			,", the first Met (atg) is amino acid %d"
			, firstMet) ;
	  else
	    vtxtPrintf (blkp
			, ", there is no Met (atg) codon in this CDS") ;
	}
      if (gmp->product &&
	  (gMW = ac_tag_table (gmp->product, "Molecular_weight", h)) &&
	  ac_has_tag (gmp->product, "Complete") &&
	  gMW && gMW->rows == 1 && gMW->cols >=2 && 
	  ac_table_float (gMW, 0, 1, -1) != -1)
	{
	  vtxtDot (blkp) ;
	  vtxtPrintf (blkp, "The calculated molecular weight of the") ; 
	  vtxtPrintf (blkp, " protein is %.1f kDa, its isoelectric point %.2g", 
		      ac_table_float (gMW, 0, 0, 0), ac_table_float (gMW, 0, 1, 0)) ; 
	}  
    }
  else
    {
      vtxtPrintf (blkp, 
		  " There is a gap estimated (not reliably)"
		  " to %d bp"
		  , lenGap) ; 
    }
  ac_free (h) ;
  return ;
}

void ficheMRNAConceptualTranslationParagraph (vTXT blkp, GMP *gmp)
{
  char *ptr ;
  vTXT buf = vtxtCreate () ;
  if (gmp->markup) vtxtMarkup (buf) ;

  ficheMRNAConceptualTranslationParagraphContent (buf, gmp) ; 
  if ((ptr = vtxtPtr (buf)))
    {
      gmpSection (blkp, gmp, "MW", "Conceptual translation, MW, pI") ;
      vtxtPrint (blkp, ptr) ;
    }
  vtxtDestroy (buf) ;
}

/* -===================================- /
 * -=  PSORT Paragraph                 =- /
 * -===================================- */	
void fichePRODUCTPsortParagraphContent (vTXT blkp, GMP *gmp, BOOL justPsortLocalization)
{
  AC_TABLE gLocalization, gDomain ;
  AC_OBJ oTmp ;
  int isProductComplete=0, ir, isComma, tt ; 
  const char	*prefix, *tTxt ;
  char	objdate[256], seqBuf[128] ;
  float percentage ; 
  Array Tbb = 0 ;
  vTXT bfr ;
  AC_HANDLE h = ac_new_handle () ;
  char *colNames[]={ "From aa", 
		     "Domain", /* in worm */
		     "Sequence",
		     0} ; 
  
  enum { COL_WHERE=0, COL_DOMAIN, COL_SEQ, COL_LAST} ; 
  
  
  if (! gmp->kantor ||
      ! (gLocalization = ac_tag_table (gmp->kantor, "Localization", h)) )
    {
      ac_free (h) ;
      return ; 
    }
  
  if (
      (ac_has_tag (gmp->product, "NH2_Complete") && ac_has_tag (gmp->product, "COOH_Complete")) ||
      ac_has_tag (gmp->product, "Complete") || (gmp->mrna && ac_has_tag (gmp->mrna, "From_prediction"))
      )
    isProductComplete=1 ; 
  
  timeShowFormat (ac_tag_date (gmp->kantor, "Psort_Date", timeNow ()), "%b %d, %Y", objdate, sizeof (objdate)) ;  
  
  vtxtPrintf (blkp, "PSORT II analysis, (K. Nakai") ; 
  gmpURL (blkp, gmp, " http://psort.nibb.ac.jp", 0) ; 
  vtxtPrintf (blkp, ") trained on yeast data") ; 
  if (0) vtxtPrintf (blkp, " and run on %s", objdate) ; 
  
  vtxtPrintf (blkp, " predicts "
	      "that the subcellular location of this %sprotein is ", 
	      isProductComplete ? "" : "partial ") ; 
  prefix="most likely" ; 
  for (ir=0 ; ir < gLocalization->rows ; ir++)
    {
      oTmp = ac_table_obj (gLocalization, ir, 0, h) ; 
      /*				tTxt = ac_name (oTmp) ; 
					if (!strcmp (tTxt, "extracellular, including cell wall"))tTxt="secreted" ; */
      tTxt = fichePSORTNameAlias (ac_name (oTmp)) ; 
      percentage = ac_table_float (gLocalization, ir, 1, 0) ; 
      if (percentage < 30)break ; 
      vtxtPrintf (blkp, "%s%s%s (%d%%)", prefix, ir>0 ? " or " : " ", tTxt, (int) (percentage)) ; 
      prefix="" ; 
    }
  
  if (ir==0)prefix="possibly " ; 
  else prefix=". Less likely possibilities are" ; 
  isComma=0 ; 
  for ( ; ir < gLocalization->rows ; ir++)
    {
      oTmp = ac_table_obj (gLocalization, ir, 0, h) ; 
      percentage = ac_table_float (gLocalization, ir, 1, 0) ; 
      vtxtPrintf (blkp, "%s%s %s (%d%%)", prefix, isComma>0 ? " or" : "", fichePSORTNameAlias (ac_name (oTmp)), (int) (percentage)) ; 
      prefix="" ; 
      isComma++ ; 
    }
  
  if (justPsortLocalization)
    return  ; 
  
  bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  vtxtPrintf (bfr, "toto") ; /* avoid zero */  
  /* initialize the table */
  tt = 0 ;
  Tbb = arrayCreate (2000, int) ;
  for (ir=0 ; colNames[ir] ; ir++)
    TBB (0, ir) = vtxtPrintf (bfr, "%s"ooo, colNames[ir]) ; 
  
  /* kantor->domain is copied as product->psort_domain except for the last column 'motif' */
  if ((gDomain = ac_tag_table (gmp->kantor, "Domain", h)))
    {
      int cntDomains=0, *indDomains, currInd, i, j, tmp ; 
      const char *last = "" ;
      const char *sSequence ; 

      indDomains= (int *) messalloc(gDomain->rows*2*sizeof (int)) ; 
      
      for (ir=0 ; ir < gDomain->rows ; ir++)
	{ /* count unique elements */
	  if (strcmp (last, ac_table_printable (gDomain, ir, 0, "")))
	    { cntDomains++ ; last = ac_table_printable (gDomain, ir, 0, "") ; } 
	  if (indDomains)indDomains[ir]=ir ; 
        }
      
      if (indDomains)
	{
	  /* sort them */
	  for (i=0 ; i < gDomain->rows-1 ; i++)
	    {
	      for (j = i+1 ; j < gDomain->rows ; j++)
		{
		  if (ac_table_int (gDomain, indDomains[i], 3, -1) > 
		      ac_table_int (gDomain, indDomains[j], 3, 0))
		    {
		      tmp=indDomains[i] ; indDomains[i]=indDomains[j] ; indDomains[j]=tmp ;  /* swap */
		    }
		}
	    }
	}
      vtxtBreak (blkp) ;
      vtxtPrintf (blkp, 
		  "\nThe following domain%s %s found:\n", _multi (gDomain->rows)
		  , _waser (gDomain->rows)) ; 
      
      for (tt = ir = 0 ; ir < gDomain->rows ; ir++)
	{
	  if (indDomains)
	    currInd = indDomains[ir] ; 
	  else 
	    currInd = ir ; 
	  oTmp = ac_table_obj (gDomain, currInd, 0, h) ; 
	  sSequence = ac_table_text (gDomain, currInd, 7, 0) ; 
	  seqBuf[0]=0 ; strncpy (seqBuf, sSequence, 128 - 1) ; 
	  seqBuf[127]=0 ; 
	  
	  tt++ ;
	  TBB (tt, COL_WHERE) = vtxtPrintf (bfr, "%d to %d"ooo
					    ,ac_table_int (gDomain, currInd, 3, 0)
					    ,ac_table_int (gDomain, currInd, 4, 0)) ;
	  {
	    const char *cp = ac_name (oTmp) ;
	    if (!strncasecmp (cp, "Nuclear_locali", 14))
	      {
		cp = "Possible nuclear localization" ;
		if (gmp->Spc == HUMAN)
		  TBB (tt, COL_DOMAIN) = vtxtPrintf (bfr, "%s"ooo, cp) ; 
		else
		  {
		    TBB (tt, COL_DOMAIN) = gmpFakeObjLink (bfr, gmp, "psort", oTmp, cp) ;
		    vtxtPrintf (bfr, ooo) ;
		  }
	      }
	    else
	      {
		TBB (tt, COL_DOMAIN) = gmpFakeObjLink (bfr, gmp, "psort", oTmp, cp) ;
		vtxtPrintf (bfr, ooo) ;
	      }
	  }
	  if (sSequence[0] && ! strstr (sSequence, "recommended to reKantorize"))
	    {
	      TBB (tt, COL_SEQ) = vtxtPrintf (bfr, "%s", seqBuf) ; 
	      if (strlen (sSequence) >= (128-1))
		vtxtPrintf (bfr, "...") ;
	      vtxtPrintf (bfr, ooo) ;
	    }
        }
      messfree (indDomains) ; 
    } 

  fichePrintSquareTable (gmp, gmp->style, Tbb, blkp, bfr, 0, 0, 0, tt + 1, 
			  COL_WHERE, COL_DOMAIN, COL_SEQ, -1) ;
  arrayDestroy (Tbb) ; 
  vtxtDestroy (bfr) ;
  ac_free (h) ;
} /* fichePRODUCTPsortParagraphContent */

void fichePRODUCTPsortParagraph (vTXT blkp, GMP *gmp)
{
  vTXT bfr = vtxtCreate () ;
  if (gmp->markup) vtxtMarkup (bfr) ;
 
  if (!gmp->product) return ; 
  
  fichePRODUCTPsortParagraphContent (bfr, gmp, 0) ; 
  
  if (vtxtPtr (bfr))
    {
      gmpSection (blkp, gmp, "Psort"
		      , "Predicted cellular localization and motifs (Psort)") ;
      vtxtPrint (blkp, vtxtPtr (bfr)) ; 
    }
  vtxtDestroy (bfr) ;
} /* fichePRODUCTPsortParagraph */

/* -===================================- /
 * -=  PFAM Paragraph                  =- /
 * -===================================- */	
static void fichePRODUCTPfamParagraphContent (vTXT blkp, GMP *gmp)
{
  AC_TABLE gPfam ;
  AC_OBJ oFam ;
  int date, ir, cntNonSignificant ; 
  char	objdate[256], linkBuf[vONLINEMAXURLSIZE] ; 
  AC_HANDLE h = 0 ;

  date = ac_tag_date (gmp->kantor, "Pfam_Date", 0) ;
  if (!date)
    return ;
  h = ac_new_handle () ;
  vtxtPrintf (blkp, "Pfam analysis (") ; 
  gmpURL (blkp, gmp, "http://pfam.xfam.org", 0) ; 
  vtxtPrintf (blkp, ")") ; 
 
  timeShowFormat (date, "%b %d, %Y", objdate, sizeof (objdate)) ; 
  vtxtPrintf (blkp, " run on %s, shows ", objdate) ; 
  
  if ((gPfam = ac_tag_table (gmp->kantor, "Pfam", h)))
    {
      int cntType, jr, sr ;
      const char *ccp, *famDescr ;
      char familyName[vSTRMAXNAME] ;
      AC_TABLE gFamGene ;
      const char *txtAccession ; 
      
      for (ir=0 ; ir < gPfam->rows ; ir++)
	{
	  oFam = ac_table_obj (gPfam, ir, 0, h) ; 
	  if (!ac_tag_printable (oFam, "Accession", 0))
	    continue ;
	  sr=ir ; 
	  
	  for (cntType=1 ; 	ir + 1 <  (gPfam->rows) && 
		!strcmp (ac_table_printable (gPfam, ir+1, 0, ""), ac_name (oFam))
		 ; ir++)
	    {
	      cntType++ ; 
	    }
	  markupStart (blkp, gmp, _TB) ; markupStart (blkp, gmp, _TR) ; markupStart (blkp, gmp, _TD) ; 
	  if (sr > 0) vtxtPrintf (blkp, "\n    There %s also ", _isare (cntType)) ; 
	  if (cntType==1)vtxtPrintf (blkp, "a ") ; 
	  else if (cntType>1)vtxtPrintf (blkp, "%d ", cntType) ; 
	  
	  ccp = ac_tag_printable (oFam, "Definition", 0) ;
	  if (ccp)
	    strcpy (familyName, ccp ? ccp : "") ;
	  if (!ccp) /* try to grab it from synonymous pfam */
	    {
	      AC_ITER ol ;
	      AC_OBJ oFamTmp ;
	      ol = ac_query_iter (gmp->db, TRUE, messprintf ("Find Pfam IS \"%s\"; > Accession ; > Quoted_in ; CLASS Pfam ; Definition", ac_name(oFam)), 0, h) ;
	      if (ol)
		{
		  while (!ccp && (oFamTmp = ac_next_obj (ol)))
		    {
		      ccp = ac_tag_printable (oFamTmp, "Definition", 0) ; 
		      if (ccp)
			strcpy (familyName, ccp ? ccp : "") ;
		      ac_free (oFamTmp) ;
		    }
		}
	    }
	  vtxtPrintf (blkp, "significant hit%s to the ", _multi (cntType)) ; 
	  
	  if ((txtAccession = ac_tag_printable (oFam, "Accession", 0)))
	    {
	      char *cq = strstr (txtAccession , ".") ; /* drop pfam version number */
	      if (cq) *cq = 0 ;
	      sprintf (linkBuf, "http://pfam.xfam.org/family/%s", txtAccession) ; 
	      gmpURL (blkp, gmp, linkBuf, familyName) ; 
	    }
	  else vtxtPrint (blkp, familyName) ; 
	  
	  
	  for (jr=sr ;  jr-sr < cntType ; jr++)
	    {
	      const char *eValue = ac_table_printable (gPfam, sr, 8, 0) ;
	      vtxtPrint (blkp, jr==sr ? "" : (jr-sr==cntType-1 ? " and" : ",") ) ; 
	      if (cntType>1)
		vtxtBreak (blkp) ;
	      vtxtPrintf (blkp, " from %d to %d, with score %d "
			  , ac_table_int (gPfam, jr, 3, 0), ac_table_int (gPfam, jr, 4, 0)
			  , (int) ac_table_float (gPfam, jr, 2, 0)) ; 
	      if (eValue && *eValue)
		vtxtPrintf (blkp, "and E = %s"
			    , eValue) ; 
	    }
	  
	  vtxtBreak (blkp) ;
	  if ((gFamGene = ac_tag_table (oFam, "Gene", h)) )
	    {
	      int nGenes=0, nTg = -1 ; 
	      nGenes = gFamGene->rows ;

	      if (nGenes == 1)
		{
		  vtxtDot (blkp) ;
		  vtxtPrintf (blkp, "No other gene in the database contains this motif") ; 
		}
	      else  if (nGenes > 1)
		{
		  int ii, jj ; 
		  AC_OBJ obj ;

		  vtxtDot (blkp) ;
		  gmpObjLink (blkp, gmp, oFam, messprintf ("%d", nGenes-1)) ;
		  
		  vtxtPrintf (blkp, " other gene%s ", _multi (nGenes-1)) ; 
		  if (nGenes < 10)
		    {		      vtxtPrintf (blkp, "/") ; 
		      for (jj=0, ii=0 ; ii < nGenes ; ii++)
			{
			  obj = ac_table_obj (gFamGene, ii, 0, h) ;
			  if (!strcmp (ac_name (obj), ac_name (gmp->gene))) continue ; 
			  if (jj) vtxtPrintf (blkp, ", ") ; 
			  gmpObjLink (blkp, gmp,  obj, 0) ;
			  jj++ ; 
			}
		      vtxtPrintf (blkp, "/ ") ; 
		    }

		  vtxtPrintf (blkp, "in the database also contain%s this motif", _verbs (nGenes-1)) ; 

		  if (gmp->Spc == WORM)
		    {
		      nTg = ac_keyset_count (ac_objquery_keyset (oFam, "follow gene ;  Transcribed_gene || Expression_title || Pattern || Expr_pattern || Has_OST || Has_cDNA_clone", h)) ;	
		      if (gmp->tg || ac_has_tag (gmp->gene, "Expression"))
			nTg-- ;
		      if (nTg == 0)
			vtxtPrintf (blkp, ", none of which is known to be expressed") ;
		      else if (nGenes > 2 && nTg == 1)
			vtxtPrintf (blkp, "and only one is known to be expressed") ;
		      else if (nTg == 1)
			vtxtPrintf (blkp, "and is known to be expressed") ;
		      else if (nTg == nGenes - 1)
			vtxtPrintf (blkp, ", all of which are known to be expressed") ;
		      else if (nTg > 0)
			vtxtPrintf (blkp, ", %d of which are known to be expressed", nTg) ;
		    }
		}
	    }	  
	  
	  markupEnd (blkp, gmp, _TD) ; markupEnd (blkp, gmp, _TR) ; markupEnd (blkp, gmp, _TB) ; 
	  
	  if (gmp->markup &&
	      (famDescr = ac_tag_printable (oFam, "Comment", 0)))
	    {
	      AC_TABLE gRef ; 
	      int ii ; 
	      
	      markupStart (blkp, gmp, _TB) ; 				 
	      markupStart (blkp, gmp, _TR) ; 
	      markupText (blkp, gmp, _TD, "    ") ; 
	      markupStart (blkp, gmp, _TD) ; 
	      vtxtBold (blkp, "[InterPro/Pfam]") ;
	      vtxtPrintf (blkp, " %s", famDescr) ; 
	      markupEnd (blkp, gmp, _TD) ; 
	      markupEnd (blkp, gmp, _TR) ; 
	      
	      if ((gRef = ac_tag_table (oFam, "Reference", h)))
		{
		  for (ii=0 ; ii < gRef->rows ; ii++)
		    {
		      if (gmp->markup)
			{ /* spaces deteriorate the presentation of the table but */
			  markupStart (blkp, gmp, _TR) ; markupStart (blkp, gmp, _TD) ; 
			  vtxtPrintf (blkp, "[%d]:", ii+1) ; 
			  markupEnd (blkp, gmp, _TD) ; markupStart (blkp, gmp, _TD) ; 
			}
		      else /* spaces are needed to please sequin validator not to recognise this as a reference */
			vtxtPrintf (blkp, "[ %d ]:", ii+1) ; 
		      markupLinkPubmed (blkp, gmp, ac_table_printable (gRef, ii, 0, ""), ac_table_printable (gRef, ii, 1, "")) ; 
		      markupEnd (blkp, gmp, _TD) ; markupEnd (blkp, gmp, _TR) ; 
		    }
		}	 
	      markupEnd (blkp, gmp, _TB) ;
	    }
	}
    }
  
  if ((cntNonSignificant = ac_tag_int (gmp->kantor, "Pfam_Non_Significant_Hits", 0)))
    {
      vtxtPrintf (blkp, "%s%s %s%d non significant Pfam hit%s"
		 , (gPfam) ? "There " :""
		 , (gPfam) ? (_isare (cntNonSignificant)) : ""
		 , (gPfam) ? "also " :""
		 , cntNonSignificant, _multi (cntNonSignificant)) ; 
    }
  else if (!gPfam)
    {
      vtxtPrintf (blkp, "that this protein does not belong to any recognized protein family") ; 
    }
  ac_free (h) ;
  return ;
}

void fichePRODUCTPfamParagraph (vTXT blkp, GMP *gmp)
{
  vTXT bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;
    
  fichePRODUCTPfamParagraphContent (bfr, gmp) ; 
  if (vtxtPtr (bfr))
    {
      gmpSection (blkp, gmp, "mrnaPfam", "Protein domains") ;
      vtxtPrint (blkp, vtxtPtr (bfr)) ; 
    }
  vtxtDestroy (bfr) ;
} /* fichePRODUCTPfamParagraph  */

/* -===================================- /
* -=  BLASTP                  =- /
* -===================================- */	

void ficheProductBlastPTableContent (vTXT blkp, GMP *gmp, int maxCol, const char *more)
{
  AC_TABLE gBlastP ;
  AC_OBJ oBl ; 
  int	ir ; 
  char linkBuf[128+vSTRMAXNAME], buf[vSTRMAXNAME], *ptr ;
  const char *titles[] = {"score", "from aa to aa", "blastP Hit Title [species] eValue links to GenBank", "from aa to aa"} ; 
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE  table = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  
  if ((gBlastP = ac_tag_table (gmp->kantor, "BlastP", h)))
    {
      const char *ccp ;
      maxCol = gBlastP->rows ;
      table = ac_empty_table (maxCol, 4, h) ;
      for (ir=0 ; ir<gBlastP->rows && ir<maxCol ;  ir++)
	{
	  oBl = ac_table_obj (gBlastP, ir, 0, h) ; 
	  ccp = vtxtPtr (txt) ;
	  vtxtPrintf (txt, "%12.0f", ac_table_float (gBlastP, ir, 2, 0)) ;
	  ccp = strnew (vtxtPtr (txt), h) ;
	  ac_table_insert_text (table, ir, 0, ccp) ;

	  vtxtClear (txt) ;
	  vtxtPrintf (txt, "%4d to %4d", ac_table_int (gBlastP, ir, 3, 0), ac_table_int (gBlastP, ir, 4, 0) ) ;
	  ccp = strnew (vtxtPtr (txt), h) ;
	  ac_table_insert_text (table, ir, 1, ccp) ;

	  vtxtClear (txt) ;
	  strcpy (buf, ac_name (oBl)) ; 
	  if ((ptr=strchr (buf, ' ')))*ptr=0 ; 
	  sprintf (linkBuf, "https://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=%s", buf) ; 
	  if ((ccp = ac_table_printable (gBlastP, ir, 7, 0)))
	    gmpURL (txt, gmp, linkBuf, ccp) ; 
	  else
	    gmpURL (txt, gmp, linkBuf, "segment of previous") ; 
	  ccp = strnew (vtxtPtr (txt), h) ;
	  ac_table_insert_text (table, ir, 2, ccp) ;

	  vtxtClear (txt) ;
	  vtxtPrintf (txt, "%4d to %4d", ac_table_int (gBlastP, ir, 5, 0), ac_table_int (gBlastP, ir, 6, 0)) ;
	  ccp = strnew (vtxtPtr (txt), h) ;
	  ac_table_insert_text (table, ir, 3, ccp) ;
	}
      if (ir < gBlastP->rows)
	{ 
	  for (int icol = 0 ; icol < 4 ; icol++)
	    ac_table_insert_text (table, ir + 1, icol, more) ;
	}
    }
  ac_table_display (blkp, table, titles, 0, 0, 0, 0, 0, 0) ;
  ac_free (h) ;
} /* ficheProductBlastPTableContent */

void fichePRODUCTBlastPParagraph (vTXT blkp, GMP *gmp)
{
  vTXT bfr = vtxtCreate () ;
  char *ptr ;
  if (gmp->markup) vtxtMarkup (bfr) ;

  if ((ptr = fichePRODUCTBlastPParagraphContent (bfr, gmp)))
    {
      gmpSection (blkp, gmp, "BlastP", "Protein homologies (BlastP results)") ;
      vtxtPrint (blkp, ptr) ; 
    } 
  vtxtDestroy (bfr) ;
} /* fichePRODUCTBlastPParagraph */

void ficheProductBlastPTableParagraph (vTXT blkp, GMP *gmp)
{
  vTXT bfr = vtxtCreate () ; 
 if (gmp->markup) vtxtMarkup (bfr) ;

  ficheProductBlastPTableContent (bfr, gmp, 0, 0) ; 
  if (vtxtPtr (bfr))
    {
      gmpSubSection (blkp, gmp, "Blastp", "Table of protein homologies") ; 
      vtxtPrint (blkp, vtxtPtr (bfr)) ; 
    }
  vtxtDestroy (bfr) ;
} /* ficheProductBlastPTableParagraph */

/***************************************************************************/
/***************************************************************************/
/* 5PRIME */	
void ficheMRNA5PrimeParagraphContent (vTXT blkp, GMP *gmp)
{
  BOOL isSl2s = FALSE  ; 
  int ir, firstMet, openL, Up_stop = ac_tag_int (gmp->product, "Up_stop", 0), cntSL=0 ; 
  int Length_5prime_UTR = ac_tag_int (gmp->mrna, "Length_5prime_UTR", 0) ;
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE gTranspliced_to = ac_tag_table (gmp->mrna, "Transpliced_to", h) ; 
  AC_TABLE gProducts = ac_tag_table (gmp->mrna, "Product", h) ;
  AC_OBJ oProduct ;
  vTXT bfr = vtxtCreate () ; 

  for (ir = 0 ; ir < gProducts->rows ; ir++)
    {
      oProduct = ac_table_obj (gProducts, ir, 0, h) ;
      if (oProduct && ac_has_tag (oProduct, "Best_product"))
	{
	  Length_5prime_UTR = ac_table_int (gProducts, ir, 1, 1) - 1 ;
	  break ;
	}
    }

  if (gmp->markup)
    {
      vtxtMarkup (bfr) ;

      vtxtPrint (bfr, "The sequence of the premessenger, with coding in cpper case and exons in alternate colors is ") ;
      gmpFakeObjLink (bfr, gmp, "DNA_mRNA:99999:99999", gmp->mrna, "here") ;
      vtxtBreak (bfr) ;

      vtxtPrint (bfr, "The sequence of spliced mRNA, with coding in upper case and exons in alternate colors is ") ;
      gmpFakeObjLink (bfr, gmp, "DNA_mRNA::", gmp->mrna, "here") ;
      
      vtxtBreak (bfr) ;

      vtxtPrint (bfr, "The sequence of the 2 kb upstream of the mRNA is ") ;
      {
	char buf1[256] ;
	AC_OBJ oCosmid ;
	AC_TABLE gCovers ;
	int  ma1, ma2 ;
	int u1, u2 ;

	oCosmid = ac_tag_obj (gmp->mrna, "Genomic_sequence", h) ;
	gCovers = ac_tag_table (gmp->mrna, "Covers", h) ;
	ma1 = ac_table_int (gCovers, 0, 2, 0) ;
	ma2 = ac_table_int (gCovers, 0, 4, 0) ;

	if (ma1 < ma2)
	  { u1 = ma1 - 2000 ; u2 = ma1 - 1 ; }
	else
	  { u1 = ma1 + 2000 ; u2 = ma1 + 1 ; }
	sprintf (buf1, "DNA:%d:%d:0", u1, u2) ;
	gmpFakeObjLink (bfr, gmp, buf1, oCosmid, "here") ;
      }

      vtxtBreak (bfr) ;



      vtxtBold (bfr, "The 5'UTR ") ;

      if (Length_5prime_UTR > 3)
	{
	  char buf1[256], buf2[256] ;

	  vtxtPrint (bfr, " contains ") ;
	  if (!gTranspliced_to)
	    vtxtPrint (bfr, "about ") ;
    
	  sprintf (buf1, "DNA_mRNA:%d:%d", 1, Length_5prime_UTR) ;
	  sprintf (buf2, "%d bp", Length_5prime_UTR) ; 
	  gmpFakeObjLink (bfr, gmp, buf1, gmp->mrna, buf2) ;
	}
      else if (gmp->Spc == WORM && Length_5prime_UTR < 0 && Length_5prime_UTR > -16) 
	vtxtPrintf (bfr, " contains%s %d bp stolen from the genome since the early yk library were by construction losing an average of 16 bp on the 5\' side. "
		    , ((Length_5prime_UTR && (!gTranspliced_to)) ? " about" : "")
		    , Length_5prime_UTR) ; 
    }
  
  /* do not export the 5' sequence
     if (sDna && Length_5prime_UTR)
     {
     vtxtSequence (bfr, acDna (sDna)) ;
    
     }
  */

  if (gTranspliced_to)
    {
      AC_KEYSET kRead ;

      kRead = ac_objquery_keyset (gmp->mrna
			      , messprintf ("Follow cDNA_clone  ;  Follow Read  ;  From_gene == \"%s\" " 
					  , ac_name (gmp->tg))
			      , 0) ;

      if (kRead && ac_keyset_count (kRead))
	cntSL = ficheREADListTransplitionStatement 
	  (bfr, gmp, gmp->gene, "The transpliced leader", 
	   kRead, 0, 1, "or ", &isSl2s) ; 
      ac_free (kRead) ; 
      
      if (cntSL)
	vtxtPrintf (bfr, 
		    "%s present in the mRNA in front of the sequence. "
		    , _isare (cntSL)) ; 
    }
 
  if (ac_has_tag (gmp->mrna, "Found5p"))
    { 
      if (ac_keyset_count (ac_objquery_keyset (gmp->mrna, ">product ; mRNA_5p_complete && best_product && !at_position_1", h))) 
	{
	  if (gmp->Spc == WORM)
	    {
	      if (!cntSL) vtxtPrintf (bfr, "This mRNA is not transpliced. ") ; 
	    }
	  else if (ac_has_tag (gmp->mrna, "Transpliced_to") && ac_has_tag (gmp->mrna, "Valid5p"))
	    vtxtPrintf (bfr, "It is defined by cap selected clones. ") ;
	}
    }      

  if (gmp->product &&
      !ac_has_tag (gmp->product, "Complete") &&
      !ac_has_tag (gmp->product, "NH2_complete") && 
      (openL = ac_tag_int (gmp->product, "Open_length", 0)) &&
      (firstMet = ac_tag_int (gmp->product, "First_Met", 0) + ac_tag_int (gmp->product, "First_ATG", 0)))
    {
      if (3*firstMet < openL)
	vtxtPrintf (bfr, 
		    "The first AUG occurs at bp %d. ", 
		    3*firstMet) ; 
      else if (openL > 0)
	vtxtPrintf (bfr, 
		    "There is no AUG codon in this CDS. ") ;
    }	
  
  
  if (Up_stop!=0)
    {
      if (gTranspliced_to)
	vtxtPrintf (bfr, 
		   "There is an in frame stop in the 5'UTR %d bp before the Met. "
		   , -Up_stop) ; 
      else
	vtxtPrintf (bfr, 
		   "The CDS is complete since there is an in frame stop in the 5'UTR %d bp before the Met. "
		   , -Up_stop) ; 
    }
  else 
    {   
      if (!ac_has_tag (gmp->mrna, "Found5p"))
	vtxtPrintf (bfr, 
		   "The mRNA may be incomplete at the 5' end: the frame is open. " ) ; 
    }
  
  if (gmp->mrna && !gmp->dna)
    gmp->dna = ac_obj_dna (gmp->mrna, gmp->h) ;
  
  if (gmp->dna && 
      Length_5prime_UTR >= SpcI[gmp->Spc].long5P)
    {
      vtxtDot (bfr) ;
      vtxtPrintf (bfr, 
		 "This 5\'UTR %d bp is among the 5%% longest "
		 "we have seen, maybe another protein is "
		 "encoded in this messenger. The 5' UTR contains "
		 , Length_5prime_UTR) ; 
      ficheCountATGCs (bfr, gmp->dna, Length_5prime_UTR) ; 
    }
  if (isSl2s)
    {
      ficheTGOperonParagraphContent (bfr, gmp, ficheOperonDistance, 1, 1, 1, TRUE) ; 
    }
  if (0)
    {
     
      AC_TABLE gPrev = gmp->tg ? ac_tag_table (gmp->tg, "Previous_gene_in_cis", h) : 0  ; 
      AC_OBJ prev  = gPrev ? ac_table_obj (gPrev, 0, 0, h) : 0 ;

      if (gPrev && gPrev->rows >= 2 && ac_table_int (gPrev, 0, 1, 10000) < 2000)
	vtxtPrintf (bfr, " This gene may be part of an operon with gene %s", ac_name (prev)) ;
      {
	vTXT buf ; 
	buf = vtxtCreate () ; 
	if (gmp->markup) vtxtMarkup (buf) ;
  
	ficheTGTitleStatement (buf, gmp, prev) ;
	if (vtxtPtr (buf))
	  vtxtPrintf (bfr, " (%s)", vtxtPtr (buf)) ; 
	vtxtDestroy (buf) ;
      }
      vtxtPrintf (bfr, ", endding %d bp upstream.",  ac_table_int (gPrev, 0, 1, 0)) ; 
      ac_free (prev) ;
    }

  if (vtxtPtr (bfr))
    {
      vtxtBreak (blkp) ;
      vtxtPrint (blkp, vtxtPtr (bfr)) ;
    }
  vtxtDestroy (bfr) ; 
  ac_free (h) ;
} /* ficheMRNA5PrimeParagraphContent */

/* -===================================- /
* -=  3PRIME                  =- /
* -===================================- */	

void ficheMRNA3PrimeParagraphContent (vTXT blkp, GMP *gmp)
{
  AC_TABLE gPolyA_Signal ; 
  int ir, ln = 0, PolyA_pos ; 
  const char    *PolyA_Signal ; 
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE gTranspliced_to = ac_tag_table (gmp->mrna, "Transpliced_to", h) ; 
  int Length_3prime_UTR = ac_tag_int (gmp->mrna, "Length_3prime_UTR", 0) ; 
  AC_TABLE gProducts = ac_tag_table (gmp->mrna, "Product", h) ;
  AC_TABLE gDna = ac_tag_table (gmp->mrna, "DNA", h) ;
  AC_OBJ oProduct ;

  for (ir = 0 ; ir < gDna->rows ; ir++)
    {
      ln = ac_table_int (gDna, ir, 1, 0) ;
      break ;
    }
  
  for (ir = 0 ; ir < gProducts->rows ; ir++)
    {
      oProduct = ac_table_obj (gProducts, ir, 0, h) ;
      if (oProduct && ac_has_tag (oProduct, "Best_product"))
	{
	  Length_3prime_UTR = ln - ac_table_int (gProducts, ir, 2, 1) ;
	  break ;
	}
    }
  
  if (!ac_has_tag (gmp->product, "Complete") && !ac_has_tag (gmp->product, "COOH_Complete"))
    vtxtPrintf (blkp, " is not reached, the gene is incomplete.") ; 
  else 
    { 
      char buf1[256], buf2[256] ;
      char *mrna_dna = ac_obj_dna (gmp->mrna, h) ;
      int mrnaLen = mrna_dna ? strlen (mrna_dna) : 0 ;

      vtxtPrint (blkp, " contains ") ;
      if (!gTranspliced_to)
	vtxtPrint (blkp, "about ") ;
     
      sprintf (buf1, "DNA_mRNA:%d:%d", mrnaLen - Length_3prime_UTR + 1,  mrnaLen) ;
      sprintf (buf2, "%d bp", Length_3prime_UTR) ; 
      gmpFakeObjLink (blkp, gmp, buf1, gmp->mrna, buf2) ;
      vtxtPrint (blkp, " followed by the polyA") ;

      gPolyA_Signal = ac_tag_table (gmp->mrna, "PolyA_Signal", h) ; 
      if (gPolyA_Signal)
	{
	  PolyA_Signal = ac_table_printable (gPolyA_Signal, 0, 0, "") ; 
	  
	  vtxtPrintf (blkp, ". The standard AATAAA polyadenylation signal ") ; 
	  
	  if (strstr (PolyA_Signal, "AATAA"))
	    {
	      PolyA_pos = ac_table_int (gPolyA_Signal, 0, 1, 0) ;
	      vtxtPrintf (blkp, "is seen%s %d bp before the polyA", (PolyA_pos!=0 ? " about" :""), PolyA_pos ) ; 
	    }
	  else if (strstr (PolyA_Signal, "Variant"))
	    {
	      char	bfr[1024] ; 
	      
	      vtxtPrintf (blkp, "does not occur, but the variant ") ; 
	      PolyA_Signal = ac_table_printable (gPolyA_Signal, 0, 1, "") ; 
	      PolyA_pos = ac_table_int (gPolyA_Signal, 0, 2, 0) ;
	      strcpy (bfr, PolyA_Signal) ; 
	      vtextUpperCase (bfr) ; 
	      vtxtPrintf (blkp, "%s is seen%s %d bp before the polyA addition site", 
			  bfr
			  , (PolyA_pos!=0 ? " about" :""), PolyA_pos ) ; 
	    }
	}
      else if (1 || gmp->Spc==WORM)
	vtxtPrintf (blkp, 
		    ". Neither the standard AATAAA polyadenylation signal "
		    "nor a variant is seen in the last 30 bp") ;  
      
      if (gmp->mrna && !gmp->dna)
	gmp->dna = ac_obj_dna (gmp->mrna, gmp->h) ;
      if (Length_3prime_UTR >= SpcI[gmp->Spc].long3P && gmp->dna)
	{
	  vtxtPrintf (blkp, 
		      ". This 3\'UTR %d bp is among the 5%% longest "
		      "we have seen, it may serve a regulatory function. It contains "
		      , Length_3prime_UTR+1) ; 
	  
	  ficheCountATGCs (blkp, gmp->dna + strlen (gmp->dna) - Length_3prime_UTR, Length_3prime_UTR) ; 
	}
    }
  ac_free (h) ;
} /* ficheMRNA3PrimeParagraphContent */


void ficheMRNA3PrimeParagraph (vTXT blkp, GMP *gmp)
{
  vTXT bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;
    
  ficheMRNA3PrimeParagraphContent (bfr, gmp) ;
  if (vtxtPtr (bfr))
    {
      vtxtBreak (blkp) ;
      vtxtBold (blkp, "The 3'UTR") ;
      vtxtPrintf (blkp, " %s", vtxtPtr (bfr)) ;
    }
  vtxtDestroy (bfr) ;
} /* ficheMRNA3PrimeParagraph */


/* -===================================- /
* -=  ANTISENS                        =- /
* -===================================- */	
int ficheTGAntisensParagraphContent (vTXT blkp, GMP *gmp, BOOL details)
{
  int nn = 0 ;
  AC_TABLE gAntisens_to ;
  AC_HANDLE h = ac_new_handle () ;

  if (ac_has_tag (gmp->tg, "gt_ag") && 
      (gAntisens_to = ac_tag_table (gmp->tg, "Antisens_to", h)))
    {
      AC_OBJ oTg ;
      int iOver, jj, ir ; 
      
      for (jj = ir = 0 ; ir < gAntisens_to->rows ; ir++)
	{	  
	  oTg = ac_table_obj (gAntisens_to, ir, 0, h) ; 
	  if (!strncmp (ac_name(oTg), "G_", 2) || !ac_has_tag (oTg, "gt_ag") || !ac_has_tag (oTg, "Antisens_to"))
	    continue ;
	  iOver = ac_table_int (gAntisens_to, ir, 1, 0) ; 
	  iOver = iOver > 0  ? iOver : -iOver ; 
	  if (iOver < 40) continue ; 
	  if (!jj++)
	    {
	      vtxtDot (blkp) ;
	      vtxtPrintf (blkp, "%d bp of this gene "
			  , iOver
			  ) ;
	      if (gmp->markup)
		gmpAction (blkp, gmp, gmp->gene, "vgene&v=2&S", "are antisense to spliced gene ") ;
	    }
	  else
	    {
	      vtxtComma (blkp) ;
	      vtxtPrintf (blkp, "%d to  "
			  , iOver
			  ) ;
	    }
	  nn++ ;
	  gmpObjLink (blkp, gmp,  oTg, 0) ; /* ficheTGBestNameStatement (blkp, oTg, 1) ;  */

	  if (details)
	    {
	      {
		vTXT bfr = vtxtCreate () ; 
		if (gmp->markup) vtxtMarkup (bfr) ;

		ficheTGTitleStatement (bfr, gmp, oTg) ; 
		if (vtxtPtr (bfr))
		  {
		    vtxtPrintf (blkp, " (%s)", vtxtPtr (bfr)) ; 
		  }
		vtxtDestroy (bfr) ;
	      }
	    }
	}
      
      if (nn && gmp->markup) /* not in asn dump mode */
	{ 
	  vtxtPrintf (blkp, ", raising the possibility of regulated alternate expression") ; 
	}
    } 
  ac_free (h) ;
  return nn ; 
} /* ficheTGAntisensParagraphContent */

/* -===================================- /
* -=  OPERON                  =- /
* -===================================- */	

static int ficheTGListOperonsDirect (DICT *dict, vTXT blkp, GMP *gmp, AC_OBJ tg, int operonDistance, const char *stopHere, int isTitles, int isStrands)
{
  AC_TABLE gPossible_operon ;
  AC_OBJ oOperon ;
  int cordDiference, totalCnt = 0, ir ; 
  BOOL myDict = FALSE ;
  AC_HANDLE h = ac_new_handle () ;

  if (stopHere == 0)
    stopHere = ac_name (tg) ; 
  
  if ((gPossible_operon = ac_tag_table (tg, "Possible_operon", h)))
    {
      if (!dict)
	{
	  dict = dictCreate (20) ;
	  myDict = TRUE ;
	}
      for (ir = 0 ; ir < gPossible_operon->rows ; ir++)
	{
	  cordDiference = ac_table_int (gPossible_operon, ir, 1, 0) ; 
	  
	  if (operonDistance < 0 && (cordDiference < operonDistance || cordDiference>0) )
	    continue ; 
	  else if (operonDistance>0 && (cordDiference>operonDistance || cordDiference < 0) )
	    continue ; 
	  
	  oOperon = ac_table_obj (gPossible_operon, ir, 0, h) ;
	  if (!dictAdd (dict, ac_name (oOperon), 0)) 	  /* do not recurse */
	    continue ;
	  
	  totalCnt += ficheTGListOperonsDirect (dict, blkp, gmp, oOperon, operonDistance, stopHere, isTitles, isStrands) ; 
	  
	  vtxtPrint (blkp, (vtxtLen (blkp) ? ", " : "")) ; 
	  ficheTGBestNameStatement (blkp, gmp, oOperon, 1) ;

	  if (isStrands)
	    /*vtxtPrintf (blkp, " (%d bp %s)", cordDiference < 0 ? -cordDiference : cordDiference, cordDiference < 0 ? "upstream" : "downstream") ;  */
	    vtxtPrintf (blkp, ", %d bp %s", cordDiference < 0 ? -cordDiference : cordDiference, cordDiference < 0 ? "upstream" : "downstream") ; 
	  if (gmp->Spc==WORM)
	    {
	      AC_KEYSET kRead = ac_objquery_keyset (oOperon, "Follow cDNA_clone  ;  Follow Read ", h) ; 
	      int which  ; 
	      BOOL isSl2s = 0  ; 

	      which = ac_keyset_count (kRead) ? ficheREADListTransplitionStatement (blkp, gmp, oOperon, ", transpliced to ", kRead, 0, 0, ", ", &isSl2s) : 0 ; 
	      ac_free (kRead) ; 
	      if (!which)
		vtxtPrintf (blkp, ", not known to be transpliced") ; 
	    }
	  
	  if (isTitles)
	    {
	      vTXT bfr = vtxtCreate () ;
	      if (gmp->markup) vtxtMarkup (bfr) ;

	      ficheTGTitleStatement (bfr, gmp, oOperon) ; 
	      if (vtxtPtr (bfr))
		{
		  vtxtPrintf (blkp, " (%s)", vtxtPtr (bfr)) ; 
		}
	      vtxtDestroy (bfr) ;
	    }
	  
	  totalCnt++ ; 
	}
    }
  if (myDict)
    dictDestroy (dict) ;

  ac_free (h) ;
  return totalCnt ; 
}

static int ficheTGOperonParagraphContent (vTXT blkp, GMP *gmp, int operonDistance, int isTitles, int isStrands, int previous, BOOL details)
{
  int which = 0, cntOperon = 0 ; 
  AC_KEYSET kRead ; 
  vTXT  ope = vtxtCreate () ;

  if (gmp->markup) vtxtMarkup (ope) ;

  if (gmp->Spc == WORM)
    {
      BOOL isSl2s = 0  ; 
      kRead = ac_objquery_keyset (gmp->tg, "Follow cDNA_clone  ;  Follow Read ", 0) ; 
      
      if (kRead && ac_keyset_count (kRead))
	which= ficheREADListTransplitionStatement (ope, gmp, gmp->gene, "This gene is transpliced to", kRead, 1, 1, "", &isSl2s) ; 
      ac_free (kRead) ; 
      if (which)
	{
	  vtxtPrintf (blkp, vtxtPtr (ope)) ;
	  vtxtClear (ope) ;
	}
      else
	vtxtPrintf (blkp, "This gene is not known to be transpliced") ;
    }

  if (gmp->Spc == WORM &&
      (cntOperon = ficheListOperons (0, ope, gmp, gmp->tg, ficheOperonDistance, 1, 1)) && 
      /*    if ((cntOperon = getOperonList (&ope, ficheOperonDistance, gmp->tg)) && */
      vtxtPtr (ope))
    {
      vtxtDot (blkp) ;
      if (which>1)
	vtxtPrintf (blkp, 
		    "It belongs to an ") ;
      else 
	vtxtPrintf (blkp, 
		       "It may belong to an ") ;
      
      vtxtBold (blkp, "operon") ;
      if (details)
	vtxtPrintf (blkp, 
		    " with gene %s"
		    , vtxtPtr (ope)) ; 
    } 

  vtxtDestroy (ope) ;
  return cntOperon ; 
}

/* -===================================- /
 * -=  CLONES ANOMALIES                =- /
 * -===================================- */	
static int ficheNewCloneAnomalyStatement (vTXT blkp, AC_OBJ oClone, char *prefix, BOOL showCloneName)
{
  AC_TABLE gAn, gComment ;
  AC_HANDLE h = ac_new_handle () ;
  int 	ii, ll, jj, isAny=0, ir ; 
  int messageOffset = vtxtMark (blkp) ;
  const char *ccp ;
  static struct  { char *nam, *repl ;  }chomper[]={
    {"Inverted_strand_dummy", "Submitted on the opposite strand" },   /* must be first */
    {"Internal_priming_on_A_rich", "Possibly primed on the genome, locally A rich (12 A in genome in the 18 bp downstream of last aligned base)"},
    {"Internal_priming", "Internal priming: the 3' end of the oligo dT primed clone lies in the CDS"},
    {"Mosaic", "Chimeric: part of this clone comes from one gene, the rest from another"}, 
    /*    {"Possible_mosaic", "Possibly chimeric: part of this clone comes from one gene, the rest from another"},  */
    {"Suspected_internal_deletion", "Suspected internal deletion"}, 
    {"Internal_capping", "This clone, from a cap selected library, nevertheless appears to be 5' incomplete."}, 
    {"Partly_inverted", "Part of the insert is inverted"}, 
    {"Possibly_unspliced", "Possibly unspliced"}, 
    {"Duplicate_clone", "duplicate clone from a PCR amplified library"}, 
    {"Incomplete_5p_Match", "incomplete 5' match: may be chimeric or rearranged"}, 
    {"Incomplete_3p_Match", "incomplete 3' match: may be chimeric or rearranged"}, 
    {"Length_anomaly", "measured length longer than reported mRNA, possible alternative"},  /* annot auto ? */   
    {"Suspected_genomic_contamination", "Suspected genomic contamination"},
    {"Ignore_this_clone", "This clone is problematic and was ignored"}, 
    {"Ignore_this_clone_automatic", "This clone was ignored"}, 
    {"Unspliced_non_coding", "Unspliced non coding"},
    {"Possibly_unspliced", "Possibly incompletely spliced"},
    { 0, 0 }
  } ; 							
  /*
    if Suspected_internal_deletion  >> Suspected internal deletion at base # (when # 
    available)
    Length_anomaly  >> measured length (.. kb) is longer than expected : possible 
  */

  gAn = ac_tag_table (oClone, "Anomalous_clone", h) ;
  if (1)
    {
      const char *tAnomaly = 0 ;
      
      for (ii=0, ll=0 ; (ii == 0) || (gAn && ii < gAn->rows) ; ii++)
	{
	  tAnomaly = ac_table_tag (gAn, ii, 0, 0) ; 
	  
	  for (jj = 0; (!jj || tAnomaly) && chomper[jj].nam ; jj++)
	    {
	      if (jj == 0)
		{
		  AC_KEYSET  ks = ac_objquery_keyset (oClone, ">read (forward && areverse) || (reverse && aforward)", 0) ;
		  int n1 = ks ? ac_keyset_count (ks) : 0 ;
		  ac_free (ks) ;
		  
		  if (!n1) continue ;
		}
	      else if (strcasecmp (tAnomaly, chomper[jj].nam))
		continue ;
	      if (jj== 1) /* "Internal_priming_on_A_rich" */
		{
		  ccp = ac_tag_printable (oClone, "Internal_priming_on_A_rich", 0) ;
		  ccp = messprintf ("Possibly primed on the genome, locally A rich (%s A in genome downstream of last aligned base)", ccp) ;
		}
	      else
		ccp = chomper[jj].repl ;
	      if (ll) vtxtPrintf (blkp, ", ") ; 
	      if (ll) ccp = gtLowerCleanUp (ccp) ;
	      if (showCloneName)
		vtxtPrintf (blkp, "%s (%s)", gtYkName (ac_name (oClone)), ccp) ;
	      else
		vtxtPrintf (blkp, "%s", ccp) ;
	      ll++ ; 
	      isAny=1 ; 
	      break ; 
	    }
	}
    }
  if ((gComment = ac_tag_table (oClone, "Comment", h)))
    {
      vtxtPrintf (blkp, " (", "") ; 
      for (ir=0 ; ir < gComment->rows ; ir++)
	{
	  vtxtPrintf (blkp, "%s%s", ir ? ", " :"", gtLowerCleanUp(ac_table_printable (gComment, ir, 0, ""))) ; 
	}
      vtxtPrintf (blkp, ")") ; 
      isAny=1 ; 
    }
  
  ac_free (h) ;
  return isAny ? messageOffset : 0 ;
}  /* ficheNewCloneAnomalyStatement */

/**************************************************************/

static int ficheMrnaAnomalousClones (vTXT blkp, GMP *gmp, AC_TABLE oClones)
{
  AC_OBJ oClone  ; 
  int ir, nn = 0  ; 

  if (oClones)
    for (ir = 0  ;  ir < oClones->rows  ;  ir++)
      {
	oClone = ac_table_obj (oClones, ir, 0, 0) ; 
	if (!nn) 
	  vtxtDot (blkp) ;
	nn += ficheNewCloneAnomalyStatement 
	  (blkp, oClone
	   , !nn ? "Anomalous clones: " : ", "
	   , TRUE
	   ) ; 
	ac_free (oClone) ;
      }
  return nn  ; 
}

/**************************************************************/
/**************************************************************/

/* -===================================- /
* -=  TG Structure                    =- /
* -===================================- */	


static void ficheNewGeneCountVariantsStatement (vTXT blkp, GMP *gmp)
{
  int nMrna=0,  nKantor = 0 ;
  AC_HANDLE h = ac_new_handle () ;

  if (gmp->gene)
    {
      nMrna = ac_keyset_count (ac_objquery_keyset (gmp->gene
						   , "{>Transcribed_gene;>mrna} SETELSE {>Genefinder;>predicted_mrna}"
						   , h)) ; 
      nKantor = ac_keyset_count (ac_objquery_keyset (gmp->gene
						     , " >Product ; >Kantor", h)) ;
    }
  vtxtBreak (blkp) ;
  if (gmp->tg)
    vtxtPrintf (blkp, "According to our analysis, this gene produces") ; 
  else
    vtxtPrintf (blkp, "According to the genome annotations, this gene produces") ; 
  if (nMrna > 1)
    vtxtPrintf (blkp, ", by alternative splicing, %d mRNAs", nMrna) ;
  else 
    vtxtPrintf (blkp, " a single mRNA") ; 
  
  if (nKantor)
    {
      vtxtPrintf (blkp, ", predicted to encode ") ; 
      if (nKantor > 1)
	vtxtPrintf (blkp, "%d distinct proteins", nKantor) ; 
      else 
	vtxtPrintf (blkp, "a single protein") ; 
    }
  ac_free (h) ;
}

/* -===================================- /
* -=  TG'S MRNAs and Products         =- /
* -===================================- */	


/* -===================================- /
* -=  TG'S MRNAs and Products         =- /
* -===================================- */	

/* -===================================- /
* -=  CLONE LIST TABLE                =- /
* -===================================- */	

static char* cleanVariantName (const char *gene, const char *text)
{
  static char buf[1024] ;
  const char* ccp = text ;
  char *cp, *cq ;

  if (gene && ! strncmp (gene, text, strlen (gene)))
    {
      ccp += strlen (gene) ; /* jump common part */
      if (*ccp == '.') ccp++ ;
    }
  strncpy (buf, ccp, 1023) ;
  

  /* we expect genename.aDec03   or for a product genename.b.aDecO3 */

  if (strlen(buf) > 5)
    {
      cq = buf + strlen (buf) - 5 ;
      cp = strstr (cq, "Dec03") ;
      if (cp)
	*cp = 0 ;
      cp = strstr (cq, "Sep04") ;
      if (cp)
	*cp = 0 ;
      cp = strstr (cq, "Nov04") ;
      if (cp)
	*cp = 0 ;
      cp = strstr (cq, "Jun05") ;
      if (cp)
	*cp = 0 ;
      cp = strstr (cq, "Aug05") ;
      if (cp)
	*cp = 0 ;
    }

  if (!buf[0])
    strncpy (buf, text, 1023) ;
  return buf ;
}

#ifdef JUNK

static int ficheCDNACloneListTable1 (vTXT blkp, GMP *gmp, AC_TABLE tcDNA_clone, DICT *estDict, BOOL cheatOnCloneLink, AC_OBJ oGeneCalledFor)
{
  int ii, icl, iss, tt, ig, im, ir, is,  jr, ali, nerr, leng, iLine ; 
  float	fVal ; 
  AC_TABLE gRead, gMrna, gGene, gTmp, gTissue ;
  AC_OBJ oRead, oMrna, oGene, oClone ;
  const char	*cp, *extName ;
  char	linkBuf[vONLINEMAXURLSIZE] ; 
  vTXT bfr ; 	
  Array Tbb = 0 ;
  char style = gmp->style ;
  AC_HANDLE h = ac_new_handle () ;
  Array clones = sortAliClones (gmp, tcDNA_clone, h) ;
  char *notToGenbank[]={"yke", "y1L", "ycL1", "y2L", "ycL2", "y4L", "ycL4", "yc", "yd", "ykm", "ykms", "mv2h", "dvp", "Vidal", "U454", "UTR", 0} ; 

  char *colNames[]={ "Clone", 
                     "Clone group",  /* in worm */
		     "Library", /* in worm */
		     "Tissue",   /* in human tissue/tissue */
		     "Measured<br/>size kb", 
		     "Properties", 
		     "Anomalies", 
		     "Link to Genbank<br/> or fasta sequence", 
		     "Sequence<br/>links to clone", 
		     "match over #bp<br/> (% length)", 
		     "# differences<br/> (% id)", 
		     "Gene and<br/>mRNA", 
		     "remark", 
		     0} ; 
  
  enum { COL_CLONE=0, COL_CLONEGROUP, COL_LIBRARY, COL_TISSUE, COL_MEASURE, COL_PROP, COL_ANOM,
	 COL_SEQ, COL_SEQ2, 
	 COL_MATCH, COL_DIFF, COL_GENELIST, COL_REMARK, COL_LAST} ; 
  
  bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;
  vtxtPrintf (bfr, "toto") ; /* avoid zero */
  /* initialize the table */
  Tbb = arrayCreate (20000, int) ;
  for (ii=0 ; colNames[ii] ; ii++)
    TBB (0, ii) = vtxtPrintf (bfr, "%s"ooo, colNames[ii]) ; 
  
  for (tt = 1, ii = 0 ; ii < arrayMax (clones) ; tt++, ii++)
    {
      icl = arr(clones, ii, CLALI).icl ;
      oClone = ac_table_obj (tcDNA_clone, icl, 0, h) ; /* we are looping along a acLIST */

      TBB (tt, COL_CLONE) =  vtxtPrintf (bfr, "%s", ac_name (oClone)) ; 
      /* gmpObjLink (bfr, gmp, oClone, 0) ; */
      vtxtPrintf (bfr, ""ooo) ;
      if (gmp->Spc == WORM)
	TBB (tt, COL_CLONEGROUP) = vtxtPrintf (bfr, "%s"ooo, ac_tag_printable (oClone, "Clone_group", "")) ; 
      TBB (tt, COL_LIBRARY) = vtxtPrintf (bfr, "%s"ooo, ac_tag_printable (oClone, "Library", "")) ; 
      if ((fVal = ac_table_float (ac_tag_table (oClone, "Length", h), 0, 1, 0)))
	TBB (tt, COL_MEASURE) = vtxtPrintf (bfr, "%.2f kb"ooo, fVal) ; 

         /* properties */
      {
	struct {char *lookFor, * whatToSay ; int checkExt ; }anomal[]={
	  {"Best__of", "recommended", 1}, 
	  /* 	  {"Specific__of", "specific", 1},  */
	  {"Complete_CDS__of", "covers entire CDS", 2}, 
	  {"Fully_sequenced", "fully sequenced", 0}, 
	  {"Resequence", "to resequence", 0}, 
	  {0, 0}} ; 
	
	
	for (iLine=0, is=0 ; anomal[is].lookFor ; is++)
	  {
	    if (anomal[is].checkExt)
	      {
		if ((gTmp = ac_tag_table (oClone, anomal[is].lookFor, h)))
		  {
		    
		    TBB (tt, COL_PROP) = vtxtPrintf (bfr, anomal[is].whatToSay) ; 
		    iLine++ ; 
		    extName = strrchr (ac_table_printable (gTmp, 0, 0, "."), '.') ; 
		    if (!extName)extName = "" ; else extName ++ ; 
		    if (extName[0])
		      {
			vtxtPrintf (bfr, " of variant .") ; 
			gmpObjLink (bfr, gmp, ac_table_obj (gTmp, 0, 0, h), extName) ; 
		      }
		    vtxtPrintf (bfr, ooo) ; 
		  }
	      }
	    else if (ac_has_tag (oClone, anomal[is].lookFor))
	      {
		
		TBB (tt, COL_PROP)=vtxtPrintf (bfr, anomal[is].whatToSay) ; 
		iLine++ ; 
		vtxtPrintf (bfr, ooo) ; 
	      }
	  }
	if (iLine)
	  vtxtPrintf (bfr, ooo) ;
      }

      /* anomalies */
      {
	int iAnom =  ficheNewCloneAnomalyStatement (bfr, oClone, 0, FALSE) ;
	if (iAnom)
	  { 	
	    TBB (tt, COL_ANOM) = iAnom ;
	    vtxtPrintf (bfr, ooo) ;
	    
	    iLine++ ; 
	  }
      }
      
      if ((gRead = ac_tag_table (oClone, "Read", h)))
	for (ir = jr = 0 ; ir < gRead->rows ; ir++)
	  {
	    oRead = ac_table_obj (gRead, ir, 0, h) ; 
	    if (!ac_has_tag (oRead, "From_gene"))
	      continue ;
	    if (estDict && !dictFind (estDict, ac_name(oRead), 0))
	      continue ;
	    if ((gTissue = ac_tag_table (oRead, "Tissue", h)))
	      {
		for (iss = 0 ; iss < gTissue->rows ; iss++)
		  {
		    cp = ac_table_printable (gTissue, iss, 0, 0) ;
		    if (!iss && cp && *cp)
		      TBB (tt, COL_TISSUE) = vtxtPrintf (bfr, "%s", cp) ;
		    else
		      vtxtPrintf (bfr, " %s", cp) ;
		  }
		vtxtPrintf (bfr, ooo) ;
	      }

	    if (jr++ > 0) tt++ ;

	    { 
	      if (cheatOnCloneLink) /* jan 2004, cheat so that the read points to the clone */
		{
		  TBB (tt, COL_SEQ2) = gmpFakeObjLink (bfr, gmp, "cDNA_clone", oClone, ac_name(oRead)) ;
		  vtxtPrintf (bfr, ooo) ;
		}
	      else
		{
		  /* check if read has an external or internal hot link */
		  for (is=0 ; notToGenbank[is] ; is++)
		    {
		      if (strstr (ac_tag_printable (oClone, "Library", ac_name(oClone)), notToGenbank[is]))
			break ; 
		    }
		  if (ac_has_tag (oRead, "Composite") || notToGenbank[is])
		    {
		      AC_OBJ oDna = ac_tag_obj (oRead, "DNA", h) ;
		      if (oDna)
			TBB (tt, COL_SEQ) = gmpFakeObjLink (bfr, gmp, "DNA::", oRead, ac_name(oRead)) ;
		      else
			TBB (tt, COL_SEQ) = vtxtPrint (bfr, ac_name (oRead)) ;
		      vtxtPrintf (bfr, ooo) ; 
		    }
		  else
		    {
		      const char *cleanNam = ac_name (oRead) ;
		      if (!strncmp(cleanNam, "GB:",3))
			cleanNam += 3 ;sprintf (linkBuf, GENBANK_LINK, cleanNam) ; 
		      TBB (tt, COL_SEQ) = gmpURL (bfr, gmp, linkBuf, ac_name (oRead)) ; 
		      vtxtPrintf (bfr, ooo) ; 
		    }
		}
	    }

	    if ((gGene = ac_tag_table (oRead, "From_gene", h)))
	      for (ig = 0 ; ig < gGene->rows ; ig++)
		{          /* report alignment quality */
		  if (ig > 0) tt++ ;
		  oGene = ac_table_obj (gGene, ig, 0, h) ;
		  leng = ac_table_int (gGene, ig, 1, 0) ;
		  ali = ac_table_int (gGene, ig, 3, 0) ;
		  nerr = ac_table_int (gGene, ig, 4, 0) ;
      
		  /* match over */
		  TBB (tt, COL_MATCH) =
		    vtxtPrintf (bfr, "%d bp (%.lf%%)"ooo, ali, 
				leng ? (ali < leng ? ali*100./leng : 100.0) : 0) ; 
		  
		  /* error */
		  if (ali)
		    {
		      if (nerr)
			TBB (tt, COL_DIFF)
			 =vtxtPrintf (bfr, "%d&nbsp;diff (%.1lf%%id)"ooo, nerr, 
				       leng ? 100.-nerr*100./leng : 1000 ) ; 
		      else 
			TBB (tt, COL_DIFF)
			 =vtxtPrintf (bfr, leng ? "no&nbsp;error (100%%id)"ooo : "") ; 
		    }

		  /* mRNAs */
		  iLine = 0 ;
		  if ((gMrna = ac_tag_table (oRead, "In_mRNA", h)))
		    {
		      for (im = 0 ; im < gMrna->rows ; im++)
			{
			  oMrna = ac_table_obj (gMrna, im, 0, h) ; 
			  if (!strstr (ac_name (oMrna), ac_name (oGene)))
			    continue ;
			 			 
			  extName = cleanVariantName (ac_name(oGeneCalledFor), ac_name(oMrna)) ;
			  if (!iLine)
			    {			      
			      TBB (tt, COL_GENELIST) = gmpObjLink (bfr, gmp, oMrna, extName) ; 
			    }
			  else if (!(iLine % 3))
			    {
			      vtxtPrint (bfr, ",<br/>") ;
			      gmpObjLink (bfr, gmp, oMrna, extName) ; 
			    }
			  else
			    {
			      vtxtPrintf (bfr, ", ") ;
			      gmpObjLink (bfr, gmp, oMrna, extName) ; 
			    }
			  iLine++ ;
			}
		      if (iLine) 
			vtxtPrintf (bfr, ooo) ;
		      
		    } 
		  else if ((gMrna = ac_tag_table (oClone, "From_gene", h)))
		    {
		      for (im = 0 ; im < gMrna->rows ; im++)
			{
			  oMrna = ac_table_obj (gMrna, im, 0, h) ; 
			  extName = cleanVariantName (ac_name(oGeneCalledFor), ac_name(oMrna)) ;
			  if (!iLine++)
			    TBB (tt, COL_GENELIST) = gmpObjLink (bfr, gmp, oMrna, extName) ; 
			  else
			    gmpObjLink (bfr, gmp, oMrna, extName) ; 
			}
		      if (iLine) 
			vtxtPrintf (bfr, "<br/>variant not shown"ooo) ;		      
		    }
		  else if (ac_has_tag (oRead, "Bad_quality"))
		    {
		      TBB (tt, COL_GENELIST) =
			vtxtPrintf (bfr, "bad quality trace"ooo) ; 
		    }		  
		}
	  }
      else if ((gMrna = ac_tag_table (oClone, "In_gene", h)))
	{
	  for (im = 0 ; im < gMrna->rows ; im++)
	    {
	      oMrna = ac_table_obj (gMrna, im, 0, h) ; 
	      extName = ac_name (oMrna) ;
	      if (!iLine++)
		TBB (tt, COL_GENELIST) = gmpObjLink (bfr, gmp, oMrna, extName) ; 
	      else
		gmpObjLink (bfr, gmp, oMrna, extName) ; 
	    }
	}
    }
  
  if (cheatOnCloneLink) /* jan 2004, cheat so that the read points to the clone */
    {
      if (gmp->Spc==WORM)
	fichePrintSquareTable (gmp, style, Tbb, blkp, bfr, 0, 0, 0, tt, 
			       /* COL_CLONE, COL_TISSUE, COL_SEQ,  */
			       COL_SEQ2, COL_TISSUE,
			       COL_MATCH, COL_DIFF, COL_CLONEGROUP, COL_LIBRARY, 
			       COL_MEASURE, COL_GENELIST, 
			       COL_PROP, COL_ANOM, -1) ; 
      else 
	fichePrintSquareTable (gmp, style, Tbb, blkp, bfr, 0, 0, 0,tt,  
			       /* COL_CLONE, COL_TISSUE, COL_SEQ,  */
			       COL_SEQ2, COL_TISSUE,
			       COL_MATCH, COL_DIFF/*, COL_CLONE_GROUP, COL_LIBRARY*/, 
			       COL_MEASURE, COL_GENELIST, 
			       COL_PROP, COL_ANOM, -1) ; 
    }
  else /* clone page */
    {
      if (gmp->Spc==WORM)
	fichePrintSquareTable (gmp, style, Tbb, blkp, bfr, 0, 0, 0, tt, 
			       COL_CLONE, COL_TISSUE, COL_SEQ,
			       /* COL_SEQ, COL_TISSUE, */
			       COL_MATCH, COL_DIFF, COL_CLONEGROUP, COL_LIBRARY, 
			       COL_MEASURE, COL_GENELIST, 
			       COL_PROP, COL_ANOM, -1) ; 
      else 
	fichePrintSquareTable (gmp, style, Tbb, blkp, bfr, 0, 0, 0,tt,  
			       COL_CLONE, COL_TISSUE, COL_SEQ,
			       /* COL_SEQ, COL_TISSUE, */
			       COL_MATCH, COL_DIFF/*, COL_CLONE_GROUP, COL_LIBRARY*/, 
			       COL_MEASURE, COL_GENELIST, 
			       COL_PROP, COL_ANOM, -1) ; 

    }
  arrayDestroy (Tbb) ;
  
  vtxtDestroy (bfr) ;
  ac_free (h) ;
  return 1 ; 
}  /* ficheCDNACloneListTable1 */

static int ficheCDNACloneListTable2 (vTXT blkp, GMP *gmp, AC_TABLE tcDNA_clone, Array clones, int halfMode, AC_OBJ oGeneCalledFor)
{
  int ir, icl, is, rX1, ali, leng, gerr, maxCols, maxRows, iLine, iClone, curRowClone, addRow ;
  float	fVal ; 
  AC_TABLE gTmp, gDna, gMrna, gConstructed_from ;
  AC_OBJ oRead, oMrna, oClone ;
  char	* extName, linkBuf[vONLINEMAXURLSIZE] ; 
  char style = gmp->style ;
  vTXT bfr ; 	
  Array Tbb = 0 ;
  char *notToGenbank[]={"yke", "y1L", "ycL1", "y2L", "ycL2", "y4L", "ycL4", "yc", "yd", "ykm", "ykms", "mv2h", "dvp", "Vidal", 0} ; 
  char *colNames[]={ "Clone",
		     "Clone<br/>links to details", 
		     "Library", 
		     "Tissue", 
		     "Measured<br/>size kb", 
		     "Properties", 
		     "Anomalies", 
		     "Link to Genbank<br/> or fasta sequence", 
		     "match over #bp<br/> (% length)", 
		     "# differences<br/> (% id)", 
		     "Length<br/>bp", 
		     "Strand", 
		     "Gene",
		     "mRNA", 
		     "Gene and<br/>mRNA", 
		     "coordinates<br/>on read", 
		     "coordinates<br/>on mRNA", 
		     "remark", 
		     0} ; 
  
  enum { COL_CLONE=0, COL_CLONE2, COL_LIBRARY, COL_TISSUE, COL_MEASURE, COL_PROP, COL_ANOM, COL_SEQ, 
	 COL_MATCH, COL_DIFF, COL_LEN, COL_STRAND, 
	 COL_GENE, COL_GENELIST, COL_GENELISTREGULAR, COL_CORDREAD, COL_CORDTRA, COL_REMARK, COL_LAST} ; 
  AC_HANDLE h = ac_new_handle () ;

  bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  for (maxRows=0, iClone=0 ; iClone < tcDNA_clone->rows ; iClone++)
    {
      AC_OBJ ocDNA_clone ;
      icl =  arr(clones, iClone, CLALI).icl ;
      ocDNA_clone = ac_table_obj (tcDNA_clone, icl, 0, h) ; /* we are looping along a acLIST */
  
      if ((gTmp = ac_tag_table (ocDNA_clone, "Read", h)))
	{
	  for (ir=0 ; ir < gTmp->rows ; ir++)
	    {
	      oRead = ac_table_obj (gTmp, ir, 0, h) ; 
	      
	      if ((gMrna = ac_tag_table (oRead, "In_mRNA", h)))
		{
		  maxRows+=gMrna->rows ; 
		  /* ocDNA_clone[iClone].pUser=assVoid ((gMrna->rows)) ;  */
		  /* clone on mRNA */
		  /*
		    oMrna = ac_table_obj (gMrna, 0, 0, h) ; 
		    oVal = acFindClass (ac_tag_table (oMrna, "Constructed_from", h), ac_name (oRead), 2) ; 
		    if (oVal)ocDNA_clone[iClone].pUser=assVoid (acInt ((void*) ((oVal-2)))) ; 
		  */
		}maxRows++ ; 
	    }
	}
      else maxRows++ ; 
    }
  
  Tbb = arrayCreate ((TbbDim+1)*20* (maxRows+2), int) ;
  
  vtxtPrintf (bfr, " ") ; 
  
  /* initialize the column titles */
  for (maxCols=0 ; colNames[maxCols] ; maxCols++)
    TBB (0, maxCols)=vtxtPrintf (bfr, "%s"ooo, colNames[maxCols]) ; 
  

  /* filling the table */
  for (curRowClone=0, iClone=0 ; iClone < arrayMax (clones) ; iClone++)
    {
      icl = arr(clones, iClone, CLALI).icl ;
      oClone = ac_table_obj (tcDNA_clone, icl, 0, h) ; /* we are looping along a acLIST */

      addRow=0 ; 
      
      if (halfMode)
	TBB (curRowClone+1, COL_CLONE2) = gmpObjLink (bfr, gmp, oClone, 0) ; 
      else
	TBB (curRowClone+1, COL_CLONE) = vtxtPrintf (bfr, "%s", ac_name (oClone)) ; 
      vtxtPrintf (bfr, ""ooo) ;
      TBB (curRowClone+1, COL_LIBRARY)=vtxtPrintf (bfr, "%s"ooo, ac_tag_printable (oClone, "Library", "")) ; 
      TBB (curRowClone+1, COL_TISSUE)=vtxtPrintf (bfr, "%s"ooo, ac_tag_printable (oClone, "Tissue", "")) ; 
      
      if ((fVal = ac_table_float (ac_tag_table (oClone, "Length", h), 0, 1, 0)))
	TBB (curRowClone+1, COL_MEASURE)=vtxtPrintf (bfr, "%.2f"ooo, fVal) ; 
      
      /* properties */
      {
	struct {char *lookFor, * whatToSay ; int checkExt ; }anomal[]={
	  {"Best__of", "recommended", 1}, 
	  {"Specific__of", "specific", 1}, 
	  {"Complete_CDS__of", "covers entire CDS", 2}, 
	  {"Fully_sequenced", "fully sequenced", 0}, 
	  {"Resequence", "to resequence", 0}, 
	  {0, 0}} ; 
	
	
	for (iLine=0, is=0 ; anomal[is].lookFor ; is++)
	  {
	    if (anomal[is].checkExt)
	      {
		if ((gTmp = ac_tag_table (oClone, anomal[is].lookFor, h)))
		  {
		    
		    TBB (curRowClone+1, COL_PROP)=vtxtPrintf (bfr, anomal[is].whatToSay) ; 
		    iLine++ ; 
		    extName = strrchr (ac_table_printable (gTmp, 0, 0, "."), '.') ; 
		    if (!extName)extName = "" ; else extName ++ ; 
		    if (extName[0])
		      {
			vtxtPrintf (bfr, " of variant .") ; 
			gmpObjLink (bfr, gmp, ac_table_obj (gTmp, 0, 0, h), extName) ; 
		      }
		    vtxtPrintf (bfr, ooo) ; 
		  }
	      }
	    else if (ac_has_tag (oClone, anomal[is].lookFor))
	      {
		
		TBB (curRowClone+1, COL_PROP)=vtxtPrintf (bfr, anomal[is].whatToSay) ; 
		iLine++ ; 
		vtxtPrintf (bfr, ooo) ; 
	      }
	  }
	if (addRow < iLine)
	  addRow=iLine ; 
      }
      
      
    /* anomalies */
      {
	int iAnom = ficheNewCloneAnomalyStatement (bfr, oClone, 0, FALSE) ;
	
	if (iAnom) 
	  {
	    TBB (curRowClone+1+iLine, COL_ANOM) = iAnom ;
	    vtxtPrintf (bfr, ooo) ;
	    if (addRow < iLine)addRow=iLine ; 
	  }
      }
      
      /* sequences */
      {
	AC_TABLE tbl ;

	if ((gTmp = ac_tag_table (oClone, "Read", h)))
	  {
	    for (iLine=0, ir=0 ; ir < gTmp->rows && ir < 50 ; ir++)
	      {
		oRead = ac_table_obj (gTmp, ir, 0, h) ; 
			
		/* lengths */
		gDna = ac_tag_table (oRead, "DNA", h) ;
		leng = ac_table_int (gDna, 0, 1, 0) ;
		TBB (curRowClone+1+iLine, COL_LEN)=vtxtPrintf (bfr, "%d"ooo, leng) ; 
	
		tbl = ac_tag_table (oRead, "From_gene", 0) ; 
		rX1 = ac_table_int (tbl, 0, 2, 0) ;
		ali = ac_table_int (tbl, 0, 3, 0) ; 
		gerr = ac_table_int (tbl, 0, 4, 0) ;
		ac_free (tbl) ;
		
		for (is=0 ; notToGenbank[is] ; is++)
		  {
		    if (!strcmp (notToGenbank[is], ac_tag_printable (oClone, "Library", "")))break ; 
		  }
		if (!notToGenbank[is])
		  {
		    const char *cleanNam = ac_name (oRead) ;
		    if (!strncmp(cleanNam, "GB:",3))
		      cleanNam += 3 ;
		    sprintf (linkBuf, GENBANK_LINK, cleanNam) ;
		    TBB (curRowClone+1+iLine, COL_SEQ)=gmpURL (bfr, gmp, linkBuf, ac_name (oRead)) ; 
		    vtxtPrintf (bfr, ooo) ; 
		  }
		else 
		  {
		    AC_OBJ oDna = ac_tag_obj (oRead, "DNA", h) ;
		    if (oDna)
		      TBB (curRowClone+1+iLine, COL_SEQ) = gmpFakeObjLink (bfr, gmp, "DNA::", oRead, ac_name(oRead)) ;
		    else
		      TBB (curRowClone+1+iLine, COL_SEQ) = vtxtPrint (bfr, ac_name (oRead)) ;
		    vtxtPrintf (bfr, ooo) ; 
		  }		
		
		/* strand */
		strcpy (linkBuf, ac_tag_printable (oRead, "Strand", "")) ; 
		vtextLowerCase (linkBuf) ; 
		TBB (curRowClone+1+iLine, COL_STRAND)=vtxtPrintf (bfr, "%s"ooo, linkBuf) ; 
		
		/* mRNAs */
		is=1 ; 
		if ((gMrna = ac_tag_table (oRead, "In_mRNA", h)))
		  {
		    int ll=0 ; AC_OBJ curmRNA ; 
		    
		    for (is=0 ; is < gMrna->rows ; is++)
		      {
			oMrna = ac_table_obj (gMrna, is, 0, h) ; 
			
			{
			  AC_OBJ oTg, oGene = 0 ;

			  if (ac_has_tag (oMrna, "From_gene"))
			    {
			      oTg = ac_tag_obj (oMrna, "From_gene", h) ;
			      oGene = oTg ? ac_tag_obj (oTg, "Gene", h) : 0 ;
			      ac_free (oTg) ;
			    }
			  else if (ac_has_tag (oMrna, "From_prediction"))
			    {
			      oTg = ac_tag_obj (oMrna, "From_prediction", h) ;
			      oGene = oTg ? ac_tag_obj (oTg, "Model_of_gene", h) : 0 ;
			      ac_free (oTg) ;
			    }
			  if (oGene)
			    {
			      TBB (curRowClone+1+iLine+is, COL_GENE) = gmpObjLink (bfr, gmp, oGene, 0) ; 
			      vtxtPrintf (bfr, ooo) ; 
			      ac_free (oGene) ;
			    }
			}

			if (is==ll)
			  {
			    AC_OBJ myTg = 0 ;
			    extName = cleanVariantName (ac_name(oGeneCalledFor), ac_name(oMrna)) ;

			    TBB (curRowClone+1+iLine+is, COL_GENELIST)
			      = gmpObjLink (bfr, gmp, oMrna, extName) ; 

			    for (ll++ ; is+ll < gMrna->rows ; ll++)
			      { /* mieg, was ll < gMrna->rows  */
				ac_free (myTg) ;
				curmRNA = ac_table_obj (gMrna, is+ll, 0, h) ; 
				if (oGeneCalledFor && 
				    (myTg = ac_tag_obj (curmRNA, "From_gene", h)) &&
				    strcmp (ac_name (oGeneCalledFor), ac_tag_printable (myTg, "Gene", "")))
				  break ; 
				extName = cleanVariantName (ac_name(oGeneCalledFor), ac_name(curmRNA)) ;				vtxtPrintf (bfr, ", ") ; 
				gmpObjLink (bfr, gmp, curmRNA, extName) ; 
			      }
			    ac_free (myTg) ;
			    vtxtPrintf (bfr, ooo) ; 
			  }
			
			
			extName = cleanVariantName (ac_name(oGeneCalledFor), ac_name(oMrna)) ;
			TBB (curRowClone+1+iLine+is, COL_GENELISTREGULAR)
			  =gmpObjLink (bfr, gmp, oMrna, extName) ; 
			vtxtPrintf (bfr, ooo) ; 
			
			/*
			  TBB (curRowClone+1+iLine+is, COL_GENELISTREGULAR)=gmpObjLink (bfr, gmp, oMrna, 0) ; 
			  vtxtPrintf (bfr, ooo) ; 
			*/
			
			if ((gConstructed_from = ac_tag_table (oMrna, "Constructed_from", h)))
			  {
			    int b1, jr ;
			    
			    for (jr = 0 ; jr < gConstructed_from->rows ; jr++)
			      if (!strcmp (ac_name(oRead), ac_table_printable (gConstructed_from, jr, 2, "")))
				{
				  
				  /* coordinates on read */
				  TBB (curRowClone+1+iLine+is, COL_CORDREAD)
				    = vtxtPrintf (bfr, "%d to %d"ooo
						  , ac_table_int (gConstructed_from, jr, 3, 0)
						  , ac_table_int (gConstructed_from, jr, 4, 0)) ;
				  
				  /* coordinates on mRNA */
				  b1 = ac_table_int (gConstructed_from, jr, 0, 0) ;
				  if (b1 < -12)  b1 = -12 ;
				  TBB (curRowClone+1+iLine+is, COL_CORDTRA)
				    = vtxtPrintf (bfr, "%d to %d"ooo
						  , b1
						  , ac_table_int (gConstructed_from, jr, 1, 0)) ;
				}
			    
			  }
		      }
		    is++ ; 
		  }
		else if (ac_has_tag (oRead, "Bad_quality"))
		  {
		    TBB (curRowClone+1+iLine, COL_GENELIST) =
		      vtxtPrintf (bfr, "bad quality trace"ooo) ; 
		  }
		
		/* match over */
		TBB (curRowClone+1+iLine, COL_MATCH) =
		  vtxtPrintf (bfr, "%d bp (%.lf%%)"ooo, ali,
			      leng ? (ali < leng ? ali*100./leng : 100.0) : 0) ; 
		
		/* error */
		if (ali)
		  {
		    if (gerr)
		      TBB (curRowClone+1+iLine, COL_DIFF)
			=vtxtPrintf (bfr, "%d&nbsp;err (%.1lf%% id)"ooo, gerr, 
				     leng ? 100.-gerr*100./leng : 1000 ) ; 
		    else 
		      TBB (curRowClone+1+iLine, COL_DIFF)
			=vtxtPrintf (bfr, leng ? "no&nbsp;error (100%% id)"ooo : "") ; 
		  }
		
		/* remarks */
		if (ac_has_tag (oRead, "Directed_sequencing"))
		  {
		    const char *txt ; 
		    
		    if ((txt = ac_tag_printable (oRead, "Transpliced_to", "")) && (strstr (txt, "SL")==txt) )
		      {
			if (txt[2]=='0')TBB (curRowClone+1+iLine, 15)=vtxtPrintf (bfr, "capped ") ; 
			else TBB (curRowClone+1+iLine, COL_REMARK)=vtxtPrintf (bfr, "SL%s ", txt+2) ; 
		      }
		    if (ac_has_tag (oRead, "PolyA_after_base"))
		      {
			if (TBB (curRowClone+1+iLine, COL_REMARK))TBB (curRowClone+1+iLine, 15)=vtxtPrintf (bfr, "AAA ") ; 
			else vtxtPrintf (bfr, "AAA ") ; 
		      }
		    if (TBB (curRowClone+1+iLine, COL_REMARK))
		      vtxtPrintf (bfr, ooo) ; 
		  }
		
		iLine+=is ; 
	      }
	  }
	if (addRow < iLine)addRow=iLine ; 
      }
      
      
      curRowClone+=addRow ; 
    }

  maxRows=1+curRowClone ; 
  
  if (halfMode==0)
    {
      GMP *gmp2 = gmpCreate (gmp->db, 0, 0, 0, 0, 0, style, 'z') ;
      gmpSection (blkp, gmp2, "SeqTab", "Sequences") ; 
      vtxtBreak (blkp) ;
      fichePrintSquareTable (gmp, style, Tbb, blkp, bfr, 0, 0, 0, maxRows, 
			     COL_SEQ, 
			     COL_GENE, 
			     COL_GENELISTREGULAR, COL_CORDREAD, COL_CORDTRA, 
			     COL_REMARK, COL_LEN, COL_STRAND, -1) ; 
      vtxtBreak (blkp) ;
      gmpDestroy (gmp2) ;
    }
  
  arrayDestroy (Tbb) ;
  
  vtxtDestroy (bfr) ;
  ac_free (h) ;

  return 1 ; 
} /* ficheCDNACloneListTable2 */



static int clAliOrder (const void *a, const void *b)
{
  const CLALI *aa = (const CLALI *)a , *bb = (const CLALI *)b ;
  if (aa->ali != bb->ali) return bb->ali - aa->ali ;
  return aa->icl - bb->icl ;
}


static Array sortAliClones (GMP *gmp, AC_TABLE tcDNA_clone, AC_HANDLE h1)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *dict = dictHandleCreate (3000, h) ;
  int ii = 0, ir, jr, n ;
  const char *ccp  ;
  AC_TABLE oTmp ;
  AC_OBJ oClone, oRead ;
  CLALI *up ;
  Array clones = arrayHandleCreate (300, CLALI, h) ;
  AC_KEYSET nms = ac_objquery_keyset (gmp->gene, ">read ; Ref_seq ; >cdna_clone", h) ;
  AC_KEYSET mrnaCoveringClones = ac_objquery_keyset (gmp->tg, ">Mrna ; >Mrna_covered_by ; >cdna_clone", h) ;
  AC_KEYSET cdsCoveringClones = ac_objquery_keyset (gmp->tg, ">Mrna ; >CDS_covered_by ; >cdna_clone", h) ;
  AC_KEYSET tilingClones = ac_objquery_keyset (gmp->tg, ">read ; mRNA_tiling ; >cdna_clone", h) ;
  for (ir = ii = 0 ; ir < tcDNA_clone->rows ; ii++, ir++)
    {
      oClone = ac_table_obj (tcDNA_clone, ir, 0, 0) ;
      up = arrayp (clones, ii, CLALI) ;
      up->icl = ii ;
      if (ac_keyset_contains (nms, oClone)) up->isNm = 1 ;
      if (ac_keyset_contains (mrnaCoveringClones, oClone)) up->isTiling = 3 ;
      else if (ac_keyset_contains (cdsCoveringClones, oClone)) up->isTiling = 2 ;
      else  if (ac_keyset_contains (tilingClones, oClone)) up->isTiling = 1 ;

      if ((oTmp = ac_tag_table (oClone, "In_mRNA", 0)))
	for (ir = 0 ; ir < oTmp->rows ; ir++)
	  {
	    ccp = ac_table_printable (oTmp, ir, 0, "") ;
	    dictAdd (dict, ccp, &n) ;
	    if (!up->mrna || lexstrcmp (ccp, up->mrna))
	      up->mrna = dictName (dict, n) ;
	  }
      ac_free (oTmp) ;
      if ((ccp = ac_tag_printable (oClone, "HINV_libs", 0)))
	{
	  dictAdd (dict, ccp, &n) ;
	  up->hinv_lib = dictName (dict, n) ;
	}
      else
	up->hinv_lib = "zzzz" ;
      
      up->ali = 0 ; /* may be we should order according to the length */
      if ((oTmp = ac_tag_table (oClone, "Read", 0))) 
	for (jr = 0 ; jr < oTmp->rows ; jr++)
	  {
	    oRead = ac_table_obj (oTmp, jr, 0, 0) ;
	    if ((ccp = ac_tag_printable (oRead, "Tissue", 0)))
	      {
		dictAdd (dict, ccp, &n) ;
		up->tissue = dictName (dict, n) ;
	      }
	    ac_free (oRead) ;
	  } 
      ac_free (oTmp) ;
      ac_free (oClone) ;
    }
  arraySort (clones, clAliOrder) ;

  ac_free (h) ;
  return clones ;
}
#endif
/**************************************************************************************/
/**************************************************************************************/

static void markupLinkPubmed (vTXT blkp, GMP *gmp, const char *id, const char *txt)
{
  gmpURL (blkp, gmp
	  , messprintf("https://www.ncbi.nlm.nih.gov/corehtml/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=%s&dopt=AbstractPlus", id)
	  , txt ? txt : id) ;
  
}

static void markupLinkAvery (vTXT blkp, GMP *gmp, const char *id)
{ 
  gmpURL (blkp, gmp, messprintf("http://elegans.swmed.edu/wli/%s", id), id) ;				
}

void ficheNewGeneBiblioChapter  (vTXT blkp, GMP *gmp, int gRif)
{
  AC_OBJ oPap ;
  AC_TABLE gPap ;
  int ir, ir1 = 0, ir2 = 0, j, nn, nRif = 0, nPubMed ;
  const char *pm, *ptr, *cit, *ccp ;
  AC_HANDLE h = ac_new_handle () ;
  vTXT bfr1 = vtxtHandleCreate (h) ;
  vTXT bfr2 = vtxtHandleCreate (h) ;
  vTXT bfr3 = vtxtHandleCreate (h) ;

  if (gmp->markup) vtxtMarkup (bfr1) ;   
  if (gmp->markup) vtxtMarkup (bfr2) ;   
  if (gmp->markup) vtxtMarkup (bfr3) ;   

  gPap = ac_tag_table (gmp->gene, "Reference", h) ;
  nRif = ac_keyset_count (ac_objquery_keyset (gmp->gene, ">Reference ; Gene_rif", h)) ;
  if (!gPap || (gRif && !nRif))
    { ac_free (h) ; return ; }
 
  /* report the global link to pubmed */

  vtxtClear (bfr1) ;
  /*   vtxtPrint (bfr1, PUBMED_LINK) ;  no %s included */
  vtxtPrintf (bfr1, PUBMED_GENEID_LINK, ac_tag_printable (gmp->gene, "GeneId", "toto" )) ; /* no %s included */
  if (0) 
    for (ir = gPap->rows - 1, ir1 = ir2 = 0, j = 0 ; ir >= 0 ; ir--)
      {
	ccp = ac_table_printable (gPap, ir, 0, 0) ;
	
	if (ccp && !strncasecmp (ccp, "pmp", 2))
	  {
	    if (ir1++) 
	      { if (j < 300) vtxtPrint (bfr1, ",") ; }
	    if (j++ < 300) vtxtPrint (bfr1, ccp+2) ; 
	  }
	else
	ir2++ ;
      }
  if (0) nPubMed = ir1 ;
  ir2 = 0 ;
  nPubMed = gPap->rows ; ;
  if (gmp->Spc == WORM && ir1)
    {
      vtxtPrintf (bfr2, "Please see %s ",  _theese(ir1)) ;
      gmpURL (bfr2, gmp, vtxtPtr (bfr1)
	      ,messprintf ("%d article%s in PubMed", ir1, _multi(ir1))
	      ) ;
      vtxtBreak (bfr2) ;
    }
  if (gRif && nRif)
    {
      vtxtPrintf (bfr2, "NCBI reference into function (RIF) annotation%s"
		  , _multi (nRif)) ;
    }
  else if (!gRif && ir2)
    vtxtPrintf (bfr2, "In addition we found %d papers for which we do not have a PubMed identifier", ir2) ;
  else
    goto done ;

  if (0)
    {
      vtxtPrintf (bfr2,"<ul>") ;  
      for (ir = gPap->rows - 1, oPap = 0 ; ir >= 0 ; ac_free (oPap), ir--)
	{
	  cit = 0 ;
	  oPap = ac_table_obj (gPap, ir, 0, h) ; 
	  if (gRif && ac_has_tag (oPap, "Gene_rif") &&
	      (pm = ac_tag_printable (oPap, "PMID", 0)) &&
	      (cit = ac_tag_printable (oPap,"Citation", 0)))
	    {
	      vtxtPrintf (bfr2,"\n<li>") ;
	      markupLinkPubmed (bfr2, gmp, pm, cit) ; 
	      if ((ptr = ac_tag_printable (oPap,"Title", 0)))
		vtxtPrintf (bfr2, " %s", ptr) ;
	      vtxtDot (bfr2) ;
	      if ((ptr = ac_tag_printable (oPap, "Gene_rif", 0)))
		{	      
		  vtxtPrintf (bfr2, "\n     <span class='rif'>") ;
		  vtxtPrintf (bfr2, " %s", ptr) ;    
		  vtxtPrintf (bfr2, "</span>\n") ; 
		}
	    }
	  else if (!gRif && !ac_has_tag (oPap, "PMID"))
	    {
	      vtxtPrintf (bfr2,"\n<li>") ;
	      
	      if (gmp->Spc == WORM && (ac_name (oPap)[0] == '['))
		{ 
		  markupLinkAvery (bfr2, gmp, ac_tag_printable (oPap, "LeonID", ac_name (oPap))) ;
		  if ((ptr = ac_tag_printable (oPap,"Title", 0)))
		    vtxtPrintf (bfr2, " %s", ptr) ;
		  continue ;
		}
	      
	      vtxtPrint (bfr2, ac_name (oPap)) ; 
	      if ((ptr = ac_tag_printable (oPap,"Title", 0)))
		vtxtPrintf (bfr2, " %s", ptr) ;
	      
	      if ((ptr = ac_tag_printable (oPap,"Citation", 0))) /* cit but no pmid, unlikely */
		vtxtPrintf (bfr2, " %s", ptr) ;
	      else
		{
		  if ((ptr = ac_tag_printable (oPap,"Journal", 0)))
		    vtxtPrintf (bfr2, " %s", ptr) ;
		  if ((nn = ac_tag_int (oPap,"Year", 0)))
		    vtxtPrintf (bfr2, " %d", nn) ;
		  if ((ptr = ac_tag_printable (oPap,"Volume", 0)))
		    vtxtPrintf (bfr2, " %s", ptr) ;
		  if (!strncmp(ac_name (oPap),"[cgc",4) && 
		      (ptr = ac_tag_printable (oPap,"Page", 0)))
		    vtxtPrintf (bfr2, " %s", ptr) ;
		}	  
	    }
	}
      vtxtPrintf (bfr2,"</ul>") ;
    }

 done:
  if (gRif)
    gmpChapter (blkp, gmp, "*BiblioRif"
		, "Bibliography and NLM 'reference into function' annotations"
		) ;
  else
    {
      if (gmp->Spc == WORM || nPubMed < 1)
	gmpChapter (blkp, gmp, "*Biblio", "Bibliography") ;
      else
	{
	  vtxtPrint (bfr3, "Bibliography:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp ") ;
	  gmpURL (bfr3, gmp, vtxtPtr (bfr1) 
		  ,messprintf ("<font color='blue'>%d article%s in PubMed</font>", nPubMed, _multi(nPubMed))
		  ) ;
	  gmpChapter (blkp, gmp, "*Biblio", vtxtPtr (bfr3)) ;
	}
    }
  vtxtPrint (blkp, vtxtPtr (bfr2)) ;
  
  if (gRif)
    gmpChapterClose (blkp, gmp, "BiblioRif", TRUE) ;
  else
    gmpChapterClose (blkp, gmp, "Biblio", TRUE) ;
  
  ac_free (h) ;
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
/*
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/  PARTS FUNCTIONS
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

/*****************************************************************************/

char *fichePrimersParagraphContent (vTXT blkp, GMP *gmp)
{
  AC_TABLE oPrimers = 0 ;
  AC_OBJ oProd = 0 ;
  int  nn1 = 0,  nbp = 0 ;
  AC_HANDLE h = ac_new_handle () ;

  if (gmp->product)
    {
      oPrimers = ac_tag_table (gmp->product, "Primers", h) ;
      nbp = ac_tag_int (oProd, "Coding_length", -1) + 3 ;
    }

  if (oPrimers && oPrimers->rows >= 2 &&  oPrimers->cols >= 2)
    {
      const char *p5, *p3, *pc5, *pc3 ;
      p5 = strnew (ac_table_printable (oPrimers,0,0,0), h) ;
      pc5 = strnew (ac_table_printable (oPrimers,0,1,0), h) ;
      p3 = strnew (ac_table_printable (oPrimers,1,0,0), h) ;
      pc3 = strnew (ac_table_printable (oPrimers,1,1,0), h) ;

      if (pc3 && pc5)
	{
	  pc3 = strstr (pc3, "T=") ;
	  pc5 = strstr (pc5, "T=") ;
	}
      if (pc3 && pc5)
	{
	  nn1++ ;
	  vtxtPrint (blkp , "Primers to amplify the predicted CDS (" ) ;
	  
	  if (nbp > 12)
	    vtxtPrintf (blkp
			, "%d bp, "
			, nbp) ;
	  vtxtPrintf (blkp
		      ,"Stop included): %s (%s, %s (%s"
		      , p5, pc5, p3, pc3
		      ) ;
	}	    
    }
  ac_free (h) ;

  return nn1 ? vtxtPtr (blkp) : 0 ;
} /* fichePrimersParagraphContent */

/*****************************************************************************/

static void fichePrimersParagraph (vTXT blkp, GMP *gmp)
{
  char *ptr ;
  vTXT bfr = vtxtCreate () ;
  if (gmp->markup) vtxtMarkup (bfr) ;

  if ((ptr = fichePrimersParagraphContent (bfr, gmp)))
    {
      gmpSection (blkp, gmp, "Primers", "Primers") ;
      vtxtPrint (blkp, ptr) ;
      vtxtBreak (blkp) ;
    }
  vtxtDestroy (bfr) ;
} /* fichePrimersParagraph */

/*****************************************************************************/

int ficheCloneParagraphContent (vTXT blkp, GMP *gmp)
{
  int ir, mypos, nClo = 0 ;
  char *acc ;
  AC_OBJ oBestClone ;
  AC_TABLE gCompleteClones, gcDNA_clone ;
  int dummy, nn1 = 0, nn2 = 0, nn3 = 0 ;
  AC_HANDLE h = ac_new_handle () ;

  DICT *dict = dictCreate (200) ;
  
  gcDNA_clone = ac_tag_table (gmp->mrna, "cDNA_clone", h) ;
  if (gcDNA_clone && gcDNA_clone->rows == 1)
    nClo = 1 ;
  else if (gcDNA_clone && gcDNA_clone->rows > 1)
    {
      AC_OBJ oClone = ac_table_obj (gcDNA_clone, 0, 0, h) ;
      /*       AC_TABLE gDNA = ac_tag_table (gmp->mrna, "DNA", h) ; */

      gCompleteClones = ac_tag_table (gmp->mrna, "Complete_CDS_clone", h) ;
      mypos = vtxtLen (blkp) ; nn2 = 0 ;
      if (gCompleteClones)
	for (ir = 0 ; nn1 < 1000 && ir < gCompleteClones->rows ;ir++)
	  {
	    oClone = ac_table_obj (gCompleteClones, ir, 0, h) ;
	    if (!dictAdd (dict, ac_name (oClone), &dummy))
	      continue ;	  
	    if (!nn2++) 
	      {
		vtxtBreak (blkp) ;
		mypos = vtxtPrintf (blkp, "Complete CDS clone%s: "
				    , gCompleteClones->rows > 1 ? "s" : "") ;
	      }
	    else  
	      vtxtPrintf (blkp, ", ") ;
	    nn1++ ;
	    acc = gtYkName (ac_name (oClone)) ;
	    gmpObjLink (blkp, gmp,  oClone, acc) ;

	    if (strlen(vtxtPtr (blkp)) - mypos > 900)
	      {
		vtxtPrintf (blkp, "...") ;
		nn1 = 1000 ;
	      }
	  }
      
      gCompleteClones = 0 ; /* ac_tag_table (gmp->mrna, "Specific_clone", h) ; */
      mypos = vtxtLen (blkp) ; nn2 = 0 ;
      if (gCompleteClones)
	for (ir = 0 ; nn1 < 1000 && ir < gCompleteClones->rows ;ir++)
	  {
	    oClone = ac_table_obj (gCompleteClones, ir, 0, h) ;
	    if (!dictAdd (dict, ac_name (oClone), &dummy))
	      continue ;	  
	    if (!nn2++) 
	      {
		vtxtBreak (blkp) ;
		mypos = vtxtPrintf (blkp, "Clone%s specific of this variant %s "
				    ,gCompleteClones->rows > 1 ? "s" : ""
				    ,gCompleteClones->rows > 1 ? "are" : "is") ;
	      }
	    else  
	      vtxtPrintf (blkp, ", ") ;
	    nn1++ ;
	    acc = gtYkName (ac_name (oClone)) ;
	    gmpObjLink (blkp, gmp,  oClone, acc) ;

	    if (strlen(vtxtPtr (blkp)) - mypos > 900)
	      {
		vtxtPrintf (blkp, "...") ;
		nn1 = 1000 ;
	      }
	  }
      
      oBestClone = ac_tag_obj (gmp->mrna, "Best_available_clone", h) ;
      if (oBestClone)
	{
	  vtxtBreak (blkp) ;
	  if (gmp->Spc == WORM)
	    {
	      vtxtPrintf (blkp, "Recommended clone (from the Kohara collection): ") ;
	    }
	  else
	    vtxtPrintf (blkp, "Recommended clone: ") ;
	  acc = gtYkName (ac_name (oBestClone)) ;
	  gmpObjLink (blkp, gmp,  oBestClone, acc) ;

	  dictAdd (dict, ac_name (oBestClone), &dummy) ;
	  nn1++ ;
	}
      
      gcDNA_clone = ac_tag_table (gmp->mrna, "cDNA_clone", h) ;
      nClo = gcDNA_clone->rows ;
      mypos = vtxtLen (blkp) ;
      if (0 &&
	  gcDNA_clone)
	for (ir = 0 ; nn1 < 1000 && ir < gcDNA_clone->rows ;ir++)
	  {
	    oClone = ac_table_obj (gcDNA_clone, ir, 0, h) ;
	    
	    if (!dictAdd (dict, ac_name (oClone), &dummy))
	      continue ;	  
	    
	    if (!nn3++) 
	      { 
		vtxtBreak (blkp) ;
		mypos = vtxtPrintf (blkp, "Other clone%s: ", ir < gcDNA_clone->rows - 1 ? "(s)" : "" ) ;
	      }
	    else if (nn1++) vtxtPrintf (blkp, ", ") ;
	    acc = gtYkName (ac_name (oClone)) ;
	    gmpObjLink (blkp, gmp,  oClone, acc) ;

	    if (strlen(vtxtPtr (blkp)) - mypos > 900)
	      {
		vtxtPrintf (blkp, "..") ;
		nn1 = 1000 ;
	      }
	  }
      if (0)
	{
	  vtxtBreak (blkp) ;
	  ficheMrnaAnomalousClones (blkp, gmp, gcDNA_clone) ;
	}
    }
  if ((gmp->style == 's' || gmp->style == 'r') &&
      vtxtPtr (blkp) && gmp->Spc == WORM)
    vtxtPrintf (blkp, " for edited clone sequences see www.wormgenes.org") ;
  dictDestroy (dict) ;
  ac_free (h) ;

  return nClo ;
} /* ficheCloneParagraphContent */

/*****************************************************************************/

typedef struct clo2libstruct { AC_OBJ clo, lib ; int idx ; } C2LIB ;
static int C2libOrder (const void *a, const void *b)
{
  const C2LIB *va = (const C2LIB *)a , *vb = (const C2LIB *)b ;
  const char *la, *lb ;

  if (va->idx != vb->idx)
    {
      return va->idx - vb->idx ;
    }
  else
    {
      la = ac_name (va->lib) ; lb = ac_name (vb->lib) ;
      return strcmp (la, lb) ;
    }
}

/* declare the lib here in the order you wish them in the output */
static DICT *ficheCloneLibDict (void)
{
  static DICT *dict = 0 ;

  if (dict)
    return dict ;
  
  dict = dictCreate (12) ;
  dictAdd (dict, "yke", 0) ;
  dictAdd (dict, "cee", 0) ;
  dictAdd (dict, "ycL1", 0) ;
  dictAdd (dict, "ycL2", 0) ;
  dictAdd (dict, "ycL4", 0) ;
  
  dictAdd (dict, "yc", 0) ;

  dictAdd (dict, "cm", 0) ;
  dictAdd (dict, "cem", 0) ;
  
  dictAdd (dict, "Barstead", 0) ;
  dictAdd (dict, "Kim", 0) ;
  dictAdd (dict, "Okkema", 0) ;
  dictAdd (dict, "cem", 0) ; 
  
  dictAdd (dict, "ykm", 0) ;
  dictAdd (dict, "ykms", 0) ;
  
  dictAdd (dict, "sperm", 0) ;

  dictAdd (dict, "ovary", 0) ;
  dictAdd (dict, "mv2h", 0) ;
  dictAdd (dict, "rare", 0) ; 
  dictAdd (dict, "yd", 0) ;

  dictAdd (dict, "cadmium", 0) ;

  return dict ;
}

static BOOL ficheCloneLibParagraphContent (vTXT blkp, GMP *gmp)
{
  int ir, jr, idx ;
  BOOL ok = FALSE ;
  AC_TABLE oTmp = 0 ;
  AC_OBJ oClone, oLib ;
  C2LIB *c2lib ;
  Array C2lib = 0 ;
  DICT *dict = ficheCloneLibDict () ; /* do not destroy */
  vTXT bfr = 0 ;
  AC_HANDLE h ;

  if (gmp->style == 'r') /* july 16, 2003 */
    return ok ;
  if (! gmp->tg)
    return ok ;
    
  h = ac_new_handle () ;
  bfr = vtxtHandleCreate (h) ;
  if (gmp->markup) vtxtMarkup (bfr) ;

  if (gmp->Spc == HUMAN)
    {
      ok = TRUE ;
      goto done ;
    }
  /* export only if the clones can be attributed to this variant */
 
  if (!gmp->pg && gmp->mrna)
    oTmp = ac_tag_table (gmp->mrna, "cDNA_clone", h) ;
  else if (gmp->gene)
    {
      oTmp = ac_tag_table (gmp->gene, "genefinder", h) ;
      if (gmp->view == 'm' && oTmp->rows > 1)
	return ok  ;
      oTmp = ac_tag_table (gmp->gene, "Has_cDNA_clone", h) ;
    }
    
  if (!oTmp || !oTmp->rows)
    return ok ;

  if (oTmp && oTmp->rows >= 1)
    {
      C2lib = arrayCreate (oTmp->rows, C2LIB) ;
      for (ir = 0 ; ir < oTmp->rows ; ir++)
	{
	  oClone = ac_table_obj (oTmp, ir, 0, h) ;
	  oLib = oClone ? ac_tag_obj (oClone, "Library", h) : 0 ;
	  dictAdd (dict, ac_name (oLib), &idx) ;
	  c2lib = arrayp (C2lib, ir, C2LIB) ;
	   c2lib->clo = oClone ;
	   c2lib->lib = oLib ;
	   c2lib->idx = idx ;
	}
      arraySort (C2lib, C2libOrder) ;
      for (ir = jr = 0 ; ir < arrayMax (C2lib) ; ir++)
	{
	  c2lib = arrp (C2lib, ir, C2LIB) ;
	  if ( ! c2lib->lib) continue ;
	  oClone = c2lib->clo ;
	  oLib = c2lib->lib ;
	  if (!ir || 
	      (
	       ac_name(oLib) != ac_name((c2lib-1)->lib) &&
	       ( strcasecmp (ac_name(oLib) , "ykms") || strcasecmp (ac_name((c2lib-1)->lib) , "ykm"))
	       )) /* mix up ykm ykms in the export */
	    {
	      jr = 0 ;
	      vtxtPrintf (bfr, "%s%s:"
			  , ir > 0 ? ";\n" : ""
			  , ac_tag_printable (oLib, "Title", ac_name (oLib))) ;
	    }
	  vtxtPrintf (bfr, "%s%s"
		      , jr++ > 0 ? ", " : " "
		      , gtYkName (ac_name (oClone))) ;
	}
    }

  if (vtxtPtr (bfr))
    {
      ok = TRUE ;
      vtxtBreak (blkp) ;
      vtxtEmptyLine (blkp, 1) ;
      vtxtBold (blkp, "Libraries: ") ;
      vtxtPrint (blkp, vtxtPtr (bfr)) ;
    }

 done:
  vtxtDestroy (bfr) ;
  arrayDestroy (C2lib) ;
  ac_free (h) ;

  return ok ;
} /* ficheCloneLibParagraphContent */

/***************************************************************************************/
/**************************  ficheNewGeneTitleParagraph ********************************/
/***************************************************************************************/

void ficheNewGeneTitleParagraph (vTXT blkp, GMP *gmp)
{
  char *ptr ;
  vTXT bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  vtxtHr (blkp, 0, 0) ;
  if ((ptr = gtGeneTitle (bfr, gmp, TRUE))) /* contains a dot. nasty when bold */
    {	
      vtxtBreak (blkp) ;
      vtxtBold (blkp, ptr) ;
      vtxtPrintf (blkp,"\n") ; /* avoid a second dot */
    }
  vtxtHr (blkp, 0, 0) ;
    
  vtxtDestroy (bfr) ; 
  return ;
}

/***************************************************************************************/
/**************************  ficheNewGeneSummaryChapter ********************************/
/***************************************************************************************/

char *ficheNewGeneAceViewSummary (vTXT blkp, GMP *gmp)
{
  BOOL isCloud ;
  
  if (gmp->tg) ficheNewGeneComplexLocusStatement (blkp, gmp) ;
  ficheNewGeneExpressionTissueStatement (blkp, gmp, TRUE) ;
  isCloud = ficheNewGeneIntronsStatement (blkp, gmp, FALSE, TRUE) ;
  if (!isCloud) ficheNewGeneAltVariantStatement (blkp, gmp, TRUE) ;
  ficheNewGeneAltFeatureStatement (blkp, gmp, TRUE, TRUE) ;
  if (gmp->tg) 
    {
      ficheTGAntisensParagraphContent (blkp, gmp, FALSE) ; 
      ficheNewGeneNmdStatement (blkp, gmp, TRUE) ;
      vtxtBreak (blkp) ;
      ficheNewGeneUorfStatement (blkp, gmp) ;
    }

  ficheNewGeneFunctionStatement (blkp, gmp) ;
  if (!isCloud) ficheNewGeneProteinStatement (blkp, gmp, TRUE) ;
  if (gmp->tg)
    {
       ficheNewGenePfamPsortStatement (blkp, gmp, isCloud, TRUE) ;
       if (!isCloud) ficheNewGeneNonGoodVariantStatement (blkp, gmp) ;
       ficheNewGeneKozakStatement (blkp, gmp) ;
       ficheNewGenePhosphositeStatement (blkp, gmp) ;
    }
  vtxtBreak (blkp) ;

  ficheNewGenePleaseQuote (blkp, gmp) ;
  return vtxtPtr (blkp) ;
} /* ficheNewGeneAceViewSummary */

/***************************************************************************************/
/* shed + refseq summary + proteome summary */
static int ficheNewGeneSummaryParagraph (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *ptr1, *ptr2, *cp ;
  vTXT buf ;
  AC_TABLE gSummary = ac_tag_table (gmp->gene, "Summary", h) ;
  const char *shedded_from_gene = ac_tag_printable (gmp->gene, "shedded_from", 0) ;
  int nn = 0, ir, refseq ;
  const char *dbSource ;

  if (gmp->Spc == WORM)
    dbSource = "Wormbase" ;
  else if (gmp->Spc == ARA)
    dbSource = "TAIR" ;
  else
    dbSource = "RefSeq" ;

  if (shedded_from_gene)
    {
      nn++ ;
      if (ac_has_tag (gmp->tg, "gt_ag"))
	vtxtPrintf (blkp,
		    "This gene contacts gene %s, but was \"shed\" from it because the two"
		    " genes do not share a single intron boundary. This \"gene\" could be"
		    " considered as an alternative variant of gene %s."
		    , shedded_from_gene
		    , shedded_from_gene
		    ) ;
      else
	vtxtPrintf (blkp,
		    "This gene contacts gene %s, but was \"shed\" from it because it does not"
		    " have standard introns: it may correspond to the alignment of unprocessed"
		    " cDNAs. You could consider this \"gene\" as an alternative variant of gene"
		    " %s, with unspliced support."
		    , shedded_from_gene
		    , shedded_from_gene
		    ) ;
    }


  for (ir = refseq = 0 ; gSummary && ir < gSummary->rows ; ir++)
    {
      ptr1 = ac_table_printable (gSummary, ir, 0, 0) ;
      if (ptr1)
	{
	  if (strstr(ptr1, messprintf ("[%s Summary",  dbSource)))
	    {
	      refseq = 1 ;
	      ptr1 += strlen("[RefSeq Summary") + 1 ;
	      if (!nn++)
		gmpSubSection (blkp, gmp, "tg_RefSeq_Summary", messprintf ("%s summary", dbSource)) ; 
	      else
		vtxtBreak (blkp) ;
	      vtxtPrint (blkp, "[") ;
	      cp = strstr (ptr1, "]") ;
	      if (cp && cp - ptr1 < 3)
		{
		  cp = ac_tag_printable (gmp->gene, "LocusLink", "") ;
		  vtxtPrint (blkp, cp) ;
		}
	    }
	  else 
	    {
	      if (!nn++)
		gmpSubSection (blkp, gmp, "tg_Manual_Summary", "Summary") ;
	    }
	  
	  vtxtPrint (blkp, gtCleanQuotes (ptr1)) ;
	}
    }

  buf = vtxtHandleCreate (h) ;
  if (gmp->markup)
    vtxtMarkup (buf) ;
  ptr2 = ficheNewGeneAceViewSummary (buf, gmp) ;
  if (ptr2 && gmp->tg)
    { 
      const char *tmpTxt ; 
      AC_KEYSET ksMrna = ac_objquery_keyset (gmp->tg, ">mrna ; gt_ag || gc_ag", h) ;
      int nMrna = ksMrna ? ac_keyset_count (ksMrna)  : 0 ;
      AC_KEYSET ksPg = ac_objquery_keyset (gmp->gene, ">Genefinder", h) ;
      int nPg =  ksPg ? ac_keyset_count (ksPg)  : 0 ;
      AC_KEYSET ksGid = ac_objquery_keyset (gmp->gene, ">geneId", h) ;
      int nGid = ksGid ? ac_keyset_count (ksGid) : 0 ;
      AC_KEYSET ksNm = ac_ksquery_keyset (ksGid, ">sequence; Ref_seq", h) ;
      int nNm = ksNm ? ac_keyset_count (ksNm)  : 0 ;
      AC_KEYSET ksGidWithNm = ac_objquery_keyset (gmp->gene
						  , ">Genefinder; >geneId_pg COUNT {>sequence ref_seq} > 0", h) ;
      int nGidWithNm = ksGidWithNm ? ac_keyset_count (ksGidWithNm) : 0 ;
      /*
	AC_KEYSET ksXm = ac_ksquery_keyset (ksGid, ">sequence; IS XM_*", h) ;
	int nXm = ksXm ? ac_keyset_count (ksXm) : 0  ;
      */
      vTXT linkBuf = vtxtHandleCreate (h) ;

      if (gmp->Spc == WORM)
	{
	  tmpTxt = ac_tag_printable (gmp->gene, "WbId", 0) ; 
	  if (tmpTxt) vtxtPrintf (linkBuf, "http://www.wormbase.org/db/get?class=Gene;name=%s", tmpTxt) ;
	}
       else if (nGid == 1)
	{
	  tmpTxt = ac_tag_printable (gmp->gene, "GeneId", 0) ; 
	  if (tmpTxt) vtxtPrintf (linkBuf, ENTREZ_LINK, tmpTxt) ; 
	}
      else if (nGid > 1)
	{
	  AC_TABLE tbl = ac_keyset_table (ksGid, 0, -1, 0, h) ;
	  int ir ;

	  vtxtPrint (linkBuf, "https://www.ncbi.nlm.nih.gov/corehtml/query.fcgi?db=gene&cmd=Retrieve&list_uids=") ;
	  for (ir = 0 ; ir < tbl->rows ; ir++)
	    vtxtPrintf (linkBuf, "%s%s"
			, ir ? "," : ""
			, ac_table_printable (tbl, ir, 0, "")
			) ;
	}
      
      if (gmp && gmp->tg)
	{
	  if (refseq)
	    {
	      refseq = 0 ;
	      vtxtBreak (blkp) ; 
	      vtxtEmptyLine (blkp, 1) ; 
	    }

	  if (nPg)
	    {	  
	      if (gmp->markup)
		vtxtPrint (blkp, "<span class='ace_summary'>") ; 
	      
	      if (nNm)
		{
		  vtxtPrintf (blkp, "%s annotates ", dbSource) ;
		  if (nGid)
		    gmpURL (blkp, gmp,  vtxtPtr (linkBuf)
			    , messprintf ("%s representative transcript%s", isOne (nNm), _multi(nNm))
			    ) ;
		  else
		    vtxtPrintf (blkp, "%s representative transcript%s"
				, isOne (nNm), _multi(nNm)) ;
		  if (nGid > 1 && nGidWithNm > 1)
		    vtxtPrintf (blkp, " from %d genes that we see as a complex", nGid) ;
		  if (nGid > 1 && nGidWithNm  == 1)
		    vtxtPrintf (blkp, " from %d predicted genes that we see see as a single gene", nGid) ;

		  if (gmp->tg && nMrna > 1)
		    {
		      const char *cp1, *cp = messprintf (">in_mrna ; COUNT {>from_gene ; IS %s} > 0"
							 , ac_protect (ac_name(gmp->tg), h)) ;
		      AC_KEYSET ks = ac_ksquery_keyset (ksNm, cp, h) ;
		      AC_TABLE tbl = ks ? ac_keyset_table (ks, 0, -1, 0, h) : 0 ;

		      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
			{
			  char *cp2 ;
			  cp1 = ac_table_printable (tbl, ir, 0, "") ;
			  cp2 = gtMrnaSuffix (ac_name (gmp->gene), cp1, h) ;
			  vtxtPrintf (blkp, "%s%s%s"
				    , ir ? (ir < tbl->rows - 1 ? ", " : " and ") : " (NM included in AceView variant"
				      , ir == 0 && tbl->rows > 1 ? "s " : ""
				      , cp2
				    ) ;
			}
		      if (ir)
			vtxtPrint (blkp, ")") ;
		    }
		}
	      else if (nPg)
		{ 
		  vtxtBreak (blkp) ;
		  vtxtPrintf (blkp, "%s predicts ", dbSource) ;

		  if (nGid && vtxtPtr (linkBuf))
		    gmpURL (blkp, gmp, vtxtPtr (linkBuf)
			    , messprintf ("%s model%s", isOne (nPg), _multi(nPg))
			    ) ;
		  else
		    vtxtPrintf (blkp, "%s model%s"
				, isOne (nPg), _multi(nPg)) ;
		  if (nGid > 1)
		    vtxtPrintf (blkp, " from %d genes", nGid) ;
		}
	  
	      
	      if (nMrna > nPg)
		vtxtPrintf (blkp, ", but %s cDNA sequences in GenBank, dbEST,  Trace and SRA, filtered against clone rearrangements, coaligned on the genome and clustered in a minimal non-redundant way by the manually supervised AceView program, support <a href=\"javascript:openAnchor ('fgene', 'Gene_compact_diagram_1')\"><font color=red> at least %d spliced variant%s</font></a>" 
			    , gmp->spci->speciesName , nMrna, nMrna > 1 ? "s" : "") ;
	      if (gmp->markup)
		vtxtPrint (blkp, "</span>") ;
	      vtxtBreak (blkp) ;
	      vtxtEmptyLine (blkp, 1) ;
	    }
	
	  else /* nPg == 0 */
	    {
	      if (gmp->markup)
		vtxtPrint (blkp, "<span class='ace_summary'>") ;
	      
	      if (strstr(dbSource, "RefSeq"))
		{
		  AC_KEYSET ksPg2 = ac_objquery_keyset (gmp->tg, ">mrna ; > Matching_genefinder", h) ;
		  int nPg2 =  ksPg2 ? ac_keyset_count (ksPg2)  : 0 ;
		  
		  dbSource = nPg2 ? "NCBI RefSeq NM" : "NCBI RefSeq NM or XM" ;
		}
	      if (nMrna > 1) 
		vtxtPrintf (blkp, "There is no %s model for this gene, but %s cDNA sequences in GenBank, coaligned on the genome and clustered in a minimal non-redundant way by the manually supervised AceView program, support at least %d alternative spliced variants" 
			    , dbSource, gmp->spci->speciesName, nMrna) ;
	      else
		vtxtPrintf (blkp, "There is no %s model for this gene, but it is supported by %s cDNA sequences in GenBank", dbSource, gmp->spci->speciesName) ;
	      if (gmp->markup)
		vtxtPrint (blkp, "</span>") ; 
	      vtxtBreak (blkp) ;
	    }
	}
      nn++ ;
      
      gmpSubSection (blkp, gmp, "tg_AceView_Summary", "AceView synopsis, each blue text links to tables and details") ;
      vtxtPrint (blkp, ptr2) ;
    }

  if (0 && nn)
    vtxtHr (blkp, 0, 0) ;

  ac_free (h) ;
  return nn ;
} /* ficheNewGeneSummaryParagraph */

/*****************/

static void ficheNewGeneEcNumberStatement (vTXT blkp, GMP *gmp)
{ 
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = ac_tag_table (gmp->gene, "EC_number", h) ;
  int ir, mx = tbl ? tbl->rows : 0 ; 

  if (mx)
    {
      vtxtBreak (blkp) ;
      vtxtBold (blkp, "EC number: ") ;
      vtxtPrintf (blkp, "This gene encodes protein%s number:", _multi(mx)) ;
      for (ir = 0 ;ir < mx ; ir++)
	vtxtPrintf (blkp, "%s %s", ir ? ",":"", ac_table_printable (tbl, ir, 0, "")) ; 
    }
  ac_free (h) ;
} /* ficheNewGeneEcNumberStatement */

/**************/

static int ficheNewGeneAliasStatement (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nAlias = 0 ;
  int ir, jj ;
  const char *cp, *prefix ;
  char buf[1024] ;
  AC_TABLE tt = 0 ;
  DICT *dict = 0 ;
  vTXT bfr = vtxtHandleCreate (h) ;  
  
  if (gmp->markup)
    {
      vtxtMarkup (bfr) ;
      vtxtPrint (bfr, "<b>Other names:</b>") ;
      vtxtPrint (bfr, " The gene is also known") ;
    }
  else
      vtxtPrint (bfr, " is also known") ;
	
  if (gmp->Spc == HUMAN || 
      gmp->Spc == RAT
      )
    {
      vtxtPrint (bfr, " as ") ;
      dict = gtGeneAliases (gmp, FALSE) ;
      nAlias = dict ? dictMax (dict) - 1 : 0 ;
      if (nAlias)
	{
	  vtxtPrintf (blkp, " %s", vtxtPtr (bfr)) ;
	  for (ir = 1 ; ir <= dictMax (dict)  ; ir++)
	    {
	      cp = dictName (dict, ir) ;
	      if (ir > 1)
		vtxtPrintf (blkp, ir == nAlias ? " or " : ", ") ; 
	      vtxtPrint (blkp, cp) ;
	    }
	}
      tt = ac_tag_table (gmp->gene, "Properties", h) ;
      if (tt && tt->rows)
	{
	  vtxtDot (blkp) ;
	  vtxtPrint (blkp, " It has been described as") ;
	  for (ir = 0 ; ir < tt->rows ; ir++)
	    vtxtPrintf (blkp, "%s %s", ir>0 ? ", " : "" , ac_table_printable (tt, ir, 0, "")) ;
	}      
    }
  else if (gmp->Spc == MOUSE)
    {
      vtxtPrint (bfr, " as ") ;
      dict = gtGeneAliases (gmp, FALSE) ;
      nAlias = dict ? dictMax (dict) - 1 : 0 ;
      if (nAlias)
	{
	  vtxtPrintf (blkp, " %s", vtxtPtr (bfr)) ;
	  for (ir = 1 ; ir <= dictMax (dict)  ; ir++)
	    {
	      cp = dictName (dict, ir) ;
	      if (ir > 1)
		vtxtPrintf (blkp, ir == nAlias ? " or " : ", ") ; 
	      vtxtPrint (blkp, cp) ;
	    }
	}
    }
  else if (gmp->Spc == ARA)
    {
      vtxtPrint (bfr, " as ") ;
      dict = gtGeneAliases (gmp, FALSE) ;
      nAlias = dict ? dictMax (dict) - 1 : 0 ;
      if (nAlias)
	{
	  vtxtPrintf (blkp, " %s", vtxtPtr (bfr)) ;
	  for (ir = 1 ; ir <= dictMax (dict)  ; ir++)
	    {
	      cp = dictName (dict, ir) ;
	      if (ir > 1)
		vtxtPrintf (blkp, ir == nAlias ? " or " : ", ") ; 
	      vtxtPrint (blkp, cp) ;
	    }
	}
      tt = ac_tag_table (gmp->gene, "Properties", h) ;
      if (tt && tt->rows)
	{
	  vtxtDot (blkp) ;
	  vtxtPrint (blkp, " It has been described as") ;
	  for (ir = 0 ; ir < tt->rows ; ir++)
	    vtxtPrintf (blkp, "%s %s", ir>0 ? ", " : "" , ac_table_printable (tt, ir, 0, "")) ;
	}      
    }
  else if (gmp->Spc == WORM)
    {
      dict = dictCreate (64) ;
      jj = 0 ;
      dictAdd (dict, ac_name(gmp->gene), 0) ;
      tt = ac_tag_table (gmp->gene, "Locus", h) ;
      for (ir = 0 ; tt && ir < tt->rows ; ir++)
	{
	  cp = ac_table_printable (tt, ir, 0, 0) ;
	  if (!dictAdd (dict, cp, 0))
	    continue ;
	  if (!jj++)
	    vtxtPrintf (blkp, " as %s", vtxtPtr (bfr)) ;
	  else
	    vtxtPrint (blkp, ",") ;
	  vtxtPrintf (blkp, " %s", cp) ;
	}
      tt = ac_tag_table (gmp->gene, "NewName", h) ;
      for (ir = 0 ; tt && ir < tt->rows ; ir++)
	{
	  cp = ac_table_printable (tt, ir, 0, 0) ;
	  if (!dictAdd (dict, cp, 0))
	    continue ;
	  if (!jj++)
	    vtxtPrint (blkp, vtxtPtr (bfr)) ;
	  else
	    vtxtPrint (blkp, ",") ;
	  vtxtPrintf (blkp, " in Wormgenes/AceView by its positional name %s", cp) ;
	}
      tt = ac_tag_table (gmp->gene, "Genefinder", h) ;
      prefix = "in Wormbase by its cosmid.number name" ;
      for (ir = 0 ; tt && ir < tt->rows ; ir++)
	{
	  char *cq ;
	  
	  cp = ac_table_printable (tt, ir, 0, 0) ;
	  if (!cp || ! strncasecmp (cp, "hmm", 3) )
	    continue ;
	  strncpy (buf, cp, 1023) ;
	  cq = buf + strlen(buf) - 1 ;
	  if (*cq >= 'a' && *cq <= 'z') *cq = 0 ;
	  cp = buf ;
	  if (!dictAdd (dict, cp, 0))
	    continue ;
	  if (!strcasecmp (ac_name(gmp->gene), cp))
	    continue ;
	  if (!jj++)
	    vtxtPrint (blkp, vtxtPtr (bfr)) ;
	  else
	    vtxtPrint (blkp, ",") ;
	  vtxtPrintf (blkp, " %s %s", prefix, cp) ;
	  prefix = "" ;
	}
      prefix = "here by its Gnomon predition.number" ;
      for (ir = 0 ; gmp->markup && tt && ir < tt->rows ; ir++)
	{
	  char *cq ;
	  
	  cp = ac_table_printable (tt, ir, 0, 0) ;
	  if (!cp || strncasecmp (cp, "hmm", 3) )
	    continue ;
	  strncpy (buf, cp, 1023) ;
	  cq = buf + strlen(buf) - 1 ;
	  if (*cq >= 'a' && *cq <= 'z') *cq = 0 ;
	  cp = buf ;
	  if (!dictAdd (dict, cp, 0))
	    continue ;
	  if (!strcasecmp (ac_name(gmp->gene), cp))
	    continue ;
	  if (!jj++)
	    vtxtPrint (blkp, vtxtPtr (bfr)) ;
	  else
	    vtxtPrint (blkp, ",") ;
	  vtxtPrintf (blkp, " %s %s", prefix, cp) ;
	  prefix = "" ;
	}
      cp = gtCloneGroup (gmp->gene, 0) ;
      if (cp && dictAdd (dict, cp, 0))
	{
	  if (!jj++)
	    vtxtPrint (blkp, vtxtPtr (bfr)) ;
	  else
	    vtxtPrint (blkp, ",") ;
	  vtxtPrintf (blkp, " in NextDB, the Nematode expression pattern database, as CE%s", cp) ;
	}
      nAlias = jj ;
    }

  if (gmp->markup)
    vtxtPrint (blkp, "</span>") ;

  dictDestroy (dict) ;

  ac_free (h) ;
  return nAlias ;
} /* ficheNewGeneAliasStatement */

/**************/

static int ficheNewDbXrefParagraphContent (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  int ir, nn = 0, nFunctional = 0 ;
  const char *ccp ;
  char *cp ;
  char linkBuf[vONLINEMAXURLSIZE] ;
  AC_ITER iter = 0 ;
  AC_OBJ clo = 0 ;
  AC_TABLE oTmp ;
  vTXT bfr = vtxtHandleCreate (h) ;

  if (gmp->markup) vtxtMarkup (bfr) ;

  if (gmp->Spc == WORM && (ccp = ac_tag_printable (gmp->gene, "Genefinder", 0)))
    {
      sprintf (linkBuf, "http://www.wormbase.org/db/get?class=Gene;name=%s",ccp);
      if (nn++) vtxtPrintf (bfr, ", ") ;
      gmpURL (bfr, gmp, linkBuf, "WormBase")  ;
    }

  if (gmp->Spc == WORM && gmp->tg && (ccp = ac_tag_printable (gmp->tg, "Clone_group", 0)))
    {
      sprintf (linkBuf, "http://nematode.lab.nig.ac.jp/db2/ShowGeneInfo.php?celk=CELK%05d",atoi(ccp+2));
      if (nn++) vtxtPrintf (bfr, ", ") ;
      gmpURL (bfr, gmp, linkBuf, "NextDB")  ;
    }

  if (gmp->Spc == WORM &&
      (1 || (ccp = gtPiano (gmp->gene))) )
    {
      sprintf (linkBuf
	       , "http://nematoda.bio.nyu.edu:8001/cgi-bin/browse/card.orf.cgi?Detailed_View=1&query=%s"
	       , ac_tag_printable (gmp->gene, "Genefinder", "*")) ;
      cp = linkBuf + strlen (linkBuf) - 1 ;
      if (cp > linkBuf && *cp >= 'a' && *cp <= 'z')
	*cp = 0 ;
      if (nn++) vtxtPrint (bfr, ", ") ;
      gmpURL (bfr, gmp, linkBuf, "RNAiDB")  ;
    }

  if (0 && /* 2007_03_06 the vidal db is totally oboslete, anyway our link via gtOst is also obsolete */
      gmp->Spc == WORM &&
      (ccp = gtOst (gmp->gene)))
    {
        sprintf (linkBuf, "http://worfdb.dfci.harvard.edu/searchallwormorfs.pl?by=name&sid=%s",ccp);
      if (nn++) vtxtPrint (bfr, ", ") ;
      gmpURL (bfr, gmp, linkBuf, "Vidal")  ;
    }

  if (0 && /*  2007_03_06 obsolete */
      gmp->Spc == WORM &&
      (ccp = gtHyman (gmp->gene)))
    {
      sprintf (linkBuf, "http://worm-srv1.mpi-cbg.de/dbScreen/") ;
      if (nn++) vtxtPrint (bfr, ", ") ;
      gmpURL (bfr, gmp, linkBuf, "Hyman movies")  ;
    }

  nFunctional = nn = 0 ;
  if (gmp->Spc == HUMAN)
    { /************************** OMIM  *****************************/
      if ((oTmp = ac_tag_table (gmp->gene, "Extern", h)))
	for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
	  {
	    const char *ccp = ac_table_printable (oTmp, ir, 0, 0) ;
	    
	    
	    if (!ccp || strncmp (ccp, "OMIM_", 5))
	      continue ;
	    { 
	      AC_OBJ myOmim = ac_table_obj (oTmp, ir, 0, h) ;
	      if ( ac_has_tag (myOmim , "OMIM_molecular"))
		continue ;
	    }
	    if (nn++)
	      vtxtPrint (bfr, ", ") ;
	    if (! nFunctional++)
	      vtxtPrint (bfr, " manual annotations from ") ;
	    sprintf (linkBuf, OMIM_LINK, ccp + 5) ;
	    gmpURL (bfr, gmp, linkBuf, ccp) ;
	  }
    }

   if (gmp->Spc == MOUSE)
    { /************************** MGI  *****************************/
      if ((oTmp = ac_tag_table (gmp->gene, "Extern", h)))
	for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
	  {
	    const char *ccp = ac_table_printable (oTmp, ir, 0, 0) ;
	    
	    if (!ccp || strncmp (ccp, "MGI_", 4))
	      continue ;

	    if (nn++)
	      vtxtPrint (bfr, ", ") ;
	    if (! nFunctional++)
	      vtxtPrint (bfr, " manual annotations from ") ;
	    sprintf (linkBuf, MGI_LINK, ccp + 4) ;
	    gmpURL (bfr, gmp, linkBuf, ccp) ;
	  }
    }

   if (0 && gmp->Spc == RAT)
    { /************************** MGI  *****************************/
      if ((oTmp = ac_tag_table (gmp->gene, "Extern", h)))
	for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
	  {
	    const char *ccp = ac_table_printable (oTmp, ir, 0, 0) ;
	    
	    if (!ccp || strncmp (ccp, "MGI_", 4))
	      continue ;

	    if (nn++)
	      vtxtPrint (bfr, ", ") ;
	    if (! nFunctional++)
	      vtxtPrint (bfr, " manual annotations from ") ;
	    sprintf (linkBuf, MGI_LINK, ccp + 4) ;
	    gmpURL (bfr, gmp, linkBuf, ccp) ;
	  }
    }

   if (gmp->Spc == ARA)
    { /************************** TAIR  *****************************/
      if ((oTmp = ac_tag_table (gmp->gene, "Extern", h)))
	for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
	  {
	    const char *ccp = ac_table_printable (oTmp, ir, 0, 0) ;
	    
	    if (!ccp || strncmp (ccp, "TAIR:", 4))
	      continue ;

	    if (nn++)
	      vtxtPrint (bfr, ", ") ;
	    if (! nFunctional++)
	      vtxtPrint (bfr, " manual annotations from ") ;
	    sprintf (linkBuf, TAIR_LINK, ccp + 5) ;
	    gmpURL (bfr, gmp, linkBuf, ccp) ;
	  } 
    }

 if (gmp->Spc == HUMAN)
    { /************************** GAD *****************************/
       AC_ITER iter = ac_objquery_iter (gmp->gene, ">Extern  GAD AND NOT AntiGad", h) ;
       AC_OBJ gad = 0 ;

       if (iter)
	 while (ac_free (gad), (gad = ac_iter_obj (iter))) 
	   {    
	     if (! nFunctional++)
	       vtxtPrint (bfr, ", manual annotations from") ;
	     /* GAD as vanished since 2018 
		const char* mygeneid = ac_tag_printable (gmp->gene, "geneid", 0) ;
	       if (nn++)
	       vtxtPrint (bfr, ", ") ;
	       sprintf (linkBuf,  "http://geneticassociationdb.nih.gov/cgi-bin/tableview.cgi?table=allview&cond=LOCUSNUM=%s"
	       , mygeneid) ;
	       gmpURL (bfr, gmp, linkBuf, messprintf (" GAD")) ;
	     */
	    break ; /* we want only one GAD */
	   }
       ac_free (gad) ; 
    }
  
  if (gmp->Spc == HUMAN)
    { /************************** HPRD *****************************/
       char buf1[256] ;

      if ((oTmp = ac_tag_table (gmp->gene, "HPRD", h)))
	for (ir=0; ir < oTmp->rows && oTmp->cols >= 1;ir++)
	  {
	    const char *cp = ac_table_printable (oTmp, ir, 0, 0) ;
	    
	    if (!cp)
	      continue ;
	    sprintf (buf1," HPRD") ;

	    if (nn++)
	      vtxtPrint (bfr, ",") ;
	    if (! nFunctional++)
	      vtxtPrint (bfr, " manual annotations from") ;
	    sprintf (linkBuf, "http://www.hprd.org/summary?protein=%s&isoform_id=%s_1"
		     , cp, cp) ;    
	    gmpURL (bfr, gmp, linkBuf, buf1) ;
	  }
     }

  if (1)
    { /************************** KEGG *****************************/
       AC_ITER iter = ac_objquery_iter (gmp->gene, ">Extern KEGG", h) ;
       AC_OBJ kegg = 0 ;
       const char* mygeneid = ac_tag_printable (gmp->gene, "geneid", 0) ;

      if (iter)
	{
	  while (ac_free (kegg), (kegg = ac_iter_obj (iter))) 
	    {
	      if (nn++)
		vtxtPrint (bfr, ",") ;
	      if (! nFunctional++)
		vtxtPrint (bfr, " manual annotations from") ;
	      
	      sprintf (linkBuf, "http://www.genome.jp/dbget-bin/show_pathway?hsa%s%s%s"
		       , ac_name(kegg)+5
		       , mygeneid ? "+" : ""
		       , mygeneid ? mygeneid : ""
		       ) ;
	      gmpURL (bfr, gmp, linkBuf, messprintf (" %s", ac_name(kegg))) ;
	    }
	  ac_free (kegg) ;
	}
    }

 
  if (1)
    { /************************** Phosphosite *****************************/
       AC_ITER iter = ac_objquery_iter (gmp->gene, ">Extern Phosphosite", h) ;
       AC_OBJ obj = 0 ;
       const char* uniprot ;

      if (iter)
	while (ac_free (obj), (obj = ac_iter_obj (iter))) 
	  {	
	    uniprot = ac_tag_printable (obj, "Uniprot", 0) ;
	    if (! uniprot)
	      continue ;
	    if (nn++)
	      vtxtPrint (bfr, ",") ;
	    if (! nFunctional++)
	      vtxtPrint (bfr, " manual annotations from") ;

	    sprintf (linkBuf, "http://www.phosphosite.org/uniprotAccAction.do?id=%s"
		     , uniprot
		     ) ;
	    gmpURL (bfr, gmp, linkBuf, " PhosphoSite") ;
	    break ; /* we want at most 1 link */
	  }
      ac_free (obj) ;
     }

 
  /******************************** SNP **************************************/
  nFunctional = 0;
  if (gmp->Spc != WORM && 
      (oTmp = ac_tag_table (gmp->gene, "GeneId", h)))
    for (ir=0; ir < oTmp->rows && oTmp->cols >= 1;ir++)
      {

	ccp = ac_table_printable (oTmp, ir, 0, 0) ;
	if (nn++)
	  vtxtPrint (bfr, ",") ;
	if (! nFunctional++)
	  vtxtPrint (bfr, " the") ;
	sprintf (linkBuf, "https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?chooseRs=all&locusId=%s", ccp) ;
	gmpURL (bfr, gmp, linkBuf, " SNP")  ;

	vtxtPrint (bfr, " view") ;
      }

  /******************************** ENTREZ **************************************/
  nFunctional = 0;
  if (gmp->Spc != WORM && 
      (oTmp = ac_tag_table (gmp->gene, "GeneId", h)))
    {
      for (ir=0; ir < oTmp->rows && oTmp->cols >= 1;ir++)
	{
	  
	  ccp = ac_table_printable (oTmp, ir, 0, 0) ;
	  if (nn++)
	    vtxtPrint (bfr, ",") ;
	  if (! nFunctional++)
	    vtxtPrint (bfr, " gene overviews from") ;

	  sprintf (linkBuf, "https://www.ncbi.nlm.nih.gov/gene/%s", ccp) ;
	  if (0 && oTmp->rows == 1) /* always show the geneid */
	    gmpURL (bfr, gmp, linkBuf, " Gene")  ;
	  else 
	    gmpURL (bfr, gmp, linkBuf, messprintf ("%s %s", ir > 0 ? "" : " Gene", ccp))  ;
	}
      
      /************************** GeneCard direct *****************************/
      if (gmp->Spc == HUMAN)
	{ 
	  /* before 2020: sprintf (linkBuf, "http://bioinformatics.weizmann.ac.il/cards-bin/carddisp?gc_id=%s", ccp); */
	  nn++ ;
	  if ((ccp = ac_tag_printable (gmp->gene, "GeneCard_id", 0)))
	    sprintf (linkBuf, "https://www.genecardsorg/cgi-bin/carddisp.pl?gc_id=%s", ccp);
	  else if ((ccp = ac_tag_printable (gmp->gene, "LocusLink", 0)))
	    sprintf (linkBuf, "https://www.genecardsorg/cgi-bin/carddisp.pl?gene=%s", ccp) ;
	  if (ccp)
	    {
	      ccp = ac_table_printable (oTmp, ir, 0, 0) ;
	      if (nn++)
		vtxtPrint (bfr, ",") ;
	      if (! nFunctional++)
		vtxtPrint (bfr, " gene overviews from") ;

	      gmpURL (bfr, gmp, linkBuf, " GeneCards")  ;
	    }
	}
    }
  /******************************** HINV **************************************/
  if (gmp->Spc == HUMAN)
    {
      int nhinv = 0 ;
      iter = ac_objquery_iter (gmp->tg, ">read HINV_cluster_id", 0) ;
      while (ac_free(clo), nhinv==0 && iter && (clo = ac_next_obj (iter)))
	{
	  ccp = ac_tag_printable (clo, "HINV_cluster_id", 0) ;
	  if (!ccp)
	    continue ;
	    
	  if (nn++)
	    vtxtPrint (bfr, ",") ;
	  if (! nFunctional++)
	    vtxtPrint (bfr, " gene overviews from") ;
	  nhinv++ ;
	  /* was */
	  if (0) sprintf (linkBuf, "http://www.jbirc.aist.go.jp/hinv/spsoup/locus_view?hix_id=%s&status=full",ccp) ;
	  /* changed 2008_03_10 following email from hinv */
          if (1) sprintf (linkBuf, "http://h-invitational.jp/hinv/spsoup/locus_view?hix_id=%s&status=full",ccp) ;
    
	  gmpURL (bfr, gmp, linkBuf, " H-Inv")  ;
	  if (0)
	    vtxtPrint (bfr, "<img src='images/btn_locusview.gif' width='85' height='16' border='0' alt='H-INV' />") ;
	  if (0)
	    vtxtPrint (bfr, "<img src='images/btn_locusview.gif' width='85' height='16' border='0'>") ;
	}
      ac_free (clo) ;
      ac_free (iter) ;
    }

  /************************** EC gene gone as of 2018 *****************************/
  nFunctional = 0 ;
  if (0 &&   
      gmp->Spc == HUMAN && gmp->tg && 
      (ccp = ac_tag_printable (gmp->gene, "LocusLink", 0)))
    {
      AC_KEYSET ks = ac_objquery_keyset (gmp->tg, ">mrna; best_in_gene ; >CDS_covered_by ; {ref_mrna} SETELSE {dna}", h) ;
      AC_TABLE tbl = ks ? ac_keyset_table (ks, 0, 1, 0, h) : 0 ;
      if (tbl && tbl->rows)
	{
	  if (nn++)
	    vtxtPrint (bfr, ",") ;
	  sprintf (linkBuf, "http://genome.ewha.ac.kr/cgi-bin/ECexpress.py?organism=human&confidence=C&query=%s"
		   , ac_table_printable (tbl, 0, 0, "")) ;
	  if (! nFunctional++)
	    vtxtPrint (bfr, " expression data from") ;
	  gmpURL (bfr, gmp, linkBuf, " ECgene")  ;
	}
    }

  /************************** GENE (gobbled UniGene in 2018)  *****************************/
  if (1 &&
      gmp->Spc != WORM &&
    (oTmp = ac_tag_table (gmp->gene, "GeneId", h)))
    {
      for (ir=0; ir < oTmp->rows && oTmp->cols >= 1;ir++)
	{
 	  
	  ccp = ac_table_printable (oTmp, ir, 0, 0) ;
	  if (nn++)
	    vtxtPrint (bfr, ",") ;
	  if (! nFunctional++)
	    vtxtPrint (bfr, " expression data from") ;

	  sprintf (linkBuf, "https://www.ncbi.nlm.nih.gov/gene/%s", ccp) ;
	  if (0 && oTmp->rows == 1) /* always show the geneid */
	    gmpURL (bfr, gmp, linkBuf, " Gene")  ;
	  else 
	    gmpURL (bfr, gmp, linkBuf, messprintf ("%s %s", ir > 0 ? "" : " Gene", ccp))  ;
	}
    }

  /************************** Unigene gone as of 2018  *****************************/
  if (0 &&
      gmp->Spc != WORM && 
      (oTmp = ac_tag_table (gmp->gene, "UniGene", h)))
    for (ir=0; ir < oTmp->rows && oTmp->cols >= 1;ir++)
      {
	ccp = ac_table_printable (oTmp, ir, 0, 0) ;
	
	if (nn++)
	  vtxtPrint (bfr, ",") ;
	if (! nFunctional++)
	  vtxtPrint (bfr, " expression data from") ;
	if (gmp->Spc == HUMAN)
	  sprintf (linkBuf, "https://www.ncbi.nlm.nih.gov/UniGene/ESTProfileViewer.cgi?uglist=Hs.%s",ccp + 3);
	else if (gmp->Spc == MOUSE)
	  sprintf (linkBuf, "https://www.ncbi.nlm.nih.gov/UniGene/ESTProfileViewer.cgi?uglist=Mm.%s",ccp + 3);
	else if (gmp->Spc == RAT)
	  sprintf (linkBuf, "https://www.ncbi.nlm.nih.gov/UniGene/ESTProfileViewer.cgi?uglist=Rn.%s",ccp + 3);
	if (oTmp->rows == 1)
	  gmpURL (bfr, gmp, linkBuf, " UniGene")  ;
	else 
	  gmpURL (bfr, gmp, linkBuf, messprintf ("%s %s", ir > 0 ? "" : "  UniGene", ccp))  ;
      }

  nFunctional = 0 ;
  if (gmp->Spc == HUMAN || gmp->Spc == MOUSE || gmp->Spc == RAT)
    {/************************** UCSC *****************************/
      AC_TABLE tt = 0 ;
      if ( (tt = ac_tag_table (gmp->gene, "IntMap", 0)) && tt->cols >= 3)
        {
	  const char *ucsctarget = "hg18" ;
	  const char *nnam = ac_table_printable (tt, 0, 0, "") ;
	  if (gmp->Spc == HUMAN)
	    ucsctarget = "hg19" ;
	  else if (gmp->Spc == MOUSE)
	    ucsctarget = "mm9" ;
	  else if (gmp->Spc == RAT)
	    ucsctarget = "rn18" ;
	  if (strstr (nnam, "_"))
	    nnam = strstr (nnam, "_")+1 ;
	  sprintf (linkBuf,
		   "http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&acembly=pack&hgt.out1=1.5x&est=dense&stsMap=hide&knownGene=squish&RefSeq=squish&mgcGenes=hide&intronEst=dense&snp126=hide&rmsk=hide&position=chr%s:%d-%d"
		   , ucsctarget
		   , nnam, ac_table_int (tt, 0, 1, 0), ac_table_int (tt, 0, 2, 0)) ;

	  if (nn++)
	    vtxtPrint (bfr, ",") ;
	  if (! nFunctional++)
	    vtxtPrint (bfr, " molecular and other annotations from") ;
	  gmpURL (bfr, gmp, linkBuf, " UCSC")  ;
        }
      ac_free (tt) ;
    }
  
  /************************** GOLD *****************************/
  if (1 && /* aug 2005, simplify the link to GOLD, obsolete */
      gmp->tg && gmp->Spc == HUMAN)
    {
      AC_KEYSET hinv = ac_objquery_keyset (gmp->tg, ">cdna_clone; HINV_libs", 0) ;

      if (ac_keyset_count (hinv) > 0)
	{
	  nn++ ;
	  sprintf (linkBuf
		   , "https://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/GOLD") ;
	  if (0)
	    {
	      vtxtBreak (bfr) ;
	      vtxtPrint (bfr,
			 "Some clones from this gene come from one of the four "
			 "large scale cDNA projects (KIAA, DKFZ, FLJ or MGC clones) and have been "
			 "subjected to the GOLD analysis"
			 ". Any feedback or comments on the "
			 ) ;
	      gmpURL (bfr, gmp, linkBuf, "draft paper") ;
	      vtxtPrint (bfr, " are welcome") ;
	    }
	  else /* august 30 2005 */
	    {
	      if (nFunctional++) vtxtPrint (bfr, ", or") ;
	      vtxtPrint (bfr, " our ") ;
	      gmpURL (bfr, gmp, linkBuf, "GOLD") ;
	      vtxtPrint (bfr, " analysis") ;
	    }
	  vtxtBreak (bfr) ;
	}
      ac_free (hinv) ;
    }

  if (gmp->Spc == HUMAN &&
      ac_has_tag (gmp->gene, "Balise") &&
      (ccp = ac_tag_printable (gmp->gene, "GeneId", 0))
      )
    {
      nn++ ;
      sprintf (linkBuf, "https://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/av.cgi?db=36a&c=geneid&q=%s", ccp) ;
      vtxtDot (bfr) ;
      vtxtPrint (bfr, "The previous AceView annotation is ") ;
      gmpURL (bfr, gmp, linkBuf, "here")  ;
    }
  
  if ((ccp = vtxtPtr (bfr)))
    {
      vtxtBreak (blkp) ;
      vtxtBold (blkp, "<b>Links to:</b> ") ;
      vtxtPrint (blkp, ccp) ;
    }
  
  ac_free (h) ;
  return nn ;
} /* ficheNewDbXrefParagraphContent */

/***************************************************************************************/
/*
  prefix= 0->void, 1->It, 2->this gene/mrna ..
  type = 1->gene, 2->tg, 3->mrna
*/

static void ficheNewGeneCoversValue (vTXT blkp, GMP *gmp, int type)
{
  int a1, a2, dx, n1 = 0, n2 = 0 ;
  AC_TABLE oTmp = 0 ;
  AC_OBJ oCosmid = 0 ;
  char buf1[1024], buf2[256] ;
  AC_HANDLE h = ac_new_handle () ;

  switch (type)
    {
    case 1:
      return ;
    case 2: 
      n1 = 2 ; n2 = 3 ;
      oTmp = ac_tag_table (gmp->tg, "Covers", h) ; 
      oCosmid = ac_tag_obj (gmp->tg, "Genomic_sequence", h) ;
      break ;
    case 3:  
      n1 = 2 ; n2 = 4 ;
      oTmp = ac_tag_table (gmp->mrna, "Covers", h) ; 
      break ;
    case 4: /* genefinder */
      n1 = 1 ; n2 = 2 ;
      oTmp = ac_tag_table (gmp->pg, "IntMap", h) ; 
      break ;
    }

  if (oTmp)
    {     
      /* the model is different for tg and mRNA ! */
      a1 = ac_table_int (oTmp, 0, n1, 0) ;
      a2 = ac_table_int (oTmp, 0, n2, 0) ; 

      dx = a2 - a1 ;
      if (dx < 0) dx = - dx ;
      dx++ ; /* plato */
      sprintf (buf2, "%.2f kb", 0.001 * dx) ;

      if (type == 3)
	{
	  sprintf (buf1, "DNA_mRNA:99999:99999:0") ;
	  gmpFakeObjLink (blkp, gmp, buf1, gmp->mrna, buf2) ;
	}
      else if (type == 2 && oCosmid)
	{
	  sprintf (buf1, "DNA:%d:%d:0", a1, a2) ;
	  gmpFakeObjLink (blkp, gmp, buf1, oCosmid, buf2) ;
	}
      else
	vtxtPrint (blkp, buf2) ;	    
    }
  ac_free (h) ;
} /* ficheNewGeneCoversValue */


static void ficheNewGenePositionOnChromosomeStatement (vTXT blkp, GMP *gmp, int type, int prefix)
{
  int a1, a2, dx ;
  AC_TABLE oTmp = 0 ;
  AC_HANDLE h = ac_new_handle () ;

  switch (type)
    {
    case 1:
      oTmp = ac_tag_table (gmp->gene, "IntMap", h) ; 
      break ;
    case 2: 
      oTmp = ac_tag_table (gmp->tg, "IntMap", h) ; 
      break ;
    case 3:  
      oTmp = ac_tag_table (gmp->tg, "IntMap", h) ; 
      break ;
    }

  if (oTmp)
    {     
      /* the model is different for tg and mRNA ! */
      a1 = ac_table_int (oTmp, 0, 1, 0) ;
      a2 = ac_table_int (oTmp, 0, 2, 0) ; 

      dx = a2 - a1 ;
      if (dx < 0) dx = - dx ;
      dx++ ; /* plato */

      if (prefix) 
	vtxtDot (blkp) ;
      vtxtPrintf (blkp, "In AceView, it covers ") ;
      ficheNewGeneCoversValue (blkp, gmp, type) ;
      vtxtPrintf (blkp
		  , ", from %d to %d (%s), on the %s strand"
		  , a1, a2
		  , dbNam
		  , a1 > a2 ? "reverse" : "direct"
		  ) ;
      
      if (!prefix)    
	vtxtPrintf (blkp, " of chromosome %s", ac_table_printable (oTmp, 0, 0, "")) ;
    }
  ac_free (h) ;
} /* ficheNewGenePositionOnChromosomeStatement */

/**********************/
/*
  This gene || unc-32
  ,also known as ZK637.8,
  maps on chromosome II at position +3.2cM (measured), +3.4 cM (interpolated).
  It covers 6210 bp, from 3002102 to 3020202, on the reverse strand.

  prefix= 0->void, 1->It, 2->this gene/mrna ..
  type = 1->gene, 2->tg, 3->mrna
*/
static int ficheNewGeneMappingSentence (vTXT blkp, GMP *gmp, int prefix)
{
  int type = 2 ; /* the only way it is now used is for a tg */
  AC_TABLE gMap ;
  AC_TABLE oTmp ;
  const char *ptr ;
  char buf1[512], buf2[512] ;
  int nn = 0 ;
  float z1 = -9999, z2 ;
  AC_HANDLE h = ac_new_handle () ;
  
  buf1[0] = buf2[0] = 0 ;
  if (prefix == 2) vtxtPrintf (blkp, "This %s %s", type == 3 ? "RNA" : "gene", ac_name (type == 3 ? gmp->mrna : gmp->gene)) ; 
  else if (prefix == 1) vtxtPrintf (blkp, "It") ;
  else ; /* no prefix needed */
  if ( gmp->Spc==WORM && 
       (gMap = ac_tag_table (gmp->gene, "Map", h)) && 
       ( z1 = ac_table_float (gMap, 0, 2, 9999)) != 9999)
    {
      vtxtPrintf (blkp, " maps on chomosome %s", ac_table_printable (gMap, 0, 0, "")) ;
      vtxtPrintf (blkp, " at position %s%4.2f (measured by recombination)" 
		  , z1 >= 0 ? "+" : "",
		  z1) ;
      nn++ ;
    }   
  if ( (gMap = ac_tag_table (gmp->gene, "InterpolatedMap", h)) && 
       (z2 = ac_table_float (gMap, 0, 1, 9999)) != 9999) 
    {
      if (!nn)
	vtxtPrintf (blkp, " maps on chomosome %s at position", ac_table_printable (gMap, 0, 0, "")) ;
      else
	vtxtPrintf (blkp, ", ") ;
      vtxtPrintf (blkp, " %s%4.2f (interpolated)" 
		  , z2 >= 0 ? "+" : "",
		  z2) ;
      nn++ ;
    }   
  

  if (!nn &&
      (ptr = ficheChromName (ac_tag_printable (gmp->gene, "Map", 0)))
      )
    { 
      nn++ ;
      vtxtPrintf (blkp, " maps on chromosome %s", ptr) ; nn++ ; 
      strncpy (buf1, ptr, 511) ;
    }
  
  
  if (!nn &&
      (ptr = ficheChromName (ac_tag_printable (gmp->gene, "IntMap", 0))) &&
      ptr && 
      strlen(ptr) < 4)
    {
      nn++ ; vtxtPrintf (blkp, " maps on chromosome %s", ptr) ; nn++ ; 
      strncpy (buf2, ptr, 511) ;
    }
  
  if ((oTmp = ac_tag_table (gmp->gene, "Cytogenetic", h)) &&
      (ptr = strnew (ac_table_printable (oTmp, 0, 0, 0), h)))
    {
      const char *ptr2 = ac_table_printable (oTmp, 0, 1, "Entrez Gene") ;
      nn++ ; vtxtPrintf (blkp, ", at %s according to %s", ptr, ptr2) ; nn++ ;
    }
  
  
  if (*buf1 && *buf2 && strcasecmp(buf1, buf2))
    {
      if (ac_has_tag (gmp->tg,"intron_boundaries"))
	vtxtPrintf (blkp,
		    ". The gene shown here is probably not the gene itself "
		    "but a close paralog since it maps on chromosome %s"
		    , buf2) ;
      else
	vtxtPrintf (blkp, ". The gene shown here is probably a pseudogene") ;
    }
  
  ficheNewGenePositionOnChromosomeStatement (blkp, gmp, type, nn) ; 

  ac_free (h) ;
  return nn ; 
} /* ficheNewGeneMappingSentence */

/***************************/

int ficheNewAceKogStatement (vTXT blkp, GMP *gmp, BOOL doTitle)
{
  AC_HANDLE h = 0, h0 = ac_new_handle() ;
  int ii, ir1, ir2, jr, kr, nn = -1 ;
  char	linkBuf[vONLINEMAXURLSIZE] ;
  const char *ccp, *ccq ;
  char *cp, *cq, *cr ;
  float score, bestScore = 0 ;
  AC_TABLE kantors = 0, akg = 0, akg2 = 0 ;
  AC_OBJ Kantor = 0 ;
  Array titles = 0 ;
  Array scores = 0 ;
  Array expects = 0 ;
  DICT *dict = 0 ;
  AC_KEYSET ks = 0 ;
  const char *targetSpecies=0, *targetAceView=0, *targetDb=0 ;
  BOOL isInteresting = FALSE ;

  if (gmp->view == 'm' && gmp->product)
    ks = ac_objquery_keyset (gmp->product, ">Kantor ; AKG", h0) ;
  else if (gmp->view != 'm' && gmp->gene)
    ks = ac_objquery_keyset (gmp->gene, ">product best_product || very_good_product ; >Kantor ; AKG", h0) ;

  if (ks) 
    kantors = ac_keyset_table (ks, 0, -1, 0, h0) ;
  for (ir2 = 0 ; ir2 < 4 ; ir2++)
    { 
      bestScore = 0 ; 
      akg = akg2 = 0 ;
      ac_free (h) ;
      h = ac_new_handle() ;
      for (ir1 = 0 ; kantors && ir1 < kantors->rows ; ir1++)
	{
	  Kantor = ac_table_obj (kantors, ir1, 0, h) ;
	  
	  switch (ir2)
	    {
	    case 0: 
	      targetSpecies = "human" ; targetAceView = "AceView gene" ; targetDb = "human" ; 
	      akg = ac_tag_table (Kantor, "AceKogHuman", h) ;
	      break ;
	    case 1: 
	      targetSpecies = "mouse" ; targetAceView = "AceView gene" ; targetDb = "mouse" ; 
	      akg = ac_tag_table (Kantor, "AceKogMouse", h) ;
	      break ;
	    case 2: 
	      targetSpecies = "C.elegans" ; targetAceView = "AceView/WormGene" ; targetDb = "worm" ; 
	      akg = ac_tag_table (Kantor, "AceKogWorm", h) ;
	      break ;
	    case 3: 
	      targetSpecies = "A.thaliana" ; targetAceView = "AceView gene" ; targetDb = "ara" ; 
	      akg = ac_tag_table (Kantor, "AceKogAra", h) ;
	      break ;
	    }
	  for (jr = 0 ; akg && jr < akg->rows ; jr++)
	    {
	      if ((ccp = ac_table_printable (akg, jr, 0, 0))) /* gene name  */
		{
		  if (!dict) dict = dictHandleCreate (4, h) ;
		  if (dictFind (dict, ccp, 0))
		    continue ;
		  dictAdd (dict, ccp, &nn) ;
		  if (ir2 == 2 && strstr (ccp, "-"))
		    isInteresting = TRUE ;
		  if (! akg2)
		    akg2 = ac_tag_table (Kantor, "AceKog", h) ;
		  if (! akg2)
		    continue ;
		  if (!titles) titles = arrayHandleCreate (16, char *, h) ;
		  if (!expects) expects = arrayHandleCreate (16, char *, h) ;
		  if (!scores) scores = arrayHandleCreate (16, float, h) ;
		  for (kr = 0 ; akg2 && kr < akg2->rows ; kr++)
		    {
		      ccq = ac_table_printable (akg2, kr, 0, "") ;
		      if (strstr (ccq, ccp))
			{
			  ccq = ac_table_printable (akg2, kr, 7, "") ;
			  cp = array (titles, nn, char*) = strnew (ccq, h) ;
			  cq = cr = strstr (cp, "[") ;
			  if (cq) 
			    {
			      *cq-- = 0 ;
			      while (cq > cp &&
				     (
				      *cq == ' ' || *cq == '.' || *cq == ','
				      )
				     )
				*cq-- = 0 ; 
			      cp = strstr (cr+1, "eValue = ") ;
			      cq = cp + strlen (cp) - 1 ;
			      if (*cq == ',') *cq = 0 ;
			      if (cp)
				array (expects, nn, char *) = strnew (cp + 9, h) ;
			    }
			  score = ac_table_float (akg2, kr, 2, 0) ;
			  array (scores, nn , float) = score ; 
			  if (bestScore < score) bestScore = score ;
			  break ;
			}
		    }
		}
	    }
	}
    
      nn = 0 ;	  
      if (dict && scores)
	{ 
	  for (ii = 1 ; ii <= dictMax (dict) ; ii++)
	    {
	      if (array (scores, ii, float) >= .88 * bestScore)
		nn++ ;
	    }
	}
      
      if (nn && dict && scores)
	{ 
	  vtxtBreak (blkp) ;
	  vtxtBold (blkp, messprintf ("The closest %s gene%s", targetSpecies, nn > 1 ? "s" : "")) ;
	  vtxtPrint (blkp, ", according to BlastP,") ;
	  if (nn == 1)
	    vtxtPrintf (blkp, " is the %s ", targetAceView) ;	  
	  else
	    vtxtPrintf (blkp, " are the %ss ", targetAceView) ;	  
	  for (ii = 1, nn = 0 ; ii <= dictMax (dict) ; ii++)
	    {
	      if (array (scores, ii, float) < .88 * bestScore)
		continue ;
	      if (nn++) vtxtPrint (blkp, ", ") ;
	      ccp = dictName (dict, ii) ;	
	      sprintf (linkBuf, "https://www.ncbi.nlm.nih.gov/IEB/Research/Acembly/av.cgi?db=%s&q=%s",targetDb, ccp);
	      gmpURL (blkp, gmp, linkBuf, ccp) ;
	      cp = titles ? array (titles, ii, char *) : 0 ;
	      if (cp)
		{
		  if (0) /* when i will have the titles truly avalilable */
		    {
		      vtxtPrint (blkp, " (") ;
		      vtxtPrintf (blkp, "%s", cp) ;
		      cp = expects ? array (expects, ii, char *) : 0 ;
		      if (cp)
			vtxtPrintf (blkp, ", e=%s", cp) ;
		      vtxtPrint (blkp, ")") ;
		    }
		  else
		    {
		      char *cpe ;
		      
		      cp = expects ? array (expects, ii, char *) : 0 ;
		      cpe = cp ? strstr(cp, ",") : 0 ;
		      if (cpe) *cpe = 0 ;
		      cpe = cp ? strstr(cp, "e") : 0 ;
		      if (cpe)
			{
			  *cpe = 0 ;
			  if (!strcmp (cp, "1"))
			    vtxtPrintf (blkp, " (e=10<sup>%s</sup>)", cpe+1) ;
			  else
			    vtxtPrintf (blkp, " (e=%s 10<sup>%s</sup>)", cp, cpe+1) ;
			}
		      else if (cp)
			vtxtPrintf (blkp, " (e=%s)", cp) ;
		    }
		}
	    }
	  if (isInteresting)
	    vtxtPrint (blkp, ", which may contain interesting functional annotation") ;
	}
      arrayDestroy (titles) ;
      arrayDestroy (expects) ;
      arrayDestroy (scores) ;
      dictDestroy (dict) ;
    }
  ac_free (h) ;
  ac_free (h0) ;
  return nn ; 
} /* ficheNewAceKogStatement */

/**************/

static char* ficheNewGeneMappingLinksAliasesParagraphContent (vTXT blkp, GMP *gmp)
{
  vtxtBold (blkp, "Map: ") ;

  vtxtPrintf (blkp, "This %sgene %s"
	      , gtIsEssentialGene (gmp->gene) ? "essential " : ""
	      , ac_name (gmp->gene)) ;
  ficheNewGeneMappingSentence (blkp, gmp, 0) ; 
  if (gmp->markup)
    ficheNewDbXrefParagraphContent (blkp, gmp) ;
 
  vtxtBreak (blkp) ;
  
  ficheNewGeneAliasStatement (blkp, gmp) ;
  ficheNewGeneEcNumberStatement  (blkp, gmp) ; 

  vtxtBreak (blkp) ;
  
  return vtxtPtr (blkp) ;
} /* ficheNewGeneMappingLinksAliasesParagraphContent */

/***************************************************************************************/

void ficheNewGeneMappingLinksAliasesSubSection (vTXT blkp, GMP *gmp)
{
  char *ptr ;
  vTXT bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  if ((ptr = ficheNewGeneMappingLinksAliasesParagraphContent (bfr, gmp)))
    {
      const char *chrom = gmp->gene ? ac_tag_printable (gmp->gene, "IntMap", 0) : 0 ;
      if (chrom && ! strncmp (chrom, "CHROMOSOME_", 11))
	chrom += 11 ;

     if (chrom && !strstr(chrom,"|"))
	gmpSubSection (blkp, gmp, "*tg_Alias_map"
		    , messprintf ("Map on chromosome %s, links to other databases and other names", chrom)
		    ) ;
      else
	gmpSubSection (blkp, gmp, "*tg_Alias_map", 
		    "Map, links to other databases and other names") ;
      vtxtPrint (blkp, ptr) ; 
      if (0) gmpChapterClose (blkp, gmp, "tg_Alias_map", TRUE) ;
    }
  vtxtDestroy (bfr) ;
} /* ficheNewGeneMappingLinksAliasesSubSection */

/***************************************************************************************/

void ficheNewMrnaStructureChapter (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  char *ptr = 0 ;
  vTXT bfr = vtxtHandleCreate (h) ; 
  vTXT bfr1 = vtxtHandleCreate (h) ; 
  if (gmp->markup) vtxtMarkup (bfr) ;
  if (gmp->markup) vtxtMarkup (bfr1) ;
  
  if (1)
    {
    gmpChapter (blkp, gmp,  "Gene_MOLECULAR", "MOLECULAR ANNOTATION of the mRNA") ;
    
    ficheNewGeneAnnotationOfVariantsTable (blkp, gmp, FALSE) ; 
    ficheNewGeneIntronsParagraph (blkp, gmp) ;
    if (0)ficheNewGenePolyAParagraph (blkp, gmp) ;
    
    gmpChapterClose (blkp, gmp, "Gene_MOLECULAR", TRUE) ;
    }

  if (0 && vtxtPrint (blkp, ptr))
    {
      const char *chrom = gmp->gene ? ac_tag_printable (gmp->gene, "IntMap", 0) : 0 ;
      if (chrom && !strstr(chrom,"|"))
	gmpChapter (blkp, gmp, "*tg_Alias_map"
		    , messprintf ("Map on chromosome %s, links to other databases and other names", chrom)
		    ) ;
      else
	gmpChapter (blkp, gmp, "*tg_Alias_map", 
		    "Map, links to other databases and other names") ;
      vtxtPrint (blkp, ptr) ; 
      gmpChapterClose (blkp, gmp, "tg_Alias_map", TRUE) ;
    }
ac_free (h) ;
} /* ficheNewGeneMappingLinksAliasesChapter */

/***************************************************************************************/

#ifdef JUNK
 junk until the titles become nice

static int ficheNewGeneProductMolecularNameStatement (vTXT blkp, GMP *gmp)
{
  int nProduct = 0 ;
  const char *ptr ;
  AC_TABLE gProd = ac_tag_table (gmp->gene, "Product", 0) ;

  /* product name */ 
  nProduct = gProd ? gProd->rows : 0 ;

  if (1)
    {
      vTXT bfr ; bfr = vtxtCreate () ;
      if (gmp->markup) vtxtMarkup (bfr) ;

      if (( ptr = gtProductMolecularName (bfr, gmp, gmp->product)) && ! strstr (ptr, "utative protein"))
	{
	  char *qptr = gtCleanUp (ptr) ;

	  vtxtDot (blkp) ;
	  if (nProduct > 1)
	    {
	      
	      vtxtPrint (blkp, "It encodes a set of ") ;
	      ptr = strstr (qptr, " member") ;
	      if (ptr && (ptr += 7) && *ptr != 's') 
		{
		  char cc = *ptr ;
		  *ptr = 0 ;
		  vtxtPrint (blkp, qptr) ;
		  vtxtPrint (blkp, "s ") ;
		  if (cc) 
		    {
		      *ptr = cc ;
		      vtxtPrint (blkp, ptr) ;
		    }
		}
	      else
		vtxtPrintf (blkp, "\"%s\"", qptr) ;
	    }
	  else
	    vtxtPrintf (blkp, "It encodes a%s %s"
			, (isVoyelle(*qptr) ? "n" : "")
			, qptr) ;
	}
      vtxtDestroy (bfr) ;
    }

  ac_free (gProd) ;
  return nProduct ;
} /* ficheNewGeneProductMolecularNameStatement */
#endif
/*****************/

static void ficheNewGenePseudogeneStatement (vTXT blkp, GMP *gmp)
{
  if (0 &&
      ac_has_tag (gmp->gene, "Pseudogene"))
    vtxtPrintf (blkp, "It is likely a pseudogene: it has no intron and the cdna clones differ a lot from the genome sequence.<br/>\n") ;
} /* ficheNewGenePseudogeneStatement */

/*****************/

static int ficheNewGeneProductPropertiesStatement (vTXT blkp, GMP *gmp)
{
  const char *ptr ;
  int ir, nn = 0 ;
  AC_TABLE gProperties = 0 ;
  
  if (gmp->Spc == WORM &&
      (gProperties = ac_tag_table (gmp->gene, "Properties", 0)))
    {
      for (ir=0 ; ir < gProperties->rows ; ir++)
	{
	  if ((ptr = ac_table_printable (gProperties, ir, 0, 0)))
	    {
	      if (!nn++)
		{
		  vtxtBreak (blkp) ;
		  vtxtBold (blkp, "Molecular properties: ") ;
		}
	      else
		vtxtBreak (blkp) ;
	      vtxtPrint (blkp, gtCleanQuotes (ptr)) ;
	    }
        }
    }
  if (nn) 
    vtxtBreak (blkp) ;
  ac_free (gProperties) ;
  return nn ;
} /* ficheNewGeneProductPropertiesStatement */

/***************************************************************************************/

static int ficheNewGenePathwaysProcessFunctionLocalizationTableNew (vTXT blkp, GMP *gmp, Array bb)
{
  AC_HANDLE h = ac_new_handle () ;
  KEYSET ksGPap = keySetHandleCreate (h) ;
  Array aa = arrayHandleCreate (256, ITRC, h) ;
  AC_TABLE tbl ; 
  vTXT bfr = vtxtHandleCreate (h) ; 
  vTXT bfrPap = vtxtHandleCreate (h) ; 
  int nGwithPap = 0, nnGwithPap = 0 ; 
  int row = 0, cols[7] = {1,2,3,4,5,0, 0} ;
  enum { COL_TYPE=0, COL_MESH, COL_EVIDENCE, COL_SOURCE, COL_GENES, COL_VOTE, COL_LAST} ; 
  const char *colNames[]={ "Type", 
			   "Description", 
			   "Evidence",
			   "Source",
			   "Related genes (* if published support)",
			   0} ; 

  /* contruct the table */
  tbl = ac_empty_table (80, 6, h) ; /* tbl->rows is a hint, not a hard limit */

  /* format the table and add the http links */
  nGwithPap = 0, nnGwithPap = 0 ;  
  vtxtClear (bfrPap) ;
  vtxtPrint (bfrPap, PUBMED_LINK) ; /* no %s included */ 

  ficheNewGeneDiseaseKeggPathwaysTableNew (tbl, gmp, aa) ;
  ficheNewGeneProcessFunctionLocalizationTableNew (tbl, gmp, bfrPap, &nGwithPap, &nnGwithPap, ksGPap, 0, aa) ;
  ficheNewGeneProcessFunctionLocalizationTableNew (tbl, gmp, bfrPap, &nGwithPap, &nnGwithPap, ksGPap, 1, aa) ;
  
  if (aa && arrayMax (aa))
    {
      int nn, row, ir, ir2 ;
      ITRC *itrc, *itrc2 ;
      int nItrc1 = arrayMax (aa), nItrc2 ;
      AC_OBJ Gene = 0 ;
      double cc, oldcc ;

      arraySort (aa, itrcOrder) ;  
      if (bb)
	for (nItrc2 = 0, ir2 = arrayMax (bb) ; nItrc2 < nItrc1 ; ir2++, nItrc2++)
	  {
	    itrc = arrp (aa, nItrc2, ITRC) ;
	    if (itrc->nn)
	      {
		itrc2 = arrayp (bb, ir2, ITRC) ;
		*itrc2 = *itrc ; 
		itrc2->Gene = 0 ;
	      }
	  }
	  
      itrc = arrp (aa, 0, ITRC) ;
      nn = itrc->nn ;
      if (nn > 1)
	{
	  row = tbl->rows ;
	  ac_table_insert_text (tbl, row, COL_MESH
			    , "<font color='red'>Genes most related through pathways, process or function (with published evidence)</font>") ;
      
	  cc = oldcc = itrc->cc ;
	  for (nItrc2 = ir = 0 ; nItrc2 < nItrc1 ; nItrc2++)
	    {
	      itrc = arrp (aa, nItrc2, ITRC) ;
	      if (2*itrc->cc >= cc) ir++ ;
	    }
	  
	  vtxtClear (bfr) ; 
	  for (nItrc2 = ir2 = 0 ; nItrc2 < nItrc1 ; nItrc2++)
	    {
	      itrc = arrp (aa, nItrc2, ITRC) ;
	      if (2*itrc->cc >= cc || FDEBUG)
		{
		  if (itrc->cc < oldcc)
		    vtxtPrint (bfr, ir2++ ? " | " : "") ;
		  else
		    vtxtPrint (bfr, ir2++ ? ", " : "") ;
		  oldcc = itrc->cc ;
		  Gene = ac_get_obj (gmp->db, "gene", name(itrc->gene), 0) ;
		  gmpObjLink (bfr, gmp, Gene, ac_name (Gene)) ;
		  ac_free (Gene) ;
	      
		  if (ir2 > 12 &&  ir > 20 && ! FDEBUG)
		    { vtxtPrint (bfr, "...") ; break ; }
		}
	    }
	  ac_table_insert_text (tbl, row, COL_GENES, vtxtPtr (bfr)) ;
	}
    }
  aa = arrayReCreate (aa, 256, ITRC) ;
  ficheNewGeneProcessFunctionLocalizationTableNew (tbl, gmp, bfrPap, &nGwithPap, &nnGwithPap, ksGPap, 2, aa) ;

   if (aa && arrayMax (aa) && bb)
    {
      int ir2 ;
      ITRC *itrc, *itrc2 ;
      int nItrc1 = arrayMax (aa), nItrc2 ;
      
      for (nItrc2 = 0, ir2 = arrayMax (bb) ; nItrc2 < nItrc1 ; ir2++, nItrc2++)
	{
	  itrc = arrp (aa, nItrc2, ITRC) ;
	  if (itrc->nn)
	    {
	      itrc2 = arrayp (bb, ir2, ITRC) ;
	      *itrc2 = *itrc ; 
	      itrc2->Gene = 0 ;
	    }
	}
    }
   if (nGwithPap > 1)
    {
      int row = tbl->rows ;
      int ir = nnGwithPap ;
      ac_table_insert_text (tbl, row, COL_MESH, "Cumulated literature") ;

      vtxtClear (bfr) ; 
      gmpURL (bfr, gmp, vtxtPtr (bfrPap)
	      , messprintf ("%d article%s", ir, _multi(ir))) ;
      ac_table_insert_text (tbl, row, COL_EVIDENCE, vtxtPtr (bfr)) ;
    }

  /* export */
  if ((row = tbl->rows))
    {
      vtxtPrint (blkp, "This section summarizes the functional aspects: pathways, processes, molecular function, enzymatic activity, or localization of the protein(s) to cell compartments. Some annotations are documented in PubMed, some are inferred. The lists of related genes with the same process or function GO annotation are reported in the last column only if they are supported by a PubMed publication") ;
      vtxtBreak (blkp) ;
      
      ac_table_display (blkp 
			, tbl, colNames
			, cols, 7
			, 0, 0, 0
			, 0
			) ;
    }

  ac_free (h) ; 

  return row ;
} /* ficheNewGeneGoLocalizationTableNew */

/*****************/

static void ficheNewGenePathwaysProcessFunctionLocalizationParagraph (vTXT blkp, GMP *gmp, Array bb)
{
  char *ptr ;
  vTXT bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  ficheNewGenePathwaysProcessFunctionLocalizationTableNew (bfr, gmp, bb) ;
  ficheNewGeneProductPropertiesStatement (bfr, gmp) ; /* WORM gene->Properties */

  if ((ptr = vtxtPtr (bfr)))
    {
      gmpSection (blkp, gmp, "tg_pathway", "Pathways, biological processes, molecular function and cellular localization (GO)") ;
      vtxtPrint (blkp, ptr) ; 
    }
  vtxtDestroy (bfr) ;
} /* ficheNewGenePathwaysProcessFunctionLocalizationParagraph */

/***************************************************************************************/

static void ficheNewGeneProductPfamMotifTableTable (AC_TABLE tbl1, GMP *gmp, Array aa)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *ccp, *txtAccession ;
  char linkBuf[2000] ;
  int ir, jr, iType, np, row ;
  AC_OBJ oPfam, gene ;
  AC_TABLE tbl = 0 ;
  AC_ITER iter ;
  AC_KEYSET otherPfam, otherGene ;
  vTXT bfr = vtxtHandleCreate (h) ;  
  vTXT bfr1 = vtxtHandleCreate (h) ;  
  enum { COL_TYPE=0, COL_MESH, COL_EVIDENCE, COL_SOURCE, COL_GENES, COL_VOTE, COL_LAST} ; 

  if (gmp->markup) vtxtMarkup (bfr) ;   
  if (gmp->markup) vtxtMarkup (bfr1) ;   
  
  iType = 0 ; 
  
  iter = ac_objquery_iter (gmp->gene, ">Pfam", h) ;
  oPfam = 0 ; vtxtClear (bfr) ;
  while (ac_free (oPfam), (oPfam = ac_iter_obj (iter)))
    {
      row = tbl1->rows ;
      if (!iType++)
	ac_table_insert_text (tbl1, row, COL_TYPE, "Motif") ;
      
      vtxtClear (bfr) ;
      vtxtPrintf (bfr, "The ") ;
      ccp = ac_tag_printable (oPfam, "Definition", ac_name (oPfam)) ;
      if ((txtAccession = ac_tag_printable (oPfam, "Accession", 0)))
	{
	  char *cq = strstr (txtAccession , ".") ; /* drop pfam version number */
	  if (cq) *cq = 0 ;
	  sprintf (linkBuf, "http://pfam.xfam.org/family/%s", txtAccession) ; 
	  gmpURL (bfr, gmp, linkBuf, ccp) ; 
	}
      else
	vtxtPrint (bfr, ccp) ; 
      
      ac_table_insert_text (tbl1, row, COL_SOURCE, "AceView") ;

      ccp = messprintf (">product best_product && good_product || very_good_product ; pfam = %s",
			ac_protect (ac_name(oPfam), h)) ;
      otherPfam = ac_objquery_keyset (gmp->gene, ccp, h) ;
      np = otherPfam ? ac_keyset_count (otherPfam) : 0 ;
      if (np > 1)
	vtxtPrintf (bfr, " domain is found in %d isoforms", np) ;
      else if (np)
	vtxtPrintf (bfr, " domain is seen in isoform") ;
      if (np)
	{
	  tbl = ac_keyset_table (otherPfam, 0, -1, 0, h) ;
	  for (ir = jr = 0 ; tbl && ir < tbl->rows ; ir++)
	    {
	      gene = ac_table_obj (tbl, ir, 0, h) ; 
	      
	      if (1 || ! ac_obj_equal (gene, gmp->gene))
		{
		  vtxtPrint (bfr, jr++ ? ", " : ": ") ;
		  gmpObjLink (bfr, gmp, gene, gtMrnaSuffix(ac_name (gmp->gene), ac_name (gene), h)) ;
		}
	      if (jr > 2 &&  tbl->rows > 8)
		{ vtxtPrint (bfr, "...") ; break ; }
	    }
	}

      if (gmp->markup || gmp->style == 'r')
	{
	  const char *ptr ;
	  
	  if (
	      (ptr = ac_tag_printable (oPfam, "Comment", 0)) ||
	      (ptr = ac_tag_printable (oPfam, "Interpro_comment", 0)) ||
	      (ptr = ac_tag_printable (oPfam, "pfam_comment", 0))
	      )
	    {
	      vtxtBreak (bfr) ;
	      vtxtBold (bfr, "[InterPro annotation]") ;
	      vtxtPrintf (bfr, "  %s", ptr) ;
	    }
	  else /* try to get it from a synonymous pfam */
	    {
	      AC_KEYSET ks = ac_objquery_keyset (oPfam, ">accession; >quoted_in; Interpro_comment", h) ;
	      AC_ITER iter2 = ac_keyset_iter (ks, TRUE, h) ;
	      AC_OBJ oPfam2 = ac_next_obj (iter2) ;
	      
	      if ((ptr = ac_tag_printable (oPfam2, "Interpro_comment", 0)))
		{
		  vtxtBreak (bfr) ;
		  vtxtBold (bfr, "<span class='explain2'>[InterPro annotation]") ;
		  vtxtPrintf (bfr, "  %s</span>", ptr) ;
		}
	      ac_free (ks) ; ac_free (oPfam2) ; ac_free (iter2) ;
	    }
	}
      ac_table_insert_text (tbl1, row, COL_MESH, vtxtPtr (bfr)) ;
      
      ac_table_insert_text (tbl1, row, COL_LAST,  "2") ;
      
      vtxtClear (bfr) ;
      vtxtClear (bfr1) ; jr = 0 ;

      if ((otherGene = ac_objquery_keyset (oPfam, ">Gene", h)))
      { 
	int nGenes, nItrc1 = arrayMax (aa), nItrc2 ; 
	KEY kU ;
	ITRC *itrc ;
	double cc = 0 ;
	BOOL show = TRUE ;
	KEY k0 = ac_obj_key (gmp->gene) ;
	AC_OBJ oG ;

	nGenes = ac_keyset_count (otherGene) ;
	if (nGenes == 1)
	  {
	    vtxtPrintf (bfr, "No other gene in the database contains this motif") ; 
	  }
	else  if (nGenes > 1)
	  {	
	    gmpObjLink (bfr, gmp, oPfam, messprintf ("%d genes: ", nGenes)) ; 
	    
	    ir = nGenes ;
	    if (ir < 2) ir = 2 ;
	    if (ir < 10000) cc = - log(((double)ir)/10000.0) ;
	    cc *= Npfam ;
	    
	    tbl = ac_keyset_table (otherGene, 0, -1, 0, h) ;
	    for (ir = 0 ; ir < tbl->rows ; ir++)
	      {
		kU = ac_table_key (tbl, ir, 0, 0) ; 
		if (kU != k0)
		  {
		    for (nItrc2 = 0 ; nItrc2 < nItrc1 ; nItrc2++)
		      {
			itrc = arrp (aa, nItrc2, ITRC) ;
			if (itrc->gene == kU)
			  { itrc->nn++ ; itrc->cc += cc  ; itrc->type = Mtype ; break ; }
		      }
		    if (nItrc2 >= nItrc1)
		      {
			itrc = arrayp (aa, nItrc1++, ITRC) ;
			itrc->gene = kU ;
			itrc->Gene = 0 ;
			itrc->nn = 1 ;
			itrc->cc = cc ;
			itrc->type = Mtype ;
		      }
		  }
		if (show)
		  {
		    oG = ac_table_obj (tbl, ir, 0, h) ;
		    if (ir) vtxtPrintf (bfr, ", ") ; 
		    gmpObjLink (bfr, gmp,  oG, ac_name (oG)) ;
		    ac_free (oG) ;
		    if (ir > 8 && tbl->rows > 12)
		      { vtxtPrintf (bfr, "...") ; break ; }
		  }
	      }
	      
	    if (gmp->Spc == WORM)
	      {
		int nTg = ac_keyset_count (ac_objquery_keyset (oPfam, "follow gene ; Transcribed_gene || Expression_title || Pattern || Expr_pattern || Has_OST || Has_cDNA_clone", h)) ;
		if (gmp->tg || ac_has_tag (gmp->gene, "Expression"))
		  nTg-- ;
		if (nTg <= 0)
		  vtxtPrintf (bfr, ", none of which is known to be expressed") ;
		else if (nGenes > 2 && nTg == 1)
		  vtxtPrintf (bfr, " and only one is known to be expressed") ;
		else if (nTg == 1)
		  vtxtPrintf (bfr, " and is known to be expressed") ;
		else if (nTg == nGenes - 1)
		  vtxtPrintf (bfr, ", all of which are known to be expressed") ;
		else if (nTg > 0)
		  vtxtPrintf (bfr, ", %d of which are known to be expressed", nTg) ;
	      }
	  }
	ac_free (tbl) ;
	ac_table_insert_text (tbl1, row, COL_GENES, vtxtPtr (bfr)) ;
      }
    }
  ac_free (oPfam) ;
  ac_free (h) ;
  return ;
} /* ficheNewGeneProductPfamMotifTableTable */

/***************************************************************************************/

static void ficheNewGeneProductPsortMotifTableTable (AC_TABLE tbl1, GMP *gmp, Array aa)
{
  AC_HANDLE h = ac_new_handle () ;
  int jj, nn, row, iProd, iMotif ;
  KEY k0 = ac_obj_key (gmp->gene) ;
  AC_OBJ oPsort, oProd ;
  AC_TABLE gPsort = 0, gProds ;
  vTXT bfr = vtxtHandleCreate (h) ;  
  enum { COL_TYPE=0, COL_MESH, COL_EVIDENCE, COL_SOURCE, COL_GENES, COL_VOTE, COL_LAST} ; 
  char **cpp, *psNN [] =
    { 
      "N_terminal_signal_domain", "an N-terminal peptide signal",

      "Transmembrane_domain", "a transmembrane domain", 
      "N_myristoylation_domain", "a N-myristoylation domain",
      "prenylation_domain", "a prenylation domain",
      "Golgi_transport_domain", "a Golgi transport domain", 
      
      "peroxisomal_domain", "a peroxisomal domain", 
      "2nd_peroximal_domain", "a second peroximal domain", 
      "vacuolar_domain", "a vacuolar domain", 
      
      "ER_retention_domain", "an ER_retention domain", 
      "Coiled_coil_region", "a coiled coil stretch", 
      "Leucine_zipper_domain", "A leucine zipper domain",
      "actin_binding_1_domain", "an actin binding domain", 
      "actin_binding_2_domain", "an actin binding domain", 
      "ATP_binding_domain", "an ATP binding domain",   

      "RNA_binding_domain", "an RNA binding domain",                           
      0, 0
    } ;

  if (gmp->markup) vtxtMarkup (bfr) ;
  
  gProds = ac_tag_table (gmp->gene, "Product", h) ;

  for (iMotif = 0, cpp = psNN ; gProds && *cpp ; iMotif += 2, cpp += 2 )
    {
      row = tbl1->rows ;
      vtxtClear (bfr) ;
      oPsort = ac_get_obj (gmp->db, "Psort", *cpp, h) ;

      for (iProd = jj = 0,  gPsort = 0 ; iProd < gProds->rows ; ac_free (oProd), ac_free (gPsort), iProd++)
	{
	  oProd = ac_table_obj (gProds, iProd, 0, h) ;
	  if (
	      (!ac_has_tag (oProd, "Best_product") && !ac_has_tag (oProd, "Very_good_product") ) ||
	      !ac_has_tag (oProd, "Good_product")
	      )
	    continue ;
	  gPsort = ac_tag_table (oProd, *cpp, h) ;
	  if (!gPsort)
	    continue ;

	  if (!jj++)
	    {
	      gmpURL (bfr, gmp, PSORT_LINK, *(cpp+1)) ;
	      vtxtPrintf (bfr, " is seen") ;
	    }
	  else
	    vtxtPrintf (bfr, ", ") ;	    
	  nn = gPsort->rows ;
	  if (nn > 1)
	    vtxtPrintf (bfr, " %d times", nn) ;
	  else if (nn == 1)
	    vtxtPrintf (bfr, " once") ;
	  if (jj == 1)
	    vtxtPrint (bfr, " in isoform ") ;
	  else
	    vtxtPrint (bfr, " in ") ;
	  gmpObjLink (bfr, gmp, oProd, gtMrnaSuffix(ac_name (gmp->gene), ac_name(oProd), h)) ;
	}
      if (!jj)
	continue ;
      ac_table_insert_text (tbl1, row, COL_MESH, vtxtPtr (bfr)) ;
      ac_table_insert_text (tbl1, row, COL_SOURCE, "AceView") ;

      ac_table_insert_text (tbl1, row, COL_LAST,  "3") ;

      { 
	int nItrc1 = arrayMax (aa), i, nItrc2 ; 
	KEY kU ;
	KEYSET ks1, ks2 ;
	ITRC *itrc ;
	double cc = 0 ;

	ks1 = query(0, messprintf ("Find Product %s", *cpp)) ;
	ks2 = ks1 ? OWQProduct2Gene (ks1, FALSE) : 0 ;
	nn = ks2 ? keySetMax (ks2) : 0 ;
	
	i = nn ;
	if (i < 2) i = 2 ;
	if (i < 10000) cc = - log(((double)i)/10000.0) ;
	cc *= Npsort ;

	for (i = 0 ; ks2 && i < keySetMax (ks2) ; i++)
	  {
	    kU = keySet (ks2, i) ; 
	    if (kU != k0)
	      {
		for (nItrc2 = 0 ; nItrc2 < nItrc1 ; nItrc2++)
		  {
		    itrc = arrp (aa, nItrc2, ITRC) ;
		    if (itrc->gene == kU)
		      { itrc->nn++ ; itrc->cc += cc  ; itrc->type = Mtype ; break ; }
		  }
		if (nItrc2 >= nItrc1)
		  {
		    itrc = arrayp (aa, nItrc1++, ITRC) ;
		    itrc->gene = kU ;
		    itrc->Gene = 0 ;
		    itrc->nn = 1 ;
		    itrc->cc = cc ;
		    itrc->type = Mtype ;
		  }
	      }
	  }
	keySetDestroy (ks1) ;
	keySetDestroy (ks2) ;
      }
      if (nn > 1)
	{  
	  vtxtClear (bfr) ;
	  gmpObjLink (bfr, gmp, oPsort, messprintf ("%s gene%s", isOne (nn) , _multi (nn))) ;

	  ac_table_insert_text (tbl1, row, COL_GENES, vtxtPtr (bfr)) ;
	}
    }

  ac_free (h) ;
  return ;
} /* ficheNewGeneProductPsortMotifTableTable */

/***************************************************************************************/

static int ficheNewGeneProductPfamPsortMotifTableNew (vTXT blkp, GMP *gmp, Array bb)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl ; 
  Array aa = arrayHandleCreate (256, ITRC, h) ;
  vTXT bfr = vtxtHandleCreate (h) ; 
  vTXT bfrPap = vtxtHandleCreate (h) ; 
  int nGwithPap = 0, nnGwithPap = 0 ; 
  int row = 0, cols[6] = {2,5,0} ;
  enum { COL_TYPE=0, COL_MESH, COL_EVIDENCE, COL_SOURCE, COL_GENES, COL_VOTE, COL_LAST} ; 
  const char *colNames[]={ "Type", 
			   "Protein domains or motifs", 
			   "Evidence",
			   "Source",
			   "Genes with same motif",
			   0} ; 

  /* contruct the table */
  tbl = ac_empty_table (80, 6, h) ; /* tbl->rows is a hint, not a hard limit */

  /* format the table and add the http links */
  nGwithPap = 0, nnGwithPap = 0 ;  
  vtxtClear (bfrPap) ;
  vtxtPrint (bfrPap, PUBMED_LINK) ; /* no %s included */ 
  if (nGwithPap > 1)
    {
      int row = tbl->rows ;
      int ir = nnGwithPap ;
      ac_table_insert_text (tbl, row, COL_MESH, "Cumulated literature") ;

      vtxtClear (bfr) ; 
      gmpURL (bfr, gmp, vtxtPtr (bfrPap)
	      , messprintf ("%d article%s", ir, _multi(ir))) ;
      ac_table_insert_text (tbl, row, COL_EVIDENCE, vtxtPtr (bfr)) ;
    }

  nGwithPap = 0, nnGwithPap = 0 ;  
  vtxtClear (bfrPap) ;
  vtxtPrint (bfrPap, PUBMED_LINK) ; /* no %s included */ 

  ficheNewGeneProductPfamMotifTableTable (tbl, gmp, aa) ;
  ficheNewGeneProductPsortMotifTableTable (tbl, gmp, aa) ;
  
  if (arrayMax (aa))
    {
      int nn, row, ir, ir2 ;
      ITRC *itrc, *itrc2 ;
      int nItrc1 = arrayMax (aa), nItrc2 ;
      double cc, oldcc ;
      AC_OBJ Gene ;
      BOOL show ;

     arraySort (aa, itrcOrder) ;  
      for (nItrc2 = 0, ir2 = arrayMax (bb) ; nItrc2 < nItrc1 ; ir2++, nItrc2++)
	{
	  itrc = arrp (aa, nItrc2, ITRC) ;
	  itrc2 = arrayp (bb, ir2, ITRC) ;
	  *itrc2 = *itrc ;
	  itrc2->Gene = 0 ;
	}
	  
      itrc = arrp (aa, 0, ITRC) ;
      nn = itrc->nn ;
      if (nn > 1)
	{
	  row = tbl->rows ;
	  ac_table_insert_text (tbl, row, COL_MESH
			    , "<font color='red'>Genes encoding proteins with most motifs in common</font>") ;
      
	  cc = oldcc = itrc->cc ;
	  for (nItrc2 = ir = 0 ; nItrc2 < nItrc1 ; nItrc2++)
	    {
	      itrc = arrp (aa, nItrc2, ITRC) ;
	      if (2*itrc->cc >= cc) ir++ ;
	    }
	  
	  vtxtClear (bfr) ; 
	  for (nItrc2 = ir2 = 0, show = TRUE ; show && nItrc2 < nItrc1 ; nItrc2++)
	    {
	      itrc = arrp (aa, nItrc2, ITRC) ;
	      if (2*itrc->cc >= cc || FDEBUG)
		{
		  if (itrc->cc < oldcc)
		    vtxtPrint (bfr, ir2++ ? " | " : "") ;
		  else
		    vtxtPrint (bfr, ir2++ ? ", " : "") ;
		  oldcc = itrc->cc ;
		  Gene = ac_get_obj (gmp->db, "gene", ac_key_name(itrc->gene), 0) ;
		  gmpObjLink (bfr, gmp, Gene, ac_name (Gene)) ;
		  ac_free (Gene) ;
		  if (ir2 > 8 &&  ir > 12 && ! FDEBUG)
		    { vtxtPrint (bfr, "...") ; show = FALSE ; }
		}
	    }
	  ac_table_insert_text (tbl, row, COL_GENES, vtxtPtr (bfr)) ;
	}
    }
  
  /* export */
  if ((row = tbl->rows))
    ac_table_display (blkp 
		      , tbl, colNames
		      , cols, 2
		      , 0, 0, 0
		      , 0
		      ) ;

  ac_free (h) ; 

  return row ;
} /* ficheNewGeneProductPfamPsortMotifTableNew */

/*****************/

static void ficheNewGeneProductPfamPsortParagraph (vTXT blkp, GMP *gmp, Array bb)
{
  char *ptr ;
  vTXT bfr = vtxtCreate () ; 

  if (gmp->markup) vtxtMarkup (bfr) ;

  ficheNewGenePseudogeneStatement  (bfr, gmp) ; 
  ficheNewGeneProductPropertiesStatement (bfr, gmp) ; /* gene->Properties */
  ficheNewGeneProductPfamPsortMotifTableNew  (bfr, gmp, bb) ;

  if ((ptr = vtxtPtr (bfr)))
    {
      gmpSection (blkp, gmp, "tg_product_pfam", "Protein domains and motifs") ;
      vtxtPrint (blkp, ptr) ; 
    }
  vtxtDestroy (bfr) ;
} /* ficheNewGeneProductPfamParagraph */

/***************************************************************************************/

static void  ficheNewGeneLocus_DescriptionStatement (vTXT blkp, GMP *gmp)
{
  int ir, jr=0 ;
  const char *ptr ;
  AC_TABLE oTmp ;
  AC_HANDLE h = ac_new_handle () ;

  /* Locus_Description */
  if ((oTmp = ac_tag_table (gmp->gene, "Locus_Description", h)))
    {
      for (ir=jr=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
	{
	  ptr = ac_table_printable (oTmp, ir, 0, 0) ;
	  if (ptr && *ptr && !strstr(ptr,"ssential gene")) /* essential gene is reported separately */
	    {
	      if (!jr++)
		{
		  vtxtDot (blkp) ;
		  if (oTmp->rows > 4)
		    vtxtPrintf (blkp, "Phenotypes and affected processes are") ;
		  else
		    vtxtPrintf (blkp, "Its phenotype is") ;
		}
	      else
		vtxtPrintf (blkp, ", ") ;
	      vtxtPrintf (blkp, " %s", gtLowerCleanUp(ptr)) ;
	    }
	}
    }
  ac_free (h) ;
} /*ficheNewGeneLocus_DescriptionStatement */

/*************/

static void  ficheNewGeneOmimStatement (vTXT blkp, GMP *gmp)
{
  int ir, jr=0 ;
  const char *ccp ;
  char linkBuf[2000] ;
  AC_TABLE oTmp ;
  AC_HANDLE h = ac_new_handle () ;

  if (0 &&  /* pas interessant, Danielle 2007_02_20 */
      (oTmp = ac_tag_table (gmp->gene, "Extern OMIM", h)))
    {
      for (ir=0;ir < oTmp->rows && oTmp->cols >= 1;ir++)
	{
	  ccp = ac_table_printable (oTmp, ir, 0, 0) ;
	  if (!jr++)
	    {
	      vtxtDot (blkp) ;
	      vtxtPrintf (blkp, "Interesting annotations are available from OMIM ") ;
	    }
	  else
	    vtxtPrintf (blkp, ", ") ;
	  sprintf (linkBuf, OMIM_LINK, ccp + 5) ; 
	  gmpURL (blkp, gmp, linkBuf, ccp + 5) ;
	}
    }
 ac_free (h) ;
} /* ficheNewGeneOmimStatement */

/***************************************************************************************/

static int ficheNewGeneLocus_PhenotypeStatement (vTXT blkp, GMP *gmp)
{
  int ir, nn = 0 ;
  AC_TABLE gLocus_Phenotype ;
  const char *txt ;
  char linkBuf[vONLINEMAXURLSIZE] ;
  int omim ;
  AC_HANDLE h = ac_new_handle () ;

  if ((gLocus_Phenotype = ac_tag_table (gmp->gene, "Locus_Phenotype", h)))
    {
      for (ir=0 ; ir < gLocus_Phenotype->rows ; ir++)
	{
	  txt = ac_table_printable (gLocus_Phenotype, ir, 0, 0) ; 
	  if (txt)
	    {
	      nn++ ; 
	      vtxtBreak (blkp) ;
	      
	      if (sscanf (txt, "[OMIM %d]", &omim) == 1 ||
		  sscanf (txt, "[MIM %d]", &omim) == 1) 
		{ /* Locus_Phenotype may look like [OMIM 25070] */
		  char cc = 0, *ptr = strstr (txt, "]") ;
		  
		  if (ptr) { ptr++ ; cc = *ptr ; *ptr = 0 ; }
		  sprintf (linkBuf, OMIM_LINK, messprintf ("%d", omim)) ;
		  gmpURL (blkp, gmp, linkBuf, txt)  ;
		  if (ptr) { txt = ptr ; *ptr = cc ; } /* jump and restore */
		  if (*txt != ' ') vtxtPrint (blkp, " ") ;
		}
	      vtxtPrint (blkp, gtCleanQuotes (txt)) ;
	      vtxtBreak (blkp) ;
	    }
	}
    }
  
  ac_free (h) ;
  return nn ;
} /*ficheNewGeneLocus_PhenotypeStatement */

/***************************************************************************************/

static vTXT ficheNewGeneFunctionLocusPhenotypeDescription (GMP *gmp, int *nnp, AC_HANDLE h0)
{
  int ir, nn = 0 ;
  AC_TABLE gLocus_description ;
  const char *ccp ;
  BOOL isEssential = FALSE ;
  vTXT txt = 0, txt1 = 0 ;
  AC_HANDLE h = 0 ;

  switch (gmp->Spc)
    {
    case ARA:
      if (ac_has_tag (gmp->gene, "Locus_description") ||
	  ac_has_tag (gmp->gene, "Locus_phenotype")
	  )
	{
	  nn = 1 ;
	  txt = vtxtHandleCreate (h0) ;
	  vtxtPrint (txt, "This gene was associated to a <a href='javascript:openAnchor (\"ffunc\",\"Gene_FUNCTION\")'>phenotype</a> by TAIR") ;
	}
      break ;
    case WORM:
      if (!ac_has_tag (gmp->gene, "Locus_description") &&
	  ac_has_tag (gmp->gene, "Locus_phenotype")
	  )
	{
	  nn = 1 ;
	  txt = vtxtHandleCreate (h0) ;
	  vtxtPrint (txt, "This gene is associated to a <a href='javascript:openAnchor (\"ffunc\",\"Gene_FUNCTION\")'>phenotype</a>") ;
	}
      else if (ac_has_tag (gmp->gene, "Locus_description"))
	{
	  txt = vtxtHandleCreate (h0) ;
	  h = ac_new_handle () ;
	  txt1 = vtxtHandleCreate (h) ;
	  nn = 0 ;
	  gLocus_description = ac_tag_table (gmp->gene, "Locus_description", h) ;
	  for (ir=0 ; ir < gLocus_description->rows ; ir++)
	    {
	      ccp = ac_table_printable (gLocus_description, ir, 0, 0) ; 
	      if (strstr (ccp, "essential gene"))
		isEssential = TRUE ;
	      else
		{
		  if (nn++)
		    vtxtComma (txt1) ;
		  vtxtPrint (txt1, gtCleanQuotes (ccp)) ;
		}
	    }
	  if (nn)
	    vtxtPrintf (txt, "This %sgene is associated to a <a href='javascript:openAnchor (\"ffunc\",\"Gene_FUNCTION\")'>phenotype</a> (%s)"
			, isEssential ? "essential " : ""
			, vtxtPtr (txt1)
			) ;
	}
      break ;
    default:
      break ;
    }
  
  ac_free (h) ;
  *nnp = nn ;
  return txt ;
} /* ficheNewGeneFunctionLocusPhenotypeDescription */

/*************/

static int ficheNewGeneRNAiStatement (vTXT blkp, GMP *gmp)
{
  const char *ptr ;
  char linkBuf[vONLINEMAXURLSIZE] ;
  int ir, jr, nn = 0, nn2 ;
  AC_TABLE gRnai, gGene, gPheno ;
  AC_OBJ oRnai, oGene ;
  AC_HANDLE h = ac_new_handle () ;

  if ((gRnai = ac_tag_table (gmp->gene, "RNAi", h)))
    for (ir=0; ir < gRnai->rows && gRnai->cols >= 1; ir++)
      {
	oRnai = ac_table_obj (gRnai, ir, 0, h) ;
	nn2 = 0 ;
	if (!strncmp("mv_",ac_name (oRnai),3))
	  continue ; 
	if ((gPheno = ac_tag_table (oRnai, "Phenotype", h)) && ac_tag_printable (oRnai, "Phenotype", 0))
	  {
	    for (jr=0; jr < gPheno->rows && gPheno->cols >= 1; jr++)
	      if ((ptr = ac_table_printable(gPheno, jr, 0, 0)))
		{
		  if (!nn2++)
		    {
		      if (!nn++)
			{
			  vtxtBreak (blkp) ; /* vtxtEmptyLine (blkp, 1) ; */
			  vtxtBold (blkp, "RNA interference results ") ;
			}
		      else
			vtxtBreak (blkp) ;
		    }
		  else
		    vtxtPrintf (blkp, ". ") ;
		  vtxtPrint (blkp, gtSetUpper(gtCleanUp(ptr))) ;
		}
	    if (strstr(ac_name (oRnai), "JA"))
	      {
		vtxtPrintf (blkp, " (by feeding genomic PCR product %s", ac_name (oRnai)) ;
		if (gmp->markup)
		  vtxtPrintf (blkp, ", shown as a red box on the graph") ;
		vtxtPrintf (blkp, ")") ;
	      }
	    else if (strstr(ac_name (oRnai), "TH"))
	      {
		vtxtPrintf (blkp, " (by injecting genomic PCR product %s", ac_name (oRnai)) ;
		if (gmp->markup)
		  vtxtPrintf (blkp, ", shown as a red box on the graph") ;
		vtxtPrintf (blkp, ")") ;
	      }
	    else if (strstr(ac_name (oRnai), "SA") || strstr(ac_name (oRnai), "FP"))
	      vtxtPrintf (blkp, " (by injecting cDNA clone %s)", ac_name (oRnai)) ;

	  }
	else if (ac_has_tag (oRnai, "No_obvious_phenotype"))
	  {
	    if (!nn2++)
	      {
		if (!nn++)
		  gmpSection (blkp, gmp, "RNAi", "RNA interference results") ;
		else
		  vtxtBreak (blkp) ;
	      }
	    ptr = ac_tag_printable (oRnai, "No_obvious_phenotype", "") ;
	    vtxtPrintf (blkp, "%s No obvious phenotype", ptr) ;
	    
	    if (strstr(ac_name (oRnai), "JA"))
	      {
		vtxtPrintf (blkp, " (by feeding genomic PCR product %s)", ac_name (oRnai)) ;
		if (gmp->markup)
		  vtxtPrintf (blkp, ", shown as a grey box on the graph") ;
	      }
	    else if (strstr(ac_name (oRnai), "TH"))
	      vtxtPrintf (blkp, " (by injecting genomic PCR product %s)", ac_name (oRnai)) ;
	    else if (strstr(ac_name (oRnai), "SA") || strstr(ac_name (oRnai), "FP"))
	      vtxtPrintf (blkp, " (by injecting cDNA clone %s)", ac_name (oRnai)) ;
	  }
	if (ac_has_tag (oRnai, "DIC_Movie_available"))
	  {
	    if (!strncasecmp (ac_name (oRnai), "TH", 2))
	      gmpURL (blkp, gmp
		      , "http://worm-srv1.mpi-cbg.de/dbScreen/"
		      , ". Movies are available on Hyman's site"
		      ) ;
	    else /* if (!strncasecmp (ac_name (oRnai), "FP", 2)) */
	      {
		sprintf (linkBuf
			 , "http://nematoda.bio.nyu.edu:8001/cgi-bin/browse/card.orf.cgi?query=%s"
			 , ac_tag_printable (gmp->gene, "Genefinder", "*")
			 ) ;
		gmpURL (blkp, gmp
			,linkBuf
			,". Movies are available on Piano's site"
			) ;	
	      }	    
	  }
	gGene = ac_tag_table (oRnai, "Gene", h) ;
	nn2 = 0 ;
	if (gGene->rows > 1)
	  {
	    vtxtPrintf (blkp, ". Warning: this double stranded RNA may also interfere with gene ") ;
	    for (jr=0; jr < gGene->rows ; jr++)
	      {
		oGene = ac_table_obj (gGene, jr, 0, h) ;
		if (!strcmp (ac_name (oGene), ac_name (gmp->gene)))
		  continue ;
		if (nn2++) vtxtPrintf (blkp, ", ") ;
		gmpObjLink (blkp, gmp, oGene, 0) ;
	      }
	  }
      }
  
  ac_free (h) ;
  return nn ;
} /* ficheNewGeneRnaiStatement */

/*************/

static int ficheNewGeneAlleleStatement (vTXT blkp, GMP *gmp)
{
  int ir, nn = 0 ;
  const char *ptr, *ptr1 ;
  AC_TABLE gAllele ;
  AC_HANDLE h = ac_new_handle () ;
  
  /* we must use our own buffer to avoid repeats 
     between allele and knoc_out_allele */
  vTXT bfr = vtxtCreate () ;
  if (gmp->markup) vtxtMarkup (bfr) ;

  if ((gAllele = ac_tag_table (gmp->gene, "Allele", h)))
    {
      for (ir=0 ; ir < gAllele->rows ; ir++)
	{
	  if ((ptr = ac_table_printable (gAllele, ir, 0, 0)))
	    {
	      if (!nn++)
		{
		  vtxtEmptyLine (blkp, 1) ;
		  vtxtBold (bfr, "Knock-out allele") ;
		  vtxtPrintf (bfr,", deletion obtained by the ") ;
		  gmpURL (bfr, gmp, "http://elegans.bcgsc.bc.ca/knockout.shtml", "Gene Knockout Consortium") ;
		}
	      else
		vtxtPrintf (bfr, " ") ;
	      vtxtPrintf (bfr, "%s", gtCleanQuotes (ptr)) ;
	    }
        }
      ac_free (gAllele) ;
    }
  
  if ((ptr1 = vtxtPtr (bfr)))
    {
      vtxtPrint (blkp, ptr1) ;
    }
  
  /* now we write directly in blkp */
  if ((gAllele = ac_tag_table (gmp->gene, "knock_out_allele", h)))
    {
      for (ir=0 ; ir < gAllele->rows ; ir++)
	{
	  /* Allele name */
	  ptr = ac_table_printable (gAllele, ir, 2, 0) ;
	  
	  if (ptr && ptr1 && strstr (ptr1, ptr))
	    continue ; /* we have a text for this allele */
	  
	  if (!ptr &&
	      !ac_table_printable (gAllele, ir, 0, 0) &&
	      !ac_table_printable (gAllele, ir, 1, 0) 
	      ) /* incomplete data */ 
	    continue ; 
	  
	  /* export something, so we need a title */
	  if (!nn++)
	    {
	      vtxtBreak (blkp) ;
	      vtxtBold (blkp, "Knock-out allele") ;
	    }
	  else
	    vtxtDot (blkp) ;
	  
	  
	  if (!ptr) /* incomplete data */ 
	    {	       
	      vtxtPrintf (blkp,", deletion seeked by the ") ;
	      gmpURL (blkp, gmp, "http://elegans.bcgsc.bc.ca/knockout.shtml", "Gene Knockout Consortium ") ;
	      continue ;
	    }
	  
	  vtxtPrintf (blkp,", deletion obtained by the ") ;
	  gmpURL (blkp, gmp, "http://elegans.bcgsc.bc.ca/knockout.shtml", "Gene Knockout Consortium ") ;
	  /* allele name */
	  vtxtPrint (blkp, gtCleanQuotes (ptr)) ;
	  
	  /* strain */
	  if ((ptr = ac_table_printable (gAllele, ir, 3, 0)))
	    {
	      vtxtPrintf (blkp, " (strain %s)", gtCleanQuotes (ptr)) ;
	    }
	  
	  /* laboratory */
	  if ((ptr = ac_table_printable (gAllele, ir, 1, 0)))
	    {
	      /* oklahoma */
	      if (!strcasecmp (ptr, "Oklahoma"))
		vtxtPrintf (blkp, " [R Barstead, Oklahoma MRF, USA]") ;
	      else if (!strcasecmp (ptr, "Vancouver"))
		vtxtPrintf (blkp, " [D Moerman, UBC Vancouver, Canada]") ;
	      else if (!strcasecmp (ptr, "Japan"))
		vtxtPrintf (blkp, " [S Mitani, Tokyo, Japan]") ;
	      else if (!strcasecmp (ptr, "Sanger"))
		vtxtPrintf (blkp, " [A Coulson, Sanger Centre, UK]") ;
	    }
        }
      ac_free (gAllele) ;
    }
  ac_free (h) ;
  
  vtxtDestroy (bfr) ;
  return nn ;
} /* ficheNewGeneAlleleStatement */

/*************/

static int ficheNewGeneStrainStatement (vTXT blkp, GMP *gmp)
{
  int ir, nn = 0 ;
  const char *ptr ;
  AC_TABLE gStrain ;
  AC_OBJ oStrain ;
  AC_HANDLE h = ac_new_handle () ;

  if ((gStrain = ac_tag_table (gmp->gene, "Strain", h)))
    {
      for (ir=0 ; ir < gStrain->rows && ir < 5 ; ir++)
	{
	  if ((oStrain = ac_table_obj (gStrain, ir, 0, h)))
	    {
	      if (!nn++)
		{
		  vtxtBreak (blkp) ; if (0) vtxtEmptyLine (blkp, 1) ;
		  vtxtBold (blkp, gStrain->rows < 2 ? "Selected strain" : "Selected strains") ;
		  if (gmp->Spc == WORM)
		    {
		      vtxtPrintf (blkp, " available from the ") ;
		      gmpURL (blkp, gmp, "http://biosci.umn.edu/CGC", "CGC") ;
		    }
		}
	      vtxtBreak (blkp) ;
	      vtxtPrint (blkp, ac_name (oStrain)) ;
	      if ((ptr = ac_tag_printable (oStrain, "Genotype", 0))) vtxtPrintf (blkp, " %s", gtCleanQuotes (ptr)) ;
	      if ((ptr = ac_tag_printable (oStrain, "Remark", 0))) vtxtPrintf (blkp, " %s", gtCleanQuotes (ptr)) ;
	    }
	}
    }
  vtxtBreak (blkp) ;
  ac_free (h) ;
  
  return nn ;
} /* ficheNewGeneStrainStatement */

/*************/

static int ficheNewGenePhenotypeParagraph (vTXT blkp, GMP *gmp)
{
  int nn = 0 ;
  char *ptr ;
  vTXT bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  if (gmp->style != 'x') /* Danielle 2007_03_16  redundant par rapport a la table */
    {
      ficheNewGeneLocus_DescriptionStatement (bfr, gmp) ; /* locus description */
      ficheNewGeneOmimStatement (bfr, gmp) ; /* OMIM */
    }
  ficheNewGeneLocus_PhenotypeStatement (bfr, gmp) ; /* locus phenotype */
  ficheNewGeneRNAiStatement (bfr, gmp) ; /* rnai */
  ficheNewGeneAlleleStatement (bfr, gmp) ; /* contains just KO alleles  */
  ficheNewGeneStrainStatement (bfr, gmp) ; /* selected strains */
  
  if ((ptr = vtxtPtr (bfr)))
    {
      if (gmp->Spc == HUMAN)
	gmpSection (blkp, gmp, "tg_Phenotype", "Phenotype") ;
      else
	gmpSection (blkp, gmp, "tg_Worm_Phenotype", "Phenotype") ;
      vtxtPrint (blkp, ptr) ; 
    }
  vtxtDestroy (bfr) ;

  return nn ;
} /* ficheNewGenePhenotypeParagraph */

/***************************************************************************************/
#ifdef JUNK
static int ficheNewGeneGoStatement (vTXT blkp, GMP *gmp)
{
  int ir,ir2, jr=0, jr2=0, jrgo = 0 ;
  const char *ptr ;
  DICT *dict = 0 ;
  AC_TABLE oTmp, oCogInfo;
  AC_OBJ oCog ;
  AC_HANDLE h = ac_new_handle () ;

  jr2 += jr ;
  jr = 0 ;
  dict = dictCreate (16) ;

  if (ac_has_tag (gmp->gene, "go_m_ace") ||
      ac_has_tag (gmp->gene, "go_b_ace") ||
      ac_has_tag (gmp->gene, "go_c_ace") 
      )
    {
      AC_TABLE gProd, gTag ;
      int nprod = 0, j2 = 0 ;
      char *cp1, *cp2, *prefix = gmp->Spc == HUMAN ? "From Entrez Proteome or GOA annotation, the" : "The" ;
      gProd = ac_tag_table (gmp->gene, "Product", h) ;
      nprod = gProd ? gProd->rows : 0 ;
      
      j2 = 0 ;
      gTag = ac_tag_table (gmp->gene, "go_m_ace", h) ;
      for (ir = 0 ; gTag && ir < gTag->rows ; ir++)
	{
	  ptr = ac_table_printable (gTag, ir, 0, 0) ;
	  if (!ptr || !dictAdd (dict, ptr, 0))
	    continue ;
	  cp1 = gtCleanUp(ptr) ;
	  cp2 = strstr (cp1, " activit") ;
	  if (cp2) *cp2 = 0 ;
	  if (!dictAdd (dict, cp1, 0))
	    continue ;

	  if (strstr(ptr,"unknown")) continue ;
	  if (!jr++)
	    {
	      vtxtDot (blkp) ;
	      vtxtPrintf (blkp, "%s %s would have"
			  , prefix
			  , nprod > 1 ? "products" : "product") ;
	    }
	  else
	    vtxtPrintf (blkp, ",") ;
	  vtxtPrintf (blkp, " %s", cp1) ;
	  if ((ptr = ac_table_printable (gTag, ir, 1, 0)) &&
	      !strstr (ptr, "NULL"))
	    vtxtPrintf (blkp, " [%s]", gtSquareCleanUp(ptr)) ;
	  j2++ ;
	}
      if (j2)
	vtxtPrintf (blkp, " activit%s", j2 > 1 ? "ies" : "y") ;
      
      
      j2 = 0 ;
      gTag = ac_tag_table (gmp->gene, "go_b_ace", h) ;
      for (ir = 0 ; gTag && ir < gTag->rows ; ir++)
	{
	  ptr = ac_table_printable (gTag, ir, 0, 0) ;
	  if (!ptr || !dictAdd (dict, ptr, 0))
	    continue ;
	  if (strstr(ptr,"unknown")) continue ;
	  if (!jr++)
	    {
	      vtxtDot (blkp) ;
	      vtxtPrintf (blkp,  "%s %s would be involved in", prefix
			  ,nprod > 1 ? "products" : "product") ;
	    }
	  else
	    vtxtPrintf (blkp,  ",%s", j2 ? "" : " would be involved in") ;
	  vtxtPrintf (blkp, " %s", gtCleanUp(ptr)) ;
	  if ((ptr = ac_table_printable (gTag, ir, 1, 0)) &&
	      !strstr (ptr, "NULL"))
	    vtxtPrintf (blkp, " [%s]", gtSquareCleanUp(ptr)) ;
	  j2++ ;
	}
      
      j2 = 0 ;
      gTag = ac_tag_table (gmp->gene, "go_c_ace", h) ;
      for (ir = 0 ; gTag && ir < gTag->rows ; ir++)
	{
	  ptr = ac_table_printable (gTag, ir, 0, 0) ;
	  if (!ptr || !dictAdd (dict, ptr, 0))
	    continue ;
	  if (strstr(ptr,"unknown")) continue ;
	  if (!jr++)
	    {
	      vtxtDot (blkp) ;
	      vtxtPrintf (blkp, "%s %s would localize in", prefix
			  ,nprod > 1 ? "products" : "product") ;
	    }
	  else
	    vtxtPrintf (blkp,  "%s", j2 ? "," : " and would localize in") ;
	  vtxtPrintf (blkp, " %s", gtCleanUp(ptr)) ;
	  if ((ptr = ac_table_printable (gTag, ir, 1, 0)) &&
	      !strstr (ptr, "NULL"))
	    vtxtPrintf (blkp, " [%s]", gtSquareCleanUp(ptr)) ;
	  j2++ ;
	}
    }
  
  jrgo = jr ;
  if (1 ||
      gmp->Spc != WORM)   /* in human case also export the PFAMs */
    { 
      jr2 += jr ;
      jr = 0 ;
    }
  
  if(!jr && /* in worm, do not duplicate */
     (ac_has_tag (gmp->gene, "go_m_pfam") ||
      ac_has_tag (gmp->gene, "go_b_pfam") ||
      ac_has_tag (gmp->gene, "go_c_pfam") ||
      ac_has_tag (gmp->gene, "go_c_psort") ||
      ac_has_tag (gmp->gene, "go_locuslink")
      )
     )
    {
      AC_TABLE gProd, gTag ;
      AC_OBJ oProd;
      int nprod = 0, j2 = 0 ;
      char *cp1, *cp2, *prefix = jrgo ? "In addition, from Pfam homology, the" : "From Pfam homology, the" ;
      
      gProd = ac_tag_table (gmp->gene, "Product", h) ;
      for (ir = 0 ; gProd && ir < gProd->rows; ir++)
	{
	  oProd = ac_table_obj (gProd, ir, 0, h) ;
	  if (ac_has_tag (oProd, "Pfam")) nprod++ ;
	  ac_free (oProd) ;
	}
      
      j2 = 0 ;
      gTag = ac_tag_table (gmp->gene, "go_m_pfam", h) ;
      for (ir = 0 ; gTag && ir < gTag->rows ; ir++)
	{
	  ptr = ac_table_printable (gTag, ir, 0, 0) ;
	  if (!ptr || !dictAdd (dict, ptr, 0))
	    continue ;
	  cp1 = gtCleanUp(ptr) ;
	  cp2 = strstr (cp1, " activit") ;
	  if (cp2) *cp2 = 0 ;
	  if (!dictAdd (dict, cp1, 0))
	    continue ;
	  if (strstr(ptr,"unknown")) continue ;

	  if (!jr++)
	    {
	      vtxtDot (blkp) ;
	      vtxtPrintf (blkp, "%s %s would have"
			  , prefix
			  , nprod > 1 ? "products" : "product") ;
	    }
	  else
	    vtxtPrintf (blkp, ",") ;
	  vtxtPrintf (blkp, " %s", cp1) ;

	  if ((ptr = ac_table_printable (gTag, ir, 1, 0)) &&
	      !strstr (ptr, "NULL"))
	    vtxtPrintf (blkp, " [%s]", gtSquareCleanUp(ptr)) ;
	  j2++ ;
	}
      if (j2)  vtxtPrintf (blkp, " activit%s", j2 > 1 ? "ies" : "y") ;
      
      
      j2 = 0 ;
      gTag = ac_tag_table (gmp->gene, "go_b_pfam", h) ;
      for (ir = 0 ; gTag && ir < gTag->rows ; ir++)
	{
	   ptr = ac_table_printable (gTag, ir, 0, 0) ;
	  if (!ptr || !dictAdd (dict, ptr, 0))
	    continue ;
	  if (strstr(ptr,"unknown")) continue ;

	  if (!jr++)
	    {
	      vtxtDot (blkp) ;
	      vtxtPrintf (blkp,  "%s %s would be involved in", prefix
			  ,nprod > 1 ? "products" : "product") ;
	    }
	  else
	    vtxtPrintf (blkp,  ",%s", j2 ? "" : " would be involved in") ;
	  vtxtPrintf (blkp, " %s", gtCleanUp(ptr)) ;
	  if ((ptr = ac_table_printable (gTag, ir, 1, 0)) &&
	      !strstr (ptr, "NULL"))
	    vtxtPrintf (blkp, " [%s]", gtSquareCleanUp(ptr)) ;
	  j2++ ;
	}
      
      j2 = 0 ;
      gTag = ac_tag_table (gmp->gene, "go_c_pfam", h) ;
      for (ir = 0 ; gTag && ir < gTag->rows ; ir++)
	{
	  ptr = ac_table_printable (gTag, ir, 0, 0) ;
	  if (!ptr || !dictAdd (dict, ptr, 0))
	    continue ;
	  if (strstr(ptr,"unknown")) continue ;

	  if (!jr++)
	    {
	      vtxtDot (blkp) ;
	      vtxtPrintf (blkp, "%s %s would localize in", prefix
			  ,nprod > 1 ? "products" : "product") ;
	    }
	  else
	    vtxtPrintf (blkp,  "%s", j2 ? "," : " and would localize in") ;
	  vtxtPrintf (blkp, " %s", gtCleanUp(ptr)) ;
	  if ((ptr = ac_table_printable (gTag, ir, 1, 0)) &&
	      !strstr (ptr, "NULL"))
	    vtxtPrintf (blkp, " [%s]", gtSquareCleanUp(ptr)) ;
	  j2++ ;
	}
      
      if (0)
	{ /* this is a duplication and not reliable for nuclear,
	   *  danielle wants it out of the gene summary as of may 26 2003 
	   */
	  j2 = 0 ;
	  gTag = ac_tag_table (gmp->gene, "go_c_psort", h) ;
	  for (ir = 0 ; gTag && ir < gTag->rows ; ir++)
	    {
	      ptr = ac_table_printable (gTag, ir, 0, 0) ;
	      if (!ptr || !dictAdd (dict, ptr, 0))
		continue ;
	      if (strstr(ptr,"unknown")) continue ;
	      
	      if (!jr++)
		{
		  vtxtDot (blkp) ;
		  vtxtPrintf (blkp, "According to PSORTII, it encodes %s which would localise in the"
			      ,nprod > 0 ? "a product" : "products") ;
		}
	      else
		vtxtPrintf (blkp,  ",%s", j2 ? "" : " and would according to PSORT localise in the") ;
	      vtxtPrintf (blkp, " %s", gtCleanUp(ptr)) ;
	      if ((ptr = ac_table_printable (gTag, ir, 1, 0)) &&
		  !strstr (ptr, "NULL"))
		vtxtPrintf (blkp, " [%s]", gtSquareCleanUp(ptr)) ;
	      j2++ ;
	    }
	}
      
      /* pfam go c ou psort title */
    } 
  
  dictDestroy (dict) ;

  /* COG */
  jr2 += jr ;
  jr = 0 ;
  if (0 &&
      ac_has_tag (gmp->gene, "Product") && (oTmp = ac_tag_table (gmp->gene, "COG", h)))
    {
      AC_TABLE gProd ;
      int nprod = 0, j2 = 0 ;
      
      gProd = ac_tag_table (gmp->gene, "Product", h) ;
      nprod = gProd->rows ;

      for (ir=0; oTmp && ir < oTmp->rows && oTmp->cols >= 1; ir++)
	{
	  oCog = ac_table_obj (oTmp, ir, 0, h) ;

	  if ((ptr = ac_tag_printable (oCog, "Title", 0)))
	    {
	      if (!strcmp (ptr, "Poorly_characterized"))
		continue ;
	      if (!jr++)
		{
		  vtxtDot (blkp) ;
		  vtxtPrintf (blkp, "From COG, the %s " 
			      , nprod > 0 ? "product has" : "products have"
			      ) ;
		}
	      else
		vtxtPrintf (blkp, ",") ;
	      vtxtPrintf (blkp, " %s"
			  , gtCleanUp(ptr)) ;
	      j2++ ;
	    }
	  if (j2)  vtxtPrintf (blkp, " activit%s", j2 > 1 ? "ies" : "y") ;
	  ac_free (oCog) ;
	}

      j2 = 0 ; 
      dict = dictCreate (16) ;
      for (ir=0;ir < oTmp->rows && oTmp->cols >= 1; ir++)
	{
	  const char *ptr2 ;
	  oCog = ac_table_obj (oTmp, ir, 0, h) ;
	  if ((oCogInfo = ac_tag_table (oCog, "Info", h)))
	    for (ir2=0;ir2 < oCogInfo->rows && oCogInfo->cols >= 2;ir2++)
	      if ((ptr = ac_table_printable (oCogInfo, ir2, 1, 0)))
		{
		  if (!dictAdd (dict, ptr, 0))
		    continue ;
		  ptr2 = ac_table_printable (oCogInfo, ir2, 0, 0) ;
		  if (!ptr2)
		    continue ;
		  if (!strcmp (ptr2, "Title"))
		    continue ;
		  if (!strcmp (ptr2, "Poorly_characterized"))
		    continue ;
		  if (!jr++)
		    vtxtPrintf (blkp, " From COG, %s involved in" 
				, nprod > 0 ? "the product is" : "the products are"
				) ;
		  else
		    vtxtPrintf (blkp,  ",%s", j2 ? "" : " involved in") ;
		  vtxtPrintf (blkp, " %s"
			      , gtCleanUp (ptr2)) ;
		  vtxtPrintf (blkp, ", more precisely %s"
			      , gtCleanUp(ptr)) ;
		  j2++ ;
		} 
	  ac_free (oCog) ;
	}
      dictDestroy (dict) ;
    }
  jr2 += jr ;
  vtxtBreak (blkp) ;
  ac_free (h) ;

  return jr2 ;
}  /* ficheNewGeneGoStatement */
#endif
/***************************************************************************************/

static BOOL ficheNewGeneDisease (vTXT blkp, GMP *gmp, Array bb)
{
  char *ptr ;
  const char *ccp ;
  int nn = 0 ;
  BOOL hasVote = FALSE ;  
  char linkBuf[2000] ;
  AC_HANDLE h = ac_new_handle () ;
  vTXT bfr = vtxtHandleCreate (h) ; 

  if (gmp->markup) vtxtMarkup (bfr) ;

  nn += ficheNewGeneDiseasePathwaysBioProcessTableNew (bfr, gmp, &hasVote, bb) ;
  if (nn)
    {
      ptr = vtxtPtr (bfr) ;
      if (gmp->Spc == MOUSE)
	{
	  AC_OBJ obj1 ;
	  gmpSection (blkp, gmp, "tg_disease", "Diseases or phenotypes (MGI)") ;
	  
	  vtxtPrint (blkp, "Look here for studies trying to associate a particular disease or phenotype to variations in sequence, expression or function of this gene") ;
	  
	  obj1 = ac_objquery_obj (gmp->gene, ">Extern MGI", h) ;
	  ccp = obj1 ? ac_name(obj1) : 0 ;
	  if (ccp && !strncmp (ccp, "MGI_", 4))
	    {
	      vtxtPrintf (blkp, ". The summary below gives hints, but may be misleading since we do not distinguish loss of function from gain of function phenotypes, and we do not report the specific genetic backgrounds (e.g. some phenotypes may be seen only in combination with other mutated genes). For more information about all aspects of the gene: phenotypes, alleles, strains, expression, etc, please see the ") ;
	      
	      sprintf (linkBuf, MGI_LINK, ccp + 4) ;
	      gmpURL (blkp, gmp, linkBuf, "Jackson Laboratory Mouse Genome Database/Informatics") ;
	      vtxtPrintf (blkp, " site") ;
	    }
	}
      else if (gmp->Spc == HUMAN)
	{
	  gmpSection (blkp, gmp, "tg_disease", "Diseases (Mesh)") ;
	  
	  vtxtPrint (blkp, "Look here for studies trying to associate a particular disease to variations in sequence, expression or function of this gene") ;
	  if (hasVote)
	    vtxtPrintf (blkp, ". Warning: read the articles because our PubMed mining only gives an indication since we cannot automatically distinguish a yes from a no answer. Also please help us by <a href='mailto:mieg@ncbi.nlm.nih.gov?subject=vote_gene%s'><font color='red'><i>voting</i></font></a>: send your documented view on the validity of the association. Authored users comments will appear in future versions and impact our lists", ac_name(gmp->gene)) ;
	}
      vtxtBreak (blkp) ;
      vtxtPrint (blkp, ptr) ; 
    }

  ac_free (h) ;
  return TRUE ;
} /* ficheNewGeneDisease */

/***************************************************************************************/
/***************************************************************************************/

static void ficheNewGeneAceKogSubSection (vTXT blkp, GMP *gmp)
{
  char *ptr ;
  vTXT bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  ficheNewAceKogStatement (bfr, gmp, TRUE) ;
  
  if ((ptr = vtxtPtr (bfr)))
    {
      gmpSubSection (blkp, gmp, "tg_acekog", "Closest AceView homologs in other species") ;
      vtxtPrint (blkp, ptr) ; 
    }
  vtxtDestroy (bfr) ;
} /* ficheNewAceKogSubSection */

/***************************************************************************************/
/***************************************************************************************/

static void ficheNewGeneTaxblastParagraph (vTXT blkp, GMP *gmp)
{
  char *ptr ;
  vTXT bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  fichePRODUCTTaxTreeParagraphContent (bfr, gmp) ; 
  if ((ptr = vtxtPtr (bfr)))
    {
      gmpSection (blkp, gmp, "tg_taxblast", "Closest GenBank protein accessions in other organisms") ;
      vtxtPrint (blkp, ptr) ; 
    }
  vtxtDestroy (bfr) ;
} /* ficheNewGeneTaxblastParagraph */

/***************************************************************************************/
/***************************************************************************************/

static void ficheNewGeneInteractionsTableTable (AC_TABLE tbl1, GMP *gmp, Array bb, BOOL with_protein)
{
  AC_HANDLE h = ac_new_handle () ;
  Array aa = 0, aaU = arrayHandleCreate (12, ITRC, h) ;
  ITRC *itrc ;
  int nItrc1 = 0, nItrc2 ;
  const char *ccr ;
  char *cp, *cq, papBuf[32400] ;
  KEY k0, k1, k2, kk ;
  int aha, ir, ir1, ir2, nn, row = 0, nGwithPap = 0, nnGwithPap = 0 ;
  AC_OBJ obj, gene2 ;
  AC_TABLE tbl, tbl2 ;
  AC_KEYSET ks2 ;
  vTXT bfr = vtxtHandleCreate (h) ;  
  vTXT bfr1 = vtxtHandleCreate (h) ; 
  vTXT bfrPap = vtxtHandleCreate (h) ;  
  enum { COL_MESH=0, COL_EVIDENCE, COL_SOURCE, COL_GENES, COL_VOTE, COL_LAST} ; 

  if (gmp->markup) vtxtMarkup (bfr) ;   
  if (gmp->markup) vtxtMarkup (bfr1) ;   
  if (gmp->markup) vtxtMarkup (bfrPap) ;

  vtxtPrint (bfrPap, PUBMED_LINK) ; /* no %s included */ 
  /* Interactions */
  /*  AC_KEYSET geneidKs = ac_objquery_keyset (gmp->gene, ">geneid", h) ; */
  k0 = ac_obj_key (gmp->gene) ;

  for (aha = 0 ; aha < 1 ; aha++) /* formelly used to loop on with_gene with_protein */
    {
      tbl = ac_tag_table (gmp->gene, with_protein ? "with_protein" : "with_gene", h) ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  obj = ac_table_obj (tbl, ir, 0, h) ;
	  ccr = ac_tag_printable (obj, "Title", ac_name(obj)) ;
	  if (!strncasecmp (ccr, "gene ", 5)) ccr += 5 ;
	  row = tbl1->rows ;
	  ac_table_insert_text (tbl1, row, COL_LAST, with_protein ? "1" : "2") ;
	  vtxtClear (bfr) ;
	  gmpObjLink (bfr, gmp, obj, ccr) ;
	  if (gmp->Spc == WORM) 
	    {
	      int nnE = 0 ;
	      kk = ac_table_key (tbl, ir, 1, 0) ;
	      for (ir1 = ir ; ir1 < tbl->rows ; ir1++)
		if (kk == ac_table_key (tbl, ir1, 1, 0))
		  {
		    if (!nnE++)
		      vtxtPrint (bfr, "<span class='explain2'>") ;
		    vtxtPrint (bfr, ac_table_printable (tbl, ir1, 1, "")) ;
		  }
	      if (nnE)	
		vtxtPrint (bfr, "</span>") ;
	      }
	  ac_table_insert_text (tbl1, row, COL_MESH, vtxtPtr (bfr)) ;
      
	  vtxtClear (bfr) ;
	  vtxtClear (bfr1) ;
	  vtxtPrint (bfr1, PUBMED_LINK) ; /* no %s included */ 

	  ks2 = ac_objquery_keyset (obj, with_protein ? ">with_protein" : "{>with_gene}SETOR{>with_protein}", h) ;
	  if (ks2 && ac_keyset_count (ks2))
	    {
	      double cc = Ninteract2 ; /* do not contribute to friends */
	      if (!aa) aa = arrayHandleCreate (64, ITRC, h) ;
	      tbl2 = ac_keyset_table (ks2, 0, -1, 0, h) ;
	      vtxtClear (bfr) ;
	      if (tbl2 && tbl2->rows)
		{ 
		  BOOL show = TRUE ;
		  vtxtPrintf (bfr, "%s gene%s", isOne (tbl2->rows) , _multi (tbl2->rows)) ;
		  for (ir2 = 0 ; tbl2 && ir2 < tbl2->rows ; ir2++)
		    {
		      k2 = ac_table_key (tbl2, ir2, 0, 0) ;
		      gene2 = ac_table_obj (tbl2, ir2, 0, h) ; 
		      for (nItrc2 = 0 ; nItrc2 < nItrc1 ; nItrc2++)
			{
			  itrc = arrp (aa, nItrc2, ITRC) ;
			  if (itrc->gene == k2)
			    { itrc->nn++ ; itrc->cc += cc ; itrc->type = with_protein ? IPtype : IGtype ; break ; }
			}
		      if (k2 != k0 && nItrc2 >= nItrc1)
			{
			  itrc = arrayp (aa, nItrc1++, ITRC) ;
			  itrc->gene = k2 ;
			  itrc->Gene = gene2 ;
			  itrc->nn = 1 ;
			  itrc->cc = cc ;
			  itrc->type = with_protein ? IPtype : IGtype ;
			}
		      if (show)
			{
			  vtxtPrint (bfr, ir ? ", " : ": ") ;
			  gmpObjLink (bfr, gmp, gene2, ac_name (gene2)) ;
			  
			  if (ir2 > 4 &&  tbl2->rows > 8)
			    { vtxtPrint (bfr, "...") ; show = FALSE ; }
			}
		    }
		}
	      vtxtDot (bfr) ;
	      ac_table_insert_text (tbl1, row, COL_GENES, vtxtPtr (bfr)) ;
	    }

	  k1 = ac_table_key (tbl, ir, 0, 0) ;
	  if (1)
	    {
	      int iU, nItrc11 = arrayMax (aaU), nItrc21 ;
	      double cc = 0 ;

	      iU = ac_keyset_count (ks2) ;
	      if (iU < 2) iU = 2 ;
	      iU = 1 ; /* Force since we do not care about 
			* eventual second interactors to weigth direct ones 
			*/
	      if (iU < 10000) cc = - log(((double)iU)/10000.0) ;
	      cc *= with_protein ? NPinteract1 : NGinteract1 ;
	      
	      for (nItrc21 = 0 ; nItrc21 < nItrc11 ; nItrc21++)
		{
		  itrc = arrp (aaU, nItrc21, ITRC) ;
		  if (itrc->gene == k1)
		    { itrc->nn++ ; itrc->cc += cc ; itrc->type = with_protein ? IPtype : IGtype ; break ; }
		}
	      if (k1 != k0 && nItrc21 >= nItrc11)
		{
		  itrc = arrayp (aaU, nItrc11++, ITRC) ;
		  itrc->gene = k1 ;
		  itrc->Gene = 0 ;
		  itrc->nn = 1 ;
		  itrc->cc = cc ;
		  itrc->type |= with_protein ? IPtype : IGtype ;
		}
	    }
	  vtxtClear (bfr) ;
	  if (gmp->Spc != WORM && ! with_protein) /* the with_protein schema does not have a ?paper in second column */
	    {
	      for (ir1 = ir, ir2 = 0 ; tbl && ir1 < tbl->rows ; ir1++)
		{
		  kk = ac_table_key (tbl, ir1, 0, 0) ;
		  if (kk != k1) break ;
		  if (!ac_table_type (tbl, ir1, 1))
		    continue ;
		  if (strncasecmp (ac_table_printable (tbl, ir1, 1 , "xx"), "pm", 2))
		    continue ;
		  strncpy (papBuf, ac_table_printable (tbl, ir1, 1, "pm1") + 2, 32000) ;
		  strcat (papBuf, ",") ;
		  if (! strstr (vtxtPtr (bfr1), papBuf))
		    { ir2++ ; vtxtPrint (bfr1, papBuf) ; }
		  if (!vtxtPtr (bfrPap) || ! strstr (vtxtPtr (bfrPap), papBuf))
		    {
		      nnGwithPap++ ;
		      vtxtPrint (bfrPap,papBuf) ;
		    }
		}
	      if (ir2)
		{ 
		  nGwithPap++ ; 
		  gmpURL (bfr, gmp, vtxtPtr (bfr1)
			  , messprintf ("%d article%s", ir2, _multi(ir2))) ;
		  ac_table_insert_text (tbl1, row, COL_EVIDENCE, vtxtPtr (bfr)) ;
		}
	    }
	 
	  vtxtClear (bfr) ;
	  if (gmp->Spc != WORM && ! with_protein)
	    {
	      k1 = ac_table_key (tbl, ir, 0, 0) ;
	      for (ir1 = ir ; ir1 < tbl->rows ; ir1++)
		{
		  kk = ac_table_key (tbl, ir1, 0, 0) ;
		  if (kk != k1) break ;
		  if (!ac_table_type (tbl, ir1, 2))
		    continue ;
		  cp = strnew (ac_table_printable (tbl, ir1, 2, ""), h) ;
		  cq = strstr (cp, ":") ;
		  if (cq) *cq = 0 ;
		  if (vtxtPtr (bfr) && strstr(vtxtPtr (bfr), cp))
		    continue ;
		  vtxtBreak (bfr) ;
		  vtxtPrint (bfr, cp) ;
		}
	      ac_table_insert_text (tbl1, row, COL_SOURCE, vtxtPtr (bfr)) ;
	    }

	  for (ir1 = ir + 1 ; ir1 < tbl->rows ; ir1++)
	    {
	      kk = ac_table_key (tbl, ir1, 0, 0) ;
	      if (kk != k1)
		break ;
	    }
	  ir = ir1 - 1 ; /* jump to the end of this interactor */
	}
    }

  if (nGwithPap > 1)
    {
      row = tbl1->rows ;
      ac_table_insert_text (tbl1, row, COL_MESH, "Cumulated literature") ;
      ir = nnGwithPap ;
      vtxtClear (bfr) ; 
      gmpURL (bfr, gmp, vtxtPtr (bfrPap)
	      , messprintf ("%d article%s", ir, _multi(ir))) ;
      ac_table_insert_text (tbl1, row, COL_EVIDENCE, vtxtPtr (bfr)) ;
    }

  if (bb && aaU && arrayMax (aaU))
    { /* register the first interactors */
      int ir2 ;
      ITRC *itrc2 ;

      arraySort (aaU, itrcOrder) ;
      if (bb)
	for (nItrc2 = 0, ir2 = arrayMax (bb) ; nItrc2 < arrayMax (aaU) ; ir2++, nItrc2++)
	  {
	    itrc = arrp (aaU, nItrc2, ITRC) ;
	    itrc2 = arrayp (bb, ir2, ITRC) ;
	    *itrc2 = *itrc ;
	    itrc2->Gene = 0 ;
	  }
    }
  
  if (aa && arrayMax (aa) && aaU && arrayMax (aaU))
    { /* we want the intersect of the first and the second interactors */
      int i ;
      KEYSET ks1, ks2, ks3 ;
      AC_OBJ Gene = 0 ;

      ks1 = keySetCreate () ;
      ks2 = keySetCreate () ;

      for (i = 0, itrc = arrp (aaU, 0, ITRC) ; i < arrayMax (aaU) ; itrc++, i++)
	keySet (ks1, i) = itrc->gene ;
      for (i = 0, itrc = arrp (aa, 0, ITRC) ; i < arrayMax (aa) ; itrc++, i++)
	keySet (ks2, i) = itrc->gene ;
      keySetSort (ks1) ; keySetCompress (ks1) ;  
      keySetSort (ks2) ; keySetCompress (ks2) ;  
      ks3 = keySetAND (ks1, ks2) ;

      if (keySetMax (ks3))
	{
	  row = tbl1->rows ;
	  ac_table_insert_text (tbl1, row, COL_MESH, "<font color='red'>Putative partners in a complex</font>") ;
	  
	  vtxtClear (bfr) ; 
	  for (i = 0 ; i < keySetMax (ks3) ; i++)
	    {
	      k1 = keySet (ks3, i) ;
	      vtxtPrint (bfr, i ? ", " : "") ;
	      Gene = ac_get_obj (gmp->db, "gene", name(k1), 0) ;
	      gmpObjLink (bfr, gmp, Gene, ac_name (Gene)) ;
	      ac_free (Gene) ;
	      
	      if (i > 10 &&  keySetMax (ks3) > 20)
		{ vtxtPrint (bfr, "...") ; break ; }
	    }
	  ac_table_insert_text (tbl1, row, COL_GENES, vtxtPtr (bfr)) ;
	}
      keySetDestroy (ks1) ;
      keySetDestroy (ks2) ;
      keySetDestroy (ks3) ;
    }
  if (aa && arrayMax (aa))
    {
      int ir2 ;
      double cc, oldcc ;

      arraySort (aa, itrcOrder) ;
      itrc = arrp (aa, 0, ITRC) ;
      nn = itrc->nn ;
      if (nn > 1)
	{
	  row = tbl1->rows ;
	  ac_table_insert_text (tbl1, row, COL_MESH, "<font color='red'>Most connected second interactors</font>") ;
	  
	  cc = oldcc = itrc->cc ;
	  for (nItrc2 = ir = 0 ; nItrc2 < nItrc1 ; nItrc2++)
	    {
	      itrc = arrp (aa, nItrc2, ITRC) ;
	      if (2*itrc->cc >= cc) ir++ ;
	    }
	  
	  vtxtClear (bfr) ; 
	  for (nItrc2 = ir2 = 0 ; nItrc2 < nItrc1 ; nItrc2++)
	    {
	      itrc = arrp (aa, nItrc2, ITRC) ;
	      if (2*itrc->cc >= cc || FDEBUG)
		{
		  if (itrc->cc < oldcc)
		    vtxtPrint (bfr, ir2++ ? " | " : "") ;
		  else
		    vtxtPrint (bfr, ir2++ ? ", " : "") ;
		  oldcc = itrc->cc ;
		  gmpObjLink (bfr, gmp, itrc->Gene, ac_name (itrc->Gene)) ;
	      
		  if (ir2 > 8 &&  ir > 12 && ! FDEBUG)
		    { vtxtPrint (bfr, "...") ; break ; }
		}
	    }
	  ac_table_insert_text (tbl1, row, COL_GENES, vtxtPtr (bfr)) ;
	}
    }
  ac_free (h) ;
  return ;
} /* ficheNewGeneInteractionsTableTable */

/***************************************************************************************/

static int ficheNewGeneInteractionsTableNew (vTXT blkp, GMP *gmp, Array bb, BOOL with_protein)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl ; 
  vTXT bfrPap = vtxtHandleCreate (h) ; 
  int row = 0, cols[6] = {1,2,3,4,0} ;
  enum { COL_MESH=0, COL_EVIDENCE, COL_SOURCE, COL_GENES, COL_VOTE, COL_LAST} ; 
  const char *colNames[]={
			   "Interacting gene or protein", 
			   "Evidence",
			   "Source",
			   "Other interactors of this interactor",
			   0} ; 
  if (with_protein)
    colNames[0] = hprintf (h, "Protein interacting with %s", ac_name(gmp->gene)) ; 
  else
    colNames[0] = hprintf (h, "Gene interacting with %s", ac_name(gmp->gene)) ; 
  if (gmp->Spc != WORM)
    colNames[3] = hprintf (h, "Other interactors of this interactor<br>(with published support)") ;
  /* contruct the table */
  tbl = ac_empty_table (80, 6, h) ; /* tbl->rows is a hint, not a hard limit */

  /* format the table and add the http links */
  vtxtClear (bfrPap) ;
  vtxtPrint (bfrPap, PUBMED_LINK) ; /* no %s included */ 
  ficheNewGeneInteractionsTableTable (tbl, gmp, bb, with_protein) ;
  /* export */
  if ((row = tbl->rows))
    ac_table_display (blkp 
		      , tbl, colNames
		      , cols, 1
		      , 0, 0, 0
		      , 0
		      ) ;

  ac_free (h) ; 

  return row ;
} /* ficheNewGeneInteractionsTableNew */

/**************/

static void ficheNewGeneInteractionsTableCaption (vTXT blkp, GMP *gmp)
{
  const char *tmpTxt ;
  return ;
  gmpCaption (blkp, gmp, "tg_interaction_caption", "Legend") ;
  vtxtPrint (blkp, 
	     "We simply display here the information collected in Entrez"
	     " Gene, but merge the various sources and the supporting articles to get a"
	     " non-redundant list"
	     ) ;
     
  tmpTxt = ac_tag_printable (gmp->gene, "GeneId", 0) ; 
  if (tmpTxt)
    {
      char	linkBuf[vONLINEMAXURLSIZE] ;
      vtxtPrint (blkp, ". More details and links to the sources are available in ") ;
      sprintf (linkBuf, ENTREZ_LINK, tmpTxt) ; 
      gmpURL (blkp, gmp, linkBuf, "Entrez") ;
    }

  vtxtPrint (blkp, 
	     ". <br>"
	     " Column 1: the interacting gene, links to the gene in AceView<br>"
	     " Column 2: evidence, links to PubMed<br>"
	     " Column 3: source of the annotation.<br>"
	     " Column 4: list of genes or proteins annotated as interacting with the"
	     " gene in column 1. To obtain the complete list of its interactors"
	     ", follow the link in column 1 and look at the 'interactions' paragraph"
	     ) ; 

} /* ficheNewGeneInteractionsTableCaption */

/**************/

static void ficheNewGeneInteractionParagraph (vTXT blkp, GMP *gmp, Array bb)
{
  char *ptr ;
  vTXT bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  ficheNewGeneInteractionsTableNew (bfr, gmp, bb, TRUE) ;
  ficheNewGeneInteractionsTableNew (bfr, gmp, bb, FALSE) ;

  if ((ptr = vtxtPtr (bfr)))
    {
      gmpSection (blkp, gmp, "tg_protein_interaction", "Interactions") ;
      vtxtPrint (blkp, ptr) ; 
      ficheNewGeneInteractionsTableCaption (blkp, gmp) ;
    }
  vtxtDestroy (bfr) ;
} /* ficheNewGeneInteractionParagraph */

/***************************************************************************************/
/***************************************************************************************/

static BOOL ficheNewGeneExpressionTissueIsBaddy (const char *ccp)
{
  const char **baddy ;
  const char *baddies[] =
    { 
      "syntheti*"
      , "partial"
      , "pooled"
      , "pcr*"
      , "mixed"
      , "mixed (*"
      , "mixed tissues"
      /* *cell?line*" */
      , 0
    } ;
  for (baddy = baddies ; *baddy ; baddy++) 
    if (pickMatch (ccp, *baddy))
      return TRUE ;
  return FALSE ;
} /* ficheNewGeneExpressionTissueIsBaddy */

/***************************************************************************************/
/**************************************************************************************/
typedef struct tissueStruct { KEY tissue ; int nn ; } TISSUE ;

static int tissueOrder (const void *va, const void *vb)
{
  const  TISSUE *up = (const TISSUE*) va ;
  const  TISSUE *vp = (const TISSUE*) vb ;
  
  if (up->tissue < vp->tissue)
    return -1 ;
  if (up->tissue > vp->tissue)
    return 1 ;
  return 0 ;
} /* tissueOrder */

static int tissueNnOrder (const void *va, const void *vb)
{
  const  TISSUE *up = (const TISSUE*) va ;
  const  TISSUE *vp = (const TISSUE*) vb ;
  int nn ;
  
  nn = up->nn - vp->nn ;
  if (nn) return - nn ; /* large nn first */
  return lexstrcmp (name(up->tissue), name(vp->tissue)) ;
} /* tissueNnOrder */

static void tissueCompress (Array tissues)
{
  TISSUE *up, *vp ;
  int i, j ;
  
  for (i = 0 ; i < arrayMax (tissues) ; i++)
    {
      up = arrp (tissues, i, TISSUE) ;
      for (j = i+1, vp = up+1 ; up->nn && j < arrayMax (tissues) ; vp++, j++)
	if (vp->tissue == up->tissue)
	  { vp->nn = 0 ; up->nn++ ; }
    }
  return ;
} /* tissueCompress */


char *ficheNewGeneExpressionTissue (vTXT blkp, GMP *gmp, AC_KEYSET clones, int nMax1, int nMax2, char showOther, int nClo)
{
  AC_HANDLE h = ac_new_handle () ;
  Array tissues = 0 ;
  TISSUE *tt ;
  int i, ir1, jr = 0, n1 = 0 ;
  AC_ITER iter = ac_keyset_iter (clones,1, h) ;
  AC_TABLE tbl1 ;
  AC_OBJ clone = 0, est = 0, lib = 0 ;
  KEY stage ;

  
  if (0 && gmp && gmp->Spc == WORM)
    return 0 ; 
  tissues = arrayHandleCreate (60, TISSUE, h) ;
  while (ac_free(clone), (clone = ac_iter_obj (iter)))
    {
      lib = ac_tag_obj (clone, "Library", h) ;
      stage = ac_tag_key (lib, "Stage", 0) ;
      if (stage)
	{
	  tt = arrayp (tissues, jr++, TISSUE) ;
	  tt->tissue = stage ;
	  tt->nn = 1 ;
	  continue ;
	}
      est = ac_tag_obj (clone, "Read", h) ;
      tbl1 = ac_tag_table (est, "Tissue", h) ;
      for (ir1 = 0 ; tbl1 && ir1 < tbl1->rows ; ir1++)
	{
	  if (!ficheNewGeneExpressionTissueIsBaddy (ac_table_printable (tbl1, ir1, 0, "")))
	    {
	      tt = arrayp (tissues, jr++, TISSUE) ;
	      tt->tissue = ac_table_key (tbl1, ir1, 0, 0) ;
	      tt->nn = 1 ;
	    }
	}
    } 
  ac_free(clone) ;
  if (arrayMax (tissues))
    {
      arraySort (tissues, tissueOrder) ;
      tissueCompress (tissues) ;
      arraySort (tissues, tissueNnOrder) ;
      
      for (i = n1 = 0 ; i < arrayMax (tissues) ; i++)
	{
	  tt = arrp (tissues, i, TISSUE) ;
	  if (!tt->nn)
	    { 
	      arrayMax (tissues) = i + 1 ;
	      break ;
	    }
	}
      for (i = n1 = 0 ; i < arrayMax (tissues) ; i++)
	{
	  tt = arrp (tissues, i, TISSUE) ;
	  if (tt->nn)
	    {
	      if (nClo > 1)
		vtxtPrintf (blkp 
			    , "%s%s (%s%s%s)"
			    , i ? ", " : ", some from "
			    , gtLowerCleanUp (name (tt->tissue))
			    , i ? "" : "seen "
			    , tt->nn > 1 ? isOne(tt->nn) : "once"
			    , i ? "" : (tt->nn > 1 ? " times" :"")
			    );
	      else if (nClo == -1)
		vtxtPrintf (blkp 
			    , "%s%s (%d)"
			    , i ? ", " : ""
			    , gtLowerCleanUp (name (tt->tissue))
			    , tt->nn
			    );
	      else
		vtxtPrintf (blkp , " %s%s", nClo ? " from " : "", gtLowerCleanUp (name (tt->tissue))) ;
	    }
	  else
	    break ;
	  if (arrayMax (tissues) > nMax2 && i > nMax1)
	    { n1 = arrayMax (tissues) - i ; break ; }
	}  
    }
  if (n1 > 0 && showOther)
    {
       if (showOther == 'm')
	 vtxtPrintf (blkp, " and %d other <a href=\"javascript:openAnchor (0,'mRNA_expression')\">tissue%s</a>", n1, _multi(n1)) ;
       else if (showOther == 'g')
	 vtxtPrintf (blkp, " and %d other <a href=\"javascript:openAnchor ('fexp', 'tg_expression')\">tissue%s</a>", n1, _multi(n1)) ;
       else
	 vtxtPrintf (blkp, " and %d other tissue%s", n1, _multi(n1)) ;
    }
  ac_free (h) ;
  return vtxtPtr (blkp) ;
} /* ficheNewGeneExpressionTissue */

/***************************************************************************************/

static int ficheNewGeneCountAnomalousClones (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_KEYSET anno = gmp->tg ? ac_objquery_keyset (gmp->tg, ">cdna_clone ; anomalous_clone", h) : 0 ;
  int nn = anno ? ac_keyset_count (anno) : 0 ;

  if (nn)
    {
      vtxtDot (blkp) ;
      if (nn > 1)
	vtxtPrintf (blkp, " We annotate  <a href=\"javascript:openAnchor ('fexp', 'tg_expression')\">structural defects or features</a> in %d cDNA clones", nn) ;
      else
	vtxtPrintf (blkp, " We annotate  <a href=\"javascript:openAnchor ('fexp', 'tg_expression')\">structural defects or features</a> in one cDNA clone") ;
    }

  ac_free (h) ;
  return nn ;
} /* ficheNewGeneCountAnomalousClones */

/***************************************************************************************/

static int ficheNewGeneExpressionTissueStatement (vTXT blkp, GMP *gmp, BOOL decorate)
{	
  const char *ptr ;
  int  level = 0, ok = 0, numReads = 0, numBuriedReads = 0, numClones = 0, numClones2 = 0 ;
  int numMrnas = 0 ; 
  AC_HANDLE h = ac_new_handle () ;
  
  /*
  if (gmp->gene)
    {
      if ((numShed = ac_keyset_count (ac_objquery_keyset (gmp->gene, ">transcribed_gene  shedded_from", h))))
	numReadsShed = ac_keyset_count (ac_objquery_keyset (gmp->gene, ">transcribed_gene  shedded_from ; > read", h)) ;
    }
  */

  vtxtBreak (blkp) ;
  if (gmp->tg)
    {
      level = 4 ;

      numMrnas = ac_keyset_count (ac_objquery_keyset (gmp->tg, "> mrna ; gt_ag || gc_ag", h)) ;
      if (gmp->Spc == WORM && gmp->gene &&
	  ac_has_tag (gmp->gene, "Has_cDNA_clone")) /* priorite au canning */
	{
	  numClones = ac_keyset_count (ac_objquery_keyset (gmp->gene, "Follow Has_cDNA_clone", h)) ;
	  numReads = 0 ;
	}
      else if (gmp->gene)
	{
	  AC_KEYSET ksr1 = ac_objquery_keyset (gmp->gene, ">transcribed_gene; >Read ; NOT Composite && NOT IS U*", h) ;
	  AC_KEYSET ksr2 = ac_ksquery_keyset (ksr1, ">buries;>buried_est;NOT Composite &&  NOT Composite && NOT IS U*", h) ;
	  /*
	    if we had access to the buried_reads->cdna_clones we would do:
	    AC_KEYSET ksc1 = ac_objquery_keyset (gmp->gene, ">transcribed_gene; >cdna_clone", h) ;
	    AC_KEYSET ksc2 = ksr2 ? ac_ksquery_keyset (ksr2, ">cdna_clone", h) : 0 ;
	    numClones = ac_keyset_count (ksc1) + (ksc2 ?  ac_keyset_count (ksc2) : 0) ;
	    numReads = ac_keyset_count (ksr1) + (ksr2 ?  ac_keyset_count (ksr2) : 0) ;
	  */
	  AC_KEYSET ksc1 = ac_objquery_keyset (gmp->gene, ">transcribed_gene; >cdna_clone; NOT X?_* && NOT IS U*", h) ; 
	  AC_KEYSET ksc2 = ac_objquery_keyset (gmp->gene, ">transcribed_gene; >cdna_clone; IS X?_* || IS U*", h) ; 
	  numReads = ac_keyset_count (ksr1) ;
	  numBuriedReads = ac_keyset_count (ksr2) ;
	  numClones = ac_keyset_count (ksc1) ;
	  numClones2 = ac_keyset_count (ksc2) ;
 	}

      if (numClones > 0)
	{
	  ok++ ;
	  if (gmp->markup && ! ac_has_tag (gmp->gene, "Cloud_gene") && numMrnas > 1)
	    vtxtPrint (blkp, "<b>Expression:</b> ") ;
	  ficheGeneExpressionLevelStatement (blkp, gmp, 0) ;
	  if (1)
	    {
	      vtxtDot (blkp) ;
	      vtxtPrint (blkp, "The ") ;
	      if (gmp->markup)
		{
		  vtxtPrintf (blkp, "<a href=\"javascript:openAnchor ('fgene', 'Gene_sequences')\">sequence</a>", ac_class (gmp->gene), ac_name (gmp->gene)) ;      
		} 
	      else
		vtxtPrint (blkp, "sequence") ;
	      vtxtPrint (blkp, " of this gene") ;
	    }
	  else
	    vtxtPrint (blkp, "It") ;
	  if (gmp->Spc == WORM)
	    {
	      if (numClones + numClones2 > 1)
		{
		  if (gmp->markup && gmp->view != 'c')
		    vtxtPrintf (blkp, " is defined by <a href=\"javascript:openAnchor"
				"('fexp', 'tg_expression')\">%d cDNA clones</a>"
				, numClones) ;
		  else
		    vtxtPrintf (blkp, " is defined by %d cDNA clones", numClones) ;
		  if (numClones2)
		    vtxtPrintf (blkp, " %s %s element%s defined by RNA-seq" 
				, numClones > 0 ? "and" :  "is defined by"
				, _isone(numClones2)
				, _multi(numClones2)
				) ;
		}
	      else
		{ 
		  AC_OBJ oClone = ac_tag_obj (gmp->tg, "cDNA_clone", h) ;
		  AC_OBJ oLib = oClone ? ac_tag_obj (oClone, "Library", h) : 0 ;

		  vtxtPrintf (blkp, " is defined by one cDNA clone") ;
		  if (oLib && (ptr = ac_tag_printable (oLib, "Title", 0)))
		    vtxtPrintf (blkp, 
				" from the %s"
				, ptr) ;
		}
	    }
	  else
	    {
	      char buff1[256] ;
		  
	      vtxtPrintf (blkp, " is defined by") ;
	      if (numReads > 1)
		{
		  sprintf (buff1, " %d GenBank accessions", numReads + numBuriedReads) ;
		  
		  if (gmp->markup)
		    vtxtPrintf (blkp, " <a href=\"javascript:openAnchor ('fexp', 'tg_expression')\">%s</a>"
				, buff1) ; 
		  else 
		    vtxtPrint (blkp, buff1) ;
		    
		  if (numBuriedReads == 0 && numReads > numClones)
		    vtxtPrintf (blkp, " from %s cDNA clone%s", _isone(numClones), _multi(numClones)) ;
		}
	      else 
		{ 
		  AC_OBJ oClone = ac_tag_obj (gmp->tg, "cDNA_clone", h) ;
		  AC_TABLE gRead = ac_tag_table (gmp->tg, "Read", h) ;
	
		  vtxtPrintf (blkp, " one GenBank accession") ;
		  vtxtPrintf (blkp, numReads > 1 ? " from one cDNA clone " : " ") ;

		  if (gRead && gRead->rows >= 1 && gmp->Spc != WORM)
		    sprintf (buff1, "%s", ac_table_printable (gRead, 0, 0, "")) ;
		  else
		    sprintf (buff1, "%s", ac_name(oClone)) ; 
		  
		  if (gmp->markup)
		    vtxtPrintf (blkp, " <a href=\"javascript:openAceViewLink ('%s', '%s', 'clones')\">%s</a>"
				, ac_class (gmp->gene), ac_name (gmp->gene), buff1) ; 
		  else 
		    vtxtPrint (blkp, buff1) ;
		}
	    }
	}
    }
  else if (gmp->pg && ac_has_tag (gmp->pg, "CDS"))
    {
      switch ((level = gtPredictedMrnaSupport (gmp->pg, &numClones)))
	{
	case 4:
	  if (numClones)
	    {
	      if (decorate)
		{
		  vtxtPrint (blkp, "The ") ;
		  if (gmp->markup)
		    {
		      vtxtPrintf (blkp, "<a href=\"javascript:openAnchor (0, 'Gene_sequences')\">sequence</a>"
				  , ac_class (gmp->gene), ac_name (gmp->gene)) ;      
		    } 
		  else
		    vtxtPrint (blkp, "sequence") ;
		  vtxtPrint (blkp, " of this gene") ;
		}
	      else
		vtxtPrint (blkp, "It ") ;
	      vtxtPrintf (blkp," is fully supported by %d cDNA clone%s" 
			  , numClones < 0 ? - numClones : numClones
			  , numClones < -1 || numClones > 1 ? "s" : ""
			  ) ;
	      if (decorate && gmp->markup && !gmp->mrna)
		vtxtPrintf (blkp,", the extent of which is shown by the turquoise bar") ;
	    }
	  break ;
	case 3:
	  if (decorate && numClones)
	    {
	      vtxtPrintf (blkp, "The existence of the gene, but not its exact ") ;
	      if (gmp->markup)
		{
		  vtxtPrintf (blkp, "<a href=\"javascript:openAnchor (0, 'Gene_sequences')\"><b>sequence</b></a>"
			      , ac_class (gmp->gene), ac_name (gmp->gene)) ;      
		} 
	      else
		vtxtPrintf (blkp, "sequence") ;
	      vtxtPrintf (blkp, ", derived here from the genome sequencing consortium annotation, is supported by %d cDNA clone%s"
			 , numClones < 0 ? - numClones : numClones
			 , numClones < -1 || numClones > 1 ? "s" : ""
			 ) ;
	      if (gmp->markup && !gmp->mrna)
		vtxtPrintf (blkp,", the extent of which is shown by the turquoise bar") ;
	    }
	  break ;
	case 2:
	  if (decorate)
	    {
	      vtxtPrintf (blkp, "The existence of the gene, but not its exact ") ;
	      if (gmp->markup)
		{
		  vtxtPrintf (blkp, "<a href=\"javascript:openAnchor (0, 'Gene_sequences')\"><b>sequence</b></a>"
			      , ac_class (gmp->gene), ac_name (gmp->gene)) ;      
		} 
	      else
		vtxtPrintf (blkp, "sequence") ;
	      vtxtPrintf (blkp, ", derived here from the genome sequencing consortium annotation, is supported by ") ;	 
	      if ((ptr = gtOst (gmp->gene)))
		{
		  char linkBuf[vONLINEMAXURLSIZE] ;
		  if (gmp->markup)
		    {
		      sprintf (linkBuf, "http://worfdb.dfci.harvard.edu/searchallwormorfs.pl?by=name&sid=%s",ptr);
		  gmpURL (blkp, gmp, linkBuf, "PCR amplification") ;
		    }
		  else
		    vtxtPrintf (blkp, "PCR amplification") ;
		  vtxtPrintf (blkp, " of the Vidal ORFeome project cDNA library") ;
		}
	      else
		vtxtPrintf (blkp, "PCR amplification of a cDNA library") ;
	    }
	  break ;
	case 1:
	  if (decorate)
	    {
	      vtxtPrintf (blkp, "The existence of the gene, but not its exact ") ;
	      if (gmp->markup)
		{
		  vtxtPrintf (blkp, "<a href=\"javascript:openAnchor (0, 'Gene_sequences')\"><b>sequence</b></a>"
			      , ac_class (gmp->gene), ac_name (gmp->gene)) ;      
		} 
	      else
		vtxtPrintf (blkp, "sequence") ;
	      vtxtPrintf (blkp, ", derived here from the genome sequencing consortium annotation, is supported by ") ;
	      vtxtPrintf (blkp, " expression data") ;
	    }
	  break ;
	case 0:	
	  if (decorate)
	    {
	      if (gmp->Spc != WORM)
		vtxtPrintf (blkp, "This predicted gene is not yet supported by cDNA or expression data") ; 
	      else
		vtxtPrintf (blkp, "It is predicted by the C.elegans Consortium (%s)"
			    " and not yet supported by cDNA or expression data"
			    , genomeRelease) ; 
	    }
	  break ;
	}
      ok++ ;
    }
   
  if (gmp->tg)
    {
      AC_KEYSET clones = ac_objquery_keyset (gmp->tg, ">cdna_clone", h) ;
      if (ac_keyset_count (clones))
	{
	  ficheNewGeneExpressionTissue (blkp, gmp, clones, 5, 8, gmp->view, ac_keyset_count (clones)) ;
	}
      ficheNewGeneCountAnomalousClones (blkp, gmp) ;
    }

  ac_free (h) ;
  return ok ;
} /* ficheNewGeneExpressionTissueStatement */

/***************************************************************************************/

static void ficheNewGeneSlStatement (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_ITER iter = gmp->tg ? ac_objquery_iter (gmp->tg, ">mrna ; Transpliced_to = SL* && Transpliced_to > SL0", h) : 0 ;
  AC_OBJ mrna = 0 ;;
  AC_TABLE tbl ;
  int k = 0, pos, ir, nn = strlen (ac_name(gmp->tg)) ;
  const char *variant = 0, *ccp = 0 ;

  while ((mrna = ac_iter_obj (iter)))
    {
      if (!strncmp (ac_name(gmp->gene), ac_name(mrna), nn) &&
	  *(ac_name(gmp->mrna) + nn) == '.')
	variant = ac_name (mrna) + nn + 1 ;
      
      if (!k)
	vtxtDot (blkp) ;

      if (variant)
	vtxtPrintf (blkp, "%s %s", k ? "," : "Variant", variant) ;
      else
	vtxtPrintf (blkp, "%s", ac_name(mrna)) ;
      vtxtPrintf (blkp, k ? " to" : " is transpliced to ") ;
      k++ ;
      tbl = ac_tag_table (mrna, "Transpliced_to", h) ;
      pos = vtxtMark (blkp) ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  ccp =  ac_table_printable (tbl, ir, 0, 0) ;
	  if (ccp && !strstr(ccp, "SL0") && (!vtxtAt (blkp, pos) || ! strstr (vtxtAt (blkp, pos), ccp)))
	    vtxtPrintf (blkp, "%s %s", ir ? "," : "", ccp) ;
	}
      ac_free (mrna) ;
    }
  ac_free (h) ;
} /*  ficheNewGeneSlStatement */

/***************************************************************************************/

static int ficheNewGeneAltVariantStatement (vTXT blkp, GMP *gmp, BOOL decorate)
{	
  int  level = 0, ok = 0, numMrnas = 0,  numMrnasWithIntron = 0 ;
  AC_KEYSET mrnasWithIntron = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  
  /*
  if (gmp->gene)
    {
      if ((numShed = ac_keyset_count (ac_objquery_keyset (gmp->gene, ">transcribed_gene  shedded_from", h))))
	numReadsShed = ac_keyset_count (ac_objquery_keyset (gmp->gene, ">transcribed_gene  shedded_from ; > read", h)) ;
    }
  */

  if (gmp->tg && gmp->gene)
    {
      level = 4 ;
       
      numMrnas = ac_keyset_count (ac_objquery_keyset (gmp->gene, ">transcribed_gene ; >mrna", h)) ;
      mrnasWithIntron = ac_objquery_keyset (gmp->gene, ">transcribed_gene ;>mrna ; gt_ag || gc_ag", h) ;
      numMrnasWithIntron = ac_keyset_count (mrnasWithIntron) ;
    }
  else if (gmp->gene && gmp->pg && ac_has_tag (gmp->pg, "CDS"))
    {
      numMrnas = ac_keyset_count (ac_objquery_keyset (gmp->gene, "Follow Genefinder", h)) ; 
      mrnasWithIntron = ac_objquery_keyset (gmp->gene, "Follow Genefinder; COUNT source_exons > 1", h) ; 
      numMrnasWithIntron = ac_keyset_count (mrnasWithIntron) ;
    }
   
  vtxtDot (blkp) ;
  vtxtPrint (blkp, " Transcription") ;
  if (numMrnas == 1)
    {
      ok++ ;
      vtxtPrintf (blkp
		  , " %s one mRNA"
		  , level ? " produces" : " would produce") ;
    }
  else if (numMrnas > 1)
    {
      ok++ ; 
      vtxtPrintf (blkp, 
		  "%s %s"
		  , level >= 4 ? " produces" : " would produce"
		  , gmp->style != 's' && (numMrnas > 4 && gmp->Spc == WORM) && gmp->tg ? "at least " : ""
		  ) ;  
      if (gmp->markup)
	{
	  char *cp ;
	  
	  if (numMrnasWithIntron > 1)
	    {
	      if (numMrnas == numMrnasWithIntron)
		cp = messprintf ("%d alternatively spliced mRNAs"
				 , numMrnas
				 ) ;
	      else 
		cp = messprintf ("%s different mRNA%s"
				 , _isone(numMrnas), _multi (numMrnas)
				 ) ;
	      vtxtPrintf (blkp, " <a href=\"javascript:openAnchor ('fmol','mRNAs')\">%s</a>", cp) ;
	    }
	  else
	    {
	      cp = "one spliced mRNA" ;
	      gmpObjLink (blkp, gmp, gmp->mrna, cp) ;
	    }
	}
      else
	vtxtPrintf (blkp, "%s different mRNA%s"
		    , _isone(numMrnas), _multi (numMrnas)) ;

      if (numMrnas != numMrnasWithIntron)
	vtxtPrintf (blkp, ", %d alternatively spliced variant%s and %d unspliced form%s"
		    , numMrnasWithIntron, _multi (numMrnasWithIntron)
		    , numMrnas - numMrnasWithIntron, _multi (numMrnas - numMrnasWithIntron)
		    ) ;
    }
  if (gmp->Spc == WORM)
    ficheNewGeneSlStatement (blkp, gmp) ;
  ac_free (h) ;
  return ok ;
} /* ficheNewGeneAltVariantStatement */

/***************************************************************************************/

static int  ficheNewGeneKozakStatement (vTXT blkp, GMP *gmp)
{	
  AC_HANDLE h = ac_new_handle () ;
  AC_KEYSET pk ;
  int nkoz = 0, ngain = 0 ;
  const char *myMet = 0 ;
  const char *pnam = 0, *myOkMet = 0 ;

  pk = ac_objquery_keyset (gmp->view == 'm' ? gmp->mrna : gmp->gene, ">product  Met && at_position_1 && First_Kozak = 1 && NEXT && First_ATG > 1 && best_product && good_product", h) ;
  if (pk && ac_keyset_count (pk))
    {
      AC_TABLE kozak ;
      AC_OBJ oProduct ;
      AC_TABLE pkt = ac_keyset_table (pk, 0, -1, 0, h) ;
      int ir, atg ;

      for (ir = 0 ; pkt && ir < pkt->rows ; ir++)
	{
	  oProduct = ac_table_obj (pkt, ir, 0, h) ;
	  kozak = ac_tag_table (oProduct, "First_Kozak", h) ;
	  if (kozak && ac_table_int (kozak, 0, 0, 0) == 1)
	    myMet = ac_table_printable (kozak, 0, 1, 0) ;
	  atg = ac_tag_int (oProduct, "First_ATG", 0) ;
	  if (!myMet || !atg)
	    continue ;
	  if (!nkoz++)
	    ngain = atg ;
	  else if (ngain > atg)
	    ngain = atg ;
	  pnam = ac_name (oProduct) ;
	  myOkMet = strnew (myMet, h) ;
	}
    }

  if (nkoz)
    {
      if (gmp->view == 'm')
	{
	  vtxtPrintf (blkp, ", annotated using as Met a ", pnam) ;
	  vtxtPrint (blkp, "<a href=\"javascript:openAnchor ('fmol','Proteins') \">Kozak-compatible</a>") ;
	  vtxtPrintf (blkp, " %s start, thereby gaining %d amino acids N-terminal to the first AUG", myOkMet, ngain) ;
	}
      else
	{
	  vtxtBreak (blkp) ;
	  if (nkoz == 1)
	    {
	      vtxtPrintf (blkp, "Isoform %s is annotated using as Met a ", pnam) ;
	      vtxtPrint (blkp, "<a href=\"javascript:openAnchor ('fmol','Proteins') \">Kozak-compatible</a>") ;
	      vtxtPrintf (blkp, " %s start, thereby gaining %d amino acids N-terminal to the first AUG", myOkMet, ngain) ;
	    }
	  else
	    {
	      vtxtPrintf (blkp, "%d isoforms are annotated using as Met a ", nkoz) ;
	      vtxtPrint (blkp, "<a href=\"javascript:openAnchor ('fmol','Proteins')\">Kozak-compatible</a>") ;
	      vtxtPrintf (blkp, " non-AUG start, thereby gaining a minimum of %d amino acids N-terminal to the first AUG", ngain) ;
	    }
	}
    }

  ac_free (h) ;
  return nkoz ;
} /* ficheNewGeneKozakStatement */

/***************************************************************************************/

static int  ficheNewGenePhosphositeStatement (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle() ;
  int nn = 0 ; 
  AC_ITER iter = 0 ;
  AC_OBJ obj = 0 ;
  char linkBuf [300] ;
  const char *modif, *uniprot ;
  
  if (gmp->gene && (iter = ac_objquery_iter (gmp->gene, ">Extern Phosphosite", h)))
    while (ac_free (obj), (obj = ac_iter_obj (iter)))
      {
	modif =  ac_tag_printable (obj, "Modif", 0) ;
	uniprot = ac_tag_printable (obj, "Uniprot", 0) ;
	if (! modif || ! uniprot)
	  continue ;
	if (!nn) 
	  {
	    vtxtDot (blkp) ;
	    vtxtPrintf (blkp, " Finally proteins from this gene may be modulated by ") ;
	    sprintf (linkBuf, "http://www.phosphosite.org/uniprotAccAction.do?id=%s", uniprot) ;
	    gmpURL (blkp, gmp, linkBuf, modif) ;
	    vtxtPrintf (blkp, ", as detailed at PhosphoSite") ;
	  }
      }
  ac_free (obj) ;
  ac_free (h) ;
  return nn ;
} /* ficheNewGenePhosphositeStatement */

/***************************************************************************************/

static int  ficheNewGeneUorfStatement (vTXT blkp, GMP *gmp)
{
  if (! gmp->gene || ! ac_has_tag (gmp->gene, "Pastille_regulation_uORF"))
    return 0 ;
  vtxtDot (blkp) ;
  vtxtPrintf (blkp, "Efficacy of translation may be reduced by the presence of a shorter translated"
	      " product  (") ;
  vtxtPrint (blkp, "<a href=\"javascript:openAnchor ('fmol','Proteins')\">uORF</a>") ;
  vtxtPrintf (blkp, ") initiating at an AUG upstream of the main open reading frame") ;
  
  {
    int nn = 0, kk = 0 ;
    AC_HANDLE h = ac_new_handle () ;
    AC_OBJ mrna = 0 ;
    AC_ITER iter = ac_objquery_iter (gmp->gene, ">product ; uORF_candidate; >mrna ; COUNT {>product; very_good_product} > 0", h) ;
    
    while ((mrna = ac_iter_obj (iter)))
      {
	if ((nn = strlen (ac_name(gmp->gene))) &&
	    !strncmp (ac_name(gmp->gene), ac_name(mrna), nn) &&
	    *(ac_name(mrna) + nn) == '.')
	  {
	    if (!kk++)
	      vtxtPrintf (blkp, " (in variant ") ;
	    else
	      vtxtPrintf (blkp, ", ") ;
	    vtxtPrintf (blkp, "%s", ac_name (mrna) + nn + 1, h) ;
	  }
	ac_free (mrna) ;
      }
    if (kk)
      vtxtPrintf (blkp, ")") ;
    ac_free (h) ;
  }

  return 1 ;
} /* ficheNewGeneUorfStatement */

/***************************************************************************************/

static int  ficheNewGeneNmdStatement (vTXT blkp, GMP *gmp, BOOL decorate)
{	
  int  nNmd = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  
  
  nNmd = ac_keyset_count (ac_objquery_keyset (gmp->tg, ">mrna ; NMD", h)) ;
  if (nNmd)
    {
      vtxtBreak (blkp) ;
      if (nNmd == 1)
	{
	  AC_OBJ nmdMrna = ac_objquery_obj (gmp->tg, ">mrna ; NMD", h) ;
	  const char *mrna = ac_name (nmdMrna) ;
	  
	  if (!strcasecmp (ac_name (gmp->gene), ac_name (nmdMrna)))
	    mrna = "" ;
	  if (strstr (ac_name (nmdMrna), ac_name (gmp->gene)))
	    mrna += strlen (ac_name (gmp->gene)) ;
	  
	  vtxtPrintf (blkp, "Note that mRNA %s was found ", mrna) ;
	  vtxtItalic (blkp, "in vivo") ;
	  vtxtPrintf (blkp, ", although it is a predicted target of") ;
	}
      else
	{
	  vtxtPrintf (blkp, "%d variants were isolated "
		      , nNmd
		      ) ;
	  vtxtItalic (blkp, "in vivo") ;
	  vtxtPrintf (blkp, ", despite the fact that they are predicted targets of") ;
	}
      vtxtPrintf (blkp
		  , " <a href=\"javascript:openAnchor ('fmol','Proteins')\">nonsense mediated mRNA decay</a> (NMD)") ;
    }

  ac_free (h) ;
  return nNmd ;
} /* ficheNewGeneNmdStatement */

/***************************************************************************************/

static int ficheNewGeneProteinStatement (vTXT blkp, GMP *gmp, BOOL decorate)
{	
  const char *prefix ;
  int ok = 0, numMrnas = 0,  numMrnasWithIntron = 0 ;
  int nGoodProduct, nCompleteProduct, nCoohCompleteProduct, nPartialProduct ;
  AC_KEYSET mrnasWithIntron = 0, products = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  
  
  if (gmp->gene)
      products = ac_objquery_keyset (gmp->gene, "{>transcribed_gene ; >mrna;>Product} SETELSE {>Product}", h) ;
  if (! products)
    goto done ;

  if (gmp->gene && gmp->tg)
    {
      numMrnas = ac_keyset_count (ac_objquery_keyset (gmp->gene, ">transcribed_gene ; >mrna", h)) ;
      mrnasWithIntron = ac_objquery_keyset (gmp->gene, ">transcribed_gene ;>mrna ; gt_ag || gc_ag", h) ;
      numMrnasWithIntron = ac_keyset_count (mrnasWithIntron) ;
    }
  else if (gmp->gene && gmp->pg && ac_has_tag (gmp->pg, "CDS"))
    {
      numMrnas = ac_keyset_count (ac_objquery_keyset (gmp->gene, "Follow Genefinder", h)) ; 
      mrnasWithIntron = ac_objquery_keyset (gmp->gene, "Follow Genefinder; COUNT source_exons > 1", h) ; 
      numMrnasWithIntron = ac_keyset_count (mrnasWithIntron) ;
    }
   
  vtxtBreak (blkp) ;
   
  nGoodProduct = ac_keyset_count (ac_ksquery_keyset (products, "(Best_product && Good_product) || very_good_product ; >kantor", h)) ;

  if (gmp->markup && numMrnas > 1 && numMrnasWithIntron && nGoodProduct)
    vtxtPrint (blkp, "<b>Protein coding potential:</b> ") ;

  /*
    numKantors = ac_keyset_count (ac_ksquery_keyset (products, "follow kantor", h)) ;
    numKantorsComplete = ac_keyset_count (ac_ksquery_keyset (products, " Complete; follow kantor", h)) ;
    numKantors5pComplete = ac_keyset_count (ac_ksquery_keyset (products, " !Complete && !COOH_Complete && NH2_Complete; follow kantor", h)) ;
    numKantors3pComplete = ac_keyset_count (ac_ksquery_keyset (products, " !Complete && COOH_Complete && !NH2_Complete; follow kantor", h)) ;
    numKantorsNotComplete = ac_keyset_count (ac_ksquery_keyset (products, " !Complete && !COOH_Complete && !NH2_Complete; follow kantor", h)) ;
  */
  

  if (!nGoodProduct)
    {
      vtxtPrintf (blkp, " No 'good' large or conserved protein is encoded by this gene, which might be considered non-coding") ;
      
      nCompleteProduct = ac_keyset_count (ac_ksquery_keyset (products, "Best_product && Complete ; >kantor", h)) ;
      nCoohCompleteProduct = ac_keyset_count (ac_ksquery_keyset (products, "Best_product && COOH_complete && ! Complete ; >kantor", h)) ;
      nPartialProduct = ac_keyset_count (ac_ksquery_keyset (products, "Best_product && ! COOH_complete ; >kantor", h)) ;
    }
  else  if (gmp->gene)
    {
      BOOL closeNeeded = FALSE ;
      int n0 = ac_keyset_count (ac_objquery_keyset (gmp->gene, ">transcribed_gene; >mrna", h)) ; 
      int n1 = ac_keyset_count (ac_objquery_keyset (gmp->gene, ">transcribed_gene;>mrna ; gt_ag || gc_ag ; COUNT {>product ; (Best_product && Good_product) || very_good_product} > 0", h)) ;
      int n2 = ac_keyset_count (ac_objquery_keyset (gmp->gene, ">transcribed_gene;>mrna ; !(gt_ag || gc_ag) ;COUNT {>product ; (Best_product && Good_product) || very_good_product} > 0", h)) ;
      if (n1 == 1 || n1 + n2 == n0) 
	vtxtPrintf (blkp, " The") ;
      if (n1 == 1)
	vtxtPrintf (blkp, " spliced", n1) ;
      if (n1 > 1)
	vtxtPrintf (blkp, " %d spliced", n1) ;
      if (n2 == 1)
	vtxtPrintf (blkp, " %s  unspliced", n1 ? "and the" : "") ;
      if (n2 > 1)
	vtxtPrintf (blkp, " %s%d unspliced", n1 ? "and " : "", n2) ;
      if (gmp->nMrna > 1)
	{
	  vtxtPrintf (blkp, " mRNAs putatively encode <a href=\"javascript:openAnchor ('fmol','Proteins')\">") ;
	  if (nGoodProduct > 1)
	    vtxtPrintf (blkp, "good proteins</a>") ;
	  else if (nGoodProduct == 1 && n1 + n2 == 1)
	    vtxtPrintf (blkp, "a good protein</a>") ;
	  else if (nGoodProduct == 1 && n1 + n2 > 1)
	    vtxtPrintf (blkp, "the same good protein</a>") ;
	}
      else
	{
	  vtxtPrintf (blkp, " mRNA putatively encodes ") ;
	  if (0)
	    gmpObjLinkAnchor (blkp, gmp, gmp->mrna, "Proteins", "a good protein") ;
	  else
	    vtxtPrintf (blkp, " <a href=\"javascript:openAnchor ('mrna','Proteins')\">a good protein</a>") ;
	}

      prefix = " (" ; 
      if (gmp->tg && nGoodProduct > 1)
	{
	  vtxtPrintf (blkp, ", altogether ") ;
	  /* pointe vers le tableau des variants */
	  if (gmp->markup)
	    {
	      closeNeeded = TRUE ;
	      vtxtPrint (blkp, "<a href=\"javascript:openAnchor ('fmol','Proteins')\">") ;
	    }
	  vtxtPrintf (blkp, "%d different isoforms", nGoodProduct) ;
	  	  
	  nCompleteProduct = ac_keyset_count (ac_ksquery_keyset (products, "((Best_product && Good_product) || very_good_product) && (Complete && at_position_1); >kantor", h)) ;
	  nCoohCompleteProduct = ac_keyset_count (ac_ksquery_keyset (products, "((Best_product && Good_product) || very_good_product) && COOH_complete && (! Complete || ! at_position_1) ; >kantor", h)) ;
	  nPartialProduct = ac_keyset_count (ac_ksquery_keyset (products, "((Best_product && Good_product) || very_good_product) && ! COOH_complete ; >kantor", h)) ;
	  	  
	  if (nCompleteProduct)
	    { vtxtPrintf (blkp, "%s%d complete", prefix , nCompleteProduct) ; prefix = ", " ; }
	  if (nCoohCompleteProduct)
	    { vtxtPrintf (blkp, "%s%d COOH complete", prefix , nCoohCompleteProduct) ; prefix = ", " ; }
	  if (nPartialProduct)
	    { vtxtPrintf (blkp, "%s%d partial", prefix , nPartialProduct) ; prefix = ", " ; }
	  if (closeNeeded)
	    vtxtPrintf (blkp, "</a>)") ;
	}
    }

 done:
  ac_free (h) ;
  return ok ;
} /* ficheNewGeneProteinStatement */

/***************************************************************************************/

static int ficheNewGenePfamPsortStatement (vTXT blkp, GMP *gmp, BOOL isCloud, BOOL fromSummary)
{	
  int nPfam = 0, nProducts ;
  int nPsort_domain = 0 ;
  AC_KEYSET products = 0 ;
  vTXT bfr = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  
  const char *prefix = " " ;
  if (gmp->gene)
    products = ac_objquery_keyset (gmp->gene, ">Product ; (best_product && good_product) || very_good_product", h) ;

  if (! products)
    goto done ;

  nProducts = ac_keyset_count (products) ;
  /* PFAM sentence */
  bfr = vtxtHandleCreate (h) ;
  if (1)
    {
      const char *txtAccession ;
      char familyName [1000] ;
      vTXT pfamBfr = vtxtHandleCreate (h) ;
      AC_OBJ oPfam = 0 ;
      AC_ITER pfamIter = ac_ksquery_iter (products, ">Pfam", h) ;

      while (ac_free (oPfam), oPfam = ac_iter_obj (pfamIter))
	{
	  if (vtxtPtr (pfamBfr) && strstr (vtxtPtr (pfamBfr), ac_name (oPfam)))
	    continue ;
	  nPfam++ ;
	  vtxtPrintf (pfamBfr, "%s", ac_name (oPfam)) ;
	  strncpy (familyName, ac_tag_printable (oPfam, "Definition", ""), 999) ; 
	  if (!*familyName)
	    {
	      /* try to get it from a synonymous pfam */
	      AC_OBJ oPfam2 = ac_objquery_obj (oPfam, ">accession; >quoted_in; Definition", h) ;
	      
	      strcpy (familyName, ac_tag_printable (oPfam2, "Definition", "")) ;
	      ac_free (oPfam2) ;
	    }
	  
	  vtxtPrint (bfr, prefix) ;
	  prefix = ", " ;
	  if (0 &&
	      gmp->markup &&
	      (txtAccession = ac_tag_printable (oPfam, "Accession", 0)))
	    {
	      char linkBuf [300] ;
	      char *cq = strstr (txtAccession , ".") ; /* drop pfam version number */
	      if (cq) *cq = 0 ;
	      sprintf (linkBuf, "http://pfam.xfam.org/family/%s", txtAccession) ; 
	      gmpURL (bfr, gmp, linkBuf, familyName) ; 
	    }
	  else if (*familyName)
	    vtxtPrint (bfr, familyName) ; 
	  else 
	    vtxtPrint (bfr, "Pfam") ; 
	}
      ac_free (oPfam) ;
    }
  if (nPfam == 1)
    {
      if (!strstr (vtxtPtr (bfr), "domain") && !strstr (vtxtPtr (bfr), "motif"))
	{
	  if (fromSummary)
	    vtxtPrint (bfr, " <a href='javascript:openAnchor (\"ffunc\",\"tg_product_pfam\") '>domain</a> [Pfam]") ;
	  else
	    vtxtPrint (bfr, " <a href='javascript:openAnchor (\"fmol\",\"Proteins\") '>domain</a> [Pfam]") ;
	}
      else
	{
	  if (fromSummary)
	    vtxtPrint (bfr, " <a href='javascript:openAnchor (\"ffunc\",\"tg_product_pfam\") '>[Pfam]</a>") ;
	  else
	    vtxtPrint (bfr, " <a href='javascript:openAnchor (\"fmol\",\"Proteins\") '>[Pfam]</a>") ;
	}
    }
  else if (nPfam > 1)
    vtxtPrint (bfr, " [Pfam]") ;

  /* Psort_domain sentence */
  if (1)
    {
      if (ac_has_tag (gmp->product, "Psort_domain"))
	{
	  int ii, n1 = 0 ;
	  char **ccp, *psNN [] =
	    { 
	      "Transmembrane_domain", "some transmembrane domains", 
	      "N_myristoylation_domain", "a N-myristoylation domain",
	      "prenylation_domain", "a prenylation domain",
	      "Golgi_transport_domain", "a Golgi transport domain", 

	      "peroxisomal_domain", "a peroxisomal domain", 
	      "2nd_peroximal_domain", "a second peroximal domain", 
	      "vacuolar_domain", "a vacuolar domain", 

	      "ER_retention_domain", "an ER_retention domain", 
	      "Coiled_coil_region", "a coiled coil stretch", 
	      "actin_binding_1_domain", "an actin binding domain", 
	      "actin_binding_2_domain", "an actin binding domain", 
	      "ATP_binding_domain", "an ATP binding domain",                              
	      0, 0
	    } ;

	  for (ii = 0, ccp = psNN ; *ccp ; ii++, ccp += 2)
	    {
	      n1 = ac_keyset_count (ac_ksquery_keyset (products
						       , messprintf ("((Best_product && Good_product) || very_good_product) && %s", *ccp)
						       , h)
				) ;
	      if (n1 && ! (vtxtPtr (bfr) && strstr (vtxtPtr (bfr), *(ccp+1))))
		{
		  vtxtPrintf (bfr, "%s%s", prefix, *(ccp+1)) ;
		  prefix = ", " ; nPsort_domain++ ;
		}
	    }
	  if (n1) vtxtPrint (bfr, "</a>") ;
	}
    }
  if (nPsort_domain)
    {
      if (fromSummary)
	vtxtPrint (bfr, " <a href='javascript:openAnchor (\"ffunc\",\"tg_product_pfam\") '>[Psort2]</a>") ;
      else
	vtxtPrint (bfr, " <a href='javascript:openAnchor (\"fmol\",\"Proteins\") '>[Psort2]</a>") ;
    }
  if (vtxtPtr (bfr))
    {
      if (isCloud)
	vtxtPrintf (blkp, ". This gene contains%s", nProducts > 2 ? " some" : "") ;
      else
	vtxtPrintf (blkp, ", %scontaining", nProducts > 2 ? "some " : "" ) ;
      if (nPfam > 1)
	{
	  if (fromSummary)
	    vtxtPrintf (blkp, " <a href='javascript:openAnchor (\"ffunc\",\"tg_product_pfam\") '>domains</a>") ;
	  else
	    vtxtPrintf (blkp, " <a href='javascript:openAnchor (\"fmol\",\"Proteins\") '>domains</a>") ;
	}

      vtxtPrintf (blkp, " %s", vtxtPtr (bfr)) ;
    }
  
  /* secreted */
  {
    int n1 = ac_keyset_count (ac_ksquery_keyset (products, "((Best_product && Good_product) || very_good_product) && Complete && Psort_title = *secreted* ; > Kantor", h)) ;
    int n2 = ac_keyset_count (ac_ksquery_keyset (products, "((Best_product && Good_product) || very_good_product) && Complete ; > Kantor", h)) ;
    if (n1)
      {
	if (n1 < n2)
	  vtxtPrintf (blkp, "; %d of the %d complete proteins appear%s to be", n1, n2, n1 > 1 ? "" : "s") ;
	else if (n1 == n2 && n1 > 1)
	  vtxtPrintf (blkp, "; the %d complete proteins appear to be", n1) ;
	else if (n1 == n2 && n1 == 1)
	  vtxtPrint (blkp, "; the complete protein appears to be") ;
	{
	  if (fromSummary)
	    vtxtPrint (blkp, " <a href='javascript:openAnchor (\"ffunc\",\"tg_pathway\") '>secreted</a>") ;
	}
      }
  }
  /* tax_common */
  {
    int n1 = ac_keyset_count (ac_objquery_keyset (gmp->gene
						  , messprintf ("COUNT {>product; complete && ((Best_product && Good_product) || very_good_product) && Tax_common_ancestor == \"%s\"} > 0 &&  COUNT {>product; ((Best_product && Good_product) || very_good_product) && Tax_common_ancestor != \"%s\"} == 0"
								, gmp->spci->speciesName, gmp->spci->speciesName )
						  , h)) ;
    if (n1 > 0)
      vtxtPrintf (blkp, ", apparently %s specific",  gmp->spci->speciesName) ;
    else
      {
	n1 = ac_keyset_count (ac_objquery_keyset (gmp->gene
						  , messprintf ("COUNT {>product; ((Best_product && Good_product) || very_good_product) && Tax_common_ancestor == \"eutheria\"} > 0 &&  COUNT {>product; ((Best_product && Good_product) || very_good_product)  && Tax_common_ancestor != \"eutheria\"  && Tax_common_ancestor != \"%s\"} == 0"
								, gmp->spci->speciesName)
						  , h)) ;
	if (n1 > 0)
	  vtxtPrintf (blkp, ", apparently mammals specific",  gmp->spci) ;
	else
	  {
	    n1 = ac_keyset_count (ac_objquery_keyset (gmp->gene
						      , messprintf ("COUNT {>product; ((Best_product && Good_product) || very_good_product) && Tax_common_ancestor == \"teleostom*\"} > 0 &&  COUNT {>product; ((Best_product && Good_product) || very_good_product) && Tax_common_ancestor != \"teleostom*\"  && Tax_common_ancestor != \"eutheria\"  && Tax_common_ancestor != \"%s\"} == 0"
								    , gmp->spci->speciesName)
						      , h)) ;
	    if (n1 > 0)
	      vtxtPrintf (blkp, ", apparently vertebrate specific",  gmp->spci) ;
	  }
      }
  }
 done:
  ac_free (h) ;
  return 1 ;
} /* ficheNewGenePfamPsortStatement */

/***************************************************************************************/

static int ficheNewGeneNonGoodVariantStatement (vTXT blkp, GMP *gmp)
{	
  const char *prefix ;
  AC_HANDLE h = ac_new_handle () ;  
  int n000 = ac_keyset_count (ac_objquery_keyset (gmp->tg, ">mrna", h)) ;
  int n00 = ac_keyset_count (ac_objquery_keyset (gmp->tg, ">mrna ;  COUNT {>product ; Good_product;} > 0", h)) ;
  int n0 = ac_keyset_count (ac_objquery_keyset (gmp->tg, ">mrna ;  COUNT {>product ; Good_product;} == 0", h)) ;
  int n1 = ac_keyset_count (ac_objquery_keyset (gmp->tg, ">mrna ; gt_ag || gc_ag ; COUNT {>product ; Good_product;} == 0", h)) ;
  int n2 = ac_keyset_count (ac_objquery_keyset (gmp->tg, ">mrna ; !(gt_ag || gc_ag) ; COUNT {>product ; Good_product;} == 0", h)) ;
  int n3 = ac_keyset_count (ac_objquery_keyset (gmp->tg, ">mrna ; COUNT {>product ; Good_product;} == 0  && COUNT { >product; Best_product && ! complete} > 0  ", h)) ;
  
  if (n0)
    {
      if (n000 > 1 && n0 == 1)
	{
	  vtxtPrintf (blkp, ". The %s mRNA variant", n00 ? "remaining " : "") ;
	  prefix = " (" ;
	  if (n1)
	    { vtxtPrintf (blkp, "%sspliced", prefix) ; prefix = ", " ; }
	  if (n2)
	    { vtxtPrintf (blkp, "%sunspliced", prefix) ; prefix = ", " ; }
	  if (*prefix == ',') prefix = "; " ;
	  if (n3)
	    { vtxtPrintf (blkp, "%spartial", prefix) ; prefix = ", " ; }
	  if (*prefix == ',' || *prefix == ';')
	    vtxtPrintf (blkp, ")") ;
	  vtxtPrint (blkp, " appears not to encode a good protein") ;
	}
      else if (n0 > 1)
	{
	  vtxtPrintf (blkp, ". The %s%d mRNA variants", n00 ? "remaining " : "", n0) ;
	  prefix = " (" ;
	  if (n1)
	    { vtxtPrintf (blkp, "%s%d spliced", prefix, n1) ; prefix = ", " ; }
	  if (n2)
	    { vtxtPrintf (blkp, "%s%d unspliced", prefix, n2) ; prefix = ", " ; }
	  if (*prefix == ',') prefix = "; " ;
	    if (n3)
	    { vtxtPrintf (blkp, "%s%d partial", prefix, n3) ; prefix = ", " ; }
	  if (*prefix == ',' || *prefix == ';')
	    vtxtPrintf (blkp, ")") ;
	  vtxtPrint (blkp, " appear not to encode good proteins") ;
	}
    }

  ac_free (h) ;
  return 1 ;
} /* ficheNewGeneNonGoodVariantStatement */

/***************************************************************************************/

static int ficheNewGeneComplexLocusStatement (vTXT blkp, GMP *gmp)
{	
  if (ac_has_tag (gmp->gene, "Complex_locus"))
    { 
      BOOL b1 = ac_has_tag (gmp->gene, "Two_product") ;
      BOOL b2 = ac_has_tag (gmp->gene, "Two_geneId") ; 
      AC_HANDLE h = ac_new_handle () ;

      int n1 = ac_keyset_count (ac_objquery_keyset (gmp->gene, ">GeneId", h)) ;
   
      vtxtBreak (blkp) ;
      if (b1)
	{
	  vtxtPrintf (blkp, "Note that this locus is complex: it appears to produce several proteins with no sequence overlap") ;
	  if (b2)
	    vtxtPrintf (blkp, ", and it merges %d different NCBI genes", n1) ;
	}
      else if (b2)
	{
	  vtxtPrintf (blkp, "Note that this locus merges %d  different NCBI genes, because some GenBank cDNA has significant sequence overlap or common intron boundaries with %s gene models", n1, n1==2 ? "both" : "several") ;
	}
      ac_free (h) ;  
    }
  return 1 ;
} /* ficheNewGeneComplexLocusStatement */

/***************************************************************************************/
/***************************************************************************************/
/* same query copied in BLY q2c to compute pastille_disease */
static AC_KEYSET ficheNewGeneDiseaseKs (AC_OBJ gene, AC_HANDLE h)
{
  return ac_objquery_keyset (gene, "{>disease} $| {>Mesh} $| {>Reference ; mesh && COUNT gene < 3 ; >Mesh  meshkey } $| {>Extern NOT AntiGad;>Disease} $| {>Extern ; KEGG_disease ; >pathway } $| {>Extern ; GAD AND NOT AntiGad ; >GAD_TITLE NOT AntiGad} $| {>Extern ; OMIM_disease; >OMIM_title } ; {!Hidden_alias_of} $| {>Hidden_alias_of} ; {! alias_of} $| {>alias_of} ;  {IS *} SETMINUS {>meshkey ; >child ; >mesh} ; !meshkey || meshkey = C* || meshkey = F*" , h) ;
} 
   
static vTXT ficheNewGeneFunctionDisease (GMP *gmp, int *np, AC_HANDLE h)
{
  AC_KEYSET ks = 0, ks1 = 0 , ks2 = 0 ;
  AC_TABLE tbl ;
  int pass, ir, jr = 0, n1 = 0 ;
  vTXT txt = 0 ;
  
  *np = 0 ; 
  
  for (pass = 0 ; pass < 2 ; pass++)
    {
      switch (pass)
	{
	case 0:
	  ks1 =  ac_objquery_keyset (gmp->gene, " {>Extern ; OMIM_disease; >OMIM_title }" , h) ;
	  ks = ks1 ;
	  break ;
	case 1:
	  ks2 =  ficheNewGeneDiseaseKs (gmp->gene, h) ;
	  if (ks1) ac_keyset_minus (ks2, ks1) ;
	  ks = ks2 ;
	  break ;
	} 
      tbl = ks ? ac_keyset_table (ks, 0, -1, 0, h) : 0 ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  if (!txt)
	    txt = vtxtHandleCreate (h) ;
	  if (ir == 0 && pass == 0) vtxtPrint (txt, "OMIM: ") ;
	  else if (ir == 0 && jr > 0 && pass == 1) vtxtPrint (txt, "; Other sources: ") ;
	  else if (ir)
	    vtxtPrint (txt, "; ") ;
	  (*np)++ ;  jr++ ;   
	  vtxtPrint (txt, ac_table_printable (tbl, ir, 0, "")) ;
	  if (pass > 0 && jr > 8 && tbl->rows - ir > 4) break ;
	}
    }
  n1 = tbl ? tbl->rows - ir - 1 : 0 ;
  if (n1 > 0)
    vtxtPrintf (txt, " and <a href='javascript:openAnchor (\"ffunc\",\"tg_disease\")'>%d other%s</a>", n1, _multi(n1)) ;
  return txt ;
} /* ficheNewGeneFunctionDisease */

/********************/

static vTXT ficheNewGeneFunctionProcess (GMP *gmp, int *np, AC_HANDLE h)
{
  AC_TABLE tbl ;
  int pass, ir, jr, jr1, n1 = 0 ;
  vTXT txt = 0 ;
  const char *ccp ;
  AC_KEYSET ks ;

  jr = jr1 = 0 ;
  for (pass = 0 ; pass < 3 ; pass++)
    {
      switch (pass)
	{
	case 0:
	  ks = ac_objquery_keyset (gmp->gene, "{>Go_b_ace}", h) ;
	  break ;
	case 1:	
	  ks = ac_objquery_keyset (gmp->gene, "{>Go_b_iea}", h) ;
	  break ;
	case 2:	
	  ks = ac_objquery_keyset (gmp->gene, "{>Go_b_pfam}", h) ;
	  break ;
	}
      tbl = ac_keyset_table (ks, 0, -1, 0, h) ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  ccp = ac_table_printable (tbl, ir, 0, "") ;
	  if (strstr (ccp, "unknown"))
	    continue ;
	  jr++ ; 
	}
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  ccp = ac_table_printable (tbl, ir, 0, "") ;
	  if (strstr (ccp, "unknown"))
	    continue ;
	  if (txt && vtxtPtr (txt) && strstr (vtxtPtr (txt), ccp))
	    { jr-- ; continue ;}
	  if (jr1 > 34 && jr > 38) 
	    break ;
	  jr1++ ;
	  
	  if (!txt)
	    txt = vtxtHandleCreate (h) ;
	  else
	    vtxtPrint (txt, ", ") ;
	  vtxtPrint (txt, ccp) ;
	}
    }
  *np = jr ;
  n1 = jr - jr1 ;
  if (n1 > 0)
    {
      if (!txt)
	txt = vtxtHandleCreate (h) ;
      vtxtPrintf (txt, " and <a href='javascript:openAnchor (\"ffunc\",\"tg_pathway\")'>%d other%s</a>", n1, _multi(n1)) ;
    }
  return txt ;
} /* ficheNewGeneFunctionProcess */

/********************/

static vTXT ficheNewGeneFunctionFunction (GMP *gmp, int *np, int *np1, AC_HANDLE h)
{
  AC_KEYSET ks ;
  AC_TABLE tbl ;
  int pass, ir, jr, jr1, n1 = 0 ;
  vTXT txt = 0 ;
  const char *ccp ;

  jr = jr1 = 0 ;
  for (pass = 0 ; pass < 3 ; pass++)
    {
      switch (pass)
	{
	case 0:
	  ks = ac_objquery_keyset (gmp->gene, "{>Go_m_ace}", h) ;
	  break ;
	case 1:	
	  ks = ac_objquery_keyset (gmp->gene, "{>Go_m_iea}", h) ;
	  break ;
	case 2:	
	  ks = ac_objquery_keyset (gmp->gene, "{>Go_m_pfam}", h) ;
	  break ;
	}
      
      tbl = ac_keyset_table (ks, 0, -1, 0, h) ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  ccp = ac_table_printable (tbl, ir, 0, "") ;
	  if (strstr (ccp, "unknown"))
	    continue ;
	  jr++ ; 
	}
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  ccp = ac_table_printable (tbl, ir, 0, "") ;
	  if (strstr (ccp, "unknown"))
	    continue ;
	  if (txt && vtxtPtr (txt) && strstr (vtxtPtr (txt), ccp))
	    { jr-- ; continue ;}
	  if (jr1 > 30 && jr > 35) 
	    break ;
	  jr1++ ;
	  
	  if (!txt)
	    txt = vtxtHandleCreate (h) ;
	  else
	    vtxtPrint (txt, ", ") ;
	  vtxtPrint (txt, ccp) ;
	}
    }
  *np = jr ; 
  *np1 = jr1 ; 
  n1 = jr - jr1 ;
  if (n1 > 0)
    {
      if (!txt)
	txt = vtxtHandleCreate (h) ;
      vtxtPrintf (txt, " and <a href='javascript:openAnchor (\"ffunc\",\"tg_pathway\")'>%d other%s</a>", n1, _multi(n1)) ;
    }
  return txt ;
} /* ficheNewGeneFunctionFunction */

/********************/

static vTXT ficheNewGeneFunctionLocalization (GMP *gmp, int *np, int *np1, AC_HANDLE h)
{
  AC_KEYSET ks ;
  AC_TABLE tbl ;
  int pass, ir, jr, jr1, n1 = 0 ;
  vTXT txt = 0 ;
  const char *ccp ;

  jr = jr1 = 0 ;
  for (pass = 0 ; pass < 4 ; pass++)
    {
       switch (pass)
	{
	case 0:
	  ks = ac_objquery_keyset (gmp->gene, "{>Go_c_ace}", h) ;
	  break ;
	case 1:	
	  ks = ac_objquery_keyset (gmp->gene, "{>Go_c_psort}", h) ;
	  break ;
	case 2:	
	  ks = ac_objquery_keyset (gmp->gene, "{>Go_c_iea}", h) ;
	  break ;
	case 3:	
	  ks = ac_objquery_keyset (gmp->gene, "{>Go_c_pfam}", h) ;
	  break ;
	}
 
      tbl = ac_keyset_table (ks, 0, -1, 0, h) ;
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  ccp = ac_table_printable (tbl, ir, 0, "") ;
	  if (strstr (ccp, "unknown"))
	    continue ;
	  jr++ ; 
	}
      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  ccp = ac_table_printable (tbl, ir, 0, "") ;
	  if (strstr (ccp, "unknown"))
	    continue ;
	  if (txt && vtxtPtr (txt) && strstr (vtxtPtr (txt), ccp))
	    { jr-- ; continue ;}
	  if (jr1 > 3 && jr > 7) 
	    break ;
	  jr1++ ;
	  
	  if (!txt)
	    txt = vtxtHandleCreate (h) ;
	  else
	    vtxtPrint (txt, ", ") ;
	  vtxtPrint (txt, ccp) ;
	}
    }
  *np = jr ; 
  *np1 = jr1 ; 
  n1 = jr - jr1 ;
  if (n1 > 0)
    {
      if (!txt)
	txt = vtxtHandleCreate (h) ;
      vtxtPrintf (txt, " and <a href='javascript:openAnchor (\"ffunc\",\"tg_pathway\")'>%d other%s</a>", n1, _multi(n1)) ;
    }
  return txt ;
} /* ficheNewGeneFunctionLocalization */

/********************/

static vTXT ficheNewGeneFunctionPathway (GMP *gmp, int *np, const char **keggp, AC_HANDLE h)
{
  AC_KEYSET ks ;
  AC_TABLE tbl ;
  int ir, n1 = 0 ;
  vTXT txt = 0 ;

  *np = 0 ; 
  ks = ac_objquery_keyset (gmp->gene, "{>Extern ; (KEGG && ! KEGG_disease) ; >pathway} SETOR {>geneid;>extern; >Biocarta_title}", h) ;
  tbl = ks ? ac_keyset_table (ks, 0, -1, 0, h) : 0 ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      if (!txt)
	txt = vtxtHandleCreate (h) ;
      else
	vtxtPrint (txt, ", ") ;
      (*np)++ ;      
      vtxtPrint (txt, ac_table_printable (tbl, ir, 0, "")) ;
      if (ir > 100 && tbl->rows > 12) break ;
    }
  n1 = tbl ? tbl->rows - ir - 1 : 0 ;
  if (n1 > 0)
    vtxtPrintf (txt, " and %d other%s", n1, _multi(n1)) ;
  if (keggp && tbl && tbl->rows == 1)
    {
      AC_ITER iter = ac_objquery_iter (gmp->gene, ">Extern ; IS kegg_* && KEGG && ! KEGG_disease", h) ;
      AC_OBJ obj = iter ? ac_iter_obj (iter) : 0 ;
      if (obj)
	{
	  *keggp = ac_name (obj) + 5 ;  
	  ac_free (obj) ;
	}
    }
  return txt ;
} /* ficheNewGenePathway */

/********************/

static vTXT ficheNewGeneFunctionInteracts (GMP *gmp, int *np, BOOL with_protein, AC_HANDLE h)
{
  AC_KEYSET ks ;
  AC_TABLE tbl ;
  int ir, n1 = 0 ;
  vTXT txt = 0 ;

  *np = 0 ; 
  if (with_protein)
    ks = ac_objquery_keyset (gmp->gene, "{>with_protein} SETELSE {>geneid ; >with_geneid; >gene LocusLink}", h) ;
  else
    ks = ac_objquery_keyset (gmp->gene, ">with_gene", h) ;
  tbl = ks ? ac_keyset_table (ks, 0, -1, 0, h) : 0 ;
  *np = tbl ? tbl->rows : 0 ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      if (!txt)
	txt = vtxtHandleCreate (h) ;
      else
	vtxtPrint (txt, ", ") ;
      vtxtPrint (txt, ac_table_printable (tbl, ir, 0, "")) ;
      if (ir > 100 && tbl->rows > 15) break ;
    }
  if (tbl && tbl->rows)
    vtextUpperCase (vtxtPtr (txt)) ;
  n1 = tbl ? tbl->rows - ir - 1 : 0 ;
  if (n1 > 0)
    vtxtPrintf (txt, " and <a href='javascript:openAnchor (\"ffunc\",\"%s\")'>%d other%s</a>"
		, with_protein ? "tg_protein_interaction" : "tg_gene_interaction" 
		, n1, _multi(n1)
		) ;
  return txt ;
} /* ficheNewGenePathwayInteracts */

/********************/

static int ficheNewGeneFunctionStatement (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn , nd = 0, npr = 0, npf = 0, npw = 0, nloc = 0, nloc_ph = 0, npf1 = 0, nloc1 = 0, nitp = 0, nitg = 0 ;
  int n1 ;
  int ngp = ac_keyset_count (ac_objquery_keyset (gmp->gene, ">Product;good_product;From_gene", h)) ;
  const char *keggp = 0 ;
  AC_KEYSET ksPaper = ac_objquery_keyset (gmp->gene, ">Reference", h) ;
  vTXT diseaseTxt = ficheNewGeneFunctionDisease (gmp, &nd, h) ;
  vTXT pathwayTxt = ficheNewGeneFunctionPathway (gmp, &npw, &keggp, h) ;
  vTXT processTxt = ficheNewGeneFunctionProcess (gmp, &npr, h) ;
  vTXT functionTxt = ficheNewGeneFunctionFunction (gmp, &npf, &npf1, h) ;
  vTXT locTxt = ficheNewGeneFunctionLocalization (gmp, &nloc, &nloc1, h) ;
  vTXT interactsPTxt = ficheNewGeneFunctionInteracts (gmp, &nitp, TRUE, h) ;
  vTXT interactsGTxt = ficheNewGeneFunctionInteracts (gmp, &nitg, FALSE, h) ;
  vTXT locPhTxt = ficheNewGeneFunctionLocusPhenotypeDescription (gmp, &nloc_ph, h) ;
  if ((gmp->Spc == ARA || gmp->Spc == WORM) && 
       (ac_has_tag (gmp->gene, "Locus_phenotype") || ac_has_tag (gmp->gene, "Locus_description"))
      )
    nloc_ph = 1 ;

  n1 = ac_keyset_count (ksPaper) ;
  nn = n1 + nd + npr + npw + nloc_ph ;
  if (nn) /* sometehing to say */
    {
      vtxtBreak (blkp) ;
      if (gmp->markup)
	vtxtPrint (blkp, "<b>Function:</b> ") ;
    }
  if (n1)
    {
      if (gmp->markup)
	{
	  AC_TABLE tbl = ac_keyset_table (ksPaper, 0, -1, 0, h) ;
	  vTXT bfr1 = vtxtHandleCreate (h) ;
	  int ir, ir1, j = 0 ;
	    
	  if (gmp->markup) vtxtMarkup (bfr1) ;   

	  vtxtPrint (bfr1, PUBMED_LINK) ; /* no %s included */
	  for (ir = tbl->rows - 1, ir1 = 0 ; ir >= 0 ; ir--)
	    {
	      const char *cp = ac_table_printable (tbl, ir, 0, 0) ;
	      if (cp && !strncasecmp (cp, "pmp", 2))
		{
		  if (ir1++) { if (j<300) vtxtPrint (bfr1, ",") ; }
		  if (j++ < 300)vtxtPrint (bfr1, cp + 2) ;
		}
	    }
	  vtxtPrintf (blkp, "There %s ", _isare(ir1)) ;
	  gmpURL (blkp, gmp, vtxtPtr (bfr1)
		  , messprintf ("%s article%s", isOne(ir1), _multi (ir1))) ;
	  vtxtPrintf (blkp, " specifically referring to this gene in PubMed") ;
	  if (n1 > ir1)
	    {
	      vtxtPrintf (blkp, ". In addition we point <a href='javascript:openAnchor (0,\"Biblio\") '>below</a> to %s abstract%s"
			  , isOne (n1 - ir1), _multi (n1 - ir1)) ;

	    }
	}		  
      else
	vtxtPrintf (blkp, "There %s %d article%s specifically referring to this gene in PubMed"
		    , n1 == 1 ? "is a" : "are"
		    , n1
		    , _multi (n1)
		    ) ;
    }
  
  if (nloc_ph)
    { 
      vtxtDot (blkp) ;
      if (locPhTxt)
	vtxtPrint (blkp, vtxtPtr (locPhTxt)) ;
    }

  if (nd + npr + npw > 0)
    {
      vtxtDot (blkp) ;
      vtxtPrintf (blkp, "Functionally, the gene has been %s "
		  , nd ? "tested for association to " : "proposed to participate in "
		  ) ;
      n1++ ;
      if (nd)
	{
	  if (gmp->markup)
	    vtxtPrintf (blkp, "%s<a href='javascript:openAnchor (\"ffunc\",\"tg_disease\")'>disease%s</a> ("
			,  nd ==1 ? "a " : ""
			, _multi (nd)
			) ;
	  else
	    vtxtPrintf (blkp, "%sdisease%s ("
			,  nd ==1 ? "a " : ""
			, _multi (nd)
			) ;
	  vtxtPrint (blkp, vtxtPtr (diseaseTxt)) ;
	  vtxtPrint (blkp, ")") ;
	}
      
      if (npw)
	{
	  if (nd && npr)
	    vtxtPrint (blkp, ", proposed to participate in ") ;
	  if (nd && !npr)
	    vtxtPrint (blkp, " and proposed to participate in ") ;
	  if (gmp->markup)
	    {
	      if (keggp)
		{
		  char linkBuf[1000] ;
		  AC_OBJ gid = ac_tag_obj (gmp->gene, "geneid", h) ;
		  if (gmp->Spc == HUMAN)
		    sprintf (linkBuf, KEGG_LINK_HUMAN, keggp, gid ? ac_name(gid) : "") ;
		  else if (gmp->Spc == MOUSE)
		    sprintf (linkBuf, KEGG_LINK_MOUSE, keggp, gid ? ac_name(gid) : "") ;
		  else if (gmp->Spc == RAT)
		    sprintf (linkBuf, KEGG_LINK_RAT, keggp, gid ? ac_name(gid) : "") ;
		  else if (gmp->Spc == ARA)
		    sprintf (linkBuf, KEGG_LINK_ARA, keggp, gid ? ac_name(gid) : "") ;
		  gmpURL (blkp, gmp, linkBuf, "pathway") ;
		}
	      else
		vtxtPrintf (blkp, "%s<a href='javascript:openAnchor (\"ffunc\",\"tg_pathway\")'>pathway%s</a>"
			    ,  npw ==1 ? "a " : ""
			    , npw > 1 ? "s" : ""
			    ) ;
	    }
	  else
	    vtxtPrintf (blkp, "%spathway%s"
			, npw==1 ? "a " : ""
			, npw > 1 ? "s" : ""
			) ;
	  vtxtPrint (blkp, " (") ;
	  vtxtPrint (blkp, vtxtPtr (pathwayTxt)) ;
	  vtxtPrint (blkp, ")") ;
	}
      if (npr)
	{
	  if (npw)
	    vtxtPrint (blkp, " and ") ;
	  else if (nd)
	    vtxtPrint (blkp, " and proposed to participate in ") ;
	  if (gmp->markup)
	    vtxtPrintf (blkp, "%s<a href='javascript:openAnchor (\"ffunc\",\"tg_pathway\")'>process%s</a> ("
			, npr ==1 ? "a " : ""
			, npr > 1 ? "es" : ""
			) ;
	  else
	    vtxtPrintf (blkp, "%sprocess%s ("
			, npr ==1 ? "a " : ""
			, npr > 1 ? "es" : ""
			) ;
	  vtxtPrint (blkp, vtxtPtr (processTxt)) ;
	  vtxtPrint (blkp, ")") ;
	}
    }
  
  if (npf1)
    {
      vtxtPrint (blkp, n1 ? ". " : "") ;
      vtxtPrint (blkp, " Proteins are expected to have molecular ") ;
      if (gmp->markup)
	vtxtPrintf (blkp, "<a href='javascript:openAnchor (\"ffunc\",\"tg_pathway\")'>function%s</a> "
		    , npf1 > 1 ? "s" : ""
		    ) ;
      else
	vtxtPrintf (blkp, "function%s "
		    , npf1 > 1 ? "s" : ""
		    ) ;
      vtxtPrint (blkp, "(") ;
      vtxtPrint (blkp, vtxtPtr (functionTxt)) ;
      vtxtPrint (blkp, ")") ;
    }

   if (nloc1)
    {
      if (npf)
	vtxtPrint (blkp, " and to") ;
      else
	{
	  vtxtPrint (blkp, n1 ? ". " : "") ;
	  vtxtPrint (blkp, " Proteins are expected to") ;
	}
      if (gmp->markup)
	vtxtPrintf (blkp, " <a href='javascript:openAnchor (\"ffunc\",\"tg_pathway\")'>localize</a> in "
		    ) ;
      else
	vtxtPrintf (blkp, " to localize in ") ;
      if (nloc1 > 1)
	vtxtPrint (blkp, "various compartments (") ;
      vtxtPrint (blkp, vtxtPtr (locTxt)) ;
      if (nloc1 > 1)
	vtxtPrint (blkp, ")") ;
    }
  
  if (nitp)
    {
      if (gmp->Spc == WORM)
	{
	  vtxtDot (blkp) ;
	  if (ngp > 1)
	    vtxtPrint (blkp, "These proteins appear to ") ;
	  else
	    vtxtPrint (blkp, "This protein appears to ") ;
	  vtxtPrint (blkp, "<a href='javascript:openAnchor (\"ffunc\",\"tg_protein_interaction\")'>interact</a> with "
		     ) ;
	  if (nitp == 1)
	    vtxtPrint (blkp, "another protein (") ;
	  else
	    vtxtPrintf (blkp, "other proteins (", nitp) ;
	}
      else
	vtxtPrintf (blkp
		    , "%s%sutative  <a href='javascript:openAnchor (\"ffunc\",\"tg_protein_interaction\")'>protein interactor%s</a> %s been described ("
		    , n1 + nd + npr + npw ? ". " : ""
		    , nitp == 1 ? "A p" : "P"
		    , _multi(nitp)
		    , nitp > 1 ? "have" : "has"
		    ) ;
	
      vtxtPrint (blkp, vtxtPtr (interactsPTxt)) ;
      vtxtPrint (blkp, ")") ;
    }

   if (nitg)
    {
      if (gmp->Spc == WORM)
	{
	  vtxtDot (blkp) ;
	  vtxtPrint (blkp, "The gene <a href='javascript:openAnchor (\"ffunc\",\"tg_protein_interaction\")'>interacts</a> with ") ;
	  if (nitg > 1)
	    vtxtPrintf (blkp, "%d other genes (", nitg) ;   
	  vtxtPrint (blkp, vtxtPtr (interactsGTxt)) ;
	  if (nitg > 1)
	    vtxtPrint (blkp, ")") ;
 	}
    }
     
  if (nd + npr + npw + nitp + nitg + nloc_ph == 0) /*  && ! ac_has_tag (gmp->gene, "Cloud_gene") */
    {
      /* NO interacts  
	 disease function and localization
      */
      vtxtDot (blkp) ;
      if (ac_has_tag (gmp->gene, "Cloud_gene"))
	{
	  vtxtPrint (blkp, "No function has yet been reported to our knowledge") ; 
	  vtxtBreak (blkp) ;
	}
      else
	{
	  vtxtPrint (blkp, "No phenotype has yet been reported to our knowledge") ; 
	  vtxtPrint (blkp, ": this gene's in vivo function is yet unknown") ;   
	  vtxtBreak (blkp) ;
	}
    }
  if (gmp->Spc == MOUSE)
    {
      AC_KEYSET ksMgi = ac_objquery_keyset (gmp->gene, ">Extern ; MGI", h) ;
      int ir, jr, nMgi = ac_keyset_count (ksMgi) ;
      vtxtBreak (blkp) ;
      if (nMgi)
	{
	  AC_TABLE mgis = ac_keyset_table (ksMgi, 0, -1, 0, h) ;
	  char linkBuf[vONLINEMAXURLSIZE] ;

	  vtxtPrintf (blkp, "Please see the Jackson Laboratory Mouse Genome Database/Informatics site ") ;
	  for (ir = jr = 0 ; ir < mgis->rows ; ir++)
	    {
	      const char *ccp = ac_table_printable (mgis, ir, 0, 0) ;
	      if (!ccp || strncmp (ccp, "MGI_", 4))
		continue ;
	      if (jr++)
		vtxtPrint (blkp, ", ") ;
	      sprintf (linkBuf, MGI_LINK, ccp + 4) ;
	      gmpURL (blkp, gmp, linkBuf, ccp) ;
	    }
	  
	  vtxtPrintf (blkp, " for in depth functional annotation of this gene") ; 
	}
      else
	vtxtPrint (blkp, "This gene is not yet annotated in MGD. ") ;
    }

  else if (gmp->Spc == ARA)
    {
      AC_KEYSET ksMgi = ac_objquery_keyset (gmp->gene, ">Extern ; TAIR", h) ;
      int ir, jr, nMgi = ac_keyset_count (ksMgi) ;
      vtxtBreak (blkp) ;
      if (nMgi)
	{
	  AC_TABLE mgis = ac_keyset_table (ksMgi, 0, -1, 0, h) ;
	  char linkBuf[vONLINEMAXURLSIZE] ;

	  vtxtPrintf (blkp, "Please see the Arabidopsis TAIR site ") ;
	  for (ir = jr = 0 ; ir < mgis->rows ; ir++)
	    {
	      const char *ccp = ac_table_printable (mgis, ir, 0, 0) ;
	      if (!ccp || strncmp (ccp, "TAIR:", 5))
		continue ;
	      if (jr++)
		vtxtPrint (blkp, ", ") ;
	      sprintf (linkBuf, TAIR_LINK, ccp + 5) ;
	      gmpURL (blkp, gmp, linkBuf, ccp) ;
	    }
	  
	  vtxtPrintf (blkp, " for in depth functional annotation of this gene") ; 
	}
      else
	vtxtPrint (blkp, "This gene is not yet annotated in TAIR. ") ;
    }

  ficheGmpG26Statement (blkp, gmp, FALSE) ; 
  ficheGmpHighlyMutableStatement (blkp, gmp, FALSE) ; 
  ficheGmpRnaEditingStatement (blkp, gmp, FALSE) ; 
  ficheGmpSelenocysteineStatement (blkp, gmp, FALSE) ; 
  ficheGmpTranslational_frameshiftStatement (blkp, gmp, FALSE) ; 

  ac_free (h) ;  
  return 1 ;
} /* ficheNewGeneFunctionStatement */

/***************************************************************************************/

static int ficheNewGeneMainSupportingClonesParagraphContent (vTXT blkp, GMP *gmp, Array tissues)
{
  AC_KEYSET clones = 0 ;
  int nClones = 0, nMrna ;
  AC_HANDLE h = 0 ;
  vTXT buf = 0 ;

  if (!gmp->tg)
    return 0 ;
  
  h = ac_new_handle () ;
  /* title of the section */
  
  clones = ac_objquery_keyset (gmp->gene, "> transcribed_gene ; >Read ; NOT composite && NOT IS U* ; >cdna_clone ", h) ;
  
  if (clones && (nClones = ac_keyset_count (clones)))
    {
      buf = vtxtHandleCreate (h) ;
      if (0)
	{
	  vtxtPrintf (buf, "<a href=\"javascript:openAceViewLink ('%s', '%s', 'clones')\">"  
		      , ac_class (gmp->gene), ac_name (gmp->gene)
		      ) ;
	  vtxtPrintf (buf, "%s cDNA clone%s</a> support%s" 
		      , _isone (nClones)
		      , _multi(nClones)
		      , _verbs(nClones)
		      ) ;
	  nMrna = ac_keyset_count (ac_objquery_keyset (gmp->tg, ">mrna", h)) ;
	  if (nMrna > 1)
	    vtxtPrintf (buf, " the  <a href=\"javascript:openAnchor ('fgene','Gene_compact_diagram_4')\">%d variants</a> of", nMrna) ;
	  vtxtPrintf (buf, " gene %s", ac_name (gmp->gene)) ;
	}
      else
	{
	  vtxtPrintf (buf, "%s cDNA clone%s support "
		      , _isone (nClones)
		      , _multi(nClones)
		      , _verbs(nClones)
		      ) ;
	  if (gmp->nMrna > 1)
	    vtxtPrintf (buf, " the %d variants of", gmp->nMrna) ;
	  vtxtPrintf (buf, " gene %s", ac_name (gmp->gene)) ;
	}
      
      gmpSubSection (blkp, gmp, "allSupportingClones", vtxtPtr (buf)) ;

      vtxtPrint (blkp, "This table helps analyze the pattern of expression of the gene") ;
      if (gmp->nMrna > 1) 
	vtxtPrint (blkp, ", the tissue, cell type or disease state specificity of the alternative variants") ;
      vtxtPrint (blkp, " and to select cDNA clones suitable for your experiments") ;
      vtxtBreak (blkp) ;

      ficheNewCloneTable (blkp, gmp, clones, 'q', 0, 0, tissues) ;
    }
  
  vtxtBreak (blkp) ;
  ac_free (h) ;

  return nClones ;
} /* ficheNewGeneMainSupportingClonesParagraphContent */

/******************/

static int ficheNewGeneExpressionParagraph (vTXT blkp, GMP *gmp, Array tissues)
{
  int nClo = 0 ;
  char *ptr ;
  vTXT bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  /* expression level */
  ficheMRNAExpressionProfileParagraphContent (bfr, gmp, FALSE, TRUE) ;
  ficheGmpYkImageStatement (bfr, gmp) ;
  ficheGmpPatternStatement (bfr, gmp) ;
  if (gmp->view == 'g' && gmp->Spc == WORM)
    {
     ficheCloneLibParagraphContent (blkp, gmp) ;
     vtxtBreak (blkp) ;
    }

  ficheGmpExpressionPatternStatement (bfr, gmp) ;

  if ((ptr = vtxtPtr (bfr)))
    {
      gmpSection (blkp, gmp, "tg_expression", "Expression") ;
      vtxtPrint (blkp, ptr) ; 
    }

  /* cdna support */
  vtxtBreak (blkp) ; 
  vtxtClear (bfr) ;
  nClo = ficheNewGeneMainSupportingClonesParagraphContent (bfr, gmp, tissues) ;

  if ((ptr = vtxtPtr (bfr)))
    {
      vtxtPrint (blkp, ptr) ; 
    }
  vtxtDestroy (bfr) ;
  return nClo ;
} /* ficheNewGeneExpressionParagraph */

/***************************************************************************************/
/***************************************************************************************/

static BOOL ficheNewGeneIntronsStatement (vTXT blkp, GMP *gmp, BOOL decorate, BOOL showTitle)
{
  int iFeet, ir, isComma, numGap, numIntrons, numAltIntrons, numMrnaWithStandardIntrons ;
  int nIntron, numGTAG, numGCAG, numATAC, numFuzzy, numFuzzy_gt_ag, numFuzzy_gc_ag, numOther, nExon  ;
  int nGoodProduct ;
  AC_TABLE  gmRNA, gSplicing ; 
  const char *tag ;
  DICT *feetDict = 0 ;
  Array feetN = 0 ;
  BOOL isCloud = FALSE ;
  AC_HANDLE h ;
  vTXT blkp1 ;

  if (!gmp->tg)
    return TRUE ;

  h = ac_new_handle () ;
  blkp1 = vtxtHandleCreate (h) ;
  if (gmp->markup)
    vtxtMarkup (blkp1) ;

  /* count all types of introns */
  gmRNA = ac_tag_table (gmp->tg, "mRNA", h) ; 
  /*   gTGProduct = ac_tag_table (gmp->tg, "Product", h) ;   */
  numMrnaWithStandardIntrons = ac_keyset_count (ac_objquery_keyset (gmp->tg, ">mrna; gt_ag || gc_ag", h)) ;
  numIntrons = numGap = numGTAG = numGCAG = numATAC = numFuzzy = numFuzzy_gt_ag = numFuzzy_gc_ag = numOther = numAltIntrons = nExon = 0 ;
  gSplicing = ac_tag_table (gmp->tg, "Splicing", h) ; 
  if (gSplicing)
    {
      for (ir=0 ; ir < gSplicing->rows ; ir++)
	{
	  tag = ac_table_printable (gSplicing, ir, 2, "") ;
	  if (strstr (tag, "Gap"))
	    numGap++ ;
	  else if (strstr (tag, "xon"))
	    nExon++ ;
	  else if (strstr (tag, "tron"))
	    {
	      numIntrons++ ;
	      if (strstr (tag, "ternativ")) numAltIntrons++ ;
	      tag = ac_table_printable (gSplicing, ir, 3, "") ;
	      if (!strcmp (tag, "gt_ag")) numGTAG++ ;
	      else if (!strcmp (tag, "gc_ag")) numGCAG++ ;
	      else if (!strcmp (tag, "at_ac")) numATAC++ ;
	      else if (!strcasecmp (tag, "fuzzy_gt_ag")) { numFuzzy_gt_ag++,  numGTAG++ ;}
	      else if (!strcasecmp (tag, "fuzzy_gc_ag")) { numFuzzy_gc_ag++ ; numGCAG++ ;}
	      else if (!strcasecmp (tag, "fuzzy")) numFuzzy++ ;
	      else
		{ 
		  numOther++ ;
		  if (!feetDict)
		    {
		      feetDict = dictHandleCreate (12, h) ;
		      feetN = arrayHandleCreate (12, int, h) ;
		    }
		  dictAdd (feetDict, tag, &iFeet) ;
		  array (feetN, iFeet, int)++ ;
		}
	    }
	}
    }

  /* 5 cas : Cloud yes
                   no   hasIntrons no   hasGoodProduct n
                                   no   hasGoodProduct y
                                   yes  hasGoodProduct n
                                   yes  hasGoodProduct y
  */
  nIntron = numGTAG+numGCAG+numATAC+numFuzzy+numOther ;
  nGoodProduct = ac_keyset_count (ac_objquery_keyset (gmp->tg, ">Product; Good_product", h)) ;
  if (ac_has_tag (gmp->gene, "Cloud_gene")) /* cloud */
    {
      vtxtDot (blkp1) ;
      vtxtPrint (blkp1, "This putative gene would produce one intronless mRNA, apparently non-coding") ;
      isCloud = TRUE ;
      if (gmp->Spc == WORM)
	ficheNewGeneSlStatement (blkp, gmp) ;
    }
  else if (!nIntron  && ! nGoodProduct) /* unspliced non-coding */
    {
      vtxtDot (blkp1) ;
      vtxtPrint (blkp1, "The gene produces one intronless mRNA, apparently non-coding") ; 
      isCloud = TRUE ;
      if (gmp->Spc == WORM)
	ficheNewGeneSlStatement (blkp, gmp) ;
    }
  else if (!nIntron  && nGoodProduct) /* unspliced coding */
    {
      vtxtDot (blkp1) ; 
      vtxtPrint (blkp1, "The gene produces one intronless mRNA, putatively encoding ") ;
      vtxtPrintf(blkp1, "%s%d protein%s%s"
		 , gmp->markup ? "<a href='javascript:openAnchor (\"fmol\",\"mrnaAnnotTable\")'>" : ""
		 , nGoodProduct, nGoodProduct > 1 ? "s" : ""
		 , gmp->markup ? "</a>" : ""
		 ) ;  
      isCloud = TRUE ;
      if (gmp->Spc == WORM)
	ficheNewGeneSlStatement (blkp, gmp) ;
    }
  else if (nIntron)
    {
      vtxtBreak (blkp1) ;
      if (showTitle && gmp->markup && numMrnaWithStandardIntrons > 1)
	{
	  vtxtPrint (blkp1, "<b>Alternative mRNA variants and regulation:</b> ") ;
	}
 
      if (! decorate)
	{
	  int ok = 0 ;

	  if (gmp->view == 'm')
	    vtxtPrintf (blkp1, "The mRNA contains") ; 
	  else
	    vtxtPrintf (blkp1, "The gene contains") ; 
	  if (nIntron)
	    {
	      vTXT buf = vtxtCreate () ;
	      
	      vtxtMarkup (buf) ;
	      vtxtPrintf (buf, " %d", nIntron) ;
	      if (nIntron > 1) 
		vtxtPrint (buf, " distinct") ;
	      if (nIntron == numGTAG)
		{ ok = 1 ; vtxtPrintf (buf, " gt-ag") ; }
	      else if (nIntron == numGCAG)
		{ ok = 1 ; vtxtPrintf (buf, " gc-ag") ; }
	      if (gmp->markup)
		vtxtPrintf (buf, " intron%s</a>", _multi (nIntron)) ;
	      else
		vtxtPrintf (buf, " intron%s", _multi (nIntron)) ;
	      if (! ok)
		{
		  char *prefix = " (" ;
		  if (numGTAG)
		    { vtxtPrintf (buf, "%s%d gt-ag", prefix, numGTAG) ; prefix = ", " ; }
		  if (numGCAG)
		    { vtxtPrintf (buf, "%s%d gc-ag", prefix, numGCAG) ; prefix = ", " ; }
		  if (numATAC)
		    { vtxtPrintf (buf, "%s%d at-ac", prefix, numATAC) ; prefix = ", " ; }
		  if (numOther)
		    { vtxtPrintf (buf, "%s%d other%s", prefix, numOther, _multi(numOther)) ; prefix = ", " ; }
		  if (numFuzzy)
		    { vtxtPrintf (buf, "%s%d fuzzy", prefix, numFuzzy) ; prefix = ", " ; }
		  vtxtPrint (buf, ")") ;
		}
	      if (gmp->markup)
		{
		  if (gmp->view == 'g' && gmRNA && gmRNA->rows > 1)
		    vtxtPrintf (blkp1, "<a href=\"javascript:openAnchor ('fmol','tg_introns')\">%s</a>", vtxtPtr (buf)) ;
		  else if (gmp->view == 'g' && gmRNA &&  gmRNA->rows == 1)
		    gmpObjLinkAnchor (blkp1, gmp, gmp->mrna, "tg_introns", vtxtPtr (buf)) ;
		  else
		    gmpObjLinkAnchor (blkp1, gmp, gmp->gene, "tg_introns", vtxtPtr (buf)) ;
		}
	      else
		vtxtPrint (blkp1, vtxtPtr (buf)) ;
	      vtxtDestroy (buf) ;
	    }
	  else
	    vtxtPrintf (blkp1, " no intron") ;
	}
      else
	{
	  vtxtDot (blkp1) ;
	  vtxtPrintf (blkp1, "The gene contains") ; 
	  if (numIntrons>0)vtxtPrintf (blkp1, " %d", numIntrons) ; 
	  else vtxtPrintf (blkp1, " no") ; 
	  vtxtPrintf (blkp1, " intron%s", _multi (numIntrons)) ; 
	  
	  if (decorate && numAltIntrons > 1)
	    {
	      if (0)
		{
		  /* this sentence is unreliable */
		  if (numIntrons == numAltIntrons)
		    vtxtPrint (blkp1, ", all of which are alternative") ;
		  else
		    vtxtPrintf (blkp1, ", %d of which %s alternative", numAltIntrons, _isare (numAltIntrons)) ; 
		}
	      if (gmp->markup)
		vtxtPrint (blkp1,
			   ". Locally, in the diagram, different introns are highlighted in different colors"
			   ". Constitutive introns correspond to a common color traversing the entire mRNA set"
			   ) ;
	    }
	}
     
      if (!decorate)
	{
	  vtxtBreak (blkp) ;
	  vtxtPrint (blkp, vtxtPtr (blkp1)) ;
	}
	
      if (decorate && numGap + numGTAG+numGCAG+numATAC+numFuzzy+numOther)
	{
	  vtxtBreak (blkp) ;
	  vtxtPrintf (blkp, "Comparison to the genome sequence shows that") ; 
	  
	  isComma=0 ; 
	  if (numGTAG - numFuzzy_gt_ag)
	    {
	      vtxtPrintf (blkp, "%s %s intron%s %s well supported and follow%s the consensual [gt-ag] rule%s"
			  , isComma ? "," : ""
			  , _isone (numGTAG - numFuzzy_gt_ag) 
			  , _multi (numGTAG - numFuzzy_gt_ag)
			  , _isare (numGTAG - numFuzzy_gt_ag) 
			  , _verbs (numGTAG - numFuzzy_gt_ag)
			  , gmp->markup ? " (pink broken line in gene diagrams, at least one clone exactly matches the genome in the 16bp bordering the splice site)": ""
			  ) ; 
	      isComma=1 ;
	    }
	  if (numGCAG - numFuzzy_gc_ag)
	    {
	      vtxtPrintf (blkp, "%s %s splice%s follow%s "
			  , isComma ? "," : ""
			  , _isone(numGCAG - numFuzzy_gc_ag)
			  , _multi (numGCAG - numFuzzy_gc_ag), _verbs (numGCAG - numFuzzy_gc_ag)) ; 
	      vtxtPrintf (blkp, "the less frequent consensus [gc-ag]%s%s%s"
			  , gmp->markup ? " (" : ""
			  , numGTAG && gmp->markup ? "also" : ""
			  , gmp->markup ? " pink broken line)": ""
			  ) ; 
	      isComma=1 ; 
	    }
	  if (numATAC)
	    {
	      vtxtPrintf (blkp, "%s %s splice%s follow%s " 
			  , isComma ? "," : ""
			  , _isone (numATAC)
			  , _multi (numATAC), _verbs (numATAC)) ; 
	      vtxtPrintf (blkp, "the U12 assosiated consensus [at-ac]%s%s%s"
			  , gmp->markup ? " (" : ""
			  , numGTAG + numGCAG && gmp->markup ? "also" : ""
			  , gmp->markup ? " pink broken line)": ""
			  ) ; 
	      isComma=1 ; 
	    }
	  if (numOther)
	    {
	      int oN, oN1 = 0 ;

	      vtxtPrintf (blkp, "%s %s splice%s "
			  , isComma ? "," : ""
			  , _isone (numOther)
			  , _multi (numOther)) ; 
	      vtxtPrintf (blkp, "%s atypical with multiple support"
			  , _isare (numOther)
			  ) ; 
	      isComma=1 ; 

	      for (iFeet = 0 ; feetN && iFeet < arrayMax (feetN) ; iFeet++)
		{
		  oN = array (feetN, iFeet, int) ;
		  if (oN < 1)
		    continue ;
		  if (oN1++) vtxtPrint (blkp, ", ") ;
		  if (oN > 1) vtxtPrintf (blkp, " %d:", oN) ;
		  vtxtPrintf (blkp, " [%s]", dictName (feetDict, iFeet)) ;
		}
	       if (gmp->markup)
		vtxtPrint (blkp," (blue broken line: might corresponding to deletion in cDNA insert, polymorphism or rearrangement or sequencing error in genome, or unusual splice site)") ;
	    }
	  if (numFuzzy_gt_ag + numFuzzy_gc_ag)
	    {
	      vtxtPrintf (blkp, "%s %s splice%s " 
			    , isComma ? "," : ""
			    , _isone (numFuzzy_gt_ag)
			    , _multi (numFuzzy_gt_ag)
			    ) ; 
	      vtxtPrintf (blkp, "%s fuzzy gt-ag"
			  , _isare (numFuzzy_gt_ag)
			  ) ; 
	      vtxtPrintf (blkp, "%s %s splice%s " 
			    , isComma ? "," : ""
			    , _isone (numFuzzy_gc_ag)
			    , _multi (numFuzzy_gc_ag)
			    ) ; 
	      vtxtPrintf (blkp, "%s fuzzy gc-ag"
			  , _isare (numFuzzy_gc_ag)
			  ) ; 
	      if (gmp->markup)
		vtxtPrint (blkp, " (straight pink line, the supporting clone differs from the genome by at least one base in the 16 bp bordering the splice site)") ;
	      isComma=1 ; 
	    }
	  if (numFuzzy)
	    {
	      vtxtPrintf (blkp, "%s %s splice%s " 
			    , isComma ? "," : ""
			    , _isone (numFuzzy)
			    , _multi (numFuzzy)
			    ) ; 
	      vtxtPrintf (blkp, "%s fuzzy"   /* non standard and fuzzy */
			  , _isare (numFuzzy)
			  ) ; 
	      if (gmp->markup)
		{
		  if (numFuzzy)
		    vtxtPrint (blkp, " (straight blue line)") ;
		  else
		    vtxtPrint (blkp, " (straight line, the supporting clone differs from the genome by at least one base in the 16 bp bordering the splice site)") ;
		}
	      isComma=1 ; 
	    }
	  if (numGap)
	    {
	      vtxtDot (blkp) ;
	      vtxtPrint (blkp
			 , "Some mRNAs also contain discontinuities or sequencing gaps") ;
	      if (gmp->markup)
		vtxtPrint (blkp, " (straight gray line)") ;

	      isComma=1 ; 
	    }
	}
    }
  ac_free (h) ;
  return isCloud ;
} /* ficheNewGeneIntronsStatement */

/***************************************************************************************/

static int ficheNewGeneAltFeatureStatement (vTXT blkp, GMP *gmp, BOOL decorate, BOOL fromSummary)
{
  AC_HANDLE h = 0 ;
  AC_OBJ oGene = gmp->gene, oMrna = 0 ;
  AC_TABLE tMrna = 0, vTable = 0 ;
  int n1, n2, n3, n0 = 0, a1, ir, jr, jj ; 
  KEYSET vv = 0 ;
  BOOL ok = FALSE ;
  
  if (!ac_has_tag (gmp->gene, "Transcribed_gene"))
    return FALSE ;
  h = ac_new_handle () ;
  tMrna = ac_tag_table (gmp->tg, "mRNA", h) ;
  n1 = ac_tag_int (oGene, "nPossiblePromotors", 0) ;
  if (!n1) /* nov 2 2003, i added a direct count more accurate of nPossiblePromotors */
    n1 = ac_tag_int (oGene, "nNonOverlappingAltFirstExons", 0) ;
  n2 = ac_tag_int (oGene, "nNonOverlappingAltLastExons", 0) ;
  n3 = 0 ;
  for (ir = jj = 0 ; tMrna && ir < tMrna->rows ; ir++)
    {
      oMrna = ac_table_obj (tMrna, ir, 0, h) ;
      vTable = ac_tag_table (oMrna, "Valid3p", h) ;
      if (!vTable)
	continue ;
      if (!vv)
	vv = keySetHandleCreate (h) ;
      a1 =  ac_table_int (tMrna, ir, 0, 0) ;
      for (jr = 0 ; vTable && jr < vTable->rows ; jr++)
	{
	  keySet (vv, jj++) = a1 + ac_table_int (vTable, jr, 0, 0) ;
	}
    }

  if (vv)
    {
      keySetCompress (vv) ;
      n3 = keySetMax (vv) ;
    }
    
  if (n1 > 1)
    {
      if (!n0++) vtxtDot (blkp) ;
      vtxtPrintf (blkp
		  , "There are %d probable <a href=\"javascript:openAnchor ('fmol','mRNA_sequences') \">alternative promotors</a>"
		  , n1) ;
    }
  
  if (n2 > 1)
    {
      if (!n0++) vtxtDot (blkp) ;
      vtxtPrintf (blkp
		  , "%s %d  non overlapping alternative last exons"
		  , n1 > 1 ? (n3 > 1 ? "," : " and") : "There are"
		  , n2) ;
    }

  if (n3 > 1)
    {
      if (!n0++) vtxtDot (blkp) ;
      vtxtPrintf (blkp
		  , "%s %d  validated  <a href=\"javascript:openAnchor ('fmol','mRNAs') \">alternative polyadenylation sites</a>"
		  , n1 > 1 || n2 > 1 ? " and" : "There are"
		  , n3) ;
    }

  if (n1 > 1 || n2 > 1 || n3 > 1) 
    { 
      if (fromSummary)
	vtxtPrintf (blkp, " (see the <a href=\"javascript:openAnchor ('fgene', 'Gene_compact_diagram_1')\">diagram</a>)");
      else
	vtxtPrintf (blkp, " (see the <a href=\"javascript:openAnchor ('fmol', 'Gene_compact_diagram_2')\">diagram</a>)");
      ok = TRUE ;
    }
  
  if (gmp->tg &&
      (
       ac_has_tag (oGene, "n5pSkippedCentralExons") ||
       ac_has_tag (oGene, "n5pSkippedCentralExons") ||
       ac_has_tag (oGene, "n3pSkippedCentralExons") ||
       ac_has_tag (oGene, "nsce") ||
       ac_has_tag (oGene, "nSkippedCentralExons") ||
       ac_has_tag (oGene, "nOverlappingCentralExons") ||
       ac_has_tag (oGene, "nRetainedIntrons")
       )
      )
    {
       if (!n0++) vtxtDot (blkp) ;
       else vtxtDot (blkp) ;
       vtxtPrintf (blkp
		   , "The mRNAs appear to differ "
		   ) ;
    }
  
  n1 = 0 ; 
  if (ac_has_tag (oGene, "n5pSkippedCentralExons")) 
    {
      if (n1++)vtxtPrintf (blkp, ", ") ;
      else vtxtPrintf (blkp, " by ") ;
      vtxtPrintf (blkp, "truncation of the 5' end") ; /* N-terminus */
    }
  
  if (ac_has_tag (oGene, "n3pSkippedCentralExons"))
    {
      if (n1++)vtxtPrintf (blkp, ", ") ;
      else vtxtPrintf (blkp, " by ") ;
      vtxtPrintf (blkp, "truncation of the 3' end") ; /* C-terminus */
    }
  
  n2 = 0 ;
  if ((n2 = ac_tag_int (oGene, "nSkippedCentralExons", 0)) ||
      ac_has_tag (oGene, "nsce"))
    {
      if (n1++)vtxtPrintf (blkp, ", ") ;
      else vtxtPrintf (blkp, " by ") ;
      if (n2 < 2)
	vtxtPrintf (blkp, "presence or absence of a <a href=\"javascript:openAnchor ('fmol','tg_introns')\">cassette exon</a>") ;
      else
	vtxtPrintf (blkp, "presence or absence of %d <a href=\"javascript:openAnchor ('fmol','tg_introns')\">cassette exons</a>", n2) ;
    }
  
  if (ac_has_tag (oGene, "nOverlappingCentralExons"))
    {
      if (n1++) vtxtPrintf (blkp, ", ") ;
      else vtxtPrintf (blkp, " by ") ;
      vtxtPrintf (blkp, "overlapping exons with different boundaries") ;
    }
  
  if (ac_has_tag (oGene, "nRetainedIntrons"))
    {
      int ni = ac_tag_int (oGene, "nRetainedIntrons", 0) ;
      if (n1++)vtxtPrintf (blkp, ", ") ;
      
      vtxtPrintf (blkp, " splicing versus retention of %s intron%s", isOne (ni), _multi(ni)) ;
    }
  if (gmp->Spc == WORM && gmp->tg)
    {
      AC_KEYSET vks = ac_objquery_keyset (gmp->tg, "cdna_clone  && ! (cdna_clone != OST*)", 0) ;
      int nv = ac_keyset_count (vks) ;
      if (nv)
	{
	  vtxtDot (blkp) ;
	  vtxtPrintf (blkp, "The gene may be incomplete, and the last 20 bp at each end may not exist in vivo, because ORFeome clones are PCR designed to match predictions") ;
	}
      ac_free (vks) ;  
  }
  if (n1) ok = TRUE ;

  ac_free (h) ;
  return ok ;
} /* ficheNewGeneAltFeatureStatement */

/******************/
#ifdef JUNK
 removed 2007_03_04 danielle
static void ficheNewGeneProductCompletionStatement (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = 0 ;
  int n = 0 ;
  int nC, nC3, nC5, nC35 ;
  AC_KEYSET products = 0 ;

  if (!gmp || !gmp->tg)
    return ;
  h = ac_new_handle () ;
  products = ac_objquery_keyset (gmp->tg, ">mrna ; gt_ag || gc_ag ; >Product ; Best_product", h) ;
  nC = ac_keyset_count (ac_ksquery_keyset (products, "COOH_Complete && NH2_Complete", h)) ;
  nC5 = ac_keyset_count (ac_ksquery_keyset (products, "COOH_Complete && !NH2_Complete", h)) ;
  nC3 = ac_keyset_count (ac_ksquery_keyset (products, "!COOH_Complete && NH2_Complete", h)) ;
  nC35 = ac_keyset_count (ac_ksquery_keyset (products, "!COOH_Complete && !NH2_Complete", h)) ;

  if (nC > 0)
    {
      if (!n++) vtxtDot (blkp) ; 
      vtxtPrintf (blkp, "%d transcript%s appear%s to be complete at both the 5' and 3' end"
		  , nC, nC > 1 ? "s" : "",  nC > 1 ? "" : "s"
		  ) ;
    }
  if (nC5 > 0)
    {
      if (!n++) 
	{ 
	  vtxtDot (blkp) ; 
	  vtxtPrintf (blkp, "%d transcript%s seem%s incomplete at the 5' end"
		      , nC5, nC5 > 1 ? "s" : "",  nC5 > 1 ? "" : "s"
		      ) ;
	}
      else
	{
	  vtxtComma (blkp) ;
	  vtxtPrintf (blkp, "%d  seem%s incomplete at the 5' end"
		      , nC5, nC5 > 1 ? "" : "s"
		      ) ;
	}
    }
  if (nC3 > 0)
    {
      if (!n++) 
	{ 
	  vtxtDot (blkp) ; 
	  vtxtPrintf (blkp, "%d transcript%s %s open at the 3' end"
		      , nC3, nC3 > 1 ? "s" : "",  nC3 > 1 ? "are" : "is"
		      ) ;
	}
      else
	{
	  vtxtComma (blkp) ;
	  vtxtPrintf (blkp, "%d %s open at the 3' end"
		      , nC3, nC3 > 1 ? "are" : "is"
		      ) ;
	}
    }
  if (nC35 > 0)
    {
      if (!n++) 
	{ 
	  vtxtDot (blkp) ; 
	  vtxtPrintf (blkp, "%d transcript%s seem%s incomplete at both ends"
		      , nC35, nC35 > 1 ? "s" : "",  nC35 > 1 ? "" : "s"
		      ) ;
	}
      else
	{
	  vtxtComma (blkp) ;
	  vtxtPrintf (blkp, "%d seem%s incomplete at both ends"
		      , nC35, nC35 > 1 ? "" : "s"
		      ) ;
	}
    }
  ac_free (h) ;
} /* ficheNewGeneProductCompletionStatement */
#endif

/******************/

static void ficheNewGenePleaseQuote (vTXT blkp, GMP *gmp)
{
  char oldStyle = gmp->markup ;
  
  if (gmp->style == 'r') gmp->markup = TRUE ;
  if (vtxtPtr (blkp))
    vtxtEmptyLine (blkp, 1) ;
  vtxtPrint (blkp,"Please quote: ") ;

  if (gmp->markup)
    vtxtPrint (blkp,"<span class='quote_aceview'>") ;
 
  if (1) gmpURL (blkp, gmp, "http://genomebiology.com/2006/7/S1/S12", "AceView: a comprehensive cDNA-supported gene and transcripts annotation, Genome Biology 2006, 7(Suppl 1):S12") ;
  else gmpURL (blkp, gmp, "Papers/2006ThierryMiegDJAceView.pdf", "AceView: a comprehensive cDNA-supported gene and transcripts annotation, Genome Biology 2006, 7(Suppl 1):S12") ;
  
  if (gmp->markup)
    vtxtPrint (blkp,"</span>") ;
  gmp->markup = oldStyle ;
} /* ficheNewGene */

/******************/

void ficheNewGeneAltFeatureChapter (vTXT blkp, GMP *gmp)
{
  char *ptr ;
  BOOL isCloud ;
  vTXT bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  isCloud = ficheNewGeneIntronsStatement (bfr, gmp, FALSE, FALSE) ;
  if (!isCloud) ficheNewGeneAltVariantStatement (bfr, gmp, TRUE) ;
  ficheNewGeneAltFeatureStatement (bfr, gmp, TRUE, FALSE) ;
  if (gmp->tg) ficheTGAntisensParagraphContent (bfr, gmp, FALSE) ; 
  if (gmp->tg) ficheNewGeneNmdStatement (bfr, gmp, TRUE) ;
  if (!isCloud) ficheNewGeneProteinStatement (bfr, gmp, TRUE) ;
  if (gmp->tg)
    {
      ficheNewGenePfamPsortStatement (bfr, gmp, isCloud, FALSE) ;
      if (!isCloud) ficheNewGeneNonGoodVariantStatement (bfr, gmp) ;
      ficheNewGeneKozakStatement (bfr, gmp) ;
    }
  if (gmp->tg) ficheNewGeneComplexLocusStatement (bfr, gmp) ;
 
 if ((ptr = vtxtPtr (bfr)))
    {
      gmpChapter (blkp, gmp, "*tg_alt_features", "Alternative mRNA variants and regulation") ;
      
      vtxtPrint (blkp, ptr) ; 
      gmpChapterClose (blkp, gmp, "tg_alt_features", TRUE) ;
    }
  vtxtDestroy (bfr) ;
} /* ficheNewGeneAltFeatureChapter */

/***************************************************************************************/

static void ficheNewGeneRegulationParagraph (vTXT blkp, GMP *gmp)
{
  char *ptr ;
  vTXT bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  if (gmp->markup && gmp->tg)
    {
      ficheTGOperonParagraphContent (bfr, gmp, ficheOperonDistance, 1, 1, 0, TRUE) ; 
      if (0) /* already explained in the abstract */
	ficheTGAntisensParagraphContent (bfr, gmp, TRUE) ; 
    }

  ficheGmpG26Statement (bfr, gmp, TRUE) ; 
  ficheGmpHighlyMutableStatement (bfr, gmp, TRUE) ; 
  ficheGmpRnaEditingStatement (bfr, gmp, TRUE) ; 
  ficheGmpSelenocysteineStatement (bfr, gmp, TRUE) ; 
  ficheGmpTranslational_frameshiftStatement (bfr, gmp, TRUE) ; 

  if ((ptr = vtxtPtr (bfr)))
    {
      gmpSection (blkp, gmp, "tg_regulation", "Regulation") ;
      vtxtPrint (blkp, ptr) ; 
    }
  vtxtDestroy (bfr) ;
} /* ficheNewGeneRegulationParagraph */

/***************************************************************************************/

void ficheNewGeneSummaryChapter (vTXT blkp, GMP *gmp)
{
  gmpChapter (blkp, gmp, "*SUMMARY", "SUMMARY") ;
  ficheNewGeneSummaryParagraph (blkp, gmp) ;  /* shed + refseq summary + proteome summary */
  ficheNewGeneMappingLinksAliasesSubSection (blkp, gmp) ;
  ficheNewGeneAceKogSubSection (blkp, gmp) ;
    gmpChapterClose (blkp, gmp, "SUMMARY", TRUE) ;
} /* ficheNewGeneSummaryChapter */

/***************************************************************************************/
/***************************************************************************************/
/* FLASH chapter */

void ficheListAndMarkUpMrnas (vTXT blkp, GMP *gmp, char type, BOOL fromGene)
{
  AC_OBJ oMrna ;
  AC_TABLE gmRNA ; 
  const char *mrna ;
  int ir, nn = 0 ;
  AC_HANDLE h = ac_new_handle () ;

  vtxtBreak (blkp) ;
  if (gmp->tg)
    {
      vtxtPrintf (blkp, "You may access from the top of the page or from here the ") ;
      gmpObjLink (blkp, gmp, gmp->gene, "gene") ;
      vtxtPrint (blkp, " and other <span class='ace_summary'>annotated mRNAs</span>: ") ;

      if ((gmRNA = ac_tag_table (gmp->tg, "mRNA", h)))
	{
	  for (ir = nn =0 ; ir < gmRNA->rows ; ir++)
	    {
	      oMrna = ac_table_obj (gmRNA, ir, 0, h) ; 
	      if (!fromGene && ac_obj_equal (gmp->mrna, oMrna))
		continue ;
	      mrna = ac_name (oMrna) ;
	      if (nn++) vtxtPrint (blkp, ", ") ;
	      gmpObjLink (blkp, gmp, oMrna, gtMrnaSuffix (ac_name (gmp->gene), mrna, h)) ;
	    }
	}
    }
  ac_free (h) ;
  return ;
} /* ficheListAndMarkUpMrnas */

/***************************************************************************************/

void ficheGenomeSummaryChapter (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = gmp->gene ? ac_tag_table (gmp->gene, "IntMap", h) : 0 ;
  int a1, a2 ;

  gmpChapter (blkp, gmp, "#*Genome_summary", "Mapping") ;
  if (tbl && tbl->cols >= 3)
    {
      AC_OBJ chrom = ac_table_obj (tbl, 0, 0, h) ;
      a1 = ac_table_int (tbl, 0, 1, 0) ;
      a2 = ac_table_int (tbl, 0, 2, 0) ;
      vtxtPrintf (blkp, "Gene %s maps on %s strand of %s%s from base %d to %d\n"
		  , ac_name (gmp->gene)  
		  , a1 < a2 ? "plus" : "minus"
		  , strncasecmp ("chr", ac_name (chrom), 3) ? "chromosome " : ""
		  , ac_name (chrom)
		  , a1
		  , a2 
		  ) ;
    }
  gmpChapterClose (blkp, gmp,  "Genome_summary", TRUE) ;
  ac_free (h) ;
  return ;
} /* ficheGenomeSummaryChapter */

/***************************************************************************************/

void ficheGenomePlotChapter (vTXT blkp, AC_DB db, GMP *gmp, BOOL isBig)
{
  AC_HANDLE h = ac_new_handle () ;
  int len = 0 ;
  char *qq, *cr ;
  AC_TABLE tbl = gmp->gene ? ac_tag_table (gmp->gene, "IntMap", h) : 0 ;

  /* obsolete, this was triggering the download 
 char *cp = messprintf ("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href=\"javascript:openAceViewAction('gene','%s','vgene&v=2&B')\"><font color='red'>Complete gene on genome diagram:</font></a>", ac_name(gmp->gene)) ; 
  vtxtPrintf (blkp, "%s", cp) ;
  */

  if (!isBig)
    gmpChapter (blkp, gmp, "#*Small_Genome_plot", "Show just this gene") ;
  else
    gmpChapter (blkp, gmp, "#*Large_Genome_plot", "Show the genomic regiont") ;

  vtxtPrintf (blkp, "This plot is new. It uses HTML5/SVG, if it does not display correctly or if the mouse hover bubbles are not nice please report the problems to mieg@ncbi.nlm.nih.gov, thank you.") ;



  if (HTML5 && tbl && tbl->cols >= 3)
    {
      AC_OBJ chrom = ac_table_obj (tbl, 0, 0, h) ;
      int a1 = ac_table_int (tbl, 0, 1, 0) ;
      int a2 = ac_table_int (tbl, 0, 2, 0) ;
      int da = 300 ;

      da = (a2 - a1)/20 ; 
      if (isBig)
	da = 3 * 20 * da ;
      a1 -= da ; a2 += da ; 

      qq = hprintf (h,   "GIF ; dimensions 3000 500 ; seqget %s -coords %d %d -view av_tg_whole ; seqdisplay ; svgdump -",  ac_protect (ac_name (chrom), h), a1, a2) ;
      cr = (char *)ac_command (db, qq, &len, h) ;  
      if (1)  cr = strchr (cr,'<') ;
      vtxtPrint (blkp, cr) ;
    }	 

    if (!isBig)
      gmpChapterClose (blkp, gmp,  "Small_Genome_plot", TRUE) ;
    else
      gmpChapterClose (blkp, gmp,  "Large_Genome_plot", TRUE) ;
  ac_free (h) ;
  return ;
} /* ficheGenomePlotChapter */

/***************************************************************************************/

void ficheNewGeneGenomeDiagramChapter (vTXT blkp, GMP *gmp)
{
  char *cp = messprintf ("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href=\"javascript:openAceViewAction('gene','%s','vgene&v=2&B')\"><font color='red'>Complete gene on genome diagram:</font></a>", ac_name(gmp->gene)) ; 
  gmpChapter (blkp, gmp, "#*Gene_genome_diagram", cp) ;

  vtxtPrint (blkp, " Please choose between the ") ;
  vtxtPrintf (blkp, "<a href=\"javascript:openAceViewAction('gene','%s','vgene&v=2&B')\">zoomable GIF version</a>."
	      , ac_name(gmp->gene)) ;
  vtxtPrint (blkp, ", and the ") ;
  vtxtPrintf (blkp, "<a href=\"javascript:openAceViewAction('gene','%s','vgene&v=2&S')\">HTML5/SVG version</a>"
		  , ac_name(gmp->gene)) ;
  vtxtBreak (blkp) ;
  vtxtPrint (blkp, "This diagram shows in true scale the gene on the genome, the mRNAs and the cDNA clones.") ;

  gmpChapterClose (blkp, gmp, "Gene_genome_diagram", TRUE) ;
} /* ficheNewGeneGenomeDiagramChapter */

/***************************************************************************************/

void ficheNewGeneMrnaDiagramChapter (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl ;
  int ir ;

  if (gmp->tg &&
      (tbl = ac_tag_table (gmp->tg, "mRNA", h)))
    {
      gmpChapter (blkp, gmp, "#*Gene_mrna_diagrams", "Annotated mRNA diagrams") ;
      
      vtxtPrint (blkp, "The mRNAs diagrams with the aligned cDNA sequence accessions"
		 " and their mismatches are available in the mRNA pages accessible from the tab at the top"
		 " of the page, or here: "
		 ) ; 
      vtxtBreak (blkp) ;
      vtxtPrint (blkp, " In Flash: ") ;
      for (ir = 0 ; ir < tbl->rows ; ir++)
	vtxtPrintf (blkp, "<a href=\"javascript:openAceViewAction('mrna','%s','vmrna&v=2&S')\">%s%s</a>"
		    , ac_table_printable (tbl, ir, 0, "")
		    , ir ? ", " : ""
		    , gtMrnaSuffix (ac_name (gmp->gene), ac_table_printable (tbl, ir, 0, ""), h)
		    ) ;
      vtxtBreak (blkp) ;
      vtxtPrint (blkp, " or in GIF: ") ;
      for (ir = 0 ; ir < tbl->rows ; ir++)
	vtxtPrintf (blkp, "<a href=\"javascript:openAceViewAction('mrna','%s','vmrna&v=2&B')\">%s%s</a>"
		    , ac_table_printable (tbl, ir, 0, "")
		    , ir ? ", " : ""
		    , gtMrnaSuffix (ac_name (gmp->gene), ac_table_printable (tbl, ir, 0, ""), h)
		    ) ;

      gmpChapterClose (blkp, gmp, "Gene_mrna_diagrams", TRUE) ;
    }
  ac_free (h) ;
  return ;
} /* ficheNewGeneMrnaDiagramChapter */

/***************************************************************************************/

void ficheNewGeneCompactDiagramChapter (vTXT blkp, GMP *gmp, int pass)
{
  AC_HANDLE h = ac_new_handle () ;
  int nMrna = gmp->gene ? ac_keyset_count (ac_objquery_keyset (gmp->gene, ">transcribed_gene ; >mrna", h)) : 0 ;
  int height = 70 + 25 * nMrna ;

  /* pass is set to distinguish the header in the cookie list */
  gmpChapter (blkp, gmp, hprintf (h, "#*Gene_compact_diagram_%d", pass), "Compact gene diagram") ;

  if (HTML5)
    {
      int len = 0 ;
      AC_DB db = ac_open_db ("local", 0) ;
      char *qq = hprintf (h, "view  -c Gene -n %s  -v DtHSEQ -svg", ac_protect (ac_name (gmp->gene), h)) ;
      vtxtPrint (blkp, (char *)ac_command (db, qq, &len, h)) ;  
    }
  else
    {
      vtxtPrintf (blkp,"\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
		  "  <!--\n" 
		  "  openAceViewElasticImage (\'%s\',\'%s\',\'%s\', %d, %d, 100) ; \n"
		  "  //-->\n"
		  "</script>\n" 
		  , "gene", ac_name (gmp->gene), "hmrna", 600, height
		  ) ;
    }  
  vtxtBreak (blkp) ;
  
  vtxtPrintf (blkp, "Alternative mRNAs are shown aligned from 5' to 3' on a virtual genome where introns have been shrunk to a minimal length. Exon size is proportional to length, intron height reflects the number of cDNAs supporting each intron, the small numbers show the support of the introns in deep sequencing (with details in mouse-over) . Introns of the same color are identical, of different colors are different. 'Good proteins' are pink, partial or not-good proteins are yellow, uORFs are green. 5' cap or3' poly A flags show completeness of the transcript.") ;

  gmpCaption (blkp, gmp, "Gene_compact_diagram_caption", "Read more...") ;

  vtxtPrintf (blkp, "Mouse over the ending of each transcript gives tissues from which the supporting cDNAs were extracted. Details on tissue of origin for each intron and exon is available from the <a href=\"javascript:openAnchor ('fmol','tg_introns')\">intron and exons table</a>.") ;
   vtxtBreak (blkp) ;
  
  vtxtPrintf (blkp, "Click on any transcript to open the specific mRNA page, to see the exact cDNA clone support and eventual SNPs and to get details on tissues, sequences, mRNA and protein annotations. Proteins supported by a single continuous cDNA sequence lead to underlining the name/ending of the variant. Names not underlined result from cDNA concatenation in the coding region and should be experimentally checked.") ;
   vtxtBreak (blkp) ;
  
   if (1) /* caption */
    {
      vtxtPrint (blkp,
		 "<font color=#007f7f >Introns</font> are depicted by broken lines; the height of the top of each "
		 "intron reflects the relative number of clones supporting this intron.  "
		 "<font color=#f700f7>]^[ A pink broken line</font> "
		 "denotes an intron with standard boundaries (gt-ag or gc-ag) "
		 "that is exactly supported (i.e. a cDNA sequence exactly matches the "
		 "genome over 16 bp, 8 on both sides of the intron).  "
		 "<font color=#f700f7>]</font> "
		 "<font color=#0000ff>^</font> "
		 "<font color=#f700f7>]</font> "
		 "<font color=#0000ff> A blue broken line</font> "
		 "denotes non-standard introns, exactly supported, but with non-standard "
		 "at-ac or any other boundaries.  "
		 "<font color=#f700f7>]-[ Pink</font> "
		 "and  "
		 "<font color=#f700f7>]</font> "
		 "<font color=#0000ff>-</font> "
		 "<font color=#f700f7>]</font> "
		 "<font color=#0000ff>blue</font> "
		 " straight lines "
		 "represent 'fuzzy' introns of the standard and non-standard types "
		 "respectively, those introns do not follow the 16 bp rule. Black straight "
		 "lines ]-[denote gaps in the alignments.  "
		 ) ;
      vtxtBreak (blkp) ;
      vtxtPrint (blkp,
		 "<font color=#007f7f >Exons:</font> "
		 "Wide filled pink areas represent putative protein coding regions, "
		 "narrow empty pink boxes represent the 5'UTR (on the left) and 3' UTR (on "
		 "the right). Flags identify validated endings: cap site on the 5' side, "
		 "polyadenylation site on the 3' side. Filled flags correspond to frequent "
		 "events while empty flags have lesser supporting cDNAs (yet all are "
		 "validated); at the 3' side, black flags are associated to the main "
		 "AATAAA signal,  "
		 "<font color=#0000ff>blue flags</font>  "
		 "to any single letter variant of the main "
		 ) ;
      vtxtPrint (blkp, ". More explanations are given in the ") ;
      gmpURL (blkp, gmp, "HelpGene.html#HelpGeneDiagram", "gene help file") ;
  
  
    }
  gmpChapterClose (blkp, gmp,  hprintf (h, "*Gene_compact_diagram_%d", pass), TRUE) ;
  ac_free (h) ;
} /* ficheNewGeneCompactDiagramChapter */

/***************************************************************************************/

void ficheNewGeneWiggleDiagramChapter (vTXT blkp, GMP *gmp, int pass)
{
  AC_HANDLE h = ac_new_handle () ;
  int height = 700 ;

  /* pass is set to distinguish the header in the cookie list */
  gmpChapter (blkp, gmp, hprintf (h, "#*Gene_wiggle_diagram_%d", pass), "RNA-seq expression profiles") ;

  if (1) /* end of the diagram proper */
    {
      vtxtPrintf (blkp,"\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
		  "  <!--\n" 
		  "  openAceViewElasticImage (\'%s\',\'%s\',\'%s\', %d ,%d, %d) ; \n"
		  "  //-->\n"
		  "</script>\n" 
		  , "gene", ac_name (gmp->gene), "wiggle", 1000, height, 100
		  ) ;
    }  
  
  vtxtPrint (blkp, "RNA-seq tags were aligned on the genome, and their density is plotted at a resolution of 5 bp. Different colors show different stages and types of RNA.") ;
  
  gmpCaption (blkp, gmp, "Gene_wiggle_diagram_caption", "Read more...") ;

  vtxtPrint (blkp, 
	     "This display integrates the short RNA-Seq tags from transcriptome projects from the Baillie, Kim and Piano laboratories, and short RNA sequences from Waterston, Fraser and Mello available in NCBI SRA or GEO and not embargoed. Thanks to all. Those RNA-Seq data were used to confirm the previous cDNA supported transcriptome and to discover new exons, new introns and new genes. In particular, close to 24,000 addition sites of trans-spliced leaders are now annotated. A main surprise is the frequent occurrence of transspliced leaders inside genes, at the acceptor site of internal exons, following introns longer than the typical 40 to 55 bp.") ;
  
 vtxtBreak (blkp) ;
 vtxtPrint (blkp, "Expression profile: Tags were aligned on the transcriptome and on the genome and their density on the genome is plotted at a resolution of 5 bp. The gene of interest is in the center of the genomic profile; the diagram shows about twice the size of the gene and at least 10 kb worth of genome. As can be seen in the top left legend, different colors show different stages and types of preparation. PolyA selected RNAs are in the green to blue register, from light green in embryos to dark green in L2, then light blue in L3 to dark blue in adults, dauers are violet and mixed stages gray.") ;
 vtxtBreak (blkp) ;
 vtxtPrint (blkp, "Small RNAs, isolated without polyA selection and sequenced with a strand specific protocol , are shown in light pink for the + strand and in light orange for the - strand (try for example <a href=\"javascript:openAceViewAction('gene','mir-39','fiche')\">mir-39</a> or <a href=\"javascript:openAceViewAction('gene','mir-*','fiche')\">mir-*</a> ). 26G (and soon 22G) endo siRNA are shown in darker colors, to reflect the recent work of the Kim and Mello groups defining new functions for endo siRNAs.") ; 
 vtxtBreak (blkp) ;
 vtxtPrint (blkp, "x and y scales: The position in the genome is given in chromosome coordinates (WS190) in the horizontal scale bar, the chromosome is the first letter/number in the top green box (e.g. region 1B30 means chromosome 1, second tile), the exact coordinate of the center of the window is in the top yellow box. For now this scale is fixed but flash images have a right-click menu to zoom in and look in more detail at the region of interest. The y scale, currently fixed at 200, shows the number of tags at each position and each stage. ") ;
 vtxtBreak (blkp) ;

 vtxtPrint (blkp, "Relationship between expression profile and annotated genes: In blue along the coordinate scale is a projection of the (AceView) genes summarizing the short and long cDNA sequences. All variants and elements are collapsed into a single drawing: on top are genes on the + strand, transcribed from left to right, in the bottom are genes on the - strand, transcribed from right to left. Clicking on a blue gene brings you to that gene's page, in effect allowing you to move along the chromosome and examine the expression profile along the genome. The projections of the Wormbase CDS models are shown in gray, for reference.") ;
 vtxtBreak (blkp) ;
 vtxtPrint (blkp, "Troubleshooting: Some genes such as tRNAs, or models with no RNA sequence support, or even bugs on our side will give a less informative view (e.g. try sup-5). To see the expression profile in the region, please scroll down to the 'neighbors diagram'. Your gene is in the center. Click on any pink gene you like (here for example ach-1), and you will now see the local profile (the small RNA sup-5 tRNA has light orange tags).") ; 
 vtxtBreak (blkp) ;
 vtxtPrint (blkp, "Other displays: the candidate new genes, transcript models and their alternative element are displayed in the 'compact gene' and 'gene neighbors' diagrams below. Click on any object to get more details. The ") ;

 vtxtPrintf (blkp, "<a href=\"javascript:openAceViewAction('gene','%s','vgene&v=2&S')\">gene on genome diagram</a>", ac_name(gmp->gene)) ;

 if (gmp->mrna)
   {
     vtxtPrintf (blkp, " and the ") ;
     vtxtPrintf (blkp, "<a href=\"javascript:openAceViewAction('mrna','%s','fiche')\">RNA diagram</a>", ac_name(gmp->mrna)) ;
   }

 vtxtPrintf (blkp, " display molecular summaries with cDNA evidence in basepair scale.  In the mRNA and gene diagrams, elements of the transcriptome supported by the new deep transcriptome data are shown highlighted in green. Detailed text and table descriptions can be found in the various pages (tabs on top, or links from the Gene summary paragraph). ") ;
 
 vtxtBreak (blkp) ;
 vtxtPrint (blkp, "This display is still experimental, thank you for any comments or wishes; please email us ") ;
 gmpURL (blkp, gmp, messprintf ("mailto:mieg@ncbi.nlm.nih.gov?subject=wiggle%s", ac_name(gmp->gene)), "mieg@ncbi.nlm.nih.gov") ;
  vtxtPrint (blkp, " also if you'd like to share some transcriptome data ") ;


 gmpChapterClose (blkp, gmp,  hprintf (h, "*Gene_wiggle_diagram_%d", pass), TRUE) ;
 ac_free (h) ;
} /* ficheNewGeneWiggleDiagramChapter */

/***************************************************************************************/

void ficheNewGeneExpressionProfileChapter (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ; 
  AC_TABLE tt = 0 ;
  
  if (gmp->Spc == HUMAN && gmp->gene  && ! ac_has_tag (gmp->gene, "Cloud_gene"))
    {
      gmpChapter (blkp, gmp, "#*RNA_seq", "RNA_seq discoveries") ;
       

      gmpSection (blkp, gmp, "Gene_expression_profile_in_primates", "Expression/conservation in primates tissues evaluated by cross-mapping to human.") ;

      vtxtPrintf (blkp, "<div class='%s' id='geneExpProfile'>\n", "shown") ;
      if (HTML5)
	{
	  int len = 0 ;
	  AC_DB db = ac_open_db ("local", 0) ;
	  char *qq = hprintf (h, "view  -c Gene -n %s  -v DtGeneExp -svg", ac_protect (ac_name (gmp->gene), h)) ;
	  vtxtPrint (blkp, (char *)ac_command (db, qq, &len, h)) ;  
	}
      else
	{
	  vtxtPrintf (blkp,"\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
		      "  <!--\n" 
		      "  openAceViewElasticImage (\'%s\',\'%s\',\'%s\', 600, 150, 100) ; \n"
		      "  //-->\n"
		      "</script>\n" 
		      , "gene", ac_name (gmp->gene)
		      , "gxp"
		      ) ;
	}

      vtxtBreak (blkp) ;
      vtxtBreak (blkp) ;

      vtxtPrintf (blkp, "RNA-seq gene expression profile across 16 selected tissues from the Non-Human Primates Reference Transcriptome  Resource (link to ") ;
      gmpURL (blkp, gmp, "http://nhprtr.org","NHPRTR project") ; 
      vtxtPrintf (blkp, ").") ;
      vtxtBreak (blkp) ;
      /* BROWN 168 40 0 -> #a82800
       * DARKRED 204 0 51 -> #cc0033
       * DARKGREEN 0 153 153 -> #009999
       * BLUE4 130 144 246 -> #8290f6
       * BLUE5  94 112 244 -> #5e70f4
       * GREEN6 13 205 141 -> #0dcd8d
       */
      vtxtPrintf (blkp, "- Primates: <font color='blue'>Apes</font> (<font color='blue'>HUM</font>: Human (Illumina BodyMap 2), <font color='blue'>CHP</font>: Chimpanzee), <font color='#cc0033'>Old World monkeys</font> (<font color='#cc0033'>PTM</font>: Pig-Tailed  Macaque, <font color='#cc0033'>JMI</font> Japanese Macaque, <font color='#cc0033'>RMI</font> Rhesus  Macaque Indian, <font color='#cc0033'>RMC</font> Rhesus  Macaque Chinese, <font color='#cc0033'>CMM</font> Cynomolgus  Macaque Mauritian, <font color='#cc0033'>CMC</font> Cynomolgus  Macaque Chinese, <font color='#a82800'>BAB</font> Olive Baboon, <font color='#a82800'>SMY</font> Sooty Mangabey); <font color='#009999'>New World  monkeys</font> (<font color='#009999'\">MST</font> common Marmoset, <font color='#009999'\">SQM</font> Squirrel  Monkey, <font color='#009999'>OWL</font> Owl Monkey); and <font color='orange'>Lemurs</font> (<font color='orange'>MLM</font> Mouse Lemur, <font color='orange'>RTL</font> Ring-Tailed  Lemur).") ;
      vtxtBreak (blkp) ;
      vtxtPrintf (blkp, "- The level for significantly expressed genes is color coded in 8 equal sized bins (light to  dark green). Light gray is for weak not-accurately measured expression (2 to 8 reads  above intergenic background); dark gray for no expression or no sequence conservation (0 read in gene)") ; 
      vtxtPrintf (blkp, ". The plot to the right shows the distribution of measured expression values  in all tissues  for <font color='#5e70f4'>all genes (blue)</font> and for  <font color='#0dcd8d'>this gene (green)</font>, in Magic index = log<sub>2</sub>(1000 sFPKM)") ;

      vtxtBreak (blkp) ;
      if ( (tt = ac_tag_table (gmp->gene, "IntMap", h)) && tt->cols >= 3)
	{
	  const char *nnam = ac_table_printable (tt, 0, 0, "") ;
          int a1 =  ac_table_int (tt, 0, 1, 0), a2 = ac_table_int (tt, 0, 2, 0) ;

	  if (a1 > a2) { int a0 = a1 ; a1 = a2 ; a2 = a0 ;}
	  if (a2 - a1 < 10000) { int a0 = (a1 + a2)/2 ; a1 = a0 - 5000 ; a2 = a0 + 5000 ; } 
	  a1 -= 5000 ; a2 += 5000 ;

	  /* add   target=\"_top\" in the <a   > to be outside aceview ! */

	  vtxtPrintf (blkp, "You may also  examine the strand-specific <a  href=\'http://www.genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=mieg&hgS_otherUserSessionName=primates_NB&position=chr%s:%d-%d\' target=\"_top\" >genome coverage plots&nbsp;on the experimental AceView/Magic hub at UCSC</a>, by tissue or by species" 
		      , nnam, a1, a2
		      ) ; 
	  vtxtPrintf (blkp, ". Tracks may be  <font color='red'>slow to load</font>; please reload if some tracks come up yellow-greenish, and thanks to UCSC for the great work!") ;
	}
      gmpCaption (blkp, gmp, "Gene_expression_profile_in_primates_caption", "Read more...") ;
      vtxtBreak (blkp) ;
      if (1) /* caption */
	{
	  vtxtPrintf (blkp, "About UCSC tracks: ") ;
	  vtxtPrintf (blkp, "you may enjoy the plots for the summed coverage over all primates' libraries (top track), summarizing 3 terabases of stranded RNA-seq. Fragments mapping on the + strand of the genome (from genes on the + strand) are red (or  dark), on minus strand blue (or light) and antisense transcribed areas are black or overlaid. The vertical scale for each track is self-adapting. Homozygous SNPs tracks are also presented") ;
	  vtxtBreak (blkp) ;

	  vtxtPrintf (blkp, "About mapping: ") ;
	  vtxtPrintf (blkp, "Primates body map RNA-seq data were stringently mapped to the human genome using the NCBI Magic pipeline. Normalized results are shown as significant FPKM (sFPKM), which includes corrections on F, K and M, computed from parameters measured directly in each RNA-seq experiment, to render the expression measures more significant and more robust to experimental biases") ;

	  vtxtPrintf (blkp, ". Only fragments with both reads mapped uniquely and over at least 80+80 bases ending with 8 exact bases on each side of each read, and facing each other in a single site or gene, are included in the computation of the sFPKM/index, in the coverage plots, and in the determination of homozygous SNPs (minimum coverage 10, minimum allele frequency 95%). But be aware that genes whose sequence evolved to become too distant from Human cannot be measured well, this bias can be appreciated in the per-species coverage plots at UCSC.") ;
	  vtxtBreak (blkp) ;
	  vtxtPrintf (blkp, "About libraries:") ;
	  vtxtPrintf (blkp, " For non-human primates, total RNA libraries used TruSeq, ribozero and the stranded UDG protocol. The human 2010 libraries used the polyA selected non-stranded protocol, with short reads (50, 75 or 50+50 bases); furthermore the insert lengths are larger in human than in the non-human primates (average insert size 187 bp in non-human primates versus 232 bp in human). These protocol differences may impact expression measures for the non polyadenylated genes (or genes with shorter or occasional polyA tails), for the pseudogenes or close gene families (specificity is reduced in humans due to shorter reads), and for the very short genes.") ;

	  vtxtBreak (blkp) ;
		      
	  if (0)  vtxtPrintf (blkp, "The SNP line indicates the number of homozyguous SNPs in each species and conservation of the sequence is color coded.") ;
	   vtxtBreak (blkp) ;
		      
	}

      gmpChapterClose (blkp, gmp, "RNA_seq", TRUE) ;
    } 
  
  ac_free (h) ;
} /* ficheNewGeneExpressionProfileChapter */

/***************************************************************************************/

void ficheNewGeneLocatorDiagramChapter (vTXT blkp, GMP *gmp, BOOL isSmall)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *ccp ;
  
  if (gmp->gene)
    {
      ccp =  ac_tag_printable (gmp->gene, "Cytogenetic", ac_tag_printable (gmp->gene, "IntMap", "")) ;
      
      if(! strncmp (ccp, "CHROMOSOME_", 11))
	ccp += 11 ;
      gmpChapter (blkp, gmp, "#*Gene_locator", messprintf ("Gene neighbors and Navigator on chromosome %s", ccp)) ;
      if (0)
	{
	  vtxtPrint (blkp,"Pink arrow: Protein coding gene, ") ;
	  vtxtPrint (blkp,"Red: spliced non protein coding. Click on a gene to navigate.<br> ") ;
	}

      vtxtPrintf (blkp, "<div class='%s' id='locatorSmall'>\n", isSmall ? "shown" : "hidden" ) ;
      if (HTML5)
	{
	  int len = 0 ;
	  AC_DB db = ac_open_db ("local", 0) ;
	  char *qq = hprintf (h, "view  -c Gene -n %s  -v DtGLOC -svg", ac_protect (ac_name (gmp->gene), h)) ;
	  vtxtPrint (blkp, (char *)ac_command (db, qq, &len, h)) ;  
	}
      vtxtBreak (blkp) ;
      vtxtPrint (blkp, "<a href=\"javascript:locatorZoom(2)\">ZOOM OUT</a>\n") ;
      if (1)
	{
	  vtxtPrintf (blkp, "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;") ;
	  vtxtPrint (blkp,"<font color='red'> D</font>:disease,") ; 
	  vtxtPrint (blkp,"<font color='brown'> C</font>:conserved,") ; 
	  vtxtPrint (blkp,"<font color='green'> I</font>:interactions,") ; 
	  vtxtPrint (blkp,"<font color='blue'> R</font>:regulation,") ; 
	  vtxtPrint (blkp,"<font color='purple'> P</font>:publications &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ") ; 
	}
      vtxtPrintf (blkp, "<a  href='LegendNeighborsDiagram.pdf' target='AceV legend'>Read more...</a>") ;
      vtxtPrint (blkp, "</div>\n") ;    

      vtxtPrintf (blkp, "<div class='%s' id='locatorBig'>\n", isSmall ? "hidden" : "shown") ;
      if (HTML5)
	{
	  int len = 0 ;
	  AC_DB db = ac_open_db ("local", 0) ;
	  char *qq = hprintf (h, "view  -c Gene -n %s  -v DtGLOCBIG -svg", ac_protect (ac_name (gmp->gene), h)) ;
	  vtxtPrint (blkp, (char *)ac_command (db, qq, &len, h)) ;  
	}
      vtxtBreak (blkp) ;
      vtxtPrint (blkp, "<a href=\"javascript:locatorZoom(1)\">ZOOM IN</a>\n") ;
      if (1)
	{
	  vtxtPrintf (blkp, "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;") ;
	  vtxtPrint (blkp,"<font color='red'>D</font>:disease,") ; 
	  vtxtPrint (blkp,"<font color='brown'>C</font>:conserved,") ; 
	  vtxtPrint (blkp,"<font color='green'>I</font>:interactions,") ; 
	  vtxtPrint (blkp,"<font color='blue'>R</font>:regulation,") ; 
	  vtxtPrint (blkp,"<font color='purple'>P</font>:publications &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ") ; 
	}
      vtxtPrintf (blkp, "<a  href='LegendNeighborsDiagram.htm' target='AceV legend'>Read more...</a>") ;
      vtxtPrint (blkp, "</div>\n") ;

 
      vtxtBreak (blkp) ;

      if (0)
	{
	  gmpCaption (blkp, gmp, "Gene_locator_caption", "Legend: The gene and its neighbors") ;
	  vtxtBreak (blkp) ;
	  vtxtPrint (blkp,
		     "A zoomable chromosomal region, centered around the gene of "
		     "interest (named in red), is shown to scale") ;
	  
	  vtxtBreak (blkp) ;
	  
	  vtxtPrint (blkp,"Each gene glyph is an arrow (in the direction of transcription) "
		     "covering the extent and strand of the clustered cDNA sequences from "
		     "GenBank and dbEST belonging specifically to the gene"
		     ) ;
	  vtxtPrint (blkp,"<br>The color of the arrow indicates the type of gene: " 
		     "<font color='#ff33cc'>pink if protein coding</font>, "
		     "<font color='red'> red if spliced and not-protein coding</font>, "
		     "<font color='blue'>blue if unspliced but relatively well expressed</font>  (putative gene) "
		     "black if neither spliced nor obviously coding and with a low level "
		     "of expression (i.e. what we call a cloud, "
		     "which appears enriched in intronic areas). "
		     "<br>"
		     "The width of the gene glyph indicates the span <font color='#993366'>level of expression</font>"
		     "<br>"
		     "Knowledge about the gene is summarized in the colored pastilles."
		     " For details see the <a href='DetailsNeighborsDiagram.htm'>help</a>"
		     ) ;
	  vtxtBreak (blkp) ;
	  vtxtPrint (blkp, 
		     "<i> This flash diagram is an experimental new development: you may use the right mouse "
		     "button to zoom, drag and print. We count on you to  report "
		     "<a href=\"mailto:mieg@ncbi.nlm.nih.gov?subject=problem_with_flash_diagram\">  any problem</a></i>"
		     ) ;
	}
      gmpChapterClose (blkp, gmp, "Gene_locator", TRUE) ;
    }
  ac_free (h) ;
} /* ficheNewGeneFlashLocatorChapter */

/***************************************************************************************/

static void ficheNewGeneBestFriends (vTXT blkp, GMP *gmp, Array bb)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  vTXT bfr = vtxtHandleCreate (h) ;
  int i, j, n, row ;
  int cc, oldcc ;
  ITRC *itrc, *ip, *jp ;
  int nItrc1 = arrayMax (bb), nItrc2 ;
  AC_OBJ Gene = 0 ;
  enum { COL_SCORE=0, COL_GENE, COL_EVIDENCE, COL_TITLE, COL_LAST} ; 
  int cols[] = {3,2,4,0} ;
  const char *colNames[] = {"Score", "Gene", "From", "Title", "Last"} ;
  const char *ccp ;

  if (gmp->markup) vtxtMarkup (bfr) ;

  arraySort (bb, itrcGeneOrder) ;
  for (i = 0 ; i < arrayMax (bb) ; i++)
    {
      ip = arrp (bb, i, ITRC) ; 
      for (j = i + 1, jp = ip + 1 ; ip->gene && ip->gene == jp->gene && j < arrayMax (bb) ; jp++, j++)
	if (jp->gene == ip->gene)
	  { 
	    jp->gene = 0 ; ip->nn += jp->nn ; jp->nn = 0 ; ip->cc += jp->cc ; jp->cc = 0 ; 
	    ip->type |= jp->type ; jp->type = 0 ; 
	  }
    }
  arraySort (bb, itrcOrder) ;
  for (n = i = 0 ; i < arrayMax(bb) && i < 8 ; i++)
    if (arr (bb, i, ITRC).nn > 1) n = 1 ; 
  itrc = arrp (bb, 0, ITRC) ;
  cc = itrc->cc ;
  if (n)
    {
      gmpSection (blkp, gmp, "#tg_best_friends", "Genes most related by function") ;
      
      /* You find this gene interesting? Learn about its friends and close relatives! */
      
      vtxtPrintf (blkp, "These are the most related genes, through combined annotation of diseases: <font color='red'>D</font>, Kegg pathways: <font color='blue'>W</font>, processes or molecular functions (GO): <font color='blue'>G</font>, cellular localization: L, domains or motifs:M and interactions: <font color='#009c3a'>I</font>. We give more weight to the more specific annotations, and to those supported by publications") ;
      vtxtBreak (blkp) ;
      vtxtPrintf (blkp, "It is just our best guess. Tell us if it helps you discover something interesting (and please ") ; 
      gmpURL (blkp, gmp, "http://genomebiology.com/2006/7/S1/S12", "cite us!") ;
      vtxtPrintf (blkp, ")") ;
      vtxtBreak (blkp) ;
      
      /* contruct the table */
      tbl = ac_empty_table (80, 6, h) ; /* tbl->rows is a hint, not a hard limit */
      oldcc = -1 ; 
      if (nItrc1 > 3) cc = arr (bb, 2, ITRC).cc ;
      for (nItrc2 = 0, row = -1 ; nItrc2 < nItrc1 ; nItrc2++)
	{
	  itrc = arrp (bb, nItrc2, ITRC) ;
	  if (itrc->cc < 1) break ;
	  if (itrc->cc < oldcc && row > 80 && ! FDEBUG) break ;
	  if (row > 50) break ; /* was 200 */
	  oldcc = itrc->cc ;
	  if ((itrc->cc > 1 && 4*itrc->cc >= cc) ||  FDEBUG)
	    {
	      row++ ;
	      ac_table_insert_text (tbl, row, COL_SCORE, messprintf ("%d", (int)itrc->cc)) ;
	      
	      vtxtClear (bfr) ;
	      Gene = ac_get_obj (gmp->db, "gene", name(itrc->gene), 0) ; 
	      gmpObjLink (bfr, gmp, Gene, ac_name (Gene)) ;
	      ac_table_insert_text (tbl, row, COL_GENE, vtxtPtr (bfr)) ;
	      
	      vtxtClear (bfr) ;
	      vtxtPrint (bfr, "<font face='courier'>") ;
	      if (itrc->type & Dtype) vtxtPrint (bfr, "<font color='red'>D</font>") ; else vtxtPrint (bfr, "-") ;
	      if (itrc->type & Wtype) vtxtPrint (bfr, "<font color='blue'>W</font>") ; else vtxtPrint (bfr, "-") ;
	      if (itrc->type & Gtype) vtxtPrint (bfr, "<font color='blue'>G</font>") ; else vtxtPrint (bfr, "-") ;
	      if (itrc->type & Ltype) vtxtPrint (bfr, "<font color='brown'>L</font>") ; else vtxtPrint (bfr, "-") ;
	      if (itrc->type & Mtype) vtxtPrint (bfr, "M") ; else vtxtPrint (bfr, "-") ;
	      if (itrc->type & (IPtype|IGtype)) vtxtPrint (bfr, "<font color='#009c3a'>I</font>") ; else vtxtPrint (bfr, "-") ;
	      /*
		if (itrc->type & Ttype) vtxtPrint (bfr, "T") ; else vtxtPrint (bfr, "-") ;
		if (itrc->type & Etype) vtxtPrint (bfr, "E") ; else vtxtPrint (bfr, "-") ;
	      */
	      vtxtPrint (bfr, "</font>") ;

	      ac_table_insert_text (tbl, row, COL_EVIDENCE, vtxtPtr (bfr)) ;

	      ccp = ac_tag_printable (Gene, "Title", 0) ;
	      if (ccp) 
		ac_table_insert_text (tbl, row, COL_TITLE, ccp) ;

	      ac_free (Gene) ;
	    }
	}
  
     /* export */
     if (row >= 0)
       ac_table_display (blkp 
			 , tbl, colNames
			 , cols, 1
			 , 0, 0, 0
			 , 0
			 ) ;
     
     vtxtBreak (blkp) ;
    }
  ac_free (h) ;
} /* ficheNewGeneBestFriends */

/***************************************************************************************/

void ficheNewGenePhenotypeFunctionChapter (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  vTXT bfr = vtxtHandleCreate (h) ; 
  vTXT bfr1 = vtxtHandleCreate (h) ; 
  Array bb = arrayHandleCreate (256, ITRC, h) ;

  if (gmp->markup) vtxtMarkup (bfr) ;   
  if (gmp->markup) vtxtMarkup (bfr1) ;   
  
  gmpChapter (bfr1, gmp, "Gene_FUNCTION", "BIOLOGICAL AND FUNCTIONAL ANNOTATION") ;

  ficheNewGeneDisease (bfr, gmp, bb) ;  /* type 1 */
  ficheNewGenePhenotypeParagraph (bfr, gmp) ; /* WORM verbose mode */

  ficheNewGeneRegulationParagraph (bfr, gmp) ; 

  ficheNewGenePathwaysProcessFunctionLocalizationParagraph (bfr, gmp, bb) ; /* type = 2 */

  ficheNewGeneProductPfamPsortParagraph (bfr, gmp, bb) ; /*  type Pfam=5  */
  ficheNewGeneInteractionParagraph (bfr, gmp, bb) ;

  if (arrayMax (bb))
    ficheNewGeneBestFriends (bfr, gmp, bb) ;
  ficheNewGeneAceKogSubSection (bfr, gmp) ;
  ficheNewGeneTaxblastParagraph (bfr, gmp) ;

  if (vtxtPtr (bfr))
    { 
      vtxtPrint (blkp, vtxtPtr (bfr1)) ; 

      vtxtPrint (blkp
		 , "In the spirit of systems biology, this chapter provides links"
		 " to all genes with similar annotations") ;
      vtxtBreak (blkp) ;
      vtxtPrint (blkp, vtxtPtr (bfr)) ;
      gmpChapterClose (blkp, gmp, "Gene_FUNCTION", TRUE) ;
    }
  else
    gmpChapterClose (bfr1, gmp, "Gene_FUNCTION", FALSE) ;

  ficheNewGeneBiblioChapter (blkp, gmp, 1) ;
  ac_free (h) ;
} /* ficheNewGenePhenotypeFunctionChapter */

/***************************************************************************************/

static int ficheNewGeneIntronsTable (vTXT blkp, GMP *gmp)
{
  vTXT bfr, bfr2 ; 
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  KEY myMrna = gmp->view == 'm' ? ac_obj_key (gmp->mrna) : 0 ;
  AC_OBJ oCosmid = ac_tag_obj (gmp->tg, "Genomic_sequence", h) ;
  AC_TABLE gCovers = ac_tag_table (gmp->tg, "Covers", h) ;
  AC_TABLE gIntMap = ac_tag_table (gmp->tg, "IntMap", h) ;
  AC_TABLE gSplicing = ac_tag_table (gmp->tg, "Splicing", h) ;
  AC_TABLE gFuzzy = ac_tag_table (gmp->tg, "Fuzzy", h) ;
  AC_TABLE gAssembled_from = ac_tag_table (gmp->tg, "Assembled_from", h) ; 
  BOOL specific, oldSpecific = FALSE ;
  int maxRows, maxCols, ir, ga1, ga2, ca1, ca2 ;
  int jt ;
  int jExon = 1, jIntron = 1 ;
  Array Tbb = 0 ;
  const char *txt ;
  const char *txt2;
  char *secColNames[]={ "Type", 
			"Length<br/>& DNA", 
			"Coordinates<br/>on mRNA", 
			"Coordinates<br/>on gene", 
			"Coordinates<br/>on genome", 
			"Supporting<br/>clones", 
			"in variant",
			"Coordinates<br/>on clone", 
			"From tissue (# clones)",
			0} ; 
  enum {COL_INTREX=0, COL_LENG, COL_CORDMRNA, COL_CORDGENE, COL_CORDGENOME, COL_SUPCLONE, COL_SUPMRNA, COL_CORDCLONE, COL_TISSUE,  COL_LAST} ; 
  
  ca1 = ac_table_int (gCovers, 0, 2, 0) ;
  ca2 = ac_table_int (gCovers, 0, 3, 0) ;

  ga1 = ac_table_int (gIntMap, 0, 1, 0) ;
  ga2 = ac_table_int (gIntMap, 0, 2, 0) ;

  bfr = vtxtHandleCreate (h) ; 
  if (gmp->markup) vtxtMarkup (bfr) ;
  bfr2 = vtxtHandleCreate (h) ; 
  if (gmp->markup) vtxtMarkup (bfr2) ;

  if (gSplicing && gAssembled_from)
    maxRows = (gSplicing->rows) * (gAssembled_from->rows) ; 
  else
    { ac_free (h) ;  return 0  ; }
  Tbb = arrayCreate ((TbbDim+1) * (maxRows+2), int) ;
  
  /* The Second Table */
  vtxtPrintf (bfr, " ") ; 
  /* initialize the column titles */
  for (maxCols=0 ; secColNames[maxCols] ; maxCols++)
    TBB(0, maxCols) = vtxtPrintf (bfr, "%s"ooo, secColNames[maxCols]) ; 
  
  /* filling the table */
  for (jt = 0, ir=0 ; ir < gSplicing->rows ; ir++)
    {
      BOOL foundMrna = FALSE ;
      int v3 = 0, v4 = 0 ;  

      ac_free (h1) ;
      h1 = ac_new_handle () ;

      specific = oldSpecific ;
      /* in mRNA */
      {
	int iMrna, jt2 = 0, iss ;
	char *extName = "" ;
	AC_TABLE gMrnas, gmSplicing ;
	AC_OBJ oMrna ;
	int a1, u1, u2, v1 , v2 ;

	gMrnas = ac_tag_table (gmp->tg, "mRNA", h) ;
	if (gMrnas)
	  for (jt2 = iMrna = 0 ; iMrna < gMrnas->rows ; iMrna++)
	    {
	      oMrna = ac_table_obj (gMrnas, iMrna, 0, h) ;
	      a1 = ac_table_int (gMrnas, iMrna, 1, 0) ; /* mrna coord in the gene */
	      
	      gmSplicing = ac_tag_table (oMrna, "Splicing", h) ;
	      if (!gmSplicing) continue ;
	      for (iss = 0 ; iss < gmSplicing->rows ; iss++)
		{
		  u1 = ac_table_int (gSplicing, ir, 0, 0) ; 
		  u2 = ac_table_int (gSplicing, ir, 1, 0) ;
		  v1 = ac_table_int (gmSplicing, iss, 0, 0) ;
		  v2 = ac_table_int (gmSplicing, iss, 1, 0) ;
		  if (v1 + a1 - 1 != u1 ||
		      v2 + a1 - 1 != u2)
		    continue ;
		  extName = cleanVariantName (ac_name (gmp->gene), ac_name(oMrna)) ;
		  if (!jt2++)
		    TBB (jt+1, COL_SUPMRNA) = gmpObjLink (bfr, gmp, oMrna,extName) ;
		  else
		    { 
		      if (jt2 % 3 == 1 && jt2 > 1)
			vtxtPrintf (bfr, ",<br/>", extName) ; 
		      else
			vtxtPrintf (bfr, ",", extName) ; 

		      gmpObjLink (bfr, gmp, oMrna, extName) ;
		    }
		  if (myMrna && myMrna == ac_obj_key (oMrna))
		    {
		      foundMrna = TRUE ; 
		      v3 = ac_table_int (gmSplicing, iss, 2, 0) ;
		      v4 = ac_table_int (gmSplicing, iss, 3, 0) ;
		    }
		  break ; /* break the iss loop */
		}
	    }
	if (jt2)
	  {
	    vtxtPrintf (bfr, ooo) ; 
	  } 
	specific = foundMrna && jt2 == 1 ? TRUE : FALSE ;
      }
      if (myMrna && !foundMrna)
	continue ; 
      jt++ ;
      if (myMrna)
	TBB (jt, COL_CORDMRNA) = vtxtPrintf (bfr, "%d to %d"ooo, v3, v4) ;
      oldSpecific = specific ;
      txt = ac_table_tag (gSplicing, ir, 2, "") ; 

      {
	int iss ; 
	struct {char *lookFor, * whatToSay ; }xonTypes[]={
	  {"Exon", "Exon"}, 
	  {"Alternative_exon", "Alternative exon"}, 
	  {"Alternative_Partial_Exon", "Alternative partial exon"},
	  {"Predicted_exon", "Predicted exon"}, 
	  {"Stolen_exon", "Inferred exon"}, 
	  {"Intron", "Intron"}, 
	  {"Predicted_intron", "Predicted intron"}, 
	  {"Alternative_intron", "Alternative intron"}, 
	  {"Stolen_intron", "Inferred intron"}, 
	  {"Gap", "Sequencing gap"}, 
	  {"ORF_Gap", "Open frame gap"}, 
	  {0, 0}} ; 
	
	
	for (iss=0 ; xonTypes[iss].lookFor ; iss++)
	  {
	    if (!strcasecmp (txt, xonTypes[iss].lookFor))break ; 
	  }
	
	/* exon */
	if (strstr (txt, "xon"))
	  {
	    TBB (jt, COL_INTREX) = vtxtPrintf (bfr, "") ;
	    if (specific) vtxtPrint (bfr, "<font color='red'>") ;
	    vtxtPrintf (bfr, "<b>%s</b>",  xonTypes[iss].whatToSay ? xonTypes[iss].whatToSay : txt) ; 
	    if (specific) vtxtPrint (bfr, "</font>") ;
	    vtxtPrintf (bfr, " %d"ooo, jExon++) ;
	  }
	/* intron */
	else if (strstr (txt, "ntron"))
	  {
	    txt2 = ac_table_printable (gSplicing, ir, 3, "") ;
	    TBB (jt, COL_INTREX) = vtxtPrintf (bfr, "") ;
	    if (specific) vtxtPrint (bfr, "<font color='red'>") ;
	    vtxtPrintf (bfr, "%s %d",  xonTypes[iss].whatToSay ? xonTypes[iss].whatToSay : txt, jIntron++) ; 
	    if (!strcasecmp (txt2, "Fuzzy"))
	      {
		int xx1 = ac_table_int (gSplicing, ir, 0, 0) ;
		int xx2 = ac_table_int (gSplicing, ir, 1, 0) ;

		vtxtPrintf (bfr, " Fuzzy ") ; 
		if (gFuzzy)
		  {
		    int zz1, zz2, ifz ;
		    for (ifz = 0 ; ifz < gFuzzy->rows ; ifz++)
		      {
			zz1 = ac_table_int (gFuzzy, ifz, 1, -1) ;
			zz2 = ac_table_int (gFuzzy, ifz, 2, -1) ;
			if (zz1 == xx1 && zz2 == xx2)
			  {			    
			    txt2 = ac_table_printable (gFuzzy, ifz, 4, 0) ;
			    if (txt2)
			      vtxtPrint (bfr, txt2) ;
			    break ;
			  }
		      }
		  }
	      }
	    else 
	      {
		if (strcmp (txt2, "gt_ag") && strcmp (txt2, "gc_ag"))
		  vtxtPrintf (bfr, " <font color='#purple'>[%c%c-%c%c]</font>", txt2[0], txt2[1], txt2[3], txt2[4]) ; 		else
		  vtxtPrintf (bfr, " [%c%c-%c%c]", txt2[0], txt2[1], txt2[3], txt2[4]) ; 
	      }
	    if (specific) vtxtPrint (bfr, "</font>") ;
	    vtxtPrint (bfr, ooo) ; 
	  }
	else {
	  TBB (jt, COL_INTREX)=vtxtPrintf (bfr, "%s", xonTypes[iss].whatToSay ? xonTypes[iss].whatToSay : txt) ; 
	  vtxtPrintf (bfr, ooo) ; 
	}
	
      }
      
      /* length & coordinate on mRNA and Gene */
      {  /* view -v fasta -c mrna -n mog-6 -xml -p 12 35 */
	char buf1[256], buf2[256] ;
	int u1, u2, xx1, xx2, len ;
	int exonColor ;

	xx1 = ac_table_int (gSplicing, ir, 0, 0) ;
	xx2 = ac_table_int (gSplicing, ir, 1, 0) ;
	len = xx2 - xx1 + 1 ;
	if (ca1 < ca2)
	  { u1 = ca1 + xx1 - 1 ; u2 = ca1 + xx2 - 1 ; }
	else
	  { u1 = ca1 - xx1 + 1 ; u2 = ca1 - xx2 + 1 ; }
	TBB (jt, COL_LENG) = vtxtPrintf (bfr, " ") ;

	if (strstr (ac_table_printable (gSplicing, ir, 2, ""), "xon"))
	  {
	    if (jExon % 2)
	      exonColor = 1 ;
	    else
	      exonColor = 2 ;
	  }
	else
	  exonColor = 0 ;
	sprintf (buf1, "DNA:%d:%d:%d", u1, u2, exonColor) ;
	sprintf (buf2, "%d&nbsp;bp", len) ;
	gmpFakeObjLink (bfr, gmp, buf1, oCosmid, buf2) ;
      
	vtxtPrintf (bfr, ooo) ; 

	if (ga1 < ga2)
	  { u1 = ga1 + xx1 - 1 ; u2 = ga1 + xx2 - 1 ; }
	else
	  { u1 = ga1 - xx1 + 1 ; u2 = ga1 - xx2 + 1 ; }

 	TBB (jt, COL_CORDGENE) = vtxtPrintf (bfr, "%d to %d"ooo, xx1, xx2) ;
 	TBB (jt, COL_CORDGENOME) = vtxtPrintf (bfr, "%d to %d"ooo, u1, u2) ;
      }
      
      /* supporting Clone */
      {
	int iEst, iEst2, jt2 = 0 ;
	AC_OBJ est ;
	int u1 = 0, u2 = 0, v1 , v2, irTbl = 0 ;
	AC_TABLE rTbl = ac_db_empty_table (gmp->db, 30, 1, h) ;

	vtxtClear (bfr2) ;
	if (strstr (txt, "xon"))
	  for (jt2 = iEst = 0 ; iEst < gAssembled_from->rows ; iEst++)
	    {
	      est = ac_table_obj (gAssembled_from, iEst, 2, h) ; 
	      u1 = ac_table_int (gSplicing, ir, 0, 0) ;
	      u2 = ac_table_int (gSplicing, ir, 1, 0) ;
	      v1 = ac_table_int (gAssembled_from, iEst, 0, 0) ;
	      v2 = ac_table_int (gAssembled_from, iEst, 1, 0) ;
	      
	      if (u1 != v1 || u2 != v2)
		continue ;
	      if (ac_has_tag (est, "Tissue"))
		ac_table_insert_type (rTbl, irTbl++, 0, &est, ac_type_obj) ;
	      /* cord on est */
	      if (jt2 == 0)  /* was jt2 < 3*/
		{
		  TBB (jt, COL_CORDCLONE) = 
		    vtxtPrintf (bfr, "%d to %d"ooo
				, ac_table_int (gAssembled_from, iEst, 3, 0) 
				, ac_table_int (gAssembled_from, iEst, 4, 0) ) ;
		  /* TBB (jt, COL_SUPCLONE) = gmpObjLink (bfr, gmp, est, 0) ;  */
		  TBB (jt, COL_SUPCLONE) = vtxtPrint (bfr, "") ;
		  vtxtPrintf (bfr2, "%s", ac_name (est)) ;
		} 
	      jt2++ ;
	    }
	if (strstr (txt, "tron"))
	  for (jt2 = iEst = 0 ; iEst < gAssembled_from->rows ; iEst++)
	    {
	      est = ac_table_obj (gAssembled_from, iEst, 2, h) ; 
	      u1 = ac_table_int (gSplicing, ir, 0, 0) ; 
	      u2 = ac_table_int (gSplicing, ir, 1, 0) ;
	      v2 = ac_table_int (gAssembled_from, iEst, 1, 0) ;
	      if (u1 != v2 + 1)
		continue ;
	       for (iEst2 = iEst + 1 ; iEst2 < gAssembled_from->rows ; iEst2++) 
		 {
		   v2 = ac_table_int (gAssembled_from, iEst2, 0, 0) ;
		   if (u2 != v2 - 1 ||
		       strcmp (ac_name (est), ac_table_printable (gAssembled_from, iEst2, 2, "")))
		     continue ;
		   if (ac_has_tag (est, "Tissue"))
		     ac_table_insert_type (rTbl, irTbl++, 0, &est, ac_type_obj) ;
		   if (jt2 == 0)  /* was jt2 < 3*/
		     {
		       TBB (jt, COL_CORDCLONE) = 
			 vtxtPrintf (bfr, "%d to %d"ooo
				     , ac_table_int (gAssembled_from, iEst, 4, 0)
				     , ac_table_int (gAssembled_from, iEst2, 3, 0)) ;
		       /* TBB (jt, COL_SUPCLONE) = gmpObjLink (bfr, gmp, est, 0) ;  */
		       TBB (jt, COL_SUPCLONE) = vtxtPrint (bfr, " ") ;
		       vtxtPrint (bfr2, ac_name (est)) ;
		     }
		   jt2++ ;
		 }
	    }
	if (vtxtPtr (bfr2))
	  {
	    char buf1[256] ;
	    
	    sprintf (buf1, "Tg_support:%d:%d", u1, u2) ;
	    if (jt2 > 1) vtxtPrintf (bfr2, "<br/>and %d other%s", jt2 - 1, _multi(jt2-1)) ;
	    gmpFakeObjLink (bfr, gmp, buf1, gmp->tg, vtxtPtr (bfr2)) ;
	    vtxtPrintf (bfr, ooo) ;
	  }
	if (irTbl)
	  {
	    AC_KEYSET rks = ac_table_keyset (gmp->db, rTbl, 0, h) ;
	    AC_KEYSET clones = ac_ksquery_keyset (rks, ">cdna_clone", h) ;
	    TBB (jt, COL_TISSUE) = vtxtPrint (bfr, "") ;
	    if (specific) vtxtPrint (bfr, "<font color='red'>") ;
	    ficheNewGeneExpressionTissue (bfr, gmp, clones, 2, 4, 'z', -1) ; 
	    if (specific) vtxtPrint (bfr, "</font>") ;
	    vtxtPrintf (bfr, ooo) ;
	  }
      }

    }
  vtxtDot (blkp) ;
  vtxtPrintf (blkp, "Sequences of exons and introns are available by clicking on the lengths in column 3") ;
  if (myMrna)
    vtxtPrintf (blkp, ". Red corresponds to exons or introns specific of mRNA variant %s", 
		cleanVariantName (ac_name (gmp->gene), ac_name(gmp->mrna))) ;
  vtxtBreak (blkp) ;
  if (myMrna)
    fichePrintSquareTable (gmp, gmp->style, Tbb, blkp, bfr, 0, 0, 0, jt + 1,
			   COL_INTREX,  COL_SUPMRNA, COL_LENG, COL_CORDMRNA, COL_CORDGENE, COL_CORDGENOME, COL_SUPCLONE, COL_TISSUE, -1) ; 
  else
    fichePrintSquareTable (gmp, gmp->style, Tbb, blkp, bfr, 0, 0, 0, jt + 1,
			   COL_INTREX, COL_SUPMRNA, COL_LENG, COL_CORDGENE, COL_CORDGENOME, COL_SUPCLONE, COL_TISSUE, -1) ; 
  
  arrayDestroy (Tbb) ;
  vtxtDestroy (bfr) ;

  if (0) vtxtPrintf (blkp,
	      "\nA clone supports an exon or an intron if it has exactly the same boundaries. "
	      "A specified intron, either typical [gt-ag] or [gc-ag] both shown in pink, "
	      "or atypical and shown in blue on the drawing, has at least one clone "
	      "exactly matching the genome over 8 bp on each side. "
	      "Some supported exons or introns may be shown, "
	      "although the corresponding variants are not displayed. "
	      "If an exon is supported by overlapping clones, "
	      "they are not listed. This is frequently the case for the last (and first) exon, "
	      "because alternative polyadenylation is so prevalent that "
	      "we have chosen to merge and show only the longest 3'UTR. "
	      "All features in the table (up to programming bugs) are supported "
	      "by mRNAs or ESTs from the public databases (DDBJ/EMBL/GenBank)."
	      ) ;
  ac_free (h1) ;
  ac_free (h) ;
  
  return 1 ; 
} /* ficheNewGeneIntronsTable */

/*****************/

static void ficheNewGeneIntronsParagraph (vTXT blkp, GMP *gmp)
{
  vTXT bfr ; 
  bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  if (gmp->view == 'g') ficheNewGeneIntronsStatement (bfr, gmp, TRUE, FALSE) ;
  ficheNewGeneIntronsTable (bfr, gmp) ;
  if (vtxtPtr (bfr))
    {
      gmpSection (blkp, gmp, "tg_introns"
		  , "Introns and exons sequence, cDNA support, and tissue") ;
      vtxtPrint (blkp, vtxtPtr (bfr)) ; 
    }
  vtxtDestroy (bfr) ;
} /* ficheNewGeneIntronsParagraph */

/***************************************************************************************/
/***************************************************************************************/
/* auto configure the titles and visible columns using the -: and 4: constructions 
 *  -> separate the table-makewr title from the web title
 */  
 
static int ficheTableConfigureColumns (const char **titles, const char **trueTitles, int *cols) 
{
  int i, j, j1, colorControl =  0 ;
  const char *cp ;
  
  for (i = 0 ; i < 60 ; i++)
    cols[i] = 0 ;
  for (i = j = 0 ; titles[i] ; i++)
    {
      cp = titles[i] ;
      trueTitles [i] = titles[i] ;
      if (*cp == '+')	
	{cp++; if (! colorControl)  colorControl= i + 1 ;}
      if (*cp == '-')	
	continue ;
      if (sscanf (cp,"%d:",&j1) == 1 && j1 >= 1 && j1 < 60)
	cols[j1-1] = i+1 ; /* this column should be displayed */      
      if ((cp = strstr (titles[i], ":")))
	trueTitles [i] = cp+ 1 ;
      if (cp && (cp = strstr (cp, "->")))
	trueTitles [i] = cp+ 2 ;
    }
  return colorControl ;
} /* ficheTableConfigureColumns */

/***************************************************************************************/
/**************************************************************************************/

static BOOL ficheNewGenePolyATableFormat (vTXT blkp, GMP *gmp, AC_TABLE tbl, int maxLine)
{
  int ir, a1, m1, x1 ;
  const char *ccp ;
  AC_HANDLE h = ac_new_handle () ;

  if (maxLine == 0)
    maxLine = tbl->rows ;
  /* otherwise use maxLine +2 so that 'more' will get trigerred */

  if (tbl) /* format the table */
    for (ir = 0 ; ir < maxLine+2 && ir < tbl->rows ; ir++)
      {
	/* 1: Position on the gene */
	m1 = ac_table_int (tbl, ir, 1, 0) ;
	a1 = ac_table_int (tbl, ir, 2, 0) ;
	ac_table_insert_text (tbl, ir, 2, messprintf ("%d", a1 + m1 -1)) ;
	
	/* 2: Positionm on mrna */
	x1 = ac_table_int (tbl, ir, 3, 0) ;
	ccp = gtMrnaSuffix (ac_name (gmp->gene), ac_table_printable (tbl, ir, 0, ""), h) ;
	ac_table_insert_text (tbl, ir, 3, messprintf ("bp %d on %s", x1, ccp)) ;

	/* 4: type */
	ccp =  ac_table_printable (tbl, ir, 7, 0) ;
	if (ccp)
	  ac_table_insert_text (tbl, ir, 5, ccp) ;
      }
  ac_free (h) ;
  return TRUE ;
} /* ficheNewGenePolyATableFormat */

/***************************************************************************************/

static void ficheNewGenePolyATable (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = 0 ;
  int colorControl = 1 ;
  const char *errorMessage = 0 ;
  /* keep in synch column enum and clumn names */
  const char *myOrder = 0 ; /* set below in the switch */
  int cols[60] ;
  const char *trueTitles [60] ;
  /* + color control column
   *  - columns are not visible
   * 1:  means see it as column 1
   */
  const char *titles[]={ "-:mRNA"
			 , "-:m1"
			 , "+1:a1->Position on the gene"
			 , "2:x1->Position on the mRNA"
			 , "3:ncl->Number of supporting accessions"
			 , "4:type->Poly-A signal"
			 , "5:amont->Distance to poly-A site"
			 , "-:variant"
			 , 0} ;
  AC_KEYSET ksTg = ac_objquery_keyset (gmp->tg, "mrna", h) ;
  vTXT bqlQ = vtxtHandleCreate (h) ;
  
  if (!ksTg) goto done ;
  /* auto configure the titles and visible columns */
  colorControl = ficheTableConfigureColumns (titles, trueTitles, cols) ;

  /* if the bqlQ or bqlQ2 queries are  modified, their mapping jj[] below must be modified */
  vtxtPrint (bqlQ,  "select m, m1, a1, x1, ncl, type, amont, variant ") ;
  vtxtPrint (bqlQ,  " from tg in @, m in tg->mrna where m#valid3p, m1 in m[1] ") ;
  vtxtPrint (bqlQ,  " , a1 in m->valid3p, x1 in a1[1], ncl in x1[1], type in ncl[1], amont in type[1]") ;
  vtxtPrint (bqlQ,  " , variant in m->variant") ;
  if (gmp->view == 'm')
    vtxtPrintf (bqlQ,  " where m like %s", ac_protect (ac_name(gmp->mrna), h)) ;
  
  /* select line ordering according to user choice */
  myOrder = "+1+2+3+4+5" ;
    
  /* contruct the bql table */
  tbl = ac_bql_table (gmp->db, vtxtPtr(bqlQ), ksTg, myOrder, &errorMessage, h) ;

  /* format the table and add the http links */
  if (tbl)
    ficheNewGenePolyATableFormat (blkp, gmp, tbl, 0) ;

  /* export */
  if (tbl && tbl->rows && tbl->cols < 18)
    ac_table_insert_text (tbl, 0, 17, "toto") ;
  if (tbl && tbl->cols)
    ac_table_display (blkp 
		      , tbl, trueTitles
		      , cols, colorControl
		      , 0, 0, 0
		      , 0
		      ) ;
 done:
  ac_free (h) ;
} /* ficheNewGenePolyATable */

/*****************/

static void ficheNewGenePolyAParagraph (vTXT blkp, GMP *gmp)
{
  vTXT bfr ; 
  bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ; 

  ficheNewGenePolyATable (bfr, gmp) ; /* deplacee dans mrna structure and expression */
  if (vtxtPtr (bfr))
    {
      gmpSection (blkp, gmp, "tg_polyA"
		  , "Validated poly A sites") ;
      vtxtPrint (blkp, vtxtPtr (bfr)) ; 
    }
  vtxtDestroy (bfr) ;
} /* ficheNewGenePolyAParagraph */

/***************************************************************************************/

void ficheNewGeneExpressionChapter (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  Array tissues = arrayHandleCreate (60, TISSUE, h) ;
  vTXT bfr = vtxtHandleCreate (h) ; 
  vTXT bfr1 = vtxtHandleCreate (h) ; 
  int nClo ;
  
  if (gmp->markup) vtxtMarkup (bfr) ;   
  if (gmp->markup) vtxtMarkup (bfr1) ;   

  /*  see if i should open the chapter */
  gmpChapter (bfr1, gmp, "*Gene_EXPRESSION", "Expression and GenBank cDNA support") ;
  nClo = ficheNewGeneExpressionParagraph (bfr, gmp, tissues) ;
  ficheNewTissuesParagraph (bfr1, gmp, tissues
			    , messprintf ("Tissues where expression was observed (from %s cDNA clone%s)"
					  , isOne(nClo), _multi(nClo))
			    , "Origin of the cDNAs, as reported in GenBank/dbEST (tissue, stage, pathological or normal) shows that the gene is expressed in "
			    ) ;

  if (vtxtPtr (bfr))
    { 
      vtxtPrint (blkp, vtxtPtr (bfr1)) ;
      vtxtPrint (blkp, vtxtPtr (bfr)) ;
      gmpChapterClose (blkp, gmp, "Gene_EXPRESSION", TRUE) ;
    }
  else
    gmpChapterClose (bfr1, gmp, "Gene_EXPRESSION_REGULATION", FALSE) ;

  ac_free (h) ;
} /* ficheNewGeneExpressionRegulationChapter */

/***************************************************************************************/
/***************************************************************************************/
/******************* ficheNewGeneMRNAsProteinsChapter ****************************/
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

static int mSizeOrder(const void *a, const void *b)
{
  const CLALI *u = (const CLALI *)a, *v = (const CLALI *)b ;
  return keySetAlphaOrder (u->oMrna, v->oMrna) ;
} /*  mSizeOrder */

/*********/

static int ficheNewGeneAnnotationOfVariantsTable (vTXT blkp, GMP *gmp, BOOL justSequences)
{
  vTXT bfr ; 
  AC_OBJ myMrna = gmp->view == 'm' ? gmp->mrna : 0 ;
  AC_TABLE gIncludes, gMrna=0, gDna, gProduct, gTGProduct = 0, gSplicing = 0 ;
  AC_OBJ oMrna = 0, oProduct = 0, oCosmid = 0 ; 
  AC_OBJ originalMrna = gmp->mrna ; 
  AC_OBJ originalProduct = gmp->product ;
  AC_OBJ originalKantor = gmp->kantor ; 

  const char *txt = 0 ;
  char *prefix = "", namebuf[256], *extName = 0 ;
  int  maxRows, maxCols, curRowClone, ir, addRow, len, mrnaLen, iLine, is, isComma ; 
  Array  Tbb = 0 ;		
  AC_KEYSET lProd=0 ; 
  Array mSize = 0 ;
  CLALI *cla ;
  AC_HANDLE h = ac_new_handle () ;
  char *promotorQuality = "" ;
  char *colNames[]={ "mRNA variant<br>and sequence", 
		     "mRNA variant", 
		     "mRNA matching the genome",
		     "Best predicted protein", 
		     "mRNA from cDNA consensus", 
		     "5' UTR",
		     "3' UTR",
		     "uORF",
		     "Upstream sequence", 
		     "Transcription<br/>unit<br/>pre-mRNA", 
		     "Downstream sequence",
		     "Overview <it>(for structural details  see the table below)</it>", 
		     "5' completeness evidence", 
		     "Sequence gap",
		     "3' completeness evidence", 
		     "# of<br/>exons", 
		     "# of<br/> clones",
		     "From tissue (no strict specificity is implied)",
		     "Domains",
		     "Predicted localization",
		     "Protein", 
		     "Completeness<br>and uniqueness",
		     "Exons<br>in CDS",
		     "Protein<br>quality",
		     "Extends from",
		     "coordinates<br/>on mRNA", 
		     "coordinates<br/>on gene",
		     "coordinates<br/>on genome",
		     "minimal set of<br>supporting clones", 
		     "Features",
		     "representative<br/>clone", 
		     0} ; 
  enum {	
    COL_TRANSCRIPT=0, COL_VARIANT, COL_VSEQ, COL_PSEQ, COL_AMSEQ,
    COL_5PUTR, COL_3PUTR, COL_UORF, COL_PROMOTOR, COL_PREMESSENGER, COL_POSTMOTOR,
    COL_OVERVIEW, COL_5PRIME, COL_SEQGAP, COL_3PRIME, 
    COL_EXON, COL_CLONE, COL_TISSUE, COL_DOMAIN, COL_LOCAL,
    COL_PROTEIN, COL_HIERARCHY, COL_EXON_CDS, COL_PROTEIN_QUALITY, COL_FROM,
    COL_CORDMRNA, COL_CORDGENE, COL_CORDGENOME, COL_SUPPORTING, COL_FEATURES,  COL_BESTCLONE, 
    COL_PF, COL_PFD,
    COL_LAST
  } ; 
  
  bfr = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (bfr) ;

  if (gmp->tg)
    {
      AC_KEYSET mKs = 0 ;
      
      oCosmid = ac_tag_obj (gmp->tg, "Genomic_sequence", h) ;
      if (myMrna)
	mKs = ac_objquery_keyset (myMrna, "from_gene", h) ;
      else if (gmp->gene)
	mKs = ac_objquery_keyset (gmp->gene, ">transcribed_gene ; > mrna", h) ;
      gMrna = mKs ? ac_keyset_table (mKs , 0, -1, 0, h) : 0 ;
      gTGProduct = gmp->gene ? ac_tag_table (gmp->gene, "Product", h) : 0 ;

      ac_free (mKs) ;
    }
  else if (gmp->pg)
    {
      gMrna = ac_tag_table (gmp->pg, "predicted_mrna", h) ; 
      gTGProduct =  gmp->gene ? ac_tag_table (gmp->gene, "Product", h) : 0 ;
    }
  if (gMrna && gTGProduct)
    maxRows= (gMrna->rows)* (gTGProduct->rows) ; 
  else
    { ac_free (h) ; return 0  ;}
  
  Tbb = arrayCreate ((TbbDim+1)* (maxRows+2), int) ;
  
  vtxtPrintf (bfr, "\n") ; 
  
  /* initialize the column titles */
  for (maxCols=0 ; colNames[maxCols] ; maxCols++)
    {
      TBB (0, maxCols)=vtxtPrintf (bfr, "%s"ooo, colNames[maxCols]) ; 
    }
  
  /* get the protein sequences */
  mSize = arrayCreate (gMrna->rows+1, CLALI) ;
  for (ir=0 ; ir < gMrna->rows ; ir++)
    {      
      cla = arrayp (mSize, ir, CLALI) ;
      cla->icl = ir ;
      cla->oMrna = ac_table_obj (gMrna, ir, 0, h) ;
    }
  arraySort (mSize, mSizeOrder) ;

  for (curRowClone=0, ir=0 ; ir < gMrna->rows ; ir++)
    {
      int proteinLength = 0 ;

      cla = arrp (mSize, ir, CLALI) ;
      oMrna= cla->oMrna ;
      gmp->mrna = oMrna ;
      gmp->kantor = 0 ;
      gSplicing = ac_tag_table (oMrna, "Splicing", h) ; 
      if (! (gProduct = ac_tag_table (oMrna, "Product", h)))
	continue ; 
      addRow=0 ; 
      
      strcpy (namebuf, ac_name (oMrna)) ;
      if (gmp->tg)
	{
	  extName = strrchr (namebuf, '.') ;  
	  if (!extName)extName = "a" ; 
	  else extName ++ ; 
	}
      else if (gmp->pg)
	{
	  extName = strrchr (namebuf, '.') ;  
	  if (extName && !strcmp (extName, ".pg"))
	    {
	      *extName = 0 ;
	      extName = strrchr (namebuf, '.') ;  
	    }
	  if (!extName)extName = "a" ; 
	  else extName ++ ; 
	}
      gDna = ac_tag_table (oMrna, "DNA", h) ;
      mrnaLen = ac_table_int (gDna, 0, 1, 0) ; 
      if (extName[0])
	{
	  TBB (curRowClone+1, COL_TRANSCRIPT)= vtxtPrintf (bfr, "") ;
	  vtxtPrintf (bfr, "Variant:") ;
	  gmpObjLink (bfr, gmp, oMrna, extName) ; 
	  {
	    char buf1[256], buf2[256] ;
	    
	    vtxtPrintf (bfr, "<br/>DNA:") ;
	    sprintf (buf1, "DNA_mRNA::") ;
	    sprintf (buf2, "%d&nbsp;bp", mrnaLen) ; 
	    gmpFakeObjLink (bfr, gmp, buf1, oMrna, buf2) ;
	  }
	 
	  if (0) /* 2007_06_08 AM this time are not nice */
	    {
	    char buf1[256], buf2[256] ;
	    AC_OBJ oAm = ac_tag_obj (oMrna, "RefSeqMaker", h) ;
	    char *am_dna = oAm ? ac_obj_dna (oAm, h) : 0 ;
	    int mrnaLen2 = am_dna ? strlen (am_dna) : 0 ;
	    if (am_dna && oAm)
	      {
		vtxtPrintf (bfr, "<br/>AM:") ;
		sprintf (buf1, "DNA::") ;
		sprintf (buf2, "%d&nbsp;bp", mrnaLen2) ; 
		gmpFakeObjLink (bfr, gmp, buf1, oAm, buf2) ;
	      }
	  }
	}
      else TBB (curRowClone+1, COL_TRANSCRIPT)=vtxtPrintf (bfr, "%d&nbsp;bp", mrnaLen) ; 
      
      vtxtPrintf (bfr, ooo) ; 



      if (oMrna)
	{
	  TBB (curRowClone+1, COL_VSEQ)= vtxtPrintf (bfr, "") ;
	  {
	    char buf1[256], buf2[256] ;
	    
	    sprintf (buf1, "DNA_mRNA::") ;
	    sprintf (buf2, "%d&nbsp;bp", mrnaLen) ; 
	    gmpFakeObjLink (bfr, gmp, buf1, oMrna, buf2) ;
	  }
	 vtxtPrintf (bfr, ooo) ; 
	}

      if (1) /* 2007_06_08 AM this time are not nice */
	{
	  char buf1[256], buf2[256] ;
	  AC_OBJ oAm = ac_tag_obj (oMrna, "RefSeqMaker", h) ;
	  char *am_dna = oAm ? ac_obj_dna (oAm, h) : 0 ;
	  int mrnaLen2 = am_dna ? strlen (am_dna) : 0 ;
	  if (am_dna && oAm)
	    { 
	      TBB (curRowClone+1, COL_AMSEQ)= vtxtPrintf (bfr, "") ;
	      sprintf (buf1, "DNA::") ;
	      sprintf (buf2, "%d&nbsp;bp", mrnaLen2) ; 
	      gmpFakeObjLink (bfr, gmp, buf1, oAm, buf2) ;
	      vtxtPrintf (bfr, ooo) ; 
	    }
	}

      
      TBB (curRowClone+1, COL_VARIANT) = gmpObjLink (bfr, gmp, oMrna, extName) ; 
      vtxtPrintf (bfr, ooo) ; 
      vtxtPrintf (bfr,"\n") ; /* protects against vtxtDot() */
     
      if (gmp->view == 'g')
	{
	  TBB (curRowClone+1, COL_OVERVIEW) = vtxtMark (bfr) ; 
	  ac_free (gmp->dna) ; 
	  if (oMrna)
	    {
	      AC_TABLE mygProduct = ac_tag_table (gmp->mrna, "Product", h) ;
	      AC_OBJ myProduct = gmp->product ;
	      
	      for (is=0 ; mygProduct && is < mygProduct->rows ; is++)
		{
		  oProduct = ac_table_obj (mygProduct, is, 0, h) ; 
		  if (! ac_has_tag (oProduct, "very_good_product") && ! ac_has_tag (oProduct, "Best_product"))
		    continue ;
		  gmp->product = oProduct ;
		  gmp->kantor = ac_tag_obj (gmp->product, "Kantor", h) ;
		  if (1 && (txt = ac_tag_printable (gmp->mrna, "Title", 0)))
		    { vtxtPrint (bfr, gtSetUpper (txt)) ; vtxtBreak (bfr) ; }
		  ficheMRNAOverviewParagraphContent (bfr, gmp) ; 
		  break ;
		}	      
	      gmp->product =  myProduct ;
	    }
	  
	  vtxtPrintf (bfr, "."ooo) ; 
	}
      for (iLine=0, is=0 ; is < gProduct->rows ; is++)
	{
	  DICT *dict = 0 ;
	  int iss ;
	  char *locExt ;
	  AC_OBJ oP2 = 0 ;
	  AC_TABLE gPep = 0 ;

	  oProduct = ac_table_obj (gProduct, is, 0, h) ; 
	  if (gmp->view == 'm' && gmp->product && ! ac_has_tag (gmp->product, "Best_product"))
	    {
	      if (!ac_obj_equal (oProduct, gmp->product))
		continue ;
	    }
	  else
	    {
	      if (! ac_has_tag (oProduct, "very_good_product") && ! ac_has_tag (oProduct, "Best_product"))
		continue ;
	    }

	  ac_free (gmp->peptide) ; gmp->product = oProduct ;
	  TBB (curRowClone+1+iLine, COL_PROTEIN) = gmpObjLink (bfr, gmp, oMrna, gtMrnaSuffix(ac_name(gmp->tg), ac_name(oProduct), h)) ;
	  vtxtPrint (bfr, ooo) ; 

	  len = 0 ;
	  if (!len)
	    len = ac_tag_int (oProduct, "Coding_length", 0)/3 ; 
	  if (!len)
	    len = ac_tag_int (oProduct, "Open_length", 0)/3 ;
	  if (! len)
	    {
	      gPep = ac_tag_table (oProduct, "Peptide", h) ;
	      if (gPep)
		{
		  len = ac_table_int (gPep, 0, 1, 0) ;
		  ac_free (gPep) ;
		}
	    }

	  proteinLength = len ;
	  vtxtPrint (bfr, ooo) ; 

	  if (len)
	    {
	      char buf1[256], buf2[256] ;
		      
	      sprintf (buf1, "PEP_Product::") ;
	      sprintf (buf2, "%d&nbsp;aa", len) ; 
	      TBB (curRowClone+1+iLine, COL_PSEQ) =  vtxtPrint (bfr, "") ;
	      gmpFakeObjLink (bfr, gmp, buf1, oProduct, buf2) ;
	      vtxtPrint (bfr, ooo) ; 
	    }


	  TBB (curRowClone+1+iLine, COL_HIERARCHY) =  vtxtPrint (bfr, "") ;
	  if (0 && gProduct->rows>1)
	    vtxtPrintf (bfr, "#%d", is+1) ; 
	  
	  if (ac_has_tag (oProduct, "NH2_Complete") &&  
	      ac_has_tag (oProduct, "at_position_1")&&  
	      ac_has_tag (oProduct, "COOH_Complete")
	      )
	    {
	      vtxtPrint (bfr, " complete<br/>") ; 
	      prefix="" ; 
	    }
	  else if (ac_has_tag (oProduct, "NH2_Complete") &&  
		   ac_has_tag (oProduct, "at_position_1") &&  
		   !ac_has_tag (oProduct, "COOH_Complete")
		   )
	    {
	      vtxtPrint (bfr, " NH2 complete<br/>") ; 
	      prefix="" ; 
	    }
	  else if (
		   (!ac_has_tag (oProduct, "NH2_Complete") || !ac_has_tag (oProduct, "at_position_1")) &&
		   ac_has_tag (oProduct, "COOH_Complete")
		   )
	    {
	      vtxtPrint (bfr, " COOH complete<br/>") ; 
	      prefix="" ; 
	    }
	  else if (
		   (!ac_has_tag (oProduct, "NH2_Complete") || !ac_has_tag (oProduct, "at_position_1")) &&
		   !ac_has_tag (oProduct, "COOH_Complete")
		   )
	    {
	      vtxtPrint (bfr, " partial<br/>") ; 
	      prefix="" ; 
	    }

	  if (0 && /* also elsewhere */
	      ac_has_tag (oProduct, "Good_product") &&  
	      ac_has_tag (oProduct, "NH2_Complete") &&
	      ac_has_tag (oProduct, "at_position_1") &&
	      ac_tag_int  (oProduct, "First_Kozak", 0) == 1
	      )
	    {
	      AC_TABLE koz = ac_tag_table (oProduct, "First_Kozak", h) ;
	      const char *ccp ;

	      if (koz && (ccp = ac_table_printable (koz, 0, 2, 0)))
		vtxtPrintf (bfr, " Starts on %s<br/>", ccp) ; 
	      prefix="" ; 
	    }
	  if (ac_has_tag (oMrna, "Gap"))
	    {
	      vtxtPrint (bfr, " mRNA gap<br/>") ; 
	      prefix="" ; 
	    }
	  
	  /* compare protein similarities */
	  /* equals */
	  if ((gIncludes = ac_tag_table (oProduct, "Identical_to", h)))
	    {
	      prefix="=" ; 
	      for (iss=0 ; iss < gIncludes->rows && iss < 1 ; iss++)
		{
		  oP2 = ac_table_obj (gIncludes, iss, 0, h) ;
		  if (ac_has_tag (oP2, "Best_product"))
		    {
		      locExt = cleanVariantName (ac_name(gmp->gene), ac_table_printable (gIncludes, iss, 0, "")) ;
		      vtxtPrintf (bfr, "%s%s", prefix, locExt) ;
		      prefix=", " ; 
		    }
		  ac_free (oP2) ;
		}
	    }
	  /* sorry but same prot seem to be repeated in included in complete/partial */
	  prefix="<br/>included in " ; 
	  if ((gIncludes = ac_tag_table (oProduct, "Included_in_complete", h)))
	    {
	      dict = dictCreate (30) ;
	      for (iss=0 ; iss < gIncludes->rows ; iss++)
		{
		  dictAdd (dict, ac_table_printable (gIncludes, iss, 0, ""), 0) ;
		  oP2 = ac_table_obj (gIncludes, iss, 0, h) ;
		  if (ac_has_tag (oP2, "Best_product"))
		    {
		      locExt = cleanVariantName (ac_name(gmp->gene), ac_table_printable (gIncludes, iss, 0, "")) ;
		      
		      vtxtPrintf (bfr, "%s%s", prefix, locExt) ; 
		      prefix=", " ; 
		    }
		  ac_free (oP2) ;
		}
	    }
	  if ((gIncludes = ac_tag_table (oProduct, "Included_in_partial", h)))
	    {
	      for (iss=0 ; iss < gIncludes->rows ; iss++)
		{
		  if (dict &&
		      dictFind (dict, ac_table_printable (gIncludes, iss, 0, ""), 0))
		      continue ;
		  oP2 = ac_table_obj (gIncludes, iss, 0, h) ;
		  if (ac_has_tag (oP2, "Best_product"))
		    {
		      locExt = cleanVariantName (ac_name(gmp->gene), ac_table_printable (gIncludes, iss, 0, "")) ;
		      vtxtPrintf (bfr, "%s%s", prefix, locExt) ;  
		      prefix=", " ; 
		    }
		  ac_free (oP2) ;
		}
	    }
	  dictDestroy (dict) ;
	  vtxtPrintf (bfr, ooo) ; 
	  
	  promotorQuality = "" ;
	  /* 5prime */
	  {
	    AC_TABLE g5prime, g5pnext ; 
	    int	slNum, len5pUTR, iss ; 
	    BOOL ok5p = FALSE ;
	    const char *ccp ;

	    len5pUTR = ac_tag_int (oProduct, "Length_5prime_UTR", 0) ;
	    if (len5pUTR > 3)
	      {
		char buf1[256], buf2[256] ;
		sprintf (buf1, "DNA_mRNA:1:%d",len5pUTR) ;
		sprintf (buf2, "%d&nbsp;bp", len5pUTR) ; 
		TBB (curRowClone+1+iLine, COL_5PUTR)=vtxtPrintf (bfr, " ") ; 	
		gmpFakeObjLink (bfr, gmp, buf1, oMrna, buf2) ;
		vtxtPrintf (bfr, ooo) ; 
	      }
	    TBB (curRowClone+1+iLine, COL_5PRIME)=vtxtPrintf (bfr, " ") ; 
	    
	    isComma=0 ; 
	    if (ac_has_tag (oMrna, "Found5p"))
	      {
		int here = vtxtMark (bfr) ;
		if ((g5prime = ac_tag_table (oMrna, "Found5p", h)))
		  {
		    if ((g5pnext = ac_tag_table (oMrna, "Transpliced_to", h)))
		      {
			if (isComma)vtxtPrintf (bfr, ", ") ; 
			for (iss=0 ; iss < g5pnext->rows ; iss++)
			  {
			    ccp = ac_table_printable (g5pnext, iss, 0, "") ; 
			    if (!sscanf (ccp, "SL%d", &slNum))
			      continue ; 
			    promotorQuality = " including Promoter" ;
			    if (slNum==0) ccp = "capped" ; 
			    if (slNum==0 && ok5p) continue ;
			    ok5p = TRUE ;
			    if (! vtxtAt (bfr, here) || !strstr (vtxtAt (bfr, here), ccp))
			      vtxtPrintf (bfr, "%s%s", _comma (iss), ccp) ; 
			  }
		      }
		    if (!ok5p && (g5pnext = ac_tag_table (oMrna, "Aggregated_5p_clones", h)))
		      {
			if (isComma)vtxtPrintf (bfr, ", ") ; 
			vtxtPrintf (bfr, "aggregated clones") ; 
			promotorQuality = " probably including promoter" ;
		      }
		  }
	      }
	    else if (len5pUTR > 3 && ac_keyset_count (ac_objquery_keyset(oMrna, ">product;good_product && best_product && up_stop", h)))
	      {
		vtxtPrintf (bfr, "5' stop") ; promotorQuality = "  possibly including promoter" ;
	      }
	    else
	      vtxtPrintf (bfr, "no evidence") ; 
	    
	    vtxtPrintf (bfr, ooo) ; 
	  }
	  /* 3prime */
	  {
	    int	len3pUTR, len5pUTR,iv ; 
	    AC_OBJ uORF ;
	    AC_TABLE valid ;

	    if ((len3pUTR = ac_tag_int (oProduct, "Length_3prime_UTR", 0)))
	      {
		char buf1[256], buf2[256] ;

		sprintf (buf1, "DNA_mRNA:%d:%d", mrnaLen - len3pUTR + 1, mrnaLen) ;
		sprintf (buf2, "%d&nbsp;bp", len3pUTR) ; 
		TBB (curRowClone+1+iLine, COL_3PUTR)=vtxtPrintf (bfr, " ") ; 	
		gmpFakeObjLink (bfr, gmp, buf1, oMrna, buf2) ;
		vtxtPrintf (bfr, ooo) ; 
	      }
	    if ((uORF = ac_objquery_obj (oMrna, ">product uORF_candidate", h)))
	      {
		char buf1[256], buf2[256] ;
		int Coding_length ;

		len5pUTR = ac_tag_int (uORF, "Length_5prime_UTR", 0) ;
		Coding_length = ac_tag_int (uORF, "Coding_length", 0) ;
		sprintf (buf1, "DNA_mRNA:%d:%d", len5pUTR + 1, len5pUTR + Coding_length) ;
		sprintf (buf2, "%d&nbsp;aa", Coding_length/3) ; 
		sprintf (buf2, " %d&nbsp;bp", Coding_length) ; 
		TBB (curRowClone+1+iLine, COL_UORF)=vtxtPrintf (bfr, " ") ; 	
		gmpFakeObjLink (bfr, gmp, buf1, oMrna, buf2) ;
		vtxtPrintf (bfr, ooo) ; 
	      }
	    
	    isComma=0 ; 
	    TBB (curRowClone+1+iLine, COL_3PRIME)=vtxtPrintf (bfr, " ") ; 
	    if (0 && ac_has_tag (oMrna, "Valid3p"))
	      {
		if (isComma)vtxtPrintf (bfr, ", ") ; 
		valid = ac_tag_table (oMrna, "Valid3p", h) ;
		vtxtPrintf (bfr, "polyA at") ; 
		for (iv = 0 ; iv < valid->rows ; iv++)
		  vtxtPrintf (bfr, "%s %d" 
			      , iv ? "," : ""
			      , ac_table_int (valid, iv, 0, 0) 
			      ) ;
		isComma=1 ; 
	      }
	    else if (ac_has_tag (oMrna, "Valid3p"))
	      {
		vtxtPrintf (bfr, "validated polyA") ; 
	      }
	    else if (len3pUTR > 3 && ac_keyset_count (ac_objquery_keyset(oMrna, ">product;good_product && best_product && down_stop", h)))
	      {
		vtxtPrintf (bfr, "3' stop") ;
	      }
	    else if (ac_keyset_count (ac_objquery_keyset(oMrna, ">product;good_product && best_product && !down_stop", h)))
	      {
		vtxtPrintf (bfr, "Partial") ;
	      }
	    else
	       {
		vtxtPrintf (bfr, "No evidence") ;
	      }
	    
	    vtxtPrintf (bfr, ooo) ; 
	  }
	  
	  /* #exons */
	  {
	    int	numExon=0, iss ; 
	    
	    if (gSplicing)
	      {
		for (iss=0 ; iss < gSplicing->rows ; iss++)
		  {
		    txt = ac_table_printable (gSplicing, iss, 4, "") ; 
		    if (strstr (txt, "xon"))numExon++ ; 
		  }
	      }
	    
	    TBB (curRowClone+1+iLine, COL_EXON)=vtxtPrintf (bfr, "%d"ooo, numExon) ; 
	  }
	  
	  /* #clones and tissues */
	  {
	    AC_TABLE gClone = ac_tag_table (oMrna,"cDNA_clone", 0) ;
	    int	numClone = gClone ? gClone->rows : 0 ;
	    AC_KEYSET clones = ac_objquery_keyset (oMrna, ">cdna_clone", h) ;

	    TBB (curRowClone+1+iLine, COL_CLONE)=vtxtPrintf (bfr, "%d"ooo, numClone) ; 

	    TBB (curRowClone+1+iLine, COL_TISSUE)=vtxtPrintf (bfr, "") ;
	    ficheNewGeneExpressionTissue (bfr, gmp, clones, 2, 2, myMrna ? gmp->view : 'z', -1) ;
	    vtxtPrintf (bfr, ooo) ;

	    ac_free (gClone) ;	    
	    ac_free (clones) ;
	  }
	  
	  /* premessenger */
	  {
	    int lenUnit ;
	    AC_TABLE gCovers = ac_tag_table (oMrna, "Covers", h) ; 

	    lenUnit = ac_table_int (gCovers, 0, 0, 0) ;
	
	    if (lenUnit && oCosmid)
	      {
		char buf1[256], buf2[256] ;
		TBB (curRowClone+1+iLine, COL_PREMESSENGER)=vtxtPrintf (bfr, " ") ; 

		sprintf (buf1, "DNA_mRNA:99999:99999") ;
		sprintf (buf2, "%d&nbsp;bp", lenUnit) ;
		gmpFakeObjLink (bfr, gmp, buf1, oMrna, buf2) ;

		vtxtPrintf (bfr, ooo) ;
	      }	    
	  }
	  /* promotor/postmotor */
	  {
	    int	ma1, ma2, lenUnit ;
	    AC_TABLE gCovers = ac_tag_table (oMrna, "Covers", h) ; 

	    lenUnit = ac_table_int (gCovers, 0, 0, 0) ;
	    ma1 = ac_table_int (gCovers, 0, 2, 0) ;
	    ma2 = ac_table_int (gCovers, 0, 4, 0) ;
	    if (lenUnit && oCosmid)
	      {
		char buf1[256], buf2[256] ;
		TBB (curRowClone+1+iLine, COL_PROMOTOR)=vtxtPrintf (bfr, " ") ; 
		if (ma1 < ma2)
		  sprintf (buf1, "DNA:%d:%d:0", ma1-2000, ma1-1) ;
		else
		  sprintf (buf1, "DNA:%d:%d:0", ma1+2000, ma1+1) ;
		sprintf (buf2, "2kb%s", promotorQuality) ;
		gmpFakeObjLink (bfr, gmp, buf1, oCosmid, buf2) ;

		vtxtPrintf (bfr, ooo) ;

		TBB (curRowClone+1+iLine, COL_POSTMOTOR)=vtxtPrintf (bfr, " ") ; 
		
		if (ma1 < ma2)
		  sprintf (buf1, "DNA:%d:%d:0", ma2+1, ma2 + 1000) ;
		else
		  sprintf (buf1, "DNA:%d:%d:0", ma2-1, ma2-1000) ;


		sprintf (buf2, "1kb") ;
		gmpFakeObjLink (bfr, gmp, buf1, oCosmid, buf2) ;

		vtxtPrintf (bfr, ooo) ;
	      }	    
	  }
	  
	  
	  /* coord on gene */
	  {
	    AC_TABLE imMrna = ac_tag_table (oMrna, "IntMap", h) ;
	    AC_TABLE imGene = ac_tag_table (gmp->gene, "IntMap", h) ;
	    int m1, m2, g1, g2, x1, x2 ;

	    if (imGene && imMrna)
	      {
		m1 = ac_table_int (imMrna, 0, 1, -999999) ;
		m2 = ac_table_int (imMrna, 0, 2, -999999) ;
		g1 = ac_table_int (imGene, 0, 1, -999999) ;
		g2 = ac_table_int (imGene, 0, 2, -999999) ;
		if (m2 != -999999 && g2 != -999999)
		  {
		    if (g1 < g2)
		      { x1 = m1 - g1 + 1 ; x2 = m2 - g1 + 1 ; }
		    else
		      { x1 = g1 - m1 + 1 ; x2 = g1 - m2 + 1 ; }
		    TBB (curRowClone+1+iLine, COL_CORDGENE)=vtxtPrintf (bfr, "%d to<br/>%d"ooo, x1, x2) ;
		    TBB (curRowClone+1+iLine, COL_CORDGENOME)=vtxtPrintf (bfr, "%d to<br/>%d"ooo, m1, m2) ;
		  }
	      }
	    ac_free (imMrna) ;
	    ac_free (imGene) ;	    
	  }

	  /* gap length */
	  {
	    int	gapLen = ac_tag_int (oMrna, "Gap_length", 0), iPr=0, iSto=0, iss ; 
	    TBB (curRowClone+1+iLine, COL_SEQGAP)=vtxtLen (bfr) ; 
	    
	    if (gapLen)vtxtPrintf (bfr, "gap:%.2lfkb<br/>", gapLen/1000.) ; 
	    if (gSplicing)for (iss=0 ; iss < gSplicing->rows ; iss++)
	      {
		txt = ac_table_printable (gSplicing, iss, 4, "") ; 
		if (!strcasecmp (txt, "Predicted_Exon"))iPr++ ; 
		else if (!strcasecmp (txt, "Stolen_Exon"))iSto++ ; 
	      }
	    if (iSto)vtxtPrintf (bfr, "%d exon%s inferred<br/>", iSto, _multi (iSto)) ; 
	    if (iPr)vtxtPrintf (bfr, "%d exon%s predicted<br/>", iPr, _multi (iPr)) ; 
	    
	    vtxtPrintf (bfr, ooo) ; 
	  }
	  
	  /* exon in CDS */
	  {
	    int x1 = ac_tag_int (oProduct, "Nb_introns_in_CDS", 0) + 1 ;
	    TBB (curRowClone+1+iLine, COL_EXON_CDS)=vtxtPrintf (bfr, "%d"ooo, x1) ;
	  }
	  /* protein quality */
	  {
	    int q = ac_tag_int (oProduct, "Quality", -999) ;
	    const char *ccp = 0 ;
	    if (ac_has_tag (oProduct, "Very_good_product"))
	      ccp = "Very good" ;
	    else if (ac_has_tag (oProduct, "Good_product"))
	      ccp = "Good" ;
	    else if (q > 0)
	      ccp = "Questionable" ;
	    else if (q > -999)
	      ccp = "Apparently non coding" ;

	    if (ccp || proteinLength)
	      { 
		char buf1[256], buf2[256] ;
		
		sprintf (buf1, "PEP_Product::") ;
		if (proteinLength)
		  sprintf (buf2, "%d&nbsp;aa", proteinLength) ;
		else
		  sprintf (buf2, " ");
		TBB (curRowClone+1+iLine, COL_PROTEIN_QUALITY)=vtxtPrint (bfr, "") ;
		gmpFakeObjLink (bfr, gmp, buf1, oProduct, buf2) ;
		if (ccp)
		  vtxtPrintf (bfr, "<br>%s", ccp) ;
		vtxtPrint (bfr, ooo) ; 
	      }
	  }
	  /* Domains/Motif */
	  {
	    int ii, n = 0, jr ;
	    AC_TABLE motifs = ac_tag_table (oProduct, "PFam", h) ;
	    AC_OBJ pfam = 0 ;
	    const char *ccp ;
	    char **cpp, *psNN [] =
	      { 
		"N_terminal_signal_domain", "N-terminal peptide signal",
		
		"Transmembrane_domain", "transmembrane domain", 
		"N_myristoylation_domain", "N-myristoylation domain",
		"prenylation_domain", "prenylation domain",
		"Golgi_transport_domain", "Golgi transport domain", 
		
		"peroxisomal_domain", "peroxisomal domain", 
		"2nd_peroximal_domain", "second peroximal domain", 
		"vacuolar_domain", "vacuolar domain", 
		
		"ER_retention_domain", "ER_retention domain", 
		"Coiled_coil_region", "coiled coil stretch", 
		"Leucine_zipper_domain", "leucine zipper domain",
		"actin_binding_1_domain", "actin binding domain", 
		"actin_binding_2_domain", "actin binding domain", 
		"ATP_binding_domain", "ATP binding domain",   
		
		"RNA_binding_domain", "RNA binding domain",                           
		0, 0
	      } ;
	    
	    for (jr = 0 ; motifs && jr < motifs->rows ; jr++)
	      {
		pfam = ac_table_obj (motifs, jr, 0, h) ;
		ccp = ac_tag_printable (pfam, "Definition", ac_name(pfam)) ;
		if (!n++)
		  TBB (curRowClone+1+iLine, COL_DOMAIN)=vtxtPrintf (bfr, "") ;
		else
		  vtxtPrint (bfr, ", ") ;
		vtxtPrint (bfr, ccp) ;
		ac_free (pfam) ;
	      }
	    ac_free (motifs) ;
	    {
	      int n1 = n ;
	      for (ii = 0, cpp = psNN ; *cpp ; ii++, cpp += 2)
		if (ac_has_tag (oProduct, *cpp))
		  {
		    if (!n1)
		      TBB (curRowClone+1+iLine, COL_DOMAIN)=vtxtPrintf (bfr, "") ;
		    else
		      vtxtPrint (bfr, n++ ? ", " : "<br>") ;
		    n1 = 1 ;
		    vtxtPrint (bfr, *(cpp+1)) ;
		  }
	    }
	    
	    if (n) 
	      vtxtPrint (bfr, ooo) ;
	  }

	  /* Localiszation/Psort_title */
	  {
	    int n = 0, jr ;
	    AC_TABLE locs = ac_tag_table (oProduct, "Psort_title", h) ;
	    const char *ccp ;

	    for (jr = 0 ; locs && jr < locs->rows ; jr++)
	      {
		ccp = ac_table_printable (locs, jr, 0, "") ;
		if (!n++)
		  TBB (curRowClone+1+iLine, COL_LOCAL)=vtxtPrintf (bfr, "") ;
		else
		  vtxtPrint (bfr, ", ") ;
		vtxtPrint (bfr, ccp) ;
	      }
	    ac_free (locs) ;

	    if (n) 
	      vtxtPrint (bfr, ooo) ;
	  }

	  /* protein features */
	  {
	    BOOL f1, f2, f3 ;
	    f1 = f2 = f3 = FALSE ;
	    { /* NMD */
	      f1 = ac_has_tag (oMrna, "NMD") && ac_has_tag (oProduct, "Best_product") ;
	    }
	    { /* multiple products */
	      AC_KEYSET ksf = ac_objquery_keyset (oProduct, ">mrna ; COUNT {> product; Product ; good_product} > 1", 0) ;
	      if (ksf && ac_keyset_count (ksf)) f2 = TRUE ;
	       ac_free (ksf) ;
	    }
	    { /* uORF */
	      AC_KEYSET ksf =  ac_objquery_keyset (oProduct, ">mRNA ; >Product ; uORF_candidate", 0) ;
	      if (ksf && ac_keyset_count (ksf)) f3 = TRUE ;
	      ac_free (ksf) ;
	    }
	    if (f1 || f2 || f3)
	      {
		TBB (curRowClone+1+iLine, COL_FEATURES) = vtxtMark (bfr) ;
		if (f1) 
		  {
		    vtxtPrint (bfr, "Predicted NMD candidate") ;
		    if (f2 || f3) vtxtBreak (bfr) ;
		  }
		if (f2)
		  {
		    vtxtPrint (bfr, "mRNA may encode multiple products") ;
		    if (f3) vtxtBreak (bfr) ;
		  }
		if (f3)
		  vtxtPrint (bfr, "Upstream uORF") ;
		vtxtPrint (bfr, ooo) ;
	      }
	  }
	  /* starts on , ends on */
	  {
	    if (ac_has_tag (oProduct, "Met") && ac_has_tag (oProduct, "at_position_1"))
	      {
		const char *myMet = "ATG" ;
		AC_TABLE kozak = 0 ;

		if (ac_tag_int (oProduct, "First_ATG", 0) == 1)
		  myMet = "ATG" ;
		else
		  {
		    kozak = ac_tag_table (oProduct, "First_Kozak", 0) ;
		    if (kozak && ac_table_int (kozak, 0, 0, 0) == 1)
		      myMet = ac_table_printable (kozak, 0, 1, "ATG") ;
		  }
		if (ac_has_tag (oProduct, "Down_stop"))
		  TBB (curRowClone+1+iLine, COL_FROM)=vtxtPrintf (bfr, "Met (%s)<br> to Stop"ooo, myMet) ; 
		else 
		  TBB (curRowClone+1+iLine, COL_FROM)=vtxtPrintf (bfr, "Met (%s)<br> to last codon"ooo, myMet) ;
		ac_free (kozak) ;
	      }
	    else
	      {
		if (ac_has_tag (oProduct, "Down_stop"))
		  TBB (curRowClone+1+iLine, COL_FROM)=vtxtPrintf (bfr, "1st codon<br> to Stop"ooo) ; 
		else 
		  TBB (curRowClone+1+iLine, COL_FROM)=vtxtPrintf (bfr, "1st<br> to last codon"ooo) ; 
	      }
	  }
	
	  /* cord on mrna */
	  {
	    int x1 = ac_table_int (gProduct, is, 1, 0) ;
	    int x2 = ac_table_int (gProduct, is, 2, 0) ;
	    if (x1 < 1) x1 = 1 ;
	    TBB (curRowClone+1+iLine, COL_CORDMRNA)=vtxtPrintf (bfr, "%d to %d"ooo, x1, x2) ;
	  }
	  
	  
	  /* main supporting Clones */
	  {
	    AC_TABLE gcovering_clone = ac_tag_table (oMrna, "CDS_covered_by", h) ; 
	    if (!gcovering_clone) gcovering_clone = ac_tag_table (oMrna, "mrna_covered_by", h) ;
	    if (!gcovering_clone) gcovering_clone = ac_tag_table (oMrna, "Complete_CDS_clone", h) ; 
	    if (gcovering_clone)
	      {
		if (gcovering_clone->rows > 1)
		  { int ir ;
		    TBB (curRowClone+1+iLine, COL_SUPPORTING) =
		      gmpObjLink (bfr, gmp, ac_table_obj (gcovering_clone, 0, 0, h), 0) ; 
		    for (ir = 1 ;ir < gcovering_clone->rows ; ir++)
		      {
			vtxtPrintf (bfr, "<br/> ") ;
			gmpObjLink (bfr, gmp, ac_table_obj (gcovering_clone, ir, 0, h), 0) ; 
		      }
		    if (0 && gcovering_clone->rows > 2)
		      vtxtPrintf (bfr, "<br/> and %d other clones", gcovering_clone->rows - 2) ;
		  }
		else 
		  TBB (curRowClone+1+iLine, COL_SUPPORTING) =
		    gmpObjLink (bfr, gmp, ac_table_obj (gcovering_clone, 0, 0, h), 0) ; 
		vtxtPrintf (bfr, ooo) ;
	      }
	  }
	  /* Best Available Clone */
	  {
	    AC_OBJ oBest_available_clone = ac_tag_obj (oMrna, "Best_available_clone", h) ; 
	    if (oBest_available_clone)
	      {
		TBB (curRowClone+1+iLine, COL_BESTCLONE)=vtxtPrintf (bfr, "%s"ooo, ac_name (oBest_available_clone)) ; 
	      }
	  }
	  
	  iLine++ ; 
	}
      if (addRow < iLine)addRow=iLine ; 
      curRowClone+=addRow ; 
    }

  gmp->mrna = originalMrna ;
  gmp->product = originalProduct ;
  gmp->kantor = originalKantor ;

  maxRows=1+curRowClone ; 
  
  /* Outputing now */
  /* SECTION sequences */
  if (gmp->tg)
    { 
      if (! justSequences) /* else it is a chapter by itself */
	
	{
	  gmpSection (blkp, gmp, "mRNA_sequences", "Sequences") ; 

	  vtxtPrint (blkp, "All sequences, including the mRNA, the premessenger or transcription unit, the 5 kb upstream and the UTRs are clickable and link to their DNA content, individual introns and exons sequences are given <a href=\"javascript:openAnchor(0,'tg_introns')\">below</a>") ;

	}
      fichePrintSquareTable (gmp, gmp->style, Tbb, blkp, bfr, 0, 0, 0, maxRows, 
				 COL_VARIANT, COL_VSEQ, 
				 COL_PSEQ, COL_5PUTR, COL_3PUTR, COL_UORF,
				 /* COL_AMSEQ, */
				 COL_PROMOTOR,COL_PREMESSENGER, COL_POSTMOTOR,
				 -1) ; 
      if (! justSequences) /* else it is a chapter by itself */
	{
	  if (!myMrna)
	    {
	      vtxtPrint (blkp, "For computer oriented people, we provide a composite "
			 ) ;
	      vtxtPrintf (blkp
			  , "<a href=\"javascript:openAceViewLink ('gene', '%s', 'fasta')\">FASTA sequence page</a>"
			  , ac_name(gmp->gene)) ;
	      if (0) gmpHelp (blkp, gmp, "help", "helpFasta") ;
	    }
	  else
	    {
	      BOOL hasStop = ac_has_tag (gmp->product, "Down_stop") ;
	      
	      if (gmp->product && !gmp->peptide)
		gmp->peptide = ac_obj_peptide (gmp->product, gmp->h) ;
	      if (gmp->peptide)
		ficheNewMrnaDecoratedPeptide (blkp, gmp, hasStop) ;
	      
	      if (!gmp->dna)
		gmp->dna = ac_obj_dna (gmp->mrna, gmp->h) ;
	      if (gmp->dna)
		ficheNewMrnaDnaDecoratedSequence (blkp, gmp, 0, 0) ;
	    }
	}
      
      if (justSequences) goto done ;
    }
  /* SECTION END */

  /* SECTION structure */
  if (gmp->tg)
    { 
      if (myMrna || gmp->nMrna == 1)
	gmpSection (blkp, gmp, "mRNAs", "mRNA structure and expression") ; 
      else if (gmp->nMrna > 1)
	gmpSection (blkp, gmp, "mRNAs", "mRNAs structure and expression") ; 
      if (!myMrna)
	{
	  vtxtPrint (blkp, "The AceView mRNA models are a non-redundant, comprehensive and curated representation of the cDNA sequence data in the public repositories (GenBank and dbEST). This table documents the physical properties of each representative mRNA, whose sequence has been corrected for sequencing errors and matched to the genome") ;
	  vtxtBreak (blkp) ;
	}
      fichePrintSquareTable (gmp, gmp->style, Tbb, blkp, bfr, 0, 0, 0, maxRows, 
			     COL_VARIANT, COL_EXON, COL_CLONE, COL_TISSUE,
			     COL_5PRIME, COL_SEQGAP, COL_3PRIME, 
			     COL_CORDGENE,COL_CORDGENOME, -1) ; 
      if (myMrna && ! justSequences)
	{
	  ficheNewMrna5PrimeParagraph (blkp, gmp) ;
	  ficheNewMrna3PrimeParagraph (blkp, gmp) ;
	}

      {
	vTXT bfr1 = vtxtCreate () ; 
	if (gmp->markup) vtxtMarkup (bfr1) ;

	ficheNewGenePolyATable (bfr1, gmp) ; 
	if (vtxtPtr (bfr1))
	  {
	    vtxtPrintf (blkp, "Validated poly A sites") ;
	    vtxtPrint (blkp, vtxtPtr (bfr1)) ; 
	  }
	vtxtDestroy (bfr1) ;
      }
      vtxtBreak (blkp) ;

    }
  /* SECTION END */

  
  /* SECTION Annotation of variants */
  if (gmp->view == 'g')
    {
      gmpSection (blkp, gmp, "mrnaAnnotTable", "Summaries of AceView transcripts and proteins") ;
      if (gMrna && gMrna->rows > 1)
	vtxtPrint (blkp, "Variants are named .a .b .c ... from longest to smallest ORF.") ;
      
      fichePrintSquareTable (gmp, gmp->style, Tbb, blkp, bfr, 0, 0, 0, maxRows, 
			     COL_VARIANT, COL_OVERVIEW, -1) ; 
      vtxtBreak (blkp) ;
    }
  /* SECTION END */

  if (gmp->Spc == WORM)
    {
      if (!gmp->tg && gmp->pg)
	{ 
	  AC_OBJ oProduct = ac_tag_obj (gmp->gene, "Product", h) ;
	  AC_OBJ oMrna = oProduct ? ac_tag_obj (oProduct, "mRNA", h) : 0 ;

	  if (oMrna)
	    {
	      vtxtBreak (blkp) ;
	      vtxtPrintf (blkp, "The %sCDS is fully annotated "
			  , gtPredictedMrnaSupport (gmp->pg, 0) < 4 ? "predicted " : ""
			  ) ;
	      gmpObjLink (blkp, gmp,  oMrna, "<b>here</b>") ; 
	    }
	}
    }
#ifdef JUNK
  if (gmp->Spc != WORM) /* not worm */
    {
      vtxtDot (blkp) ;
      vtxtPrintf (blkp,
		  "\nWarning: we annotate only one open reading frame (ORF) per mRNA, "
		  "choosing the longest, and deriving its sequence from " 
		  ) ;
      
      if (!ac_has_tag (gmp->gene, "use_am"))
	vtxtPrintf (blkp,
		    "the underlying genome. "
		    "If there is an error in the genome, a better ORF "
		    "may be derived from the cDNA consensus sequence. " 
		    "It is also possible "
		    );
      else
	vtxtPrintf (blkp,
		    "the cDNA consensus. " 
		    "It is possible "
		    ) ;
      
      vtxtPrintf (blkp,
		  "that the cell uses another frame, or makes more than one product per mRNA. "
		  "The ORF we annotate on each mRNA is shown as a broad "
		  "solid pink area on the drawing. An open reading frame that does not "
		  "cover most of the standard gt-ag or gc-ag intron boundaries "
		  "(both drawn in pink, blue being reserved for atypical splice sites) is in our "
		  "opinion suspicious. If you are interested in the gene, we recommend that "
		  "you reanalyse yourself all these possibilities using the sequences given "
		  "<a href=\"javascript:openAceViewLink ('%s', '%s', 'fasta')\"><b>here</b></a>, "
		  "in particular the AceView reference sequences, which represent the "
		  "consensus of cDNA sequences guided by the genome sequence"
		  , ac_class (gmp->gene), ac_name (gmp->gene)
		  ) ;
    }
#endif /* JUNK */

  /* SECTION Proteins */
  gmpSection (blkp, gmp, "Proteins", "Predicted protein properties") ; 
  if (!myMrna)
    vtxtPrint (blkp
	       , "This table allows to see at a glance from the last"
	       " column if an isoform has its exonic structure fully supported"
	       " by a single clone (the variant identifier a, b, c under such mRNA is underlined in the gene diagrams),"
	       " or if it requires concatenation"
	       " of two or more cDNA clones (identifier not underlined).") ;
  vtxtBreak (blkp) ;

  if (gmp->Spc == WORM)
    fichePrintSquareTable (gmp, gmp->style, Tbb, blkp, bfr, 0, 0, 0, maxRows, 
			   COL_PROTEIN,  COL_PROTEIN_QUALITY, COL_EXON_CDS, COL_DOMAIN, COL_LOCAL, COL_HIERARCHY, COL_FROM, COL_CORDMRNA, COL_FEATURES, COL_SUPPORTING, COL_BESTCLONE, -1) ; 
  else
    fichePrintSquareTable (gmp, gmp->style, Tbb, blkp, bfr, 0, 0, 0, maxRows, 
			   COL_PROTEIN,  COL_PROTEIN_QUALITY, COL_EXON_CDS, COL_DOMAIN,  COL_LOCAL, COL_HIERARCHY, COL_FROM, COL_CORDMRNA, COL_FEATURES, COL_SUPPORTING, -1) ; 
  vtxtBreak (blkp) ;
  /* SECTION END */

  vtxtBreak (blkp) ;

  /* free the Table for next Usage */
 done:
  arrayDestroy (Tbb) ;
  arrayDestroy (mSize) ;
  vtxtDestroy (bfr) ;
  ac_free (lProd) ; 

  ac_free (h) ;

  return 1 ; 
} /* ficheNewGeneAnnotationOfVariantsTable */

/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

void ficheNewGeneMOLECULESChapter (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  vTXT bfr = vtxtHandleCreate (h) ; 
  vTXT bfr1 = vtxtHandleCreate (h) ; 

  if (gmp->markup) vtxtMarkup (bfr) ;   
  if (gmp->markup) vtxtMarkup (bfr1) ;   
  /*  see if i should open the chapter */
  gmpChapter (bfr1, gmp,  "Gene_MOLECULAR", "MOLECULAR ANNOTATION OF ACEVIEW mRNAs") ;
  
  /* mesure sl1, protein, localization ancetres */
  /* proprietes des mrnas completeness, 5'length, 3'length, representative clones */
  /* proteins, completeness, length, pI MW, hiearchy representative clone -> best candidate clone */

  /* this stuff contains 3 section in one function,uuuurk ! */
  if (1) 
    ficheNewGeneAnnotationOfVariantsTable (bfr, gmp, FALSE) ;
  if (0)    
    ficheNewGeneCountVariantsStatement (bfr, gmp) ;

  ficheNewGeneIntronsParagraph (bfr, gmp) ;
  if (0) ficheNewGenePolyAParagraph (bfr, gmp) ;
  
  if (vtxtPtr (bfr))
    { 
      vtxtPrint (blkp, vtxtPtr (bfr1)) ; 
      vtxtPrint (blkp, vtxtPtr (bfr)) ; 
      gmpChapterClose (blkp, gmp, "Gene_MOLECULAR", TRUE) ;
    }
  else
    gmpChapterClose (bfr1, gmp, "Gene_MOLECULAR", FALSE) ;

  ac_free (h) ;
} /* ficheNewGeneMOLECULESChapter */

/***************************************************************************************/

void ficheNewGeneSequenceChapter (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  vTXT bfr = vtxtHandleCreate (h) ; 
  vTXT bfr1 = vtxtHandleCreate (h) ; 

  if (gmp->markup) vtxtMarkup (bfr) ;   
  if (gmp->markup) vtxtMarkup (bfr1) ;   
  /*  see if i should open the chapter */
  gmpChapter (bfr1, gmp,  "*Gene_sequences", "Sequences: click on the numbers to get the DNA") ;
  
  /* this stuff contains 3 sectiona in one function,uuuurk ! */
  ficheNewGeneAnnotationOfVariantsTable (bfr, gmp, TRUE) ;
 
  if (vtxtPtr (bfr))
    { 
      vtxtPrint (blkp, vtxtPtr (bfr1)) ; 
      vtxtPrint (blkp, vtxtPtr (bfr)) ; 
      gmpChapterClose (blkp, gmp, "Gene_sequences", TRUE) ;
    }
  else
    gmpChapterClose (bfr1, gmp, "Gene_sequences", FALSE) ;

  ac_free (h) ;
} /* ficheNewGeneSequenceMOLECULESChapter */

/***************************************************************************************/
/********************* Mrna Chapter ****************************************************/
/***************************************************************************************/

BOOL ficheNewMrnaSummaryChapter (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  
  /* we have to use SUMMARY as first para bacause of main.js */ 
  if (gmp->variant && gmp->gene && ac_keyset_count (ac_objquery_keyset(gmp->gene,">transcribed_gene;>mrna",h)) > 1)
    gmpChapter (blkp, gmp, "*SUMMARY", messprintf ("mRNA summary: %s, variant %s"
						   , ac_name(gmp->gene), gmp->variant)) ; 
  else
    gmpChapter (blkp, gmp, "*SUMMARY", "mRNA summary") ;

  ficheMRNAOverviewParagraphContent (blkp, gmp) ;
  vtxtBreak (blkp) ;
  if (gmp->markup)
    {
      vtxtPrint (blkp, "<span class='ace_summary'>") ; 
      vtxtPrint (blkp, "To mine knowledge about the gene, please click the ") ;
      vtxtPrintf (blkp, "<a href='javascript:openAceViewAction (\"gene\", \"%s\", \"fgene\")'>"
		  , ac_name(gmp->gene)) ;
      vtxtPrintf (blkp,"'Gene Summary'</a>") ;
      if (!ac_has_tag (gmp->gene, "Cloud_gene"))
	{
	  vtxtPrint (blkp, " or the ") ;
	  vtxtPrintf (blkp, "<a href='javascript:openAceViewAction (\"gene\", \"%s\", \"ffunc\")'>"
		      , ac_name(gmp->gene)) ;
	  vtxtPrintf (blkp,"'Function and related genes'</a>") ;
	}
      vtxtPrint (blkp, " tab at the top of the page") ;
      vtxtPrint (blkp, "</span>\n") ;
    }
  gmpChapterClose (blkp, gmp, "SUMMARY", TRUE) ;
  
  ac_free (h) ;
  return TRUE ;
} /* ficheNewMrnaSummaryChapter */

/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

void ficheNewMrnaFlashDiagramChapter (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  int dnaLn = 0, height ;
  AC_TABLE tbl = ac_tag_table (gmp->mrna, "DNA", h) ;
  dnaLn = tbl ? ac_table_int (tbl, 0, 1, 200) : 200 ;
  height = 200 + dnaLn/10 ;
  if (height > 800) height = 800 ;
  if (1)
    {
      gmpChapter (blkp, gmp, "#*mRNA_diagram", "mRNA diagram") ;
      
      vtxtPrint (blkp,
		 "The diagram shows in pink (the wide colored part is the best predicted protein"
		 ", the thin empty parts are the UTRs), from 5' (top) to 3' (bottom), the AceView "
		 "mRNA reconstructed from the GenBank cDNA sequences (in black to the "
		 "right)."
		 ) ;
      if (ac_objquery_obj (gmp->mrna, "  COUNT {>product good_product } > 0 ;>Product ; uORF_candidate", h))
	vtxtPrint (blkp,
		   "<a href ='https://www.ncbi.nlm.nih.gov/pubmed/15489325'>"
		   " 5' uORFs"
		   "</a>"
		   " are annotated in green"
		   ) ;
      vtxtPrint (blkp,
		 " Protein annotations are displayed to the left"
		 ) ;
 
      if (0) ficheListAndMarkUpMrnas (blkp, gmp, 'S', FALSE) ;
      vtxtBreak (blkp) ;

      /* the diagram proper */
      vtxtEmptyLine (blkp, 1) ;

      if (HTML5)
	{
	  int len = 0 ;
	  AC_DB db = ac_open_db ("local", 0) ;
	  char *cr, *qq = 0 ;



	  qq = hprintf (h,   "GIF ; dimensions 3000 300 ; seqget -class mrna %s -view av_mrna_whole ; seqdisplay ;  svgdump -",  ac_protect (ac_name (gmp->mrna), h)) ;
	  cr = (char *)ac_command (db, qq, &len, h) ;  
	  if (0) cr = strchr (cr,'<') ;
	  vtxtPrint (blkp, cr) ;

	  ac_db_close (db) ;
	}
      else
	{
	  vtxtPrintf (blkp,"\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
		      "  <!--\n" 
		      "  openAceViewImage (\'%s\',\'%s\',\'%s\', %d) ;\n "
		      "  //-->\n"
		      "</script>\n" 
		      , "mrna", ac_name (gmp->mrna), "vmrna"
		      , height
		      ) ;
	}
      vtxtEmptyLine (blkp, 1) ;
      /* end of the diagram proper */

      gmpCaption (blkp, gmp, "MrnaDiagramCaption", "Read more...") ;
      vtxtPrint (blkp,
		 "<font color=#007f7f >On the pink mRNA</font> are indicated, as wide boxes, the predicted protein(s), "
		 "colored according to their protein coding score (pink for good, blue for "
		 "less good, yellow for dubious),"
		 ) ;
      if (ac_objquery_obj (gmp->mrna, "  COUNT {>product good_product } > 0 ;>Product ; uORF_candidate", h))
	vtxtPrint (blkp,
		   "<a href ='https://www.ncbi.nlm.nih.gov/pubmed/15489325'>"
		   " 5' uORFs"
		   "</a>"
		   " are annotated in green"
		   ) ;
      vtxtPrint (blkp,
		 ", the untranslated regions (UTRs, thin and "
		 "empty), and the eventual validated 5' and 3' ends (flags down or up)."
		 " Introns have been spliced out, leaving triangles behind"
		 ) ;
      vtxtBreak (blkp) ;
      if (gmp->Spc == WORM)
	{
	  vtxtPrint (blkp,
		     "<font color=#007f7f >The colored curves,</font> superimposed on the diagram, indicate the number "
		     "of sequences supporting each base, as seen in RNA deep sequencing experiments. Each color "
		     "represents a different stage (for the worm) or tissue (for human, mouse..). The scale is "
		     "arbitrary, but it is often easy to spot that some exons or introns are better or less expressed "
		     "than others, or are mostly expressed in a different tissue. For worm the color code is "
		     "Embryo,L1,L2,L3,L4,Adult: light green evolving to dark blue, dauer: violet, mixed: gray."
		     ) ;
	  vtxtBreak (blkp) ;
	}
      vtxtPrint (blkp,
		 "<font color=#007f7f >To the right of the mRNA are the GenBank sequences</font> matching and defining "
		 "the transcript. Any mismatch with the local genome sequence is color "
		 "coded (red for single base variation, blue for single base insertion or "
		 "deletion). Single nucleotide polymorphisms or errors in the genome "
		 "sequence appear as a line across multiple cDNAs. The sign 0 at a 5' end "
		 "indicates a putative cap-site (start of transcription), circles at the "
		 "3' end poly A sequences. Some spurious internal priming and other "
		 "defects in the cDNAs are annotated with a red dot below the sequence "
		 "(but the flags on the mRNA pink diagram are our best guess at validated "
		 "5' and 3' ends). Sequences highlighted in turquoise are the RefSeq NM, "
		 "pink is the consensus of the cDNAs best matching the genome, pale yellow "
		 "were submitted on the opposite (probably wrong) strand"
		 ) ;
      vtxtPrint (blkp, ". More explanations are given in the ") ;
      gmpURL (blkp, gmp, "HelpmRNA.html#HelpmRNADiagram", "mRNA help file") ;
      vtxtPrint (blkp, ". Mouse over provides details, clicking even more. Scale is in bases or kilobases (k)") ;

      vtxtBreak (blkp) ;

      gmpChapterClose (blkp, gmp, "*mRNA_diagram", TRUE) ;
    }
  ac_free (h) ;
} /* ficheNewMrnaFlashDiagramChapter */

/***************************************************************************************/
/***************************************************************************************/
/**************************************************************/

static void ficheNewMrnaSupportingClonesParagraphContent (vTXT blkp, GMP *gmp, Array tissues)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_KEYSET clones = ac_objquery_keyset (gmp->mrna, ">cDNA_clone", h) ; 
  int nClo = clones ? ac_keyset_count (clones) : 0 ;
  vTXT bfr = vtxtHandleCreate (h) ; 
  char *cp ;
  if (gmp->markup) vtxtMarkup (bfr) ;

  if (nClo)
    {
      gmpSubSection (blkp, gmp, "allSupportingClones", messprintf ("%d selected cDNA clone%s support%s this mRNA"
								   , nClo, _multi(nClo), _verbs(nClo)
								   )
		     ) ;

      if (nClo > 1)
	{
	  /*  ficheCloneParagraphContent (bfr, gmp) ;
	   *  best available clone etc, probably obsolete
	   */
	  
	  vtxtPrint (blkp, "This table helps analyze the pattern of expression of the mRNA") ;
	  if (gmp->nMrna > 1) 
	    vtxtPrint (blkp, ", the tissue, cell type or disease state specificity of the alternative variants") ;
	  vtxtPrint (blkp, " and to select cDNA clones suitable for your experiments") ;
	  vtxtBreak (blkp) ;
	}
  

      switch (gmp->style)
	{
	case 'x':
	  if (0 &&
	      nClo > 24)
	    {
	      cp = hprintf (h, "<a href=\"javascript:openAceViewLink ('%s', '%s', 'clones')\"><i>more</i></a>", ac_class (gmp->mrna), ac_name (gmp->mrna)) ; 
	      ficheNewCloneTable (blkp, gmp, clones, 'q', 12, cp, 0) ;
	    }
	  else
	    ficheNewCloneTable (blkp, gmp, clones, 'q', 0, 0, tissues) ;
	  vtxtBreak (blkp) ;
	      if (vtxtPtr (bfr))
		vtxtPrint (blkp, vtxtPtr (bfr)) ; 
	      break ;
	case 'r':
	  break ;
	case 's':
	      ficheCloneLibParagraphContent (bfr, gmp) ;
	      if (vtxtPtr (bfr))
		vtxtPrint (blkp, vtxtPtr (bfr)) ; 
	      break ;
	}
      vtxtBreak (blkp) ;
    }

  ac_free (h) ;
} /* ficheNewMrnaSupportingClonesParagraphContent */

/***************************************************************************************/

BOOL ficheNewMrnaExpressionCloneSupportChapter (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  Array tissues = arrayHandleCreate (60, TISSUE, h) ;
  int nmrna, nclo ;
  vTXT bfr = vtxtHandleCreate (h) ;
  vTXT bfr1 = vtxtHandleCreate (h) ;

  if (gmp->markup) vtxtMarkup (bfr) ;
  if (gmp->markup) vtxtMarkup (bfr1) ;

  nmrna = gmp->gene ? ac_keyset_count (ac_objquery_keyset (gmp->gene, ">transcribed_gene;>mrna", h)) : 0 ;
  nclo = gmp->gene && nmrna > 1 ? ac_keyset_count (ac_objquery_keyset (gmp->gene, ">transcribed_gene;>cdna_clone", h)) : 0 ;

  gmpChapter (bfr1, gmp, "*mRNA_expression", "Expression and GenBank cDNA support") ;
  ficheNewMrnaSupportingClonesParagraphContent (bfr, gmp, tissues) ;
  vtxtPrint (blkp, vtxtPtr (bfr1)) ;
  ficheNewTissuesParagraph (blkp, gmp, tissues
			    , "Tissues where expression was observed (from cDNA clones)"
			    , hprintf (h, "GenBank cDNAs whose sequence support the %s variant%s were found in "
				       , ac_name(gmp->mrna)
				       , nmrna > 1 && nclo > 5 && arrayMax (tissues) > 5 ?
				       " (no strict specificity is implied, since some partial sequences might match other variants as well)"
				       : ""
				       )
			    ) ;
  vtxtPrint (blkp, vtxtPtr (bfr)) ;
	     
  gmpChapterClose (blkp, gmp, "*mRNA_expression", TRUE) ;

  ac_free (h) ;
  return TRUE ;
} /* ficheNewMrnaCloneChapter */

/***************************************************************************************/

BOOL ficheNewTgSupportChapter (vTXT blkp, GMP *gmp, AC_KEYSET ks1, const char *title) 
{
  AC_HANDLE h = ac_new_handle () ;
  Array tissues = arrayHandleCreate (60, TISSUE, h) ;
 vTXT bfr = vtxtHandleCreate (h) ; 
  vTXT bfr1 = vtxtHandleCreate (h) ; 
   
  if (gmp->markup) vtxtMarkup (bfr) ;   
  if (gmp->markup) vtxtMarkup (bfr1) ;   
  
  gmpChapter (bfr1, gmp, "clone_table", title) ;
  ficheNewCloneTable (bfr, gmp, ks1, 'q', 0, 0, tissues) ;
  ficheNewTissuesParagraph (bfr1, gmp, tissues
			    , "Tissues where expression was observed (from cDNA clones)"
			    , "Origin of the cDNAs, as reported in GenBank/dbEST (tissue, stage, pathological or normal) indicates that this feature is expressed in "
			    ) ;
 if (vtxtPtr (bfr))
    { 
      vtxtPrint (blkp, vtxtPtr (bfr1)) ;
      vtxtPrint (blkp, vtxtPtr (bfr)) ;
      gmpChapterClose (blkp, gmp, "clone_table", TRUE) ;
    }
  else
    gmpChapterClose (bfr1, gmp, "clone_table", FALSE) ;
  
  ac_free (h) ;
  return TRUE ;
} /* ficheNewTgSupportChapter */

/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

static BOOL ficheNewMrnaAmConsensusSentence (vTXT blkp, GMP *gmp)
{
  BOOL ok = FALSE ;
  int nerr = 0 ;

  if (ac_has_tag (gmp->mrna, "From_AM"))
    {
      nerr = ac_tag_int (gmp->mrna, "Tiling_error", 0) ;
      if (nerr) vtxtDot (blkp) ;
      switch (nerr)
	{
	case 0:
	  break ;
	case 1:
	  vtxtPrintf (blkp, "The mRNA sequence representing the cDNA clone consensus differs from the genome in one position") ;
	  break ;
	case 2:
	  vtxtPrintf (blkp, "The mRNA sequence representing the cDNA clone consensus differs from the genome in two positions") ;
	  break ;
	case 3:
	  vtxtPrintf (blkp, "The mRNA sequence representing the cDNA clone consensus differs from the genome in three positions") ;
	  break ;
	default:
	  vtxtPrintf (blkp, "The mRNA sequence representing the cDNA clone consensus differs from the genome in %d positions", nerr) ;
	  break ;
	} 
      ok = TRUE ;
    }
  else if (ac_has_tag (gmp->mrna, "RefSeqMaker") )
    {
      nerr = ac_tag_int (gmp->mrna, "Tiling_error", -1) ;
      vtxtDot (blkp) ;
      if (nerr > 0)
	{
	  vtxtPrintf (blkp, "We annotate here the sequence derived from the genome, although the best path through the available clones differs from it in %d position%s, which may be SNPs, sequencing errors, or mutations"
		      , nerr, nerr > 1 ? "s" : "") ;
	  
	  vtxtPrintf (blkp, " (link to the ") ;
	  vtxtPrintf (blkp, "<a href='javascript:openAnchor (0,\"mrna_AMsequence\")'>cDNA derived sequence consensus</a>)") ;
	}
      else if (nerr == 0)
	vtxtPrintf (blkp, "The cDNA clones consensus validates the genome sequence") ;
      ok = TRUE ;
    }

  return ok ;
} /* ficheNewMrnaAmConsensusSentence */

/***************************************************************************************/

static void ficheNewMrnaLinkToOtherProductsSentence (vTXT blkp, GMP *gmp, char type)
{
  AC_OBJ product ;
  AC_TABLE products ;
  const char *ccp ;
  int ir, nn = 0 ;
  AC_HANDLE h = ac_new_handle () ;

  if (! gmp->product || ac_has_tag (gmp->product, "Best_product"))
    {
      products = gmp->mrna ? ac_tag_table (gmp->mrna, "Product", h) : 0 ;
      for (ir = nn = 0 ; products &&  products->rows > 1 && ir < products->rows ; ir++)
	{
	  product = ac_table_obj (products, ir, 0, h) ;
	  if (ac_obj_equal (product, gmp->product))
	    continue ;
	  if (!nn++)
	    vtxtPrintf (blkp,	
			"%d minor protein%s, annotated here ("
			, products->rows - 1
			, _multi(products->rows -1)
			) ;
	  else
	    vtxtPrint (blkp, ", ") ;
	  if (ac_has_tag (product, "uORF_candidate"))
	    vtxtPrint (blkp, "uORF ") ;
	  ccp = ac_name(product) ;
	  if (strstr (ccp, ac_name(gmp->gene)))
	    ccp += strlen (ac_name(gmp->gene)) ;
	  gmpObjLink (blkp, gmp, product, ccp) ;
	}
      if (nn)
	vtxtPrint (blkp,	
		   ", in pale blue or yellow) can conceptually be deduced from this mRNA sequence "
		   "and might possibly be translated in vivo") ;
    }
  ac_free (h) ;
  if (nn)
    vtxtBreak (blkp) ;
} /* ficheNewMrnaLinkToOtherProductsSentence */

/***************************************************************************************/

BOOL ficheNewMrnaProteinChapter (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  vTXT bfr = vtxtHandleCreate (h) ;  
  vTXT bfr1 = vtxtHandleCreate (h) ;  

  if (gmp->markup) vtxtMarkup (bfr) ;
  if (gmp->markup) vtxtMarkup (bfr1) ;

  gmpChapter (bfr1, gmp, "*Protein", "Protein annotation") ;
  vtxtPrint (bfr,	
	     "Warning:  The protein sequence is predicted from the mRNA sequence; "
	     "confirmation of the protein's existence and extent will require more "
	     "protein sequence data" 
	     ) ;
  ficheNewMrnaAmConsensusSentence (bfr, gmp) ;  
  vtxtBreak (bfr) ;

  ficheNewMrnaLinkToOtherProductsSentence (bfr, gmp, 'S') ;

  ficheMRNAConceptualTranslationParagraph (bfr, gmp) ; 
  if (gmp->markup)
    fichePrimersParagraph (bfr, gmp) ; 
  fichePRODUCTPsortParagraph (bfr, gmp) ; 
  fichePRODUCTPfamParagraph (bfr, gmp) ; 
  if (gmp->markup)
    fichePRODUCTBlastPParagraph (bfr, gmp) ; 
    
  fichePRODUCTTaxTreeParagraph (bfr, gmp) ; 

  if (vtxtPtr (bfr))
    { 
      vtxtPrint (blkp, vtxtPtr (bfr1)) ;
      vtxtPrint (blkp, vtxtPtr (bfr)) ;
      gmpChapterClose (blkp, gmp, "Protein", TRUE) ;
    }
  else
    gmpChapterClose (bfr1, gmp, "Protein", FALSE) ;

  ac_free (h) ;
  return TRUE ;
} /* ficheNewMrnaProteinChapter */

/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************/
/* 5PRIME */
static void ficheNewMrna5PrimeLengthSentence (vTXT blkp, GMP *gmp, int Length_5prime_UTR)
{
  BOOL isTranspliced_to = ac_has_tag (gmp->mrna, "Transpliced_to") ; 
  
  vtxtBreak (blkp) ;
  vtxtPrint (blkp, "The 5'UTR ") ; 

  if (Length_5prime_UTR > 0)
    {
      char buf1[256], buf2[256] ;
      
      vtxtPrint (blkp, " contains ") ;
      if (!isTranspliced_to)
	vtxtPrint (blkp, "about ") ;
      
      sprintf (buf1, "DNA_mRNA:%d:%d:0", 1, Length_5prime_UTR) ;
      sprintf (buf2, "%d bp", Length_5prime_UTR) ; 
      gmpFakeObjLink (blkp, gmp, buf1, gmp->mrna, buf2) ;
    }
  else if (gmp->Spc == WORM && Length_5prime_UTR < 0 && Length_5prime_UTR > -16) 
    vtxtPrintf (blkp, " contains%s %d bp stolen from the genome since the early yk library were by construction losing an average of 16 bp on the 5\' side. "
		, ((Length_5prime_UTR && (!isTranspliced_to)) ? " about" : "")
		, Length_5prime_UTR) ;
}

/**************************/

static void ficheNewMrna5PrimeParagraphContent (vTXT blkp, GMP *gmp, int Length_5prime_UTR)
{

  int firstMet, openL, Up_stop = ac_tag_int (gmp->product, "Up_stop", 0), cntSL=0 ; 
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE gTranspliced_to = ac_tag_table (gmp->mrna, "Transpliced_to", h) ; 
  BOOL isSl2s = FALSE ;
  vTXT bfr = vtxtCreate () ; 

  if (gmp->markup) vtxtMarkup (bfr) ;

  if (gTranspliced_to)
    {
      AC_KEYSET kRead ;

      kRead = ac_objquery_keyset (gmp->mrna
			      , messprintf ("Follow cDNA_clone  ;  Follow Read  ;  From_gene == \"%s\" " 
					  , ac_name (gmp->tg))
			      , 0) ;

      if (kRead && ac_keyset_count (kRead))
	cntSL = ficheREADListTransplitionStatement 
	  (bfr, gmp, gmp->gene, "The transpliced leader", 
	   kRead, 0, 1, "or ", &isSl2s) ; 
      ac_free (kRead) ; 

      if (cntSL)
	vtxtPrintf (bfr, 
		   "%s present in the mRNA in front of the sequence. "
		   , _isare (cntSL)) ; 
    }
  
  if (ac_has_tag (gmp->mrna, "Found5p"))
    {
      if (gmp->Spc == WORM)
	{
	  if (!cntSL) vtxtPrintf (bfr, "This mRNA is not transpliced. ") ; 
	}
      else if (ac_has_tag (gmp->mrna, "Transpliced_to") && ac_has_tag (gmp->mrna, "Valid5p") &&
	       ac_keyset_count (ac_objquery_keyset (gmp->mrna, ">product ; mRNA_5p_complete && best_product && !at_position_1", h)) == 0
	       )
	vtxtPrintf (bfr, "It is defined by cap selected clones. ") ;
    }      

  if (gmp->product &&
      !ac_has_tag (gmp->product, "Complete") &&
      !ac_has_tag (gmp->product, "NH2_complete") && 
      (openL = ac_tag_int (gmp->product, "Open_length", 0)) &&
      (firstMet = ac_tag_int (gmp->product, "First_Met", 0) + ac_tag_int (gmp->product, "First_ATG", 0)))
    {
      if (3*firstMet < openL)
	vtxtPrintf (bfr, 
		    "The first AUG occurs at bp %d. ", 
		    3*firstMet) ; 
      else if (openL > 0)
	vtxtPrintf (bfr, 
		    "There is no AUG codon in this CDS. ") ;
    }	
  
  
  if (Up_stop!=0)
    {
      vtxtPrintf (bfr, 
		  "There is an in frame stop in the 5'UTR %d bp before the Met. "
		  , -Up_stop) ; 
    }
  else 
    {   
      if (!ac_has_tag (gmp->mrna, "Found5p"))
	vtxtPrintf (bfr, 
		   "The mRNA may be incomplete at the 5' end: the frame is open. " ) ; 
    }
  
  if (gmp->mrna && !gmp->dna)
    gmp->dna = ac_obj_dna (gmp->mrna, gmp->h) ;
  
  if (gmp->dna && 
      Length_5prime_UTR >= SpcI[gmp->Spc].long5P)
    {
      vtxtDot (bfr) ;
      vtxtPrintf (bfr, 
		 "This 5\'UTR %d bp is among the 5%% longest "
		 "we have seen, it may play a regulatory role. The 5' UTR contains "
		 , Length_5prime_UTR+1) ; 
      ficheCountATGCs (bfr, gmp->dna, Length_5prime_UTR) ; 
    }
  if (isSl2s)
    {
      ficheTGOperonParagraphContent (bfr, gmp, ficheOperonDistance, 1, 1, 1, TRUE) ; 
    }
  if (0)
    {
     
      AC_TABLE gPrev = gmp->tg ? ac_tag_table (gmp->tg, "Previous_gene_in_cis", h) : 0  ; 
      AC_OBJ prev  = gPrev ? ac_table_obj (gPrev, 0, 0, h) : 0 ;

      if (gPrev && gPrev->rows >= 2 && ac_table_int (gPrev, 0, 1, 10000) < 2000)
	vtxtPrintf (bfr, " This gene may be part of an operon with gene %s", ac_name (prev)) ;
      {
	vTXT bfr2 ; 
	bfr2 = vtxtCreate () ; 
	if (gmp->markup) vtxtMarkup (bfr2) ;

	ficheTGTitleStatement (bfr2, gmp, prev) ;
	if (vtxtPtr (bfr2))
	  vtxtPrintf (bfr, " (%s)", vtxtPtr (bfr)) ; 
	vtxtDestroy (bfr2) ;
      }
      vtxtPrintf (bfr, ", endding %d bp upstream.",  ac_table_int (gPrev, 0, 1, 0)) ; 
      ac_free (prev) ;
    }

  if (vtxtPtr (bfr))
    {
      vtxtPrint (blkp, vtxtPtr (bfr)) ;
    }
  vtxtDestroy (bfr) ; 
  ac_free (h) ;
} /* fichenewMrna5PrimeParagraphContent */

/* -===================================- /
* -=  3PRIME                  =- /
* -===================================- */	

static void ficheNewMrna3PrimeParagraphContent (vTXT blkp, GMP *gmp, int Length_3prime_UTR)
{
  AC_TABLE gPolyA_Signal ; 
  int ir, ln = 0, PolyA_pos ; 
  const char    *PolyA_Signal ; 
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE gTranspliced_to = ac_tag_table (gmp->mrna, "Transpliced_to", h) ; 
  AC_TABLE gProducts = ac_tag_table (gmp->mrna, "Product", h) ;
  AC_TABLE gDna = ac_tag_table (gmp->mrna, "DNA", h) ;
  AC_OBJ oProduct ;

  for (ir = 0 ; ir < gDna->rows ; ir++)
    {
      ln = ac_table_int (gDna, ir, 1, 0) ;
      break ;
    }
  
  for (ir = 0 ; ir < gProducts->rows ; ir++)
    {
      oProduct = ac_table_obj (gProducts, ir, 0, h) ;
      if (oProduct && ac_has_tag (oProduct, "Best_product"))
	{
	  Length_3prime_UTR = ln - ac_table_int (gProducts, ir, 2, 1) ;
	  break ;
	}
    }
  
  if (!ac_has_tag (gmp->product, "Complete") && !ac_has_tag (gmp->product, "COOH_Complete"))
    vtxtPrintf (blkp, " is not reached, the gene is incomplete.") ; 
  else 
    { 
      char buf1[256], buf2[256] ;
      char *mrna_dna = ac_obj_dna (gmp->mrna, h) ;
      int mrnaLen = mrna_dna ? strlen (mrna_dna) : 0 ;

      vtxtPrint (blkp, " contains ") ;
      if (!gTranspliced_to)
	vtxtPrint (blkp, "about ") ;
     
      sprintf (buf1, "DNA_mRNA:%d:%d:0", mrnaLen - Length_3prime_UTR + 1,  mrnaLen) ;
      sprintf (buf2, "%d bp", Length_3prime_UTR) ; 
      gmpFakeObjLink (blkp, gmp, buf1, gmp->mrna, buf2) ;
      vtxtPrint (blkp, " followed by the polyA") ;

      gPolyA_Signal = ac_tag_table (gmp->mrna, "PolyA_Signal", h) ; 
      if (gPolyA_Signal)
	{
	  PolyA_Signal = ac_table_printable (gPolyA_Signal, 0, 0, "") ; 
	  
	  vtxtPrintf (blkp, ". The standard AATAAA polyadenylation signal ") ; 
	  
	  if (strstr (PolyA_Signal, "AATAA"))
	    {
	      PolyA_pos = ac_table_int (gPolyA_Signal, 0, 1, 0) ;
	      vtxtPrintf (blkp, "is seen%s %d bp before the polyA", (PolyA_pos!=0 ? " about" :""), PolyA_pos ) ; 
	    }
	  else if (strstr (PolyA_Signal, "Variant"))
	    {
	      char	bfr[1024] ; 
	      
	      vtxtPrintf (blkp, "does not occur, but the variant ") ; 
	      PolyA_Signal = ac_table_printable (gPolyA_Signal, 0, 1, "") ; 
	      PolyA_pos = ac_table_int (gPolyA_Signal, 0, 2, 0) ;
	      strcpy (bfr, PolyA_Signal) ; 
	      vtextUpperCase (bfr) ; 
	      vtxtPrintf (blkp, "%s is seen%s %d bp before the polyA", 
			  bfr
			  , (PolyA_pos!=0 ? " about" :""), PolyA_pos ) ; 
	    }
	}
      else if (1 || gmp->Spc==WORM)
	vtxtPrintf (blkp, 
		    ". Neither the standard AATAAA polyadenylation signal "
		    "nor a variant is seen in the last 30 bp") ;  
      
      if (gmp->mrna && !gmp->dna)
	gmp->dna = ac_obj_dna (gmp->mrna, gmp->h) ;
      if (Length_3prime_UTR >= SpcI[gmp->Spc].long3P && gmp->dna)
	{
	  vtxtPrint (blkp, 
		     ". This 3\'UTR is among the 5% longest "
		     "we have seen, it may serve a regulatory function. It contains "
		     ) ; 
	  
	  ficheCountATGCs (blkp, gmp->dna + strlen (gmp->dna) - Length_3prime_UTR, Length_3prime_UTR) ; 
	}
    }
  ac_free (h) ;
} /* ficheNew3PrimeParagraphContent */

/***************************************************************************************/

static void ficheNewMrna5PrimeParagraph (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE gProducts = ac_tag_table (gmp->mrna, "Product", h) ;
  AC_OBJ oProduct ;
  int ir ;
  int Length_5prime_UTR = 0 ;
 
  for (ir = 0 ; ir < gProducts->rows ; ir++)
    {
      oProduct = ac_table_obj (gProducts, ir, 0, h) ;
      if (oProduct && ac_has_tag (oProduct, "Best_product"))
	{
	  Length_5prime_UTR = ac_table_int (gProducts, ir, 1, 1) - 1 ;
	  break ;
	}
    }

  if (Length_5prime_UTR >= 3)
    {  
      ficheNewMrna5PrimeLengthSentence (blkp, gmp, Length_5prime_UTR) ;
      vtxtDot (blkp) ;
      ficheNewMrna5PrimeParagraphContent (blkp, gmp, Length_5prime_UTR) ;
    }

  ac_free (h) ;
} /*  ficheNewMrna5PrimeParagraph */


/***************************************************************************************/

static void ficheNewMrna3PrimeParagraph (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  vTXT bfr= vtxtHandleCreate (h) ; 
  int Length_3prime_UTR = ac_tag_int (gmp->mrna, "Length_3prime_UTR", 0) ; 
  AC_TABLE gProducts = ac_tag_table (gmp->mrna, "Product", h) ;
  AC_TABLE gDna = ac_tag_table (gmp->mrna, "DNA", h) ;
  AC_OBJ oProduct ;
  int ir, ln = 0 ;

  for (ir = 0 ; ir < gDna->rows ; ir++)
    {
      ln = ac_table_int (gDna, ir, 1, 0) ;
      break ;
    }
  
  for (ir = 0 ; ir < gProducts->rows ; ir++)
    {
      oProduct = ac_table_obj (gProducts, ir, 0, h) ;
      if (oProduct && ac_has_tag (oProduct, "Best_product"))
	{
	  Length_3prime_UTR = ln - ac_table_int (gProducts, ir, 2, 1) ;
	  break ;
	}
    }
  

  if (gmp->markup) vtxtMarkup (bfr) ;

  ficheNewMrna3PrimeParagraphContent (bfr, gmp, Length_3prime_UTR) ;
  if (vtxtPtr (bfr))
    {
      vtxtBreak (blkp) ;
      vtxtPrint (blkp, "The 3'UTR") ;
      vtxtPrintf (blkp, " %s", vtxtPtr (bfr)) ; 
      vtxtBreak (blkp) ;
    }

  ac_free (h) ;
} /*  ficheNewMrna3PrimeParagraph */

/***************************************************************************************/
/***************************************************************************************/
#define ExonColor0 "<font color=#000000>"
#define ExonColor1 "<font color=#0000ff>"
#define ExonColor2 "<font color=#006666>"
#define ExonColor3 "<font color=#FF0000>"

/***************************************************************************************/

void ficheNewMrnaDecoratedPeptide (vTXT blkp, GMP *gmp, BOOL hasStop)
{
  AC_HANDLE h = ac_new_handle () ;
  char *pep = strnew (gmp->peptide, h) ;
  vTXT bfr = vtxtHandleCreate (h) ;
  char *cp, cc, buf[2] ;
  const char *ccq ;
  int a1, a2, a3, a11, a22, jExon, ir, jj, p1, p2 ;
  AC_TABLE gSplicing = ac_tag_table (gmp->mrna, "Splicing", h) ;
  AC_TABLE gProducts = ac_tag_table (gmp->mrna, "Product", h) ;
 
  buf[1] = 0 ;
  for (ir = p1 = p2 = 0 ; gProducts && ir < gProducts->rows ; ir++)
    if (ac_obj_key (gmp->product) == ac_table_key (gProducts, ir, 0, 0))
      { 
	p1 = ac_table_int (gProducts, ir, 1, 0) ; 
	p2 = ac_table_int (gProducts, ir, 2, 0) ; 
      }
  if (!p2) goto done ;
  if (gmp->markup) vtxtMarkup (bfr) ;
  
  cp = strchr (pep, '*') ;
  if (cp)
    cp[1] = 0 ;
  vtxtPrint (bfr, pep) ;
  vtextUpperCase (vtxtPtr (bfr)) ; /* bring to upper case */
  
  vtxtBreak (bfr) ;
  cp = pep = vtxtPtr (bfr) ;
  vtxtPrintf (blkp, ">Conceptual protein %s, %d aa, with exons in alternate blue green colors and red amino acids encoded at junctions", ac_name (gmp->product) , strlen (cp)) ;
 /* do it yourself with colors */
  {
    vtxtPrint (blkp, "<font face='courier' color=#204040>\n<pre AA_START>\n") ;
    
    for (jExon = ir = 0, jj = 150 ; gSplicing && ir < gSplicing->rows ; ir++)
      {
	a1 = ac_table_int (gSplicing, ir, 2, 0) ; /* unspliced coords */
	a2 = ac_table_int (gSplicing, ir, 3, 0) ;
	ccq = ac_table_printable (gSplicing, ir, 4, "") ;
	if (!strstr (ccq, "xon"))
	  continue ;
	if (p2 && a1 > p2)
	  break ;
	if (p2 && a2 > p2)
	  a2 = p2 ;
	jExon++ ;
	if (a2 < p1)
	  continue ;
	if ((jExon-1) % 2)
	  vtxtPrint (blkp, ExonColor1) ;
	else
	  vtxtPrint (blkp, ExonColor2) ;
	if (a1 < p1)
	  a1 = p1 ;
	cp = pep + (a1 - p1)/3 ;
	/* jj is the number of char available on the current line */
	a11 = p1 + 3 * ((int)(a1 - p1 + 2)/3) ;
	a22 = a2 - (a2 + 1 - p1)%3 ;
	for (a3 = a11 + jj ; a11 <= a22 ; a11 = a3, a3 += 150)
	  {
	    cp = pep + (a11 - p1)/3 ;
	    if (a3 > a22 + 1)
	      {
		a3 = a22 + 1 ;
		jj = a3 - a11 + 150 - jj ;
	      }
	    else
	      jj = 150 ;
	    cc = *(pep + (a3 - p1)/3) ;
	    *(pep + (a3 - p1)/3) = 0 ;
	    vtxtPrint (blkp,  cp) ;
	    *(pep + (a3 - p1)/3) = cc ;
	    if (jj < 150)
	      break ;
	    vtxtPrint (blkp, "<AA>\n") ;
	  }
	if (a22 < a2)
	  {
	    vtxtPrint (blkp, "</font>") ;
	    vtxtPrint (blkp, ExonColor3) ;
	    buf[0] = *(pep + (a2 - p1)/3) ;
	    vtxtPrint (blkp, buf) ;
	    jj += 3 ;
	  }
	jj = 150 - jj ; /* to be used when starting the next exon */
	if (jj <= 0) jj = 150 ;
	vtxtPrint (blkp, "</font>") ;
      }
    if (hasStop)
      {
	buf[0] = '*' ; buf[1] = 0 ;
	vtxtPrint (blkp, buf) ;
      }
    vtxtPrint (blkp, "</pre AA_END>\n</font>\n<br/>\n") ;
  }
 done:
  ac_free (h) ;
} /* ficheNewMrnaDecoratedPeptide */

/***************************************************************************************/

static void ficheNewMrnaPeptideParagraph (vTXT blkp, GMP *gmp)
{
  BOOL hasStop = ac_has_tag (gmp->product, "Down_stop") ;
  
  if (gmp->product && !gmp->peptide)
    gmp->peptide = ac_obj_peptide (gmp->product, gmp->h) ;
  if (gmp->peptide)
    {
      gmpSection (blkp, gmp, "mRNA_Protein_sequence"
		  , messprintf ("Protein sequence: %d aa", strlen(gmp->peptide))) ;
      if (gmp->style != 'x')
	vtxtPeptide (blkp, gmp->peptide, hasStop) ;
      else
	ficheNewMrnaDecoratedPeptide (blkp, gmp, hasStop) ;
    }
}  /* ficheNewPeptideParagraph */

/***************************************************************************************/
/***************************************************************************************/

char *ficheMrnaGenomeDNA (vTXT blkp, AC_DB db, AC_OBJ oMrna, char style, BOOL isPromotor)
{
  int dx = 100 ; /* dna line length */
  AC_HANDLE h = handleCreate () ;
  const char *ccp ;
  char *cp, cc, *sDna, *title ;
  int ir, a1, a2, a3, jExon, jj ;
  int	ma1, ma2, len, u1, u2 ;
  AC_TABLE gCovers = ac_tag_table (oMrna, "Covers", h) ; 
  AC_OBJ oCosmid = ac_tag_obj (oMrna, "Genomic_sequence", h) ;
  AC_TABLE gCoding = ac_tag_table (oMrna,"Coding", h) ;
  AC_TABLE gSplicing = ac_tag_table (oMrna,"Splicing", h) ;

  ma1 = ac_table_int (gCovers, 0, 2, 0) ;
  ma2 = ac_table_int (gCovers, 0, 4, 0) ;
  if (isPromotor)
    { u1 = ma1 - 2000 ; u2 = ma1-1 ; }
  else
    { u1 = ma1 ; u2 = ma2 ; }
  sDna = oCosmid ? ac_zone_dna (oCosmid, u1, u2, h) : 0 ;
  
  if (sDna && (len = strlen (sDna)))
    {
      GMP *gmp = gmpCreate (db, 0, 0, 0, 0, 0, style, 'z') ;
      vtextLowerCase (sDna) ;
      /* set to upper case the coding part */
      for (ir = 0 ; gCoding && ! isPromotor && ir < gCoding->rows ; ir++)
	{
	  a1 = ac_table_int (gCoding, ir, 0, 0) ; /* unspliced coords */
	  a2 = ac_table_int (gCoding, ir, 1, 0) ;
	  ccp = ac_table_printable (gCoding, ir, 4, "") ;
	  if (*ccp < 'A') continue ; 
	  cp = sDna + a1 - 2 ;
	  while (cp++ && a1++ <= a2)
	    *cp = ace_upper (*cp) ;
	}
      if (isPromotor)
	title = hprintf (h, ">Genomic promoter of mRNA %s region 2000 bp upstream of the first base<br/>", ac_name (oMrna)) ;
      else
	title = hprintf (h, ">Genomic premessenger of mRNA %s, %d&nbsp;bp, with coding in upper case and exons in alternate colors and introns in black<br/>", ac_name (oMrna) , strlen (sDna)) ;
      if (style != 'x' || isPromotor)
	{
	  vtxtSequence (blkp, sDna) ;
	}
      else /* do it yourself with colors */
	{
	  vtxtPrintf (blkp, "<br/>\n<font color=#007f7f><b>%s </b></font>", title) ;
	  if (0)
	    { /* BUG MOI est mal defini */
	      gmpObjLink (blkp, 0 /* should be gmp */, oMrna
			 ,"<img SRC='images/arrowup.gif' border='0' width='14' height='14'></a>"
			  ) ;
	    }
	  /*	    gmpHelp (blkp, gmp, jmp->file, jmp->help) ; */
	  vtxtPrint (blkp, "<br/>\n") ;
	  vtxtPrint (blkp, "<font face='courier' color=#204040>\n<pre DNA_START>\n") ;
	  for (jExon = ir = 0, jj = dx ; gSplicing && ir < gSplicing->rows ; ir++)
	    {
	      a1 = ac_table_int (gSplicing, ir, 0, 0) ; /* unspliced coords */
	      a2 = ac_table_int (gSplicing, ir, 1, 0) ;
	      ccp = ac_table_printable (gSplicing, ir, 4, "") ;
	      if (strstr (ccp, "xon"))
		{
		  if ((jExon++) % 2)
		    vtxtPrint (blkp, ExonColor1) ;
		  else
		    vtxtPrint (blkp, ExonColor2) ;
		}
	      else
		vtxtPrint (blkp, ExonColor0) ;
	      cp = sDna + a1 - 1 ;
	      for (a3 = a1 + jj ; a1 < a2 ; a1 = a3, a3 += dx)
		{
		  cp = sDna + a1 - 1 ;
		  if (a3 > a2 + 1)
		    {
		      a3 = a2 + 1 ;
		      jj = a3 - a1 + dx - jj ;
		    }
		  else
		    jj = dx ;
		  cc = *(sDna + a3 - 1) ;
		  *(sDna + a3 - 1) = 0 ;
		  vtxtPrint (blkp,  cp) ;
		  *(sDna + a3 - 1) = cc ;
		  if (jj < dx)
		    break ;
		  vtxtPrint (blkp, "<DNA>\n") ;
		}	      
	      jj = dx - jj ; /* to be used when starting the next exon */
	      if (jj <= 0) jj = dx ;
	      
	      vtxtPrint (blkp, "</font>") ;
	    }
	  vtxtPrint (blkp, "</pre DNA_END>\n<br/>\n") ;
	}
      gmpDestroy (gmp) ;
    }
  ac_free (h) ;
  return vtxtPtr (blkp) ;
} /* ficheMrnaGenomeDNA */

/***************************************************************************************/

static void ficheNewPremessengerDnaParagraph (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  int n, ma1, ma2 ;
  AC_TABLE gCovers ;
	
  gCovers = ac_tag_table (gmp->mrna, "Covers", h) ;
  ma1 = ac_table_int (gCovers, 0, 2, 0) ;
  ma2 = ac_table_int (gCovers, 0, 4, 0) ;
  
  n = ma2 > ma1 ? ma2 - ma1 + 1 : ma1 - ma2 + 1 ;
  if (n)
    gmpSection (blkp, gmp, "Premessenger_sequence", messprintf ("Premessenger sequence: %d bp", n)) ;

  ficheMrnaGenomeDNA (blkp, gmp->db, gmp->mrna, gmp->style, FALSE) ;

  ac_free (h) ;
} /* ficheNewPremessengerDnaParagraph */

/***************************************************************************************/

static void ficheNewPromotorDnaParagraph (vTXT blkp, GMP *gmp)
{
  if (gmp->mrna && ac_has_tag (gmp->mrna, "Found5p"))
    gmpSection (blkp, gmp, "Promotor_sequence", "Putative promoter region: 2kb upstream of the mRNA") ;
  else if (gmp->product && ac_has_tag (gmp->product, "Good_product") && 
	   ac_has_tag (gmp->product, "Best_product") && ac_has_tag (gmp->product, "NH2_complete"))
    gmpSection (blkp, gmp, "Promotor_sequence", "Putative promoter region: 2kb upstream of the mRNA") ;
  else
    gmpSection (blkp, gmp, "Promotor_sequence", "Genome sequence: 2kb upstream of the mRNA") ;
  ficheMrnaGenomeDNA (blkp, gmp->db, gmp->mrna, gmp->style, TRUE) ;
} /* ficheNewPromotorDnaParagraph */

/***************************************************************************************/
/***************************************************************************************/

void ficheNewMrnaDnaDecoratedSequence (vTXT blkp, GMP *gmp, int x1, int x2)
{
  int dx = 100 ; /* dna line length */
  AC_HANDLE h = ac_new_handle () ;
  char *sDna = ac_obj_dna (gmp->mrna, h) ;
  vTXT bfr = vtxtHandleCreate (h) ;
  char *cp, cc ;
  const char *ccq ;
  int a1, a2, a3, jExon, ir, jj ;
  int ln = strlen (sDna) ;
  int i5 = ac_tag_int (gmp->mrna, "Length_5prime_UTR", 0) ; 
  int i3 = ac_tag_int (gmp->mrna, "Length_3prime_UTR", 0) ; 
  AC_TABLE gSplicing = ac_tag_table (gmp->mrna, "Splicing", h) ;
  AC_TABLE gProducts = ac_tag_table (gmp->mrna, "Product", h) ;
  AC_OBJ oProduct ;

  if (gmp->markup) vtxtMarkup (bfr) ;
 
  for (ir = 0 ; ir < gProducts->rows ; ir++)
    {
      oProduct = ac_table_obj (gProducts, ir, 0, h) ;
      if (oProduct && ac_has_tag (oProduct, "Best_product"))
	{
	  i5 = ac_table_int (gProducts, ir, 1, 1) - 1 ;
	  i3 = ln - ac_table_int (gProducts, ir, 2, 1) ;
	  break ;
	}
    }
  if (x2 > ln) x2 = ln ;
  vtxtPrint (bfr, sDna) ;
  vtextUpperCase (vtxtPtr (bfr)) ; /* bring to upper case */
  /* now bring the UTR to lower case */
  cp = vtxtPtr (bfr) ;
  if (i5 > 0) while (i5--) { *cp = ace_lower (*cp) ; cp++ ; }
  cp = vtxtPtr (bfr) + strlen (vtxtPtr (bfr)) - 1 ;
  if (i3 > 0) while (i3--) { *cp = ace_lower (*cp) ; cp-- ; }
  
  
  cp = sDna = vtxtPtr (bfr) ;
  if (x2) *(cp+x2) = 0 ;
  if (x1) cp += x1 - 1 ; /* must be shifted after setting x2 */
  
  if (gmp->style != 'x')
    vtxtSequence (blkp, sDna) ;
  else /* do it yourself with colors */
    {
      if (x2)
	vtxtPrintf (blkp, ">mRNA %s, bp %d to %d with coding in upper case and exons in alternate colors<br/>", ac_name (gmp->mrna) , x1, x2) ;
      else
	vtxtPrintf (blkp, ">mRNA %s, %d bp with coding in upper case and exons in alternate colors<br/>", ac_name (gmp->mrna) , strlen (sDna)) ;

      vtxtPrint (blkp, "<font face='courier' color=#204040>\n<pre DNA_START>\n") ;
      
      for (jExon = ir = 0, jj = dx ; gSplicing && ir < gSplicing->rows ; ir++)
	{
	  a1 = ac_table_int (gSplicing, ir, 2, 0) ; /* unspliced coords */
	  a2 = ac_table_int (gSplicing, ir, 3, 0) ;
	  ccq = ac_table_printable (gSplicing, ir, 4, "") ;
	  if (!strstr (ccq, "xon") && !strstr(ccq,"Gap"))
	    continue ;
	  if (x2 && a1 > x2)
	    break ;
	  if (x2 && a2 > x2)
	    a2 = x2 ;
	  if (a2 < x1)
	    continue ;
	  if ((jExon++) % 2)
	    vtxtPrint (blkp, ExonColor1) ;
	  else
	    vtxtPrint (blkp, ExonColor2) ;
	  if (a1 < x1)
	    a1 = x1 ;
	  cp = sDna + a1 - 1 ;
	  /* jj is the number of char available on the current line */
	  for (a3 = a1 + jj ; a1 <= a2 ; a1 = a3, a3 += dx)
	    {
	      cp = sDna + a1 - 1 ;
	      if (a3 > a2 + 1)
		{
		  a3 = a2 + 1 ;
		  jj = a3 - a1 + dx - jj ; /* line length written */
		}
	      else
		jj = dx ;
	      cc = *(sDna + a3 - 1) ;
	      *(sDna + a3 - 1) = 0 ;
	      vtxtPrint (blkp,  cp) ;
	      *(sDna + a3 - 1) = cc ;
	      if (jj < dx)
		break ;
	      vtxtPrint (blkp, "<DNA>\n") ;
	    }	      
	  jj = dx - jj ; /* to be used when starting the next exon */
	  if (jj <= 0) jj = dx ;
	  vtxtPrint (blkp, "</font>") ;
	}
      vtxtPrint (blkp, "</pre DNA_END>\n</font>\n<br/>\n") ;
    }
  ac_free (h) ;
} /* ficheNewMrnaDnaDecoratedSequence */

/***************************************************************************************/

static void ficheNewMrnaDnaParagraph (vTXT blkp, GMP *gmp)
{
  if (!gmp->dna)
    gmp->dna = ac_obj_dna (gmp->mrna, gmp->h) ;

  if (gmp->dna)
    { 
       gmpSection (blkp, gmp, "mRNA_sequence"
		  , messprintf ("mRNA sequence: %d bp derived from the genome", strlen (gmp->dna))) ;
       ficheNewMrnaDnaDecoratedSequence (blkp, gmp, 0, 0) ;
    }
} /* ficheNewDnaParagraph */

/***************************************************************************************/

static void ficheNewMrnaDna_AMParagraph (vTXT blkp, GMP *gmp)
{
  char *dna ;
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ oAm = ac_tag_obj (gmp->mrna, "RefSeqMaker", h) ;
  
  if ((oAm = ac_tag_obj (gmp->mrna, "RefSeqMaker", h)) &&
      (dna = ac_obj_dna (oAm, h)))
    {
      gmpSection (blkp, gmp, "mRNA_AM_sequence"
		  , messprintf ("mRNA sequence: %d bp derived from the cDNA clones best matching the genome", strlen (dna))) ;
      
      if (!gmp->dna)
	gmp->dna = ac_obj_dna (gmp->mrna, gmp->h) ;
      if (gmp->dna && !strcasecmp (dna, gmp->dna))
	vtxtPrint (blkp, "Identical to above") ;
      else
	vtxtSequence (blkp, dna) ;
    }
  
  ac_free (h) ;
}  /* ficheNewMrnaDna_AMParagraph */

/***************************************************************************************/

BOOL ficheNewMrnaSequenceSubChapter (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  vTXT bfr = vtxtHandleCreate (h) ;  

  if (gmp->markup) vtxtMarkup (bfr) ;
  ficheNewMrnaPeptideParagraph (bfr, gmp) ; 
  ficheNewMrnaDnaParagraph (bfr, gmp) ; 
  if (0) ficheNewMrnaDna_AMParagraph (bfr, gmp) ; 
  ficheNewPromotorDnaParagraph (bfr, gmp) ;
  if (0) ficheNewPremessengerDnaParagraph (bfr, gmp) ;
  if (vtxtPtr (bfr))
    {
      vtxtPrint (blkp, vtxtPtr (bfr)) ; 
      vtxtBreak (blkp) ;
    }
  ac_free (h) ;

  return TRUE ;
} /* ficheNewMrnaSequenceSubChapter */

/***************************************************************************************/

BOOL ficheNewMrnaSequenceChapter (vTXT blkp, GMP *gmp)
{
  gmpChapter (blkp, gmp, "*mRNA_structure", "Sequences") ; 

  ficheNewMrnaSequenceSubChapter (blkp, gmp) ;
  gmpChapterClose (blkp, gmp, "mRNA_structure", TRUE) ;
  
  return TRUE ;
} /* ficheNewMrnaStructureSequenceChapter */

/***************************************************************************************/
/***************************************************************************************/
/* called from fichegraph.c */
/* style = 's', view = ''m' */
char *ficheNewMrnaSubmissionComment (vTXT blkp1, GMP *gmp)
{
  BOOL debug = FALSE ;
  
  vTXT buf2 = vtxtCreate () ;
  vTXT blkp = vtxtCreate () ;

  if (debug) vtxtPrint (blkp, " TOTO1 ") ;
  ficheNewGeneSummaryParagraph (blkp, gmp) ; /* shed + refseq summary + proteome summary */

  if (debug) vtxtPrint (blkp, " TOTO2 ") ;
  if (ficheNewGeneAliasStatement (buf2, gmp))
    {
      if (0) gmpSection (blkp, gmp, "tg_Alias_map", "Alias names") ;
      vtxtPrintf (blkp, "This %sgene %s"
		  , gtIsEssentialGene (gmp->gene) ? "essential " : ""
		  , ac_name (gmp->gene)) ;
      vtxtPrint (blkp, vtxtPtr (buf2)) ;
      ficheNewAceKogStatement (blkp, gmp, TRUE) ;
    }
  if (debug) vtxtPrint (blkp, " TOTO3 ") ;
  ficheNewGeneExpressionParagraph (blkp, gmp, 0) ; 

  if (debug) vtxtPrint (blkp, " TOTO4 ") ;
  gmpSection (blkp, gmp, "mRNA_summary", "mRNA summary") ;
  if (gmp->variant)
    vtxtPrintf (blkp, 
		" The report below describes variant %s", gmp->variant) ; 
  vtxtBreak (blkp) ;
  ficheMRNAOverviewParagraphContent (blkp, gmp) ;

  if (debug) vtxtPrint (blkp, " TOTO5 ") ;
  /*   ficheNewGeneEncodeProductPfamParagraph (blkp, gmp) ; */

  if (debug) vtxtPrint (blkp, " TOTO6 ") ;
  ficheNewGenePhenotypeParagraph (blkp, gmp) ;

  if (debug) vtxtPrint (blkp, " TOTO7 ") ;
  ficheNewGeneInteractionParagraph (blkp, gmp, 0) ;
  
  if (debug) vtxtPrint (blkp, " TOTO8 ") ;
  ficheNewGeneRegulationParagraph (blkp, gmp) ; 

  vtxtEmptyLine (blkp, 1) ;
 

  vtextCollapseSpaces (vtxtPtr (blkp)) ;
  
  vtxtPrint (blkp1, vtxtPtr (blkp)) ;
  vtxtDestroy (blkp) ;

  return vtxtPtr (blkp1) ;
} /* ficheNewMrnaSubmissionComment */

/***************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/

static BOOL ficheNewCloneTableFormat (vTXT blkp, GMP *gmp, AC_TABLE tbl, int maxLine, Array tissues)
{
  int i, ir, irg, is, a1, a2, da, x1, x2, xx ;
  float zz ;
  const char *ccp ;
  AC_HANDLE h = 0, h1 = ac_new_handle () ;
  AC_OBJ gene, mrna = 0, clone = 0, est = 0 ;
  AC_TABLE gTmp ;
  vTXT buf = vtxtHandleCreate (h1) ;
  vTXT linkBuf = vtxtHandleCreate (h1) ;
  char *notToGenbank[]={"yke", "y1L", "ycL1", "y2L", "ycL2", "y4L", "ycL4", "yc", "yd", "ykm", "ykms", "mv2h", "dvp", "Vidal", "Exelixis", "TI_", 0} ; 
  char *libNam[6], *libNam2[6] ;

  if (gmp->markup)
    vtxtMarkup (buf) ;
  libNam2 [0] = "<a href=\"http://www.kazusa.or.jp/huge/clone.req\">KIAA</a>" ;
  libNam2 [1] = "<a href=\"http://www.rzpd.de\">DKFZ</a>" ;
  libNam2 [2] = "<a href=\"http://www.nbrc.nite.go.jp/e/hflcdna2-e.html\">FLJ</a>" ;
  libNam2 [3] = "<a href=\"https://mgc.nci.nih.gov\">MGC</a>" ;
  
  libNam [0] = "KIAA" ;
  libNam [1] = "DKFZ" ;
  libNam [2] = "FLJ" ;
  libNam [3] = "MGC" ;
  libNam [4] = 0 ;

  if (maxLine == 0)
    maxLine = tbl->rows ; 
  /* otherwise use maxLine +2 so that 'more' will get trigerred */

  if (tbl) /* format the table */
    for (ir = 0 ; ir < maxLine+2 && ir < tbl->rows ; ir++)
      { 
	ac_free (h) ;
	h = ac_new_handle () ;

	/* est->genbank */
	clone = ac_table_obj (tbl, ir, 0, h) ; 
	est = ac_table_obj (tbl, ir, 1, h) ; 
	if (!clone || !est)
	  continue ;
	/* 9:Anomalies voir code + internal priming + internal priming on A rich genome 12A/13
	 * manque inverted internal priming + internal priming on A rich genome 12A/13"
	 */
	vtxtClear (buf) ; 
	ficheNewCloneAnomalyStatement (buf, clone, 0, FALSE) ;
	irg = 0 ;
	if (gmp->tg)
	  {
	    AC_KEYSET ksg = ac_objquery_keyset (clone
						, messprintf (">from_gene ;>gene ;; ! IS \"%s\"", ac_name (gmp->gene))
						, h) ; 
	    AC_TABLE tblg = ksg ? ac_keyset_table (ksg, 0, -1, 0, h) : 0 ;

	    for (irg = 0 ; tblg && irg < tblg->rows ; irg++)
	      {
		vtxtComma (buf) ;
		vtxtPrintf (buf, "%s"
			    , irg ? "" : "also hits gene "
			    ) ;
		gmpObjLink (buf, gmp, ac_table_obj (tblg, irg, 0, h), 0) ;
	      }
	  }
	if (vtxtPtr (buf))
	  ac_table_insert_text (tbl, ir, 14, vtxtPtr (buf)) ;
	else
	  ac_table_insert_text (tbl, ir, 14, 0) ;

	vtxtClear (buf) ; 
	/* check if read has an external or internal hot link */
	for (is=0 ; notToGenbank[is] ; is++)
	  {
	    if (strstr (ac_tag_printable (clone, "Library", ac_name(clone)), notToGenbank[is]))
	      break ;
	  } 
	if (ac_has_tag (est, "Composite") || notToGenbank[is])
	  {
	    AC_OBJ oDna = ac_tag_obj (est, "DNA", h) ;
	    if (oDna)
	      gmpFakeObjLink (buf, gmp, "DNA::", est, ac_name(est)) ;
	    else
	      vtxtPrint (buf, ac_name (est)) ;
	  }
	else
	  {
	    const char *cleanNam = ac_name (est) ;
	    if (!strncmp(cleanNam, "GB:",3))
	      cleanNam += 3 ;
	    vtxtClear (linkBuf) ; 
	    vtxtPrintf (linkBuf, GENBANK_LINK, cleanNam) ; 
	    gmpURL (buf, gmp, vtxtPtr (linkBuf), ac_name (est)) ;
	  }
	if (irg)
	  {
	    vtxtBreak (buf) ;
	    gmpObjLink (buf, gmp, clone, "<i>matches multiple genes</i>") ;
	  }
	if (vtxtPtr (buf))
	  ac_table_insert_text (tbl, ir, 1, vtxtPtr (buf)) ;

	/* 2 tissue, DATA;DATA not available in aql but exists in bql implemented as obj=>tag */
	vtxtClear (buf) ; 
	{
	  AC_KEYSET ksg = ac_objquery_keyset (clone
					      , ">read ; > tissue"
					      , h) ; 
	  AC_TABLE tblg = ksg ? ac_keyset_table (ksg, 0, -1, 0, h) : 0 ;
	  int nTissue = tissues ? arrayMax(tissues) : 0 ;
	  TISSUE *tt ;

	  a1 =  ac_table_int (tbl, ir, 4, 1000) ;  
	  for (irg = 0 ; tblg && irg < tblg->rows ; irg++)
	    {
	      if (tissues && !ficheNewGeneExpressionTissueIsBaddy (ac_table_printable (tblg, irg, 0, "")))
		{
		  tt = arrayp (tissues, nTissue++, TISSUE) ;
		  tt->tissue = ac_table_key (tblg, irg, 0, 0) ;
		  tt->nn = 1 ;
		}
	      vtxtPrintf (buf, "%s"
			  , irg ? ", " : ""
			  ) ;
	      if (a1 < 150)
		vtxtPrintf (buf, "<font color='red'>%s</font>", ac_table_printable (tblg, irg, 0, "")) ;
	      else
		vtxtPrintf (buf, "%s", ac_table_printable (tblg, irg, 0, "")) ;
	    }
	  }
	if (vtxtPtr (buf))
	  ac_table_insert_text (tbl, ir, 2, vtxtPtr (buf)) ;
	else
	  ac_table_insert_text (tbl, ir, 2, "") ;

	/* "3: mRNA */
	mrna = ac_table_obj (tbl, ir, 3, h) ; 
	vtxtClear (buf) ; 
	if (mrna)
	  {
	    if (gmp->tg)
	      gmpObjLink (buf, gmp, mrna, gtMrnaSuffix(ac_name(gmp->tg), ac_name(mrna), h)) ;
	    else
	      gmpObjLink (buf, gmp, mrna, ac_name(mrna)) ;
	  }
	else
	  {
	    gene = ac_table_obj (tbl, ir, 18, h) ;
	    if (gene)
	      {
		vtxtPrint (buf, "Gene ") ;
		gmpObjLink (buf, gmp, gene, ac_name(gene)) ;
		vtxtPrintf (buf, ", variant not shown") ;
	      }
	    else
	      vtxtPrintf (buf, "not aligned") ;
	  }
	ac_table_insert_text (tbl, ir, 3, vtxtPtr (buf)) ;
	
	if (mrna)
	  {
	    /* "4:From bp to bp in mRNA" */
	    a1 =  ac_table_int (tbl, ir, 4, 0) ;  
	    a2 =  ac_table_int (tbl, ir, 5, 0) ;
	    da = 0 ;
	    if (a1 < 1 && a1 > -8) { da = 1 - a1 ; a1 = 1 ; } /* be pudibond */
	    vtxtClear (buf) ;
	    vtxtPrintf (buf, "%d to %d",  a1, a2) ;
	    ac_table_insert_text (tbl, ir, 5, vtxtPtr (buf)) ;
	  
	    /* "5:From bp to bp in accession" */
	    x1 =  ac_table_int (tbl, ir, 6, 0) ;  
	    x2 =  ac_table_int (tbl, ir, 7, 0) ;  
	    if (x1 < x2) x1 += da ;
	    else x2 -= da ;	      
	    vtxtClear (buf) ;
	    if (gmp->Spc != WORM && x1 > 30 && x1 < x2) 
	      vtxtPrintf (buf, "<font color='red'>%d</font>",  x1) ;
	    else
	      vtxtPrintf (buf, "%d",  x1) ;
	    if (gmp->Spc != WORM && x2 > 30 && x1 > x2) 
	      vtxtPrintf (buf, " to <font color='red'>%d</font>",  x2) ;
	    else
	      vtxtPrintf (buf, " to %d", x2) ;
	    ac_table_insert_text (tbl, ir, 7, vtxtPtr (buf)) ;
	  
	    /* "6:ali->Accession match over (% length) */
	    x1 =  ac_table_int (tbl, ir, 8, 0) ;
	    x2 =  ac_table_int (tbl, ir, 9, 0) ;  
	    if (x1 == 0) x1 = 1 ;
	    xx = 100*x2/x1 ;
	    vtxtClear (buf) ;
	    vtxtPrintf (buf, "%d/%d<br>(%d %%)",  x2, x1, xx) ;
	    ac_table_insert_text (tbl, ir, 9, vtxtPtr (buf)) ;
	
	    /* "7:error->Base differences relative to genome (% id) // 23 diff<br>(%.1f %%id)" */
	    x2 =  ac_table_int (tbl, ir, 10, 0) ;  
	    zz = 1000 - 1000*x2/x1 ; zz /= 10.0 ;
	    vtxtClear (buf) ;
	    if (x2 == 0)
	      vtxtPrintf (buf, " %d diff<br>(100 %%id)",  x2) ;
	    else
	      vtxtPrintf (buf, " %d diff<br>(%.1f %%id)",  x2, zz) ;
	    ac_table_insert_text (tbl, ir, 10, vtxtPtr (buf)) ;
	  }
	
	/* 8:Features (one feature per line)
	 * refseq  / tiling clone 
	 * mapped in more than one gene
	 * SL1
	 * orfeome tags
	 * protein coding modifs: exactly encodes the predicted protein OR encodes the 
	                          predicted protein with amino acid variations [list ...]"
	 *
	 : fully sequenced + mapped in more than one gene + tiling clone + refseq 
	 + exactly encodes the predicted protein OR encodes the predicted protein
	 with amino acid variations [list ...]"
	*/

	vtxtClear (buf) ;
	{
	  struct {char *tag, * whatToSay ; int showData ; } prop []={
	    /* {"Best__of", "recommended", 1},  */
	    /* 	  {"Specific__of", "specific", 1},  */
	    {"Complete_CDS__of", "covers entire CDS", 2}, 
	    {"Fully_sequenced", "fully sequenced", 0}, 
	    {"Resequence", "to resequence", 0}, 
	    {0, 0}} ; 
	  
	  if (ac_table_printable (tbl, ir, 15, 0))
	    {
	      vtxtComma (buf) ;
	      vtxtPrintf (buf, "RefSeq") ;  
	    }

	  if (ac_table_printable (tbl, ir, 16, 0))
	    {
	      vtxtComma (buf) ;
	      vtxtPrintf (buf, "tiling clone") ;
	    }
	  if (gmp->Spc == HUMAN )
	    {
	      if ((ccp = ac_table_printable (tbl, ir, 17, 0)))
		{
		  for (i = 0 ; libNam[i] ; i++)
		    if (!strcasecmp (ccp, libNam[i]))
		      {
			vtxtComma (buf) ;
			vtxtPrintf (buf, "available from %s", libNam2[i]) ;
			break ;
		      }
		}
	    }


	  for (is = 0 ; prop[is].tag ; is++)
	    {
	      switch (prop[is].showData)
		{
		case 0:
		  if (ac_has_tag (clone, prop[is].tag))
		    {
		      vtxtComma (buf) ;
		      vtxtPrintf (buf, prop[is].whatToSay) ;  
		    }
		  break ;
		case 1:
		case 2:
		  if ((gTmp = ac_tag_table (clone, prop[is].tag, h)))
		    {
		      vtxtComma (buf) ;
		      vtxtPrintf (buf, prop[is].whatToSay) ; 
		      ccp = strrchr (ac_table_printable (gTmp, 0, 0, "."), '.') ; 
		      if (!ccp)
			ccp = "" ; 
		      else 
			ccp++ ; 
		      if (ccp[0])
			{
			  vtxtPrintf (buf, " of variant .") ; 
			  gmpObjLink (buf, gmp, ac_table_obj (gTmp, 0, 0, h), ccp) ; 
			}
		    }
		  break ;
		}
	    }
	}
	/* Transplicing */
	if ((ccp = ac_tag_printable (est, "Transpliced_to", "")))
	  {
	    vtxtComma (buf) ;
	    if (!strcmp (ccp, "SL0")) 
	      vtxtPrint (buf, "capped") ; 
	    else if (!strncmp (ccp, "SL", 2))
	      vtxtPrint (buf, ccp) ;
	  }
	if (ac_has_tag (est, "PolyA_after_base"))
	  { 
	    vtxtComma (buf) ;
	    vtxtPrint (buf, "AAA ") ; 
	  }

	/* export the features */
	if (vtxtPtr (buf))
	  ac_table_insert_text (tbl, ir, 13, vtxtPtr (buf)) ;
	else
	  ac_table_insert_text (tbl, ir, 13, "") ;

	/* 10: variations */
	{ 
	  /* prod2 in r->covers_product, var in prod2[type] , */
	  AC_OBJ prod = ac_table_obj (tbl, ir, 11, h) ;
	  KEY estKey = ac_obj_key (est) ;
	  int jr ;
	  
	  ac_table_insert_text (tbl, ir, 12, 0) ;

	  gTmp = ac_tag_table (prod, "covered_by", h) ;
	  for (jr = 0 ; gTmp && jr < gTmp->rows ; jr++)
	    {
	      if (estKey != ac_table_key (gTmp, jr, 0, 0))
		continue ;
	      if (ac_table_int (gTmp, jr, 2, -1) == 0 && ac_table_int (gTmp, jr, 4, -1) == 0)
		ac_table_insert_text (tbl, ir, 12, "exact") ;
	      else if ((ccp =  ac_table_printable (gTmp, jr, 5, 0)))
		ac_table_insert_text (tbl, ir, 12, ccp) ;
	    }
	}
      }
  ac_free (h) ;
  ac_free (h1) ;
  return TRUE ;
} /* ficheNewCloneTableFormat */

/*******************/

int ficheNewCloneTable (vTXT blkp, GMP *gmp, AC_KEYSET clones, char orderBy, int maxLine, const char *more, Array tissues)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tbl = 0, tbl2 = 0 ;
  int colorControl = 1 ;
  const char *errorMessage = 0 ;
  /* keep in synch column enum and clumn names */
  const char *myOrder = 0 ; /* set below in the switch */
  int cols[60] ;
  const char *trueTitles [60] ;
  /* + color control column
   *  - columns are not visible
   * 1:  means see it as column 1
   */
  const char *titles[]={ "-:Clone"
			 , "1:cDNA accession<br><i>Links to the sequence</i>"
			 , "2:Tissue->Tissue<br><br><i>most 5' clones<br>in <font color='red'>red</font>"
			 , "+3:Match<br>mRNA"
			 , "-:a1"
			 , "4:From bp<br>to bp<br>in mRNA"
			 , "-:x2->x1"
			 , "5:From bp<br>to bp<br>in accs."
			 , "-:len to be aligned"
			 , "9:ali->Accession<br>match over<br>(% length)" /* 345 bp<br>(%d %%) */
			 , "10:error->Base differences<br> relative to genome<br>(% identity)" /* 23 diff<br>(%.1f %%id) */
			 , "-:covers_product"
			 , "6:Variation->Clone encodes<br> complete protein<br>(with AA variation)"
			 , "7:Inverted->Features"
			 , "8:Anomalies->Anomalies<br>detected by<br>AceView"
			 , "-:Ref_seq"
			 , "-:tiling_clone"
			 , "-:hinv_libs"
			 , "-:gene"
			 , "-:other_gene"
			 , 0} ;
  vTXT bqlQ = vtxtHandleCreate (h) ;
  vTXT bqlQ2 = vtxtHandleCreate (h) ;
  
  /* if the bqlQ or bqlQ2 queries are  modified, their mapping jj[] below must be modified */
  vtxtPrint (bqlQ,  "select c,r,c,m,a1,a2, x1,x2,len,ali,err,prod,c,inv,ano,ref_seq, tiling, hlib, gene") ;
  vtxtPrint (bqlQ,  " from c in @active:1, m in c->in_mrna where exists m, tg in m->from_gene, gene in tg->gene ") ;
  if (gmp->gene)
    vtxtPrintf (bqlQ,  " where gene == %s ",  ac_protect (ac_name (gmp->gene), h)) ;
  vtxtPrint (bqlQ,  ", r in c->read, t in r->tissue, a1 in m->constructed_from, a2 in a1[1], r1 in a2[1] where r1 == r, x1 in r1[1], x2 in x1[1], tg1 in r->from_gene where tg1 == tg, len in tg1[1], ali in len[2], err in ali[1], inv in r#inverted , ano in c#anomalous_clone, prod in m->product where (not prod or prod#best_product)") ;
  vtxtPrint (bqlQ,  ",  ref_seq in r#ref_seq, tiling in r#mRNA_tiling, hlib in c->hinv_libs ") ;
 
  /* prod2 in r->covers_product, var in prod2[type] , */

  if (1)
    {
      vtxtPrint (bqlQ2,  "select c,r,c,inv,ano,ref_seq, hlib, gene") ;
      vtxtPrint (bqlQ2,  " from c in @active:1 where not c#in_mrna, r in c->read, t in r->tissue, inv in r#inverted , ano in c#anomalous_clone,  ref_seq in r#ref_seq,  hlib in c->hinv_libs, tg in c->from_gene, gene in tg->gene") ;
    }
  else
    {
      vtxtPrint (bqlQ2,  "select c,r,c,inv,ano,ref_seq, hlib, gene") ;
      vtxtPrint (bqlQ2,  " from c in @active:1 where not c#in_mrna, r in c->read, t in r->tissue, inv in r#inverted , ano in c#anomalous_clone,  ref_seq in r#ref_seq[0],  hlib in c->hinv_libs, tg in c->from_gene, gene in tg->gene") ;
    }
  
  /* auto configure the titles and visible columns using the -: and 4: constructions 
     *  -> separate the table-makewr title from the web title
     */  
  colorControl = ficheTableConfigureColumns (titles, trueTitles, cols) ;
  
  if (gmp->Spc == WORM)
    titles[1] = "Sequence" ;
  
  /*
    accession// tissue// match over #bp<br>(%len) // nb diff (%identite) // gene// mrna// 
    coords on reasd// coords on mrna 
    // properties (tilingclon e or refseq or 
    fully sequences
    error on product
    // anomalies 
    */
  
  /* select line ordering according to user choice */
  switch ((int)orderBy)
    {
    default:
      myOrder = "+4+5+1+2" ; 
      break ;
    }
  
  /* contruct the table bql syntax */
  if (0) printf ("bql -active %s\n", vtxtPtr(bqlQ)) ;
  tbl = ac_bql_table (gmp->db, vtxtPtr(bqlQ), clones, myOrder, &errorMessage, h) ;
  tbl2 = ac_bql_table (gmp->db, vtxtPtr(bqlQ2), clones, 0, &errorMessage, h) ;
  
   if (tbl2) 
    {  /* cumulate the 2 tables */
      int i, ir, jr ;
      AC_OBJ obj ;
      const char *ccp ;
      int jj[] = { 0, 1, 2, 13, 14, 15, 17, 18 } ;

      if (! tbl)
	{
	  tbl = tbl2 ; tbl2 = 0 ;
	}
      if (tbl2)
	for (ir = 0, jr = tbl->rows  ; ir < tbl2->rows ; jr++, ir++)
	  {
	    /* clone read */
	    for (i = 0 ; i < 2 ; i++)
	      {
		obj = ac_table_obj (tbl2, ir, i, h) ; 
		if (obj) 
		  ac_table_insert_type (tbl, jr, jj[i], &obj, ac_type_obj) ;
	      }
	    /* t, inv,ano, ref_seq, hlib */
	    for (i = 2 ; i < 7 ; i++)
	      {
		ccp = ac_table_printable (tbl2, ir, i, 0) ; 
		if (ccp)
		  ac_table_insert_text (tbl, jr, jj[i], ccp) ;
	      }
	    /* gene */
	    for (i = 7 ; i < 8 ; i++)
	      {
		obj = ac_table_obj (tbl2, ir, i, h) ; 
		if (obj) 
		  ac_table_insert_type (tbl, jr, jj[i], &obj, ac_type_obj) ;
	      }
	  }
    }      

  /* format the table and add the http links */
  if (tbl)
    ficheNewCloneTableFormat (blkp, gmp, tbl, maxLine, tissues) ;

  /* export */
  if (tbl && tbl->rows && tbl->cols < 18)
    ac_table_insert_text (tbl, 0, 17, "toto") ;
  if (tbl && tbl->cols)
    ac_table_display (blkp 
		      , tbl, trueTitles
		      , cols, colorControl
		      , 0, maxLine, more
		      , 0
		      ) ;
  else if (errorMessage)
    vtxtPrint (blkp, errorMessage) ;

  ac_free (h) ;    
  return 1 ; 
} /* ficheNewCloneTable */

/**************************************************************************************/

static void ficheNewTissuesSentence (vTXT blkp, GMP *gmp, Array tissues, const char *subtitle)
{
  TISSUE *tt ;
  int i ;

  /* tisues of these clones */
  if (arrayMax (tissues))
    {
      arraySort (tissues, tissueOrder) ;
      tissueCompress (tissues) ;
      arraySort (tissues, tissueNnOrder) ;
      vtxtPrint (blkp, subtitle) ;

      for (i = 0 ; i < arrayMax (tissues) ; i++)
	{
	  tt = arrp (tissues, i, TISSUE) ;
	  if (tt->nn)
	    vtxtPrintf (blkp 
			, "%s%s (%s%s%s)"
			, i ? ", " : ""
			, gtLowerCleanUp (name (tt->tissue))
			, i ? "" : "seen "
			, tt->nn > 1 ? isOne (tt->nn) : "once"
			, i ? "" : (tt->nn > 1 ? " times" : "") 
			) ;
	  else
	    break ;
	}  
      vtxtBreak (blkp) ;
    }
} /* ficheNewTissuesSentence */

/**************************************************************************************/

static void ficheNewTissuesParagraph (vTXT blkp, GMP *gmp, Array tissues, const char *title, const char *subtitle)
{
  /* tisues of these clones */
  if (arrayMax (tissues))
    {
      gmpSubSection (blkp, gmp, "Tissues", title) ;
      ficheNewTissuesSentence (blkp, gmp, tissues, subtitle) ;
    }
} /* ficheNewTissuesParagraph */

/**************************************************************************************/

void ficheNewCloneParagraph (vTXT vtxt, AC_DB db, AC_KEYSET clones, char style, char orderBy, Array tissues)
{
  GMP *gmp = gmpCreate (db, 0, 0, 0, 0, 0, style, 'z') ;
  vTXT blkp = vtxtCreate () ; 
  if (gmp->markup) vtxtMarkup (blkp) ;

  ficheNewCloneTable (blkp, gmp, clones, orderBy, 0, 0, tissues) ;
  if (vtxtPtr (blkp))
    {
      if (0) gmpSection (vtxt, gmp, "cloneListTable", "Clone Description") ; 
      vtxtPrintf (vtxt, "%s", vtxtPtr (blkp)) ; 
      vtxtBreak (vtxt) ;
    }
  vtxtDestroy (blkp) ;
  gmpDestroy (gmp) ;
} /* ficheNewCloneParagraph */

/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
