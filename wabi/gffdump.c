#include "../wac/ac.h"
#include "../wh/dict.h"
#include "../wh/vtxt.h"
#include "../wh/regular.h"
#include "../wh/regular.h"
#include "../wh/utils.h"
#include <stdio.h>
#include "../wh/dna.h"

/* Dumps transcripts from an AceDB db conforming to new 'mrna' schema as GFF (GTF) to standard out.
*/


/* identify the best product */
static AC_OBJ gffBestProduct (AC_OBJ mrna, BOOL *startp, BOOL *stopp, BOOL *cdsp
			    , int *p1p, int *p2p, int otherVeryGoodProduct, AC_HANDLE h)
{
  int i, nvg = 0, ok ;
  AC_OBJ oProduct, bestProduct = 0 ;
  AC_TABLE gProduct ;

  gProduct = ac_tag_table (mrna, "Product", h) ;
  if (gProduct)
    {
      for (i = 0; gProduct && i < gProduct->rows; i++) 
	{
	  ok = 0 ;
	  oProduct = ac_table_obj (gProduct, i, 0, h) ;
	  if (otherVeryGoodProduct && ! ac_has_tag (oProduct, "Best_product") && 
	      ac_has_tag (oProduct, "Very_good_product"))
	    nvg++ ;
	  if (otherVeryGoodProduct > 0 && otherVeryGoodProduct == nvg)
	    { ok = TRUE ; nvg++ ; }
	  
	  if (!otherVeryGoodProduct && ac_has_tag (oProduct, "Best_product") && 
	      ac_has_tag (oProduct, "Good_product"))
	    ok = TRUE ;
	  if (ok)
	    
	    {
	      bestProduct = oProduct ;
	      *p1p =  ac_table_int (gProduct, i, 1, -1) ;
	      *p2p =  ac_table_int (gProduct, i, 2, -1) ;
	      *cdsp = TRUE ;
	      *startp = ac_has_tag (bestProduct, "NH2_complete");
	      *stopp = ac_has_tag (bestProduct, "Down_stop");
	      break ;
	    }
	}
   }
  return  bestProduct ;
} /* gffBestProduct */

/**************************************************************************/

static int gffDumpMrna (AC_OBJ mrna, BOOL isQual, int otherVeryGoodProduct, int debug) 
{
  BOOL isGencode = FALSE ; /* dedicated dump for Gencode, april 2005 */
  BOOL EXP = isGencode ; /* to add names long */
  enum { UNKNOWN, EXON, INTRON, GAP, ORFGAP } tag ;
  int x1, x2, a1, a2, c1, c2, p1, p2, y1, y2, type ;
  char map0[1000] ;
  int ir, jr, ok = 0 ;
  int prodQual ;
  BOOL start, stop, cds, mrnaOk ;
  AC_TABLE gMap, gSplicing, gCov ;
  AC_OBJ geneBox, tg,  bestProduct ;
  const char *cp ;
  const char *map, *gene_id, *transcript_id, *product_id ;
  AC_HANDLE h = ac_new_handle () ;
  DICT *dict = dictHandleCreate (100, h) ;
  vTXT gNam = vtxtHandleCreate (h) ;
  vTXT tNam = vtxtHandleCreate (h) ;
  vTXT mNam = vtxtHandleCreate (h) ;
  vTXT pNam = vtxtHandleCreate (h) ;

  start = stop = cds = FALSE ; p1 = p2 = -1 ;
  bestProduct = gffBestProduct (mrna, &start, &stop, &cds, &p1, &p2, otherVeryGoodProduct, h) ;
  if (otherVeryGoodProduct && !bestProduct)
    goto done ;
    
  dictAdd (dict, "toto", 0) ;
  transcript_id = ac_name(mrna);
  product_id = bestProduct ? ac_name(bestProduct) : 0 ;
  tg = ac_tag_obj (mrna, "From_gene", h) ;
  if (!tg) goto done ;
  geneBox = ac_tag_obj (tg, "Gene", h) ;
  if (!geneBox) goto done ;
  gene_id = ac_name (geneBox) ;
  
  gMap = ac_tag_table (mrna, "IntMap", h);
  cp = ac_table_printable (gMap, 0, 0, "XXX") ;
  strcpy (map0, cp) ; map = map0 ;
  
  if (1)
    {
      const char *ccp ;
      
      if (!strncmp(transcript_id,"G_t",3))
	{
	  ccp = transcript_id + 4 ;
	  if (*ccp == '_' || *(++ccp) == '_')
	    transcript_id = ccp + 1 ;
	}
      if (product_id && !strncmp(product_id,"G_t",3))
	{
	  ccp = product_id + 4 ;
	  if (*ccp == '_' || *(++ccp) == '_')
	    product_id = ccp + 1 ;
	}
      if (!strncmp(gene_id,"G_t",3))
	{
	  ccp = gene_id + 4 ;
	  if (*ccp == '_' || *(++ccp) == '_')
	    gene_id = ccp + 1 ;
	}
      if (!strncmp(map,"c_",2))
	{
	  map += 2 ;
	}
    }

  if (isGencode &&
      ! ac_has_tag (mrna, "gt_ag") &&
      ! ac_has_tag (mrna, "gc_ag"))
  {
    AC_KEYSET ks1 = ac_objquery_keyset (mrna, ">cdna_clone ", h) ;
    AC_KEYSET ks2 = ac_ksquery_keyset (ks1, ">read ; from_gene:6 < 3", h) ;
    AC_KEYSET ks3 = ac_ksquery_keyset (ks1, ">read ; from_gene:6 < 4", h) ;
    int y1, y2, y3 ;

    y1 = ac_keyset_count (ks1) ;
    y2 = ac_keyset_count (ks2) ;
    y3 = ac_keyset_count (ks3) ;

    if (y1 > 6 && 2*y2 < y1)
      printf ("PSEUDO1 %d %d %s\n", y1, y2, ac_name (mrna)) ;
    if (y1 > 6 && 10*y3 < 8*y1)
      printf ("PSEUDO2 %d %d %s\n", y1, y2, ac_name (mrna)) ;
  }
  
  if (isGencode)
    {  
      vtxtPrintf (gNam, "Acv_%s", gene_id + 2) ;
      vtxtPrintf (mNam, "Acv_%s", transcript_id + 2) ;
      if (product_id)
	vtxtPrintf (pNam, "Acv_%s", product_id + 2) ;
    }
  else
    {
      vtxtPrintf (gNam, "%s", gene_id) ;
      vtxtPrintf (mNam, "%s", transcript_id) ;
      if (product_id)
	vtxtPrintf (pNam, "%s", product_id) ;
    }

  /*
    if (ac_has_tag (geneBox, "Spliced_gene"))
    vtxtPrint (gNam, "Spliced_gene") ;
    else if (ac_has_tag (geneBox, "Single_exon_gene"))
    vtxtPrint (gNam, "Single_exon_gene") ;
  */

  if (ac_has_tag (geneBox, "Cloud") || ac_has_tag (geneBox, "Cloud_gene"))
    vtxtPrint (tNam, "Cloud") ;
  else
    vtxtPrint (tNam, "cDNA_supported") ;

  y1 = ac_table_int (gMap, 0, 1, 1) ;
  y2 = ac_table_int (gMap, 0, 2, 1) ;
  if (debug) printf ("y1 = %d   y2 = %d\n", y1, y2);



  /* how good is the transcript */
  
  if (product_id)
    {
      if (ac_has_tag (bestProduct, "Complete"))
	{ prodQual  = 2 ; if (EXP) vtxtPrintf (pNam, "_Complete") ; }
      else if (ac_has_tag (bestProduct, "NH2_Complete"))
	{ prodQual  = 4 ; if (EXP) vtxtPrintf (pNam, "_NH2_Complete") ; }
      else if (ac_has_tag (bestProduct, "COOH_complete"))
	{ prodQual  = 6 ; if (EXP) vtxtPrintf (pNam, "_COOH_complete") ; }
      else
	{ prodQual  = 8 ; if (EXP) vtxtPrintf (pNam, "_Partial") ; }
      
      gCov = ac_tag_table (bestProduct, "Covered_by", h);
      for (ir = 0 ; gCov && ir < gCov->rows ; ir++)
	{
	  if (ac_table_int (gCov, ir, 2, -1) == 0 &&
	      ac_table_int (gCov, ir, 4, -1) == 0)
	    {
	      if (EXP) vtxtPrintf (pNam, "_%s", ac_table_printable (gCov, ir, 0, "")) ;
	      prodQual  -= 1 ;
	      break ;
	    }      
	}
    }
  else
    prodQual  =  0 ;

  /* quality of the mrna
     find read covering all introns
  */
  gSplicing = ac_tag_table (mrna, "Splicing", h) ;

  /* sanity check */
  mrnaOk = TRUE ;

  for (ir = 0 ; gSplicing && ir < gSplicing->rows ; ir++)
    {  
      a1 = ac_table_int (gSplicing, ir, 0, 1) ;
      if (ir && a1 != a2 + 1)
	mrnaOk = FALSE ;
      a2 = ac_table_int (gSplicing, ir, 1, 1) ;
    }
  if (! mrnaOk)
    goto done ;

  mrnaOk = FALSE ;
  if (EXP && gSplicing &&  gSplicing->rows > 0)
    {
      int xx1, xx2 ; /* limite des introns */
      
      xx1 = ac_table_int (gSplicing, 3, 0, -1) ;
      ir = gSplicing->rows - 1 ;
      xx2 = ac_table_int (gSplicing, 2, ir, -1) ;
    
      gCov =  ac_tag_table (mrna, "Constructed_from", h) ;
      for (ir = 0 ; gCov && ir < gCov->rows ; ir++)
	{
	  if (ac_table_int (gCov, 0, ir, 999999) < xx1 &&
	      (ac_table_int (gCov, 1, ir, -999999) > xx2))
	    {
	      vtxtPrintf (mNam, "_1:%s", ac_table_printable (gCov, ir, 2, "")) ;
	      mrnaOk = TRUE ; 
	      break ;
	    }
	}
    }
  if (EXP && !mrnaOk)
    {
      gCov = ac_tag_table (mrna, "Mrna_covered_by", h);
      if (gCov && gCov->rows > 0)
	{
	  vtxtPrintf (mNam, "_%d:", gCov->rows) ;
	  for (ir = 0 ; gCov && ir < gCov->rows && ir < 4 ; ir++)
	    vtxtPrintf (mNam, "%c%s:", ir == 0 ? ':' : '/', 
			ac_table_printable (gCov, ir, 0, "")) ;

	}
    }

  
  
  if (debug) printf ("p1 = %d   p2 = %d\n", p1, p2);
  
  for (ir = 0 ; gSplicing && ir < gSplicing->rows ; ir++)
    {  
      a1 = ac_table_int (gSplicing, ir, 0, 1) ;
      a2 = ac_table_int (gSplicing, ir, 1, 1) ;
      x1 = ac_table_int (gSplicing, ir, 2, 1) ; 
      x2 = ac_table_int (gSplicing, ir, 3, 1) ;
      cp = ac_table_printable (gSplicing, ir, 4, "") ;
      tag = UNKNOWN ;
      if (strstr (cp, "xon"))
	tag = EXON ;
      else if (strstr (cp, "tron"))
	tag = INTRON ;
      else if (strstr (cp, "ORF_Gap"))
	tag = ORFGAP ;
      else if (strstr (cp, "Gap"))
	tag = GAP ; /* should be GAP */

      cp = ac_table_printable (gSplicing, ir, 5, "") ;
      dictAdd (dict, cp, &type) ;
      /* store all the coordinates in arrays... */
    
      /* cheat with orf gaps for the exports */
      if (tag == EXON)
	for (jr = ir+1 ; jr < gSplicing->rows ; jr++)
	  {  
	    if (strstr (ac_table_printable (gSplicing, jr, 4, ""), "ORF_Gap") ||
		strstr (ac_table_printable (gSplicing, jr, 4, ""), "xon")
		)
	      {
		ir = jr ; 
		a2 = ac_table_int (gSplicing, jr, 1, 1) ;  
		x2 = ac_table_int (gSplicing, jr, 3, 1) ;
	      }
	    else
	      break ;
	  }
   
      if (tag == ORFGAP)
	printf ("# ORFGAP\n") ;

      if (debug) printf ("a1 = %d   a2 = %d  x1 = %d x2 = %d\n", a1, a2, x1, x2) ;
      if (product_id && ! isQual && start && tag == EXON && x1 <= p1 && x2 >= p1)  /* start_codon */
	{
	  c1 = x1 + p1 - x1 ; /* first base of start codon, in mrna coords */
	  c2 = c1 + 2 ;
	  /* no-plato: gff is 1-based, same as aceview */
	  if (y1 <= y2)
	    printf("%s\tAceView\tstart_codon\t%d\t%d\t.\t+\t"
		   , map, y1 + c1 - 1, y1 + c2 - 1) ;
	  else
	    printf("%s\tAceView\tstart_codon\t%d\t%d\t.\t-\t"
		   , map, y1 - c2 + 1, y1 - c1 + 1) ; 

	  printf("0\tgene_id \"%s\"; transcript_id \"%s\"; x1 \"%d\"; x2 \"%d\"; exon_number \"%s\"; "
		   , vtxtPtr (gNam), vtxtPtr (mNam), c1, c2, dictName (dict, type) ) ;
	  if (product_id)
	    printf(" protein_id \"%s\";"
		   , vtxtPtr (pNam));
	  printf("\n") ;
	}
	      
      if (product_id && ! isQual && (tag == EXON || tag == GAP) && x2 >= p1 && x1 <= p2) /* CDS */
	{
	  int p22 = stop ? p2  : p2 ; /* include the stop */
	  c1 = a1 < p1 ? p1 : a1 ; 
	  c2 = a2 < p22 ? a2 : p22 ; 
	  if (c1 <= c2)
	    {
	      ok |= 2 ;
	      if (y1 <= y2)
		printf("%s\tAceView\tCDS\t%d\t%d\t.\t+\t"
		       , map, y1 + c1 - 1, y1 + c2 - 1) ; /* , prodQual */
	      else
		printf("%s\tAceView\tCDS\t%d\t%d\t.\t-\t"
		       , map, y1 - c2 + 1, y1 - c1 + 1) ; /* , prodQual */
	      printf("0\tgene_id \"%s\"; transcript_id \"%s\"; x1 \"%d\"; x2 \"%d\"; exon_number \"%s\";"
		     /* , frame just export a . here */
		     , vtxtPtr (gNam), vtxtPtr (mNam), 	c1, c2, dictName (dict, type));      
	      printf(" protein_id \"%s\";"
		     , vtxtPtr (pNam)) ;

	      printf("\n") ;
	    }
	}
      if (!otherVeryGoodProduct && ! isQual && (tag == ORFGAP || tag == EXON || tag == GAP)) /* exon */
	{
	  c1 = a1 ;
	  c2 = a2 ;
	  if (c1 <= c2)
	    {
	      ok |= 1 ;
	      if (y1 <= y2)
		printf("%s\tAceView\texon\t%d\t%d\t.\t+\t"
		       , map, y1 + c1 - 1, y1 + c2 - 1) ; /*, mrnaQual */
	      else
		printf("%s\tAceView\texon\t%d\t%d\t.\t-\t"
		       , map, y1 - c2 + 1, y1 - c1 + 1) ; /*, mrnaQual */
	      printf(".\tgene_id \"%s\"; transcript_id \"%s\"; x1 \"%d\"; x2 \"%d\"; exon_number \"%s\";"
		     /* , frame just export a . here */
		     , vtxtPtr (gNam), vtxtPtr (mNam)
		     , x1, x2
		     , dictName (dict, type));		      
	      printf("\n") ;
	    }
	}

      if (product_id && ! isQual && stop &&  tag == EXON && x1 <= p2 && x2 >= p2)  /* stop_codon */
	{
	  c1 = x1 + p2 - x1 - 2 ; /* first base of start codon, in mrna coords */
	  c2 = c1 + 2 ;
	  /* pla-to: gff is x1:zero-based, x2:1-based.  base 1 is called 0,1 */
	  if (y1 <= y2)
	    printf("%s\tAceView\tstop_codon\t%d\t%d\t.\t+\t"
		   , map, y1 + c1 - 1, y1 + c2 - 1) ;
	  else
	    printf("%s\tAceView\tstop_codon\t%d\t%d\t.\t-\t"
		   , map, y1 - c2 + 1, y1 - c1 + 1) ;
	  printf("0\tgene_id \"%s\"; transcript_id \"%s\"; x1 \"%d\"; x2 \"%d\"; exon_number \"%s\"; "
		 , vtxtPtr (gNam), vtxtPtr (mNam), c1, c2, dictName (dict, type)) ;
	  if (product_id)
	    printf(" protein_id \"%s\";"
		 , vtxtPtr (pNam));
	  printf("\n") ;
	}

      if (!otherVeryGoodProduct && ! isQual &&  tag == INTRON) /* intron */
	{
	  c1 = a1 ;
	  c2 = a2 ;
	  if (c1 <= c2)
	    {
	      if (y1 <= y2)
		printf("%s\tAceView\tintron\t%d\t%d\t.\t+\t"
		      , map, y1 + c1 - 1, y1 + c2 - 1) ;
	      else
		printf("%s\tAceView\tintron\t%d\t%d\t.\t-\t"
		       , map, y1 - c2 + 1, y1 - c1 + 1) ;
	      printf("0\tgene_id \"%s\"; transcript_id \"%s\"; intron_feet \"%s\"; intron_length \"%s\";"
		     /* , frame just export a . here */
		     , vtxtPtr (gNam), vtxtPtr (mNam)
		     , dictName (dict, type)
		     , ac_table_printable (gSplicing, ir, 7, "") 
		     );		      
	      printf("\n") ;
	    }
	}
    }

 done:
  ac_free (h) ;
  return ok ;
}

/**************************************************************************/
/* look if the mrna is cloud-like */
static BOOL gffGoodMrna (AC_OBJ mrna)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL ok = FALSE ;
  AC_KEYSET ks ;

  if (ac_has_tag (mrna, "gc_ag") || ac_has_tag (mrna, "gt_ag"))
    { ok = TRUE ; goto done ; }
  
  ks = ac_objquery_keyset (mrna, ">cdna_clone", h) ;
  if (ac_keyset_count (ks) > 5)
    { ok = TRUE ; goto done ; }
  /* ks1 = ac_ksquery_keyset (ks, "IS nm_* || IS nr_*", h) ; */

  ks = ac_objquery_keyset (mrna, ">product ; good_product", h) ;
  if (ac_keyset_count (ks) > 0)
    { ok = TRUE ; goto done ; }
  
 done:
  ac_free (h) ;
  return ok ;
}  /* gffGoodMrna */

/**************************************************************************/
/* get the alignment of the exon in spl format */

static Array gffSplGetCosmidDna (AC_OBJ Tg, int *c1p, int *c2p, AC_HANDLE h)
{
  Array cosmidDna = 0 ;
  AC_TABLE gCovers = ac_tag_table (Tg, "Covers", h) ;
  AC_OBJ Cosmid ;
  KEY cosmid = 0 ;
  
  Cosmid = ac_tag_obj (Tg, "Genomic_sequence", h) ;
  cosmid = ac_obj_key (Cosmid) ;
  *c1p = ac_table_int (gCovers, 0, 2, 0) ;
  *c2p = ac_table_int (gCovers, 0, 3, 0) ;
  cosmidDna = dnaGet (cosmid) ;
  
  return cosmidDna ;
} /* gffSplGetCosmidDna */

/**************************************************************************/

static char* gffSplIntronFoot (char *buf, Array cosmidDna, Array cosmidDnaR, int c1, int c2, int a1, int a2)
{
  int u1, u2 ;
  if (c1 < c2) 
    {
      u1 = c1 + a1 - 1 ; u2 = c1 + a2 - 1 ; 
      buf[0] = ace_upper(dnaDecodeChar [(int)arr(cosmidDna, u1 - 1, char)]) ;
      buf[1] = ace_upper(dnaDecodeChar [(int)arr(cosmidDna, u2 - 1, char)]) ;
    }
  else 
    {
      c1 = arrayMax(cosmidDnaR) - c1 + 1 ;
      u1 = c1 + a1 - 1 ; u2 = c1 + a2 - 1 ; 
      buf[0] = ace_upper(dnaDecodeChar [(int)arr(cosmidDnaR, u1 - 1, char)]) ;
      buf[1] = ace_upper(dnaDecodeChar [(int)arr(cosmidDnaR, u2 - 1, char)]) ;
    }
  buf[2] = 0 ;

  return buf ;
}  /* gffSplIntronFoot */

/**************************************************************************/
/* per cent of A in given exon */
static int gffSplAcontent (Array cosmidDna, Array cosmidDnaR, int c1, int c2, int a1, int a2)
{
  Array aa ;
  const char *ccp ;
  int nA, nn ;
  int u, u1, u2 ;
  if (c1 < c2) 
    {
      aa = cosmidDna ;
    }
  else 
    {
      c1 = arrayMax(cosmidDnaR) - c1 + 1 ;
      aa = cosmidDnaR ;
    }
  u1 = c1 + a1 - 1 ; u2 = c1 + a2 - 1 ;
  for (nA = 0, nn = 0, u = u1, ccp = arrp (aa, u-1,char) ; u <= u2 ; ccp++, nn++, u++)
    if (*ccp == A_) nA++ ;
  if (nn == 0) nn = 1 ;
  nA *= 100 ;
  return nA/nn ;
} /* gffSplAcontent */

/**************************************************************************/
/* get the alignment of the exon in spl format */

static int gffSplGetAli (Array rDna, int x1, int *x2p
			 , Array cosmidDna, Array cosmidDnaR
			 , int a1, int *a2p, int c1, int c2
			 , Array err
			 , int *nAlip, int *nExactp, double *eQualp, vTXT eBuf
			 , BOOL isNm, BOOL futureExon
			 , AC_HANDLE h
			 )
{
  int ii, nErr, nHits = 0, nAli, nAliMin ;
  int a2 = *a2p, x2 = *x2p ;
  int errMax ;
  int u1, u2, oldx, dx ;
  BOOL isDown ;

  /* accepts biologists coordinates (1-based, extremities included) 
   * always give *a1p < *a2p
   * Array aceDnaDoubleTrackErrors (Array  dna1, int *x1p, int *x2p, BOOL isDown,
                                      Array dna2, Array dna2R, int *a1p, int *a2p, 
                                      int *NNp, Array err, int maxJump, int maxError, BOOL doExtend)
  */
    
  if (c1 < c2) { u1 = c1 + a1 - 1 ; u2 = c1 + a2 - 1 ; }
  else { u1 = c1 - a1 + 1 ; u2 = c1 - a2 + 1 ; }
  if (u1 > u2) /* swap */
    {
      Array aa ;
   
      aa = cosmidDna ; cosmidDna = cosmidDnaR ; cosmidDnaR = aa ;
      u1 = arrayMax(aa) - u1 + 1 ;
      u2 = arrayMax(aa) - u2 + 1 ;
    }
  isDown = (x1 < x2 ? TRUE : FALSE) ;
  if (isDown && isNm && !futureExon) /* clip the terminal A */
    {
      char *cp ;
      for (cp = arrp (rDna, x2 - 1, char) ; *cp == A_ ; cp--)
	{ x2-- ; a2-- ; u2-- ; x2p-- ; a2p-- ;}
    }
  nErr = 0 ; nAli = 0 ;
  nAliMin = a1 < a2 ? a2 - a1 + 1 : a1 - a2 + 1 ;
  dx = isDown ? x2 - x1 + 1 : x1 - x2 + 1 ;
  if (nAliMin < dx) nAliMin = dx ;
  if (isDown)
    oldx = x1 - 2 ; /* just before - Plato */
  else
    oldx = x1 + 1 - 1 ; /* just after - Plato */
  aceDnaDoubleTrackErrors (rDna, &x1, &x2, isDown /* bio coordinates  */
			   , cosmidDna, cosmidDnaR, &u1, &u2
			   , 0, err, 8, 0, FALSE, 0) ;

  errMax = err ? arrayMax(err) : 0 ;
  if (errMax)
    {
      A_ERR *ep ;
      int nR = 0, nI = 0, nD = 0, nExact = 0 ;

      for (ii = 0, ep = arrp (err, ii, A_ERR) ; ii < errMax ; ii++, ep++)
	{
	  switch (ep->type)
	    {
	    case ERREUR:
	      if (nI) { vtxtPrint (eBuf, "I") ; if (nI > 1)  vtxtPrintf (eBuf, "%d", nI) ; nI = 0 ; }
	      if (nD) { vtxtPrint (eBuf, "D") ; if (nD > 1)  vtxtPrintf (eBuf, "%d", nD) ; nD = 0 ; }
	      dx = isDown ? ep->iShort - oldx - 1 : oldx - ep->iShort - 1 ;
	      if (dx > 0)
		{
		  if (nR) { vtxtPrint (eBuf, "R") ; if (nR > 1)  vtxtPrintf (eBuf, "%d", nR) ; nR = 0 ; }
		  vtxtPrintf (eBuf, "M") ; if(dx > 1) vtxtPrintf (eBuf, "%d", dx) ;
		  nExact += dx ;
		}
	      oldx = ep->iShort ;
	      nR++ ;
	      nErr++ ;
	      break ;
	    case TROU:
	    case TROU_DOUBLE:
	    case TROU_TRIPLE:
	      if (nR) { vtxtPrint (eBuf, "R") ; if (nR > 1)  vtxtPrintf (eBuf, "%d", nR) ; nR = 0 ; }
	      if (nD) { vtxtPrint (eBuf, "D") ; if (nD > 1)  vtxtPrintf (eBuf, "%d", nD) ; nD = 0 ; }
	      dx = isDown ? ep->iShort - oldx - 1 : oldx - ep->iShort - 1 ;
	      if (dx > 0)
		{
		  if (nI) { vtxtPrint (eBuf, "I") ; if (nI > 1)  vtxtPrintf (eBuf, "%d", nI) ; nI = 0 ; }
		  vtxtPrintf (eBuf, "M") ; if(dx > 1) vtxtPrintf (eBuf, "%d", dx) ;
		  nExact += dx ;
		}
	      oldx = isDown ? ep->iShort - 1 : ep->iShort + 1 ;
	      switch (ep->type)
		{
		case TROU:
		  nErr++ ;
		  nI++ ;
		  break ;
		case TROU_DOUBLE:
		  nErr += 2 ;
		  nI += 2 ;
		  break ;
		case TROU_TRIPLE:
		  nErr += 2 ;
		  nI += 2 ;
		  break ;
		default:
		  break ;
		}
	      break ;
	    case INSERTION: 
	    case INSERTION_DOUBLE:
	    case INSERTION_TRIPLE:
	      if (nI) { vtxtPrint (eBuf, "I") ; if (nI > 1)  vtxtPrintf (eBuf, "%d", nI) ; nI = 0 ; }
	      if (nR) { vtxtPrint (eBuf, "R") ; if (nR > 1)  vtxtPrintf (eBuf, "%d", nR) ; nR = 0 ; }
	      dx = isDown ? ep->iShort - oldx - 1 : oldx - ep->iShort - 1 ;
	      if (dx > 0)
		{
		  if (nD) { vtxtPrint (eBuf, "D") ; if (nD > 1)  vtxtPrintf (eBuf, "%d", nD) ; nD = 0 ; }
		  vtxtPrintf (eBuf, "M") ; if(dx > 1) vtxtPrintf (eBuf, "%d", dx) ;
		  nExact += dx ;
		}
	      oldx = isDown ? ep->iShort + 1 : ep->iShort - 1 ;
	      switch (ep->type)
		{
		case INSERTION: 
		  nErr++ ;
		  nAli++ ;
		  nD++ ;
		  break ;
		case INSERTION_DOUBLE:
		  nErr += 2 ;
		  nAli += 2 ;
		  nD += 2 ;
		  break ;
		case INSERTION_TRIPLE:
		  nErr += 3 ;
		  nAli += 3 ;
		  nD += 3 ;
		  break ;
		default:
		  break ;
		}
	      break ;
	    default:
	      break ;
	    }
	}
      if (nR) { vtxtPrint (eBuf, "R") ; if (nR > 1)  vtxtPrintf (eBuf, "%d", nR) ; nR = 0 ; }
      if (nD) { vtxtPrint (eBuf, "D") ; if (nD > 1)  vtxtPrintf (eBuf, "%d", nD) ; nD = 0 ; }
      if (nI) { vtxtPrint (eBuf, "I") ; if (nI > 1)  vtxtPrintf (eBuf, "%d", nI) ; nI = 0 ; }

      dx = isDown ? x2 - oldx - 1 : oldx - x2 - 1 ;
      if (dx > 0)
	{
	  vtxtPrintf (eBuf, "M") ; if(dx > 1) vtxtPrintf (eBuf, "%d", dx) ;
	  nExact += dx ;
	}
      nAli += nExact ;
      *nExactp = nExact ;
    }
  else
    { vtxtPrintf (eBuf, "M%d", nAliMin) ; *nExactp = nAliMin ; }
  if (nAli < nAliMin) nAli = nAliMin ;
  *nAlip = nAli ; 
  *eQualp = (nAli - (double)nErr)/nAli ;
  
  return nHits ;
} /* gffSplGetAli */

/**************************************************************************/

static int gffSplDumpTg (AC_DB  db, AC_OBJ Tg, BOOL debug)
{
  AC_HANDLE h = ac_new_handle () ;
  Array err = 0, rDna = 0, cosmidDna = 0, cosmidDnaR = 0 ;
  int qual, ok = 0, ok1, tg1, tg2, a1, a2, b1, b2, c1, c2, x1, x2, ir, jr, jr2, nerr, nExact, nAli ;
  int nExactTotal, nTotal, chromNumber ;
  char *rNam, *chrNam, *chromNam, footBuf[3] ; ;
  const char *ccp ;
  BOOL isDown, isReadDown ;
  AC_OBJ Rr, Chrom ;
  AC_TABLE gDb, gMap, gAF, gReads ;
  KEY rr ;
  static int nCluster = 0 ;
  double eQual ;
  vTXT eBuf = vtxtHandleCreate (h) ;
  BOOL oldExon, isNm, futureExon ;

  gMap = ac_tag_table (Tg, "IntMap", h) ;
  if (!gMap) goto done ;
  chromNam = hprintf (h, "%s", ac_table_printable (gMap, 0, 0, "toto")) ;
  Chrom = ac_get_obj (db, "Sequence", chromNam, h) ;
  if (!Chrom) goto done ;
  gDb = ac_tag_table (Chrom, "Database", h) ;
  ccp = ac_table_printable (gDb, 0, 1, 0) ;
  if (ccp)
    chrNam = hprintf (h, "%s", ccp) ;
  else
    chrNam = hprintf (h, "%s.1", ac_name(Chrom)) ;
  if (strstr(chromNam, "Y")) chromNumber = 24 ;
  else if (strstr(chromNam, "X")) chromNumber = 23 ;
  else if (strstr(chromNam, "M")) chromNumber = 25 ;
  else if (sscanf (chromNam, "%d", &ir) == 1)
    chromNumber = ir ;
  else
    chromNumber = 99 ;
  tg1 = ac_table_int (gMap, 0, 1, 1) ;
  tg2 = ac_table_int (gMap, 0, 2, 1) ;
  isDown = tg1 < tg2 ? TRUE : FALSE ;
  
  cosmidDna = gffSplGetCosmidDna (Tg, &c1, &c2, h) ;
  cosmidDnaR = dnaCopy (cosmidDna) ;
  reverseComplement (cosmidDnaR) ;

  gReads = ac_tag_table (Tg, "Read", h);
  gAF = ac_tag_table (Tg, "Assembled_from", h);
  if (!gAF || ! gReads) goto done ;
  for (ir = 0 ; ir < gReads->rows ; ir++)
    {
      ok1 = 0 ;
      ac_free (rDna) ; oldExon = FALSE ;
      qual = ac_table_int (gReads, ir, 5, 0) ;
      if (qual >= 17)
	continue ; /* drop doubtful alignments */
      rr = ac_table_key (gReads, ir, 0, 0) ;
      if (!rr) continue ;
      Rr = ac_table_obj (gReads,ir,0, h) ;
      isNm = ac_has_tag (Rr, "Ref_seq") ;
      gDb = ac_tag_table (Rr, "Database", h) ;
      ccp = ac_table_printable (gDb, 0, 1, 0) ;
      if (ccp)
	rNam = hprintf (h, "%s", ccp) ;
      else
	rNam = hprintf (h, "%s.1", ac_key_name(rr)) ;
      nExactTotal = 0 ; nTotal = 0 ;
      for (jr = 0 ; jr < gAF->rows ; jr++)
	{
	  if (rr != ac_table_key (gAF, jr, 2, 0)) 
	    continue ;
	  if (0 && strcmp(ac_key_name(rr), "NM_007100"))
	    continue ;
	  if (!ok1) { ok1 = 1 ; nCluster++ ;}
	  a1 = ac_table_int (gAF, jr, 0, 1) ;
	  a2 = ac_table_int (gAF, jr, 1, 1) ;
	  x1 = ac_table_int (gAF, jr, 3, 1) ;
	  x2 = ac_table_int (gAF, jr, 4, 1) ;
	  isReadDown = x1 < x2 ? TRUE : FALSE ;
	  nerr =  ac_table_int (gAF, jr, 5, 0) ;
	  vtxtClear (eBuf) ;
	  futureExon = FALSE ;
	   for (jr2 = jr + 1 ; ! futureExon && jr2 < gAF->rows ; jr2++)
	     {
	       if (rr == ac_table_key (gAF, jr2, 2, 0)) 
		 futureExon = TRUE ;
	     }
	  if (1 || nerr)
	    {
	      if (!rDna)
		rDna = dnaGet (rr) ;
	      if (!err)
		err = arrayHandleCreate (32, A_ERR, h) ;
	      gffSplGetAli (rDna, x1, &x2, cosmidDna, cosmidDnaR, a1, &a2, c1, c2, err, &nAli, &nExact, &eQual, eBuf, isNm, futureExon, h) ;
	      if (0 && arrayMax(err) && !nerr)
		printf ("\t****nerr==0/%d****", arrayMax(err)) ;
	    }
	  else
	    {
	      nAli = a2 - a1 + 1 ; 
	      nExact = a2 - a1 + 1 ; 
	      eQual = 1 ; 
	      vtxtPrintf (eBuf, "M%d", nAli) ;
	    }
	  printf("%s%02d%08d", isReadDown ? "+" : "-", chromNumber, nCluster) ;
	  /* nali should be the length of the ali + the - for indels */
	  printf ("\t%s\t%s\t%g\t%d", rNam, chrNam, eQual, nAli) ;
	  if (isDown) { b1 = tg1 + a1 - 1 ; b2 = tg1 + a2 - 1 ; }
	  else { b1 = tg1 - a1 + 1 ; b2 = tg1 - a2 + 1 ; }
	  printf ("\t%d\t%d\t%d\t%d\t", x1, x2, b1, b2) ;
	  nExactTotal += nExact ;
	  nTotal += nAli ;
	  if (oldExon) 
	    printf ("%s", gffSplIntronFoot (footBuf, cosmidDna, cosmidDnaR, c1, c2, a1-2, a1-1)) ;
	  printf ("<exon>") ; 
	  oldExon = TRUE ;
	  if (futureExon) 
	    printf ("%s", gffSplIntronFoot (footBuf, cosmidDna, cosmidDnaR, c1, c2, a2+1, a2+2)) ;
	  printf ("\t%s", vtxtPtr (eBuf) ? vtxtPtr (eBuf) : " ") ;
	  if(0)nExactTotal *= 100 ;
	  if (nTotal == 0) nTotal = 1 ;
	  if(0)nExactTotal /= nTotal ;
	  if (0 && !futureExon)
	    printf ("\t%d\t%d\t%d", gffSplAcontent (cosmidDna, cosmidDnaR, c1, c2, a1, a2), nExactTotal, nTotal) ;
	  printf ("\n") ;
	}
      ok += ok1 ;
    }
 done:
  ac_free (rDna) ;
  ac_free (cosmidDna) ;
  ac_free (cosmidDnaR) ;
  ac_free (h) ;
  return ok ;
} /* gffSplDumpTg */

/**************************************************************************/
/**************************************************************************/
/* minimal gtf dump of the predicted genes to be used by magic metadata.tcsh */
static int gffDumpPg (AC_DB  db, AC_OBJ Pg, BOOL debug)
{
  AC_HANDLE h = ac_new_handle () ;
  int ok = 0, ir, jr ;
  int a1, a2, x1, x2, u1, u2, pass, nIntrons ;
  BOOL isDown ;
  AC_TABLE iMap, exons ;
  const char *ccp ;

  iMap = ac_tag_table (Pg, "IntMap", h) ;
  if (! (iMap && iMap->rows > 0 && iMap->cols >= 3))
    goto done ;
  exons = ac_tag_table (Pg, "Source_Exons", h) ;
  if (! (exons && exons->rows > 0 && exons->cols >= 2))
    goto done ;

  a1 = ac_table_int (iMap, 0, 1, 0) ;
  a2 = ac_table_int (iMap, 0, 2, 0) ;
  isDown = a1 < a2 ? TRUE : FALSE ;

  /* count the introns */
  nIntrons = 0 ;
  for (ir = 1 ; ir < exons->rows ; ir++)
    if (ac_table_int (exons, ir, 0, 0) > ac_table_int (exons, ir - 1, 1, 0) + 1)
      nIntrons++ ;

  for (pass = 1 ; pass >= 0 ; pass--)
    {
      for (ir = jr = 0 ; ir < exons->rows ; ir++, jr++)
	{
	  x1 = ac_table_int (exons, ir, 0, 0) ;
	  x2 = ac_table_int (exons, ir, 1, 0) ;
	  ccp = ac_table_printable (exons, ir, 2, 0) ;
	  if (x1 < 0 || x2 < 0 || x1 > x2)
	    {
	      messerror ("x1 = %d < 0 or x2 = %d  < 0 or x1 < x2  in %s \n", x1, x2, ac_name (Pg)) ;
	      continue ;
	    } 
	  if (pass && (! ccp || (strcasecmp (ccp, "Exon") && strcasecmp (ccp, "CDS"))))
	    continue ;

	  if ( pass == 0) /* fuse CDS/UTR belonging so same exon */
	    {
	      while (ir < exons->rows - 1 &&
		  x2 + 1 == ac_table_int (exons, ir + 1, 0, 0)
		  )
		{ x2 = ac_table_int (exons, ir + 1, 1, 0) ; ir++ ; }
	    }
	  if (isDown)
	    { u1 = a1 + x1 - 1 ; u2 = a1 + x2 - 1 ; } /* acedb conventions */
	  else
	    { u1 = a1 - x1 + 1 ; u2 = a1 - x2 + 1 ; } /* acedb conventions */
	  
	  if (isDown) ;   /* gtf conventions */
	  else
	    { int u0 = u1 ; u1 = u2 ; u2 = u0 ; }   /* gtf conventions */
	  
	  printf ("%s\t%s%s\t%s\t%d\t%d\t.\t%c\t.\t"
		  , ac_table_printable (iMap, 0, 0,"unmapped")
		  , nIntrons > 0 ? "spliced"  : "single_exon"
		  , pass == 0 ?  "_transcript" : "_protein_coding"
		  , pass == 0 ? "exon" : "CDS"
		  , u1, u2
		  , isDown ? '+' : '-'
		  ) ;
	  
	  ok = 0 ;
	  ccp = ac_tag_printable (Pg, "Model_of", 0) ;
	  if (ccp)
	    printf (" gene_id \"%s\";", ccp) ; 
	  if (ccp && ! strncmp (ccp, "X__",3))
	    printf (" gene_name \"%s\";", ccp+3) ; 
	    
	  ccp = ac_name (Pg) ;
	  printf (" transcript_id \"%s\";", ccp) ; 
	  if (ccp && ! strncmp (ccp, "_",1))
	    printf (" transcript_name \"%s\";", ccp+1) ; 
	  ccp = ac_tag_printable (Pg, "NM_id", 0) ;
	  if (ccp)
	    printf (" NM_id \"%s\";", ccp) ; 
	  ccp = ac_tag_printable (Pg, "GeneId_pg", 0) ;
	  if (ccp)
	    printf (" NCBI_GeneId \"%s\";", ccp) ; 
	  
	  if (pass == 0)
	    printf (" exon_number \"%d\";", jr + 1) ;
	  
	  printf ("\n") ;
	}
    }
  ok = 1 ;

 done:
  ac_free (h) ;
  return ok ;
}
/**************************************************************************/
/**************************************************************************/

static void usage (void)
{
  char * usage = "gffdump [acedb_dir | host:port] [-no_cloud] [-multiProduct] [-qual] [-filter filter (i.e. \'IS *3\')]\n" ;
  
  fprintf (stderr, "Usage: %s\n", usage);
  fprintf (stderr, "   acedb_dir || host : the acedb database\n") ; 
  fprintf (stderr, "   -no_cloud : reject cloud_genes \n") ; 
  fprintf (stderr, "   -multiProduct : dump very_good non best products\n") ; 
  fprintf (stderr, "   -qual  : exports Bad_quality for all problematic entries\n") ; 
  fprintf (stderr, "   -filter <query> : i.e. IS a* (any acedb query is ok) \n") ;  
  fprintf (stderr, "   -spl : dump in NCBI_splign (kapustin) format \n") ;  
  fprintf (stderr, "   -pg [-nmid] : dump the predicted_genes [with nmid] \n") ;  
  
  exit (1) ;
} /* usage */

/**************************************************************************/

int main (int argc, const char * argv[]) 
{
  AC_HANDLE h ;
  AC_DB  db;
  AC_ITER  tgs ;
  AC_TABLE mrnas ;
  AC_OBJ tg, mrna;
  const char *dbName = "(missing dabase argument)" ;
  const char *err = 0 ;
  int debug = 0, ir, np = 0, nm = 0, ng = 0, nm1, np1, pass, ok ;
  const char *filter = "" ; /* "IS * && COUNT {>product ;peptide:2 >= 100} > 0"; */
  char *qq ;
  BOOL isPg, hasNm, multiProduct, noCloud = FALSE, isQual = FALSE, isSpl = FALSE ;
  
  
  if (argc <2)
    usage () ;
  
  dbName = argc>=2 ? argv[1] : "." ;
  getCmdLineOption (&argc, argv, "-filter", &filter) ;
  noCloud =  getCmdLineOption (&argc, argv, "-no_cloud", 0) ;  
  isSpl =  getCmdLineOption (&argc, argv, "-spl", 0) ;  
  multiProduct =  getCmdLineOption (&argc, argv, "-multiProduct", 0) ; 
  isPg = getCmdLineOption (&argc, argv, "-pg", 0) ;  
  hasNm = getCmdLineOption (&argc, argv, "-nmid", 0) ;  
  
  isQual = getCmdLineOption (&argc, argv, "-qual", 0) ;
  if (getCmdLineOption (&argc, argv, "-h", 0)) usage () ;
  if (getCmdLineOption (&argc, argv, "-help", 0)) usage () ;
  if (getCmdLineOption (&argc, argv, "--help", 0)) usage () ;
 
  if (argc <2)
    usage () ;
  
  if (!( db = ac_open_db (dbName, &err)))
    {
      printf("ERROR: Cannot open the database %s, %s sorry\n\n", dbName, err ? err : "");
      exit (1) ;
    }

  if (isPg && hasNm)
    qq = "find Predicted_gene model_of && NM_id " ;
  else if (isPg)
        qq = "find Predicted_gene model_of " ;
  else if (noCloud)
    qq = "find gene ! cloud_gene ; >transcribed_gene";
  else
    qq = "find transcribed_gene";
  fprintf (stderr, "// %s\n", qq) ;
  if (filter)
    qq = messprintf ("%s ; %s", qq, filter) ;

  tgs = ac_dbquery_iter (db, qq, 0) ;
  tg = 0 ;
  while (ac_free (tg), tg = ac_iter_obj (tgs))
    {
      ng++ ; nm1 = 0 ;
      if (isPg)
	{
	  nm += gffDumpPg (db, tg, debug) ;
	  continue ; 
	}
      if (isSpl)
	{
	  nm += gffSplDumpTg (db, tg, debug) ;
	  continue ; 
	}
      h = ac_new_handle () ;
      mrnas = ac_tag_table (tg, "mrna", h) ; 
      for (pass = 0 ; pass < 2 ; pass++)
	for (ir = 0 ; mrnas && ir < mrnas->rows ; ir++)
	  {
	    mrna = ac_table_obj (mrnas, ir, 0, h) ;
	    if (! pass && ! ac_has_tag (mrna, "best_in_gene"))
	      continue ;
	    if (pass && ac_has_tag (mrna, "best_in_gene"))
	      continue ;
	    if (1 &&   /* 0 => dump all trancript even bad ones == SEQC */ 
		nm1 && noCloud && !gffGoodMrna (mrna))
	      continue ; 
	    else
	      {
		np1 = 0 ;
		ok = 2 ;
		while (ok & 0x2) /* we may need to export yet another product */
		  {
		    ok = gffDumpMrna (mrna, isQual, np1, debug) ;
		    if (ok & 0x1) nm1 = 1 ; /* exons exported */
		    if (ok & 0x2) np1++ ; /* CDS exported */
		    if (! multiProduct)
		      break ;
		  }
		if (nm1) nm++ ;
		np += np1 ;
	      }
	  }
      ac_free (h) ;
    }
  
  fprintf (stderr, "// Exported %d products, %d mrnas in %d genes)\n", np, nm, ng) ;

  ac_db_close (db) ; 
  return 0 ;
}

/**************************************************************************/
/**************************************************************************/
