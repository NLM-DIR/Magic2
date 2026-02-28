#include "../wac/ac.h"
#include "../wfiche/gtitle.h"
#include "pick.h"

static BOOL debug = FALSE ;
static const char *gtGene2Locus (AC_OBJ oGene) ;
static char *gtProductBaseTitle (vTXT blkp, GMP *gmp, AC_OBJ oProduct, AC_OBJ oGF, BOOL isGene, BOOL *isPutativep) ;
static char *gtGeneDescriptorTag (vTXT blkp, GMP *gmp) ;
static char* gtJavaScriptProtect (const char *ptr, AC_HANDLE h) ;

/********************************************************************/
/********************************************************************/
/************************* JUMPS ************************************/
/* section is the offset in the table, it allows to identify (say to open/close) a section by number
 * the number of subsection should not go over 1000
 */
typedef struct { int used ; int level; BOOL isOpen; int index; char nam [64], txt[1024], file[64], help[64] ; }  JUMPPOINTMENU ;
static char  jumpPointGeneNam[1024], jumpPointMrnaNam[1024] ;
static JUMPPOINTMENU  jumpPointMenu [] =
{
  /* a level 3 marks nothing but returns its help pointer */
  /* a level 2 marks level 1 above and level 9 under */
  /* a level in range 10-100 is javascript */
  /* a level in range > 100 is a help on click submenu */

  /* a hash in call to gmpSection makes the title RED: i.e.  gmpSection (blkp, gmp, "#tg_disease",..) */
  /***********************/
  /****** gene context ***/
  /***********************/
  /*** gene GENE ****/
  { 0 , -1 , TRUE , 0 , "SUMMARY" , "summary" , "", "" } ,
  { 0 , 3 , TRUE , 0 , "tg_RefSeq_Summary" , "", "HelpGene", "" } ,
  { 0 , 3 , TRUE , 0 , "tg_Manual_Summary" , "", "HelpGene" , "" } ,
  { 0 , 3 , TRUE , 0 , "tg_AceView_Summary" , "", "HelpGene", "" } ,
  { 0 , 2 , TRUE , 0 , "Gene_genome_diagram" , "" , "", "" } ,
  { 0 , 2 , TRUE , 0 , "RNA_seq" , "" , "", "" } ,
  { 0 , 2 , TRUE , 0 , "Gene_expression_profile_in_primates" , "gxpProfile", "", "" } ,
  { 0 , 2 , FALSE , 0 , "Gene_expression_profile_in_primates_caption", "", "", "" } ,
  { 0 , 2 , TRUE , 0 , "Primate_wiggles" , "", "", "" } ,
  { 0 , 2 , FALSE , 0 , "Primate_wiggles_caption", "", "", "" } ,
  { 0 , 2 , FALSE , 0 , "Gene_mrna_diagrams" , "" , "", "" } ,
  { 0 , 2 , TRUE , 0 , "Gene_compact_diagram_1" , "" , "", "" } ,
  { 0 , 2 , FALSE , 0 , "Gene_compact_diagram_2" , "" , "", "" } ,
  { 0 , 2 , FALSE , 0 , "Gene_compact_diagram_3" , "" , "", "" } ,
  { 0 , 2 , FALSE, 0 , "Gene_compact_diagram_4" , "" , "", "" } ,
  { 0 , 2 , FALSE , 0 , "Gene_compact_diagram_5" , "" , "", "" } ,
  { 0 , 2 , FALSE , 0 , "Gene_compact_diagram_6" , "" , "", "" } ,
  { 0 , 2 , FALSE , 0 , "Gene_compact_diagram_caption", "", "", "" } ,
  { 0 , 2 , TRUE , 0 , "Gene_sequences" , "", "HelpGene", "" } ,
  { 0 , 2 , TRUE , 0 , "Gene_wiggle_diagram_1" , "" , "", "" } ,
  { 0 , 2 , FALSE , 0 , "Gene_wiggle_diagram_caption", "", "", "" } ,
  { 0 , 2 , TRUE , 0 , "Gene_locator" ,  "Locator", "HelpGene", "" } ,
  { 0 , 2 , FALSE , 0 , "Gene_locator_caption", "", "", "" } ,
  { 0 , 2 , FALSE , 0 , "tg_Alias_map" , "Map<br>and<br>Links" , "HelpGene", "HelpAlias" } ,
  { 0 , 2 , TRUE , 0 , "Biblio" ,  "Bibliography" , "HelpGene", "" } ,
  /*** gene MOLECULES ***/
  { 0 , 1 , FALSE , 0 , "Gene_MOLECULAR" , "Molecules" , "HelpGene", "" } ,
  { 0 , 2 , TRUE , 0 , "tg_alt_features" , "" , "HelpGene", "HelpAltFeature" } ,
  { 0 , 2 , FALSE , 0 , "mRNAs" , "" , "HelpGene", "HelpTranscripts" } ,
  { 0 , 2 , FALSE , 0 , "mrnaAnnotTable" , "" , "HelpGene", "HelpmrnaAnnotTable" } ,
  { 0 , 2 , FALSE , 0 , "tg_introns" , "" , "HelpGene", "HelpIntronsExons" } ,
  { 0 , 2 , FALSE , 0 , "tg_promotors" , "" , "HelpGene", "HelpIntronsExons" } ,
  { 0 , 2 , FALSE , 0 , "tg_polyA" , "" , "HelpGene", "HelpIntronsExons" } ,
  { 0 , 2 , FALSE , 0 , "mRNA_sequences" , "" , "HelpGene", "HelpTranscripts" } ,
  { 0 , 2 , FALSE , 0 , "Shed_variants" , "" , "HelpGene", "HelpShed" } ,
  { 0 , 2 , FALSE , 0 , "Proteins" , "" , "HelpGene", "HelpProteins" } ,
  { 0 , 2 , FALSE , 0 , "Introns_and_exons" , "" , "" } ,
  { 0 , 2 , TRUE , 0 ,  "tg_dna" , "" , "" } ,
  { 0 , 2 , TRUE , 0 ,  "tg_pep" , "" , "" } ,
  
  /*** gene mrna list ***/ 
  { 0 , 2 , TRUE , 0 , "gene_mrna_list" , "mRNAs" , "", "" } ,

  /*** gene EXPRESSION  ****/
  { 0 , 1 , TRUE , 0 , "Gene_EXPRESSION" , "Expression Tissue" , "HelpGene", "" } ,
  { 0 , 2 , TRUE , 0 , "tg_expression" , "" , "HelpGene", "HelpExpression" } ,
  { 0 , 2 , TRUE , 0 , "tg_cdna_support" , "" , "HelpGene", "HelpExpression" } ,
  { 0 , 2 , TRUE , 0 , "tg_expression_profile" , "" , "HelpGene", "HelpExpressionProfile" } ,

  /*** gene FUNCTION and REGULATION ****/
  { 0 , 1 , TRUE , 0 , "Gene_FUNCTION" , "Function" , "HelpGene", "" } ,
  { 0 , 2 , TRUE , 0 , "tg_disease" , "" , "HelpGene", "HelpGO" } ,
  { 0 , 2 , TRUE , 0 , "tg_pathway" , "" , "HelpGene", "HelpGO" } ,
  { 0 , 2 , TRUE , 0 , "tg_Worm_Phenotype" , "" , "HelpGene", "HelpPhenotype" } ,
  { 0 , 2 , TRUE , 0 , "tg_Phenotype" , "" , "HelpGene", "HelpPhenotype" } ,
  { 0 , 2 , TRUE , 0 , "tg_regulation" , "" , "HelpGene", "HelpRegulation" } ,
  { 0 , 2 , TRUE , 0 , "tg_product", "" , "HelpGene", "HelpProducts" } ,
  { 0 , 2 , TRUE , 0 , "tg_product_psort", "" , "HelpGene", "HelpProducts" } ,
  { 0 , 2 , TRUE , 0 , "tg_product_pfam", "" , "HelpGene", "HelpProducts" } ,
  { 0 , 2 , TRUE , 0 , "tg_acekog" , "" , "HelpGene", "HelpInteractions" } ,
  { 0 , 2 , TRUE , 0 , "tg_taxblast" , "" , "HelpGene", "HelpInteractions" } ,
  { 0 , 2 , TRUE , 0 , "tg_protein_interaction" , "" , "HelpGene", "HelpInteractions" } ,
  { 0 , 2 , TRUE , 0 , "tg_gene_interaction" , "" , "HelpGene", "HelpInteractions" } ,
  { 0 , 2 , TRUE , 0 , "tg_interaction_caption", "", "", "" } ,
  { 0 , 2 , TRUE , 0 , "tg_best_friends", "", "", "" } ,
  { 0 , 2 , TRUE , 0 , "BiblioRif",  "RIF" , "HelpGene", "" } ,
  /*** gene PROTEINS ***/
  { 0 , 1 , TRUE , 0 , "GENE_PROTEINS" ,  "Protein" , "HelpGene", "" } ,
  
  /*** gene clones ***/ 
  { 0 , 1 , FALSE , 0 , "cDNA_clones" ,  "cDNA clones" , "HelpGene", "" } ,
  { 0 , 81 , FALSE , 0 , "allSupportingClones" , "" , "HelpGene", "HelpAllClones" } ,
  { 0 , 2 , FALSE , 0 ,  "tgSupportingClones" , "" , "HelpGene", "HelpMainClones" } ,
  { 0 , 2 , FALSE , 0,  "NMSupportingClones" , "" , "HelpGene", "helpMainClones" } ,

  /***********************/
  /*** clone page *********/
  /***********************/
  /*** clone summary ***/
  { 0 , 1 , TRUE , 0 , "_cloneSummary" ,  "Gene  page" , "HelpGene", "" } ,
  { 0 , 1 , FALSE , 0 , "clone_AceView_Summary" , "<b>AceView</b> Summary" , "HelpGene", "" } ,
  { 0 , 1 , TRUE , 0 , "clone_table" , "Clone table" , "HelpGene", "" } ,

  /***********************/
  /*** mRNA page *********/
  /***********************/
  /*** mrna summary ***/
  { 0 , -1 , TRUE , 0 , "mRNA_summary" ,  "" , "HelpmRNA", "" } ,
  /*** mrna diagram ***/
  { 0 , 2 , TRUE , 0 , "mRNA_diagram" , "Diagram" , "", "" } ,
  { 0 , 2 , FALSE , 0 , "MrnaDiagramCaption", "", "", "" } ,
  /*** mrna expression and clone support ***/
  { 0 , 2 , FALSE , 0 , "mRNA_expression" ,  "mRNA Expression" , "HelpmRNA", "" } ,
  /*** protein annotation ***/
  { 0 , 2 , FALSE, 0 , "Protein" ,  "Protein" , "HelpmRNA", "helpmRNASummary" } ,
  { 0 , 2 , FALSE , 0 , "MW" ,"" , "HelpmRNA", "" } ,
  { 0 , 2 , FALSE , 0 , "Primers" , "" ,  "HelpmRNA", "HelpPrimer" } ,
  { 0 , 2 , FALSE , 0 , "Psort" , "" , "HelpmRNA", "helpPredictedCellularLocalization" } ,
  { 0 , 2 , FALSE , 0 , "mrnaPfam" , "" , "HelpmRNA", "helpProteinFamily" } ,
  { 0 , 2 , FALSE , 0 , "BlastP" , "" , "HelpmRNA", "helpProteinHomologies" } ,
  { 0 , 2 , FALSE , 0 , "Taxonomy" , "" , "HelpmRNA", "helpLineageAndClosestHomologs" } ,

  /*** mRNA structure and sequence ***/ 
  { 0 , 2 , FALSE , 0 , "mRNA_structure" ,  "Structure" , "HelpmRNA", "helpmRNAstructure" } ,
  { 0 , 2 , FALSE , 0 , "mRNA_struct" , "" , "HelpmRNA", "helpmRNAstructure" } ,
  { 0 , 2 , FALSE , 0 , "CDS_struct" , "" , "HelpmRNA", "helpmRNAstructure" } ,
  { 0 , 2 , FALSE , 0 , "mRNA_Protein_sequence", "", "HelpmRNA", "" } ,
  { 0 , 2 , FALSE , 0 , "mRNA_sequence" , "", "HelpmRNA", "" } ,
  { 0 , 2 , FALSE , 0 , "mRNA_AM_sequence" , "", "HelpmRNA", "" } ,
  { 0 , 2 , FALSE , 0 , "Promotor_sequence", "" , "HelpmRNA", "" } ,
  { 0 , 2 , FALSE , 0 , "Premessenger_sequence" , "" , "HelpmRNA", "" } ,
  { 0 , 2 , FALSE , 0 , "Tissues" , "" , "HelpmRNA", "" } ,
  /*
    { 0 , 2 , TRUE , 0 , "mRNA_promotor" ,  "mRNA promotor" , "HelpmRNA", "" } ,
    { 0 , 2 , TRUE , 0 , "mRNA AM promotor" ,  "mRNA AM promotor" , "HelpmRNA", "" } ,
  */

  /*** genome svg plot */
  { 0 , 2 , TRUE , 0 , "Genome_summary" , "" , "HelpGenomeSummary", "" } ,
  { 0 , 2 , FALSE , 0 , "Large_Genome_plot" , "" , "HelpGenomePlot", "" } ,
  { 0 , 2 , TRUE , 0 , "Small_Genome_plot" , "" , "HelpGenomePlot", "" } , 

  /***********************/
  /****** diagrams *******/
  /***********************/
  { 0 , 1 , TRUE , 0 , "", "Diagram" , "" } ,
  { 0 , 85, TRUE , 0 , "DiagToggle" ,  "<b>Show/Hide</b> the clickable diagram on the right side of the screen" , "HelpJan", "" } ,
  { 0 , 86, TRUE , 0 , "DiagFull" ,  "Blows the clickable diagram to <b>full screen</b>, with zoom/move and more details" , "HelpJan", "" } ,
  /* hack: put 102 above 101 so that 102 does not trigger 101 */
  { 0 , 102 , TRUE , 0 , "mrna_DiagHelp" ,  "Help on the mRNA diagram" , "HelpmRNA", "helpGraphicalmRNA" } ,
  { 0 , 101 , TRUE , 0 , "tg_DiagHelp" ,  "Help on the gene diagram" , "" } ,
  
  /***********************/
  { 0 , 0, FALSE, 0, "", "", ""}
  /***********************/
  /***********************/
} ;

void gmpJumpPointInit (void) 
{
  int ii ;
  JUMPPOINTMENU *jmp = jumpPointMenu ;

  for (ii=1, jmp = jumpPointMenu ; jmp->level ; ii++, jmp++)
    {
      jmp->used = 0 ; 
      jmp->index = ii ;
    }
  jumpPointGeneNam[0] = 0 ;
  jumpPointMrnaNam[0] = 0 ;
}

/*****************************************************************/

void gmpJumpPointDestroy (void)
{
  return ;
}

/*****************************************************************/

static void showJP (void) /* for debugging */
{
  JUMPPOINTMENU *jmp ;
 
  for (jmp = jumpPointMenu ; jmp->level ; jmp++)  /* go to last */
    if (jmp->used)
      printf ("index=%d level=%3d: used=%d %s\n", jmp->index, jmp->level, jmp->used, jmp->nam) ;
}

void gmpJumpPointShow (vTXT blkp, GMP *gmp, BOOL showNow, BOOL verbose)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL needFfunc = FALSE ;
  const char *bulle ;
  int nCommand = 0 ;
  vTXT commandBuf = vtxtHandleCreate (h) ;
  
  if (0) showJP () ; /* for computer happiness */

  vtxtPrint (commandBuf,  "var cannedCommands = new Array (\n\"toto\"") ;
  
  /* create a special division that will be dynamically transfered to the tabs frame */
  vtxtPrint (blkp, "<div id='999'>\n") ;
  {
    vtxtPrint (blkp, "<table width=\"98%%\" border=2 bordercolor='blue'>") ;
    {
      bulle = " Table of contents, aims, and hints on how to use this resource" ;
      vtxtPrint (blkp, "  <tr VALIGN=TOP bgcolor=\"#d5d5ff\">\n") ; 
      vtxtPrintf(blkp, "    <td><a title='%s'><font color='black'>?</font></a></td>", bulle) ;
      if (gmp->gene)
	{
	  bulle = "All we learnt about the gene; diagrams, maps, outside links, PubMed" ;
	  vtxtPrintf(blkp, "    <td>") ;
	  nCommand++ ;
	  if (gmp->view == 'g' && gmp->requestedChapter == 1)
	    {
	      vtxtPrint (commandBuf, ", \"openAceViewChapter(0)\"") ;
	      vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)' title='%s'><font color='red'>Gene Summary</font></a>"
			  , nCommand, bulle) ;
	    }
	  else if (gmp->view == 'g')
	    {
	      vtxtPrintf (commandBuf, ", \"openAceViewAction ('gene', '%s', 'fgene')\""
			  , gtJavaScriptProtect(ac_name(gmp->gene), h)) ;
	      vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)' title='%s'>Gene Summary</a>"
			  , nCommand, bulle) ;
	    }
	  else
	     {
	      vtxtPrintf (commandBuf, ", \"openAceViewAction ('gene', '%s', 'fgene')\""
			      , gtJavaScriptProtect(ac_name(gmp->gene), h)) ;
	      vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)' title='%s'>Gene Summary</a>"
			  , nCommand, bulle) ;
	    }

	  vtxtPrint (blkp, "</td>\n") ;
	}
      if (gmp->gene)
	{
	  bulle = "Vertical diagram in true scale, with aligned cDNAs and neighbours." ; 
	  nCommand++ ;
	  vtxtPrintf(blkp, "    <td>") ;
	  vtxtPrintf (commandBuf, ", \"openAceViewAction ('gene', '%s','vgene&v=2&S')\""
		      , gtJavaScriptProtect(ac_name(gmp->gene), h)) ;
	  vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)' title='%s'>Gene on genome </a>", nCommand, bulle) ;
	  vtxtPrint (blkp, "</td>\n") ;
	}

      if (gmp->mrna)
      {
	AC_KEYSET ksm = gmp->gene ? ac_objquery_keyset (gmp->gene, ">transcribed_gene;>mrna", h) :  ac_objquery_keyset (gmp->mrna, "IS *", h) ;
	AC_TABLE tbl = ksm ? ac_keyset_table (ksm, 0, -1, 0, h) : 0 ;
	int ir ;
	AC_OBJ mrna = 0 ;
	char *cq ;
	
	bulle = "Diagram, structure, sequences, putative proteins, cDNA clones, tissues, accessions, ..." ;

	if (tbl->rows > 1)
	  {
	    for (ir = 0 ; tbl && ir < tbl->rows ; ir++, ac_free (mrna))
	      {
		if (!ir)
		  vtxtPrintf (blkp, "    <td><a title='%s'><font color=black>mRNA:</font></a>", bulle) ;
		else
		  vtxtPrint (blkp, ", ") ;
		mrna = ac_table_obj (tbl, ir, 0, h) ;
		cq = gtMrnaSuffix (ac_name (gmp->gene), ac_name(mrna), h) ;
		if (gmp->view == 'm' && ac_obj_equal (mrna, gmp->mrna))
		  {
		    nCommand++ ;
		    vtxtPrint (commandBuf, ", \"openAceViewChapter(0)\"") ;
		    vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)'><font color=red size=+2>%s</font></a>"
				, nCommand, cq) ;
		  }
		else
		  {
		     nCommand++ ;
		     vtxtPrintf (commandBuf, ", \"openAceViewLink ('mrna', '%s')\""
			      , gtJavaScriptProtect(ac_name(mrna), h)) ;
		     vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)' title='%s'>%s</a>"
				 , nCommand, bulle, cq) ;
		  }
	      }
	  }
	else
	  {
	    vtxtPrint (blkp, "    <td>") ;
	    if (gmp->view == 'm')
	      {
		nCommand++ ;
		vtxtPrint (commandBuf, ", \"openAceViewChapter(0)\"") ;
		vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)' title='%s'><font color=red size=+2>%s</font></a>"
			    , nCommand, bulle, "mRNA") ;
	      }
	    else
	      {
		nCommand++ ;
		vtxtPrintf (commandBuf, ", \"openAceViewLink ('mrna', '%s')\""
			    , gtJavaScriptProtect(ac_name(gmp->mrna), h)) ;
		vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)' title='%s'>%s</a>"
			    , nCommand, bulle, "mRNA") ;
	      }
	  }
	  vtxtPrint (blkp, "</td>\n") ;
      }

      if (gmp->gene /* && gmp->view == 'g' */)
	{
	  if (gmp->nMrna > 1)
	    {
	      bulle = "mRNAs, sequences, completeness, alt.introns and exons, proteins, promoters, polyA sites, cDNAs..." ;
	      nCommand++ ;
	      vtxtPrintf(blkp, "    <td>") ;
	      if (gmp->requestedChapter == 2) 
		{
		  vtxtPrint (commandBuf, ", \"openAceViewChapter(0)\"") ;
		  vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)'><font color='red' title='%s'>Alternative mRNAs features, proteins, introns, exons, sequences</font></a>", nCommand, bulle) ;
		}
	      else
		{
		  vtxtPrintf (commandBuf, ", \"openAceViewAction ('gene', '%s', 'fmol')\""
			      , gtJavaScriptProtect(ac_name(gmp->gene), h)) ;
		  vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)' title='%s'>Alternative mRNAs features, proteins, introns, exons, sequences</a>"
			      , nCommand, bulle) ;
		}
	      vtxtPrint (blkp, "</td>\n") ;
	    }
	  
	  vtxtPrintf(blkp, "    <td>") ;
	  bulle = "Accessions aligned in each mRNA, tissues, cDNA clones annotations, choose..." ;
	  nCommand++ ;
	  if (gmp->requestedChapter == 3) 
	    {
	      vtxtPrint (commandBuf, ", \"openAceViewChapter(0)\"") ;
	      vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)'><font color='red' title='%s'>Expression Tissue</font></a>"
			  , nCommand, bulle) ;
	    }
	  else
	    {
	      vtxtPrintf (commandBuf, ", \"openAceViewAction ('gene', '%s', 'fexp')\""
			      , gtJavaScriptProtect(ac_name(gmp->gene), h)) ;
	      vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)' title='%s'>Expression Tissue</a>"
			  , nCommand, bulle) ;
	    }
	  vtxtPrint (blkp, "</td>\n") ;
	  
	  if (gmp->requestedChapter == 4 ||
	      gmp->Spc == WORM ||
	      ac_has_tag (gmp->gene, "Pastille_disease") ||
	      ac_has_tag (gmp->gene, "Pastille_conserved") ||
	      ac_has_tag (gmp->gene, "Interacts") ||
	      ac_has_tag (gmp->gene, "Extern") ||
	      ac_has_tag (gmp->gene, "GO_b_pfam") ||
	      ac_has_tag (gmp->gene, "GO_b_ace") ||
	      ac_has_tag (gmp->gene, "GO_c_pfam") ||
	      ac_has_tag (gmp->gene, "GO_c_ace") ||
	      ac_has_tag (gmp->gene, "GO_c_psort") ||
	      ac_has_tag (gmp->gene, "Pathway") ||
	      ac_has_tag (gmp->gene, "Pfam") ||
	      (gmp->product && ac_has_tag (gmp->product, "AKG"))
	      )
	    needFfunc = TRUE ;
	  else
	    {
	      vTXT bb = vtxtHandleCreate (h) ;
	      ficheNewGenePhenotypeFunctionChapter (bb, gmp) ;
	      if (vtxtPtr (bb))
		needFfunc = TRUE ;
	    }
	  if (needFfunc)
	    {	
	      bulle = "Disease, Pathways, Process, Localization, Motifs, RIF, Related genes, Systems biology..." ;      
	      vtxtPrintf(blkp, "    <td>") ; 
	      nCommand++ ;
	      if (gmp->requestedChapter == 4) 
		{
		  vtxtPrint (commandBuf, ", \"openAceViewChapter(0)\"") ;
		  vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)' title='%s'><font color='red'>Function, regulation, related genes </font></a>"
			      , nCommand, bulle) ;		  
		}
	      else
		{
		  vtxtPrintf (commandBuf, ", \"openAceViewAction ('gene', '%s', 'ffunc')\""
			      , gtJavaScriptProtect(ac_name(gmp->gene), h)) ;
		  vtxtPrintf (blkp, "<a href='javascript:fireCommand (%d)' title='%s'>Function, regulation, related genes </a>", nCommand, bulle) ;
		}
	      if (ac_has_tag (gmp->gene, "Pastille_disease"))
		{
		  /* graphColor (CERISE) ; */
		  vtxtPrintf (blkp, "<font color='red'>D</font>") ;
		}
	      if (ac_has_tag (gmp->gene, "Pastille_conserved"))
		{
		  /* graphColor (PURPLE) ; */
		  vtxtPrintf (blkp, "<font color='purple'>C</font>") ; 
		}
	      if (ac_has_tag (gmp->gene, "Interacts"))
		{
		/* graphColor (DARKGREEN) ; */
		  vtxtPrintf (blkp, "<font color='#009c3a'>I</font>") ; 
		}
	      /* graphColor (BLACK) ; */

	      vtxtPrint (blkp, "</td>\n") ;
	    }
	}

      vtxtPrint (blkp, "</tr>\n") ;
     }
    vtxtPrint (blkp, "</table>") ;
  }
  vtxtPrint (blkp,"</div>\n") ;

  /* create a special division that will be dynamically transfered to the tabs frame */

  if (nCommand)
    {
      vtxtPrint (commandBuf, ") ;\n") ;
      vtxtPrint (blkp, "<div id='executor' onclick='jsExecutor()'>\n") ;
      
      vtxtPrint (blkp, "\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
		 "  <!--\n"
		 ) ;
      vtxtPrint (blkp, vtxtPtr (commandBuf)) ;
      vtxtPrint (blkp,"  //-->\n"
		 "</script>\n"
		 ) ;	   
      
      vtxtPrint (blkp, "</div>\n") ;
    }

  if (verbose)
    {
      vtxtEmptyLine (blkp, 1) ;
        
      if (! ac_has_tag (gmp->gene, "Cloud_gene"))
	{
	  vtxtPrint (blkp, "To <span class='ace_summary'>mine knowledge</span> about the gene, please click the ") ;
	  vtxtPrintf (blkp, "<span class='ace_summary'>'Gene Summary'</span>") ;
	  
	  vtxtPrint (blkp, " or the <span class='ace_summary'>'Function, regulation, related genes '</span>") ;
	  vtxtPrint (blkp, " tab at the top of the page. The ") ;
	  vtxtPrintf (blkp, "<a href='javascript:openAceViewAction (\"gene\", \"%s\", \"fgene\")'>"
		      , ac_name(gmp->gene)) ;
	  vtxtPrintf (blkp,"'Gene Summary'</a>") ;
	  vtxtPrint (blkp, " page includes all we learnt about the gene, functional annotations of neighboring genes, maps, links to other sites and the bibliography. The ") ;
	  vtxtPrintf (blkp, "<a href='javascript:openAceViewAction (\"gene\", \"%s\", \"ffunc\")'>"
		      , ac_name(gmp->gene)) ;
	  vtxtPrintf (blkp,"'Function, regulation, related genes '</a>") ;
	  vtxtPrintf (blkp," page includes Diseases (D), Pathways, GO annotations, conserved domains (C), interactions (I) reference into function, and pointers to all genes with the same functional annotation") ;
	}
      else
	{
	  vtxtPrint (blkp, "To <span class='ace_summary'>mine knowledge</span> about the gene, please click the ") ;
	  vtxtPrintf (blkp, "<a href='javascript:openAceViewAction (\"gene\", \"%s\", \"fgene\")'>"
		      , ac_name(gmp->gene)) ;
	  vtxtPrintf (blkp,"'Gene Summary'</a>") ;
	  
	  vtxtPrint (blkp, " tab at the top of the page which includes functional annotations of neighboring genes, maps, links to other sites") ;
	}
      
      vtxtBreak (blkp) ;
      
      if (gmp->nMrna > 1) 
	{
	  vtxtPrint (blkp, "To <span class='ace_summary'>compare alternative variants</span>, their summarized annotations, predicted proteins, introns and exons, or to access any sequence, click the ") ;
	  vtxtPrintf (blkp, "<a href='javascript:openAceViewAction (\"gene\", \"%s\", \"fmol\")'>"
		      , ac_name(gmp->gene)) ;
	  vtxtPrintf (blkp,"'Alternative mRNAs features'</a>") ;
	  
	  vtxtPrint (blkp, " tab. To see a <span class='ace_summary'>specific mRNA variant</span> diagram, sequence and annotation, click the variant name in the ") ;
	  vtxtPrintf (blkp, "<a href='javascript:openAceViewAction (\"mrna\", \"%s\", \"fiche\")'>"
		      , ac_name(gmp->mrna)) ;
	  vtxtPrintf (blkp,"'mRNA'</a> tab") ;
	}
      else
	{
	  vtxtPrint (blkp, "To see the <span class='ace_summary'>mRNA diagram, sequence and annotation</span>,  click the ") ;
	  vtxtPrintf (blkp, "<a href='javascript:openAceViewAction (\"mrna\", \"%s\", \"fiche\")'>"
		      , ac_name(gmp->mrna)) ;
	  vtxtPrintf (blkp,"'mRNA'</a> tab") ;
	}
      
      vtxtPrint (blkp, ". To examine <span class='ace_summary'>expression data</span> from all cDNAs clustered in this gene by AceView, click the ") ;
      vtxtPrintf (blkp, "<a href='javascript:openAceViewAction (\"gene\", \"%s\", \"fexp\")'>"
		  , ac_name(gmp->gene)) ;
      vtxtPrintf (blkp,"'Expression tissue'</a>") ;
      
      vtxtBreak (blkp) ;
      
      vtxtEmptyLine (blkp, 1) ;
      
      vtxtPrint (blkp, "If you know more about this gene, or found errors, please ") ;
      vtxtPrint (blkp, "<a href=\"mailto:mieg@ncbi.nlm.nih.gov\">share</a>") ;
      vtxtPrint (blkp, " your knowledge. Thank you ! ") ;
    }

  if (showNow) /* show now, otherwise, main.js will show after  receiving the remainder file */
    vtxtPrint (blkp, "\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
	       "  <!--\n"
	       "    aceViewOnLoad() ;\n "
	       "  //-->\n"
	       "</script>\n"
	       ) ;
  
  ac_free (h) ;
  return ;
}

/*****************************************************************/

static JUMPPOINTMENU *gmpJumpPoint (vTXT blkp, GMP *gmp, const char *header)
{
  JUMPPOINTMENU *jmp = jumpPointMenu ;
  
  if (!header || !gmp->markup)
    return 0 ;

  if (gmp->view == 'g' && gmp->gene)
    strncpy (jumpPointGeneNam, ac_name (gmp->gene), 1023) ;
  if (gmp->view == 'm')
    strncpy (jumpPointMrnaNam, ac_name (gmp->mrna), 1023) ;

  for (jmp = jumpPointMenu ; jmp->level ; jmp++)
    if (!strcmp (jmp->nam, header))
      {
	jmp->used |= 1 ; /* |= 2: chapterClose, |= 2: blockStart */
	return jmp ;
      }
  return 0 ;
}

/********************************************************************/
/* the index can be used to open/close a section */
int gmpJumpPointIndex (char *header)
{
  JUMPPOINTMENU *jmp = jumpPointMenu ;
  
  if (!header)
    return 0 ;

  for (jmp = jumpPointMenu ; jmp->level ; jmp++)
    if (!strcmp (jmp->nam, header))
      return jmp->index ;

  return 0 ;
}

/********************************************************************/
/********************************************************************/
#ifdef JUNK
static BOOL isVoyelle (char c1)
{
  int c = ace_lower (c1) ;
  switch (c)
    {
    case 'a':    case 'e':    case 'i':    case 'o':    case 'u':
      return TRUE ;
    }
  return TRUE ;
}
#endif
/********************************************************************/

static char* gtDoCleanUp (const char *ptr, int doLower)
{
  char *cp, *cq ; 
  int nn, n, doIt = 1 ;
  static char *buf = 0 ;

  if (!ptr || !*ptr) return "" ;

  cp = buf ; /* we do not lose ptr if ptr results from a previous call */
  buf = strnew (ptr, 0) ;
  if (cp) messfree (cp) ;

  cp = buf - 1 ;
  
  while (*++cp)
    switch (*cp)
      {
      case '_': *cp = ' ' ; break ;
      case '\"': *cp = '\'' ; break ;
      } 

  /* clean initial spaces */
  cp = cq = buf ;
  while (*cp == ' ') cp++ ;
  n = cp - cq ;
  if (n)
    while ((*cq++ = *cp++)) ;

  /* clean terminal punctuation */
  cp = buf + strlen (buf) ;
  cp-- ;
  while (cp > buf && (*cp == ' ' || *cp == ',' || *cp == '.')) *cp-- = 0 ;
  
  /* clear initial upper case, if first word is totally lower case */
  doIt = 1 ;
  
  cp = buf ; nn = *cp ? 1 : 0 ;
  while (doIt == 1)
    switch (*++cp)
      {
      case 0:
      case ' ':
      case ',':
      case ';':
      case '.':
      case '-':
      case '+':
      case '_':
      case '/':
	doIt = 2 ;
	break ;
      default:
	nn++ ;
	if (*cp != ace_lower (*cp))
	  doIt = 0 ;
	break ;
      }
  if (doLower && nn > 1 && doIt == 2)
    { cp = buf ; *cp = ace_lower (*cp) ; } 

  if (doLower == 2)
    {
      cp = cq = buf ;
      while (*cq)
	{
	  cp = cq - 1 ; doIt = 0 ;
	  while (*++cp)
	    {
	      if (!*cp || *cp == ' ' || *cp == ',')
		break ;
	      if (*cp != ace_upper (*cp))
		doIt = 1 ;
	    }
	  if (doIt)
	    {
	      cp = cq - 1 ;
	      while (*++cp)
		if (!*cp || *cp == ' ' || *cp == ',')
		  break ;
		else
		  *cp = ace_lower (*cp) ;
	    }
	  cq = cp ;
	  while (*cq == ' ' || *cq == ',') cq++ ;
	}
    }
  if (doLower == 3)
    { /* renove surrounding [] if they exist */
      cp = cq = buf ;
      if (*cp == '[')
	{
	  while (*cq++) ;
	  cq-- ;
	  if (*cq == ']')
	    {
	      *cq = 0 ;
	      cq = cp + 1 ;
	      while ((*cp++ = *cq++)) ;
	    }
	}
    }
  
  return buf ;
} /* gtDoCleanUp */

/****/
char* gtCleanQuotes (const char *ptr) /* keep all upper/lower */
{
  return gtDoCleanUp (ptr, 0) ;
}
/****/
char* gtCleanUp (const char *ptr)
{   /* change first letter to upper except if full Up*/
  return gtDoCleanUp (ptr, 1) ;
}
/****/
char* gtLowerCleanUp (const char *ptr)
{ /* change first letter to lower except if some other letter is Up */
  return gtDoCleanUp (ptr, 2) ;
}
/****/
char* gtSquareCleanUp (const char *ptr)
{   /* change first letter to upper except if full Up*/
  return gtDoCleanUp (ptr, 3) ;
}

/****/
static char* gtJavaScriptProtect (const char *cq, AC_HANDLE h)
{   /* change first letter to upper except if full Up*/
  char cc, *cp, *cp0 = halloc (3 * strlen (cq), h) ;
  char b[3];

  cp = cp0 ; 
  while ((cc = *cq++))
    {
      if (isalnum (cc) || cc ==  '_' || cc == ',' || cc == '.' || cc == '-' )
	*cp++ = cc ;
      else
	{
	  sprintf(b,"%02x",cc);
	  *cp++ = '%';
	  *cp++ = b[0];
	  *cp++ = b[1];
	}
    }
  *cp++ = 0 ;      
  return cp0 ;
}

/********************************************************************/
/* set initial upper case if first word is totally lower case */
char* gtSetUpper (const char *buf)
{
  char *cp ;
  int doIt = 1 ;
  static char *buf1 = 0 ;

  if (buf && *buf)
    {
       messfree (buf1) ;
      buf1 = strnew (buf, 0) ;
      cp = buf1 - 1;
      while (doIt == 1)
	switch (*++cp)
	  {
	  case 0:
	  case ' ':
	  case ',':
	  case ';':
	  case '.':
	    doIt = 2 ;
	    break ;
	  default:
	    if (*cp != ace_lower (*cp))
	      doIt = 0 ;
	    break ;
	  }
      if (doIt == 2)
	{ cp = buf1 ; *cp = ace_upper (*cp) ; } 
    }
  return buf1 ;
}

/********************************************************************/
/* set whole buf to upper case */
char* gtSetAllUpper (const char *buf)
{
  char *cp ;
  static char *buf1 = 0 ;

  if (buf && *buf)
    {
      messfree (buf1) ;
      buf1 = strnew (buf, 0) ;
      cp = buf1 - 1;
      while (*++cp)
	*cp = ace_upper (*cp) ;
    }
  return buf1 ;
}

/********************************************************************/
/* renormalize y1L, yc, yg, GGB:, GB: names */
char* gtYkName (const char *old)
{
  static char buff[200] ;
  char *acc = buff ;
  const char *ptr  = old ;
    /* renormalize yc yd y1L  to yk */
 if (*ptr=='y')
    {
      if (*(ptr+1)=='k')
	strcpy (buff, old) ;
      else if (*(ptr+1) == 'c' || *(ptr+1)=='d')
	{
	  strcpy(buff,"yk") ;
	  strcat(buff,ptr+2) ;
	}
      else if (*(ptr+2) == 'L')
	{
	  strcpy(buff,"yk") ;
	  strcat(buff,ptr+3) ;
	  acc = buff ;
	}
    }
  else
    {
      if (!strncmp (old,"GGB:",4)) 
	strcpy (buff, old+4) ;
      else if (!strncmp (old,"GB:",3)) 
	strcpy (buff, old+3) ;
      else
	strcpy (buff, old) ;
    }
  
  return acc ;
} /* gtYkName */

/***********************/

BOOL gtIsPseudogene (AC_OBJ oGene)
{
  BOOL ok = FALSE ;
  AC_OBJ oPg = 0, oMethod = 0 ;

  if (!ac_has_tag (oGene,"Transcribed_gene") &&
      (oPg = ac_tag_obj (oGene, "Genefinder", 0)) &&
      (oMethod = ac_tag_obj (oPg,"Method", 0)) &&
      strstr (ac_name(oMethod),"seudogene"))
    ok = TRUE ;
  ac_free (oPg) ; ac_free (oMethod) ;

  return ok ;
} /* gtIsPseudogene */

/***********************/

BOOL gtIsCloudGene (AC_OBJ oGene)
{
  BOOL ok = FALSE ;

  if (ac_has_tag (oGene,"Cloud_gene"))
    ok = TRUE ;

  return ok ;
} /* gtIsCloudGene */

/*********/

BOOL gtIsEssentialGene (AC_OBJ oGene)
{
  BOOL ok = FALSE ;
  AC_TABLE tbl = 0 ;
  const char *ptr ;
  int ir ;

  if ((tbl = ac_tag_table (oGene, "Locus_description", 0)))
    for (ir = 0 ; ir < tbl->rows && tbl->cols >= 1; ir++)
      {
	ptr = ac_table_printable (tbl, ir, 0, 0) ;
	if (ptr && *ptr && strstr(ptr,"ssential gene")) 
	  ok = TRUE ;
      }

  ac_free (tbl) ;
  return ok ;
} /* gtIsEssentialGene */

/*********/

static BOOL gtIsPrecursor (AC_OBJ oProduct, AC_OBJ oGF)
{
  BOOL ok = FALSE ;
  AC_TABLE tbl = 0 ;
  const char *ptr ;
  int ir ;

  if (
      (oProduct && (tbl = ac_tag_table (oProduct, "Kantor_title", 0)))||
      (oGF &&  (tbl = ac_tag_table (oGF, "has_est_kantor_title", 0)))
      )
    for (ir = 0 ; ir < tbl->rows && tbl->cols >= 1; ir++)
      {
	ptr = ac_table_printable (tbl, ir, 0, 0) ;
	if (ptr && strstr(ptr,"likely secreted or extracellular")) 
	  ok = TRUE ;
      }
  
  ac_free (tbl) ;
  return ok ;
} /* gtIsPrecursor */

/*********/

static BOOL gtHasLocus (AC_OBJ oGene)
{
  return  oGene && 
    (ac_has_tag (oGene, "Locus")) ?
    TRUE : FALSE ;
} /* gtHasLocus */

/*********/

static BOOL gtHasEst (AC_OBJ oGene, AC_OBJ oGF)
{
  return  
    (oGene && ac_has_tag (oGene, "Transcribed_gene")) || 
    (oGene && ac_has_tag (oGene, "Has_cDNA_clone")) || 
    (oGF && ac_has_tag (oGF, "has_est")) ?
    TRUE : FALSE ;
} /* gtHasEst */ 

/*********/

static BOOL gtHasOst (AC_OBJ oGene)
{
  AC_HANDLE h = handleCreate () ;
  BOOL ok = FALSE ;
  AC_TABLE tbl, oPh ;
  AC_OBJ oR ;
  int ir, jr ;
  const char *ptr ;
  
  if (0) return FALSE ; /* not yet public */
  
  if (oGene && (tbl = ac_tag_table (oGene, "RNAi", h)))
    for (ir = 0 ; ir < tbl->rows && tbl->cols >= 1; ir++)
      {
	oR = ac_table_obj (tbl, ir, 0, h) ;
	if (!oR) continue ;
	if (strncmp("mv_", ac_name(oR), 3))
	  continue ; 
	if ((oPh = ac_tag_table (oR, "Phenotype", h)))
	  for (jr = 0 ; !ok && jr < oPh->rows && oPh->cols >= 1; jr++)
	    if ((ptr = ac_table_printable (oPh, jr, 0, 0)) &&
		!strcasecmp (ptr, "amplified"))
	      ok = TRUE ;
      }
  ac_free (h) ;
  
  return ok ;
} /* gtHasOst */

/*********/

static BOOL gtHasExpression (AC_OBJ oGene)
{
  AC_HANDLE h = 0 ;
  BOOL ok = FALSE ;
  AC_TABLE tbl, gg ;
  AC_OBJ oR ;
  int ir ;
  
  if (!oGene)
    return FALSE ;
  
  if (ac_has_tag (oGene, "Locus_phenotype") || 
      ac_has_tag (oGene, "Locus_description") || 
      ac_has_tag (oGene, "Strain") || 
      ac_has_tag (oGene, "Expr_pattern") || 
      ac_has_tag (oGene, "Pattern")
      )
    return TRUE ;

  h = handleCreate () ;
  if (oGene && (tbl = ac_tag_table (oGene, "RNAi", h)))
    for (ir = 0 ; !ok && ir < tbl->rows && tbl->cols >= 1; ir++)
      {
	oR = ac_table_obj (tbl, ir, 0, h) ;
	if (!oR) continue ;
	if (!strncmp("mv_",ac_name(oR),3))
	  continue ; 
	gg = ac_tag_table (oR, "Gene", h) ;
	if (gg && gg->rows < 2 && ac_has_tag (oR, "Phenotype"))	 
	  ok = TRUE ;
      }

  ac_free (h) ;
  return ok ;
} /* gtHasExpression */

/*********/

BOOL gtIsExperimentalGene (AC_OBJ oGene, AC_OBJ oGF)
{
  return gtHasEst (oGene, oGF) ;
} /* gtGeneIsExperimental */

/********************************************************************/

char* gtCloneGroup (AC_OBJ oGene, BOOL isLink)
{
  AC_HANDLE h = handleCreate () ;
  AC_OBJ oClo, oTg ;
  AC_TABLE gClg, gTg ;
  const char *ccp = 0 ;
  int nn, ir1, ir2, nMax, nHit ;
  static char buf[12] ;

  oClo = ac_tag_obj (oGene, "Clone_group", h) ;
  if (!oClo && (gTg = ac_tag_table (oGene, "Transcribed_gene", h)))
    for(nMax = ir1 = 0 ; ir1 < gTg->rows; ir1++)
      {
	oTg = ac_table_obj (gTg, ir1, 0, h) ;
	gClg = oTg ? ac_tag_table (oTg, "Clone_group", h) : 0 ;
	if (gClg)
	  for(ir2 = 0 ; ir2 < gClg->rows; ir2++)
	    {
	      nHit = ac_table_int (gClg,ir2,1, -1) ;
	      if (nHit > nMax)
		{ nMax = nHit ; oClo = ac_table_obj (gClg,ir2, 0, h) ; }
	    }
      }
  if (!oClo && (gTg = ac_tag_table (oGene, "Has_cdna_clone", h)))
    for(nMax = -1,  ir1 = 0 ; ir1 < gTg->rows; ir1++)
      {
	oTg = ac_table_obj (gTg, ir1, 0, h) ;
	gClg = oTg ? ac_tag_table (oTg, "Clone_group", h) : 0 ;
	if (gClg)
	  for(ir2 = 0 ; !oClo && ir2 < gClg->rows; ir2++)
	    {
	      oClo = ac_table_obj (gClg,ir2, 0, h) ; 
	    }
      }
  *buf = 0 ;
  if (oClo && (ccp = ac_name (oClo)) && (nn = atoi(ccp+2)))
    {
      if (isLink)
	sprintf (buf, "CELK%05d",nn) ;
      else
	sprintf (buf, "YK%d",nn) ;
    }

  ac_free (h) ;
  return *buf ? buf : 0 ;
} /* gtCloneGroup */

/********************************************************************/

char* gtGeneCard (AC_OBJ oGene)
{
  static char buf[128] ;
  const char *ccp ;
  char *cp ;

  if ((ccp = ac_tag_printable (oGene, "GeneCard_id", 0)))
    {
      sprintf (buf, "GeneCard:%s", ccp) ;
      cp = buf+9 ;
      while (*cp) { *cp = ace_upper(*cp) ; cp++ ; }
      return buf ;
    }

  return 0 ;
} /* gtGeneCard */

/********************************************************************/

char* gtOst (AC_OBJ oGene)
{
  AC_HANDLE h = handleCreate () ;
  AC_TABLE tbl ;
  AC_OBJ oRnai ;
  const char *ccp = 0, *phe ;
  static char buf[64] ;
  int ir ;
  
  if ((tbl = ac_tag_table (oGene, "RNAi", h)))
    for(ir = 0 ; ir < tbl->rows && tbl->cols >= 1 ; ir++)
      {
	ccp = ac_table_printable (tbl, ir, 0, "") ;
	if (!strncmp ("mv_", ccp, 3))
	  {
	    strncpy (buf, ccp+3, 63) ;
	    oRnai = ac_table_obj (tbl, ir, 0, h) ;
	    phe = ac_tag_text (oRnai, "Phenotype", 0) ;
	    if (phe && *phe)
	      { ac_free (h) ; return buf ; }
	  }
      }
  ac_free (h) ;
  return 0 ;
} /* gtOst */

/********************************************************************/

char* gtHyman (AC_OBJ oGene)
{
   AC_HANDLE h = handleCreate () ;
  AC_TABLE tbl ;
  AC_OBJ oR ;
  static char buf[64] ;
  int ir ;
  
  if ((tbl = ac_tag_table (oGene, "RNAi", h)))
    for(ir = 0 ; ir < tbl->rows && tbl->cols >= 1 ; ir++)
      {
	oR = ac_table_obj (tbl, ir, 0, h) ; 
	strncpy (buf, ac_name (oR) + 3, 63) ;
	if (!strncasecmp (ac_name (oR), "TH", 2) && ac_has_tag (oR, "DIC_Movie_available"))
	  { ac_free (h) ; return buf ; }
      } 
  ac_free (h) ;
  return 0 ;
} /* gtHyman */

/********************************************************************/

char* gtPiano (AC_OBJ oGene)
{
  AC_HANDLE h = handleCreate () ;
  AC_TABLE tbl ;
  AC_OBJ oR ;
  int ir ;
  
  if ((tbl = ac_tag_table (oGene, "RNAi", h)))
    for(ir = 0 ; ir < tbl->rows && tbl->cols >= 1 ; ir++)
      {
	oR = ac_table_obj (tbl, ir, 0, h) ; 
	if (strncasecmp (ac_name (oR), "FP", 2) && 
	    ac_has_tag (oR, "DIC_Movie_available"))
	  { ac_free (h) ; return "FP" ; }
      } 
  ac_free (h) ;
  return 0 ;
} /* gtPiano */

/***********************/

DICT *gtGeneAliases (GMP *gmp, BOOL isLink)
{
  AC_HANDLE h = handleCreate () ;
  AC_TABLE tLoc, tGeneId, tGenefinder, tLocLink, tNewNameOld ;
  int ir ;
  const char *ccp ;
  char *cp, buf[1000]   ;
  AC_OBJ oGene = gmp->gene, gf ;
  DICT *dict = dictCreate (12) ;

  tLoc = ac_tag_table (oGene,"Locus", h) ;
  tLocLink = ac_tag_table (oGene,"LocusLink", h) ;
  tGeneId = ac_tag_table (oGene,"GeneId", h) ;
  tGenefinder = ac_tag_table (oGene,"Genefinder", h) ;
  tNewNameOld = ac_tag_table (oGene,"NewNameOld", h) ;

  dictAdd (dict, ac_name (oGene), 0) ;

  if (tLoc)
    for (ir = 0 ; ir < tLoc->rows ; ir++)
      if ((ccp = ac_table_printable (tLoc, ir, 0, 0)))
	dictAdd (dict, ccp, 0) ;

  if (tLocLink)
    for (ir = 0 ; ir < tLocLink->rows ; ir++)
      if ((ccp = ac_table_printable (tLocLink, ir, 0, 0)) &&
	  strncmp (ccp, "LOC", 3) &&
	  strncmp (ccp, "OTTHUM", 6) &&
	  strncmp (ccp, "LocusID_", 8))
	dictAdd (dict, ccp, 0) ;

  if (gmp->Spc != WORM && tGeneId)
    for (ir = 0 ; ir < tGeneId->rows ; ir++)
      if ((ccp = ac_table_printable (tGeneId, ir, 0, 0)) &&
	  strncmp (ccp, "LOC", 3) &&
	  strncmp (ccp, "OTTHUM", 6) &&
	  strncmp (ccp, "LocusID_", 8))
	dictAdd (dict, messprintf ("LOC%s", ccp), 0) ;

  if (tNewNameOld)
    for (ir = 0 ; ir < tNewNameOld->rows ; ir++)
      if ((ccp = ac_table_printable (tNewNameOld, ir, 0, 0)))
	dictAdd (dict, ccp, 0) ;

  if (tGenefinder)
    for(ir = 0 ; ir < tGenefinder->rows; ir++)
      { 
	strcpy (buf, ac_table_printable (tGenefinder,ir,0, "")) ;
	if (gmp->Spc == WORM && !isLink)
	  {
	    cp = buf + strlen(buf) - 1 ;
	    if (*cp >= 'a' && *cp <= 'z') *cp = 0 ;
	  }
	if (gmp->Spc != WORM &&
	    strstr (buf, ac_name(gmp->gene)))
	  continue ; /* throw away _CACNA1I_1,  _CACNA1I_2 */
      
	if (gmp->Spc != WORM &&
	    *buf == '_')
	  {
	    char buf1[100] ;
	    BOOL ok = TRUE ;
	    int i ;

	    strncpy (buf1, buf+1,  99) ;
	    cp = strstr (buf1, "_") ;
	    if (cp) *cp = 0 ;
	    for (i = 1 ; ok && i <= dictMax (dict) ; i++)
	      if (!strcmp (buf1, dictName (dict, i)))
		ok = FALSE ; 
	    if (!ok)
	      continue ; /* throw away _BAGE_4 */
	  }

	if (*buf == '_')
	  dictAdd (dict, buf+1, 0) ;
	else if (*buf)
	  dictAdd (dict, buf, 0) ;
      }
  else
    {
      tGenefinder = ac_tag_table (gmp->tg,"Matching_genefinder_gene", h) ;
      if (tGenefinder)
	for(ir = 0 ; ir < tGenefinder->rows; ir++)
	  { 
	    gf = ac_table_obj (tGenefinder,ir,0, h) ;
	    if (!gf) continue ;
	    ccp = ac_tag_printable (gf, "GeneId_pg", 0) ;
	    if (!ccp)
	      ccp = ac_tag_printable (gf, "LocusId", 0) ;
	    if (!ccp)
	      continue ;
	    cp = strnew (ccp, h) ;
	    if (!isLink)
	      {
		char *cq = cp + strlen(cp) - 1 ;
		if (*cq >= 'a' && *cq <= 'z') *cq = 0 ;
	      }
	    if (*cp == '_') cp++ ;
	    if (!strncasecmp(cp,"LocusId:",7)) /* a hack to compensate a bad parsing in human-34 */
	      dictAdd (dict, messprintf("GeneId:%s",cp+7), 0) ;
	    else
	      dictAdd (dict, messprintf("GeneId:%s",cp), 0) ;
	  }
    }

  if (!isLink &&
      gmp->Spc == WORM &&
      (ccp = ac_tag_printable (oGene,"NewName",0))) /* do not show in ficheasn */
    dictAdd (dict, ccp, 0) ;

  if ((ccp = gtCloneGroup (oGene, isLink)))
    dictAdd (dict, ccp, 0) ;
  
  if (gmp->Spc == HUMAN && (ccp = gtGeneCard (oGene)))
    dictAdd (dict, ccp, 0) ;

  if (!dictMax (dict))
    dictDestroy (dict) ;

  ac_free (h) ;
  return  dict ;
} /* gtGeneAliases  */

/********************************************************************/
/********************************************************************/

static void gtMolecularWeight (vTXT blkp, AC_OBJ oProduct, AC_OBJ oGF, BOOL showPi)
{
  float pI = -200.0, mw = 0 ;
  AC_OBJ oP = 0, oM = 0 ;
  AC_TABLE tbl = 0 ;
  AC_HANDLE h = handleCreate () ;

  showPi = 0 ;
  if (oProduct)
    oP = oProduct ;
  else if (oGF &&
	   (oM = ac_tag_obj (oGF, "Predicted_mrna", h)))
    oP = gtMrna2Product (oM, h) ;
  
  if (oP &&
      (tbl = ac_tag_table (oP, "Molecular_weight", h)) &&
      tbl->cols >= 2 && 
      (mw = ac_table_float(tbl, 0, 0, 0)) && 
      (pI = ac_table_float(tbl, 0, 1, 0)) &&
      pI > -100
      ) 
    {
      if (showPi)
	vtxtPrintf (blkp, " (%.1f kDa, pI %.1f)", mw, pI) ;
      else
	vtxtPrintf (blkp, " (%.1f kDa)", mw) ;
    }
  ac_free (h) ;
} /* gtMolecularWeight */

/********************************************************************/

static char* gtGeneBaseTitle (vTXT blkp, GMP *gmp, BOOL isPutative, BOOL isShedded)
{
  AC_HANDLE h = handleCreate () ;
  AC_OBJ oPg = 0 ;
  BOOL isComplex = FALSE ;
  AC_TABLE tbl ;
  int ir, nn = 0 ;

  if (gtIsEssentialGene (gmp->gene))
    { nn++ ; vtxtPrintf (blkp, " essential") ; }
  else if (!gtHasEst (gmp->gene, gmp->pg) && 
	   !gtHasOst (gmp->gene) && 
	   !gtHasExpression (gmp->gene) &&
	   (oPg = ac_tag_obj (gmp->gene, "GeneFinder", h)) &&
	   ac_has_tag (oPg, "CDS") /* meaningless for tRNA etc */
	   )
    { nn++ ; vtxtPrintf (blkp, " predicted") ; isPutative = FALSE ; }

  if (gtIsCloudGene (gmp->gene))
    { nn++ ; vtxtPrintf (blkp, " cloud") ; }

  if (ac_has_tag (gmp->gene, "Complex_locus"))
    {
      isPutative = FALSE ; vtxtPrintf (blkp, " complex locus ") ; isComplex = TRUE ;
    }
  if (isPutative)
    { /* show only if nothing better and we did not already say 'predicted' */
      if (isShedded)
	vtxtPrintf (blkp, " variant ") ;
      else if (gtIsCloudGene (gmp->gene))
	vtxtPrintf (blkp, " gene ") ;
      else if (ac_has_tag (gmp->gene, "Single_exon_gene") && ac_has_tag (gmp->gene, "Pastille_coding"))
	vtxtPrint (blkp, "  single exon coding gene ") ;
      else if (ac_has_tag (gmp->gene, "Single_exon_gene"))
	vtxtPrint (blkp, "  gene ") ;
      else  if (ac_has_tag (gmp->gene, "Pastille_spliced_non_coding"))	   
	vtxtPrintf (blkp, " spliced non-coding gene ") ;
      else  if (ac_has_tag (gmp->gene, "Spliced_gene") && ac_has_tag (gmp->gene, "Pastille_coding"))	   
	vtxtPrintf (blkp, " spliced coding gene ") ;

      if (!gtIsEssentialGene (gmp->gene) &&
	  ac_has_tag (gmp->gene, "Locus_description"))
	vtxtPrint (blkp, "with phenotype ") ;
    }
  else if (!isComplex)
    vtxtPrint (blkp, " gene ") ;
  
  if (gmp->markup && (tbl = ac_tag_table (gmp->gene, "GeneId", h)))
    {
      char linkBuf[4096] ;
      vTXT buf = vtxtHandleCreate (h) ;

      for (ir = 0 ; ir  < tbl->rows ; ir++)
	vtxtPrintf (buf, "%s%s", ir ? "," : "", ac_table_printable (tbl, ir, 0, "")) ;
      sprintf (linkBuf, ENTREZ_LINK, vtxtPtr (buf)) ;

      gmpURL (blkp, gmp, linkBuf, ac_name(gmp->gene)) ;
    }
  else
    vtxtPrint (blkp, ac_name(gmp->gene)) ; 
  

  ac_free (h) ;
  return vtxtPtr (blkp) ;
} /* gtGeneBaseTitle */

/********************************************************************/
/* type could be Genetic or molecular */

static void gtOneGeneClassDescriptor (vTXT blkp, AC_OBJ oGene, AC_OBJ oLocus, const char *ptr, int nn, BOOL isGene)
{
  char *cp ;
  int len ;

  vtxtPrintf (blkp, "%s%s", nn > 0 ? ", " : "", ptr) ;
  if (oLocus && 
      !(isGene && !strcasecmp(ac_name(oGene), ac_name (oLocus))))
    {
      len = vtxtLen (blkp) ;
      vtxtPrintf (blkp, " %s", ac_name (oLocus)) ; 
      cp = vtxtPtr (blkp) + len ;
      if (!isGene) vtextUpperCase (cp) ;
    }
}

/********/

static char *gtGeneClassDescriptor (vTXT blkp, AC_OBJ oGene, const char *shortKantorTitle, char *descriptor, BOOL isGene, BOOL *isPutativep)
{
  AC_HANDLE h = handleCreate () ; 
  AC_TABLE tbl ;
  AC_OBJ oLocus, oGeneClass ;
  int jj, ir, nn = 0, n0, n00, n2, n4, n5 ;
  const char *ptr ;
  DICT *dict = 0 ;
  BOOL isMol = 0, doClean = FALSE, isEssential = gtIsEssentialGene (oGene) ;
  const char *molecDescription ;

  if (isGene || ! ac_has_tag (oGene, "Complex_locus"))
    molecDescription = descriptor ? descriptor : shortKantorTitle ; 
  else /* mrna or product from complex gene */
    molecDescription = (! *isPutativep && shortKantorTitle) ? shortKantorTitle : descriptor ;

  if (molecDescription)
    {
      dict = dictHandleCreate (16, h) ;
      dictAdd (dict, molecDescription, 0) ;
    }
  vtxtClear (blkp) ;

  for (n0 = n00 = n2 = n4 = n5 = -1, jj = 0 ; jj < 6 ; jj++)
    {
      if ((tbl = ac_tag_table (oGene, "Locus", h)))
	for (ir = 0 ; ir < tbl->rows && tbl->cols >= 1; ir++)
	  {
	    oLocus = ac_table_obj (tbl, ir, 0, h) ;
	    if (!(oGeneClass = ac_tag_obj (oLocus, "Gene_class", h)))
	      continue ;
	    isMol = ac_has_tag (oGeneClass, "Molecular") ;
	    ptr = ac_tag_printable (oGeneClass, "Description", 0) ;
	    if (isMol && descriptor)
	      continue ; /* would be redundant */
	    if (!dict)
	      dict = dictHandleCreate (16, h) ;
	    if (dictFind (dict, ac_name(oGeneClass), 0))
	      continue ;
	    if (dictFind (dict, ptr, 0))
	      continue ;
	    switch (jj)
	      { 
	      case 0: /* true molec or wb id */
		if (n0 < 0 && isEssential && isMol && ptr &&
		    strstr (ac_name(oGene), ac_name(oLocus))
		    )
		  { 
		    if (*isPutativep && doClean)
		      { vtxtClear (blkp) ; nn = 0 ; }
		    *isPutativep = FALSE ;
		    dictAdd (dict, ac_name(oGeneClass), 0) ;
		    n0 = ir ; gtOneGeneClassDescriptor (blkp, oGene, oLocus, ptr, nn++, isGene) ;
		  } 
		break ;
	      case 1:  /* any other molec */
		if (n0 < 0 && isEssential && isMol && ptr)
		  { 
		    if (*isPutativep && doClean)
		      { vtxtClear (blkp) ; nn = 0 ; }
		    *isPutativep = FALSE ;
		    dictAdd (dict, ac_name(oGeneClass), 0) ;
		    n0 = ir ; gtOneGeneClassDescriptor (blkp, oGene, oLocus, ptr, nn++, isGene) ;
		  }
		break ;
	      case 2: /* true genetic */
		if (n2 < 0 && !isMol && ptr &&
		    strstr (ac_name(oGene), ac_name(oLocus)))
		  { 
		    if (*isPutativep && doClean)
		      { vtxtClear (blkp) ; nn = 0 ; }
		    *isPutativep = FALSE ;
		    dictAdd (dict, ac_name(oGeneClass), 0) ;
		    n2 = ir ; gtOneGeneClassDescriptor (blkp, oGene, oLocus, ptr, nn++, isGene) ;
		  } 
		break ;
	      case 3: /* other genetic */
		if (n2 != ir && !isMol && ptr)
		  { 
		    if (*isPutativep && doClean)
		      { vtxtClear (blkp) ; nn = 0 ; }
		    *isPutativep = FALSE ;
		    dictAdd (dict, ac_name(oGeneClass), 0) ;
		    gtOneGeneClassDescriptor (blkp, oGene, oLocus, ptr, nn++, isGene) ;
		  } 
		break ;
	      case 4: /* exact molec !essential */
		if (n4 < 0 && !isEssential && isMol && ptr &&
		    strstr (ac_name(oGene), ac_name(oLocus)))
		  { 
		    if (*isPutativep && doClean)
		      { vtxtClear (blkp) ; nn = 0 ; }
		    *isPutativep = FALSE ;
		    dictAdd (dict, ac_name(oGeneClass), 0) ;
		    n4 = ir ; gtOneGeneClassDescriptor (blkp, oGene, oLocus, ptr, nn++, isGene) ;
		  } 
		break ;
	      case 5: /* autres molec */
		if (n0 != ir  && n4 != ir && isMol && ptr)
		  { 
		    if (*isPutativep && doClean)
		      { vtxtClear (blkp) ; nn = 0 ; }
		    *isPutativep = FALSE ;
		    dictAdd (dict, ac_name(oGeneClass), 0) ;
		    n5 = ir ; gtOneGeneClassDescriptor (blkp, oGene, oLocus, ptr, nn++, isGene) ;
		  }
		break ;
	      }	      
	  }
      switch (jj)
	{
	case 1: /* true molec or wb id */
	  if (n0 < 0 && isEssential && molecDescription)
	    { 
	      n0 = -1 ;
	      n00 = 1 ;
	      if (*isPutativep)
		doClean = TRUE ;
	      *isPutativep = FALSE ;
	      if (1 && isGene) 
		vtxtPrintf (blkp, "%sencoding "
			    , vtxtPtr (blkp) ? ", " : ""
			    ) ; /* , isVoyelle (*molecDescription) ? "n" : "" */
	      else if (nn) 
		vtxtPrintf (blkp, "%s"
			    , vtxtPtr (blkp) ? ", " : " ") ;
	      gtOneGeneClassDescriptor (blkp, 0, 0, molecDescription, 0, isGene) ;
	      nn++ ;
	    }
	  break ;
	case 5: /* autres molec */ 
	  if (n00 < 0 && n4 < 0 && n5 < 0 && molecDescription)
	    {
	      n4 = -1 ; 
	      if (*isPutativep)
		doClean = TRUE ;
	      *isPutativep = FALSE ;
	      if (1 && isGene) 
		vtxtPrintf (blkp, "%sencoding "
			    , vtxtPtr (blkp) ? ", " : ""
			    ) ;
	      else if (nn) 
		vtxtPrintf (blkp, "%s"
			    , vtxtPtr (blkp) ? ", " : " ") ;
	      gtOneGeneClassDescriptor (blkp, 0, 0, molecDescription, 0, isGene) ;
	      nn++ ;
	    }
	  break ;
	}
    }
  ac_free (h) ;

  return vtxtPtr (blkp) ;
} /* gtGeneClassDescriptor */

/********************************************************************/

static BOOL gtBaseGeneDescriptor (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = handleCreate () ;
  BOOL ok = FALSE ;
  AC_TABLE tbl ;
  AC_OBJ oCog ;
  int nn = 0, ir ;
  const char *ptr ;
  
  if ((tbl = ac_tag_table (gmp->gene, "Locus_description", h)))
    for (ir = 0 ; ir < tbl->rows && tbl->cols >= 1; ir++)
      {
	ptr = ac_table_printable (tbl, ir, 0, 0) ;
	if (ptr && *ptr && !strstr(ptr,"ssential gene")) /* essential gene is reported separately */
	  {
	    ok = TRUE ;
	    vtxtPrintf (blkp, "%s%s", 
			nn++ > 0 ? "; " : "" , gtCleanUp(ptr)) ;
	  }
      }
  
  ok = gtGeneDescriptorTag (blkp, gmp) ? TRUE : FALSE ;

  if (!ok &&
      (ptr = ac_tag_text(gmp->gene, "WB_id", 0)))
    { 
      ok = TRUE ; 
      vtxtPrintf (blkp,  "%s%s", 
		  nn++ > 0 ? "; " : "" , gtCleanUp(ptr)) ;
    }
   
  if (!ok &&
      (oCog = ac_tag_obj (gmp->gene, "COG", h)) &&
      (ptr = ac_tag_text (oCog, "Title", 0))
      )
    { 
      ok = TRUE ; 
      vtxtPrintf (blkp,  "%s%s", 
		 nn++ > 0 ? "; " : "" , gtCleanUp(ptr)) ;
    }

  ac_free (h) ;
  return ok ;
} /* BaseGeneDescriptor */

/********************************************************************/

static char *gtBaseGeneFunction (vTXT blkp, GMP *gmp, AC_OBJ oProduct, BOOL doQuote)
{
  AC_HANDLE h = handleCreate () ; 
  AC_TABLE oIP, tbl, oCogInfo;
  AC_OBJ oCog, oPfam ;
  int nn=0, nn2, ir, jr ;
  const char *ptr, *ccp ;
  char *ptr2 ;
  char *qq1, *qq = doQuote ? "\"" : "" ;

  if ((tbl = ac_tag_table (gmp->gene, "Locus_description", h)))
    for (ir = 0 ; ir < tbl->rows && tbl->cols >= 1; ir++)
      {
	ptr = ac_table_printable(tbl, ir, 0, 0) ;
	if (ptr && *ptr && !strstr(ptr,"ssential gene")) /* essential gene is reported separately */
	  vtxtPrintf (blkp, "%s%s%s%s", 
		      nn++ > 0 ? ", " : "" , qq, gtSetUpper(gtCleanUp(ptr)),qq) ;
      }

  if (oProduct && /* better to get those in the product */
      (tbl = ac_tag_table (oProduct, "Pfam", h)))
    for (ir = 0 ; ir < tbl->rows; ir++)
      {
	oPfam = ac_table_obj (tbl, ir, 0, h) ;
	if ((oIP = ac_tag_table (oPfam, "Interpro", h)))
	  for (jr=0 ;jr < oIP->rows && oIP->cols >= 4;jr++) 
	    if ((ptr = ac_table_printable (oIP, jr, 1, 0)))
	      {
		vtxtPrintf (blkp, "%s%s %s  "
			    , nn++ > 0 ? ", " : "" , qq 
			    ,gtSetUpper(gtCleanUp(ptr))
			    ) ;
		ccp = ac_table_printable (oIP, jr, 3, "") ;
		vtxtPrintf (blkp, "%s [Pfam/Interpro/GO]%s"
			    ,ccp
			    ,qq) ;
	      }
      }
  
  if ((tbl = ac_tag_table (gmp->gene, "COG", h)))
    for (ir = 0 ; ir < tbl->rows && tbl->cols >= 1; ir++)
      {
	oCog = ac_table_obj (tbl, ir, 0, h) ;
	if ((oCogInfo = ac_tag_table (oCog, "Info", h)))
	  for (jr=0 ;jr < oCogInfo->rows && oCogInfo->cols >= 2;jr++)
	    if ((ptr = ac_table_printable (oCogInfo, jr, 1, 0)))
	      vtxtPrintf (blkp, "%s%s%s [COG annotation]%s", 
			  nn++ > 0 ? ", " : "" , qq, gtSetUpper(gtCleanUp(ptr)),qq) ;
      }
  
  if ((tbl = ac_tag_table (gmp->gene, "Locus_Phenotype", h)))
    for (ir = 0 ; ir < tbl->rows && tbl->cols >= 1; ir++)
      {
	if ((ptr = ac_table_printable(tbl, ir, 0, 0)))
	  {
	    vtxtPrintf (blkp, "%s %s",
			nn++ > 0 ? ", " : "" , qq) ;
	    vtxtPrintf (blkp, gtSetUpper(gtCleanUp(ptr))) ;
	    ptr2 = vtxtPtr (blkp) ; 
	    ptr2 += strlen(ptr2) - 1 ;
	    if (*ptr2 == '\n') *ptr2 = ' ' ;
	    vtxtPrint (blkp,  qq) ;
	  }
      }
  if ((tbl = ac_tag_table (gmp->gene, "Allele", h)))
    for (ir = 0 ; ir < tbl->rows && tbl->cols >= 1; ir++)
      {
	ptr = ac_table_printable (tbl, ir, 0, 0) ;
	if (ptr)
	  {
	    vtxtPrintf (blkp, "%s %sAllele %s%s", 
			nn++ > 0 ? ", " : "" , qq, gtSetUpper(gtCleanUp(ptr)),qq) ;
	  }
      }
  if ((tbl = ac_tag_table (gmp->gene, "RNAi", h)))
    for (ir = 0 ; ir < tbl->rows && tbl->cols >= 1; ir++)
      {
	oCog = ac_table_obj (tbl, ir, 0, h) ;
	qq1 = qq ; nn2 = 0 ;
	if (!strncmp("mv_", ac_name(oCog), 3))
	  continue ; 
	if ((oCogInfo = ac_tag_table (oCog, "Phenotype", h)) && ac_tag_text (oCog, "Phenotype", 0))
	  for (jr=0 ;jr < oCogInfo->rows && oCogInfo->cols >= 1; jr++)
	    if ((ptr = ac_table_printable (oCogInfo, jr, 0, 0)))
	      {
		nn2++ ;
		vtxtPrintf (blkp, "%s%s%s", 
			    nn++ > 0 ? ", " : "" , qq1, gtSetUpper(gtCleanUp(ptr))) ;
		if (ac_has_tag (oCog, "DIC_Movie_available"))
		  vtxtPrintf (blkp, "%s %s", 
			      ". Movies of embryo development in RNA interference experiments "
			      "are available at http://worm-srv1.mpi-cbg.de/dbScreen. ") ;
		qq1 = "" ; /* concatene al the pheno of the same rnai */
	      }
	if (nn2) vtxtPrint (blkp,  qq) ;
      }
  
  ac_free (h) ;
  return vtxtPtr (blkp) ;
} /* BaseGeneFunction */

/********************************************************************/

char *gtGeneName (vTXT blkp, GMP *gmp)
{
  vtxtClear (blkp) ;
  if (debug) vtxtPrint (blkp, "GeneName::") ;
  vtxtPrint (blkp,  ac_name (gmp->gene)) ;
  
  return vtxtPtr (blkp) ;
} /* gtGeneName */

/********************************************************************/

static char *gtSheddedGeneTitle (vTXT blkp, GMP *gmp)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_OBJ shedFromGene = gmp->tg ? ac_tag_obj (gmp->tg, "Shedded_from", h) : 0 ;

  if (shedFromGene)
    {
      vtxtPrint (blkp, ", contacting gene ") ;
      gmpObjLink (blkp, gmp, shedFromGene, ac_name (shedFromGene)) ;
    }
  ac_free (h) ; 

  return vtxtPtr (blkp) ;
} /* gtSheddedGeneTitle */

/********************************************************************/

char *gtGeneTitle (vTXT blkp, GMP *gmp, BOOL showSpecies)
{
  char *ptr1, *ptr2 = 0 ;
  BOOL isPutative = FALSE, isShedded = FALSE ;
  vTXT buf1 = vtxtCreate () ;
  vTXT buf2 = vtxtCreate () ;
  
  vtxtClear (blkp) ;
  
  if (debug) vtxtPrint (blkp, "GeneTitle::") ;
  
  if (showSpecies)
    vtxtPrint (blkp,  gmp->spci->speciesName) ;
  
  if (gmp->tg && ac_has_tag (gmp->tg, "Shedded_from"))
    isShedded = TRUE ;
  else
    ptr2 = gtProductBaseTitle (buf2, gmp, 0, 0, TRUE, &isPutative) ;
  if (!ptr2 || isPutative)
    {
      gtGeneBaseTitle (blkp, gmp, TRUE, isShedded) ;
    }
  else
    {
      ptr1 = gtGeneBaseTitle (buf1, gmp, FALSE, FALSE) ;
      if (ptr1 && strstr (ptr1, "encoding")) 
	/* title imported diretly from locuslink */
	vtxtPrint (blkp,  ptr1) ;      
      else 
	vtxtPrintf (blkp, "%s, %s", ptr1, ptr2) ;   
    }
 
  if (! strstr (vtxtPtr (blkp), "seudogene") &&
      (
       ac_has_tag (gmp->gene, "Pseudogene") ||
       (
	gmp->pg && 
	!strcasecmp (ac_tag_printable (gmp->pg, "Method", ""), "Pseudogene")
	)
       )
      )
    vtxtPrint (blkp, ", probable pseudogene") ;
  if (isShedded)
    gtSheddedGeneTitle (blkp, gmp) ;
  vtxtPrint (blkp, ".") ;
  gtSetUpper (vtxtPtr(blkp)) ;

  vtxtDestroy (buf1) ;
  vtxtDestroy (buf2) ;

  return vtxtPtr (blkp) ;
}  /* gtGeneTitle */

/********************************************************************/

char *gtGeneDescriptor (vTXT blkp, GMP *gmp)
{
  vtxtClear (blkp) ;
  if (debug) vtxtPrint (blkp, "GeneDescriptor::") ;
  gtBaseGeneDescriptor (blkp, gmp) ;
 
  gtSetUpper (vtxtPtr(blkp)) ;
  return vtxtPtr (blkp) ;
} /* gtGeneDescritor */

/********************************************************************/

char *gtGeneFunctionalName (vTXT blkp, GMP *gmp)
{
  AC_OBJ oProduct = gtGene2Product (gmp->gene, 0) ;
  char *cp = gtProductName (blkp, gmp, oProduct) ;
  ac_free (oProduct) ;
  return cp ;
} /* gtGeneFunctionalName */

/********************************************************************/

char *gtGeneFunctionalDescriptor (vTXT blkp, GMP *gmp, AC_OBJ oProduct, BOOL doQuote)
{
  AC_OBJ omm = 0, opp = 0 ;
  vtxtClear (blkp) ;
  
  omm = ac_tag_obj (oProduct,"Predicted_mRNA", 0) ;
  if (omm) opp = gtMrna2Product (omm, 0) ;
  if (opp) oProduct = opp ;

  gtBaseGeneFunction (blkp, gmp, oProduct, doQuote) ;
  gtSetUpper (vtxtPtr(blkp)) ;

  ac_free (omm) ;
  ac_free (opp) ;
  return vtxtPtr (blkp) ;
} /* gtGeneFunctionalDescriptor */

/********************************************************************/
/**********************************************************************/

static char *gtGeneIdTitle  (vTXT blkp, AC_OBJ oGene)
{
  AC_HANDLE h = handleCreate () ;
  AC_OBJ obj = 0 ;
  const char *ptr ;
  AC_ITER iter = 0 ;
  int n = 0 ;

  iter = ac_objquery_iter (oGene, "{>GeneId} SETELSE {>LocusLink} ; Title", h) ;
  
  while (ac_free (obj), iter && (obj = ac_iter_obj (iter)))
    {
      if ((ptr = ac_tag_text (obj, "Title", 0)))
	vtxtPrintf (blkp, "%s%s", n++ ? " and " : "", ptr) ;
    }
  ac_free (h) ;
  return vtxtPtr (blkp) ;
} /* gtGeneidTitle */

/********************************************************************/
/* 19 may 2003,
 * we decide to disregard gene->product_title
 * and to transfer the kantor->hand.title into
 * using the concatenation of all gene->descriptor
 *
 */
static char *gtGeneDescriptorTag (vTXT blkp, GMP *gmp)
{
  AC_TABLE tbl = ac_tag_table (gmp->gene, "Descriptor", 0) ;
  int nn = 0, ir ;
  const char *ptr = 0 ;

  if (gmp->style == 's' &&
      (tbl = ac_tag_table (gmp->gene, "Genbank_product_descriptor", 0)))
    for (ir = 0 ; ir < tbl->rows && tbl->cols >= 0; ir++)
      if ((ptr = ac_table_printable (tbl, ir, 0, 0)) && *ptr)
	{
	  if (strstr (ptr, "Summary:"))
	    continue ;
	  vtxtPrintf (blkp, "%s%s", 
		      nn++ > 0 ? "; " : "" , gtCleanUp(ptr)) ;
	}

  if (! nn &&
      (tbl = ac_tag_table (gmp->gene, "Descriptor", 0)))
    for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
      if ((ptr = ac_table_printable (tbl, ir, 0, 0)) && *ptr)
	{
	  if (strstr (ptr, "Summary:"))
	    continue ;
	  vtxtPrintf (blkp, "%s%s", 
		      nn++ > 0 ? "; " : "" , gtCleanUp(ptr)) ;
	}
  
  ac_free (tbl) ;
  return vtxtPtr (blkp) ;
} /* gtProductDescriptorTitle */

/********************************************************************/
/* kantor tile (short or long or replaced by molecular gene class descriptor */
static char *gtProductBaseTitle (vTXT blkp, GMP *gmp, AC_OBJ oProduct, AC_OBJ oGF, BOOL isGene, BOOL *isPutativep)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL ok = FALSE ;
  int nn = 0 ;
  const char *ptr, *skt = 0 ;
  char *gclTitle = 0 ;
  char *descriptor = 0 ;
  vTXT buf, buf1 ;

  if (oProduct)
    {
      skt = ac_tag_text (oProduct, "Short_kantor_title", 0) ;
    }
  else if (oGF)
    {  
      AC_OBJ oPredictedMrna =  ac_tag_obj (oGF, "Predicted_mRNA", h) ;
      AC_OBJ oPredictedProduct =  oPredictedMrna ? ac_tag_obj (oPredictedMrna, "Product", h) : 0 ;

      skt = oPredictedProduct ? strnew (ac_tag_printable (oPredictedProduct, "Short_kantor_title", 0), h) : 0 ; 
    }
  else if (gmp->gene)
    {  /* select the best product, i.e. the first one, preferably real */
      oProduct = gmp->product ;
      if (isPutativep && oProduct &&
	  (ptr = ac_tag_printable (oProduct, "Blastp_title",0)) &&
	  strstr (ptr, "rotein of unknown function"))
	*isPutativep = TRUE ;
      else
	skt = oProduct ? strnew (ac_tag_printable (oProduct, "Short_kantor_title", 0), h) : 0 ; 
    }
  else
    { ac_free (h) ; return NULL ; }

  if (0 && skt && !strcasecmp (skt, "putative protein"))
    skt = 0 ;
  buf = vtxtHandleCreate (h) ;
  buf1 = vtxtHandleCreate (h) ;

  if (gmp->gene)
    {
      descriptor = gtGeneDescriptorTag (buf1, gmp) ;
      
      if (! descriptor && ! ac_has_tag (gmp->gene, "Shedded_from") &&
	  gmp->Spc != WORM)
	descriptor = gtGeneIdTitle (buf1, gmp->gene) ;

      gclTitle = gtGeneClassDescriptor (buf, gmp->gene, skt, descriptor, isGene, isPutativep) ;
    }

  if (gclTitle)
    {
      ok = TRUE ;      
      if (nn++) vtxtPrint (blkp, ", ") ;
      vtxtPrint (blkp,  gtCleanUp(gclTitle)) ;
    } 
  
  if (!ok)
    {
      if (oGF)
	vtxtPrint (blkp,  ac_name(oGF)) ;
      else if (oProduct)
	vtxtPrint (blkp,  ac_name(oProduct)) ;
      ok = TRUE ;
    }

  if (!isGene &&
      (oProduct || (oGF && ac_has_tag (oGF, "Has_exact_est"))))
    gtMolecularWeight (blkp, oProduct, oGF, TRUE) ;

  if ((ptr = vtxtPtr (blkp))  &&
      !strstr (ptr,  "precursor") &&
      !strstr (ptr,  "secreted") &&
      gtIsPrecursor (oProduct, oGF))
    vtxtPrint (blkp, " precurssor") ;

  if (!isGene &&
      (ptr = gtGene2Locus (gmp->gene)))
    vtxtPrintf (blkp, " (%s)", ptr) ;

  ac_free (h) ;

  return ok ? vtxtPtr (blkp) : 0 ;
} /* gtProductBaseTitle */

/********/

static const char *gtGene2Locus (AC_OBJ oGene)
{
  AC_TABLE tbl ;
  int ir ;

  if ((tbl = ac_tag_table (oGene, "Locus", 0)))
    for (ir = 0 ; ir < tbl->rows && tbl->cols >= 1; ir++)
      {
	if (strstr (ac_name(oGene), ac_table_printable (tbl, ir, 0, "")))
	  { ac_free (tbl) ; return ac_name (oGene) ; }
      }
  ac_free (tbl) ;
  return ac_name (oGene) ;
}

/*********/

AC_OBJ gtPredictedMrna2Product (AC_OBJ oGF, AC_HANDLE h)
{
  AC_OBJ oProd = 0, oMrna = ac_tag_obj (oGF,"Predicted_mrna", 0) ;
  if (oMrna)
    {
      oProd = gtMrna2Product (oMrna, h) ;
      ac_free (oMrna) ;
    }
  return oProd ;
}

/*********/

AC_OBJ gtMrna2Product (AC_OBJ oMrna, AC_HANDLE h)
{
  AC_TABLE tbl ;
  AC_OBJ oProduct = 0, oBestProd = 0 ;
  int n, ir, nmax ;
  
  if (!oMrna)
    return 0 ;

  tbl = ac_tag_table (oMrna, "Product", 0) ;
  for (ir=nmax=0 ; tbl && ir < tbl->rows ; ir++)
    {
      oProduct = ac_table_obj (tbl,ir,0, h) ;
      if (ac_has_tag (oProduct, "Best_product"))
	{ oBestProd = oProduct ; break ; }
     else
       ac_free (oProduct) ;
    }

  if (! oBestProd)
    {
      for (ir=nmax=0 ; tbl && ir < tbl->rows ; ir++)
	{
	  n = ac_table_int (tbl,ir,2, 0) - ac_table_int (tbl,ir,1, 0) ;
	  if (!oBestProd || n > nmax)
	    { nmax = n ; ac_free (oBestProd) ; oBestProd = ac_table_obj (tbl,ir,0, h) ;}
	}
    }
  ac_free (tbl) ;
  return oBestProd ;
} /* gtMrna2Product */

/*********/

AC_OBJ gtGene2Product (AC_OBJ oGene, AC_HANDLE h)
{
  AC_TABLE gProd ;
  AC_OBJ oProd, oBestProd = 0 ;
  int n, ir, nmax = 1 ;

  gProd = ac_tag_table (oGene, "Product", 0) ;
  if (!gProd)
    return NULL;
  for (ir=nmax=0 ; ir < gProd->rows ; ir++)
    {
      oProd = ac_table_obj (gProd, ir, 0, h) ;
      n = 0 ;
      if (ac_has_tag (oProd, "Best_product"))
	n = ac_tag_int (oProd, "Open_length", 0) ;
      if (n > nmax)
	{ ac_free (oBestProd) ; nmax = n ; oBestProd = oProd ;}
      else
	ac_free (oProd) ;
    }
  ac_free (gProd) ;
  return oBestProd ;
} /* gtMrna2Product */

/*********/

static void gtMrnaCompleteness (vTXT blkp, AC_OBJ oMrna, BOOL showDetails)
{
  AC_OBJ oProduct = gtMrna2Product(oMrna, 0) ;
  int isComplete=0 ;
  int isCompleteCDS=0 ;
  
  if (ac_has_tag(oMrna,"Complete") || (ac_has_tag(oMrna,"Found5p") && ac_has_tag(oMrna,"Found3p") ))
    isComplete=1;
  if(ac_has_tag(oProduct,"NH2_Complete") && ac_has_tag(oProduct,"COOH_Complete"))
    isCompleteCDS=1;
  
  if(isComplete && isCompleteCDS)
    vtxtPrint (blkp,"complete mRNA");
  /* kim wants to discard the correct but to deatilled next paragraph */
  else if (!showDetails)
    vtxtPrint (blkp,"mRNA");
  else
    {
      /* this is all correct, just banished from refseq. we may use it for genbank */
      if(isCompleteCDS)
	vtxtPrint (blkp,"complete CDS mRNA");
      
      else if (ac_has_tag(oProduct,"NH2_Complete") && !ac_has_tag(oProduct,"COOH_Complete") )
	vtxtPrint (blkp,"partial mRNA, 3' incomplete");
      
      else if (!ac_has_tag(oProduct,"NH2_Complete") && ac_has_tag(oProduct,"COOH_Complete") )
	if(FALSE)vtxtPrint (blkp,"possibly 5' incomplete mRNA");
	else vtxtPrint (blkp,"mRNA");
      
      else 
	vtxtPrint (blkp,"partial mRNA, internal fragment");
    }
  ac_free (oProduct) ;

} /* gtMrnaCompleteness */

/*********/

static char *gtMrnaBaseTitle (vTXT blkp, GMP *gmp, BOOL showCompleteness)
{
  BOOL isPutative = FALSE ;
  AC_OBJ oProd = 0 ;
  char *ptr = 0 ;
  vTXT buf = vtxtCreate () ;

  oProd = gtMrna2Product(gmp->mrna, 0) ;
  if (gtProductBaseTitle (buf, gmp, oProd, 0, 0, &isPutative))
    ptr = vtxtPtr (buf) ;

  if (gmp->gene &&
      ptr &&
      !strstr (ptr, "ssential") && /* do not repeat essential */
      gtIsEssentialGene (gmp->gene))
    vtxtPrint (blkp, "essential ") ;

  if (gtIsCloudGene (gmp->gene))
    vtxtPrint (blkp, "cloud ") ;
  if (ptr)
    vtxtPrintf (blkp, " %s ", ptr) ;

  if (gmp->variant) 
    vtxtPrintf (blkp, "alternative variant %s, ", gmp->variant) ;
  
  gtMrnaCompleteness (blkp, gmp->mrna, showCompleteness) ;
  
  ac_free (oProd) ;
  vtxtDestroy (buf) ;
  return vtxtPtr (blkp) ;
} /* gtMrnaBaseTitle */

/********************************************************************/

static void gtMrnaMarkUpBaseTitle (vTXT blkp, GMP *gmp)
{
  char *ptr, *cp ;
  AC_OBJ oLocusLink, lid ;
  AC_HANDLE h = handleCreate () ;
  vTXT buf = vtxtHandleCreate (h) ;
  
  ptr = gtMrnaBaseTitle (buf, gmp, FALSE) ;

  if (ac_has_tag (gmp->gene, "Shedded_from"))
    vtxtPrint (blkp,   ptr) ; 
  else if ((lid = ac_tag_obj (gmp->gene, "GeneId", h)))
    {
      char linkBuf[4096] ;

      oLocusLink = ac_tag_obj (gmp->gene, "LocusLink", h) ;
      if (!oLocusLink) oLocusLink = gmp->gene ;
      sprintf (linkBuf, ENTREZ_LINK, ac_name(lid))  ; 
      
      if ((cp = strstr (ptr, ac_name(oLocusLink))))
	{
	  *cp = 0 ;
	  vtxtPrint (blkp,  ptr) ;
	  gmpURL (blkp, gmp, linkBuf, ac_name(oLocusLink)) ;
	  vtxtPrint (blkp,   cp + strlen (ac_name(oLocusLink))) ;
	}
      else
	{
	  gmpURL (blkp, gmp, linkBuf, ac_name(oLocusLink)) ;
	  vtxtPrintf (blkp, ": %s",  ptr) ; 
	}
    }
  else if ((lid = ac_tag_obj (gmp->gene, "LocusId", h)))
    {
      char linkBuf[4096] ;

      oLocusLink = ac_tag_obj (gmp->gene, "LocusLink", h) ;
      if (!oLocusLink) oLocusLink = gmp->gene ;
      sprintf (linkBuf, ENTREZ_LINK, ac_name(lid))  ; 
      
      if ((cp = strstr (ptr, ac_name(oLocusLink))))
	{
	  *cp = 0 ;
	  vtxtPrint (blkp,  ptr) ;
	  gmpURL (blkp, gmp, linkBuf, ac_name(oLocusLink)) ;
	  vtxtPrint (blkp,   cp + strlen (ac_name(oLocusLink))) ;
	}
      else
	{
	  gmpURL (blkp, gmp, linkBuf, ac_name(oLocusLink)) ;
	  vtxtPrintf (blkp, ": %s",  ptr) ; 
	}
    }
  else
    vtxtPrint (blkp,   ptr) ; 
  ac_free (h) ;

   return ;
} /* gtMrnaMarkUpBaseTitle */

/********************************************************************/

char *gtMrnaName (vTXT blkp, GMP *gmp)
{
  vtxtClear (blkp) ;
  if (debug) vtxtPrint (blkp, "MrnaName::") ;
  gtMrnaBaseTitle (blkp, gmp, 0) ;
 
  gtSetUpper (vtxtPtr(blkp)) ; 
  vtxtPrint (blkp, ".") ;
  return vtxtPtr (blkp) ;
} /* gtMrnaName */

/********************************************************************/

char *gtMrnaSuffix (const char *geneName, const char *mrnaName, AC_HANDLE h) 
{
  char *cp1, *cp = strnew (mrnaName, h) ;

  if (strstr (cp, geneName) == cp) 
    cp += strlen (geneName) ;
  cp1 = strstr(cp, "-unspliced") ;
  if (cp1) *cp1 = 0 ;
  if (strlen (cp) > 5)
    { /* cut the date suffix */
      char *cp2 = cp + strlen (cp) - 5 ;
      *cp2 = 0 ;
    }
  if (cp1)
    { 
      cp1 = cp + strlen(cp) ;
      *cp1++ = '-' ; *cp1++ = 'u' ; *cp1 = 0 ; 
    }
  if (!*cp)
    {
      cp = strnew (mrnaName, h) ;
      if (strlen (cp) > 5)
	{ /* cut the date suffix */
	  char *cp2 = cp + strlen (cp) - 5 ;
	  *cp2 = 0 ;
	}
    }
  return cp ;
} /* gtMrnaName */

/********************************************************************/

char *gtMrnaTitle (vTXT blkp, GMP *gmp)
{
  vtxtClear (blkp) ;
  if (debug) vtxtPrint (blkp, "MrnaTitle::") ;
  if (1) vtxtPrintf (blkp, "%s ", gmp->spci-> speciesName) ;

  if (gmp->markup)
    gtMrnaMarkUpBaseTitle (blkp, gmp) ;
  else
    gtMrnaBaseTitle (blkp, gmp, gmp->markup) ;
 
  gtSetUpper (vtxtPtr(blkp)) ;
  return vtxtPtr (blkp) ;
} /* gtMrnaTitle */

/********************************************************************/
/********************************************************************/

char *gtProductName (vTXT blkp, GMP *gmp, AC_OBJ oProduct)
{ 
  BOOL isPutative = FALSE ;

  vtxtClear (blkp) ;
  if (debug) vtxtPrint (blkp, "ProductName::") ;
  gtProductBaseTitle (blkp, gmp, oProduct, 0, 0, &isPutative) ;
  /*
    cp = vtxtPtr (blkp) ;
  if (cp) while (*++cp) if (*cp == '(') *(cp-1)=0 ; 
  cut the mW and pI */
  if (0) vtxtPrint (blkp, ".") ; /* causes kans to add a [c.elegans]. after the existing . */

  gtSetUpper (vtxtPtr(blkp)) ;
  return vtxtPtr (blkp) ;
} /* gtProductName */

/********************************************************************/

char *gtProductMolecularName (vTXT blkp, GMP *gmp, AC_OBJ oProduct)
{ 
  char *ptr ;
  const char *skt ;
  vtxtClear (blkp) ;

  if (!gtGeneDescriptorTag (blkp, gmp) &&
      oProduct &&
      (skt =  ac_tag_text (oProduct, "Short_kantor_title", 0)))
    vtxtPrint (blkp,  skt) ;

  ptr = vtxtPtr (blkp) ;
  if (ptr)
    gtCleanUp (ptr) ;
  return vtxtPtr (blkp) ;
} /* gtProductName */

/********************************************************************/

char *gtProductTitle (vTXT blkp, GMP *gmp)
{
  BOOL isPutative = FALSE ;

  vtxtClear (blkp) ;
  if (debug) vtxtPrint (blkp, "ProductTitle::") ;
  vtxtPrintf (blkp, "%s, ", ac_name(gmp->product)) ;
  gtProductBaseTitle (blkp, gmp, gmp->product, 0, 0, &isPutative) ;
  vtxtPrintf (blkp, " [%s]", gmp->spci-> speciesName) ; 
  vtxtPrint (blkp, ".") ;

  gtSetUpper (vtxtPtr(blkp)) ;
  return vtxtPtr (blkp) ;
} /* gtProductTitle */

/********************************************************************/
/********************************************************************/
/* abs(*nclp) is set to the # of clones, and is negattive is mainly yk */
int gtPredictedMrnaSupport (AC_OBJ oGF, int *nclp)
{
  AC_HANDLE h = handleCreate () ;
  int n = 0, ir, ncl, nyk ;
  AC_TABLE gCl ;
  AC_OBJ oGene ;
  oGene = ac_tag_obj (oGF, "Model_of_gene", h) ;
  
  if (nclp)
    {
      *nclp = 0 ;
      if ((gCl = ac_tag_table (oGene, "Has_cDNA_clone", h)))
	{
	  for (nyk = ncl = ir = 0 ; ir < gCl->rows ; ir++)
	    {
	      ncl++ ;
	      if (*(ac_table_printable (gCl, ir, 0, "")) == 'y')
		nyk++ ;
	    }
	  if (2*nyk > ncl) ncl = -ncl ;
	  *nclp = ncl ;
	}
    }
  if (ac_has_tag (oGF, "has_exact_est"))
    n = 4 ;
  else if (gtHasEst(oGene, oGF))
    n = 3 ;
  else if (gtHasOst(oGene))
    n = 2 ;
  else if (gtHasExpression(oGene))
    n = 1 ;
  else
    n = 0 ;

  ac_free (h) ;
  return n ;
} /* gtMrnaSupport */

/*********/
#ifdef JUNK
may 19, 2003 seems useless
static void gtMrnaShowSupport (vTXT blkp, AC_OBJ oGF)
{
  AC_OBJ oGene ;
  oGene = ac_tag_obj (oGF, "Model_of_gene") ;

  if (ac_has_tag (oGF, "has_exact_est"))
    vtxtPrint (blkp, ", fully supported by cDNA") ;
  else if (gtHasEst(oGene, oGF))
    vtxtPrint (blkp, ", with cDNA support") ;
  else if (gtHasOst(oGene))
    vtxtPrint (blkp, ", with RT-PCR support") ;
  else if (gtHasExpression(oGene))
    vtxtPrint (blkp, ", supported by expression data") ;
  else
    vtxtPrint (blkp, ", not yet supported by cDNA or RT-PCR") ;

  return ;
} /* gtMrnaShowSupport */
#endif
/*********/

static char *gtPredictedMrnaBaseTitle (vTXT blkp, GMP *gmp, AC_OBJ oGF, BOOL showDetails)
{
  AC_HANDLE h = handleCreate () ;
  BOOL isPutative = FALSE ;
  AC_OBJ oNN, oGene = ac_tag_obj (oGF, "Model_of_gene", h) ;
  const char *gName = oGene ? ac_name (oGene) : 0 ;
  char *ptr = 0 ;
  char *gfName = gtGfName (oGF) ;
  vTXT buf = vtxtCreate () ;

  ptr = gtProductBaseTitle (buf, gmp, 0, oGF, 0, &isPutative) ;
  if (oGene &&
      ptr &&
      !strstr (ptr, "ssential") &&
      gtIsEssentialGene (oGene))
    vtxtPrint (blkp, "essential ") ;

  if (gtIsCloudGene (oGene))
    vtxtPrint (blkp, "cloud ") ;

  if (ptr)
    vtxtPrint (blkp,  ptr) ;
  else
    {
      vtxtPrint (blkp,  gfName) ;
      if (!strcmp (ac_name(oGF), gName) &&
	  (oNN = ac_tag_obj (oGene, "NewName", h)))
	gName = ac_name (oNN) ;
    }
  
  if (gmp && gmp->variant)
    vtxtPrintf (blkp, ", alternative variant %s", gmp->variant) ;
  
  if (!showDetails)
    {
      if (!gtHasLocus (oGene) && !gtHasEst(oGene, oGF) &&  !gtHasOst (oGene) && !gtHasExpression (oGene))
	vtxtPrint (blkp, ", predicted mRNA") ;
      else if (ac_has_tag (oGF, "has_exact_est"))
	vtxtPrint (blkp, ", mRNA") ;
      else
	vtxtPrint (blkp, ", mRNA") ;
    }
  else
    {
      if (!gtHasLocus (oGene) && !gtHasEst(oGene, oGF) &&  !gtHasOst (oGene) && !gtHasExpression (oGene))
	vtxtPrint (blkp, ", predicted CDS") ;
      else if (ac_has_tag (oGF, "has_exact_est"))
	vtxtPrint (blkp, ", complete CDS") ;
      else
	vtxtPrint (blkp, ", CDS") ;
    }
  
  ac_free (h) ;
  vtxtDestroy (buf) ;
  return vtxtPtr (blkp) ;
}  /* gtPredictedMrnaBaseTitle */

/********************************************************************/

char *gtPredictedMrnaName (vTXT blkp, GMP *gmp, AC_OBJ oGF)
{
  vtxtClear (blkp) ;
  if (debug) vtxtPrint (blkp, "PredictedMrnaName::") ;
  gtPredictedMrnaBaseTitle (blkp, gmp, oGF, 0) ;

  vtxtPrint (blkp, ".") ;
  gtSetUpper (vtxtPtr(blkp)) ;
  return vtxtPtr (blkp) ;
} /* gtPredictedMrnaName */

/********************************************************************/

char *gtPredictedMrnaTitle (vTXT blkp, GMP *gmp)
{
  vtxtClear (blkp) ;
  if (debug) vtxtPrint (blkp, "PredictedMrnaTitle::") ;
  if (1) vtxtPrintf (blkp, "%s ", gmp->spci->speciesName) ; 

  gtPredictedMrnaBaseTitle (blkp, gmp, gmp->pg, 0) ;
 
  vtxtPrint (blkp, ".") ;
  gtSetUpper (vtxtPtr(blkp)) ;
  return vtxtPtr (blkp) ;
} /* gtMrnaTitle */

/********************************************************************/
/********************************************************************/
/* name, activity ->function */
char *gtPredictedProductName (vTXT blkp, GMP *gmp, AC_OBJ oGF, BOOL insist)
{
  BOOL isPutative = FALSE ;

  vtxtClear (blkp) ;
  if (debug) vtxtPrint (blkp, "PredictedProductName::") ;
  
  if (insist && !gtPredictedMrnaSupport (oGF, 0))
    vtxtPrint (blkp, "predicted CDS, ") ;
  gtProductBaseTitle (blkp, gmp, 0, oGF, 0, &isPutative) ;

  if (0) vtxtPrint (blkp, ".") ; /* causes kans to add a [c.elegans]. after the existing . */

  gtSetUpper (vtxtPtr(blkp)) ;
  return vtxtPtr (blkp) ;
} /* gtPredictedProductName */

/********************************************************************/

char *gtPredictedProductTitle (vTXT blkp, GMP *gmp)
{ 
  BOOL isPutative = FALSE ;

  vtxtClear (blkp) ;
  if (debug) vtxtPrint (blkp, "PredictedProductTitle::") ;
  if (!gtPredictedMrnaSupport (gmp->pg, 0))
    vtxtPrint (blkp, "predicted CDS") ;
  vtxtPrintf (blkp, "%s (%s), ", ac_name(gmp->pg), gtGene2Locus (gmp->gene)) ;
  gtProductBaseTitle (blkp, gmp, 0, gmp->pg, 0, &isPutative) ;

  vtxtPrintf (blkp, " [%s]", gmp->spci-> speciesName) ; 
  vtxtPrint (blkp, ".") ;
  
  gtSetUpper (vtxtPtr(blkp)) ;
  return vtxtPtr (blkp) ;
} /* gtPredictedProductTitle */

/********************************************************************/

char *gtGfName (AC_OBJ oGF) 
{
  AC_HANDLE h = handleCreate () ; 
  AC_OBJ oGene, oNM = 0 ;
  AC_TABLE gTag ;
  const char *variant = 0 ;
  static char *old = 0 ;
  static vTXT buf = 0 ;

  if (old && !strcmp(ac_name(oGF), old))
    return vtxtPtr (buf) ;

  messfree (old) ; vtxtDestroy (buf) ;
  oGene = ac_tag_obj (oGF, "Model_of_gene", h) ;
  gTag = ac_tag_table (oGene, "Genefinder", h) ;
  buf = vtxtCreate () ;

  if (gTag && gTag->rows > 1)
    {
      oNM = ac_tag_obj (oGF, "NM_id", h) ;
      variant = oNM ? ac_tag_printable (oNM, "Variant", 0) : 0 ;
    }

  if (variant)
    vtxtPrintf (buf, "%s.%s", ac_name(oGene), variant) ;
  else
    vtxtPrint (buf,  ac_name(oGF)) ;
  
  old = strnew (ac_name (oGF), 0) ;

  ac_free (h) ;
  return vtxtPtr (buf) ;
}  /* gtGfName */

/********************************************************************/

char *gtReportTitles (AC_DB db, AC_OBJ oProduct)
{
  AC_OBJ oMrna = ac_tag_obj (oProduct, "Mrna", 0), oGF = 0 ;
  char *ptr ;
  GMP *gmp = gmpCreate (db, 0, 0, 0, 0, oProduct, 's', 'm') ;
  vTXT buf = vtxtCreate () ;
  vTXT blk = vtxtCreate () ;

  if (ac_has_tag (oMrna, "From_gene"))
    {
      if ((ptr = gtMrnaName (blk, gmp)))
	vtxtPrintf (buf, "Mrna %s\nTitle \"%s\"\n\n", ac_name(oMrna), ptr) ;
      if ((ptr = gtProductName (blk, gmp, oProduct)))
	vtxtPrintf (buf, "Product %s\nTitle \"%s\"\n\n", ac_name(oProduct), ptr) ;
    }
  else if ((oGF = ac_tag_obj (oMrna, "From_prediction", 0)))
    {
      if ((ptr = gtPredictedMrnaName (blk, gmp, oGF)))
	vtxtPrintf (buf, "Mrna %s\nTitle \"%s\"\n\n", ac_name(oMrna), ptr) ;
      if ((ptr = gtPredictedProductName (blk, gmp, oGF, TRUE)))
	vtxtPrintf (buf, "Product %s\nTitle \"%s\"\n\n", ac_name(oProduct), ptr) ; 
    }

  ptr = strnew (vtxtPtr (buf), 0) ;
  ac_free (oMrna) ;
  ac_free (oGF) ;
  gmpDestroy (gmp) ;
  vtxtDestroy (blk) ;
  vtxtDestroy (buf) ;

  return ptr ;
} /* gtReportTitles */

/********************************************************************/
/********************************************************************/

char *gtPredictedTrnaTitle (vTXT blkp, GMP *gmp, BOOL showSpecies, char **typep, AC_HANDLE h)
{
  AC_OBJ oGF = gmp->pg ;
  const char *ptr, *type ;

  vtxtClear (blkp) ;
  if (debug) vtxtPrint (blkp, "PredictedTrnaTitle::") ;
  if (showSpecies) vtxtPrintf (blkp, "%s ", gmp->spci-> speciesName) ;
 
  ptr = ac_tag_printable (gmp->gene, "Descriptor", ac_name(oGF)) ;

  vtxtPrint (blkp,  gtSetUpper (gtCleanUp(ptr))) ;
  if ((ptr = gtGene2Locus(gmp->gene)))
    vtxtPrintf (blkp, " (%s)", ptr) ;

  if (ac_has_tag (oGF, "Pseudogene"))
    type = "pseudogene" ;
  else if (ac_has_tag (oGF, "tRNA"))
    type = "tRNA" ;
  else
     {
       type = ac_tag_printable (oGF, "Transcript", 0) ;
       if (!type)
	 type = ac_tag_printable (oGF, "Method", "RNA") ;
     }
  if (typep) *typep = strnew (type, h) ;

  vtxtPrintf (blkp, ", %s", type) ;

  vtxtPrint (blkp, ".") ;
  gtSetUpper (vtxtPtr(blkp)) ;

  return vtxtPtr (blkp) ;
} /* gtMrnaTitle */

/********************************************************************/
/********************************************************************/
/********************** gmp constructors ***************************************/
/* gene:mrna:protein, style = s|r|x, view = g:gene|m:mrna|b:blast */

GMP *gmpCreate (AC_DB db, AC_OBJ oGene, AC_OBJ oTg, AC_OBJ oMrna, AC_OBJ oGF, AC_OBJ oProduct, char style, char view)
{
  AC_HANDLE h = handleCreate () ;
  GMP *gmp = (GMP *) messalloc (sizeof (GMP)) ;

  gmp->h = h ;

  gmp->db = db ;
  gmp->style = style ;
  gmp->pictureType = 'B' ; /* gif by default, tthe other possibility is 'S' for flash */
  gmp->view = view ;
  gmp->Spc = ficheDBCreatureID (db) ;
  gmp->spci = SpcI + gmp->Spc ;
  if (style == 'x') gmp->markup = TRUE ;
  if (view == 'z') return gmp ;
  if (oMrna)
    {
      gmp->mrna = oMrna ;

      if ((gmp->tg = ac_tag_obj (oMrna, "From_gene", h))) 
	gmp->gene = ac_tag_obj (gmp->tg, "Gene", h) ;
      else if ((gmp->pg = ac_tag_obj (oMrna, "From_prediction", h)))
	gmp->gene = ac_tag_obj (gmp->pg, "Model_of_gene", h)  ;

      gmp->product = gtMrna2Product (gmp->mrna, h) ;
    }
  else if (oTg)
    {
      gmp->tg = oTg ;

      gmp->mrna =  ac_tag_obj (oTg, "mRNA", h) ;
      gmp->gene = ac_tag_obj (gmp->tg, "Gene", h) ;
      if (gmp->mrna)
	gmp->product = gtMrna2Product (gmp->mrna, h) ;
    }
  else if (oGene)
    {
      AC_TABLE gTg ;
      int ir ;

      gmp->gene = oGene ;
      gTg = ac_tag_table (oGene, "Transcribed_gene", h) ;
      for (ir = 0 ; gTg && ir < gTg->rows ; ir++)
	{
	  gmp->tg = ac_table_obj (gTg, ir, 0, h) ;
	  if (gTg->rows > 1 && ac_has_tag (gmp->tg, "Shedded_from"))
	    ac_free (gmp->tg) ;
	  else
	    break ;
	}
      for (ir = 0 ; gTg && !gmp->tg && ir < gTg->rows ; ir++)
	{
	  gmp->tg = ac_table_obj (gTg, ir, 0, h) ;
	  if (ac_has_tag (gmp->tg, "gt_ag"))
	    break ;
	}
      for (ir = 0 ; gTg && !gmp->tg && ir < gTg->rows ; ir++)
	{
	  gmp->tg = ac_table_obj (gTg, ir, 0, h) ;
	  break ;
	}
      ac_free (gTg) ;
      if (!gmp->tg)
	gmp->pg = ac_tag_obj (oGene, "Genefinder", h);
      gmp->product = gmp->gene ? gtGene2Product (gmp->gene, h) : 0 ;
      gmp->mrna = gmp->product ? ac_tag_obj (gmp->product, "Mrna", h) : 0 ;
    }

  else if (oGF)
    {
      gmp->pg = oGF ;
      gmp->gene = ac_tag_obj (gmp->pg, "Model_of_gene", h) ;
      gmp->mrna = ac_tag_obj (oGF, "Predicted_mrna", h) ; 
      gmp->product = gtMrna2Product (gmp->mrna, h) ;
    }

  else if (oProduct)
    {
      gmp->product = oProduct ;
      gmp->mrna = ac_tag_obj (gmp->product, "mRNA", h) ;
      gmp->tg = gmp->mrna ? ac_tag_obj (gmp->mrna, "From_gene", h) : 0 ;
      gmp->gene = gmp->tg ? ac_tag_obj (gmp->tg, "Gene", h) : 0 ;
    }

  gmp->useAm = gmp->mrna ? ac_has_tag (gmp->mrna, "From_AM") : FALSE ;
  gmp->kantor = gmp->product ? ac_tag_obj (gmp->product, "Kantor", h) : 0 ;
  if (gmp->mrna && gmp->gene)
    {
      int nn, nMrna = 0 ;
      AC_KEYSET ks = ac_objquery_keyset (gmp->gene, ">transcribed_gene ; >mrna", 0) ;

      nMrna = ks ? ac_keyset_count (ks) : 0 ;
      ac_free (ks) ;
      
      if (nMrna > 1 &&
	  (nn = strlen (ac_name(gmp->gene))) &&
	  !strncmp (ac_name(gmp->gene), ac_name(gmp->mrna), nn) &&
	  *(ac_name(gmp->mrna) + nn) == '.')
	gmp->variant = strnew (ac_name (gmp->mrna) + nn + 1, h) ;
    }
  else if (gmp->mrna && gmp->tg)
    {
      int nn, nMrna = 0 ;
      AC_KEYSET ks = ac_objquery_keyset (gmp->tg, ">mrna", 0) ;

      nMrna = ks ? ac_keyset_count (ks) : 0 ;
      ac_free (ks) ;
      if (nMrna > 1 &&
	  (nn = strlen (ac_name(gmp->gene))) &&
	  !strncmp (ac_name(gmp->gene), ac_name(gmp->mrna), nn) &&
	  *(ac_name(gmp->mrna) + nn) == '.')
	gmp->variant = strnew (ac_name (gmp->mrna) + nn + 1, h) ;
    }

  if (gmp->gene)
    {
      AC_KEYSET ks = ac_objquery_keyset(gmp->gene,">transcribed_gene;>mrna", 0) ;  
      gmp->nMrna = ac_keyset_count (ks) ;
      ac_free (ks) ;
    }

  if (!oProduct && !oTg && !gmp->gene)
    {
      ac_free (gmp->h) ;
      messfree (gmp) ; /* resets gmp=0 */
    }
  return gmp ;
}

void gmpDoDestroy (GMP *gmp)
{
  if (gmp)
    {
      ac_free (gmp->h) ;
      messfree (gmp) ;
    }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

BOOL gmpHelp (vTXT blkp, GMP *gmp, char *file, const char *title)
{
  if (gmp->markup)
    {
      vtxtPrintf (blkp
		  , "<a href='%s.html#%s'><span class='explain'> ?</span></a>\n"
		  , file, title) ;
      return TRUE ;
    }
  return FALSE ;
} /* gmpHelp */

/******************************************************************************/
/******************************************************************************/
/* declare a yellow line */
BOOL gmpBlockInit (vTXT blkp, GMP *gmp, BOOL showIt, int page)
{
  if (! gmp->superDiv) /* initialize the system */
    {
      gmp->superDiv = page+10000 ;
      if (showIt) /* to see the yellow line */
	{
	  vtxtPrintf (blkp, "<span class='hh0'><a  class='showing'  id='%d' toggle='none'>"
		     "TABLE OF CONTENTS / OPEN CLOSE ALL PARAGRAPHS"
		     "</a></span><br/>\n"
		      , gmp->superDiv
		     ) ;
	}
      else /* just declare the div, and access it via javascript:openAllChapters */
	vtxtPrintf (blkp, "<div  class='showing' id='%d'></div>\n", gmp->superDiv) ;
      vtxtPrintf (blkp, "\n<div class='hiding' id='allChapters'></div>\n") ;
    }
  return TRUE ;
} /* gmpBlockInit */

/******************************************************************************/
/* notice that the block is not apparent in the chapter
 * so it can be started asynchroneously
 */
BOOL gmpBlockStart (vTXT blkp, GMP *gmp, const char *header)
{
  JUMPPOINTMENU *jmp ;
  /* open the chapter blocker */
  if (gmp->blockerDiv)
    return FALSE ;

  /* ATTENTION the _ in  id='_%s' toggle='_%s' is needed to distinguish these
     from <a href="#" name="%s">
     otherwise stupid internet explorer is totally confused
  */
  if (header && (jmp = gmpJumpPoint (blkp, gmp, header)))
    {
      if (gmp->superDiv < 2) /* block open */
	{
	  /* so that the summary is open */
	  vtxtPrintf (blkp, "\n<div class='showing' toggle='allChapters' id='c_%s'></div>\n"
		      , header) ;
	  vtxtPrintf (blkp, "\n<div class='shown' toggle='c_%s'>\n", header) ;
	}
      else  /* block closed */
	{
	  vtxtPrintf (blkp, "\n<div class='hiding' toggle='allChapters' id='c_%s'></div>\n"
		      , header) ;
	  vtxtPrintf (blkp, "\n<div class='hidden' toggle='c_%s'>\n", header) ;
	}
      jmp->used |= 4 ; 
    }
  gmp->blockerDiv = 1 ;
  return TRUE ;
} /* gmpBlockStart */

/******************************************************************************/
/* close open sections, captions. register */
BOOL gmpBlockClose (vTXT blkp, GMP *gmp, const char *header)
{ 
  if (gmp->markup && gmp->blockerDiv > 0)
    {
      vtxtPrintf (blkp, "</div c_%s>\n", header) ;  /* close the chapter blocker */
      gmp->blockerDiv = 0 ;
    }
  return TRUE ;
} /* gmpBlockClose */

/******************************************************************************/
/* close open sections, captions. register */
BOOL gmpChapterClose (vTXT blkp, GMP *gmp, const char *header, BOOL used)
{
  JUMPPOINTMENU *jmp ;
  if (gmp->markup && used && gmp->captionDiv + gmp->sectionDiv > 0)
    {
      while (gmp->captionDiv-- > 0) /* is a caption is open, close it */
	vtxtPrint (blkp, "\n      </div caption>") ;
      while (gmp->sectionDiv-- > 0) /* if a section is open, close it */
	vtxtPrint (blkp, "\n    </div section>") ;  
      vtxtPrint (blkp, "\n") ;
    } 
  if (header && (jmp = gmpJumpPoint (blkp, gmp, header)))
    { if (used) jmp->used |= 0x2 ; else jmp->used = 0 ;}
  gmp->captionDiv = gmp->sectionDiv = 0 ;   
  return TRUE ; 
} /* gmpChapterClose  */

/******************************************************************************/

BOOL gmpImportRemainder (vTXT blkp, GMP *gmp, AC_OBJ obj)
{
  if (gmp->markup && obj)
    {
      /* temporary table of content division */
      vtxtPrintf (blkp, "<div  id=\'%d\'>"   
		  , 100001
		  ) ; 
      
      vtxtPrintf (blkp, "<font color='red'> Please read the abstract !<br>The report should appear soon."
		  " The delay depends both on the speed of your connection and on the load of our server,"
		  " but we hope that it does not exceed 30s, please <a href=\"mailto:mieg?subject=problem loading the second part of the fiche of gene %s\">email</a> if it is much slower, thanks."
		  "</font>"
		  , ac_name(gmp->gene)
		  ) ; 
      
      vtxtPrintf (blkp,"\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
		  "  <!--\n" 
		  "  openAceViewChapter (\'%s\',\'%s\', \'%d\') ;\n "
		  "  //-->\n"
		  "</script>\n" 
		  , ac_class (obj), ac_name (obj), 100000
		  ) ;
      
      vtxtPrintf (blkp, "</div>\n") ;
    }
  return TRUE ;
} /* gmpImportRemainder */

/******************************************************************************/

BOOL gmpChapter2 (vTXT blkp, GMP *gmp, const char *header, const char *title, const char *bubble) 
{
  JUMPPOINTMENU *jmp;
  BOOL makeDiv = FALSE ;
  BOOL isDiag = FALSE ;

  if (*header=='#') { isDiag = TRUE ; header++ ;}
  if (*header=='*') { makeDiv = TRUE ; header++ ;}
  vtxtBreak (blkp) ;
  if (gmp->markup)
    {
      jmp = gmpJumpPoint (blkp, gmp, header) ;
      vtxtPrint (blkp, "\n") ;
      gmp->chapterHeader = "toto" ;
      vtxtPrintf (blkp
		  , "<a href='#' name='%s'></a>"
		  , header) ;
      {
	gmp->chapterDiv = jmp->index ;
	gmp->chapterHeader = header ;
	/* open the chapter */
	if (bubble)
	  vtxtPrintf (blkp, "\n  <span class='%s'><a  class='%s' toggle='%d' id='_%s' title='%s'>%s</a></span>\n"
		      , isDiag ? "hhD" : "hh1"
		      , jmp->isOpen ? "showing" : "hiding"
		      , gmp->superDiv 
		      , header 
		      , bubble 
		      , title
		      ) ;
	else
	  vtxtPrintf (blkp, "\n  <span class='%s'><a  class='%s' toggle='%d' id='_%s'>%s</a></span>\n"
		      , isDiag ? "hhD" : "hh1"
		      , jmp->isOpen ? "showing" : "hiding"
		      , gmp->superDiv 
		      , header 
		      , title
		      ) ;
	  

	if (1) /* show the up arrow */
	  vtxtPrint (blkp, "<a href='#top'><img SRC='images/arrowup.gif' border=0 width=14 height=14 ALT='back to top'></a></span>") ;
	if (jmp && *(jmp->help))  /* show 'explain' if available */
	  {
	    vtxtPrintf (blkp, "<span class='%s' toggle='_%s'>\n"
			, ! jmp || jmp->isOpen ? "shown" : "hidden"  
			, header
			) ;
	    gmpHelp (blkp, gmp, jmp->file, jmp->help) ; 
	    vtxtPrint (blkp, "</span>") ;
	  }
	  vtxtPrintf (blkp, "<br/>\n") ;	
	  
	  if (0) /* this trxt would appear at the top of the chapter */
	    {
	      vtxtPrint (blkp, "toto") ;
	      vtxtBreak (blkp) ;
	    }
	  if (makeDiv)
	    {
	      /* open a div toggled by the section and its parent chapter/superChapter 
	       * and initialised as open or closed according to the menu
	       */
	      
	      vtxtPrintf (blkp, "<div class='%s' toggle='_%s'>\n"
			  , ! jmp || jmp->isOpen ? "shown" : "hidden"  
			  , header
			  ) ;
	      gmp->sectionDiv = 1 ;
	    }
      }
    }
  else
    vtxtPrintf (blkp, "\n%s\n", title) ;
  return TRUE ;
} /* gmpChapter2 */

/******************************************************************************/

BOOL gmpChapter (vTXT blkp, GMP *gmp, const char *header, const char *title)
{
  return gmpChapter2 (blkp, gmp, header, title, 0) ;
} /* gmpChapter */

/******************************************************************************/

BOOL gmpSection2 (vTXT blkp, GMP *gmp, char *header, const char *title, const char *subTitle)
{
  BOOL isDiag = FALSE ;
  JUMPPOINTMENU *jmp ;
  
  if (*header=='#') { isDiag = TRUE ; header++ ;}
  vtxtBreak (blkp) ;
  if (gmp->markup)
    {
      while (gmp->captionDiv-- > 0) /* is a caption is open, close it */
	vtxtPrint (blkp, "\n      </div caption>\n") ;
      gmp->captionDiv = 0 ;
      while (gmp->sectionDiv-- > 0) /* if a section is open, close it */
	vtxtPrint (blkp, "\n    </div section>") ;  
      gmp->sectionDiv = 0 ;
      vtxtPrint (blkp, "\n") ;

      jmp = gmpJumpPoint (blkp, gmp, header) ;
      if (jmp && title)
	{
	  vtxtPrintf (blkp
		      , "<a href='#' name='%s'></a>"
		      , header) ;
	  vtxtPrintf (blkp, "<span class='%s'><a  class='%s' toggle='_%s' id='_%s'>%s</a></span>\n"
		      , isDiag ? "hh2D" : "hh2"
		      , jmp->isOpen ? "showing" : "hiding"
		      , gmp->chapterHeader ? gmp->chapterHeader : "toto"
		      , header 
		      , title
		      ) ;
	  if (1) /* show the up arrow */
	    {
	      vtxtPrintf (blkp, "<span class='closed' toggle='_%s'>" 
			  , gmp->chapterHeader ? gmp->chapterHeader : "toto"
			  ) ;
	      vtxtPrint (blkp, "<a href='#top'><img SRC='images/arrowup.gif' border=0 width=14 height=14 ALT='back to top'></a>") ;
	      vtxtPrintf (blkp, "</span>\n") ;
	    }
	  if (jmp && *(jmp->help))  /* show 'explain' if available */
	    {
	      vtxtPrintf (blkp, "<span class='%s' toggle='_%s'>\n"
			  , ! jmp || jmp->isOpen ? "shown" : "hidden"
			  , header
			  ) ;
	      gmpHelp (blkp, gmp, jmp->file, jmp->help) ;
	      vtxtPrint (blkp, "</span>") ;
	    }
	  vtxtPrintf (blkp, "<br/>\n") ;
	  if (subTitle)
	    {
	      vtxtPrintf (blkp, "%s", subTitle) ;
	      vtxtBreak (blkp) ;
	    }

	  /* open a div toggled by the section
	   * and initialised as open or closed according to the menu
	   */
	  vtxtPrintf (blkp, "<div class='%s'  toggle='_%s'>\n"
		      , ! jmp || jmp->isOpen ? "shown" : "hidden"
		      , header
		      ) ;
	  gmp->sectionDiv = 1 ;  
	}
    }
  else
    if (title) vtxtPrintf (blkp, "\n%s\n", title) ;
  return TRUE ;
} /* gmpSection2 */

/******************************************************************************/

BOOL gmpSection (vTXT blkp, GMP *gmp, char *header, const char *title)
{
  return gmpSection2 (blkp, gmp, header, title, 0) ;
} /* gmpSection */

/******************************************************************************/
/* caption: a title a help and an embedded open/close */
BOOL gmpCaption (vTXT blkp, GMP *gmp, char *header, const char *title)
{
  JUMPPOINTMENU *jmp ;
  
  vtxtBreak (blkp) ;
  if (gmp->markup)
    {
      while (gmp->captionDiv-- > 0) /* is a caption is open, close it */
	vtxtPrint (blkp, "\n      </div caption>\n") ;
      gmp->captionDiv = 0 ;
      vtxtPrint (blkp, "\n") ;

      jmp = gmpJumpPoint (blkp, gmp, header) ;
      if (jmp && title)
	{
	  vtxtPrintf (blkp
		      , "<a href='#' name='%s'></a>"
		      , header) ;
 
	  vtxtPrintf (blkp, "<span class='hhcaption'><a class='hiding' id='_%s'>%s</a>\n"
		      , header
		      , title
		      ) ;
	  vtxtPrintf (blkp, "</span><br/>\n") ;

	  /* open a div toggled by the caption and but not its parent chapter/superChapter 
	   * and initialised as open or closed according to the menu
	   */
	  vtxtPrintf (blkp, "<div class='%s' toggle='_%s'>\n"
		      , ! jmp || jmp->isOpen ? "shown" : "hidden"
		      , header
		      ) ;
	  gmp->captionDiv = 1 ; /* a caption is part of the section div */
	}
    }
  else if (title) 
    vtxtPrintf (blkp, "\n%s\n", title) ;
  return TRUE ;
} /* gmpCaption */

/* a title and a help but no go up arrow and no open close */
BOOL gmpSubSection (vTXT blkp, GMP *gmp, char *header, const char *title)
{
  vtxtBreak (blkp) ;
  if (gmp->markup)
    {
      JUMPPOINTMENU *jmp = gmpJumpPoint (blkp, gmp, header) ;
      vtxtPrintf (blkp
		  , "<a href='#' name='%s'></a>"
		  , header) ;
 
      if (title)
	{
	  vtxtPrintf (blkp, "<span class='hh3'>%s</span>\n", title) ;
	  if (jmp && *(jmp->help))
	    gmpHelp (blkp, gmp, jmp->file, jmp->help) ;
	  vtxtPrint (blkp, "<br/>\n") ;
	}
    }
  else
    if (title) 
      vtxtPrintf (blkp, "\n%s\n", title) ;
  return TRUE ;
} /* gmpSubSection */

BOOL gmpImage (vTXT blkp, GMP *gmp, AC_OBJ obj)
{
  if (gmp->markup)
    {
      vtxtPrint (blkp, "<script language=\'JavaScript\' type=\'text/javascript\'>\n") ;
      vtxtPrint (blkp, "  <!--\n") ;
      vtxtPrintf (blkp,
		  "ficheObjectImage(\"%s\",\"%s\");\n"
		  ,ac_class(obj), ac_name(obj)
		  ) ;
      vtxtPrint (blkp, "  //-->\n") ;
      vtxtPrint (blkp, "</script>\n") ;
    }
  return TRUE ;
}  /* gmpImage */

int gmpURL (vTXT blkp, GMP *gmp, char *link, const char *text) 
{
  int nn ;
  
  if (gmp->markup)
    nn = vtxtPrintf (blkp,"<a href=\"%s\" target=\"_top\">%s</a>",link, text ? text : link) ;
  else
    nn = vtxtPrint (blkp,  text ? text : link) ;
  
  return nn ;
} /* gmpURL */  

BOOL gmpAction (vTXT blkp, GMP *gmp, AC_OBJ obj, const char *action, const char *text)
{
  if (gmp->markup)
    {
      vtxtPrintf (blkp, "<a href=\"javascript:openAceViewAction ('%s', '%s', '%s') ; \" >%s</a>"
			 , ac_class(obj), ac_name (obj), action
			 , text) ;
    }
  else
    vtxtPrint (blkp, text) ;
  return TRUE ;
}  /* gmpAction */

static int gmpFakeObjLinkAnchor (vTXT blkp, GMP *gmp, const char *cl, AC_OBJ obj
				 , const char *anchor, const char *text) 
{
  int nn = 0 ;
  
  if (obj && !text) text = ac_name (obj) ;
  if (obj && gmp->markup)
    {
      if (!cl && obj) cl = ac_class(obj) ;
      if (obj && anchor)
	nn = vtxtPrintf (blkp, "<a href=\"javascript:openAceViewLinkAnchor ('%s', '%s', '%s', '%s') ; \" >"
			 , cl, ac_name (obj), anchor
			 , strncasecmp (cl, "dna", 3) ? (strcasecmp (cl, "gene") ? "Fiche" : "fgene") : "DNA") ;
      else if (obj)
	nn = vtxtPrintf (blkp, "<a href=\"javascript:openAceViewLink('%s', '%s', '%s') ; \" >"
			 , cl, ac_name (obj)
			 , strncasecmp (cl, "dna", 3) ? (strcasecmp (cl, "gene") ? "Fiche" : "fgene") : "DNA") ;
      else if (anchor)
	nn = vtxtPrintf (blkp, "<a href=\"#%s\">", anchor) ;
      if (text) vtxtPrint (blkp,  text) ;
      vtxtPrint (blkp,"</a>") ;
    }
  else
    nn = text ? vtxtPrint (blkp,  text) : 0 ;
  
  return nn ;
} /* gmpFakeObjLink */  

int gmpFakeObjLink (vTXT blkp, GMP *gmp, const char *cl, AC_OBJ obj, const char *text) 
{
  return gmpFakeObjLinkAnchor (blkp, gmp, cl, obj, 0, text) ;
} /* gmpFakeObjLink */  
int gmpObjLink (vTXT blkp, GMP *gmp, AC_OBJ obj, const char *text) 
{
  return gmpFakeObjLinkAnchor (blkp, gmp, 0, obj, 0, text) ;
} /* gmpObjLink */  

int gmpObjLinkAnchor (vTXT blkp, GMP *gmp, AC_OBJ obj, const char *anchor, const char *text) 
{
  return gmpFakeObjLinkAnchor (blkp, gmp, 0, obj, anchor, text) ;
} /* gmpObjLink */  

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

