#ifndef  BIOLOG_H_DEF
#define BIOLOG_H_DEF

/* #include    <acec.h> */
#include    "../wac/ac.h"
#include   "vtxt.h"
#include "dict.h"
#include "keyset.h"

/* sentence composition macros */
#define _verbs(v_src)    ((v_src)>1 ? "" : "s") 
#define _multi(v_src)    ((v_src)>1 ? "s" : "") 
#define _isare(v_src)    ((v_src)>1 ? "are" : "is") 
#define _theese(v_src)    ((v_src)>1 ? "these" : "this") 
#define _isone(v_src)    isOne(v_src)   /* one, 2, 3, 4 ... */
#define _haves(v_ind)    ((v_ind)>1 ? "have" : "has") 
#define _comma(v_ind)    ((v_ind)==0 ? "" : ", ") 
#define _waser(v_ind)    ((v_ind)>1 ? "were" : "was") 
#define _dodoes(v_ind)    ((v_ind)>1 ? "do" : "does") 
#define _dodoes(v_ind)    ((v_ind)>1 ? "do" : "does") 

/* markup */
enum    {_sMarkup=0,
	 _PN,_PN2,_PB,_SN,_SN2,_SB,_TB,_TR,_TD,_SE,_BR,_LI,_UL
};					   

void ficheSectionizeFiche (char * ficheComments,char ** titleFromFiche, char **commentSection);
extern char genomeRelease[];

typedef struct {
    char * BlastBuf,*PfamBuf,*PsortBuf,* TaxBuf,* ExpasyBuf;
    }rawFICHE;

void    kantorRawCallectRaw(vTXT blkp,char * keyname,char * sequence,rawFICHE * praw,char * todo,vTXT  log);
void 	kantorRawRemoveFiles(char * keyname,char * todo);

int kantorParseBLAST(vTXT blkp,vTXT  hintBlk,char * blastBuf,char* separ);
int kantorParseTAXBLAST(vTXT blkp,char * taxblastBuf);
int kantorParsePFAM(vTXT blkp,vTXT  hintBlk,char * pfamBuf);

char * kantorGenerateTitle(char * dbName,char * keyname);

char * swormView(vTXT blkp,char * dbName,char * clsName,char * objName,char style, char *view, char *params);


enum {
	HUMAN=0,
	WORM,
        ARA,
        DROSO,
	COLI,
        MOUSE,
        RAT,
	SPECIESCNT
};

#define MESH_LINK "https://www.nlm.nih.gov/cgi/mesh/2007/MB_cgi?term=%s"
#define COG_LINK "https://www.ncbi.nlm.nih.gov/cgi-bin/COG/palox?%s"
#define OLD_GO_LINK "http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=details&depth=1&query=%s"
#define GO_LINK "http://www.ebi.ac.uk/ego/DisplayGoTerm?id=%s"
#define PFAM_LINK "http://ebi.ac.uk/interpro/entry/pfam/%s"
#define PSORT_LINK "http://psort.hgc.jp/%s"
#define KEGG_LINK_HUMAN "http://www.genome.jp/dbget-bin/show_pathway?hsa%s+%s"
#define KEGG_LINK_MOUSE "http://www.genome.jp/dbget-bin/show_pathway?mmu%s+%s"
#define KEGG_LINK_RAT "http://www.genome.jp/dbget-bin/show_pathway?rno%s+%s"
#define KEGG_LINK_ARA "http://www.genome.jp/dbget-bin/show_pathway?ath%s+%s"
#define BIOCARTA_LINK  "http://www.biocarta.com/pathfiles/%sPathway.asp"
#define OMIM_LINK "https://www.ncbi.nlm.nih.gov/omim/%s"
#define TAIR_LINK "http://www.arabidopsis.org/servlets/TairObject?type=locus&name=%s"
#define MGI_LINK "http://www.informatics.jax.org/searches/accession_report.cgi?id=MGI:%s"
#define MGI_ALLELE_LINK "http://www.informatics.jax.org/searches/allele_report.cgi?markerID=MGI:%s"

#define PUBMED_LINK      "https://www.ncbi.nlm.nih.gov/pubmed/"
#define PUBMED_GENEID_LINK      "https://www.ncbi.nlm.nih.gov/pubmed?linkname=gene_pubmed&from_uid=%s"

#define GENBANK_LINK "https://www.ncbi.nlm.nih.gov/nucest/%s"
#define ENTREZ_LINK "https://www.ncbi.nlm.nih.gov/gene/%s"
#define RAT_LINK "https://www.ncbi.nlm.nih.gov/gene/%s"
   

typedef struct tagSPECIESINFO {
  char * speciesName,* speciesGenomicName;	
  int		long5P,long3P;
}SPECIESINFO;

extern	SPECIESINFO SpcI[];
extern	int Spc;

typedef struct { AC_DB db ; AC_OBJ chrom, gene, tg, mrna, pg, est, clone, product, kantor, oNM ; int nMrna ; char *dna, *peptide, *variant ; char style, view, pictureType ; BOOL useAm, markup ; SPECIESINFO *spci ; int Spc, captionDiv, sectionDiv, chapterDiv, blockerDiv, superDiv, requestedChapter ; AC_HANDLE h ; const char *chapterHeader ; } GMP ;

GMP *gmpCreate (AC_DB db, AC_OBJ oGene, AC_OBJ oTg, AC_OBJ oMrna, AC_OBJ oGF, AC_OBJ oProduct, char style, char view) ;
void gmpDoDestroy (GMP *gmp) ;
#define gmpDestroy(_gmp) (gmpDoDestroy(_gmp), (_gmp)=0)

char *fAsnGenerateMRNA (vTXT blkp, GMP *gmp, char * ficheComments, int isHeader) ;
char * fAsnGenerateGene (AC_DB db, AC_OBJ oGene, int isDetail, char style);
char * ficheChromosomeDump (vTXT vtxt, AC_DB db, AC_OBJ oMap, char style);


int fogicTaxTree (vTXT blkp, GMP *gmp, AC_OBJ gTax_tree, int level);
void fogicTaxCrit (vTXT blkp, GMP *gmp, AC_OBJ gTax_tree, int level, int * cnt);

char * kantorGetTitles (vTXT blkp, char * kantorName, int isComplete);
char *kantorReportTitles (char *productName) ;
void ficheMRNA5PrimeParagraphContent (vTXT blk, GMP *gmp) ;

int ficheDBCreatureID (AC_DB db);

void ficheGmpMapDescriptionParagraph (vTXT blkp, GMP *gmp, int type) ;


int ficheMrnaTitleParagraph (vTXT blkp, GMP *gmp) ;
int ficheProductTitleParagraphContent (vTXT blk, GMP *gmp);
int ficheTGAntisensParagraphContent (vTXT blkp, GMP *gmp, BOOL details) ;
char * fogicArguableCommonAncestors (const char * namKantor);
void ficheGmpMainSupportingClones  (vTXT blkp, GMP *gmp) ;
int ficheCloneParagraphContent (vTXT blkp, GMP *gmp) ;
void ficheGmpFunctionPart (vTXT blkp, GMP *gmp) ;
char *fichePrimersParagraphContent (vTXT blkp, GMP *gmp) ;

void ficheMRNATranscriptParagraph (vTXT blkp, GMP *gmp);
void ficheMRNAConceptualTranslationParagraph (vTXT blkp, GMP *gmp);
void fichePRODUCTPsortParagraph (vTXT blkp, GMP *gmp);
void fichePRODUCTPfamParagraph (vTXT blkp, GMP *gmp);
void fichePRODUCTBlastPParagraph (vTXT blkp, GMP *gmp);
void fichePRODUCTTaxTreeParagraph (vTXT blkp, GMP *gmp) ;
void ficheMRNA5PrimeParagraph (vTXT blkp, GMP *gmp);
void ficheMRNA3PrimeParagraph (vTXT blkp, GMP *gmp);
void ficheMRNASplicingParagraph (vTXT blkp, GMP *gmp);

void fichePEPTIDEParagraph (vTXT blkp, GMP *gmp) ;
void ficheDNAParagraph (vTXT blkp, GMP *gmp);
void ficheTGOverviewParagraph (vTXT blkp, GMP *gmp);

void ficheGenomeSummaryChapter (vTXT blkp, GMP *gmp);
void ficheGenomePlotChapter (vTXT blkp, AC_DB db, GMP *gmp, BOOL isBig);

void ficheGENEFunctionParagraph (vTXT blkp, GMP *gmp);
void ficheGENENeighborhoodParagraph (vTXT blkp, GMP *gmp);
char *ficheGeneBiologyPart (vTXT blkp, GMP *gmp);
void ficheGeneOverviewParagraph (vTXT blkp, GMP *gmp);
void ficheGENEmrnaListParagraph (vTXT blkp, GMP *gmp);
void ficheGENEIntronExonListParagraph (vTXT blkp, GMP *gmp);

void ficheProductBlastPTableParagraph (vTXT blkp, GMP *gmp);
void ficheProductBlastPTableContent (vTXT blkp, GMP *gmp, int maxCol, const char *more) ;
void ficheGENEPhenotypeAndDescriptorParagraphContent (vTXT blkp, GMP *gmp, int whatToOut);
void ficheGENERegulationParagraphContent (vTXT blkp, GMP *gmp) ;
void ficheMRNAConceptualTranslationParagraphContent (vTXT blkp, GMP *gmp);
void fichePRODUCTPsortParagraphContent (vTXT blkp, GMP *gmp, BOOL justPsortLocalisation);
char *fichePRODUCTBlastPParagraphContent (vTXT blkp, GMP *gmp) ;
char *ficheTAXClosestAncestorStatement (GMP *gmp) ;
char *fichePRODUCTTaxTreeParagraphContent (vTXT blkp, GMP *gmp) ;

void ficheMRNA3PrimeParagraphContent (vTXT blkp, GMP *gmp) ;
void ficheMRNAGeneDescriptionPart (vTXT blkp, GMP *gmp, int isMRNA);
void ficheMRNASupportingClonesPart (vTXT blkp, GMP *gmp) ;
void ficheMRNAProteinAnalysisPart (vTXT blkp, GMP *gmp) ;
void ficheMRNAStructurePart (vTXT blkp, GMP *gmp) ;
void ficheMRNASequencesDescriptionPart (vTXT blkp, GMP *gmp) ;
void ficheTGBestNameStatement (vTXT blkp, GMP *gmp, AC_OBJ tg, int isLocalLink) ;
char *ficheMrnaGenomeDNA (vTXT blkp, AC_DB db, AC_OBJ oMrna, char style, BOOL isPromotor) ;

/* new stuff 2004 */
void ficheNewGeneTitleParagraph (vTXT blkp, GMP *gmp) ;
char *ficheNewGeneAceViewSummary (vTXT blkp, GMP *gmp) ;
void ficheNewGeneSummaryChapter (vTXT blkp, GMP *gmp) ;

void ficheNewGeneExpressionProfileChapter (vTXT blkp, GMP *gmp) ;
void ficheNewGeneLocatorDiagramChapter (vTXT blkp, GMP *gmp, BOOL isSmall) ;
void ficheNewGeneGenomeDiagramChapter (vTXT blkp, GMP *gmp) ;
void ficheNewGeneCompactDiagramChapter (vTXT blkp, GMP *gmp, int pass) ;
void ficheNewGeneWiggleDiagramChapter (vTXT blkp, GMP *gmp, int pass) ;
void ficheNewGeneSequenceChapter (vTXT blkp, GMP *gmp) ;
/* void ficheNewGeneCompactDiagramLink (vTXT blkp, GMP *gmp) ; */
/* void ficheNewGeneWiggleDiagramLink (vTXT blkp, GMP *gmp) ; */
void ficheNewGeneMrnaDiagramChapter (vTXT blkp, GMP *gmp) ;

void ficheNewGeneMappingLinksAliasesChapter (vTXT blkp, GMP *gmp) ;
void ficheNewGenePhenotypeFunctionChapter (vTXT blkp, GMP *gmp) ;
void ficheNewGeneExpressionChapter (vTXT blkp, GMP *gmp) ;
void ficheNewGeneAltFeatureChapter (vTXT blkp, GMP *gmp) ;
void ficheNewGeneMOLECULESChapter (vTXT blkp, GMP *gmp) ;
void ficheNewGeneBiblioChapter  (vTXT blkp, GMP *gmp, int gRif) ;
char *ficheNewGeneExpressionTissue (vTXT blkp, GMP *gmp, AC_KEYSET clones, int nMax1, int nMax2, char showOther, int nClo) ;

BOOL ficheNewMrnaSummaryChapter (vTXT blkp, GMP *gmp) ;
void ficheNewMrnaStructureChapter (vTXT blkp, GMP *gmp) ;
void ficheNewMrnaFlashDiagramChapter (vTXT blkp, GMP *gmp) ;
BOOL ficheNewMrnaExpressionCloneSupportChapter (vTXT blkp, GMP *gmp) ;
BOOL ficheNewTgSupportChapter (vTXT blkp, GMP *gmp, AC_KEYSET ks1, const char *title) ;
BOOL ficheNewMrnaProteinChapter (vTXT blkp, GMP *gmp) ;
BOOL ficheNewMrnaSequenceChapter (vTXT blkp, GMP *gmp) ;
char *ficheNewMrnaSubmissionComment (vTXT blkp, GMP *gmp) ;
void ficheNewMrnaDecoratedPeptide (vTXT blkp, GMP *gmp, BOOL hasStop) ;
void ficheListAndMarkUpMrnas (vTXT blkp, GMP *gmp, char type, BOOL fromGene) ;

int ficheNewAceKogStatement (vTXT blkp, GMP *gmp, BOOL doTitle) ;

int ficheNewCloneTable (vTXT blkp, GMP *gmp, AC_KEYSET clones, char orderBy, int maxLine, const char *more, Array tissues) ;
void ficheNewCloneParagraph (vTXT vtxt, AC_DB db, AC_KEYSET clones, char style, char orderBy, Array tissues) ;

void ficheGoldCloneParagraph (vTXT blkp, GMP *gmp) ;
#endif  /* BIOLOG_H_DEF */


