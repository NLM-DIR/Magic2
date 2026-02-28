/*
 * authors: Danielle and Jean Thierry-Mieg, NCBI, 
 * March 2012
 * geneindex.cmain_clone_ncbi.ace   
 *   Input: ace file containing run/gene RNA-seq counts
 */

#include "regular.h"
#include <ac.h>
#include <bigarray.h>
#include <acedna.h>
#include "bitset.h"
#include "topology.h"
#include <math.h>

#define LEMING 0
#define NVIRTUALSTRATA 4

static double FIX = 0 ;

typedef struct dcStruct { 
  float index, variance, divariance, indexMin, indexMax ; int N ; int isLow ;
  double seqs, tags, kb ;
} DC ;
typedef struct bigGeneStruct { int gene ; float kb ;} BG ;
typedef struct ddcStruct { 
  double nReads, nReadsOk ; 
  long int nerr ; int a2g, partial, orphan, badTopo, multi, multi2 ; 
} DDC ;
typedef struct pcStruct { 
  int gene, ok ; double index1, index2, dIndex, *idx, var1, var2, dvar, sig, d2dw, overlap, danielle, AUC, geneid ; 
} PC ;
typedef struct rwStruct { int run ; double alpha, beta1, beta2 ; } RW ;  /* Run-Weighted struct,  used to construct iterated groups */
typedef struct endpointStruct { int run1, run2, run3, run4 ; }  ENDPOINT ;
typedef struct compareStruct { 
  AC_HANDLE h ;
  int compare ; 
  KEYSET runs ; 
  KEYSET predicts ; 
  KEYSET ratios ;
  Array shape ;
  BOOL correlation, autosort, show_all_runs, private ;
  BOOL samplePairing, expressionProfile, snpProfile, showAllHistos ;
  BOOL compareSNP ;
  BOOL compareINTRON ;
  BOOL who_is_who ;
  BOOL profile ;
  BOOL mixing_ratios ;
  BOOL done ;
  int vRun1, vRun2 ;
  int nDEG[3] ;
  int nDEG185[3] ;
  int nGroupDEG[3] ;
  int nGroupDEG185[3] ;
  float fdr[3] ; 
  int fdrThreshold[3] ; 
  float fdr1[201], fdr2[201] ;
  int nVirtualStrata ; /* number of virtual strata actually generated */
  int forceIbest ;
  KEYSET genePlus[3] ;
  Array noise, scores[3] ; /* array of differential scores and noise for all genes */
} COMPARE ;

typedef struct rcStruct { 
  int run, runid, sample, machine, title, otherTitle, sortingTitle, sortingTitle2, nRuns ; 
  int Genes_with_index,  Genes_with_reliable_index,  Genes_with_one_reliable_index, Genes_with_index_over_8, Genes_with_index_over_10, Genes_with_index_over_12, Genes_with_index_over_15, Genes_with_index_over_18, Genes_with_index_over_21, Genes_with_index_over_24, Genes_with_index_over_27, Genes_with_index_over_30 ;
  int captured_genes_with_index,  captured_genes_with_reliable_index,  captured_genes_with_one_reliable_index, captured_genes_with_index_over_8, captured_genes_with_index_over_10, captured_genes_with_index_over_12, captured_genes_with_index_over_15, captured_genes_with_index_over_18, captured_genes_with_index_over_21, captured_genes_with_index_over_24, captured_genes_with_index_over_27, captured_genes_with_index_over_30 ;
  double Low_index, zeroIndex, NA_crossover_index, geneTags, genomicKb, intergenicKb, intergenicDensity, targetKb, prunedTargetKb, bigGenesKb, capturedBigGenesKb, anyTargetKb, X2, alpha, beta1, beta2, fineTune ; 
  double seqs, tags, kb ;
  double captured_seqs, captured_tags, captured_kb ;

  double nReads, nReadsOk ; 
  long int nerr ; int a2g, partial, orphan, badTopo, multi, multi2 ; 
  double CnReads, CnReadsOk ; 
  long int Cnerr ; int Ca2g, Cpartial, Corphan, CbadTopo, Cmulti, Cmulti2 ; 

  Array aa, daa, aag, antiAa ;
  KEYSET runs, compare_to, titration ;
  Array bigGenes, capturedBigGenes ;
  BOOL isSNP, cumulDone, addCounts, isAny, solid, hasData, avoid, PolyA_selected_RNA, isPairedEnd, selectedVariance, snpCompare ;
  int indexDone, fragmentLength, private, isSubLib ;
  int accessibleLength5kb ;  /* max gene length effectivelly sequenced, evaluated on transcripts > 8kb, usually limited by the 3p biais */
  int accessibleLength8kb ;  /* max gene length effectivelly sequenced, evaluated on transcripts > 8kb, usually limited by the 3p biais */
  int isMA ; /* probe length */
  int pA ; /* 0, 1, 2 overrides gx->pA  */
  BOOL DGE ; /* probe length */
  int groupLevel ;
  int capture ;
  float minIntronSupport ;
  float accepted_tags ;
  float shift_expression_index ;
  float abaque[31] ;
} RC ;

typedef struct gcStruct { 
  int title, affy, geneId, geneModel, intronGene, intronMrna, capturedIntronGene, nmid, gene_type, fromGene, fromTranscript, length, boxLength, chrom, a1, a2, intronType, intronSharedD, intronSharedA, alias ; 
  BOOL isSNP, notIsTranscriptA, targeted, captured, hasGoodProduct, isGood ;
  Array geneGroup ;
  KEYSET captures, capturesTouch ;
  int Group_level ;
  int gMix ;
  BOOL hasWeights ; /* group has weight, we may use it as a predictor (not yet done) compute a log proba log(exp( w X) / (  1 + exp (w X))) */
  int HGNC_family ;
  int intronAv, intronRefSeq, hasKnownDonorAcceptor ;
  float minIndex, maxIndex, averageIndex, tags, score ; 
  float weight ; /* weight of a gene in the gene group */
  float score1, score2, bonus, av1, av2 ;     /* half scores and zones when analyzing the comparative histos */
  float intrinsicSigma ; /* inherited from SEQC_main: indicates a gene hard to quantify */
  int ratioa, ratiob ;
  int ix1, ix2, iy1, iy2 ;
  int za1, za2, zb1, zb2 ;
  char strand, isIntron, isLow1, isLow2 ; 
} GC ;


typedef struct ggStruct { 
  int chrom, a1, a2 ;
  KEYSET genes ;
} GG ;

typedef struct gggStruct { 
  int gene ;
  float weight ;
} GGG ;

typedef struct scalStruct { 
  int N ;
  double X, Y ;
} SCAL ;

typedef struct gxStruct { 
  AC_HANDLE h ; 
  AC_DB db ;           /* ACEDB database containing the semantics of the experiment */
  const char *title ;
  const char *deepFileName ;
  const char *outFileName ;
  const char *runListFileName ;
  const char *runAceFileName ;
  const char *maskFileName ;
  const char *referenceGenome ;
  const char *chromAliasFileName ;
  const char *stableGenesFileName ;
  const char *markerSelectFileName ;
  const char *markerRejectFileName ;
  const char *splitMrnaFileName ;
  const char *select ;
  const char *method ;
  const char *target_class ; /* import nh_Ali for that target_class */
  const char *seedGene ;
  const char *genePlusFileName ;
  const char *geneMinusFileName ;
  const char *geneGroupFileName ;
  const char *geneClustersFileName ;
  const char *snpFileName ; 
  int isMA ;
  int pA ;
  BOOL DGE ;
  int geneSelectionMethod ;

  BOOL gzi, gzo ;
  BOOL showAnyRun ;
  BOOL unique ;
  BOOL anti ; /* count the Anti_run */
  BOOL fineTune, iterate, snpEval, snpCompare ;
  BOOL phenoGene ;
  BOOL targeted ; /* only export targeted genes */
  const char *captured ; /* only export genes with tag "Capture <captured>" */
  BOOL TGx ; /* TGx strategy, use discrete scale */
  int digital, digitalCorrelation ;
  int normalize ;
  int localIteration ;
  Array runs ;           /* array of RC */
  DICT *runDict ;
  DICT *captureDict ;
  Array genes ;           /* array of GC */
  Array geneIds ;           /* array of GC */
  DICT *geneDict ;
  DICT *geneIdDict ;
  DICT *maskDict ;
  DICT *chromDict ;
  DICT *chromAliasDict ;
  DICT *compareDict ;
  Array splitMrnaArray ;
  DICT *splitMrnaDict ;
  Associator splitMrnaAss ;
  KEYSET chromAliasKs ;
  KEYSET geneLns ;      /* gene lengths */   
  Stack info ;
  DICT *affyDict ;
  DICT *userGivenGenePlusMinusDict ;
  KEYSET affy2gene ;
  KEYSET lowVarianceGenes ;
  KEYSET userGivenGenePlus, userGivenGeneMinus,  userGivenGenePlusMinus ;
  Array signatureGenes, signatureGenesSorted, genePlus, endpoints ;
  int malus, minIndex, seaLevel, removeLimit, wall ;
  int maxGroupLevel ;
  int maxGeneGroupLevel ;
  int runAny, runAnySolid, runAnyILM, isRunMax, nHistoRun ;         /* global counts */
  int subsample ; /* max number of Million reads per run */
  float ratio_bound, minVariance ;
  float genomeLengthInKb ;
  int histo_shifting, maxGenePlus ;
  BOOL isINTRON, isTranscript, isSNP, isMRNAH, isMicroRNA, keepIndex, medianCenter, noHisto, exportDiffGenes, hasSelectedVariance, skipEmptyGenes, deepTsf ;
  BOOL hasRunId, hasRunSample, hasMachine, hasRunTitle, hasRunOtherTitle, hasRunSortingTitle, hasRunSortingTitle2, hasGroup, hasCompare, hasTitration ;
  BOOL hasGeneAffy, hasGeneId, hasGeneModel, hasGeneNm, hasGeneType, hasGeneTitle, hasGeneAlias, hasGeneLength, hasGeneBoxLength,hasGeneChrom, hasIntronType, hasIntronShared,hasKnownDonorAcceptor,hasFromGene, hasFromTranscript, hasRefSeqAv, hasGoodProduct, hasCapture, hasCaptureTouch ;
  BOOL hasAccessibleLength ;
  ACEOUT aoFineTune, aoAUC2 ;

  BOOL gene2deep2index, intronGeneMix, hasDaa ;
  BOOL correlation, compare_to, samplePairing, expressionProfile, snpProfile, showAllHistos, html, who_is_who, mixing_ratios ;
  float threshold, minFoldChange, minWhoScore ;
  const char *export ;
  const char *htmlSpecies ;

  Array XY, seedRuns ;
  Array geneClusters ;
  double gauss5[21] ;
  double gauss10[21] ;
  double gauss20[41] ;
  double gauss30[61] ;
  double gauss40[81] ;
  double AUC, MCC, AUC2, MCC2 ;
  const char *isLowTitle[3] ;
  Array degHistos ;
  Array baddyBatchies ; /* raise the threshold for genes coming in the random strata */
  Array baddyBatchiesAll ; /* register all random strata values, so we can later subtract the average before we compute the FDR histograms */

  Array gza ; /* used to sort the exported index tables */
  Array compares ; /* classe compare: protocol classe for titration etc */
  Array pairScores ;  /* memorize the pair score when computing the .compare files, use them in the .profile files */
  Array pValues ;  /* memorize the pValue of the pair score when computing the .compare files, use them in the .profile files */
  Array pairs2cmp ;
  DICT *pairDict ;  /* memorize the compare corresponding to a pair of runs */
  DICT *pairScoresDict ;  /* memorize the pair score when computing the .compare files, use them in the .profile files */
  KEYSET pairScores2compare ;
  KEYSET targeted_genes ;
  Array diffGenes ;
} GX ;


static int pcOrder (const void *a, const void *b) ;
static Array gxAllRunsCosineUsingGenesRun1Run2 (GX *gx, int run1, int run2, KEYSET genePlus, KEYSET geneMinus, Array sp, Array sm, AC_HANDLE h) ;
static BOOL gxRunsCosineUsingGenes (GX *gx, KEYSET genes, Associator assPlus, Associator assMinus, Array sp, Array sm, int run1, int run2, double *alphap, int *nnp, int lowAll, BOOL isDiff) ;
static int gxExportComparison (GX *gx, COMPARE *compare, Array rws, KEYSET genePlus, KEYSET geneMinus, int run1, int run2, BOOL histo, BOOL wantGeneGroup) ;
static void gxExportSignatures (GX *gx, int run1, int run2, int run3, int run4, int pass) ;
static Array gxTrueFalse_AUC_MCC (GX *gx, ACEOUT ao, Array rws, int runa, int runb, int runa1, int runb1, int pass) 
  ;
static int gxDghInitializeBaddyBatchies (GX *gx, AC_HANDLE h) ;
static void gxDghRegisterScore (GX *gx, int iVV, PC *pc) ;
static void gxDghRegisterHisto (GX *gx, int run1, int run2, int nRun1, int nRun2) ;
static void gxDghExportHistos (GX *gx, COMPARE *compare, int run1, int run2, int isPm) ;
static int gxExportDiffGenes (GX *gx, COMPARE *compare, Array pp, int run1, int run2, KEYSET vA, KEYSET vB, int isPm) ;
static void gxOneGroupCumul (GX *gx, RC *rc, int myGene) ;
static BOOL gxSnpEvaluate (GX *gx, int gene, const char *snpName) ;
static void gxExportTableHeader  (GX *gx, ACEOUT ao, int type) ;

static void gxDghShiftScores (GX *gx, COMPARE *cmp) ;
static void gxDoRegisterComparedGenes (GX *gx, int gene, int iCmp, int iRun, int jRun, float iIndex, float jIndex, float diffScore, float pValue) ; 
static BOOL gxMakeComparativeHistos (GX *gx, int gene
				       , int run1, int run2, KEYSET vA, KEYSET vB
				       , double *hh0s, double *hh1s, float *chi2p
				       , int *nR1p, int *nR2p, int *nHeterop, int *nHomop, int *nObsp
				       , double *thresholdp, double *Ginip, double *pValuep
				       , BOOL isVirtual) ;



/*************************************************************************************/
/*************************************************************************************/

typedef struct splitMrnaStruct { 
  int gene, mrna, gMix
    , x1, x2
    ;
} SPLITMRNA ;

/*************************************************************************************/
/*************************************************************************************/

typedef struct gzaStruct { int gene ; float index ; } GZ ;
static int gzOrder (const void *a, const void *b)
{
  const GZ *up = (const GZ *)a, *vp = (const GZ *)b ;
  float z = up->index - vp->index ;

  if (z > 0) return -1 ;
  if (z < 0) return 1 ;
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/

static int endPointOrder (const void *a, const void *b)
{
  const ENDPOINT *up = (const ENDPOINT *)a, *vp = (const ENDPOINT *)b ;
  int n ;

  n = up->run1 - vp->run1 ; if (n) return n ;
  n = up->run2 - vp->run2 ; if (n) return n ;
  n = up->run3 - vp->run3 ; if (n) return n ;
  n = up->run4 - vp->run4 ; if (n) return n ;

  return 0 ;
} /* endpointOrder */

/*************************************************************************************/

static double hhShow (double *hh, BOOL show)
{
  double cumul = 0 ;
  if (hh)
    {
      for (int i = 0 ; i <256 ; i++)
	cumul += hh[i] ;
      if (show)
	for (int i = 0 ; i <256 ; i++)
	  if (hh[i]>0) 
	    fprintf (stderr, "%d\t%.2g\n", i, hh[i]) ;
      fprintf (stderr, "cumul\t%.2g\n", cumul) ;
    }
  return cumul ;
}

/*************************************************************************************/
/*************************************************************************************/

static void gxParseInit (GX *gx)
{
  gx->geneLns = arrayHandleCreate (1000, float, gx->h) ;
  gx->info = stackHandleCreate (10000, gx->h) ;
  gx->runs = arrayHandleCreate (500, RC, gx->h) ;
  gx->genes = arrayHandleCreate (50000, GC, gx->h) ;
  gx->geneIds = arrayHandleCreate (50000, GC, gx->h) ;
  gx->runDict = dictHandleCreate (500, gx->h) ;
  gx->captureDict = dictHandleCreate (32, gx->h) ;
  if (gx->maskFileName) gx->maskDict = dictHandleCreate (500, gx->h) ;
  if (gx->chromAliasFileName)
    {
      gx->chromAliasDict = dictHandleCreate (500, gx->h) ;
      gx->chromAliasKs = keySetHandleCreate (gx->h) ;
    }
  gx->geneDict = dictHandleCreate (50000, gx->h) ;
  gx->chromDict = dictHandleCreate (50, gx->h) ;
  gx->geneIdDict = dictHandleCreate (50000, gx->h) ;
  gx->affyDict = dictHandleCreate (50000, gx->h) ;
  gx->affy2gene = keySetHandleCreate (gx->h) ;
  gx->signatureGenes = arrayHandleCreate (64, KEYSET, gx->h) ;
  gx->signatureGenesSorted = arrayHandleCreate (64, KEYSET, gx->h) ;
  gx->genePlus = arrayHandleCreate (64, KEYSET, gx->h) ;
  gx->compares = arrayHandleCreate (64, COMPARE, gx->h) ;
  gx->compareDict = dictHandleCreate (64, gx->h) ;

  if (1)
    {  /* initialise the capture dict in the desired order */
      char *cq, *cp = getenv ("CAPTURES") ;
      while (cp)
	{
	  cq = strchr (cp, ' ') ;
	  if (cq) *cq = 0 ;
	  if (strlen (cp))
	    dictAdd (gx->captureDict, cp, 0) ;
	  cp = cq ;
	  if (cp) cp++ ;
	}
    }
  
  pushText (gx->info, "_NA_") ;
 
  if (1)
    {
      RC *rc ;
      gx->maxGroupLevel = 1 ; /* always allow _sumOfAll */
      dictAdd (gx->runDict, "_SumOfAllReadsInProject", &gx->runAny) ;
      rc = arrayp (gx->runs, gx->runAny, RC) ;
      rc->addCounts = TRUE ;
      rc->isAny = TRUE ;
      rc->groupLevel = 1 ;
      rc->run = gx->runAny ;  
      rc->sample = rc->title = rc->runid = stackMark (gx->info) ;
      pushText (gx->info, "SumOfAllReadsInProject") ;
      if (0)
	{
	  pushText (gx->info, "_AnyRunILM") ;
	  dictAdd (gx->runDict, "_AnyRun", &gx->runAnyILM) ;
	  rc = arrayp (gx->runs, gx->runAnyILM, RC) ;
	  rc->addCounts = TRUE ;
	  rc->isAny = TRUE ;
	  rc->run = gx->runAny ; 
	  rc->private = 0 ;
	  rc->sample = rc->title = rc->runid = stackMark (gx->info) ;
	  
	  pushText (gx->info, "_AnyRunSolid") ;
	  dictAdd (gx->runDict, "_AnyRun", &gx->runAnySolid) ;
	  rc = arrayp (gx->runs, gx->runAny, RC) ;
	  rc->addCounts = TRUE ;
	  rc->isAny = TRUE ;
	  rc->run = gx->runAny ;  
	  rc->sample = rc->title = rc->runid = stackMark (gx->info) ;
	  pushText (gx->info, "_AnyRun") ;
	}
    }

  if (1)
    {
      int i ;
      double z = 0 ;
      double *gauss5 = gx->gauss5 ;
      
      gauss5[0] = 1 ;
      for (i = 1 ; i <= 10 ; i++)
	gauss5[i] = exp (- i*i/50.0) ;   /* 2 sigma^2 = 18, sigma = 3 (sqrt_two-fold) */
      for (z = 1, i = 1 ; i <= 10 ; i++)
	z += 2 * gauss5[i] ;
      for (i = 0 ; i <= 10 ; i++)
	gauss5[i] /= z ;
      if (0)
	for (i = 0 ; i <= 10 ; i++)
	  fprintf (stderr, "Gauss\t%d\t%.2f\n", i, gauss5[i]) ;
    }
  if (1)
    {
      int i ;
      double z = 0 ;
      double *gauss10 = gx->gauss10 ;
      
      gauss10[0] = 1 ;
      for (i = 1 ; i <= 10 ; i++)
	gauss10[i] = exp (- i*i/50.0) ;   /* 2 sigma^2 = 50, sigma = 5 (sqrt_two-fold) */
      for (z = 1, i = 1 ; i <= 10 ; i++)
	z += 2 * gauss10[i] ;
      for (i = 0 ; i <= 10 ; i++)
	gauss10[i] /= z ;
      if (0)
	for (i = 0 ; i <= 10 ; i++)
	  fprintf (stderr, "Gauss\t%d\t%.2f\n", i, gauss10[i]) ;
    }
  if (1)
    {
      int i ;
      double z = 0 ;
      double *gauss20 = gx->gauss20 ;
      
      gauss20[0] = 1 ;
      for (i = 1 ; i <= 20 ; i++)
	gauss20[i] = exp (- i*i/200.0) ;  /* 2 sigma^2 = 200, sigma = 10 (twofold) */
      for (z = 1, i = 1 ; i <= 20 ; i++)
	z += 2 * gauss20[i] ;
      for (i = 0 ; i <= 20 ; i++)
	gauss20[i] /= z ;
      if (0)
	for (i = 0 ; i <= 20 ; i++)
	  fprintf (stderr, "Gauss\t%d\t%.2f\n", i, gauss20[i]) ;
    }
  if (1)
    {
      int i ;
      double z = 0 ;
      double *gauss30 = gx->gauss30 ;
      
      gauss30[0] = 1 ;
      for (i = 1 ; i <= 30 ; i++)
	gauss30[i] = exp (- i*i/450.0) ;  /* 2 sigma^2 = 450, sigma = 15 (about 3 fold) */
      for (z = 1, i = 1 ; i <= 30 ; i++)
	z += 2 * gauss30[i] ;
      for (i = 0 ; i <= 30 ; i++)
	gauss30[i] /= z ;
      if (0)
	for (i = 0 ; i <= 30 ; i++)
	  fprintf (stderr, "Gauss\t%d\t%.2f\n", i, gauss30[i]) ;
    }
  if (1)
    {
      int i ;
      double z = 0 ;
      double *gauss40 = gx->gauss40 ;
      
      gauss40[0] = 1 ;
      for (i = 1 ; i <= 40 ; i++)
	gauss40[i] = exp (- i*i/800.0) ;  /* 2 sigma^2 = 800, sigma = 20 (fourfold) */
      for (z = 1, i = 1 ; i <= 40 ; i++)
	z += 2 * gauss40[i] ;
      for (i = 0 ; i <= 40 ; i++)
	gauss40[i] /= z ;
      if (0)
	for (i = 0 ; i <= 40 ; i++)
	  fprintf (stderr, "Gauss\t%d\t%.2f\n", i, gauss40[i]) ;
    }
  hhShow (0, FALSE) ; /* for compiler happiness */

} /* gxParseInit */

/*************************************************************************************/
/* parse list of runs, to impose the order of the columns in the outputs */
static int gxParseRunList (GX *gx, const char* fileName)
{
  AC_HANDLE h = ac_new_handle () ; 
  ACEIN ai ; 
  const char *ccp ;

  ai = aceInCreate (fileName, gx->gzi, h) ;
  if (! ai)
    messcrash ("cannot open input file -runList %s", fileName ? fileName : "stdin") ;
  aceInSpecial (ai, "\"\n\t") ;

  while (aceInCard (ai)) 
    if ((ccp = aceInWord (ai)))
      if (*ccp != '/')
	dictAdd (gx->runDict, ccp, 0) ;
  ac_free (h) ;

  return dictMax (gx->runDict) ;
} /* gxParseRunList  */

/*************************************************************************************/
/* parse list of transcripts that should be excluded */
static int gxParseMask (GX *gx, const char* fileName)
{
  AC_HANDLE h = ac_new_handle () ; 
  ACEIN ai ; 
  char *cp, *cq ;

  ai = aceInCreate (fileName, FALSE, h) ;
  if (! ai)
    messcrash ("cannot open input file -mask  %s", fileName ? fileName : "stdin") ;
  aceInSpecial (ai, "\"\n\t") ;

  while (aceInCard (ai)) 
    if ((cp = aceInWord (ai)))
      {
	if (0)  /* we do not want to mask the whole gene by default, just the masked transcripts or maked gene if explicitly given */
	  { cq = strchr (cp, '.') ; if (cq) *cq = 0 ; }
	dictAdd (gx->maskDict, cp, 0) ;
      }
  ac_free (h) ;

  return dictMax (gx->maskDict) ;
} /* gxParseMask */

/*************************************************************************************/

static const char *NB_alias_double_runs_rename (const char *old)
{
  char *news[] = { "Rhs963", "Rhs967", "Rhs969", "Rhs971", "Rhs973", "Rhs975", "Rhs977", "Rhs979", "Rhs981", "Rhs1217", "Rhs1219", "Rhs1221", "Rhs1223", "Rhs1225", "Rhs1227", 0 } ;
  char *olds[] = { "Rhs964", "Rhs968", "Rhs970", "Rhs972", "Rhs974", "Rhs976", "Rhs978", "Rhs980", "Rhs982", "Rhs1218", "Rhs1220", "Rhs1222", "Rhs1224", "Rhs1226", "Rhs1228", 0 } ;

  int i ;

  for (i = 0 ; olds[i] ;  i++)
    if (! strcmp (old, olds[i]))
      return news[i] ;
  return old ;
} /* NB_alias_double_runs_rename */

/*************************************************************************************/
/* parse list of chrom alias => rename on the fly the IntMap names */
static int gxParseChromAlias (GX *gx, const char* fileName)
{
  AC_HANDLE h = ac_new_handle () ; 
  ACEIN ai ; 
  const char *ccp ;
  int chrom, chrom2 ;

  ai = aceInCreate (fileName, FALSE, h) ;
  if (! ai)
    messcrash ("cannot open input file -chromAlias  %s", fileName ? fileName : "stdin") ;
  aceInSpecial (ai, "\"\n\t") ;

  while (aceInCard (ai)) 
    if ((ccp = aceInWord (ai)))
      {
	if (*ccp == '#') continue ;
	dictAdd (gx->chromAliasDict, ccp, &chrom) ;
	if ((ccp = aceInWord (ai)))
	  dictAdd (gx->chromAliasDict, ccp, &chrom2) ;
	keySet (gx->chromAliasKs, chrom) = chrom2 ;
      }
  ac_free (h) ;

  return dictMax (gx->chromAliasDict) ;
} /* gxParseChromAlias  */

/*************************************************************************************/

static int gxSplitMrnaParse (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  int gene, mrna, x1, x2, gMix ;
  const char *ccp ;
  DICT *dict = 0 ;
  Array aa  = 0 ;
  ACEIN ai = aceInCreate (gx->splitMrnaFileName, 0, h) ;
  SPLITMRNA *up ;
  Associator ass ;
  GC* gc ;

  if (! ai)
    messcrash ("Sorry, i cannot find the -split_mRNAs file : %s", gx->splitMrnaFileName) ;
  aa = gx->splitMrnaArray = arrayHandleCreate (10000, SPLITMRNA, gx->h) ;
  dict = gx->splitMrnaDict = dictHandleCreate (100, gx->h) ;
  ass = gx->splitMrnaAss = assHandleCreate (gx->h) ;

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#')
	continue ;
      dictAdd (gx->geneDict, ccp, &gene) ;

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp)
	continue ;
      dictAdd (dict, ccp, &mrna) ;
      
      x1 = x2 = -1 ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &x1) ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &x2) ;

      if (x2 < 0 || x1 < 0)
	messcrash ("FATAL ERROR: In -split_mRNAs, missing coordinates in mRNA %s, line %d, file %s"
		   , dictName (dict, mrna) 
		   , aceInStreamLine (ai) 
		   , gx->splitMrnaFileName
		   ) ;
      
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp)
	continue ;

      x1 = x2 = 0 ; /* the coordinates where taken care of in bestali */

      aceInStep (ai, '\t') ; /* jump col 5 */
      ccp = aceInWord (ai) ;
      if (! ccp)
	continue ;
      gMix = gene ;
      if (strcmp (ccp, dictName (gx->geneDict, gene)))
	{
	  aceInStep (ai, '\t') ;
	  ccp = aceInWord (ai) ;
	  if (! ccp)
	    continue ;
	  dictAdd (gx->geneDict, ccp, &gMix) ; /* name(oldNam) */
	  if (0) 
	    dictAdd (gx->geneDict, messprintf ("%s(%s)", dictName (gx->geneDict, gene), ccp), &gMix) ;
	}
      
      assMultipleInsert (ass, assVoid(mrna), assVoid (nn)) ;
      up = arrayp (aa, nn++, SPLITMRNA) ;
      up->gene = gene ;
      up->gMix = gMix ;
      if (x1) 
	{
	  up->x1 = x1 ;
	  up->x2 = x2 ;
	}
      gc = arrayp (gx->genes, gene, GC) ; /* make room */
      gc->gMix = gMix ;
    }

  ac_free (h) ;
  return nn ;
} /* gxSplitMrnaParse */

/*************************************************************************************/
/* cleanup gene(gene) into gene */
static const char *gxSplitMrnaAlias (GX *gx, char *nam)
{
  char *cp = strchr (nam, '(') ;
  char *cq = strchr (nam, ')') ;

  if (strstr (nam, "AKT3_test"))
    return 0 ;
  if (cp && cq)
    {
      *cp = *cq = 0 ;
      if (strcasecmp (nam, cp+1))
	{ *cp = '(' ; *cq = ')' ; } /* restore */
      /* else we are happy to remove the redundant gene(gene) */

      return nam ;
    }
  else
    {
      GC *gc = 0 ;
      int gene = 0 ;

      if (dictFind (gx->geneDict , nam, &gene)
	  && gene < arrayMax (gx->genes)
	  && (gc = arrp (gx->genes, gene, GC))
	  && gc->gMix
	  )
	return dictName (gx->geneDict, gc->gMix) ;
    }
  return nam ;
} /* gxSplitMrnaAlias */

/*************************************************************************************/
/* parse the ace file run->deep and run->Length */
static int gxAceParse (GX *gx, const char* fileName,BOOL metaData)
{
  AC_HANDLE h = ac_new_handle () ; 
  ACEIN ai ; 
  Array aa, daa ;
  int nn = 0, gene = 0, geneId = 0, run = 0, affy = 0, compare = 0, nLine = 0, cover, count, wild ;
  int geneOld = 0, geneTrue = 0 ;
  float index = 0, seqs, tags, kb, nReads, nReadsOk, nerr, a2g, partial, orphan, badTopo, multi, multi2 ;
  const char *ccp ;
  DC *dc = 0 ;  
  DDC *ddc = 0 ;
  RC *rc = 0 ;
  GC *gc = 0 ;
  BOOL isOut, isRun, isAli, isGene, isGeneId, isIntron, isCompare, badSnp = FALSE ;
  /*   Stack buf = stackHandleCreate (1000,h) ; */
  DICT *markerSelectDict = 0 ;
  DICT *markerRejectDict = 0 ;
  BOOL wantDaa = gx->export && !gx->subsample && strchr (gx->export, 'v') ? TRUE : FALSE ;

  if (gx->markerSelectFileName)
    {
      ai = aceInCreate (gx->markerSelectFileName, gx->gzi, h) ;
      markerSelectDict = dictHandleCreate (10000, h) ;
      while (aceInCard (ai)) 
	{ 
	  ccp = aceInWord (ai) ;
	  if (ccp) dictAdd (markerSelectDict, ccp, 0) ;
	}

      ac_free (ai) ;
    }
  if (gx->markerRejectFileName)
    {
      ai = aceInCreate (gx->markerRejectFileName, gx->gzi, h) ;
      markerRejectDict = dictHandleCreate (10000, h) ;
      while (aceInCard (ai)) 
	{ 
	  ccp = aceInWord (ai) ;
	  if (ccp) dictAdd (markerRejectDict, ccp, 0) ;
	}

      ac_free (ai) ;
    }


  ai = aceInCreate (fileName, gx->gzi, h) ;
  if (! ai)
    messcrash ("cannot open input file -deep %s", fileName ? fileName : "stdin") ;
  aceInSpecial (ai, "\"\n\t") ;

  isOut = TRUE ;
  isRun = isGene = isIntron = isGeneId = isAli = isCompare = FALSE ; run = gene = 0 ;
  while (aceInCard (ai)) 
    { /* parse the ace file */
      if (! metaData && gx->isSNP)
	{
	  ccp = aceInPos (ai) ;
	  if (! ccp || strstr (ccp, "Incompatible_strands") || strstr (ccp, "N_rich") )
	    badSnp = TRUE ;
	}

      ccp = aceInWord (ai) ; 
      if (! (++nLine % 10000000))
	fprintf (stderr, "%d lines %d accepted records\n", nLine, dictMax(gx->geneDict)) ;
      if (ccp && ccp[0] == '/' && ccp[1] == '/') 
	continue ;
      if (! metaData && gx->isSNP)
	{
	  /*  new format, september 2023 : parse a tsf file exported from TSNP_DB script/snp3.tcsh */
	  float f = -10 ;
	  if (! ccp || *ccp == '#')
	    continue ;
	  dictAdd (gx->geneDict, ccp, &gene) ;
	  gc = arrayp (gx->genes, gene, GC) ; /* make room */
	  isGene = TRUE ; isOut = FALSE ;
	  
	  ccp = aceInWord (ai) ;
	  run = 0 ; isRun = FALSE ; isOut = FALSE ;
	  if (! gx->runListFileName)
	    dictAdd (gx->runDict, ccp, &run) ;
	  else
	    dictFind (gx->runDict, ccp, &run) ;
	  if (run > 0)
	    isRun = TRUE ;
	  else
	    continue ;

	  ccp = aceInWord (ai) ; /* format, trash */

	  if (! aceInFloat (ai, &f))
	    continue ;
	  cover = count = 0 ;
	  aceInInt (ai, &cover) ;
	  aceInInt (ai, &count) ;
	  aceInInt (ai, &wild) ;
          if (cover < count + wild) 
	    cover = count + wild ;
          
	  if (1)  /* dromadaire 2016_11_24 */
	    {
	      if (f < .4) 
		{ f = 0 ; wild += count ; count = 0 ; }
	      if (f > .9) 
		{ f = 1 ; count += wild ; wild = 0 ; }
	    }
	  if (cover >= 20)
	    {
	      rc = arrayp (gx->runs, run, RC) ;
	      rc->isSNP = TRUE ;
	      rc->hasData = TRUE ;
	      aa = rc->aa ;
	      if (! aa)
		aa = rc->aa = arrayHandleCreate (50000, DC, gx->h) ;
	      
	      /* register */
	      
	      if (gx->snpEval || gx->snpCompare )
		{  /* evaluate if the score of this snp, export and reinitialize */
		  int runMax = arrayMax (gx->runs) ;

		  geneTrue = gene ;
		  if (geneOld && geneTrue != geneOld)
		    {
		      if (nn && ! badSnp)
			{
			  int groupLevel ;
			  for (groupLevel = -1 ; groupLevel <= gx->maxGroupLevel ; groupLevel++)
			    for (run = 1, rc = arrayp (gx->runs, run, RC) ; run < runMax - 1 ; rc++, run++)
			      if ((! rc->runs || rc->addCounts) && rc->groupLevel == groupLevel)
				gxOneGroupCumul (gx, rc, 1) ;
			  gxSnpEvaluate (gx,1, dictName(gx->geneDict, geneOld)) ;
			}
		      /* reinitialize the aa arrays */
		      for (run = 1, rc = arrayp (gx->runs, run, RC) ; run < runMax - 1 ; rc++, run++)
			{
			  aa = rc->aa ;
			  rc->cumulDone = FALSE ;
			  if (! aa)
			    aa = rc->aa = arrayHandleCreate (5, DC, gx->h) ;
			  else if(arrayMax (aa))
			    aa = rc->aa = arrayReCreate (aa, 5, DC) ;
			}
		      badSnp = FALSE ;
		    }
		  geneOld = geneTrue ;
		  gene = 1 ;
		}
	      nn++ ;

	      if (1)
		{
		  /* register so that gxSnpEvaluate has something to evaluate */
		  dc = arrayp (aa, gene, DC) ;
		  dc->seqs = count ;
		  dc->tags = cover ;
		  dc->kb = count ; /* 2016_02_10 was genotype */
		  index = 0 + 100.0 * ((float)count)/((float)cover +.0001) ; 
		  dc->indexMin = dc->indexMax = dc->index = index + 1000 ;
		  gc->tags++ ; /* to insure the SNP is exported */
		  gc->isGood = TRUE ;
		}
	    }
	  continue ;
	}
      if (! ccp || ! *ccp)
	{
	  isOut = TRUE ;
	  isRun = isGene = isGeneId = isIntron = isAli = isCompare = FALSE ; run = gene = compare = 0 ;
	  continue ;
	}
      if (isOut && ! isRun && ! isGene && ! isGeneId && ! isAli && !strcasecmp (ccp, "Compare"))
	{
	  ccp = aceInWord (ai) ;
	  if (! ccp)
	    continue ;
	  dictAdd (gx->compareDict, ccp, &compare) ;
	  isCompare = TRUE ; isOut = FALSE ;
	  continue ;
	}
      if (isCompare)
	{
	  COMPARE *up = arrayp (gx->compares, compare, COMPARE) ;
	  up->compare = compare ;
	  if (!strcasecmp (ccp, "Runs"))
	    {
	      int run1 = 0 ;
	      
	      ccp = aceInWord (ai) ;
	     
	      run1 = 0 ;
	      if (! gx->runListFileName)
		dictAdd (gx->runDict, ccp, &run1) ;
	      else
		dictFind (gx->runDict, ccp, &run1) ;
	      if (run1 > 0)
		{
		  int i ;
		  KEYSET ks = up->runs ;
	      
		  if (! ks)
		    ks = up->runs = keySetHandleCreate (gx->h) ;

		  i = keySetMax (ks) ;
		  aceInInt (ai, &i) ;
		  keySet (ks, i) = run1 ; /* ordered list of runs */
		  if (aceInInt (ai, &i))
		    {
		      if (! up->ratios)
			up->ratios = keySetHandleCreate (gx->h) ;
		      keySet (up->ratios, run1) = i ; /* a priori ratio */
		    }
		}

	      continue ;
	    }
	  if (!strcasecmp (ccp, "Mixing_ratios"))
	    {
	      gx->mixing_ratios = up->mixing_ratios = TRUE ;
	      continue ;
	    }
	  if (!strcasecmp (ccp, "Correlation") || !strcasecmp (ccp, "Correlation_analysis"))
	    {
	      up->correlation = TRUE ;
	      continue ;
	    }
	  if (!strcasecmp (ccp, "Force_threshold"))
	    {
	      int th = 0 ;
	      if (aceInInt (ai, &th) && th > 0 && th < 200)
		up->forceIbest = th ;
	      continue ;
	    }
	  if (!strcasecmp (ccp, "Covariance") || !strcasecmp (ccp, "Covariance_analysis") ||
	      !strcasecmp (ccp, "Correlation") || !strcasecmp (ccp, "Correlation_analysis"))
	    {
	      up->correlation = TRUE ;
	      continue ;
	    }
	  if (!strcasecmp (ccp, "Who_is_who"))
	    { 
	      up->who_is_who = TRUE ;
	      gx->who_is_who = TRUE ;
	      gx->samplePairing = up->samplePairing = TRUE ;
	      if ( gx->isSNP)
		{
		  up->compareSNP = TRUE ;
		}
	      continue ;
	    }
	  if (!strcasecmp (ccp, "SNP"))
	    {
	      up->compareSNP = TRUE ;
	      fprintf (stderr, "## compare %s compareSNP\n", dictName(gx->compareDict,up->compare)) ;
	      continue ;
	    }
	  if (!strcasecmp (ccp, "INTRON"))
	    {
	      up->compareINTRON = TRUE ;
	      fprintf (stderr, "## compare %s compareINTRON\n", dictName(gx->compareDict,up->compare)) ;
	      continue ;
	    }
	  if (!strcasecmp (ccp, "Profile"))
	    {
	      gx->expressionProfile = up->expressionProfile = TRUE ;
	      continue ;
	    }
	  if (!strcasecmp (ccp, "Show_all_histos"))
	    {
	      gx->showAllHistos = up->showAllHistos = TRUE ;
	      continue ;
	    }
	  if (!strcasecmp (ccp, "Classify_by_homozygous_SNPs"))
	    {
	      gx->snpProfile = up->snpProfile = TRUE ;
	      continue ;
	    }
	  if (!strcasecmp (ccp, "Autosort"))
	    {
	      up->correlation = TRUE ;
	      up->autosort = TRUE ;
	      continue ;
	    }
	  if (!strcasecmp (ccp, "Sort_by_sorting_titles"))
	    {
	      up->correlation = TRUE ;
	      up->autosort = FALSE ;
	      continue ;
	    }
	  if (!strcasecmp (ccp, "Show_all_runs"))
	    {
	      up->correlation = TRUE ;
	      up->show_all_runs = TRUE ;
	      continue ;
	    }
	  if (!strcasecmp (ccp, "Predicts"))
	    {
	      int i = 0 ;

	      ccp = aceInWord (ai) ;
	      if (! ccp)
		continue ;
	      dictAdd (gx->compareDict, ccp, &i) ;
	      if (! up->predicts)
		up->predicts = keySetHandleCreate (gx->h) ;
	      keySet (up->predicts, keySetMax(up->predicts)) = i ;
	      array (gx->compares, i, COMPARE).compare = i ;
	      if (! gx->aoAUC2)
		gx->aoAUC2 = aceOutCreate (gx->outFileName, ".predictions.txt", 0, gx->h) ;
	      continue ;
	    }
	  if (0 &&                           /* 2019_05_05 remove Titration and shape : used in SEQC to list titrating genes */
	      !strcasecmp (ccp, "Shape"))
	    {
	      int i = 0 ;
	      float x ;
	      Array s = up->shape ;
	      BOOL ok1 = TRUE ;

	      if (! s)
		{
		  s = up->shape = arrayHandleCreate (8, float, gx->h) ;
		}

	      ccp = aceInPos (ai) ;
	      if (!ccp) 
		continue ;
	      if (
		  (ccp[0] >= '0' && ccp[0] <= '9') || 
		  (ccp[1] >= '0' && ccp[1] <= '9') || 
		  (ccp[1] && ccp[2] >= '0' && ccp[2] <= '9')
		  )
		{
		  while (aceInFloat (ai, &x))
		    array (s, i++, float) = x ;  
		}
	      else
		{
		  ccp-- ;
		  while (*++ccp)
		    {
		      switch (*ccp)
			{
			case '+': x = .1 ; break ;
			case '-': x = - .1 ; break ;
			case '=': x = 0 ; break ;
			default: ok1 = FALSE ; break ;
			}
		      array (s, i++, float) = x ;  
		    }
		}
	      
	      if (ok1 && i > 1)
		gx->hasTitration = TRUE ;
	      else
		{
		  arrayDestroy (up->shape) ;
		}
	    }

	  continue ;
	}
      if ( isOut &&
	   (
	    (gx->isTranscript && !strcasecmp (ccp, "Transcript")) || 
	    (gx->isTranscript && !strcasecmp (ccp, "mRNA")) || 
	    (!gx->isINTRON && !gx->isTranscript && !strcasecmp (ccp, "Gene")) || 
	    (gx->isINTRON && !strcasecmp (ccp, "Intron")) ||
	    (gx->isINTRON && !strcasecmp (ccp, "Gene"))
	    )
	   )
	{
	  if (!strcasecmp (ccp, "Gene")) 
	    {
	      char *cp = aceInWord (ai) ;
	      ccp = gxSplitMrnaAlias (gx, cp) ;
	    }
	  else
	    ccp = aceInWord (ai) ;
	  if (ccp && strncmp (ccp, "G_Any_", 6))
	    {
	      char *cr ;
	      if (0 && ! strcmp (ccp, "MIF(SLC2A11andMIF)"))
		invokeDebugger() ;
	      if (gx->maskDict &&  dictFind (gx->maskDict, ccp, 0))
		continue ;
	     if (gx->keepIndex && gx->userGivenGenePlusMinusDict && ! dictFind (gx->userGivenGenePlusMinusDict, ccp, 0))
	       continue ;

	      if (markerSelectDict && ! dictFind(markerSelectDict, ccp, 0))
		continue ;
	      if (markerRejectDict && dictFind(markerRejectDict, ccp, 0))
		continue ;

	     if (dictFind (gx->affyDict, ccp, &affy))
		gene = keySet (gx->affy2gene, affy) ;
	      if (! gene)
		{
		  dictAdd (gx->geneDict, ccp, &gene) ;
		}
	      gc = arrayp (gx->genes, gene, GC) ; /* make room */
	      if (gc->gMix && gc->gMix != gene)
		{
		  gene = gc->gMix ;
		  gc = arrayp (gx->genes, gene, GC) ; /* make room */
		}
	      isGene = TRUE ; isOut = FALSE ;

	      if (gx->isMRNAH && gx->isTranscript && ! strstr(ccp, ".a"))
		gc->notIsTranscriptA = TRUE ;

	      if (strstr(ccp,":"))
		gc->isSNP = TRUE ;
	      if ((gx->isINTRON || gx->intronGeneMix) &&
		  (cr = strstr (ccp, "__")) && strchr (cr+3,'_')
		  )
		{ 
		  isIntron =  TRUE ; gc->isIntron = 1 ; 
		}
	    }
	  continue ;
	}

      /* put this clause at the top to accelerate the code since Run_U or nU is most frequent line of input */
      if (ccp && isGene && gx->isINTRON && 
	  !strcasecmp (ccp, "de_uno")
	  )
	{
	  ccp = aceInWord (ai) ;
	  if (! ccp) continue ;

	  run = 0 ;
	  if (! gx->runListFileName)
	    dictAdd (gx->runDict, ccp, &run) ;
	  else
	    dictFind (gx->runDict, ccp, &run) ;
	  if (run <= 0)
	    continue ;

	  rc = arrayp (gx->runs, run, RC) ;
	  rc->cumulDone = TRUE ;
	  aa = rc->aa ;
	  daa = rc->daa ;
	  if (! aa)
	    aa = rc->aa = arrayHandleCreate (50000, DC, gx->h) ;
	  index = seqs = tags = kb = nReads = nReadsOk = nerr = a2g = orphan = partial = badTopo = multi = multi2 = 0 ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &tags) ;
	  nn++ ;
	  dc = arrayp (aa, gene, DC) ;
	  ddc = 0 ;
	  dc->seqs += tags ;
	  dc->tags += tags ;
	  gc->tags += tags ;
	  rc->geneTags += tags ;
	  if (gc->capturedIntronGene)
	    rc->captured_tags += tags ;
	  gc->isGood = TRUE ;
	  if (tags && ! kb) kb = tags/10.0 ; /* assume length = 100bp */
	  dc->kb += tags ;
	  

	  if (gc->intronGene)
	    { /* also register the introns counts of that gene */
	      dc = arrayp (aa, gc->intronGene, DC) ;
	      ddc = 0 ;
	      if (tags > dc->tags)
		{
		  dc->seqs = tags ;
		  dc->tags = tags ;
		  gc->isGood = TRUE ;
		  if (tags && ! kb) kb = tags/10.0 ; /* assume length = 100bp */
		}
	      dc->kb += tags ;
	    }
	  if (tags) rc->hasData = TRUE ;
	  continue ;
	}
      /* put this clause at the top to accelerate the code since Run_U or nU is most frequent line of input */
      if (ccp && isGene && 
	  (
	   (gx->unique && !strcasecmp (ccp, "Run_U")) || 
	   (!gx->unique && !strcasecmp (ccp, "Run_nU")) ||
	   (gx->anti && !strcasecmp (ccp, "Anti_run")) 
	   )
	  )
	{
	  BOOL anti = FALSE ;
	  if (!strcasecmp (ccp, "Anti_run")) 
	    anti = TRUE ;
	  /* get the count for this gene/run */
	  ccp = aceInWord (ai) ;
	  if (! ccp) continue ;

	  if (0) ccp = NB_alias_double_runs_rename (ccp) ;
	  run = 0 ;
	  if (! gx->runListFileName)
	    dictAdd (gx->runDict, ccp, &run) ;
	  else
	    dictFind (gx->runDict, ccp, &run) ;
	  if (run <= 0)
	    continue ;

	  rc = arrayp (gx->runs, run, RC) ;
	  if (anti)
	    {
	      aa = rc->antiAa ;
	      daa = 0 ;
	    }
	  else
	    {
	      aa = rc->aa ;
	      daa = rc->daa ;
	    }
	  if (! aa)
	    {
	      if (anti)
		{
		  aa = rc->antiAa = arrayHandleCreate (50000, DC, gx->h) ;
		}
	      else
		{
		  aa = rc->aa = arrayHandleCreate (50000, DC, gx->h) ;
		  if (wantDaa && ! gx->isTranscript && ! gx->isSNP && ! gx->isINTRON && gx->export && strchr (gx->export, 'a') )
		    {
		      daa = rc->daa = arrayHandleCreate (50000, DDC, gx->h) ;
		    }
		}
	    }

	  /* get the counts */
	  index = seqs = tags = kb = nReads = nReadsOk = nerr = a2g = orphan = partial = badTopo = multi = multi2 = 0 ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &index) ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &seqs) ;
	  aceInStep (ai,'\t') ; ccp = aceInWord(ai) ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &tags) ;
	  aceInStep (ai,'\t') ; ccp = aceInWord(ai) ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &kb) ;
	  aceInStep (ai,'\t') ; ccp = aceInWord(ai) ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &nReads) ;
	  aceInStep (ai,'\t') ; ccp = aceInWord(ai) ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &nReadsOk) ;
	  aceInStep (ai,'\t') ; ccp = aceInWord(ai) ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &nerr) ;
	  aceInStep (ai,'\t') ; ccp = aceInWord(ai) ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &a2g) ;
	  aceInStep (ai,'\t') ; ccp = aceInWord(ai) ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &partial) ;
	  aceInStep (ai,'\t') ; ccp = aceInWord(ai) ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &orphan) ;
	  aceInStep (ai,'\t') ; ccp = aceInWord(ai) ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &badTopo) ;
	  aceInStep (ai,'\t') ; ccp = aceInWord(ai) ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &multi) ;
	  aceInStep (ai,'\t') ; ccp = aceInWord(ai) ;
	  aceInStep (ai,'\t') ; aceInFloat (ai, &multi2) ;

	  nn++ ;
	  dc = arrayp (aa, gene, DC) ;
	  ddc = daa ? arrayp (daa, gene, DDC) : 0 ;
	  if (gx->keepIndex)
	    {
	      dc->index = 1000 + index ;
	      if (ccp && ! strcmp (ccp, "NA"))
		dc->isLow = 1 ;
	    }

	  if (0 && ! strncmp(dictName (gx->geneDict, gene), "X__A", 4))
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      char *sample = rc->sample ? stackText (gx->info, rc->sample) : 0 ; 
	      
	      if (sample && strstr (sample, "Fatigue"))
		{
		  seqs *= 2 ; tags *= 2 ; kb *= 2 ;
		}
	    }
	  if (0 && ! strncmp(dictName (gx->geneDict, gene), "X__B", 4))
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      char *sample = rc->sample ? stackText (gx->info, rc->sample) : 0 ; 
	      
	      if (sample && ! strstr (sample, "Fatigue"))
		{
		  seqs *= 2 ; tags *= 2 ; kb *= 2 ;
		}
	    }

	  if (tags && ! kb) kb = tags/50.0 ; /* assume length = 20bp */
	  dc->seqs += seqs ;
	  dc->tags += tags ;
	  gc->tags += tags ;
	  rc->geneTags += tags ;
	  
	  rc->nReads += nReads ;
	  rc->nReadsOk += nReadsOk ;
	  rc->kb += kb ;
	  rc->nerr += nerr ;
	  rc->a2g += a2g ;
	  rc->partial += partial ;
	  rc->orphan += orphan ;
	  rc->badTopo += badTopo ;
	  rc->multi += multi ;
	  rc->multi2 += multi2 ;

	  if (rc->capture && gc->captures)
	    {
	      KEYSET ks = gc->captures ;
	      int i ;
	      for (i = 0 ; i < keySetMax (gc->captures) ; i++)
		if (rc->capture == keySet (ks, i))
		  {
		    rc->captured_seqs += seqs ;
		    rc->captured_tags += tags ;
		    rc->captured_kb += kb ;


		    rc->CnReads += nReads ;
		    rc->CnReadsOk += nReadsOk ;
		    rc->Cnerr += nerr ;
		    rc->Ca2g += a2g ;
		    rc->Cpartial += partial ;
		    rc->Corphan += orphan ;
		    rc->CbadTopo += badTopo ;
		    rc->Cmulti += multi ;
		    rc->Cmulti2 += multi2 ;



		    break ;
		  }
	    }
	  gc->isGood = TRUE ;
	  dc->kb += kb ;

	  if (ddc)
	    {
	      if (nReads)
		gx->hasDaa = TRUE ;
	      ddc->nReads += nReads ;
	      ddc->nReadsOk += nReadsOk ;
	      ddc->nerr += nerr ;
	      ddc->a2g += a2g ;
	      ddc->partial += partial ;
	      ddc->orphan += orphan ;
	      ddc->badTopo += badTopo ;
	      ddc->multi += multi ;
	      ddc->multi2 += multi2 ;
	    }

	  if (tags) rc->hasData = TRUE ;
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Run"))
	{
	  ccp = aceInWord (ai) ;
	  isRun = FALSE ; isOut = FALSE ;
	  if (ccp) 
	    {
	      run = 0 ; isRun = FALSE ;
	      if (! gx->runListFileName)
		dictAdd (gx->runDict, ccp, &run) ;
	      else
		dictFind (gx->runDict, ccp, &run) ;
	      if (run > 0)
		isRun = TRUE ;
	    }
	  continue ;
	}
      if (ccp && isRun && !strcasecmp (ccp, "Avoid"))
	{
	  RC *rc = arrayp (gx->runs, run, RC) ;
	  rc->avoid = TRUE ;
	  continue ;
	}
      if (ccp && isOut && !isRun && !strcasecmp (ccp, "Ali"))
	{
	  ccp = aceInWord (ai) ;
	  isAli = FALSE ; isOut = FALSE ;
	  if (ccp) 
	    {
	      run = 0 ; isRun = FALSE ;
	      if (! gx->runListFileName)
		dictAdd (gx->runDict, ccp, &run) ;
	      else
		dictFind (gx->runDict, ccp, &run) ;
	      if (run > 0)
		isAli = TRUE ; 
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "SEQC_sigma"))
	{
	  float z = 0 ;
	  if (aceInFloat (ai, &z))
	    {
	      gc = arrayp (gx->genes, gene, GC) ;
	      if (1) gc->intrinsicSigma = z ;  /* dromadaire */
	      gc->isGood = TRUE ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Targeted"))
	{
	  if (isGene && gx->targeted)
	    {
	      gc = arrayp (gx->genes, gene, GC) ;
	      gc->targeted = TRUE ;
	      gc->isGood = TRUE ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Shift_expression_index"))
	{
	  float z = 0 ;
	  if (aceInFloat (ai, &z))
	    {
	      if (0 && isRun)
		{
		  rc = arrayp (gx->runs, run, RC) ;
		  rc->shift_expression_index = z ;
		}
	    }
	  continue ;
	}
      if (ccp && !gx->isINTRON && !strcasecmp (ccp, "Capture_touch"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp)
	    {
	      int capture = 0 ;
	      gx->hasCaptureTouch = TRUE ;

	      dictAdd (gx->captureDict, ccp, &capture) ;
	      if (isGene)
		{
		  KEYSET ks ;
		  gc = arrayp (gx->genes, gene, GC) ;
		  ks = gc->capturesTouch ;
		  if (!ks)
		    ks = gc->capturesTouch = keySetHandleCreate (gx->h) ;
		  keySetInsert (ks, capture) ;
		}
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Capture"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp)
	    {
	      int capture = 0 ;
	      gx->hasCapture = TRUE ;

	      if (!gx->targeted_genes)
		gx->targeted_genes = keySetHandleCreate (gx->h) ;
      	      dictAdd (gx->captureDict, ccp, &capture) ;
	      if (isGene)
		{
		  KEYSET ks ;
		  gc = arrayp (gx->genes, gene, GC) ;
		  if (gx->captured && ! strcmp (ccp,gx->captured))
		    { gc->captured = TRUE ; gc->isGood = TRUE ; }
		  ks = gc->captures ;
		  if (!ks)
		    ks = gc->captures = keySetHandleCreate (gx->h) ;
		  keySetInsert (ks, capture) ;
		  keySet (gx->targeted_genes, capture)++ ;
		}
	      else if (isRun)
		{
		  rc = arrayp (gx->runs, run, RC) ;
		  rc->capture = capture ;
		}
	    }
	  continue ;
	}
     if (ccp && !strcasecmp (ccp, "GeneId"))
	{
	  int old ;
	  ccp = aceInWord (ai) ;
	  if (ccp && isGene)
	    {
	      gc = arrayp (gx->genes, gene, GC) ;
	      gc->isGood = TRUE ;
	      old = gc->geneId ;
	      if (old && strstr(stackText (gx->info,old),ccp)) continue ;
	      if (old) ccp = messprintf ("%s;%s",stackText (gx->info,old),ccp) ;
	      gc->geneId = stackMark (gx->info) ;
	      pushText (gx->info, ac_unprotect (ccp,h)) ;
	      gx->hasGeneId = TRUE ;
	    }
	  else if (isOut && ccp)
	    {
	      dictAdd (gx->geneIdDict, ccp, &geneId) ;
	      isGeneId = TRUE ; isOut = FALSE ;
	    }
	  continue ;
	}

     if (ccp  && isGene && gx->isINTRON && !strcasecmp (ccp, "Gene"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp)
	    {
	      gc = arrayp (gx->genes, gene, GC) ;
	      dictAdd (gx->geneDict, ccp, &(gc->intronGene)) ;
	      /*
		gc->geneModel = stackMark (gx->info) ;
		pushText (gx->info, ac_unprotect (ccp,h)) ;
	      */
	      gx->hasFromGene = TRUE ;
	    }
	  continue ;
	}

     if (ccp  && isGene && gx->isINTRON && !strcasecmp (ccp, "In_mRNA"))
	{
	  gc = arrayp (gx->genes, gene, GC) ;
	  gc->intronMrna++ ;
	  continue ;
	}

     if (ccp  && isGene && gx->isINTRON && !strcasecmp (ccp, "Known_donor"))
	{
	  ccp = aceInWord (ai) ;
	  gc = arrayp (gx->genes, gene, GC) ;
	  gc->hasKnownDonorAcceptor |= 1 ;
	  gx->hasKnownDonorAcceptor= TRUE ;
	  continue ;
	}

     if (ccp  && isGene && gx->isINTRON && !strcasecmp (ccp, "Known_acceptor"))
	{
	  gc = arrayp (gx->genes, gene, GC) ;
	  gc->hasKnownDonorAcceptor |= 2 ;
	  gx->hasKnownDonorAcceptor= TRUE ;
	  continue ;
	}

     if (ccp && gx->isTranscript  && !strcasecmp (ccp, "Model_of"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp)
	    {
	      if (!strncmp(ccp, "X__",3))ccp+=3 ;
	      gc = arrayp (gx->genes, gene, GC) ;
	      gc->geneModel = stackMark (gx->info) ;
	      pushText (gx->info, ac_unprotect (ccp,h)) ;
	      gx->hasGeneModel = TRUE ;
	    }
	  continue ;
	}

      if (! isAli && ! isGene && ! isGeneId && ! isRun)
	continue ;

      if (ccp && isGene && !strcasecmp (ccp, "Length"))
	{
	  int ln = -1 ;
	  if (aceInInt (ai, &ln) && ln > 0)
	    {
	      gc = arrayp (gx->genes, gene, GC) ;
	      gc->isGood = TRUE ;
	      gc->length = ln ;
	      gx->hasGeneLength = TRUE ;
	    }
	  continue ;
	}

      if (ccp && isIntron &&
	  (
	   !strcasecmp (ccp, "gt_ag") ||
	   !strcasecmp (ccp, "ct_ac") ||
	   !strcasecmp (ccp, "gc_ag") ||
	   (
	    !strcasecmp (ccp, "Other") && (ccp = aceInWord(ai))
	    )
	   )
	  )
	{
	  gc = arrayp (gx->genes, gene, GC) ;
	  if (gc->isIntron)
	    {
	      int t = 0 ;
	      
	      dictAdd (gx->affyDict, ccp, &t) ;
	      gc = arrayp (gx->genes, gene, GC) ;
	      gc->intronType = t ;
	      gx->hasIntronType = TRUE ;
	    }
	  continue ;
	}

      if (ccp && isIntron && 
	   !strcasecmp (ccp, "RefSeq") 
	  )
	{
	  gc = arrayp (gx->genes, gene, GC) ;
	  if (gc->isIntron)
	    {
	      gx->hasRefSeqAv = TRUE ;
	      aceInInt (ai, &(gc->intronRefSeq)) ;
	    }
	  continue ;
	}
      if (ccp && isIntron &&
	   !strcasecmp (ccp, "mRNA") 
	  )
	{
	  gc = arrayp (gx->genes, gene, GC) ;
	  if (gc->isIntron)
	    {
	      gx->hasRefSeqAv = TRUE ;
	      aceInInt (ai, &(gc->intronAv)) ;
	    } 
	  continue ;
	}

      if (ccp && isAli && !strcasecmp (ccp, "Intergenic_density"))
	{ 
	  float gkb ;
	  if (aceInFloat (ai, &gkb)) 
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->intergenicDensity = gkb ;
	    }
	  continue ;
	}

      if (ccp && isAli && !strcasecmp (ccp, "Intergenic"))
	{ 
	  float gkb ;
	  if (aceInFloat (ai, &gkb)) 
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->intergenicKb = gkb ;
	    }
	  continue ;
	}
      if (ccp && isAli && !strcasecmp (ccp, "Candidate_introns"))
	{ 
	  int i ;
	  char *cx = 0 ;
	  float known = 0, new = 0, minS = 0 ;
	  for (i = 0 ; i < 12 ; i++)
	    if (!aceInWord (ai))
	      break ;
	  if (i == 12 
 	      && aceInFloat (ai, &known) 
	      && (cx = aceInWord (ai))
	      && ! strcasecmp (cx, "New_support")
	      && aceInFloat (ai, &new) 
	      && (cx = aceInWord (ai))
	      && ! strcasecmp (cx, "New_minS")
	      && aceInFloat (ai, &minS) 
	      ) 
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->tags = known + new ;
	      rc->minIntronSupport = minS ;
	    }
	  continue ;
	}
      if (1 && ccp && isAli && !strcasecmp (ccp, "Accessible_length_5kb"))
	{ 
	  int aln = 0 ;
	  if (aceInInt (ai, &aln)) 
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->accessibleLength5kb = aln ;
	      gx->hasAccessibleLength = TRUE ;
	    }
	  continue ;
	}
      if (1 && ccp && isAli && !strcasecmp (ccp, "Accessible_length_8kb"))
	{ 
	  int aln = 0 ;
	  if (aceInInt (ai, &aln)) 
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->accessibleLength8kb = aln ;
	      gx->hasAccessibleLength = TRUE ;
	    }
	  continue ;
	}

      if (ccp && isAli && !strcasecmp (ccp, "Accepted"))
	{
	  float raw = 0 ;

	  ccp = aceInWord (ai) ; 
	  if (! ccp) continue ;
	  ccp = aceInWord (ai) ; 
	  if (! ccp) continue ;
	  if (aceInFloat (ai, &raw))
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->accepted_tags = raw ;
	    }
	  continue ;
	}

      if (ccp && isAli && !strcasecmp (ccp, "h_Ali"))
	{
	  int i,j ;
	  float gkb ;

	  ccp = aceInWord (ai) ;
	  if (! ccp || strcmp (ccp, "Z_genome")) continue ;
	  
	  for (i = j = 0 ; i < 4 ; i++)
	     if ((ccp = aceInWord (ai))) j++ ;
	  if (j < 4) continue ;
	  if (aceInFloat (ai, &gkb))
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->genomicKb = gkb ;
	    }
	  continue ;
	}

      if (ccp && isAli && !strcasecmp (ccp, "nh_Ali"))
	{
	  int i,j ;
	  float tags = 0, gkb = 0 ;
	  RC *rc = arrayp (gx->runs, run, RC) ;

	  ccp = aceInWord (ai) ;
	  if (! ccp) continue ;

	  if (gx->target_class && !strcmp (ccp, gx->target_class))
	    {	  
	      for (i = j = 0 ; i < 2 ; i++)
		if ((ccp = aceInWord (ai))) j++ ;
	      if (j < 2) continue ;
	      if (aceInFloat (ai, &tags))
		rc->tags = tags ;

	      for (i = j = 0 ; i < 1 ; i++)
		if ((ccp = aceInWord (ai))) j++ ;
	      if (j < 1) continue ;
	      if (aceInFloat (ai, &gkb))
		rc->targetKb = gkb ;
	    }
	  if (!strcasecmp (ccp, "Any"))
	    {	  
	      for (i = j = 0 ; i < 4 ; i++)
		if ((ccp = aceInWord (ai))) j++ ;
	      if (j < 4) continue ;
	      gkb = 0 ;
	      if (aceInFloat (ai, &gkb))
		rc->anyTargetKb = gkb ;
	    }
	  continue ;
	}
      
      if (ccp && !strncasecmp (ccp, "PolyA", 5))   /* 2014_06_26 was PolyA_selected_RNA */
	{
	  if (isRun)
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->PolyA_selected_RNA = TRUE ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Paired_end"))  
	{
	  if (ccp && isRun)
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->isPairedEnd = TRUE ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "DGE"))  
	{
	  if (ccp && isRun)
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->DGE = TRUE ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "polyA"))  
	{
	  if (ccp && isRun)
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->pA = 2 ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Total_RNA"))  
	{
	  if (ccp && isRun)
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->pA = 1 ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Variance"))  
	{
	  if (ccp && isRun)
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      gx->hasSelectedVariance = rc->selectedVariance = TRUE ;
	    }
	  continue ;
	}

      if (ccp && !strcasecmp (ccp, "Sample"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp && isRun)
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->sample = stackMark (gx->info) ;
	      pushText (gx->info, ac_unprotect (ccp,h)) ;
	      gx->hasRunSample = TRUE ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "RunId"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp && isRun && !strstr (ccp, "NULL") && ! (strlen(ccp) > 100))
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->runid = stackMark (gx->info) ;
	      pushText (gx->info, ac_unprotect (ccp, h)) ;
	      gx->hasRunId = TRUE ;
	    }
	  continue ;
	}
      if (ccp && (!strcasecmp (ccp, "Affy") || !strcasecmp (ccp, "MicroArray")))
	{
	  int old ;
	  ccp = aceInWord (ai) ;
	  if (ccp && isGene)
	    {
	      int affy ;
	      gc = arrayp (gx->genes, gene, GC) ;
	      old = gc->affy ;
	      if (old && strstr(stackText (gx->info,old),ccp)) continue ;
	      dictAdd (gx->affyDict, ccp, &affy) ;
	      if (old) ccp = messprintf ("%s;%s",stackText (gx->info,old),ccp) ;
	      gc->affy = stackMark (gx->info) ;
	      pushText (gx->info, ac_unprotect (ccp, h)) ;
	      gx->hasGeneAffy = TRUE ;
	      keySet (gx->affy2gene, affy) = gene ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "From_gene"))
	{  
	  int old ;
	  ccp = aceInWord (ai) ;
	  if (ccp && gene && (gx->isINTRON || gx->isTranscript || gx->isMA))
	    {
	      gc = arrayp (gx->genes, gene, GC) ;
	      gc->isGood = TRUE ;
              old = gc->fromGene ;
	      if (old && strstr(stackText (gx->info,old),ccp)) continue ;
	      if (old) ccp = messprintf ("%s;%s",stackText (gx->info,old),ccp) ;
	      if (gc->isIntron)
		{
		  dictAdd (gx->geneDict, ccp, &(gc->intronGene)) ;
		  gc->hasKnownDonorAcceptor |= 4 ;
		  gx->hasKnownDonorAcceptor= TRUE ;
		}
	      else
		{
		  gc->fromGene = stackMark (gx->info) ;
		  pushText (gx->info, ac_unprotect (ccp,h)) ;
		}
	      gx->hasFromGene = TRUE ;
	    } 
	  continue ;
	}

      if (ccp && !strcasecmp (ccp, "From_transcript"))
	{  
	  int old ;
	  ccp = aceInWord (ai) ;
	  if (ccp && gene && (gx->isINTRON || gx->isTranscript || gx->isMA))
	    {
	      gc = arrayp (gx->genes, gene, GC) ;
              old = gc->fromTranscript ;
	      if (old && strstr(stackText (gx->info,old),ccp)) continue ;
	      if (old) ccp = messprintf ("%s;%s",stackText (gx->info,old),ccp) ;
	      gc->fromTranscript = stackMark (gx->info) ;
	      pushText (gx->info, ac_unprotect (ccp,h)) ;
	      gx->hasFromTranscript = TRUE ;
	    } 
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "NM_id"))
	{
	  int old ;
	  ccp = aceInWord (ai) ;
	  if (ccp && isGene)
	    {
	      gc = arrayp (gx->genes, gene, GC) ;
	      old = gc->nmid ;
	      if (old && strstr(stackText (gx->info,old),ccp)) continue ;
	      if (old) ccp = messprintf ("%s;%s",stackText (gx->info,old),ccp) ;
	      gc->nmid = stackMark (gx->info) ;
	      pushText (gx->info, ac_unprotect (ccp,h)) ;
	      gx->hasGeneNm = TRUE ;
	    }
	  if (ccp && isGeneId)
	    {
	      gc = arrayp (gx->geneIds, geneId, GC) ;
	      old = gc->nmid ;
	      if (old && strstr(stackText (gx->info,old),ccp)) continue ;
	      if (old) ccp = messprintf ("%s;%s",stackText (gx->info,old),ccp) ;
	      gc->nmid = stackMark (gx->info) ;
	      pushText (gx->info, ac_unprotect (ccp,h)) ;
	      gx->hasGeneNm = TRUE ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Alias"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp && isGene)
	    {
	      int old ;
	      
	      gc = arrayp (gx->genes, gene, GC) ;
	      old = gc->alias ;
	      if (old && strstr(stackText (gx->info,old),ccp)) continue ;
	      if (old) ccp = messprintf ("%s;%s",stackText (gx->info,old),ccp) ;
	      gc->alias = stackMark (gx->info) ;
	      pushText (gx->info, ac_unprotect (ccp,h)) ;
	      gx->hasGeneAlias = TRUE ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Good_product"))
	{
	  if (isGene)
	    {
	      gc = arrayp (gx->genes, gene, GC) ;
	      gc->hasGoodProduct = gx->hasGoodProduct = TRUE ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "IntMap"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp && isGene)
	    {
	      int chrom ;

	      gc = arrayp (gx->genes, gene, GC) ;
	      gc->isGood = TRUE ;
	      if (gx->chromAliasDict && 
		  dictFind (gx->chromAliasDict, ccp, &chrom) &&
		  keySet (gx->chromAliasKs, chrom)
		  )
		ccp = dictName (gx->chromAliasDict, keySet (gx->chromAliasKs, chrom)) ;
	      dictAdd (gx->chromDict, ccp, &chrom) ;
	      gc->chrom = chrom ;
	      gx->hasGeneChrom = TRUE ;
	      gc->a1 = gc->a2 = 0 ; gc->strand = '+' ;
	      if (aceInInt (ai, &(gc->a1)))
		aceInInt (ai, &(gc->a2)) ;
	      if (gc->a1 > gc->a2)
		{ 
		  int a0 = gc->a1 ; 
		  gc->a1 = gc->a2 ; 
		  gc->a2 = a0 ; 
		  gc->strand = '-' ;
		}
	      if (gc->a1 < gc->a2)
		{
		  gc->boxLength = gc->a2 - gc->a1 + 1 ;
		  gx->hasGeneBoxLength = TRUE ; 
		  if (0) /* this was 1 between at least sept 2021 and may 27 2022, a bug, because we only had length and not yet boxLength */
		    { 
		      gc->length = gc->a2 - gc->a1 + 1 ;
		      gx->hasGeneLength = TRUE ; 
		    }
		}
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Group_level")) 
	{
	  if (gx->geneGroupFileName && isRun)
	    {
	      int x = 0 ;
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      if (aceInInt (ai, &x))
		rc->groupLevel = x ;
	      if (x > gx->maxGroupLevel)
		gx->maxGroupLevel = x ;
	    }
	  continue ;
	}

      if (ccp && (!strcasecmp (ccp, "Group") ||  !strcasecmp (ccp, "Sublibrary_of")))
	{
	  BOOL isSubLib = FALSE ;
	  if (!strcasecmp (ccp, "Sublibrary_of")) 
	    {
	      if (isRun)
		{
		  RC *rc = arrayp (gx->runs, run, RC) ;
		  isSubLib = TRUE ; 
		  rc->isSubLib = 1 ; /* do not export these groups */
		  rc->groupLevel = -1 ;
		  rc->private = 1 ;
		}
	    }
	  ccp = aceInWord (ai) ;
	  if (ccp && isRun)
	    {
	      RC *rc ;
	      int run2 = 0 ;

	      if (! gx->runListFileName)
		dictAdd (gx->runDict, ccp, &run2) ;
	      else
		dictFind (gx->runDict, ccp, &run2) ;
	      if (run2 > 0)
		{
		  rc = arrayp (gx->runs, run2, RC) ;
		  if (! rc->runs)
		    rc->runs = keySetHandleCreate (gx->h) ;
		  keySet (rc->runs, keySetMax (rc->runs)) = run ;
		  gx->hasGroup = TRUE ;
		  if (isSubLib || gx->snpEval)
		    rc->addCounts = TRUE ;
		  if (! rc->groupLevel && ! isSubLib)
		    rc->groupLevel = 1 ;
		} 
	    }
	  continue ;
	}

      if (ccp && !strcasecmp (ccp, "Private")) /* do not export these groups */
	{
	  if (isRun)
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->private = 2 ;
	    }
	  continue ;
	}

      if (ccp && !strcasecmp (ccp, "Title"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp)
	    {
	      if (isRun)
		{
		  RC *rc = arrayp (gx->runs, run, RC) ;
		  rc->title = stackMark (gx->info) ;
		  pushText (gx->info, ac_unprotect (ccp,h)) ;
		  gx->hasRunTitle = TRUE ;
		}
	      if (isGene)
		{
		  gc = arrayp (gx->genes, gene, GC) ;
		  gc->title = stackMark (gx->info) ;
		  gc->isGood = TRUE ;
		  pushText (gx->info, ac_unprotect (ccp,h)) ;  
		  gx->hasGeneTitle = TRUE ;
		}
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Gene_type"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp)
	    {
	      if (isGene)
		{
		  gc = arrayp (gx->genes, gene, GC) ;
		  gc->isGood = TRUE ;
		  gc->gene_type = stackMark (gx->info) ;
		  pushText (gx->info, ac_unprotect (ccp,h)) ;  
		  gx->hasGeneType = TRUE ;
		}
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Sorting_title"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp)
	    {
	      if (isRun)
		{
		  RC *rc = arrayp (gx->runs, run, RC) ;
		  rc->sortingTitle = stackMark (gx->info) ;
		  pushText (gx->info, ac_unprotect (ccp,h)) ;
		  gx->hasRunSortingTitle = TRUE ;
		}
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Sorting_title_2"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp)
	    {
	      if (isRun)
		{
		  RC *rc = arrayp (gx->runs, run, RC) ;
		  rc->sortingTitle2 = stackMark (gx->info) ;
		  pushText (gx->info, ac_unprotect (ccp,h)) ;
		  gx->hasRunSortingTitle2 = TRUE ;
		}
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Other_title"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp)
	    {
	      if (isRun)
		{
		  RC *rc = arrayp (gx->runs, run, RC) ;
		  rc->otherTitle = stackMark (gx->info) ;
		  pushText (gx->info, ac_unprotect (ccp,h)) ;
		  gx->hasRunOtherTitle = TRUE ;
		}
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "Machine"))
	{
	  ccp = aceInWord (ai) ;
	  if (ccp && isRun)
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->machine = stackMark (gx->info) ;
	      pushText (gx->info, ac_unprotect (ccp,h)) ;
	      gx->hasMachine = TRUE ;
	    }
	  continue ;
	}
      if (ccp && !strcasecmp (ccp, "SOLiD"))
	{
	  if (isRun)
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->solid = TRUE ;
	    }
	  continue ;
	}

      if (ccp && !strcasecmp (ccp, "Add_counts"))
	{
	   if (ccp && isRun)
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->addCounts = TRUE ;
	    }
	   continue ;
	}
      if (ccp && !strcasecmp (ccp, "Average_fragment_length"))
	{
	  int fLength = 0 ;
	  if (aceInInt (ai, &fLength) && fLength > 30 && isRun)
	    {
	      RC *rc = arrayp (gx->runs, run, RC) ;
	      rc->fragmentLength = fLength ;
	    }
	  continue ;
	}
      if (ccp && isGene && !strcasecmp (ccp, "Validated_U"))
	{
	  /* get the count for this gene/run */
	  ccp = aceInWord (ai) ;
	  if (! ccp) continue ;

	  continue ; /* 2014_02_02 : we never use ->val so i stop reading it */
#ifdef JUNK
	  run = 0 ;
	  if (! gx->runListFileName)
	    dictAdd (gx->runDict, ccp, &run) ;
	  else
	    dictFind (gx->runDict, ccp, &run) ;
	  if (run <= 0)
	    continue ;

	  rc = arrayp (gx->runs, run, RC) ;
	  aa = rc->aa ;
	  if (! aa)
	    aa = rc->aa = arrayHandleCreate (50000, DC, gx->h) ;

	  /* get the counts */
	  val = 0 ;
	  aceInStep (ai,'\t') ; aceInInt (ai, &val) ;

	  if (val)
	    {
	      dc = arrayp (aa, gene, DC) ;
	      dc->val = val ;
	    }
#endif
	  continue ;
	}
      if (ccp && isGene && (!strcasecmp (ccp, "Micro_array")))
	{
	  /* get the index for this gene/micro-array */
	  ccp = aceInWord (ai) ;
	  if (! ccp) continue ;

	  run = 0 ;
	  if (! gx->runListFileName)
	    dictAdd (gx->runDict, ccp, &run) ;
	  else
	    dictFind (gx->runDict, ccp, &run) ;
	  if (run <= 0)
	    continue ;

	  rc = arrayp (gx->runs, run, RC) ;
	  rc->isMA = gx->isMA ? gx->isMA : 20 ;
	  rc->hasData = TRUE ;
	  aa = rc->aa ;
	  if (! aa)
	    aa = rc->aa = arrayHandleCreate (50000, DC, gx->h) ;

	  /* get the counts */
	  aceInStep (ai,'\t') ; aceInFloat (ai, &index) ;

	  nn++ ;
	  dc = arrayp (aa, gene, DC) ;
	  dc->tags = 1 ;
	  dc->indexMin = dc->indexMax = dc->index = index + 1000 ;
	  gc->tags++ ; /* to insure the gene is exported */
	  gc->isGood = TRUE ;
	  continue ;
	}
    }
  fprintf (stderr, "// %d lines %d accepted records\n", nLine, dictMax(gx->geneDict)) ; 
  if (gx->endpoints)
    {
      int i ;

      if (0) /* switch in run order, this seems silly*/
	for (i = 0 ; i < arrayMax (gx->endpoints) ; i++)
	  {
	    ENDPOINT *ep = arrp (gx->endpoints, i, ENDPOINT) ;
	    if (ep->run1 > ep->run2) /* switch */
	      {
		run = ep->run1 ; ep->run1 = ep->run2 ; ep->run2 = run ;
		run = ep->run3 ; ep->run3 = ep->run4 ; ep->run4 = run ;
	    }
	  }
      arraySort (gx->endpoints, endPointOrder) ;
      arrayCompress (gx->endpoints) ;
     }

  if (gx->compares)
    {
      int i, j, iCompare = 0 ;
      COMPARE *compare, *compare2 ;
      ENDPOINT *endpoint ;

      gx->hasCompare = FALSE ;
      /* keep happy few */
      for (j = iCompare = 0 ; iCompare < arrayMax (gx->compares) ; iCompare++)
	{
	  compare = arrp(gx->compares, iCompare, COMPARE) ;
	  if (! compare->runs)
	    continue ;    
	  if (gx->isSNP && ! compare->compareSNP)
	    keySetDestroy (compare->runs) ;
	  if (gx->isINTRON && ! compare->compareINTRON)
	    keySetDestroy (compare->runs) ;
	}
      
      /* compress the list of runs keeping the order */
      for (iCompare = 0 ; iCompare < arrayMax (gx->compares) ; iCompare++)
	{
	  compare = arrp(gx->compares, iCompare, COMPARE) ;
	  if (! compare->runs)
	    continue ;

	  for (i = j = 0 ; i < keySetMax (compare->runs) ; i++)
	    {
	      run = keySet (compare->runs, i) ;
	      if (run) keySet (compare->runs, j++) = run ;
	    }
	  keySetMax (compare->runs) = j ;
	  if (j < 1 || (j < 2 && ! compare->snpProfile))
	    keySetDestroy (compare->runs) ;
	  else
	    gx->hasCompare = TRUE ; 
	}
      
      /* create the old fashioned enppoint structure */
      for (iCompare = 0 ; iCompare < arrayMax (gx->compares) ; iCompare++)
	{
	  compare = arrp(gx->compares, iCompare, COMPARE) ;
	  if (! compare->runs || keySetMax (compare->runs) != 2)
	    continue ;
	  if (! gx->endpoints)
	    gx->endpoints = arrayHandleCreate (16, ENDPOINT, gx->h) ;
	  endpoint = arrayp (gx->endpoints, arrayMax(gx->endpoints), ENDPOINT) ;
	  endpoint->run1 = keySet (compare->runs, 0) ;
	  endpoint->run2 = keySet (compare->runs, 1) ;

	  if (compare->predicts)
	    {
	      int i2 = 0, i3 = 0 ;
	      
	      for (i2 = 0 ; i2 < keySetMax (compare->predicts) ; i2++)
		{
		  compare2 = arrp(gx->compares, keySet(compare->predicts, i2), COMPARE) ;
		  if (compare2->runs && keySetMax (compare2->runs) == 2)
		    {
		      if (i3++)
			{
			  endpoint = arrayp (gx->endpoints, arrayMax(gx->endpoints), ENDPOINT) ;
			  *endpoint = *(endpoint - 1) ;
			}
		      endpoint->run3 = keySet (compare2->runs, 0) ;
		      endpoint->run4 = keySet (compare2->runs, 1) ;
		    }
		}
	    }
	}
    }


  if (gx->endpoints)
    {
      arraySort (gx->endpoints, endPointOrder) ;
      arrayCompress (gx->endpoints) ;
    }

  if ((gx->snpEval  || gx->snpCompare) && nn && ! badSnp && geneTrue)
    {	
      int runMax = gx->runs ? arrayMax (gx->runs) : 0 ;
      int groupLevel ;
      for (groupLevel = -1 ; groupLevel <= gx->maxGroupLevel ; groupLevel++)
	for (run = 1, rc = arrayp (gx->runs, run, RC) ; run < runMax - 1 ; rc++, run++)
	  if ((! rc->runs || rc->addCounts) && rc->groupLevel == groupLevel)
	    gxOneGroupCumul (gx, rc, 1) ;
      gxSnpEvaluate (gx,1, dictName(gx->geneDict, geneTrue)) ;
    }
  gx->hasGeneNm = TRUE ;

  if (metaData && gx->hasGroup) 
    { /* check that sublibraries do not contain runs */
      int runMax = arrayMax (gx->runs) ;

       for (run = 0, rc = arrayp (gx->runs, 0, RC) ; run < runMax ; rc++, run++)
	 {
	   if (rc && rc->run && rc->runs && rc->isSubLib)
	     messcrash ("Run %s is both a  group and a sublibrary_of"
			, dictName (gx->runDict, rc->run)) ;
	 }
       /* check that sublibraries not part of a group, other than the library */
       for (run = 0, rc = arrayp (gx->runs, 0, RC) ; run < runMax ; rc++, run++)
	 {
	   int i, j = 0, k = 0 ;
	   if (rc->run  && rc->runs)
	     {
	       KEYSET ks1 = rc->runs ;
	       keySetSort (ks1) ;
	       keySetCompress (ks1) ;

	       for (i = j = k = 0 ; i < keySetMax (ks1) ; i++)
		 {
		   int run2 = keySet (ks1, i) ;
		   RC *rc2 = arrayp (gx->runs, run2, RC)  ;
		   if (rc2->isSubLib) j++ ;
		   k++ ;
		 }
	       if (j && j != k)
		  messcrash ("Run %s contains both sublibrary_of and normal runs"
			     , dictName (gx->runDict, rc->run)) ;
	       if (j)
		 rc->addCounts = TRUE ;
	     }
	 }
    }
  
  if (metaData)
    {
      int run, runMax = arrayMax (gx->runs) ;
      for (run = 0, rc = arrayp (gx->runs, 0, RC) ; run < runMax ; rc++, run++)
	rc->run = run ;
    }
  if (metaData && gx->hasGroup) 
    { /* 2016_03_02 
       * if A is a non additive supergroup containing runs  r1, r2 and groups g1 g2
       * expand g1, g2 into more runs, and remove gi, unless gi is additive 
       * so that we do not compute the average of an average
       * and we do not double count the overlapping subgroups
       * but we directly average the leaf runs
       *
       * Note that additive groups never contain overlapping subgroups
       * and never contain non-additive subgroups
       */
      int run, runMax = arrayMax (gx->runs) ;
      int level, maxLevel = gx->maxGroupLevel ;
      RC *rc ;
      
      for (level = 1 ; level <= maxLevel ; level++)
 	{
	  for (run = 0, rc = arrayp (gx->runs, 0, RC) ; run < runMax ; rc++, run++)
	    {
	      if (rc->run && rc->runs && rc->groupLevel == level)
		{
		  KEYSET ks1 = rc->runs ;
		  KEYSET ks2 = keySetHandleCreate (gx->h) ;
		  RC *rc2 ;
		  int i, j, run2 ;
		  for (i = j = 0 ; j < keySetMax (ks1) ; j++)
		    {
		      run2 = keySet (ks1, j) ;
		      rc2 =  arrayp (gx->runs, run2, RC) ;
		      if (rc2 && rc2->runs && ! rc2->addCounts)
			{
			  int k ;
			  for (k = 0 ; k < keySetMax (rc2->runs) ; k++)
			    keySet (ks2, i++) = keySet (rc2->runs, k) ;
			}
		      else
			keySet (ks2, i++) = run2 ;
		    }
		  keySetSort (ks2) ;
		  keySetCompress (ks2) ;
		  rc->runs = ks2 ;
		  keySetDestroy (ks1) ;
		}
	    }
	}
    }

  ac_free (h) ;
  fprintf (stderr, "// %s Found %d %s, %d runs or groups"
	   , timeShowNow () 
	   , dictMax (gx->geneDict)
	   , gx->isINTRON ? "introns" : (gx->isTranscript ? "transcripts" : "genes")
	   , dictMax (gx->runDict)
	   ) ; 
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, ", %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
  }

  return nn ;
} /* gxAceParse */

/*************************************************************************************/
/* parse the procaryote tsf file run->deep and run->Length */
static int gxDeepTsfParse (GX *gx, const char* fileName,BOOL metaData)
{
  AC_HANDLE h = ac_new_handle () ; 
  ACEIN ai ; 
  Array aa ;
  int nn = 0, gene, run, a1, a2 ;
  float tags, kb ;
  const char *ccp ;
  DC *dc = 0 ;  
  RC *rc = 0 ;
  GC *gc = 0 ;

  /*   Stack buf = stackHandleCreate (1000,h) ; */
  DICT *markerSelectDict = 0 ;
  DICT *markerRejectDict = 0 ;

  if (gx->markerSelectFileName)
    {
      ai = aceInCreate (gx->markerSelectFileName, gx->gzi, h) ;
      markerSelectDict = dictHandleCreate (10000, h) ;
      while (aceInCard (ai)) 
	{ 
	  ccp = aceInWord (ai) ;
	  if (ccp) dictAdd (markerSelectDict, ccp, 0) ;
	}

      ac_free (ai) ;
    }
  if (gx->markerRejectFileName)
    {
      ai = aceInCreate (gx->markerRejectFileName, gx->gzi, h) ;
      markerRejectDict = dictHandleCreate (10000, h) ;
      while (aceInCard (ai)) 
	{ 
	  ccp = aceInWord (ai) ;
	  if (ccp) dictAdd (markerRejectDict, ccp, 0) ;
	}

      ac_free (ai) ;
    }


  ai = aceInCreate (fileName, gx->gzi, h) ;
  if (! ai)
    messcrash ("cannot open input file -deep %s", fileName ? fileName : "stdin") ;
  aceInSpecial (ai, "\"\n\t") ;

  while (aceInCard (ai)) 
    { /* parse the tsf file */
      ccp = aceInWord (ai) ; /* gene */
      if (!ccp || *ccp == '#' || *ccp == '/') continue ;

      if (gx->maskDict &&  dictFind (gx->maskDict, ccp, 0))
	continue ;
      if (gx->keepIndex && gx->userGivenGenePlusMinusDict && ! dictFind (gx->userGivenGenePlusMinusDict, ccp, 0))
	continue ;
      if (markerSelectDict && ! dictFind(markerSelectDict, ccp, 0))
	continue ;
      if (markerRejectDict && dictFind(markerRejectDict, ccp, 0))
	continue ;

      dictAdd (gx->geneDict, ccp, &gene) ;
      gc = arrayp (gx->genes, gene, GC) ; /* make room */

      ccp = aceInWord (ai) ; /* run */
      if (! ccp) continue ;
      
      run = 0 ;
      if (! gx->runListFileName)
	dictAdd (gx->runDict, ccp, &run) ;
      else
	dictFind (gx->runDict, ccp, &run) ;
      if (run <= 0)
	continue ;

      rc = arrayp (gx->runs, run, RC) ;
      rc->cumulDone = TRUE ;
      aa = rc->aa ;
      if (! aa)
	aa = rc->aa = arrayHandleCreate (50000, DC, gx->h) ;


      ccp = aceInWord (ai) ; /* format  ignore */
      ccp = aceInWord (ai) ; /* chrom  ignore */
      aceInInt (ai, &a1) ;   /* gene coordinates ignore */
      aceInInt (ai, &a2) ;   /* gene coordinates ignore */

      aceInFloat (ai, &kb) ; /* kb */
      nn++ ;
      dc = arrayp (aa, gene, DC) ;
      tags = 10 * kb ; /* assume 100 bases */
      dc->seqs += tags ;
      dc->tags += tags ;
      dc->kb += kb ;

      rc->geneTags += tags ;
      rc->nReads += tags ;
      rc->nReadsOk += tags ;
      rc->kb += kb ;

      gc->tags += tags ;
      if (gc->capturedIntronGene)
	rc->captured_tags += tags ;
      gc->isGood = TRUE ;
      if (tags) rc->hasData = TRUE ;

    }
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, ", %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
  }

  return nn ;
} /* gxDeepTsfParse */

/*************************************************************************************/

static void gxSubsample (GX *gx)
{
  int run, runMax = gx->runs ? arrayMax (gx->runs) : 0 ;

  for (run = 1 ; run < runMax ; run++)
    {
      int iAA = 0 ;
      float zMin ;
      double z, targetKb = 0, oldTags ;

      RC *rc = arrayp (gx->runs, run, RC) ;

      if (rc->geneTags < 1000.0 * gx->subsample)
	{
	  if (0)
	    {
	      int n = rc->aa ? arrayMax (rc->aa) : 0 ;
	      rc->tags = rc->seqs = rc->targetKb = 1 ;
	      if (n)
		rc->aa = arrayCreate (n, DC) ;
	    }
	  if (0) fprintf (stderr, "#SUBSAMPLING %d %s %.0f %.0f  keep\n", gx->subsample, dictName (gx->runDict, run), rc->tags, rc->geneTags) ;
	  continue ;
	}
      zMin = gx->subsample * 1000.0 / rc->geneTags ;
      oldTags = rc->geneTags ;
      rc->geneTags = 0 ;

      for (iAA = 0 ; iAA < 2 ; iAA++)
	{
	  Array aa = iAA ? rc->aa : rc->antiAa ;
	  
	  if (aa && arrayMax (aa))
	    {
	      int gene, gMax = arrayMax (aa) ;
	      DC *dc ;
	      GC *gc ; 
	      
	      for (gene = 0, dc = arrp (aa, 0, DC) ; gene < gMax ; dc++, gene++)
		{
		  int j = dc->tags, k = 0 ;

		  while (j--)
		    if (randfloat() < zMin)
		      k++ ;
		  z = dc->tags ? (1.0 * k)/dc->tags : 0 ;
		  dc->seqs *= z ;
		  dc->kb *= z ;

		  gc =  arrp (gx->genes, gene, GC) ; /* make room */
		  gc->tags = gc->tags - dc->tags + k ;
		  dc->tags = k ;
		  rc->geneTags += dc->tags ;
		  targetKb += dc->kb ;
		}
	    }
	}
      if (0) fprintf (stderr, "#SUBSAMPLING %d %s %.0f %.0f %.0f\n", gx->subsample, dictName (gx->runDict, run), rc->tags, oldTags, rc->geneTags) ;
      z = rc->targetKb > 0 ? targetKb / rc->targetKb : 0 ;
      rc->targetKb = targetKb ;
      rc->anyTargetKb *= z ;
      rc->genomicKb *= z ;
      rc->intergenicKb *= z ;
      rc->intergenicDensity *= z ;
    }

  return ;
} /* gxSubsample */

/*************************************************************************************/
/*************************************************************************************/

#ifdef JUNK
/* useful for debugging */
static void gxHackZeroCounts (GX *gx)
{
  int gene, run, runG ;
  RC *rc ;
  DC *dc ;
  Array aa ;
  KEYSET runs ;
  int i, runMax = arrayMax (gx->runs) ;

  for (run = 1, runG = 0 ; !runG && run < runMax ; run++)
    if (! strcmp  (dictName (gx->runDict, run), "NB_male")) /* Stg1_MYCN_not_amplified_117 */
      runG = run ;
  
  if (! runG ||
      ! (rc = arrayp (gx->runs, runG, RC)) ||
      ! (runs = rc->runs) ||
      ! keySetMax (runs)
      )
    return ;
  
  if (1 && dictFind (gx->geneDict, "X__ACTG1", &gene))
    {
      for (i = 0 ; i < keySetMax (runs) ; i+=1)
	{
	  run = keySet (runs, i) ;
	  rc = arrayp (gx->runs, run, RC) ;
	  aa = rc->aa ;
	  if (aa) 
	    {
	      dc = arrayp (aa, gene, DC) ;
	      if (dc)
		dc->tags = dc->kb = dc->seqs = 0 ;
	    }
	}
      return ;
    }    

  for (i = 0 ; i < keySetMax (runs) ; i++)
    {
      int hackNullCount = ( (i % 2 == 0) ? 200 : 400) ;

      run = keySet (runs, i) ;
      rc = arrayp (gx->runs, run, RC) ;
      aa = rc->aa ;

      if (aa)
	for (gene = 1, dc = arrayp (aa, gene, DC) ; gene < arrayMax (aa) && gene < hackNullCount ; gene++, dc++) 
	  dc->tags = dc->kb = dc->seqs = 0 ;
    }
}
#endif

/*************************************************************************************/

static void gxGlobalCounts (GX *gx)
{
  int gene, run, geneMax ;
  RC *rc ;
  DC *dc ;
  GC *gc ;
  double kb, t, t2, z, z2, zbig, zE, tE ;
  int removeLimit = gx->removeLimit, pass ;
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = 0 ; 
  Array aa ;
  int n, runMax = arrayMax (gx->runs) ;

  /* To compute the cumul avoid Globin
   * eliminate from the cumul the genes contributing over 1/removeLimit  default 1/50 =2% (was 1 per cent, or .1%) 
   */
  
  gx->isRunMax = 0 ;
  geneMax = arrayMax (gx->genes) ;
  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
    {
      int nGeneWithGeneId = 0 ;
      rc->run = run ;
      if (rc->runs)
	continue ;
      gx->isRunMax++ ;
      aa = rc->aa ;
      if (! aa) 
	continue ;
      dc = arrayp (aa, dictMax(gx->geneDict), DC) ; /* regularize the size of the aa table in all runs */

      kb = rc->prunedTargetKb = rc->targetKb ; rc->bigGenesKb = 0 ; rc->capturedBigGenesKb = 0 ;
      if (rc->bigGenes)  arrayMax (rc->bigGenes) = 0 ;
      if (rc->capturedBigGenes)  arrayMax (rc->capturedBigGenes) = 0 ;
      /* 2 passes :: first count how many kb are in genes with a geneId, then flag the big genes */
      for (pass = 0, kb = 0 ; pass < 2 ; pass++)
	{
	  t = t2 = tE = z = z2 = zE = zbig = nGeneWithGeneId = 0 ;
	  for (n = 0, gene = 1, dc = arrayp (aa, gene, DC), gc = arrayp (gx->genes, gene, GC) ; gene < arrayMax (aa) ; gene++, gc++, dc++) 
	    {
	      if (gc->geneGroup)
		continue ;
	      if (0 && rc->capture && gc->captures && keySetFind (gc->captures, rc->capture, 0))
		continue ;
	      if (!strncmp (dictName(gx->geneDict,gene), "ERCC-", 5))
		{
		  zE += dc->kb ;
		  tE += dc->tags ;
		}
	      else if (pass == 0 || removeLimit * dc->kb < kb)  /* dromadaire 50 seems best, 1000 3000 is terrible */
		{                   /* the removeLimit is only applied in the second pass (pass == 1) */
		  z += dc->kb ;
		  t += dc->tags ;
		  if (gene < geneMax && 
		      ((gx->hasGoodProduct && gc->hasGoodProduct) || gc->geneId)
		      )
		    { 
		      nGeneWithGeneId++ ;
		      z2 += dc->kb ;
		      t2 += dc->tags ;
		    }
		  n++ ;
		}
	      else
		{ 
		  BG *bg ;
		  zbig += dc->kb ;
		  if (! ao) 
		    {
		      if (gx->isINTRON)
			ao = aceOutCreate (gx->outFileName, ".big_intron.txt", FALSE, h) ;
		      else
			ao = aceOutCreate (gx->outFileName, ".big_gene.txt", FALSE, h) ;
		      aceOutDate (ao, "###", gx->title) ;
		      if (gx->isINTRON)
			aceOutf (ao, "# Intron") ;
		      else if (gx->isSNP)
			aceOutf (ao, "# SNV") ;
		      else if (gx->isTranscript)
			aceOutf (ao, "# mRNA") ;
		      else
			aceOutf (ao, "# Gene") ;
		      aceOutf (ao, "\tRun\tSupporting reads") ;
		      if (gx->hasCapture)
			{
			  if (gx->isINTRON)
			    aceOutf (ao, "\tGene") ;
			  aceOutf (ao, "\tCaptured by\t Touched by") ;
			  aceOutf (ao, "\n") ;
			}
		    }
		  if (gx->isINTRON && ! gc->isIntron)
		    continue ;
		  if (! rc->bigGenes) 
		    rc->bigGenes = arrayHandleCreate (32, BG, gx->h) ;
		  if (gx->hasCapture && ! rc->capturedBigGenes) 
		    rc->capturedBigGenes = arrayHandleCreate (32, BG, gx->h) ;
		  bg = arrayp (rc->bigGenes, arrayMax (rc->bigGenes), BG) ;
		  bg->gene = gene ; bg->kb = dc->kb ;
		  aceOutf (ao, "%s\t%s\t%g", dictName(gx->geneDict,gene), dictName (gx->runDict, run), dc->tags) ;
		  if (gx->hasCapture)
		    {
		      KEYSET ks = gc->captures ;
		      int i ;
		      char *sep = "\t" ;
		      aceOutf (ao, "\t") ;
		      if (gx->isINTRON)
			{
			  KEY i2g = gc->intronGene ;
			  GC *gc2 = i2g ? arrayp (gx->genes, i2g, GC) : 0 ;
			  aceOutf (ao, "\t%s", i2g ? dictName(gx->geneDict,i2g) : "") ;
			  ks = gc2 ? gc2->captures : 0 ;
			}
		      
		      if (ks)
			for (i = 0 ; i < keySetMax (ks) ; i++)
			  {
			    int capture = keySet (ks, i) ;
			    if (capture && capture == rc->capture)
			      {
				BG *bgc = arrayp (rc->capturedBigGenes, arrayMax (rc->capturedBigGenes), BG) ;
				bgc->gene = gene ; bgc->kb = dc->kb ;
				rc->capturedBigGenesKb += dc->kb ;
			      }
			    if (capture)
			      {
				const char *ccp = dictName (gx->captureDict, capture) ;
				if (ccp && (1 || ccp[2] == 0)) 
				  {
				    aceOutf (ao, "%s%s", sep, ccp) ;
				    sep = "," ;
				  }
			      }
			  }
		    
		      if (gx->hasCaptureTouch)
			{
			  KEYSET ks = gc->captures ;
			  KEYSET ksT = gc->capturesTouch ;
			  int i ;
			  char *sep = "\t" ;
			  if (gx->isINTRON)
			    {
			      KEY i2g = gc->intronGene ;
			      GC *gc2 = i2g ? arrayp (gx->genes, i2g, GC) : 0 ;
			      ks = gc2 ? gc2->captures : 0 ;
			      ksT = gc2 ? gc2->capturesTouch : 0 ;
			    }
			  if (ksT)
			    for (i = 0 ; i < keySetMax (ksT) ; i++)
			      {
				int capture = keySet (ksT, i) ;
				if (capture && ! (ks && keySetFind (ks, capture, 0)))
				  {
				    const char *ccp = dictName (gx->captureDict, capture) ;
				    if (ccp && ccp[2] == 0) 
				      {
					aceOutf (ao, "%s%s", sep, ccp) ;
					sep = "," ;
				      }
				  }
			      }
			}
		    }
		  aceOutf (ao, "\n") ;
		}
	    }
	  if (3 * (z + zbig) < rc->targetKb) /* drop the true counts, if we are using a very reduced list */
	    {
	      kb = rc->prunedTargetKb = rc->targetKb ; 
	    }
	  else
	    {
	      if (nGeneWithGeneId > 10000)
		{
		  kb = z2 ;
		}
	      else
		{
		  if (2 * zE > z) {z += zE ; t += tE ; } /* beware of pure ERCC runs */
		  kb = z ; t2 = t ;
		}
	      if (pass)
		{
		  rc->prunedTargetKb = kb ; /* rather normalize via the genes WITH geneid, if we have them available */
		  rc->tags = t2 ;
		  rc->bigGenesKb = zbig ;
		}
	    }
	  if (0)
	    fprintf (stderr, "prunedTargetK=%lf intergenicKb=%ff\n", rc->prunedTargetKb, rc->intergenicKb* gx->seaLevel/3) ;
	  if (0)
	    rc->prunedTargetKb += rc->intergenicKb * gx->seaLevel / 3.0 ;
	}
    }  
  ac_free (h) ;
  return ;
} /* gxGlobalCounts */

/*************************************************************************************/
/* once eigen genes are declared and their tags are counted, they behave normally everywhere else */
static void gxGeneGroupCount (GX *gx)
{
  int level, geneMax = arrayMax (gx->genes) ;
  int gene, runMax = arrayMax (gx->runs) ;
  GC *gc ;

  for (level = 1 ; level <= gx->maxGeneGroupLevel ; level++)
    for (gene = 1, gc = arrp (gx->genes, gene, GC) ; gene < geneMax ; gc++, gene++)
      {
	if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
	  continue ;
	if (gc->Group_level == level)
	  {
	    Array aa, gggg = gc->geneGroup ;
	    int run, j, jMax = arrayMax (gggg) ;
	    GC *gc1 ;
	    RC *rc ;
	    BOOL hasLength = gc->length ? TRUE : FALSE ; /* may have been provided as .ace data */
	    if (gx->targeted && ! gc->targeted)
	      continue ;
	    if (0)
	      if (gx->captured && ! gc->captured)
		continue ;
	    if (! hasLength)
	      gc->length =  0 ;
	    gc->tags = 0 ;
	    for (j = 0 ; j < jMax ; j++)
	      {
		DC *dc, *dc1 ;
		int gene1 = arr (gggg, j, GGG).gene ;
		if (gene1 && gene1 < keySetMax (gx->genes))
		  {
		    gc1 = arrp (gx->genes, gene1, GC) ;
		    if (0)
		      if (gx->captured && ! gc1->captured)
			continue ;
		    if (gx->targeted && ! gc1->targeted)
		      continue ;
		    if (! hasLength) gc->length += gc1->length ;  
		    gc->tags += gc1->tags ;
		    gc->targeted |= gc1->targeted ;
		    gc->captured |= gc1->captured ;
		    for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++)  
		      {
			if (rc->runs)
			  continue ;
			aa = rc->aa ;
			if (! aa || gene > arrayMax (aa) || gene1 > arrayMax (aa)) 
			  continue ;
			dc = arrayp (aa, gene, DC) ;
			dc1 = arrayp (aa, gene1, DC) ;
			dc->seqs += dc1->seqs ;
			dc->tags += dc1->tags ;
			dc->kb += dc1->kb ;
		      }
		  }
	      }
	  }
      }
} /* gxGeneGroupCount */

/*************************************************************************************/
/*************************************************************************************/
/* In each group of runs, cumulate the counts */ 
static void gxOneGroupCumul (GX *gx, RC *rc, int myGene)
{
  AC_HANDLE h = ac_new_handle () ;
  int run, run2, i, gene, nRuns ;
  int PolyA_selected_RNA = 0 ;
  double accessibleLength5kb = 0 ;
  double accessibleLength8kb = 0 ;
  RC *rc2 ;
  DC *dc, *dc2 ;
  DC *antiDc, *antiDc2 ;
  DDC *ddc, *ddc2 ;
  Array aa, daa, aa2, daa2 ;
  Array antiAa, antiAa2 ;
  BitSet bbLow = myGene ? 0 : bitSetCreate (dictMax (gx->geneDict) + 1, h) ;
  double Kdelta = 0, Kshift = 1010 ;  /* Knuth algo, see keyset.c, 1010 is a good guess at the average value */

  if (rc->runs && !rc->cumulDone)
    {
      run = rc->run ;
      rc->cumulDone = TRUE ;
      keySetSort (rc->runs) ;
      keySetCompress (rc->runs) ;
      if (bbLow) bitSetReCreate (bbLow, dictMax (gx->geneDict) + 1) ;
      nRuns = 0 ;
      rc->nRuns = keySetMax (rc->runs) ;
      if (gx->isMA) rc->isMA = gx->isMA ;
      aa = rc->aa ; daa = rc->daa ; antiAa = rc->antiAa ;
      rc->seqs = rc->tags = rc->geneTags = rc->prunedTargetKb = 0 ; 
      rc->genomicKb = rc->intergenicKb = rc->targetKb = rc->anyTargetKb = 0 ;
      rc->Genes_with_one_reliable_index = 0 ;
      rc->PolyA_selected_RNA = 0 ; 
      PolyA_selected_RNA = 0 ;  
      accessibleLength5kb = rc->accessibleLength5kb = 0 ; 
      accessibleLength8kb = rc->accessibleLength8kb = 0 ; 
      if (! rc->addCounts)
	{
	  rc->zeroIndex = 0 ;
	  rc->Low_index = 0 ;
	}

      if (! aa)
	{
	  if (0)
	    fprintf (stderr, ".... Creating group cumul level %d run %d %s\n"
		     , rc->groupLevel
		     , run
		     , dictName (gx->runDict, run)
		     ) ;
	  aa = rc->aa = arrayHandleCreate (50000, DC, gx->h) ;
	  if (! gx->isTranscript && ! gx->isSNP && ! gx->isINTRON)
	    daa = rc->daa = arrayHandleCreate (50000, DDC, gx->h) ;
	}
      else
	{
	  aa = rc->aa = arrayReCreate (aa, arrayMax (aa), DC) ;
	  if (! gx->isTranscript && ! gx->isSNP  && ! gx->isINTRON)
	    daa = rc->daa = arrayReCreate (daa, arrayMax (daa), DDC) ;
	}
      if (gx->anti)
	{ 
	  if (! antiAa)
	    antiAa = rc->antiAa = arrayHandleCreate (50000, DC, gx->h) ;
	  else
	    antiAa = rc->antiAa = arrayReCreate (antiAa, arrayMax (antiAa), DC) ; 
	}
      for (i = 0 ; i < keySetMax (rc->runs) ; i++)
	{
	  int gMax ;
	  run2 = keySet (rc->runs, i) ;
	  rc2 = arrayp (gx->runs, run2, RC) ;
	  aa2 = rc2->aa ;
	  daa2 = rc2->daa ;
	  antiAa2 = rc2->antiAa ;
	  if (! aa2)
	    continue ;
	  
	  
	  nRuns++ ;
	  rc->nReads += rc2->nReads ;
	  rc->seqs += rc2->seqs ;
	  rc->tags += rc2->tags ;
	  rc->geneTags += rc2->geneTags ;

	  rc->nReadsOk += rc2->nReadsOk ;
	  rc->kb += rc2->kb ;
	  rc->nerr += rc2->nerr ;
	  rc->a2g += rc2->a2g ;
	  rc->partial += rc2->partial ;
	  rc->orphan += rc2->orphan ;
	  rc->badTopo += rc2->badTopo ;
	  rc->multi += rc2->multi ;
	  rc->multi2 += rc2->multi2 ;


	  rc->targetKb += rc2->targetKb ;
	  rc->anyTargetKb += rc2->anyTargetKb ;
	  rc->genomicKb += rc2->genomicKb ;
	  rc->intergenicKb += rc2->intergenicKb ;
	  rc->prunedTargetKb += rc2->prunedTargetKb ;
	  rc->bigGenesKb += rc2->bigGenesKb ;
	  rc->capturedBigGenesKb += rc2->capturedBigGenesKb ;
	  PolyA_selected_RNA += rc2->PolyA_selected_RNA ;
	  accessibleLength5kb += rc2->accessibleLength5kb * rc2->tags ;
	  accessibleLength8kb += rc2->accessibleLength8kb * rc2->tags ;
	  
	  if (! rc->addCounts)
	    {
	      rc->zeroIndex += rc2->zeroIndex ;
	      rc->Low_index += rc2->Low_index ;
	    }
	  gMax = myGene && myGene < arrayMax (aa2) ? myGene + 1 : arrayMax (aa2) ;
	  for (gene = myGene ? myGene : 1, dc2 = arrayp (aa2, gene, DC) ; gene < gMax ; gene++, dc2++) 
	    {
	      GC *gc = gene < arrayMax (gx->genes) ? arrp (gx->genes, gene, GC) : 0 ;
	      if (! gc || ! dc2->index) 
		continue ;
	      if (gx->targeted && (! gc || ! gc->targeted))
		continue ;
	      if (0)
		if (gx->captured && (! gc || ! gc->captured))
		  continue ;
	      if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
		continue ;

	      dc = arrayp (aa, gene, DC) ; 
	      antiDc = antiAa ? arrayp (antiAa, gene, DC) : 0 ; 
	      antiDc2 = antiAa2 ? arrayp (antiAa2, gene, DC) : 0 ; 
	      ddc = daa ? arrayp (daa, gene, DDC) : 0 ;
	      ddc2 = daa2 ? arrayp (daa2, gene, DDC) : 0 ;
	      if (gx->isSNP)
		{
		  if (dc2->tags)
		    {
		      dc->seqs += dc2->seqs ; /* count */
		      dc->tags += dc2->tags ; /* cover */
		      dc->index += dc2->index ;
		      if (dc->kb  == 0 || dc->indexMin > dc2->index)
			dc->indexMin = dc2->index ;
		      if (dc->indexMax < dc2->index)
			dc->indexMin = dc2->index ;
		      dc->kb   ++ ; /* count the runs in which the SNP is measured */
		    }
		}
	      else
		{
		  if (nRuns == 1)
		    {
		      dc->indexMin = dc2->index ;
		      dc->variance = 0 ;
		      dc->index = 0 ;
		      if (antiDc && antiDc2)
			{
			  antiDc->indexMin = antiDc2->index ;
			  antiDc->variance = 0 ;
			  antiDc->index = 0 ;
			}
		    }
		  if (dc->indexMin > dc2->index)
		    dc->indexMin = dc2->index ;
		  if (dc->indexMax < dc2->index)
		    dc->indexMin = dc2->index ;
		  
		  dc->seqs += dc2->seqs ;
		  dc->tags += dc2->tags ;
		  dc->kb   += dc2->kb ;  
		  
		  if (antiDc && antiDc2)
		    {
		      if (antiDc->indexMin > antiDc2->index)
			antiDc->indexMin = antiDc2->index ;
		      if (antiDc->indexMax < antiDc2->index)
			antiDc->indexMin = antiDc2->index ;
		      antiDc->seqs += antiDc2->seqs ;
		      antiDc->tags += antiDc2->tags ;
		      antiDc->kb   += antiDc2->kb ;  
		    }
		  if (0 && gene == 9009 && run == 194)
		    fprintf (stderr, "Gene %d:%s run %d:%s run2 %d:%s tags %g %g\n"
			     , gene, dictName(gx->geneDict,gene)
			     , run,  dictName(gx->runDict, run)
			     , run2,  dictName(gx->runDict, run2)
			     , dc2->tags, dc->tags
			     )
		      ;  
		  if (ddc && ddc2)
		    {
		      ddc->nReads += ddc2->nReads ;
		      ddc->nReadsOk += ddc2->nReadsOk ;
		      ddc->nerr     += ddc2->nerr ;  
		      ddc->a2g     += ddc2->a2g ;  
		      ddc->partial  += ddc2->partial ;  
		      ddc->orphan   += ddc2->orphan ;  
		      ddc->badTopo   += ddc2->badTopo ;  
		      ddc->multi   += ddc2->multi ;  
		      ddc->multi2   += ddc2->multi2 ;  
		    }
		  
		  if (! dc2->index)
		    {
		      fprintf (stderr, "gxOneGroupCumul, uninitialized index group=%s (level %d)  run=%s (level %d)gene=%s (g=%d)\n"
			       , dictName (gx->runDict,run), rc ? rc->groupLevel : -99
			       , dictName (gx->runDict,run2), rc2 ? rc2->groupLevel : -99
			       , dictName (gx->geneDict,gene), gene
			       ) ; 
		      fprintf (stderr, "May be it is a synchro problem  try :\n    rm tmp/GENERUNS/*.ace\n") ;
		      fprintf (stderr, "then  rerun MAGIC g4\n") ;
		      messcrash ("Group cumul  uninitialized index group") ;
		      gene = arrayMax (aa2) ;
		      invokeDebugger() ;
		      continue ;
		    }
		  Kdelta = dc2->index - Kshift - dc->index ; 
		  dc->index += Kdelta/nRuns ;  /* 2013_12_09 always use the index of the gene, Knuth rolling average */
		  dc->variance += Kdelta * Kdelta  * (nRuns - 1)/nRuns ;
		  dc->isLow += dc2->isLow ;
		  if (bbLow && ! dc2->isLow) bitSet (bbLow, gene) ;
		  if (0 && 
		      (
		       (!gx->isTranscript && gene == 55690) ||
		       (gx->isTranscript && gene == 142153)
		       ) && (run > 1)
		      )
		    {
		      printf ("GrCumul\t%s\t%s\t%s\t%d\t%.4f\t%.4f\t%.2f\t%.2f\t%f\n"
			      , dictName (gx->geneDict,gene), dictName (gx->runDict,run)
			      , dictName (gx->runDict,run2) , nRuns 
			      , dc2->index - 1000, dc->index
			      , dc2->tags, dc->tags/nRuns  
			      , dc->variance/nRuns
			      ) ;
		    }
		}
	    }
	}
      if (nRuns)
	{
	  int gMax = myGene && myGene < arrayMax (aa) ? myGene + 1 : arrayMax (aa) ;
	  for (gene = myGene ? myGene : 1, dc = arrayp (aa, gene, DC) ; gene < gMax; gene++, dc++) 
	    {
	      GC *gc = gene < arrayMax (gx->genes) ? arrp (gx->genes, gene, GC) : 0 ;
	      if (gx->targeted && (! gc || ! gc->targeted))
		continue ;
	      if (0)
		if (gx->captured && (! gc || ! gc->captured))
		  continue ;
	      ddc = daa ? arrayp (daa, gene, DDC) : 0 ;
	      if (gx->isSNP)
		{
		  if (rc->addCounts)
		    {
		      int count = dc->seqs ;
		      int cover = dc->tags ;
		      
		      if (cover >= 10 && (myGene + dc->kb) > nRuns/2)
			{
			  float index = 0 + 20.0 * ((float)count)/((float)cover +.0001) ; 
			  dc->index = index + 1000 ;
			}
		      else
			dc->tags = dc->seqs = dc->kb = dc->index = 0 ;
		    }
		  else
		    {
		      if (dc->kb > nRuns/2) 
			{
			  dc->index /= dc->kb ;
			}
		      else
			dc->tags = dc->seqs = dc->kb = dc->index = 0 ;
		    }
		}
	      else
		{
		  if (1 || !rc->addCounts || gx->keepIndex)
		    {
		      
		      dc->index += Kshift ; /* is already normalized by Knuth  in add count mode the index will be recomputed */
		      if (nRuns <= 1)
			dc->variance = -999 ; 
		      else
			{
			  /* up to 2019_12_03  we have 2 errors
			   * we were dividing by N, which underestimates the variance for small N
			   * and we seem to have an algorithm error when implementing Knuth
			   */
			  Kdelta =  dc->variance/nRuns ;         /* 2016?      population variance */
			  Kdelta =  dc->variance/(nRuns - 1) ;   /* 2019_12_03 sample estimated variance */
			  dc->variance = Kdelta ;   
			      /* standard-deviation: sigma =  sqrt (Kdelta) */
			      /* coef variation on log normal data: CV = sqrt (2^variance - 1)  (CF wikipedia) */
			}
		      if (0 &&  gene == 55690)
			printf ("Final %s idx=%f var=%f sig=%f\n"
				, dictName(gx->geneDict, gene)
				, dc->index-1000, Kdelta, dc->variance) ;
		      
		      if (0 &&  /* 2016_03_28 : do NOT divide the counts ; was take the count averages */
			  ! rc->addCounts
			  )  
			{ 
			  if (nRuns)
			    {
			      dc->seqs  /= nRuns ;
			      dc->tags  /= nRuns ;
			      dc->kb    /= nRuns ;  
			      
			      if (ddc)
				{
				  ddc->nReads    /= nRuns ;
				  ddc->nReadsOk  /= nRuns ;
				  ddc->nerr      /= nRuns ;   
				  ddc->a2g       /= nRuns ;   
				  ddc->partial   /= nRuns ;
				  ddc->orphan    /= nRuns ; 
				  ddc->badTopo   /= nRuns ;  
				  ddc->multi     /= nRuns ;  
				  ddc->multi2    /= nRuns ;
				}
			    }
			}
		      if (! rc->addCounts)
			dc->isLow = (5 * dc->isLow >= 4 * nRuns) ? 1 : 0 ; /* was 3 : 2 modified 2014_06_04 */
		      if (rc->addCounts)
			dc->index  = 0 ; /* will be computed later */ 
		    }
		}
	    }
	  
	  rc->PolyA_selected_RNA = PolyA_selected_RNA > nRuns/2 ? TRUE : FALSE ;
	  rc->accessibleLength5kb = rc->tags ? accessibleLength5kb / rc->tags : 0 ; 
	  rc->accessibleLength8kb = rc->tags ? accessibleLength8kb / rc->tags : 0 ; 

	  if (! rc->addCounts)
	    {
	      rc->Genes_with_one_reliable_index = bitSetCount (bbLow) ;

	      rc->zeroIndex /= nRuns ;
	      rc->Low_index /= nRuns ;

	      if (0) /* 2016_03_28 : do NOT divide the counts ; was take the count averages */
		{
		  if (nRuns)
		    {
		      rc->targetKb         /= nRuns ;
		      rc->anyTargetKb      /= nRuns ; 
		      rc->prunedTargetKb   /= nRuns ; 
		      rc->bigGenesKb       /= nRuns ; 
		      rc->genomicKb        /= nRuns ;
		      rc->intergenicKb     /= nRuns ;
		    }
		}
	    }
	}
      else
	{ arrayDestroy (rc->aa) ; arrayDestroy (rc->daa) ; }
    }
  
  fflush (stdout) ;
  ac_free (h) ;
  return ;
} /* gxOneGroupCumul */

/*************************************************************************************/
/*************************************************************************************/

static float gxComputeOneIndex (int run, int gene, GX *gx, DC *dc, int *isLowp)
{
  float index = 0 ;  /* depends on total counts */
  double damper = 1000 ; /* damper in bp */
  float average_tag_ln ;
  double noise_bp = 0, noise_tags = 0, bonusTags ;
  double logDeux ;
  double z, abp, bp, genomicKb ;
  int ln, ln0 ;
  int accessibleLength ;
  /* static int debugDamper = 0 ; */
  logDeux = log((double)2.0) ;
  RC *rc ;
  GC *gc = gene < arrayMax (gx->genes) ? arrp (gx->genes, gene, GC) : 0 ;
  BOOL isLow = FALSE ;
  double wall = gx->wall > 1 ? gx->wall : 1 ;
  BOOL isIntron = gx->isINTRON || (gc && gc->isIntron) ;

  if (0)
    {
      if (0 || ! strncmp (dictName(gx->geneDict, gene),"AAK1.", 5) ||
	  (
	   0 || ! strcmp (dictName(gx->runDict, run),"HR_TERTlow_MYCNosLow_46") || ! strcmp (dictName(gx->runDict, run),"testTit")
	   )
	  )
	invokeDebugger () ;
    }
  *isLowp = FALSE ;
  if (! strcmp (dictName (gx->runDict,run), "Rhs1603"))
    rc = 0 ; /* hack to create a break point in gdb */
  if (! strcmp (dictName (gx->runDict,run), "ChronicFatigueGroup2"))
    rc = 0 ; /* hack to create a break point in gdb */

  rc = arrayp (gx->runs, run, RC) ;

  if (rc->isPairedEnd)   /* mieg 2014_09_07 */
    wall *= 2 ;  

  if (rc->isSNP)
    {       
      return dc->index - 1000 ;
    }
  if (gx->keepIndex && ! rc->runs)
    {   
      if (! rc->zeroIndex ||  dc->index - 1000 < rc->zeroIndex)
	rc->zeroIndex = dc->index - 1000 ;
      return dc->index - 1000 ;
    }

  if (! isIntron && ! gx->target_class && rc->tags < 90000)  /* 2022_06 allow 100k, useful for subsampling studies */
    return -999 ;
  if (rc->tags < 1000 && ! gx->target_class) 
    return -999 ;

  average_tag_ln  = 1000 * rc->prunedTargetKb / (1 + rc->tags) ;

  ln = ln0 = gc ? gc->length : 0 ;	    

  if (rc->PolyA_selected_RNA)
    { /* case -pA polyA selected */
      if (rc->accessibleLength8kb > 0)
	accessibleLength = rc->accessibleLength8kb ;
      else if (rc->accessibleLength5kb > 0)
	accessibleLength = rc->accessibleLength5kb ;
      else
	accessibleLength = 3000 ;

      if (ln < 20) ln = ln0 = 1000 ; /* synchronize with  baSaturatedGeneLengths */
    }
  else
    { 
      if (rc->accessibleLength8kb > 5000)
	accessibleLength = rc->accessibleLength8kb ;
      else if (rc->accessibleLength5kb > 5000)
	accessibleLength = rc->accessibleLength5kb ;
      else
	accessibleLength = 5000 ;

      if (ln < 20) ln = ln0 = 1000 ; /* synchronize with  baSaturatedGeneLengths */
    }
  if (gc && gc->geneGroup)
    accessibleLength *= arrayMax (gc->geneGroup) ;

  if (ln > accessibleLength) ln = accessibleLength ;
  if (gx->isMRNAH && gc && ! gc->geneGroup  && gc->notIsTranscriptA && ln > 400) ln = 400 ;

  if (gx->DGE || rc->DGE) 
    ln = 4000 ;
  if (0) ln = 1000 ;
  if (gx->isMicroRNA) { ln = 1000 ;  average_tag_ln = 35 ; }
  if (average_tag_ln < 35)    /* we cannot usefully map shorter reads */
    average_tag_ln = 35 ;

  if (gx->isMA)
    {
      ln = ln0 = average_tag_ln = gx->isMA ;
    }

  /* 2012_06_03, add a new correction, 
   * estimate the new exonic + intronic as z% of the annotated exon, this number of reads is removed from
   * the genomic read to estimate the truly genomic random reads
   * we measured z=5.8% (3.6% new exon or junction + 2.2 intronic, in Brain 
   * and 4.3% in UHR (2.6% in new exon and junctions and 1.7% in introns) 
   * this removes most of the genomic reads in good experiments, but not much in bad ones
   * we could then use a higher value of the sea level which now just affects the bad experiments
   */

  if (LEMING || isIntron)
    { noise_bp =  noise_tags = genomicKb = 0 ; wall = 4 ; }
  else
    {
      z = 5.0/100.0 ;          /* dromadaire */
      genomicKb = rc->genomicKb - z * rc->targetKb ;
      if (genomicKb < 0) 
	genomicKb = 0 ;
      genomicKb = gx->seaLevel * genomicKb / 3 ;  
      if (rc->DGE || gx->DGE || rc->pA == 1 || (!rc->pA && gx->pA == 1)) genomicKb /= 10 ;      /* total, 90% go to introns and should not be attributed to genomic contamination */

      if (rc->intergenicKb)
	genomicKb = rc->intergenicKb * gx->seaLevel/3 ; /* was 0.4 seaLevel * rc->intergenicKb, sellevel default == 3 */
      if (rc->intergenicKb > rc->prunedTargetKb)  /* 2020_01_27 we are facing a DNA-seq experiment */
	genomicKb *= 2 ;
      noise_bp = ln0 * genomicKb/gx->genomeLengthInKb ;  /* here use the true length of the sequence, do not staturate at 3kb  */
     
      noise_tags = noise_bp/average_tag_ln ;
      if (noise_tags < 10)
	noise_bp = 3 * noise_tags ;
      else if (noise_tags > 500)
	noise_bp = 1.2828 * noise_tags ;
      else
	noise_bp = noise_tags  + 6.3245 * sqrt (noise_tags) ;
      noise_bp *= average_tag_ln ;
      noise_tags = noise_bp/average_tag_ln ;
    }

  /* mieg: 2011_03_14 we try to substract the genomic noise
   * note that we should measure the effective length rather than rely on annotations
   * CF Double stranded or star-seq normalisation is supposed to yield a more uniform coverage
   * and may also influence the index is some way or other
   */
  index = -999 ;

  bp = 1000 * dc->kb - noise_bp ;  
  abp = 1000 * (rc->prunedTargetKb + rc->intergenicKb * gx->seaLevel/3.0) ; /* 2020_01_27 add the intergenic tags to the denominator */
  bonusTags = dc->tags - noise_tags - wall ;

  /****************** new method 2013_04_13
   * we separate damping the (quantized)  number of tags 
   * from the denominator, normalizing the size of the genes and the run
   * for large counts, the formula is not modified
   * for zero counts, the low value is around 7 for 60M reads (Rhs844)
   */

  if (! (LEMING && isIntron)) /* new method 2013_04_13*/
    {
      double tt ;
      int lf ;

      if (bonusTags < 0)
	{ 
	  isLow = TRUE ;
	  *isLowp = TRUE ;
	}
      else
	*isLowp = FALSE ;

      if (gx->isMA)
	lf = 0 ;
      else if (rc->fragmentLength)
	lf = rc->fragmentLength ;
      else
	lf = 300 ; /* assumed length of the fragment, we should have measured it */ 

      z = 1.0E12/ ((ln - lf > lf ? ln - lf : lf) * abp) ;
      
      if (0)  /* dromadaire 2015_06_27, i revert to if (0) */
	{
	  /* dromadaire, method 2013_05_20 to 2013_07_31, gave bad predictive power for fatigue using damper=3+noise */
	  bp = 1000 * dc->kb ;  
	  tt = bp/average_tag_ln ; z *= average_tag_ln ;
	  if (tt < 0) tt = 0 ;
	  damper = 25 + noise_tags ; /* 2013_05_20   before we substracted noise_bp from bp and had damper=5 */
	  /* dromadaire, it seems we used damper = 5 + noise_tags 2015_01_31 with good results */
	  if (1 && gx->isMRNAH && gc && gc->notIsTranscriptA) damper =  4 + noise_tags ;
	}
      else  /* excellent change: OS and Fav d not change, EF HREF HROS go up one point */
	{ /* 2013_07_31 revert to previous method */
	  bp = 1000 * dc->kb - noise_bp ;  
	  tt = bp/average_tag_ln ; z *= average_tag_ln ;
	  if (tt < 0) tt = 0 ; 
	  damper = 5 ;
	  if (1) damper = 3 ; /* SEQC2 test 2021_09_03 */
	  damper = 4 ; /* 2023_09_21, so if tt=0 log((tt+sqrt(damper +tt^2)/2)=log(1)=0 */
	}

      if(gx->isMA) damper = 9 ; /* because we count a local coverage, not a number of reads */
      index = (log(( tt + sqrt (damper  + tt * tt))/((double)2.0)) + log(z))/logDeux ;	 
      if (0 && gene == 224 && dc->tags == 0)
	fprintf (stderr, "%s\t%s\tnoise_tags=%.2f\tdamper=%.2f\tindex=%.2f\n"
		 , dictName(gx->geneDict, gene)
		 , dictName(gx->runDict, run)
		 ,  noise_tags, damper, index
		 ) ;

      if (gx->malus) 
	index -= gx->malus ;
      if (! rc->zeroIndex)
	{
	  /* mieg 2023_09_21: rc->zero cannot depend on the particular gene 
	   * i choose 2kb
	   */
	  z = 1.0E12/ (2000 * abp) ; z *= average_tag_ln ; 
	  rc->zeroIndex = (log(sqrt(damper)/2.0) + log(z))/logDeux ; 
	}
      if (0 && !strcmp (dictName(gx->geneDict,gene),"Lgals3"))
	{
	  printf ("run %s Ln=%d  index=%.2f dc->kb=%g rc->prunedTargetKb=%g rc->targetKb=%g\n", dictName(gx->runDict, run),ln, index, dc->kb, rc->prunedTargetKb, rc->targetKb) ;
	}
      index += rc->shift_expression_index ;
      return index ;
    }
  else if (0 && gc->intronGene && isIntron)    /* relative read counts */
    { 
      DC *dc2 = arrp (rc->aa, gc->intronGene, DC) ;
      isLow = dc->tags < 10 ? TRUE : FALSE ; 
      index = 10 + log((4.0 + dc->tags)/(4+dc2->tags))/log(2.0) ;
    }
  else if (isIntron)    /* Naive Leming Shi NB Index */
    { 
      isLow = dc->tags < 10 ? TRUE : FALSE ; 
      index = 10 + log(1.0 + dc->tags)/log(2.0) ; 
    }
  else if (LEMING)   /* Naive Leming Shi NB Index */
    { 
      isLow = dc->tags < 10 ? TRUE : FALSE ; 
      index = log((1.0 + dc->tags)*1000000.0/(1+rc->tags))/log(2.0) ;
    }

  if (0 && !strcmp (dictName(gx->geneDict,gene),"Lgals3"))
    {
      printf ("run %s Ln=%d  index=%.2f dc->kb=%g rc->prunedTargetKb=%g rc->targetKb=%g\n", dictName(gx->runDict, run),ln, index, dc->kb, rc->prunedTargetKb, rc->targetKb) ;
    }
  *isLowp = isLow ;
  index += rc->shift_expression_index ;
  return index ;
} /* gxComputeOneIndex */

/*************************************************************************************/
/* given the counts, compute the index for each run/gene */
static int gxComputeAllIndex (GX *gx)
{
  int run, gene, nn = 0, iln = 0, pass ;
  DC *dc ;
  RC *rc ;
  GC *gc ;
  Array a, aa, lns = 0, indexes = 0 ;
  DICT *runDict = gx->runDict ;
  int runMax = arrayMax (gx->runs) ;
  int geneMax = arrayMax (gx->genes) ;
  /* int removeLimit = gx->removeLimit ; */
  double kb ;

  a = arrayCreate (dictMax (gx->geneDict), float) ;

  rc = arrayp (gx->runs, gx->runAny, RC) ;
  if (! rc->runs)
    {
      int i ;
      RC *rc2 ;

      rc->runs = keySetHandleCreate (gx->h) ;
      for (i = 0, run = 1, rc2 = arrayp (gx->runs, run, RC) ; run < runMax ; rc2++, run++) 
	if (! rc2->runs)
	  keySet (rc->runs, i++) = run ;
    }

  for (pass = -1 ; pass <= gx->maxGroupLevel ; pass++)
    {
      for (run = 1, rc = arrayp (gx->runs, run, RC) ; run < runMax ; rc++, run++) 
	{
	  if (rc->indexDone) continue ; 
	  if (pass != rc->groupLevel) continue ; 

	  rc->indexDone = 1 ;
	  if (rc->runs)
	    gxOneGroupCumul (gx, rc, 0) ; /* compute the tag counts of all groups + index of non-additive */
	  if (pass > 0 && ! rc->addCounts) 
	    continue ;
	  aa = rc->aa ;
	  if (! aa) 
	    continue ;
	  lns = arrayReCreate (lns, geneMax, int) ; iln = 0 ;
	  indexes = arrayReCreate (indexes, geneMax, int) ;
	  for (gene = 1, dc = arrayp (aa, 1, DC), kb = 0 ; gene < arrayMax (aa) ; gene++, dc++) 
	    {
	      gc = gene < geneMax ? arrayp (gx->genes, gene, GC) : 0 ;
	      if (!gc)
		continue ;
	      if (gx->targeted && ! gc->targeted)
		continue ;
	      if (0)
		if (gx->captured && (! gc->captured))
		  continue ;
	      if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
		continue ;
	      /* compute and register the index */

	      nn++ ;
	      dc->indexMin = dc->indexMax = dc->index = 1000 + gxComputeOneIndex (run, gene, gx, dc, &(dc->isLow)) ;
	      if (gx->medianCenter)
		array (indexes, gene, float) = dc->index ;
	      if (strncmp (dictName(gx->geneDict,gene), "ERCC-", 5))
		kb += dc->kb ;
	      if (gc) /* rc->title && strstr (stackText (gx->info, rc->title), "Liver")) */
		{
		  /* gc->minIndex is used as a lower bound when searching for DEG and in gxCorrelation */
		  if (gc->minIndex == 0)
		    gc->minIndex = gc->maxIndex = dc->index ;
		  if (gc->minIndex > dc->index )
		    gc->minIndex = dc->index ;
		  if (gc->maxIndex < dc->index )
		    gc->maxIndex = dc->index ;
		}
	      if (! dc->isLow)
		{
		  array (lns, iln++, int) = gc->length ;
		}
	    }
      	  if (gx->medianCenter)
	    {
	      float dz = 0 ;
	      arraySort (indexes, floatOrder) ;
	      dz = array(indexes, arrayMax(indexes)/2, float) - 1010 ;
	      for (gene = 1, dc = arrayp (aa, 1, DC) ; gene < arrayMax (aa) ; gene++, dc++) 
		{
		  GC *gc = arrayp (gx->genes, gene, GC) ; 
		  if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
		    continue ;
		  dc->indexMin = dc->indexMax = dc->index = dc->index - dz ;
		}
	      fprintf(stderr, "// MEDCENTER\t%s\t%.2f\n", dictName(gx->runDict,run), dz) ;
	    }
	  /* find the median length and compute the abaque */
	  if (iln > 0)
	    {
	      int oldLn, iDummy, nDummy = 0 ;
	      DC dcDummy ;
	      float indexDummy ;
	      float average_tag_ln ;
	      
	      arraySort (lns, intOrder) ;
	      gene = 1 ;
	      gc = gene < geneMax ? arrayp (gx->genes, gene, GC) : 0 ;
	      oldLn = gc->length ;
	      dc = arrayp (aa, 1, DC) ;
	      dcDummy = *dc ;
	      
	      average_tag_ln  = 1000 * rc->prunedTargetKb / (1 + rc->tags) ;
	      gc->length = 32 * average_tag_ln ;
	      /* alternative : gc->length = arr (lns , iln/2, int) ; */

	      for (iDummy = nDummy = 0 ; iDummy <= 30 ; iDummy++)
		{
		  if (! iDummy) nDummy = 0 ;
		  else if (iDummy == 1) nDummy = 1 ;
		  else nDummy <<= 1 ;
		  dcDummy.tags = nDummy ;
		  dcDummy.kb = average_tag_ln * nDummy/1000.0 ;
		  indexDummy = gxComputeOneIndex (run, gene, gx, &dcDummy, &(dcDummy.isLow)) ;
		  rc->abaque[iDummy] = indexDummy ;
		}
	      gc->length = oldLn ;
	    }      
	}
    }
  arrayDestroy (lns) ;

  for (run = 1, rc = arrayp (gx->runs, run, RC) ; run < runMax ; rc++, run++) 
    if (rc->nRuns && rc->indexDone < 2)
      {
	int ia, ira, runa, NX, NY, nRuns ;
	RC *rca = 0 ;
	DC* dca ;
	Array a = 0, aaa ;
	double z, X, X2, X3, X4, Y, Y2, L2, L3, L4 ;
	double C2[2], C3[4], C4[6], c2, Kshift = 1010 ;
	int isLow, isHigh ;

	/* compute the binomial coeff, needed for the L-coef, they depend on n: number of runs in the sample */
        C2[0] = C2[1] = C3[0] = C3[1] = C3[2] = C3[3] = C4[0] = C4[1] = C4[2] = C4[3] = C4[4] = C4[5] = 1 ;
	rc->indexDone = 2 ;
	aaa = rc->aa ;
	if (! aaa) 
	  continue ;
	dc = arrayp (aaa, dictMax (gx->geneDict) - 1, DC) ; /* make room */
	for (gene = 1, dc = arrayp (aaa, 1, DC) ; gene < arrayMax (aaa) ; gene++, dc++) 
	  if (1)
	    {
	      GC *gc = arrayp (gx->genes, gene, GC) ; 
	      if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
		continue ;

	      isLow = isHigh = 0 ;
	      NX = 0 ; X = X2 = X3 = X4 = 0 ;
	      NY = 0 ; Y = Y2 = 0 ;
	      L2 = L3 = L4 = 0 ;
	      a = arrayReCreate (a, dictMax (gx->runDict), float) ; ia = 0 ;
	      nRuns = rc->nRuns ;
	      for (ira = 0 ; ira < nRuns ; ira++)
		{
		  runa = keySet (rc->runs, ira) ;
		  rca = arrp (gx->runs, runa, RC) ;
		  aa = rca->aa ;
		  if (! aa) continue ;
		  array (a, ira, float) = 1000 + rca->zeroIndex ; /* make room, ensure we have a measure for every run */
		  dca = arrayp (aa, gene, DC) ;
		  if (dca && dca->index > 1000)
		    z = dca->index + rca->fineTune ;
		  else
		    {
		      z = rca->zeroIndex + rca->fineTune + 1000 ;  /* not yet available at this stage of the code, so == 0 */
		    }
		  z -= Kshift ;
		  array (a, ira, float) = z ; 
		
		  /* the average and variance discard the NA values */
		  if (0 && dca->isLow) /* 2013_11_18, since we now compute low index with run-dependant plateau, we should include all runs */
		    isLow++ ;
		  else
		    {
		      isHigh++ ;
		      NY++ ;
		      Y += z ;
		      Y2 += z * z ;
		    }
		}
	      arraySort (a, floatOrder) ;
	      for (ira = 1 ; ira <= nRuns ; ira++)  /* number 1 to nRuns, easier for the binomial coef */
		{
		  z = array (a, ira - 1, float) ;
		  if (z > -100)
		    {
		      /* kurtosis = mu4/sigma^4 - 3  see wikipedia 
		       * where muN = 1/n sum(i=1 to n) (xi - xmoyen)^N
		       * and by definition sigma^2 = mu2

		       * if mN = 1/n sum (i =1 to n ) xi^N
		       * we have  
		       *        xmoyen = m1,  mu1 = 0, 
		       *        mu2 = m2 - m1^2
		       *        mu3 = m3 - 3 m1 m2 + 2 m1^3
		       *        mu4 = m4 - 4 m1 m3 + 6 m1^2 m2 - 3 m1^4
		       *
		       * Kurtosis, we must substract 3 for good additivity properties
		       *   define the Nth-cumulan kN
		       *   formula unknown
		       *   the Kurtosis, or standardized 4th cumulant is defined as
		       *     Kurtosis = gamma2 = k4/k2^2 = mu4/mu2^2 - 3
                       *   a good estimator g2 of gamma2, for a sample of size n 
		       *   is given by the formula (wikipedia:kurtosis
		       *     g2 = {(n-1)/(n-2)(n-3)} {(n+1)mu4/mu2^2 - 3 (n-1)}
		       *
		       *
		       *  L moments, first we order the values, then the order statistics of 6,9,3,8
		       *     are xi = 3,6,8,9   counting from 1, i.e. x1=3 x4=9
		       *  In a normal sample of size 6, we have 97% = chance to have the true mean in [x2,x5]
		       *  for 4 point we 2 * 1/2^4 chance of having all point on same size of true mean = 1/8 > 5%
		       *  for 5 point we 2 * 1/2^5 chance of having all point on same size of true mean = 1/16 > 5%
		       *  For size 5 or lower we should use min-max, for 6 use [x2,x5]
		       *     we now define the estimator of N-th L-moment as
		       *    lN = {1/N} { N! (n - N)!/n! } { sum(i=1 to n) [z(i,n;N) xi ]  }
		       *  The important idea is that we do not compute x^2 or x^3... so we are not 
		       *  over influnced by the outliers, rather we use a binomial weight from Pascal triangle
                       *    z(i,n;1) =  1   hence l1 = {1/n} { sum xi } == m1 the conventional average
		       *    z(i,n;2) =  C(1,i-1) - C(1,n-i)     where C(a,b) is the binomial b!/a! (b-a)!
		       *    z(i,n;3) =  C(2,i-1) - 2 C(1,i-1)C(1,n-i) + C(2,n-i)
		       *    z(i,n;4) =  C(3,i-1) - 3 C(2,i-1)C(1,n-i) + 3 C(1,i-1)C(2,n-i) - C(3,n-i)
		       *
		       *  L-moment ratios tauN = lN/l2 belong to [-1,1]  -> L-skewness and L-kurtosis tau4 is in [-1/4,1]	
		       *  The coefficent of L-variation (i.e. L-CV) is defined as tau1 = l2/l1 belongs to [0,1] == Gini coef
		       */
		      NX++ ; X += z ; X2 += z * z ; X3 += z * z * z ; X4 += z * z * z * z ;
		      if (ira == 1) 
			c2 = 1 - 1/(nRuns - 1.0) ;
		      else if (ira == nRuns) 
			c2 = 1/(nRuns - 1.0) - 1 ;
		      else
			c2 =  1/(ira - 1.0) - 1/(nRuns - ira) ;
		      L2 += c2 * z ;

		      L3 += (C3[0] - 2*C3[1]*C3[2] + C3[3]) * z ;
		      L4 += (C4[0] - 3*C4[1]*C4[2] + 3*C4[3]*C4[4] - C4[5]) * z ;
		    }
		}
	      if (nRuns > 2)
		L2 = L2 / (2.0 * nRuns * (nRuns - 2.0)) ;

	      ia = rc->nRuns - 1 ;
	      dc->indexMin = array (a, ia/10, float) ;
	      dc->indexMax = array (a, ia - ia/10, float) ;
	      
	      if (NX > 1)
		{ X /= NX ; X2 /= NX ; }
	      if (NY > 1)
		{ Y /= NY ; Y2 /= NY ; }

	      /*
		variance is now computed in groupCumul
		using a numerically stable algorithm
	      */
	      if (NX > 1 && isHigh > isLow)
		{
		  z = X2 - X*X ;
		  dc->divariance = z > 0.00001 ? sqrt (z) : 0 ;
		}
	      else
		dc->divariance = -999 ;
	      dc->N = NX ;

	      if (dc->index < dc->indexMin) dc->indexMin = dc->index ;
	      if (dc->index > dc->indexMax) dc->indexMax = dc->index ;

	      /* look for a dromadaire
	       * find the point which splits the data in 2 parts of minimal combined variance
	       *  retain this minimal variance
	       */
	      if (1 && NX > 20 && 100 * NX >= 95 * rc->nRuns && isHigh > isLow)
		{
		  float split, av, sigma, av1, sigma1, av2, sigma2, sigmaCombined ;

		  diVariance (a, &split, &av, &sigma, &av1, &sigma1, &av2, &sigma2, &sigmaCombined) ;
		  dc->divariance = sigma > 0 ? sigmaCombined * sigmaCombined/sigma : 0 ;
		}
	    }
      }
  fprintf (stderr, "// %s computed %d index for %d genes in %d runs", 
	   timeShowNow (), nn, dictMax (gx->geneDict), dictMax (runDict)) ;
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, ", %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
  }

  arrayDestroy (a) ;
  arrayDestroy (indexes) ;
  arrayDestroy (lns) ;
  return nn ;  
} /* gxComputeAllIndex */

/*************************************************************************************/
/* given the counts, compute the average index of each gene across all runs */
static int gxComputeAverageIndex (GX *gx)
{
  int run, gene, nn = 0 ;
  DC *dc, *dca ;
  RC *rc ;
  GC *gc ;
  Array a, aa ;
  int runMax = arrayMax (gx->runs) ;
  int geneMax = arrayMax (gx->genes) ;

  a = arrayCreate (dictMax (gx->geneDict), DC) ;

  for (run = 1, rc = arrayp (gx->runs, run, RC) ; run < runMax ; rc++, run++) 
    {
      if (! rc->indexDone) continue ; 
      if (rc->runs)
	continue ;
      aa = rc->aa ;
      if (! aa) 
	continue ;
      nn++ ;
      for (gene = 1, dc = arrayp (aa, 1, DC) ; gene < arrayMax(aa) ; gene++, dc++) 
	{
	  GC *gc = arrayp (gx->genes, gene, GC) ; 
	  if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
	    continue ;
	  dca = arrayp (a, gene, DC) ;
	  dca->index += dc->index ;
	}
    }

  if (nn > 0)
    for (gene = 1, dca = arrayp (a, 1, DC), gc =  arrayp (gx->genes, gene, GC) ; gene < geneMax ; gene++, dca++, gc++) 
      {
	if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
	  continue ;
	gc->averageIndex = dca->index / nn ;
      }

  arrayDestroy (a) ;
  return nn ;  
} /* gxComputeAverageIndex */

/*************************************************************************************/
/* parse list of gene groups: groupName, followed on the same line by gene names 
 *  to exploit the macro-genes created at the Broad
 */
static int gxParseGeneGroup (GX *gx, const char* fileName)
{
  AC_HANDLE h = ac_new_handle () ; 
  ACEIN ai ; 
  const char *ccp ;
  int n = 0 ;
  BOOL inGG = FALSE ;
  GC *gc, *myGgg ;
  KEYSET ks ;
  
  ai = aceInCreate (fileName, FALSE, h) ;
  if (! ai)
    messcrash ("cannot open input file -geneGroup %s", fileName ? fileName : "stdin") ;
  aceInSpecial (ai, "\"\n\t") ;
  
  while (aceInCard (ai)) 
    if ((ccp = aceInWord (ai)))
      {
	if (! ccp || ! *ccp)
	  inGG = FALSE ;
	if (! strcasecmp (ccp, "Gene_group"))
	  {
	    int group ;
	    ccp = aceInWord (ai) ;
	    if (! ccp || ! *ccp)
	      continue ;
	    inGG = TRUE ;
	    dictAdd (gx->geneDict, ccp, &group) ;
	    gc = arrayp (gx->genes, group, GC) ;
	    if (! gc->Group_level)
	      gc->Group_level = 1 ;
	    if (! gx->maxGeneGroupLevel)
	      gx->maxGeneGroupLevel = 1 ;
	    myGgg = gc ;
	    if (! gc->geneGroup)
	      ks = gc->geneGroup = arrayHandleCreate (32, GGG, gx->h) ;
	    n = 0 ;
	    gc->isGood = TRUE ;
	    continue ;
	  }
	
	if (! inGG)
	  continue ;
	else if (! strcasecmp (ccp, "Title"))
	  {
	    ccp = aceInWord (ai) ;
	    if (ccp) 
	      {
		gc->title = stackMark (gx->info) ;
		pushText (gx->info, ac_unprotect (ccp,h)) ;
	      }
	  }
	else if (! strcasecmp (ccp, "GroupId"))
	  {
	    ccp = aceInWord (ai) ;
	    if (ccp) 
	      {
		gc->geneId = stackMark (gx->info) ;
		pushText (gx->info, ac_unprotect (ccp,h)) ;
	      }
	  }
	else if (! strcasecmp (ccp, "HGNC_family"))
	  {
	    ccp = aceInWord (ai) ;
	    if (ccp) 
	      {
		gc->HGNC_family = stackMark (gx->info) ;
		pushText (gx->info, ac_unprotect (ccp,h)) ;
	      }
	  }
	else if (! strcasecmp (ccp, "Genes"))
	  {
	    char *cp = aceInWord (ai) ;
	    ccp = gxSplitMrnaAlias (gx, cp) ;
	    if (ccp) 
	      {
		GGG *ggg ;
		int gene ;
		float w = 1 ;
		int n1 = arrayMax (ks) ;
		char *cp = messprintf ("X__%s", ccp) ;
		n++ ;
		if (! dictFind (gx->geneDict, cp, &gene))
		  dictAdd (gx->geneDict, ccp, &gene) ;
		ggg = arrayp (ks, n1, GGG) ;
		ggg->gene = gene ;
		ggg->weight = w = 1 ;
		aceInNext (ai) ;
		if (aceInFloat (ai, &w))
		  {
		    ggg->weight = w ;
		    myGgg->hasWeights = TRUE ; 
		  }
	      }
	  }
 	else if (! strcasecmp (ccp, "Group_level"))
	  {
	    int level = 0 ;
	    aceInNext (ai) ;
	    if (aceInInt (ai, &level))
	      {
		gc->Group_level = level ;
		if (gx->maxGeneGroupLevel < level)
		  gx->maxGeneGroupLevel = level ;
	      }
	  }
      }
  ac_free (h) ;

  return n ;
} /* gxParseGeneGroup  */

/*************************************************************************************/
static void gxPrepareGeneClusters (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai ;
  Array aa ;
  GG *gg ;
  GC *gc ;
  int nn = 0, gene, chrom , a1, a2, ngg1 = 0, ngg2 = 0 ;
  const char *ccp ;

  aa = gx->geneClusters = arrayHandleCreate (1000, GG, gx->h) ;
  ai = aceInCreate (gx->geneClustersFileName, FALSE, h) ;
  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (! ccp)
	continue ;
      dictAdd (gx->chromDict, ccp, &chrom) ;
      a1 = a2 = 0 ;
      aceInStep (ai, '\t') ; aceInInt(ai, &a1) ;
      aceInStep (ai, '\t') ; aceInInt(ai, &a2) ;
      if (a1 && a2)
	{
	  gg = arrayp (aa, nn++, GG) ;
	  gg->chrom = chrom ;
	  gg->a1 = a1 ;
	  gg->a2 = a2 ;
	  gg->genes = keySetHandleCreate (gx->h) ;
	  ngg1++ ;
	}
    }

  for (gene = 1, gc = arrp (gx->genes, gene, GC) ; gene < dictMax(gx->geneDict) ; gc++, gene++)
    {
      int i, chrom = gc->chrom, a1 = gc->a1 ;
      if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
	continue ;
      for (i = 0, gg = arrp (aa, i, GG) ; i < nn ; i++, gg++)
	if (gg->chrom == chrom && gg->a1 < a1 && gg->a2 > a1)
	  {
	    int j = keySetMax (gg->genes) ;
	    keySet (gg->genes, j) = gene ;
	    ngg2++ ;
	    break ;
	  }
    }

  fprintf (stderr, "Constructed %d gene groups, using a total of %d genes\n", ngg1, ngg2) ;
  ac_free (h) ;
} /* gxPrepareGeneClusters */

/*************************************************************************************/

static void gxComputeGeneClusterIndex (GX *gx)
{
  int i, ii, gene, run, ng ;
  DC *dc, *dcg ;
  GG *gg ;
  RC *rc ;
  Array aa, aag ;
  int geneClusterMax = arrayMax(gx->geneClusters) ;
  int runMax = arrayMax (gx->runs) ;

  for (run = 1, rc = arrayp (gx->runs, run, RC) ;  geneClusterMax > 0 && run < runMax ; rc++, run++) 
    {
      if (! rc->indexDone) continue ; 
      if (rc->runs)
	continue ;
      aa = rc->aa ;
      if (! aa) 
	continue ;
      aag = rc->aag ;
      if (! aag)
	aag = rc->aag = arrayHandleCreate (geneClusterMax, RC, gx->h) ;

      for (ii = 0, gg = arrp(gx->geneClusters, 0, GG) ; ii < geneClusterMax ; ii++, gg++)
	{
	  dcg = arrayp (aag, ii, DC) ;
	  memset (dcg, 0, sizeof(DC)) ;
	  ng = gg->genes ? keySetMax (gg->genes) : 0 ;
	  if (ng > 0)
	    {
	      for (i = 0 ; i < ng ; i++)
		{
		  gene = keySet (gg->genes, i) ;
		  dc = arrayp (aa, gene, DC) ;
		  dcg->index += dc->index ;
		  dcg->seqs += dc->seqs ;
		  dcg->tags += dc->tags ;
		}
	      dcg->index /= ng ;
	    }
	}
    }
} /* gxComputeGeneClusterIndex  */

/*************************************************************************************/

static void gxExportGeneClusterIndex (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, run ;
  GG *gg ;
  DC *dc ;
  RC *rc ;
  int geneClusterMax = arrayMax(gx->geneClusters) ;
  int runMax = arrayMax (gx->runs) ;
  ACEOUT ao = aceOutCreate (gx->outFileName, ".geneClusters.ace",  gx->gzo, h) ;
  aceOutDate (ao, "//", gx->title) ;

  for (ii = 0, gg = arrp (gx->geneClusters, 0, GG)  ; ii < geneClusterMax ; gg++, ii++)
    {
      BOOL ok = FALSE ;

      for (run = 1, rc = arrayp (gx->runs, run, RC) ; run < runMax ; rc++, run++) 
	{
	  if (rc->private || ! rc->aag || ii >= arrayMax (rc->aag))
	    continue ;

	  dc = arrp(rc->aag, ii, DC) ;
	  if (! dc->tags)
	    continue ;

	  if (! ok)
	    {
	      ok = TRUE ;
	      aceOutf (ao, "Gene Zone%s__%d_%d // ng=%d\n"
		       , dictName (gx->chromDict, gg->chrom)
		       , gg->a1, gg->a2
		       , gg->genes ? keySetMax(gg->genes) : 0
		       ) ;
	    }

	  aceOutf (ao, "Run_U %s  %.2f %.2f seqs %.2f tags %.2f kb\n"
		   , dictName (gx->runDict, run)
		   , dc->index - 1000
		   , dc->seqs, dc->tags, dc->kb
		   ) ;
	}
      if (ok)
	aceOut (ao, "\n") ;
    }
  ac_free (h) ;
} /* gxExportGeneClusterIndex  */

/*************************************************************************************/
/*************************************************************************************/
/* given the index in all runs, 
 * construct the sum of all runs and compute its index (done in gxGlobalCout
 * then compute the linear regression of each run towards the any run, limited to well expressed genes with low variance
 * then shift each run accordingly and iterate
 */

/* select the most stable well expressed genes */
static float gxFineTuneMaxVariance (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, gene ;
  float maxVar = 0 ;
  DC *dc ;
  RC *rc = arrayp (gx->runs, gx->runAny, RC) ;
  Array aa = rc->aa 
    , vv = arrayHandleCreate (dictMax (gx->geneDict) + 1, float, h) ;

  if (gx->stableGenesFileName)
    {
      ACEIN ai = aceInCreate (gx->stableGenesFileName, FALSE, h) ;
      const char *ccp ;

      gx->lowVarianceGenes = keySetHandleCreate (gx->h) ;
      while (aceInCard (ai))
	{
	  char *cp = aceInWord (ai) ;
	  if (! ccp || *ccp == '#') continue ;
	  ccp = gxSplitMrnaAlias (gx, cp) ;
	  dictAdd (gx->geneDict, ccp, &gene) ;
	  keySet (gx->lowVarianceGenes, nn++) = gene ;
	}
    }

  else
    {
      for (gene = 1, dc = arrayp (aa, 1, DC) ; gene <= dictMax (gx->geneDict) ; gene++, dc++)
	if (dc->index >= 1015 && dc->variance >= 0)
	  array (vv, nn++, float) = dc->variance ;
      arraySort (vv, floatOrder) ;
      maxVar = arr (vv, arrayMax (vv)/3, float) ;
      
      gx->lowVarianceGenes = keySetHandleCreate (gx->h) ;
      for (nn = 0, gene = 1, dc = arrayp (aa, 1, DC) ; gene <= dictMax (gx->geneDict) ; gene++, dc++)
	if (dc->index >= 1015 && dc->variance <= maxVar && dc->variance >= 0)
	  keySet (gx->lowVarianceGenes, nn++) = gene ;
    }

  keySetSort (gx->lowVarianceGenes) ;
  ac_free (h) ;
  return maxVar ;
} /* gxFineTuneMaxVariance */

/**************************/
/* adjust one run, on genes with index > 15 and variance < maxVar */
static void gxFineTuneOneRun (GX *gx, int pass, int run0, int run, float maxVar)
{
  AC_HANDLE h = ac_new_handle () ;
  RC *rc0 = arrayp (gx->runs, run0, RC) 
    , *rc = arrayp (gx->runs, run, RC) ;
  DC *dc0, *dc ;
  POINT2D *pp ;
  Array aa0 = rc0->aa, aa = rc->aa ;
  Array ww = arrayHandleCreate (dictMax (gx->geneDict) + 1, POINT2D, h) ;
  ACEOUT ao = gx->aoFineTune ;
  int nn = 0, gene ;
  static int firstPass = 0 ;
  float X, Y ;

  if (! ao)
    {
      ao = gx->aoFineTune = aceOutCreate (gx->outFileName, ".fineTune.txt", FALSE, gx->h) ;
      aceOutDate (ao, "###", gx->title) ;
      aceOut (ao, "## Linear regression towards the global run: run(gene) = a*global(gene) + b ;  r,w = usual stats\n") ;
      aceOut (ao, "#Pass\tRun\tRunId\tTitle\tGenes\ta\tb\tdb\tr\tw\n") ;
    }
  if (rc->runs)
    rc->indexDone = 1 ; /* force future recalculation of the variance */
  if (! rc->aa || run == run0)
      goto done ;

  
  if (0) /* adjust the index to have 20,000 genes at index 10, or set 10k genes at index 13 */
    {
      /* this method does not work, testes 2013_02_17 on Neuroblastoma data 
       * the AUC predictors on the test set  all fell by 0.1 to .4 
       * not much, but significant
       */
      KEYSET histo = keySetHandleCreate (h) ;
      int i, n ;
      float dx ;

      /* histo of the index */
      for (gene = 1,  dc = arrayp (aa, 1, DC) ; gene <= dictMax (gx->geneDict) ; gene++, dc++)
	{
	  GC *gc = arrayp (gx->genes, gene, GC) ; 
	  if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
	    continue ;
	  
	  n = 100.0*(dc->index - 1000.0) ; if (n < 0) n = 0 ;
	  keySet (histo,n)++ ;
	}
      /* cumulated histo */
      for (n = 0, i = keySetMax (histo) - 1 ; i>=0 ; i--)
	{
	  n += keySet (histo, i) ;
	  if (n >= 10000) break ;
	}
      if (i > 1000 && i < 1600) /* adjust but not too much */
	{
	  dx = (1300 - i)/100.0 ;
	  aceOutf (ao, "%s\t%.1f\n", dictName (gx->runDict, run), dx) ;
	  for (gene = 1,  dc = arrayp (aa, 1, DC) ; gene <= dictMax (gx->geneDict) ; gene++, dc++)
	     {
	       GC *gc = arrayp (gx->genes, gene, GC) ; 
	       if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
		 continue ;
	       dc->index += dx ;
	       if (dc->index < 0) 
		 dc->index = 0 ;
	     }
	}
    }

  else
    {
      /* construct a dot plot of all good genes */

      for (X = Y = 0, nn = 0, gene = 1, dc0 = arrayp (aa0, 1, DC) , dc = arrayp (aa, 1, DC) ; gene <= dictMax (gx->geneDict) ; gene++, dc0++, dc++)
	if (dc->index >= 1015 && dc->index > rc->zeroIndex + 2 && keySetFind (gx->lowVarianceGenes, gene, 0))
	  {
	    pp = arrayp( ww, nn++, POINT2D) ; 
	    pp->x = dc0->index - 1000 ; pp->y = dc->index + rc->fineTune - 1000 ;
	    X += pp->x ; Y += pp->y ;
	    if (0) { pp->x = pp->y = nn ; }
	    if (nn < 10 && firstPass < 1)
	      fprintf (stderr, "%s\t%s\t%s\t%s\t%s\t%.1f\t%.1f\n"
		       , dictName (gx->runDict, run0)
		       , dictName (gx->runDict, run)
		       , rc->runid ?  stackText (gx->info, rc->runid) : dictName (gx->runDict, run)
		       , rc->title ?  stackText (gx->info, rc->title) : dictName (gx->runDict, run)
		       , dictName (gx->geneDict, gene)
		       , pp->x
		       , pp->y
		       ) ;
	  }
      aceOutf (ao, "Using %d genes: ", nn) ;
  
      /* adjust */
      {
	double a = 1, b = 0, r = 0, w = 0 ;
	linearRegression (ww, &a, &b, &r, &w) ;
	b = nn ? (Y - X)/nn : 0 ;
	if (1) rc->fineTune -= b ;    /* dromadaire */
	aceOutf (ao, "%d\t%s\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\n"
		 , pass, dictName (gx->runDict, run)
		 , rc->runid ?  stackText (gx->info, rc->runid) : dictName (gx->runDict, run)
		 , rc->title ?  stackText (gx->info, rc->title) : dictName (gx->runDict, run)
		 , nn, a, rc->fineTune, -b, r, w
		 ) ;
	if (firstPass < 1)
	  fprintf (stderr, "%d\t%s\t%d\t%f\t%f\t%f\t%f\n"
		   , pass, dictName (gx->runDict, run), nn, a, b, r, w
		   ) ;
      }
    }
  firstPass++ ;
 done:
  ac_free (h) ;
  return ;
} /* gxFineTuneOneRun */

/**************************/

static void gxFineTune (GX *gx)
{
  int pass, run, runMax = dictMax (gx->runDict) + 1 ;
  float maxVar = 0 ;

  for (pass = 0 ; runMax > 2 && pass < 3 ; pass++)
    {
      fprintf (stderr, "// %s: fineTune pass %d\n", timeShowNow (), pass) ;
      maxVar = gxFineTuneMaxVariance (gx) ;
      for (run = 1 ; run < runMax ; run++) 
	gxFineTuneOneRun (gx, pass, gx->runAny, run, maxVar) ;
      gxComputeAllIndex (gx) ; /* will just recompute the variance of gx->runAny */
    }
  fprintf (stderr, "// %s: fineTune done\n", timeShowNow ()) ;
  return ;
} /* gxFineTune */

/*************************************************************************************/
/*************************************************************************************/
/* test the rankGene method */

typedef struct rankStruct { int gene ; double index, Z ; } RANK ;

static int rankGeneOrder (const void *a, const void *b)
{
  const RANK *up = (const RANK *)a, *vp = (const RANK *)b ;
  int n ;

  n = up->gene - vp->gene ; if (n) return n ;
  return 0 ;
} /* rankGeneOrder */


static int rankZOrder (const void *a, const void *b)
{
  const RANK *up = (const RANK *)a, *vp = (const RANK *)b ;
  float x ;

  x = up->Z - vp->Z ; return (x > 0) ? 1 : (x < 0 ? -1 : 0) ;

} /* rankZOrder */


static void gxRankGene (GX *gx)
{
  int gene, run, runMax = dictMax (gx->runDict) + 1 ;
  double Z = 0 ;
  Array aa, bb = 0 ;
  DC *dc ;
  RC *rc ;
  RANK *rank ;

  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
    {
      aa = rc->aa ;
      if (! aa || rc->runs)
	continue ;
      arrayDestroy (bb) ;
      /* copy and number the lines */
      bb = arrayCreate (arrayMax (aa), RANK) ;
      for (gene = arrayMax(aa) - 1, dc = arrayp (aa, gene, DC), rank = arrayp (bb, gene, RANK) ; gene >= 0 ; dc--, rank--, gene--)
	{
	  rank->gene = gene ;
	  rank->index = dc->index ;
	}
      /* sort */
      arraySort (bb, rankZOrder) ;
      /* compute the Z */
      for (gene = arrayMax(aa) - 1, Z = 0, rank = arrayp (bb, gene, RANK) ; gene >= 0 ; gene--, rank--)
	{
	  Z += rank->index ;	  
	  rank->Z = Z ;
	}
      arraySort (bb, rankGeneOrder) ;
      for (gene = arrayMax(aa) - 1,  dc = arrayp (aa, gene, DC), rank = arrayp (bb, gene, RANK) ; gene >= 0 ; gene--, dc--, rank--)
	{
	  dc->index = 1000 + 20.0 * (1 - rank->Z/Z) ;
	}      
    }

  arrayDestroy (bb) ;

  fprintf (stderr, "// %s: rankGene done\n", timeShowNow ()) ;
  return ;
} /* gxRankGene */

/*************************************************************************************/
/*************************************************************************************/
/* construct a group with 20% of the sample where seedGene is most extreme */
static void gxConstructSeedGroup (ACEOUT ao, GX *gx, Array pp, int runa, int runb)
{
  RC *rc ;
  PC *pc ;
  int i, nRuns, run ;

  run = runa ;
  rc = arrayp (gx->runs, run, RC) ;
  nRuns = arrayMax (pp) ;
  rc->runs = keySetHandleCreate (gx->h) ;
  for (i = 0, pc = arrayp (pp, arrayMax(pp)-1, PC); i < nRuns/5 ; i++, pc--)
    {
      if (pc->gene)
	{
	  keySet (rc->runs, i) = pc->gene ; /* horrible hack */  
	  aceOutf (ao, "%s\t%s\t%.2f\n", dictName (gx->runDict, keySet (rc->runs,i)), dictName (gx->runDict, runa), pc->dIndex) ;
	}
    }
  keySetSort (rc->runs) ;
  keySetCompress (rc->runs) ;
  rc->nRuns = keySetMax (rc->runs) ;
  rc->aa = arrayHandleCreate (50000, DC, gx->h) ;   
  rc->compare_to = arrayHandleCreate (50000, DC, gx->h) ; 
  keySet (rc->compare_to, 0) = runb ;
  
  /* register the runs of my seedGroup_high seedGroup_low */
  run = runb ;
  rc = arrayp (gx->runs, run, RC) ;
  nRuns = arrayMax (pp) ;
  rc->runs = keySetHandleCreate (gx->h) ;
  for (i = 0, pc = arrayp (pp, 0, PC); i < nRuns/5 ; i++, pc++)
    {
      if (pc->gene)
	{
	  keySet (rc->runs, i) = pc->gene ; /* horrible hack */  
	  aceOutf (ao, "%s\t%s\t%.2f\n", dictName (gx->runDict, keySet (rc->runs,i)), dictName (gx->runDict, runb), pc->dIndex) ;
	}
    }
  keySetSort (rc->runs) ;
  keySetCompress (rc->runs) ;
  rc->nRuns = keySetMax (rc->runs) ;
  rc->aa = arrayHandleCreate (50000, DC, gx->h) ; 
  rc->compare_to = arrayHandleCreate (50000, DC, gx->h) ; 
  rc->compare_to = keySetHandleCreate (gx->h) ;
  keySet (rc->compare_to, 0) = runa ;
  
} /* gxConstructSeedGroup */

/*************************************************************************************/
/* construct a group with 20% of the sample where seedGene is most extreme */
static void gxSeedGroup (GX *gx, int pass)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = 0 ;
  int gene, run, runa, runb, runMax = arrayMax (gx->runs), nRuns ;
  RC *rc ;
  DC *dc ;
  PC *pc ;
  Array aa, pp ;
  
  if (dictFind (gx->geneDict, gx->seedGene, &gene))
    {
      ao = aceOutCreate (gx->outFileName, messprintf (".SG%d_seededOn_%s.txt", pass, gx->seedGene), FALSE, h) ;
      aceOutDate (ao, "###", gx->title) ;
      aceOutf (ao, "## 20%% runs with highest then lowest expression of the %s\n", gx->seedGene) ;
      aceOut (ao, "# Run\tLow/High\tIndex\n") ;
      pp = arrayHandleCreate (dictMax (gx->geneDict), PC, h) ;
      for (nRuns = 0, run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	{
	  if (rc->runs || rc->private)
	    continue ;
	  aa = rc->aa ;
	  if (! aa || rc->runs)
	    continue ;   
	  dc = arrayp (aa, gene, DC) ;
	  {
	    pc = arrayp (pp, nRuns++, PC) ;
	    pc->gene = run ;         /* horrible hack , i am lazy, i prefer to reuse the PC type */
	    pc->dIndex = dc->index > 200 ? dc->index + rc->fineTune - 1000 : dc->index - 1000 ; /* horrible hack */
	  }
	}
      arraySort (pp, pcOrder) ;    /* horrible hack */  

      /* register the runs of my seedGroup_high seedGroup_low */
      dictAdd (gx->runDict, messprintf("SG%d_%s_high", pass, gx->seedGene),&runb) ;
      dictAdd (gx->runDict, messprintf("SG%d_%s_low", pass, gx->seedGene),&runa) ;

      gxConstructSeedGroup (ao, gx, pp, runa, runb) ;
    }
  ac_free (h) ;
  gxComputeAllIndex (gx) ;  /* will just do the newly created groups */  
  return ;
} /* gxSeedGroup */


/*************************************************************************************/
/*************************************************************************************/
/* for each quantitative phenotype available per run, probably via run->sample
 * int the training phase
 * for each gene
 *   take that set of runs
 *   sort them according to each gene vai the index of the gene in that run
 *   take the top and bottom quarter of the run list
 *   average the phenotype in that quarter, or maybe take the median
 *   or may be multiply the delta phenotype around the median phenotype
 *   by the delta gene index for that run around the median index of the whole population
 *   this gives a differential signal for the gene
 *   the original idea was to take the difference of the area under the Kaplan Mayer curves
 *   but it may be the same
 *
 * in the prediction phase,
 *   for each run
 *   we sum the (centered) gene expression times the gene phenoIndex
 *   and this yield a predicted phenotype value for that run
 *   we classify the runs by this parameter and compute the true kaplan mayer of the top and low quarters
 */
 
typedef struct phenoStruct { int run, pheno, unit ; float z ; } PHENO ;
#ifdef JUNK
static void gxPhenoGenesOne (GX *gx, int pheno)
{
  AC_HANDLE h = ac_new_handle () ;

  /* get a list of runs for that pheno with their pheno value 
   * for each gene sort those runs
   * sort a struct (run, pheno, geneindex) by index
   * compute the gene index in a struct (gene, index) by pheno
   * use that in predictions
   */

  ac_free (h) ;
  return ;
} /* gxPhenoGenesOne */
#endif

static void gxPhenoGenes (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  int ir, nn = 0 ; 
  int runMax = arrayMax (gx->runs) ;
  int run, pheno, unit ;
  AC_OBJ Run, Sample ;
  AC_TABLE tbl ;
  DICT *phenoDict = dictHandleCreate (64, h) ;
  PHENO *ph ;
  Array phenos = arrayHandleCreate (10 * runMax, PHENO, h) ;
  RC *rc ;

  if (gx->db)
    {
      for (run = 1, rc = arrp (gx->runs, run, RC) ; run < runMax - 1 ; rc++, run++)
	{
	  Run = ac_get_obj (gx->db, "Run",  dictName (gx->runDict, run), h) ;
	  if (Run)
	    {
	      AC_HANDLE h1 = ac_new_handle () ;
	      Sample = ac_tag_obj (Run, "Sample", h1) ;
	      if (Sample)
		{
		  tbl = ac_tag_table (Sample, "Phenotype", h) ;
		  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
		    {
		      float z = ac_table_float (tbl, ir, 2, -999999) ;
		      if (z > -999998) /* found a value */
			{
			  const char *tag = ac_table_printable (tbl, ir, 0, "X") ;
			  dictAdd (phenoDict, tag, &pheno) ;
			  tag = ac_table_printable (tbl, ir, 1, "X") ;
			  dictAdd (phenoDict, tag, &unit) ;
			  ph = arrayp (phenos, nn++, PHENO) ;
			  ph->run = run ; ph->pheno = pheno ; ph->unit = unit ; ph->z = z ;
			}
		    }
		}
	      ac_free (h1) ;
	    }
	}
    }
  
  /* gxDrawKaplanMayer (gx, phenos) ; */
  
  ac_free (h) ;
  return ;
} /* gxPhenoGenes */

/*************************************************************************************/
/*************************************************************************************/
/* compute for each run the minimal index */
static void gxAllGeneCorrelations (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  RC *rc ;
  int i, j, run, runMax = arrayMax (gx->runs) ;
  int gMax = dictMax (gx->geneDict) + 1 ;
  int ok, N ;
  float minLevel = 0 ;
  double X, Y, X2, Y2, XY ;
  ACEOUT ao = aceOutCreate (gx->outFileName, ".gene_correl.txt", gx->gzo, h) ;
  DICT *dict = gx->geneDict ;
  
  aceOutDate (ao, "###", "Correlation of all expressed genes") ;
  aceOutf (ao, "# Gene 1\tGene 2\tCorrelation coefficient <xy>/sqrt(<x2><y2>)\n") ;
  
  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++)
    if (minLevel < rc->NA_crossover_index)
      minLevel = rc->NA_crossover_index ;
  minLevel += 1 ;

  if (0 && gMax > 1000) gMax = 1000 ;

  for (i = 1 ; i < gMax ; i++)
    {
      GC *gc1 = arrp (gx->genes, i, GC) ;
      if (1)
	if (gx->captured && gc1->captured)
	  continue ;
      if (gx->targeted && ! gc1->targeted)
	continue ;
		  
      if (gc1->maxIndex + 1000 < 10) 
	continue ;
      if (0 && strcmp (dictName(dict, i), "X__ACTBL2"))
	continue ;
      for (j = i+1 ; j < gMax ; j++)
	{
	  GC *gc2 = arrp (gx->genes, j, GC) ;
	  if (gc2->maxIndex + 1000 < 10) 
	    continue ;
	  if (gx->targeted && ! gc2->targeted)
	    continue ;
	  if (1)
	    if (gx->captured && ! gc2->captured)
	      continue ;
	  if (0 && strcmp (dictName(dict, j), "X__ACTR3BP2"))
	    continue ;
	  ok = N = 0 ; X = X2 = Y =Y2 = XY = 0 ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      double z1, z2 ;
	      DC *dc1, *dc2 ;
	      Array aa = rc->aa ;
	      
	      if (! aa || j >= arrayMax (aa))
		continue ;
	      
	      dc1 = arrayp (aa, i, DC) ;
	      dc2 = arrayp (aa, j, DC) ;
	      
	      z1 = dc1->index - 1000 ;
	      z2 = dc2->index - 1000 ;
	      if (dc1->isLow || z1 < minLevel) z1 =  minLevel ;
	      if (dc2->isLow || z2 < minLevel) z2 =  minLevel ;
	      if (1)
		{
		  N++ ;
		  X += z1 ; X2 += z1 * z1 ;
		  Y += z2 ; Y2 += z2 * z2 ; XY += z1 * z2 ;
		}
	      if (z1 + z2 > 2 * minLevel + 1)
		ok++ ;
	    }
	  
	  if (ok > 10) /* number of expressed values */
	    {
	      double x = X2/N - X*X/(N * N) ;
	      double y = Y2/N - Y*Y/(N * N) ;
	      double z = XY/N - X*Y/(N * N) ;
	      if (x * y == 0) {x = y = z = 1 ; }
	      if (1)
		{
		  z = z / sqrt (x*y) ;
		  if (z > .4 || z < -.4) /* export */ 
		    {
		      aceOutf (ao, "%s\t%s\t%.2f\n"
			       , dictName (dict, i) 
			       , dictName (dict, j) 
			       , z
			       ) ;
		      /* symmetrize the matrix */
		      aceOutf (ao, "%s\t%s\t%.2f\n"
			       , dictName (dict, j) 
			       , dictName (dict, i) 
			       , z
			       ) ;
		      if(0) aceOutf (ao, "N=%d X=%f X2=%f Y=%f Y2=%f XY=%f\n", N,X,X2,Y,Y2,XY) ;
		    }
		}
	    }
	}
    }

  ac_free (h) ;
} /*  gxAllGeneCorrelations */

/*************************************************************************************/
/*************************************************************************************/
/* compute for each run the minimal index */
static void gxNumberOfGenesPerRunAboveGivenIndexDo (GX *gx, BOOL show, BOOL isCapture)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, j3, run, gene, n, n2, n3, ngid, maxIndex = 0, pass, pass2 ;
  Array aa, hhh, hhh2, hhhMA, hhhGID, hh = 0, hh2 = 0, hh3 = 0, hh4 = 0, hhMA = 0, hhGID = 0 ;
  DC *dc ;
  GC *gc ;
  RC *rc ;
  int runMax = arrayMax (gx->runs) ;
  
  hhh = arrayHandleCreate (100 * arrayMax (gx->runs), int, h) ;
  hhh2 = arrayHandleCreate (100 * arrayMax (gx->runs), int, h) ;
  hhhMA = arrayHandleCreate (100 * arrayMax (gx->runs), int, h) ;
  hhhGID = arrayHandleCreate (100 * arrayMax (gx->runs), int, h) ;

  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
    {
      double lowestIndex = 9999999 ;

      aa = rc->aa ;
      if (! aa || ! arrayMax (aa))
	continue ;
      if (! show && 
	  (
	   (! isCapture && rc->Genes_with_index) ||
	   (isCapture && rc->captured_genes_with_index) ||
	   (isCapture && ! rc->capture)
	   )
	  )
	continue ;

      hh = arrayReCreate (hh, 100, int) ;
      hh2 = arrayReCreate (hh2, 100, int) ;
      hh3 = arrayReCreate (hh3, 5000, int) ;
      hh4 = arrayReCreate (hh4, 5000, int) ;
      hhMA = arrayReCreate (hhMA, 5000, int) ;
      hhGID = arrayReCreate (hhGID, 5000, int) ;

      /* compute the histogram of index per run */
      for (gene = j3 = 1, dc = arrayp (aa, 1, DC), gc = arrayp (gx->genes, gene, GC) ; gene < arrayMax (aa) ; gene++, dc++, gc++) 
	if (! gc->geneGroup && (gx->keepIndex || dc->tags > 0))
	  {
	    /* hh : histogram step 1 of all values
	       hh2: histogram step 1 of reliable values
	       hh2: histogram step 1 of reliable values for genes with affy
               hh3: histogram step 1/100 of reliable values
               hh4: histogram step 1/100 of low values
	    */
	    if (gx->targeted && ! gc->targeted)
	      continue ;
	    if (isCapture)
	      {
		BOOL ok = FALSE ;
		if (rc->capture)
		  {
		    KEYSET ks = gc->captures ;
		    int i ;
		    for (i = 0 ; ks && ! ok && i < keySetMax (ks) ; i++)
		      if (rc->capture == keySet (ks, i))
			ok = TRUE ;
		  }
		if (! ok)
		  continue ;
	      }
	    if (0)
	      if (gx->captured && ! gc->captured)
		continue ;
	    
	    n = dc->index + .001 - 1000 ;
	    if (! dc->tags || n<0) n = 0 ;
	    if (gx->isSNP || n > 0)
	      array (hh, n, int)++ ;
	    if (! dc->isLow && (gx->isSNP || n) > 0 && dc->index + .001 - 1000 < lowestIndex)
	      lowestIndex = dc->index + .001 - 1000 ;
	    if (dc->isLow && n >= 1) n = 1 ;
	    if (gx->isSNP || n >= 1)
	      {
		array (hh2, n, int)++ ;
		gc = arrp (gx->genes, gene, GC) ;
		if (gc->affy)
		  array (hhMA, n, int)++ ;
		if (gc->geneId)
		  array (hhGID, n, int)++ ;
	      }
	    n = 100 * (dc->index + .001 - 1000) ;
	    if (! dc->tags || n<0) n = 0 ;
	    
	    if (! dc->isLow)
	      array (hh3, n, int)++ ;
	    else
	      array (hh4, n, int)++ ;
	  }
      
      /* compute the cumulated hh and hh2, histo */
      if (maxIndex < arrayMax (hh))
	maxIndex = arrayMax (hh) ;
      for (n = n2 = n3 = ngid = 0, i = arrayMax (hh) - 1 ; i > 0 ; i--)
	{
	  n += arr (hh, i, int) ;
	  arr (hh, i, int) = n ;
	  array (hhh, 100 * run + i, int) = n ;
	  
	  n2 += array (hh2, i, int) ;
	  arr (hh2, i, int) = n2 ;
	  array (hhh2, 100 * run + i, int) = n2 ;
	  
	  n3 += array (hhMA, i, int) ;
	  arr (hhMA, i, int) = n3 ;
	  array (hhhMA, 100 * run + i, int) = n3 ;

	  ngid +=  array (hhGID, i, int) ;
	  arr (hhGID, i, int) = ngid ;
	  array (hhhGID, 100 * run + i, int) = ngid ;
	}
      
      /* compute the 3% quantile = low index */  
      for (n3 = 0, i = 0 ; i < arrayMax (hh3) ; i++)
	n3 += arr (hh3, i, int) ;
      for (n = 0, i = 0 ; i < arrayMax (hh3) ; i++)
	{
	  n += arr (hh3, i, int) ;
	  if (100 * n > 3 * n3)
	    break ;
	}
      if (! isCapture)
	rc->Low_index = i/100.0 ; 
      /* look for the position where there are as many low index above and high index below 
       * to do that we compute the cumul of the NA hh4 from above and the cumul of the not NA hh3 from below
       * and we stop when they cross
       */
      for (n = 0, i = arrayMax (hh4) - 1 ; i >= 0 ; i--)
	{
	  n += arr (hh4, i, int) ;
	  arr (hh4, i, int) = n ;
	}
      for (n = 0, i = 0 ; i < arrayMax (hh3) ; i++)
	{
	  n += arr (hh3, i, int) ;
	  arr (hh3, i, int) = n ;
	  if (n > array(hh4, i, int))
	    break ;	      
	}
      if (! isCapture)
	{
	  rc->NA_crossover_index = i/100.0 ; 
	  
	  rc->Genes_with_index = array (hh, 1, int) ;
	  rc->Genes_with_reliable_index = array (hh2, 2, int) ;
	  rc->Genes_with_index_over_8 = array (hh2, 8, int) ;
	  rc->Genes_with_index_over_10 = array (hh2, 10, int) ;
	  rc->Genes_with_index_over_12 = array (hh2, 12, int) ;
	  rc->Genes_with_index_over_15 = array (hh2, 15, int) ;
	  rc->Genes_with_index_over_18 = array (hh2, 18, int) ;
	  rc->Genes_with_index_over_21 = array (hh2, 21, int) ;
	  rc->Genes_with_index_over_24 = array (hh2, 24, int) ;
	  rc->Genes_with_index_over_27 = array (hh2, 27, int) ;
	  rc->Genes_with_index_over_30 = array (hh2, 30, int) ;
	}
      else
	{
	  rc->captured_genes_with_index = array (hh, 1, int) ;
	  rc->captured_genes_with_reliable_index = array (hh2, 2, int) ;
	  rc->captured_genes_with_index_over_8 = array (hh2, 8, int) ;
	  rc->captured_genes_with_index_over_10 = array (hh2, 10, int) ;
	  rc->captured_genes_with_index_over_12 = array (hh2, 12, int) ;
	  rc->captured_genes_with_index_over_15 = array (hh2, 15, int) ;
	  rc->captured_genes_with_index_over_18 = array (hh2, 18, int) ;
	  rc->captured_genes_with_index_over_21 = array (hh2, 21, int) ;
	  rc->captured_genes_with_index_over_24 = array (hh2, 24, int) ;
	  rc->captured_genes_with_index_over_27 = array (hh2, 27, int) ;
	  rc->captured_genes_with_index_over_30 = array (hh2, 30, int) ;
	}
    }
  if (! isCapture)
    {
      int run, run2, nRuns ;
      RC *rc2 ;

      for (run = 1, rc = arrayp (gx->runs, run, RC) ; run < runMax ; rc++, run++) 
	if (rc->runs && !rc->addCounts)
	  {
	    nRuns = 0 ; rc->NA_crossover_index  = 0 ;
	    for (i = 0 ; i < keySetMax (rc->runs) ; i++)
	      {
		run2 = keySet (rc->runs, i) ;
		rc2 = arrayp (gx->runs, run2, RC) ;
		if (! rc2->aa)
		  continue ;
		nRuns++ ;
		rc->NA_crossover_index += rc2->NA_crossover_index ;
	      }
	    
	    if (nRuns)
	      {
		rc->NA_crossover_index /= nRuns ;
	      }
	  }
    }

   if (!isCapture && show)
    for (pass = 1 ; pass < 4 ; pass++)
      {
	char *suffix = "" ;
	  Array hhhX = 0 ;

	switch (pass)
	  {
	  case 0: 
	    hhhX = hhh2 ;
	    break ;
	  case 1:
	    hhhX = hhh2 ;
	    suffix = "dropLow" ; /* danielle no longer likes the suffix it 2014_02_04 */
	    suffix = "any" ;
	    break ;
	  case 2:
	    hhhX = hhhMA ;
	    suffix = "withMA" ;
	    if (! gx->hasGeneAffy)
	      continue ;
	    break ;
	  case 3:
	    hhhX = hhhGID ;
	    suffix = "withGeneId" ;
	    if (! gx->hasGeneId)
	      continue ;
	    break ;
	  }
	
	{
	  ACEOUT ao ;
	  Array nn1 = arrayHandleCreate (runMax+1, int, h) ;
	  Array nn2 = arrayHandleCreate (runMax+1, int, h) ;
	  
	  if (pass == 0)
	    ao = aceOutCreate (gx->outFileName
			  , ".expression_index_histogram.txt"
			  , FALSE, h) ;
	  else
	    ao = aceOutCreate (gx->outFileName
				    , hprintf (h, ".Expression_index_histograms.%s.txt", suffix)
			       , FALSE, h) ;
	  
	  aceOutDate (ao, "###", gx->title) ;
	  aceOutf (ao, "\n## Table 1 : Identification and alignment statistics of the runs and groups. The first column applies to  the sum of all runs in the project") ;
	  aceOutf (ao, "## Table 2: Histogram of the number of genes at a given index. The index is equal to the log2(sFPKM) + 10), index 10 corresponds to sFPKM 1. The sFPKM ressembles the FPKM with mutiple corrections described in Genome Biology 2015, 16:133\n") ;
	  aceOutf (ao, "## Table 3: Corresponding cumulated histogram\n") ;

	  aceOutf (ao, "\n## Table 1:\n") ; 
	  gxExportTableHeader (gx, ao, 51) ;
	  
	  if (pass && gx->seedGene) 
	    {
	      aceOut (ao, "\nClassifier") ;
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private) continue ;
		  aa = rc->aa ;
		  if (! aa || ! arrayMax (aa))
		    continue ;
		  aceOutf (ao, "\t%.2f", rc->alpha) ;
		}
	    }
	  
	  maxIndex = (maxIndex + 9)/10 ;
	  maxIndex *= 10 ;
	  
	  for (pass2 = 0 ; pass2 < 2 ; pass2++)
	    {
	      aceOutf (ao, "\n## Table %d: %s\n"
		       , 2 + pass2
		       , pass2 == 0 ? "Histogram of the number of genes at a given index" : "Cumulated histogram: number of genes above a given index"
		       ) ;
	      gxExportTableHeader  (gx, ao, 1000) ;
	      for (i = gx->isSNP ? 0 : 1 ; i <= maxIndex ; i++)
		{
		  if (i == 1 && !  gx->isSNP)
		    aceOutf (ao, "\nTouched") ;
		  else if (i >= 2 && !  gx->isSNP)
		    aceOutf (ao, "\n%d", i) ;
		  else if (i == 0 && gx->isSNP)
		    aceOutf (ao, "\n+/+") ;
		  else if (i == 10 && gx->isSNP)
		    aceOutf (ao, "\nm/+") ;
		  else if (i == 20 && gx->isSNP)
		    aceOutf (ao, "\nm/m") ;
		  else
		    aceOutf (ao, "\n%d", i) ;
		  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		    { 
		      int z ;
		      
		      if (rc->private) continue ;
		      aa = rc->aa ;
		      if (! aa || ! arrayMax (aa))
			continue ;
		      z = array (hhhX, 100 * run + i, int)  - (pass2 == 0 ? array (hhhX, 100 * run + i + 1, int) : 0) ;
		      aceOutf (ao, "\t%d", z) ; 
		      if (i >  (gx->isSNP ? 0 : 1 ))
			{
			  z = z  - (pass2 > 0 ? array (hhhX, 100 * run + i + 1, int) : 0) ;
			  array (nn1, run, int) += z ;
			  array (nn2, run, int) += z * i ;
			}
		    }
		}
	      aceOutf (ao, "\n") ;
	    }
	  
	  aceOutf (ao, "\nNumber of expressed genes") ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    if (! rc->private)
	      aceOutf (ao, "\t%d", array (nn1, run, int) ) ;
	  aceOutf (ao, "\nAverage index of these genes") ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    if (! rc->private)
	      aceOutf (ao, "\t%.1f", array (nn2, run, int)/(.0001 + array (nn1, run, int)) ) ;
	  aceOutf (ao, "\n") ;
	}
      }

  arrayDestroy (hh) ;
  arrayDestroy (hh2) ;
  arrayDestroy (hh3) ;
  arrayDestroy (hh4) ;



  ac_free (h) ;
  return ;
} /* gxNumberOfGenesPerRunAboveGivenIndexDo */

static void gxNumberOfGenesPerRunAboveGivenIndex (GX *gx, BOOL show)
{
  gxNumberOfGenesPerRunAboveGivenIndexDo (gx, show, FALSE) ;
  if (gx->hasCapture)
    gxNumberOfGenesPerRunAboveGivenIndexDo (gx, show, TRUE) ;
  return ;
} /* gxNumberOfGenesPerRunAboveGivenIndex */

/*************************************************************************************/
/*************************************************************************************/
/* compute for each run the minimal index */
static KEYSET gxParsePlusMinusFile (GX *gx, BOOL pm)
{
  AC_HANDLE h = ac_new_handle () ;
  KEYSET ks ;
  int k, n = 0 ;
  const char *ccp ;
  const char *fnam = pm ? gx->genePlusFileName : gx->geneMinusFileName ;
  ACEIN ai = aceInCreate (fnam, FALSE, h) ;

  ks = keySetHandleCreate (gx->h) ;
  if (! gx->userGivenGenePlusMinusDict)
    gx->userGivenGenePlusMinusDict = dictHandleCreate (500, gx->h) ;

  while (aceInCard (ai))
    {
      char *cp ;
      while ((cp = aceInWord (ai)))
	{
	  if (strcmp(cp, "Minus") && strcmp(cp, "Plus") && strcmp(cp, "genes:"))
	    {
	      ccp = gxSplitMrnaAlias (gx, cp) ;
	      dictAdd (gx->geneDict, ccp, &k) ;
	      keySet (ks, n++) = k ;
	      dictAdd (gx->userGivenGenePlusMinusDict, ccp, 0) ;
	    }
	}
    }
  if (n > 1)
    {
      keySetSort (ks) ;
      keySetCompress (ks) ;
    }
  fprintf (stderr, "Found %d genes in %s\n"
	   , keySetMax (ks)
	   , fnam
	   ) ;
  ac_free (h) ;
  return ks ;
} /* gxParsePlusMinusFile */

/*************************************************************************************/
/*************************************************************************************/
/* compute for each run the number of titrating genes */
typedef struct ttStruct { int iCompare, rMax, rr[60] ; float shape[60] ; BitSet *bbb ; Array pcs ; Array hh ; } TT ;
/*************************************************************************************/

static int rwOrder (const void *a, const void *b)
{
  const RW *up = (const RW *)a, *vp = (const RW *)b ;
  int n ; double x ;

  x = up->alpha - vp->alpha ; if (x > 0) return 1 ; if (x < 0) return -1 ;
  n = up->run - vp->run ; if (n) return n ;

  return 0 ;
} /* rwOrder */

/*************************************************************************************/

const char *noX__ (const char *ccp, BOOL doIt)
{
  if (doIt && ccp && !strncmp(ccp, "X__", 3)) ccp += 3 ;
  return ccp ;
} /* noX__ */

/*************************************************************************************/
#ifdef JUNK
/*************************************************************************************/

static int ttOrder (const void *a, const void *b)
{
  const TT *up = (const TT *)a, *vp = (const TT *)b ;
  int i, n ;

  for (i = 0 ; i < up->rMax || i < vp->rMax ; i++) 
    {
      n = up->rr[i] - vp->rr[i] ;
      if (n) return n ;
    }
  return 0 ;
} /* ttOrder */

/*************************************************************************************/

static int pcsOrder (const void *a, const void *b)
{
  const PC *up = (const PC *)a, *vp = (const PC *)b ;
  int  n ; double x ;

  x = up->sig - vp->sig ; if (x > 0) return -1 ; if (x < 0) return 1 ;
  n = up->gene - vp->gene ; if (n) return n ;

  return 0 ;
} /* pcsOrder */

/*************************************************************************************/

static void gxTitration (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  int i,j, ii = 0, gene, n, maxGene = dictMax (gx->geneDict), pass, maxTG = 0, nTG ;
  Array aa, tt ;
  BitSet bb ;
  DC *dc ;
  PC *pc ;
  RC *rc ;
  TT *tp ;
  int runMax = arrayMax (gx->runs), rMax = 0 , maxIndex = 0 ;
  COMPARE *compare ;
  int iCompare ;
  int iMax = gx->compares ? keySetMax (gx->compares) : 0 ;

  /* register in a single structure all titration series */
  tt = arrayHandleCreate (64, TT, h) ;
  for (ii = iCompare = 0 ; iCompare < iMax ; iCompare++)
    {
      BOOL ok = TRUE ;
      int run2 = 0 ;
      RC *rc2 ;

      compare = arrp (gx->compares, iCompare, COMPARE) ;
      if (! compare->shape || ! compare->runs)
	continue ;
      if (gx->isINTRON && ! compare->compareINTRON)
	continue ;
      if (gx->isSNP && ! compare->compareSNP)
	 continue ;

      for (i = 0 ; ok && i < keySetMax (compare->runs) && i < 60 ; i++)
	{
	  run2 = keySet (compare->runs, i) ;
	  rc2 = run2 < runMax ? arrayp (gx->runs, run2, RC) : 0 ;
	  if (! rc2 || ! rc2->aa)
	    ok = FALSE ;
	}
      if (! ok)
	continue ;

      tp = arrayp (tt, ii++, TT) ;
      tp->iCompare = iCompare ;
      tp->hh = arrayHandleCreate (32, int, h) ; /* histo of titrating genes */
      tp->pcs = arrayHandleCreate (1000, PC, h) ;
      for (i = 0 ; i < keySetMax (compare->runs) && i < 60 ; i++)
	{
	  tp->rr[i] = keySet (compare->runs, i) ;
	  tp->shape[i] = compare->shape ? array(compare->shape, i, float) : .1 ;
	}
      tp->rMax = i ;
      if (i > rMax) rMax = i ;

      /* create the mirror image */
      tp = arrayp (tt, ii++, TT) ;
      tp->iCompare = - iCompare ;
      tp->hh = arrayHandleCreate (32, int, h) ; /* histo of titrating genes */
      tp->pcs = arrayHandleCreate (1000, PC, h) ;
      tp->rMax = i ;
      for (j = 0 ; j < i ; j++)
	{
	  tp->rr[j] = (tp-1)->rr[j] ;
	  tp->shape[j] = - ((tp-1)->shape[j]) ;
	}
    }
  arraySort (tt, ttOrder) ;
  arrayCompress (tt) ;

  /* search for the titrating genes */ 
  for (ii = 0, tp = arrp (tt, ii, TT) ; ii < arrayMax (tt) ; tp++, ii++)
    {
      double minIndex ;
      for (minIndex = 2, i = 0 ; i < tp->rMax ; i++)
	{
	  rc = arrayp (gx->runs, tp->rr[i], RC) ; 
	  if (rc->zeroIndex > minIndex)
	    minIndex = rc->zeroIndex ;
	  if (rc->NA_crossover_index > minIndex)
	    minIndex = rc->NA_crossover_index ;
	}

      for (nTG = 0, gene = 1 ; gene <= maxGene ; gene++)
	{ 
	  int neq = 0, iFirstIndex ;
	  double index, lastIndex = 0 ;
	  BOOL ok = TRUE ;  /* strict titration */
	  BOOL isLow, lastLow = FALSE ;
	  double sig = 0, idx[tp->rMax] ;
	  GC *gc = arrayp (gx->genes, gene, GC) ;
	  if (gx->targeted && ! gc->targeted)
	    continue ;
	  if (1)
	    if (gx->captured && ! gc->captured)
	      continue ;

	  for (i = 0 ; i < tp->rMax ; i++)
	    {
	      rc = arrayp (gx->runs, tp->rr[i], RC) ; 
	      aa = rc->aa ;
	      if (! aa || ! arrayMax (aa))
		{ ok = FALSE ; break ; }
	      dc = arrayp (rc->aa, gene, DC) ;
	      index = dc->index + rc->fineTune - 1000 ;
	      idx[i] = index ;
	      isLow = dc->isLow ;
	      if (1 && index < minIndex) index = minIndex ;
	      if (1 && index < 2 ) index = 2 ;
	      if (i > 0)
		{
		  float d = (tp->shape)[i-1] ;
		  
		  if (isLow && lastLow)
		    {
		      if (d != 0)
			neq++ ;
		    }
		  else if (! isLow && lastLow)
		    { 
		      if (d <= 0) ok = FALSE ; 
		    }
		  else if (isLow && ! lastLow)
		    {
		      if (d >= 0) ok = FALSE ; 
		    }
		  else if (index > lastIndex)
		    {
		      if (d < 0) ok = FALSE ; 
		      if (index > lastIndex + .1 && d == 0)
			ok = FALSE ;
		      sig += (index - lastIndex) ;
		    }
		  else if (index < lastIndex)
		    {
		      if (d > 0) ok = FALSE ; 
		      if (index < lastIndex - .1 && d == 0)
			ok = FALSE ;
		      sig -= (index - lastIndex) ;
		    }
		  else 
		    {
		      if (d != 0)
			neq++ ;
		    }
		}
	      lastLow = isLow ;
	      lastIndex = index ;
	    }
	 
	  if (neq < 2 &&  ok && sig >  gx->minFoldChange)
	    {
	      float dx = idx[0] - idx[tp->rMax - 1] ;

	      if (dx < 0) dx = -dx ;
	      iFirstIndex = 10 * dx ;

	      array (tp->hh, iFirstIndex, int)++ ;
	      pc = arrayp (tp->pcs, arrayMax(tp->pcs), PC) ;
	      pc->gene = gene ; 
	      pc->sig = sig ;
	      pc->idx = halloc(tp->rMax * sizeof(double), h) ;
	      memcpy (pc->idx, idx, sizeof(idx)) ;
	      nTG++ ;
	    }
	}
      if (nTG > maxTG) maxTG = nTG ;
    }
  /* export the list of titrating genes */
  for (ii = 0, tp = arrp (tt, ii, TT) ; ii < arrayMax (tt) ; tp++, ii++)
    arraySort (tp->pcs, pcsOrder) ;
  if (1)
    {
      int BBMAX=dictMax (gx->runDict) + 100 ;
      ACEOUT ao = aceOutCreate (gx->outFileName, ".TitratingGenes.txt", FALSE, h) ;
      double x ;
      int jj ;
      int nbb, nbbMax = 200 ; /* limit the number of titration quadratic comparisons */

      aceOutDate (ao, "###", gx->title) ;
      /* export the run names in as many lines as needed, one column per titration group */
      aceOutf (ao, "\nCompare") ;
      for (ii = 0, tp = arrp (tt, ii, TT) ; ii < arrayMax (tt) ; tp++, ii++)
	{
	  i = tp->iCompare ;
	  if (i > 0)
	    aceOutf (ao, "\t%s Up", dictName (gx->compareDict, i)) ;
	  else
	    aceOutf (ao, "\t%s Down", dictName (gx->compareDict, - i)) ;
	}
      for (i = 0 ; i < rMax ; i++)
	{
	  aceOutf (ao, "\nRun %d", i + 1) ;
	  for (ii = 0, tp = arrp (tt, ii, TT) ; ii < arrayMax (tt) ; tp++, ii++)
	    aceOutf (ao, "\t%s", i < tp->rMax ? dictName (gx->runDict, tp->rr[i]) : "") ;
	  if (i == 0)
	    {
	      for (nbb = ii = 0 ; ii < arrayMax (tt) ; ii++)
		for (jj = ii+1 ; jj < arrayMax (tt) ; jj++)
		  {
		    if (nbb++ > nbbMax) continue ;
		    if (1 || ii >= arrayMax(tt)/2 || jj < arrayMax(tt)/2)
		      {
			tp = arrp (tt, ii, TT) ;
			if (tp->iCompare > 0)
			  aceOutf (ao, "\t%s Up",  dictName (gx->compareDict, tp->iCompare)) ;
			else
			  aceOutf (ao, "\t%s Down",  dictName (gx->compareDict, - tp->iCompare)) ;
			tp = arrp (tt, jj, TT) ;
			if (tp->iCompare > 0)
			  aceOutf (ao, "/%s Up",  dictName (gx->compareDict, tp->iCompare)) ;
			else
			  aceOutf (ao, "/%s Down",  dictName (gx->compareDict, - tp->iCompare)) ;
		      }
		  }
	    }
	}
      /* position the ppc on the first gene of each serie */
       for (jj = 0, tp = arrp (tt, jj, TT) ; jj < arrayMax (tt) ; tp++, jj++)
	 {
	   /* initialise the bitsets */
	   tp->bbb = halloc (BBMAX * sizeof(BitSet), h) ;
	   for (ii = 0 ; ii < BBMAX && ii < arrayMax (tt) ; ii++)
	     tp->bbb[ii] = bitSetCreate (dictMax (gx->geneDict), h) ;
	 }
      /* export the gene names */
      for (i = 0 ; i < maxTG ; i++)
	{
	  aceOutf (ao, "\n%d", i) ;
	  for (ii = 0, tp = arrp (tt, ii, TT) ; ii < arrayMax (tt) ; tp++, ii++)
	    {
	      if (i < arrayMax (tp->pcs))
		{
		  pc = arrp (tp->pcs, i, PC) ;
		  aceOutf (ao, "\t%s", noX__ (dictName (gx->geneDict, pc->gene), !gx->isINTRON)) ;

		  if (1)
		    {
		      int ix ;

		      for (ix = 0 ; ix < tp->rMax ; ix++)
			{
			  aceOutf(ao, "%c%.2f"
				  , ix ? ',' : '('
				  , pc->idx[ix]
				  ) ;
			}
		      for (ix = 1 ; ix < tp->rMax ; ix++)
			{
			  aceOutf(ao, "%c%.2f"
				  , ix > 1 ? ',' : ';'
				  , pc->idx[ix] - pc->idx[ix-1]
				  ) ;
			}
		      aceOut (ao, ")") ;  
		    }
		  bitSet (tp->bbb[ii], pc->gene) ;
		}
	      else
		aceOut (ao, "\t") ;
	    }

	  for (nbb = 0, ii = 0 ; ii < arrayMax (tt) ; ii++)
	    for (jj = ii+1 ; jj < arrayMax (tt) ; jj++)
	      {
		if (nbb++ > nbbMax) continue ;
		if (1 || ii >= arrayMax(tt)/2 || jj < arrayMax(tt)/2)
		  {
		    if (tp->bbb && arrayExists(tp->bbb[ii]))
		      {
			bb = bitSetCopy (tp->bbb[ii], 0) ;
			x = 100.0 * bitSetAND (bb, tp->bbb[jj])/(i+1) ;
			aceOutf (ao, "\t%.1f", x) ;
			bitSetDestroy (bb) ;
		      }
		    else
		      aceOutf (ao, "\tX") ;
		  }
	      }
	}

      	aceOut (ao, "\n") ;
    }
  
  /* compute the cumulated histo */
  for (ii = 0, tp = arrp (tt, ii, TT) ; ii < arrayMax (tt) ; tp++, ii++)
   {
     for (n = 0, i = arrayMax (tp->hh) - 1 ; i >= 0 ; i--)
       {
	 n += arr (tp->hh, i, int) ;
	 arr (tp->hh, i, int) = n ;
       }
     if (maxIndex <  arrayMax (tp->hh))
       maxIndex = arrayMax (tp->hh) ;
   }

  for (pass = 0 ; pass < 1 ; pass++) /* pass == 0 : strict,  1 : sloped linear fit */
    {
      ACEOUT ao = aceOutCreate (gx->outFileName, ".NumberOfTitratingGeneAboveGivenIndex.txt", FALSE, h) ;
      const char *ccp ;
    
      aceOutDate (ao, "###", gx->title) ;
      aceOutf (ao, "\n## Histogram of the number of titrating genes at least once above a given differential index: %s", gx->title ? gx->title  : "") ;
      aceOut (ao, "\n## The first column gives the minimal difference d of index between the first and last run using a log2 scale") ;
      aceOutf (ao, "\n## The remaining columns give the number of genes for which the index follow the requested shape"
	      , rMax, rMax) ; 

      aceOutf (ao, "\nCompare") ;
      for (ii = 0, tp = arrp (tt, ii, TT) ; ii < arrayMax (tt) ; tp++, ii++)
	{
	  i = tp->iCompare ;
	  if (i > 0)
	    aceOutf (ao, "\t%s Up", dictName (gx->compareDict, i)) ;
	  else
	    aceOutf (ao, "\t%s Down", dictName (gx->compareDict, - i)) ;
	}
      /* export the run names in as many lines as needed, one column per titration group */
      for (i = 0 ; i < rMax ; i++)
	{
	  aceOutf (ao, "\nRun %d", i + 1) ;
	  for (ii = 0, tp = arrp (tt, ii, TT) ; ii < arrayMax (tt) ; tp++, ii++)
	    aceOutf (ao, "\t%s", i < tp->rMax ? dictName (gx->runDict, tp->rr[i]) : "") ;
	}  
      if (gx->hasRunId) 
	{
	  for (i = 0 ; i < rMax ; i++)
	    {
	      aceOutf (ao, "\nRunId %d", i + 1) ;
	      for (ii = 0, tp = arrp (tt, ii, TT) ; ii < arrayMax (tt) ; tp++, ii++)
		{
		  ccp = 0 ;
		  if (i < tp->rMax)
		    {
		      rc = arrayp (gx->runs, tp->rr[i], RC) ; 
                      ccp = rc->runid ? stackText (gx->info, rc->runid) : 0 ;
		    }
		  aceOutf (ao, "\t%s", ccp ? ccp : "") ;
		}  
	    }
	}
      if (gx->hasRunSample) 
	{
	  for (i = 0 ; i < rMax ; i++)
	    {
	      aceOutf (ao, "\nRunId %d", i + 1) ;
	      for (ii = 0, tp = arrp (tt, ii, TT) ; ii < arrayMax (tt) ; tp++, ii++)
		{
		  ccp = 0 ;
		  if (i < tp->rMax)
		    {
		      rc = arrayp (gx->runs, tp->rr[i], RC) ; 
                      ccp = rc->sample ? stackText (gx->info, rc->sample) : 0 ;
		    }
		  aceOutf (ao, "\t%s", ccp ? ccp : "") ;
		}  
	    }
	}
      if (gx->hasRunTitle) 
	{
	  for (i = 0 ; i < rMax ; i++)
	    {
	      aceOutf (ao, "\nRunId %d", i + 1) ;
	      for (ii = 0, tp = arrp (tt, ii, TT) ; ii < arrayMax (tt) ; tp++, ii++)
		{
		  ccp = 0 ;
		  if (i < tp->rMax)
		    {
		      rc = arrayp (gx->runs, tp->rr[i], RC) ; 
                      ccp = rc->title ? stackText (gx->info, rc->title) : 0 ;
		    }
		  aceOutf (ao, "\t%s", ccp ? ccp : "") ;
		}  
	    }
	}
 
      for (i = 0 ; i < maxIndex ; i++)
	{
	  aceOutf (ao, "\n%.1f", i/10.0) ;
	  for (ii = 0, tp = arrp (tt, ii, TT) ; ii < arrayMax (tt) ; tp++, ii++)
	    aceOutf (ao, "\t%d", array (tp->hh, i, int)) ;
	}
      aceOutf (ao, "\n") ;
    }
  ac_free (h) ;
  return ;
} /* gxTitration */
#endif
/*************************************************************************************/
/*************************************************************************************/

static int pcOrder (const void *a, const void *b)
{
  const PC *up = (const PC *)a, *vp = (const PC *)b ;
  double n ;

  n = up->dIndex - vp->dIndex ; 
    if (n > 0) return -1 ;
    if (n < 0) return 1 ;

  return 0 ;
}

/*************************************************************************************/

static int pcDvarOrder (const void *a, const void *b)
{
  const PC *up = (const PC *)a, *vp = (const PC *)b ;
  double n ;

  n = up->dvar - vp->dvar ; 
  if (n > 0) return -1 ;
  if (n < 0) return 1 ;

  n = up->index1 + up->index2 - vp->index1 - vp->index2 ;
  if (n > 0) return -1 ;
  if (n < 0) return 1 ;

  n = up->dIndex - vp->dIndex ; 
  if (n > 0) return -1 ;
  if (n < 0) return 1 ;

  return 0 ;
} /* pcDvarOrder */

/*************************************************************************************/
/* sort the most significant genes */
static int pcSigOrder (const void *a, const void *b)
{
  const PC *up = (const PC *)a, *vp = (const PC *)b ;
  double n, nid ;

  n = up->sig - vp->sig ; 
  if (n > 0) return -1 ;
  if (n < 0) return 1 ;

  nid = up->geneid - vp->geneid ;
  if (nid > 0) nid = 1 ;
  if (nid < 0) nid = -1 ;
  if (nid && up->geneid && vp->geneid) nid = 0 ;

  n = 3 * nid + up->index1 + up->index2 - vp->index1 - vp->index2 + (up->dIndex - vp->dIndex)/10.0 ;
  if (n > 0) return -1 ;
  if (n < 0) return 1 ;


  n = up->dIndex - vp->dIndex ; 
  if (n > 0) return -1 ;
  if (n < 0) return 1 ;

  return 0 ;
} /* pcSigOrder */

/*************************************************************************************/
/* Evaluate the score of the double histogram */
/* copied from the awk script: histoselect.awk */
static double  doScoreDoubleHisto (double sx0, double sy0, double dx, double dy, double ratio_bound)
{
  double dzu, z = ratio_bound ; 

  if (sx0 < 3 * sy0) dx -= sx0 ;
  if (1)
    {  /* methode additive, count area only if above 2 fold. linear interpolation in [1.5, 2]  */
      z = 1 ;
      if (dy > 0)
	{ 
	  double m = ratio_bound ;
	  z = dx/dy ;
	  if (z < m)
	    z = 2.2 * z - (2.2 * m - 1) ;
	  if (z < 0)
	    z = 0 ;
	  if (z > 1)
	    z = 1 ; 
	}
    }
  else
    {    /* methode multiplicative */
      if (dy > 0) 
	z = dx/dy ;
      if (z < 3)
	z = 2 * z - 3 ; 
      if (z > 0) 
	{
	  if (z > ratio_bound) 
	z = ratio_bound ;
	} 
    }
  if (z < 0) z = 0 ;

  dzu = (dx - dy) * z ; 

  return dzu > 0 ? dzu : 0 ;
}

static double  scoreDoubleHisto (double sx0, double sy0, double dx, double dy, float *ratio_boundp)
{
  double R, z, bestz = 0 ;
  *ratio_boundp = 1 ;
  for (R = 1 ; R <= 5 ; R++)
    {
      z = doScoreDoubleHisto (sx0, sy0, dx, dy, R) ;
      if (z > bestz)
	{ bestz = z ; *ratio_boundp = R ; }
    }
  return bestz ;
}

/*************************************************************************************/
/*************************************************************************************/

static COMPARE *gxPairCmp (GX *gx, int run1, int run2)
{
  char buf[128] ;
  int k ;
  COMPARE *cmp ;

  if (run1 > run2)
    { int run0 = run1 ; run1 = run2 ; run2 = run0 ;}

  sprintf (buf, "virtual_%d_%d", run1, run2) ;
  dictAdd (gx->pairDict, buf, &k) ;
  cmp = arrayp (gx->pairs2cmp, k, COMPARE) ;
  cmp->compare = k ;
  cmp->vRun1 = run1 ;
  cmp->vRun2 = run2 ;
  cmp->private = TRUE ;

  return cmp ;
} /* gxPairCmp */

/*************************************************************************************/

static void gxPairScoreRegister (GX *gx, int iCompare, int run1, int run2, int gene, float score, float pValue)
{
  int rrg,  sign = 1 ;
  char buf[128] ;
  int k ;
  static int kk = 0 ;
  if (score > 200) score = 200 ;
  if (score < -200) score = -200 ;
  if (run1 > run2)
    { 
      int run0 = run1 ; run1 = run2 ; run2 = run0 ; sign = -1 ;
    }

  if (gx->pairScores2compare)
    {
      sprintf (buf, "%d_%d_000", run1, run2) ;
      dictAdd (gx->pairScoresDict, buf, &rrg) ;
      keySet (gx->pairScores2compare, rrg)  = iCompare ;

      sprintf (buf, "%d_%d_%d", run1, run2, gene) ;
      k = dictAdd (gx->pairScoresDict, buf, &rrg) ;
      if (gene && ! k)
	{
	  if (kk++ < 20)
	    {
	      fprintf (stderr, "New score for run1=%d run2=%d gene = %d olds=%g news=%g  %s %s %s\n"
		       , run1, run2, gene
		       , array (gx->pairScores, rrg, float) - 1000
		       , sign * score
		       , dictName (gx->runDict, run1)
		       , dictName (gx->runDict, run2)
		       , dictName (gx->geneDict, gene)
		       ) ;
	      invokeDebugger () ;
	    }
	}
      array (gx->pairScores, rrg, float) = sign * score + 1000 ;
      array (gx->pValues, rrg, float) = pValue ;
    }
} /* gxPairScoreRegister */

/**************************************************************/
/**************************************************************/
/*************************************************************************************/

static float gxScoreHistos (GX *gx, int gene, RC *rc1, RC *rc2, int run1, int run2, double *b, double *c, float ch2
			    , BOOL isVirtual, double threshold, int *sensep)
{
  float ratio_bound  = gx->ratio_bound ;
  int jjShift = gx->histo_shifting ;  
  GC *gc = arrp (gx->genes, gene, GC) ;  
  float boni[2], scores[2], score, scoreA[2], scoreB[2] ;
  int ratioA[2], ratioB[2] ;
  int iiix1[2], iiix2[2], iiiy1[2], iiiy2[2], sens[2] ; 
  int aa1[2], aa2[2], bb1[2], bb2[2] ;
  int i, iMin ;
  int debug = 0 ;
  int raBest1 = 0, raBest2 = 0, rbBest1 = 0, rbBest2 = 0 ; 
  BOOL killScore = FALSE ;
  double  aBest1 = 0, aBest2 = 0, bBest1 = 0, bBest2 = 0, av1 = 0, av2 = 0, av1s[2], av2s[2], 
    sx0 = 0, sy0 = 0, sx, sy, sx1, sy1, mx, my, dza, dzb, 
    x, y, dx, dy, ixa, iya, ixn, iyn, jxa, jya, jxn, jyn ;
  double peakx, peaky, dp, bonus ;
  int a1 = 0, a2 = 0, b1 = 0, b2 = 0 ;
  int ix1 = 0, ix2 = 0, jx1 = 0, jx2 = 0 ;
  int iy1 = 0, iy2 = 0, jy1 = 0, jy2 = 0 ;
  int jjj, ijjx, ijjy, ijjxa, ijjxb, ijjya, ijjyb ;
  /* normalize = FALSE ;   done inside gxMakeComparativeHisto normalize = gx->normalize */

  if (debug > 0 
      && strstr (dictName(gx->runDict, run1), "l")
      && strstr (dictName(gx->geneDict, gene), "zzzzAldh3a1")
      )
    { debug = 1 ; invokeDebugger () ; }
  if (! strncmp (dictName(gx->geneDict, gene), "zzzzCyp1a",5))
    debug = 1 ;
  if (threshold < rc1->NA_crossover_index)
    threshold = rc1->NA_crossover_index ;
  if (threshold < rc2->NA_crossover_index)
    threshold = rc2->NA_crossover_index ;
  threshold -= 0 ;  /* 2015_01_28, was zero 2014_03_01 was 1.5 ; 2013_08_08 was 2 dromadaire */
  if (0 && gx->isMA)
    threshold = 8 ;
  if (gx->threshold && threshold < gx->threshold)
    threshold = gx->threshold ;
  if (gx->isSNP || gc->isSNP)
    threshold = 0 ;
 
  /* initialize the loop */
  iMin = 10 * threshold ;                        /* v1 j'utilise cette condition */
  iMin = (gx->digital ? 0 : 10 * threshold)  ;   /* v2, je retablit cette condition */
  iMin = 0 ; /* 2015_06_6 why did we play with that ? */

  /* compute the total surface and the max of the continuous zone (above index 2) */
  for (sx = 0, i = 0 ; i < 256 ; i++)
    sx += b[i] ;
  for (sy = 0, i = 0 ; i < 256 ; i++)
    sy += c[i] ;
  if (! gx->isSNP)
    { sx0 = b[0] ; sy0 = c[0] ; }

  for (mx = 0, i = 00 ; i < 256 ; i++)
    if (b[i] > mx) mx = b[i] ;
  for (my = 0, i = 00 ; i < 256 ; i++)
    if (c[i] > my) my = c[i] ;
  
  for (jjj = 0 ; jjj < 2 ; jjj++)
    {
      sens[jjj] = 0 ;
      scores[jjj] = 0 ;
      scoreA[jjj] = 0 ;
      scoreB[jjj] = 0 ;
      ratioA[jjj] = ratioB[jjj] = 0 ;
      boni[jjj] = 0 ;
      av1s[jjj] = 0 ;
      av2s[jjj] = 0 ;
      iiix1[jjj] = 0 ; iiix2[jjj] = 0 ; 
      iiiy1[jjj] = 0 ; iiiy2[jjj] = 0 ; 
      aa1[jjj] = 0 ; aa2[jjj] = 0 ; 
      bb1[jjj] = 0 ; bb2[jjj] = 0 ; 
    }

  for (jjj = 0 ; jjj < 2 ; jjj++)
    {
      peakx = 0; peaky = 0; dp = 0; bonus = 0 ;

      aBest1 = 0 ; aBest2 = 0 ; bBest1 = 0 ; bBest2 = 0 ;
      raBest1 = 0 ; raBest2 = 0 ; rbBest1 = 0 ; rbBest2 = 0 ;
      a1 = 0 ; a2 = 0 ; b1 = 0 ; b2 = 0 ;
      dx = 0 ; dy = 0 ; 
      sx1 = sy1 = 0 ;
      ijjxa = 0 ; ijjxb = 0 ; ijjya = 0 ; ijjyb = 0 ;
      ixa = 0 ; ixn = 1 ; iya = 0 ; iyn = 1 ;
      jxa = 0 ; jxn = 1 ; jya = 0 ; jyn = 1 ;
      ix1 = 0 ; ix2 = 0 ; iy1 = 0 ; iy2 = 0 ;
      jx1 = 0 ; jx2 = 0 ; jy1 = 0 ; jy2 = 0 ;
      
      
      /* look for the best zone from 0 to i 
       * collect at the same time the position of the dominant peaks
       */
      /* true average */
      for (i = 0 ; i < 256 ; i++)
	{
	   x = b[i] ; y = c[i] ; 
	   if (x > 0 && x > mx/300 ) 	{ jxa += i * x ; jxn += x ; }
	   if (y > 0 && y > my/300 ) 	{ jya += i * y ; jyn += y ; }
	}
	  
      for (i = 0 - jjShift ; i < 256 + jjShift ; i++)
	{ 
	  if (0)
	    { x = b[i] ; y = b[i] ; }
	  else
	    {
	      ijjx = i + jjShift * (1 - jjj) ;
	      ijjy = i + jjShift * jjj ;
	      
	      if (ijjx < 0 || ijjx > 255) x = 0 ;
	      else x = b[ijjx] ;

	      if (ijjy < 0 || ijjy > 255) y = 0 ;
	      else y = c[ijjy] ;
	    }
	  
	  dx += x ; dy += y ;
	  
	  if (i < iMin) continue ;
	  if (i == iMin) {  x = dx ; y = dy ; }
	  
	  /* average of significant region */
	  if (x > 1.5 * y && x > mx/3 ) { ixa += ijjx * dx ; ixn += dx ; sx1 += x ; if(! ix1) ix1 = ijjx ; ix2 = ijjx ; }
	  if (y > 1.5 * x && y > my/3) 	{ iya += ijjy * dy ; iyn += dy ; sy1 += y ; if(! iy1) iy1 = ijjy ; iy2 = ijjy ; }

	  dza = scoreDoubleHisto (sx0, sy0, dx,dy, &ratio_bound) ;
	  if (dza > aBest1) { aBest1 = dza ; a1 = ijjx ; raBest1 = 100.0 * dx / (dx + dy + .00001) ; }
	  
	  dzb = scoreDoubleHisto (sy0, sx0, dy,dx, &ratio_bound) ;
	  if (dzb > bBest1) { bBest1 = dzb ; b1 = ijjy ; rbBest1 = 100.0 * dy / (dx + dy + .00001) ; }
	  if (x + y > 0 &&
	      ( debug > 1 ||
		( strstr (dictName(gx->geneDict, gene), "ItgamZZ") > 0)
		) 
	      )
	    fprintf (stderr, "+++ %s i=%d  x=%.2f y=%.2f dx=%.1f dy=%.1f dza=%.1f dzb=%.1f aB1=%.1f aB2=%.1f bB1=%.1f bB2=%.1f ijjx=%d x=%.1f jxa=%.1f jxn=%.1f av=.%1f\n", dictName(gx->geneDict, gene), i, x,y,dx,dy,dza,dzb,aBest1, aBest2, bBest1,bBest2,ijjx,x,jxa,jxn,jxa/(-.01+jxn)) ;
	}
      

      /* look for the best zone from na  to i */
      dx = 0 ; dy = 0 ;
      for (i = 255 + jjShift ; i >= 0 - jjShift ; i = i - 1)
	{ 
	  if (0)
	    { x = b[i] ; y = b[i] ; }
	  else
	    {
	      ijjx = i + jjShift * (1 - jjj) ;
	      ijjy = i + jjShift * jjj ;
	      
	      if (ijjx < 0 || ijjx > 255) x = 0 ;
	      else x = b[ijjx] ;

	      if (ijjy < 0 || ijjy > 255) y = 0 ;
	      else y = c[ijjy] ;
	    }

	  dx += x ; dy += y ;
	  if (i < iMin) continue ;
	  
	  if (x > 1.5 * y && x > mx/3 ) { if(! jx2) jx2 = ijjx ; jx1 = ijjx ; }
	  if (y > 1.5 * x && y > my/3) 	{ if(! jy2) jy2 = ijjy ; jy1 = ijjy ; }
	  
	  dza = scoreDoubleHisto (0,0, dx,dy, &ratio_bound) ;
	  if (dza > aBest2) { aBest2 = dza ; a2 = ijjx ; raBest2 = 100.0 * dx / (dx + dy + .00001) ; }
	  
	  dzb = scoreDoubleHisto (0, 0, dy,dx, &ratio_bound) ;
	  if (dzb > bBest2) { bBest2 = dzb ; b2 = ijjy ; rbBest2 = 100.0 * dy / (dx + dy + .00001) ; }

	  if ((debug > 1 && i%5 == 0) ||
	      ( strstr (dictName(gx->geneDict, gene), "ItgamZZ") > 0)
	      ) 
	    fprintf (stderr, "--- %s i=%d  x=%.2f y=%.2f dx=%.1f dy=%.1f dza=%.1f dzb=%.1f aBest1=%.1f aBest2=%.1f bBest1=%.1f bBest2=%.1f\n", dictName(gx->geneDict, gene), i, x,y,dx,dy,dza,dzb,aBest1, aBest2, bBest1,bBest2) ;
	} 
      
      dx = 0 ; dy = 0 ;
      if (aBest1 > aBest2)  { dx = ix2 - ix1 ;} /* D8  measure the zone width and add a bonus */
      if (aBest1 < aBest2)  { dx = jx2 - jx1 ;}
      if (bBest1 > bBest2)  { dy = iy2 - iy1 ;}
      if (bBest1 < bBest2)  { dy = jy2 - jy1 ;}
      
       score = 0 ;
      if (aBest1 >= aBest2 && bBest2 >= bBest1) 
	{
	  sens[jjj] = 1 ;
	  scores[jjj] = aBest1 + bBest2 + bonus ;
	  scoreA[jjj] = aBest1 ; 	  
	  scoreB[jjj] = bBest2 ; 
	  ratioA[jjj] = raBest1 ; 	  
	  ratioB[jjj] = rbBest2 ; 
	  iiix1[jjj] = ix1 ; iiix2[jjj] = ix2 ; 
	  iiiy1[jjj] = jy1 ; iiiy2[jjj] = jy2 ; 
	  aa1[jjj] = 0 ; bb1[jjj] = b2 > 1 ? b2 : 1 ;
	  aa2[jjj] = a1 ; bb2[jjj] = 256 ;
	}
      if (aBest1 <= aBest2 && bBest2 <= bBest1)
	{
	  sens[jjj] = 2 ;
	  scores[jjj] = aBest2 + bBest1 + bonus ;
	  scoreA[jjj] = aBest2 ; 	  
	  scoreB[jjj] = bBest1 ; 
	  ratioA[jjj] = raBest2 ; 	  
	  ratioB[jjj] = rbBest1 ; 
	  iiix1[jjj] = jx1 ; iiix2[jjj] = jx2 ; 
	  iiiy1[jjj] = iy1 ; iiiy2[jjj] = iy2 ; 
	  aa1[jjj] = a2 > 1 ? a2 : 1 ; bb1[jjj] = 0 ;
	  aa2[jjj] = 256 ; bb2[jjj] = b1 ;
	}

      if (0)
	{
	  bonus = 0 ;

	  if (1)
	    {
	      peakx = ixa/ixn ;  peaky = iya/iyn ; dp = peaky - peakx ; 
	      if (dp > 0) 
		bonus = peaky ;
	      else
		bonus = peakx ;
	    }
	  if (0)  /*  danielle */
	    {
	      if (100 * scoreA[jjj] > 15 * sx && 100 * scoreB[jjj] > 10 * sy)
		bonus += .5 * (dp - 20)/10 * (sx + sy)/85 ;
	      if (100 * scoreB[jjj] > 15 * sy && 100 * scoreA[jjj] > 10 * sx)
		bonus += .5 * ((dp - 20)/10) * (sx + sy)/85 ;
	    }

	  if (0 && 10 * sx1 > sx && 10 * sy1 > sy)
	    {
	      /*  peak distance in significant differential positions */
	      peakx = ixa/ixn ;  peaky = iya/iyn ; dp = peaky - peakx ; if (dp < 0) dp = -dp ;
	      if (dp > 40) dp = 40 ;
	      bonus += 1 * 0.5 * 0.5 * dp * (sx + sy)/100 ;
	      
	      /* peak distance in significant not necesseraly differential positions */
	      peakx = jxa/jxn ;  peaky = jya/jyn ; dp = peaky - peakx ; if (dp < 0) dp = -dp ;
	      if (dp > 40) dp = 40 ;
	      dp -= 10 ;  
              bonus +=  7 * (dp/40)  * (sx + sy)/100 ; 
	      /* fprintf(stderr, "# %s Bonus= %.1f delta=%.1f = %.1f | %.1f\n", dictName(gx->geneDict, gene),  bonus, dp, peakx, peaky) ; */
	    }

	  if (bonus < 0) bonus = 0 ;
	  if (0 && bonus >  0) fprintf (stderr, "##dp\t%s\t%.1f\n", dictName(gx->geneDict, gene), dp) ;
	}

      av1s[jjj] = jxn ? jxa/jxn : 0 ;  
      av2s[jjj] = jyn ? jya/jyn : 0 ;  
      boni[jjj] = bonus ;
      scores[jjj] *= (100.0 + bonus)/100.0 ;

      if (debug > 0)
        printf ("==+ %s\tscore=%.1f\tsens %d\t\tA %.1f::%d:%d//%.1f::%d:%d//%.1f::%d:%d\tB %.1f::%d:%d//%.1f::%d:%d//%.1f::%d:%d\n"
		, dictName(gx->geneDict, gene)
		,scores[jjj], jjj, scoreA[jjj], ijjxa,ijjxb,aBest1,ix1,ix2,aBest2,jx1,jx2,scoreB[jjj],ijjya,ijjyb,bBest1,iy1,iy2,bBest2,jy1,jy2) ;

      if (debug > 0 ||
	  ( strstr (dictName(gx->geneDict, gene), "zzWFDC2") > 0 && strstr (dictName(gx->runDict, run1),"HR") > 0)
	  )
	fprintf (stderr, "==+ %s\tscore=%.1f\tsens %d\t\tA %.1f::%d:%d//%.1f::%d:%d//%.1f::%d:%d\tB %.1f::%d:%d//%.1f::%d:%d//%.1f::%d:%d\n"
		 , dictName(gx->geneDict, gene)
		 ,scores[jjj], jjj, scoreA[jjj], ijjxa,ijjxb,aBest1,ix1,ix2,aBest2,jx1,jx2,scoreB[jjj],ijjya,ijjyb,bBest1,iy1,iy2,bBest2,jy1,jy2) ;
    }    
  if (scores[0] < scores[1]) /* take the least favorable score */
    jjj = 0 ;
  else
    jjj = 1 ;
  /* TOTOTOTO */
  killScore = FALSE ;
  if ( 
      (sens[0] == 1 && av2s[jjj] < 10 * threshold) || /* a is left of b */
      (sens[0] == 2 && av1s[jjj] < 10 * threshold)    /* b is left of a */
      )
   killScore = TRUE ;  /* mieg 2023_09_25 kill the score if the high peak is below the threshold */
  /* killScore = FALSE ; */
  if (! killScore && sens[0] == sens[1])
    {
      score = scores[jjj] ;
      bonus = boni[jjj] ;
      aBest1 = scoreA[jjj] ;
      bBest1 = scoreB[jjj] ;
      raBest1 = ratioA[jjj] ;
      rbBest1 = ratioB[jjj] ;
      av1 = av1s[jjj] ;
      av2 = av2s[jjj] ;
      ix1 = iiix1[jjj] ; 
      ix2 = iiix2[jjj] ;   
      iy1 = iiiy1[jjj] ;      
      iy2 = iiiy2[jjj] ;
      a1 = aa1[jjj] ; a2 = aa2[jjj] ; b1 = bb1[jjj] ; b2 = bb2[jjj] ;
      /*
	bonus = 0 ;
	if (ix2 - ix1 > 20 && aBest1 > 6 * (sx + sy)/80)
	bonus += aBest1/2 ;
	bonus = (ix2 - ix1) + (iy2 - iy1)
      */

      /* reextend the grey zone : start from the current limits 
       * grow as long as R ratio is not reached and diminish the score
       */
      if (1) /* zoneGrise */
	{
	  double bonusA = 0, bonusB = 0 ;
	  if (sens[0] == 1) /* a is left of b */
	    {
	      for (bonus = 0, ijjx = a2 ; ijjx > a1 ; ijjx--)
		{
		  x =  b[ijjx] ; y = c[ijjx] ;
		  if (x > y * ratio_bound)
		    {
		      if (c[ijjx + 1] == 0) ijjx++ ;
		      if (c[ijjx + 1] == 0) ijjx++ ;
		      break ;
		    }
		  if (x > y) bonus += x - y ;
		}
	      bonusA = bonus ; bonus = 0 ;
	      a2 = ijjx - jjShift ;
	      if (a1 == 0 && c[0] * ratio_bound > b[0]) a1 = 1 ;
	      if(a2 <= 0) a2 = 1 ;
	      for (ijjx = b1 ; ijjx < b2 ; ijjx++)
		{
		  x =  b[ijjx] ; y = c[ijjx] ;
		  if (ijjx > 10 * threshold && y > x * ratio_bound)
		    {
		      if (ijjx - 1 > 0 && b[ijjx - 1] == 0) ijjx-- ;
		      if (ijjx - 1 > 0 && b[ijjx - 1] == 0) ijjx-- ;
		      break ;
		    }
		  if (y > x) bonus += y - x ;
		}
	      b1 = ijjx ;
	      aBest1 -= bonusA ; bBest1 -= bonus ;
	      score -= (bonusA + bonus) ;
	      
	      if (score < 0) score = 0 ;
	    }
	  
	  else  /* b is left of a */
	    {
	      for (bonus = 0, ijjx = b2 ; ijjx > b1 ; ijjx--)
		{
		  x =  b[ijjx] ; y = c[ijjx] ;
		  if (y > x * ratio_bound)
		    {
		      if (b[ijjx + 1] == 0) ijjx++ ;
		      if (b[ijjx + 1] == 0) ijjx++ ;
		      break ;
		    }
		  if (y > x) bonus += y - x ;
		}
	      b2 = ijjx ;
	      bonusB = bonus ; bonus = 0 ;
	      if(b2 <= 0) b2 = 1 ;
	      if (b1 == 0 && b[0] * ratio_bound > c[0]) b1 = 1 ;
	      for (ijjx = a1 ; ijjx < a2 ; ijjx++)
		{
		  x =  b[ijjx] ; y = c[ijjx] ;
		  if (ijjx > 10 * threshold && x > y * ratio_bound)
		    {
		      if (ijjx - 1 > 0 && c[ijjx - 1] == 0) ijjx-- ;
		      if (ijjx - 1 > 0 && c[ijjx - 1] == 0) ijjx-- ;
		      break ;
		    }
		  if (x > y) bonus += x - y ;
		}
	      a1 = ijjx ;
	      aBest1 -= bonus ; bBest1 -= bonusB ;
	      score -= (bonusB + bonus) ;

	      if (score < 0) score = 0 ;
	    }
	}

    }
  else
    score = 0 ;

  if (score <= 0 || aBest1 <= 0 || bBest1 <= 0)   
    {
      jjj = -1 ; 
      score = 0 ;
      bonus = 0 ;
      aBest1 = 0 ;
      bBest1 = 0 ;
      raBest1 = 0 ;
      rbBest1 = 0 ;
      av1 = av2 = 0 ;
      ix1 = 0 ;
      ix2 = 0 ;
      iy1 = 0 ;
      iy2 = 0 ;
      a1 = a2 = b1 = b2 = 0 ;
    }
  if (debug > 0 || 
      ( strstr (dictName(gx->geneDict, gene), "ItgamZZ"))
      ) 
    {
      fprintf (stderr, "=== %s\tscore=%.1f\tsens %d\t\tA %.1f::%d:%d ra=%d\tB %.1f::%d:%d  rb=%d\n"
	       , dictName(gx->geneDict, gene)
	       ,score, jjj, aBest1, ix1, ix2, raBest1, bBest1, iy1, iy2, rbBest1
	       ) ;
      for (i = 0 ; i < 256 ; i++)
	fprintf (stderr, "...%s\t%d\t%.1f\t%.1f\n"
		 ,  dictName(gx->geneDict, gene)
		 , i, b[i], c[i]
		 ) ;
    }      

  ch2 /= 2 ;  /* 2014_06_07 */
  if (ch2 > 30) ch2 = 30 ;
  if (ch2 < 0) ch2 = 0 ;
  ch2 = 0 ; /* dromadaire, this boost may be a false good idea 2014_06_08 */

  if (killScore)
    {
      gc->score = 0 ;
      gc->score1 = 0 ;
      gc->score2 = 0 ;
    }
  else
    {
      gc->score = (score + ch2  < 200  ? score + ch2 : 200) ;
      gc->score1 = (aBest1 + ch2/2 < 100 ? aBest1 + ch2/2 : 100) ;
      gc->score2 = (bBest1 + ch2/2 < 100 ? bBest1 + ch2/2 : 100) ;
    }
  gc->ratioa = raBest1 ;
  gc->ratiob = rbBest1 ;
  gc->bonus = bonus ;
  gc->ix1 = ix1 ; gc->ix2 = ix2 ;
  gc->iy1 = iy1 ; gc->iy2 = iy2 ;
  
  gc->av1 = av1 ; gc->av2 = av2 ;  gc->isLow1 = 0 ;  gc->isLow2 = 0 ;
  gc->za1 = a1 ; gc->za2 = a2 ;  
  gc->zb1 = b1 ; gc->zb2 = b2 ;

  if (rc1 && rc1->aa && rc2 && rc2->aa && arrayMax(rc1->aa) > gene && arrayMax(rc2->aa) > gene)
    {
      DC *dc ;
      dc = arrp (rc1->aa, gene, DC) ;
      gc->av1 = dc->index ; gc->isLow1 = dc->isLow ? (dc->tags == 0 ? 2 : 1) : 0 ;
      dc = arrp (rc2->aa, gene, DC) ;
      gc->av2 = dc->index ; gc->isLow2 = dc->isLow ? (dc->tags == 0 ? 2 : 1) : 0 ;
    }

  if (sensep) *sensep = sens[0] ;

  return score > 0 ? score : 0 ;
}  /* gxScoreHistos */

static void gxScoreByOverlap (void)
{
  messcrash ("gxScorebyOverlap is obsolete code") ;
#ifdef JUNK
	  int i, nx, old, x0, x1 ;
		  double z, s0, s1, s0a, s1a, sum0, sum1, *a0p, *a1p, mx0, mx1, sdiff ;
		  double hh0s[256], hh1s[256] ;
		  /*  int danRegion0 = 20, danRegion1 = 70, danRegion2 = 80, danRegion3 = 120 ; */

		  double dan00, dan01, dan10, dan11, danDiff, danStep, danScore ;
  int danRegion0 = 10, danRegion1 = 35, danRegion2 = 40, danRegion3 = 60 ; /* since we smooth evey .2, not .1 */
  double z, s0, s1, s0a, s1a, sum0, sum1, *a0p, *a1p, mx0, mx1, sdiff ;
 /* obsolete code */

  /* compute the 2 areas */
  /* find the extent starting at index 5 over 1/10 of max */
  for (s0a = s1a = sum0 = sum1 = sdiff = 0, x0 = x1 = i = nx = old = 0, a0p = hh0s, a1p = hh1s ; i < 256 ; a0p++, a1p++, i++)
    {
      sum0 += *a0p ; sum1 += *a1p ;
      if (*a0p + *a1p > mx0 + mx1) 
	{
	  if (x0 == 0) x0 = i ;
	  x1 = i ;
	}
      if (*a0p < *a1p) 
	{
	  s0a += *a1p - *a0p ;
	  sdiff +=  *a1p - *a0p ;
	  if (2 * (*a0p) < *a1p && *a1p > mx0)
	    {
	      s0a += *a1p - *a0p ; 
	      if (old == -1) nx++ ; old =  1 ; 
	    }
	  if (5 * (*a0p) < *a1p && *a1p > mx0)
	    s0a += *a1p - *a0p ; 
	  
	  /* danielle's loaded oreilles */
	  if (i <= danRegion0)                         /* below index 1, cumulate the values */
	    { dan10 += *a0p; dan11 += *a1p ; }
	  else if (i == danRegion1)
	    {
	      if (dan11 > 10 * dan10 && dan11 > mx1) 
		danDiff += 3  * (dan11 - dan10) * (1.0 + dan11 / (dan11 + dan10)) ; /* was 1.5 in D3 */
	    }
	  else
	    {
	      if (*a1p > 0) 
		danDiff += (*a1p - *a0p) * (1.0 + (*a1p)/(*a0p + *a1p)) ;
	      if ((*a1p > 3 * (*a0p))  && 
		  (*a1p >  mx1/10) &&
		  (*a0p < mx0/1)
		  )
		{
		  if (i > danRegion3)
		    danStep += 2.0 * (1.0 + (*a1p)/(*a0p + *a1p)) ;
		  else if (i > danRegion2)
		    danStep += 2.0 * (1.0 + (*a1p)/(*a0p + *a1p)) ;
		  else if (i > danRegion1)
		    danStep += 1.5 * (1.0 + (*a1p)/(*a0p + *a1p)) ;
		}
	    }
	}
      else
	{
	  s1a += *a0p - *a1p ;  
	  sdiff += *a0p - *a1p ;
	  if (2 * (*a1p) < *a0p && *a0p > mx1)
	    {
	      s1a += *a0p - *a1p ; 
	      if (old == 1) nx++ ; old =  -1 ; 
	    }
	  if (5 * (*a1p) < *a0p && *a0p > mx1)
	    s1a += *a0p - *a1p ;
	  
	  /* danielle's loaded oreilles */
	  if (i <= danRegion0)                         /* below index 1, cumulate the values */
	    { dan00 += *a0p; dan01 += *a1p ; }
	  else if (i == danRegion1)
	    {
	      if (dan00 > 10 * dan01 && dan00 > mx0) 
		danDiff += 3.0 * (dan00 - dan01) * (1.0 + dan00 / (dan00 + dan01)) ;
	    }
	  else
	    {
	      if (*a0p > 0) 
		danDiff += (*a0p - *a1p) * (1.0 + (*a0p)/(*a0p + *a1p)) ;
	      if ((*a0p > 3 * (*a1p))  && 
		  (*a0p >  mx0/10) &&
		  (*a1p < mx1/1)
		  )
		{
		  if (i > danRegion3)
		    danStep += 2.0 * (1.0 + (*a0p)/(*a0p + *a1p)) ;
		  else if (i > danRegion2)
		    danStep += 2.0 * (1.0 + (*a0p)/(*a0p + *a1p)) ;
		  else if (i > danRegion1)
		    danStep += 1.5 * (1.0 + (*a0p)/(*a0p + *a1p)) ;
		}
	    }
	}
      
      if (0 &&
	  ! strcmp (dictName(gx->geneDict, pc->gene), "EPHA5") &&
	  ! strcmp (dictName(gx->runDict, run1), "Favorable_NB_178")
	  )
	{
	  if (0) fprintf (stderr, "# EPHA5 i=%d x=%.1f y=%.1f danDiff=%.1f dnaStep=%.1f mx0=%.1f, mx1=%.1f\n"
			  , i, *a0p, *a1p
			  , danDiff, danStep, mx0,mx1) ;
	}
    }
  if (0)
    {  /* danielle's method step + oreilles - cross-penalty */
      if (x0 < 50) x0 = 50 ;
      if (x1 < x0) x1 = x0 ;
      pc->overlap = x1 - x0 + sdiff * 230/(1 + sum0 + sum1) - 30 * nx ;
    }

    {  /* danielle's method step + oreilles - cross-penalty */
      danDiff *= 1 ;
      danStep *= 2 ; /* method D3. was *= 1 in D2, may be we need to double again since we smooth every .2 */
      pc->danielle = danDiff + danStep - 100 * nx ;
    }

#endif
  return ;
} /* gxScoreByOverlap */

/*************************************************************************************/
#ifdef JUNK
static float gxScoreOneGeneHisto (GX *gx, int gene, int run1, int run2, float *pValue)
{
  AC_HANDLE h = ac_new_handle () ;
  int sense, nRun1 = 0, nRun2 = 0 ;
  RC *rc1, *rc2 ;
  double threshold = 0 ;
  BOOL isVirtual = FALSE ;
  float ch2 = 0, score = 0 ;
  double hh0s[256], hh1s[256] ;
  KEYSET vA = 0, vB = 0, vA0, vB0 ;
  
  vA0 = keySetHandleCreate (h) ;
  vB0 = keySetHandleCreate (h) ;
  
  rc1 = arrp (gx->runs, run1, RC) ;
  rc2 = arrp (gx->runs, run2, RC) ; 
  
  keySet (vA0, 0) = run1 ;
  keySet (vB0, 0) = run2 ;
  
  vA =  rc1->runs ? rc1->runs : vA0 ;
  vB =  rc2->runs ? rc2->runs : vB0 ;
  
  if (1)
    {
      memset (hh0s, 0, sizeof (hh0s)) ;
      memset (hh1s, 0, sizeof (hh1s)) ;
    }

  gxMakeComparativeHistos (gx, gene, run1, run2, vA, vB, hh0s, hh1s
			   , &ch2, &nRun1, &nRun2, 0, 0, 0
			   , &threshold, 0, 0, isVirtual
			   ) ;
  
  score = gxScoreHistos (gx, gene, rc1, rc2, run1, run2, hh0s, hh1s, ch2, isVirtual, threshold, &sense) ;
  
  ac_free (h) ;
  return score ;
} /* gxScoreOneGeneHisto */
#endif
/*************************************************************************************/
/* filter on minimal delta and sig */
static void gxRuns2geneFilter (GX *gx, PC *pc, float minFoldChange) 
{
  GC *gc ;

  pc->ok = 0 ;
  if (gx->genePlusFileName)
    { 
      if (pc->dIndex >= minFoldChange) 
	pc->ok = 1 ; 
      return ;
    }
  
  /* short genes fluctuate too much to be predictive */
  gc = arrayp (gx->genes, pc->gene, GC)  ;
  if (! gc || (gc->length > 1 && gc->length < 150)) 
    return ;
  
  if (pc->dIndex > minFoldChange && pc->index1 >= 3 && 
      (pc->overlap > .5 || pc->d2dw > .5 || pc->sig > .5)) 
    pc->ok = 1 ; 
} /* gxRuns2geneFilter */

/*************************************************************************************/
/* extract the 'nks' most significant genes */
static void gxRuns2geneGetBest (GX *gx, COMPARE *compare, int isPm)
{
  AC_HANDLE h = ac_new_handle () ;

  int fdrThreshold = compare->fdrThreshold[isPm] ;
  PC *pc ;
  int nMax ;
  KEYSET genePlus = compare->genePlus[isPm] ;
  Array scores = compare->scores[isPm] ;
  int n, gene, geneMax = arrayMax (scores) ;
  Array pp = arrayHandleCreate (geneMax, PC, h) ;

  if (! genePlus)
    genePlus = compare->genePlus[isPm] = keySetHandleCreate (compare->h) ;

  for (gene = n = 0 ; gene < geneMax ; gene++)
    {
      float score = arr (scores, gene, float) ;
      if (score +.5 >= fdrThreshold)
	{
	  GC *gc = arrp (gx->genes, gene, GC) ;
	  if (1 || ! gc->geneGroup)
	    {
	      pc = arrayp (pp, n++, PC) ;
	      pc->gene = gene ;
	      pc->sig = score ;
	    }
	}
    }
  nMax = n ;
  arraySort (pp, pcSigOrder) ;

  if (n)
    for (n = 0, pc = arrp (pp, 0, PC) ; n < nMax ; n++, pc++)
      keySet (genePlus, n) = pc->gene ;

  ac_free (h) ;
  return ;
} /*  gxRuns2geneGetBest */

/*************************************************************************************/
/*************************************************************************************/
/* export the list of the top NN genes for each run and group */
static void gxExportTop (GX *gx, int NN)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (gx->outFileName, ".topGeneList.txt", FALSE, h) ;
  Array gza = arrayHandleCreate (dictMax(gx->geneDict) + 1, GZ, h) ;
  KEYSET ks ;
  int gene, run, runMax = dictMax(gx->runDict) + 1 ;
  Array aa, kss  = arrayHandleCreate (runMax, Array, h) ;
  int i ;
  RC *rc ;
  DC *dc ;
  GC *gc ;
  GZ *gz ;

  gxExportTableHeader (gx, ao, 51) ;

  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
    {
      if (rc->private)
	continue ;
      aa = rc->aa ;
      if (aa)
	{
	  gza = arrayReCreate (gza, dictMax(gx->geneDict) + 1, GZ) ;
	  for (gene = 1, gc = arrp (gx->genes, gene, GC) ; gene <= dictMax (gx->geneDict) ; gc++, gene++)
	    {
	      if (gc->tags < 1) 
		continue ;
	  
	      gz = arrayp (gza, gene, GZ) ;
	      gz->gene = gene ;
	      dc = arrayp (aa, gene, DC) ;
	      gz->index = dc->index ;
	    }
	  arraySort (gza, gzOrder) ;
	  ks = array(kss, run, Array) = arrayHandleCreate (NN, KEY, h) ;
	  for (i = 0 ; i < NN ; i++)
	    keySet (ks, i) = array (gza, i, GZ).gene ;
	}
    }
  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
    if (! rc->private && rc->aa) 
       aceOutf (ao, "\t%s", dictName (gx->runDict, run)) ;
  aceOut (ao, "\n") ;
  for (i = 0 ; i < NN ; i++)
    {
      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	{
	  if (rc->private) 
	    continue ;
	  ks = array(kss, run, Array) ;
	  if (ks)
	    { 
	      gene = keySet (ks, i) ;
	      aceOutf (ao, "\t%s", gene ? noX__ (dictName (gx->geneDict, gene), !gx->isINTRON) : "") ;  
	    }
	}
      aceOut (ao, "\n") ;
    }

  ac_free (h) ;
  return ;
}  /* gxExportTop */

/*************************************************************************************/
/*************************************************************************************/

/* compute the AUC based on the class 0 or 1 declared in gz->gene 
 * U, called the Matt-Whitney-Wicoxon statistics, is by definition the the unormalized AUC
 * auc = U / n1*n2, where n1, n2 are the sizes of the 2 groups, so n1*n2 is the surface of the rectangle
 */
static double gzAUC (Array gza, double *Ginip, double *pValuep) 
{
  GZ *gz ;
  int ir, irMax = arrayMax (gza) ;
  int nn, U, w, n1, n2 ;
  double auc ;

  for (nn = n1 = w = n2 = U = 0, gz = arrp (gza, 0, GZ), ir = 0 ; ir < irMax ; ir++, gz++)
    {
      if (gz->gene == 0)
	{ w = 1 ; n2++ ; nn++; } 
      else
	{ w = 0 ; n1++ ; U += n2 ; nn++ ;  }
    }

  if (w == 1) /*  register the last value */
    n1++ ; 
    
   auc = 0 ; 
   if (n1 * n2 > 0) 
     {
       auc = U ; /* cast to double */
       auc /= (n1 * n2) ; 
     }
   if (Ginip)
     *Ginip = 100 * (2 * auc - 1) ;

   if (pValuep)
     wilcoxon (U, n1, n2, pValuep, 0) ;   /* compute the pValue via the Wilcoxon */

   return 100 * auc ;
} /* gzAUC */

/*************************************************************************************/

static BOOL gxMakeComparativeHistos (GX *gx, int gene
				       , int run1, int run2, KEYSET vA, KEYSET vB
				       , double *hh0s, double *hh1s, float *chi2p
				       , int *nR1p, int *nR2p, int *nHeterop, int *nHomop, int *nObsp
				       , double *thresholdp, double *Ginip, double *pValuep
				       , BOOL isVirtual)
{
  AC_HANDLE h = ac_new_handle () ;
  RC *rc, *rc1, *rc2 ;
  GZ *gz ;
  GC *gc = arrp (gx->genes, gene, GC) ;
  double z, zT, hh0[256], hh1[256], th = 0 ;
  int minIndex, run, i, jj, nr, nr1, nth = 0, nT, gzn = 0 ;
  int nL[2], nH[2] ;  /* low or high values in run1 or 2 => chi2 */
  KEYSET vAB = 0 ;
  Array gza = arrayHandleCreate (arrayMax(gx->runs), GZ, h) ;
  Array gzb = arrayHandleCreate (arrayMax(gx->runs), GZ, h) ;
  BOOL normalize = gx->normalize ;
  float shift = gx->histo_shifting/10.0 ;    

  rc1 = arrayp (gx->runs, run1, RC); 
  rc2 = arrayp (gx->runs, run2, RC); 

  minIndex = (rc1->Low_index > rc2->Low_index ? rc1->Low_index : rc2->Low_index) ;
  if (minIndex < 2) minIndex = 2 ;
  minIndex = 2 ;  /* 2015_02_06, use always, since we now trust the low index range */

  if (nHomop) *nHomop = 0 ;
  if (nHeterop) *nHeterop = 0 ;
  if (nObsp) *nObsp = 0 ;
  if (chi2p) *chi2p = 0 ;
  if (pValuep) *pValuep = 0 ;

  memset (hh0, 0, sizeof(hh0)) ;
  memset (hh1, 0, sizeof(hh1)) ;
  memset (hh0s, 0, sizeof(hh0)) ;
  memset (hh1s, 0, sizeof(hh1)) ;
  for (jj = 0 ; jj < 2 ; jj++)
    {
      nT = 0 ; zT = 0 ;
      vAB = jj == 0 ? vA : vB ;
      if (jj == 0) 
	*nR1p = vA ? keySetMax (vA) : 0 ;
      else
	*nR2p = vB ? keySetMax (vB) : 0 ;

      nL[jj] = nH[jj] = 0 ;
      for (nr = nr1 = 0 ; nr1 < (vAB ? keySetMax (vAB) : 1) ; nr1++)
	{
	  run = keySet (vAB, nr1) ;
	  rc = arrayp (gx->runs, run, RC); 
	  if (rc->avoid) continue ;
	  if (! rc->aa || (! rc->runs && ! rc->hasData)) continue ;  /* skip empty runs */
	  {
	    Array aa = rc->aa ;
	    DC *dc = aa ? arrayp (aa, gene, DC) : 0 ;
	    double z = dc ? dc->index + rc->fineTune - 1000 : 0 ;
	    int iz ;
	    
	    gz = arrayp (gza, gzn,GZ) ;
	    gz->gene = jj ; gz->index = z - (jj ? shift : -shift) ;
	    gz = arrayp (gzb, gzn,GZ) ;
	    gz->gene = jj ; gz->index = z + (jj ? shift : -shift) ;
            gzn++ ;

	    if (dc->isLow) 
	      nL[jj]++ ;
	    else
	      if (z > rc->NA_crossover_index + 2)
		nH[jj]++ ;
	      
	    if (gx->isSNP)
	      {
		if (jj == 0) 
		  *nR1p = dc->kb ;
		else
		  *nR2p = dc->kb ;
	      }		

	    th += rc->NA_crossover_index ;
	    nth++ ;

	    if (gx->isSNP && ! dc->tags)  /* SNP case: seqs=count, tags=cover, do not measured (no coverage) */
	      continue ;

	    if (0)
	      {
		iz = 5 * (z + .1) ;  /* could be 0.05 and 10 below replacing 5 */
		if (iz < 5 * minIndex) iz = 0  ; 
	      }
	    else
	      {
		iz = 10 * (z + .05) ;  /* could be 0.05 and 10 below replacing 5 */
		if (iz < 10 * minIndex) iz = 0  ; 
	      }
	    nr++ ;
	    if (0 && gene == 10524)
	      {
		nT++ ; zT += z ;
		printf("MkHval %s %s   %.2f %.2f\n", dictName (gx->geneDict,gene), dictName (gx->runDict,run), z, zT/nT) ;
	      }
	    if (nObsp && dc && dc->tags > 10) 
	      {
		(*nObsp)++ ;
		if (nHomop && z > 18) (*nHomop)++ ;
		else if (nHeterop && z > 8) (*nHeterop)++ ;
	      }

	    if (1 && ! strncmp (dictName(gx->geneDict, gene),"kloychabu", 8) &&
		( ! strcmp (dictName(gx->runDict, run) , "Rhs5325") 
		  )
		)
	      invokeDebugger () ;

	    if (1) /* perform a rough spreading, followed by finer smoothing below */
	      {
		if (iz < 0) iz = 0 ;
		if (iz > 245) iz = 245 ;
		if (! gx->isSNP || iz == 0 || iz == 245) 
		  {
		    if (jj == 0) hh0[iz]++ ; else hh1[iz]++ ;
		  }
		else if (dc->tags > 10 && gc->intrinsicSigma < 1) /* 2014_01_31, seems silly to over smooth */
		  {
		    if (jj == 0) hh0[iz]++ ; else hh1[iz]++ ;
		  }
		else if (dc->tags > 8 && gc->intrinsicSigma < 2)
		  {
		    int iz1 = iz > 10 ? iz - 10 : 0 ;
		    int iz2 = iz < 235 ? iz + 10 : 245 ;
		    if (jj == 0)
		      {
			hh0[iz]+= .5 ; hh0[iz1]+= .25 ;  hh0[iz2]+= .25 ;
		      }
		    else
		      {
			hh1[iz]+= .5 ; hh1[iz1]+= .25 ;  hh1[iz2]+= .25 ;
		      }
		  }
		else if (dc->tags > 6 && gc->intrinsicSigma < 3)
		  {
		    int iz1 = iz > 15 ? iz - 15 : 0 ;
		    int iz2 = iz < 230 ? iz + 15 : 245 ;
		    if (jj == 0)
		      {
			hh0[iz]+= .5 ; hh0[iz1]+= .25 ;  hh0[iz2]+= .25 ;
		      }
		    else
		      {
			hh1[iz]+= .5 ; hh1[iz1]+= .25 ;  hh1[iz2]+= .25 ;
		      }
		  }
		else 
		  {
		    int iz1 = iz > 20 ? iz - 20 : 0 ;
		    int iz2 = iz < 225 ? iz + 20 : 245 ;
		    if (jj == 0)
		      {
			hh0[iz]+= .5 ; hh0[iz1]+= .25 ;  hh0[iz2]+= .25 ;
		      }
		    else
		      {
			hh1[iz]+= .5 ; hh1[iz1]+= .25 ;  hh1[iz2]+= .25 ;
		      }
		  }
	      }
	  }
	}
      
      gx->nHistoRun = nr ;
      if (nr > 0)  /* normalize  and smooth */
	{
	  int j ;
	  
	  z = 100.0 / nr ;
	  if (! normalize) z = 1 ;

	  if (jj == 0)
	    {
	      float intrinsic_sigma = gc->intrinsicSigma ;
              if (! intrinsic_sigma) 	      /* dromadaire 2015_06_27 : for genes with unknown intrinsicSigma use value 1 not .5 */
		intrinsic_sigma = gx->isMA ? .5 : .5 ; /* using 1 lower EF, does not change HR */
	      for(i = 0 ; i < 256 ; i++)
		hh0[i] *= z ; 
	      hh0s[0] = 0 * hh0[0] ;
	      if (0) 
		{
		  for (i = 0 ; i <= 255 ; i++)
		    hh0s[i] += hh0[i] ;
		}
	      else if (gc->intrinsicSigma <= .5)
		{
		  for (i = 0 ; i <= 255 ; i++)
		    for (j = -10 ; j <= 10 ; j++)
		      hh0s[i] += hh0[(i+j < 0 ? 0 : ( i+j < 255 ? i+j : 255))] * gx->gauss5[j>=0 ? j : -j] ;
		}
	      else if (gc->intrinsicSigma < 1) 
		{
		  for (i = 0 ; i <= 255 ; i++)
		    for (j = -10 ; j <= 10 ; j++)
		      hh0s[i] += hh0[(i+j < 0 ? 0 : ( i+j < 255 ? i+j : 255))] * gx->gauss10[j>=0 ? j : -j] ;
		}
	      else if (gc->intrinsicSigma < 2)
		{
		  for (i = 0 ; i <= 255 ; i++)
		    for (j = -20 ; j <= 20 ; j++)
		      hh0s[i] += hh0[(i+j < 0 ? 0 : ( i+j < 255 ? i+j : 255))] * gx->gauss20[j>=0 ? j : -j] ;
		}
	      else if (gc->intrinsicSigma < 3)
		{
		  for (i = 0 ; i <= 255 ; i++)
		    for (j = -30 ; j <= 30 ; j++)
		      hh0s[i] += hh0[(i+j < 0 ? 0 : ( i+j < 255 ? i+j : 255))] * gx->gauss30[j>=0 ? j : -j] ;
		}
	      else 
		{
		  for (i = 0 ; i <= 255 ; i++)
		    for (j = -40 ; j <= 40 ; j++)
		      hh0s[i] += hh0[(i+j < 0 ? 0 : ( i+j < 255 ? i+j : 255))] * gx->gauss40[j>=0 ? j : -j] ;
		}
	    }
	  if (jj == 1)
	    {
	      for (i = 0 ; i < 256 ; i++)
		hh1[i] *= z ;
	      hh1s[0] = 0 *  hh1[0] ;
	      if (0)
		{
		  for (i = 0 ; i <= 255 ; i++)
		    hh1s[i] += hh1[i] ;
		}
	      else if (gc->intrinsicSigma <= .5)
		{
		  for (i = 0 ; i <= 255 ; i++)
		    for (j = -10 ; j <= 10 ; j++)
		      hh1s[i] += hh1[(i+j < 0 ? 0 : ( i+j < 255 ? i+j : 255))] * gx->gauss5[j>=0 ? j : -j] ;
		}
	      else if (gc->intrinsicSigma < 1)
		{
		  for (i = 0 ; i <= 255 ; i++)
		    for (j = -10 ; j <= 10 ; j++)
		      hh1s[i] += hh1[(i+j < 0 ? 0 : ( i+j < 255 ? i+j : 255))] * gx->gauss10[j>=0 ? j : -j] ;
		}
	      else if (gc->intrinsicSigma < 2)
		{
		  for (i = 0 ; i <= 255 ; i++)
		    for (j = -20 ; j <= 20 ; j++)
		      hh1s[i] += hh1[(i+j < 0 ? 0 : ( i+j < 255 ? i+j : 255))] * gx->gauss20[j>=0 ? j : -j] ;
		}
	      else if (gc->intrinsicSigma < 3)
		{
		  for (i = 0 ; i <= 255 ; i++)
		    for (j = -30 ; j <= 30 ; j++)
		      hh1s[i] += hh1[(i+j < 0 ? 0 : ( i+j < 255 ? i+j : 255))] * gx->gauss30[j>=0 ? j : -j] ;
		}
	      else 
		{
		  for (i = 0 ; i <= 255 ; i++)
		    for (j = -40 ; j <= 40 ; j++)
		      hh1s[i] += hh1[(i+j < 0 ? 0 : ( i+j < 255 ? i+j : 255))] * gx->gauss40[j>=0 ? j : -j] ;
		}
	    }
	}
    }
  if (1)
    {
      double u = 0, v = 0 ;
      z = 100.0 ;
      if (! normalize) z = 10 * nr ;

      for(i = 0 ; i < 256 ; i++)
	{
	  u += hh0s[i] ;
	  v += hh1s[i] ;
	}
      
      if (u > 0)
	for(i = 0 ; i < 256 ; i++)
	  hh0s[i] *= z/u ;
      if (v > 0)
	for(i = 0 ; i < 256 ; i++)
	  hh1s[i] *= z/v ;
      if (gene == 50331)
	fprintf (stderr, "### %s histo area %.2g %.2g\n", dictName (gx->geneDict, gene), u, v) ;
    }

  if (0 && gene == 10524)
    {
      float xt0 = 0.1, xt1 = 0.1 ; float zt0 = 0, zt1 = 0 ;
      for(i = 0 ; i < 256 ; i++)
	{
	  zt0 += i *  hh0s[i] ; xt0 +=   hh0s[i] ;
	  zt1 += i *  hh1s[i] ; xt1 +=   hh1s[i] ;
	  printf ("MkHhs %s %s %s     %d %.1f %.1f    %.1f %.1f\n"
		  , dictName (gx->geneDict,gene)
		  , dictName (gx->runDict,run1)
		  , dictName (gx->runDict,run2)
		  , i, hh0s[i], hh1s[i], zt0/xt0, zt1/xt1) ;
	}
    }

  if (thresholdp)
    *thresholdp = nth ? th/nth : 0 ;
  if (chi2p)   /* find the extreme genes high (>threshold+2) in nH cases, isLow in nL cases */
    chi2(nL[0], nL[1], nH[0], nH[1], chi2p) ;

  arraySort (gza, gzOrder) ;
  
  if (Ginip)
    { /* select worse scenario */
      double p1 = 1, p2 = 1, g1 = 0, g2 = 0 ;

      arraySort (gza, gzOrder) ;
      arraySort (gzb, gzOrder) ;
      gzAUC (gza, &g1, &p1) ;
      gzAUC (gzb, &g2, &p2) ;
      
      if (g1 * g1 < g2 * g2)
	*Ginip = g1 ;
      else
	*Ginip = g2 ;
      if (p1 >= 0 && p1 <= p2) /* 1 is best, select 2 */
	{
	  *Ginip = g2 ;
	  if (pValuep) *pValuep = p2 ;
	}
      if (p2 >= 0 && p2 < p1) /* 2 is best, select 1 */
	{
	  *Ginip = g1 ;
	  if (pValuep) *pValuep = p1 ;
	}
    }

  ac_free (h) ; 
  return TRUE ;
} /* gxMakeComparativeHistos */

/*************************************************************************************/
/* select genes above a given differential and export them 
 * optionally count the FDR on a set of virtual randomizedly stratified control classes
 */
static Array gxDoCompare (GX *gx, COMPARE *compare
			     , int run1, int run2, KEYSET vA, KEYSET vB
			     , int isVirtual, int isPm
			     )
{
  AC_HANDLE h = ac_new_handle () ;
  int k, gene, ipp, aaMax ;
  DC *dc1, *dc2 ;
  RC *rc1 = 0, *rc2 = 0 ;
  PC *pc ;
  GC *gc ;
  RC reserveRc1, reserveRc2 ;
  Array aa1, aa2, pp = 0 ;
  double ix1, ix2, dIndex, dvar0 ;
  float minFoldChange ;
  float minIndex = 0 ;

  if (isVirtual)
    {
      rc1 = arrp (gx->runs, run1, RC) ;
      if (rc1->run != run1) messcrash ("mixup in gxDoCompare run1=%d:%s != rc1->run=%d:%s"
				       , run1, dictName (gx->runDict, run1)
				       , rc1->run, dictName (gx->runDict, rc1->run)
				       ) ;
      reserveRc1 = *rc1 ;
      rc1->cumulDone = FALSE ; /* rc1->addCounts = FALSE ; */
      rc1->fineTune = 0 ;
      rc1->aa = rc1->daa = rc1->aag = 0 ;
      rc1->runs = vA ;

      rc2 = arrp (gx->runs, run2, RC) ;
      if (rc2->run != run2) messcrash ("mixup in gxDoCompare run2=%d:%s != rc2->run=%d:%s"
				       , run2, dictName (gx->runDict, run2)
				       , rc2->run, dictName (gx->runDict, rc2->run)
				       ) ;
      reserveRc2 = *rc2 ;
      rc2->cumulDone = FALSE ; /* rc2->addCounts = FALSE ; */
      rc2->fineTune = 0 ;
      rc2->aa = rc2->daa = rc2->aag = 0 ;
      rc2->runs = vB ;

      gxOneGroupCumul (gx, rc1, 0) ; 
      gxOneGroupCumul (gx, rc2, 0) ; 
    }

  rc1 = arrayp (gx->runs, run1, RC) ;
  aa1 = rc1->aa ;
  
  rc2 = arrayp (gx->runs, run2, RC) ;
  aa2 = rc2->aa ;
  if (!aa1 || ! aa2)
    return 0 ;

  minIndex = rc1->NA_crossover_index ;
  if (minIndex < rc2->NA_crossover_index)
    minIndex = rc2->NA_crossover_index ; 

  ipp = 0 ;
  pp = arrayHandleCreate (10000, PC, h) ;
  aaMax = arrayMax (aa1) ; 
  if (aaMax > arrayMax (aa2))
    aaMax = arrayMax (aa2) ;
	  
  minFoldChange = gx->minFoldChange ;
  if (minFoldChange < 1)
    {
      float u1, u2 ;
      u1 = rc1->runs ? keySetMax (rc1->runs) : 0 ; 
      u2 = rc2->runs ? keySetMax (rc2->runs) : 0 ; 
      if (u1 < 4 || u2 < 4)
	minFoldChange = 1.0 ;
      else if (u1 + u2 < 8 &&	minFoldChange < 1 &&
	       minFoldChange < 1.5 - 1.1 * (u1+u2)/8.0
	       )
	minFoldChange = 1.5 - 1.1 * (u1+u2)/8.0 ;  /* interpolate down to .4 at or over 50 */
    }
  
  for (gene = 1, dc1 = arrayp (aa1, 1, DC), dc2 = arrayp (aa2, 1, DC) ; gene < aaMax ; gene++, dc1++, dc2++) 
    {
      gc = arrayp (gx->genes, gene, GC)  ;

      if (gx->targeted && ! gc->targeted)
	continue ;
      if (0)
	if (gx->captured && ! gc->captured)
	  continue ;
      if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
	continue ;
      
      if (! gc || (gc->length > 0 && gc->length < 150))   /* short genes fluctuate too much to be predictive */
	continue ;
      if (0 && gc->geneGroup) /* 2019_05_30 */
	continue ;
      
      ix1 = dc1->index + rc1->fineTune - 1000 ; if (! gx->isSNP && ix1 < minIndex) ix1 = minIndex ;
      ix2 = dc2->index + rc2->fineTune - 1000 ; if (! gx->isSNP && ix2 < minIndex) ix2 = minIndex ;
      dIndex = ix1 - ix2 ;
      if (dc1->isLow)    /* 2014_06_04 22h54 */
	continue ;
      /*
	if (1 && ! strcasecmp (dictName (gx->geneDict, gene), "UKv4_A_23_P2258"))
	invokeDebugger () ;
	if (0 && ! strcasecmp (dictName (gx->geneDict, gene), "NBPF8.2.voAug10-unspliced"))
	invokeDebugger () ;
	if (!strcmp (dictName(gx->geneDict, gene), "EEF2K.aAug10:4206:A2G:GATGCTaGGTGCG:GATGCTgGGTGCG") && ! strcmp (dictName(gx->runDict, run1), "Unfavorable_NB_91"))
	invokeDebugger () ;
      */
      
      if (ix1 <= 0 || dIndex < minFoldChange)
	continue ;
      
      if (! strncasecmp (dictName (gx->geneDict, gene), "ERCC-", 5))
	continue ;
      if (! strncasecmp (dictName (gx->geneDict, gene), "ERCC_vector", 11))
	continue ;
      
      gc = arrp (gx->genes, gene, GC) ;
      
      pc = arrayp (pp, ipp++, PC) ;
      pc->gene = gene ;
      pc->index1 = ix1 ; pc->index2 = ix2 ; pc->dIndex = dIndex ;
      pc->dvar = pc->d2dw = pc->sig = pc->overlap = pc->danielle = 0 ;
      
#ifdef JUNK
      if (0) /* divariance */
	{
	  pc->var1 = dc1->divariance ;  pc->var2 = dc2->divariance ; 
	  
	  if (                    /* dc->N is the number of significant values used to compute the divariance */
	      (                  
	       dc2->N > 0 && 
	       (rc1->nRuns < 40 || rc2->nRuns < 40) &&
	       dc1->indexMin < dc2->indexMax
				 ) ||
	      (1 * rc1->nRuns > 4 * dc1->N)
				  ) ;
	  else
	    {
	      if (dc2->N)
		{
		  if (pc->var1 < 0 && pc->var2 >= 0) pc->dvar = 2 * pc->var2 ;
		  else if (pc->var1 >= 0 && pc->var2 < 0) pc->dvar = 2 * pc->var1 ;
		  else if (pc->var1 < 0 && pc->var2 < 0) pc->dvar = 0 ;
		  else
		    pc->dvar = pc->var1 + pc->var2 ; /* sqrt (pc->var1 * pc->var1 + pc->var2 * pc->var2) */
		  
		  pc->d2dw = rc1->nRuns - dc1->N ;
		  if (0 && dc1->indexMin > dc2->indexMax) /* disfavor crossing histograms */
		    pc->d2dw += 10000 ;
		}
	      else
		pc->dvar = pc->var1 ;
	    }
	}
#endif
      if (1)  /* 2012_08_13: Danielle's new method based on the union and intersection of the histograms of expression */
	{
	  /* for a given gene, loop on all patients and compute the hitogram in the plus minus groups
	   * norm the histo by dividing by the number of runs in each group
	   * measure the surface of the union and intersect histo = histo of max/min of the 2 values
	   * sig = 1 - intersect/union will vary between 1: disconnected histo and zero identical histos
	   */
	  int i, nRun1 = 0, nRun2 = 0 ;
	  double z, s0, s1, *a0p, *a1p, mx0, mx1, Gini = 0 ;
	  float ch2 = 0 ;
	  double hh0s[256], hh1s[256] ;
	  /*  int danRegion0 = 20, danRegion1 = 70, danRegion2 = 80, danRegion3 = 120 ; */
	  double threshold = 0 ;
	  gxMakeComparativeHistos (gx, gene, run1, run2, vA, vB, hh0s, hh1s
				   , &ch2, &nRun1, &nRun2, 0, 0, 0
				   , &threshold, &Gini, 0, isVirtual
				   ) ;
	  
	  /* find the max in the co &pValuentinuous part */
	  for (s0 = s1 = 0, mx0 = mx1 = 0, i = 1, a0p = hh0s, a1p = hh1s ; i < 256 ; a0p++, a1p++, i++)
	    {
	      if (*a0p > mx0) mx0 = *a0p ;
	      if (*a1p > mx1) mx1 = *a1p ;
	      
	      if (*a0p < *a1p) 
		{ s0 += *a0p ; s1 += *a1p ; }
	      else
		{ s0 += *a1p ; s1 += *a0p ; }
	    }
	  z = 1 - s0/s1 ; /* old overlap method */
	  pc->overlap = z ;

	  if (gx->geneSelectionMethod >= 7)     /* 2012_09_8 */
	    {
	      int sense ;
	      
	      pc->danielle = gxScoreHistos (gx, gene, rc1, rc2, run1, run2, hh0s, hh1s, ch2, isVirtual, threshold, &sense) ;
	      if (0 && pc->danielle > 0 && sense != 2) /* TOTOTOTO */
		pc->danielle = 0 ;
	      if (gx->geneSelectionMethod == 8)
		pc->danielle = 2 * Gini ; /* Gini == 100 is a perfect score */
	      
	    }
	  else if (0)
	    gxScoreByOverlap() ;
	  
	  pc->geneid = gc->geneId ;
	  if (0 && ! strcasecmp (dictName (gx->geneDict, gene), "UKv4_A_23_P2258"))
	    invokeDebugger () ;
	}     
    }

  if (! ipp)
    goto done ;
  /* compute the median of all divariances */
  if (0)
    {
      arraySort (pp, pcDvarOrder) ;
      pc = arrayp (pp, (arrayMax (pp)/2), PC) ;
      dvar0 = pc->dvar ? .66 * pc->dvar : 1.0 ;

      for (ipp = 0, pc = arrp (pp, 0, PC) ; ipp < arrayMax (pp) && pc->dvar <= 0 ; ipp++, pc++) ;
      if (ipp < arrayMax (pp)/2)
	{
	  pc = arrayp (pp, (arrayMax (pp) + ipp)/2, PC) ;
	  dvar0 = pc->dvar ? .66 * pc->dvar : 1.0 ;
	}
      else
	dvar0 = 1 ;

      /* damp the variance by adding the median and recompute the significativity */
      if (gx->geneSelectionMethod < 7)
	for (ipp = 0, pc = arrp (pp, 0, PC) ; ipp < arrayMax (pp) ; ipp++, pc++)
	  {
	    double z = 0 ;
	    int nr = arrayMax (vA) ; /* nb of runs in this group */
	    
	    /* each missing value contributes a 3*dvar0/nr */
	    if (pc->d2dw >= 10000) { z = dvar0 ; pc->d2dw -= 10000 ; }
	    z += dvar0 * ( 1 + 1 * 3 * pc->d2dw/nr) + (pc->dvar ? pc->dvar : dvar0) ;
	    
	    z = z > 0 ? pc->dIndex/z : 0 ;
	    
	    pc->d2dw = z ;
	  }
    }

  for (ipp = k = 0, pc = arrp (pp, 0, PC) ; ipp < arrayMax (pp) ; ipp++, pc++)
    {
      switch (gx->geneSelectionMethod)
	{
	case 1: pc->sig = pc->d2dw ; break ;
	case 2: pc->sig = 1 ; break ;
	case 3: pc->sig = pc->overlap ; break ;
	case 4: pc->sig = pc->d2dw + pc->overlap - .1 ; break ;
	case 5: pc->sig = pc->danielle ; break ;
	case 7: pc->sig = pc->danielle ; break ;
	case 8: pc->sig = pc->danielle ;  break ; /* use the Gini coefficient */
	}
      
      gxRuns2geneFilter (gx, pc, minFoldChange) ;
      if (pc->ok)
	{
	  if (pc->sig > 180) k++ ;
	  if (0 && pc->gene == 44002)
	    fprintf (stderr, "Gene %d run1=%d run2=%d isPm=%d isVirtual=%d sig=%.2f\n", pc->gene, run1, run2, isPm, isVirtual, pc->sig) ;
	  gxDghRegisterScore (gx, isPm ? isPm - 1 : isVirtual - 1 + 2, pc) ;
	}
      /*       pc->ok = 0 ; pc->sig = 0 ; */
    }
  
  fprintf (stderr, "run1=%d run2=%d IVV=%d isPm=%d pp->max = %d  ok=%d\n", run1, run2, isVirtual, isPm, arrayMax (pp), k) ;
  if (isPm == 2)
    compare->noise = arrayHandleCreate (dictMax(gx->geneDict) +1, float, compare->h) ;

  if (isPm)
    {  
      int nRun1 = rc1 && rc1->runs ? arrayMax (rc1->runs) : 0 ;
      int nRun2 = rc2 && rc2->runs ? arrayMax (rc2->runs) : 0 ;
      
      for (int i = 1 ; i<= 2 ; i++)
	if (! compare->scores[i])
	  compare->scores[i] = arrayHandleCreate (dictMax(gx->geneDict) +1, float, compare->h) ;

      arraySort (pp, pcSigOrder) ;

      if (isPm == 2)
	{
	  gxDghShiftScores (gx, compare) ;
	  gxDghRegisterHisto (gx, run2, run1, nRun2, nRun1) ;
	  gxDghExportHistos (gx, compare, run2, run1, 1) ;
	  gxDghExportHistos (gx, compare, run1, run2, 2) ;

	  gxRuns2geneGetBest (gx, compare, 1) ;
	  gxRuns2geneGetBest (gx, compare, 2) ;

	  if (0) /* should be used in actual pair compare not in the virtual case */
	    gxExportDiffGenes (gx, compare, pp, run1, run2, vA, vB, isPm) ;
	}
    }

 done:    
  if (isVirtual)
    {
      rc1 = arrp (gx->runs, run1, RC) ;
      arrayDestroy (rc1->aa) ;
      arrayDestroy (rc1->daa) ;
      arrayDestroy (rc1->aag) ;
      *rc1 = reserveRc1 ;
      if (rc1->run != run1) messcrash ("mixup in gxDoCompare run1=%d:%s != rc1->run=%d:%s"
				       , run1, dictName (gx->runDict, run1)
				       , rc1->run, dictName (gx->runDict, rc1->run)
				       ) ;
      
      rc2 = arrp (gx->runs, run2, RC) ;
      arrayDestroy (rc2->aa) ;
      arrayDestroy (rc2->daa) ;
      arrayDestroy (rc2->aag) ;
      *rc2 = reserveRc2 ;
      if (rc2->run != run2) messcrash ("mixup in gxDoCompare run2=%d:%s != rc2->run=%d:%s"
				       , run2, dictName (gx->runDict, run2)
				       , rc2->run, dictName (gx->runDict, rc2->run)
				       ) ;
    }
  
  ac_free (h) ;
  return pp ;  
} /* gxDoCompare */

/*************************************************************************************/
/*************************************************************************************/
typedef struct dghStruct { KEY run1, run2 ; int iVV, nRun1, nRun2 ; double histo[201] ;} DGH ;

/*************************************************************************************/
#ifdef JUNK
static int dghOrder (const void *a, const void *b)
{
  const DGH *up = (const DGH *)a, *vp = (const DGH *)b ;
  int n ;

  n = up->iVV - vp->iVV ; if (n) return n ;
  n = up->run1 - vp->run1 ; if (n) return n ;
  n = up->run2 - vp->run2 ; if (n) return n ;


  return 0 ;
} /* dghOrder */
#endif

/*************************************************************************************/
/*
  baddybatchy array, sum up the score of each gene seemingly significant in the virtual strata iVV > 0
*/
static int gxDghInitializeBaddyBatchies (GX *gx, AC_HANDLE h)
{
  Array bb = gx->baddyBatchies ;
  Array bbAll = gx->baddyBatchiesAll ;
  int nVV = 2 + NVIRTUALSTRATA ; /* desired number of virtual strata */ 

  if (!bb)
    {
      bb = gx->baddyBatchies = arrayHandleCreate (arrayMax (gx->genes), float, gx->h) ;
      bbAll = gx->baddyBatchiesAll = arrayHandleCreate (nVV * arrayMax (gx->genes), float, gx->h) ;
    }
  else
    {
      bb = gx->baddyBatchies = arrayReCreate (bb, arrayMax (gx->genes), float) ;
      bbAll = gx->baddyBatchiesAll = arrayReCreate (bbAll, nVV * arrayMax (gx->genes), float) ;
    }
  return 0 ;
} /* gxDghInitializeBaddyBatchies */

/*************************************************************************************/
/* iVV = 0 : signal  run1 > run2,  
 * iVV = 1 : signal  run1 < run2,  
 * iVV >= 2 : random strata
 */
static void gxDghRegisterScore (GX *gx, int iVV, PC *pc)
{
  Array bb = gx->baddyBatchies ;
  Array bbAll = gx->baddyBatchiesAll ;
  int nVV = NVIRTUALSTRATA + 2 ; /* desired number of virtual strata */ 

  if (!bbAll)
    messcrash ("Bad call to gxDghRegisterScore") ;

  array (bbAll, nVV * pc->gene + iVV, float) = pc->sig ;
  if (iVV >= 2)
    array (bb, pc->gene, float) += pc->sig ;

  return ;
} /* gxDghRegisterScore */

/******************************************/
/* for each gene and each random strata, remove the gene noise value 
 * then recompute the signal and the noise histos
 */
static void gxDghShiftScores (GX *gx, COMPARE *cmp)
{
  Array bb = gx->baddyBatchies ;
  Array bbAll = gx->baddyBatchiesAll ;
  int gene, gMax = arrayMax (bb) ;
  int iVV, nVV = 2 + NVIRTUALSTRATA ; /* desired number of virtual strata */ 



  for (gene = 0 ; gene < gMax ; gene++)
    {
      float z, shift ;

      shift = array (bb, gene, float)/nVV ;
      if (0) shift = 0 ;
      array (cmp->noise, gene, float) = shift ;

      for (iVV = 0 ; iVV < nVV ; iVV++)
	{
	  z = array (bbAll, nVV * gene + iVV, float) ;
	  z -= shift ;
	  if (z < 0) z = 0 ;
	  array (bbAll, nVV * gene + iVV, float) = z ;

	  if (iVV < 2)
	    array (cmp->scores[iVV+1], gene, float) = z ;
	}
    }

  return ;
} /* gxDghShiftScores */

/*************************************************************************************/


static void gxDghRegisterHisto (GX *gx, int run1, int run2, int nRun1, int nRun2)
{
  Array bbAll = gx->baddyBatchiesAll ;
  int nVV = 2 + NVIRTUALSTRATA ; /* desired number of virtual strata */ 
  int kk = 0 ;

  if (!bbAll)
    messcrash ("Bad call to gxDghRegisterHisto") ;

  for (int iVV = 0 ; iVV < nVV ; iVV += 1)
    {
      int jj = arrayMax (gx->degHistos) ;
      DGH *dgh = arrayp (gx->degHistos, jj, DGH) ;
      int gMax = arrayMax (gx->genes) ;

      dgh->run1 = run1 ;
      dgh->run2 = run2 ;
      dgh->nRun1 = nRun1 ;
      dgh->nRun2 = nRun2 ;
      
      dgh->iVV = iVV ;
      double *hh = dgh->histo ;

      for (int i = 0 ; i <= 200 ; i++)
	hh[i] = 0 ;

      /* compute the histogram of gene index rounding by int value */
      for (int g = 0 ; g < gMax ; g++)
	{
	  GC *gc = arrayp (gx->genes, g, GC)  ;
	  if (! gc->geneGroup)
	    {
	      float z = array (bbAll, nVV * g + iVV, float) ;
	      int n = z  + .49 ;
	      
	      if (n > 200) n = 200 ;
	      if (n > 0)
		hh[n]++ ;
	      if (1 && run1 == 15 && run2 == 16 && n >= 156 && iVV == 1)
		{
		  kk++ ;
		  if (kk < 16)
		    fprintf(stderr, "==== score %d  kk=%d %s\n", n, kk, dictName (gx->geneDict, g)) ;
		}
	    }
	}
    }
  fprintf(stderr, "==== score 185  %d genes \n", kk) ;
  return ;
} /* gxDghRegisterHisto */

/******************************************/
/* adjust the FDR threhshold */
static void gxDghExportHistos (GX *gx, COMPARE *compare, int run1, int run2, int isPm)
{
  AC_HANDLE h = ac_new_handle () ;
  Array degHistos = arrayHandleCreate (arrayMax(gx->degHistos)+1, DGH, h) ;
  static ACEOUT ao = 0 ;
  DGH *dgh ;
  int i, ii, iiMax, jj, nn ;
  int  nRun1 = 0, nRun2 = 0 ;
  int  nCRun1 = 0, nCRun2 = 0 ;
  BOOL hasControls = FALSE ;
  double z, zn, GAUSS ;
  double cE[201] ;   /* smoothed experiment cumul */
  double aE[201] ;   /* smoothed experiment */
  static double aExp1[201] ;   /* copy of smoothed experiment 1 */
  double cE0[201] ;  /* cumul before smoothing */
  static double cX[201] ;   /* smoothed average resampling cumul */
  static double aX[201] ;   /* smoothed average resampling */
  static double aX0[201] ;  /* before smoothing */
  float fdr[201] ;
  static int myDEG[3] ;
  static int myDEG185[3] ;
  float *myFdr = compare->fdr ;
  int *myFdrThreshold = compare->fdrThreshold ;

  /* keep both orientation */ 
  compare->nVirtualStrata = 0 ;
  if (arrayMax(gx->degHistos) < 1)
    goto done ;
  else
    {  /* make space */
      /* copy */
      for (ii = jj = 0, dgh = arrp (gx->degHistos, ii, DGH) ; ii < arrayMax (gx->degHistos) ; ii++, dgh++)
	if (1) /* dgh->run1 == run1 && dgh->run2 == run2) */
	  {
	    DGH *dgh1 = arrayp (degHistos, jj++, DGH) ;
	    *dgh1 = *dgh ;
	    
	    if (dgh->iVV >= 2)
	      {
		compare->nVirtualStrata++ ;
		hasControls = TRUE ;
		nCRun1 = dgh->nRun1 ;
		nCRun2 = dgh->nRun2 ;
	      }
	  }
    }
  
  /* clean up */
  if (isPm == 1)
    for (int i = 0 ; i < 3 ; i++)
      {
	myFdr[i] = 0 ;
	myFdrThreshold[i] = 0 ;
      }
  
  if (1)
    {
      RC *rc1 = arrayp (gx->runs, run1, RC) ;
      RC *rc2 = arrayp (gx->runs, run2, RC) ;
      nRun1 = rc1 && rc1->runs ? arrayMax (rc1->runs) : 0 ;
      nRun2 = rc2 && rc2->runs ? arrayMax (rc2->runs) : 0 ;
    }

  if (isPm == 1)
    {
      ao = aceOutCreate (gx->outFileName
			 , messprintf(".%s.%s.%s.DEG_FDR_selection.txt"
				      , dictName (gx->pairDict, compare->compare)
				      , dictName (gx->runDict, run1)
				      , dictName (gx->runDict, run2)
				      )
			 , FALSE, 0) ;
      aceOutDate (ao, "###", gx->title) ;
      
      
      aceOutf (ao, 
	       "## Histogram of the number of differentially expressed gene per score in different experimental conditions.\n" 
	       "## For each experiment, and each method: RNA_seq AceView, RefSeq, EBI, genes or transcripts ; Agilent microarray, RNA-seq at probe locations\n" 
	       "## the histogram of the number of genes or transcripts or probe with a given differential score is presented. The score are in the range [0,200].\n"
	       "## In the same conditions, the test is computed by randomly resampling %d times the same samples in 2 pseudo-phenotypic classes alpha and beta,\n"
	       "## attributing the same proportion of A and B phenotypes to the alpha and beta classes, so the alpha-beta split is random and orthogonal to the A-B split.\n" 
	       "## and can be used to estimate the noise that occurs independently of the experimental conditions\n"
	       "## The control histogram corresponds to the average of the %d resamplings.\n"
	       "## The false discovery rate (FDR) in each bin is estimated by dividing the number of controls by the number of observed differential genes.\n" 
	       "## The estimator (Theta) is used to find the region with lowest relative noise\n"
	       "### %s\n"
	       , compare->nVirtualStrata , compare->nVirtualStrata
	       , dictName (gx->pairDict, compare->compare)
	       ) ;
    }

  /* cE  smoothed experiment cumul 
   * aE  smoothed experiment 
   * cE0 cumul before smoothing 
   * cX  smoothed average resampling cumul 
   * aX  smoothed average resampling 
   * aX0 before smoothing 
   */

  iiMax = arrayMax (degHistos) ;

  memset (fdr, 0, sizeof (fdr)) ;
  memset (cE, 0, sizeof (cE)) ;
  memset (aE, 0, sizeof (aE)) ;
  memset (cE0, 0, sizeof (cE0)) ;
  if (isPm == 1)
    {
      memset (cX, 0, sizeof (cX)) ;
      memset (aX, 0, sizeof (aX)) ;
      memset (aX0, 0, sizeof (aX0)) ;
    }

  for (nn = ii = 0 ; ii < iiMax ; ii++)
    {
      dgh = arrp (degHistos, ii, DGH)  ;

      if (dgh->iVV == isPm - 1)   /* cumul of the experiment */
	{
	  for (zn = 0, i = 200 ; i >= 0 ; i--)
	    {
	      z = dgh->histo[i] ;
	      zn += z ;
	      cE0[i] = zn ;     
	      if (0 && run2 == 15 && run1 == 16 && i >= 185 && dgh->iVV == 1)
		fprintf(stderr, "++++++ h[%d]=%.0f cumul=%.0f\n", i, z, zn) ;
	    }
	}

      if (isPm == 1 && dgh->iVV >= 2)       /* sum of the resampling */
	{
	  nn++ ;
	  for (i = 0 ; i <= 200 ; i++)
	    aX0[i] += dgh->histo[i] ;
	}
    }

  if (isPm == 1 && nn)   /* average the resampling */
    {
      for (i = 0 ; i <= 200 ; i++)
	aX0[i] /= nn ;
    }

  /* smooth , this multiplies by GAUSS all values, without changing the integral */
  GAUSS = 64 ;
  if (GAUSS == 64)
    for (ii = 0 ; ii < arrayMax (degHistos) ; ii++)
      {
	double xx[210] ;
	memset (xx, 0, sizeof(xx)) ;
	
	dgh = arrp (degHistos, ii, DGH)  ;

	if (isPm == 1 && dgh->iVV == 1)
	  continue ;
	if (isPm == 2 && dgh->iVV != 1)
	  continue ;

	for (i = 0 ; i <= 200 ; i++)
	  xx[i+5] = dgh->histo[i] ;
	for (i = 0 ; i <= 200 ; i++)
	  dgh->histo[i] = xx[i-3+5] + 6 * xx[i-2+5] + 15 * xx[i-1+5] + 20 * xx[i+5] + 15 * xx[i+1+5] + 6 * xx[i+2+5] + xx[i+3+5] ;
	/* rescale the reminders so that the integral is invariant */
	dgh->histo[0]   += 22 * xx[0+5] ;
	dgh->histo[200] += 22 * xx[200+5] ;
	dgh->histo[1]   +=  7 * xx[1+5] ;
	dgh->histo[199] +=  7 * xx[199+5] ;
	dgh->histo[2]   +=  1 * xx[2+5] ;
	dgh->histo[198] +=  1 * xx[198+5] ;
	for (i = 0 ; i <= 200 ; i++)
	  dgh->histo[i] /= 64.0 ;
      }

  for (nn = ii = 0 ; ii < iiMax ; ii++)
    {
      dgh = arrp (degHistos, ii, DGH)  ;

      if (dgh->iVV == isPm - 1)   /* smoothed experiment */
	{
	  for (zn = 0, i = 200 ; i >= 0 ; i--)
	    {
	      z = dgh->histo[i] ;
	      zn += z ;
	      cE[i] = zn ;    /* smoothed cumul */
	      aE[i] = z ;     /* smoothed values */
	      if (0 && run2 == 15 && run1 == 16 && i >= 185 && dgh->iVV == 1)
		fprintf(stderr, "++++++ h[%d]=%.0f smoothed cumul=%.0f\n", i, z, zn) ;	    }
	}

      if (isPm == 1 && dgh->iVV >= 2)       /* smoothed sum of the resampling */
	{
	  nn++ ;
	  for (i = 0 ; i <= 200 ; i++)
	    aX[i] += dgh->histo[i] ;
	}
    }

  if (isPm == 1)
    {
      for (i = 0 ; i <= 200 ; i++)
	aExp1[i] = aE[i] ;

      /* average the smoothed resampling */
      for (i = 0 ; i <= 200 ; i++)
	aX[i] /= nn ;
      /* cumul the smoothed resampling */
      for (zn = 0, i = 200 ; i >= 0 ; i--)
	{
	  z = aX[i] ;
	  zn += z ;
	  cX[i] = zn ;
	}
    }


  /* export */
  /* export the coordinates */ 
  if (isPm == 1)
    {
      aceOutf (ao, "# Compare\tGene over expressed in\trelative to\tNumber of sample in 1st group\tNumber of samples in 2nd group\tDEG_score") ;
      for (i = 0 ; i <= 200 ; i++)
	aceOutf (ao, "\t%d", i) ;
      aceOut (ao, "\n") ;
    }

  /* export the basics */
  if (ao && isPm)
    {
      double best = 0 ; 
      int iBest = 200 ;

      /* export the cumul of the experiment */
      aceOutf (ao, "%s", dictName (gx->pairDict, compare->compare)) ;
      aceOutf (ao, "\t%s", dictName (gx->runDict, run1)) ;
      aceOutf (ao, "\t%s", dictName (gx->runDict, run2)) ;
      aceOutf (ao, "\t%d\t%d", nRun1, nRun2) ;
      aceOutf (ao, "\tSmoothed cumul of experiment") ;
      for (i = 0 ; i <= 200 ; i++)
	aceOutf (ao, "\t%.2f", cE0[i]) ;
      aceOut (ao, "\n") ;
      
      /* export the smoothed control */
      aceOutf (ao, "%s", dictName (gx->pairDict, compare->compare)) ;
      if (isPm == 1)
	{
	  aceOutf (ao, "\t%s", dictName (gx->runDict, run1)) ;
	  aceOutf (ao, "\t%s", dictName (gx->runDict, run2)) ;
	  aceOutf (ao, "\t%d\t%d", nCRun1, nCRun2) ;
	  aceOutf (ao, "\tSmoothed cumul of controls") ;
	}
      else
	{
	  aceOutf (ao, "\t%s", dictName (gx->runDict, run2)) ;
	  aceOutf (ao, "\t%s", dictName (gx->runDict, run1)) ;
	  aceOutf (ao, "\t%d\t%d", nCRun2, nCRun1) ;
	  aceOutf (ao, "\tSmoothed cumul of controls") ;
	}
      for (i = 0 ; i <= 200 ; i++)
	aceOutf (ao, "\t%.2f", cX[i]) ;
      aceOut (ao, "\n") ;
      
      /* export the smoothed cumul of the experiment */
      aceOutf (ao, "%s", dictName (gx->pairDict, compare->compare)) ;
      aceOutf (ao, "\t%s", dictName (gx->runDict, run1)) ;
      aceOutf (ao, "\t%s", dictName (gx->runDict, run2)) ;
      aceOutf (ao, "\t%d\t%d", nRun1, nRun2) ;
      aceOutf (ao, "\tSmoothed cumul of experiment") ;
      for (i = 0 ; i <= 200 ; i++)
	aceOutf (ao, "\t%.2f", cE[i]) ;
      aceOut (ao, "\n") ;
      
      /* export the FDR */
      aceOutf (ao, "%s", dictName (gx->pairDict, compare->compare)) ;
      aceOutf (ao, "\t%s", dictName (gx->runDict, run1)) ;
      aceOutf (ao, "\t%s", dictName (gx->runDict, run2)) ;
      aceOutf (ao, "\t%d\t%d", nRun1, nRun2) ;
      aceOutf (ao, "\tFDR") ;

      for (i = 0 ; i <= 200 ; i++)
	{
	  int v = cE[i] ;
	  double z, w = cX[i] ;
	  if (v == 0) v = 1 ;
	  z = w/v ;
	  if (z > 1) z = 1 ;
	  if (ao)
	    aceOutf (ao, "\t%.2g", z) ;
	  
	  if (! hasControls) z = -.01 ;
	  fdr[i] = z ;
	}
      aceOutf (ao, "\n") ;

      /* sex genes and other genuine baddies with very high randomized score 
       * are deduced from the FDR histo by shiftHisto 
       */
      /* export the TDR */

      if (ao)
	{
	  aceOutf (ao, "%s", dictName (gx->pairDict, compare->compare)) ;
	  aceOutf (ao, "\t%s", dictName (gx->runDict, run1)) ;
	  aceOutf (ao, "\t%s", dictName (gx->runDict, run2)) ;
	  aceOutf (ao, "\t%d\t%d", nRun1, nRun2) ;
	  aceOutf (ao, "\tTheta") ;
	}
      for (iBest = 200, i = 1 ; i <= 200 ; i++)
	{ 
	  int v = cE[i] ; 
	  int w = cX[i] ; 
	  int vw = v - w ; 

	  if (vw < 0) vw = 0 ;
	  if (vw <= 0) z = 0 ;
	  else if (w < .1) z = 10*vw ;
	  else z = vw/w ;
	  if (! hasControls) z = 1 ;
	  if (ao) 
	    aceOutf (ao, "\t%.1f", z) ;
	  if (z > best && fdr[i] < .7 && (i > 40 || fdr[i] < .1))
	    { best = z ; iBest = i ; if (fdr[i] < .001) break ; }
	}
      aceOut (ao, "\n") ;
      /* move up again untill the local noise reaches 20% of the local signal */
      for (i = iBest ; i>0 ; i--)
	if (
	    (fdr[i] > .045 && (fdr[i] > .4 || 100 * cX[i] >   20 * cE[i])) ||
	    100 * cX[i] >   30 * cE[i]
	    )
	  break ;
      iBest = i ;
      /* but if we already have over 10% FDR move up again and increase the number of genes until the FDR moves up by another 10% */
      if (fdr[iBest]> 0.2)
	{
	  for (i = iBest ; i>0 ; i--)
	    {
	      if (100 * cX[i] >  30 * cE[i])
		break ;
	      if (fdr[i] > fdr[iBest] + 0.10 || 2*fdr[i] > 3*fdr[iBest])
		break ;
	    }
	  i++ ;
	}
      /* look for another higher minimum */
      iBest = i < 200 ? i : 200 ;
      if (iBest > 80)
	{
	  float z = fdr[iBest] ;
	  for (i = iBest ; z < .5 && i>100 ; i--)
	    {
	      if (i < 150 && 100 * cX[i] >   60 * cE[i])
		break ;
	      if (i < 100 && 100 * cX[i] >   30 * cE[i])
		break ;
	      if (fdr[i] < z)
		{ iBest = i ; z = fdr[i] ; }
	    }
	  i = iBest ;
	}
      /* move right a little */
      iBest = i ;
      if (fdr[iBest] > .05)
	{
	  float z, z0 = fdr[iBest] ;
	  for (i = iBest ; i < 180 && i < iBest + 20 ; i++)
	    {
	      z = fdr[i] ;
	      if (z <= .05)
		{ iBest = i ; break ; }
	      if (z < z0)
		{ iBest = i ; z0 = z ; }
	    }
	  i = iBest ;
	}
      /* move right a little */
      if (fdr[i] > .05)
	{ 
	  if (i < 150) i += 5 ; 
	  else if (i < 160) i += 4 ; 
	  else if (i < 170) i += 3 ; 
	  else if (i < 180) i += 2 ;
	  else if (i < 195) i++ ;
	}
      if (! cE[iBest]) /* no gene */
	i = 200 ;
      
      if (! hasControls) i = 180 ;
      iBest = i ;
      
      if (compare->forceIbest)
	iBest = compare->forceIbest ;
      
      myFdr [isPm] = fdr[iBest] ;
      myFdrThreshold[isPm] = iBest ;

      if (isPm == 1)
	{
	  for (i = 0 ; i <= 200 ; i++)
	    compare->fdr1[i] = fdr[i] ;
	}
      else
	{
	  for (i = 0 ; i <= 200 ; i++)
	    compare->fdr2[i] = fdr[i] ;
	}
  
      if (isPm == 2 &&  myFdr [isPm - 1] > 0 &&
	  myFdr [isPm] > myFdr [isPm - 1] &&
	  myFdrThreshold[isPm] <  myFdrThreshold[isPm - 1]
	  )
	{
	  for (i = myFdrThreshold[isPm] ; i < myFdrThreshold[isPm - 1] ; i++)
	    if (fdr[i] <= myFdr [isPm - 1])
	      break ;
	  iBest = i ;
	  myFdr [isPm] = fdr[iBest] ;
	  myFdrThreshold[isPm] = iBest ;
	}
      
      if (ao)
	{
	  aceOutf (ao, "%s", dictName (gx->pairDict, compare->compare)) ;
	  aceOutf (ao, "\t%s", dictName (gx->runDict, run1)) ;
	  aceOutf (ao, "\t%s", dictName (gx->runDict, run2)) ;
	  aceOutf (ao, "\t%d\t%d", nRun1, nRun2) ;
	}
      if (cE[0] == 0)
	{ iBest = 0 ; fdr[0] = 0 ; } 
      myFdr[isPm] = fdr[iBest] ;
      myFdrThreshold [isPm] = iBest ;
      myDEG[isPm] = cE0[iBest] ;
      myDEG185[isPm] = iBest < 185 ? cE0[185] : cE0[iBest] ;

      if (0 && run2 == 15 && run1 == 16 && isPm == 2)
	{
	  for (i = 200 ; i >= 185 ; i--)
	    fprintf(stderr, "+/+/+/ cE[%d]=%.0f myDEg=%d iBest=%d\n", i, cE[i], myDEG[isPm], iBest) ;
	}
      if (1)
	{
	  float f =  100*fdr[iBest] ;
	  if (iBest == 200)
	    {
	      myDEG[isPm] = 0 ;
	      myDEG185[isPm] = 0 ;
	      f = 0 ;
	    }
	  if (ao)
	    aceOutf (ao, "\tSelected threshold\t%d\tDEG\t%d\t%%FDR\t%.3g\tDEG185\t%d\n\n", iBest, myDEG[isPm], f, myDEG185[isPm]) ;	  
	}
    }

  if (ao && isPm == 2)
    {
      /* export the smoothed experiments */
      aceOutf (ao, "%s", dictName (gx->pairDict, compare->compare)) ;
      aceOutf (ao, "\t%s", dictName (gx->runDict, run2)) ;
      aceOutf (ao, "\t%s", dictName (gx->runDict, run1)) ;
      aceOutf (ao, "\t%d\t%d", nRun2, nRun1) ;
      aceOutf (ao, "\tSmoothed experiment >") ;
      for (i = 0 ; i <= 200 ; i++)
	aceOutf (ao, "\t%.2f", aExp1[i]) ;
      aceOut (ao, "\n") ;

      aceOutf (ao, "%s", dictName (gx->pairDict, compare->compare)) ;
      aceOutf (ao, "\t%s", dictName (gx->runDict, run1)) ;
      aceOutf (ao, "\t%s", dictName (gx->runDict, run2)) ;
      aceOutf (ao, "\t%d\t%d", nRun1, nRun2) ;
      aceOutf (ao, "\tSmoothed experiment <") ;
      for (i = 0 ; i <= 200 ; i++)
	aceOutf (ao, "\t%.2f", aE[i]) ;
      aceOut (ao, "\n\n") ;
    }

  if (hasControls && ao && isPm == 2) /* export the smoothed controls */
    {      
      /* export the average of the controls */
      aceOutf (ao, "%s", dictName (gx->pairDict, compare->compare)) ;
      aceOutf (ao, "\t%s", dictName (gx->runDict, run2)) ;
      aceOutf (ao, "\t%s", dictName (gx->runDict, run1)) ;
      aceOutf (ao, "\t%d\t%d", nRun1, nRun2) ;
      aceOutf (ao, "\tAverage of the %d randomized resamplings", compare->nVirtualStrata) ;
      for (i = 0 ; i <= 200 ; i++)
	aceOutf (ao, "\t%.2f", aX[i]) ;
      aceOut (ao, "\n") ;

      for (ii = 0 ; ii < arrayMax (degHistos) ; ii++)
	{
	  dgh = arrp (degHistos, ii, DGH)  ;
	  if (dgh->iVV >= 2)
	    {
	      aceOutf (ao, "%s", dictName (gx->pairDict, compare->compare)) ;
	      aceOutf (ao, "\t%s", dictName (gx->runDict, run2)) ;
	      aceOutf (ao, "\t%s", dictName (gx->runDict, run1)) ;
	      aceOutf (ao, "\t%d\t%d", nCRun2, nCRun1) ;
	      
	      aceOutf (ao, "\tResampling %d", dgh->iVV - 1) ;
	      for (i = 0 ; i <= 200 ; i++)
		aceOutf (ao, "\t%.2f", dgh->histo[i]/GAUSS) ;
	      aceOut (ao, "\n") ;
	    }
	}
    }

 done:

 if (isPm)
    {	
      COMPARE *cmp = compare ;
      if (! cmp->h)
	cmp->h = ac_new_handle () ;
      cmp->nDEG[isPm] = myDEG[isPm] ;
      cmp->nDEG185[isPm] = myDEG185[isPm] ;
    }

  /* now that both directions have been exported we clean up */
  if (isPm == 2)
    arrayMax (gx->degHistos) = 0 ;
  if (isPm == 2)
    ac_free (ao) ;
  ac_free (h) ;
} /* gxDghExportHistos */

/*************************************************************************************/
/*************************************************************************************/
/* select genes above a given differential and export them 
 * in virtual case, run1 == run2 == 0, but vA and vB are given in extension
 */
static BOOL gxCompareOrControl (GX *gx, COMPARE *compare, int run1, int run2
				, KEYSET vA, KEYSET vB
				, int isVirtual)
{
  /* select the characteristic genes */

  /* select the characteristic genes */
  if (isVirtual)
    gxDoCompare (gx, compare, run1, run2, vA, vB, isVirtual, 0) ; /* virtual */
  if (! isVirtual)
    {
      gxDoCompare (gx, compare, run1, run2, vA, vB, 0, 1) ; /* pp */
      gxDoCompare (gx, compare, run2, run1, vB, vA, 0, 2) ; /* pm */
    }
  return TRUE ;
}  /* gxCompareOrControl */

/*************************************************************************************/
/* given 2 groups A and B , contruct virtual groups vA and vB 
 * such that vA, vB have the same cardinal as A and B
 *           the proportion of A and B in vA and vB are equal and equal to A/B
 * The goups are constructed using the by orbiting around a prime
 */

static BOOL gxMakeVirtualGroups (GX *gx, int run1, int run2, KEYSET vA, KEYSET vB, int iVV)
{
  RC *rc1 = arrayp (gx->runs, run1, RC) ;
  RC *rc2 = arrayp (gx->runs, run2, RC) ;
  KEYSET A = rc1->runs ;
  KEYSET B = rc2->runs ;
  KEYSET AorB = 0 ;
  int nA = A ? keySetMax (A) : 1 ;
  int nB = B ? keySetMax (B) : 1 ;
  int nAB = nA + nB ;
  int i, j, zz, n ;
  BitSet bb = bitSetCreate(nAB, 0) ;
  BOOL ok = FALSE ;
  BOOL debug = FALSE ;

  if (nA <= 1 || nB <= 1) /* otherwise we destroy the rc->aa table, i.e. the run counts when crating the virtual rc */
    return FALSE ;
  
  ok = TRUE ;
  if (1)
    {
      int x, y, z, t ; /* we want (x,y) in vA, (z,t) in vB
			* where (x,z) is from A and (y,t) is from B
			* and all ratio are equal x/y = z/t => x/z = y/t = (x+y)/(z+t)=nA/nB
			* and x+z <= nA , y+t <= nB , 
			* and all numbers as as high as we can
			*/
  
      if (0)
	{  /* in this method vA and vB are balanced but have the same size */
	  x = z = nA/2 ; y = t = nB/2 ;
	}
      else
	{ /* in this method we try to have vA vB in the proportions of A B, say 5 versus 20 */
	  int nnA = nA <= nB ? nA : nB ;   /* nnA is the smaller group */
	  int nnB = nA <= nB ? nB : nA ;
	  int tA ;
	  int k = 1 + nnB / nnA ;  /* always >= 2 */
	  int X = 0, Y, Z, T ;
	  float bestDelta = 10000000, delta, r, r0 ;
	  r0 = nnB ; r0 /= nnA ;  /* size ratio we want to optimize size vB/vA relative to B/A */
	  for (tA = 2 * (nnA/2) ; tA <= nnA + 1 || tA <= 1.2 * nnA ; tA++)  /* we may prefer a small vA group slightly larger than nnA */
	    {
	      for (x = 1 ; x <= 1 + tA/k  && x < tA ; x++)    /* vA = {x,y}   vB = {z,t} */
		{
		  int k0, k1, k2 ;
		  
		  y = tA - x ; /* complete group vA */
		  k1 = (nnA - x)/x ; /* maximal multiplicator */
		  k2 = (nnB - y)/y ; /* maximal multiplicator */
		  k0 = k1 < k2 ? k1 : k2 ;
		  if (k0 < 1) continue ;
		  z = k0 * x ; t = k0 * y ;  /* the 2 groups are balanced */
		  if (x + z <= nnA && y + t <= nnB) /* always TRUE */
		    {
		      r = z + t ; r /= (x + y) ;
		      delta = (r - r0) * ( r - r0) ;
		      if (delta <= bestDelta) /* optimize the size of the groups */
			{ 
			  X = x ; Y = y ; Z = z ; T = t ; bestDelta = delta ;
			}
		    }
		}
	    }
	  if (X == 0) return FALSE ;
	  if (nA <= nB)
	    { x = X ; y = Y ; z = Z ; t = T ; }
	  else
	    { x = T ; y = Z ; z = Y ; t = X ; } /* swith back the A/B groups */
	}

      if (x * y * z * t < 1) goto done ;
      /* we have the desired partition, now we pick randomly */
      /* u = (x+y) * (x+y - 1) * (z+t) * (z+t - 1) / 2 ;  total number of splits */
      if (1)
	{
	  for (zz = 0 ; zz < 2 ; zz++)
	    {
	      AorB = (zz == 0 ? A : B) ;
	      int n1 = (zz == 0 ? nA : nB) ;
	      int n2 = (zz == 0 ? x : y) ;
	      int delta ;
	      int primes[] = {11,23,17,13,29,23,31,37,41,0} ; 
	      /* find delate prime relative to n1 */
	      for (i = 0 ; primes[i] ; i++)
		{
		  delta = primes[i] ;
		  if (n1 - delta * (n1/delta) > 0)
		    break ;
		  delta = 1 ;
		}
	      bitSetMINUS (bb, bb) ;
	      for (i = 0, j = arrayMax (vA) ; i < n2 ; i++)
		{
		  int k1, k  =  ((randint () % n1) + n1) % n1 ;
		  for (k1 = 0 ; k1 < n1 ; k1++)
		    {
		      n = (k + k1 * delta) % n1 ;
		      if (! bit (bb, n))
			{
			  bitSet (bb, n) ;
			  keySet (vA, j++) = keySet (AorB, n) ;
			  if (debug) fprintf (stderr, "vA zz=%d n=%d k=%d %s\n",zz,n,k,dictName(gx->runDict,(keySet (AorB, n)))) ;
			  break ;
			}
		    }
		}
	      n2 = (zz == 0 ? z : t) ; 
	      for (i = 0, j = arrayMax (vB) ; i < n2 ; i++)
		{
		  int k1, k  =  ((randint () % n1) + n1) % n1 ;
		  for (k1 = 0 ; k1 < n1 ; k1++)
		    {
		      n = (k + k1 * delta) % n1 ;
		      if (! bit (bb, n))
			{
			  bitSet (bb, n) ;
			  keySet (vB, j++) = keySet (AorB, n) ;
			  if (debug) fprintf (stderr, "vB zz=%d n=%d k=%d %s\n",zz,n,k,dictName(gx->runDict,(keySet (AorB, n)))) ;
			  break ;
			}
		    }
		}
	    }
	}
    }
  keySetSort (vA) ;  keySetSort (vB) ; keySetCompress (vA) ; keySetCompress (vB) ;
  
  if (debug)
    {  /* we must verify that the vA vB groups have the correct size and proportions */
      KEYSET ksa, ksb ;

      ksa = keySetAND (A, vA) ;
      ksb = keySetAND (B, vB) ;

      fprintf (stderr, "gxMakeVirtualGroups %s:%d  %s:%d ---> vA:%d=%d+%d  vB:%d=%d+%d\t%s\n"
	       , dictName(gx->runDict, run1), nA
	       , dictName(gx->runDict, run2), nB
	       , keySetMax (vA), keySetMax(ksa),  keySetMax (vA) - keySetMax(ksa)
	       , keySetMax (vB),  keySetMax (vB) - keySetMax(ksb), keySetMax(ksb)
	       , timeShowNow ()
	       ) ;
      if (1) 
	{
	  fprintf (stderr, "------A :") ;
	  for (i=0 ; i < 10 && i < keySetMax(A) ; i++)
	    fprintf (stderr, " %s", dictName(gx->runDict, keySet (A,i))) ;
	  fprintf (stderr, "\n------B :") ;
	  for (i=0 ; i < 10 && i < keySetMax(B) ; i++)
	    fprintf (stderr, " %s",  dictName(gx->runDict, keySet (B,i))) ;
	  fprintf (stderr, "\n------vA :") ;
	  for (i=0 ; i < 10 && i < keySetMax(vA) ; i++)
	    fprintf (stderr, " %s", dictName(gx->runDict, keySet (vA,i))) ;
	  fprintf (stderr, "\n------vB :") ;
	  for (i=0 ; i < 10 && i < keySetMax(vB) ; i++)
	    fprintf (stderr, " %s",  dictName(gx->runDict, keySet (vB,i))) ;
	  fprintf (stderr, "\n") ;
	}
      keySetDestroy (ksa) ;
      keySetDestroy (ksb) ;
    }
 done:
  bitSetDestroy (bb) ;
  return ok ;
} /* gxMakeVirtualGroups */

/*************************************************************************************/
/* clean up */
static void gcCleanUp (GX *gx)
{
  int gene, geneMax ;
  GC *gc ;
  
  for (gene = 0, gc = arrp (gx->genes, 0, GC), geneMax = arrayMax (gx->genes) ; gene < geneMax ; gc++, gene++)
    {
      gc->score = gc->score1 = gc->score2 = 0 ; 
      gc->za1 = gc->za2 = gc->zb1 = gc->zb2 = 0 ;
    }
} /* gcCleanUp */

/*************************************************************************************/
/* select genes above a given differential and export them 
 * optionally count the FDR on a set of virtual randomizedly stratified control classes
 */
static BOOL gxCompareAndControl (GX *gx, COMPARE *compare)
{
  AC_HANDLE h = ac_new_handle () ;
  int nVV = NVIRTUALSTRATA ; /* desired number of virtual strata */ 
  BOOL ok = FALSE ;
  int run1 = compare->vRun1 ;
  int run2 = compare->vRun2 ;
  RC *rc1 = arrayp (gx->runs, run1, RC) ;
  RC *rc2 = arrayp (gx->runs, run2, RC) ;
  KEYSET vA, vB ; /* virtual groups */
  
  vA = keySetHandleCreate (h) ;
  vB = keySetHandleCreate (h) ;

  compare->nVirtualStrata = 0 ;
  arrayDestroy (gx->baddyBatchies) ;
  gxDghInitializeBaddyBatchies (gx, h) ;

  if (! rc1->addCounts && ! rc2->addCounts  /* an addCount group is atomic like a run, it makes no sense to shuffle its runs */
      )
    {
      int iVV ;
      
      for (iVV = 0 ; iVV < nVV ; iVV++)
	{
	  vA = keySetReCreate (vA) ;
	  vB = keySetReCreate (vB) ;
	  if (gxMakeVirtualGroups (gx, run1, run2, vA, vB, iVV))
	    {
	      compare->nVirtualStrata++ ;
	      gcCleanUp (gx) ;
	      gxCompareOrControl (gx, compare, run1, run2, vA, vB, iVV+1) ;
	    }	
	  else
	    break ;
	}
    }

  /* finaly compare the actual run
   * by now the FDR histograms are available
   * rc1, rc2 may have been reallocated via the creation of virtual runs
   */

  if (rc1->runs)
    vA = rc1->runs ;
  else
    { vA = keySetReCreate (vA) ; keySet (vA, 0) = run1 ; }
  if (rc2->runs)
    vB = rc2->runs ;
  else
    { vB = keySetReCreate (vB) ; keySet (vB, 0) = run2 ; }

  ok = gxCompareOrControl (gx, compare, run1, run2, vA, vB, 0) ;

  ac_free (h) ;
  return ok ;
} /* gxCompareAndControl */

/*************************************************************************************/
/*************************************************************************************/
/* evaluate the score of a single SNP and export it */
static BOOL gxSnpEvaluate (GX *gx, int gene, const char *snpName)
{
  int ii, scoreBin[21], nHomo = 0, nHetero = 0, nObs = 0 ;
  double bestScore = 0 ;
  DC *dc1, *dc2 ;
  static int firstPass = 0 ;

  memset(scoreBin, 0, sizeof(scoreBin)) ;

  if (gx->snpEval)
    {
      int compare ;
      COMPARE *up ;

      for (compare = 0 ; compare < arrayMax (gx->compares) ; compare++)
	{
	  up = arrayp (gx->compares, compare, COMPARE) ;
	  if (up->compareSNP && up->runs)
	    {
	      int ii, jj, iMax = arrayMax (up->runs) ;
	      for (ii = 0 ; ii < iMax - 1 ; ii++)
		{
		  int run1 = keySet (up->runs, ii) ;
		  RC *rc1 = arrayp (gx->runs, run1, RC) ;
		  rc1->cumulDone = FALSE ;
		  if (!rc1->aa || ! arrayMax (rc1->aa))
		    {
		      if (1) fprintf (stderr, "# No measure for compareSNP %s run:%s\n"
				      , dictName(gx->compareDict, up->compare)
				      , dictName(gx->runDict, run1)
				      ) ;
		      continue ;
		    }
		  for (jj = ii+1 ; jj < iMax ; jj++)
		    {
		      int run2 = keySet (up->runs, jj) ;
		      RC *rc2 = arrayp (gx->runs, run2, RC) ; 
		      rc2->cumulDone = FALSE ;
		      if (!rc2->aa || ! arrayMax (rc2->aa))
			continue ;
		      if (ii < jj)
			{
			  double s1, s2, c1, c2, g1, g2, m1, m2 ;
			  
			  if (! firstPass)
			    {
			      firstPass++ ;
			      printf ("# SNV\tRun\tTitle\tGenotype\t%% variant\tCoverage\tVariant\tReference\tRun\tTitle\tGenotype\t%% variant\tCoverage\tVariant\tReference\tchi2 (10.83->1/1000, 15->1/10,000, 20->1/100,000, 24->1/1,000,000\n") ;
			    }
			  
			  dc1 = arrp (rc1->aa, 1, DC) ;  
			  c1 = dc1->tags ;
			  if (c1 <10) 
			    continue ;
			  m1 = dc1->seqs ;
			  s1 = 100.0*m1/(c1+.0001) ;
			  g1 = dc1->kb ;
			  if (g1 == 0 && m1 > 0) g1 = .3 ;
			  if (g1 == 1 && s1 < 40) g1 = .3 ;
			  if (g1 == 1 && s1 > 65) g1 = 1.3 ;
			  if (g1 == 2 && c1 - m1 > 0) g1 = 1.3 ;
			  

			  dc2 = arrp (rc2->aa, 1, DC) ;  
			  c2 = dc2->tags ;
			  if (c2 <10) 
			    continue ;
			  m2 = dc2->seqs ;
			  s2 = 100.0*m2/(c2+.0001) ;
			  g2 = dc2->kb ;
			  if (g2 == 0 && m2 > 0) g2 = .3 ;
			  if (g2 == 1 && s2 < 40) g2 = .3 ;
			  if (g2 == 1 && s2 > 65) g2 = 1.3 ;
			  if (g2 == 2 && c2 - m2 > 0) g2 = 1.3 ;
	
			  if (chi2 (dc1->seqs, dc2->seqs, dc1->tags - dc1->seqs, dc2->tags - dc2->seqs, 0)  >= 0x4 /* risk 5/100, 0x4=1/100, 0x8 = 1/1000 */
			      /*	  ||
			       (g1 == 0 && g2 >=1) || 
			       (g2 == 0 && g1 >=1) || 
			       (g1 == 2 && g2 <=1) || 
			       (g2 == 2 && g1 <=1)  || 
			       (m1 == 0 && c1 >= 10 && g2 == .5) || 
			       (m2 == 0 && c2 >= 10 && g1 == .5) || 
			       (m1 == 2 && c1 >= 10 && g2 == 1.5) || 
			       (m2 == 2 && c2 >= 10 && g1 == 1.5)
			       */
			       )  
			    {
			      float cz = 0 ;
			      int im1=m1, iw1=c1-m1, im2=m2,iw2=c2-m2 ;
			      chi2 (im1, iw1, im2, iw2, &cz) ;

			      if (cz > 10.83) /*  10.83 ~ 1/1000  // 15: proba < 1/10,000 */
				printf ("%s\t%s\t%s\t%.1f\t%.1f\t%.0f\t%.0f\t%.0f\t%s\t%s\t%.1f\t%.1f\t%.0f\t%.0f\t%.0f\t%.3g\n"
					, snpName
					, dictName (gx->runDict, run1)
					, rc1->title ? stackText (gx->info, rc1->title) : "-" 
					, g1, s1, c1, m1, c1 - m1
					, dictName (gx->runDict, run2)
					, rc2->title ? stackText (gx->info, rc2->title) : "-" 
					, g2, s2, c2, m2, c2 - m2 
					, cz
					) ;
			    }
			}
		    }
		}
	    }
	}
      return TRUE ;
    }

  if (gx->snpEval && ((bestScore > 10 && nHomo + nHetero > 5) || nHomo > 5))
    {
      printf ("%s\t%.1f\t%d\t%d\t%d\t#", snpName, bestScore, nObs - nHetero - nHomo, nHetero, nHomo) ;
      for (ii = 0 ; ii < 21 ; ii++)
	printf ("\t%d", scoreBin[ii]) ;
      printf ("\n") ;
    }

  return bestScore > 50 ? TRUE : FALSE ;
} /* gxSnpEvaluate */

/*************************************************************************************/
/*************************************************************************************/
/* symmetrize */
static int keySetMerge (KEYSET ks1, KEYSET ks2, Associator assPlus, Associator assMinus)
{
  int i, j, gene ;

  assReCreate (assPlus) ;
  assReCreate (assMinus) ;
  for (i = j = 0 ; i < keySetMax (ks1) ; i++, j++)
   {
     gene = keySet (ks1, i) ;
     assInsert (assPlus, assVoid(gene), assVoid(i+1)) ;
   }
  for (i = 0 ; i < keySetMax (ks2) ; i++, j++)
    {
      gene = keySet (ks2, i) ;
      assInsert (assMinus, assVoid(gene), assVoid(i+1)) ;
    }

  return j ;
} /* keySetMerge */

/*************************************************************************************/
/* select genes above a given differential and export them */
static int gxCompare_to (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ; 
  int iCompare, nn = 0, run1, run2 ;
  RC *rc1, *rc2 ;
  Array rws = 0 ;
  KEYSET genePlus = 0, geneMinus = 0 ;
  Array sp, sm ;

  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, "// gxCompare_to start :  %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
  }

  sp = arrayHandleCreate (100, SCAL, h) ;
  sm = arrayHandleCreate (100, SCAL, h) ;

  gxExportSignatures (gx, 0, 0, 0, 0, -1) ; /* initialise the signature file */

  for (iCompare = 1 ; gx->compares && iCompare < arrayMax(gx->compares) ; iCompare++)
    {
      COMPARE *cmp, *compare = arrayp (gx->compares, iCompare, COMPARE) ;
      if (! compare->runs || keySetMax (compare->runs) != 2)
	continue ;
      if (gx->isINTRON && ! compare->compareINTRON)
	continue ;
      if (gx->isSNP && ! compare->compareSNP)
	continue ;
      if (compare->private)
	continue ;

      /* TOTOTOTO */
       run1 = keySet (compare->runs, 0) ;
       run2 = keySet (compare->runs, 1) ;
	  
       cmp = gxPairCmp (gx, run1, run2) ;
       if (run1 == cmp->vRun1)
	 {
	   genePlus = cmp->genePlus[1] ;
	   geneMinus = cmp->genePlus[2] ;
	 }
       else
	 {
	   genePlus = cmp->genePlus[2] ;
	   geneMinus = cmp->genePlus[1] ;
	 }

       if (! cmp)
	 continue ;
       if (! cmp->done)
	 continue ;
       rc1 = arrayp (gx->runs, run1, RC) ;
       if (! rc1->aa || rc1->groupLevel < 0) continue ;
       rc2 = arrayp (gx->runs, run2, RC) ;
       if (! rc2->aa || rc2->groupLevel < 0) continue ;
       
       if (1)
	 {
	   /* genePlus and geneMinus are sorted by importance */
	   /* order the runs */
	   sp = arrayReCreate (sp, keySetMax (genePlus), SCAL) ;
	   sm = arrayReCreate (sm, keySetMax (geneMinus), SCAL) ;

	   /* all histos should be registered probably in sp, sm, passed to cosines in case we use digital == 3, and passed to export comparison, so we do not recompute them 
	      double hh0s[256], hh1s[256] ;
	      gxMakeComparativeHistos (gx, gene, run1, run2, rc1->runs, rc2->runs, hh0s, hh1s
	      , 0, &nRun0, &nRun1, 0, 0, 0
	      , 0
	      , &Gini,  &pValue
	      , FALSE
	      ) ;
	  */
	       
	   rws = gxAllRunsCosineUsingGenesRun1Run2 (gx, run1, run2, genePlus, geneMinus, sp, sm, h) ; /* rw: {run,alpha} */
	   if (! rws || ! arrayMax (rws)) 
	     continue ;

	   if (0) /* there is another call to gxTRueFalse elsewhere if compare->predicts */
	     gxTrueFalse_AUC_MCC (gx, 0, rws, run1, run2, 0, 0, 0) ; /* computes the AUC, Gini, pValue and in each run sets rc->alpha, beta1, beta2 */
	   gxExportSignatures (gx, run1, run2, 0, 0, 0) ;      /* exports in a single file for all compare, the functions alpha, beta1, beta2  */
	   if (1)  /* beta coeff 2014_11_13 : gives the alpha classification of all runs showing if the signature works on training, if applicable also scores the predicts->compare */
	     {
	       gxExportComparison (gx, compare, rws, genePlus, geneMinus, run1, run2, FALSE, FALSE) ; /* beta: alpha, beta1,beta2 plots, contrib of each gene to alpha */
	       gxExportComparison (gx, compare, rws, genePlus, geneMinus, run1, run2, TRUE, FALSE) ;  /* histo: histos of all genes used in the signatures */
	       gxExportComparison (gx, compare, rws, genePlus, geneMinus, run1, run2, TRUE, TRUE) ;  /* histo: histos of all pathways used in the signatures */
	     }
	 }
    }
  
  gxExportSignatures (gx, 0, 0, 0, 0, -2) ; /* close the signature file */

  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, "// gxCompare_to done :  %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
  }
  ac_free (h) ;
  return nn ;  
} /* gxCompare_to */
 
/*************************************************************************************/
/*************************************************************************************/

static void gxExportTableHeaderLegend (GX *gx, ACEOUT ao, const char* title, int type)
{
  aceOutf (ao, "\n#%s", title) ;
  if(type == 51) ;
  else if (type == 41)
    {
      if (gx->method) aceOut (ao, "\t#Method") ;
      if (gx->hasRunId) aceOut (ao, "\t#RunId") ;
      if (gx->hasRunSample) aceOut (ao, "\t#Sample") ;
      if (gx->hasRunTitle) aceOut (ao, "\t#Title") ;	
      if (type == 22 || type == 32) aceOut (ao, "\t#Run") ;
      aceOutf (ao, "\t#%s", title) ;
    }
  else if (type == 52) ;
  else 
    {
      if (gx->method) aceOut (ao, "\t#Method") ;
      if (gx->hasGeneId) aceOut (ao, "\t#NCBI GeneId") ;
      if (0 && gx->hasGeneNm) aceOut (ao, "\t#RefSeq transcript Id") ;
      if (gx->hasGeneLength) aceOut (ao, "\t#Length") ;
      if (gx->hasKnownDonorAcceptor) aceOut (ao, "\t#Known/New donor acceptor") ;
      if (gx->hasFromGene) aceOut (ao, "\t#From gene") ;
      if (gx->hasIntronType) aceOut (ao, "\t#Type") ;
      if (gx->isINTRON) aceOut (ao, "\t#In AceView") ;
      if (gx->hasFromTranscript) aceOut (ao, "\t#From transcript") ;
      if (gx->hasGeneAffy) aceOut (ao, "\t#Microarray") ;
      if (gx->hasCapture) aceOut (ao, "\t#Captured by3") ;
      if (gx->hasCaptureTouch) aceOut (ao, "\t#Touched by") ;
      if (gx->hasGeneChrom) aceOut (ao, "\t#Chromosome\t#Strand\t#from base\t#to base") ;
      if (gx->hasGeneType) aceOut (ao, "\t#Type") ;
      if (gx->hasGeneTitle) aceOut (ao, "\t#Title") ;
      aceOut (ao, "\t#Measured object") ;
      if (gx->htmlSpecies) aceOut (ao, "\t#Link") ;
      aceOutf (ao, "\t#%s", title) ;
    }
} /* gxExportTableHeaderLegend  */

/*************************************************************************************/

static int bigGeneOrder (const void *a, const void *b)
{
  const BG *up = (const BG *)a, *vp = (const BG *)b ;
  float z = up->kb - vp->kb ;

  if (z > 0) return -1 ;
  if (z < 0) return 1 ;
  return up->gene - vp->gene ;
}

static void gxBigGenesDo (ACEOUT ao, GX *gx, RC *rc, const char *target, BOOL isAce, BOOL isCapture)
{
  KEYSET bigG = isCapture ?  rc->capturedBigGenes : rc->bigGenes ;

  if (target)
    aceOutf (ao, "-D High_genes %s%s\n"
	     , isCapture ? "Captured_" : ""
	     , target
	     ) ;
  
  if (bigG && keySetMax (bigG))
    {
      AC_HANDLE h = ac_new_handle () ;
      int i ;
      BG *bg ;
  
      arraySort (bigG, bigGeneOrder) ;
      arrayCompress (bigG) ;
      for (i = 0, bg = arrp (bigG, i, BG) ; i < arrayMax (bigG) ; bg++, i++)
	{
	  if (! bg->gene) continue ;
	  if (isAce)
	    aceOutf (ao, "High_genes %s%s \t%s %.2f Mb\n"
		     , isCapture ? "Captured_" : ""
		     , target
		     , ac_protect (noX__ (dictName (gx->geneDict, bg->gene), !gx->isINTRON), h)
		     , bg->kb/1000 
		     ) ;
	  else
	    aceOutf (ao, "%s%s (%.2f Mb)"
		     , (i ? ", " : "")
		     , noX__ (dictName (gx->geneDict, bg->gene), !gx->isINTRON)
		     ,  bg->kb/1000 
		     ) ;
	}
      ac_free (h) ;
    }
}  /* gxBigGenesDo */

/********/

static void gxBigGenes (ACEOUT ao, GX *gx, RC *rc, const char *target, BOOL isAce)
{
  if (isAce)
    aceOutf (ao, "-D High_genes %s\n", target) ;
  gxBigGenesDo (ao, gx, rc, target, isAce, FALSE) ;
  if (isAce && gx->hasCapture && rc->capturedBigGenes)
    gxBigGenesDo (ao, gx, rc, target, isAce, TRUE) ;
} /* gxBigGenes */

/*************************************************************************************/

static void gxExtendedRunIds (ACEOUT ao, GX *gx, RC *rc, char* sep, int pass)
{
  if (rc->runid)
    aceOutf (ao, "%s%s", sep,  stackText (gx->info, rc->runid)) ;
  else if (rc->runs && rc->groupLevel < 1 && pass < 2)
    {
      int i, j ;
      for (i = j = 0 ; i < keySetMax (rc->runs) ; i++)
	{
	  RC *rc2 = arrp (gx->runs, keySet (rc->runs, i), RC) ;
	  if (rc2->aa && ! rc2->isAny)
	    {
	      gxExtendedRunIds (ao, gx, rc2, sep, pass + 1) ;
	      sep = ", " ;
	    }
	}
    }
}  /* gxExtendedRunIds */

/*************************************************************************************/
/* export a table of index, tags, kb, variabce ... 
 * type == 3: just the groups
 */
static void gxExportTableHeader  (GX *gx, ACEOUT ao, int type)
{
  int run ;
  RC *rc ;
  int runMax = dictMax(gx->runDict) + 1, pass ;
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao2 = ! gx->isINTRON && ! gx->keepIndex && type == 41 ? aceOutCreate (gx->outFileName, ".abaque.txt", FALSE, h) : 0 ;
  int  doNotExportMetaData = 0 ;
  int iDoubleTitle, doubleTitle = (type == 5 ? 2 : 1) ;
  const char *GM = gx->isTranscript ? "Transcripts" : (gx->isMA ? "Probes" : "Genes") ;

  if (type >= 1000)
    { type = 51 ; doNotExportMetaData = 1 ; }
  else
    aceOutDate (ao, "###", gx->title) ;

  if (runMax > arrayMax(gx->runs))
    runMax = arrayMax(gx->runs) ;

  if (ao2) aceOutf (ao2, "## %s\t%s\tfile=%s", timeShowNow(), gx->title ? gx->title : "", aceOutFileName (ao2)) ;
  
  if (gx->isINTRON)
    {
      if (gx->unique) aceOutf (ao, "\n# Exon junctions and their unique support per group and run.") ;
      else aceOutf (ao, "\n# Exon junctions and their support per group and run.") ;
      if (gx->unique) aceOutf (ao, "\n# Only uniquely mapped reads contribute to this table ") ;
      aceOutf (ao, "\n# On the indicated chromosome, a1 and a2 are the coordinates of the first and last base of the intron ") ;
      aceOutf (ao, "in the reference genome %s", gx->referenceGenome) ;
      aceOutf (ao, "\n# For example in a gt_ag intron, a1 is the coordinate of the g of gt_, a2 of the g of _ag.") ;
      aceOutf (ao, "\n# If a1 < a2, the intron is on the + strand, if a1 > a2 on the - strand.") ;
      aceOutf (ao, "\n# All coordinates are 1-based, i.e. the first base of the chromosome is numbered 1.") ;
      aceOutf (ao, "\n# For each group and each run, we count the number of reads with best match ") ;
      aceOutf (ao, "across the exon-junction") ;
      if (1)
	{
	  aceOutf (ao, ", covering at least 8 nucleotides of each bordering exon.") ;
	  aceOutf (ao, "\n# Considering that MAGIC clips the RNA alignments ") ;
	  aceOutf (ao, "right before mismatches closer than 8 nucleotides from the edge of the alignment, all reads supporting an intron ") ;
	  aceOutf (ao, "are either locally exact or extend further than 8 nucleotides inside the bordering exons.") ;
	}
      aceOutf (ao, "\n# The RefSeq and EST columns show the number of Genbank/dbEST supporting accessions.") ; 
      aceOutf (ao, "\n# The RNA_seq support column shows the number of reads supporting the intron in all RNA-seq analyzed so far in MAGIC.") ; 

      if (0) aceOutf (ao, "\n# \n# Chromosome\ta1\ta2") ;
      if (0) aceOutf (ao, "\tStrand\tType\tIntron Length\tGene\tRefSeq\tEST\tDeep RNA-seq\tAll runs in project") ;
    }

  for (pass = doNotExportMetaData  ; pass < 2 ; pass++)
    {
      gxExportTableHeaderLegend (gx, ao, "Run", type) ;
      if (ao2) gxExportTableHeaderLegend (gx, ao2, "Run", type) ;
      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	{
	  if (! rc->run) rc->run = run ;
	  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
	    {
	      aceOutf (ao, "\t%s", dictName (gx->runDict, rc->run)) ;
	      if (doubleTitle == 2)
		aceOutf (ao, " %s", iDoubleTitle == 0 ? " index" : " Coeffficent of variation Variance") ;
	      if (ao2) aceOutf (ao2, "\t%s", dictName (gx->runDict, rc->run)) ;
	    }
	}

      if (doNotExportMetaData)
	continue ;
      
      if (gx->hasRunId)
	{
	  gxExportTableHeaderLegend (gx, ao, "Detailed RunIds", type) ;
	  if (ao2) gxExportTableHeaderLegend (gx, ao2, "Detailed RunIds", type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{   
		  aceOut (ao, "\t") ;
		  gxExtendedRunIds (ao, gx, rc, "", 0) ;
		  
		  if (ao2)
		    {
		      aceOut (ao2, "\t") ;
		      gxExtendedRunIds (ao2, gx, rc, "", 0) ;
		    }
		}
	    }
	}
      if (gx->hasRunSample)
	{
	  gxExportTableHeaderLegend (gx, ao, "Sample", type) ;
	  if (ao2) gxExportTableHeaderLegend (gx, ao2, "Sample", type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  if (rc->sample && ! rc->runs)
		    aceOutf (ao, "\t%s", stackText (gx->info, rc->sample)) ;
		  else
		    aceOutf (ao, "\t%s", dictName (gx->runDict, rc->run)) ;
		  if (ao2)
		    {
		      if (rc->sample && ! rc->runs)
			aceOutf (ao2, "\t%s", stackText (gx->info, rc->sample)) ;
		      else
			aceOutf (ao2, "\t%s", dictName (gx->runDict, rc->run)) ;
		    }
		}
	    }
	}
      
      if (gx->hasRunTitle)
	{
	  gxExportTableHeaderLegend (gx, ao, "Title", type) ;
	  if (ao2) gxExportTableHeaderLegend (gx, ao2, "Title", type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  if (rc->title)
		    aceOutf (ao, "\t%s", stackText (gx->info, rc->title)) ;
		  else
		    aceOutf (ao, "\t%s", dictName (gx->runDict, rc->run)) ;
		  if (ao2)
		    {
		      if (rc->title)
			aceOutf (ao2, "\t%s", stackText (gx->info, rc->title)) ;
		      else
			aceOutf (ao2, "\t%s", dictName (gx->runDict, rc->run)) ;
		    }
		}
	    }
	}
      
      if (gx->hasRunSortingTitle)
	{
	  gxExportTableHeaderLegend (gx, ao, "Sorting_title", type) ;
	  if (ao2) gxExportTableHeaderLegend (gx, ao2, "Sorting_title", type) ;
	  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
	    {
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {
		      if (rc->sortingTitle)
			aceOutf (ao, "\t%s", stackText (gx->info, rc->sortingTitle)) ;
		      else
			aceOut (ao, "\t") ;
		      if (ao2)
			{
			  if (rc->sortingTitle)
			    aceOutf (ao2, "\t%s", stackText (gx->info, rc->sortingTitle)) ;
			  else
			aceOut (ao2, "\t") ;
			}
		    }
		}
	    }
	}
      if (gx->hasRunSortingTitle2)
	{
	  gxExportTableHeaderLegend (gx, ao, "Sorting_title_2", type) ;
	  if (ao2) gxExportTableHeaderLegend (gx, ao2, "Sorting_title_2", type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  if (rc->sortingTitle2)
		    aceOutf (ao, "\t%s", stackText (gx->info, rc->sortingTitle2)) ;
		  else
		    aceOut (ao, "\t") ;
		  if (ao2)
		    {
		      if (rc->sortingTitle2)
			aceOutf (ao2, "\t%s", stackText (gx->info, rc->sortingTitle2)) ;
		      else
			aceOut (ao2, "\t") ;
		    }
		}
	    }
	}
      if (gx->hasRunOtherTitle)
	{
	  gxExportTableHeaderLegend (gx, ao, "Also known as", type) ;
	  if (ao2) gxExportTableHeaderLegend (gx, ao2, "Also known as", type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  if (rc->otherTitle)
		    aceOutf (ao, "\t%s", stackText (gx->info, rc->otherTitle)) ;
		  else
		    aceOut (ao, "\t") ;
		  if (ao2)
		    {
		      if (rc->otherTitle)
			aceOutf (ao2, "\t%s", stackText (gx->info, rc->otherTitle)) ;
		      else
			aceOut (ao2, "\t") ;
		    }
		}
	    }
	}
      
      if (pass == 1) /* STOP repeating the header */
	{
	  aceOut (ao, "\n") ;
	  if (ao2) aceOut (ao2, "\n") ;

	  continue ;
	}
      
      if (gx->hasMachine)
	{
	  gxExportTableHeaderLegend (gx, ao, "Machine", type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  if (rc->machine)
		    aceOutf (ao, "\t%s", stackText (gx->info, rc->machine)) ;
		  else
		    aceOut (ao, "\t") ;
		}
	    }
	}
      
      if (1 || ! gx->isINTRON)
	{
	  gxExportTableHeaderLegend (gx, ao, "Average zero count (NE) index", type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  aceOutf (ao, "\t%.2f", rc->zeroIndex - FIX) ;
		}
	    }
	}
      
      if (1 || ! gx->isINTRON)
	{
	  int iab, nab = 0 ;
	  for (iab = 0 ; iab <= 30 ; iab++)
	    {
	      if (! iab) nab = 0 ;
	      else if (iab == 1) nab = 1 ;
	      else nab <<= 1 ;
	      if (nab < 32)
		{
		  if (0) gxExportTableHeaderLegend (gx, ao, messprintf("%.2f fold", nab/32.0), type) ;
		  if (ao2) gxExportTableHeaderLegend (gx, ao2, messprintf("%.2f fold", nab/32.0), type) ;
		}
	      else
		{
		  if (0) gxExportTableHeaderLegend (gx, ao, messprintf("%d fold", nab/32), type) ;
		  if (ao2) gxExportTableHeaderLegend (gx, ao2, messprintf("%d fold", nab/32), type) ;
		}
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {
		      if (0) aceOutf (ao, "\t%.2f", rc->abaque[iab]) ;
		      if (ao2) aceOutf (ao2, "\t%.2f", rc->abaque[iab]) ;
		    }
		}
	    }
	}
      
      if (1 || ! gx->isINTRON)
	{
	  gxExportTableHeaderLegend (gx, ao, "Average Non Accurate (NA) measurable index", type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  aceOutf (ao, "\t%.2f", rc->NA_crossover_index - FIX) ;
		}
	    }
	}
      
      if (gx->isINTRON)
	{
	  gxExportTableHeaderLegend (gx, ao, "Exon junction with support", type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  aceOutf (ao, "\t%d", rc->Genes_with_index) ;
		}
	    }
	}
      if (1)
	{
	  if (1)
	    {
	      gxExportTableHeaderLegend (gx, ao, hprintf (h, "High %s " , GM), type) ;
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {
		      aceOut (ao, "\t") ;
		      if (rc->bigGenes) 
			gxBigGenes (ao, gx, rc, 0, FALSE) ;
		      
		      if (ao2)
			{
			  aceOut (ao2, "\t") ;
			  if (rc->bigGenes) 
			    gxBigGenes (ao, gx, rc, 0, FALSE) ;
			}
		    }
		}
	    }
	  
	  if (! gx->keepIndex)
	    {
	      gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s touched, including Non  Accurate (NA) values" , GM), type) ;
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {
		      aceOutf (ao, "\t%d", rc->Genes_with_index) ;
		    }
		}
	      gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s significantly expressed (for a group: in at least one run)" , GM), type) ;
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {
		      aceOutf (ao, "\t%d", rc->runs && ! rc->addCounts ? rc->Genes_with_one_reliable_index : rc->Genes_with_reliable_index) ;
		    }
		}
	      
	      gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s significantly expressed (for a group: in the majority of the runs)" , GM), type) ;
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {
		      aceOutf (ao, "\t%d", rc->Genes_with_reliable_index) ;
		    }
		}
	    }
	  gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s with index >= 8" , GM), type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  aceOutf (ao, "\t%d", rc->Genes_with_index_over_8) ;
		}
	    }

	  gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s with index >= 10" , GM), type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  aceOutf (ao, "\t%d", rc->Genes_with_index_over_10) ;
		}
	    }

	  gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s with index >= 12" , GM), type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  aceOutf (ao, "\t%d", rc->Genes_with_index_over_12) ;
		}
	    }

	  gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s with index >= 15" , GM), type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  aceOutf (ao, "\t%d", rc->Genes_with_index_over_15) ;
		}
	    }

	  gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s with index >= 18" , GM), type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  aceOutf (ao, "\t%d", rc->Genes_with_index_over_18) ;
		}
	    }

	  gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s with index >= 21" , GM), type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  aceOutf (ao, "\t%d", rc->Genes_with_index_over_21) ;
		}
	    }

	  gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s with index >= 24" , GM), type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  aceOutf (ao, "\t%d", rc->Genes_with_index_over_24) ;
		}
	    }

	  gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s with index >= 27" , GM), type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  aceOutf (ao, "\t%d", rc->Genes_with_index_over_27) ;
		}
	    }

	  gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s with index >= 30" , GM), type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  aceOutf (ao, "\t%d", rc->Genes_with_index_over_30) ;
		}
	    }

	  gxExportTableHeaderLegend (gx, ao, "Accessible transcript length, reflectig 3\' bias", type) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
	      for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		{
		  if (rc->accessibleLength8kb > 0)
		    aceOutf (ao, "\t%d", rc->accessibleLength8kb) ;
		  else if (rc->accessibleLength5kb > 0)
		    aceOutf (ao, "\t%d", rc->accessibleLength5kb) ;
		  else
		    aceOut (ao, "\t") ;
		}
	    }

	  if (! gx->keepIndex)
	    {
	      char *MB = (gx->isINTRON ? "Junction support" : "MB") ;
	      gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s aligned in this run, including non-unique" , "reads"), type) ;
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {
		      aceOutf (ao, "\t%.0f", rc->tags ) ;
		    }
		}

	      gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s aligned in these genes, including non-unique" , MB), type) ;
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {	
		      aceOutf (ao, "\t%.1f", rc->targetKb/1000 ) ;
		    }
		}

	      gxExportTableHeaderLegend (gx, ao, hprintf (h, "denominator: %s %sin genes with geneId minus high genes" , MB, gx->unique ? "uniquely aligned " : ""), type) ;
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {
		      aceOutf (ao, "\t%.1f", rc->prunedTargetKb/1000 ) ;
		    }
		}

	      gxExportTableHeaderLegend (gx, ao, hprintf (h, "%s aligned in intergenic, including non-unique" , MB), type) ;
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {
		      aceOutf (ao, "\t%.1f", rc->intergenicKb/1000 ) ;
		    }
		}

	      gxExportTableHeaderLegend (gx, ao, hprintf (h, "%% %s aligned in all genes, including high genes" , MB), type) ;
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {
		      if (rc->anyTargetKb > 0) 
			aceOutf (ao, "\t%.2f", 100.0 * rc->targetKb/rc->anyTargetKb ) ;
		      else
			aceOut (ao, "\tNA") ;
		    }
		}

	      gxExportTableHeaderLegend (gx, ao, hprintf (h, "%% %s used in denominator: genes with geneId minus high genes" , MB), type) ;
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {
		      if (rc->anyTargetKb > 0) 
			aceOutf (ao, "\t%.2f", 100.0 * rc->prunedTargetKb/rc->anyTargetKb ) ;
		      else
			aceOut (ao, "\tNA") ;
		    }
		}

	      gxExportTableHeaderLegend (gx, ao, hprintf (h, "%% %s aligned in the high genes" , MB), type) ;
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {
		      if (rc->anyTargetKb > 0) 
			aceOutf (ao, "\t%.2f", 100.0 * rc->bigGenesKb/rc->anyTargetKb ) ;
		      else
			aceOut (ao, "\tNA") ;
		    }
		}

	      gxExportTableHeaderLegend (gx, ao, hprintf (h, "%% Intergenic" , GM), type) ;
	      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
		{
		  if (rc->private || ! rc->aa || (type == 3 && run > 1 && ! rc->runs) || (type == 5 && ! rc->selectedVariance)) continue ;
		  for (iDoubleTitle = 0 ; iDoubleTitle < doubleTitle ; iDoubleTitle++)
		    {
		      if (rc->anyTargetKb > 0) 
			aceOutf (ao, "\t%.2f", 100.0 * rc->intergenicKb/rc->anyTargetKb ) ;
		      else
			aceOut (ao, "\tNA") ;
		    }
		}
	    }
	}
      if (0) aceOutf (ao, "\tMaximum observed in the top run\tPercentage of the total in the top run\tNumber of runs collecting 25%% of the total\t50%%\t75%%\t100%%") ;
      gxExportTableHeaderLegend (gx, ao, messprintf ("High %s,  collecting over %.0f%% of all reads", GM, 100.0 / (gx->removeLimit ? gx->removeLimit : 1000000)), type) ;
      if (ao2) gxExportTableHeaderLegend (gx, ao2, messprintf ("High %s,  collectingsy over %.0f%% of all reads", GM, 100.0 / (gx->removeLimit ? gx->removeLimit : 1000000)), type) ;

      if (0) gxExportTableHeaderLegend (gx, ao, "Gene", type) ;
      aceOut (ao, "\n") ;
      if (ao2) aceOut (ao2, "\n") ;
    }

  ac_free (h) ;

  return ;
} /* gxExportTableHeader */

/*************************************************************************************/
/* export a table of index, tags, kb, variance ... */
/* order the genes by maximal run index, stored in GZ structure */

static int gxExportTable (GX *gx, int type)
{
  AC_HANDLE h = ac_new_handle () ; 
  ACEOUT ao ;
  int gene, run, igz ;
  Array aa ;
  DC *dc ;
  DDC *ddc ;
  GC *gc ;
  RC *rc ;
  int runMax = arrayMax (gx->runs) ;
  Array allZ = arrayHandleCreate (runMax, float, h) ;
  double z ;
  double ln2 = log(2.0) ;
  BOOL wantGeneGroup = FALSE ;
  const char *types[] = { "expression_index", (! gx->isINTRON ? "reads_aligned_per_gene" : "reads_aligned_per_intron"), "kb_aligned_per_gene", "CV", "sFPKM"
			  , "index_CV", "zebre",  (! gx->isINTRON ? "reads_aligned_per_gene_per_million_raw_reads" : "reads_aligned_per_intron_per_million_raw_reads"), "X", "X", "header"    /* 5->10 */
			  , "nReads", "nReadsOk", "nerr", "a2g", "partial", "orphan" /* 11->16 */
			  , "badtopo", "multi", "multi2"   /* 17->19 */
			  , "X", "X", "X", "X", "X", "X"    /* 20->25 */
  } ;
  Array gza = 0 ;
  GZ *gz ;

  if (type >= 1000)
    {
      type -= 1000 ;
      wantGeneGroup = TRUE ;
    }
  /* export a table one line per gene, one run or group per column  */
  if (type < 0 || type > 19) messcrash ("bad type %d > 19 || < 0 in gxExportTable (GX *gx, int type)", type) ;
  ao = aceOutCreate (gx->outFileName, messprintf (".%s%s.txt", (wantGeneGroup ? "pathways.":""),types[type]), type == 10 ? FALSE :  gx->gzo, h) ;
  gxExportTableHeader (gx, ao, type) ;
  if (type == 10)
    goto done ;
  /* order all genes according to the maximal index among the exported runs */
 
  if (! gx->gza)
    {
      gx->gza = gza = arrayHandleCreate (dictMax(gx->geneDict) + 1, GZ, gx->h) ;
      for (gene = 1, gc = arrp (gx->genes, gene, GC) ; gene <= dictMax (gx->geneDict) ; gc++, gene++)
	{
	  if ( (gx->targeted && !gc->targeted && ! gx->isINTRON) || ( gx->skipEmptyGenes && gc->tags < 1 && ! gc->geneId)) /* force exportation of a gene by setting geneid */
	    continue ;	 
	  if (1)
	    if (gx->captured && ! gc->captured)
	      continue ;
 
	  if (gx->lowVarianceGenes && ! keySetFind (gx->lowVarianceGenes, gene, 0))
	    continue ;
	  if (! gc->isGood)
	    continue ;
	  if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
	    continue ;

	  gz = arrayp (gza, gene, GZ) ;
	  gz->gene = gene ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	    {
	      if (type == 3 && run > 1 && ! rc->runs)
		continue ;
	      if (rc->private)
		continue ;
	      aa = rc->aa ;
	      if (aa)
		{
		  dc = arrayp (aa, gene, DC) ;
		  z = dc->index > 200 ? dc->index + rc->fineTune - 1000 : dc->index ;
		  if (! dc->tags && z == 0) { dc->isLow = TRUE ; z = rc->zeroIndex ; } 
		  if (z > gz->index) gz->index = z ;
		}
	    }
	}
      gx->gza = gza ;
    }
  gza = gx->gza ;
  arraySort (gza, gzOrder) ;
  
  for (igz = arrayMax (gza), gz = arrp (gza, 0, GZ) ; igz > 0 ; igz--, gz++)
    {
      float zSum = 0, zMax2 = -9999 ;
      gene = gz->gene ;
      if (! gene)
	continue ;
      gc = arrp (gx->genes, gene, GC) ;
      if (gc->gMix && gc->gMix != gene)
	continue ;
      if (gx->targeted && ! gc->targeted)
	continue ;
      if (1)
	if (gx->captured && ! gc->captured)
	  continue ;

      if (! gc->isGood)
	continue ;
		  

      if (wantGeneGroup && ! gc->Group_level)
	continue ;
      if (! wantGeneGroup && gc->Group_level)
	continue ;
      if (gx->lowVarianceGenes && ! keySetFind (gx->lowVarianceGenes, gene, 0))
	continue ;
 
      if (type == 6) /* need to check on all genes on each run to eliminate the gene */
	{
	  BOOL ok = FALSE ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC); ! ok && run < runMax ; rc++, run++) 
	    if (rc->aa  && rc->antiAa)
	      {
		DC *antiDc =  arrayp (rc->antiAa, gene, DC) ;
		DC *dc     =  arrayp (rc->aa,     gene, DC) ;
		double zi = dc ? dc->tags : 0 ;
	 	double antiZi = antiDc ? antiDc->tags : 0 ;
		double z = (antiZi > 20 ?  antiZi/(zi + antiZi) : 0) ;
		if (z > .9830) ok = TRUE ;
	      }
	  if (! ok)
	    continue ;
	}
      
      aceOutf (ao, "\n%s", noX__(dictName (gx->geneDict, gene), !gx->isINTRON )) ;

      if (gx->method) aceOutf (ao, "\t%s", gx->method) ;
      if (gx->hasGeneId) 
	{
	  if (gc->geneGroup)
	    aceOutf (ao, "\tGroup") ;
	  else
	    aceOutf (ao, "\t%s", gc->geneId ? stackText (gx->info, gc->geneId) : "") ;
	}

 
      if (0 && gx->hasGeneNm)
	{
	  aceOutf (ao,  "\t") ;
	  if (gc->nmid)
	    {
	      char *buf = strnew (stackText (gx->info, gc->nmid), 0) ;
	      
	      if (buf && *buf)
		aceOutf (ao, "%s", buf) ;

	      ac_free (buf) ;
	    }
	  else if (gc->geneId)
	    {
	      int gg = 0 ;
	      char *cp, *cq, *cr ;
	      char *buf = strnew (stackText (gx->info, gc->geneId), 0) ;
	      int jj = 0 ;
	      
	      cp = cq = buf ;
	      while (cq)
		{
		  cq = strstr (cp, ";") ;
		  if (cq) *cq++ = 0 ; 
		  if (dictFind (gx->geneIdDict, cp, &gg))
		    {
		      GC *gcid = arrp (gx->geneIds, gg, GC) ;
		      cr = stackText (gx->info, gcid->nmid) ;
		      
		      if (cr)
			aceOutf (ao, "%s%s"
				 , jj++ ? ";" : ""
				 , cr
				 ) ;
		    }
		  cp = cq ;
		}
	      ac_free (buf) ;
	    }
	}

      if (gx->hasGeneLength) aceOutf (ao, "\t%d", gx->isMA ? gx->isMA : gc->length) ;
      if (0) aceOutf (ao, "\thasFromGene=%d isINTRON=%d", gx->hasFromGene, gx->isINTRON) ; 
      if (gx->hasKnownDonorAcceptor)
	{
	  if (gc->intronMrna)
	    aceOutf (ao, "\tKnown") ;
	  else
	    switch (gc->hasKnownDonorAcceptor)
	      {
	      case 0: aceOutf (ao, "\tNew") ; break ;
	      case 1:  aceOutf (ao, "\tNew_D") ; break ;
	      case 2:  aceOutf (ao, "\tNew_A") ; break ;
	      case 3:  aceOutf (ao, "\tNew_AD") ; break ;
	    }
	}
      if (gx->hasFromGene)
	{
	  if (gx->isINTRON && gc->isIntron)
	    aceOutf (ao, "\t%s", gc->intronGene ? dictName (gx->geneDict, gc->intronGene) : "") ;
	  else
	    aceOutf (ao, "\t%s", gc->fromGene ? noX__(stackText (gx->info, gc->fromGene), TRUE) : "") ;
	}
      if (gx->hasFromTranscript) aceOutf (ao, "\t%s", gc->fromTranscript ? noX__(stackText (gx->info, gc->fromTranscript), TRUE) : "") ;
      if (gx->hasIntronType)
	{
	  if (gx->isINTRON && gc->isIntron && gx->affyDict)
	    aceOutf (ao, "\t%s", gc->intronType ? dictName (gx->affyDict, gc->intronType) : "") ;
	}
      if (gx->isINTRON && gc->isIntron)
	aceOutf (ao, "\t%d", gc->intronMrna) ;
      if (gx->hasGeneAffy) aceOutf (ao, "\t%s", gc->affy ? stackText (gx->info, gc->affy) : "") ;
      if (gx->hasCapture)
	{
	  KEYSET ks = gc->captures ;
	  int i ;
	  char *sep = "" ;
	  aceOutf (ao, "\t") ;
	  if (gx->isINTRON)
	    {
	      KEY i2g = gc->intronGene ;
	      GC *gc2 = i2g ? arrayp (gx->genes, i2g, GC) : 0 ;
	      ks = gc2 ? gc2->captures : 0 ;
	    }

	  if (ks)
	    for (i = 0 ; i < keySetMax (ks) ; i++)
	      {
		int capture = keySet (ks, i) ;
		if (capture)
		  {
		    const char *ccp = dictName (gx->captureDict, capture) ;
		    if (ccp && (1 || ccp[2] == 0)) 
		      {
			aceOutf (ao, "%s%s", sep, ccp) ;
			sep = "," ;
		      }
		  }
	      }
	}
      if (gx->hasCaptureTouch)
	{
	  KEYSET ks = gc->captures ;
	  KEYSET ksT = gc->capturesTouch ;
	  int i ;
	  char *sep = "" ;
	  aceOutf (ao, "\t") ;
	  if (gx->isINTRON)
	    {
	      KEY i2g = gc->intronGene ;
	      GC *gc2 = i2g ? arrayp (gx->genes, i2g, GC) : 0 ;
	      ks = gc2 ? gc2->captures : 0 ;
	      ksT = gc2 ? gc2->capturesTouch : 0 ;
	    }
	  if (ksT)
	    for (i = 0 ; i < keySetMax (ksT) ; i++)
	      {
		int capture = keySet (ksT, i) ;
		if (capture && ! (ks && keySetFind (ks, capture, 0)))
		  {
		    const char *ccp = dictName (gx->captureDict, capture) ;
		    if (ccp && ccp[2] == 0) 
		      {
			aceOutf (ao, "%s%s", sep, ccp) ;
			sep = "," ;
		      }
		  }
	      }
	}
      if (gx->isINTRON)
	{
	  char *cp, *cq, buf[500], *zchrom = 0 ;
	  int z1 = 0, z2 = 0  ;
	  memcpy (buf, dictName (gx->geneDict, gene), 499) ;
	  
	  if (0)
	    {
	      sscanf (buf, "%s__%d_%d", zchrom, &z1, &z2) ;
	      aceOutf (ao, "\t%s\t%s\t%d\t%d"
		       , zchrom ? zchrom : "-"
		       , z2 > z1 ? "-" : "+"
		       , z1, z2
		       ) ;
	    }
	  else
	    {
	      cp = buf ; cq = strstr (cp, "__") ;
	      if (cq) *cq = 0 ;
	      aceOutf (ao, "\t%s", cp ? cp : "") ;
	      cp = 0 ; 
	      if (cq)
		{
		  cp = cq + 2 ;
		  cq = strstr (cp, "_") ;
		  if (cq) *cq = 0 ;
		  z1 = atoi (cp) ;
		}

	      cp = 0 ; 
	      if (cq)
		{
		  cp = cq + 1 ;
		  cq = strstr (cp, "_") ;
		  if (cq) *cq = 0 ;
		  z2 = atoi (cp) ;
		}

	      aceOutf (ao, "\t%s\t%d\t%d"
		       , z2 > z1 ? "-" : "+"
		       , z1, z2
		       ) ;
	    }
	}
      else if (gx->hasGeneChrom) 
	{
	  if (gc->chrom)
	    aceOutf (ao, "\t%s\t%c\t%d\t%d", dictName (gx->chromDict, gc->chrom), gc->strand, gc->a1, gc->a2) ;
	  else
	    aceOutf (ao, "\t\t\t\t") ;
	}
      if (gx->hasGeneType) aceOutf (ao, "\t%s", gc->gene_type ? stackText (gx->info, gc->gene_type) : "") ;
      if (gx->hasGeneTitle) aceOutf (ao, "\t%s", gc->title ? stackText (gx->info, gc->title) : noX__ (dictName (gx->geneDict, gene), !gx->isINTRON)) ;
      aceOutf (ao, "\t%s", noX__(dictName (gx->geneDict, gene), !gx->isINTRON)) ;
      if (gx->method) aceOutf (ao, ":%s", gx->method) ;

      if (gx->htmlSpecies) 
	aceOutf (ao, "\thttps://www.ncbi.nlm.nih.gov/AceView/db=%s&q=%s"
		 , gx->htmlSpecies
		 , noX__(dictName (gx->geneDict, gene), TRUE)
		 ) ;
      aceOutf (ao, "\t%s", noX__(dictName (gx->geneDict, gene), !gx->isINTRON )) ;
      
      zMax2 = -9999 ;
      zSum = 0 ;
      arrayMax (allZ) = 0 ;
      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	{
	  double z = -100, z1 ;
	  long int zi ;

	  float zMax = -9999 ;
	  if (type == 3 && run > 1 && ! rc->runs)
	    continue ;
	  if (rc->private)
	    continue ;
	  if (type == 5 && ! rc->selectedVariance)
	    continue ;
	  aa = rc->aa ;
	  if (! aa  || (type == 6 && ! rc->antiAa))
	    continue ;
	  if (1)
	    {
	      dc = arrayp (aa, gene, DC) ;
	      ddc = rc->daa ? arrayp (rc->daa, gene, DDC) : 0 ;
	      switch (type)
		{
		case 0:
		  z = dc->index > 200 ? dc->index + rc->fineTune - 1000 : dc->index ;
		  if (! dc->tags && z == 0) { dc->isLow = TRUE ; z = rc->zeroIndex ; } 
		  aceOutf (ao, "\t%s%.2f", dc->isLow ?  (dc->tags ? "00" : "000")  : "",  z - FIX) ;
		  if (z - FIX > zMax) zMax = z-FIX ;
		  break ;
		case 4:  /* sFPKM */
		  z = dc->index > 200 ? dc->index + rc->fineTune - 1000 : dc->index ;
		  if (! dc->tags && z == 0) { dc->isLow = TRUE ; z = rc->zeroIndex ; } 
		  {
		    double z1 ;
		    z = z1 = exp (z*log(2.0))/1000.0 ; 
		    aceOutf (ao, "\t%s%.2f", dc->isLow ? (dc->tags ? "00" : "000")  : "",  z1) ;
		    if (z1 > zMax) zMax = z1 ;
		  }
		  break ;
		case 1:
		  z = zi = dc->tags ;
		  aceOutf (ao, "\t%s%ld",  dc->isLow ?  (dc->tags ? "00" : "000")  : "", zi) ;
		  if (zi > zMax) zMax = zi ;
		  break ;
		case 7: /* reads aligned per_million raw reads */
		  z = dc->tags ; z = 1000000.0 * z / (rc->accepted_tags ? rc->accepted_tags : 1000000) ;
		  zi = z ;
		  aceOutf (ao, "\t%s%.3f",  dc->isLow ?  (dc->tags ? "00" : "000")  : "", z) ;
		  if (z > zMax) zMax = z ;
		  break ;
		case 2:
		  z = dc->kb ;
		  aceOutf (ao, "\t%.2f", z) ;
		  if (z > zMax) zMax = z ;
		  break ;
		case 3:
		  z = dc->variance ;  
		  /* dc->CV = sqrt (2^variance - 1) , where variance = sigma^2 , (CF wikipedia) */
		  z = sqrt (exp (z *ln2) - 1.0) ;
		  z1 = dc->index > 200 ? dc->index + rc->fineTune - 1000 : dc->index ;
		  if (z < 0) 
		    aceOutf (ao, "\tNA") ;
		  else
		    {
		      aceOutf (ao, "\t%s%.3f", dc->isLow ?  (dc->tags ? "0000" : "00000")  : (z1 < 10 ? "000" : (z1 < 12 ? "00" : "")), z) ;
		      if (z > zMax) zMax = z ;
		    }
		  break ;
		case 5:
		  z = z1 = dc->index > 200 ? dc->index + rc->fineTune - 1000 : dc->index ;
		  if (! dc->tags) { dc->isLow = TRUE ; z = rc->zeroIndex ; } 
		  aceOutf (ao, "\t%s%.2f", dc->isLow ?  (dc->tags ? "00" : "000")  : "",  z - FIX) ;
		  z = dc->variance ;
		  /* dc->CV = sqrt (2^variance - 1) , where variance = sigma^2 , (CF wikipedia) */
		  z = sqrt (exp (z *ln2) - 1.0) ;
		  if (z < 0) 
		    aceOutf (ao, "\tNA") ;
		  else
		    {
		      aceOutf (ao, "\t%.3f", z)  ;
		      if (z > zMax) zMax = z ;
		    }
		  break ;
		case 6:
		  zi = dc->tags ;
		  {
		    Array antiAa = rc->antiAa ;
		    DC *antiDc =  (antiAa ? arrayp (antiAa, gene, DC) : 0) ;
		    long int antiZi = (antiDc ? antiDc->tags : 0) ;
		    float z = (zi + antiZi > 10 ?  (100.0 * antiZi)/(zi + antiZi + .0001) : -10) ;
		    if (1) aceOutf (ao, "\t%ld:%ld::%.1f", zi,antiZi, z) ;
		    if (0) aceOutf (ao, "\t%.1f", z) ;
		    if (z > zMax) zMax = z ;
		  }
		  break ;
		case 14:
		  if (ddc)
		    {
		      z = ddc->a2g ;
		      aceOutf (ao, "\t%.2f", z) ;
		      if (z > zMax) zMax = z ;
		    }
		  else
		    aceOutf (ao, "\tNA") ;
		  break ;
		case 15:
		  if (ddc)
		    {
		      z = ddc->partial ;
		      aceOutf (ao, "\t%.2f", z) ;
		      if (z > zMax) zMax = z ;
		    }
		  else
		    aceOutf (ao, "\tNA") ;
		  break ;
		case 16:
		  if (ddc)
		    {
		      z = ddc->orphan ;
		      aceOutf (ao, "\t%.2f", z) ;
		      if (z > zMax) zMax = z ;
		    }
		  else
		    aceOutf (ao, "\tNA") ;
		  break ;
		}
	    }
	  else	if (1 || gx->runListFileName)
	    aceOutf (ao, "\t") ;
	  if (run >= 2 && rc->groupLevel == 0)
	    {
	      if (z > 0)
		{
		  zSum += z ;
		  array (allZ, arrayMax(allZ), float) = z ;
		}
	      if (zMax > zMax2)
		zMax2 = zMax ;
	    }
	}
      aceOutf (ao, "\t%.2f", (zMax2  > -100.0) ? zMax2 : -100.0) ; 
      aceOutf (ao, "\t%.2f", (zSum > 0 && zMax2 > 0) ? (100 * zMax2/zSum) : -100.0) ; 
      arraySort (allZ, floatOrder) ;
      if (1)
	{ /* export the number of runs needed to reach 25% 50% 75% of the total */
	  float s = 0, z25 = zSum/4.0, z50 = 2.0 * z25, z75 = 3.0 * z25, z100 = .9999 * zSum ;
	  int ix = arrayMax (allZ) ;
	  int i, k, i25, i50, i75, i100 ;
	  i25 = i50 = i75 = i100 = 0 ;
	  for(k = 1, i = ix - 1 ; i >= 0 ; k++, i--)
	    {
	      s += arr (allZ, i, float) ;
	      if (s >= z25 && ! i25) i25 = k ;
	      if (s >= z50 && ! i50) i50 = k ;
	      if (s >= z75 && ! i75) i75 = k ;
	      if (s >= z100 && ! i100) { i100 = k ; break ; }
	    }
	  aceOutf (ao, "\t%d\t%d\t%d\t%d", i25, i50, i75, i100) ;
	}
    }

 done:
  aceOut (ao, "\n") ;
  ac_free (h) ;

  return 0 ;
} /* gxExportTable */

/*************************************************************************************/
/* export a gene ace file including groups counts and index computed on the fly */
static int gxExportGeneAceFileSummary (GX *gx, int type)
{
  AC_HANDLE h = ac_new_handle () ; 
  ACEOUT ao ;
  RC *rc ;
  int run, runMax = arrayMax (gx->runs) ;
  Array aa = 0 ;
  const char *target = gx->target_class ;

  if (! target || ! *target) target = "X" ;
  if (strchr ( target, '_')) 
    target = strchr ( target, '_') + 1 ;
  ao = aceOutCreate (gx->outFileName, ".withIndex.ace", gx->gzo, h) ;
  aceOutDate (ao, "//", gx->title) ;

  if (! gx->isINTRON && ! gx->isTranscript)
    {
      const char *target2 = target ;
      char subsamp[256] ;

      memset (subsamp, 0, sizeof(subsamp)) ;
      if (gx->subsample)
	sprintf (subsamp, "_%dk", gx->subsample) ;

      if (gx->captured)
	target2 = hprintf (h, "%s.%s", target, gx->captured) ;
      
      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	{
	  aa = rc->aa ;
	  if (!aa)
	    continue ;

	  aceOutf (ao, "Ali%s \"%s\"\n", subsamp, dictName (gx->runDict, run)) ;
	  if (rc->Genes_with_index)
	    aceOutf (ao, "Zero_index %s %.1f\nLow_index  %s %.1f\nCross_over_index  %s %.1f\nGenes_touched  %s %d\nGenes_with_index %s  %d\nGenes_with_index_over_10 %s  %d\nGenes_with_index_over_12  %s %d\nGenes_with_index_over_15 %s  %d\nGenes_with_index_over_18  %s %d\nGenes_with_index_over_21  %s %d\nGenes_with_index_over_24  %s %d\nGenes_with_index_over_27  %s %d\nGenes_with_index_over_30  %s %d\n"
		     , target2 , rc->zeroIndex
		     , target2 , rc->Low_index
		     , target2 , rc->NA_crossover_index
		     , target2 , rc->Genes_with_index
		     , target2 , rc->Genes_with_reliable_index
		     , target2 , rc->Genes_with_index_over_10
		     , target2 , rc->Genes_with_index_over_12
		     , target2 , rc->Genes_with_index_over_15
		     , target2 , rc->Genes_with_index_over_18
		     , target2 , rc->Genes_with_index_over_21
		     , target2 , rc->Genes_with_index_over_24
		     , target2 , rc->Genes_with_index_over_27
		     , target2 , rc->Genes_with_index_over_30
		     ) ;
	  
	      aceOutf (ao, "Reads_in_genes %s %.0f\n"
		       , target2, rc->nReads
		       ) ;
	      aceOutf (ao, "Reads_in_pairs_in_genes %s %.0f\n"
		       , target2, rc->nReadsOk
		       ) ;
	      aceOutf (ao, "Orphan_reads_in_genes %s %d\n"
		       , target2, rc->orphan 
		       ) ;
	      aceOutf (ao, "Multi_aligned_reads_in_genes %s %d\n"
		       , target2, rc->multi
		       ) ;
	      aceOutf (ao, "Mismatches_in_genes %s %ld\n"
		       , target2, rc->nerr
		       ) ;
	      aceOutf (ao, "A2G_in_genes %s %d\n"
		       , target2, rc->a2g
		       ) ;

	  if (rc->captured_genes_with_index)
	    {	      
	      aceOutf (ao, "Targeted_genes %s %d\nGenes_touched  Captured_%s %d\nGenes_with_index Captured_%s  %d\nGenes_with_index_over_10 Captured_%s  %d\nGenes_with_index_over_12  Captured_%s %d\nGenes_with_index_over_15 Captured_%s  %d\nGenes_with_index_over_18  Captured_%s %d\nGenes_with_index_over_21  Captured_%s %d\nGenes_with_index_over_24  Captured_%s %d\nGenes_with_index_over_27  Captured_%s %d\nGenes_with_index_over_30  Captured_%s %d\n-D High_genes  Captured_%s \n"
		       , dictName (gx->captureDict, rc->capture) , keySet (gx->targeted_genes, rc->capture)
		       , target2 , rc->captured_genes_with_index
		       , target2 , rc->captured_genes_with_reliable_index
		       , target2 , rc->captured_genes_with_index_over_10
		       , target2 , rc->captured_genes_with_index_over_12
		       , target2 , rc->captured_genes_with_index_over_15
		       , target2 , rc->captured_genes_with_index_over_18
		       , target2 , rc->captured_genes_with_index_over_21
		       , target2 , rc->captured_genes_with_index_over_24
		       , target2 , rc->captured_genes_with_index_over_27
		       , target2 , rc->captured_genes_with_index_over_30
		       , target2
		       ) ;
	      aceOutf (ao, "Reads_in_genes Captured_%s %.0f\n"
		       , target2, rc->CnReads
		       ) ;
	      aceOutf (ao, "Reads_in_pairs_in_genes Captured_%s %.0f\n"
		       , target2, rc->CnReadsOk
		       ) ;
	      aceOutf (ao, "Orphan_reads_in_genes Captured_%s %d\n"
		       , target2, rc->Corphan 
		       ) ;
	      aceOutf (ao, "Multi_aligned_reads_in_genes Captured_%s %d\n"
		       , target2, rc->Cmulti
		       ) ;
	      aceOutf (ao, "Mismatches_in_genes Captured_%s %ld\n"
		       , target2, rc->Cnerr
		       ) ;
	      aceOutf (ao, "A2G_in_genes Captured_%s %d\n"
		       , target2, rc->Ca2g
		       ) ;

	      aceOutf (ao, "nh_Ali Captured_%s %.0f seqs %.0f tags %.0f kb_aligned %.2f bp\n"
		       , target2
		       , rc->captured_seqs
		       , rc->captured_tags
		       , rc->captured_kb
		       , 1000*rc->captured_kb/(1+rc->captured_tags)
		       ) ;
	    }
	  if (rc->runs && rc->groupLevel > 0 && ! rc->addCounts)
	    aceOutf (ao, "Genes_expressed_in_at_least_one_run %s %d\n",
		     target2, rc->Genes_with_one_reliable_index) ;
	  aceOutf (ao, "Mb_aligned %f\n",  rc->anyTargetKb/1000) ;
	  aceOutf (ao, "Mb_in_genes %s %f\n", target2, rc->targetKb/1000) ;
	  aceOutf (ao, "Mb_in_genes_with_GeneId_minus_high_genes  %s %f\n", target2, rc->prunedTargetKb/1000) ;
	  aceOutf (ao, "Mb_in_high_genes %s %f\n", target2, rc->bigGenesKb/1000) ;
	  if (rc->captured_kb)
	    aceOutf (ao, "Mb_in_high_genes Captured_%s %f\n", target2, rc->capturedBigGenesKb/1000) ;
	  if (rc->captured_kb/1000)
	    aceOutf (ao, "Mb_in_genes Captured_%s %f\n", target2, rc->captured_kb/1000) ;
	  aceOutf (ao, "Intergenic  %f kb\n", rc->intergenicKb) ;
	  if (rc->bigGenes)
	    gxBigGenes (ao, gx, rc, target2, TRUE) ;
	  
	  aceOutf (ao, "\n") ;
	}
 
      if (! gx->isINTRON && gx->isTranscript)
	for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	  {
	    if (rc->private)
	      continue ;
	    aa = rc->aa ;
	    if (!aa)
	      continue ;
	    aceOutf (ao, "Ali \"%s\"\n", dictName (gx->runDict, run)) ;
	    aceOutf (ao, "mRNA_with_index %s %d\nmRNA_expressed_in_at_least_one_run %s %d\n\n"
		     , target2, rc->Genes_with_reliable_index
		     , target2, rc->runs && ! rc->addCounts ?  rc->Genes_with_one_reliable_index : rc->Genes_with_reliable_index
		     ) ;
      }
    }
  ac_free (h) ;
  return 1 ;
}  /* gxExportGeneAceFileSummary */

/*************************************************************************************/
/* export a gene ace file including groups counts and index computed on the fly */
static int gxExportGeneAceFile (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ; 
  ACEOUT ao ;
  int gene, run, nn = 0 ;
  Array aa, daa ;
  DC *dc ;  DDC *ddc = 0 ;
  GC *gc ;
  RC *rc ;
  int runMax = arrayMax (gx->runs) ;
  int geneMax = arrayMax (gx->genes) ;

  geneMax = dictMax (gx->geneDict) ;
  ao = aceOutCreate (gx->outFileName, ".ace", gx->gzo, h) ;
  aceOutDate (ao, "//", gx->title) ;

  for (gene = 1, gc = arrp (gx->genes, gene, GC) ; gene <= geneMax ; gene++, gc++) 
    {
      if (gx->targeted && ! gc->targeted)
	continue ;
      if (0)
	if (gx->captured && ! gc->captured)
	  continue ;
      
      if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
	continue ;

      if (! gc->isGood)
	continue ;
		  
      aceOutf (ao, "%s %s\n", gx->isINTRON ? "Intron" : (gx->isTranscript ? "Transcript" : (gc->geneGroup ? "pathways" : "Gene")), ac_protect (dictName (gx->geneDict, gene), h)) ;
      for (run = 1, rc = arrayp (gx->runs, run, RC); run < runMax ; rc++, run++) 
	{
	  aa = rc->aa ;
	  daa = rc->daa ;
	  if (!aa)
	    continue ;
	  if (rc->private)
	    continue ;
	  dc = arrayp (aa, gene, DC) ;
	  ddc = daa ? arrayp (daa, gene, DDC) : 0 ;
	  if (dc->tags >= 0 && dc->index > -900)
	    {
	      nn++ ;
	      if (rc->runs)
		{
		  aceOutf (ao, "Group_%s  %s %.2f %.2lf seqs %.2lf tags %.2lf kb"
			   , gx->unique ? "U" : "nU"
			   , dictName (gx->runDict, run)
			   , dc->index > 200 ? dc->index + rc->fineTune - 1000 : dc->index - 1000
			   , dc->seqs, dc->tags, dc->kb
			   ) ;

		  if (ddc)
		    {
		      aceOutf (ao, " %.2lf reads %.2lf compRead"
			       , ddc->nReads, ddc->nReadsOk
			       ) ;
		      aceOutf (ao, " %ld err  %d a2g"
			       , ddc->nerr, ddc->a2g
			       ) ;
		      aceOutf (ao, " %d partial %d orphan %d badTopo"
			       , ddc->partial, ddc->orphan, ddc->badTopo
			       ) ;
		      aceOutf (ao, " %d Multi %d AmbStrand"			      
			       , ddc->multi, ddc->multi2
			       ) ;
		    }
		  else
		    aceOutf (ao, " 0 reads 0 compRead 0 err 0 a2g 0 partial 0 orphan 0 badTopo 0 Multi 0 AmbStrand") ;

		  aceOutf (ao, " 0 %s"
			   , dc->isLow ? (dc->tags ? "NA" : "NE") : "ok"
			   ) ;
		  if (! dc->isLow && dc->variance >= 0)
		    aceOutf (ao, " %.4f sigma", sqrt(dc->variance)) ;
		  aceOut (ao, "\n") ;
		}
	      else
		{
		  if (0 && gx->isINTRON) continue ;

		  aceOutf (ao, "Run_%s %s %.2f %.2f seqs %.2f tags %.2f kb "
			   , gx->unique ? "U" : "nU"
			   , dictName (gx->runDict, run)
			   , dc->index > 200 ? dc->index + rc->fineTune - 1000 : dc->index - 1000
			   , dc->seqs, dc->tags, dc->kb
			   
			 ) ;
		  if (ddc)
		    {
		      aceOutf (ao, " %.2lf reads %.2lf compRead"
			       , ddc->nReads, ddc->nReadsOk
			       ) ;
		      aceOutf (ao, " %ld err  %d a2g"
			       , ddc->nerr, ddc->a2g
			       ) ;
		      aceOutf (ao, " %d partial %d orphan %d badTopo"
			       , ddc->partial, ddc->orphan, ddc->badTopo
			       ) ;
		      aceOutf (ao, " %d Multi %d AmbStrand"			      
			       , ddc->multi, ddc->multi2
			       ) ;
		    }
		  else
		    aceOutf (ao, " 0 reads 0 compRead 0 err 0 a2g 0 partial 0 orphan 0 badTopo 0 Multi 0 AmbStrand") ;
		  aceOutf (ao, " 0 %s\n"
			   , dc->isLow ? (dc->tags ? "NA" : "NE") : "ok"
			 ) ;
		}
	    }
	}
      aceOut (ao, "\n") ;
    }

  ac_free (h) ;
  return nn ;
}  /* gxExportGeneAceFile */

/*************************************************************************************/

int rcAlphaOrder (const void *a, const void *b)
{
  const RC *up = (const RC *)a, *vp = (const RC *)b ;
  double z ;

  z = up->alpha - vp->alpha ;
    if (z > 0) return 1 ;
    if (z < 0) return -1 ;

  return up->run - vp->run ;
} /* rcAlphaOrder */

/*************************************************************************************/
/* check if the variance depends of the signal */
static int gxExportMeanversusVariancePlot (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ; 
  ACEOUT ao = aceOutCreate (gx->outFileName, ".signal2variance.txt", FALSE, h) ;
  int N, gene, run ;
  Array aa ;
  DC *dc ;
  RC *rc ;
  int runMax = arrayMax (gx->runs) ;
  int level ;
  double X, X2, z, w ;
  
  /* export a table one line per gene, one run or group per column */
  aceOutDate (ao, "###", gx->title) ;
  gxExportTableHeader (gx, ao, 100) ;
  
  for (level = 0 ; level < 30 ; level++)
    {
      aceOutf (ao, "\n%d", level) ;
      for (run = 1, rc = arrayp (gx->runs, run, RC) ; run < runMax ; rc++, run++) 
	{
	  if (rc->private)
	    continue ;
	  w = 0 ;
	  aa = rc->aa ;
	  if (aa)
	    {
	      N = 0 ; X = X2 = 0 ;  
	      if (level > 0)
		{
		  for (gene = 0, dc = arrayp (aa, gene, DC)  ; gene < arrayMax (aa) ; dc++, gene++)
		    {
		      z = dc->index > 200 ? dc->index + rc->fineTune - 1000 : dc->index - 1000 ;
		      if (z <= 0) z = 0 ;
		      if (z >= level && z < level+1)
			{
			  N++ ; X += z ; X2 += z * z ;
			}
		    }
		  if (N)
		    w = X2/N - X*X/(N * N) ;
		}
	      
	      aceOutf (ao, "\t%.2f", w) ;
	    }
	}	
    }

  aceOut (ao, "\n") ;
  ac_free (h) ;

  return 0 ;
} /* gxExportMeanversusVariancePlot */

/*************************************************************************************/
/* if ks1 ks2 are distinct, return ks1 ks2
 * if ks2 is included in ks1 return ks2 and {ks1 minus ks2} and vice versa
 * otherwise return ks1 minus intersect and ks2 minus intersect
 */
static void keySetExclude (KEYSET *ks1p, KEYSET *ks2p, KEYSET ks1, KEYSET ks2, AC_HANDLE h)
{
  KEYSET a, b ;

  a = keySetMINUS (ks1, ks2) ;
  b = keySetMINUS (ks2, ks1) ;
  if (keySetMax (a) && keySetMax(b))
    {
      *ks1p = arrayHandleCopy (a, h) ;
      *ks2p = arrayHandleCopy (b, h) ;
    }
  else if (keySetMax(a))
    {
      *ks1p = arrayHandleCopy (a, h) ;
      *ks2p = arrayHandleCopy (ks2, h) ;
    }
  else 
    {
      *ks1p = arrayHandleCopy (ks1, h) ;
      *ks2p = arrayHandleCopy (b, h) ;
    }
  keySetDestroy (a) ;
  keySetDestroy (b) ;
} /* keySetExclude */

/*************************************************************************************/
/*************************************************************************************/
/* count the true and false positives and the AUC and optimize the MCC
 * assuming the training groups run1 run2 are true, and respecting the alpha order in runs
 * we then construct the beta1 FN/N and beta2 FP/P coeficient of certitude
< */
static Array gxTrueFalse_AUC_MCC (GX *gx, ACEOUT ao, Array rws, int run1, int run2, int run3, int run4, int pass)
{
  AC_HANDLE h = ac_new_handle () ; 
  int ir, irMax = arrayMax (rws), nn = 0 ;
  int s, w, u = 0, f = 0 ;
  int bestjj, jj, bestTpTn ;
  int a[irMax], u2f[irMax+1] ;
  int ttp[irMax], ttn[irMax], tfp[irMax], tfn[irMax], sp[irMax], sn[irMax], accu[irMax], jj2ir[irMax], ir2jj[irMax] ;
  double bestmcc, mcc[irMax], b, bonus[irMax] ;
  float alpha[irMax], da = 0 ;
  RW *rw, *rwa ;
  RC *rc ;
  double AUC, beta, xx[irMax], yy[irMax], beta1Min, beta1Max, beta2Min, beta2Max, tbeta1[irMax], tbeta2[irMax] ;
  BOOL emptyAlpha = TRUE ;
  BOOL debug = FALSE ;
  KEYSET ks1a = 0 ;
  KEYSET ks2a = 0 ;
  KEYSET ks1 = 0 ;
  KEYSET ks2 = 0 ;

  memset (xx, 0, sizeof(xx)) ;
  memset (yy, 0, sizeof(yy)) ;
  if (run1 && run2)
    {
      RC *rc1 = 0, *rc2 = 0 ;

      rc1 = arrayp (gx->runs, run1, RC) ;
      rc2 = arrayp (gx->runs, run2, RC) ;
  
      if (rc1) ks1 = rc1->runs ;
      if (rc2) ks2 = rc2->runs ;
    }

  keySetExclude (&ks1a, &ks2a, ks1, ks2, h) ;

  if (0 && strstr(dictName(gx->runDict, run1), "avorable_Tr") > 0)
    debug = TRUE ;
   
  memset (a, 0, sizeof (a)) ;
  memset (u2f, 0, sizeof (u2f)) ;
  memset (ttp, 0, sizeof (ttp)) ;
  memset (ttn, 0, sizeof (ttn)) ;
  memset (tfp, 0, sizeof (tfp)) ;
  memset (tfn, 0, sizeof (tfn)) ;
  memset (tbeta1, 0, sizeof (tbeta1)) ;
  memset (tbeta2, 0, sizeof (tbeta2)) ;
  memset (jj2ir, 0, sizeof (jj2ir)) ;
  memset (ir2jj, 0, sizeof (ir2jj)) ;
  memset (alpha, 0, sizeof (alpha)) ;
  memset (sp, 0, sizeof (sp)) ;
  memset (sn, 0, sizeof (sn)) ;
  memset (accu, 0, sizeof (accu)) ;
  memset (mcc, 0, sizeof (mcc)) ;
  memset (bonus, 0, sizeof (bonus)) ;

  for (nn = w = f = s = 0, rw = arrp (rws, 0, RW), ir = 0 ; ir < arrayMax (rws) ; ir++, rw++)
    {
      rc = arrp (gx->runs, rw->run, RC) ;
      if (rc->avoid) continue ;
      if (ks2a && keySetFind (ks2a, rw->run, 0))
	{ w = 1 ; f++ ;                       nn++ ; a[nn] = 1 ;  alpha[nn] = rw->alpha ; jj2ir[nn] = ir ; ir2jj[ir] = nn + 1 ; } 
      if (ks1a && keySetFind (ks1a, rw->run, 0))
	{ w = 0 ; u++ ; s += f ; u2f[u] = f ; nn++ ; a[nn] =-1 ;  alpha[nn] = rw->alpha ; jj2ir[nn] = ir ; ir2jj[ir] = nn + 1 ; }
    }

  emptyAlpha = (nn < 1 || alpha[nn-1] == alpha[0]) ? TRUE : FALSE ;
  if (debug)
    invokeDebugger () ;
  
  if (w == 1)
    {  /*  register the last value */
      u++ ; u2f[u] = f ; 
    }
  u2f[u+1] = f ;
   beta = s ; 
   if (u*f > 0) beta /= (u*f) ; 
   AUC = 100.0 * beta ;
   if (emptyAlpha ) AUC = s = 0 ;

   /* drawing of the ROC curve */
   
   aceOutf (ao, "s=%d uf=%d u=%d f=%d AUC=%.2f\t",s,u*f,u,f,AUC) ;
   aceOutf (ao, "\nFalse positive counts\tFalse positive counts") ; 
   for (int i = 1 ; i <=  u ; i++)
     aceOutf (ao, "\t%d", u2f[i+1]-u2f[i]) ;

   aceOutf (ao, "\nTrue positive counts\t%s", "title") ;
   for (int i = 1 ; i <= u; i++)
     aceOutf (ao, "\t%d", u2f[i]);

   aceOutf (ao, "\nFalse positive rate\tFalse positive rate") ; 
   for (int x = 0 ; x <= 100; x++)
     aceOutf (ao, "\t%.2f", x/100.0) ;
   aceOutf (ao, "\nTrue positive rate\t%s", "title") ;
   u2f[u+1] = u2f[u] ;
   for (int x = 0 ; x <= 100 ; x++)
     {
       float y = 0  ; 
       int i = u*x/100.0 + 1 ;
       if (x == 100)
	 y = 1 ;
       else
	 {
	   float y1 = u2f[i] ;
	   float dy = u2f[i+1] - y1 ;
	   float dx = u * x / 100.0 - i + 1 ;

	   y = y1 ;
	   if (dx>0)
	     y = y + dx * dy / 100.0 ;
	   if (i < -2)
	     aceOutf (ao, "ERROR dx=%f x=%d u=%d ux=%d i=%d dy=%.2f y1=%.2f y=%.2f\n",dx,x,u,u*x,i,dy,y1,y);
	   if (u2f[u] > 0) y = y / u2f[u];
	 }
       aceOutf (ao, "\t%.3f",y);
     }
   aceOutf (ao, "\n");

   /*  look for the threshold with optimal MCC Matthews correlation coefficient */
   bestjj = 0 ; bestmcc = 0 ; bestTpTn = 100000000 ; da = 0 ;
   for (jj = 1 ; jj <= nn ; jj++)
     {
       int i, tp = 0, fp = 0, tn = 0,fn = 0 ;
       /* TP all positive below jj
	* FP all positive above jj
	* TN all positive above jj
	* FN all positive below jj
	*/
       
       for (i = 0 ; i < jj ; i++)
	 {
	   if (a[i] == 1) tp++ ;
	   if (a[i] == -1) fp++ ;
	 }
       for (i = jj ; i <= nn ; i++)
	 {
	   if (a[i] == 1) fn++ ;
	   if (a[i] == -1) tn++ ;
	 }
       ttp[jj] = tp ;
       ttn[jj] = tn ;
       tfp[jj] = fp ;
       tfn[jj] = fn ;
       sn[jj] = 0 ; if (tp>0)   sn[jj] = 100.0*tp/(1.0*tp+fn) ;
       sp[jj] = 0 ; if (tn>0)   sp[jj] = 100.0*tn/(1.0*tn+fp) ;
       /*
	 # Youden's J statistics = intercept of ROC curve with anti diagonal  = sn + sp -1  # not recommended
	 # Gini coefficint = 2 * AUC - 1
	 ppv[jj] = 0 ;  if (tp + fp > 0) tp/(tp + fp) ;   # positive predictive value or Precision
	 npv[jj] = 0 ;  if (tn + fn > 0) tn/(tn + fn) ;   # negative predictive value
	 fdr[jj] = 0 ;  if (fp + tp > 0) fp/(fp + tp) ;   # FDR false discovery rate
	 fpr[jj] = 0 ;  if (fp + tn > 0) fp/(fp + tn) ;   # false positive rate or fall-out
       */
       accu[jj] = 0 ; if (tn + tp > 0) accu[jj] = 100.0*(tn + tp)/(1.0*tn+tp+fn+fp) ;
       mcc[jj] = (tp+fp) ;
       mcc[jj] *= (tp+fn) ;
       mcc[jj] *= (tn+fp) ;
       mcc[jj] *= (tn+fn) ;
       if (mcc[jj] > 0) mcc[jj] = (tp*tn - fp*fn)/sqrt(mcc[jj]) ;
       b = tp + fp  - f ; if (b < 0) b = -b ; b /= sqrt(f) ;
       bonus[jj] = .3 - b/20 ;
       if (mcc[jj]+bonus[jj] > bestmcc) { bestjj = jj ; bestmcc = mcc[jj] + bonus[jj] ; bestTpTn = ttp[jj] + ttn[jj] ; da = alpha[jj] - alpha[jj-1] ; }
     }

   if (bestTpTn) {}  /* for compiler happiness */
#ifdef JUNK
   /* little adjustment looking for a larger alpha gap */
     for (jj = bestjj + 10 ; jj > bestjj - 10 && jj > 1 ; jj--)
       {
	 if (jj >= nn) continue ;
	 if (ttp[jj] + ttn[jj] >= bestTpTn - 1 && a[jj-1] == 1 && a[jj] == -1 && (alpha[jj] - alpha[jj-1]) > da)
	   { da = alpha[jj] - alpha[jj-1] ; bestjj = jj ; break ; }
     }
     beta = bestjj > 0 ?  (alpha[bestjj] + alpha[bestjj-1])/2 : 0 ;
     {
       int i;
       for (rw = arrp (rws, 0, RW), i = 0 ; i < arrayMax (rws) ; i++, rw++) 
	 rw->alpha -= 0*beta ; /* we do not need to shift alpha, since now we shift beta */
     }
#else
     if (da) {}  /* for compiler happiness */
#endif

   gx->AUC = AUC ;
   gx->MCC = 100 * mcc[bestjj] ;
   gx->AUC2 = gx->MCC2 = 0 ;

   /* initialize */
   beta1Min = beta2Min =  1000 ;
   beta1Max = beta2Max = -1000 ;

   for (rw = arrp (rws, 0, RW), ir = 0 ; ir < arrayMax (rws) ; ir++, rw++)
     {  
       double z ;
       rw->beta1 = rw->beta2 = -1000 ;
       
       /* for the runs in the training set, we locate them on the ROC curve */
       jj = ir2jj[ir] ; /* 1 + true_jj */
       if (jj > nn + 1)
	 messcrash ("jj=%d > nn=%d in gxTrueFalse_AUC_MCC", jj, nn) ;
       if (jj > 0) 
	 {
	   int jj3 ;
	   jj-- ; jj3 = jj - 1 ;
	   while (++jj3 <= nn)
	     {
	       z = ttp[jj3] + tfp[jj3] ; 
	       if (z > 0)
		 {
		   tbeta2[jj] = rw->beta2 = ttp[jj3] / z ; /* TP/P [jj] */
		   break ;
		 }
	     }
	   if (beta2Min > rw->beta2) beta2Min = rw->beta2 ;
	   if (beta2Max < rw->beta2) beta1Max = rw->beta2 ;

	   z = ttn[jj] + tfn[jj] ;       tbeta1[jj] = rw->beta1 = z>1 ? ttn[jj] / z : 0 ; /* TN/N [jj] */
	   
	   if (beta1Min > rw->beta1) beta1Min = rw->beta1 ;
	   if (beta1Max < rw->beta1) beta1Max = rw->beta1 ;
	 }
     }


   /* interpolation of the unknown alpha */
   for (rwa = 0, rw = arrp (rws, 0, RW), ir = 0 ; ir < arrayMax (rws) ; ir++, rw++)
     {
       int jr ;
       RW *rwb = 0 ;

       if (rw->beta1 > -999)
	 { rwa = rw ; continue ; }
       if (rwa == 0) /* stabilize on minimal value */
	 {
	   rw->beta1 = beta1Max ;
	   rw->beta2 = beta2Min ;
	   continue ;
	 }
       for (rwb = rw, jr = ir ; jr < arrayMax (rws) ; jr++, rwb++)
	 if (rwb->beta1 > -999)
	   {
	     rw->beta1 = rwa->beta1 ;
	     rw->beta2 = rwa->beta2 ;
	     if (rwb->alpha != rwa->alpha)
	       {
		 rw->beta1 += (rw->alpha - rwa->alpha) * (rwb->beta1 - rwa->beta1)/(rwb->alpha - rwa->alpha) ;
		 rw->beta2 += (rw->alpha - rwa->alpha) * (rwb->beta2 - rwa->beta2)/(rwb->alpha - rwa->alpha) ;
	       }
	     rwa = rw ;
	     break ;
	   }
       if (rw->beta1 < -999)
	 { rw->beta1 = beta1Min ; rw->beta2 = beta2Max ; }
     }
	 
   /* reinject in original runs */
   for (rw = arrp (rws, 0, RW), ir = 0 ; ir < arrayMax (rws) ; ir++, rw++)
     {
       RC *rc1 = arrp (gx->runs, rw->run, RC) ;
       rc1->beta1 = rw->beta1 ;
       rc1->beta2 = rw->beta2 ;
     }

   if (emptyAlpha) 
     {
       beta = 0 ;
     }
   if (debug)
     {
       jj = bestjj ;
       fprintf (stderr, "\n\\n\ntMCC sort Bestj=%d AUC=%.2f %s \t%s\n"
		, jj, 100*beta
		, dictName(gx->runDict, run1)
		, dictName(gx->runDict, run2)
		);
        fprintf (stderr, "TP=%d TN=%d FP=%d FN=%d "
		, ttp[jj],ttn[jj], tfp[jj], tfn[jj]
		);
        fprintf (stderr, "SENS=%d SP=%d ACCURACY=%d MCC=%.1f\n"
		,100*sn[jj], 100*sp[jj], 100*accu[jj], 100*mcc[jj]
		);
       if (1)
	 for (jj = 1 ; jj < nn ; jj++)
	   printf ("MCC jj=%d a=%d tp=%d tn=%d fp=%d fn=%d SENS=%d SP=%d ACCURACY=%d MCC=%.1f SN+SP=%d SN*SP=%d xx=%.2f yy=%.2f\tfx+uy=%.2f\talpha=%.2f\tbeta1=%.2f\tbeta2=%.2f\n"
		  ,jj, a[jj]
		  ,ttp[jj],ttn[jj],tfp[jj],tfn[jj]
		  ,sn[jj], sp[jj], accu[jj],mcc[jj], sn[jj]+sp[jj],sn[jj]*sp[jj]
		  ,xx[jj],yy[jj], f*xx[jj]+u*yy[jj]
		  , alpha[jj], tbeta1[jj], tbeta2[jj]
		  );
     }

   if (run3 && run4) /* in the same coordinate system construct the AUC MCC of a validation by stratification group */
     {
       RC *rc3 = arrp (gx->runs, run3, RC) ;
       RC *rc4 = arrp (gx->runs, run4, RC) ;
       KEYSET ks3a, ks3 = rc3 ? rc3->runs : 0 ;
       KEYSET ks4a, ks4 = rc4 ? rc4->runs : 0 ;
       int tp, fp, tn, fn ;
       double z ;

       if (ks3 && ks4)
	 {
	   keySetExclude (&ks3a, &ks4a, ks3, ks4, h) ;
	   for (nn = w = f = u = s = tp = tn = fp = fn = 0, rw = arrp (rws, 0, RW), ir = 0 ; ir < arrayMax (rws) ; ir++, rw++)
	     {
	       if (keySetFind (ks1, rw->run, 0) || keySetFind (ks2, rw->run, 0)) /* exclude the training set from the test set */
		 continue ;
	       rc = arrp (gx->runs, rw->run, RC) ;
	       if (rc->avoid) continue ;
	       if (keySetFind (ks4a, rw->run, 0))
		 { w = 1 ; f++ ; if(rw->beta1 > rw->beta2) tp++ ; else fn++ ;} 
	       if (keySetFind (ks3a, rw->run, 0))
		 { w = 0 ; u++ ; s += f ; if(rw->beta1 < rw->beta2) tn++ ; else fp++ ;}
	     }
	   if (w == 1)
	     {  /*  register the last value */
	       u++ ; u2f[u] = f ; 
	     }
	   beta = s ; 
	   if (u*f > 0) beta /= (u*f) ; 
	   gx->AUC2 = 100.0 * beta ;

	   z = (tp+fp) ;
	   z *= (tp+fn) ;
	   z *= (tn+fp) ;
	   z *= (tn+fn) ;
	   gx->MCC2 = (z > 0) ? 100 * (tp*tn - fp*fn)/sqrt(z) : 0 ;
	   aceOutf (ao, "%s %s\tmethod=%d tp=%d fp=%d tn=%d fn=%d z=%.1f MCC2=%.1f AUC=%.1f AUC2=%.1f\titer%d\t%s\t%s\n"
		    , dictName(gx->runDict,run1)
		    , dictName(gx->runDict,run2)
		    , gx->geneSelectionMethod
		    , tp, fp, tn, fn, z
		    , gx->MCC2, gx->AUC, gx->AUC2
		    , pass
		    , dictName(gx->runDict,run3)
		    , dictName(gx->runDict,run4)
		    ) ;
	   if (gx->aoAUC2)
	     aceOutf (gx->aoAUC2, "%s %s\tmethod=%d tp=%d fp=%d tn=%d fn=%d z=%.1f MCC2=%.1f AUC=%.1f AUC2=%.1f\titer%d\t%s\t%s\n"
		      , dictName(gx->runDict,run1)
		      , dictName(gx->runDict,run2)
		      , gx->geneSelectionMethod
		      , tp, fp, tn, fn, z
		      , gx->MCC2, gx->AUC, gx->AUC2
		      , pass
		      , dictName(gx->runDict,run3)
		      , dictName(gx->runDict,run4)
		      ) ;
	 } 
     } 
   ac_free (h) ;
  return 0 ;
} /* gxTrueFalse_AUC_MCC */

/*************************************************************************************/

static void gxShowGeneHisto (ACEOUT ao, GX* gx, double *hh, int nRun) 
{
  int i, j ;
  double x, s ;
  
  for (s = 0, i = 0 ; i < 256 ; i++)
    s += hh[i] ;
  for (i = 0 ; i < 260 ; i+= 5)
    {
      x = hh[i] ;
      if (1)
	{
	  for (x = 0, j = -2 ; j <= 2 ; j++)
	    if (i+j >= 0 && i+j < 256) x += hh[i+j] ;
	  if (0 && i == 240)
	    for (j = 243 ; j < 256 ; j++)
	      x +=  hh[i+j] ;
	}
      aceOutf (ao, "\t%.2f", x) ;
    }
  aceOutf (ao, "\t%.2f\t%d", s, nRun) ;
} /* gxShowGeneHisto */

/*************************************************************************************/
/* The run are sorted in rows by their alpha coefficient and export for each run the index of the most characteristic genes */
static int gxExportComparison (GX *gx, COMPARE *compare, Array rws, KEYSET genePlus, KEYSET geneMinus, int run1, int run2, BOOL histo, BOOL wantGeneGroup)
{
  AC_HANDLE h = ac_new_handle () ; 
  ACEOUT ao = 0 ;
  const char *species =  gx->htmlSpecies ? gx->htmlSpecies : "species" ;
  int ii, jj, jpp, gene, genePass, pass, betaPass ;
  int run3 = 0, run4 = 0 ;
  vTXT compareWhat = vtxtHandleCreate (h) ;
  RW *rw ;
  RC *rc ;
  DC *dc ;
  GC *gc ;
  double z ;
  KEYSET genePM ;

  if (!run3)
    { run3 = run1 ; run4 = run2 ; }

  ao = aceOutCreate (gx->outFileName
		     , hprintf (h, ".%s.%s.%s.txt"
				, dictName (gx->compareDict, compare->compare)
				, histo ? "histo" : "beta"
				, wantGeneGroup ? "pathways" : "genes"
				)
		     ,  FALSE, h
		     ) ;
  aceOutDate (ao, "###", gx->title) ;
  aceOutf (ao, "## %s : This classifier is based on the comparison of the groups \t%s\t%s\n"
	   , dictName (gx->compareDict, compare->compare)
	   , dictName (gx->runDict, run1)
	   , dictName (gx->runDict, run2)
	   ) ;
  aceOutf (ao, "## Parameters:  %s  %s R= %.1f %s"
	   , gx->normalize ? "The area of each histogram is normalized to 100" : "The area of each histogram represents the actual number of samples"
	   , gx->geneSelectionMethod == 7 ? "" : messprintf ("geneSelectionMethod=%d",  gx->geneSelectionMethod)
	   , gx->ratio_bound
	   , wantGeneGroup ? "Pathways" : "Genes"
	   ) ;
  
  if (1)
    {
      int nMax = gx->maxGenePlus ;
      int np = 0, nm = 0, iMax = keySetMax (genePlus) ;
      for (int i = 0 ; i < iMax && np < nMax ; i++)
	{
	  int gene = keySet (genePlus, i) ;
	  gc = arrp (gx->genes, gene, GC) ;  
	  if (wantGeneGroup && ! gc->geneGroup)
	    continue ;
	  if (! wantGeneGroup && gc->geneGroup)
	    continue ;
	  np++ ;
	}
      iMax = keySetMax (geneMinus) ;
      for (int i = 0 ; i < iMax && nm < nMax ; i++)
	{
	  int gene = keySet (geneMinus, i) ;
	  gc = arrp (gx->genes, gene, GC) ;  
	  if (wantGeneGroup && ! gc->geneGroup)
	    continue ;
	  if (! wantGeneGroup && gc->geneGroup)
	    continue ;
	  nm++ ;
	}

      aceOutf (ao, "## We required maxGenes\t%d\tand used\t%d\tplus genes and\t%d\tminus genes\n"
	       , gx->maxGenePlus
	       , np, nm
	       ) ;
    }

  aceOutf (ao, "\n##AUC:\t%s\t%.1f\tMCC\t%.1f"
	   , dictName (gx->compareDict, compare->compare)
	   , gx->AUC, gx->MCC
	   ) ;
  {
    const char *ccp = "" ;
    switch (gx->digital)
      {
      case 0:
	ccp = "produit scalaire" ;
	break ;
      case 1:
	ccp = "weigthed classifier" ;
	break ;
      case 2:
	ccp = "fully digital" ;
	break ;
      case 3:
	ccp = "Mobymatic" ;
	break ;
      case 4:
	ccp = "SNP, count by classes" ;
	break ;
      case 5:
	ccp = "produit scalaire avec seuil" ;
	break ;
      }
    aceOutf (ao, "\t%s\n", ccp) ;
  }

  if (compare->predicts)
    for (ii = 0 ; ii < keySetMax (compare->predicts) ; ii++)
      {
	int iCompare2 = keySet (compare->predicts, ii) ;
	COMPARE *compare2 = arrayp (gx->compares, iCompare2, COMPARE) ;
	if (! compare2->runs || keySetMax (compare2->runs) != 2)
	  continue ;
	run3 = keySet (compare2->runs, 0) ;
	run4 = keySet (compare2->runs, 1) ;
	gx->AUC2 = gx->MCC2 = -10 ;
	/* side effect, the last run3 run4 pair is used below when reporting the validating phenotype line may be the AUC shoud be written then */
	gxTrueFalse_AUC_MCC (gx, ao, rws, run1, run2, run3, run4, 1) ;
	fprintf (stderr, "AUC2: %g, MCC2:%g\n", gx->AUC2, gx->MCC2) ;
	aceOutf (ao, "##AUC2:\t%s\tpredicts\t%s\t%.1f\tMCC2\t%.1f\n"
		 , dictName (gx->compareDict, compare->compare)
		 , dictName (gx->compareDict, compare2->compare) 
		 , gx->AUC2, gx->MCC2
		 ) ;
      }

  rc = arrayp (gx->runs, run1, RC) ;
  aceOutf (ao, "Plus genes: overexpressed in\t%s\t%d\truns used", dictName (gx->runDict, run1),rc->runs ? arrayMax(rc->runs): 1); 
  for (jj = 0 ; genePlus && jj < keySetMax (genePlus) ; jj++)
    {
      int gene = keySet (genePlus, jj) ;
      gc = arrp (gx->genes, gene, GC) ;  
      if (wantGeneGroup && ! gc->geneGroup)
	continue ;
      if (! wantGeneGroup && gc->geneGroup)
	continue ;

      aceOutf (ao, "\t%s", gene ? noX__ (dictName (gx->geneDict, gene), !gx->isINTRON) : "") ;
    }
  rc = arrayp (gx->runs, run2, RC) ;
  aceOutf (ao, "\nMinus genes: overexpressed in\t%s\t%d\truns used", dictName (gx->runDict, run2), rc->runs ? arrayMax(rc->runs): 1); 
  for (jj = 0 ; geneMinus && jj < keySetMax (geneMinus) ; jj++)
    {
      int gene = keySet (geneMinus, jj) ;
      gc = arrp (gx->genes, gene, GC) ;  
      if (wantGeneGroup && ! gc->geneGroup)
	continue ;
      if (! wantGeneGroup && gc->geneGroup)
	continue ;

      aceOutf (ao, "\t%s", gene ? noX__ (dictName (gx->geneDict, gene), !gx->isINTRON) : "") ;
    }
 
  for (genePass = 0 ; genePass < 2 ; genePass++)
    {
      genePM = genePass == 0 ? genePlus : geneMinus ;
      rc = arrayp (gx->runs, genePass == 0 ? run1 : run2, RC) ;
      if (! rc->runs || ! arrayMax (rc->runs))
	continue ;
      aceOutf (ao, "\n\n\nGenes:") ;
      vtxtClear (compareWhat) ;
      {
	RC *rc1 = rc ;
	RC *rc2 =  arrayp (gx->runs, genePass == 0 ? run2 : run1, RC) ;
	vtxtPrintf (compareWhat, "overexpressed in\t%s relative to %s using %d versus %d runs"
		    , dictName (gx->runDict, (genePass == 0 ? run1 : run2))
		    , dictName (gx->runDict, (genePass == 0 ? run2 : run1))
		    , (rc1->runs ? arrayMax(rc->runs) : 0)
		    , (rc2->runs ? arrayMax(rc->runs) : 0)
		    ) ;
      }
      if (!histo && genePass == 0)
	{
	  aceOutf (ao, "\n##\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tRun:"); 
	  for (pass = 0 ; pass < 2 ; pass++)
	    for (rw = arrp (rws, 0, RW), ii = 0 ; ii < arrayMax (rws) ; ii++, rw++) 
	      {
		rc = arrp (gx->runs, rw->run, RC) ;
		if (! rc->run || ! rc->aa || rc->private) continue ;
		if ((pass == 0 && rc->runs) || (pass == 1 && ! rc->runs))
		  continue ;
		aceOutf (ao, "\t%s", dictName (gx->runDict, rc->run)) ;
	      }
	  aceOutf (ao, "\n##\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tTitle:") ;
	  for (pass = 0 ; pass < 2 ; pass++)
	    for (rw = arrp (rws, 0, RW), ii = 0 ; ii < arrayMax (rws) ; ii++, rw++) 
	      {
		rc = arrp (gx->runs, rw->run, RC) ;
		if (! rc->run || ! rc->aa || rc->private) continue ;
		if ((pass == 0 && rc->runs) || (pass == 1 && ! rc->runs))
		  continue ;
		aceOutf (ao, "\t%s",  rc->title ? stackText (gx->info, rc->title) : "") ;
	      }
	  aceOutf (ao, "\n##\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tRunId:") ;
	  for (pass = 0 ; pass < 2 ; pass++)
	    for (rw = arrp (rws, 0, RW), ii = 0 ; ii < arrayMax (rws) ; ii++, rw++) 
	      {
		rc = arrp (gx->runs, rw->run, RC) ;
		if (! rc->run || ! rc->aa || rc->private) continue ;
		if ((pass == 0 && rc->runs) || (pass == 1 && ! rc->runs))
		  continue ;
		aceOutf (ao, "\t%s",  rc->runid ? stackText (gx->info, rc->runid) : "") ;
	      }
	  
	  aceOutf (ao, "\n##\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tSample:") ;
	  for (pass = 0 ; pass < 2 ; pass++)
	    for (rw = arrp (rws, 0, RW), ii = 0 ; ii < arrayMax (rws) ; ii++, rw++) 
	      {
		rc = arrp (gx->runs, rw->run, RC) ;
		if (! rc->run || ! rc->aa || rc->private) continue ;
		if ((pass == 0 && rc->runs) || (pass == 1 && ! rc->runs))
		  continue ;
		aceOutf (ao, "\t%s",  rc->sample ? stackText (gx->info, rc->sample) : "") ;
	      }
	  
	  if (gx->hasRunSortingTitle)
	    {
	      aceOutf (ao, "\n##\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tSorting title:") ;
	      for (pass = 0 ; pass < 2 ; pass++)
		for (rw = arrp (rws, 0, RW), ii = 0 ; ii < arrayMax (rws) ; ii++, rw++) 
		  {
		    rc = arrp (gx->runs, rw->run, RC) ;
		    if (! rc->run || ! rc->aa || rc->private) continue ;
		    if ((pass == 0 && rc->runs) || (pass == 1 && ! rc->runs))
		      continue ;
		    aceOutf (ao, "\t%s",  rc->sample ? stackText (gx->info, rc->sortingTitle) : "") ;
		  }
	    }
	  if (gx->hasRunSortingTitle2)
	    {
	      aceOutf (ao, "\n##\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tSorting title 2:") ;
	      for (pass = 0 ; pass < 2 ; pass++)
		for (rw = arrp (rws, 0, RW), ii = 0 ; ii < arrayMax (rws) ; ii++, rw++) 
		  {
		    rc = arrp (gx->runs, rw->run, RC) ;
		    if (! rc->run || ! rc->aa || rc->private) continue ;
		    if ((pass == 0 && rc->runs) || (pass == 1 && ! rc->runs))
		      continue ;
		    aceOutf (ao, "\t%s",  rc->sample ? stackText (gx->info, rc->sortingTitle2) : "") ;
		  }
	    }
	  if (gx->hasRunOtherTitle)
	    {
	      aceOutf (ao, "\n##\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tOther title:") ;
	      for (pass = 0 ; pass < 2 ; pass++)
		for (rw = arrp (rws, 0, RW), ii = 0 ; ii < arrayMax (rws) ; ii++, rw++) 
		  {
		    rc = arrp (gx->runs, rw->run, RC) ;
		    if (! rc->run || ! rc->aa || rc->private) continue ;
		    if ((pass == 0 && rc->runs) || (pass == 1 && ! rc->runs))
		      continue ;
		    aceOutf (ao, "\t%s",  rc->sample ? stackText (gx->info, rc->otherTitle) : "") ;
		  }
	    }
	  	  
	  if (run1 && run2)
	    {
	      RC *rc1 = arrayp (gx->runs, run1, RC) ;
	      RC *rc2 = arrayp (gx->runs, run2, RC) ;
	      KEYSET ks1a = 0, ks1 = rc1 ? rc1->runs : 0 ;
	      KEYSET ks2a = 0, ks2 = rc2 ? rc2->runs : 0 ;
	      
	      if (ks1 && ks2)
		{
		  keySetExclude (&ks1a, &ks2a, ks1, ks2, h) ;
		  aceOutf (ao, "\n##\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tTraining_phenotype:"); 
		  for (pass = 0 ; pass < 2 ; pass++)
		    for (rw = arrp (rws, 0, RW), ii = 0 ; ii < arrayMax (rws) ; ii++, rw++) 
		      {
			int tt = 0 ;
			
			if (ks1a && keySetFind (ks1a, rw->run, 0)) tt = -1 ;
			if (ks2a && keySetFind (ks2a, rw->run, 0)) tt = 1 ;
			rc = arrp (gx->runs, rw->run, RC) ;
			if (! rc->run || ! rc->aa) continue ;
			if ((pass == 0 && rc->runs) || (pass == 1 && ! rc->runs))
			  continue ;
			aceOutf (ao, "\t%d", tt) ;
		      }
		}
	    }
	  
	  if (run3 && run4)
	    {
	      RC *rc3 = arrp (gx->runs, run3, RC) ;
	      RC *rc4 = arrp (gx->runs, run4, RC) ;
	      KEYSET ks3a = 0, ks3 = rc3 ? rc3->runs : 0 ;
	      KEYSET ks4a = 0, ks4 = rc4 ? rc4->runs : 0 ;
	      
	      if (ks3 && ks4)
		{
		  keySetExclude (&ks3a, &ks4a, ks3, ks4, h) ;
		  
		  aceOutf (ao, "\n##\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tValidation_phenotype:"
			   , dictName (gx->runDict, run3)
			   , dictName (gx->runDict, run4)
			   ); 
		  for (pass = 0 ; pass < 2 ; pass++)
		    for (rw = arrp (rws, 0, RW), ii = 0 ; ii < arrayMax (rws) ; ii++, rw++) 
		      {
			int tt = 0 ;
			
			if (ks3 && keySetFind (ks3, rw->run, 0)) tt = -1 ;
			if (ks4 && keySetFind (ks4, rw->run, 0)) tt = 1 ;
			rc = arrp (gx->runs, rw->run, RC) ;
			if (! rc->run || ! rc->aa || rc->private) continue ;
			if ((pass == 0 && rc->runs) || (pass == 1 && ! rc->runs))
			  continue ;
			aceOutf (ao, "\t%d", tt) ;
		      }
		}
	    }
	  
	  aceOutf (ao, "\n##\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tAlpha score"); 
	  for (pass = 0 ; pass < 2 ; pass++)
	    for (rw = arrp (rws, 0, RW), ii = 0 ; ii < arrayMax (rws) ; ii++, rw++) 
	      {
		rc = arrp (gx->runs, rw->run, RC) ;
		if (! rc->run || ! rc->aa || rc->private) continue ;
		if ((pass == 0 && rc->runs) || (pass == 1 && ! rc->runs))
		  continue ;
		aceOutf (ao, "\t%.2f", rw->alpha) ;
	      }
	  
	  for (betaPass = 0 ; betaPass < 2 ; betaPass++)
	    {
	      aceOutf (ao, "\n##\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%s:"
		       , betaPass ? "TP/P" : "TN/N"
		       ); 
	      for (pass = 0 ; pass < 2 ; pass++)
		for (rw = arrp (rws, 0, RW), ii = 0 ; ii < arrayMax (rws) ; ii++, rw++) 
		  {
		    rc = arrp (gx->runs, rw->run, RC) ;
		    if (! rc->run || ! rc->aa || rc->private) continue ;
		    if ((pass == 0 && rc->runs) || (pass == 1 && ! rc->runs))
		      continue ;
		    aceOutf (ao, "\t%.2f", 100.00 * (betaPass ? rw->beta2 : rw->beta1) ) ;
		  }
	    }
	}    
      
      aceOutf (ao, "\n#Species\tComparison\t\tOverexpressed in\tOrdinal\tDifferential element\tNoise score in %d random resamplings\tMagic differential score, noise substracted\tGini coefficient\t Wilcoxon Mann Whitney p-value", compare->nVirtualStrata); 
      aceOutf (ao, "\tLog2(fold change)\tHigh Average index\tLow Average index") ;    
      aceOutf (ao, "\tGene\tGeneId\tGene type\tGene descriptor") ; 
	  
      aceOutf (ao, "\tCoordinates%s %sChromosome strand : from to"
	       , gx->referenceGenome ? " on " : ""
	       , gx->referenceGenome ? gx->referenceGenome : ""
	       ) ;
      if (histo)
	{
	  int i ;

	  /*
	    const char *ccp = "Gene" ;	  
	    if (gx->isTranscript) ccp = "Transcript" ;
	    if (gx->isINTRON) ccp = "Intron" ;
	    if (gx->isSNP) ccp = "SNV" ;
	  */

	  aceOutf (ao, "\tHalf score\tzone start-end\tpeak start-end\tDifferential element") ;

	  for (i = 0 ; i <= 260 ; i += 5)
	    aceOutf (ao, "\t%.1f", i/10.0) ;
	  aceOutf (ao, "\tTotal\tNumber of runs") ;
	}
      else
	 aceOutf (ao, "\tDifferential element") ;

      int nMax = gx->maxGenePlus ;
      for (jj = jpp = 0 ; genePM && jj < keySetMax (genePM) && jj < nMax ; jj++)
	{
	  double hh0s[256], hh1s[256] ;
	  int nRun0 = 0, nRun1 = 0 ;
	  RC *rc1 = arrp (gx->runs, run1, RC) ;
	  RC *rc2 = arrp (gx->runs, run2, RC) ;
	  GC *gc2 ;
	  double pValue = -1, Gini = 0 ;
	  
	  /* int chi2 (int a1, int a2, int b1, int b2, float *c2p) */
	  gene = keySet (genePM, jj) ;	
	  gc = arrp (gx->genes, gene, GC) ;  
	  if (! gene || ! rc1->runs || ! rc2->runs) 
	    continue ;
	  if (wantGeneGroup && ! gc->geneGroup)
	    continue ;
	  if (! wantGeneGroup && gc->geneGroup)
	    continue ;
	  gc2 = gc ;
	  if (gx->hasFromGene)
	    {
	      const char *gNam = gc->fromGene ? stackText (gx->info, gc->fromGene) : 0 ;
	      if (gNam)
		{
		  int fg = 0 ;
		  if (dictFind (gx->geneDict, gNam, &fg))
		    gc2 = arrp (gx->genes, fg, GC) ;
		}
	    }
	  Gini = 0 ; pValue = -1 ;
	  gxMakeComparativeHistos (gx, gene, run1, run2, rc1->runs, rc2->runs, hh0s, hh1s
				   , 0, &nRun0, &nRun1, 0, 0, 0
				   , 0 /* , &threshold */
				   , &Gini,  &pValue
				   , FALSE
				   ) ;

	  if (1 && gene == 50331)
	    {
	      double u = 0, v = 0 ;
	      for(int i = 0 ; i < 256 ; i++)
		{
		  u += hh0s[i] ;
		  v += hh1s[i] ;
		}
	      fprintf (stderr, "###=====  %s histo area %.2g %.2g\n", dictName (gx->geneDict, gene), u, v) ;
	    }

	  aceOutf (ao, "\n%s\t%s\t%s\t%d", species,  dictName (gx->compareDict, compare->compare), vtxtPtr (compareWhat), ++jpp) ;
	  aceOutf (ao, "\t%s", noX__ (dictName (gx->geneDict, gene), !gx->isINTRON)) ;
	  {
	    Array scores ;
	    COMPARE *cmp = gxPairCmp (gx, run1, run2) ;
	    int isPm = (run1 == cmp->vRun1 ? 1 : 2) ;
	    float noise = 0, score = 0 ;

	    if (genePass)
	      isPm = 3 - isPm ;
	    scores= cmp->scores [isPm] ;
	    score = array (scores, gene, float) ;
	    noise = array (cmp->noise, gene, float) ;
	    aceOutf (ao, "\t%.1f\t%.1f", noise, score) ;
	  
	    if (1)
	      {
		float iIndex = rc1->aa ? (arr (rc1->aa, gene, DC)).index : 0 ;
		float jIndex = rc1->aa ? (arr (rc2->aa, gene, DC)).index : 0 ;
		gxDoRegisterComparedGenes (gx, gene, compare->compare, run1, run2, iIndex -1000, jIndex -1000, score, pValue) ; 
	      }
	  }
	  aceOutf (ao, "\t%.1f\t%.1g", Gini, pValue) ;
	  aceOutf (ao, "\t%.2f\t%s%.2f\t%s%.2f"
		   , gc->av1 - gc->av2
		   , gx->isLowTitle[(int)gc->isLow1], gc->av1 - 1000
		   , gx->isLowTitle[(int)gc->isLow2], gc->av2 - 1000
		   ) ;

	  if (gx->hasFromGene)
	    {
	      const char *gNam = gc->fromGene ? stackText (gx->info, gc->fromGene) : 0 ;
	      aceOutf (ao, "\t%s", gNam ? gNam : "") ;
	    }
	  else
	    aceOutf (ao, "\t%s", noX__ (dictName (gx->geneDict, gene), !gx->isINTRON)) ;

	  if (gc2->geneGroup)
	    aceOutf (ao, "\tGroup") ;
	  else
	    aceOutf (ao, "\t%s", gc2->geneId ? stackText (gx->info, gc2->geneId) : "") ;

	  aceOutf (ao, "\t%s", gc2->gene_type ? stackText (gx->info, gc2->gene_type) : "") ;
	  aceOutf (ao, "\t%s", gc2->title ? stackText (gx->info, gc2->title) : "") ;
	  if (gc->chrom)
	    aceOutf (ao, "\t%s%c:%d-%d", dictName (gx->chromDict, gc->chrom), gc->strand, gc->a1, gc->a2) ;
	  else
	    aceOutf (ao, "\t") ;

	  if (histo)
	    {

	      if (genePass == 0)
		{
		  aceOutf (ao, "\t%.2f\t%d-%d\t%d-%d"
			   , gc->score1, gc->za1, gc->za2, gc->ix1, gc->ix2
			   ) ;
		  aceOutf (ao, "\t%s:%d/.", noX__ (dictName (gx->geneDict, gene), !gx->isINTRON), (int)(gc->score1 + .49)) ;
		  
		  gxShowGeneHisto (ao, gx, hh1s, nRun1) ;
		  
		  aceOutf (ao, "\n\t\t\t\t%d", jpp) ;
		  aceOutf (ao, "\t\t\t\t\t\t\t\t") ;

		  if (gx->hasFromGene)
		    {
		      const char *gNam = gc->fromGene ? stackText (gx->info, gc->fromGene) : 0 ;
		      aceOutf (ao, "\t%s", gNam ? gNam : "") ;
		    }
		  else
		    aceOutf (ao, "\t%s", noX__ (dictName (gx->geneDict, gene), !gx->isINTRON)) ;

		  aceOutf (ao, "\t\t\t\t\t%.2f\t%d-%d\t%d-%d"
			   ,  gc->score2, gc->zb1, gc->zb2, gc->iy1, gc->iy2
			   ) ;
		  aceOutf (ao, "\t%s:./%d", noX__ (dictName (gx->geneDict, gene), !gx->isINTRON), (int)(gc->score2 + .49)) ;
		  
		  gxShowGeneHisto (ao, gx, hh0s, nRun0) ;
		}
	      else
		{
		  aceOutf (ao, "\t%.2f\t%d-%d\t%d-%d"
			   , gc->score2, gc->zb1, gc->zb2, gc->iy1, gc->iy2
			   ) ;
		  aceOutf (ao, "\t%s:./%d", noX__ (dictName (gx->geneDict, gene), !gx->isINTRON), (int)(gc->score2 + .49)) ;
		  
		  gxShowGeneHisto (ao, gx, hh1s, nRun1) ;
		  aceOutf (ao, "\n\t\t\t\t%d", jpp) ;
		  aceOutf (ao, "\t\t\t\t\t\t\t\t") ;

		  if (gx->hasFromGene)
		    {
		      const char *gNam = gc->fromGene ? stackText (gx->info, gc->fromGene) : 0 ;
		      aceOutf (ao, "\t%s", gNam ? gNam : "") ;
		    }
		  else
		    aceOutf (ao, "\t%s", noX__ (dictName (gx->geneDict, gene), !gx->isINTRON)) ;

		  aceOutf (ao, "\t\t\t\t\t%.2f\t%d-%d\t%d-%d"
			   ,  gc->score1, gc->za1, gc->za2, gc->ix1, gc->ix2
			   ) ;
		  aceOutf (ao, "\t%s:./%d", noX__ (dictName (gx->geneDict, gene), !gx->isINTRON), (int)(gc->score1 + .49)) ;
		  
		  gxShowGeneHisto (ao, gx, hh0s, nRun0) ;
		}
	    }
	  else
	    {
	      aceOutf (ao, "\t%s", noX__ (dictName (gx->geneDict, gene), !gx->isINTRON)) ;
	      for (pass = 0 ; pass < 2 ; pass++)
		for (rw = arrp (rws, 0, RW), ii = 0 ; ii < arrayMax (rws) ; ii++, rw++) 
		  {
		    rc = arrp (gx->runs, rw->run, RC) ;
		    if (! rc->run || ! rc->aa || rc->private) continue ;
		    if ((pass == 0 && rc->runs) || (pass == 1 && ! rc->runs))
		      continue ;
		    
		    dc = arrp (rc->aa, gene, DC) ;
		    z =  dc->index > 200 ? dc->index + rc->fineTune - 1000 : dc->index - 1000 ;
		    if (z > 0) aceOutf (ao, "\t%.2f", z) ;
		    else aceOut (ao, "\t") ;
		  }
	    }
	}
      aceOutf (ao, "\n\n") ;

    }

  aceOutf (ao, "\n\n") ;		     

  ac_free (h) ;
  return 0 ;
} /* gxExportComparison */

/*************************************************************************************/
/*************************************************************************************/
/* for each gene, if compare contains N groups, show N histo, limiting to minFoldChange between the most extreme groups */
static void gxOneShowAllHistos (GX *gx, int iCompare, int gMax)
{
  AC_HANDLE h = ac_new_handle () ;
  int gene, run ;
  RC *rc ;
  GC *gc ;
  DC *dc ;
  COMPARE *compare = arrp(gx->compares, iCompare, COMPARE) ;
  int geneMax = arrayMax (gx->genes) ;

  ACEOUT ao = aceOutCreate (hprintf (h, "%s.%s", gx->outFileName, dictName (gx->compareDict, iCompare)), ".allHistos.txt", FALSE, h) ;

  aceOutDate (ao, "###", gx->title) ;
  aceOutf (ao, "## All expression histograms for %s, for each gene  minFoldChange %f between the most extreme groups\n",  dictName (gx->compareDict, iCompare), gx->minFoldChange) ;
  aceOutf (ao, "# Gene\tGroup\tTitle") ;
  for (gene = 0 ; gene <= 260 ; gene += 5) 
    aceOutf (ao, "\t%.1f", gene/10.0) ;
  aceOut (ao, "\n") ;
  for (gene = 1, gc = arrayp (gx->genes, gene, GC) ; gene < geneMax ; gene++, gc++) 
    {
      /* check that the gene is expressed somewhere and varies */
      float minIndex = 999999, maxIndex = -999999 ;

      if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
	continue ;
      for (run = 0 ; run < keySetMax (compare->runs) ; run++)
	{
	   rc = arrayp (gx->runs, run, RC) ;
	   if (! rc->aa)
	     continue ;
	   dc = arrayp (rc->aa, gene, DC) ;
	   if (minIndex > dc->index && dc->index > 100) minIndex = dc->index ;
	   if (maxIndex < dc->index && dc->index > 100) maxIndex = dc->index ;
	}
      if (maxIndex < minIndex + gx->minFoldChange/2)
	continue ;
      
      /* export each expression histogram */
      for (run = 0 ; run < keySetMax (compare->runs) ; run++)
	{
	  double hh1s[256], hh2s[256] ; /*  threshold = 0, pValue = 1 ; */
	  int run1, run2, nRun1, nRun2 ;
	  RC *rc1, *rc2 ;

	  run1 = run2 = keySet (compare->runs, run) ;
	  rc1 = rc2 = arrayp (gx->runs, run1, RC) ;
	  if (! rc1->aa)
	    continue ;
	  if (run + 1 < keySetMax (compare->runs)) 
	    {
	      run2 =  keySet (compare->runs, run + 1) ;
	      rc2 = arrayp (gx->runs, run2, RC) ;
	      if (rc2->aa)
		run++ ;   /*  kill 2 birds with one stone */
	      else
		{ 
		  run2 = run1 ; rc2 = rc1 ;
		}
	    }

	  gxMakeComparativeHistos (gx, gene, run1, run2, rc1->runs, rc2->runs, hh1s, hh2s
				   , 0, &nRun1, &nRun2, 0, 0, 0
				   , 0, 0, 0 /* , &threshold, &Gini, &pValue */
				   , FALSE
				   ) ;
	  aceOutf (ao, "%s\t%s", noX__ (dictName (gx->geneDict, gene), !gx->isINTRON), dictName (gx->runDict, run1)) ;
	  aceOutf (ao, "\t%s:%s" , noX__ (dictName (gx->geneDict, gene), !gx->isINTRON), rc1->title ?  stackText (gx->info, rc1->title) : dictName (gx->runDict, run1)) ;
	  gxShowGeneHisto (ao, gx, hh1s, nRun1) ;
	  aceOut (ao, "\n") ;
	  
	  if (1 && run2 != run1) /* kill 2 birds with one stone */
	    {
	      aceOutf (ao, "%s\t%s", noX__ (dictName (gx->geneDict, gene), !gx->isINTRON), dictName (gx->runDict, run2)) ;
	      aceOutf (ao, "\t%s:%s" ,  noX__ (dictName (gx->geneDict, gene), !gx->isINTRON), rc2->title ?  stackText (gx->info, rc2->title) : dictName (gx->runDict, run2)) ;
	      gxShowGeneHisto (ao, gx, hh2s, nRun2) ;
	      aceOut (ao, "\n") ;
	    }
	}
    }
  ac_free (h) ;
  return ;
}  /* gxOneShowAllHistos */

/*************************************************************************************/

static int gxShowAllHistos (GX *gx)
{
  int iCompare = 0 ;
  COMPARE *compare ;

   for (iCompare = 0 ; iCompare < arrayMax (gx->compares) ; iCompare++)
     {
       compare = arrp(gx->compares, iCompare, COMPARE) ;
       if (gx->isINTRON && ! compare->compareINTRON)
	 continue ;
       if (gx->isSNP && ! compare->compareSNP)
	 continue ;
       if (compare->showAllHistos && compare->runs && keySetMax(compare->runs))
	 gxOneShowAllHistos (gx, iCompare, keySetMax(compare->runs)) ;
     }
   return 0 ;
} /* gxShowAllHistos */

/*************************************************************************************/
/*************************************************************************************/
/* Export the info on all diff genes */
static int gxExportDiffGenes (GX *gx, COMPARE *compare, Array pp, int run1, int run2, KEYSET vA, KEYSET vB, int isPm)
{
  static AC_HANDLE h = 0 ;
  static ACEOUT ao = 0 ;
  static vTXT compareWhat = 0 ;
  int i, ipp, jpp, gene, pass ;
  PC *pc ;
  GC *gc ;
  const char *species =  gx->htmlSpecies ? gx->htmlSpecies : "species" ;

  if (isPm == FALSE || !h || !ao)
    {
      if (h) { ac_free (h) ; h = 0 ; ao = 0 ; }
      h = ac_new_handle () ; 
      compareWhat = vtxtHandleCreate (h) ;
      fprintf (stderr, "=== compare %s fdr=%.2f thr=%d\n"
	       , dictName (gx->compareDict, compare->compare)
	       , compare->fdr[isPm]
	       , compare->fdrThreshold[isPm]
	       ) ;
      if (1)
	{
	  ao = aceOutCreate (gx->outFileName
			     , hprintf (h, ".%s.compare.txt"
					, dictName (gx->compareDict, compare->compare)
					)
			     , FALSE, h
			     ) ;
	  aceOutDate (ao, "###", gx->title) ;
	  aceOutf (ao, "## Parameters:  Minimal score selected by optimizing the FDR on %d randomized trials, at FDR = %.2f%% and threshold = %d. %s  %s R= %.2f"
		   ,  compare->nVirtualStrata 
		   , 100 * compare->fdr[isPm] , compare->fdrThreshold[isPm] 
		   , gx->normalize ? "The area of each histogram is normalized to 100" : "The area of each histogram represents the actual number of samples"
		   , gx->geneSelectionMethod == 7 ? "" : messprintf ("geneSelectionMethod=%d",  gx->geneSelectionMethod)
		   , gx->ratio_bound
		   ) ;
	  
	  vtxtPrintf (compareWhat, "overexpressed in\t%s relative to %s using %d versus %d measures"
		      , dictName (gx->runDict, run1)
		      , dictName (gx->runDict, run2)
		      , keySetMax(vA)
		      , keySetMax(vB)
		      ) ;
	  aceOutf (ao, "\n##Comparison: %s :: (private keys run1=%d, run2=%d) "
		   , vtxtPtr (compareWhat)
		   , run1, run2
		   ) ;
	  
	  aceOutf (ao, "\n#Species\tComparison\t\tOverexpressed in\tOrdinal\tDifferential elementl\tNoise score in %d random resamplings\tMagic differential score, noise substracted\t%% FDR\tGini coefficient\t Wilcoxon Mann Whitney p-value", compare->nVirtualStrata); 
	  aceOutf (ao, "\tLog2(fold change)\tHigh Average index\tLow Average index") ;    
	  aceOutf (ao, "\tGene\tGeneId\tGene type\tGene descriptor") ; 
	  aceOutf (ao, "\tCoordinates%s %sChromosome strand : from to"
		   , gx->referenceGenome ? " on " : ""
		   , gx->referenceGenome ? gx->referenceGenome : ""
		   ) ;
	  
	  aceOutf (ao, "\tHalf score\tZone start-end\tpeak start-end") ;
	  aceOutf (ao, "\tDifferential element") ;
	  
	  for (i = 0 ; i <= 260 ; i += 5)
	    aceOutf (ao, "\t%.1f", i/10.0) ;
	  aceOutf (ao, "\tTotal\tRuns") ;
	}
    }
  else
    aceOutf (ao, "\n## Parameters:  Minimal score selected by optimizing the FDR on %d randomized trials, %% at FDR = %.2f and3 threshold = %d"
	     ,  compare->nVirtualStrata 
	     , 100 * compare->fdr[isPm] , compare->fdrThreshold[isPm] 
	     ) ;

    
  vtxtClear (compareWhat) ;
  vtxtPrintf (compareWhat, "overexpressed in\t%s relative to %s using %d versus %d runs"
	      , dictName (gx->runDict, run1)
	      , dictName (gx->runDict, run2)
	      , keySetMax(vA)
	      , keySetMax(vB)
	      ) ;
  
  for (ipp = jpp = 0 ; ipp < arrayMax (pp) ; ipp++)
    {
      double Gini = 0, hh0s[256], hh1s[256], pValue = 1 ;
      int nRun0 = 0, nRun1 = 0 ;
      pc = arrp (pp, ipp, PC) ;
      if (pc->ok < 1) continue ;
      gene = pc->gene ;
      gc = arrp (gx->genes, gene, GC) ;  
      if (gx->targeted && ! gc->targeted)
	continue ;
      if (1)
	if (gx->captured && ! gc->captured)
	  continue ;
		  
      if (! gene) 
	continue ;

      float score = pc->sig ;
      if (! pc->ok || score <= compare->fdrThreshold[isPm]) 
	continue ;
     
      if (1 && ! strcmp (dictName (gx->geneDict, gene), "22__36651046_36653128"))
	fprintf (stderr, "%s %s %s %g\n", dictName (gx->geneDict, gene), dictName (gx->runDict, run1), dictName (gx->runDict, run2), gc->score) ;
      gxPairScoreRegister (gx, compare->compare, run1, run2, gene, score, pValue) ;
      
      for (pass = 0 ; pass < 2 ; pass++)
	{
	  if (ao)
	    {
	      if (pass == 0)
		aceOutf (ao, "\n%s\t%s\t%s\t%d.1", species,  dictName (gx->compareDict, compare->compare), vtxtPtr (compareWhat), ++jpp) ;
	      else
		aceOutf (ao, "\n%s\t%s\t%s\t%d.2", species,  dictName (gx->compareDict, compare->compare), vtxtPtr (compareWhat), jpp) ;
	      
	      aceOutf (ao, "\t%s", noX__(dictName (gx->geneDict, gene), !gx->isINTRON)) ;
	    }
	  if (ao)
	    {
	      int s =  gc->score ;
	      float x = 0 ; 
	      
	      /* Average score in the  NVIRTUALSTRATA random resamplings */ 
	      aceOutf (ao, "\t%.2f\t%.2f", -10, pc->sig) ;

	      if (s > 200) s = 200 ;
	      if (s < -200) s = -200 ;

	      if (s > 0 && s <= 200) 
		x = 100 * compare->fdr[isPm] ;   /* fdr[s] ;  TOTOTOTO */
	      if (x < 1) 
		aceOutf (ao, "\t< 1") ;
	      else
		aceOutf (ao, "\t%.2f", x) ;
	      
	      aceOutf (ao, "\t%.2f\t%.2g", Gini, pValue) ;
	      
	      if (pass == 0)
		aceOutf (ao, "\t%.2f\t%s%.2f\t%s%.2f"
			 , gc->av1 - gc->av2
			 , gx->isLowTitle[(int)gc->isLow1], gc->av1 - 1000
			 , gx->isLowTitle[(int)gc->isLow2], gc->av2 - 1000
			 ) ;
	      else
		aceOutf (ao, "\t\t\t") ;
	      
	      /* this gene column is mandatory and useful if merging the GENE and MRNAH and INTRON and SNP and REARRENGEMENTS  tables, to group the gene and all its transcripts, its introns etc */
	      if (pass >= 0)
		{
		  GC *gc2 = gc ;
		  
		  if (gx->hasFromGene)
		    {
		      const char *gNam = gc->fromGene ? stackText (gx->info, gc->fromGene) : 0 ;
		      if (gNam)
			{
			  int fg = 0 ;
			  if (dictFind (gx->geneDict, gNam, &fg))
			    gc2 = arrp (gx->genes, fg, GC) ;
			  aceOutf (ao, "\t%s", gNam) ;
			}
		      else
			aceOutf (ao, "\t") ;
		    }
		  else
		    aceOutf (ao, "\t%s", noX__(dictName (gx->geneDict, gene), ! gx->isINTRON)) ;
		  
		  if (gc2->geneGroup)
		    aceOutf (ao, "\tGroup") ;
		  else
		    aceOutf (ao, "\t%s", gc2->geneId ? stackText (gx->info, gc2->geneId) : "") ;
		  aceOutf (ao, "\t%s", gc2->gene_type ? stackText (gx->info, gc2->gene_type) : "") ;
		  aceOutf (ao, "\t%s", gc2->title ? stackText (gx->info, gc2->title) : "") ;
		  if (gc->chrom)
		    aceOutf (ao, "\t%s%c:%d-%d", dictName (gx->chromDict, gc->chrom), gc->strand, gc->a1, gc->a2) ;
		  else
		    aceOutf (ao, "\t") ;
		}
	      else
		{
		  aceOutf (ao, "\t\t\t\t\t") ;
		}
	      
	      if (1)
		{
		  if ((pass == 0 && keySet(compare->runs, 0) == run1) || (pass == 1 && keySet(compare->runs, 0) != run1))
		    {
		      aceOutf (ao, "\t%.2f\t%d-%d\t%d-%d"
			       , gc->score1, gc->za1, gc->za2, gc->ix1, gc->ix2
			       ) ;
		      
		      aceOutf (ao, "\t%s:%d/.", noX__(dictName (gx->geneDict, gene), !gx->isINTRON), (int)(gc->score1 + .49)) ;
		      gxShowGeneHisto (ao, gx, hh0s, nRun0) ;
		      
		    }
		  else
		    {
		      aceOutf (ao, "\t%.2f\t%d-%d\t%d-%d"
			       , gc->score2, gc->zb1, gc->zb2, gc->iy1, gc->iy2
			       ) ;
		      
		      aceOutf (ao, "\t%s:%d/.",  noX__ (dictName (gx->geneDict, gene), !gx->isINTRON), (int)(gc->score2 + .49)) ;
		      
		      gxShowGeneHisto (ao, gx, hh1s, nRun1) ;
		    }
		} 
	    }
	}
    }

  if (ao)
    aceOutf (ao, "\n\n\n\n") ;
  
  if (isPm == 2)
    {
      ac_free (h) ;
      h = 0 ; ao = 0 ; 
    }

  return 0 ;
} /* gxExportDiffGenes */

/*************************************************************************************/

static void gxExportSignatures (GX *gx, int run1, int run2, int run3, int run4, int pass)
{
  static ACEOUT ao = 0 ;
  static AC_HANDLE h = 0 ;
  RC *rc, *rc1, *rc2 ;
  int run ;

  return ; /* 2019_05_31  : obsolete */
  if (pass == -2) /* close */
    { aceOutf (ao, "\n") ; ao = 0 ; ac_free (h) ; return ; }

  if (pass == -1) /* open */
    {
      h = ac_new_handle () ;
      ao = aceOutCreate (gx->outFileName, ".newSampleClassificationBySignature.txt", FALSE, h) ;
      aceOutDate (ao, "###", gx->title) ;
      gxExportTableHeader (gx, ao, 32) ;
      return ;
    }

  rc1 = arrayp (gx->runs, run1, RC) ;
  rc2 = arrayp (gx->runs, run2, RC) ;

  aceOutf (ao, "\nRatio:%d %s/%s", pass, dictName (gx->runDict, run1), dictName (gx->runDict, run2)) ;
  if (0 && run3 && run4 && run3 != run1)
    aceOutf (ao, "//%s/%s", dictName (gx->runDict, run3), dictName (gx->runDict, run4)) ;
  aceOutf (ao, "//AUC:%.1f MCC:%.1f", gx->AUC, gx->MCC) ;
  if (0) aceOutf (ao, " /AUC2:%.1f MCC2:%.1f ",gx->AUC2, gx->MCC2) ;

  if (gx->method) aceOutf (ao, "\t%s", gx->method) ;

  if (gx->hasRunId) 
    aceOutf (ao, "\t%s/%s"
	     , rc1->runid ? stackText (gx->info, rc1->runid) : dictName (gx->runDict, run1)
	     , rc2->runid ? stackText (gx->info, rc2->runid) : dictName (gx->runDict, run2)
	     ) ;
  if (gx->hasRunSample) 
    aceOutf (ao, "\t%s/%s"
	     , rc1->sample ? stackText (gx->info, rc1->sample) : dictName (gx->runDict, run1)
	     , rc2->sample ? stackText (gx->info, rc2->sample) : dictName (gx->runDict, run2)
	     ) ;
  if (gx->hasRunTitle) 
    aceOutf (ao, "\t%s/%s"
	     , rc1->title ? stackText (gx->info, rc1->title) : dictName (gx->runDict, run1)
	     , rc2->title ? stackText (gx->info, rc2->title) : dictName (gx->runDict, run2)
	     ) ;


  if (1)
    {
      aceOutf (ao, "\nAlpha score:%d %s/%s", pass, dictName (gx->runDict, run1), dictName (gx->runDict, run2)) ;
      for (run = 1, rc = arrayp (gx->runs, run, RC) ; run < arrayMax (gx->runs) ; rc++, run++) 
	aceOutf (ao, "\t%.2f", rc->alpha) ;	
    }
  if (1)
    {
      aceOutf (ao, "\nTN/N percentage:%d %s/%s", pass, dictName (gx->runDict, run1), dictName (gx->runDict, run2)) ;
      for (run = 1, rc = arrayp (gx->runs, run, RC) ; run < arrayMax (gx->runs) ; rc++, run++) 
	aceOutf (ao, "\t%.2f", 100.0 * rc->beta1) ;
      aceOutf (ao, "\nTP/P percentage:%d %s/%s", pass, dictName (gx->runDict, run1), dictName (gx->runDict, run2)) ;
      for (run = 1, rc = arrayp (gx->runs, run, RC) ; run < arrayMax (gx->runs) ; rc++, run++) 
	aceOutf (ao, "\t%.2f", 100.0 * rc->beta2) ;
    }
  aceOutf (ao, "\n") ;
  return ;
} /* gxExportSignatures */

/*************************************************************************************/
/* compute alpha for 'run' by comparing the 'genes'  to 'run1' and 'run2' 
 * run is any run or group and we compute the correl (run1, run2)
 */
/* compute the scalar product of run1 and run2 limited to genes
 */
static BOOL gxRunsCosineUsingGenes (GX *gx, KEYSET genes, Associator assPlus, Associator assMinus, Array sp, Array sm, int run1, int run2, double *alphap, int *nnp, int lowAll, BOOL isDiff) 
{
  RC *rc1, *rc2 ;
  GC *gc ;
  DC *dc1, *dc2 ;
  Array aa1, aa2 ;
  int low ;
  int geneMax, gene ;
  double minLevel = lowAll ;
  int N ;
  double x, y, z, zs, w, X, Y, X2, Y2, XY ;
  BOOL ok ;
  void *vp = 0, *vm = 0 ;
  int gp = 0, gm = 0 ;
  SCAL *scalp ;
  int doubleLowSnp = 0 ;
  int doubleHighSnp = 0 ;
  rc1 = arrayp (gx->runs, run1, RC) ;
  rc2 = arrayp (gx->runs, run2, RC) ;

  rc1->run = run1 ;
  aa1 = rc1->aa ;
  if (! aa1) 
    return FALSE ;
  
  aa2 = rc2->aa ;
  if (! aa2) 
    return FALSE ;

  N = 0 ; X = X2 = Y = Y2 = XY = 0 ;  
  if (run1 == run2)
    { N = arrayMax (aa1) ; goto done ; }
  minLevel = rc1->NA_crossover_index ;
  if (minLevel < rc2->NA_crossover_index)
    minLevel = rc2->NA_crossover_index ;
  
  if (minLevel < lowAll)
    minLevel = lowAll ;

  geneMax = arrayMax (aa1) ;
  if (geneMax > arrayMax (aa2))
    geneMax = arrayMax (aa2) ;
  gc = arrayp (gx->genes, geneMax - 1, GC) ; /* make room */

  for (gene = 1, gc = arrayp (gx->genes, gene, GC), dc1 = arrayp (aa1, gene, DC), dc2 = arrayp (aa2, gene, DC)  ;
       gene < geneMax ; gc++, dc1++, dc2++, gene++)
    {
      ok = FALSE ;
      if (gx->targeted && ! gc->targeted)
	continue ;
      if (1)
	if (gx->captured && ! gc->captured)
	  continue ;
      if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
	continue ;

      if (! gc->isGood)
	continue ;
		  
      if (! genes && gc->maxIndex - gc->minIndex < 1)
	continue ;
      if (dc1->isLow && dc2->isLow)  /* 2013_08_31 ATTENTION, this may alter the correlation and the compare_to functions */
	continue ;
      ok = FALSE ;
      vp = vm = 0 ;
      if ( gx->digital == 1 ||
	   gx->digital == 4 ||
	  (genes && keySetFind (genes, gene, 0)) ||  
	  (assPlus && assFind (assPlus, assVoid(gene), &vp)) ||
	  (assMinus && assFind (assMinus, assVoid(gene), &vm))
	  )
	{
	  float seuil ;
	  gp = assInt (vp) ; gm =assInt (vm) ; /* ordinal of gene in plus and minus sets */
	  x = dc1->index + rc1->fineTune - 1000 ;
	  y = dc2->index + rc2->fineTune - 1000 ;
	  switch (gx->digital)
	    {
	    case 0: /* systeme classique, produit scalaire des index  used by gxCorrelation */
	      low = 0 ;
	      if (x < minLevel) { low++ ; x = minLevel ; }
	      if (y < minLevel) { low++ ; y = minLevel ; }
	      if (low < 2)
		ok = TRUE ;
	      if (isDiff)
		{ x -= gc->averageIndex ; y-= gc->averageIndex ; }
	      break ; 

	    case 5: /* produit scalaire avec seuil index  used by gxCorrelation */
	      /* much less  good than digital == 1 for NB gene prediction */
	      low = 0 ;
	      if (x < minLevel) { low++ ; x = minLevel ; }
	      if (y < minLevel) { low++ ; y = minLevel ; }
	      if (low < 2)
		ok = TRUE ;
	      seuil = 0.8 ;
	      if (seuil < 2 * gc->intrinsicSigma) seuil =  gc->intrinsicSigma ;
	      if (seuil < 2 * gx->minFoldChange) seuil = gx->minFoldChange ;
	      if (x < y + seuil && x > y - seuil) { x = y = (x + y)/2 ; }
	      if (isDiff)
		{ x -= gc->averageIndex ; y-= gc->averageIndex ; }
	      break ; 

	    case 1: /* systeme weighted digital */
	      z = 0 ; x *= 10 ;
	      zs = gc->score > 40 ? gc->score - 40 : 0 ; /* 2013_02_12 */
	      zs = zs * zs ;   /* dromadaire */
	      if (x >= gc->za1 && x <= gc->za2 ) z = - zs ;
	      if (x >= gc->zb1 && x <= gc->zb2 ) z =   zs ;
	      x = z ;
	      z = 0 ; y *= 10 ;
	      if (y >= gc->za1 && y <= gc->za2 ) z = - zs ;
	      if (y >= gc->zb1 && y <= gc->zb2 ) z =   zs ;
	      y = z ;
	      ok = TRUE ;
	      break ;

	    case 2: /* systeme fully digital */
	      z = 0 ; x *= 10 ;
	      zs = gc->score > 40 ? gc->score - 40 : 0 ; /* 2013_02_12 */
	      zs = zs * zs ;
	      if (x >= gc->za1 && x <= gc->za2 ) z = - 1 ;
	      if (x >= gc->zb1 && x <= gc->zb2 ) z =   1 ;
	      x = z ;
	      z = 0 ; y *= 10 ;
	      if (y >= gc->za1 && y <= gc->za2 ) z = - 1 ;
	      if (y >= gc->zb1 && y <= gc->zb2 ) z =   1 ;
	      y = z ;
	      ok = TRUE ;
	      break ;

	    case 3: /* Mobymatic, 2015_03, not yet tested as of 2016_06 */
	      z = 0 ; x *= 10 ;
	      /*
	      {
		float  a, b, c ;
		a = 100 * hhs[0, x] ; b = 100 * hhs[1, x] ; c = (a + b) / 100 ; 
		if (c > 0) 
		  { a /= c ; b /= c ; }
		if (a > b + 40) z = - (a - b) ;
		if (b > a + 40) z = + (b - a) ;
	      }
	      x = z/100 ;
	      z = 0 ; y *= 10 ;
	      {
		float  a, b, c ;
		a = 100 * hhs[0, y] ; b = 100 * hhs[1, y] ; c = (a + b) / 100 ; 
		if (c > 0) 
		  { a /= c ; b /= c ; }
		if (a > b + 40) z = - (a - b) ;
		if (b > a + 40) z = + (b - a) ;
	      }
	      y = z/100 ;
	      ok = TRUE ;
	      */
	      break ;
	    case 4: /* SNP who_is_who */
	      low = 0 ;
	      if (x < 0 || y < 0)
		break ;
	      if (x < 4) { low++ ; x = 0 ; }
	      else if (x > 18) { x = 20 ; }

	      if (y < 4) { low++ ; y = 0 ; }
	      else if (y > 18) { y = 20 ; }

	      if (0 && gc->averageIndex > 10) /* beware or genomic errors */
		{ x = 20 - x ; y = 20 - y ; }
	      if (low < 2)
		ok = TRUE ;
	      else
		doubleLowSnp++ ;
	      if (x > 18 && y > 18)
		doubleHighSnp++ ;
	      if (isDiff)
		{ x -= gc->averageIndex ; y -= gc->averageIndex ; }
	      break ; 

	    case 4444: /* this version seems bugged SNP, count by classes removing sampling errors : 2016_06_07 */
	      if (x >= 0 && y >= 0)  /* both runs should have measured the SNP */
		{
		  low = 0 ;
		  if (x < 4) { low++ ; }
		  if (y < 4) { low++ ; }
		  if (low < 2000)   
		    ok = TRUE ;   /* wild/wild does not contribute to the scalar product */
		  if (ok)   /* 0: using allelic frequencies gives correl aroung 70% */
		    {      /* 1: using classes gives correl around 96% */ 
		      if (x > 18) x = 20 ;
		      else if (x < 2) x = 0 ;
		      else x = 10 ;
		      if (y > 18) y = 20 ;
		      else if (y < 2) y = 0 ;
		      else y = 10 ;
		      if (0)
			{
			  ok = FALSE ; /* discard all hetero, only keep 1 or 2 home */ 
			  if (x == 20 || y == 20) ok = TRUE ;
			  if (x == 10 || y == 10) ok = FALSE ;
			}
		    }
		  if (1)
		    { /* 2016_06 */
		      if (isDiff)
			{ x -= gc->averageIndex ; y-= gc->averageIndex ; }
		    }
		  else
		    { /* 2016_11_24 */
		      if (0 && gc->averageIndex > 10) /* beware or genomic errors */
			{ x = 20 - x ; y = 20 - y ; }
		      if (x + y == 0) 
			ok = FALSE ;
		      else
			ok = TRUE ;
		      if (x == 0)
			{
			  if (y ==  0) { x =   0 ; y =   0 ; }
			  if (y == 10) { x = -10 ; y =  10 ; }
			  if (y == 20) { x = -14 ; y =  14 ; }
			}
		      else if (x == 10)
			{
			  if (y ==  0) { x =  10 ; y = -10 ; }
			  if (y == 10) { x =  14 ; y =  14 ; }
			  if (y == 20) { x = -10 ; y =  10 ; }
			}
		      else if (x == 20)
			{
			  if (y ==  0) { x =  14 ; y = -14 ; }
			  if (y == 10) { x =  10 ; y = -10 ; }
			  if (y == 20) { x =  14 ; y =  14 ; }
			}
		    }
		  break ;
		}
	      break ;
	    }
	}
      if (ok)
	{
	  N++ ; X += x ; X2 += x * x ; Y += y ; Y2 += y * y ; XY += x*y ;
	  if (gp && sp)
	    { 
	      scalp = arrayp (sp, gp, SCAL) ;
	      scalp->N = 1 ;
	      scalp->X = x ;
	      scalp->Y = y ;
	    }
	  if (gm && sm)
	    { 
	      scalp = arrayp (sm, gm, SCAL) ;
	      scalp->N = 1 ;
	      scalp->X = x ;
	      scalp->Y = y ;
	    }
	}
    }
  x = y = w = 0 ;
  if (0)
    fprintf (stderr, "#########  doubleHighSnp = %d  run1 %s  run2 %s\n"
	    , doubleHighSnp
	    , dictName (gx->runDict, run1) 
	    , dictName (gx->runDict, run2) 
	    ) ;
  if (0)
    fprintf (stderr, "#########  doubleLowSnp = %d\n",  doubleLowSnp) ;
  if (doubleHighSnp)
    {
      if (doubleHighSnp > doubleLowSnp)
	doubleHighSnp = doubleLowSnp ;
      N += doubleHighSnp ; /* we are effectively adding lots of points at zero zero 
			    * to insure that a few points at 1,0 and 0,1 do not create
			    * a negative correlation 
			    */ 

    }
 done:
  
  if (N)
    {  /* computig N*N may overflow the int range */
      x = X/N ; x = x * x ; x = X2/N - x ; 
      y = Y/N ; y = y * y ; y = Y2/N - y ;
      w = X/N ; w *= Y/N ; w = XY/N - w ;
      if (x * y > 0)
	w = w / sqrt (x*y) ;
      else
	w = 1 ;
 
      if (0)
	fprintf (stderr, "#########  N=%d, X=%.2g Y=%.2g X2=%.2g Y2=%.2g XY=%.2g x=%.2g y=%.2g w=%.2g\n",N, X,Y,X2,Y2,XY,x,y,w) ;
      *alphap = w ;
      if (nnp) 
	*nnp = N ;
      return TRUE ;
    }
  *alphap = -3 ;
  return FALSE ;
} /* gxRunsCosineUsingGenes */

/*************************************************************************************/
/* compute for all runs R  cosinus (R,run1) - cosinus (R,run2) */
static Array gxAllRunsCosineUsingGenesRun1Run2 (GX *gx, int run1, int run2, KEYSET genePlus, KEYSET geneMinus, Array sp, Array sm, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  int jj, run, nn = 0 ;
  double w1 = 0, w2 = 0 ;
  RC *rc = 0 ;
  RC *rc1 = arrayp (gx->runs, run1, RC) ;
  RC *rc2 = arrayp (gx->runs, run2, RC) ;
  Array rws = 0 ;
  RW *rw ;
  int runMax = arrayMax (gx->runs) ;
  int nnn[] = {1000000, 0, 2,10,20,40,80,160,320,600,0} ;
  Associator assPlus = assHandleCreate (h) ;
  Associator assMinus = assHandleCreate (h) ;

  for (jj = 0 ; nnn[jj] ; jj++)
    {
      keySetMerge (genePlus, geneMinus, assPlus, assMinus) ;   /* will sort compress */
      /* rws = array */
      if (rc1->aa && rc2->aa && !rc1->isSubLib && ! rc2->isSubLib)
	{
	  rws = arrayHandleCreate (arrayMax (gx->runs), RW, h0) ;
	  for (run = 1, rc = arrayp (gx->runs, run, RC) ; run < runMax - 1 ; rc++, run++)
	    {
	      if (! rc->aa) continue ;
	      gxRunsCosineUsingGenes (gx, 0, assPlus, assMinus, sp, sm, run, run1, &w1, 0, 0, FALSE) ;
	      gxRunsCosineUsingGenes (gx, 0, assPlus, assMinus, sp, sm, run, run2, &w2, 0, 0, FALSE) ;
	      
	      rw = arrayp (rws, nn++, RW) ;
	      rw->run = run ;
	      rc->alpha = rw->alpha = 100 * (w1 - w2) ;
	    }	  
	  arraySort (rws, rwOrder) ;  
	}
    }
  ac_free (h) ;
  return rws ;
} /* gxAllRunsCosineUsingGenesRun1Run2 */

/*************************************************************************************/

static void gxOneCorrelationExport (GX *gx, int iCompare, int lowAll, int NR1
				   , Array aa, KEYSET runs1, KEYSET runs2, KEYSET rParent, BOOL isDiff)
{
  AC_HANDLE h = ac_new_handle () ;
  int irun1, irun2, run1, run2 ;
  RC *rc1, *rc2 ;
  float w ;
  COMPARE *compare = arrp (gx->compares, iCompare, COMPARE) ;
  ACEOUT ao = aceOutCreate (gx->outFileName, messprintf ("%s.Correlation.%s.MinIndex%d.txt"
							 , isDiff ? ".Differential" : ""
							 , dictName (gx->compareDict, iCompare)
							 , lowAll)
			    , FALSE, h) ;
  aceOutDate (ao, "###", gx->title) ;
  aceOutf (ao, "## %s.Correlation.%s.MinIndex%d.txt\n"
	   , gx->outFileName ? gx->outFileName : "stdout"
	   , dictName (gx->compareDict, iCompare)
	   , lowAll) ;
  aceOutf (ao, "## The lines are in the order imposed by the user in Compare->runs:2%s.\n"
	   , compare->autosort ? ", the columns are sorted automatically" : "the columns are ordered by their sorting_title"
	   ) ;
  
  {
    /* export the header reordered according to runs2 */
    Array runsOk = gx->runs ;
    
    gx->runs = arrayCreate (arrayMax(gx->runs), RC) ;
    for (irun2 = 0 ; irun2 < keySetMax(runs2) ; irun2++)
      {
	run2 =  keySet (runs2, irun2) ;
	rc1 = arrp (runsOk, run2, RC) ; rc1->run = run2 ;
	rc2 = arrayp (gx->runs, irun2+1, RC) ;
	*rc2 = *rc1 ;
      }
    if (1) gxExportTableHeader (gx, ao, 52) ;
    arrayDestroy (gx->runs) ;
    gx->runs = runsOk ;
  }
  
  if (rParent)
    {
      KEYSET closestKs = keySetHandleCreate (h) ;
      int parent, closest = 0 ;
      aceOutf (ao, "\n#Parent") ;
      for (irun2 = 0 ; irun2 < keySetMax (runs2) ; irun2++)
	{
	  run2 = keySet (runs2, irun2) ;
	  parent = keySet (rParent, run2) ;
	  aceOutf (ao, "\t%s", parent ? dictName (gx->runDict, parent) : "") ;
	}
      aceOutf (ao, "\n#Closest") ;
      for (irun2 = 0 ; irun2 < keySetMax (runs2) ; irun2++)
	{
	  float bestW = -100 ;
	  run2 = keySet (runs2, irun2) ;
	  for (irun1 = 0 ; irun1 < keySetMax (runs1) ; irun1++)
	    {
	      run1 = keySet (runs1, irun1) ;
	      w = array(aa, run2 * NR1 + run1, float) ;
	      if (run1 != run2 && w > bestW)
		{
		  bestW = w ;
		  keySet (closestKs, run2) = run1 ;
		}
	    }
	  closest = keySet (closestKs, run2) ;
	  aceOutf (ao, "\t") ;
	  if (closest)
	    aceOutf (ao, "%s", dictName (gx->runDict, closest)) ;
	}
      aceOutf (ao, "\n#Score to closest") ;
      for (irun2 = 0 ; irun2 < keySetMax (runs2) ; irun2++)
	{
	  run2 = keySet (runs2, irun2) ;
	  closest = keySet (closestKs, run2) ;
	  aceOutf (ao, "\t") ;
	  if (closest)
	    {
	      w = array(aa, run2 * NR1 + closest, float) ;
	      aceOutf (ao, "%.2f", w) ;
	    }
	}
      aceOutf (ao, "\n#Score to parent") ;
      for (irun2 = 0 ; irun2 < keySetMax (runs2) ; irun2++)
	{
	  run2 = keySet (runs2, irun2) ;
	  parent = keySet (rParent, run2) ; 
	  if (parent)
	    {
	      w = array(aa, run2 * NR1 + parent, float) ;
	      aceOutf (ao, "\t%.2f", w) ;
	    }
	  else 
	    aceOut (ao, "\t") ;
	}
 
   }

  aceOutf (ao, "\n# Title") ;
  for (irun2 = 0 ; irun2 < keySetMax (runs2) ; irun2++)
    {
      run2 = keySet (runs2, irun2) ;
      rc2 = arrayp (gx->runs, run2, RC) ;
      if (rc2->title)
	aceOutf (ao, "\t%s", stackText (gx->info, rc2->title)) ;
      else
	aceOutf (ao, "\t%s", dictName (gx->runDict, rc2->run)) ;
    }

  for (irun1 = 0 ; irun1 < keySetMax (runs1) ; irun1++)
    {
      run1 = keySet (runs1, irun1) ;
      rc1 = arrayp (gx->runs, run1, RC) ;
      
      aceOutf (ao, "\n%s", dictName (gx->runDict, run1)) ;

      for (irun2 = 0 ; irun2 < keySetMax (runs2) ; irun2++)
	{
	  run2 = keySet (runs2, irun2) ;
	  w = array(aa, run2 * NR1 + run1, float) ;
	  if (0) printf ("YYYYY\t%d\t%s\t%d\t%s\t%.2f\n", run1, dictName(gx->runDict, run1),  run2, dictName(gx->runDict, run2), w) ;
	  if (1 || w > 0)
	    aceOutf (ao, "\t%.2f", w) ;
	  else
	    aceOutf (ao, "\t") ;
	}
      if (rc1->title)
	aceOutf (ao, "\t%s", stackText (gx->info, rc1->title)) ;
      else
	aceOutf (ao, "\t%s", dictName (gx->runDict, rc1->run)) ;
    }
  
  aceOut (ao, "\n") ;
  
  ac_free (h) ;
  return ;
} /*  gxOneCorrelationExport */

/*************************************************************************************/
/* try to sort the table using the chronoligical ordering */
typedef struct csStruct { int irun ; int run ; int nn ; float weight ; unsigned char word[1024] ; } CS ;

static int csOrder (const void *a, const void *b)
{
  const CS *up = (const CS*)a ;
  const CS *vp = (const CS*)b ;
  
  return  up->weight - vp->weight ; 
}

static int csWordOrder (const void *a, const void *b)
{
  const CS *up = (const CS*)a ;
  const CS *vp = (const CS*)b ;
  
  return  strcmp ((char *)up->word, (char *)vp->word) ;
}

static void gxOneCorrelationSortWeightColumnsOld (int NR1, Array cc1, Array cc3, Array kk1, KEYSET runs1, KEYSET runs2)
{
  int irun1, irun2, run1, run2 ;
  float w ;
  CS *cs ;
  
  kk1 = arrayReCreate (kk1, keySetMax(runs1), CS) ;
  for (irun1 = 0 ; irun1 < keySetMax (runs1) ; irun1++)
    {
      run1 = keySet (runs1, irun1) ;
      cs = arrayp (kk1, irun1, CS) ;
      cs->irun = irun1 ; cs->run = run1 ; cs->weight = 0 ;
      for (irun2 = 0 ; irun2 < keySetMax (runs2) ; irun2++)
	{
	  run2 = keySet (runs2, irun2) ;
	  if (arrayMax (cc3) > run2 * NR1 + run1)
	    {
	      w = arr (cc3, run2 * NR1 + run1, int) ;
	      if (w > 0) cs->weight += 1000.0 * irun2 ;
	    }
	  if (arrayMax (cc1) > run2 * NR1 + run1)
	    {
	      w = arr (cc1, run2 * NR1 + run1, int) ;
	      if (w > 0) cs->weight += 1 * irun2 ;
	    } 
	}
    }
} /* gxOneCorrelationSortWeightColumnsOld */

/****************/

static void gxOneCorrelationSortWeightLinesOld (int NR1, Array cc1, Array cc3, Array kk2, KEYSET runs1, KEYSET runs2)
{
  int irun1, irun2, run1, run2 ;
  float w ;
  CS *cs ;
  
  kk2 = arrayReCreate (kk2, keySetMax(runs2), CS) ;
  for (irun2 = 0 ; irun2 < keySetMax (runs2) ; irun2++)
    {
      run2 = keySet (runs2, irun2) ;
      cs = arrayp (kk2, irun2, CS) ;
      cs->irun = irun2 ; cs->run = run2 ; cs->weight = 0 ;
      for (irun1 = 0 ; irun1 < keySetMax (runs1) ; irun1++)
	{
	  run1 = keySet (runs1, irun1) ;
	  if (arrayMax (cc3) > run2 * NR1 + run1)
	    {
	      w = arr (cc3, run2 * NR1 + run1, int) ;
	      if (w > 0) cs->weight += 1000.0 * irun1 ;
	    }
	  if (arrayMax (cc1) > run2 * NR1 + run1)
	    {
	      w = arr (cc1, run2 * NR1 + run1, int) ;
	      if (w > 0) cs->weight += 1 * irun1 ;
	    }
	}
    }
} /* gxOneCorrelationSortWeightLinesOld */

/****************/

static void gxOneCorrelationChronoSortOld (GX *gx, int iCompare, int NR1, Array aa, KEYSET runs1, KEYSET runs2, int lowAll)
{
  AC_HANDLE h = ac_new_handle () ;
  Array cc1 = arrayHandleCreate (arrayMax(aa), float, h) ;
  Array cc3 = arrayHandleCreate (arrayMax(aa), float, h) ;
  float w ;
  int i, j, pass ;
  int irun1, irun2, run1, run2 ;
  CS *cs ;
  Array kk1 = 0, kk2 = 0 ;
  
  for (irun2 = 0 ; irun2 < keySetMax (runs2) ; irun2++)
    {
      Array bb = arrayCreate (keySetMax(runs1), CS) ;
      run2 = keySet (runs2, irun2) ; 
      for (j = irun1 = 0 ; irun1 < keySetMax (runs1) ; irun1++)
	{
	  run1 = keySet (runs1, irun1) ;
	  w = arr (aa, run2 * NR1 + run1, float) ;
	  if (w > 0)
	    {
	      cs = arrayp (bb, j++, CS) ;
	      cs->run = run1 ;
	      cs->weight = - 1000.0 * w ;
	    }
	}
      arraySort (bb, csOrder) ;
      for (i = 0 ; i<j && i < 3 ; i++)
	{
	  cs = arrp (bb, i, CS) ;
	  if (cs->weight < 0)
	    {
	      if (i < 1)
		array (cc1, run2 * NR1 + cs->run, float) = 1 ;
	      if (i < 3)
		array (cc3, run2 * NR1 + cs->run, float) = 1 ;
	    }
	}
      arrayDestroy (bb) ;
    }
  
  /* chronological ordering applied after thresholding */
  for (pass = 0 ; pass < 100 ; pass++)
    {
      int ok = 0 ;
      /* weight columns */
      kk1 = arrayReCreate (kk1, keySetMax(runs1), CS) ;
      gxOneCorrelationSortWeightColumnsOld (NR1, cc1, cc3, kk1, runs1, runs2) ;
      /* order columns */
      arraySort (kk1, csOrder) ;
      for (i = 0, cs = arrp(kk1, 0, CS) ; i < arrayMax(runs1) ;cs++,  i++)
	{
	  if (keySet (runs1, i) != cs->run) 
	    ok++ ;
	  keySet (runs1, i) = cs->run ;
	}
      
      /* weight lines */
      kk2 = arrayReCreate (kk2, keySetMax(runs2), CS) ;
      gxOneCorrelationSortWeightLinesOld (NR1, cc1, cc3, kk2, runs1, runs2) ;
      /* order lines */
      arraySort (kk2, csOrder) ;
      for (i = 0, cs = arrp(kk2, 0, CS) ; i < arrayMax(runs2) ;cs++,  i++)
	{
	  if (keySet (runs2, i) != cs->run) 
	    ok++ ;
	  keySet (runs2, i) = cs->run ;
	}
      if (0) 
	fprintf (stderr, "---  gxOneCorrelationSort %s pass %d minIndex %d  flips %d : %s\n"
		 , dictName (gx->compareDict, iCompare)
		 , pass, lowAll, ok
		 , timeShowNow() ) ;

      if (ok == 0) break ;
    }
  ac_free (h) ;
  return ;
} /* gxOneCorrelationChronoSortOld */

/****************/

static void gxOneCorrelationTopoChronoSort (GX *gx, int iCompare, int NR1, Array aa, KEYSET runs1, KEYSET runs2, int lowAll)
{
  AC_HANDLE h = ac_new_handle () ;
  int i, iMax = keySetMax (runs1) ;
  int j, jMax = keySetMax (runs2) ;
  KEYSET lines = keySetHandleCopy (runs1, h) ;
  KEYSET cols = keySetHandleCopy (runs2, h) ;
  KEYSET ks ;
  Array bb = arrayHandleCreate (iMax * jMax, double, h) ; 
  array (bb, iMax * jMax - 1, double) = 0 ;

  for (i = 0 ; i < iMax ; i++)
    {
      int run1 = keySet (runs1, i) ;
      for (j = 0 ; j < jMax ; j++)
	{
	  int run2 = keySet (runs2, j) ;
	  float w = arr (aa, run2 * NR1 + run1, float) ;
	  if (w > 20)
	    arr (bb, i * jMax + j, double) = w ;
	}
    }
  topoChronoOrder (bb, lines, cols) ;

  /* reorder runs1 as a function of cols which contains indices, not runs */
  if (1)
    {
      ks = keySetHandleCopy (runs1, h) ;
      for (i = 0 ; i < iMax ; i++)
	keySet (runs1, i) = keySet (ks, keySet (lines, i)) ;
    }
  if (1)
    {
      ks = keySetHandleCopy (runs2, h) ;
      for (j = 0 ; j < jMax ; j++)
	keySet (runs2, j) = keySet (ks, keySet (cols, j)) ;
    }

  ac_free (h) ;
  return ;
} /* gxOneCorrelationTopoChronoSort */

/****************/
/* sort strictly in each group defined by beeig closest to that line, in score decreasing for that line */
#ifdef JUNK

static void gxOneCorrelationSort (GX *gx, int iCompare, int NR1, Array aa, KEYSET runs1, KEYSET runs2)
{
  AC_HANDLE h = ac_new_handle () ;
  Array bb1 = arrayHandleCreate (keySetMax(runs1), CS, h) ;
  Array bb2 = arrayHandleCreate (keySetMax(runs2), CS, h) ;
  float w ;
  int i, j ;
  int irun1, irun2, run1, run2 ;
  char *word ;
  CS *cs ;
  KEYSET ks ;

  for (irun2 = 0 ; irun2 < keySetMax (runs2) ; irun2++)
    {
      cs = arrayp (bb2, irun2, CS) ;
      cs->run = irun2 ;
      bb1 = arrayCreate (keySetMax(runs1), CS) ;
      run2 = keySet (runs2, irun2) ; 
      for (j = irun1 = 0 ; irun1 < keySetMax (runs1) ; irun1++)
	{
	  run1 = keySet (runs1, irun1) ; /* runs1 are the groups (colored lines of the cavariance table) */
	  w = arr (aa, run2 * NR1 + run1, float) ;
	  if (w > 0)
	    {
	      cs = arrayp (bb1, j++, CS) ;
	      cs->run = irun1 ;
	      cs->weight = - 1000.0 * w ;
	    }
	}
      arraySort (bb1, csOrder) ;
      for (i = 0 ; i<j && i < 64 ; i++)
	{
	  cs = arrp (bb1, i, CS) ;
	  irun1 = cs->run ;
	  word = array(bb2, irun2, CS).word ;
	  
	  if (i == 0)
	    {   /* 2014_08)24: intercale le score du premier avant la liste ordinale des autres */
	      float c1 = - cs->weight / 1000 ;
	      int c ;
	      if (c1 < 0) c1 = 0 ;
	      if (c1 > 100) c1 = 100 ;
	      c = c1 ;
	      word[0] = irun1 < 254 ? irun1 + 1 : 255 ;
	      word[1] = 101 - c ;
	    }
	  else
	    word[i+1] = irun1 < 254 ? irun1 + 1 : 255 ;
	}
      arrayDestroy (bb1) ;
    }
  arraySort (bb2, csWordOrder) ;
  ks = keySetCopy (runs2) ;
  for (i = 0, cs = arrp (bb2, i, CS) ; i < keySetMax (runs2) ; i++, cs++)
    {
      int run =  keySet (ks, cs->run) ; 
      keySet (runs2, i) = run ;
      if (0)
	{
	  RC *rc = arrayp (gx->runs,  run, RC) ;

	  fprintf (stderr, "\n---- %s", dictName(gx->runDict, run)) ;
	  fprintf (stderr, "\t %s", rc->title ?  stackText (gx->info, rc->title) : dictName(gx->runDict, run)) ;
	  for (j = 0 ; j < 8 ; j++)
	    {
	      word = cs->word ;
	      if (word[j] >= 0)
		fprintf (stderr, "\t%d", word[j]-1) ;
	      else
		break ;
	    }
	}
    }

  keySetDestroy (ks) ;

  ac_free (h) ;
  return ;
} /* gxOneCorrelationSort */
#endif

/****************/
/* sort inside each group by a combination of weigths */
static void gxOneCorrelationSort2 (GX *gx, int iCompare, int NR1, Array aa, KEYSET runs1, KEYSET runs2)
{
  AC_HANDLE h = ac_new_handle () ;
  Array bb1 = arrayHandleCreate (keySetMax(runs1), CS, h) ;
  Array bb2 = arrayHandleCreate (keySetMax(runs2), CS, h) ;
  float w ;
  int i, j, pass ;
  int irun1, irun2, run1, run2 ;
  unsigned char *word ;
  CS *cs ;
  KEYSET ks ;

  for (irun2 = 0 ; irun2 < keySetMax (runs2) ; irun2++)
    {
      cs = arrayp (bb2, irun2, CS) ;
      cs->run = irun2 ;
      bb1 = arrayCreate (keySetMax(runs1), CS) ;
      run2 = keySet (runs2, irun2) ; 
      for (j = irun1 = 0 ; irun1 < keySetMax (runs1) ; irun1++)
	{
	  run1 = keySet (runs1, irun1) ; /* runs1 are the groups (colored lines of the cavariance table) */
	  w = arr (aa, run2 * NR1 + run1, float) ;
	  if (w > 0)
	    {
	      cs = arrayp (bb1, j++, CS) ;
	      cs->run = irun1 ;
	      cs->weight = - 1000.0 * w ;
	    }
	}
      arraySort (bb1, csOrder) ;
      for (i = 0 ; i<j && i < 64 ; i++)
	{
	  cs = arrp (bb1, i, CS) ;
	  irun1 = cs->run ;
	  word = array(bb2, irun2, CS).word ;
	  
	  if (i == 0)
	    {   /* 2014_08)24: intercale le score du premier avant la liste ordinale des autres */
	      float c1 = - cs->weight / 1000 ;
	      int c ;
	      int j1 ;
	      float damper = 1 ;
	      if (j >= 1) c1 = 0 ;
	      for (j1 = 1 ; j1 < j ; j++)
		{
		  damper *= 1.3 ;
		  c1 -=  (cs+1)->weight / (1000 * damper) ;
		}
	      if (c1 < 0) c1 = 0 ;
	      if (c1 > 100) c1 = 100 ;
	      c = c1 ;
	      word[0] = irun1 < 254 ? irun1 + 1 : 255 ;
	      word[1] = 101 - c ;
	    }
	  else
	    word[i+1] = irun1 < 254 ? irun1 + 1 : 255 ;
	}
      arrayDestroy (bb1) ;
    }
  arraySort (bb2, csWordOrder) ;
  ks = keySetCopy (runs2) ;
  for (i = 0, cs = arrp (bb2, i, CS) ; i < keySetMax (runs2) ; i++, cs++)
    {
      int run =  keySet (ks, cs->run) ; 
      keySet (runs2, i) = run ;
      if (0)
	{
	  RC *rc = arrayp (gx->runs,  run, RC) ;

	  fprintf (stderr, "\n---- %s", dictName(gx->runDict, run)) ;
	  fprintf (stderr, "\t %s", rc->title ?  stackText (gx->info, rc->title) : dictName(gx->runDict, run)) ;
	  for (j = 0 ; j < 8 ; j++)
	    {
	      word = cs->word ;
	      if (word[j] >= 0)
		fprintf (stderr, "\t%d", word[j]-1) ;
	      else
		break ;
	    }
	}
    }
  keySetDestroy (ks) ;

  /* last phase: sort in a frist block the runs that reach .60
   *  then the remainder
   */
  ks = keySetCopy (runs2) ;
  for (i = 0, pass = 0 ; pass < 2 ; pass++)
    {
      int irun1, irun2 ;

      for (irun2 = 0 ; irun2 < keySetMax (ks) ; irun2++)
	{
	  float bestW = -100 ;
	  run2 = keySet (ks, irun2) ;
	  for (irun1 = 0 ; irun1 < keySetMax (runs1) ; irun1++)
	    {
	      run1 = keySet (runs1, irun1) ;
	      w = array(aa, run2 * NR1 + run1, float) ;
	      if (w > bestW)
		{
		  bestW = w ;
		}
	    }
	  if (
	      (pass == 0 && bestW >= .6) ||
	      (pass == 1 && bestW < .6)
	      )
	    keySet (runs2, i++) = run2 ;
	}
    }
  keySetDestroy (ks) ;

  ac_free (h) ;
  return ;
} /* gxOneCorrelationSort2 */

/*************************************************************************************/

static void gxOneCorrelation (GX *gx, int iCompare, int lowAll, BOOL isDiff)
{
  AC_HANDLE h = ac_new_handle () ;
  Array aa = 0 ;
  KEYSET runs1, runs2, rParent = 0 ;
  int irun, n, gene, run1, run2, low ;
  RC *rc1, *rc2 ;
  GC *gc ;
  COMPARE *compare = arrp(gx->compares, iCompare, COMPARE) ;
  KEYSET genes = keySetHandleCreate (h) ;
  double w ;
  int NR1 ;
  int geneMax = arrayMax (gx->genes) ;
  int oldDigital = gx->digital ;
  
  gx->digital = gx->digitalCorrelation ;  /* dromadaire 2016_06_08 : prodit scalair avec seuil */
  
  runs1 = keySetHandleCreate (h) ;    /* the groups (compare->runs) or colored lines, use to classify the runs */
  runs2 = keySetHandleCreate (h) ;    /* the runs we want to classify */
  for (n = gene = 0, gc = arrayp (gx->genes, gene, GC) ; gene < geneMax ; gc++, gene++)
    {
      if (gc->maxIndex - gc->minIndex < 2)
	continue ;
      if (gc->maxIndex < 8 || gc->maxIndex < lowAll)
	continue ;
      keySet (genes, n++) = gene ;
    } 
  keySetSort(genes) ; keySetCompress (genes) ;
  
  /* the list runs1 is the list of runs or groups given in the compare */
  for (n = 0, irun = 0 ; irun < keySetMax (compare->runs) ; irun++)
    {
      run1 = keySet (compare->runs, irun) ;
      rc1 = arrp (gx->runs, run1, RC) ;
      if (! rc1->private)
	{
	  keySet (runs1, n++) = run1 ;
	}
    }

  rParent = keySetHandleCreate (h) ;    /* the parent of the runs we want to classify */
  for (n = irun = 0 ; irun < keySetMax (runs1) ; irun++)
    {
      run1 = keySet (runs1, irun) ; 
      keySet (runs2, n++) = run1 ; 
      rc1 = arrp (gx->runs, run1, RC) ;
      if (rc1->groupLevel > 0)
	keySet (rParent, run1) = run1 ;
    }
  if (compare->show_all_runs)  /* in runs2, we want all runs */
    {
      for (run2 = 0 ; run2 < arrayMax (gx->runs) ; run2++)
	{
	  rc2 = arrp (gx->runs, run2, RC) ;
	  if (rc2->run)
	    if (! rc2->private && ! rc2->isSubLib && run2 != gx->runAny && rc2->groupLevel == 0)
	      {
		if (gx->isSNP || ( rc2->Genes_with_index_over_10 > 4000 &&  rc2->Genes_with_index_over_18 < 800))
		  keySet (runs2, n++) = run2 ;
	      }
	}
    }
  /* do not yet sort-compress runs2, because in the granpParent search, we want to start from runs1 list */
  if (1)  /* we want in runs2 the same runs1 and their filiation */
    { 
      /* copy list 1 */
      /* expand recursivelly */
      for (irun = 0 ; irun < keySetMax (runs2) ; irun++)
	{
	  int j, grandParent ;
	  run1 = keySet (runs2, irun) ;
	  rc1 = arrp (gx->runs, run1, RC) ;
	  grandParent = keySet (rParent, run1) ;
 	  if (rc1->runs && ! rc1->addCounts)
	    for (j = 0 ; j < keySetMax (rc1->runs) ; j++)
	      {
		run2 = keySet (rc1->runs, j) ;
		if (grandParent && !keySet (rParent, run2))
		  keySet (rParent, run2) = grandParent ;
		rc2 = arrp (gx->runs, run2, RC) ;
		if (! rc2->private && ! rc2->isSubLib)
		  { 
		    if (gx->isSNP || ( rc2->Genes_with_index_over_10 > 4000 &&  rc2->Genes_with_index_over_18 < 800))
		      {
			if (! compare->show_all_runs) 
			  keySet (runs2, n++) = run2 ;
		      }
		  }
	      }
	}
    }

  if (1) /* 2016_12_14 : remove the groups from the column */
    {
      for (n = irun = 0 ; irun < keySetMax (runs2) ; irun++)
	{
	  run2 = keySet (runs2, irun) ; 
	  rc2 = arrp (gx->runs, run2, RC) ;
	  if (rc2->groupLevel == 0)
	    keySet (runs2, n++) = run2 ;
	}
    }
  keySetSort (runs2) ; keySetCompress (runs2) ;

  NR1 = keySetMax (gx->runs) ;        /* not relevant, used as place older to construct square matrices */
  aa = arrayHandleCreate (NR1 * (1+dictMax(gx->runDict)), float, h) ;
  array (aa, NR1 * (1+dictMax(gx->runDict) - 1), float) = 0 ; /* make room */
  for (low = 0, irun = 0 ; irun < keySetMax (runs1) ; irun++)
    {
      int j ;
      run1 = keySet (runs1, irun) ;
      rc1 = arrayp (gx->runs, run1, RC) ;

      for (j = 0 ; j < keySetMax (runs2) ; j++)
	{
	  run2 = keySet (runs2, j) ;
	  rc2 = arrp (gx->runs, run2, RC) ;

	  low = rc1->NA_crossover_index > rc2->NA_crossover_index ? rc1->NA_crossover_index : rc2->NA_crossover_index ;
	  if (low < lowAll) low = lowAll ;

	  if (gxRunsCosineUsingGenes (gx, genes, 0,0,0,0, run1, run2, &w, 0, low, isDiff))
	    {
	      array (aa, run2 * NR1 + run1, float) = 100 * w ;
	      if (0) 
		fprintf (stderr, "XXXX\t%d\t%s\t%d\t%s\t%.2f\n", run1, dictName(gx->runDict, run1),  run2, dictName(gx->runDict, run2), w) ;
	    } 
	}
    }
  
  /* do not sort the groups, the user must sort them using the numbers */
  if (compare->autosort)
    {
      if (1) gxOneCorrelationTopoChronoSort (gx, iCompare, NR1, aa, runs1, runs2, lowAll) ; /* sort the groups (and the runs) */
      if (0) gxOneCorrelationChronoSortOld (gx, iCompare, NR1, aa, runs1, runs1, lowAll) ; /* sort the groups (and the runs) */
      if (0) gxOneCorrelationSort2 (gx, iCompare, NR1, aa, runs1, runs2) ;  /* sort the runs in a more simple way */
    }
  else if (keySetMax(runs1) == keySetMax (runs2))
    {
      KEYSET ks1 = keySetCopy (runs1) ;
      KEYSET ks3, ks2 = keySetCopy (runs2) ;
      keySetSort (ks1) ;
      keySetSort (ks2) ;
      ks3 = keySetMINUS (ks1, ks2) ;
      if (keySetMax (ks3) == 0)     /* only happens if there are no groups in runs1, since groups are always removed from runs */ 
	{
	  int i ;
	  for (i = 0 ; i < keySetMax (ks2) ; i++)
	    keySet (runs2, i) = keySet (runs1, i) ;
	}
      keySetDestroy (ks1) ;
      keySetDestroy (ks2) ;
      keySetDestroy (ks3) ;
    }
  gxOneCorrelationExport (gx, iCompare, lowAll, NR1, aa, runs1, runs2, rParent, isDiff) ;

  gx->digital = oldDigital ;  
  ac_free (h) ;
  return ;
} /*  gxOneCorrelation */

/*************************************************************************************/

static void gxCorrelation (GX *gx, int lowAll, BOOL isDiff)
{
  int iCompare = 0 ;
  COMPARE *compare ;

   for (iCompare = 0 ; iCompare < arrayMax (gx->compares) ; iCompare++)
     {
       compare = arrp(gx->compares, iCompare, COMPARE) ;
       if (gx->isINTRON && ! compare->compareINTRON)
	 continue ;
       if (gx->isSNP && ! compare->compareSNP)
	 continue ;

       if (compare->correlation && compare->runs && keySetMax(compare->runs))
	 gxOneCorrelation (gx, iCompare, lowAll, isDiff) ;
     }
} /* gxCorrelation */

/*************************************************************************************/
/*************************************************************************************/

typedef struct gtStruct { int gene, length, intronGene, intronType, intronRefSeq, intronAv, colMax, colMin, geneGroup, capture, intronSharedD, intronSharedA ; float min, max, delta, score, score_1_2, index[100]; int tags[100] ; char shape[100], sigShape[100] ; BOOL isIntron ; } GT ;

static int gtOrder (const void *a, const void *b)
{
  const GT *up = (const GT *)a, *vp = (const GT *)b ;
  float z = up->score - vp->score ;

  if (0)
    {
      if (up->score >= 150 && vp->score < 150) return -1 ;
      if (vp->score >= 150 && up->score < 150) return 1 ;
      
      if (vp->score >= 150)
	{
	  z = up->colMax - vp->colMax ;
	  if (z > 0) return 1 ;
	  if (z < 0) return -1 ;
	  z = up->colMin - vp->colMin ;
	  if (z > 0) return 1 ;
	  if (z < 0) return -1 ;
	}
    }
  if (0)
    {
      z = up->score - vp->score ;
      if (z > 0) return -1 ;
      if (z < 0) return 1 ;
      
      z = up->colMax - vp->colMax ;
      if (z > 0) return 1 ;
      if (z < 0) return -1 ;
      z = up->colMin - vp->colMin ;
      if (z > 0) return 1 ;
      if (z < 0) return -1 ;
    }

  if (1)
    {
      /* signif shape, then score sum, then score max */
      int n ;
      const char *cp = up->sigShape ;
      const char *cq = vp->sigShape ;
      
      n = up->geneGroup - vp->geneGroup ; if (n) return n ;

      if (up->score && ! vp->score)  return -1 ;
      if (! up->score && vp->score)  return 1 ;
      if (cp && !cq) return -1 ;
      if (cq && !cp) return 1 ;


      n = strcmp (cp, cq) ; if (n) return n ;

      z = up->score - vp->score ;
      if (z > 0) return -1 ;
      if (z < 0) return 1 ;
    }
  return 0 ;
} /* gtOrder */

/*************************************************************************************/
/*************************************************************************************/

typedef struct dgStruct { 
  int gene, iCompare, iRun, jRun ; 
  float iIndex, jIndex, diffScore, pValue ; 
} DG ;

static int dgOrder (const void *a, const void *b)
{
  const DG *up = (const DG *)a, *vp = (const DG *)b ;
  int g = up->gene - vp->gene ;
  float z = up->diffScore - vp->diffScore ;

  if (g) return g ;
  return z > 0 ? -1 : 1 ;  /* high value first */
  return 0 ;
} /* dgOrder */

/*************************************************************************************/

static void gxDoRegisterComparedGenes (GX *gx, int gene, int iCmp, int iRun, int jRun, float iIndex, float jIndex, float diffScore, float pValue) 
{
  if (! gx->diffGenes)
    gx->diffGenes = arrayHandleCreate (1000, DG, gx->h) ;
  if (gene) /* register */
    {
      DG *dg = arrayp (gx->diffGenes, arrayMax(gx->diffGenes), DG) ;
      dg->gene = gene ;
      dg->iCompare = iCmp ;
      dg->diffScore = diffScore ;
      dg->pValue = pValue ;
      dg->iRun = iRun ;
      dg->jRun = jRun ;
      dg->iIndex = iIndex ;
      dg->jIndex = jIndex ;
    }
  return ;
} /* gxRegisterComparedGenes  */

/*************************************************************************************/

static void gxRegisterComparedGenes (GX *gx, GT *gt, int iCmp, int iRun, int jRun, float diffScore, float pValue)
{
  if (gt) /* register */
    {
      gxDoRegisterComparedGenes (gx, gt->gene, iCmp, iRun, jRun,  gt->index[iRun], gt->index[jRun], diffScore, pValue) ;
    }
  return ;
} /* gxRegisterComparedGenes  */

/*************************************************************************************/

static void gxExportComparedGenes (GX *gx)
{
  if (gx->diffGenes && arrayMax (gx->diffGenes)) /* export and clean up */
    {
      AC_HANDLE h = ac_new_handle () ;
      int oldGene = 0 ;
      int ii, iMax = arrayMax ( gx->diffGenes) ;
      int jj, jMax = dictMax (gx->compareDict) ;
      DG *dg ;
      ACEOUT ao = aceOutCreate (gx->outFileName, ".differential_genes.ace", gx->gzo, h) ;
      DICT *dict = dictHandleCreate (256, h) ;
      KEYSET ks = keySetHandleCreate (h) ;

      /* we must be sure that we eliminate old compare values in all genes */
     for (jj = 1 ; jj <= jMax ; jj++)
	{  /* protect the compare names once for all genes */
	  dictAdd (dict, ac_protect (dictName (gx->compareDict, jj), h), &ii) ;
	  keySet (ks, jj) = ii ;
	}

      arraySort (gx->diffGenes, dgOrder) ;
      arrayCompress (gx->diffGenes) ;
      for (ii = 0, dg = arrp (gx->diffGenes, 0, DG) ; ii < iMax ; ii++, dg++)
	{
	  if (dg->gene != oldGene)
	    {
	      oldGene = (dg->gene) ;
	      aceOutf (ao, "\nGene %s\n", 
		       ac_protect (noX__ (dictName (gx->geneDict, dg->gene), !gx->isINTRON), h)
		     ) ;
	      for (jj = 1 ; jj <= jMax ; jj++)
		if (keySet (ks, jj))
		  aceOutf (ao, "-D Compare %s\n", 
			   dictName (dict, keySet (ks,jj))
			   ) ;
	    }
	  if (dg->iCompare)
	    aceOutf (ao, "Compare %s %.2f %s %.2f %s %.2f\n"
		     , dictName (dict, keySet (ks, dg->iCompare))
		     , dg->diffScore
		     , dictName (gx->runDict, dg->iRun)
		     , dg->iIndex /* exp(dg->iIndex * log(2.0))/1000.0 */
		     , dictName (gx->runDict, dg->jRun)
		     , dg->jIndex /* exp(dg->jIndex * log(2.0))/1000.0 */
		     ) ;
	}
      aceOut (ao, "\n") ;
      ac_free (h) ;
     }
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, "// gxExportComparedGenes done :  %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
  }
} /* gxExportComparedGenes */

/*************************************************************************************/
/*************************************************************************************/

static void gxShowFpkm (ACEOUT ao, float z)
{
  char *prefix = "" ;
  
  if (z < -8000)
    { z += 10000 ; prefix = "NE:" ; }
  else if (z < -800)
    { z += 1000 ; prefix = "NA:" ; }
  aceOutf (ao, "\t%s%.2f", prefix, exp (z*log(2.0))/1000.0) ;
} /* gxShowFpkm */

/*************************************************************************************/

static void gxOneExpressionProfile (GX *gx, int iCompare, int pass
				    , BOOL showPvalue	
				    , BOOL wantGeneGroup
				    )
{
  AC_HANDLE h = ac_new_handle () ;
  Array gts = 0, aa = 0 ;
  int irun, gene, run, i, line, n1 = 0, n2 = 0, n3 = 0 ;
  int  n_diff_scores, n_diff_scores1, n_diff_scores2 ;
  float sum_score, sum_score1, sum_score2 ;
  RC *rc ;
  GC *gc ;
  DC *dc ;
  COMPARE *compare = arrp(gx->compares, iCompare, COMPARE) ;
  int geneMax = arrayMax (gx->genes) ;
  double threshold = 0 ;
  GT *gt ;
  ACEOUT ao = 0 ;
  int rMax = compare->runs ? keySetMax(compare->runs) : 0 ;
  KEYSET bs = keySetHandleCreate (h) ;
  const char *oldSigShape ;
  double ln2 = log(2.0) ;
  int nBADcapt = 0, nGOOD1capt = 0, nGOOD2capt = 0, nUNIQUEcapt = 0, nSCORE1capt = 0, nSCORE2capt = 0 ;
  float sumGoodScore1capt = 0, sumGoodScore2capt = 0, sumUniqueScorecapt = 0  ;
  int nBAD = 0, nGOOD1 = 0, nGOOD2 = 0, nUNIQUE = 0, nSCORE1 = 0, nSCORE2 = 0, nCAPTURED = 0 ;
  float sumGoodScore1 = 0, sumGoodScore2 = 0, sumUniqueScore = 0 ;

  switch (pass + 3 * (showPvalue ? 1 : 0))
    {
    case 0: /* all profiles */
      ao = aceOutCreate (hprintf (h, "%s.%s", gx->outFileName, dictName (gx->compareDict, iCompare)), wantGeneGroup ? ".score.pathways.profiles.txt" : (! gx->isINTRON ? ".score.genes.profiles.txt" : ".score.introns.profiles.txt"), FALSE, h) ; 
      break ;
    case 1: /* DEG   max > 0 */
      ao = aceOutCreate (hprintf (h, "%s.%s", gx->outFileName, dictName (gx->compareDict, iCompare)), (! gx->isINTRON ? ".score.genes.differential_profiles.txt" : ".score.introns.differential_profiles.txt"), FALSE, h) ; 
      break ;
    case 2: /* TOP DEG max >= 185 */
			 ao = aceOutCreate (hprintf (h, "%s.%s", gx->outFileName, dictName (gx->compareDict, iCompare)), (! gx->isINTRON ? ".score.genes.top_differential_profiles.txt" : ".score.introns.top_differential_profiles.txt"), FALSE, h) ; 
     break ;

    case 3: /* all profiles */
      ao = aceOutCreate (hprintf (h, "%s.%s", gx->outFileName, dictName (gx->compareDict, iCompare))
			 , wantGeneGroup ? ".pValue.pathways.profiles.txt" : (! gx->isINTRON ? ".pValue.genes.profiles.txt" : ".pValue.introns.profiles.txt")
			 , FALSE, h) ; 
      break ;
    case 4: /* DEG   max > 0 */
			 ao = aceOutCreate (hprintf (h, "%s.%s", gx->outFileName, dictName (gx->compareDict, iCompare)), (! gx->isINTRON ? ".pValue.genes.differential_profiles.txt" : ".pValue.introns.differential_profiles.txt"), FALSE, h) ; 
      break ;
    case 5: /* TOP DEG max >= 185 */
			 ao = aceOutCreate (hprintf (h, "%s.%s", gx->outFileName, dictName (gx->compareDict, iCompare)), (! gx->isINTRON ? ".pValue.genes.top_differential_profiles.txt" : ".pValue.introns.top_differential_profiles.txt"), FALSE, h) ; 
     break ;
    }

  if (rMax > 100) rMax = 100 ;
  
  aceOutDate (ao, "###", gx->title) ;
  aceOut (ao, "#Order") ; n1 = 1 ;
  if (gx->isINTRON) { aceOut (ao, "\t#Intron (chrom__from_to)") ; n1++ ;}
  else if (gx->isTranscript) { aceOut (ao, "\t#Transcript") ; n1++ ;}
  else  { aceOutf (ao, "\t#%s", wantGeneGroup ? "Gene-group" : "Gene") ; n1++ ;}
  if (gx->hasGeneLength) { aceOut (ao, "\t#Length") ; n1++; }
  if (gx->hasGeneId) { aceOut (ao, "\t#NCBI GeneId") ; n1++ ; } 
  if (gx->isINTRON && gx->hasFromGene) { aceOut (ao, "\t#Gene") ; n1++ ; } 
  if (gx->hasGeneModel) { aceOut (ao, "\t#Gene") ; n1++ ; } 
  if (gx->hasIntronType) { aceOut (ao, "\t#Intron type") ; n1++ ; } 
  if (gx->hasKnownDonorAcceptor){ aceOut (ao, "\t#New (known donor acceptor)") ; n1++ ; } 
  if (gx->hasIntronShared) { aceOut (ao, "\t#Shared donor\tShared acceptor (New Av RefSeq)") ; n1+=2 ; } 
  if (gx->hasRefSeqAv) { aceOut (ao, "\t#Number of AceView transcripts including this intron") ; n1+=1 ; } 
  if (gx->hasCapture) { aceOut (ao, "\t#Captured by") ; n1++ ; } 
  if (gx->hasCaptureTouch) { aceOut (ao, "\t#Touched by") ; n1++ ; } 
  if (0 && gx->hasGeneNm) { aceOut (ao, "\t#RefSeq transcript Id") ; n1++ ; } 
  if (gx->hasGeneAlias) { aceOut (ao, "\t#Alias") ; n1++ ; } 
  if (gx->hasGeneChrom) { aceOut (ao, "\t#Chromosome\t#Strand\t#from base\t#to base") ; n1 += 4 ; } 
  if (gx->hasGeneType) { aceOut (ao, "\t#Type") ; n1++ ; } 
  if (gx->hasGeneTitle) { aceOut (ao, "\t#Title") ; n1++ ; } 
  
  if (gx->isINTRON) { aceOut (ao, "\t#Intron") ; n1++ ;}
  else if (gx->isTranscript) { aceOut (ao, "\t#Transcript") ; n1++ ;}
  else  { aceOutf (ao, "\t#%s", wantGeneGroup ? "Gene-group" : "Gene") ; n1++ ;}

  for (irun = 0 ; irun < rMax ; irun++)
    aceOutf (ao, "\tExpression index %s", dictName (gx->runDict, keySet (compare->runs, irun))) ;

  aceOut (ao, "\tMaximal expression at\tMinimal expression at\tShape\tSignificant shape") ; n2 = 4 ;
  aceOut (ao, "\tGene") ; n2++ ;

  if (rMax <= 15)
    {
      int jrun ;
      for (irun = 0 ; irun < rMax ; irun++)
	for (jrun = irun + 1 ; jrun < rMax ; jrun++)
	  {
	    aceOutf (ao, "\tDifferential score if %s > %s"
		     , dictName (gx->runDict, keySet (compare->runs, irun))
		     , dictName (gx->runDict, keySet (compare->runs, jrun))
		     ) ;
	    aceOutf (ao, "\tDifferential score if %s < %s"
		     , dictName (gx->runDict, keySet (compare->runs, irun))
		     , dictName (gx->runDict, keySet (compare->runs, jrun))
		     ) ;
	  }

    }
  else
    {
      for (irun = 0 ; irun < rMax - 1 ; irun++)
	{
	  int jrun = irun + 1 ;
	  aceOutf (ao, "\tDifferential score if %s > %s"
		   , dictName (gx->runDict, keySet (compare->runs, irun))
		   , dictName (gx->runDict, keySet (compare->runs, jrun))
		   ) ;
	  aceOutf (ao, "\tDifferential score if %s > %s"
		   , dictName (gx->runDict, keySet (compare->runs, jrun))
		   , dictName (gx->runDict, keySet (compare->runs, irun))
		   ) ;
	}
    }

  aceOutf (ao, "\tMaximal differential score (using consistency across the individual samples defining  a group)") ; n3 = 1 ;
  aceOutf (ao, "\tSum of the differential scores of the odd columns") ; n3++ ;
  aceOutf (ao, "\tSum of the differential scores of the even columns") ; n3++ ;
  aceOutf (ao, "\tSum of the differential scores") ; n3++ ;
  aceOutf (ao, "\tNumber of comparisons with differential scores") ; n3++ ;
  aceOutf (ao, "\tInconsistencies: odd and even") ; n3++ ;

  aceOutf (ao, "\t") ; n3++ ;
  for (irun = 0 ; irun < rMax ; irun++)
    {
      if (gx->isINTRON)
	aceOutf (ao, "\tSupport %s", dictName (gx->runDict, keySet (compare->runs, irun))) ;
      else
	aceOutf (ao, "\tsFPKM %s", dictName (gx->runDict, keySet (compare->runs, irun))) ;
    }
  aceOutf (ao, "\tMax sFPKM\tMin sFPKM\tFold change") ;
   
  aceOutf (ao, "\tIndex Max\tIndex Min\tLog2(fold change)") ;

  if (0)
    {
      aceOut (ao, "\n#Title") ;  
      for (i = 1 ; i < n1 ; i++)
	aceOut (ao, "\t") ;
      for (irun = 0 ; irun < rMax ; irun++)
	{
	  rc = arrp (gx->runs, keySet (compare->runs, irun), RC) ;
	  aceOutf (ao, "\t%s", rc->title ? stackText (gx->info, rc->title) : "" ) ;
	}
      for (i = 0 ; i < n2 + n3 + 6 ; i++)
	aceOut (ao, "\t") ;
    }
  aceOut (ao, "\n#Number of DEGs, maximal FDR and differential score threshold") ;
  
  for (i = 1 ; i < n1 ; i++)
    aceOut (ao, "\t") ;
  for (irun = 0 ; irun < rMax ; irun++)
    aceOutf (ao, "\t") ;
  for (i = 0 ; i < n2 ; i++)
    aceOut (ao, "\t") ;
  if (rMax <= 15)
    {
      int jrun ;
      for (irun = 0 ; irun < rMax ; irun++)
	for (jrun = irun + 1 ; jrun < rMax ; jrun++)
	  {
	    int run1 = keySet (compare->runs, irun) ;
	    int run2 = keySet (compare->runs, jrun) ;
	    COMPARE *cmp = gxPairCmp (gx, run1, run2) ;

	    if (cmp)
	      {
		int dr = (run1 > run2 ? 0 : 1) ;
		aceOutf (ao, "\t%d DEG >= 185, %d DEG > 0, at FDR %.2f%% and threshold %d"
			 , cmp->nDEG185[2 - dr]
			 , cmp->nDEG[2 - dr]
			 , 100 * cmp->fdr[2 - dr]
			 , cmp->fdrThreshold[2 - dr] 
			 ) ;
		aceOutf (ao, "\t%d DEG >= 185, %d DEG > 0, at FDR %.2f%% and threshold %d"
			 , cmp->nDEG185[1 + dr]
			 , cmp->nDEG[1 + dr]
			 , 100 * cmp->fdr[1 + dr]
			 , cmp->fdrThreshold[1 + dr]
			 ) ;
	      }
	    else
	      aceOutf (ao, "\t\t") ;
	  }
    }
  else
    {
      for (irun = 0 ; irun < rMax - 1 ; irun++)
	{
	  int run1 = keySet (compare->runs, irun) ;
	  int run2 = keySet (compare->runs, irun + 1) ;
	  COMPARE *cmp = gxPairCmp (gx, run1, run2) ;
	  if (1)
	    {
	      int dr = (run1 > run2 ? 0 : 1) ;
	      aceOutf (ao, "\t%d DEG >= 185, %d DEG > 0, at FDR %.2f%% and threshold %d"
		       , cmp->nDEG185[2 - dr]
		       , cmp->nDEG[2 - dr]
		       , 100 * cmp->fdr[2 - dr]
		       , cmp->fdrThreshold[2 - dr] 
		       ) ;
	      aceOutf (ao, "\t%d DEG >= 185, %d DEG > 0, at FDR %.2f%% and threshold %d"
		       , cmp->nDEG185[1 + dr]
		       , cmp->nDEG[1 + dr]
		       , 100 * cmp->fdr[1 + dr]
		       , cmp->fdrThreshold[1 + dr]
		       ) ;
	    }
	  else
	    aceOutf (ao, "\t\t") ;
	}
    }

  aceOut (ao, "\n#Group number") ;
  
  for (i = 1 ; i < n1 ; i++)
    aceOut (ao, "\t") ;
  for (irun = 0 ; irun < rMax ; irun++)
    aceOutf (ao, "\t%d", 1 + irun) ;
  for (i = 0 ; i < n2 ; i++)
    aceOut (ao, "\t") ;
  if (rMax <= 15)
    {
      int jrun ;
      for (irun = 0 ; irun < rMax ; irun++)
	for (jrun = irun + 1 ; jrun < rMax ; jrun++)
	  aceOutf (ao, "\t %d>%d\t %d<%d", 1 + irun, 1 + jrun, 1 + irun, 1 + jrun) ;
    }
  else
    {
      for (irun = 0 ; irun < rMax - 1 ; irun++)
	aceOutf (ao, "\t %d>%d\t %d<%d", 1 + irun, 2 + irun, 1 + irun, 2 + irun) ;
    }
  for (i = 0 ; i < n3 ; i++)
    aceOut (ao, "\t") ;
  for (irun = 0 ; irun < rMax ; irun++)
    aceOutf (ao, "\t%d", 1 + irun) ;
 
  aceOut (ao, "\n") ;

  gts = arrayHandleCreate (geneMax + 1, GT, h) ;
  array (gts, geneMax, GT).gene = 0 ; /* make room */ ;
  n1 = 0 ;
  for (i = irun = 0 ; irun < rMax ; irun++)
    {
      run = keySet (compare->runs, irun) ;
      rc = arrp (gx->runs, run, RC) ;
      if (! rc->aa || rc->private)
	continue ;

      if (! i++)
	threshold = rc->NA_crossover_index ;
      if (threshold < rc->NA_crossover_index)
	threshold = rc->NA_crossover_index ;

      aa = rc->aa ;
      for (gene = 1, gc = arrayp (gx->genes, gene, GC), dc = arrayp (aa, gene, DC) ;
       gene < geneMax ; gc++, dc++, gene++)
	{
	  if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
	    continue ;
	  if (! gc->isGood)
	    continue ;
	  if (!gx->isINTRON && gx->targeted && !gc->targeted)
	    continue ;
	  if (!gx->isINTRON && gx->captured && ! gc->captured)
	    continue ;
	  gt = arrayp (gts, gene, GT) ;
	  if (! gt->gene)
	    n1++ ;
	  gt->gene = gene ;
	  gt->intronGene = gc->intronGene ;
	  gt->intronType = gc->intronType ;
	  gt->intronAv = gc->intronAv ;
	  gt->intronSharedD = gc->intronSharedD ;
	  gt->intronSharedA = gc->intronSharedA ;
	  gt->length = gc->length ;
	  gt->isIntron = gc->isIntron ;
	  gt->geneGroup = gc->geneGroup ? 1 : 0 ;
	  gt->index[irun] = dc->index  - 1000 ; 
	  if (dc->isLow)
	    gt->index[irun] -= (dc->tags ? 1000 : 10000) ;
	  gt->tags[irun] = dc->tags ;
	    
	} 
    }
  fprintf (stderr, "# ===== n1 = %d\n", n1) ;
  for (gene = 1 ; gene < geneMax ; gene++)
    {
      double z0, z1 ;
      gt = arrayp (gts, gene, GT) ;
      if (! gt->gene) 
	continue ;
      z0 = gt->min = gt->max = gt->index[0] ; 
      if (z0 <  threshold) z0 =  threshold ;
      gt->colMax = gt->colMin = 0 ;
      for (irun = 1 ; irun < rMax ; irun++)
	{
	  char cc ;
	  z1 = gt->index[irun] ;
	  if (z1 <  threshold) z1 = threshold ;
	  
	  if (gt->min > z1) { gt->colMin = irun ; gt->min = z1 ; }
	  if (gt->max < z1) { gt->colMax = irun ; gt->max = z1 ; }

	  if (1)
	    {
	      if (z1 > z0 + .4)
		cc = '/' ;
	      else if (z1 < z0 - .4)
		cc = '\\' ;
	      else if (z1 == threshold)
		cc = '_' ;
	      else 
		cc = '~' ;
	      gt->shape[irun - 1] = cc ;
	    }
	  z0 = z1 ;
	}
      /*
	if (gt->min < 10) gt->min = 10 ;
	if (gt->max < 10) gt->max = 10 ;
      */
      gt->delta = gt->max - gt->min ;
      gt->score = -10 + gt->delta ;
      if (0 && gt->delta < .4) gt->delta = 0 ;
    }
  n1 = 0 ;
  for (gene = 1 ; gene < geneMax ; gene++)
    {
      int irun ;
      int nS = 1, nUp = 0, nDown = 0 ;
      gt = arrayp (gts, gene, GT) ;
      gt->score = 0 ;
      if (! gt->gene) 
	continue ;
      gc = arrayp (gx->genes, gt->gene, GC) ;
      if (wantGeneGroup == FALSE && gc->geneGroup)
	continue ;
      if (wantGeneGroup == TRUE && ! gc->geneGroup)
	continue ;
      if (!gx->isINTRON && (gx->targeted && !gc->targeted))
	continue ;
      if (gx->isINTRON && !gc->tags && ! gc->intronRefSeq && ! gc->intronAv)
	continue ;
      if (1)
	if (gx->captured && ! gc->captured)
	  continue ;

      if (1 && ! strcmp (dictName (gx->geneDict, gt->gene), "22__36651046_36653128"))
	invokeDebugger () ;
      if (rMax <= 15)
	{
	  int jrun ;
	  for (irun = 0 ; irun < rMax ; irun++)
	    for (jrun = irun + 1 ; jrun < rMax ; jrun++)
	      { 
		int run1 = keySet (compare->runs, irun) ;
		int run2 = keySet (compare->runs, jrun) ;
		COMPARE *cmp = gxPairCmp (gx, run1, run2) ;
		int isPm = run1 == cmp->vRun1 ? 1 : 2 ;
		Array scores1 = cmp ? cmp->scores[1] : 0 ;
		Array scores2 = cmp ? cmp->scores[2] : 0 ;
		float z1 = scores1 ? array (scores1, gt->gene, float) : 0 ;
		float z2 = scores2 ? array (scores2, gt->gene, float) : 0 ;
		float z = z1 > z2 ? z1 : z2 ;
		gt->sigShape[nS] = '~' ;
		
		if (
		    (z1 > z2 && z1 + .49 >= cmp->fdrThreshold[1]) ||
		    (z2 > z1 && z2 + .49 >= cmp->fdrThreshold[2])
		    )
		  {
		    gt->score += z ;
		    if (z > 0)
		      {
			if ((z1 - z2)*(2*isPm  - 3) < 0)
			  {
			    nDown++ ;
			    gt->sigShape[nS] = '\\' ;
			  }
			else
			  {
			    nUp++ ;
			    gt->sigShape[nS] = '/' ;
			  } 
		      }
		  }
		if(nS < 98) nS++ ;
	      }
	}
      else
	{
	  for (irun = 0 ; irun < rMax - 1 ; irun++)
	    {
	      int jrun = irun + 1 ;
	      int run1 = keySet (compare->runs, irun) ;
	      int run2 = keySet (compare->runs, jrun) ;
	      COMPARE *cmp = gxPairCmp (gx, run1, run2) ;
	      int isPm = run1 == cmp->vRun1 ? 1 : 2 ;
	      Array scores1 = cmp ? cmp->scores[1] : 0 ;
	      Array scores2 = cmp ? cmp->scores[2] : 0 ;
	      float z1 = scores1 ? array (scores1, gt->gene, float) : 0 ;
	      float z2 = scores2 ? array (scores2, gt->gene, float) : 0 ;
	      float z = z1 > z2 ? z1 : z2 ;
	      gt->sigShape[nS] = '~' ;

	      if (
		  (z1 > z2 && z1 + .49 >= cmp->fdrThreshold[1]) ||
		  (z2 > z1 && z2 + .49 >= cmp->fdrThreshold[2])
		  )
		{
		  gt->score = z ;
		  if (z > 0)
		    {
			if ((z1 - z2)*(2*isPm  - 3) < 0)
			{
			  nDown++ ;
			  gt->sigShape[nS] = '\\' ;
			}
		      else
			{
			  nUp++ ;
			  gt->sigShape[nS] = '/' ;
			}
		    }
		}
	      if(nS < 99) nS++ ;
	    }
	}
      gt->sigShape[nS] = 0 ;
      if (nUp && ! nDown)
	gt->sigShape[0] = '/' ;
      else if (! nUp && nDown)
	gt->sigShape[0] = '\\' ;
      else if (nUp && nDown)
	gt->sigShape[0] = 'x' ;
      else
	gt->sigShape[0] = 'z' ;
    }

  /* export */
  arraySort (gts, gtOrder) ; 
  oldSigShape = 0 ;
  for (gene = line = 0 ; gene < arrayMax (gts) ; gene++)
    {
      gt = arrayp (gts, gene, GT) ;
      if (! gt->gene)
	continue ;
      if (gx->isINTRON && ! gt->isIntron)
	continue ;
      if (wantGeneGroup != gt->geneGroup)
	continue ;
      switch (pass)
	{
	case 1:
	  if (gt->score <= 0) 
	    continue ;
	  break ;
	case 2:
	  if (gt->score < 185) 
	    continue ;
	  break ;
	default:
	  break ;
	}
      if (
	  (oldSigShape && strcmp (gt->sigShape, oldSigShape)) ||
	  (gene > 1 && gt->geneGroup != (gt-1)->geneGroup)
	  )
	aceOut (ao, "\n\n\n\n\n") ;
      oldSigShape = gt->sigShape ;

      if (keySet(bs, gt->gene))
	{
	  fprintf (stderr, "Double export of gene %s\n", dictName(gx->geneDict, gt->gene)) ;
	  invokeDebugger () ;
	}
      keySet (bs, gt->gene) = gene ;
      gc = arrayp (gx->genes, gt->gene, GC) ;
      aceOutf (ao, "%d\t%s" , ++line, noX__(dictName (gx->geneDict, gt->gene), !gx->isINTRON)) ;
      
      if (gx->hasGeneLength)
	aceOutf (ao, "\t%d", gt->length) ;

      if (gx->hasGeneId)
	{ 
	  if (gc->geneGroup)
	    aceOutf (ao, "\tGroup") ;
	  else
	    aceOutf (ao, "\t%s", gc->geneId ? stackText (gx->info, gc->geneId) : "AAA") ;
	}
      if (gx->isINTRON && gx->hasFromGene)
	{ 
	  if (gc->geneGroup)
	    aceOutf (ao, "\tGroup") ;
	  else 
	    aceOutf (ao, "\t%s", gc->intronGene ? dictName (gx->geneDict, gc->intronGene) : "YYY") ;
	}
      else if (gx->hasGeneModel)
	{ 
	  if (gc->geneGroup)
	    aceOutf (ao, "\tGroup") ;
	  else
	    aceOutf (ao, "\t%s", gc->geneModel ? stackText (gx->info, gc->geneModel) : "XXX") ;
	}
      if (gx->hasIntronType)
	{
	  if (gx->isINTRON && gt->isIntron && gx->affyDict)
	    aceOutf (ao, "\t%s", gt->intronType ? dictName (gx->affyDict, gt->intronType) : "") ;
	}
      if (gx->hasKnownDonorAcceptor)
	{
	  switch (gc->hasKnownDonorAcceptor)
	    {
	    case 0: aceOutf (ao, "\tNew") ; break ;
	    case 1:  aceOutf (ao, "\tNew_KnownD") ; break ;
	    case 2:  aceOutf (ao, "\tNew_KnownA") ; break ;
	    case 3:  aceOutf (ao, "\tNew_KnownAD") ; break ;
	    default: aceOutf (ao, "\tKnown") ; break ;
	    }
	}
      if (gx->hasIntronShared)
	{
	  if (gx->isINTRON && gt->isIntron && gx->affyDict)
	    aceOutf (ao, "\t%s", gt->intronSharedD ? dictName (gx->affyDict, gt->intronSharedA) : "") ;
	  if (gx->isINTRON && gt->isIntron && gx->affyDict)
	    aceOutf (ao, "\t%s", gt->intronSharedA ? dictName (gx->affyDict, gt->intronSharedA) : "") ;
	}
      if (gx->hasRefSeqAv)
	{
	  if (gx->isINTRON && gc->isIntron)
	    aceOutf (ao, "\t%d", gt->intronAv) ;
	}
      if (gx->hasCapture)
	{
	  KEYSET ks = gc->captures ;
	  int i ;
	  char *sep = "" ;
	  aceOutf (ao, "\t") ;
	  if (gx->isINTRON)
	    {
	      KEY i2g = gt->intronGene ;
	      GC *gc2 = i2g ? arrayp (gx->genes, i2g, GC) : 0 ;
	      ks = gc2 ? gc2->captures : 0 ;
	    }
	  if (ks)
	    for (i = 0 ; i < keySetMax (ks) ; i++)
	      {
		int capture = keySet (ks, i) ;
		if (capture)
		  {
		    const char *ccp = dictName (gx->captureDict, capture) ;
		    if (ccp && ccp[2] == 0) 
		      {
			aceOutf (ao, "%s%s", sep, ccp) ;
			sep = "," ;
		      }
		  } 
		if (! gt->capture)
		  {
		    for (irun = 0 ; irun < rMax ; irun++)
		      {
			run = keySet (compare->runs, irun) ;
			rc = arrp (gx->runs, run, RC) ;
			if (! rc->aa || rc->private)
			  continue ;
			if (rc->capture == capture)
			  {
			    gt->capture = capture ;
			    nCAPTURED++ ;
			    break ;
			  }
		      }
		  }
	      }
	}

      if (gx->hasCaptureTouch)
	{
	  KEYSET ks = gc->capturesTouch ;
	  KEYSET ksT = gc->capturesTouch ;
	  int i ;
	  char *sep = "" ;
	  aceOutf (ao, "\t") ;
	  if (gx->isINTRON)
	    {
	      KEY i2g = gt->intronGene ;
	      GC *gc2 = i2g ? arrayp (gx->genes, i2g, GC) : 0 ;
	      ks = gc2 ? gc2->captures : 0 ;
	      ksT = gc2 ? gc2->capturesTouch : 0 ;
	    }
	  if (ksT)
	    for (i = 0 ; i < keySetMax (ksT) ; i++)
	      {
		int capture = keySet (ksT, i) ;
		if (capture && ! (ks && keySetFind (ks, capture, 0)))
		  {
		    const char *ccp = dictName (gx->captureDict, capture) ;
		    if (ccp && ccp[2] == 0) 
		      {
			aceOutf (ao, "%s%s", sep, ccp) ;
			sep = "," ;
		      }
		  } 
	      }
	}
 
      if (0 && gx->hasGeneNm)
	{
	  aceOutf (ao,  "\t") ;
	  if (gc->nmid)
	    {
	      char *buf = strnew (stackText (gx->info, gc->nmid), 0) ;
	      
	      if (buf && *buf)
		aceOutf (ao, "%s", buf) ;

	      ac_free (buf) ;
	    }
	  else if (gc->geneId)
	    {
	      int gg = 0 ;
	      char *cp, *cq, *cr ;
	      char *buf = strnew (stackText (gx->info, gc->geneId), 0) ;
	      int jj = 0 ;
	      
	      cp = cq = buf ;
	      while (cq)
		{
		  cq = strstr (cp, ";") ;
		  if (cq) *cq++ = 0 ; 
		  if (dictFind (gx->geneIdDict, cp, &gg))
		    {
		      GC *gcid = arrp (gx->geneIds, gg, GC) ;
		      cr = stackText (gx->info, gcid->nmid) ;
		      
		      if (cr)
			aceOutf (ao, "%s%s"
				 , jj++ ? ";" : ""
				 , cr
				 ) ;
		    }
		  cp = cq ;
		}
	      ac_free (buf) ;
	    }
	}

      if (gx->hasGeneAlias)  aceOutf (ao, "\t%s", gc->alias ? stackText (gx->info, gc->alias) : "") ;
      if (gx->hasGeneChrom) 
	{
	  if (gc->chrom)
	    aceOutf (ao, "\t%s\t%c\t%d\t%d", dictName (gx->chromDict, gc->chrom), gc->strand, gc->a1, gc->a2) ;
	  else
	    aceOutf (ao, "\t\t\t\t") ;
	}
      if (gx->hasGeneType) aceOutf (ao, "\t%s", gc->gene_type ? stackText (gx->info, gc->gene_type) : "") ;
      if (gx->hasGeneTitle) 
	{
	  aceOutf (ao, "\t%s", gc->title ? stackText (gx->info, gc->title) : noX__(dictName (gx->geneDict, gt->gene), !gx->isINTRON)) ;
	  if (gx->hasGeneAlias && gc->alias) /* try to concatenate the title of the alias */
	    {
	      int gene2 = 0 ;
	      if (dictFind (gx->geneDict, stackText (gx->info, gc->alias), &gene2) && 
		  gene2 > 0 && gene2 < arrayMax (gx->genes)
		  )
		{
		  GC *gc2 = arrayp (gx->genes, gene2, GC) ;
		  if (gc2->title)
		    aceOutf (ao, ". %s",  stackText (gx->info, gc2->title)) ;
		}
	    }	  
	} 

      aceOutf (ao, "\t%s", noX__(dictName (gx->geneDict, gt->gene), !gx->isINTRON)) ;
      gt->max = -10000 ; gt->min = 10000 ;
      for (irun = 0 ; irun < rMax ; irun++)
	{
	  float z =  gt->index[irun] ;
	  char *prefix = "" ;

	  if (z < -8000)
	    { z += 10000 ; prefix = "NE:" ; }
	  else if (z < -800)
	    { z += 1000 ; prefix = "NA:" ; }
	  if (z > gt->max) 
	    { gt->max = z ; gt->colMax = irun ; }
	  if (z < gt->min) 
	    { gt->min = z ; gt->colMin = irun ; }
	  aceOutf (ao, "\t%s%.2f", prefix, z - FIX) ;
	}

      aceOutf (ao, "\t%d\t%d\t%s\t%s"
	       , gt->colMax + 1, gt->colMin + 1
	       , gt->shape
	       , gt->sigShape + 1
	       ) ;
      if (1)
	{
	  const char *ccp = 0 ;
	  ccp = noX__(dictName (gx->geneDict, gt->gene), !gx->isINTRON) ;
	  if (!ccp && gc->title) ccp = stackText (gx->info, gc->title) ;
	  if (!ccp && gc->Group_level) ccp = noX__(dictName (gx->geneDict, gt->gene), !gx->isINTRON) ;
	  if (!ccp) ccp =  noX__(dictName (gx->geneDict,  gt->gene), !gx->isINTRON) ;
	  aceOutf (ao, "\t%s", ccp ? ccp : "") ;
	}

      sum_score = sum_score1 = sum_score2 = 0 ;
      n_diff_scores = n_diff_scores1 = n_diff_scores2 = 0 ;

      if (0 && ! strcmp (dictName (gx->geneDict, gt->gene), "22__36651046_38087262"))
	invokeDebugger () ;
      float bestScore = -100 ;
      if (rMax <= 15)
	{
	  int jrun ;
	  for (irun = 0 ; irun < rMax ; irun++)
	    for (jrun = irun + 1 ; jrun < rMax ; jrun++)
	      { 
		float z ;
		float z1 = gt->index[irun] ;
		float z2 = gt->index[jrun] ;
		int isLow = 0 ;
		
		if (z1 < -8000)
		  { z1 += 10000 ; isLow++ ; }
		else if (z1 < -800)
		  { z1 += 1000 ; isLow++ ; }
		if (z2 < -8000)
		  { z2 += 10000 ; isLow += 2 ; }
		else if (z2 < -800)
		  { z2 += 1000 ; isLow += 2 ; }
		

		z = z1 - z2 ;
		if (z < 0)  z = -z ; 
 
		if (
		    z <  gx->minFoldChange ||
		    isLow == 3 ||   /* both low */
		    (isLow == 1 && z1 > z2) ||   /* low run cannot win */
		     (isLow == 2 && z1 < z2)      /* low run cannot win */
		    )
		  aceOut (ao, "\t\t") ;
		else
		  {
		    int type = 0 ;
		    int run1 = keySet (compare->runs, irun) ;
		    int run2 = keySet (compare->runs, jrun) ;
		    COMPARE *cmp = gxPairCmp (gx, run1, run2) ;
		    int isPm = run1 == cmp->vRun1 ? 1 : 2 ;
		    Array scores1 = cmp ? cmp->scores[1] : 0 ;
		    Array scores2 = cmp ? cmp->scores[2] : 0 ;
		    float z1 = scores1 ? array (scores1, gt->gene, float) : 0 ;
		    float z2 = scores2 ? array (scores2, gt->gene, float) : 0 ;
		    float z ;
		    float pValue = 10 ;		    
		
		    if (z1 < 0) z1 = 0 ;
		    if (z1 > 200) z1 = 200 ;
		    if (z2 < 0) z2 = 0 ;
		    if (z2 > 200) z2 = 200 ;
		    z = z1 > z2 ? z1 : z2 ;

		    if (
			(z1 > z2 && z1 + .49 >= cmp->fdrThreshold[1]) ||
			(z2 > z1 && z2 + .49 >= cmp->fdrThreshold[2])
			)
		      {
			if (z1 > z2)
			  {
			    type += isPm ;
			    pValue = cmp->fdr1[(int)z1] ;
			  }
			else
			  {
			    type += 3 - isPm ;
			    pValue = cmp->fdr2[(int)z2] ;
			  }
		      }
		      switch (type)
			{
			case 1:
			  if (z > bestScore)
			    bestScore = z ;
			  sum_score1 += z ;
			  n_diff_scores1++ ;
			  sum_score += z ;
			  n_diff_scores++ ;
			  
			  if (showPvalue)
			    aceOutf (ao, "\t%.2g\t", pValue) ;
			  else
			    aceOutf (ao, "\t%.2f\t", z) ;
			  gxRegisterComparedGenes (gx, gt, iCompare, irun, jrun, z, pValue) ;
			  break ;
			case 2:
			  if (z > bestScore)
			    bestScore = z ;
			  sum_score2 += z ;
			  n_diff_scores2++ ;
			  sum_score += z ;
			  n_diff_scores++ ;
		
			  if (showPvalue)
			    aceOutf (ao, "\t\t%.2g", pValue) ;
			  else
			    aceOutf (ao, "\t\t%.2f", z) ;
			  gxRegisterComparedGenes (gx, gt, iCompare, irun, jrun, z, pValue) ;
			  break ;
			case 3:
			  messcrash ("type = 3, both scores cannot be positive") ;
			  break ;
			default:
			  aceOut (ao, "\t\t") ;
		    }
		  }
	      }
	}
      else
	{
	  for (irun = 0 ; irun < rMax - 1 ; irun++)
	    {
	      int jrun = irun + 1 ;
	      
	      if (1)
		{
		  float z ;
		  float z1 = gt->index[irun] ;
		  float z2 = gt->index[jrun] ;
		  int isLow = 0 ;
		  
		  if (z1 < -8000)
		    { z1 += 10000 ; isLow++ ; }
		  else if (z1 < -800)
		    { z1 += 1000 ; isLow++ ; }
		  if (z2 < -8000)
		    { z2 += 10000 ; isLow += 2 ; }
		  else if (z2 < -800)
		    { z2 += 1000 ; isLow += 2 ; }
		  

		  z = z1 - z2 ;
		  if (z < 0)  z = -z ; 
		  
		  if (
		      z <  gx->minFoldChange ||
		      isLow == 3 ||   /* both low */
		      (isLow == 1 && z1 > z2) ||   /* low run cannot win */
		      (isLow == 2 && z1 < z2)      /* low run cannot win */
		      )
		    aceOut (ao, "\t\t") ;
		  
		  else
		    {
		      int type = 0 ;
		      int run1 = keySet (compare->runs, irun) ;
		      int run2 = keySet (compare->runs, jrun) ;
		      COMPARE *cmp = gxPairCmp (gx, run1, run2) ;
		      int isPm = run1 == cmp->vRun1 ? 1 : 2 ;
		      Array scores1 = cmp ? cmp->scores[1] : 0 ;
		      Array scores2 = cmp ? cmp->scores[2] : 0 ;
		      float z1 = scores1 ? array (scores1, gt->gene, float) : 0 ;
		      float z2 = scores2 ? array (scores2, gt->gene, float) : 0 ;
		      float z = 0 ;
		      float pValue = -10 ;		    
		      
		      if (z1 < 0) z1 = 0 ;
		      if (z1 > 200) z1 = 200 ;
		      if (z2 < 0) z2 = 0 ;
		      if (z2 > 200) z2 = 200 ;
		      z = z1 > z2 ? z1 : z2 ;

		      if (
			  (z1 > z2 && z1 + .49 >= cmp->fdrThreshold[1]) ||
			  (z2 > z1 && z2 + .49 >= cmp->fdrThreshold[2])
			  )
			{
			  if (z1 > z2)
			    {
			      type += isPm ;
			      pValue = cmp->fdr1[(int)z1] ;
			    }
			  else
			    {
			      type += 3 - isPm ;
			      pValue = cmp->fdr2[(int)z2] ;
			    }
			}
		      switch (type)
			{
			case 1:
			  if (z > bestScore)
			    bestScore = z ;
			  sum_score1 += z ;
			  n_diff_scores1++ ;
			  sum_score += z ;
			  n_diff_scores++ ;
			  
			  if (showPvalue)
			    aceOutf (ao, "\t%.2g\t", pValue) ;
			  else
			    aceOutf (ao, "\t%.2f\t", z) ;
			  gxRegisterComparedGenes (gx, gt, iCompare, irun, jrun, z, pValue) ;
			  break ;
			case 2:
			  if (z > bestScore)
			    bestScore = z ;
			  sum_score2 += z ;
			  n_diff_scores2++ ;
			  sum_score += z ;
			  n_diff_scores++ ;
		
			  if (showPvalue)
			    aceOutf (ao, "\t\t%.2g", pValue) ;
			  else
			    aceOutf (ao, "\t\t%.2f", z) ;
			  gxRegisterComparedGenes (gx, gt, iCompare, irun, jrun, z, pValue) ;
			  break ;
			case 3:
			  messcrash ("type = 3, both scores cannot be positive") ;
			  break ;
			default:
			  aceOut (ao, "\t\t") ;
			}
		    }
		}
	    }
	}

     
      /* best differential score */
      if (bestScore > -30)
	aceOutf (ao, "\t%.2f"
		 , bestScore
		 ) ;
      else
	aceOut (ao, "\t0") ;
      
      if (sum_score1 > -30)
	aceOutf (ao, "\t%.2f", sum_score1) ;
      else
	aceOut (ao, "\t") ;
      if (sum_score2 > -30)
	aceOutf (ao, "\t%.2f", sum_score2) ;
      else
	aceOut (ao, "\t") ;
      if (sum_score > -30)
	aceOutf (ao, "\t%.2f", sum_score) ;
      else
	aceOut (ao, "\t") ;

      aceOutf (ao, "\t%d", n_diff_scores) ;
      if (gt->capture)
	{
	  if (n_diff_scores1 && n_diff_scores2)
	    {
	      nBADcapt++ ;
	      aceOutf (ao, "\t%s",  n_diff_scores1 && n_diff_scores2 ? "Inconsistent" : "") ;
	    }
	  else
	    {
	      aceOut (ao, "\t") ;
	      if (n_diff_scores1 > 1)
		{ nGOOD1capt++ ; sumGoodScore1capt += sum_score1 ; nSCORE1capt += n_diff_scores1 ; }
	      if (n_diff_scores2 > 1)
		{ nGOOD2capt++ ; sumGoodScore2capt += sum_score2 ; nSCORE2capt += n_diff_scores2 ;}
	      if (n_diff_scores1 + n_diff_scores2 == 1)
		{ nUNIQUEcapt++ ; sumUniqueScorecapt += sum_score1 + sum_score2 ; }
	    }
	}
      else
	{
	  if (n_diff_scores1 && n_diff_scores2)
	    {
	      nBAD++ ;
	      aceOutf (ao, "\t%s",  n_diff_scores1 && n_diff_scores2 ? "Inconsistent" : "") ;
	    }
	  else
	    {
	      aceOut (ao, "\t") ;
	      if (n_diff_scores1 > 1)
		{ nGOOD1++ ; sumGoodScore1 += sum_score1 ; nSCORE1 += n_diff_scores1 ; }
	      if (n_diff_scores2 > 1)
		{ nGOOD2++ ; sumGoodScore2 += sum_score2 ; nSCORE2 += n_diff_scores2 ;}
	      if (n_diff_scores1 + n_diff_scores2 == 1)
		{ nUNIQUE++ ; sumUniqueScore += sum_score1 + sum_score2 ; }
	    }
	}
      aceOut (ao, "\t") ;
      if (gx->isINTRON)
	for (irun = 0 ; irun < rMax ; irun++)
	  aceOutf (ao, "\t%d", gt->tags[irun]) ;
      else
	for (irun = 0 ; irun < rMax ; irun++)
	  gxShowFpkm (ao, gt->index[irun] - FIX) ;
      gxShowFpkm (ao, gt->max - FIX) ;
      gxShowFpkm (ao, gt->min - FIX) ;
      aceOutf (ao, "\t%.2f", exp ((gt->max - gt->min) * ln2)) ;

      aceOutf (ao, "\t%.2f\t%.2f\t%.2f"
	       , gt->max - FIX, gt->min - FIX, gt->max - gt->min
	       ) ;

      if (gx->hasGeneAlias && gc->alias)  aceOutf (ao, "/%s", stackText (gx->info, gc->alias)) ;
      aceOut (ao, "\n") ;
    }

  aceOutf(ao, "\n\n## Statistics on %s", dictName (gx->compareDict, iCompare)) ;
  aceOutf(ao, "\tNumber of non targeted DEG genes\tAverage differential score\tAverage number of differential measures per gene") ;
  aceOutf(ao, "\tNumber of DEG genes among the %d targeted genes\tAverage differential score\tAverage number of differential measures per gene", nCAPTURED) ;
  aceOutf(ao, "\tNumber of DEG genes\tAverage differential score\tAverage number of differential measures per gene") ;

  aceOutf(ao, "\n# Good DEG with at least 2 consistent measures") ;
  aceOutf(ao, "\t%d\t%.2f\t%.2f",nGOOD1+nGOOD2, (sumGoodScore1+sumGoodScore2)/(.000001+nGOOD1+nGOOD2), (nSCORE1+nSCORE2)/(.000001+nGOOD1+nGOOD2)) ;
  aceOutf(ao, "\t%d\t%.2f\t%.2f",nGOOD1capt+nGOOD2capt, (sumGoodScore1capt+sumGoodScore2capt)/(.000001+nGOOD1capt+nGOOD2capt), (nSCORE1capt+nSCORE2capt)/(.000001+nGOOD1capt+nGOOD2capt)) ;
  aceOutf(ao, "\t%d\t%.2f\t%.2f"
	  , nGOOD1+nGOOD2+nGOOD1capt+nGOOD2capt
	  , (sumGoodScore1+sumGoodScore2+sumGoodScore1capt+sumGoodScore2capt)/(.000001+nGOOD1+nGOOD2+nGOOD1capt+nGOOD2capt)
	  , (nSCORE1+nSCORE2+nSCORE1capt+nSCORE2capt)/(.000001+nGOOD1+nGOOD2+nGOOD1capt+nGOOD2capt)
	  ) ;

  aceOutf(ao, "\n# Good DEG with at least 2 consistent measures A > B") ;
  aceOutf(ao, "\t%d\t%.2f\t%.2f",nGOOD2, sumGoodScore2/(.000001+nGOOD2), nSCORE2/(.000001+nGOOD2)) ;
  aceOutf(ao, "\t%d\t%.2f\t%.2f",nGOOD2capt, sumGoodScore2capt/(.000001+nGOOD2capt), nSCORE2capt/(.000001+nGOOD2capt)) ;
  aceOutf(ao, "\t%d\t%.2f\t%.2f"
	  , nGOOD2+nGOOD2capt
	  , (sumGoodScore2+sumGoodScore2capt)/(.000001+nGOOD2+nGOOD2capt)
	  , (nSCORE2+nSCORE2capt)/(.000001+nGOOD2+nGOOD2capt)
	  ) ;

  aceOutf(ao, "\n# Good DEG with at least 2 consistent measures A < B") ;
  aceOutf(ao, "\t%d\t%.2f\t%.2f",nGOOD1, sumGoodScore1/(.000001+nGOOD1), nSCORE1/(.000001+nGOOD1)) ;
  aceOutf(ao, "\t%d\t%.2f\t%.2f",nGOOD1capt, sumGoodScore1capt/(.000001+nGOOD1capt), nSCORE1capt/(.000001+nGOOD1capt)) ;
  aceOutf(ao, "\t%d\t%.2f\t%.2f"
	  , nGOOD1+nGOOD1capt
	  , (sumGoodScore1+sumGoodScore1capt)/(.000001+nGOOD1+nGOOD1capt)
	  , (nSCORE1+nSCORE1capt)/(.000001+nGOOD1+nGOOD1capt)
	  ) ;

  aceOutf(ao, "\n# DEG with only one measure") ;
  aceOutf(ao, "\t%d\t%.2f\t%d",nUNIQUE, sumUniqueScore/(.000001+nUNIQUE), 1) ; 
  aceOutf(ao, "\t%d\t%.2f\t%d",nUNIQUEcapt, sumUniqueScorecapt/(.000001+nUNIQUEcapt), 1) ; 
  aceOutf(ao, "\t%d\t%.2f\t%d",nUNIQUE+nUNIQUEcapt, (sumUniqueScore+sumUniqueScorecapt)/(.000001+nUNIQUE+nUNIQUEcapt), 1) ; 

  aceOutf(ao, "\n# Inconsistent DEG\t%d\t\t\t%d\t\t\t%d\t\t", nBAD, nBADcapt,nBAD+nBADcapt) ;
  aceOutf(ao, "\n") ;

  ac_free (h) ;
  return ;
}  /* gxOneExpressionProfile */

/*************************************************************************************/

static BOOL gxPairEvaluate (GX *gx, COMPARE *compare)
{
  BOOL ok = FALSE ;

  fprintf (stderr, "....... PairEvaluate A :: %s %s  compare=%s\n"
	   , dictName (gx->runDict, compare->vRun1)
	   , dictName (gx->runDict, compare->vRun2)
	   , dictName (gx->pairDict, compare->compare)
	   ) ;
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, ": %s %d Mb\n", timeShowNow(), mx) ;
  }

  ok = gxCompareAndControl (gx, compare) ;
  compare->done = TRUE ;
  
  return ok ;
} /* gxPairEvaluate */

/*************************************************************************************/

static void  gxCreatePrivatePairs (GX *gx, int iCompare)
{
  COMPARE *cmp, *compare = arrayp(gx->compares, iCompare, COMPARE) ;
  int irun, jrun, run1, run2 ;
  KEYSET runs = compare->runs ;
  int rMax = runs ? keySetMax(runs) : 0;  

#ifdef JUNK
  RC *rc ;
  double threshold = 8.5 ; 

  for (irun = 0 ; irun < rMax ; irun++)
    {
      run = keySet (compare->runs, irun) ;
      rc = arrp (gx->runs, run, RC) ;
      if (! rc->aa || rc->private)
	continue ;
     if (threshold < rc->NA_crossover_index)
	threshold = rc->NA_crossover_index ;
    }
#endif

  for (irun = 0 ; irun < rMax ; irun++)
    for (jrun = irun + 1 ; jrun < rMax ; jrun++)
      { 
	if (compare->profile && rMax > 15 && jrun != irun +1)
	  continue ;
	if (1)
	  {
	    run1 = keySet (runs, irun) ;
	    run2 = keySet (runs, jrun) ;
	  }
	else
	  {  /* hack, what happens if we switch, does it affect the results ? 2023_10_5 
	     * it seems quasi reproducible, except random 
	     * does not affect the same sequence to the  same runs
	     */
	    run1 = keySet (runs, rMax - irun - 1) ;
	    run2 = keySet (runs, rMax - jrun - 1) ;
	  }
	cmp = gxPairCmp (gx, run1, run2) ;
	if (! cmp->done)
	  gxPairEvaluate (gx, cmp) ;
     }
} /* gxCreatePrivatePairs  */
 
/*************************************************************************************/

static void gxExpressionProfile (GX *gx, int pass)
{
  AC_HANDLE h = ac_new_handle () ;
  int iCompare = 0 ;
  COMPARE *compare ;
  DICT *allPairs = dictHandleCreate (200, h) ;

  if (pass == 0)
    {
      for (iCompare = 0 ; iCompare < arrayMax (gx->compares) ; iCompare++)
	{ 
	  char buf[1024] ;
	  
	  compare = arrp(gx->compares, iCompare, COMPARE) ;
	  if (gx->isINTRON && ! compare->compareINTRON)
	    continue ;
	  if (gx->isSNP && ! compare->compareSNP)
	    continue ;
	  if (compare->runs && keySetMax (compare->runs) == 2)
	    {
	      int run1 = keySet (compare->runs, 0) ;
	      int run2 = keySet (compare->runs, 1) ;
	      if (run1 > run2)
		{ int run0 = run1 ; run1 = run2 ; run2 = run0 ; } 
	   sprintf (buf, "%d_%d", run1, run2) ;
	   
	   dictAdd (allPairs, buf, 0) ;
	    }
	}
    }

   /* ATTENTION, these functions may increase and realloc gx->compares */
   if (pass == 0)
     {
       for (iCompare = 1 ; iCompare < arrayMax (gx->compares) ; iCompare++)
	 {
	   compare = arrp(gx->compares, iCompare, COMPARE) ;
	   if (gx->isINTRON && ! compare->compareINTRON)
	     continue ;
	   if (gx->isSNP && ! compare->compareSNP)
	     continue ;
	   if (compare->private)
	     continue ;
	   gxCreatePrivatePairs (gx, iCompare) ;
	 }
     }
   else
     {
       int v ;
       int vMax = 1 ; /* Danielle 2019_05_30 only wants all profiles (v==0) and no longer diff_profiles (v==1) or top_diff (v==2) */
       for (v = 0 ; v < vMax ; v++)
	 {               /* pass , showPvalue, wantGeneGroup */
	   for (iCompare = 0 ; iCompare < arrayMax (gx->compares) ; iCompare++)
	     {
	       compare = arrp(gx->compares, iCompare, COMPARE) ;
	       if (! compare->private && compare->expressionProfile)
		 {
		   gxOneExpressionProfile (gx, iCompare, v, TRUE, FALSE) ;
		   gxOneExpressionProfile (gx, iCompare, v, FALSE, FALSE) ;
		   if (gx->maxGeneGroupLevel)
		     {
		       gxOneExpressionProfile (gx, iCompare, v, TRUE, TRUE) ;
		       gxOneExpressionProfile (gx, iCompare, v, FALSE, TRUE) ;
		     }
		 }
	     }
	 }
     }
   {
     int mx ;
     messAllocMaxStatus (&mx) ;   
     fprintf (stderr, "// gxExpressionProfile done :  %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
   }
   ac_free (h) ;
} /* gxExpressionProfile */

/*************************************************************************************/
/*************************************************************************************/

static void gxMixingRatios (GX *gx)
{
  AC_HANDLE h = ac_new_handle () ;
  int iCompare = 0 ;
  COMPARE *compare ;
  
   for (iCompare = 0 ; iCompare < arrayMax (gx->compares) ; iCompare++)
     { 
       compare = arrp(gx->compares, iCompare, COMPARE) ;
       if (! compare->mixing_ratios)
	 continue ;
       if (keySetMax (compare->runs) < 3)
	 continue ;
       if (keySetMax (compare->ratios) < 3)
	 continue ;
     }

   ac_free (h) ;
} /* gxMixingRatios */

/*************************************************************************************/
/*************************************************************************************/

static void  gxOneSamplePairingExport (ACEOUT ao, GX *gx, CS *cs1, CS *cs2, int pass)
{
  switch (pass)
    {
    case 0:
      aceOutf (ao, "\t%.2f", - cs2->weight/10000.0) ;
      break ;
    case 1:
      if (- cs2->weight/10000.0 > 0)
	{
	  int run1 = cs1->run ;
	  RC *rc1 = arrayp (gx->runs, run1, RC) ;
	  int run2 = cs2->run ;
	  RC *rc2 = arrayp (gx->runs, run2, RC) ;
	  
	  aceOutf (ao, "\t%s::%s"
		   , rc1->title ? stackText (gx->info, rc1->title) : dictName (gx->runDict, run1)
		   , rc2->title ? stackText (gx->info, rc2->title) : dictName (gx->runDict, run2)
		   ) ;
	}
      else
	aceOut (ao, "\t") ;
      break ;
    case 2:
      aceOutf (ao, "\t%d", cs2->nn) ;
      break ;
    }
}  /* gxOneSamplePairingExport */

/*************************************************************************************/

typedef struct { int run1, run2, nn ; float score ; } GOSP ;
static int gospOrder (const void *a, const void *b)
{
  const GOSP *up = (const GOSP *)a, *vp = (const GOSP *)b ;
  float z = up->score - vp->score ;
  int n ;

  if (z > 0) return -1 ; /* bigScore first */
  if (z < 0) return 1 ;
  n =  up->run1 - vp->run1 ; if (n) return n ;
  n =  up->run2 - vp->run2 ; if (n) return n ;
  return 0 ;
}

/*************************************************************************************/

static BOOL gxWhoIsWhoOneChronoOrdering (GX *gx, Array aa, float delta, float z, Array runs1, Array runs2, int iCompare, COMPARE *compare)
{ 
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (gx->outFileName, messprintf (".%s.who_is_who.clustered.txt", dictName (gx->compareDict, iCompare)), FALSE, h) ;
  Array bb ;
  KEYSET lines, cols ;
  int i, j, iMax, jMax ;
  float minZ = gx->minWhoScore ;

  aceOutDate (ao, "###", gx->title) ;
  aceOutf (ao, "## Pairs of runs, probably coming from the same individual, are sorted\n") ;
  aceOutf (ao, "## The comparison is based on the correlation of substitution SNPs seen as high in at least 2 runs\n") ;
  aceOutf (ao, "## The differential correlation is given for each pair, and the table is sorted by chronological ordering\n") ;

  iMax = arrayMax (gx->runs) ;
  jMax = arrayMax (gx->runs) ;
  if (0 && iMax > 8) iMax = jMax = 8 ;  /* hack to debug more easily when needed */
  bb = arrayHandleCreate (iMax * jMax, double, h) ;
  array (bb, iMax * jMax - 1, double) = 0 ;
  for (i = 0 ; i < arrayMax (aa) ; i++)
    {
      GOSP *up = arrayp (aa, i, GOSP) ;
      int run1 = up->run1 ;
      int run2 = up->run2 ;
      float z1 = (up->score + delta) * z ;
      if (z1 < minZ) z1 = 0 ;
      if (run1 >= iMax || run2 >= iMax)
	continue ;
      array (bb, run1 * iMax + run1, double) = 100.0 ;
      array (bb, run2 * iMax + run2, double) = 100.0 ;
     
      array (bb, run1 * iMax + run2, double) = z1 ;
      array (bb, run2 * iMax + run1, double) = z1 ;
    }

  lines = keySetHandleCreate (h) ;
  cols = keySetHandleCreate (h) ;
  keySet (lines, iMax - 1) = 0 ;
  keySet (cols, jMax - 1) = 0 ;
  topoChronoOrder (bb, lines, cols) ;

  /* show */

  if (gx->hasRunTitle)
    {
      for (j = 0 ; j < jMax ; j++)
	if (keySet (cols, j) &&  (0 || array (bb, j * jMax + j, double) > 0))
	  {
	    RC *rc = arrayp (gx->runs, keySet (cols, j), RC) ;
	    aceOutf (ao, "\t%s", rc->title ? stackText (gx->info, rc->title) : "") ;
	  }
      aceOutf (ao, "\n") ;
    }
  for (j = 0 ; j < jMax ; j++)
    if (keySet (cols, j) &&  (0 || array (bb, j * jMax + j, double) > 0))
      aceOutf (ao, "\t%s", dictName (gx->runDict, keySet (cols, j))) ;
  for (i = 0 ; i < iMax ; i++)
    {
      if (keySet (lines, i) &&  (0 || array (bb, i * jMax + i, double) > 0))
	{
	  aceOutf (ao, "\n%s", dictName (gx->runDict, keySet (lines, i))) ;
	  for (j = 0 ; j < jMax ; j++) 
	    if (keySet (cols, j) &&  (0 || array (bb, j * jMax + j, double) > 0))
	      {
		float z = array (bb, i * jMax + j, double) ;
		if (z > 0)
		  aceOutf (ao, "\t%.2f", z) ;
		else
		  aceOut (ao, "\t") ;
	      }
	  if (gx->hasRunTitle) 
	    {
	      RC *rc = arrayp (gx->runs, keySet (lines, i), RC) ;
	      aceOutf (ao, "\t%s", rc->title ? stackText (gx->info, rc->title) : "") ;
	    }
	}
    }
  aceOut (ao, "\n") ;

  ac_free (h) ;
  return TRUE ;
} /* gxWhoIsWhoOneChronoOrdering */

/*************************************************************************************/

static BOOL gxOneSamplePairing (ACEOUT ao, GX *gx, int iCompare, COMPARE *compare, BOOL isDiff)
{
  AC_HANDLE h = ac_new_handle () ;
  Array hh = 0, aa = 0 ;
  int lowAll = 12 ;
  int runMax = arrayMax (gx->runs) ;
  int geneMax = arrayMax (gx->genes) ;
  int oldDigital = gx->digital ;
  int ii = 0 ;
  int irun1, irun2, n, gene, pass, run1, run2, low ;
  Array runs1 = arrayHandleCreate (runMax, CS, h), runs2 = 0 ;
  RC *rc1, *rc2 ;
  GC *gc ;
  CS* cs1, *cs2 ;
  KEYSET genes = keySetHandleCreate (h) ;
  double w ;
  gx->digital = gx->digitalCorrelation ;  
  
  if (0) fprintf (stderr, "######### gx->digitalCorrelation=%d\n", gx->digital) ;
  /* create a searchable keyset */
  for (n = 0, low = 0, irun1 = 0 ; irun1 < keySetMax (compare->runs) ; irun1++)
    {
      run1 = keySet (compare->runs, irun1) ;
      rc1 = arrayp (gx->runs, run1, RC) ;

      cs1 = arrayp (runs1, n++, CS) ;     /* twin runs have consecutive numbers */
      cs1->run = run1 ; cs1->irun = irun1 ;
      if (rc1->aa && lowAll < rc1->NA_crossover_index)
	lowAll = rc1->NA_crossover_index ;
    }
  if (gx->isSNP)
    lowAll = 4 ;
  else
    for (n = gene = 0, gc = arrayp (gx->genes, gene, GC) ; gene < geneMax ; gc++, gene++)
      {
	if (runMax > 10 && gc->maxIndex - gc->minIndex < 2)
	  continue ;
	if (gc->maxIndex < lowAll)
	  continue ;
	keySet (genes, n++) = gene ;
      } 
  keySetSort(genes) ; keySetCompress (genes) ;  
  runs2 = arrayHandleCopy (runs1, h) ;

  /* for each run compare->runs */


  if (compare->who_is_who)
    {
      ii = arrayMax (runs1) * arrayMax (runs2) ;
      hh = arrayHandleCreate (ii, float, h) ;
      aa = arrayHandleCreate (ii, GOSP, h) ;
    }

  for (ii = 0, low = 0, irun1 = 0 ; irun1 < arrayMax (runs1) ; irun1++)
    {
      cs1 = arrayp (runs1, irun1, CS) ;  
      run1 = cs1->run ;
      rc1 = arrayp (gx->runs, run1, RC) ;

      if (rc1->private || rc1->isSubLib) continue ;
      lowAll = 12 ;

      /* compute the correlation with all other compare->runs */
      for (irun2 = 0 ; irun2 < arrayMax (runs2) ; irun2++)
	{
	  cs2 = arrayp (runs2, irun2, CS) ;  
	  if (1 && cs1->irun == cs2->irun)
	    {
	      cs2->weight = - 1000000 ;
	      cs2->nn = arrayMax (gx->genes) ;
	      continue ;
	    }
	  else
	    {
	      int nn = 0 ;
	      run2 = cs2->run ;
	      rc2 = arrayp (gx->runs, run2, RC) ;

	      if (rc2->private || rc2->isSubLib) continue ;
	     
	      if (gx->isSNP) 
		{ low = 0 ; isDiff = FALSE ; }
	      else
		{
		  low = rc1->NA_crossover_index > rc2->NA_crossover_index ? rc1->NA_crossover_index : rc2->NA_crossover_index ;
		  if (low < lowAll) low = lowAll ;
		}
	      w = 0 ;
	      gxRunsCosineUsingGenes (gx, genes, 0,0,0,0, run1, run2, &w, &nn, low, isDiff) ;
	      if(0 && (irun1+irun2<12 || 100 * w > 50)) fprintf (stderr, "Pairing %d\t%.1f\t%s\t%s\n",nn,100*w, dictName (gx->runDict, cs1->run), dictName (gx->runDict, cs2->run)) ; 
	      if (nn >= 20)  /* dromadaire minimal number of measured SNP */
		{
		  cs2->weight = - 1000000.0 * w ; /* minus, because we reuse csOrder which is increasing, sorry */
		  cs2->nn = nn ;
		}
	      else
		cs2->nn = 0 ;
	    }
	}
      /* sort an export the twin and the 5 best */
      if (1) arraySort (runs2, csOrder) ;
      if (ao && ! compare->who_is_who)
	{
	  aceOutf (ao, "%s\t%s", dictName (gx->runDict, cs1->run), rc1->title ? stackText (gx->info, rc1->title) : "") ;
	  
	  for (pass = 0 ; pass < 3 ; pass++)
	    {
	      if (0)
		{
		  /* look for the twin */
		  for (irun2 = 0, cs2 = arrayp (runs2, irun2, CS)  ; irun2 < arrayMax (runs2) ; cs2++, irun2++)
		    if (cs2->irun == cs1->irun + (cs1->irun % 2 == 0 ? 1 : -1))
		      gxOneSamplePairingExport (ao, gx, cs1, cs2, pass) ;
		}
	      /* now export the 20 best */
	      for (irun2 = 0, cs2 = arrayp (runs2, irun2, CS)  ; irun2 < 20  ; cs2++, irun2++)
		{
		  if (cs2->nn > 0 &&  irun2 < arrayMax (runs2))
		    gxOneSamplePairingExport (ao, gx, cs1, cs2, pass) ; 
		  else
		    aceOut (ao, "\t") ;
		}
	      aceOut (ao, pass == 2 ? "\n" : "\t") ;	  
	    }
	}
      else if (aa)
	{
	  for (irun2 = 0, cs2 = arrayp (runs2, irun2, CS)  ; irun2 < arrayMax (runs2) ; cs2++, irun2++)
	    if ( cs1->run && cs2->run && cs1->run < cs2->run &&  - cs2->weight/10000.0 > -230 && cs2->nn > 0)
	      {
		GOSP *up = arrayp (aa, ii, GOSP) ;
		
		up->run1 = cs1->run ;
		up->run2 = cs2->run ;
		array (hh, ii, float) = up->score = - cs2->weight/10000.0 ;
		up->nn = cs2->nn ;
		ii++ ;
	      }
	}
    }
  
  if (compare->who_is_who)
    {
      /* rescale the SNP score to 0 100
       * the median is often at delta = -30
       * we want to move this to zero and keep 100 at 100
       *  score -> (score + delta) * 100/(100+delta)
       */
      float z, delta ;
      ACEOUT ao = aceOutCreate (gx->outFileName, messprintf (".%s.who_is_who.txt", dictName (gx->compareDict, iCompare)), FALSE, h) ;

      aceOutDate (ao, "###", gx->title) ;
      aceOutf (ao, "## Pairs of runs, probably coming from the same individual, are listed\n") ;
      aceOutf (ao, "## The comparison is based on the correlation of substitution SNPs seen as high in at least 2 runs\n") ;
      aceOutf (ao, "# Run 1\tRun 2\tDifferential correlation\tNumber of SNPs\tTitle 1\tTitle 2\n") ;



      arraySort (aa, gospOrder) ;
      arraySort (hh, floatOrder) ;
      ii = arrayMax (hh) / 2 ;
      delta = - array (hh, ii, float) ;
      if (0 && delta < 0) delta = 0 ;
      z = 100 / (100 + delta) ;

      for (ii = 0 ; ii < arrayMax (aa) ; ii++)
	{
	  GOSP *up = arrayp (aa, ii, GOSP) ;
	  int run1 = up->run1 ;
	  RC *rc1 = arrayp (gx->runs, run1, RC) ;
	  int run2 = up->run2 ;
	  RC *rc2 = arrayp (gx->runs, run2, RC) ;
		
	  aceOutf (ao, "%s\t%s"
		   , dictName (gx->runDict, run1)     
		   , dictName (gx->runDict, run2)     
		   ) ;
	  aceOutf (ao, "\t%.2f", (up->score + delta) * z) ;
	  aceOutf (ao, "\t%d", up->nn) ;
	  aceOutf (ao, "\t%s\t%s"
		   , rc1->title ? stackText (gx->info, rc1->title) : dictName (gx->runDict, run1)
		   , rc2->title ? stackText (gx->info, rc2->title) : dictName (gx->runDict, run2)
		   ) ;
	  aceOutf (ao, "\n") ;
	}

      gxWhoIsWhoOneChronoOrdering (gx, aa, delta, z, runs1, runs2, iCompare, compare) ;
    }
  gx->digital = oldDigital ;  
  ac_free (h) ;
  return TRUE ;
} /* gxOneSamplePairing */

/*************************************************************************************/
/* find the correlation between all runs in the compare->runs
 * and list the 5 best matches
 */
static void gxSamplePairing (GX *gx)
{
  int n, iCompare = 0 ;
  COMPARE *compare ;

  if (!gx->compares || ! arrayMax (gx->compares))
    return ;
  
  for (iCompare = 0, compare = arrp(gx->compares, iCompare, COMPARE) ; iCompare < arrayMax (gx->compares) ; compare++, iCompare++)
    {
       if (gx->isINTRON && ! compare->compareINTRON)
	 continue ;
       if (gx->isSNP && ! compare->compareSNP)
	 continue ;
       if (0 && !gx->isSNP && compare->compareSNP)
	 continue ;

      if (compare->who_is_who && compare->runs && keySetMax(compare->runs))
	 {
	   gxOneSamplePairing (0, gx, iCompare, compare, FALSE) ;
	 }
      else if (! gx->isSNP && compare->samplePairing && compare->runs && keySetMax(compare->runs))
	 {
	   AC_HANDLE h = ac_new_handle () ;
	   ACEOUT ao = aceOutCreate (gx->outFileName, messprintf (".%s.sample_pairing.txt", dictName (gx->compareDict, iCompare)), FALSE, h) ;
	   aceOutDate (ao, "###", gx->title) ;
	   aceOutf (ao, "# Run\tTitle") ;
	   for (n = 1 ; n<=20 ; n++)
	     aceOutf (ao, "\tClosest %d", n) ;
	   aceOut (ao,"\t") ;
	   for (n = 1 ; n<=20 ; n++)
	     aceOutf (ao, "\tClosest %d", n) ;
	   for (n = 1 ; n<=20 ; n++)
	     aceOutf (ao, "\tClosest %d", n) ;
	   aceOut (ao,"\n") ;

	   gxOneSamplePairing (0, gx, iCompare, compare, TRUE) ;
	   gxOneSamplePairing (0, gx, iCompare, compare, FALSE) ;

	   ac_free (h) ;
	 }
     }
} /* gxSamplePairing */

/*************************************************************************************/
/*************************************************************************************/
/* classification by signatures
 * Starting from a compare
 * we select in each compare->runs (called the training set) the most characteristic
 * SNPs or genes, defined as the top N which are higly expressed, or m/m in each
 * training run (or subgroup) and the least possible number of other training runs
 * Then we count for all runs of the project the number of characteristic gene/SNP
 * they share with each training run, then for each run we sort the training
 * runs by decreasing order and then sort all runs by that signature
 */

typedef struct gxsnpStruct { 
  AC_HANDLE h ;
  BitSet isHigh ;      /* SNP high, indeex by s * runMax + run */
  BitSet isMeasured ;  /* SNP measured, indeex by s * runMax + run */
  Array sharedMatrix ;      /* class characteristic SNPs also high in run     : relevant numerator */
  Array measuredMatrix ;  /* class characteristic SNPs also measured in run : relevant denominator */

  BigArray snps ;
  DICT *snpDict ;
  KEYSET iRuns, jRuns, iRunsSorted, iCounts ;
} GXSNP ;

typedef struct snpStruct { int nm, nm2, nh, nh2, ok ; } SNP ;

static GXSNP *gxSnpInit (GX *gx)
{
  int ri ;
  GXSNP *gxsnp = halloc (sizeof(GXSNP), gx->h) ;

  gxsnp->h = ac_new_handle () ;
  gxsnp->snps = bigArrayHandleCreate (1000000, SNP, gxsnp->h) ;
  gxsnp->snpDict = dictHandleCreate (100000, gxsnp->h) ;

  gxsnp->iCounts = keySetHandleCreate (gxsnp->h) ;
  gxsnp->iRuns = keySetHandleCreate (gxsnp->h) ;
  gxsnp->jRuns = keySetHandleCreate (gxsnp->h) ;
  gxsnp->isHigh = bitSetCreate (1000000, gxsnp->h) ;
  gxsnp->isMeasured = bitSetCreate (1000000, gxsnp->h) ;

  for (ri = dictMax (gx->runDict) ; ri >= 1 ; ri--)
    keySet (gxsnp->iRuns, ri) = ri ;

  return gxsnp ;
} /* gxSnpInit */

/*************************************************************************************/
/* parse the snps in the format generated in directory SNPH by the snp.c program */
static void gxSnpParse (GX *gx, GXSNP *gxsnp)
{
  AC_HANDLE h  = ac_new_handle () ;
  ACEIN ai = aceInCreate (gx->snpFileName, gx->gzi, h) ;
  const char *ccp ;
  int mx, count, cover, wild, s, run, ki ;
  float genotype = 0 ;
  Stack buf = stackHandleCreate (256, h) ;
  SNP *snp ;
  BitSet HH = gxsnp->isHigh ;
  BitSet MM = gxsnp->isMeasured ;
  BigArray snps = gxsnp->snps ;
  DICT *snpDict = gxsnp->snpDict ;
  KEYSET iCounts = gxsnp->iCounts ;
  int riMax = dictMax (gx->runDict) + 1 ;
  int rMax = riMax ; 
  int nnn0 = 0, nnn1 = 0 ;

  messAllocMaxStatus (&mx) ;  
  fprintf (stderr, "... gxSnpParse       \t%s\tmax memory %d Mb\tstart\n", timeShowNow(), mx) ;
  aceInSpecial (ai, "\"\n\t") ;
  while (aceInCard (ai))
    { 
      ccp = aceInWord (ai) ;
      if (! ccp) continue ;
      if (*ccp =='#') continue ;
      stackClear (buf) ;
      pushText (buf, ccp) ;
      catText (buf, ":") ;
       aceInStep (ai,'\t') ; ccp = aceInWord (ai) ; if (! ccp) continue ;
      catText (buf, ccp) ;
      catText (buf, ":") ;
       aceInStep (ai,'\t') ; ccp = aceInWord (ai) ; if (! ccp) continue ;
      catText (buf, ccp) ;
      catText (buf, ":") ;

      aceInStep (ai,'\t') ; ccp = aceInWord (ai) ; if (! ccp) continue ; catText (buf, ccp) ;  catText (buf, ":") ;
      aceInStep (ai,'\t') ; ccp = aceInWord (ai) ; if (! ccp) continue ; catText (buf, ccp) ;  catText (buf, ":") ;
      
      aceInStep (ai,'\t') ; ccp = aceInWord (ai) ; if (! ccp) continue ;
      run = 0 ;
      if (! gx->runListFileName)
	dictAdd (gx->runDict, ccp, &run) ;
      else
	dictFind (gx->runDict, ccp, &run) ;
      if (run <= 0)
	continue ;
       aceInStep (ai,'\t') ; if (! aceInFloat (ai, &genotype))
	continue ;
    
       aceInStep (ai,'\t') ; ccp = aceInWord (ai) ; /* frequency */ if (! ccp) continue ;
      
      cover = count = wild = 0 ;
      aceInStep (ai,'\t') ; aceInInt (ai, &cover) ;
      aceInStep (ai,'\t') ; aceInInt (ai, &count) ;
      aceInStep (ai,'\t') ; aceInInt (ai, &wild) ;
      
      if (cover < 10)
	continue ;
      if (0 && (100 * count < 95 * cover && 100 * wild < 95 * cover))
	continue ;
 
      {
	int i = 8 , dummy ;
	while (i--)
	  { aceInStep (ai,'\t') ; aceInInt (ai, &dummy) ; }  /* jump the counts per strand and the error rates */
	 aceInStep (ai,'\t') ; ccp = aceInWord (ai) ; if (! ccp) continue ;
	 if (! strcmp (ccp, "Incompatible_strands")) continue ;   /* we certainly need that */
	 if (0 && strcmp (ccp, "Both_strands")) continue ;    /* this fileter is optional */
      }
      dictAdd (snpDict, stackText(buf,0), &s) ;
      snp = bigArrayp (snps, s, SNP) ;
      ki = s * rMax + run ;
      if (bit (MM, ki)) /* do not double count if the datafile has repetitions */
	continue ; 
      nnn0++ ; 
      bitSet (MM, ki) ; snp->nm++ ; 
      if (100 * count >= 95 * cover)
	{
	  bitSet (HH, ki) ; snp->nh++ ; 
	  keySet (iCounts, run)++ ;
	  nnn1++ ;
	}
    }
  messAllocMaxStatus (&mx) ;  
  fprintf (stderr, "... gxSnpParse       \t%s\tmax memory %d Mb\tfound %d measured, %d high SNPs\n", timeShowNow(), mx, nnn0, nnn1) ;
  ac_free (h) ;
} /* gxSnpParse */

/*************************************************************************************/

static void gxSnpSelect (GXSNP *gxsnp, BOOL unSelect)
{
  int s, sMax = dictMax (gxsnp->snpDict) + 1 ;
  SNP *snp ;
  BigArray snps = gxsnp->snps ;

  for (s = 1, snp = bigArrp (snps, s, SNP) ; s < sMax ; snp++, s++)
    {
      if (unSelect)
	{
	  if (snp->ok == 1)
	    snp->ok = 0 ;
	}
     }
} /* gxSnpSelect */

/*************************************************************************************/
/* for each SNP in jRuns , count in how many other jRuns they appear and select the N most specific */
static void gxSnpCountSpecificSNPs (ACEOUT ao, GX *gx, GXSNP *gxsnp, int iCompare, int NN)
{
  AC_HANDLE h = ac_new_handle () ;
  BitSet HH = gxsnp->isHigh ;
  BitSet MM = gxsnp->isMeasured ;
  KEYSET iRuns = gxsnp->iRuns ;
  int riMax = keySetMax (iRuns) ; 
  KEYSET jRuns = gxsnp->jRuns ;
  int rjMax = keySetMax (jRuns) ; 
  int rMax = riMax ;
  int rj, s, k, mx ;
  int n, n0, n1, n0M ;
  int layer, layersH[rjMax+1], layersM[rjMax+1], nnn1 = 0  ;
  SNP *snp ;
  BigArray snps = gxsnp->snps ;
  int sMax = dictMax (gxsnp->snpDict) + 1 ;

  rMax = riMax = keySetMax (iRuns) ; 
  rjMax = keySetMax (jRuns) ; 
 
  for (n0 = 0, s = 1, snp = bigArrp (snps, s, SNP) ; s < sMax ; snp++, s++)
    {
      /* count how nany time the snp is covered or is high in the COMPARE specific jRuns */
      snp->nh2 = snp->nm2 = 0 ;
      for (rj = 0 ; rj < rjMax ; rj++)
	{ 
	  int run =   keySet (jRuns, rj) ;
	  int k = s * rMax + run ;
	  if (bit(HH, k)) 
	    { n0++ ; snp->nh2++ ; }
	  if (bit(MM, k)) 
	    { n0++ ; snp->nm2++ ; }
	}
    }
 
  aceOutf (ao, "## %s\t%s\tfile=%s\n", timeShowNow(), dictName (gx->compareDict, iCompare), aceOutFileName (ao)) ;
  aceOut (ao, "## For each run, number of SNPs present in n of these runs, the n=1 are fully specific.\n") ;
  aceOut (ao, "# Run\tTitle") ;
  for (rj = 1 ; rj <= rjMax ; rj++)
    aceOutf (ao, "\tmutant in %d", rj) ;
  aceOutf (ao, "\tTotal number of homozygous SNPs in this run\tSelected number\tTotal selection\t\t") ;
  for (rj = 1 ; rj <= rjMax ; rj++)
    aceOutf (ao, "\tm/m or +/+ positions meassured in %d classes ", rj) ;
  aceOutf (ao, "\tTotal") ;

  for (rj = 0 ; rj < rjMax ; rj++)
    { 
      int run =   keySet (jRuns, rj) ;
      RC *rc = arrp (gx->runs, run, RC) ;
      if (rc->private) continue ;
      aceOutf (ao, "\n%s (%d)", dictName (gx->runDict, run),keySet (gxsnp->iCounts, run) ) ;
      aceOutf (ao, "\t%s",  rc->title ?  stackText (gx->info, rc->title) : dictName (gx->runDict, run)) ;
      if (keySet (gxsnp->iCounts, run) < 10)
	continue ;
      memset (layersH, 0, sizeof(layersH)) ;
      memset (layersM, 0, sizeof(layersM)) ;
      for (n0 = 0, s = 1, snp = bigArrp (snps, s, SNP) ; s < sMax ; s++, snp++)
	{
	  if (snp->nm) 
	    {
	      k = s * rMax + run ; 
	      if (bit(HH, k))
		{
		  layersH[snp->nh2]++ ; /* layers = number of snps seen in 1 to rjMax jRuns */  
		  if (snp->ok)
		    n0++ ;
		}
	      if (bit(MM, k))
		{
		  layersM[snp->nm2]++ ; /* layers = number of snps seen in 1 to rjMax jRuns */  
		}
	    }
	} 
   

      for (n = n0, n1 = 0, layer = 1 ; layer <= rjMax ; layer++)
	  for (s = 1, snp = bigArrp (snps, s, SNP) ; s < sMax && (1 || n < NN || (layer == 1 && n1 < 1000 && 5 * n1 < NN)) ; s++, snp++)
	    {
	      k = s * rMax + run ; 
	      if (snp->nh2 == layer && 
		  /* 100 * snp->nm2 > 80 * rjMax && */
		  bit(HH, k)
		  )
	 	{ n++ ; if (layer == 1) n1++ ; if (! snp->ok) { snp->ok = 1 ; nnn1++ ; }}
	    } 	
   
      for (n0 = 0, k = 1 ; k <= rjMax ; k++)
	{
	  aceOutf (ao, "\t%d", layersH[k]) ;
	  n0 += layersH[k] ;
	} 
      aceOutf (ao, "\t%d\t%d\t%d\t\t", n0, n, nnn1) ;
      for (n0M = 0, k = 1 ; k <= rjMax ; k++)
	{
	  aceOutf (ao, "\t%d", layersM[k]) ; 
	  n0M += layersM[k] ;
	} 
      aceOutf (ao, "\t%d", n0M) ;

      if (0) fprintf (stderr, "..... gxSnpCountSpecificSNPs selected %d SNPs in run %d :: %s   -> total %d\n"
		      , n
		      , rj
		      , dictName (gx->runDict, run)
		      , nnn1
		      ) ; 
    }
  aceOutf (ao, "\n\n\n") ;
  messAllocMaxStatus (&mx) ;  
  fprintf (stderr, "... gxSnpCountSpecific\t%s\tmax memory %d Mb\tselected %d SNPs in the union of all runs\n", timeShowNow(), mx, nnn1) ;
  ac_free (h) ;
} /* gxSnpCountSpecificSNPs */

/*************************************************************************************/
/* count for all runs how many of its SNPs are characteritic of each training runs  */
static void gxSnpCountCharacters (GX *gx,  GXSNP *gxsnp)
{
  BitSet HH = gxsnp->isHigh ;
  BitSet MM = gxsnp->isMeasured ;
  int ki, kj, ri, rj, rij, s ;
  int sMax = dictMax (gxsnp->snpDict) + 1 ;
  KEYSET iRuns = gxsnp->iRuns ;
  int riMax = keySetMax (iRuns) ; 
  KEYSET jRuns = gxsnp->jRuns ;
  int rjMax = keySetMax (jRuns) ; 
  int rMax = riMax ;
  char bi, bj ;
  SNP *snp ;
  BigArray snps = gxsnp->snps ;

  Array shared  ;           /* class characteristic SNPs also high in run     : relevant numerator */
  Array measured  ;       /* class characteristic SNPs also measured in run : relevant denominator */

  shared = gxsnp->sharedMatrix  = arrayHandleCreate (rMax * rMax, int, gxsnp->h) ;
  measured = gxsnp->measuredMatrix  = arrayHandleCreate (rMax * rMax, int, gxsnp->h) ;

  for (s = 1, snp = bigArrp (snps, s, SNP) ; s < sMax ; s++, snp++)
    if (1 || snp->ok)
      for (rj = 0 ; rj < rjMax ; rj++)
	{
	  for (ri = 0 ; ri < riMax ; ri++)
	    {
	      ki = s * rMax + ri ; kj = s * rMax + keySet(jRuns,rj) ;
	      rij = ri * rjMax + rj ;
	      if (bit(MM,ki) && bit(HH,kj)) 
		{ 
		  array(measured, rij, int)++ ;
		  bi = bit(HH,ki) ? 1 : 0 ; bj = bit(HH,kj) ? 1 : 0 ; 
		  if (bi == bj) array(shared, rij, int)++ ;
		}
	    }
	}
} /* gxSnpCountCharacters */

/*************************************************************************************/

static void gxSnpSortByCharacters (GX *gx,  GXSNP *gxsnp)
{  
  CS *cs ;
  RC *rci ; 
  int runi, ri, rj, rij, j, mx ;
  unsigned char *word ;
  Array shared = gxsnp->sharedMatrix  ;           /* class characteristic SNPs also high in run     : relevant numerator */
  Array measured = gxsnp->measuredMatrix  ;   /* class characteristic SNPs also measured in run : relevant denominator */
  KEYSET iRuns = gxsnp->iRuns ;
  int riMax = keySetMax (iRuns) ; 
  KEYSET jRuns = gxsnp->jRuns ;
  int rjMax = keySetMax (jRuns) ; 
  Array bb2 = arrayCreate (riMax, CS) ;
  KEYSET iRunsSorted = gxsnp->iRunsSorted = keySetHandleCreate (gxsnp->h) ;

   for (ri = 1 ; ri < riMax ; ri++)
    {
      Array bb1 ;
   
      bb1 = arrayCreate (rjMax, CS) ;
      runi =  keySet (iRuns,ri) ;
      rci = arrp (gx->runs,runi, RC) ;
      if (rci->private) continue ;
      if (keySet (gxsnp->iCounts, runi) < 10)
	continue ;

      cs = arrayp (bb2, ri, CS) ;
      cs->run = ri ;
      
      for (rj = 0 ; rj < rjMax ; rj++)
	{
	  /* rjj = keySet(jRuns,rj) * rjMax + rj ;   */
	  rij = ri * rjMax + rj ; 
	  cs = arrayp (bb1, rj, CS) ;
	  cs->run = rj ;
	  {
	    float h = array(shared, rij, int) ;
	    float m = array(measured, rij, int) ;

	    if (m < 1)
	      cs->weight = 0 ;
	    else
	      cs->weight = - 1000.0 * h/m ; 
	    if (0) fprintf (stderr, "....ri=%d rj=%d %s\t%s\t h=%.0f\tm=%.0f\ts=%.1f\n"
			    , ri, rj
			    , dictName (gx->runDict, keySet (iRuns, ri))
			    , dictName (gx->runDict, keySet (jRuns, rj))
			    , h, m, cs->weight
			    ) ;
	  }
	}
      arraySort (bb1, csOrder) ;
      word = array(bb2, ri, CS).word ;
      for (j = 0 ; j < rjMax && j < 64 ; j++)
	{
	  cs = arrp (bb1, j, CS) ;
	  rj = cs->run ;
	  if (j == 0)
	    {   /* 2014_08)24: intercale le score du premier avant la liste ordinale des autres */
	      int c, x ;
	      float c1 = - cs->weight / 10 ;
	      if (c1 < 0) c1 = 0 ;
	      if (c1 > 100) c1 = 100 ;
	      c = c1 ;
              x = 2 * rj + 257 ;  
	      word[0] =  (x >> 1) & 0xff ;
	      word[1] =  (x ) & 0xff ;
	      word[2] = 101 - c ;
	    }
	  else
	    {
	      int x = 2 * rj + 257 ; 
	      word[2*j + 3] =  (x >> 1) & 0xff ;
	      word[2*j + 4] =  (x ) & 0xff ;
	    }
	}
      word[2*j+3] = 0 ;
     arrayDestroy (bb1) ;
    } 
  arraySort (bb2, csWordOrder) ;
  
  for (ri = 0, cs = arrp (bb2, 0, CS) ; ri < riMax ; ri++, cs++)
    keySet (iRunsSorted, ri) = cs->run ;

  messAllocMaxStatus (&mx) ;  
  fprintf (stderr, "... gxSnpSort        \t%s\tmax memory %d Mb\tsorted %d values\n", timeShowNow(), mx, arrayMax (bb2)) ;

  arrayDestroy (bb2) ;

  return ;
} /* gxSnpSortByCharacters */

/*************************************************************************************/

static void gxSnpExport (ACEOUT ao, GX *gx, GXSNP *gxsnp, int iCompare)
{ 
  AC_HANDLE h = ac_new_handle () ;
  int ri, ri2, rj, rij, run, pass, mx, nn1 = 0 ;
  KEYSET iCounts = gxsnp->iCounts ;
  KEYSET iRunsSorted = gxsnp->iRunsSorted ; 
  int riMax = keySetMax ( gxsnp->iRuns) ;
  KEYSET iRuns = gxsnp->iRuns ;
  KEYSET jRuns = gxsnp->jRuns ;
  int rjMax = keySetMax (jRuns) ; 
  Array shared = gxsnp->sharedMatrix  ;           /* class characteristic SNPs also high in run     : relevant numerator */
  Array measured = gxsnp->measuredMatrix  ;   /* class characteristic SNPs also measured in run : relevant denominator */

  RC *rc ;

  for (pass = 0 ; pass < 4 ; pass++)
    {
      switch (pass)
	{
	case 0 :
	  aceOutf (ao, "# %s\tPercentage of class characteristic SNPs seen as m/m versus +/+ in each run\t %s\tfile=%s\n", timeShowNow(), dictName (gx->compareDict, iCompare), aceOutFileName (ao)) ;
	  break ;
	case 1:
	  aceOutf (ao, "# %s\tExact fractio  of class characteristic SNPs seen as m/m versus +/+ in each run\t %s\tfile=%s\n", timeShowNow(), dictName (gx->compareDict, iCompare), aceOutFileName (ao)) ;
	  break ;
	case 2:
	  aceOutf (ao, "# %s\tNumber of class characteristic SNPs seen as m/m versus +/+ in each run\t %s\tfile=%s\n", timeShowNow(), dictName (gx->compareDict, iCompare), aceOutFileName (ao)) ;
	  break ;
	case 3:
	  aceOutf (ao, "# %s\tNumber of class characteristic SNPs covered at least 10 times in each run\t %s\tfile=%s\n", timeShowNow(), dictName (gx->compareDict, iCompare), aceOutFileName (ao)) ;
	  break ;
	}
      aceOut (ao, "# Run\tRun") ;
      for (ri = 0 ; ri < riMax ; ri++)
	{
	  ri2 = keySet(iRunsSorted,ri) ;
	  run = keySet(iRuns,ri2) ; 
	  if (!run)
	    continue ;
	  if (keySet(iCounts, run) < 1)
	    continue ; 
	  rc = arrp (gx->runs, run, RC) ;
	  if (rc->private) continue ;
	  if (keySet (gxsnp->iCounts, run) < 10)
	    continue ;

	  aceOutf (ao, "\t%s", dictName (gx->runDict, run)) ;
	}
      aceOut (ao, "\n# Run\tSample") ;
      for (ri = 0 ; ri < riMax ; ri++)
	{
	  ri2 = keySet(iRunsSorted,ri) ;
	  run = keySet(iRuns,ri2) ; 
	  if (!run)
	    continue ;
	  if (keySet(iCounts, run) < 1)
	    continue ;
	  rc = arrp (gx->runs, run, RC) ;
	  if (rc->private) continue ;
	  if (keySet (gxsnp->iCounts, run) < 10)
	    continue ;

	  aceOutf (ao, "\t%s",  rc->sample ?  stackText (gx->info, rc->sample) : dictName (gx->runDict, run)) ;
	}

      aceOut (ao, "\n# Run\tSorting_title") ;
      for (ri = 0 ; ri < riMax ; ri++)
	{
	  ri2 = keySet(iRunsSorted,ri) ;
	  run = keySet(iRuns,ri2) ; 
	  if (!run)
	    continue ;
	  if (keySet(iCounts, run) < 1)
	    continue ;
	  rc = arrp (gx->runs, run, RC) ;
	  if (rc->private) continue ;
	  if (keySet (gxsnp->iCounts, run) < 10)
	    continue ;

	  aceOutf (ao, "\t%s",  rc->sortingTitle ?  stackText (gx->info, rc->sortingTitle) : dictName (gx->runDict, run)) ;
	}
   
      aceOut (ao, "\n# Run\tSorting_title_2") ;
      for (ri = 0 ; ri < riMax ; ri++)
	{
	  ri2 = keySet(iRunsSorted,ri) ;
	  run = keySet(iRuns,ri2) ; 
	  if (!run)
	    continue ;
	  if (keySet(iCounts, run) < 1)
	    continue ;
	  rc = arrp (gx->runs, run, RC) ;
	  if (rc->private) continue ;
	  if (keySet (gxsnp->iCounts, run) < 10)
	    continue ;

	  aceOutf (ao, "\t%s",  rc->sortingTitle2 ?  stackText (gx->info, rc->sortingTitle2) : dictName (gx->runDict, run)) ;
	}
   
      aceOut (ao, "\n# Run\tOther_title") ;
      for (ri = 0 ; ri < riMax ; ri++)
	{
	  ri2 = keySet(iRunsSorted,ri) ;
	  run = keySet(iRuns,ri2) ; 
	  if (!run)
	    continue ;
	  if (keySet(iCounts, run) < 1)
	    continue ;
	  rc = arrp (gx->runs, run, RC) ;
	  if (rc->private) continue ;
	  if (keySet (gxsnp->iCounts, run) < 10)
	    continue ;

	  aceOutf (ao, "\t%s",  rc->otherTitle ?  stackText (gx->info, rc->otherTitle) : dictName (gx->runDict, run)) ;
	}
       
      aceOut (ao, "\n# Run\tTitle") ;
      for (ri = 0 ; ri < riMax ; ri++)
	{
	  ri2 = keySet(iRunsSorted,ri) ;
	  run = keySet(iRuns,ri2) ; 
	  if (!run)
	    continue ;
	  if (keySet(iCounts, run) < 1)
	    continue ;
	  rc = arrp (gx->runs, run, RC) ;
	  if (rc->private) continue ;
	  if (keySet (gxsnp->iCounts, run) < 10)
	    continue ;

	  aceOutf (ao, "\t%s",  rc->title ?  stackText (gx->info, rc->title) : dictName (gx->runDict, run)) ;
	}
         
      for (rj = 0 ; rj < rjMax ; rj++)
	{
	  int runi, runj =   keySet (jRuns, rj) ;
	  rc = arrp (gx->runs, runj, RC) ;
	  if (!runj)
	    continue ;
	  if (keySet(iCounts, runj) < 10)
	    continue ;
	  if (rc->private) continue ;
	  aceOutf (ao, "\n%s", dictName (gx->runDict, runj)) ;
	  aceOutf (ao, "\t%s",  rc->title ?  stackText (gx->info, rc->title) : dictName (gx->runDict, runj)) ;
	  
	  for (ri = 0 ; ri < riMax ; ri++)
	    {
	      float m, h ;
	      RC *rci ;

	      ri2 = keySet(iRunsSorted,ri) ; 
	      runi = keySet(iRuns,ri2) ; 
	      if (!runi)
		continue ;
	      if (keySet (iCounts, runi) < 10)
		continue ;
	      rci = arrp (gx->runs, runi, RC) ;
	      if (rci->private) continue ;

	      rij = ri2 * rjMax + rj ;
	      m = array(measured, rij, int) ; 
	      h = array(shared, rij, int) ; 
	      if (m < 1) m = 1 ;
	      switch (pass)
		{
		case 0:
		  aceOutf (ao, "\t%.6f", m > 0 ? 100.0 * h/m : 0) ;
		  nn1++ ;
		  break ;
		case 1:
		  aceOutf (ao, "\t%.0f:%.0f",  h, m) ;
		  break ;
		case 2:
		  aceOutf (ao, "\t%.0f",  h) ;
		  break ;
		case 3:
		  aceOutf (ao, "\t%.0f",  m) ;
		  break ;
		}
	    }
	}
      aceOutf (ao, "\n") ;
      aceOutf (ao, "\n") ;
      aceOutf (ao, "\n") ;
    }  
  messAllocMaxStatus (&mx) ;  
  fprintf (stderr, "... gxSnpExport      \t%s\tmax memory %d Mb\texported %d values in %s\n", timeShowNow(), mx,  nn1, aceOutFileName (ao)) ;
  ac_free (h) ;
} /* gxSnpExport */

/*************************************************************************************/

static void gxTestTroisPoints (ACEOUT ao, GX *gx, GXSNP *gxsnp, int iCompare)
{   
  typedef struct  { int run1, run2, run3,  nn, n[8] ;} T3P ;
  AC_HANDLE h = ac_new_handle () ;
   BitSet HH = gxsnp->isHigh ;
   BitSet MM = gxsnp->isMeasured ;
   BigArray snps = gxsnp->snps ;
   int sMax = keySetMax (snps) ; 
   KEYSET iRuns = gxsnp->iRuns ;
   int riMax = keySetMax (iRuns) ; 
   KEYSET jRuns = gxsnp->jRuns ;
   int rjMax = keySetMax (jRuns) ; 
   int rMax = riMax ;
   int mx, t3Max = 0 ;
   int s, rj1, rj2, rj3, kj1, kj2, kj3, run1, run2, run3, b1, b2, b3 ;
   SNP *snp ;
   RC *rc ;
   T3P *t3p ;
   Array t3ps = arrayHandleCreate (rjMax * rjMax * rjMax, T3P, h) ;

   for (s = 1, snp = bigArrp (snps, s, SNP) ; s < sMax ; s++, snp++)
     {
       if (! snp->ok) continue ;

       for (rj1 = 0 ; rj1 < rjMax ; rj1++)
	 {
	   run1 =  keySet (jRuns,rj1) ;
	   rc = arrp (gx->runs,run1, RC) ;
	   if (rc->private) continue ;
	   if (keySet (gxsnp->iCounts, run1) < 10)
	     continue ;
	   kj1 = s * rMax + keySet(jRuns,rj1) ;
	   if (! bit(MM,kj1))continue ;
	   b1 = bit (HH, kj1)  ? 1 : 0 ; 

	   for (rj2 = rj1 + 1 ; rj2 < rjMax ; rj2++)
	     {
	       run2 =  keySet (jRuns,rj2) ;
	       rc = arrp (gx->runs,run2, RC) ;
	       if (rc->private) continue ;
	       if (keySet (gxsnp->iCounts, run2) < 10)
		 continue ;
	       kj2 = s * rMax + keySet(jRuns,rj2) ;
	       if (! bit(MM,kj2))continue ;
	       b2 = bit (HH, kj2)  ? 1 : 0 ; 

	       for (rj3 = rj2 + 1 ; rj3 < rjMax ; rj3++)
		 {
		   run3 =  keySet (jRuns,rj3) ;
		   rc = arrp (gx->runs,run3, RC) ;
		   if (rc->private) continue ;
		   if (keySet (gxsnp->iCounts, run3) < 10)
		     continue ;
		   kj3 = s * rMax + keySet(jRuns,rj3) ;
		   if (! bit(MM,kj3))continue ;
		   
		   t3p = arrayp (t3ps, rj1 * rjMax * rjMax + rj2 * rjMax + rj3, T3P) ;
		   t3p->run1 = run1 ; t3p->run2 = run2 ; t3p->run3 = run3 ; 
		   t3p->nn++ ;
		   b3 = bit (HH, kj3) ? 1 : 0 ; 
		   t3p->n[b1 + 2*b2 + 4*b3]++ ;
		 }
	     }
	 }
     }
   
   aceOutf (ao, "\n\n\nThree points tests: for each triplet of runs, given in columns 1, 2, 3, select all SNPs measured in all 3 runs and report their genotype as a Venn diagram\n") ;
   aceOutf (ao, 
	    "In the A-B-C test, the run  with the highest number of specific alleles, say +mm or m++, is most distant from the 2 others\n"
	    "the relevant sums are given in the columns called A,B,C special, for example B special = A ~B C    +   ~A B ~C\n"
	    "Columns 1,2,3 give the run name, 4,5,6 give the special percentages, 7-15 give the absolute counts for the Venn diagram, 16-24 the relative Venn diagram,\n"
	    "The denominator, use to compute the percentages is the total number of SNPs measured in all 3 runs, the sum of the Venn diagram\n"
	    "These triplets allow to construct recursivelly the philogenic tree following the simple rule that if A is distant from the BC pair\n"
	    "A and BC sjould sit on different branches of the full tree.\n"
	    "The tree is consistant if it respects all triplets, or at least all triplets with high contrast.\n"
	    ) ;
   
   t3Max = arrayMax (t3ps) ;
   aceOutf (ao, "A          \tB          \tC          \t\tA++ or +BC\t+B+ or A+C\t\"++C\" or AB+\t\tTotal\t\t+++\t\tA++\t+B+\t\"++C\"\t\tAB+\tA+C\t+BC\t\t\tABC\t\t\t\t+++\t\tA++\t+B+\t\"++C\"\t\tAB+\tA+C\t\"+BC\"\t\tABC") ;
   for (s = 0, t3p = arrp (t3ps, 0, T3P) ; s < t3Max ; t3p++, s++)
     {
       double z = t3p->nn/100.0 ;
       RC *rc ;

       if (t3p->nn < 1) continue ;
       rc = arrp (gx->runs, t3p->run1, RC) ;
       if (rc->sortingTitle)
	 aceOutf (ao, "\n%s", stackText (gx->info, rc->sortingTitle)) ;
       else
	 aceOutf (ao, "\n%s", dictName (gx->runDict, t3p->run1)) ;

        rc = arrp (gx->runs, t3p->run2, RC) ;
       if (rc->sortingTitle)
	 aceOutf (ao, "\t%s", stackText (gx->info, rc->sortingTitle)) ;
       else
	 aceOutf (ao, "\t%s", dictName (gx->runDict, t3p->run2)) ;

       rc = arrp (gx->runs, t3p->run3, RC) ;
       if (rc->sortingTitle)
	 aceOutf (ao, "\t%s", stackText (gx->info, rc->sortingTitle)) ;
       else
	 aceOutf (ao, "\t%s", dictName (gx->runDict, t3p->run3)) ;

       aceOutf (ao, "\t\t%.2f\t%.2f\t%.2f"
		, t3p->n[1]/z + t3p->n[6]/z, t3p->n[2]/z + t3p->n[5]/z, t3p->n[4]/z + t3p->n[3]/z
		) ;

       aceOutf (ao, "\t\t%d\t\t%d\t\t%d\t%d\t%d\t\t%d\t%d\t%d\t\t%d"
		, t3p->nn
		, t3p->n[0]
		, t3p->n[1], t3p->n[2], t3p->n[4]
		, t3p->n[3], t3p->n[5], t3p->n[6]
		, t3p->n[7]
		) ;
       aceOutf (ao, "\t\t\t\t%.2f\t\t%.2f\t%.2f\t%.2f\t\t%.2f\t%.2f\t%.2f\t\t%.2f"
		, t3p->n[0]/z
		, t3p->n[1]/z, t3p->n[2]/z, t3p->n[4]/z
		, t3p->n[3]/z, t3p->n[5]/z, t3p->n[6]/z
		, t3p->n[7]/z
		) ;
     }
   aceOutf (ao, "\n") ;

   messAllocMaxStatus (&mx) ;  
   fprintf (stderr, "... gxTestTroisPoints\t%s\tmax memory %d Mb\ttested %d triplets\n", timeShowNow(), mx, t3Max) ;
   ac_free (h) ;
   return ;
}  /* gxTestTroisPoints */

/*************************************************************************************/

static void gxSnpCorrelation (GX *gx) 
{
  COMPARE *compare ;
  int iCompare,  mx ;
  int NN = 1000000 ;
  GXSNP *gxsnp ;

  gxsnp = gxSnpInit (gx) ;

  messAllocMaxStatus (&mx) ;   
  fprintf (stderr, "... gxSnpCorrelation  \t%s\tmax memory %d Mb\tstart\n", timeShowNow(), mx) ;

  gxSnpParse (gx, gxsnp) ;
  
  for (iCompare = 1, compare = arrp (gx->compares, iCompare, COMPARE) ; 
       iCompare <= dictMax(gx->compareDict) ; iCompare++, compare++)
    {
      if (gx->isINTRON && ! compare->compareINTRON)
	continue ;
      if (gx->isSNP && ! compare->compareSNP)
	continue ;
      if (compare->snpProfile && compare->runs && keySetMax (compare->runs) > 0)
	{
	  AC_HANDLE h = ac_new_handle () ;
	  ACEOUT ao = aceOutCreate (gx->outFileName
				    , messprintf (".%s.SNP_profile.txt", dictName (gx->compareDict, iCompare))
				    , FALSE, h) ;
	  aceOutDate (ao, "###", gx->title) ;
	  
	  messAllocMaxStatus (&mx) ;   
	  fprintf (stderr, "... gxSnpCorrelation  \t%s\tmax memory %d Mb\tstart\t%s\n", timeShowNow(), mx, dictName (gx->compareDict, iCompare)) ;
	  
	  {
	    int rj ;
	    keySetMax (gxsnp->jRuns) = 0 ;
	    for (rj = 0 ; rj < keySetMax (compare->runs) ; rj++)
	      keySet (gxsnp->jRuns, rj) = keySet (compare->runs, rj) ;
	  }
	  
	  gxSnpCountSpecificSNPs (ao, gx, gxsnp, iCompare, NN) ;
	  gxSnpCountCharacters (gx, gxsnp)  ;
	  gxSnpSortByCharacters (gx,  gxsnp) ;
	  gxSnpExport (ao, gx, gxsnp, iCompare) ;
	  gxTestTroisPoints (ao, gx, gxsnp, iCompare) ;
	  /* unselect to restore the data as from the parser */
	  gxSnpSelect (gxsnp, TRUE) ;
	  ac_free (h) ;
	}
    }
  messAllocMaxStatus (&mx) ;   
  fprintf (stderr, "... gxSnpCorrelation  \t%s\tmax memory %d Mb\tdone\n", timeShowNow(), mx) ;

  ac_free (gxsnp->h) ;
  ac_free (gxsnp) ;
  
  messAllocMaxStatus (&mx) ;   
  fprintf (stderr, "... gxSnpCorrelation  \t%s\tmax memory %d Mb\tafter freeing gxSnp\n", timeShowNow(), mx) ;
} /*  gxSnpCorrelation */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
  fprintf  (stderr,
	    "// geneindex: analyse gene expression index tables\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, March 2012, mieg@ncbi.nlm.nih.gov\n"
	    "// Purpose\n"
	    "// Given a table of read counts, compute the run index, cluster in groups\n"
	    "// Analyse the table to create heat map, PCA, clusters of various sorts\n"
	    "//\n"
	    "// Syntax:\n"
	    "// geneindex [options]\n"
	    "// The options are all optional and may be specified in any order\n"
	    "//   -o output_file_prefix\n"
	    "//      all exported files will be named output_file_prefix.action"
            "//   -gzo : gzip all output files\n"
	    "//   -gzi : gunzip stdin input\n"
	    "//      all input files called .gz are automatically gunzipped\n"
	    "//   -export [aitmbv] : export expression tables one line per gene\n"
	    "//      a: gene ace file, including groups and indexes computed on the fly\n"
	    "//      i: index, t:tags, m:tags per M raw reds, b:kilobp, v:group variances\n"
	    "//   -malus int\n"
	    "//      lower all expression index by this constant, useful to adjust microarrays\n"
	    "//   -correlation -minIndex <int> : \n"
	    "//   -covariance -minIndex <int> : \n"
	    "//      export a table of run comparison by correlation\n"
	    "//      in the scalar product values below minIndex are raised to minIndex\n"
	    "//   -compare : export differential for all run->Compare_to found in -runAce file\n"
	    "//     -showAnyRun : in compare table show all runs, not just the relevant ones\n"
	    "//     -minFoldChange : default .4 do not consider differentail genes below this limit\n"
	    "//     -exportDiffGenes : export all genes with at least minFoldChange\n"
	    "//   -geneClusters fileName : consider 'meta' gene for example clusters on chromosome regions\n"
	    "//       the file has 3 columns chrom coord coord, all genes in that region make the group\n"
	    "//   -runList\n"
	    "//      list of run names, impose the order of the runs in the outputs, reject other runs\n" 
	    "//   -runAce aceFileName\n"
	    "//        optional ace file, gives genes info, and  groups, sample, machine, title info of each runs\n"
	    "//      -mask fileName : mask out genes with names listed in this file\n"
	    "//      -chromAlias file : rename the chromosome in the ace fileI\n"
	    "//   -deepGene file_name [-u | -nu] [-stableGenes file] [-mask file] [-keepIndex]\n"
	    "//      acefile with gene->deep index and counts\n"
	    "//      -u [default] or -nu selects Run_U or Run_nU counts in the ace file\n"
	    "//      -keepIndex accept the index provided in the file, rather than compute from the counts\n"
	    "//      -stableGenes fileName : fine tune the index by linear regression around this list\n"
	    "//        format 2 words per line: chrom-behind-IntMap-tag  new-chrom-name\n"
	    "//   -deepTranscript file_name\n"
	    "//      acefile with Transcript->deep  counts\n"
	    "//   -deepSNP file_name\n"
	    "//      snp table 20 columns, changes normalize default to FALSE\n"
	    "//     -selectSNP filename\n"
	    "//      keep only SNPs listed in that single column file like: PARN.bAug10:34G2T\n"
	    "//   -snpEval\n"
	    "//      Export all SNPs in -deepSNP file, differential for any pair of runs or groups\n"
	    "//      declared in Compare objects with tag Apply SNP. Export in SNPH/*/$MAGIC.snp.sorted.differential\n"
	    "//   -snpCompare\n"
	    "//      Reexport all SNPs in -deepSNP file, differential for any pair of runs\n"
	    "//      declared in a compare object with tag Apply SNP\n"
	    "//   -snpMinWhoScore int [default 20]\n"
	    "//      min score dispalyed in who_is_who.cluster table\n"
	    "//   -snpCorrelation file_name\n"
	    "//   -snpCovariance file_name\n"
	    "//      file is in SNPH/.snp format\n"
	    "//   -normalize -ratio_bound <number> -iterate [0-2]  : deep cases\n"
	    "//   -pA | -total : mandatory\n"
	    "//     in pA case, saturate lengths to 3kb, in total to 10kb\n"
	    "//   -MA : micro-array data\n"
	    "//   -microRNA : micro-RNA data which was not aligned but directly counted as genes\n"
	    "//   -digital <int> : default 1, weighting of each gene when computing cosine of 2 runs\n"
	    "//      0: use the gene index, 1: +-score, 2: fully digital, +-1 according to index zone, 3:mobymatic, 4:SNP by classes\n"
	    "//   -selectMarker filter\n"
	    "//     keep only genes|transcripts|snps (according to -deepGene|SNP...) whose name are listed\n"
	    "//   -rejectMarker filter\n"
	    "//     reject genes|transcripts|snps (according to -deepGene|SNP...) whose name are listed\n"
	    "//   -skipEmptyGenes : skip genes or mRNAS with zero count\n"
	    "//   -captured name:  only export genes with the tag Capture <name>\n"
	    "//     the default is to export all genes, i.e. present in the target fasta file used by the aligner\n"
	    "//   -targeted :  only export genes with the tag targeted\n"
	    "//     the default is to export all genes, i.e. present in the target fasta file used by the aligner\n"
	    "//   -minVariance : only consider genes with some minimal variance in the largest group\n"
	    "//   -genePlus fp -geneMinus fm\n"
	    "//     Compare all groups declared as  compare_to in -runAce, using these 2 lists of genes\n"
	    "//     fp and fm files list genes, one per line, named as in the -deepGene file [or -deepIntron]\n"
	    "//   -seedGene gene_name\n"
	    "//     construct a pair of groups with the 1/5 most extre expressots\n"
	    "//   -seaLevel int -removeLimit int : modify the normalization of the index\n"
	    "//      seaLevel [default 3] : do not call genes low relative to intergenic level\n"
	    "//      removeLimit [default 50] : reads in genes gathering over 1/limit of all reads\n"
	    "//        do not contribute to the denominator in the normalization of the gene index\n"
	    "//   -referenceGenome text : name of genome release i.e. GRCh38\n"
	    "//   -genomeLengthInKb float : Length of the genome in kilobases, default : 3 10^6\n"
	    "//   -subsample int : desired max number of thousands reads per run\n"
	    "//      In each gene random subsample to keep a total of n kilo aligned read per run\n"
	    "// Caveat:\n"
	    "//   Lines starting with '#' are considered as comments and dropped out\n"
	    ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  geneindex -help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/

int main (int argx, const char **argv)
{
  GX gx ;
  AC_HANDLE h = 0 ;

  FIX = log(1000.0)/log(2.0) ;  /*substract from the magic index to convert back to RPKM */
  FIX = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&gx, 0, sizeof (GX)) ;
  gx.h = h ;
  getCmdLineOption (&argx, argv, "-o", &gx.outFileName) ;


  if (gx.outFileName && gx.outFileName[0] == 'X')
    {
      wilcoxon (-99,-99,-99,0,0) ;
      exit (0) ;
    }

  if (gx.outFileName)
    {
      /* the existence of the .done file proves that the code worked correctly */
      AC_HANDLE h1 = ac_new_handle () ;
      ACEOUT ao = aceOutCreate (gx.outFileName, ".done", FALSE, h1) ;
      const char *cp ;
      if (!ao)
	messcrash ("cannot create file %s.done", gx.outFileName) ;
      cp = aceOutFileName (ao) ;
      ac_free (h1) ;
      if (unlink (cp) == -1)
	messcrash ("Failed to remove check file %s (%s)",
		  cp, messSysErrorText());
    }

  gx.geneSelectionMethod = 7 ; /* 7: histo with gray zone,  8: Gini */ 

  /* defaults */
  gx.digital = 1 ;    
       /* 
	* 0: use the index, 
	* 1: +-score, 
	* 2: +-1 (fully digital) *
	* 3:mobymatic, 
	* 4: SNP count by classes , 
	* 5: like 0 with threshold, (not good for NB predictions) 
	*/
  gx.digitalCorrelation = 0 ;  /* we could try 5, but this is not yet tested */
  gx.normalize = TRUE ;
  gx.localIteration = FALSE ;

  gx.maxGenePlus = 5000 ;
  gx.maxGenePlus = 60 ;


  gx.ratio_bound = 1.0 ; /* was 4,  2.5 ; */
  gx.histo_shifting = 1 ;  /* unit is dixieme de point d'index */
  
  gx.minFoldChange = .4 ;  /* seems best  was 1.0   dromadaire */
  gx.pA = 2 ;
  if (getCmdLineBool (&argx, argv, "-pA"))
    gx.pA = 2 ;
  if (getCmdLineBool (&argx, argv, "-total"))
    gx.pA = 1 ;

  
  /* MICRO ARRAY PARAMETERS */
  if (getCmdLineInt (&argx, argv, "-MA", &(gx.isMA)))
    {
      if (1)
	{
	  gx.histo_shifting = 2 ; 
	  gx.digital = 0 ;    /* 0: use the plain index, 1: +-score, 2: +-1 (fully digital) */
	  gx.normalize = TRUE ;
	}
      gx.pA = 1 ; /* divide by 10 the evaluation of the genomic contamination */
    }
  /* DGE : like in Helicos Digital Gene Expression (only capture the 3prime end or a fixed position 
   *    usually set at the Run level in rc->DGE
   *    globally set in deepSNP and deepIntron mode
   */
  if ((gx.DGE = getCmdLineBool (&argx, argv, "-DGE")))
    {
      if (1)
	{
	  gx.histo_shifting = 2 ; 
	  gx.digital = 0 ;    /* 0: use the plain index, 1: +-score, 2: +-1 (fully digital) */
	  gx.normalize = TRUE ;
	}
      gx.pA = 1 ;
    }
  
  getCmdLineInt (&argx, argv, "-local", &gx.localIteration);
  /* optional arguments */

  if (argx == 1 ||
      getCmdLineBool (&argx, argv, "-h") ||
      getCmdLineBool (&argx, argv, "-help") ||
      getCmdLineBool (&argx, argv, "--help")
      )
    usage (0) ;

  /* input output redirections */
  gx.gzi = getCmdLineBool (&argx, argv, "-gzi");
  gx.gzo = getCmdLineBool (&argx, argv, "-gzo");
  gx.unique = TRUE ;
  if (getCmdLineBool (&argx, argv, "-nu")) 
    gx.unique = FALSE ;
  if (getCmdLineBool (&argx, argv, "-u")) 
    {
      if (! gx.unique)
	messcrash ("Sorry options -u and -nu are incompatible, you cannot set both at the same time") ;
      gx.unique = TRUE ;
    }

  gx.isMicroRNA = getCmdLineBool (&argx, argv, "-microRNA");
  gx.phenoGene = getCmdLineBool (&argx, argv, "-phenoGene");
  gx.fineTune = getCmdLineBool (&argx, argv, "-fineTune");
  gx.iterate = getCmdLineBool (&argx, argv, "-iterate");
  gx.showAnyRun = getCmdLineBool (&argx, argv, "-showAnyRun");
  gx.html = getCmdLineBool (&argx, argv, "-html");
  gx.keepIndex = getCmdLineBool (&argx, argv, "-keepIndex");
  gx.medianCenter = getCmdLineBool (&argx, argv, "-medianCenter");
  gx.noHisto = getCmdLineBool (&argx, argv, "-noHisto"); 
  gx.skipEmptyGenes = getCmdLineBool (&argx, argv, "-skipEmptyGenes"); 
  gx.exportDiffGenes = getCmdLineBool (&argx, argv, "-exportDiffGenes");
  gx.TGx = getCmdLineBool (&argx, argv, "-TGx");

  gx.subsample = 0 ;
  if (1)
    {
      const char *ssBuf = 0 ;
      if (getCmdLineOption (&argx, argv, "-subsample", &ssBuf))
	{
	  char cc = 0 ;
	  sscanf (ssBuf, "%d%c", &gx.subsample, &cc) ;
	  switch ((int)cc)
	    {
	    case 'k' : break ;
	    case 'M' : gx.subsample *= 1000 ; break ;
	    case 'G' : gx.subsample *= 1000000 ; break ;
	    default: usage ("-subsample xx requires <int>[kMG]\n") ; break ;
	    }
	}
      fprintf (stderr , "sub=%d\n", gx.subsample) ;
    }
  if (gx.subsample && gx.keepIndex)
    messcrash ("Sorry options -subsample a nd -keepIndex are incompatible") ;

  gx.targeted = getCmdLineBool (&argx, argv, "-targeted");
  getCmdLineOption (&argx, argv, "-captured", &gx.captured);

  gx.intronGeneMix = getCmdLineBool (&argx, argv, "-intronGeneMix");
  getCmdLineOption (&argx, argv, "-deepGene", &gx.deepFileName) ;
  gx.deepTsf = getCmdLineBool (&argx, argv, "-deepTsf") ;
  getCmdLineOption (&argx, argv, "-snpCorrelation", &gx.snpFileName) ;
  getCmdLineOption (&argx, argv, "-snpCovariance", &gx.snpFileName) ;
  if (getCmdLineOption (&argx, argv, "-deepIntron", &gx.deepFileName))
    {
      gx.isINTRON = gx.hasGeneChrom = TRUE ; gx.DGE = TRUE ; gx.pA = 1 ; gx.digital = 1 ;
      gx.skipEmptyGenes = TRUE ;
    }
  if (getCmdLineOption (&argx, argv, "-deepTranscript", &gx.deepFileName))
    gx.isTranscript = TRUE ;
  if (getCmdLineOption (&argx, argv, "-deepSNP", &gx.deepFileName))
    { 
      gx.isSNP = TRUE ; gx.keepIndex = TRUE ; gx.DGE = TRUE ; gx.pA = 1 ; gx.digital = 1 ;
      gx.digitalCorrelation = 4 ;  /* we could try 5, but this is not yet tested */ 
      gx.skipEmptyGenes = TRUE ;
    }
  /* optionally modified */
  getCmdLineInt (&argx, argv, "-digital", &gx.digital);
  if (gx.digital < 0 || gx.digital > 5)
    usage (messprintf ("Error: parameter -digital %d out of range, the value should be 0,1,2,3 or 4", gx.digital)) ;

  getCmdLineOption (&argx, argv, "-selectMarker", &gx.markerSelectFileName) ;
  getCmdLineOption (&argx, argv, "-rejectMarker", &gx.markerRejectFileName) ;
  gx.snpEval = getCmdLineBool (&argx, argv, "-snpEval") ;
  if (gx.snpEval)
    gx.normalize = FALSE ;
  gx.snpCompare = getCmdLineBool (&argx, argv, "-snpCompare") ;
  if (gx.snpCompare)
    gx.normalize = FALSE ;

  if (getCmdLineOption (&argx, argv, "-htmlSpecies", &gx.htmlSpecies))
    {
      if (! strcasecmp (gx.htmlSpecies,"hs"))
	gx.htmlSpecies = "human" ;
      if (! strcasecmp (gx.htmlSpecies,"ec"))
	gx.htmlSpecies = "E.coli" ;
    }
  
  getCmdLineInt (&argx, argv, "-normalize", &gx.normalize);

  gx.referenceGenome = "NCBI_37" ;
  getCmdLineOption (&argx, argv, "-referenceGenome", &gx.referenceGenome) ;
  getCmdLineOption (&argx, argv, "-runList", &gx.runListFileName) ;
  getCmdLineOption (&argx, argv, "-runAce", &gx.runAceFileName) ;
  getCmdLineOption (&argx, argv, "-mask", &gx.maskFileName) ;
  getCmdLineOption (&argx, argv, "-chromAlias", &gx.chromAliasFileName) ;
  getCmdLineOption (&argx, argv, "-stableGenes", &gx.stableGenesFileName) ;
  getCmdLineOption (&argx, argv, "-geneGroup", &gx.geneGroupFileName) ;
  getCmdLineOption (&argx, argv, "-method", &gx.method) ;
  getCmdLineOption (&argx, argv, "-target_class", &gx.target_class) ;
  getCmdLineOption (&argx, argv, "-export", &gx.export) ;
  if (gx.export && strchr (gx.export, 'z')) gx.anti = TRUE ; /* parse the Anti_run */
  getCmdLineOption (&argx, argv, "-seedGene", &gx.seedGene) ;
  getCmdLineOption (&argx, argv, "-genePlus", &gx.genePlusFileName) ;
  getCmdLineOption (&argx, argv, "-geneMinus", &gx.geneMinusFileName) ;
  gx.minIndex = 0 ; getCmdLineInt (&argx, argv, "-minIndex", &gx.minIndex) ;
  gx.seaLevel = 3 ; getCmdLineInt (&argx, argv, "-seaLevel", &gx.seaLevel) ;
  gx.removeLimit = 50 ;  getCmdLineInt (&argx, argv, "-removeLimit", &gx.removeLimit) ;
  gx.wall = 4 ;  getCmdLineInt (&argx, argv, "-wall", &gx.wall) ;

  getCmdLineInt (&argx, argv, "-malus", &gx.malus) ; 
  getCmdLineInt (&argx, argv, "-maxGene", &gx.maxGenePlus) ; 
  getCmdLineFloat (&argx, argv, "-ratio_bound", &gx.ratio_bound) ; 
  getCmdLineFloat (&argx, argv, "-minVariance", &gx.minVariance) ; 
  getCmdLineFloat (&argx, argv, "-threshold" , &gx.threshold) ;
  getCmdLineFloat (&argx, argv, "-minFoldChange" , &gx.minFoldChange) ;
  gx.minWhoScore = 5 ;
  getCmdLineFloat (&argx, argv, "-snpMinWhoScore" , &gx.minWhoScore) ;
  getCmdLineInt (&argx, argv, "-histo_shifting", &gx.histo_shifting) ; 
  getCmdLineOption (&argx, argv, "-title", &gx.title) ;
  getCmdLineOption (&argx, argv, "-geneClusters", &gx.geneClustersFileName) ;
  getCmdLineOption (&argx, argv, "-split_mRNAs", &(gx.splitMrnaFileName)) ;
  gx.correlation = 
    getCmdLineBool (&argx, argv, "-correlation") || 
    getCmdLineBool (&argx, argv, "-covariance") ;
  gx.compare_to = getCmdLineBool (&argx, argv, "-compare") ;
  gx.isMRNAH = getCmdLineBool (&argx, argv, "-MRNAH") ;

  gx.genomeLengthInKb = 3000000 ;
  getCmdLineFloat (&argx, argv, "genomeLengthInKb", &gx.genomeLengthInKb) ;
  if (getCmdLineBool (&argx, argv, "-test"))
    {
      exit(0) ;
    }
  if ((gx.genePlusFileName || gx.geneMinusFileName) && ! (gx.genePlusFileName && gx.geneMinusFileName))
    usage ("-genePlus gene_list -geneMinus gene_list parameters come in pairs, sorry") ;
  if (gx.genePlusFileName)
    {
      if (! filName (gx.genePlusFileName, 0, "r"))
	usage (messprintf ("-genePlus expects a file name with one gene per line: Sorry, I cannot open file %s", gx.genePlusFileName)) ;
      if (! filName (gx.geneMinusFileName, 0, "r"))
	usage (messprintf ("-geneMinus expects a file name with one gene per line: Sorry, I cannot open file %s", gx.geneMinusFileName)) ;
    }

  if (argx > 1)
    usage (messprintf ("COMMAND LINE ERROR: Unknown argument %s,", argv[1])) ;

  fprintf (stderr, "// geneindex start: %s\n", timeShowNow()) ;
  /* Input is an ace file Run->Deep */

  gxParseInit (&gx) ;

  if (gx.runListFileName)     /* impose the order in which they will be reported */
    gxParseRunList (&gx, gx.runListFileName) ;
  if (gx.splitMrnaFileName)
      gxSplitMrnaParse (&gx) ;
  if (gx.maskFileName)
    gxParseMask (&gx, gx.maskFileName) ;
  if (gx.chromAliasFileName)
    gxParseChromAlias (&gx, gx.chromAliasFileName) ;
  if (gx.runAceFileName)
    gxAceParse (&gx, gx.runAceFileName,TRUE) ;


  gx.isLowTitle[0] = "" ;
  gx.isLowTitle[1] = "NA/" ;
  gx.isLowTitle[2] = "NE/" ;


  if (gx.snpEval)
    {
      int compare ;
      COMPARE *up ;

      for (compare = 0 ; compare < arrayMax (gx.compares) ; compare++)
	{
	  up = arrayp (gx.compares, compare, COMPARE) ;
	  if (up->compareSNP && up->runs)
	    {
	      int ii, iMax = arrayMax (up->runs) ;
	      for (ii = 0 ; ii < iMax ; ii++)
		{
		  int run1 = keySet (up->runs, ii) ;
		  RC *rc1 = arrayp (gx.runs, run1, RC) ;
		  rc1->snpCompare = TRUE ; /* will help gxOneGroupCumul on next snp */
		}
	    }
	}
    }

  if (gx.genePlusFileName)
    {
      gx.userGivenGenePlus = gxParsePlusMinusFile (&gx, TRUE) ;
      gx.userGivenGeneMinus = gxParsePlusMinusFile (&gx, FALSE) ;
      gx.userGivenGenePlusMinus = keySetOR (gx.userGivenGenePlus, gx.userGivenGeneMinus) ;
    }

  if (gx.snpFileName)
    {
      gxSnpCorrelation (&gx) ;
    }
  else if (gx.snpEval || gx.snpCompare)
    {
      gxAceParse (&gx, gx.deepFileName,FALSE) ;
    }
  else if (gx.deepFileName)
    {

      if (gx.geneGroupFileName)
	gxParseGeneGroup (&gx, gx.geneGroupFileName) ;
      if (gx.deepTsf)
	{
	  if (! gxDeepTsfParse (&gx, gx.deepFileName,FALSE))
	    messcrash ("Cannot parse the .tsf file %s given as intput, expecting run->Deep counts\n"
		   , gx.deepFileName ? gx.deepFileName : "stdin" ) ;
	}
      else
	{
	  if (! gxAceParse (&gx, gx.deepFileName,FALSE))
	    messcrash ("Cannot parse the .ace file %s given as intput, expecting run->Deep counts\n"
		       , gx.deepFileName ? gx.deepFileName : "stdin" ) ;
	}
      if (gx.subsample)
	gxSubsample (&gx) ;
      if (gx.geneClustersFileName)
	gxPrepareGeneClusters (&gx) ;
      /*  gxHackZeroCounts (&gx) ; */
      gxGlobalCounts (&gx) ;


      if (gx.geneGroupFileName)
	gxGeneGroupCount (&gx) ;

      gxComputeAllIndex (&gx) ; 
      gxComputeAverageIndex (&gx) ; 
      if (gx.geneClustersFileName)
	{ 
	  gxComputeGeneClusterIndex (&gx) ; 
	  gxExportGeneClusterIndex  (&gx) ; 
	}
      if (0) gxRankGene (&gx) ; /* use the rank of the gene as a definition of the index */
      if (0 && gx.fineTune)    /* dromadaire  adjusting across heterogeneous samples seems very detrimental 
				* it is also likely that in gxFine tune the call to computeAllIndex is messing 
				* with the index of the additive and non additive groups
				*/
	gxFineTune (&gx) ; /* fine tune the global linear correlation of all runs */
      if (gx.genePlusFileName)
	{
	  goto done ;
	}
      /* crucial: sets NA_crossover_index and a minimal high gene expression */
      gxNumberOfGenesPerRunAboveGivenIndex (&gx, FALSE) ;
      if (dictMax(gx.compareDict) > 0)
	{
	  gx.degHistos = arrayHandleCreate (100, DGH, gx.h) ;
	  gx.pairs2cmp = arrayHandleCreate (1000, COMPARE, gx.h) ;
	  gx.pairDict = dictHandleCreate (1000, gx.h) ;
	  gx.pairScoresDict = dictHandleCreate (1000000, gx.h) ;
	  gx.pairScores = arrayHandleCreate (1000000, float, gx.h) ;
	  gx.pValues = arrayHandleCreate (1000000, float, gx.h) ;
	  gx.pairScores2compare = keySetHandleCreate (gx.h) ;
	  gxExpressionProfile (&gx, 0) ;  /* CreatePrivateCompare */
	  for (int i = 1 ; i <= dictMax (gx.compareDict) ; i++)
	    dictAdd (gx.pairDict, dictName (gx.compareDict, i), 0) ;
	}


      if (gx.seedGene) 
	gxSeedGroup (&gx, 1) ;

      if (0) /* this can generate a huge table which take hours */
	{
	  gxAllGeneCorrelations (&gx) ;
	  goto done ;
	}
      
      if (gx.phenoGene)
	gxPhenoGenes (&gx) ;

      gxExportGeneAceFileSummary (&gx, 0) ; /* count genes above given index */
      if (0 && gx.export) gxExportTable (&gx, 10) ; /* header */
      if (gx.export && strchr (gx.export, 'a')) gxExportGeneAceFile (&gx) ; /* ace file */
      if (gx.export && strchr (gx.export, 'i')) gxExportTable (&gx, 0) ;    /* gene index */
      if (gx.export && strchr (gx.export, 't')) gxExportTable (&gx, 1) ; /* tags */
      if (gx.export && strchr (gx.export, 'b')) gxExportTable (&gx, 2) ; /* kb */
      if (gx.export && strchr (gx.export, 'v')) gxExportTable (&gx, 3) ; /* variance */
      if (gx.export && strchr (gx.export, 'i')) gxExportTable (&gx, 4) ; /* sFPKM */
      if (gx.export && strchr (gx.export, 'v') && gx.hasSelectedVariance) gxExportTable (&gx, 5) ; /* variance */
      if (gx.export && strchr (gx.export, 'm')) gxExportTable (&gx, 7) ; /* tags per million raw reads*/

      if (gx.maxGeneGroupLevel && gx.export && strchr (gx.export, 'i')) gxExportTable (&gx, 1000) ; /* gene_group index */
      if (gx.maxGeneGroupLevel && gx.export && strchr (gx.export, 't')) gxExportTable (&gx, 1001) ; /* gene_group index */
      if (gx.maxGeneGroupLevel && gx.export && strchr (gx.export, 'b')) gxExportTable (&gx, 1002) ; /* gene_group index */
      if (gx.maxGeneGroupLevel && gx.export && strchr (gx.export, 'm')) gxExportTable (&gx, 1007) ; /* gene_group index */

      if (0 && gx.export && strchr (gx.export, 'z')) gxExportTable (&gx, 6) ; /* zebre:antisense */

#ifdef JUNK
      if (gx.hasDaa) gxExportTable (&gx, 11) ; /* nReads */
      if (gx.hasDaa) gxExportTable (&gx, 12) ; /* nReadsOk */
      if (gx.hasDaa) gxExportTable (&gx, 13) ; /* nerr */
      if (gx.hasDaa) gxExportTable (&gx, 17) ; /* badtopo */
      if (gx.hasDaa) gxExportTable (&gx, 18) ; /* multi */
      if (gx.hasDaa) gxExportTable (&gx, 19) ; /* multi2 */

      if (gx.hasDaa) gxExportTable (&gx, 14) ; /* a2g */
      if (gx.hasDaa) gxExportTable (&gx, 15) ; /* partial */
      if (gx.hasDaa) gxExportTable (&gx, 16) ; /* orphan */
#endif

      if (gx.export || gx.correlation || gx.compare_to) gxExportTop (&gx, 500) ; /* export the top genes of all runs and groups */
      ac_free (gx.gza) ;

      if (1 && gx.correlation)
	{
	  if (0) gxCorrelation (&gx, 5, TRUE) ;
	  if (0) gxCorrelation (&gx, 8, FALSE) ;
	  if (0) gxCorrelation (&gx, 8, TRUE) ; 

	  if (1) gxCorrelation (&gx, 10, FALSE) ;
	  if (1) gxCorrelation (&gx, 10, TRUE) ;
	  if (1) gxCorrelation (&gx, 12, FALSE) ;
	  if (1) gxCorrelation (&gx, 12, TRUE) ; 
	}
      if (1 && gx.compare_to) 
	{
	  gxCompare_to (&gx) ; 
	  
	  if (gx.showAllHistos)
	    gxShowAllHistos (&gx) ; /* will register */
	  if ((1 || gx.isSNP) && gx.who_is_who)  /* 2019_05_05 remove other outputs */
	    gxSamplePairing (&gx) ;
	  gxExportComparedGenes (&gx) ;
	}
      if (gx.expressionProfile && gx.compare_to)
	gxExpressionProfile (&gx, 1) ;  /* use registered value */
      if (gx.mixing_ratios && gx.compare_to)
	gxMixingRatios (&gx) ;  /* use registered value */

      if (0) 
	gxExportMeanversusVariancePlot (&gx) ;
    }

 done:

  fprintf (stderr, "// gx.digital = %d\n", gx.digital) ;
  {
    /* the existence of the .done file proves that the code worked correctly */
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, "// done: %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
    if (gx.outFileName)
      {
	ACEOUT ao = aceOutCreate (gx.outFileName, ".done", FALSE, h) ;
	aceOutf (ao,  "// done: %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
      }
  }
  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
} /* main */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
