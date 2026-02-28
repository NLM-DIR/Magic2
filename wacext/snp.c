/* snp.c
 * Authors: Jean Thierry-Mieg, NCBI, 
 * 2010:2014
 * SNP search in the output of the MAGIC suite
*/

/* #define ARRAY_CHECK     */

#include "../wac/ac.h"
#include "bitset.h"
#include "keyset.h"
#include "query.h"
#include "dict.h"
#include "aceio.h"
#include <wiggle.h>

/* static int SNPDEBUG=0 ; */
static void usage (char *message) ;

typedef struct hitStruct { 
  int tag , score     /* identifier of the short sequence tag, score of the ali */
    , mult            /* multiplicity of the tag */
    , ln              /* length to be aligned (minus the vectors and pA) */
    , ali             /* length of the alignment (measured on the tag) */
    , x1, x2
    , class           /* the first letter of the target class, defining the hierarchy */
    , target_class    /* entry in targetDict, A_mito ... Z_genome */
    , gene            /* entry in targetDict */
    , unicity         /* number of targets in that class, negative if we change strand */
    , target          /* entry in targetDict */
    , a1, a2, nN, nErr
    , errTag          /* entry in errDict, errors in tag coordinates and orientation */
    , errTarget       /* entry in errDict, errors in target coordinates and orientation */
    , prefix, suffix  /* overhangs, read in the tag, on the strand extending the alignment */
    , targetPrefix, targetSuffix  /* 30bp, up and downstream, read on the target in the orientation of the tag */
    , dPair           /* distance to pair */
    , isRead2         /* 1 if a read< */
    ;
} HIT ;

typedef struct snpStruct {  
  AC_HANDLE h ; BOOL gzi, gzo, smooth ; 
  const char *inFileName, *outFileName, *selectFileName, *title, *fastaFileName, *targetGeneFileName, *newIntronsFileName ;
  const char *remapFileName1, *remapFileName2, *runQualityFileName, *runListFileName ;
  const char *run, *target_class, *snpListFileName, *snpValFileName,*project, *Reference_genome ;
  BOOL mito, qual, BRS_detect, BRS_Detect, BRS_count, pool, ventilate, solid, BRS_merge, mergeExportPopulationTable, aliExtend , phasing ;
  BOOL vcfName ;
  int ooFrequency, editSequence ;
  Stack targetDna ; /* actual dna ofd the targets */
  Array target2dna, targetDnaLength, targetGene ; /* array of offset in targetDna */
  int minCover, minFrequency, minMutant, minAli, minAliPerCent,  Remove_inserts_shorter_than ;
  BOOL hits2BRS ;
  int mult, n, NN, snp ;
  int eeLastSub, eeLastInsert ;
  BOOL unique, db_intersect, remap, db_pheno, db_report, db_VCF ;
  BOOL vcfTable, frequencyTable, prevalenceTable, frequencyHisto, dbCount ;
  BOOL dropMultiplicity ; /* count each sequence only one, dropping the #multiplicity parameter */
  BOOL differential ; /* limit to Variant with tag differential */
  BOOL dropMonomodal ;
  int analyzeEditedHits ;
  int anyRunMeasured, anyRunDropped, anyRunPure, anyRunHigh, anyRunMid, anyRunLow, anyRunWild, anyRunNA ;
  Array target2geneAtlas ;
  Array target2exonAtlas ;
  Array selectZone ;
  Array runQuality ;
  Array aliSnps ;
  Array errCounts ;
  Array errorsP, dnas ;
  BigArray snps ;
  Array aaa, aaaP ;  /* array of bigArrays of BRS struct */
  enum { STRATEGY_ZERO=0, STRATEGY_RNA_SEQ, STRATEGY_EXOME, STRATEGY_GENOME} strategy ;
  KEYSET selectZoneIndex ;
  KEYSET runs ;
  Stack snipnet ;
  DICT *selectDict ; 
  DICT *ventilationDict ;
  DICT *chromDict ;
  DICT *targetDict ;
  DICT *snpDict ;
  DICT *snpNamDict ;
  const char *selected8kbFileName ;
  DICT *selected8kbDict ;
  BitSet selected8kbBitSet ;
  DICT *geneDict ;
  Stack aliExtendStack ;
  DICT *runDict ;
  DICT *target_classDict ;
  DICT *eeDict ;
  DICT *oldSnps ;
  AC_DB db ;
  BitSet run_sample ;
  ACEIN ai ; ACEOUT ao ; 
} SNP ;

static int snpPrettyNames (SNP *snp, AC_OBJ snpObj
			    , vTXT gTxt, vTXT cTxt, vTXT pTxt
			    , int *gMapp, int *gPosp
			    , char *gType, char *cType, vTXT pType
			    , vTXT location
			    , vTXT geneboxes, vTXT avgeneboxes, vTXT gsnippet, vTXT snippet, vTXT pSnippet
			   ) ;

typedef struct brsStruct { int pos, type, bp, bm, rp, rm, sp, sm, tp, tm ; short ee ; } BRS ;
typedef struct brs2Struct { int bp, bm, rp, rm, sp, sm, tp, tm ; int oldSnp ; char base ; } BRS2 ;

/*************************************************************************************/
/* parse the runQuality file
 * Two columns: position in read (sequencing machine cycle)
 * average fastqQuality at this position: -10 log10 (error proba)
 * all counts are mutiplied by the coefficient proposed here
 */

static int snpParseRunQualityFile (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai ;
  int part, x, nn[10], xMax = 0, nnn ;
  float q ;
  Array aa = 0 ;

  memset (nn, 0, sizeof(nn)) ;
  aa = snp->runQuality = arrayHandleCreate (200, int, snp->h) ;
  ai = aceInCreate (snp->runQualityFileName, 0, h) ;
  while (aceInCard (ai))
    {
      if (aceInInt (ai, &part) && aceInStep(ai, '\t') && aceInInt (ai, &x) && aceInStep(ai, '\t') && aceInFloat (ai, &q))
	{
	  if (q < 10) q = 10 ;
	  if (q > 30) q = 30 ;
	  q = 100  - 5 * (30 - q) ;
	  array (aa, part+10*x, int) = q ;
	  if (x > xMax) xMax = x ;
	  nn[part]++ ;
	}
      if (part > 9) 
	messcrash ("Sorry snpParseRunQualityFile can only handle 9 parts for a run, the first column says %d line %d in file %s"
		   , part, aceInStreamLine (ai), snp->runQualityFileName ) ;
    }

  x = arrayMax (aa) ;
  if (! x) 
    arrayDestroy (snp->runQuality)  ;
  else if (!nn[0]) /* create in column 0 the average of all the columns */
    {
      for (nnn = part = 0 ; part < 10 ; part++)
	if (nn[part])
	  nnn++ ;
      if (! nnn) nnn = 1 ;
       for (x = 0 ; x <= xMax ; x++)
	{
	  for (q = part = 0 ; part < 10 ; part++)
	    if (nn[part])
	      q +=  array (aa, part+10*x, int) ; 
	  array (aa, 10*x, int) = q /nnn ; 
	}
    }
  ac_free (h) ;
  return x ;
} /* snpParseRunQualityFile */

/*************************************************************************************/
/* parse the coordinates of genes mapping on the mito or genome target */
/* typically using WIGGLEREMAP directory */
typedef struct gtStruct { int gene, a1, a2 ; BOOL isUp ; } GT ;
typedef struct zoneStruct { int zone, a1, a2, vent ; } ZONE ;

static int snpParseTargetGeneFile (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *ccp = snp->targetGeneFileName + strlen (snp->targetGeneFileName) - 3 ;
  ACEIN ai ;
  int nn = 0, line = 0 ;
  Array aa = 0 ;
  GT *gt ;

  snp->targetGene = aa = arrayHandleCreate (1000, GT, snp->h) ;
  snp->geneDict = dictHandleCreate (100, snp->h) ;

  ai = aceInCreate (snp->targetGeneFileName, 0, h) ;
  aceInSpecial (ai, "\t\n") ;

  nn = 0 ;
  while (aceInCard (ai)) /* will close pp->targetFile */
    {
      line++ ;
      ccp = aceInWord (ai) ;
      if (ccp && *ccp == '#') continue ;
      gt = arrayp (aa, nn++, GT) ;

      dictAdd (snp->geneDict, ccp, &(gt->gene)) ;
      aceInInt (ai, &(gt->a1)) ;
      aceInInt (ai, &(gt->a2)) ;
      if (gt->a1 > gt->a2)
	{ int a = gt->a1 ; gt->a1 = gt->a2 ; gt->a2 = a ; gt->isUp = TRUE ; }
      nn++ ; /* success */
    }
  arrayMax (aa) = nn ;

  ac_free (h) ;

  return nn ;
} /* snpParseTargetGeneFile */

/*************************************************************************************/

static int snpParseFastaFile (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  char *cp, *cq ;
  ACEIN ai ;
  int state = 0, nn = 0, line = 0, mrna = 0 ;
  long  nbp = 0 ;
  Stack s ;

  snp->targetDna = s = stackHandleCreate (1000000, snp->h) ;
  snp->target2dna = arrayHandleCreate (1024, int, snp->h) ;
  snp->targetDnaLength = arrayHandleCreate (1024, int, snp->h) ;

  pushText (s, "toto") ; /* avoid zero */
  ai = aceInCreate (snp->fastaFileName, FALSE, h) ;
  aceInSpecial (ai, "\n") ;

  state = nn = 0 ;
  while (aceInCard (ai)) /* will close pp->targetFile */
    {
      line++ ;
      cp = aceInWord (ai) ;
      if (cp && *cp == '#') continue ;
      if (cp) switch (state)
	{
	case 1:
	case 2:
	  if (*cp != '>')
	    {
	      nbp += strlen (cp) ;
	      if (state == 1)
		{
		  state = 2 ;
		  pushText (s, cp) ;
		}
	      else
		catText (s, cp) ;
	      array (snp->targetDnaLength, mrna, int) = nbp ;
	      break ;
	    }
      /* fall thru to new sequence */
	case 0: /* expecting   >target_name */ 
	  if (*cp == '>')
	    {
	      if (*cp != '>' || ! *(cp+1))
		{
		  fprintf (stderr, "// bad character in file %s line %d, expecting a '>', got %s\n"
			   , snp->fastaFileName, line, cp) ;
		  return FALSE ;
		}
	      cq = strstr (cp, "|Gene|") ; if (cq) *cq = 0 ;
	      cq = strstr (cp, "|GENE|") ; if (cq) *cq = 0 ;
	      if (! strncmp (cp, ">MRNA:",6)) cp += 6 ;
	      if (! strncmp (cp, ">",1)) cp += 1 ;
	      /* in case -select, only collect the relevant sequences in the fasta file */
	      if ( ! snp->selectDict || dictFind (snp->selectDict, cp, 0))
		{
		  dictAdd (snp->targetDict, cp, &mrna) ;
		  state = 1 ; nn++ ;
		  array (snp->target2dna, mrna, int) = stackMark (s) ;
		}
	      else
		state = 0 ;
	    }
	  break ;
	}
    }
  ac_free (h) ;
  fprintf (stderr, "// %s snpParseFastaFile selected %d sequences %ldbp\n", timeShowNow(), nn, nbp) ; 
  return nn ;
} /* snpParseFastaFile */

/*************************************************************************************/
/*************************************************************************************/
#define NN 32
typedef struct ssStruct { int pos, count ; short snp, qual, strand ; } SS ;

static int ssOrder (const void *va, const void *vb)
{
  const SS *a = (const SS *)va, *b = (const SS *)vb ;
  int n ;

  n = a->pos - b->pos ; if (n) return n ;
  n = a->strand - b->strand ; if (n) return n ;
  n = a->snp - b->snp ; if (n) return n ;
  n = a->qual - b->qual ; if (n) return n ;

  return 0 ;
} /* ssOrder */

/*************************************************************************************/
#define NN 32
typedef struct ssmStruct { int run, variant, target, pos, leftShift 
			     , snp, count, coverR, coverS, coverT
			     ;} SSM ;

/*************************************************************************************/

typedef struct remapStruct { int target_class, gene, target, x1, x2, chrom, a1, a2 ; } RMP ;

/*****************/

static int atlasOrder1 (const void *va, const void *vb)
{
  const RMP *a = (const RMP *)va, *b = (const RMP *)vb ;
  int n ;
  n = a->target_class - b->target_class ; if (n) return n ;
  n = a->target - b->target ; if (n) return n ;
  n = a->x1 - b->x1 ; if (n) return n ;
  n = a->x2 - b->x2 ; if (n) return n ;
  n = a->chrom - b->chrom ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;

  return 0 ;
} /* atlasOrder1 */

/*****************/

static int atlasOrder2 (const void *va, const void *vb)
{
  const RMP *a = (const RMP *)va, *b = (const RMP *)vb ;
  int n ;
  n = a->chrom - b->chrom ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  n = a->target_class - b->target_class ; if (n) return n ;
  n = a->target - b->target ; if (n) return n ;
  n = a->x1 - b->x1 ; if (n) return n ;
  n = a->x2 - b->x2 ; if (n) return n ;

  return 0 ;
} /* atlasOrder1 */

/*************************************************************************************/
/* a tampon in bp, to allocate the snp to the promotor region of the gene */
#define PROMOTOR_DELTA 500

static int snpCreateAtlas (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai ;
  Array atlas = 0, geneAtlas = 0, map, geneMap ;
  const char *ccp ;
  RMP *up ;
  int nn = 0, x1, x2, a1, a2, target, chrom, gene, target_class ;

  x1 = x2 = a1 = a2 = 0 ;
  snp->target2geneAtlas = arrayHandleCreate (4, Array, snp->h) ;
  snp->target2exonAtlas = arrayHandleCreate (4, Array, snp->h) ;
  if (snp->remapFileName1)
    ai = aceInCreate (snp->remapFileName1, FALSE, h) ;
  else
    ai = aceInCreate (snp->remapFileName2, FALSE, h) ;
  aceInSpecial (ai, "\t\n") ;

  while (ai && aceInCard (ai))
    {
      ccp = aceInWord (ai) ; /* target_class */
      if (! ccp || *ccp == '#')
	continue ;
      if (snp->target_class && strcmp (ccp, (snp->target_class)))
	continue ;
      ccp = "any" ; /* ignore target class for the moment */
      dictAdd (snp->target_classDict, ccp, &target_class) ;
      atlas = array (snp->target2exonAtlas, target_class, Array) ;
      geneAtlas = array (snp->target2geneAtlas, target_class, Array) ;
      if (! atlas)
	{
	  atlas = arrayHandleCreate (100000, Array, snp->h) ;
	  array (snp->target2exonAtlas, target_class, Array) = atlas ;
	}
      if (! geneAtlas)
	{
	  geneAtlas = arrayHandleCreate (10000, Array, snp->h) ;
	  array (snp->target2geneAtlas, target_class, Array) = geneAtlas ;
	}
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* target */
      if (! ccp || *ccp == '#')
	continue ; 
      dictAdd (snp->targetDict, ccp, &target) ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &x1))  /* mRNA exon coordinates */
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &x2))  /* mRNA exon coordinates */
	continue ;
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* chrom */
      if (! ccp || *ccp == '#')
	continue ; 
      dictAdd (snp->chromDict, ccp, &chrom) ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &a1))  /* mRNA exon coordinates */
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &a2))  /* mRNA exon coordinates */
	continue ;
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; /* chrom */
      if (! ccp || *ccp == '#')
	continue ;
      dictAdd (snp->targetDict, ccp, &gene) ;
      if (1)
	{
	  map = array (atlas,  snp->remapFileName1 ? target : chrom, Array) ;
	  if (! map)
	    map = array (atlas, snp->remapFileName1 ? target : chrom, Array) = arrayHandleCreate (8, RMP, snp->h) ;
	  geneMap = array (geneAtlas,  snp->remapFileName1 ? target : chrom, Array) ;
	  if (! geneMap)
	    geneMap = array (geneAtlas, snp->remapFileName1 ? gene : chrom, Array) = arrayHandleCreate (8, RMP, snp->h) ;
	  nn++ ;
	  up = arrayp (map, arrayMax (map), RMP) ;
	  up->target = target ;
	  up->chrom = chrom ;
	  if (snp->remapFileName1 || a1 < a2)
	    {
	      up->x1 = x1 ; up->x2 = x2 ;
	      up->a1 = a1 ; up->a2 = a2 ;
	    }
	  else
	    {
	      up->x1 = x2 ; up->x2 = x1 ;
	      up->a1 = a2 ; up->a2 = a1 ;
	    }

	  up =  arrayMax (geneMap) ? arrayp (geneMap, arrayMax (geneMap) - 1, RMP) : 0 ;
	  if (! up || up->target != gene)
	    up = arrayp (geneMap, arrayMax (geneMap), RMP) ;
	  up->target = gene ;
	  up->chrom = chrom ;
	  if (a1 < a2)
	    {
	      a1 -= PROMOTOR_DELTA ;
	      if (! up->a1 || a1 < up->a1) up->a1 = a1 ; 
	      if (a2 > up->a2) up->a2 = a2 ;
	      up->x1 = 1 ; up->x2 = up->a2 - up->a1 + 1 ;
	    }
	  else
	    {
	      a1 += PROMOTOR_DELTA ;
	      if (! up->a1 || a2 < up->a1) up->a1 = a2 ; 
	      if (a1 > up->a2) up->a2 = a1 ;
	      up->x2 = 1 ; up->x1 = up->a2 - up->a1 + 1 ;
	    }
	}
    }

  for (a1 = 0 ; atlas && a1 < arrayMax (atlas) ; a1++)
    {
      map = array (atlas, a1, Array) ;
      if (map)
	arraySort (map, snp->remapFileName1 ? atlasOrder1 : atlasOrder2) ;
    }
  for (a1 = 0 ; geneAtlas && a1 < arrayMax (geneAtlas) ; a1++)
    {
      map = array (geneAtlas, a1, Array) ;
      if (map)
	{
	  arraySort (map, snp->remapFileName1 ? atlasOrder1 : atlasOrder2) ;
	  arrayCompress (map) ;
	}
    }
   
  fprintf (stderr, "snpCreateAtlas found %d exons in file %s\n"
	   , nn
	   , snp->remapFileName1 ? snp->remapFileName1 :  snp->remapFileName2
	   ) ;

  ac_free (h) ;
  return nn ;
} /* snpCreateAtlas */

/*************************************************************************************/
/* 
 * We import all the relevant data from the ZZ database
 * Remap the transcript variants into genome coordinates
 */

static BOOL snpRemap1Do (SNP *snp, int mrna, int x1, int *chromp, int *a1p, int *strandp)
{
  int ii, target_class = 0 ;
  RMP *up ;
  Array atlas, map ;
  
  dictAdd (snp->target_classDict, "any", &target_class) ;
  if (snp->target2exonAtlas && 
      target_class < arrayMax (snp->target2exonAtlas) && 
      (atlas =  arr (snp->target2exonAtlas, target_class, Array)) &&
      mrna < arrayMax(atlas) &&
      (map = array (atlas, mrna, Array))
      )
    {
      for (ii = 0, up = arrp (map, 0, RMP) ; ii < arrayMax (map) ; ii++, up++)
	{
	  if (up->x1 <= x1 && up->x2 >= x1)
	    {
	      if (up->a1 < up->a2)
		{ *a1p = up->a1 + x1 - up->x1 ; *strandp = 1 ; *chromp = up->chrom ; return TRUE ; }
	      else
		{ *a1p = up->a1 - x1 + up->x1 ; *strandp = - 1 ; *chromp = up->chrom ; return TRUE ; }
	    }
	  else if (x1 < up->x1)
	    break ;
	}
    }
  return FALSE ;
} /* snpRemap1Do */

/*************************************************************************************/
/* scan the VariantDB acedb database
 * add the remap info 
 * Remap the transcript variants into genome coordinates
 */
static int snpRemap1 (SNP *snp)
{
  AC_HANDLE  h1 = 0, h = ac_new_handle () ;
  AC_ITER iter ;
  AC_OBJ variant = 0 ;
  AC_TABLE mrnaTable = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  int pos, a1, strand, chrom, mrna, nn1 = 0, nn2 = 0 ;
  const char *mm ;
  const char *errors = 0 ;
  
  if (snp->db)
    {
      iter = ac_query_iter (snp->db, TRUE, "find variant mRNA && ! IntMap", 0, h) ;
      while (ac_free (variant), ac_free (h1), variant = ac_iter_obj (iter))
	{
	  h1 = ac_new_handle () ;
	  nn1++ ;
	  mrnaTable = ac_tag_table (variant, "mRNA", h1) ;
	  if (mrnaTable)
	    {
	      pos = ac_table_int (mrnaTable, 0, 1, 0) ;
	      mm = ac_table_printable (mrnaTable, 0, 0, 0) ;
	      if (mm && dictFind (snp->targetDict, mm, &mrna) && snpRemap1Do (snp, mrna, pos, &chrom, &a1, &strand))
		{
		  nn2++ ;
		  vtxtPrintf (txt, "Variant %s\n", ac_protect (ac_name (variant), h1)) ;
		  vtxtPrintf (txt, "IntMap %s %d %d\n\n", dictName (snp->chromDict, chrom), a1, strand) ; 
		}
	    }
	}
      ac_parse (snp->db, vtxtPtr (txt), &errors, 0, h) ; 
    }
  fprintf(stderr, "snpRemap1 found %d variants remapped %d\n", nn1, nn2) ;

  if (errors && *errors) fprintf(stderr, "snpRemap parsing error %s\n", errors) ;
  ac_free (h1) ;
  ac_free (h) ;
  return nn2 ;
} /* snpRemap1 */

/*************************************************************************************/
/* 
 * We import all the relevant data from the ZZ database
 *   Remap the genome variants into transcript coordinates
 *   Remap to the genebox without giving precise coordinates
 */

static BOOL snpRemap2geneBoxDo (SNP *snp, int chrom, int pos, int *JJp, int *geneBoxp, int *x1p, int *strandp)
{
  int ii, target_class = 0 ;
  RMP *up ;
  Array atlas, map ;
  
  dictAdd (snp->target_classDict, "any", &target_class) ;
  if (snp->target2geneAtlas && 
      target_class < arrayMax (snp->target2geneAtlas) && 
      (atlas =  arr (snp->target2geneAtlas, target_class, Array)) &&
      chrom < arrayMax(atlas) &&
      (map = array (atlas, chrom, Array))
      )
    {
      for (ii = *JJp, up = arrp (map, ii, RMP) ; ii < arrayMax (map) ; ii++, up++)
	{
	  if (up->a1 <= pos && up->a2 >= pos)
	    {
	      if (up->x1 < up->x2)
		{ *x1p = up->x1 + pos - up->a1 - PROMOTOR_DELTA ; *strandp = 1 ; *geneBoxp = up->target ; *JJp = ii + 1 ; return TRUE ; }
	      else
		{ *x1p = up->x1 - pos + up->a1  - PROMOTOR_DELTA ; *strandp = - 1 ; *geneBoxp = up->target ; *JJp = ii + 1 ;  return TRUE ; }
	    }
	  else if (pos < up->a1)
	    break ;
	}
    }
  return FALSE ;
} /* snpRemap2Do */

static BOOL snpRemap2Do (SNP *snp, int chrom, int pos, int *JJp, int *mrnap, int *x1p, int *strandp)
{
  int ii, target_class = 0 ;
  RMP *up ;
  Array atlas, map ;
  
  dictAdd (snp->target_classDict, "any", &target_class) ;
  if (snp->target2exonAtlas && 
      target_class < arrayMax (snp->target2exonAtlas) && 
      (atlas =  arr (snp->target2exonAtlas, target_class, Array)) &&
      chrom < arrayMax(atlas) &&
      (map = array (atlas, chrom, Array))
      )
    {
      for (ii = *JJp, up = arrp (map, ii, RMP) ; ii < arrayMax (map) ; ii++, up++)
	{
	  if (up->a1 <= pos && up->a2 >= pos)
	    {
	      if (up->x1 < up->x2)
		{ *x1p = up->x1 + pos - up->a1 ; *strandp = 1 ; *mrnap = up->target ; *JJp = ii + 1 ; return TRUE ; }
	      else
		{ *x1p = up->x1 - pos + up->a1 ; *strandp = - 1 ; *mrnap = up->target ; *JJp = ii + 1 ;  return TRUE ; }
	    }
	  else if (pos < up->a1)
	    break ;
	}
    }
  return FALSE ;
} /* snpRemap2Do */

/*************************************************************************************/
/* scan the VariantDB acedb database
 * add the remap info 
 *   Remap the genome variants into transcript coordinates\n"
 */
static int snpRemap2 (SNP *snp)
{
  AC_HANDLE  h1 = 0, h = ac_new_handle () ;
  AC_ITER iter ;
  AC_OBJ variant = 0 ;
  AC_TABLE intMapTable = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  int pos, pos2, x1 = 0, strand, chrom, mrna = 0, geneBox = 0, nn1 = 0, nn2 = 0, JJ=0 ;
  const char *chromNam ;
  const char *errors = 0 ;
  
  iter = ac_query_iter (snp->db, TRUE, "find variant IntMap && ! geneBox", 0, h) ;
  while (ac_free (variant), ac_free (h1), variant = ac_iter_obj (iter))
    {
      h1 = ac_new_handle () ;
      nn1++ ;
      intMapTable = ac_tag_table (variant, "IntMap", h1) ;
      if (intMapTable)
	{
	  chromNam = ac_table_printable (intMapTable, 0, 0, 0) ;
	  pos = ac_table_int (intMapTable, 0, 1, 0) ;
	  pos2 = ac_table_int (intMapTable, 0, 2, 0) ;
	  if (pos2 == 1 || pos2 > pos) pos2 = 1 ;
	  else if (pos2 == -1 || pos2 < pos) pos2 = 1 ;
	  
	  if (pos && dictFind (snp->chromDict, chromNam, &chrom))
	    {
	      JJ = 0 ;
	      if (! ac_has_tag (variant, "mRNA"))
		while (snpRemap2Do (snp, chrom, pos, &JJ, &mrna, &x1, &strand))
		  {
		    nn2++ ;
		    vtxtPrintf (txt, "Variant %s\n", ac_protect (ac_name (variant), h1)) ;
		    vtxtPrintf (txt, "mRNA %s %d %d\n\n", dictName (snp->targetDict, mrna), x1, strand * pos2) ; 
		  }
	      JJ = 0 ;
	      while (snpRemap2geneBoxDo (snp, chrom, pos, &JJ, &geneBox, &x1, &strand))
		{
		  nn2++ ;
		  vtxtPrintf (txt, "Variant %s\n", ac_protect (ac_name (variant), h1)) ;
		  vtxtPrintf (txt, "GeneBox %s %d %d\n\n", dictName (snp->targetDict, geneBox), x1, strand * pos2) ; 
		}
	    }
	}
    }
  fprintf(stderr, "snpRemap2 found %d variants remapped %d\n", nn1, nn2) ;
  if (vtxtPtr (txt))
    {
      ACEOUT ao = aceOutCreate ("toto", 0, 0, h) ;
      aceOut (ao,  vtxtPtr (txt)) ;
      ac_free (ao) ;
    }
  ac_parse (snp->db, vtxtPtr (txt), &errors, 0, h) ; 
  if (*errors) fprintf(stderr, "snpRemap parsing error %s\n", errors) ;
  ac_free (h1) ;
  ac_free (h) ;
  return nn2 ;
} /* snpRemap2 */

/*************************************************************************************/

static void showSsm (Array aa)
{
  int ii ;
  SSM *ssm ;
  
  for (ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      ssm = arrp (aa, ii, SSM) ; 
      fprintf (stderr, "%d:\t%d\tsnp:%d\t%d\tR:%d\tS:%d\tT:%d\trun=%d\n"
	       , ii
	       , ssm->pos
	       , ssm->snp
	       , ssm->count
	       , ssm->coverR
	       , ssm->coverS
	       , ssm->coverT
	       , ssm->run
	       ) ;
    }
  
  if (0) showSsm (0) ;
} /* showSsm */
  
/*************************************************************************************/
#ifdef JUNK
static int ssmOrderBySnp (const void *va, const void *vb)
{
  const SSM *a = (const SSM *)va, *b = (const SSM *)vb ;
  int n ;

  n = a->target - b->target ; if (n) return n ;
  n = a->pos - b->pos ; if (n) return n ;
  n = a->snp - b->snp ; if (n) return n ;
  n = a->run - b->run ; if (n) return n ;

  return 0 ;
} /* ssmOrderBySnp */
#endif
/*************************************************************************************/

static int ssmDeepTableReportOrder (const void *va, const void *vb)
{
  const SSM *a = (const SSM *)va, *b = (const SSM *)vb ;
  int n ;

  n = a->target - b->target ; if (n) return n ;
  n = a->pos - b->pos ; if (n) return n ;
  n = a->variant - b->variant ; if (n) return n ;
  n = a->run - b->run ; if (n) return n ;

  return 0 ;
} /* ssmDeepTableReportOrder */

/*************************************************************************************/

static void snpExportErrors (SNP* snp, int nAli, long int nBpAli, Array errorsP, DICT *eeDict)
{
  int i ;
  SS *ee ;

  /* export the total counts */

  aceOutf (snp->ao, "nAli\t%d\n", nAli) ;
  aceOutf (snp->ao, "nBpAli\t%ld\n", nBpAli) ;

  if (errorsP) /* export the error profile */
    {
      for (i = 0 ; i < arrayMax (errorsP) ; i++)
	{
	  ee = arrp (errorsP, i, SS) ;
	  if (ee->count)
	    aceOutf(snp->ao, "ERR\t%s\t%d\t%d\n", dictName(eeDict, ee->snp), ee->qual, ee->count/100) ;
	}
    }
} /* snpExportErrors */

/*************************************************************************************/

static DICT *eeDictCreate (SNP *snp, AC_HANDLE h)
{
  DICT *eeDict = dictHandleCreate (10000, h) ;

  /* force wild type single letters fisrt */
  dictAdd (eeDict, "toto", 0) ;
  dictAdd (eeDict, "a", 0) ;
  dictAdd (eeDict, "t", 0) ;
  dictAdd (eeDict, "g", 0) ;
  dictAdd (eeDict, "c", 0) ;
  dictAdd (eeDict, "n", 0) ;
  dictAdd (eeDict, "W", 0) ;
  
  dictAdd (eeDict, "a>t", 0) ;
  dictAdd (eeDict, "a>g", 0) ;
  dictAdd (eeDict, "a>c", 0) ;
  dictAdd (eeDict, "a>n", 0) ;
  
  dictAdd (eeDict, "t>a", 0) ;
  dictAdd (eeDict, "t>g", 0) ;
  dictAdd (eeDict, "t>c", 0) ;
  dictAdd (eeDict, "t>n", 0) ;
  
  dictAdd (eeDict, "g>a", 0) ;
  dictAdd (eeDict, "g>t", 0) ;
  dictAdd (eeDict, "g>c", 0) ;
  dictAdd (eeDict, "g>n", 0) ;
  
  dictAdd (eeDict, "c>a", 0) ;
  dictAdd (eeDict, "c>t", 0) ;
  dictAdd (eeDict, "c>g", 0) ;
  dictAdd (eeDict, "c>n", &(snp->eeLastSub)) ;
  
  dictAdd (eeDict, "+a", 0) ;
  dictAdd (eeDict, "+t", 0) ;
  dictAdd (eeDict, "+g", 0) ;
  dictAdd (eeDict, "+c", 0) ;
  
  dictAdd (eeDict, "*+a", 0) ;
  dictAdd (eeDict, "*+t", 0) ;
  dictAdd (eeDict, "*+g", 0) ;
  dictAdd (eeDict, "*+c",  &(snp->eeLastInsert)) ;
  
  dictAdd (eeDict, "-a", 0) ;
  dictAdd (eeDict, "-t", 0) ;
  dictAdd (eeDict, "-g", 0) ;
  dictAdd (eeDict, "-c", 0) ;
  
  dictAdd (eeDict, "*-a", 0) ;
  dictAdd (eeDict, "*-t", 0) ;
  dictAdd (eeDict, "*-g", 0) ;
  dictAdd (eeDict, "*-c", 0) ;

  return eeDict ;
} /* eeDictCreate */

/*************************************************************************************/
/* 19028711:cc>oo,19028712:cg>oo  should become 19028711:oo2,19028712:oo3,19028713:oo1 */
static int snpCountSnpsOO (char *snpBuf, char *xsnpBuf, int s1, int s2)
{
  char buf[1000] ;
  char xbuf[1000] ;
  char *cp, *cq, *cr, *cnew = snpBuf ; ;
  char *cpx, *cqx, *crx, *cnewx = xsnpBuf ; ;
  int i, j, k, nn = 0, nnx = 0 ;
  int pos[100], noo = 0 ;
  BOOL hasMore ;

  strncpy (buf+1,snpBuf, 999) ; buf[0] = ',' ;
  strncpy (xbuf+1,xsnpBuf, 999) ; xbuf[0] = ',' ;
  /* split on , */

  *cnew = *cnewx = 0 ;
  cq = buf ;
  hasMore = TRUE ;
  while (hasMore)
    {
      cp = cq + 1 ;
      cq = strchr (cp+1, ',') ; /* next block or zero */
      if (cq) 
	{ hasMore = TRUE ; *cq = 0 ; }
      else
	hasMore = FALSE ;
      cr = strstr (cp, "oo") ;
      if (! cr)  /* just copy the buffer */
	{
	  if (nn++) *(cnew-1) = ',' ;
	  while ((*cnew++ = *cp++)) ;
	}
      else
	{
	  cr = strchr(cp, ':') ;
	  *cr = 0 ;
	  pos[noo++] = i = atoi(cp) ;
	}
    }
  cqx = xbuf ;
  hasMore = TRUE ;
  while (hasMore)
    {
      cpx = cqx + 1 ;
      cqx = strchr (cpx+1, ',') ; /* next block or zero */
      if (cqx)
	{ hasMore = TRUE ; *cqx = 0 ; }
      else
	hasMore = FALSE ;

      crx = strstr (cpx, "oo") ;
      if (! crx)  /* just copy the buffer */
	{
	  if (nnx++) *(cnewx-1) = ',' ;
	  while ((*cnewx++ = *cpx++)) ;
	}
      else
	{
	  crx = strchr(cpx, ':') ;
	  *crx = 0 ;
	}
    }
  if (noo) /* export the oo positions */
    {
      for (i = 0 ; i < noo ; i++)
	if (pos[i]) /* export the positions only once */
	  {
	    k = 6 ;
	    for (j = i+1 ; j < noo ; j++)
	      {
		if (pos[i] == pos[j]) pos[j] = 0 ;
		if (pos[i] == pos[j] - 1) k = 2 ;
	      }
	    for (j = i-1 ; j > 0 ; j--)
	      {
		if (pos[i] == pos[j] + 1) k |= 0x1 ;
	      }
	
	    /* count only within the acceptable zone */
	    if (pos[i] >= s1 && pos[i] <= s2) 
	      {
		cr = messprintf ("%d:o>o%d", pos[i],k & 0x3) ;
		if (nn++) *(cnew-1) = ',' ;
		while ((*cnew++ = *cr++)) ;
	      }
	    if (k & 0x4 && pos[i] + 1 >= s1 && pos[i] + 1 <= s2) 
	      {
		cr = messprintf ("%d:o>o%d", pos[i] + 1,1) ;
		if (nn++) *(cnew-1) = ',' ;
		while ((*cnew++ = *cr++)) ;
	      }
	    if (cnew > snpBuf + 900) break ;
	  }
    }
  cp = snpBuf + strlen(snpBuf) ; *++cp = 0 ;
  return noo ;
} /* snpCountSnpsOO */

/*************************************************************************************/
/* BRST cover
 * R repeat or roche, confirmed letter count, i.e. we cross over to a different letter on both side
 * S subtitution or solexa, we cover the letter
 * T trasition or lifeTech incoherent with the reference
 * In each case we require that the cover is 8bp on each side
 * so a would be mutation would have been detected and not clipped
 * Otherwise, we would artificially increase the cover of the wild type
 * For example with 36-mers detecting an heterozygote mutation
 * 16 positions are lost so the rate apparent rate would not be 1/2
 * but 20/20+36 = .35 
 */

static int runQ (Array runQuality, int runQualityPrefix, int x)
{
  int q = 100 ;

  if (x > 0 && runQuality)
    {
      int  iq = runQualityPrefix + 10 * x ;

      if (runQualityPrefix >= 0  && iq < arrayMax (runQuality))
	q = arr (runQuality, iq, int) ;
      else
	q = 0 ;

    }
  return q ;
} /* runQ */

static int snpHits2BRS (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  int leftOverhang, rightOverhang ; 
  int leftROverhang, rightROverhang ; 
  int leftOverhangEffectif, rightOverhangEffectif ;
  int ne = 0, nme = 0, mult, unicity ;
  int nn = 0, i, *ip, ali, toBeAligned ;
  int a1, a2, r1, r2, s1, s2, dummy, x, x1, x2, q ;
  int pos = 0, sn, target, targetZone = 0, nAli = 0, qual, count, nN ;
  int oldmrna, oldpos, oldsn, strand = 0, iSelect = 0, zzSelect = 0, zzda ;
  long int nBpAli = 0 ;  
  BOOL     isDown = TRUE ;
  Array errorsP = arrayHandleCreate (1000, SS, h) ;
  SS *ss, *ee ;
  ACEIN ai = snp->ai ;
  ACEOUT ao = snp->ao ;
  Array aosU = arrayHandleCreate (64,ACEOUT,h) ;
  Array aosNU = arrayHandleCreate (64,ACEOUT,h) ;
  const char *ccp ;
  char *cp, *cq, *cr, *cpx, *cqx, *crx, snpBuf[1000], xsnpBuf[1000], targetBuf[1000], cc ;
  DICT *eeDict = eeDictCreate (snp, h) ;
  DICT *posDict = 0 ;
  DICT *zzDict = dictHandleCreate (10000, h) ;
  Array posDictArray = arrayHandleCreate (1000, DICT*,h) ;
  char posBuffer[4096] ;
  int posIndex, target_class ;
  Array cover, covers, snps ;
  Array coversDownArray = arrayHandleCreate (1000, Array, h) ;
  Array coversUpArray = arrayHandleCreate (1000, Array, h) ;
  char tagBuf[4095] ;  char newTagBuf[4095] ;
  int nmembl, memsize, tc ;
  int oldXSnp = 0, oldXRun = 0, oldXTarget = 0 ; /* old exported value */
  ZONE *zp = 0, *zzp = 0 ;
  BOOL hasMore = FALSE ;

  BOOL hasPair = FALSE ;  /* we autodetect the pairs */
  BOOL goodPair = TRUE ;
  
  KEYSET zzDas = keySetHandleCreate (h) ;
  KEYSET targetZone2target = keySetHandleCreate (h) ;
  Array snpsArray ;
  Array snpsUpArray = arrayHandleCreate (64, Array, h) ;
  Array snpsDownArray = arrayHandleCreate (64, Array, h) ;
  int snpTotal = 0 ;
  int ventU = 0, ventNU = 0, ventUF = 0, ventNUF = 0 ;
  int runQualityPrefix = 0 ;
  int  Remove_inserts_shorter_than = snp->Remove_inserts_shorter_than ;
  Array runQuality = snp->runQuality ;
  BOOL dropMultiplicity = snp->dropMultiplicity ;

  runQuality = 0 ;
  if (snp->ventilate) aceInSpecial (ai, "\n") ;
  memset (tagBuf, 0, sizeof(tagBuf)) ;

  int ooE1, ooE2 = 0, ooE3 = 0 ;
  dictAdd (eeDict, "o>o1", &ooE1) ;
  dictAdd (eeDict, "o>o2", &ooE2) ;
  dictAdd (eeDict, "o>o3", &ooE3) ;

  leftROverhang = rightROverhang = 8 ;
  if (snp->strategy == STRATEGY_RNA_SEQ)
    {
      leftOverhang = rightOverhang = 8 ;
    }
  else
    {
      leftOverhang = rightOverhang = 4 ;
    }

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (!ccp || *ccp == '#')
	continue ;
      if (runQuality) 
	{
	  const char *ccq = ccp + strlen(ccp) - 1 ;
	  int i = 2 ;

	  runQualityPrefix = 1 ;
	  if (*ccq == '>') runQualityPrefix = 1 ;
	  else if (*ccq == '<') runQualityPrefix = 2 ;
	  else if (0) /* this is very dangerous if the users follow a different naming convention */
	    {
	      while (ccq > ccp && *ccq != '/' && i--) ccq-- ;
	      runQualityPrefix = (*ccq == '/') ? atoi (ccq+1) : 0 ;
	      if (runQualityPrefix < 0 || runQualityPrefix > 9)
		continue ; /* reject this probe */
	    }
	}
      if (! strncmp (tagBuf, ccp, 4095)) continue ;  /* count each tag only once */
      strncpy (newTagBuf, ccp, 4095) ;
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; /* score */
      aceInStep (ai, '\t') ; aceInInt (ai, &mult) ; /* multiplicity */

      count = dropMultiplicity ? 1 : mult ; /* 2014_12_12: new option, was mult since 2011_06_13, was 1 before */

      /* check partial : score is per pair but ali and toBeAligned are per read
       * so the constraint on >minAliPerCent applies independantly to each read
       *    2017_09_14, after FDA-MOD challenge, we relax and accept ali > 100 even is < snp->minAliPerCent
       */
      ali = toBeAligned = 0 ;
      aceInStep (ai, '\t') ; aceInInt (ai, &toBeAligned) ;        /* clipped length */
      aceInStep (ai, '\t') ; aceInInt (ai, &ali) ;  /* length aligned */
      if (ali < snp->minAli) continue ;
      if (ali < 140 &&                    /* 100 per read translates into 200 for the pair which seems secure */
	  100 * ali < snp->minAliPerCent * toBeAligned)
	continue ;

      aceInStep (ai, '\t') ; aceInInt (ai, &x1) ;   /* x1 tag coord */
      aceInStep (ai, '\t') ; aceInInt (ai, &x2) ;   /* x2 tag coord */

      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; /* target_class */
      if (!ccp)
	continue ;
      if (snp->mito && strcmp (ccp, "A_mito"))
	continue ;
      if (snp->target_class && ! strstr (ccp, snp->target_class ))
	continue ;
      dictAdd (snp->target_classDict, ccp, &target_class) ; tc = *ccp ;
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; /* gene ignore */
      aceInStep (ai, '\t') ; aceInInt (ai, &unicity) ;        /* unicity, ignore */
      if (snp->unique && unicity > 1)
	continue ;
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; /* mrna */
      if (! ccp)
	continue ;
      if (!strncmp (ccp, "MRNA:",5)) ccp+= 5 ;

      strncpy (targetBuf, ccp, 999) ; /* mrna */
      dictAdd (snp->targetDict, targetBuf, &target) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &a1) ;
      aceInStep (ai, '\t') ; aceInInt (ai, &a2) ;

      /* ATTENTION in ventilate case we will count the aceInWord till we reach the pairing info */

      isDown = TRUE ;
      if (a1 > a2)
	{ 
	  isDown = FALSE ;
	  i = a1 ; a1 = a2 ; a2 = i ; 
	  i = x1 ; x1 = x2 ; x2 = i ; 
	}

      /*
	2011_03_24
	It is quasi impossible to remap before constructing the SNPs because
	we cannot work per base, we have to work per segment

	if (snp->target2exonAtlas && ! snpRemapSegment (aa, snp, target, a1, a2)))
	continue ;
      */

      zzp = 0 ; zzda = zzSelect = 0 ; /* if known zone, declare the coordinates inside the zone */
      if (snp->selectDict && !(snp->ventilate && tc < 'C'))
	{
	  BOOL ok = FALSE ;
	  int zone, iSelectMax =  snp->selectZone ? arrayMax (snp->selectZone) - 1 : 0 ;

	  tc=0 ;
	  if (dictFind (snp->selectDict, targetBuf, &zone))
	    {
	      if (snp->selectZone)
		{
		  iSelect = keySet (snp->selectZoneIndex, zone) ;
		  zp = arrp (snp->selectZone, iSelect, ZONE) ;
		  while (iSelect > 0 && (zp->zone > zone || zp->a1 > a2)) { zp-- ; iSelect-- ; }
		  while (iSelect < iSelectMax  && (zp->zone < zone || (zp->a2 && zp->a2 < a1))) { zp++ ; iSelect++ ; }
		  if (zp->zone == zone && (zp->a1 + zp->a2 == 0 || (zp->a1 <= a2 && zp->a2 >= a1) || (zp->a1 <= a1 && zp->a2 >= a2)))
		    { 
		      ok = TRUE ; 
		      if (zp->a1 + zp->a2 > 0) 
			{ 
			  zzp = zp ;  zzSelect = iSelect ; zzda = zp->a1 ; 
			  keySet (snp->selectZoneIndex, zone) = iSelect ;
			}
		    }
		}
	      else
		ok = TRUE ;
	    }
	  if (!ok)
	    continue ;
	}
      strcpy (tagBuf, newTagBuf) ;

      /* do not compute the SNP, rather reexport the line in N-files according to the zone */
      if (snp->ventilate)
	{
	  char *cpn = 0 ;

	  if (! zzp && !tc) continue ;
	  if (unicity > 2) continue ;

	  if (1)
	    {
	      int i, deltaPair = 0 ;
	      
	      /* position on col 22 = pair fate */
	      for (i = 14 ; i<= 21 ; i++)
		{ aceInStep (ai, '\t') ; aceInWord (ai) ; }
	      aceInStep (ai, '\t') ; aceInInt (ai, &deltaPair) ;
	      
	      if (deltaPair > 30 || deltaPair < -30) /* at least one correctly annotated pair was found */
		hasPair = TRUE ;
	      goodPair = TRUE ;
	      if (hasPair && deltaPair < 0 && deltaPair > NON_COMPATIBLE_PAIR && deltaPair != -2 && deltaPair != -5) /* synchronize to hack with wiggle.c and bestali.c*/
		goodPair = FALSE ;
	      if ( Remove_inserts_shorter_than && deltaPair <  Remove_inserts_shorter_than && deltaPair > - Remove_inserts_shorter_than)
		goodPair = FALSE ;
#ifdef JUNK
                synchronize with snp.c
		int mateAliDirect  = 0, mateToAliDirect  = 0,  mateScoreDirect = 0 ; /* in case they are available in cols 23 24 */
		if (aceInStep (ai, '\t') &&  aceInInt (ai, &mateScore) && aceInStep (ai, '\t') &&  aceInInt (ai, &mateAli) && aceInStep (ai, '\t') &&  aceInInt (ai, &mateToAli) )
		  {
		    mateAli = mateAliDirect ; mateScore = mateScoreDirect ; mateToAli = mateToAliDirect ; /* if I find them directly that is ok */
		  }
#endif
	    }
	  if (hasPair && ! goodPair)
	    continue ;
	  

	  /* check unicity to ventilate into u or nu files */
	  if (unicity < 2)
	    {
	      ventU++ ;
	      ao = array (aosU, zzp ? zzp->vent : 10000 + (tc - '0' + 1), ACEOUT) ;
	      if (!ao)
		{
		  ventUF++ ;
		  if (zzp) cpn = messprintf(".%s.hits.u", dictName (snp->ventilationDict, zzp->vent)) ;
		  else if (tc == '0') cpn = ".SpikeIn.hits.u" ;
		  else if (tc == 'A') cpn = ".mito.hits.u" ;
		  else if (tc == 'B') cpn = ".rrna.hits.u" ;
		  ao = array (aosU,  zzp ? zzp->vent : 10000 + (tc - '0' + 1), ACEOUT) = aceOutCreate (snp->outFileName, cpn, snp->gzo, h) ;
		}
	    }
	  else
	    {
	      ventNU++ ;
	      ao = array (aosNU, zzp ? zzp->vent : 10000 + (tc - '0' + 1), ACEOUT) ;
	      if (!ao)
		{
		  ventNUF++ ;
		  if (zzp) cpn =  messprintf(".%s.hits.nu", dictName (snp->ventilationDict, zzp->vent)) ;
		  else if (tc == '0') cpn = ".SpikeIn.hits.nu" ;
		  else if (tc == 'A') cpn = ".mito.hits.nu" ;
		  else if (tc == 'B') cpn = ".rrna.hits.nu" ;
		  ao = array (aosNU,  zzp ? zzp->vent : 10000 + (tc - '0' + 1), ACEOUT) = aceOutCreate (snp->outFileName, cpn, snp->gzo, h) ;
		}
	    }
	    
	  aceInCardBack (ai) ;
	  aceInCard (ai) ;
	  cp = aceInPos (ai) ;
	  for (i = 1, cq = cp ; i <= 17 ; i++)
	    if (cq) cq = strchr (cq+1, '\t') ;
	  if (cq) 
	    {
	      *cq++ = '\n' ;
	      *cq = 0 ;
	      aceOut (ao, cp) ;
	    }
	  continue ;
	}

      dictAdd (zzDict, messprintf("%s_._%d", targetBuf, iSelect), &targetZone) ; 
      keySet (targetZone2target, targetZone) = target ;

      nN = 0 ; aceInStep (ai, '\t') ; aceInInt (ai, &nN) ;   
      aceInStep (ai, '\t') ; aceInInt (ai, &dummy) ;        /* nerror, ignore */

      snpsArray = isDown ? snpsDownArray : snpsUpArray ;
      snps = array (snpsArray, targetZone, Array) ;
      if (! snps) 
	array (snpsArray, targetZone, Array) = snps = arrayHandleCreate (1000, SS, h) ;

      covers = isDown ? coversDownArray : coversUpArray ;
      cover = array (covers, targetZone, Array) ;
      if (! cover)
	cover = array (covers, targetZone, Array) = arrayHandleCreate (2 * (a2 + 1 - zzda), int, h) ;

      keySet (zzDas, targetZone) = zzda ;

      /* S estimate */
      leftOverhangEffectif = leftOverhang ;
      rightOverhangEffectif = rightOverhang ;
      r1 = a1 ; r2 = a2 ;
     
      if (x1 < x2)
	{
	  if (x1 <=  leftOverhang)
	    leftOverhangEffectif = leftOverhang  - x1 + 1 ;
	  else
	    leftOverhangEffectif = 0 ;
	  a1 += leftOverhangEffectif ;
	  if (toBeAligned - x2 <  rightOverhang)
	    rightOverhangEffectif = rightOverhang - (toBeAligned - x2) ; 
	  else
	    rightOverhangEffectif =  0 ;
	  a2 -= rightOverhangEffectif ;
	}
      else
	{
	  if (x2 <=  leftOverhang)
	    leftOverhangEffectif = leftOverhang  - x2 + 1 ;
	  else
	    leftOverhangEffectif = 0 ;
	  a1 += leftOverhangEffectif ;
	  if (toBeAligned - x1 <  rightOverhang)
	    rightOverhangEffectif = rightOverhang - (toBeAligned - x1) ; 
	  else
	    rightOverhangEffectif = 0 ;
	  a2 -= rightOverhangEffectif ;
	}
 
      /* R case, we never detect Indel closer than 8, so we must reconsider the overhang */
     if (leftROverhang > leftOverhang || rightROverhang > rightOverhang)
       {
	 leftOverhangEffectif = leftROverhang ;
	 rightOverhangEffectif = rightROverhang ;
	 
	 if (x1 < x2)
	   {
	     if (x1 <=  leftROverhang)
	       leftOverhangEffectif = leftROverhang  - x1 + 1 ;
	     else
	       leftOverhangEffectif = 0 ;
	     r1 += leftOverhangEffectif ;
	     if (toBeAligned - x2 <  rightROverhang)
	       rightOverhangEffectif = rightROverhang - (toBeAligned - x2) ; 
	     else
	       rightOverhangEffectif =  0 ;
	     r2 -= rightOverhangEffectif ;
	   }
	 else
	   {
	     if (x2 <  leftROverhang)
	       leftOverhangEffectif = leftROverhang  - x2 + 1 ;
	     else
	       leftOverhangEffectif = 0 ;
	     r1 += leftOverhangEffectif ;
	     if (toBeAligned - x1 <  rightROverhang)
	       rightOverhangEffectif = rightROverhang - (toBeAligned - x1) ; 
	     else
	       rightOverhangEffectif =  0 ;
	     r2 -= rightOverhangEffectif ;
	   }
       }
     else
       {
	 r1 = a1 ; r2 = a2 ; /* inherit S values */
       }

      /* because the insert is presented on the next base and there is a also a deficit of 1 for the deletions */
      /* R cover assymetry, wait for a letter change before counting */
      if (snp->target2dna && snp->targetDna && target < arrayMax(snp->target2dna))
	{
	  i = arr (snp->target2dna, target, int) ;
	  cp = stackText (snp->targetDna, i + r1 - 1) ; cc = *cp ;
	  while (*cp++ == cc) { r1++ ; }

	  if (1) /* 2015_03_15  this makes the R count flat over the repeat, possibly better if the read is on negative strand ? */
	    {
	      i = arr (snp->target2dna, target, int) ;
	      cp = stackText (snp->targetDna, i + r2 - 1) ; cc = *cp ;
	      while (cp[0] == cp[1]) { cp-- ; r2-- ; }
	    }
	}

      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;      /* error in tag coords */
      strncpy (xsnpBuf, ccp ? ccp : " ", 999) ;  
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; /* error in target coords */
      if (ccp)  /* complex correction, we now hope that it is optimal */
	{
	  if ((strstr (ccp, ">oo"))) /* we may need to reclip 2 more bases ; */
	    {
	      int noo = 0, oo[256], x = 0, dx ;
	      const char *ccq = ccp, *ccr ;

	      /* find all the oo coordinates in read space */
	      for (ccq = ccp ; noo < 255 && (ccq = strstr (ccq, ">oo")) ; ccq++)
		{
		  for (ccr = ccq - 1 ; *ccr != ':' && ccr > ccp ; ccr--) ;
		  if (*ccr == ':')
		    {
		      for (--ccr ; *ccr >= '0' && *ccr <= '9' && ccr >= ccp ; ccr--) ;
		      for (x = 0, ccr++ ; *ccr != ':' ; ccr++)
			x = 10 * x + (*ccr - '0') ;
		      oo[noo++] = x ;
		    }
		}
	      
	      /* start from the first coord of the alignment and eat 2 more bases for each oo */
	      i = dx = 0 ;
	      while (i < noo && a1 + dx >= oo[i])
		{ i++ ; dx += 2 ; }
	      a1 += dx ;
	      i = dx = 0 ;
	      while (i < noo && r1 + dx >= oo[i])
		{ i++ ; dx += 2 ; }
	      r1 += dx ;

	      /* start from the first coord of the alignment and eat 2 more bases for each oo */
	      i = noo - 1 ; dx = 0 ;
	      while (i >= 0 && a2  - dx <= oo[i])
		{ i-- ; dx += 2 ; }
	      a2 -= dx ;
	      i = noo - 1 ; dx = 0 ;
	      while (i >= 0 && r2  - dx <= oo[i])
		{ i-- ; dx += 2 ; }
	      r2 -= dx ;
	    }
	}
      if (Remove_inserts_shorter_than)
	{ 
	  int i, deltaPair = 0 ;
	      
	  /* position on col 22 = pair fate */
	  for (i = 18 ; i<= 21 ; i++)
	    { aceInStep (ai, '\t') ; aceInWord (ai) ; }
	  aceInStep (ai, '\t') ; aceInInt (ai, &deltaPair) ;
	  if ( Remove_inserts_shorter_than && deltaPair <  Remove_inserts_shorter_than && deltaPair > - Remove_inserts_shorter_than)
	    continue ;
	}
      if (zzp && r1 > zzp->a2) r1 = zzp->a2 ;
      if (zzp && r2 > zzp->a2) r2 = zzp->a2 ;
      r1 -= zzda ; r2 -= zzda ; if (r1 < 1) r1 = 1 ; if (r2 < 0) r2 = 0 ;
      strand = isDown ? 1 : -1 ;
      for (i = r2, ip = arrayp (cover, 2*i + 0, int), x = x2 + strand * (r2 + zzda - a2) ; i >= r1 ; ip -= 2, i--, x -= strand)
	(*ip) += count * 100 ; /* 100 -> runQ (runQuality, runQualityPrefix, x) but we did not compensate for indels */
      ip =  arrayp (cover, 0, int) ;
      if (i > 0 && (*ip == 0 || *ip > i)) *ip = i ;

      /* S cover, base i is protected by at least the overhang */
      s1 = a1 ;
      s2 = a2 ; 
      if (zzp && s1 > zzp->a2) s1 = zzp->a2 ;
      if (zzp && s2 > zzp->a2) s2 = zzp->a2 ;
      s1 -= zzda ; s2 -= zzda ; if (s1 < 1) s1 = 1 ; if (s2 < 0) s2 = 0 ;
      /* run in the direction of the read to minimize the shifts in quality coordinates if there is an indel */
  
      for (i = s2, ip = arrayp (cover, 2*i + 1, int) , x = x2 + strand * (s2 + zzda - a2) ; i >= s1 ; ip -= 2, i--, x -= strand)
	(*ip) += count * 100 ; /* 100 -> runQ (runQuality, runQualityPrefix, x) ; */
      ip =  arrayp (cover, 0, int) ;
      if (i > 0 && (*ip == 0 || *ip > i)) *ip = i ;

      s1 += zzda ; s2 += zzda ;
      if (! ccp)
	continue ;
      if (0) memset (snpBuf, 0, 999) ;
      strncpy (snpBuf, ccp, 998) ; /* snip */

      if (ali > 0)
	{ nBpAli += ali ; ali = 1 ; }
      else
	{ nBpAli -= ali ; ali = -1 ; }

      nAli ++ ;
      cp = snpBuf + strlen(snpBuf) ; *++cp = 0 ;
      if (strcmp (snpBuf, "-"))
	{
	  if (strstr(snpBuf,"oo")) 
	    snpCountSnpsOO (snpBuf, xsnpBuf, s1, s2) ; 
	  cp = snpBuf ; cc = *cp ;
	  cpx = xsnpBuf ; hasMore = TRUE ;
	  while (hasMore)
	    { /* break on ',' */
	      cq = cp ; cqx = cpx ;
	      while (*cq && *cq != ',') cq++ ;
	      cc = *cq ; 
	      if (cc) 
		{ hasMore = TRUE ; *cq = 0 ; }
	      else
		hasMore = FALSE ;

	      
	      while (*cqx && *cqx != ',') cqx++ ;
	      *cqx = 0 ;

	      cr = cp ; pos = 0 ;
	      while (*cr && *cr != ':') 
		pos = 10 * pos + (*cr++ - '0') ;
	      *cr++ = 0 ;
	      cp = cr ;
	      while (*cr && *cr != '#')  cr++ ;

	      crx = cpx ; x = 0 ;
	      while (*crx && *crx != ':') 
		x = 10 * x + (*crx++ - '0') ;
	      *crx++ = 0 ;
	      cpx = crx ;
	      while (*crx && *crx != '#')  crx++ ;

	      qual = 100 ;
	      if (*cr == '#')
		{
		  *cr++ = 0 ;
		  qual = 0 ;
		  while (*cr>='0' && *cr <= '9')
		    { qual = 10 * qual + (*cr - '0') ; cr++ ; }
		}
	      if (strstr(cp,":"))
		invokeDebugger() ;
	      if (*cp == '*') cp++ ; 
	      if (1) { char *cp1 = cp ; while ((*cp1)) { if(*cp1 == 'r') *cp1 = 'a' ; if (*cp1 == 'y') *cp1 = 't' ; cp1++ ; }}
	      dictAdd (eeDict, cp, &sn) ;
	      cp = cq + 1 ;
	      cpx = cqx + 1 ;

	      q = 100 ; /* 100 -> runQ (runQuality, runQualityPrefix, x) ; */
	      ne += q ;
	      nme += count * q ;

	      if (pos < a1  || pos > a2 )
		{ /* 2010_11_23: we compensate a glitch
		   * If the errors were starting with a double delete, clipalign
		   * was reporting the error exactly at the edge of the alignment 
		   * skip those cases
		   */
		  ne-- ; nme -= count * q ;
		}
	      else if (!zp ||  (zp->a1 + zp->a2 == 0) || (pos >= zp->a1 && pos <= zp->a2)) /* acceptable pos */
		{
		  /* accumulate en place, so RAM needed is independant of the number of reads */
		  sprintf (posBuffer, "%c%d.%d", strand > 0 ? 'p' : 'm', pos - zzda,sn) ;
		  posDict = array (posDictArray, targetZone, DICT*) ;
		  if (! posDict)
		    posDict = array (posDictArray, targetZone, DICT*) = dictHandleCreate (100, h) ;
		  dictAdd (posDict, posBuffer, &posIndex) ;
		  snpsArray = strand > 0 ? snpsDownArray : snpsUpArray ;
		  snps = array (snpsArray, targetZone, Array) ;
		  if (! snps)
		    snps = array (snpsArray, targetZone, Array) = arrayCreate (1000, SS) ;
		  ss = arrayp (snps, posIndex, SS) ;
		  ss->pos = pos - zzda ;
		  ss->snp = sn ;
		  if (sn < 0) invokeDebugger() ;
		  if (ss->pos == 73126786) invokeDebugger() ;
		  ss->qual = qual ; 
		  if (! ss->count) snpTotal++ ;
		  if (sn == ooE1 || sn == ooE2 || sn == ooE3)
		    ss->count +=  q * count ;
		  else
		    ss->count +=   q * count ;
		  ss->strand = strand ;
		}
	    }
	}
    }

  nmembl =  messAllocMaxStatus (&memsize) ;
  if (! snp->ventilate)
    fprintf (stderr, "%s snps parsed, found %d snp, dictMax(eeDict)=%d, messalloc %d Mb in %d blocks\n"
	     , timeShowNow(), snpTotal, dictMax(eeDict), memsize, nmembl) ;
  for (strand = -1 ; strand < 2 ; strand += 2)
    {
      snpsArray = strand > 0 ? snpsDownArray : snpsUpArray ;
	
      for (targetZone = 0 ; targetZone < arrayMax (snpsArray) ; targetZone++)
	{	      
	  snps = array (snpsArray, targetZone, Array) ;
	  if (snps)
	    arraySort (snps, ssOrder) ;
	}
    }
  nmembl =  messAllocMaxStatus (&memsize) ;
  fprintf (stderr, "%s snps parsed and sorted, found %d snp, messalloc %d Mb in %d blocks\n",  timeShowNow(), snpTotal, memsize, nmembl) ;

  /* construct the quality profile for each type of error */
  /* in the same loop rescale down by 100, because of snpQ(runQuality) */
  for (strand = 1 ; strand >= -1 ; strand -= 2)
    {
      snpsArray = strand > 0 ? snpsDownArray : snpsUpArray ;
      
      for (targetZone = 0 ; targetZone < arrayMax (snpsArray) ; targetZone++)
	{	      
	  snps = array (snpsArray, targetZone, Array) ;
	  if (! snps) continue ;
	  
	  arrayp (snps, arrayMax (snps), SS)->snp = 0 ; /* so we loop once more */
	  for (i = oldmrna = oldpos = oldsn = 0 ; i < arrayMax (snps) ; i++)
	    { 
	      ss = arrp (snps, i, SS) ;
	      ss->count /= 1 ; /* rescale down by 100, because of snpQ(runQuality) */
	      if (ss->count)
		{
		  int bsnp = 0 ;
		  char buf[12] ;
		  char *b1, bb ;
		  const char *b2 ;
		  if (ss->strand > 0)
		    bsnp = ss->snp ;
		  else
		    {
		      bsnp = 0 ;
		      b1 = buf ;
		      b2 = dictName (eeDict, ss->snp) - 1 ;
		      while (*++b2)
			{
			  switch (*b2)
			    {
			    case 'a': bb = 't' ; break ;
			    case 't': bb = 'a' ; break ;
			    case 'g': bb = 'c' ; break ;
			    case 'c': bb = 'g' ; break ;
			    default: bb = *b2 ; break ;
			    }
			  *b1++ = bb ;
			}
		      *b1 = 0 ;
		      if (*buf) 
			dictAdd (eeDict, buf, &bsnp) ;
		      if (bsnp == 1000)
			invokeDebugger() ;
		    }
		  if (bsnp < 0) invokeDebugger() ;
		  if (bsnp)
		    {
		      ee = arrayp (errorsP, 256*bsnp + ss->qual, SS) ;
		      ee->count += ss->count ;
		      ee->qual = ss->qual ;
		      ee->snp = bsnp ;
		    }
		}
	    }
	}
    }
  nmembl =  messAllocMaxStatus (&memsize) ;
  fprintf (stderr, "%s max errorsP = %d dictMax(eeDict)=%d, messalloc %d Mb in %d blocks\n",  timeShowNow(), arrayMax(errorsP), dictMax(eeDict), memsize, nmembl) ;
  if (0)
    for (i = 1 ; i<= dictMax(eeDict) ; i++)
      fprintf (stderr, "%d::%s\n",i,dictName(eeDict, i)) ;
  if (1) /* abandon the strand and quality and condensate the snps */
    for (strand = 1 ; strand >= -1 ; strand -= 2)
      {
	int i, j, k ;
	SS *ss, *ss2, *ss3 ;
	
	snpsArray = strand > 0 ? snpsDownArray : snpsUpArray ;
	
	for (targetZone = 0 ; targetZone < arrayMax (snpsArray) ; targetZone++)
	  {	      
	    snps = array (snpsArray, targetZone, Array) ;
	    if (! snps) continue ;
	
	    for (i = k = 0, ss = ss3 = arrp (snps, 0, SS) ; i < arrayMax (snps) ; ss++, i++)
	      { 
		if (ss->count == 0) continue ;
		pos = ss->pos ; sn = ss->snp ; strand = ss->strand ;
		for (j = i + 1, ss2 = ss + 1;  j < arrayMax (snps) && ss2->pos == pos && ss2->strand == strand && ss2->snp == sn ; ss2++, j++)
		  {
		    ss->count += ss2->count ;
		    ss2->count = 0 ;
		  }
		if (! snp->qual) 
		  ss->qual = 100 ; /* abandon the strand and quality */
		if (k<i) *ss3 = *ss ;
		k++ ; ss3++ ;
	      }
	    arrayMax (snps) = k ;
	  }
      }

  
  snpExportErrors (snp, nAli, nBpAli, errorsP, eeDict) ;

  aceOutf (ao, "## Attention: the BRST counts are multiplied by up to 100: the a posteriori quality\n") ;
  aceOutf (ao, "# Target\tCoordinate\tBase (W=wild_type)\tStrand\tRun\tB\tR\tS\tT\n") ;
  
  if (1) /* at each position report the errors and then the wild type */
    {
      int *ncp = 0, iss, ii, nS, nW ;
      SS *ss0 ;
      const char *bb ;
      char wild ;
      int strand, noo1, noo2, noo3 ;
      char strandChar ;

      for (strand = 1 ; strand >= -1 ; strand -= 2)
	{
	  covers = strand > 0 ? coversDownArray : coversUpArray ;
	  snpsArray = strand > 0 ? snpsDownArray : snpsUpArray ;

	  for (targetZone = 0 ; targetZone < arrayMax (covers) ; targetZone++)
	    {
	      cover = arr (covers, targetZone, Array) ;
	      if (! cover || ! arrayMax (cover)) continue ;
	      
	      snps = array (snpsArray, targetZone, Array) ;
	      if (! arrayExists (snps))
		continue ;

	      target = keySet (targetZone2target, targetZone) ;
	      zzda = keySet (zzDas, targetZone) ;

	      /* jump the non S-cover bases */
	      pos = arr (cover, 0, int) ;  ncp = arrp (cover, 1 + 2 * pos, int) ; 
	      ss0 = arrp (snps, 0, SS), iss = 0 ;
	      for ( ; pos <= arrayMax (cover)/2 ; ncp += 2, pos++) 
		{
		  nS = 0 ;
		  /* position at this coordinate */
		  while (iss > 0 && (iss >=arrayMax (snps)  || ss0->pos > pos))
		    { iss-- ; ss0-- ;}
		  for ( ; iss < arrayMax (snps) && (ss0->pos < pos) ; ss0++, iss++) ;
 
		  noo1 = noo2 = noo3 = 0 ; /* zero the T count */
		  if (ss0->pos == pos)
		    {  
		      /* count the oo */
		      for (ss = ss0, ii = iss ; (ii < arrayMax (snps)) && ss->pos == pos && ss->strand == strand ; ss++, ii++) 
			{
			  /*
			    if (pos + zzda == 2933832) fprintf(stderr, "****** 2933832, ii=%d ss->count=%d snp=%d\n", ii, ss->count, ss->snp) ;
			  */
			  if (ss->count && ss->snp == ooE1)
			    {
			      noo1 = ss->count ;
			    }
			  if (ss->count && ss->snp == ooE2)
			    {
			      noo2 = ss->count ;
			    }
			  if (ss->count && ss->snp == ooE3)
			    {
			      noo3 = ss->count ;
			    }
			  
			    if (pos + zzda == 220958625)
			    fprintf(stderr, "****** 1,2,3 = %d # %d # %d\n", noo1, noo2, noo3) ;
			  
			}
		      /* report the errors */
		      for (ss = ss0, ii = iss ; (ii < arrayMax (snps)) && (ss->pos == pos) && ss->strand == strand ; ss++, ii++) 
			if (ss->count)
			  {
			    if (ss->snp == ooE1 || ss->snp == ooE2 || ss->snp == ooE3)
			      continue ;
			    bb = dictName (eeDict, ss->snp) ;
			    if (strstr (bb, ">n")) { noo3++ ; continue ; }
			    aceOutf (ao, "%s\t%d\t"
				     , target == oldXTarget ? "~" : dictName (snp->targetDict, target) 
				     , pos + zzda
				     ) ;
			    oldXTarget = target ;
			    
			    if (oldXSnp != ss->snp && bb[0] == '+')
			      {
				int iiw, wild ;
				wild = 'W' ;
				if (snp->target2dna && target < arrayMax (snp->target2dna) &&
				    (iiw = arr (snp->target2dna, target, int)) &&
				    iiw + pos < stackMark (snp->targetDna))
				  wild = ace_lower (*(stackText (snp->targetDna, iiw + pos + zzda - 1))) ;
				aceOutf (ao, "%c", wild) ;
			      }

			    aceOutf (ao, "%s"
				     , oldXSnp == ss->snp ? "~" : bb
				     ) ;
			    oldXSnp = ss->snp ;

			    strandChar =  strand > 0 ? '+' : '-' ;
			    aceOutf (ao, "\t%c\t%s" 
				     , strandChar
				     , oldXRun == 99999 ? "~" : snp->run
				     ) ;
			    oldXRun = 99999 ;

			    /* 
			       if (bb[1] == '>') type = 'S' ;
			       else if (bb[0] == '*') type = 'R' ;
			       else type = 'T' ;
			    */
			    aceOutf (ao, "\t%d\t%d\t%d\t%d"
				     , ss->count
				     , *(ncp - 1)   /* - noo R cover */
				     , *(ncp)       /* - noo S cover */
				     , noo2 + noo3   /* T cover */
				     ) ;
			     /* a substitution or an indel is a non-support of the wild type */
			    nS += ss->count ;
			    
			    aceOutf (ao, "\n") ;
			  }
		    }
		  /* export the wild type */
		  nW = *ncp > nS + noo1 + noo2 + noo3 ? *ncp - nS - noo1 - noo2 - noo3 : 0 ;
		  if (nW == 0)
		    continue ;
		  wild = 'W' ;
		  if (snp->target2dna && target < arrayMax (snp->target2dna) &&
		      (ii = arr (snp->target2dna, target, int)) &&
		      ii + pos < stackMark (snp->targetDna))
		    wild = ace_lower (*(stackText (snp->targetDna, ii + pos + zzda - 1))) ;
		  aceOutf (ao, "%s\t%d"
			   , target == oldXTarget ? "~" : dictName (snp->targetDict, target) 
			   , pos + zzda
			   ) ;
		  oldXTarget = target ;
		  aceOutf (ao, "\t%c%s\t%c"
			   , wild
			   , snp->qual ? ":100" : ""
			   , strand > 0 ? '+' : '-'
			   ) ;
		  oldXSnp = 999999 ;

			    aceOutf (ao, "\t%s"
				     , oldXRun == 99999 ? "~" : snp->run
				     ) ;
			    oldXRun = 99999 ;

		  
		  aceOutf (ao, "\t%d\t%d\t%d\t%d\n"
			   , nW 
			   , *(ncp - 1), *(ncp), noo2 + noo3  /* R S T cover */
			   ) ;
		}
	    }
	}
    }
  
  ac_free (h) ;

  if (snp->ventilate)
    { 
      fprintf (stderr, "%s File %s ventilated into %d u lines in %d files\n", timeShowNow(), aceInFileName (ai), ventU, ventUF) ;
      fprintf (stderr, "%s File %s ventilated into %d nu lines in %d files\n", timeShowNow(), aceInFileName (ai), ventNU, ventNUF) ;
    }
  else
    fprintf (stderr, "%s Found %d errors (%d including multiplicities) in %ld bp in %d alignments\n", timeShowNow(), ne/100, nme/100, nBpAli, nAli) ;

  return nn ;
} /* snpCountSNPs */
/*  run --select TUTU/toto.zone2 --runQuality tmp/SNP/AmishExome.Quality_profile.txt --unique --count --strategy $Strategy --minAliPerCent 90 --target_class Z_genome --run Ghs99 -i TUTU/toto3.hit */
/*************************************************************************************/
/* Danielle does not like name starting with + or - */

static char *cleanSnpName (const char * ccp, int strand, int *type)
{
  static char buf[250] ;
  char buf0[250], *cp = buf ;

  if (! ccp) return 0 ;
  strcpy (buf0, ccp) ;
  if (strand < 0)
    {
      cp = buf0 - 1 ;
      while (*++cp)
	switch ((int)*cp)
	  {
	  case 'a': *cp = 't' ; break ;
	  case 't': *cp = 'a' ; break ;
	  case 'g': *cp = 'c' ; break ;
	  case 'c': *cp = 'g' ; break ;
	  }
    }

  cp = buf ; ccp = buf0 ;
  *type = 'S' ; 
  if (*ccp == '*')
    { *type = 'R' ; *cp++ = *ccp++ ; }
  if (*ccp == '+')
    { strcpy (cp, "Ins ") ; cp += 4 ; }
  if (*ccp == '-')
    { strcpy (cp, "Del ") ; cp += 4 ; }
  while (*ccp == '+' || *ccp == '-')
    ccp++ ;
  strcpy (cp, ccp) ;
  cp = buf ;
  while (*cp++)
    if (*cp == ':') *cp = 0 ;


  cp = buf - 1 ;
  while (*++cp)
    switch ((int)*cp)
      {
      case 'a': *cp = 'A' ; break ;
      case 't': *cp = 'T' ; break ;
      case 'g': *cp = 'G' ; break ;
      case 'c': *cp = 'C' ; break ;
      }
  
  return buf ;
} /* cleanSnpName */

/*************************************************************************************/
#ifdef JUNK

/* construct a signature for the quadruplet */
static int zz (SNP *snp, int run)
{
  static DICT *dict = 0 ;
  int nn ;

  if (! dict)
    dict = dictHandleCreate (1000, snp->h) ;
  dictAdd (dict, messprintf("%d", run), &nn) ;
  return nn ;
} /* zz */

/*************************************************************************************/

static long int snpCoalesce (BigArray snps)
{
  long int ii, jj, nn ;
  SSM *ssm, *ssm2, *ssm3  ;

  /* accumulate the counts in identical SNP */
  bigArraySort (snps, ssmOrderBySnp) ;
  for (ii = nn = 0,  ssm3 = ssm = bigArrp (snps, 0, SSM) ; ii < bigArrayMax (snps) ; ssm++, ii++)
    {
      if (ssm->count > 0)
	{
	  /* accumulate */
	  for (jj = ii+1 ,  ssm2 = ssm+1  ; 
	       jj < bigArrayMax (snps) && ssm2->target == ssm->target && ssm2->pos == ssm->pos && ssm2->snp == ssm->snp && ssm2->run == ssm->run ; 
	       ssm2++, jj++)
	    {
	      /* divide by 100, because the counts were scale up by 100 to use the qualities */
	      ssm->count += ssm2->count/100 ;
	      ssm2->count = 0 ;
	      ssm->coverR += ssm2->coverR/100 ;
	      ssm->coverS += ssm2->coverS/100 ;
	      ssm->coverT += ssm2->coverT/100 ;
	    }
	  if (nn < jj) *ssm3 = *ssm ;
	  nn++ ; ssm3++ ;
	}
    }
  bigArrayMax (snps) = nn  ;
  return nn ;
} /* snpCoalesce */

/*************************************************************************************/

static long int snpBRSparseFile (SNP *snp, BOOL forgetRepeats, AC_HANDLE h)
{
  int j, strand = 0, nLines = 0 ;
  int qual, count ;
  int oldXMrna = 0, oldXPos = 0, oldXStrand = 0 , oldXSnp = 0, oldXRun = 0 ;
  long int nAli = 0, nBpAli = 0 ;
  ACEIN ai = snp->ai ;
  BigArray snps ;
  Array  dna = 0 ; 
  SSM *ssm ;
  SS *ee ;
  const char *ccp ;
  long int insideCover = 0 ;
  int inside = 0, insideStart = 0, insideLast = 0, insideTarget = 0 ;
  int minCover = 100 * snp->minCover ;
  ACEOUT aoW = aceOutCreate (snp->outFileName, ".BV", snp->gzo, h) ;
  aceInSpecial (ai, "\n") ;

  snp->run_sample = bitSetCreate (1000, h) ;
  snp->eeDict = eeDictCreate (snp, h) ;
  dictAdd (snp->eeDict, "toto", 0) ; /* avoid 1 */
  snp->errorsP = arrayHandleCreate (1000, SS, h) ;
  snps = snp->snps = bigArrayHandleCreate (10000, SSM, h) ;

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (!ccp || *ccp == '#')
	continue ;
      if (!strcmp (ccp, "nAli"))
	{ 
	  aceInStep (ai, '\t') ;
	  if (aceInInt (ai, &j)) 
	    nAli += j ;
	  continue ; 
	}
      if (!strcmp (ccp, "nBpAli"))
	{ 
	  aceInStep (ai, '\t') ;
	  if (aceInInt (ai, &j)) 
	    nBpAli += j ;
	  continue ; 
	}
      if (!strcmp (ccp, "ERR"))
	{
	  int sn ;

	  if (! (ccp = aceInWord (ai))) continue ; /* type */
	  dictAdd (snp->eeDict, ccp, &sn) ; 
	  aceInStep (ai, '\t') ;
	  if (aceInInt (ai, &qual))
	    { 
	      aceInStep (ai, '\t') ;
	      if (aceInInt (ai, &count))
		{
		  if (qual > 256) messcrash ("line %s, quality = %d should not be > 256, sorry", aceInStreamLine (ai), qual) ;
		  ee = arrayp (snp->errorsP, 256 * sn + qual, SS) ;
		  ee->count += count ;
		  ee->qual = qual ;
		  ee->snp = sn ;
		}
	    }
	  continue ;
	}
      /* actual snp */
      ssm = bigArrayp (snps, bigArrayMax (snps), SSM) ;
      if (*ccp == '~' && oldXMrna) ssm->target = oldXMrna ;
      else  
	{
	  dictAdd (snp->targetDict, ccp, &(ssm->target)) ;
	  if (snp->selected8kbDict &&
	      dictFind (snp->selected8kbDict, ccp, 0))
	    bitSet (snp->selected8kbBitSet, ssm->target) ;
	  dna = array (snp->dnas, ssm->target, Array) ;
	  if (! dna)
	    dna = array (snp->dnas, ssm->target, Array) = arrayHandleCreate (1000, char, snp->h) ;
	}
      oldXMrna =  ssm->target ;

      aceInStep (ai, '\t') ;
      ccp = aceInPos (ai) ;
      if (*ccp == '~')
	{
	  aceInWord (ai) ; /* gobble it up */
	  ssm->pos = oldXPos ;
	}
      else if (aceInInt (ai, &(ssm->pos)))
	{  
	  if (ssm->pos == 0)
	    ssm->pos = oldXPos ;
	  oldXPos = ssm->pos ;
	}
      else
	{
	  messerror ("Expecting a position in column 2 , line %d, in file %s"
		     , aceInStreamLine (ai), aceInFileName (ai)) ;
	  continue ;
	}   

      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp) 
	{
	  messerror ("Expecting a SNP in column 3 , line %d, empty column, in file %s"
		     , aceInStreamLine (ai), aceInFileName(ai)) ;
	  continue ;
	}
      if (*ccp == '~' && oldXSnp) ssm->snp = oldXSnp ;
      else if (*ccp == '-' && ccp[1] == 0) 
	{ bigArrayMax (snps)-- ; continue ; } /* we are out of the genome reference */
      else if (0 && ccp [1])
	{ bigArrayMax (snps)-- ; continue ; } /* HACK: we only keep the ref and kill the variants */
      else 
	{
	  char c, c1 ;
	  int i3 ;

	  if (forgetRepeats && *ccp == '*') /* forget the * because it is different in fasta and csfasta (Solid) */
	    ccp++ ;
	  c = *ccp ;
	  if (dna)
	    {
	      if (c == '-')  /* expecting --at go and fetch the a */
		{
		  for (i3 = 1 ; c == '-' && i3 < 4 ; i3++)
		    c = *(ccp+i3) ;
		}
	      switch ((int)c)
 		{
		case 'a':
		case 't':
		case 'g':
		case 'c':
		  array (dna, ssm->pos, char) = c ;
		  break ;
		}
	    }
	  c = *ccp ;
	  ssm->snp = 0 ;
	  switch ((int) c)
	    {
	    case 0:
	      break ;
	    case '+':            /* expecting ++at  */
	    case '-':            /* expecting -a  */
	      dictAdd (snp->eeDict, ccp, &(ssm->snp)) ;
	      break ;
	    default:
	      c1 = *(ccp+1) ;
	      switch ((int) c1)
		{
		case 0:          /* expecting a  */
		case ':':        /* expecting a:34  */
		  break ;
		case '>':        /* expecting a>g  */
		  dictAdd (snp->eeDict, ccp, &(ssm->snp)) ;
		  break ;
		case '+':        /* expecting a++gt where a is wild type  */
		case '-':
		  dictAdd (snp->eeDict, ccp+1, &(ssm->snp)) ;
		  break ;
		}
	      break ;
	    }
	}
      
      oldXSnp = ssm->snp ;
      
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (!ccp)
	{
	  messerror ("Expecting a strand info in column 4 , line %d, empty column, in file %s"
		     , aceInStreamLine (ai), aceInFileName (ai)) ;
	  continue ; 
	}
      if (*ccp == '+') strand = 1 ; 
      else if (*ccp == '-') strand = -1 ;
      else if (*ccp == '~') strand = oldXStrand ;
      else
	{
	  messerror ("Expecting a strand info in column 4 , line %d, got %s, in file %s"
		     , aceInStreamLine (ai), ccp, aceInFileName (ai)) ;
	  continue ;
	}
      oldXStrand = strand ;	      
      
      /* from now on ssm->snp is negative for the reverse strand */
      ssm->snp = (1 + ssm->snp) * strand ;
      
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ; 
      if (! ccp)
	{
	  messerror ("Expecting a run in column 5, line %d, empty column, in file %s"
		     , aceInStreamLine (ai), aceInFileName (ai)) ;
	  continue ;
	}
      if (*ccp == '~' && oldXRun) ssm->run = oldXRun ;
      else  dictAdd (snp->runDict, ccp, &(ssm->run)) ;
      oldXRun = ssm->run ;

      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &(ssm->count)))
	{
	  messerror ("Expecting a count integer in column 9, line %d, empty column, in file %s"
		   , aceInStreamLine (ai), aceInFileName (ai)) ;
	  continue ;
	}

      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &(ssm->coverR)))
	{
	  messerror ("Expecting R coverage integer in column 10, line %d, empty column, in file %s"
		   , aceInStreamLine (ai), aceInFileName (ai)) ;
	  continue ;
	}

      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &(ssm->coverS))) /* S */ 
	{
	  messerror ("Expecting S coverage integer in column 11, line %d, empty column, in file %s"
		   , aceInStreamLine (ai), aceInFileName (ai)) ;
	  continue ;
	}

      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &(ssm->coverT))) /* T */ 
	{
	  messerror ("Expecting T coverage integer in column 12, line %d, empty column, in file %s"
		   , aceInStreamLine (ai), aceInFileName (ai)) ;
  	  continue ;
	}
    
      if (0) bitSet (snp->run_sample,  zz (snp, ssm->run)) ;
      nLines++ ;

      if (! (nLines % 5000000)) /* every 5M lines, condense the table */
	snpCoalesce (snps) ;
      
      if (ssm->snp * ssm->snp == 1)
	switch (inside)
	  {
	  case 1:
	    if (ssm->coverS < minCover ||  insideTarget != ssm->target)
	      {
		inside = 0 ; 
		aceOutf (aoW, "%s\t%d\t%d\t%ld\n"
			 , dictName (snp->targetDict, ssm->target)
			 , insideStart,  insideLast
			 , insideLast >= insideStart ? insideCover/(100 * (insideLast - insideStart + 1)) : 0
			 ) ;
		insideCover = 0 ;
	      }
	    else
	      {
		insideCover += ssm->coverS ; insideLast = oldXPos ;
		break ;
	      }
	    /*  fall thru, this may happen if we changed target */
	  case 0 :
	    if (ssm->coverS >= minCover)
	      {
		inside = 1 ; insideStart = oldXPos ;  insideLast = oldXPos ; insideTarget = ssm->target ; insideCover = ssm->coverS ;
	      }
	    break ;
	  }      
    } /* aceInCard */
	
  snpCoalesce (snps) ;

  if (inside == 1)
    aceOutf (aoW, "%s\t%d\t%d\t%ld\n"
	     ,  dictName (snp->targetDict, insideTarget)
	     , insideStart, insideLast
	     , insideCover/(100 * (insideLast - insideStart + 1))
	     ) ;

  fprintf (stderr, "// Parsed %ld lines in input file : %s\n", bigArrayMax (snps), timeShowNow ()) ;
  return bigArrayMax (snps) ;
} /* snpBRSparseFile */
#endif
/*************************************************************************************/

static int snpParseSnpList (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  int n, nn = 0, target = 0, sn, pos ;
  ACEIN ai ;
  char *cp, *cq ;
  BRS2 *brs2 ;
  DICT *dict ;
  DICT *targetDict = snp->targetDict ;
  BigArray aap ;

  dict = snp->oldSnps = dictHandleCreate (100000, snp->h) ;
  ai = aceInCreate (snp->snpListFileName, 0, h) ;
  if (! ai)
    messcrash ("Cannot open file -snp_list %s",  snp->snpListFileName) ;

  while (aceInCard (ai))
    {

      cp = aceInWord (ai) ;
      if (!cp || *cp == '#')
	continue ;
      if (! strcasecmp (cp, "Variant"))
	continue ;

      n = strlen(cp) ; if (n > 2047) messcrash ("bad variant column 3 file %s, line %d, name  too long > 2048 char : %s\n"
					     , aceInFileName (ai), aceInStreamLine (ai), cp
					     ) ;
      cq = strstr (cp, ":") ;
      if (! cq)
	messcrash ("bad variant (no : between the sequence name and the position in column 3 file %s, line %d : %s\n"
		   , aceInFileName (ai), aceInStreamLine (ai), cp
		   ) ;
      dictAdd (dict, cp, &sn) ;
      *cq++ = 0 ;
      dictAdd (targetDict, cp, &target) ;
      aap = array (snp->aaaP, target, BigArray) ;
      if (! aap)
	aap = array (snp->aaaP, target, BigArray) = bigArrayCreate (10000, BRS2) ;

      cp = cq ;
      cq = strchr (cq, ':') ; 
      if (cq) *cq++ = 0 ;
      if (sscanf (cp, "%d", &pos) != 1 || pos < 1)
	messcrash ("bad variant (no : between the position and the type in column 3 in file %s, line %d : %s\n"
		   , aceInFileName (ai), aceInStreamLine (ai), dictName (dict, sn)
		   ) ;
      brs2 = bigArrayp (aap, pos, BRS2) ;
      brs2->oldSnp = sn ;
    }
  ac_free (h) ;

  nn = dictMax (dict) ;
  fprintf (stderr, "%s found %d known SNPs in file %s\n", timeShowNow(), nn, snp->snpListFileName) ;
  return nn ; 
} /* snpParseSnpList */

/*************************************************************************************/
#ifdef JUNK

static BOOL  snpBRS2snpName (SNP *snp, SSM *ssm, char* namBuf, char *typeBuf, char *Abuf, char *Bbuf, char *tagBuf, int *deltap, const char *requestedName) 
{
  int sn, sign = 0, d = 0 ;
  BOOL isRepeat = FALSE ;
  char buf0[1024], *cp , hookBase = 'n' ;
  const char *ccp ;

  sn = ssm->snp ; if (sn < 0) sn = - sn ; sn-- ;
  ccp = sn ? dictName (snp->eeDict, sn) : "W" ;
  
  if (requestedName)  /* deduce all data from the previous name */
    {
      char *cp, *cq ;
      
      if (SNPDEBUG)
	fprintf(stderr, "X2 %s\n", ccp) ;
      sprintf (namBuf, "%s", requestedName) ;
      strncpy (buf0, requestedName, 1000) ;
      cp = strchr (buf0, ':') ;
      cp = strchr (cp+1, ':') ;
      sprintf (typeBuf, "%s", cp+1) ;
      strncpy (buf0, typeBuf, 1000) ;
      cp = buf0 ;
      cq = strchr (cp, ':') ;
      
      if (SNPDEBUG)
	fprintf(stderr, "X3 %s\n", cp) ;
      if (! strncmp (cp, "Sub", 3)) 
	{
	  *deltap = 0 ;
	  Abuf[0] = tagBuf[0] = cp[4] ;
	  tagBuf[1] = '>' ;
	  Bbuf[0] = tagBuf[2] = cp[6] ;
	  Abuf[1] = Bbuf[1] = tagBuf[3] = 0 ;
	  isRepeat = FALSE ;
	  if (SNPDEBUG)
	    fprintf(stderr, "X6 tagBuf=%s\n", tagBuf) ;
	}
      else if (cq && ! strncmp (cp, "Del", 3)) 
	{
	  char *cr ;
	  int i = 0 ;

	  *cq++ = 0 ;
	  if (strchr (cp, '/'))
	    isRepeat = TRUE ;
	  cr = strchr (cq, ':') ;
	  if (cr) *cr++ = 0 ;
      if (SNPDEBUG)
	fprintf(stderr, "X4 cp=%s cq=%s cr=%s\n", cp, cq, cr) ;
	  cq++ ;  /* jump :anchorBase */
	  cr = Abuf ;
	  while (++i < 12 && *cq && *cq == ace_upper (*cq))
	    *cr++ = *cq++ ;
	  *cr = 0 ;
	  *deltap = strlen (Abuf) ;
	  Bbuf[0] = '.' ; Bbuf[1] = 0 ;
	  cp = tagBuf ;
	  if (isRepeat)
	    *cp++ = '*' ;
	  *cp++ = '-' ;
	  strcpy (cp, Abuf) ;
	  if (SNPDEBUG)
	    fprintf(stderr, "X6 tagBuf=%s\n", tagBuf) ;
  	}
      else if (cq && ! strncmp (cp, "Ins", 3)) 
	{
	  char *cr ;
	  int i = 0 ;

	  *cq++ = 0 ;
	  if (strchr (cp, '/'))
	    isRepeat = TRUE ;
	  cr = strchr (cq, ':') ;
	  if (cr) *cr++ = 0 ;
	  if (SNPDEBUG)
	    fprintf(stderr, "X4 cp=%s cq=%s cr=%s\n", cp, cq, cr) ;
	  cq = cr;  /* move to the interesting part */
	  cq++ ;  /* jump :anchorBase */
	  cr = Bbuf ;
	  while (++i < 12 && *cq && *cq == ace_upper (*cq))
	    *cr++ = *cq++ ;
	  *cr = 0 ;
	  *deltap = strlen (Bbuf) ;
	  Abuf[0] = '.' ; Abuf[1] = 0 ;
	  cp = tagBuf ;
	  if (isRepeat)
	    *cp++ = '*' ;
	  *cp++ = '+' ;
	  strcpy (cp, Bbuf) ;
      if (SNPDEBUG)
	fprintf(stderr, "X6 tagBuf=%s\n", tagBuf) ;
	  
	}
      cp = Abuf - 1 ; while (*++cp) *cp = ace_lower(*cp) ;
      cp = Bbuf - 1 ; while (*++cp) *cp = ace_lower(*cp) ;
      cp = tagBuf - 1 ; while (*++cp) *cp = ace_lower(*cp) ;
      
      if (SNPDEBUG)
	fprintf(stderr, "X7 Abuf=%s Bbuf=%s tagBuf=%s namBuf=%s delta=%d\n", Abuf, Bbuf, tagBuf, namBuf, *deltap) ;

      if (Abuf[11]) messcrash ("Abuff 11") ;
      if (Bbuf[11]) messcrash ("Bbuff 11") ;
      return isRepeat ;
    }

  if (sn || SNPDEBUG)
    fprintf(stderr, "Y1 %s\n", ccp) ;
  if (sn)
    {
      strcpy (buf0, ccp) ;
      cp = buf0 - 1 ;
      while (*++cp)
	switch ((int)*cp)
	  {
	  case 'a': *cp = 'A' ; break ;
	  case 't': *cp = 'T' ; break ;
	  case 'g': *cp = 'G' ; break ;
	  case 'c': *cp = 'C' ; break ;
	  case '2': *cp = '>' ; break ;
	  }
      
      typeBuf[0] = 0 ;
      if (buf0[1] == '>')
	{
	  sprintf (typeBuf, "Sub:%c:%c", buf0[0], buf0[2]) ;
	  Abuf[0] = buf0[0] ; Abuf[1] = 0 ; 
	  Bbuf[0] = buf0[2] ; Bbuf[1] = 0 ; 
	}
      
      *deltap = 0 ; 
      ccp = buf0 ;
      
      if (*ccp == '*')
	{ ccp++ ; isRepeat = TRUE ; }
      
      /* we have ++aa if from BRS file and InsAA if from requested list */
      if (*ccp == '+')
	{ strcpy (typeBuf, "Ins") ; sign = 1 ; }
      else if (*ccp == '-')
	{  strcpy (typeBuf, "Del") ; sign = -1 ;  }
      else if (! strncmp (ccp, "Ins", 3))
	{   ccp += 3 ; d = strlen (ccp) ; strcpy (typeBuf, "Ins") ; sign = 1 ;}
      else if (! strncmp (ccp, "Del", 3))
	{ ccp += 3 ; d = strlen (ccp) ;  strcpy (typeBuf, "Del") ; sign = -1 ;  }
      while (*ccp == '+' || *ccp == '-')
	{ d++ ; ccp++ ; }
      *deltap = sign * d ;  /* number of bases added or substracted */
      
      if (sign >= 1)  {  Abuf[0] = '.' ; Abuf[1] = 0 ; strncpy (Bbuf, ccp, 3) ; Bbuf[3] = 0 ; }
      else if (sign <= -1)  {  Bbuf[0] = '.' ; Bbuf[1] = 0 ; strncpy (Abuf, ccp, 3) ; Abuf[3] = 0 ; }
      else { Abuf[0] = tagBuf[0] = ccp[0] ;  Abuf[1] = 0 ; tagBuf[1] = '>' ; Bbuf[0] = tagBuf[2] = ccp[2] ;  Bbuf[1] = 0 ; tagBuf[3] = 0 ; }
      if (Abuf[11]) messcrash ("Abuff 11") ;
      cp = Abuf ; 
      while (*cp) 
	{ *cp = ace_lower(*cp) ; cp++ ; }
      if (Abuf[11]) messcrash ("Abuff 11") ;
      cp = Bbuf  ; 
      while (*cp) 
	{ *cp = ace_lower(*cp) ; cp++ ; }
      cp = tagBuf ; 
      while (*cp) 
	{
	  *cp = ace_lower(*cp) ;
	  cp++ ;
	}
      
      
      if (sign == 1)
	{
	  Array dna = array (snp->dnas, ssm->target, Array) ;
	  char s1[1024], s2[1024], slide[32] ;
	  
	  s1[0] = s2[0] = 'n' ; 
	  s1[1] = s2[1] = 0 ; 
	  slide[0] = 0 ;
	  /* roll back */
	  if (dna && ssm->pos > 1 && ssm->pos - 2 < arrayMax (dna))
	    {
	      int i, j, dx = 0 ;
	      while (arr (dna, ssm->pos - 1 - dx, char)  == ace_lower(buf0[d + ((10000 * d - dx + d-1) % d)]))
		dx++ ;
	      ssm->leftShift = dx ;
	      hookBase = arr (dna, ssm->pos - dx - 1, char) ;
	      s1[0] = s2[0] = hookBase ;
	      i = j = 1 ;
	      for (j = 0 ; j < d && j < 1000 ; j++)
		s2[j+1] = ace_upper(Bbuf[(1000*d +j - dx)%d]) ;
	      j++ ; /* because of the hookBase */
	      for (i = 1 ; i < dx + 2 && i < 1000 && j < 1000 ; j++, i++)
		s1[i] = s2[j] = arr (dna, ssm->pos - dx + i - 1, char) ;
	      s1[i] = s2[j] = 0 ;
	      if (1)
		{
		  int i ;
		  char buf[d+1] ;
		  char *cp = tagBuf ;
		  
		  if (dx) 
		    { 
		      isRepeat = TRUE ;
		      *cp++ = '*' ;
		    }
		  *cp++ = '+' ;
		  for (i = 0 ; i < d ; i++)
		    *cp++ = buf[i] = Bbuf[(1000*d +i - dx)%d] ;
		  *cp++ = buf[i] = 0 ;
		  memcpy (Bbuf, buf, d+1) ;
		  for (i = 0 ; i < d ; i++)
		    buf0[d+i] = ace_upper(buf[i]) ;
		  if (dx) sprintf (slide, "/%d", d+dx) ;
		}
	    }
	  if (d == 1 && ! slide[0])
	    sprintf (typeBuf, "Ins%s:%s:%s", slide, s1, s2) ;
	  else 
	    sprintf (typeBuf, "Ins_%d%s:%s:%s", d, slide, s1, s2) ;
	}
      else if (sign == -1)
	{
	  Array dna = array (snp->dnas, ssm->target, Array) ;
	  char s1[1024], s2[1024], slide[32] ;
	  int dx = 0 ;
	  
	  s1[0] = s2[0] = 'n' ; 
	  s1[1] = s2[1] = 0 ; 
	  slide[0] = 0 ;
	  /* roll back */
	  if (dna && ssm->pos > 1 && ssm->pos  < arrayMax (dna))
	    {
	      int i, j, dy ;
	      int a = ssm->pos - 1 ;
	      
	      while (arr (dna, a - dx, char)  == Abuf[(10000 * d - dx + d - 1) % d]) 
		dx++ ;
	      ssm->leftShift = dx ;
	      hookBase = arr (dna, ssm->pos - dx - 1, char) ;
	      dy = 0 ; 
	      while (arr (dna, ssm->pos - dx + dy + d, char)  == Abuf[(10000 * d - dx + dy ) % d]) 
		dy++ ;
	      s1[0] = s2[0] = hookBase ;
	      i = j = 1 ;
	      for (i = 0 ; i < d && i < 1000 ; i++)
		s1[i+1] = ace_upper(Abuf[i]) ;
	      i++ ; /* because of the hookBase */
	      for (j = 1 ; j < dy + 2 && i < 1000 && j < 1000 ; j++, i++)
		s1[i] = s2[j] = arr (dna, ssm->pos - dx + d + j - 1, char) ;
	      s1[i] = s2[j] = 0 ;
	      if (1)
		{
		  int i ;
		  char *cp = tagBuf ;
		  
		  if (dy)
		    {
		      isRepeat = TRUE ;
		      *cp++ = '*' ;
		    }
		  *cp++ = '-' ;
		  for (i = 0 ; i < d ; i++)
		    *cp++ = Abuf[i] = arr (dna, ssm->pos - dx + i, char) ;
		  *cp++ = Abuf[i] = 0 ;
      if (Abuf[11]) messcrash ("Abuff 11") ;
		  if (dy) sprintf (slide, "/%d", d+dy) ;
		}
	    }
	  if (d == 1 && !slide[0])
	    sprintf (typeBuf, "Del%s:%s:%s", slide, s1, s2) ;
	  else 
	    sprintf (typeBuf, "Del_%d%s:%s:%s", d, slide, s1, s2) ;
	}
      sprintf (namBuf, "%s:%d:%s", dictName (snp->targetDict, ssm->target),  ssm->pos - ssm->leftShift, typeBuf) ;
    }
  else
    sprintf (namBuf, "%s:%d:", dictName (snp->targetDict, ssm->target),  ssm->pos) ;
      if (Abuf[11]) messcrash ("Abuff 11") ;
      if (Bbuf[11]) messcrash ("Bbuff 11") ;

  return isRepeat ;  
} /* snpBRS2snpName */
#endif

/*************************************************************************************/
/* from wolfram: p = proba (mean - alpha*sigma < x < mean + alpha*sigma)
 *  1-p     alpha
 * 5/100 = 1.95996 sigma
 * 1/100 = 2.57583 sigma
 * 10-3  = 3.29053 sigma
 *
 * sigma = sqrt (1/N sum(x - <x>)^2)
 *    or we replace N by N-1 if <x> is estimated from the distrib 
 *
 *  0x1  : not_low  Pooled data, the a priori proba is {0 to 100%, step 1}
 *  0x2  : not intermediate
 *  0x4  : not_high
 *
 *  0x01  : not_ww   Non pooled data, the a priori proba is {0, 50, 100%}
 *  0x02  : not wm
 *  0x04  : not_mm
 */

#define MX 110
static int snpNotLow (int count, int cover, double *zp, double *z1p, double *z2p, int limit)
{
  static int firstPass = 0 ;
  double cumul, cumulLow, cumulMid, cumulHigh ;
  static int isLow[MX] ;
  static int isHigh[MX] ;
  static int isMid1[MX] ;
  static int isMid2[MX] ;
  static int isww[MX] ;
  static int ismm[MX] ;
  static int iswm1[MX] ;
  static int iswm2[MX] ;
  double alpha, beta, risk ;
  double z, dz, z1, z2 ;
  int i, pp, result = 0 ;

  if (cover < 0) messcrash ("cover < 0") ;
  if (count < 0) messcrash ("count < 0") ;
  /* alpha if the usual 1.96 of 1.96sigma, beta is used as 3/n for risk 5% */
  risk = 5.0/100 ; alpha = 1.95996 ;  beta = 3 ;
  risk = 1.0/1000 ; alpha = 3.29053 ;  beta = 5 ;
  risk = 1.0/100 ; alpha = 2.57583 ;  beta = 4 ;

#define CC(_n,_p) cc[MX*(_n)+(_p)]
#define CCC(_pp,_n,_p) ccc[MX*MX*(_pp)+MX*(_n)+(_p)]

  if (! firstPass)
    {
      firstPass = 1 ;
      /* compute Pascal's C(n,p) coefficients */
      int n, p ;
      double c, cc[MX * MX], ccc[101 * MX * MX], z, z1, z2, z21 ;

      memset(isLow,0, sizeof(isLow)) ;
      memset(isMid1,0, sizeof(isMid1)) ;
      memset(isMid2,0, sizeof(isMid2)) ;
      memset(isHigh,0, sizeof(isHigh)) ;

      memset(isww,0, sizeof(isww)) ;
      memset(iswm1,0, sizeof(iswm1)) ;
      memset(iswm2,0, sizeof(iswm2)) ;
      memset(ismm,0, sizeof(ismm)) ;

      memset(cc,0, sizeof(cc)) ;
      memset(ccc,0, sizeof(ccc)) ;
      
      for (n = 0 ; n < MX ; n++)
	{
	  c = CC(n,0) = CC(n,n) = 1 ;
	  for (p = 1 ; 2*p <= n ; p++)
	    CC(n,p) = CC(n, n - p) = c = c * (n- p + 1) / p ;

	  if (0 && n < 10)
	    {
	      fprintf(stderr, "Pascal n=%d ::", n) ;
	      for (p = 0 ; p <= n ; p++)
		fprintf(stderr, "\t%g", CC(n,p)) ; 
	      fprintf(stderr, "\n") ;
	    }
	}
      
      /* compute for every proba pp the binomial distribution */
      for (pp = 0 ; pp <= 100 ; pp++)
	{
	  z1 = pp/100.0 ; z2 = 1 - z1 ; z21 = z2 > 0 ? z1/z2 : 1 ;
	  for (n = 0 ; n < MX ; n++)
	    {
	         /* compute z2^n */
	      for (i = 1, z = 1 ; i <= n ; i++) 
		z *= z2 ;
	      CCC(pp,n,0) = z ; 
	      if (0) fprintf (stderr, "proba(%d,%d | %d) = %g\n", 0,n,pp,CCC(pp,n,0)) ;
	      for (p = 1 ; p < n ; p++)
		{
		  z *= z21 ;
		  CCC(pp,n,p) = z * CC(n,p) ;
		  if (0) fprintf (stderr, "proba(%d,%d | %d) = %g\n", p,n,pp,CCC(pp,n,p)) ;
		}
	         /* compute z2^n */
	      for (i = 1, z = 1 ; i <= n ; i++) 
		z *= z1 ;
	      CCC(pp,n,n) = z ;
	      if (0) fprintf (stderr, "proba(%d,%d | %d) = %g\n", n,n,pp,CCC(pp,n,n)) ;
	    }
	}
	
      /* cumulated distrib
       * foreach choice of p,q
       * I want to know if 
       * sum(pp=0 to 10) ccc(pp,n,p) exceeds or the risk
       * if not then i can conclude no_low
       * so for a given n, the first such value for p is the
       * threshold we want to find
       */
	 
      /* pooled samples : low high mid */
      for (n = 0 ; n < MX ; n++)
	{
	  for (p = 0 ; p <= n ; p++)
	    {      
	      cumul = cumulLow = cumulMid = 0 ;
	      for (pp = 0 ; pp <= 100 ; pp++)
		{
		  cumul += CCC(pp,n,p) ;
		  if (pp <= 5) cumulLow += CCC(pp,n,p) ;
		  if (pp <= 95) cumulMid += CCC(pp,n,p) ;
		}
	      if (cumulLow > risk * cumul) /* this value of p is too small to be not_low */
		isLow[n] = p + 1 ;
	      if (cumulMid > risk * cumul) /* this value of p is too small to be not_mid upper bound */
		isMid2[n] = p + 1 ;
	      if (0) fprintf (stderr, "=== n=%d p=%d  cumulLow=%g  cumul=%g  isLow=%d\n", n, p, cumulLow, cumul, isLow[n]) ;
	    }
	  for (p = n  ; p >= 0 ; p--)
	    {      
	      cumul = cumulMid = cumulHigh = 0 ;
	      for (pp = 0 ; pp <= 100 ; pp++)
		{
		  cumul += CCC(pp,n,p) ;
		  if (pp >= 95) cumulHigh += CCC(pp,n,p) ;
		  if (pp > 5) cumulMid += CCC(pp,n,p) ;
		}
	      if (cumulHigh > risk * cumul) /* this value of p is too big to be not_high */
		isHigh[n] = p - 1 ;
	      if (cumulMid > risk * cumul) /* this value of p is too big to be not_mid */
		isMid1[n] = p - 1 ;
	    }
	}
    
      /* individual samples : ww wm mm */
      for (n = 0 ; n < MX ; n++)
	{
	  for (p = 0 ; p <= n ; p++)
	    {      
	      cumul = cumulLow = cumulMid = 0 ;
	      for (pp = 0 ; pp <= 100 ; pp+=1)
		{
		  if ((pp > 2 && pp < 49) || (pp > 51 && pp < 98)) continue ;
		  cumul += CCC(pp,n,p) ;
		  if (pp <= 2) cumulLow += CCC(pp,n,p) ;
		  if (pp <= 97) cumulMid += CCC(pp,n,p) ;
		}
	      if (cumulLow > risk * cumul) /* this value of p is too small to be not_low */
		isww[n] = p + 1 ;
	      if (cumulMid > risk * cumul) /* this value of p is too small to be not_mid upper bound */
		iswm2[n] = p + 1 ;
	      if (0) fprintf (stderr, "=== n=%d p=%d  cumulLow=%g  cumul=%g  isLow=%d\n", n, p, cumulLow, cumul, isww[n]) ;

	    }
	  for (p = n  ; p >= 0 ; p--)
	    {      
	      cumul = cumulMid = cumulHigh = 0 ;
	      for (pp = 0 ; pp <= 100 ; pp+=1)
		{
		  if ((pp > 2 && pp < 49) || (pp > 51 && pp < 98)) continue ;
		  cumul += CCC(pp,n,p) ;
		  if (pp >= 98) cumulHigh += CCC(pp,n,p) ;
		  if (pp >= 3) cumulMid += CCC(pp,n,p) ;
		}
	      if (cumulHigh > risk * cumul) /* this value of p is too big to be not_high */
		ismm[n] = p - 1 ;
	      if (cumulMid > risk * cumul) /* this value of p is too big to be not_mid */
		iswm1[n] = p - 1 ;
	    }
	}
    
      if (0)
	{
	  fprintf (stderr, "n\tLow\tM1\tM2\tHigh\tz1\tz2(low) at risk %g\n", risk) ;
	  for (n = 1 ; n < MX ; n++)
	    {
	      fprintf (stderr, "%d\t%d\t%d\t%d\t%d"
		       , n
		       , isLow[n]
		       , isMid1[n]
		       , isMid2[n]
		       , isHigh[n]
		       ) ;
	      if (n)
		{
		  p = isLow[n] ;

		  z = (p)/(1.0*n) ; if (z < 0) z = 0 ;if (z > 1) z = 1 ;
		  dz = sqrt(z * (1-z)/n) ;
		  fprintf (stderr, "\t%.2f\t%.2f", 100 * (z - alpha * dz), 100 * (z + alpha * dz)) ;
		}
	      fprintf (stderr, "\n") ;
	    }
	}
    }

  if (cover == 0)
    { z = dz = 0 ; }
  else
    {
      z = (double)count/cover ;
      if (z > 1) z = 1 ;
      dz = sqrt(z * (1-z)/cover) ;
    }
  /* limits of confidence interval at chosen risk */
  z1 = z - alpha * dz ;
  z2 = z + alpha * dz ;
  if (dz == 0 && cover > 0 && z == 1)
    z1 = 1 - beta/cover ;
  if (dz == 0 && cover > 0 && z == 0)
    z2 = 0 + beta/cover ;
  if (z1 < 0) z1 = 0 ;
  if (z2 > 1) z2 = 1 ;
  
  if (zp)
    {
      *zp  = 100 * z  ;
      *z1p = 100 * z1 ;
      *z2p = 100 * z2 ;
    }
  if (cover < 4 || cover < limit) return 0 ;
  
  /* return resonable guesses */
  result = 0 ;
  if (count == 0)   
    result = 0x2 | 0x4 | 0x20 | 0x40 ;  /* wild type */
  if (count == 1)
    {
      if (cover >= 10)
	result |= 0x2 | 0x4 | 0x20 | 0x40 ;  /* wild type */
      else
	result |= 0x4 | 0x40 ;  /* not high */
    }
  if (cover - count == 0)   
    result |= 0x2 | 0x1 | 0x20 | 0x10 ;  /* mm */
  if (cover - count == 1)
    {
      if (cover >= 10)
	result |= 0x2 | 0x1 | 0x20 | 0x10 ;  /* mm */
      else
	result |= 0x1 | 0x10 ;  /* not low */
    }
  if (cover - count == 2)
    {
      if (cover >= 12)
	result |= 0x2 | 0x1 | 0x20 | 0x10 ;  /* mm */
      else
	result |= 0x1 | 0x10 ;  /* not low */
    }
  if (cover >= MX)
    {
      /*we  compute limits of the normal Gaussian interval */
      if (100 * z1 > 2.5)
	result |= 0x1 ;   /* not_low (lower bound < 1.5%) */
      if (0 && z > beta/cover)
	result |= 0x1 ;   /* not_low (z > upper variation from zero) this does not work we get a case in NBE 30/50,000 = 0.1% counted as not_low */
      if (100 * z2 < 95)
	result |= 0x4 ;   /* not_high */
      if (100 * z1 > 65 || 100 * z2 < 35)
	result |= 0x2 ;   /* not_mid */

      if (100 * z1 > 2.5)
	result |= 0x10 ;  /* not_ww */
      if (0 && z > beta/cover)
	result |= 0x10 ;   /* not_low (z > upper variation from zero) */
      if (100 * z2 < 95)
	result |= 0x40 ;   /* not_mm */
      if (100 * z1 > 65 || 100 * z2 < 35)
	result |= 0x20 ;   /* not_wm */
    }
  else if (cover < MX)
    {  
      if (0)
	{
	  if (count >= isLow[cover])  /* not_low */
	    result |= 0x1 ;
	  if (count <= isHigh[cover])  /* not_high */
	    result |= 0x4 ;
	  if (count <= isMid1[cover] || count >= isMid2[cover])  /* not_mid */
	    result |= 0x2 ;
	}
      /* we adopt the discrete values */
      if (count >= isww[cover])  /* not_low */
	result |= 0x10 | 0x1 ;
      if (count <= ismm[cover])  /* not_high */
	result |= 0x40 | 0x4 ;
      if (count <= iswm1[cover] || count >= iswm2[cover])  /* not_mid */
	result |= 0x20 | 0x2 ;
    }

  return result ;
} /* snpNotLow */

/*************************************************************************************/

static void snpNotLowShow (ACEOUT ao, int count, int called, double *zp, double *z1p, double *z2p, int limit)
{
  float dosage = -1 ;
  int dosage2 = -1 ;


  int notLow = snpNotLow (count, called, zp, z1p, z2p, limit) ;

  if ((notLow & 0x01) && (notLow & 0x02) && (notLow & 0x04)) dosage = (*zp < 40 ? 0.5 : 1.5) ; 
  else if ((notLow & 0x02) && (notLow & 0x04)) dosage2 = 0 ;
  else if ((notLow & 0x01) && (notLow & 0x04)) dosage2 = 1 ;
  else if ((notLow & 0x01) && (notLow & 0x02)) dosage2 = 2 ;
  else if (notLow & 0x01) dosage = 1.7 ;
  else if (notLow & 0x04) dosage = 0.3 ;
  
  if (dosage + dosage2 < -1.9) 
    aceOut (ao, "\t-") ;
  else if (dosage2 < 0)
    aceOutf (ao,  "\t%.1f", dosage) ;
  else
    aceOutf (ao, "\t%d", dosage2) ;
  
  return ;
} /* snpNotLowShow */

/*************************************************************************************/
#ifdef JUNK
static int ssmOrderByC2a (const void *va, const void *vb)
{
  const SSM *a = (const SSM *)va, *b = (const SSM *)vb ;
  int n, s1, s2 ;

  n = a->target - b->target ; if (n) return n ;
  n = a->pos - b->pos ; if (n) return n ;
  s1 = a->snp ; if (s1<0) s1 = -s1 ;
  s2 = b->snp ; if (s2<0) s2 = -s2 ;
  n = a->run - b->run ; if (n) return n ;
  n = s1 - s2 ; if (n) return n ; /* get wild type first */
  n = a->snp - b->snp ; if (n) return -n ; /* strand before antistrand */

  return 0 ;
} /* ssmOrderByC2a */

/*************************************************************************************/
/*************************************************************************************/

static BOOL snpIsCompatibleToGenome (SNP *snp, SSM *ssm, int leftShift, const char *type, const char *Abuf)
{
  Array  dna = snp->dnas ? array (snp->dnas, ssm->target, Array) : 0 ;
  const char *ccp ;

  if (dna)
    {
      ccp = arrp (dna, ssm->pos - leftShift, char) ;

      if (! strncmp (type, "Ins", 3))
	return TRUE ;
      else if (! strncmp (type, "Del", 3))
	{
	  while (*Abuf)
	   {
	     if (*ccp && ace_lower (*ccp) != ace_lower(*Abuf))
	       return FALSE ;
	     ccp++ ; Abuf++ ;
	   }
	}
      else if (*ccp && ace_lower (*ccp) != *Abuf)
	return FALSE ;
    }
  return TRUE ;
} /* snpIsCompatibleToGenome */

/*************************************************************************************/
/* find the counts on other strand at same position, same snp or wild
 * return FALSE is the corresponding point is above (already exported)
 * or if the requested snp does not match the genome
 * divide by 100 all reported counts to compensate the runQ a posteriori quality factor
 *
 * for deletions the wild count is R 
 * for insertions the wild count is R (like del)  of previous base if I insert the same base tttgg-> tttTgg
 * or local S (like for a substitution) if I inser something else      tttgg -> tttAgg
 */
static BOOL snpCountOtherStrand (SNP *snp, SSM *ssm0, int sn, long int ii0, int delta, int leftShift, const char *typeBuf, const char *Abuf
				 , int *cover, int *calledp, int *calledm, int *mp, int *mm, int *wp, int *wm
				 , int *oo1p, int *oo1m, int *oo2p, int *oo2m, int pass)
{
  SSM *ssm, *ssm2 ;
  int coverp = 0, coverm = 0 ;
  int pos = ssm0->pos, target = ssm0->target, run = ssm0->run ;
  long int ii, ii2, iMax = bigArrayMax (snp->snps) ;
  int b = 0, donep, donem ;

  if (sn < 0) sn = -sn ; 
  *cover = *calledp = *calledm = *mp = *mm = *wp = *wm = *oo1p = *oo1m = *oo2p = *oo2m = 0 ;

  if (! snpIsCompatibleToGenome (snp, ssm0, leftShift, typeBuf, Abuf))
    return FALSE ;

  /* move up */
  for (ssm = ssm0 - 1, ii = ii0 - 1 ; ii >= 0 && ssm->pos == pos && ssm->target == target && ssm->run == run; ssm--, ii--)
    if (ssm->snp == sn || ssm->snp == -sn) return FALSE ; /* so we export only once */
  /* if solid, we need the previous oo value */
  if (snp->solid)
    for (ssm2 = ssm, ii2 = ii ; ii2 >= 0 && ssm2->pos >= pos - 1 && ssm2->target == target && ssm2->run == run; ssm2--, ii2--)
      {
	if (ssm2->snp > 0) *oo1p = ssm2->coverT ;
	if (ssm2->snp < 0) *oo1m = ssm2->coverT ;
      }

  if (delta > 0) /* we need to know what base was inserted to compare it to the previous base */
    {
      const char* ccp = dictName (snp->eeDict, sn - 1) ;

      while (*ccp == '*') ccp++ ;    
      while (*ccp == '+') ccp++ ;
      if (!strncmp (ccp, "Ins", 3)) ccp += 3 ;
      b = ace_lower(*ccp) ;
    }
  if (delta < 0) /* we need to know what base was inserted to compare it to the previous base */
    {
      const char* ccp = dictName (snp->eeDict, sn - 1) ;

      while (*ccp == '*') ccp++ ;    
      while (*ccp == '-') ccp++ ;
      if (!strncmp (ccp, "Del", 3)) ccp += 3 ;
      b = ace_lower(*ccp) ;
    }

  donep = donem = 0 ;
  for (ii++, ssm++ ; ii < iMax && ssm->pos == pos && ssm->target == target && ssm->run == run ; ssm++, ii++) 
    {
      /* with the else, the counts are ok even if we question on a wild-type ssm0 */

      if (delta == 0)       /* substitution: mutant, wild type and coverage are directly available */
	{
	  if (ssm->snp >= 1)  { coverp = ssm->coverS ;   *oo2p = ssm->coverT ; }
	  if (ssm->snp <= -1)  { coverm = ssm->coverS ;  *oo2m = ssm->coverT ; }
	  
	  if (ssm->snp >= 1) *calledp += ssm->count ;
	  else *calledm += ssm->count ;
	  
	  if (ssm->snp == sn)        *mp = ssm->count ;
	  else if (ssm->snp == -sn)  *mm = ssm->count ; 
	  
	  if (ssm->snp == 1)        *wp = ssm->count ;
	  else if (ssm->snp == -1)  *wm = ssm->count ;
	  
	  /* the insertions also support the wild type */
	  if (ssm->snp > snp->eeLastSub) /* we are on the plus strand */
	    {
	      const char* ccp = dictName (snp->eeDict, ssm->snp - 1) ;
	      
	      while (*ccp == '*') ccp++ ;    
	      if (*ccp == '+')
		{
		  *wp += ssm->count ;
		}
	    }
	  if (ssm->snp <  - snp->eeLastSub) /* we are on the minus strand */
	    {
	      const char* ccp = dictName (snp->eeDict, -ssm->snp - 1) ;
	      
	      while (*ccp == '*') ccp++ ; 
	      if (*ccp == '+')
		{
		  *wm += ssm->count ;
		}
	    }
	}
      else if (delta < 0)    /* deletion: mutant is directly available, cover = coverR of this letter
			      *  locate the stretch [r1,r2] of letter identical to r0 in the reference, r0 is somewhere inside [r1, r2]
			      *  wild = cases with correct number of letters = cover - any delete in [r1, r2] - any insert at r2+1
			      */
	{
	  if (! donep && ssm->snp == sn)
	    {
	      int r1, r2 ;
	      Array dna = array (snp->dnas, ssm->target, Array) ;
	      r1 = r2 = ssm->pos ;
	      if (dna)
		{
		  /* locate the stretch [r1,r2] in the reference with the same base b as the base beeing deleted */
		  const char *ccp = arrp(dna, ssm->pos - 2, char) ;
		  while (*ccp-- == b) r1-- ;
		  ccp = arrp(dna, ssm->pos + 1, char) ;
		  while (*ccp++ == b) r2++ ;
		}
	      
	      /* count in r1, r2 the number of cases with indels on the required stretch, add the inserts of the same base at the next position */
	      
	      for (ssm2 = ssm, ii2 = ii ; ii2>=0 && ssm2->pos >= r1 && ssm2->target == target && ssm2->run == run ; ssm2--, ii2--) ;
	      if (ssm2->pos < r1) { ssm2++ ; ii2++ ; }
	      for ( ; ii2 < iMax && ssm2->pos <= r2 + 1 ; ssm2++, ii2++)
		{
		  if (ssm2->pos == r2 + 1)
		    { /* deduce the insertions of the same base from the wild count */
		      if (ssm2->snp >  snp->eeLastSub) /* we are on the plus strand */
			{
			  const char* ccp = dictName (snp->eeDict, ssm2->snp - 1) ;
			  
			  while (*ccp == '*') ccp++ ;    
			  if (*ccp == '+')
			    {
			      while (*ccp == '+') ccp++ ;
			      if (!strncmp (ccp, "Ins", 3)) ccp += 3 ;
			      if (b == ace_lower(*ccp))
				*wp -= ssm2->count ;
			    }
			}
		      if (ssm2->snp <  - snp->eeLastSub) /* we are on the minus strand */
			{
			  const char* ccp = dictName (snp->eeDict, -ssm2->snp - 1) ;
			  
			  while (*ccp == '*') ccp++ ; 
			  if (*ccp == '+')
			    {
			      while (*ccp == '+') ccp++ ;
			      if (b == ace_lower(*ccp))
				*wm -= ssm2->count ;
			    }
			}
		      
		      continue ;
		    }
		  if (ssm2->snp > snp->eeLastSub || ssm2->snp < - snp->eeLastSub) 
		    {
		      if (ssm2->pos > r1)
			{
			  if (ssm2->snp > 0) *wp -= ssm2->count ;
			  else *wm -=  ssm2->count ;
			}
		      if (ssm2->pos == r1) /* only accept a deletion of base b. not an insert of a previous base */
			{
			  const char* ccp = dictName (snp->eeDict, (ssm2->snp > 0 ? ssm2->snp - 1 : -ssm2->snp - 1)) ;
			  
			  while (*ccp == '*') ccp++ ; 
			  if (*ccp == '-') /* in principle we do not need to check which base it is */
			    {
			      while (*ccp == '-') ccp++ ;
			      if (b == ace_lower(*ccp))
				{
				  if (ssm2->snp > 0) *wp -= ssm2->count ;
				  if (ssm2->snp < 0) *wm -= ssm2->count ;
				}
			    }
			}
		    }
		  if (ssm2->snp == sn && ssm2->pos == pos)
		    {
		      if (sn > 0) *mp = ssm2->count ;
		      else  *mm = ssm2->count ;
		    }
		  if (ssm2->snp == - sn && ssm2->pos == pos)
		    {
		      if (sn < 0) *mp = ssm2->count ; 
		      else  *mm = ssm2->count ; 
		    }
		  if (ssm2->snp > 0) { coverp = ssm2->coverR ;  *oo2p = ssm->coverT ; }  /* for cover, we autaumatically catch the rigthmost base b */
		  if (ssm2->snp < 0) { coverm = ssm2->coverR ;  *oo2m = ssm->coverT ; } /* for cover, we autaumatically catch the rigthmost base b */
		}
	      *wp += coverp ;
	      *wm += coverm ;
	      *calledp = coverp ; 
	      *calledm = coverm ;
	      donep = 1 ;
	    }
	  if (pass == 1 && ! donep && ssm->pos == pos && (ssm->snp == 1 ||  ssm->snp == -1))
	    {
	      if (ssm->snp == 1)   { coverp = ssm->coverR ;  *calledp += ssm->coverR ; *wp = ssm->coverR ; *oo2p = ssm->coverT ; }
	      if (ssm->snp == -1)  { coverm = ssm->coverR ;  *calledm += ssm->coverR ;  *wm = ssm->coverR ; *oo2m = ssm->coverT ; }
	      donep = 1 ;
	    }
	}
      else if (delta > 0)    /* insertion: mutant is directly available, 
			      *  locate the stretch [r1,r2] of letter identical to the inserted letter in the reference, most often r0 = r2 + 1
			      * if the inserted letter is different from the letter before, cover = soverS[r0]
			      * if the inserted letter is identical to the previou one (duplication) cover = coverR[r2]
			      * where [r1,r2] is the stretch of identical letters
			      *  wild = cases with correct number of letters = cover - any insert or delete in [r1, r2] - any insert at r2+1
			      */
	{
	  if (! donep && ssm->snp == sn)
	    {
	      int r1, r2 ;
	      Array dna = array (snp->dnas, ssm->target, Array) ;
	      r1 = r2 = ssm->pos - 1 ;
	      if (dna)
		{
		  /* insertion at r0: mutant is directly available
		   * if dup, insertion of same letter as previous one
		   *   locate that strectch [r1,r2], cover = coverR of [r1,r2] (should be a constant), wild = cover - del anywhere in [r1,r2] - any inset at r0
		   * if insert at r0, instertion of a ditinct letter
		   *   cover = coverS at r0, wild = cover - any del at r0 - any insert at r0
		   * if a coverR - delete here or before in TTT stratch - insert at nect base */
		  
		  /* locate the stretch [r1,r2] in the reference with the same base b as the base beeing deleted */
		  const char *ccp = arrp(dna, ssm->pos - 1, char) ;
		  if (* ccp == b)
		    { 
		      /* this insertion is really a duplication, we must study the [r1,r2] stretch */
		      int i2 = ii ;
		      while (i2-- > 0 && *--ccp == b) r1-- ;
		    }
		  else
		    {
		      /* this is a new letter, counts are directly available */
		      r1 = r2 = ssm->pos ;
		      for (ssm2 = ssm, ii2 = ii ; ii2>=0 && ssm2->pos >= r1 && ssm2->target == target && ssm2->run == run ; ssm2--, ii2--) ; 
		      if (ssm2->pos < r1) { ssm2++ ; ii2++ ; }
		      for ( ; ii2 < iMax && ssm2->pos <= r2  ; ssm2++, ii2++)
			{
			  if (ssm2->snp == sn)   /* direct support for this insert */
			    {
			      if (sn > 0) { *mp = ssm2->count ; *oo2p = ssm->coverT ; }
			      else { *mm = ssm2->count ; *oo2m = ssm->coverT ; }
			    }
			  else if (ssm2->snp == - sn)
			    {
			      if (sn > 0) { *mm = ssm2->count ; *oo2m = ssm->coverT ; }
			      else { *mp =  ssm2->count ; *oo2p = ssm->coverT ; }
			      
			    }
			  if (ssm2->snp > 0) { coverp = ssm2->coverS ;  *oo2p = ssm->coverT ; }  /* for cover, we autaumatically catch the rigthmost base b */
			  if (ssm2->snp < 0) { coverm = ssm2->coverS ;  *oo2m = ssm->coverT ; } /* for cover, we autaumatically catch the rigthmost base b */
			  
			  if (ssm2->snp > snp->eeLastSub) *wp -=  ssm2->count ;	  /* any other indel must be deduced from the wild type */
			  else if (ssm2->snp < - snp->eeLastSub) *wm -=  ssm2->count ;
			}
		      *wp += coverp ;
		      *wm += coverm ;
		      *calledp = coverp ; 
		      *calledm = coverm ;
		      donep = 1 ;
		      continue ; /* this simple case is resolved */
		    }	
		}
	      
	      /* count in r1, r2 the number of cases with indels on the required stretch, add the inserts of the same base at the next position */
	      
	      for (ssm2 = ssm, ii2 = ii ; ii2>=0 && ssm2->pos >= r1 && ssm2->target == target && ssm2->run == run ; ssm2--, ii2--) ;
	      if (ssm2->pos < r1) { ssm2++ ; ii2++ ; }
	      r2 = pos ;
	      for ( ; ii2 < iMax && ssm2->pos <= r2 ; ssm2++, ii2++)
		{
		  if (ssm2->pos == pos)
		    { 
		      if (ssm2->snp == sn)   /* direct support for this insert */
			{
			  if (sn > 0) { *mp = ssm2->count ; *oo2p = ssm->coverT ; }
			  else { *mm = ssm2->count ; *oo2m = ssm->coverT ; }
			}
		      else if (ssm2->snp == - sn)
			{
			  if (sn > 0) { *mm = ssm2->count ; *oo2m = ssm->coverT ; }
			  else { *mp =  ssm2->count ; *oo2p = ssm->coverT ; }
			}
		      
		      /* deduce the insertiton of the same base from the wild count */
		      if (ssm2->snp >  snp->eeLastSub) /* we are on the plus strand */
			{
			  const char* ccp = dictName (snp->eeDict, ssm2->snp - 1) ;
			  
			  while (*ccp == '*') ccp++ ;    
			  if (*ccp == '+')
			    {
			      while (*ccp == '+') ccp++ ;
			      if (b == ace_lower(*ccp))
				*wp -= ssm2->count ;
			    }
			}
		      if (ssm2->snp <  - snp->eeLastSub) /* we are on the minus strand */
			{
			  const char* ccp = dictName (snp->eeDict, -ssm2->snp - 1) ;
			  
			  while (*ccp == '*') ccp++ ; 
			  if (*ccp == '+')
			    {
			      while (*ccp == '+') ccp++ ;
			      if (b == ace_lower(*ccp))
				*wm -= ssm2->count ;
			    }
			}
		      
		      continue ;
		    }
		  if (ssm2->snp > snp->eeLastSub || ssm2->snp < - snp->eeLastSub) /* count the insert or the deletes in the strech as non wild */
		    {
		      if (ssm2->pos > r1)  /* inserts and deletes both count as non wild-type number of letters */
			{
			  if (ssm2->snp > snp->eeLastSub ) *wp -= ssm2->count ;
			  if (ssm2->snp < - snp->eeLastSub ) *wm -= ssm2->count ;
			}
		      if (ssm2->pos == r1) /* only accept a deletion of base b. not an insert of a previous base */
			{
			  const char* ccp = dictName (snp->eeDict, (ssm2->snp > 0 ? ssm2->snp - 1 : -ssm2->snp - 1)) ;
			  
			  while (*ccp == '*') ccp++ ; 
			  if (*ccp == '-') /* in principle we do not need to check which base it is */
			    {
			      while (*ccp == '-') ccp++ ;
			      if (b == ace_lower(*ccp))
				{
				  if (ssm2->snp > 0) *wp -= ssm2->count ;
				  if (ssm2->snp < 0) *wm -= ssm2->count ;
				}
			    }
			}
		    }
		  if (ssm2->snp > 0) coverp = ssm2->coverR ;   /* for cover, we autaumatically catch the rigthmost base b */
		  if (ssm2->snp < 0) coverm = ssm2->coverR ;   /* for cover, we autaumatically catch the rigthmost base b */
		}
	      *wp += coverp ;
	      *wm += coverm ;
	      *calledp = coverp ; 
	      *calledm = coverm ;
	      donep = 1 ;
	    }
	  if (pass == 1 && ! donep && ssm->pos == pos && (ssm->snp == 1 ||  ssm->snp == -1))
	    {
	      for (ssm2 = ssm, ii2 = ii ; ii2 >= 0 && ssm2->pos >= pos - 1 ; ssm2--, ii2--) 
		if (ssm2->pos ==  pos - 1)
		  {
		    if (ssm2->snp == 1)   { coverp = ssm2->coverR ;  *calledp = ssm2->coverR ; *wp = ssm2->coverR ; *oo2p = ssm2->coverT ; }
		    if (ssm2->snp == -1)  { coverm = ssm2->coverR ;  *calledm = ssm2->coverR ;  *wm = ssm2->coverR ; *oo2m = ssm2->coverT ; }
		  }
	      donep = 1 ;
	    }
	  
	}
    }


  if (coverp < 0) coverp = 0 ;
  if (coverm < 0) coverm = 0 ;
  if (*calledp < 0) *calledp = 0 ;
  if (*calledm < 0) *calledm = 0 ;
  if (*mp < 0) *mp = 0 ;
  if (*mm < 0) *mm = 0 ;
  if (*wp < 0) *wp = 0 ;
  if (*wm < 0) *wm = 0 ;

  if (*oo1p < 0) *oo1p = 0 ;
  if (*oo2p < 0) *oo2p = 0 ;
  if (*oo1m < 0) *oo1m = 0 ;
  if (*oo2m < 0) *oo2m = 0 ;

  *cover = coverp + coverm ; 
  
  *cover /= 100 ; *calledp  /= 100 ; *calledm  /= 100 ; 
  *mp  /= 100 ; *mm  /= 100 ; *wp  /= 100 ; *wm  /= 100 ; 
  *oo1p  /= 100 ; *oo1m  /= 100 ; *oo2p  /= 100 ; *oo2m  /= 100 ; 

  return *cover > 0 ? TRUE : FALSE ;
} /* snpCountOtherStrand */

/*************************************************************************************/

static int snpBRS2snpExport (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = snp->ao ;
  BigArray snps = snp->snps ;
  DICT *oldSnps = snp->oldSnps ;
  DICT *targetDict = snp->targetDict ;
  SSM *ssm ;
  /*
   * The call 
     ~/ace/bin.ICC_centos7/snp -BRS_count -run RNA_ROCR1_A_1 -strategy RNA_seq -i tmp/SNP_BRS/RNA_ROCR1_A_1/zoner.5.BRS.u.gz -snp_list Global_SNP_LIST/SnpA2R2.zoner.5.Variant.list | grep DDR2.aAug10 | grep 1336
     * crashes if i declare Abuf[12], Bbuf[12] and also for Abuf[24], Bbuf[24] but not if 24 12, mieg 2023_06_03
     * no idea why , valgrind gives a hint
     */

  char namBuf[2048], typeBuf[1024], Abuf[12], Bbuf[12], tagBuf[1000] ;
  int i, ii, nn = 0, cover = 0, mutant, calledp, calledm, wild, delta = 0 ;
  int mp, mm, wp, wm, oo1p, oo1m, oo2p, oo2m ;
  double z, z1, z2 ;
  float cz2 ;
  int iOldSnp = 0 ;
  long int coveredPos[1001] ;
  BitSet bb1, bb[128] ;
  /*   KEYSET goodRuns = keySetHandleCreate (h) ; */
  int prefix = 6 ;   /* how many base we want to export on each side of the event in the SNPnet */
  char bufD[256], bufS[256], *ss ;
  Array dna ;
  BOOL showTitle = TRUE ;
  /* variables used to compute the cumulated histogram of the transcript coverage, in order to evaluate the 3' bias */
  int wCountTranscripts = 0, wKeptTranscripts = 0 ;
  /*  
      KEYSET wCover = keySetHandleCreate (h) ; 
      int wCumul = 0 ;
      int oldTarget = 0, wMaxPos = 0 , wMaxValue = 0, wCountTranscripts = 0 ;
  */
  KEYSET wGlobalCover = keySetHandleCreate (h) ;
  KEYSET wGlobalShiftedCover = keySetHandleCreate (h) ;
  int wStep = 10 ;
  ACEOUT wao = 0 ;
  BOOL wDebug2 = FALSE ;

  memset (tagBuf, 0, sizeof (tagBuf)) ;
  memset (typeBuf, 0, sizeof(typeBuf)) ;
  memset (namBuf, 0, sizeof(namBuf)) ;

  if (0 &&  /* the wiggles are transfered back to bestali->mrnaSupport */
      ! oldSnps &&  snp->strategy == STRATEGY_RNA_SEQ)  /* detect mode, RNA_seq mode, compute the cover histogram to allow a better evaluation of the 3' bias */
    {
      wao = aceOutCreate (snp->outFileName, ".3pHisto.txt", FALSE, h) ;
      
      aceOutDate (wao, "###", snp->run) ;
      aceOutf (wao, "## Number of contributing transcripts\t%d\tCumul\t%ld\n", wCountTranscripts) ;
      aceOutf (wao,
	       "## The maxima of the contributing transcripts, usually representing the major polyA addition site, is further than 5kb from the 5' end of the annotated transcript.\n"
	       "## This implies that the histogram below 5kb represents the 3' biais, and above 5kb a commbination of the 3'biaias and the prevalence of very long transcripts\n"
	       "## The coverage plots were piled up, aligning the maxima at position. The x coordinates run 3' to 5' and represent the distance to the 3' end\n"
	       "# Run\tDistance from 3' end"
	       ) ;
      for (i = 0 ; i <= 16000 ; i+=10)
	aceOutf (wao, "\t%d", i) ;
    }

  memset (coveredPos, 0, sizeof (coveredPos)) ;

  if (showTitle)
    {
      aceOutDate (ao, "###", snp->project) ;
      aceOutf (ao, 
	       "##  1: Target :: Reference chromosome or transcript.\n"
	       "##  2: Pos :: Position in the reference of the hook base, \n"
	       "##      -all our coordinates start at 1 (not zero).\n"
	       "##      -all bases correspond to the plus strand of the reference\n"
	       "##      -for a substitution,the position of the modified base is given\n" 
	       "##      -for an indel, the position of the hook (last conserved) base is given\n" 
	       "##  3: Type :: Sub,Del,Ins : reference : modified seq, starting on hook base.\n"
	       "##  4: Reference :: the local reference sequence.\n"
	       "##  5: Variation :: the local variant sequence:\n"
	       "##      -We show, in upper case, or as -, the modified letters\n"
	       "##      -For a substitution, the modified base sits at the declared position\n"
	       "##      -For an insertion or a deletion, the first - sits at the declared position\n"
	       "##      -Usually, we show, in lower case, %d bases on each side of the event\n"
	       "##       but for an insertion or a deletion, the situation is more complex.\n"
	       "##       Consider the deletion-GA 2000:tGAGAc -> tGA--c,\n"
	       "##       where 2000 is the position of the second G in the reference.\n"
	       "##       The same deletion could as accurately be called \n"
	       "##       deletion-GA 1998:t--GAc or deletion-AG 1999:tG--Ac \n"
	       "##       So we show in upper case the whole GAGA sliding region\n"
	       "##  6: Run :: Name of the experiment as declared in the Meta-database\n"
	       "##  7: Dose :: Allelic dosage: two doses denote an homozyguous m/m variant,\n"
	       "##       one dose, an m/+ heterozyguous, zero dose, a +/+ wild type\n"
	       "##       When the counts are low or the frequency far from 0/50/100%%\n"
	       "##       we say 0.5 or 1.5 to indicate our incertitude\n"
	       "##  8: Fraction :: The percentage of reads supporting the variant\n"
	       "##  9: Cover :: The total number of unambiguous reads covering the event.\n"
	       "##      -For an indel the read must bridge across the variation\n"
	       "##       For example if we delete an A in a  gAAAAt stretch,\n"
	       "##       we count the reads bridging from g to t across the AAAA\n"
	       "##       or from g to c in delete GA 2000:tGAGAc -> tGA--c\n"
	       "## 10: V= :: The total number of reads supporting the variant\n"
	       "## 11: R= :: The total number of reads supporting the reference\n"
	       "##      -The coverage may be bigger than the sum V+R if some reads\n"
	       "##       support a different variation at the same position\n"
	       "## 12: V+ :: The total number of reads supporting the variant\n"
	       "##       aligning on the plus strand of the reference\n"
	       "## 13: R+ :: The total number of reads supporting the plus\n"
	       "##       strand of the reference\n"
	       "## 14: V- :: The total number of reads supporting the variant\n"
	       "##       aligning on the minus strand of the reference\n"
	       "## 15: R- :: The total number of reads supporting the minus\n"
	       "##       strand of the reference\n"
	       "## 16: Amb :: The number of ambiguous or minor allele reads at this position\n"  
	       "##       by minor we mean at a frequency lower that the %d%% threshold for this file\n"
	       "## 17: N+ :: The total number of ambiguous bases at this position\n"
	       "##       in reads aligning on the plus strand of the reference\n" 
	       "## 18: N- :: The total number of ambiguous bases at this position\n"
	       "##       in reads aligning on the minus strand of the reference\n" 
	       , prefix
	       , snp->minFrequency
	       ) ;
      
      if (snp->solid)
	aceOut (ao,
		"##      -For SOLiD sequencing, we give 2 numbers for each strand, \n"
		"##      for example gTc 2:3  1:0, to indicate on the plus strand\n"
		"##      two reads not matching the gT transition, one not matchingn Tc\n"
		"##      and on the minus strand, one read not matching gT\n"
		) ;
      aceOut (ao, 
	      "## 19: Type :: Visible on both strands, on one or Strand Incompatible\n"
	      "## 20: Chi2 :: The chi2 of the counts when comparing the 2 strands\n"
	      "## 21: N_rich :: 40%% N at this position, or on either side in Solid\n"
	      "## 22: Dose on plus strand\n"
	      "## 23: Dose on minus strand\n"
	      ) ;
    }

  aceOut (ao, "#\n#Target\tPosition\tType\tReference\tVariation\tRun\tDose\t%%\tCover\tV=\tR=\tV+\tR+\tV-\tR-\tN+\tN-\tAmb\tType\tChi2\tN_rich\tDose+\tDose-\n") ;
  
  if (oldSnps)
    {
      for (ii = 0 ; ii < 1 ; ii++)
	bb[ii] = bitSetCreate (1024, h) ;
    }

  if (!snps || ! bigArrayMax(snps))
    goto done ;

  bigArraySort (snps, ssmOrderByC2a) ;
  for (ii = 0, ssm = bigArrp (snps, 0, SSM) ; ii < bigArrayMax (snps) ; ssm++, ii++)
    {
      /*       keySet (goodRuns, ssm->run) = 1 ; */

#ifdef JUNK
      if (wao &&
	  (ssm->snp == 1 || ssm->snp == -1)                   /* just count at the wild positions */
	  ) 
	{
	  if (ssm->target != oldTarget)
	    {
	      if (wCover && wMaxPos > 0 &&  keySetMax (wCover) > 8000)
		wCountTranscripts++ ;

	      if (wCover &&  keySetMax (wCover) > 8000)
		{
		  int i, j, nn, subCumul = 0 ;

		  /* find the max of the histo under the condition that we eat at most 10% of the total */
		  wMaxPos = keySetMax (wCover) - 1 ; wMaxValue = 0 ;
		  for (nn = 0, i =  wMaxPos ; i >= 0 ; i--)
		     {
		       j = keySet (wCover, i) ;
		       nn += j ;
		       if (10 * nn < wCumul && j >  wMaxValue)
			 { wMaxPos = i ; wMaxValue = j ; subCumul = nn ; }
		     }

		  if (wDebug) /* debugging */
		    {
		      aceOutf (wao, "\n%s\t%s", snp->run, dictName (targetDict, oldTarget)) ;
		      for (i = 0 ; i <= wMaxPos && i <= 16000 ; i += wStep)
			aceOutf (wao, "\t%d",  keySet (wCover, wMaxPos - i)) ;
		    }

		  if (wMaxPos >= 8000)
		    {
		      wKeptTranscripts++ ;   

		      /* register the shifted smoothed histo */
		      for (i = 0 ;  i <= wMaxPos && i <= 16000 ; i += wStep)
			keySet (wGlobalShiftedCover, i) += keySet (wCover, wMaxPos - i) ; /* count backwards, starting at the max, which is probably close tot he 3' end of the transcript */
		      
		      if (wDebug2) /* debugging */
			{
			  wMaxPos = keySetMax (wCover) - 1 ;
			  for (i = 0 ;  i <= wMaxPos && i <= 16000 ; i += wStep)
			    keySet (wGlobalCover, i) += keySet (wCover, wMaxPos - i) ; /* count backwards, starting at the max, which is porbbaly close tot he 3' end of the transcript */
			}
		    }
		  if (wDebug2) /* debugging */
		    {
		      aceOutf (wao, "\n%s raw shift=%d cumul=%d subCumul=%d\t%s", snp->run, keySetMax (wCover) - wMaxPos, wCumul, subCumul, dictName (targetDict, oldTarget)) ;
		      wMaxPos = keySetMax (wCover) - 1 ;
		      for (i = 0 ;  i <= wMaxPos && i <= 16000 ; i+= wStep)
			aceOutf (wao, "\t%d",  keySet (wCover, wMaxPos - i)) ;
		    }

		}
	      wMaxPos = wCumul = wMaxValue = 0 ;
	      wCover = keySetReCreate (wCover) ;
	      oldTarget = ssm->target ;

	      if (wao && ssm->target && strstr(dictName (targetDict, ssm->target), ".a"))
		wMaxPos = 1 ; /* otherwise the histo will not be computed */ 
	      if (snp->selected8kbBitSet && bitt (snp->selected8kbBitSet, ssm->target))
		wMaxPos = 1 ; /* force select */		
	    }
	  if (wMaxPos)
	    {
	      int myCover = ssm->coverS/100 ;
	      keySet (wCover, ssm->pos) +=  myCover ;
	      wCumul +=  myCover ;
	      /*
		myCover = keySet (wCover, ssm->pos) ;
		if (myCover >  wMaxValue)
		{  wMaxValue = myCover ;  wMaxPos = ssm->pos ; }
	      */
	    }	  
	}
#endif

      if (ssm->snp == 1 || ssm->snp == -1) continue ;
      memset (Abuf, 0, sizeof(Abuf)) ; memset (Bbuf, 0, sizeof(Bbuf)) ;
      snpBRS2snpName (snp, ssm, namBuf, typeBuf, Abuf, Bbuf, tagBuf, &delta, 0) ;

      /* report the SNPs,
       * detect mode :  if they are over threshold 
       * count mode : if they are over requested
       */
      if (oldSnps)  /* count mode */
	{
	  if (dictFind (oldSnps, namBuf, &iOldSnp))
	    bitSet (bb[0], iOldSnp) ;                                   /* requested position */
	  else
	    continue ;
	} 

      if (! snpCountOtherStrand (snp, ssm, ssm->snp, ii, delta, ssm->leftShift, typeBuf, Abuf, &cover, &calledp, &calledm, &mp, &mm, &wp, &wm, &oo1p, &oo1m, &oo2p, &oo2m, 0))
	continue ;
      if (Abuf[11]) messcrash ("Abuff 11") ;
      mutant = mp + mm ; wild = wp + wm ;
      
      if (0)
	{  /* 2015_03_15, obsolete */
	  int mycover = mutant + wild ;
	  if (mycover > 1000) mycover = 1000 ;
	  coveredPos[mycover]++ ;
	}

      iOldSnp = 0 ;

      if (! oldSnps) /* detect mode */
       {
	 if (snp->minFrequency * (calledp + calledm) <= 100 * mutant &&
	     snp->minMutant <= mutant &&
	     snp->minCover <= (calledp + calledm)
	     )  ;                                                 /*  above thresholds  */
	 else
	   continue ;
       }

      if (1)
	{
	  aceOutf (ao, "%s\t%d\t%s"
		   , dictName (targetDict, ssm->target), ssm->pos - ssm->leftShift
		   , typeBuf
		   ) ;
	  dna = array (snp->dnas, ssm->target, Array) ;
	  memset (bufD, 0, sizeof (bufD)) ;
	  memset (bufS, 0, sizeof (bufS)) ;
	  if (dna)
	    {
	      int i, j ;
	      char *cp ;

	      for (j = 0, i = ssm->pos - ssm->leftShift - prefix ; j < prefix && i < 0 ; j++, i++)  
		bufD[j] = 'N' ;
	      for ( ; j < prefix && i < arrayMax(dna) ; j++, i++)  
		bufD[j] = arr (dna, i, char) > 0 ? ace_upper(arr (dna, i, char)) : 'N' ;
	      bufD[j] = 0 ;
	      if (delta == 0) i++ ;
	      else if (delta < 0) i -= delta ; 
	      /* if delta > 0) do not move */
	      for (j = 0 ; j < prefix && i < arrayMax(dna) ; j++, i++)  
		bufS[j] =  arr (dna, i, char) > 0 ? ace_upper(arr (dna, i, char)) : 'N' ;
	      for (; j < prefix ; j++) 
		bufS[j] = 'N' ;
	      bufS[j] = 0 ;
	      if (delta < 0)
		{
		  cp = Abuf ;
		  for (j = 0 ;  ace_lower (bufS[j]) == ace_lower (*cp) ; j++, cp = (cp < Abuf - delta - 1 ? cp + 1 : Abuf))
		    bufS[j] = ace_lower(bufS[j]) ;
		}
	      if (delta > 0)
		{
		  cp = Bbuf ;
		  for (j = 0 ; ace_lower (bufS[j]) == ace_lower (*cp) ; j++, cp = (cp < Bbuf + delta - 1 ? cp + 1 : Bbuf))
		    bufS[j] = ace_lower(bufS[j]) ;
		}
	    }
	  {
	    const char *mmm = "---" ;
	    aceOutf (ao, "\t%s%s%s\t%s%s%s"
		     , bufD, delta <= 0 ? Abuf : mmm + 3 - delta , bufS
		     , bufD, delta >= 0 ? Bbuf : mmm + 3 + delta , bufS
		     ) ;
	  }
	  aceOutf (ao, "\t%s", dictName (snp->runDict, ssm->run)) ;
	}
      if (1)
	{
	  snpNotLowShow (ao, mutant, calledp + calledm, &z, &z1, &z2, 0) ;
	  aceOutf (ao, "\t%.1f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d"
		   , z
		   , calledp + calledm, mutant, wild
		   , mp, wp, mm, wm
		   , cover - calledp - calledm > 0 ? cover - calledp - calledm : 0
		   ) ;
	  if (snp->solid) /* solid */
	    aceOutf (ao, "\t%d:%d\t%d:%d"
		     , oo1p, oo2p
		     , oo1m, oo2m
		     ) ;
	  else
	    aceOutf (ao, "\t%d\t%d"
		     , oo2p
		     , oo2m
		     ) ;
	}
      
      if (1)
	{
	  ss = "Both_strands" ; cz2 = 0 ;
	  if (mp + wp + mm + wm < snp->minCover) ss = "-" ;
	  else if (2 * (mp + wp) < snp->minCover) ss = "Minus_strand_only" ;
	  else if (2 * (mm + wm) < snp->minCover) ss = "Plus_strand_only" ;
	  else
	    { 
	      double x, y ;
	      chi2 (mp, wp, mm, wm, &cz2) ;
	      x = mp/(1.0 + mp + (double)wp) ;
	      y = mm/(1.0 + mm + (double)wm) ;
	      if (cz2 > 100 && (x > y + .5 || x < y - .5)) ss = "Incompatible_strands" ; /* extremely high and clear contradiction */ 
	    }

	  aceOutf (ao, "\t%s\t%.1f", ss, cz2) ;
	}
      aceOutf (ao, "\t%s"
	       , 100 * (oo1p + oo1m + oo2p + oo2m) > 40 * cover ? "N_rich" : "-"
	       ) ;
       if (1) /* dosage on plus and minus strand */
	{                             
	  snpNotLowShow (ao, mp, calledp, &z, &z1, &z2, 0) ;
	  snpNotLowShow (ao, mm, calledm, &z, &z1, &z2, 0) ;
	}

      aceOutf (ao, "\n") ;
    }

   if (1 && oldSnps)
    {
      long int jj, kk ;
      int jOldSnp ;
      char buf[1000], *cp ;
      const char *ccp ;
      Associator missed2old = 0 ;
      const void *vp ;
      DICT *mDict = 0 ;

      for (kk = 0 ; kk < 1 /* keySetMax(goodRuns) */ ; kk++)
	{
	  /*
	    if (snp->keepSample && kk == 0)
	    continue ;
	    if (snp->keepSample && ! keySet(goodRuns, kk))
	    continue ;
	    if (! snp->keepSample && kk > 0)
	    continue ;
	  */
	  mDict = dictHandleCreate (10000, h) ; /* missed snp */
	  missed2old = assReCreate (missed2old) ;
	  for (jj = 1 ; jj <= dictMax (oldSnps) ; jj++)
	    {
	      if (bitt (bb[kk], jj))  /* this snp was reported */
		continue ;
	      	      
	      strcpy (buf, dictName (oldSnps, jj)) ;
	      cp = strchr(buf,':') ;
	      cp = strchr(cp+1,':') ;
	      *cp = 0 ; /* we keep   target:position */
	      dictAdd (mDict, buf, &iOldSnp) ; 	 
	      assMultipleInsert (missed2old, assVoid (iOldSnp), assVoid (jj)) ;
	    }
	  /* multiple insert car si il y a 2 mutations au meme endroit,  missed2old n'en retient que une ? */

	  /* now we scan again the table */
	  bb1 = bitSetCreate (1024, h) ;
	  for (ii = 0, ssm = bigArrp (snps, 0, SSM) ; ii < bigArrayMax (snps) ; ssm++, ii++)
	    {
	      int iBucket = 0 ;
	      Array bucket = 0 ;
	      if (1 && ssm->snp != 1 && ssm->snp != -1) continue ; /* just check the wild type */
	      if (0 && ssm->snp) continue ; /* just check the wild type */
	      memset (Abuf, 0, sizeof(Abuf)) ; memset (Bbuf, 0, sizeof(Bbuf)) ;
	      snpBRS2snpName (snp, ssm, namBuf, typeBuf, Abuf, Bbuf, tagBuf, &delta, 0) ;
	      if (Abuf[11]) messcrash ("Abuff 11") ;
	      cp = strchr(namBuf,':') ;
	      cp = strchr(cp+1,':') ;
	      *cp = 0 ; /* we keep   target:position */
	      if (! dictFind (mDict, namBuf, &jOldSnp))
		continue ; 

	      if (bitt (bb1, jOldSnp)) /* export only once */
		continue ;
	      bitSet (bb1, jOldSnp) ; 

	      bucket = 0 ; iBucket = 0 ; SNPDEBUG=1 ;
	      while (assFindNext (missed2old, assVoid (jOldSnp), &vp, &bucket, &iBucket)) 
		  {
		    int sn = 0 ;
		    jj = assInt (vp) ;
		    ccp = dictName (oldSnps, jj) ;

		    if (iBucket >= 2 && ! bucket)
		      ;
		    /* nnn++ ; printf("nnn=%d\n",nnn) ; */
		    memset (Abuf, 0, sizeof(Abuf)) ; memset (Bbuf, 0, sizeof(Bbuf)) ; memset (tagBuf, 0, sizeof(tagBuf)) ;
		    snpBRS2snpName (snp, ssm, namBuf, typeBuf, Abuf, Bbuf, tagBuf, &delta, ccp) ;

		    dictAdd (snp->eeDict, tagBuf, &sn) ;  sn++ ;
		    if (! snpCountOtherStrand (snp, ssm, sn, ii, delta, ssm->leftShift, typeBuf, Abuf, &cover, &calledp, &calledm, &mp, &mm, &wp, &wm, &oo1p, &oo1m, &oo2p, &oo2m, 1))
		      continue ;
		    if (Abuf[11]) messcrash ("Abuff 11") ;

		    mutant = mp + mm ; wild = wp + wm ;
		    iOldSnp = 0 ;
		    
		    if (1)
		      {
			aceOutf (ao, "%s\t%d\t%s"
				 , dictName (targetDict, ssm->target), ssm->pos - ssm->leftShift
				 , typeBuf
				 ) ;
			dna = array (snp->dnas, ssm->target, Array) ;
			bufS[0] = bufD[0] = 0 ;
			if (dna)
			  {
			    int i, j ;

			    for (j = 0, i = ssm->pos - ssm->leftShift - prefix ; j < prefix && i < 0 ; j++, i++)  
			      bufD[j] = 'N' ;
			    for ( ; j < prefix && i < arrayMax(dna) ; j++, i++)  
			      bufD[j] = arr (dna, i, char) > 0 ? ace_upper(arr (dna, i, char)) : 'N' ;
			    bufD[j] = 0 ;
			    if (delta == 0) i++ ;
			    else if (delta < 0) i -= delta ; 
			    /* if delta > 0) do not move */
			    for (j = 0 ; j < prefix && i < arrayMax(dna) ; j++, i++)  
			      bufS[j] =  arr (dna, i, char) > 0 ? ace_upper(arr (dna, i, char)) : 'N' ;
			    for (; j < prefix ; j++) 
			      bufS[j] = 'N' ;
			    bufS[j] = 0 ;
			    if (delta < 0)
			      {
				cp = Abuf ;
				for (j = 0 ;  ace_lower (bufS[j]) == ace_lower (*cp) ; j++, cp = (cp < Abuf - delta - 1 ? cp + 1 : Abuf))
				  bufS[j] = ace_lower(bufS[j]) ;
			      }
			    if (delta > 0)
			      {
				cp = Bbuf ;
				for (j = 0 ; ace_lower (bufS[j]) == ace_lower (*cp) ; j++, cp = (cp < Bbuf + delta - 1 ? cp + 1 : Bbuf))
				  bufS[j] = ace_lower(bufS[j]) ;
			      }
			  }
			{
			  const char *mmm = "---" ;
			  aceOutf (ao, "\t%s%s%s\t%s%s%s"
				   , bufD, delta <= 0 ? Abuf : mmm + 3 - delta , bufS
				   , bufD, delta >= 0 ? Bbuf : mmm + 3 + delta , bufS
				   ) ;
			}
			aceOutf (ao, "\t%s", dictName (snp->runDict, ssm->run)) ;

			snpNotLowShow (ao, mutant, calledp + calledm, &z, &z1, &z2, 0) ;
			
			aceOutf (ao, "\t%.1f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d"
				 , z
				 , calledp + calledm, mutant, wild
				 , mp, wp, mm, wm
				 , cover - calledp - calledm
				 ) ;
			if (snp->solid && oo1p + oo2p) /* solid */
			  aceOutf (ao, "\t%d:%d"
				   , oo1p, oo2p
				   ) ;
			else
			  aceOutf (ao, "\t%d", oo2p) ;
			if (snp->solid && oo1m + oo2m) /* solid */
			  aceOutf (ao, "\t%d:%d"
				   , oo1m, oo2m
				   ) ;
			else
			  aceOutf (ao, "\t%d", oo2m) ;

			if (1)
			  {
			    char *ss = "Both_strands" ; cz2 = 0 ;
                            if (mp + wp + mm + wm < snp->minCover) ss = "-" ;
			    else if (2 * (mp + wp) < snp->minCover) ss = "Minus_strand_only" ;
			    else if (2 * (mm + wm) < snp->minCover) ss = "Plus_strand_only" ;
			    else
			      { 
				double x, y ;
				chi2 (mp, wp, mm, wm, &cz2) ;
				x = mp/(1.0 + mp + (double)wp) ;
				y = mm/(1.0 + mm + (double)wm) ;
				if (cz2 > 100 && (x > y + .5 || x < y - .5)) ss = "Incompatible_strands" ; /* extremely high  and clear contradiction */ 
			      }
			    aceOutf (ao, "\t%s\t%.1f", ss, cz2) ;
			  } 
			aceOutf (ao, "\t%s"
				 , 100 * (oo1p + oo1m + oo2p + oo2m) > 40 * cover ? "N_rich" : "-"
				 ) ;
			snpNotLowShow (ao, mp, calledp, &z, &z1, &z2, 0) ;
			snpNotLowShow (ao, mm, calledm, &z, &z1, &z2, 0) ;
			aceOutf (ao, "\n") ;
		      }
		  }
	    }
	}
      assDestroy (missed2old) ;
    }
  
  /* export the global histo */
  if (wKeptTranscripts)
    {
      long int wCumul = 0 ;
      
      for (i = 0 ; i < keySetMax (wGlobalCover) && i < 16000 ; i+= wStep)
	wCumul += keySet (wGlobalCover, i) ;
      aceOutf (wao, "\n%s\tKept %d/%d transcripts", snp->run, wKeptTranscripts, wCountTranscripts) ;
      for (i = 0 ; i < keySetMax (wGlobalShiftedCover) && i < 16000 ; i+= wStep)
	aceOutf (wao, "\t%d", keySet (wGlobalShiftedCover, i)) ;

      if (wDebug2)
	{
	  aceOutf (wao, "\n%s\tKept %d/%d raw transcripts", snp->run, wKeptTranscripts, wCountTranscripts) ;
	  for (i = 0 ; i < keySetMax (wGlobalCover) && i < 16000 ; i+= wStep)
	    aceOutf (wao, "\t%d", keySet (wGlobalCover, i)) ;
	}
      aceOutf (wao, "\n") ;
    }

 done:
  ac_free (h) ;

  return nn ;
} /* snpBRS2snpExport */

/*************************************************************************************/

static int snpBRS2snp (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;

  snp->dnas = arrayHandleCreate (8, Array, h) ;
  snpBRSparseFile (snp, FALSE, h) ;
  snpBRS2snpExport (snp) ;

  ac_free (h) ;
  return nn ;
} /* snpBRS2snp */
#endif

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* scan the VariantDB acedb database
 * associate the SNPs and the phenotypes
 */
static int snpPhenotype (SNP *snp)
{
  AC_HANDLE  h1 = 0, h = ac_new_handle () ;
  typedef struct phenoStruct { KEY run, sample, p1, p2 ; float z ; } PHENO ;
  int nn = 0, nnnn = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  AC_DB db = snp->db ;
  AC_ITER iter ;
  AC_TABLE dr = 0, phenoTable = 0, phenoTable2 = 0 ;
  AC_OBJ variant = 0 ;
  Array phenos = 0, myPhenos = 0 ;
  PHENO *pheno, *myPheno ;
  Associator phenoAss = 0 ;
  DICT *phenoDict = 0 ;
  const char *errors = 0 ;
  int i, ir, jr, phenoKey ;
  KEY run, genotype ;
  char *genotypes[] =   { "ww", "wx", "wm", "mx", "mm", 0 } ;
  int genotypeScore[] = { -100,  50 , 100 , 150 , 200 , 0 } ;
  KEY genotypeTags [12] ;

  /* init the genoypeTags */
  {
    char **cpp ;
    for (i = 0, cpp = genotypes ; *cpp ; i++, cpp++)
      genotypeTags[i] = str2tag(*cpp) ;
  }

  /* create an array of phenotype measurements per run */
  phenoTable = ac_bql_table (db, "select r,p1,p2,z from r in class \"run\" where r#is_run, s in r->sample, p1 in s->Phenotype, p2 in p1[1], z in p2[1]", 0, 0, &errors, h) ;

  if (phenoTable && phenoTable->rows)
    {
      float ztot = 0 ;

      phenos = arrayHandleCreate (phenoTable->rows, PHENO, h) ;
      myPhenos = arrayHandleCreate (256, PHENO, h) ;
      phenoDict = dictHandleCreate (256, h) ;
      phenoAss = assHandleCreate (h) ;
      for (ir = jr = 0 ; ir < phenoTable->rows ; ir++)
	{
	  run = ac_table_key (phenoTable, ir, 0, 0) ;
	  if (! run) continue ;

	  assInsert (phenoAss, assVoid(run), assVoid(jr)) ;
	  pheno = arrayp (phenos, jr++, PHENO) ;
	  pheno->run = run ;
	  pheno->p1 =  ac_table_key (phenoTable, ir, 1, 0) ;
	  pheno->p2 =  ac_table_key (phenoTable, ir, 2, 0) ;
	  pheno->z = ac_table_float (phenoTable, ir, 3, -999999999) ;
	  if (pheno->z < -999999998) 
	    jr-- ;
	  else
	    ztot += pheno->z ;
	}
      arrayMax (phenos) = jr ; /* because we may have cancelled a line */

      /* reexport again but without the run so that the phenoDict is correctly sorted */
	phenoTable2 = ac_bql_table (db, "select p1,p2,z from r in class \"run\" where r#is_run, s in r->sample, p1 in s->Phenotype, p2 in p1[1], z in p2[1]  where z", 0, 0, &errors, h) ; 
      for (ir = 0 ; ir < phenoTable2->rows ; ir++)
	{
	  dictAdd (phenoDict
		   , messprintf ("%d_%d"
				 , ac_table_key (phenoTable2, ir, 0, 0)
				 , ac_table_key (phenoTable2, ir, 1, 0)
				 )
		   , &phenoKey) ;
	}      
      /* center the pheno */
      {
	Array aan = arrayCreate (dictMax(phenoDict), int) ;
	Array aaz = arrayCreate (dictMax(phenoDict), float) ;
	
	for (ir = 0, pheno = arrp (phenos, ir, PHENO) ; jr && ir < jr ; pheno++, ir++)
	  {
	    if (! dictFind (phenoDict
			    , messprintf ("%d_%d", pheno->p1, pheno->p2)
			    , &phenoKey))
	      messcrash ("Cannot locate %s-> %s %s in phenoDict", name(pheno->run), name(pheno->p1), name(pheno->p2));
	    array (aaz, phenoKey, float) += pheno->z ;
	    array (aan, phenoKey, int) ++ ;
	  }
	for (ir = 0, pheno = arrp (phenos, ir, PHENO) ; jr && ir < jr ; pheno++, ir++)
	  {
	    dictFind (phenoDict
		      , messprintf ("%d_%d", pheno->p1, pheno->p2)
		      , &phenoKey) ;
	    pheno->z -= array (aaz, phenoKey, float)/array (aan, phenoKey, int) ;
	  }
	arrayDestroy (aan) ;
	arrayDestroy (aaz) ;
      }
    }

  iter = ac_query_iter (db, TRUE, "find variant FQ ", 0, h) ;

  while (ac_free (variant), ac_free (h1), vtxtClear (txt), variant = ac_iter_obj (iter))
    {
      nnnn++ ;
      if (0 && nnnn > 20)
	break ;
      h1 = ac_new_handle () ; 
      vtxtPrintf (txt, "Variant %s\n-D Phenotype\n", ac_protect (ac_name (variant), h1)) ;
	
      dr = ac_tag_table (variant, "FQ", h1) ;
      if (! dr || ! dr->rows) continue ;

      /* compute the convolution between the phenotype and the mutant/wild-type counts */
      myPhenos = arrayReCreate (myPhenos, 256, PHENO) ;      
      for (ir = 0 ; phenoDict && ir < dr->rows ; ir++)
	{
	  void *vp ;
	  KEY *kp ;
	  int score ;
	  
	  genotype = ac_table_key (dr, ir, 0, 0) ;
	  if (genotype)
	    {
	      for (i = score = 0, kp = genotypeTags ; *kp ; kp++, i++)
		if (*kp == genotype)
		  { score = genotypeScore [i] ; break ; }
	      if (score)
		{
		  run = ac_table_key (dr, ir, 1, 0) ;	
	  
		  if (run && assFind (phenoAss, assVoid(run), &vp))
		    {
		      for (jr = assInt (vp), pheno = arrp (phenos, jr, PHENO) ; jr < arrayMax (phenos) && pheno->run == run ; pheno++, jr++)
			{
			  if (dictFind (phenoDict, messprintf ("%d_%d", pheno->p1, pheno->p2), &phenoKey))
			    {
			      myPheno = arrayp (myPhenos, phenoKey, PHENO) ;
			      myPheno->p1 = pheno->p1 ;
			      myPheno->p2 = pheno->p2 ;
			      myPheno->z += score * pheno->z ;
			    }
			}
		    }
		}
	    }
	}
      
      /* export the phenotypes */
      if (1 && arrayMax (myPhenos))
	{
	  for (jr = 0, myPheno = arrp (myPhenos, 0, PHENO); jr < arrayMax (myPhenos); myPheno++, jr++)
	    if (myPheno->p1)
	      vtxtPrintf(txt, "Phenotype \"%s\" \"%s\" %f\n", ac_key_name(myPheno->p1),  ac_key_name(myPheno->p2), myPheno->z) ;
	  
	  /* edit this variant */
	  vtxtPrintf (txt, "\n") ;
	  ac_parse (db, vtxtPtr (txt), &errors, 0, h1) ; 
	  if (errors && *errors)
	    messerror (errors) ;
	}
    }
  ac_free (h) ;
  return nn ;
} /* snpPhenotype */
 
/*************************************************************************************/
/*************************************************************************************/
/* extend the snpnet using the genome
 */
/*************************************************************************************/
/* returns sliding and edits the buffers */
static BOOL snpIsSliding (char *typ, char *buf1, char *buf2, BOOL isRNA, BOOL isDown, int *vcfPosp, char *vcfBuf)
{
  int i, j, ddx = strlen(typ) - 3 ;
  BOOL isSliding = FALSE ;
  char *cp, *cq, *bb[3] ;
  
  bb[0] = buf1 ; bb[1] = buf2 ; bb[2] = typ + 3 ;
  for (i = 0 ; i < 3 ; i++)
    for (cp = bb[i] ; *cp ; cp++)
      {
	if (isRNA)
	  *cp = ace_upper (rnaDecodeChar[(int)dnaEncodeChar[(int) *cp]]) ;
	else
	  *cp = ace_upper (dnaDecodeChar[(int)dnaEncodeChar[(int) *cp]]) ;	
      }
  if (vcfBuf)
    vcfBuf[0] = 0 ;

  if (! isDown)
    { /* complement buf1 and buf2 */
      char cc1, cc2 ;
      for (i = 0 ; i < 2 ; i++)
	{
	  int j, ln = strlen(bb[i]) ;
	  cp = bb[i] ;
	  for (j = 0 ; j < (ln+1)/2 ; j++)
	    {
	      cc1 = cp[j] ; cc2 = cp[ln - 1 - j] ;
	      if (isRNA)
		{
		  cp[j] = ace_lower(rnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)cc2])]) ;
		  cp[ln - 1 - j] = ace_lower(rnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)cc1])]) ;
		}
	      else
		{
		  cp[j] = ace_upper(dnaDecodeChar[(int)dnaEncodeChar[(int)cc2]]) ;
		  cp[ln - 1 - j] = ace_upper(dnaDecodeChar[(int)dnaEncodeChar[(int)cc1]]) ;
		}
	    }
	}
    }
  
  if (typ[1] == '>' || typ[1] == '2')
    {
      typ[1] = '>' ;
      if (isDown)
	{
	  if (isRNA)
	    {
	      typ[0] = ace_lower(rnaDecodeChar[(int)dnaEncodeChar[(int)typ[0]]]) ;
	      typ[2] = ace_lower(rnaDecodeChar[(int)dnaEncodeChar[(int)typ[2]]]) ;
	    }
	  else
	    {
	      typ[0] = ace_upper(dnaDecodeChar[(int)dnaEncodeChar[(int)typ[0]]]) ;
	      typ[2] = ace_upper(dnaDecodeChar[(int)dnaEncodeChar[(int)typ[2]]]) ;
	    }
	}
      else
	{
	  if (isRNA)
	    {
	      typ[0] = ace_lower(rnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)typ[0]])]) ;
	      typ[2] = ace_lower(rnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)typ[2]])]) ;
	    }
	  else
	    {
	      typ[0] = ace_upper(dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)typ[0]])]) ;
	      typ[2] = ace_upper(dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)typ[2]])]) ;
	    }
	}

      if (vcfBuf)
	{
	  vcfBuf[0] = ace_upper (typ[0]) ;
	  vcfBuf[1] = ' ' ;
	  vcfBuf[2] = ace_upper (typ[2]) ; 
	  vcfBuf[3] = 0 ;
	}
    }
  
  if (!strncmp (typ, "Del", 3))
    { 
      int vcfPos = vcfPosp ? *vcfPosp - 1 : 0 ;
      
      cp = buf1 ; cq = buf2 ;
      while (ace_upper (*cp) == ace_upper(*cq))
	{ cp++ ; cq++ ; } 
      cp-- ; cq-- ;
      
      j = 0 ;
      while (cp > buf1 && ace_upper (cp[0]) == ace_upper (cp[ddx]))
	{ cp-- ; vcfPos-- ; j++ ; }
      if (j >= ddx)
	isSliding = TRUE ; 
      
      if (vcfBuf)
	{
	  *vcfPosp = vcfPos ;
	  for (i = 0 ; i < ddx + 1 && i < 60 ; i++) 
	    vcfBuf[i] = ace_upper (cp[i]) ;
	  vcfBuf[i++] = ' ' ;
	  vcfBuf[i++] = ace_upper (cp[0]) ;
	  vcfBuf[i++] = 0 ;
	}
    } 
  else if (!strncmp (typ, "Ins", 3))
    {
      int vcfPos = vcfPosp ? *vcfPosp - 1 : 0 ;
      
      cp = buf1 ; cq = buf2 ;
      while (ace_upper (*cp) == ace_upper(*cq))
	{ cp++ ; cq++ ; } 
      cp-- ; cq-- ;
      
      j = 0 ;
      while (cq > buf2 && ace_upper (cq[0]) == ace_upper (cq[ddx]))
	{ cq-- ; vcfPos-- ;  j++ ; }
      if (j >= ddx)
	isSliding = TRUE ; 
       
      if (vcfBuf)
	{
	  *vcfPosp = vcfPos ;
	  vcfBuf[0] = vcfBuf[2] = ace_upper (*cq) ;
	  vcfBuf[1] = ' ' ;
	  for (i = 1 ; i < ddx + 1 && i < 60 ; i++) 
	    vcfBuf[2+i] = ace_upper (cq[i]) ;
	  vcfBuf[2+i] = 0 ;
	}
    }
	      
  /* fix the upper lower case */
  cp = buf1 ; cq = buf2 ;
  while (*cp && *cq && ace_upper (*cp) == ace_upper(*cq))
    { 
      *cp = ace_upper (*cp) ; cp++ ;
      *cq = ace_upper (*cq) ; cq++ ;
    }
  if (!strncmp (typ, "Del", 3))
    {
      for (i = 0 ; *cp && i < ddx ; cp++, i++)
	*cp = ace_lower (*cp) ;
    }
  else if (!strncmp (typ, "Ins", 3))
    {
      for (i = 0 ; *cq && i < ddx ; cq++, i++)
       *cq = ace_lower (*cq) ;
    }
  else if (*cp && *cq)
    {
      *cp = ace_lower (*cp) ; cp++ ;
      *cq = ace_lower (*cq) ; cq++ ;
    }
  while (*cp)
    {
      *cp = ace_upper (*cp) ; cp++ ;
    }
  while (*cq)
    {
      *cq = ace_upper (*cq) ; cq++ ;
    }
  
  return isSliding ;
} /* snpIsSliding */

/*************************************************************************************/
#ifdef JUNK
static int snpExtendGenomicSnpNet  (SNP *snp)
{
  AC_HANDLE  h1 = 0, h = ac_new_handle () ;
  int n = 0, nn = 0, nnnn = 0 ; ;
  vTXT txt = vtxtHandleCreate (h) ;
  AC_DB db = snp->db ;
  AC_ITER iter ;
  AC_TABLE map, intMap ;
  AC_OBJ variant = 0, seq = 0 ;
  KEY chrom, oldChrom = 0 ;
  const char  *dna1, *dna = 0 ;
  char buf1[1024], buf2[1024], typBuf[2048] ; /* typBuf must accpomodate long multisub */
  const char *typ ;
  int a0, a11, a22, dx1, dx2 ;
  const char *errors = 0 ;
  BOOL isDown = TRUE ;

  if (1)
    iter = ac_query_iter (db, TRUE, "find variant", 0, h) ;  
  else
    iter = ac_query_iter (db, TRUE, "find variant IntMap && (IS \"NC_045512:26550:Del_21:gTTGAAGAGCTTAAAAAGCTCC:g\")", 0, h) ;
  while (ac_free (variant), ac_free (h1), vtxtClear (txt), variant = ac_iter_obj (iter))
    {
      intMap = ac_tag_table (variant, "IntMap", h1) ;      
      BOOL ok = FALSE ;
      
      memset (typBuf, 0, sizeof (typBuf)) ;
      nnnn++ ;
      if (0) fprintf (stderr, "---- %d\n", nnnn) ;
      if (ac_has_tag (variant, "Found_in_mRNA"))
	{
	  map = ac_tag_table (variant, "mRNA", h1) ;
	  if (intMap && ac_table_int (intMap, 0, 2, 1) == -1) 
	    isDown = FALSE ;
	}
      else
	map = intMap ;
      if (! map || map->cols < 3) 
	continue ;
      chrom = ac_table_key (map, 0, 0, 0) ;
      if (chrom != oldChrom)
	{
	  if (! chrom) continue ; 
	  if (ac_has_tag (variant, "Found_in_mRNA"))
	    seq = ac_tag_obj (variant, "mRNA", h1) ;
	  else
	    seq = ac_get_obj (db, "Sequence", ac_table_printable (map, 0, 0, 0), h1) ;
	  if (! seq) continue ;
	  ac_free (dna) ;
	  dna = ac_obj_dna (seq, h) ; 
	  n =  dna ? strlen (dna) : 0 ;
	  oldChrom = chrom ;
	}
      if (! dna)
	continue ;
      n =  strlen (dna) ;
      a0 = ac_table_int (map, 0, 1, 0) ;
      a11 = a0 - 20 ; if (a11 < 1) a11 = 1 ;
      if (a0 < 1)
	continue ;
      a22 = a0 + 1020 ;
      if (n < a22) a22 = n ;
      if (a22 <= a11)
	continue ;
      dna1 = dna + a11 - 1 ;
      strncpy (buf1, dna1, a0 - a11) ; dx1 = a0 - a11 ;
      strncpy (buf2, dna1, a0 - a11) ; dx2 = a0 - a11 ;
      
      ok = FALSE ;
      typ = ac_tag_printable (variant, "Transition", 0) ;
      if (! typ) typ = ac_tag_printable (variant, "Transversion", 0) ;
      if (typ)
	{
	  ok = TRUE ;
	  typBuf[0] = buf1[a0 - a11] = ace_upper (typ[0]) ; dx1++ ;
	  typBuf[2] = buf2[a0 - a11] = ace_upper (typ[2]) ; dx2++ ;
	  typBuf[1] = '2' ;
	  typBuf[3] = 0 ;
	  if (a22 > a11 + dx1 + 20) a22 = a11 + dx1 + 1 ;
	  memcpy (buf1 + dx1, dna1 + dx1, a22 - a11 - dx1) ; 
	  memcpy (buf2 + dx2, dna1 + dx1, a22 - a11 - dx2) ; 
	  dx1 += a22 - a11 - dx1 ; dx2 += a22 - a11 - dx2 ;
	}
      if (! ok)
	{
	  typ = ac_tag_printable (variant, "Multi_substitution", 0) ;
	  if (typ)
	    {
	      AC_TABLE t1 = ac_tag_table (variant, "Multi_substitution", h1) ;

	      ok = TRUE ;
	      if (t1 && t1->cols >= 3 && ac_table_printable (t1,0,1,0) && ac_table_printable (t1,0,2,0))
		{
		  const  char *ccp ;
		  int tbn = 0 ;
		  ccp = ac_table_printable (t1,0,1,0)  ;
		  n = strlen (ccp) ;
		  if (n + dx1 > 500) n = 500 - dx1 ;
		  if (n > 0)
		    { memcpy (buf1 + dx1, ccp, n) ; dx1 += n ; memcpy (typBuf, ccp, n) ; tbn = n ; }
		  typBuf[tbn++] = '2' ;
		  ccp = ac_table_printable (t1,0,2,0)  ;
		  n = strlen (ccp) ;
		  if (n + dx2 > 1000) n = 100 - dx2 ;
		  if (n > 0)
		    { memcpy (buf2 + dx2, ccp, n) ;  dx2 += n ; memcpy (typBuf+tbn, ccp, n) ; tbn += n ; } 
		  typBuf[tbn++] = 0 ;		 

		  n = 20 ;
		  if (a22 < a11 + dx1 + n) n = a22 - a11 - dx1 ;
		  if (n>0)
		    {
		      memcpy (buf1 + dx1, dna1 + dx1, n) ; 
		      memcpy (buf2 + dx2, dna1 + dx1, n) ; 
		      dx1 += n ; dx2 += n ;
		    }
		}
	    }
	}
      if (! ok)
	{
	  typ = ac_tag_printable (variant, "Single_deletion", 0) ;
	  if (typ)
	    {
	      ok = TRUE ;
	      int i, dx = 1 ;
	      
	      memcpy (buf1 + dx1, dna1 + dx1, 1) ;                  /* because the variant name reports a la VCF the base before the deletion */
	      memcpy (buf2 + dx2, dna1 + dx1, 1) ; dx1++ ; dx2++ ;  /* because the variant name reports a la VCF the base before the deletion */

	      strcpy (typBuf, "Del") ;
	      for (i = 0 ; i < dx ; i++)
		typBuf[3+i] = ace_upper (typ[i+3]) ;
	      typBuf[3+i] = 0 ;
	      for (i = 0 ; i < dx ; i++, dx1++)
		buf1[dx1] = ace_upper (dna1[dx1]) ;
	      n = 20 ;
	      if (a22 < a11 + dx1 + n) n = a22 - a11 - dx1 ;
	      if (n>0)
		{
		  memcpy (buf1 + dx1, dna1 + dx1, n) ;
		  memcpy (buf2 + dx2, dna1 + dx1, n) ; 
		  dx1 += n ; dx2 += n ;
		}
	    }
	}
      if (! ok)
	{
	  int dx = ac_tag_int (variant, "Multi_deletion", 0) ;
	  if (dx && dx < 1000)
	    {
	      ok = TRUE ;
	      int i ;
	      
	      memcpy (buf1 + dx1, dna1 + dx1, 1) ;                  /* because the variant name reports a la VCF the base before the deletion */
	      memcpy (buf2 + dx2, dna1 + dx1, 1) ; dx1++ ; dx2++ ;  /* because the variant name reports a la VCF the base before the deletion */

	      strcpy (typBuf, "Del") ;
	      for (i = 0 ; i < dx ; i++, dx1++)
		{
		  buf1[dx1] = ace_upper (dna1[dx1]) ;
		  typBuf[3+i] = buf1[dx1] ;
		}
	      typBuf[3+i] = 0 ;
	      n = 20 ;
	      if (a22 < a11 + dx1 + n) n = a22 - a11 - dx1 ;
	      if (n>0)
		{
		  memcpy (buf1 + dx1, dna1 + dx1, n) ;
		  memcpy (buf2 + dx2, dna1 + dx1, n) ; 
		  dx1 += n ; dx2 += n ;
		}
	    }
	}
      if (! ok)
	{
	  typ = ac_tag_printable (variant, "Single_insertion", 0) ;
	  if (typ)
	    {
	      int i, dx = 1 ;
	      ok = TRUE ;
	      
	      memcpy (buf1 + dx1, dna1 + dx1, 1) ;                  /* because the variant name reports a la VCF the base before the deletion */
	      memcpy (buf2 + dx2, dna1 + dx1, 1) ; dx1++ ; dx2++ ;  /* because the variant name reports a la VCF the base before the deletion */

	      strcpy (typBuf, "Ins") ;
	      for (i = 0 ; i < dx ; i++)
		typBuf[3+i] = ace_upper (typ[i+3]) ;
	      typBuf[3+i] = 0 ;
	      for (i = 0 ; i < dx ; i++, dx2++)
		buf2[dx2] = ace_upper (typ[i + 3]) ;
	      n = 20 ;
	      if (a22 < a11 + dx1 + n) n = a22 - a11 - dx1 ;
	      if (n>0)
		{
		  memcpy (buf1 + dx1, dna1 + dx1, n) ;
		  memcpy (buf2 + dx2, dna1 + dx1, n) ; 
		  dx1 += n ; dx2 += n ;
		}
	    }
	}
      if (! ok)
	{
	  int dx = ac_tag_int (variant, "Multi_insertion", 0) ;
	  if (dx)
	    {
	      AC_TABLE t1 = ac_tag_table (variant, "Multi_insertion", h1) ;
	      if (t1 && t1->cols >= 2)
		{
		  const char *ccp = ac_table_printable (t1,0,1,0) ;

		  memcpy (buf1 + dx1, dna1 + dx1, 1) ;                  /* because the variant name reports a la VCF the base before the deletion */
		  memcpy (buf2 + dx2, dna1 + dx1, 1) ; dx1++ ; dx2++ ;  /* because the variant name reports a la VCF the base before the deletion */

		  if (ccp)
		    {
		      if (dx <= strlen (ccp))
			{
			  dx = strlen (ccp) ; 
			  if (dx > 1000 + dx2)
			    dx = 1000 - dx2 ;
			  if (dx > 0)
			    {
			      int i ;
			      
			      ok = TRUE ;
			      strcpy (typBuf, "Ins") ;
			      for (i = 0 ; i < dx ; i++)
				typBuf[3+i] = ace_upper (ccp[i]) ;
			      typBuf[3+i] = 0 ;
			      for (i = 0 ; i < dx ; i++, dx2++)
				buf2[dx2] = ace_upper (ccp[i]) ;
			      n = 20 ;
			      if (a22 < a11 + dx1 + n) n = a22 - a11 - dx1 ;
			      if (n>0)
				{
				  memcpy (buf1 + dx1, dna1 + dx1, n) ;
				  memcpy (buf2 + dx2, dna1 + dx1, n) ; 
				  dx1 += n ; dx2 += n ;
				}
			    } 
			}
		    }
		}
	    }
	}
      if (! ok)
	continue ;
      buf1[dx1] = buf2[dx2] = 0 ;

      /* edit this variant */	  
      /* ATTENTION 2018_01_12 if we start from_mrna and we are on negative strand this is completely FAUX
       * we need to complement everything, but how do we deal with introns ?
       * also the vcfPos should be a genomic position, not a0 which is a mrna position
       * we should start from a0 on the mrna, slide it on the mrna, then remap it to the genome
       *
       * if we start from the genome, we slide it and then remap it to the mrna
       * but it may have jumped an intron
       * in both cases, on the negative strand, we must slide in opposite directions on mrna and on genome
       */ 
      if (1)
	{ 
	  int vcfPos = intMap ? ac_table_int (intMap, 0, 1, 0) : 0 ;
	  char vcfBuf[128] ;
	  BOOL sliding = FALSE ;

	  vtxtPrintf (txt, "Variant %s\n", ac_protect (ac_name (variant), h1)) ;
	  vtxtPrintf(txt, "Reference_genomic_sequence  %s\n", buf1) ;
	  vtxtPrintf(txt, "Observed__genomic_sequence   %s\n", buf2) ;
	  vtxtPrintf(txt, "Typ %s\n", typBuf) ;
	  
	  sliding =  snpIsSliding (typBuf, buf1, buf2, FALSE, isDown, &vcfPos, vcfBuf) ;
	  vtxtPrintf(txt, "VCF %d %s\n", vcfPos, vcfBuf) ;
	  if (sliding)
	    {
	      vtxtPrintf(txt, "In_repeat\n") ;

	      if (typBuf[4] && !strncmp (typBuf, "Del", 3))
		vtxtPrintf(txt, "Dim%c\nIn_repeat\nDiminution\n", typBuf[3]) ;
	      if (typBuf[4] && !strncmp (typBuf, "Ins", 3))
		vtxtPrintf(txt, "Dup%c\nIn_repeat\nDuplication\n", typBuf[3]) ;
	    }
	  else
	    vtxtPrintf(txt, "-D In_repeat\n") ;
	  
	  vtxtPrintf(txt, "\n") ;
	  ac_parse (db, vtxtPtr (txt), &errors, 0, h1) ;  
	}
      vtxtClear (txt) ;
      nn++ ;
    }
  ac_free (h1) ;
  ac_free (h) ;
  return nn ;
} /* snpExtendGenomicSnpNet */
#endif

static int snpCompatibleStrands (int cover, int a1, int a2, int b1, int b2, double *zp, double *z1p, double *z2p, int *chiFlagp, int limit)
{
  int z3 = 0 ;
  
  /*
   * z1 = snpNotLow (a1, a1+a2, 0, 0, 0, 4) ;
   * z2 = snpNotLow (b1, b1+b2, 0, 0, 0, 4) ;
   */
  z3 = snpNotLow (a1+b1, cover, zp, z1p, z2p, limit) ;

  if (chiFlagp)
    *chiFlagp = chi2 (a1,a2,b1,b2, 0) ;
  
  return z3 ;    
} /* snpCompatibleStrands */

/*************************************************************************************/
/*************************************************************************************/
/* parse a (collection of) .snp tables, report per 6-mer the frequency of oo errors */
static int snpOoFrequency (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, i, nb = 0 ;
  int cc, oligo = 0 ;  /* 32 bits = 8*4 so we can store an 8-mer */
  ACEIN ai = snp->ai ;
  ACEOUT ao = snp->ao ;
  const char *ccp ;
  BOOL oldDown = TRUE, isDown = TRUE ;
  int B,R,S,T ;
  int LN = snp->ooFrequency ;
  int MXX = (1 << (2 * LN)) ;
  int pos = 0, oldPos = 0, nw[LN], noo[LN] ;
  int mask, nw3D[MXX], no3D[MXX], nw3U[MXX], no3U[MXX] ;
  
  mask = (1 << (2*LN)) - 1 ;

  memset (nw, 0, sizeof (nw)) ;
  memset (noo, 0, sizeof (noo)) ;


  memset (nw3D, 0, sizeof (nw3D)) ;
  memset (no3D, 0, sizeof (no3D)) ;
  memset (nw3U, 0, sizeof (nw3U)) ;
  memset (no3U, 0, sizeof (no3U)) ;

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (! ccp || !strncmp (ccp, "ERR", 3))
	continue ;
      if (*ccp != '~') /* new target, clear the buffers */
	goto ooCleanUp ;
      if (! aceInInt (ai, &pos))
	goto ooCleanUp ;
      ccp = aceInWord (ai) ;  /* the base or the snp */
      if (! ccp)
	goto ooCleanUp ;
      if (ccp[1]) /* a snp */
	continue ;
      /* else ccp is the wild type base */
      if (nb && pos != oldPos + 1)  /* is it truly the next base */
	goto ooCleanUp ;
      switch ((int)*ccp)
	{
	case 'a': cc = 0x00 ; break ;
	case 't': cc = 0x01 ; break ;
	case 'g': cc = 0x02 ; break ;
	case 'c': cc = 0x03 ; break ;
	default: goto ooCleanUp ;
	}
      
      ccp = aceInWord (ai) ;  /* strand */
      if (! ccp)
	goto ooCleanUp ;
      switch ((int)*ccp)
	{
	case '+': isDown = TRUE ; break ;
	case '-': isDown = FALSE ; break ;
	default: goto ooCleanUp ;
	}
      if (nb && isDown != oldDown)
	 goto ooCleanUp ;
      oldDown = isDown ; oldPos = pos ;
      nb++ ; oligo <<= 2 ; oligo |= cc ; oligo &= mask ;
      ccp = aceInWord (ai) ;  /* run, drop it */
      if (! aceInInt (ai, &B) || ! aceInInt (ai, &R) || ! aceInInt (ai, &S) || ! aceInInt (ai, &T))
	goto ooCleanUp ;
      for (i = 0 ; i < LN - 1 ; i++)
	{ nw[i] = nw[i+1] ; noo[i] = noo[i+1] ; }
      nw[LN - 1] = B ; noo[LN - 1] = T ;
      if (nb >= LN) 
	{
	  if (isDown)
	    {
	      nw3D[oligo] += nw[LN/2 - 1] ;   /* i have read this base as the third base of the hexamer */
	      no3D[oligo] += noo[LN/2 - 1] ;
	    }
	  else
	    {
	      nw3U[oligo] += nw[LN/2 - 1] ;   /* i have read this base as the third base of the hexamer */
	      no3U[oligo] += noo[LN/2 - 1] ;
	    }
	}
      continue ;

    ooCleanUp:
      oldDown = isDown ;
      oldPos = pos ;
      nb = 0 ; oligo = 0 ;
    }

  /* report */
  aceOutDate (ao, "###", snp->project) ;
  aceOutf (ao, "# Frequency of incoherent transitions between base %d and %d a %d-mer in Solid data\n", LN/2, LN/2+1,LN) ;
  aceOutf (ao, "Oligo\tThe central ransition is coherent read on strand +\tThe central transition is incoherent read on strand +\tFrequency of incoherent transition read on strand +") ;
  aceOutf (ao, "\tOligo\tThe central transition is coherent read on strand -\tThe central transition is incoherent read on strand -\tFrequency of incoherent transition read on strand -") ;
  aceOutf (ao, "\tOligo\tThe central transition is coherent read on any strand\tThe central transition is incoherent read on any strand\tFrequencyof incoherent transition read on any strand") ;
  for (oligo = 0 ; oligo < MXX ; oligo++)
    {
      char oBuf[LN+1], *cp ;
      float z ;
      int oligoR = 0, ccR = 0 ;

      for (i = 0, cp = oBuf ; i < LN ; cp++, i++)
	{
	  cc = (oligo >> (2*LN - 2 - 2 * i)) & 0x3 ;
	  switch (cc)
	    {
	    case 0x00: *cp = 'a' ; ccR = 0x01 ; break ;
	    case 0x01: *cp = 't' ; ccR = 0x00 ; break ;
	    case 0x02: *cp = 'g' ; ccR = 0x03 ; break ;
	    case 0x03: *cp = 'c' ; ccR = 0x02 ; break ;
	    }
	  oligoR |= (ccR << (2*i)) ;
	}
      *cp = 0 ;
      z = nw3D[oligo] + no3D[oligo] + nw3U[oligoR] + no3U[oligoR] ; if (z < 1) continue ;
      z = nw3D[oligo] + no3D[oligo] ; if (z < 1) z = 1 ;
      aceOutf (ao, "\n%s\t%d\t%d\t%.1f",
	       oBuf, nw3D[oligo], no3D[oligo], 100.0 * no3D[oligo]/z 
	       ) ;
      z = nw3U[oligoR] + no3U[oligoR] ; if (z < 1) z = 1 ;
      aceOutf (ao, "\t%s\t%d\t%d\t%.1f",
	       oBuf, nw3U[oligoR], no3U[oligoR], 100.0 * no3U[oligoR]/z 
	       ) ;
      z = nw3D[oligo] + no3D[oligo] + nw3U[oligoR] + no3U[oligoR] ; if (z < 1) z = 1 ;
      aceOutf (ao, "\t%s\t%d\t%d\t%.1f",
	       oBuf, nw3D[oligo] + nw3U[oligoR], no3D[oligo] + no3U[oligoR], 100.0 * (no3D[oligo] + no3U[oligoR])/z 
	       ) ;
    }
  aceOutf (ao, "\n") ;

  ac_free (h) ;
  return nn ;
} /* snpOoFrequency */

/*************************************************************************************/
/*************************************************************************************/
/* for the snps in -snp_list and the genome in -fasta, export n bases around the modified sequence */
typedef struct sneStruct { int classe, chrom1, chrom2, a1, a2, type, typeM, typeW ; BOOL isDown1, isDown2 ; } SNE ;

/*************************************************************************************/

static int sneOrder (const void *va, const void *vb)
{
  const SNE *a = (const SNE *)va, *b = (const SNE *)vb ;
  int n ;

  n = a->classe - b->classe ;  if (n) return n ;
  n = a->chrom1 - b->chrom1 ; if (n) return n ;
  n = a->chrom2 - b->chrom2 ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  n = a->isDown1 - b->isDown1 ; if (n) return n ;
  n = a->isDown2 - b->isDown2 ; if (n) return n ;
  n = a->type - b->type ; if (n) return n ;

  return 0 ;
} /* sneOrder */

/*************************************************************************************/

static int snpParseTranslocations (SNP *snp, Array aa, const char **cls) 
{
  AC_HANDLE h = 0 ;
  ACEIN ai = 0 ;
  int nn = arrayMax (aa) ;
  int i, classe = 0, chrom1, chrom2, a1, a2 ;
  BOOL isDown1 = TRUE, isDown2 = TRUE ;
  const char *ccp, **clp ;
  SNE *up ; ;

  if (! snp->newIntronsFileName)
    return 0 ;
  h = ac_new_handle () ;
  ai = aceInCreate (snp->newIntronsFileName, FALSE, h) ;
  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (! ccp)
	continue ;
      classe = 0 ;
      for (clp = cls, i = 0 ; !classe && *clp ; i++, clp++)
	if (! strcmp (ccp, *clp))
	  classe = i ;
      if (! classe) 
	continue ;

      /* first chromosome */
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp)
	continue ;
      if (!dictFind (snp->targetDict, ccp, &chrom1))
	continue ;
      
       /* first orientation */
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp)
	continue ;
      else if (! strcmp (ccp, "Forward"))
	isDown1 = TRUE ;
      else if (! strcmp (ccp, "Reverse"))
	isDown1 = FALSE ;
      else
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &a1))
	continue ;
      
      
      /* second chromosome */
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp)
	continue ;
      if (!dictFind (snp->targetDict, ccp, &chrom2))
	continue ;
      
       /* second orientation */
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp)
	continue ;
      else if (! strcmp (ccp, "Forward"))
	isDown2 = TRUE ;
      else if (! strcmp (ccp, "Reverse"))
	isDown2 = FALSE ;
      else
	continue ;
      aceInStep (ai, '\t') ;
      if (! aceInInt (ai, &a2))
	continue ;

      if (isDown1)
	{
	  up = arrayp (aa, nn++, SNE) ;
	  up->chrom1 = chrom1 ;
	  up->chrom2 = chrom2 ;
	  up->a1 = a1 - 1 ;
	  up->a2 = a2 + 1 ;
	  up->type = 0 ;
	  up->classe = classe ; 
	  up->isDown1 = isDown1 ;
	  up->isDown2 = isDown2 ;
	}
      else /* exchange */
	{
	  up = arrayp (aa, nn++, SNE) ;
	  up->chrom1 = chrom2 ;
	  up->chrom2 = chrom1 ;
	  up->a1 = a2 - 1 ;
	  up->a2 = a1 + 1 ;
	  up->type = 0 ;
	  up->classe = classe ; 
	  up->isDown1 = isDown2 ;
	  up->isDown2 = isDown1 ;
	}   
      if (classe == 9) /* try the reciprocal translocation */
	{
	  SNE *vp = arrayp (aa, nn++, SNE) ;
	  vp->chrom1 = up->chrom2 ;
	  vp->chrom2 = up->chrom1 ;
	  vp->a1 = up->a2 - 2 ;
	  vp->a2 = up->a1 + 2 ;
	  vp->type = 0 ;
	  vp->classe = classe ; 
	  vp->isDown1 = up->isDown2 ;
	  vp->isDown2 = up->isDown1 ;
	}
    }
  
  ac_free (h) ;	
  return nn ;
} /* snpParseTranslocations */

/*************************************************************************************/

static int snpExportEditedSequence (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  int dx, delta = snp->editSequence ;
  int a1, a2, b1, ii, nn = 0, nok = 0, target, x, type ;
  int snpId = 0, oldChrom1 = 0, oldA1 = 0 ;
  const char *ccp ;
  char bufW[100000], bufM[100000], *cp, *cq, *dna1, *dna2 ;
  DICT *typeDict = dictHandleCreate (100, h) ;
  Array aa = arrayHandleCreate (20000, SNE, h) ;
  SNE *up ;
  ACEOUT ao, aoSNP = aceOutCreate (snp->outFileName, ".snp.fasta", snp->gzo, h) ;
  ACEOUT aoBRK = aceOutCreate (snp->outFileName, ".brk.fasta", snp->gzo, h) ;

  const char *cls[] = {"SNP", "MICROINSERTION", "INSERTION","MICRODELETION", "DELETION", "INTRON", "MICRODUPLICATION", "DUPLICATION", "INVERSION", "TRANSLOCATION", 0} ;
  const char *cls2[] = {"SNV", "Ins", "INS", "Del", "DEL", "ITR", "Dup", "DUP", "INV", "BRK", 0} ;

  /* collec the SNPs and create a sorted array */

  memset (bufM, 0, sizeof(bufM)) ;
  memset (bufW, 0, sizeof(bufW)) ;
  snpParseSnpList (snp) ;  /* fills dict: snp->oldSnps */
  for (ii = 1 ; snp->oldSnps && ii <= dictMax (snp->oldSnps) ; ii++)
    {
      target = type = 0 ;
      ccp = dictName (snp->oldSnps, ii) ;
      strncpy (bufM, ccp, 999) ;
      for (cp = bufM ; *cp && *cp != ':' ; cp++) ;
      if (*cp == ':') 
	{
	  *cp = 0 ; 
	  if (!dictFind (snp->targetDict, bufM, &target))
	    continue ;
	}
      for (x = 0, cp++ ; *cp && *cp >= '0' && *cp <= '9' ; cp++)
	x = 10 * x + (*cp - '0') ;
      if (*cp == ':') cp++ ; /* optional ':' */
      if (*cp)
	dictAdd (typeDict, cp, &type) ;
      if (target && type)
	{
	  up = arrayp (aa, nn++, SNE) ;
	  up->chrom1 = target ; 
	  up->chrom2 = target ; 
	  up->a1 = x ;
	  up->a2 = x ;
	  up->isDown1 = TRUE ;
	  up->isDown2 = TRUE ;
	  up->a1 = x - 1 ;  /* last base before the error */
	  if (! strncmp (cp, "Ins", 3))
	    {
	      up->classe = 1 ;
	      up->a2 = x ;
	      dictAdd (typeDict, cp+3, &(up->typeM)) ;
	      dictAdd (typeDict, "-",  &(up->typeW)) ;
	    }
	  else if (! strncmp (cp, "Del", 3))
	    {
	      up->classe = 3 ;
	      up->a2 = x + strlen(cp+3) ;
	      dictAdd (typeDict, cp+3, &(up->typeW)) ;
	      dictAdd (typeDict, "-",  &(up->typeM)) ;
	    }
	  else if (cp[1] == '2' || cp[1] == '>')
	    {
	      up->a2 = x + 1 ; /* first base after the error */
	      up->classe = 0 ;
	      dictAdd (typeDict, cp+2,  &(up->typeM)) ;
	      cp[1] = 0 ;
	      dictAdd (typeDict, cp,  &(up->typeW)) ;
	    }
	  else
	    {
	      up->classe = 99 ;
	    }	  
	}
    }

  /* parse the candidate translocations */
  snpParseTranslocations (snp, aa, cls) ; 

  arraySort (aa, sneOrder) ;
  
  for (ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      up = arrp (aa, ii, SNE) ;
      memset (bufM, 'n', 2 * delta + 1) ;
      memset (bufW, 'n', 2 * delta + 1) ;
      bufM[2 * delta + 1] = 0 ;
      bufW[2 * delta + 1] = 0 ;

      /* collect the wild type DNA from the fasta file */
      dna1 = stackText (snp->targetDna, array (snp->target2dna, up->chrom1, int)) ;
      dna2 = stackText (snp->targetDna, array (snp->target2dna, up->chrom2, int)) ;
      a1 = up->a1 - delta ; a2 = up->a1 + delta ;
      b1 = 0 ; 
      if (a1 < 0) { b1 = - a1 ; a1 = 0 ; }
      if (a2 > array (snp->targetDnaLength, up->chrom1, int))
	a2 = array (snp->targetDnaLength, up->chrom1, int) ;
      memcpy (bufW, dna1 + a1 + b1, a2 - a1 + 1);

      /* copy the 5' part of bufM from bufW upstream from the SNP */
      memcpy (bufM, bufW, delta) ;
      /* insert the SNP */
      dx = 0 ; cq = bufM + delta ;
      switch (up->classe)
	{
	case 0: /* SNP */
	  ccp = dictName (typeDict, up->typeM) ;
	  *cq++ = ace_upper((int)*ccp) ;
	  bufW[delta] = ace_upper((int)bufW[delta]) ;
	  dx = 1 ;
	  {
	    char buf3[4] ;
	    int k ;
	    for (k=0 ; k < 3 ; k++)
	      buf3[k] = ace_upper((int)bufW[delta - 1 + k]) ;
	    buf3[3] = 0 ;
	    dictAdd (typeDict, buf3, &(up->typeW)) ;
	    buf3[1] = ace_upper((int)*ccp) ;
	    dictAdd (typeDict, buf3, &(up->typeM)) ;
	  }
	  break ;
	case 1: /* MICROINSERTION */ 
	case 2: /* INSERTION */ 
	  dictAdd (typeDict, "-",  &(up->typeW)) ;
	  ccp = dictName (typeDict, up->typeM) ;
	  while (*ccp)
	    {
	      *cq++ = ace_upper((int)*ccp++) ; 
	      dx++ ;
	    } 
	  break ;
	case 3: /* MICRODELETION */ 
	case 4: /* DELETION */
	case 5: /* INTRON */
	  dictAdd (typeDict, "-",  &(up->typeM)) ;
	  {
	    char buf3[16] ;
	    int k = up->a2 - up->a1 + 1 ;

	    if (! up->typeW)
	      {
		if (k < 13)
		  {
		    strncpy (buf3, dna1, k) ;
		    buf3[k] = 0 ;
		    dictAdd (typeDict, buf3, &(up->typeW)) ;
		  }
		else
		  {
		    strncpy (buf3, dna1, 9) ;
		    buf3[9] = '+' ;  buf3[10] = '+' ; buf3[11] = '+' ;  buf3[12] = 0 ; 
		    dictAdd (typeDict, buf3, &(up->typeW)) ;
		  }
	      }
	    k = up->a2 - up->a1 - 1 ;
	    if (k > delta) k = delta ;
	    if (k > 0)
	      while (k--)
		bufW[delta+k] = ace_upper((int)bufW[delta+k]) ;
	  }
	  break ;
	case 6: /* MICRODUPLICATION */ 
	case 7: /* DUPLICATION */ 
	  dictAdd (typeDict, "-",  &(up->typeW)) ;
          break ;
	case 8: /* INVERSION */
	case 9: /* TRANSLOCATION */
	default:
	  dictAdd (typeDict, "-",  &(up->typeM)) ;
	  dictAdd (typeDict, "-",  &(up->typeW)) ;
	  {
	    int k = delta ;
	    while (k--)
	      bufW[delta+k] = ace_upper((int)bufW[delta+k]) ;
	  }
	  break ;
	}

      /* collect the 3' part of the buffer with the sequence behind the SNP */
      dna2 = stackText (snp->targetDna, array (snp->target2dna, up->chrom2, int)) ;
      if (up->isDown2)
	memcpy (cq, dna2 + up->a2 - 1, delta - dx + 1) ;
      else
	{
	  int i ;
	  const char *ccq = dna2 + up->a2 - 1 ;
	  for (i = 0 ; i < delta - dx ; i++)
	    *cq++ = dnaDecodeChar [(int)complementBase(dnaEncodeChar[(int)*ccq--])];
	  *cq = 0 ;
	}

      /* export */
      switch (up->classe)
	{
	case 99: 
	  ao = 0 ;
	  break ;
	case 0: 
	case 1: 
	case 3: 
	case 6:
	  ao = aoSNP ;
	  break ;
	default:
	  ao = aoBRK ;
	  break ;
	}

      /* we give an identifier to quasi contiguous SNPs */
      if (up->chrom1 != oldChrom1 || up->a1 < oldA1 - delta || up->a1 > oldA1 + delta)
	snpId++ ;
      oldChrom1 = up->chrom1 ; oldA1 = up->a1 ;

      if (up->typeW  && strstr (dictName (typeDict, up->typeW), "nnnnnn")) continue ;
      if (up->typeM  && strstr (dictName (typeDict, up->typeM), "nnnnnn")) continue ;
      if (ao)
	{
	  /* export the mutant */
	  aceOutf (ao, ">M:%s:%d:%s:%d:%s:%s.%s|Gene|%s.%d\n%s\n"
		   , dictName(snp->targetDict, up->chrom1)
		   , up->a1
		   , dictName(snp->targetDict, up->chrom2)
		   , up->a2
		   , cls2[up->classe]
		   , up->typeW ? dictName (typeDict, up->typeW) : "-"
		   , up->typeM ? dictName (typeDict, up->typeM) : "-"
		   , snp->run, snpId
		   , bufM
		   ) ;
	  /* as a control, export also the wild type */
	  aceOutf (ao, ">W:%s:%d:%s:%d:%s:%s.%s|Gene|%s.%d\n%s\n"
		   , dictName(snp->targetDict, up->chrom1)
		   , up->a1
		   , dictName(snp->targetDict, up->chrom2)
		   , up->a2
		   , cls2[up->classe]
		   , up->typeW ? dictName (typeDict, up->typeW) : "-"
		   , up->typeM ? dictName (typeDict, up->typeM) : "-"
		   , snp->run, snpId
		   , bufW
		   ) ;
	}
    }
  
  fprintf (stderr, "Exported %d edited sequences of length %d", nok, 2 * delta + 1 ) ;
  ac_free (h) ;

  return nok ;
} /*  snpExportModifiedSequence */

/*************************************************************************************/
/*************************************************************************************/

typedef struct aesStruct { int snp, target, aceNam, vcfNam, a1, da, m, mm, mp, w, wm, wp, other, otherp, otherm, oo1m, oo1p, oo2m, oo2p, snipnet1, snipnet2 ; char mb, mb2, mb3, wb, wb2, wb3, type[8] ; int phasingScore[8], phase ; } AES ;



static int aesOrder (const void *va, const void *vb)
{
  const AES *a = (const AES *)va, *b = (const AES *)vb ;
  int n ;

  n = a->target - b->target ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->snp - b->snp ; if (n) return n ;

  return 0 ;
} /* aesOrder */

/*************************************************************************************/
/* parse a snp table exported previously by detect or count */
static int snpAliExtendParseSnpList (SNP *snp, Array raw)
{
  AC_HANDLE h = ac_new_handle () ;
  Array aa = arrayHandleCreate (100000, AES, snp->h) ;
  int nn = 0, target, s, a1, aceNam = 0, vcfNam = 0 ;
  ACEIN ai ;
  AES *aes ;
  char *cp ;
  DICT *targetDict = snp->targetDict  ;
  DICT *snpDict = snp->snpDict  ;
  DICT *snpNamDict = snp->snpNamDict ;
  char buf[1000] ;
  ai = aceInCreate (snp->snpListFileName, 0, h) ;
  if (! ai)
    messcrash ("Cannot open file -snp_list %s",  snp->snpListFileName) ;
  
  if (!snpDict)
    snpDict = snp->snpDict = dictHandleCreate (100000, snp->h) ;
 if (!snpNamDict)
    snpNamDict = snp->snpNamDict = dictHandleCreate (100000, snp->h) ;


  snp->aliSnps = aa ;
  while (aceInCard (ai))
    {
      cp = aceInWord (ai) ;
      if (cp && *cp == '#')
	continue ;
      if (!cp)
	{ aceInStep(ai, '\t') ; cp = aceInWord (ai) ; }
      if (!cp || *cp == '#')
	continue ;
    if (0 && ! dictFind (targetDict, cp, &target))
	{
	  if (snp->phasing)
	    {  /* optionally skip the snp number in the RESULT file */
	       aceInStep(ai, '\t') ;
	       cp = aceInWord (ai) ;
	       if (!cp || ! dictFind (targetDict, cp, &target))
		 continue ;
	    }
	  else
	    continue ; /* only use the targets for which we have a fasta file */
	}
      dictAdd (targetDict, cp, &target) ; 
      aceInStep(ai, '\t') ;
      if (!aceInInt(ai, &a1))
	continue ;
      aceInStep(ai, '\t') ;
      cp = aceInWord (ai) ;
      if (!strcmp (cp, "+") || !strcmp (cp, "-")) /* jump the optional strand column */
	{
	  aceInStep(ai, '\t') ;
	  cp = aceInWord (ai) ;
   	}

      aes = arrayp (aa, nn++, AES) ;
      aes->target = target ;
      aes->a1 = a1 ;
      if (cp[1] == '2' || cp[1] == '>')
	{
	  aes->wb = dnaEncodeChar [(int)cp[0]] ;
	  aes->mb = dnaEncodeChar [(int)cp[2]] ;
	  aes->type[0] = dnaDecodeChar [(int)aes->wb] ;
	  aes->type[1] = '>' ;
	  aes->type[2] = dnaDecodeChar [(int)aes->mb] ;
	  aes->type[3] = 0 ;
	  aes->da = 0 ;
	}
      else if (! strncasecmp (cp, "Del", 3))
	{
	  aes->type[0] = '-' ;
	  aes->mb = 0 ;
	  aes->wb = dnaEncodeChar [(int)cp[3]] ;
	  aes->type[1] = dnaDecodeChar [(int)aes->wb] ;
	  aes->da = - (strlen(cp) - 3) ;  /* deletion */
	  if (aes->da < -1)
	    {
	      aes->wb2 = dnaEncodeChar [(int)cp[4]] ;
	      aes->type[2] = dnaDecodeChar [(int)aes->wb2] ;
	    }
	  if (aes->da < -2)
	    {
	      aes->wb3 = dnaEncodeChar [(int)cp[5]] ;
	      aes->type[3] = dnaDecodeChar [(int)aes->wb3] ;
	    }
	  aes->type[1 - aes->da]  = 0 ;
	}
      else if (! strncasecmp (cp, "Ins", 3))
	{
	  aes->type[0] = '-' ;
	  aes->wb = 0 ;
	  aes->mb = dnaEncodeChar [(int)cp[3]] ;
	  aes->type[1] = dnaDecodeChar [(int)aes->mb] ;
	  aes->da = (strlen(cp) - 3) ;  /* deletion */
	  if (aes->da > 1)
	    {
	      aes->mb2 = dnaEncodeChar [(int)cp[4]] ;
	      aes->type[2] = dnaDecodeChar [(int)aes->mb2] ;
	    }
	  if (aes->da > 2)
	    {
	      aes->mb3 = dnaEncodeChar [(int)cp[5]] ;
	      aes->type[3] = dnaDecodeChar [(int)aes->mb3] ;
	    }
	  
	  aes->type[1 + aes->da]  = 0 ;
	}

      sprintf (buf, "%s:%d:%s", dictName(snp->targetDict,aes->target), aes->a1, aes->type) ;
      dictAdd (snpDict, buf, &s) ;
      aes->snp = s ; 
      /* parse the aceName */
      aceInStep(ai, '\t') ;  cp = aceInWord (ai) ;
      dictAdd (snpNamDict, cp, &(aceNam)) ; 
      /* parse the vcfName */
      aceInStep(ai, '\t') ;  cp = aceInWord (ai) ;
      dictAdd (snpNamDict, cp, &(vcfNam)) ;
      /* parse the genomic snpnet */
      if (! snp->phasing)
	{
	  aceInStep(ai, '\t') ;  cp = aceInWord (ai) ;
	  if (cp && *cp != '-' && ! strstr(cp, "NNNN"))
	    {
	      aes->snipnet1 = stackMark (snp->snipnet) ;
	      pushText (snp->snipnet, cp) ;
	    }
	  aceInStep(ai, '\t') ; cp = aceInWord (ai) ;
	  if (cp && *cp != '-' && ! strstr(cp, "NNNN"))
	    {
	      aes->snipnet2 = stackMark (snp->snipnet) ;
	      pushText (snp->snipnet, cp) ;
	    }
	}
      aes = arrayp (raw, s, AES) ;
      aes->target = target ;
      aes->a1 = a1 ;
      aes->aceNam = aceNam ;
      aes->vcfNam = vcfNam ;
    }

  arraySort (aa, aesOrder) ;

  ac_free (h) ;	
  return nn ;
} /* snpAliExtendParseSnpList */

/*************************************************************************************/

static int snpAliExtendGetHits (SNP *snp, ACEIN ai, Array hits, int *runQualityPrefixp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  char buf[1024] ;
  HIT *up ;
  char *cp, *cq ;
  BOOL isRead2 = FALSE ;

  buf[0] = 0 ;
  hits = arrayReCreate (hits, 32, HIT) ;
  *runQualityPrefixp = 1 ;
  stackClear (snp->aliExtendStack) ;
  pushText (snp->aliExtendStack, "toto") ;
  while (aceInCard (ai))
    {
      cp = aceInWord (ai) ;
      if (! cp || *cp == '#') 
	continue ;
      cq = cp + strlen (cp) - 1 ;
      isRead2 = FALSE ;
      switch ((int)(*cq))
	{
	case '<':
	  isRead2 = TRUE ;
	  /* fall through */
	case '>':
	  *cq = 0 ;
	  break ;
	}
      if (nn && buf[0] && strcmp (cp, buf))
	{
	  aceInCardBack (ai) ;
	  break ;
	}
      strncpy (buf, cp, 999) ;
      up = arrayp (hits, nn++, HIT) ; aceInStep(ai, '\t') ;
      aceInInt (ai, &(up->score)) ; aceInStep(ai, '\t') ;
      aceInInt (ai, &(up->mult)) ; aceInStep(ai, '\t') ;
      aceInInt (ai, &(up->ln)) ; aceInStep(ai, '\t') ;
      aceInInt (ai, &(up->ali)) ; aceInStep(ai, '\t') ;
      aceInInt (ai, &(up->x1)) ; aceInStep(ai, '\t') ;
      aceInInt (ai, &(up->x2)) ; aceInStep(ai, '\t') ;
      up->isRead2 = isRead2 ;
      cp = aceInWord (ai) ;
      if (cp && snp->target_class && strcmp (cp, snp->target_class))
	{ nn-- ; continue ; } aceInStep(ai, '\t') ;

      if (snp->runQuality)
	{
	  const char *ccq = buf + strlen(buf) - 1 ;
	  int i = 2 ;

	  if (*ccq == '>') *runQualityPrefixp = 1 ;
	  else if (*ccq == '<') *runQualityPrefixp = 2 ; 
	  else if(0)
	    {
	      while (ccq > buf && *ccq != '/' && i--) ccq-- ;
	      *runQualityPrefixp = (*ccq == '/') ? atoi (ccq+1) : 0 ;
	      if (*runQualityPrefixp < 0 || *runQualityPrefixp > 9)
		continue ; /* reject this probe */
	    }
	}
      if (up->ali < 140 && 
	  100 * up->ali < snp->minAliPerCent * up->ln)
	continue ;
      cp = aceInWord (ai) ;  aceInStep(ai, '\t') ;/* gene, drop */
      aceInInt (ai, &(up->unicity)) ; aceInStep(ai, '\t') ;
      cp = aceInWord (ai) ; /* target */
      if (!cp || ! dictFind (snp->targetDict, cp, &up->target))
	{ nn-- ; continue ; } aceInStep(ai, '\t') ;
      aceInInt (ai, &(up->a1)) ; aceInStep(ai, '\t') ;
      aceInInt (ai, &(up->a2)) ; aceInStep(ai, '\t') ;
      aceInInt (ai, &(up->nN)) ; aceInStep(ai, '\t') ;
      aceInInt (ai, &(up->nErr)) ;  aceInStep(ai, '\t') ;

      cp = aceInWord (ai) ;
      if (cp && strcmp (cp, "-"))
	{
	  up->errTag = stackMark (snp->aliExtendStack) ;
	  pushText (snp->aliExtendStack, cp) ;
	}
      else
	up->errTag = 0 ;

      aceInStep(ai, '\t') ; cp = aceInWord (ai) ;
      if (cp && strcmp (cp, "-"))
	{
	  up->errTarget = stackMark (snp->aliExtendStack) ;
	  pushText (snp->aliExtendStack, cp) ;
	}
      else
	up->errTarget = 0 ;

      aceInStep(ai, '\t') ; cp = aceInWord (ai) ;
      if (cp && strcmp (cp, "-"))
	{
	  up->prefix = stackMark (snp->aliExtendStack) ;
	  pushText (snp->aliExtendStack, cp) ;
	}
      else
	up->prefix = 0 ;

      aceInStep(ai, '\t') ; cp = aceInWord (ai) ;
      if (cp && strcmp (cp, "-"))
	{
	  up->suffix = stackMark (snp->aliExtendStack) ;
	  pushText (snp->aliExtendStack, cp) ;
	}
      else
	up->suffix = 0 ;
      aceInStep(ai, '\t') ; cp = aceInWord (ai) ; /* target prefix, drop it */
     aceInStep(ai, '\t') ; cp = aceInWord (ai) ; /* target suffix, drop it */
      up->dPair = 0 ;
      aceInStep(ai, '\t') ; aceInInt (ai, &(up->dPair)) ;
      if (up->dPair < 0 && up->dPair > -20)
	break ; /* do not use mate pairs for phasing is not a good pair */
      if (up->unicity > 1)
	break ; /* do not use mate pairs because w e may mix up 2 mapping of same read, but the single read remains ok */
   }

  arrayMax (hits) = nn ;
  if (! nn)
    ac_free (ai) ;

  ac_free (h) ;	
  return nn ;
} /* snpAliExtendGetHits */

/*************************************************************************************/

static int snpAlignExtendCheck (SNP *snp, Array err, HIT *up
				, Array tagZone, int y1, int y2
				, Array targetZone, int a0, int a1, int a2, int a3, int b1, int b2
				, Array aliSnps, int jaes
				, Array badAes
				, int runQualityPrefix
				, BOOL isDown
				, float *supportp
				) 
{
  A_ERR *ep ;
  int iep, xx1, xx2, bb1, bb2 ;
  int a1a, a1b ; /* feet to bridge around an indel */
  int iaes, nbad = 0, x ;
  AES *aes ;

  bb1 = b1 ; bb2 = b2 ;
  xx1 = y1 ; xx2 = y2 ;
  err = aceDnaDoubleTrackErrors (tagZone, &xx1, &xx2, TRUE,
				 targetZone, 0, &bb1, &bb2,
				 0, err, 3, 0, FALSE, 0) ;
  
  /* compare the errors with the relevant list of snps and 
   * increment the support for  the wild type or for the mutation
   */
  
  for (iaes = jaes, aes = arrp (aliSnps, iaes, AES) ; iaes < arrayMax (aliSnps)  && aes->a1 <= a2 ; aes++, iaes++)
    {
      int isM = 0 ;
      
      if (aes->target != up->target)
	break ;
      if (a0 > aes->a1 || a3 < aes->a1)
	continue ;
      if (aes->da <= 0 &&  badAes && aes->wb != arr (targetZone, aes->a1 - (a1 - b1) - 1, char))
	continue ;
      a1a = a1b = aes->a1 ;
      if (aes->da > 0) /* insertion */
	{
	  char cm = aes->mb, *cp = arrp (targetZone, aes->a1 - (a1 - b1) - 1, char) ;
	  cp -= aes->da ;
	  while (cm == *cp && a1a > a0) { cp -= aes->da ; a1a -= aes->da ; }
	  a1a -= aes->da ;
	  if (a1a - (a1 - b1) < b1 || a1b - (a1 - b1) > b2)
	    continue ;
	}
      if (aes->da < 0) /* deletion */
	{
	  char cm = aes->wb, *cp = arrp (targetZone, aes->a1 - (a1 - b1) - 1, char) ;
	  while (cm == *cp && a1a > a0) { cp += aes->da ; a1a += aes->da ; }

	  a1b -= aes->da ;
	  if (a1a - (a1 - b1) < b1 || a1b - (a1 - b1) > b2)
	    continue ;
	}
      if (
	  (a1a - (a1 - b1) >= b1 && a1a - (a1 - b1) < bb1) ||
	  (a1b - (a1 - b1) <= b2 && a1b - (a1 - b1) > bb2)
	  )
	isM = -1 ;
      x = -1 ;
      if (arrayMax (err))
	for (ep = arrp (err, 0, A_ERR), iep = 0 ; isM == 0 && iep < arrayMax (err) ; iep++, ep++)
	  {
	    if (aes->a1 - (a1 - b1) - 1 == ep->iLong)
	      {
		isM = -1 ;
		switch (ep->type)
		  {
		  case ERREUR:  
		    if (aes->da == 0 &&
			aes->mb == ep->baseShort)
		      isM = 1 ;
		    break ;
		  case TROU:
		    if (aes->da == -1 && 
			arr(targetZone,ep->iLong,char) == aes->wb
			)
		      isM = 1 ;
		    break ;
		  case TROU_DOUBLE: 
		    if (aes->da == -2 && 
			arr(targetZone,ep->iLong,char) == aes->wb 
			&& arr(targetZone,ep->iLong+1,char) == aes->wb2
			) 
		      isM = 1 ; 
		    break ;
		  case TROU_TRIPLE:  
		    if (aes->da == -3 && 
			arr(targetZone,ep->iLong,char) == aes->wb && 
			arr(targetZone,ep->iLong+1,char) == aes->wb2  && 
			arr(targetZone,ep->iLong+2,char) == aes->wb3
			) 
		      isM = 1 ; 
		    break ;
		  case INSERTION: 
		    if (aes->da == 1 && 
			ep->baseShort == aes->mb
			) 
		      isM = 1 ; 
		    break ;
		  case INSERTION_DOUBLE:   
		    if (aes->da == 2 && 
			ep->baseShort == aes->mb && 
			arr(tagZone,ep->iShort+1,char) == aes->mb2
			) 
		      isM = 1 ;
		    break ;
		  case INSERTION_TRIPLE:   
		    if (aes->da == 3 && 
			ep->baseShort == aes->mb  && 
			arr(tagZone,ep->iShort+1,char) == aes->mb2 && arr(tagZone,ep->iShort+2,char) == aes->mb3
			) 
		      isM = 1 ;
		    break ;
		  default: 
		    break ;
		  }
		if (isM == 1) 
		  x = ep->baseShort ;
	      }
	    else if (a1b - (a1 - b1) - 1 < ep->iLong && (iep == 0 || (ep-1)->iLong < a1b - (a1 - b1) - 1 - 3))
	      { isM = 0 ; break ; }
	    else if (a1b - (a1 - b1) - 1 > ep->iLong + 3 && iep == arrayMax (err) - 1)
	      { isM = 0 ; break ; }
	  }
      
      if (badAes)
	{
	  if (supportp) messcrash ("(badAes && supportp) are incompatible options in snpAlignExtendCheck") ;
	  x = 20 ; /* we do not yet know the correct position */
	  if (isM == 1) 
	    {
	      double z = up->mult * runQ (snp->runQuality, runQualityPrefix, x) ; /* attention we did not compensate for indels */
	      aes->m += z ;	 
	      if (isDown) aes->mp += z ;
	      else aes->mm += z ;
	    }
	  else if (isM == 0)
	    {
	      double z = up->mult * runQ (snp->runQuality, runQualityPrefix, x) ; /* attention we did not compensate for indels */
	      aes->w += z ;
	      if (isDown) aes->wp += z ;
	      else aes->wm += z ;
	    }
	  if (isM == -1)  /* unresolved cases */
	    keySet (badAes, nbad++) = iaes ;
	}
      else
	{
	  int x = 20 ; 
	  if (!supportp) messcrash ("(!badAes && !supportp) are incompatible options in snpAlignExtendCheck") ;
	  *supportp = up->mult * runQ (snp->runQuality, runQualityPrefix, x) ; /* attention we did not compensate for indels */
	  return isM ;
	}
    }

  return nbad ;
} /* snpAlignExtendCheck */

/*************************************************************************************/

static int snpAliExtendAnalyse (SNP *snp, Array hits, int runQualityPrefix)
{
  AC_HANDLE h = ac_new_handle () ;
  HIT *up ;
  int i, ii, ns, nn = 0, nbad = 0 ;
  int iaes, jaes ;
  int a0, a1, a2, a3, da, aL ;
  int b0, b1, b2 ;
  int x1, x2, ddx ;
  int y1, y2 ;
  float support = 0 ;

  const char *ccp = 0 ;
  const char *targetDna = 0 ;
  static Array err = 0, tagZone = 0, targetZone = 0 ;
  Array aliSnps = snp->aliSnps ;  
  BOOL isDown ;
  AES *aes ;
  static KEYSET badAes = 0 ;

  if (! aliSnps || ! arrayMax (aliSnps))
    return 0 ;

  for (ii = 0, up = arrp (hits, 0, HIT) ; ii < arrayMax (hits) ; ii++, up++)
    {
      nbad = 0 ; badAes = keySetReCreate (badAes) ;

      /* switch to the positive strand of the target */
      isDown = TRUE ;
      if (up->a2 < up->a1)
	{
	  int a, x ;

	  isDown = FALSE ;
	  a = up->a1 ; up->a1 = up->a2 ; up->a2 = a ;
	  /*
	   * we DO NOT reverse the tag coordinates, because we will stay in target coords
	   * but we switch the prefix and the suffix
	   * we do not complement since they are always oriented outwards
	   */
	  x = up->prefix ; up->prefix = up->suffix ; up->suffix = x ;
	}
      a1 = up->a1 - (up->prefix ? strlen(stackText (snp->aliExtendStack, up->prefix)) : 0) ;
      a2 = up->a2 + (up->suffix ? strlen(stackText (snp->aliExtendStack, up->suffix)) : 0) ;
 
      /* check for candidate SNPs in the area */
      ns = 0 ;
      for (iaes = 0, jaes = -1, aes = arrp (aliSnps, 0, AES) ; iaes < arrayMax (aliSnps) ; aes++, iaes++)
	{
	  if (aes->a1 > a2 || aes->a1 < a1) 
	    continue ;
	  if (aes->target != up->target)
	    continue ;
	  ns++ ;
	  if (jaes < 0) jaes = iaes ;
	}
      if (! ns)
	continue ;
      /* recover the correct zone of the target DNA */
      targetDna = 0 ;
      targetDna = stackText (snp->targetDna, array (snp->target2dna, up->target, int)) ;
      if (! targetDna)
	continue ;
      aL = strlen (targetDna) ; 

      /* a.. are the coordinates in the target
       * b.. are there projections in targetZone (so between 1 and maybe 200)
       * x.. are the coordinates on the read
       * y.. are their projections in the tag zone
       */

      a1 = up->a1 ; a2 = up->a2 ; da = a2 - a1 + 1 ;
      a0 = a1 > 50 ? a1 - 50 : 1 ;
      a3 = a2 + (1000 < aL ? a2 + 1000 : aL) ; 

      /* dx =  dy = da ; */

      b1 = 51 ; b2 = b1 + (a2 - a1) ;
      b0 = b1 - (a1 - a0) ; 

      x1 = up->x1 ; x2 = up->x2 ;
      /*
	b3 = b1 + (a3 - a0) ;

	dx = x2 - x1 + 1 ;
	x0 = b0 ; 
	x3 = b3 ;
      */

      y1 = 51 ; y2 = y1 + (x2 - x1) ;
      /* y0 = y1 - (x1 - x0) ;  y3 = y1 + (x3 - x0) ; */

      /* create the buffers */
      targetZone = arrayReCreate (targetZone, da + 202, char) ;
      array (targetZone, da + 101, char) = 0 ;
      memset (arrp(targetZone, 0, char), 'n', da + 100) ;

      tagZone = arrayReCreate (tagZone, da + 202, char) ;
      array (tagZone, da + 101, char) = 0 ;
      memset (arrp(tagZone, 0, char), 'n', da + 100) ;

     /* at this point we are in atgc, not in A_ T_ G_ C_ */
      memcpy (arrp(targetZone, b0 - 1, char), targetDna + a0 - 1, a3 - a0 + 1) ;
    
      /* recover the zone of the tag matching the target */
      /* copy the target */
      memcpy (arrp(tagZone, 0, char), arrp (targetZone, 0, char), da + 101) ;

      /* most likely the prefix/suffix are null or were clipped because they differ from the target */
      /* add the prefix */ 
      ccp = up->prefix ? stackText (snp->aliExtendStack, up->prefix) : 0 ;
      if (ccp)
	{
	  int ln = strlen (ccp), ln1 = 0 ; 

	  if (ln > 30) ln = 30 ; 
	  i = 0 ; ccp-- ;
	  while (ln > 0 && b1 - i >= 1 && *++ccp)
	    {
	      char cc = dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)*ccp])] ;
	      arr (tagZone, b1 - i - 2, char) = cc ; ln1++ ; i++ ;
	    }
	  b1 -= ln1 ; a1 -= ln1 ; x1 -= ln1 ; y1 -= ln1 ;
	}      
      /* add the suffix */
      ccp = up->suffix ? stackText (snp->aliExtendStack, up->suffix) : 0 ;
      if (ccp)
	{
	  int ln = strlen (ccp) ; 

	  if (ln > 30) ln = 30 ;
	  memcpy (arrp (tagZone, b2, char), ccp, ln) ;
	  arr (tagZone, b2 + ln, char) = 0 ;
	  b2 += ln ; a2 += ln ; x2 += ln ; y2 += ln ;
	}
      
      /* add the differences between the tag and the target */
      ddx = 0 ;
      ccp = up->errTarget ? stackText (snp->aliExtendStack, up->errTarget) : 0 ;
      if (ccp)
	{  /* scan the errors and edit the tag sequence accordingly, possibly edit dx in case of indel */
	  char cc, *cp3, *cp2, *cp1 = strnew (ccp, 0) ;
	  int x ;
	  
	  while (cp1)
	    {
	      cp2 = strchr (cp1, ',') ;
	      if (cp2) *cp2 = 0 ;
	      cp3 = strchr (cp1, ':') ;
	      if (!cp3) break ;
	      *cp3++ = 0 ; x = atoi(cp1) ;
	      if (x < a1 || x  > a2) break ;
	      while (*cp3 == '*') cp3++ ;
	      if (cp3[1] == '>') /* substitution, in target coordinates */
		{
		  cc = ace_lower (cp3[0]) ;
		  if (cp3[0] != ace_lower(array(targetZone, x - (a1 - b1) - 1, char)))
		    break ;
		  cc = cp3[2] ;
		  array(tagZone, x - (a1 - b1) - 1, char) = cc ;
		}
	      else if (cp3[3] == 'o') ; /* oo transition in SOLiD ignore */
	      else if (cp3[0] == '-' && cp3[1] == '-' && cp3[2] == '-')   /* triple deletion */
		{
		  char * cp ;
		  
		  cc = ace_lower (cp3[3]) ;
		  if (cc != ace_lower(array(targetZone, x - (a1 - b1) - 1, char)))
		    break ;
		  cp = arrp (tagZone, x - (a1 - b1) - 1, char) ;
		  while (*cp) { *cp = *(cp+3) ; cp++ ; }
		  ddx -= 3 ;
		}
	      else if (cp3[0] == '-' && cp3[1] == '-')   /* double deletion */
		{
		  char * cp ;

		  cc = ace_lower (cp3[2]) ;
		  if (cc != ace_lower(array(targetZone, x - (a1 - b1) - 1, char)))
		    break ;
		  cp = arrp (tagZone, x - (a1 - b1) - 1, char) ;
		  while (*cp) { *cp = *(cp+2) ; cp++ ; }
		  ddx -= 2 ;
		}
	      else if (cp3[0] == '-')   /* deletion */
		{
		  char * cp ;
		  
		  cc = ace_lower (cp3[1]) ;
		  if (cc != ace_lower(array(targetZone, x - (a1 - b1) - 1, char)))
		    break ;
		  cp = arrp (tagZone, x - (a1 - b1) - 1, char) ;
		  while (*cp) { *cp = *(cp+1) ; cp++ ; }
		  ddx-- ;
		}
	      else if (cp3[0] == '+' && cp3[1] == '+' && cp3[2] == '+')   /* triple insertion */
		{
		  char *cp, *cq ;
		  
		  cc = ace_lower (cp3[3]) ;
		  cp = arrp (tagZone, x - (a1 - b1) - 1, char) ;
		  for (cq = cp ; *cq ; cq++) ;
		  cq++ ; cq++ ; *(cq+1) = 0 ;
		  while (cq > cp) { *cq = *(cq-3) ; cq-- ; }
		  *cp = cc ; *(cp+1) = ace_lower (cp3[4]) ; *(cp+2) = ace_lower (cp3[5]) ;
		  ddx += 3 ;
		}
	      else if (cp3[0] == '+' && cp3[1] == '+')   /* double insertion */
		{
		  char *cp, *cq ;
		  
		  cc = ace_lower (cp3[2]) ;
		  cp = arrp (tagZone, x - (a1 - b1) - 1, char) ;
		  for (cq = cp ; *cq ; cq++) ;
		  cq++ ; *(cq+1) = 0 ;
		  while (cq > cp) { *cq = *(cq-2) ; cq-- ; }
		  *cp = cc ; *(cp+1) = ace_lower (cp3[3]) ;
		  ddx += 2 ;
		}
	      else if (cp3[0] == '+')   /* insertion */
		{
		  char *cp, *cq ;

		  cc = ace_lower (cp3[1]) ;
		  cp = arrp (tagZone, x - (a1 - b1) - 1, char) ;
		  for (cq = cp ; *cq ; cq++) ;
		  *(cq+1) = 0 ;
		  while (cq > cp) { *cq = *(cq-1) ; cq-- ; }
		  *cp = cc ;
		  ddx++ ;
		}
	      
	      
	      cp1 = cp2 ? cp2 + 1 : 0 ;
	    }
	}
      arrayMax (tagZone) += ddx ;
      
      /* encode the target and the tag and align, extending in both directions */
      dnaEncodeArray (targetZone) ;  dnaEncodeArray (tagZone) ;
      
      err = arrayReCreate (err, 32, A_ERR) ;
      nbad = snpAlignExtendCheck (snp, err, up,  tagZone, y1, y2, targetZone, a0, a1, a2, a3, b1, b2
				  , aliSnps, jaes, badAes, runQualityPrefix, isDown, 0) ;

      if (nbad) /* try to realign on a modified genome and see if the error goes away */
	{
	  Array newTargetZone = 0 ;

	  int ibad, isM ;
	  for (ibad = 0 ; ibad < nbad ; ibad++)
	    {
	      int i, dda = 0 ;
	      char *cp, *cp3 ;

	      aes = arrp (aliSnps, keySet (badAes, ibad), AES) ; 	 
	      if (a0 > aes->a1 || a3 < aes->a1)
		continue ;
	      if (aes->da <= 0 &&  aes->wb != arr (targetZone, aes->a1 - (a1 - b1) - 1, char))
		continue ;
	      newTargetZone = arrayCopy (targetZone) ;

	      cp3 = aes->type ;
	      if (cp3[1] == '2')
		arr (newTargetZone, aes->a1 - (a1 - b1) - 1, char) = aes->mb ;
	      else if (!strncmp(cp3,"Del",3))   /* single deletion */
		{
                  dda = strlen(cp3) - 3 ;
		  for (i = aes->a1 - (a1 - b1) - 1, cp = arrp (newTargetZone, i, char) ; i < arrayMax(newTargetZone) - dda ; cp++, i++)
		    *cp = *(cp + dda) ;
		  dda = -dda ;
		}
	      else if (!strncmp(cp3,"Ins",3))   /* single deletion */
		{
                  dda = strlen(cp3) - 3 ;
		  for (i = arrayMax(newTargetZone) - 1 , cp = arrp (newTargetZone, i, char) ; i >= aes->a1 - (a1 - b1) - 1 + dda ; cp--, i--)
		    *cp = *(cp - dda) ;
		  arr (newTargetZone, aes->a1 - (a1 - b1) - 1, char) = aes->mb ;

		  if (dda >= 2) arr (newTargetZone, aes->a1 - (a1 - b1) - 1 + 1, char) = aes->mb2 ;
		  if (dda >= 3) arr (newTargetZone, aes->a1 - (a1 - b1) - 1 + 2, char) = aes->mb3 ;
		}
	      isM = snpAlignExtendCheck (snp, err, up,  tagZone, y1, y2
					 , newTargetZone, a0, a1, a2 + dda, a3 + dda, b1 + dda, b2 + dda
					 , aliSnps, keySet (badAes, ibad), 0, runQualityPrefix, isDown
					 , &support) ;
	      if (isM == 0) /* the edition of the genome fixed the error */
		{
		  aes->m += support ;	 
		  if (isDown) aes->mp += support ;
		  else aes->mm += support ;
		}
	      else if (isM == 1) /* the edition of the genome fixed the error */
		{
		  aes->w += support ;
		  if (isDown) aes->wp += support ;
		  else aes->wm += support ;
		}
	      else if (isM == -1) 
		{
		  aes->other += support ;
		  if (isDown) aes->otherp += support ;
		  else aes->otherm += support ;
		}
	      arrayDestroy (newTargetZone) ;
	    }
	}      
    }   /* loop on all hits */

  ac_free (h) ;	
  return nn ;
} /* snpAliExtendAnalyse */

/*************************************************************************************/

static int snpAliExtendReport (SNP *snp, ACEOUT ao)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, iaes ;
  AES *aes ;
  Array aliSnps = snp->aliSnps ;
  int cover, m, w, mp, mm, wp, wm, otherp, otherm, calledp, calledm ;
  int oo1p = 0, oo2p = 0, oo1m = 0, oo2m = 0 ;
  double z, z1, z2 ;

  for (iaes = 0, aes = arrp (aliSnps, 0, AES) ; iaes < arrayMax (aliSnps) ; aes++, iaes++)
    {
      m = aes->m/100 ; w = aes->w/100 ; otherp = (aes->otherp)/100 ; otherm = (aes->otherm)/100 ;
      calledp = (aes->mp + aes->wp + aes->otherp)/100 ;
      calledm = (aes->mm + aes->wm + aes->otherm)/100 ;
      if (calledp + calledm < snp->minCover) continue ;

      mm = aes->mm/100 ; wm = aes->wm/100 ;
      mp = aes->mp/100 ; wp = aes->wp/100 ;
      oo1p = aes->oo1p/100 ; oo1m = aes->oo1m/100 ;
      oo2p = aes->oo2p/100 ; oo2m = aes->oo2m/100 ;
      cover = (aes->m + aes->w + aes->other + aes->oo1p + aes->oo1m + aes->oo2p + aes->oo2m)/100 ;

      aceOutf (ao, "%s\t%d\t%s\t%s\t%s\t%s"
	       , dictName (snp->targetDict, aes->target)
	       , aes->a1, aes->type
	       , aes->snipnet1 ? stackText (snp->snipnet, aes->snipnet1) : "-"
	       , aes->snipnet2 ? stackText (snp->snipnet, aes->snipnet2) : "-"
	       , snp->run ? snp->run : "-"
	       ) ;

      if (1)
	{
	  snpNotLowShow (ao, m, calledp + calledm, &z, &z1, &z2, 0) ;
	  aceOutf (ao, "\t%.1f\t%d\t%d\t%d\t%d\t%d\t%d\t%d"
		   , z
		   , cover, m, w
		   , mp, wp, mm, wm
		   ) ;
	}
      
      if (1)
	{
	  if (snp->solid) /* solid */
	    aceOutf (ao, "\t%d:%d\t%d:%d"
		     , oo1p, oo2p
		     , oo1m, oo2m
		     ) ;
	  else
	    aceOutf (ao, "\t%d\t%d"
		     , oo1p + oo2p
		     , oo1m + oo2m
		     ) ;
	}
      if (1)
	{
	  const char *ss = "Both_strands" ; 
	  float cz2 = 0 ;
	  if (mp + wp + mm + wm < snp->minCover) ss = "-" ;
	  else if (2 * (mp + wp) < snp->minCover) ss = "Minus_strand_only" ;
	  else if (2 * (mm + wm) < snp->minCover) ss = "Plus_strand_only" ;
	  else
	    { 
	      double x, y ;
	      chi2 (mp, wp, mm, wm, &cz2) ;
	      x = mp/(1.0 + mp + (double)wp) ;
	      y = mm/(1.0 + mm + (double)wm) ;
	      if (cz2 > 100 && (x > y + .5 || x < y - .5)) ss = "Incompatible_strands" ; /* extremely high  and clear contradiction */ 
	    }

	  aceOutf (ao, "\t%s\t%.1f", ss, cz2) ;
	}
      aceOutf (ao, "\t%s"
	       , 100 * (oo1p + oo1m + oo2p + oo2m) > 40 * cover ? "N_rich" : (10 * (otherp + otherm) > cover ? "Several_forms" : "-")
	       ) ;
      
      snpNotLowShow (ao, mp, calledp, &z, &z1, &z2, 0) ;
      snpNotLowShow (ao, mm, calledm, &z, &z1, &z2, 0) ;
      aceOutf (ao, "\n") ;
     }
  ac_free (h) ;	
  return nn ;
} /* snpAliExtendReport */

/*************************************************************************************/

static int snpAliExtend (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  int runQualityPrefix = 0 ;
  Array hits = arrayHandleCreate (32, HIT, h) ;
  ACEIN ai = aceInCreate (snp->inFileName, snp->gzi, snp->h) ;
  ACEOUT ao = aceOutCreate (snp->outFileName, snp->unique ? ".extend.u.txt" : ".extend.nu.txt", snp->gzo, snp->h) ;
  
  snp->snipnet = stackHandleCreate (1000000, h) ;
  snp->aliExtendStack = stackHandleCreate (1024, h) ;
  while (snpAliExtendGetHits (snp, ai, hits, &runQualityPrefix))
    snpAliExtendAnalyse (snp, hits, runQualityPrefix) ;

  snpAliExtendReport (snp, ao) ;

  ac_free (h) ;	
  return nn ;
} /* snpAliExtend */

/*************************************************************************************/
/*************************************************************************************/
/* snp1 and snp2 are offsets in the  snp->aliSnps Array */
typedef struct phasingStruct { int snp1, snp2, a1, a2, mm, ww, mw, wm, m1, w1, m2, w2 ; } PHS ;

static int phsOrder (const void *va, const void *vb)
{
  const PHS *a = (const PHS *)va, *b = (const PHS *)vb ;
  int n ;

  n = a->snp1 - b->snp1 ; if (n) return n ;
  n = a->snp2 - b->snp2 ; if (n) return n ;

  return 0 ;
} /* phsOrder */

/*************************************************************************************/

static long int snpPhasingCompress (BigArray ph)
{
  if (bigArrayMax (ph) > 1)
    {
      PHS *xp, *yp ;
      long int iph, jph = 0 , iMax = bigArrayMax (ph) ;
      
      bigArraySort (ph, phsOrder) ;
      
      for (iph = 1, jph = 1, yp = bigArrp (ph, 0, PHS), xp = yp + 1 ; iph < iMax ; iph++, xp++)
	{
	  if (xp->snp1 == yp->snp1 &&
	      xp->snp2 == yp->snp2
	      )
	    {
	      yp->mm += xp->mm ;
	      yp->mw += xp->mw ;
	      yp->wm += xp->wm ;
	      yp->ww += xp->ww ;
	      yp->m1 += xp->m1 ;
	      yp->m2 += xp->m2 ;
	      yp->w1 += xp->w1 ;
	      yp->w2 += xp->w2 ;
	    }
	  else
	    {
	      yp++; jph++ ;
	      *yp = *xp ;
	    }
	}
      bigArrayMax (ph) = jph ;
    }
  return bigArrayMax (ph) ;
} /* snpPhasingCompress  */

/*************************************************************************************/
static int snpPhasingExport (SNP *snp, BigArray ph, Array raw)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEOUT ao = aceOutCreate (snp->outFileName, ".phasing", snp->gzo, h) ;
  long int da = 0, iph = 0, iMax = bigArrayMax (ph) ;
  PHS *up ;
  DICT *dict = snp->snpDict ;
  DICT *phaseDict = dictHandleCreate (10000, h) ;
  int pass, mm, ww, mw, wm, m1, w1, m2, w2, M1, M2, W1, W2, pb, line ;

  mm = ww = mw = wm = m1 = w1 = m2 = w2 =  M1 = W1 = M2 = W2 = pb = line = 0 ;

  aceOutDate (ao, "###", "Phasing info") ;
  aceOutf (ao, "# SNP1 ....\t SNP2 ....\tmm\tww\tmw\twm\tm1\tw1\tm2\tw2\tam\taw\tbm\tbw\tdistance\tk>3\ta/a or b/b\n") ;

  for (pass = pb = 0 ; pass < 3 ; pb = 0, pass++)
    {
      if (pass == 1)
	{ /* reinitialize the phasingScore counters */
	  int ii, iMax = arrayMax (raw) ;

	  for (ii = 0 ; ii < iMax ; ii++)
	    {
	      int a,b,c ;
	      AES *ap = arrayp (raw, ii, AES) ;
	      a = ap->phasingScore[0] ;
	      b = ap->phasingScore[1] ;
	      c = ap->phasingScore[2] ;
	      
	      if (a == 0 && c > 0 && b >= 0)
		ap->phasingScore[3] =  ap->phasingScore[2] ;
	      ap->phasingScore[0] = 0 ;
	      ap->phasingScore[1] = 0 ;
	      ap->phasingScore[2] = 0 ;
	      ap->phasingScore[6] = 0 ;
	    }
	  continue ;
	}

      for (iph = 0, up = bigArrp (ph, 0, PHS) ; iph < iMax ; iph++, up++)
	{
	  AES *ap = arrayp (raw, up->snp1, AES) ;
	  AES *bp = arrayp (raw, up->snp2, AES) ;
	  int cis, trans, ok ;
	  BOOL isH ;
	  char *kk ;
	  
	  cis = up->mm +  up->ww ;
	  trans = up->mw +  up->wm ;
	  kk = "" ; ok = 5 ;
	  if (cis * trans > 0)
	    {
	      kk = "**" ; pb++ ; ok = 2 ;		
	      if (
		  (20 * cis < trans) || 
		  (20 * trans < cis)
		  ) 
		{ kk = "*" ; ok = 1 ; }
	    }
	  if ( ok == 5 &&
	       cis + trans >= 3
	       )
	    { ok = 0 ; kk = "++" ; }
	  if (ap->phasingScore[3] + bp->phasingScore[3])
	    { kk = "***" ; ok = 5 ; }

	  isH = (ap->m >= 10 * ap->w) || (bp->m >= 10 * bp->w) || (ap->w >= 10 * ap->m) || (bp->w >= 10 * bp->m) ;

	  if (isH)
	    ok = 5 ;
	  if (ok < 3)
	    {
	      (ap->phasingScore[ok])++ ;
	      (bp->phasingScore[ok])++ ;
	      if (ok == 0)
		{
		  (ap->phasingScore[6]) += cis + trans ;
		  (bp->phasingScore[6]) += cis + trans ;
		}	   
	    }
	  if (pass == 0)
	    continue ;
	  
	  if (isH)
	    continue ;

	  if (pass == 2 && ok == 0) /* phasing */
	    {
	      int p1 = ap->phase ;
	      int p2 = bp->phase ;
	      
	      if (1)
		{
		  if (p1 && ! p2) { p2 = (cis > trans ? p1 : -p1) ; }
		  else if (p2 && ! p1) { p1 = (cis > trans ? p2 : -p2) ; }
		  else if (p1 == 0 &&  p2 == 0)
		    {
		      char *buf = strnew (dictName (dict, up->snp1), 0) ;
		      char *cr = strchr (buf, ':') ;
		      if (cr) *cr = '_' ;
		      cr = strchr (buf, ':') ;
		      if (cr) *cr = 0 ;	
		      dictAdd (phaseDict, buf, &p1) ;
		      p2 = (cis > trans ? p1 : -p1) ;
		    }
		}
	      ap->phase = p1 ;
	      bp->phase = p2 ;
	    }

	  aceOutf (ao, "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\td=%d\t%s%s\t%d\t%d\n"
		   , dictName (dict, up->snp1)
		   , dictName (dict, up->snp2)
		   , up->mm
		   , up->ww
		   , up->mw
		   , up->wm
		   , up->m1, up->w1
		   , up->m2, up->w2
		   , ap->m, ap->w
		   , bp->m, bp->w
		   , bp->a1 - ap->a1
		   , kk
		   , isH ? "H" : ""
		   , ap->phase
		   , bp->phase
		   ) ;
	  line++ ;
	  mm += up->mm ;
	  ww += up->ww ;
	  mw += up->mw ;
	  wm += up->wm ;
	  m1 += up->m1 ;
	  w1 += up->w1 ;
	  m2 += up->m2 ;
	  w2 += up->w2 ;
	  da += bp->a1 - ap->a1 ; 
	  M1 += ap->m ;
	  W1 += ap->w ;
	  M2 += bp->m ;
	  W2 += bp->w ;

	}
    }
  if (line == 0) line = 1 ;
  aceOutf (ao, "# cumul ....\tcumul ....\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%ld\t%d\n"
	   , mm
	   , ww
	   , mw
	   , wm
	   , m1, w1
	   , m2, w2
	   , M1, W1
	   , M2, W2
	   , da / line
	   , pb
	   ) ;
  aceOutf (ao, "# average ..\taverage ..\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%ld\t%d\n"
	   , mm/ line
	   , ww/ line
	   , mw/ line
	   , wm/ line
	   , m1/ line, w1/ line
	   , m2/ line, w2/ line
	   , M1/ line, W1/ line
	   , M2/ line, W2/ line
	   , da / line
	   , 0
	   ) ;

  ac_free (ao) ;
  ao = aceOutCreate (snp->outFileName, ".phasing_score", snp->gzo, h) ;
  if (ao)
    {
      int ii, iMax = arrayMax (raw) ;
      
      aceOutDate (ao, "###", "Phasing success and arrors") ;
      aceOutf (ao, "# SNP\tPhases correctly with n other SPS using at least 3 molecules\tMild problems\tMore that 5%% cis/trans or vice versa\tOnly problematic links\n") ;

      for (ii = 0 ; ii < iMax ; ii++)
	{
	  int a,b,c,d ;
	  AES *ap = arrayp (raw, ii, AES) ;
	  a = ap->phasingScore[0] ;
	  b = ap->phasingScore[1] ;
	  c = ap->phasingScore[2] ;
	  d = ap->phasingScore[3] ;

	  if (a + b + c + d)
	    aceOutf (ao, "%s\t%d\t%d\t%d\t%d\n"
		     , dictName (dict, ii)
		     , a, b, c, d
		     ) ;
	}

      ac_free (ao) ;
    }

  ac_free (ao) ;
  ao = aceOutCreate (snp->outFileName, ".phases.ace", snp->gzo, h) ;
  if (ao)
    {
      int ii, iMax = arrayMax (raw) ;
      DICT *snpNamDict = snp->snpNamDict ;

      aceOutDate (ao, "//", "Phasing") ;

      for (ii = 0 ; ii < iMax ; ii++)
	{
	  AES *ap = arrayp (raw, ii, AES) ;
	  int phase = ap->phase ;
	  int strand = 1 ;
	  int d ;
	  int nam = ap->vcfNam ; /* sould be aceNam */

	  d = ap->phasingScore[3] ;
	  if (phase)
	    {
	      if (phase < 0)
		{ strand = 2 ; phase = -phase ; }
	      
	      aceOutf (ao, "Variant \"%s\"\nPhase %d \"%s\" %d\n\n"
		       , snpNamDict && nam ? dictName (snpNamDict, nam) : dictName(dict, ii)
		       , strand
		       , dictName (phaseDict, phase)
		       , ap->phasingScore[6] 
		       ) ;
	    }
	  else if (d)
	    aceOutf (ao, "Variant \"%s\"\nPhase %d\n\n"
		     , snpNamDict && nam ? dictName (snpNamDict, nam) : dictName(dict, ii)
		     , -d
		     ) ;

	}

      ac_free (ao) ;
    }

  ac_free (h) ;
  return 1 ;
} /* snpPhasingExport */

/*************************************************************************************/
/* plug the list of mismatches, known from the snp file, in the gg array */
static void snpPhasingInterpretOneMissmatch (SNP *snp, HIT *up, Array gg, Array raw, char *cp)
{
  char type[8] ;
  DICT *dict = snp->snpDict ;
  int s ;
  char *cq ;
  
  cq = strchr (cp, ':') ;
  if (cq > cp)
    {
      int a1 = 0 ;

      *cq = 0 ;
      sscanf (cp, "%d", &a1) ;
      cp = cq + 1 ;
      if (a1 > 0)
	{
	  char buf[256] ;
	  if (cp[0] == '*') cp++ ;
	  strncpy (type, cp, 7) ;
	  type[7] = 0 ; 

	  sprintf (buf, "%s:%d:%s", dictName(snp->targetDict,up->target), a1, type) ;
	  if (dictFind (dict, buf, &s))
	    {
	      AES *ap = arrayp (raw, s, AES) ; 
	      PHS *phs ;
	      int ig ;
	      BOOL ok = FALSE ;

	      for (ig = 0 ; ! ok &&  ig < arrayMax (gg) ; ig++)
		{
		  phs = arrp (gg, ig, PHS) ;
		  if (ap->a1 == phs->a1)
		    ok = TRUE ; 
		}
	      if (! ok)
		{
		  phs = arrayp (gg, arrayMax (gg), PHS) ;
		  phs->snp1 = s ;
		  phs->a1 = ap->a1 ;
		  phs->m1 = up->mult ; phs->w1 = 0 ;
		  ap->m += up->mult ;
		}
	    }
	}
    }
  return ;
} /* snpPhasingInterpretOneMissmatch  */

/*************************************************************************************/
/* register all the snps seen in that HIT (easy) and all the wild type covered (less easy) */
static BOOL snpPhasingInterpretOneHit (SNP *snp, HIT *up, Array gg, Array raw)
{
  int a1 = up->a1, a2 = up->a2 ;
  AES *aes ;
  int target = up->target ;
  static int iSnp = 0 ;
  Stack s = snp->aliExtendStack ;
  Array snps = snp->aliSnps ;
  int iSnpMax = arrayMax (snps) ;
    
  if (a1 > a2)
    { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }
  if (up->errTarget)
    {
      char *cp, *cq ;

      cp = stackText (s, up->errTarget) ;
      while ((cq = strchr (cp, ',')))
	{
	  *cq = 0 ;
	  snpPhasingInterpretOneMissmatch (snp, up, gg, raw, cp) ;
	  cp = cq + 1 ;
	}
      snpPhasingInterpretOneMissmatch (snp, up, gg, raw, cp) ;
    }
  /* position back */
  if (iSnp) iSnp-- ;
  aes = arrayp (snps, iSnp, AES) ;
  while (iSnp > 0 && aes->target >= target && aes->a1 >= a1)
    { aes-- ; iSnp-- ; }
  /* scan for releant snps */
  while (iSnp < iSnpMax && aes->target < target)
    { aes++ ; iSnp++ ; }
  while (iSnp < iSnpMax && aes->target == target && aes->a1 <= a2 - 8)
    {
      if ( aes->a1 >= a1 + 8 )
	{
	  int ig ;
	  PHS *phs ;
	  BOOL ok = FALSE ;
	  /* check if we have seen that snp in our fragment */
	  for (ig = 0 ; ! ok &&  ig < arrayMax (gg) ; ig++)
	    {
	      phs = arrp (gg, ig, PHS) ;
	      if (aes->a1 == phs->a1)
		ok = TRUE ; 
	    }
	  if (!ok)
	    {
	      AES *ap = arrayp (raw, aes->snp, AES) ; 
	      phs = arrayp (gg, arrayMax (gg), PHS) ;
	      
	      phs->snp1 = aes->snp ;
	      phs->a1 = aes->a1 ;
	      phs->w1 = up->mult ; phs->m1 = 0 ;
	      ap->w += up->mult ;
	    }
	}
      aes++ ; iSnp++ ;
    }
  return TRUE ;
} /* snpPhasingInterpretOneHit */

/*************************************************************************************/

static void snpPhasing (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  Array hits = arrayHandleCreate (32, HIT, h) ;
  Array gg = arrayHandleCreate (32, PHS, h) ;
  Array raw = arrayHandleCreate (100000, AES, h) ;
  BigArray ph = bigArrayHandleCreate (100000, PHS, h) ;
  ACEIN ai = aceInCreate (snp->inFileName, snp->gzi, h) ;
  int kk = 0, runQualityPrefix = 0 ;
  long int iph = 0 ;

  snp->aliExtendStack = stackHandleCreate (100000, snp->h) ;
  snpAliExtendParseSnpList (snp, raw) ;  /* fills array: snp->aliSnps */
  arraySort (snp->aliSnps, aesOrder) ; 

  /* get all hits for a given read, we should do it per clone */
  while (snpAliExtendGetHits (snp, ai, hits, &runQualityPrefix))
    {
      int ii, jj ;
      HIT *up ;
      PHS *ap, *bp ;
      /* for all hits of the same clone, check if we cover as wild type or mutant the given snps */
      gg = arrayReCreate (gg, 32, PHS) ;
      for (ii = 0, up = arrp (hits, ii, HIT) ; ii < arrayMax (hits) ; up++, ii++)
	snpPhasingInterpretOneHit (snp, up, gg, raw) ;
      /* register all the phasing info */
      arraySort (gg, phsOrder) ;

      for (ii = 0, ap = arrp (gg, ii, PHS) ; ii < arrayMax (gg) ; ap++, ii++)
	{
	  if (ii && ap->snp1 == (ap-1)->snp1)
	    continue ;
	  for (jj = ii + 1, bp = ap + 1 ; jj < arrayMax (gg) ; bp++, jj++)
	    {
	      PHS *xp, *yp ;
	      PHS *wp ;
	      
	      if (bp->snp1 == (bp-1)->snp1)
		continue ;
	      
	      wp = bigArrayp (ph, iph++, PHS) ;
	      kk++ ;
	      if (ap->snp1 < bp->snp1)
		{
		  xp = ap ; yp = bp ;
		}
	      else
		{
		  xp = bp ; yp = ap ;
		}
	      wp->snp1 = xp->snp1 ;
	      wp->snp2 = yp->snp1 ;
	      if (xp->m1 && yp->m1)
		wp->mm += xp->m1 ;
	      else if (xp->w1 && yp->w1)
		wp->ww += xp->w1 ;
	      else if (xp->w1 && yp->m1)
		wp->wm += xp->w1 ;
	      else if (xp->m1 && yp->w1)
		wp->mw += xp->m1 ;
	      wp->m1 += xp->m1 ;
	      wp->w1 += xp->w1 ;
	      wp->m2 += yp->m1 ;
	      wp->w2 += yp->w1 ;
	    }
	  if (0 && kk > 1000000)
	    { iph = snpPhasingCompress (ph) ; kk = 0 ; }
	}
    }
  snpPhasingCompress (ph) ;
  snpPhasingExport (snp, ph, raw) ;

  ac_free (h) ;
} /* snpPhasing */

/*************************************************************************************/
/*************************************************************************************/
typedef struct aehStruct { int target, score, chrom1, pos1, chrom2, pos2, cover, mp, mm, wp, wm
			     , bad, oo1p, oo1m, oo2p, oo2m, multip, multim
			     , type, typeW, typeM, ss, genotype
			     ; } AEH ;

/*************************************************************************************/

static int aehOrder (const void *va, const void *vb)
{
  const AEH *a = (const AEH *)va, *b = (const AEH *)vb ;
  int n ;

  n = a->chrom1 - b->chrom1 ; if (n) return n ;
  n = a->pos1 - b->pos1 ; if (n) return n ;
  n = a->chrom2 - b->chrom2 ; if (n) return n ;
  n = a->pos2 - b->pos2 ; if (n) return n ;
  n = a->type - b->type ; if (n) return n ;

  return 0 ;
} /* aehOrder */

/*************************************************************************************/
/* we need a delayed registration because a probe matching both
 * the reference and the mutant should not count in either
 *
 * In particular this avoid the problem of counting the support for indels
 * reads that do not read the base following the indel will match wild and mutant 
 * and be eliminated
 */
static void snpAnalyzeEditedHitsRegisterPreviousProbe (Array aa, Array bb, DICT *phasingDict, Array phasing, int delta, BOOL isMulti)
{
  int ii, jj, kill,  mult ;
  AEH *up, *vp ;
  PHS *phs ;

  /* for each target, we check if we support both the mutant and the wild type 
   * at the same score, if so we kill both
   * else we kill the lower score
   */
  for (ii = 0, up = arrp (bb, 0, AEH) ; ii < arrayMax (bb) ; up++, ii++)
    {
      kill = 0 ;
      if (up->target)
	{
	  for (jj = ii + 1, vp = up + 1 ; jj < arrayMax (bb) ; vp++, jj++)
	    if (up->target == vp->target)
	      {
		if (up->score == vp->score)
		  { kill = 1 ; vp->target = 0 ; }
		if (up->score > vp->score)
		  vp->target = 0 ;
		if (up->score < vp->score)
		  kill = 1 ;
	      }
	}
      if (kill) up->target = 0 ;
    }
	
  /* we check also if we support several distinct positions and do the phasing */
  if (! isMulti)
    {
      AEH *uu, *vv ;
      int nph = 0 ;

      for (ii = 0, up = arrp (bb, 0, AEH) ; ii < arrayMax (bb) ; up++, ii++)
	{
	  if (up->target)
	    {
	      for (jj = ii + 1, vp = up + 1 ; jj < arrayMax (bb) ; vp++, jj++)
		if (vp->target)
		  {
		    if (up->chrom1 !=  vp->chrom1)
		      continue ;
		    if (up->pos1 == vp->pos1)
		      continue ;
		    if (up->pos1 < vp->pos1)
		      { uu = up ; vv = vp ; }
		    else
		      { uu = vp ; vv = up ; }
		    if (uu->pos1 + delta < vv->pos1)
		      continue ;
		    dictAdd (phasingDict, messprintf ("%d:%d", uu->target, vv->target), &nph) ;
		    phs = arrayp (phasing, nph, PHS) ;
		    phs->snp1 = uu->target ;
		    phs->snp2 = vv->target ;
		    mult = uu->mp + uu->mm + uu->wp + uu->wm ;

		    if (uu->mp + uu->mm && vv->mp + vv->mm) phs->mm += mult ;
		    else if (uu->wp + uu->wm  &&  vv->wp + vv->wm) phs->ww += mult ;
		    else if (uu->wp + uu->wm  &&  vv->mp + vv->mm) phs->wm += mult ;
		    else if (uu->mp + uu->mm  &&  vv->wp + vv->wm) phs->mw += mult ;
		  }
	    }
	}
    }

  for (ii = 0 ; ii < arrayMax (bb) ; ii++)
    {
      up = arrp (bb, ii, AEH) ;
      if (! up->target)
	continue ;

      vp = arrayp (aa, up->target, AEH) ;
      if (isMulti)
	{
	  vp->multip += up->mp + up->oo1p + up->oo2p ;
	  vp->multim += up->mm + up->oo1m + up->oo2m ;
	}
      else
	{
	  /* we integrate the previous count */
	  up->mp += vp->mp ;
	  up->mm += vp->mm ;
	  up->wp += vp->wp ;
	  up->wm += vp->wm ;
	  
	  up->bad += vp->bad ;
	  up->oo1p += vp->oo1p ;
	  up->oo2p += vp->oo2p ;
	  up->oo1m += vp->oo1m ;
	  up->oo2m += vp->oo2m ;
	  
	  /* we clean up the hacked tags */
	  up->ss = up->genotype = 0 ;
	  
	  /* the we copy back the whole record */
	  *vp = *up ;
	}
    }
  bb = arrayReCreate (bb, 100, AEH) ;
} /* snpAnalyzeEditedHitsRegisterPreviousProbe */

/*************************************************************************************/

static int snpAnalyzeEditedHits (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, nok = 0 ;
  int delta =  snp->analyzeEditedHits ; /* central position */
  DICT *targetDict = dictHandleCreate (1000000, h) ;
  DICT *chromDict = dictHandleCreate (100, h) ;
  ACEIN ai = snp->ai ;
  ACEOUT ao = 0 ;
  char probeBuf[1000], pErr[1000], tErr[1000], wm ;
  const char *ccp ;
  char *cp, targetBuf[1000] ;
  int score, mult, toBeAligned, ali, x1, x2, tmult, a1, a2, nN, nErr, target, chrom1, chrom2, pos1, pos2 ;
  int type, typeW, typeM, oo1, oo2 ;
  AEH *up ;
  BOOL isDown, isSolid = FALSE ;
  Array aa = arrayHandleCreate (100000, AEH, h) ;
  Array bb = arrayHandleCreate (100, AEH, h) ;
  Array phasing = arrayHandleCreate (100, PHS, h) ;
  DICT *phasingDict = dictHandleCreate (100, h) ;
  BOOL isMulti = FALSE ;

  memset (probeBuf, 0, sizeof(probeBuf)) ;
  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ; if (!ccp) continue ;
      if (*ccp == '#') continue ;
      if (! strcmp (ccp, "zzzzz"))
	{ isMulti = TRUE ; continue ; }
      if (probeBuf[0] && strcmp (probeBuf, ccp) && arrayMax(bb))
	snpAnalyzeEditedHitsRegisterPreviousProbe (aa, bb, phasingDict, phasing, delta, isMulti) ;
      strncpy (probeBuf, ccp, 999) ;
      aceInStep (ai, '\t') ;  aceInInt (ai, &score) ;
      aceInStep (ai, '\t') ;  aceInInt (ai, &mult) ;
      aceInStep (ai, '\t') ;  aceInInt (ai, &toBeAligned) ;
      aceInStep (ai, '\t') ;  aceInInt (ai, &ali) ;
      if (ali < 140 && 100 * ali < snp->minAliPerCent * toBeAligned)
	continue ;
      aceInStep (ai, '\t') ;  aceInInt (ai, &x1) ;
      aceInStep (ai, '\t') ;  aceInInt (ai, &x2) ;
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; if (!ccp || strcmp (ccp, "s_snp")) continue ; 
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; if (!ccp) continue ;  /* gene */
      aceInStep (ai, '\t') ;  aceInInt (ai, &tmult) ;
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; if (!ccp) continue ; 
      wm = *ccp ; 
      dictAdd (targetDict, ccp+2, &target) ;

      strncpy (targetBuf,  dictName (targetDict, target), 999) ;

      for (cp = targetBuf ; *cp != ':'  && *cp ; cp++) ;
      *cp++ = 0 ;
      dictAdd (chromDict, targetBuf, &chrom1) ;

      for (pos1 = 0 ; *cp != ':'  && *cp ; cp++)
	pos1 = 10 * pos1 + (*cp - '0') ;
      *cp++ = 0 ;

      for (ccp = cp ; *cp != ':'  && *cp ; cp++) ;
      *cp++ = 0 ;
      dictAdd (chromDict, ccp, &chrom2) ;

      for (pos2 = 0 ; *cp != ':'  && *cp ; cp++)
	pos2 = 10 * pos2 + (*cp - '0') ;
      *cp++ = 0 ;

      for (ccp = cp ; *cp != ':'  && *cp ; cp++) ;
      *cp++ = 0 ;
      dictAdd (chromDict, ccp, &type) ;

      for (ccp = cp ; *cp != '.'  && *cp ; cp++) ;
      *cp++ = 0 ;
      dictAdd (chromDict, ccp, &typeW) ;

      for (ccp = cp ; *cp != ':'  && *cp ; cp++) ;
      *cp++ = 0 ;
      dictAdd (chromDict, ccp, &typeM) ;

      if (! strncmp (dictName (chromDict, type), "SN", 2))
	{
	  cp = messprintf("%c>%c", (dictName (chromDict, typeW))[1], (dictName (chromDict, typeM))[1]) ;
	  dictAdd (chromDict, cp, &type) ;
	}

      aceInStep (ai, '\t') ;  aceInInt (ai, &a1) ;
      aceInStep (ai, '\t') ;  aceInInt (ai, &a2) ;
      if (a1 < a2) { isDown = TRUE ; }
      else { int a0 = a1 ; a1 = a2 ; a2 = a0 ; isDown = FALSE ; }
      if (a1 > delta - 1 || a2 < delta + 1) continue ;
      aceInStep (ai, '\t') ;  aceInInt (ai, &nN) ;
      aceInStep (ai, '\t') ;  aceInInt (ai, &nErr) ;
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; if (!ccp) continue ; strncpy (pErr, ccp, 999) ;
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; if (!ccp) continue ; strncpy (tErr, ccp, 999) ;

      /************* count the oo ********/
      
      oo1 = oo2 = 0 ;
      if (strstr (tErr, "oo"))
	{
	  char *cp1 = tErr, *cp2 ;
	  int x ;

	  isSolid = TRUE ;
	  while (*cp1)
	    {    /* parse    2345:ag>oo */
	      for (cp2 = cp1 ; *cp2 != 'o' && *cp2 ; cp2++) ;
	      if (*cp2 == 'o' && *(cp2 - 4) == ':')
		{
		  cp2 -= 5 ;
		  while (*cp2 >= '0' && *cp2 <= '9') cp2-- ;
		  for (x = 0, cp2++ ; *cp2 != ':' ; cp2++)
		    x = 10 * x + *cp2 - '0' ;
	
		  if (x == delta - 1) /* the transition to the left of the snp is dubious */
		    oo1 += mult ;
		  if (x == delta)     /* the transition to the right of the snp is dubious */
		    oo2 += mult ;

		  for (cp1 = cp2 + 1 ; *cp1 != ':' && *cp1 ; cp1++) ; /* find the next error position */
		}
	      else
		break ;
	    }
	}
      
      ii = isDown ? 0 : 1 ;
      switch (wm)
	{
	case 'M': ii += 0 ; break ;
	case 'W': ii += 2 ; break ;
	default: continue ;
	}

      up = arrayp (bb, arrayMax(bb), AEH) ;
      up->target = target ;
      up->score = score ;
      up->chrom1 = chrom1 ;
      up->chrom2 = chrom2 ;
      up->pos1 = pos1 ;
      up->pos2 = pos2 ;
      up->type = type ;
      up->typeW = typeW ;
      up->typeM = typeM ;

      /* ATTENTION: if the count types are edited one must 
	 must synchronize
	 snpAnalyzeEditedHitsRegisterPreviousProbe 
      */
      if (oo1 + oo2 == 0)
	{
	  switch (ii)
	    {
	    case 0: up->mp += mult ; break ;
	    case 1: up->mm += mult ; break ;
	    case 2: up->wp += mult ; break ;
	    case 3: up->wm += mult ; break ;
	    }
	}
      else
	{
	  up->bad = 1 ;
	  switch (ii)
	    {
	      /* forget the oo on the mutant side thay should never occur */
	    case 2: up->ss = 1 ; up->oo1p += oo1 ; up->oo2p += oo2 ; break ;
	    case 3: up->ss = 1 ; up->oo1m += oo1 ; up->oo2m += oo2 ; break ;
	    }
	}
    }
  
  if (bb && arrayMax(bb))
    snpAnalyzeEditedHitsRegisterPreviousProbe (aa, bb, phasingDict, phasing, delta, isMulti) ;
 
  if (! arrayMax(aa))
    messcrash ("No data read from the input file, sorry") ;
  /* sort aa must come AFTER exporting the phasing */
  arraySort (aa, aehOrder) ;
  /* duplicate the wild type corresponding to the same positions */
  for (ii = 0, up = arrp (aa, 0, AEH) ; ii < arrayMax (aa) ; ii++, up++)
    {
      int jj, wp = 0, wm = 0, cover = 0 ;
      AEH *vp ;

      chrom1 = up->chrom1 ; chrom2 = up->chrom2 ;
      pos1 = up->pos1 ; pos2 = up->pos2 ;

      for (jj = ii, vp = up ; jj < arrayMax (aa) && vp->chrom1 == chrom1 && vp->chrom2 == chrom2 && vp->pos1 == pos1 && vp->pos2 == pos2 ; jj++, vp++)
	{ wp += vp->wp ; wm += vp->wm ; cover += vp->mp + vp->mm + vp->wp + vp->wm ; }
      for (jj = ii, vp = up ; jj < arrayMax (aa) && vp->chrom1 == chrom1 && vp->chrom2 == chrom2 && vp->pos1 == pos1 && vp->pos2 == pos2 ; jj++, vp++)
	{ vp->wp = wp ; vp->wm = wm ; vp->cover = cover ; }
      ii = jj - 1 ; up = vp - 1 ;
    }
  
  /* analyze the genotypes */
  for (ii = 0, up = arrp (aa, 0, AEH) ; ii < arrayMax (aa) ; ii++, up++)
    {
      int chiFlag = 0, notLow = 0 ;
      double z = 0, z1 = 0, z2 = 0 ;
      char *ss = "-" ;

      /*
	int m, w ;
	m = up->mp + up->mm ; 
	w = up->wp + up->wm ; 
      */

      if (! up->cover || ! up->target) continue ;

      /* count the baddy_too_many_hits which align on the same spot */
      /* deal with the tags that map in several places */
      /* deal with the phasing */
      /* size of the pseudo exons */
      /* y a t-il encore des overhangs locaux alors que les mutants sont integres */

      notLow = snpCompatibleStrands (up->cover, up->mp, up->wp, up->mm, up->wm, &z, &z1, &z2, &chiFlag, snp->minCover) ;
      ss = cp = "-" ;
      switch (chiFlag)
	{
	case 0:
	  ss = 0 ;
	  if (2 * (up->mp + up->wp) < snp->minCover) ss = "Minus_strand_only" ;
	  else if (2 * (up->mm + up->wm) < snp->minCover) ss = "Plus_strand_only" ;
	  break ;
	case 0x1:
	  ss = "Both_strands" ;
	  break ;
	case 0x2:
	  ss = "Incompatible_strands_A_5_100" ;
	  break ;
	case 0x4:
	  ss = "Incompatible_strands_B_1_100" ;
	  break ;
	case 0x8:
	  ss = "Incompatible_strands_C_1_1000" ;
	  break ;
	}
      dictAdd (chromDict, ss, &(up->ss)) ;
      
      if ((notLow & 0x01) && (notLow & 0x02) && (notLow & 0x04))
	cp = z < 40 ? "intermediate" : "m/?" ;
      else if ((notLow & 0x02) && (notLow & 0x04)) 
	cp = "+/+" ;
      else if ((notLow & 0x01) && (notLow & 0x04)) 
	cp = "m/+" ;
      else if ((notLow & 0x01) && (notLow & 0x02)) 
	cp = "m/m" ;
      else if (notLow & 0x01) 
	cp = "m/?" ;
      else if (notLow & 0x02) 
	cp = z < 40 ? "intermediate" : "m/?" ;
      else if (notLow & 0x04) 
	cp = "+/?" ;
      else
	cp = "-" ;
      
      dictAdd (chromDict, cp, &(up->genotype)) ;
    }

  /* export the results */
  ao = aceOutCreate (snp->outFileName, ".countEdited", snp->gzo, h) ;
  aceOutDate (ao, "###", snp->project) ;
  aceOutf (ao, "## Number of reads supporting the SNP or the wild type, per strand\n") ;
  aceOutf (ao, "## The position is one based and indicate the first modified base, [mw]+- refers to strand +/-\n") ;
  aceOutf (ao, "## In case of break points, columns 4 and 5 refer to the other half of the breakpointThe position is one based and indicate the first modified base, [mw]+- refers to strand +/-\n") ;

  aceOutf (ao, "# Chromosome\tPosition\tType\tReference or second chrom\tVariant or second position\tRun\tDosage\t%% mutant\tCover\tMutant\tWild\tm+\tw+\tm-\tw-\tN\tN+\tN-\tcompatibility\tchi2\t-\tmultiply mapped+\tmultiply mapped-\t-\n") ;
  for (ii = 0, up = arrp (aa, 0, AEH) ; ii < arrayMax (aa) ; ii++, up++)
    {
      int m = up->mp + up->mm ; 
      int w = up->wp + up->wm ; 
      int notLow, chiFlag = 0 ;
      char *ss, *cp ;
      const char *ccp ;
      float dose ;
      double z, z1, z2 ;

      if (up->cover && up->target)
	{
	  /* dosage */
	  /* chi2 */
	  notLow = snpCompatibleStrands (up->cover, up->mp, up->wp, up->mm, up->wm, &z, &z1, &z2, &chiFlag, snp->minCover) ;
	  ss = cp = "-" ;
	  switch (chiFlag)
	    {
	    case 0:
	      if (2 * (up->mp + up->wp) < snp->minCover) ss = "Minus_strand_only" ;
	      else if (2 * (up->mm + up->wm) < snp->minCover) ss = "Plus_strand_only" ;
	      break ;
	    case 0x1:
	      ss = "Both_strands" ;
	      break ;
	    case 0x2:
	      ss = "Incompatible_strands_A_5_100" ;
	      break ;
	    case 0x4:
	      ss = "Incompatible_strands_B_1_100" ;
	      break ;
	    case 0x8:
	      ss = "Incompatible_strands_C_1_1000" ;
	      break ;
	    }
	  if ((notLow & 0x01) && (notLow & 0x02) && (notLow & 0x04))
	    { dose = z < 40 ? .5 : 1.5 ; cp = z < 40 ? "intermediate" : "m/?" ; }
	  else if ((notLow & 0x02) && (notLow & 0x04)) 
	    { dose = 0 ; cp = "+/+" ; }
	  else if ((notLow & 0x01) && (notLow & 0x04)) 
	    { dose = 1 ; cp = "m/+" ; }
	  else if ((notLow & 0x01) && (notLow & 0x02)) 
	    { dose = 2 ; cp = "m/m" ; }
	  else if (notLow & 0x01) 
	    { dose = 1 ; cp = "m/?" ; }
	  else if (notLow & 0x02) 
	    { dose =  z < 40 ? .5 : 1.5 ; cp = z < 40 ? "intermediate" : "m/?" ; }
	  else if (notLow & 0x04) 
	    { dose = .5 ; cp = "+/?" ; }
	  else
	    { dose = -1 ; cp = "-" ; }
	  
	  if (up->type)
	    {
	      ccp = dictName (chromDict, up->type) ;
	      aceOutf (ao, "%s\t%d\t%s%s%s"
		       , dictName (chromDict, up->chrom1)
		       , up->pos1 + ( ccp[1] == '>' || !strncmp (ccp, "Del", 3) || !strncmp (ccp, "Ins", 3) ? 1 : 0)
		       , ccp
		       , !strncmp (ccp, "Del", 3) &&  up->typeW ? dictName (chromDict, up->typeW)  : ""
		       , !strncmp (ccp, "Ins", 3) && up->typeM ? dictName (chromDict, up->typeM)  : ""
		       ) ;

	      if (ccp &&
		  ( ccp[1] == '>' || !strncmp (ccp, "Del", 3) || !strncmp (ccp, "Ins", 3))
		  )
		{
		  /* export the snipnet */
		  aceOutf (ao, "\t%s\t%s"
			 , up->typeW ? dictName (chromDict, up->typeW) : "-"
			 , up->typeM ? dictName (chromDict, up->typeM) : "-" 
			 ) ;
		}
	      else  /* BRK point export the coordinate on the other chrom */
		{
		  aceOutf (ao, "\t%s\t%d"
			   , dictName (chromDict, up->chrom2)
			   , up->pos2
			   ) ;
		}
	      
	      aceOutf (ao, "\t%s", snp->run) ;
	      if (dose == 0) aceOutf (ao, "\t0") ;
	      else if (dose == 1) aceOutf (ao, "\t1") ;
	      else if (dose == 2) aceOutf (ao, "\t2") ;
	      else if (dose == -1) aceOutf (ao, "\t-") ;
	      else aceOutf (ao, "\t%.1f", dose) ;

	      aceOutf (ao, "\t%.1f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d"
		       , z
		       , up->cover, m, w, up->mp, up->wp, up->mm, up->wm
		       , up->bad
		       ) ;

	      if (isSolid)
		aceOutf (ao, "\t%d:%d\t%d:%d"
			 , up->oo1p, up->oo2p, up->oo1m, up->oo2m
			 ) ;
	      else
		aceOutf (ao, "\t%d\t%d"
			 , up->oo1p + up->oo2p, up->oo1m, up->oo2m
			 ) ;

	      aceOutf (ao, "\t%s\t%.1f\t%s\t%d\t%d\t%s"
		       , ss, z2
		       , 100 * (up->oo1p + up->oo1m + up->oo2p + up->oo2m) > 40 * up->cover ? "OO" : "-"
		       , up->multip, up->multim
		       , 100 * (up->multip + up->multim) > 40 * up->cover ? "REPEATED" : "-"
		       ) ;
	    }
	  aceOutf (ao, "\n") ;
	}
    }
  ac_free (ao) ;


  /* export the statistics */
  if (1)
    {
      int tMut = 0, tSnp = 0 ;
      int ntransitions = 0, ntransversions = 0 ;
      int tntransitions = 0, tntransversions = 0 ;
      int ns = 0, ni = 0, nd = 0, no = 0, u ;
      int tns = 0, tni = 0, tnd = 0, tno = 0 ;
      int _Ins, _Del, _Other, i, ii ;
      const char *typ ;
      char buf[1024], **typep ;
      KEYSET ks = keySetHandleCreate (h) ;
      KEYSET mks = keySetHandleCreate (h) ;
      DICT *tDict = dictHandleCreate (100, h) ;
      char *types[] = { "A>G", "T>C", "G>A", "C>T"
			, "A>T", "T>A", "G>C", "C>G"
			, "A>C", "T>G", "G>T", "C>A"
			, "Ins A", "Ins T", "Ins G", "Ins C"
			, "Del A", "Del T", "Del G", "Del C"
			, 0  } ;

      dictAdd (chromDict, "Ins", &_Ins) ;
      dictAdd (chromDict, "Del", &_Del) ;
      for (typep = types ; *typep ; typep++)
	dictAdd (tDict, *typep, &ii) ;
      dictAdd (tDict, "Other", &_Other) ;
      ao = aceOutCreate (snp->outFileName, ".stats.txt", snp->gzo, h) ;

      for (ii = 0, up = arrp (aa, 0, AEH) ; ii < arrayMax (aa) ; ii++, up++)
	{
	  int m = up->mp + up->mm ; 
	  int w = up->wp + up->wm ; 
	  int cover = m + w ;
	  
	  if (!(cover >= 4 && up->target && 100.0 * m >= 20 * cover && m >= 3))
	    continue ;
	  
	  tMut += m ;
	  tSnp++ ;
	
	  typ = dictName (chromDict, up->type) ;
	  if(up->type == _Del)
	    {
	      nd++ ; tnd += m ; 
	      sprintf (buf, "Del %s", dictName (chromDict, up->typeW)) ;
	      dictAdd (tDict, buf, &i) ;
	      keySet (ks, i)++ ;
	      keySet (mks, i) += m ;
	    }
	  else if(up->type == _Ins)
	    {
	      ni++ ; tni += m ; 
	      sprintf (buf, "Ins %s", dictName (chromDict, up->typeM)) ;
	      dictAdd (tDict, buf, &i) ;
	      keySet (ks, i)++ ;
	      keySet (mks, i) += m ;
	    }
	  else if (typ[1] == '>')
	    {
	      ns++ ; tns += m ;
	      dictAdd (tDict, dictName (chromDict, up->type), &i) ;
	      keySet (ks, i)++ ;
	      keySet (mks, i) += m ;
	      if (i <= 4)
		{
		  ntransitions++ ;
		  tntransitions += m ;		  
		}
	      else
		{
		  ntransversions++ ;
		  tntransversions += m ;		  
		}
	    }
	  else
	    {
	      ns++ ; tns += m ;
	      keySet (ks, _Other)++ ;
	      keySet (mks, _Other) += m ;
	    }
	}

      aceOutDate (ao, "###", snp->project) ;
      aceOutf (ao, "## Number of variants validated with frequency above 20%% and coverage 4 or more\n") ;
      aceOutf (ao, "# Run\tAny variants\tSubstitutions\tTransitions\tTransversions\tInsertions\tDeletions\tOther") ;
      for (typep = types ; *typep ; typep++)
	aceOutf (ao, "\t%s", *typep) ;

      aceOutf (ao, "\t\t\tRun\tReads supporting the variants\tSubstitutions\tTransitions\tTransversions\tInsertions\tDeletions\tOther") ;
      for (typep = types ; *typep ; typep++)
	aceOutf (ao, "\t%s", *typep) ;

      aceOutf (ao, "\t\t\tRun\tPrevalence of variants\tSubstitutions\tTransitions\tTransversions\tInsertions\tDeletions\tOther") ;
      for (typep = types ; *typep ; typep++)
	aceOutf (ao, "\t%s", *typep) ;

      aceOutf (ao, "\t\t\tRun\tAverage number of reads supporting each type of variant\tSubstitutions\tTransitions\tTransversions\tInsertions\tDeletions\tOther") ;
      for (typep = types ; *typep ; typep++)
	aceOutf (ao, "\t%s", *typep) ;


      /* counts */
      aceOutf(ao, "\n%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d", snp->run, tSnp, ns, ntransitions, ntransversions, ni, nd, no) ;
      for (typep = types ; *typep ; typep++)
	{
	  dictFind (tDict, *typep, &ii) ;
	  aceOutf(ao, "\t%d", keySet (ks, ii)) ;
	}

      /* support */
      aceOutf(ao, "\t\t\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d", snp->run, tMut, tns, tntransitions, tntransversions, tni, tnd, tno) ;
      for (typep = types ; *typep ; typep++)
	{
	  dictFind (tDict, *typep, &ii) ;
	  aceOutf(ao, "\t%d", keySet (mks, ii)) ;
	}

      /* percent of the type of SNP among all SNPs */
      aceOutf(ao, "\t\t\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f", snp->run, 100.0*tSnp/tSnp, 100.0*ns/tSnp, 100.0*tntransitions/tSnp, 100.0*tntransversions/tSnp, 100.0*ni/tSnp, 100.0*nd/tSnp, 100.0*no/tSnp);
      for (typep = types ; *typep ; typep++)
	{
	  dictFind (tDict, *typep, &ii) ;
	  aceOutf(ao, "\t%.2f", 100.0*keySet (mks, ii)/(double)tSnp) ;
	}

      /* average support per type of SNP */
      if (tSnp == 0) tSnp = 1 ; 
      if (ns == 0) ns = 1 ; 
      if (ni == 0) ni = 1 ; 
      if (nd == 0) nd = 1 ;
      aceOutf(ao, "\t\t\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", snp->run, tSnp/(double)tSnp, tns/(double)ns,  tntransitions/(double)ns,  tntransversions/(double)ns, tni/(double)ni, tnd/(double)nd);
      for (typep = types ; *typep ; typep++)
	{
	  dictFind (tDict, *typep, &ii) ;
	  u = keySet (ks, ii) ; if (u == 0) u = 1 ;
	  aceOutf(ao, "\t%.2f", keySet (mks, ii)/(double)u) ;
	}
      
      aceOutf(ao, "\n") ;
    }

  ac_free (h) ;

  return nok ;
} /* snpAnalyzeEditedHits */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/* parse an acedb table with 8 columns generated by  */
static int snpDeepTableReport (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0, ii = 0, jj, strand = 0, nLines = 0 ;
  int run, variant ;
  int allStrands = 0 ;
  ACEIN ai = snp->ai ;
  ACEOUT ao = snp->ao ;
  Array snps = arrayHandleCreate (10000, SSM, h) ;
  Array cumulBelow = arrayHandleCreate (1000, int, h) ;
  Array cumulOver = arrayHandleCreate (1000, int, h) ;
  Array cumulMiddle = arrayHandleCreate (1000, int, h) ;
  Array cumulAbsent = arrayHandleCreate (1000, int, h) ;
  Array cumulTotal = arrayHandleCreate (1000, int, h) ;
  DICT *runDict = dictHandleCreate (200,h) ;
  DICT *eeDict = eeDictCreate (snp, h) ;
  DICT *vDict = dictHandleCreate (100000, h) ;
  BOOL ok = FALSE ;
  SSM *ssm, *ssm2 ;
  float z ;
  char bufSnp[100] ;
  const char *ccp, *ccq ;

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ; /* Variant name, use to group the results */
      if (!ccp || *ccp == '#')
	continue ;
      if ( *ccp == '/' &&  *(ccp+1)== '/')
	continue ;

      dictAdd (vDict, ccp, &variant) ;
      /* actual snp */
      ssm = arrayp (snps, arrayMax (snps), SSM) ;
      ssm->variant = variant ;

      ccp = aceInWord (ai) ; /* mRNA */
      dictAdd (snp->targetDict, ccp, &(ssm->target)) ;
      
      aceInInt (ai, &(ssm->pos)) ; /* position */
	  
      ccp = aceInWord (ai) ;  /* type */
      strcpy (bufSnp, ccp) ;
      if (bufSnp[1] == '2' && bufSnp[3] == 0) bufSnp[1] = '>' ;
      dictAdd (eeDict, bufSnp, &(ssm->snp)) ;

      ccp = aceInWord (ai) ; /* parent run  */
      dictAdd (runDict, ccp, &(ssm->run)) ;

      ccq = dictName(runDict, ssm->run) ;
      ccq += strlen(ccq) - 2 ;
      ccp = aceInWord (ai) ; /* parent experiment stranding   */
      (ssm->snp)++ ;
      if (!strcasecmp(ccq, ".r") ||  !strcasecmp (ccp,"Reverse")) { allStrands |= 0x2 ; ssm->snp *= -1 ; }
      else { allStrands |= 0x1  ; }

      if (! aceInInt (ai, &(ssm->count)))
	messcrash ("Expecting a count in column 7, line %d, empty column, in file %s"
		   , aceInStreamLine (ai), aceInFileName (ai)) ;
      if (! aceInInt (ai, &(ssm->coverS)))
	messcrash ("Expecting coverage integer in column 8, line %d, empty column, in file %s"
		   , aceInStreamLine (ai), aceInFileName (ai)) ;
      nLines++ ;
    } /* aceInCard */

  if (1)
    {
      aceOutDate (ao, "###", snp->title) ;
      aceOutf (ao, "# SNPs across %d runs at threshold %d%% and minimal count %d and minimal coverage %d\n"
	       , dictMax (runDict)
	       , snp->minFrequency
	       , snp->minMutant
	       , snp->minCover
	       ) ;
      aceOutf (ao, "#Target") ;
      if (snp->target2dna)   
	aceOutf (ao, "\tsequence") ;
      if (snp->targetGene)   
	aceOutf (ao, "\tGene") ;
      aceOutf (ao, "\tPosition\tStrand\tSNP") ;
      
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	{
	  aceOutf (ao, "\t%s Variant",  dictName (runDict, run)) ;
	  aceOutf (ao, "\t%s Coverage",  dictName (runDict, run)) ;
	  aceOutf (ao, "\t%s %%",  dictName (runDict, run)) ;
	}
      
      aceOutf (ao, "\tAbsent\t<=5%%\t5-90%%\t>=90%%") ;
      aceOutf (ao, "\n") ;
    }
  
  arraySort (snps, ssmDeepTableReportOrder) ;
  for (ii = nn = 0 ; ii < arrayMax (snps) ; ii++)
    {
      ssm = arrp (snps, ii, SSM) ;
      /* select the snp/positions to be exported */
      ok = FALSE ;
      variant = ssm->variant ; 
      for (jj = ii,  ssm2 = ssm ;
	   ! ok && jj < arrayMax (snps) && ssm2->variant == variant ;
	   ssm2++, jj++)
	{
	  if (strstr (dictName(eeDict, ssm->snp > 0 ? ssm->snp : - ssm->snp), ">n")) 
	    continue ;

	  if (snp->minCover + snp->minMutant + snp->minFrequency > 0)
	    {
	      if (snp->minCover > ssm2->coverS)
		continue ;
	      if (ssm2->coverS * snp->minFrequency > 100 * ssm2->count)
		continue ;
	      if (snp->minMutant > ssm2->count)
		continue ;
	      ok = TRUE ;
	    }
	  else
	    ok = TRUE ;
	}

           /* export the selected positions */
      if (ok)
	for (strand = 1 ; strand > -2 ; strand -= 2)
	  { /* export a report */
	    if (0) showSsm (snps) ;
	    int itg = 0 ;
	    int nAbsent, nAbove, nMiddle, nBelow ;

	    if (allStrands == 0x1 && strand != 1)
	      continue ;
	    if (allStrands == 0x2 && strand != -1)
	      continue ;
	      
	    
	    nAbsent = nAbove = nMiddle = nBelow = 0 ;
	    aceOutf (ao, "%s"
		     , ssm->target ? dictName (snp->targetDict, ssm->target) : "NULL"
		     ) ;
	    if (snp->target2dna)
	      {
		if (snp->targetDna && ssm->target < arrayMax(snp->target2dna))
		  {
		    int w = 8 ;
		    int x = arr (snp->target2dna, ssm->target, int) ;
		    int r1 = x + ssm->pos - 1 - w ;
		    int r2 = x + ssm->pos - 1 ;
		    int r3 = x + ssm->pos - 1 + w + 1 ;
		    char *cp, *cq ;
		    char *cp1 = stackText (snp->targetDna, r1) ; 
		    char *cp2 = stackText (snp->targetDna, r2) ; 
		    char *cp3 = stackText (snp->targetDna, r3) ; 
		    char buf[2*w + 2] ;
		    
		    buf[2*w + 1] = 0 ;
		    for (cp = cp1, cq = (strand == 1 ? buf : buf + 2*w)  ; cp < cp2 ; cp++, cq += strand)
		      {
			if (strand == -1)
			  *cq = ace_upper (dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int) *cp])]) ;
			else
			  *cq = ace_upper (*cp) ;
		      }
		    for ( ; cp <= cp2 ; cp++, cq += strand)
		      {
			if (strand == -1)
			  *cq =  ace_lower (dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int) *cp])]) ;
			else
			  *cq = ace_lower (*cp) ;
		      } 
		    for ( ; cp < cp3 ; cp++, cq += strand)
		      {
			if (strand == -1)
			  *cq = ace_upper (dnaDecodeChar[(int)complementBase(dnaEncodeChar[(int) *cp])]) ;
			else
			  *cq = ace_upper (*cp) ;
		      }
		    
		    aceOutf (ao, "\t%s", buf) ;
		  }
		else
		  aceOut (ao, "\tnnnnnnnnNnnnnnnnn") ;
	      }
	    if (snp->targetGene)
	      {
		/* itg is initialised as 0 outside of the ii loop */
		GT *gt ;
		int pos = ssm->pos ;
		BOOL nonCoding = FALSE ;

		for (gt = arrp (snp->targetGene, itg, GT) ; itg > 0 && gt->a1 > pos ; itg--, gt--) ;
		for ( ; itg < arrayMax (snp->targetGene) &&  gt->a2 < pos ; itg++, gt++) ;
		if ( itg < arrayMax (snp->targetGene) && gt->a1 <= pos && gt->a2 >= pos)
		  {
		    nonCoding = FALSE ;
		    aceOutf (ao, "\t%c%s:%d"
			     , gt->isUp ? '-' : '+'
			     , dictName (snp->geneDict, gt->gene)
			     , gt->isUp ? gt->a2 - pos + 1 : pos - gt->a1+1
			     ) ;
		    if (snp->targetDna && ssm->target < arrayMax(snp->target2dna))
		      {
			extern char codonMito (const char*) ;
			int x = arr (snp->target2dna, ssm->target, int) ;
			int r1 = x + pos - 1 ; /* base b subject to snp */
			char *cp1 = stackText (snp->targetDna, r1) ; 
			int i, dx ;
			char buf[4] ;
			char AA1, AA2 = 'X' ;
			int sn = ssm->snp > 0 ? ssm->snp : - ssm->snp ;
			const char *mysnp = dictName (eeDict, sn) ;
			
			if (! strncasecmp (dictName (snp->geneDict, gt->gene), "trna", 4) ||
			    ! strncasecmp (dictName (snp->geneDict, gt->gene), "rrn", 3)
			    )
			  nonCoding = TRUE ;
			    
			if (gt->isUp)
			  {
			    dx = (gt->a2 - pos) % 3 ; /* zero if b is the first base of the codon */
			    for (cp1 += dx, i = 0 ; i<3 ; cp1--, i++)
			      buf[i] = complementBase((int)dnaEncodeChar[(int) *cp1]) ;
			  }
			else
			  {
			    dx = (pos - gt->a1) % 3 ; /* zero if b is the first base of the codon */
			    for (cp1 -= dx, i = 0 ; i<3 ; cp1++, i++)
			      buf[i] = dnaEncodeChar[(int) *cp1] ;
			  }
			buf[3] = 0 ; /* block the codon */
			AA1 = snp->mito ? codonMito (buf) : codon (buf) ;
			if (mysnp[1] == '>')
			  {
			    if (gt->isUp)
			      {
				buf[dx] = complementBase(dnaEncodeChar[(int)mysnp[2]]) ;
			      }
			    else
			      {
				buf[dx] = dnaEncodeChar[(int)mysnp[2]] ;
			      }
			    AA2 = snp->mito ? codonMito (buf) : codon (buf) ;
			  }
			if (! nonCoding)
			  {
			    if (AA1 == AA2)
			      aceOutf (ao, "(==%c%d)", AA1
				       , gt->isUp ? 1 + (gt->a2 - pos)/3 : 1 + (pos - gt->a1)/3
				       ) ;
			    else
			      aceOutf (ao, "(%c%d%c)", AA1, 1 + (pos - gt->a1)/3, AA2) ;
			  }
		      }
		  }
		else
		  aceOut (ao, "\txxx") ;
	      }
	    
	    aceOutf (ao, "\t%s\t%d\t%s"
		     , strand > 0 ? "+" : "-"
		     , ssm->pos
		     , dictName (eeDict,  ssm->snp > 0 ? ssm->snp : - ssm->snp)
		     ) ;
	    
	    
	    for (run = 1 ; run <= dictMax (runDict) ; run++)
	      {			
		ok = 0 ;
		for (jj = ii ,  ssm2 = ssm ;
		     !ok && jj < arrayMax (snps) && ssm2->variant == variant ;
		     ssm2++, jj++)
		  {
		    /* export the relevant strand */
		    if (strand * ssm2->snp < 0  || run != ssm2->run)
		      continue ;

		    aceOutf (ao, "\t%d\t%d",  ssm2->count, ssm2->coverS) ;
		    z = ssm2->coverS ;
		    if (z == 0) z = 1 ;
		    z = 100 * ssm2->count/z ;
		    
		    if (snp->minCover + snp->minFrequency > 0)
		      {
			if (ssm2->coverS < snp->minCover)
			  z = -10 ;
		      }

		    aceOutf (ao, "\t%.1f",  z) ;
		    ok = 1 ; 
		    if (z >= 90) { nAbove++ ; array (cumulOver, run, int)++ ; }
		    else if (z <= -9) { nAbsent++ ;array (cumulAbsent, run, int)++ ; }
		    else if (z <= 5) { nBelow++ ;array (cumulBelow, run, int)++ ; }
		    else {  nMiddle++ ; array (cumulMiddle, run, int)++ ; }
		  }
		if (! ok) /* we found no support, report the observed cases */
		  {
		    aceOutf (ao, "\t0\t0\t-10") ;
		    nAbsent++ ; array (cumulAbsent, run, int)++ ;
		  }
	      }
	    aceOutf (ao, "\t%d\t%d\t%d\t%d", nAbsent, nBelow, nMiddle, nAbove) ;
	    aceOut (ao, "\n") ;
	  }
      for (jj = ii,  ssm2 = ssm ;
	   jj < arrayMax (snps) &&  ssm2->variant == variant ;
	   ssm2++, jj++) ;
      ii = jj - 1 ;
    }

  if (1)
    {
      aceOutf (ao, "p>=90\t\t\t\t\t") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	{
	  aceOutf (ao, "\t0\t0\t%d", array (cumulOver, run, int)) ;
	  array (cumulTotal, run, int) += array (cumulOver, run, int) ;
	}
      
      aceOutf (ao, "\np<=5\t\t\t\t\t") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	{
	  aceOutf (ao, "\t0\t0\t%d", array (cumulBelow, run, int)) ;
	  array (cumulTotal, run, int) += array (cumulBelow, run, int) ;
	}
      aceOutf (ao, "\nIntermediate\t\t\t\t\t") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	{
	  aceOutf (ao, "\t0\t0\t%d", array (cumulMiddle, run, int)) ;
	  array (cumulTotal, run, int) += array (cumulMiddle, run, int) ;
	}
      aceOutf (ao, "\nUndecidable\t\t\t\t\t") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	{
	  aceOutf (ao, "\t0\t0\t%d", array (cumulAbsent, run, int)) ;
	  array (cumulTotal, run, int) += array (cumulAbsent, run, int) ;
	}

      aceOutf (ao, "\nTotal\t\t\t\t\t") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	{
	  aceOutf (ao, "\t0\t0\t%d", array (cumulTotal, run, int)) ;
	}
      aceOutf (ao, "\n") ;
    }
  

  ac_free (h) ;
  return nn ;
} /* snpDeepTableReport */

/*************************************************************************************/

static int selectZoneOrder (const void *va, const void *vb)
{
  const ZONE *a = (const ZONE *)va, *b = (const ZONE *)vb ;
  int n ;

  n = a->zone - b->zone ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  n = a->vent - b->vent ; if (n) return n ;

  return 0 ;
} /* selectZoneOrder */

/*************************************************************************************/

static void snpPrepareZoneSelection (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *ccp ; int n = 0, zone ;
  Array aa ;
  ACEIN ai = aceInCreateFromFile (snp->selectFileName, "r", 0, h) ; 
  DICT *selectDict, *ventDict ;
  ZONE *zp ;
  KEYSET ks ;
  int dummy ;

  snp->selectDict = selectDict = dictHandleCreate (1000, snp->h) ;
  snp->ventilationDict = ventDict = dictHandleCreate (1000, snp->h) ;
  dictAdd (ventDict, "1", &dummy) ;
  snp->selectZone = aa = arrayHandleCreate (1000, ZONE, snp->h) ;
  snp->selectZoneIndex = ks = keySetHandleCreate (snp->h) ;
  if (ai)
    {
      while (aceInCard (ai))
	{
	  if ((ccp = aceInWord (ai)))
	    {
	      dictAdd (selectDict, ccp, &zone) ;
	      zp = arrayp (aa, n++, ZONE) ;
	      zp->zone = zone ; 
	      zp->a1 = 1 ; zp->a2 = 1<<30 ; zp->vent = dummy ;
	      if (aceInInt (ai, &zp->a1) &&  aceInInt (ai, &zp->a2) && (ccp = aceInWord (ai)))
		dictAdd (ventDict, ccp, &zp->vent) ;
	    }
	}
    }
  else
    usage (messprintf ("Cannot open the -select file : %s", snp->selectFileName)) ;

  arraySort (aa, selectZoneOrder) ;
  for (zp = arrp (aa, 0, ZONE), n = 0 ; n < arrayMax (aa) ; n++, zp++)
    {
      if (! keySet (ks, zp->zone)) 
	keySet (ks, zone) = n ;
    }

  if (! arrayMax (aa))
    usage (messprintf("The -select file is empty, I quit : %s", snp->selectFileName)) ;
  fprintf (stderr, "%s Found %d target, %d zones in the -select file: %s\n", timeShowNow(), dictMax(selectDict), arrayMax (aa), snp->selectFileName) ;
  ac_free (h) ;
  return ;
} /* snpPrepareZoneSelection */

/*************************************************************************************/
/* given a collection of snp txt file, find the intersect in all samples
 */
typedef struct wiStruct { float w ; int i ; } WI ;

static int wiOrder (const void *a, const void *b) 
{
  const WI *up = (const WI *)a, *vp = (const WI *)b ;
  
  if (up->w < vp->w) return -1 ;
  if (up->w > vp->w) return 1 ;
  return up->i - vp->i ;
} /* wiOrder */

/***************/

static void snpChronoOrder (Array aa, KEYSET aa1) 
{
  AC_HANDLE h = ac_new_handle () ;
  int iMax = keySetMax (aa1) ;
  int i, j, i2, j2, pass ;
  float w ;
  KEYSET col = keySetHandleCreate (h) ;
  KEYSET lin = keySetHandleCreate (h) ;
  Array ww = arrayHandleCreate (iMax, WI, h) ;

  for (j = 0 ; j < iMax ; j++)
    {
      keySet (col, j) = keySet (aa1, j) ;
      keySet (lin, j) = keySet (aa1, j) ;
    }

  for (pass = 0 ; pass < 10 ; pass++)
    {
      /* sort by weight of the lines */
      for (i = 0 ; i < iMax ; i++)
	{
	  i2 = keySet (lin, i) ;
	  for (w = 0, j = 0 ; j < iMax ; j++)
	    w += j * arr (aa, i2 * iMax + keySet (col,j), float) ;
	  array (ww, i, WI).w = w ;
	  array (ww, i, WI).i = i2 ;
	}
      arraySort (ww, wiOrder) ;
      for (i = 0 ; i < iMax ; i++)
	keySet (lin, i) = array (ww, i, WI).i ;

      /* sort by weight of the columns */
      for (j = 0 ; j < iMax ; j++)
	{
	  j2 = keySet (col, j) ;
	  for (w = 0, i = 0 ; i < iMax ; i++)
	    w += i * arr (aa, j2 + iMax * keySet (lin, i), float) ;
	  array (ww, j, WI).w = w ;
	  array (ww, j, WI).i = j2 ;
	}
      arraySort (ww, wiOrder) ;
      for (j = 0 ; j < iMax ; j++)
	keySet (col, j) = array (ww, j, WI).i ;
    }
  for (j = 0 ; j < iMax ; j++)
     keySet (aa1, j) = keySet (col, j) ;
  ac_free (h) ;
}

/***********/

static void snpIntersect (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = snp->ai ;
  ACEOUT ao = snp->ao ;
  DICT *runDict = dictHandleCreate (300, h) ;
  DICT *snpDict = dictHandleCreate (1000000, h) ;
  DICT *sampleDict = dictHandleCreate (300, h) ;
  const char *ccp ;
  const char *errors = 0 ;
  int pos, run, sn, ii, ir, jr, irMax, n ;
  vTXT txt = vtxtHandleCreate (h) ;
  BitSet bb, bbMM, bb1, bb2 ;
  Array bbMMs = arrayHandleCreate (200, BitSet, h) ;
  Array bbs = arrayHandleCreate (200, BitSet, h) ;
  Array aa = arrayHandleCreate (200, float, h) ;
  KEYSET histo = keySetHandleCreate (h) ;
  KEYSET histoMM = keySetHandleCreate (h) ;
  KEYSET histoCover = keySetHandleCreate (h) ;
  KEYSET intersect = keySetHandleCreate (h) ;
  KEYSET uni = keySetHandleCreate (h) ;
  KEYSET aa1 = keySetHandleCreate (h) ; 
  KEYSET covers = keySetHandleCreate (h) ;
  KEYSET run2sample = keySetHandleCreate (h) ;
  double u ;
  int cover ;
  float dosage, fq ;
  
  aceOutDate (ao, "###", snp->title) ;

  if (snp->db && snp->project)
    {
      char *qq = messprintf ("find run ; project == \"%s\"", snp->project) ;
      AC_KEYSET ks = ac_dbquery_keyset (snp->db, qq, h) ;
      AC_TABLE tbl = ac_keyset_table (ks, 0, 0, TRUE, h) ;
      const char *ccp ;
      int ir ;

      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  ccp = ac_table_printable (tbl, ir, 0, 0) ;
	  dictAdd (runDict, ccp, 0) ;
	}
    }

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#') 
	continue ;
      vtxtClear (txt) ;
      vtxtPrintf (txt, "%s:",ccp) ;

      if (!aceInInt (ai, &pos))
	continue ;

      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#') 
	continue ;
      vtxtPrintf (txt, "%d%s", pos, ccp) ;

      ccp = aceInWord (ai) ;  /* snipnet 1 */
      ccp = aceInWord (ai) ;  /* snipnet 2 */
      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#') 
	continue ;
      if (snp->db && snp->project && ! dictFind (runDict, ccp, &run))
	continue ;
      dictAdd (runDict, ccp, &run) ;

      dictAdd (snpDict, vtxtPtr (txt), &sn) ;

      dosage = 0 ;
      aceInFloat (ai, &dosage) ;   /* genotype */

      if (! aceInFloat (ai, &fq) || fq <  snp->minFrequency)
	continue ;
      if (!aceInInt (ai, &cover) || cover < 4 || cover < snp->minCover)
	continue ;
      keySet (covers, sn)++ ;
      if (dosage < .1)
	continue ;
      bb = array (bbs, run, BitSet) ;
      if (! bb)
	bb = array (bbs, run, BitSet) = bitSetCreate (1000000, h) ;
      bitSet (bb, sn) ;
      bbMM = array (bbMMs, run, BitSet) ;
      if (! bbMM)
	bbMM = array (bbMMs, run, BitSet) = bitSetCreate (1000000, h) ;
      if (dosage > 1.9)
	bitSet (bbMM, sn) ;
    }

  /* make sure all BitSets are large enough so we use bit in the rest of the code */
  sn =  dictMax (snpDict) + 1 ;
  for (ii = 1 ; ii < arrayMax (bbs) ; ii++)
    {
      bb = arr (bbs, ii, BitSet) ;
      if (! bb)
	bb = array (bbs, run, BitSet) = bitSetCreate (dictMax (snpDict), h) ;
      bitUnSet (bb, sn) ;
    }
  for (ii = 1 ; ii < arrayMax (bbMMs) ; ii++)
    {
      bbMM = arr (bbMMs, ii, BitSet) ;
      if (! bbMM)
	bbMM = array (bbMMs, run, BitSet) = bitSetCreate (dictMax (snpDict), h) ;
      bitUnSet (bbMM, sn) ;
    }
  
  /* count of how many sample see a given snp */
  for (sn = 1 ; sn <= dictMax (snpDict) ; sn++)
    {
      keySet (histoCover, keySet (covers, sn))++ ;
      for (n = 0, ii = 1 ; ii < arrayMax (bbs) ; ii++)
	{
	  bb = arr (bbs, ii, BitSet) ;
	  if (bb && bit (bb, sn)) n++ ;
	}
      keySet (histo, n) += 1 ;
      for (n = 0, ii = 1 ; ii < arrayMax (bbMMs) ; ii++)
	{
	  bbMM = arr (bbMMs, ii, BitSet) ;
	  if (bbMM && bit (bbMM, sn)) n++ ;
	}
      keySet (histoMM, n) += 1 ;
    }

  /* a few variants (the breakpoints ?) seem to be present twice in the files, we fix that */
  for (n = ii = 0 ; n > dictMax (runDict) && n < keySetMax (histoCover) ; n++)
    ii +=  keySet (histoCover, n) ;
  keySet (histoCover, dictMax (runDict)) += ii ;


  if (snp->db)
    {
      AC_TABLE tbl = 0 ;
      const char *ccp ;
      vTXT txt = vtxtHandleCreate (h) ;
      int ir, ii, n ;
      
      vtxtPrintf (txt, "table = ") ;
      vtxtPrintf (txt, "colonne 1 \n class Run \n\n") ;
      vtxtPrintf (txt, "colonne 2 \nFrom 1\nClass Sample\nTag Sample\n\n") ;
      tbl = ac_tablemaker_table (snp->db, vtxtPtr (txt), 0, ac_tablemaker_text , 0 , 0, &errors, h) ;

      for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
	{
	  ccp = ac_table_printable (tbl, ir, 0, 0) ;
	  if (ccp && dictFind (runDict, ccp, &ii))
	    {
	      ccp = ac_table_printable (tbl, ir, 1, 0) ;
	      if (ccp)
		{
		  dictAdd (sampleDict, ccp, &n) ;
		  keySet (run2sample, ii) = n ;
		}
	    }
	}
    }

  aceOutf (ao, "Number of samples") ;
  for (n = 0 ; n <= dictMax (runDict) ; n++) 
    aceOutf (ao, "\t%d", n) ;
  aceOutf (ao, "\tTotal\nNumber of variants measured in n samples") ;
  for (n = ii = 0 ; n <= dictMax (runDict) && n < keySetMax (histoCover) ; n++) 
    {
      ii += keySet (histoCover, n) ;
      aceOutf (ao, "\t%d", keySet (histoCover, n)) ;
    }
  aceOutf (ao, "\nNumber of variants seen as m/m or m/+ in n samples") ;
  for (n = ii = 0 ; n < keySetMax (histo) ; n++) 
    {
      ii += keySet (histo, n) ;
      aceOutf (ao, "\t%d", keySet (histo, n)) ;
    }
  aceOutf (ao, "\nNumber of variants seen as m/m in n samples") ;
  for (n = ii = 0 ; n < keySetMax (histoMM) ; n++) 
    {
      ii += keySet (histoMM, n) ;
      aceOutf (ao, "\t%d", keySet (histoMM, n)) ;
    }
  aceOutf (ao, "\n") ;


  irMax = arrayMax (bbs) ;
  for (ir = 1 ; ir < irMax ; ir++)
    {
      bb1 = arr (bbs, ir, BitSet) ;
      if (! bb1)
	bb1 = array (bbs, run, BitSet) = bitSetCreate (dictMax (snpDict), h) ;
      for (jr = 1 ; jr < irMax ; jr++)
	{
	  bb2 = arr (bbs, jr, BitSet) ;
	  if (! bb2)
	    bb2 = array (bbs, run, BitSet) = bitSetCreate (dictMax (snpDict), h) ;
	   for (sn = 1 ; sn <= dictMax (snpDict) ; sn++)
	    {
	      if (1)
		{
		  if (bit (bb1, sn) && bit (bb2, sn)) 
		    keySet (intersect, ir * irMax + jr)++ ;
		  if (bit (bb1, sn) || bit (bb2, sn)) 
		    keySet (uni, ir * irMax + jr)++ ;
		}
	    }
	}
    }

  aceOutf (ao, "\n\nNumber of variants seen as mutant in either run") ;
  for (jr = 1 ; jr <= irMax ; jr++) aceOutf (ao, "\t") ;
  aceOutf (ao, "\t\tNumber of variants seen in either run") ;

  aceOut (ao, "\nRun") ;  if (snp->db) aceOut (ao, "\t") ;
  for (jr = 1 ; jr < irMax ; jr++) aceOutf (ao, "\t%s", dictName (runDict, jr)) ;
  aceOut (ao, "\t\tRun") ;  if (snp->db) aceOut (ao, "\t") ;
  for (jr = 1 ; jr < irMax ; jr++) aceOutf (ao, "\t%s", dictName (runDict, jr)) ;

  if (snp->db)
    {
      aceOutf (ao, "\n\tSample") ;
      for (jr = 1 ; jr < irMax ; jr++)
	{
	  n =  keySet (run2sample,jr) ; 
	  aceOutf (ao, "\t%s", n ? dictName (sampleDict, n) : "") ;
	}
      aceOutf (ao, "\t\t\tSample") ;
      for (jr = 1 ; jr < irMax ; jr++) 
	{
	  n =  keySet (run2sample,jr) ; 
	  aceOutf (ao, "\t%s", n ? dictName (sampleDict, n) : "") ;
	}
    }

  for (ir = 1 ; ir < irMax ; ir++)
    {
      aceOutf (ao, "\n%s", dictName (runDict, ir)) ;
      if (snp->db)
	{
	  n =  keySet (run2sample,keySet (aa1,ir)) ; 
	  aceOutf (ao, "\t%s", n ? dictName (sampleDict, n) : "") ;
	}
      for (jr = 1 ; jr < irMax ; jr++)
	aceOutf (ao, "\t%d", keySet (uni, ir * irMax + jr)) ;
      aceOutf (ao, "\t\t%s", dictName (runDict, ir)) ;
      if (snp->db)
	{
	  n =  keySet (run2sample,keySet (aa1,ir)) ; 
	  aceOutf (ao, "\t%s", n ? dictName (sampleDict, n) : "") ;
	}
      for (jr = 1 ; jr < irMax ; jr++)
	{
	  u = keySet (intersect, ir * irMax + ir) ;
	  if (u == 0) u =  1 ;
	  aceOutf (ao, "\t%.1f", 100.0 * keySet (uni, ir * irMax + jr)/u) ;
	}
    }


  aceOutf (ao, "\n\nNumber of variants seen in both runs\t") ;
  for (jr = 1 ; jr < irMax ; jr++) aceOutf (ao, "\t") ;
  aceOutf (ao, "\t\tNumber of variants seen in both runs\t") ;

  aceOut (ao, "\nRun") ;  if (snp->db) aceOut (ao, "\t") ;
  for (jr = 1 ; jr < irMax ; jr++) aceOutf (ao, "\t%s", dictName (runDict, jr)) ;
  aceOut (ao, "\t\tRun") ;  if (snp->db) aceOut (ao, "\t") ;
  for (jr = 1 ; jr < irMax ; jr++) aceOutf (ao, "\t%s", dictName (runDict, jr)) ;

  if (snp->db)
    {
      aceOutf (ao, "\n\tSample") ;
      for (jr = 1 ; jr < irMax ; jr++)
	{
	  n =  keySet (run2sample,jr) ; 
	  aceOutf (ao, "\t%s", n ? dictName (sampleDict, n) : "") ;
	}
      aceOutf (ao, "\t\t\tSample") ;
      for (jr = 1 ; jr < irMax ; jr++) 
	{
	  n =  keySet (run2sample,jr) ; 
	  aceOutf (ao, "\t%s", n ? dictName (sampleDict, n) : "") ;
	}
    }

  for (ir = 1 ; ir < irMax ; ir++)
    {
      aceOutf (ao, "\n%s", dictName (runDict, ir)) ;  
      if (snp->db)
	{
	  n =  keySet (run2sample,keySet (aa1,ir)) ; 
	  aceOutf (ao, "\t%s", n ? dictName (sampleDict, n) : "") ;
	}
      for (jr = 1 ; jr < irMax ; jr++)
	aceOutf (ao, "\t%d", keySet (intersect, ir * irMax + jr)) ;
      aceOutf (ao, "\t\t%s", dictName (runDict, ir)) ;
      if (snp->db)
	{
	  n =  keySet (run2sample,keySet (aa1,ir)) ; 
	  aceOutf (ao, "\t%s", n ? dictName (sampleDict, n) : "") ;
	}
      for (jr = 1 ; jr < irMax ; jr++)
	{
	  int ni, nu ;

	  u = keySet (intersect, ir * irMax + ir) ;
	  ni = keySet (intersect, ir * irMax + jr) ;
	  if (u == 0) u =  1 ;
	  aceOutf (ao, "\t%.1f", 100.0 * ni/u) ;

	  u = nu = keySet (uni, ir * irMax + jr) ; 
	  /* complement of set theoretic distance */
	  array (aa, ir * irMax + jr, float) = nu > 0 ? 100.0 * (ni/u) : 0 ;

	  /* new idea: intersect divided by smaller */
	  u = keySet (intersect, ir * irMax + ir) ; 
	  nu = keySet (intersect, jr * irMax + jr) ; 
	  if (u > nu) u = nu ; 
	  if (u == 0) u =  1 ;
	  array (aa, ir * irMax + jr, float) = 100.0 * (ni/u) ;
	}
    }

  aceOutf (ao, "\n\n") ;

  /* export the best pairs */ 
  for (ir = 0 ; ir < irMax ; ir++) 
    {
      for (jr = 1 ; jr < irMax ; jr++)
	{
	}
    }
  /* reorder the aa matrix of relative intersects and reexport it */
  for (ir = 0 ; ir < irMax ; ir++) 
    keySet (aa1, ir) = ir ;
  snpChronoOrder (aa, aa1) ;
  if (0)
    aceOutf (ao, "\n\nSet theoretic distance between reordered samples") ;
  if (1)
    aceOutf (ao, "\n\nIntersect/smaller category") ;
  aceOutf (ao, "\nRun") ;
  if (snp->db)
    aceOutf (ao, "\tSample") ;
  for (jr = 1 ; jr < irMax ; jr++) aceOutf (ao, "\t%s", dictName (runDict, keySet (aa1,jr))) ;
  if (snp->db)
    {
      aceOutf (ao, "\nRun\tSample") ;
      for (jr = 1 ; jr < irMax ; jr++)
	{
	  n =  keySet (run2sample,keySet (aa1,jr)) ; 
	  aceOutf (ao, "\t%s", n ? dictName (sampleDict, n) : "") ;
	}
    }

  for (ir = 1 ; ir < irMax ; ir++)
    {
      aceOutf (ao, "\n%s", dictName (runDict, keySet (aa1,ir))) ; 
      if (snp->db)
	{
	  n =  keySet (run2sample,keySet (aa1,ir)) ; 
	  aceOutf (ao, "\t%s", n ? dictName (sampleDict, n) : "") ;
	}
      for (jr = 1 ; jr < irMax ; jr++)
	{
	  aceOutf (ao, "\t%.1f", array (aa, keySet (aa1,ir) * irMax + keySet (aa1,jr), float)) ;
	}
    }

  aceOutf (ao, "\n\n") ;


  ac_free (h) ;
} /* snpIntersect */

/*************************************************************************************/
/*************************************************************************************/
/* merge files exported by the -count2genome option */
typedef struct snpReportStruct {int sn, chrom, pos, type, sq1, sq2, run, cover, m, w, mp, wp, mm, wm, bad, multip, multim, oo1p, oo2p, oo1m, oo2m ; float frequency, dose ; } SNR ;

static int snrOrder (const void *a, const void *b)
{
  const SNR *up = (const SNR *)a, *vp = (const SNR *)b ;
  int n ;

  n = up->sn - vp->sn ; if (n) return n ;
  n = up->run - vp->run ; if (n) return n ;

  return 0 ;
} /* snrOrder */

/*************************************************************************************/

static BOOL snpParseOneRecord (ACEIN ai, SNR *sp, vTXT txt, DICT *runDict, DICT *snpDict, Stack snipnet, BOOL *solidp)
{
  const char *ccp ;

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#') 
	continue ;
      vtxtClear (txt) ;
      memset (sp, 0, sizeof(SNR)) ;
      dictAdd (runDict, ccp, &(sp->chrom)) ;
      vtxtPrintf (txt, "%s:",ccp) ;
      
      if (!aceInInt (ai, &(sp->pos)) || sp->pos <= 0)
	continue ;

      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#') 
	continue ;
      dictAdd (runDict, ccp, &(sp->type)) ;
      vtxtPrintf (txt, "%d%s", sp->pos, ccp) ;

      ccp = aceInWord (ai) ;
      if (ccp) { sp->sq1 = stackMark (snipnet) ; pushText (snipnet, ccp) ; vtxtPrintf (txt, ".%s", ccp) ;} 
      ccp = aceInWord (ai) ;
      if (ccp) { sp->sq2 = stackMark (snipnet) ; pushText (snipnet, ccp) ; vtxtPrintf (txt, ".%s", ccp) ;} 

      /* in case of BRK points the snipnet contains the coords of the second chromosome */
      dictAdd (snpDict, vtxtPtr (txt), &(sp->sn)) ;

      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#') 
	continue ;
      dictAdd (runDict, ccp, &(sp->run)) ;

      ccp = aceInWord (ai) ; /* dosage drop it */
      sp->dose = (ccp && *ccp != '-' ? atof (ccp) : -1) ;
      aceInFloat (ai, &(sp->frequency)) ;
      aceInInt (ai, &(sp->cover)) ;
      aceInInt (ai, &(sp->m)) ;
      aceInInt (ai, &(sp->w)) ;
      aceInInt (ai, &(sp->mp)) ;
      aceInInt (ai, &(sp->wp)) ;
      aceInInt (ai, &(sp->mm)) ;
      aceInInt (ai, &(sp->wm)) ;
      aceInInt (ai, &(sp->bad)) ;
      ccp = aceInWord (ai) ; 
      if (ccp)
	{
	  char *cq = strstr (ccp, ":") ;
	  if (cq) 
	    { *solidp = TRUE ; sscanf (ccp, "%d:%d", &(sp->oo1p), &(sp->oo2p)) ; }
	  else
	    sscanf (ccp, "%d", &(sp->oo2p)) ;
	}
      ccp = aceInWord (ai) ; 
      if (ccp)
	{
	  char *cq = strstr (ccp, ":") ;
	  if (cq) 
	    { *solidp = TRUE ; sscanf (ccp, "%d:%d", &(sp->oo1m), &(sp->oo2m)) ; }
	  else
	    sscanf (ccp, "%d", &(sp->oo2m)) ;
	}
      ccp = aceInWord (ai) ;  /* chi2 strand flag, drop it */
      ccp = aceInWord (ai) ;  /* oo flag, drop it */

      aceInInt (ai, &(sp->multip)) ;
      aceInInt (ai, &(sp->multim)) ;
      return TRUE ;
    }
  return FALSE ;
} /* snpParseOneRecord */

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

static DICT *myTargetDict ;

static int targetOrder (const void *a, const void *b)
{
  const int *up = (const int *)a, *vp = (const int *)b ;
  const char *tu = dictName (myTargetDict, *up) ;
  const char *tv = dictName (myTargetDict, *vp) ;

  return lexstrcmp (tu, tv) ;
  return 0 ;
} /* snrOrder */

/*************************************************************************************/

static KEYSET snpBRS_sort_target (SNP *snp, AC_HANDLE h)
{
  Array ttt = arrayHandleCreate (dictMax(snp->targetDict)+1, int, h) ;
  int ii = dictMax (snp->targetDict) ; 
  
  while (ii--)
    array (ttt, ii, int) = ii + 1 ;
  
  myTargetDict = snp->targetDict ;
  arraySort (ttt, targetOrder) ;
  
  return ttt ;
} /* snpBRS_sort_target */

/*************************************************************************************/
/*************************************************************************************/
/* expect as input a concatenation of BRS files
 * merge them into a single file
 *
 * we parse the 7 columns into a data structure, 
 * merge on the fly
 * reexport at the end in the same format
 */
static int snpBRSMergeOrder (const void *va, const void *vb)
{
  const BRS *a = (const BRS *)va, *b = (const BRS *)vb ;
  int n ;

  n = a->pos - b->pos ; if (n) return n ;
  n = a->ee - b->ee ; if (n) return n ;

  return 0 ;
} /*  snpBRSMergeOrder */

/*************************************************************************************/

static int snpBRSMergeExport (SNP *snp, ACEOUT ao, BigArray aa, int target)
{
  BRS *brs ;
  int nn = 0 ;
  long int i, iMax = bigArrayMax (aa) ;
  const char *snpRun = snp->run ;
  const char *targetName = dictName (snp->targetDict, target) ;
  DICT *eeDict = snp->eeDict ;
  char *strands = "+-" ;

  for (int pass = 0 ; pass < 2 ; pass++)
    for (i = 0, brs = bigArrp (aa, i, BRS) ; i < iMax ; i++, brs++)
      {
	if (brs->bp == -1)
	  continue ;
	if (brs->ee == -1) 
	  continue ;
	if (pass == 1 && ! brs->bm)
	  continue ;
	if (pass == 0 && ! brs->bp)
	  continue ;
	aceOutf (ao, "%s\t%d\t%s\t%c\t%s"
		 , nn ? "~" : targetName
		 , brs->pos
		 , dictName (eeDict, brs->ee)
		 , strands[pass]
		 , nn ? "~" : snpRun 
	       ) ;
	nn++ ;
	if (pass == 1)
	  aceOutf (ao, "\t%d\t%d\t%d\t%d\n"
		   , brs->bm, brs->rm, brs->sm, brs->tm
		   ) ;
	else
	  aceOutf (ao, "\t%d\t%d\t%d\t%d\n"
		   , brs->bp, brs->rp, brs->sp, brs->tp
		   ) ;
      }
  return nn ;
} /* snpBRSMergeExport */

/*************************************************************************************/

static void snpBRSParseCompress (BigArray aa)
{
  register BRS *up, *vp, *wp ;
  register long int ii, i, j, iMax ;

  bigArraySort (aa, snpBRSMergeOrder) ;
  iMax = bigArrayMax (aa) ;
  for (i = ii = 0, up = wp = bigArrp (aa, 0, BRS) ; i < iMax ; i++, up++)
    {
      if (up->bp == -1) continue ;
      int pos = up->pos ;
      short ee = up->ee ;

      for (j = i + 1, vp = up + 1 ;
	   j < iMax && vp->pos == pos && vp->ee == ee ;
	   j++, vp++)
	{
	  up->bp += vp->bp ;
	  up->bm += vp->bm ;
	  up->rp += vp->rp ;
	  up->rm += vp->rm ;
	  up->sp += vp->sp ;
	  up->sm += vp->sm ;
	  up->tp += vp->tp ;
	  up->tm += vp->tm ;
	  vp->bp = -1 ;
	}
      if (ii < i) { *wp = *up ; }
      ii++ ; wp++ ;
    }	
  bigArrayMax (aa) = ii ;
  return ;					\
}  /* snpBRSParseCompress */

/*************************************************************************************/

static int snpBRSparse (SNP *snp)
{
  BRS *brs ;
  ACEIN ai = snp->ai ;
  DICT *eeDict = snp->eeDict ;
  DICT *targetDict = snp->targetDict ;

  BigArray aa = 0 ;  /* ee B array of the current target */
  Array aaa = snp->aaa ;
  BigArray aaP = 0 ;  /* RS array of the current target */
  Array aaaP = 0 ;
  Array errCounts = snp->errCounts ;

  int ee, x, state = 0, oldTarget = 0, oldTargetDummy = 0 ;
  int nnFiles = 0 ; /* count the BRS files */
  int nnIn = 0 ;    /* global line counts */
  long int nAli = 0 ;  /* global read count */
  long int nBpAli = 0 ; /* global base count */

  const char *ccp ;
  
  if (snp->BRS_detect || snp->BRS_count)
    aaaP = snp->aaaP ;

  aceInSpecial (ai, "\"\n\t") ;
  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#')
	continue ;
      else if (! strcmp (ccp, "nAli"))  /* overall read count */
	{
	  nnFiles++ ;
	  state = 0 ;
	  aceInStep (ai,'\t') ; 
	  if (aceInInt (ai, &x))
	    nAli += x ;
	}
      else if (! strcmp (ccp, "nBpAli"))  /* overall base count */
	{
	  state = 0 ;
	  aceInStep (ai,'\t') ; 
	  if (aceInInt (ai, &x)) 
	    nBpAli += x ;
	}
      else if (! strcmp (ccp, "ERR"))     /* overall error count */
	{
	  state = 0 ;
	  aceInStep (ai,'\t') ;  ccp = aceInWord (ai) ;
	  dictAdd (eeDict, ccp, &ee) ;
	  
	  aceInStep (ai,'\t') ; aceInInt (ai, &x) ; /* dummy 100, drop it */ 
	  aceInStep (ai,'\t') ;
	  if (aceInInt (ai, &x)) 
	    array (errCounts, ee, long int) += x ;
	  continue ;
	}

      if (state == 0 && aa)
	{
	  snpBRSParseCompress (aa) ;
	  aa = aaP = 0 ; state = oldTarget = oldTargetDummy = 0 ;
	}

      if (state == 0 && *ccp == '~')  /* target not specified */
	continue ;

      if (*ccp != '~')
	{
	  int target ;

	  state = 1 ;
	  dictAdd (targetDict, ccp, &target) ;
	  if (target != oldTarget)
	    {
	      if (aa)
		{
		  snpBRSParseCompress (aa) ;
		}
	      aa = array (aaa, target, BigArray) ;
	      if (! aa)
		aa = array (aaa, target, BigArray) = bigArrayHandleCreate (10000, BRS, snp->h) ;
	      if (aaaP)
		{
		  aaP = array (aaaP, target, BigArray) ;
		  if (! aaP)
		    aaP = array (aaaP, target, BigArray) = bigArrayHandleCreate (10000, BRS2, snp->h) ;
		}
	      oldTarget = target ;
	    }
	}
      

      if (state == 1)
	{
	  BRS B ;
	  char base, strand ;

	  aceInStep (ai,'\t') ; if (! aceInInt (ai, &B.pos)) continue ;
	  aceInStep (ai,'\t') ; if (! (ccp = aceInWord (ai))) continue ;
	  B.type = 0 ;
	  if (ccp[1] == 0)
	    B.type = 0 ;
	  else if (ccp[1] == '>')
	    B.type = 1 ;
	  else if (ccp[0] == '+' || ccp[1] == '+') 
	    B.type = 2 ;
	  else if (ccp[0] == '-' || (ccp[0] == '*' && ccp[1] == '-')) 
	    B.type = 3 ;
	  else
	    messcrash ("Bad snp descriptor %s in column 3 line %d of BRS file %s\n"
		       , ccp
		       , aceInStreamLine (ai), aceInFileName (ai)
		       ) ;
	  dictAdd (eeDict, ccp, &ee) ; B.ee = (short) (ee & 0xffff) ; 
	  base = ccp[0] ;
	  aceInStep (ai,'\t') ; if (! (ccp = aceInWord (ai))) continue ;
	  strand = *ccp ;
	  B.bm = 0 ;
	  aceInStep (ai,'\t') ; if (! aceInWord (ai)) continue ;
	  aceInStep (ai,'\t') ; if (! aceInInt (ai, &B.bp)) continue ;
	  aceInStep (ai,'\t') ; if (! aceInInt (ai, &B.rp)) continue ;
	  aceInStep (ai,'\t') ; if (! aceInInt (ai, &B.sp)) continue ;
	  aceInStep (ai,'\t') ; if (! aceInInt (ai, &B.tp)) continue ;

	  brs = bigArrayp (aa, bigArrayMax(aa), BRS) ;
	  brs->pos = B.pos ;
	  brs->type = B.type ;
	  brs->ee = B.ee ;
	  if (strand == '+')
	    { 
	      brs->bp = B.bp ; 
	      brs->rp = B.rp ; 
	      brs->sp = B.sp ; 
	      brs->tp = B.tp ; 
	    }
	  else
	    { 
	      brs->bm = B.bp ; 
	      brs->rm = B.rp ; 
	      brs->sm = B.sp ; 
	      brs->tm = B.tp ; 
	    }

	  if (aaP && ! B.type)
	    {
	      BRS2 *brs2 = bigArrayp (aaP, B.pos, BRS2) ;
	      if (strand == '+')
		{
		  brs2->bp += B.bp ;
		  brs2->rp += B.rp ;
		  brs2->sp += B.sp ;
		  brs2->tp += B.tp ;
		}
	      else
		{
		  brs2->bm += B.bp ;
		  brs2->rm += B.rp ;
		  brs2->sm += B.sp ;
		  brs2->tm += B.tp ;
		}
	      brs2->base = base ;
	    }
	  nnIn++ ;
	}      
    }
  if (aa)
    snpBRSParseCompress (aa) ; 

  fprintf (stderr, "//  snpParse found %d files %d lines %ld ali %ld bases in file %s\n"
	   , nnFiles, nnIn, nAli, nBpAli
	   , aceInFileName (ai)
	   ) ;
  return nnIn ;
} /* snpBRSparse */

/*************************************************************************************/
/* insure that to each snp corresponds a wild type line */
static int snpBRScompleteOne (SNP *snp, BigArray aa, BigArray aaP)
{
  long ii, iiMax = bigArrayMax (aa) ;
  BRS *brs ;
  BRS2 *brs2 ;
  DICT *eeDict = snp->eeDict ;

  for (ii = 0, brs = bigArrp (aa, 0, BRS) ; ii < iiMax ; ii++, brs++)
    {
      if (brs->bp == -1) continue ;
      if (brs->ee < 1) continue ;
      brs2 = bigArrayp (aaP, brs->pos, BRS2) ;
      if (! brs2->base)
	{
	  const char *ee = dictName (eeDict, brs->ee) ;
	  while (*ee == '*') ee++ ;
	  while (*ee == '-') ee++ ;
	  brs2->base = ace_lower (ee[0]) ;
	}
      if (brs->sp && ! brs2->sp)
	{
	  brs2->rp = brs->rp ;
	  brs2->sp = brs->sp ;
	  brs2->tp = brs->tp ;
	}
      if (brs->sm && ! brs2->sm)
	{
	  brs2->rm = brs->rm ;
	  brs2->sm = brs->sm ;
	  brs2->tm = brs->tm ;
	}
    }
  return 0 ;
}

/*************************************************************************************/

static int snpBRScomplete (SNP *snp)
{
  int nInc = 0 ;
  BigArray aa = 0 ;  /* ee B array of the current target */
  BigArray aaP = 0 ;  /* RS array of the current target */
  Array aaa = snp->aaa ;  /* ee B array of the current target */
  Array aaaP = snp->aaaP ;  /* ee B array of the current target */
  int target, targetMax = arrayMax (aaa) ;
  
  /* check each target */
  for (target = 1 ; target < targetMax ; target++)
    {
      aa = array (aaa, target, BigArray) ;
      aaP = array (aaaP, target, BigArray) ;
      if (aa && bigArrayMax (aa))
	{
	  nInc += snpBRScompleteOne (snp, aa, aaP) ;
	}
    }

  return nInc ;
} /* snpBRScomplete */

/*************************************************************************************/
/*************************************************************************************/

static int snpExportBRS (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  DICT *eeDict = snp->eeDict ;
  int nOut = 0 ;
  Array errCounts = snp->errCounts ;
  BigArray aa = 0 ;  /* ee B array of the current target */
  Array aaa = snp->aaa ;  /* ee B array of the current target */
  int ee ;
  long int nAli = 0  ;  /* gloabal read count */
  long int nBpAli = 0 ; /* gloabal base count */
  ACEOUT ao = aceOutCreate (snp->outFileName, ".snp", snp->gzo, h) ;
  aceOutDate (ao, "##", snp->title ? snp->title : "BRS file ") ; 

  /* export start  */
  /* export global counts */
  aceOutf (ao, "nAli\t%ld\nnBpAli\t%ld\n", nAli, nBpAli) ;
  for (ee = 1 ; ee <= dictMax (eeDict) ; ee++)
    {
      long int z = array (errCounts, ee, long int) ;
      if (z > 0)
	aceOutf (ao, "ERR\t%s\t100\t%ld\n"
		 , dictName (eeDict, ee)
		 , z
		 ) ;
    }
  
  /* export the caption */
  aceOutf (ao, "## Attention: the BRST counts are multiplied by up to 100: the a posteriori quality\n") ;
  aceOutf (ao, "# Target\tCoordinate\tBase (W=wild_type)\tStrand\tRun\tB\tR\tS\tT\n") ;
  
  /* export each target */
  {
    Array ttt = snpBRS_sort_target (snp, h) ;

    for (int ii = 0 ; ii < arrayMax (ttt) ; ii++)
      {
	int target = array (ttt, ii, int) ;
	aa = array (aaa, target, BigArray) ;
	if (aa && bigArrayMax (aa))
	  {
	    nOut += snpBRSMergeExport (snp, ao, aa, target) ;
	  }
      }
  }
  
  aceOut (ao, "\n") ;
  return nOut ;
} /* snpExportBRS */

/*************************************************************************************/

static void snpBRSMerge (SNP *snp)
{
  int nIn = snpBRSparse (snp) ;
  int nOut = snpExportBRS (snp) ;
  
  int mx = 0 ;
  messAllocMaxStatus (&mx) ;   
  fprintf (stderr, "... snpBRSMerge       \t%s\tmax memory %d Mb imported %d  BRS lines exported %d\n", timeShowNow(), mx, nIn, nOut) ;

  return ;
} /* snpBRSMerge */

/*************************************************************************************/
/******************************* snpBRSMerge end *************************************/
/*************************************************************************************/

static void snpBRS2snpsExportHeader (SNP *snp, ACEOUT ao, int prefix)
{
  BOOL showTitle = TRUE ;

  if (showTitle)
    {
      aceOutDate (ao, "###", snp->project) ;
      aceOutf (ao, 
	       "##  1: SnpName :: target:Position:type:Reference:Variant\n"
	       "##      -we use VCF conventions\n"
	       "##      -the target is a chromosome or transcript,\n"
	       "##      -all our coordinates start at 1 (not zero),\n"
	       "##      -the type is either Sub, Ins or Del,\n"
	       "##      -for a substitution, x is the position of the modifed base,\n"
	       "##      -for an indel, x is the position of the conserved hook base left of the variation\n"
	       "##      -all sliding indels are pushed left in the orientation of the target,\n"
	       "##      -all bases correspond to the plus strand of the reference,\n"
	       "##      -example t:100:Sub:A:C, base 100 A of target t is modified as C,\n"
	       "##               t:100:Del:aG:a, base 101 G of target t is deleted,\n"
	       "##               t:100:Ins:a:aG, base G is inserted between A100 and G101.\n"
	       "##  2: Run name.\n"
	       "##  3: 6i, format descriptor announcing 6 additive integers.\n"
	       "##  4,5,6: Coverage, Mutant, Reference counts seen on plus strand of the reads.\n"
	       "##  7,8,9: Coverage, Mutant, Reference counts seen on minus strand of the reads.\n"
	       ) ;
      
      if (0 && snp->solid)
	aceOut (ao,
		"##      -For SOLiD sequencing, we give 2 numbers for each strand, \n"
		"##      for example gTc 2:3  1:0, to indicate on the plus strand\n"
		"##      two reads not matching the gT transition, one not matchingn Tc\n"
		"##      and on the minus strand, one read not matching gT\n"
		) ;
    }
  
  aceOut (ao, "#\n#SNP\tRun\tFormat\tCover +\tMutant +\tReference +\tCover -\tMutant -\tReference -\n") ;
  
  return ;
} /* snpBRS2snpsExportHeader */

/*************************************************************************************/

static char *snpBRSname (SNP *snp, int target, BRS *brs, BRS2 *brs2, char *buf)
{
  char *type = "XXX" ;
  int k, k2, kk, x, pos = brs->pos ;
  const char *ee = dictName (snp->eeDict, brs->ee) ;
  char Abuf[10], Bbuf[10] ;
  BRS2 *b ;

  x = pos ;
  switch (brs->type)
    {
    case 1: /* sub */
      type = "Sub" ;
      Abuf[0] = ace_upper(brs2->base) ;  Abuf[1] = 0 ;
      Bbuf[0] = ace_upper(ee[2]) ;  Bbuf[1] = 0 ;
      break ;
    case 2: /* ins */
      type = "Ins" ;
      if (*ee == '*') ee++ ;
      if (ee[1] == '+') ee++ ;
      while (*ee && *ee == '+') ee++ ;
      kk = 0 ; while (ee[kk]) kk++ ;
      for (x = pos - 1, b = brs2 - 1, k = kk - 1 ; x > 1 && b->base == ee[k] ; b--, x--, k =(kk + k - 1) % kk) 
	;
      Abuf[0] = ace_lower(b->base) ;  Abuf[1] = 0 ;
      Bbuf[0] = ace_lower(b->base) ; 
      for (k2 = 0 ; k2 < kk ; k2++) 
	Bbuf[k2+1] = ace_upper (ee[(k+kk+k2+1)%kk]) ; 
      Bbuf[k2+1] = 0 ;
      break ;
    case 3: /* del */
      type = "Del" ;
      if (*ee == '*') ee++ ;
      while (*ee && *ee == '-') ee++ ;
      kk = 0 ; while (ee[kk]) kk++ ;
      for (x = pos - 1, b = brs2 - 1, k = kk - 1 ; x > 1 && b->base == ee[k] ; b--, x--, k = (kk + k - 1) % kk) 
	;
      Abuf[0] = ace_lower(b->base) ;  
      for (k = 1 ; k <= kk ; k++) 
	Abuf[k] = ace_upper (b[k].base) ; 
      Abuf[k] = 0 ;
      Bbuf[0] = ace_lower(b->base) ;  Bbuf[1] = 0 ;
      break ;
    }
  
  sprintf (buf, "%s:%d:%s:%s:%s"
	   , dictName (snp->targetDict, target)
	   , x
	   , type
	   , Abuf, Bbuf
	   ) ;
  
  return buf ;
} /* snpBRSname */

/*************************************************************************************/

static int snpBRS2tsfExportLoop (SNP *snp, ACEOUT ao, BigArray aa, BigArray aaP, int target) 
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  long ii, iiMax = arrayMax (aa) ;
  DICT *oldSnps = snp->oldSnps ;
  int minMutant = snp->minMutant ;
  int minCover = snp->minCover ;
  int minFrequency = snp->minFrequency ;
  char *snpNam, buf[1024] ;
  BitSet oldDone = 0 ;
  snpNam = buf ; 

  if (snp->BRS_count)
    oldDone = bitSetCreate (dictMax (oldSnps)+1, h) ;
  for (ii = 0 ; ii < iiMax ; ii++)
    {
      BRS *brs = bigArrp (aa, ii, BRS) ;
      if (brs && brs->type) /* type 0:wildtype, 1:sub, 2:ins, 3:del */
	{
	  int pos = brs->pos ;
	  int old ;

	  if (pos > 0 && pos < arrayMax (aaP))
	    {
	      BRS2 *brs2 = bigArrp (aaP, pos, BRS2) ;
	      int covp = brs->type == 1 ? brs2->sp : brs2->rp  ;
	      int covm = brs->type == 1 ? brs2->sm : brs2->rm  ;
	      int cov = covp + covm ;
	      int bb = brs->bp + brs->bm ;

	      if (snp->BRS_count)
		{
		  if (brs->type < 2 && ! brs2->oldSnp)
		    continue ;
		  snpBRSname (snp, target, brs, brs2, snpNam) ;
		  if (! dictFind (oldSnps, snpNam, &old))
		    continue ;
		  if (bitt (oldDone, old))
		    continue ;
		  bitSet (oldDone, old) ;
		}
	      else if (snp->BRS_detect)
		{
		  if (bb < 200 * minMutant ||
		      cov < 100 * minCover ||
		      100 * bb < minFrequency * cov
		      ) 
		    continue ;
		  snpBRSname (snp, target, brs, brs2, snpNam) ;
		}
	      else
		continue ;

	      nn++ ;
	      aceOutf (ao, "%s\t%s\t6i\t%d\t%d\t%d\t%d\t%d\t%d\n"
		       , snpNam
		       , snp->run
		       , covp/100
		       , brs->bp/100
		       , brs2->bp/100
		       , covm/100
		       , brs->bm/100
		       , brs2->bm/100
		       ) ;
	    }
	}
    }
  
  iiMax = oldSnps ? dictMax (oldSnps) : 0 ;
  if (iiMax > nn+1 )
    {
      const char *targetNam = dictName (snp->targetDict, target) ;
      int tLen = strlen (targetNam) ;
      char buf[tLen+12] ;
      
      sprintf (buf, "%s:%%d:%%c:", targetNam) ;
      for (ii = 1 ; ii <= iiMax ; ii++)
	{
	  if ( ! bitt (oldDone, ii))
	    {
	      int pos = 0, pos2 ;
	      char type = 0 ;
	      const char *snpNam = dictName (oldSnps, ii) ;
	      if (sscanf (snpNam, buf, &pos, &type) != 2)
		continue ;
	      if (pos < 1 || pos > bigArrayMax (aaP))
		continue ;
	      if (type == 'S')
		pos2 = pos ;
	      else if (type == 'I')
		pos2 = pos + 1 ;
	      else if (type == 'D')
		pos2 = pos + 1 ;
	      else
		continue ;
	      BRS2 *brs2 = bigArrp (aaP, pos, BRS2) ;
	      if (brs2->sp + brs2->sm > 0)
		aceOutf (ao, "%s\t%s\t6i\t%d\t%d\t%d\t%d\t%d\t%d\n"
			 , snpNam
			 , snp->run
			 , pos2 > pos ? brs2->rp/100 : brs2->sp/100
			 , 0
			 , brs2->bp/100
			 , pos2 > pos ? brs2->rm/100 : brs2->sm/100
			 , 0
			 , brs2->bm/100
			 ) ;
	    }
	}
    }

  ac_free (h) ;
  return nn ;
} /* snpBRS2tsfExportLoop */

/*************************************************************************************/

static int snpBRS2tsfExport (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nnOut = 0 ;
  int prefix = 6 ;   /* how many base we want to export on each side of the event in the SNPnet */
  Array aaa = snp->aaa ;  /* ee B array of the current target */
  Array aaaP = snp->aaaP ;  /* ee B array of the current target */
  Array ttt = snpBRS_sort_target (snp, h) ;
  ACEOUT ao = aceOutCreate (snp->outFileName
			    , snp->BRS_count ? ".snp_count.tsf" : ".snp_detect.tsf"
			    , snp->gzo
			    , h) ;

  aceOutDate (ao, "##", snp->title ? snp->title : "tsf file ") ; 
  snpBRS2snpsExportHeader (snp, ao, prefix) ;
    
  for (int ii = 0 ; ii < arrayMax (ttt) ; ii++)
    {
      int target = array (ttt, ii, int) ;
      BigArray aa = array (aaa, target, BigArray) ;
      BigArray aaP = array (aaaP, target, BigArray) ;
      if (aa && bigArrayMax (aa))
	nnOut += snpBRS2tsfExportLoop (snp, ao, aa, aaP, target) ;
    }

  
  ac_free (h) ;
  return nnOut ;
} /* snpBRS2tsfExport */

/*************************************************************************************/

static void snpBRS2tsf (SNP *snp)
{
  int nnIn = snpBRSparse (snp) ;
  int nnInc = snpBRScomplete (snp) ;
  int nnOut = snpBRS2tsfExport (snp) ;
  int mx = 0 ;
  
  messAllocMaxStatus (&mx) ;   
  fprintf (stderr, "... snpBRS2snps       \t%s\tmax memory %d Mb parsed %d:%d BRS lines eported %d snps\n", timeShowNow(), mx, nnIn, nnInc, nnOut) ;

  return ;
} /* snpBRS2tsf */

/*************************************************************************************/
/******************************* snpBRS2snps end *************************************/
/*************************************************************************************/
#ifdef JUNK
static void snpMerge (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = snp->ai ;
  ACEOUT ao = snp->ao ;
  DICT *runDict = dictHandleCreate (300, h) ;
  DICT *snpDict = dictHandleCreate (1000000, h) ;
  Stack snipnet = stackHandleCreate (1000000, h) ;
  SNR s, *sp ;
  Array aa = arrayHandleCreate (1000000, SNR, h) ;
  vTXT txt = vtxtHandleCreate (h) ;
  char *ss, *cp ;
  double z, z1, z2 ;
  float zchi2 = 0 ;
  int ii ;
  BOOL solid = FALSE ;


  aceOutf (ao, "## %s\n", timeShowNow ()) ;
  if (snp->title)
    aceOutf (ao, "## %s\n", snp->title) ;

  while (snpParseOneRecord (ai, &s, txt, runDict, snpDict, snipnet, &solid))
    {
      sp = arrayp (aa, s.sn, SNR) ; 
      sp->sn = s.sn ;
      sp->chrom = s.chrom ;
      sp->pos = s.pos ;
      sp->type = s.type ;
      if (s.sq1) sp->sq1 = s.sq1 ;
      if (s.sq2) sp->sq2 = s.sq2 ;
      sp->run = 0 ;
      sp->cover += s.cover ;
      sp->bad += s.bad ;
      sp->m += s.m ;
      sp->w += s.w ;
      sp->mp += s.mp ;
      sp->wp += s.wp ;
      sp->mm += s.mm ;
      sp->wm += s.wm ;
      sp->bad += s.bad ;
      sp->multip += s.multip ;
      sp->multim += s.multim ;
      sp->oo1p += s.oo1p ;
      sp->oo2p += s.oo2p ;
      sp->oo1m += s.oo1m ;
      sp->oo2m += s.oo2m ;
    }

      /*  2       41384390        G>A     CATGGAgGAAACT   CATGGAaGAAACT   Ghs103  2       100.0   14      14      0       3       0       11 00       0:0     0:0 */
    
  for (ii = 0 ; ii < arrayMax (aa) ; ii++)
    {
      sp = arrp (aa, ii, SNR) ;
      if (! sp->sn)
	continue ;
      if (sp->cover < snp->minCover)
	continue ;

      /* dosage */
      /* chi2 */
      z = sp->cover ? 100.0 * (sp->mp + sp->mm)/(double)sp->cover : 0 ;
      ss = cp = "-" ;

      if (z  < snp->minFrequency)
	continue ;
      aceOutf (ao, "%s\t%d\t%s\t%s\t%s\t%s"
	       , dictName (runDict, sp->chrom)
	       , sp->pos
	       , dictName (runDict, sp->type)
	       , sp->sq1 ? stackText (snipnet, sp->sq1) : "-"
	       , sp->sq2 ? stackText (snipnet, sp->sq2) : "-"
	       , snp->run ? snp->run : "-"
	       ) ;

      ss = "Both_strands" ; z2 = 0 ; zchi2 = 0 ;
      if (sp->mp + sp->wp + sp->mm + sp->wm < snp->minCover) ss = "-" ;
      else if (2 * (sp->mp + sp->wp) < snp->minCover) ss = "Minus_strand_only" ;
      else if (2 * (sp->mm + sp->wm) < snp->minCover) ss = "Plus_strand_only" ;
      else
	{ 
	  double x, y ;
	  chi2 (sp->mp, sp->wp, sp->mm, sp->wm, &zchi2) ;
	  x = sp->mp/(1.0 + sp->mp + (double)sp->wp) ;
	  y = sp->mm/(1.0 + sp->mm + (double)sp->wm) ;
	  if (z2 > 100 && (x > y + .5 || x < y - .5)) ss = "Incompatible_strands" ; /* extremely high  and clear contradiction */ 
	}

      snpNotLowShow (ao, sp->mp + sp->mm, sp->cover, &z, &z1, &z2, 0) ;

      aceOutf (ao, "\t%.1f", z) ;
      
      aceOutf (ao, "\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d"
	       , sp->cover
	       , sp->m, sp->w
	       , sp->mp, sp->wp
	       , sp->mm, sp->wm
	       , sp->bad
	       ) ;
      if (solid)
	aceOutf (ao, "\t%d:%d\t%d:%d"
		 , sp->oo1p, sp->oo2p
		 , sp->oo1m, sp->oo2m
		 ) ;
      else
	aceOutf (ao, "\t%d\t%d"
		 , sp->oo2p
		 , sp->oo2m
		 ) ;

      aceOutf (ao, "\t%s\t%.1f\t%s"
	       , ss, zchi2
	       , 100 * (sp->oo1p + sp->oo1m + sp->oo2p + sp->oo2m) > 40 * sp->cover ? "N_rich" : "-"
	       ) ;
      if (0) /* obsolete */
	aceOutf (ao, "\t%d\t%d\t%.1f\t%s"
	       , sp->multip, sp->multim
	       , z2
	       , 100 * (sp->multip + sp->multim) > 5 * sp->cover ? "REPEATED" : "-"
	       ) ;

      snpNotLowShow (ao, sp->mp, sp->mp + sp->wp, &z, &z1, &z2, 0) ;
      snpNotLowShow (ao, sp->mm, sp->mm + sp->wm, &z, &z1, &z2, 0) ;
      aceOutf (ao, "\n") ;
    }

  ac_free (h) ;
} /* snpMerge */
#endif
/*************************************************************************************/
typedef struct runInfoStruct { int machine, sample, system1, system2, tissue1, tissue2 ; } RI ;

static int snpMergeExportPopulationTable (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = snp->ai ;
  ACEOUT ao = 0, aos[12] ;
  DICT *runDict = dictHandleCreate (300, h) ;
  DICT *snpDict = dictHandleCreate (1000000, h) ;
  DICT *infoDict = dictHandleCreate (3000, h) ;
  Stack snipnet = stackHandleCreate (1000000, h) ;
  SNR s, *sp, *sq ;
  Array aa = arrayHandleCreate (1000000, SNR, h) ;
  vTXT txt = vtxtHandleCreate (h) ;
  int nn = 0, aaMax, ii, run, oldSn = 0, iType ;
  BOOL solid = FALSE ;
  KEYSET rM = keySetHandleCreate (h) ;
  char *outNam[] = {"genotypes", "counts", 0 } ;
  int nHomo,  nHetero, nWild ;
  int mp, mm, wp, wm ;
  const char *ccp  ;
  RI *info ;
  Array infos = arrayHandleCreate (300, RI, h) ;

  /* PREORDER THE runs according to the complementary info */
  aceInSpecial (ai, "\t\"\n") ;
  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (! ccp) continue ;
      if (! strcmp (ccp, "ZZZZZ"))
	break ;
      dictAdd (runDict, ccp, &run) ;
      info = arrayp (infos, run, RI) ;
      
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp) continue ;
      dictAdd (infoDict, ccp, &(info->machine)) ;
   
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp) continue ;
      dictAdd (infoDict, ccp, &(info->sample)) ;
      
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp) continue ;
      dictAdd (infoDict, ccp, &(info->system1)) ;
      
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp) continue ;
      dictAdd (infoDict, ccp, &(info->tissue1)) ;
      
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp) continue ;
      dictAdd (infoDict, ccp, &(info->tissue2)) ;
      
      aceInStep (ai, '\t') ;
      ccp = aceInWord (ai) ;
      if (! ccp) continue ;
      dictAdd (infoDict, ccp, &(info->system2)) ;
    }      

  while (snpParseOneRecord (ai, &s, txt, runDict, snpDict, snipnet, &solid))
    {
      sp = arrayp (aa, nn++, SNR) ; 
      *sp = s ;
      if (sp->dose > 0)
	keySet (rM, sp->run) = 1;
    }
  arraySort (aa, snrOrder) ;

  for (iType = 0 ; outNam[iType] ; iType++)
    {
      ao = aos[iType] = aceOutCreate (snp->outFileName, messprintf(".%s.txt", outNam[iType]), snp->gzo, h) ;
      aceOutf (ao, "## %s\n", timeShowNow ()) ;
      if (snp->title)
	aceOutf (ao, "## %s\n", snp->title) ;

      aceOutf (ao, "## Distribution of %s of the variants observed in the population\n"
	       , outNam[iType]) ;
      aceOutf (ao, "\n# Chromosome\tPosition\tType\tReference\tVariant\tHomozygous reference\tHeterozygous\tHomozygous variant\tFrequency\tCover\tVariant\tReference\tVar strand+\tRef strand+\tVar strand-\tRef strand-") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	if (keySet (rM, run))
	  aceOutf (ao, "\t%s", dictName (runDict, run)) ;
      aceOutf (ao, "\n# Machine\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	if (keySet (rM, run))
	  {
	    info = arrayp (infos, run, RI) ;
	    ccp = info->machine ? dictName (infoDict, info->machine) : "-" ;
	    aceOutf (ao, "\t%s", ccp) ;
	  }
      aceOutf (ao, "\n# Sample\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	if (keySet (rM, run))
	  {
	    info = arrayp (infos, run, RI) ;
	    ccp = info->sample ? dictName (infoDict, info->sample) : "-" ;
	    aceOutf (ao, "\t%s", ccp) ;
	  }
      aceOutf (ao, "\n# System\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	if (keySet (rM, run))
	  {
	    info = arrayp (infos, run, RI) ;
	    ccp = info->system1 ? dictName (infoDict, info->system1) : "-" ;
	    aceOutf (ao, "\t%s", ccp) ;
	  }
      aceOutf (ao, "\n# Subsystem\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	if (keySet (rM, run))
	  {
	    info = arrayp (infos, run, RI) ;
	    ccp = info->system2 ? dictName (infoDict, info->system2) : "-" ;
	    aceOutf (ao, "\t%s", ccp) ;
	  }
      aceOutf (ao, "\n# Tissue\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	if (keySet (rM, run))
	  {
	    info = arrayp (infos, run, RI) ;
	    ccp = info->tissue1 ? dictName (infoDict, info->tissue1) : "-" ;
	    aceOutf (ao, "\t%s", ccp) ;
	  }
      aceOutf (ao, "\n# Sub tissue\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	if (keySet (rM, run))
	  {
	    info = arrayp (infos, run, RI) ;
	    ccp = info->tissue2 ? dictName (infoDict, info->tissue2) : "-" ;
	    aceOutf (ao, "\t%s", ccp) ;
	  }
      aceOutf (ao, "\n# Run\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t") ;
      for (run = 1 ; run <= dictMax (runDict) ; run++)
	if (keySet (rM, run))
	  {
	    ccp = dictName (runDict, run) ;
	    aceOutf (ao, "\t%s", ccp) ;
	  }
    }

  aaMax = arrayMax (aa) ;
  sp = arrayp (aa, aaMax, SNR) ; /* guarantee a null record at the end */
  for (ii = 0 ; ii < aaMax ; ii++)
    {
      sp = arrp (aa, ii, SNR) ;
      if (! sp->sn || sp->sn == oldSn ) continue ;
      oldSn = sp->sn ;
      /* prevalence in the population */
      nHomo = nHetero = nWild = 0 ;
      mp = mm = wp = wm = 0 ;
      for (sq = sp ; sq->sn == oldSn ; sq++)
	if (sq->run && keySet (rM, sq->run))
	  {
	    if (0 && sq->dose > .4 && sq->dose < 1.7) sq->dose = 1 ;
	    if (sq->dose == 2) nHomo++ ;
	    else if (sq->dose == 0) nWild++ ;
	    else if (sq->dose < 0) ;
	    else nHetero++ ;

	    mp += sq->mp ; wp += sq->wp ;
	    mm += sq->mm ; wm += sq->wm ;	    
	  }

      /* identify the variant */ 
      for (iType = 0 ; iType < 3 && outNam[iType] ; iType++)
	{
	  int cover = mp + mm + wp + wm ;
	  int coverp = mp + wp ;
	  int coverm = mm + wm ;

	  ao = aos[iType] ;
	  if (cover == 0) cover = 1 ;
	  if (coverp == 0) coverp = 1 ;
	  if (coverm == 0) coverm = 1 ;

	  aceOutf (ao, "\n%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.1f\t%d\t%d\t%d\t%d\t%d\t%d\t%d"
		   , dictName (runDict, sp->chrom)
		   , sp->pos
		   , dictName (runDict, sp->type)
		   , sp->sq1 ? stackText (snipnet, sp->sq1) : "-"
		   , sp->sq2 ? stackText (snipnet, sp->sq2) : "-"
		   , nWild, nHetero, nHomo
		   , 100.0 * (mp+mm)/cover, 100.0 * mp/coverp , 100.0 * mm/coverm, cover, mp+mm, wp+wm
		   , mp,wp,mm,wm
		   ) ;
	}

      for (run = 1, sq = sp ; run <= dictMax (runDict) ; run++)
	if (keySet (rM, run))
	  {
	    BOOL ok1 = FALSE ;
	    /* do not restart at sq = sp: rely on the idea that aa is ordered by sp->run */
	    for (sq = sp ; ! ok1 && sq->sn == oldSn  ; sq++)
	      if (sq->run == run)
		{
		  for (iType = 0 ; iType < 2 && outNam[iType] ; iType++)
		    {
		      ao = aos[iType] ;
		      switch (iType)
			{
			case 0: 
			  if (sq->dose == 0) aceOutf (ao, "\t0") ;
			  else if (sq->dose == 1) aceOutf (ao, "\t1") ;
			  else if (sq->dose == 2) aceOutf (ao, "\t2") ;
			  else if (sq->dose < 0) aceOutf (ao, "\t") ;
			  else aceOutf (ao, "\t%.1f", sq->dose) ; 
			  break ;
			case 1:
			  aceOutf (ao, "\t%d\t%d\t%d\t%d"
				   , sq->mp, sq->wp
				   , sq->mm, sq->wm
				   ) ;
			  break ;
			}	
		    }
		  ok1 = TRUE ;
		}
	    if (! ok1)	
	      {
		for (iType = 0 ; iType < 2 && outNam[iType] ; iType++)
		  {
		    ao = aos[iType] ;
		    aceOutf (ao, "\t") ;
		  }
	      }
	  }
    }
  aceOutf (ao, "\n") ;
  for (iType = 0 ; iType < 3 && outNam[iType] ; iType++)
    {
      ao = aos[iType] ;
      if (ao) aceOutf (ao, "\n") ;
      ac_free (aos[iType]) ;
    }

  return nn ;
} /* snpMergeExportPopulationTable */

/*************************************************************************************/
/*************************************************************************************/
/* Get the ordered PREORDER THE runs according to the complementary info */
static  AC_TABLE snpProjectRuns (SNP *snp,  AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_KEYSET ks = 0 ;
  AC_TABLE tt = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;
  const char *errors = 0 ;
  char *cp ;

  cp =  hprintf (h, "Find run project == \"%s\" AND ((union_of && SNP && add_counts) OR (sublibraries) OR (Is_run && ! Sublibrary_of))",  snp->project) ;
  ks =  ac_dbquery_keyset (snp->db, cp, h) ;
  if (! ks || ! ac_keyset_count (ks))
    messcrash (" snpProjectRuns no run found using query=%s", cp) ;

  vtxtPrintf (txt, "table = ") ;
  vtxtPrintf (txt, "colonne 1 \n class Run \n\n") ;
  vtxtPrintf (txt, "colonne 2 \nFrom 1\nClass Text\nTag Sorting_title\n\n") ;
  vtxtPrintf (txt, "colonne 3 \nFrom 1\nClass Text\nTag Title\n\n") ;
  vtxtPrintf (txt, "colonne 4 \nFrom 1\nClass Sample\nTag Sample\n\n") ;

  tt = ac_tablemaker_table (snp->db, vtxtPtr (txt), ks, ac_tablemaker_text , 0 , "2+1", &errors, h0) ;

  if (errors)
    messcrash (" snpProjectRuns table maker error = %s\nerr = %s", vtxtPtr (txt), errors) ;
  if (! tt)
    messcrash (" snpProjectRuns no table found using table =%s", vtxtPtr (txt)) ;
  if (! tt->rows)
    messcrash (" snpProjectRuns no run found using table = %s", vtxtPtr (txt)) ;
  fprintf (stderr, "snpProjectRuns  found %d runs\n", tt ? tt->rows : 0) ;
  ac_free (h) ;
  
  return tt ;
} /* snpProjectRuns */

/*************************************************************************************/
/* cumulate in each variant, for each group with SNP_RNA/DNA tag the global counts
 * per gene dose and per strand according to the notes in wspec/models.wrm 
 */
static void snpDbCount (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = 0 ;
  const char *errors = 0 ;
  char oldNam[1024], newNam[1024] ;
  vTXT txt = vtxtHandleCreate (h) ;
  
  if (! (ai = aceInCreate (snp->inFileName, snp->gzi, h)))
    messcrash ("-db_count cannot open file : %s\n", snp->inFileName) ;
 
  memset (newNam, 0, sizeof(newNam)) ;
  memset (oldNam, 0, sizeof(oldNam)) ;

  /* scan the file,
   * reject unkown variants
   * cumulate all counts related to the same snp
   * store in the database
   */
  while (aceInCard (ai))
    {
      /* recognize the snp name */
      memset (newNam, 0, sizeof(newNam)) ;
      
      if (strcmp (oldNam, newNam))
	{   /* store the previous counts */
	  if (vtxtPtr (txt))
	    {
	      ac_parse (snp->db, vtxtPtr (txt), &errors, 0, h) ; 
	      vtxtClear (txt) ;
	    }
	}
      strcpy (oldNam, newNam) ;
    }

  
  if (vtxtPtr (txt)) /* store the last results */
    ac_parse (snp->db, vtxtPtr (txt), &errors, 0, h) ; 

  ac_free (h) ;
  return ;
} /* snpDbCount */

/*************************************************************************************/

static BOOL snpPleaseDropMonomodal (int nww, int nlow, int nwm, int nhigh, int npure)
{
  int nm = nww + nlow + nwm + nhigh + npure ;
  if (10 * nlow > nm && nlow > 5 &&
      ! (nww > nlow + 10 && nwm + npure > nlow + 10) &&
      ! (npure > nhigh)
      )
    return TRUE ;
  return FALSE ;
} /* snpDropMonomodal */

/*************************************************************************************/

static BOOL snpFrequencyTableExportOneLine (SNP *snp, ACEOUT aoF, ACEOUT aoC
					    , int max
					    , Array ff, Array ff1, Array ff2
					    , int oldTarget, int oldPos, int oldType
					    , AC_TABLE groups, AC_TABLE groups2runs
					    , KEYSET r2ir, DICT *runDict
					    , int mp, int mm, int wp, int wm
					    )
{
  int i, ir, jr, pass ;
  float z, zMax ;
  KEY key ;
  int nm = 0, nww = 0, nwm = 0, npure = 0, nlow = 0, nhigh = 0, nNA = 0 ;
  BOOL isDrop = FALSE ;
  ACEOUT ao = 0 ;

  for (i = 0 ; i < max ; i++)
    {
      z = array (ff, i, float) - 1000 ;
      nm++ ; /* measured */
      if (z < 0) { nNA++ ; nm-- ; }
      else if (z < 5) nww++ ; /* wild */
      else if (z < 20) nlow++ ;
      else if (z < 80) nwm++ ;
      else if (z < 95) nhigh++ ;
      else npure++ ;
    }
  if (nlow + nwm + nhigh + npure == 0)
    return TRUE ;
  if (nwm + nhigh + npure == 0)
    isDrop = TRUE ;
  else if (snp->dropMonomodal)
    isDrop = snpPleaseDropMonomodal (nww, nlow, nwm, nhigh, npure) ;

  if (! isDrop)
    { /* evaluate the anyRun statistics */
      snp->anyRunMeasured++ ;
      if (npure) snp->anyRunPure++ ;
      else if (nhigh) snp->anyRunHigh++ ;
      else if (nwm) snp->anyRunMid++ ;
      else if (nlow) snp->anyRunLow++ ;
      else if (nww) 
	snp->anyRunWild++ ;
      else if (nNA) snp->anyRunNA++ ;    
     }
  else
    {
      snp->anyRunDropped++ ;
    }

  for (pass = 0 ; pass < 2 ; pass++)
    {
      ao = pass == 0 ? aoF : aoC ;
      aceOutf (ao, "\t%d\t%d\t%d\t%d", mp, mm, wp, wm) ;
      aceOutf (ao, "\t%s%d", (isDrop ? "+" : ""), chi2 (mp, mm, wp, wm, &z)) ; /* type : isDrop + sign allows sorting but does not affect excell */
      aceOutf (ao, "\t%.1f", z) ;
      aceOutf (ao, "\t%d\t%d\t%d\t%d\t%d\t%d\t%d", nm, nNA, nww, nlow, nwm, nhigh, npure) ;
      if (1) /* allelic frequency in the population (before 2020_01_17) */
	aceOutf (ao, "\t%.2f", nm > 20 ? (.2 * nlow + 1.0 * nwm + 1.8 * nhigh + 2 * npure )/ (0.02 * nm) : -10.0) ;
      if (1) /* average frequency 2020_01_17 */
	aceOutf (ao, "\t%.2f", mp + mm + wp + wm >= 20 ? (mp + mm)/ (0.01 * (mp + mm + wp + wm)) : -20.0) ;
      aceOutf (ao, "\t") ;
      for (ir = jr = 0 ; groups && ir < groups->rows ; ir++)
	{
	  nm = nNA = nww = nlow = nwm = nhigh = npure = 0 ;
	  key = ac_table_key (groups, ir, 0, 0) ;
	  for (jr = 0 ; key && groups2runs && jr < groups2runs->rows ; jr++)
	    {
	      int run ;
	      if (ac_table_key (groups2runs, jr, 0, 0) == key && dictFind (runDict, ac_table_printable (groups2runs, jr, 1, "-"), &run))
		{
		  i =  keySet(r2ir, run) ;
		  z = array (ff, i, float) - 1000 ;
		  nm++ ; /* measured */
		  if (z < 0) { nNA++ ; nm-- ; }
		  else if (z < 5) nww++ ; /* wild */
		  else if (z < 20) nlow++ ;
		  else if (z < 80) nwm++ ;
		  else if (z < 95) nhigh++ ;
		  else npure++ ;
		}
	    }
	  aceOutf (ao, "\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t", nm, nNA, nww, nlow, nwm, nhigh, npure, nm > 0 ? (.2 * nlow + 1.0 * nwm + 1.8 * nhigh + 2 * npure )/ (0.02 * nm) : 0) ;
	}
      zMax = -10 ;
      for (i = 0 ; i < max ; i++)
	{
	  z = array (ff, i, float) - 1000 ; 
	  if (z > zMax) zMax = z ;
	}
      aceOutf (ao, "\t") ;
      aceOutPercent (ao, zMax) ;
      for (i = 0 ; i < max ; i++)
	{
	  if (pass == 0)
	    {
	      z = array (ff, i, float) - 1000 ; 
	      aceOutf (ao, "\t") ;
	      if (z  >= 0 ) aceOutPercent (ao, z) ;
	      else  aceOutf (ao, "-%d", snp->minCover) ;
	    }
	  else
	    {
	      int z1 = array (ff1, i, int) ;
	      int z2 = array (ff2, i, int) ;
	      aceOutf (ao, "\t") ;
	      if (z2  >= 0 ) aceOutf(ao, "%d|%d", z1, z2) ;
	    }
	}
    }
  return isDrop ;
}  /* snpFrequencyTableExportOneLine */

/*************************************************************************************/

static void snpPrettyCoordinatesName (vTXT txt, int a1, int a2, char *type, BOOL isRNA, BOOL isSliding)
{
  int ddx =  strlen (type) - 3 ;
  
  if (a2 > a1 || a2 == 1) a2 = 1 ;
  else a2 = -1 ;

  if (ddx == 0 && type[1] == '>')
    vtxtPrintf (txt, "%d%s", a1, type) ;
  else if (! strncasecmp (type,  "Ins", 3))
    {
      if (isSliding)
	{ 
	  if (ddx == 1)
	    vtxtPrintf (txt,"%ddup%s", a1 - 1, type + 3) ;
	  else
	    vtxtPrintf (txt,"%d_%ddup%s", a1 - ddx, a1 - 1, type + 3) ;
	}
      else
	vtxtPrintf (txt,"%d_%dins%s", a1 - 1, a1, type + 3) ;
    }
  else if (! strncasecmp (type,  "Del", 3))
    {
      char *t = isSliding ? "dim" : "del" ;
      if (ddx == 1)
	vtxtPrintf (txt,"%d%s%s", a1, t, type + 3) ;
      else
	vtxtPrintf (txt,"%d_%d%s%s", a1, a1 + ddx - 1, t, type + 3) ;
    }
} /* snpPrettyCoordinatesName */

/*************************************************************************************/
/* cosntruct the 3 standard names of a snp 
 * return exonic | protein_changing : 0x1 | 0x2 
 */
static int snpPrettyNames (SNP *snp, AC_OBJ snpObj
			    , vTXT gTxt, vTXT cTxt, vTXT pTxt
			    , int *gMapp, int *gPosp
			    , char *gType, char *cType, vTXT pType
			    , vTXT location
			    , vTXT geneboxes, vTXT avgeneboxes, vTXT gsnippet, vTXT snippet, vTXT pSnippet
			   )
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE imap = 0, tt1 = 0, tt2 = 0,  ttg1 = 0, ttg2 = 0, tt3 = 0,  ttm = 0 ;
  int ir, a1, a2, ga1, ga2 = 0, type = 0 ;
  const char *ccp, *product ;
  BOOL isIntronic = TRUE ;
  BOOL isSliding = FALSE ;
  KEY key ;

  vtxtClear (gTxt) ;
  vtxtClear (cTxt) ;
  vtxtClear (pTxt) ;
  vtxtClear (location) ;
  vtxtClear (avgeneboxes) ;
  vtxtClear (geneboxes) ;
  vtxtClear (gsnippet) ;
  vtxtClear (snippet) ;
  vtxtClear (pSnippet) ;
  cType[0] = 0 ;
  vtxtClear (pType) ;
 
  strcpy (cType, gType) ; /* as given in the snp object */

  /* genome name */ 
  *gMapp = 0 ; *gPosp = 0 ;
  imap = ac_tag_table (snpObj, "IntMap", h) ;
  if (imap && imap->cols >= 3)
    {
      int x ;
      ga1 = ac_table_int (imap, 0, 1, 0) ;
      ga2 = ac_table_int (imap, 0, 2, 0) ;
      dictAdd (snp->targetDict, ac_table_printable (imap, 0, 0, "-"), gMapp) ;
      *gPosp =  ga1 ;
      if (1)
	{
	  AC_OBJ map = ac_table_obj (imap, 0, 0, h) ;
	  BOOL hasTitle = FALSE ;

	  ccp = 0 ;
	  if (map)
	    ccp = ac_tag_printable (map, "Title", 0) ;
	  if (ccp) 
	    {
	      hasTitle = TRUE ;
	      vtxtPrintf (gTxt, "%s", ccp) ;
	    }
	  else
	    {
	      ccp = ac_table_printable (imap, 0, 0, "-") ;
	      
	      vtxtPrintf (gTxt, "%s", ccp) ;
	    }

	  vtxtPrintf (gTxt, ":g.", ccp) ;
	  
	  /* in-repeat applies to the genomic repeat */
	  isSliding = ac_has_tag (snpObj, "In_repeat") ;
	  snpPrettyCoordinatesName (gTxt, ga1, ga2, gType, FALSE, isSliding) ;

	  if (hasTitle || snp->Reference_genome)
	    {
	      vtxtPrintf (gTxt, "(") ;
	      if (hasTitle)
		vtxtPrintf (gTxt, "%s", ac_name(map)) ;
	      if (hasTitle && snp->Reference_genome) 
		vtxtPrintf (gTxt, ",") ;
	      if (snp->Reference_genome)
		vtxtPrintf (gTxt, "%s", snp->Reference_genome) ;
	      vtxtPrintf (gTxt, ")") ;
	    }
	}
     /*
	Potential_splice_disruption Near_donor    // -16, -1, 1,2 3, 5,, ne pas remplir ces tags detailles, pvalue < 1/100, and > 20% splice disrupting (2 stars), add spl at the end of the c.name
	Near_acceptor //  -2, -1  Rivas et al., Effect of predicted protein truncating SNVs on transcriptome. Science 348:666-669, May 2015.
      */
      if ((x = ac_tag_int (snpObj, "Near_donor", 0)))
	vtxtPrintf (gTxt, ".spl(?donor%d)", x) ;
      if ((x =ac_tag_int (snpObj, "Near_acceptor",0)))
	vtxtPrintf (gTxt, ".spl(?donor%d)", x) ;
    }

  /* cDNA name */
  ttg1 = ac_tag_table (snpObj, "Reference_genomic_sequence", h) ;
  ttg2 = ac_tag_table (snpObj, "Observed__genomic_sequence", h) ;
  tt1 = ac_tag_table (snpObj, "Reference_sequence", h) ;
  tt2 = ac_tag_table (snpObj, "Observed__sequence", h) ;

  vtxtPrintf (location, "Intergenic") ;
  ttm = ac_tag_table (snpObj, "GeneBox", h) ;
  if (ttm)
    {
      int x = -999999, y ;
      for (ir = 0 ; ir < ttm->rows ; ir++)
	{
	  vtxtPrintf (geneboxes, "%s%s", ir ? ", " : "", ac_table_printable (ttm, ir, 0, "")) ;
	  y = ac_table_int (ttm, ir, 1, -999999) ;
	  if (y > x) x = y ;
	}
      if (x > 0) 
	{  
	  isIntronic = TRUE ;
	  vtxtClear (location) ; vtxtPrintf (location, "Intronic") ; 
	}
      else if (x < 0 && x > -99999)
	{ vtxtClear (location) ; vtxtPrintf (location, "Promotor region") ; }
    }
  ttm = ac_tag_table (snpObj, "mRNA", h) ;
  if (ttm)
    {
      BOOL first = TRUE ;
       for (ir = 0 ; ir < ttm->rows && ir < 1 ; ir++)
	 {
	   AC_OBJ mrna = ac_table_obj (ttm, ir, 0, h) ;
	   if (mrna)
	     {
	       char *cq = vtxtPtr (geneboxes) ;
	       ccp = ac_tag_printable (mrna, "Gene", 0) ;

	       if (ccp && (! cq || ! strstr (cq, ccp)))
		 {
		   if (first)
		     {
		       first = FALSE ;
		        vtxtClear (geneboxes) ;
			vtxtPrintf (geneboxes, "%s", ccp) ;
		     }
		   else
		     vtxtPrintf (geneboxes, ", %s", ccp) ;
		   if (1)
		     {
		       AC_OBJ gene = ac_tag_obj (mrna, "Gene", h) ;
		       if (gene)
			 {
			   ccp = ac_tag_printable (gene, "Title", 0) ;
			   if (ccp)
			     vtxtPrintf (geneboxes, ": %s", ccp) ;
			   ac_free (gene) ;
			 }
		     }
		 }
	       ac_free (mrna) ;
	     }
	 }
    }

  vtxtClear (gsnippet) ;
  if (ttg1 && ttg2)  /* by default use the gsnippet in the genome orientation */
    {
      const char *ccp = ac_table_printable (ttg1, 0, 0, 0) ;
      if (ccp) 
	{
	  vtxtPrintf (gsnippet, "%s", ccp) ;
	  ccp = ac_table_printable (ttg2, 0, 0, 0) ;
	  if (ccp) vtxtPrintf (gsnippet," > %s",  ccp) ;
	}
    }
  vtxtClear (snippet) ;
  if (tt1 && tt2)  /* by default use the snippet in the given orientation */
    {
      const char *ccr ;
      const char *ccp = ac_table_printable (tt1, 0, 0, "") ; 
      const char *ccq = ac_table_printable (tt2, 0, 0, "") ;
      int i, k1 = strlen (ccp) ;
      int k2 = strlen (ccq) ;
      char *cp, buf1 [k1+7] ;
      char *cq, buf2 [k2+7] ;
      int ddx = strlen (cType) - 3 ;

      if (ace_lower (cType[0]) == 'd')
	ddx *= -1 ;
      buf1[0] = buf2[0] = 0 ;
      if (ga2 >= 0)
	{
	  ccr = ccp ;
	  i = 0 ;
	  while (*ccr)
	    {
	      char c = rnaDecodeChar[(int)dnaEncodeChar[(int)*ccr]] ;
	      buf1[i++] = c ; ccr++ ;
	    }
	  buf1[i++] = 0 ;
	  ccr = ccq ;
	  i = 0 ;
	  while (*ccr)
	    {
	      char c = rnaDecodeChar[(int)dnaEncodeChar[(int)*ccr]] ;
	      buf2[i++] = c ; ccr++ ;
	    }
	  buf2[i++] = 0 ;
	}
      else  /* complement to obtain the g orientation */
	{
	  ccr = ccp + strlen (ccp)  ;
	  i = 0 ;
	  while (--ccr >= ccp)
	    {
	      char c = rnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)*ccr])] ;
	      buf1[i++] = c ;
	    }
	  buf1[i++] = 0 ;

	  ccr = ccq + strlen (ccq)  ;
	  i = 0 ;
	  while (--ccr >= ccq)
	    {
	      char c = rnaDecodeChar[(int)complementBase(dnaEncodeChar[(int)*ccr])] ;
	      buf2[i++] = c ;
	    }
	  buf2[i++] = 0 ;
	}	  
      cp = buf1 ; cq = buf2 ;
      while (ace_lower(*cp) == ace_lower(*cq))
	{
	  *cp = ace_upper(*cp) ;
	  *cq = ace_upper(*cq) ;
	  cp++ ; cq++ ;
	}
      if (ddx == 0)
	{
	  *cp = ace_lower(*cp) ;
	  *cq = ace_lower(*cq) ;
	  cp++ ; cq++ ;
	}
      else if (ddx > 0)
	{
	  for (i = 0 ; i < ddx ; i++)
	    { *cq = ace_lower(*cq) ; cq++ ; }
	}
      else if (ddx < 0)
	{
	  for (i = 0 ; i < -ddx ; i++)
	    { *cp = ace_lower(*cp) ; cp++ ; }
	}
      while (*cp)
	{
	  *cp = ace_upper(*cp) ;
	  cp++ ; 
	}
      while (*cq)
	{
	  *cq = ace_upper(*cq) ;
	  cq++ ;
	}
      vtxtPrintf (snippet,"%s > %s",  buf1, buf2) ;
      isSliding = snpIsSliding (cType, buf1, buf2, TRUE, TRUE, 0, 0) ;
    }
  
  if (ttm && ttm->rows > 0)
    {
      int ddx = strlen (cType) - 3 ;
      vtxtClear (location) ; vtxtPrintf (location, "Exonic") ;
      type |= 0x1 ;
      if (tt2 && tt2->rows >= 1 && tt2->cols >= 3)
	{
	  key = ac_table_key (tt2, 0, 2, 0) ;
	  if (key)
	    for (ir = 0 ; ir < ttm->rows ; ir++)
	      {
		if ( key == ac_table_key (ttm ,ir, 0, 0))
		  {  
		    a1 =  ac_table_int (ttm,ir, 1, 0) ;
		    a2 =  ac_table_int (ttm,ir, 2, 1) ;
		    vtxtPrintf (cTxt, "%s:r.", ac_table_printable (tt2, 0, 2, "")) ;
		    snpPrettyCoordinatesName (cTxt, a1, a2, cType, TRUE, isSliding) ;
		    /* translated snippet */
		    ccp = ac_table_printable (tt2, 0, 1, 0) ;
		    if (ccp)
		      {
			char *cp, *cq, *ccp2 = hprintf (h, "%s", ccp) ;
			char *ccp1 = hprintf (h, "%s", ac_table_printable (tt1, 0, 1, 0)) ;
		
			/* fix the upper lower case */
			cp = ccp1 ; cq = ccp2 ;
			while (*cp && *cq && ace_upper (*cp) == ace_upper(*cq))
			  { 
			    *cp = ace_upper (*cp) ; cp++ ;
			    *cq = ace_upper (*cq) ; cq++ ;
			  }
			if (ddx == 3)
			  {
			    if (!strncmp (cType, "Del", 3)) 
			      { *cp = ace_lower (*cp) ; cp++ ; }
			    if (!strncmp (cType, "Ins", 3)) 
			      { *cq = ace_lower (*cq) ; cq++ ; }
			    while (*cp && *cq && ace_upper (*cp) == ace_upper(*cq))
			      { 
				*cp = ace_upper (*cp) ; cp++ ;
				*cq = ace_upper (*cq) ; cq++ ;
			      }
			  }
			else
			  {  
			    while (*cp)
			      {
				*cp = ace_upper (*cp) ; cp++ ;
			      }
			    while (*cq)
			      {
				*cq = ace_upper (*cq) ; cq++ ;
			      }
			  }
			cq = strchr (ccp1, 'x') ; if (cq) cq[1] = 0 ;
			cq = strchr (ccp2, 'x') ; if (cq) cq[1] = 0 ;
			vtxtPrintf (pSnippet, "%s > %s", ccp1, ccp2) ;
		      }
		    else
		      vtxtPrintf (pSnippet, "####") ;
		    break ;
		  }
	      }
	}
      else if (ttm->rows >= 1)
	{
	  vtxtPrintf (cTxt, "%s:r.", ac_table_printable (ttm,0, 0, "?")) ; 
	  a1 =  ac_table_int (ttm,0, 1, 0) ;
	  a2 =  ac_table_int (ttm,0, 2, 1) ;
	  snpPrettyCoordinatesName (cTxt, a1, a2, cType, TRUE, isSliding) ;
	}
    }

  /* protein name */
  product = tt2 ? ac_table_printable (tt2, 0, 2, 0) : 0 ;
  if (product)
    { 
      int xx = ac_table_int (tt2, 0, 4, -1) ;
      if (xx >= 1)
	{
	  if ((tt3 = ac_tag_table (snpObj, "Synonymous", h)))
	    { 
	      const char *ccp = ac_table_printable (tt3, 0, 0, 0) ;
	      
	      vtxtClear (location) ; vtxtPrintf (location, "Coding"); 	 

	      vtxtClear (pTxt) ; vtxtPrintf (pTxt, "p.%s", ccp ? ccp : "=") ;
	      vtxtClear (pType) ;
	      vtxtPrint (pType, "Synonymous") ;
	      if (ccp)
		{
		  char buf[4] ;
		  memcpy (buf, ccp, 3) ; buf[3] = 0 ;
		  vtxtPrintf (pType, " %s", ccp) ;
		}
	    }
	  
	  else if ((tt3 = ac_tag_table (snpObj, "Length_variation", h))) 
	    {
	      int ir ;
	      
	      type |= 0x2 ;
	      for (ir = 0 ; ir < tt3->rows ; ir++)
		{
		  const char *ccp = ac_table_printable (tt3, ir, 0, 0) ;
		  
		  if (ccp && strcasecmp (ac_table_printable (tt3, ir, 0, ""), "Extension"))
		    {
		      vtxtClear (pType) ; vtxtPrint (pType, ccp) ;
		      vtxtClear (location) ; vtxtPrintf (location, "Coding");
		      ccp = ac_table_printable (tt3, ir, 1, 0) ;
		      if (ccp) vtxtPrintf (pTxt, "p.%s", ccp) ;
		      break ;
		    }
		}
	    }
	  
	  else if ((tt3 = ac_tag_table (snpObj, "AA_substitution", h))) 
	    { 
	      type |= 0x2 ;
	      vtxtClear (location) ; vtxtPrintf (location, "Coding") ; 
	      ccp = ac_table_printable (tt3, 0, 0, "") ;
	      vtxtClear (pType) ; vtxtPrintf (pType, "AA_sub %s", ccp) ;
	      vtxtPrintf (pTxt, "p.%s", ccp) ;
	    } 
	  else
	    {
	      vtxtClear (location) ; vtxtPrintf (location, "Coding");
	      vtxtClear (pType) ; vtxtPrint (pType, ac_table_printable (tt3, 0, 1, "-")) ;
	    }
	}
      else
	{
	  if (ac_has_tag (snpObj, "Met1_gained"))
	    {
	      ccp = ac_tag_printable (snpObj, "Met1_gained", 0) ;
	      vtxtClear (location) ; vtxtPrintf (location, "Coding") ; 
	      if (ccp)
		{
		  vtxtClear (pType) ; vtxtPrintf (pType, "Met1_gained %s", ccp) ;
		  vtxtPrintf (pTxt, "p.%s", ccp) ;
		}
	    }
	  else
	    { 
	      vtxtClear (pType) ; vtxtPrintf (pType, "5' UTR %s", cType) ; 
	      vtxtClear (location) ; vtxtPrintf (location, "5' UTR") ;	
	    }
	}
    }
  else if (ac_has_tag (snpObj, "Non_coding_transcript")) 
    {  
      vtxtClear (pType) ; vtxtPrintf (pType, "Non_coding_transcript %s", cType) ; 
      vtxtClear (location) ;vtxtPrintf (location, "Non_coding_transcript") ; 
    }
  else if (ac_has_tag (snpObj, "UTR_5prime")) 
    { 
      vtxtClear (pType) ; vtxtPrintf (pType, "5' UTR %s", cType) ;
      vtxtClear (location) ;vtxtPrintf (location, "5' UTR") ; 
    }
  else if (ac_has_tag (snpObj, "UTR_3prime"))
    { 
      vtxtClear (pType) ; vtxtPrintf (pType, "3' UTR %s", cType) ; 
      vtxtClear (location) ; vtxtPrintf (location, "3' UTR") ; 
    }
  else if (isIntronic)
    { 
      vtxtClear (pType) ; vtxtPrintf (pType, "Intronic") ;
    }
  
  if (ac_has_tag (snpObj, "Potential_splice_disruption") )
    vtxtPrint (location, "_potential_splice_disruption") ;

  /* protect the future printf */
  if (! vtxtPtr (gTxt)) vtxtPrintf (gTxt, " ") ;
  
  if (1)
    { 
      const char *errors = 0 ;
      const char *ccp ;
      vTXT txt = vtxtHandleCreate (h) ;
      
      vtxtPrintf (txt, "Variant %s\n", ac_protect (ac_name(snpObj), h)) ;
      ccp = vtxtPtr (gTxt) ;
      if (ccp && strlen (ccp) > 2)
	vtxtPrintf (txt, "gName %s\n",  ac_protect (ccp, h)) ;
      ccp = vtxtPtr (cTxt) ;
      if (ccp && strlen (ccp) > 2)
	vtxtPrintf (txt, "rName %s\n",  ac_protect (ccp, h)) ;
       ccp = vtxtPtr (pTxt) ;
      if (ccp && strlen (ccp) > 2)
	 vtxtPrintf (txt, "pName %s\n",  ac_protect (ccp, h)) ;
      vtxtPrintf (txt, "\n") ;
      ac_parse (snp->db, vtxtPtr (txt), &errors, 0, h) ;  
      if (errors && *errors) 
	fprintf(stderr, "snpPrettyNames parsing error %s\n", errors) ;
    }

  ac_free (h) ;
  return type ;
} /* snpPrettyNames */

/*************************************************************************************/

static AC_OBJ snpGetObj (SNP *snp, const char *targetName, int pos, const char *typeName, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  char *cp = 0 ;
  int n ;
  AC_OBJ obj = 0 ;
  vTXT txt = vtxtHandleCreate (h) ;

  vtxtPrintf (txt, "%s:%d:%s"
	      , targetName
	      , pos
	      , typeName
	      ) ;
  cp = vtxtPtr (txt) ;
  n = strlen (cp) ; if (n > 2) cp += n - 2 ;
  if (*cp == '>') *cp = '2' ;
  else cp = 0 ;

  if (snp->db && txt &&  vtxtPtr (txt))
    obj = ac_get_obj (snp->db, "Variant", vtxtPtr (txt), h0) ;

  ac_free (h) ;
  return obj ;
} /* snpGetObj */

 /*************************************************************************************/

static void snpStLouisData (SNP *snp, AC_OBJ snpObj, ACEOUT ao, int maxNHLBI) 
{
  AC_HANDLE h = ac_new_handle () ;
  int ir ;
  AC_TABLE tt1 = 0 ;
  
  if (snpObj)
    {
      tt1 = ac_tag_table (snpObj, "Saint_Louis", h) ;
      if (tt1)
	{
	  for (ir = 0 ; ir <  maxNHLBI ; ir++)
	    aceOutf (ao, "\t%s", ac_table_printable (tt1, ir, 1, "-")) ;
	}
      else
	{
	  for(ir = 0 ; ir < maxNHLBI ; ir++)
	    aceOutf (ao, "\t-") ;
	}
   }
  else
    {
      for(ir = 0 ; ir <  maxNHLBI ; ir++)
	aceOutf (ao, "\t-") ;
    }

  ac_free (h) ;
  return ;
} /* snpStLouisData */

/*************************************************************************************/

static void snpFrequencyTableHeader (ACEOUT ao, int maxNHLBI,  AC_TABLE groups, const char *title)
{
  int ir ;

  aceOutf (ao, "\n##\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t") ;
  for (ir = 0 ; ir < maxNHLBI ; ir++) aceOut (ao, "\t") ; 
  for (ir = 0 ;  groups && ir < groups->rows; ir++)
    aceOutf (ao, "\t\t\t\t\t\t\t\t\t") ;
  aceOutf (ao, "\t%s", title) ;
}

/*************************************************************************************/

/* input snp-count table (per snp, counts of m/w on both strands 
 * db -> MetaDB : gives runs/groups and runs metadata
 * export a table of frequency, one column per run. one line per snp
 */ 
static int snpFrequencyTable (SNP *snp, BOOL preRun)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = snp->ai ;
  ACEOUT ao, aoF = 0, aoC = 0 ;
  DICT *eeDict = eeDictCreate (snp, h) ;  /* a,t,g,c and all error types in canonical order */
  DICT *selectDict = snp->selectDict ;
  DICT *targetDict = dictHandleCreate (1000, h) ;
  DICT *runDict = dictHandleCreate (1000, h) ;
  int nn = 0 ;
  int ir, run, gPos, line = 0 , pass ;
  int gMap ;
  int type, pos, target, oldType = 0, oldPos = 0, oldTarget = 0 ;
  AC_TABLE runs =  snpProjectRuns (snp, h) ;
  KEYSET r2ir = keySetHandleCreate (h) ;
  BOOL new = TRUE ;
  const char *ccp ;
  const char *errors = 0 ;
  int max = runs ? runs->rows + 1 : 1 ;
  int prettyType = 0 ;
  int mp = 0, mm = 0, wp = 0, wm = 0 ; /* wild types and mutants on both strands */ 
  int snpNA[max], snpMeasured[max], snpDropped[max], snpWild[max], snpLow[max], snpMid[max], snpHigh[max], snpPure[max] ;
  int        snpExAny[max], snpExWild[max], snpExLow[max], snpExMid[max], snpExHigh[max], snpExPure[max] ; /* Exonic */
  int        snpPcAny[max], snpPcWild[max], snpPcLow[max], snpPcMid[max], snpPcHigh[max], snpPcPure[max] ; /* protein changing */
  int snpPcAnyAny = 0, snpPcWildAny = 0, snpPcLowAny = 0, snpPcMidAny = 0, snpPcHighAny = 0, snpPcPureAny = 0 ;
  int snpExAnyAny = 0,snpExWildAny = 0, snpExLowAny = 0, snpExMidAny = 0, snpExHighAny = 0, snpExPureAny = 0 ; /* Exonic */

  Array ff = 0 ; 
  Array ff1 = 0 ; 
  Array ff2 = 0 ; 
  float z, cover, mutant ;
  float minCover = snp->minCover ;
  Array hh, histos = arrayHandleCreate (100, Array, h) ;
  vTXT txt = vtxtHandleCreate (h) ;
  vTXT gTxt = vtxtHandleCreate (h) ;
  vTXT cTxt = vtxtHandleCreate (h) ;
  vTXT pTxt = vtxtHandleCreate (h) ;
  vTXT geneboxes = vtxtHandleCreate (h) ;
  vTXT avgeneboxes = vtxtHandleCreate (h) ;
  vTXT location = vtxtHandleCreate (h) ;
  vTXT gsnippet = vtxtHandleCreate (h) ;
  vTXT snippet = vtxtHandleCreate (h) ;
  vTXT pSnippet = vtxtHandleCreate (h) ; 
  vTXT pType = vtxtHandleCreate (h) ;
  char gType[32], cType[32] ;

  int  maxNHLBI = 0 ;
  AC_OBJ snpObj = 0 ;
  AC_TABLE groups = 0, groups2runs = 0 ;
  vTXT stlTxt = 0 ;

  memset (snpNA, 0, sizeof (snpNA)) ;
  memset (snpMeasured, 0, sizeof (snpMeasured)) ;
  memset (snpDropped, 0, sizeof (snpDropped)) ;
  memset (snpWild, 0, sizeof (snpWild)) ;
  memset (snpLow, 0, sizeof (snpLow)) ;
  memset (snpMid, 0, sizeof (snpMid)) ;
  memset (snpHigh, 0, sizeof (snpHigh)) ;
  memset (snpPure, 0, sizeof (snpPure)) ;

  memset (snpExAny, 0, sizeof (snpExAny)) ;
  memset (snpExWild, 0, sizeof (snpExWild)) ;
  memset (snpExLow, 0, sizeof (snpExLow)) ;
  memset (snpExMid, 0, sizeof (snpExMid)) ;
  memset (snpExHigh, 0, sizeof (snpExHigh)) ;
  memset (snpExPure, 0, sizeof (snpExPure)) ;
  
  memset (snpPcAny, 0, sizeof (snpPcAny)) ;
  memset (snpPcWild, 0, sizeof (snpPcWild)) ;
  memset (snpPcLow, 0, sizeof (snpPcLow)) ;
  memset (snpPcMid, 0, sizeof (snpPcMid)) ;
  memset (snpPcHigh, 0, sizeof (snpPcHigh)) ;
  memset (snpPcPure, 0, sizeof (snpPcPure)) ;
  
  if (1)  /* dedicated code to add info imported from the StLouis server of the NHLBI project, see dir RESULTS/StLouis */
    {
      AC_HANDLE h1 = ac_new_handle () ;
      AC_ITER iter = ac_dbquery_iter (snp->db, "Find Variant Saint_louis", h1) ;
      AC_OBJ obj1 = 0 ;
      AC_TABLE tt1 = 0 ;

      if (iter && (obj1 = ac_iter_obj (iter)) && (tt1 = ac_tag_table (obj1, "Saint_Louis", h1)))
	{
	  stlTxt = vtxtHandleCreate (h) ;
	  maxNHLBI = tt1->rows ;
	  for (ir = 0 ; ir < tt1->rows ; ir++)
	    vtxtPrintf (stlTxt, "\t%s", ac_table_printable (tt1, ir, 0, "-")) ;
	}
      ac_free (obj1) ;
      ac_free (h1) ;
  }
  /*
    AC_TABLE groups = ac_bql_table (snp->db, hprintf (h,"select st,   g from g in class \"run\" where  g#runs and g->project == \"%s\" and not  g#add_counts, st in g->sorting_title", snp->project), 0, 0, 0, h) ;
    AC_TABLE groups2runs = ac_bql_table (snp->db, hprintf (h,"select g,r from g in class \"run\" where g#runs and g->project == \"%s\" and not g#add_counts, r in g->runs", snp->project), 0, 0, 0, h) ;
  */

  if (snp->dbCount)
    snpDbCount (snp) ;
  if (! snp->frequencyHisto && ! snp->frequencyTable)
    goto done ;
  vtxtClear (txt) ;
  vtxtPrintf (txt, "colonne 1 \n class Run \n Condition Runs && SNP && ! Add_counts && project == \"%s\"\n", snp->project) ;
  vtxtPrintf (txt, "colonne 2 \nFrom 1\nClass Text\nTag Sorting_title\n\n") ;
  groups = ac_tablemaker_table (snp->db, vtxtPtr (txt), 0, ac_tablemaker_text , 0 , "2+1", &errors, h) ;

  vtxtClear (txt) ;
  vtxtPrintf (txt, "colonne 1 \n class Run \n Condition Runs && ! Add_counts && project == \"%s\"\n", snp->project) ;
  vtxtPrintf (txt, "colonne 2 \nFrom 1\nClass Run\nTag Runs\n\n") ;
  groups2runs = ac_tablemaker_table (snp->db, vtxtPtr (txt), 0, ac_tablemaker_text , 0 , 0, &errors, h) ;


  if (! preRun && snp->frequencyTable)
    for (pass = 0 ; pass < 2 ; pass++)
      {
	switch (pass)
	  {
	  case 0:
	    ao = aoF = aceOutCreate (snp->outFileName, ".snp_list_with_allele_frequency_per_sample.txt", snp->gzo, h) ;
	    break ;
	  case 1:
	    ao = aoC = aceOutCreate (snp->outFileName, ".snp_list_with_allele_counts_per_sample.txt", snp->gzo, h) ;
	    break ;
	  }
	    
	aceOutDate (ao, "###", snp->project) ;
	if (maxNHLBI)
	  {
	    aceOutf (ao, "## The NHLBI columns were downloaded for comparison from\n") ;
	    aceOutf (ao, "## Citation: Exome Variant Server, NHLBI GO Exome Sequencing Project (ESP), Seattle, WA (URL: http://evs.gs.washington.edu/EVS/) downloaded  December 10, 2014\n") ;
	    aceOutf (ao, "## Acknoledgment for publication:: The authors would like to thank the NHLBI GO Exome Sequencing Project and its ongoing studies which produced and provided exome variant calls for comparison: the Lung GO Sequencing Project (HL-102923), the WHI Sequencing Project (HL-102924), the Broad GO Sequencing Project (HL-102925), the Seattle GO Sequencing Project (HL-102926) and the Heart GO Sequencing Project (HL-103010).\n") ;
	  }
	aceOutf (ao, "# Chromosome\tPosition\tStrand\tMagic id\tVCF id\tGenome SNV type on plus strand of genome\tSNiPpets on plus strand of the genome\tSNiPpets on transcribed strand if applicable\tGenome SNV\tGene area\tAceView gene area\tcDNA SNV\tcDNA SNV type on plus strand of transcript\tProtein SNV\tProtein SNV type\tProtein snippets\tLocation relative to transcripts\tFull SNV name") ;
	if (maxNHLBI) aceOut (ao, vtxtPtr (stlTxt)) ;
	
	if (0)
	  { /* eliminate these 2 columns, add magic name and vcf name as columns 1 and 2 */
	    aceOutf (ao, "\tRS\tDifferential between") ;
	  }
	
	aceOutf (ao, "\tVariant counts on the plus strand in the whole cohort\tVariant counts on the minus strand in the whole cohort\tReference counts on the plus strand in the whole cohort\tReference counts on the minus strand in the whole cohort\tCompatibility\tchi2") ;
	aceOutf (ao, "\tMeasured in\tNA\tReference (<5% variant)\tLow (5-20%)\tMid (20-80%)\tHigh(80-95%)\tPure variant (>95%)\tAllele frequency in cohort\tAverage allele frequency") ;
	for (ir = 0 ; groups && ir < groups->rows; ir++)
	  {
	    ccp = ac_table_printable ( groups, ir, 0, "") ;
	    aceOutf (ao, "\t\t%s Measured in\t%s NA\t%s Reference (<5% variant)\t%s Low (5-20%)\t%s Mid (20-80%)\t%s High (80-95%)\t%s Pure variant (>95%)\tAllele frequency in %s"
		     , ccp
		     , ccp
		     , ccp
		     , ccp
		     , ccp
		     , ccp
		     , ccp
		     , ccp
		     ) ;
	  }
	aceOutf (ao, "\tRun\tMaximal allele frequency in any run") ;
	for (ir = 0 ; ir < runs->rows ; ir++)
	  aceOutf (ao, "\t%s", ac_table_printable (runs, ir, 0, "-")) ;
	
	snpFrequencyTableHeader (ao, maxNHLBI, groups, "Sorting title\tAny run") ;
	for (ir = 0 ; ir < runs->rows ; ir++)
	  aceOutf (ao, "\t%s", ac_table_printable (runs, ir, 1, "-")) ;
	
	snpFrequencyTableHeader (ao, maxNHLBI, groups, "Title\tAny run") ;
	for (ir = 0 ; ir < runs->rows ; ir++)
	  aceOutf (ao, "\t%s", ac_table_printable (runs, ir, 2, "-")) ;
	
	snpFrequencyTableHeader (ao, maxNHLBI, groups, "Sample\tAny sample") ;
	for (ir = 0 ; ir < runs->rows ; ir++)
	  aceOutf (ao, "\t%s", ac_table_printable (runs, ir, 3, "-")) ;
	aceOutf (ao, "\n") ;
      }

  dictAdd (runDict, "_zero_", 0) ;
  for (ir = 0 ; ir < runs->rows ; ir++)
    {
      ccp = ac_table_printable (runs, ir, 0, "_zero_") ;
      dictAdd (runDict, ccp, &run) ;
      keySet (r2ir, run) = ir ;
    }
  
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ; if (! ccp || *ccp == '#') continue ;
      vtxtClear (txt) ;
      vtxtClear (gTxt) ;
      vtxtClear (cTxt) ;
      vtxtClear (pTxt) ;
      vtxtClear (geneboxes) ;
      vtxtClear (avgeneboxes) ;
      vtxtClear (gsnippet) ;
      vtxtClear (snippet) ;
      vtxtClear (pSnippet) ;
      vtxtClear (pType) ;
      if (0 && line++ > 20) break ;
      gType[0] = 0 ;
      cType[0] = 0 ;
      
      if (selectDict)
	{
	  int zone, iSelect ;
	  
	  if (! dictFind (selectDict, ccp, &zone)) 
	    continue ;
	  dictAdd (targetDict, ccp, &target) ;
	  aceInStep (ai, '\t') ; if (! aceInInt (ai, &pos)) continue ;
	  
	  if (snp->selectZone)
	    {
	      ZONE *zp ;
	      int iSelectMax = arrayMax (snp->selectZone) ;
	      
	      iSelect = keySet (snp->selectZoneIndex, zone) ;
	      zp = arrp (snp->selectZone, iSelect, ZONE) ;
	      while (iSelect > 0 && (zp->zone > zone || zp->a2 > pos)) { zp-- ; iSelect-- ; }
	      while (iSelect < iSelectMax  && (zp->zone < zone || (zp->a1 && zp->a2 < pos))) { zp++ ; iSelect++ ; }
	      if (zp->zone == zone && (zp->a1 + zp->a2 == 0 || (zp->a1 <= pos && zp->a2 >= pos))) ;
	      else continue ;
	    }
	}
      else
	{
	  dictAdd (targetDict, ccp, &target) ;
	  aceInStep (ai, '\t') ; if (! aceInInt (ai, &pos)) continue ;
	}
      
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; if (! ccp || *ccp == '#') continue ;
      dictAdd (eeDict, ccp, &type) ;
      
      if (target != oldTarget || pos != oldPos || type != oldType)
	{
	  if (snp->frequencyTable)
	    {
	      if (oldTarget) 
		{
		  BOOL isDrop = snpFrequencyTableExportOneLine (snp, aoF, aoC, max - 1, ff, ff1, ff2, oldTarget, oldPos, oldType, groups, groups2runs, r2ir, runDict, mp, mm, wp, wm) ;
		  if (! isDrop)
		    {  
		      float zPc = -10, zEx = -10 ;
		      for (ir = 0 ; ir < runs->rows ; ir++)
			{
			  float z = array (ff, ir, float) - 1000 ;
			  
			  if (z >= 0) snpMeasured[ir]++ ;
			  
			  if (z < 0) snpNA[ir]++ ;
			  else if (z < 5) snpWild[ir]++ ;
			  else if (z < 20) snpLow[ir]++ ;
			  else if (z < 80) snpMid[ir]++ ;
			  else if (z < 95) snpHigh[ir]++ ;
			  else snpPure[ir]++ ;
			  
			  if (prettyType & 0x2) /* ProteinChanging */
			    {
			      if (zPc < z) zPc = z ;
			      if (z >= 0) snpPcAny[ir]++ ;
			      if (z < 0) ;
			      else if (z < 5) snpPcWild[ir]++ ;
			      else if (z < 20) snpPcLow[ir]++ ;
			      else if (z < 80) snpPcMid[ir]++ ;
			      else if (z < 95) snpPcHigh[ir]++ ;
			      else snpPcPure[ir]++ ;
			    }
			  if (prettyType & 0x1) /* Exonic */
			    {
			      if (zEx < z) zEx = z ;
			      if (z >= 0)  snpExAny[ir]++ ;
			      if (z < 0) ;
			      else if (z < 5) snpExWild[ir]++ ;
			      else if (z < 20) snpExLow[ir]++ ;
			      else if (z < 80) snpExMid[ir]++ ;
			      else if (z < 95) snpExHigh[ir]++ ;
			      else snpExPure[ir]++ ;
			    }
			  
			}
		      if (prettyType & 0x2) /* ProteinChanging */
			{
			  z = zPc ; 
			  snpPcAnyAny++ ;
			  if (z < 0) ;
			  else if (z < 5) snpPcWildAny++ ;
			  else if (z < 20) snpPcLowAny++ ;
			  else if (z < 80) snpPcMidAny++ ;
			  else if (z < 95) snpPcHighAny++ ;
			  else snpPcPureAny++ ;
			}
		      
		      if (prettyType & 0x1) /* Exonic */
			{
			  z = zEx ;
			  snpExAnyAny++ ;
			  if (z < 0) ;
			  else if (z < 5) snpExWildAny++ ;
			  else if (z < 20) snpExLowAny++ ;
			  else if (z < 80) snpExMidAny++ ;
			  else if (z < 95) snpExHighAny++ ;
			  else snpExPureAny++ ;
			}
		    }
		  else 
		    {
		      for (ir = 0 ; ir < runs->rows ; ir++)
			{
			  float z = array (ff, ir, float) - 1000 ;
			  if (z >= 0) 
			    snpDropped[ir]++ ;
			  else
			    snpNA[ir]++ ;
			}
		    }
		}

	      prettyType = 0 ;
	      ccp =  dictName (targetDict, target) ;
	      strncpy (gType, dictName (eeDict, type), 31) ;
	      ac_free (snpObj) ; 
	      snpObj = snpGetObj (snp , ccp, pos, gType, h) ;
	      oldTarget = 0 ;
	      if (snp->differential && ! ac_has_tag (snpObj, "Differential"))
		continue ;
	      if (aoF)
		for (pass = 0 ; pass < 2 ; pass++)
		  {
		    ao = pass ? aoF : aoC ;
		    prettyType = snpPrettyNames (snp, snpObj, gTxt, cTxt, pTxt, &gMap, &gPos, gType, cType, pType, location, geneboxes, avgeneboxes, gsnippet, snippet, pSnippet) ;
		    
		    aceOutf (ao, "\n") ;
		    
		    aceOutf (ao, "%s\t%d\t+\t%s",  ac_tag_printable (snpObj, "IntMap", " "), gPos, gType) ;
		    aceOutf (ao, "\t%s", ac_name(snpObj)) ;
		    if (1)
		      {
			AC_TABLE vcf = ac_tag_table (snpObj, "VCF", h) ;
			
			if (vcf && vcf->cols >= 3 &&
			    (ccp = ac_tag_printable (snpObj, "IntMap", 0))
			    ) 
			  aceOutf (ao, "\t%s:%s:%s:%s"
				   , ccp
				   , ac_table_printable (vcf, 0, 0, "")
				   , ac_table_printable (vcf, 0, 1, "")
				   , ac_table_printable (vcf, 0, 2, "")
				   ) ;
			else
			  aceOutf (ao, "\t") ;
			ac_free (vcf) ;
		      }
		    
		    aceOutf (ao, "\t%s", vtxtPtr (gsnippet) ? vtxtPtr (gsnippet) : "") ; /* genomic snippet */
		    aceOutf (ao, "\t%s", vtxtPtr (snippet) ? vtxtPtr (snippet) : "") ; /* snippet */
		    aceOutf (ao, "\t%s", vtxtPtr (gTxt)) ; /* genomic g.name */
		    aceOutf (ao, "\t%s", vtxtPtr (geneboxes) ? vtxtPtr (geneboxes) : "") ;
		    aceOutf (ao, "\t%s", vtxtPtr (avgeneboxes) ? vtxtPtr (avgeneboxes) : "") ;
		    aceOutf (ao, "\t%s", vtxtPtr (cTxt) ? vtxtPtr (cTxt) : "") ; /* cDNA c.name */
		    aceOut (ao, "\t") ;
		    
		    aceOutf (ao, "%s", cType) ;          /* cDNA type, possibly the complement of the type */
		    aceOutf (ao, "\t%s", vtxtPtr (pTxt) ? vtxtPtr (pTxt) : "") ; /* protein p.Name */
		    aceOutf (ao, "\t%s", vtxtPtr (pType) ? vtxtPtr (pType) : "") ;
		    aceOutf (ao, "\t%s", vtxtPtr (pSnippet) ?  vtxtPtr (pSnippet) : "") ;
		    /* composite name */
		    aceOutf (ao, "\t%s", vtxtPtr (location)) ;
		    aceOutf (ao, "\t%s", vtxtPtr (gTxt)) ; /* genomic g.name */ 
		    if (vtxtPtr (cTxt))
		      aceOutf (ao, ", %s", vtxtPtr (cTxt)) ;
		    if (vtxtPtr (pTxt))
		      aceOutf (ao, ", %s", vtxtPtr (pTxt)) ;
		    if (vtxtPtr (pSnippet))
		      aceOutf (ao, " ; %s", vtxtPtr (pSnippet)) ;
		  }
	    }
	  
	  ff = arrayReCreate (ff, max, float) ; wp = wm = mp = mm = 0 ;
	  ff1 = arrayReCreate (ff1, max, int) ; 
	  ff2 = arrayReCreate (ff2, max, int) ; 
	  oldTarget = target ; oldPos = pos ; oldType = type ; new = TRUE ;
	  nn++ ;
	}
      
      if (aoF && snp->frequencyTable && new) 
	{  
	  if (maxNHLBI)
	    {
	      snpStLouisData (snp, snpObj, aoF, maxNHLBI) ;
	      snpStLouisData (snp, snpObj, aoC, maxNHLBI) ;
	    }
	  if (0)
	    {
	      aceOutf (aoF, "\t-") ; /* RS number */
	      aceOutf (aoC, "\t-") ; /* RS number */
	    }
	}
      if (0 && aoF && snp->frequencyTable && new) 
	{ 
	  AC_TABLE dd = ac_tag_table ( snpObj, "Differential", h)  ;
	  ao = aoF ;
	  aceOutf (ao, "\t") ;
	  if (dd && dd->rows && dd->cols >= 2)
	    for (ir = 0 ; ir < dd->rows ; ir++)
	      {
		AC_HANDLE h1 = ac_new_handle () ;
		AC_OBJ run1, run2 ;

		aceOutf (ao, "%s%s_"
			 , ir ? ", " : ""
			 , ac_table_printable (dd, ir, 0, "") /* differential type */
			 ) ;

		run1 = ac_table_obj (dd, ir, 1, h) ;
		run2 = ac_table_obj (dd, ir, 2, h) ;

		aceOutf (ao, "%s", ac_table_printable (dd, ir, 1, "")) ;
		if (ac_has_tag (run1, "Sorting_title"))
		  aceOutf (ao, "(%s)", ac_tag_printable (run1, "Sorting_title", "-")) ;
		aceOutf (ao, "/%s", ac_table_printable (dd, ir, 2, "")) ;
		if (ac_has_tag (run2, "Sorting_title"))
		  aceOutf (ao, "(%s)", ac_tag_printable (run2, "Sorting_title", "-")) ;

		ac_free (h1) ;
	      }
	  else
	    aceOutf (ao, "-") ;
	}

      new = FALSE ;
      
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;  /* jump the 2 snippets */
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;  /* jump the 2 genomic snippets */
      aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; 
      
      if (! ccp || *ccp == '#') continue ;
      if (! dictFind (runDict, ccp, &run)) continue ;
      aceInStep (ai, '\t') ;  if (! aceInFloat (ai, &z)) continue ;
      aceInStep (ai, '\t') ;  if (! aceInFloat (ai, &z)) continue ;
      aceInStep (ai, '\t') ;  if (! aceInFloat (ai, &cover)) continue ;

      ir =  keySet(r2ir, run) ;
      if (cover < minCover) 
	continue ; 

      aceInStep (ai, '\t') ;  if (! aceInFloat (ai, &mutant)) continue ;
      if (mutant > cover) mutant = cover ;
      if (snp->frequencyHisto)
	{
	  hh = array (histos, ir, Array) ;
	  if (! hh)
	    hh = array (histos, ir, Array) = arrayHandleCreate (101, double, h) ;
	  utSmoothHisto (hh, mutant, cover, 1) ;
	  { int zi = .4 + 100.0 * mutant / (cover * 1.0) ;
	    if (cover > 0) array (hh, 120 + zi, double) += 1 ;
	  }
	}
      if (1)
	{
	  float z =  100.0 * mutant/((float)cover) ;
	  array (ff, ir, float) = 1000 + z ;
	  array (ff1, ir, int) = mutant ;
	  array (ff2, ir, int) = cover ;
	}
      
      {
	int x ;
	aceInStep (ai, '\t') ;  if (! aceInInt (ai, &x)) continue ; /* wild type */
	aceInStep (ai, '\t') ;  if (! aceInInt (ai, &x)) continue ; mp += x ;
	aceInStep (ai, '\t') ;  if (! aceInInt (ai, &x)) continue ; wp += x ;
	aceInStep (ai, '\t') ;  if (! aceInInt (ai, &x)) continue ; mm += x ;
	aceInStep (ai, '\t') ;  if (! aceInInt (ai, &x)) continue ; wm += x ;
      }
    }
  if (snp->frequencyTable)
    {
      float z ;
      int maxMeasured = 0 ;

      for (ir = 0 ; ir < runs->rows ; ir++)
	if (maxMeasured < snpMeasured[ir]) 
	  maxMeasured = snpMeasured[ir] ;
      if ( oldTarget) 
	{
	  BOOL isDrop = snpFrequencyTableExportOneLine (snp, aoF, aoC, max - 1, ff, ff1, ff2, oldTarget, oldPos, oldType, groups, groups2runs, r2ir, runDict, mp, mm, wp, wm) ;  
	  if (! isDrop)
	    {  
		      float zPc = -10, zEx = -10 ;
		      for (ir = 0 ; ir < runs->rows ; ir++)
			{
			  float z = array (ff, ir, float) - 1000 ;
			  
			  if (z >= 0) snpMeasured[ir]++ ;

			  if (z < 0) snpNA[ir]++ ;
			  else if (z < 5) snpWild[ir]++ ;
			  else if (z < 20) snpLow[ir]++ ;
			  else if (z < 80) snpMid[ir]++ ;
			  else if (z < 95) snpHigh[ir]++ ;
			  else snpPure[ir]++ ;
			  
			  if (prettyType & 0x2) /* ProteinChanging */
			    {
			      if (zPc < z) zPc = z ;
			      if (z >= 0) snpPcAny[ir]++ ;
			      if (z < 0) ;
			      else if (z < 5) snpPcWild[ir]++ ;
			      else if (z < 20) snpPcLow[ir]++ ;
			      else if (z < 80) snpPcMid[ir]++ ;
			      else if (z < 95) snpPcHigh[ir]++ ;
			      else snpPcPure[ir]++ ;
			    }
			  if (prettyType & 0x1) /* Exonic */
			    {
			      if (zEx < z) zEx = z ;
			      if (z >= 0)  snpExAny[ir]++ ;
			      if (z < 0) ;
			      else if (z < 5) snpExWild[ir]++ ;
			      else if (z < 20) snpExLow[ir]++ ;
			      else if (z < 80) snpExMid[ir]++ ;
			      else if (z < 95) snpExHigh[ir]++ ;
			      else snpExPure[ir]++ ;
			    }
			  
			}
		      if (prettyType & 0x2) /* ProteinChanging */
			{
			  z = zPc ; 
			  snpPcAnyAny++ ;
			  if (z < 0) ;
			  else if (z < 5) snpPcWildAny++ ;
			  else if (z < 20) snpPcLowAny++ ;
			  else if (z < 80) snpPcMidAny++ ;
			  else if (z < 95) snpPcHighAny++ ;
			  else snpPcPureAny++ ;
			}

		      if (prettyType & 0x1) /* Exonic */
			{
			  z = zEx ;
			  snpExAnyAny++ ;
			  if (z < 0) ;
			  else if (z < 5) snpExWildAny++ ;
			  else if (z < 20) snpExLowAny++ ;
			  else if (z < 80) snpExMidAny++ ;
			  else if (z < 95) snpExHighAny++ ;
			  else snpExPureAny++ ;
			}
		    }
		  else 
		    {
		      for (ir = 0 ; ir < runs->rows ; ir++)
			{
			  float z = array (ff, ir, float) - 1000 ;
			  if (z >= 0) 
			    snpDropped[ir]++ ;
			  else
			    snpNA[ir]++ ;
			}
		    }
	}
      if (aoF)
	for (pass = 0 ; pass < 2 ; pass++)
	  {
	    float x ;
	    ao = pass ? aoF : aoC ;
	    
	    aceOutf (ao, "\n") ;
	    
	    /* export the global counts */
	    snpFrequencyTableHeader (ao, maxNHLBI, groups, hprintf (h, "All sites (variant in at least one sample) covered at least %d times ", snp->minCover )) ;
	    aceOutf (ao, "\t%d", snp->anyRunMeasured + snp->anyRunDropped) ;
	    for (ir = 0 ; ir < runs->rows ; ir++)
	      aceOutf (ao, "\t%d", snpMeasured[ir] + snpDropped[ir]) ;
	    snpFrequencyTableHeader (ao, maxNHLBI, groups, hprintf (h, "Sites not measurable (less than %d reads)", snp->minCover)) ;
	    aceOutf (ao, "\t%d", snp->anyRunNA) ;
	    for (ir = 0 ; ir < runs->rows ; ir++)
	      aceOutf (ao, "\t%d", snpNA[ir]) ;
	    snpFrequencyTableHeader (ao, maxNHLBI, groups, "Rejected monomodal sites (likely sequencing or mapping error)") ;
	    aceOutf (ao, "\t%d", snp->anyRunDropped) ;
	    for (ir = 0 ; ir < runs->rows ; ir++)
	      aceOutf (ao, "\t%d", snpDropped[ir]) ;
	    snpFrequencyTableHeader (ao, maxNHLBI, groups, "Well measured sites") ;
	    aceOutf (ao, "\t%d",  snp->anyRunMeasured) ;
	    for (ir = 0 ; ir < runs->rows ; ir++)
	      aceOutf (ao, "\t%d", snpMeasured[ir]) ;
	    
	    x = snp->anyRunMeasured > 0 ? 100.0 / snp->anyRunMeasured : 0 ;
	    snpFrequencyTableHeader (ao, maxNHLBI, groups, "% Reference (<5% variant)") ;
	    aceOutf (ao, "\t%.2f", x * snp->anyRunWild) ;
	    for (ir = 0 ; ir < runs->rows ; ir++)
	      { z = snpMeasured[ir] ? snpMeasured[ir] : 1 ; aceOutf (ao, "\t%.2f", 100.0 * snpWild[ir]/z) ; }
	    snpFrequencyTableHeader (ao, maxNHLBI, groups, "% Low (5-20% variant)") ;
	    aceOutf (ao, "\t%.2f", x * snp->anyRunLow) ;
	    for (ir = 0 ; ir < runs->rows ; ir++)
	      { z = snpMeasured[ir] ? snpMeasured[ir] : 1 ;aceOutf (ao, "\t%.2f", 100.0 * snpLow[ir]/z) ; }
	    snpFrequencyTableHeader (ao, maxNHLBI, groups, "% Mid (20-80% variant)") ;
	    aceOutf (ao, "\t%.2f", x * snp->anyRunMid) ;
	    for (ir = 0 ; ir < runs->rows ; ir++)
	      { z = snpMeasured[ir] ? snpMeasured[ir] : 1 ;aceOutf (ao, "\t%.2f", 100.0 * snpMid[ir]/z) ; }
	    snpFrequencyTableHeader (ao, maxNHLBI, groups, "% High (80-95% variant)") ;
	    aceOutf (ao, "\t%.2f", x * snp->anyRunHigh) ;
	    for (ir = 0 ; ir < runs->rows ; ir++)
	      { z = snpMeasured[ir] ? snpMeasured[ir] : 1 ;aceOutf (ao, "\t%.2f", 100.0 * snpHigh[ir]/z) ; }
	    snpFrequencyTableHeader (ao, maxNHLBI, groups, "% Pure variant (>95% variant)") ;
	    aceOutf (ao, "\t%.2f", x * snp->anyRunPure) ;
	    for (ir = 0 ; ir < runs->rows ; ir++)
	      { z = snpMeasured[ir] ? snpMeasured[ir] : 1 ;aceOutf (ao, "\t%.2f", 100.0 * snpPure[ir]/z) ; }
	    if (1)
	      {
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "% Variant (20-100%)" ) ;
		aceOutf (ao, "\t%.2f", x * snp->anyRunMid + x * snp->anyRunHigh + x * snp->anyRunPure) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  { z = snpMeasured[ir] ? snpMeasured[ir] : 1 ; aceOutf (ao, "\t%.2f", 100.0 * (snpMid[ir] + snpHigh[ir] + snpPure[ir] )/z) ; }
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "Variant (20-100%)") ;
		aceOutf (ao, "\t%d", snp->anyRunPure + snp->anyRunMid + snp->anyRunHigh) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpPure[ir] +  snpMid[ir] + snpHigh[ir]) ;
	      }
	    if (1)
	      {
		
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "Reference (<5% variant)") ;
		aceOutf (ao, "\t%d", snp->anyRunWild) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpWild[ir]) ;
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "Low (5-20%)") ;
		aceOutf (ao, "\t%d", snp->anyRunLow) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpLow[ir]) ;
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "Mid (20-80%)") ;
		aceOutf (ao, "\t%d", snp->anyRunMid) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpMid[ir]) ;
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "High (80-95%)") ;
		aceOutf (ao, "\t%d", snp->anyRunHigh) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpHigh[ir]) ;
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "Pure variant (>95%)") ;
		aceOutf (ao, "\t%d", snp->anyRunPure) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpPure[ir]) ;
	      }
	    
	    if (0)
	      { /* very wrong in posterior code */
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "Heterozygosity index") ;
		aceOutf (ao, "\t%.2f", snp->anyRunMeasured ? (100.0 * (snp->anyRunLow + snp->anyRunMid + snp->anyRunHigh)/(1.0*snp->anyRunMeasured)) : 0.0) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  {	      
		    if (snpMeasured[ir] > 0) 
		      {
			float z = snpLow[ir] + snpMid[ir] + snpHigh[ir] ;
			z = 100.0 * z / snpMeasured[ir] ;
			aceOutf (ao, "\t%.2f", z) ;
		      }
		    else
		      aceOutf (ao, "\t%.2f", 0) ;
		  }
	      }
	    
	    
	    /********************/ 
	    
	    snpFrequencyTableHeader (ao, maxNHLBI, groups, "Protein changing Sites") ;
	    aceOutf (ao, "\t%d", snpPcAnyAny) ;
	    for (ir = 0 ; ir < runs->rows ; ir++)
	      aceOutf (ao, "\t%d", snpPcAny[ir]) ;
	    snpFrequencyTableHeader (ao, maxNHLBI, groups, "Protein changing Reference (<5% variant)") ;
	    aceOutf (ao, "\t%d", snpPcWildAny) ;
	    for (ir = 0 ; ir < runs->rows ; ir++)
	      aceOutf (ao, "\t%d", snpPcWild[ir]) ;
	    
	    if (1)
	      {
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "Protein changing Low (5-20% AF)") ;
		aceOutf (ao, "\t%d", snpPcLowAny) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpPcLow[ir]) ;
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "Protein changing Mid (20-80% AF)") ;
		aceOutf (ao, "\t%d", snpPcMidAny) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpPcMid[ir]) ;
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "Protein changing High (80-95% AF)") ;
		aceOutf (ao, "\t%d", snpPcHighAny) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpPcHigh[ir]) ;
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "Protein changing Pure variant (>95% AF)") ;
		aceOutf (ao, "\t%d", snpPcPureAny) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpPcPure[ir]) ;
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "Protein changing variant (20-95% AF)") ;
		aceOutf (ao, "\t%d", snpPcPureAny) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpPcMid[ir] + snpPcHigh[ir] + snpPcPure[ir]) ;
	      }
	    if (0)
	      {
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "Protein changing Intermediate frequency (5-95%)") ;
		aceOutf (ao, "\t%d", snpPcLowAny + snpPcMidAny + snpPcHighAny) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpPcLow[ir]+snpPcMid[ir]+snpPcHigh[ir]) ;
	      }
	    
	    /********************/ 
	    
	    snpFrequencyTableHeader (ao, maxNHLBI, groups, "exonic Sites") ;  /* exonic, not Exonic, to influence later call to sort */
	    aceOutf (ao, "\t%d", snpExAnyAny) ;
	    for (ir = 0 ; ir < runs->rows ; ir++)
	      aceOutf (ao, "\t%d", snpExAny[ir]) ;
	    snpFrequencyTableHeader (ao, maxNHLBI, groups, "exonic Reference (<5% variant)") ;
	    aceOutf (ao, "\t%d", snpExWildAny) ;
	    for (ir = 0 ; ir < runs->rows ; ir++)
	      aceOutf (ao, "\t%d", snpExWild[ir]) ;
	    if (1)
	      {
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "exonic Low (5-20%)") ;
		aceOutf (ao, "\t%d", snpExLowAny) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpExLow[ir]) ;
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "exonic Mid (20-80%)") ;
		aceOutf (ao, "\t%d", snpExMidAny) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpExMid[ir]) ;
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "exonic High (80-95%)") ;
		aceOutf (ao, "\t%d", snpExHighAny) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpExHigh[ir]) ;
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "exonic Pure variant (>95%)") ;
		aceOutf (ao, "\t%d", snpExPureAny) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpExPure[ir]) ;
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "exonic variant (20-95%)") ;
		aceOutf (ao, "\t%d", snpExPureAny) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpExMid[ir] + snpExHigh[ir] + snpExPure[ir]) ;
	      }
	    if (0)
	      {
		snpFrequencyTableHeader (ao, maxNHLBI, groups, "exonic Intermediate frequency (5-95%)") ;
		aceOutf (ao, "\t%d", snpExLowAny + snpExMidAny + snpExHighAny) ;
		for (ir = 0 ; ir < runs->rows ; ir++)
		  aceOutf (ao, "\t%d", snpExLow[ir] + snpExMid[ir] + snpExHigh[ir]) ;
	      }
	    
	    /********************/ 
	    aceOutf (ao, "\n") ;
	  }
    }
  
  if (snp->frequencyHisto)
    {
      int zi, dzi = 0 ;
      ao = aceOutCreate (snp->outFileName, ".SNP_frequency_histos.txt", FALSE, h) ;
      aceOutDate (ao, "###", snp->project) ;
      aceOutf (ao, "# Run") ;
      for (ir = 0 ; ir < runs->rows ; ir++)
	aceOutf (ao, "\t%s", ac_table_printable (runs, ir, 1, "-")) ;
      aceOutf (ao, "\n# Sorting title") ;
      for (ir = 0 ; ir < runs->rows ; ir++)
	aceOutf (ao, "\t%s", ac_table_printable (runs, ir, 0, "-")) ;
       aceOutf (ao, "\n# Title") ;
      for (ir = 0 ; ir < runs->rows ; ir++)
	aceOutf (ao, "\t%s", ac_table_printable (runs, ir, 2, "-")) ;
       aceOutf (ao, "\n# Sample") ;
      for (ir = 0 ; ir < runs->rows ; ir++)
	aceOutf (ao, "\t%s", ac_table_printable (runs, ir, 3, "-")) ;
      for (zi = 0 ; zi <= 220 ; zi++)
	{
	  if (zi > 100 && zi < 119)  
	    aceOutf (ao, "\n") ;
	  else if (zi == 119) 
	    { aceOutf (ao, "\n# Raw histogram, no smoothing") ; dzi = 120 ; }
	  else
	    {
	      aceOutf (ao, "\n%d", zi - dzi) ;
	      for (ir = 0 ; ir < runs->rows ; ir++)
		{ 
		  hh = array (histos, ir, Array) ;
		  if (hh)
		    aceOutf (ao, "\t%.1f", array (hh, zi, double)) ;
		  else
		    aceOutf (ao, "\t-") ;
		}
	    }
	}
      aceOutf (ao, "\n") ;
    }

 done:
  ac_free (h) ;
  return nn ;
} /*  snpFrequencyTable */

/*************************************************************************************/
/*************************************************************************************/

static BOOL snpVcfCounts (AC_OBJ s, const char *tag, const char *run, int *m, int *w,  int *cover)
{
  AC_HANDLE h = ac_new_handle () ;BOOL ok = FALSE ;
  AC_TABLE tt = ac_tag_table (s, tag, h) ;
  
  *m = *w = 0 ;
  if (cover) *cover = 0 ;
  if (tt && tt->rows && tt->cols >= 3)
    {
      int ir ;
      for (ir = 0 ; ir < tt->rows ; ir++)
	{
	  const char *ccp = ac_table_printable (tt, ir, 0, 0) ;
	  if (ccp && !strcasecmp (ccp, run))
	    break ;
	}
      if (ir < tt->rows)
	{
	  ok = TRUE ;
	  *m = ac_table_int (tt, ir, 1, 0) ;
	  *w = ac_table_int (tt, ir, 2, 0) ;
   	  if (cover)
	    {
	      *cover = ac_table_int (tt, ir, 3, 0) ;
	      if (*cover < *m + *w)
		*cover = *m + *w ;
	    }
	}
    }
  return ok ;
} /* snpVcfCounts */

/*************************************************************************************/

static void snpVcfRunCounts (vTXT txt, SNP *snp, AC_OBJ s, const char *run)
{
  int mutant, wild, cover, mm, mp, wm, wp ;
  float frequency = 0 ;

  snpVcfCounts (s, "nsCounts", run, &mutant, &wild, &cover) ;
  snpVcfCounts (s, "fCounts", run, &mp, &wp, 0) ;
  snpVcfCounts (s, "rCounts", run, &mm, &wm, 0) ;
  
  if (cover)
    { 
      AC_HANDLE h = ac_new_handle () ;
      AC_TABLE tt = ac_tag_table (s, "Phase", h) ;
      int strand =  tt ? ac_table_int (tt, 0, 0, 0) : 0 ;
      const char *genotype, *phaseName = tt ? ac_table_printable (tt, 0, 1, ".") : "." ;

      frequency =  100.0 * mutant/cover ;
      if (strand)
	{
	  if (frequency > 80)
	    genotype = "1|1" ;
	  else if (frequency < 20)
	    genotype = "0|0" ;
	  else if (strand == 1)
	    genotype = "1|0" ;
	  else if (strand == 2)
	    genotype = "0|1" ;
	}
      else
	{
	  if (frequency > 80)
	    genotype = "1/1" ;
	  else if (frequency < 20)
	    genotype = "0/0" ;
	  else 
	    genotype = "1/0" ;
	}

      vtxtPrintf (txt, "\t%s:%.3f:%d:%d,%d:%d,%d:%d,%d:%s"
		  , genotype
		  , frequency/100
		  , cover
		  , wild, mutant
		  , wp, mp
		  , wm, mm
		  , phaseName
		  ) ;

      ac_free (h) ;
    }
  else
    {
      vtxtPrintf (txt, "\t.:.:.:.:.") ;
    }
} /*  snpVcfRunCounts */

/*************************************************************************************/

static void snpVcfTableHeader (SNP *snp, ACEOUT ao, KEYSET runs)
{
  AC_HANDLE h = ac_new_handle () ;
  AC_TABLE tt = 0 ;
  int ir ;
  const char *errors = 0 ;

  aceOutf (ao, "%s", 
	   "##fileformat=VCFv4.2\n"
	   "##fileDate=20170912\n"
	   "##source=MAGIC_NCBI.v2017_09_11\n"
	   "##reference=GRCh37.p10_Primary_Assembly_GCF_000001405.22_https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.22\n"
	   ) ;
  tt = ac_bql_table (snp->db , "select m, t, ln from m in ?map, ln in m->length, t in m->title where t", 0, 0, &errors, h) ;
  for (ir = 0 ; tt && ir < tt->rows ; ir++)
    {
      const char *ccp = ac_table_printable (tt, ir, 0, "") ;
      
      if (ccp && strlen (ccp) < 3)
	{
	  aceOutf (ao, "##contig=<ID=%s,chromosome=chr%s,RefSeqID=%s,length=%d>\n"
		   , ac_table_printable (tt, ir, 0, "") 
		   , ac_table_printable (tt, ir, 0, "") 
		   , ac_table_printable (tt, ir, 1, "") 
		   , ac_table_int (tt, ir, 2, 0) 
		   ) ;
	}
    }
  aceOutf (ao, "%s", 
	   "##ALT=<ID=NON_REF,Description=\"Represents a particular alternative allele at this location, a single alternative per line\">\n"
	   "##FILTER=<ID=Reject_M,Number=0,Type=Flag,Description=\"monomodal\">\n"
	   "##FILTER=<ID=Reject_P,Number=0,Type=Flag,Description=\"Not consistently phasable\">\n"
	   "##INFO=<ID=In_51bp_genomic_repeat,Number=1,Type=Integer,Description=\"The modified base is at the centre of a 51 bp segment repeated N times on the either strand of the genome.\">\n"
	   "##INFO=<ID=In_101bp_genomic_repeat,Number=1,Type=Integer,Description=\"The modified base is at the centre of a 101 bp segment repeated N times on the either strand of the genome.\">\n"
	   "##INFO=<ID=MOD,Number=0,Type=Flag,Description=\"This variant appears to be MOD specific, it is not detected by MAGIC or by GIAB in sample NA12878.\">\n"
	   "##INFO=<ID=Monomodal,Number=0,Type=Flag,Description=\"This variant is hard to read, possibly it sits in a genome compression or is generated by a systematic sequencer error: the histogram of the observed frequency of this variant in hundreds of samples is monomodal in the intermediate range, without clear homozygotes.\">\n"
	   "##INFO=<ID=Phasing,Number=1,Type=Integer,Description=\"Number of heterozygous SNVs consistently phased with this  heterozygous SNV in cis or in trans,  by at least 3 paired end fragments with less than 5% contradictions\">\n"
	   "##INFO=<ID=Non_phasable,Number=1,Type=Integer,Description=\"Number of heterozygous SNVs with inconsistent phasing with this SNV: the SNVs are seemingly seen associated in both cis and trans. Either there is a copy number variation and more than 2 copies, or the SNV is an artefact\">\n"
	   "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
	   "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth, only read pairs aligning uniquely and over their quasi-entire length contribute ton SNV calling\">\n"
	   "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Variant allele frequency for the ref and alt alleles in the order listed\">\n"
	   "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n"
	   "##FORMAT=<ID=ADP,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed, aligned on strand plus of the genome; this is used to estimate strand bias\">\n"
	   "##FORMAT=<ID=ADM,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed, aligned on strand minus of the genome; used to estimate strand bias\">\n"
	   "##FORMAT=<ID=PS,Number=1,Type=String,Description=\"Phase information: all heterozygous alleles in this phase group are measured to be in cis or trans relative to each other, in at least three connecting pairs of reads\">\n"
	   "##phasing=partial\n"
	   ) ;
  
  aceOutf (ao, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT") ;
    
  if (snp->runs)
    {
      DICT *dict = snp->runDict ;
      KEYSET runs = snp->runs ;
      int ir, irMax = keySetMax (runs) ; 
      
      for (ir = 0 ; ir < irMax ; ir++)
	{
	  int run = keySet (runs, ir) ;
	  aceOutf (ao, "\t%s", dictName (dict, run)) ;
	}
    }

  aceOutf (ao, "\n") ;
  
  return ;
} /* snpVcfTableHeader  */

/*************************************************************************************/
/* parse RunsSortedList */
static KEYSET snpVcfTableGetRuns (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (snp->runListFileName, FALSE, h) ;
  int ir = 0 ;

  snp->runs = keySetHandleCreate (snp->h) ;
  while (aceInCard (ai))
    {
      char *cp = aceInWord (ai) ;
      int run ;

      if (!cp || ! *cp || *cp == '#' || *cp == '/')
	continue ;

      dictAdd (snp->runDict, cp, &run) ;
      keySet (snp->runs, ir++) = run ;
    }
  ac_free (ai) ;
  ac_free (h) ;

  return snp->runs ;
} /* snpVcfTableGetRuns */

/*************************************************************************************/

static void snpVcfTable (SNP *snp)
{
  AC_HANDLE h1 = 0, h = ac_new_handle () ;
  ACEOUT ao = 0 ;
  int nn = 0 ;
  AC_DB db = snp->db ;
  AC_OBJ s = 0 ;
  AC_ITER iter = ac_query_iter (db, TRUE, hprintf (h, "find variant project == %s", snp->project), 0, h) ;
  vTXT txt = vtxtHandleCreate (h) ;
  /*
    const char **runp, *runs[] = {"Magic11", 0, "MagicWG2", "GIAB", "7_18.FDA.Walia", "8_24.Atgenomix.conf30", "8_30.NYU_Langone", "7_27.NIH_Lack", "8_31.Hard", "Richa", 0} ;
  */
  KEYSET runs = snpVcfTableGetRuns (snp) ;
  
  ao = aceOutCreate (snp->outFileName, ".vcf", snp->gzo, h) ;
  snpVcfTableHeader (snp, ao, runs) ;
  
  while (ac_free (s), ac_free (h1), vtxtClear (txt), s = ac_iter_obj (iter))
    { 
      h1 = ac_new_handle () ;
      char *cq, *cp = strnew (ac_name(s), h1) ;
      
      vtxtClear (txt) ;
      /* chrom */
      cq = strchr (cp, ':') ; if (! cq) continue ; 
      *cq = 0 ; vtxtPrintf (txt, "%s", cp) ; cp = cq + 1 ; 
      if (! *cp) continue ;
      /* pos */
      cq = strchr (cp, ':') ; if (! cq) continue ; 
      *cq = 0 ; vtxtPrintf (txt, "\t%s", cp) ; cp = cq + 1 ; 
      if (! *cp) continue ;
      vtxtPrintf (txt, "\t.") ; /* RS id */
      /* reference */
      cq = strchr (cp, ':') ; if (! cq) continue ; 
      *cq = 0 ; bufferToUpper (cp) ; vtxtPrintf (txt, "\t%s",cp ) ; cp = cq + 1 ; 
      if (! *cp) continue ;
      /* variant */
      bufferToUpper (cp) ; vtxtPrintf (txt, "\t%s", cp) ; cp = cq + 1 ; 
      
      /*
      if (! snpVcfCounts (s, "nsCounts", snp->run, &mutant, &wild, &cover))
	continue ;
      if (cover < 10) continue ;
      if (mutant == 0) continue ;
      frequency =  100.0 * mutant/cover ;
      if (frequency < 20) continue ;
      */

      vtxtPrintf (txt, "\t.") ; /* quality */
      if (1)
	{
	  int phase = ac_tag_int (s, "Phase", 0) ;
	  
	  if (ac_has_tag (s, "Rejected"))
	    vtxtPrintf (txt, "\tReject_M") ; /* filter */
	  else if (phase < 0)
	    vtxtPrintf (txt, "\tReject_P") ; /* filter */ 
	  else if (0 && ac_has_tag (s, "Genome_51mer_repeats"))
	    vtxtPrintf (txt, "\tReject_R") ; /* filter */ 
	  else
	     vtxtPrintf (txt, "\tPASS") ; /* filter */
	}
      
      vtxtPrintf (txt, "\t") ;
      if (1)
	{  /* info */
	  char *sep = "" ; 
	  int x ;
	  int phase = ac_tag_int (s, "Phase", 0) ;
	  const char *coding ;
	  
	  if (0 && ! ac_has_tag (s, "MagicWG2") && ! ac_has_tag (s, "GIAB"))
	    { vtxtPrintf (txt, "%sMOD", sep) ; sep = ";" ; }
	  x = ac_tag_int (s, "Genome_51mer_repeats", 0) ;
	  if (x > 0)
	    { vtxtPrintf (txt, "%sIn_51bp_genomic_repeat=%d", sep, x) ; sep = ";" ; }
	  x = ac_tag_int (s, "Genome_101mer_repeats", 0) ;
	  if (x > 0)
	    { vtxtPrintf (txt, "%sIn_101bp_genomic_repeat=%d", sep, x) ; sep = ";" ; }
	  if (ac_has_tag (s, "Rejected"))
	    { vtxtPrintf (txt, "%sMonomodal", sep) ; sep = ";" ; }
	  if (phase > 0)
	    { 
	      AC_TABLE tt = ac_tag_table (s, "Phase", h1) ;
	      vtxtPrintf (txt, "%sPhasing=%d", sep, ac_table_int (tt, 0, 2, 0)) ; sep = ";" ; 
	      ac_free (tt) ;
	    }
	  if (phase < 0)
	    { vtxtPrintf (txt, "%sNon_phasable=%d",sep, -phase) ; sep = ";" ; }
	  coding= ac_tag_printable (s, "pName", 0) ;
	  if (0 && coding)
	    { 
	      /* 
		 AC_TABLE tt1 = ac_tag_table (s, "Reference_sequence", h1) ;
		 AC_TABLE tt2 = ac_tag_table (s, "Observed__sequence", h1) ;
		 vtxtPrintf (txt, "Coding=\"") ; sep = ";" ; 
	      */	    
	    }
	}
	   
      /* INFO : proprietes du SNP, distance to nearest neighbour, difficult region, nb de 51mers du genome, delins insert sub deletes dup dim strand discordant  */
      vtxtPrintf (txt, "\tGT:AF:DP:AD:ADP:ADM:PS") ;
      /*
	genotype = "." ;
	     if (frequence > 90) genotype = "1/1" ;
	     else  if (frequence < 5) genotype = "1/1" ;
	     else  genotype = "1/0" ;
	     vtxtPrintf (txt, "\t%s", genotype) ;
      */
        
      if (snp->runs)
	{
	  DICT *dict = snp->runDict ;
	  KEYSET runs = snp->runs ;
	  int ir, irMax = keySetMax (runs) ; 
	  
	  for (ir = 0 ; ir < irMax ; ir++)
	    {
	      int run = keySet (runs, ir) ;
	      const char* runName = dictName (dict, run) ;
	      
	      snpVcfRunCounts (txt, snp, s, runName) ;
	    }
	}      
      if (0) vtxtPrintf (txt, "\t%s", ac_has_tag (s, "GIAB_easy") ? "EASY" : "DIFFICULT") ;
      
      vtxtPrintf (txt, "\n") ;
      
      nn++ ;
      aceOut (ao, vtxtPtr (txt)) ;
    }
  
  
  fprintf (stderr, "// Exported %d snp in file %s\n", nn, aceOutFileName (ao)) ;
  ac_free (h) ;
  return ;
} /*  snpVcfTable */

/*************************************************************************************/
/*************************************************************************************/

typedef struct gpStruct {int gene, nCells, nSnps, type[7], nPcCells, nPcSnps, pcType[7] ; int *runType ; float *runScore ; float *runPcScore ; } GP ;
static int snpPrevalenceTable (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = snp->ai ;
  ACEOUT ao = 0 ;
  int line = 0, nn = 0, gene = 0, run = 0, snvTypeCol = 0, alFreqCol = 0 ;
  KEYSET col2run = keySetHandleCreate (h) ; 
  AC_TABLE runs =  snpProjectRuns (snp, h) ;
  int runMax = 0, typeMax = 6 ;
  GP *gp ;
  Array genes = arrayHandleCreate (50000, GP, h) ;
  DICT *geneDict = dictHandleCreate (50000, h) ;
  DICT *runDict = dictHandleCreate (500, h) ;

  /*
   * lister les comptages par gene dans une struct
   * sort
   * cumul par gene
   * report 
   */

  if (1)
    {
      int ir ;
      const char *ccp ;
      for (ir = 0 ; ir < runs->rows ; ir++)
	{
	  ccp = ac_table_printable (runs, ir, 0, 0) ;
	  if (ccp)
	    dictAdd (runDict, ccp, 0) ;
	}
    }
  aceInSpecial (ai, "\n") ;
  
  /* parse the title line */
  while (aceInCard (ai))
    {
      line++ ;
      if (line == 2) /* grab the list of runs */
	{
	  int i ; char cutter = 0 ;
	  const char *ccp = aceInWordCut (ai, "\t", &cutter) ; 
	  if (! ccp || strcmp (ccp, "# Chromosome"))
	    {
	      fprintf (stderr, "ERROR : Aborting  snpPrevalenceTable while parsing line %d of file %s", line, aceInFileName (ai)) ;
	      fprintf (stderr, "ERROR : Expecting column 1 line 2 of input to be \"# Chromosome\"\n") ;
	      goto done ;
	    }
	  for (i = 2 ;  ; i++)
	    {
	      cutter = 0 ;
	      ccp =  aceInWordCut (ai, "\t", &cutter) ;
	      if (!ccp && ! cutter) break ;
	      if (!strcmp (ccp, "SNV type"))
		{ snvTypeCol = i ; continue ; }
	      if (!strcmp (ccp, "Allele frequency in cohort"))
		{ alFreqCol = i ; continue ; }

	      if (ccp && dictFind (runDict, ccp, &run))
		{
		  keySet (col2run, i) = run ;
		  if (runMax <= run)
		    runMax = run + 1 ;
		}
	    }
	  i = keySetMax (col2run) ;
	  fprintf (stderr, "Recognised %d run names in file %s\n", i, aceInFileName (ai)) ;
	  if (i < 1) goto done ;
	  if (! snvTypeCol)
	    messcrash ("In line 2 of file %s , there was no column called \"SNV type\" needed to distinguish protein changing SNPs", aceInFileName (ai)) ;
	  break ;
	}
    }
  if (line < 2)
    goto done ;
  
  while (aceInCard (ai))
    {
      int i; char cutter = 0 ;
      const char *ccp = aceInWord (ai) ;
      BOOL isProteinChanging = FALSE ;
      BOOL isMajorAllele = FALSE ;
      line++ ;
      if (! ccp || *ccp == '#') continue ;

	/* eliminate monomodal snps */
	ccp = aceInPos (ai) ;
	if (snp->dropMonomodal && ccp)
	  {
	    AC_HANDLE h3 = ac_new_handle () ;
	    ACEIN card = aceInCreateFromText (ccp, 0, h3) ;
	    int i, t = 0, kk[8], kkk = 0 ;
	    float z = 0 ;
	    aceInSpecial (card, "\n") ;
	    aceInCard (card) ;
	    
	    memset (kk, 0, sizeof(kk)) ;
	    for (i = 1 ; i < keySetMax (col2run) ; i++)
	      {
		ccp = aceInWordCut (card, "\t", 0) ;
		if (! ccp || *ccp == '#') continue ;

		/*
		  if (i>=1 && i <= 14 && line>=5923 && line <= 5927)
		  fprintf(stderr,"..%d.%d...%s\n",line,i, ccp) ;
		  if (line > 5930) exit(0) ;
		*/
		t = 0 ;
		if (i == alFreqCol)
		  {
		    if (sscanf(ccp,"%f",&z) == 1 && z > 50)
		      isMajorAllele = TRUE ;
		    continue ;
		  }
		if (i == snvTypeCol)
		  {
		    /* fprintf(stderr,"..%d.%d...%s\n",line,i, ccp) ; */
		    if (!strcmp (ccp, "Synonymous") ||
			!strcmp (ccp, "Non_coding_transcript") ||
			strstr (ccp, "UTR") ||
			!strcmp (ccp, "Frame_preserving_indel") ||
			! strcmp (ccp, "-")
			)
		      continue ;
		    if (!strncmp (ccp, "AA_sub", 6) ||
			!strncmp (ccp, "Met1", 4) ||
			!strncmp (ccp, "Ter", 3) ||
			!strncmp (ccp, "Frame", 5)  ||
			!strncmp (ccp, "Extension", 9)  ||
			!strncmp (ccp, "AA_to_Stop", 10) ||
			!strncmp (ccp, "Stop_to_AA", 10)
			)
		      { isProteinChanging = TRUE ; /* fprintf(stderr,"....%s\n",ccp) ; */ }
		  }
		if (keySet (col2run, i))
		  {
		    if (sscanf(ccp,"%f",&z) == 1)
		      {
			if (z >= 95) 
			  { t = 5 ;  }
			else  if (z >= 80) 
			  { t = 4 ;  }
			else  if (z > 20) 
			  { t = 3 ;  }
			else  if (z > 5) 
			  t = 2 ;
			else  if (z >= 0) 
			  t = 1 ;
			else
			  t = 0 ;
		      }
		    kk[t]++ ; if (t) kkk++ ;
		  }
	      }
	    ac_free (h3) ;
	    
	    /* if monomodal, jump this snp */
	    if (snpPleaseDropMonomodal (kk[1], kk[2], kk[3], kk[4], kk[5]))
	      continue ;
	    
	    /*  this would reject pure wild type
	      if (kk[2] > kk[3] && kk[3] > kk[4] && kk[4] > kk[5])
	      continue ;
	    */
	  }

	for (i = 1 ; i <= 9 ; i++)
	  { aceInWordCut (ai, "\t", &cutter) ; }
	i = 10 ; ccp = aceInWordCut (ai, "\t", &cutter) ; /* gene */
	if (! ccp) continue ;
	dictAdd (geneDict, ccp, &gene) ;
	gp = arrayp (genes, gene, GP) ;
	gp->gene = gene ; /* we need this if we later sort the table */
	if (0 && ! gp->runType) gp->runType = (int *)halloc (runMax * typeMax * sizeof (int), h) ;
	if (! gp->runScore) 
	  {
	    gp->runScore = (float *)halloc (runMax * sizeof (float), h) ;
	    for (i = 0 ; i < runMax ; i++)
	      gp->runScore[i] = -10 ; 
	  }
 	if (! gp->runPcScore) 
	  {
	    gp->runPcScore = (float *)halloc (runMax * sizeof (float), h) ;
	    for (i = 0 ; i < runMax ; i++)
	      gp->runPcScore[i] = -10 ; 
	  }
	gp->nSnps++ ;
	if (isProteinChanging)	
	  gp->nPcSnps++ ;
	for (i = 11 ;  i < keySetMax (col2run) ; i++)
	  {
	    /* get value z in each run 
	     * split value into a dosage
	     * foreach group to which the run belongs, cumul the counts per dosage
	     */
	    char *ccp = aceInWordCut (ai, "\t", &cutter) ; /* ignore */

	    run = keySet (col2run, i) ;
	    if (ccp && run)
	      {
		float z ; int t = 0 ; 
		
		if (sscanf (ccp, "%f", &z) == 1)
		  {
		    if (isMajorAllele && z >= 0)
		      z = 100 - z ;
		    
		    if (z >= 0 && gp->runScore[run] == -10) gp->runScore[run] = 0 ;
		    if (z >= 0 && isProteinChanging && gp->runPcScore[run] == -10) gp->runPcScore[run] = 0 ;
		    if (z >= 95) 
			{
			  t = 5 ; gp->runScore[run] += 2 ;   
			  if (isProteinChanging)
			      gp->runPcScore[run] += 2 ;   
			}
		      else  if (z >= 80) 
			{  
			  t = 4 ; gp->runScore[run] += 1.5 ; 
			  if (isProteinChanging)
			    gp->runPcScore[run] += 1.5 ;   
			}
		      else  if (z > 20) 
			{ 
			  t = 3 ; gp->runScore[run] += 1.0 ; 
			  if (isProteinChanging)
			    gp->runPcScore[run] += 1.0 ;   
			}
		      else  if (z > 5) 
			t = 2 ;
		      else  if (z >= 0) 
			t = 1 ;
		      else
		       t = 0 ;
		    }
		
		  gp->nCells++ ;
		  gp->type[t]++ ;
		  if (isProteinChanging)
		    {
		      gp->nPcCells++ ;
		      gp->pcType[t]++ ;
		    }
		  if (t > typeMax) messcrash ("t = %d > typeMax = %d\n", t, typeMax) ;
		  /* gp->runType [run * typeMax + t]++ ; */
		}
	  }
      }
  
  ao = aceOutCreate (snp->outFileName, ".snp_prevalence_per_gene.txt", snp->gzo, h) ;
  aceOutDate (ao, "###", snp->project) ;
  aceOutf (ao, "# Gene") ;
  /*
  const char *types[] = { "%s NA"
			  , "%s Reference (<5% variant)"
			  , "%s Low (5-20%)"
			  , "%s Mid (20-80%)"
			  , "%s High (80-95%)"
			  , "%s Pure variant (>95%)"
  } ;
  */
  aceOut (ao, "\tSNP sites and types\tObservations\tNA\tReference (<5% variant)\tLow (5-20%)\tMid (20-80%)\tHigh (80-95%)\tPure variant (>95%)\tMid to Pure (20-100%)") ;
  aceOut (ao, "\t\t%%NA\t%%Reference (<5% variant)\t%Low (5-20%)\t%Mid (20-80%)\t%High (80-95%)\t%Pure variant (>95%)\t% Mid to Pure (20-100%)") ;
    aceOutf (ao, "\t") ;
  {
    aceOutf (ao, "\tRuns with no measured site") ;
    aceOutf (ao, "\tRuns with no mutant allele") ;
    aceOutf (ao, "\tRuns with at least one mutant allele") ;  /* new */
    aceOutf (ao, "\tRuns with one mutant allele") ;
    aceOutf (ao, "\tRuns with 2 mutant alleles") ;
    aceOutf (ao, "\tRuns with 3 or more mutant alleles") ;
    aceOutf (ao, "\tRuns with 10 or more mutant alleles") ;
  }

  aceOutf (ao, "\t\tProtein changing SNP sites and types\tProtein changing observations\tNA\tProtein changing Reference (<5% variant)\tProtein changing Low (5-20%)\tProtein changing Mid (20-80%)\tProtein changing High (80-95%%)\tProtein changing Pure variant (>95%)\tProtein changing Mid to Pure (20-100%%)") ; /* new */
  aceOutf (ao, "\t\t%% NA\t%% Protein changing Reference (<5% variant)\t%%Protein changing Low (5-20%)\t%Protein changing %Mid (20-80%)\t%%Protein changing High (80-95%%)\t%%Protein changing Pure variant (>95%)\t%%Protein changing Mid to Pure (20-100%%)") ; /* new */
    aceOutf (ao, "\t") ;
  {
    aceOutf (ao, "\tRuns with no measured protein changing site") ;
    aceOutf (ao, "\tRuns with no observed protein changing mutant allele") ; /* new */
    aceOutf (ao, "\tRuns with at least one protein changing mutant allele") ;
    aceOutf (ao, "\tRuns with one protein changing mutant allele") ;
    aceOutf (ao, "\tRuns with 2 protein changing mutant alleles") ;
    aceOutf (ao, "\tRuns with 3 or more protein changing mutant alleles") ;
    aceOutf (ao, "\tRuns with 10 or more protein changing mutant alleles") ;
  }
  aceOutf (ao, "\t") ;
  {
    int j, run ;
    for (j = 0 ; j < keySetMax(col2run) ; j++)
      {
	run = keySet (col2run, j) ;
	if (run)
	  aceOutf (ao, "\t%s", dictName(runDict, run)) ;
      }
  }
  
  if (1)
    {
      int i, j, t, k, pass ;
      for (i = 0, gp = arrp (genes, i, GP) ; i < arrayMax (genes) ; i++, gp++)
	{
	  if (! gp->nCells)
	    continue ;
	  if (gp->type [3] + gp->type [4] + gp->type [5] < 1)
	    continue ;
	    
	  aceOutf (ao, "\n%s", dictName (geneDict, gp->gene)) ;
	  /*
	    chrom a1 a2  CF geneexpressionindex table
	    type
            titre
	  */
	  for (pass = 0 ; pass < 2 ; pass++)
	    {
	      if (pass == 0)
		{
		  aceOutf (ao, "\t%d\t%d", gp->nSnps, gp->nCells) ;
		  for (t = 0  ; t < 6 ; t++)
		    aceOutf (ao, "\t%d", (int)gp->type [t]) ;
		  aceOutf (ao, "\t%d", gp->type [3] + gp->type [4] + gp->type [5]) ;
		  aceOutf (ao, "\t") ;
		  for (t = 0  ; t < 6 ; t++)
		    aceOutf (ao, "\t%.1f", (float)gp->type [t]*100.0/gp->nCells) ;
		  aceOutf (ao, "\t%.1f",(float)(gp->type [3] + gp->type [4] + gp->type [5])*100.0/gp->nCells ) ;
		}
	      else
		{
		  int nc = gp->nPcCells ;

		  aceOutf (ao, "\t") ;
		  if (nc == 0) nc = 1 ;
		  aceOutf (ao, "\t%d\t%d", gp->nPcSnps, gp->nPcCells) ;
		  for (t = 0  ; t < 6 ; t++)
		    aceOutf (ao, "\t%d", (int)gp->pcType [t]) ;
		  aceOutf (ao, "\t%d", gp->pcType [3] + gp->pcType [4] + gp->pcType [5]) ;
		  aceOutf (ao, "\t") ;
		  for (t = 0  ; t < 6 ; t++)
		    aceOutf (ao, "\t%.1f", (float)gp->pcType [t]*100.0/nc) ;
		  aceOutf (ao, "\t%.1f",(float)(gp->pcType [3] + gp->pcType [4] + gp->pcType [5])*100.0/nc) ;
		}
	      aceOutf (ao, "\t") ;
	      for (k = -1 ; k < 6 ; k++)
		{
		  int n = 0 ;
		  for (j = 0 ; j < keySetMax(col2run) ; j++)
		    {
		      run = keySet (col2run, j) ;
		      if (run)
			{
			  float f ;
			  
			  if (pass == 0) 
			    f = gp->runScore ? gp->runScore[keySet (col2run,j)] : -10 ; 
			  else
			    f = gp->runPcScore ? gp->runPcScore[keySet (col2run,j)] : -10 ; 
			  switch (k)
			    {
			    case -1:
			      if (f == -10)
				n++ ;
			      break ;
			    case 0:
			      if (f == 0)
				n++ ;
			      break ;
			    case 1:  /* at least 1 */
			      if (f > 0)
				n++ ;
			      break ;
			    case 2:  /* exactly 1 */
			      if (f >= 1 && f < 1.8)
				n++ ;
			      break ;
			    case 3:  /* exactly 2 */
			      if (f >= 1.8 && f < 3)
				n++ ;
			      break ;
			    case 4: /* 3 or more */
			      if (f >= 3)
				n++ ;
			      break ;
			    case 5:   /* 10 or more */
			      if (f >= 10)
				n++ ;
			      break ;
			    }
			}
		    }
		  aceOutf (ao, "\t%d", n) ;
		}
	    }
	  aceOutf (ao, "\t") ;
	  for (j = 0 ; j < keySetMax(col2run) ; j++)
	    {
	      run = keySet (col2run, j) ;
	      if (run)
		aceOutf (ao, "\t%.1f", gp->runScore[keySet (col2run,j)]) ;
	    }
	}
    }
  aceOutf (ao, "\n") ;
  
 done:
  ac_free (h) ;
  return nn ;
} /*  snpPrevalenceTable */

/*************************************************************************************/
/*************************************************************************************/

static int snpParseSelected8kbTargetFile (SNP *snp)
{
  AC_HANDLE h = ac_new_handle () ;
  int nn = 0 ;
  const char *ccp ;
  ACEIN ai = aceInCreate (snp->selected8kbFileName, 0, h) ;

  if (! ai)
    messcrash ("Sorry, i cannot find the -selected8kbList file : %s", snp->selected8kbFileName) ;

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (ccp && *ccp != '#')
	dictAdd (snp->selected8kbDict, ccp, &nn) ;
    }

  ac_free (h) ;
  return nn ;
} /* snpParseSelected8kbTargetFile */

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (message)
    {
      fprintf (stderr, "// program snp : %s\n%s\n// Try snp -h or --help\n", timeShowNow(), message) ;
      exit (2) ;
    }

  fprintf  (stderr,
	    "// SNP.c:\n"
	    "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2009, mieg@ncbi.nlm.nih.gov\n"
	    "//   This successive phases of this program are automatically invoked by the MAGIC SNP pipeline\n"
	    "// Purpose\n"
	    "//   Analyze the mismatches explicitly reported in the best .hits files exported by the MAGIC aligner\n" 
	    "//   in view of identifying reliable sequence variants, SNPs or short indels\n"
	    "//   and exporting large detailed counting tables analogous to .vcf files as well as concise\n"
	    "//   user friendly allele frequency tables, fully annotated with metadata on runs and SNVs\n"
	    "//   When relevant, the geneindex program can then perform a SNV based differential analysis\n"
	    "//   a covariance analysis or consctruct classifiers, or predictors\n" 
	    "//\n"
	    "//   The program works in 3 successive phases, automatically invoked by the MAGIC pipeline\n"
	    "//   Phase 1: Import a set of .hits[.gz] files and create a .BRS[.gz] files, optionally gzipped\n"
	    "//            The BRS file counts at each position the support for the reference allele and all variants\n"
	    "//              carefully counting as support for indels only the reads bridging over repeated bases\n"
            "//   Phase 2: Import a BRS file and export a .snp file, one line per SNP\n"
            "//            2a: detect in the BRS file all SNPs satisfying the quality thresholds\n"
	    "//            2b: list the union of all SNPs, detected by phase 2a in any run\n" 
	    "//            2c: count in the BRS file the support in all runs of all SNPs present in this list\n"
	    "//                in this way, any SNP detected in even just one sample, is quantified in all samples, including homozygous reference\n"
	    "//   Phase 3: Import a .snp file and analyze it\n"
	    "//            3a: export a .ace file to be loaded in a VariantDB acedb database aware of gene structure\n"
	    "//            3b: remap the mRNA SNVs on the genome and vice-versa (using VariantDB)\n"
	    "//            3d: report various user-oriented tables and histograms, of biological imterest\n"
	    "// Input/output\n"
	    "//   -i input_file  [-gzi] : default stdin,  if the file is called .gz or if -gzi is specified, the input is automatically gunzipped\n"
	    "//   -o output_file [-gzo] : default stdout, if -gzo, the output files are gzipped\n"
	    "//            each output file is named as 'output_file'\n"
	    "//            followed by action specific endings constructed by the program, often an action produces several files [optionally .gz if -gzo is requested]\n"
	    "//\n"
	    "// Phase 1 actions: Transform the .hits files into .BRS files\n"
	    "//   -run run_name : name of the run, mandatory \n"
	    "//   -ventilate -select selectFile\n"
	    "//      ventilate the .hits files into a collection of files named by zone defined in  selectFile\n"
	    "//      This phase is useful to keep file and memory requirements under control\n"
	    "//      the selectFile has 4 tab separated columns : t x1 x2 zone\n"
	    "//      t is a target name, column 11 in the .hits file: usually a chrom or a transcript\n"
	    "//      x1 x2 are coordinates on this target giving the extent of the named zone\n"
	    "//      This program allows to select a zone (targeted exome experiment)\n"
	    "//      it can also be used to 'transpose' the hits file, organized by sequencing fragments\n"
	    "//      into a set of say 20 to 200 geographic zones, or gene sets zones, before using -count\n"
	    "//   -hits2BRS -run run_name -fasta f -strategy s [-runQuality f] [-select selectFile] [-minAli n] [-minAliPerCent n] [-dropMultiplicity]:\n"
	    "//       report in BRST format the coverage and errors found in a .hits file exported by clipalign\n"
            "//     -strategy [ Exome | Genome | RNA_seq ]    default:RNA_seq\n"
	    "//        should match the choice in clipalign, influences the coverage counts\n"
	    "//     -fasta : file f should contain, in fasta format. the DNA of the target \"reference\" sequences\n"
	    "//     -select : only report hits against the selected zones \n"
	    "//     -minAli : discard hits shorter than the specified limit in bp\n"
	    "//     -minAliPercent :  discard reads not aligning on a n%% of their clipped lenght and ali < 140\n"
	    "//       the optional -runQuality file is described below\n"
	    "//       The output, BRST format gives for each position in each target the reference and variant counts\n"
	    "//        distinguishing B (base count) R(repeat count) S(substitution count) T(transition count)\n"
	    "//        These for numbers are needed to correctly estimate the wild type counts of Insertions/deletions\n"
	    "//        and the transition counts in the SOLiD transition based sequencing\n"
	    "//     -dropMultiplicity : count only once identical fragment, rather that using the #multiplicity of the reads\n"
	    "//   -runQuality file_name :  Optional filtering on averaged positional fastq quality\n"
	    "//        3 column file: p x q\n"
	    "//        where p is 0 or the part 1 or 2 of a paired end read, named /1 /2 in column 1\n"
	    "//        where x is the position in the read, i.e. the sequencing machine cycle\n"
	    "//        and q the average fastq quality at this position: -10 log10(mismatch probability)\n"
	    "//        In MAGIC, the original fastq quality coefficients are judged unreliable and discarded\n" 
	    "//        rather the averaged q is estimated automatically after alignment and stored in MetaDB\n"
	    "//        effect: the numbers reported by -count are multiplied by 100 - 5 * (30 - q)/30, bounded to [0,100]\n"
	    "//   -BRS_merge -run group_name :  export a single data file combining all the entries\n"
	    "//      The input should be a concatenation of BRST files exported using the -hits2BRS option\n"
	    "//\n"
	    "// Phase 2 actions: Transform in 2 passes, BRS_detect, BRS_count, the .BRS files into .snp files\n"
	    "//    -run run_name : name of the run, mandatory \n"
	    "//    -BRS_detect  : Phase 2a, list all SNP which pass the following thresholds\n"
	    "//    -BRS_Detect  : optional, export immediatly in the tsf format defined below\n"
	    "//      The input file should be a [concatenation of] BRS[.gz] file created by hits2BRS\n"
	    "//      -minFrequency f : [default 20] f is a integer number between 0 and 100\n"
	    "//      -minMutant m : [default 4] m is an integer >= 0\n"
	    "//      -minCover  c : [default 10] c is an integer > 0\n"
	    "//        Only export variants seen at least m times, with coverage at least c and frequency at least f%%\n"
	    "//   -BRS_count -snp_list f -run runName : Phase 2b, export in .tsf format the counts for a requested list of SNPs\n"
	    "//      The input file should be a [concatenation of] BRS[.gz] file created by hits2BRS\n"
	    "//      The mandatory runName will be included as column 2 in each exported line\n" 
	    "//      The mandatory file f contains the list exported by BRS_detect\n"
	    "//      The listed snps are exported if covered (c>0) even if there is no mutant (m=0)\n"
	    "//      -solid : work in transition-space, report the uncorrected transition errors on both side of the SNP as (oo1:oo2)\n" 
	    "//      The tsf format has 9 columns: snpName runName 6i cp mp rp cm mm rm  (cover mutant reference for both strands)\n"
	    "//\n"
	    "// Phase 3 actions: Analyse the .snp files\n"
	    "//   The analyses rely on a the existence of an acedb VariantDB.$zone database, aware of genes, transcripts and coding structures\n"
	    "//     and containing a copy of the metadata of the runs, originally hand constructed in MetaDB\n"
	    "//     The parameter -db points to this database, we recommend one database per zone, allowing parallelization\n"
	    "//   -db_remap2genome tmp/METADATA/mrnaRemap.gz  -db ACEDB\n"
	    "//      Remap the transcript variants into genome coordinates\n"
	    "//   -db_remap2genes tmp/METADATA/mrnaRemap.gz -db ACEDB\n"
	    "//      Remap the genome variants into transcript coordinates\n"
	    "//   -db_count -i count_file -db ACEDB  : scan the input file and adjust in the ACEDB database the variant->population and ->strand counts\n"
	    "//   Finally various kinds of user friendly reports, if -db is documented, the tables will show run->sample, title and other info\n"
	    "//   -db_frequencyTable  -i count_file  -db ACEDB -project project [-select selectFile] [-minCover c] [-Reference_genome] [-dropMonomodal]\n"
	    "//       Report a long table, one SNP per line, one run per column, and several histograms\n"
	    "//       Use minCover and selectFile to restrict the report to a zone or a cover threshold\n"
	    "//       Use -Reference_genome to qualify the g.name in the exported table\n"
	    "//   -db_frequencyHisto   -db ACEDB -project project [-select selectFile] [-minCover c] \n"
	    "//       Report a smoothed prevalence histo, for each run: how many SNPs are seen at 0,...100 percent, using inverse binomial stats\n"
	    "//       Use minCover and selectFile to restrict the report to a zone or a cover threshold\n"
	    "//   -db_prevalenceTable   -db ACEDB -project project [-select selectFile] [-minCover c] [-Reference_genome] [-dropMonomodal] \n"
	    "//       Report a table, one GENE per line, one run per column, and several histograms\n"
	    "//       Use minCover and selectFile to restrict the report to a zone or a cover threshold\n"
	    "//       Use -Reference_genome to qualify the g.name in the exported table\n"
	    "//   -db_vcfTable -db ACEDB -project p -run_list list -o f [-gzo] \n"
	    "//       Export variants in VCF format\n"
	    "//         selects variant with tag project $P\n"
	    "//                                  nsCount $run in run_list\n"
	    "//         writes in file f.vcf[.gz]\n"
	    "//   -db_report           -db ACEDB : report general statistics, obsolete\n"
	    "//   -db_intersect        -db ACEDB : analyze how many SNPs are seen in how many samples\n"
	    "//   -db_pheno            -db ACEDB : get the phenotype of the variants\n"
	    "//        defined as the average of the absolute value of the deviation from the median of the phenotype of each run\n"
	    "//        weighted by the observed prevalence p of the SNP in that run (p = 0 for a wild type, p = 100 for a homozygous SNP)\n"
	    "//        this gives an indication of which SNP is correlated to which phenotype\n"
	    "//   -mergeExportPopulationTable : export a large table of the genotypes one column per run\n"
	    "//      The input should be a concatenation of .snp files exported in phase 2\n"

	    "//\n"
	    "// Other actions:\n"
	    "//   -ooFrequency 6 : use 2,4 or 6\n"
	    "//      Report the frequency of SOLiD transition errors at the center of each possible 6-mers\n"
	    "//\n"
	    "//   SNP confirmation and phasing:  \n"
	    "//   editSequence and countEdited are useful to confirm SNP counts, but do not work well for paired-end sequencing\n"
	    "//   they will at some point be superseded by aliExtend\n"
	    "//   -editSequence delta [-newIntrons newIntrons_file_name]: requires -fasta and -snp_list\n"
	    "//      Export a fasta file of edited segment: the SNP with delta bases on each side\n"
	    "//      plus fused fragments covering the candidate rearrangments listed in the newIntron file\n"
	    "//   -countEdited [-minAliPerCent <int>] [-runQuality f] :\n"
	    "//      count and analyze the hits on the fasta file exported by -editSequence\n"
	    "//      optionally reject reads not aligning on a sufficient fraction of their clipped lenght\n" 
	    "//   -aliExtend -snp_list filename [-target_class class]\n"
	    "//      draft program attempting to reanalyze the hits files given on stdin trying to locally remap the same reads\n"
	    "//      on the locally combinatorial variations of the genome implied by the snp_list i.e. 3:45321_A>C 3:45667_InsA \n"

	    //\n"
	    "// Filters\n"
	    "//   -select targetFile : only report SNPs (re)mapping into these targets\n"
	    "//      The first column selects the target, optionally coordinates in columns 2 and 3 restrict to a zone\n"
	    "//      This zone selection is required in -ventilate and optional in the phase 3 reports\n"
	    "//   -mito : select just the mitochondria\n"
	    "//   -help : export this on line help (-h or --help are also recognized)\n"
	    "// Caveat:\n"
	    "//   Lines starting with '#' in the input files are considered as comments and dropped out\n"
	    ) ;

	    /* obsolete	    "//   -qual : report the solexa quality of each base, default drop it\n" */


  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  SNP snp ;
  AC_HANDLE h = 0 ;
  const char *ccp, *dbName = 0 ;
 
  freeinit () ; 
  messErrorInit (argv[0]) ;

  if (0)
    {
    int i=1, line=1 ;
    ACEIN ai, aj = aceInCreate (argv[1], 0, h) ;
    aceInSpecial (aj, "\n") ;
    aceInCard (aj) ;
    ccp = aceInPos (aj) ;
    ai = aceInCreateFromText (ccp, 0, h) ; 
    aceInSpecial (ai, "\n") ;
    aceInCard (ai) ;
    ccp = aceInWordCut (aj, "\t", 0) ; 
    fprintf(stderr,"..%d.%d...%s\n",line,i++, ccp) ;
     ccp = aceInWordCut (aj, "\t", 0) ; 
    fprintf(stderr,"..%d.%d...%s\n",line,i++, ccp) ;
     ccp = aceInWordCut (aj, "\t", 0) ; 
    fprintf(stderr,"..%d.%d...%s\n",line,i++, ccp) ;
 
    i=11 ;
   ccp = aceInWordCut (ai, "\t", 0) ; 
    fprintf(stderr,"..%d.%d...%s\n",line,i++, ccp) ;
    ccp = aceInWordCut (ai, "\t", 0) ;
    fprintf(stderr,"..%d.%d...%s\n",line,i++, ccp) ;
    ccp = aceInWordCut (ai, "\t", 0) ;
    fprintf(stderr,"..%d.%d...%s\n",line,i++, ccp) ;
    ccp = aceInWordCut (ai, "\t", 0) ;
    fprintf(stderr,"..%d.%d...%s\n",line,i++, ccp) ;
    exit (1) ;
  }

  h = ac_new_handle () ;
  memset (&snp, 0, sizeof (SNP)) ;
  snp.h = h ;

  snp.runDict = dictHandleCreate (100, h) ;

  /* optional arguments */

  if (getCmdLineBool (&argc, argv, "-h") ||
      getCmdLineBool (&argc, argv, "-help") ||
      getCmdLineBool (&argc, argv, "--help")
      )
    usage (0) ;

  snp.smooth = getCmdLineBool (&argc, argv, "-smooth") ;
  snp.smooth = TRUE ;

  getCmdLineOption (&argc, argv, "-db", &dbName) ;

  snp.BRS_merge = getCmdLineBool (&argc, argv, "-BRS_merge") ;

  snp.minCover = 10 ;  snp.minMutant = 4 ;  snp.minFrequency = 20 ; /* defaults for BRS2snp */

  snp.mergeExportPopulationTable = getCmdLineBool (&argc, argv, "-mergeExportPopulationTable") ;
  snp.ventilate = getCmdLineBool (&argc, argv, "-ventilate") ;
  getCmdLineInt (&argc, argv, "-minCover", &snp.minCover) ;
  getCmdLineInt (&argc, argv, "-minMutant", &snp.minMutant) ;
  getCmdLineInt (&argc, argv, "-minFrequency", &snp.minFrequency) ;
  getCmdLineInt (&argc, argv, "-minAli", &snp.minAli) ;
  getCmdLineInt (&argc, argv, "-minAliPerCent", &snp.minAliPerCent) ;
  getCmdLineOption (&argc, argv, "-i", &snp.inFileName) ;
  getCmdLineOption (&argc, argv, "-o", &snp.outFileName) ;
  getCmdLineOption (&argc, argv, "-select", &snp.selectFileName) ;
  getCmdLineOption (&argc, argv, "-title", &snp.title) ;
  getCmdLineOption (&argc, argv, "-run", &snp.run) ;
  getCmdLineOption (&argc, argv, "-project", &snp.project) ;
  getCmdLineOption (&argc, argv, "-Reference_genome", &snp.Reference_genome) ;
  if (snp.Reference_genome)
    {
      char *cp = strchr (snp.Reference_genome, '.') ; if (cp) *cp = 0 ; 
    }
  getCmdLineOption (&argc, argv, "-run_list", &snp.runListFileName) ;
  getCmdLineOption (&argc, argv, "-snp_list", &snp.snpListFileName) ;
  getCmdLineOption (&argc, argv, "-newIntrons", &snp.newIntronsFileName) ;
  snp.gzi = getCmdLineBool (&argc, argv, "-gzi") ;
  snp.gzo = getCmdLineBool (&argc, argv, "-gzo") ;
  snp.mito = getCmdLineBool (&argc, argv, "-mito") ;
  snp.dropMultiplicity = getCmdLineBool (&argc, argv, "-dropMultiplicity") ;
  snp.unique = getCmdLineBool (&argc, argv, "-unique") ;
  snp.differential = getCmdLineBool (&argc, argv, "-differential") ;
  snp.dropMonomodal = getCmdLineBool (&argc, argv, "-dropMonomodal") ;
  snp.solid = getCmdLineBool (&argc, argv, "-solid") ;
  snp.dbCount = getCmdLineBool (&argc, argv, "-db_count") ;
  snp.vcfTable = getCmdLineBool (&argc, argv, "-db_vcfTable") ;
  getCmdLineOption (&argc, argv, "-db_val", &snp.snpValFileName) ;
  snp.frequencyTable = getCmdLineBool (&argc, argv, "-db_frequencyTable") ;
  snp.prevalenceTable = getCmdLineBool (&argc, argv, "-db_prevalenceTable") ;
  snp.frequencyHisto = getCmdLineBool (&argc, argv, "-db_frequencyHisto") ;
  snp.qual = getCmdLineBool (&argc, argv, "-qual") ;
  snp.hits2BRS = getCmdLineBool (&argc, argv, "-hits2BRS") ;
  getCmdLineInt (&argc, argv, "-Remove_inserts_shorter_than", &(snp.Remove_inserts_shorter_than)) ;

  snp.aliExtend = getCmdLineBool (&argc, argv, "-aliExtend") ;
  snp.phasing = getCmdLineBool (&argc, argv, "-phasing") ;
  getCmdLineInt (&argc, argv, "-ooFrequency", &(snp.ooFrequency)) ;
  getCmdLineInt (&argc, argv, "-editSequence", &(snp.editSequence)) ;
  getCmdLineInt (&argc, argv, "-countEdited", &(snp.analyzeEditedHits)) ;

  snp.db_report = getCmdLineBool (&argc, argv, "-db_report") ;

  snp.BRS_detect = getCmdLineBool (&argc, argv, "-BRS_detect") ;
  snp.BRS_Detect = getCmdLineBool (&argc, argv, "-BRS_Detect") ;
  snp.BRS_count = getCmdLineBool (&argc, argv, "-BRS_count") ;
  /* .keepSample = getCmdLineBool (&argc, argv, "-keepSample") ; */
  getCmdLineOption (&argc, argv, "-fasta", &(snp.fastaFileName)) ;
  getCmdLineOption (&argc, argv, "-targetGene", &snp.targetGeneFileName) ;
  getCmdLineOption (&argc, argv, "-selected8kbList", &snp.selected8kbFileName) ;

  getCmdLineOption (&argc, argv, "-db_remap2genome", &snp.remapFileName1) ;
  getCmdLineOption (&argc, argv, "-db_remap2genes", &snp.remapFileName2) ;
  snp.db_VCF = getCmdLineBool (&argc, argv, "-db_VCF") ;
  snp.db_intersect = getCmdLineBool (&argc, argv, "-db_intersect") ;
  snp.db_pheno = getCmdLineBool (&argc, argv, "-db_pheno") ;

  getCmdLineOption (&argc, argv, "-target_class", &snp.target_class) ;
  getCmdLineOption (&argc, argv, "-runQuality", &snp.runQualityFileName) ;
  snp.pool = getCmdLineBool (&argc, argv, "-pool") ;

 if (dbName)
    {
      const char *errors ;
      snp.db = ac_open_db (dbName, &errors);
      if (! snp.db)
	messcrash ("Failed to open db %s, error %s", dbName, errors) ;
    }
  snp.strategy = STRATEGY_RNA_SEQ ; /* default */
  if (getCmdLineOption (&argc, argv, "-strategy", &ccp))
    {
      if (! strcasecmp (ccp, "Exome") ||
	  ! strcasecmp (ccp, "Genome")
	  )
	{
	  snp.strategy = STRATEGY_EXOME ;
	  snp.strategy = STRATEGY_RNA_SEQ ; /* default */
	}
      else if (! strcasecmp (ccp, "RNA_seq"))
	{
	  snp.strategy = STRATEGY_RNA_SEQ ;
	}
      else
       {
	 fprintf (stderr, "-strategy %s , should be Exome, Genome or RNA_seq, sorry", ccp) ;
	 exit (1) ;
       }
    }

  fprintf (stderr, "// start: %s\n", timeShowNow()) ;
  if (argc > 1)
    usage (messprintf("sorry unknown arguments %s", argv[1])) ;

  cleanSnpName(0,0,0); /* for compiler happiness */
  
  /*   snpOpenInputFile (&snp) ; */
  snp.ai = aceInCreate (snp.inFileName, snp.gzi, snp.h) ;
  aceInSpecial (snp.ai,"\t\n") ;

  if (snp.selectFileName)
    snpPrepareZoneSelection (&snp) ;
 
  snp.runDict = dictHandleCreate (64, snp.h) ; 
  snp.targetDict = dictHandleCreate (10000, snp.h) ;
  snp.target_classDict = dictHandleCreate (4, snp.h) ;
  snp.chromDict = dictHandleCreate (100, snp.h) ;
  snp.eeDict = eeDictCreate (&snp, snp.h) ;  /* a,t,g,c and all error types in canonical order */
  snp.aaa = arrayHandleCreate (1024, BigArray, snp.h) ; /* Array of all aa in target */
  snp.aaaP = arrayHandleCreate (1024, BigArray, snp.h) ; 
  snp.errCounts = arrayHandleCreate (1024, long int, snp.h) ;


  if (snp.selected8kbFileName)
    {
      snp.selected8kbDict = dictHandleCreate (100, snp.h) ;
      snp.selected8kbBitSet = bitSetCreate (300000, snp.h) ;
      snpParseSelected8kbTargetFile (&snp) ;
    }
  if (snp.runQualityFileName)
    snpParseRunQualityFile (&snp) ;

  if (snp.targetGeneFileName)
    snpParseTargetGeneFile (&snp) ;
  if (snp.fastaFileName)
    snpParseFastaFile (&snp) ;
  if (snp.remapFileName1 || snp.remapFileName2)
    snpCreateAtlas (&snp) ;

  if (snp.hits2BRS || snp.ventilate)
    {
      if (! snp.run)
	messcrash ("-hits2BRS/ventilate options require -run run_name") ;
      if (snp.remapFileName1 || snp.remapFileName2)
	messcrash ("It is too complex to -count and -remap at the same time") ;
      snp.ao = aceOutCreate (snp.outFileName, 0, snp.gzo, snp.h) ;
      snpHits2BRS (&snp) ;
    }
  else if (snp.BRS_merge)
    snpBRSMerge (&snp) ;
  else if (snp.remapFileName1)
    { snpRemap1 (&snp) ;  snpRemap2 (&snp) ; }
  else if (snp.remapFileName2)
    snpRemap2 (&snp) ;
  else if (snp.BRS_detect)
    {
      if (snp.snpListFileName)
	{
	  fprintf (stderr, "FATAL ERROR; Wrong parameters, option -BRS_detect detects all SNPs passing the thresholds, it is incompatible with -snp_list\nTry snp --help") ;
	  exit (1) ;
	}
      snpBRS2tsf (&snp) ;
    }
  else if (snp.BRS_count)
    {
      if (snp.snpListFileName)
	snpParseSnpList (&snp) ;
      else
	{
	  fprintf (stderr, "FATAL ERROR; Wrong parameters, option -BRS_count requires option -snp_list\nTry snp --help") ;
	  exit (1) ; 
	}
      snpBRS2tsf (&snp) ;
    }
  else if (snp.db_pheno)
    {
      snpPhenotype (&snp) ;
    }
  else if (snp.db_report)
    {
      snpDeepTableReport (&snp) ;
    }
  else if (snp.ooFrequency)
    {
      if (snp.ooFrequency != 2 && snp.ooFrequency != 4 && snp.ooFrequency != 6)
	messcrash ("The ooFrequency oligo length should be 2, 4 or 6, sorry") ;
      snpOoFrequency (&snp) ;
    }
  else if (snp.editSequence)
    {
      if (snp.editSequence < 0)
	messcrash ("-editSequence should be a positive number = the zone to export on each size of the SNP\n") ;
      if (! snp.snpListFileName)
	messcrash ("-editSequence requires -snp_list and will export the DNA around these SNPs, sorry") ;
      if (! snp.fastaFileName)
	messcrash ("-editSequence requires -fasta parameter pointing to the reference sequence, sorry") ;
      snpExportEditedSequence (&snp) ;
    }
  else if (snp.phasing)
    { 
      if (! snp.run)
	messcrash ("-phasing requires -run run_name, sorry") ;

      if (! snp.snpListFileName)
	messcrash ("-phasing requires -snp_list and will flag SNPs validated by phasing") ;
       snpPhasing (&snp) ;
    }
  else if (snp.aliExtend)
    {
      if (! snp.run)
	messcrash ("-aliExtend requires -run run_name, sorry") ;
      if (! snp.fastaFileName)
	messcrash ("-aliExtend requires -fasta parameter pointing to the reference sequence, sorry") ;
      if (! snp.snpListFileName)
	messcrash ("-aliExtend requires -snp_list and will locally realign on these variations of the same target, sorry") ;
      snpParseFastaFile (&snp) ;
      snpAliExtendParseSnpList (&snp, 0) ;  /* fills array: snp->aliSnps */
      snpAliExtend (&snp) ;
    }
  else if (snp.db_intersect)
     snpIntersect (&snp) ;
  else if (snp.mergeExportPopulationTable)
    snpMergeExportPopulationTable (&snp) ;
  else if (snp.frequencyTable)
    {
      if (! snp.db) 
	usage ("-db_frequencyTable requires -db MetaDB") ;
      if (! snp.project) 
	usage ("-db_frequencyTable requires -project projectName") ;
      if (! snp.inFileName)
	usage ("-db_frequencyTable requires -i snp_file_name, the file is parsed twice hence stdin cannot be used") ;
      snpFrequencyTable (&snp, FALSE) ;
      /*
      ac_free (snp.ai) ;
      snp.ai = aceInCreate (snp.inFileName, snp.gzi, snp.h) ;
      aceInSpecial (snp.ai,"\t\n") ; 
      snpFrequencyTable (&snp, FALSE) ;
      */
    }
  else if (snp.vcfTable)
    {
      if (! snp.db) 
	usage ("-db_vcfTable requires -db MetaDB") ;
      if (! snp.runListFileName)
	messcrash ("-db_vcfTable requires -run_list, the list of runsto be exported") ;
      if (0 && ! snp.run) 
	usage ("-db_vcfTable requires -run runName") ;
      if (! snp.project) 
	usage ("-db_frequencyTable requires -project projectName") ;
      snpVcfTable (&snp) ;
    }
   else if (snp.prevalenceTable)
    {
      if (! snp.db) 
	usage ("-db_prevalenceTable requires -db MetaDB") ;
      if (! snp.project) 
	usage ("-db_prevalenceTable requires -project projectName") ;
      snpPrevalenceTable (&snp) ;
    }
  else if (snp.frequencyHisto)
    {
      if (! snp.db) 
	usage ("-db_frequencyHisto requires -db MetaDB") ;
      if (! snp.project) 
	usage ("-db_frequencyHisto requires -project projectName") ;
      snpFrequencyTable (&snp, FALSE) ;
    }
   else if (snp.dbCount)
    {
      if (! snp.db) 
	usage ("-db_count requires -db MetaDB") ;
      if (! snp.project) 
	usage ("-db_count requires -project projectName") ;
      if (! snp.inFileName && ! snp.snpValFileName)
	usage ("-db_count requires -i snp_file_name OR -db_val snp2.val.txt file list") ;
      snpFrequencyTable (&snp, FALSE) ;
    }
 else if (snp.analyzeEditedHits)
    {
      if (! snp.run)
	messcrash ("-editSequence requires -run run_name, sorry") ;
      if (snp.analyzeEditedHits < 10)
	messcrash ("-editSequence value (now %d) should be at least 10, it should be the middle of the analyzed reads"
		   , snp.analyzeEditedHits) ;
      snpAnalyzeEditedHits (&snp) ;
    }
  fprintf (stderr, "// done: %s\n", timeShowNow()) ;
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   fprintf (stderr, "// max memory %d Mb\n", mx) ;
  }
  ac_free (snp.db) ;
  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

