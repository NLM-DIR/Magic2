/*
 * sortalign : RNA aligner

 * Created April 18, 2025
 * In collaboration with Greg Boratyn, NCBI

 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality

 * The algorithm extract words from the union of the targets and from the reads creating 2 long tables
 * Both tables are sorted
 *  sorting the targets is a preprocessing, so the time does not count
 *  sorting the reads is done datablock by datablock
 * The two sorted tables are 'merged' which is local in memory
 *  exporting a sequential table.
 * That table is sorted per read
 *  alignments are extended using the magicblast jumper algorithm
 * Finally auxiliary tables, like introns, or wiggles, are exported on demand

 * The code is heavilly parallelized using the channel paradigm borrowed from the Go language
 */


/*
  #define ARRAY_CHECK
  #define MALLOC_CHECK
*/
/* Examples:  sortalign -h */

/* edit the version only if you edit the way the index is constructed */ 
#define INDEXVERSION "2026_02_7"

/***********************************************************************************/

#ifdef USEGPUzzzz
  #define NAGENTS 12
  #define NBLOCKS 12
#else
  #define YANN 
  #define NAGENTS 40
  #define NBLOCKS 40
#endif

/***********************************************************************************/

#include "ac.h"
#include "channel.h"
#include <dirent.h>
#include <zlib.h>
#include <stdatomic.h>
#include "../wsra/sra_read.h"
#include "sa.common.h"



typedef struct nodeStruct { double x ; CHAN *cx, *cy, *cu, *cv, *done ; int k ; } NODE ;

typedef enum {FASTA=1, FASTQ, FASTC, RAW, SRA, SRACACHE, SRACACHE1, SRACACHE2, INTRONS, GFF} DnaFormat ;
typedef struct targetClassStruct {
  char targetClass ; /* single char a-z, A-Z */
  int priority ;
  DnaFormat format ;
  const char *fileName ;
  Array cws ;
} TC ;
		  
typedef struct runClassStruct {
  int run ;   /* index in pp->runDict */
  BOOL pairedEnd ;
  DnaFormat format ;
  int jump5r1, jump5r2 ; /* jump bases at the 5prime end of read1 and read2 */
  const char *fileName1 ;
  const char *fileName2 ;
  atomic_int lane ;
  atomic_int laneDone ;
  ACEOUT aoSam ;
} RC ;

typedef struct aStruct {
  char a1L[32];
  char a1R[32];
  char a2L[32];
  char a2R[32];
} ADAPTORS ;


#define SLMAX 1000
#define LETTERMAX 1000
#define OVERHANGMAX 600   /* 4*30*5  overhang: 32 bases max, starting on base atgc, 5 counts atgcn per position */

typedef struct PSDstruct {
  int run ;   /* index in p.runDict */
  int nFiles ;
  long int nBase1, nBase2 ;
  long int nPairs, nReads ;
  Array lengthDistribution ;
  int minReadLength, maxReadLength ;
  long int letterProfile1[5 * LETTERMAX] ;
  long int letterProfile2[5 * LETTERMAX] ;
  long int ATGCN[5] ;

} PSD ;

typedef struct runStatStruct {
  PSD p ;
  long int lowEntropy, tooShort ;
  long int lowEntropyBases, tooShortBases ;
  long int nCompatiblePairs, nIncompatiblePairs, n2ChromsPairs, nOrphans, nCirclePairs ;
  long int nPairsAligned, nBaseAligned1, nBaseAligned2 ;
  long int cds, utr, intronic, intergenic ;
  long int nMultiAligned[11] ;
  long int nReadsAlignedPerTargetClass[256] ;
  long int nBasesAlignedPerTargetClass[256] ;
  long int nPerfectReads ;
  long int nAlignments ;
  long int nClippedPolyA ;
  long int nClippedPolyT ;
  long int nClippedSL ;
  long int nClippedAdaptor1L ;
  long int nClippedAdaptor2L ;
  long int nClippedAdaptor1R ;
  long int nClippedAdaptor2R ;
  long int nTello[12] ;
  long int nTelloRead ;
  long int SlRead ;
  long int nClippedSls[SLMAX] ;
  long int nSupportedIntrons ;
  long int nIntronSupportPlus ;
  long int nIntronSupportMinus ;
  int gt_ag_Support ;
  int ct_ac_Support ;
  float intronStranding ;
  long int wiggleCumul ; /* in million bases */
  long int wiggleLCumul ; /* in million bases */
  long int wiggleRCumul ; /* in million bases */
  Array insertLengthDistribution ;
  ADAPTORS adaptors ;
  long int nN, nErr, nMID ;
  int accessibleLength5kb ;
  int accessibleLength8kb ;
  long int overhangL1 [OVERHANGMAX] ; /* read 1 left overhang  */
  long int overhangL2 [OVERHANGMAX] ; /* read 2 left overhang */
  
  long int overhangR1 [OVERHANGMAX] ; /* read 1 right overhang  */
  long int overhangR2 [OVERHANGMAX] ; /* read 2 right overhang */
  
  int GF[256], GR[256] ; /* number of reads aligned per target_class on Forward/Reverse strand */
  Array errors ;  /* substitutions, insertions, deletions counts */
  /* coverage of long transcripts ? */
} RunSTAT ;
		  
typedef struct bStruct {
  AC_HANDLE h ;
  int run, lane, readerAgent ;
  RC rc ;
  DICT *dict, *errDict ;
  mytime_t start, stop ;
  int gpu ;
  BigArray dnaCoords ;   /* offSets of the dna in the globalDna array */
  Array dnas ;           /* Array of const char Arrays */
  Array dnasR ;          /* Their reverse complement, only computed for the genome */
  Array quals ;
  BigArray globalDna ;   /* concatenation of all sequence DNAs separated by blocks of nnn */
  BigArray globalDnaR ;  /* concatenation of all reverse sequences, in the same order */
  BigArray hits ;        /* BigArray of read<->genome hits */
  BigArray *cwsN ;         /* BigArray of codeWords */
  long unsigned int nSeqs ;  /* number of sequences in bloc */
  long unsigned int length ; /* cumulated number of bases */
  long unsigned int nerr ;   /* cumulated number of errors */
  long unsigned int nAli ;   /* cumulated number of alignments */
  long unsigned int aliDx ;  /* cumulated aligned read length */
  long unsigned int aliDa ;  /* cumulated genome coverage */

  char *gzBuffer ;
  
  /*   BitSet isAligned ; */
  BigArray aligns ; /* final alignments */  
  RunSTAT runStat ;
  Array cpuStats ;
  Array errors ;
  Array wiggles ;
  Array wigglesL ;
  Array wigglesR ;
  Array wigglesP ;
  Array wigglesNU ;
  Array confirmedPolyAs ;
  Array confirmedSLs ;
  Array confirmedIntrons ;
  Array doubleIntrons ;
  BOOL isGenome ;
  int isRna ; /* 2: user defined RNA, -2: user defined DNA, 1: autodefined RNA, -1 autodefined DNA */

  int step, skips0, skips1, skips2, skips3, skips4, skipsFound, skipsNotFound ;
  long int genomeLength ;
  long int nSLsSupport ;
  long int nPolyASupport ;
  long int nIntronSupportPlus ;
  long int nIntronSupportMinus ;    
  vTXT txt1, txt2 ; /* a pair of reusable txt buffer */ 
} BB ;  

typedef struct pStruct {
  AC_HANDLE h ;
  AC_HANDLE bamHandle ;
  BOOL debug, gzi, gzo ;
  BOOL createIndex ;
  BOOL noJump ;
  BOOL align ;
  BOOL justStats ;
  BOOL wiggle ;
  BOOL wiggleEnds ;
  BOOL introns ;
  BOOL snps ;
  BOOL ignoreIntronSeeds ;
  BOOL isWorm ;
  const char *runName ;
  const char *inFileName ;
  const char *inConfigFileName ;
  const char *outFileName ;
  const char *indexName ;
  const char *tFileName ;
  const char *tConfigFileName ;
  const char *tFileBinaryCwsName ;
  const char *tFileBinaryDnaName ;
  const char *tFileBinaryDnaRName ;
  const char *tFileBinaryIdsName ;
  const char *tFileBinaryCoordsName ;

  Array SLs ;
  const char *rawAdaptor1R, *rawAdaptor2R ;
  ADAPTORS adaptors ;
  /* Agents:
     R read parser
     G genome parser
     L load regulator,
     C code words,
     S sort words,
     M match words
     O order hits
     A align in detail
     E export ali
  */
  CHAN *fpChan ; /* RC, Run Config -> readParser */
  CHAN *npChan ; /* int, number of BB emitted for each inFile */
  CHAN *gmChan ; /* Genome is ready, signals the Matcher */
  CHAN *plChan ; /* Parser sends a BB to the load regulator */
  CHAN *lcChan ; /* Load regulator sends a BB to the word Coder */
  CHAN *csChan ; /* Coder sends words to the Sorter */
  CHAN *smChan ; /* Sorter sends words to the Matcher */
  CHAN *moChan ; /* Matcher needs seeds to be Ordered */
  CHAN *oaChan ; /* Ordered matches  sent to the Aligner */
  CHAN *aeChan ; /* Aligner results to be Exported */
  CHAN *doneChan ; /* return to main program */
  CHAN *wwChan ; /* Wiggler */
  CHAN *wwDoneChan ; /* Wiggler */


  BB bbG ;  /* genome or genes target */
  Array tArray ;
  Array runStats ;
  Array runLanes ;
  Array runLanesDone ;
  Array geneBoxes ;  /* gene coordinates per chromosome */
  Array geneCounts ; /* expression counts per chromosomes */
  BigArray knownIntrons ;
  Array wiggles ;
  Array wiggleCumuls ;
  Array wigglesL ;
  Array wiggleLCumuls ;
  Array wigglesR ;
  Array wiggleRCumuls ;
  Array wigglesP ;
  Array wigglesNU ;
  Array cdss ;
  Array utrs ;
  Array intronics ;
  Array intergenics ;
  Array confirmedPolyAs ;
  Array confirmedSLs ;
  Array confirmedIntrons ;
  Array doubleIntrons ;
  BOOL fasta, fastq, fastc, raw, solid, sra, sraCaching ;
  BOOL sam, bam, hitsFormat, exportSamSequence, exportSamQuality, qualityFactors ;
  BOOL strand, antiStrand ;
  BOOL isDna, isRna ;
  int bonus[256] ;
  DICT *runDict ;
  DICT *targetClassDict ;
  DICT *geneDict ;
  DICT *mrnaDict ;
  DICT *wiggleFileDict ;
  const char *method ;
  int run ;
  int nFiles ;  /* number of input sequence files */
  int agent ;   /* instance of the agent */
  int nIndex ;  /* number of parallel target indexes */
  int nBlocks ;
  int blocMaxBases ; /* max number of bases read as one bloc */
  int iStep ;        /* default 2, take a read seed every iStep */
  int tStep ;        /* default 3, take a target seed every tStep */
  int maxTargetRepeats ;
  int tMaxTargetRepeats ;
  int seedLength ;
  int maxIntron ;
  int BMAX ; /* max number of bases in a block, default 200M */
  int errCost ;
  int errMax ;       /* (--align case) max number of errors in seed extension */
  int minScore, minAli, minAliPerCent ;
  int minLength, minEntropy ;
  int errRateMax ;       /* (--align case) max number of errors in seed extension */
  int OVLN ;
  int gpu ;
  int maxSraGb ; /* max number of Gigabases in each SRA download, 0 : no max */
  BOOL sraOutFormatPE ; /* default: 4 lines per pair (>id1, atgc, >id2. atgc */
  BOOL deduplicate ;
  BOOL sraOutFormatPEQ ; /* 8 lines per pair (twice fastq) */
  BOOL sraOutFormatFasta ; /* 2 files per pair, fasta format */
  BOOL sraOutFormatFastq ; /* 2 files per pair, fastq format */
  long int wiggleCumul ; /* in million bases */
  long int cds, utr, intronic, intergenic ;
  BOOL splice ;
  long int nRawReads, nRawBases, genomeLength ;
  int wiggle_step ;
  float *runStranding ;
} PP ;

#ifdef JUNKINCOMMON

#define NSHIFTEDTARGETREPEATBITS 8

typedef struct codeWordsStruct {
  unsigned int seed ; /* 32 bits = 16 bases, 2 bits per base */
  int nam ; /* index in readDict or chromDict << 1 | (0x1 for minus words) */
  int pos ;  /* bio coordinate of first letter of seed */
  unsigned int intron ;
} __attribute__((aligned(16))) CW ;

typedef struct hitStruct {
  unsigned int read ;  /* index in readDict */
  unsigned int chrom ; /* index in chromDict << 1 | (0x1 if minus strand) */
  unsigned int a1 ;  /* bio coordinates on chrom (base 1) */
  unsigned int x1 ;  /* bio coordinate on read */
} __attribute__((aligned(16))) HIT ;
#endif
typedef struct countChromStruct {
  float weight ;
  int seeds, chrom ;
  int seed1, seed2, seed4, seed8, seed16, seed32 ;
  int i1, i2, x1, x2 ;
  unsigned int a1 ;  /* bio coordinates on chrom (base 1) */
  unsigned int a2 ;  /* bio coordinate on read */
}  __attribute__((aligned(32))) COUNTCHROM ;

typedef struct alignStruct {
  int read ;
  int targetClass ;
  int chrom ;
  int a0, x0 ;  /* coordinate of the hit before extension */
  int a1, a2 ;  /* bio position in chrom, a1 < a2 on minus strand */
  int x1, x2 ;  /* bio position in read, x1 < x2 always */
  int chain, chainX1, chainX2, chainA1, chainA2 ;
  int id, previous, next ;
  int ali, chainAli, score, chainScore ;
  int pairScore, mateChrom, mateA1, mateA2, pairLength ;
  int nN, nErr, nMID, chainErr, chainMID ;
  int nTargetRepeats ;
  int nChains ;
  int readLength ;
  int leftClip, rightClip ;
  int pA ;
  int errShort, errLong ; /* bb->dict */
  int leftOverhang, rightOverhang ; /* bb->dict */
  int donor, acceptor ;
  Array errors ;
} __attribute__((aligned(32))) ALIGN ;

typedef struct geneBoxStruct {
  int chrom ;
  int a1, a2 ; 
  int gene ; /* index in pp->geneDict */
  int mrna ;
  char flag ;
  char strand ;
  char friends ;
} __attribute__((aligned(32))) GBX ;
		  
typedef struct geneStruct {
  int chrom ;
  int a1, a2 ; 
  int gene ; /* index in pp->geneDict */
} __attribute__((aligned(16))) GENE ;
		  
typedef struct intronStruct {
  int run ;
  int mrna ;
  int chrom ;
  int n, nR, a1, a2 ;
  char feet[6] ;
} __attribute__((aligned(32))) INTRON ;
		  
typedef struct doubleIntronStruct {
  int run ;
  int chrom ;
  int n, nR, a1, a2, b1, b2 ;
  char feet1[6] ;
  char feet2[6] ;
} __attribute__((aligned(32))) DOUBLEINTRON ;

typedef struct polyAStruct {
  int chrom ;  /* chrom indicates the strand */
  int run ;
  int a1 ; /* first A in the orientation of the chrom */
  int n ;
} __attribute__((aligned(32))) POLYA ;

typedef struct SLsStruct {
  int chrom ;  /* chrom indicates the strand */
  int run ; /* run << 4 | SLnumber (0..14 always < 15) */
  int a1 ; /* first A in the orientation of the chrom */
  int n ;
} __attribute__((aligned(32))) SLS ;

typedef struct cpuStatStruct {
  char nam[32] ;
  int agent ;
  int nB ; /* blocs treated */
  long int n ;
  clock_t tA ; /* time time active */
} CpuSTAT ;


#define step1 256
#define step2 512
#define step3 1024
#define step4 4096
/* was 256 1024 4096 16384 */

#define NTARGETREPEATBITS 5


#define mstep1 255
#define mstep2 510
#define mstep3 1020
#define mstep4 4080
/* was 255 1020 4080 16320 */

#include <pthread.h>
#include <time.h>

typedef struct timespec TMS ;

/**************************************************************/

/* sa.main.c */
void saUsage (char *message, int argc, const char **argv) ; 
void saCpuStatRegister (const char *nam, int agent, Array cpuStats, clock_t t1, clock_t t2, long int n) ;
  
/* sa.config.c */
void saConfigIsIndexAccessible (PP *pp) ;
void saConfigIsOutDirWritable (PP *pp) ;
int saConfigCheckTargetIndex (PP *pp) ; 
Array saConfigGetRuns (PP *pp, Array runStats) ;    
/* saSetGetAdaptors:
 * set=0: get, 1: set, 2: hard set non rewritable
 * isRna: >=0: favor introns, search polyAs, export gene expression
 *        <0: none of the above, disfavor deletions.
 */
BOOL saSetGetAdaptors (int set, int *isRnap, ADAPTORS *aa, int run) ;

/* sa.gff.c */
long int saGffParser (PP *pp, TC *tc) ;
void saGffBinaryParser (PP *pp) ;
long int saIntronParser (PP *pp, TC *tc) ; 

/* sa.sort.c */
int saSort (BigArray aa, int type) ; /* return usedGPU */
/* sa.gpusort.c */


/* sa.sam.c */
int saSamExport (ACEOUT ao, ACEOUT aoe, const PP *pp, BB *bb) ;       
ACEOUT saSamCreateFile (const PP *pp, BB *bb, BOOL isError, AC_HANDLE h) ; 

/* sa.introns.c */
void saIntronsOptimize (BB *bb, ALIGN *vp, ALIGN *wp, Array dnaG)  ; 
void saPolyAsExport (PP *pp, Array aaa) ;
void saPolyAsCumulate (PP *pp, BB *bb) ;
void saSLsExport (PP *pp, Array aaa) ;
void saSLsCumulate (PP *pp, BB *bb) ;
void saIntronsExport (PP *pp, Array aaa) ;
void saDoubleIntronsExport (PP *pp, Array aaa) ;
void saIntronsCumulate (PP *pp, BB *bb) ;
void saDoubleIntronsCumulate (PP *pp, BB *bb) ;
int saSupportedIntrons (const PP *pp, int run) ;
void saIntronStranding (PP *pp, Array aa) ;

/* sa.targetIndex.c */
Array saTargetParseConfig (PP *pp) ;
void saTargetIndexCreate (PP *pp) ;
void saTargetIndexGenomeParser (const void *vp) ;
void saDictMapRead (DICT *dict, const char *fNam) ;
void saDictMapWrite (DICT *dict, const char *fNam) ;

/* sa.seeds.c */
int saCodeIntronSeeds (PP *pp, BB *bbG) ;
void saCodeSequenceSeeds (const PP *pp, BB *bb, int step, BOOL isTarget) ;
  
/* sa.sequenceParser.c */
void saSequenceParse (const PP *pp, RC *rc, TC *tc, BB *bb, int isGenome) ;
int saSequenceParseSraDownload (const PP *pp, const char *sraID) ;
void saSequenceParseGzBuffer (const PP *pp, BB *bb) ;
void saSequenceDeduplicate (const PP *pp, BB *bb) ;

/* sa.uringSequenceParser.c */
void saUringSequenceParser (const PP *pp, RC *rc, TC *tc, BB *bb) ;

/* sa.wiggle */
void saWiggleCumulate (const PP *pp, BB *bb) ;
void saWiggleExport (PP *pp, int nAgents) ;

/* sa.stats */
void saRunStatExport (const PP *pp, Array runStats) ;
void saCpuStatExport (const PP *pp, Array stats) ;
void saCpuStatCumulate (Array aa, Array a) ;
void saRunStatsCumulate (int run, PP *pp, BB *bb) ;
BOOL saReadAdaptors (ADAPTORS *adaptors, RunSTAT *up, BOOL coded) ;
  
/* sa.align */
void saAlign (const void *vp) ;
void saAlignDo (const PP *pp, BB *bb) ;
int saAlignOrder (const void *va, const void *vb) ;

/* sa.utils */
int saBestNumactlNode (void) ;

/* sa.tests */
void saCreateRandomGenome (PP *pp, int nMb) ;

/**************************************************************/
/**************************************************************/
/**************************************************************/
