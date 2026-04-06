/* sa.common.h
 * Declarations common to code for GPU and CPU
 *
 */

#define NSHIFTEDTARGETREPEATBITS 8

typedef struct codeWordsStruct {
  unsigned int seed ; /* 32 bits = 16 bases, 2 bits per base */
  unsigned int nam ; /* index in readDict or chromDict << 1 | (0x1 for minus words) */
  unsigned int pos ;  /* bio coordinate of first letter of seed */
  unsigned int intron ;
} __attribute__((aligned(16))) CW ;


typedef struct hitStruct {
  unsigned int read ;  /* index in readDict */
  unsigned int chrom ; /* index in chromDict << 1 | (0x1 if minus strand) */
  unsigned int a1 ;  /* bio coordinates on chrom (base 1) */
  unsigned int x1 ;  /* bio coordinate on read */
} __attribute__((aligned(16))) HIT ;

typedef struct seedMatchStruct {
  unsigned int readSeed ;    /* place holder, set to zero */
  unsigned int read ;        /* index in bb->dict << 1 | (0x1 for minus words) */
  unsigned int x1 ;          /* bio coordinate of first letter of seed in read */
  unsigned int readFlags ;   /* copied from read index */
  
  unsigned int targetSeed ;  /* place holder, set to zero */
  unsigned int target ;      /* index in pp->bbG.dict << 1 | (0x1 for minus chromosome strand) */
  unsigned int a1 ;          /* bio coordinate of first letter of seed in target */
  unsigned int targetFlags ; /* copied from target index */
} __attribute__((aligned(32))) SEEDMATCH ;

