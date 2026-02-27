/*
 * sa.tabular.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created February 26, 2026

 * This code is public.

 * This code handles the MagicBlast tabular format
 * Given an alignment in internal format (ALI structures and _ERR) it exports a tabular file
 */

#include "sa.h"

/*********************************************************************/
/*********************************************************************/

/* export  Tabular
 * given an Array with all hits corresponding to a fragment

 * The tabular format has 25 columns
 */

void saExportTabular (const PP *pp, BB *bb, Array aa)
{
  int ii, iMax = arrayMax (aa) ;
  ALIGN *up = iMax ? arrp (aa, 0, ALIGN) : 0 ;
  ACEOUT ao = 0 ;
  
  for (ii = 0 ; ii < iMax ; ii++, up++)
    {
      /* 1: read identifier */
      aceOut (ao, dictName (bb->dict, up->read >> 1)) ; 
      /* 2: target identifier */
      aceOut (ao, "\t") ;
      aceOut (ao, dictName (pp->bbG.dict, up->chrom >> 1)) ; 
      /* 3: percent identity of the alignment */
      float z = 100 ;
      aceOut (ao, "\t") ;
      aceOutPercent (ao, z) ;
      /* 4 5 6 not used */
      /* 7, 8 ali start/stop on read */
      aceOut (ao, "\t\t\t\t") ;
      aceOutInt (ao, up->x1) ;
      aceOut (ao, "\t") ;
      aceOutInt (ao, up->x2) ;
      /* 9, 10 ali start/stop on target */
      aceOut (ao, "\t") ;
      aceOutInt (ao, up->a1) ;
      aceOut (ao, "\t") ;
      aceOutInt (ao, up->a2) ;
      /* 11, 12 not used */
      /* 13 alignment score */
      aceOut (ao, "\t\t\t") ;
      aceOutInt (ao, up->chainScore) ;
      /* 14/15  read/chrom strands */
      aceOut (ao, up->x1 < up->x2 ? "\t+" : "\t-") ;
      aceOut (ao, up->a1 < up->a2 ? "\t+" : "\t-") ;
      /* 16 read length */
      aceOut (ao, "\t") ;
      aceOutInt (ao, up->x2 - up->x1 + 1) ;
      /* 17 BTOP string */
      aceOut (ao, "\t") ;
      aceOutInt (ao, up->x2 - up->x1 + 1) ;
      /* 18 target multiplicity */
      aceOut (ao, "\t") ;
      aceOutInt (ao, 1) ; /* nb of exported chains */
      /* 19 not used */
      /* 20 chain */
      aceOut (ao, "\t\t") ;
      aceOutInt (ao, up->chain) ; /* nb of exported chains */
      /* 21 complement of left overhang */
      aceOut (ao, "\t") ;
      aceOut (ao, up->leftOverhang ? dictName (bb->dict, up->leftOverhang) : "-") ;
      /* 22 right overhang */
      aceOut (ao, "\t") ;
      aceOut (ao, up->rightOverhang ? dictName (bb->dict, up->rightOverhang) : "-") ; 
      /* 23 chrom of other read of pair */ 
      aceOut (ao, "\t") ;
      aceOut (ao, "-") ; /* right */
      /* 24 start pos of mate alignment (negative if divergent)  */
      aceOut (ao, "\t") ;
      aceOutInt (ao, 0) ; /* right */
      /* 25 chain score */
      aceOut (ao, "\t") ;
      aceOutInt (ao, up->chainScore) ; /* right */
      /*********************DONE*******************/
    }
      
  return ;
}
