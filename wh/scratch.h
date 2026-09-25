#ifndef SCRATCH_DEF
#define SCRATCH_DEF
#include "ac.h"

/* Transparent temporary disk storage of BigArrays
 *
 * The stored arrays are numbered sequencially
 * They must be of homogeneous TYPE (inforced weakly via aa->size)
 * The disk space is released by calling ac_free (scr) or ac_free (h)
 *
 * int k = scratchPut (scr, aa) : returns a chronological integer
 * int kMax = scratchCount (scr) : return the number of Put
 * scratchGet (scr, aa, k, TYPE, h)
 *    k > kMax is forbidden
 *    if (aa == 0) Get allocates a new array of TYPE on h 
 *    if (aa != 0) Get checks aa->size, reuses aa and adjust aa->max
 */

typedef struct scrStruct *SCR ;
SCR   scratchCreate (const char *dirName, AC_HANDLE h) ;
int   scratchCount (SCR scr) ;
int   scratchPut (SCR scr, BigArray aa) ;
BigArray uScratchGet (SCR scr, int k, BigArray aa, int typeSize, AC_HANDLE h) ;
#define scratchGet(scr,k,aa,TYPE,h) uScratchGet((scr),(k),aa,sizeof(TYPE),(h))

#endif
