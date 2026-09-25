#ifndef SCRATCH_DEF
#define SCRATCH_DEF
#include "ac.h"

typedef struct scrStruct *SCR ;
SCR   scratchCreate (const char *dirName, AC_HANDLE h) ;
int   scratchArrayCount (SCR scr) ;
int   scratchArrayPut (SCR scr, Array aa) ;

Array scratchArrayDoGet (SCR scr, int k, Array aa, int typeSize, AC_HANDLE h) ;
#define scratchArrayGet(scr,k,aa,TYPE,h) scratchArrayDoGet((scr),(k),aa,sizeof(TYPE),(h))

#endif
