/*  Last edited: Apr 25 15:14 1995 (srk) */

/* $Id: disk_.h,v 1.2 2017/06/06 20:22:00 mieg Exp $ */

 
             /*********************************************/
             /* disk_.h                                   */
             /* disk block universal structure            */
             /* replicated in all classes                 */
             /*********************************************/
 
#ifndef DEF_DISK__H 
#define DEF_DISK__H

#ifndef DEF_DISK
#include "disk.h"
#endif
/* lscpu   on Linux gives the cache size of the hardware */
/* sysctl  on Mac */
#define BLOC_SIZE 1024      /*number of char in a block*/
extern int _DISQUE_BLOC_SIZE ;
typedef struct  blockheader
  { DISK disk,                 /* where */
         nextdisk ;
    int session ;              /* when */
    KEY key ;                  /* who */
    char type ;                /* how */
  } BLOCKHEADER ;


#endif


