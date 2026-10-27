/*  File: channel_.h
 *  channels interface : derived from the GO language channels
 *  Author: Danielle and Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov) Mark Sienkiewicz (sienkiew@stsci.edu)
 *  Copyright (C) J Thierry-Mieg, 2015
 * -------------------------------------------------------------------
 * Created: mars 2015, mieg
 *   2026: revised. channelSelect and the remote_channel extension were removed.
 *-------------------------------------------------------------------
 */

#ifndef CHANNEL__H
#define CHANNEL__H

#include "wego.h"
#include <pthread.h>

/*************************************************************************/
/* Private implementation :
 */

/* Please use these macros
 * macros are used because the allow to pass around C-types as parameters
 */
#define channelCreate(_cMax,_type,_h) uChannelCreate((_cMax),sizeof(_type),(_h))

#define channelGet(_chan,_vp,_type) (uChannelMultiGet((_chan),(_vp),sizeof(_type),1,TRUE))
#define channelGive(_chan,_vp,_type) (*(_type*)uChannelGet((_chan),(_vp),sizeof(_type)))
#define channelPut(_chan,_vp,_type) uChannelMultiPut((_chan),(_vp),sizeof(_type),1,TRUE)

#define channelMultiGet(_chan,_vp,_max,_type) (uChannelMultiGet((_chan),(_vp),sizeof(_type),(_max),TRUE))
#define channelMultiPut(_chan,_vp,_max,_type) uChannelMultiPut((_chan),(_vp),sizeof(_type),(_max),TRUE)

#define channelTryGet(_chan,_vp,_type) uChannelMultiGet((_chan),(_vp),sizeof(_type),1,FALSE)
#define channelTryPut(_chan,_vp,_type) uChannelMultiPut((_chan),(_vp),sizeof(_type),1,FALSE)

#define channelGetMaxWait(_chan,_vp,_type,_maxSeconds) (uChannelGetMaxWait((_chan),(_vp),sizeof(_type),(_maxSeconds)))
#define channelGetmaxWait channelGetMaxWait  /* old, misspelled name */

 /*   please do not call these functions directly, they are private but in C they
 *   need to be exposed to allow the compiler to understand the public macros.
 */

CHAN *uChannelCreate (int cMax, int size, AC_HANDLE h) ;

void *uChannelGet (CHAN *c, void *vp, int size) ;

int uChannelMultiGet (CHAN *c, void *vp, int size, int max, BOOL wait) ;
int uChannelMultiPut (CHAN *c, void *vp, int size, int max, BOOL wait) ;

int uChannelGetMaxWait (CHAN *c, void *vp, int size, double maxSeconds) ;

/*************************************************************************/
/*************************************************************************/
/* Private details of the Channel struture,
 * please do not use these variables directly, only use the public interface in channel.h
 *
 * The records live in a ring buffer of cMax slots:
 *   head  = slot of the oldest record (next one to be read)
 *   count = number of records currently stored, 0 <= count <= cMax
 * Everything is protected by mutex.
*/

typedef struct channelStruct1 CHAN1 ;

struct channelStruct1 {
  AC_HANDLE h ;
  void *magic ;             /* check that channel exists */
  int size ;                /* size of the records manipulated by this channel */
  int cMax ;                /* max number of unprocessed recors, over this limit channelPut waits */
  int head ;                /* slot number, modulo cMax, of the next record to be read by channelGet */
  int count ;               /* number of records stored in the channel */
  int getWaiting ;          /* number of threads blocked in notEmpty */
  int putWaiting ;          /* number of threads blocked in notFull */
  long nPut, nGet ;         /* global counters: records written and read */
  int nSources ;            /* used by channelAddSources and by channelCloseSource */
  char title[24] ;
  BOOL isClosed ;
  BOOL debug ;
  pthread_mutex_t mutex ;     /* mutex, to be locked before using this channel */
  pthread_cond_t notEmpty ;  /* signal here to activate waiting channelGet calls */
  pthread_cond_t notFull ;   /* signal here to activate waiting channelPut calls */
} ;

struct channelStruct { CHAN1 c ; unsigned char *vp ; } ;

#endif
