/*  File: session_.h
 *  Author: Fred Wobus (fw@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1998
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 11 11:20 1998 (fw)
 * Created: Mon Nov 23 17:21:38 1998 (fw)
 *-------------------------------------------------------------------
 */


#ifndef _SESSION__H
#define _SESSION__H

#ifdef ACEDB5

typedef struct block /* sizeof(block) must be < BLOC_SIZE */
  { 
    DISK gAddress ;
    int mainRelease ,  subDataRelease, subCodeRelease , subSubCodeRelease ;
    int session ;
    char dbName[32] ; /* added as of release 4.3 */
    int magic ; /* added feb 98: add between read and write */
  } BLOCK;   /* super block info */

#else  /* !ACEDB5 */

#define BLKMX (BLOC_SIZE - sizeof(BLOCKHEADER) - sizeof(DISK)\
	       -sizeof(int) -sizeof(int) -sizeof(int) - sizeof(int)\
	       - (32*sizeof(char) -sizeof(int)))

typedef struct block /* sizeof(block) must be < BLOC_SIZE */
  { BLOCKHEADER  h ;
    DISK gAddress ;
    int mainRelease ,  subDataRelease, subCodeRelease ;
    int byteSexIndicator;
    char dbName[32] ; /* added as of release 4.3 */
    char c[BLKMX] ;
    int magic ; /* added feb 98: add between read and write */
  } BLOCK;   /*the transfer unit between disk and cache*/

#endif /* !ACEDB5 */

void swapSuperBlock(BLOCK *bp) ;

/******* session tree **********/
typedef struct STstructure *STP ;
typedef struct STstructure {
  int number;
  KEY key, father, ancester ;
  STP sta ;
  BOOL isDestroyed, isReadlocked, isPermanent ;
  char *user;
  char *date;
  char *title; 
  int generation, x, y, len;
} ST ;

/*********** function priate to the session package ********/

BOOL I_own_the_process(void);

Array sessionTreeCreate(BOOL IsDeadSessionShown); /* array of ST */
void sessionTreeDestroy(Array sessionTree);

void sessionStart(KEY fromSession);
BOOL sessionDoClose (void) ;	/* returns TRUE if 
				   modifications had been saved */
void sessionInitialize (void);	/* set start time of session
				   and init userKey */

VoidRoutine sessionChangeRegister (VoidRoutine func);
/************************************************************/

#endif /* _SESSION__H */
