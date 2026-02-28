/*  File: call.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: Header file for message system to allow calls by name
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 19 11:06 1998 (fw)
 * * Nov  3 16:15 1994 (mieg): callCdScript, first cd to establish 
     the pwd  of the command, needed for ghostview etc.
 * Created: Mon Oct  3 14:57:16 1994 (rd)
 *-------------------------------------------------------------------
 */

/* $Id: call.h,v 1.8 2015/08/11 22:15:21 mieg Exp $ */


#ifndef DEF_CALL_H
#define DEF_CALL_H
 
#include "regular.h"
#include "aceio.h"

typedef int MESSAGERETURN ;
typedef void (*CallFunc)(va_list args) ;

BOOL call (char *name, ...) ;
int callScript (const char *script, const char *args) ;
int callCdScript (const char *dir, const char *script, const char *args) ; 
FILE* callScriptPipe (const char *script, const char *args) ;
FILE* callCdScriptPipe (const char *dir, const char *script, const char *args) ;
BOOL externalAsynchroneCommand (const char *command, const char *parms,
                                void *look, void(*g)(FILE *f, void *lk)) ;
void externalFileDisplay (const char *title, FILE *f, Stack s) ;
void externalPipeDisplay (const char *title, FILE *f, Stack s) ;
void acedbMailComments(void) ;
void externalCommand (const char* command) ;
int callSystem (const char *cmd) ; /* replacement for 'system' recommended by solaris */
BOOL callExists (const char *name) ;

/* bi directional pipe, missing in SUSE Linux */
int callPipe (const char *command, const char *inBuf, unsigned int inSize, ACEOUT fo) ;

#endif

