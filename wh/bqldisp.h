/*  File: bqldisp.h
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: public header for bqldispsheet operations
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 17:21 1998 (fw)
 * Created: Thu Aug  6 13:40:24 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: bqldisp.h,v 1.7 2007/03/21 21:03:14 mieg Exp $ */
  
#ifndef DEF_BQLDISP_H  
#define DEF_BQLDISP_H

#include "acedb.h"
void bqldCreate (void) ;
BOOL bqlDisplay (KEY new, KEY old, BOOL isOldGraph) ;

#endif /* DEF_BQLDISP_H */
 
