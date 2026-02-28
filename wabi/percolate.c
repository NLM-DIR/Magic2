/*  File: percolate.c
 *  Author: Danielle et Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg. 2009
 *-------------------------------------------------------------------
 * This file is part of the AceView extension to 
        the ACEDB genome database package, written by
 *         Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *        Jean Thierry-Mieg (CRBM du CNRS, France) mieg@ncbi.nlm.nih.gov
 *
 * Description:
   Align cDNA
 * Exported functions:
 * HISTORY:
 * Created: Thu Sep 14 2009 (mieg)
 *-------------------------------------------------------------------
 */

/* %W% %G% */
#include "acedb.h"
#include "percolate.h"
#include "bs.h"
#include "query.h"

/*******************************************************************/
/*******************************************************************/
typedef struct pcStrut { int a1 ; } PC ;

static BOOL cdnaPercolateGetIntMap (OBJ obj, int *g1p, int *g2p)
{
  BOOL ok = FALSE ;

  return ok ;
} /* cdnaPercolateGetIntMap  */

/*******************************************************************/

static Array cdnaPercolateGetTgExons (KEY tg) 
{ 
  Array aa = arrayCreate (1000, PC) ;
  int g1, g2 ;
  OBJ Tg = bsCreate (tg) ;

  cdnaPercolateGetIntMap (Tg, &g1, &g2) ;
  
  bsDestroy (Tg) ;
  return aa ;
} /* cdnaPercolateGetTgExons */

/*******************************************************************/

static Array cdnaPercolateGetOneMrnaExons (KEY mrna, Array aa) 
{ 
  int g1, g2 ;
  OBJ Mrna = bsCreate (mrna) ;

  cdnaPercolateGetIntMap (Mrna, &g1, &g2) ;

  bsDestroy (Mrna) ;
  return aa ; 
} /* cdnaPercolateGetOneMrnaExons */

/*******************************************************************/

static Array cdnaPercolateGetMrnaExons (KEY tg) 
{ 
  Array aa = arrayCreate (1000, PC) ;
  int ii ;
  KEYSET mrnas = queryKey (tg, ">Mrna") ;

  for (ii = 0; ii < keySetMax (mrnas) ; ii++)
    cdnaPercolateGetOneMrnaExons (keySet (mrnas, ii), aa) ;

  return aa ;
} /* cdnaPercolateGetMrnaExons */

/*******************************************************************/

static int cdnaPercolateTg (KEY tg)
{
  int nn = 0 ;
  Array tgx, mmx ;

  messout ("Percolating gene %s", name (tg)) ;
  
  /* gather the corrdinate and level of support of all exons
   * and introns in the coordinate system of the tg
   */

  tgx = cdnaPercolateGetTgExons (tg) ;
  mmx = cdnaPercolateGetMrnaExons (tg) ;

  /*
    ests = cdnaPercolateCreatePseudoEsts (tgx, mmx) ;
    cdnaPercolateExportEsts (ests) ; 
  */
  arrayDestroy (tgx) ;
  arrayDestroy (mmx) ;
  return nn ;
} /* cdnaPercolateTg */

/*******************************************************************/

int cdnaPercolate (KEYSET tgs)
{
  int ii, jj ;
  KEY tg ;

  for (ii = jj = 0 ; ii < keySetMax (tgs) ; ii++)
    {
      tg = keySet (tgs, ii) ;
      if (cdnaPercolateTg (tg))
	keySet (tgs, jj++) = tg ;
    }
  keySetMax (tgs) = jj ;
  
  messout ("Percolate constructed pseudo ESTs, in %d genes", jj) ;
  return jj ;
} /* cdnaPercolate */

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
