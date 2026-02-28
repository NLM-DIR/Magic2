/*  File: spreaddata.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:	part of the spreadDisp package
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 22 17:49 1998 (fw)
 * * Aug  4 10:55 1992 (mieg): Working on graph output
 * Created: Thu Jun 11 22:04:35 1992 (mieg)
 *-------------------------------------------------------------------
 */

/* @(#)sprddata.c	1.23 10/29/97 */

#include "acedb.h"
#include "ac.h"

#include "display.h"
#include "spreaddisp_.h"

#include "lex.h"
#include "pick.h"
#include "bs.h"
#include "systags.h"
#include "sysclass.h"
#include "query.h"
#include "plot.h"
#include "freeout.h"

#include <ctype.h>

static void spreadPick(int box, double x, double y) ;
static void spreadSelect(SPREAD spread, int box) ;
static void spreadKeyBoard (int k) ;
static void spreadDispDestroy (void);

extern BOOL bundleDisplay (KEYSET ks, KEYSET ks2, Stack s) ;
extern Array pos2col ; 

/*****************************************/

magic_t GRAPH2SPREAD_ASSOC = "GRAPH2SPREAD";

SPREAD currentSpread (char *caller)
{
  /* find and verify SPREAD struct on active graph */
  SPREAD spread = 0 ;

  if (!(graphAssFind(&GRAPH2SPREAD_ASSOC, &spread)))
    messcrash("%s() could not find SPREAD on graph", caller);
  if (!spread)
    messcrash("%s() received NULL SPREAD pointer", caller);
  if (spread->magic != &SPREAD_MAGIC)
    messcrash("%s() received non-magic SPREAD pointer", caller);
  
  return spread; 
} /* currentSpread */


void spreadDispCreate (BOOL oldGraph)
{
  SPREAD spread ;

  if (!oldGraph && !displayCreate (DtSpreadSheet)) 
    return ;

  graphRegister (DESTROY, spreadDispDestroy) ;

  spread = spreadCreate() ;
  spread->graph = graphActive() ;
  if (oldGraph)
    graphAssRemove (&GRAPH2SPREAD_ASSOC) ;
  graphAssociate (&GRAPH2SPREAD_ASSOC, spread) ;
  spread->sortColonne = 1 ;
  spreadDefineColonne(TRUE) ;

  return;
} /* spreadDispCreate */

/*****************************************/

static void spreadDispDestroy (void)
{    
  Graph mapG ;
  SPREAD spread = currentSpread("spreadDispDestroy");
  
  mapG = spread->mapGraph ;

  spreadDestroy(spread) ;

  if (graphExists(mapG))
    {
      graphActivate(mapG) ;
      graphDestroy() ;
    }

  return;
} /* spreadDispDestroy */

/**********************************************************/
/**********************************************************/

static void spreadExportKeySet (void)
{
  KEYSET ks, ks2 = 0 ; 
  Graph  g = graphActive () ;
  COL *c ;
  int i, col, col2, j = 0, j2 = 0 ;
  static int nExport = 0 ; /* to number the exported lists */
  Array l, t ; 
  KEY key ;
  SPREAD spread = currentSpread("spreadExportKeySet") ;

  if (!arrayMax(spread->colonnes))
    return ;

  col = spread->activeColonne ;
  c = arrayp(spread->colonnes,col,COL) ;

  switch (c->type)
    {
    case 'D':
    case 'P':
    case 'k':
      break ;
    default:
      messout ("The active colonne does not contain keys or DNA or Peptide") ;
      return ;
    }
  switch (c->showType)
    {
    case SHOW_MULTI:
      messout ("Sorry, you cannot export a multi-valued column") ;
      return ;
    default:
      break ;
    }

  t = spread->tableau ;
  ks = keySetCreate () ;
  if (c->type == 'D' || c->type == 'P') 
    {
      col2 = c->from - 1 ;
      /*       c2 = arrayp(spread->colonnes,col2,COL) ; */
      ks2 = keySetCreate () ;
      
      for (i = 0; i < arrayMax(spread->tableau) ; i++)
	if (
	    ( l = array(t, i, Array)) &&
	    ( key = arr(array(t, i, Array) , col, SPCELL).u.k)
	    )
	  {
	    Array a = array(t, i, Array) ;
	    keySet(ks, j++)  = key ;	
	    key = arr(a , col2, SPCELL).u.k ;
	    if (!key)
	      messcrash ("pb avec le keyset des sequences") ;
	    keySet(ks2, j2++)  = key ;	
	  }
    }
  else 
    {
      for (i = 0; i < arrayMax(spread->tableau) ; i++)
	if ( 
	    ( l = array(t, i, Array)) &&
	    ( key = arr(array(t, i, Array) , col, SPCELL).u.k) 
	    )
	  if (!arr(l, col, SPCELL).empty)
	    keySet(ks, j++)  = key ;	
    }
  if (!arrayMax(ks))
    {
      messout ("The active colonne is empty, sorry") ;
      keySetDestroy(ks) ;
      return ;
    }
  switch (c->type)
    {
    case 'D': 
      bundleDisplay (ks, ks2, c->dnaStack) ;
      graphRetitle("Exported DNA") ;
      keySetDestroy (ks) ;
			keySetDestroy (ks2) ;
      break ;

    case 'P': 
      bundleDisplay (ks, ks2, c->pepStack) ;
      graphRetitle("Exported Peptides") ;
      keySetDestroy (ks) ;
			keySetDestroy (ks2) ;
      break ;

    case 'k':
      keySetSort(ks) ;
      keySetCompress(ks) ;
      displayCreate(DtKeySet) ;
      graphRetitle(messprintf ("%d:%s", ++nExport, *c->subtitleBuffer ? c->subtitleBuffer : "Exported keyset")) ;
      keySetShow (ks,0) ; ks = 0 ; /* will destroy ks */
      keySetSelect () ;
      graphActivate (g) ;
    }

  return;
} /* spreadExportKeySet */

/*****************************************************/

void spreadImportKeySet (void)
{
  KEYSET ks = 0, ks2 = 0 ;
  void *vp ;
  COL *c ;
  SPREAD spread = currentSpread("spreadImportkeySet") ;

  if (!arrayMax(spread->colonnes))
    spreadInitColonne(spread) ;

  if (!spreadCheckConditions(spread))
     return ;

  if (!arrayMax(spread->colonnes))
    { messout ("// First define the table") ;
    return ;
    }

  c = arrp(spread->colonnes, 0, COL) ;
  if (!c->classe)
    { messout ("First choose the class of column 1 ") ;
      return ;
    }
  
  if (*c->conditionBuffer &&
      ! condCheckSyntax(messprintf(" %s", c->conditionBuffer)))
    { messout ("First fix the syntax error in column 1's condition") ;
      return ;
    }

  if (!keySetActive(&ks, &vp))
    { messout("First select a keySet, thank you") ;
      return ;
    }
 
  ks2 = spreadFilterKeySet (spread, ks) ;

  spread->isActiveKeySet = TRUE ;
  spreadRecomputeKeySet(spread, ks2, 0, 0, 0) ;
  spread->isActiveKeySet = FALSE ;
  keySetDestroy (ks2) ;
  spreadDisplayData(spread) ;

  return;
} /* spreadImportKeySet */

/************************************************/

void spreadSwitchColonnes (void)
{
  int c1, c2, p1, p2 ;
  COL *c ;
  SPREAD spread = currentSpread("spreadSwitchColonnes") ;

  c1 = spread->activeColonne ;
  p1 = arrayMax(spread->pos2col) ;  /* search the corresponding
				       active position */
  while (p1--)
    if (c1 == arr(spread->pos2col, p1, int) )
      break ;
  p2 = p1 ;
  while(TRUE)
    {
      p2++ ;   /* The position i ll switch with p1 */
      if (p2 >= arrayMax(spread->pos2col)) 
	return ;

      c2 =  arr(spread->pos2col, p2, int) ;
      c = arrp(spread->colonnes, c2, COL) ;
      
      if (!c->hidden)
	break ;
    }
                  /* Do switch */
  arr(spread->pos2col, p1, int) = c2 ;   
  arr(spread->pos2col, p2, int) = c1 ;
	
  spreadReorder(spread) ;
  spreadDisplayData(spread) ;

  return;
} /* spreadSwitchColonnes */
/*****************/
/* mhmp 01.04.99 pour smooth */
static float cleanScale(float scale)
{
  int inter;
  inter = 1 ;
  if (scale < 1) {
    inter = 1./scale ;
    inter = utMainRoundPart (inter) ;
    return 1./inter ;
  }
  else scale = 1;
  return scale ;
}

/*****************/

static BOOL spreadHistoFillArray (SPREAD spread, int col, Array a, 
				  float *xminp, float *xmaxp, float *scalep, int *xfacp)
{
  int n, n1, xfac = 1 ;
  int step = 1 ;
  /* mhmp 14.09.01 step arrondi */
  int debut ;
  float x = 0, scale = 0, xmin = 0, xmax = 0, d = 0 ;
  Array tt ;
  COL *c1 ;

  c1 = arrp( spread->colonnes, col, COL) ;
  tt = spread->tableau ;
  n = arrayMax (tt) ;
  switch (c1->realType)
    {
    case 'c': /*mhmp 22.02.99 + 'c' */
      scale = 0 ; xmin = 1; xmax = 0 ;
      while (n--) 
	{ 
	  x = (arr (arr(spread->tableau, n, Array), col, SPCELL)).u.i ;
	  if (xmin > xmax) xmin = xmax = x ;
	  else { if (xmin > x) xmin = x ; if(xmax < x) xmax = x ; } 
	}
      if (xmax - xmin < 10000) scale = 1.0 ;
      else scale = 10000.0 / (utArrondi(xmax) - utArrondi(xmin)) ;
      scale = cleanScale (scale) ;/* mhmp 01.04.99 pour smooth */
      if (scale != 1){
	step = utArrondi(1./scale) ;
	if(!messQuery(messprintf("%d-stepping recommended\n Etes-vous d'accord?",step)))
	  scale = 1 ;
      }
      n = arrayMax(tt) ;
      while (n--) 
	{
	  x = arr (arr(spread->tableau, n, Array), col, SPCELL).u.i ;
          n1 = (x - xmin ) * scale ;
	  array(a, n1, int)++ ;
	}
      break ;
    case 'i':
    case 'f':
      scale = 0 ; xmin = 1; xmax = 0 ;
      while (n--) 
	{ 
	  x = (arr (arr(spread->tableau, n, Array), col, SPCELL)).z ;
	  if (xmin > xmax) xmin = xmax = x ;
	  else { if (xmin > x) xmin = x ; if(xmax < x) xmax = x ; } 
	}
      if (xmax - xmin < 0.0001) xmax = xmin + 0.0001;     
      d = xmax - xmin ;
      while (d < 10)
	{
	  xfac *= 10 ;
	  d *= xfac ;
	}
      xmin *= xfac ;
      xmax *= xfac ;
      scale = 10000.0 / (utArrondi(xmax) - utArrondi(xmin)) ;
      scale = cleanScale (scale) ;
      if (scale != 1){
	step = utArrondi(1./scale) ;
	if(!messQuery(messprintf("%d-stepping recommended\n Etes-vous d'accord?",step)))
	  scale = 1 ;
      }
      n = arrayMax(tt) ;
      while (n--)
	{
	  x = arr (arr(spread->tableau, n, Array), col, SPCELL).z ;
	  x *= xfac ;
	  n1 = (x - xmin ) * scale ;
	  if (n1 >= 0)
	    array(a, n1, int)++ ;
	}
      break ;
    default:
      return FALSE ;
    }
  /* mhmp 14.09.01 step arrondi */
  if (step > 1){
    debut = xmin ;
    if (xmin >= 0) {
      debut = (debut/step) * step ;
      xmin = debut ;
    }
    else {
      debut = (-debut/step) * step + step ;
      xmin = - debut ;
    }
  }

  *xminp = xmin ; *xmaxp = xmax ; *scalep = scale ; *xfacp = xfac ;
  return TRUE ;
}

static void spreadHistoCreate (void)
{
  Array a = 0 ;
  int col, xfac ;
  float xmin, xmax, scale ;
  COL* c1 ;
  char *cp, *cq ;
  SPREAD spread = currentSpread("spreadHistoCreate") ;

  col = spread->activeColonne ;
  c1 = arrp( spread->colonnes, col, COL) ;
  cp = spread->titleBuffer ;
  cq = c1->subtitleBuffer ;
  a = arrayCreate (1000, int) ;
  if (spreadHistoFillArray (spread, col, a, &xmin, &xmax, &scale, &xfac))
    plotShiftedHisto (cp,cq, a, xmin, xmax, scale, xfac) ;
  else
    {
      messout ("First select a numeric column") ;
      arrayDestroy (a) ;
    }
  return  ;
}

/*****************/

static void spread2DHistoCreate (void)
{
  Array xy = 0, tt = 0 ;
  int n, i, col0 = 0, col1 = 0, col2 = 0 ;
  COL *c0, *c1, *c2 ;
  char *cp=0, *cp1=0, *cp2=0 ;
  SPCELL *sp0, *sp1, *sp2 ;
  POINT2D *pp ;   /* a 2 dim point */
  SPREAD spread = currentSpread("spreadHistoCreate") ;

  col0 = spread->activeColonne ;
  c0 = arrp( spread->colonnes, col0, COL) ; /* mhmp 20.06.02 col1->col0 */ 
  cp = spread->titleBuffer ;
  c1 = c2 = 0 ;
  switch (c0->realType)
    {
    case 'i':  case 'c': case 'f':
      c1 = c0 ; col1 = col0 ; /* c0 = 0 ;  mhmp 20.06.02 */
      if (col0) col0-- ;  /* a wild guess */
      break ;
    case 'k':  /* we found a column of keys, search now for 2 col of numbers */
      break ;
    default:
      goto usage ;
      break ;
    }
  if (!c1)
    {
      for (i = 0 ; i < arrayMax(pos2col) ; i++)
	if (col0 ==  arr(pos2col,i, int))
	  break ;
      for (i++; i < arrayMax(pos2col) ; i++)
	{
	  col1 =  arr(pos2col,i, int) ;
	  c1 = arrp( spread->colonnes, col1, COL) ;
	  if (c1->hidden)
	    continue ;
	  cp1 = c1->subtitleBuffer ;
	  switch (c1->realType)
	    {
	    case 'i':  case 'c': case 'f':
	      break ;
	    default:
	      c1 = 0 ;
	      break ;
	    }
	  if (c1) break ;
	}
    }
  if (c1)
    {
      for (i = 0 ; i < arrayMax(pos2col) ; i++)
	if (col1 ==  arr(pos2col,i, int))
	  break ;
      for (i++; i < arrayMax(pos2col) ; i++)
	{
	  col2 =  arr(pos2col,i, int) ;
	  c2 = arrp( spread->colonnes, col2, COL) ;
	  if (c2->hidden)
	    continue ;
	  cp2 = c2->subtitleBuffer ;
	  switch (c2->realType)
	    {
	    case 'i':  case 'c': case 'f':
	      break ;
	    default:
	      c2 = 0 ;
	      break ;
	    } 
	  if (c2) break ;
	}
    }
  if (!c1 || !c2)
    goto usage ;

  cp1 = c1->subtitleBuffer ;

  tt = spread->tableau ;
  n = arrayMax (tt) ;
 
  xy = arrayCreate (n, POINT2D) ;
  while (n--)
    {
      pp = arrayp (xy, n, POINT2D) ;
      pp->x = pp->y = 0 ;
      sp0 = arrp(arr(spread->tableau, n, Array), col0, SPCELL) ;
      sp1 = arrp(arr(spread->tableau, n, Array), col1, SPCELL) ;
      sp2 = arrp(arr(spread->tableau, n, Array), col2, SPCELL) ;
      
      if (!sp0->empty) 
	switch (c0->realType)
	  {
	  case 'k': 
	    pp->k = sp0->u.k ;
	    break ;
	  }
      if (!sp1->empty) 
	switch (c1->realType)
	  {
	  case 'i': case 'c':
	    pp->x = sp1->z ;
	    break ;
	  case 'f':
	    pp->x = sp1->z ;
	    break ;
	  }
      if (!sp2->empty) 
	switch (c2->realType)
	  {
	  case 'i': case 'c':
	    pp->y = sp2->z ;
	    break ;
	  case 'f':
	    pp->y = sp2->z ;
	    break ;
	  }
    }

  plot2d (cp, cp1, cp2, xy) ;
  return ;

 usage:
  messout ("%s\n%s",
	   "First select a column of keys, followed by 2 col of numbers",
	   " or a col of numbers followed by another one") ;
  return ;
}

/*************************************************************************************/
/*************************************************************************************/
#define NNh 100

static double logfac (int n)
{
  int i ;
  double z ;

  if (n < 100)
    for (i = 1, z = 1 ; i <= n ; i++)
      z += log (i) ;
  else /* Stirling formula */
    z = n * log(n) - n ;
  return z ;
}

/*************************************************************************************/

static double *sprdHistoSmoothDistrib (int m, int n)
{
  static AC_HANDLE h = 0 ;
  int i, j, k, m1 ;
  static Array aaa = 0 ;
  static double *abig = 0 ;
  double z, zmax, *aa, NNz = NNh, p, q, t ;

  if (!h) 
    h = ac_new_handle () ;

  if (m > n) { i = m ; m = n ; n = i ; }
  if (n > 10000 || 400*m < n || 400 * (n-m) < n)
    {
      if (!abig) 
	abig = halloc ((NNh+1)*sizeof (double), h) ;
      aa = abig ;
      z = .5 + 100.0*m/(double)n ;
      i = z ;
      memset (aa, 0, (NNh+1)*sizeof (double)) ; aa[i] = t = 1 ;
    }
  else if (n > 400)
    {
      if (!abig) 
	abig = halloc ((NNh+1)*sizeof (double), h) ;
      aa = abig ;
      aa[0] = aa[NNh] = t =  0 ;
      for (i = 1, zmax = 0 ; i < NNh ; i++)
	{
	  p = i/NNz ; q = 1 - p ;
	  z = 0 ;
	  z += m * log (p) ;
	  z += (n-m) * log (q) ;
	  z += logfac(n) ;
	  z -= logfac(m) ;
	  z -= logfac(n-m) ;
	  aa[i] = z ;
	  if (i == 1 || z > zmax) zmax = z ;
	}
      for (i = 0 ; i <= NNh ; i++)  
	{
	  aa[i] -= zmax ;
	  aa[i] = exp (aa[i]) ;
	  t += aa[i] ;
	}
      for (i = 0 ; i <= NNh ; i++)  	  
	aa[i] /= t ;
    }
  else
    {
      if (! aaa)
	aaa = arrayCreate (NNh * NNh, double*) ;
      k = m + (n-1) * (n-1) ;
      aa = array (aaa, k, double *) ;
      if (! aa)
	{
	  aa = array (aaa, k, double *) =  halloc ((NNh+1)*sizeof (double), h) ;
	  for (i = 0, t = 0 ; i <= NNh ; i++)
	    {
	      /* Cmn = n!/m! (n-m)! p^m q^(n-m) */
	      p = i/NNz ; q = 1 - p ;
	      if (m > n/2)
		{ z = p ; p = q ; q = z ; m1 = n - m ; }
	      else
		m1 = m ;
	      for (z = 1, j = 1 ; j <= m1 ;j++)
		z = z * p * (n + 1 - j)/j ; 
	      for (j = 1 ; j <= n - m1 ; j++)
		z = z * q ;
	      aa[i] = z ; 
	      t += z ;
	    } 
	  for (i = 0 ; t> 0 && i <= NNh ; i++)  
	    aa[i] /= t ;
	}
    }
  
  return aa ;
}  /* sprdHistoSmoothDistrib */

/*****************/

static void spreadSmoothHistoCreate (void)
{
  Array a, xy = 0, tt = 0 ;
  int n, nerr = 0, i, col0 = 0, col1 = 0, col2 = 0, *ip ;
  COL *c0, *c1, *c2 ;
  double *aa ;
  char *cp2=0 ;
  SPCELL *sp0, *sp1, *sp2 ;
  POINT2D *pp ;   /* a 2 dim point */
  SPREAD spread = currentSpread("spreadSmoothHistoCreate") ;

  col0 = spread->activeColonne ;
  c0 = arrp( spread->colonnes, col0, COL) ; /* mhmp 20.06.02 col1->col0 */ 
  c1 = c2 = 0 ;
  switch (c0->realType)
    {
    case 'i':  case 'c': case 'f':
      c1 = c0 ; col1 = col0 ; /* c0 = 0 ;  mhmp 20.06.02 */
      if (col0) col0-- ;  /* a wild guess */
      break ;
    case 'k':  /* we found a column of keys, search now for 2 col of numbers */
      break ;
    default:
      goto usage ;
      break ;
    }
  if (!c1)
    {
      for (i = 0 ; i < arrayMax(pos2col) ; i++)
	if (col0 ==  arr(pos2col,i, int))
	  break ;
      for (i++; i < arrayMax(pos2col) ; i++)
	{
	  col1 =  arr(pos2col,i, int) ;
	  c1 = arrp( spread->colonnes, col1, COL) ;
	  if (c1->hidden)
	    continue ;
	  switch (c1->realType)
	    {
	    case 'i':  case 'c': case 'f':
	      break ;
	    default:
	      c1 = 0 ;
	      break ;
	    }
	  if (c1) break ;
	}
    }
  if (c1)
    {
      for (i = 0 ; i < arrayMax(pos2col) ; i++)
	if (col1 ==  arr(pos2col,i, int))
	  break ;
      for (i++; i < arrayMax(pos2col) ; i++)
	{
	  col2 =  arr(pos2col,i, int) ;
	  c2 = arrp( spread->colonnes, col2, COL) ;
	  if (c2->hidden)
	    continue ;
	  cp2 = c2->subtitleBuffer ;
	  switch (c2->realType)
	    {
	    case 'i':  case 'c': case 'f':
	      break ;
	    default:
	      c2 = 0 ;
	      break ;
	    } 
	  if (c2) break ;
	}
    }
  if (!c1 || !c2)
    goto usage ;

  tt = spread->tableau ;
  n = arrayMax (tt) ;
 
  xy = arrayCreate (n, POINT2D) ;
  while (n--)
    {
      pp = arrayp (xy, n, POINT2D) ;
      pp->x = pp->y = 0 ;
      sp0 = arrp(arr(spread->tableau, n, Array), col0, SPCELL) ;
      sp1 = arrp(arr(spread->tableau, n, Array), col1, SPCELL) ;
      sp2 = arrp(arr(spread->tableau, n, Array), col2, SPCELL) ;
      
      if (!sp0->empty) 
	switch (c0->realType)
	  {
	  case 'k': 
	    pp->k = sp0->u.k ;
	    break ;
	  }
      if (!sp1->empty) 
	switch (c1->realType)
	  {
	  case 'c':
	    pp->x = sp1->u.i ;
	    break ;
	  case 'i':
	  case 'f':
	    pp->x = sp1->z ;
	    break ;
	  }
      if (!sp2->empty) 
	switch (c2->realType)
	  {
	  case 'c':
	    pp->y = sp2->u.i ;
	    break ;
	  case 'i':
	  case 'f':
	    pp->y = sp2->z ;
	    break ;
	  }
    }
  
  a = arrayCreate (100, int) ;

  for (n = nerr = 0 ; n < arrayMax (xy) ; n++)
    {
      pp = arrayp (xy, n, POINT2D) ;
      if (pp->x < 0 || pp->y < 0 || pp->x > pp->y)
	{ nerr++ ; continue ; }
      aa = sprdHistoSmoothDistrib (pp->x, pp->y) ;
      for (i = 0 ; i <= NNh ; i++)
	{
	  ip = arrayp (a, i, int) ;
	  (*ip) += 1000 * aa[i] ;
	}
    }
  ac_free (xy) ;
  plotShiftedHisto ("Smoothed histogram",cp2, a, 0, 100, 1, 1) ;
  return ;

 usage:
  messout ("%s\n%s",
	   "First select a column of keys, followed by 2 col of numbers",
	   " or a col of numbers followed by another one") ;
  return ;
}

/*****************/

static void spreadAsciiDumpDo (SPREAD spread)
{
  static char sep = '\t' ; 
  FILE *f ;
  int level ;

  f = filqueryopen (0, 0, "txt", "w", 
		    "Where should I ascii-dump this table ?") ;
  if (!f)
    return ;
 
  level = freeOutSetFile (f) ;
  spreadDoDump(spread, sep, 0, TRUE) ; /* style 0, like 'a' but no format line */
  freeOutClose (level) ;
  filclose(f) ;

  return;
} /* spreadAsciiDumpDo */

/*****************/

static void spreadAsciiDump(void)
{
  SPREAD spread = currentSpread("spreadAsciiDump") ;
  spread->showTitle = FALSE ;
  spreadAsciiDumpDo (spread) ;

  return ;
} /* spreadAsciiDump */

/*****************/

static void spreadAsciiDumpT(void)
{
  SPREAD spread = currentSpread("spreadAsciiDump") ;
  spread->showTitle = TRUE ;
  spreadAsciiDumpDo (spread) ;

  return ;
} /* spreadAsciiDumpT */

/*****************/

static void spreadHtmlDump(void)
{
  static char sep = '\t' ; 
  FILE *f ;
  int level ;
  SPREAD spread = currentSpread("spreadAsciiDump") ;
  BOOL old ;
  
  f = filqueryopen (0, 0, "htm", "w", 
		    "Where should I html-dump this table ?") ;
  if (!f)
    return ;
 
  level = freeOutSetFile (f) ;
  freeOut ("<html>\n<body bgcolor=white>\n<h2>") ;
  if (*spread->titleBuffer) freeOutf ("<h2>\n%s\n</h2>", spread->titleBuffer) ;
  spreadDumpLegend (spread, 'x') ;
  
  freeOut ("<p>\n<table border=1>\n") ;

  old = spread->showTitle ;
  spread->showTitle = TRUE ;
  spreadDoDump(spread, sep, 'x', TRUE) ; 
  spread->showTitle= old ;

  freeOut ("</table>\n</body>\n</html>\n") ;
  freeOutClose (level) ;
  filclose(f) ;

  return;
} /* spreadHtmlDump */

/**********************************************************/

static void shouldWeDestroy(void)
{
  if (messQuery("Do you really want to quit the Table_Maker ?"))
    graphDestroy() ;

  return;
}


/*mhmp 21.04.98 pour avoir graphClear correct */
static void spreadDefineColonne2(void)
{
  spreadDefineColonne(TRUE) ; 
  /* spreadDefineColonne(TRUE) ; */

  return;
} /* spreadDefineColonne2 */


/************/

static MENUOPT spreadMenu[] =
  {
   {shouldWeDestroy, "Quit"},
   {help, "Help"},
   {graphPrint, "Print"},
   {spreadExportKeySet, "Export column"},
   {spreadAsciiDump, "Export data in tab delimited format"},
   {spreadAsciiDumpT, "Export data in tab delimited, with title"},
   {spreadHtmlDump, "Export data as html table"},
   {spreadDefineColonne2,"Edit Query"},
   {spreadSwitchColonnes, "Switch columns"},
   {spreadHistoCreate, "Histo"},
   {spreadSmoothHistoCreate, "Smooth Histo"},
   {spread2DHistoCreate, "2D-Plot"},
   {spreadMapCreate, "Comparative map"},
   {0, 0} 
   } ;

/*****************/

void spreadDisplayData(SPREAD spread)
{
  int i , j , j1, jj, maxCol , box, color = 0, col = 4 , line , nn = 0 , oldFormat = 0 ;
  COL *c ;
  Array lineArray , oldLineArray = 0 ; 
  SPCELL *su ;
  BSunit *u, *uOld = 0 ;
  char *cp ;
  char small[130] ;
  char timeBuf[25] ;

  BOOL isAnySubtitle = FALSE ;
  float y2 ;
  int maxLines = 1000 ;

  graphUnMessage () ;
  pos2col = spread->pos2col ;
  
  graphActivate(spread->graph) ;
  graphPop() ;
  graphClear() ;
  graphRegister (PICK, spreadPick) ;
  graphRegister (KEYBOARD, spreadKeyBoard) ;
  graphHelp("Table_Maker_Data") ;
  graphMenu (spreadMenu) ;

  if(arrayExists(spread->tableau) && arrayMax(spread->tableau))
    maxLines = arrayMax(spread->tableau);
  if (maxLines > 1000)
    {
      graphMessage (
"Attention, the graphic display is limited to 1000 lines, but the full table can be queried and printed") ;
      maxLines = 1000;
    }

  line = 1 ;

  if (stackExists(spread->comments) && !stackAtEnd(spread->comments))
    { stackCursor(spread->comments,0) ;
      while ((cp = stackNextText(spread->comments)))
	{ graphText(cp, 2, line++) ;
	}
      ++line ;
    }
  
#ifdef EDITBUTTON
  spread->editModeBox = graphButton("Edit", spreadEditButton, 1, 1.2) ;
  graphBoxDraw(spread->editModeBox,BLACK, 
	       spread->editMode ? LIGHTBLUE : WHITE ) ;
#endif

  i = graphBoxStart () ;
  graphButtons (spreadMenu, 1, line, 60) ;
  graphBoxEnd () ;
  graphBoxDim (i, 0, 0, 0, &y2) ;
  line = y2 + 2 ;

  if (*spread->titleBuffer)
    { graphText(spread->titleBuffer, 2, line) ;
      line += 2 ;
    }

  spread->activeBox = 0 ;
  maxCol = arrayMax(spread->colonnes) ;
	spread->tableauBox = graphBoxStart() ;
  if (arrayExists(spread->tableau) && arrayMax(spread->tableau))
    { spread->flags =   /* duplicated in sprdDoDump */
	arrayReCreate(spread->flags, arrayMax(spread->tableau), char) ;
      array(spread->flags, arrayMax(spread->tableau) - 1 , char) = 0 ;
				    
      for (j = 0, col = 4; j < maxCol; j++)
	{
	  j1 = arr(pos2col,j, int) ;
	  c = arrp(spread->colonnes,j1,COL) ;
	  if (!c->hidden)
	    { graphBoxStart() ; /* mhmp 20.11.02 */
				if (*c->subtitleBuffer)
				{ graphText(c->subtitleBuffer, col, line) ;
				  isAnySubtitle = TRUE ;
				}
			  col += c->width ;
				graphBoxEnd() ;
	    }
	}
      if (isAnySubtitle)
	line += 2 ;

      for (i = 0 ; i < maxLines ; i++)
	{ lineArray = arr(spread->tableau, i, Array) ;
	  if (oldLineArray && !spreadOrder(&oldLineArray, &lineArray))
	    { arr(spread->flags, i, char) = 1 ; continue ; }
	          /* Useful if some colonnes are hidden */
	  nn++ ;
	  oldLineArray = lineArray ;
	  for(j = jj = 0,  col = 4; j < maxCol; j++)
	    {
	      j1 = arr(pos2col,j, int) ;
	      c = arrp(spread->colonnes,j1,COL) ;
	      if (!c->hidden)
		{ 
		  su = arrp(lineArray, j1,SPCELL) ;
		  u = &(su->u) ;
		box = graphBoxStart(); 
		if (!jj++) 
		  {
		    if (!uOld || u->k != uOld->k) color = 1 - color ;
		    uOld  = u ;
		  }
		/* create box even in empty case, otherwise fools Pick */
		if (!arr(lineArray, j1,SPCELL).empty) 
		  {
		    if (c->showType == SHOW_MULTI)
		      {
			if (stackExists(c->text) && u->i)
			  graphText(stackText(c->text, u->i), col, line) ;
		      }
		    else switch (c->type)
		      {
		      case 0:
			break ;
		      case 'c':
			graphText(messprintf("%d", u->i), col, line) ;
			break ;
		      case 'i': 
			graphText(messprintf("%lld", (long long int)su->z), col, line) ;
			break ;
		      case 'f':
			{
			  double zf = su->z > 0 ? su->z : - su->z ;
			  long long int izf = su->z ;
			  if (izf > 0 && zf - izf < ACE_FLT_RESOLUTION)
			    graphText(messprintf("%lld", su->z > 0 ? izf : -izf), col, line) ;
			  else
			    graphText(messprintf("%6lg", su->z), col, line) ;
			}		
			break ;
		      case 'd':
			if (u->time)
			  graphText(timeShow (u->time, timeBuf, 25), col, line) ;
			break ;
		      case 'k': case 'K': case 'n': case 'b':
			if (iskey(u->k) == 2)
			  oldFormat = graphTextFormat(BOLD) ;
			if (iskey(u->k))
			  graphText(name(u->k), col, line) ;
			if (iskey(u->k) == 2)
			  graphTextFormat(oldFormat) ;
			break ;
		      case 't':
			if (stackExists(c->text) && u->i)
			  graphText(stackText(c->text, u->i), col, line) ;
			break ;
		      case 'D':
			if (stackExists(c->dnaStack) && u->k) {
			  memset(small, 0, 130) ;
			  strncpy(small, stackText(c->dnaStack, u->k), 127) ;
			  if (strlen(stackText(c->dnaStack, u->k)) > 127)
			    strcat(small, "..") ;
			  graphText(small, col, line) ;
			}
			break ;
		      case 'P':
			if (stackExists(c->pepStack) && u->k) {
			  memset(small, 0, 130) ;
			  strncpy(small, stackText(c->pepStack, u->k), 127) ;
			  if (strlen(stackText(c->pepStack, u->k)) > 127)
			    strcat(small, "..") ;
			  graphText(small, col, line) ;
			}
			break ;
		      default:
			messcrash("Unknown type in spreadDisplayData") ;
		      }
		  }
		graphBoxEnd() ;
		if (color)
		  graphBoxDraw (box, BLACK, PALEBLUE) ;
		col += c->width ;
		}
	    }
	  ++line ;
	}
    }	  
  graphBoxEnd() ;
  if(!nn || maxLines ==  arrayMax(spread->tableau)) 
    graphText(messprintf("%d lines", nn), 35, 0.0) ;
  else
    graphText(messprintf("%d lines, 1000 shown",arrayMax(spread->tableau) ), 35, 0.0) ;
  graphTextBounds (col + 8, line + 3) ;  

  i = 0 ;
  for(j = 0 ,  col = 4; j < maxCol; j++)
    {
      j1 = arr(pos2col,j, int) ;
      c = arrp(spread->colonnes,j1,COL) ;
      if (!c->hidden)
	i++ ;
    }
  spread->numberDisplayedColonnes = i ;
  graphRedraw() ;

  if (graphExists(spread->mapGraph))
    { spreadMapCreate() ;
      graphActivate(spread->graph) ;
      graphPop() ;
    }

  return;
} /* spreadDisplayData */

/*****************/
/*****************/

static void spreadSelect(SPREAD spread, int box)
{
  int 
    j, j1, n = box - spread->tableauBox - 1 ,
    max = spread->numberDisplayedColonnes ,
    maxCol = arrayMax(spread->colonnes) ;
  COL *c ;
  SPCELL *su ;
  BSunit u ;
  char timeBuf[25] ;
  char *cp ;

  if (n >= max)
    n = n - max ; /* mhmp 20.11.02 */
  if (spread->activeBox > spread -> tableauBox)
    graphBoxDraw (spread->activeBox, BLACK, WHITE) ;
   
  spread->activeBox = box ;
  graphBoxDraw (spread->activeBox, BLACK, RED) ;

  j1 = n / max ;
  for(j = 0 , cp = arrp(spread->flags, 0, char) ; 
      j < arrayMax(spread->tableau); cp++, j++)
    if (!*cp)
      { if (!j1--)
	  { spread->activeLine = arr(spread->tableau, j, Array) ;
	    break ;
	  }
      }
  n = n%max ;
  for(j = 0 ; j < maxCol; j++)
    {
      j1 = arr(pos2col,j, int) ;
      c = arrp(spread->colonnes,j1,COL) ;
      if (!c->hidden)
	{ if (!n--)
	    { spread->activeColonne  = j1 ;  /* may be j */
	      break ;
	    }
	}
    }
	  
  c = arrp( spread->colonnes, spread->activeColonne, COL) ;
  su = arrp (spread->activeLine, spread->activeColonne, SPCELL) ;
  u = su->u ;
  
   if (c->showType == SHOW_MULTI)
     graphPostBuffer (stackExists(c->text) && u.i ? stackText(c->text, u.i) : " " ) ;
   else switch (c->type)
    {
    case 'b':
      break ;
    case 'k': case'n': case 'K':
      graphPostBuffer (name(u.k)) ;
      break ;
    case 'i':
      graphPostBuffer (messprintf("%lld", (long long int) su->z)) ;
      break ;
    case 'f':
      graphPostBuffer (messprintf("%lg", su->z)) ;
      break ;
    case 'd':
      graphPostBuffer (timeShow (u.time, timeBuf, 25)) ;
      break ;
    case 't':
      graphPostBuffer (stackExists(c->text) && u.i ? stackText(c->text, u.i) : " " ) ;
      break ;
    default:
      break ;
    }

  return;
} /* spreadSelect */

void spreadSelectFromMap(SPREAD spread, int line, int colonne)
{
  int 
    j, j1, j2,
    max = spread->numberDisplayedColonnes ,
    maxCol = arrayMax(spread->colonnes) ;
  COL *c ;
  char *cp ;

  for(j = 0 , j1 = 0, cp = arrp(spread->flags, 0, char) ; 
      j < arrayMax(spread->tableau); cp++, j++)
    if (!*cp)
      { if (j == line)
	  break ;
	j1++ ;
      }
  if (j == arrayMax(spread->tableau))
    return ; /* failure */
  
  for(j = 0 ; j < maxCol; j++)
    {
      j2 = arr(pos2col,j, int) ;
      c = arrp(spread->colonnes,j2,COL) ;
      if (!c->hidden)
	{ if (j2 == colonne)
	    { spreadSelect(spread, max*j1 + j  +  spread->tableauBox + 1 ) ;
	      return ;
	    }
	}
    }

  return;
} /* spreadSelectFromMap */

/*
static void spreadEditButton(void)
{
  SPREAD spread = currentSpread ("spreadEdit") ;

  spread->editMode = spread->editMode  ? FALSE : TRUE ;
  graphBoxDraw(spread->editModeBox,BLACK, 
	       spread->editMode ? LIGHTBLUE : WHITE ) ;
}
*/
static Graph sprdAddKeyGraph = 0 ;
static BOOL isAddingKey = FALSE ;
static void sprdAddKey1 (KEY key)
{ 
  BSunit *u ;
  SPREAD spread = currentSpread ("sprdAddkey1") ;
  
  u = &(arr(spread->activeLine, spread->activeColonne, SPCELL).u) ; 
  if (class(key) != class(u->k))
    { messout("The selected object is not in the correct class, sorry") ;
      return ;
    }
  if (u->k != key)
    { u->k = key ;
      spreadReorder(spread) ;
      spreadDisplayData(spread) ;
    }

  return;
} /* sprdAddKey1 */


static void sprdAddKey (KEY key)
{ 
  graphActivate(sprdAddKeyGraph) ;
  graphRegister (MESSAGE_DESTROY, 0) ;
  graphUnMessage() ; /* this if i came from short cut sprdPick */
  displayUnBlock() ;
  isAddingKey = FALSE ;
  sprdAddKey1(key) ;

  return;
} /* sprdAddKey */


static void sprdAddKeyByName1 (void)
{
  char *cp ;
  KEY kk ;
  BSunit *u ;
  SPREAD spread = currentSpread ("sprdAddkey1") ;
  
  u = &(arr(spread->activeLine, spread->activeColonne, SPCELL).u) ; 
  if (!messPrompt("New value ", name(u->k),"w"))
    return ;
  cp = freeword() ;
		    
  if (!lexword2key(cp, &kk, class(u->k)))
    { if (messQuery("Unknown object, do you want to create it ?"))
	lexaddkey(cp, &kk, class(u->k)) ;
      else
	return ;
    }
  u->k = kk ;
  
  spreadReorder(spread) ;
  spreadDisplayData(spread) ;
  
  return;
} /* sprdAddKeyByName1 */


static void sprdAddKeyByName (void)
{ 
  graphActivate(sprdAddKeyGraph) ;
  isAddingKey = FALSE ;
  displayUnBlock() ;
  sprdAddKeyByName1() ;
  
  return;
} /* sprdAddKeyByName */


static void spreadFollow (SPREAD spread, int box)
{
  COL* c = arrp( spread->colonnes, spread->activeColonne, COL) ;
  SPCELL *su = arrp(spread->activeLine, spread->activeColonne, SPCELL) ;
  BSunit *u = &(su->u) ;
  BOOL touched = FALSE ;

  if (spread->editMode)
    switch (c->type)
      {
      case 'b':
	u->k = u->k ? 0 : c->tag ;
	touched = TRUE ;
	break ;
      case 'k': case'n': case 'K':
	isAddingKey = TRUE ;
	sprdAddKeyGraph = graphActive() ;
	displayBlock (sprdAddKey, "Or cancel, you will then be prompted to type the name."
		      "If you press the TAB key you will autocomplete "
		      "or get a list of possible entries.") ;
	graphRegister (MESSAGE_DESTROY, sprdAddKeyByName) ;
	return ;
      case 'i':
        if (messPrompt("New value ", messprintf("%lld",(long long int)su->z),"i"))
	  freeint(&(u->i)) ;
	touched = TRUE ;
	break ;
      case 'f':
        if (messPrompt("New value ", messprintf("%lg",su->z),"f"))
	  freefloat(&(u->f)) ;
	touched = TRUE ;
	break ;
      case 't':
        if (messPrompt("New value ",
		       stackExists(c->text) && u->i ? stackText(c->text, u->i) : "" , 
		       "t"))
	  { u->i = stackMark(c->text) ;
	    pushText(c->text, freepos()) ;
	  }
	touched = TRUE ;
	break ;
      default:
	break ;
      }
  else
    if (arr(spread->activeLine, spread->activeColonne, SPCELL).parent)
      display (arr(spread->activeLine, spread->activeColonne, SPCELL).parent, 
	       arr(spread->activeLine, spread->activeColonne, SPCELL).grandParent
	       , /* spread->activeColonne ? (char*)0 : mieg: 2005_10_14: always use tree display */ 
	       TREE) ;
  if(touched)
    { 
      spreadReorder(spread) ;
      spreadDisplayData(spread) ;
    }  

  return;
} /* spreadFollow */

static void spreadKeyBoard (int k)
{ return ;
}

static void spreadPick (int box, double x, double y)
{
  SPREAD spread = currentSpread("spreadPick") ;
  
  if (box == 0)			/* picked background */
    return ;

  if (box > spread->tableauBox) 
    {
      int n = box - spread->tableauBox - 1 ,
	max = spread->numberDisplayedColonnes ;
      if (isAddingKey)
	{ 
	  int j, j1, maxCol = arrayMax(spread->colonnes) ;
	  COL *c = 0 ;
	  char *cp ;
	  Array aa = 0 ;
	  if (n >= max) {/* mhmp 20.11.02 */
	    n = n - max ;
	    j1 = n / max ;
	    for(j = 0 , cp = arrp(spread->flags, 0, char) ; 
		j < arrayMax(spread->tableau); cp++, j++)
	      if (!*cp)
		{ 
		  if (!j1--)
		    { 
		      aa = arr(spread->tableau, j, Array) ;
		      break ;
		    }
		}
	    n = n%max ;
	    for(j = 0 ; j < maxCol; j++)
	      { 
		j1 = arr(pos2col,j, int) ;
		c = arrp(spread->colonnes,j1,COL) ;
		if (!c->hidden && !n--)
		  break ;
	      }
	    
	    if (ace_lower(c->type) == 'k' || c->type == 'n')
	      sprdAddKey ((arr(aa, j1, SPCELL).u).k) ;
	  }
	  return ;
	}
      
      
      if (box == spread->activeBox) /* a second hit - follow it */
	{
	  if (n >= max)
	    spreadFollow (spread, box) ;
	}
      else
	spreadSelect (spread, box) ;
    }
  else switch (box)		/* a control box */
    {
    default: break ;
    }
  
  return;
} /* spreadPick */

/******************************************************/


/******************************************************/
/******************************************************/
 
