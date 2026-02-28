#include "../wac/ac.h"
#include "mytime.h"

#include "../whooks/systags.h"
#include "../whooks/tags.h"
#include "freeout.h"
#include "query.h"

extern AC_OBJ ac_key2obj (AC_DB db, KEY key, BOOL fillIt, AC_HANDLE handle) ;
typedef struct{ KEY gene, map ; int a1, a2; BOOL isUp ; } GENE ;
typedef struct{ KEY g1, g2 ; int a1, a2, b1, b2, dx ; BOOL isUp ; } GENE2 ;
static AC_DB DB = 0 ;
static vTXT txt = 0 ;

/***************************************************************/

static void swapA1 (Array genes)
{
  int j1 = arrayMax(genes) ;
  GENE *hh ;

  while (j1--)
    {
      hh = arrayp(genes, j1, GENE) ;
      if (hh->a1 > hh->a2)
	{ 
	  int tmp = hh->a1 ; 
	  hh->a1 = hh->a2 ; hh->a2 = tmp ; 
	  hh->isUp = ! hh->isUp ; 
	}
    }
}

/***************************************************************/

static int geneOrder (const void *a, const void *b) 
{
  const GENE *ga = (const GENE *)a, *gb = (const GENE *)b ;

  if (ga->map != gb->map)
    return ga->map - gb->map ;

  return 
    ga->a1 != gb->a1 ?
    ga->a1 - gb->a1 :
    ga->a2 - gb->a2 ;
}

/***************************************************************/

static int gene2Order (const void *a, const void *b) 
{
  const GENE2 *ga = (const GENE2 *)a, *gb = (const GENE2 *)b ;

  if (ga->g1 != gb->g1)
    return ga->g1 < gb->g1 ? -1 : 1 ;

  return 
    - ga->dx + gb->dx ;     /* longest overlap first */
}

/***************************************************************/

static void crossOverlaps (KEY g1, KEY g2)
{
  vtxtPrintf (txt,"Transcribed_gene \"%s\"\nOverlaps \"%s\"\n\n"
	     , ac_key_name(g1), ac_key_name(g2)
	     ) ;
} /* crossOverlaps */

/******/

static void crossTo_be_fused_with (KEY g1, KEY g2)
{
  vtxtPrintf (txt,"Transcribed_gene \"%s\"\nTo_be_fused_with \"%s\"\n\n"
	    , ac_key_name(g1), ac_key_name(g2)
	    ) ;
} /* crossTo_be_fused_with */

/******/

static void crossAntisens (KEY g1, KEY g2, int dx, BOOL isCoding)
{
  char *tag = isCoding ? "Coding_Antisens_to" : "Antisens_to" ;

  if (g1 == g2)
    return ;
  vtxtPrintf (txt,"Transcribed_gene \"%s\"\n%s \"%s\" %d\n\n"
	    , tag, ac_key_name(g1), ac_key_name(g2), dx
	    ) ;
  vtxtPrintf (txt,"Transcribed_gene \"%s\"\n%s \"%s\" %d\n\n"
	    , tag, ac_key_name(g2), ac_key_name(g1), dx
	    ) ;
} /* crossAntisens */

/******/

static void crossOperon (KEY g1, KEY g2, int dx)
{
  char *mm = "Possible_operon" ;

  if (!g1 || !g2 || g1 == g2)
    return ;
  vtxtPrintf (txt,"Transcribed_gene \"%s\"\n%s \"%s\" %d\n\n"
	    , ac_key_name(g1), mm, ac_key_name(g2), dx
	    ) ;
  vtxtPrintf (txt,"Transcribed_gene \"%s\"\n%s \"%s\" %d\n\n"
	    , ac_key_name(g2), mm, ac_key_name(g2), -dx
	    ) ;
}

/***************************************************************/
/***************************************************************/

typedef struct exonStruct { int x1, x2 ; KEY tag ; } EXON ;

static int exonOrder (const void *va, const void *vb)
{
  const EXON *up = (const EXON *)va, *vp = (const EXON *)vb ;

  if (up->x1 != vp->x1)
    return up->x1 - vp->x1 ;
  else
    return up->x2 - vp->x2 ;
} /* exonOrder */

/******/

static void exonFuse (Array exons) 
{
  int i, j, y1, y2  ;
  EXON *up, *vp ;

  /* fuse exons with common part , happens in case of alternative exons of a tg */
  for (i = 0, up = arrp (exons, i, EXON) ; i < arrayMax(exons) ; i++, up++)
    { 
      if (!up->tag) continue ;
      for (j = i + 1, vp = arrp (exons, j, EXON) ; j < arrayMax(exons) ; j++, vp++)
	{ 
	  if (!vp->tag) continue ;
	  y1 = up->x1 > vp->x1 ? up->x1 : vp->x1 ;
	  y2 = up->x2 < vp->x2 ? up->x2 : vp->x2 ;
	  if (y1 <= y2) /* fuse */
	    { vp->tag = 0 ; if (up->x2 < vp->x2) up->x2 = vp->x2 ; }
	}
    }
  
  /* keep happy few */
  for (i = j = 0, up = vp = arrp (exons, i, EXON) ; i < arrayMax(exons) ; i++, up++)
    { 
      if (up->tag)
	{
	  if (up != vp) *vp = *up ;
	  vp++ ; j++ ;
	}
    }
  arrayMax (exons) = j ;
} /* exonFuse */

/***************************************************************/
/* type 0: no exons, 1: tg */
static void checkStop (GENE *gh, GENE *hh, int type1, int type2, vTXT txt)
{
  KEY g1, g2, product ;
  const char *mm, *isExon ;
  AC_OBJ G1 = 0, G2 = 0, Product = 0 ;
  BOOL isUp1 ;
  int i, j, dx, tmp, a1, a2, b1, b2, x1, x2, z1, z2, tag1=2, tag2=2 ;
  int u1, u2 ;
  AC_TABLE units = 0 ;
  Array exons1 = 0, exons2 = 0 ;
  EXON *up, *vp ;
  int firstMet = 0 ;
  char pId[2], pId5[3], pId3[3] ;
  AC_HANDLE h = ac_new_handle() ;
  
  pId[0]=0; pId[1]=0;
  pId5[0]='5' ; pId5[1]=0; pId5[2]=0;
  pId3[0]='3' ; pId3[1]=0; pId3[2]=0;
  
  z1 = z2 = 0 ;
  g1 = gh->gene ; a1 = gh->a1 ; a2 = gh->a2 ; isUp1 = gh->isUp ;
  g2 = hh->gene ; b1 = hh->a1 ; b2 = hh->a2 ;
  
  if (a1 > a2)
    return ;
  
  x1 = a1 > b1 ? a1 : b1 ; /* x1 x2 = the intersect */
  x2 = a2 < b2 ? a2 : b2 ;
  if (x1 > x2)
    return ;
  
  switch (type1)
    {
    case 1: /* transcribed gene */ 
    case 2: /* predicted gene */
    case 3: /* mrna */
    case 4: /* coding mrna */
      G1 = ac_key2obj (DB, g1, TRUE, h) ;
      switch (type1)
	{
	case 1:	mm = "Splicing" ; tag1 = 2 ; break ;
	case 2: mm = "Source_exons" ; break ;
	case 3:	mm = "Splicing" ; tag1 = 4 ; break ;
	case 4:	
	  mm = "Coding" ; tag1 = 4 ; 
	  units = ac_tag_table (G1, "Product", h) ;
	  for (i = 0; i < units->rows ; i++)
	    {
	      product = ac_table_key (units, i, 0, 0) ;
	      if (keyFindTag (product, str2tag ("Best_product")))
		pId[0] = pId3[1] = pId5[1] = 'A'+ i ;
	      break ;
	    }
	  break ;
	}
      exons1 = arrayHandleCreate (60, EXON, h) ;
      units = ac_tag_table (G1, mm, h) ;
      for (i = j = 0; i < units->rows ; i++)
	{
	  u1 = ac_table_int (units, i, 0, 0) ;
	  u2 = ac_table_int (units, i, 1, 0) ;
	  z1 += u2 - u1 + 1 ;
	  if (!isUp1)
	    {
	      u1 += a1 - 1 ;
	      u2 += a1 - 1 ;
	    }
	  else
	    {
	      tmp = u1 ;
	      u1 = a2 - u2 + 1 ;
	      u2 = a2 - tmp + 1 ;
	    }
	  if (!type2)
	    {
	      if (u1 < x1) u1 = x1 ;
	      if (u2 > x2) u2 = x2 ;
	    } 
	  isExon = ac_table_printable (units, i, tag1, "toto") ;
	  if (type1 == 2 || 
	      ((type1 == 1 || type1 == 3) && strstr (isExon, "xon")) ||
	      (
	       type1 == 4 && 
	       strstr (isExon, pId) && 
	       !strstr (isExon, pId3) &&
	       !strstr (isExon, pId5)
	       ))
	    { 
	      if (j>0 && up->x2 ==  u1 - 1)
		up->x2 = u2 ;
	      else
		{
		  up = arrayp (exons1, j++, EXON) ;
		  up->tag = 1 ; up->x1 = u1 ; up->x2 = u2 ;
		}
	    }
	}
      break ;
    }
  switch (type2)
    {
    case 1: /* transcribed gene */ 
    case 2: /* predicted gene */
    case 3: /* mrna */
    case 4: /* coding mrna */
      G2 = ac_key2obj (DB, g2, TRUE, h) ;
      switch (type2)
	{
	case 1:	mm = "Splicing" ; tag2 = 2 ; break ;
	case 2: mm = "Source_exons" ; break ;
	case 3:	mm = "Splicing" ; tag2 = 4 ; break ;
	case 4:	
	  mm = "Coding" ; tag2 = 4 ; 
	  units = ac_tag_table (G2, "Product", h) ;
	  for (i = 0; units && i < units->rows ; i++)
	    {
	      product = ac_table_key (units, i, 0, 0) ;
	      if (keyFindTag (product, str2tag ("Best_product")))
		pId[0] = pId3[1] = pId5[1] = 'A'+i;
	      if (ac_has_tag (Product, "Best_product"))
		{
		  pId[0] = pId3[1] = pId5[1] = 'A'+ i ;
		  /*  p1 = ac_table_int (units, i, 1, 0) ;
		   * p2 = ac_table_int (units, i, 2, 0) ;  
		   * product in mrna coord 
		   */
		  if ((Product = ac_table_obj (units, i, 0, h)) &&
		      (
		       (firstMet = ac_tag_int (Product, "First_Met", 0))|| 
		       (firstMet = ac_tag_int (Product, "First_ATG", 0))
		       ) &&
		      firstMet > 1
		      )
		    firstMet = 3 * (firstMet - 1) ;
		  else
		    firstMet = 0 ;
		}
	      break ;
	    }
	  break ;
	}
      exons2 = arrayHandleCreate (60, EXON, h) ;
      units = ac_tag_table (G2, mm, h) ;
      for (i = j = 0; units && i < units->rows ; i++)
	{
	  u1 = ac_table_int (units, i, 0, 0) ;
	  u2 = ac_table_int (units, i, 1, 0) ;
	  z2 += u2 - u1 + 1 ;
	  if (!isUp1)
	    {
	      u1 += b1 - 1 ;
	      u2 += b1 - 1 ;
	    }
	  else
	    {
	      tmp = u1 ;
	      u1 = b2 - u2 + 1 ;
	      u2 = b2 - tmp + 1 ;
	    }
	  if (!type2)
	    {
	      if (u1 < x1) u1 = x1 ;
	      if (u2 > x2) u2 = x2 ;
	    } 
	  isExon = ac_table_printable (units, i, tag2, "toto") ;
	  if (type2 == 2 || 
	      ((type2 == 1 || type2 == 3) && strstr (isExon, "xon")) ||
	      (
	       type2 == 4 && 
	       strstr (isExon, pId) && 
	       !strstr (isExon, pId3) &&
	       !strstr (isExon, pId5)
	       ))
	    { 
	      if (j>0 && up->x2 ==  u1 - 1)
		up->x2 = u2 ;
	      else
		{
		  up = arrayp (exons2, j++, EXON) ;
		  up->tag = 1 ; up->x1 = u1 ; up->x2 = u2 ;
		}
	    }
	}
      break ;
    }
  
  if (exons1) { arraySort(exons1, exonOrder) ; exonFuse (exons1) ; }
  if (exons2) { arraySort(exons2, exonOrder) ; exonFuse (exons2) ; }
  
  if (type2 == 4 && exons2 && arrayMax(exons2) && 
      ac_has_tag (G2, "gap_length"))
    /*  hasGap2 = TRUE ; */
  
  dx = 0 ;
  if (exons1 && exons2 && arrayMax(exons1) && arrayMax(exons2))
    {
      /* shift to the pg stop */
      i = arrayMax(exons1) - 1;
      up = arrp (exons1, i, EXON) ;
      /* relocate x1,x2 to the actual exons */
      x2 = up->x2 ;
      
      for (j = dx = 0 ; j < arrayMax(exons2) ; j++)
	{ 
	  vp = arrp (exons2, j, EXON) ;
	  if (vp->x1 <= x2 && vp->x2 >= up->x2)
	    {
	      dx += x2 - vp->x1 + 1 ;
	      vtxtPrintf(txt, "mRNA %s\n", ac_protect (ac_name (G2), h)) ;
	      vtxtPrintf (txt, "Stop_of  %s %d\n\n" 
			  , ac_protect (ac_name (G1), h)
			  , dx ) ;
	      break ;
	    }
	  else
	    dx += vp->x2 - vp->x1 + 1 ;
	}
    }
  ac_free (h) ;
} /* checkStop */

/***************************************************************/
/* type 0: no exons, 1: tg */
static int checkExonOverlap (GENE *gh, GENE *hh, int type1, int type2, int *typep,
			     int *z1p, int *z2p, int *exonSamep, int *exonDiffp)
{
  KEY g1, g2, product ;
  const char *mm, *isExon ;
  AC_OBJ G1 = 0, G2 = 0, Product = 0 ;
  BOOL isUp1 ;
  int type = 0, i, j, dx, ddx, dy, dy1, dy2, tmp, a1, a2, b1, b2, x1, x2, y1, y2, tag1=2, tag2=2 ;
  int u1, u2 ;
  AC_TABLE units = 0 ;
  Array exons1 = 0, exons2 = 0 ;
  EXON *up, *vp ;
  int firstMet = 0 ;
  int exonSame = 0, exonDiff = 0 ;
  char pId[2], pId5[3], pId3[3] ;
  AC_HANDLE h = ac_new_handle() ;

  pId[0]=0; pId[1]=0;
  pId5[0]='5' ; pId5[1]=0; pId5[2]=0;
  pId3[0]='3' ; pId3[1]=0; pId3[2]=0;
  
  g1 = gh->gene ; a1 = gh->a1 ; a2 = gh->a2 ; isUp1 = gh->isUp ; /* map1 = gh->map ; */
  g2 = hh->gene ; b1 = hh->a1 ; b2 = hh->a2 ; /* isUp2 = hh->isUp ;  map2 = hh->map ; */
  
  *z1p = *z2p = 0 ;
  if (a1 > a2)
    return 0 ;

  x1 = a1 > b1 ? a1 : b1 ; /* x1 x2 = the intersect */
  x2 = a2 < b2 ? a2 : b2 ;
  if (x1 > x2)
    return 0 ;

  switch (type1)
    {
    case 1: /* transcribed gene */ 
    case 2: /* predicted gene */
    case 3: /* mrna */
    case 4: /* coding mrna */
      G1 = ac_key2obj (DB, g1, TRUE, h) ;
      switch (type1)
	{
	case 1:	mm = "Splicing" ; tag1 = 2 ; break ;
	case 2: mm = "Source_exons" ; break ;
	case 3:	mm = "Splicing" ; tag1 = 4 ; break ;
	case 4:	
	  mm = "Coding" ; tag1 = 4 ; 
	  units = ac_tag_table (G1, "Product", h) ;
	  for (i = 0; i < units->rows ; i++)
	    {
	      product = ac_table_key (units, i, 0, 0) ;
	      if (keyFindTag (product, str2tag ("Best_product")))
		pId[0] = pId3[1] = pId5[1] = 'A'+ i ;
	      break ;
	    }
	  break ;
	}
      exons1 = arrayHandleCreate (60, EXON, h) ;
      units = ac_tag_table (G1, mm, h) ;
      for (i = j = 0; i < units->rows ; i++)
	{
	  u1 = ac_table_int (units, i, 0, 0) ;
	  u2 = ac_table_int (units, i, 1, 0) ;
	  *z1p += u2 - u1 + 1 ;
	  if (!isUp1)
	    {
	      u1 += a1 - 1 ;
	      u2 += a1 - 1 ;
	    }
	  else
	    {
	      tmp = u1 ;
	      u1 = a2 - u2 + 1 ;
	      u2 = a2 - tmp + 1 ;
	    }
	  if (!type2)
	    {
	      if (u1 < x1) u1 = x1 ;
	      if (u2 > x2) u2 = x2 ;
	    } 
	  isExon = ac_table_printable (units, i, tag1, "toto") ;
	  if (type1 == 2 || 
	      ((type1 == 1 || type1 == 3) && strstr (isExon, "xon")) ||
	      (
	       type1 == 4 && 
	       strstr (isExon, pId) && 
	       !strstr (isExon, pId3) &&
	       !strstr (isExon, pId5)
	       ))
	    { 
	      if (j>0 && up->x2 ==  u1 - 1)
		up->x2 = u2 ;
	      else
		{
		  up = arrayp (exons1, j++, EXON) ;
		  up->tag = 1 ; up->x1 = u1 ; up->x2 = u2 ;
		}
	    }
	}
      break ;
    }
  switch (type2)
    {
    case 1: /* transcribed gene */ 
    case 2: /* predicted gene */
    case 3: /* mrna */
    case 4: /* coding mrna */
      G2 = ac_key2obj (DB, g2, TRUE, h) ;
      switch (type2)
	{
	case 1:	mm = "Splicing" ; tag2 = 2 ; break ;
	case 2: mm = "Source_exons" ; break ;
	case 3:	mm = "Splicing" ; tag2 = 4 ; break ;
	case 4:	
	  mm = "Coding" ; tag2 = 4 ; 
	  units = ac_tag_table (G2, "Product", h) ;
	  for (i = 0; units && i < units->rows ; i++)
	    {
	      product = ac_table_key (units, i, 0, 0) ;
	      if (keyFindTag (product, str2tag ("Best_product")))
		pId[0] = pId3[1] = pId5[1] = 'A'+i;
	      if (ac_has_tag (Product, "Best_product"))
		{
		  pId[0] = pId3[1] = pId5[1] = 'A'+ i ;
		  /* p1 = ac_table_int (units, i, 1, 0) ; */
		  /*  p2 = ac_table_int (units, i, 2, 0) ;   product in mrna coord */
		  if ((Product = ac_table_obj (units, i, 0, h)) &&
		      (
		       (firstMet = ac_tag_int (Product, "First_Met", 0))|| 
		       (firstMet = ac_tag_int (Product, "First_ATG", 0))
		       ) &&
		      firstMet > 1
		      )
		    firstMet = 3 * (firstMet - 1) ;
		  else
		    firstMet = 0 ;
		}
	      break ;
	    }
	  break ;
	}
      exons2 = arrayHandleCreate (60, EXON, h) ;
      units = ac_tag_table (G2, mm, h) ;
      for (i = j = 0; units && i < units->rows ; i++)
	{
	  u1 = ac_table_int (units, i, 0, 0) ;
	  u2 = ac_table_int (units, i, 1, 0) ;
	  *z2p += u2 - u1 + 1 ;
	  if (!isUp1)
	    {
	      u1 += b1 - 1 ;
	      u2 += b1 - 1 ;
	    }
	  else
	    {
	      tmp = u1 ;
	      u1 = b2 - u2 + 1 ;
	      u2 = b2 - tmp + 1 ;
	    }
	  if (!type2)
	    {
	      if (u1 < x1) u1 = x1 ;
	      if (u2 > x2) u2 = x2 ;
	    } 
	  isExon = ac_table_printable (units, i, tag2, "toto") ;
	  if (type2 == 2 || 
	      ((type2 == 1 || type2 == 3) && strstr (isExon, "xon")) ||
	      (
	       type2 == 4 && 
	       strstr (isExon, pId) && 
	       !strstr (isExon, pId3) &&
	       !strstr (isExon, pId5)
	       ))
	    { 
	      if (j>0 && up->x2 ==  u1 - 1)
		up->x2 = u2 ;
	      else
		{
		  up = arrayp (exons2, j++, EXON) ;
		  up->tag = 1 ; up->x1 = u1 ; up->x2 = u2 ;
		}
	    }
	}
      break ;
    }

  if (exons1) { arraySort(exons1, exonOrder) ; exonFuse (exons1) ; }
  if (exons2) { arraySort(exons2, exonOrder) ; exonFuse (exons2) ; }
  
  if (type2 == 4 && exons2 && arrayMax(exons2) && 
      ac_has_tag (G2, "gap_length"))
  
  dx = 0 ;
  if (exons1 && exons2)
    {
      /* shift to the Met */
      if (!isUp1 && firstMet)
	{ vp = arrp (exons2, 0, EXON) ; vp->x1 += firstMet ; }
      if (isUp1 && firstMet)
	{ vp = arrp (exons2, arrayMax(exons2) - 1, EXON) ; vp->x2 -= firstMet ; }
      /* relocate x1,x2 to the actual exons */
      up = arrp (exons1, 0, EXON) ;
      vp = arrp (exons2, 0, EXON) ;
      x1 = up->x1 > vp->x1 ? up->x1 : vp->x1 ;
      up = arrp (exons1, arrayMax(exons1) - 1, EXON) ;
      vp = arrp (exons2, arrayMax(exons2) - 1, EXON) ;
      x2 = up->x2 < vp->x2 ? up->x2 : vp->x2 ;
      type = 0 ;
      /* look who is longer at both ends */
      for (i = 0; i < arrayMax(exons1) ; i++)
	{ 
	  up = arrp (exons1, i, EXON) ;
	  ddx = 0 ;
	  if (up->x2 >= up->x1)
	    {
	      for (j = 0; j < arrayMax(exons2) ; j++)
		{ 
		  vp = arrp (exons2, j, EXON) ;
		  if (vp->x2 >= vp->x1)
		    {
		      y1 = up->x1 > vp->x1 ?  up->x1 : vp->x1 ;
		      y2 = up->x2 < vp->x2 ?  up->x2 : vp->x2 ;
		      if (y1 <= y2)
			ddx += y2 - y1 + 1 ;
		    }
		  /* assymetric because genefinder always supposed to be complete */
		  if (!isUp1 && i == 0 && j == 0)
		    {
		      if (up->x1 ==  vp->x1) type |= 0x20 ;
		      else if (up->x1 >  vp->x1) type |= 0x4 ;  /* shorter */
		      else if (up->x1 <  vp->x1) type |= 0x8 ;  /* longer */
		    }
		  if (isUp1 && i == 0 && j == 0)
		    {
		      if (up->x1 ==  vp->x1) type |= 0x10 ;
		      else if (up->x1 >  vp->x1) type |= 0x1 ;  /* shorter */
		      else if (up->x1 <  vp->x1) type |= 0x2 ;  /* longer */
		    }
		  if (!isUp1 && i == arrayMax(exons1) - 1 && j == arrayMax(exons2) - 1)
		    {
		      if (up->x2 ==  vp->x2) type |= 0x10 ;
		      else if (up->x2 <  vp->x2) type |= 0x1 ; /* shorter */
		      else if (up->x2 >  vp->x2) type |= 0x2 ; /* longer */
		    }
		  if (isUp1 && i == arrayMax(exons1) - 1 && j == arrayMax(exons2) - 1)
		    {
		      if (up->x2 ==  vp->x2 - firstMet) type |= 0x20 ;
		      else if (up->x2 <  vp->x2) type |= 0x4 ; /* shorter */
		      else if (up->x2 >  vp->x2) type |= 0x8 ; /* longer */
		    }
		}
	    }
	  dx += ddx ;
	}
      /* look for a difference in overlapping region */
      for (dy = i = j = 0; type && i < arrayMax(exons1) ; i++)
	{ 
	  up = arrp (exons1, i, EXON) ;
	  ddx = 0 ;
	  if (up->x1 < x1) up->x1 = x1 ;
	  if (up->x2 > x2) up->x2 = x2 ;
	  
	  if (up->x2 >= up->x1)
	    {
	      for (j = 0; j < arrayMax(exons2) ; j++)
		{ 
		  vp = arrp (exons2, j, EXON) ;
		  if (vp->x1 < x1) vp->x1 = x1 ;
		  if (vp->x2 > x2) vp->x2 = x2 ;
		  if (vp->x2 >= vp->x1)
		    {
		      y1 = up->x1 > vp->x1 ?  up->x1 : vp->x1 ;
		      y2 = up->x2 < vp->x2 ?  up->x2 : vp->x2 ;
		      if (y1 <= y2)
			{
			  ddx = y2 - y1 + 1 ;
			  if (i > 0 && j > 0 && 
			      i < arrayMax(exons1) - 1 && 
			      j < arrayMax(exons2) - 1
			      )
			    {
			      if (ddx != up->x2 - up->x1 + 1 ||
				  ddx != vp->x2 - vp->x1 + 1)
				{
				  if (y1 > x1 && y2 < x2) exonDiff++ ;
				}
			      else
				{
				  if (y1 > x1 && y2 < x2) exonSame++ ;
				}
			    }
			  dy += ddx ;
			}
		    }
		}
	    }
	}
      /* look for an exon above or below known end  */
      for (dy1 = j = 0 ; j < arrayMax(exons1) ; j++)
	{ 
	  up = arrp (exons1, j, EXON) ;
	  if (up->x2 >= up->x1)
	    {
	      y1 = up->x1  > x1 ? up->x1 : x1 ;
	      y2 = up->x2  < x2 ? up->x2 : x2 ;
	      if (y1 <= y2)
		dy1 += y2 - y1 + 1 ;
	    }
	}
      for (dy2 = j = 0 ; j < arrayMax(exons2) ; j++)
	{ 
	  vp = arrp (exons2, j, EXON) ;
	   if (vp->x2 > vp->x1)
	    {
	      y1 = vp->x1 > x1 ? vp->x1 : x1 ;
	      y2 = vp->x2 < x2 ? vp->x2 : x2 ;
	      if (y1 <= y2)
		dy2 += y2 - y1 + 1 ;
	    }
	}
      if (dy == dy1 && dy == dy2) /* perfect agreement in the intersect */
	{
	  if ((type & 0x30) == 0x30) /* agree at both ends */
	    type |= 0x800 ;  /* identical */
	  else if ((type & 0x28) && (type & 0x12)) /* longer at both ends */
	    type |= 0x400 ;   /* includes */
	  else if ((type & 0x24) && (type & 0x11)) /* shorter at both ends */
	    type |= 0x200 ;  /* included */
	}
      else
	type |= 0x100 ;
    }
  else if (exons1)
    {
      for (i = 0; i < arrayMax(exons1) ; i++)
	{ 
	  up = arrp (exons1, i, EXON) ;
	  if (up->x2 > up->x1)
	    dx += up->x2 - up->x1 + 1 ;
	}
    }
  else if (exons2)
    { 
      for (i = 0; i < arrayMax(exons2) ; i++)
	{ 
	  up = arrp (exons2, i, EXON) ; 
	  if (up->x2 > up->x1)
	    dx += up->x2 - up->x1 + 1 ;
	}
    }
  else
    {
      dx = x2 - x1 + 1 ;
    }

  if (typep) *typep = type ;
  if (exonSamep) *exonSamep = exonSame ;
  if (exonDiffp) *exonDiffp = exonDiff ;

  ac_free (h) ;
  return dx ;
}

/***************************************************************/

static void searchOverlappingTr (Array trGenes)
{
  GENE *hh1, *hh2 ;
  KEY g1, g2, map1, map2 ;
  int j1, j2, dx, a1, a2, b1, b2, z1, z2 ;

  for (j1 = 0 ; j1 < arrayMax(trGenes) - 1 ; j1++)
    {
      hh1 = arrayp(trGenes, j1, GENE) ;
      g1 = hh1->gene ; a1 = hh1->a1 ; a2 = hh1->a2 ; map1 = hh1->map ;

      for (j2 = j1 + 1 ; j2 < arrayMax(trGenes) ; j2++)
	{
	  hh2 = arrayp(trGenes, j2, GENE) ;
	  g2 = hh2->gene ;  b1 = hh2->a1 ; b2 = hh2->a2 ; map2 = hh2->map ;
	  if (map1 != map2) break ; 

	  dx = ( b2 < a2 ? b2 : a2) -  (b1 > a1 ? b1 : a1) ;	
	  /* da = a2 - a1 ; db = b2 - b1 ; */
 
	  if (dx > 0 && hh1->isUp != hh2->isUp)
	    {  
	      dx = checkExonOverlap (hh1, hh2, 4, 4, 0, &z1, &z2, 0, 0) ;
	      if (dx > 20)
		crossAntisens (g1, g2, dx, TRUE) ;
	    }
	}
    }
}

/***************************************************************/

static void searchOverlappingTG (Array tGenes)
{
  GENE *hh1, *hh2 ;
  KEY g1, g2, map1, map2 ;
  int j1, j2, dx, da, db, a1, a2, b1, b2, z1, z2 ;
  BOOL  firstDownInCis ;

  for (j1 = 0 ; j1 < arrayMax(tGenes) - 1 ; j1++)
    {
      hh1 = arrayp(tGenes, j1, GENE) ;
      g1 = hh1->gene ; a1 = hh1->a1 ; a2 = hh1->a2 ; map1 = hh1->map ;
      firstDownInCis = TRUE ;

      for (j2 = j1 + 1 ; j2 < arrayMax(tGenes) ; j2++)
	{
	  hh2 = arrayp(tGenes, j2, GENE) ;
	  g2 = hh2->gene ;  b1 = hh2->a1 ; b2 = hh2->a2 ; map2 = hh2->map ;
	  if (map1 != map2) break ;
	  dx = ( b2 < a2 ? b2 : a2) -  (b1 > a1 ? b1 : a1) ;	
	  da = a2 - a1 ; db = b2 - b1 ;
 
	  if (hh1->isUp == hh2->isUp &&
	      firstDownInCis &&
	      b1 > a2)
	    {
	      if (!hh1->isUp)
		crossOperon (g1, g2, b1 - a2) ;
	      else
		crossOperon (g2, g1, b1 - a2) ;
	      firstDownInCis = FALSE ;
	      break ;
	    } 

	  if (hh1->isUp == hh2->isUp && (3*dx > da || 3*dx > db))
	    {	
	      crossOverlaps (g1, g2) ;
	      dx = checkExonOverlap (hh1, hh2, 1, 1, 0, &z1, &z2, 0, 0) ;
	      if (dx > 20)
		{
		  if (da > db || (da == db  && g1 < g2))
		    crossTo_be_fused_with (g1, g2) ;
		  else
		    crossTo_be_fused_with (g2, g1) ;
		}
	    }
	      
	  else if (dx > 0 && hh1->isUp != hh2->isUp)
	    {   
	      crossOverlaps (g1, g2) ;
	      dx = checkExonOverlap (hh1, hh2, 1, 1, 0, &z1, &z2, 0, 0) ;
	      if (dx > 20)
		crossAntisens (g1, g2, dx, FALSE) ;
	    }
	  if (!firstDownInCis && b1 > a2)
	    break ;
	}
    }
}

/*********/

static void crossPg2Tg (KEY pg, KEY tg)
{
  vtxtPrintf (txt,"Sequence \"%s\"\nMatching_transcribed_gene \"%s\"\n\n"
	      , ac_key_name(pg), ac_key_name(tg)
	      ) ;
} /* crossPg2Tg */

/*********/
#ifdef JUNK

static void compatibility (GENE *ph, GENE *hh, vTXT txt)
{
  AC_HANDLE h = ac_new_handle () ;
  KEY pg = ph->gene, tg = hh->gene, tprod = 0 ;
  int ii, jj, jj0, a1 = ph->a1, a2 = ph->a2, b1 = hh->a1, b2 = hh->a2 ;
  GENE *PGene = 0, *TGene = 0 ;
  BOOL ok = FALSE ;
  BOOL pUp = ph->isUp ;
  BOOL tUp = hh->isUp; 
  int x, y1, y2, z1, z2 ;
  AC_TABLE punits = 0, tunits = 0 ;

  if (1)
    {
      punits = ac_tag_table (PGene, "Source_Exons", h) ;
      tunits = ac_tag_table (TGene, "Derived_sequence", h) ;

      tprod = 0 ; jj0 = 0 ; jj = -1 ;
      for (ii = 0 ; ii < punits->rows ; ii++)
	{
	  y1 = ac_table_int (punits, ii, 0, 0) ;
	  y2 = ac_table_int (punits, ii, 1, 0) ;

	  if (pUp)
	    { y1 = a2 - y1 + 1 ; y2 = a2 - y2 + 1 ; }
	  else
	    { y1 += a1 - 1 ; y2 += a1 - 1 ; }
	  /* now search that exon in the tg */
	retry:
	  for (ok = FALSE, jj = jj0 ; jj < tunits->rows ; jj++)
	    {
	      z1 = ac_table_int (tunits, ii, 1, 0) ;
	      z2 = ac_table_int (tunits, ii, 2, 0) ;

	      x = ac_table_int (tunits, ii, 4, 0) ;

	      if (x != 1 + ii/2) continue ;
	      if (tprod && tprod != ac_table_key (tunits, ii, 0, 0))
		{ jj0 = jj ; tprod = ac_table_key (tunits, ii, 0, 0) ; goto retry ; }
	      if (tUp)
		{ z1 = b2 - z1 + 1 ; z2 = b2 - z2 + 1 ; }
	      else
		{ z1 += b1 - 1 ; z2 += b1 - 1 ; }
	      if ( z2 - z1 == y2 - y1) 
		/* just the diff because of the float in map model */
		/* z1 == y1 && z2 == y2) */
		{ ok = TRUE ; jj0 = jj ; tprod = ac_table_key (tunits, ii, 0, 0) ; break ; }
	    }
	  if (!ok)
	    goto abort ;
	}
      /* so now all exons of tg match the exons of pg */
      /* i should test the reciprocal */
      if (jj >= 0 &&
	  jj < arrayMax(tunits) && 
	  (ac_table_key (punits, jj, 0, 0) == tprod &&
	  (ac_table_key (punits, jj, 3, 0) == str2tag("exon"))
	goto abort ;
      /* success */  
      vtxtPrintf (txt,"Gene \"%s\"\nIdentical_to_genefinder \"%s\"\n\n"
		  , ac_key_name (tg), ac_key_name (pg) ) ;
    }
  
abort:
  ac_free (h) ;
}
#endif

/*********/

static void matchProduct_genefinder (Array pGenes, Array products)
{
  AC_HANDLE h = ac_new_handle () ;
  GENE *ph, *hh ;
  int j1, jh ;
  Array hits = arrayCreate (12,GENE) ;

  if (! arrayMax (pGenes) || !arrayMax (products))
    return ;
  for (j1 = jh = 0 ; j1 < arrayMax(pGenes) ; j1++)
    {
      ph = arrayp(pGenes, j1, GENE) ;
      hh = arrayp (hits, jh++, GENE) ;
      hh->gene = ph->gene ; 
      if (! ph->isUp) hh->a1 = ph->a2 ;
      else hh->a1 = ph->a1 ;
      hh->a2 = hh->a1 + 1 ; /* dummy */
      hh->isUp = FALSE ;  /* pg */
    }
  for (j1 = 0 ; j1 < arrayMax(products) ; j1++)
    {
      ph = arrayp(products, j1, GENE) ;
      hh = arrayp (hits, jh++, GENE) ;
      hh->gene = ph->gene ; 
      if (! ph->isUp) hh->a1 = ph->a2 ;
      else hh->a1 = ph->a1 ;
      hh->a2 = hh->a1 + 2 ; /* dummy */
      hh->isUp = TRUE ;  /* product */
    }

  arraySort (hits, geneOrder) ;
  for (jh = 0 ; jh < arrayMax(hits) ; jh++)
    {
      ph = arrayp(hits, jh, GENE) ;
      for (j1 = jh-1 ; j1 >= 0 ; j1--) 
	{
	  hh = arrayp(hits, j1, GENE) ;
	  if (hh->a1 < ph->a1)
	    break ;
	}
      for ( ; j1 < arrayMax(hits) ; j1++)
	{
	  hh = arrayp(hits, j1, GENE) ;
	  if (hh == ph)
	    continue ;
	  if (hh->a1 < ph->a1)
	    continue ;
	  if (hh->a1 > ph->a1)
	    break ;
	  if (! ph->isUp  && hh->isUp)
	    vtxtPrintf (txt, "Product %s\nSame_stop_as_model %s\n\n"
			, ac_protect (ac_key_name(hh->gene), h)
			, ac_protect (ac_key_name(ph->gene), h)
			) ;
	  /* vtxtPrintf (txt, "") ; inutile ce truc ? */
	  if (ph->isUp && hh->isUp )
	    vtxtPrintf (txt, "Product %s\nSame_stop_as %s\n\n"
			, ac_protect (ac_key_name(hh->gene), h)
			, ac_protect (ac_key_name(ph->gene), h)
			) ;
	}
    }
  ac_parse (DB, vtxtPtr (txt), 0, 0, h) ;
  ac_free (h) ;
} /* matchProduct_genefinder */

/***************************/

static int  matchTg_genefinder (Array pGenes, Array tGenes)
{
  GENE *ph, *hh ;
  KEY pg, tg, map1, map2 ; BOOL isUp ;
  int j1, j2, jh, dx, a1, a2, b1, b2, z1, z2 ;
  int ok = 0 ;
  Array hits = arrayCreate (12,GENE2) ;
  GENE2 *mm, *mm1 ;

  if (! arrayMax (pGenes) || !arrayMax (tGenes))
    return 0 ;
  for (j1 = jh = 0 ; j1 < arrayMax(pGenes) ; j1++)
    {
      ph = arrayp(pGenes, j1, GENE) ;
      pg = ph->gene ; a1 = ph->a1 ; a2 = ph->a2 ; isUp = ph->isUp ; map1 = ph->map ;
      for (j2 = 0 ; j2 < arrayMax(tGenes) ; j2++)
	{
	  hh = arrayp(tGenes, j2, GENE) ;
	  if (isUp != hh->isUp) continue ;
	  tg = hh->gene ;  b1 = hh->a1 ; b2 = hh->a2 ; map2 = hh->map ;
	  if (map1 < map2)
	    break ;
	  if (map1 != map2)
	    continue ;
	  if (b1 > a2) break ;
	  dx = ( b2 < a2 ? b2 : a2) -  (b1 > a1 ? b1 : a1) ;	
	  /* da = a2 - a1 ; db = b2 - b1 ; */
	  if (dx > 20)
	    {
	      dx = checkExonOverlap (ph, hh, 2, 1, 0, &z1, &z2, 0, 0) ;
	      if  (dx > 20)
		{
		  ok++ ;
		  mm = arrayp (hits, jh++, GENE2) ;
		  mm->g1 = pg ; mm->g2 = tg ;
		  mm->dx = dx ;
		  mm->isUp = isUp ;
		  mm->a1 = ph->a1 ; mm->a2 = ph->a2 ;
		  mm->b1 = hh->a1 ; mm->b2 = hh->a2 ;
		}
	    }
	}
    }
  arraySort (hits, gene2Order) ;
  for (jh = 0, mm1 = 0 ; jh < arrayMax (hits) ; jh++)
    {
      mm = arrp (hits, jh, GENE2) ;
      if (mm1 &&
	  mm->g1 == mm1->g1 &&
	  mm->dx < mm1->dx)
	continue ;
      j1 = jh ; mm1 = mm ;
      crossPg2Tg (mm->g2, mm->g2) ;
    }
  arrayDestroy (hits) ;
  return ok ;
}

/*********/

static void crossPgPg (KEY pg1, KEY pg2)
{
  KEY gene1 = keyGetKey (pg1, str2tag ("Model_of_gene")) ;
  KEY gene2 = keyGetKey (pg2, str2tag ("Model_of_gene")) ;

  if (gene1 && !gene2)
    vtxtPrintf (txt,"Sequence \"%s\"\nModel_of_gene \"%s\"\n\n"
		, ac_key_name (pg2), ac_key_name (gene1)
		) ;
  if (gene2 && !gene1)
    vtxtPrintf (txt,"Sequence \"%s\"\nModel_of_gene \"%s\"\n\n"
		, ac_key_name (pg1), ac_key_name (gene2)
		) ;
} /* crossPgPg */

/*********/

static int  matchGenefinder_genefinder (Array pGenes)
{
  GENE *ph, *hh ;
  KEY pg1, pg2, map1, map2 ; BOOL isUp ;
  int j1, j2, dx,  a1, a2, b1, b2, z1, z2 ;
  int ok = 0 ;

  if (! arrayMax (pGenes))
    return 0 ;
  for (j1 = 0 ; j1 < arrayMax(pGenes) ; j1++)
    {
      ph = arrayp(pGenes, j1, GENE) ;
      pg1 = ph->gene ; a1 = ph->a1 ; a2 = ph->a2 ; isUp = ph->isUp ; map1 = ph->map ;
      for (j2 = j1+1 ; j2 < arrayMax(pGenes) ; j2++)
	{
	  hh = arrayp(pGenes, j2, GENE) ;
	  if (isUp != hh->isUp) continue ;
	  pg2 = hh->gene ;  b1 = hh->a1 ; b2 = hh->a2 ; map2 = hh->map ;
	  if (map1 < map2)
	    break ;
	  if (map1 != map2)
	    continue ;
	  if (b1 > a2) break ;
	  dx = ( b2 < a2 ? b2 : a2) -  (b1 > a1 ? b1 : a1) ;	
	  /* da = a2 - a1 ; db = b2 - b1 ; */
	  if (dx > 20)
	    {
	      dx = checkExonOverlap (ph, hh, 2, 2, 0, &z1, &z2, 0, 0) ;
	      if  (dx > 20)
		{
		  ok++ ;
		  crossPgPg (pg1, pg2) ;
		}
	    }
	}
    }
  return ok ;
}

/***************************************************************/

static void crossPGeneRnai (KEY pg, KEY rnai)
{
  KEY gene = keyGetKey (pg, str2tag ("Model_of_gene")) ;
  
  if (gene)
    vtxtPrintf (txt,"Gene \"%s\"\nRNAi \"%s\"\n\n"
		, ac_key_name (gene), ac_key_name (rnai)
	      ) ;
} /* crossPGeneRnai */

/*********/

static int matchPGene_Rnai (Array pGenes, Array rnaiGenes)
{
  GENE *ph, *hh ;
  KEY pg, rnai, map1, map2 ;
  int j1, j2, dx, a1, a2, b1, b2, z1, z2 ;
  int ok = 0 ;

  if (! arrayMax (pGenes) || !arrayMax (rnaiGenes))
    return 0 ;

  printf ("// Attempting to match %d predicted_genes to %d rnai\n", arrayMax(pGenes), arrayMax(rnaiGenes)) ;
  if (arrayMax(rnaiGenes))
    for (j1 = 0 ; j1 < arrayMax(pGenes) ; j1++)
    {
      ph = arrayp(pGenes, j1, GENE) ;
      pg = ph->gene ; a1 = ph->a1 ; a2 = ph->a2 ; map1 = ph->map ;
      for (j2 = 0 ; j2 < arrayMax(rnaiGenes) ; j2++)
	{
	  hh = arrayp(rnaiGenes, j2, GENE) ;
	  /* if (isUp != hh->isUp) continue ; DOES NOT MATTER fro RNAi */
	  rnai = hh->gene ;  b1 = hh->a1 ; b2 = hh->a2 ; map2 = hh->map ;
	  if (map1 < map2)
	    break ;
	  if (map1 != map2)
	    continue ;
	  if (b1 > a2) break ;
	  dx = ( b2 < a2 ? b2 : a2) -  (b1 > a1 ? b1 : a1) ;	
	  /*  da = a2 - a1 ; db = b2 - b1 ; */
	  if (dx > 25)
	    {
	      dx = checkExonOverlap (ph, hh, 2, 0, 0, &z1, &z2, 0, 0) ;
	      if  (dx > 25)
		{ 
		  ok++ ;
		  crossPGeneRnai (pg, rnai) ;
		}
	    }
	}
    }
  return ok ;
} /* matchPGene_Rnai */

/***************************************************************/

static void crossTGeneRnai (KEY tg, KEY rnai)
{
  KEY gene =  keyGetKey (tg, str2tag ("Gene")) ;

  if (gene && rnai)
    vtxtPrintf (txt,"Gene \"%s\"\nRNAi \"%s\"\n\n"
		, ac_key_name (gene), ac_key_name (rnai)
	      ) ;
} /* crossTGeneRnai */

/*********/

static int matchTGene_Rnai (Array tGenes, Array rnaiGenes)
{
  GENE *ph, *hh ;
  KEY tg, rnai, map1, map2 ; 
  int j1, j2, dx,  a1, a2, b1, b2, z1, z2 ;
  int ok = 0 ;

  printf ("// Attempting to match %d genes to %d rnai\n", arrayMax(tGenes), arrayMax(rnaiGenes)) ;
  if (! arrayMax (tGenes) || !arrayMax (rnaiGenes))
    return 0 ;

  for (j1 = 0 ; j1 < arrayMax(tGenes) ; j1++)
    {
      ph = arrayp(tGenes, j1, GENE) ;
      tg = ph->gene ; a1 = ph->a1 ; a2 = ph->a2 ; map1 = ph->map ;
      for (j2 = 0 ; j2 < arrayMax(rnaiGenes) ; j2++)
	{
	  hh = arrayp(rnaiGenes, j2, GENE) ;
	  /* if (isUp != hh->isUp) continue ; DOES NOT MATTER fro RNAi */
	  rnai = hh->gene ;  b1 = hh->a1 ; b2 = hh->a2 ; map2 = hh->map ;
	  if (map1 < map2)
	    break ;
	  if (map1 != map2)
	    continue ;
	  if (b1 > a2) break ;
	  dx = ( b2 < a2 ? b2 : a2) -  (b1 > a1 ? b1 : a1) ;	
	  /* da = a2 - a1 ; db = b2 - b1 ; */
	  if (dx > 25)
	    {
	      dx = checkExonOverlap (ph, hh, 1, 0, 0, &z1, &z2, 0, 0) ;
	      if  (dx > 25)
		{
		  ok++ ;
		  crossTGeneRnai (tg, rnai) ;
		}
	    }
	}
    }
  return ok ;
} /* matchTGene_Rnai */

/***************************************************************/

static void crossTrRnai (KEY tr, KEY rnai)
{
  KEY tg =  keyGetKey (tr, str2tag ("From_gene")) ;

  if (tg)
    crossTGeneRnai (tg, rnai) ;
} /* crossTrRnai */

/*********/

static int matchTr_Rnai (Array trGenes, Array rnaiGenes)
{
  GENE *ph, *hh ;
  KEY tr, rnai, map1, map2 ; 
  int j1, j2, dx, a1, a2, b1, b2, z1, z2 ;
  int ok = 0 ;

  printf ("// Attempting to match %d genes to %d rnai\n", arrayMax(trGenes), arrayMax(rnaiGenes)) ;
  if (! arrayMax (trGenes) || !arrayMax (rnaiGenes))
    return 0 ;

  for (j1 = 0 ; j1 < arrayMax(trGenes) ; j1++)
    {
      ph = arrayp(trGenes, j1, GENE) ;
      tr = ph->gene ; a1 = ph->a1 ; a2 = ph->a2 ; map1 = ph->map ;
      for (j2 = 0 ; j2 < arrayMax(rnaiGenes) ; j2++)
	{
	  hh = arrayp(rnaiGenes, j2, GENE) ;
	  /* if (isUp != hh->isUp) continue ; DOES NOT MATTER fro RNAi */
	  rnai = hh->gene ;  b1 = hh->a1 ; b2 = hh->a2 ; map2 = hh->map ;
	  if (map1 < map2)
	    break ;
	  if (map1 != map2)
	    continue ;
	  if (b1 > a2) break ;
	  dx = ( b2 < a2 ? b2 : a2) -  (b1 > a1 ? b1 : a1) ;	
	  /* da = a2 - a1 ; db = b2 - b1 ; */
	  if (dx > 25)
	    {
	      dx = checkExonOverlap (ph, hh, 1, 0, 0, &z1, &z2, 0, 0) ;
	      if  (dx > 25)
		{
		  ok++ ;
		  crossTrRnai (tr, rnai) ;
		}
	    }
	}
    }
  return ok ;
} /* matchTr_Rnai */

/***************************************************************/

static void compatibilityTr (GENE *ph, GENE *hh, vTXT txt)
{
  AC_HANDLE h = ac_new_handle () ;
#ifdef JUNK

  int nx, ny, ii, jj, jj0, a1 = ph->a1, a2 = ph->a2, b1 = hh->a1, b2 = hh->a2 ;
  Array tunits = 0 ;
  Array punits = 0 ;
  BOOL pUp = ph->isUp ;
  BOOL tUp = hh->isUp; 
  int y1, y2, z1, z2 ;
  KEY tr = hh->gene, type,  _Exon ;


  _Exon = aceTag ("Exon",0) ;

  tunits = aceArrayCreate (80, AceUnit) ; 
  punits = aceArrayCreate (80, AceUnit) ; 

  TR = aceOpenKey(tr, 0) ; /* do not write lock tr , this would prevent the XREFs when updating pg */
  if (PG && TR)
    {
      aceGetArray (PG, str2tag("Source_Exons"), punits, 2) ;
      aceGetArray (TR, str2tag("Splicing"), tunits, 3) ;

      jj0 = 0 ; jj = -1 ;
      for (nx = 0, ny = 0, ii = 0 ; ii < arrayMax(punits) ; nx++, ii += 2)
	{
	  up = arrp (punits, ii, AceUnit) ;
	  y1 = up[0].i ; y2 = up[1].i ;
	  if (pUp)
	    { y1 = a2 - y1 + 1 ; y2 = a2 - y2 + 1 ; }
	  else
	    { y1 += a1 - 1 ; y2 += a1 - 1 ; }
	  /* now search that exon in the tg */
	  z1 = z2 = -1 ;
	  for (;ny < arrayMax(tunits)/3;)
	  {
	    vp = arrp (tunits, 3*ny, AceUnit) ;
	    ny++; /* increment here because of subsequent break  */
	    z1 = vp[0].i ; z2 = vp[1].i ; type = vp[2].k ;
	    if (type == _Exon)
	      break ;
	  }
	  if (tUp)
	    { z1 = b2 - z1 + 1 ; z2 = b2 - z2 + 1 ; }
	  else
	    { z1 += b1 - 1 ; z2 += b1 - 1 ; }
	  if ((nx > 0 && z1 != y1) || 
	      (nx == 0 && !pUp && z1 > y1) ||
	      (nx == 0 && pUp && z1 < y1) ||
	      z2 != y2)
	    goto abort ;
	}
      for (;ny < arrayMax(tunits)/3; )
	{
	  vp = arrp (tunits, 3*ny, AceUnit) ;
	  ny++; /* increment here because of subsequent break  */
	  z1 = vp[0].i ; z2 = vp[1].i ; type = vp[2].k ;
	  if (type == _Exon)
	    goto abort ;  /* one exon too much */
	}
    
      /* success */
       if (aceAddTag (PG, str2tag("Identical_to_mRNA")))
	 aceAddKey (PG, tr) ;
    }
 abort:
#endif

  ac_free (h) ;
}

/*********/

static void compatibilityPcr (GENE *ph, GENE *hh, vTXT txt)
{
#ifdef JUNK
  int ii, jj, jj0, a1 = ph->a1, a2 = ph->a2, b1 = hh->a1, b2 = hh->a2;
  Array tunits = 0 ;
  Array punits = 0 ;


  AceUnit *vp ;
  AceInstance TR = 0 ;

  BOOL pUp = ph->isUp ;

  int z1, z2, z3, z3max, isPcr5 = 0, isPcr3 = 0 ;
  int pStart, pStop, tStop ;
  KEY tr = hh->gene, type,  _Exon,  _Intron ;

  acePushContext () ;

  _Exon = aceTag ("Exon",0) ;
  _Intron = aceTag ("Intron",0) ;
  tunits = aceArrayCreate (80, AceUnit) ; 
  punits = aceArrayCreate (80, AceUnit) ; 

  jj0 = 0 ; jj = -1 ;
  /*  goto tr coordinates */
  if (pUp)
    {
      pStart = b2 - a2 + 1 ;
      pStop =  b2 - a1 + 1 ; 
      tStop = b2 - b1 + 1 ;
    }
  else
    {
      pStart = a1 - b1 + 1 ;
      pStop = a2 - b1 + 1 ;
      tStop = b2 - b1 + 1 ;
    }

      /* 1:INCOMPLET 2:ABSENT 3:INEXACT 4:EXACT */
      /* now search first exon in the tg */
  TR = aceOpenKey(tr, 0) ; /* do not write lock tr , this would prevent the XREFs when updating pg */
  aceGetArray (TR, str2tag("Splicing"), tunits, 5) ;
  z3max = 0 ;
  for (ii = 0 ; ! z3max && ii < arrayMax(tunits) ; ii += 5)
    {
      vp = arrp (tunits, ii, AceUnit) ;
      z1 = vp[0].i ; z2 = vp[1].i ; type = vp[4].k ; z3  = vp[3].i ;
      if (type == _Exon) z3max = z3 ;
    }
   for (ii = isPcr5 = 0 ; !isPcr5 && ii < arrayMax(tunits) ; ii += 5)
    {
      vp = arrp (tunits, ii, AceUnit) ;
      z1 = vp[0].i ; z2 = vp[1].i ; type = vp[4].k ; z3  = vp[3].i ;
      
      if (z1 == pStart && type == _Exon && z3 == z3max)
	isPcr5 = 4 ;  /* EXACT */
      else if (
	       (z1 == pStart && type != _Intron && z3 != z3max) ||
	       (z1 < pStart && z2 > pStart  && type != _Intron)
	       )
	{
	  if (aceGotoTag (TR, str2tag("Found5p")))
	    isPcr5 = 3 ;  /* INEXACT starts inside */
	  else
	    isPcr5 = 1 ;  /* INCOMPLET tg may start here or elsewhere */
	}
      else if (z1 < pStart && z2 > pStart && z2 < pStart +12 && type == _Intron)
	isPcr5 = 3 ;  /* INEXACT oligo is inside intron */
      else if (z1 < pStart && z2 > pStart +12 && type == _Intron)
	isPcr5 = 2 ;  /* ABSENT  oligo in tg intron */
    }
  if (!isPcr5)   /* outside the loop because of utr parts */
    {
      if (aceGotoTag (TR, str2tag("Found5p")))
	isPcr5 = 2 ;  /* ABSENT tg known to stop upstream */
      else
	isPcr5 = 1 ;  /* INCOMPLET tg may start here or downstream */
      
    }
  z3max = 0 ;
  for (ii = 0 ; ii < arrayMax(tunits) ; ii += 5)
    {
      vp = arrp (tunits, ii, AceUnit) ;
      z1 = vp[0].i ; z2 = vp[1].i ; type = vp[4].k ; z3  = vp[3].i ;
      if (type == _Exon) z3max = z3 ;
    }
  for (ii = isPcr3 = 0 ; !isPcr3 && ii < arrayMax(tunits) ; ii += 5)
    {
      vp = arrp (tunits, ii, AceUnit) ;
      z1 = vp[0].i ; z2 = vp[1].i ; type = vp[4].k ; z3  = vp[3].i ;
      
      if (z2 == pStop && type == _Exon && z3 == z3max)
	isPcr3 = 4 ;  /* exact */
      else if (
	       (z2 == pStop && type != _Intron && z3 != z3max)  || /* i.e. also UTR */
	       /* INEXACT stops on wrong exon */
	       (z1 <= pStop && z2 > pStop && type != _Intron )
	       )
	{
	  if (aceGotoTag (TR, str2tag("Found3p")))
	    isPcr3 = 3 ; /* INEXACT stops inside */
	  else
	    isPcr3 = 1 ;  /* INCOMPLET tg may start here or downstream */
	}
      else if (z2 < pStop  && z2 > pStop - 12 && type == _Intron)
	isPcr3 = 3 ;  /* INEXACT oligo is inside intron */
      else if (z1 < pStop && z1 > pStop - 12 && z2 > pStop && type == _Intron)
	isPcr3 = 3;  /* INEXACT  oligo in tg intron */
      else if (z1 <= pStop - 12 && z2 > pStop && type == _Intron)
	isPcr3 = 2 ;  /* ABSENT  oligo in tg intron */
    }
  if (!isPcr3)   /* outside the loop because of utr parts */
    {
      if (aceGotoTag (TR, str2tag("Found3p")))
	isPcr3 = 2 ;  /* ABSENT tg known to stop upstream */
      else
	isPcr3 = 1 ;  /* INCOMPLET tg may start here or downstream */
      
    }
  acePopContext (TRUE) ;   /* needed BEFORE we edit the obj, which is ridiculous */
  if (isPcr5 || isPcr3)
    {
      aceAddTag (PG, str2tag("Matching_mRNA")) ;
      aceAddKey (PG, tr) ;
      if (acePushType (PG))
	{
	  switch (isPcr5)
	    {
	      case 4: 
		aceAddTag (PG, str2tag("mExact")) ;
		break ;
	      case 3: 
		aceAddTag (PG, str2tag("mInexact")) ;
		break ;
	      case 2: 
		aceAddTag (PG, str2tag("mAbsent")) ;
		break ;
	      case 1: 
		aceAddTag (PG, str2tag("mNoInfo")) ;
		break ;
	    }
	  switch (isPcr3)
	    {
	      case 4: 
		aceAddTag (PG, str2tag("sExact")) ;
		break ;
	      case 3: 
		aceAddTag (PG, str2tag("sInexact")) ;
		break ;
	      case 2: 
		aceAddTag (PG, str2tag("sAbsent")) ;
		break ;
	      case 1: 
		aceAddTag (PG, str2tag("sNoInfo")) ;
		break ;
	    }

	  acePopType (PG) ;
	}
    }
  /* success 
     aceAddTag (PG, str2tag("Matching_transcript"))Identical_to_genefinder")))
     aceAddKey (PG, tr) ;
     pushdobj add Start 1 2 3 4 5
  */
  
#endif
 
}

/*********/

static int matchTr_genefinder (Array pGenes, Array trGenes)
{
  GENE *ph, *hh, *ah ;
  KEY pg, map1, map2 ; 
  BOOL isUp ;
  int j1, j2, dx, dxa, dxb, a1, a2, b1, b2, iaa, z1a, z2a, z1b, z2b ;
  int ok = 0, bestType = 0, exonSame = 0, exonDiff = 0,
    exonSamea = 0, exonDiffa = 0, exonSameb = 0, exonDiffb = 0 ;
  AC_OBJ PG = 0 ;
  static Array aa = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  AC_HANDLE h1 = ac_new_handle () ;
  vTXT txt = vtxtHandleCreate (h1) ;
  vTXT stopTxt = vtxtHandleCreate (h1) ;

  if (! arrayMax (pGenes) || !arrayMax (trGenes))
    return 0 ;

  for (j1 = 0 ; j1 < arrayMax(pGenes) ; j1++)
    {  
      ph = arrayp(pGenes, j1, GENE) ;
      pg = ph->gene ;
      ac_free (h) ; 
      h = ac_new_handle () ;
      PG = ac_key2obj (DB, pg, TRUE, h) ;
      if (! PG)
	continue ;
      vtxtPrintf (txt,"Sequence %s\n-D Matching_mRNA\n"
		  , ac_protect (ac_name (PG), h)
		  ) ;
      bestType = 0 ;

      aa = arrayReCreate (aa, 12, GENE) ;
      iaa = bestType = 0 ;

      a1 = ph->a1 ; a2 = ph->a2 ; isUp = ph->isUp ; map1 = ph->map ;
      for (j2 = 0 ; j2 < arrayMax(trGenes) ; j2++)
	{
	  hh = arrayp(trGenes, j2, GENE) ;
	  if (isUp != hh->isUp) continue ;
	  /* tr = hh->gene ; */  b1 = hh->a1 ; b2 = hh->a2 ; map2 = hh->map ;
	  if (map1 < map2)
	    break ;
	  if (map1 != map2)
	    continue ;
	  if (b1 > a2 + 100) 
	    break ;
	  dx = ( b2 < a2 ? b2 : a2) -  (b1 > a1 ? b1 : a1) ;	
	  /* da = a2 - a1 ; db = b2 - b1 ; */
	  if (dx > 0) /* they touch, we must examine if some exons are common */
	    {
	      int type = 0, typea = 0, typeb = 0 ;
	      checkStop (ph, hh, 2, 3, stopTxt) ;
	      ok++ ; exonSame = exonDiff = 0 ;
	      dxa = checkExonOverlap (ph, hh, 2, 3, &typea, &z1a, &z2a, &exonSamea, &exonDiffa) ;
	      dxb = checkExonOverlap (ph, hh, 2, 4, &typeb, &z1b, &z2b, &exonSameb, &exonDiffb) ;
	      if (
		  (3 *dxa > z1a || 3 * dxa > z2b) ||
		  (3 *dxb > z1b || 3 * dxb > z2b) 
		  )
		{
		  if (typea > typeb)
		    {
		      type = typea ;
		      exonSame = exonSamea ;
		      exonDiff = exonDiffa ;
		    }
		  else
		    {
		      type = typeb ;
		      exonSame = exonSameb ;
		      exonDiff = exonDiffb ;
		    }
		  type = (type << 16) + exonSame ;
		  ah = arrayp (aa, iaa++, GENE) ;
		  ah->gene = hh->gene ;
		  ah->a1 = type ;
		  ah->a2 = exonDiff ;
		  if (bestType < type)
		    bestType = type ;
		}
		
	      compatibilityTr (ph, hh, txt) ; 
	      compatibilityPcr (ph, hh, txt) ;
	    }
	}
      for (iaa = 0 ; iaa < arrayMax (aa) ; iaa++)
	{
	  int type ;
	  char *cp ;

	  ah = arrayp (aa, iaa, GENE) ;
	  if (ah->a1 == bestType ||
	      ((ah->a1 & 0x4000000) && (bestType & 0x4000000))
	      ) ; /* accept several 'includes': a long prediction can match several mrna */
	  else
	    continue ; 
	  type = (ah->a1) >> 16 & 0xffff;
	  exonSame = ah->a1 & 0xffff ;
	  exonDiff = ah->a2 ;
	  cp = messprintf ("Matching_mRNA %s "
			  , ac_protect (ac_key_name(ah->gene), h)
			  ) ;
	  if (type & 0x800) vtxtPrintf (txt, "%s Identical\n", cp) ;
	  if (type & 0x400) vtxtPrintf (txt, "%s Includes\n", cp) ;
	  if (type & 0x200) vtxtPrintf (txt, "%s Included\n", cp) ;
	  if (type & 0x100) vtxtPrintf (txt, "%s Different\n", cp) ;

	  if (type & 0x20) vtxtPrintf (txt, "%s Same_5p\n", cp) ;
	  if (type & 0x10) vtxtPrintf (txt, "%s Same_3p\n", cp) ;

	  if (type & 0x8) vtxtPrintf (txt, "%s Longer_5p\n", cp) ;
	  if (type & 0x4) vtxtPrintf (txt, "%s Shorter_5p\n", cp) ;
	  if (type & 0x2) vtxtPrintf (txt, "%s Longer_3p\n", cp) ;
	  if (type & 0x1) vtxtPrintf (txt, "%s Shorter_3p\n", cp) ;

	  if (exonSame)
	    vtxtPrintf (txt,"Matching_mRNA %s Identical_internal_exon %d\n"
			, ac_protect (ac_key_name(ah->gene), h)
			, exonSame
			) ;
	  if (exonDiff)
	    vtxtPrintf (txt,"Matching_mRNA %s Different_internal_exon %d\n"
			, ac_protect (ac_key_name(ah->gene), h)
			, exonDiff
			) ;
	}
      vtxtPrint (txt, "\n") ;	
    }
  
  if (0) ac_parse (DB, vtxtPtr (txt), 0, 0, h) ;
  if (1) ac_parse (DB, vtxtPtr (stopTxt), 0, 0, h) ;
  
  ac_free (h) ;
  ac_free (h1) ;
  return ok ;
} /* matchTr_genefinder */

/***************************************************************/

static void crossMissingReads (KEY tg, AC_OBJ r, AC_HANDLE h)
{
  vtxtPrintf (txt, "Transcribed_gene %s\nMissing_read %s\n\n"
	      , ac_protect (ac_key_name (tg), h)
	      , ac_protect (ac_name (r), h)
	    ) ;
} /* crossMissingReads  */

/******/

static void searchMissingRead (Array tgenes)
{
  KEY g1 ;
  GENE *hh ;
  int j1 ; char *cp ;
  AC_ITER iter = 0 ;
  AC_OBJ r = 0 ;
  AC_HANDLE h = ac_new_handle () ;

  for (j1 = 0 ; j1 < arrayMax(tgenes) - 1 ; j1++)
    {
      hh = arrayp(tgenes, j1, GENE) ;
      g1 = hh->gene ; 
      cp = messprintf ("Finf Transcribed_gene %s ; FOLLOW cDNA_clone ; FOLLOW Read NOT From_gene = %s"
		       , ac_key_name(g1)
		       , ac_key_name(g1)
		       ) ;
      iter = ac_dbquery_iter (DB, cp, h) ;
      while (ac_free (r), r = ac_iter_obj (iter))
	{
	  crossMissingReads (g1, r, h) ;
	}
      ac_free (iter) ;
    }
  ac_free (h) ;
} /* searchMissingRead */

/***************************************************************/

static void searchCtfFile (void)
{
#ifdef JUNK
  KEYSet ks = aceQuery (0, "Find EST yk*", 0) ;
  int i, n ; 
  KEY *kp, key,  _ctf= aceTag ("CTF_File",0), _scf= aceTag ("SCF_File",0) ;
  AceInstance ai = 0 ;
  char *cp ;
  FILE *f ;

  i = keySetMax (ks) ;
  printf ("Found %d est\n", i) ;
  if (!i || !_scf) return ;
  kp = arrp (ks, 0, KEY) - 1 ;
  while (kp++, i--)
    {
      key = *kp ;
      ai = aceOpenKey(*kp, 0) ;
      if (ai &&  aceGotoTag (ai, _ctf))
	{
	  if (!aceGotoChild (ai) || !aceGetText (ai, &cp))
	    printf ("NoFileName\n") ;
	  else
	    {
	      f = filopen (messprintf("CTF/%s",cp),"", "r") ;
	      if (!f)
		printf ("NoFile %s %s\n", ac_key_name(*kp), cp) ;
	      else
		{ 
		  if (!fseek (f, 0L, 2))
		    { 
		      n  = ftell (f) ;
		      if (n > 10000)
			{
			  printf ("OK  %s %s %d ", ac_key_name(*kp), cp, n) ;  
			  if (aceGotoTag (ai, _scf) && aceGotoChild (ai) && aceGetText (ai, &cp))
			     printf (" %s ", cp) ;
			  printf ("\n") ;
			}

		      else if (!n < 10)
			printf ("NoLength  %s %s %d\n", ac_key_name(*kp), cp, n) ;
		      else
			printf ("NoGood  %s %s %d\n", ac_key_name(*kp), cp, n) ;
		    }
		  else
		    printf ("NoSeekEnd  %s %s\n", ac_key_name(*kp), cp) ;
		  filclose (f) ;
		}
	    }
	}
      else
	printf ("NoTag %s %s %s\n", ac_key_name(*kp), ac_key_name (_ctf), aceErrorMessage (0)) ;
      aceFree (ai) ;
    }
#endif
  fprintf (stderr, "This function searchCTFfile is obsolete\n") ;
} /* searchCtfFile */

/***************************************************************/

static void searchScfFile (void)
{
#ifdef JUNK
  KEYSET ks = aceQuery (0, "Find EST yk* SCF_File AND DNA", 0) ;
  int i, n ; 
  KEY *kp, key,  _scf= aceTag ("SCF_File",0) ;
  AceInstance ai = 0 ;
  char *cp ;
  FILE *f ;

  i = keySetMax (ks) ;
  printf ("Found %d est\n", i) ;
  if (!i || !_scf) return ;
  kp = arrp (ks, 0, KEY) - 1 ;
  while (kp++, i--)
    {
      key = *kp ;
      ai = aceOpenKey(*kp, 0) ;
      if (ai &&  aceGotoTag (ai, _scf))
	{
	  if (!aceGotoChild (ai) || !aceGetText (ai, &cp))
	    printf ("NoText %s\n", ac_key_name(*kp)) ;
	  else
	    {
	      f = filopen (messprintf("SCF/%s",cp),"", "r") ;
	      if (!f)
		printf ("NoFile %s %s\n", ac_key_name(*kp), cp) ;
	      else
		{ 
		  if (!fseek (f, 0L, 2))
		    { 
		      n  = ftell (f) ;
		      if (n > 10000)
			printf ("OK  %s %s %d\n", ac_key_name(*kp), cp, n) ;
		      else if (!n < 10)
			printf ("NoLength  %s %s %d\n", ac_key_name(*kp), cp, n) ;
		      else
			printf ("NoGood  %s %s %d\n", ac_key_name(*kp), cp, n) ;
		    }
		  else
		    printf ("NoSeekEnd  %s %s\n", ac_key_name(*kp), cp) ;
		  filclose (f) ;
		}
	    }
	}
      else
	printf ("NoTag %s %s %s\n", ac_key_name(*kp), ac_key_name (_scf), aceErrorMessage (0)) ;
      aceFree (ai) ;
    }
#endif
  fprintf (stderr, "This function searchSCFfile is obsolete\n") ;
} /* searchScfFile */

/***************************************************************/

static void showExons (Array exons)
{
  EXON *up ;
  int i ;
  
  if (0 && arrayExists(exons))
    for (i = 0 ; i < arrayMax(exons) ; i++)
    {
      up = arrp (exons, i, EXON) ;
      printf("%d:\t%10d%10d  %d\n", i, up->x1, up->x2, up->tag) ;
    }
} /* showExons */

/***************************************************************/

static void compareGenes (int mode)
{
  GENE *hh = 0 ;
  int ii, j1 ;
  Array tGenes=0, pGenes=0, trGenes=0, gGenes = 0, rnaiGenes = 0 , products ;
  KEYSET ksTg = 0, ksPg = 0, ksTr = 0, ksRnai = 0, ksProduct ;
  KEY tg, pg, map ;
  OBJ obj ;
  float a1, a2; int x1, x2 ;

  ksTg = query(0, "FIND Transcribed_gene ") ;
  ksTr = query(0, "FIND mRNA   ; from_gene") ;  /*  IS 1F417.b OR IS G_yk323h11.b */
  ksPg = query(0, "FIND Predicted_gene ") ; /*  IS Y74C9A.3 OR IS Y74C9A.2 */
  ksRnai = query(0, "FIND RNAi  SMAP ") ;
  ksProduct =  query(0, "FIND Product IntMap  ") ;

  /* printf ("selected %s\n", ac_key_name(cosmid)) ;  */
  /*
    AC_HANDLE h = ac_new_handle () ;
    AC_KEYSET ks = ac_dbquery_keyset (snp->db, "Find  Transcribed_gene", h) ;
    printf ("found %d lines", ac_keyset_count (ks)) ;
  */

  if (mode == 0 || mode == 41 || mode == 2 || mode == 3)
    {
      pGenes = arrayCreate (arrayMax(ksPg), GENE) ; j1 = 0 ;
      for (ii = 0 ; ii < arrayMax(ksPg) ; ii++)
	{
	  pg = keySet (ksPg, ii) ;
	  obj = bsCreate (pg) ;
	  if (!obj)
	    continue ;
	  bsGoto (obj, 0) ;
	  if (bsGetKey (obj,str2tag("IntMap"), &map) &&
	      bsGetData (obj, _bsRight, _Int, &x1) &&
	      bsGetData (obj, _bsRight, _Int, &x2) )
	    {
	      if (!hh || hh->gene != pg)
		hh = arrayp(pGenes, j1++, GENE) ;
	      hh->gene = pg ;
	      hh->map = map ;
	      hh->a1 = x1 ;
	      hh->a2 = x2 ; 
	      bsDestroy (obj) ;
	    }  
	  if (!obj)
	    continue ;
	  if (bsGetKey (obj, _Map, &map) &&
	      !strncasecmp("CHROMO",ac_key_name(map), 6) &&
	      bsPushObj (obj) &&
	      bsGetData (obj, _Left, _Float, &a1) &&
	      bsGetData (obj, _Right, _Float, &a2) )
	    {
	      hh = arrayp(pGenes, j1++, GENE) ;
	      hh->gene = pg ;
	      hh->map = map ;
	      hh->a1 = a1 ;
	      hh->a2 = a2 ;
	    }
	  bsDestroy (obj) ;
	}
    }
 
  if (mode == 0 || mode == 1 || mode == 2 || mode == 4 || mode == 5)
    {
      tGenes = arrayCreate (arrayMax(ksTg), GENE) ; j1 = 0 ;
      for (ii = 0 ; ii < arrayMax(ksTg) ; ii++)
	{
	  tg = keySet (ksTg, ii) ;
	  obj = bsCreate (tg) ;
	  if (!obj)
	    continue ;
	  if (bsGetKey (obj,str2tag("IntMap"), &map) &&
	      bsGetData (obj, _bsRight, _Int, &x1) &&
	      bsGetData (obj, _bsRight, _Int, &x2) )
	    {
	      if (!hh || hh->gene != tg)
		hh = arrayp(tGenes, j1++, GENE) ;
	      hh->gene = tg ;
	      hh->map = map ;
	      hh->a1 = x1 ;
	      hh->a2 = x2 ;
	      bsDestroy (obj) ;
	    }
	  if (!obj)
	    continue ;
	  
	  bsGoto (obj, 0) ;
	  if (bsGetKey (obj, _Map, &map) &&
	      !strncasecmp("CHROMO",ac_key_name(map), 6) &&
	      bsPushObj (obj) &&
	      bsGetData (obj, _Left, _Float, &a1) &&
	      bsGetData (obj, _Right, _Float, &a2) )
    {
	      hh = arrayp(tGenes, j1++, GENE) ;
	      hh->gene = tg ;
	      hh->map = map ;
	      hh->a1 = a1 ;
	      hh->a2 = a2 ;
	    }   
	  bsDestroy (obj) ;
	}
    }

  if (mode == 0 || mode == 1 || mode == 3 || mode == 31)
    {
      trGenes = arrayCreate (arrayMax(ksTr), GENE) ; j1 = 0 ;
      for (ii = 0 ; ii < arrayMax(ksTr) ; ii++)
	{
	  tg = keySet (ksTr, ii) ;
	  obj = bsCreate (tg) ;
	  if (!obj)
	    continue ;
	  if (bsGetKey (obj,str2tag("IntMap"), &map) &&
	      bsGetData (obj, _bsRight, _Int, &x1) &&
	      bsGetData (obj, _bsRight, _Int, &x2) )
	    {
	      if (!hh || hh->gene != tg)
		hh = arrayp(trGenes, j1++, GENE) ;
	      hh->gene = tg ;
	      hh->map = map ;
	      hh->a1 = x1 ;
	      hh->a2 = x2 ;
	      bsDestroy (obj) ;
	    }
	  if (!obj)
	    continue ;
	  
	  bsGoto (obj, 0) ;
	  if (bsGetKey (obj, _Map, &map) &&
	      !strncasecmp("CHROMO",ac_key_name(map), 6) &&
	      bsPushObj (obj) &&
	      bsGetData (obj, _Left, _Float, &a1) &&
	      bsGetData (obj, _Right, _Float, &a2) )
	    {
	      hh = arrayp(trGenes, j1++, GENE) ;
	      hh->gene = tg ;
	      hh->map = map ;
	      hh->a1 = a1 ;
	      hh->a2 = a2 ;
	    }   
	  bsDestroy (obj) ;
	}
    }

#ifdef JUNK
  KEYSET ksGene = query(0, "FIND Gene  SMAP ") ;

  gGenes = aceArrayCreate (arrayMax(ksGene), GENE) ; j1 = 0 ;
  printf ("// Found %d gene\n", ksGene ? keySetMax(ksGene) : 0) ;
  for (ii = 0 ; ii < arrayMax(ksGene) ; ii++)
    {
      tg = keySet (ksGene, ii) ;
      obj = bsCreate (tg) ;
      if (!obj)
	continue ;
      if (bsGetKey (obj,str2tag("IntMap"), &map) &&
	  bsGetData (obj, _bsRight, _Int, &x1) &&
	  bsGetData (obj, _bsRight, _Int, &x2) )
	{
	  if (!hh || hh->gene != tg)
	    hh = arrayp(gGenes, j1++, GENE) ;
	  hh->gene = tg ;
	  hh->map = map ;
	  hh->a1 = x1 ;
	  hh->a2 = x2 ;
	  bsDestroy (obj) ;
	}
      if (!obj)
	continue ;

      bsGoto (obj, 0) ;
      if (bsGetKey (obj, _Map, &map) &&
	  !strncasecmp("CHROMO",ac_key_name(map), 6) &&
	  bsPushObj (obj) &&
	  bsGetData (obj, _Left, _Float, &a1) &&
	  bsGetData (obj, _Right, _Float, &a2) )
	{
	  hh = arrayp(gGenes, j1++, GENE) ;
	  hh->gene = tg ;
	  hh->map = map ;
	  hh->a1 = a1 ;
	  hh->a2 = a2 ;
	}   
      bsDestroy (obj) ;
    }
#endif
  
  if (mode == 0 || mode == 3)
    {
      products = arrayCreate (arrayMax(ksProduct), GENE) ; j1 = 0 ;
      printf ("// Found %d products\n", ksProduct ? keySetMax(ksProduct) : 0) ;
      for (ii = 0 ; ii < arrayMax(ksProduct) ; ii++)
	{
	  tg = keySet (ksProduct, ii) ;
	  obj = bsCreate (tg) ;
	  if (!obj)
	    continue ;
	  if (bsGetKey (obj,str2tag("IntMap"), &map) &&
	      bsGetData (obj, _bsRight, _Int, &x1) &&
	      bsGetData (obj, _bsRight, _Int, &x2) )
	    {
	      if (!hh || hh->gene != tg)
		hh = arrayp(products, j1++, GENE) ;
	      hh->gene = tg ;
	      hh->map = map ;
	      hh->a1 = x1 ;
	      hh->a2 = x2 ;
	      bsDestroy (obj) ;
	    }
	  if (!obj)
	    continue ;
	}
    }

  if (mode == 0 || mode == 4 || mode == 31 || mode == 41)
    {
      rnaiGenes = arrayCreate (arrayMax(ksRnai), GENE) ; j1 = 0 ;
      printf ("// Found %d rnai\n", ksRnai ? keySetMax(ksRnai) : 0) ;
      for (ii = 0 ; ii < arrayMax(ksRnai) ; ii++)
	{
	  tg = keySet (ksRnai, ii) ;
	  obj = bsCreate (tg) ;
	  if (!obj)
	    continue ;
	  if (bsGetKey (obj,str2tag("IntMap"), &map) &&
	      bsGetData (obj, _bsRight, _Int, &x1) &&
	      bsGetData (obj, _bsRight, _Int, &x2) )
	    {
	      if (!hh || hh->gene != tg)
		hh = arrayp(rnaiGenes, j1++, GENE) ;
	      hh->gene = tg ;
	      hh->map = map ;
	      hh->a1 = x1 ;
	      hh->a2 = x2 ;
	      bsDestroy (obj) ;
	    }
	  if (!obj)
	    continue ;
	  
	  bsGoto (obj, 0) ;
	  if (bsGetKey (obj, _Map, &map) &&
	      !strncasecmp("CHROMO",ac_key_name(map), 6) &&
	      bsPushObj (obj) &&
	      bsGetData (obj, _Left, _Float, &a1) &&
	      bsGetData (obj, _Right, _Float, &a2) )
	    {
	      hh = arrayp(rnaiGenes, j1++, GENE) ;
	      hh->gene = tg ;
	      hh->map = map ;
	      hh->a1 = a1 ;
	      hh->a2 = a2 ;
	    }   
	  bsDestroy (obj) ;
	}
    }
 
  /* Example of accessing the dna 
     ac_obj_dna() ;
  */


  if (tGenes) { swapA1 (tGenes) ;    arraySort (tGenes, geneOrder) ; }
  if (pGenes) { swapA1 (pGenes) ;    arraySort (pGenes, geneOrder) ; }
  if (trGenes) { swapA1 (trGenes) ;   arraySort (trGenes, geneOrder) ; }
  if (gGenes) { swapA1 (gGenes) ;    arraySort (gGenes, geneOrder) ; }
  if (rnaiGenes) { swapA1 (rnaiGenes) ; arraySort (rnaiGenes, geneOrder) ; }
  if (products) { swapA1 (products) ; arraySort (products, geneOrder) ; }

  switch (mode)
    {
    case 0:
      /*  done in makemrna.c
	searchOverlappingTG (tGenes) ;
      */
      matchGenefinder_genefinder (pGenes) ;
      matchTg_genefinder (pGenes, tGenes) ;
      matchTr_genefinder (pGenes, trGenes) ;
      matchTGene_Rnai (tGenes, rnaiGenes) ; 
      matchPGene_Rnai (pGenes, rnaiGenes) ;
      matchProduct_genefinder (pGenes, products) ;
      break ;
    case 1:
      searchOverlappingTG (tGenes) ; 
      searchOverlappingTr (trGenes) ; 
      break ;
    case 2:
      matchTg_genefinder (pGenes, tGenes) ; 
      break ;
    case 3:
      matchTr_genefinder (pGenes, trGenes) ; 
      matchProduct_genefinder (pGenes, products) ;
      break ;
    case 31:
      matchTr_genefinder (trGenes, rnaiGenes) ; 
      break ;
    case 4:
      matchTGene_Rnai (tGenes, rnaiGenes) ;
      break ;
    case 41:
      matchPGene_Rnai (pGenes, rnaiGenes) ;
      break ;
    case 5:
      /* matchGene_Gene (gGenes, gGenes) ; */
      break ;
    case 6:
      searchMissingRead (tGenes) ;
      break ;
    default:
      fprintf (stderr, "Usage gene2gene 1::Tg-tg,  2::pg-Tg 3::pg-Tr,  4::tg-RNAi, 41::pg-RNAi 5::gg-gg 6::MissingReads,  10::all\n") ;
      break ;
    }
} /* compareGenes */

/***************************************************************/


static void loopOnAllGenome (int mode)
{
  if (mode == 21) { searchScfFile () ; return ; }
  if (mode == 22) { searchCtfFile () ; return ; }
  compareGenes (mode) ;
}

/***************************************************************/

int main(int argc, char **argv)
{ 
  int mode = 0 ; /* 0: execute all modules */
  const char *error = 0 ;
  AC_HANDLE h = 0 ;

  freeinit () ;
  freeOutInit () ;
  switch (argc)
    {
    case 4:
    case 3:
      mode = atoi(argv[2]) ;
      /* fall thru */
    case 2:
      break ;
    default:
      printf ("usage: gene2gene ACEDB mode\n") ;
      printf ("mode 1::Tg-tg,  2::pg-Tg 3::pg-mRNA, 31::mRNA-RNAi, 4::tg-RNAi, 41::pg-RNAi 5::gg-gg 6::MissingReads,  0::all\n") ;
      exit (1) ;
    }

  messcrash ("This code is obsolete, use gene2chrom2 or tacembly->assembly->cdna_21 / cdna_50") ;
  showExons(0) ; /* to please the compiler */
  if (0) matchTr_Rnai (0, 0) ;  /* to plese the compiler */
  if ((DB = ac_open_db (argv[1], &error)))
    {
      const char *error ;
      fprintf(stderr, "// gene2gene mode %d, start:%s\n", mode, timeShowNow()) ;
      
      if (!ac_parse (DB, "// gene2gene2", &error, 0, h))
	fprintf(stderr, "// Sorry, I cannot write in this databaser: %s\n", error) ;
      else
	{
	  txt = vtxtHandleCreate (h) ;
	  loopOnAllGenome (mode) ;
	  fprintf (stderr, "// gene2gene created a buffer of %ld char\n", vtxtMark (txt)) ;
	  if (0) fprintf (stderr, "%s\n", vtxtPtr (txt)) ;
	  if (vtxtPtr (txt) && !ac_parse (DB, vtxtPtr (txt), &error, 0, h))
	    fprintf(stderr, "// Parse error: %s\n", error) ;
	}
      ac_db_close (DB) ;
      ac_free (h) ;
      printf("// gene2gene done:  %s\n", timeShowNow()) ;
    }
  else
    {
      fprintf (stderr, "Sorry, I cannot open database %s:  %s\n", argv[1], error) ;
      printf ("Cannot open database: %s\n", error) ;
    }

  exit (0) ;
}

/***************************************************************/
/***************************************************************/
