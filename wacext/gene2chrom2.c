 
#include "ac.h"
#include "mytime.h"
#include "dict.h"
#include "vtxt.h"
 
static void usage (void) ;

static BOOL debug = FALSE, 
  showIntMap = FALSE, 
  showInt2gMap = TRUE,
  chromOnly = TRUE,
  exportJunction = FALSE ;

static DICT *Maps, *Genes ;
typedef struct { int gene, map, overlap ; int a1, a2; double gMap ; char gMapName[4]; BOOL isUp; } GENE ;
extern double intMap2Gmap (int chrom, double X1plusX2half) ;

/***************************************************************/

static BOOL locateGene (AC_OBJ Gene, int mode, Array genes)
{
  AC_OBJ ai = 0, parent, child = 0 ;
  int ii, a1, a2, x1, x2 ;
  BOOL isUp = FALSE, pleaseBreak ;
  GENE *nn ;
  AC_TABLE trs = 0 ;
  AC_HANDLE h ; 
  
  if (!Gene)
    return 0 ;
  h = ac_new_handle () ;
  a1 = a2 = x1 = x2 = 0 ;

  if (mode == 11) /* Transcript */
    {
      if (!Gene || 
	  !(parent = ac_tag_obj (Gene, "Genomic_sequence",h)))
	{ 
	  /* printf ("// no parent in \"%s\"\n", ac_name (Gene)) ; */
	  ac_free (h) ;
	  return 0 ;
	}
      
      ai = parent ;
      a1 = -1 ;
      trs = ac_tag_table (parent, "Transcripts", h) ;
      for (ii = 0 ; trs && ii < trs->rows ; ii++)
	  {
	    if (!strcmp (ac_name (Gene), ac_table_printable (trs, ii, 0, "")))
	      {
		a1 = x1 = ac_table_int (trs, ii, 1, -1) ;
		a2 = x2 = ac_table_int (trs, ii, 2, -1) ;
		break ;
	      }
	  }
      
      if (a1 == -1) 
	{
	  
	  printf ("// transcrip \"%s\" absent from parent \"%s\"\n", 
		  ac_name (Gene), ac_name (parent)) ;
	  ac_free (h) ;
	  return  0 ;
	}
      if (debug) printf("// transcript \"%s\" seq \"%s\" %d %d\n",
			ac_name (Gene), ac_name (parent), a1, a2) ;
      

      child = parent ;
    }
  else if (mode == 1 || mode == 12 || mode == 13 || mode == 14 || mode == 22)  /* Transcribed gene, mrna, genes */
    {
      if (!Gene || 
	  !(parent = ac_tag_obj (Gene, "Genomic_sequence",h)))
	{ 
	  /* printf ("// no parent in \"%s\"\n", ac_name (Gene)) ; */
	  ac_free (h) ;
	  return 0 ;
	}
      
      ai = parent ;
      a1 = -1 ;
      trs = ac_tag_table (parent
			  , mode == 1 ? "Transcribed_gene" : (mode == 12 ? "Mrnas" : ( mode == 13 ? "RNAi" : (mode == 14 ? "OST" : "Genes")))
			  , h) ;
      for (ii = 0 ; trs && ii < trs->rows ; ii++)
	  {
	    if (!strcmp (ac_name (Gene), ac_table_printable (trs, ii, 0, "")))
	      {
		a1 = x1 = ac_table_int (trs, ii, 1, -1) ;
		a2 = x2 = ac_table_int (trs, ii, 2, -1) ;
		break ;
	      }
	  }
      
      if (a1 == -1) 
	{
	  
	  printf ("// Gene \"%s\" absent from parent \"%s\" (#rows = =%d)\n", 
		  ac_name (Gene), ac_name (parent), trs ? trs->rows : 0) ;
	  ac_free (h) ;
	  return  0 ;
	}
      if (debug) printf("// Gene \"%s\" seq \"%s\" %d %d\n",
			ac_name (Gene), ac_name (parent), a1, a2) ;
      

      child = parent ;
    }
  else
    ai = child = Gene ;

  parent = ac_tag_obj (ai, "Source",h) ;
  
  pleaseBreak = FALSE ;
  while (parent)
    {
      ai = parent ;
      x1 = -1 ;
      trs = ac_tag_table (parent, "Subsequence", h) ;
      for (ii = 0 ; trs && ii < trs->rows ; ii++)
	{
	  if (!strcmp (ac_name (child), ac_table_printable (trs, ii, 0, "")))
	    {
	      x1 = ac_table_int (trs, ii, 1, -1) ;
	      x2 = ac_table_int (trs, ii, 2, -1) ;
	      break ;
	    }
	}
      if (x1 == -1)
	{ /* if the child is localised in a map with same name as parent, ok */
	  trs = ac_tag_table (child, "IntMap", h) ;
	  if (trs && trs->rows == 1 && trs->cols == 3 &&
	      !strcmp (ac_name (parent), ac_table_printable (trs, 0, 0, "")))
	    { 
	      x1 = ac_table_int (trs, 0, 1, -1) ;
	      x2 = ac_table_int (trs, 0, 2, -1) ;
	      pleaseBreak = TRUE ;
	    }
	}
      if (x1 == -1)
	{
	  printf("// error \"%s\" absent_from \"%s\"\n", ac_name (child), ac_name (parent)) ;
	  ac_free (h) ;
	  return 0 ;
	}
      if (!strcmp (ac_name (Gene), ac_name (child)))
	{
	  a1 = x1 ; a2 = x2 ;
	}
      else if (x1 < x2) { a1 += x1 - 1 ; a2 += x1 - 1 ; }
      else { a1 = x1 - a1 + 1 ; a2 = x1 - a2 + 1 ; }
      if (debug) printf("// Gene \"%s\" seq \"%s\" %d %d\n",
			ac_name (Gene), ac_name (parent), a1, a2) ;
      child = parent ;
      parent = 0 ;
      if (!pleaseBreak)
	parent = ac_tag_obj (ai, "Source",h) ;
    }


  /* register */
  if (a1 || a2)
    {
      int chrom = -1 ;
      const char *cp = 0, *gMapName = 0 ;

      if (a1 > a2) { x1 = a1 ; a1 = a2 ; a2 = x1 ; isUp = TRUE ;}
      else isUp = FALSE ;
      
      nn = arrayp (genes, arrayMax(genes), GENE) ;
      
      dictAdd (Genes, ac_name (Gene), &(nn->gene)) ;
      dictAdd (Maps, ac_name (child), &(nn->map)) ;
      if ((cp = ac_tag_printable (Gene, "Overlap_right", 0)))
	dictAdd (Genes, cp, &(nn->overlap)) ;
      nn->a1 = a1 ;
      nn->a2 = a2 ;
      nn->isUp = isUp ;
      
      cp = dictName (Maps, nn->map) ;
      if (0 && !strncmp(cp,"CHROMOSOME_",11))
	cp += 11 ;
      if (!strcmp (dictName (Maps, nn->map),"I") ||
	  (cp && !strcmp (cp,"I")))
	{ chrom = 1 ; gMapName = "I" ; }
      if (!strcmp (dictName (Maps, nn->map),"II") ||
	  (cp && !strcmp (cp,"II")))
	{ chrom = 2 ; gMapName = "II" ; }
      if (!strcmp (dictName (Maps, nn->map),"III") ||
	  (cp && !strcmp (cp,"III")))
	{ chrom = 3 ; gMapName = "III" ; }
      if (!strcmp (dictName (Maps, nn->map),"IV") ||
	  (cp && !strcmp (cp,"IV")))
	{ chrom = 4 ; gMapName = "IV" ; }
      if (!strcmp (dictName (Maps, nn->map),"V") ||
	  (cp && !strcmp (cp,"V")))
	{ chrom = 5 ; gMapName = "V" ; }
      if (!strcmp (dictName (Maps, nn->map),"X") ||
	  (cp && !strcmp (cp,"X")))
	{ chrom = 6 ; gMapName = "X" ; }
      
      if (chrom > -1)
	{
	  nn->gMap =
	    intMap2Gmap (chrom, (double)(a1+a2)/2.0) ;
	  strcpy (nn->gMapName, gMapName) ;
	}
      else
	nn->gMapName[0] = 0 ;
      }
  ac_free (h) ;
  return TRUE ;
} /* locateGene */

/***************************************************************/

static int gapOrder (const void *a, const void *b)
{
  const GENE *ga = (const GENE*)a, *gb = (const GENE *)b ;
  
  if (ga->map != gb->map)
    return ga->map - gb->map ;
  if (ga->a1 < gb->a1)
    return ga->a1 - gb->a1 ;
  return ga->a2 - gb->a2 ;
}

/***************************************************************/

static Array getGaps (AC_DB db)
{
  Array gaps = arrayCreate (100, GENE) ;
  GENE *gg ;
  AC_ITER iter =  0 ;
  int jj = 0, a1, a2 ;
  AC_OBJ gap = 0, map = 0;
  AC_TABLE gMap ;
  AC_HANDLE h = ac_new_handle () ;

  iter = ac_dbquery_iter (db, "FIND Sequence Is_Gap", h) ;
  while (ac_free (gap), gap = ac_next_obj (iter))
    {
      jj++ ;
      printf ("// Gap  \"%s\"\n" , ac_name (gap)) ;
      gMap = ac_tag_table (gap, "Intmap", h) ;
      if (gMap && 
	  (map = ac_table_obj (gMap, 0, 0, h)) &&
	  (a1 = ac_table_int (gMap, 0, 1, 0)) &&
	  (a2 = ac_table_int (gMap, 0, 2, 0))
	  )
	{
	  gg = arrayp (gaps, jj++, GENE) ;
	  dictAdd (Genes, ac_name(gap), &(gg->gene)) ;
	  dictAdd (Maps, ac_name(map), &(gg->map)) ;
	  gg->a1 = a1 ;
	  gg->a2 = a2 ;
	  gg->isUp = FALSE ;
	}
      ac_free (gMap) ;
    }
  printf ("// Found %d gaps\n", jj) ;
  
  arraySort (gaps, gapOrder) ;
  for (jj = 0 ; jj < arrayMax(gaps) ; jj++)
    {
      gg = arrayp (gaps, jj, GENE) ;
      printf ("// Map \"%s\"\tGap \"%s\"\t\t%d\t%d\n",
	     dictName (Maps, gg->map),  dictName (Genes, gg->gene), gg->a1, gg->a2) ;
    } 
  printf ("\n\n") ;
  ac_free (h) ;
  return gaps ;
}

/****************************************************************/

static BOOL getExact (GENE *nn, int zz, Array gaps) 
{
  GENE *gg ;
  int ii = arrayMax (gaps) ;
  int z1 = -1, z2 = -1 ;
  int map = nn->map ;

  gg = arrp (gaps, 0, GENE) - 1 ;
  /* z1, z2 will be a gapless interval including zz */
  while (gg++, ii--)
    {
      if (map != gg->map)
	continue ;

      /*  gotMap = TRUE ; */
      if (gg->a2 < zz) 
	z1 = gg->a2 ;

      if (gg->a1 > zz) 
	{
	  z2 = gg->a2 ;
	  break ;
	}      
    }
  /*
    printf ("// gene \"%s\" \t gotMap \t %d map \"%s\" z1=%d z2=%d zz=%d   a1 = %d a2 = %d\n",
	    ac_name (nn->gene), gotMap, ac_name (map),z1, z2, zz, nn->a1, nn->a2) ;
  */
  if (z2 == -1) z2 = nn->a1 + 1000 ;
  if (z1 < nn->a1 && z2 > nn->a2)
    return TRUE ;
  return FALSE ;
}

/****************************************************************/

static void newName (AC_DB db, Array genes)
{
  DICT *dict = dictCreate(arrayMax(genes) + 1) ;
  DICT *sdict = dictCreate(200) ;
  GENE *nn, *sn ;
  const char *ccp ;
  int a1, a2, i, n, mega, bis, section, chrom = 0 ; char code, *cp, *cq ;
  int zz, zero[7], trueZero[7] ;
  Array sections = arrayCreate (100,GENE) ;
  Array gaps = 0 ;
  BOOL isExact = FALSE ;
  BOOL s = TRUE ;  /* to find them */

  for (i=0 ; i < 7 ; i++)
    { zero[i] = 0 ; trueZero[i] = 0 ; }
  
  for (i=0 ; i < arrayMax(genes) ; i++)
    {
      nn = arrp(genes, i, GENE) ;
      
      if (!strcmp(dictName (Genes, nn->gene), s ? "dpy-5" : "1G0"))   /* G_YK1557  dpy-5 */
	trueZero[1] = nn->isUp ? -nn->a1 : nn->a2 ;
      if (!strcmp(dictName (Genes, nn->gene), s ? "dpy-10" : "2H0"))   /* G_YK3231  dpy-10 not at trueZero on rd's map !*/
	trueZero[2] = nn->isUp ? -nn->a1 : nn->a2 ;
      if (!strcmp(dictName (Genes, nn->gene), s ? "unc-32" : "3K1"))     /* G_YK91" unc-32 wrongly attributed on rd's */
	trueZero[3] = nn->isUp ? -nn->a1 : nn->a2 ;
      if (!strcmp(dictName (Genes, nn->gene), s ? "dpy-13" : "4F1"))   /* G_YK1902   dpy-13 */
	trueZero[4] = nn->isUp ? -nn->a1 : nn->a2 ;
      if (!strcmp(dictName (Genes, nn->gene), s ? "dpy-11" : "5G992"))   /* dpy-11, danielle->5H1, vrai dpy-11 = 5G992 */
	trueZero[5] = nn->isUp ? -nn->a1 : nn->a2 ;
      if (!strcmp(dictName (Genes, nn->gene), s ? "dpy-6" : "XJ70"))   /* on danielle-> XJ1,   vrai dpy-6 = XJ70 */
	trueZero[6] = nn->isUp ? -nn->a1 : nn->a2 ;
    }   
  
  for (i=0 ; i < 7 ; i++)
    {
      zero[i] = trueZero[i] ;
      if (trueZero[i] <= 0) 
	{  
	  trueZero[i] *= -1 ;
	  zero[i] = trueZero[i] ; zero[i]  %= 1000000 ;  
	  if (zero[i] > 0) zero[i] -= 1000000 ;
	}
      else 
	{ 
	  zero[i] = trueZero[i] - 1000 ; 
	  zero[i]  %= 1000000 ; 
	  if (zero[i] > 0) zero[i] -= 1000000 ; 
	}
    }
  gaps = getGaps (db) ;
  for (i=0 ; i < arrayMax(genes) ; i++)
    {
      nn = arrp(genes, i, GENE) ;

      if (nn->isUp) { a1 = nn->a2 ; a2 = nn->a1 ;}
      else { a1 = nn->a1 ; a2 = nn->a2 ; }

      if (!chromOnly || !strncasecmp("chr",dictName (Maps, nn->map),3))
	{
	  if (showIntMap)
	    printf("\nTranscribed_gene \"%s\"\nIntMap \"%s\" %d %d \n",
		    dictName (Genes, nn->gene), dictName (Maps, nn->map), a1,a2);
	  else
	    printf("\nTranscribed_gene \"%s\"\nMap \"%s\" Left %d\nMap \"%s\" Right %d \n",
		   dictName (Genes, nn->gene), dictName (Maps, nn->map), a1,  dictName (Maps, nn->map), a2) ;
	}
      if (strlen(dictName (Maps, nn->map)) > 11 && !strncmp("CHROMOSOME_",dictName (Maps, nn->map),11))
	{
	  zz = 0 ;

	  ccp = dictName (Maps, nn->map) + 11 ; cq = 0 ; chrom = 0 ; 
	  if (!strcmp(ccp,"I"))         { zz = zero[1]; cq = "1" ; chrom = 1 ; }
	  else if (!strcmp(ccp,"II"))   { zz = zero[2]; cq = "2" ; chrom = 2 ; }
	  else if (!strcmp(ccp,"III"))  { zz = zero[3]; cq = "3" ; chrom = 3 ; }
	  else if (!strcmp(ccp,"IV"))   { zz = zero[4]; cq = "4" ; chrom = 4 ; }
	  else if (!strcmp(ccp,"V"))    { zz = zero[5]; cq = "5" ; chrom = 5 ; }
	  else if (!strcmp(ccp,"X"))    { zz = zero[6]; cq = "X" ; chrom = 6 ; }

	  isExact = chrom ? getExact (nn, trueZero[chrom], gaps) : FALSE ;

	  nn->a1 -= zz ; nn->a2 -= zz ;

	  if (cq && nn->a1 && nn->a2)
	    {
	      if (!nn->isUp)   /* gene is a1 --->  a2 */
		{ 
		  n = 2 * (nn->a2/2000) + 1 ;                /* chose an odd number */
		  if (n > 0 && 1000 * n > nn->a2) n -= 2 ;   /* should be above 3p end */
		  if (1000 * n < nn->a1 &&                   /* but if outside */
		      nn->a1 - 1000 * n + 500 > 1000*(n+2) - nn->a2)   /* best to take closest to inside, but favor shift downstream */
		    n += 2 ;		    
		}
	      else         /* gene is a1 <---  a2 */
		{ 
		  n = 2 * ((nn->a1+1999)/2000)  ;            /* chose an even number below 3p */
		  if (n > 0 && 
		      1000 * n > nn->a2 &&                   /* but if outside */	
		      1000 * n - nn->a2 + 500 > nn->a1 - 1000 * (n-2))   /* best to take closest to to inside, but favor shift downstream */
		    n -= 2 ;		    
		}
              if (zero[chrom]) n += 1000 ;
	      mega = n/1000 ;

	      n -= 1000 * mega ;
	      
	      code = 'A' + mega ;  if (zero[chrom]) code-- ;
	      cp = messprintf("%s%c",cq,code) ;
	      printf("Section \"%s\"\n",cp) ;
	      dictAdd (sdict, cp, &section) ;
	      sn = arrayp (sections, section,GENE) ;
	      sn->gene = section ; sn->map = nn->map ; sn->a1 = mega ; sn->a2 = chrom ;  
	      if (zero[chrom]) sn->a1 -= 2 ;
	      cp = messprintf("%s%c%d",cq,code,n) ; bis = 1 ;
	      while (dictFind (dict,cp,0))
		{ cp = messprintf("%s%c%d_%d",cq,code,n,bis++) ; }
	      dictAdd (dict, cp, 0) ;
	      printf("DefinitiveName \"%s\"\n",cp) ;
	      if (isExact)
		 printf("DefinitiveName\n") ;
	    }
	}	
      printf("\n") ;
    }
  
  for (i = 0 ; i < arrayMax(sections) ; i++)
    {
      sn = arrp (sections, i, GENE) ;

      a1 = (sn->a1)* 1000000 ;
      if (zero[sn->a2]) a1 += 1000000 + zero[sn->a2] ; 
      a2 = a1 + 999999 ;
      if (a1 < 1) a1 = 1 ;
      printf("Section \"%s\"\nMap \"%s\" Left %d\nMap \"%s\" Right %d\n\n",
	     dictName(sdict,sn->gene), dictName (Maps, sn->map), a1,  dictName (Maps, sn->map), a2) ;
    }
  for (i=1 ; i < 7 ; i++)
    { printf("// zero[%d] = %d\n",i,zero[i]) ; }
  
  arrayDestroy (sections) ;
  arrayDestroy (gaps) ;
  dictDestroy (dict) ;
  dictDestroy (sdict) ;
}

/***************************************************************/

static BOOL loopOnAllGenome (int mode, AC_DB db)
{
  int n, ii, jj, g1, g2, g3 ;
  AC_ITER iter = 0 ;
  AC_OBJ Gene ;
  Array genes = 0 ;
  GENE *gg, *gh, *gk ;

  if (exportJunction)
    mode += 100 ;

  switch (mode)
    {
    case 1:
      iter = ac_dbquery_iter (db, "FIND Transcribed_gene",0) ;
      break ;
    case 11:
      iter = ac_dbquery_iter (db, "FIND Transcript",0) ;
      break ;
    case 12:
      iter = ac_dbquery_iter (db, "Find Mrna",0) ;
      break ;
    case 22:
      iter = ac_dbquery_iter (db, "FIND Gene",0) ;
      break ;
    case 2:
      iter = ac_dbquery_iter (db, "{FIND Predicted_gene} SETOR {find sequence  Gene_model && Source}",0) ;
      break ;
    case 13:
      iter = ac_dbquery_iter (db, "Find Rnai ; SMAP ",0) ;
      break ;
    case 14:
      iter = ac_dbquery_iter (db, "Find OST ; SMAP",0) ;
      break ;
    case 3: case 103:
      iter = ac_dbquery_iter (db, "FIND Sequence Genomic || Parts || GS_tiling",0) ;
      break ;
    case 4: case 104:
      iter = ac_dbquery_iter (db, "FIND Genome_sequence",0) ;
      break ;
    case 5: case 105:
      iter = ac_dbquery_iter (db, "FIND Section ; FOLLOW Sequence",0) ;
      break ;
    case 6: case 106:
      iter = ac_dbquery_iter (db, "find sequence is_bac",0) ;
      break ;
    default:
      messcrash ("Unrecognized mode %d",  mode) ;
      usage() ;
      break ;
    }
  
  genes = arrayCreate (30000,GENE) ;
  Gene = 0 ;
  while (ac_free (Gene), (Gene = ac_next_obj (iter)))
    {
      locateGene (Gene, mode, genes) ;
      if (debug)  break ;
    }
  printf("// Mapped %d objects\n", genes ? arrayMax(genes) : 0) ;
  if (arrayMax(genes))
    switch (mode)
      {
      case 1:
	newName(db, genes) ; 
	break ;
      case 2: case 3: case 4: case 5: case 6:  /* locate a Sequence */
	for (n = 0 ; n < arrayMax(genes) ; n++)
	  {
	    gg = arrayp (genes, n, GENE) ;
	    if (!chromOnly || !strncasecmp("chr",dictName (Maps, gg->map),3))
	      {
		int a1, a2 ;
		
		if (gg->isUp) { a1 = gg->a2 ; a2 = gg->a1 ; }
		else { a1 = gg->a1 ; a2 = gg->a2 ; }
		if (showIntMap)
		  printf("Sequence \"%s\"\nIntMap \"%s\" %d %d \n\n",
			 dictName (Genes, gg->gene), dictName (Maps, gg->map), a1, a2) ;
		else
		  printf("Sequence \"%s\"\nMap \"%s\" Left %d\nMap \"%s\" Right %d \n\n",
			 dictName (Genes, gg->gene), 
			 dictName (Maps, gg->map), a1,  
			 dictName (Maps, gg->map), a2) ;
		if (showInt2gMap && *gg->gMapName)
		  printf("Sequence \"%s\"\nFull_name (%s;%.2f)\n\n",
			 dictName (Genes, gg->gene), gg->gMapName, (float) gg->gMap) ;
	      } 
	  }
	break ;
      case 11: /* locate a Transcript */
	for (n = 0 ; n < arrayMax(genes) ; n++)
	  {
	    gg = arrayp (genes, n, GENE) ;
	    if (!chromOnly || !strncasecmp("chr", dictName (Maps, gg->map),3))
	      {
		int a1, a2 ;
		
		if (gg->isUp) { a1 = gg->a2 ; a2 = gg->a1 ; }
		else { a1 = gg->a1 ; a2 = gg->a2 ; }
		if (showIntMap)
		  printf("Transcript \"%s\"\nIntMap \"%s\" %d %d \n\n",
			 dictName (Genes, gg->gene), dictName (Maps, gg->map), a1, a2) ;
		else
		  printf("Transcript \"%s\"\nMap \"%s\" Left %d\nMap \"%s\" Right %d \n\n",
			 dictName (Genes, gg->gene), 
			 dictName (Maps, gg->map), a1,  
			 dictName (Maps, gg->map), a2) ;
	      } 
	  }
	break ;
      case 12: /* locate a Mrna */
	for (n = 0 ; n < arrayMax(genes) ; n++)
	  {
	    gg = arrayp (genes, n, GENE) ;
	    if (!chromOnly || !strncasecmp("chr",dictName (Maps, gg->map),3))
	      {
		int a1, a2 ;
		
		if (gg->isUp) { a1 = gg->a2 ; a2 = gg->a1 ; }
		else { a1 = gg->a1 ; a2 = gg->a2 ; }
		if (showIntMap)
		  printf("Mrna \"%s\"\nIntMap \"%s\" %d %d \n\n",
			 dictName (Genes, gg->gene), dictName (Maps, gg->map), a1, a2) ;
		else
		  printf("Mrna \"%s\"\nMap \"%s\" Left %d\nMap \"%s\" Right %d \n\n",
			 dictName (Genes, gg->gene), 
			 dictName (Maps, gg->map), a1,  
			 dictName (Maps, gg->map), a2) ;
	      } 
	  }
	break ;
      case 13: /* locate a RNAi  */
	for (n = 0 ; n < arrayMax(genes) ; n++)
	  {
	    gg = arrayp (genes, n, GENE) ;
	    if (!chromOnly || !strncasecmp("chr",dictName (Maps, gg->map),3))
	      {
		int a1, a2 ;
		
		if (gg->isUp) { a1 = gg->a2 ; a2 = gg->a1 ; }
		else { a1 = gg->a1 ; a2 = gg->a2 ; }
		if (showIntMap)
		  printf("RNAi \"%s\"\nIntMap \"%s\" %d %d \n\n",
			 dictName (Genes, gg->gene), dictName (Maps, gg->map), a1, a2) ;
		else
		  printf("RNAi \"%s\"\nMap \"%s\" Left %d\nMap \"%s\" Right %d \n\n",
			 dictName (Genes, gg->gene), 
			 dictName (Maps, gg->map), a1,  
			 dictName (Maps, gg->map), a2) ;
	      } 
	  }
	break ;
      case 14: /* locate OST */
	for (n = 0 ; n < arrayMax(genes) ; n++)
	  {
	    gg = arrayp (genes, n, GENE) ;
	    if (!chromOnly || !strncasecmp("chr",dictName (Maps, gg->map),3))
	      {
		int a1, a2 ;
		
		if (gg->isUp) { a1 = gg->a2 ; a2 = gg->a1 ; }
		else { a1 = gg->a1 ; a2 = gg->a2 ; }
		if (showIntMap)
		  printf("OST \"%s\"\nIntMap \"%s\" %d %d \n\n",
			 dictName (Genes, gg->gene), dictName (Maps, gg->map), a1, a2) ;
		else
		  printf("OST \"%s\"\nMap \"%s\" Left %d\nMap \"%s\" Right %d \n\n",
			 dictName (Genes, gg->gene), 
			 dictName (Maps, gg->map), a1,  
			 dictName (Maps, gg->map), a2) ;
	      } 
	  }
	break ;
      case 22: /* Gene */
	for (n = 0 ; n < arrayMax(genes) ; n++)
	  {
	    gg = arrayp (genes, n, GENE) ;
	    if (!chromOnly || !strncasecmp("chr",dictName (Maps, gg->map),3))
	      {
		int a1, a2 ;
		
		if (gg->isUp) { a1 = gg->a2 ; a2 = gg->a1 ; }
		else { a1 = gg->a1 ; a2 = gg->a2 ; }
		if (showIntMap)
		  printf("Gene \"%s\"\nIntMap \"%s\" %d %d \n\n",
			 dictName (Genes, gg->gene), dictName (Maps, gg->map), a1, a2) ;
		else
		  printf("Gene \"%s\"\nMap \"%s\" Left %d\nMap \"%s\" Right %d \n\n",
			 dictName (Genes, gg->gene), 
			 dictName (Maps, gg->map), a1,  
			 dictName (Maps, gg->map), a2) ;
		if (showInt2gMap && *gg->gMapName)
		  printf("Gene \"%s\"\nInterpolatedMap \"%s\" %.2f\n\n",
			 dictName (Genes, gg->gene), gg->gMapName, (float) gg->gMap) ;
	      } 
	  }
	break ;
      case 103:   /* junctions */
      case 104:
	for (n = 0 ; n < arrayMax(genes) ; n++)
	  {
	    gg = arrayp (genes, n, GENE) ;
	    
	    if (gg->overlap)
	      {
		
		for (ii = 0 ; ii < arrayMax(genes) ; ii++)
		  {
		    gh = arrayp (genes, ii, GENE) ;
		    if (gh->gene == gg->overlap && gh->map && gh->map == gg->map)
		      {
			int c1, c2 ;

			g1 = gg->gene ; g2 = gh->gene ;
			c1 = c2 = gg->a1 ;
			if (c1 > gg->a2) c1 = gg->a2 ;
			if (c1 > gh->a1) c1 = gh->a1 ;
			if (c1 > gh->a2) c1 = gh->a2 ;

			if (c2 < gg->a2) c2 = gg->a2 ;
			if (c2 < gh->a1) c2 = gh->a1 ;
			if (c2 < gh->a2) c2 = gh->a2 ;

			if (c2 - c1 < 20000 &&
			    (g3 = gh->overlap))
			  {
			    for (jj = 0 ; jj < arrayMax(genes) ; jj++)
			      {
				gk = arrayp (genes, jj, GENE) ;
				if (gk->gene == g3 && gk->map && gk->map == gg->map)
				  {

				    if (c1 > gk->a1) c1 = gk->a1 ;
				    if (c1 > gk->a2) c1 = gk->a2 ;
				    
				    if (c2 < gk->a1) c2 = gk->a1 ;
				    if (c2 < gk->a2) c2 = gk->a2 ;

				    printf("Sequence \"%s_%s_%s\"\nJunction\nIntMap \"%s\" %d %d\nParts \"%s\"\nParts \"%s\"\nParts \"%s\"\n\n",
					   dictName (Genes, g1), 
					   dictName (Genes, g2), 
					   dictName (Genes, g3),
					   dictName (Maps, gg->map), c1, c2,
					   dictName (Genes, g1), 
					   dictName (Genes, g2), 
					   dictName (Genes, g3)
					   ) ;
				    printf ("Sequence \"%s\"\nSubsequence \"%s_%s_%s\" %d %d\n\n",
					    dictName (Maps, gg->map),
					    dictName (Genes, g1), 
					    dictName (Genes, g2), 
					    dictName (Genes, g3),
					    c1, c2) ;
				    break ;
				  }
			      }
			  }
			else
			  {
			    printf("Sequence \"%s_%s\"\nJunction\nIntMap \"%s\" %d %d\nParts \"%s\"\nParts \"%s\"\n\n",
				   dictName (Genes, g1), 
				   dictName (Genes, g2), 
				   dictName (Maps, gg->map), c1, c2,
				   dictName (Genes, g1), 
				   dictName (Genes, g2)
				   ) ;
			    printf ("Sequence \"%s\"\nSubsequence \"%s_%s\" %d %d\n\n",
				    dictName (Maps, gg->map),
				    dictName (Genes, g1), 
				    dictName (Genes, g2), 
				    c1, c2) ;
			  }
			break ;
		      }
		  }
	      }	       
	  }	
      }
  arrayDestroy (genes) ;
  
  return TRUE ;
}


/***************************************************************/
/***************************************************************/
 
static void usage (void)
{ fprintf (stderr, "Usage: gene2chrom [options] $ACEDB\n"
	   "  version 2, 5 july 1999\n"
           "  exports coordinates in top most sequence parent\n"
	   "\t-tg:  export transcribed_gene new names and coordinates\n"
	   "\t-tr:  export transcripts coordinates\n"
	   "\t-mr:  export mrna coordinates\n"
	   "\t-gg:  export genes coordinates\n"
	   "\t-pg:  export predicted_gene coordinates\n"
	   "\t-rnai:  export rnai coordinates\n"
	   "\t-ost:  export ost coordinates\n"
	   "\t-gs:  export \'sequence genomic\' coordinates\n"
           "\t-GS:  export \'Genome_sequence\' coordinates\n"
           "\t-bac:  export \'sequence is_bac\' coordinates\n"
           "\t-xs:  export \'section->sequence\' coordinates\n"
	   "\t-i: export IntMap values (to avoid float rounding)\n"
	   "\t-G: export interpolated gmap\n"
	   "\t-chrom: (default) export only if topmost sequence is called chr*\n"
	   "\t-any: export whatever the name of the topmost sequence\n"
	   "\t-j: export gs or GS junctions based on overlap info\n"
	   "\t-test: does some test program\n"
	   );
}

/***************************************************************/

int main(int argc, char **argv)
{ 
  AC_DB db = 0 ;
  int mode = 0 ;
  Maps = dictCreate (256) ;
  Genes = dictCreate (25000) ;

  for (--argc, ++argv ; argc > 1 ; --argc, ++argv)
    if (!strcmp (*argv, "-h"))
      usage () ;
    else if (!strcmp (*argv, "-help"))
      usage () ;
    else if (!strcmp (*argv, "---help"))
      usage () ;
    else if (!strcmp (*argv, "-tg"))
      mode = 1 ;
    else if (!strcmp (*argv, "-tr"))
      mode = 11 ;
    else if (!strcmp (*argv, "-mr"))
      mode = 12 ;
    else if (!strcmp (*argv, "-gg"))
      mode = 22 ;
    else if (!strcmp (*argv, "-pg"))
      mode = 2 ;
    else if (!strcmp (*argv, "-rnai"))
      mode = 13 ;
    else if (!strcmp (*argv, "-ost"))
      mode = 14 ;
    else if (!strcmp (*argv, "-gs"))
      mode = 3 ;
    else if (!strcmp (*argv, "-GS"))
      mode = 4 ;
    else if (!strcmp (*argv, "-xs"))
      mode = 5 ;
    else if (!strcmp (*argv, "-bac"))
      mode = 6 ;
    else if (!strcmp (*argv, "-chrom"))
      chromOnly = TRUE ;
    else if (!strcmp (*argv, "-any"))
      chromOnly = FALSE ;
    else if (!strcmp (*argv, "-i"))
      showIntMap = TRUE ;
    else if (!strcmp (*argv, "-j"))
      exportJunction = TRUE ;
    else if (!strcmp (*argv, "-G"))
      showInt2gMap = TRUE ;
    else if (!strcmp (*argv, "-test"))
      mode = -1 ;
    else
      { fprintf (stderr, "Unrecognized option %s\n", *argv) ;
        usage() ;
	exit (-1) ;
      }
  
  if (!argc)
      { fprintf (stderr, "Database location MUST be specified\n") ;
        usage() ;
	exit (-1) ;
      }

  {
    const char *error = 0 ;
    if (!(db = ac_open_db (*argv, &error)))
      messcrash ("Failed to open db %s, error %s", argv[1], error) ;
  }

  printf("// gene2chrom mode  start:%s\n", timeShowNow()) ;

  loopOnAllGenome(mode, db) ;

  ac_db_close (db) ;   /* savesession */
  printf("// gene2chrom mode = %d done:  %s\n", mode, timeShowNow()) ;
  return 0 ;
}

/***************************************************************/
/***************************************************************/
