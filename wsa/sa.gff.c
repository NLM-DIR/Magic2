/*
 * sa.gff.c

 * This module is part of the sortalign package
 * A new RNA aligner with emphasis on parallelisation by multithreading and channels, and memory locality
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg and Greg Boratyn, NCBI/NLM/NIH
 * Created April 18, 2025

 * This is public.


 * This module implements all operations related
 * the parsing the gff annotations 

 */
/* #define ARRAY_CHECK  */

#include "sa.h"

/**************************************************************/
/**************************************************************/

static void showGBX (Array aa, int start)
{
  if (aa)
    {
      GBX *up = aa && start < arrayMax (aa) ? arrp (aa, start, GBX) : 0 ;
      if (up)
	for (int i = start ; i < arrayMax (aa) && i < start + 50 ; up++, i++)
	  fprintf (stderr, "%d :: gene:%d m:%d flag:%x strand:%x chr:%d a=%d__%d\n"
		   , i
		   , up->gene, up->mrna, up->flag, up->strand
		   , up->chrom, up->a1, up->a2
		   ) ;
    }
  return ;
} /* showGBX */

/**************************************************************/

static int exonMrnaOrder (const void *va, const void *vb)
{
  const GBX *up = va ;
  const GBX *vp = vb ;
  int n ;
  n = up->gene - vp->gene ; if (n) return n ;
  n = up->mrna - vp->mrna ; if (n) return n ;
  n = up->chrom - vp->chrom; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  return 0 ;
} /* exonMrnaOrder */

/**************************************************************/

static int exonGeneOrder (const void *va, const void *vb)
{
  const GBX *up = va ;
  const GBX *vp = vb ;
  int n ;
  n = up->gene - vp->gene ; if (n) return n ;
  n = up->chrom - vp->chrom; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->flag - vp->flag ; if (n) return -n ;
  return 0 ;
} /* exonGeneOrder */

/**************************************************************/

static int exonA1Order (const void *va, const void *vb)
{
  const GBX *up = va ;
  const GBX *vp = vb ;
  int n ;
  n = up->chrom - vp->chrom; if (n) return n ;
  n = up->a1 - vp->a1 ; if (n) return n ;
  n = up->a2 - vp->a2 ; if (n) return n ;
  n = up->strand - vp->strand ; if (n) return n ;
  n = up->gene - vp->gene ; if (n) return n ;
  n = up->mrna - vp->mrna ; if (n) return n ;
  return 0 ;
} /* exonMrnaOrder */

/**************************************************************/
/**************************************************************/
/* memory map the annotated introns and gene boxes */
void saGffBinaryParser (PP *pp)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL READONLY = TRUE ;
  DICT *geneDict = pp->geneDict = dictHandleCreate (10000, pp->h) ;
  DICT *mrnaDict = pp->mrnaDict = dictHandleCreate (10000, pp->h) ;

  if (0)
    {
      char *fNam = "/mrnaNames.sortali" ;
      char *cp = filName (pp->indexName, fNam, "rb") ;
      if (cp)
	{
	  saDictMapRead (mrnaDict, cp) ;
	  fprintf (stderr, "saGffParser found %d mrna identifiers in file %s\n", dictMax (mrnaDict), cp) ;
	}
    }
  
  if (1)
    { /* used in the introns and gene expression reports */
      char *fNam = "/geneNames.sortali" ;
      char *cp = filName (pp->indexName, fNam, "rb") ;
      if (cp)
	{
	  saDictMapRead (geneDict, cp) ;
	  fprintf (stderr, "saGffParser found %d gene identifiers in file %s\n", dictMax (geneDict), cp) ;
	}
    }

  if (1) /*annotated intron coordinates and gene name */
    {
      char *fNam = "/knownIntrons.sortali" ;
      char *cp = filName (pp->indexName, fNam, "rb") ;
      if (cp)
	{
	  pp->knownIntrons = bigArrayMapRead (cp, GBX, READONLY, pp->h) ;
	  long int nnI = bigArrayMax (pp->knownIntrons) ;
	  fprintf (stderr, "saGffParser found %ld introns coordinates in file %s\n", nnI, cp) ;
	}
    }

  if (pp->wiggle)
    { /* gene expression records */
      long int nnE = 0 ;
      BigArray geneBoxes = 0 ;
      char *fNam = "/geneBoxes.sortali" ;
      char *cp = filName (pp->indexName, fNam, "rb") ;
      if (cp)
	{
	  geneBoxes = bigArrayMapRead (cp, GBX, READONLY, h) ;
	  nnE = bigArrayMax (geneBoxes) ;
	  fprintf (stderr, "saGffParser found %ld gene fragments coordinates in file %s\n", nnE, cp) ;
	}

      if (nnE) /* split per chromosomes */
	{
	  GBX *up, *vp ;
	  long int ii ;
	  pp->geneBoxes = arrayHandleCreate (dictMax (pp->bbG.dict) + 1, Array, pp->h) ;
	  for (ii = 0, up = bigArrp (geneBoxes, 0, GBX) ; ii < nnE ; ii++, up++)
	    {
	      int chrom = up->chrom ;
	      Array chromBoxes = array (pp->geneBoxes, chrom, Array) ;
	      if (! chromBoxes)
		chromBoxes = array (pp->geneBoxes, chrom, Array) = arrayHandleCreate (10000, GBX, pp->h) ;
	      vp = arrayp (chromBoxes, arrayMax (chromBoxes), GBX) ;
	      *vp = *up ;	      
	    }
	}
    }
} /* saGffBinaryParser */

/**************************************************************/
/* transform a set of possibly overlapping gene segments into a clean bords-francs set */
static int refineGeneExons (Array bb, Array aa, KEYSET borders)
{
  int ii, jj, k, kk, n, aMax = arrayMax (aa), bMax ;
  GBX *up, *vp, *wp ;

  keySetMax (borders) = kk = 0 ;
  arrayMax (bb) = bMax = 0 ;
  
  arraySort (aa, exonGeneOrder) ;
  arrayCompress (aa) ;
  aMax = arrayMax (aa) ;
  for (ii = jj = 0, up = vp = wp = arrp (aa, 0, GBX) ; ii < aMax ; ii++, up++)
    {
      if (up->gene == vp->gene && up->chrom == vp->chrom && up->a1 == vp->a1 && up->a2 == vp->a2 && up->flag < vp->flag)
	continue ;
      if (wp < up) *wp = *up ;
      vp = wp ;
      wp++ ; jj++ ;
    }
  aMax = arrayMax (aa) = jj ;
      
  for (ii = bMax = 0, up = arrp (aa, 0, GBX) ; ii < aMax ; ii++, up++)
    {
      for (n = 1, jj = ii + 1, vp = up + 1 ; jj < aMax && vp->gene == up->gene ; jj++, vp++)
	n++ ;
      if (n == 1) /* export */
	{
	  wp = arrayp (bb, bMax++, GBX) ;
	  *wp = *up ;
	}
      else
	{
	  /* locate all borders */
	  for (k = kk = 0, vp = up ; k < n ; vp++, k++)
	    {
	      keySet (borders, kk++) = vp->a1 ;
	      keySet (borders, kk++) = vp->a2 + 1 ;
	    }
	  keySetSort (borders) ;
	  keySetCompress (borders) ;
	  kk = keySetMax (borders) ;
	  
	  /* split the vp segments on all included borders */
	  for (jj = ii, vp = up ; jj < ii + n ; jj++, vp++)
	    for (k = 0 ; k < kk ; k++)
	      {
		int a = keySet (borders, k) ;
		if (vp->a1 < a && vp->a2 >= a - 1)
		  {
		    wp = arrayp (bb, bMax++, GBX) ;
		    *wp = *vp ;
		    wp->a2 = a - 1 ;
		    vp->a1 = a ;
		  }
	      }
	}
      ii += n - 1 ; up += n - 1 ; /* shift to the end of this gene */
    }

  arraySort (bb, exonGeneOrder) ;
  /* now all segments for a given gene are redundant saucisson non overlapping slices */
  for (ii = kk = 0, up = wp = arrp (bb, 0, GBX) ; ii < bMax ; ii++, up++)
    {
      char bestFlag = up->flag ;
      for (n = 1, jj = ii + 1, vp = up + 1 ; jj < bMax && vp->gene == up->gene && vp->a1 == up->a1 && vp->a2 == up->a2 ; vp++, jj++)
	{
	  n++ ;
	  if (vp->flag > bestFlag) bestFlag = vp->flag ;
	}
      if (wp < up) *wp = *up ;
      wp->flag = bestFlag ;
      kk++ ; wp++ ;

      ii += n - 1 ; up += n - 1 ;
    }
  /* now all segements for a given gene are non-redundant best flag sausisson slices */
  /* they are never ambiguous since we work one gene at a time */
  arrayMax (bb) = bMax = kk ;

  return bMax ;  
} /* refineGeneExons */

/**************************************************************/
/* transform a set of genome wide segments into a clean set analysing ambiguities */
static int refineGeneBoxes (Array bb, Array aa, KEYSET borders)
{
  int ii, jj, k0, k, kk, n, aMax = arrayMax (aa), bMax ;
  GBX *up, *vp, *wp ;

  keySetMax (borders) = kk = 0 ;
  arrayMax (bb) = bMax = 0 ;
  
  arraySort (aa, exonA1Order) ;
  arrayCompress (aa) ;
  aMax = arrayMax (aa) ;
  for (ii = k = k0 = 0, up = arrp (aa, 0, GBX) ; ii < aMax ; ii++, up++)
    {
      for (n = 1, jj = ii + 1, vp = up + 1 ; jj < aMax && vp->chrom == up->chrom ; jj++, vp++)
	n++ ;
      if (n == 1) /* export */
	{
	  wp = arrayp (bb, bMax++, GBX) ;
	  *wp = *up ;
	}
      else
	{
	  /* locate all borders */
	  for (k = kk = 0, vp = up ; k < n ; vp++, k++)
	    {
	      keySet (borders, kk++) = vp->a1 ;
	      keySet (borders, kk++) = vp->a2 + 1 ;
	    }
	  keySetSort (borders) ;
	  keySetCompress (borders) ;
	  kk = keySetMax (borders) ;
	  /* split the vp segments on all inluded borders */
	  for (jj = ii, vp = up ; jj < ii + n ; jj++, vp++)
	    {
	      k = k0 ;
	      while (k > 0 && keySet (borders, k) > vp->a1)
		k-- ;
	      while (k < kk && keySet (borders, k) < vp->a1)
		k++ ;
	      k0 = k ;
	      
	      for (k = k0 ; k < kk ; k++)
		{
		  int a = keySet (borders, k) ;
		  
		  if (vp->a1 < a)
		    {
		      if (vp->a2 >= a - 1)
			{
			  wp = arrayp (bb, bMax++, GBX) ;
			  *wp = *vp ;
			  wp->a2 = a - 1 ;
			  vp->a1 = a ;
			}
		      else
			break ;
		    }
		}
	    }
	}
      ii += n - 1 ; up += n - 1 ; /* shift to the end of this chromosome */
    }
  
  arraySort (bb, exonA1Order) ;
  /* now all segments for a given chromosome are redundant saucisson non overlapping slices */
  for (ii = kk = 0, up = wp = arrp (bb, 0, GBX) ; ii < bMax ; ii++, up++)
    {
      GBX *bestVp = up ;
      char bestFlag = 0 ;
      int n1 = 0, n2 = 0 ;
      
      for (n = 0, jj = ii, vp = up ; jj < bMax && vp->chrom == up->chrom && vp->a1 == up->a1 && vp->a2 == up->a2 && vp->strand == up->strand ; vp++, jj++)
	{
	  n++ ;
	  switch (vp->flag)
	    {
	    case 1: n1++ ; break ;
	    case 2: /* fall through */
	    case 4: n2++ ; break ;
	    }
	  if (vp->flag > bestFlag) { bestFlag = vp->flag ; bestVp = vp ; }
	}
      
      bestVp->flag |= 0x10 ; /* master box per strand, only one contributing to sponge */
      switch (bestFlag)
	{
	case 1: /* intron */
	  if (n1 > 1) /* keep only one non attributed intron */
	    {
	      for (jj = ii, vp = up ; jj < ii + n ; jj++, vp++)
		{
		  if (vp != bestVp) vp->flag = 0 ;
		  else  vp->gene = 0 ;
		}
	    }
	  break ;
	case 2:
	case 4:
	  if (n2 > 1) /* count ambiguous genes, discard introns */
	    {
	      for (jj = ii, vp = up ; jj < ii + n ; jj++, vp++)
		{
		  if (vp->flag >= 2)  { vp->friends = n2 < 256 ? n2 : 255 ; vp->flag |= 0x8 ; } /* ambiguous */
		  else  vp->flag = 0 ;
		}
	    }
	  break ;
	}
      
      for (jj = ii, vp = up ; jj < ii + n ; jj++, vp++)      
	if (vp->flag)
	  {
	    if (wp < vp) *wp = *vp ;
	    kk++ ; wp++ ;
	  }
      ii += n - 1 ; up += n - 1 ;

    }
  /* now all segments for a given chromosome strand are ready */
  arrayMax (bb) = bMax = kk ;

  return bMax ;  
} /* refineGeneBoxes */

/**************************************************************/

long int saGffParser (PP *pp, TC *tc)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (tc->fileName, 0, h) ;
  long int nnE = 0, nnI = 0, nnGE = 0, nnME = 0, nnB = 0 ;
  int chrom = 0, gene = 0, mrna = 0 ;
  int a1, a2 ;
  int line = 0 ;
  BB *bbG = &(pp->bbG) ;
  const int twoMb = (0x1 << 21) ;
  BOOL isGff = FALSE ;
  BOOL isGtf = FALSE ;
  BOOL isIntrons = FALSE ;
  DICT *chromDict = bbG->dict ;
  DICT *geneDict = pp->geneDict ;
  DICT *mrnaDict = pp->mrnaDict ;
  Array geneExons = 0 ;
  Array geneBoxes = 0 ;
  Array geneBoxes2 = 0 ;
  Array mrnaExons = 0 ;
  KEYSET borders = 0 ;
  char strand = 0 ;
  
  if (! geneDict) geneDict = pp->geneDict = dictHandleCreate (1 << 15, pp->h) ;
  if (! mrnaDict) mrnaDict = pp->mrnaDict = dictHandleCreate (1 << 15, pp->h) ;

  if (! pp->knownIntrons)   pp->knownIntrons = bigArrayHandleCreate (100000, GBX, pp->h) ;
  if (! pp->geneBoxes) pp->geneBoxes = arrayHandleCreate (2 * dictMax (chromDict) + 1, Array, pp->h) ;

  geneExons = arrayHandleCreate (100000, GBX, h) ;
  borders = keySetHandleCreate (h) ;
  geneBoxes = arrayHandleCreate (100000, GBX, h) ;
  geneBoxes2 = arrayHandleCreate (100000, GBX, h) ;
  mrnaExons = arrayHandleCreate (100000, GBX, h) ;
  
  if (1)
    {
      const char *cp = tc->fileName ;
      int i = strlen (cp) ;
      if (i >= 8 && !strcmp (cp + i - 8, ".introns")) isIntrons = TRUE ;
      else if (i >= 11 && !strcmp (cp + i - 11, ".introns.gz")) isIntrons = TRUE ;
      else if (i >= 4 && !strcmp (cp + i - 4, ".gff")) isGff = TRUE ;
      else if (i >= 5 && !strcmp (cp + i - 5, ".gff3")) isGff = TRUE ;
      else if (i >= 4 && !strcmp (cp + i - 4, ".gtf")) isGtf = TRUE ;
      else if (i >= 7 && !strcmp (cp + i - 7, ".gff.gz")) isGff = TRUE ;
      else if (i >= 8 && !strcmp (cp + i - 8, ".gff3.gz")) isGff = TRUE ;
      else if (i >= 7 && !strcmp (cp + i - 7, ".gtf.gz")) isGtf = TRUE ;

      if (isIntrons) fprintf (stderr, "Parsing .intron file %s\n", cp) ;
      else if (isGff) fprintf (stderr, "Parsing gff file %s\n", cp) ;
      else if (isGtf) fprintf (stderr, "Parsing gtf file %s\n", cp) ;
      else messcrash ("\n\tOption I in -T tConfig_file_name expects a gff or gtf file\n\tCannot parse file %s found in tConfig file %s\n"
		      , cp, pp->tConfigFileName ) ;
    }

  if (pp->seedLength && pp->seedLength < 16)
    messcrash ("\nSorry, to study the introns defined in file %s,\n the seed length must be at least 16, not %d\n"
	       , aceInFileName (ai)
	       , pp->seedLength
	       ) ;

  if (isIntrons)
    return saIntronParser (pp, tc) ;
  

  nnI = bigArrayMax (pp->knownIntrons) ;
    
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      int type ;
      char *cp = aceInWord (ai) ;
      line++ ;
      if (! cp || *cp == '#' || *cp == '%' || *cp == '/')
	continue ;
      chrom = 0 ; a1 = a2 = 0 ;
      
      if (! dictFind (bbG->dict, hprintf (h, "G.%s", cp), &chrom))
	messcrash ("\nUnknown chromosome name %s, not matching the G targets, at line %d of file %s\n"
		   , cp
		   , line
		   , aceInFileName (ai)
		   ) ;


      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;  /* method, drop it */
      
      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;  /* type */
      type = 0 ;
      if (! strcasecmp (cp, "exon")) type = 1 ;
      if (! strcasecmp (cp, "CDS")) type = 2 ;
      if (! type) continue ;
  

      aceInStep (ai, '\t') ;
      aceInInt (ai, &a1) ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &a2) ;
      if (!a1 || !a2 || a1 == a2)
	continue ;

      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;  /* some kind of score or a dot, drop it */
      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;  /* strand */
      if (! cp)
	continue ;
      else if (*cp == '+')
	strand = 0x0 ;
      else if (*cp == '-')
	strand = 0x1 ;
      else
	continue ;
      
      aceInStep (ai, '\t') ;
      cp = aceInWord (ai) ;  /* probably . (some kind of score ?) drop it */
      if (! cp)
	continue ;
      
      aceInStep (ai, '\t') ;
      char cutter = 0 ;
      cp = aceInWordCut (ai, "\t", &cutter) ;  /* the metadata */
      if (! cp)
	continue ;
      mrna = 0 ;
      if (type >= 1) /* identify the transcript */
	{ /* expect : transcript_id "NM_001364814.1"; */
	  char *cq = strstr (cp, "transcript_id") ;
	  if (cq)
	    {
	      AC_HANDLE h1 = ac_new_handle () ;
	      ACEIN aj = aceInCreateFromText (cq + 13, 0, h1) ;
	      aceInCard (aj) ;
	      char *cr = aceInWord (aj) ;
	      if (cr)
		dictAdd (mrnaDict, cr, &mrna) ;
	      ac_free (h1) ;
	    }
	}
      gene = 0 ;
      if (type >= 1) /* identify the gene */
	{ /* expect : transcript_id "NM_001364814.1"; */
	  char *cq = strstr (cp, "gene_id") ;
	  if (cq)
	    {
	      AC_HANDLE h1 = ac_new_handle () ;
	      ACEIN aj = aceInCreateFromText (cq + 7, 0, h1) ;
	      aceInCard (aj) ;
	      char *cr = aceInWord (aj) ;
	      if (cr)
		dictAdd (geneDict, cr, &gene) ;
	      ac_free (h1) ;
	    }
	}
      
      if (type == 1 && gene) 
	{
	  GBX *up = arrayp (geneExons, nnGE++, GBX) ;
	  up->gene = gene ;
	  up->strand = strand ;
	  up->flag = 2 ;   /* exon : CDS | UTR, not yet determined */
	  up->chrom = chrom ;
	  up->a1 = a1 ;
	  up->a2 = a2 ;
	}

      if (type == 2 && gene)
	{
	  GBX *up = arrayp (geneExons, nnGE++, GBX) ;
	  up->gene = gene ;
	  up->strand = strand ;
	  up->flag = 4 ;  /* CDS */
	  up->chrom = chrom ;
	  up->a1 = a1 ;
	  up->a2 = a2 ;
	}

      if (type >= 1 && gene && mrna)
	{
	  GBX *up = arrayp (mrnaExons, nnME++, GBX) ;
	  up->gene = gene ;
	  up->mrna = mrna ;
	  up->strand = strand ;
	  up->flag = 2 ;	  
	  up->chrom = (chrom << 1) | strand ; 
	  up->a1 = a1 ;
	  up->a2 = a2 ;
	}
    }

  if (nnME) /* create the mRNA specific corresponding introns */
    {
      GBX *up, *vp, *wp ;
      int ii, jj ;
      
      arraySort (mrnaExons, exonMrnaOrder) ;
      for (up = arrp (mrnaExons, 0, GBX), ii = 0 ; ii < nnME ; ii++, up++)
	if (up->mrna)  
	  {
	    int n = 0 ;
	    for (vp = up + 1, jj = ii + 1 ; jj < nnME && vp->gene == up->gene && vp->mrna == up->mrna && vp->strand == up->strand ; vp++, jj++) 
	      {
		int chrom = up->chrom ;
		int a1 = vp[-1].a2 + 1 ;
		int a2 = vp[0].a1 - 1 ;
		int da = a2 - a1 + 1 ;
		if (da < 10)
		  continue ;
		if (0)
		  fprintf (stderr, "Intron %ld : %d %d\n", nnI, a1, a2) ;
		n++ ;
		wp = arrayp (geneExons, nnGE++, GBX) ;
		wp->gene = up->gene ;
		wp->flag = 1 ;   /* Intron */
		wp->strand = up->strand ;
		wp->chrom = chrom >> 1 ;
		wp->a1 = a1 ;
		wp->a2 = a2 ;
		
		if (da >= twoMb) /* our intron_seed format only allows 21 bits for the intron length */ 
		  continue ;

		wp = bigArrayp (pp->knownIntrons, nnI++, GBX) ;
		wp->chrom = chrom ;
		wp->mrna = 0 ;
		wp->flag = 1 ;
		wp->gene = up->gene ;
		wp->a1 = a1 ;
		wp->a2 = a2 ;
		
		vp->mrna = 0 ; /* avoid looping */
	      }
	    up->mrna = 0 ;
	    up += n - 1 ; ii += n - 1 ;
	  }
      bigArraySort (pp->knownIntrons, exonA1Order) ;
      bigArrayCompress (pp->knownIntrons) ;
      nnI = bigArrayMax (pp->knownIntrons) ;
      
      /* memory map the known introns */
      if (nnI)
	{
	  char *fNam = hprintf (h, "%s/knownIntrons.sortali", pp->indexName) ;
	  bigArrayMapWrite (pp->knownIntrons, fNam) ;
	  fprintf (stderr, "saGffParser exported %ld knownIntrons\n", bigArrayMax (pp->knownIntrons)) ;
	}
      ac_free (mrnaExons) ;
    }
  
  
  if (nnGE) /* refine the CDS/exons/introns segments, one gene at a time */
    {
      nnB = refineGeneExons (geneBoxes, geneExons, borders) ;  /* create clean gene saucisson */
      nnB = refineGeneBoxes (geneBoxes2, geneBoxes, borders) ; /* create geneome wide clean saucisson */
      
      if (nnB) /* construct and memory map the final gene boxes */
	{
	  GBX *up, *vp ;
	  int ii ;
	  BigArray geneBoxes3 = bigArrayHandleCreate (nnB, GBX, h) ;
	  arraySort (geneBoxes2, exonA1Order) ;
	  for (ii = 0, nnE = 0,  up = arrp (geneBoxes2, 0, GBX) ; ii < nnB ; ii++, up++)
	    {
	      vp = bigArrayp (geneBoxes3, nnE++, GBX) ;
	      *vp = *up ;
	    }
	  
	  char *fNam = hprintf (h, "%s/geneBoxes.sortali", pp->indexName) ;
	  bigArrayMapWrite (geneBoxes3, fNam) ;
	  fprintf (stderr, "saGffParser exported %ld gene box segments\n", bigArrayMax (geneBoxes3)) ;
	}
    }

  if (geneDict && dictMax (geneDict))
    {
      char *fNam = hprintf (h, "%s/geneNames.sortali", pp->indexName) ;
      saDictMapWrite (geneDict, fNam) ;
      fprintf (stderr, "saGffParser exported %d gene identifiers\n", dictMax (geneDict)) ;
    }
  if (0 && mrnaDict && dictMax (mrnaDict))
    {
      char *fNam = hprintf (h, "%s/mrnaNames.sortali", pp->indexName) ;
      saDictMapWrite (mrnaDict, fNam) ;
      fprintf (stderr, "saGffParser exported %d mrna identifiers\n", dictMax (mrnaDict)) ;
    }

  fprintf (stderr, "+++++++ Found %ld boxes %ld introns in file %s\n", nnE, nnI, tc->fileName) ;
  ac_free (h) ;

  return nnE + nnI ;
} /* gffParser */

/**************************************************************/

long int saIntronParser (PP *pp, TC *tc)
{
  AC_HANDLE h = ac_new_handle () ;
  long int nn = 0 ;
  BB *bbG = &(pp->bbG) ;
  ACEIN ai = aceInCreate (tc->fileName, 0, h) ;
  int chrom= 0 ;
  int a1, a2, da ;
  int line = 0 ;
  const int twoMb = (0x1 << 21) ;
  GBX *up ;
  
  if (pp->seedLength < 16)
    messcrash ("\nSorry, to study the introns defined in file %s,\n the seed length must be at least 16, not %d\n"
	       , aceInFileName (ai)
	       , pp->seedLength
	       ) ;
  if (! pp->knownIntrons)
    pp->knownIntrons = bigArrayHandleCreate (100000, GBX, pp->h) ;

  nn = bigArrayMax (pp->knownIntrons) ;
  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      char *cp = aceInWord (ai) ;
      line++ ;
      if (! cp || *cp == '#')
	continue ;
      chrom = 0 ; a1 = a2 = 0 ;

      if (! dictFind (bbG->dict, hprintf (h, "G.%s", cp), &chrom))
	messcrash ("\nUnknown chromosome name %s, not matching the G targets, at line %d of file %s\n"
		   , cp
		   , line
		   , tc->fileName
		   ) ;

      aceInStep (ai, '\t') ;
      aceInInt (ai, &a1) ;
      aceInStep (ai, '\t') ;
      aceInInt (ai, &a2) ;
      if (!a1 || !a2 || a1 == a2)
	continue ;
      da = (a1 < a2 ? a2 - a1 + 1 : a1 - a2 + 1) ;
      if (da >= twoMb) /* our format only allows 21 bits for the intron length */ 
	continue ;
      up = bigArrayp (pp->knownIntrons, nn++, GBX) ;
      if (a1 < a2)
	{
	  up->chrom = chrom << 1 ; 
	  up->a1 = a1 ;
	  up->a2 = a2 ;
	}
      else
	{
	  up->chrom = (chrom << 1) | 0x1 ;
	  up->a1 = a2 ;
	  up->a2 = a1 ;
	}
	
    }

  bigArraySort (pp->knownIntrons, exonA1Order) ;
  bigArrayCompress (pp->knownIntrons) ;
  
  /* memory map the known introns */
  if (pp->knownIntrons && bigArrayMax (pp->knownIntrons))
    {
      char *fNam = hprintf (h, "%s/knownIntrons.sortali", pp->indexName) ;
      bigArrayMapWrite (pp->knownIntrons, fNam) ;
      fprintf (stderr, "genomeCreateBinary exported %ld known_introns\n", bigArrayMax (pp->knownIntrons)) ;
    }

  showGBX (0, 0) ; /* for compiler happiness */
  fprintf (stderr, "+++++++ Found %ld introns in file %s\n", nn, tc->fileName) ;
  ac_free (h) ;

  return nn ;
} /* saIntronParser */

/**************************************************************/
/**************************************************************/
