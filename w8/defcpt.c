/*  File: defcpt.c
 *  Author: Ulrich Sauvage (ulrich@kaa.crbm.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Dec 11 17:31 1998 (fw)
 *	-	static defCptOut initialized in defCptOpen() rather than file scope
 *	-	return value for defCptDoJoinDiagonal
 * * Jul 24 14:19 1996 (ulrich) 
 * * Apr 16 19:01 1996 (ulrich)
 * Created: Tue Aug 30 13:09:34 1994 (ulrich)
 *-------------------------------------------------------------------
 */

/* $Id: defcpt.c,v 1.64 2020/05/30 16:50:34 mieg Exp $ */

#include "acedb.h"
#include "freeout.h"

#ifndef NON_GRAPHIC
#include "display.h"
#include "chrono.h"
#include "plot.h"

#define DEFCPTGET(name) \
  DEFCPT look ; \
  if (!graphAssFind (&DEFCPTMAG, &look)) \
    messcrash ("defcpt-graph not found in %s", name) ; \
  if (look->magic != DEFCPTMAG) \
    messcrash ("%s received a wrong defcpt-pointer", name)

#endif /* !NON_GRAPHIC */

static int DEFCPTMAG  =  727830152 ;  /* repeated in blyDnaGet */

#include "command.h"		/* for tStatus */
#include "lex.h"
#include "bs.h"
#include "a.h"
#include "dna.h"
#include "cdna.h"
#include "sysclass.h"
#include "classes.h"
#include "systags.h"
#include "tags.h" 
#include "pick.h"
#include "topology.h"
#include "query.h"
#include "bump.h"
#include "interval.h"
#include "dnaalign.h"
#include "fingerp.h"
#include "mytime.h"
#include "session.h"
#include "parse.h"
#include "ctf.h"
#include "../wabi/chromorepeats.h"
#include "acembly.h"  
 
static BOOL DEBUG = FALSE ;    /* additional printouts */
extern void monDnaForget(DEFCPT look, KEY key) ;

typedef struct { int g1, g2 ; } GENE_PAIR ;
typedef struct { KEY key ; int tai ; } KEY_INT ;
typedef struct { KEY gene, ch ; float p ; } GMAPST ;
typedef struct { int marker, first ; float pos ; } SORT_DAT ;

static DEFCPT nonGraphicLook = 0 ;
static DEFCPT myLook = 0 ;

static int nmrid = 0 ;
static Associator assDefLook = 0 ;
static Stack diffaction = 0 ;

static FILE *defCptOut = 0 ;

static void defCptOrderBySubclones (DEFCPT look) ;
static void lookArrayDestroy (DEFCPT look) ;
static void uLocalDoDestroy (DEFCPT look) ;
static void defCptRestart (DEFCPT look) ;

#ifndef NON_GRAPHIC
static void localDestroy (void) ;
static void defcomputeDraw (DEFCPT look, char *text) ;
static Graph myGraph = 0 ;
#endif /* !NON_GRAPHIC */

static int nblook = 0 ;

#define localDoDestroy(x) ((x) ? uLocalDoDestroy(x), x=0, TRUE : FALSE)

static FREEOPT mkaction[] =
{ { 82, "Ace.mbly"},
    {'q', "Quit [<filename>] : Quits ace.mbly [and reports on file]"},
    {'l', "Load -all | -active | name : Loads sequences, discarding previous set (ex: load a1*)"},
    {'A', "add query : Adds in result of the query (ex: Add >?Clone; >Cloning_Vector)"},
    {'n', "newSCF : Adds new_scf reads to current assembly"},
    {'b', "Basecall : Replace external by acembly base call for tace-active set"},   
    {'g', "get n : get n oligos per sequence (ex: get 12, sort 90, assemble 8)"},
    {'s', "sort d : sorts results of get to limit 50 <d<=100"},
    {'a', "assemble e : assemble results of sort at error rate e %%"},
    {'j', "join [diagonal] : Joins at low strigency preassembled contigs"},
    {'J', "Juxtapose : Juxtapose the contigs into an assembly"},
    {'F', "fix : Recomputes the consensus from the elementary traces"},
    {'i', "insert : Tries to assemble into the same contigs all other reads"},
    { 10, "extend : Extend the reads beyond the present contig"},
    {'L', "Align <reference> : Align on wild type (i.e. on a reference sequence)"},
    {'V', "Vector_Clipping [-f] [-o out.ace] : look for vector [-f forces reevaluation] in loaded or active set, if no out.ace, modifies the database"},
    {'Q', "BadQuality [-n 6 -d 100] : Labels sequences with more than 6 n in best 100 bp < 600 bp  in loaded or active set"},
    {'C', "Clip_on [Excellent, Good, Fair, Tile, Consolidate] : Evaluates quality, {then  move clips]"},
    { 21, "Order_by_Size"},
    { 22, "Order_by_Subclones"},
    { 23, "Make_Subclones : Identify reads differing only by endding as from same subclone"},
    { 24, "Make_subseq [-dna] prefix subprefix size overlap [-worm]: split active set into subsequeneces of source:name"},
    { 25, "Make_dna_cover -n 20 [-justIntron] [-worm] [-pg] [-o filename]: \n"
      "export a tab delimited file of all tags of length n present in each active gene/tg/mrna,\n"
      "in intron case they are left padded with N so the exon junction is right of base n,\n"
      "the fragments are called gene.[XIAN] X=exon, I=intron A=polyA n=SLn with n=1,2,3,4..."},
    { 13, "autoedit"} ,
/*   18, "TrainNN",  */
    {'R', "Rename name : Saves as name, destroy any previous sequence called name or name.#"},
    {'S', "save_as [<name>] : Saves as name, or by incrementing (ex: a5 --> saved as a6)"},
    {1001, "cDNA_1 [chrom].. align all est on all [or a list of chromosomes]"},
    {1002, "cDNA_2 [chrom].. align active set of est est on all, do not destroy previous genes [or a list of chromosomes]"},
    {1003, "cDNA_3 align est on one sequence"},
    {1021, "cDNA_21 [chrom].. create link_ and gap_ on all [or a list of chromosomes]"},
    {3001, "cDNA_e1 [chrom].. align all est on all [or a list of chromosomes]"},
    {3002, "cDNA_e2 align active est on all [or a list of chromosomes]"},
    {1301, "cDNA_p1 -f filename  [-s step]  [-c comb] create word frequency table of active sequence set"},
    {1302, "cDNA_p2 -f frequencyfile  -l n [-t | -o outfile] preassign all reads [just tag] with at least n hits  to active sequence set"},
    {1005, "cDNA_5 [-f file] [-alter] [-nonSliding] [-coord] [-short 1|2] export confirmed introns"},
    {1006, "cDNA_6 [-f file] export exon to search oligos"},
    {1061, "cDNA_virtual_tile create virtual tiles for super_long predicted genes"},
    {1071, "cDNA_71 [-locally] realign all transcribed_genes"},
    {1073, "cDNA_73 [-locally]  [-fuse_locally] [-clean_killed_mRNA] [-split_cloud] [-repeats] [-rubber]  realign active set of transcribed_genes"},
  {1074, "Percolate active set of tg: create pseudo-EST percolating through the most favorable set of exons"}, 
  {1075, "Tag_shedded_genes Tag the principal and shedded active transcribed_genes"},
  {1076, "cDNA_76 Split in active set of transcribed genes  the XG_ exons on Xends_ signals"},
  {1077, "cDNA_77 Kill premRNA echoes on other starnd\n"},
  {1031, "cDNA_31 allocate all cDNAs to best gene"},
  {1032, "cDNA_32 allocate active cDNA to best gene"},
  {1038, "cDNA_bury bury cdna_clones in mrna>200 clones in active list of Tg"},
  {1039, "cDNA_39 report alignment quality of active est keyset"},
  {1047, "Translate_est translate est and create product <est>.est"},
  {1048, "Transfer_pg [name] Transfer the predicted genes to the mRNA class"},
  {1049, "Cluster_pg Cluster alternative predicted genes"},
  {1050, "cDNA_50 Operons overlaps and comparison to predicted genes"},
  {1051, "SlideVector  active ESTs till cdna differs from genomic"},
  {1052, "SlidePolyA  active ESTs till cdna differs from genomic, and count terminal A"},
  {1053, "PolyAtag  active mRNAs : tag Internal priming and label trustable polyA sites"},
  {1062, "Align_Sage  [-on_mrna]  align active set of  sage on all genomic sequences"},
  {1063, "Align_Primers [-on_mrna] align active set of primers on all genomic [mRNA] sequences"},
  {1064, "Align_RNAi [-f outfile.ace] double align the active Rnai set based on its allready aligned Primers"},
  {1311, "Genebox_all_primers  [-f file] export primer-pairs for all geneboxes"},
  {1313, "Genebox_primers  [-f file] export primer-pairs for active set of geneboxes"},
  {1060, "cDNA_60 export active set of sequence in CTF format"},
  {1160, "zone_export -o outfile  export active set gene as fasta zones for probaalign"},
  {1070, "chromoRepeats analyse all repeats in active sequence set"},
  {1080, "cDNA_80 create cDNA_clone info, for active set"},
  {1081, "Duplicate_clones containing several reads, for active set"},
  {1090, "cDNA_Fiche Run fiche on all genes present active set"},
  {1091, "cDNA_Kantor Read the kantor info in active set of mrna"},
  {1092, "cDNA_Compare [-coding] Compare the [coding] exons/introns in active set of tg"},
  {1095, "cDNA_Tiling add tiling path to the active set of mrna"},
  {1096, "Gene_intron add intron class info in active set of genes"},
  {1097, "Gene_probe_walls add probe walls (to split probe-sets) in active set of genes"},
  {1098, "Intron_hierarchy look for introns including another one"},
  {1101, "cDNA_FlipAllGenes Flips all transcribed_genes if need be"},
  {1102, "cDNA_FlipGenes Flips the active transcribed_genes if need be"},
  {1111, "cDNA_Flag_suspected_internal_deletion [ -ignore] Flags unique clones defining a weird alternative intron"},
  {1131, "cDNA_Rename [-f file | -g [-c chrom]] rename active set of genes using funny phonems, -g:go do rename, -c -f select newname"},
  {1121, "RNA_editing flag active set of transcribed_genes looking for rna editing"},

  {1200, "saucisse find most frequent words"},
  {1210, "Tiling_regions (gifacembly mode only) precompute active set of tiling regions"},
  {2001, "Align_1 align all reads on all reference sequences"},
  {2002, "Align_2 align tace active reads on reference sequence"},
  {2003, "Align_3 align tace active sequences on all genomic sequences"},

    {'T', "tace [-f <file>y] : execute tace commands until ^D"},
    {'r', "Status: reports on memory usage"},
    {16, "Report : Reports present length of contigs"},
     {'B', "Statistics : number the differences between sequences and choosen consensus"}
/*
      case 'v':
	key = dnaAlignDoMakeSuperLink (ksNew, 0) ;
	freeOut (messprintf ("\n// made Assembly %s\n", name (key))) ;
	break ;
*/
} ;

/***************************************************************/
/**************************    Display    **********************/
/***************************************************************/

#define TESCOLO \
   if (!arrayExists(look->colOrder)) \
      { int i = look->nd, *ip ; \
	look->colOrder = arrayCreate(i, int) ; \
	ip = arrayp(look->colOrder, i - 1, int) ; \
	while(i--) \
	  *ip-- = i ; \
      }

#define TESLIN \
   if (!arrayExists(look->linOrder)) \
      { int i = look->nm, *ip ; \
	look->linOrder = arrayCreate(i, int) ; \
	ip = arrayp(look->linOrder, i - 1, int) ; \
	while(i--) \
	  *ip-- = i ; \
      }

/***************************************************************/

static int mardefOrder(const void *a, const void *b)
{
  return
    ((*(const KEY*)a) & F_WHO) > ((*(const KEY*)b) & F_WHO)  ?
      1 : ((*(const KEY*)a) & F_WHO) == ((*(const KEY*)b) & F_WHO) ? 0 : - 1 ;
}

/***************************************************************/

typedef struct { int marker, nn ; } intPair ;

#ifndef NON_GRAPHIC
static int intPairOrder(const void *a, const void *b)
{
  return
    ((const intPair *) a)->nn - ((const intPair *) b)->nn ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static Array defCptTranspMat(Array donnee)
{ Array result = 0 ;
  int *ip, i = arrayMax(donnee) ;

  result = arrayCreate(i, int) ;
  array(result, i - 1, int) = 0 ;
  ip = arrp(donnee, i - 1, int) + 1 ;
  while(ip--, i--)
    arr(result, *ip, int) = i ;
  return result ;
} /* defCptTranspMat */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defTreeSortM(void)
{ KEYSET ksm, ksd ;
  KEY *keyp, *keyq ;
  int i, j, n, nm, k = 0, imax, jmax, kmax, *ip, who ;
  Array tm, newMarkers = 0, defOrder = 0 ;
  DEFCPTGET("defTreeSortM") ;

  tm = arrayCreate(look->nm, int) ;
  array(tm, look->nm - 1, int) = 0 ;
  defOrder = defCptTranspMat(look->colOrder) ;
  imax = look->nd ;
  ip = arrp(look->colOrder, 0, int) ;
  for (i = 0 ; i < imax ; ip++, i++)
    { ksm = arr(look->marInDef, *ip, KEYSET) ;
      jmax = keySetMax(ksm) ;
      keyp = arrp(ksm, 0, KEY) - 1 ;
      newMarkers = arrayReCreate(newMarkers, 32, intPair) ;
      nm = 0 ;
      while(keyp++, jmax--)
	{ if ((*keyp & F_NO) || arr(tm, *keyp & F_WHO, int) || (*keyp & F_V_FLAG))
	    continue ;
	  who = *keyp & F_WHO ;
	  ksd = arr(look->defInMar, who, KEYSET) ;
	  kmax = keySetMax(ksd) ;
	  keyq = arrp(ksd, 0, KEY) - 1 ;
	  n = 0 ;
	  while(keyq++, kmax--)
	    if (!(*keyq & F_NO) && (j = arr(defOrder, *keyq & F_WHO, int)) > i &&
		!(*keyq & F_V_FLAG))
	      n += j - i ;  /* on compte le poids des restant sur la ligne et non le nombre d'elements restant sur la ligne*/
	  arrayp(newMarkers, nm, intPair)->marker = who ;
	  arrp(newMarkers, nm++, intPair)->nn = n ;
	  arr(tm, who, int) = 1 ;
	}
      if (nm)
	arraySort(newMarkers, intPairOrder) ;
      for (j = 0 ; j < nm ; j++)
	arr(look->linOrder, k++, int) = arrp(newMarkers, j, intPair)->marker ;
    }
            /* add the last ones in */
  ip = arrp(tm, 0, int) ;
  imax = look->nm ;
  for (i = 0 ; i < imax ; ip++, i++)
    if (!(*ip))
      arr(look->linOrder, k++, int) = i ;
  arrayDestroy(tm) ;
  arrayDestroy(newMarkers) ;
  arrayDestroy(defOrder) ;
  look->mapStatus |= F_SOMAR ;
  look->mapStatus &= ~F_SHOWM ;
  defcomputeDraw (look, 0) ;
} /* defTreeSortM */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static int defCptSortOrder(const void *a, const void *b)
{ int i ;

  i = (((const SORT_DAT *)a)->pos > ((const SORT_DAT *)b)->pos) ? 1 : 0 ;
  if (!i)
    i = (((const SORT_DAT *)a)->pos < ((const SORT_DAT *)b)->pos) ? -1 : 0 ;
  if (!i)
    i = (((const SORT_DAT *)b)->first - ((const SORT_DAT *)a)->first) ; /* ordre decroissant car on veut le premier d'abord */
  return i ;                  /* was a)->pos in last test ??? */
} /* defCptSortOrder */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptSortMarkers(void)
{ KEYSET *ksmp ;
  KEY *keyp ;
  int i, j, n, nm, k, imin, imax, jmax, *ip, who ;
  Array tm, newMarkers = 0, defOrder = 0 ;
  DEFCPTGET("defCptSortMarkers") ;

  TESCOLO ;
  TESLIN ;
  tm = arrayCreate(look->nm, int) ;
  array(tm, look->nm - 1, int) = 0 ;
  defOrder = defCptTranspMat(look->colOrder) ;
  newMarkers = arrayCreate(look->nm, SORT_DAT) ;
  ksmp = arrp(look->defInMar, 0, KEYSET) ;
  imax = look->nm ;
  nm = 0 ;
  for (i = 0 ; i < imax ; ksmp++, i++)
    { if (keySet(*ksmp, 0) & F_V_FLAG)
	continue ;
      jmax = keySetMax(*ksmp) ;
      keyp = arrp(*ksmp, 0, KEY) - 1 ;
      n = 0 ;
      k = 0 ;
      imin = look->nd ;
      while(keyp++, jmax--)
	{ if ((*keyp & F_NO) || (*keyp & F_V_FLAG))
	    continue ;
	  who = *keyp & F_WHO ;
	  j = arr(defOrder, who, int) ;
	  n += j ;
	  k++ ;
	  if (j < imin)
	    imin = j ;
	}
      if (k)
	{ arrayp(newMarkers, nm, SORT_DAT)->marker = i ;
	  arrp(newMarkers, nm, SORT_DAT)->first = imin ;
	  arrp(newMarkers, nm++, SORT_DAT)->pos = n / k ;
	  arr(tm, i, int) = 1 ;
	}
    }
  if (nm)
    arraySort(newMarkers, defCptSortOrder) ;
  k = 0 ;
  for (j = 0 ; j < nm ; j++)
    arr(look->linOrder, k++, int) = arrp(newMarkers, j, SORT_DAT)->marker ;
            /* add the last ones in */
  ip = arrp(tm, 0, int) ;
  imax = look->nm ;
  for (i = 0 ; i < imax ; ip++, i++)
    if (!(*ip))
      arr(look->linOrder, k++, int) = i ;
  arrayDestroy(tm) ;
  arrayDestroy(newMarkers) ;
  arrayDestroy(defOrder) ;
  look->mapStatus |= F_SOMAR ;
  look->mapStatus &= ~F_SHOWM ;
  defcomputeDraw (look, 0) ;
} /* defCptSortMarkers */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static int defCptWhoSup(Array donnee, Array result)
{ KEYSET *ksp ;
  int i, cc = 0 ;
  char *cp ;

  i = arrayMax(donnee) ;
  array(result, i - 1, char) = 0 ;
  cp = arrp(result, i - 1, char) + 1 ;
  ksp = arrp(donnee, i - 1, KEYSET) + 1 ;
  while(ksp--, cp--, i--)
    { 
      if (keySetMax(*ksp))
	{
	  if (keySet(*ksp, 0) & (((unsigned int)248) << 24)) /* flag autre que F_ZERO */
	    { cc++ ;
	      *cp = 1 ;
	    }
	}
      else   /* pas de marqueurs pour la def */
	{ cc++ ;
	  *cp = 1 ;
	}
    }
  return cc ;
} /* defCptWhoSup */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static int dupCostYes(DEFCPT look)
{ int cost = - look->nd, i, imax, j, *ip, who, costmax = look->nd ;
  Array inter = 0, defSup = 0, marSup = 0 ;
  KEYSET *ksp ;
  KEY *keyp ;
  char *cp ;

  defSup = arrayCreate(look->nd, char) ;
  cost += defCptWhoSup(look->defInMar, defSup) ;
  marSup = arrayCreate(look->nm, char) ;
  defCptWhoSup(look->marInDef, marSup) ;
  inter = arrayCreate(look->nd, char) ;
  array(inter, look->nd, char) = 0 ;
  imax = look->nm ;
  ip = arrp(look->linOrder, 0, int) - 1 ;
  while(ip++, imax--)
    { 
      if (arr(defSup, *ip, char))
	{ costmax-- ;
	  continue ;
	}
      ksp = arrp(look->defInMar, *ip, KEYSET) ;
      j = keySetMax(*ksp) ;
      keyp = arrp(*ksp, 0, KEY) ;
      who = *keyp & F_WHO ;
      for (i = 0 ; i < look->nd ; i++)
	{ 
	  if (i == who)
	    { 
	      if (arr(marSup, who, char))
		continue ;
	      if (*keyp & F_NO)
		{ if (arr(inter, who, char))
		    { cost += 2 ; /* pour que les inconnues comptent comme 1/2 dans le cost */
		      arr(inter, who, char) = 0 ;
		    }
		}
	      else /* F_YES */
		arr(inter, who, char) = 1 ;
	      if (--j)
		{ keyp++ ;
		  who = *keyp & F_WHO ;
		}
	    }
	  else if (arr(inter, i, char) && !arr(marSup, i, char))
	    { cost++ ; /* inco donc 1/2 */
	      arr(inter, i, char) = 0 ;
	    }
	}
    }
  i = look->nd ;
  cp = arrp(inter, 0, char) ;
  while(i--)
    if (*cp++)
      cost +=2 ;
  arrayDestroy(defSup) ; 
  arrayDestroy(marSup) ; 
  arrayDestroy(inter) ;
  return (50 * cost) / costmax ; /* le cost est le % de trou dans chaque def */
} /* dupCostYes */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static int dupCostNo(DEFCPT look)
{ int cost = - look->nd, i, imax, j, *ip, who, costmax = look->nd ;
  Array inter = 0, defSup = 0, marSup = 0 ;
  KEYSET *ksp ;
  KEY *keyp ;
  char *cp ;

  defSup = arrayCreate(look->nd, char) ;
  cost += defCptWhoSup(look->defInMar, defSup) ;
  marSup = arrayCreate(look->nm, char) ;
  defCptWhoSup(look->marInDef, marSup) ;
  inter = arrayCreate(look->nd, char) ;
  array(inter, look->nd - 1, char) = 0 ;
  imax = look->nm ;
  ip = arrp(look->linOrder, 0, int) - 1 ;
  while(ip++, imax--)
    { 
      if (arr(defSup, *ip, char))
	{ costmax-- ;
	  continue ;
	}
      ksp = arrp(look->defInMar, *ip, KEYSET) ;
      j = keySetMax(*ksp) ;
      keyp = arrp(*ksp, 0, KEY) - 1 ;
      while(keyp++, j--)
	{ who = *keyp & F_WHO ;
	  if (arr(marSup, who, char))
	    continue ;
	  if (*keyp & F_NO)
	    { if (arr(inter, who, char))
		{ cost++ ;
		  arr(inter, who, char) = 0 ;
		}
	    }
	  else /* F_YES */
	    arr(inter, who, char) = 1 ;
	}
    }
  i = look->nd ;
  cp = arrp(inter, 0, char) ;
  while(i--)
    if (*cp++)
      cost++ ;
  arrayDestroy(defSup) ; 
  arrayDestroy(marSup) ; 
  arrayDestroy(inter) ;
  return (100 * cost) / costmax ;
} /* dupCostNo */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static int dupCostED(DEFCPT look)
{ int cost = - look->nm, imax, j, *ip, who, costmax = look->nm ;
  Array inter = 0, defSup = 0, marSup = 0 ;
  KEYSET *ksp ;
  KEY *keyp ;

  defSup = arrayCreate(look->nd, char) ;
  defCptWhoSup(look->defInMar, defSup) ;
  marSup = arrayCreate(look->nm, char) ;
  cost += defCptWhoSup(look->marInDef, marSup) ;
  inter = arrayCreate(look->nm, char) ;
  array(inter, look->nm - 1, char) = 0 ;
  imax = look->nd ;
  ip = arrp(look->colOrder, 0, int) - 1 ;
  while(ip++, imax--)
    { 
      if (arr(marSup, *ip, char))
	{ costmax-- ;
	  continue ;
	}
      ksp = arrp(look->marInDef, *ip, KEYSET) ;
      j = keySetMax(*ksp) ;
      keyp = arrp(*ksp, 0, KEY) - 1 ;
      while(keyp++, j--)
	{ who = *keyp & F_WHO ;
	  if (arr(defSup, who, char))
	    continue ;
	  if ((*keyp & F_PERE) && !arr(inter, who, char))
	    { cost++ ;
	      arr(inter, who, char) = 1 ;
	    }
	  else if ((*keyp & F_MERE) && arr(inter, who, char))
	    { cost++ ;
	      arr(inter, who, char) = 0 ;
	    }
	}
    }
  arrayDestroy(defSup) ; 
  arrayDestroy(marSup) ; 
  arrayDestroy(inter) ;
  return (100 * cost) / costmax ;
} /* dupCostED */
#endif /* !NON_GRAPHIC */

/***************************************************************/
/* attention : Pour le moment on ne tient pas compte des brokenLink
   (sur les paires que l'on sait assemblees) et une def ne contenant
   que des NO (rare ?) compte pour cost -1 et non pas 0 */

#ifndef NON_GRAPHIC
static int dupCost(DEFCPT look)
{ if (look->method == M_EDWARD)
    return dupCostED(look) ;
  if (!(look->whatDis & F_E_YN) || (look->whatDis & 1)) 
    return dupCostYes(look) ;   /* only YES data or unknown outside */
  return dupCostNo(look) ;
} /* dupCost */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static MENUOPT defMapGraphMenu[] = {
  { graphDestroy, "Quit"},
  { help, "Help"},
  { graphPrint, "Print"},
  { NULL, NULL  }} ;

/***************************************************************/

static void defMapDisplay(void)
{
  BUMP bp = bumpCreate(100, 0) ;
  int i, j, imax, line = 4, *ip, cost ;
  char *cp ;
  int x1 ;
  float y1 ;
  KEYSET *ksp ;
  KEY *keyp, kludge ;
  Array defOrder = 0, marOrder = 0 ;
  DEFCPTGET("defMapDisplay") ;

  messStatus("please wait for Display") ;
  TESCOLO ;
  TESLIN ;

  if (look->method == M_FINGPR) /* Finger Print ; dessin du gel */
    { fpClearDisplay() ;
      i = arrayMax(look->colOrder) ;
      ip = arrp(look->colOrder, 0, int) - 1 ;
      while(ip++, i--)
	fpDisplay(keySet(look->def, *ip)) ;
    }

  defOrder = defCptTranspMat(look->colOrder) ;
  if (!(look->method == M_EDWARD) && !(look->mapStatus & F_SOMAR))
    defCptSortMarkers() ;
  marOrder = defCptTranspMat(look->linOrder) ;
  if (graphActivate(look->defMapGraph))
    graphPop() ;
  else
    { look->defMapGraph = graphCreate(TEXT_FULL_SCROLL, "Interval/Locus Map", 0.0,0.25,0.4,0.65) ;
      graphMenu(defMapGraphMenu) ;
      graphHelp("Deficiency_map") ;
    }
  graphClear() ;

  if (look->selectedMap)
    graphText(messprintf("Map : %s", name(look->selectedMap)), 1., .2) ;

  cost = dupCost(look) ;
  graphText(messprintf("Cost : %d", cost), 1, 1) ;
  imax = look->nd ;
  ip = arrp(look->colOrder, 0, int) ;
  for (i = 0 ; i < imax ; ip++, i++)
    { cp = name(keySet(look->def, *ip)) ;
      x1 = 0 ; y1 = 2 * i ;
      bumpItem(bp, 1, strlen(cp) + 2.0, &x1, &y1) ;
      graphText(cp, 2 * i + 4, line + x1) ;
    }
  line += bumpMax(bp) ;
  graphBoxStart() ;
  line += 2 ;
  imax = look->nm ;
  x1 = look->nd ;
  ip = arrp(look->linOrder, 0, int) ;
  for (i = 0 ; i < imax ; ip++, i++)
    {
      if (class(keySet(look->mar, *ip)))
	graphText(name(keySet(look->mar, *ip)), 2 * x1 + 8, i + line) ;
      else if (look->method == M_ASSEMB)
	graphText(dnaAlignDecodeOligo(keySet(look->mar, *ip)), 2 * x1 + 8, i + line) ;
    }
  i = look->nd ;                  /* dessin */
  ip = arrp(defOrder, 0, int) - 1 ;
  ksp = arrp(look->marInDef, 0, KEYSET) - 1 ;
  while(ip++, ksp++, i--)
    { j = keySetMax(*ksp) ;
      keyp = arrp(*ksp, 0, KEY) - 1 ;
      if (keySet(*ksp, 0) & F_V_FLAG)
	kludge = 1 ;
      else
	kludge = 0 ;
      while(keyp++, j--)
	{ x1 = arr(marOrder, *keyp & F_WHO, int) ;
	  switch((*keyp & F_FLAG) | kludge)
	    {
	    case 0: case F_YES:
	      graphText("#", 2 *(*ip) + 4, line + .2 + x1) ;
	      break ;
	    case F_NO:
	      if ((look->whatDis & 7) != 1)
		graphText("-", 2 *(*ip) + 4, line + x1) ;
	      break ;
	    default:
	      switch(*keyp & F_PEME)
		{
		case 0: case F_YES:
		  graphText("+", 2 *(*ip) + 4, line - .1 + x1) ;
		  break ;
		case F_NO:
		  if ((look->whatDis & 7) != 1)
		    graphText("_", 2 *(*ip) + 4, line - .2 + x1) ;
		  break ; 
		}
	      break ;
	    }
	}
    }
  graphBoxEnd() ;
  look->mapStatus |= F_SHOWM ;
  bumpDestroy(bp) ;
  arrayDestroy(defOrder) ;
  arrayDestroy(marOrder) ;
  graphTextBounds(2 * look->nd + 48, look->nm + line + 5) ;
  graphRedraw() ;
  defcomputeDraw (look, 0) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/
#ifndef NON_GRAPHIC
static void defCptHisto(char *text1, char *text2)
{ Array h1, h2, h3, h4, h5, defOrder ;
  int i, j, k, n, max, *ip ;
  KEYSET *ksp ;
  KEY *keyp ;
  DEFCPTGET("defCptHisto") ;

  i = look->nd ;              /* Histo des nombres de marqueurs par def */
  h1 = arrayCreate(i, int) ;
  h3 = arrayCreate(15, int) ; /* Histo des nombres de def par nombre de marqueur qu'ils contiennent */
  ip = arrayp(h1, i - 1, int) ;
  ksp = arrp(look->marInDef, i - 1, KEYSET) + 1 ;
  while(ksp--, i--)
    { j = keySetMax(*ksp) ;
      k = 0 ;
      keyp = arrp(*ksp, 0, KEY) - 1 ;
      while(keyp++, j--)
	if (!(*keyp & F_NO))
	  k++ ;
      *ip-- = k ;
      array(h3, k, int)++ ;
    }
  i = look->nm ;              /* Histo des nombres de defs par marqueur */
  h2 = arrayCreate(i, int) ;
  ip = arrayp(h2, i - 1, int) ;
  ksp = arrp(look->defInMar, i - 1, KEYSET) + 1 ;
  while(ksp--, i--)
    { j = keySetMax(*ksp) ;
      k = 0 ;
      keyp = arrp(*ksp, 0, KEY) - 1 ;
      while(keyp++, j--)
	if (!(*keyp & F_NO))
	  k++ ;
      *ip-- = k ;
    }
  TESCOLO ;
  h4 = arrayCreate(25, int) ;
  h5 = arrayCreate(25, int) ;
  defOrder = defCptTranspMat(look->colOrder) ;
  i = look->nm ;
  ksp = arrp(look->defInMar, 0, KEYSET) - 1 ;
  while(ksp++, i--)
    { k = look->nd ;
      max = - 1 ;
      j = keySetMax(*ksp) ;
      keyp = arrp(*ksp, j - 1, KEY) + 1 ;
      while(keyp--, j--)
	{ if (!(*keyp & F_NO))
	    {
	      if ((n = arr(defOrder, *keyp & F_WHO, int)) < k)
		k = n ;
	      if (n > max)
		max = n ;
	    }
	}
      array(h4, max - k, int)++ ;
      j = keySetMax(*ksp) ;
      while(keyp++, j--)
	if (!(*keyp & F_NO))
	  array(h5, arr(defOrder, *keyp & F_WHO, int) - k, int)++ ;
    }
  arrayDestroy(defOrder) ;
  plotHisto(messprintf("Number of %s per %s", text2, text1), h1) ;
  plotHisto(messprintf("Number of %s per %s", text1, text2), h2) ;
  plotHisto(messprintf("Number of %s per number of %s inside", text1, text2), h3) ;
  plotHisto(messprintf("Number of %s per interval between first & last %s", text2, text1), h4) ;
  plotHisto(messprintf("Number of %s at x %s from their first appearance", text2, text1), h5) ;
  look->mapStatus |= F_HIST ;
  defcomputeDraw (look, 0) ;
  /* note that h1, h2, h3, h4, and h5 will 
     be destroyed later by plotHisto */
}

/***************************************************************/

static void defCptRmHisto(void)
{
  DEFCPTGET("defCptRmHisto") ;
  
  plotHistoRemove() ;
  look->mapStatus &= ~F_HIST ;
  defcomputeDraw (look, 0) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/
/*******************   Filter   Functions   ********************/
/***************************************************************/

#ifndef NON_GRAPHIC
static int myArrayOrder(const void *a, const void *b)
{ return *(const int *)a - *(const int *)b ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
#define myabs(__x) ((__x) < 0 ? (-(__x)) : (__x))

static void defCptDoSuppChim(DEFCPT look)
{ KEYSET *ksp, *ksq, ks = 0 ;
  KEY *keyp, *keyq ;
  int i, j, k, n, m, max = look->nd, max2, nb = 0, *ip ;
  Array h1, h2, defOrder ;

  h1 = arrayCreate(look->nd, int) ;
  h2 = arrayCreate(look->nd, int) ;
  array(h2, look->nd - 1, int) = 0 ;
  defOrder = defCptTranspMat(look->colOrder) ;
  ksp = arrp(look->marInDef, 0, KEYSET) ;
  for (i = 0 ; i < look->nd ; ksp++, i++)
    { j = keySetMax(*ksp) ;
      keyp = arrp(*ksp, 0, KEY) - 1 ;
      k = 0 ;
      ks = keySetReCreate(ks) ;
      while(keyp++, j--)
	{ ksq = arrp(look->defInMar, *keyp & F_WHO, KEYSET) ;
	  n = keySetMax(*ksq) ;
	  keyq = arrp(*ksq, 0, KEY) - 1 ;
	  while(keyq++, n--)
	    keySet(ks, k++) = *keyq & F_WHO ;
	}
      keySetSort(ks) ;
      keySetCompress(ks) ;
      m = k = keySetMax(ks) ;
      j = 0 ;
      keyq = arrp(ks, 0, KEY) - 1 ;
      while(keyq++, k--)
	{
	  if ((n = myabs(arr(defOrder, i, int) - arr(defOrder, *keyq, int))) > 20)
	    array(h1, n, int)++ ;
	  j += n ;
	}
      arr(h2, i, int) = j / m ;
    }
  arraySort(h2, myArrayOrder) ;
  plotHisto("nouveau", h1) ;
  plotHisto("nouveau - bis", h2) ;
  arrayDestroy(defOrder) ;
  if (look->manEntry)
    { if (messPrompt("Limit Value", "", "i"))
	{ freeint(&max2) ;
	  if (max2 < max && max2 > 10)
	    max = max2 ;
	}
    }
  else max = look->dlimit ;
  ip = arrp(h2, 0, int) - 1 ;
  i = look->nd ;
  ksp = arrp(look->marInDef, 0, KEYSET) - 1 ;
  while(ip++, ksp++, i--)
    if (*ip >= max)
      { nb++ ;
	keySet(*ksp, 0) |= F_NON_DIAG ;
      }
  if (diffaction)
    catText(diffaction, messprintf("filter %d // -> supp %d seq\n", max, nb)) ;
  look->mapStatus = 1 ;
  look->mapStatus |= F_HIST ;
  arrayDestroy(h1) ;
  arrayDestroy(h2) ;
  defcomputeDraw (look, messprintf("I suppress %d sequences", nb)) ;
}
#endif /* !NON_GRAPHIC */

#ifndef NON_GRAPHIC
static void defCptSuppChim(void)
{ DEFCPTGET ("") ;
 defCptDoSuppChim(look) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptSupOneFlag(DEFCPT look, KEY flag)
{ int i, j ;
  KEYSET *ksp ;
  KEY unflag = ~flag, *keyp ;

  i = look->nd ;
  ksp = arrp(look->marInDef, 0, KEYSET) - 1 ;
  while(ksp++, i--)
    keySet(*ksp, 0) &= unflag ;
  i = look->nm ;
  ksp = arrp(look->defInMar, 0, KEYSET) - 1 ;
  while(ksp++, i--)
    { keySet(*ksp, 0) &= unflag ;
      if ((look->whatDis & F_E_TFLAG) && (flag & F_NON_DIAG))
	{ j = keySetMax(*ksp) - 1 ;
	  keyp = arrp(*ksp, 0, KEY) ;
	  while(keyp++, j--)
	    *keyp &= unflag ;
	}
    }
  if ((look->whatDis & F_E_TFLAG) && (flag & F_NON_DIAG))
    look->whatDis &= ~F_E_TFLAG ;
  look->mapStatus &= ~(F_SHOWM | F_SORT | F_SOMAR) ;
  defcomputeDraw (look, 0) ;
} /* defCptSupOneFlag */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptSupSing(void)
{ int i ;
  KEYSET *ksp ;
  DEFCPTGET("defCptSupSing") ;

  i = look->nd ;
  ksp = arrp(look->marInDef, 0, KEYSET) - 1 ;
  while(ksp++, i--)
    if (keySetMax(*ksp) < 2)
      keySet(*ksp, 0) |= F_SINGLE ;
  look->whatDis |= F_E_ZFLAG ;
  look->mapStatus &= ~(F_SHOWM | F_SORT | F_SOMAR) ;
  defcomputeDraw (look, 0) ;
} /* defCptSupSing */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defTreeSupr(void)
{
  static float xx = .95 ;
  float x1 ;
  int i,j, nn = 0 , n , n0 = 0 ;
  Array nb = 0 ;
  double my , cur, cumul, u0  ;
  KEYSET *ksp ;
  KEY *keyp ;
  DEFCPTGET("defTreeSupr") ;

  if (!messPrompt
      ("Confidence Level (.5 ... 1)", messprintf("%f",xx), "f"))
    return ;
  freefloat(&x1);
  if (x1 > 1 || x1 < .5 )
    { messout("Please, choose a value (%f) between 0.5 and 1", x1) ;
      return ;
    }
  nb = arrayCreate(look->nm, int) ;
  xx = x1 ;
  i = look->nm ;  /* compte les def par marqueur */
  ksp = arrp(look->defInMar, 0, KEYSET) - 1 ;
  while(ksp++, i--)
    { n = 0 ;
      j = keySetMax(*ksp) ;
      keyp = arrp(*ksp, 0, KEY) ;
      if (*keyp & F_SUP_OVER)
	*keyp &= ~F_SUP_OVER ;
      while (j--)
	if (!(*keyp++ & F_NO))
	  n++, nn++ ;
      array(nb, i, int) = n ;
    }
  my = nn / (double)look->nm ;
  cumul = u0 = exp(-my) ; 
  cur = 1 ;
  while(cumul < xx)
    { cur = cur * my / ++n0 ;
      cumul += cur * u0 ;
    }
  i = look->nm ;  /* compte les def par marqueur */
  while(i--)
    if (arr(nb, i, int) > n0)
      keySet(arr(look->defInMar, i, KEYSET), 0) |= F_SUP_OVER ;
  arrayDestroy(nb) ;
  look->whatDis |= F_E_ZFLAG ;
  look->mapStatus &= ~(F_SHOWM | F_SORT | F_SOMAR) ;
  defcomputeDraw (look, 0) ;
} /* defTreeSupr */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defTreeNonDiag(void)
{ int i, j = 0, min, nmax ;
  Array defOrder ;
  KEYSET *ksp ;
  KEY *keyp ;
  DEFCPTGET("non diag") ;

  if (!messPrompt
      ("Please, choose a width", messprintf("%i",j), "i"))
    return ;
  freeint(&nmax) ;
  TESCOLO ;
  defOrder = defCptTranspMat(look->colOrder) ;
  i = look->nm ;
  ksp = arrp(look->defInMar, 0, KEYSET) - 1 ;
  while(ksp++, i--)
    { j = keySetMax(*ksp) ;
      keyp = arrp(*ksp, 0, KEY) - 1 ;
      min = look->nd ;
      while(keyp++, j--)
	{ if (*keyp & F_NON_DIAG)
	    *keyp &= ~F_NON_DIAG ;
	  if (!(*keyp & F_NO) && arr(defOrder, *keyp & F_WHO, int) < min)
	    min = arr(defOrder, *keyp & F_WHO, int) ;
	}
      j = keySetMax(*ksp) ;
      while(keyp--, j--)
	if (!(*keyp & F_NO) && (arr(defOrder, *keyp & F_WHO, int) >= min + nmax))
	  *keyp |= F_NON_DIAG ;
    }
  look->whatDis |= F_E_TFLAG ;
  look->mapStatus &= ~(F_SHOWM | F_SORT | F_SOMAR) ;
  defcomputeDraw (look, 0) ;
} /* defTreeNonDiag */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defSupFlag(void)
{ DEFCPTGET ("defSupFlag") ;

  if (look->whatDis & F_E_FLAG)
    { defCptSupOneFlag(look, F_V_FLAG) ;
      look->whatDis &= ~F_E_FLAG ;
      look->mapStatus &= ~(F_SHOWM | F_SORT | F_SOMAR) ;
    }
  defcomputeDraw (look, 0) ;
} /* defSupFlag */
#endif /* !NON_GRAPHIC */

/***************************************************************/
/**************************    Divers    ***********************/
/***************************************************************/

#ifndef NON_GRAPHIC
static int gmapOrder(const void *a, const void *b)
{ float x ;

  x = ((const GMAPST *) a)->p - ((const GMAPST *) b)->p ;
  return
    x < 0 ? - 1 : (x == 0 ? 0 : 1) ;
} /* gmapOrder */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptMapStan(DEFCPT look)
{ Array gmapOrdered ;
  int i, j, *ip ;
  GMAPST *gmsp ;

  gmapOrdered = arrayCopy(look->gmap) ;
  arraySort(gmapOrdered, gmapOrder) ;
  i = look->nm ;
  look->linOrder = arrayCreate(i, int) ;
  gmsp = arrp(gmapOrdered, i - 1, GMAPST) + 1 ;
  ip = arrayp(look->linOrder, i - 1, int) ;
  while(gmsp--, i--)
    {
      if (keySetFind(look->mar, gmsp->gene, &j))
	*ip-- = j ;
      else
	messcrash("defCptMapStan, gene : %s", name(gmsp->gene)) ;
    }
  arrayDestroy(gmapOrdered) ;
} /* defCptMapStan */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptGetPos(DEFCPT look)
{ int i = arrayMax(look->mar) ;
  float p ;
  OBJ obj ;
  GMAPST *gmsp ;
  KEY *keyp ;

  look->gmap = arrayCreate(i, GMAPST) ;
  gmsp = arrayp(look->gmap, i - 1, GMAPST) + 1 ;
  keyp = arrp(look->mar, i - 1, KEY) + 1 ;
  while(gmsp--, keyp--, i--)
    { gmsp->gene = *keyp ;
      obj = bsCreate(*keyp) ;
      if (obj)
	{ if (bsFindKey(obj, _Map, look->selectedMap))
	    { gmsp->ch = look->selectedMap ;
	      bsPushObj(obj) ;
	      if (bsGetData(obj, _Position, _Float, &p))
		gmsp->p = p ;
	    }
	  bsDestroy(obj) ;
	}
    }
  if (look->nm) defCptMapStan(look) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

static int defCptDoTree(DEFCPT look)
{ float lambda, mu ;
  TREE_DEF *seg ;
  int i ;
  char timeBuf[25] ;


  if (!look->def || keySetMax (look->def) < 2)
    { fprintf (defCptOut, " // Nothing to sort") ;
      return 0 ;
    }

  arrayDestroy (look->colOrder) ;
  arrayDestroy (look->linOrder) ;
  TESCOLO ;
  TESLIN ;

  if (arrayExists(look->maillon))
    { i = arrayMax(look->maillon) ;
      seg = arrp(look->maillon, i - 1, TREE_DEF) + 1 ;
      while(seg--, i--)
	{ keySetDestroy(seg->ks) ;
	  keySetDestroy(seg->x) ;
	}
    }
  arrayDestroy(look->maillon) ;
  if (look->manEntry)
    { freeforcecard(look->distance) ;
      freefloat(&lambda) ;
      if (lambda > 1 || lambda <= 0)
	{ messout("Please, choose a distance limit value (%f) between 0 and 1", lambda) ;
	ici:
	  if (!messPrompt("distance limit (0 ... 1)", messprintf("%f", lambda), "f"))
	    return 0 ;
	  freefloat(&mu) ;
	  if (mu > 1 || mu <= 0)
	    { messout("Please, choose a value (%f) between 0 and 1", mu) ;
	      goto ici ;
	    }
	  lambda = mu ;
	  sprintf(look->distance, "%f", lambda) ;
	}
      look->dlimit = (int)(100 * lambda) ;
    }
  messStatus("I Sort the data") ;
  if ((i = intCptTree(look)))
    {
      look->mapStatus |= (F_SORT | 3) ;
      if (look->maillon && arrayMax(look->maillon))
	look->mapStatus &= ~F_ASSSEG ;
      else arrayDestroy(look->maillon) ;
      look->mapStatus &= ~(F_SHOWM | F_SOMAR) ;
      if (diffaction)
	catText(diffaction, messprintf("sort %d // -> Find %d seg ; %s\n"
				       , look->dlimit, i, timeShow(timeNow(), timeBuf, 25))) ;
    }
  else
    { if (look->display)
	messout("no distances < max") ;
      look->mapStatus = ((look->mapStatus & ~F_BOUT) | 6) ;
    }
#ifndef NON_GRAPHIC
  defcomputeDraw (look, messprintf("I find %d segment", i)) ;
#endif /* !NON_GRAPHIC */
  return i ;
}

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptTree(void)
{ DEFCPTGET("defCptTree") ;

  defCptDoTree(look) ;
} /* defCptTree */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptTreatUnk(void)
{ int i, j ;
  char buf[140] ;
  DEFCPTGET("treat unk") ;

  j = look->whatDis & F_UNK_T ;
  sprintf(buf, "How do you want to consider unknown datas ?\n%s\n%s\n%s\n%s\n%s",
	  "unconsidered (0)", "outside (1)", "1/2 outside (2)",
	  "1/2 in - 1/2 outside (3)", "proportionally (4)") ;
  if (!messPrompt(buf, messprintf("%i", j), "i"))
    return ;
  freeint(&i) ;
  if (i > 4 || i < 0)
    { messout("You are entering a bad value. Try again") ;
      return ;
    }
  look->whatDis &= ~F_UNK_T ;
  look->whatDis |= i ;
  look->mapStatus &= ~(F_SORT | F_SOMAR) ;
  defcomputeDraw (look, 0) ;
} /* defCptTreatUnk */
#endif /* !NON_GRAPHIC */

/***************************************************************/
/********************      Ace - mbly     **********************/
/***************************************************************/
#ifndef NON_GRAPHIC
void defCptAddSeqIn(KEY link, KEY key)
{ DEFCPT look = defCptGetLook (link) ;

  dnaAlignAddSeqIn(look, key) ;
}
#endif  /* !NON_GRAPHIC */

/***************************************************************/
/********************      Def  -  Dup     *********************/
/***************************************************************/

#ifndef NON_GRAPHIC
static void defDatum(KEY gene, Array cc, KEYSET dd)
{ OBJ obj = bsCreate(gene) ;
  int i, n = arrayMax(cc), dummy ;
  KEY def ;
  LINK *ca ;
  Array aa = 0 ;

  if(!obj) return ;
  aa  = arrayReCreate(aa, 30, BSunit) ;
  if( bsFindTag(obj, _Inside) &&
     bsFlatten(obj, 2, aa))
    for (i = 1 ; i < arrayMax(aa) ; i += 2)
      if (def = arr(aa, i, BSunit).k, keySetFind(dd, def, &dummy))
	{ ca = arrayp(cc, n++, LINK) ;
	  ca->a = def ;
	  ca->b = gene ;
	  ca->type = F_YES ;
	  ca->group = 0 ;
	}
  bsGoto(obj, 0) ;
  arrayMax(aa) = 0 ;
  if( bsFindTag(obj,_Outside) &&
     bsFlatten(obj, 2, aa))
    for (i = 1 ; i < arrayMax(aa) ; i += 2)
      if (def = arr(aa, i, BSunit).k, keySetFind(dd, def, &dummy))
	{ ca = arrayp(cc, n++, LINK) ;
	  ca->a = def ;
	  ca->b = gene ;
	  ca->type = F_NO ;
	  ca->group = 0 ;
	} 
  bsDestroy(obj) ;
  arrayDestroy(aa) ;
} /* defDatum */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static BOOL defCptDoParseDB(DEFCPT look, KEYSET map2def, char *cp)
{ Array cc = 0, aa, histo, vertex ;
  LINK *c = 0 ;
  int n, box , line = 15, i, iLocus ;
  KEY map, locus ; 
  OBJ Map , Locus ; 
  KEYSET map2genes = 0, defs, genes ;

  if (map2def)
    map2genes = query(map2def, ">STS") ;
  else if (!lexword2key(cp, &map, _VMap) || !(Map = bsCreate(map)))
    { messout ("Unknown map, I can t proceed") ;
      return FALSE ;
    }
  else
    { map2genes = keySetCreate() ;
      map2def = keySetCreate() ;
      look->selectedMap = map ;
      aa = arrayCreate(100, BSunit) ;
      if (bsFindTag (Map, _Contains))
	bsFlatten(Map, 2, aa) ;
      for (iLocus = 1 ; iLocus < arrayMax(aa) ; iLocus += 2)
	{ locus = arr(aa, iLocus, BSunit).k ;
	  if ((Locus = bsCreate(locus)))
	    { 
	      if (bsFindKey(Locus, _Map, map) &&
		  bsPushObj(Locus))
		{
		  if (bsFindTag(Locus, _Position))
		    keySet(map2genes, keySetMax(map2genes)) = locus ;
		  else if (bsFindTag(Locus, _Left))
		    keySet(map2def, keySetMax(map2def)) = locus ;
		}
	      bsDestroy(Locus) ;
	    }
	}
      arrayDestroy(aa) ;
    }
  keySetSort(map2genes) ;
  keySetCompress(map2genes) ;
  keySetSort(map2def) ;
  keySetCompress(map2def) ;
  if (keySetMax(map2genes) && keySetMax(map2def))
    defcomputeDraw (look, messprintf("Map : %s%s%d %s, %d %s", 
			      look->selectedMap ? name(look->selectedMap) : "",
			      look->selectedMap ? ", " : "",
			      keySetMax(map2def), className(keySet(map2def, 0)),
			      keySetMax(map2genes), className(keySet(map2genes, 0)))) ;
  cc = arrayCreate(50, LINK) ;
  for (i = 0; i < keySetMax(map2genes) ; i++)
    defDatum(keySet(map2genes, i), cc, map2def) ;
  keySetDestroy(map2genes) ;
  keySetDestroy(map2def) ;
  defs  = keySetCreate() ;
  genes  = keySetCreate() ;
  n = arrayMax(cc) ;
  if (n) c = arrp(cc, 0, LINK) - 1 ;
  while(c++, n--)
    { keySetInsert(defs, c->a) ;
      keySetInsert(genes, c->b) ;
    }
  box = graphBoxStart() ;
  if (keySetMax(defs) && keySetMax(genes))
    graphText(messprintf(" %d data,  %d %s tested against %d %s",
			 arrayMax(cc), keySetMax(defs), className(keySet(defs, 0)),
			 keySetMax(genes), className(keySet(genes, 0))), 3, line++) ;
  graphText("I will now assemble these into contigs, please wait",
	    3,line++) ;
  graphBoxDraw(box, BLACK,WHITE) ;
  n = arrayMax(cc) ;
  vertex = arrayCreate(n/2, VERTEX) ;
  look->nbContig = topoConnectedComponents(cc, vertex) ;
  arrayMax(cc) = n ;
  box = graphBoxStart() ;
  graphText(messprintf("found %d contigs", look->nbContig), 2, line++) ;
  graphBoxDraw(box, BLACK,WHITE) ;
  histo = arrayCreate(20,int) ;
  n = arrayMax(vertex) ;
  while(n--)
    array(histo, arr(vertex,n,VERTEX).group, int)++ ;
  plotHisto("Genes or Def per  contig", histo) ;
  look->mapStatus |= F_HIST ;
  look->dataArray = cc ;
  if (!arrayMax(vertex))
    { arrayDestroy(look->dataArray) ;
      lookArrayDestroy(look) ;
      arrayDestroy(vertex) ;
      keySetDestroy(defs) ;
      keySetDestroy(genes) ;
      return FALSE ;
    }
  arrayDestroy(vertex) ;
  keySetDestroy(defs) ;
  keySetDestroy(genes) ;
  return TRUE ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static BOOL defCptParseDB(DEFCPT look)
{ KEYSET bb = 0, ks = 0 ;
  void *dummy ;
  char *cp = 0 ;

  if (look->nbContig)
    messcrash("look pas reinitialise correctement") ;
  if (keySetActive(&bb, &dummy) && keySetMax(bb))
    { ks = query(bb, "CLASS YAC") ;
      if (!keySetMax(ks))
	keySetDestroy(ks) ;
    }
  graphActivate(look->defMapCtlGraph) ;
  graphPop() ;
  if (!ks)
    { if (messPrompt("Give the name of the Map you want to analyse","","w"))
	cp = freeword() ;
      else
	return FALSE ;
    }
  return defCptDoParseDB(look, ks, cp) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptDoGetData(DEFCPT look, int t)
{ 
  int i, j, n, box, line, no = 0 ;
  Array cc ; 
  LINK *c ;
  KEYSET *ksp ;
  KEY *keyp ;

  lookArrayDestroy(look) ;
  cc  = look->dataArray ;
  if(!cc)
    messcrash("defCptGetData receives a null Array") ;
  look->def = keySetCreate() ;
  look->mar = keySetCreate() ;
  n = arrayMax(cc) ;
  while(n--)
    { c = arrp(cc, n, LINK) ;
      if (t == c->group)
	{ keySetInsert(look->def, c->a) ;
	  keySetInsert(look->mar, c->b) ;
	}
    }
  look->nd = keySetMax(look->def) ;
  look->nm = keySetMax(look->mar) ;
  i = look->nd ;
  look->marInDef = arrayCreate(i, KEYSET) ;
  array(look->marInDef, i - 1, KEYSET) = 0 ;
  ksp = arrp(look->marInDef, 0, KEYSET) ;
  while(i--)
    *ksp++ = keySetCreate() ;
  i = look->nm ;
  look->defInMar = arrayCreate(look->nm, KEYSET) ;
  array(look->defInMar, i - 1, KEYSET) = 0 ;
  ksp = arrp(look->defInMar, 0, KEYSET) ;
  while(i--)
    *ksp++ = keySetCreate() ;
  n = arrayMax(cc) ;
  while(n--)
    { c = arrp(cc, n, LINK) ;
      if (t == c->group && c->type)
	{ 
	  if (keySetFind(look->def, c->a, &i) && keySetFind(look->mar, c->b, &j))
	    { ksp = arrp(look->marInDef, i, KEYSET) ;
	      keySet(*ksp, keySetMax(*ksp)) = j | (c->type) ;
	      ksp = arrp(look->defInMar, j, KEYSET) ;
	      keySet(*ksp, keySetMax(*ksp)) = i | (c->type) ;
	      if (c->type == F_NO)
		no = 1 ; /* pour savoir si que des in ou in et out */
	    }
	  else
	    messcrash("Cant find %d -> %s\n%d -> %s", i, name(c->a), j, name(c->b)) ;
	}
    }
  i = look->nd ;
  ksp = arrp(look->marInDef, 0, KEYSET) - 1 ;
  while(ksp++, i--)
    arraySort(*ksp, mardefOrder) ;
  i = look->nm ;
  ksp = arrp(look->defInMar, 0, KEYSET) - 1 ;
  while(ksp++, i--)
    arraySort(*ksp, mardefOrder) ;
  defCptGetPos(look) ;
  look->mapStatus = (1 | (look->mapStatus & F_HIST)) ;
  if (no)
    look->whatDis = (F_E_YN | 1) ;
  else /* que des yes => suppression des flag YES pour shortDistances */
    { i = arrayMax(look->marInDef) ;
      ksp = arrp(look->marInDef, 0, KEYSET) - 1 ;
      while(ksp++, i--)
	{ j = keySetMax(*ksp) ;
	  keyp = arrp(*ksp, 0, KEY) - 1 ;
	  while(keyp++, j--)
	    *keyp &= F_WHO ;
	}
      i = arrayMax(look->defInMar) ;
      ksp = arrp(look->defInMar, 0, KEYSET) - 1 ;
      while(ksp++, i--)
	{ j = keySetMax(*ksp) ;
	  keyp = arrp(*ksp, 0, KEY) - 1 ;
	  while(keyp++, j--)
	    *keyp &= F_WHO ;
	}
      look->whatDis = 1 ; /* que des Yes : unknown out (1) ou rien (0) ? */
    }
  if (look->selectedMap)
    defcomputeDraw (look, messprintf("Map : %s", name(look->selectedMap))) ;
  else defcomputeDraw (look, 0) ;
  line = 14 ;
  box = graphBoxStart () ;
  graphText(messprintf("Contig :  %d", t), 1, line);
  line += 2 ;
  if (look->nm && look->nd)
    graphText(messprintf("%d %s, %d %s",
			 look->nd, className(keySet(look->def, 0)),
			 look->nm, className(keySet(look->mar, 0))), 
	      12, line++) ;
  graphBoxDraw (box, BLACK, WHITE) ;
  look->Line = 8;
  if(!look->nm)
    { lookArrayDestroy(look) ;
      return ;
    }
} /* defCptDoGetData */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptDoExportOrder(DEFCPT look, FILE *f, char buf[81])
{ char *buf2 ;
  int i, j, m, n, x, imax ;
  Array marOrder ;
  KEYSET *ksp ;
  KEY *keyp ;

  switch(look->method)
    { case M_EDWARD:
	imax = look->nd ;
	buf2 = className(keySet(look->def, 0)) ;
	for (i = 0 ; i < imax ; i++)
	  { j = arr(look->colOrder, i, int) ;
	    if (!(keySet(arr(look->marInDef, j, KEYSET), 0) & F_V_FLAG))
	      fprintf(f, "%s %s\nMap %s Position %d\n\n", buf2, name(keySet(look->def, j)), buf, i) ;
	  }
	look->method = M_FIN ;
#ifndef NON_GRAPHIC
	defcomputeDraw (look, 0) ;
#endif /* !NON_GRAPHIC */
	break ;
      default:
	imax = look->nm ;
	buf2 = className(keySet(look->mar, 0)) ;
	for (i = 0 ; i < imax ; i++)
	  { j = arr(look->linOrder, i, int) ;
	    fprintf(f, "%s %s\nMap %s Position %d\n\n", buf2, name(keySet(look->mar, j)), buf, i) ;
	  }
	buf2 = className(keySet(look->def, 0)) ;
	marOrder = defCptTranspMat(look->linOrder) ;
	ksp = arrp(look->marInDef, 0, KEYSET) ;
	for (i = 0 ; i < look->nd ; ksp++, i++)
	  { j = keySetMax(*ksp) ;
	    keyp = arrp(*ksp, 0, KEY) - 1 ;
	    x = look->nm ;
	    m = - 1 ;
	    while(keyp++, j--)
	      { if (!(*keyp & F_NO))
		  { if ((n = arr(marOrder, *keyp & F_WHO, int)) < x)
		      x = n ;
		    if (n > m)
		      m = n ;
		  }
	      }
	    if (m != -1)
	      { fprintf(f, "%s %s\nMap %s Left %d\n", buf2, name(keySet(look->def, i)), buf, x) ;
		fprintf(f, "Map %s Right %d\n\n", buf, m) ;
	      }
	  }
	break ;
      }
} /* defCptDoExportOrder */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptExportOrder(void)
{ FILE *f ;
  KEY myMap, vMap ;
  char *tmpName  = 0 ;
  char *cp, buf[81], filname[FIL_BUFFER_SIZE] ;
  DEFCPTGET("defCptExportOrder") ;

  TESLIN ;
  strcpy(filname, "best.tree") ;
  if (!(messPrompt("Please, give a name for the exported map", "Test", "w")))
    return ;
  cp = freeword() ;
  strncpy(buf, cp, 80) ;
  lexaddkey (buf, &myMap, _VMap) ;
  f = filtmpopen(&tmpName,"w") ;
  if (!f)
    return ;
  defCptDoExportOrder(look, f, buf) ;
  filclose(f) ;
  if ((f = filopen (tmpName, 0, "r")))
    parseFile (f, 0, 0) ;
  filtmpremove (tmpName) ;
  if (lexReClass(myMap, &vMap,_VvMap) &&
      !lexlock (vMap))
    arrayKill (vMap) ;
  display (myMap, 0, 0) ;
} /* defCptExportOrder */
#endif

/***************************************************************/

static int defCptContigOrder(const void *a, const void *b)
{ return
    ((const KEY_INT*)b)->tai - ((const KEY_INT*)a)->tai ;
}

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptAllContigs(DEFCPT look)
{ int i, j = 0 ;
  char *cp, buf[81], *buf2, *tmpName = 0 ;
  KEY link, myMap, vMap ;
  Array contig = 0 ;
  FILE *f ;
  KEY_INT *kip ;

  if (!messPrompt("Please, give a name for The Maps", "Test", "w"))
    { look->mapStatus |= F_NO_CONTIG ;
      defcomputeDraw (look, 0) ;
      return ;
    }
  cp = freeword() ;
  strncpy(buf, cp, 80) ;
  f = filtmpopen(&tmpName, "w") ;
  if (!f)
    return ;
  look->display = FALSE ;
  lexaddkey(messprintf("_Assembly.%s", buf), &link, _VMap) ;
  contig = arrayCreate(look->nbContig, KEY_INT) ;
  arrayp(contig, look->nbContig - 1, KEY_INT)->key = 0 ;
  kip = arrp(contig, 0, KEY_INT) ;
  for (i = 1 ; i <= look->nbContig ; i++)
    { defCptDoGetData(look, i) ;
      defCptDoTree(look) ;
      defCptSortMarkers() ;
      buf2 = messprintf("%s.%d", buf, i) ;
      lexaddkey(buf2, &myMap, _VMap) ;
      kip->key = myMap ;
      (kip++)->tai = look->nm - 1 ;
      defCptDoExportOrder(look, f, buf2) ;
    }
  look->display = TRUE ;
  arraySort(contig, defCptContigOrder) ;
  i = look->nbContig ;
  kip = arrp(contig, 0, KEY_INT) - 1 ;
  while(kip++, i--)
    { if (kip->tai)
	{ fprintf(f, "Map %s\nMap %s Left %d\n", name(kip->key), name(link), j) ;
	  j += kip->tai ;
	  fprintf(f, "Map %s Right %d\n\n", name(link), j) ;
	  j += 20 ;
	}
    }
  filclose(f) ;
  arrayDestroy(contig) ;
  if (filopen(tmpName, 0, "r"))
    parseFile(f, 0, 0) ;
  filtmpremove(tmpName) ;
  if (lexReClass(link, &vMap, _VvMap) && !lexlock(vMap))
    arrayKill(vMap) ;
  display(link, 0, 0) ;
  look->mapStatus |= F_NO_CONTIG ; /* a voir ou il faut retomber */
  defcomputeDraw (look, 0) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptGetData(void)
{ int t = 0 ;
  DEFCPTGET("defCptGetData") ;

  look->mapStatus |= ~F_NO_CONTIG ;
  messStatus("Get Data") ;
  if (!look->dataArray && !defCptParseDB(look))
    { messout("I have no data to proceed") ;
      defcomputeDraw (look, 0) ;
      return ;
    }
  if (look->nbContig == 1)
    t = 1 ;
  else if(!messPrompt("Def contig number ? (0 for All)","1","i") || 
	  !freeint(&t) || t > look->nbContig)
    { look->mapStatus |= F_NO_CONTIG ;
      if (t > look->nbContig)
	messout(messprintf("there are only %d contigs", look->nbContig)) ;
      if (look->selectedMap)
	defcomputeDraw (look, messprintf("Map : %s", name(look->selectedMap))) ;
      else defcomputeDraw (look, 0) ;
      return ;
    }
  if (t)
    defCptDoGetData (look, t) ;
  else
    defCptAllContigs(look) ;
} /* defCptGetData */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptDefDupHisto(void)
{ defCptHisto("Deficiency", "marker") ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/
/********************    Finger   Print    *********************/
/***************************************************************/
#ifndef NON_GRAPHIC
static void defCptGetClones(void) /* provisoirement pas de version tace de fp */
{
  extern void pmapCptGetData(Array marInDef, Array defInMar,
			     KEYSET clones, KEYSET bands) ;
  DEFCPTGET("defCptGetClones") ;

  messStatus("GetClones") ;
  look->def = keySetCreate() ;
  look->mar = keySetCreate() ;
  look->marInDef = arrayCreate(100, KEYSET) ;
  look->defInMar = arrayCreate(100, KEYSET) ;
  pmapCptGetData(look->marInDef, look->defInMar, look->def, look->mar) ;
  look->nm = keySetMax(look->mar) ;
  if (!look->nm)
    { lookArrayDestroy(look) ;
      return ;
    }
  look->nd = keySetMax(look->def) ;
  defCptGetPos(look) ;
  look->mapStatus = 1 ;
  look->whatDis = 1 ;
  defcomputeDraw (look, messprintf("Found %d Clones and %d Bands", look->nd, look->nm)) ;
} /* defCptGetClones */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptFingPrHisto(void)
{ defCptHisto("Clone", "Band") ;
}
#endif /* pas de version tace pour le finger print */
/***************************************************************/
/********************  Sequence Alignment  *********************/
/***************************************************************/

static void defCptNbOligoChoice (DEFCPT look)
{ int i = 0 ;

  if ((i = keySetMax (look->def)))
    i = 1000 / i ;
  if (i < 2) i = 2 ;
  if (i > 60) i = 60 ;
  look->whatDis = (unsigned char)i ;
  sprintf (look->nboligo, "%d", i) ;
}

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptDoInsertBad(DEFCPT look)
{ dnaAlignInsertBad(look) ;
  if (diffaction)
    catText(diffaction, messprintf("insert // -> %d seq\n", keySetMax(look->def))) ;
  defcomputeDraw (look, messprintf("Found %d sequences", keySetMax(look->def))) ;
  defCptNbOligoChoice (look) ;
}
#endif /* !NON_GRAPHIC */

#ifndef NON_GRAPHIC
static void defCptInsertBad(void)
{ DEFCPTGET("insert Bad") ;

  defCptDoInsertBad(look) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

static void defCptDoFixSegConsensus(DEFCPT look, BOOL new)
{ KEYSET maille ;
  char timeBuf[25] ;

  messStatus("Fix consensus") ;
  if (new || !new)  /* ulrich's fix */
    {
      maille = look->def ;
      look->def = keySetCreate() ;
      dnaAlignFixSegConsensus(look, maille) ;
      keySetDestroy(maille) ;
    }
  else
    {
    }
  if (diffaction)
    catText(diffaction, messprintf("fix // %s\n", timeShow(timeNow(), timeBuf, 25))) ;
  look->mapStatus &= ~F_BOUT ;
  look->mapStatus |= (F_ASSSEG | 4 | F_IS_FIXED) ;
#ifndef NON_GRAPHIC
  defcomputeDraw (look, 0) ;
#endif /* !NON_GRAPHIC */
} /* defCptDoFixSegConsensus */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptFixSegConsensus(void)
{ DEFCPTGET("fix seg") ;

  defCptDoFixSegConsensus(look, FALSE) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/
/***************************************************************/

static void defCptDoTrainNNonKeySet (KEYSET aa)
{ int i ; 
  KEY key ;
  char timeBuf[25] ;

  messStatus("AutoEdit") ;
  if (!aa)
    return ;
  for (i = 0 ; i < keySetMax(aa) ; i++)
    { key = keySet(aa, i) ;
      if (dnaReClass (key, &key))
	nnContigTrain (key);
    }
  nnContigTrain (0) ; /* close export file */
  if (diffaction)
    catText(diffaction, messprintf("TrainNN // TrainNN contigs ; %s\n",
				   timeShow(timeNow(), timeBuf, 25))) ;
}

/***************************************************************/
/***************************************************************/

static int defCptDoPatchKeySet (KEYSET aa)
{ int i, n, nn = 0 ; 
  KEY key ;
  char timeBuf[25] ;

  messStatus("AutoEdit") ;
  if (!aa) 
    return 0 ;
  for (i = 0 ; i < keySetMax(aa) ; i++)
    { 
      key = keySet(aa, i) ;
      n = 0 ;
      if (dnaReClass (key, &key) && 
	  !baseCallRedoBaseCall (key, &n)) /* FALSE si F4 */
	{ nn += n ; break ; }
      nn += n ;
    }
  if (diffaction)
    catText(diffaction, messprintf("autoedit // autoedit contigs ; %s\n", 
				   timeShow(timeNow(), timeBuf, 25))) ;
  return nn ;
}

/***************************************************************/
#ifdef JUNK
static int defCptDoFlagRnaEditingKeySet (KEYSET aa)
{ int i, n, nn = 0 ; 
  KEY key ;
  char timeBuf[25] ;

  messStatus("Flag RNA editing") ;
  if (!aa)
    return 0 ;
  for (i = 0 ; i < keySetMax(aa) ; i++)
    { 
      key = keySet(aa, i) ;
      n = 0 ;
      if (dnaReClass (key, &key) && 
	  !baseCallRedoBaseCall (key, &n)) /* FALSE si F4 */
	{ nn += n ; break ; }
      nn += n ;
    }
  if (diffaction)
    catText(diffaction, messprintf("autoedit // autoedit contigs ; %s\n", 
				   timeShow(timeNow(), timeBuf, 25))) ;
  return nn ;
}
#endif
/***************************************************************/

static int defCptDoBaseCall (KEYSET aa)
{ int i, n, nn = 0 ; 
  KEY key ;

  messStatus("defCptDoBaseCall") ;
  if (!aa)
    return 0 ;
  for (i = 0 ; i < keySetMax(aa) ; i++)
    { 
      key = keySet(aa, i) ;
      n = 0 ;
      if (dnaReClass (key, &key) &&
	  !baseCallRedoBaseCall (key, &n)) /* FALSE si F4 */
	  { nn += n ; break ; }
      nn += n ;
    }
  return nn ;
}

/***************************************************************/
#ifndef NON_GRAPHIC
/***************************************************************/

static void defCptPatchKeySet (void)
/* graphic */
{ DEFCPTGET ("defCptPatchKeySet") ;

  defCptDoPatchKeySet (look->def) ;
}

/***************************************************************/

#ifdef THIS_CODE_IS_NEVER_USED

static void defCptPatchActiveKeySet (void)
/* graphic */
{
  KEYSET aa ;
  int nn ;
  void *dummy ;
  DEFCPTGET("defCptPatchActiveKeySet") ;

  if (!isWriteAccess ())	
    { messout("Sorry, you do not have write access") ;
      return ;
    }

  if (!keySetActive(&aa, &dummy))
    { messout("First select a keyset containing sequences") ;
      return ;
    }
  nn = defCptDoPatchKeySet (aa) ;
  defcomputeDraw (look, messprintf ("Performed %d auto-editions", nn)) ;
}


#endif /* THIS_CODE_IS_NEVER_USED */

/***************************************************************/

static void defCptChoosePaire(KEY key)
/* graphic */
{ DEFCPTGET("Choose Paire") ;

  array(look->dataArray, arrayMax(look->dataArray), KEY) = key ;
  keySetSort(look->dataArray) ;
  keySetCompress(look->dataArray) ;
  if (arrayMax(look->dataArray) != 2)
    { displayRepeatBlock () ;
      return ;
    }
  if (arrayMax (look->dataArray) != 2)
    messcrash("Impossible Max pour dataArray") ;
  dnaAlignAsmbPaire (look, arr(look->dataArray, 0, KEY),
		    arr(look->dataArray, 1, KEY)) ;
  arrayDestroy (look->dataArray) ;
  defcomputeDraw(look, 0) ;
}

/***************************************************************/

static void defCptPaireDestroy(void)
/* graphic */
{ DEFCPTGET ("Paires Destroy") ;

  displayUnBlock () ;
  arrayDestroy(look->dataArray) ;
}

/***************************************************************/

static void defCptAsmbPaireActive (void)
/* graphic */
{ DEFCPTGET ("Asmb Paire") ;

  if (!isWriteAccess ())	
    { messout("Sorry, you do not have write access") ;
      return ;
    }
  if (look->dataArray)
    messcrash("look->dataArray should not exit") ;
  look->dataArray = arrayCreate(2, KEY) ;
  graphRegister(MESSAGE_DESTROY, defCptPaireDestroy) ;
  displayBlock (defCptChoosePaire, "Pick on sequence") ;
}

/***************************************************************/

void defCptMakeSuperLink (void)
/* graphic */
{ KEY link ;
  KEYSET ks ;

  if (!isWriteAccess ())	
    { messout("Sorry, you do not have write access") ;
      return ;
    }
  if (!keySetActive (&ks, 0))
    { messout ("First select a keySet containing sequences") ;
      return ;
    }
  link = dnaAlignDoMakeSuperLink (ks, 0) ;
  display(link, 0, 0) ;
}
/***************************************************************/

static void defCptFixActKeyset (void)
/* graphic */
{ DEFCPTGET("fix keyset") ;

  dnaAlignFixActKeyset () ;
  look->mapStatus |= F_IS_FIXED ;
  defcomputeDraw (look, 0) ;
}

/***************************************************************/
#endif /* !NON_GRAPHIC */
/***************************************************************/

static void defCptDoExtendReadsExt (DEFCPT look)
{
  char timeBuf[25] ;

  messStatus("Extend Reads (ext)") ;
  abiFixExtend (look->link, look->def) ;
  look->mapStatus = (F_ASSSEG | 5 | (look->mapStatus & F_HIST)) ;
/* defCptDoFixSegConsensus (look) ; */
  if (diffaction)
    catText(diffaction, messprintf("extend  // %s ; extend reads ext\n", timeShow(timeNow(), timeBuf, 25))) ;
#ifndef NON_GRAPHIC
  defcomputeDraw (look, 0) ;
#endif /* !NON_GRAPHIC */
} /* defCptDoExtendReadsExt */

#ifndef NON_GRAPHIC
static void defCptExtendReadsExt (void)
{ DEFCPTGET("extend read") ;

  defCptDoExtendReadsExt (look) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

static void defCptDoGetOligos(DEFCPT look)
{ int i ;
  KEYSET *ksp ;
  char timeBuf[25] ;

  if (look->manEntry)
    { freeforcecard(look->nboligo) ;
      freeint(&i) ;
      if (i > 0 && i < 256)
	look->whatDis = (unsigned char)i ;
    }
  if (!look->def)
    look->def = keySetCreate () ;
  look->mar = keySetReCreate(look->mar) ;
  if (arrayExists(look->marInDef))
    { if ((i = arrayMax(look->marInDef)))
	{ ksp = arrp(look->marInDef, 0, KEYSET) - 1 ;
	  while(ksp++, i--)
	    keySetDestroy(*ksp) ;
	}
    }
  look->marInDef = arrayReCreate(look->marInDef, 100, KEYSET) ;
  if (arrayExists(look->defInMar))
    { if ((i = arrayMax(look->defInMar)))
	{ ksp = arrp(look->defInMar, 0, KEYSET) - 1 ;
	  while(ksp++, i--)
	    keySetDestroy(*ksp) ;
	}
    }
  look->defInMar = arrayReCreate(look->defInMar, 100, KEYSET) ;
  messStatus("Get Sequences") ;
  dnaAlignGetData(look) ;
  if ((look->nd = keySetMax(look->def)) < 2)
    { messout("GetOligo: Only %d sequence, i can't get", look->nd) ;
      if (!look->mapStatus) /* i e c'est le premier tour et l'utilisateur */
	{ keySetDestroy(look->def) ; /* a charge 1 seq. => retour au depart */
	  keySetDestroy(look->mar) ;
	}
#ifndef NON_GRAPHIC
      defcomputeDraw (look, 0) ;
#endif /* !NON_GRAPHIC */
      return ;
    }
  if(!(look->nm = keySetMax(look->mar)))
    { messout("I haven't found any oligo") ; /* idem */
      if (!look->mapStatus) /* attention look->mapStatus ne doit jamais etre nul en dehors du depart */
	{ keySetDestroy(look->def) ;
	  keySetDestroy(look->mar) ;
	}
#ifndef NON_GRAPHIC
      defcomputeDraw (look, 0) ;
#endif /* !NON_GRAPHIC */
      return ;
    }
  if (diffaction)
    catText(diffaction, messprintf("get %d // %s\n", look->whatDis, timeShow(timeNow(), timeBuf, 25))) ; /* pour faire le fichier des actions */
  arrayDestroy(look->colOrder) ;
  arrayDestroy(look->linOrder) ;
  TESCOLO ;
  TESLIN ;
  look->mapStatus &= ~F_BOUT ;
  look->mapStatus |= (F_HIST | 1) ;
#ifndef NON_GRAPHIC
  defcomputeDraw (look, 0) ;
#endif /* !NON_GRAPHIC */
} /* defCptDoGetOligos */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptGetOligos(void)
{ DEFCPTGET("get sequences") ;

  defCptDoGetOligos(look) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptGetMoreSequences(void)
{
  int i, j ;
  char buf[1000] ;
  DEFCPTGET("more oligo") ;

  sprintf(buf, "How many oligos do you want") ;
  j = (look->whatDis + 2) ;
 ici:
  if (messPrompt(buf, messprintf("%d", j), "i"))
    { freeint(&i) ;
      if (i < 0 || i > 256)
	{ sprintf(buf, "Bad value (%d) ; please try again", i) ;
	  goto ici ;
	}
      look->whatDis = (unsigned char)i ;
      sprintf(look->nboligo, "%d", i) ;
    }
  else return ;
  look->manEntry = FALSE ; /* kludge pour imposer un nb oligos a GetOligos */
  defCptDoGetOligos(look) ;
  look->manEntry = TRUE ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/
/*  probablement inutile
static void defCptMakeOneSegment(void)
{ int i, j, k = 0, imax, jmax ;
  TREE_DEF maille, *ksmp ;

  if (arrayMax(look->maillon) == 1)
    { messout("This work is already done") ;
      return ;
    }
  imax = arrayMax(look->maillon) ;
  maille.ks = keySetCreate() ;
  maille.x = keySetCreate() ;
  ksmp = arrp(look->maillon, 0, TREE_DEF) ;
  for (i = 0 ; i < imax ; ksmp++, i++)
    { jmax = keySetMax(ksmp->ks) ;
      for (j = 0 ; j < jmax ; j++)
	{ keySet(maille.ks, k) = keySet(ksmp->ks, j) ;
	  keySet(maille.x, k++) = keySet(ksmp->x, j) ;
	}
      keySetDestroy(ksmp->ks) ;
      keySetDestroy(ksmp->x) ;
    }
  look->maillon = arrayReCreate(look->maillon, 1, TREE_DEF) ;
  array(look->maillon, 0, TREE_DEF) = maille ;
  if (diffaction)
    catText(diffaction, "one\n") ;
}
*/
/***************************************************************/

static BOOL defCptFixErrorRate(DEFCPT look)
{ int mu ;

  if (look->manEntry)
    { freeforcecard(look->choix) ;
      freeint(&mu) ;
      if (mu > 50 || mu < 0)
	{ if (mu < 0)
	    messout("Please, choose a rate (%d) between 0 and 50", mu) ;
	  if (mu > 50)
	    messout("Are you sure to want %d %% error ; please confirm", mu) ;
	ici:
	  if (!messPrompt("error rate", messprintf("%d", mu), "i"))
	    return FALSE ;
	  freeint(&mu) ;
	  if (mu < 0)
	    { messout("Please, choose a value (%d) between 0 and 50", mu) ;
	      goto ici ;
	    }
	  sprintf(look->choix, "%d", mu) ;
	}
      look->taux = mu ;
    }
  return TRUE ;
}

/***************************************************************/

static void defCptEndAssfunc(DEFCPT look, KEYSET ks)
{ look->mapStatus = (F_ASSSEG | 5 | (look->mapStatus & F_HIST)) ;
#ifndef NON_GRAPHIC
  defcomputeDraw (look, messprintf("Found %d contigs, rejected %d Baddies",
				   keySetMax(look->def), look->rejected ?
				   keySetMax (look->rejected) : 0)) ;
#endif /* !NON_GRAPHIC */
  if (keySetMax(look->def))
    keySetDestroy(ks) ;
  else
    { keySetDestroy(look->def) ;
      look->def = ks ;
    }
  keySetDestroy(look->mar) ;
  defCptNbOligoChoice (look) ;
}

/***************************************************************/

static void defCptDoSegments(DEFCPT look)
{ KEYSET defbis = 0 ;
  char timeBuf[25] ;

  messStatus("Assembling") ;
  if (!defCptFixErrorRate(look))
    return ;
  defbis = look->def ;
  look->def = keySetCreate() ;
  if (!look->maillon || !arrayMax (look->maillon))
    { messout ("sorry no assembly") ;
      goto abort ;
    }
  dnaAlignCptSegments(look) ;

  if (diffaction)
    catText(diffaction, messprintf("assemble %d // -> Found %d contigs, rejected %d Baddies ; %s\n",
				   look->taux, keySetMax(look->def),
				   look->rejected ?
				   keySetMax(look->rejected) : 0, 
				   timeShow(timeNow(), timeBuf, 25))) ;
 abort:
  defCptEndAssfunc(look, defbis) ;
}

/***************************************************************/

static Array defCptCollectDummy (KEY join, KEYSET jj, int niveau)
{ KEY key ;
  BSunit *u ;
  OBJ Join = 0 ;
  Array units = 0, new = 0, daughter = 0 ;
  int p1, p2, x1 , x2,  dummy, i, i1, j = 0 ;

  if (dnaReClass(join, &join))
    Join = bsCreate(join) ;
  if (!Join)
    return 0 ;
  units = arrayCreate (90, BSunit) ;
  new = arrayCreate (90, BSunit) ;
  if (bsFindTag (Join, _Previous_contig) &&
      bsFlatten (Join, 3, units))
    { for (i = 0 ; i < arrayMax(units) ; i += 3)
	{ u = arrp(units,i,BSunit) ;
	  key = u[0].k ;
	  p1 = u[1].i ; p2 = u[2].i ;
	  if (keySetFind (jj, key, &dummy))
	    { daughter = defCptCollectDummy (key, jj, niveau + 1) ;
	      if (daughter) for (i1 = 0 ; i1 < arrayMax(daughter) ; i1 += 3)
		{ u = arrp(daughter,i1,BSunit) ;
		  array(new, j++, BSunit).k = u[0].k ;
		  if (p1 < p2)
		    { x1 = p1 + u[1].i - 1 ; x2 = p1 + u[2].i - 1 ; }
		  else
		    { x1 = p1 - u[1].i + 1 ; x2 = p1 - u[2].i + 1 ; }		    
		  array(new, j++, BSunit).i = x1 ;
		  array(new, j++, BSunit).i = x2 ;
		}
	      arrayDestroy (daughter) ;
	    }
	  else
	    { array(new, j++, BSunit).k = key ;
	      array(new, j++, BSunit).i = p1 ;
	      array(new, j++, BSunit).i = p2 ;
	    }
	  
	}
    }
  arrayDestroy (units) ;
  bsDestroy (Join) ;
  if (niveau)
    return new ;
  Join = bsUpdate(join) ;
  if (Join && bsAddArray (Join, _Previous_contig, new, 3))
    bsSave (Join) ;
  else
    bsDestroy (Join) ;
  arrayDestroy (new) ;
  return 0 ;    
}

/***************************************************************/

static int defCptDoJoinDiagonal (DEFCPT look)
{ int i, j = 0, nj = 0, fl, imax ;
  KEY key1, key2, join, link ;
  KEYSET ks = 0, ks1 = 0, ks2 = 0, ks3 = 0 ;
  Array subseq = 0, flag = 0 ;
  OBJ obj = 0 ;
  mytime_t tim ;
  char *cp ;

  messStatus ("Joining") ;
  if (!look || !look->link)
    return 0 ;
  if (class (look->link) == _VSequence)
    cp = dnaAlignNewName (name (look->link)) ;
  else
    cp = dnaAlignNewName (0) ;
  subseq = arrayCreate (20, BSunit) ;
  if (!(obj = bsCreate (look->link)) || !bsFindTag (obj, _Subsequence) ||
      !bsFlatten (obj, 1, subseq) || (arrayMax (subseq) < 2))
    goto abort ;
  bsDestroy (obj) ;
  i = arrayMax (subseq) ;
  flag = arrayCreate (i, int) ;
  array (flag, i-1, int) = 0 ;
  for (i = 1 ; i < arrayMax (subseq) ; i++)
    { key1 = arrp (subseq, i - 1, BSunit)->k ;
      key2 = arrp (subseq, i, BSunit)->k ;
      fl = arr (flag, i - 1, int) + 2 * arr (flag, i, int) ;
      if (!key1 || !key2)
	continue ;
/* verification qu'il y a au moins un subclone en commun */
      ks1 = queryKey (key1, ">Assembled_from ; >Subclone") ;
      ks2 = queryKey (key2, ">Assembled_from ; >Subclone") ;
      ks3 = keySetAND (ks1, ks2) ;
      imax = keySetMax (ks3) ;
      keySetDestroy (ks1) ;
      keySetDestroy (ks2) ;
      keySetDestroy (ks3) ;
      if (!imax)
	continue ;
/* fin de verif */

      join = dnaAlignAsmbPaire21 (key1, key2, 0, 0, 7, fl, cp) ;
      if (!join)
	join = dnaAlignAsmbPaire21 (key1, key2, 0, 0, 4, fl, cp) ;
      if (!join)
	join = dnaAlignAsmbPaire21 (key1, key2, 0, 0, -1, fl, cp) ;
      if (!join)
	continue ;
      nj++ ;
      arrp (subseq, i - 1, BSunit)->k = 0 ;
      arrp (subseq, i, BSunit)->k = join ;
      arr (flag, i, int) = 1 ;
    }
  if (!nj)
    goto abort ;
  lexaddkey (cp, &link, _VSequence) ;
  tim = timeParse ("today") ;
  if ((obj = bsUpdate (link)))
    if (bsAddKey (obj, _Derived_from, look->link))
      bsAddData (obj, _bsRight, _DateType, &tim) ;
  bsSave (obj) ;
  defCptChangeLook (look, look->link, link) ;
  ks = keySetCreate () ;
  look->def = keySetReCreate (look->def) ;
  for (i = 0 ; i < arrayMax (subseq) ; i++)
    if ((join = arrp (subseq, i, BSunit)->k))
      { 
	if (arr (flag, i, int))
	  dnaAlignFixContig (look->link, join) ;
	else if (!dnaAlignCopyContig (join, &join, cp, 0))
	  continue ;
	keySet (ks, j) = join ;
	dnaSubClass (join, &join) ;
	keySet (look->def, j++) = join ;
      }
  alignToolsAdjustLink (look->link, ks, 0) ;
  keySetDestroy (ks) ;
  keySetSort (look->def) ;
 abort:
  messfree (cp) ;
  bsDestroy (obj) ;
  arrayDestroy (subseq) ;
  return nj ;
}

/***************************************************************/

static int defCptDoJoinSegments(DEFCPT look)
{ KEY key, key1, join ;
  int i, j, max, nj = 0 ;
  Array defOrd = 0 ;
  OBJ obj = 0 ;
  KEYSET def, aa = keySetCreate() , jj = keySetCreate () ;
  char timeBuf[25] ;

  messStatus ("Joining") ;
  def = look->def ;
  if (!def || !arrayMax(def))
    return 0 ;
  if (!look->tour)		/* premier passage */
    alignToolsDestroy_Segs (look->id) ;
  look->tour++ ;
  defOrd = arrayCreate (arrayMax (def), KEY_INT) ;
  j = 0 ; nj = 0 ;
  for (i = 0 ; i < keySetMax(def) ; i++)
    { 
      key1 = keySet(def, i) ;
      dnaReClass(key1, &key1) ;
      keySet(def, i) = key1 ;
      if ((obj = bsCreate (key1)) && bsGetKey (obj, _DNA, &key) &&
	  bsGetData (obj, _bsRight, _Int, &max))
	{ arrayp (defOrd, j, KEY_INT)->key = key1 ;
	  arrp (defOrd, j++, KEY_INT)->tai = max ;
	}
      bsDestroy (obj) ;
    }
  arraySort (defOrd, defCptContigOrder) ;
  for (i = 0 ; i < arrayMax(defOrd) ; i++)
    { key = arrp (defOrd, i, KEY_INT)->key ;
      if (!key)
	continue ;
      if (arrp (defOrd, i, KEY_INT)->tai < 2000) /* end the loop */
	break ;
      for (j = i + 1 ; j < arrayMax(defOrd) ; j++)
      { key1 = arrp (defOrd, j, KEY_INT)->key ;
	if (!key1)
	  continue ;
	join = dnaAlignAsmbPaire21 (key, key1, look->id, look->tour, 7, 0, 0) ;
	if (!join)
	  continue ;
	keySet (jj, nj++) = join ;
	dnaReClass(join, &join) ;
	key = arrp (defOrd, i, KEY_INT)->key = join ;
	arrp (defOrd, j, KEY_INT)->key = 0 ;
	j = i ; /* reiterate the inner j loop */
      }
    }
  /* collect the contigs */
  aa = keySetCreate () ; j = 0 ;
  for (i = 0 ; i < arrayMax(defOrd) ; i++)
    { key = arrp (defOrd, i, KEY_INT)->key ;
      if (key)
	{ 
	  dnaSubClass(key, &key) ;
	  keySet (aa, j++) = key ;
	}
    }
  arrayDestroy (defOrd) ;

  if (look->mapStatus && arrayMax(def))
    { if (nj)
	look->mapStatus = (F_ASSSEG | 5 | (look->mapStatus & F_HIST)) ;
#ifndef NON_GRAPHIC
      defcomputeDraw (look, 
	 messprintf(" // I made %d joins -> %d contigs", nj, keySetMax (look->def))) ;
#endif /* !NON_GRAPHIC */
    }
  if (diffaction)
    catText(diffaction, messprintf("join // -> Made %d joins ; %s\n",
				   nj, timeShow(timeNow(), timeBuf, 25))) ;
  keySetSort(aa) ;
  keySetCompress(aa) ;
  keySetDestroy (def) ;
  look->def = arrayCopy (aa) ;
  if (look->display)
    dnaDispGraph (aa, look->tour) ;
  for (i = 0 ; i < arrayMax(aa) ; i++)
    defCptCollectDummy (keySet (aa, i), jj, 0) ;
  keySetDestroy (jj) ;
  if (!nj) /* no join */
    look->tour-- ; /*en reponse au ++ du haut de la fonction */
  return nj ;
}

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptJoinSegments(void)
{ DEFCPTGET("Joining Segments") ;
  defCptDoJoinSegments (look) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptSegments(void)
{ DEFCPTGET("Compute Segments") ;

  defCptDoSegments(look) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/
/* tri les segments qui restent dans l'ordre decroissant de la taille
   et essaye de les assembler dans cet ordre */
#ifndef NON_GRAPHIC
static void defCptDoTryAll(DEFCPT look)
{ KEYSET defbis ;
  char timeBuf[25] ;

  messStatus("Assembling all") ;
  if (!defCptFixErrorRate(look))
    return ;
  defbis = look->def ;
  look->def = keySetCreate() ;
  dnaAlignTryAll (look) ;
  if (diffaction)
    catText(diffaction, messprintf("tryall %d // -> Found %d contigs, rejected %d Baddies ; %s\n",
				   look->taux, keySetMax(look->def),
				   look->rejected ? 
				   keySetMax (look->rejected) : 0,
				   timeShow(timeNow(), timeBuf, 25))) ;
  defCptEndAssfunc(look, defbis) ;
} /* defCptDoTryAll */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptTryAll(void)
{ DEFCPTGET("Try All") ;

  defCptDoTryAll(look) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptDoNewTryPaires(DEFCPT look)
{ KEYSET defbis ;
  char timeBuf[25] ;

  messStatus("Assembling Paires") ;
  if (!defCptFixErrorRate(look))
    return ;
  defbis = look->def ;
  look->def = keySetCreate() ;
  dnaAlignNewTryPaires(look) ;
  if (diffaction)
    catText(diffaction, messprintf("newtrypaires %d // -> Found %d contigs, rejected %d Baddies ; %s\n",
				   look->taux, keySetMax(look->def),
				   look->rejected ? 
				   keySetMax (look->rejected) : 0,
				   timeShow(timeNow(), timeBuf, 25))) ;
  defCptEndAssfunc(look, defbis) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptNewTryPaires(void)
{ DEFCPTGET ("Try Paires") ;

  defCptDoNewTryPaires(look) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

BOOL defCptReport (KEY link, Array tgdContig, int tglimit, int *nTGcontig, int *nTGreads,
		   int *nTGbase, int limit, int *nGcontig, int *nGreads, int *nGbase, 
		   int *nPcontig, int *nPreads, int *nPbase)
{ OBJ obj = 0 ;
  Array units = 0, unit2 = 0 ;
  BSunit *u ;
  int i, j = 0, max, nbr, ntgc = 0, ntgr = 0, ntgb = 0, ngc = 0, ngr = 0, ngb = 0 ;
  int npc = 0, npr = 0, npb = 0 ;
  KEY dnakey ;

  units = arrayCreate (30, BSunit) ;
  if (!(obj = bsCreate (link)) || !bsFindTag (obj, _Subsequence) ||
      !bsFlatten (obj, 1, units) || !(i = arrayMax (units)))
    { arrayDestroy (units) ;
      bsDestroy (obj) ;
      return FALSE ;
    }
  bsDestroy (obj) ;
  u = arrp (units, 0, BSunit) - 1 ;
  while (u++, i--)
    {
      if (!(obj = bsCreate (u->k)))
	continue ;
      if (bsGetKey (obj, _DNA, &dnakey) && 
	  bsGetData (obj, _bsRight, _Int, &max))
	{ unit2 = arrayReCreate (unit2, 200, BSunit) ;
	  if (bsFindTag (obj, _Assembled_from) && bsFlatten (obj, 1, unit2))
	    nbr = arrayMax (unit2) ;
	  else nbr = 1 ; /* should not happen */
	  if (max > tglimit)
	    { array (tgdContig, j++, int) = max ;
	      array (tgdContig, j++, int) = nbr ;
	      ntgc++ ;
	      ntgr += nbr ;
	      ntgb += max ;
	    }
	  else if (max > limit) 
	    { ngc++ ;
	      ngr += nbr ;
	      ngb += max ;
	    }
	  else
	    { npc++ ;
	      npr += nbr ;
	      npb += max ;
	    }
	}
      bsDestroy (obj) ;
    }
  arrayDestroy (units) ;
  arrayDestroy (unit2) ;
  *nTGcontig = ntgc ;
  *nTGreads = ntgr ;
  *nTGbase = ntgb ;
  *nGcontig = ngc ;
  *nGreads = ngr ;
  *nGbase = ngb ;
  *nPcontig = npc ;
  *nPreads = npr ;
  *nPbase = npb ;
  return TRUE ;
}

/***************************************************************/

static void defCptPreReport (DEFCPT look)
{ int i, ntgc, ntgr, ntgb, ngc, ngr, ngb, npc, npr, npb, nrr ;
  Array tgdContig ;
  KEYSET ks ;

  if (!(look->link))
    { if (look->def)
	fprintf (defCptOut, "%d Contigs", keySetMax (look->def)) ;
      else
	fprintf (defCptOut, "No sequences loaded") ;
      return ;
    }
  tgdContig = arrayCreate (30, int) ;
  if (defCptReport (look->link, tgdContig, 3000, &ntgc, &ntgr, &ntgb, 1000,
		    &ngc, &ngr, &ngb, &npc, &npr, &npb))
    { fprintf (defCptOut, "\n// Resultat :\n") ;
      fprintf (defCptOut, "// Assembly %s\n", name (look->link)) ;
      if (arrayMax (tgdContig))
	{ fprintf (defCptOut, "// Contig > 3 kb : %d\n", arrayMax (tgdContig) / 2) ;
	  fprintf (defCptOut, "//         Total : %d bases, %d reads\n", ntgb, ntgr) ;
	  for (i = 0 ; i < arrayMax (tgdContig) ; i += 2)
	    fprintf (defCptOut, "//      %d bases en %d reads\n", arr(tgdContig, i, int),
		     arr (tgdContig, i + 1, int)) ;
	}
      ks = query (0,"Find Read") ; nrr = keySetMax(ks) ; keySetDestroy (ks) ;
      fprintf (defCptOut, "// Contig between 3 and 1 kb : %d\n", ngc) ;
      fprintf (defCptOut, "//             Total : %d bases, %d reads\n", ngb, ngr) ;
      fprintf (defCptOut, "// Contig <  1 kb : %d\n", npc) ;
      fprintf (defCptOut, "//             Total : %d bases, %d reads\n", npb, npr) ;
      fprintf (defCptOut, "// Total number of contigs            : %8d\n", ntgc+ngc+npc) ;
      fprintf (defCptOut, "// Total number of incorporated reads : %8d/%d present in database\n", ntgr+ngr+npr, nrr) ;
      fprintf (defCptOut, "// Total number of bases              : %8d\n", ntgb+ngb+npb) ;
    }
  arrayDestroy (tgdContig) ;
}

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptSaveAs(void)
{
  char timeBuf[25] ;
  char *cq = 0 ;
  DEFCPTGET("defCptSaveAs") ;

  dnaAlignSaveAs(look, &cq) ;
  look->step = 0 ;
  if (cq)
    {
      if (diffaction)
	catText(diffaction, messprintf("save %s // save consensus as ; %s", cq, timeShow(timeNow(), timeBuf, 25))) ;
      messfree(cq) ;
    } /* faut-il stocker le nom du consensus dans le game ou non ? */
} /* defCptSaveAs */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptAssHisto(void)
{ defCptHisto("Sequence", "Oligo") ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/
/********************   Genes  Alignment   *********************/
/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptFilter(DEFCPT look, int limit)
{ int i, ng ;
  KEYSET *kstmp ;

  ng = arrayMax(look->defInMar) ;
  i = arrayMax(look->marInDef) ;
  kstmp = arrp(look->marInDef, 0, KEYSET) - 1 ;
  while(kstmp++, i--)
    if ((100 * keySetMax(*kstmp)) / ng < limit)
      keySet(*kstmp, 0) |= F_ISOLE ;
} /* defCptFilter */
#endif /* !NON_GRAPHIC */

/***************************************************************/

/* dans ce cas les def sont des alleles de genes et les genes sont les
animaux ou on a teste les alleles comme issus du pere ou de la mere */

#ifndef NON_GRAPHIC
static void defCptDoGetGenes(DEFCPT look, KEYSET aa)
{ KEYSET tmp, tmq, *tmr ;
  KEY key, newkey ;
  int i, j, k, limit ;
  OBJ obj ;
  KEY _LocusP, _LocusM ;

  lexaddkey ("LocusP", &_LocusP, 0) ;
  lexaddkey ("LocusM", &_LocusM, 0) ;

  /* aa is a set of Loci */  
  look->mar = query(aa, "{> GameteP } $| {> GameteM}") ;
  if (!keySetMax(look->mar))
    { keySetDestroy(look->mar) ;
      messout("First (2) select a keyset containing gametes") ;
      return ;
    }
  messStatus("Get Genes") ;
  chrono("Matrice") ;
  tmp = query(look->mar, messprintf("> LocusP")) ;
  tmq = query(look->mar, messprintf("> LocusM")) ;
  look->def = keySetOR(tmp, tmq) ;
  keySetDestroy(tmp) ;
  keySetDestroy(tmq) ;
  if (!keySetMax(look->def))
    { keySetDestroy(look->mar) ;
      keySetDestroy(look->def) ;
      messout("First (3) select a keyset containing gametes") ;
      return ;
    }
  keySetSort(look->mar) ;
  keySetCompress(look->mar) ;
  keySetSort(look->def) ;
  keySetCompress(look->def) ;
  look->nd = keySetMax(look->def) ;
  look->nm = keySetMax(look->mar) ;
  i = look->nd ;
  look->marInDef = arrayCreate(i, KEYSET) ;
  array(look->marInDef, i - 1, KEYSET) = 0 ;
  tmr = arrp(look->marInDef, 0, KEYSET) - 1 ;
  while(tmr++, i--)
    *tmr = keySetCreate() ;
  i = look->nm ;
  look->defInMar = arrayCreate(i, KEYSET) ;
  array(look->defInMar, i - 1, KEYSET) = 0 ;
  tmr = arrp(look->defInMar, 0, KEYSET) - 1 ;
  while(tmr++, i--)
    *tmr = keySetCreate() ;
  for (j = 0 ; j < look->nm ; j++)
    { key = keySet(look->mar, j) ;
      tmq = arr(look->defInMar, j, KEYSET) ;
      if ((obj = bsCreate(key)))
	{ k = F_PERE ;
	  if (bsGetKey(obj, _LocusP, &newkey))
	    do { keySetFind(look->def, newkey, &i) ;
		 tmp = arr(look->marInDef, i, KEYSET) ;
		 keySet(tmp, keySetMax(tmp)) = (j | k) ;
		 keySet(tmq, keySetMax(tmq)) = (i | k) ;
	       } while(bsGetKey(obj, _bsDown, &newkey)) ;
	  k = F_MERE ;
	  if (bsGetKey(obj, _LocusM, &newkey))
	    do { keySetFind(look->def, newkey, &i) ;
		 tmp = arr(look->marInDef, i, KEYSET) ;
		 keySet(tmp, keySetMax(tmp)) = (j | k) ;
		 keySet(tmq, keySetMax(tmq)) = (i | k) ;
	       } while(bsGetKey(obj, _bsDown, &newkey)) ;
	  bsDestroy(obj) ;
	}
    }
  i = look->nd ;
  tmr = arrp(look->marInDef, 0, KEYSET) - 1 ;
  while(tmr++, i--)
    arraySort(*tmr, mardefOrder) ;
  i = look->nm ;
  tmr = arrp(look->defInMar, 0, KEYSET) - 1 ;
  while(tmr++, i--)
    arraySort(*tmr, mardefOrder) ;
  arrayDestroy (look->colOrder) ;
  arrayDestroy (look->linOrder) ;
  TESCOLO ;
  TESLIN ;
  chronoReturn() ;
 ici:
  if (messPrompt("unknown limit", "", "i"))
    { freeint(&limit) ;
      if (limit < 0 || limit > 100)
	{ messout("please choose a value between 0 and 100") ;
	  goto ici ;
	}
      defCptFilter(look, 100 - limit) ;
    }
  chrono("tree") ;
  look->dlimit = 100 ;
  intCptTree(look) ;
  chronoReturn() ;
  if (look->defMapCtlGraph)
    graphActivate(look->defMapCtlGraph) ;
  defCptSortMarkers() ;
  defMapDisplay() ;
  look->mapStatus = 1 ;
  defcomputeDraw (look, 0) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptGetGenes(void)
{ KEYSET aa ;
  void *dummy ;
  DEFCPTGET("Get Genes") ;

  if (!keySetActive(&aa, &dummy))
    { messout("First (1) select a keyset containing gametes") ;
      return ;
    }
  defCptDoGetGenes(look, aa) ;
}
#endif

/***************************************************************/

KEY lastAssembly (void) 
{ KEY last = 0 ;
  int i ;
  KEYSET ks = query (0, ">?Assembly NOT IS _*") ;

  i = keySetMax (ks) ;
  if (i) last = keySet (ks, i - 1) ;
  keySetDestroy (ks) ;
  return last ;
}

/***************************************************************/

static void defCptJuxtapose (DEFCPT look)
{ int n ;
  KEYSET ks0, ks1, ks2, ks3, ks4 ;
  KEY key ;

  ks0 = keySetCreate () ;
  for (n = 0 ; n < keySetMax(look->def) ; n++)
    { 
      dnaReClass (keySet(look->def, n), &key) ;
      keySet (ks0, n) = key ; 
    }
  ks1 = query (ks0, "{NOT Subsequence } $| {> Subsequence}") ;
  keySetDestroy (look->def) ;

  ks2 = query (ks1, "!Assembled_from AND DNA") ;
  ks3 = query (ks1, "Assembled_from") ;
  ks4 = keySetOR (ks2, ks3) ;
  keySetSort (ks4) ;
  keySetCompress (ks4) ;
  contigLengthSort (ks4) ;
  look->def = ks4 ; ks4 = 0 ;
  if (!look->link)
    lexaddkey ("Phrap", &look->link, _VSequence) ;
  alignToolsAdjustLink (look->link, look->def, 0) ;
  keySetDestroy (ks0) ;
  keySetDestroy (ks1) ;
  keySetDestroy (ks2) ;
  keySetDestroy (ks3) ;
}

/***************************************************************/
/****************************************************************/

BOOL defCptChromosomeOrigin (int *zero, int *trueZero, KEY *maps)
{
  int ii, a1, a2 ;
  char *cp ;
  KEY map, gene = 0 ;
  KEYSET ks = 0 ;
  OBJ obj = 0 ;

  for (ii = 0 ; ii < 7 ; ii++)
    { 
      zero[ii] = 0 ; trueZero[ii] = 0 ; maps[ii] = 0 ;
      map = 0 ; cp = 0 ;
      switch (ii)
	{
	case 0:
	  continue ;
	case 1:
	  /*  G_YK1557 1G0" dpy-5, f27c1.8, mail ann rose 14 avril 99   1G0:6:5.3Mb(jan2001) */
	  cp = "(IS dpy-5 || IS  f27c1.8 || IS 1G0)" ;
	  break ;
	case 2:
	  /* G_YK3231 2H0  dpy-10 t14b4.7  not Zero on rd's map !    2H0:7:6.7Mb(jan2001) */
	  cp = "(IS dpy-10 || IS t14b4.7 || IS t14b4.7a || IS 2H0)" ;
	  break ;
	case 3:
	  /* G_YK91 3K1  unc-32,         3K1:9:8.9Mb (jan2001) */
	  cp = "(IS unc-32 || IS  zk637.8a || IS 3K1)" ;
	  break ;
	case 4:
	  /* G_YK1902 4F1 dpy-13 f30b5.1                            4F1:5:4.2Mb(jan2001) */
	  cp = "(IS dpy-13 || IS  f30b5.1 || IS 4F1)" ;
	  break ;
	case 5:
	  /* G_YK3192 5H1 , no knownsequence, danielle's bet   5H1:7:6.7Mb(jan2001) */
	  /* faux de 8k, vrai dpy-11 un peu upstream dpy-11 ==  f46e10.9 == 5G992 */
	  cp = "(IS  G_YK3192 || IS  f46e10.2 || IS 5H1)" ;
	  break ;
	case 6:  /* X */
	  /* G_YK1353 XJ1 dpy-6 f41d9.1 g_yk1353 (why? not         XJ1:9:8.4Mb(jan2001) ???*/
	  cp = "(IS dpy-6 || IS f41d9.1 || IS XJ1)" ;
	  break ;
	}
      ks = query (0, messprintf ("{Find newname %s ; >gene ; IntMap} $/ {locus %s; >gene ; IntMap} $/ {Find locus %s; >Transcribed_gene ; IntMap}  $/ {Find sequence %s ; IntMap} $/ {Find locus %s; >sequence ; IntMap}",cp,cp,cp,cp,cp)) ;
      if (keySetMax(ks))
	{
	  gene = keySet (ks,0) ;
	  if ((obj = bsCreate (gene)) )
	    {
	      if (bsGetKey (obj, str2tag("IntMap"), &map) &&
		  bsGetData (obj, _bsRight, _Int, &a1) &&
		  bsGetData (obj, _bsRight, _Int, &a2)) ;
	      else map = 0 ;
	      bsDestroy (obj) ;
	    }

	}
      keySetDestroy (ks) ;
      if (!map)
	continue ;
      trueZero[ii] = a1 > a2 ? -a1 : a2 ;

      if (trueZero[ii] <= 0) 
	{  
	  trueZero[ii] *= -1 ;
	  zero[ii] = trueZero[ii] ; zero[ii]  %= 1000000 ;  
	  if (zero[ii] > 0) zero[ii] -= 1000000 ;
	}
      else 
	{ 
	  zero[ii] = trueZero[ii] - 1000 ; 
	  zero[ii]  %= 1000000 ; 
	  if (zero[ii] > 0) zero[ii] -= 1000000 ; 
	}
      if (ii == 3)
	zero[ii] -= 1000000 ;
      maps[ii] = map ;
      printf ("%d: Chromosome %s gene %s   zero = %d truezero = %9d\n", ii, name(map), name(gene), zero[ii], trueZero[ii]) ;

    }
  return TRUE ;
}

/***************************************************************/
/* droso chromo are labelled X Xh, 2L 2h 2R, 3L 3h 3R, 4, U, Yh */
static int drosoOrder (const void *va, const void *vb)
{
  KEY a = *(const KEY*)va, b = *(const KEY*)vb ;
  char *na = name(a), *nb = name(b) ;

  if (*na != *nb)
    return (na - nb) ;
  if (!*(na+1) && *(nb+1)) return -1 ;
  if (*(na+1) && !*(nb+1)) return 1 ;
  if (*(na+1) == 'L') return -1 ;
  if (*(nb+1) == 'L') return 1 ;
  if (*(na+1) == 'R') return 1 ;
  if (*(nb+1) == 'R') return -1 ;
  return 0 ;
}

/***************************************************************/
/* split into pieces */
static void defCptMakeSubsequences (DEFCPT look)
{ 
  int size = 0, overlap = -1, i, ii, jj, jjMax, nn, myZero, myTrueZero ; 
  int zero[7], trueZero[7] ;
  char *buffer = 0 , *cp = 0, *prefix = 0, *prefix2 = 0, *source = 0, *subsource = 0;
  KEY seq, maps[7] ; 
  Array dna = 0 ;
  BOOL withDna = FALSE, isWorm = FALSE, isDroso = FALSE, isAra = FALSE ;
  Stack s = stackCreate (10000000) ;
  char nnbuf[512], nnbuf2[512] ;
  
  ii = keySetMax (look->taceActif) ;
  cp = freeword () ; 
  if (cp && !strcmp(cp,"-dna"))
    { withDna = TRUE ; cp = freeword() ; }
  if (cp) prefix = strnew (cp, 0) ;
  cp = freeword () ; 
  if (cp) prefix2 = strnew (cp, 0) ;
  freeint(&size) ;
  freeint(&overlap) ;
  cp = freeword () ; 
  if (cp  && !strcmp(cp,"-worm")) isWorm = TRUE ;
  if (cp  && !strcmp(cp,"-ara")) isAra = TRUE ;
  if (cp  && !strcmp(cp,"-droso")) isDroso = TRUE ;
  
  if (!ii || size < 1000 || !prefix || !prefix2 || overlap < 0 || overlap > size/3)
    {
      messout ("%s\n%s\n%s",
	       "Please provide 2 prefixes (c t), one for the future source sequence",
	       "another for the future subsequences",
	       "then a size (50000: > 1000) and overlap (200: > size/3)") ;
      goto done ;
    } 
  
  if (isWorm || isDroso || isAra)
    {
      size = 110000 ; overlap = 10000 ;
      defCptChromosomeOrigin (zero, trueZero, maps) ;
    }
  

  printf ("// %d sequences size = %d overlap = %d prefix = %s prefix2 = %s\n",
	  ii, size, overlap, prefix, prefix2) ;
  buffer = messalloc (size + 1) ;

  if (isDroso)
    arraySort (look->taceActif, drosoOrder) ;
  for (ii = 0 ; ii < keySetMax(look->taceActif) ; ii++)
    {
      seq = keySet(look->taceActif, ii) ;
      dna = dnaGet(seq) ;
      if (!dna)
	{
	  messerror ("No dna in defCptMakeSubsequences %s", name(seq)) ;
	  continue ;
	}
      dnaDecodeArray (dna) ;
      jj =  dna ? arrayMax(dna) : 0 ;
      printf ("Sequence %s length %d\n", name(seq), jj) ;

      myZero = myTrueZero = 0 ;
      for (i = 0 ; i < 7 ; i++)
	if (!strcmp (name(maps[i]), name (seq)))
	  {
	    myZero = zero[i] ;
	    myTrueZero = trueZero[i] ;
	    break ;
	  }

      if (isWorm)
	{
	  if (!strcmp (name(seq),"CHROMOSOME_I")) cp = "1" ;
	  else if (!strcmp (name(seq),"CHROMOSOME_II")) cp = "2" ;
	  else if (!strcmp (name(seq),"CHROMOSOME_III")) cp = "3" ;
	  else if (!strcmp (name(seq),"CHROMOSOME_IV")) cp = "4" ;
	  else if (!strcmp (name(seq),"CHROMOSOME_V")) cp = "5" ;
	  else if (!strcmp (name(seq),"CHROMOSOME_X")) cp = "X" ;
	  else continue ;

	  source = strnew (messprintf("%s", name(seq)), 0) ;
	  subsource = strnew (messprintf("%s",cp), 0) ;
	}
      else if (isAra)
	{
	  if (!strcmp (name(seq),"CHR1")) cp = "1" ;
	  else if (!strcmp (name(seq),"CHR2")) cp = "2" ;
	  else if (!strcmp (name(seq),"CHR3")) cp = "3" ;
	  else if (!strcmp (name(seq),"CHR4")) cp = "4" ;
	  else if (!strcmp (name(seq),"CHR5")) cp = "5" ;
	  else continue ;

	  source = strnew (messprintf("%s", name(seq)), 0) ;
	  subsource = strnew (messprintf("%s",cp), 0) ;
	}  
      else if (isDroso)
	{
	  cp = name(seq) ;
	  source = strnew (messprintf("CH%c", *cp), 0) ;
	  subsource = strnew (messprintf("%c",*cp), 0) ;
	}
      else
	{
	  source = strnew (messprintf("%s_%s",prefix,name(seq)), 0) ;
	  subsource = strnew (messprintf("%s_%s",prefix2,name(seq)), 0) ;
	}

      jjMax =  arrayMax(dna) ;
      for (jj = myZero, nn = 0 ; jj < jjMax ; nn++, jj += size)
	{
	  int j1, j2, dj1 ;

	  j1 = jj + 1 ; j2 = jj + size ;
	  if (j2 > jjMax)
	    j2 = jjMax ;
	  if (j1 < 1) { dj1 = 1 - j1 ; j1 = 1 ; }
	  else dj1 = 0 ;
	  if (j2 < 1) j2 = 1 ;
	  if (jj == myZero || j2 - j1 > overlap)
	    {
	      if (isWorm)
		{
		  sprintf(nnbuf, "%s%c%d", subsource, 'A' + nn/10, (nn%10)*10) ;
		  sprintf(nnbuf2, "%s%c%d", subsource, 'A' + (nn+1)/10, ((nn+1)%10)*10) ;
		}
	      else
		{
		  sprintf(nnbuf, "%s_%d", subsource, nn) ;
		  sprintf(nnbuf2, "%s_%d", subsource, nn + 1) ;
		}
	      catText (s, messprintf ("Sequence %s\nSubsequence %s %d %d\n\n",
				      source, nnbuf, j1, j2)) ;
	      catText (s, messprintf ("Sequence %s\ngenomic\n\n",
				      nnbuf)) ;
	      if (j2 < jjMax)
		catText (s, messprintf ("Sequence %s\nOverlap_right %s %d\n\n",
					nnbuf, nnbuf2, size + 1 - overlap - dj1)) ;
	      if (withDna)
		{
		  catText (s, messprintf ("DNA %s\n", nnbuf)) ;
		  memset(buffer, 0, size+1) ;
		  cp = arrp(dna, j1 - 1, char) ;
		  memcpy(buffer, cp, j2 - j1 + 1) ;
		  catText (s, buffer) ;
		  catText (s,"\n\n") ;
		}
	    }
	  jj -= overlap ;
	}
      arrayDestroy (dna) ;
      if ( !(nn %20))
	{  parseBuffer (stackText(s,0), 0) ; stackClear (s) ; }
    }
  parseBuffer (stackText(s,0), 0) ;

 done:
  messfree (prefix) ;
  messfree (prefix2) ;
  messfree (source) ;
  messfree (subsource) ;
  messfree (buffer) ;
  stackDestroy (s) ;
}

/***************************************************************/
/***************************************************************/
/*  BOOL SPAN  
 * TRUE: export a set of nn-mer without repetitions staring at every position,  spanning all the mrnas, labelled as X or I
 *       this file shoul be used to compute the mapability of the transcriptome by realigning it on the fasta file of all mrnas
 *
 * FALSE: export the exons, then the introns, then SLs (if worm) and then polyA 
 *        this file should be used to find solexa tags supporting the exons,introns, polyA and SL elements of the transcripts
 */
typedef struct coverFlagStruct {KEY map, mrna ; KEYSET mrnas ; int dna, dm, a1, a2, pos, sl ; char type ;} CF ;
static void defCptMakeDnOneCoverDict (KEY mrna, int nn,  DICT *dict, Array flags, BOOL span, BOOL isWorm, BOOL justIntron, BOOL isPg, AC_HANDLE h)
{
  CF *cf ;
  KEY Nmrna ;
  OBJ Mrna = bsCreate (mrna) ;
  Array dna = 0 ;
  Array aa = arrayCreate (128, BSunit) ;
  BSunit *uu ;
  int a1, a2, m1, m2, dm, amax, g1, g2, ii, i, j, n, iSL, u1, u2 ;
  char *cp, *cq, *cr, cc, buf18[nn+19], bufSL[nn+23], bufI[2*nn+20] ;
  const char *ccp ;
  char *SL[] = {"toto"
		, "GGTTTAATTACCCAAGTTTGAG"
		, "GGTTTTAACCCAGTTACTCAAG"
		, "GGTTTTAACCCAGTTAACCAAG"
		, "GTTTTAACCCAGTTTAACCAAG" /* SL4: left G clipped */
		, "NGGTTTTAACCCAGTTACCAAG"
		, "GGTTTAAAACCCAGTTACCAAG"
		, "GGTTTTAACCCAGTTAATTGAG"
		, "GGTTTTTACCCAGTTAACCAAG"
		, "GGTTTATACCCAGTTAACCAAG"
		, "GTTTTAACCCAAGTTAACCAAG"
		, "NGGTTTTAACCAGTTAACTAAG" /* SL11: left G clipped */
		, "GGTTTTAACCCATATAACCAAG"
		, 0

  } ;

  memset (buf18, 'A', nn+18) ;
  buf18[nn+18] = 0 ;
  g1 = g2 = 0 ;
  if (bsGetArray (Mrna, str2tag("IntMap"), aa, 3))
    {
      uu = arrp (aa, 0, BSunit) ;
      /*  map = uu[0].k ; */
      g1 = uu[1].i ;
      g2 = uu[2].i ;
    }
  dna = dnaGet (mrna) ;
  if (!dna) goto done ;
  dnaDecodeArray (dna) ;
  amax = arrayMax (dna) ;
  if (1 && !isPg && bsGetArray (Mrna, str2tag("Splicing"), aa, 6))
    for (ii = 0 ; ii < arrayMax (aa) ; ii += 6) 
      {
	uu = arrp (aa, ii, BSunit) ;
	a1 = uu[0].i ; a2 = uu[1].i ;
	m1 = uu[2].i ; m2 = uu[3].i ;
	if (strstr(name(uu[4].k), "tron")) /* export the introns */
	  {
	    ccp = uu[5].s ;
	    for (i = 0 ; i < (span ? nn-1 : 1) ; i++)
	      {
		cr = bufI ;
		memset (bufI, 0, sizeof(bufI)) ;
		u1 = m1 - nn + i ; 
		if (u1 < 0 && span && i < nn-2) continue ;
		if (u1 < 0) 
		  {
		    u1 = -u1 ;
		    memset(cr, 'N', u1) ;
		    cr += u1 ;
		    u1 = 0 ;
		  }
		dm = m1 - u1 ;
		cp = arrp (dna, u1, char) ;
		if (span)
		  { u2 = u1 + nn -1 ;  if (u2 > amax && i > 0) continue ;}
		else
		  { u2 = m1 + nn - 1 ; }
		if (u2 > amax) u2 = amax ;
		if (u2 - u1 < nn - 1) continue ;
		cq = arrp (dna, u2 , char) ;
		if (g1 < g2) { u1 = g1 + a1 - 1 ; u2 = g1 + a2 - 1 ; }
		else { u1 = g1 - a1 + 1 ; u2 = g1 - a2 + 1 ; }
		memcpy (cr, cp, cq - cp + 1) ;
		dictAdd (dict, bufI, &n) ; 
		cf = arrayp (flags, n, CF) ;
		if (cf->mrna == mrna) ;
		else if (!cf->mrna) 
		  { cf->mrna = mrna ; cf->pos = m1 ; }
		else
		  {
		    if (! cf->mrnas) 
		      {
			cf->mrnas = arrayHandleCreate (32, KEY, h) ;
			keySet (cf->mrnas, 0) = cf->mrna ;
		      }
		    keySetInsert (cf->mrnas, mrna) ;
		  }
		cf->a1 = u1 ; cf->a2 = u2 ; cf->type |= 0x4 ; cf->dm = nn ; cf->dna = n ;
		if (!strcmp (ccp, "gt_ag"))
		  cf->type |= 0x20 ;
		if (!strcmp (ccp, "gc_ag"))
		  cf->type |= 0x10 ;
	      }
	  }
	else if (!justIntron && strstr(name(uu[4].k), "xon")) /* export the exons */
	  {
	    for (i = m1 - 1 ; i <= (span ? m2 - nn : m1 - 1) && i <= amax - nn  ; i++)
	      {
		cp = arrp (dna, i, char) ;
		if (strlen(cp) < nn) continue ;
		if (span) 
		  cq = cp + nn ;
		else
		  cq = arrp (dna, m2, char) ;
		cc = *cq ; *cq = 0 ;
		dictAdd (dict, cp, &n) ; 
		cf = arrayp (flags, n, CF) ;
		if (cf->mrna == mrna) ;
		else if (!cf->mrna) 
		  { cf->mrna = mrna ; cf->pos =  i - m1 + a1 ; }
		else
		  {
		    if (! cf->mrnas) 
		      {
			cf->mrnas = arrayHandleCreate (32, KEY, h) ;
			keySet (cf->mrnas, 0) = cf->mrna ;
		      }
		    keySetInsert (cf->mrnas, mrna) ;
		  }
		if (g1 < g2) { u1 = g1 + a1 - 1 + (i - m1 + 1) ; u2 = (span ? u1 + nn - 1 : g1 + a2 - 1) ; }
		else { u1 = g1 - a1 + 1 - (i - m1 + 1) ; u2 = (span ?  u1 - nn + 1 : g1 - a2 + 1) ;}
		cf->a1 = u1 ; cf->a2 = u2 ; cf->type |= 0x8 ; cf->dna = n ;
		*cq = cc ;
	      }	
	  }
      }


  if (isPg && bsGetArray (Mrna, str2tag("Source_Exons"), aa, 3))
    for (ii = 0 ; ii < arrayMax (aa) ; ii += 3) 
      {
	/* cumulate the previous exons to find the position in pg coordinates */
	for (m1 = 1,  i = 0 ; i < ii ; i += 3)
	  {
	    uu = arrp (aa, i, BSunit) ;
	    m1 += uu[1].i - uu[0].i + 1 ;  
	  }
	uu = arrp (aa, ii, BSunit) ;
	a1 = uu[-2].i + 1 ; a2 = uu[0].i - 1 ;  
        m2 = m1 + uu[1].i - uu[0].i ;
	Nmrna = keyGetKey (mrna, str2tag("NM_id")) ;
	if (! Nmrna) Nmrna = mrna ;
	if (ii > 0) /* export the introns */
	  {
	    for (i = 0 ; i < (span ? nn-1 : 1) ; i++)
	      {
		u1 = m1 - nn + 1 + i ; 
		if (u1 < 0 && span && i < nn-2) continue ;
		if (u1 < 0) u1 = 0 ;
		dm = m1 - u1 ;
		cp = arrp (dna, u1, char) ;
		if (span)
		  { u2 = u1 + nn -1 ;  if (u2 > amax && i > 0) continue ;}
		else
		  { u2 = m1 + nn - 1 ; }
		if (u2 > amax) u2 = amax ;
		if (u2 - u1 < nn - 1) continue ;
		cq = arrp (dna, u2 , char) ;
		if (g1 < g2) { u1 = g1 + a1 - 1 ; u2 = g1 + a2 - 1 ; }
		else { u1 = g1 - a1 + 1 ; u2 = g1 - a2 + 1 ; }
		cc = *cq ; *cq = 0 ;
		dictAdd (dict, cp, &n) ; 
		cf = arrayp (flags, n, CF) ;
		if (cf->mrna == Nmrna) ;
		else if (!cf->mrna) 
		  { cf->mrna = Nmrna ; cf->pos = m1 ; }
		else
		  {
		    if (! cf->mrnas) 
		      {
			cf->mrnas = arrayHandleCreate (32, KEY, h) ;
			keySet (cf->mrnas, 0) = cf->mrna ;
		      }
		    keySetInsert (cf->mrnas, Nmrna) ;
		  }
		cf->a1 = u1 ; cf->a2 = u2 ; cf->type |= 0x4 ; cf->dm = dm ; cf->dna = n ;
		*cq = cc ;
	      }
	  }
	
	if (!justIntron) /* export the exons */
	  {
	    for (i = m1 - 1 ; i <= (span ? m2 - nn : m1 - 1) && i <= amax - nn  ; i++)
	      {
		cp = arrp (dna, i, char) ;
		if (strlen(cp) < nn) continue ;
		if (span) 
		  cq = cp + nn ;
		else
		  cq = arrp (dna, m2, char) ;
		cc = *cq ; *cq = 0 ;
		dictAdd (dict, cp, &n) ; 
		cf = arrayp (flags, n, CF) ;
		if (cf->mrna == Nmrna) ;
		else if (!cf->mrna) 
		  { cf->mrna = Nmrna ; cf->pos =  i - m1 + a1 ; }
		else
		  {
		    if (! cf->mrnas) 
		      {
			cf->mrnas = arrayHandleCreate (32, KEY, h) ;
			keySet (cf->mrnas, 0) = cf->mrna ;
		      }
		    keySetInsert (cf->mrnas, Nmrna) ;
		  }
		if (g1 < g2) { u1 = g1 + a1 - 1 + (i - m1 + 1) ; u2 = (span ? u1 + nn - 1 : g1 + a2 - 1) ; }
		else { u1 = g1 - a1 + 1 - (i - m1 + 1) ; u2 = (span ?  u1 - nn + 1 : g1 - a2 + 1) ;}
		cf->a1 = u1 ; cf->a2 = u2 ; cf->type |= 0x8 ; cf->dna = n ;
		*cq = cc ;
	      }	
	  }
      }
  /* export the polyA */
  if (!justIntron && !span && bsGetArray (Mrna, str2tag("Valid3p"), aa, 2))
    for (ii = 0 ; ii < arrayMax (aa) ; ii += 2) 
      {
	uu = arrp (aa, ii, BSunit) ;
	a1 = uu[0].i ;
	m1 = uu[1].i ;
	for (i = 1 ; i < (span ? 20 : 2) ; i++)
	    {
	      if (i < 0) continue ;
	      cp = arrp (dna, m1 -  nn +  i, char) ;
	      strncpy (buf18, cp, nn) ;
	      cp = buf18 + nn - i ;
	      memset (cp, 'A', 18) ;
	      cp = buf18 ;
	      if (span) 
		buf18[nn] = 0 ;
	      else
		buf18[nn+18] = 0 ;
	      dictAdd (dict, cp, &n) ; 
	      cf = arrayp (flags, n, CF) ;
	      if (cf->mrna == mrna) ;
	      else if (!cf->mrna) 
		{ cf->mrna = mrna ; cf->pos = m1 ; }
	      else
		{
		  if (! cf->mrnas) 
		    {
		      cf->mrnas = arrayHandleCreate (32, KEY, h) ;
		      keySet (cf->mrnas, 0) = cf->mrna ;
		    }
		  keySetInsert (cf->mrnas, mrna) ;
		}
	      if (g1 < g2) { u1 = g1 + a1 - 1 ; u2 = u1 + 1 ; }
	      else { u1 = g1 - a1 + 1 ; u2 = u1 - 1 ; }
	      cf->a1 = u1 ; cf->a2 = u2 ; cf->type |= 0x2 ; cf->dm = nn - i ; cf->dna = n ;
	    }	
      }
  /* export the SL1 */
  if (!justIntron && !span && isWorm)
    for (iSL = 1 ; SL[iSL] ; iSL++)
      {	
	strncpy (bufSL, SL[iSL], 22) ;
	for (i = 0 ; i < 100 ; i++)
	  {
	    if (i < amax - nn)
	      {
		cp = arrp (dna, i, char) ;
		if (strlen(cp) < nn) continue ;
		if (i>1 && strncmp(cp-2,"ag", 2)) continue ;
		strncpy (bufSL+22, cp, nn) ;
		bufSL[22+nn-8] = 0 ;
		for (j = 0 ; j < (span ? 14 : 1) ; j++) /* at least 8 letters from the SL */
		  {
		    cp = bufSL + j ;
		    cq = cp + (span ? nn : 22+nn) ;
		    cc = *cq ; *cq = 0 ;
		    dictAdd (dict, cp, &n) ; 
		    cf = arrayp (flags, n, CF) ;
		    if (cf->mrna == mrna) ;
		    else if (!cf->mrna) 
		      { cf->mrna = mrna ; cf->pos = 23 - i ; }
		    else
		      {
			if (! cf->mrnas) 
			  {
			    cf->mrnas = arrayHandleCreate (32, KEY, h) ;
			    keySet (cf->mrnas, 0) = cf->mrna ;
			  }
			keySetInsert (cf->mrnas, mrna) ;
		      }
		    if (g1 < g2) { u1 = g1 + i - 1 ; u2 = u1 + 1 ; }
		    else { u1 = g1 - i + 1 ; u2 = u1 - 1 ; }
		    cf->a1 = u1 ; cf->a2 = u2 ; cf->type |= 0x1 ; cf->sl = iSL ; cf->dm = j ; cf->dna = n ;
		    *cq = cc ;
		  }
	      }
	  }
      }	
 done:
  arrayDestroy (aa) ;
  arrayDestroy (dna) ;
  bsDestroy (Mrna) ;
  return ;
} /* defCptMakeDnOneCoverDict */

/***************************************************************/

static void defCptMakeDnaOneCover (KEY tg, int nn, BOOL span, BOOL isWorm, BOOL justIntron, BOOL isPg)
{
  int i, ii ;
  char *type, buf[1000], *cq ;
  const char *ccg, *ccm ;
  KEY gene, map, mrna ;
  KEYSET geneIds, geneIds2, geneIds3, ks = queryKey (tg, isPg ? "CLASS Predicted_gene" : "CLASS Transcribed_gene ; >mrna") ;
  DICT *dict = 0 ;
  Array flags = 0 ;
  CF *cf ;
  AC_HANDLE h = ac_new_handle () ;

  dict = dictHandleCreate (1000, h) ;
  flags = arrayHandleCreate (1000, CF, h) ;
  gene = keyGetKey (tg, str2tag(isPg ? "Model_of_gene" : "Gene")) ;
  map = keyGetKey (tg, str2tag("IntMap")) ;
  if (isPg)
    geneIds = queryKey (tg, ">GeneId_pg") ;
  else
    geneIds = gene ? queryKey (gene, ">GeneId") : keySetCreate () ;

  for (ii = 0 ; ii < keySetMax (ks) ; ii++)
    {
      mrna = keySet(ks, ii) ;
      defCptMakeDnOneCoverDict (mrna, nn, dict, flags, span, isWorm, justIntron, isPg, h) ;
    }

  for (ii = 0 ; ii < arrayMax (flags) ; ii++)
    {
      cf = arrp (flags, ii, CF) ;
      if (cf->type & 0x8) type = "X" ;
      else if (cf->type & 0x4)
	{
	  if (cf->type & 0x20)
	    type = "Igt_ag" ;
	  else if (cf->type & 0x10)
	    type = "Igc_ag" ;
	  else
	    type = "I" ;
	}
      else if (cf->type & 0x2) type = "A" ;
      else if (cf->type & 0x1) type = messprintf ("SL%d", cf->sl) ;
      else continue ;
      ccg = gene ? name (gene) : name (tg) ;
      ccm = name (cf->mrna) ;
      if (isPg && *ccg == '_') 
	{ 
	  strncpy (buf, ccg + 1, 999) ;
	  ccg = buf ;
	  cq = buf + strlen(ccg) - 2 ;
	  if (*cq == '_' && *(cq+1)>= '1' && *(cq+1)<= '9')
	    *cq = 0 ;
	}
      if (isPg && *ccm == '_') ccm++ ;
      freeOutf ("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%d\t%s"
		, ccg
		, type
		, name(map), cf->a1, cf->a2, cf->dm
		, dictName(dict, cf->dna)
		, cf->pos
		, ccm
		) ;
      
      for (i = 0 ; cf->mrnas && i < keySetMax (cf->mrnas) ; i++)
	{
	  mrna= keySet (cf->mrnas, i) ;
	  if (mrna != cf->mrna)
	    {
	      ccm = name (mrna) ;
	      if (isPg && *ccm == '_') ccm++ ;
	      freeOutf (",%s", name (mrna)) ;
	    }
	}
      freeOutf ("\t") ;
      geneIds3 = geneIds ;
      if (!isPg)
	{
	  if (cf->mrnas)
	    geneIds2 = query (cf->mrnas, "COUNT geneid == 1 ; >GeneId") ;
	  else
	    geneIds2 = queryKey (cf->mrna, "COUNT geneid == 1 ; >GeneId") ;
	  if (keySetMax (geneIds2))
	    geneIds3 = geneIds2 ;
	}
      for (i = 0 ; i < keySetMax (geneIds3) ; i++)
	freeOutf ("%s%s", i>0 ? "," : "", name (keySet(geneIds3, i))) ;
      keySetDestroy (geneIds2) ;
      freeOutf ("\n") ;
    }
  
  keySetDestroy (geneIds) ;
  keySetDestroy (ks) ;
  ac_free (h) ;
  return ;
} /* defCptMakeDnaOneCover */

/***************************************************************/
/* export all tags covering each gene, once per gene */
static void defCptMakeDnaCover (KEYSET ks0)
{
  char *cp = 0 ;
  int nn = 0 ;
  KEYSET ks = 0 ;
  int ii, mx = 0, level = 0 ;
  FILE *f = 0 ;
  BOOL span = FALSE, isWorm = TRUE, justIntron = FALSE, isPg = FALSE ;

  while ((cp=freeword()))
    {
       if (!f && !strcmp (cp, "-o"))
	{
	  f = 0 ;
	  if (!(cp = freeword()))
	    messout ("Missing out_file_name: Usage make_dna_cover -n 35 [-o out_file_name] [-span]") ;
	  else if (! (f = filopen (cp, 0, "w")))
	    messout ("Cannot open the out_file %s", cp) ;
	  if (!f)	    
	    goto done ;
	}
       else if (!nn && !strcmp (cp, "-n"))
	 {
	   if (! freeint (&nn) || (nn < 8))
	     {
	       messout ("Missing -n value or value < 8: Usage make_dna_cover -n 35 [-o out_file_name] [-span] [-justIntron]") ;
	       goto done ;
	     }
	 }
       else if (!span && !strcmp (cp, "-span"))
	 span = TRUE ;
       else if (!span && !strcmp (cp, "-worm"))
	 isWorm = TRUE ;
       else if (!span && !strcmp (cp, "-pg"))
	 isPg = TRUE ;
       else if (!span && !strcmp (cp, "-justIntron"))
	 justIntron = TRUE ;
       messout("%s",cp) ;
    }
  if (!nn) 
    {
      messout ("Missing tag length -n : Usage make_dna_cover -n 35 [-o out_file_name] [-span]") ;
      goto done ;
    }
  ks = query (ks0, "{CLASS Gene; >Transcribed_gene} SETOR {CLASS Transcribed_gene} SETOR {CLASS predicted_gene}") ;
  mx = keySetMax (ks) ;
  if (!mx)
    {
      messout ("No gene or transcribed_gene or mRNA in the active keyset") ;
      goto done ;
    }

  messout ("Start defCptMakeDnaOneCover nn=%d span=%s", nn, span?"TRUE":"FALSE") ;

  if (f) level = freeOutSetFile (f) ;
  for (ii = 0 ; ii < mx ; ii++)
    defCptMakeDnaOneCover (keySet (ks, ii), nn, span, isWorm, justIntron, isPg) ;

 done:
  if (level) freeOutClose (level) ;
  if (f) filclose (f) ;
  keySetDestroy (ks) ;
  return ;
} /* defCptMakeDnaCover */

/***************************************************************/
#ifndef NON_GRAPHIC
static void defCptMakeStatErreur (void)
{ KEYSET aa ;
  void *dummy ;

  if (!keySetActive(&aa, &dummy))
    { messout("First select a keyset containing sequences") ;
      return ;
    }
  statisticsMakeErreur (aa) ;
}

/***************************************************************/

static void defCptMakeStatBarreaux (void)
{ KEYSET aa ;
  void *dummy ;

  if (!keySetActive(&aa, &dummy))
    { messout("First select a keyset containing sequences") ;
      return ;
    }
  statisticsCountGroup (aa) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/
/********************    Menu   options    *********************/
/***************************************************************/

extern long int assBounce , assFound , assNotFound , assInserted , assRemoved  ;

static void defCptDoReadAction(DEFCPT look, int level)
{ int i, j, n, dt = 0 , dt1 = 0 ;
  char *cp, *cq = 0 ;
  KEY opt, key ;
  OBJ obj = 0 ;
  KEYSET ks = 0 ;
  FILE *f = 0 ;
  BOOL subclone ;
  mytime_t tim, tim1 ;
  extern BOOL isInteractive ;
  BOOL oldInteractive = isInteractive ;
  char timeBuf[25] ;

  if (!defCptOut)
    defCptOut = stdout ; 

  isInteractive = TRUE ;
  tim  = tim1 = timeNow () ;
  if (diffaction)
    catText(diffaction, messprintf("Action %%%%1* // %s\n", timeShow(tim, timeBuf, 25))) ;

  look->manEntry = FALSE ;
  acemblyInit() ;
  cDNAAlignInit () ;
  while(TRUE)
    { subclone = FALSE ;
      if (freelevelselect(level, &opt, mkaction))
	{ cp = freekey2text(opt, mkaction) ;
	  while (cp && *cp && *cp != ' ')
	    fputc (*cp++, defCptOut) ;
	  switch(opt)
	    {
	    case (KEY)(-1):
	      goto done ;
	    case 1001: case 1002: case 1003:
	    case 1011: case 1012: case 1013:
	    case 1005: case 1006:
	    case 1021: case 1022: case 1023: case 1024:
	    case 1047: case 1048: case 1049: 
	    case 1050: case 1051: case 1052: case 1053: 
	    case 1062: case 1063: case 1064:
	    case 1031: case 1032: case 1038: case 1039:
	    case 1061:
	    case 1071: case 1073: case 1074: case 1075: case 1076: case 1077: case 1080: case 1081:
	    case 1090: case 1091: case 1092: case 1095: case 1096: case 1097: case 1098:
	    case 1101: case 1102: case 1111: case 1121: case 1131:
            case 1301: case 1302: case 1311: case 1313: 
	      fMapcDNADoSelect (opt - 1000, look->actif, look->taceActif) ;
	      break ;
	    case 1060:
	      ctfTraceExportKeySet (look->actif) ;
	      break ;
	    case 1160:
	      geneZoneExport (look->taceActif) ;
	      break ;
	    case 1070:
	      chromoRepeatsKeyset (look->actif) ;
	      break ;
	    case 1200:
	      saucisseTest (look->taceActif) ;
	      break ;
	    case 1210:
#ifndef NON_GRAPHIC
	      htilePrecompute (look->taceActif) ;
#endif
	      break ;
	    case 2001: case 2002: case 2003:
	    case 3001: case 3002:
	      fMapcDNADoSelect (opt - 1000, look->actif, look->taceActif) ;
	      break ;
            case 'q':
	      ks = query (0, "FIND Sequence __*") ;
	      for (i = 0 ; i < keySetMax(ks) ; i++)
	      if ((obj = bsUpdate(keySet(ks, i))))
		bsKill (obj) ;
	      keySetDestroy (ks) ;
	      ks = query (0, "FIND DNA __*") ;
	      for (i = 0 ; i < keySetMax(ks) ; i++)
		arrayKill (keySet(ks, i)) ;
	      keySetDestroy (ks) ;
	      goto done ;
	    case 19: subclone = TRUE ;   /* charge les dogs, pas au point */
	      /* fall thru to load */
	    case 'l':
	      freenext () ; cp = freepos() ;
	      if (!*cp)
		{ fprintf (defCptOut, " // Load requires a sequence name. I quit") ;
		  goto done ;
		}
	      else
		{ KEY last ;
		  if (!strncmp (cp, "-last", 5) &&
		      (last = lastAssembly ())) 
		    cp = name (last) ;
		  else if (!strncmp (cp, "-all", 4))
		    cp = " CLASS Read " ; /* SCF_File AND NOT Assembled_from AND NOT Vector AND NOT IS _*  */
		  fprintf (defCptOut, "  %s  ", cp) ;
		}
	      cq = (char *)messalloc(1256) ;
	      if (!strncasecmp (cp, "-active", 7))
		sprintf (cq, "-active") ;
	      else if (!subclone)
		sprintf (cq,">? Sequence %s", cp) ;
	      else
		sprintf (cq,">? Clone IS %s ; >Read", cp) ;
              printf (" // load got: ##%s##\n", cq) ;
	      if (look->def) /* donc load au milieu du fichier => stop lecture */
		{ fprintf (defCptOut, " // restarting ") ;
		  defCptRestart (look) ;
		}
	      i = dnaAlignLoad (look, cq) ;
	      messfree (cq) ;
	      if (i)
		fprintf (defCptOut, "  // %d sequences loaded ", i) ;
	      else
		fprintf (defCptOut, " // No Sequence, sorry") ;
	      break ;
	    case 'A':  /* add */
	      i = look->def ? keySetMax (look->def) : 0 ;
	      if (!i)
		break ;
	      cp = freepos() ;
	      fprintf (defCptOut, " %s ", cp) ;
	      i = dnaAlignLoad (look, cp) ;
	      fprintf (defCptOut, "  // %d sequences added ", i) ;
	      break ;
	    case 'n':  /* newScf */
	      i = dnaAlignLoad (look, 0) ; /* 0 => baseCallNewScf() */
	      fprintf (defCptOut, "  // %d new scf reads added ", i) ;
	      break ;
	    case 'g': /* get */
	      if (!look->def || !keySetMax (look->def))
		break ;
	      if (!freeint(&i))
		i = 2 ;
	      fprintf (defCptOut, " %d ", i) ;
	      look->whatDis = (unsigned char)i ;
	      defCptDoGetOligos(look) ;
	      if (!look->def) /* || keySetMax(look->def) < 2) */
		goto done ;
	      break ;
	    case 's': /* sort */
	      if (!look->def || keySetMax (look->def) < 2) /* no sort for less than 2 seq */
		break ;
	      if (!freeint(&i)) /* dans le fichier il faut 100*distancemax !! */
		i = 100 ;
	      fprintf (defCptOut, " %d // %d sequences", i, keySetMax(look->def)) ;
 	      look->dlimit = i ;
	      i = defCptDoTree(look) ;
	      if (i > 1)
		fprintf (defCptOut, "  // Planted %d trees", i) ;
	      else
		fprintf (defCptOut, "  // Planted %d tree", i) ;
	      break ;
	    case 'a': /* assemble */
	      if (!look->def || !keySetMax (look->def))
		break ;
	      freeint(&i) ;
	      fprintf (defCptOut, " %d ", i) ;
	      look->taux = i ;
	      defCptDoSegments(look) ;
	      break ;
	    case 'j': /* join */
	      if (!look->def || !keySetMax (look->def))
		break ;
	      if ((cp = freeword ()) && !strncasecmp (cp, "diag", 4))
		{ i = defCptDoJoinDiagonal (look) ;
		  fprintf (defCptOut, " diagonal") ;
		}
	      else
		i = defCptDoJoinSegments (look) ;
	      fprintf (defCptOut, "  // I made %d joins ", i) ;
	      break ;
	    case 'J': /* juxtapose */
	      if (!look->def || !keySetMax (look->def))
		break ;
	      defCptJuxtapose (look) ;
	      break ;
	    case 'F': /* Fix */
	      if (!look->def || !keySetMax (look->def))
		break ; 
	      if ((cp = freeword ()) && !strncasecmp (cp, "-new", 4))
		defCptDoFixSegConsensus(look, TRUE) ;
	      else
		defCptDoFixSegConsensus(look, FALSE) ;
	      break ;
	    case 'i':
	      freenext () ; cp = freepos () ;
	      if (!*cp)
		ks = query (0, ">?Read IS *") ;
	      else if (!strncasecmp (cp, "-subclone", 9))
		{ cp += 9 ;
		  ks = query (0, messprintf (">?Subclone ; >Read %s", cp)) ;
		}
	      else
		ks = query (0, messprintf (">?Read %s", cp)) ;
	      if (look->link) dnaAlignDoReInsertLoners (look->link, ks) ;
	      keySetDestroy (ks) ;
	      break ;

/*	    case 7:
	      if (!look->def || !keySetMax (look->def))
		break ;
	      printf("I filter the data\n") ;
	      i = 0 ; freeint(&i) ;
	      fprintf (defCptOut, " %d ", i) ;
	      j = look->dlimit ;
	      look->dlimit = i ;
	      defCptDoSuppChim(look) ;
	      look->dlimit = j ;
	      break ;
*/

	    case 'V':
	      ks = look->def ;
	      if (!ks || !keySetMax(ks))
		ks = look->actif ;
	      if (!ks || !keySetMax(ks))
		break ; 
	      {
		int levelOut = 0 ;
		BOOL force = FALSE ;
		FILE *ff  = 0 ;
		Stack vStack = 0 ;
		while ((cp = freeword()))
		  {
		    if (!strcmp (cp, "-f"))
		      force = TRUE ;
		    if (!strcmp (cp, "-o") &&
			(cp = freeword ()))
		      {
			if (*cp == '-')
			  freeback () ;
			else
			  {
			    if ((ff = filopen (cp, 0, "w")))
			      levelOut = freeOutSetFile (ff) ;
			  }
		      }
		  }
		if (!ff)
		  {
		    vStack = stackCreate (100000) ;
		    levelOut = freeOutSetStack (vStack) ;
		  }
		i = trackVector (ks, 0, force) ;
		if (vStack)
		  {
		    freeOutClose (levelOut) ;
		    parseBuffer (stackText (vStack, 0), 0) ;
		    stackDestroy (vStack) ;
		  }
		if (ff)
		  {
		    freeOutClose (levelOut) ;
		    filclose (ff) ;
		  }
	      }
	      fprintf (defCptOut, "  // %d vector clips found",i) ;
	      break ;
	    case 'Q':
	      ks = look->def ;
	      if (!ks || !keySetMax(ks))
		ks = look->actif ;
	      if (!ks || !keySetMax(ks))
		break ; 
	      if (!bsIsTagInClass(_VSequence, str2tag("Bad_quality")))
		fprintf (defCptOut, "  // Bad_quality is not declared in the Sequence model, sorry") ;
	      else
		{ 
		  i = trackBadQuality (ks) ;
		  fprintf (defCptOut, "  // Labelled %d bad quality sequences",i) ;
		}
	      break ;
	    case 'C': /* clip_on */
	      if (!look->def || !keySetMax (look->def))
		break ;
	      j = 0 ;
	      if ((cp = freeword()))
		switch (ace_upper(*cp))
		  { 
		  case 'G': case 'F': case 'E':
		    key = (KEY)(*cp) ;
		    fprintf (defCptOut, " %s ", cp) ; j = 0 ;
		    i = baseCallUnclipKeySet (look->link, look->def, key, &j) ;
		    fprintf (defCptOut, " // %d clips moved by an average of %d bases ",
			     i, i > 0 ? j/i : 0 ) ;
		    break ;
		  case 'M':
		    fprintf (defCptOut, " %s ", cp) ;
		    i = 0 ; freeint (&i) ; j = 0 ;
		    fprintf (defCptOut, "%d ", i) ;
		    n = baseCallClipContig2Max (look->link, i, &j) ;
		    fprintf (defCptOut, " // %d clips moved by an average of %d bases ",
			     n, n > 0 ? j/n : 0 ) ;
		    break ;
		  case 'T': /* Tile */
		  case 'D': /* Double Tile */
		    key = (KEY)(*cp) ;
		    fprintf (defCptOut, " %s ", cp) ;
		    j = 0 ;
		    i = baseCallTileContigs (look->def, key, &j) ;
		    fprintf (defCptOut, " // %d clips moved by an average of %d bases ",
			     i, i > 0 ? j/i : 0 ) ;
		    break ;
		  case 'C': /* Consolidate */
		    key = (KEY)(*cp) ;
		    fprintf (defCptOut, " %s ", cp) ;
		    j = 0 ;
		    i = abiFixDouble (look->link, look->def, &j) ;
		    fprintf (defCptOut, " // %d clips moved by an average of %d bases ",
			     i, i > 0 ? j/i : 0 ) ;
		    break ;
		  default:
		    baseCallUnclipKeySet (look->link, look->def, 0, &j) ;
		    fprintf (defCptOut, "%s // Usage to move Clips: Clip_on [Excellent | Good | Fair]", cp) ;
		  }
	      else
		{ i = baseCallUnclipKeySet (look->link, look->def, 0, &j) ;
		  fprintf (defCptOut, " // Clips Evaluated") ;
		}
	      break ;
	    case 10:
	      if (!look->def || !keySetMax (look->def))
		break ;
	      defCptDoExtendReadsExt (look) ;
	      break ;
	    case 'L': /* align <target> */    
	      if (!look->def || !keySetMax (look->def))
		break ;   
	      freenext () ; cp = freepos () ;
	      if (!*cp)
		{ fprintf (defCptOut, "// Usage: Align <reference_sequence>") ;
		  break ;
		}
	      if (!lexword2key(cp, &key, _VDNA) || /* get a seq with dna */
		  !lexword2key(cp, &key, _VSequence))  
		{ fprintf (defCptOut, "// Sorry, unknown reference_sequence %s", cp) ;
		  break ;
		}  
	      doAssembleAllTraces (key, look->def, 'n') ;
	      break ;
	    case 21:  /* order by size */
	      if (!look->link)
		break ;
	      ks = queryKey (look->link, ">Subsequence") ;
	      contigLengthSort (ks) ;
	      alignToolsAdjustLink (look->link, ks, 0) ;
	      keySetDestroy (ks) ;
	      break ;
	    case 22: /* order by subclone */
	      if (look->link)
		  defCptOrderBySubclones (look) ;
	      break ;
	    case 23: /* make subclones */
	      if (!look->def || !keySetMax (look->def))
		break ;
	      baseCallMakeSubclones (look->def) ;
	      break ;
	    case 24: /* make subsequences */
	      if (!look->taceActif || !keySetMax (look->taceActif))
		break ;
	      defCptMakeSubsequences (look) ;
	      break ;
	    case 25: /* make_dna_cover nn filename */
	      if (!look->taceActif || !keySetMax (look->taceActif))
		break ;
	      defCptMakeDnaCover (look->taceActif) ;
	      break ;
	    case 13:
	      if (!look->def || !keySetMax (look->def))
		break ;
	      i = defCptDoPatchKeySet (look->def) ;
	      fprintf (defCptOut, " // %d auto-editions performed ", i) ;
	      break ;
	    case 'b':   /* base call */
	      if (!look->taceActif || !keySetMax (look->taceActif))
		break ;
	      i = defCptDoBaseCall (look->taceActif) ;
	      fprintf (defCptOut, " // %d traces called ", i) ;
	      break ;
	    case 18:
	      if (!look->def || !keySetMax (look->def))
		break ;
	      defCptDoTrainNNonKeySet (look->def) ;
	      break ;
	    case 'S': /* Save_as */
	      if (!look->def || !keySetMax (look->def))
		break ;
	      look->step = 0 ;
	      cp = freeword () ; 
	      if (cp)
		{ fprintf (defCptOut, " %s ", cp) ;
		  cq = strnew(cp, 0) ;
		  dnaAlignSave (look, cq, TRUE) ; /* pas a moi => peut changer */
		}
	      else
		{ cq = dnaAlignSaveDefault (look) ;
		  fprintf (defCptOut, " %s ", cq) ;
		}
	      messfree(cq) ;
	      sessionDoSave (TRUE) ;
	      break ;
	    case 'R': /* Rename */
	      if (!look->def || !keySetMax (look->def))
		break ;
	      look->step = 0 ;
	      cp = freeword () ; 
	      if (!cp)
		break ;
	      fprintf (defCptOut, " %s ", cp) ;
	      cq = strnew(cp, 0) ;
	      { KEYSET ksq2, ksq = 
		  query (0, messprintf
			 ("{FIND Sequence IS %s || IS %s.*}", cq, cq)) ;
	        ksq2 = query (ksq, "> DNA") ;
		i = keySetMax(ksq) ;
		while (i--)
		  { if ((obj = bsUpdate(keySet(ksq,i))))
		    bsKill (obj) ;
		  }
		keySetDestroy (ksq) ;	
		i = keySetMax(ksq2) ;
		while (i--)
		  { monDnaForget (look, keySet(ksq2,i)) ;
		    arrayKill (keySet(ksq2,i)) ;
		  }
		keySetDestroy (ksq2) ;
	      }
	      dnaAlignSave (look, cq, TRUE) ; /* pas a moi => peut changer */
	      messfree(cq) ;
	      sessionDoSave (TRUE) ;
	      break ;
	    case 'T':  /* tace */
	      f = stdin ;
	      if ((cp = freeword())) /* -f filename */
		{ 
		  if (!strncmp(cp, "-f", 2))
		    if (!(cp = freeword()) || /* filename */
			!(f = fopen(cp,"r")))
		      { fprintf (defCptOut 
				 , "Sorry: I cannot find file \"%s\"\n"
				 , cp ? cp : ""
				 ) ;
			break ;
		      }
		}
	      freenext () ;  /* get parameters */
	      cq = strnew(freepos(), 0) ;
	      i = isInteractive ;
	      commandExecute (freesetfile (f, cq), FALSE, TRUE, defCptOut, 0, 3,0) ;
	      messfree (cq) ;
	      isInteractive = i ;
	      break ;
	    case 'r':
	      tStatus () ;
	      break ;
	    case 16:
	      defCptPreReport (look) ;
	      break ;
	    case 'B':
	      if ((cp = freeword()) && lexword2key (cp, &key, _VSequence))
		statisticsDoMakeErreur (key, look->def) ;
	      else
		statisticsMakeErreur (look->def) ;
	      break ;
	    default:
	      fprintf (defCptOut, " // Unknown command, sorry") ;
	      break ;
	    }
	  if (DEBUG)
	    { mysize_t aUsed, aMade, aAlloc, aReal ;
	      int nmessalloc, messalloctot, dnaMade, dnaTotal ;
	      nmessalloc = messAllocStatus (&messalloctot) ;
	      monDnaReport (&dnaMade, &dnaTotal) ;
	      arrayStatus (&aMade, &aUsed, &aAlloc, &aReal) ;
	      fprintf (defCptOut, " // nMessAlloc %d totMessAlloc %d kb nArr %ld totArr %ld kb, nMondna %d totMonDna %d kb\n",
		       nmessalloc, messalloctot/1024, (long)aMade, (long)aAlloc/1024, dnaMade, dnaTotal/1024 ) ;
	      fprintf (stderr, " // %ld assInserted, %ld assRemoved, %ld assFound, %ld assNotFound, %ld assBounce\n",
		       assInserted, assRemoved, assFound, assNotFound, assBounce) ;
	    }
	  timeDiffSecs (tim, timeNow (), &dt) ;
	  timeDiffSecs (tim1, timeNow (), &dt1) ;
	  fprintf (defCptOut, " // + %d = %d seconds\n", dt1, dt) ;
	  tim1 = timeNow () ;
	}
    }
 done:
  timeDiffSecs (tim, timeNow (), &dt) ;
  fprintf (defCptOut," // end of ace-mbly  %s, total time %d seconds\n\n",
	   timeShow(timeNow(), timeBuf, 25), dt) ;
  freeclose(level) ;
  look->manEntry = TRUE ;
  isInteractive = oldInteractive ;
  return ;
}

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptReadAction(void)
{ static char fileName[FIL_BUFFER_SIZE] , dirName[DIR_BUFFER_SIZE] ;
  FILE *ff ;
  char *param = 0 ;
  int level ;
  DEFCPTGET("defCptReadAction") ;

  strcpy(dirName, sessionFilName("", 0, 0)) ;
  ff = filqueryopen(dirName, fileName, "smb", "r",
		    "Command File ?") ;
  if (!ff)
    { 
      defcomputeDraw (look, 0) ;
      return ;
    }
  level = freesetfile(ff, param) ;
  defCptDoReadAction(look, level) ;
} /* defCptReadAction */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void defCptWriteAction(void)
{ FILE *ff ;
  static char fileName[FIL_BUFFER_SIZE], dirName[DIR_BUFFER_SIZE] ;

  if (!diffaction)
    return ;
  strcpy(dirName, sessionFilName("", 0, 0)) ;
  if (!(ff = filqueryopen(dirName, fileName, "smb", "w", 
			 "Name of Ace-mbly session ?")))
    return ;
  fprintf(ff, "%s", stackText(diffaction, 0)) ;
  filclose(ff) ;
} /* defCptWriteAction */
#endif /* !NON_GRAPHIC */

/***************************************************************/

#ifndef NON_GRAPHIC
static void pleaseAssemble(void)
{ DEFCPTGET ("Assemble") ;

  look->method = M_ASSEMB ;
  look->whatDis = 2 ; /* in this case number of oligos */
  look->taux = 2 ;
  sprintf(look->nboligo, "%d", (int)2) ;
  sprintf(look->choix, "%d", (int)3) ;
  *look->param = 0 ;
  diffaction = stackCreate(1024) ;
  graphActivate(look->defMapCtlGraph) ;
  graphHelp("Sequence_Assembly") ;
  graphRetitle("Sequence Assembly") ;
  defcomputeDraw (look, 0) ;
}

/***************************************************************/

static void pleaseFingerPrint(void)
{ DEFCPTGET ("FingPr") ;

  look->method = M_FINGPR ;
  graphHelp("Finger_Print_Mapping") ;
  graphRetitle("Finger Print Mapping") ;
  defcomputeDraw (look, 0) ;
}

/***************************************************************/

static void pleaseDefDup(void)
{ DEFCPTGET ("Def-Dup") ;

  look->method = M_DEFDUP ;
  /* look->whatDis initialise dans defCptGetData */
  graphHelp("Genetic_Mapping_by_Markers_and_Segments") ;
  graphRetitle("Genetic Mapping by Markers and Segments") ;
  defcomputeDraw (look, 0) ;
}

/***************************************************************/

static void pleaseAlignment(void)
{ DEFCPTGET ("Align") ;

  look->method = M_EDWARD ;
  graphHelp("Alignment_of_Genes") ;
  graphRetitle("Alignment of Genes") ;
  defcomputeDraw (look, 0) ;
}

/***************************************************************/

static MENUOPT qhpMenu[] = {
  { graphDestroy, "Quit"},
  { help, "Help"},
  { NULL, NULL }} ;

static MENUOPT mapOrAssembleMenu[] = {
  { pleaseAssemble, "Sequence Assembly"},
  { pleaseFingerPrint, "FingerPrint"},
  { pleaseDefDup, "Def/Dup genetic mapping"},
  { pleaseAlignment, "Gene Alignment"},
  { NULL, NULL }} ;
#endif

/***************************************************************/

#ifndef NON_GRAPHIC
static void defcomputeDraw (DEFCPT look, char *text)
{ int line = 3 ;

  if (!graphActivate(look->defMapCtlGraph))
    return ;
  graphPop() ;
  graphClear() ;
  graphText(text, 2, 2) ;
  graphButtons (qhpMenu, 0.25, 0.25, 60) ;
  switch (look->method)
    {
    case 0:
      graphButtons (mapOrAssembleMenu, 0.25, 2.25, 60) ;
      break ;
    case M_ASSEMB:
      switch(look->mapStatus & F_BOUT)
	{
	case 6:
	  goto ici ;
	case 5:
	  graphButton ("Fix Present Consensus", defCptFixSegConsensus,  3, 11.25) ;
	case 4:
	  if (look->rejected && keySetMax (look->rejected))
	    graphButton ("Re Insert Bad Sequences", defCptInsertBad, 25, 11.25) ;
	  goto ici ;
	case 3:case 2:
	  if (!(look->mapStatus & F_SHOWM))
	    graphButton ("Show Matrix",defMapDisplay, 18.25,5.25) ;
	  graphButton ("More Oligos", defCptGetMoreSequences, 23.25,7.25) ;
/*	  graphButton ("Make One Segment", defCptMakeOneSegment, 35.63, 7.25) ; */
	  graphButton ("Try All Assembly", defCptTryAll, 35.63, 7.25) ;
	  graphButton ("Try Spec Paires", defCptNewTryPaires, 25, 16.5) ;
	  graphButton ("Filter", defCptSuppChim, 53, 7.25) ;
	  if (!(look->mapStatus & F_ASSSEG))
	    { graphButton ("3: Assemble Sequences", defCptSegments, .25, 9.25) ;
	      graphText("Error Rate : ", 30, 9.25) ;
	      look->choixBox = graphTextEntry(look->choix, 4, 47, 9.25, 0) ;
	    }
	case 1:
	  if ((look->mapStatus & F_BOUT) < (unsigned char)3)
	    { graphButton ("2: Sort Sequences", defCptTree, 0.25,5.25) ;
	      graphText("Distance limit : ", 22, 5.25) ;
	      look->disBox = graphTextEntry(look->distance, 5, 43, 5.25, 0) ;
	      graphButton ("2.bis: Join Sequences", defCptJoinSegments, 52, 9.25) ;
	    }
/*        defCptChoixBox = graphButton ("Adjust Parameters", defCptChoix, 3.25,7.25) ;
	  graphBoxDraw(defCptChoixBox, BLACK, defchoix ? GREEN : WHITE) ; */
	  graphButton ("Histograms", defCptAssHisto, 28.25, 1.25) ;
	case 0:
	  if (!look->mapStatus)
	    { graphButton ("Read a File", defCptReadAction, 0.25, 6.25) ;
	      graphText("Name of Sequences to load :", 13, 6.25) ;
	      look->paramBox = graphTextEntry(look->param, 25, 40, 6.25, 0) ;
	      graphButton ("1: Get KeySet of Sequences", defCptGetOligos, 0.25,3.25) ;
	      graphText("Number of Oligos : ", 30, 3.25) ;
	      look->nboligBox = graphTextEntry(look->nboligo, 4, 50, 3.25, 0) ;
	    }
	ici:
	  if ((look->mapStatus & F_BOUT) > (unsigned char)3 && keySetMax(look->def) > 1)
	    { graphButton ("4.bis: Join Sequences", defCptJoinSegments, 52, 9.25) ;
	      graphButton ("4: Iterate the assembly", defCptGetOligos, 0.25,3.25) ;
	      graphText("Number of Oligos : ", 30, 3.25) ;
	      look->nboligBox = graphTextEntry(look->nboligo, 4, 50, 3.25, 0) ;
	    }	      
	  line = 12 ;
	  graphButton("Save Game", defCptWriteAction, 13, .25) ;
	  graphText ("Other tools", .25, line += 2) ;
	  /* graphButton ("Fit Sequences To Target", oldAssembleAllTraces, 2, line += 1.5) ; */
	  graphButton ("Fit Sequences To Target", newAssembleAllTraces, 2, line) ;
	  graphButton ("Make Statistics", defCptMakeStatErreur, 15, line + 1.5) ;
	  graphButton ("Count Group", defCptMakeStatBarreaux, 55, line + 1.5) ;
	  /*	  graphButton ("Check Oligos", statisticsTestOligo, 55, line + 3) ;*/
	  graphButton ("Fix Active Keyset", defCptFixActKeyset, 48, line) ;
	  graphButton ("Autoedit", defCptPatchKeySet, 36.5, line) ;
	  graphButton ("Cpt One Paire", defCptAsmbPaireActive, 2, line + 1.5) ;
	  if (look->method & 3)
	    graphButton ("Make SuperLink", defCptMakeSuperLink, 35, line + 1.5) ;
	  if (look->mapStatus & F_IS_FIXED)
	    { graphButton ("Save As", defCptSaveAs, 53, line + 1.5) ;
	    }	      
	  if ((look->mapStatus & F_BOUT) /*  == 4  */)
	    graphButton ("Extend Reads", defCptExtendReadsExt,  23, 12.75) ;
	  if (look->mapStatus & F_HIST)
	    graphButton("Clean Histograms", defCptRmHisto, 40.50, 1.25) ;
	  break ;
	}
      break ;
    case M_FINGPR:
      if (!look->mapStatus)
	{ graphButton("Get Clones", defCptGetClones, 2, 2) ;
	  break ;
	}
      if (!(look->mapStatus & F_SORT))
	graphButton("Sort Segment", defCptTree, 2, 3.5) ;
      graphText("Filter Functions", .25, 5) ;
      graphButton("Suppress overpresent bands", defTreeSupr, 2, 6) ;
      if (look->mapStatus & F_SORT)
	graphButton("Suppress off diagonal bands", defTreeNonDiag, 2, 7.5) ;
      if (look->whatDis & F_E_FLAG)
	graphButton("use all Data", defSupFlag, 33, 7.5) ;
      graphText("Display Functions", .25, 9) ;
      if (!(look->mapStatus & F_SHOWM))
	graphButton("Show Matrix", defMapDisplay, 2, 10) ;
      graphButton("Histograms", defCptFingPrHisto, 15, 10) ;
      if (look->mapStatus & F_HIST)
	graphButton("Clean Histograms", defCptRmHisto, 27, 10) ;
      break ;
    case M_DEFDUP:
      if (!look->mapStatus)
	{ graphButton ("Select a Map", defCptGetData, 2, 4) ;
	  break ;
	}
      if (look->mapStatus & F_NO_CONTIG)
	{ graphButton("Select a Contig", defCptGetData, 2, 4) ;
	  break ;
	}
      graphButton ("Change Contig", defCptGetData, 2, 4) ;
      if (!(look->mapStatus & F_SORT))
	graphButton("Sort Segment", defCptTree, 2, 5.5) ;
      else
	{ graphButton("Show Map", defCptExportOrder, 2, 5.5) ;
	  if (!(look->mapStatus & F_SOMAR))
	    { graphButton("Sort Marker", defCptSortMarkers, 50, 5.5) ;
	      graphButton("Old Sort Marker", defTreeSortM, 46, 12) ;
	    }
	}
      if (look->whatDis & F_E_YN)
	graphButton("Treatment of Unknown", defCptTreatUnk, 25, 5.5) ;
      graphText("Filter Functions", .25, 7) ;
      graphButton("Suppress trivial segments", defCptSupSing, 2, 8) ;
      graphButton("Suppress overpresent markers", defTreeSupr, 33, 8) ;
      if (look->mapStatus & F_SORT)
	graphButton("Suppress off diagonal markers", defTreeNonDiag, 2, 9.5) ;
      if (look->whatDis & F_E_FLAG)
	graphButton("use all Data", defSupFlag, 33, 9.5) ;
      graphText("Display Functions", .25, 11) ;
      if (!(look->mapStatus & F_SHOWM))
	graphButton("Show Matrix", defMapDisplay, 2, 12) ;
      graphButton("Histograms", defCptDefDupHisto, 15, 12) ;
      if (look->mapStatus & F_HIST)
	graphButton("Clean Histograms", defCptRmHisto, 27, 12) ;
      break ;
    case M_EDWARD:
      if (!(look->mapStatus))
	{ 
	  graphButton("Get Data", defCptGetGenes, 0.25, 3.25) ;
	}
      else if (look->mapStatus == 1)
	graphButton("Ace-Dump Allele Order", defCptExportOrder, 2, 3.5) ;
      break ;
    default:
      break ;
    }
  graphTextBounds (100, 17) ;
  graphRedraw() ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/
/* destruction des tableaux de look sauf dataArray (pour changer de contig) */
/* ? pour knownPairs et knownOrder */
static void lookArrayDestroy(DEFCPT look)
{ int i ;
  KEYSET *ksp ;
  TREE_DEF *seg ;

  look->nd = look->nm = 0 ;
  arrayDestroy (look->linOrder) ;
  arrayDestroy (look->colOrder) ;
  if (arrayExists (look->maillon))
    { if ((i = arrayMax (look->maillon)))
	{ seg = arrp (look->maillon, i - 1, TREE_DEF) + 1 ;
	  while(seg--, i--)
	    { keySetDestroy (seg->ks) ;
	      keySetDestroy (seg->x) ;
	    }
	}
      arrayDestroy (look->maillon) ;
    }
  if (arrayExists(look->defInMar))
    { if ((i = arrayMax(look->defInMar)))
	{ ksp = arrp(look->defInMar, 0, KEYSET) - 1 ;
	  while(ksp++, i--)
	    keySetDestroy(*ksp) ;
	}
      arrayDestroy (look->defInMar) ;
    }
  if (arrayExists (look->marInDef))
    { if ((i = arrayMax(look->marInDef)))
	{ ksp = arrp(look->marInDef, 0, KEYSET) - 1 ;
	  while(ksp++, i--)
	    keySetDestroy(*ksp) ;
	}
      arrayDestroy (look->marInDef) ;
    }
  arrayDestroy (look->gmap) ;
  keySetDestroy (look->def) ;
  keySetDestroy (look->mar) ;
  keySetDestroy (look->rejected) ;
}

/***************************************************************/

static void defCptRestart (DEFCPT look)
{ lookArrayDestroy(look) ;
  intrinsicTreeDestroy() ;
  dnaAlignDestroy(look) ;
#ifndef NON_GRAPHIC
  plotHistoRemove() ;
  if (graphActivate(look->defMapGraph))
    graphDestroy() ;
  if (graphActivate(look->defTreeGraph))
    graphDestroy() ;
#endif
  look->assDnaGet = assBigCreate (10000) ;
}

/***************************************************************/

static void uLocalDoDestroy (DEFCPT look)
{ char *kp = 0 ;
#ifndef NON_GRAPHIC
  Graph gActive = graphActive () ; /* some destruction affect Graphs ; we want keep the active one */
#endif

  if (!look || !look->magic)
    return ;
  if (look->magic != DEFCPTMAG)
    messcrash ("Destroy received a wrong look pointer") ;
  nblook-- ;
  look->magic = 0 ;
  lookArrayDestroy(look) ;
  intrinsicTreeDestroy() ;
  dnaAlignDestroy (look) ;
  keySetDestroy (look->actif) ;
  arrayDestroy(look->dataArray) ;
  arrayDestroy(look->knownPairs) ;
  arrayDestroy(look->knownOrder) ;
  if (assDefLook && look->link)
    assRemove (assDefLook, kp + look->link) ;
#ifndef NON_GRAPHIC
  plotHistoRemove () ;
  if (graphActivate (look->defMapGraph))
    graphDestroy() ;
  if (graphActivate (look->defTreeGraph))
    graphDestroy() ;
  graphActivate (gActive) ;	/* restore context */
#endif /* !NON_GRAPHIC */
  messfree (look) ;
  if (diffaction)
    stackDestroy(diffaction) ;
}

/***************************************************************/

#ifndef NON_GRAPHIC
static void localDestroy(void)
{ DEFCPTGET (look) ;

  localDoDestroy (look) ;
}
#endif /* !NON_GRAPHIC */

/***************************************************************/
/***************************************************************/
/* attention ne sortir de tace que en detruisant look (localDestroy) */

static DEFCPT defComputeNewLook (void)
{ DEFCPT look = (DEFCPT)messalloc(sizeof(struct DEFCPTSTUFF)) ;
  nblook++ ;
  look->magic = DEFCPTMAG ;
  look->id = ++nmrid ;
  look->Line = 8 ;
  look->dlimit = 100 ;
  look->manEntry = FALSE ;
  look->display = FALSE ;
  look->assDnaGet = assBigCreate (10000) ;
  sprintf(look->distance, "%f", (float)1.) ;
  return look ;
}

/***************************************************************/

void defComputeTace (int level, KEYSET ks)
     /* called by command.c (for tacembly) */
{ 
  char timeBuf[25] ;

  if (!sessionGainWriteAccess()) /* try to grab it */
    { 
      /* may occur is somebody else grabbed it */
      messout ("Sorry, you cannot gain write access") ;
      return ;
    }
  nonGraphicLook = defComputeNewLook () ;
  if (!defCptOut)
    defCptOut = stdout ;
  fprintf (defCptOut, " // start ace-mbly; %s\n", 
	   timeShow(timeNow(), timeBuf, 25)) ;
  nonGraphicLook->taceActif = ks ;
  nonGraphicLook->actif = keySetExists (ks) ? query (ks, "CLASS Sequence || CLASS mrna") : 0 ;
  defCptDoReadAction (nonGraphicLook, level) ;
  localDoDestroy (nonGraphicLook) ;

  freeclose (level) ;
} /* defComputeTace */

/***************************************************************/
#ifndef NON_GRAPHIC

static void defCptPick(int box)
{ DEFCPTGET ("pick") ;

  if (look->nboligBox && box == look->nboligBox)
    graphTextEntry(look->nboligo, 0, 0, 0, 0) ;
  else if (look->disBox && box == look->disBox)
    graphTextEntry(look->distance, 0, 0, 0, 0) ;
  else if (look->choixBox && box == look->choixBox)
    graphTextEntry(look->choix, 0, 0, 0, 0) ;
  else if (look->paramBox && box == look->paramBox)
    graphTextEntry(look->param, 0, 0, 0, 0) ;
  else
    defcomputeDraw (look, 0) ;
}

void defCompute (void)
{ DEFCPT look = 0 ;

  if (graphActivate (myGraph))
    { if (!graphAssFind (&DEFCPTMAG, &look)) 
	messcrash ("graph not found in %s", name) ; 
      if (look->magic != DEFCPTMAG) 
	messcrash ("%s received a wrong pointer", name) ;
      defcomputeDraw (look, 0) ;
      myLook = look ;
      return ;
    }
  look = defComputeNewLook () ;
  look->display = TRUE ;
  look->manEntry = TRUE ;
  look->defMapCtlGraph = myGraph =
    graphCreate (TEXT_SCROLL, "Mapping and Assembly",
				      0.60,0.0,0.65,0.27) ;
  graphRegister(DESTROY, localDestroy) ;
  graphRegister(PICK, defCptPick) ;
  graphAssociate (&DEFCPTMAG , look) ;
  graphMenu(qhpMenu) ;
  graphHelp("Mapping and Assembly") ;
  graphColor (BLACK) ;
  defcomputeDraw (look, 0) ;
  myLook = look ;
}

/***************************************************************/

void defCptOpen (KEY link)
{

  /* Initialize here rather than at file level, for WIN32;
	 However, where does this go during WinAce (graphical) execution? */
  if (!defCptOut)
    defCptOut = stdout ; 

  defCompute () ; /* reopen the window */
  if (myLook->link == link)
    return ;

  if (myLook->link && myLook->link != link &&
      !messQuery
      (messprintf
       ("You are already assembling %s, %s %s",
	name(myLook->link), 
	"should we restart on ", name(link))))
      return ;
  defCptChangeLook (myLook, myLook->link, link) ;
  myLook->def = queryKey (link, ">Subsequence ; >DNA") ;
  pleaseAssemble () ;
}

/***************************************************************/

   /* call with com == 0 to start up the command window */
KEY defCptExecuteCommand (KEY link, KEYSET def, KEYSET actif, char *com, char *param)
{ int level = 0 ;
  DEFCPT look = 0 ;
  KEY key = 0 ;

  if (!link)
    link = 1 ; /* to be associated */
  if (com && *com)
    level = freesettext(com, param) ;
/* else   ? level n'etait pas initialise et pourquoi le else ? 
    look->display = TRUE ; */

  look = defCptGetLook (link) ;
  lookArrayDestroy(look) ; /* restart */
  intrinsicTreeDestroy() ;
  if (def)
    look->def = def ;
  else if (actif)
    look->actif = actif ;
  else
    { look->def = queryKey (link, ">Subsequence ; >DNA") ;
/*    look->rejected = queryKey (link, ">Rejected ; >DNA") ; */
    }

  if (level)
    { key = look->link ;
      if ((look->def && arrayMax (look->def)) || (look->actif && arrayMax (look->actif)))
	defCptDoReadAction(look, level) ;
      if (key == 1)
	{ key = look->link ;
	  defCptDestroyLook (key) ;
	}
      else if (key != look->link)
	{ defCptChangeLook (look, key, look->link) ; /* normally done before */
	  key = look->link ;
	}
      return key ;
    }
  else
    return 0 ;
}
#endif
/***************************************************************/
/***************************************************************/
/***************************************************************/

static void defCptGiveDirection (KEY key1, KEY key2, KEY tag, int *sens1, int *sens2)
{ int i1, i2, x1, x2 ;
  Array dog1 = 0, dog2 = 0 ;
  OBJ obj1 = 0, obj2 = 0 ;
  BSunit *u1, *u2 ;

  dog1 = arrayCreate (20, BSunit) ;
  obj1 = bsCreate (key1) ;
  if (!obj1 || !bsFindTag (obj1, tag) || !bsFlatten (obj1, 2, dog1))
    { *sens1 = 0 ; *sens2 = 0 ;
      goto abort ;
    }
  dog2 = arrayCreate (20, BSunit) ;
  obj2 = bsCreate (key2) ;
  if (!obj2 || !bsFindTag (obj2, tag) || !bsFlatten (obj2, 2, dog2))
    { *sens1 = 0 ; *sens2 = 0 ;
      goto abort ;
    }
  x1 = 0 ;
  x2 = 0 ;
  i1 = arrayMax (dog1) / 2 ;
  u1 = arrp (dog1, 0, BSunit) ;
  while (i1--)
    { i2 = arrayMax (dog2) / 2 ;
      u2 = arrp (dog2, 0, BSunit) ;
      while (i2--)
	{ if (u1->k == u2->k)
	    { x1 += (u1+1)->i ;
	      x2 += (u2+1)->i ;
	    }
	  u2 += 2 ;
	}
      u1 += 2 ;
    }
  *sens1 = x1 > 0 ? 1 : (x1 < 0 ? -1 : 0) ; /* sens 1 puis 2 ils doivent se faire front */
  *sens2 = x2 < 0 ? 1 : (x2 > 0 ? -1 : 0) ;
 abort:
  bsDestroy (obj1) ;
  bsDestroy (obj2) ;
  arrayDestroy (dog1) ;
  arrayDestroy (dog2) ;
}

/***************************************************************/

static BOOL defCptPolarFirst (TREE_DEF ksm, Array bilan, KEY tag, int *dir)
{ int i, max, x1, x2 ;
  KEY key1, key2 ;

  i = keySetMax (ksm.ks) ;
  max = arrayMax (bilan) ;
  if (i < 2)
    { array (bilan, max++, BSunit).k = keySet (ksm.ks, 0) ;
      array (bilan, max, BSunit).i = 1 ;
      return FALSE ;
    }
  key1 = keySet (ksm.ks, 0) ;
  key2 = keySet (ksm.ks, 1) ;
  for (i = 1 ; i < keySetMax (ksm.ks) && keySet (ksm.x, i) ; i++) ;
  if (i < keySetMax (ksm.ks))
    key2 = keySet (ksm.ks, i) ;
  defCptGiveDirection (key1, key2, tag, &x1, &x2) ;
  *dir = x1 ;
  return TRUE ;
}

/***************************************************************/

static void defCptPolarNext (TREE_DEF ksm, Array bilan, KEY tag, int fath, int dir)
{ int i, j, k, max, xmax, dirf, dir1, dir2, x2 ;
  KEY key1, key2, done = 1 << 24 ;
  BOOL pass = TRUE ;

  i = keySetMax (ksm.ks) ;
  key1 = keySet (ksm.ks, fath) ;
  xmax = keySet (ksm.x, fath) ;
  dirf = dir ;
  k = fath ;
 again:
  for (j = k + 1 ; j < i ; j++)
    { x2 = keySet (ksm.x, j) ;
      if (x2 == done)
	continue ;
      if (x2 < xmax)
	{ if (keySet (ksm.x, k) != done)
	    { max = arrayMax (bilan) ;
	      array (bilan, max++, BSunit).k = key1 ;
	      array (bilan, max, BSunit).i = dirf ;
	      keySet (ksm.x, k) = done ;
	    }
	  return ;
	}
      key2 = keySet (ksm.ks, j) ;
      defCptGiveDirection (key1, key2, tag, &dir1, &dir2) ;
      if (!dirf)
	dirf = dir1 ;
      if (pass && (x2 > xmax))
	{ if (dirf * dir1 < 0)
	    defCptPolarNext (ksm, bilan, tag, j, - dir2) ;
	  continue ;
	}
      if (keySet (ksm.x, k) != done)
	{ max = arrayMax (bilan) ;
	  array (bilan, max++, BSunit).k = key1 ;
	  array (bilan, max, BSunit).i = dirf ;
	  keySet (ksm.x, k) = done ;
	  j = k ;
	  pass = FALSE ;
	  continue ;
	}
      if (x2 > xmax)
	{ defCptPolarNext (ksm, bilan, tag, j, dir2) ;
	  continue ;
	}
      k = j ;
      dirf = dir2 ;
      key1 = key2 ;
    }
  if (keySet (ksm.x, k) != done)
    { max = arrayMax (bilan) ;
      array (bilan, max++, BSunit).k = key1 ;
      array (bilan, max, BSunit).i = dirf ;
      keySet (ksm.x, k) = done ;
      if (pass)
	{ pass = FALSE ;
	  goto again ;
	}
    }
}

/***************************************************************/

static void defCptPolarity (DEFCPT look, KEY sub, KEY tag)
{ int i, dir ;
  TREE_DEF ksm ;
  Array bilan = 0 ;

  i = arrayMax (look->maillon) ;
  if (!i)
    return ;
  bilan = arrayCreate (2*keySetMax (look->def), BSunit) ;
  for (i = 0 ; i < arrayMax (look->maillon) ; i++)
    { ksm = arr (look->maillon, i, TREE_DEF) ;
      if (defCptPolarFirst (ksm, bilan, tag, &dir))
	defCptPolarNext (ksm, bilan, tag, 0, dir) ;
    }
  alignToolsAdjustLink (look->link, 0, bilan) ;
}

/***************************************************************/
/* attention look->def est ici un keySet de Sequence */
static BOOL defCptDoOrderKeySet (DEFCPT look, KEY tag, KEY sub)
{ int i, j, k, nmd ;
  Array temp = 0, *ap ;
  OBJ obj = 0 ;
  KEY *keyp ;
  KEYSET *ksp, ks, resu ;
  BSunit *up ;

  if (!keySetExists (look->def) || !(i = keySetMax (look->def)))
    return FALSE ;
  keySetSort (look->def) ; /* normalement fait dans la fonction appelante ? */
  keySetCompress (look->def) ;
  i = keySetMax (look->def) ;
  temp = arrayCreate (i, Array) ;
  array (temp, i - 1, Array) = 0 ;
  ap = arrp (temp, 0, Array) - 1 ;
  keyp = arrp (look->def, 0, KEY) - 1 ;
  look->mar = keySetCreate () ;
  while (ap++, keyp++, i--)
    { *ap = arrayCreate (20, BSunit) ;
      if ((obj = bsCreate (*keyp)) && bsFindTag (obj, tag) && 
	  bsFlatten (obj, 1, *ap))
	{ j = arrayMax (*ap) ;
	  up = arrp (*ap, 0, BSunit) - 1 ;
	  k = keySetMax (look->mar) ;
	  while (up++, j--)
	    keySet (look->mar, k++) = up->k ;
	}
      bsDestroy (obj) ;
    }
  keySetSort (look->mar) ;
  keySetCompress (look->mar) ;
  if (!keySetMax (look->mar))
    { i = arrayMax (temp) ;
      ap = arrp (temp, 0, Array) - 1 ;
      while (ap++, i--)
	arrayDestroy (*ap) ;
      arrayDestroy (temp) ;
      keySetDestroy (look->mar) ;
      return FALSE ;
    }
  look->nd = keySetMax (look->def) ;
  look->nm = keySetMax (look->mar) ;
  look->whatDis = 1 ;
  look->method = M_DEFDUP ;
  i = look->nd ;
  look->marInDef = arrayCreate(i, KEYSET) ;
  array(look->marInDef, i - 1, KEYSET) = 0 ;
  i = look->nm ;
  look->defInMar = arrayCreate(i, KEYSET) ;
  array(look->defInMar, i - 1, KEYSET) = 0 ;
  ksp = arrp(look->defInMar, 0, KEYSET) ;
  while(i--)
    *ksp++ = keySetCreate() ;
  ap = arrp (temp, 0, Array) ;
  ksp = arrp(look->marInDef, 0, KEYSET) ;
  for (i = 0 ; i < look->nd ; ksp++, ap++, i++)
    { *ksp = keySetCreate () ;
      nmd = 0 ;
      j = arrayMax (*ap) ;
      if (j) /* mieg */
	{ up = arrp (*ap, 0, BSunit) - 1 ;
	  while (up++, j--)  
	    { if (!keySetFind (look->mar, up->k, &k))
		messcrash ("Le marqueur doit exister a ce niveau") ;
	      keySet (*ksp, nmd++) = k ;
	      ks = arr (look->defInMar, k, KEYSET) ;
	      keySet (ks, keySetMax (ks)) = i ;
	    }
	}
      arrayDestroy (*ap) ;
    }
  arrayDestroy (temp) ;
  i = look->nd ;
  ksp = arrp(look->marInDef, 0, KEYSET) - 1 ;
  while(ksp++, i--)
    arraySort(*ksp, mardefOrder) ;
  i = look->nm ;
  ksp = arrp(look->defInMar, 0, KEYSET) - 1 ;
  while(ksp++, i--)
    arraySort(*ksp, mardefOrder) ;
  arrayDestroy (look->colOrder) ;
  arrayDestroy (look->linOrder) ;
  TESCOLO ;
  TESLIN ;
  if (!intCptTree (look))
    return FALSE ;
  resu = keySetCopy (look->def) ;
  for (i = 0 ; i < look->nd ; i++)
    keySet (look->def, i) = keySet (resu, arr (look->colOrder, i, int)) ;
  keySetMax (look->def) = i ;
  if (look->link)
    defCptPolarity (look, sub, tag) ;
  keySetDestroy (resu) ;
  return TRUE ;
}

/***************************************************************/
/* restart look */
BOOL defCptOrderKeySet (KEYSET aTrier, KEY tag, KEY link, KEY sub)
{ BOOL end = FALSE, isCreate = FALSE ;
  KEY i = 1 ;
  char *kp = 0 ;
  DEFCPT look = 0 ;

  if (!keySetExists (aTrier) || !keySetMax (aTrier))
    return FALSE ;
  keySetSort (aTrier) ; /* normalement fait dans la fonction appelante ? */
  keySetCompress (aTrier) ;
  if (!link && assDefLook)
    while (assFind (assDefLook, kp + i, &look))
      i++ ;
  else i = link ;
  if (assDefLook && !assFind (assDefLook, kp + i, &look))
    isCreate = TRUE ;
  look = defCptGetLook (i) ;
  lookArrayDestroy (look) ; /* restart look */
  look->def = aTrier ;
  look->link = link ;
  if (defCptDoOrderKeySet (look, tag, sub))
    end = TRUE ;
  look->def = 0 ;
  if (isCreate)
    { look->link = i ;
      localDoDestroy (look) ;
    }
  else
    lookArrayDestroy (look) ;
  return end ;
}

/***************************************************************/

static void findMultipletsInContigs (KEY contig, KEY subcl, KEY mult)
{ OBJ Seq, Contig = bsUpdate (contig) ;
  int x1, x2, lng, old ;
  KEY seq, dog ;
  BSMARK mark = 0, mark2 = 0 ;
  
  if (!Contig)
    return ;
  if (bsFindTag (Contig, mult))
    bsRemove (Contig) ;

  if (bsGetKey (Contig, _Assembled_from, &seq))
    do
      { if ((Seq = bsCreate (seq)))
	  { if (bsGetKey (Seq, subcl, &dog))
	      { lng = 0 ;
		mark = bsMark (Contig, mark) ;
		if (bsGetData (Contig, _bsRight, _Int, &x1) &&
		    bsGetData (Contig, _bsRight, _Int, &x2))
		  lng = x2 - x1 ;
		bsAddKey (Contig, mult, dog) ;
		if (lng)
		  { mark2 = bsMark (Contig, mark2) ;
		    old = 0 ;
		    if (bsGetData (Contig, _bsRight, _Int, &old))
		      lng += old ; /* so conflicting directions average out */
		    bsGoto (Contig, mark2) ;
		    bsAddData (Contig, _bsRight, _Int, &lng) ;
		  }
		bsGoto (Contig, mark) ;
	      }
	    bsDestroy (Seq) ;
	  }
      } while (bsGetKey (Contig, _bsDown, &seq)) ;
  bsSave (Contig) ;
  bsMarkFree (mark) ;
  bsMarkFree (mark2) ;
}

/***************************************************************/

static void defCptMoveLocally (DEFCPT look, KEY tag, KEY sub)
{ int i, p1, p2, d1, d2, s1, s2, max = 0 ;
  KEY key1, key2 ;
  Array order = 0, unit = 0 ;
  OBJ Link = 0 ;
  BSunit *u ;

  Link = bsCreate (look->link) ;
  if (!Link)
    return ;
  unit = arrayCreate (90, BSunit) ;
  if (!bsFindTag (Link, sub) || !bsFlatten (Link, 3, unit) || !arrayMax (unit))
    goto abort ;
  key1 = arrp (unit, 0, BSunit)->k ;
  p1 = arrp (unit, 1, BSunit)->i ;
  p2 = arrp (unit, 2, BSunit)->i ;
  order = arrayCreate (90, BSunit) ;
  for (i = 3 ; i < arrayMax (unit) ; i += 3)
    { u = arrp (unit, i, BSunit) ;
      key2 = u[0].k ;
      d1 = u[1].i ;
      d2 = u[2].i ;
      defCptGiveDirection (key1, key2, tag, &s1, &s2) ;
      if ((s1 * (p2 - p1) < 0) && (s2 * (d2 - d1) < 0))
	{ array (order, max++, BSunit).k = key2 ;
	  array (order, max++, BSunit).i = d2 - d1 ;
	  continue ;
	}
      array (order, max++, BSunit).k = key1 ;
      array (order, max++, BSunit).i = p2 - p1 ;
      key1 = key2 ; p1 = d1 ; p2 = d2 ;
    }
  array (order, max++, BSunit).k = key1 ;
  array (order, max++, BSunit).i = p2 - p1 ;
  alignToolsAdjustLink (look->link, 0, order) ;
  arrayDestroy (order) ;
 abort:
  arrayDestroy (unit) ;
  bsDestroy (Link) ;
}

/***************************************************************/
/* attention change le contenu du look ; normalement appelee par execute command 
   avec un look vide a part look->def */
static void defCptOrderBySubclones (DEFCPT look)
{ int i ;

  keySetDestroy (look->def) ;
  if (look->link)
    look->def = queryKey (look->link, ">Subsequence") ;
  i = keySetMax (look->def) ;
  if (!i) return ;
  while (i--)
    findMultipletsInContigs (keySet (look->def, i), _Subclone, _Multiplet) ;
  if (defCptDoOrderKeySet (look, _Multiplet, _Subsequence))
    defCptMoveLocally (look, _Multiplet, _Subsequence) ;
}

/***************************************************************/

void defCptForget (KEY link, KEY key)
{ DEFCPT look = 0 ;
  char *kp = 0 ;

  if (!link)
    { if (assDefLook)
	while (assNext (assDefLook, &kp, &look))
	  dnaAlignForget (look, key) ;
      if (nonGraphicLook)
	dnaAlignForget (nonGraphicLook, key) ;
      if (myLook)
	dnaAlignForget (myLook, key) ;
    }
  else if (assDefLook && assFind (assDefLook, kp + link, &look))
    dnaAlignForget (look, key) ;
}

/***************************************************************/

void defCptChangeLook (DEFCPT look, KEY oldLink, KEY newLink)
{ char *kp = 0 ;
  KEY link ;

  if (!assDefLook)
    { assDefLook = assCreate () ;
      isTimeStamps = FALSE ;
    }
  if (oldLink)
    link = oldLink ;
  else
    link = look->link ;
  if (link) /* voir le cas pour le look static de la fonction */
    assRemove (assDefLook, kp + link) ;
  look->link = newLink ;
  assInsert (assDefLook, kp + newLink, look) ;
}

/***************************************************************/
/* Si look->link = 0 association avec la cle 1 qui ne peut etre un link (ancien look statique) */
DEFCPT defCptGetLook (KEY link)
{ DEFCPT look = 0 ;
  char *kp = 0 ;

  if (!assDefLook)
    { assDefLook = assCreate () ;
      isTimeStamps = FALSE ;
    }
  if (!link)
    { if (nonGraphicLook)
	look = nonGraphicLook ;
#ifndef NON_GRAPHIC
      else
	graphAssFind (&DEFCPTMAG, &look) ;
#endif      
      if (!look)
	messcrash ("Look non initialise") ;
    }
  else if (!assFind (assDefLook, kp + link, &look))
    { look = defComputeNewLook () ;
      assInsert (assDefLook, kp + link, look) ;
      look->link = link ;
    }
  if (look->magic != DEFCPTMAG)
    messcrash ("%s received a wrong pointer in defCptGetLook", "") ;
  return look ;
}

/***************************************************************/

void defCptDestroyLook (KEY link)
{ DEFCPT look = 0 ;
  char *kp = 0 ;

  if (!assDefLook)
    { isTimeStamps = FALSE ;
      return ;
    }
  if (link && assFind (assDefLook, kp + link, &look) && look)
    { 
#ifndef NON_GRAPHIC
      if (myLook && myLook->link == link) /* i e same look as defcpt's window */
	{ if (myLook != look)
	    messerror ("Should be the same look") ;
	  return ;
	}
#endif
      if (link != look->link)
	messerror ("Pb dans le look -> voir pour le destroy") ;
      if (myLook == look) myLook = 0 ;
      localDoDestroy (look) ;
    }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
 
