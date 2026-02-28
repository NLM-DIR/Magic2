/*  File: hseqdisp.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@ncbi.nlm.nih.gov
 *
 * Description: Display of the genetic map
 * Exported functions:
       hseqDisplay
 * HISTORY:

 * Created: July 15 2002 (mieg)
 *-------------------------------------------------------------------
 */

/*
#define ARRAY_CHECK
#define MALLOC_CHECK
*/

#include "acedb.h"
#include "bump.h"
#include "wac/ac.h"
#include "wfiche/gtitle.h"
#include "display.h"
#include "pick.h"
#include "session.h"
#include "systags.h"
#include "tags.h"
#include "classes.h"
#include "sysclass.h"
#include "bs.h"
#include "lex.h"
#include "a.h"
#include "graph.h"
#include "hseqdisp.h"
#include "query.h"
#include "cdna.h"

static AC_DB ac_db = 0 ;

typedef struct HMapStruct {
  magic_t *magic;        /* == HMap_MAGIC */
  AC_HANDLE h ;
  Graph graph ;

  float mag, yCurrent, leftMargin, offset ; 
  int graphWidth, graphHeight ;
  int max ;

  BOOL hideHeader ;

  int   activeBox ;
  int	messageBox ;
  char  messageText[128] ;

  int starBox ;
  char starBuffer [1000] ; /* a place to write the * when i touch an mrna */

  Array segs ;      	/* array of SEG's from the obj */
  Array introns, intronSupport ;
  Array boxIndex ;    /* box -> int (index into look->segs) */
  int lastTrueSeg ;	/* arrayMax(look->segs) after Convert */

  void (*mapDraw) (void) ;
} HMap ;

typedef struct HSeqStruct {
  magic_t *magic;        /* == HSeq_MAGIC */
  AC_HANDLE h ;
  AC_DB db ;
  HMap *map ;
  int type ; /* 1: gene, 2:mrna, 3:product */
  BOOL onlyNm, onlyTouched, allProbes, showNoSignalProbes, showOnlyExactProbes, showAmbiguousProbes ;
  KEY intMap, key, gene, tg, mrna, product, from ;
  int a1, a2 ; /* int map coords of the outside object */
  unsigned int flag ;
  int maxIntronClone ;
  DICT *tissues ;
  Array menu ;
} *HSeq;

static int  MAGGIE = -1 ;

static BOOL 
oldOnlyNm = FALSE
, oldOnlyTouched = FALSE
  , oldAllProbes = FALSE
  , oldShowNoSignalProbes = TRUE
  , oldShowOnlyExactProbes = FALSE
  , oldShowAmbiguousProbes = FALSE
  , oldHideHeader = FALSE 
  ;
typedef enum { Tzero, Tgene, Tmrna
	       , MrnaExon, VeryGoodCodingExon, GoodCodingExon, PoorCodingExon, uORFExon, Utr5, Utr3
	       , MrnaIntron, MrnaIntronFuzzy, GAP
	       , Valid5p, V5Capped, V5SL, V5Aggregate, V5stop
	       , Valid3p, V3aataaa, V3variant
	       , MrnaProbe, GENE
} HSTYPE  ;
char *hseqType[] = { "Tzero", "Tgene", "Tmrna"
		     , "MrnaExon", "VeryGoodCodingExon", "GoodCodingExon", "PoorCodingExon", "uORFExon", "Utr5", "Utr3"
		     , "MrnaIntron", "MrnaIntronFuzzy", "GAP"
		     , "Valid5p", "V5Capped", "V5SL", "V5Aggregate", "V5stop"
		     , "Valid3p", "V3aataaa", "V3variant"
		     , "MrnaProbe", "GENE"
} ;

#define NTAGS 8            /* how many supporting projects */
#define NTAGS_OUTSIDER 1    /* do not count the primates */
static  char *titleTags[] = { "UHR pooled cells", "C", "D", "Brain", "Blood", "Neuroblastoma", "Other", "Primates bodymap", 0} ;	

typedef struct
  { 
    KEY key, mrna ;
    int col, bgcol ;
    int a0, a1, a2, b1, b2, x1, x2 , nClo, nCloNM, nTagsA, nTagsC, nTagsD, nTagsB, nTagsBlood, nTagsNB, nTagsO , nTagsP, ln, bumpy ;
    int nTags[NTAGS] ;
    int sls, tissue ; /* dict  (hseq->tissues) */
    HSTYPE type, subType ;
    char txtType[12] ;
    int widgetBox, tableBox ;
    BOOL exactHit ;
    unsigned int flag ;
    float a,c,d,b,e,f ;
  } SEG ;

static BOOL hasTagCount = FALSE ;

static magic_t HSeq_MAGIC = "HSeq";
static magic_t HMap_MAGIC = "HMap";
static BOOL hseqConvert (HSeq look) ;

/************************************************************/
/************************************************************/
#define HMAP2GRAPH(_hMap,_x) (((_x) - (_hMap)->offset) * (_hMap)->mag + (_hMap)->leftMargin)

/************************************************************/

static HMap *hMapCurrent (char *caller)
{
  HMap *map = 0 ;

  if (!graphAssFind (&HMap_MAGIC,&map))
    messcrash("%s() could not find HMap on graph", caller);
  if (!map)
    messcrash("%s() received NULL HSeq pointer", caller);
  if (map->magic != &HMap_MAGIC)
    messcrash("%s() received non-magic HMap pointer", caller);
  
  return map ;
} /* hMapCurrent */

/************************************************************/

static HSeq currentHSeq (char *caller)
{
  HSeq hseq;

  if (!graphAssFind (&HSeq_MAGIC,&hseq))
    messcrash("%s() could not find HSeq on graph", caller);
  if (!hseq)
    messcrash("%s() received NULL HSeq pointer", caller);
  if (hseq->magic != &HSeq_MAGIC)
    messcrash("%s() received non-magic HSeq pointer", caller);
  
  return hseq;
} /* currentHSeq */

/************************************************************/

static void hMapResize (void)
{
  HMap *map = hMapCurrent ("hmapPick") ;
  HSeq look = currentHSeq("hseqDraw") ;

  pickRememberDisplaySize (DtHSEQ) ;

  graphFitBounds (&map->graphWidth,&map->graphHeight) ;
  map->leftMargin = 2 ;
  map->mag = ( (map->graphWidth - map->leftMargin -4))/map->max ;
  map->offset = 0 ;
  if (look->from)
    {
      map->mag = ( (map->graphWidth - map->leftMargin -4))/10000 ;
      map->offset = look->from ;
      look->from = 0 ;
    }
  map->mapDraw () ;
} /* hMapResize */

/************************************************************/

static void hMapKbd (int k)
{ 
  HMap *map = hMapCurrent ("hmapPick") ;

  if (map)
    switch (k)
      {
      case LEFT_KEY :
	map->offset -= map->graphWidth/(4 * map->mag) ;
	break ;
      case RIGHT_KEY :
	map->offset += map->graphWidth/(4 * map->mag) ;
	break ;
      case UP_KEY :
	map->mag *= 1.414 ;
	break ;
      case DOWN_KEY :
	map->mag /= 1.414 ;
	break ;
      case HOME_KEY :
	hMapResize () ;
	return ;
      default:
	return ;
      }
  map->mapDraw () ;
}

/************************************************************/

static void hMapFollow (HMap *map, int box)
{
  SEG *seg ;
  int iSeg ;
  if (box > 0 &&
      box < arrayMax(map->boxIndex) &&
      (iSeg = arr (map->boxIndex, box, int)) &&
      iSeg > 0 && 
      (seg = arrp(map->segs, iSeg, SEG)))
    display (seg->key, 0, TREE) ;
} /* hMapFollow */

/************************************************************/

static void hMapSelect (HMap *map, int box)
{
  SEG *seg = 0 ;
  int iSeg ;

  if (map->activeBox > 0 &&
      map->activeBox < arrayMax(map->boxIndex) && 
      (iSeg = arr (map->boxIndex, map->activeBox, int)) &&
      iSeg > 0 &&
      (seg = arrp(map->segs, iSeg, SEG)))
    {
      if (seg->widgetBox > 0 &&
	  seg->widgetBox < arrayMax(map->boxIndex))
	graphBoxDraw (seg->widgetBox, seg->col, seg->bgcol) ;
      if (seg->tableBox > 0 &&
	  seg->tableBox < arrayMax(map->boxIndex))
	graphBoxDraw (seg->tableBox, BLACK, WHITE) ;
    }
  
  map->activeBox = 0 ; seg = 0 ;
  if (box > 0 &&
      box < arrayMax(map->boxIndex) &&
      (iSeg = arr (map->boxIndex, box, int)) &&
      iSeg > 0 &&
      (seg = arrp(map->segs, iSeg , SEG)))
    { 
      map->activeBox = box ;
      graphBoxDraw (map->activeBox, BLACK, CYAN) ;
      if (seg->widgetBox > 0 &&
	  seg->widgetBox < arrayMax(map->boxIndex))
	graphBoxDraw (seg->widgetBox, BLACK, RED) ;
      if (seg->tableBox > 0 &&
	  seg->tableBox < arrayMax(map->boxIndex))
	graphBoxDraw (seg->tableBox, BLACK, CYAN) ;
     }

  if (seg && map->starBox)
    {
      if (map->activeBox && map->starBox)
	{
	  int x1, x0 = HMAP2GRAPH (map, seg->a1) ;

	  x1 = map->graphWidth ;
	  if (x0 > x1)
	    x0 = 1 ;
	  if (x1 > 999) x1 = 999 ;
	  memset (map->starBuffer, ' ', x1) ;
	  map->starBuffer [x0 - 1] = '*' ;
	  map->starBuffer [x1] = 0 ;
	}
      else
	map->starBuffer [0] = 0 ;

      graphBoxDraw (map->starBox, CERISE, -1) ;
    }

  if (map->activeBox && map->messageBox)
    { 
      strncpy (map->messageText, name (seg->key), 100) ;
      strcat (map->messageText, messprintf(" %d  %d", seg->x1, seg->x2)) ;
      
      graphBoxDraw (map->messageBox, BLACK, PALEBLUE) ;
    }
} /* hMapSelect */

/*************************************/

static void hMapPick (int box,  double x , double y) 
{
  HMap *map = hMapCurrent ("hmapPick") ;

  if (!box)
    return ;
  if (box < arrayMax (map->boxIndex) && arr(map->boxIndex,box,int)) /* a SEG */
    {
      if (box == map->activeBox) /* a second hit - follow it */
	hMapFollow (map, box) ;
      else
	hMapSelect (map, box) ;
    }

  return;
} /* hMapPick */

/************************************************************/

static HMap *hMapCreate (void *look, AC_HANDLE handle, void (*mapDraw) (void))
{
  HMap *map = (HMap *) halloc (sizeof (struct HMapStruct), handle) ;
  map->magic = &HMap_MAGIC ;
  map->h = handle ; /* do not destroy, belongs to calling routine */
  map->boxIndex = arrayHandleCreate (64,int, handle) ;
  map->segs = arrayHandleCreate (64, SEG, handle) ;
  map->introns = arrayHandleCreate (64, SEG, handle) ;
  map->intronSupport = arrayHandleCreate (64, SEG, handle) ;
  map->activeBox = 0 ;
  map->lastTrueSeg = arrayMax(map->segs) ;
  map->mapDraw = mapDraw ;

  map->hideHeader = oldHideHeader ;

  return map ;
} /* hMapCreate */

/************************************************************/

static void hMapGraphAssociate (HMap *map)
{
  graphRegister (KEYBOARD, hMapKbd) ;
  graphRegister (PICK, hMapPick) ;
  graphAssRemove (&HMap_MAGIC) ;
  graphAssociate (&HMap_MAGIC, map) ; /* attach map to graph */
  graphRegister (RESIZE, hMapResize) ;
  map->graph = graphActive() ;

  hMapResize () ; /* redraws */
  return ;
} /* hMapCreate */

/************************************************************/
/************************************************************/

static int  hseqSegOrder(const void *a, const void *b)  /* for arraySort() call */ 
{
  const SEG *sa = (const SEG*)a, *sb = (const SEG*) b ;
  int n ;

  n = sa->type - sb->type ;
  if (n) return n ;
  if (sa->mrna != sb->mrna)
    {
      n = sa->bumpy - sb->bumpy ;
      if (n) return n ;
      n = sa->a0 - sb->a0 ; 
      if (n) return n ;
      n = lexstrcmp(name (sa->mrna), name (sb->mrna)) ;
      if (n) return n ;
    }
  n = sa->a1 - sb->a1 ;
  if (n) return n ;
  n = sa->a2 - sb->a2 ;
  if (n) return n ;
  return sa->key - sb->key ;
}

/************************************************************/

static int  heatOrder(const void *a, const void *b)  /* for arraySort() call */ 
{
  const SEG *sa = (const SEG*)a, *sb = (const SEG*) b ;
  int n ;

  n = sa->x2 - sb->x2 ;
  if (n) return n ;
  return sa->key - sb->key ;
}

/************************************************************/
/***************** Registered routines *********************/

static void hseqDestroy (void)
{
  HSeq look = currentHSeq("hseqDestroy") ;

  if (look)
    {
      ac_free (look->h) ;
      
      graphAssRemove (&HSeq_MAGIC) ; /* detach look from graph */
      
      look->magic = 0 ;
      ac_free (look) ;
    }

  return;
} /* hseqDestroy */

/*************************************************************/

static void hseqRecalculate(void)
{
  HSeq look = currentHSeq("hseqRecalculate") ;
  
  if (look)
    {
      look->map->segs = arrayReCreate (look->map->segs,64, SEG) ;
      hseqConvert (look) ;
      look->map->mapDraw () ;
    } 
}

/*************************************/
/*************************************************************/
/************************************************************/

static int hseqExonsOrder (const void *a, const void *b)
{
  const SEG *sa = (const SEG *)a, *sb = (const SEG *)b ;
  int n ;

  n = sa->a1 - sb->a1 ;
  if (n) return n ;
  n = sa->a2 - sb->a2 ;
  return n ;
} /* hseqExonsOrder */

/************************************************************/

static int hseqShowSegs (Array segs)
{
  SEG *ss ;
  int ii;

  if (segs)
    for (ii = 0, ss = arrp (segs, 0, SEG) ; ii < arrayMax (segs) ; ss++, ii++)
      {
      printf ("%3d\t%12s\tbump=%d\ta=\t%5d %5d\tb=\t%5d %5d\tx=\t%5d %5d\tnclo=%5d\t%s\n"
	      , ii, hseqType[ss->type], ss->bumpy, ss->a1, ss->a2, ss->b1, ss->b2, ss->x1, ss->x2
	      , ss->nClo, name(ss->mrna)
	      ) ;
      }
  return 0 ;
} /* hseqExonsOrder */

/************************************************************/
static AC_OBJ  hseq_key_obj (HSeq look, KEY key)
{
  return ac_get_obj (look->db, className (key), name(key), look->h) ;
}

/************************************************************/

static BOOL hseqMrnaConvert2 (HSeq look, AC_OBJ mrna, int a0, AC_HANDLE h)
{
  BOOL ok = TRUE, exactOk ;
  AC_TABLE tbl, tbsl, tblx, spl = ac_tag_table (mrna, "Splicing", h) ;
  int ln, a1, a2, b1, b2, p1, p2, nClo, xx = 0, jr, ir, iSeg = arrayMax (look->map->segs)  ;
  int j1, jProduct, mrnaLength ;
  SEG *seg ;
  KEY probeKey = 0, product ;
  AC_OBJ Product = 0, Probe = 0 ;
  AC_KEYSET pks = 0 ;
  char buf [1000] ;
  const char *ccp ;
  BOOL isGood = FALSE ;
  BOOL isBest = FALSE ;
  BOOL isVeryGood = FALSE ;
  BOOL isUorf = FALSE ;

  hseqShowSegs (0) ; /* for compiler happiness */

  tbl = ac_tag_table (mrna, "Product", h) ; 
  for (jProduct = 0 ; jProduct < 3 ; jProduct++) 
    {
      isGood = isBest = isVeryGood = isUorf = FALSE ;
      /* draw up to 3 products per mRNA */
      j1 = p1 = p2 = 0 ; product = 0 ;
      for (ir = 0 ; tbl && tbl->cols > 2 && ir < tbl->rows ; ir++)
	{
	  isGood = isBest = isVeryGood = FALSE ;
	  Product = ac_table_obj (tbl, ir, 0, h) ; 
	  if (ac_has_tag (Product, "uORF_candidate") || ac_has_tag (Product, "Very_good_product") || ac_has_tag (Product, "Best_product"))
	    {
	      p1 = ac_table_int (tbl, ir, 1, 0) ; 
	      p2 = ac_table_int (tbl, ir, 2, 0) ; 
	      isVeryGood = ac_has_tag (Product, "Very_good_product")  || (ac_has_tag (Product, "Best_product") && ac_has_tag (Product, "Good_product")) ;
	      isGood = ac_has_tag (Product, "Good_product") ;
	      isBest = ac_has_tag (Product, "Best_product") ;
	      isUorf = !isGood && ! isBest && ac_has_tag (Product, "uORF_candidate") ;
	    }
	  product = ac_obj_key (Product) ;
	  ac_free (Product) ;
	  if (j1 == jProduct && (isVeryGood || isUorf ||isBest ))
	    break ;
	  j1++ ;
	}
      if (jProduct && j1 != jProduct) /* we found the first second ... very good product */
	continue ;
      if (isGifDisplay && jProduct && !isVeryGood && !isBest && !isUorf)
	continue ;
      if (isVeryGood && !isBest)
	isVeryGood = FALSE ;
      for (ir = 0 ; spl && ir < spl->rows ; ir++)
	{
	  if (!ir && !jProduct)
	    {
	      AC_KEYSET ks = ac_objquery_keyset (mrna, "! valid5p ; >Product ; good_product && best_product && up_stop", h) ;
	      if (0 && ac_keyset_count (ks)) 
		{
		  a1 = ac_table_int (spl, ir, 0, 0) ; 
		  b1 = 1 ;

		  seg = arrayp (look->map->segs, iSeg++, SEG) ;
		  seg->mrna = seg->key = ac_obj_key (mrna) ;
		  seg->type = Valid5p ;
		  seg->subType = V5stop ;
		  seg->b1 = b1 ; 
		  seg->b2 = b1 ;
		  seg->a1 = a0 + a1 - 1 ;
		  seg->a2 = a0 + a2 - 1 ; /* genomic coordinates */
		}
	      ac_free (ks) ;
	    }

	  if (!jProduct && strstr (ac_table_printable (spl, ir, 4, ""), "Gap"))
	    {
	      a1 = ac_table_int (spl, ir, 0, 0) ; 
	      a2 = ac_table_int (spl, ir, 1, 0) ; 
	      b1 = ac_table_int (spl, ir, 2, 0) ; 
	      b2 = ac_table_int (spl, ir, 3, 0) ; 
	      ln = a2 - a1 + 1 ;
	      seg = arrayp (look->map->introns, arrayMax(look->map->introns), SEG) ;
	      seg->mrna = seg->key = ac_obj_key (mrna) ;
	      seg->type = MrnaIntron ;
	      seg->subType = GAP ;
	      seg->b1 = b1 ; 
	      seg->b2 = b2 ;
	      seg->a0 = a0 ;
	      seg->a1 = a0 + a1 - 1 ;
	      seg->a2 = a0 + a2 - 1 ; /* genomic coordinates */
	      seg->ln = ln ;
	      xx += 100 ;
	    }
	  else if (!jProduct && strstr (ac_table_printable (spl, ir, 4, ""), "tron"))
	    {
	      a1 = ac_table_int (spl, ir, 0, 0) ; 
	      a2 = ac_table_int (spl, ir, 1, 0) ; 
	      b1 = ac_table_int (spl, ir, 2, 0) ; 
	      b2 = ac_table_int (spl, ir, 3, 0) ; 
	      ln = a2 - a1 + 1 ;
	      seg = arrayp (look->map->introns, arrayMax(look->map->introns), SEG) ;
	      seg->mrna = ac_obj_key (mrna) ;
	      seg->type = MrnaIntron ;
	      ccp = ac_table_printable (spl, ir, 5, 0) ; 
	      if (strstr (ccp, "uzzy"))
		{
		  ccp += 6 ;
		  seg->subType = MrnaIntronFuzzy ;
		}
	      else
		{
		  seg->subType = MrnaIntron ;
		}  
	      if (strcmp (ccp,"gt_ag") && strcmp (ccp,"gc_ag"))
		seg->flag = 1 ;
	      strcpy (seg->txtType+1, ccp) ;
	      seg->txtType[0] = '[' ;
	      seg->txtType[3] = '-' ;
	      seg->txtType[6] = ']' ;
	      
	      seg->b1 = b1 ; 
	      seg->b2 = b2 ;
	      seg->a0 = a0 ;
	      seg->a1 = a0 + a1 - 1 ;
	      seg->a2 = a0 + a2 - 1 ; /* genomic coordinates */
	      seg->ln = ln ;
	      xx += 100 ;
	    }
	  else if (strstr (ac_table_printable (spl, ir, 4, ""), "xon"))
	    {
	      a1 = ac_table_int (spl, ir, 0, 0) ; 
	      a2 = ac_table_int (spl, ir, 1, 0) ; 
	      b1 = ac_table_int (spl, ir, 2, 0) ; 
	      b2 = ac_table_int (spl, ir, 3, 0) ; 
	      ln = a2 - a1 + 1 ;
	      if (p1 > b1)
		{
		  xx = p1 <= b2 ? p1 - b1 : b2 - b1 + 1 ;
		  seg = arrayp (look->map->segs, iSeg++, SEG) ;
		  seg->mrna = seg->key = ac_obj_key (mrna) ;
		  seg->type = MrnaExon ;
		  seg->subType = Utr5 ;
		  seg->b1 = b1 ; 
		  seg->b2 = b1 + xx - 1 ;
		  seg->a0 = a0 ;
		  seg->a1 = a0 + a1 - 1 ;
		  seg->a2 = a0 + a1 - 1 + xx - 1 ; /* genomic coordinates */
		  seg->ln = ln ;
		  b1 = p1 ;
		  a1 = a1 + xx ;
		}
	      if (b1 <= b2 && b1 <= p2)
		{
		  xx = p2 < b2 ? b2 - p2 : 0 ;
		  seg = arrayp (look->map->segs, iSeg++, SEG) ;
		  seg->mrna = ac_obj_key (mrna) ;
		  seg->key = product ;
		  seg->type = MrnaExon ;
		  seg->subType = isVeryGood ? VeryGoodCodingExon : (isGood ? GoodCodingExon : (isUorf ? uORFExon : PoorCodingExon)) ;
		  if (isUorf) seg->subType = uORFExon  ;
		  seg->a0 = a0 ;
		  seg->a1 = a0 + a1 - 1 ;
		  seg->a2 = a0 + a2 - 1 - xx ; /* genomic coordinates */
		  seg->b1 = b1 ; 
		  seg->b2 = b2 - xx ;  
		  seg->ln = isUorf ? p2 - p1 + 1 : ln ;
		}
	      if (p2 < b2)
		{
		  xx = b1 < p2 ? b2 - p2 : b2 - b1 + 1;
		  seg = arrayp (look->map->segs, iSeg++, SEG) ;
		  seg->mrna = seg->key = ac_obj_key (mrna) ;
		  seg->type = MrnaExon ;
		  seg->subType = Utr3 ;
		  seg->a0 = a0 ;
		  seg->a2 = a0 + a2 - 1 ;
		  seg->a1 = a0 + a2 - 1 - xx + 1 ; /* genomic coordinates */
		  seg->b2 = b2 ;
		  seg->b1 = b2 - xx + 1 ; 
		  seg->ln = ln ;
		  b2 = p2 ;
		  a2 = a2 - xx ;
		}
	    }
	}
    }

  mrnaLength = ac_tag_int (mrna , "Covers", 0) ;
  tbl = ac_tag_table (mrna, "Valid3p", h) ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      a1 = ac_table_int (tbl, ir, 0, 0) ; 
      b1 = ac_table_int (tbl, ir, 1, 0) ; 

      if (a1 > mrnaLength)
	continue ;

      nClo = ac_table_int (tbl, ir, 2, 0) ; /* nb suppoting clones */

      seg = arrayp (look->map->segs, iSeg++, SEG) ;
      seg->mrna = seg->key = ac_obj_key (mrna) ;
      seg->type = Valid3p ;
      ccp = ac_table_printable (tbl, ir, 3, "") ;
      if (ccp && !strcmp (ccp, "AATAAA"))
	{ seg->subType = V3aataaa ; strcpy (seg->txtType, ccp) ; }
      else
	seg->subType = V3variant ;

      seg->b1 = b1 ; 
      seg->b2 = b1 ;
      seg->a1 = a0 + a1 - 1 ;
      seg->a2 = a0 + a1 - 1 ; /* genomic coordinates */
      seg->nClo = nClo ;
    }
  tbl = ac_tag_table (mrna, "Valid5p", h) ;
  /* disregard a valid5p if protein was constructed from AA==1 
   * it means that the protein was computed before we decided
   * that the mrna is 'aggregated' or 'capped'
   */
  if (tbl && ac_keyset_count (ac_objquery_keyset (mrna, ">product ; mRNA_5p_complete && best_product && !at_position_1", h))) 
    ac_free (tbl) ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      int isl = 0 ;

      a1 = ac_table_int (tbl, ir, 0, 0) ; 
      b1 = ac_table_int (tbl, ir, 1, 0) ; 
      nClo = ac_table_int (tbl, ir, 2, 0) ; /* nb suppoting clones */
      ccp = ac_table_printable (tbl, ir, 3, "") ;

      seg = arrayp (look->map->segs, iSeg++, SEG) ;
      seg->mrna = seg->key = ac_obj_key (mrna) ;
      seg->type = Valid5p ;
      tbsl = ac_tag_table (mrna, "Transpliced_to", h) ;
      for (isl = 0 ; tbsl && isl < tbsl->rows ; isl++)
	{
	  ccp = ac_table_printable (tbsl, isl, 0, "") ;
	  /* get table of clones right of SL0, get list of tissues */
	  if (!strcasecmp (ccp, "SL0"))
	    {
	      if (seg->subType != V5SL)
		seg->subType = V5Capped ;
	    }
	  else if (!strncasecmp (ccp, "SL", 2))
	    {
	      seg->subType = V5SL ;
	      if (seg->sls )
		{
		  if (strstr (dictName (look->tissues, seg->sls), ccp))
		    continue ;
		  ccp = messprintf ("%s, %s", dictName (look->tissues, seg->sls), ccp) ;
		}
	      if (! look->tissues)
		look->tissues = dictHandleCreate (64, look->h) ;
	      dictAdd (look->tissues, ccp, &seg->sls) ;
	    }
	  else if (ac_has_tag (mrna, "Aggregated_5p_clones"))
	    {
	      if (seg->subType != V5SL && seg->subType != V5Capped)
		seg->subType = V5Aggregate ;
	    }
	}
      
      seg->b1 = b1 ; 
      seg->b2 = b1 ;
      seg->a1 = a0 + a1 - 1 ;
      seg->a2 = a0 + a1 - 1 ; /* genomic coordinates */
      seg->nClo = nClo ;
    }

  tblx = ac_tag_table (mrna, "Probe_exact_hit", h) ;
  tbl = ac_tag_table (mrna, "Probe_hit", h) ;

  strcpy (buf, "IS * ") ;
  if (! look->showAmbiguousProbes)
    strcat (buf, messprintf (" && Confirmed_gene = \"%s\" ", name(look->gene))) ;
  if (! look->showNoSignalProbes)
    strcat (buf, " && Signal ") ;
 
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      probeKey = ac_table_key (tbl, ir, 0, 0) ; 
      Probe = ac_table_obj (tbl, ir, 0, 0) ; 
      if (ac_has_tag (Probe, "mRNA_many_hits"))
	continue ;
      p1 = ac_table_int (tbl, ir, 1, -999) ; 
      p2 = ac_table_int (tbl, ir, 2, -999) ; 
      if (p2 == -999)
	continue ;
      exactOk = FALSE ;
      for (jr = 0 ; tblx && ! exactOk && jr < tblx->rows ; jr++)
	if (ac_table_key (tblx, jr, 0, 0) == probeKey)
	  exactOk = TRUE ;
      if (look->showOnlyExactProbes && !exactOk)
	continue ;
      Probe = ac_table_obj (tbl, ir, 0, 0) ;      
      pks = ac_objquery_keyset (Probe, buf, h) ;
      if (ac_keyset_count (pks))
	{
	  seg = arrayp (look->map->segs, iSeg++, SEG) ;
	  seg->mrna = 0 ; /* would prevent arrayCompress ac_obj_key (mrna) ; */
	  seg->type = MrnaProbe ;
	  seg->key = ac_table_key (tbl, ir, 0, 0) ;  
	  seg->exactHit = exactOk ;
	  /* seg->x1 = p1 ;  seg->x2 = p2 ; */
	  for (jr = 0 ; spl && jr < spl->rows ; jr++)
	    {
	      if (strstr (ac_table_printable (spl, jr, 4, ""), "xon"))
		{
		  a1 = ac_table_int (spl, jr, 0, 0) ; 
		  a2 = ac_table_int (spl, jr, 1, 0) ; 
		  b1 = ac_table_int (spl, jr, 2, 0) ; 
		  b2 = ac_table_int (spl, jr, 3, 0) ; 
		  if (p1 <= b2 && p2 >= p1)
		    {
		      seg->a1 = a0 + a1 - 1 + p1 - b1 ; 
		      seg->a2 = a0 + a1 - 1 + p2 - b1 ; 
		      break ;
		    }
		}
	    }
	  Probe = ac_table_obj (tbl, ir, 0, 0) ;
	  seg->a = ac_tag_float (Probe, "Signal_A", -9999) ;
	  seg->c = ac_tag_float (Probe, "Signal_C", -9999) ;
	  seg->d = ac_tag_float (Probe, "Signal_D", -9999) ;
	  seg->b = ac_tag_float (Probe, "Signal_B", -9999) ;
	  seg->e = ac_tag_float (Probe, "Signal_E", -9999) ;
	  seg->f = ac_tag_float (Probe, "Signal_F", -9999) ;
	  {
	    KEY mycol ;

	    if ((mycol =  ac_tag_key (Probe, "Colour", 0)))
	      seg->col = WHITE + mycol - _WHITE ;
	    if ((mycol =  ac_tag_key (Probe, "Background_colour", 0)))
	      seg->bgcol = WHITE + mycol - _WHITE ;
	  }
	}
      ac_free (pks) ;
      ac_free (Probe) ;
    }
  return ok ;
} /* hseqMrnaConvert2 */

/************************************************************/

static BOOL  hseqTgConvert2 (HSeq look, AC_OBJ tg, int g1, int g2)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL ok = TRUE ;
  AC_TABLE tbl, introns ;
  AC_OBJ clo, est, mrna ;
  AC_KEYSET aks = 0 ;
  vTXT vtissues = 0 ;
  int a1 = 0, m1, ir, jr, iSeg, maxClo = 0 ;
  KEY type ;
  int x1, x2, y1, y2, n, nCloNM,  nTagsA, nTagsC, nTagsD, nTagsB, nTagsBlood , nTagsNB, nTagsO , nTagsP ;
  int iTags, nTags[NTAGS] ;

  SEG *seg ;

  tbl = ac_tag_table (tg, "IntMap", h) ;
  if (tbl && tbl->cols > 2)
    {
      a1 = ac_table_int (tbl, 0, 1, 0) ; 
    }
  
  if (g1 < g2)
    a1 = a1 - g1 + 1 ;
  else
    a1 = g1 - a1 + 1 ; 

  /*
    introns =  ac_aql_table (look->db, "select m, a1, a2, ntA, ntC, ntD, ntB, ntBlood, ntNB, nto, ntp  from tg in @active:1, ii in tg->Intron, m in ii->IntMap, a1 in m[1], a2 in a1[1],  ntA in ii->UHR, ntC in ii->Seqc_C, ntD in ii->Seqc_D, ntB in ii->Brain, ntBlood in ii->Blood, ntNB in ii->NB,  nto in ii->Other_sample, ntp in ii->Primates"
			   , ac_objquery_keyset (tg, "IS *", h), 0, 0, h) ;
 
    tbl = ac_aql_table (look->db, "select m, a1, a2 from tg in @active:1, m in tg->mrna, a1 in m[1], a2 in a1[1]"
		      , ac_objquery_keyset (tg, "IS *", h), "+2+3+1", 0, h) ;
  */
  introns =  ac_bql_table (look->db, "select m, a1, a2, ntA, ntC, ntD, ntB, ntBlood, ntNB, nto, ntp  from tg in @active:1, ii in tg->Intron, m in ii->IntMap, a1 in m[1], a2 in a1[1],  ntA in ii->UHR, ntC in ii->Seqc_C, ntD in ii->Seqc_D, ntB in ii->Brain, ntBlood in ii->Blood, ntNB in ii->NB,  nto in ii->Other_sample, ntp in ii->Primates"
			   , ac_objquery_keyset (tg, "IS *", h), 0, 0, h) ;
  tbl = ac_bql_table (look->db, "select m, a1, a2 from tg in @active:1, m in tg->mrna, a1 in m[1], a2 in a1[1]"
		      , ac_objquery_keyset (tg, "IS *", h), "+2+3+1", 0, h) ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      mrna = ac_table_obj (tbl, ir, 0, h) ; 
      if (look->onlyTouched && ! ac_has_tag (mrna, "Probe_hit"))
	continue ;
      if (look->onlyNm)
	{
	  aks = ac_objquery_keyset (mrna, ">cdna_clone ; IS NM*", h) ;
	  if (!ac_keyset_count (aks))
	    continue ;
	}
      m1 = ac_table_int (tbl, ir, 1, 0) ; 
      hseqMrnaConvert2 (look, mrna, a1 + m1 - 1, h) ;
    }

  iSeg = arrayMax (look->map->intronSupport)  ;
  tbl = ac_tag_table (tg, "Intron_boundaries", h) ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    {
      int jt = 0, di = 0 ;
      AC_TABLE reads = ac_db_empty_table (look->db, 64, 1, h) ;

      type = ac_table_key (tbl, ir, 0, 0) ;
      if (type == _Other) di = 1 ;
      x1 = ac_table_int (tbl, ir, 2+di, 0) ;
      x2 = ac_table_int (tbl, ir, 3+di, 0) ;
      n = nCloNM = nTagsA = nTagsC = nTagsD = nTagsB = nTagsBlood = nTagsNB = nTagsO = nTagsP = 0 ;
      for (iTags = 0 ; iTags < NTAGS ; iTags++)
	nTags[iTags] = 0 ;
      for (jr = ir ; jr < tbl->rows ; jr++)
	{
	  y1 = ac_table_int (tbl, jr, 2+di, 0) ;
	  y2 = ac_table_int (tbl, jr, 3+di, 0) ;
	  if (y1 != x1 || y2 != x2) break ;
	  n++ ;
	  clo = ac_table_obj (tbl, jr, 4+di, h) ;
	  est = ac_tag_obj (clo, "Read", h) ;
	  if (ac_has_tag (est, "Ref_Seq"))
	    nCloNM++ ;
	  if (est && ac_has_tag (est, "Tissue"))
	    ac_table_insert_type (reads, jt++, 0, &est, ac_type_obj) ;	    
	}
      for (jr = 0 ; introns && jr < introns->rows ; jr++)
	{
	  y1 = ac_table_int (introns, jr, 1, 0) ; /* intmap coord */
	  y2 = ac_table_int (introns, jr, 2, 0) ;
	  if (g1 < g2 && y1 < y2)
	    { y1 = y1 - g1 + 1 - a1 + 1 ; y2 = y2 - g1 + 1 - a1 + 1 ; }
	  else if (g1 > g2 && y1 > y2)
	    { y1 = g1 - y1 + 1 - a1 + 1 ; y2 = g1 - y2 + 1 - a1 + 1 ; }
	  else
	    continue ;	  
	  if (y1 == x1 && y2 == x2)
	    {
	      nTagsA = ac_table_int (introns, jr, 3, 0) ;
	      nTagsC = ac_table_int (introns, jr, 4, 0) ;
	      nTagsD = ac_table_int (introns, jr, 5, 0) ;
	      nTagsB = ac_table_int (introns, jr, 6, 0) ;

	      nTagsBlood = ac_table_int (introns, jr, 7, 0) ;
	      nTagsNB= ac_table_int (introns, jr, 8, 0) ;
	      nTagsO = ac_table_int (introns, jr, 9, 0) ;
	      nTagsP = ac_table_int (introns, jr, 10, 0) ;
	      if (nTagsB + nTagsA  + nTagsBlood  + nTagsNB + nTagsO  + nTagsP > 0)
		hasTagCount = TRUE ; /* in this database we have a tag count */

	      for (iTags = 0 ; iTags < NTAGS ; iTags++)  
		{
		  nTags[iTags] = ac_table_int (introns, jr, iTags + 3, 0) ;
		  if (nTags[iTags] > 0) 
		    hasTagCount = TRUE ; /* in this database we have a tag count */
		}
	      break ;
	    }
	}
      ir += n-1 ; /* those lines are consumed */
      seg = arrayp (look->map->intronSupport, iSeg++, SEG) ;
      seg->a1 = a1 + x1 - 1 ;
      seg->a2 = a1 + x2 - 1 ;
      seg->type = MrnaIntron ;
      seg->nClo = n ;
      seg->nCloNM = nCloNM ;

      seg->nTagsA = nTagsA ;
      seg->nTagsC = nTagsC ;
      seg->nTagsD = nTagsD ;
      seg->nTagsB = nTagsB ;

      seg->nTagsBlood = nTagsBlood ;
      seg->nTagsNB= nTagsNB ;
      seg->nTagsO = nTagsO ;
      seg->nTagsP = nTagsP ;

       for (iTags = 0 ; iTags < NTAGS ; iTags++)  
	 seg->nTags[iTags] = nTags[iTags] ;

      if (maxClo < n) maxClo = n ;
      if (jt)
	{
	  if (!vtissues)
	    vtissues = vtxtHandleCreate (h) ;
	  if (!look->tissues)
	    look->tissues = dictHandleCreate (64, look->h) ;
	  vtxtClear (vtissues) ;
	  /* 2011_01_25 : wrong reads is an AC_TABLE but should be an AC_KEYSET 
	    if ((ccp = ficheNewGeneExpressionTissue (vtissues, 0, reads, 2, 2, 0, seg->nClo)))
	    dictAdd (look->tissues, ccp, &seg->tissue) ;
	  */
	}
    }
  if (look->maxIntronClone < maxClo)
    look->maxIntronClone = maxClo ;

  ac_free (h) ;
  return ok ;
} /* hseqTgConvert2 */

/************************************************************/

static BOOL  hseqGeneConvert (HSeq look)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL ok = FALSE ;
  KEY map, gene = look->gene ;
  AC_OBJ Gene = 0 ;
  AC_TABLE tbl ;
  int ir, a1 = 0, a2 = 0, iSeg = 0, xxGap = 50 ;
  SEG *seg, *ss ;
  Array exons = arrayHandleCreate (64, SEG, h) ;
  int ii, iTags, xx = 0 ;

  Gene = hseq_key_obj (look, gene) ;
  if (!Gene)
    goto abort ;

  tbl = ac_tag_table (Gene, "IntMap", h) ;
  if (tbl && tbl->cols > 2)
    {
      seg = arrayp (look->map->segs, iSeg++, SEG) ;
      seg->key = gene ;
      seg->type = Tgene ;
      map = ac_table_key (tbl, 0, 0, 0) ;
      a1 = ac_table_int (tbl, 0, 1, 0) ; 
      a2 = ac_table_int (tbl, 0, 2, 0) ; 
      if (! look->intMap)
	{
	  look->intMap = map ;
	  look->a1 = a1 ;
	  look->a2 = a2 ;
	  seg->a1 = 1 ;
	  seg->a2 = a1 < a2 ? a2 - a1 : a1 - a2 ;
	}
      else if (map == look->intMap)  
	{
	  if (look->a1 < look->a2)
	    { seg->x1 = a1 - look->a1 + 1 ; seg->x2 = a2 - look->a1 + 1 ; }
	  else
	    { seg->x1 = look->a1 - a1 + 1 ; seg->x2 = look->a1 - a2 + 1 ; }
	}
      arrayMax (look->map->segs) -= 1 ; /* destroy this box after ! */
    }
  else
    goto abort ;
  tbl = ac_tag_table (Gene, "Transcribed_gene", h) ;
  look->maxIntronClone = 0 ;
  for (ir = 0 ; tbl && ir < tbl->rows ; ir++)
    hseqTgConvert2 (look, ac_table_obj (tbl, ir, 0, h), a1, a2) ;
  /* rationalize the introns and their support */
  for (iSeg = 0, seg = arrayp (look->map->introns, 0, SEG) ; iSeg < arrayMax (look->map->introns) ; iSeg++, seg++)
    {
      if (seg->type != MrnaIntron)
	continue ;
      a1 = seg->a1 ; a2 = seg->a2 ;
      for (ii = 0, ss = arrayp (look->map->intronSupport, 0, SEG) ; ii < arrayMax (look->map->intronSupport) ; ss++, ii++)
	{
	  if (ss->a1 == a1 && ss->a2 == a2)
	    { 
	      seg->nClo = ss->nClo ; seg->nCloNM = ss->nCloNM ;
	      seg->nTagsA = ss->nTagsA ;  seg->nTagsB = ss->nTagsB ; 
	      seg->nTagsC = ss->nTagsC ;  seg->nTagsD = ss->nTagsD ; 
	      seg->nTagsBlood = ss->nTagsBlood ;  seg->nTagsNB = ss->nTagsNB ; seg->nTagsO = ss->nTagsO ; seg->nTagsP = ss->nTagsP ; 

	      for (iTags = 0 ; iTags < NTAGS ; iTags++)
		seg->nTags[iTags] = ss->nTags[iTags] ; 
	      seg->tissue = ss->tissue ; break ;
	    }
	}
    }
  
  /* now squish the coordinates to remove the introns */
  /* 1) collect all exons */
  for (iSeg = ii = 0, seg = arrayp (look->map->segs, 0, SEG) ; iSeg < arrayMax (look->map->segs) ; iSeg++, seg++)
    {
      if (seg->type != MrnaExon)
	continue ;
      ss = arrayp (exons, ii++, SEG) ;
      ss->a1 = seg->a1 ;
      ss->a2 = seg->a2 ;
    }
  
  /* 2) merge exons touching or intersecting with each other */
  arraySort (exons, hseqExonsOrder) ;
  for (iSeg = 0, seg = arrayp (exons, 0, SEG) ; iSeg < arrayMax (exons) ; iSeg++, seg++)
    {
      if (seg->key)
	continue ;
      for (ii = iSeg + 1, ss = seg + 1 ; ii < arrayMax (exons) ; ii++, ss++)
	{
	  if (ss->a1 <= seg->a2 + 1)
	    {
	      ss->key = 1 ;
	      if ( ss->a2 > seg->a2)
		seg->a2 = ss->a2 ;
	      if ( ss->x2 > seg->x2)
		seg->x2 = ss->x2 ; /* closest end of previous exon */
	    }
	}
    }
  /* 3) compress */
  for (iSeg = ii = 0, ss = seg = arrayp (exons, 0, SEG) ; iSeg < arrayMax (exons) ; iSeg++, seg++)
    {
      if (!seg->key)
	{
	  if (ii < iSeg) *ss = *seg ;
	  ii++ ; ss++ ;
	}
    }
  arrayMax (exons) = ii ;
  /* 4) squish the coordinates */
  xxGap = 2000 ;  /* was 50, this is to make room to write the number of tags supporting the introns */
  for (xx = 1, iSeg = 0, seg = arrayp (exons, 0, SEG) ; iSeg < arrayMax (exons) ; iSeg++, seg++)
    {
      if (iSeg > 0 && seg->a1 > (seg-1)->a2 + 1)
	{ 
	  ss = arrayp (look->map->segs, arrayMax (look->map->segs), SEG) ;
	  ss->type = GAP ;
	  ss->a1 = (seg-1)->a2 + 1 ;
	  xx += xxGap ; 
	  ss->a2 = seg->a1 - 1 ;
	}
      seg->x1 = xx ;
      xx += seg->a2 - seg->a1 ;
      seg->x2 = xx ;
      xx++ ;
    }
  /* 4) rewrite the coordinates of each segment */
  for (iSeg = ii = 0, seg = arrayp (look->map->segs, 0, SEG) ; iSeg < arrayMax (look->map->segs) ; 
       iSeg++, seg++)
    {
      switch (seg->type)
	{
	case  MrnaProbe:
	case  MrnaExon:
	case Valid5p:
	case Valid3p:
	  for (ii = 0, ss = arrayp (exons, 0, SEG) ; ii < arrayMax (exons) ; ii++, ss++)
	    {
	      if (seg->a1 <= ss->a2 && seg->a2 >= ss->a1)
		{
		  seg->x1 = ss->x1 + seg->a1 - ss->a1 + 1 ;
		  seg->x2 = ss->x1 + seg->a2 - ss->a1 + 1 ;
		}
	      if (seg->a2 < ss->a1 -1)
		break ;
	    }
	  break ;
	default:
	  break ;
	}
    }
  for (iSeg = ii = 0, seg = arrp (look->map->introns, 0, SEG) ; iSeg < arrayMax (look->map->introns) ; 
       iSeg++, seg++)
    {
      int ok = 0 ;

      switch (seg->type)
	{
	case MrnaIntron:
	  for (ii = 0, ss = arrayp (exons, 0, SEG) ; ok < 2 && ii < arrayMax (exons) ; ii++, ss++)
	    {
	      if (seg->a1 >= ss->a1 - 1 &&
		  seg->a1 <=  ss->a2 + 1)
		{ ok++ ; seg->x1 = ss->x1 + seg->a1 - ss->a1 ; }
	      if (seg->a2 >= ss->a1 - 1 &&
		  seg->a2 <=  ss->a2 + 1)
		{ ok++ ; seg->x2 = ss->x1 + seg->a2 - ss->a1 ; }
	    }
	  if (ok < 2)
	    seg->type = 0 ;
	  break ;
	default:
	  break ;
	}
    }

  /* 5) find the new max */
  look->map->max = 0 ;
  for (iSeg = 0, seg = arrayp (look->map->segs, 0, SEG) ; iSeg < arrayMax (look->map->segs) ; 
       iSeg++, seg++)
    if (seg->type != Tgene && seg->x2 > look->map->max)
      look->map->max = seg->x2 ;
  look->map->max *= 1.1 ;

  ok = TRUE ;
 abort:
  ac_free (h) ;
  return ok ;
} /* hseqGeneConvert */

/************************************************************/

static BOOL  hseqTgConvert (HSeq look)
{
  look->tg = look->key ;
  look->key = look->gene = keyGetKey (look->tg, _Gene) ;
  return hseqGeneConvert (look) ;
}  /* hseqTgConvert */

/************************************************************/

static BOOL  hseqMrnaConvert (HSeq look)
{
  look->mrna = look->key ;
  look->key = look->tg = keyGetKey (look->mrna, _From_gene) ;
  return hseqTgConvert (look) ;
} /* hseqMrnaConvert */

/************************************************************/

static BOOL hseqProductConvert (HSeq look)
{
  look->product = look->key ;
  look->key = look->mrna = keyGetKey (look->product, _mRNA) ;
  return hseqMrnaConvert (look) ;		  
}  /* hseqProductConvert */

/************************************************************/

static BOOL hseqConvert (HSeq look)
{
  BOOL ok = FALSE ;

  if (class(look->key) == _VGene)
    { look->gene = look->key ; look->type = 1 ; }
  else if (class(look->key) == _VTranscribed_gene)
    { look->mrna = look->key ; look->type = 2 ;}
  else if (class(look->key) == _VmRNA)
    { look->mrna = look->key ; look->type = 3 ;}
  else if (class(look->key) == _VmProduct)
    { look->product = look->key ; look->type = 4 ; }
  else 
    goto abort ;
  
  switch (look->type)
    {
    case 1:
      ok = hseqGeneConvert (look) ;
      break ;
    case 2:
      ok = hseqTgConvert (look) ;
      break ;
    case 3:
      ok = hseqMrnaConvert (look) ;
      break ;
    case 4:
      ok = hseqProductConvert (look) ;
      break ;
    default:
      break ;
    }
  arraySort (look->map->segs,  hseqSegOrder) ;
  arrayCompress (look->map->segs) ;
 abort:
  return ok ;
}

/************************************************************/

static void hseqAllMrna (void)
{
  HSeq look = currentHSeq("hseqAllMrna") ;
  oldOnlyNm = look->onlyNm = FALSE ; 
  oldOnlyTouched =look->onlyTouched = FALSE ; 
  hseqRecalculate () ;
} /* hseqAllMrna */

/************************************************************/

static void hseqOnlyTouched (void)
{
  HSeq look = currentHSeq("hseqDestroy") ;
  oldOnlyTouched = look->onlyTouched = TRUE ;
  oldOnlyNm = look->onlyNm = FALSE ; 
  hseqRecalculate () ;
}

/************************************************************/

static void hseqJustNm (void)
{
  HSeq look = currentHSeq("hseqJustNm") ;
  oldOnlyNm = look->onlyNm = TRUE ;
  oldOnlyTouched = look->onlyTouched = FALSE ; 
  hseqRecalculate () ;
} /* hseqJustNm */

/************************************************************/

static void hseqNoSignalProbes (void)
{
  HSeq look = currentHSeq("hseqNoSignalProbes") ;
  oldShowNoSignalProbes = look->showNoSignalProbes = TRUE ;
  hseqRecalculate () ;
} /* hseqNoSignalProbes */

/************************************************************/

static void hseqJustSignalProbes (void)
{
  HSeq look = currentHSeq("hseqJustSignalProbes") ;
  oldShowNoSignalProbes = look->showNoSignalProbes = FALSE ;
  hseqRecalculate () ;
} /* hseqJustSignalProbes  */

/************************************************************/

static void hseqOnlyExactProbes (void)
{
  HSeq look = currentHSeq("hseqOnlyExactProbes") ;
  oldShowOnlyExactProbes = look->showOnlyExactProbes = TRUE ;
  hseqRecalculate () ;
} /* hseqOnlyExactProbes */

/************************************************************/

static void hseqAlsoInexactProbes (void)
{
  HSeq look = currentHSeq("hseqOnlyExactProbes") ;
  oldShowOnlyExactProbes = look->showOnlyExactProbes = FALSE ;
  hseqRecalculate () ;
} /*  hseqAlsoInexactProbes */

/************************************************************/

static void hseqAmbiguousProbes (void)
{
  HSeq look = currentHSeq("hseqAmbiguousProbes") ;
  oldShowAmbiguousProbes = look->showAmbiguousProbes = TRUE ;
  hseqRecalculate () ;
} /* hseqAmbiguousProbes */

/************************************************************/

static void hseqNoAmbiguousProbes (void)
{
  HSeq look = currentHSeq("hseqNoAmbiguousProbes") ;
  oldShowAmbiguousProbes = look->showAmbiguousProbes = FALSE ;
  hseqRecalculate () ;
}

/************************************************************/

static void hseqHideHeader (void)
{
  HSeq look = currentHSeq("hseqNoAmbiguousProbes") ;
  oldHideHeader = look->map->hideHeader = ! look->map->hideHeader ;
  hseqRecalculate () ;
}

/************************************************************/

static void hseqMenuSelect (HSeq look)
{
  int ii = 0 ;
  MENUOPT *mm ;
  
  arrayDestroy (look->menu) ;
  look->menu = arrayHandleCreate (30, MENUOPT, look->h) ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = graphDestroy ; mm->text = "Quit" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = help ; mm->text = "Help";

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = graphPrint ; mm->text = "Print" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = displayPreserve ; mm->text = "Preserve" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = hseqHideHeader ; mm->text = look->map->hideHeader ? "Show header" : "Hide header" ;

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = hseqRecalculate ; mm->text ="Recalculate" ;

  if (look->onlyNm || look->onlyTouched)
    {
      mm = arrayp (look->menu, ii++, MENUOPT) ;
      mm->f = hseqAllMrna ; mm->text = "All mRNAs" ;
    }

  if (!look->onlyTouched)
    {
      mm = arrayp (look->menu, ii++, MENUOPT) ;
      mm->f = hseqOnlyTouched ; mm->text = "Restrict to touched mRNAs" ;
    }
  if (! look->onlyNm)
    {
      mm = arrayp (look->menu, ii++, MENUOPT) ;
      mm->f = hseqJustNm ; mm->text = "Restrict to NM" ;
    }

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  if (look->showNoSignalProbes)
    {
      mm->f = hseqJustSignalProbes ; mm->text = "Restrict to probes with signal" ;
    }
  else
    {
      mm->f = hseqNoSignalProbes ; mm->text = "Add probes without signal" ;
    }

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  if (look->showAmbiguousProbes)
    {
      mm->f = hseqNoAmbiguousProbes ; mm->text = "Restrict to gene specific probes" ;
    }
  else
    {
      mm->f = hseqAmbiguousProbes ; mm->text = "Show non gene specific probes" ; 
    }

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  if (look->showOnlyExactProbes)
    {
      mm->f = hseqAlsoInexactProbes ; mm->text = "Add approximate probe hits" ;
    }
  else
    {
      mm->f =  hseqOnlyExactProbes ; mm->text = "Restrict to exact hits" ;
    }

  mm = arrayp (look->menu, ii++, MENUOPT) ;
  mm->f = 0 ; mm->text = 0;


} /* hseqMenuSelect */

/************************************************************/

static int hseqDrawAcdbSeg (HSeq look, SEG *seg, float bottomProbe, float dSignal, int ii, float *heightp)
{
  float  s1, s2, s, sm, scale, x, y ;
  int col, box = 0, box2, ix, iy, dx = 10, dy = 5, margin = 12 ;
  int nbOfCells = (look->map->graphWidth - margin)/dx ;

  if (MAGGIE) dy += 2 ;
  ix = ii % nbOfCells ;
  iy = ii/nbOfCells ;

  x = ix * dx + margin ;
  y = bottomProbe + iy * dy ;
  if (*heightp < (iy+1) * dy)
    *heightp = (iy+1) * dy ;
  if (! ix)
    {
      if (MAGGIE)
	{
	  graphBoxStart () ;
	  graphText ("C1", 1, y + 0) ;
	  graphText ("C2", 1, y + 1) ;
	  graphText ("C3", 1, y + 2) ;
	  graphText ("D1", 1, y + 3) ;
	  graphText ("D2", 1, y + 4) ;
	  graphText ("D3", 1, y + 5) ;
	  graphBoxEnd () ;
	}
      else
	{
	  graphBoxStart () ;
	  graphText ("Sample A", 1, y + 0) ;
	  graphText ("Sample C", 1, y + 1) ;
	  graphText ("Sample D", 1, y + 2) ;
	  graphText ("Sample B", 1, y + 3) ;
	  graphBoxEnd () ;
	}
	
    }

  s1 = +99999 ;
  if (seg->a > 0 && seg->a < s1) s1 = seg->a ;
  if (seg->c > 0 && seg->c < s1) s1 = seg->c ;
  if (seg->d > 0 && seg->d < s1) s1 = seg->d ;
  if (seg->b > 0 && seg->b < s1) s1 = seg->b ;
  if (s1 == 99999)
   goto done ;

  box = graphBoxStart () ;
  
  s2 = -99999 ;
  if (seg->a > 0 && seg->a > s2) s2 = seg->a ;
  if (seg->c > 0 && seg->c > s2) s2 = seg->c ;
  if (seg->d > 0 && seg->d > s2) s2 = seg->d ;
  if (seg->b > 0 && seg->b > s2) s2 = seg->b ;

  if (!MAGGIE && seg->a > 0 && seg->c > 0 && seg->d > 0 && seg->b > 0)
    {
      col = 0 ;
      if (seg->a >= seg->c &&  seg->c >= seg->d && seg->d >= seg->b)
	col = RED ;
      if (seg->a <= seg->c &&  seg->c <= seg->d && seg->d <= seg->b)
	col = GREEN ;
      if (col)
	{
	  graphColor (col) ;
	  graphRectangle (x -.2, y - .2, x + 9.2, y + 4.2) ;
	  graphRectangle (x -.1, y - .1, x + 9.1, y + 4.1) ;
	  graphColor (BLACK) ;
	}
    }
  if (!MAGGIE && ! seg->exactHit)
    graphText ("*", x + 8, y + 3.0) ;
  sm = (s1 + s2) / 2 ;
  scale = 4.0/(dSignal > 0 ? dSignal : 1) ;
  /* draw the 4 diagrams */

  if (MAGGIE)
    {
      if (seg->a > 0)
	{
	  s = 2.0 + scale * (seg->a - sm) ;
	  box2 = graphBoxStart () ;
	  graphRectangle (x+5, y + 0, x + 5 + s, y + .6 + 0) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box2, seg->col, seg->bgcol) ;
	  graphBoxSetPick (box2, FALSE) ;
	}
      if (seg->b > 0)
	{
	  s = 2.0 + scale * (seg->b - sm) ;
	  box2 = graphBoxStart () ;
	  graphRectangle (x+5, y + 1, x + 5 + s, y + .6 + 1) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box2, seg->col, seg->bgcol) ;
	  graphBoxSetPick (box2, FALSE) ;
	}
      if (seg->c > 0) 
	{
	  s = 2.0 + scale * (seg->c - sm) ;
	  box2 = graphBoxStart () ;
	  graphRectangle (x+5, y + 2, x + 5 + s, y + .6 + 2) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box2, seg->col, seg->bgcol) ;
	  graphBoxSetPick (box2, FALSE) ;
	}
      if (seg->d > 0) 
	{
	  s = 2.0 + scale * (seg->d - sm) ;
	  box2 = graphBoxStart () ;
	  graphRectangle (x+5, y + 3, x + 5 + s, y + .6 + 3) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box2, seg->col, seg->bgcol) ;
	  graphBoxSetPick (box2, FALSE) ;
	}
      if (seg->e > 0)
	{
	  s = 2.0 + scale * (seg->e - sm) ;
	  box2 = graphBoxStart () ;
	  graphRectangle (x+5, y + 4, x + 5 + s, y + .6 + 4) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box2, seg->col, seg->bgcol) ;
	  graphBoxSetPick (box2, FALSE) ;
	}
      if (seg->f > 0)
	{
	  s = 2.0 + scale * (seg->f - sm) ;
	  box2 = graphBoxStart () ;
	  graphRectangle (x+5, y + 5, x + 5 + s, y + .6 + 5) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box2, seg->col, seg->bgcol) ;
	  graphBoxSetPick (box2, FALSE) ;
	}
    }
  else
    {
      if (seg->a > 0)
	{
	  s = 2.0 + scale * (seg->a - sm) ;
	  box2 = graphBoxStart () ;
	  graphRectangle (x+5, y + 0, x + 5 + s, y + .6 + 0) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box2, seg->col, seg->bgcol) ;
	  graphBoxSetPick (box2, FALSE) ;
	}
      if (seg->c > 0) 
	{
	  s = 2.0 + scale * (seg->c - sm) ;
	  box2 = graphBoxStart () ;
	  graphRectangle (x+5, y + 1, x + 5 + s, y + .6 + 1) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box2, seg->col, seg->bgcol) ;
	  graphBoxSetPick (box2, FALSE) ;
	}
      if (seg->d > 0) 
	{
	  s = 2.0 + scale * (seg->d - sm) ;
	  box2 = graphBoxStart () ;
	  graphRectangle (x+5, y + 2, x + 5 + s, y + .6 + 2) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box2, seg->col, seg->bgcol) ;
	  graphBoxSetPick (box2, FALSE) ;
	}
      if (seg->b > 0)
	{
	  s = 2.0 + scale * (seg->b - sm) ;
	  box2 = graphBoxStart () ;
	  graphRectangle (x+5, y + 3, x + 5 + s, y + .6 + 3) ;
	  graphBoxEnd () ;
	  graphBoxDraw (box2, seg->col, seg->bgcol) ;
	  graphBoxSetPick (box2, FALSE) ;
	}
    }

  /* write out the numbers */
  box2 = graphBoxStart () ;
  if (MAGGIE)
    {
      if (seg->a > 0)
	graphText (messprintf ("%.2f", seg->a), x, y + 0) ; 
      if (seg->b > 0) 
	graphText (messprintf ("%.2f", seg->b), x, y + 1) ; 
      if (seg->c > 0) 
	graphText (messprintf ("%.2f", seg->c), x, y + 2) ; 
      if (seg->d > 0)
	graphText (messprintf ("%.2f", seg->d), x, y + 3) ; 
      if (seg->e > 0)
	graphText (messprintf ("%.2f", seg->e), x, y + 4) ; 
      if (seg->f > 0)
	graphText (messprintf ("%.2f", seg->f), x, y + 5) ; 
    }
  else
    {
      float sa = 0 ;
      if (seg->a > 0)
	graphText (messprintf ("%.2f", seg->a - sa), x, y + 0) ; 
      if (seg->c > 0) 
	graphText (messprintf ("%.2f", seg->c - sa), x, y + 1) ; 
      if (seg->d > 0) 
	graphText (messprintf ("%.2f", seg->d - sa), x, y + 2) ; 
      if (seg->b > 0)
	graphText (messprintf ("%.2f", seg->b - sa), x, y + 3) ; 
    }
  graphBoxEnd () ;
  graphBoxSetPick (box2, FALSE) ;
  graphBoxDraw (box2, seg->col, seg->bgcol) ;

  graphBoxEnd () ;
done:
  return box ;
} /* hseqDrawAcdbSeg */

/************************************************************/

static int hseqDrawAcdbSegs (HSeq look,float *heightp, float dSignal)
{
  int box, iSeg, iNN, y0 ;
  SEG *seg ;

  y0 = *heightp ;
  for (iSeg = iNN = 0, seg = arrp (look->map->segs, 0, SEG); iSeg < arrayMax(look->map->segs) ; iSeg++, seg++)
    {
      switch (seg->type)
	{
	case MrnaProbe: 	  
	  if ((box = hseqDrawAcdbSeg (look, seg, y0, dSignal, iNN, heightp)))
	    {
	      iNN++ ;
	      array (look->map->boxIndex, box, int) = iSeg ;
	      seg->tableBox = box ;
	    }
	  break ;
	default:
	  break ;
	}
    }
  return iNN ;
} /* hseqDrawAcdbSeg */

/************************************************************/

static void hseqDrawHeatCircles (HSeq look, float *yoffsetp)
{
  float y, s, a1, oldA1 = 0, oldy = 0, radius, dyMax, yScale ; 
  int jj = 0, col = BLACK, iSeg ;
  Array heatArray = arrayCreate (256, SEG) ;
  SEG *ss, *seg ;
  
  yScale = 2;
  dyMax = 3 ;
  *yoffsetp += 2 + yScale * dyMax ;
  if (1) graphLine (look->map->leftMargin - 2, *yoffsetp, look->map->graphWidth - 5, *yoffsetp) ; 
  
  graphColor (LAVANDER) ;
  for (iSeg = 0, seg = arrp (look->map->segs, 0, SEG); iSeg < arrayMax(look->map->segs) ; iSeg++, seg++)
    {
      if (seg->type != MrnaProbe || seg->a < -1000 || seg->b < -1000)
	continue ;
      y = seg->a - seg->b ;
      if (MAGGIE) 
	y = seg->a + seg->b + seg->c - seg->d - seg->e - seg->f ;
      if (y > dyMax) y = dyMax ;
      if (y < -dyMax) y = -dyMax ;
      y *= yScale ;
      a1 = HMAP2GRAPH (look->map, seg->a1) ;
      if (oldA1 == 0)
	oldA1 = a1 - 2 ;
      graphLine (oldA1,  *yoffsetp  - oldy, a1,  *yoffsetp - y) ;
      oldA1 = a1 ;
      oldy = y ;
      ss = arrayp (heatArray, jj++, SEG) ;
      *ss = *seg ;
      ss->x1 = 100 * y ;
      ss->x2 = y > 0 ? ss->x1 : - ss->x1 ;
    }
  if (oldA1 > 0)
    graphLine (oldA1,  *yoffsetp + 5  - oldy, oldA1 + 2, *yoffsetp + 5) ;
  
  arraySort (heatArray, heatOrder) ;
  
  for (iSeg = 0, seg = arrp (heatArray, 0, SEG); iSeg < arrayMax(heatArray) ; iSeg++, seg++)
    {
      if (seg->type != MrnaProbe || seg->a < -1000 || seg->b < -1000)
	continue ;
      y = seg->a - seg->b ;
      if (MAGGIE) 
	y = seg->a + seg->b + seg->c - seg->d - seg->e - seg->f ;
      if (y > 3) y = 3 ;
      if (y < -3) y = -3 ;
      a1 = HMAP2GRAPH (look->map, seg->a1) ;
      s = y > 0 ? seg->a : seg->b ;
      if (MAGGIE) s += 5 ;
      if (s < 7)
	col = GRAY ;
      else if (s < 8)
	 col = VIOLET ;
      else if (s < 9)
	col = LIGHTBLUE ;
      else if (s < 9.5)
	col = LIGHTCYAN ;
      else if (s < 10)
	col =LIGHTGREEN ;
      else if (s < 10.5)
	col = BROWN ;
      else if (s < 11)
	col = ORANGE ;
      else if (s < 12)
	col = LIGHTRED ;
      else if (s < 13)
	col = RED ;
      else 
	col = CERISE ;
      graphColor (col) ;
      if (
	  (seg->a > seg->c &&  seg->c > seg->d && seg->d > seg->b) ||
	  (seg->a < seg->c &&  seg->c < seg->d && seg->d < seg->b)
	  )
	{
	  radius = .9 ; ;
	  graphFillArc (a1, *yoffsetp - y * yScale, radius, 0, 360) ;
	}
      else
	graphCircle (a1, *yoffsetp - y * yScale, .6) ;
    }
  graphColor (BLACK) ;
  arrayDestroy (heatArray) ;
  *yoffsetp += 2 + yScale * dyMax ;
} /* hseqDrawHeatCircles */

/************************************************************/

static void hseqDrawSignalPlot (HSeq look, float *yoffsetp)
 {
  float x, y, a1, a2, da, dyMax, yScale, xScale ; 
  int i, jj = 0, col = BLACK, iSeg ;
  float olda1[6], oldy[6] ;
  SEG *seg ;
  BOOL isLog = TRUE ;

  if (isLog)
    {
      yScale = 15/1000.0;  /* was dyMax = 15 ; yScale = 1 ; */
      dyMax = 800 ;
    }
  else
    {
      yScale = 1 ;
      dyMax = 15 ;
    }
  a1 = a2 = 0 ;
  for (iSeg = 0, jj = 0, seg = arrp (look->map->segs, 0, SEG); iSeg < arrayMax(look->map->segs) ; iSeg++, seg++)
    {
      if (seg->type != MrnaProbe || seg->a < -1000 || seg->b < -1000)
	continue ;
      if (! jj++) { a1 = seg->a1 ; a2 = seg->a2 ;}
      if (seg->a1 < a1) a1 = seg->a1 ;
      if (seg->a2 > a2) a2 = seg->a2 ;
    }
  if (a1 == a2)
    return ;
  da = a2 - a1 ;
  a1 -= da/30 ; a2 += da/30 ;
  da = a2 - a1 ;
  
  (*yoffsetp) += 2 ;
  graphLine (HMAP2GRAPH (look->map,a1), *yoffsetp, 6, *yoffsetp + 5) ;
  graphLine (HMAP2GRAPH (look->map,a1) - .5, *yoffsetp, 6  - .5, *yoffsetp + 5) ;
  graphLine (HMAP2GRAPH (look->map,a2), *yoffsetp, look->map->graphWidth - 5, *yoffsetp + 5) ;
  graphLine (HMAP2GRAPH (look->map,a2) + .5, *yoffsetp, look->map->graphWidth - 4.5, *yoffsetp + 5) ;

  xScale = (look->map->graphWidth - 6  - 5)/da ;
  (*yoffsetp) += 5 ;

  (*yoffsetp) += 2 + yScale * dyMax ;
  graphLine (6 - 2, *yoffsetp, look->map->graphWidth - 5, *yoffsetp) ; 
  graphLine (6 - 2, *yoffsetp, 6  - 2, *yoffsetp - yScale * dyMax) ;
  if (! isLog)
    for (i = 0 ; i < dyMax ; i += 2)
      graphText (messprintf("%3d",i), 6 - 5, *yoffsetp - yScale * i - .5) ;
  else
    for (i = 0 ; i < dyMax ; i += 100)
      graphText (messprintf("%3d",i), 6 - 5, *yoffsetp - yScale * i - .5) ;
  graphColor (PALEVIOLET) ;
  if (! isLog)
    for (i = 2 ; i < dyMax ; i += 2)
      graphLine (6 - 2, *yoffsetp - yScale * i, look->map->graphWidth - 5, *yoffsetp - yScale * i) ; 
  else
    for (i = 100 ; i < dyMax ; i += 100)
      graphLine (6 - 2, *yoffsetp - yScale * i, look->map->graphWidth - 5, *yoffsetp - yScale * i) ; 
   
  for (iSeg = 0, jj = -1, seg = arrp (look->map->segs, 0, SEG); iSeg < arrayMax(look->map->segs) ; iSeg++, seg++)
    {
      if (seg->type != MrnaProbe || seg->a < -1000 || seg->b < -1000)
	continue ;
      jj++ ;
      if (0) /* show all lines */
	for (i = 0 ; i < (MAGGIE ? 5 : 4)  ; i++)
	  {
	    switch (i)
	      {
	      case 0: y = seg->a ; col = BLUE ; break ;
	      case 1: y = seg->b ; col = MAGGIE ? BLUE : RED ; break ;
	      case 2: y = seg->c ; col = MAGGIE ? BLUE : LIGHTBLUE ; break ;
	      case 3: y = seg->d ; col = MAGGIE ? RED : LIGHTMAGENTA   ; break ;
	      case 4: y = seg->e ; col = RED ; break ;
	      case 5: y = seg->f ; col = RED ; break ;
	      }
	    graphColor (col) ;
	    x = 6 + xScale * (seg->a1 - a1) ; 
	    if (y < 1) y = 1 ;
	    if (y > dyMax) y = dyMax ;
	    if (MAGGIE) y = exp (y * log(2)) ;
	    y = *yoffsetp - y * yScale ;
	    if (jj)
	      graphLine (olda1[i], oldy[i], x, y) ;
	    graphCircle (x, y,.5) ;
	    olda1[i] = x ; oldy[i] = y ;
	  }
      else /* show error bars */
	{
	  double z, z2, w, w2, a, b,c,d,e,f ;
	  float width = graphLinewidth (-1) ;
	  graphLinewidth (.2) ;
	  a = exp (seg->a * log(2)) ;
	  b = exp (seg->b * log(2)) ;
	  c = exp (seg->c * log(2)) ;
	  d = exp (seg->d * log(2)) ;
	  e = exp (seg->e * log(2)) ;
	  f = exp (seg->f * log(2)) ;
	  z = (a + b + c)/3.0 ;
	  w = (d + e + f)/3.0 ;
	  z2 = (a - z)*(a - z) +  (b - z)*(b - z) + (c - z)*(c - z) ;
	  w2 = (d - w)*(d - w) +  (e - w)*(e - w) + (f - w)*(f - w) ;
	  /* divide by (ddl - 1) since we measured the average */
	  z2 = sqrt (z2/2) ; w2 = sqrt (w2/2) ; 
	  x = 6 + xScale * (seg->a1 - a1) ; 
	  graphColor (BLUE) ;
	  graphLine (x, *yoffsetp - (w + w2) * yScale, x, *yoffsetp - (w - w2) * yScale) ;
	  graphLine (x - .6, *yoffsetp - (w + w2) * yScale, x + .6, *yoffsetp - (w + w2) * yScale) ;
	  graphLine (x - .6, *yoffsetp - (w - w2) * yScale, x + .6, *yoffsetp - (w - w2) * yScale) ;
	  graphFillArc (x, *yoffsetp - w * yScale, .6, 0, 360) ;
	  x += .2 ;
	  graphColor (RED) ;
	  graphLine (x, *yoffsetp - (z + z2) * yScale, x, *yoffsetp - (z - z2) * yScale) ;
	  graphLine (x - .6, *yoffsetp - (z + z2) * yScale, x + .6, *yoffsetp - (z + z2) * yScale) ;
	  graphLine (x - .6, *yoffsetp - (z - z2) * yScale, x + .6, *yoffsetp - (z - z2) * yScale) ; 
	  graphFillArc (x, *yoffsetp - z * yScale, .6, 0, 360) ;
	  graphLinewidth (width) ;
	}
    }
  if (MAGGIE) /* hard coded lines */
    {
      graphColor (MIDBLUE) ;
      for (x = 100 ; x < 400 ; x += 257)
	graphLine (6 + xScale * x, *yoffsetp - 700 * yScale, 6 + xScale * x, *yoffsetp) ;
    }
  graphColor (BLACK) ;
  /* draw scale */
  if (da > 10)
    {
      int i, dx1, dx2 ;
      dx1 = utMainRoundPart (da/10.0) ;
      dx2 = utMainPart (dx1/5.0) ;
      for (i = 0 ; dx2 * i < da ; i++)
	graphLine (6 + xScale * i * dx2, *yoffsetp, 6 + xScale * i * dx2, *yoffsetp + .3) ; 
      for (i = 0 ; dx1 * i < da ; i++)
	{
	  graphLine (6 + xScale * i * dx1, *yoffsetp, 6 + xScale * i * dx1, *yoffsetp + .8) ; 
	  graphText (messprintf ("%d", i * dx1), 5.5 + xScale * i * dx1, *yoffsetp + 1) ;
	}
      
    }
  

  *yoffsetp += 4 ;
 } /* hseqDrawSignalPlot */

/************************************************************/
/* little uitility to right flush a number */
static char *uTab(int w, int n)
{
  int i ;
  static char buf[256] ; 

  memset (buf, 0, 256) ;
  sprintf(buf, "%d", n) ;
  i = strlen (buf) ;
  if (i < w)
    {
      memset (buf,  ' ', w - i) ;
      sprintf(buf + w - i, "%d", n) ;
    }
  return buf ;
} /* uTab */

/************************************************************/

static void hseqDrawValidEnds (HSeq look, SEG *seg0, int iSeg, float y)
{
  AC_HANDLE h = ac_new_handle () ;
  int box, ii = iSeg + 1, iTags ;
  SEG *seg = seg0 + 1 ;
  float u1 ;
  KEY mrna = seg0->mrna ;
  vTXT vtxt = vtxtHandleCreate (h) ;

  for (ii = iSeg + 1, seg = seg0 + 1 ; ii < arrayMax (look->map->segs) ; ii++, seg++)
    switch (seg->type)
      {
      case MrnaExon: continue ;
      case Valid5p:
	if (seg->mrna != mrna) continue ;
	box = graphBoxStart () ;
	u1 = HMAP2GRAPH (look->map, seg->x1) ;
	if (seg->subType == V5stop)
	  graphColor (BLUE) ;
	graphLine (u1, y - 1, u1, y) ; 
	if (seg->subType == V5Capped || seg->subType == V5SL)
	  graphFillRectangle (u1, y - 1, u1 + .4, y - .6) ;
	else
	  graphRectangle (u1, y - 1, u1 + .4, y - .7) ;
	graphColor (BLACK) ;
	graphBoxEnd () ;
	switch (seg->subType)
	  {
	    /* seg->nClo = valid5p[2] */
	  case V5Capped: /* noir plein Transpliced_to      SL0 */
	    {
	      vtxtClear (vtxt) ;
	      vtxtPrintf (vtxt, "capped 5' end, %d accession%s\n"
			  , seg->nClo, _multi(seg->nClo)
			  ) ;

	      if (seg->tissue) vtxtPrintf (vtxt, "%s\n", dictName (look->tissues, seg->tissue)) ;

	      if (1)
		{
		  if (seg->nTagsA + seg->nTagsC  +seg->nTagsD + seg->nTagsB  + seg->nTagsBlood  + seg->nTagsNB + seg->nTagsO  + seg->nTagsP > 0)
		    vtxtPrintf (vtxt, "%s RNA-seq supporting reads\n", uTab (8, seg->nTagsB + seg->nTagsA  + seg->nTagsBlood  + seg->nTagsNB + seg->nTagsO)) ;
		  
		  if (seg->nTagsA) vtxtPrintf (vtxt, "%s UHR pooled cells\n",  uTab (8, seg->nTagsA)) ;
		  if (seg->nTagsC) vtxtPrintf (vtxt, "%s C\n",  uTab (8, seg->nTagsC)) ;
		  if (seg->nTagsD) vtxtPrintf (vtxt, "%s D\n",  uTab (8, seg->nTagsD)) ;
		  if (seg->nTagsB) vtxtPrintf (vtxt, "%s Brain\n",  uTab (8, seg->nTagsB)) ;
		  if (seg->nTagsBlood) vtxtPrintf (vtxt, "%s Blood\n",  uTab (8, seg->nTagsBlood)) ;
		  if (seg->nTagsNB) vtxtPrintf (vtxt, "%s Neuroblastoma\n",  uTab (8, seg->nTagsNB)) ;
		  if (seg->nTagsO) vtxtPrintf (vtxt, "%s Other\n",  uTab (8, seg->nTagsO)) ;
		  if (seg->nTagsP) vtxtPrintf (vtxt, "(Also %d primates bodymap)\n", seg->nTagsP) ;
		  vtxtPrintf (vtxt, "\n") ;
		}
	      if (1)
		{
		  int n = 0 ;
		  for (iTags = 0 ; iTags < NTAGS ; iTags++)
		    n += seg->nTags[iTags] ;
		  if (n)
		    vtxtPrintf (vtxt, "%s RNA-seq supporting reads\n", uTab (8,n)) ;
		  for (iTags = 0 ; iTags < NTAGS ; iTags++)
		    if (seg->nTags[iTags])
		      vtxtPrintf (vtxt, "%s %s\n"
				  ,  uTab (8, seg->nTags[iTags])
				  , titleTags[iTags]
				  ) ;
		  vtxtPrintf (vtxt, "\n") ;
		}
	      graphBubbleInfo (box, seg->mrna, className(seg->mrna), vtxtPtr (vtxt)) ;
	    }
	    break ;
	  case V5SL: /* noir plein Transpliced_to      SL1 */
	    {
	      vtxtClear (vtxt) ;
	      vtxtPrintf (vtxt, "Transpliced 5' end, %d accession%s\n"
			  , seg->nClo, _multi(seg->nClo)
			  ) ;
	
	      if (seg->tissue) vtxtPrintf (vtxt, "%s\n", dictName (look->tissues, seg->tissue)) ;
	      {
		if (seg->nTagsA + seg->nTagsC  +seg->nTagsD + seg->nTagsB + seg->nTagsBlood  + seg->nTagsNB + seg->nTagsO  + seg->nTagsP > 0)
		  vtxtPrintf (vtxt, "%s RNA-seq supporting reads\n", uTab (8, seg->nTagsB + seg->nTagsA  + seg->nTagsBlood  + seg->nTagsNB + seg->nTagsO)) ;
		if (seg->nTagsA) vtxtPrintf (vtxt, "%s UHR pooled cells\n",  uTab (8, seg->nTagsA)) ;
		if (seg->nTagsC) vtxtPrintf (vtxt, "%s C\n",  uTab (8, seg->nTagsC)) ;
		if (seg->nTagsD) vtxtPrintf (vtxt, "%s D\n",  uTab (8, seg->nTagsD)) ;
		if (seg->nTagsB) vtxtPrintf (vtxt, "%s Brain\n",  uTab (8, seg->nTagsB)) ;

		if (seg->nTagsBlood) vtxtPrintf (vtxt, "%s Blood\n",  uTab (8, seg->nTagsBlood)) ;
		if (seg->nTagsNB) vtxtPrintf (vtxt, "%s Neuroblastoma\n",  uTab (8, seg->nTagsNB)) ;
		if (seg->nTagsO) vtxtPrintf (vtxt, "%s Other\n",  uTab (8, seg->nTagsO)) ;
		if (seg->nTagsP) vtxtPrintf (vtxt, "(Also %d primates bodymap)\n", seg->nTagsP) ;
		vtxtPrintf (vtxt, "\n") ;
	      }
	      {
		int n = 0 ;
		for (iTags = 0 ; iTags < NTAGS ; iTags++)
		  n += seg->nTags[iTags] ;
		if (n)
		  vtxtPrintf (vtxt, "%s RNA-seq supporting reads\n", uTab (8,n)) ;
		for (iTags = 0 ; iTags < NTAGS ; iTags++)
		  if (seg->nTags[iTags])
		    vtxtPrintf (vtxt, "%s %s\n"
				,  uTab (8, seg->nTags[iTags])
				, titleTags[iTags]
				) ;
		vtxtPrintf (vtxt, "\n") ;
	      }
	      	      
	      graphBubbleInfo (box, seg->mrna, className(seg->mrna), vtxtPtr (vtxt)) ;
	    }
	    break ;
	  case V5Aggregate: /* noir creux mrna->Aggregated_5p_clones */
	    graphBubbleInfo (box, seg->mrna, className(seg->mrna)
			     , "5' end, aggregated 5'clones"
			     ) ;
	    break ;
	  case V5stop: /* bleu creux  */
	    graphBubbleInfo (box, seg->mrna, className(seg->mrna)
			     , "Putative 5' end, stop in UTR"
			     ) ;
	    break ;
	  default:
	    break ;
	  }
	break ;
      case Valid3p:
	if (seg->mrna != mrna) continue ;
	box = graphBoxStart () ;
	u1 = HMAP2GRAPH (look->map, seg->x1) ;
	if (seg->subType == V3variant)
	  graphColor (BLUE) ;
	graphLine (u1, y - 1, u1, y) ; 
	if (seg->nClo > 5)
	  graphFillRectangle (u1, y - 1, u1 - .4, y - .6) ;
	else
	  graphRectangle (u1, y - 1, u1 - .4, y - .7) ;
	graphColor (BLACK) ;
	graphBoxEnd () ;
	if (seg->nClo >= 0)
	switch (seg->subType)
	  {
	  case V3aataaa: /* noir plein */
	    graphBubbleInfo (box, seg->mrna, className(seg->mrna)
			     , messprintf("Validated 3' end, %d accession%s"
					  , seg->nClo
					  , _multi(seg->nClo)
					  ) 
			     ) ;
	    break ;
	  case V3variant: /* noir plein */
	    graphBubbleInfo (box, seg->mrna, className(seg->mrna)
			     , messprintf("Validated 3' end, %d accession%s"
					  , seg->nClo 
					  , _multi(seg->nClo)
					  ) 
			     ) ;
	    break ;
	  default:
	    break ;
	  }
	break ;
       default:
	 break ;
      }
  ac_free (h) ;
} /* hseqDrawValidEnds */

/************************************************************/
static char *isOne (int nn)
{
  static char buf[12] ;

  if (nn == 0) return "zero" ;
  if (nn == 1) return "one" ;
  sprintf (buf, "%d", nn) ;
  return buf ;
}

static int numIntronCol = 29 ;
static KEY intronCol[] = { PALEYELLOW , PALEGREEN, PALECYAN,  PALEMAGENTA, PALEVIOLET , 
			   RED1, PALEBLUE, PALEGRAY, GREEN1, PALEORANGE, BLUE1, 
			   PALERED, 
			   GREEN2, RED2, BLUE2,
			   GREEN3, RED3, BLUE3,
			   LIGHTCYAN, LIGHTORANGE, LIGHTMAGENTA, LIGHTVIOLET, 
			   GREEN4, RED4, BLUE4,
			   LIGHTRED, LIGHTGREEN, LIGHTGRAY, LIGHTBLUE
} ;
static void hseqDrawGene (HSeq look)
{
  AC_HANDLE h = ac_new_handle () ;
  int ii, iy, iymax, iSeg, iBound = 0, olda2 = 0, pass = 0, acdbHeight = 0, nClo ;
  SEG *seg, *oldSeg = 0 ;
  KEY oldMrna = 0 ;
  int box, iTags ; 
  float topMrna = 0, bottomProbe = 0 ;
  float x, y, y00 = 0, dy = 0, ddy = .6, boxHeight = .4, dSignal = 0, dyClo ;
  float u1, u2, oldu2 = 0, dummy;
  BOOL inProbes = FALSE ;
  char *pgmQ[]  = {"ABI",    "AFX",  "AGL"    ,  "GEH",     "ILM",   "NCI",  "EPP",     "GEX",  "QGN",  "TAQ", 0} ;
  int pgmcolQ[] = { ORANGE, YELLOW, PALEORANGE, VIOLET, LIGHTCYAN, MIDBLUE, PALEBLUE, PALERED, PALEGREEN, CERISE, 0} ;
  char *pgmJ[]  = {"Jun2",    "Jun1",  "Jun3", 0} ;
  char *pgmJ2[]  = {"b0805_9681",    "b0805_616",  "b0805_11137", 0} ;
  int pgmcolJ[] = {PALEBLUE, PALERED, PALEGREEN, 0} ;
  char **pgm, **pgm2  ;
  int *pgmcol ;
  BUMP bumper = 0 ;
  vTXT vtxt = vtxtHandleCreate (h) ;
  KEYSET boundaries = keySetHandleCreate (h) ;

  if (MAGGIE) { pgm = pgmJ ; pgm2 = pgmJ2 ; pgmcol = pgmcolJ ; }
  else { pgm = pgmQ ; pgm2 = pgmQ ; pgmcol = pgmcolQ ; }

  y = look->map->yCurrent ; iymax = -1 ;
  y += 2 ; /* top margin, room for the bubbles */
  dy = 2.3 ;
  for (pass = 0 ; pass < 3 ; pass++)
    {
      bumpDestroy (bumper) ;/*  iNN = 0 ; */
      bumper = bumpCreate (12,  0) ;
      for (iSeg = 0, seg = arrp (look->map->segs, 0, SEG); iSeg < arrayMax(look->map->segs) ; iSeg++, seg++)
	{
	  if (!pass)
	    seg->widgetBox = seg->tableBox = 0 ;
	  u1 = HMAP2GRAPH (look->map, seg->x1 - .5) ;
	  u2 =  HMAP2GRAPH (look->map, seg->x2 + .5) ;
	  if (u1 < 1) u1 = 1 ;
	  if (u2 >= look->map->graphWidth)
	    u2 = look->map->graphWidth - 1 ;
	  if (u1 > u2)
	    continue ;
	  iy = 0 ;
	  switch (seg->type)
	    {
	    case MrnaProbe:
	      if (pass == 0)
		continue ; 
	      else if (pass == 1)
		{
		  if (! inProbes || iymax == -1) /* first mrnaProbe */
		    {
		      y += dy ;
		    }
		  if (0 && seg->mrna != oldMrna) 
		    {
		      y += dy ;
		    } 
		  oldMrna = seg->mrna ;
		  if (! seg->col)
		    seg->col = BLACK ;
		  if (! seg->bgcol)
		    {
		      char *cp = name(seg->key) ;
		      char **cq = pgm - 1 ;
		      int *pc = pgmcol - 1 ;
		      
		      while (++pc, *++cq)
			if (!strncmp (cp, *cq, strlen(*cq)))
			  { seg->bgcol = *pc ; break ; }
		    }
		  if (! seg->bgcol)
		    seg->bgcol = GREEN ; 
		  boxHeight = .4 ;
		  inProbes = TRUE ;
		  dummy = u1 ; iy = 0 ;
		  bumpItem (bumper, 1, u2 - u1, &iy, &dummy) ;
		  if (iymax < iy) iymax = iy ;
		  if (bottomProbe < y + iy + .4 + boxHeight)
		    bottomProbe = y + iy + .4 + boxHeight ;
		  {
		    float s1, s2 ;
		    s1 = +99999 ;
		    if (seg->a > 0 && seg->a < s1) s1 = seg->a ;
		    if (seg->c > 0 && seg->c < s1) s1 = seg->c ;
		    if (seg->d > 0 && seg->d < s1) s1 = seg->d ;
		    if (seg->b > 0 && seg->b < s1) s1 = seg->b ;
		    
		    s2 = -99999 ;
		    if (seg->a > 0 && seg->a > s2) s2 = seg->a ;
		    if (seg->c > 0 && seg->c > s2) s2 = seg->c ;
		    if (seg->d > 0 && seg->d > s2) s2 = seg->d ;
		    if (seg->b > 0 && seg->b > s2) s2 = seg->b ;

		    if (s1 < 1000 && s2 > -1000 && s2 - s1 > dSignal)
		      dSignal = s2 - s1 ;
		  }
		}
	      else if (pass == 2)
		{
		  if (seg->a < 0 &&
		      seg->c < 0 &&
		      seg->d < 0 &&
		      seg->b < 0)
		    continue ;
		  break ;
		}		  
	      break ;
	    case MrnaExon: 
	      if (pass)
		continue ;
	      if (iymax >= -1) /* first mrna after probes */
		{
		  if (iymax >= 0) y += iymax ;
		  iymax = -2 ;
		  y00 = y ;
		}
	      y = y00 + dy * seg->bumpy ;
	      if (seg->mrna == oldMrna && !inProbes)
		{
		  graphColor (MAGENTA) ;
		  if (u1 > oldu2 + .1)
		    {
		      int is ;
		      SEG *ss = 0, *ssok = 0 ;
		      char *cp = 0 ;
		      int oldColor = -1, bColor = -1 ;	 
		      
		      if (oldSeg && look->map->introns)
			{
			  for (is = 0, ss = arrp (look->map->introns, 0, SEG) ; is < arrayMax (look->map->introns) ; is++, ss++) 
			    if (ss->type == MrnaIntron && ss->a2 == seg->a1 - 1 && ss->a1 == oldSeg->a2 + 1)
			      { ssok = ss ; break ; }
			}
		      
		      array (look->map->boxIndex, box = graphBoxStart(), int) = iSeg ;
		      
		      nClo = ssok ? ssok->nClo : 1 ;
		      dyClo = 1.0 ;
		      if (look->maxIntronClone > 1)
			{
			  double v, u = log ((double) .8 * look->maxIntronClone) ;
			  u = v = exp (u/8.0) ;
			  dyClo = .1 ;
			  while (nClo > v) { dyClo += .14 ; v = u*v ; }
			}
		      if (ssok && ssok->subType != GAP)
			{
			  {
			    int z;
			    void *vp ;
			    
			    z =  ((int)(171 * oldu2) % 30269) ^ (((int)(172 * u1)) % 30307) ; /* some large int */
			    if (z == 0) z = 1 ; /* unlucky ! */
			    if (graphAssFind (assVoid (z), &vp))
			      bColor = assInt (vp) ;
			    else
			      {
				bColor = intronCol [seg->bumpy ? (seg->bumpy % (numIntronCol - 1) + 1) : 0 ] ;
				graphAssociate (assVoid(z), assVoid(bColor)) ;
			      }
			  }
			  
			  
			  if (ssok->flag)
			    oldColor = graphColor (BLUE) ;
			  if (ssok->subType == MrnaIntronFuzzy)
			    graphLine (oldu2, y + .4, u1, y + .4) ;
			  else
			    {
			      graphLine (oldu2, y + ddy, (oldu2 + u1)/2, y +.2 - dyClo) ;
			      graphLine ((oldu2 + u1)/2, y +.2 - dyClo, u1, y + ddy) ;
			    }
			  if (hasTagCount)
			    {
			      float dd = 1.5, d1=0 ;
			      int n = ssok->nClo + ssok->nTagsA + ssok->nTagsC +  ssok->nTagsD + ssok->nTagsB + ssok->nTagsBlood + ssok->nTagsNB + ssok->nTagsO ;
			      n = ssok->nClo ;
			      for (iTags = 0 ; iTags < NTAGS - NTAGS_OUTSIDER ; iTags++)
				n += ssok->nTags[iTags] ;

			      {
				char *tagText = messprintf("%d", n) ;
				float *tHp, tH[] = {.8,.7,.6,0,.5,.4,0.0} ;
				float lTT = strlen (tagText) ;
				BOOL ok = FALSE ;

				for (tHp = tH ; !ok && *tHp ; tHp++)
				  if (u1 - oldu2 > *tHp * lTT + dd)
				    {
				      float oldTH = graphTextHeight (*tHp) ; 
				      int col = graphColor (BLACK) ;
				      graphText (tagText, (u1+oldu2)/2 - *tHp * (lTT+d1)/2.0, y) ;
				      graphTextHeight (oldTH) ; 
				      graphColor (col) ;
				      ok = TRUE ;
				    }

				tagText = messprintf(".%dk", n/100) ;
				lTT = strlen (tagText) ;
				if (n >= 100 && n < 1000)
				  for (tHp = tH ; !ok && *tHp ; tHp++)
				    if (u1 - oldu2 > *tHp * lTT + dd)
				      {
					float oldTH = graphTextHeight (*tHp) ; 
					int col = graphColor (BLACK) ;
					graphText (tagText, (u1+oldu2)/2 - *tHp * (lTT+d1)/2.0, y) ;
					graphTextHeight (oldTH) ; 
					graphColor (col) ;
					ok = TRUE ;
				      }

				tagText = messprintf("%dk", n/1000) ;
				lTT = strlen (tagText) ;
				if (n >= 1000)
				  for (tHp = tH ; !ok && *tHp ; tHp++)
				    if (u1 - oldu2 > *tHp * lTT + dd)
				      {
					float oldTH = graphTextHeight (*tHp) ; 
					int col = graphColor (BLACK) ;
					graphText (tagText, (u1+oldu2)/2 - *tHp * (lTT+d1)/2.0, y) ;
					graphTextHeight (oldTH) ; 
					graphColor (col) ;
					ok = TRUE ;
				      }

				tagText = messprintf(".%dM", n/100000) ;
				lTT = strlen (tagText) ;
				if (n >= 100000 && n < 10000000)
				  for (tHp = tH ; !ok && *tHp ; tHp++)
				    if (u1 - oldu2 > *tHp * lTT + dd)
				      {
					float oldTH = graphTextHeight (*tHp) ; 
					int col = graphColor (BLACK) ;
					graphText (tagText, (u1+oldu2)/2 - *tHp * (lTT+d1)/2.0, y) ;
					graphTextHeight (oldTH) ; 
					graphColor (col) ;
					ok = TRUE ;
				      }

				tagText = messprintf("%dM", n/1000000) ;
				lTT = strlen (tagText) ;
				if (n >= 1000000)
				  for (tHp = tH ; !ok && *tHp ; tHp++)
				    if (u1 - oldu2 > *tHp * lTT + dd)
				      {
					float oldTH = graphTextHeight (*tHp) ; 
					int col = graphColor (BLACK) ;
					graphText (tagText, (u1+oldu2)/2 - *tHp * (lTT+d1)/2.0, y) ;
					graphTextHeight (oldTH) ; 
					graphColor (col) ;
					ok = TRUE ;
				      }
			      }
			    }
			  if (oldColor > -1)
			    graphColor (oldColor) ;
		
			  vtxtClear (vtxt) ; 
			  vtxtPrintf (vtxt, "%d bp %s %s"
				      , ssok->ln
				      , ssok->txtType
				      , ssok->subType == MrnaIntronFuzzy ? "fuzzy intron" : "intron"
				      ) ;
			  {
			    int n = 0, n1 = 0 ;
			    for (iTags = 0 ; iTags < NTAGS  - NTAGS_OUTSIDER  ; iTags++)
			      n += ssok->nTags[iTags] ;
			    for (; iTags < NTAGS ; iTags++)
				n1 += ssok->nTags[iTags] ;
			    if (ssok->nClo ==  ssok->nCloNM && n + n1 == 0)
			      vtxtPrintf (vtxt, "\n only supported by %s RefSeq model%s"
					  , isOne(ssok->nClo), _multi(ssok->nClo)
					  ) ;
			    else
			      {
				vtxtPrintf (vtxt, "\n") ;
				if (0)
				  {
				    if (ssok->nClo - ssok->nCloNM > 0)
				      vtxtPrintf (vtxt, " %d accession%s\n"
						  , ssok->nClo - ssok->nCloNM, _multi(ssok->nClo - ssok->nCloNM) 
						  ) ;
				    if ( ssok->nCloNM > 0)
				      vtxtPrintf (vtxt, " %d RefSeq%s\n"
						  , ssok->nCloNM, _multi(ssok->nCloNM) 
						  ) ;
				  }
				else
				  {
				    vtxtPrintf (vtxt, " %d GenBank accession%s\n"
						, ssok->nClo, _multi(ssok->nClo)
						) ;
				  }
			      }
			    
			    if (ssok->tissue) vtxtPrintf (vtxt, "%s\n", dictName (look->tissues, ssok->tissue)) ;
			    if (n > 0)
			      vtxtPrintf (vtxt, "%d RNA-seq supporting reads\n", n) ;
			    for (iTags = 0 ; iTags < NTAGS - NTAGS_OUTSIDER; iTags++)
			      if (ssok->nTags[iTags])
				vtxtPrintf (vtxt, "%s %s\n"
					    ,  uTab (8, ssok->nTags[iTags])
					    , titleTags[iTags]
					    ) ;
			    for ( ; iTags < NTAGS; iTags++)
			      if (ssok->nTags[iTags])
				vtxtPrintf (vtxt, "(also %d %s)\n"
					    , ssok->nTags[iTags]
					    , titleTags[iTags]
					    ) ;

			  vtxtPrintf (vtxt, "\n") ;
			  graphBubbleInfo (box, ssok->mrna, className(ssok->mrna), vtxtPtr (vtxt)) ;
			  }
			}
		      else if (ssok)
			{
			  int oldColor = graphColor (LIGHTGRAY) ;
			  graphLine (oldu2, y + .4, u1, y + .4) ;
			  graphColor (oldColor) ; 
			  cp = messprintf ("Sequence gap %d bp", ssok->ln) ;
			}
		      keySet (boundaries, iBound++) = (seg-1)->x2 ;
		      keySet (boundaries, iBound++) = seg->x1 ;
		      graphBoxEnd () ;
		      if (bColor != -1)
			graphBoxDraw (box, -1, bColor) ;
		      if (cp)
			graphBubbleInfo (box, seg->key, "mRNA", cp) ;
		    }
		  else if (seg->b1 > (seg - 1)->b2 + 1)
		    {
		      int boxi = graphBoxStart () ;
		      graphLine (oldu2 - .6, y - .6, oldu2 + .6, y - .6) ;
		      graphLine (oldu2 - .6, y - .6, oldu2, y) ;
		      graphLine (oldu2 + .6, y - .6, oldu2, y) ; 
		      graphBoxEnd () ;
		      graphBubbleInfo (boxi, seg->key, "mRNA"
				       , messprintf ("%d bp possibly an intron"
						     , seg->ln 
						     )
			       ) ;
		    }
		  graphColor (BLACK) ;
		}
	      oldMrna = seg->mrna ; 
	      seg->col = MAGENTA ;

	      if (seg->subType == Utr5 || seg->subType == Utr3)
		{
		  boxHeight = .25 ; 
		  seg->bgcol = WHITE ;
		}
	      else if (seg->subType == VeryGoodCodingExon)
		{ 
		  boxHeight = .4 ;
		  seg->bgcol = PALEMAGENTA ;
		}
	      else if (seg->subType == GoodCodingExon)
		{ 
		  boxHeight = .4 ;
		  seg->bgcol = PALEVIOLET ;
		}
	      else if (seg->subType == uORFExon)
		{ 
		  boxHeight = .4 ;
		  seg->bgcol = PALEGREEN ;
		}
	      else
		{
		  boxHeight = .4 ;
		  seg->bgcol = PALEYELLOW ;
		}
	      if (!inProbes && (seg+1)->mrna != seg->mrna)
		{
		  AC_HANDLE h1 = ac_new_handle () ;
		  KEYSET ks ;
		  const char *ccp ;
		  int ln,  myTissue = 0 ;
		  AC_OBJ Mrna = ac_get_obj (look->db, "mrna", name(seg->mrna), h1) ;
		  AC_KEYSET clones = ac_objquery_keyset (Mrna, ">cdna_clone", h1) ;

		  vtxtClear (vtxt) ;
		  if ((ccp = ficheNewGeneExpressionTissue (vtxt, 0, clones, 3, 3, 0, ac_keyset_count (clones))))
		    {
		      if (!look->tissues)
			look->tissues = dictHandleCreate (64, look->h) ;
		      dictAdd (look->tissues, ccp, &myTissue) ;
		    }
		  
		  vtxtClear (vtxt) ; 
		  if (! topMrna)
		    topMrna = y ;
		  ccp = strstr (name (seg->key), name (look->key)) ;
		  if (ccp)
		    ccp += strlen (name (look->key)) ;
		  {
		    char buf [40], *cq ;
		    if (ccp && *ccp++ == '.') ;
		    else ccp = name (seg->key) ;
		    strncpy (buf, ccp, 39) ;
		    cq = strstr (buf, "-unspliced") ;
		    if (cq) *(cq+2) = 0 ;
		    /* vtxtPrint (vtxt, buf) ; */
		    vtxtPrint (vtxt, 1+ gtMrnaSuffix (name (look->key), name(seg->mrna), h)) ;
		  }

		  box = graphBoxStart () ;
		  {
		    float oldTH = graphTextHeight (1.6) ;
		    graphText (vtxtPtr (vtxt), u2+.3, y - .0) ;
		    ln = vtxtMark (vtxt) ;
		    if (ac_has_tag (Mrna,  "Well_supported"))
		      graphLine (u2+.3, y + 1.0, u2+1.3, y + 1.0) ; 
		    graphTextHeight (oldTH) ;
		  }

		  {
		    vTXT myRefSeq = 0 ;
		    
		    ks = queryKey (seg->mrna, ">cdna_clone ; >Read Ref_Seq && IS NM_*") ;
		    if (keySetMax(ks))
		      graphText ("[NM]", u2 + ln + 1, y - .1) ;
		    else
		      {
			keySetDestroy (ks) ;
			ks = queryKey (seg->mrna, ">cdna_clone ; >Read Ref_Seq && IS NM_*") ;
			if (keySetMax(ks))
			  graphText ("[NR]", u2 + ln + 1, y - .1) ;
		      }
		    graphBoxEnd () ;
		    if (ks && keySetMax (ks))
		      {
			int i ;

			myRefSeq = vtxtHandleCreate (h1) ;
			for (i = 0 ; i < keySetMax(ks) ; i++)
			  vtxtPrintf (myRefSeq, ", %s", name(keySet(ks, i))) ;
		      }
		    graphBubbleInfo (box, seg->mrna, className(seg->mrna)
				     , messprintf("%d accession%s%s%s"
						  , ac_keyset_count (clones) 
						  , _multi(ac_keyset_count (clones))
						  , myRefSeq ? vtxtPtr (myRefSeq) : ""
						  , myTissue ? dictName (look->tissues, myTissue) : ""
						  ) 
				     ) ;
		  }
		  keySetDestroy (ks) ;
		  ac_free (h1) ;
		} 
	      if (inProbes)
		{
		  y += 1.3 ; 
		  graphText (name (seg->key), 1, y) ;
		}
	      inProbes = FALSE ; 
	      if (iSeg>0 && 
		  (seg->mrna != (seg - 1)->mrna ||
		   (seg->b1 > (seg - 1)->b2 + 1 && seg->x1 > (seg - 1)->x2 + 1)
		   )
		  )
		keySet (boundaries, iBound++) = seg->a1 ;
	      if (iSeg < arrayMax (look->map->segs) &&
		  (seg->mrna != (seg + 1)->mrna || seg->x2 + 1 < (seg + 1)->x1 )
		  /* || seg->subType != (seg+1)->subType */
		  )
		keySet (boundaries, iBound++) = seg->a2 ;
	      oldu2 = u2 ; olda2 = seg->a2 ; oldSeg = seg ;
	      if (iSeg < arrayMax (look->map->segs) &&
		  ( seg->mrna != (seg + 1)->mrna ||
		    (seg + 1)->type != MrnaExon
		    )
		  )
		hseqDrawValidEnds (look, seg, iSeg, y+.15) ;
	      break ;
    case GAP:
	      continue ;
	    default: continue ;
	      seg->col = BLACK ;
	      seg->bgcol = BLUE ;
	      break ;
	    }

	  switch (pass)
	    {
	    case 0: /* mrna */
	      array (look->map->boxIndex, box = graphBoxStart(), int) = iSeg ;
	      graphRectangle (u1, y + iy + .4 - boxHeight , u2, y + iy + .4 + boxHeight) ;
	      graphBoxEnd () ;
	      graphBoxDraw (box, seg->col, seg->bgcol) ;
	      graphBubbleInfo (box, seg->key, className(seg->key)
			       , messprintf ("%d bp %s"  /* , %d accessions */
					     , seg->ln 
					     , seg->subType == uORFExon ? "uORF " : "exon"
					     /*  , seg->nClo */
					     ) 
			       ) ;
	      break ;
	    case 1: /* probe widgets */
	      array (look->map->boxIndex, box = graphBoxStart(), int) = iSeg ;
	      graphRectangle (u1, y + iy + .4 - boxHeight , u2, y + iy + .4 + boxHeight) ;
	      graphBoxEnd () ;
	      graphBoxDraw (box, seg->col, seg->bgcol) ;
	      seg->widgetBox = box ;
	      graphBubbleInfo (box, seg->key, "Probe"
			       , messprintf ("%s toto", name(seg->key)) 
			       ) ;
	      break ;
	    case 2:
	      /*
		box = hseqDrawAcdbSeg (look, seg, bottomProbe + 18, dSignal, iNN++, &acdbHeight) ;
		array (look->map->boxIndex, box, int) = iSeg ;
		seg->tableBox = box ;
	      */
	      break ;
	    }
	}
    }


 look->map->yCurrent += y ;

 look->map->starBuffer [0] = 0 ;
 look->map->starBox = graphBoxStart () ;
 if (bottomProbe > 0)
   graphTextPtr (look->map->starBuffer, 1, bottomProbe + 1, look->map->graphWidth) ; 
 graphBoxEnd () ;


 if (iBound && topMrna < bottomProbe)
   {
     keySetSort (boundaries) ;
     keySetCompress (boundaries) ;
     graphColor (LIGHTBLUE) ;
     for (iBound = olda2 = 0 ; iBound < keySetMax (boundaries) ; iBound++)
       {
	 ii = keySet (boundaries, iBound) ;
	 if (ii <= olda2 +1)
	   continue ;
	 olda2 = ii ;
	 x = HMAP2GRAPH (look->map, ii) ;
	 graphLine (x, topMrna, x, bottomProbe) ;
       }
     graphColor (GRAY) ;
     for (iSeg = 0, seg = arrp (look->map->segs, 0, SEG); iSeg < arrayMax(look->map->segs) ; iSeg++, seg++)
       {
	 if (seg->type != GAP)
	   continue ;
	 graphRectangle (seg->a1, topMrna, seg->a2, bottomProbe) ;
       }
     graphColor (BLACK) ;
   }

 if (!MAGGIE && bottomProbe > 0)
   {
     hseqDrawHeatCircles (look, &bottomProbe) ;
     hseqDrawAcdbSegs (look, &bottomProbe, dSignal) ;
   }

 if (MAGGIE && bottomProbe > 0)
   {
     bottomProbe += acdbHeight ;
     hseqDrawSignalPlot (look, &bottomProbe) ;
     if (0)
       hseqDrawAcdbSegs (look, &bottomProbe, dSignal) ;
   }

 if (0)
   {
     graphLine (0, bottomProbe + 0.2, look->map->graphWidth, bottomProbe+ 0.2) ;
     graphText ("A signal", 2, bottomProbe + .8) ;
     graphLine (0, bottomProbe + 2.4, look->map->graphWidth, bottomProbe + 2.4) ;
     graphText ("C signal", 2, bottomProbe + 2.8) ;
     graphLine (0, bottomProbe + 4.6, look->map->graphWidth, bottomProbe + 4.6) ;
     graphText ("D signal", 2, bottomProbe + 4.8) ;
     graphLine (0, bottomProbe + 6.8, look->map->graphWidth, bottomProbe + 6.8) ;
     graphText ("B signal", 2, bottomProbe + 6.8) ;
     graphLine (0, bottomProbe + 9.0, look->map->graphWidth, bottomProbe + 9.0) ;
   }

 if (bottomProbe > 0)
   { 
     char **cq = pgm - 1 ;
     int *pc = pgmcol - 1 ;
     
     for (ii = 0, cq = pgm2, pc = pgmcol ; *cq ; ii++, cq++, pc++)
       {
	 box = graphBoxStart () ;
	 graphText (*cq, 30 + (MAGGIE ? 12 : 6) * ii, 3.5) ;
	 graphBoxEnd () ;
	 graphBoxDraw (box, BLACK, *pc) ;
       }
   }
 look->map->yCurrent += 0.1 * dy ;
 if (isGifDisplay) /* adapt the height of the gif/flash screen */
   {
     svgGraphResize (look->map->graphWidth, look->map->yCurrent) ;
     swfGraphResize (look->map->graphWidth, look->map->yCurrent) ;
     graphFitBounds (&look->map->graphWidth,&look->map->graphHeight) ;
   }
 graphColor (ORANGE) ;
 if (0) graphText ("uuu", 1, look->map->graphHeight) ; /* force graphHeight */
 graphColor (BLACK) ;
 bumpDestroy (bumper) ;
 ac_free (h) ;
} /* hseqDrawGene */

/************************************************************/
typedef struct HSEQBPstruct { KEY mrna ; int bumpy ;} HSEQBP ;

static int HSEQBPorder (const void *a, const void *b)  /* for arraySort() call */ 
{
  const HSEQBP *sa = (const HSEQBP*)a, *sb = (const HSEQBP*) b ;

  return  keySetAlphaOrder (&(sa->mrna), &(sb->mrna)) ;
} /* HSEQBPorder */

/************************************************************/
/* bump the mrna and their names */
static void hseqDrawBumpMrnas (HSeq look)
{
  BUMP bumper = bumpCreate (120,  0) ;
  KEY old ;
  SEG *seg ;
  int bumpy, x1, x2, iSeg, jSeg, ii, iiaa = 0, jj ;
  float u1, u2 ;
  Array aa = arrayCreate (64, HSEQBP) ;
  Array bumpy2 = arrayCreate (64, int) ;
  HSEQBP *h ;

  for (iSeg = 0, seg = arrayp (look->map->segs, 0, SEG) ; iSeg < arrayMax (look->map->segs) ; 
       iSeg++, seg++)
    if (seg->type == MrnaExon && seg->mrna)
      {
	x1 = x2 = seg->x1 ;
	old = seg->mrna ;
	for (jSeg = iSeg ; jSeg < arrayMax (look->map->segs) && seg->type == MrnaExon && seg->mrna == old ; jSeg++, seg++)
	  x2 = seg->x2 ;
	seg = arrayp (look->map->segs, iSeg, SEG) ;
	bumpy = 0 ;
	u1 = HMAP2GRAPH (look->map, x1) ;
	u2 = HMAP2GRAPH (look->map, x2) ;
	bumpItem (bumper, 1, (u2 - u1) + 6, &bumpy, &u1) ;
	h = arrayp (aa, iiaa++, HSEQBP) ;
	h->mrna= seg->mrna ; h->bumpy = bumpy + 1 ;
	for (jSeg = iSeg ; jSeg < arrayMax (look->map->segs) && seg->type == MrnaExon && seg->mrna == old ; jSeg++, seg++)
	  {
	    seg->bumpy = bumpy + 1 ;
	  }
	iSeg = jSeg - 1 ; seg-- ;
      }
  arraySort (aa, HSEQBPorder) ;
  for (ii = jj = 0 ; ii < arrayMax (aa) ; ii++)
    {
      h = arrp (aa, ii, HSEQBP) ;
      if (!array (bumpy2, h->bumpy, int))
	array (bumpy2, h->bumpy, int) = ++jj ;
    }
  for (iSeg = 0, seg = arrayp (look->map->segs, 0, SEG) ; iSeg < arrayMax (look->map->segs) ; 
       iSeg++, seg++)
    {
      if (seg->bumpy)
	seg->bumpy = array (bumpy2, seg->bumpy, int) - 1 ;
    }

  bumpDestroy (bumper) ;
  arrayDestroy (aa) ;
  arrayDestroy (bumpy2) ;
  arraySort (look->map->segs, hseqSegOrder) ;
} /*hseqDrawBumpMrnas */

/************************************************************/

static void hseqDrawScale (HSeq look)
{
  float u1, u2, u3, du, y1, y2 ;
  int sBox, i, sc = utMainRoundPart (look->map->max/10) ;
  int sc1 =  utMainPart (sc/5), sc3 ;
  
  sBox = graphBoxStart () ;
  u1 = .7 ;
  du = HMAP2GRAPH (look->map, sc) - HMAP2GRAPH (look->map, 0) ;
  u2 = u1 + du ;
  y2 = look->map->graphHeight - .07 ;
  graphLine (u1, y2, u2, y2) ;
  y1 = y2 - .4 ;
  graphLine (u1, y1, u1, y2) ; graphLine (u2, y1, u2, y2) ;
  y1 = y2 - .2 ;
  du = HMAP2GRAPH (look->map, sc1) - HMAP2GRAPH (look->map, 0) ;
  for (i = 1 ; i * sc1 < sc ; i++)
    { 
      u3 = u1 + i * du ;
      if (2*i*sc1 == sc)
	{
	  graphLine (u3, y2 - .4, u3, y2) ;
	  sc3 = i*sc1 ;
	  while (sc3 >= 1000)
	    sc3 /= 1000 ;
	  graphText (messprintf ("%d", sc3), u3 - .3 , y1 - 1.2) ;
	}
      else
	graphLine (u3, y1, u3, y2) ;
    }
  if (sc >= 1000000)
    graphText (messprintf ("%dMb", sc/1000000), u2 - .5 , y1 - 1.2) ;
  else if (sc >= 1000)
    graphText (messprintf ("%dkb", sc/1000), u2 - .5 , y1 - 1.2) ;
  else
    graphText (messprintf ("%dbp", sc), u2 - .5 , y1 - 1.2) ;
  graphText ("0",  u1 - .4, y1 - 1.2) ;
  graphText (" ",  1, y2 - .6) ;
  graphBoxEnd () ;
  graphBoxDraw (sBox, BLACK, WHITE) ;
} /* hseqDrawScale */

/************************************************************/

static int comaFormat (char *buf1, int x)
{
  int i, j, k ;
  char buf[64] ; 
  
  sprintf(buf1, "%d" , x) ;
  
  for (j = k = 0,  i = strlen(buf1) ; i > 0 ; i--)
    {
      buf[j++] = buf1[i - 1] ;
      if ((++k % 3) == 0)
	buf[j++] = ',' ;
    }
  buf[j] = 0 ;
  if (buf[j-1] == ',') buf[j-1] = 0 ;
  for (j = 0,  i = strlen(buf) ; i > 0 ; i--)
    buf1[j++] = buf[i-1] ;
  buf1[j] = 0 ;

  return j ;
} /* comaFormat */

/************************************************************/

static void hseqDraw (void)
{
  float xx = 1 ;
  float yy = 2 ;
  char *cp ;
  const char *ccp ;
  HSeq look = currentHSeq("hseqDraw") ;

  hseqMenuSelect (look) ;
  graphClear () ;

  graphMenu (arrp(look->menu, 0, MENUOPT)) ;
  cp = messprintf("Gene %s", name(look->key)) ;
  graphText (cp, xx, yy) ;  xx +=  strlen (cp) + 2 ;
  graphText ("5'", xx, yy) ;  xx += 3 ;
  /* draw arrow */
  { 
    float *fp, dy = .2, y = yy + 2.5*dy, u1 = xx, ua = xx + 6, u2 = xx + 8 ;
    Array pp = arrayCreate (16, float) ;
    
    fp = arrayp (pp, 15, float) ; /* make room */
    fp = arrp (pp, 0, float) ;
    *fp++ = u1 ; *fp++ = y - dy ;
    *fp++ = ua ; *fp++ = y - dy ;
    *fp++ = ua ; *fp++ = y - 3 * dy ;
    *fp++ = u2 ; *fp++ = y ;
    *fp++ = ua ; *fp++ = y + 3 * dy ;
    *fp++ = ua ; *fp++ = y + dy ;
    *fp++ = u1 ; *fp++ = y + dy ;
    *fp++ = u1 ; *fp++ = y - dy ;
    graphLineSegs (pp) ;
    arrayDestroy (pp) ;
    xx = u2 + 1 ;
  }
  graphText ("3'", xx, yy) ;  xx += 3 ;
  ccp = name(look->intMap) ;
  if (0 && strstr(ccp, "CHROMOSOME_")) ccp += strlen ("CHROMOSOME_") ;
  
  {
    char buf1[64], buf2[64] ;
    comaFormat (buf1, look->a1) ;
    comaFormat (buf2, look->a2) ;
    graphText (messprintf ("encoded on %s strand of chromosome %s from %s to %s",
			   look->a1 < look->a2 ? "plus" : "minus"
			   , ccp
			   , buf1, buf2
			   )
	       , xx, yy) ;
  }

  look->map->messageBox = graphBoxStart () ;
  graphTextPtr (look->map->messageText, 1, 1, look->map->graphWidth) ;
  graphBoxEnd () ;

  look->map->yCurrent = 3.5 ;

  arrayMax(look->map->boxIndex) = 0 ;
  hseqDrawBumpMrnas (look) ;
  hseqDrawGene (look) ;
  hseqDrawScale (look) ;
  graphRedraw () ;
}

/************************************************************/
/*************************************************************/

BOOL hseqDisplay (KEY key, KEY from, BOOL isOldGraph)
{
  HSeq  oldlook = 0, look ;
  KEY gene = 0, mrna = 0, product = 0 ;
  AC_HANDLE handle = 0 ;
		/* centre on parent object if poss */
  if (!key)
    return FALSE ;

  cDNAAlignInit () ;
  if (!ac_db)
    ac_db = ac_open_db (0,0) ;

  if (
      isOldGraph && 
      graphAssFind (&HSeq_MAGIC, &oldlook) &&
      oldlook->magic == &HSeq_MAGIC &&
      oldlook->key == key &&
      oldlook->map &&
      graphActivate (oldlook->map->graph)
      )
    { 
      oldlook->from = from ;
      hseqDraw () ;
      return TRUE ;
    }
  if (oldlook)
    hseqDestroy () ;

  handle = handleCreate () ;
  look = (HSeq) halloc (sizeof (struct HSeqStruct), 0) ;  
  look->magic = &HSeq_MAGIC;

  look->db = ac_db ;          /* cosmetic, since we are inside xace */
  look->h = handle ;
  look->key = key ;
  look->from = from ;
  look->gene = gene ;
  look->mrna = mrna ;
  look->product = product ;

  look->onlyNm = oldOnlyNm ;
  look->onlyTouched = oldOnlyTouched ;
  look->allProbes = oldAllProbes ;
  look->showNoSignalProbes = oldShowNoSignalProbes ;
  look->showOnlyExactProbes = oldShowOnlyExactProbes ;
  look->showAmbiguousProbes = oldShowAmbiguousProbes ;

  if (MAGGIE == -1)
    {
      KEYSET ks = query (0, "find probe IS Jun_*") ;
      if (keySetMax (ks))
	MAGGIE = 1 ;
      else
	MAGGIE = 0 ;
      keySetDestroy (ks) ;
    }

  look->map = hMapCreate (look, look->h, hseqDraw) ;
  if (!hseqConvert (look) || ! (isOldGraph || displayCreate (DtHSEQ)))
    { 
      handleDestroy (handle) ;
      ac_free (look) ;
      display (key,from,TREE) ;
      return FALSE ;
    }
    
  graphRegister (DESTROY, hseqDestroy) ;

  graphRetitle (name(key)) ;

  graphAssociate (&HSeq_MAGIC, look) ; /* attach look to graph */
  hMapGraphAssociate (look->map) ; /* redraws */

  return TRUE ;
} /* hseqDisplay */

/*************************************************************/
/*************************************************************/
/*************************************************************/
