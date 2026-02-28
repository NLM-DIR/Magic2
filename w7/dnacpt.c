/*  File: dnacpt.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Apr 22 19:06 1999 (rd)
 * * Nov 19 09:15 1998 (edgrif): make char* to TINT_<colours> unsigned.
 * * Jul 16 12:26 1998 (edgrif): Introduced new private fmap_.h, really
 *      I would like to just use the public header fmap.h, but
 *      dnacptDumpSegs accesses the LookStruct of fmap so I'm using the
 *      internal header for now.
 *      This file illustrates why headers/function prototypes should be
 *      used, it calls fMapDumpSegs (look, fil) but fMapDumpSegs only
 *      takes one argument...the look struct...sigh...
 * * Feb 18 11:42 1996 (mieg)
 * * Feb  7 19:23 1993 (rd): changed TotalLength to work on Sequence 
 	keyset,	and set the sequence Length value if WriteAccess
 * * Jan 15 17:50 1992 (mieg): splice site consensus
 * * Oct 24 17:45 1991 (mieg): wobble on first letter
 * Created: Thu Oct 24 17:45:14 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: dnacpt.c,v 1.26 2020/05/30 16:50:33 mieg Exp $ */

#include "acedb.h"
#include "graph.h"
#include "sysclass.h"
#include "systags.h"
#include "tags.h"
#include "classes.h"
#include "a.h"
#include "bs.h"
#include "dna.h"
#include "bitset.h"
#include "peptide.h"
#include "lex.h"
#include "plot.h"
#include "pick.h"
#include "fmap_.h"
#include "session.h"
#include "restriction.h"
#include "query.h"
#include "call.h"
#include "peptide.h"
#include "parse.h"
#include "table.h"
#include "freeout.h"
/* #include "dnasearch.h" */
#include "dnaalign.h"

 /* we may change _VMotif in the future */
#define MYMOTIF  _VMotif
#define BOX_LENGTH 40

typedef struct DNACPTSTUFF
{ void *magic;        /* == &MAGIC*/
  Graph graph, dnacptGraph ;
  int restrictionBox ;
  char restriction [BOX_LENGTH+ 1] ;
  char restriction2 [BOX_LENGTH+ 1] ;
  int curr, line, frame, mmn, mmch ;
  BOOL useKeySet, complement, amino ;
  KEYSET cutOnce, cutTwice, cutThrice,cutMore, cutNone, enz ;
  int menuState, useKeySetButtonBox, aminoButtonBox ;
} *DNACPT ;

static int MAGIC_CPT = 93452 ;
static int MAGIC_CODON_USAGE = 93453 ;
static Graph lastDnacptGraph = 0 ;
static KEYSET  enzSet = 0 ;

#define graphBoxBox(_box) { \
	       float _x1, _y1, _x2, _y2 ; \
	       graphBoxDim (_box, &_x1, &_y1, &_x2, &_y2) ; \
	       graphRectangle (_x1 - 0.4, _y1 - 0.2, _x2 + 0.4, _y2 + 0.25) ; \
		}	   

#include "display.h"	/* must come after DNACPT definition */

#define STATUS_PROTEIN 0x02 /* from fmapdisp.c */

#define DNACPTGET(name)     DNACPT dnacpt ;\
                                       \
                          if (!graphAssFind (&MAGIC_CPT, &dnacpt)) \
		            messcrash ("(%s) can't find graph",name) ; \
                          if (dnacpt->magic != &MAGIC_CPT) \
                            messcrash ("(%s) received a wrong pointer",name)  ;\
                          displayPreserve()  

void dnaAnalyse (void) ;
void dnaCptFindRepeats(void) ;
static void dnacptDisplay(void) ;
static void dnacptPick (int k, double x) ;
static int restrictionCompletion(char *cp, int len) ;
/* static void dnacptSearchFrancois (void) ; */
static void dnacptEstMatch (void)  ;

/* #define CODE_FILE */
/* writes out coding sequences found in "Codon Usage" option to coding.seq */

#ifdef CODE_FILE
static FILE *codeFile = 0, *exonFile = 0 ;
#endif

static TABLE *table = 0 ;
static KEY tabKey = 0 ;

/************************************************************************/
/*************************************************************************/
                 /*dnaRepaint has only white paint now! */
void dnaRepaint(Array colors)
{ char *cp ;
  mysize_t mysize = (mysize_t) sizeof(char) ;

  if (arrayMax (colors))
    { cp = arrp(colors, 0, char) ;
      memset (cp, WHITE, mysize) ;
    }
}

/************************************************************************/

void dnacptClear (void)
{ Array dna, colors ;
  KEY seqKey ;
  void *look ;
  DNACPTGET("dnacptClear") ;

  if (dnacpt && !dnacpt->useKeySet &&
      fMapActive (&dna, &colors, &seqKey, &look))
    dnaRepaint(colors) ; 

  dnacptDisplay () ;
}

/************************************************************************/

static void dnacptLengths (void)
{
  Array dna ;
  register int i, n =0, t = 0 ;
  Array ll = arrayCreate(50,int) ;
  KEYSET kSet ;
  void *look ;
  int length ;
  DNACPTGET("dnacptLengths") ;

  dnacptDisplay() ;

  if (dnacpt->useKeySet)
    { 
      if (!keySetActive (&kSet, &look)) 
	{ messout ("First select a KeySet containing sequences") ;
	return ;
	}
      for (i = 0 ; i < keySetMax(kSet) ; i++)
	if ((dna = dnaGet(keySet(kSet,i))))
	  { length = arrayMax (dna) ;
	  n += length ;
	  t++ ;
	  array (ll,arrayMax(dna)/1000,int)++ ;
	  arrayDestroy(dna) ;
	  }
      graphText (messprintf ("Found %d sequences, total length %d",t,n),
		 4,dnacpt->line++) ;
      graphTextBounds (80,dnacpt->line += 2) ;
      graphRedraw () ;
      plotHisto ("Sequences lengths in kilobase", ll) ;
    }
  else
    messout ("This only works on a keyset of sequences") ;
}

/*************************************************************************/

static void dnacptFingerPrint(void)
{
  Array dna ;  Array colors ; KEY key ;
  void *look ;
  Array hind3, sau3a , fp ; 
  int i, j, t = 0, from, to, origin ;
  char *tmpFil = 0 ;
  FILE *f = filtmpopen (&tmpFil,"w") ;
  DNACPTGET("dnacptFingerPrint") ;

  if(!fMapActive(&dna,&colors,&key,&look))
    { messout("Please, first select a dna window") ;
      return ;
    }
  fMapFindZone (look, &from, &to, &origin) ;

  /* HindIII == AAGCTT */
  /* Sau3AI = GATC */

          /* note that pickMatch returns pos + 1 or zero */

  hind3 = arrayCreate(30,int) ;
  sau3a = arrayCreate(30,int) ;
  fp = arrayCreate(30,int) ;

  dnacptFingerPrintCompute(dna, from, to, colors, hind3, sau3a, fp);

  graphText(messprintf("Found %d hind3 sites",t = arrayMax(hind3)),
	    2, dnacpt->line++) ;

  for(i=0;i<t;i++)
   { graphText(messprintf("(%d)",arr(hind3,i,int)),
	       3, dnacpt->line ++) ;
     if (f) fprintf(f, "(%d)\n",arr(hind3,i,int)) ;
   }
  dnacpt->line += 3 ;

	
  graphText(messprintf("Position and length of the finger printing bands"),
	    2, dnacpt->line++) ;
  if (f)
    fprintf(f, "\n\n\nPosition and length of the finger printing bands\n"),

  dnacpt->line += 3 ;
  graphText(messprintf("%d Ordered Bands",t = arrayMax(fp)),
	    2, dnacpt->line++) ;
  if (f) fprintf(f, "\n\n%d Ordered Bands\n",t = arrayMax(fp)) ;
  for(i=0;i<t;i++)
   { graphText(messprintf("(%d)",arr(fp,i,int)),
	  3, dnacpt->line++) ;
     if (f) fprintf(f, "(%d)\n",arr(fp,i,int)) ;
   }
  
  graphTextBounds (80,dnacpt->line += 2) ;
  graphRedraw() ;
  
  if(fp)
    { i = arrayMax(fp) ;
      if(i)
	{ Array histo = arrayCreate(2000,int) ;
	  while(i--)
	    { j = arr(fp,i,int) ;
	      if (j)
		array(histo,j,int)++ ;
	    }
	  plotHisto
	    (messprintf("%s Finger Print",
			name(key)), 
	     histo) ;
	}
    }

  arrayDestroy(fp) ;
  arrayDestroy(hind3) ;
  arrayDestroy(sau3a) ;    
  if (f) filclose(f) ;
  fMapReDrawDNA (look) ;
}

void dnacptFingerPrintCompute (Array dna, int from, int to,
			       Array colors, Array hind3, Array sau3a,
			       Array fp)
{
  int i , j , n , t = 0 , nfp ,
     h, u, s, ih, is ;
  register unsigned char *cp , *cp0 ;
  BOOL lastIsHind3 = FALSE ;
  static char h3[]={'*',A_,A_,G_,C_,T_,T_,'*',0} ; /* RD initialisation */
  static char s3a[]={'*',G_,A_,T_,C_,'*',0} ; /* RD initialisation */

  /* HindIII == AAGCTT */
  /* Sau3AI = GATC */

          /* note that pickMatch returns pos + 1 or zero */

      /* Find all hind3 sites */
  t = 0 ;
  n = arrayMax(dna) ;
  if (n >= to) n = to-1 ;
  if (colors) dnaRepaint(colors) ;
  cp = cp0 = arrp(dna,from, unsigned char) ;
  while ((i = pickMatch((char *)cp, h3)))
    { 
      if (i>n)
	break ;
      cp += i ; i = cp - cp0 - 1 ;
      array(hind3,t++,int) = i - from ;
      j = 6 ;
      if (colors) while(j--)
	arr(colors,from + i + j,char) |= TINT_LIGHTGREEN ;
    }


      /* Find all sau3a sites */
  t = 0 ;
  cp = cp0 = arrp(dna,from, unsigned char) ;
  while ((i = pickMatch((char *)cp , s3a)))
    { 
      if(i>n)
	break ;
      cp += i ; i = cp - cp0 - 1 ;
      array(sau3a,t++,int) = i - from ;
      j = 4 ;
      if (colors) while(j--)
	arr(colors,from + i + j,char) |= TINT_CYAN ;
    }

	
  /* Now we measure from hind3 towards sau3a sites */
  nfp = 0;

  ih = (long)arrayMax(hind3) - 1 ;
  is = (long)arrayMax(sau3a) - 1 ;
  u  = arrayMax(dna) ;
  
  h = ih>=0 ? array(hind3,ih--,int) : -1 ;
  s = is>=0 ? array(sau3a,is--,int) : 0 ;
  lastIsHind3 = FALSE ;
  while(u > 0 || s > 0 || h > 0)
    if(s > 0 && h < s) /* I am on a sua3a site */
      { if(lastIsHind3) array(fp,nfp++,int) = u - s ;
	
	u = s ;
	s = is>=0 ?  array(sau3a,is--,int) : 0 ;
	lastIsHind3 = FALSE ;
      }
    else /* I am on a hind3 site */
      { array(fp,nfp++,int) = u - h ;

	u = h ;
	h = ih>=0 ? array(hind3,ih--,int) : -1 ;
	lastIsHind3 = TRUE ;
      }

  arraySort(fp,intOrder) ;
}

/***************************/

static void dnacptTranslate (Array dna, Array protein, int from, int length, int frame, char *tt)
{
  char *cp, *pp ; int  n = length/3 ;
  protein = arrayReCreate(protein, n + 3, char) ;

  cp = arrp(dna,from,char)  - 3 ;
  arrayMax(protein) = n ;
  pp = arrp(protein,0,char) - 1 ;
  
  while(cp +=3, pp++, n-- )
    *pp = frame >=0 ? e_codon(cp, tt) : e_reverseCodon(cp, tt) ;
}

/*************************************************************/

static void dnacptColorMatches (Array sites, int n, int col,
			        Array colors, int len, int gap)
{ int  t, i, j, c ;
  Site *s ;
  unsigned char mycolor[5] , *cp = mycolor ;

  *cp++ = TINT_RED ;
  *cp++ =  TINT_LIGHTGRAY ;
  *cp++ =  TINT_MAGENTA ;
  *cp++ =  TINT_CYAN ;
  *cp++ =  TINT_LIGHTGREEN ;

  for (t = n, s = arrp(sites, t, Site) ; t < arrayMax(sites) ; s++, t++)
				/* not the extremities */
    { if (!s->type) continue ;
      i = s->i ;
      if (i < 0 || i + gap*(len - 1) >= arrayMax(colors))
	{ /* messerror ("dnacptColorMatches over flow") ;*/ continue ;}
      j = len ; 
      switch (s->type)
	{
	case 1: c = TINT_LIGHTGREEN ; break ; /* for up and down oligos */
	case 2: c = TINT_CYAN ; break ;
	default: c = mycolor[col % 5]; break ;
	}
      if (gap > 1)
	c |= 128 ;
      while (j--)
	arr(colors, i + gap*j, char) |= c ;
    }
}

/*****************************************************************/
static int END_MARK = 0 ;
/* use frame = -1 reversed template, 
                0 direct = original template, 
		1 palindrome template 
*/
static int dnacptMatch (Array dna, int frame,
			Array sites, char *template, 
			int from, int length, int mark, 
			int maxError, int maxn)
{
  register int  n = arrayMax(dna), t = 0, i ;
  register unsigned char *cp, *cp0 ;
  int templateLength = strlen(template) ;
  char matchString[256] ;

  if (templateLength > 255)
    { messerror ("templateLength %d > 255 in dnacptMatch, can t proceed",
		 templateLength) ;
    return 0 ;
    }

  matchString[templateLength] = 0 ;
  
  if (from < 0)
    from = 0 ;
  if (length < 0)
    length = 0 ;
  if (n > from + length)
    n = from + length ;
  if (from > n)
    from = n ;

  cp = arrp(dna,from, unsigned char) ;
  cp0 = arrp(dna,0, unsigned char) ;

           /* note that dnacptMatch returns pos + 1 or zero */
  t = arrayMax(sites)  ;
  if (!t)
    { arrayp(sites,t,Site)->i = from ;
      arrayp(sites,t++,Site)->mark = 0 ;
    }
  else
    t-- ; /* cancel the previous END */

  while ((i = dnaPickMatch (cp, length, (unsigned char *)template, maxError, maxn)))
    { 
      if (i > n)
	break ;
      cp += i ; length -= i ;
      arrayp(sites,t,Site)->i = cp - cp0 - 1 ;
      arrayp(sites,t,Site)->type = frame + 2 ; 
			/* bit 1 = -ve frame, bit 2 = +ve frame */
      arrayp(sites,t++,Site)->mark = mark ;

      if (table)		/* put matches in the table */
	{ int j = tabMax(table,0) ;
	  int k ;

	  tableKey(table,j,0) = tabKey ;
	  if (frame == -1)
	    { tableInt(table,j,1) = cp - cp0 + templateLength - 1 ;
	      tableInt(table,j,2) = cp - cp0 ;
	      for (k = 0 ; k < templateLength ; ++k)
		matchString[k] = dnaDecodeChar[(int)complementBase(*(cp+templateLength-k-2))] ;
	      tableSetString(table,j,3,matchString) ;
	    }
	  else
	    { tableInt(table,j,1) = cp - cp0 ;
	      tableInt(table,j,2) = cp - cp0 + templateLength - 1 ;
	      for (k = 0 ; k < templateLength ; ++k)
		matchString[k] = dnaDecodeChar[(int)*(cp+k-1)] ;
	      tableSetString(table,j,3,matchString) ;
	    }
	}
    }

  arrayp(sites,t,Site)->i = n ; /* last fragment */
  arrayp(sites,t++,Site)->mark = END_MARK ;

/*
  printf ("dnacptMatch found %d sites for frame %d template %s\n",
          t - old, frame, dnaDecodeString(template)) ;
*/

  return arrayMax(sites) ;
}

/***************************************************************/

static void dnacptAminoMatch (Array dna, int frame, 
			      Array sites, char *template, 
			      int from, int length, int mark, char *tt)
{ register int  t = 0 , i ;
  register char *cp, *cp0, *cr ;
  static Stack s = 0 ;
  static Array protein = 0 ;

  from = from + (frame + 6)%3 ;
  if (from + length > arrayMax(dna))
    length = arrayMax(dna) - from ;
  if (length <= 0)
    return ;
  length = 3*(length/3) ;
  if (length < 3)
    return ;
  protein = arrayReCreate(protein,length,char) ;
  dnacptTranslate(dna, protein, from, length, frame, tt) ;

  cp = cp0 = arrp(protein, 0 , char) ;
  array(protein,length/3,char) = 0 ;  /* neccesary in pickMatch */

          /* note that pickMatch returns pos + 1 or zero */
  t = arrayMax(sites)  ;
  if (!t)
    { arrayp(sites,t,Site)->i = from ;
      arrayp(sites,t++,Site)->mark = 0 ;
    }
  else
    t-- ; /* cancel the previous END */

  s = stackReCreate(s, 20) ;
  pushText(s, "*") ;
  catText(s, template) ;
  catText(s, "*") ;
  if (frame < 0)  /* reverse */
    { cr = stackText(s, strlen(template)) ;
      while ((*cr-- = *template++)) ;
      *++cr = '*' ;
    }
  template = stackText(s, 0) ;

  while ((i = pickMatch(cp, template)))
    { 
      if (i>length/3)
	break ;
      cp += i ; 
      arrayp(sites,t,Site)->i = 3 * (cp - cp0 - 1) + from ;
      arrayp(sites,t,Site)->type = frame >= 0 ? 2 : 1 ;
      arrayp(sites,t++,Site)->mark = mark ;
    }
  arrayp(sites,t,Site)->i = from + length ; /* last fragment */
  arrayp(sites,t++,Site)->mark = END_MARK ;
  if (frame == 2)
    arrayDestroy(protein) ;
}

/****************************************************************/

void dnacptMultipleMatch (Array sites, 
			  Array dna, Array protein, int frame,
			  Array colors, Stack s, Stack sname, 
			  int from, int length, 
			  BOOL amino, int mmch, int mmn, char *tt)
{
  int n, col = 5; /* depends on defn of TINTS, start with magenta and go up */
  char *cp, *cq, *cr, *cs; 
  int mark ;
  
  /*
  int beg = 0 ;
    if (amino)
      beg = (frame + 6) % 3 ;
  */
  stackCursor(s, 0) ; stackCursor(sname, END_MARK) ;
 
  while ((cp = stackNextText(s)))
    { stackNextText(sname) ;
      mark = stackPos(sname) ;
      n = arrayMax(sites) ;
      if (amino)
	{
	  dnacptAminoMatch (dna, frame, sites, cp, 
			    from, length, mark, tt) ;
	  if (colors && arrayMax(sites) > n )
	    dnacptColorMatches (sites, n ? n-1 : 0, frame>=0 ? 6:7,
				colors, strlen(cp), 3);
	}
      else
	{  /* reverse the template */
	  cr = messalloc(strlen(cp)+ 1) ;
	  cq = cr ; cs = cp ;
	  while(*cs) cs++ ;
	  while(--cs >= cp) *cq++ = *cs ;
	  *cq = 0 ; cq = cr ;
	  while(*cq = complementBase(*cq), *++cq) ; /* complement */
	  if (!strcmp(cp, cr))
	    { if (frame == 0) /* flag as palindrome on +ve search */ 
		dnacptMatch (dna, 1, sites, cp, 
			     from, length, mark, mmch, mmn) ; 
	                      /* do not search in negative frame */
	      if (colors && arrayMax(sites) > n)
		dnacptColorMatches (sites, n ? n - 1 : 0 , col--, colors,
				    strlen(cp), 1) ;
	    }
	  else
	    { dnacptMatch(dna, frame, sites, frame == -1 ? cr : cp, 
			  from, length, mark, mmch, mmn) ; 
	      if (colors && arrayMax(colors) > 10 && arrayMax(sites) >n )
		dnacptColorMatches (sites, n ? n-1 : 0, frame == 0 ? 6:7,
				    colors, strlen(cp), 1) ;
	    }
	  messfree(cr) ;
	}
      
      if (colors && arrayMax(sites) > n ) 
	dnacptColorMatches (sites, n ? n - 1 : 0 , ++col,  colors, 
			    strlen(cp), amino ? 3 : 1 ) ;
    }
}

/*************************************************************************/

static void dnacptXnone (void)
{
  KEYSET ks ;
  DNACPTGET ("dnacptNone") ;

  ks = dnacpt->cutNone ;
  if (keySetExists(ks) && keySetMax(ks))
    keySetNewDisplay(keySetCopy(ks),"Enzymes cutting none") ;
}

static void dnacptXonce (void)
{
  KEYSET ks ;
  DNACPTGET ("dnacptOnce") ;

  ks = dnacpt->cutOnce ;
  if (keySetExists(ks) && keySetMax(ks))
    keySetNewDisplay(keySetCopy(ks),"Enzymes cutting once") ;
}

static void dnacptXtwice (void)
{
  KEYSET ks ;
  DNACPTGET ("dnacptTwice") ;

  ks = dnacpt->cutTwice ;
  if (keySetExists(ks) && keySetMax(ks))
    keySetNewDisplay(keySetCopy(ks),"Enzymes cutting twice") ;
}

static void dnacptXthrice (void)
{
  KEYSET ks ;
  DNACPTGET ("dnacptThrice") ;

  ks = dnacpt->cutThrice ;
  if (keySetExists(ks) && keySetMax(ks))
    keySetNewDisplay(keySetCopy(ks),"Enzymes cutting thrice") ;
}

static void dnacptXmore (void)
{
  KEYSET ks ;
  DNACPTGET ("dnacptMore") ;

  ks = dnacpt->cutMore ;
  if (keySetExists(ks) && keySetMax(ks))
    keySetNewDisplay(keySetCopy(ks),"Enzymes cutting more") ;
}

static MENUOPT dnacptXbuttons[] =  /* others */
  {
 { dnacptXnone, "Nowhere"},
 { dnacptXonce, "Once"},
 { dnacptXtwice, "Twice"},
 { dnacptXthrice, "Thrice"},
 { dnacptXmore, "Often"},
 { 0, 0 }
   } ;

/*************************************************************************/

int siteOrder (const void *a, const void *b)
{ return ((const Site*)a)->i - ((const Site*) b)->i ;
}

static void dnacptShowMatches(KEY seqKey, Array sites, Stack sname, int origin)
{ register int  t, i, box, nrow, n ;
  Site *site ;
  Array aa ; char *cp, **cpp ; float dl ;
  KEY key ; Stack s2 = 0 ;
  DNACPTGET("dnacptShowMatches") ;

  arraySort(sites, siteOrder) ;
  t = arrayMax(sites) ;
  
  if (dnacpt->enz) 
    { graphText ("Enzymes cutting:", 1,dnacpt->line - .6) ; dnacpt->line += 1.8 ;
      box = graphBoxStart() ;
      graphButtons (dnacptXbuttons, 1, dnacpt->line, 40) ; 
      graphBoxEnd() ;
      graphBoxDim(box,0,0,0,&dl) ; 
      dnacpt->line = dl + 2.2 ;
    }
  graphText(messprintf("Found %d sites %s in %s",t - 2, 
		       dnacpt->restriction, name(seqKey)), 
	    1,dnacpt->line++) ;
  if (t<=2)
    return ;

  dnacpt->line++ ;
  nrow = (t+1)/2 ;
  i = arrayMax(sites) ; site = arrp(sites,0,Site)  - 1 ; 
  while(site++, i--) 
    if (strlen(stackText(sname,site->mark)) > 9)
      { nrow = t+1 ; break ; }  /* just one col for display */
  site = arrp(sites,0,Site) ;
  graphText (messprintf("         %5d %s",site->i  - origin + 1, stackText(sname,site->mark)),
	     1, dnacpt->line) ;
  for(i=1, site = arrp(sites,1,Site) ;i<t;i++, site++)
    graphText (messprintf ("+ %4d = %5d %s ",
			   site->i - (site-1)->i ,
			   site->i - origin + 1 , /* No ZERO in biology */
			   stackText(sname,site->mark)),
	       (i >= nrow) ? 27 : 1, dnacpt->line + i%nrow) ;
  dnacpt->line += nrow + 1 ;

  aa = arrayCreate (200, int) ;
  t = arrayMax(sites) ;
  site = arrp(sites,0,Site) - 1 ;
  while (site++, t--)
    array(aa, site->mark,int)++ ;
 
  graphFitBounds (&mapGraphWidth,&mapGraphHeight) ;
  s2 = stackCreate (300) ;

  stackCursor (sname, 0) ;
  dnacpt->cutNone = keySetReCreate(dnacpt->cutNone) ; n = 0 ;
  while (i = stackPos(sname), (cp=stackNextText(sname)))
    if (!array(aa, i, int))
      {
	if (lexword2key(cp, &key, MYMOTIF))
	  keySet(dnacpt->cutNone, n++) = key ;
	catText (s2, stackText(sname,i)) ;
	catText (s2, "  ") ;
      }
  if (stackMark(s2))
    {
      graphText ("Sites not found", 1,  dnacpt->line++) ;
      for (cpp = uBrokenLines (stackText(s2,0), mapGraphWidth - 2) ; *cpp ; ++cpp)
	graphText(*cpp, 2, dnacpt->line++) ;
    }

  stackCursor (sname, 0) ;
  stackClear (s2) ;
  dnacpt->cutOnce = keySetReCreate(dnacpt->cutOnce) ; n = 0 ;
  while (i = stackPos(sname), (cp=stackNextText(sname)))
    if (array(aa, i, int) == 1 &&
	strcmp(cp,"Origin") && 	strcmp(cp,"End"))
      { 
	if (lexword2key(cp, &key, MYMOTIF))
	  keySet(dnacpt->cutOnce, n++) = key ;
	catText (s2, stackText(sname,i)) ;
	catText (s2, ": ") ;
	t = arrayMax(sites) ;
	site = arrp(sites,0,Site) - 1 ;
	while (site++, t--)
	  if (i == site->mark)
	    { catText(s2, messprintf(" %d",  site->i - origin + 1)) ; break ; }
	catText (s2, "  ") ;
      }
  if (stackMark(s2))
    {
      graphText ("Sites found once", 1, dnacpt->line++) ;
      for (cpp = uBrokenLines (stackText(s2,0), mapGraphWidth - 2) ; *cpp ; ++cpp)
	graphText(*cpp, 2, dnacpt->line++) ;
    }

  stackCursor (sname, 0) ; 
  stackClear (s2) ;
  dnacpt->cutTwice = keySetReCreate(dnacpt->cutTwice) ; n = 0 ;
  while (i = stackPos(sname), (cp=stackNextText(sname)))
    if (array(aa, i, int) == 2)
      { 
	if (lexword2key(cp, &key, MYMOTIF))
	  keySet(dnacpt->cutTwice, n++) = key ;
	catText (s2, stackText(sname,i)) ;
	catText (s2, ": ") ;
	t = arrayMax(sites) ;
	site = arrp(sites,0,Site) - 1 ;
	while (site++, t--)
	  if (i == site->mark)
	    catText(s2, messprintf(" %d",  site->i - origin + 1)) ;
	catText (s2, "  ") ;
      }
  if (stackMark(s2))
    {
      graphText ("Sites found twice", 1, dnacpt->line++) ;
      for (cpp = uBrokenLines (stackText(s2,0), mapGraphWidth - 2) ; *cpp ; ++cpp)
	graphText(*cpp, 2, dnacpt->line++) ;
    }

  stackCursor (sname, 0) ; 
  stackClear (s2) ;
  dnacpt->cutThrice = keySetReCreate(dnacpt->cutThrice) ; n = 0 ;
  while (i = stackPos(sname), (cp=stackNextText(sname)))
    if (array(aa, i, int) == 3)
      { 
	if (lexword2key(cp, &key, MYMOTIF))
	  keySet(dnacpt->cutThrice, n++) = key ;
	catText (s2, stackText(sname,i)) ;
	catText (s2, ": ") ;
	t = arrayMax(sites) ;
	site = arrp(sites,0,Site) - 1 ;
	while (site++, t--)
	  if (i == site->mark)
	    catText(s2, messprintf(" %d",  site->i - origin + 1)) ;
	catText (s2, "  ") ;
      }
  if (stackMark(s2))
    {
      graphText ("Sites found thrice", 1, dnacpt->line++) ;
      for (cpp = uBrokenLines (stackText(s2,0), mapGraphWidth - 2) ; *cpp ; ++cpp)
	graphText(*cpp, 2, dnacpt->line++) ;
    }

  stackCursor (sname, 0) ;
  stackClear (s2) ;
  dnacpt->cutMore = keySetReCreate(dnacpt->cutMore) ; n = 0 ;
  while (i = stackPos(sname), (cp=stackNextText(sname)))
    if (array(aa, i, int) > 3)
      { 
	if (lexword2key(cp, &key, MYMOTIF))
	  keySet(dnacpt->cutMore, n++) = key ; 
	if (stackMark(s2))
	  catText (s2, "\n") ;
	catText (s2, messprintf("%d * %s: ", array(aa, i, int),stackText(sname,i))) ;
	t = arrayMax(sites) ;
	site = arrp(sites,0,Site) - 1 ;
	while (site++, t--)
	  if (i == site->mark)
	    catText(s2, messprintf(" %d",  site->i - origin + 1)) ;
      }
  if (stackMark(s2))
    {
      graphText ("Sites found more often", 1, dnacpt->line++) ;
      for (cpp = uBrokenLines (stackText(s2,0), mapGraphWidth - 2) ; *cpp ; ++cpp)
	graphText(*cpp, 2, dnacpt->line++) ;
    }
  graphLine (0,dnacpt->line,mapGraphWidth,dnacpt->line) ; dnacpt->line++ ;
  arrayDestroy (aa) ;
  stackDestroy (s2) ;
}

/*************************************************************************/

static KEY matchSeqTag = 0 ;

static void getMatchSeqTag (void)
{
  if (!matchSeqTag)
#ifdef ACEDB4
    lexaddkey ("Match_sequence", &matchSeqTag, 0) ;
#else
    lexaddkey ("Site", &matchSeqTag, 0) ;
#endif
}
 

Stack dnacptAnalyseRestrictionBox(char *text, BOOL amino, Stack *snp)
{
  OBJ obj ;
  KEY key ;
  char *cname ,*cp;
  Stack s = 0 ;
  static Stack s1 = 0 , s2 = 0 ;
  BOOL error = FALSE ;
  int n , i ;
  KEY _VOligo =  pickWord2Class("Oligo") ;

  getMatchSeqTag () ;		/* make sure matchSeqTag is set */
 
  if(!enzSet && (!text || !*text))
    { messout("First type in a restriction site") ;
      return 0 ;
    }
 
  s = stackCreate(80) ; 
  s1  = stackReCreate(s1, 80) ;
  s2  = stackReCreate(s2, 80) ;
  *snp = stackReCreate(*snp, 50) ;
  pushText(*snp,"Origin") ;
  END_MARK  =  stackMark(*snp) ;
  pushText(*snp,"End") ;

  if (!text && enzSet)
    { i = keySetMax(enzSet) ; n = 0 ;
      while (i--)
	{ key = keySet(enzSet,i) ;
	  if (class(key) == MYMOTIF
	      && (obj = bsCreate(key)))
	    { if (bsGetData(obj,matchSeqTag,_Text,&cname))
		{ n++ ; pushText(s, cname) ; pushText(*snp,name(key)) ;}
	      bsDestroy(obj) ;
	    }
	  else if (_VOligo && class(key) == _VOligo
	      && (obj = bsCreate(key)))
	    { if (bsGetData(obj,_Sequence,_Text,&cname))
		{ n++ ; pushText(s, cname) ; pushText(*snp,name(key)) ;}
	      bsDestroy(obj) ;
	    }
	}
    }
  else
    { stackTokeniseTextOn(s1, text, ";") ; /* tokeniser */
      /* Now replace known enzyme names by their sequence */
      /* and produce a clean message */
      while ((cp = stackNextText(s1)))
	{
	  if(!amino 
	     && lexword2key(cp,&key,MYMOTIF)
	     && (obj = bsCreate(key)))
	    { catText(s2, name(key)) ;
	      catText(s2, " ;  ") ;
	      if (bsGetData(obj,matchSeqTag,_Text,&cname))
		{ pushText(s, cname) ;  pushText(*snp,name(key)) ;}
	      bsDestroy(obj) ;
	    }
	  else
	    {
	      pushText(s, cp) ;  pushText(*snp,cp) ;
	      catText(s2, cp) ;
	      catText(s2, " ;  ") ;
	    }
	}
      memset(text, 0, BOX_LENGTH) ;  
      strncpy(text, stackText(s2, 0), 39) ;
    }
 
        /* Verify that we have sequences */
        /* go to lower case */
  stackCursor(s, 0) ;
  while ((cp = stackNextText(s)))
    { 
      --cp ;
      while(*++cp)
	 *cp = ace_lower(*cp) ;
    }
       /* if dna, try to encode */
  stackCursor(s, 0) ;
  if (!amino)
    while ((cp = stackNextText(s)))
      { 
	n = strlen(cp) ;
	dnaEncodeString(cp) ;
	if (n > strlen(cp))
	  error = TRUE ;
      }
  
  if (error)
    { messout 
	("%s %s %s %s %s",
	 "Type names of known motifs, or valid dna sequences,",
	 "separated by semi columns. For example:",
	 "BamHI ; aatgga ; mksw\n",
	 "Use the tab button to autocomplete the enzymedisplaycrenames",
	 "or use one the enzyme list buttons") ;
      stackDestroy (s) ;
    }
  else
    stackCursor (s, 0) ;

  return s ;
}
 
/**********************************/

static void dnacptRestriction (void)
{
  Array dna, colors ;
  KEY  seqKey ; KEYSET  matchingSeqs = 0, ksa = 0 ;
  void *look ;
  Array sites = 0 ;
  Stack s , sname = 0 ; 
  int from = 0, length, to , origin;
  int ii, frame, level;
  Stack s1 = 0 ;
  char *tt ;
  DNACPTGET("dnacptRestriction") ;
  
  if (!graphCheckEditors(dnacpt->graph,0))
    return ;
  if (strcmp(dnacpt->restriction, dnacpt->restriction2))
    keySetDestroy (dnacpt->enz) ;
  /* the zeroing of dnacpt->enz is complicated,
     but this is to get a nice interface for the user and to avoid
     changing geldisp.c
     mieg: may 97 
     */
  if (dnacpt->enz)
    dnacptDisplay() ;
  enzSet = dnacpt->enz ; /* static communication with geldisp !, sorry */
  s = dnacptAnalyseRestrictionBox (dnacpt->enz ? 0 : dnacpt->restriction, 
				   dnacpt->amino, &sname) ;
  enzSet = 0 ;
  strcpy(dnacpt->restriction2,dnacpt->restriction) ;
  if (!s)
    { dnacptDisplay() ;
      return ;
    }

  if (dnacpt->useKeySet)
    { KEY seq ; KEYSET kSet ; void *look ;
      int nn = 0 , nms = 0 ;

      if (!keySetActive(&kSet, &look)) 
	{ messout("First select a KeySet containing sequences") ;
	  return ;
	}
      ksa = keySetAlphaHeap (kSet, keySetMax(kSet)) ; 
      if (!keySetMax(ksa))
	{ messout("First select a KeySet containing sequences") ;
	  keySetDestroy (ksa) ;
	  return ;
	}
      matchingSeqs = keySetCreate() ;		     
      table = tableCreate (64, "kiis") ;
      tableSetColumnName (table, 0, "seq") ;
      tableSetColumnName (table, 1, "start") ;
      tableSetColumnName (table, 2, "end") ;
      tableSetColumnName (table, 3, "match") ;

      for (ii = 0 ; ii < keySetMax(ksa) ; ii++)
	if ((dna = dnaGet(seq=keySet(ksa,ii))))
	  { 
	    nn++ ;
	    tt = pepGetTranslationTable (seq, 0) ;
	    sites = arrayReCreate(sites, 20, Site) ;
	    tabKey = keySet(ksa,ii) ;
	    if (dnacpt->amino)
	      for (frame = -3 ; frame < 3 ; frame++)
		dnacptMultipleMatch (sites, dna, 0, frame, 
				     0, s, sname, 
				     0, arrayMax(dna), 
				     TRUE, dnacpt->mmch, dnacpt->mmn, tt) ;
	    else
	      for (frame = -1 ; frame < 1 ; frame++)
		dnacptMultipleMatch (sites, dna, 0, frame,
				     0, s, sname, 
				     0, arrayMax(dna),
				     FALSE, dnacpt->mmch, dnacpt->mmn, tt) ;
	    if (arrayMax(sites) > 2)
	      { 
		dnacptShowMatches (seq, sites, sname, 1) ;
		keySet(matchingSeqs,nms++) = seq ;
	      } 
	    arrayDestroy(dna) ;
	  }
      keySetDestroy (ksa) ;
      arrayDestroy(sites) ;
      keySetSort(matchingSeqs) ;
      keySetCompress(matchingSeqs) ;
      graphText (messprintf ("I scanned %d sequences", nn),
		 2, dnacpt->line++) ;
      if (dnacpt->amino)
	graphText (messprintf ("in all six frames", nn),
		   2, dnacpt->line++) ;
      { Graph old = graphActive() ;
	displayCreate(DtKeySet) ;
	graphRetitle("Matching sequences") ;
	keySetShow(matchingSeqs, 0) ;
      matchingSeqs = 0 ;
	graphActivate(old) ;
      }
      keySetDestroy (matchingSeqs) ; /* if empty destroy here */

      { FILE *f ;
	if ((f = filqueryopen (0, 0, "match", "w", "File for match list")))
	  { level = freeOutSetFile (f) ;
	    tableOut (table, '\t', 0) ;
	    filclose (f) ;
	    freeOutClose(level) ;
	  }
      }
      s1 = stackCreate(1000) ;
      level = freeOutSetStack (s1) ; 
      tableOut (table, '\t', 0) ;
      freeOutClose(level) ;
      externalFileDisplay (messprintf ("%d matches", table ? tableMax(table) : 0),  0, s1) ;
      stackDestroy (s1) ;
      tableDestroy (table) ;
      table = 0 ;
    }
  else
    { if (!fMapActive (&dna, &colors, &seqKey, &look))
	{ messout ("First select a DNA window "
		   "or the 'KeySet' button") ;
	  return ;
	}
      dnaRepaint(colors) ;
      fMapFindZone(look, &from, &to, &origin) ;
      length = to - from ;
      tt = pepGetTranslationTable (seqKey, 0) ;

      sites = arrayReCreate(sites, 20, Site) ;
      if (dnacpt->amino) 
	for (frame = -3 ; frame < 3 ; frame++)
	  dnacptMultipleMatch (sites, dna, 0, frame, colors, s, sname,
			       from, length, TRUE, dnacpt->mmch, dnacpt->mmn, tt) ;
      else
	for (frame = -1 ; frame < 1 ; frame++)
	  dnacptMultipleMatch (sites, dna, 0, frame, colors, s, sname,
			       from, length, FALSE, dnacpt->mmch, dnacpt->mmn, tt) ;
      dnacptShowMatches(seqKey,sites,sname, origin) ;
	
      { Graph old = graphActive() ;
	fMapRegisterSites (look, sites, sname) ; /* sites and sitenames will be destroyed there */
	mapColSetByName ("Summary bar", TRUE) ;
	fMapDraw(look,0) ;
	sname = 0 ; sites = 0 ;
	graphActivate (old) ;
      }
      arrayDestroy(sites) ;
    }
  stackDestroy(s) ;
  graphTextBounds (80, ++dnacpt->line) ;
  graphPop() ;
  graphRedraw () ;
}

static void dnacptRestriction2 (char *text)
{  dnacptRestriction () ; } /* to please graphCompletionEntry */

static void dnacptSelectEnz(void)
{ OBJ obj ;
  int i , n ; char  *cname ;
  KEY key ;
  KEY _VOligo =  pickWord2Class("Oligo") ;
  KEYSET kSet ;
  void *look ;
  DNACPTGET("enzSet") ;

  getMatchSeqTag () ;		/* make sure matchSeqTag is set */

  if (!keySetActive(&kSet, &look)) 
    { messout("First select a keySet (containing Restriction enzymes)") ;
      return ;
    }
  keySetDestroy (dnacpt->enz) ;
  dnacpt->enz = keySetCreate() ;
  i = keySetMax(kSet) ; n = 0 ;
  while (i--)
    { key = keySet(kSet,i) ;
      if (class(key) == MYMOTIF
	  && (obj = bsCreate(key)))
	{ if (bsGetData(obj,matchSeqTag,_Text,&cname))
	    keySet(dnacpt->enz, n++) = key ;
	  bsDestroy(obj) ;
	}
      else if (_VOligo && class(key) == _VOligo
	  && (obj = bsCreate(key)))
	{ if (bsGetData(obj,_Sequence,_Text,&cname))
	    keySet(dnacpt->enz, n++) = key ;
	  bsDestroy(obj) ;
	}
    }
  if (!n)
    { messout("The selected keySet does not contain any Oligo, Enzyme or other Motif, sorry") ;
      keySetDestroy(dnacpt->enz) ;
      return ;
    }
  arraySort(dnacpt->enz,keySetAlphaOrder) ;
  sprintf(dnacpt->restriction,"%d sites", n) ; 
  strcpy(dnacpt->restriction2,dnacpt->restriction) ;
  graphBoxDraw(dnacpt->restrictionBox, -1, -1);
  graphCompletionEntry (restrictionCompletion, dnacpt->restriction,0,0,0,0) ; /* make it acive */
}

/***********************************************************/
/***********************************************************/

static void dnacptEmblDump (void)
{
  Array dna, colors ; 
  void  *look ;
  KEY	key ;
  KEYSET kSet ;
  extern void emblDumpKey (KEY key) ;
  extern void emblDumpKeySet (KEYSET ks) ;

  DNACPTGET("dnacptEmblDump") ;
  
  if (dnacpt->useKeySet)
    {
      if (keySetActive (&kSet, &look))
	emblDumpKeySet (kSet) ;
      else
	messout("First select a keySet containing sequences") ;
    }
  else
    {
      if (fMapActive (&dna,&colors,&key,&look))
	emblDumpKey (key) ;
      else
	messout ("First select a DNA window or the "
		 "'KeySet' button") ;
    }

  return;
} /* dnacptEmblDump */

void dnacptFastaDump(void)
{
  Array dna ;
  Array colors ;
  KEY  seqKey ;
  int from = 0, to, origin ;
  FILE *fil ;
  DNACPTGET("dnacptFastaDump") ;
  
  if (dnacpt->useKeySet)
    {
      KEYSET kSet ;
      void *look ; 
            
      if (!keySetActive(&kSet, &look)) 
	{ messout("First select a keySet containing sequences") ;
	  return ;
	}
      if ((fil = dnaFileOpen ()))
	{ dnaDumpFastAKeySet (kSet, fil, 0) ;
	  filclose (fil) ;
	}
    }
  else
    {
      OBJ obj ;
      KEY titleKey = 0, father = 0 ;
      void *look ;
      
      if (!fMapActive (&dna, &colors, &seqKey, &look))
	{ messout ("First select a DNA window or the "
		   "'KeySet' button") ;
	  return ;
	}
      fMapFindZone (look, &from, &to, &origin) ;
      
      if (!(fil = dnaFileOpen()))
	return ;
      if (fMapFindZoneFather (look, from, to, &father, &origin))
	{  
	  dnaDumpFastA (dna, from, to-1, 
			messprintf ("%s %d-%d", name(father), 
				    origin+1, origin + to - from), 
			fil, 0) ;
	}
      else if ((obj = bsCreate(seqKey)))
	{ if (bsGetKey (obj, _Title, &titleKey))
	    dnaDumpFastA (dna, from, to-1, 
			  messprintf ("%s %d-%d %s", name(seqKey), 
				      from+1, to, name(titleKey)), 
			  fil, 0) ;
	  else
	    dnaDumpFastA (dna, from, to-1, 
			  messprintf ("%s %d-%d", name(seqKey), 
				      from+1, to), 
			  fil, 0) ;
	  bsDestroy (obj) ;
	}
      else
	dnaDumpFastA (dna, from, to-1, 
		      messprintf ("%s %d-%d", name(seqKey), 
				  from+1, to), 
		      fil, 0) ;

      filclose (fil) ;
    }
}

/*****************************************************************/
/*****************************************************************/

void dnacptMakeTestSequences(void)
{
  static char dnaDecodeSpec[] =
    { 'a', 't', 'g', 'c', 'a', 't', 'g', 'c', 'n', 'n', 'n', 'n' } ;
  static char fname[FIL_BUFFER_SIZE], dname[DIR_BUFFER_SIZE] ;
  Array dna ;
  Array colors ;
  KEY  seqKey ;
  void *look = 0 ;
  char *cp ;
  int from = 0, to, origin, ll, cover, y, u, max, nn, i, x ;
  int offset, taux, type ;
  unsigned char *ip ;
  FILE *fil ;
  char cq[12] ;
  DNACPTGET("dnacptFastaDump") ;
  
  if (dnacpt->useKeySet)
    { messout("Can't do that from keySet") ;
      return ;
    }
  
  if(!fMapActive(&dna,&colors,&seqKey, &look))
    { messout ("First select a DNA window or the "
	       "'KeySet' button") ;
      return ;
    }
  fMapFindZone (look, &from, &to, &origin) ;
   
  if (!messPrompt("Nameprefix, length, coverage and error rate","a 300 10 2","wiiiz"))
	return ;

  cp = freeword() ;
  strncpy(cq, cp, 8) ; 
  freenext() ;
  freeint(&ll) ; freenext() ;  freeint (&cover) ; freeint(&taux) ;
  if (ll < 40)
    {
      messout("length %d too short") ;
      return ;
    }
  if (taux > 20 || taux < 0)
    {
      messout("please choose a value for error rate between 0 and 20") ;
      return ;
    }
  if (taux)
    {
      taux = 200 / taux ;
      offset = randint() % taux ;
    }
  else offset = 10 ;
  
  if (!(fil = filqueryopen(dname, fname, "ace", "w", "Export")))
    return ;
  
  if(ll < 100) ll = 100 ;
  max = to - from ;

  dnaRepaint(colors) ; 
  nn = max*cover / ll ;
  
  while(nn--)
    { x = randint() ;
      x <<= 8 ; x += randint() ;
      if (x < 0)
	x = - x ;      
      x %= max ;
      y = randint() % 2 ;
      fprintf(fil, "\n\nDNA : %s", cq) ;

      if (y)
	{ u = x ;	 
	}
      else
	{ u = x - ll ;
	  if (u < 1 ) u = 1 ;
	}
      while (u)
	{ fprintf(fil, "%c", '0' + u % 10) ;
	  u /= 10 ;
	}
      fprintf(fil, "%c\n", y ? 'D' : 'R') ;
      
      if (y)
	{ i = ll ;
	  if (x + i > max)
	    i = max - x ;
	  for (cp = arrp(dna, from + x, char), ip = arrp(colors, from + x, unsigned char) ;
	       i ;
	       i--, cp++, ip++)
	    { 
	      if (!offset)
		{
		  offset = randint() % taux ;
		  type = randint() % 16 ;
		  switch (type)
		    {
		    case 0:case 1:case 2:case 3:case 8:case 9:case 10:case 11:
		      putc(dnaDecodeSpec[type], fil) ;
		      break ;
		    case 4:case 5:case 6: case 7:
		      cp-- ;
		      putc(dnaDecodeSpec[type], fil) ;
		      break ;
		    default:
		      break ;
		    }
		}
	      else
		{
		  if (taux)
		    offset-- ;
		  *ip = TINT_CYAN ;
		  putc((int)dnaDecodeChar[(int)*cp], fil) ;
		}
	    }
	}
      else
	{ i = ll ;
	  if (x - i <= 0)
	    i = x ;
	  for (cp = arrp(dna, from + x, char) , ip = arrp(colors, from + x, unsigned char) ;
	       i ;
	       i--, cp--, ip-- )
	    { 
	      if (!offset)
		{
		  offset = randint() % taux ;
		  type = randint() % 16 ;
		  switch (type)
		    {
		    case 0:case 1:case 2:case 3:case 8:case 9:case 10:case 11:
		      putc(dnaDecodeSpec[type], fil) ;
		      break ;
		    case 4:case 5:case 6: case 7:
		      cp++ ;
		      putc(dnaDecodeSpec[type], fil) ;
		      break ;
		    default:
		      break ;
		    }
		}
	      else
		{
		  if (taux)
		    offset-- ;
		  *ip = TINT_LIGHTGREEN ;
		  putc((int)dnaDecodeChar[(int)complementBase(*cp)], fil) ;
		}
	    }
	}
     } 
      
  ip = arrp(colors, 0, unsigned char) ;
  i = arrayMax(colors) ; x = 0 ;
  *ip = 0 ;
  while(ip++, i--)
    if (*ip && !(*(ip - 1))) x++ ;
  messout("I exported %d contigs", x) ;

  fMapReDrawDNA (look) ;
  filclose (fil) ;
}

/*************************************************************************/
/**************************************************************************/

static void dnaAG_TC(void)
{
  register int i , j ;
  int this = 0, old = 0 ;
  int n ;
  Array plot ;
  Array dna, colors ; 
  KEY key ;
  char *cp , keep ;
  void *look ;
  static char fileName[FIL_BUFFER_SIZE],dirName[DIR_BUFFER_SIZE] ;
  static  FILE *f = 0 ;
  DNACPTGET("dnaAG_TC") ;

  if(!fMapActive(&dna,&colors,&key,&look))
    { messout("First select a dna window") ;
      return ;
    }
  plot = arrayCreate(50,int) ;
  dnacptDisplay() ;
  dnaRepaint(colors) ;
  cp = arrp(dna,0,char) ;
  n =arrayMax(dna) ;
  
  strncpy(fileName, name(key), 22) ;
  f= filqueryopen (dirName, fileName, "agct", "w","Where to export the agct results ?") ;
  if(f)
    fprintf(f,"\n\n %s\n",name(key)) ;

  for(i = j = 0 ; i<n ; cp++, i++ )
    { if(*cp & R_)  /* R is A|G */
	this = 1 ;
      if(*cp & Y_)  /* Y is T|C */
	  this = 2 ;
      if(this == old)
	j++;
      else
	{ if(j>12)
	    { array(plot,j,int)++ ;
	      if(f)
		fprintf(f,"%d pos %d\n  ",j, i-j) ;
	      keep = *(cp+2);
	      *(cp+2) = 0 ;
	      if(f)
		fprintf(f," %s\n",
			dnaDecodeString(cp-j-(i>1 ? 2 : i)));
	      *(cp+2) = keep ;
	      while(j--)
		arr(colors,i-j-1,char) = RED ;
	    }
	  j = 1 ;
	}
      old = this ;
      this = 0 ;
    }

  if(f)
    filclose(f) ;

  /*  plotHisto("Sites found",plot) ;
      rplaced by arrayDestroy : do not do both */
  arrayDestroy(plot) ;

  fMapReDrawDNA (look) ;
}

/*************************************************************************/
/*************************************************************************/
/**********   Consensus splicing sequence ********************************/
   /* accumulates in sites5 (resp sites), a bit of dna around from (resp to) 
      complementing if from > to
      */
      
static int dnacptAddSplice (Array dna, Array  site5, Array site3, 
			    int from, int to, int lengthEx, int lengthIn, char *filter)
{
  int i = from ; char *cp ; int *ip ; int n = lengthEx + lengthIn + 1 ;
  char *f = filter ;

  if (from < to)
    {
      cp = arrp(dna, 0, char) + from + 1 ;
      if (!((*cp) == *f) || ! (*(cp+1) == *(f+1)))
	return 0 ;
      cp = arrp(dna, 0, char) + to - 2 ;
      if (!((*cp) == *(f+3)) || ! (*(cp+1) == *(f+4)))
	return 0 ;
      for(i = from - lengthEx, cp = arrp(dna, 0, char) + i , ip = arrayp(site5, 0, int);
	  i <= from + lengthIn && i < arrayMax(dna); cp ++, i++, ip++ )
	{ 
	  if (i<0) continue ;
	  switch(*cp)
	    {
	    case A_ : (*ip)++ ; break ;
	    case C_ : (*(ip + n))++ ; break ;
	    case G_ : (*(ip + 2*n))++ ; break ;
	    case T_ : (*(ip + 3*n))++ ; break ;
	    default: (*(ip + 4*n))++ ; break ;
	    }
	}
      for(i = to - lengthIn, cp = arrp(dna, 0, char) + i , ip = arrp(site3, 0, int);
	  i <= to + lengthEx && i < arrayMax(dna); cp ++, i++, ip++ )
	{ if (i<0) continue ;
	  switch(*cp)
	    {
	    case A_ : (*ip)++ ; break ;
	    case C_ : (*(ip + n))++ ; break ;
	    case G_ : (*(ip + 2*n))++ ; break ;
	    case T_ : (*(ip + 3*n))++ ; break ;
	    default: (*(ip + 4*n))++ ; break ;
	    }
	}
    }
  else
    {
      cp = arrp(dna, 0, char) + from ;
      if (!((*cp) && complementBase((*f))) || ! (*(cp+1) && complementBase((*(f+1)))))
	return 0 ;
      cp = arrp(dna, 0, char) + to ;
      if (!((*cp) && complementBase((*(f+3)))) || ! (*(cp+1) && complementBase((*(f+4)))))
	return 0 ;
      for(i = from + lengthEx , cp = arrp(dna, 0, char) + i , ip = arrp(site5, 0, int);
	  i >= from - lengthIn  && i >= 0 ; cp--, i--, ip++ )
	{ 
	  if (i >= arrayMax(dna)) continue ;
	  switch(*cp)
	    {
	    case T_ : (*ip)++ ; break ;
	    case G_ : (*(ip + n))++ ; break ;
	    case C_ : (*(ip + 2*n))++ ; break ;
	    case A_ : (*(ip + 3*n))++ ; break ;
	    default: (*(ip + 4*n))++ ; break ;
	    }
	}
      for(i = to + lengthIn , cp = arrp(dna, 0, char) + i , ip = arrp(site3, 0, int);
	  i >= to - lengthEx  && i >= 0; cp--, i--, ip++ )
	{ if (i >= arrayMax(dna)) continue ;
	  switch(*cp)
	    {
	    case T_ : (*ip)++ ; break ;
	    case G_ : (*(ip + n))++ ; break ;
	    case C_ : (*(ip + 2*n))++ ; break ;
	    case A_ : (*(ip + 3*n))++ ; break ;
	    default: (*(ip + 4*n))++ ; break ;
	    }
	}
    }
  return 1 ;
}

/*************************************************************************/

    /* display the cumulted splice sites in the current gmapcpt window */
static void  dnacptShowSpliceConsensus (DNACPT dnacpt, char *title, 
					Array  site,int  lengthEx, int lengthIn) 
{ 
  int n = lengthIn + lengthEx + 1 , i , j ;
  graphText(title, 5, dnacpt->line += 3) ;
  for (j=0; j<5; j++)
    for(i=0; i< n; i++)
      graphText(messprintf("%7d", array(site, n*j + i, int)),
		7 + 7*i, dnacpt->line + 2 + j) ;
  graphText("A", 2, dnacpt->line + 2) ;
  graphText("C", 2, dnacpt->line + 3) ;
  graphText("G", 2, dnacpt->line + 4) ;
  graphText("T", 2, dnacpt->line + 5) ;
  graphText("X", 2, dnacpt->line + 6) ;

  if(*title == '5')
    { graphText("--- intron ---> ", 7 * lengthEx + 12, dnacpt->line + 1) ;
      graphLine(7*lengthEx + 13.05, dnacpt->line + 1.5, 
		7*lengthEx + 13.05, dnacpt->line + 7.5) ;
    }
  else
    { int p = 7 * lengthIn - 14 ;
      if (p<0) p = 1 ;
      graphText("--- intron --->", p , dnacpt->line + 1 ) ;
      graphLine(7*lengthIn + 7.05, dnacpt->line + 1.5, 
		7*lengthIn + 7.05, dnacpt->line + 7.5) ;
    }
  dnacpt->line += 8 ;

  graphText ("Relative usage", 2,  dnacpt->line) ;
  dnacpt->line += 2 ;

  
  graphRedraw() ;
}

/*************************************************************************/
    /* Gets the dna of the parent sequence 
     */  
static Array  dnacptGetSpliceSitesDna(KEY  gene, int *a1p, int *a2p, BOOL tg)
{
  int a1 = 0, a2 = 0 ;
  OBJ Gene = 0 , Parent = 0 ;
  KEY parent = 0 ;
  Array dna = 0 ;
  BOOL ok = FALSE ;

  Gene = bsCreate (gene) ;
  if (!tg && Gene)
   {
     bsGetKey (Gene, _Source, &parent) ;
     if (parent && (Parent = bsCreate (parent)) &&
	 bsFindKey (Parent, _Subsequence, gene) &&
	 bsGetData (Parent, _bsRight, _Int, &a1) &&
	 bsGetData (Parent, _bsRight, _Int, &a2) &&
	 (dna = dnaGet (parent)) &&
	 a1 >= 0 && a2 >= 0 &&
	 a1 < arrayMax(dna) && a2 < arrayMax(dna))
       ok = TRUE ;
   }
  else if (tg && Gene)
   {
     bsGetKey (Gene, str2tag ("Genomic_sequence"), &parent) ;
     if (parent && (Parent = bsCreate (parent)) &&
	 bsFindKey (Parent, str2tag ("Transcribed_gene"), gene) &&
	 bsGetData (Parent, _bsRight, _Int, &a1) &&
	 bsGetData (Parent, _bsRight, _Int, &a2) &&
	 (dna = dnaGet (parent)) &&
	 a1 >= 0 && a2 >= 0 &&
	 a1 < arrayMax(dna) && a2 < arrayMax(dna))
       ok = TRUE ;
   }
  
  if (!ok)
    arrayDestroy (dna) ;
  else if (a1 > a2)
    {
      a1 = (long)arrayMax(dna) + 1 - a1 ;
      a2 = (long)arrayMax(dna) + 1 - a2 ;
      reverseComplement (dna) ;
    }
  a1-- ;  /* plato */
  
  bsDestroy (Parent) ;
  bsDestroy (Gene) ;
  *a1p = a1 ; *a2p = a2 ;
  return dna ;
}

/*************************************************************************/
    /* Gets the sites recursivelly 
     *  taking orientation into account
     */  
static int  dnacptGetSpliceSites(Array sites, KEY  seq, int from, int to, BOOL tg, BOOL isAlter)
{ int i, n, nfound = 0, x1, x2 ;
  OBJ Seq = bsCreate(seq) ;
  Array exons = arrayCreate(8, BSunit) ;

  if(Seq && !tg && Seq && bsFindTag(Seq, _Source_Exons) &&
     bsFlatten(Seq, 2, exons))
    { dnaExonsSort (exons) ;
      for (i = 1; i + 1  < arrayMax(exons) ; i += 2) /* introns are between exons */
	{ 
	  x1 = arr(exons, i, BSunit).i - 1 ;
	  x2 =arr(exons, i + 1, BSunit).i - 1 ;

	  n = arrayMax(sites) ;
	  array(sites, n, int) = from + x1 ;
	  array(sites, n + 1, int) = from + x2 ;
	  nfound = 1 ;
	}
    }
  else if(Seq && tg && Seq && bsFindTag(Seq, str2tag ("Splicing")) &&
	  bsFlatten(Seq, 3, exons))
    { 
      for (i = 0; i  < arrayMax(exons) ; i += 3) /* introns are tagged in column 3 */
	{ 
	  int  n = arrayMax(sites) ;
	  KEY tt = arr(exons, i + 2, BSunit). k ;
	  if ((!isAlter && tt == str2tag ("Intron")) ||
	      (isAlter && strstr (name(tt), "ntron") && strstr (name(tt), "ltern")))
	    { 
	      x1 = arr(exons, i, BSunit).i  - 2 ;
	      x2 = arr(exons, i + 1, BSunit).i ;

	      n = arrayMax(sites) ;
	      array(sites, n, int) = from + x1 ;
	      array(sites, n + 1, int) = from + x2 ;
	      nfound = 1 ;
	    }
	}
    }

  arrayDestroy(exons) ;
  bsDestroy(Seq) ;
  return nfound ;
}

/*************************************************************************/
  /* Given tha active key set
   * look for dnasequences whose descriptor is in keyset
   * in each, look recursively for the splice sites,
   *  accumulate them and display.
   * by starting from _VDNA, we are garanteed not to double count.
   */
static void dnacptDoSpliceConsensus(BOOL tg)
{
  KEY gene = 0 ;
  void *look ;
  char *cp, *cq, filter[6], filterDecoded[6] ;
  int lengthEx = 8, lengthIn = 12 , i , ii, nn = 0 , nintrons = 0, a1, a2 ;
  KEYSET kSet = 0, ks2 = 0; 
  BOOL isalter = FALSE ;
  Array dna = 0, sites = 0,  site5 = 0, site3 = 0;
  DNACPTGET("dnacptSpliceConsensus") ;
      
  if(!dnacpt->useKeySet)
    { messout("I just know how to do that from a keySet") ;
      return ;
    }
  if (!keySetActive(&kSet, &look)) 
    { messout("First select a keySet containing genes") ;
      return ;
    }

  {
    if (!messPrompt("concensus","gt_ag alt","w")) 
      return ;
    cp = cq = freeword() ;
    if (strlen (cp) == 5 &&
	(*cp != 'a' || *cp != 't'  || *cp != 'g'  || *cp != 'c')  &&
	(cp++) &&
	(*cp != 'a' || *cp != 't'  || *cp != 'g'  || *cp != 'c')  &&
	(cp++) &&
	(*cp == '_') &&
	(cp++) &&
	(*cp != 'a' || *cp != 't'  || *cp != 'g'  || *cp != 'c')  &&
	(cp++) &&
	(*cp != 'a' || *cp != 't'  || *cp != 'g'  || *cp != 'c')  
	)
      {
	strcpy (filter, cq) ;
	strcpy (filterDecoded, cq) ;
	cp = filter ;
	*cp =  dnaEncodeChar[(int)(*cp)] ; cp++ ;
	*cp =  dnaEncodeChar[(int)(*cp)] ; cp++ ;cp++ ;
	*cp =  dnaEncodeChar[(int)(*cp)] ; cp++ ;
	*cp =  dnaEncodeChar[(int)(*cp)] ;
      }
    else
      {
	messout("I only understand 5 lower cases atgcn like gc_ag") ;
	return ;
      }
    if (freeword())
      isalter = TRUE ;
  }
    
  ks2 = keySetCreate() ;
  sites = arrayCreate(100, int) ;
  site5 =arrayCreate(5*(lengthIn  + lengthEx + 1), int) ; /* 5 is atgc + ambiguous */
  site3 =arrayCreate(5*(lengthIn + lengthEx), int) ;     /* 5 is atgc + ambiguous */
  array(site5,5*(lengthIn  + lengthEx + 1) - 1 , int) = 0 ;
  array(site3,5*(lengthIn  + lengthEx + 1) - 1 , int) = 0 ;
  for (ii = 0 ; ii < keySetMax(kSet) ; ii++)
    {
      gene = keySet (kSet, ii) ;
      dna = dnacptGetSpliceSitesDna (gene, &a1, &a2, tg) ;
      if (dna)
	{
	  nn++ ;
	  keySetInsert(ks2,gene) ;
	  arrayMax(sites) = 0 ;
	  dnacptGetSpliceSites(sites, gene,  a1, a2, tg, isalter) ;
	  for (i=0 ; i<arrayMax(sites) ; i+= 2)
	     nintrons += dnacptAddSplice(dna, site5, site3, 
			    arr(sites,i,int), arr(sites, i+1, int), lengthEx, lengthIn, filter ) ;
	  arrayDestroy (dna) ;
	}
    }
        
  if (!nn)
    { messout("First select a keySet containing genes") ;
      return ;
    }
  graphText(messprintf("I scanned %d genes, containing %d %s %s introns",
		       nn, nintrons, filterDecoded
		       , isalter ? "alternative" : "non alternative"), 
	    2,dnacpt->line += 2) ;
  dnacptShowSpliceConsensus(dnacpt, "5\'consensus",site5, lengthEx, lengthIn) ;
  dnacptShowSpliceConsensus(dnacpt, "3\'consensus",site3, lengthEx, lengthIn) ;
  nn = (lengthIn + lengthEx + 1) * 7 + 3 ;
  if (nn < 80) nn = 80 ;
  graphTextBounds (nn,dnacpt->line += 2) ;

  arrayDestroy(site3) ;
  arrayDestroy(site5) ;
  arrayDestroy(sites) ;
  keySetDestroy(kSet) ;
  keySetShow(ks2, look) ;
}

static void dnacptSpliceConsensus(void)
{ dnacptDoSpliceConsensus(FALSE) ; }
static void dnacptConfirmedSpliceConsensus(void)
{ dnacptDoSpliceConsensus(TRUE) ; }

/*************************************************************************/
/*************************************************************************/
/******************   Codon usage         ********************************/
   /* this is a static because it follows the usual UCGA display of the
      genetic code with codon in 0, 65 ; rather than our usual A_, T_, G_, C_
      internal convention, all ambiguos codons are treated as 64.
      */

static char* codonFullName[64] = {
  "Phe","Phe","Leu","Leu","Ser","Ser","Ser","Ser","Tyr","Tyr","***","***","Cys","Cys","***","Trp",
  "Leu","Leu","Leu","Leu","Pro","Pro","Pro","Pro","His","His","Gln","Gln","Arg","Arg","Arg","Arg",
  "Ile","Ile","Ile","Met","Thr","Thr","Thr","Thr","Asn","Asn","Lys","Lys","Ser","Ser","Arg","Arg",
  "Val","Val","Val","Val","Ala","Ala","Ala","Ala","Asp","Asp","Glu","Glu","Gly","Gly","Gly","Gly" } ;

/**********/
  /* Gives the per thousandth of usage of synonymous codons */
static void dnacptRenormaliseUsage(Array usage)
{
  int i, j ;
  int equiv[64] , tot[64] ;
  char name[5] ;
 
  for(i = 0 ; i < 64; i++)
    equiv[i] = -1 ;

  for(i = 0 ; i < 64; i++)
    if(equiv[i] == -1)
      { tot[i] = 0 ;
	strcpy(name, codonFullName[i]) ;
	for (j = i; j <64 ; j++)
	if(!strcmp(name, codonFullName[j]))
	  { tot[i] += arr(usage,j, int) ;
	    equiv[j] = i ;
	  }
      }

  for(i = 0 ; i < 64; i++)
    if(tot[equiv[i]])
      arr(usage, i, int) = 
	(1000 * arr(usage, i, int)) / tot[equiv[i]] ;
}

/**********/

typedef struct LOOKusageStruct
{ void *magic ;
  Graph graph ;
  Array usage ;
  int myUsage ;
} *LOOKusage  ;

static void dnacptShowUsageDestroy(void)
{ LOOKusage look = 0 ;

  if (graphAssFind (&MAGIC_CODON_USAGE, &look) &&
      look->magic == &MAGIC_CODON_USAGE )
    { if (look->myUsage == arrayExists(look->usage))
	arrayDestroy(look->usage) ; 
      look->magic = 0 ;
      messfree(look) ;
    }
}

/********/

static int aaUse(Array usage, char *cp, int total)
{
  int  subtotal = 0 , i ;
 
  if(arrayMax(usage) >= 64){
    for(i = 0 ; i < 64; i++)
      if (codonFullName[i] && !strcmp(cp, codonFullName[i]))
	subtotal += arr(usage,i, int) ;
    
    return subtotal ;
  }
  else return 0;
}

/********/

static void showAA(char *text, int n, int x , int  y, int total)
{
  if (total)
    graphText(messprintf("%s%5.1f%%", text, (100.0 * n)/total), x, y) ;
  else
    graphText(messprintf("%s  -", text), x, y) ;
}

/********/

static int dnacptShowSomePeptide(Array usage, int total,
				 char *cp, int x, int y)
{ int i, n = 0, nn = 0 ; char buf[4] ;
  for (i = 0 ; *cp ; i++, cp+= 3)
    { 
      strncpy(buf,cp,3) ; buf[3] = 0 ;
      n = aaUse(usage, buf, total) ; nn += n ;
      showAA(buf, n, 8 + 11*i, y , total ) ; 
    }
  return nn ;
}

/********/

static int  dnacptShowPeptideUsage(Array usage, int y, int total) 
{ int  n, subTotal = 0 ; 

  graphText("Amino acids usage:", 2, y) ; y += 2 ;
  
  n = dnacptShowSomePeptide(usage, total, "LysArgHis", 8, y + 2 ) ;
  showAA("Basic", n, 4, y + 1 , total) ;
  subTotal += n ;
 
  n = dnacptShowSomePeptide(usage, total, "AspGlu", 8, y + 4 ) ;
  showAA("Acidic", n, 4, y + 3 , total) ;
  subTotal += n ;

  n = dnacptShowSomePeptide(usage, total, "AsnGlnCysMetSerThr", 8, y + 6 ) ;
  showAA("Neutral", n, 4, y + 5 , total) ;
  subTotal += n ;
  showAA("Hydrophilic", subTotal, 1, y  , total) ;
  subTotal = 0 ;

  n = dnacptShowSomePeptide(usage, total, "GlyAlaValProLeuIle", 8, y + 10 ) ;
  showAA("Aliphatic", n, 4, y + 9 , total) ;
  subTotal += n ;

  n = dnacptShowSomePeptide(usage, total, "PheTyrTrp", 8, y + 12 ) ;
  showAA("Aromatic", n, 4, y + 11, total) ;
  subTotal += n ;
  showAA("Hydrophobic", subTotal, 1, y + 8 , total) ;
  subTotal = 0 ;

  return 21 ;  /* Overall number of lines used. */
}

/********/

static void dnacptShowCodonUsage (int nn, int ncds, Array usage)
{
  int n = arrayMax(usage), total = 0, i, j, k, codon ;
  int y = 4 ;
  static char name[] = {'U','C','A','G'} ;
  char buffer[40] ;
  Graph old = graphActive () ;
  LOOKusage look = (LOOKusage) messalloc(sizeof(struct LOOKusageStruct)) ;

  /* displayCreate(DtCodons) ; -g Text_FIT missing in acembly wspec, so use graphCreate() */
  graphCreate(TEXT_FIT,"Codon Usage",.1, .5, .7, .85) ;
  graphHelp("Codon_usage") ;
  graphPop() ;
  graphRegister(DESTROY, dnacptShowUsageDestroy) ;
  graphRegister(RESIZE, graphRedraw) ;
  graphAssociate (&MAGIC_CODON_USAGE, look) ;
  look->magic = &MAGIC_CODON_USAGE ;
  look->usage = usage ;
  look->myUsage = arrayExists(usage) ;
  look->graph = graphActive() ;

  graphText (messprintf ("I scanned %d sequences, containing %d coding sequences",
			 nn, ncds), 2, 2) ;

  n = arrayMax(usage) ;
  while (n--)
    total += arr (usage, n, int) ;
    
  y += dnacptShowPeptideUsage(usage, y, total) ;

  graphText ("Codon usage.", 2, y - 4) ;
  /* graphText("Percentage of codon occurence.", 2, y - 4) ; */

  graphText (messprintf ("Total %d codons, %d stops", total, 
		 arr(usage,10,int) + arr(usage,11,int) + arr(usage,14,int)),
	     40, y - 4) ;
  i = array(usage, 64, int) ;
  if (i)
    graphText (messprintf ("%d ambiguous codons.", i), 10, y - 3) ;
  n = arrayMax(usage) ;
  if(n > 65)
    { /* if there are codons uncomplete  */
      i = array(usage, 65, int) ;
      if (i)
	graphText (messprintf ("%d uncomplete codons.", i), 40, y - 3) ;
    }
  for (i = 0 ; i < 5; i++)
    graphLine(5, y + 8 * i, 73, y + 8*i) ;
  for (i = 0 ; i < 5; i++)
    graphLine(5 + 17 *i , y , 5 + 17 * i, y + 32 ) ;
  
  for(i=0; i<4; i++)
    graphText(messprintf("%c", name[i]), 13 + 17 * i, y  -1.3) ;
  for(i=0; i<4; i++)
    graphText(messprintf("%c", name[i]), 3, y + 4 + 8 * i) ;

  dnacptRenormaliseUsage(usage) ;
  for(i = 0; i<4 ; i++)
    for ( j = 0 ; j < 4 ; j++)
      for ( k = 0 ; k < 4 ; k++)
	{ codon = 16*i + 4*j + k ;
	  buffer[0] = name[i] ;
	  buffer[1] = name[j] ;
	  buffer[2] = name[k] ;
	  buffer[3] = buffer[4] = ' ' ;
	  buffer[5] = 0 ;
	  strcat(buffer, codonFullName[codon]) ; /* 5, 6, 7 */
	  buffer[8] = 0 ;
	  strcat(buffer, messprintf("%5.1f%%", arr(usage,codon,int)/ 10.0)) ;
	  graphText(buffer, 7 + 17*j, y + 2 + 8*i + k + (k/2) ) ;
	}
  graphRedraw() ;
  graphActivate (old) ;
}

/*************************************************************************/
   /* accumulates in usage the codon used in dna messenger */      
static void dnacptAddCds2 (Array usage,  Array dna)
{
  int i = arrayMax(dna), codon = 0 , frame = 0 ;
  char *cp = arrp(dna, 0, char) ;
  
  while (i--)
    { 
      if (codon != 64)
	switch (*cp++)
	  {
	  case T_: break ;
	  case C_: codon += 1 ; break ;
	  case A_: codon += 2 ; break ;
	  case G_: codon += 3 ; break ;
	  default: codon = 64 ; break ;
	  }

      frame++ ;
      frame %= 3 ;
      if (frame)
	{ if (codon != 64) codon <<= 2 ; }
      else 
	{ array(usage, codon, int)++ ;
/*   Stops should not occur 
	  if (codon == 10 || codon == 11 || codon == 14 )
	     invokeDebugger() ;
*/
	  codon = 0 ;
	}
    }    
  if (frame)
    array(usage, 65, int)++ ;
}

/****************************************************************************/

static void dnacptWordUsage(void)
{
  Array  histo = arrayCreate (1000, int) ;
  Array words = 0 ;
  KEYSET kSet ; void *look ;
  int i, nn ;
  unsigned char *cp = 0 ;
  const char *ccp ;
  unsigned long nb = 0, nw = 0 ;
  static char fname[FIL_BUFFER_SIZE] = { 0 } , dname[DIR_BUFFER_SIZE] = {'.', 0 } ;
  char buf [2 * FIL_BUFFER_SIZE + 3 ] ;
  FILE *fil ;
  DNACPTGET("dnacptWordUsage") ;
  
  if(!dnacpt->useKeySet)
    { messout("I just know how to do that from a keySet") ;
    return ;
    }
  if (!keySetActive(&kSet, &look)) 
    { messout("First select a keySet containing sequences") ;
    return ;
    }
  if (!fname[0])
    strcpy (fname, "wordFrequency.15") ; 
  if (!(fil = filqueryopen(dname, fname, "15", "w", "Export the 500 Mb word frquency table")))
    ccp = 0 ;
  else
    {
      ccp = buf ;
      strcpy (buf, messprintf ("%s/%s", dname, fname)) ;
    }
	    
  words = dnaGetWordUsage (kSet, 1, &nb, &nw, TRUE, ccp, 0) ;
  i = arrayMax(words) ;
  cp = arrp (words, 0, unsigned char) - 1;
  while (cp++, i--)
    {
      nn = (int) (*cp & 0xf) ;
      if (nn > 0)
	array(histo, nn, int)++ ;
      nn = (int) (((*cp) >> 4) & 0xf) ;
      if (nn > 0)
	array(histo, nn, int)++ ;
    }
  if (arrayMax(histo))
    plotHisto (messprintf ("y = #15mer seen x times, %d base", nb) ,histo) ;
  else
    arrayDestroy (histo) ;
  if (!nw)
    { messout("First select a keySet containing sequences") ;
      return ;
    }

  arrayDestroy (words) ;
}

/****************************************************************************/
   /* construct the plus and back messengers  dna1 and dna2 */

static int dnacptAddCds (Array dna, Array cds, Array usage,KEY dnaKey)
{ int ncds = 0, from, to, i, j ;
  static Array gene = 0, exon = 0 ; 
  int n = 0, ne = 0 ;
  char *cp ;
#ifdef CODE_FILE
  KEY seq = 0 ;
#endif
  int kExon ;

  gene = arrayReCreate(gene, 1000, char) ;
  exon = arrayReCreate(exon, 300, char) ;

  for (i=0 ; i < arrayMax(cds) ; i+= 2)
    {
      from = arr(cds,i,int) ; to = arr(cds, i+1, int) ;

      if (from == -1)   /* begin of seq */
	{ 
#ifdef CODE_FILE
	  seq = to ;
#endif
	  kExon = 1;
	  n = 0 ;
	  arrayMax(gene) = 0 ;
	  continue ;
	}
      
      if (from == -2)   /* end of seq */
	{ ++ncds ; 
	  dnacptAddCds2 (usage, gene) ;
#ifdef CODE_FILE
	  dnaDecodeArray (gene) ;
	  array(gene, arrayMax(gene), char) = 0 ;
	  fprintf (codeFile, "%s %s\n", name(to), arrp(gene,0,char)) ;
#endif	     
	  continue ;
	}

      if(from >= arrayMax(dna)){
	return 0;
      }
      cp = arrp(dna, from, char) ; /* copy dna */
      if (from <= to)
	for (j = from ; j <= to ; j++, cp++)
	  { if (ne || !(n%3))	/* keep in frame */
	      array(exon,ne++,char) = *cp ;
	    array(gene,n++,char) = *cp ;
	  }
      else
	for (j = from ; j >= to ; j--, cp--)
	  { if (ne || !(n%3))
	      array(exon,ne++,char) = complementBase(*cp) ;
	    array(gene,n++,char) = complementBase(*cp) ;
	  }

#ifdef CODE_FILE
      dnaDecodeArray (exon) ;
      array(exon, arrayMax(exon), char) = 0 ;
      fprintf (exonFile, "%s_%d %s\n", name(seq), kExon, arrp(exon,0,char)) ;
#endif	     
      arrayMax(exon) = 0 ;
      ne = 0 ;
      ++kExon ;
    }

  return ncds ;
}

/*************************************************************************/
    /* Gets the coding sequences recursivelly 
     *  taking orientation into account
     */  
static BOOL  dnacptGetCds(Array cds, KEY  seq, int from, int to)
{ int i, min, max, n ;
  BOOL result = FALSE ;
  OBJ Seq = bsCreate(seq) ;
  Array exons = arrayCreate(8, BSunit) ;
  Array subSeq = arrayCreate(9, BSunit) ;
  if (!Seq)
    return FALSE ;
  if (bsFindTag(Seq, _CDS))  /* seq is by itself coding sequences */
    { 
      min = 1 ; max = 0 ;	/* max 0 means go to the end of the spliced sequence */
      if (bsGetData(Seq, _bsRight, _Int, &min))
	bsGetData(Seq, _bsRight, _Int, &max) ;
      n = arrayMax(cds) ;
      array(cds, n++  , int) = -1 ;
      array(cds, n++ , int) = seq ; /* begin of seq */
      if (!bsFindTag(Seq, _Source_Exons) || !bsFlatten(Seq, 2, exons))
	{ array(exons, 0, BSunit).i = 1 ;  /* Failure: Make a fake entry */
	  array(exons, 1, BSunit).i = (from < to) ? to - from + 1 : from - to + 1 ;
	} 
      dnaExonsSort (exons) ;
      for (i = 0; i < arrayMax(exons) ; i += 2, n += 2 ) /*  exons themselves */
	{ result = TRUE ; 
	  if (from < to)        /* -1 because no zero in biology */
	    {
	      array(cds, n, int) = from + arr(exons,i,BSunit).i - 1 ;
	      array(cds, n + 1, int) = from + arr(exons,i + 1 ,BSunit).i - 1 ;
	    }
	  else
	    {
	      array(cds, n, int) = from - arr(exons,i,BSunit).i + 1 ;
	      array(cds, n + 1, int) = from - arr(exons,i + 1 ,BSunit).i + 1 ;
	    }
	}
      array(cds, n  , int) = -2 ;
      array(cds, n + 1 , int) = seq ; /* end of seq */
    }
  else  /* recurse though the sub sequences */
    {
      if (bsFindTag(Seq, _Subsequence) && bsFlatten(Seq, 3, subSeq))
	for (i = 0 ; i < arrayMax(subSeq) ; i += 3)
	  {
	    if (from < to)        /* -1 because no zero in biology */
	      dnacptGetCds(cds, arr(subSeq,i,BSunit).k , 
			   from + arr(subSeq,i+1,BSunit).i - 1,
			   from + arr(subSeq,i+2,BSunit).i - 1 ) ;
	    else
	      dnacptGetCds(cds, arr(subSeq,i,BSunit).k ,
			   from - arr(subSeq,i+1,BSunit).i + 1 ,
			   from - arr(subSeq,i+2,BSunit).i + 1 ) ;
	  }
    }
  arrayDestroy(subSeq) ;
  arrayDestroy(exons) ;
  bsDestroy(Seq) ;
  return result ;
}

/*************************************************************************/
  /* Given the active key set
   * look for dnasequnces whose descriptor is in keyset
   * in each, look recursively for the coding sequences,
   *  accumulate the codon usage and display.
   * by starting from _VDNA, we are garanteed not to double count.
   */

static void dnacptCodonUsage(void)
{
  Array dna ; KEY seq = 0 , dnaKey = 0 ;
  KEYSET kSet , ks2 ; void *look ;
  Array usage , cds ;
  int   ii,  n, nn = 0 , ncds = 0 ;
  DNACPTGET("dnacptCodonUsage") ;
      
  if(!dnacpt->useKeySet)
    { messout("I just know how to do that from a keySet") ;
      return ;
    }
  if (!keySetActive(&kSet, &look)) 
    { messout("First select a keySet containing sequences") ;
      return ;
    }

#ifdef CODE_FILE
  codeFile = fopen ("coding.seq","w") ;
  exonFile = fopen ("exons.seq","w") ;
#endif

  ks2 = keySetCreate() ;
  cds = arrayCreate(50,int) ;
  usage = arrayCreate(66, int) ;

  for (ii = 0 ; ii < keySetMax(kSet) ; ii++)
    if ((dna = dnaGet(seq=keySet(kSet,ii))))
      { nn++ ;
      arrayMax(cds) = 0 ;
      dnacptGetCds (cds, seq,  0, arrayMax(dna) - 1) ;
      n = dnacptAddCds(dna, cds, usage, dnaKey) ;
      if (n)
	{ ncds += n ;
	keySetInsert(ks2,seq) ;
	}
      arrayDestroy(dna) ;
      }
    

#ifdef CODE_FILE
  filclose (codeFile) ;
  filclose (exonFile) ;
#endif
      
  if (!nn)
    { messout("First select a keySet containing sequences") ;
      return ;
    }

  if(arrayMax(usage)  > 0)
    dnacptShowCodonUsage(nn, ncds, usage) ;
  else
    messout("No codon information found");

  /*   arrayDestroy(usage) ; destroyed in show.. */
  arrayDestroy(cds) ;
  /* 
     keySetDestroy(kSet) ;
     keySetShow(ks2, look) ;
     */
}

/*****************************************/
/*****************************************/

static void dnacptAminoTrue(void)
{ DNACPTGET("dnacptAminoTrue") ;
  dnacpt->amino = TRUE ;
  if (dnacpt->aminoButtonBox)
    {
      graphBoxDraw (dnacpt->aminoButtonBox, BLACK, LIGHTRED) ;
      graphBoxDraw (dnacpt->aminoButtonBox-1, BLACK, WHITE) ;
    }
}

static void dnacptAminoFalse(void)
{ DNACPTGET("dnacptAminoFalse") ;
  dnacpt->amino = FALSE ;
  if (dnacpt->aminoButtonBox)
    {
      graphBoxDraw (dnacpt->aminoButtonBox, BLACK, WHITE) ;
      graphBoxDraw (dnacpt->aminoButtonBox-1, BLACK, LIGHTRED) ;
    }
}

static void dnacptUseKeySet (void)
{ DNACPTGET("dnacptUseKeySet") ;
  if (dnacpt->useKeySet) return ;
  dnacpt->useKeySet = TRUE ;
  if (dnacpt->useKeySetButtonBox)
    {
      graphBoxDraw (dnacpt->useKeySetButtonBox, BLACK, LIGHTRED) ;
      graphBoxDraw (dnacpt->useKeySetButtonBox-1, BLACK, WHITE) ;
    }
}

void dnacptDontUseKeySet (void)	/* called in fmapdisplay */
{ DNACPTGET("dnacptDontUseKeySet") ;
  if (!dnacpt->useKeySet) return ;
  dnacpt->useKeySet = FALSE ;
  if (dnacpt->useKeySetButtonBox)
    {
      graphBoxDraw (dnacpt->useKeySetButtonBox, BLACK, WHITE) ;
      graphBoxDraw (dnacpt->useKeySetButtonBox-1, BLACK, LIGHTRED) ;
    }
}

/*****************************************/

static void dnacptDestroy (void)
{    
  DNACPTGET("dnacptDestroy") ;
  
  dnacpt->magic = 0 ;
  if (dnacpt->dnacptGraph == lastDnacptGraph)
    lastDnacptGraph = 0 ;
  keySetDestroy (dnacpt->cutNone) ;
  keySetDestroy (dnacpt->cutOnce) ;
  keySetDestroy (dnacpt->cutTwice) ;
  keySetDestroy (dnacpt->cutThrice) ;
  keySetDestroy (dnacpt->cutMore) ;
  keySetDestroy (dnacpt->enz) ;
  messfree (dnacpt) ;
}

/*******************************************/

extern int ksetClassComplete(char *text, int len, int classe) ;
static int restrictionCompletion(char *cp, int len)
{
  return ksetClassComplete(cp, len, MYMOTIF) ;
}

/**********************************************************/

static void dnacptPick (int k, double x)
{    
  DNACPTGET("dnacptPick") ;
  
  if (!k)
    return ;

  else  if (k == dnacpt->restrictionBox)
    graphCompletionEntry (restrictionCompletion, dnacpt->restriction,0,0,0,0) ;

  graphActivate(dnacpt->graph) ;
}

/**********************************************************/

static void dnacptAllEnz (void)
{    
  KEYSET ks = query(0,"Find Motif  (Restriction OR Restriction_enzyme)") ;
  DNACPTGET("dnacptAllEnz") ;

  if (keySetMax(ks))
   {
     graphEntryDisable () ;
     keySetNewDisplay (ks, "All restriction enzymes") ;
     dnacptSelectEnz () ;
   }     
}

/**********************************************************/

static void dnacptNrEnz (void)
{    
  KEYSET ks = query(0,"Find Motif (Restriction OR Restriction_enzyme) AND NOT redundant") ;
  DNACPTGET("dnacptNrEnz") ;

  if (!keySetMax(ks))
    { messout("Many different enzymes are equivalent (isishizomeres) "
	      "You may want to define a non redundant subset.\n "
	      "To do this, use the edit button of keyset window  to add the tag "
	      "redundant to some retriction enzymes") ; ;
    dnacptAllEnz() ;
    }
  else
   { 
     graphEntryDisable () ;
     keySetNewDisplay (ks, "Non_redundant  restriction enzymes") ; 
     dnacptSelectEnz () ;
   }
}

/**********************************************************/

static void dnacptLabEnz (void)
{    
  KEYSET ks = query(0,"Find KeySet Lab_Enzymes ; NEIGHBOURS") ;
  DNACPTGET("dnacptLabEnz") ;

  if (!keySetMax(ks))
    { messout("You do not have a keyset \"Lab_enzymes\" defined. "
	      "Select some enzymes and Use the Save_as menu entry "
	      "from the keyset window") ;
    dnacptAllEnz() ;
    }
  else
   {
     graphEntryDisable () ;
     keySetNewDisplay (ks, "Lab_enzymes") ;
     dnacptSelectEnz () ;
   }
}

/**********************************************************/

static void dnacptCopy (void)
{ lastDnacptGraph = 0 ; dnaAnalyse(); }

/**********************************************************/

extern void gelDisplay (void) ;

static void dnacptDumpSegs (void) ;
extern void fMapOspInit (void) ;
extern void fMapBlastInit (void) ;

static void setState0(void)
{
  DNACPTGET("setState") ;
  dnacpt->menuState = 0 ; dnacptDisplay () ;
}

static void setState1(void)
{
  DNACPTGET("setState") ;
  dnacpt->menuState = 1 ; dnacptDisplay () ;
}

static void setState2(void)
{
  DNACPTGET("setState") ;
  dnacpt->menuState = 2 ; dnacptDisplay () ;
}

static void setState3(void)
{
  DNACPTGET("setState") ;
  dnacpt->menuState = 3 ; dnacptDisplay () ;
}

static MENUOPT dnacptMenuBG[] = /* window background */
  {
 { graphDestroy, "Quit" },
 { help, "Help" },
 { graphPrint,"Print" },
 /* give choice of various states */
 { setState0, "Main-Menu"},
 { dnacptDisplay, "Clear" },
 { 0, 0 }
   } ;

static MENUOPT dnacptMenuTrivia[] = 
  {
 { graphDestroy, "Quit" },
 { help, "Help" },
 { graphPrint,"Print" },
 { dnacptCopy, "Copy" },
 { dnacptDisplay, "Clear" },
 { 0, 0 }
   } ;

static MENUOPT dnacptMenu0[] = 
  {
 /* give choice of various states */
 { setState1, "Restriction Analysis"},
 { setState2, "Export"},
 { setState3, "Other Tools"},
 { 0, 0 }
   } ;

static MENUOPT dnacptMenuE[] =  /* external programs */
  {
 { dnaCptFindRepeats, "Repeats" },
 { fMapBlastInit, "Blast" },
 { fMapOspInit, "Osp" },
 { 0, 0 }
   } ;

static MENUOPT dnacptMenu1[] =  /* restriction enzymes */
  {
 { gelDisplay, "Gels" },
 /*  { dnacptRestDistrib, "Motif distribution" }, */
 { setState0, "Main-Menu"}, 
 { 0, 0 }
   } ;

static MENUOPT dnacptMenu10[] =  /* restriction enzymes */
  {
 { setState0, "Main-Menu"}, 
 { 0, 0 }
   } ;

static MENUOPT dnacptMenuEnz[] =  /* restriction enzymes */
  {
 { dnacptNrEnz, "NR Enzymes" },
 { dnacptLabEnz, "Lab Enzymes" },
 { dnacptSelectEnz, "Selected KeySet" },
 { 0, 0 }
   } ;

static MENUOPT dnacptMenu2[] =  /* exports */
  {
 { dnacptFastaDump, "Write/print/mail in Fasta format" },
 { dnacptEmblDump, "Export annotated sequence in EMBL format" },
 { 0, 0 }
   } ;

static MENUOPT dnacptMenuOthers[] =  /* others */
  {
 { dnacptLengths, "Sequence lengths" },
 { dnacptCodonUsage, "Codon usage" },
 /* { dnacptSearchFrancois, "dnaSearch (FML)"}, */
 { dnacptWordUsage, "Word usage" },
 { dnacptEstMatch, "EST/genome match"}, 
 { dnacptSpliceConsensus, "Splice consensus (predicted genes)" },
 { dnacptConfirmedSpliceConsensus, "Splice consensus (transcribed genes)" },
 { dnacptFingerPrint, "Worm Finger print" },
 { dnaAG_TC,"poly R/poly Y" },
 { dnacptMakeTestSequences, "Make pieces" },
 { dnacptDumpSegs, "Dump segs" },
 { 0, 0 }
   } ;

/*************************************************************/

static BOOL MMcheck (int n)
{ if (n < 0 || n > 8) { messout ("Allowed mismatch range [0,8]") ;  return FALSE ;}
  return TRUE ;
}

static int dnacptDisplayTopPart(DNACPT dnacpt, MENUOPT* mm, int line)
{ int box ; 
  float dl = line ;

  graphFitBounds (&mapGraphWidth,&mapGraphHeight) ;

  if (mm)
    { box = graphBoxStart() ; 
    graphButtons (mm, .5, line, mapGraphWidth) ; 
    graphBoxEnd() ; 
    graphBoxDim(box,0,0,0,&dl) ; 
    }
  
  line = dl + 1.2 ; graphLine(0,line, mapGraphWidth, line) ; line += 1.2 ;

  graphText ("Apply to:", 1, line) ;
  box = graphButton ("DNA Window active zone", dnacptDontUseKeySet, 10, line - .5) ;
  if (!dnacpt->useKeySet) graphBoxDraw (box, BLACK, PALEBLUE) ;
  line += 1.2 ;
  box = graphButton ("KeySet of sequences",dnacptUseKeySet, 10, line) ;
  if (dnacpt->useKeySet) graphBoxDraw (box, BLACK, PALEBLUE) ;
  dnacpt->useKeySetButtonBox = box ;

  line += 2.2 ;
  graphLine(0,line, mapGraphWidth, line) ;
  line += 1.2 ;

  return line ;
}

static int dnacptDisplayEnzymes(DNACPT dnacpt, int line)
{ int box ; 
  float dl = 1 ;

  graphText ("Site(s): ", 0.5, line) ; 
  box = graphButton ("DNA", dnacptAminoFalse, 18.5, line - .3) ;
  if (!dnacpt->amino) graphBoxDraw (box, BLACK, PALEBLUE) ;
  box = graphButton ("AA",dnacptAminoTrue, 23.5,line - .3) ;
  if (dnacpt->amino) graphBoxDraw (box, BLACK, PALEBLUE) ;
  dnacpt->aminoButtonBox = box ;
  graphButton ("help", help, 28, line - .3) ;
  /*
    if (dnacpt->complement)
    { box = graphBoxStart() ;
    graphText("Complementary sequence", 35.,line) ;
    graphBoxEnd() ;
    graphBoxDraw(box,BLACK,RED) ;
    }
    */
  
  line += 1.2 ;

  if (!dnacpt->amino) 
    graphText ("ex: hindIII ; atgtnnnat ; wska  (Tab key)", 0.5, line) ; 
  else
    { graphText ("ex: HVRDH ; PG?YQ ; (?: any letter; *: any word)", 0.5, line++) ; 
      graphText ("    Use biomotif for complex patterns", 0.5, line) ;
    }
  line += 1.95 ;
      
  dnacpt->restrictionBox =
    graphCompletionEntry (restrictionCompletion, 
			  dnacpt->restriction, 
			  BOX_LENGTH, 2.5, line, 
			  dnacptRestriction2) ; 
  line += 2.2;
  
  box = graphBoxStart() ;
  graphButtons (dnacptMenuEnz, .5, line, mapGraphWidth) ;
  graphBoxEnd() ;
  graphBoxDim(box,0,0,0,&dl) ; 
  line = dl + 2.2 ;

  graphIntEditor ("Max unknown n", &dnacpt->mmn, 0.5, line, MMcheck) ; line += -1.0 ;
  graphIntEditor ("Max mismatch", &dnacpt->mmch, 0.5, line, MMcheck) ; line += .8 ;

  box = graphButton ("Clear",dnacptDisplay,26.2, line - 1.3) ;
 
  box = graphButton ("Search",dnacptRestriction,26, line + .5) ;
  graphBoxBox(box) ;

  line += 2.9 ;
  graphLine(0,line, mapGraphWidth, line) ;
  line += .8 ;
  

  graphCompletionEntry (restrictionCompletion, dnacpt->restriction,0,0,0,0) ; /* make it acive */

  return line ;
}

/*************************************************************/

static void dnacptDisplay(void)
{ int box ; 
  float dl, line = 1 ;
  DNACPTGET("dnacptDisplay") ;
 
  graphActivate (dnacpt->graph) ;
  graphClear () ;
  graphMenu(dnacptMenuBG) ;   

  graphFitBounds (&mapGraphWidth, &mapGraphHeight) ;

  box = graphBoxStart() ;
  graphButtons (dnacptMenuTrivia, .5, line, mapGraphWidth) ;
  graphBoxEnd() ;
  graphBoxDim(box,0,0,0,&dl) ;
  line = dl + 1.5 ;
      
  switch(dnacpt->menuState)
    {
    case 0:
      graphText("Built in tools", .5, line) ;
      line += 1.5 ;
      box = graphBoxStart() ;
      graphButtons (dnacptMenu0, .5, line, mapGraphWidth) ;
      graphBoxEnd() ;
      graphBoxDim(box,0,0,0,&dl) ;
      line = dl + .5 ;
      
      graphText("External analysis programs", .5, line) ;
      line += 1.5 ;
      box = graphBoxStart() ;
      graphButtons (dnacptMenuE, .5, line, mapGraphWidth) ;
      graphBoxEnd() ;
      graphBoxDim(box,0,0,0,&dl) ;
      line = dl + 2.2 ;
      break ;
    case 1:
      line = dnacptDisplayTopPart(dnacpt, dnacptMenu1, line) ;
      line = dnacptDisplayEnzymes (dnacpt, line) ;
      break ;
    case 2:
      line = dnacptDisplayTopPart(dnacpt, dnacptMenu10, line) ;
      box = graphBoxStart() ;
      graphButtons (dnacptMenu2, .5, line, mapGraphWidth) ;
      graphBoxEnd() ;
      graphBoxDim(box,0,0,0,&dl) ;
      line = dl + 2.2 ;
      break ;
    default:
      line = dnacptDisplayTopPart(dnacpt, dnacptMenu10, line) ;
      box = graphBoxStart() ;
      graphButtons (dnacptMenuOthers, .5, line, mapGraphWidth) ;
      graphBoxEnd() ;
      graphBoxDim(box,0,0,0,&dl) ;
      line = dl + 2.2 ;
      break ;
    }
  line += 2 ;
  dnacpt->line = line ;
  graphTextBounds (80,line) ;
  graphRedraw() ;
  pickRememberDisplaySize("DtDnaTool") ;
}

/****************************************************************/
/****************************************************************/
/********************  public routines   ************************/

void dnaAnalyse(void)
{ 
  DNACPT dnacpt ;

  if(graphActivate(lastDnacptGraph))
    { 
      graphPop() ;
      return ;
    }

  lastDnacptGraph = displayCreate (DtDnaTool) ;
  if (!lastDnacptGraph)
      return ;

  dnacpt=(DNACPT)messalloc(sizeof(struct DNACPTSTUFF));
  dnacpt->magic = &MAGIC_CPT ;
  dnacpt->menuState = 0 ;
  dnacpt->graph = lastDnacptGraph ; /* provision for multi windows */

  graphRegister (DESTROY,dnacptDestroy) ;
  graphRegister (PICK,dnacptPick) ;
  /*  graphRegister (RESIZE, dnacptDisplay) ;: no good, clears the results ! */

  graphAssociate (&MAGIC_CPT, dnacpt);

  strcpy(dnacpt->restriction,"") ;
  dnacpt->amino = FALSE ;

  dnacptDisplay() ;
}

/*************************************************************/
/*************************************************************/
#ifdef JUNK
 This finction is quite cryptic, i comment it out, mieg, may 97			 

typedef struct { KEY key ; int mark , nseg, length ; } RD ;

static Stack rdStack = 0 ;

static Array dnacptGetRD(void)
{ int n = arrayMax(enzSet) ;
  Array rd = arrayCreate(n, RD) ;
  RD *r ; KEY key ; char *cp ; OBJ obj ;
  
  getMatchSeqTag () ;		/* make sure matchSeqTag is set */
 
  rdStack = stackReCreate(rdStack, 50) ;

  while(n--)
    { r = arrayp(rd, n, RD) ;
      key = keySet(enzSet,n) ;
      r->key = key ;
      r->mark = stackMark(rdStack) ;
      if ((obj = bsCreate(key)))
	{ if(bsGetData(obj, matchSeqTag, _Text, &cp))
	    { pushText(rdStack, cp) ;
	      dnaEncodeString(stackText(rdStack,r->mark)) ;
	    }
	  bsDestroy(obj) ;
	}
      r->nseg = r->length = 0 ;
    }
  return rd ;
}
  

static void dnacptShowDistribution(int nseq, int nbases, Array rd)
{ int i = arrayMax(rd) ; RD* r ;
  DNACPTGET("dnacptShowDistrib") ;
  
  graphActivate(dnacpt->graph) ;
  graphText(messprintf("Searched %d sequences, %d bases", nseq, nbases) ,
	    1, dnacpt->line++) ;

  for(i=0; i< arrayMax(rd); i++)
    { 
      r = arrp(rd, i, RD) ;
      graphText(messprintf ("%d segments,", r->nseg), 3,  dnacpt->line) ; 
      graphText(messprintf(" average length %d", (r->nseg ? r->length/r->nseg : 0)), 
		15,  dnacpt->line) ; 
      graphText(messprintf(" expected length %d", 1<< 2*strlen(stackText(rdStack,r->mark))),
		40,  dnacpt->line) ; 
      graphText(messprintf(" Site %s ", dnaDecodeString(stackText(rdStack,r->mark))),
		62,  dnacpt->line++) ; 
    }
  graphTextBounds (80,dnacpt->line += 2) ;
  graphRedraw() ;
}

static void dnacptRestDistrib(void)
{
  Array dna , sites ; KEY seq = 0; void *look ; 
  int i, ii, nseq = 0 , nbases = 0 ; KEYSET kSet ;
  Array rd ;
  RD *r ;
  DNACPTGET("dnacptRestDistrib") ;
      
  if (!graphCheckEditors(dnacpt->graph,0))
    return ;
  if (!dnacpt->useKeySet)
    { messout("I only know how to do that from a keySet") ;
      return ;
    }
  if (!keySetActive(&kSet, &look)) 
    { messout("First select a keySet containing sequences") ;
      return ;
    }

  if (!enzSet || !arrayMax(enzSet))
    { messout ("First select a enzSet") ;
      return ;
    }
  rd = dnacptGetRD() ;
  
  for (ii = 0 ; ii < keySetMax(kSet) ; ii++)
    if ((dna = dnaGet(seq=keySet(kSet,ii))))
      { nbases += arrayMax(dna) ; nseq++ ;
      
      i = keySetMax(rd) ;
	  while(i--)
	    { r = arrp(rd, i, RD) ;
	    
	    sites = arrayReCreate(sites, 20, Site) ;
	    dnacptMatch(dna, 1, sites, stackText(rdStack, r->mark),
			  0, arrayMax(dna), 0, 0, 0) ;
	    if (arrayMax(sites) > 3)
	      { r->length += 
		  arr(sites, arrayMax(sites) - 2, Site).i - 
		      arr(sites, 1, Site).i ;
	      r->nseg += arrayMax(sites) - 3 ;
	      }
	    }
      }

      
  if (!nseq)
    { messout("First select a keySet containing sequences") ;
      return ;
    }

  dnacptShowDistribution(nseq, nbases, rd) ;
}
#endif

/*****************************************/
 
static void dnacptPhrapRepeats(char target)
{
  KEYSET aa = 0, bb = 0;
  void *dummy ;
  FILE *f ;
  char *tmpName  = 0 ;
  static BOOL firstPass = TRUE ;
  Array dna, colors ; 
  KEY seqKey = 0 ;
  void *look = 0 ;
  DNACPT dnacpt = 0 ;
                               
  graphAssFind (&MAGIC_CPT, &dnacpt) ;

  if (dnacpt && dnacpt->useKeySet)
    {
      if (!keySetActive(&aa, &dummy))
	{ messout("First select a keyset containing sequences") ;
	  return ;
	}
      bb = query (aa, "CLASS Sequence ; DNA ") ;
      if (!keySetMax(bb))
	{ messout("First select a keyset containing sequences") ;
	  return ;
	}
    }
  else
    { if (!fMapActive (&dna, &colors, &seqKey, &look))
	{ messout ("First select a DNA window "
		   "or the 'KeySet' button") ;
	  return ;
	}
      dnaRepaint(colors) ;
      
      bb = queryKey (seqKey, ">Subsequence") ;
      if (!keySetMax(bb))
	keySet (bb,0) = seqKey ;
     }

  f = filtmpopen (&tmpName, "w") ;
  if (!f)
    goto abort ;
  dnaDumpFastAKeySet (bb, f, 0) ;
  filclose (f) ;

  if (firstPass)
    { firstPass = FALSE ;
      messout ("wscripts/findrepeats \nContributed by Jeremy Parsons (wash U)") ;
    }
  messStatus ("Looking for Repeats") ;

  callScript("findrepeats", messprintf(" %s &", tmpName)) ; 
#ifdef JUNK
  note that we are not really able to reimport the results
  if (!callScript("findrepeats", tmpName)) ;  /*  &&
      (f = filopen(tmpName,"ace","r")))
    { parseFile(f,0,0) ;
      filremove(tmpName,"ace") ;
    }                            trying to display in ace, not too good */
  filtmpremove(tmpName) ;

  if (look)
    { fMapPleaseRecompute (look) ;  /* set look->pleaseRecompute */
      fMapDraw(look,0) ;
    }
#endif

 abort:
  keySetDestroy(aa) ;
  keySetDestroy(bb) ;
}

void dnaCptFindRepeats(void)
{ dnacptPhrapRepeats('r') ;
}

void dnacptPhrap(void)
{ dnacptPhrapRepeats('p') ;
}

/***************************************************************/

static void dnacptDumpSegs (void)
{ 
  void *dummy = 0 ;
  KEYSET kSet = 0 ;
  FILE *fil = 0 ;
  int level ;
  extern void fMapDumpSegsKeySet (KEYSET) ;
  DNACPTGET("dnacptDumpSegs") ;
      
  if (!dnacpt->useKeySet)
    { messout ("I only know how to do that from a keySet") ;
      return ;
    }
  if (!keySetActive(&kSet, &dummy)) 
    { messout ("First select a keySet containing sequences") ;
      return ;
    }
  if (!(fil = filqueryopen (0, 0, "info", "w", "File to write")))
    return ;
  
  level = freeOutSetFile (fil) ;
  fMapDumpSegsKeySet (kSet) ;
  freeOutClose (level) ;

  filclose (fil) ;
}

/***************************************************************/
/***************************************************************/
#undef DNASEARCH 
#ifdef DNASEARCH

static void dnacptSearchFrancoisDisplay (KEY seqKey, DnaSearch_Result *dr, DNACPT dnacpt, int from, int origin, Array colors)
{
  int i,j; 
  float score = 0 ;
  
  graphText( messprintf("%s: Found %d hits", name (seqKey), dr->occurences),
	     2, dnacpt->line++) ;
  dnacpt->line++ ;
  for (i = 0; dnacpt->line < 1000 && i < dr->itemcount; i++)
    {
      if (dr->items[i].type < 0) 
	{
	  graphText( messprintf("     %5d..%5d  score %g",
				dr->items[i].x1 + from - origin + 1,
				dr->items[i].x2 + from - origin + 1,
				dr->items[i].params.score),
		     2, dnacpt->line++) ;
	  score = dr->items[i].params.score;
	} 

      if (colors && dr->items[i].type >= 0)
	{
	   graphText( messprintf("  detail: %5d..%5d match query %d..%d",
		 dr->items[i].x1,
		 dr->items[i].x2,
		 dr->items[i].params.details.p1,
		 dr->items[i].params.details.p2),
		   2, dnacpt->line++) ;

	  if (dr->items[i].params.details.info>=0)
	    {
	      int index = dr->items[i].params.details.info;
	      char cc, *cp, *cq, *s=dr->items[i].params.details.InfoToString(dr->infos+index);
	      
	      cp = s ; 
	      while (*cp)
		{
		  cq = cp ;
		  while (*cq && *cq != '\n') cq++ ;
		  cc = *cq ; *cq = 0 ;
		  graphText (cp, 2, dnacpt->line++) ;
		  cp = cq ;
		  if (cc) cp++ ;
		}
	      if (s) free(s);
	  }
	}
      
      if (colors)
	for (j = dr->items[i].x1; j <= dr->items[i].x2; j++) 
	  {
	    if (dr->items[i].type == 0)
	      { /* approximate string */
		if (score > -0.1) arr(colors, j + from, char) |= TINT_RED ;
		else arr(colors, j + from, char) |= TINT_CYAN ;
	      }
	    else if (dr->items[i].type < 0) 
	      arr(colors, j + from, char) |= TINT_LIGHTGREEN;
	  }
    }
  dnacpt->line ++ ;
  DnaSearch_PrintResult (stdout, dr, 1) ;
}
 
static void dnacptSearchFrancois(void)
{
  KEY seqKey = 0 ;
  char *cp ;
  static char *motif = 0 ;
  int from, to, origin ;
  Array dna = 0 ;
  Array colors = 0 ;
  KEYSET kSet = 0 ;
  BOOL isFmap = TRUE ;
  void *look = 0 ;
#if !defined(ALPHAJUNKH) 
  DNACPTGET("dnacptSearchFrancois") ;
      
  if(!dnacpt->useKeySet)
    { 
      if(!fMapActive(&dna,&colors,&seqKey,&look))
	{ 
	  messout("First select a dna window") ;
	  return ;
	}
      isFmap = TRUE ; 
      fMapFindZone (look, &from, &to, &origin) ;
      if (to > arrayMax (dna))
	to = arrayMax(dna) ;
      if (from < 0)
	from = 0 ;
    }
  else
    { 
      if (!keySetActive (&kSet, &look)) 
	{ messout ("First select a KeySet containing sequences") ;
	return ;
	}
    }

  if (messPrompt ("seq motif ", motif ? motif : "atttagtc", "tz"))
    {
      messfree (motif) ;
      cp = freeword() ;
      motif = strnew(cp,0) ;


      if (dnacpt->useKeySet || (isFmap && dna && from + 1 < to))
	{
	  DnaSearch_Result *dr = 0 ;
	  DnaSearch_Query dq = DnaSearch_CompileQuery(motif, 1) ;
	  if (dq)
	    {
	      if (dnacpt->useKeySet)
		{ 
		  int ii ;
		  Array kdna = 0 ;
		  KEYSET rKset = keySetCreate () ;
		  
		  for (ii = 0 ; ii < keySetMax(kSet) ; ii++)
		    if ((kdna = dnaGet (keySet (kSet, ii))))
		      {
			DnaSearch_InitSearch (dq, 1) ;
			DnaSearch_SearchBlock (dq, arrp(kdna, 0, char), arrayMax(kdna)) ;
			dr = DnaSearch_TerminateSearch (dq) ;
			if (dr->occurences > 0)
			  {
			    keySetInsert (rKset, keySet (kSet, ii)) ;
			    dnacptSearchFrancoisDisplay (keySet (kSet, ii), dr, dnacpt, 0, 0, 0) ; 
			  }		
			DnaSearch_DestroyResult (dr) ;

			arrayDestroy (kdna) ;
		      }
		  keySetNewDisplay (rKset, motif) ; /* will destroy rKset */
		}
	      else
		{
		  DnaSearch_InitSearch (dq, 1) ;
		  DnaSearch_SearchBlock (dq, arrp(dna, from, char), to - from) ;
		  dr = DnaSearch_TerminateSearch (dq) ;
		  dnaRepaint(colors);
		  dnacptSearchFrancoisDisplay (seqKey, dr, dnacpt, from, origin, colors) ; 
		  DnaSearch_DestroyResult (dr) ;
		}
	      DnaSearch_DestroyQuery (dq) ;
	    }
	  else
	    {
	      int iposbad ;
	      char *posbad = messalloc(dnasearch_parseerrorinfo.last + 2) ;

	      memset (posbad, 32, dnasearch_parseerrorinfo.last) ;
	      for (iposbad = dnasearch_parseerrorinfo.first; iposbad <= dnasearch_parseerrorinfo.last; iposbad++)
		posbad[iposbad] = '^' ;
	      posbad[dnasearch_parseerrorinfo.last + 1] = 0 ;
	      graphText ("stupid query", 2, dnacpt->line++) ;
	      graphText (motif, 2, dnacpt->line++);
	      graphText (posbad, 2, dnacpt->line++);
	      graphText (messprintf ("%s", dnasearch_parseerrorinfo.msg), 2, dnacpt->line++);

	      dnacpt->line++ ;
	      messfree (posbad) ;
	    }
	}
    }
  graphRedraw () ;
  if (!dnacpt->useKeySet)
    fMapReDrawDNA (look) ;
  /* do not arrayDestroy (dna) ; i borrowed it */
#endif /* ALPHA INTEL */
}

#endif  /* FML */

/* find the longest match between est and genome in each est
 *  and report a histogram
 */
static int dnacptFindOneLongestMatchInTg (KEY est, KEY tg)
{
  int i, ii, a1, a2, x1, x2, x02, b, b1, best = 0, covers, from, to, nN = 0 ;
  char *dummy ;
  static Array dna = 0, dnaEst = 0 , err = 0 ;
  OBJ Tg = 0 ;
  static Array units = 0 ;
  static KEY oldEst = 0, oldTg = 0, _af = 0 ;
  BSunit *up ;
  KEY cosmid ;
  A_ERR *ep ;
  BOOL isCosmidReversed = TRUE ;

  if (!_af)
    _af = str2tag ("Assembled_from") ;

  if (!est) /* clean up */
    { 
      arrayDestroy (err) ;
      arrayDestroy (dna) ; 
      arrayDestroy (dnaEst) ;
      arrayDestroy (units) ;
      return 0 ;
    }
  
  if (est != oldEst) /* clean up */
    arrayDestroy (dnaEst) ;
  
  if (tg != oldTg) /* clean up */
    arrayDestroy (dna) ; 
  
  if (!dnaEst && (!dna || arrayMax (dna) > 1))
    { dnaEst = dnaGet (est) ; }
  if (!dnaEst)
    return 0 ;
  
  if (!err)
    err = arrayCreate (32, A_ERR) ;
  if ((Tg = bsCreate (tg)))
    {
      covers = from = to = -1 ;
      if (!dna &&
	  bsGetKey (Tg, str2tag ("Genomic_sequence"), &cosmid) &&
	  bsGetData (Tg, str2tag ("covers"), _Int, &covers) &&
	  bsGetData (Tg, _bsRight, _Text, &dummy) &&
	  bsGetData (Tg, _bsRight, _Int, &from) &&
	  bsGetData (Tg, _bsRight, _Int, &to))
	{
	  dna = dnaGet (cosmid) ; isCosmidReversed = FALSE ;
	  if (from < 1 || from > arrayMax (dna) ||
	      to < 1 || to > arrayMax (dna))
	    arrayMax (dna) = 1 ; /* block study of this gene */
	}
      if (dna && arrayMax (dna) > 1)
	{	
	  units = arrayReCreate (units, 4000, BSunit) ;
	  if (bsGetArray (Tg, _af, units, 5))
	    for (ii = 0 ; ii < arrayMax (units) ; ii += 5)
	      {
		up = arrp (units, ii, BSunit) ;
		if (up[2].k != est)
		  continue ;
		b = 0 ;
		a1 = up[0].i ; a2 = up[1].i ;
		x1 = up[3].i ; x2 = up[4].i ;
		if (from < to)
		  { a1 = from + a1 - 1 ; a2 = from + a2 - 1 ; }
		else
		  { a1 = from - a1 + 1 ; a2 = from - a2 + 1 ; }
		if (a1 >= 1 && a2 <= arrayMax (dna) &&
		    x1 >= 1 && x2 <= arrayMax (dnaEst) &&
		    x2 >= 1 && x2 <= arrayMax (dnaEst))
		  {
		    if (x1 < x2 && a1 < a2)
		      {
			if (isCosmidReversed) 
			  { isCosmidReversed = FALSE ; reverseComplement (dna) ; }
		      }
		    else if (x1 > x2 && a1 > a2)
		      {
			int a0, x0 ;
			x0 = x1 ; x1 = x2 ; x2 = x0 ;
			a0 = a1 ; a1 = a2 ; a2 = a0 ;

			if (isCosmidReversed) 
			  { isCosmidReversed = FALSE ; reverseComplement (dna) ; }
		      }
		    else if (x1 < x2 && a1 > a2)
		      {
			if (! isCosmidReversed) 
			  { isCosmidReversed = TRUE ; reverseComplement (dna) ; }
			a1 = (long)arrayMax (dna) - a1 + 1 ;
			a2 = (long)arrayMax (dna) - a2 + 1 ;
		      }
		    else if (x1 > x2 && a1 < a2)
		      {
			int a0, x0 ;
			x0 = x1 ; x1 = x2 ; x2 = x0 ;
			a0 = a1 ; a1 = a2 ; a2 = a0 ;

			if (! isCosmidReversed) 
			  { isCosmidReversed = TRUE ; reverseComplement (dna) ; }
			a1 = (long)arrayMax (dna) - a1 + 1 ;
			a2 = (long)arrayMax (dna) - a2 + 1 ;
		      }
		    x02 = x2 ;
		    err = aceDnaTrackErrors (dnaEst, x1, &x2, dna, a1, &a2
					     , &nN, err, 6, 6, FALSE, 0, FALSE) ;
		    if (! arrayMax (err))
		      {
			b1 = x1 < x2 ? x2 - x1 + 1 : x1 - x2 + 1 ;
			b = a2 - a1 + 1 ;
			if (b > b1) b = b1 ;
		      }
		    else if (x2 >= x02 - 10) /* ignore cases where tracking failed */
		      {
			for (i = 0, ep = arrp (err, 0, A_ERR) ; i < arrayMax (err) ; i++, ep++)
			  {
			    b1 = ep->iLong - a1 + 1 ;
			    a1 = ep->iLong ;
			    if (b1 > b) b = b1 ;
			  }
			b1 = a2 - a1 + 1 ;
			if (b1 > b) b = b1 ;
		      }
		  }
		if (best < b) best = b ;
	      }      
	}
      bsDestroy (Tg) ;
    }

  /* do not destroy the static arrays */
  if (best && best < 10) invokeDebugger () ;
  return best ;
}/*  dnacptFindOneLongestMatchInTg  */

/*********************************************************************/

static int dnacptFindOneLongestMatch (KEY est)
{
  int i, b, best = 0 ;
  KEY tg ;
  KEYSET genes = 0 ;

  if (!est) /* clean up */
    { 
      dnacptFindOneLongestMatchInTg (0, 0) ;
    }
  else
    {
      genes = queryKey (est, ">From_gene") ;
      for (i = 0 ; i < keySetMax (genes) ; i++)
	{
	  tg = keySet (genes, i) ;
	  b =  dnacptFindOneLongestMatchInTg (est, tg) ;
	  if (b > best) best = b ;
	}
      keySetDestroy (genes) ;
    }

  return best ;
}/*  dnacptFindOneLongestMatch  */

/*********************************************************************/
/* loop per gene, much faster in object access */

void dnacptEstMatch (void)
{
  Array histo = arrayCreate (1000, int) ; /* no est > 32000 */
  int best, ii, jj ;
  KEY est, est2 ;
  KEYSET ests = 0, ks = 0 ;
  BitSet bb = bitSetCreate (lexMax (_VSequence), 0) ;
  void *look ;
  DNACPTGET("dnacptestMatch") ;

  if(!dnacpt->useKeySet)
    { messout("I just know how to do that from a keySet") ;
      return ;
    }
  if (!keySetActive(&ks, &look)) 
    { messout("First select a keySet containing sequences") ;
      return ;
    }

  for (ii = 0 ; ii < keySetMax (ks) ; ii++)
    {
      est = keySet (ks, ii) ;
      if (!bitt (bb, KEYKEY (est)))
	{
	  ests = queryKey (est, ">from_gene ; > read ") ;
	  for (jj = 0 ; jj < keySetMax (ests) ; jj++)
	    {
	      est2 = keySet (ests, jj) ;
	      if (1) est2 = est ;
	      if (!bitt (bb, KEYKEY (est2)))
		{
		  bitSet (bb, KEYKEY (est2)) ;
		  if ((best = dnacptFindOneLongestMatch (est2)))
		    {
		      (array (histo, best, int))++ ;
		      if (1 && best < 30) printf ("%5d %s\n", best, name(est)) ;
		    }
		}
    if (1) break ;
	    }
	  keySetDestroy (ests) ;
	}
    }

  dnacptFindOneLongestMatch (0) ; /* cleanUp */
  plotHisto ("Longest EST/genome match of each EST", histo) ; /* will destroy histo */

  bitSetDestroy (bb) ;

  return ;
}  /* dnacptEstMatch */

/*********************************************************************/
/*********************************************************************/
