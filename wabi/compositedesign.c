#include "ac.h"
#include "cdna.h"
#include "makemrna.h"

/* #define ARRAY_CHECK */
/* #define MALLOC_CHECK */

/**********************************************************************************/
/**********************************************************************************/
/* This strategy is design to handle deep sequencing
 * We construct an oriented graph where each exon is a vertex
 * Each intronis an arc weighted by the number and samples of the supporting tags
 * The introns are ordered by weight
 * For each intron in decreasing order extend it into a path as follows
 *   If the next intron is not yet in a path, extend to the following exon connected by the heaviest intron
 *   If it is attach to the lowest path arriving at this exon from the other side
 *   unless the way i arrive on the exon is much smaller than the way the best path exits the exon
 *   or enough full paths have been constructed
 * Stop when all introns are used or enough paths have been constructed 
 */
typedef struct compositeDesignStruct { Array exons, covering ; BitSet bbIn, bbOut ; KEYSET reads, introns ; int c1, c2 ; AC_HANDLE h ; } DS ;
typedef struct dsVertexStruct { int a1, a2, x1, x2, donor, acceptor, cover, ln, end, path, type, score, nn, clipable, tested ; long unsigned int flag ; KEY est, cDNA_clone, endKey, slKey, SL ; } DSX ;

/**********************************************************************************/

static void showDsFlags (SC *sc, Array aa, long unsigned int flag)
{
  DSX *up ;
  int ii ;
  long unsigned int un = 1 ;
  long unsigned int mask = (un << 32) - 1 ;
  if (aa)
    {
      for (ii = 0, up = arrp (aa, 0, DSX) ; ii < arrayMax (aa) ; ii++, up++)
	{
	  if ((up->flag & flag) != flag)
	    continue ;

	  if (0) if (up->a1 > 6000 || up->a1 < 6300) continue ;
	  
	  fprintf( stderr, "%d:: %d %d\t %d %d\t d=%d a=%d c=%d s=%d e=%d path=%d\t"
		   , ii, up->a1, up->a2, sc->a1 + up->a1 - 1, sc->a1 + up->a2 - 1
		   , up->donor, up->acceptor, up->cover, up->score, up->end, up->path
		   ) ;
	  if (flag) fprintf( stderr, "flags %lu:%lu\t", up->flag & mask, up->flag >> 32) ;
	  if (gGene & up->type) fprintf (stderr, "Gene ") ;
	  if (gGap & up->type) fprintf (stderr, "Gap ") ;
	  if (gLink & up->type) fprintf (stderr, "Link ") ;
	  if (gGhost & up->type) fprintf (stderr, "Ghost ") ;
	  if (gFuseToGhost & up->type) fprintf (stderr, "FuseToGhost ") ;
	  if (gDroppedGene & up->type) fprintf (stderr, "DroppedGene ") ;
	  if (gSuspect & up->type) fprintf (stderr, "Suspect ") ;
	  if (gMicro & up->type) fprintf (stderr, "Micro-") ;
	  if (gX & up->type) fprintf (stderr, "Exon ") ;
	  if (gI & up->type) fprintf (stderr, "Intron ") ;
	  if (gJ & up->type) fprintf (stderr, "Fuzzy ") ;
	  if (gS & up->type) fprintf (stderr, "SL ") ;
	  if (gS0 & up->type) fprintf (stderr, "SL0 ") ;
	  if (gReal5p & up->type) fprintf (stderr, "r5p ") ;
	  if (g5 & up->type) fprintf (stderr, "X5 ") ;
	  if (g3 & up->type) fprintf (stderr, "X3 ") ;
	  if (gD & up->type) fprintf (stderr, "Debut ") ;
	  if (gDF & up->type) fprintf (stderr, "DebutFuzzy ") ;
	  if (gCompleteCDS & up->type) fprintf (stderr, "CDS ") ;
	  if (gF & up->type) fprintf (stderr, "Fin ") ;
	  if (gFF & up->type) fprintf (stderr, "FinFuzzy ") ;
	  if (gReal3p & up->type) fprintf (stderr, "r3p ") ;
	  if (gA & up->type) fprintf (stderr, "polyA ") ;
	  if (gB & up->type) fprintf (stderr, "Alter ") ;
	  if (gStolen & up->type) fprintf (stderr, "Stolen ") ;
	  if (gPredicted & up->type) fprintf (stderr, "Predicted ") ;

	  if (up->clipable) fprintf (stderr, " *") ;
	  
	  fprintf (stderr, "\n") ;
	}
    }
  return ;
} /* showDsFlags */

static void showDs (SC *sc, Array aa)
{
  return showDsFlags (sc, aa, 0) ;
} /* showDs */

/**********************************************************************************/

static int dsScoreOrder (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  n = a->score - b->score ; if (n) return -n ; /* large scores first */
  n = a->type - b->type ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  
  return 0 ;
} /* dsScoreOrder */

/**********************************************************************************/

static int dsEstOrder (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  n = a->est - b->est ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  n = a->type - b->type ; if (n) return n ;
  
  return 0 ;
} /* dsScoreOrder */

/**********************************************************************************/

static int sFlagsOrder (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  if (a->flag > b->flag)  return -1 ; /* large scores first */
  if (a->flag < b->flag)  return 1 ; /* large scores first */
  n = a->a1 - b->a1 ; if (n) return n ;
  
  return 0 ;
} /* sFlagsOrder */

/**********************************************************************************/

static int dsCoverOrder (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  n = a->cover - b->cover ; if (n) return -n ; /* large scores first */
  n = a->type - b->type ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  
  return 0 ;
} /* dsCoverOrder */

/**********************************************************************************/

static int dsEndOrder (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  n = a->end - b->end ; if (n) return -n ; /* large cover first */
  n = a->cover - b->cover ; if (n) return -n ; /* large cover first */
  n = a->type - b->type ; if (n) return n ;
  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  
  return 0 ;
} /* dsEndOrder */

/**********************************************************************************/

static int dsA1Order (const void *va, const void *vb)
{
  DSX *a = (DSX *)va, *b = (DSX *)vb ;
  int n ;

  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  n = a->type - b->type ; if (n) return n ;
  
  return 0 ;
} /* dsA1rder */

/**********************************************************************************/

static int hitA1Order (const void *va, const void *vb)
{
  HIT *a = (HIT *)va, *b = (HIT *)vb ;
  int n ;

  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  n = a->est - b->est ; if (n) return n ;
  
  return 0 ;
} /* dsA1rder */

/**********************************************************************************/

static int smrnaA1Order (const void *va, const void *vb)
{
  SMRNA *a = (SMRNA *)va, *b = (SMRNA *)vb ;
  int n ;

  n = a->a1 - b->a1 ; if (n) return n ;
  n = a->a2 - b->a2 ; if (n) return n ;
  
  return 0 ;
} /* smrnaA1rder */

/**********************************************************************************/

static int ssHappyFew (Array ss)
{
  int iMax = arrayMax (ss) ;
  if (iMax)
     {
       int ii, jj ;
       DSX *ssp, *ssp1 ;

       arraySort (ss, dsA1Order) ;
       for (ii = jj = 0, ssp1 = ssp = arrp (ss, 0, DSX) ; ii < iMax ; ii++, ssp++)
	 {
	   if (ssp->cover > 0)
	     {
	       if (jj < ii)
		 *ssp1 = *ssp ;
	       jj++ ; ssp1++ ;
	     }
	 }
       iMax = arrayMax (ss) = jj ;
     }

  arraySort (ss, dsA1Order) ;
  return iMax ;
} /* ssHappyFew */

/**********************************************************************************/
/* Rationalize the exon flags */
static void mrnaDesignCleanExons (SC *sc, Array ss)
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL debug = FALSE ;
  int ii ;
  DSX *up ;
  int iiMax = ssHappyFew (ss) ;


  /* create and boost some scores */
  for (ii = 0, up = arrp (ss, 0, DSX) ; ii < iiMax ; ii++, up++)
    {
      up->nn = ii ;
      up->score = up->cover ;
      if (up->type & gI) up->score *= 3 ;
      if (up->type & gA) up->score *= 3 ;
    }

  if (iiMax == 1)   
    {       /* eliminate weak short cloud exons */
      up = arrayp (ss,0, DSX) ;
      int da = up->a2 > up->a1 ? up->a2 - up->a1 : up->a1 - up->a2 ;
      if (da < 50 && up->score < 10)
	iiMax = arrayMax (ss) = 0 ;
    }

  /* if we have   2 introns with same acceptor
   * the exon bit going to the second donor
   * should have the support of the short intron
   */
  if (iiMax)
    for (ii = 0, up = arrp (ss, 0, DSX) ; ii < iiMax ; ii++, up++)
      {
	if (up->a2 < up->a1)
	  invokeDebugger () ;
	if (up->type & gI)
	  {
	    DSX *vp ;
	    int jj, ok2 = 0 ;
	    
	    for (jj = ii + 1, vp = up + 1 ; jj < iiMax && vp->a2 < up->a2 && vp->a1 == vp[-1].a2 + 1   ; jj++, vp++)
	      if (vp->type & gI)
		{ ok2 = TRUE ; break ; }
	    if (ok2 && vp->a1 < up->a1 + 30)
	      for (int i = 1 ; ii + i < jj ; i++)
		{
		  if (up[i].score < vp->score)
		    up[i].score = vp->score ;
		}
	  }
      }

  /* polyA are saved */
  if (0 && iiMax)
    for (ii = iiMax - 1, up = arrp (ss, ii, DSX) ; ii >= 0 ; ii--, up--)
      {
	if (up->type & gA)
	  { 
	    DSX *vp ;
	    int i, jj, bestScore = up->score ;
	    int b2 = up->a1 - 1 ;
	    for (i = ii - 1, vp = arrp (ss, i, DSX) ; i >= 0 ; vp--, i--)
	      {
		if (vp->type & gI)
		  {
		    if (vp->a1 == b2 + 1)
		      break ;
		    continue ;
		  }
		if (vp->a2 != b2) break ;
		b2 = vp->a1 - 1 ;
		if (vp->score > bestScore) bestScore = vp->score ;
	      }
	    jj = i ;
	    for (i = jj + 1, vp = arrp (ss, i, DSX) ; i <= ii ; i++, vp++)
	      if (vp->type & gX)
		vp->score = bestScore ;
	  }
      }  
  
  /* overlapping elements are killed if below 1% */
  iiMax = ssHappyFew (ss) ;
  if (iiMax)
    {
      int rho = 100 ; /* rejection rate */
      /* exact overlaps of exons introns */
      for (ii = 0, up = arrp (ss, 0, DSX) ; ii < iiMax - 1 ; ii++, up++)
	{
	  DSX *vp = up + 1 ;
	  if (vp->a1 == up->a1 && vp->a2 == up->a2)
	    {
	      if (rho * vp->cover < up->cover)
		vp->cover = vp->score = 0 ;
	      if (rho * up->cover < vp->cover)
		up->cover = up->score = 0 ;
	    }
	}

      iiMax = ssHappyFew (ss) ;
      /* generalized cassettes */
      for (ii = 0, up = arrp (ss, 0, DSX) ; ii < iiMax - 1 ; ii++, up++)
	{
	  int jj, kk ;
	  DSX *vp, *wp ;
	  for (jj = ii + 1, vp = up + 1 ; jj < iiMax && vp->a1 == up->a1 ; jj++, vp++)
	    if (up->cover && vp->cover)
	      { /* vp is necesarilly longer than up try to extend up */
		int b2 = up->a2 ;
		int cover2 = up->cover ;
		int ns = 0 ;
		KEYSET cs = 0 ;
		for (kk = jj + 1, wp = vp + 1 ; kk  < iiMax && wp->a1 < b2 + 2 ; kk++, wp++)
		  {
		    if (wp->cover && wp->a1 == b2 + 1 && wp->a2 <= vp->a2)
		      {
			if (wp->cover > cover2 && ! (wp->type & (gReal5p | gReal3p)))
			  cover2 = wp->cover ;
			b2 = wp->a2 ;
			if (! cs) cs = keySetCreate () ;
			keySet (cs, ns++) = kk ;
		      }
		  }
		if (b2 == vp->a2)
		  { /* cassette found */
		    if (cover2 > rho * vp->cover)
		      vp->cover = vp->score = 0 ;
		    if (rho * cover2 < vp->cover)
		      {
			up->cover = up->score = 0 ;
			for (int m = 0 ; m < keySetMax (cs) ; m++)
			  {
			    wp = arrp (ss, keySet(cs, m), DSX) ;
			    wp->cover = wp->score = 0 ;
			  }
		      }
		  }
		keySetDestroy (cs) ;
	      }
	}
    }

  iiMax = ssHappyFew (ss) ;
  /* number the segments */
  if (iiMax)
    {
      int a1, a2 ;

      a1 = a2 = arr (ss, 0, DSX).a1 ;
      for (ii = 0, up = arrp (ss, 0, DSX) ; ii < iiMax ; ii++, up++)
	{
	  up->nn = ii ;
	  a2 = up->a2 ;
	}
      if (a2 < a1 + 20)
	arrayMax (ss) = 0 ;
    }
  if (debug)
    {
      fprintf (stderr, "mrnaDesignCleanExons\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignCleanExons done\n") ;
    }
  ac_free (h) ;
   return ;
} /* mrnaDesignCleanExons */

/**********************************************************************************/
/* get the X Composite and their counts */
static BOOL mrnaDesignGetGraphElements (DS *ds, S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  BOOL debug = TRUE ; /* tyty TRUE ; */
  HIT *up, *up1 ;
  DSX *ssp ;
  int i1, ii, jj, iMax = 0, ir, iss = 0 ;
  int nXI = 0, nXE = 0, nXA = 0, nXSL = 0, nXends = 0 ;
  BOOL isDown ;
  Array units ;
  BSunit *uu ;
  OBJ Est = 0 ;
  static int nnnn = 0 ;
  nnnn++ ;
  AC_HANDLE h = ac_new_handle () ;
  Array ss = arrayHandleCreate (256, DSX, ds->h) ;
  KEYSET ks1 = keySetHandleCreate (h) ;
  KEYSET ks2 = keySetHandleCreate (h) ;
  KEYSET ks5 = keySetHandleCreate (h) ;
  KEYSET ks3 = keySetHandleCreate (h) ;
  int iiMax = arrayMax ( s2m->plainHits) ;
  int ii1Max = arrayMax ( gmrna->hits) ;
  BOOL hasGoody = FALSE ;
  int aMax = sc->a2 - sc->a1 ;
  if (aMax < 0) aMax = - aMax ;
  /* grab the introns in the plainHits EST alignments */
  isDown = (sc->a1 < sc->a2) ? TRUE : FALSE ;
  units = arrayHandleCreate (12, BSunit, h) ;

  for (ii = 0, up = arrp (s2m->plainHits, ii, HIT) ; ii < iiMax ; ii++, up++)
    for (i1 = 0, up1 = arrp (gmrna->hits, 0, HIT) ; i1 < ii1Max ; i1++, up1++)
      {
	int b1, b2 ;

	if (up->a1 > up->a2 - 10) continue ;
	if (sc->a1 < sc->a2)
	  {
	    b1 = up1->a1 + sc->a1 - 1 ;
	    b2 = up1->a2 + sc->a1 - 1 ;
	    if (b2 < 0 || b1 > sc->a2) continue ;
	  }
	else
	  {
	    b2 = - up1->a1 + sc->a1 + 1 ;
	    b1 = - up1->a2 + sc->a1 + 1 ;
	    if (b2 < 0 || b1 > sc->a1) continue ;
	  }

	
	if (
	    ((up1->type & gI) && (up->type & gI) && up->a1 == b1 && up->a2 == b2) ||
	    ((up1->type & gX) && (up->type & gX) && up->a1 < b2 && up->a2 > b1)
	    )
	    {
	      BOOL flipped = up->x1 > up->x2 ? TRUE : FALSE ;
	      int a1 = isDown ? up->a1 - sc->a1 + 1 : sc->a1 - up->a2 + 1 ;
	      int a2 = isDown ? up->a2 - sc->a1 + 1 : sc->a1 - up->a1 + 1 ;
	      int xdx = up->x2 - up->x1 ;
	      int cdx = up->clipEnd - up->clipTop ;
	      if (cdx < 0) cdx = - cdx ;
	      if (xdx < 0) xdx = - xdx ;
	      
	      if (0 && debug)
		{
		  fprintf (stderr, "......getSS ii=%d a1=%d a2=%d %s \n", ii, a1, a2, name (up->est)) ;
		  showDs (sc, ss) ;
		  fprintf (stderr, "......getSS done\n") ;
		}
	      if ( ! strncmp (name(up->est), "XI_", 3) ||
		   ! strncmp (name(up->est), "XW_", 3) ||
		   ! strncmp (name(up->est), "XY_", 3)
		   )
		{
		  if (! keyFindTag (up->est, _Composite))
		    continue ;
		  KEYSET ks = queryKey (up->est, "(IS XY_* || IS XW_*) && (ct_ac || other || small_deletion) && (ct_ac > 1 || other > 1 || (ct_ac && other)) && ! composite > 50") ;
		  if (keySetMax (ks))
		    {
		      for (int i = 0 ; i < keySetMax (ks) ; i++)
			{
			  KEY est = keySet (ks, i) ;
			  OBJ Est = bsUpdate (est) ;
			  if (Est)
			    {
			      if (bsFindTag (Est, str2tag ("Is_read")))
				bsRemove (Est) ;
			    }
			  bsSave (Est) ;
			}
		      keySetDestroy (ks) ;
		      continue ;
		    }
		  keySetDestroy (ks) ;
		  if (up->type & gX)
		    { /* reclip the exon borders */
		      int dx = up->x2 - up->x1 - 18 ;
		      if (dx > 0)
			{
			  if (up->x1 == 1)
			    {
			      up->x1 += dx ;
			      if (up->a1 < up->a2) up->a1 += dx ;
			      else up->a1 -= dx ;
			    }
			  else if (up[1].est != up->est)
			    {
			      up->x2 -= dx ;
			      if (up->a1 < up->a2) up->a2 -= dx ;
			      else up->a2 += dx ;
			    }
			  a1 = isDown ? up->a1 - sc->a1 + 1 : sc->a1 - up->a2 + 1 ;
			  a2 = isDown ? up->a2 - sc->a1 + 1 : sc->a1 - up->a1 + 1 ;
			}
		    }
		    
		  if ((Est = bsCreate (up->est)))
		    {
		      int s = 0 ;
		      KEY intron = 0 ;
		      if (bsGetArray (Est, _Composite, units, 1))
			for (int jj = 0 ; jj < 1 && jj < arrayMax (units) ; jj+=1)	
			  {
			    uu = arrp (units, jj, BSunit) ;
			    s =  uu[0].i ; 
			  }
		      if (s && bsGetArray (Est, _Intron, units, 1))
			for (int jj = 0 ; jj < 1 && jj < arrayMax (units) ; jj+=1)	
			  {
			    uu = arrp (units, jj, BSunit) ;
			    intron = uu[0].k ;
			    if (keyFindTag (intron, _Other))
			      s /= 5 ;
			    if (0 && ! strncmp (name(up->est), "XW_", 3))
			      s *= 3 ; /* il y a plus de bruit car on n'a pas selectionne */
			  }
		      
		      if (s > 0 && a1 >= 1 && a2 < aMax)
			{
			  DSX *ssp = arrayp (ss, iss++, DSX) ;
			  ssp->a1 = a1 ;
			  ssp->a2 = a2 ;
			  hasGoody = TRUE ;
			  s *= 1 ;
			  if (up->type & gI)
			    { nXI++; ssp->type = gI ; ssp->cover = s ; }
			  else if (up->type & gX)
			    { nXE++; ssp->type = gX ; ssp->cover = s ; ssp->clipable = 2 ; }
			}
		      
		      bsDestroy (Est) ;
		    }
		  continue ;
		}
	      
	      /* register the exon support and  shrink the higher expressed segments */
	      if ( ! strncmp (name(up->est), "XG_", 3) ||
		   ! strncmp (name(up->est), "XH_", 3)
		   )
		{
		  int iss0 = iss, b1, b2 ;
		  if (10 * xdx < 9 * cdx) continue ;
		  if (! keyFindTag (up->est, _Composite))
		    continue ;
		  if (! strcmp (name(up->est), "XG_U00096.3__179215_176175"))
		    invokeDebugger () ;
		  if ((Est = bsCreate (up->est)))
		    {
		      if (bsGetArray (Est, _Composite, units, 3) && arrayMax (units) >= 3)
			{
			  if (up->reverse == flipped)
			    for (ir = 0 ; ir < arrayMax (units) ; ir += 3)
			      {
				uu = arrp (units, ir, BSunit) ;
				b1 = a1 + uu[0].i - 1 ; 	
				b2 = a1 + uu[1].i - 1 ;

				if (b1 < 1) b1 = 1 ;
				if (b2 >= aMax) b2 = aMax - 1 ;
				if (b1 > b2 - 10) continue ;
				ssp = arrayp (ss, iss++, DSX) ;
				ssp->type = gX ;
				ssp->a1 = b1 ; ssp->a2 = b2 ; ssp->cover = uu[2].i ;
				if (ssp->cover > 100)
				  hasGoody = TRUE ;
				if (ir == 0 || ir == arrayMax (units) - 3)
				  ssp->clipable = 1 ;
				if (0 && debug)
				  fprintf(stderr, "%s %d %d %d %s %s \n", name(up->est), ssp->a1, ssp->a2, ssp->cover, up->reverse ? "reverse" : "forward", isDown ? "Down":  " Up") ;
			      }
			  else
			    for (ir = arrayMax (units) - 3 ; ir >= 0 ; ir -= 3)
			      {
				uu = arrp (units, ir, BSunit) ;
				b1 = a2 - uu[1].i + 1 ; 	
				b2 = a2 - uu[0].i + 1 ;
				
				if (b1 < 1) b1 = 1 ;
				if (b2 >= aMax) b2 = aMax - 1 ;
				if (b1 > b2 - 10) continue ;
				ssp = arrayp (ss, iss++, DSX) ;
				ssp->a1 = b1 ; ssp->a2 = b2 ; ssp->cover = uu[2].i ;
				if (ssp->cover > 100)
				  hasGoody = TRUE ;
				if (ir == 0 || ir == arrayMax (units) - 3)
				  ssp->clipable = 1 ;
				if (0 && debug)
				  fprintf(stderr, "----%s %d %d %d %s %s \n", name(up->est), ssp->a1, ssp->a2, ssp->cover, up->reverse ? "reverse" : "forward", isDown ? "Down":  " Up") ;
			      }
			}
		      
		      bsDestroy (Est) ;
		      
		      /* shrink the higher expressed segments */
		      if ( ! strncmp (name(up->est), "XG_", 3) || ! strncmp (name(up->est), "XH_", 3))
			{
			  int ii ;
			  for (ii = iss0 ; ii < iss - 1 ; ii++)
			    {
			      DSX *ssq ;
			      ssp = arrayp (ss, ii, DSX) ;
			      ssq = ssp + 1 ;
			      if (ssp->a1 < ssp->a2 && ssp->a2 == ssq->a1 - 1)
				{
				  if (ssp->cover < ssq->cover && ssq->a2 > ssq->a1 + 10)
				    { ssp->a2 += 10 ; ssq->a1 += 10 ; }
				  if (ssp->cover > ssq->cover && ssp->a2 > ssp->a1 + 10)
				    { ssp->a2 -= 10 ; ssq->a1 -= 10 ; }
				}
			      if (ssp->a1 > ssp->a2 && ssp->a2 == ssq->a1 + 1)
				{
				  messcrash ("a1>a2 should never happen") ;
				  if (ssp->cover < ssq->cover && ssq->a1 > ssp->a1 + 10)
				    { ssp->a2 -= 10 ; ssq->a1 -= 10 ; }
				  if (ssp->cover > ssq->cover && ssq->a2 > ssq->a1 + 10)
				    { ssp->a2 += 10 ; ssq->a1 += 10 ; }
				}
			    }
			  if (0)
			    {
			      ssp = arrayp (ss, iss0, DSX) ;
			      ssp->a1 += (ssp->a1 < ssp->a2 ? +10 : -10) ;
			      ssp = arrayp (ss, iss - 1, DSX) ;
			      ssp->a2 += (ssp->a1 < ssp->a2 ? -10 : +10) ;
			      if (ssp->a1 > ssp->a2) ssp->cover = 0 ;
			    }
			  for (ii = iss0 ; ii < iss ; ii++)
			    {
			      ssp = arrayp (ss, ii, DSX) ;
			      if (0 && debug)
				fprintf(stderr, "... %s %d %d %d \n", name(up->est), ssp->a1, ssp->a2, ssp->cover) ;
			    }
			}
		      continue ;
		    }
		}
	      
	      /* register the predicted cds support */
	      if (! strncmp (name(up->est), "Xcds_", 5))
		{
		  if (! keyFindTag (up->est, _Composite))
		    continue ;
		  if (10 * xdx < 9 * cdx) continue ;
		  if (a1 < 1) a1 = 1 ;
		  if (a2 >= aMax) a2 = aMax - 1 ;
		  if (a1 < a2 - 10)
		    {
		      if ((Est = bsCreate (up->est)))
			{
			  int s = 12 ;
			  nXE++ ;
			  ssp = arrayp (ss, iss++, DSX) ;
			  ssp->a1 = a1 ;
			  ssp->a2 = a2 ; 
			  ssp->cover = s ;
			  ssp->type = gX | gCompleteCDS  | gFF ; /* was | gFF */
			  ssp->clipable = 1 ;
			}
		      bsDestroy (Est) ;
		    }
		  continue ;
		}
	      
	      /* register the polyA */
	      if (! strncmp (name(up->est), "XA_", 3) ||
		  ! strncmp (name(up->est), "XSL", 3) ||
		  ! strncmp (name(up->est), "Xends_", 6)
		  )
		{
		  if (! keyFindTag (up->est, _Composite))
		    continue ;
		  if (10 * xdx < 9 * cdx) continue ;
		  if ((Est = bsCreate (up->est)))
		    {
		      int s = 0 ;
		      if (bsGetArray (Est, _Composite, units, 1))
			for (int jj = 0 ; jj < 1 && jj < arrayMax (units) ; jj+=1)	
			  {
			    uu = arrp (units, jj, BSunit) ;
			    s =  uu[0].i ; 
			  }
		      if (s)
			{
			  if (! strncmp (name(up->est), "Xends_ELF_", 9) && ! isDown)
			    continue ;
			  if (! strncmp (name(up->est), "Xends_ERF_", 9) && ! isDown)
			    continue ;
			  if (! strncmp (name(up->est), "Xends_ERR_", 9) &&  isDown)
			    continue ;
			  if (! strncmp (name(up->est), "Xends_ELR_", 9) &&  isDown)
			    continue ;
			  
			  hasGoody = TRUE ;
			  BOOL clipable = FALSE ;
			  int type = 0 ;
			  KEY slKey = 0 ;
			  
			  if (! strncmp (name(up->est), "XA_", 3))
			    {
			      nXA++ ;
			      a1 = a2 ;
			      type = gX | gA | g3 ;
			    }
			  else if (! strncmp (name(up->est), "XSL", 3))
			    {
			      nXSL++ ;
			      slKey = up->est ;
			      a2 = a1 ;
			      type = gX | gS | g5 | gReal5p ;
			    }
			  else
			    {
			      clipable = TRUE ;
			      
			      if (! strncmp (name(up->est), "Xends_ELF.", 10) && isDown) 
				{ nXends++ ; type = (gReal5p | gDF | gX) ; if (a2 > a1 + 30) { a1 += 0 ; a2 -= 0*30 ; }}
			      else if (! strncmp (name(up->est), "Xends_ERF.", 10) && isDown) 
				{ nXends++ ; type = (gReal3p | gFF | gX) ; if (a2 < a1 - 30) { a2 -= 0 ; a2 += 0*30 ; }} 
			      else if (! strncmp (name(up->est), "Xends_ERR.", 10) && ! isDown) 
				{ nXends++ ; type = (gReal5p | gDF | gX) ; if (a2 < a1 - 30) { a1 += 0 ; a2 += 0*30 ; }} 
			      else if (! strncmp (name(up->est), "Xends_ELR.", 10) && ! isDown) 
				{ nXends++ ; type = (gReal3p | gFF | gX) ; if (a2 > a1 + 30) { a1 -= 0 ; a2 -= 0*30 ; }}
			    }

			  if (a1 < a2)
			    {
			      if (a1 < 1) a1 = 1 ;
			      if (a2 >= aMax) a2 = aMax - 1 ;
			    }
			  if (a1 < a2)
			    {
			      ssp = arrayp (ss, iss++, DSX) ;
			      ssp->a1 = a1 ;
			      ssp->a2 = a2 ;
			      ssp->cover = s ;
			      ssp->end = s ;
			      ssp->clipable = clipable ;
			      ssp->type = type ;
			      ssp->slKey = slKey ;
			      ssp->endKey = up->est ;
			      ssp->cDNA_clone = up->cDNA_clone ;
			    }
			}
		      bsDestroy (Est) ;
		    }
		  continue ;
		}
	    }
      }

  if (debug)
    {
       fprintf (stderr, "mrnaDesignGetGraphElements\n") ;
       showDs (sc, ss) ;
       fprintf (stderr, "mrnaDesignGetGraphElements  done\n") ;
    }
  if (! hasGoody)
    { ac_free (h) ; return FALSE ;}
  arraySort (ss, dsA1Order) ;
  arrayCompress (ss) ;
  if (debug)
    {
       fprintf (stderr, "mrnaDesignGetGraphElements\n") ;
       showDs (sc, ss) ;
       fprintf (stderr, "mrnaDesignGetGraphElements  done\n") ;
       fprintf (stderr, "s2m->plainHits\n") ;
       showHits (s2m->plainHits) ;
       fprintf (stderr, "s2m->plainHits  done\n") ;
       fprintf (stderr, "gmrna->hits\n") ;
       showHits (gmrna->hits) ;
       fprintf (stderr, "gmrna->hits  done\n") ;
       fprintf (stderr, "ss\n") ;
       showDs (sc, ss) ;
       fprintf (stderr, "ss  done\n") ;
     }

  iMax = ssHappyFew (ss) ;
  if (iMax) /* kill overlapping introns below 1/1000 */
     {
       for (jj = 0, ssp = arrp (ss, jj, DSX) ; jj < iMax - 1 ; ssp++, jj++)
         {
	   int i ;
	   DSX *ssp1 ;

	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gI)
	     {
	       int b2 = ssp->a2 ;
	       for (i = 1, ssp1 = ssp + i ; i + jj < iMax && ssp1->a1 <= b2 ; i++, ssp1++)
		 {
		   if (ssp1->type & gI)
		     {
		       if (ssp->cover / 1000 > ssp1->cover)
			 ssp1->cover = 0 ;
		       if (ssp1->cover / 1000 > ssp->cover)
			 ssp->cover = 0 ;
		     }
		 }
	     }
	 }
     }
  
  iMax = ssHappyFew (ss) ;
  if (! iMax)
    goto done ;

if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements kill overlapping introns\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements kill overlapping introns done\n") ;
    }
    /* clip end elements overapping an intron */
  for (int ii = 0 ; ii < iMax ; ii++)
    {
      ssp = arrp (ss, ii, DSX) ;
      if (ssp->type & gReal5p)
	{
	  for (int jj = 0 ; jj < iMax ; jj++)
	    {
	      DSX *ssp1 = arrp (ss, jj, DSX) ;
	      if (ssp1->type & gI)
		if (ssp1->a1 > ssp->a1 + 20 && ssp1->a1 < ssp->a2)
		  ssp->a2 = ssp1->a1 - 1 ;
	    }
	}
      if (ssp->type & gReal3p)
	{
	  for (int jj = 0 ; jj < iMax ; jj++)
	    {
	      DSX *ssp1 = arrp (ss, jj, DSX) ;
	      if (ssp1->type & gI)
		if (ssp1->a2 > ssp->a1 && ssp1->a2 < ssp->a2 - 20)
		  ssp->a1 = ssp1->a2 + 1 ;
	    }
	}
    }
    
  /*  register all boundaries */
  keySet (ks1, 0) = 0 ;
  keySet (ks2, 0) = 0 ;
  for (int ii = 0, m = 1, n = 1 ; ii < iMax ; ii++)
    {
      ssp = arrp (ss, ii, DSX) ;
      keySet (ks1, m++) = ssp->a1 ;         /* begin */
      keySet (ks2, n++) = ssp->a2 ;         /* end */
      keySet (ks1, m++) = ssp->a2 + 1 ;     /* virtual begin */
      keySet (ks2, n++) = ssp->a1 - 1 ;     /* virtual end */
      if (ssp->type & gI)
	{
	  keySet (ks5, m++) = ssp->a2 + 1 ;     /* exon acceptor */
	  keySet (ks3, n++) = ssp->a1 - 1 ;     /* exon donor */
	}
    }
  keySetSort (ks1) ;
  keySetSort (ks2) ;
  keySetSort (ks5) ;
  keySetSort (ks3) ;
  keySetCompress (ks1) ;
  keySetCompress (ks2) ;
  keySetCompress (ks5) ;
  keySetCompress (ks3) ;
  
  /* split when needed */
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements before split\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements before split done\n") ;
    }
  if (1)
    {
      int m = 0, n = 0, b1 ;
      int mMax = keySetMax (ks1), nMax = keySetMax (ks2) ;
      Array ss1 = arrayHandleCreate (2*ii + 1 , DSX, ds->h ) ;
      int j5Max = keySetMax (ks5) ;
      int j3Max = keySetMax (ks3) ;
      
      for (int ii = 0, jj = 0 ; ii < iMax ; ii++)
	  {
	    int i4 = 0 ;
	    ssp = arrp (ss, ii, DSX) ;
	    if (ssp->type & gX)
	      {
		while (m > 0 && keySet (ks1, m) >= ssp->a1)
		  m-- ;
		while (n > 0 && keySet (ks2, n) >= ssp->a1)
		  n-- ;
		while (m < mMax  && keySet (ks1, m) < ssp->a1)
		  m++ ;
		while (n < nMax && keySet (ks2, n) < ssp->a1)
		  n++ ;
		b1 = ssp->a1 ;
		while (n < nMax && keySet (ks2, n) <= ssp->a2)
		  {
		    if (1 ||  /* otherwise in 0_923 je cree un trou avec l'exon au dessus */
			i4++ || (! ssp->clipable || keySet (ks2, n)  > b1 + 10))
		      {
			DSX *ssp1 = arrayp (ss1, jj++, DSX) ;
			*ssp1 = *ssp ;
			ssp1->a1 = b1 ; ssp1->a2 = keySet (ks2, n) ;
		      }
		    b1 = keySet (ks2, n)  + 1 ;
		    n++ ;
		  }
		if (n == nMax) n-- ;
		if (b1 < ssp->a2 && (! ssp->clipable || b1 < ssp->a2 - 10))
		  {
		    DSX *ssp1 = arrayp (ss1, jj++, DSX) ;
		    *ssp1 = *ssp ;
		    ssp1->a1 = b1 ; ssp1->a2 = ssp->a2 ;
		  }
	      }
	    else if (ssp->type & gI)
	      {
		DSX *ssp1 = arrayp (ss1, jj++, DSX) ;
		*ssp1 = *ssp ;
	      }	    
	  }
      arrayDestroy (ss) ;
      ss = ss1 ;
      iMax = arrayMax (ss) ;
      
      for (int ii = 0, j3 = 0, j5 = 0 ; ii < iMax ; ii++)
	{
	  ssp = arrayp (ss1, ii, DSX) ;
	  if (ssp->type & gX)
	    {
	      int b1 = ssp->a1 ;
	      int b2 = ssp->a2 ;
	      
	      while (j5 > 0 && keySet (ks5, j5) > b1) j5-- ;
	      while (j5 < j5Max - 1 && keySet (ks5, j5) < b1) j5++ ;
	      if (keySet (ks5, j5) == b1) ssp->type |= g5 ;
	      
	      while (j3 > 0 && keySet (ks3, j3) > b2) j3-- ;
	      while (j3 < j3Max - 1 && keySet (ks3, j3) < b2) j3++ ;
	      if (keySet (ks3, j3) == b2) ssp->type |= g3 ;
	    }
	}
    }

  if (0 && debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after split\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after split done\n") ;
    }
	      
  arraySort (ss, dsA1Order) ;
	    
  if (0 && debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after split sort\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after split sort done\n") ;
    }
	      
  /* at each position choose the maximal value */
  if (1)
    {
      int jj = 0, gIX = gI | gX ;

      for (int ii = 0 ; ii < iMax ; ii++)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  
	  int type = ssp->type & gIX ;
	  int b1 = ssp->a1 ;
	  int b2 = ssp->a2 ;
	  DSX *ssp1 ;
	  
	  for (jj = ii + 1, ssp1 = ssp + 1 ; jj < iMax && (ssp1->type & gIX) == type && ssp1->a1 == b1 && ssp1->a2 == b2 ; jj++, ssp1++)
	    {
	      if ((ssp->type & gX) && (ssp1->type & (gDF | gFF)))
		ssp1->cover = ssp->cover ;
	      ssp->type |= ssp1->type ;
	      if (ssp->cover < ssp1->cover)
		ssp->cover = ssp1->cover ;
	      if (ssp->end < ssp1->end)
		{ ssp->end = ssp1->end ; ssp->endKey = ssp1->endKey ; }
	      if (ssp->flag < ssp1->flag)
		ssp->flag = ssp1->flag ;
	      if (! ssp->est) ssp->est = ssp1->est ;
	      if (! ssp->cDNA_clone) ssp->cDNA_clone = ssp1->cDNA_clone ;
	      if (! ssp->slKey) ssp->slKey = ssp1->slKey ;
	      if (! ssp->endKey) ssp->endKey = ssp1->endKey ;
	      ssp1->cover = 0 ;
	      ssp1->flag = 0 ;
	    }
	  ii = jj - 1 ;
	}
    }

  iMax = ssHappyFew (ss) ;
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after top score\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after top score done\n") ;
    }

 
  if (iMax > 1) /* merge continuous exons */
     {
       for (jj = 0, ssp = arrp (ss, jj, DSX) ; jj < iMax - 1 ; ssp++, jj++)
	 {
	   if (ssp->a1 + 12 < ssp->a2 && (ssp->type & g5))
	     {
	       if ((ssp[1].type & gX) && ! (ssp[1].type & g5) && ssp->cover < ssp[1].cover && ssp[1].a1 == ssp->a2 + 1)
		 ssp->cover = ssp[1].cover ; 
	     }
	 }
       for (jj = 1, ssp = arrp (ss, jj, DSX) ; jj < iMax ; ssp++, jj++)
	 {
	   if (ssp->a1 +12 < ssp->a2 && (ssp->type & g3) && ssp[-1].a2 + 1 == ssp->a1)
	     {
	       if ((ssp[-1].type & gX) && ! (ssp[-1].type & g3) && ssp->cover < ssp[-1].cover)
		 ssp->cover = ssp[-1].cover ; 
	     }
	 }
       for (jj = 0, ssp = arrp (ss, jj, DSX) ; jj < iMax - 1 ; ssp++, jj++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if ((ssp->type & gX) && ! (ssp->type & (gS | gA | g5 | g3)))
	     {
	       int i ;
	       int b1 = ssp->a2 + 1 ;
	       DSX *ssp1 ;
	       for (i = jj + 1, ssp1 = arrp (ss, i, DSX) ; ssp1->a1 == b1 && ssp->type == ssp1->type && i < iMax ; ssp1++, i++)
		 {
		   if (ssp->cover < ssp1->cover)
		     ssp->cover = ssp1->cover ;
		   if (ssp->end < ssp1->end)
		     { ssp->end = ssp1->end ; ssp->endKey = ssp1->endKey ; }
		   ssp->a2 = ssp1->a2 ;
		   ssp1->cover = 0 ;
		   b1 = ssp->a2 + 1 ;
		 }
	     }
	 }
     }

  iMax = ssHappyFew (ss) ;
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after merge exons \n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after merge exons done\n") ;
    }

 
  if (iMax) /* absorb in exons 2% of intron support */
     {
       for (jj = 0, ssp = arrp (ss, 0, DSX) ; jj < iMax ; jj++, ssp++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gI)
	     {
	       DSX *ssp1 ;
	       int i, c = ssp->cover/50 ;
	       int b1 = ssp->a1 ;
	       int b2 = ssp->a2 ;
	       if (c > 0)
		 {
		   for (i = -1, ssp1 = ssp - 1 ; ssp1->a1 == b1 && jj + i >= 0 ; ssp1--, i--)
		     if ((ssp1->type & gX) && ! (ssp1->type & (gS | gA | g5 | g3)))
		       ssp1->cover -= c ;
		   for (i = 1, ssp1 = ssp + 1 ; ssp1->a2 <= b2 && jj + i < iMax ; ssp1++, i++)
		     if ((ssp1->type & gX) && ! (ssp1->type & (gS | gA | g5 | g3)))
		       ssp1->cover -= c ;
		 }
	     }
	 }
     }

  iMax = ssHappyFew (ss) ;
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after deduce in exon 3%% of intron score\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after deduce in exon 3%% of intron score done\n") ;
    }


  if (iMax) /* absorb leaking exons at donor or acceptor site */
     {
       for (jj = 0, ssp = arrp (ss, 0, DSX) ; jj < iMax ; jj++, ssp++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gI)
	     {
	       DSX *ssp1 ;
	       int i ;
	       int b1 = ssp->a1 ;
	       int b2 = ssp->a2 ;

	       for (i = -1, ssp1 = ssp - 1 ; ssp1->a1 == b1 && jj + i >= 0 ; ssp1--, i--)
		 if ((ssp1->type & gX) && ssp1->a2 < b1 + 12 && ! (ssp1->type & (gA | g3)))
		   {
		     int k = 1 ;
		     while (ssp1[k].type & gI)
		       k++ ;
		     if ((ssp1[k].type & gX) && ssp1[k].a1 == ssp1->a2 + 1)
		       ssp1->cover = ssp1[k].cover ;
		     else
		       ssp1->cover = 0 ;
		   }
	       for (i = 1, ssp1 = ssp + 1 ; ssp1->a1 <= b2 && jj + i < iMax ; ssp1++, i++)
		 if ((ssp1->type & gX) && ssp1->a1 > b2 - 12 && ! (ssp1->type & (gS | g5) ))
		   {
		     int k = 1 ;
		     while (ssp1[-k].type & gI)
		       k++ ;
		     if ((ssp1[-k].type & gX) && ssp1[-k].a2 == ssp1->a1 - 1)
		       ssp1->cover = ssp1[-k].cover ;
		     else
		       ssp1->cover = 0 ;
		   }
	     }
	 }
     }

  /* keep happy few */
  iMax = ssHappyFew (ss) ;
  
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after absorb leaking exons\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after absorb leaking exons done\n") ;
    }
  
  
  /* regularize the real5p */
  if (1)
    {
      /* SLs cannot differ by 1 base, since they are acceptor sites */
      for (int ii = 0 ; ii < iMax - 1  ; ii++)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  DSX *ssp1 = ssp + 1 ;
	  
	  if (ssp->type & gS)
	    if ((ssp1->type & gS) && (ssp1->a1 == ssp->a1 + 1))
	      {
		if (ssp->end >= ssp1->end)
		  ssp1->type &= ~gS ;
		else
		  ssp->type &= ~gS ;
	      }
	}
      for (int ii = 0 ; ii < iMax  ; ii++)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->type & gS)
	    {
	      DSX *ssp1 ;
	      int jj, b1 = ssp->a2 + 1 ;
	      for (jj = ii + 1, ssp1 = ssp + 1 ; jj < iMax ; jj++, ssp1++)
		{
		  if (ssp1->type & gI) continue ;
		  if ((ssp1->type & gX ))
		    {
		      if (ssp1->a1 < b1 + 30 && ssp1->a1 > ssp->a2 + 1)
			ssp->a2 = ssp1->a1 - 1 ;
		      break ;
		    }
		}
	    }
	}

      for (int ii = 0 ; ii < iMax - 1 ; ii++)
	{ /* clean up stretches of r5p */
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->type & gReal5p)
	    {
	      DSX *ssp1 ;
	      int jj, b2, k ;
	      
	      /* discard gReal5p leaking above gS */
	      for (int pass = 0 ; pass < 2 ; pass++)
		for (jj = ii + 1, b2 = ssp->a2 + 1, ssp1 = ssp + 1, k = 1 ; jj < iMax && ssp1->a1 == b2 ; jj++, ssp1++, k++)
		  {
		    if (ssp1->type & gI)
		      {
			if (pass == 0)
			  b2 = ssp1->a2 + 1 ;
			continue ;
		      }
		    if (! (ssp1->type & gReal5p)) break ;
		    if (ssp1->type & gS)
		      {
			for (int m = 0 ; m < k ; m++)
			  if (! (ssp1[-m-1].type & gS))
			    ssp1[-m-1].type &= ~ gReal5p ;
			break ;
		      }
		    b2 = ssp1->a2 + 1 ;
		    if (0 && b2 > ssp->a1 + 30) break ;
		}
	    }
	}

      for (int pass = 0 ; pass < 2 ; pass++)
	for (int ii = 0 ; ii < iMax - 1 ; ii++)
	  { /* clean up stretches of r5p */
	    ssp = arrp (ss, ii, DSX) ;
	    if (ssp->type & gReal5p)
	      {
		DSX *ssp1 ;
		int jj, b2, k ;
		
		for (jj = ii + 1, b2 = ssp->a2 + 1, ssp1 = ssp + 1, k = 1 ; jj < iMax && ssp1->a1 == b2 ; jj++, ssp1++, k++)
		  {
		    if (ssp1->type & gI)
		      {
			if (pass == 0)
			  b2 = ssp1->a2 + 1 ;
			continue ;
		      }
		    if (! (ssp1->type & gReal5p)) break ;
		    if (ssp1->cover > ssp->cover)
		      for (int m = 0 ; m < k ; m++)
			ssp[m].cover = ssp1->cover ;
		    if (ssp1->end > ssp->end)
		      { ssp->end = ssp1->end ; ssp->endKey = ssp1->endKey ; }
		    ssp1->type &= ~gReal5p ;
		    b2 = ssp1->a2 + 1 ;
		  }
	      }
	  }


      for (jj = 1 ; jj < iMax ; jj++)
	{ /* drop exons above real5p */
	  ssp = arrp (ss, jj, DSX) ;
	  if (ssp->cover && (ssp->type & (gReal5p | gS)))
	    {  /* look for continuous exons above the real5p with no or real5p above */
	      int ii ;
	      BOOL ok = TRUE ;
	      for (ii = jj - 1 ; ok && ii >= 0 ; ii--)
		{ /* if they are weak */
		  DSX* ssp1 = arrp (ss, ii, DSX) ;
		  if (ssp->type & (gReal5p | gS))
		    ok = FALSE ;
		  if (ii == 0 || ssp1[-1].a2 + 1 != ssp1->a1)
		    break ;
		  if (ssp1->cover > ssp->cover / 2)
		    ok = FALSE ;
		}
	      if (ok)
		for ( ; ii < jj ; ii++)
		  {
		    DSX* ssp1 = arrp (ss, ii, DSX) ;
		    ssp1->cover = 0 ;
		  }
	    }
	}

      iMax = ssHappyFew (ss) ;
      for (jj = 1 ; jj < iMax ; jj++)
	{ /* drop internal exons above real5p */
	  ssp = arrp (ss, jj, DSX) ;
	  if (ssp->type & gReal5p)
	    {
	      DSX *ssp1 ;
	      int k = 0, jj, bestCover = ssp->cover, b2 = ssp->a1 - 1 ;
	      
	      for (jj = ii - 1, ssp1 = ssp - 1 ; jj >= 0 && ssp1->a2 == b2 ; jj--, ssp1--)
		{
		  if (ssp1->type & gI) { if (k > 0) { k = 0; break ; } continue ;}
		  if (ssp1->type & (gA | gReal5p | gS | gReal3p)) {k = 0 ; break ; }
		  if (ssp1->cover > bestCover / 2) {k = 0 ; break ; }
		  b2 = ssp1->a1 - 1 ; k++ ;
		}
	      for (jj = ii - 1, ssp1 = ssp - 1 ; k > 0 && jj >= 0 && ssp1->a2 == b2 ; jj--, ssp1--)
		{
		  if (ssp1->type & gI) continue ;
		  b2 = ssp1->a1 - 1 ;
		  ssp1->cover = 0 ; k-- ;
		}
	    }
	}
      iMax = ssHappyFew (ss) ;
      for (int ii = iMax -1 ; ii > 0 ; ii--)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (0 &&  (ssp->type & gReal5p))
	    {
	      DSX *ssp1, *ssp0 = ssp ;
	      int jj, bestCover = ssp->cover, b2 = ssp->a1 - 1 ;

	      for (jj = ii - 1, ssp1 = ssp - 1 ; jj >= 0 && ssp1->a2 == b2 ; jj--, ssp1--)
		{
		  if (ssp1->type & gI) continue ;
		  if (ssp1->type & (gA | gReal3p)) break ;
		  if (ssp1->cover < bestCover / 100)
		    {
		      for ( ; jj >= 0 && ssp1->a2 == b2 && (ssp1->type & gReal5p) ; jj--, ssp1--)
			{
			  if (ssp1->type & gI) continue ;
			  b2 = ssp1->a1 - 1 ;
			  ssp1->type &= ~ gReal5p ;
			}
		      break ;
		    }
		  else
		    {
		      if (ssp1->cover > bestCover) bestCover = ssp1->cover ;
		      ssp0->type &= ~gReal5p ;
		      ssp1->type |= gReal5p ;
		      ssp0 = ssp1 ;
		      b2 = ssp1->a1 - 1 ;
		    }
		}
	    }
	}
    }
  
  /* regularize the real3p */
  if (1)
    {
      /* rehook the 1 bp polyA to previous exons */
      for (int ii = 0 ; ii < iMax  ; ii++)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->type & (gA | gReal3p))
	    {
	      DSX *ssp1 ;
	      int jj, b2 = ssp->a1 - 1 ;
	      for (jj = ii - 1, ssp1 = ssp - 1 ; jj >= 0 ; jj--, ssp1--)
		{
		  if (ssp1->type & gI) continue ;
		  if ((ssp1->type & gX ))
		    {
		      if (ssp1->a2 > b2 - 30 && ssp->a1 > ssp1->a2 + 1)
			ssp->a1 = ssp1->a2 + 1 ;
		      break ;
		    }
		}
	    }
	}


      for (int ii = 0 ; ii < iMax - 1 ; ii++)
	{ /* clean up stretches of r3p */
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->type & gReal3p)
	    {
	      DSX *ssp1 ;
	      int jj, b2, k ;
	      
	      for (int pass = 0 ; pass < 2 ; pass++)
		{
		  for (jj = ii + 1, b2 = ssp->a2 + 1, ssp1 = ssp + 1, k = 1 ; jj < iMax && ssp1->a1 == b2 ; jj++, ssp1++, k++)
		    {
		      if (ssp1->type & gI)
			{
			  if (pass == 0)
			    b2 = ssp1->a2 + 1 ;
			  continue ;
			}
		      if (! (ssp1->type & gReal3p)) break ;
		    }
		  
		  for (int m = 0 ; m < k -1 ; m++)
		    ssp[m].type &= ~gReal3p ;
		  break ;
		}
	    }
	}

      for (int ii = 0 ; ii < iMax - 1 ; ii++)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->type & gReal3p)
	    {
	      DSX *ssp1, *ssp2 = ssp ;
	      int jj, b2 = ssp->a2 + 1, jjA = ii ;

	      /* search for a quasi contiguous polyA */
	      for (jj = ii + 1, ssp1 = ssp + 1 ; jj < iMax && ssp1->a1 <= b2+10 ; jj++, ssp1++)
		{
		  if (ssp1->type & gI) continue ;
		  if (ssp1->type & (gS | gReal5p))
		    break ;
		  if (0 && (ssp1->type & gReal3p) && (ssp1->end > 3 * ssp->end)) /* tyty */
		    jjA = jj ;
		  if (ssp1->type & gA)
		    { jjA = jj ; ssp2 = ssp1 ; }
		  if (ssp1->type & gReal3p)
		    { jjA = jj ; ssp2 = ssp1 ; }
		  if ((ssp1->type & gCompleteCDS) && ssp1->a2 < b2 + 30)
		    { jjA = jj ; ssp2 = ssp1 ; }
		  if (ssp1->a2 > ssp2->a2 + 50) break ;
		  b2 = ssp1->a2 + 1 ;
		}
	      /* raise their cover */
	      b2 = ssp->a2 + 1 ;
	      int bestEnd = ssp->end ;
	      KEY bestKey = ssp->endKey ;
	      for (jj = ii , ssp1 = ssp ; jj <= jjA ; jj++, ssp1++)
		{
		  if (ssp1->type & gI) continue ;
		  if (ssp1->cover < ssp->cover)
		    ssp1->cover = ssp->cover ;
		  if (ssp1->a1 > b2)
		    ssp1->a1 = b2 ;
		  if (ssp1->end < bestEnd)
		    { ssp1->end = bestEnd ; ssp1->endKey = bestKey ; }
		  if (ssp1->type & gReal3p)
		    if (ssp1->end > bestEnd)
		      { bestEnd = ssp1->end ; bestKey = ssp1->endKey ; }

		  b2 = ssp1->a2 + 1 ;
		  if (jj < jjA)
		    ssp1->type &= ~gReal3p ;
		  else
		    ssp1->type |= gReal3p ;
		}
	    }
	}
    }
    
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after regularize r5p/r3p\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after regularize r5p/r3p done\n") ;
    }

  if (iMax) /* merge SL in next exon */
     {
       for (jj = 0, ssp = arrp (ss, jj, DSX) ; jj < iMax - 1 ; ssp++, jj++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if ((ssp->type & (gReal5p | gS)) && ! (ssp->type & (gA | gReal3p)))
	     {
	       int i ;
	       DSX *ssp1 ;
	       for (i = jj + 1, ssp1 = arrp (ss, i, DSX) ; i < iMax ; ssp1++, i++)
		 {
		   if (ssp1->type & gI) continue ;
		   if (ssp1->cover <= 0) continue ;
		   if (ssp1->a1 == ssp->a2 + 1)
		     {
		       ssp1->type |= g5 ;
		     }
		   else
		     break ;
		 }
	     }
	 }
     }

  if (iMax) /* equalize the score of SL clusters */
     {
       for (jj = 0, ssp = arrp (ss, jj, DSX) ; jj < iMax ; ssp++, jj++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gS)
	     {
	       int i, k = 0, cMax = ssp->cover ;
	       DSX *ssp1 ;
	       int b1 = ssp->a2 + 1 ;
	       for (i = jj + 1, ssp1 = arrp (ss, i, DSX) ; i < iMax ; ssp1++, i++)
		 {
		   if (ssp1->type & gI)
		     continue ;
		   if (ssp1->a1 > b1)
		     break ;
		   if (ssp1->a1 > ssp->a1 + 40)
		     break ;
		   if (ssp1->type & gS)
		     {
		       k = i ;
		       if (ssp1->cover > cMax)
			 cMax = ssp1->cover ;
		     }
		   b1 = ssp1->a2 + 1 ;
		 }
	       for (i = jj ; i <= k ; i++)
		 {
		   ssp1 = arrp (ss, i, DSX) ;
		   if (ssp1->type & gI)
		     continue ;
		   ssp1->cover = cMax + i - jj ; /* by computing the lower one first we insure we do both */
		 }
	     }
	 }
     }

  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after merge SL\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after merge SL done\n") ;
    }

  if (0 && iMax) /* merge polyA in previous exon */
    {
      for (jj = iMax - 1, ssp = arrp (ss, jj, DSX) ; jj > 0 ; jj--, ssp--)
	{
	  int jj1 = jj ;
	  if (ssp->cover <= 0)   continue ;
	  if ((ssp->type & gX) && (ssp->type & (gA)) && ! (ssp->type & (g5 | gS)))
	    {
	      int i ;
	      int c = ssp->cover ;
	      DSX *ssp1 ;

	      for (i = jj - 1, ssp1 = arrp (ss, i, DSX) ; i >= 0 ; i--, ssp1--)
		{
		  if (ssp1->type & gI) continue ;
		  if (ssp1->cover <= 0) continue ;

		  if (ssp1->a2 >= ssp->a1 - 120)
		    {
		      jj1 = i ;
		      if (c < ssp1->cover)
			c = ssp1->cover ;
		    }
		  if (ssp1->type & gReal3p)
		    {
		      ssp->type |= gReal3p ;
		      ssp->cover = c ;
		      break ;
		    }
		}
	      for (i = jj - 1, ssp1 = arrp (ss, i, DSX) ; i >= jj1 ; i--, ssp1--)
		{
		  if (ssp1->type & gI) continue ;
		  if (ssp1->cover <= c) ssp1->cover = c ;
		  if (ssp->type & gReal3p)
		    ssp1->type &= ~gReal3p ;
		}
	      jj = jj1 ;
	    }
	}
    }
  
  /* keep happy few */
  iMax = ssHappyFew (ss) ;
  
  if (debug)
    {
      fprintf (stderr, "mrnaDesignGetGraphElements after merge same, SL, pA\n") ;
      showDs (sc, ss) ;
      fprintf (stderr, "mrnaDesignGetGraphElements after merge same, SL, pA done\n") ;
    }


  if (iMax) /* clip exons above first strong signal */
     {
       int bestCover = 0 ;
       for (jj = 0, ssp = arrp (ss, 0, DSX) ; jj < iMax ; jj++, ssp++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gI)
	     continue ;
	   if (ssp->cover > bestCover)
	     bestCover = ssp->cover ;
	   if (2 * ssp->cover < bestCover)
	     break ;
	   if (ssp->type & (gS |gReal5p))
	     break ;
	   if (ssp->type & gReal3p)
	     break ;
	 }
       if (ssp->type & (gS | gReal5p))
	 {
	   /* I found a real5p, kill all exons above */
	   for (int i = 0 ; i < jj ; i++)
	     if (ssp[-i-1].a2 < ssp->a1 && (10*ssp[-i - 1].cover <= bestCover || (ssp[-i-1].type & (gDF | gFF))))
	       ssp[-i - 1].cover = 0 ;
	     else
	       break ;
	   if (! (ssp->type & (gS | g5)))
	     {
	       ssp->a1 -= 0 ;  /* was 15 move back to original position */
	       if (ssp->a1 < 1)
		 ssp->a1 = 1 ;
	     }
	 }
     }
  if (iMax) /* clip exons after last strong signal */

    {
       int bestCover = 0 ;
       for (jj = iMax - 1, ssp = arrp (ss, jj, DSX) ; jj > 0 ; jj--, ssp--)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gI)
	     break ;
	   if (ssp->type & gReal5p)
	     break ;
	   if (ssp->type & gReal3p)
	     break ;
	   if (2 * ssp->cover < bestCover)
	     break ;
	   if (ssp->cover > bestCover)
	     bestCover = ssp->cover ;
	 }
       if (ssp->type & gReal3p)
	 {
	   /* I found a real3p, kill all exons below */
	   DSX *ssp1 ;
	   int i ;
	   for (i = iMax - 1, ssp1 = arrp (ss, i, DSX) ; i > jj ; i--, ssp1--)
	     {
	       if (ssp1->type & gA)
		 break ;
	       else
		 ssp1->cover = 0 ;
	     }
	 }
     }
  if (iMax) /* clip exons after last polyA */
    {
       int bestCover = 0 ;
       for (jj = iMax - 1, ssp = arrp (ss, jj, DSX) ; jj > 0 ; jj--, ssp--)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gI)
	     break ;
	   if (ssp->type & gA)
	     break ;
	   if (ssp->type & gReal3p)
	     break ;
	   if (2 * ssp->cover < bestCover)
	     break ;
	   if (ssp->cover > bestCover)
	     bestCover = ssp->cover ;
	 }
       if (2 * ssp->cover >= bestCover)
	 if (ssp->type & (gA | gReal3p))
	   {
	     /* I found a polyA, kill all exons below */
	     for (int i = jj + 1 ; i < iMax ; i++)
	       ssp[i - jj].cover = 0 ;
	   }
    }
  
  iMax = ssHappyFew (ss) ;
  if (iMax)  /* kill intron ending below last exon */
    {
      int b2 = 0 ;
      for (ii = iMax - 1 ; ii > 0 ; ii--)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->type & gX)
	    { b2 = ssp->a2 ; break ; }
	}
      for (ii = iMax - 1 ; ii > 0 ; ii--)
	{
	  ssp = arrp (ss, ii, DSX) ;
	  if (ssp->a2 > b2 && (ssp->type & gI))
	    ssp->cover = 0 ;
	}
    }      
  
  iMax = ssHappyFew (ss) ;
  if (0 &&    /* bad idea check worm_WS294 gene chr3__0_692 */
      iMax) /* create a cut over the SL */
    {
       for (ii = 0, ssp = arrp (ss, 0, DSX) ; ii < iMax ; ii++, ssp++)
         {
	   if (ssp->cover <= 0)   continue ;
	   if (ssp->type & gS)
	     {
	       int k, b2 = ssp->a1 - 1 ;
	       for (k = 1 ; k < ii ; k++)
		 {
		   DSX *ssp1 = ssp - k ;
		   if (ssp1->type & gI) continue ;
		   if (ssp1->a2 != b2) break ;
		   if (! (ssp1->type & gS) && ssp1->cover < ssp->cover/30)
		     { ssp1->cover = 0 ; break ; }
		   b2 = ssp1->a1 - 1 ;
		 }
	     }
	 }
    }

 done:
  /* keep happy few */
  iMax = ssHappyFew (ss) ;
  mrnaDesignCleanExons (sc, ss) ;
  ds->exons = ss ;
  
  if (0 || debug)
     {
       fprintf (stderr, "mrnaDesignGetGraphElements final exons/introns\n") ;
       showDs (sc, ds->exons) ;
       fprintf (stderr, "mrnaDesignGetGraphElements final exons/introns done\n") ;
     }

   ac_free (h) ;
   return ds->exons && arrayMax (ds->exons) ? TRUE : FALSE ;
} /* mrnaDesignGetGraphElements */

/**********************************************************************************/
/* extend the path Down via the best supported route */
static  BOOL N200 = FALSE ;  

static int mrnaDesignExtendDownRaw (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn0, int nIntron, int nStart, int nStop, int bestScore, int maxCover, long unsigned int flag)
{ 
  int length = 0 ;
  int a2, ii, nn, score = 0 ;
  int sMax = arrayMax (segs) ;
  DSX *up, *vp, *wp = 0 ;
  up = arrp (segs, nn0, DSX) ;
  a2 = up->a2 + 1 ;

  for (vp = up, ii = nn0 ; ii >= 0 ; ii--, vp--)
    {
      if (vp->donor == path)
	{
	  if (vp->type & gI)
	    break ;
	  if (bestScore < vp->score)
	    bestScore = vp->score ;
	}
    }

  for (nn = 0, vp = up + 1, ii = nn0 + 1 ; ii < sMax && vp->a1 <= a2 ; ii++, vp++)
    {
      if (vp->a1 != a2) continue ;
      if ((vp->flag & flag) != flag) continue ;
      if (N200 && ! (vp->type &  (gA | gS | gReal3p | gReal5p)) && 200 * vp->score < bestScore && 200 * vp->end < bestScore) continue ;
      if (score < vp->score) { score = vp->score ; nn = vp->nn ; wp = vp ; }
    }

  if (! score)
    return length ;

  if (0 && nStop && 
      (wp->type & (gS | gS0 | gReal5p))
      )
    return length ;
  
  if (0 && ! (wp->type & (gA | gReal3p)))
    {
      if (nStop)
	{
	  if (up->score / 10 > score)
	    return length ;
	}
      if (up->score / 100 > score)
	return length ;

      if (score < 0 * bestScore / 100)
	return length ;
    }

  length = vp->a2 - vp->a1 + 1 ;
  keySet (ks, arrayMax(ks)) = nn ;
  
  wp->donor = path ;
  
  if (wp->type & gI)
    { nIntron++ ; nStop = 0 ; }
  if (wp->type & gS)
    nStart++ ;
  if (wp->type & (gA | gReal3p))
    nStop++ ;
  score = ((wp->type & gI) ? wp->score/3 :  wp->score) ;
  if (score > bestScore)
    bestScore = score ;
  if (wp->cover > maxCover)
    maxCover = wp->cover ;
  if (0) fprintf (stderr, "Path=%d  %d:%d >> %d:%d\n", path, up->a1, up->a2, wp->a1, wp->a2) ;
  length += mrnaDesignExtendDownRaw (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, bestScore, maxCover, flag) ;

  return length ;
} /* mrnaDesignExtendDownRaw */

/***********/

static int mrnaDesignExtendDownLoopCDS (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn0, int nIntron, int nStart, int nStop, int bestScore, int maxCover, long unsigned int flag)
{ 
  int length = 0 ;
  int a2, ii, nn, score = 0 ; 
  int sMax = arrayMax (segs) ;
  DSX *up, *vp, *wp = 0 ;
  
  up = arrp (segs, nn0, DSX) ;
  a2 = up->a2 + 1 ;
  if (1) /* first try to chain the CDS */
    for (nn = 0, vp = up + 1, ii = nn0 + 1 ; ii < sMax && vp->a1 <= a2 ; ii++, vp++)
      {
	if (vp->a1 != a2) continue ;
	if ((vp->flag & flag) != flag) continue ;
	if (N200 && ! (vp->type &  (gA | gS | gReal3p | gReal5p)) && 200 * vp->score < bestScore && 200 * vp->end < bestScore) continue ;
	if (! (vp->type &  (gCompleteCDS)))
	  continue ;
	if (score < vp->cover) { score = vp->cover ; nn = vp->nn ; wp = vp ; }
      }
  if (! score)
    return   mrnaDesignExtendDownRaw (sc, ds, segs, ks, path, nn0, nIntron, nStart, nStop, bestScore, maxCover, flag) ;

  length = vp->a2 - vp->a1 + 1 ;
  keySet (ks, arrayMax(ks)) = nn ;
  wp->donor = path ;
  if (score > bestScore)
    bestScore = score ;
  if (wp->cover > maxCover)
    maxCover = wp->cover ;
  length += mrnaDesignExtendDownLoopCDS (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, bestScore, maxCover, flag) ;

  return length ;
} /* mrnaDesignExtendDownLoopCDS */

/***********/

static int mrnaDesignExtendDown (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn, int nIntron, int nStart, int nStop, int useCDS, int bestScore, int maxCover, long unsigned int flag)
{
  if (useCDS)
    return mrnaDesignExtendDownLoopCDS (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, bestScore, maxCover, flag) ;

  return mrnaDesignExtendDownRaw (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, bestScore, maxCover, flag) ;
} /* mrnaDesignExtendDown */

/**********************************************************************************/
/* extend the path Up via the best supported route */
static int mrnaDesignExtendUp (SC *sc, DS *ds, Array segs, KEYSET ks, int path, int nn0
			       , int nIntron, int nStart, int nStop, int useCDS, int bestScore, int maxCover, long unsigned int flag)
{ 
  int length = 0 ;
  int a1, ii, nn = 0, score = 0 ;
  int sMax = arrayMax (segs) ;
  DSX *up, *vp, *wp = 0 ;
  
  up = arrp (segs, nn0, DSX) ;
  if (up->nn != nn0) messcrash ("nn != nn0 in  mrnaDesignExtendUp") ;
  a1 = up->a1 - 1 ;

  if (useCDS)
    for (nn = ii = 0, vp = arrp (segs, 0, DSX) ; ii < sMax && vp->a1 <= a1 ; ii++, vp++)
      {
	if (vp->a2 != a1) continue ;
	if ((vp->flag & flag) != flag) continue ;
	if (! (vp->type &  (gCompleteCDS)))
	  continue ;
	if (N200 && ! (vp->type &  (gA | gS | gReal3p | gReal5p)) && 200 * vp->score < bestScore && 200 * vp->end < bestScore) continue ;
	if (score < vp->score) { score = vp->score ; nn = vp->nn ; wp = vp ; }
    }
  if (! score)
    {
      useCDS = 0 ;
      for (nn = ii = 0, vp = arrp (segs, 0, DSX) ; ii < sMax && vp->a1 <= a1 ; ii++, vp++)
	{
	  if (vp->a2 != a1) continue ;
	  if ((vp->flag & flag) != flag) continue ;
	  if (N200 && ! (vp->type &  (gA | gS | gReal3p | gReal5p)) && 200 * vp->score < bestScore && 200 * vp->end < bestScore) continue ;
	  if (score < vp->score) { score = vp->score ; nn = vp->nn ; wp = vp ; }
	}
    }

  if (! score)
    return length ;

  if (0 && (up->type & (gS | gReal5p)) && up->score > wp->score/2)
    return length ;

  if (0 && nStart &&
      (wp->type & (gA | gReal3p))
      )
    return length ;
    
  if (! (up->type & gI) && ! (wp->type & (gS | gS0 | gReal5p)))
    {
      if (up->type & gReal5p)
	{
	  if (up->score/ 10 > score)
	    return length ;
	}
      if (up->score / 100 > score)
	return length ;

      if (score < bestScore / 100)
	return length ;
    }

  /* find local exon score */
  if (0 && nStart && ! (wp->type & (gS | gReal5p | gCompleteCDS)))
    if (score < 0 * bestScore / 100)
      return length ;
    
  if (! useCDS)
    {
      if (0 && nStart && (wp->type & ( gA | gReal3p) && ! (wp->type & gCompleteCDS)))
	return 0 ;
      if (0 && nIntron && (wp->type & (gA | gReal3p)))
	return 0 ;
    }
  
  length = wp->a2 - wp->a1 + 1 ;
  keySet (ks, arrayMax(ks)) = nn ;
  
  wp->donor = path ;
  if (! useCDS)
    {
      if (wp->type & gI)
	{ nIntron++ ; nStart = 0 ; }
      if (wp->type & (gS | gReal5p))
	nStart++ ;
      if (wp->type & (gA | gReal3p))
	nStop++ ;
    }
  if (wp->cover > bestScore)
    bestScore = wp->cover ;
  if (wp->cover > maxCover)
    maxCover = wp->cover ;
  if (0) fprintf (stderr, "Path=%d  %d:%d << %d:%d\n", path, up->a1, up->a2, wp->a1, wp->a2) ;
  length += mrnaDesignExtendUp (sc, ds, segs, ks, path, nn, nIntron, nStart, nStop, useCDS, bestScore, maxCover, flag) ;

  return length ;
} /* mrnaDesignExtendUp */

/*********************************************************************/
/* get name of smallest Clone Group or smallest clone */
static void compositDesignSelectCosmid (SMRNA *tgp)
{
  KEY cosmid = tgp->cosmid ;
  int a1 = tgp->a1, a2 = tgp->a2 ;
  
  KEYSET ks = 0 ;
  KEY part, chrom = 0 ;
  OBJ Cosmid = 0 ;
  int ii, b1, b2, c1 = 0, c2 = 0, g1 = 0, g2 = 0 ;
  /* construct g1, g2 the IntMap coords of the tgene */
  if ((Cosmid = bsCreate (cosmid)))
    {
      if (bsGetKey (Cosmid, _IntMap, &chrom) &&
          bsGetData (Cosmid, _bsRight, _Int, &c1) &&
          bsGetData (Cosmid, _bsRight, _Int, &c2)) 
	{} ;
      bsDestroy (Cosmid) ;
    }
  if (c1 && c2)
    {
      if (c1 < c2) { g1 = c1 + a1 - 1 ; g2 = c1 + a2 - 1 ; }
      else { g1 = c1 - a1 + 1 ; g2 = c1 - a2 + 1 ; }
    }

  
  /* i try to see if my genes fits inside one of the parts */
  if (g1 && g2)
    {
      ks = queryKey (cosmid, ">Parts") ;
      for (ii = 0 ; ii < keySetMax (ks) ; ii++)
        {
          part = keySet (ks, ii) ;
          if ((Cosmid = bsCreate (part)))
            {
              if (bsGetKey (Cosmid, _IntMap, &chrom) &&
                  bsGetData (Cosmid, _bsRight, _Int, &b1) &&
                  bsGetData (Cosmid, _bsRight, _Int, &b2))
		{} ;
              bsDestroy (Cosmid) ;
            }
          if (b1 < b2 && b1 <= g1 && b1 <= g2 && b2 >= g1 && b2 >= g2)
            { a1 = g1 - b1 + 1 ; a2 = g2 - b1 + 1 ; cosmid = part ; break ; }
          if (b1 > b2 && b2 <= g1 && b2 <= g2 && b1 >= g1 && b1 >= g2)
            { a1 = b1 - g1 + 1 ; a2 = b1 - g2 + 1 ; cosmid = part ; break ; }
        }
      keySetDestroy (ks) ;
    }

  tgp->cosmid = cosmid ;
  tgp->a1 = a1 ;
  tgp->a2= a2 ;
  tgp->chrom = chrom ;
  tgp->g1 = g1 ;
  tgp->g2 = g2 ;
} /* compositDesignSelectCosmid (MRNA *tgp) */

/**********************************************************************************/

static int mrnaDesignExport (S2M *s2m, SC *sc, DS *ds, Array segs, KEYSET ks, int path, Array smrnas)
{
  int a1 = 0, ii,jj,  nn, maxScore = 0 ;
  DSX *up, *up2 ;
  SMRNA *smrna ;
  HIT *vp = 0 ;
  BOOL debug = FALSE ;
  AC_HANDLE h = ac_new_handle () ;
  
  Array exions = arrayHandleCreate (keySetMax(ks), DSX, h) ;
  smrna = arrayp(smrnas, path - 1, SMRNA) ;
  arrayDestroy (smrna->dnas) ;
  arrayDestroy (smrna->introns) ;
  arrayMax (smrnas) = path ;
  if (! smrna->hits)
    smrna->hits = arrayHandleCreate (keySetMax(ks), HIT, s2m->h) ;
  else
    smrna->hits = arrayReCreate (smrna->hits, keySetMax(ks), HIT) ;
  smrna->introns = keySetHandleCreate (ds->h) ;
  smrna->clones = keySetHandleCreate (s2m->h) ;
  keySetSort (ks) ; /* now the elements are in a1 order and correctly chained */
  for (ii = 0, maxScore = 100 ; ii < keySetMax (ks) ; ii++)
    {
      nn = keySet (ks, ii) ; 
      up = arrp (segs, nn, DSX) ;
      if (up->score > maxScore)
	maxScore = up->score ;
      up2 = arrayp (exions, ii, DSX) ;
      *up2 = *up ;
    }
   if (0 && debug)
     {
       fprintf (stderr, "mrnaDesignExport path %d exons/introns\n", path) ;
       showDs (sc, exions) ;
       fprintf (stderr, "mrnaDesignExport exons/introns done\n") ;
     }
  
   if (0 || debug)
     {
       fprintf (stderr, "### mrnaDesignExport path %d exons/introns\n", path) ;
       showDs (sc, exions) ;
       fprintf (stderr, "### mrnaDesignExport path %d exons/introns done\n", path) ;
     }

   for (vp = 0, ii = jj = 0 ; ii < keySetMax (ks) ; ii++)
    {
      up = arrp (exions, ii, DSX) ;
      if (!ii && (up->type & gI))
	continue ;
      if (! up->score)
	continue ;
      if (!vp)
	a1 = smrna->a1 = up->a1 ; 
      smrna->a2 = up->a2 ; 
      if (!vp || ((up->type & (gX | gI)) != (vp->type & (gX | gI))))
	{	
	  vp = arrayp (smrna->hits, jj++, HIT) ;
	  vp->a1 = up->a1 - a1 + 1 ;
	}
      vp->a2 = up->a2 - a1 + 1 ;
      if (up->type & gReal5p)
	vp->clipTop = up->end ;
      if (up->type & gReal3p)
	vp->clipEnd = up->end ;
      if (up->type & gReal3p)
	{
	  HIT *pA ;
	  if (! smrna->pA)
	    smrna->pA = arrayHandleCreate (8, HIT, ds->h) ;
	  pA = arrayp (smrna->pA, arrayMax (smrna->pA), HIT) ;
	  pA->a2 = up->a2 ; 
	  pA->type = gReal3p ; pA->est = up->endKey ; pA->cDNA_clone = up->cDNA_clone ;
	  if (up->cDNA_clone)
	    {
	      KEY est = keyGetKey (up->cDNA_clone, _Read) ;
	      pA->zone = keyGetInt (est, _Composite) ;
	    }
	}
      if (up->type & gA)
	{
	  HIT *pA ;
	  if (! smrna->pA)
	    smrna->pA = arrayHandleCreate (8, HIT, ds->h) ;
	  pA = arrayp (smrna->pA, arrayMax (smrna->pA), HIT) ;
	  pA->a2 = up->a2 ; 
	  pA->type = gA ; pA->est = up->endKey ; pA->cDNA_clone = up->cDNA_clone ;
	  if (up->cDNA_clone)
	    {
	      KEY est = keyGetKey (up->cDNA_clone, _Read) ;
	      pA->zone = keyGetInt (est, _Composite) ;
	    }
	}
      if (up->type & gS)
	{
	  HIT *pA ;
	  if (! smrna->pA)
	    smrna->pA = arrayHandleCreate (8, HIT, ds->h) ;
	  pA = arrayp (smrna->pA, arrayMax (smrna->pA), HIT) ;
	  pA->a1 = up->a1 - smrna->a1 + 1 ;
	  pA->type = gS ;
	  pA->est = keyGetKey (up->slKey, _Transpliced_to) ;
	  if (pA->est)
	    pA->cDNA_clone = keyGetKey (up->slKey, _cDNA_clone) ;
	  else
	    lexaddkey ("SLz", &pA->est, _VMotif) ;
	  if (up->cDNA_clone)
	    {
	      KEY est = keyGetKey (up->cDNA_clone, _Read) ;
	      pA->zone = keyGetInt (est, _Composite) ;
	    }
	}
      if (vp->zone < up->cover) vp->zone = up->cover ;
      vp->type |= ( up->type  & (~ gFF)) ;
      if (jj > 1)
	vp->type &= (~ (gS | gS0 | gD | gReal5p | g5 | g3 | gDF)) ;
    }
   if (smrna->pA)
     {
       for (ii = 0 ; ii < arrayMax (smrna->pA) ; ii++)
	 {
	   HIT *pA = arrayp (smrna->pA, ii, HIT) ;
	   if (pA->type == gReal3p)
	     pA->a2 = pA->a2 - smrna->a2 ;
	   if (pA->type == gA)
	     {
	       pA->a2 = pA->a2 - smrna->a2 ;
	       if (pA->a2 < - 100 || pA->a2 > 20)
		 pA->type = 0 ;
	     }
	 }
     }

  for (ii = 0 ; ii < jj - 1 ; ii++)
    {
      vp = arrayp (smrna->hits, ii, HIT) ;
      if (vp->type & gI)
	{
	  SMRNA mm ;
	  KEY intron ;
	  if (sc->a1 < sc->a2)
	    {
	      mm.a1 = vp->a1 + sc->a1 + a1 - 2 ;
	      mm.a2 = vp->a2 + sc->a1 + a1 - 2 ;
	    }
	  else
	    {
	      mm.a1 = - vp->a1 + sc->a1 - a1 + 2 ;
	      mm.a2 = - vp->a2 + sc->a1 - a1 + 2 ;
	    }
	  mm.cosmid = s2m->cosmid ;
	  compositDesignSelectCosmid (&mm) ;
	  char *iName = messprintf ("%s__%d_%d", name (mm.chrom), mm.g1, mm.g2) ;
	  lexaddkey (iName, &intron, _VIntron) ;
	  keySet (smrna->introns, keySetMax (smrna->introns)) = intron ;
	}
    }
  if (debug)
    showHits (smrna->hits) ;

  ac_free (h) ;
  return path ;
} /* mrnaDesignExport G_t_NC_003281.10_1_4561 */

/**********************************************************************************/

static BOOL mrnaDesignIsNewPath (Array ss, Array ksPaths, KEYSET ks0, int path, int tested, int xStart, int xStop, AC_HANDLE h)
{
  BOOL ok = TRUE ;
  int PM = 0 ;
  KEYSET ks = keySetHandleCopy (ks0, h) ;
  int iMax = keySetMax (ks) ;
  keySetSort (ks) ;

  if (path > PM && iMax && xStart)
    {
      int ii = keySet (ks, 0) ;
      DSX *ssp = arrp (ss, ii, DSX) ;
      if (xStart != ssp->a1)
	ok = FALSE ;
    }
  if (path > PM && iMax && xStop)
    {
      int ii = keySet (ks, iMax - 1) ;
      DSX *ssp = arrp (ss, ii, DSX) ;
      if (xStop != ssp->a2)
	ok = FALSE ;
    }

  if (1)
    {
      int minCover = 0, maxCover = 0, maxC = 0 ;
      int nI = 0 ;
      int nEnds = 0 ;
      
      for (int i = 0 ; i < iMax ; i++)
	{
	  int jj = keySet (ks, i) ;
	  DSX *ssp = arrp (ss, jj, DSX) ;
	  if (ssp->type & (gS | gReal3p | gReal5p))
	    nEnds++ ;
	  if (ssp->cover > maxC) maxC = ssp->cover ;
	  if (ssp->type & gI)
	    {
	      nI++ ;
	      if (! maxCover) minCover = ssp->cover ;
	      if (ssp->cover > maxCover) maxCover = ssp->cover ;
	      if (ssp->cover < minCover) minCover = ssp->cover ;
	    }
	}
      if (path > PM && 200 * minCover < 1 * maxCover)
	ok = FALSE ;
      if (nI + nEnds == 0 && maxC < 1000)
	ok = FALSE ;
    }
  
  for (int ii = 0 ; ok && ii < arrayMax (ksPaths) ; ii++)
    {
      KEYSET ks1 = array (ksPaths, ii, KEYSET) ;
      ok = FALSE ;
      if (! ks1 || keySetMax (ks1) != iMax)
	ok = TRUE ;
      else
	{
	  for (int i = 0 ; !ok && i < iMax ; i++)
	    if (keySet (ks, i) != keySet (ks1, i))
	      ok = TRUE ;
	}
    }

  /* flag the paths */
  for (int i = 0 ; i < iMax ; i++)
    {
      int ii = keySet (ks, i) ;
      DSX *up = arrp (ss, ii, DSX) ;
      if (ok) up->path = path ;
      up->tested = tested ;
      up->donor = 0 ;
    }
  if (ok)
    array (ksPaths, arrayMax(ksPaths), KEYSET) = ks ;
  else
    keySetDestroy (ks) ;
  return ok ;
} /* mrnaDesignIsNewPath */

/**********************************************************************************/

static int mrnaDesignFlagDescendants (Array ss, int ii, long unsigned int flag, int maxRatio, int maxI)
{
  int jj, n = 0 ;
  int iMax = arrayMax (ss) ;
  DSX *vp, *up = arrp (ss, ii, DSX) ;
  int a1 = up->a2 + 1 ;
  
  up->flag |= flag ;
  if ((up->type & gI) && up->cover > maxI)
    maxI = up->cover ;
  for (jj = ii + 1, vp = up + 1 ; jj < iMax ; vp++, jj++)
    {
      if (vp->a1 > a1)
	break ;
      if (vp->flag & flag)
	continue ;
      if ((vp->type & gI) && maxI > maxRatio * vp->cover)
	continue ;
      if (vp->a1 == a1)
	mrnaDesignFlagDescendants (ss, jj , flag, maxRatio, maxI) ;
      n++ ;
    }

  return n ;
} /* mrnaDesignFlagDescendants */

/**********************************************************************************/

static int mrnaDesignFlagAscendants (Array ss, int ii, long unsigned int flag, int maxRatio, int maxI)
{
  int n = 0 ;
  int jj ;
  DSX *vp, *up = arrp (ss, ii, DSX) ;
  int a2 = up->a1 - 1 ;
  
  up->flag |= flag ; 
  if ((up->type & gI) && up->cover > maxI)
    maxI = up->cover ;

  for (jj = ii - 1, vp = up - 1 ; jj >= 0 ; vp--, jj--)
    {
      if (vp->flag & flag)
	continue ;
      if ((vp->type & gI) && maxI > maxRatio * vp->cover)
	continue ;
      if (vp->a2 == a2)
	mrnaDesignFlagAscendants (ss, jj , flag, maxRatio, maxI) ;
      n++ ;
    }
  return n ;
} /* mrnaDesignFlagAscendants */

/**********************************************************************************/

static void mrnaDesignSetEndShades (KEYSET light, KEYSET dark)
{
  int nLight = keySetMax (light) ;
  int nDark = keySetMax (dark) ;

  if (nLight + nDark)
    {
      AC_HANDLE h = ac_new_handle () ;
      vTXT txt = vtxtHandleCreate (h) ;
      const char *errors = 0 ;
      AC_DB db = ac_open_db (0, &errors) ; /* local database, cannot fail */
      
      keySetSort (light) ; keySetCompress (light) ;
      keySetSort (dark) ; keySetCompress (dark) ;
      for (int i = 0 ; i < nLight ; i++)
	{
	  const char *nam = name(keySet(light, i)) ;
	  vtxtPrintf (txt, "Sequence %s\n", nam) ;
	  if (! strncmp (nam, "Xends_ELF", 9) || ! strncmp (nam, "Xends_ERR", 9))
	    vtxtPrintf (txt, "Colour PALECYAN\n") ;
	  if (! strncmp (nam, "Xends_ERF", 9) || ! strncmp (nam, "Xends_ELR", 9))
	    vtxtPrintf (txt, "Colour PALEYELLOW\n") ;
	  vtxtPrintf (txt, "\n") ;
	}
      for (int i = 0 ; i < nDark ; i++)
	{
	  const char *nam = name(keySet(dark, i)) ;
	  vtxtPrintf (txt, "Sequence %s\n", nam) ;
	  if (! strncmp (nam, "Xends_ELF", 9) || ! strncmp (nam, "Xends_ERR", 9))
	    vtxtPrintf (txt, "Colour CYAN\n") ;
	  if (! strncmp (nam, "Xends_ERF", 9) || ! strncmp (nam, "Xends_ELR", 9))
	    vtxtPrintf (txt, "Colour YELLOW\n") ;
	  vtxtPrintf (txt, "\n") ;
	}
      ac_parse (db, vtxtPtr (txt), &errors, 0, h) ;
      
      ac_free (h) ;
    }
  return ;
} /*  mrnaDesignSetEndShades */

/**********************************************************************************/

static int mrnaDesignFindStartEndPairs (Array ss, Array ss2, Array sFlags, Array starts, Array stops)
{
  BOOL debug = FALSE ;
  AC_HANDLE h = ac_new_handle () ;
  DSX *up, *up2, *vp = 0, *wp = 0, *sFlag = 0 ;
  int i, iMax = arrayMax (ss), iStart = 0, iStop = 0, startCover = 0, stopCover = 0 ;
  int nFlags = 0 ;
  long unsigned int un = 1 ;
  KEYSET light = keySetHandleCreate (h) ;
  KEYSET dark = keySetHandleCreate (h) ;
  int nLight = 0, nDark = 0 ;
  
  for (i = 0, up = arrp (ss, 0, DSX) ; i < iMax ; i++, up++)
    {
      if (! (up->type & gS) && (up->type & gReal5p) && up->cover > 20 * up->end)
	{ up->type &= ~gReal5p ; keySet (light, nLight++) = up->endKey ; }
      if ((up->type & gReal3p) && up->cover > 20 * up->end)
      	{ up->type &= ~gReal3p ; keySet (light, nLight++) = up->endKey ; }

    }      
  
  for (i = 0, up2 = arrp (ss2, 0, DSX) ; i < iMax ; i++, up2++)
    {
      up = arrp (ss, up2->nn, DSX) ;
      up->path = i ; /* temporary */
    }

  for (i = 0, up = arrp (ss, 0, DSX), vp = 0 ; i < iMax ; i++, up++)
    {
      if (up->type & (gS | gReal5p))
	{
	  if (!vp || vp->a2 < up->a1 + 1)
	    {
	      vp = arrayp (starts, iStart++, DSX) ;
	      *vp = *up ;
	      vp->nn = up->path ;
	      vp->path = i ;
	      vp->flag = 0 ;
	      vp->type = up->type & gS ;
	      startCover += vp->end ;
	    }
	  vp = up ;
	}
    }

  for (i = iMax - 1, up = arrp (ss, i, DSX), wp = 0 ; i >= 0 ; i--, up--)
    {
      if (up->type & gReal3p)
	  {
	    if (! wp || wp->a1 > up->a2 + 1)
	      {
		wp = arrayp (stops, iStop++, DSX) ;
		*wp = *up ;
		wp->nn = up->path ;
		wp->path = i ;
		wp->flag = 0 ;
		wp->type = 0 ;
		stopCover += wp->end ;
	    }
	    wp = up ;
	  }
    }

  for (i = 0, up2 = arrp (ss2, 0, DSX), vp = 0 ; i < iMax ; i++, up2++)
    { up2->path = 0 ; up2->flag = 0 ; }

  /* normalize */
  if (startCover > 0)
    for (int i = 0 ; i < iStart ; i++)
      {
	DSX *vp = arrayp (starts, i, DSX) ;
	vp->end *=  10000.01 / startCover ;
	if (vp->type & gS)
	  vp->end += 10000 ;
	else
	  vp->type = 0 ;
	vp->flag = 0 ;
	if ((vp->flag & gS) || vp->end > 1000)
	  keySet (dark, nDark++) = vp->endKey ;
	else
	  keySet (light, nLight++) = vp->endKey ;
      }
  if (stopCover > 0)
    for (int i = 0 ; i < iStop ; i++)
      {
	DSX *wp = arrayp (stops, i, DSX) ;
	wp->type = 0 ; wp->flag = 0 ;
	wp->end *=  10000.01 / stopCover ;
	if (wp->end > 1000)
	  keySet (dark, nDark++) = wp->endKey ;
	else
	  keySet (light, nLight++) = wp->endKey ;
      }
  mrnaDesignSetEndShades (light, dark) ;
  /* sort */
  iStart = ssHappyFew (starts) ;
  iStop = ssHappyFew (stops) ;
  if (iStart > 1) arraySort (starts, dsEndOrder) ;
  if (iStop > 1) arraySort (stops, dsEndOrder) ;
  if (iStart > 15) iStart = arrayMax (starts) = 15 ;
  if (iStop > 15) iStop = arrayMax (stops) = 15 ;
  if (0 || debug)
    {
      for (int i = 0 ; i < iStart ; i++)
	{
	  vp = arrayp (starts, i, DSX) ;
	  fprintf (stderr, "..start %lu : a1=%d startCover %d %% %ld\n"
		   , un << i, vp->a1, vp->end, (long) startCover * vp->end / 100 
		   ) ;
	}
      for (int i = 0 ; i < iStop ; i++)
	{
	  wp = arrayp (stops, i, DSX) ;
	  fprintf (stderr, "..stop %lu : a2=%d stopCover %d %% %ld\n"
		   , un << i, wp->a2, wp->end, (long) stopCover * wp->end / 100
		   ) ;
	}
    }
  
  /* flag the segs and all their descendants or ascendants */
  for (int i = 0 ; i < iStart ; i++)
    {
      long unsigned int flag = un << i ;
      vp = arrayp (starts, i, DSX) ;
      mrnaDesignFlagDescendants (ss, vp->path, flag, 200, 0) ;
    }
  for (int i = 0 ; i < iStop ; i++)
    {
      long unsigned int flag = un << (i + 32) ;
      vp = arrayp (stops, i, DSX) ;
      mrnaDesignFlagAscendants (ss, vp->path, flag, 200, 0) ;
    }

  /* label the good starts (strongest or saves a new intron or a new stop) */
  long unsigned int allFlag = 1 ;
  if (iStart)
    {
      long unsigned int flag2 = 1 ;
      vp = arrp (starts, 0, DSX) ;
      vp->type |= 1 ; allFlag |= flag2 ; vp->flag = flag2 ; 
      for (int ii = 1 ; ii < iStart ; ii++)
	{ /* regarder si on decouvre de nouveaux introns ou stops */
	  BOOL ok = FALSE ;

	  flag2 = un << ii ;
	  vp = arrp (starts, ii, DSX) ;
	  if ((vp->type & gS) || vp->end > 6000)  
	    ok = TRUE ;
	  else if (vp->end > 0)
	    ok = TRUE ;
	  else if (! (vp->flag & allFlag))
	    {
	      int i ;
	      for (i = 0, up = arrp (ss, 0, DSX) ; ! ok && i < iMax ; i++, up++)
		if ((up->flag & flag2) && ! (up->flag & allFlag))
		  {
		    if (up->type & gI)
		      { /* this starts yields a new intron */
			ok = TRUE ;
		      }
		    if (vp->end > 100 && (up->type & gX))
		      { /* this starts yields a new exon */
			ok = TRUE ;
		      }
		  }
	    }
	  if (ok)
	    { vp->type |= 1 ; allFlag |= flag2 ; vp->flag = flag2 ; }
	}
    }
  if (iStop)
    {
      long unsigned int flag2 = un << (32) ;
      wp = arrp (stops, 0, DSX) ;
      wp->type = 1 ; allFlag |= flag2 ; wp->flag = flag2 ;
      for (int ii = 1 ; ii < iStop ; ii++)
	{ /* regarder si on decouvre de nouveaux introns ou starts */
	  BOOL ok = FALSE ;

	  flag2 = un << (32 + ii) ;
	  wp = arrp (stops, ii, DSX) ;
	  if (wp->end > 6000)
	    ok = TRUE ;
	  else if (! (wp->flag & allFlag))
	    {
	      int i ;
	      for (i = 0, up = arrp (ss, 0, DSX) ; ! ok && i < iMax ; i++, up++)
		if ((up->flag & flag2) && ! (up->flag & allFlag))
		  {
		    if (up->type & gI)  /* this starts yields a new intron */
		      ok = TRUE ;
		    if (wp->end > 100 && up->type & (gX))  /* this starts yields a new exon */
		      ok = TRUE ;
		}
	    }
	  if (ok)
	    { wp->type = 1 ; allFlag |= flag2 ; wp->flag = flag2 ;}
	}
    }
  
  /* create the flag pairs */
  sFlag = arrayp (sFlags, 0, DSX) ; /* default, but do not set nFlags++ */
  long unsigned int flag3 = 0 ;
  if (iStart)
    {
      for (int ii = 0 ; ii < iStart ; ii++)
	{
	  vp = arrp (starts, ii, DSX) ;
	  if (vp->type)
	    { /* foreach good start, select the highest matching stop */
	      BOOL ok = FALSE ;
	      long unsigned int flag = un << (ii) ;
	      for (int j = 0 ; ! ok && j < iStop ; j++)
		{
		  wp = arrp (stops, j, DSX) ;
		  if (wp->a2 > vp->a1)
		    {
		      up = arrp (ss, wp->path, DSX) ;
		      if (up->flag & flag)
			{
			  long unsigned int flag2 = un << (32 + j) ;
			  flag3 |= flag2 ;
			  if (debug)
			    fprintf (stderr, "Flag %d, a1 = %d, a2 = %d, score = %d / %d, flag = %lu / %lu\n"
				     , nFlags, vp->a1, wp->a2, vp->end, wp->end, flag, flag2 >> 32) ; 
			  sFlag = arrayp (sFlags, nFlags++, DSX) ;
			  sFlag->flag = flag | flag2 ;
			  if (vp->type & gS) sFlag->nn = vp->nn + 1 ;
			  sFlag->nn = vp->nn + 1 ;
			  ok = TRUE ;
			}
		    }
		}
	      if (! ok)
		{
		  if (debug)
		    fprintf (stderr, "Flag %d, a1 = %d, a2 = open, score = %d, flag = %lu\n"
			     , nFlags, vp->a1, vp->end, flag) ; 
		  sFlag = arrayp (sFlags, nFlags++, DSX) ;
		  sFlag->flag = flag ;

		}
	    }
	  else if (vp->end > 10)
	    { /* we are in a descendant start, select its closest stop */
	      int b2 = 0 ;
	      int bestj = -1 ;
	      long unsigned int flag = un << (ii) ;
	      for (int j = 0 ; j < iStop ; j++)
		{
		  wp = arrp (stops, j, DSX) ;
		  if (wp->a2 > vp->a1)
		    {
		      if (wp->end > vp->end/2)
			{
			  up = arrp (ss, wp->path, DSX) ;
			  if ((up->flag & flag) && (!b2 || up->a2 < b2))
			    {
			      b2 = up->a2 ;
			      bestj = j ;
			    }
			}
		    }
		}
	      if (bestj >= 0)
		{
		  long unsigned int flag2 = un << (32 + bestj) ;
		  wp = arrp (stops, bestj, DSX) ;
		  flag3 |= flag2 ;
		  if (debug)
		    fprintf (stderr, "Descendant flag %d, a1 = %d, a2 = %d, score = %d / %d, flag = %lu\n"
			     , nFlags, vp->a1, wp->a2, vp->end, wp->end, flag) ; 
		  sFlag = arrayp (sFlags, nFlags++, DSX) ;
		  sFlag->flag = flag | flag2 ;
		}
	      else
		{
		  if (debug)
		    fprintf (stderr, "Descendant flag %d, a1 = %d, a2 = open, score = %d, flag = %lu\n"
			     , nFlags, vp->a1, vp->end, flag) ; 
		  sFlag = arrayp (sFlags, nFlags++, DSX) ;
		  sFlag->flag = flag ;
		}
	    }
	}
    }
  if (iStop)
    {
      for (int j = 0 ; j < iStop ; j++)
	{
	  wp = arrp (stops, j, DSX) ;
	  long unsigned int flag2 = un << (32 + j) ;
	  if (! (flag2 & flag3)) /* new stop */
	    if (wp->type || wp->end > 10)
	      {
		/* select the closest start */
		int b1 = 0 ;
		int besti = -1 ;
		for (int i = 0 ; i < iStart ; i++)
		{
		  vp = arrp (starts, i, DSX) ;
		  if (vp->a1 < wp->a2)
		    {
		      if (vp->end > wp->end/2)
			{
			  up = arrp (ss, vp->path, DSX) ;
			  if ((up->flag & flag2) && (!b1 || up->a1 > b1))
			    {
			      b1 = up->a1 ;
			      besti = i ;
			    }
			}
		    }
		}
	      if (besti >= 0)
		{
		  long unsigned int flag = un << (besti) ;
		  vp = arrp (starts, besti, DSX) ;
		  flag3 |= flag2 ;
		  if (debug)
		    fprintf (stderr, "New stop Flag %d, a1 = %d, a2 = %d, score = %d, flags %lu / %lu\n"
			     , nFlags, vp->a1, wp->a2, wp->end, flag, flag2 >> 32) ;
		  sFlag = arrayp (sFlags, nFlags++, DSX) ;
		  sFlag->flag = flag | flag2 ;
		  sFlag->nn = wp->nn + 1 ;
		}
	      else
		{
		  if (debug)
		    fprintf (stderr, "New stop Flag %d, a1 = open, a2 = %d, score = %d, flag = %lu >> 32\n"
			     , nFlags, wp->a2, wp->end, flag2) ; 
		  sFlag = arrayp (sFlags, nFlags++, DSX) ;
		  sFlag->flag = flag2 ;
		  sFlag->nn = wp->nn + 1 ;
		}
	      }
	}
    }
  if (0) arraySort (sFlags, sFlagsOrder) ;
  for (int pass = 0 ; nFlags && pass < 2 ; pass++) /* check if we captured all introns and all cds */
    {
      DSX *up, *up2 ;
      int ii ;
      int myFlag = (pass == 0 ? gI : gCompleteCDS) ;
      long unsigned int un = 1 ;
      long unsigned int mask = (un << 32) - 1 ;
      for (ii = 0, up2 = arrp (ss2, 0, DSX) ; ii < iMax ; ii++, up2++)
	{
	  up = arrp (ss, up2->nn, DSX) ;
	  if (up->type & myFlag) 
	    { /* is this intron is part of a flag pair */
	      BOOL ok = FALSE ;
	      for (int i = 0 ; ! ok && i < nFlags ; i++)
		{
		  sFlag = arrp (sFlags, i, DSX) ;
		  flag3 = sFlag->flag ;
		  if ((up->flag & flag3) == flag3)
		    ok = TRUE ;
		}
	      if (! ok)
		{ /* create a new flag pair */
		  long unsigned int flag1 = 0, flag2 = 0 ;
		  if (up->flag & mask) /* take the best start */
		    {
		      flag1 = 1 ;
		      for (int i = 0 ; i < 32 ; i++)
			if ((up->flag >> i) & 1)
			  break ;
			else
			  flag1 <<= 1 ;
		    }
		  else
		    flag1 = 0 ;
		  if (up->flag >> 32) /* take the best stop */
		    {
		      long unsigned int un = 1 ;
		      flag2 = un << 32 ;
		      for (int i = 0 ; i < 32 ; i++)
			if ((up->flag >> (32+i)) & 1)
			  break ;
			else
			  flag2 <<= 1 ;
		    }
		  else
		    flag2 = 0 ;
		  if (debug)
		    fprintf (stderr, "New pair Flag %d : %lu/%lu inherited from the intron %d/%d\n"
			     , nFlags, flag1, flag2 >> 32, up->a1, up->a2
			     ) ;
		  sFlag = arrayp (sFlags, nFlags++, DSX) ;
		  sFlag->nn = ii + 1 ;
		  sFlag->flag = flag1 | flag2 ;
		}
	    }
	}
      arrayCompress (sFlags) ;
      nFlags = arrayMax (sFlags) ;
    }

  return nFlags ; 
} /* mrnaDesignFindStartEndPairs */

/**********************************************************************************/

static int mrnaDesignFindPaths (S2M *s2m, SC *sc, DS *ds, Array smrnas)
{
  BOOL debug = TRUE ;
  int iFlag, path, path0, tested, nIntron, nStart, nStop, useCDS ;
  int eeMax = ssHappyFew (ds->exons) ;
  Array segs, segs2 ;
  DSX *vp, *vp2, *sFlag ; ;
  KEYSET ks = keySetHandleCreate (ds->h) ;
  Array ksPaths = arrayHandleCreate (20, KEYSET, ds->h) ;
  long unsigned int flag = 0 ;
  Array sFlags = arrayHandleCreate (32, DSX, ds->h) ;
  Array starts = arrayHandleCreate (32, DSX, ds->h) ;
  Array stops = arrayHandleCreate (32, DSX, ds->h) ;
  
  arraySort (ksPaths, dsScoreOrder) ; /* for computer happiness */

  
  if (1 && debug)
     {
       fprintf (stderr, "mrnaDesignFindPaths start\n") ;
       showDs (sc, ds->exons) ;
       fprintf (stderr, "mrnaDesignFindPaths start done\n") ;
     }
  
  /* sort a copy by cover */
  segs = ds->exons ;
  segs2 = arrayHandleCopy (segs, ds->h) ;
  arraySort (segs2, dsCoverOrder) ;

  if (arrayMax (segs))
    mrnaDesignFindStartEndPairs (segs, segs2, sFlags, starts, stops) ;

  int ii ;
  for (ii = 0, vp = arrp (segs, 0, DSX) ; ii < eeMax ; ii++, vp++)
    vp->path = 0 ;
  
  /* construct all paths
   * recursion on start-end pairs specified by flags
   */
  for (path = path0 = 0, iFlag = 0 ; iFlag < arrayMax (sFlags) ; path0 = path, iFlag++)
    for (useCDS = 0, sFlag = arrp (sFlags, iFlag, DSX) ; useCDS >=0 ; useCDS--)
      while (1)
	{
	  int ii2 = 0, bestIi2 = -1 ;
	  flag = sFlag->flag ;
	  vp2 = arrp (segs2, 0, DSX) ;
	  int bigScore = 0 ;
	  int maxCover = 0 ;
	  long unsigned int un = 1 ;
	  long unsigned int mask = (un << 32) - 1 ;
	  long unsigned int flag1 = flag & mask ;
	  long unsigned int flag2 = flag ^ flag1 ;
	  int iStart = -1, iStop = -1 ;
	  int xStart = 0, xStop = 0 ;
	  
	  if (flag1)
	    {
	      int i = 0 ;
	      long unsigned int flag3 = flag1 ;
	      for (i = 0 ; flag3 ; flag3 >>=1, i++)
		if (flag3 & 1)
		  break ;
	      DSX *wp = arrp (starts, i, DSX) ;
	      iStart = wp->path ;
	      xStart = wp->a1 ;
	    }
	  if (flag2)
	    {
	      int i = 0 ;
	      long unsigned int flag3 = flag2 >> 32 ;
	      for (i = 0 ; flag3 ; flag3 >>=1, i++)
		if (flag3 & 1)
		  break ;
	      DSX *wp = arrp (stops, i, DSX) ;
	      iStop = wp->path ;
	      xStop = wp->a2 ;
	    }
	  
	  for (ii2 = 0, vp2 = arrp (segs2, 0, DSX) ; ii2 < eeMax ; ii2++, vp2++)
	    {
	      /* find highest score for that path */
	      vp = arrp (segs, vp2->nn, DSX) ;
	      if ((vp->flag & flag) != flag)
		continue ;
	      if (vp->cover > bigScore)
		bigScore = vp->cover ;
	    }
	  if ((! flag1 || ! flag2)&&  sFlag->nn)
	    {
	      vp2 = arrp (segs2, sFlag->nn - 1, DSX) ;
	      vp = arrp (segs, vp2->nn, DSX) ;
	      bestIi2 = sFlag->nn - 1 ;
	      maxCover = vp->cover ;
	      sFlag->nn = 0 ;
	      /* sFlag->nn = 0 ;  tyty */
	    }
	  else
	    for (ii2 = 0, vp2 = arrp (segs2, 0, DSX) ; ii2 < eeMax ; ii2++, vp2++)
	      {
		/* find highest scoring seg not yet incorporated in a path */
		vp = arrp (segs, vp2->nn, DSX) ;
		if (vp->tested == iFlag + 1)
		  continue ;
		if ((vp->flag & flag) != flag)
		  continue ;
		if (! vp->score)
		  continue ;
		if (0 && vp->path && vp2->nn != iStart && vp2->nn != iStop)
		  continue ;
		if (vp->path > path0)  /* already used for this pair of sflags */
		  continue ;
		if ((vp->type & gX) && vp->score < bigScore / 100)
		  continue ;
		if ((vp->type & gI) && vp->score < bigScore / 500)
		  continue ;
		if (vp2->nn != iStart && (vp->type & gDF))
		  continue ;
		if (vp2->nn != iStop && (vp->type & gFF))
		  continue ;
		if (0 && useCDS && ! (vp->type &  gCompleteCDS))
		  continue ;
		/* do not start on a mini exon */
		if (vp->a2 < vp->a1 + 10 && (vp->type & gX) && (!(vp->type & (gS | gReal5p | gA | gReal3p))))
		  continue ;
		if (vp->cover > maxCover)
		  {
		    maxCover = vp->cover ;
		    bestIi2 = ii2 ;
		  }
	      }
	  
	  if (bestIi2 < 0)
	    break ;
	  ii2 = bestIi2 ;

	  vp2 = arrp (segs2, ii2, DSX) ;
	  vp = arrp (segs, vp2->nn, DSX) ;
	  vp->tested = tested = iFlag +1 ; /* avoid zero */

	  path++ ;
	  if (0 || debug)
	    {
	      fprintf (stderr, "+++New path %d, Flag %d = %lu/%lu, start on %d %d score=%d maxCover=%d\n", path, iFlag, flag & mask, flag >> 32, vp->a1, vp->a2, vp->score, maxCover) ;
	    }
	  ks = keySetReCreate (ks) ;
	  keySet (ks, 0) = vp->nn ;
	  vp->donor = path ;
	  if (1) /* ! (vp->type & gB) */
	    {
	      int bestScore = ((vp->type & gI) ? vp->score/3 :  vp->score) ;
	      nIntron = (vp->type & gI) ? 3 : 0 ;
	      nStart = 1 ;
	      nStop = (vp->type & (gA | gReal3p)) ? 1 : 0 ;
	      mrnaDesignExtendDown (sc, ds, segs, ks, path, vp->nn, nIntron, nStart, nStop, useCDS, bestScore, maxCover, flag) ;
	      nStart = (vp->type & gS) ? 1 : 0 ;
	      nStop = 1 ;
	      mrnaDesignExtendUp (sc, ds, segs, ks, path, vp->nn, nIntron, nStart, nStop, useCDS, bestScore, maxCover, flag) ;
	      if (mrnaDesignIsNewPath (segs, ksPaths, ks, path, tested, xStart, xStop, ds->h))
		mrnaDesignExport (s2m, sc, ds, segs, ks, path, smrnas) ;
	      else
		path-- ;
	    }
	  if(useCDS) /* we only want the best CDS guy */
	    break ;
	}

  return path ;
} /* mrnaDesignFindPaths */

/**********************************************************************************/
/**********************************************************************************/

/* s2m->>compositeDesignCovering was created by  mrnaDesignUsingCompositeStrategy below
 * the final mRNAs have been contructed by makemrna.c
 * the task is to reintrocuce the SL. polyA etc flags
 * and the quantitative covering of all elements
 */
static void mrnaDesignSetOneCompletenessFlag (S2M *s2m, SC* sc, HIT *up, SMRNA *smrna, Array estHits, BOOL isTop)
{
  Array covering = s2m->compositeDesignCovering ;
  int iss, issMax = covering ? arrayMax (covering) : 0  ;
  DSX *ssp ;
  unsigned int topFlags = gS | gS0 | gReal5p ;
  unsigned int bottomFlags = gReal3p | gA ;
  unsigned int cFlags = isTop ? topFlags : bottomFlags ;
  int a1 = up->a1, a2 = up->a2 ;
  if (!(up->type & gX) || !(up->type & cFlags))
    return ;
  up->type &= ~cFlags ;

  if (issMax)
    for (iss = 0, ssp = arrp (covering, 0, DSX) ; iss < issMax ; iss++, ssp++)
      {
	if (ssp->a1 > a2)
	  break ;
	if (ssp->a2 < a1)
	  continue ;
	if (isTop && ssp->a1 >= a1 && ssp->a1 <= a1 + 10)
	  {
	    up->type |= (ssp->type & cFlags) ;
	  }
	if (!isTop && ssp->a2 >= a2 - 10 && ssp->a2 <= a2)
	  {
	    up->type |= (ssp->type & cFlags) ;
	  }
      }

  return ;
} /*  mrnaDesignSetOneCompletenessFlag */

/**********************************************************************************/

static void mrnaDesignSetMrnas (S2M *s2m, SC* sc, DS *ss, SMRNA *gmrna, Array smrnas)
{
  HIT *up ;
  int iii, j ;
  SMRNA *smrna ;

  return ;
  
  for (iii = 0 ; iii >-1 && iii < arrayMax(smrnas) ; iii++) 
    {
      smrna = arrp (smrnas, iii, SMRNA) ;
      if (smrna->hits && arrayMax (smrna->hits))      
	{
	  /* search the flags of the first exon */
	  for (j = 0, up = arrp (smrna->hits, 0, HIT) ; j < 1 && j < arrayMax(smrna->hits);  up++, j++)
	    mrnaDesignSetOneCompletenessFlag (s2m, sc, up, smrna, gmrna->estHits, TRUE) ;
	  
	  /* search the polyA flags of the last exon */
	  for (j = arrayMax(smrna->hits) - 1 ; j >= 0 && (up = arrp (smrna->hits, j, HIT)) ; j = -1)
	    mrnaDesignSetOneCompletenessFlag (s2m, sc, up, smrna, gmrna->estHits, FALSE) ;
	}
    }
  
  return ;
} /*  mrnaDesignSetMrnas */

/*********************************************************************/
/*********************************************************************/
/* in each mrna, already annotated for completeness, annotate the ORFs */

void mrnaDesignSavePepAndKantor (KEY product, Array pep)
{
  int checkSum = hashArray (pep) ;
  char buf[64] ;
  KEY kantor = 0 ;
  OBJ Product = 0, Kantor = 0 ;
  
  if (checkSum > 0)
    sprintf (buf, "P%d", checkSum) ;
  else
    sprintf (buf, "N%d", - checkSum) ;
  
  lexaddkey (buf, &kantor, _VKantor) ;
  if ((Product = bsUpdate (product)))
    {
      bsAddKey (Product, _Kantor, kantor) ;
      bsSave (Product) ;
    }
  if ((Kantor = bsUpdate(kantor)))
    {
      bsAddKey (Kantor, _Representative_product, product) ;
      bsSave (Kantor) ;
    }
  peptideStore (kantor, pep) ;
  peptideStore (product, pep) ;
} /* mrnaDesignSavePepAndKantor */

/**********************************************************************************/

static int mrnaDesignOrfQualityOrder (const void *a, const void *b)
{
  const ORFS *sa = (const ORFS*) a, *sb = (const ORFS*) b ;

  int n ;
  n = sa->quality - sb->quality ; if (n) return -n ;
  n = sa->p1 - sb->p1 ; if (n) return n ;
  return 0 ;
} /* mrnaDesignOrfQualityOrder */

/*********************************************************************/
/* notice the inversion function wants single, the tag is multiple */
static int mrnaPleaseUseLeucine (void)
{
  static int style = -1 ;

  if (style == -1)
    {
      KEYSET ks = query (0, "Find clone main_Clone && UseLeucine") ;
      style = keySetMax (ks) ? 1 : 0 ;
      
      keySetDestroy (ks) ;
    }
  if (style == -1)
    style = 0 ; /* default */

  return style ;
}

/**********************************************************************************/

static void mrnaDesignSetOneMrnaOrf (S2M *s2m, SC* sc, DS *ss, SMRNA *gmrna, Array smrnas, SMRNA *smrna)
{
  AC_HANDLE h = ac_new_handle () ;
  Array gDna, mDna, introns, orfs, pep ;
  int ii, jj, iMax = arrayMax (smrna->hits), dnaMax = 0, nOrf = 0, nI = 0 ;
  int uOrf = -1, uOrfStart = 0  ;
  KEYSET x2a ;
  int offset = 0 ;
  BOOL useLeu = mrnaPleaseUseLeucine () ;
  HIT *up, *up2 ;
  ORFS *orf ;
  int bestOrf = -1, bestOrfLn = 0, maxStop = 0 ;
  char *cp, ntgType[8] ;
  char *translationTable = pepGetTranslationTable (s2m->cosmid1 ? s2m->cosmid1 : s2m->cosmid, 0) ;
  int dnaRab = 9 ;
  strcpy (ntgType, "...NNN.") ;
    
  /* grab the genomic DMNA */
  if (sc->isUp)
    { gDna = s2m->dnaR ; offset = arrayMax (gDna) - sc->a1 + smrna->a1 - 1 ; }
  else
    { gDna = s2m->dnaD  ; offset = sc->a1 + smrna->a1  - 2 ; }
  
  /* extract the spliced DNA */
  up2 = arrp (smrna->hits, iMax - 1, HIT) ;
  smrna->orfs = orfs = arrayHandleCreate (32, ORFS, ss->h) ;
  introns = arrayHandleCreate (32, int, ss->h) ;
  smrna->mdna = mDna = arrayHandleCreate (up2->a2 + dnaRab + 2, char, ss->h) ;
  x2a = smrna->x2a = keySetHandleCreate (ss->h) ;
  cp = arrayp (mDna, 0, char) ;
  for (ii = dnaMax = 0, jj = 1, up = arrp (smrna->hits, 0, HIT) ; ii < iMax ; ii++, up++)
    {
      if (up->type & gX)
	{
	  int ln = up->a2 - up->a1 + 1 ;
	  for (int i = 0 ; i < ln + dnaRab + 30 ; i++)
	    keySet (x2a, jj + i) = up->a1 + i ;
	  jj += ln ;
	  if (up->type & (gS | gReal5p))
	    smrna->isNH2Complete = TRUE ;
	  if (up->type & (gA | gReal3p))
	    smrna->isCOOHComplete = TRUE ;
	  if (up->a2 > up2->a2)
	    messcrash ("Bad coordinates in mrnaDesignSetOneMrnaOrfs") ;
	  memcpy (cp, arrp (gDna, offset + up->a1 - 1, char), ln + dnaRab) ;
	  dnaMax += ln ;  /* the 3 extra base allow for a stop and will be overwritten */
	  cp += ln ;
	}
      else if (up->type & gI)
	{
	  array (introns, nI, int) = dnaMax ;
	  nI++ ;
	}
    }
  cp += dnaRab ; *cp++ = 0 ; *cp++ = 0 ; /* double zero terminate */
  dnaMax += dnaRab ;
  arrayMax (mDna) = dnaMax ;
  /* translate in all 3 frames and search the longet ORFS */
  pep = arrayHandleCreate (1024, char, h) ;
  for (int frame = 0 ; frame < 3 ; frame++)
    {
      int iBase, start = -3 + frame, lastStop = -3 + frame, firstATG = -1, iPep = 0, leu = -1 , dx = 0 ;
      if (! smrna->isNH2Complete)
	start = frame ;
      BOOL isLeu = FALSE ;
      BOOL isMet = FALSE ;
      for (iBase = frame, cp = arrp (mDna, iBase, char) ; iBase <= dnaMax  ; iBase += 3, cp += 3)
	{
	  char tt = 0 ; /* not a peptide code */
	  if (iBase < dnaMax - 3)
	    {
	      tt = e_codon (cp, translationTable) ;
	      array (pep, iPep++, char) = tt ;
	    }
	  switch ((int) tt)
	    {
	    case 'M': /* start on Met */
	      if (firstATG == -1)
		firstATG = iBase ;
	      if (start < 0)
		{ start = iBase ; isMet = TRUE ; }
	      else if (isLeu && start + MINI_LEU2MET > iBase)
		{ start = iBase ; isLeu = FALSE ; isMet = TRUE ;}
	      break ;
	    case 'L': /* start on Leu */
	      if (useLeu &&
		  start < 0 &&
		  leu == -1 &&
		  *(cp+2) == G_ && 
                  ((*(cp-3) == G_ && *(cp+3) == G_) || (*(cp-3) == A_)) &&
                  (*(cp) == A_ || *(cp+1) == T_)
                  )
		{
                  ntgType[0] = ace_lower (dnaDecodeChar[(int) *(cp-3)]) ;
                  ntgType[1] = '.' ;
                  ntgType[2] = '.' ;
                  ntgType[3] = ace_upper (dnaDecodeChar[(int) *cp]) ;
                  ntgType[4] = ace_upper (dnaDecodeChar[(int) *(cp+1)]) ;
                  ntgType[5] = ace_upper (dnaDecodeChar[(int) *(cp+2)]) ;
                  ntgType[6] = *(cp+3) == G_ ? 'g' : '.' ;
                  ntgType[7] = 0 ;
		  
		  start = iBase ; isLeu = TRUE ; 
		  leu = iBase ;
		}
	      break ;
	    case 0:
	      if (0 && smrna->isCOOHComplete) /* do not register a 3' open ORF */
		break ;
	      dx = 0 ;
	      /* fall thru */
	    case '*': /* stop */
	      if (isLeu && iBase - start < MINI_LEU2END)
		{
		  isLeu = FALSE ; start = firstATG ;
		  isMet = firstATG > 0 ? TRUE : FALSE ;
		}
	      if (start >= 0 && iBase > start + 90)
		{ /* register */
		  orf = arrayp (orfs, nOrf++, ORFS) ;
		  orf->frame = frame + 1 ;
		  orf->isLeu = isLeu ;
		  orf->leu = leu ;
		  orf->hasStop = (tt == '*' ? TRUE : FALSE) ;
		  orf->hasStart = isLeu || isMet ;
		  orf->smrna = smrna ;
		  array (pep, iPep++, char) = 0 ;
		  if (start < uOrfStart)
		    { uOrfStart = start ; uOrf = nOrf - 1 ; }
		  orf->p1 = start + 1 ;  /* human coordinates in the mRNA */
		  orf->p2 = iBase + 3 ;
		  orf->m1 = arr (x2a, orf->p1, int) ;
		  orf->m2 = arr (x2a, orf->p2, int) ;
		  if (lastStop >= 0)
		    orf->upStop = lastStop ;
		  int dk = (start - lastStop - 3)/3 ;
		  array (pep, iPep++,char) = 0 ;
		  if (dk)
		    for (int i = dk ; i < iPep ; i++)
		      arr (pep, i-dk, char) = arr (pep, i, char) ;
		  arrayMax (pep) = iPep - 1 - dk + (tt == '*' ? -1 : 0) ;
		  pepEncodeArray (pep) ;
		  orf->pep = arrayHandleCopy (pep, ss->h) ;
		  orf->codingLn = orf->p2 - orf->p1 + 1 + dx - 3 ;
		  orf->openLn = orf->p2 - lastStop - 3 + dx - 3 ;
		  if (firstATG >= 0) orf->firstATG = 1 + (firstATG - start)/3 ;
		  if (smrna->isCOOHComplete && ! orf->hasStop)
		    {

		      if (orf->p2 - orf->p1 + 1 > bestOrfLn + 100 )
			{
			  bestOrfLn = orf->p2 - orf->p1 + 1 ; bestOrf = nOrf - 1 ;
			}
		    }
		  else
		    {
		      if (orf->p2 - orf->p1 + 1 > bestOrfLn)
			{
			  bestOrfLn = orf->p2 - orf->p1 + 1 ; bestOrf = nOrf - 1 ;
			}
		    }
		  orf->weight = pepWeight (orf->pep)/1000.0 ; /* kDalton */
		  orf->pI = pepPI (orf->pep) ;
		  orf->quality = orf->codingLn ;
		  for (int i = 0 ; i < arrayMax (introns) ; i++)
		    {
		      int x = arr (introns, i, int) ;
		      if (x > orf->p1 && x < orf->p2)
			orf->nIntronsInsideCDS++ ;
		      else
			orf->nIntronsOutsideCDS++ ;
		    }
		}
	      start = -3 + frame ; /* behind a stop i must find a Met */
	      firstATG = -1 ;
	      lastStop = iBase ;
	      if (bestOrf == nOrf && lastStop > maxStop)
		maxStop = lastStop ;
	      isLeu = FALSE ;
	      isMet = FALSE ;
	      leu = -1 ;
	      iPep = 0 ;
	      dx = 0 ;
	      pep = arrayHandleCreate (1024, char, h) ;
	      break ;
	    }
	}
    }
  if (dnaRab)
    {  /* reclip */
      int dx = dnaMax - maxStop - 3 ;
      if (dx < 0) dx = 0 ;
      if (! maxStop)
	;
      else if (dx >= dnaRab) /* reclip the rab */
	arrayMax (mDna) = dnaMax - dnaRab ; 
      else
	{
	  int dy = dnaRab - dx ; /* needed Rab */
	  arrayMax (mDna) = dnaMax - dx ;
	  up = arrp (smrna->hits, iMax -1, HIT) ;
	  up->a2 += dy ; /* extend the mrna exon */
	  if (smrna->pA)
	    {
	      for (int i = 0 ; i < arrayMax (smrna->pA) ; i++)
		{
		  HIT *pA = arrayp (smrna->pA, i, HIT) ;
		  if (pA->type == gReal3p || pA->type == gA)
		    pA->a2 -= dy ;
		}
	    }
	}
    }
  
  if (bestOrf >= 0)
    arrp(orfs, bestOrf, ORFS)->isBest = TRUE ;
  if (uOrf >= 0 && bestOrfLn > 360 && uOrf != bestOrf)
    arrp(orfs, uOrf, ORFS)->isUorf = TRUE ;
  
  ac_free (h) ;
  return ;
} /*  mrnaDesignSetOneMrnaOrfs */

/*********************************************************************/

void mrnaDesignSaveOneMrnaCoding (SMRNA *mm, OBJ Mrna) 
{
  AC_HANDLE h = ac_new_handle () ;
  Array hits = arrayHandleCreate (32, DSX, h) ;
  KEYSET ks = keySetHandleCreate (h) ;
  ORFS *orf ;   DSX *up ;
  int ii, jj, iMax = arrayMax (mm->orfs) ;
  const KEY _Coding = str2tag ("Coding") ;
  KEYSET x2a = mm->x2a ;
  
  /* collect orf extremities */
  for (ii = jj = 0 ; ii < iMax ; ii++)
    {
      orf = arrp (mm->orfs, ii, ORFS) ;
      keySet (ks, jj++) = orf->p1 ;
      keySet (ks, jj++) = orf->p2 + 1 ;
    }
  keySet (ks, jj++) = arrayMax (mm->mdna) + 1 ;
  /* add the intron boundaries */
  for (int i = 0, x = 0 ; i < keySetMax (x2a) ; i++)
    {
      int y = keySet(x2a, i) ;
      if (y > x + 1)
	keySet (ks, jj++) = i ;
      x = y ;
    }
    
  keySetSort (ks) ;
  keySetCompress (ks) ;
  /* create subsegs, with flags */
  for (ii = jj = 0 ; ii < iMax && ii < 25 ; ii++)
    {
      int x = 1 ;
      orf = arrp (mm->orfs, ii, ORFS) ;
      for (int j = 0 ; j < keySetMax (ks) ; j++)
	{
	  int y = keySet (ks, j) ;
	  up = arrayp (hits, jj++, DSX) ;
	  up->a1 = x ; up->a2 = y - 1 ;
	  if (y <= orf->p1)
	    up->donor = 1 ;
	  else if (x >= orf->p2)
	    up->end = 1 ;
	  up->type = ii ;
	  x = y ;
	}
    }
  /* sort merge export */
  arraySort (hits, dsA1Order) ; 

  if (arrayMax (hits))
    for (ii = 0 ; ii < arrayMax (hits) ; ii++) 
      {
	up = arrp (hits, ii, DSX) ;
	int x1 = up->a1, x2 = up->a2 ;
	int a1 = keySet(x2a, x1) ;
	int a2 = keySet(x2a, x2) ;
	char flag[128] ;
	int k = 0 ;
	bsAddData (Mrna, _Coding, _Int, &a1) ;
	bsAddData (Mrna, _bsRight, _Int, &a2) ;
	bsAddData (Mrna, _bsRight, _Int, &x1) ;
	bsAddData (Mrna, _bsRight, _Int, &x2) ;
	for (ii = ii ;  ii < arrayMax (hits) ; ii++, up++)
	  {
	    if (up->a1 == x1 && up->a2 == x2)
	      {
		if (up->donor)
		  flag[k++] = '5' ;
		if (up->end)
		  flag[k++] = '3' ;
		flag[k++] = 'A' + up->type ;
		flag[k++] = ' ' ;
	      }
	    else
	      { ii-- ; up-- ; break ; }
	  }
	flag[k++] = 0 ;
	bsAddData (Mrna, _bsRight, _Text, flag) ;
      }
  else
    {
      HIT *up ;
      int ii, iMax = arrayMax (mm->hits) ;
      for (ii = 0, up = arrp (mm->hits, 0, HIT) ; ii < iMax ; ii++, up++)
	{
	  if (up->type & gX)
	    {
	      bsAddData (Mrna, _Coding, _Int, &up->a1) ;
	      bsAddData (Mrna, _bsRight, _Int, &up->a2) ;
	      bsAddData (Mrna, _bsRight, _Int, &up->x1) ;
	      bsAddData (Mrna, _bsRight, _Int, &up->x2) ;
	      bsAddData (Mrna, _bsRight, _Text, "3A") ;
	    }
	}
    }
  ac_free (h) ;
  return ;
} /* mrnaDesignSaveOneMrnaCoding */

/**********************************************************************************/

static void mrnaDesignSaveOneMrnaDna (SMRNA *mm)
{
  AC_HANDLE h = ac_new_handle () ;
  KEY mrna = mm->tr ;
  const char *mrnaName = name (mrna) ;

  if (mm->mdna)
    {
      Array dna = arrayHandleCopy (mm->mdna, h) ;
      KEY dnaKey = 0 ;
      int dnaMax = arrayMax (mm->mdna) ;
      char *dnaName = hprintf(h, "mRNA:%s", mrnaName) ;
      lexaddkey (dnaName, &dnaKey, _VDNA) ;
      OBJ Mrna = bsUpdate (mm->tr) ;
      bsAddKey (Mrna, _DNA, dnaKey) ;
      bsAddData (Mrna, _bsRight, _Int, &dnaMax) ;
      bsSave (Mrna) ;
      dnaStoreDestroy (dnaKey, dna) ; /* calls bsUpdate(seqObj(dnaKey)) */
    }
  ac_free (h) ;
} /* mrnaDesignSaveOneMrnaDna */

/**********************************************************************************/
/* export introns, in TG or in Mrna */
static void mrnaDesignSaveOneIntrons (KEYSET introns, OBJ obj)
{
  if (introns && keySetMax (introns))
    {
      KEYSET ks ;
      int k ;
      
      arraySort (introns, keySetAlphaOrder) ;
      keySetCompress (introns) ;
      for (int ii = 0 ;  ii < keySetMax (introns) ; ii++)
	bsAddKey (obj, _Intron, keySet (introns, ii)) ;
      ks = query (introns, "gt_ag") ;
      k = keySetMax (ks) ;
      if (k)
	bsAddData (obj, _gt_ag, _Int, &k) ;
      keySetDestroy (ks) ;

      ks = query (introns, "gc_ag") ;
      k = keySetMax (ks) ;
      if (k)
	bsAddData (obj, _gc_ag, _Int, &k) ;
      keySetDestroy (ks) ;

      ks = query (introns, "ct_ac") ;
      k = keySetMax (ks) ;
      if (k)
	bsAddData (obj, _ct_ac, _Int, &k) ;
      keySetDestroy (ks) ;

      ks = query (introns, "at_ac") ;
      k = keySetMax (ks) ;
      if (k)
	bsAddData (obj, str2tag("at_ac"), _Int, &k) ;
      keySetDestroy (ks) ;

      ks = query (introns, "Other") ;
      k = keySetMax (ks) ;
      if (k)
	{
	  bsAddData (obj, _Other, _Text, "nn_nn") ;
	  bsAddData (obj, _bsRight, _Int, &k) ;
	}
      keySetDestroy (ks) ;
    }
} /* mrnaDesignSaveOneIntrons */

/**********************************************************************************/
/* set in each genes the tags  clones reads assembled_from */
static void mrnaDesignSaveOneEstHits (Array estHits, OBJ obj)
{
  int ii, iiMax = arrayMax ( estHits) ;
  DSX *up ;
  
  /* Read cDNA_clone */
  if (iiMax)
    {
      arraySort (estHits, dsEstOrder) ;
      for (ii = 0, up = arrp (estHits, ii, DSX) ; ii < iiMax ; ii++, up++)
	{
	  bsAddKey (obj, _Read, up->est) ;
	  bsAddData (obj, _bsRight, _Int, &(up->x1)) ;
	  bsAddData (obj, _bsRight, _Int, &(up->x2)) ;
	  bsAddData (obj, _bsRight, _Int, &(up->ln)) ;
	  bsAddData (obj, _bsRight, _Int, &(up->score)) ;
	  bsAddData (obj, _bsRight, _Int, &(up->x1)) ;
	  bsAddData (obj, _bsRight, _Int, &(up->x2)) ;
	  
	  bsAddKey (obj, _cDNA_clone, up->cDNA_clone) ;
	}
      
      /* Assembled_from */
      arraySort (estHits, dsA1Order) ;
      for (ii = 0, up = arrp (estHits, ii, DSX) ; ii < iiMax ; ii++, up++)
	{
	  int k = 0 ;      
	  bsAddData (obj, _Assembled_from, _Int, &up->a1) ;
	  bsAddData (obj, _bsRight, _Int, &up->a2) ;
	  bsAddKey (obj, _bsRight, up->est) ;
	  bsAddData (obj, _bsRight, _Int, &up->x1) ;
	  bsAddData (obj, _bsRight, _Int, &up->x2) ;
	  bsAddData (obj, _bsRight, _Int, &k) ;
	}
    }
} /* mrnaDesignSaveOneEstHits */

/**********************************************************************************/

static void mrnaDesignSaveOneMrnaSplicing (SMRNA *mm, OBJ Mrna)
{
  /* splicing */
  if (mm->hits)
    {
      HIT *up ;
      int i, a0, a1, a2, b1 = 1, b2, nX = 0, da, ngt_ag = 0 ;
      int iMax = arrayMax (mm->hits) ;

      bsAddTag (Mrna, str2tag("Best_in_gene")) ;
      for (i = 0, up = arrp (mm->hits, i, HIT), a0 = up->a1 ; i < iMax ; up++, i++)
	{
	  if (up->clipTop)
	  {
	    int k = 1 ;
	    bsAddData (Mrna, str2tag ("Aggregated_5p_clones"), _Int, &up->clipTop) ;
	    bsAddData (Mrna, str2tag("Valid5p"), _Int, &k) ;
	    bsAddData (Mrna, _bsRight, _Int, &k) ;
	    bsAddData (Mrna, _bsRight, _Int, &up->clipTop) ;
	  }
	  
	  a1 = up->a1 - a0 + 1 ;
	  a2 = up->a2 - a0 + 1 ;
	  da = a2 - a1 + 1 ;
	  if (up->type & gX)
	    b2 = b1 + da - 1 ;
	  else
	    { b1 = b1 - 1 ; b2 = b1 + 1 ; }
	  bsAddData (Mrna, _Splicing, _Int, &a1) ;
	  bsAddData (Mrna, _bsRight, _Int, &a2) ;
	  bsAddData (Mrna, _bsRight, _Int, &b1) ;
	  bsAddData (Mrna, _bsRight, _Int, &b2) ;
	  if (up->type & gX)
	    {
	      char *cp = messprintf ("%d", ++nX) ;
	      bsAddKey (Mrna, _bsRight, _Exon) ;
	      bsAddData (Mrna, _bsRight, _Text, cp) ;
	      b1 = b2 + 1 ;
	    }
	  else
	    {
	      bsAddKey (Mrna, _bsRight, _Intron) ;
	      bsAddData (Mrna, _bsRight, _Text, "gt_ag") ;
	      ngt_ag++ ;
	      b1 = b2 ;
	    }
	  bsAddData (Mrna, _bsRight, _Text, "Length") ;
	  bsAddData (Mrna, _bsRight, _Int, &da) ;
	  bsAddData (Mrna, _bsRight, _Text, "bp") ;
	  bsAddData (Mrna, _bsRight, _Int, &up->zone) ;
	  bsAddData (Mrna, _bsRight, _Text, "cover") ;
	}
      bsAddData (Mrna, str2tag ("Nb_exons"), _Int, &nX) ;
      if (mm->pA)
	for (int i = 0 ; i < arrayMax (mm->pA) ; i++)
	  {
	    HIT *pA = arrp (mm->pA, i, HIT) ;
	    if (pA->type & gA)
	      {
		bsAddData (Mrna, str2tag("PolyA_found"), _Int, &pA->a2) ;
		if (pA->cDNA_clone)
		  {
		    bsAddKey (Mrna, _bsRight, pA->cDNA_clone) ;
		    bsAddData (Mrna, _bsRight, _Int, &pA->zone) ;
		  }
	      }
	    if ((pA->type & gS) && pA->a1 == 1)
	      {
		bsAddKey (Mrna, str2tag ("Transpliced_to"), pA->est) ;
		if (pA->cDNA_clone)
		  {
		    bsAddKey (Mrna, _bsRight, pA->cDNA_clone) ;
		    bsAddData (Mrna, _bsRight, _Int, &pA->zone) ;
		  }
	      }
	    if (pA->type & gReal3p)
	      {
		bsAddData (Mrna, str2tag ("Aggregated_3p_clones"), _Int, &pA->zone) ;
		bsAddData (Mrna, _bsRight, _Int, &pA->a2) ;
	      }
	  }
      
      mrnaDesignSaveOneEstHits (mm->estHits, Mrna) ;
      mrnaDesignSaveOneIntrons (mm->introns, Mrna) ;
    }
} /* mrnaDesignSaveOneMrnaSplicing */

/**********************************************************************************/

static void mrnaDesignSaveOneMrnaEstHits (SMRNA *mm, OBJ Mrna)
{
  /* clones */
  if (mm->estHits)
    {
      int ii, iMax = arrayMax (mm->estHits) ;
      DSX *up ;
      for (ii = 0, up = arrp (mm->estHits, 0, DSX) ; ii < iMax ; ii++, up++)
	bsAddKey (Mrna, _cDNA_clone, up->cDNA_clone) ;
    }
} /* mrnaDesignSaveOneMrnaEstHits */

/**********************************************************************************/

static void mrnaDesignSetOneMrnaEstHits (S2M *s2m, SC* sc, DS *ss, SMRNA *mm)
{
    Array plainHits = s2m->plainHits ;
  Array estHits ;
  HIT *up ;
  DSX *dsxp ;
  int g1 = mm->a1, g2 = mm->a2 ;
  int ii, jj, iMax = arrayMax (plainHits) ;
  estHits = mm->estHits = arrayHandleCreate (128, DSX, ss->h) ;

  return ; /* ATTTENTION je ne sais pas qui sont les reads du transcript */
  
  // Read: x1 x2 (read coord) length err qual longet exact
  // Assembled_from a1 a1 (tg coords) Est x1 x2 nerr
  for (ii = jj = 0, up = arrp (plainHits, 0, HIT) ; ii < iMax ; up++, ii++)
    {
      if (up->a1 < g2 && up->a2 > g1)
	{
	  dsxp = arrayp (estHits, jj++, DSX) ;
	  dsxp->est = up->est;
	  dsxp->cDNA_clone = up->cDNA_clone ;
	  dsxp->a1 = g1 > up->a1 ? 1 : up->a1 - g1 + 1 ;
	  dsxp->a2 = g2 < up->a2 ? g2 - g1 + 1 : up->a2 - g1 + 1 ;
	  dsxp->a1 = up->a1 - g1 + 1 ;
	  dsxp->a2 = up->a2 - g1 + 1 ;
	  if (up->x1 < up->x2)
	    {
	      dsxp->x1 = g1 > up->a1 ? up->x1 + g1 - up->a1 : up->x1 ;
	      dsxp->x2 = g2 < up->a2 ? up->x2 + g2 - up->a2 : up->x2 ;
	    }
	  else
	    {
	      dsxp->x1 = g1 > up->a1 ? up->x1 - g1 + up->a1 : up->x1 ;
	      dsxp->x2 = g2 < up->a2 ? up->x2 - g2 + up->a2 : up->x2 ;
	    }
	  dsxp->x1 = up->x1 ;
	  dsxp->x2 = up->x2 ;
	}
    }
  return ; 
} /* mrnaDesignSetOneMrnaEstHits */

/**********************************************************************************/
/* set in each genes the tags  clones reads assembled_from */
static void mrnaDesignSetOneGeneEstHits (S2M *s2m, SC* sc, DS *ss, SMRNA *tgp)
{
  Array plainHits = s2m->plainHits ;
  Array estHits ;
  HIT *up ;
  DSX *dsxp ;
  BOOL isDown = TRUE ;
  int a1 = tgp->a1, a2 = tgp->a2 ;
  int dx, dy ;
  int ii, jj, iMax = arrayMax (plainHits) ;
  estHits = tgp->estHits = arrayHandleCreate (128, DSX, ss->h) ;
  // Read: x1 x2 (read coord) length err qual longet exact
  // Assembled_from a1 a1 (tg coords) Est x1 x2 nerr
  if (a1 > a2) {a1 = tgp->a2 ; a2 = tgp->a1 ; isDown = FALSE ; }
  for (ii = jj = 0, up = arrp (plainHits, 0, HIT) ; ii < iMax ; up++, ii++)
    {
      if ((up->type & gX) && up->a1 < a2 && up->a2 > a1)
	{
	  dx = dy = 0 ;
	  dsxp = arrayp (estHits, jj++, DSX) ;
	  dsxp->est = up->est;
	  if (!strncmp (name (up->est), "XI_", 3))
	    {
	      if (up->x1 < 20) dx = up->x1 - 1 ; 
	      else dy = 16 ;
	      }
	  if (!strncmp (name (up->est), "XI_", 3))
	    {
	      if (up->x1 < 20) dx = up->x1 - 1 ; 
	      else if (up->x2 < 16) dy = up->x2 - 1 ;
	    }
	  dsxp->cDNA_clone = up->cDNA_clone ;
	  if (isDown)
	    {
	      if (up->a1 < up->a2)
		{
		  dsxp->a1 = up->a1 - a1 + 1 - dx ;
		  dsxp->a2 = up->a2 - a1 + 1 + dy ;
		  dsxp->x1 = up->x1 - dx ;
		  dsxp->x2 = up->x2 + dy ;
		}
	      else
		{
		  messerror ( "bad orientation in %s", name(up->est)) ;
		  dsxp->a1 = up->a2 - a1 + 1 ;
		  dsxp->a2 = up->a1 - a1 + 1 ;
		  dsxp->x1 = up->x2 ;
		  dsxp->x2 = up->x1 ;
		}
	    }
	  else
	    {
	      if (up->a1 > up->a2)
		{
		  messerror ( "bad orientation in %s", name(up->est)) ;
		  dsxp->a1 = - up->a1 + a2 + 1 ;
		  dsxp->a2 = - up->a2 + a2 + 1 ;
		  dsxp->x1 = up->x1 ;
		  dsxp->x2 = up->x2 ;
		}
	      else
		{
		  dsxp->a1 = - up->a2 + a2 + 1 - dx ;
		  dsxp->a2 = - up->a1 + a2 + 1 + dy ;
		  dsxp->x1 = up->x2 - dx ;
		  dsxp->x2 = up->x1 + dx ;
		}
	    }
	}
    }
  return ; 
} /* mrnaDesignSetOneGeneEsmatHits */

/**********************************************************************************/

static void mrnaDesignSaveOneOrf (SMRNA *mm, ORFS *orf, int iOrf, OBJ Mrna)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *mrnaName = name (mm->tr) ;
  KEY methodProduct ;
  KEY product = 0 ;
  OBJ Product = 0 ;
  const char *productName = iOrf ? hprintf(h, "%s%d", mrnaName, iOrf) : mrnaName ;
  int dnaMax = arrayMax (mm->mdna) ;
  
  lexaddkey ("PRODUCT", &methodProduct, _VMethod) ;
  lexaddkey (productName, &product, _VmProduct) ;
  mrnaDesignSavePepAndKantor (product, orf->pep) ;
  if ((Product = bsUpdate (product)))
    {
      bsAddData (Product, str2tag("Molecular_weight"), _Float, &orf->weight) ;
      bsAddData (Product, _bsRight, _Float, &orf->pI) ;
      bsAddData (Product, str2tag("Expasy"), _Float, &orf->weight) ;
      bsAddData (Product, _bsRight, _Float, &orf->pI) ;
      if (mm->gene) bsAddKey (Product, _From_gene, mm->gene) ;
      if (orf->hasStart && orf->hasStop)
	bsAddTag (Product, _Complete) ;
      if (mm->isNH2Complete)
	bsAddTag (Product, str2tag ("mRNA_5p_complete")) ;
      if (orf->hasStart)
	{
	  bsAddTag (Product, str2tag ("at_position_1")) ;
	  if (orf->firstATG > 0)
	    bsAddData (Product, str2tag ("First_ATG"), _Int, &orf->firstATG) ;
	  if (orf->upStop)
	    { int k = orf->p1 - orf->upStop - 1 ; bsAddData (Product, str2tag ("Up_stop"), _Int, &k) ; }
	  if (orf->leu >= 0)
	    {
	      int n = (orf->p1 - orf->leu)/3 ;
	      if (n < 0) n -= 1 ; else n++ ;
	      if (n <= 1)
		{
		  bsAddData (Product, _First_Kozak, _Int, &n) ;
		  bsAddData (Product, _bsRight, _Text, orf->ntgType) ;
		}
	    }		  
	}
      int k ;
      if (orf->hasStop)
	{
	  k = orf->p2 - orf->p1 - 1 ;
	  bsAddData (Product, str2tag ("Down_stop"), _Int, &k) ;
	}
      if (orf->isBest)
	bsAddTag (Product, str2tag ("Best_product")) ;
      if (orf->codingLn >= 300)
	bsAddTag (Product, str2tag ("Good_product")) ;
      if (orf->codingLn >= 600)
	bsAddTag (Product, str2tag ("Very_good_product")) ;
      if (orf->isUorf)
	bsAddTag (Product, str2tag ("uORF_candidate")) ;
      bsAddData (Product, str2tag ("Frame"), _Int, &orf->frame) ;
      bsAddData (Product, str2tag ("Open_length"), _Int, &orf->openLn) ;
      bsAddData (Product, _bsRight, _Text, "bp") ;
      k = orf->openLn/3 ;
      bsAddData (Product, _bsRight, _Int, &k) ;
      bsAddData (Product, _bsRight, _Text, "aa") ;
      if (orf->p1 > 1) { orf->p1-- ; bsAddData (Product, str2tag ("Length_5prime_UTR"), _Int, &orf->p1) ; orf->p1++ ; }
      bsAddData (Product, str2tag ("Coding_length"), _Int, &orf->codingLn) ;
      bsAddData (Product, _bsRight, _Text, "bp") ;
      k = orf->codingLn/3 ;
      bsAddData (Product, _bsRight, _Int, &k) ;
      bsAddData (Product, _bsRight, _Text, "aa") ;
      
      if (orf->p2 < dnaMax) { int k = dnaMax - orf->p2 ; bsAddData (Product, str2tag ("Length_3prime_UTR"), _Int, &k) ; }
      if (orf->nIntronsOutsideCDS) bsAddData (Product, str2tag ("Nb_Introns_outside_CDS"), _Int, &orf->nIntronsOutsideCDS) ;
      if (orf->nIntronsInsideCDS) bsAddData (Product, str2tag ("Nb_Introns_in_CDS"), _Int, &orf->nIntronsInsideCDS) ;	    
      if (0) bsAddKey (Product, _mRNA, mm->tr) ;
      
      {
	int k = 1 ;
	bsAddData (Product, str2tag ("Source_exons"), _Int, &k) ;
	k = orf->p2 - orf->p1 + 1 ;
	bsAddData (Product, _bsRight, _Int, &k) ;
	bsAddKey (Product, _bsRight, _Exon) ;
      }
      bsAddKey (Product, _Method, methodProduct) ;
    }	  
  bsSave (Product) ;
  if (1) /* (Mrna = bsUpdate (mrna) */
    {
      bsAddKey (Mrna, _Product, product) ;
      bsAddData (Mrna, _bsRight, _Int, &orf->p1) ;
      bsAddData (Mrna, _bsRight, _Int, &orf->p2) ;
      bsAddData (Mrna, _bsRight, _Int, &orf->m1) ;
      bsAddData (Mrna, _bsRight, _Int, &orf->m2) ;
      /* bsSave (Mrna) ; */
    }
  if (orf->isBest)
    {
      int k ;
      bsAddData (Mrna, str2tag ("Longest_ORF"), _Int, &orf->openLn) ;
      bsAddData (Mrna, _bsRight, _Text, "bp") ;
      k = orf->openLn/3 ; bsAddData (Mrna, _bsRight, _Int, &k) ;
      bsAddData (Mrna, _bsRight, _Text, "aa") ;
      bsAddData (Mrna, str2tag ("Longest_CDS"), _Int, &orf->codingLn) ;
      bsAddData (Mrna, _bsRight, _Text, "bp") ;
      k = orf->codingLn/3 ; bsAddData (Mrna, _bsRight, _Int, &k) ;
      bsAddData (Mrna, _bsRight, _Text, "aa") ;
      k = orf->p1 - 1 ;
      if (k > 0) bsAddData (Mrna, str2tag ("Length_5prime_UTR"), _Int, &k) ;
      k = arrayMax (mm->mdna) - orf->p2 ;
      if (k > 0) bsAddData (Mrna, str2tag ("Length_3prime_UTR"), _Int, &k) ;
    }

  return ;
} /* mrnaDesignSaveOneOrf */

/**********************************************************************************/

static void mrnaDesignSaveOneMrna (SMRNA *tg, SMRNA *mm)
{
  Array orfs = mm->orfs ;
  OBJ Mrna = 0 ;
  mm->composite = TRUE ;
  mrnaDesignSaveOneMrnaDna (mm) ;
  Mrna = bsUpdate (mm->tr) ;
  mrnaDesignSaveOneMrnaSplicing (mm, Mrna) ;
  arraySort (orfs, mrnaDesignOrfQualityOrder)  ;
  mrnaDesignSaveOneMrnaCoding (mm, Mrna) ;
  if (tg) mrnaDesignSaveOneMrnaEstHits (tg, Mrna) ; /* ATTENTION i cheat here */
  if (orfs)
    {
      for (int iOrf = 0 ; iOrf < arrayMax (mm->orfs) ; iOrf++)
	{
	  ORFS *orf = arrayp (orfs, iOrf, ORFS) ;
	  if (orf->isBest || orf->isUorf || orf->codingLn >= 240)
	    mrnaDesignSaveOneOrf (mm, orf, iOrf, Mrna) ;
	}
      arrayDestroy (mm->orfs) ; /* incompatible with makemrna which uses sanother type of orfs */
    }
  bsSave (Mrna) ;
  mrnaAddKantorInfo (mm->tr) ;

  return ;
} /* mrnaDesignSaveMrnas */

/**********************************************************************************/

static void mrnaDesignSetOrfs (S2M *s2m, SC* sc, DS *ss, SMRNA *gmrna, Array smrnas)
{
  for (int ii = 0 ; ii < arrayMax (smrnas) ; ii++)
    {
      SMRNA *mm = arrp (smrnas, ii, SMRNA) ;
      if (mm->hits)
	mrnaDesignSetOneMrnaOrf (s2m, sc, ss, gmrna, smrnas, mm) ;
    }
  return ; 
} /*  mrnaDesignSetOrfs */

/*********************************************************************/
/* get name of smallest Clone Group or smallest clone */
static KEY mrnaDesignNameGeneBySection (SMRNA *tgp)
{
  KEY gg, chrom = tgp->chrom ;
  int g1, g2 ;
  int B = 1 ; /* attention, we must not produce name conflicts */
  int M = 1000000, K = 1000 ;
  int nM = 0 ; /* megabase position */
  int nK = 0 ; /* kilobase position */
  int nB = 0 ; /* B position */

  /* construct g1, g2 the IntMap coords of the tgene */

  g1 = tgp->g1 ; g2 = tgp->g2 ;
      
  if (g1 < g2)
    {
      nM = g2 / M ;
      int dg = g2 - g1 + 1 ;
      nK = (g2 - M *nM) / 1000 ;
      if (dg < 3000)
	{
	  B = 10 ;
	  nB = (g2 - M * nM - K * nK) / B ;
	  if (nB % 2 == 0) nB-- ; /* should be odd */
	}
      else
	{
	  B = 1000 ;
	  if (nK % 2 == 0) nK-- ; /* should be odd */
	}
    }
  else
    {
      nM = g2 / M ;
      int dg = g1 - g2 + 1 ;
      if (dg < 3000)
	{
	  B = 10 ;
	  nK = (g2 - M *nM) / K ;
	  nB = (g2 - M * nM - K * nK + B - 1) / B ; /* +B-1 because we want to be dowmstream of g2 */
	  if (nB % 2 == 1) nB++ ; /* should be even */
	}
      else
	{
	  B = 1000 ;
	  nK = (g2 - M *nM + K - 1) / K ;  /* +K-1 because we want to be dowmstream of g2 */
	  if (nK % 2 == 1) nK++ ; /* should be even */
	}
    }
  if (B < 1000)
    {
      lexaddkey (messprintf ("G_%s__%d.%d.%02d", name(chrom), nM, nK, nB), &gg, _VTranscribed_gene) ;
    }
  else
    {
      lexaddkey (messprintf ("G_%s__%d.%d", name(chrom), nM, nK), &gg, _VTranscribed_gene) ;
    }
    
  return gg ;
} /* mrnaDesignNameGeneBySection */

/*********************************************************************/

static KEY mrnaDesignNameGeneByLocusLink (SMRNA *tgp)
{
  KEY gene = 0, LL = 0 ;

  if (LL)
    {
      int k = 0 ;
      lexaddkey (name(LL), &gene, _VTranscribed_gene) ;
      while (keyFindTag (gene, _Splicing))
	{
	  k++ ;
	  lexaddkey (messprintf("%s.%d", name(LL), k), &gene, _VTranscribed_gene) ;
	}
    }
  else
    return  mrnaDesignNameGeneBySection (tgp) ;
  return gene ;
} /* mrnaDesignNameGeneByLocusLink */

/*********************************************************************/

KEY mrnaDesignGetGeneName (SMRNA *tgp)
{
  int n ;
  KEY cGroup = 0, gene = 0 ;
  int type = mrnaPleaseNameBy () ;

  /* first find a candidate gene name */
  switch (type)
    { 
    case 2:    /* LocusLink name */
      type = 1 ;
      return mrnaDesignNameGeneByLocusLink (tgp) ;        
    case 1:
    default:   /* worm case */
      return mrnaDesignNameGeneBySection (tgp) ;
    }
  
  /* iterate to avoid mixing doubles */
  for (n = 0 ; TRUE ; n++)       /* iterate on gene_name candidate */
    {
      BOOL isNew = TRUE ;
      
      switch (type)
        { 
        case 1:    /* name by section */
          if (n == 0)
            lexaddkey (messprintf("%s", name(cGroup)), &gene, _VTranscribed_gene) ;
          else
            lexaddkey (messprintf("%s_%d", name(cGroup), n), &gene, _VTranscribed_gene) ;
          break ;
        default:   /* worm case */
          if (n == 0)
            {
              if (strncmp("G_",name(cGroup),2))
                lexaddkey (messprintf("G_%s", name(cGroup)), &gene, _VTranscribed_gene) ;
              else
                lexaddkey (messprintf("%s", name(cGroup)), &gene, _VTranscribed_gene) ;
            }
          else
            { 
              if (strncmp("G_",name(cGroup),2))
                lexaddkey (messprintf("G_%s_%d", name(cGroup), n), &gene, _VTranscribed_gene) ;
              else
                lexaddkey (messprintf("%s_%d", name(cGroup), n), &gene, _VTranscribed_gene) ;
            }
          break ;
        }
      if (keyFindTag (gene, _Splicing)) /* name already in use */
        {
          /* printf(" spl:%s(%d/%d) ", name(gene), i, j1) ; */
          isNew = FALSE ;
        } 
      /* printf(" ok:%s(j1=%d) ", name(gene), j1) ; */
      /* success this gene name is good */
      if (isNew)
        break ;
    }
  
  return gene ;
} /* mrnaDesignGetGeneName */

/**********************************************************************************/
/* cluster the mrnas in genes, A, B, A__B and name these genes and hence the mrnnas */
static void mrnaDesignSetGenes (S2M *s2m, SC* sc, DS *ss, SMRNA *gmrna, Array smrnas)
{
  AC_HANDLE h = ac_new_handle () ;
  /* search transitivelly for contacts
   * start with one gene per mrna
   * check for contact, and merge the genes
   *
   * later we will need to recognize the A__B genes
   */

  Array tgs, smrnas2 ;
  SMRNA *smrna, *tgp ;
  int ii, jj, iMax = arrayMax (smrnas) ;
  
  gmrna->tgs = tgs = arrayHandleCreate (12, SMRNA, ss->h) ;
  gmrna->genes = keySetHandleCreate (s2m->h) ;
  /* number the smrnas */
  if (iMax)
    {
      for (ii = 0, smrna = arrp (smrnas, ii, SMRNA) ; ii < iMax ; smrna++, ii++)
	smrna->tr = ii ;
      /* sort them by coordinates */
      smrnas2 = arrayHandleCopy (smrnas, h) ;
      arraySort (smrnas2, smrnaA1Order) ;
      
      /* scan and cluster in a tg the overlapping mrnas */
      int a1 = 0, a2 = 0 ;
      for (ii = 0, jj = -1, a2 = -999, smrna = arrp (smrnas2, ii, SMRNA) ; ii < iMax ; smrna++, ii++)
	{
	  SMRNA *smrna0 = arrp (smrnas, smrna->tr, SMRNA) ;
	  if (!a1) a1 = smrna->a1 ;
	  if (smrna->a1 > a2) jj++ ;
	  if (smrna->a2 > a2) a2 = smrna->a2 ;
	  tgp = arrayp (tgs, jj, SMRNA) ;
	  if (! tgp->mrnas) tgp->mrnas = arrayHandleCreate (12, SMRNA, ss->h) ;
	  array (tgp->mrnas, arrayMax (tgp->mrnas), SMRNA) = *smrna0 ;
	  smrna0->orfs = 0 ; /* incompatibility with makemrna */
	}
      jj = arrayMax (tgs) ;
      /* name the genes and the mrnas */
      
      for (ii = 0, tgp = arrp (tgs, 0, SMRNA) ; ii < jj ; ii++, tgp++)
	{
	  KEY tg = 0, mrna = 0 ;
	  const char *tgName ;
	  Array tg2m = tgp->mrnas ;
	  int mMax = arrayMax (tg2m), a1 = 0, a2 = 0, iHit ;
	  
	  /* find the min/max mrna->a1 */
	  a1 = a2 = 0 ;
	  for (int m = 0 ; m < mMax ; m++)
	    {
	      SMRNA *mm = arrp (tgp->mrnas, m, SMRNA) ;
	      if (m == 0 || mm->a1 < a1) a1 = mm->a1 ;
	      if (m == 0 || mm->a2 > a2) a2 = mm->a2 ;
	    }
	  
	  /* shift the mrnas coordinates to be relative to the tg coordinates */
	  for (int m = 0 ; m < mMax ; m++)
	    {
	      SMRNA *mm = arrp (tgp->mrnas, m, SMRNA) ;
	      mm->a1 = mm->a1 - a1 + 1 ;
	      mm->a2 = mm->a2 - a1 + 1 ;
	    }
	  /*  cosmid position of the tg */
	  tgp->cosmid = s2m->cosmid ;
	  if (sc->a1 < sc->a2)
	    {
	      tgp->a1 = sc->a1 + a1 - 1 ;
	      tgp->a2 = sc->a1 + a2 - 1 ;
	    }
	  else
	    {
	      tgp->a1 = sc->a1 - a1 + 1 ;
	      tgp->a2 = sc->a1 - a2 + 1 ;
	    }
	  mrnaDesignSetOneGeneEstHits (s2m, sc, ss, tgp) ;
	  
	  /* best cosmid coordinate */
	  compositDesignSelectCosmid (tgp) ;
	  
	  tgp->gene = tg = mrnaDesignGetGeneName (tgp) ;
	  
	  keySet (gmrna->genes, ii) = tg ;
	  tgName = name(tg) ;
	  
	  /* name the mrna */
	  for (int m = 0 ; m < mMax ; m++)
	    {
	      char *mrnaName ;
	      SMRNA *mm = arrp (tgp->mrnas, m, SMRNA) ;
	      mm->gene = tg ;
	      if (m <= 21)
		mrnaName = hprintf (h, "%s.%c", tgName, 'a'+ m) ;
	      else
		mrnaName = hprintf (h, "%s.v%d", tgName, m - 21) ;
	      lexaddkey (mrnaName, &mrna, _VmRNA) ;
	      mm->tr = mrna ;
	    }
	  
	  /* register the mrna->hits */
	  iHit = 0 ;
	  tgp->hits = arrayHandleCreate (64, HIT, ss->h) ;
	  tgp->introns = keySetHandleCreate (ss->h) ;
	  tgp->composite = TRUE ;
	  for (int m = 0, kk = 0 ; m < mMax ; m++)
	    {
	      SMRNA *mm = arrp (tgp->mrnas, m, SMRNA) ;
	      mrnaDesignSetOneMrnaEstHits (s2m, sc, ss, mm) ;
	      mm->composite = TRUE ;
	      
	      for (int i = 0 ; i < arrayMax (mm->hits) ; i++)
		{
		  HIT *up = arrp (mm->hits, i, HIT) ;
		  HIT *vp = arrayp (tgp->hits, iHit++, HIT) ;
		  vp->a1 = mm->a1 + up->a1 - 1 ;
		  vp->a2 = mm->a1 + up->a2 - 1 ;
		  vp->type = up->type ;
		  vp->zone = up->zone ;
		}
	      
	      if (mm->introns)
		{
		  arraySort (tgp->introns, keySetAlphaOrder) ;
		  for (int j = 0 ; j < keySetMax (mm->introns) ; j++)
		    keySet (tgp->introns, kk++) = keySet (mm->introns, j) ;
		}
	    }
	  arraySort (tgp->introns, keySetAlphaOrder) ;
	  keySetCompress (tgp->introns) ;
	  arraySort (tgp->hits, hitA1Order) ;
	  arrayCompress (tgp->hits) ;
	}
    }

  ac_free (h) ;
  
  return ;
} /*  mrnaDesignSetGenes */

/**********************************************************************************/
/**********************************************************************************/
/* cluster the mrnas in genes, A, B, A__B and name these genes and hence the mrnnas */
static void mrnaDesignSplitXGH (S2M *s2m, SC* sc, DS *ss, SMRNA *gmrna, Array smrnas)
{
  /* search transitivelly for contacts
   * start with one gene per mrna
   * check for contact, and merge the genes
   *
   * later we will need to recognize the A__B genes
   */
  /*
  Array tgs, smrnas2 ;
  SMRNA *smrna, *tgp ;
  int ii, jj, iMax = arrayMax (smrnas) ;
  */
  return ;
}  /* mrnaDesignSplitXGH */

/**********************************************************************************/

static void mrnaDesignSaveGenes (S2M *s2m, SC* sc, DS *ss, SMRNA *gmrna, Array smrnas)
{
  SMRNA *tgp ;
  Array tgs = gmrna->tgs ;
  int ii, iMax = arrayMax (tgs) ;

  if (iMax)
    {
      for (ii = 0, tgp = arrp (tgs, 0, SMRNA) ; ii < iMax ; ii++, tgp++)
	{
	  OBJ TG, Cosmid ;
	  Array mrnas = tgp->mrnas ;
	  int jj, jMax = arrayMax (mrnas) ;
	  KEY cosmid = tgp->cosmid ;
	  KEY chrom = tgp->chrom ;
	  int g1 = tgp->g1, g2 = tgp->g2 ;
	  int da, a1 = tgp->a1, a2 = tgp->a2 ;
	  
	  
	  da = (a1 < a2 ? a2 - a1 + 1 : a1 - a2 + 1) ;
	  
	  if ((TG = bsUpdate (tgp->gene)))
	    {
	      bsAddKey (TG, str2tag ("Genomic_sequence"), cosmid) ;
	      bsAddData (TG, _Covers, _Int, &da) ;
	      bsAddData (TG, _bsRight, _Text, "bp from") ;
	      bsAddData (TG, _bsRight, _Int, &a1) ;
	      bsAddData (TG, _bsRight, _Int, &a2) ;
	      if (chrom)
		{
		  bsAddKey (TG, _IntMap, chrom) ;
		  bsAddData (TG, _bsRight, _Int, &g1) ;
		  bsAddData (TG, _bsRight, _Int, &g2) ;
		}
	      mytime_t now = timeNow () ;
	      bsAddData (TG, _Date, _DateType, &now) ;
	      if (tgp->hits)
		{
		  int nI = 0, nX = 0 ;
		  int i, iMax = arrayMax (tgp->hits) ;
		  HIT *up ;
		  for (i = 0, up = arrp (tgp->hits, 0, HIT) ; i < iMax ; i++, up++)
		    {
		      if (up->type & gX) nX++ ;
		      if (up->type & gI) nI++ ;
		    }
		  if (nX) bsAddData (TG, str2tag ("Nb_possible_exons"), _Int, &nX) ;
		  if (nI) bsAddData (TG, str2tag ("Nb_confirmed_introns"), _Int, &nI) ;
		}
	      int k = 1 ;
	      /* Splicing */
	      for (jj = 0 ; jj < arrayMax (tgp->hits) ; jj++)
		{
		  BOOL ok = TRUE ;
		  HIT *up = arrp (tgp->hits, jj, HIT) ;
		  int da = up->a2 - up->a1 + 1 ;
		  bsAddData (TG, _Splicing, _Int, &up->a1) ;
		  bsAddData (TG, _bsRight, _Int, &up->a2) ;
		  if (up->type & gX)
		    {
		      char *cp = messprintf ("%d", k++) ;
		      bsAddKey (TG, _bsRight, _Exon) ;
		      bsAddData (TG, _bsRight, _Text, cp) ;
		    }
		  else
		    {
		      int ig1, ig2 ;
		      char *iName ;
		      KEY intron = 0 ;
		      
		      bsAddKey (TG, _bsRight, _Intron) ;
		      if (g1 < g2)
			{
			  ig1 = g1 + up->a1 - 1 ;
			  ig2 = g1 + up->a2 - 1 ;
			}
		      else
			{
			  ig1 = g1 - up->a1 + 1 ;
			  ig2 = g1 - up->a2 + 1 ;
			}
		      iName = messprintf ("%s__%d_%d", name (chrom), ig1, ig2) ;
		      if (lexword2key (iName, &intron, _VIntron))
			{
			  if (keyFindTag (intron, _gt_ag))
			    bsAddData (TG, _bsRight, _Text, "gt_ag") ;
			  else if (keyFindTag (intron, _gc_ag))
			    bsAddData (TG, _bsRight, _Text, "gc_ag") ;
			  else if (keyFindTag (intron, str2tag("at_ac")))
			    bsAddData (TG, _bsRight, _Text, "at_ac") ;
			  else if (keyFindTag (intron, _Other))
			    {
			      const char *feet = keyGetText (intron, _Other) ;
			      ok = FALSE ;
			      if (0) bsAddData (TG, _bsRight, _Text, "Other") ; /* cannot put 2 texts here */
			      bsAddData (TG, _bsRight, _Text, *feet ? feet : "nn_nn") ;
			    }
			}
		      else
			bsAddData (TG, _bsRight, _Text, "gt_ag") ;
		    }
		  if (ok) bsAddData (TG, _bsRight, _Text, "Length") ;
		  bsAddData (TG, _bsRight, _Int, &da) ;
		  bsAddData (TG, _bsRight, _Text, "bp") ;
		  bsAddData (TG, _bsRight, _Int, &up->zone) ;
		  bsAddData (TG, _bsRight, _Text, "cover") ;
		}
	      for (jj = 0 ; jj < jMax ; jj++)
		{
		  SMRNA *mm = arrp (mrnas, jj, SMRNA) ;
		  
		  bsAddKey (TG, _mRNA, mm->tr) ;
		  bsAddData (TG, _bsRight, _Int, &mm->a1) ;
		  bsAddData (TG, _bsRight, _Int, &mm->a2) ;
		  
		  if (chrom)
		    {
		      OBJ Mrna = bsUpdate (mm->tr) ;
		      if (Mrna)
			{
			  int b1, b2, da ; /* mm coords in cosmid */
			  int gm1, gm2 ;   /* mm coords in chromosome */
			  
			  if (g1 < g2)
			    {
			      b1 = a1 + mm->a1 - 1 ;
			      b2 = a1 + mm->a2 - 1 ;
			      gm1 = g1 + mm->a1 - 1 ;
			      gm2 = g1 + mm->a2 - 1 ;
			    }
			  else
			    {
			      b1 = a1 - mm->a1 + 1 ;
			      b2 = a1 - mm->a2 + 1 ;
			      gm1 = g1 - mm->a1 + 1 ;
			      gm2 = g1 - mm->a2 + 1 ;
			    }
			  da = mm->a2 - mm->a1 + 1 ;
			  bsAddKey (Mrna, str2tag ("Genomic_sequence"), cosmid) ;
			  bsAddData (Mrna, _Covers, _Int, &da) ;
			  bsAddData (Mrna, _bsRight, _Text, "bp from") ;
			  bsAddData (Mrna, _bsRight, _Int, &b1) ;
			  bsAddData (Mrna, _bsRight, _Text, "to") ;
			  bsAddData (Mrna, _bsRight, _Int, &b2) ;
			  
			  mm->chrom = chrom ;
			  mm->g1 = gm1 ; mm->g2 = gm2 ;
			  
			  bsAddKey (Mrna, _IntMap, chrom) ;
			  bsAddData (Mrna, _bsRight, _Int, &gm1) ;
			  bsAddData (Mrna, _bsRight, _Int, &gm2) ;
			  bsSave (Mrna) ;
			}
		    }
		}
	      
	      mrnaDesignSaveOneEstHits (tgp->estHits, TG) ;
	      mrnaDesignSaveOneIntrons (tgp->introns, TG) ;
	      
	      bsSave (TG) ;
	    }
	  if ((Cosmid = bsUpdate (cosmid)))
	    {
	      bsAddKey (Cosmid, _Transcribed_gene, tgp->gene) ;
	      bsAddData (Cosmid, _bsRight, _Int, &a1) ;
	      bsAddData (Cosmid, _bsRight, _Int, &a2) ;
	      
	      for (jj = 0 ; jj < jMax ; jj++)
		{
		  SMRNA *mm = arrp (mrnas, jj, SMRNA) ;
		  int b1, b2 ;
		  
		  if (a1 < a2)
		    {
		      b1 = a1 + mm->a1 - 1 ;
		      b2 = a1 + mm->a2 - 1 ;
		    }
		  else
		    {
		      b1 = a1 - mm->a1 + 1 ;
		      b2 = a1 - mm->a2 + 1 ;
		    }
		  
		  bsAddKey (Cosmid, str2tag("mRNAs"), mm->tr) ;
		  bsAddData (Cosmid, _bsRight, _Int, &b1) ;
		  bsAddData (Cosmid, _bsRight, _Int, &b2) ;
		}
	      bsSave (Cosmid) ;
	    }
	  for (jj = 0 ; jj < jMax ; jj++)
	    {
	      SMRNA *mm = arrp (mrnas, jj, SMRNA) ;
	      mrnaDesignSaveOneMrna (jj == 0 ? tgp : 0, mm) ;
	      mm->composite = TRUE ;
	      if (0) mrnaTiling (sc, mm->tr, FALSE) ;
	    }
	}
    }
  return ;
} /*  mrnaDesignSaveGenes */

/**********************************************************************************/
/**********************************************************************************/

BOOL mrnaDesignUsingCompositeStrategy (S2M *s2m, SC* sc, SMRNA *gmrna, Array smrnas)
{
  AC_HANDLE h = 0 ;
  DS ds ;
  int i ;
  HIT *up ;
  BOOL ok = FALSE ;

  if (arrayMax (smrnas) < 1) 
    return FALSE; 

  /* check if we are using composite reads */
  if (arrayMax (s2m->geneHits))
    for (i = 0, up = arrp (s2m->geneHits, 0, HIT), ok = FALSE ; !ok && i < arrayMax (s2m->geneHits) ; i++, up++)
      if (keyFindTag (up->est, _Composite))
	ok = TRUE ;
  if (! ok)
    return FALSE ;

  memset (&ds, 0, sizeof (ds)) ;
  ds.h = h = ac_new_handle () ;
  ds.exons = arrayHandleCreate (100, DSX, h) ;
  ds.introns = arrayHandleCreate (100, DSX, h) ;
  ds.reads = keySetHandleCreate (h) ;

  ok = FALSE ;
  if (mrnaDesignGetGraphElements (&ds, s2m, sc, gmrna, smrnas))
    {
      ok = TRUE ;
      smrnas = arrayReCreate (smrnas, 12, SMRNA) ;
      mrnaDesignFindPaths (s2m, sc, &ds, smrnas) ;
      mrnaDesignSetOrfs (s2m, sc, &ds, gmrna, smrnas) ;
      mrnaDesignSetMrnas (s2m, sc, &ds, gmrna, smrnas) ;
      mrnaDesignSetGenes (s2m, sc, &ds, gmrna, smrnas) ;
      mrnaDesignSplitXGH (s2m, sc, &ds, gmrna, smrnas) ;

      mrnaDesignSaveGenes (s2m, sc, &ds, gmrna, smrnas) ;
      mrnaAnalyseNeighbours (gmrna->genes) ;
    }
  else
    arrayMax (smrnas) = 0 ;
  ac_free (h) ;
  return ok ;
} /*  mrnaDesignUsingCompositeStrategy */

/**********************************************************************************/
/**********************************************************************************/
#ifdef JUNK
void mrnaDesignAbsorbInIntron (void) 
{
  const char *errors = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  AC_DB db = ac_open_db (0, 0) ; /* local database, cannot fail */
  const char *qqG = "select chr,a1,a2,x,x1,x2,c from x in ?Sequence where x like \"XH_*\" or x like \"XG_*\", chr in x->Intmap, a1 in chr[1], a2 in chr[2], x1 in x->composite, x2 in x1[1], c in x2[1]" ;
  const char *qqI = "select chr,a1,a2,c from x in ?Sequence where x like \"XI_*\" or x like \"XJ_*\", chr in x->Intmap, a1 in chr[1], a2 in chr[2], c in x->composite" ;
  const char *qqA = "select chr,a1,a2,x,c from x in ?Sequence where x like \"XH_*\" or x like \"XG_*\", chr in x->Intmap, a1 in chr[1], a2 in chr[2], c in x->composite" ;

  AC_TABLE aI = ac_bql_table (db, qqI, 0, 0, &errors, h) ;
  AC_TABLE aG = ac_bql_table (db, qqG, 0, 0, &errors, h) ;
  AC_TABLE aA = ac_bql_table (db, qqA, 0, 0, &errors, h) ;

  int iI = 0, iIMax = aI->rows ;
  int iA = 0, iAMax = aA->rows ;
  int iG = 0, iGMax = aG->rows ;

  for (iI = 0 ; iI < iIMax ; iI++) 
    {
      KEY iChr = ac_table_key (aI, iI, 0, 0) ;
      int i1 = ac_table_int (aI, iI, 1, 0) ;
      int i2 = ac_table_int (aI, iI, 2, 0) ;
      int iC = ac_table_int (aI, iI, 3, 0) ;

      continue ;
      
      for (iA ; iA > 0 ; iA--)
	{
	  KEY aChr = ac_table_key (aA, iA, 0, 0) ;
	  int a1 = ac_table_int (aA, iA, 1, 0) ;
	  int a2 = ac_table_int (aA, iA, 2, 0) ;
	  if (aChr < iChr || a2 < i1)
	    break ;
	}
      for (iA ; iA < iAMax ; iA++)
	{
	  KEY aX = ac_table_key (aA, iA, 0, 0) ;
	  KEY aChr = ac_table_key (aA, iA, 1, 0) ;
	  int a1 = ac_table_int (aA, iA, 2, 0) ;
	  int a2 = ac_table_int (aA, iA, 3, 0) ;
	  int aC = ac_table_int (aA, iA, 4, 0) ;
	  if (aChr > iChr || a1 > i2)
	    break ;
	}
    }
  ac_free (h) ;
  return ;
} /* mrnaDesignAbsorbInIntron */
#endif
/**********************************************************************************/
/**********************************************************************************/

BOOL mrnaDesignCutOnePreMrna (KEY gene, KEY green, KEYSET ks)
{
  BOOL debug = FALSE ;
  BOOL new = FALSE ;
  int nn = ks ? keySetMax (ks) : 0 ;
  const char *errors = 0 ;
  AC_HANDLE h = ac_new_handle () ;
  AC_DB db = ac_open_db (0, 0) ; /* local database, cannot fail */
  AC_OBJ Gene = ac_get_obj (db, "Transcribed_gene", name(gene), h) ;
  const char *qqMap = "select chr, a1, a2 from chr in  @->IntMap, a1 in chr[1], a2 in a1[1]" ;
  AC_TABLE map = ac_obj_bql_table (Gene, qqMap, 0, &errors, h) ;
  AC_OBJ Green = ac_get_obj (db, "Sequence", name(green), h) ;
  AC_TABLE greenMap = ac_obj_bql_table (Green, qqMap, 0, &errors, h) ;
  const char *qq2 = "select a1, a2, r, x1, x2 from a1 in @->Assembled_from, a2 in a1[1], r in a2[1], x1 in r[1], x2 in x1[1]" ;
  AC_TABLE afrom = ac_obj_bql_table (Gene, qq2, 0, &errors, h) ;
  AC_TABLE covers = 0 ;
  const char *qqGreen2Gene = "select green, gene, cov, a1, a2, c1, c2 from green in @, gene in  green->From_gene, cov in gene->covers, a1 in gene->assembled_from, a2 in a1[1], r in a1[2] where r == green, c1 in r[1], c2 in r[2]"  ;
  AC_TABLE green2gene = ac_obj_bql_table (Green, qqGreen2Gene, "-4", &errors, h) ;
  int ir, g1 = 0, g2 = 0, gr1 = 0, gr2 = 0, z1 = 0, z2 = 0, dz = 0, xStop = 0, xStart = 0 ;
  const char *rStart = "Xends_ELF" ;
  const char *rStop  = "Xends_ERF" ;
  int nBreaks = 0 ;
  int nGoodBreaks = 0 ;
  KEYSET breaks = keySetHandleCreate (h) ;
  KEYSET bCovers = keySetHandleCreate (h) ;
  Array tgStart = arrayHandleCreate (12, int, h) ;
  Array tgStop = arrayHandleCreate (12, int, h) ;
  char *dna = ac_obj_dna (Green, h) ;
  int d2 = strlen (dna) ;
  vTXT txt = vtxtHandleCreate (h) ;
    
  for (ir = 0 ; map && ir < map->rows ; ir++)
    {
      g1 = ac_table_int (map, ir, 1, 0) ;
      g2 = ac_table_int (map, ir, 2, 0) ;
    }
  if (g1 > g2)
    {
      rStart = "Xends_ERR" ;
      rStop  = "Xends_ELR" ;
    }

  for (ir = 0 ; map && ir < map->rows ; ir++)
    {
      gr1 = ac_table_int (greenMap, ir, 1, 0) ;
      gr2 = ac_table_int (greenMap, ir, 2, 0) ;
    }

  /* get green offset in tg cooordinates */
  if (! afrom || ! map)
    goto done ;
  for (ir = 0 ; ir < afrom->rows ; ir++)
    if (! strcmp (name (green), ac_table_printable (afrom, ir, 2, "toto")))
      { /* green coords in tg */
	z1 = ac_table_int (afrom, ir, 0, 0) ;
	z2 = ac_table_int (afrom, ir, 1, 0) ;
	int z11 = ac_table_int (afrom, ir, 3, 0) ;
	dz = z1 - z11 ;
	break ;
      }
  if (! z1 || ! z2)
    goto done ;

  /* locate potential breakpoints in tg coordinates */
  int cStart = 0, cStop = 0 ;
  for (ir = 0 ; ir < afrom->rows ; ir++)
    {
      if (! strncmp (rStop, ac_table_printable (afrom, ir, 2, "toto"), 9))
	{ /* green coords in tg */
	  AC_OBJ z = ac_table_obj (afrom, ir, 2, h) ;
	  xStop = ac_table_int (afrom, ir, 1, 0) ;
	  cStop = ac_tag_int (z, "Composite", 0) ;
	  ac_free (z) ;
	}
      if (! strncmp ("XSL", ac_table_printable (afrom, ir, 2, "toto"), 3))
	{ /* green coords in tg */
	  AC_OBJ z = ac_table_obj (afrom, ir, 2, h) ;
	  xStart = ac_table_int (afrom, ir, 0, 0) ;
	  cStart = ac_tag_int (z, "Composite", 0) ;
	  ac_free (z) ;
	  if (xStop && xStop < xStart + 30)
	    {
	      int bb = (xStop + xStart) / 2 ;
	      if (bb > z1 + 10 && bb < z2 - 10)
		{
		  keySet (breaks, nBreaks) = bb - dz ;
		  keySet (bCovers, nBreaks) = -100 ;
		  nBreaks++ ;
		}
	      xStop = 0 ;
	    }
	}
      if (! strncmp (rStart, ac_table_printable (afrom, ir, 2, "toto"), 9))
	{ /* green coords in tg */
	  AC_OBJ z = ac_table_obj (afrom, ir, 2, h) ;
	  xStart = ac_table_int (afrom, ir, 0, 0) ;
	  cStart = ac_tag_int (z, "Composite", 0) ;
	  ac_free (z) ;
	  if (xStop && xStop < xStart + 30)
	    {
	      int bb = (xStop + xStart) / 2 ;
	      if (bb > z1 + 10 && bb < z2 - 10)
		{
		  keySet (breaks, nBreaks) = bb - dz ;
		  keySet (bCovers, nBreaks) = (cStop + cStart) / 2 ;
		  nBreaks++ ;
		}
	      xStop = 0 ;
	    }
	}
    }

  
  const char *qq1 = "select x1, x2, c from x1 in @->Composite, x2 in x1[1], c in x2[1]" ;
  covers = ac_obj_bql_table (Green, qq1, 0, &errors, h) ;
  if (covers && nBreaks)
    { /* check if the candidate breakpoint falls in an exon cover dip */
      int ir = 0 ;

      for (int ib = 0 ; ib < nBreaks ; ib++)
	{
	  int u1, u2, c, bb = keySet (breaks, ib) ;
	  for (; ir >= 0 ; ir--)
	    {
	      u1 = ac_table_int (covers, ir, 0, 0) ;
	      if (u1 <= bb)
		break ;
	    }
	  if (ir < 0) ir = 0 ;
	  for (; ir < covers->rows ; ir++)
	    {
	      u2 = ac_table_int (covers, ir, 1, 0) ;
	      c = ac_table_int (covers, ir, 2, 0) ;
	      if (u2 >= bb)
		break ;
	    }
	  if (ir >= covers->rows)
	    { ir-- ; continue ; }
	  u1 = ac_table_int (covers, ir, 0, 0) ;
	  c = ac_table_int (covers, ir, 2, 0) ;	      
	  if (u1 <= bb && u2 >= bb)
	    {
	      int ok = 0 ;
	      int bc = keySet (bCovers, ib) ;
	      if (bc == -100) ok = 3 ;
	      if (! ok || 4 * c < bc)
		{
		  for (int j = 1 ; j < 4  && j + ir < covers->rows ; j++)
		    if (1 * c < ac_table_int (covers, ir + j, 2, 0))
		      { ok = 1 ; break ; }
		  for (int j = 1 ; j < 4  && ir - j >= 0 ; j++)
		    if (1 * c < ac_table_int (covers, ir - j, 2, 0))
		      { ok += 2 ; break ; }
		}
	      if (ok != 3)
		keySet (breaks, ib) = 0 ;  /* drop it */
	      else
		nGoodBreaks++ ;
	    }
	}
    }

  /* get green offset in other genes */
  for (int itg = 0 ; itg < green2gene->rows ; itg++)
      { /* green coords in tg */
	KEY tg = ac_table_key (green2gene, itg, 1, 0) ;
	int cov = ac_table_int (green2gene, itg, 2, 0) ;
	int u1 = ac_table_int (green2gene, itg, 3, 0) ;
	int u2 = ac_table_int (green2gene, itg, 4, 0) ;
	int v1 = ac_table_int (green2gene, itg, 5, 0) ;
	int v2 = ac_table_int (green2gene, itg, 6, 0) ;

	fprintf (stderr , "################ %s %d %d %d %d %d\n", name(tg), cov, u1, u2, v1, v2) ;
	array (tgStart, itg, int) = v1 - u1 + 1 ;
	array (tgStop, itg, int) = v1 - u1 + 2 + cov - 1 ;
	if (itg > 0 && array (tgStop, itg - 1, int) < array (tgStart, itg, int) + cov/3)
	  {
	    int bb = (array (tgStop, itg - 1, int) + array (tgStart, itg, int))/2 ;
	    keySet (breaks, nBreaks++) = bb ;
	    nGoodBreaks++ ;
	  }
      }

  /* remove green from overlapping genes */
  for (int itg = 0 ; itg < green2gene->rows - 1 ; itg++)
    { /* green coords in tg */
      KEY tg = ac_table_key (green2gene, itg, 1, 0) ;
      if (  /* itg includes itg+1 */
	    d2 > array (tgStart, itg, int) &&
	    1 < array (tgStop, itg, int) &&
	    d2 > array (tgStart, itg+1, int) &&
	    1 < array (tgStop, itg+1, int) &&
	    array (tgStop, itg, int) >  array (tgStop, itg + 1, int)
	    )
	vtxtPrintf (txt, "Transcribed_gene %s\n-D Read %s\n-D cdna_clone %s\n\n"
		    , name (tg), name (green)  , name (green)
		    ) ;
    }
  
 done:
  nBreaks = 0 ;
  if (nGoodBreaks)
    { /* compress */
      int jb = 0, bb0 = 10 ;
      
      keySetSort (breaks) ;
      for (int ib = 0 ; ib < keySetMax(breaks) ; ib++)
	{
	  int bb = keySet (breaks, ib) ;
	  if (bb && bb > bb0 && bb < d2 - 10)
	    { keySet (breaks, jb++) = bb ; bb0 = bb + 10 ; }
	}
      keySet (breaks, jb++) = d2 ;
      nGoodBreaks = nBreaks = keySetMax (breaks) = jb ;
    }


  if (nBreaks)
    {
      keySet (breaks, nBreaks++) = d2 ;

      vtxtPrintf (txt, "Sequence %s\n-D Is_read\n-D From_gene\n-D In_mRNA\n\n", name (green)) ;
      vtxtPrintf (txt, "cDNA_clone %s\n-D From_gene\n-D In_mRNA\n\n", name (green)) ;
      vtxtPrintf (txt, "-D DNA %s\n\n", name (green)) ;
      vtxtPrintf (txt, "-D cDNA_clone %s\n\n", name (green)) ;
      vtxtPrintf (txt, "-D Sequence %s\n\n", name (green)) ;

      /* remove green from all genes */
      for (int itg = 0 ; itg < green2gene->rows ; itg++)
	{ /* green coords in tg */
	  KEY tg = ac_table_key (green2gene, itg, 1, 0) ;
	  int u1 = ac_table_int (green2gene, itg, 3, 0) ;
	  int u2 = ac_table_int (green2gene, itg, 4, 0) ;
	  
	  vtxtPrintf (txt, "Transcribed_gene %s\n-D Assembled_from %d %d %s\n\n", name(tg), u1, u2, name(green)) ;
	}

      if (covers)
	{
	  int d11 = 1, d12 ;
	  for (int ib = 0 ; ib < nBreaks ; ib++)
	    {
	      char *nam = hprintf (h, "%s.%d", name (green), ib+1) ;
	      int bb = keySet (breaks, ib) ;
	      
	      if (bb < d11)
		continue ;
	      if (bb > d2)
		continue ;
	      fprintf (stderr, "Breaking %s at position %d\n", name (green), bb) ;
	      
	      d12 = bb ;
	      if (d12 < 1)
		d12 = 1 ;
	      if (d12 < d2)
		d12 -= 10 ;
	      /* char *dnaZone = ac_zone_dna (Green, d11, d12, h) ; */
	      dna[d12]=0 ;
	      vtxtPrintf (txt, "DNA %s\n%s\n\n", nam, dna + d11 - 1) ;
	      
	      for (int itg = 0 ; itg < green2gene->rows ; itg++)
		{ /* green coords in tg */
		  KEY tg = ac_table_key (green2gene, itg, 1, 0) ;
		  if (d12 > array (tgStart, itg, int) &&
		      d11 < array (tgStop, itg, int)
		      )
		    {
		      if (itg < green2gene->rows - 1 &&
			  d12 > array (tgStart, itg+1, int) &&
			  d11 < array (tgStop, itg+1, int) &&
			  array (tgStop, itg, int) >  array (tgStop, itg + 1, int) /* itg includes itg+1 */
			  )
			;
		      else
			{
			  vtxtPrintf (txt, "cDNA_clone %s\nRead %s\nFrom_gene %s\n\n", nam, nam, name (tg)) ;
			  int u1 = d11 - array (tgStart, itg, int) + 1 ;
			  int u2 = d12 - array (tgStart, itg, int) + 1 ;
			  vtxtPrintf (txt, "Transcribed_gene %s\nRead %s\nAssembled_from %d %d %s %d %d\n\n", name(tg), nam, u1, u2, nam, d11 - d11 + 1, d12 - d11 + 1) ;
			}
		    }
		}
	      vtxtPrintf (txt, "Sequence %s\nIs_read\nForward\nColour GREEN3\n", nam) ;
	      if (map->rows)
		vtxtPrintf (txt, "IntMap %s %d %d\n"
			    , ac_table_printable (map, 0, 0, 0)
			    , gr1 < gr2 ? gr1 + d11 - 1 : gr1 - d11 + 1
			    , gr1 < gr2 ? gr1 + d12 - 1 : gr1 - d12 + 1
			    ) ;
	      
	      for (int ir = 0; ir < covers->rows ; ir++)
		{
		  int u1 = ac_table_int (covers, ir, 0, 0) ;
		  int u2 = ac_table_int (covers, ir, 1, 0) ;
		  int c = ac_table_int (covers, ir, 2, 0) ;
		  if (u1 > d12)
		    break ;
		  if (u1 < d11) u1 = d11 ;
		  if (u2 > d12) u2 = d12 ;
		  if (u1 <= u2)
		    vtxtPrintf (txt, "Composite %d %d %d\n", u1 - d11 + 1, u2 - d11 + 1, c) ;
		}
	      vtxtPrintf (txt, "\n") ;
	      d11 = bb + 10 ;
	    }
	}
    }
  if (vtxtPtr (txt))
    {
      ac_parse (db, vtxtPtr (txt), &errors, 0, h) ;
      if (debug) fprintf (stderr, "%s\n\n", vtxtPtr (txt)) ;
    }
  if (ks)
    keySet (ks, nn++) = green ;
  ac_free (h) ;
  ac_free (db) ;
  
  return new ;
} /* mrnaDesignCutOnePreMrna */

/**********************************************************************************/
/* cut the premra XG/XH exons reads on strong termination signals */

static void mrnaDesignCutGenePreMrna (KEY gene)
{
  KEYSET greens = queryKey (gene, ">Read ; Is_Read && (IS XG_* || IS XH_*)") ;
  for (int i = 0 ; i < keySetMax (greens) ; i++)
    {
      KEY green = keySet (greens, i) ;
      mrnaDesignCutOnePreMrna (gene, green, 0) ;
    }
  return ;
} /* mrnaDesignCutGenePreMrna */

/**********************************************************************************/
/* cut the premra XG/XH exons reads on strong termination signals */
void mrnaDesignCutAllPreMrna (KEYSET tgs)
{
  KEYSET genes = query (tgs, "CLASS Transcribed_genes") ;

  if (genes)
    for (int i = 0 ; i < keySetMax (genes) ; i++)
      mrnaDesignCutGenePreMrna (keySet (genes, i)) ;
  keySetDestroy (genes) ;
  
  return ;
} /* mrnaDesignCutAllPreMrna */
 
/**********************************************************************************/
/**********************************************************************************/

int mrnaDesignCleanEcho (void)
{
  AC_HANDLE h = ac_new_handle () ;
  const char *errors = 0 ;
  AC_DB db = ac_open_db (0, &errors) ; /* local database, cannot fail */
  KEY _ps = str2tag ("ProcessStatus") ;  
  const char *qqxg = "select r, chr, a1, a2, s, x1, x2, c from r in ?Sequence where r like \"XG_*\" && !r#ProcessStatus , chr in r->IntMap, a1 in chr[1], a2 in chr[2] where -a2 < 330000, s in r->strand, x1 in r->Composite, x2 in x1[1], c in x2[1]" ;
  AC_TABLE xxg = ac_bql_table (db, qqxg, 0, 0, &errors, h) ;
  const char *qqxxi = "select r, chr, a1, a2, s, c from r in ?Sequence where r#composite && r like \"XI_*\", I in r->intron, chr in I->IntMap, a1 in chr[1], a2 in chr[2] where -a2 < 330000, s in r->strand, c in r->Composite" ;
  AC_TABLE xxi = ac_bql_table (db, qqxxi, 0, 0, &errors, h) ;
  const char *qqxxs = "select r, chr, a1, a2, s, c from r in ?Sequence where r#composite && !r#ProcessStatus && ! r like \"XG_*\" && ! r like \"XH_*\", chr in r->IntMap, a1 in chr[1], a2 in chr[2] where -a2 < 330000, s in r->strand, c in r->Composite" ;
  AC_TABLE xxs = ac_bql_table (db, qqxxs, 0, 0, &errors, h) ;
  KEYSET maps = keySetHandleCreate (h) ;
  vTXT txt = vtxtHandleCreate (h) ;
  int ccc = 0, nKilled = 0 ;
  KEY oldxg = 0 ;
      
  if (xxg) for (int ir = 0 ; ir < xxg->rows ; ir++)
    {
      KEY map = ac_table_key (xxg, ir, 1, 0) ;
      if (map)
	keySetInsert (maps, map) ;
    }
  for (int ii = 0 ; ii < keySetMax (maps) ; ii++)
    {
      AC_HANDLE h1 = ac_new_handle () ;
      KEY map = keySet (maps, ii) ;
      Array pm, plus = arrayHandleCreate (100000, int, h1) ;
      Array minus = arrayHandleCreate (100000, int, h1) ;

      /* construct a wiggle for each strand */
      for (int ir = 0 ; xxi && ir < xxi->rows ; ir++)
	{
	  KEY map1 = ac_table_key (xxi, ir, 1, 0) ;
	  if (map == map1)
	    {
	      int a1 = (ac_table_int (xxi, ir, 2, 0) + 4) / 10 ;
	      int a2 = (ac_table_int (xxi, ir, 3, 0) + 4) / 10 ;
	      KEY strand = ac_table_key (xxi, ir, 4, 0) ;
	      int c = ac_table_int (xxi, ir, 5, 0) ;

	      if (strand == _Reverse) { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }
	      pm = a1 < a2 ? plus : minus ;
	      if (a1 > a2) { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }
	      for (int a = a1 + 1 ; a < a2 ; a++)
		{
		  int c1 = array (pm, a, int) ;
		  array (plus, a, int)  = c + c1 ; /* introns are cumulative */
		  array (minus, a, int)  = c + c1 ; /* introns are cumulative */
		}
	    }
	}

      for (int ir = 0 ; xxg && ir < xxg->rows ; ir++)
	{
	  KEY map1 = ac_table_key (xxg, ir, 1, 0) ;
	  if (map == map1)
	    {
	      int b1, b2, a1 = ac_table_int (xxg, ir, 2, 0) ;
	      int a2 = ac_table_int (xxg, ir, 3, 0) ;
	      int x1 = ac_table_int (xxg, ir, 5, 0) ;
	      int x2 = ac_table_int (xxg, ir, 6, 0) ;
	      int c = ac_table_int (xxg, ir, 7, 0) ;

	      if (x1 > x2) continue ;
	      if (a1 <= a2) { b1 = a1 + x1 - 1 ; b2 = a1 + x2 - 1 ; }
	      else { b1 = a1 - x1 + 1 ; b2 = a1 - x2 + 1 ; }
	      if (b1 > b2) { int b0 = b1 ; b1 = b2 ; b2 = b0 ; }
	      b1 = (b1 + 4) / 10 ; b2 = (b2 + 4) / 10 ;
	      pm = a1 < a2 ? plus : minus ;
	      for (int a = b1 + 1 ; a < b2 ; a++)
		{
		  int c1 = array (pm, a, int) ;
		  if (c1 < c) array (pm, a, int)  = c ;
		}
	    }
	}

      /* now clean up */
      for (int ir = 0 ; xxs && ir < xxs->rows ; ir++)
	{
	  KEY map1 = ac_table_key (xxs, ir, 1, 0) ;
	  if (map == map1)
	    {
	      int a1 = (ac_table_int (xxs, ir, 2, 0) + 4) / 10 ;
	      int a2 = (ac_table_int (xxs, ir, 3, 0) + 4) / 10 ;
	      KEY strand = ac_table_key (xxs, ir, 4, 0) ;
	      int n = 0, c = ac_table_int (xxs, ir, 5, 0) ;
	      long int ccL = 0 ;
	      
	      if (strand == _Reverse) { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }
	      pm = a1 > a2 ? plus : minus ; /* opposite strand wiggle */
	      if (a1 > a2) { int a0 = a1 ; a1 = a2 ; a2 = a0 ; }
	      for (int a = a1 + 1 ; a < a2 ; a++)
		{
		  int c1 = array (pm, a, int) ;
		  ccL += c1 ; n++ ;
		}
	      if (n) ccL /= n ; /* average premessenger of other strand */
	      ccL /= 100 ;
	      if (ccL > c)
		ccc = 0 ;
	      else
		ccc = c - ccL ;
	      if (ccc < 0) ccc = 0 ;
	      vtxtPrintf (txt, "Sequence %s\nProcessStatus  %d\n-D Composite\n", ac_table_printable (xxs, ir, 0, "toto"), c) ;
	      if (ccc > 0)
		vtxtPrintf (txt, "Composite %d\n\n", ccc) ;
	      else
		{
		  nKilled++ ;
		  vtxtPrintf (txt, "-D Is_read\nColour BROWN\n\n", c) ;
		}
	    }
	}

      ccc = 0 ;
      oldxg = 0 ;
      for (int ir = 0 ; xxg && ir < xxg->rows ; ir++)
	{
	  KEY xg = ac_table_key (xxg, ir, 0, 0) ;
	  KEY map1 = ac_table_key (xxg, ir, 1, 0) ;
	  int cc = 0 ;
	  
	  if (xg && map == map1 && ! keyFindTag (xg, _ps))
 	    {
	      int b1, b2, n = 0, a1 = ac_table_int (xxg, ir, 2, 0) ;
	      int a2 = ac_table_int (xxg, ir, 3, 0) ;
	      int x1 = ac_table_int (xxg, ir, 5, 0) ;
	      int x2 = ac_table_int (xxg, ir, 6, 0) ;
	      int c = ac_table_int (xxg, ir, 7, 0) ;

	      if (!x1 || x1 > x2) continue ;
	      if (a1 < a2) { b1 = a1 + x1 - 1 ; b2 = a1 + x2 - 1 ; }
	      else { b1 = a1 - x1 + 1 ; b2 = a1 - x2 + 1 ; }
	      if (b1 > b2) { int b0 = b1 ; b1 = b2 ; b2 = b0 ; }
	      b1 = (b1 + 4) / 10 ; b2 = (b2 + 4) / 10 ;
	      pm = a1 > a2 ? plus : minus ;
	      for (int a = b1 + 1 ; a < b2 ; a++)
		{
		  if (a < 0) invokeDebugger () ;
		  int c1 = array (pm, a, int) ;
		  cc += c1 ; n++ ;
		}
	      if (n) cc /= n ; /* average premessenger of other strand */
	      cc = c - cc / 100 ; /* new value */
	      if (cc < 0) cc = 0 ;
	      if (xg != oldxg)
		{
		  if (oldxg && ! ccc)
		    {
		      nKilled++ ;
		      vtxtPrintf (txt, "-D Is_read\nColour BROWN\n") ;
		    }
		  oldxg = xg ;
		  ccc = 0 ;
		  vtxtPrintf (txt, "\nSequence %s\nProcessStatus\n-D Composite\n", name (xg)) ;
		}
	      if (cc > 0 && x1 <= x2)
		{
		  ccc += cc ;
		  vtxtPrintf (txt, "Composite %d %d %d\n", x1, x2, cc) ;
		}
	    }
	}
      if (oldxg && ! ccc)
	{
	  nKilled++ ;
	  vtxtPrintf (txt, "-D is_read\nColour BROWN\n") ;
	}
      vtxtPrintf (txt, "\n") ;
    }

  if (vtxtPtr (txt))
    ac_parse (db, vtxtPtr (txt), &errors, 0, h) ;

  fprintf (stderr, "mrnaDesignCleanEcho kill %d reads in %d maps\n", nKilled, keySetMax (maps)) ;
  ac_free (h) ;

  return nKilled ;
} /* mrnaDesignCleanEcho */

/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
