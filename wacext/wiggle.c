/*
 * Create and combine wiggles, i.e. tag density profiles
 * Several formats are recognized as input or output
 *   probealignbest hit plain and binary output
 *   three formats from UCSC
 *   the AW binary format genberated by the present program
 * Create letter profiles from a fasta or fastc file
 */


#include "wiggle.h"
static void usage (const char *error) ;

/*************************************************************************************/
/*************************************************************************************/

static int sxMapOrder (const void *va, const void *vb)
{
  const SXMAP *a = (const SXMAP *) va, *b = (const SXMAP *) vb ;
  int x ;

  x = a->map - b->map ; if (x) return x ;
  x = a->x1 - b->x1 ; if (x) return x ;
  x = a->x2 - b->x2 ; if (x) return x ;

  x = a->remap - b->remap ; if (x) return x ;
  x = a->a1 - b->a1 ; if (x) return x ;
  x = a->a2 - b->a2 ; if (x) return x ;
 
  return 0 ;
}

/*************************************************************************************/
/* expect 6 columns, 
 * target x1 x2 chrom a1 a2
 * target could be a cosmid or an NM
 *   target a1 a2 is as found in the hits file
 *   is is remapped to chrom a1 a2
 */ 
static int sxGetMap (WIGGLE *sx) 
{
  int nn = 0, target_map, map, x1, x2, remap, a1, a2 ;
  char *cp, buf[10000] ;
  AC_HANDLE h = 0 ;
  ACEIN ai = 0 ;
  SXMAP *sxmap, *sxmap2 ;

  sx->mapDict = dictHandleCreate (1024, sx->h) ;
  sx->map2remap = keySetHandleCreate (h) ;
  sx->target_mapDict = dictHandleCreate (32, sx->h) ;
  sx->remapDict = dictHandleCreate (1024, sx->h) ;
  if (! sx->mapFileName) 
    return 0 ;
  dictAdd (sx->mapDict, "______toto", 0) ;
  dictAdd (sx->remapDict, "______toto", 0) ;
  dictAdd (sx->target_mapDict, "______toto", 0) ;
  h = ac_new_handle () ;

  ai = aceInCreate ( sx->mapFileName, sx->gzi, h) ;
  if (! ai)
    usage ("Cannot open the map file") ;

  aceInSpecial (ai, "\t\n") ;
  sx->targets = arrayHandleCreate (10000, SXMAP, sx->h) ;

  while (ai && aceInCard (ai))
    {
      map = remap = x1 = x2 = a1 = a2 = 0 ;
      cp = aceInWord (ai) ; /* target: av pg HINV ... */
      if (!cp) 
	continue ;
      dictAdd (sx->target_mapDict, cp, &target_map) ;

      cp = aceInWord (ai) ; /* mrna name */
      if (!cp) 
	continue ;
      if (!strncmp (cp, "MRNA:", 5)) cp += 5 ;

      dictAdd (sx->mapDict, messprintf ("%s:%s",dictName (sx->target_mapDict, target_map),cp), &map) ;
      if (aceInInt (ai, &x1))
	{
	  if (! aceInInt (ai, &x2))
	    messcrash ("missing coordinates in column 3 in %s line %d", 
		       sx->mapFileName, aceInStreamLine (ai)) ;

	  cp = aceInWord (ai) ;
	  if (cp) 
	    {
	      strncpy (buf, cp, 9999) ;
	      if (! aceInInt (ai, &a1) || ! aceInInt (ai, &a2))
		messcrash ("missing coordinates in column 5 or 6 in %s line %d", 
			   sx->mapFileName, aceInStreamLine (ai)) ;
	      /*
		if ((sx->strand && a1 > a2) || (sx->antistrand && a1 < a2)) 
		continue ;
	      */ 
	      if (buf[0] && sx->sxxChromosome && strcasecmp (sx->sxxChromosome, buf)) continue ;
	      dictAdd (sx->remapDict, buf, &remap) ;
	    }
	  else
	    { /* duplicate column 1,2,3 */
	      const char * ccp = dictName (sx->mapDict, map) ;
	      if (ccp && sx->sxxChromosome && strcasecmp (sx->sxxChromosome, ccp)) continue ;
	      dictAdd (sx->remapDict, ccp, &remap) ;
	      a1 = x1 ; a2 = x2 ;	      
	    }
	}
      else
	{
	  const char *ccp = dictName (sx->mapDict, map) ;
	  if (ccp && sx->sxxChromosome && strcasecmp (sx->sxxChromosome, ccp)) continue ;
	  dictAdd (sx->remapDict, ccp, &remap) ;
	  a1 = x1 = 1 ; a2 = x2 = 1<<30 ;
	}
	  
      if (x1 > x2)
	{ 
	  int x ;

	  x = x1 ; x1 = x2 ; x2 = x ; 
	  x = a1 ; a1 = a2 ; a2 = x ;
	}

      if (
	  (a1 < a2 && a2 - a1 != x2 - x1) || 
	  (a1 > a2 && a1 - a2 != x2 - x1)
	  )
	{
	  if (0) messerror ("In %s x2-x1=(%d - %d) = %d  !=  a2-a1 = (%d-%d) = %d", sx->mapFileName, x2, x1, x2 - x1, a2, a1, a2 - a1) ;
	  continue ;
	}
      if (remap) keySet (sx->map2remap, map) = remap ;
      sxmap = arrayp (sx->targets, nn++, SXMAP) ;
      sxmap->map = map ;
      sxmap->remap = remap ;
      sxmap->target_map = target_map ;

      sxmap->x1 = x1 ;
      sxmap->x2 = x2 ;
      sxmap->a1 = a1 ;
      sxmap->a2 = a2 ;
    }
  arraySort (sx->targets, sxMapOrder) ;

  /* discard overlapping segments */
  if (nn > 1)
    {
      int i ;
      for (i = 1, sxmap = arrp (sx->targets, i, SXMAP), sxmap2 = sxmap - 1 ; i < nn ; sxmap++, sxmap2++, i++)
	{
	  if (sxmap2->map == sxmap->map && sxmap2->x2 >= sxmap->x1)
	    sxmap2->x2 = sxmap->x1 - 1 ;
	}
    }
  ac_free (h) ;
  fprintf (stderr, "// wiggle -remap : sxTargetWiggleGetMap done %s\n", timeShowNow()) ;
  return nn ;
} /* sxTargetWiggleGetMap */

/*************************************************************************************/
/*************************************************************************************/

static void sxGetSelectionDo (ACEIN ai, DICT *dict)
{
  const char *ccp ;

  while (aceInCard (ai))
    {
      ccp = aceInWord (ai) ;
      if (ccp) dictAdd (dict, ccp, 0) ;
    }
} /* sxGetSelectionDo */

/***********/

static void sxGetSelection (WIGGLE *sx)
{
  ACEIN ai ;

  if (sx->selectFileName)
    {
      if (sx->in != BHIT)
	messcrash ("Sorry, option -select is only compatible with -I BHIT, since otherwise the identifiers are not available") ;
      sx->selectDict = dictHandleCreate (1000000, sx->h) ;
      ai = aceInCreateFromFile (sx->selectFileName, "r", 0, sx->h) ;
      sxGetSelectionDo (ai, sx->selectDict) ;
      ac_free (ai) ;
    }
  if (sx->rejectFileName)
    {
      if (sx->in != BHIT)
	messcrash ("Sorry, option -select is only compatible with -I BHIT, since otherwise the identifiers are not available") ;
      sx->rejectDict = dictHandleCreate (1000000, sx->h) ;
      ai = aceInCreateFromFile (sx->rejectFileName, "r", 0, sx->h) ;
      sxGetSelectionDo (ai, sx->rejectDict) ;
      ac_free (ai) ;
    }
} /* sxGetSelection */

/***********/

static void sxVentilate (WIGGLE *sx)
{
  AC_HANDLE h = 0 ;
  ACEIN ai = sx->ai ;
  ACEOUT ao = 0, eo = 0, po = 0 ;
  Array aos = arrayHandleCreate (64,ACEOUT,h) ;
  Array eos = arrayHandleCreate (64,ACEOUT,h) ;
  Array pos = arrayHandleCreate (64,ACEOUT,h) ;
  char *ccp, probeBuf[1024], gBuf[1024],  mBuf[1024], target_class = 0 , buf[10001], wall1[256], wall2[256], *fr, *Efr ;
  char oldProbeBuf[1024] ;
  int chrom, score, tc, tc2, a1, a2, e1, e2, p1, p2, x1, x2, ali, toali, map, mult, tmult, firstX1 = 0 ;
  int oldScore = 0,  oldX1 = 0, oldX2 = 0, oldA1 = 0, oldA2 = 0 ;
  int oldMap = 0, oldtc = 0, tcGenome = 0 ;
  /* int  mateAli = 0, mateToAli = 0,  oldAli = 0, oldToAli = 0  ; */
  int nN = 0, nErr = 0 ;
  int ifr, iEfr ;
  float weight ;
  int newProbe ;
  BOOL unique = sx->unique ;
  BOOL non_unique = sx->non_unique ;
  BOOL partial = sx->partial ;
  BOOL noPartial = sx->noPartial ;
  BOOL hierarchic = sx->hierarchic ;
  BOOL firstRead ;
  BOOL inTranscript = FALSE ;
  BOOL isPartial = FALSE ;
  int lastScore = 0, lastAli = 0, lastToAli = 0 ;
  int EndLength = 30 ;
  int chainNumber, chainFrom, chainTo, mateScore, mateAli , mateToAli ;

  probeBuf[0] = 0 ;   oldProbeBuf[0] = 0  ;  
  gBuf[0] = 0 ; mBuf[0] = 0 ;
  
  dictAdd (sx->target_mapDict, "Z_genome", &tcGenome) ; 

  aceInSpecial (ai, "\n") ;
  while (aceInCard (ai))
    {
      ccp = aceInPos (ai) ;
      if (! ccp || *ccp == '#') continue ;
      strncpy (buf, ccp, 10000) ;
      
      ccp = aceInWord (ai) ;
      if (! ccp || *ccp == '#') continue ;

      newProbe = 0 ; 
      firstRead = TRUE ;
      if (sx->pair)
	{
	  int n = strlen (ccp) ;
	  if (ccp[n-1] == '<')
	    firstRead = FALSE ;
	}
      if ( strcmp (ccp, probeBuf))
	{
	  newProbe = 1 ;
	  strcpy (oldProbeBuf, probeBuf) ;
	  oldScore = lastScore = 0 ;
	  /* oldAli = lastAli ; 	  oldToAli = lastToAli ; */
	  lastScore = lastAli = lastToAli = 0 ;
	  strcpy (probeBuf, ccp) ;
	  
	  oldtc = oldMap = gBuf[0] = 0 ;
 
#ifdef JUNK
	  mateAli = mateToAli = 0 ; 
	  if (sx->pair)
	    {
	      int n = strlen (ccp) ;
	      if (! strncmp (ccp, oldProbeBuf, n-1) && 
		  (ccp[n-1] == '<' || ccp[n-1] == '>')  &&
		  (oldProbeBuf[n-1] == '<' || oldProbeBuf[n-1] == '>')
		  )
		{
		  mateAli = oldAli ;
		  mateToAli = oldToAli ;
		}
	    }
#endif	      
	}
	  
 
      strncpy (probeBuf, ccp, 999) ;
      if (! aceInStep (ai, '\t') || ! aceInInt (ai, &score)) continue ;
      if (! aceInStep (ai, '\t') || ! aceInInt (ai, &mult)) continue ;
      if (! aceInStep (ai, '\t') || ! aceInInt (ai, &toali)) continue ;
      if (! aceInStep (ai, '\t') || ! aceInInt (ai, &ali)) continue ;

      if (! aceInStep (ai, '\t') || ! aceInInt (ai, &x1)) continue ;
      if (! aceInStep (ai, '\t') || ! aceInInt (ai, &x2)) continue ;
      if (! aceInStep (ai, '\t') || ! (ccp = aceInWord(ai)) || *ccp > 'Z') continue ;

      if (sx->target_class && ccp && strcmp (ccp, sx->target_class))
	continue ;
      target_class = *ccp ;
      if (hierarchic && target_class == 'Z')
	continue ;
      dictAdd (sx->target_mapDict, ccp, &tc) ; 
      dictAdd (sx->remapDict, ccp, &chrom) ; 

      /* gene */
      if (! aceInStep (ai, '\t') || ! (ccp = aceInWord(ai))) continue ;
      /* newGene =  strcmp (ccp, gBuf) ; */
      strncpy (gBuf, ccp, 1000) ;

      if (! aceInStep (ai, '\t') || ! aceInInt (ai, &tmult)) continue ;
      if (! aceInStep (ai, '\t') || ! (ccp = aceInWord(ai))) continue ;

      /* newMrna =   newGene ? TRUE : strcmp (ccp, mBuf) ; */
      strncpy (mBuf, ccp, 1000) ;

      /* ok, ccp == the chromosome or the mRNA */
      map = tc2 = 0 ;	  

      if (0 && target_class != 'B') /* hack to do just the rrna */
	continue ;
      inTranscript = FALSE ;
      switch (target_class)
	{
	case 'A': /* mito */ 
          dictAdd (sx->mapDict, "mito", &map) ;
	  break ;
	case '0': /* SpikeIn */ 
          dictAdd (sx->mapDict, ccp, &map) ;
	  break ;
	case 'z': /* genome */
	case 'Z': /* genome */
	  if (strchr(ccp, '|') || strlen (ccp)>100) continue ;
	  dictAdd (sx->remapDict, ccp, &chrom) ;
	  dictAdd (sx->mapDict, ccp, &map) ;
	  break ;
	case 'B': /* ribosomal rrna now integrated in the mrnaRemap.gz table is treated as RefSeq and AceView */ 
	default:
	  if (keySetMax (sx->map2remap))
	    {
	      if (! dictFind (sx->mapDict, messprintf("%s:%s", dictName (sx->target_mapDict, tc), ccp), &map))
		continue ;
	      inTranscript = TRUE ;
	      chrom = keySet (sx->map2remap, map) ;
	      if (! chrom)
		continue ;
	      tc2 = 1 ;
	    }
	  else
	    dictAdd (sx->mapDict, ccp, &map) ;
	  break ;
	} 

      /* check unicity */
      if (tmult != 1 && tmult != -2 && unique) 
	continue ;
      if ((tmult == 1 || tmult == -2) && non_unique) 
	continue ;

      if (! aceInStep (ai, '\t') || ! aceInInt (ai, &a1)) continue ;
      if (! aceInStep (ai, '\t') || ! aceInInt (ai, &a2)) continue ;
      if (! aceInStep (ai, '\t') || ! aceInInt (ai, &nN)) continue ;
      if (! aceInStep (ai, '\t') || ! aceInInt (ai, &nErr)) continue ;

      /************* restore natural order ******************/
      /* x1<x2 is the convention, this should never be necessary, do it just in case */
      if (x1 > x2)  
	{ int x0 = x1 ; x1 = x2 ; x2 = x0 ; x0 = a1 ; a1 = a2 ; a2 = x0 ; }

      /************* remove exon overlap *********************/
      /* hierarchic clean up before we map to the wiggles */
      if (oldtc && target_class != oldtc)
	continue ;

      /* clip exon overlap when both exon leak in the intron */
      if (0 && oldX1 && map == oldMap)  
	{
	  int dx = oldX2 - x1 + 1 ;

	  if (dx > 0)  /* overlap on the read */
	    {  /* check the overlap exists also in the target */
	      if ((a2 - a1)*(oldA2 - oldA1) < 0)  /* opposite strands */
		dx = 0 ;
	      if (a1 < a2 && a1 > oldA2)
		dx = 0 ;   /* we already have a hole */
	      if (a1 > a2 && a1 < oldA2)
		dx = 0 ;
	    }
	  if (dx > 0 && oldX2 > x1) /* drop */
	    continue ;
	    
	  if (dx > 0) /* reclip */
	    {	    
	      x1 += dx ;
	      if (a1 < a2) a1 += dx ;
	      else  a1 -= dx ;
	    }
	}
      /* memorize BEFORE rejecting on the basis of the strand */ 
      oldX1 = x1 ; oldX2 = x2 ; oldA1 = a1; oldA2 = a2 ; oldScore = score ; oldMap = map ; oldtc = target_class ;

      /***********************************/
      /* e1 e2 will be used to position the start of the reads on the template */
      e1 = e2 = 0 ;
      if (strcmp (probeBuf, oldProbeBuf) &&
	  (x1 < 10  || (wall1[0] && wall1[0] == ace_upper(wall1[0])))
	  ) 
	{
	  e1 = a1 ; e2 = a2 ; 
	  if (e1 < e2 && e2 > e1 + EndLength) e2 = e1 + EndLength ;
	  if (e1 > e2 && e2 < e1 - EndLength) e2 = e1 - EndLength ;
	}
      /***********************************/
      /* 2024_01_12: in non stranded case, trust the strand of the annotated gene */
      if (! firstRead) {int dummy = a1 ; a1 = a2 ; a2 = dummy ;}
      if (inTranscript && sx->strand && a1 > a2) continue ;
      if (inTranscript && sx->antistrand && a1 < a2) continue ;
      if (sx->antistrand) {int dummy = a1 ; a1 = a2 ; a2 = dummy ;}
      if (inTranscript && sx->ventilate && ! sx->strand && ! sx->antistrand)
	{
	  if (a1 > a2) {int dummy = a1 ; a1 = a2 ; a2 = dummy ;}
	}

      if (sx->maxErr && nErr > sx->maxErr && 100 * nErr > sx->maxErrRate * ali)
	continue ;

      /************** ENDS ***************/
      strcpy (oldProbeBuf, probeBuf) ;

      /*  2015_06_12 : getting this right is amazingly difficult !
       * E [LR] {FR} : E = ends, LR = left or right end of the genebox in the orientation of the genome
       *                          FR = the gene is on the forward or the reverse strand of the genome
       *  this applies to reads mapped directly to the genome
       *  for reads mappaed to a gene, the names stay if the gene maps to the plus strand
       *  both letters LR  and FR are exchanged if the gene maps to the negative starnd 
       */
      if (e1 < e2) 
	{
	  iEfr = 0 ; Efr = "EL" ;
	  if (a1<a2)  { iEfr = 0 ; Efr = "ELF" ; }
	  else            { iEfr = 1 ; Efr = "ELR" ; }
	}
      else
	{ 
	  int e0 = e1 ; e1 = e2 ; e2 = e0 ;
	  iEfr = 2 ; Efr = "ER" ; 
	  if (a1>a2)  { iEfr = 3 ; Efr = "ERR" ; }
	  else            { iEfr = 2 ; Efr = "ERF" ; }
	}
      

      /* missmatches and overhangs, reported in col 14 to 21 are ignored in this program */

      if (1)
	{
	  int mateAliDirect  = 0, mateToAliDirect  = 0,  mateScoreDirect = 0 ; /* in case they are available in cols 23 24 */
	  
	  aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; /* list of missmatches in read frame */
	  aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; /* list of missmatches in target frame */
	  aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; if (! ccp) ccp = "-" ; strncpy (wall1, ccp, 120) ;
	  aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; if (! ccp) ccp = "-" ; strncpy (wall2, ccp, 120) ;
	  aceInStep (ai, '\t') ; ccp = aceInWord (ai) ;	/* target prefix */
	  aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; /* target suffix */
 
	  chainFrom = chainTo = 0 ;
	  /* keep only matched pairs */
	  {
	    int deltaPair = 0 ;
	    aceInStep (ai, '\t') ;  aceInInt (ai, &deltaPair) ; aceInStep (ai, '-') ; 
	      
	    /*                synchronize with snp.c */
	    aceInStep (ai, '\t') ; aceInInt (ai, &mateScore) ; aceInStep (ai, '-') ; 
	    aceInStep (ai, '\t') ; aceInInt (ai, &mateAli) ; aceInStep (ai, '-') ; 
	    aceInStep (ai, '\t') ;  aceInInt (ai, &mateToAli) ; aceInStep (ai, '-') ; 
	    
	    mateAli = mateAliDirect ; mateScore = mateScoreDirect ; mateToAli = mateToAliDirect ; /* if I find them directly that is ok */
	    
	    aceInStep (ai, '\t') ; ccp = aceInWord (ai) ; 
	    if (ccp && !strcmp(ccp, "chain"))
	      {
		aceInStep (ai, '\t') ;  aceInInt (ai, &chainNumber) ; aceInStep (ai, '-') ; 
		aceInStep (ai, '\t') ;  aceInInt (ai, &chainFrom) ; aceInStep (ai, '-') ; 
		aceInStep (ai, '\t') ;  aceInInt (ai, &chainTo) ; aceInStep (ai, '-') ; 
	      }
	  }
	}

      if (tmult <= 0) tmult = 1 ;
      weight = mult/((float)tmult) ;
 
      /************** partials ***************/
      /* check partial */
      p1 = p2 = 0 ;
      
      /* quality control */
      lastAli = ali ; lastToAli = toali ; lastScore = score ;
      if (newProbe) firstX1 = x1 ;
      if (x1 < 30) toali -= x1 ;
      isPartial = FALSE ;
      if (partial && ! noPartial)  /* may happen, the partial option also serves to export the ends */
	{                        /* the partial are just a supplementary info */
	  if ( (ali > 25 || ali > sx->minAliLength || 100 * ali > toali * sx->minAliRate) &&  ali < toali - 30)
	    {
	      if (
		  (! chainTo || x2 == chainTo) &&
		  ! (wall2[0] && wall2[0] == ace_upper(wall2[0])) &&
		  (ali + firstX1 < toali -30) 
		  )
		{ p2 = a2 ; p1 = p2 - (a1 < a2 ? 30 : -30) ; isPartial = TRUE ;}
	      if (
		  (! chainFrom || x1 == chainFrom) &&
		  ! (wall1[0] && wall1[0] == ace_upper(wall1[0])) &&
		  x1 > 30
		  )
		{ p1 = a1 ; p2 = p1 + (a1 < a2 ? 30 : -30) ; isPartial = TRUE ; }
	    }
	}

      if (a1 < a2)
	{ ifr = 0 ; fr = "f" ; }
      else
	{ 
	  int a0 = a1 ; a1 = a2 ; a2 = a0 ; ifr = 1 ; fr = "r" ; 
	  int p0 = p1 ; p1 = p2 ; p2 = p0 ; 
	}
	  
      ao = eo =  po = 0 ;
      int pass = 0 ;
      if (! partial)
       {
	 ao = array (aos, ifr + 4 * tc2 + 8 * pass + 16 * chrom, ACEOUT) ;
	 if (! ao)
	   {
	     ao = aceOutCreate (sx->outFileName
				, messprintf(".%s.%s.%s.minerr%d.%s"
					     , (tc2  ? "remapped" : "mapped")
					     , dictName(sx->remapDict, chrom)
					     , fr, (pass ? sx->minErrRate : 0)
					     , (sx->out == BHIT ? "hits" : (unique ? "u.BG" :"nu.BG") )
					     )
				, sx->gzo, h
				) ;
	     array (aos, ifr + 4 * tc2 + 8 * pass + 16 * chrom, ACEOUT) = ao ;
	     aceOutf (ao,"# %s\n", timeShowNow()) ;
	   }
       }
      if (partial && e1)
	{
	  eo = array (eos, iEfr + 4 * tc2 + 8 * pass + 16 * chrom, ACEOUT) ; 
	  if (! eo)
	    {
	      eo = 
		aceOutCreate (sx->outFileName
			      , messprintf(".%s.%s.%s.minerr%d.%s.%s"
					   , (tc2  ? "remapped" : "mapped")
					   , dictName(sx->remapDict, chrom)
					   , Efr, (pass ? sx->minErrRate : 0)
					   , (non_unique ? "nu" : "u")               /* by defult they are calld .u */ 
					   , (sx->out == BHIT ? "hits" : "BG")
					   )
			      , sx->gzo, h
			      ) ;
	      array (eos, iEfr + 4 * tc2 + 8 * pass + 16 * chrom, ACEOUT) = eo ;
	      aceOutf (eo,"# %s\n", timeShowNow()) ;
	      aceOutf (eo, "trackName type=bedGraph\n") ;
	    }
	}

      if (partial && p1 && isPartial)
	{
	  po = array (pos, ifr + 4 * tc2 + 8 * pass + 16 * chrom, ACEOUT) ;
	  if (!po)
	    {
	      po = 
		aceOutCreate (sx->outFileName
			      , messprintf(".%s.%s.%s.minerr%d.%s"
					   , (tc2  ? "remapped" : "mapped")
					   , dictName(sx->remapDict, chrom)
					   , fr, (pass ? sx->minErrRate : 0)
					   , (sx->out == BHIT ? "pp.hits" : "pp.BG")
					   )
			      , sx->gzo, h
			      ) ;
	      array (pos, ifr + 4 * tc2 + 8 * pass + 16 * chrom, ACEOUT) = po ;
	      aceOutf (po,"# %s\n", timeShowNow()) ;
	      aceOutf (po, "trackName type=bedGraph\n") ;
	    }
	}
      
      if (sx->out == BHIT)
	{
	  if (ao) 
	    { 
	      aceOutf (ao, "%s\t%d\t%d\t%f\t\t",  dictName(sx->mapDict, map), a1, a2, weight) ; 
	      aceOutf (ao, "%s\n", buf) ;
	    }
	  if (eo) 
	    { 
	      aceOutf (eo, "%s\t%d\t%d\t%f\t\t",  dictName(sx->mapDict, map), e1, e2, weight) ; 
	      aceOutf (eo, "%s\n", buf) ;
	    }
	  if (po) 
	    { 
	      aceOutf (po, "%s\t%d\t%d\t%f\t\t",  dictName(sx->mapDict, map), p1, p2, weight) ; 
	      aceOutf (po, "%s\n", buf) ;
	    }
	}
      if (sx->out == BG)
	{
	  /* a1 in BG is zero based */
	  a1-- ;
	  if (ao && a1 >= 0) aceOutf (ao, "%s\t%d\t%d\t%f\n",  dictName(sx->mapDict, map), a1, a2, weight) ; 
	  
	  e1-- ;
	  if (eo && e1 >= 0) aceOutf (eo, "%s\t%d\t%d\t%f\n",  dictName(sx->mapDict, map), e1, e2, weight) ; 
	  
	  p1-- ;
	  if (po && p1 >= 0) aceOutf (po, "%s\t%d\t%d\t%f\n",  dictName(sx->mapDict, map), p1, p2, weight) ; 
	}
    }
  /* close the multiple output files */
  for (tc = 0 ; tc < arrayMax (aos) ; tc++)
    {
      ao = array (aos, tc, ACEOUT) ; array (aos, tc, ACEOUT) = 0 ;
      ac_free (ao) ;
    }
  for (tc = 0 ; tc < arrayMax (eos) ; tc++)
    {
      eo = array (eos, tc, ACEOUT) ; array (eos, tc, ACEOUT) = 0 ;
      ac_free (eo) ;
    }
   for (tc = 0 ; tc < arrayMax (pos) ; tc++)
    {
      po = array (pos, tc, ACEOUT) ; array (pos, tc, ACEOUT) = 0 ;
      ac_free (po) ;
    }
  
  tc = oldScore ; /* for compailer happiness */
  sleep (3) ; /* we seem to get incomplete .gz files in this export */
  ac_free (h) ; 
} /* sxVentilate */

/***********/

static void sxStrandShift (WIGGLE *sx)
{
  AC_HANDLE h = 0 ;
  WIGGLE sxf, sxr ; 
  ACEOUT ao ;

  /* prepare 2 WIGGLE structures */
  memcpy (&sxf, sx, sizeof (WIGGLE)) ;
  memcpy (&sxr, sx, sizeof (WIGGLE)) ;
  sxf.h = h ;
  sxr.h = h ;

  ao = aceOutCreate (sx->outFileName, ".strand_shift.txt", sx->gzo, h) ;

  sxf.dict = dictHandleCreate (100000, h) ;
  sxf.aaa = arrayHandleCreate (100, Array, h) ;
  
  sxf.ai = aceInCreate (sx->strandShift_f, sx->gzi, h) ;
  if (! sxf.ai) messcrash ("Cannot file open %s, sorry", sx->strandShift_f) ;

  sxr.dict = dictHandleCreate (100000, h) ;
  sxr.aaa = arrayHandleCreate (100, Array, h) ;
  
  sxr.ai = aceInCreate (sx->strandShift_r, sx->gzi, h) ;
  if (! sxr.ai) messcrash ("Cannot file open %s, sorry", sx->strandShift_r) ;

  /* parse the 2 wiggles */
  sxWiggleParse (&sxf, 0, 0) ;
  sxWiggleParse (&sxr, 0, 0) ;

  /* compare */
  {
    Array aaf, aar ;
    int i, ii, jj, iMax, k, iLimitMax = 0, iaaa ;
    int dx, dx0, dxmax = 1+sx->strandShift_max/sx->out_step ;
    double z, uu, vv, u1, v1, uv[dxmax], u2[dxmax], nu ;
    WIGGLEPOINT *z1p, *z2p ;
    int *limitp, limits[] = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 10000000, -1} ;
    long int cumul[100000], any[100], both[100], just[100] ;

    memset (cumul, 0, sizeof(cumul)) ;
    memset (both, 0, sizeof(both)) ;
    memset (just, 0, sizeof(just)) ;
    memset (any, 0, sizeof(any)) ;

    nu = uu = vv = u1 = v1 = 0 ; 
    for (dx = 0 ; dx < dxmax ; dx++)
      uv[dx] = u2[dx] = 0 ;

    for (iaaa = 0 ; iaaa < arrayMax (sxf.aaa) ; iaaa++)
      {
	if (0 && iaaa != 2) continue ;
	aaf = array (sxf.aaa, iaaa, Array) ;  /* i do not know why it has to be 2 */
	aar = array (sxr.aaa, iaaa, Array) ;
	ii = aaf ? arrayMax (aaf) : 0 ; 
	jj = aar ? arrayMax (aar) : 0 ;
	iMax = ii < jj ? ii : jj ; iMax -= dxmax ;
	if (iMax <= 0) continue ;
	
	z1p = arrp(aaf,0,WIGGLEPOINT) ; z2p = arrp(aar,0,WIGGLEPOINT) ;
	dx0 = (z2p->x - z1p->x)/sx->out_step ;
	if (dx0 > 0) 
	  { z1p += dx0 ; iMax -= dx0 ; }
	if (dx0 < 0) 
	  { z2p += dx0 ; iMax += dx0 ; }
	for (ii = 0 ; ii < iMax ; z1p++, z2p++, ii++)
	  { 
	    z = z1p->y ; uu += z*z ; u1 += z ;  nu++ ;
	    if (z > 0) 
	      for (dx = 0 ; dx < dxmax ; dx++)
		{ u2[dx] += z * ((z1p+dx)->y) ; uv[dx] += z * ((z2p+dx)->y) ; }
	    z = z2p->y ; vv += z*z ; v1 += z ;
	    z = z1p->y +  z2p->y ; 
	    for (i = 0, limitp = limits ; *limitp > -1 ; i++, limitp++)
	      if (z >= *limitp)
		{
		  /* at various thresholds, compute the number of position at a given strand percentage */
		  k = .49 + 100 * z1p->y/z ; 
		  cumul[k + 200 * i] += sx->out_step ;
		  if (i > iLimitMax) iLimitMax = i ;
		  if ((k >= 30 && k <= 70) || (z1p->y >= 5 && z2p->y >= 5))
		    both[i] +=  sx->out_step ;
		  if (k < 1 || k > 99) just[i] +=  sx->out_step ;
		  any[i] +=  sx->out_step ;

		}
	  }
      }
    if (0)
      {  /* we cannot substract the mean values, they are too close to zero, sorry */
	uu -= u1 * u1 / nu ; vv -= v1 * v1 / nu ;
	z = sqrt (uu * vv) ;
	for (dx = 0 ; dx < dxmax ; dx++)
	  {
	    uv[dx] = (uv[dx] - u1 * v1 / nu) / z ;
	    u2[dx] = (u2[dx] - u1 * u1 / nu) / uu ;
	  }
      }
    else
      {
	z = sqrt (uu * vv) ;
	for (dx = 0 ; dx < dxmax ; dx++)
	  {
	    uv[dx] = (uv[dx]) / z  ;
	    u2[dx] = (u2[dx]) / uu ;
	  }
      }	  

    aceOutf (ao, "# %s\n", timeShowNow()) ;
    aceOutf (ao, "# Autocorrelation of the wiggle on the top strand, used as control\n") ;
    aceOutf (ao, "# Trans correlation of the 2 strand showing the average lag of the minus strand wiggle\n") ;
    aceOutf (ao, "# Distance\tCis autocorrelation\tTrans correlation\t\n") ;
    for (dx = 0 ; dx < dxmax ; dx++)
      aceOutf (ao, "%d\t%g\t%g\n", sx->out_step * dx, u2[dx], uv[dx]) ;
    aceOutf (ao, "\n\n") ;

    ao = aceOutCreate (sx->outFileName, ".strand_coverage_per_threshold.txt", sx->gzo, h) ;
    aceOutf (ao, "# %s\n", timeShowNow()) ;
    aceOutf (ao, "# Histogram above various thresholds of the percentage of coverage on the plus strand\n") ;
    aceOutf (ao, "# Percent +") ;
    for (i = 0, limitp = limits ; *limitp > -1 && i <= iLimitMax ; i++, limitp++)
      aceOutf (ao, "\tthr %d", *limitp) ;
    for (k = 0 ; k <= 100 ; k++)
      {
	aceOutf (ao, "\n%d", k) ;
	for (i = 0, limitp = limits ; *limitp > -1 && i <= iLimitMax ; i++, limitp++)
	  aceOutf (ao, "\t%ld", cumul[k + 200 * i]) ;
      }
    aceOutf (ao, "\nTotal") ;
    for (i = 0, limitp = limits ; *limitp > -1 && i <= iLimitMax ; i++, limitp++)
      aceOutf (ao, "\t%ld", any[i]) ;
    aceOutf (ao, "\nOnly_read_on_one_strand") ;
    for (i = 0, limitp = limits ; *limitp > -1 && i <= iLimitMax ; i++, limitp++)
      aceOutf (ao, "\t%ld", just[i]) ;
    aceOutf (ao, "\nReliably_read_on_both_strands") ;
    for (i = 0, limitp = limits ; *limitp > -1 && i <= iLimitMax ; i++, limitp++)
      aceOutf (ao, "\t%ld", both[i]) ;
    /* 
       any
       Only read on single strand (0 or 100/any)
       Reliably read on both strand (30% to 70% or at least 5 on each strand)
    */
    aceOutf (ao, "\n\n") ;
  }

  ac_free (h) ;
  return ;
} /* sxStrandShift */

/* this method does not work well, if gives a compression rate 
 *  on 60% or 70% relative to BF.gz
 */

/**********************************************************/
/*************  predictor/compressor, modelled on CTF *****/
/**********************************************************/

static Array wiggle_compress_decorrelate (WIGGLE *sx, Array aa)
{
  int ii ;
  WIGGLEPOINT *wp ;
  int u1, u2, u3, u4, z , *zp ;
  int predictionMode = sx->BF_predictor ;
  int iMax = arrayMax (aa) ;
  Array a = arrayHandleCreate (iMax, int, sx->h) ;

  u1 = u2 = u3 = u4 = 0 ; 
  array (a, iMax -1, int) = 0 ; /* make room */
  for (ii = 0, wp = arrp (aa, ii, WIGGLEPOINT), zp = arrp (a, 0, int) ; ii < iMax ; ii++, wp++, zp++) 
    {
      float y = wp->y ;
      int u0 = y + .0001 ;
      switch (predictionMode)
	{
	case 0: z = 0 ;
	case 1: z = u1 ; break ; /* predict flat, transmit derivative */
	case 2: z = 2*u1 - u2 ; break ; /* predict line trans dd2 */
	case 3: z = 3*u1 - 3*u2 + u3 ; break ; /* predict parabole */
	case 4: z = 4*u1 - 6*u2 + 4*u3 - u4; break ; /* overpredict ! */
	}
      u4 = u3 ; u3 = u2 ; u2 = u1 ; u1 = u0 ;
      *zp = u0 - z ;
    } 
  return a ;
} /* wiggle_compress_decorrelate */

/**********************************************************/

static void wiggle_compress_write_one (WIGGLE *sx, Array aa)
{
  ACEOUT ao = sx->ao ;

  if (1)
    {
      Array a = wiggle_compress_decorrelate (sx, aa) ;
      int i, iMax = arrayMax (a) ;
      int *ip = arrp (a, 0, int) ;

      if (1)
	{
	  for (i = 0 ; i < iMax ; i++, ip++)
	    {
	      if (  (*ip > -126 && *ip < 126))
		{
		  char cc = *ip ;
		  aceOutf (ao, "%c", cc) ;
		}
	      else if (  (*ip > -126 * 256 && *ip < 126 * 256))
		{
		  short cc = *ip ;
		  aceOutf (ao, "%c%c%c", 127,cc/256,cc & 0xff) ;
		}
	      else
		aceOutf (ao, "%d\n", *ip) ;
	    }
	}
      else
       aceOutBinary (ao, ip, iMax * sizeof(int)) ;
    }
} /* wiggle_compress_write_one */

/**********************************************************/

static void wiggle_compress_write (WIGGLE *sx)
{
  Array aa ;
  int remap ;
  
  for (remap = 0 ; remap < arrayMax (sx->aaa) ; remap++)
    {
      aa =  arr (sx->aaa, remap, Array) ;
      if (arrayExists (aa) && arrayMax (aa))
	wiggle_compress_write_one (sx, aa) ;
    }
  return ;
} /* wiggle_compress_write */

/**********************************************************/

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char *error)
{
  fprintf (stderr, "// Usage: wiggle \n") ;
  fprintf (stderr, "// Example:  wiggle -I BHIT -O BV -o foo/bar -mapon sections.map\n") ;  
  fprintf (stderr,
	   "// author: mieg@ncbi.nlm.nih.gov\n"
	   "//  no guarantee, no restriction whatsoever\n"
	   "// Input\n"
	   "//   -i input_file : [default: stdin] data file to analyze\n"
	   "//   -I input_format : [default: BV = tab-delimited hits] format of the input file, defined below\n"
	   "// Output\n"
	   "//   -o output_file_prefix : [default: stdout] prefix for the name of the processed sequence\n"
	   "//   -O output_format : [mandatory] format of the output file, defined below\n"
	   "// \n"
	   "// gzip (recommended)\n"
	   "//   -gzi -gzo : gzi to decompress the input file, gzo to compress the output file\n"
	   "// \n"
	   "// Input formats\n"
	   "//   In BF,BV,BG formats, a pure tab delimited file without header is acceptable\n"
#ifdef JUNK
	   "//   AW : AceView Wiggle format, compressed, machine independant, binary format\n"
#endif
	   "//   BF : one column fixed step UCSC Bed format : gives the height at every idx position\n"
	   "//     -target <name> , the name of the target, i.e. chrom7, MYH9\n"
	   "//     -in_step <int> , default 10, distance between successive input values in 'y' format\n"
	   "//     -in_span <int> , default in_step, coverage of each data point\n"
	   "//     -in_x1 <int> , default 1, 1-based position of the first input value in 'y' format\n" 
	   "//     These parameters are reset on each line labelled\n"
	   "//     fixedStep chrom=<name> start=<int> step=<int>\n"
	   "//   BV : two columns variable step UCSC Bed format : position 'tab' value\n"
	   "//     -target <name> , the name of the target, i.e. chrom7, MYH9\n"
	   "//     -in_span <int> , default 1, coverage of each data point\n"
	   "//     This parameter is reset on each line labelled\n"
	   "//     variableStep chrom=<name> span=<int>\n"
	   "//   BG : UCSC tab delimited bedGraph format : target x1 x2 value \n"
	   "//     In BG format, x1 is zero-based but x2 is 1 based, so the orientation is always known\n"
	   "//\n"
	   "//   BHIT :13 columns tab delimited table (COUNT/manip/tissue/*.hits): \n"
	   "//         as exported by clipalign and by bestali\n" 
	   "//         tag,score,tag-multiplicity, clipped_length,aligned_length,x1,x2\n"
	   "//         target_class,gene,unicity,target,a1,a2\n"
	   "//         additional columns are ignored, lines starting with # are ignored\n"
	   "//     tag is a short sequence identifier, contributing to the wiggle\n"
	   "//     score is the score of this hit\n"
	   "//     multiplicity  is the number of time this exact tag has been sequenced\n"
	   "//     length and ali are the clipped length of the tag, and the length of the alignment\n"
	   "//     x1,x2 are the 1-based coordinates of the alignment on the tag (ignored)\n"
	   "//     class is the class of the target, i.e. mito, genome, av, RefSeq, EBI, HINV\n"
	   "//     gene is an otional gene identifier\n"
	   "//     unicity is the number of genes or genome position at equal score for this tag\n"
	   "//     a1,a2 are the 1-based coordinates on the target covered by the tag\n"
	   "//\n"
	   "//     By default the hits on the 2 strands of the target are combined\n"
	   "//     unless you specify\n" 
	   "//        -strand : only consider in the hit file the hits on the positive strand\n"
	   "//        -antistrand : only consider the hits on the negative strand\n"
/* not implemented, sorry
	   "//        -stranded : requires -o option, create 2 stranded output files called .f and .r\n" 
*/
	   "//     By default all the tags in the hit file are combined\n"
	   "//     unless you specify\n"
	   "//        -select sfileName : only consider in the hit file the sequences listed there.\n"
	   "//        -reject rfileName : reject from the hit file the sequences listed there.\n"
	   "//     -target_class className : only keep the lines with the correct class in column 8\n"
	   "//     -unique : applies to BHIT input, only consider unique hits within the target class: col10 == 1\n"
	   "//     -nu |-non_unique : idem but only consider non-unique hits within the target class: col10 > 1\n"
	   "//     -force_unique : reset mult =1 before checing unicity\n"
	   "//     -partial : idem but only consider partial or orphan or incompatible pairs\n"
	   "//         partial implies unique, i.e. forget the multi-aligned partial reads\n" 
           "//     -noPartial : partial are treated as complete\n"
	   "//         defaul in nu case so that u + pp + nu ar additive\n"
	   "//         useful for nanopore or pacbio runs where all ali are partial\n"
	   "//     -delta <int>: \n"
	   "//         in this case all hits are shortened by delta at both ends \n"
	   "//         this can be used to create a 2*delta gap around unsupported introns\n"
	   "//\n"
	   "// Output formats\n"
	   "//     [-trackName name] , overruled by target names specifed in the input file"
	   "//   BF : one column fixed step UCSC Bed format : gives the height at every idx position\n"
	   "//     -out_step <int> , default 10, distance between successive input values in 'y' format\n"
	   "//     -BF_compressor <int n>:  n = 1,2,3 or 4\n"
	   "//     -BF_predictor <int n>:  n = 1,2,3 or 4\n"
	   "//        degree of the polynomial predictor used to compress the BF format\n"
	   "//   BV : two columns variable step UCSC Bed format : position 'tab' value\n"
	   "//     -out_span <int> , default 1, coverage of each data point\n"
	   "//     variableStep chrom=<name> span=<int>\n"
	   "//   BG : UCSC tab delimited bedGraph format : target x1 x2 value \n"
	   "//     In BG format, x1 is zero-based but x2 is 1 based, so the orientation is always known\n"
	   "//   AM : "
	   "//     the wiggle is exported in .ace format as class:mRNA\\nWiggle 'name' coord number-of-tags\n"
	   "//   AG : -feature name\n"
	   "//     the wiggle is exported in .ace format as Sequence target\\nFeature 'name' coord coord number-of-tags 'trackName'\n"
	   "//   Count : \n"
	   "//     coverage statistics are exported on stdout\n"
	   "//   -cumul : \n"
	   "//     coverage statistics are exported on  <-out>.cumul\n"
           "//   -peaks : export the peaks above mincover with area above 1.5 * minCover * width\n"
	   "//     peaks are exported on <-out>.peaks as a 4 columns table\n"
	   "//     target(e.g. the target of the wiggle) x1 x2 max area\n"
	   "//   -multiPeaks <int n> -wiggle1 <file1> -wiggle2 <file2> -swiggle1 <sfile1> -swiggle2 <sfile2> -stranding <float s> -minCover <int c>\n"
	   "//     file 1/2 are wiggle files (of type -I, optionally .gzipped) representing sense and antisense alignments\n"
	   "//     stranding [default 100] is the percentage of reads mapping to the top strand, for example 98.5, or equivalently 1.5\n"
	   "//     the idea is first to substract a fraction (2*s)/100 of file2 from file1,\n"
	   "//     to screen the artefact that a fraction a large antisense peak of a not perfectly stranded protocol will leak to the sense strand,\n"
	   "//     then to export coverons representing region covered at least c times, as stepwise-plateaux with granularity n-fold\n" 
	   "//     This option is typically run first of w1=forward wiggle w2=reverse wiggle, then w1=R w2=F, yieding clean strand specific exons\n"
	   "//     If the control stranded wiggles swiggle[12] are provided, they control the final strand ratio of the (non stranded) wiggles\n"
	   "//     In that case the stranding should be that of the control wiggles\n"
	   "//   -transcriptsEnds <file name prefix ff> -I <type> [-gzi] -stranding <int> -minCover <int c> : detect transcripts-ends\n"
	   "//     Analyze the files ff.Exx.type[.gz], where Exx = ELF, ELR, ERF, ERR, representing the left/rigth ends of genes on the F/R strands\n"
	   "//     Export the corresponding 4 types of coverons, with coverage at least an optimized fraction of c, and high E[L/R]x contrast\n" 
	   "//     The same scaled sreening strategy is applied in -multiPeaks and -transcriptsEnds options\n"  
	   "//   -lengthCoverage : export the local aligned length\n"
	   "//     Requires -I BH | BG, the local length is available only if the input is in BG or BH format.\n"
	   "//     Export a wiggle every step base where each read contibutes its aligned length\n"
	   "//     Dividing by the local coverage yields the local averaged aligned length\n"
	   "//     Dips are expected close to difficulties (Repeated regions, SNPs, deletions, RNA editing ...)\n"
	   "//     Since numbers are exported, they are additive and the average can be computed for a group of files\n"
	   /*************
           "//   -wiggleRatio <int L> -wiggleRatioDamper <int D> -wiggle1 <file1> -wiggle2 <file2>peaks :\n"
	   "//     file1 and file2 are 2 wiggle files of type -I\n"
	   "//     D (default 100) damps the ratio to avoid detecting fluctuations of small values\n" 
	   "//     the program computes at every outstep the log of the damped ratio log(a+D/b+D)\n"
	   "//     peaks above +-L are exported on <-out>.peaks as a 4 columns table\n"
	   "//     target(e.g. the target of the wiggle) x1 x2 max area ratio(saturated at 1000)\n"
	   **************/
	   
	   "//\n"   
	   "// -maxCover number : discard positions above this coverage  \n"
	   "// -minCover number : discard positions below this coverage  \n"
	   "//     example: -minCover 100 -O BV will export a shadow of all 'exons' above 100\n"
	   "// -scale float : multiply all hights by scale\n"
	   "//\n"   
	   "// -mapon map_file -wiggleDir name \n"
	   "//    Remap the hits on the given targets, ignore the remainder\n"
	   "//       map_file has 1 or 3 columns: only export the hits corresponding to these targets\n"
	   "//       map_file has 7 columns: target_class target x1 x2 chrom a1 a2, hits to target segment x1,x2 are remapped to chrom a1,a2,\n"
	   "//       To map on the genome the hits to each exon of a spliced mRNA give one line per exon.\n"
	   "//    -mapOnChrom chrom : if specified, keep only the lines in map_file where column4 == chrom\n" 
           "//    -wiggleDir dirName:  the wiggle is exported in wiggle format in the externalDirectory,\n"
	   "//          the file is named name.slx.\n"
	   "//\n"
	   "// -mrnaWiggle exon_coordinates_file : requires -o option\n"
	   "//       Export a wiggle per target,\n"
	   "//       Target has 1 column of mrna names, corresponding to col 1 of the hits file\n"
	   "//       hits has 3 columns: map a1 a2, matching the appelations in target\n"
	   "//       A value is exported each time there is a change in coverage\n"
	   "//       the wiggle is exported in .ace format as class:mRNA tag:Solexa \n"
	   "//            i.e.: Solexa 'trackName'   coord-on-mRNA     number-of-tags\n"
	   "// -gauss sigma : sigma beeing an integer, in bp.\n"
	   "//    Smooth the wiggle by convolution with a Gaussian  exp(-x^2/2*sigma^2)\n"
	   "//    this is mathematically equivalent to filtering out the high frequencies in Fourier space\n"
	   "// -strand_shift max -ssf file1.f -ssr file2.r : auto-correlation of forward and reverse wiggles\n"
	   "//    In ChIP-seq for example the reads on the 2 strands are shifted by the effective length of the library\n"
	   "//    This function reports the correlation (cosine in the scalar product) for shifts up to max bp\n"
	   "// -ventilate [-hierarchic] [-minAliRate number -minAliLength number -naxSnp number -maxErrRate number] [-pair fragmentLength]\n"
	   "//    Implies -I BHIT -O BHIT, requires -o, splits a HIT file in one BG file per chromosome\n" 
	   "//    in -hierarchic mode, a single target per target_class is exported , excluding the genome\n"
	   "//    reject ali shorter than minAliLength [default 70] and below minAliRate [default 0]\n"
	   "//    reject ali shorter above maxErr and above maxErrRate [default 0 and 0%%]\n"
	   "//    pair rejects all reads which are not pair compatible\n"
	   );
  if (error)
    fprintf (stderr, "// ERROR: %s\n", error) ;
  exit (1) ;
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  WIGGLE sx ;
  const char *ccp ;
  char iType[13] ;
  AC_HANDLE h = 0 ;

  freeinit () ; 
  messErrorInit (argv[0]) ;

  h = ac_new_handle () ;
  memset (&sx, 0, sizeof (WIGGLE)) ;
  sx.h = h ;
  
  if (getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;
  
  
  sx.ventilate = getCmdLineBool(&argc, argv, "-ventilate") ;
  sx.cumul = getCmdLineBool(&argc, argv, "-cumul") ; 
  sx.peaks = getCmdLineBool (&argc, argv, "-peaks") ;
  if (getCmdLineInt (&argc, argv, "-multiPeaks", &(sx.multiPeaks)))
    sx.peaks = TRUE ;
  if (getCmdLineOption(&argc, argv, "-transcriptsEnds", &(sx.transcriptsEndsFileName)))
    sx.peaks = TRUE ;
  sx.minAliLength = 100 ; sx.minAliRate = 50 ; /* defaults */
  getCmdLineInt(&argc, argv, "-minAliRate", &(sx.minAliRate)) ;
  getCmdLineInt(&argc, argv, "-minAliLength", &(sx.minAliLength)) ;
  
  sx.maxErr = sx.maxErrRate = 0 ; /* defaults */
  getCmdLineInt(&argc, argv, "-maxErr", &(sx.maxErr)) ;
  getCmdLineInt(&argc, argv, "-maxErrRate", &(sx.maxErrRate)) ;
  getCmdLineInt(&argc, argv, "-minErrRate", &(sx.minErrRate)) ;

  getCmdLineInt(&argc, argv, "-BF_compressor", &(sx.BF_compressor)) ;
  getCmdLineInt(&argc, argv, "-BF_predictor", &(sx.BF_predictor)) ;

  getCmdLineInt(&argc, argv, "-pair", &(sx.pair)) ;
  sx.wiggleScale1 = 1.0 ; sx.wiggleScale2 = 0 ;
  if (getCmdLineFloat(&argc, argv, "-stranding", &(sx.wiggleScale2)))
    {
      float x = sx.wiggleScale2 ;
      if (x > 20 && x < 80)
	sx.wiggleScale2 = 1 ; /* add up the 2 strand */
      else
	sx.wiggleScale2 = -2 * (100 - x) / 100.0 ;
    }

  /* strand shift */
  getCmdLineInt(&argc, argv, "-strand_shift", &(sx.strandShift_max)) ;
  getCmdLineOption (&argc, argv, "-ssf", &(sx.strandShift_f));
  getCmdLineOption (&argc, argv, "-ssr", &(sx.strandShift_r));

  /* input output selection */

  getCmdLineOption (&argc, argv, "-i", &(sx.inFileName)) ;
  getCmdLineOption (&argc, argv, "-o", &(sx.outFileName)) ;

  if (getCmdLineOption (&argc, argv, "-I", &ccp))
    {
      strncpy (iType, ccp, 12) ;
      if (! sxCheckFormat ("I", &sx.in, ccp, sx.fileType))
	usage (messprintf ("Unknown format in option -I %s", ccp)) ;
      if (sx.in == COUNT)
	usage (messprintf ("Count is not an acceptable input format in option -I %s", ccp)) ;
    }
  else
    usage ("missing input format option -I") ;
  if (getCmdLineOption (&argc, argv, "-O", &ccp))
    {
      if (!sxCheckFormat ("O", &sx.out, ccp, sx.fileType))
	usage (messprintf ("Unknown format in option -O %s", ccp)) ;
    }
  else if (sx.ventilate)
    {
      sxCheckFormat ("O", &sx.out, "BHIT", sx.fileType) ;
    }
  else if (!sx.strandShift_max)
    usage ("missing output format option -O") ;


  /* specific input options relative to 'wy' imput format */
  getCmdLineInt(&argc, argv, "-in_step", &(sx.in_step));
  getCmdLineInt(&argc, argv, "-in_span", &(sx.in_span));
  getCmdLineInt(&argc, argv, "-in_x1", &(sx.in_step));
  

  sx.delta = 0 ;
  
  switch (sx.out)
    {
    case BHIT : 
    case BG:
    case AF:  /* allele frequency */
    case COUNT:
      sx.out_step = 10 ;
      break ;
    default: /* BF, BV, AM, AG, AW  */
      sx.out_step = 10 ;
      break ;
    }
  getCmdLineInt(&argc, argv, "-out_step", &(sx.out_step)) ;
  getCmdLineInt(&argc, argv, "-delta", &(sx.delta));
  getCmdLineInt(&argc, argv, "-out_span", &(sx.out_span));

  sx.maxCover = ACEDB_MAXINT ;
  sx.minCover = ACEDB_MAXINT ; sx.minCover *= -1 ; /* in 2 steps to cast to float */
  getCmdLineFloat(&argc, argv, "-minCover", &(sx.minCover));
  getCmdLineFloat(&argc, argv, "-maxCover", &(sx.maxCover));

  getCmdLineInt(&argc, argv, "-gauss", &(sx.gauss));

  getCmdLineOption (&argc, argv, "-title", &sx.title) ;

  getCmdLineOption (&argc, argv, "-select", &sx.selectFileName) ;
  getCmdLineOption (&argc, argv, "-reject", &sx.rejectFileName) ;

  sx.strand = getCmdLineBool (&argc, argv, "-strand") ;
  sx.antistrand = getCmdLineBool (&argc, argv, "-antistrand") ;
  if (0)sx.stranded = getCmdLineBool (&argc, argv, "-stranded") ;

  getCmdLineOption (&argc, argv, "-nmWiggle", &sx.nmWiggle) ;
  sx.mrnaWiggle = getCmdLineBool (&argc, argv, "-mrnaWiggle");
  sx.unique = getCmdLineBool (&argc, argv, "-unique");
  sx.forceUnique = getCmdLineBool (&argc, argv, "-force_unique");
  sx.non_unique = (getCmdLineBool (&argc, argv, "-nu") ||  getCmdLineBool (&argc, argv, "-non_unique") ) ;
  sx.partial = (getCmdLineBool (&argc, argv, "-partial")) ;
  sx.noPartial = getCmdLineBool (&argc, argv, "-noPartial") ;
  if (sx.partial && ! sx.ventilate && sx.out != BG)
    usage ("option -partial is only compatible with -o BG") ;
  if (sx.unique && sx.non_unique)
    usage ("Sorry: options -unique and -non_unique are incompatible") ;
  if (0 && sx.partial && sx.noPartial)
    usage ("Sorry: options -partial and -noPartial are incompatible") ;
  if (sx.partial && sx.non_unique)
    usage ("Sorry: options -partial and -non_unique are incompatible") ;
  if (sx.ventilate && !sx.unique && !sx.non_unique && !sx.partial)
    usage ("Sorry: options -ventilate requires -unique or -non_unique or -partial") ;
  sx.gzi = getCmdLineBool (&argc, argv, "-gzi");
  sx.gzo = getCmdLineBool (&argc, argv, "-gzo");
  if (sx.ventilate && sx.gzo)
    usage ("Sorry: options -ventilate is incompatible with -gzo, it would saturate the computer") ;

  getCmdLineOption (&argc, argv, "-mapon", &(sx.mapFileName));
  getCmdLineOption (&argc, argv, "-mapOnChrom", &(sx.sxxChromosome));
  getCmdLineOption (&argc, argv, "-wiggleDir", &sx.wiggleDir) ;
  sx.hierarchic = getCmdLineBool (&argc, argv, "-hierarchic") ;
  sx.lengthCoverage = getCmdLineBool (&argc, argv, "-lengthCoverage") ;

  /* ACEDB Sequence->Feature .ace format */
  getCmdLineOption (&argc, argv, "-feature", &sx.wiggleFeature) ;
  getCmdLineOption (&argc, argv, "-trackName", &(sx.trackName));
  getCmdLineOption (&argc, argv, "-target_class", &(sx.target_class));
  getCmdLineOption (&argc, argv, "-min_target_class", &(sx.min_target_class));

  /* wiggle ratio detection of differential peaks */
  getCmdLineInt (&argc, argv, "-wiggleRatio", &(sx.wiggleRatio)) ;
  sx.wiggleRatioDamper = 5 ;
  getCmdLineInt (&argc, argv, "-wiggleRatioDamper", &(sx.wiggleRatioDamper)) ;
  getCmdLineOption (&argc, argv, "-wiggle1", &(sx.wiggleFileName1)) ;
  getCmdLineOption (&argc, argv, "-wiggle2", &(sx.wiggleFileName2)) ;
  getCmdLineOption (&argc, argv, "-swiggle1", &(sx.swiggleFileName1)) ;
  getCmdLineOption (&argc, argv, "-swiggle2", &(sx.swiggleFileName2)) ;
 
  if (getCmdLineOption (&argc, argv, "-strategy", &(sx.strategy)) &&
      ! strcasecmp (sx.strategy, "RNA_seq"))
    sx.RNA_seq = TRUE ;
 
  if (argc > 1)
    usage (messprintf ("unknown option %s\n", argv[1]));

  /**********/

  fprintf (stderr, "// start: %s\n", timeShowNow()) ;

  sx.ai = aceInCreate (sx.inFileName, sx.gzi, h) ;
  aceInSpecial (sx.ai,"\t\n") ;
  sx.ao = aceOutCreate (sx.outFileName, messprintf(".%s", sx.fileType), sx.gzo, h) ;
 
  if (0) aceOutf (sx.ao, "// %s\n", timeShowNow()) ;
  
  /****** actual work ********/
  if (sx.strandShift_max)
    sxStrandShift (&sx) ;
  else if (sx.ventilate)
    {
      sxGetMap (&sx) ;
      sxVentilate (&sx) ;
    }
  else
    {
      sx.dict = dictHandleCreate (100000, h) ;
      sx.aaa = arrayHandleCreate (100, Array, h) ;
  
      sxGetSelection (&sx) ;
      sxGetMap (&sx) ;
      
      ac_free (sx.ai) ;
      if (sx.transcriptsEndsFileName)  /* transcriptsEnds */
	{
	  int pass ;

	  for (pass = 0 ; pass < 4 ; pass++)
	    {
	      const char *elf, *elr, *erf, *err ;
	      /* z = ax + by   ==   b (a/b x + y) */
	      /* deduct the other end from the good end */
	      switch (pass)
		{
		case 0: /* get ELF */
		  elf = "ELF" ; erf = "ERF" ; 
		  elr = "ELR" ; err = "ERR" ;
		  break ;
		case 1: /* get ERF */
		  erf = "ELF" ; elf = "ERF" ; 
		  err = "ELR" ; elr = "ERR" ;
		  break ;
		case 2: /* get ELR */
		  elf = "ELR" ; erf = "ERR" ; 
		  elr = "ELF" ; err = "ERF" ;
		  break ;
		case 3: /* get ERR */
		  erf = "ELR" ; elf = "ERR" ; 
		  err = "ELF" ; elr = "ERF" ;
		  break ;
		}
	      sx.aoPeaks = aceOutCreate (hprintf (h, "%s.%s", sx.outFileName, elf), ".transcriptsEnds", sx.gzo, h) ;
	      /* add or remove the scaled opposite strand "other end"*/
	      if (sx.wiggleScale2)
		{
		  sx.ai = aceInCreate (hprintf(h,"%s.%s.%s%s",sx.transcriptsEndsFileName, err, iType, sx.gzi ? ".gz" : ""), sx.gzi, h) ;
		  if (sx.ai)
		    {
		      aceInSpecial (sx.ai,"\t\n") ;
		      sxWiggleParse (&sx, 0, 0) ;      /* parse x */
		      ac_free (sx.ai) ;

		      if (sx.wiggleScale2 != 1)
			sxWiggleScale (&sx,  sx.wiggleScale2) ;
		    }
		}

	      /* add the good strand "other end" */
	     sx.ai = aceInCreate (hprintf(h,"%s.%s.%s%s",sx.transcriptsEndsFileName, erf, iType, sx.gzi ? ".gz" : ""), sx.gzi, h) ;
	      if (sx.ai)
		{
		  aceInSpecial (sx.ai,"\t\n") ;
		  sxWiggleParse (&sx, 0, 0) ;      /* parse x */
		  ac_free (sx.ai) ;
		}
	      if (sx.wiggleScale1 != 1)
		sxWiggleScale (&sx,  sx.wiggleScale1) ;
	      if (sx.wiggleScale2 < 0)
		sxWiggleFloor (&sx, 0) ;
	      /* sxWiggleShift (&sx, sx.wiggleRatioDamper) ; */
	      sxWiggleCopy (&sx) ;    /* push the "other end" on the wiggle stack */
	      sx.aaa = 0 ;  sx.aaa = arrayHandleCreate (100, Array, h) ; sxGetMap (&sx) ;

	      /* add or remove the scaled "good end" opposite strand */
	      if (sx.wiggleScale2)
		{
		  sx.ai = aceInCreate (hprintf(h,"%s.%s.%s%s",sx.transcriptsEndsFileName, elr, iType, sx.gzi ? ".gz" : ""), sx.gzi, h) ;
		  if (sx.ai)
		    {
		      aceInSpecial (sx.ai,"\t\n") ;
		      sxWiggleParse (&sx, 0, 0) ;      /* parse x */
		      ac_free (sx.ai) ;
		    }
		  if (sx.wiggleScale2 != 1)
		    sxWiggleScale (&sx,  sx.wiggleScale2) ;
		}

	      /* add the good end good strand */
	     sx.ai = aceInCreate (hprintf(h,"%s.%s.%s%s",sx.transcriptsEndsFileName, elf, iType, sx.gzi ? ".gz" : ""), sx.gzi, h) ;
	      if (sx.ai)
		{
		  aceInSpecial (sx.ai,"\t\n") ;
		  sxWiggleParse (&sx, 0, 0) ;      /* parse x */
		  ac_free (sx.ai) ;
		}
	      if (sx.wiggleScale1 != 1)
		sxWiggleScale (&sx,  sx.wiggleScale1) ;
	      if (sx.wiggleScale2 < 0)
		sxWiggleFloor (&sx, 0) ;
	      
	      sxWiggleEndRatio (&sx) ;

	      sxWiggleExport (&sx) ;
	    }
 	}
      else /* if (!sx.transcriptends)  multiPeaks */
	{
	  /* z = ax + by   ==   b (a/b x + y) */
	  if (sx.swiggleFileName1 || sx.swiggleFileName2)
	    {
	      Array Aaaa = 0, Baaa = 0, Alpha = 0 ;
	      Array aaa1 = 0 ;
	      AC_HANDLE h = 0 ;

	      if (! sx.swiggleFileName2)
		usage ("Missing parameter: -swiggle1 requires -swiggle2") ;
	      if (! sx.swiggleFileName1)
		usage ("Missing parameter: -swiggle2 requires -swiggle1") ;
	      if (! sx.wiggleFileName1)
		usage ("Missing parameter: -swiggle1 requires -wiggle1") ;
	      if (! sx.wiggleFileName2)
		usage ("Missing parameter: -swiggle1 requires -wiggle2") ;

	      /* get a, b  = stranded wiggles, smooth them over sigma = 50 bases
	       * compute alpha = a + 10 / a + b + 20
	       *          beta = b + 10 / a + b + 20
	       *    => alpha + beta = 1, so we only need to compute alpha
	       * then read the non stranded wiggle x, y. 
	       * compute z = x + y
	       *         x' = alpha * z 
	       *         y' =  beta * z 
	       * use x', y' in place of x, y
	       *  since we just want to use x', 
	       *  we can compute in place x, then z = x + y, then x' = alpha * z
	       */
 
	      h = ac_new_handle () ;
	      if (1) /* parse a */
		{
		  sx.ai = aceInCreate (sx.swiggleFileName1, sx.gzi, h) ;
		  sxWiggleParse (&sx, 0, 0) ;      /* parse x */
		  ac_free (sx.ai) ;

		  if (sx.gauss)
		    sxWiggleGauss (&sx) ; 
		  Aaaa = sx.aaa ; sx.aaa = arrayHandleCreate (arrayMax (Aaaa), Array, h) ;
		}
	      if (1) /* parse b */
		{
		  sx.ai = aceInCreate (sx.swiggleFileName2, sx.gzi, h) ;
		  aceInSpecial (sx.ai,"\t\n") ;
		  sxWiggleParse (&sx, 0, 0) ;      /* parse x */
		  ac_free (sx.ai) ;

		  if (sx.gauss)
		    sxWiggleGauss (&sx) ; 
		  Baaa = sx.aaa ; sx.aaa = arrayHandleCreate (arrayMax (Aaaa), Array, h) ;
		}
	      if (1) /* construct the damped ratio alpha */
		Alpha = sxWiggleStrandedRatio (&sx, Aaaa, Baaa, sx.wiggleRatioDamper, h) ; /* => a + 10 / a + b + 10 */
	      if (1) /* parse x */
		{
		  sx.ai = aceInCreate (sx.wiggleFileName1, sx.gzi, h) ;
		  aceInSpecial (sx.ai,"\t\n") ;
		  sxWiggleParse (&sx, 0, 0) ;      /* parse x1 */
		  ac_free (sx.ai) ;
		}
	      if (1) /* preserve the first wiggle and parse again */
		{
		  aaa1 = sx.aaa ;
		  sx.aaa = arrayHandleCreate (arrayMax (Aaaa), Array, h) ;
		  sx.ai = aceInCreate (sx.wiggleFileName1, sx.gzi, h) ;
		  aceInSpecial (sx.ai,"\t\n") ;
		  sxWiggleParse (&sx, 0, 0) ;      /* parse x1 */
		  ac_free (sx.ai) ;
		}
	      if (1) /* parse y, and add into z = x + y */
		{
		  sx.ai = aceInCreate (sx.wiggleFileName2, sx.gzi, h) ;
		  aceInSpecial (sx.ai,"\t\n") ;
		  sxWiggleParse (&sx, 0, 0) ;      /* parse x2 -> x1+x2 */
		  ac_free (sx.ai) ;    
		}
	      if (1) /* compute x' */
		{
		  /* x' = alpha*(x1+x2) OR x1 if (alpha undefined) */
		  sxWiggleMultiplyLocally (&sx, Alpha, aaa1) ;      /* ratio */
		}
	      sxWiggleExport (&sx) ;
	   
	      sx.aaa = 0 ;
	      ac_free (h) ;
	    }
	  else if (sx.wiggleFileName1)
	    { 
	      ac_free (sx.ai) ;
	      sx.ai = aceInCreate (sx.wiggleFileName1, sx.gzi, h) ;
	      aceInSpecial (sx.ai,"\t\n") ;
 
	      sxWiggleParse (&sx, 0, 0) ;      /* parse x */
	      if (sx.wiggleScale2 && sx.wiggleFileName2)
		{
		  float scale1 = sx.wiggleScale1 ;
		  float scale2 = sx.wiggleScale2 ;
		  
		  if (sx.wiggleScale2) 
		    scale1 /= sx.wiggleScale2 ;   
		  if (scale1 != 1)
		    sxWiggleScale (&sx, scale1) ;   /* obtain a/b x */
		  ac_free (sx.ai) ;
		  sx.ai = aceInCreate (sx.wiggleFileName2, sx.gzi, h) ;
		  aceInSpecial (sx.ai,"\t\n") ;
		  sxWiggleParse (&sx, 0, 0) ;     /* add y, obtain a/b x + y */
		  if (scale2 != 1)
		    sxWiggleScale (&sx, scale2) ; /* scale again, obtain ax + by */
		  if (scale2 < 0)
		    sxWiggleFloor (&sx, 0) ;
		}
	      ac_free (sx.ai) ;
	    }
	  else  /* single wiggle */
	    { 
	      sx.ai = aceInCreate (sx.inFileName, sx.gzi, h) ;
	      aceInSpecial (sx.ai,"\t\n") ;
	      sxWiggleParse (&sx, 0, 0) ;
	      ac_free (sx.ai) ;
	    }
	  if (sx.aaa)
	    {
	      sxWiggleGauss (&sx) ; 
	      if (sx.in == BF && sx.BF_predictor)
		wiggle_compress_write (&sx) ;
	      else
		sxWiggleExport (&sx) ;
	    }
	}
    }
      
  /****** clean up *********/
  
  fprintf (stderr, "// done: %s\n", timeShowNow()) ;
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   fprintf (stderr, "// max memory %d Mb\n", mx) ;
  }

  ac_free (sx.h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

 
