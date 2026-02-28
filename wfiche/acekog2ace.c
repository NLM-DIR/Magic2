#include "ac.h"
#include "mytime.h"

#define PREFIX(b,s) ( strncmp(b,s,sizeof(s)-1) == 0 )

static char *getPrefix (ACEIN fi, char *prefix)
{
  char *cp ;
  
  while (aceInCard (fi))
    {
      cp = aceInPos (fi) ;
      if (!cp || strncmp (cp, prefix, strlen(prefix)))
	continue ;
      return cp + strlen(prefix) ; /* the rest of the line */
    }
  return 0 ;
}

static int aceKog2aceParse (ACEIN fi, BOOL isN, DICT *geneDict, Array titles)
{
  int nKantor = 0, subjectTitle, gene, oldGene = -1, link, product ;
  int a1, a2, x1, x2,  b1, y1, eValue, nHit ;
  float bitScore = 0, bestScore = 0 ;
  char *cp, *cq, species ;
  Stack s = stackCreate (1000) ;
  int spec[256], tag[256] ;
  double expect ;

  pushText (s, "-") ;
  for (a1 = 0 ; a1 < 256 ; a1++)
    spec[a1] = 0 ;
  spec['w'] = stackMark (s) ;
  pushText (s, "Caenorhabditis elegans") ;

  spec['a'] = stackMark (s) ;
  pushText (s, "Arabidopsis thaliana") ;

  spec['h'] = stackMark (s) ;
  pushText (s, "Homo sapiens") ;

  spec['m'] = stackMark (s) ;
  pushText (s, "Mus musculus") ;

  spec['r'] = stackMark (s) ;
  pushText (s, "Rattus norvegicus") ;

  spec['d'] = stackMark (s) ;
  pushText (s, "Drosophila melanogaster") ;

  spec['t'] = stackMark (s) ;
  pushText (s, "Tardigrade") ;

  tag['w'] = stackMark (s) ;
  pushText (s, isN ? "AceKogNWorm" : "AceKogWorm") ;

  tag['a'] = stackMark (s) ;
  pushText (s, isN ? "AceKogNAra" : "AceKogAra") ;

  tag['t'] = stackMark (s) ;
  pushText (s, isN ? "AceKogNTardigrade" : "AceKogTardigrade") ;

  tag['h'] = stackMark (s) ;
  pushText (s, isN ? "AceKogNHuman" : "AceKogHuman") ;

  tag['m'] = stackMark (s) ;
  pushText (s, isN ? "AceKogNMouse" : "AceKogMouse") ;

  tag['r'] = stackMark (s) ;
  pushText (s, isN ? "AceKogNRat" : "AceKogRat") ;

  tag['d'] = stackMark (s) ;
  pushText (s, isN ? "AceKogNDroso" : "AceKogDroso") ;
 lao: 

  bestScore = 0 ;
  cp = getPrefix (fi, "<b>Query=</b>") ;
  if (!cp)
    {
      stackDestroy (s) ;
      return nKantor ;
    }

  if (!cp || !*cp++ || !*cp)
    goto lao ;
  nKantor++ ;
  if (isN)
    {
      printf ("\nmRNA \"%s\"\n", cp+5) ;
      printf("AceKogN_Date %s\n", timeShowNow ());
      if (0) printf ("-D AKGN\n") ;/* no because we merge several species */
    }
  else
    {
      printf ("\nKantor \"%s\"\n", cp) ;
      printf("AceKog_Date %s\n", timeShowNow ());
      printf("Kantor_Date %s\n", timeShowNow ());
      if (0) printf ("-D AKG\n") ;/* no because we merge several species */
    }
  nHit = 0 ;
  oldGene = -1 ;

 nextHit:
  cp = getPrefix (fi, "><a name =") ;
  if (!cp) goto lao ;

 /* one letter species code */
  cp = strstr (cp, "</a>") ;
  if (!cp) goto lao ;
  cp += 4 ;
  species = *cp++ ; 

 /* first word of title is gene name */
  cq = strstr (cp, "|") ;
  if (!cq) goto lao ;
  *cq++ = 0 ;
  dictAdd (geneDict, cp, &gene) ;
  cp = cq ;

 /* next word of title is locuslink/newname name */
  cq = strstr (cp, "|") ;
  if (!cq) goto lao ;
  *cq++ = 0 ;
  dictAdd (geneDict, cp, &link) ;
  cp = cq ;

  /* next word of title is product name */
  cq = strstr (cp, "|") ;
  if (!cq) goto lao ;
  if (!strncmp (cp-3,".pg|",4)) goto lao ;
  *cq++ = 0 ;
  product = stackMark (s) ; 
  pushText (s, aceInProtect (fi, cp)) ;
  cp = cq ;

  /* rest if full title on several lines */
  subjectTitle = stackMark (s) ; 
  pushText (s, aceInUnprotect (fi, cp)) ;
  if (aceInCard (fi) &&
      (cp = aceInPos (fi)) &&
      (cp = strstr (cp, "SET\\/")))
    {
      cp += 4 ; *cp = ' ' ;
      cq = strstr (cp, "query") ;
      if (cq ) *cq = 0 ;
      catText (s, cp) ;
    }

  /* Score =  214 bits (545), Expect = 1e-55 */
  cp = getPrefix (fi, "Score") ;
  if (! aceInWord (fi)  || /* Score */
      ! aceInWord (fi)  || /* = */
      ! aceInFloat (fi, &bitScore)  ||  /* 214 */
      ! aceInWord (fi) ||
      ! aceInWord (fi) ||
      ! aceInWord (fi) ||
      ! aceInWord (fi) ||
      ! (cp = aceInWord (fi)))
    {
      fprintf (stderr, "// bad score line %d: %s\n", aceInStreamLine (fi), cp) ;
      goto lao ;
    }
  eValue = stackMark (s) ;
  pushText (s, cp) ; /* very often out of range of a float */
  expect = 1 ;
  if (strstr (cp, "e-")) 
    expect = 0 ;
  else
    sscanf (cp, "%lf", &expect) ;
  aceInCard (fi) ;  aceInCard (fi) ;  aceInCard (fi) ;
  b1 = y1 = -1 ;
  while (1) 
    {
      /* Query: 135 IESSRDLLHRIKDEVGAPGIVVGVSVDGKEVWSEGLGYADVENRVPCKPETVMRIASISK 194 */
      if (!aceInCard (fi)) break ;
      cp = aceInWord (fi) ;
      if (!cp)
	continue ;
      if (strcmp (cp, "Query:"))
	break ;
      if (!aceInInt (fi, &a1) || !aceInWord (fi) || !aceInInt (fi, &a2))
	break ;
      if (b1 == -1) b1 = a1 ;
      /* I+ +++L+       G PG+ + VS+DGK VW  G GYA++E+   C  ++VMRIASISK */
      aceInCard (fi) ;  /* jump that line */
      
      /* Sbjct: 115 IKKAKELVETTMAIQGIPGLSIAVSLDGKMVWKSGFGYANLESFAKCTGDSVMRIASISK 294 */
      if (!aceInCard (fi)) break ;
      cp = aceInWord (fi) ;
      if (!cp || strcmp (cp, "Sbjct:"))
	break ;
      if (!aceInInt (fi, &x1) || !aceInWord (fi) || !aceInInt (fi, &x2))
	break ;
      if (y1 == -1) y1 = x1 ; 
      if (!aceInCard (fi)) break ;
    }

  if (gene &&  gene != oldGene)
    if (y1 != -1 && 100 * bitScore > 88 * bestScore 
	&& (!isN || (isN && expect < .001))      
	)      /* success export one hit */
      {
	oldGene = gene ;
	printf ("%s %s %s"
		, stackText (s, tag[(int)species])
		, dictName (geneDict, gene)
		, dictName (geneDict, link)
		) ;

	if (titles && link < arrayMax (titles))
	  {
	    char *cp = 0 ;
	    
	    cp = arr (titles, link, char *) ;
	    if (cp && *cp)
	      printf (" \"%s\"",  cp) ;
	  }
	printf ("\n") ;
	
	printf ("AceKog%s %s \"%s\" %g %d %d %d %d %s\n"
		, isN ? "N" : ""
		, stackText (s, product)
		, stackText (s, isN ? tag[(int)species] : spec[(int)species])
		, bitScore, b1, a2, y1, x2
		, aceInProtect (fi, messprintf ("%s [%s] eValue = %s",
						stackText (s, subjectTitle) 
						, stackText (s, spec[(int)species])
						, stackText (s, eValue)
						)
				)
		) ;
	if (! nHit || bitScore > bestScore)
	  bestScore = bitScore ;
	nHit++ ;
      }
  if (nHit >= 36)
    goto lao ;
  goto nextHit ;
} /* aceKog2aceParse */

/******************************************************************/

static void usage (void)
{
  fprintf (stderr, "// Usage acekog2ace < blast.out >! blast.ace \n") ;
  exit (1) ;
} /* usage */

/******************************************************************/

static int aceKog2aceGetTitles (const char *titleFileName, DICT *geneDict, Array titles, AC_HANDLE h0)
{
  AC_HANDLE h = ac_new_handle () ;
  ACEIN ai = aceInCreate (titleFileName, 0, h) ;
  char cutter ;
  int gene ;
  int nn = 0 ;

  if (ai)
    {
      aceInSpecial (ai, "\n") ;
      while (aceInCard (ai))
	{
	  char *cp = aceInWordCut (ai, "\t", &cutter) ; /* gene */
	  if (! cp || *cp == '#' || cutter != '\t')
	    continue ;
	  cp = aceInWordCut (ai, "\t", &cutter) ; /* geneNewName of LocusLinkId */
	  if (! cp || *cp == '#' || cutter != '\t')
	    continue ;
	  if (dictAdd (geneDict, cp, &gene))
	    nn++ ;
	  cp = aceInWordCut (ai, "\t", &cutter) ; /* title */
	  array (titles, gene, char *) = strnew (cp, h0) ;
	}
    }
  ac_free (h) ;
  fprintf (stderr, "acekog2ace parsed %d titles\n", nn) ;
  return nn ;
} /* aceKog2aceGetTitles */

/******************************************************************/
/* split the stdin file on Query= */
static int aceKog2aceSplit (ACEIN fi, BOOL isN,  const char *titleFileName)
{
  AC_HANDLE h = ac_new_handle () ;
  Array titles =  0 ;
  Stack s = stackHandleCreate (10000, h) ;
  ACEIN splitFi = 0 ;
  int nn = 0, line = 0 ;
  char *cp, *prefix = "<b>Query=" ;
  DICT *geneDict = dictHandleCreate (50000, h) ;

  if (titleFileName)
    {
      titles = arrayHandleCreate (20000, char *, h) ;
      aceKog2aceGetTitles (titleFileName, geneDict, titles, h) ;
    }
  while (TRUE) 
    {
      if (aceInCard (fi))
	cp = aceInPos (fi) ;
      else
	cp = 0 ;
      line++ ;
      if (!cp || !strncmp (cp, prefix, strlen(prefix)))
	{
	  if (stackMark (s))
	    {
	      splitFi = aceInCreateFromText (stackText (s, 0), 0, 0) ;
	      nn += aceKog2aceParse (splitFi, isN, geneDict, titles) ;
	      aceInDestroy (splitFi) ;
	    }
	  stackClear (s) ;
	}
      if (cp)
	{ catText (s, cp) ; catText (s, "\n") ; }
      else
	break ;
    }

  ac_free (h) ;
  return nn ;
}

/******************************************************************/

int main (int argc, const char **argv)
{
  const char *titleFileName = 0 ;

  getCmdLineOption (&argc, argv, "-t", &titleFileName) ;
  if (argc == 1 || (argc == 2 && !strcmp(argv[1],"-n")))
    {
      ACEIN fi = aceInCreateFromStdin (0, 0, 0) ;
      int nn = aceKog2aceSplit (fi,argc==2 ? TRUE : FALSE, titleFileName) ;
      fprintf (stderr, "// Parse %d blast outputs\n", nn) ;
      aceInDestroy (fi) ;
      printf ("\n") ;
    }
  else
    usage () ;
    
  return 0 ;
}

/******************************************************************/
/******************************************************************/
