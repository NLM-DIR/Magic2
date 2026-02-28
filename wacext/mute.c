/*
  Version 1 Yann Thierry-Mieg 990803
modified but not tested by Jean, sept 2 2003, to use AceC
*/

#include "../wac/ac.h"
#include "array.h"
#include "regular.h"
#include "dna.h"
#include "session.h"
#include "peptide.h"
#include "../wacext/enzlist.h"

static void usage (void) ;
static Array getResuArray ( char * seq , int sizeseq, BOOL debug)  ;
 
/***************************************************************/
/***************************************************************/
/* ACTION */

static void go (Array dna, BOOL ismotif , int pos , char *aatarget , char *motif ,char *tmotif, char *mutseq, BOOL debug)
{
  Array  enzymes=enzInit(), prot=0 , mutdna=0 , resuArray=0 ;
  char targetBuffer[4], *targetprot=0, *startseq=0, *resuseq=0, *resuseq0 = 0, *resuseqn=0, *cp=0, *cq=0 , *cp1=0 , *cq1=0 , enztoaddl [25], enztoaddr [25] , buffer[64] ;
  int i , iEnz, j , iResu, i8, i9, imin, imax, pos3 ;
  int AAtestlength, bestDiffer ;
  int premiere_fois=0, difposl =0 , difposr =0 , differ=0 ,enzoverlaps=0 ,enzoverlape=0  ;
  char c1, c2 ;
  


  ENZ *enz ;
  FILE *ff = filopen ( mutseq , 0 , "w" );
 
  if (!ff) messcrash ("no file") ;
  prot = arrayCreate(arrayMax(dna)/3 +1 ,char);
  
  for (i=0, j=0; i < arrayMax(dna) - 2 ; i += 3)
    array (prot,j++,char)= e_codon(arrp(dna, i, char), 0) ; /* A_ T_ G_ C_ format to peptide */
  
  array (prot,j,char)= 0 ; arrayMax(prot) = j ; /* ensure terminal zero */
  
  /*find pos and target from motif given if  $ismotif=TRUE*/
  if (ismotif) {
    int i1,j1;
    int gagne=0 , posmotif=0 ;
  
    for ( cp = motif ; *cp ; cp++) 
      {
	if ((pepEncodeChar[(int)*cp] < 0) || ( *cp == '*' )) {
	  fprintf (stdout,"Initial motif is not a correct peptide.<br/>Acceptable codes are A=Alanine, C=Cysteine, D=Aspartate, E=Glutamate, F=Phenylalanine, G=Glycine, H=Histidine, I=Isoleucine, K=Lysine, L=Leucine , M=Methionine, N=Asparagine, P=Proline , Q=Glutamine, R=Arginine, S=Serine, T=Threonine, V=Valine, W=Tryptophan, X=Tyrosine .<br/><hr></hr>");
	  exit(-1);
	}
      }
    
    for ( posmotif = 0 , cp = motif , cq = tmotif ; *cp == *cq && *cp ; posmotif++ , cp++ , cq++ );
    if ( !*cp)
      {
	printf ("<br/><hr></hr>Initial and target motif are identical<hr></hr><br/>");
	exit (-1);
      }
    else
      {	
	targetBuffer[0]=*cp;
	targetBuffer[1]='/';
	targetBuffer[2]=*cq ;
	targetBuffer[3]=0 ;
	aatarget = targetBuffer ;
      }

    
    for (i1=0, gagne=0; gagne < 2 && i1 < arrayMax(prot) - strlen(motif); i1++) 
      {	
	cp = arrp(prot,i1,char); 
	cq = motif;
	j1 = strlen(motif) ;
	while (j1--)
	  if (*cp++ != *cq++)
	    break ;
	if (j1  == -1)
	  { gagne++ ; pos = i1 + posmotif + 1 ;}
      }

    switch (gagne)
      {
      case 0:
	printf ("<hr></hr>Motif was not found in proteic sequence. Your dna sequence may not be in the right reading frame.\n<br/><hr></hr>");
	exit(-1);
	break ;
      case 1:
	printf ("Found motif %s at position %d ;mutation target %s <p>",motif,pos,aatarget);
	break ;
      default:
	printf ("<br/><hr></hr>Motif %s was present more than once in sequence. Give a longer motif please.<hr></hr><br/>",motif);
	exit (-1);
	break ;
      }
  }







  
  /*testing ARGV for consistency*/ 
  c1 = *aatarget ;
  c2 = array (prot, (pos-1), char) ;
  if (c1 != c2)
    {
     /* usage ();*/
      fprintf (stdout,"Amino Acid %c not at position %d.\nThe AA"
	       "present in this position is %s\n",
	       *aatarget,pos,pepName[(int)c2]);
      exit(-1);
    };
  
  
  if (pepEncodeChar[(int)*(aatarget+0)] < 0 ||   /* should be of the form A/R */
      *(aatarget+1) != '/' ||
      pepEncodeChar[(int)*(aatarget+2)] < 0 ||
      *(aatarget+3))
    {   
      fprintf (stdout,"<hr></hr>Your target motif is not a correct Peptide."
	       "It must have the same length as initial peptide,"
	       "and one mismatch with same.<br/>Acceptable codes are " 
	       "A=Alanine, C=Cysteine, D=Aspartate, E=Glutamate, F=Phenylalanine, G=Glycine, H=Histidine, I=Isoleucine, K=Lysine, L=Leucine , M=Methionine, N=Asparagine, P=Proline , Q=Glutamine, R=Arginine, S=Serine, T=Threonine, V=Valine, W=Tryptophan, X=Tyrosine .<br/><hr></hr>");
      exit(-1);
    };
  
  /****************************/
  /* run  -f /var/tmp/mutecgi.dna.28921 -motif ALMQ -tmotif ALRQ -mutseq /var/tmp/toto1224 */

  targetprot = messalloc( 256) ;
  startseq = messalloc ( 256) ;
  resuseq0 = messalloc( 256) ;
 

  for (iEnz=0 ;iEnz < arrayMax(enzymes); iEnz++) {
    enz = arrayp(enzymes, iEnz, ENZ) ;
    AAtestlength=(enz->length + 2)/3;  /* + 2 pas plus 1 */
    

    if ( (pos-1-AAtestlength < 0) || (pos + AAtestlength > arrayMax (prot)) )  
	break ; /*no test if outside range*/ 
    
    /*load target proteic sequence and initial dna sequence*/
    memset (startseq, 0, 256) ;
    memset (targetprot, 0, 256) ;

    memcpy (targetprot , arrp(prot,pos-1-AAtestlength,char) , 2*AAtestlength + 1) ;
    *(targetprot + AAtestlength) = *(aatarget+2) ;
    memcpy (startseq, arrp(dna, 3*(pos-1-AAtestlength),char), 3*(2*AAtestlength + 1)) ;
    

   /* printf ("\ntarget prot: %s,AAtestl= %d,\t",targetprot,AAtestlength);*/
    /******************************/
    
    
    /*printf ("\ntesting enzyme %s for fit",enz->name);*/
    for (j= 3*AAtestlength + 1 - (enz->length); j < 3*AAtestlength+3 ; j++) 
      {  
	memset (resuseq0, 0, 256) ;
	memcpy (resuseq0, startseq,  6*AAtestlength+3) ;
	memcpy (resuseq0 + j, enz->dna, enz->length) ;
	resuseq0[6*AAtestlength+3] = 0; 
	if (j + enz->length >= 6*AAtestlength+3)
	  messcrash ("bad length") ;

	/*    Taking Care of enzymes with ambiguous bases 
	 *  left and right of target mutation : take bases 3 by 3 
	 *  and test for compatibility with
	 * starting sequence; if compatible replace with startseq bases
	 */

	for ( cp = startseq , cq = resuseq0 , i=0 ; 
	      *cp && i < 3*AAtestlength  ; i+=3, cp+=3 , cq +=3 ) 
	  if ( ( *cp & *cq ) && 
	       (*(cp + 1) & *(cq + 1 )) && 
	       ( *(cp + 2)  & *( cq + 2 ) )
	       ) {
	    *cq = *cp ;
	    *(cq + 1) = *(cp + 1) ;
	    *(cq + 2) = *(cp + 2) ;
	  }
	    
	    
	for ( cp = startseq + 3*AAtestlength + 3 , 
		cq = resuseq0  + 3*AAtestlength + 3 , i=0 ; 
	      *cp && i < 3*AAtestlength  ; i+=3, cp+=3 , cq +=3 ) 
	  if ( ( *cp & *cq ) && 
	       (*(cp + 1) & *(cq + 1 )) && 
	       (*(cp + 2)  & *( cq + 2 ) )
	       ) {
	    *cq = *cp ;
	    *(cq + 1) = *(cp + 1) ;
	    *(cq + 2) = *(cp + 2) ;
	  }
	
	/********** Most of the n have been taken care of ; 
	 * Merry go round on ambiguities left 
	 ********/

	if (debug) printf ("Enzyme %s :\n <br/>", enz->name );
	resuArray = getResuArray ( resuseq0 , 6 * AAtestlength + 3, FALSE );
	if (debug)
	  printf ("max(resuArray = %d \n<br/>", arrayMax(resuArray)) ;
	
	bestDiffer = 100000 ;
	for (iResu = 0 ; iResu < arrayMax ( resuArray ) ; iResu++) {
	  resuseq = arr ( resuArray , iResu , char * );

	  
	  cp = targetprot ;
	  cq = resuseq ;
	  
	  for (i = 0 ; i < 2*AAtestlength + 1 ; i++, cp++, cq +=3) 
	    if (e_codon(cq, 0) != *cp)
	      break;
	
	  if (i == 2*AAtestlength+1) 
	    { 
	      int i1, j1;
	      int gagne = 0 ;
	    
	      for (i1 = 0; gagne < 1 && i1 < arrayMax(dna) - enz->length; i1++) 
		{
		  cp=arrp(dna,i1,char); cq=enz->dna;
		  j1 = enz->length ;
		  while (j1--)
		    if (*cp++ != *cq++) 
		      break ;
		  if (j1 == -1)
		    gagne++ ; 
		}
	      if (gagne < 1) {
		if (premiere_fois==0) {
		  printf ("\n<p><B>Initial peptide:</B> ");
		  for (i = 3*(pos-4) ; i < 3*(pos+4); i += 3) 
		    printf ("%c", e_codon(arrp(dna, i, char), 0));
		  printf ("  translated from  :");
		  for (i = 3*(pos-4) ; i < 3*(pos+4); i ++) 
		    printf ("%c", dnaDecodeChar[(int)arr(dna, i, char)]);
		  printf ("\n<BR/><B>Mutated peptide:</B> %s\t\n<p>",targetprot);
		  printf ("The enzymes listed will only cut on the mutated sequence,"
			  "they will not cut the original sequence.\n<p>");
		  printf ("<TABLE><TR><TH>Enzyme</TH><TH>Site</TH><TH>"
			  "Resulting sequence</TH><TH>Codes for</TH><TH>"
			  "Extra Bp for cut </TH><TH>Oligos must overlap"
			  "</TH><TH>Mismatches from to</TH><TH>Nb of bp changes"
			  "</TH><TH>Seq added to oligos</TH></TR>\n");
		  
		  
		  
		  premiere_fois++;
		}
		
		
		for ( cp1 = startseq , cq1 = resuseq , difposl = 0 ;
		      *cp1 == *cq1 && difposl <  6*AAtestlength+3 ; 
		      cp1++ , cq1++ , difposl++) ;
		
		for ( cp1 = startseq + 6*AAtestlength+2 , cq1 = resuseq + 6*AAtestlength+2 , 
			difposr = 6*AAtestlength+2 ; 
		      *cp1 == *cq1 && difposr > 0 ; 
		      cp1-- , cq1-- , difposr-- ) ;
		
		for ( cp1 = startseq , cq1 = resuseq , differ = 0 , i = 0;
		      i < 6*AAtestlength+3 ; 
		      cp1++ , cq1++ , i++)  if ( *cp1 != *cq1 ) differ++;


		if (differ < bestDiffer)
		  {
		    bestDiffer = differ ;
		    
		    /* zone de travail */
		    pos3 = 3 * pos ;
		    imin = pos3 - 20 ;
		    if (imin < 0) imin = 0 ;
		    imax = pos3 + 20 ;
		    if (imax > arrayMax(dna))
		      imax = arrayMax(dna) ;
		    
		    /* copy the base dna */
		    for (i8 = imin, cp1 = buffer ; i8 < imax ; i8++, cp1++)
		      *cp1 = arr (dna, i8, char) ;
		    *cp1 = 0 ;
		    
		    /* MUTATE */
		    i8 = 3 * (pos-1-AAtestlength) - imin ;
		    i9 = 3 * (2*AAtestlength + 1) ; /* longueur */
		    strncpy (buffer + i8, resuseq, i9) ;
		    
		    /* copie les bouts utiles */
		    enzoverlaps = 3*(pos-1-AAtestlength) + j - enz->overhang;
		    enzoverlape = 3*(pos-1-AAtestlength) + j + enz->length + enz->overhang;
		    
		    for (i8 = enzoverlaps - imin , i9 = 0 ; i9 < difposr + 1  - j + enz->overhang ; i9++)
		      if (buffer[i8+i9] == arr(dna, imin + i8 + i9, char))
			enztoaddl[i9] =  ace_upper(dnaDecodeChar[(int)buffer[i8+i9]]) ;
		      else
			enztoaddl[i9] =  ace_lower(dnaDecodeChar[(int)buffer[i8+i9]]) ;
		    enztoaddl [i9] = 0 ;
		    
		    for (i8 = enzoverlape - imin - 1, i9 = 0 ; 
			 i9 < j + enz->length + enz->overhang - difposl ; i9++)  
		      if (buffer[i8-i9] == arr(dna, imin + i8 - i9, char))
			enztoaddr[i9] =  ace_upper(dnaDecodeChar[(int)complementBase(buffer[i8 - i9])]) ;
		      else
			enztoaddr[i9] =  ace_lower(dnaDecodeChar[(int)complementBase(buffer[i8 - i9])]) ;
		    enztoaddr [i9] = 0 ;
		    
		    
		    /*
		      printf ("\n<b><p>YO!!</b>Enzyme %s,<br/>of sequence %s,fits"
		      "target mutation and is absent from sequence\n<br/>", 
		      enz->name,enz->seq);
		      printf ("Found match at position :<br/> \n%s ==>  %s <br/>",
		      dnaDecodeString(resuseq),targetprot);
		      */
		    
		    
		    /*
		      resuseqn = messalloc(256) ;
		      strcpy(resuseqn, resuseq) ;
		      */
		    resuseqn = strnew (resuseq,0); 
		    for ( cp1 = startseq , cq1 = resuseqn ; *cp1 ; cp1++ ,cq1++ ) 
		      if (*cp1 != *cq1 ) *cq1 = N_ ;
		    
		    mutdna = dnaCopy (dna);
		    if (arrayMax(mutdna) < 3*(pos-1-AAtestlength) +  6*AAtestlength+3)
		      messcrash ("mutdna is too short") ;
		    memcpy ( arrp ( mutdna  ,3*(pos-1-AAtestlength) , char ) , resuseqn , 6*AAtestlength+3 );
		    if ( ff ) dnaDumpFastA (mutdna ,0 , arrayMax ( mutdna ) , enz->name , ff , 0 );
		    fprintf ( ff ,"\n");
		    
		    messfree (resuseqn) ;
		    arrayDestroy (mutdna) ;

		    printf ("<TR><TD>%s</TD><TD>%s</TD><TD>%s</TD><TD>%s</TD><TD>%d</TD><TD>%d-%d</TD><TD>%d-%d </TD><TD> %d</TD><TD>%s==%s</TD></TR>\n",
			    enz->name , enz->seq , dnaDecodeString(resuseq) , targetprot , enz->overhang ,
			    enzoverlaps ,enzoverlape  ,
			    3*(pos-1-AAtestlength) + difposl -1 , 3*(pos-1-AAtestlength) + difposr +1 , differ,
			    enztoaddl, enztoaddr);
		    /* int enzresu = iEnz; */
		  } /* bestDiffer */
	      }
	      else 
		break ;/* fin de la boucle resuArray ; 
			  printf("<br/>\nEnzyme %s fits target mutation but is present\n\n<br/>",enz->name);*/
	    }
	} /*******end resuArray *****/
      }
  }
  if ( premiere_fois != 0 ) printf ("</TABLE><P>\n");
  /* printf ("Enzyme selected: %s",arrayp(enzymes, enzresu , ENZ)->name) ;*/
  
  
  
  /*  printf ("\n\n target: %s \n",pepName[(int)*(aatarget+2)]);*/
  
  filclose (ff);
  return ;
}  /* go */


/***************************************************************/
/*          subfunctions                                       */
/***************************************************************/

static Array getResuArray ( char * seq , int sizeseq, BOOL debug) 
{
  Array rr= arrayCreate (64, char *);
  char *cp ,  *buff , *buff2 , cc;
  int  i=0 , j =0 ,i1 =0 , nrr = 0, nr1 = 0 , nr2 = 0 ;

  array (rr, nrr++, char*) = messalloc(sizeseq+1) ;
  nr1 = 0 ; nr2 = 1 ;

  for ( cp = seq, j = 0 ; * cp && j <= sizeseq ; cp ++, j++ ) 
    {
      
      
      switch ( *cp ) 
	{
	case A_:
	case T_:
	case G_:
	case C_: 
	case 0:
	  for (i = nr1 ; i < nr2 ; i++)
	    { 
	      buff = arr (rr, i, char*) ;
	      buff[j] = *cp ;
	    }
	  break;
	default:
	  /* modifie les valeurs nr1 nr2 ;*/
	  for (i1 = 0 ; i1 < 4 ; i1++)
	    {
	      
	      cc = 0x1 << i1 ;
	      if (cc & *cp) /* base acceptable */
		{
		  for (i = nr1 ; i < nr2 ; i++)
		    { 
		      buff = arr (rr, i, char*) ;
		      buff2 = array (rr, nrr++, char*) = messalloc (sizeseq+1) ;
		      memcpy (buff2, buff, j) ;
		      buff2[j] = cc ;
		    }
		}
	    }

	  nr1 = nr2 ;
	  nr2 = nrr ;

	  break ;
	}
    }
  
  
  
  
  
  for (j = 0, i = nr1 ; i < nr2 ; j++, i++)
    {
      if (i != j)
	{
	  messfree (arr(rr, j, char*)) ;
	  arr(rr, j, char*) = arr (rr, i, char*) ;
	  arr (rr, i, char*) = 0 ;
	}
    }
  for (i = j ; i < arrayMax(rr) ; i++)
    messfree (arr(rr, i, char*)) ;
  
  arrayMax(rr) = j ;

  if (debug)
    for (i = 0 ; i < arrayMax(rr) ; i++)
      {
	
	printf ("getResArray i = %3d %s :: ", i, 
		dnaDecodeString(seq)) ;
	printf ("%s \n<br/>", 
		dnaDecodeString(arr(rr, i, char*))) ;
      }

  return rr ;
}

/***************************************************************/
/*  utilitaires */
/***************************************************************/

static Array getDnaFromAcedb (char *seqName, char *acedb) 
{
  AC_DB db ;
  AC_OBJ seq ;
  int j ;
  Array dna = 0 ;
  char *dnaTxt = 0 ;
  const char *error = 0 ;

  if (!acedb)
    {
      fprintf (stdout, "Database location MUST be specified\n") ;
      usage() ;
      exit (-1) ;
    }

  db = ac_open_db (acedb, &error) ;
  if (!db)
    exit (1) ;

  seq = seqName ? ac_get_obj (db, "Sequence", seqName, 0) : 0 ;
  if (seq)
    dnaTxt = ac_obj_dna (seq, 0) ;   /*  dna is in  A_ T_ G_ C_ format */

  if (dnaTxt && (j = strlen (dnaTxt)))
    {
      dna = arrayCreate (j, char) ;
      array (dna, j, char) = 0 ; arrayMax(dna) = j ; /* ensure zero terminus */
      memcpy (arrp (dna, 0, char), dnaTxt, j) ;
      dnaEncodeArray (dna) ;  /* go to A_ T_ G_ C_ format */
    }

  aceQuit(FALSE) ;   /* savesession */
  return dna ;
} /* getDnaFromAcedb */

/***************************************************************/

static Array getDnaFromFile (char *filName)
{
  FILE *f = fopen (filName, "r") ;
  Array dna = 0 ;
  int line = 0, j = 0, pos = 0 ;
  char *cp, buffer[65000] ;

  if (!f)
    {
      fprintf (stdout, "<hr></hr>Cannot open file %s\n<HR></HR>", filName) ;
      usage() ;
      exit (-1) ;
    }

  dna = arrayCreate (64000, char) ;

  while (fgets(buffer, 64000, f))
    {
      cp = buffer ;
      if (!line++ && *cp == '>')
	continue ;
      pos = 0 ; /* count char in line */
      while (*cp)
	{
	  switch (*cp)
	    {
	    case 'u':
	    case 'U':
	      *cp = 't' ;
	      /* fall thru */
	    case 'a':
	    case 't':
	    case 'g':
	    case 'c':
	    case 'A':
	    case 'T':
	    case 'G':
	    case 'C':
	      pos++;
	      array (dna, j++, char) = *cp ;
	      break ;
	    case '\n':
	    case 13:
	    case ' ':
	    case '0':
	    case '1':
	    case '2':
	    case '3':
	    case '4':
	    case '5':
	    case '6':
	    case '7':
	    case '8':
	    case '9':
	      pos++;
	      break ;
	    default:
	      {
		fprintf (stdout, "Line %d, position %d, wrong base #%c#(%d)\n", 
			 line-1, pos+1, *cp, *cp) ;
		printf ("<HR></HR>Your pasted sequence is not dna. You must paste a  DNA sequence in FASTA or plain text format<P><HR></HR>");
		exit (-1) ;
	      }
	    }
	  cp++ ;
	}
    }
  array (dna, j, char) = 0 ; arrayMax(dna) = j ; /* ensure zero terminus */
  dnaEncodeArray (dna) ;  /* go to A_ T_ G_ C_ format */
  return dna ;
}


/***************************************************************/
/************************* interface ***************************/
/***************************************************************/

static void usage (void)
{ fprintf (stdout, "Usage: mute [options] $ACEDB\n"
	   "  version 31 august 1999\n"
           "  Directed mutagenesis\n"
	   "\t-f: fasta file\n"
	   "\t-s: sequence_name\n"
	   "\t-acedb: database directory\n"
	   "\t-pos: position of the mutation on proteic sequence\n"
	   "\t\tstarts from 1\n"
	   "\t-target: wild AA/target AA \n"
	   "\t\tie M/N  for Met->Asn\n"
	   "\t-motif: wild motif"
	   "\t-tmotif: mutation desired in motif"
	   ) ;
}

/***************************************************************/

int main(int argc, char **argv)
{ 
  Array dna = 0 ;
  char *seqName = 0, *acedb = 0, *filName = 0 , *mutseq = 0, *aatarget=0, *motif=0, *tmotif=0;
  int pos=0;
  char test[12] ;
  BOOL debug = FALSE ;

  freeinit() ;

  if (debug)
    {
      strcpy (test,"gtmkac") ; 
      dnaEncodeString (test) ;
      getResuArray (test,6,1) ;
    }

 /*debug**** for (i=0;i<argc;i++) printf ("%d: %s  ",i,argv[i]);*/
  for (--argc, ++argv ; argc > 1 ; --argc, ++argv)
    if (!strcmp (*argv, "-s"))
      { if (--argc) seqName = *(++argv) ; }
    else if (!strcmp (*argv, "-f"))
      { if (--argc)  filName  = *++argv ; }
    else if (!strcmp (*argv, "-acedb"))
      {  if (--argc) acedb  = *++argv ; }
    else if (!strcmp (*argv, "-pos"))
      {  if (--argc) pos  = atoi(*++argv) ; }
    else if (!strcmp (*argv, "-target"))
      {  if (--argc) aatarget  = *++argv ; }
    else if (!strcmp (*argv, "-motif"))
      {  if (--argc) motif  = *++argv ; }
    else if (!strcmp (*argv, "-tmotif"))
      {  if (--argc) tmotif  = *++argv ; }
    else if (!strcmp (*argv, "-mutseq"))
      {  if (--argc) mutseq  = *++argv ; }
    else
      { fprintf (stdout, "Unrecognized option %s\n", *argv) ;
        usage() ;
	exit (-1) ;
      }

  if (filName)
    dna =  getDnaFromFile (filName);

  else  if (seqName)
    dna = getDnaFromAcedb (seqName, acedb) ;
  else {printf ("\nno sequence specified\n");usage ();exit(-1);}
  if ( pos && aatarget ) {
    if (pos > 0) {
      go (dna,FALSE, pos, aatarget,0,0,mutseq, debug); 
    } else {
      usage();fprintf(stdout,"Position starts on 1");exit(-1);
    }
  } else if (motif && tmotif) {
    if ( strlen(motif) == strlen (tmotif))
      go (dna,TRUE,0,0,motif,tmotif,mutseq, debug);
    else fprintf (stdout, " <hr></hr>Initial and target motif must have same length.");
  } else {  
    fprintf (stdout, "Not enough command line arguments\n") ;
    usage() ;
    exit (-1) ;
  }
  
  if (0)
    {
      printf ("// date:\t %s\n", timeShowNow()) ;
      printf ("// please send editions and comments to yann@beta.crbm.cnrs-mop.fr\n<p>") ;
    }
  
  
/*  printf("\n// done (if the file is empty, it means that all test succeeded)\n") ;*/
  return 0 ;
}

/***************************************************************/
/***************************************************************/















