#include    <vstd.h>  
#include "../wfiche/biolog.h"
#include    "../wfiche/horribletitle.h"

#define topSentencesForSimilarities 10
#define	minScorePercentageJunkBlast 65
#define minSentencesJunkBlast 0

#define	maxWordsPerSentence 1000 /* NOT a tuning parameter need for the size of sentence arrays */
#define	percentFollowingAllowedToReduceScore	80 /* next word in the sentence is allowed to reduce the score this% much */
#define	percentSimilarityAllowedToReduceScore	80 /* next word in the sentence is allowed to reduce the score this% much */
#define	percentAllowedToDropKernel	0 /* the sentence containing kernel phrases are allowed to have this much of the score of the phrase */
#define	percentDropInHalfanumeric  66 /* if the best is halfanumeric we take the folowing if its score is no less than this much  */

#define	maxFinalSentences	30	/* how many final sentences */
#define	scoreDifferenceOutputPercentage 10 /* sentences with scores this% less than the max are not printed*/

#define	averageFollowScores 0

/* 
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ classification determination functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

static char * specialTreatedSymbols="[]()"; /* special symbols which are not a part of words even if they are written without separating space */
static char * sentenceBreakerSymbols=";";/* breaks the sentences */
static char * breakWordsContain="|"; /* breaks the sentences, also canceles the word it is in */
static char * skipWordsContain="~=:{}"; /* this words are just skipped during the scanning */

static char * skippedUnitedNumbers= /* measurment units which are usually accompanied with numbers */
	"kda "_0
	"kd "_0
	"k "_0
/*	"long "_00 */
	_0;

static char * skippedVariantsAfter=
	"isoform"_0
	"variant"_0
	"in chromosome"_0
	"CG"_0
	_0;

static char * skippedVariantsBothSides=
	"subunit"_0
	_0;

static char * skippedOnlyWithClassifiers=
	"chain"_0
	_0;

static int titleIsClassif(char * src)
{
	int pos,i;
	static char * greekLet="alpha"_0"beta"_0"gamma"_0"delta"_0"epsilon"_00;
	static char * otherClassifs="medium"_0"middle"_0"large"_0"small"_00;
	

	if(isdigit((int)src[0]))	 /* start with number */
		return (!src[1] || isspace((int)src[1])) ? 1 : 2 ;

	if(isalpha((int)src[0]) && (!src[1] || isspace((int)src[1])))	 /* single latin letter */
		return 1; /* is classif */

	if(strchr("ivx",src[0]) )
{ /* roman number */
		for(i=0;src[i] && !isspace((int)src[i]);i++) {
			if(!strchr("ivx",src[i]))break;
		}
		if(!src[i] || isspace((int)src[i]))return 1;
	}

	if((pos=vstrIsOnList(src,greekLet,0,1))!=-1)
{
		if(!src[pos] || isspace((int)src[pos]))return 1; /* is classif */
		return 2;
	}
	if((pos=vstrIsOnList(src,otherClassifs,0,1))!=-1)
{
		if(!src[pos] || isspace((int)src[pos]))return 1; /* is classif */
		return 2;
	}

	return 0;

}


static int titleCheckIfClassification(char * src,int pos,int len)
{
	int p=pos;
	static char * classif="0123456789ivx";

	while(isspace((int)src[p]) && (p-pos)<len-1)p++;
	if(isdigit((int)src[p]))return 1;
	else if(strchr(classif,src[p]))return 1;
	return 0;
}

/* 
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ hit Info Cleanup functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

static char * notTheOnlyFinalWord=
	"protein"_0
	"like protein"_0
	"factor"_0
	"gene"_0
	"homolog"_0 
	"chromosome"_0
	_0;
static char * donotStartTheFinalWord=
	"like"_0
	_0;


static char * spaceSymbols="\t ,";
static char * newlineSymbols="\n";

static char * unrepeateableSymbols="-";
static char * badHitsContain=
	"(xp_"_0 
	"(xm_"_0
	_00;
static char * badHitsStartWith=
        "PREDICTED:"_0
        "similar to "_0
	_00;
static char * skipSentencesContain=
/*	"contains"_0
	"similarity"_0
	"similar to"_0
	"comes from"_0 */
/*	"data source"_0 */
	"expressed sequence"_0 
	"contains ests"_0
	"comes from"_0
	"unnamed"_0
	"supported by genscan"_0
	"sm protein containing"_0
	"come from this gene"_0 
	"predicted using"_0 
	_00;

static char * skipWordsWithNotation =
	"riken cdna"_0
	"cdna"_0
	"bcdna"_0
	"clone"_0
	"conserved hypothetical protein"_0
/*	"protein"_0*/
	_00;

static char * titleFilterPhrases = 
	"caenorhabditis elegans"_0
	"human and rodent"_0
	"human"_0"homo sapiens"_0
	"drosophila melanogaster"_0
	"fruit fly"_0
	"mus musculus"_0" - mouse"_0
	"rat"_0
	"pig roundworm"_0
	" -"_0
	"related to"_0
	"weakly similar to"_0
	"similar to"_0
	"contains similarity to"_0
/*	"contains a"_0			  */
	"putative"_0
	"unknown"_0
	"precursor"_0
	"hypothetical"_0
	"has homology with"_0
	"arabidopsis thaliana"_0
	"fission yeast"_0
/*	" like"_0*/
	" and"_0
/*	" for"_0 */
	" to the"_0 
	" of the"_0 
	" or"_0
	" to"_0 
	"supported by"_0
	_0;

static char * skipSingleWordSentences=
	"orf"_0
	"est"_0
	"protein"_0
	"factor"_0
	"gene"_0
	_00;

static char * skipCertainSentences=
	"gene product"_0
	"expressed protein"_0
	"conserved protein"_0
	"test protein"_0
	_00;

static int titleCountDigAlf(char * src,int * iDig, int * iAlf,int * iDot, int * iDas)
{
	int l;
	if(iDig)*iDig=0;
	if(iAlf)*iAlf=0;
	if(iDot)*iDot=0;
	if(iDas)*iDas=0;
	for(l=0;src[l] && !isspace((int)src[l]);l++)
{
		if(iDig)*iDig+=isdigit((int)src[l]) ? 1 : 0;
		if(iDot)*iDot+=src[l]=='.' ? 1 : 0;
		if(iDas)*iDas+=strchr("-_/",src[l]) ? 1 : 0 ;
		if(iAlf)*iAlf+=isalpha((int)src[l]) ? 1 : 0;
	}
	return l;
}

static int titleCleanSpeciesSpecificJunk(char * src, int len, char * species)
{
  int i,l,k=0,sp,isSpace=1,iDig,iAlf,iDot,iDas,kk,modif;
  /*
    used to be 
    static char *cp, * spcN[]={
    "caenorhabditis elegans"_00,
    "human"_0"homo sapiens"_00,
    "arabidopsis thaliana"_00,
    "drosophila melanogaster"_0"oikopleura dioica"_00,
    "mus musculus"_00
  };
  */
  char *spcN [SPECIESCNT] ;

  for (i=0; i< SPECIESCNT ; i++) /* mieg may 16, 2003, to stay in synch with species enum */
    {
      spcN[i] = SpcI[i].speciesName ;
    } 
  spcN[i] = 0 ;
  /* find the current species */
  for(sp=0;sp<vDim(spcN);sp++)
    {
      if(vstrIsOnList(species,spcN[sp],0,1)!=-1)break;
    }
  /*if(sp==vDim(spcN))return len;*/
  
  for(i=0;i<len;i++)
    {
      
      /* mark beginnings of words */
      if(isspace((int)src[i]))
	isSpace=1;
      else if(isSpace)
	{ 
	  isSpace=0;
	  
	  /* count digits dots and alphabeticals */
	  l=titleCountDigAlf(src+i,&iDig,&iAlf,&iDot,&iDas);
	  
	  switch (sp) {
	  case 0: /* C elegans */
	    if(l>=5 && iDig && iDot && iAlf)i+=l;
	    break;
	  case 1: /* human */
	    {
	      modif=0;
	      if(iDig && ( 
			  !strncasecmp(src+i,"cgi",3) ||
			  !strncasecmp(src+i,"mgc",3) ||
			  !strncasecmp(src+i,"ebi",3) ||
			  !strncasecmp(src+i,"kiaa",4) || 
			  !strncasecmp(src+i,"pro",3) ||
			  !strncasecmp(src+i,"flj",3) ||
			  !strncasecmp(src+i,"xp",2) ||
			  !strncasecmp(src+i,"xm",2) 
			  ))
		{
		  i+=l;
		  modif=1;
		}
	      if(modif)
		{
		  while(isspace((int)src[i]))i++; /* pass the white space and the word "protein" */
		  if(!strncasecmp(src+i,"protein",7))i+=7;
		}
	      break;	
	    }
	    
	  case 2: /* drosophila */
	    modif=0;
	    /*
	      if(iDig && (
	      (src[i]=='c' && src[i+1]=='g' ) || 
	      (src[i]=='r' && src[i+1]=='e' ) ))
	      {
	      i+=l;
	      while(isspace((int)src[i]))i++;
	      if(!strncasecmp(src+i,"gene product",12))i+=12;
	      }*/
	    if(iDig && ( 
			!strncasecmp(src+i,"cg",2) ||
			!strncasecmp(src+i,"re",2) 
			))
	      {
		i+=l;
		modif=1;
	      }
	    if(modif)
	      {
		while(isspace((int)src[i]))i++; /* pass the white space and the word "protein" */
		if(!strncasecmp(src+i,"gene product",12))i+=12;
	      }
	    
	    break;
	  case 3: /* mouse */
	    modif=0;
	    if(iDig && ( 
			!strncasecmp(src+i,"flj",3) ||
			!strncasecmp(src+i,"cg",2) 
			))
	      {
		i+=l;
		modif=1;
	      }else if(!strncasecmp(src+i,"fis",3))
		{
		  i+=3;
		}
	    if(modif)
	      {
		while(isspace((int)src[i]))i++; /* pass the white space and the word "protein" */
		if(!strncasecmp(src+i,"protein",7))i+=7;
	      }
	    
	    if(!strncasecmp(src+i,"mouse",4))
	      i+=l;
	    else if(!strncasecmp(src+i,"data source",11))
	      i+=11;
	    for(kk=0;kk<l;kk++)if(strchr(":~",*(src+i+kk)))break; /* contains bad symbol  ? */
	    if(kk<l)i+=l;
	    
	    break;	
	  default:
	    modif=0;
	    if(iDig && ( 
			!strncasecmp(src+i,"ebi",3) 
			))
	      {
		i+=l;
		modif=1;
	      }
	    if(modif)
	      {
		while(isspace((int)src[i]))i++; /* pass the white space and the word "protein" */
		if(!strncasecmp(src+i,"protein",7))i+=7;
	      }
	    break;
	  }
	  isSpace=1;
	}
      
      if(i!=k)src[k]=src[i];
      k++;
    }
  return k;
}

static char * titlePretreatJunkyBlast(char * hitInfo)
{
  /* pretreat the input */
  char ch,*lastStc,*lastWrd,*ptr,prv='\n',*nxt,cr;
  int  i,k,l,is=0,iDig,iAlf,iDot,iDas,ll,badHit=0,ip;
  int isSpecies=0,iStc;
  char species[vSTRMAXNAME];
  int maxScore=0,curScore=0,prvSpace=1;
  
  
  /* 
     hitInfo 		the buffer to clean up
     hitInfo+i 		current reading position from the buffer 
     hitInfo+k 		current writing position in the buffer 
     ch				the current character to be kept in the buffer
     prv				the last meaningful character read from the buffer
     isSpecies 		flag shows if current position is inside of species brakets
     lastWrd			pointer to the last meaningfull word written to the buffer 
     lastStc 		pointer to the last meaningfull sentence written to the buffer 
     
  */
  
  sscanf( hitInfo,"%d",&maxScore);
  
  /* clean bad words and sentences */
  for(iStc=0,k=0,lastWrd=lastStc=hitInfo,i=0;(ch=hitInfo[i]);i++)
    {
      
      if(prv=='\n')
	{ /* checkfor every new line check if it is NM_ or PREDICTED or similar to etc */
	  char *cp = strstr (hitInfo+i, " ") + 1 ; /* jump the score */

	  badHit=vstrSearchSubstring(hitInfo+i,badHitsContain,1,"\n"_00,1) ? 1 : 0 ;
	  badHit|=(cp && cp == vstrSearchSubstring(cp,badHitsStartWith,1,"\n"_00,1))? 1 : 0 ;
	  
	}
      
      ch=tolower(ch);
      if(strchr("/",ch))
	{hitInfo[i]=ch=' ';}
      
      if(!isSpecies)
	{
	  if(strchr(breakWordsContain,ch) || strchr(sentenceBreakerSymbols,ch))
	    { 
				/*  If we need to cancel the word we just scanned, 
				    reset the write position before the last scanned word.
				    Skip the rest of the sentence by seting the reading 
				    position to the end of the line. However, if there 
				    is an "sp|" hit within the rest of this line, the reading 
				    position should be set to scan it
				    after canceling what was scanned before on this line.
				    Since we might already be here because of '|' 
				    in "sp|" we should substract 3 while checking for 'sp|'.
				    Also scan for the species to the endo f the sentence.
				    To keep the score of the sentence we set the lastWord 
				    and writing position (k) one word after the beginning 
				    of the sentence. And the reading position is set 
				    to skip the "sp|" itself;
				*/
	      if(!strchr(sentenceBreakerSymbols,ch))
		k=lastWrd-hitInfo;
	      for(i++;hitInfo[i] && hitInfo[i]!='\n' ;i++)
		{
		  cr=hitInfo[i];
		  if(!strncasecmp(hitInfo+i-3,"sp|",3))
		    { 
		      
		      lastWrd=vstrSkipWords(lastStc,1,0);
		      k=lastWrd-hitInfo;
		      i=vstrSkipWords(hitInfo+i,1,0)-hitInfo;
		      
		      isSpecies=0;species[0]=0;
		      break;
		    }
		  if(cr=='(' || cr=='[')
		    { isSpecies++;is=0;species[0]=0;}
		  else if(cr==')' || cr==']')
		    {isSpecies--;if(isSpecies<0)isSpecies=0;species[is]=0;}
		  else if(isSpecies)species[is++]=tolower(cr);
		}
	    }
	  else if (!strncmp (hitInfo+i, "eValue=", 7)) 
	    {                   /* skip score */
	      char *cp = hitInfo + i ;
	      l= 0 ;
	      while (*cp != 0 && *cp != '\n') { cp++ ; l++ ; }
	      i+=l;
	      k=lastWrd-hitInfo;
	    }
	  else if(strchr(skipWordsContain,ch))
	    { 
				/* 	Skip the word by setting the writing position before 
					the last word. Move  the reading position to a space 
					or a special character */
	      k=lastWrd-hitInfo;
	      for(i++; hitInfo[i] && !strchr(specialTreatedSymbols,hitInfo[i]) && !isspace((int)(hitInfo[i]));i++);
	    }
	  else if((l=titleCountDigAlf(hitInfo+i,&iDig,&iAlf,&iDot,&iDas))>=8 &&
		  l < 16 && iAlf && iDig 
		  /* && (iDot || !iDas)   unremoved (may 21), removed may 16, 2003 mieg */
		  )
	    {
				/* skip the number+letter identifiers without dash by setting 
				   the reading position after it and the writing position to the 
				   previous word */
	      if (0)
		{
		  char *cp =  hitInfo + i, cc ;
		  cc = *(cp+l) ; *(cp+l)=0 ;
		  printf("Supressing %s\n",cp) ;
		  *(cp+l)=cc ;
		}
	      i+=l;
	      k=lastWrd-hitInfo;
	    }
	  else if((l=vstrIsOnList(hitInfo+i,skippedUnitedNumbers,0,1))!=-1  && prvSpace)
	    {
				/* 	If the current position points to a measurment unit (if the unit name 
					is not accidentaly a part of another word prvSpace==true). 
					Skip it by setting the reading position after the measurment 
					unit and by rolling back the writing position to cancel some 
					spaces and digits wich migh be accompaniing the unit itself. 
					The lastWord position should be set back to account for canceled 
					numbers.
				*/
	      i+=l;
	      for(; k>0 && isspace((int)(hitInfo[k-1]));k--);
	      for(; k>0 && (isdigit((int)(hitInfo[k-1])) || hitInfo[k-1]=='.');k--);
	      lastWrd=hitInfo+k;
	    }
	  
	  else if(	(l=vstrIsOnList(hitInfo+i,skippedVariantsAfter,0,1))!=-1 ||
			(l=vstrIsOnList(hitInfo+i,skippedVariantsBothSides,0,1))!=-1 || 
			(l=vstrIsOnList(hitInfo+i,skippedOnlyWithClassifiers,0,1))!=-1
			)
	    {
				/* 	If the current position points to an unimportand subclassing 
					statement skip it by setting the reading position after it. 
					In case if it is followed by a number or a classifier word,
					skipt it also. In case if the classificator can be before the 
					subclassing word, roll back the writing position to remove 
					it.
				*/
	      
	      ip=i;i+=l;ll=i;
	      for(; hitInfo[i] && isspace((int)(hitInfo[i]));i++); /* skip spaces */
	      if(titleIsClassif(hitInfo+i))
		{
		  for(; hitInfo[i] && !isspace((int)(hitInfo[i]));i++); /* skip spaces */
		}
	      else { 
		if(	(l=vstrIsOnList(hitInfo+ip,skippedVariantsBothSides,0,1))!=-1 || 
			(l=vstrIsOnList(hitInfo+ip,skippedOnlyWithClassifiers,0,1))!=-1
			)
		  {	
		    /* rolling back to skip the scanned spaces and the classifcator */
		    for(l=k; l>0 && isspace((int)(hitInfo[l-1]));l--); 
		    for(; l>0 && !isspace((int)(hitInfo[l-1]));l--);
		    if(titleIsClassif(hitInfo+l))
		      {
			k=l;
			lastWrd=hitInfo+k;
		      }else if((vstrIsOnList(hitInfo+ip,skippedOnlyWithClassifiers,0,1))!=-1)
			{
			  ll=ip;	
			}
		  }
		i=ll;
	      }
	      /*				continue;*/
	    }
	  
	  else if((l=vstrIsOnList(hitInfo+i,skipWordsWithNotation,0,1))!=-1 && prvSpace)
	    {
				/*  Skip word with following notations after checking 
				    if it is not an accidental part of the previous word. */
	      while(isspace((int)hitInfo[i+l]))l++;
	      ll=titleCountDigAlf(hitInfo+i+l,&iDig,0,0,0);
	      if(iDig)i+=l+ll;
	    }
	  
	  if((l=vstrIsOnList(hitInfo+i,titleFilterPhrases,0,1))!=-1 && !isalpha((int)hitInfo[i+l]))
	    {
				/* completely skip some phrases out */
	      i+=l;/* continue; */
	    }
	}
      
      ch=tolower(hitInfo[i]);
      
      /*if((ch=='(' &&  (!i || isspace((int)hitInfo[i-1]) ) ) || ch=='[')
	{ */
      if(ch=='('  || ch=='[')
	{ 
	  /* whatever comes after this ... is a name of a species 
	     actually it is incremental, because some annotators 
	     nest parenthesis in the species names */
	  isSpecies++;is=0;species[0]=0;
	  continue;
	}					 
      /* else if( (ch==')' && (!hitInfo[i+1]  || isspace((int)hitInfo[i+1])) ) || ch==']')
	 {*/
      else if(ch==')'  || ch==']')
	{
	  /* that's it, species name has just finished */
	  isSpecies--;if(isSpecies<0)isSpecies=0;
	  species[is]=0;
	  ch=tolower(hitInfo[i]);
	  continue;
	}
      
      if(isSpecies && !strchr(newlineSymbols,ch))
	{
	  /* the species name is accumulating here */
	  species[is++]=ch;continue;
	}
      
      
      /*if(k!=i)*/hitInfo[k]=ch;
      
      prvSpace=0;
      if(strchr(spaceSymbols,ch))
	{
	  /* 	we remember the position of last accepted words 
		to be able to cancel things after we read 
		and realize how bad they were */
			lastWrd=hitInfo+k; 
			/*			if(isspace((int)prv))k--;
						else */ch=hitInfo[k]=' ';
			prvSpace=1;
	}
      else if(strchr(newlineSymbols,ch))
	{
	  /* finita la comedia:
	     1. 	Check if the last character endline repeats ? 
	     That means that the string already is empty.
	     2. Zero terminate species string just in case it wasn't 
	     finished properly before.
	     3. Skip sentence containing some nasty staff
	     4. Skip hits of C.Elegans to self 
	     5. proceed to cleaning species specific junk 
	     
	  */
	  if(ch==prv)k--; 
	  isSpecies=0;species[is]=0; /* finita la comedia. */
	  
	  ptr=lastStc+((*lastStc=='\n') ? 1 : 0 );
	  
	  /* last cleanup of bad and unhealthy content */
	  if(badHit || vstrSearchSubstring(ptr,skipSentencesContain,1,"\n"_00,1))
	    {
	      k=lastStc-hitInfo;
	    }
	  /*else if(!strcasecmp(species,"caenorhabditis elegans"))
	    {
	    k=lastStc-hitInfo;
	    } */
	  else {
	    l=titleCleanSpeciesSpecificJunk(lastStc,k-(lastStc-hitInfo)+1,species);
	    k=lastStc+l-hitInfo-1;
	  }
	  
	  
	  if(k!=lastStc-hitInfo)
	    { 
				/*  if the last cleanup didn't decide to remove the whole 
				    sentence we do some more filtering. 
				    First we analyze how many alhpabeticals, digits 
				    and other messy symbols are there.
				    Then we carefully skip the score not to discard
				    as a junk, then we discard all negligible score hits. 
				*/
	      
	      l=titleCountDigAlf(ptr,&iDig,&iAlf,&iDot,&iDas);
	      if(iDig && !iAlf)
		{ 
		  sscanf(ptr,"%d",&curScore);if(!maxScore)maxScore=curScore;
		  ptr=vstrSkipWords(ptr,1,0);
		}
	      
	      if(curScore<maxScore*minScorePercentageJunkBlast/100 && iStc>=minSentencesJunkBlast) /* too small score */
		k=lastStc-hitInfo; 
	      else {
		/*  After determining the pointer to the next word 
		    filter if the whole leftover is the score itself.
		    Then in case if it is a single word sentence 
		    check if it's a notation with a dot or its one 
		    of our uneligible single-word sentences. 
		    And then check if the leftover is one of known 
		    ugly criminal letters. */
		
		if(ptr && ptr<hitInfo+k) nxt=vstrSkipWords(ptr,1,0); else nxt=0;
		if(	(!ptr || ptr>=hitInfo+k) ) /* only score ? */
		  k=lastStc-hitInfo; 
		
		else if( (!nxt || nxt>=hitInfo+k) ) {/* single word ? */
		  l=titleCountDigAlf(ptr,&iDig,&iAlf,&iDot,&iDas);
		  if( (iDig && iDot) ||  /* notation name */
		      (vstrIsOnList(ptr,skipSingleWordSentences,0,1)!=-1 ) )
		    { /* bad single word */
		      k=lastStc-hitInfo; 
		    }
		}
		else if( (vstrIsOnList(ptr,skipCertainSentences,0,1)!=-1 ) )
		  { /* bad sentence */
		    k=lastStc-hitInfo;						
		  }
	      }
	    }
	  if(k!=lastStc-hitInfo)iStc++;
	  
	  /* if it passed so much filtering, clean it up from 
	     empty space. Finalize the sentence by \n,
	     remember newly found lastWrd and lastStc.
	     Close the species and keep a note that the last 
	     character wasn't a space  */
	  for(;k>0 && isspace((int)hitInfo[k-1]) ;k--); 
	  hitInfo[k]='\n';
	  lastStc=hitInfo+k;
	  lastWrd=hitInfo+k; 
	  species[0]=0;
	  prvSpace=1;
	}

      
      if(ch!=prv || !strchr(unrepeateableSymbols,ch)) /* get rid of unrepeateable symbols */
	{prv=ch;k++;}
    }
  hitInfo[k]=0; /* finita la comedia */
  
return hitInfo;
}


/* 
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ kernel phrase functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

static int titleJoinIntArraysByAnd(int * dst,int mk, int * p1, int i1, int * p2,int i2)
{
  int i,j,k=0;
  
  for(i=0;i<i1 && k<mk;i++)
    {
      for(j=0;j<i2 && k<mk;j++)
	{
	  if(p1[i]==p2[j])
	    {dst[k++]=p1[i];break;}
	}
    }
  return k;
}

static void titlePrintDic(vSTR * out,vMEG * dic1,int * ind)
{
  int i1,me,maxScore=0;
  
  for(i1=0;i1<dic1->cnt && i1<maxFinalSentences;i1++)
    {
      
      if(ind)me=ind[i1]; else me=i1;
      
      if(i1==0)maxScore=vtextTreeWordScore(dic1,me);
      if(*((char * )vmegId(dic1,me)) && vtextTreeWordScore(dic1,me)>=(scoreDifferenceOutputPercentage*maxScore/100))
	vstrPrintf(out,"%d %d %s\n",vtextTreeWordScore(dic1,me),vtextTreeWordCnt(dic1,me),vmegId(dic1,me));
    }
}


static void titleFollowWord(vMEG * dicStc,vMEG * dicWrd,int startWrd,int * stc,int stcdim,int * wrds,int iPos,vMEG * dicRes, int maxWrds)
{
  vMEX * ref,* substc;
  int i,j,score,cntcrs,crs[maxWordsPerSentence],scorecrs,subAny=0,curWrd=startWrd,iStc,allcnt,allcntcrs;
  
  ref=vtextTreeWordRef(dicWrd,curWrd, vTEXTTREELNKNXT); /* my next links */
  score=0;allcnt=0;
  for(j=0;j<stcdim;j++)
    {
      score+=vtextTreeWordScore(dicStc,stc[j])*vtextTreeWordCnt(dicStc,stc[j]);
      allcnt+=vtextTreeWordCnt(dicStc,stc[j]);
    }
  if(averageFollowScores && allcnt)score/=allcnt;
  
  wrds[iPos++]=startWrd;
  
  for(i=0;i<ref->cnt;i++)
    {
      curWrd=vmexArr( ref, i, int);
      
      if(curWrd==vTEXTTREEENDWRD)break;
      for(j=0;j<iPos && j < maxWrds;j++)if(wrds[j]==curWrd)break;
      if(j<iPos)continue;
      
      substc=vtextTreeWordRef(dicWrd,curWrd, vTEXTTREELNK); /* my sentences */
      
      /* get cross section of all sentences containing both words */
      cntcrs=titleJoinIntArraysByAnd(crs,vDim(crs),stc,stcdim,(int *)(substc->buf),substc->cnt);
      
      scorecrs=0;allcntcrs=0;
      for(j=0;j<cntcrs;j++)
	{
	  scorecrs+=vtextTreeWordScore(dicStc,crs[j])*vtextTreeWordCnt(dicStc,crs[j]);
	  allcntcrs+=vtextTreeWordCnt(dicStc,crs[j]);
	}
      if(averageFollowScores && allcntcrs)scorecrs/=allcntcrs;
      
      /* see if this doesn't decrease the score too much */
      if(iPos==1 || scorecrs>=percentFollowingAllowedToReduceScore*score/100)
	{
	  /* check if this word can be the first */
	  /*			if(!vstrSearchSubstring((char *)vmegId(dicWrd,curWrd),notTheFirstWord,1,_00,1))
				{*/
	  titleFollowWord(dicStc,dicWrd,curWrd,crs,cntcrs,wrds,iPos,dicRes,maxWordsPerSentence);
	  subAny=1;
	  /*			}*/
	}
    }
  
  /* the finale ... remember the current sentence  */
  if(!subAny)
    { 
      vSTR buf;
      vTEXTTREEWRD lnk;
      int	iAdd;
      vMEX * ref;
      
      vSet0(buf);vSet0(lnk);lnk.score=score;lnk.cnt=allcnt;
      
      for(i=1;i<iPos && i < maxWrds;i++) /* compose the destination sentence */
	vstrPrintf(&buf,"%s%s",i>1 ? " ":"" ,(char * )vmegId(dicWrd,wrds[i]));
      
      if(vstrPtr(&buf))
	{
	  /* check if the adding string is not someones substring */
	  for(iStc=0;iStc<dicRes->cnt;iStc++)if(!strcmp((char *)vmegId(dicRes,iStc),vstrPtr(&buf)))break;
	  
	  if(iStc==dicRes->cnt)
	    { /* not a part of the previous sentence */
				/* if(vmegFind(dicStc,0,vstrPtr(&buf)))
				   { }*/
	      vmegAdd(dicRes,&iAdd,vstrPtr(&buf),&lnk,sizeof(lnk)); /* add it to dictionary */
	      ref=vtextTreeWordRef(dicRes,iAdd,vTEXTTREELNK);
	      for(i=1;i<iPos && i < maxWrds;i++) /* add the links to it's words */
		vtextTreeLinkAttach(ref,wrds[i]);
	    }
	}
      
      vstrEmpty(&buf);
    }
}

/* 
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ similarities functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


static void titleJoinSimilarities (vMEG * dicSim,char * st1, char * st2,char * sBuf,int simcnt)
{
	int i,* pij,len,iSim,l;
	
	/* fill with spaces */
	for(i=0;i<simcnt;i++)sBuf[i]=' ';
	sBuf[simcnt-1]=0;

	for(iSim=0;iSim<dicSim->cnt;iSim++)
{
		pij=(int * )vmegId(dicSim,iSim);
		len=*((int *)(vmegData(dicSim,iSim)));
		
		/* filter half word similarities if it starts not at the beginning of words */
		if( ( pij[0]>0 && isalpha((int)st1[pij[0]]) && isalpha((int)st1[pij[0]-1]) ) || 
			( pij[1]>0 && isalpha((int)st2[pij[1]]) && isalpha((int)st2[pij[1]-1])) )
{
			while (st1[pij[0]]!=' ' && len>0)
{ pij[0]++;pij[1]++;len--;} /* adjust the lenght and start */
			*((int *)(vmegData(dicSim,iSim)))=len; 
		}

		/* filter half word similarities if it ends not at the beginning of words */
		if( ( st1[pij[0]+len] && st1[pij[0]+len]!=' ' && st1[pij[0]+len-1]!=' ') && 
			( st2[pij[1]+len] && st2[pij[1]+len]!=' ' && st2[pij[1]+len-1]!=' ') )
{
			while (!isdigit((int)(st1[pij[0]+len-1])) && st1[pij[0]+len-1]!=' ' && len>0)
{ len--;} /* adjust the lenght and start */
			*((int *)(vmegData(dicSim,iSim)))=len; 
		}

		if(st1[pij[0]] && strchr("/+",st1[pij[0]]))
{ /* extend if it was punctuation mark */
			while (pij[0]>0 && st1[pij[0]]!=' ')
{ pij[0]--;pij[1]--;len++;} /* adjust the lenght and start */
		}

		if(!len)continue;

		if(len<4 && !titleCheckIfClassification(st1,pij[0],len) )continue;

		/* the similarity itself */
		for(l=0;l<len;l++)sBuf[pij[0]+l]=st1[pij[0]+l];
	}

}


static int titleFindSimilarities(vMEG* dicStc,int * stc, int stccnt,vMEG * dicRes)
{
	int i,j,iRes,/*maxScore,score,*/cnt1,cnt2,scr1,scr2;
	char * s1,*s2,sBuf[vTEXTMAXLINE]/*,st1[vTEXTMAXLINE],st2[vTEXTMAXLINE]*/;
	vTEXTTREEWRD lnk;
	vMEG dicSim;
		
	vSet0(lnk);
	vSet0(dicSim);
	
	for(i=0;i<stccnt;i++)
{
		s1=(char *)vmegId(dicStc,stc[i]);
		cnt1=vtextTreeWordCnt(dicStc,stc[i]);
		scr1=vtextTreeWordScore(dicStc,stc[i]);

		for(j=i;j<stccnt;j++)
{
			s2=(char *)vmegId(dicStc,stc[j]);
			cnt2=vtextTreeWordCnt(dicStc,stc[j]);
			scr2=vtextTreeWordScore(dicStc,stc[j]);

			/*if(s1==s2 && cnt1==1)
{sBuf[0]=0; continue; }*/
			vmegEmpty(&dicSim);
			vtxtCompareSentences(s1,s2,&dicSim);

			titleJoinSimilarities (&dicSim,s1,s2,sBuf,vDim(sBuf));
			vstrFindReplaceSymbols(sBuf,sBuf,0,vSTRBLANK," ",0,1,1);
			vstrCleanEnds(sBuf,sBuf,0,vSTRBLANK,1);
			
			if(i==j)lnk.cnt=(cnt1-1)*cnt1/2;
			else lnk.cnt=cnt1*cnt2;
			if(!sBuf[0]  || !lnk.cnt)continue;
			lnk.score=lnk.cnt*(scr1+scr2)/2;
			
			if(!vmegAdd(dicRes,&iRes,sBuf,&lnk,sizeof(lnk)))
{ /* add it to dictionary */
				vtextTreeWordScore(dicRes,iRes)+=lnk.score; /* if it was there before, just add it */
				vtextTreeWordCnt(dicRes,iRes)+=lnk.cnt; 
			}

			
		}
	}
	vmegEmpty(&dicSim);
	return dicRes->cnt;
}

/* 
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ sentence refinement functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/

static void titleRefineOrder(vMEG * dicWrd, vMEG * dicStc, char * dst,int * cwrd, int cntwrd,char * caseBuf)
{
  int i,j,iMinNum,l,ii,tt,allCaps,allLets;
  int wrds[maxWordsPerSentence];
  char buf1[vTEXTMAXLINE],buf2[vTEXTMAXLINE],* ptr,*totu;
  
  for(j=0;j<cntwrd;j++)
    {
      sprintf(buf1," %s ",(char * )vmegId(dicStc,cwrd[j]));
      
      for(i=0;i<dicWrd->cnt;i++)
	{
	  sprintf(buf2," %s ",(char * )vmegId(dicWrd,i));
	  
	  if(!(ptr=strstr(buf1,buf2)))break; /* doesn't exist in this sentence */
	  wrds[i]=ptr-buf1;
	}
      if(i==dicWrd->cnt)break; /* found a sentence with all words */
    }
  
  if(j!=cntwrd)
    {
      for(i=0;i<dicWrd->cnt;i++)
	{
	  iMinNum=0;
	  for(j=0;j<dicWrd->cnt;j++)
	    {
	      if(wrds[j]<wrds[iMinNum])iMinNum=j;
	    }
	  ptr=(char * )vmegId(dicWrd,iMinNum);l=strlen(ptr);
	  if(caseBuf)
	    {
	      for(ii=0;caseBuf[ii];ii++)
		{
		  if( !strncasecmp(caseBuf+ii,ptr,l) && 
		      !isalpha((int)caseBuf[ii+l]) && 
		      (l==0 || !isalpha((int)caseBuf[ii-1])) )
		    {
		      
		      for(allLets=0,allCaps=0,tt=0,totu=caseBuf+ii;tt<l;totu++,tt++)
			{
			  if( (caseBuf[ii+tt])>='A' && (caseBuf[ii+tt])<='Z' )allCaps++;
			  if( ((caseBuf[ii+tt])>='A' && (caseBuf[ii+tt])<='Z' )  ||
			      ((caseBuf[ii+tt])>='a' && (caseBuf[ii+tt])<='z' ) )allLets++;
			  
			}
		      
		      if (l < vTEXTMAXLINE - 1)
			{ strncpy(buf1,caseBuf+ii,l);ptr=buf1;buf1[l]=0; }
		      if(l<4 || allCaps!=allLets)
			break;
		    }
		}
	    }
	  dst+=sprintf(dst,"%s%s",i ? " " : "",ptr);
	  wrds[iMinNum]=0x7fffffff;
	}
    }else {
      for(i=0;i<dicWrd->cnt;i++)
	{
	  ptr=(char * )vmegId(dicWrd,i);l=strlen(ptr);
	  if(caseBuf)
	    {
	      for(ii=0;caseBuf[ii];ii++)
		{
		  if( !strncasecmp(caseBuf+ii,ptr,l) && 
		      !isalpha((int)caseBuf[ii+l]) && 
		      (l==0 || !isalpha((int)caseBuf[ii-1])) )
		    {
		      
		      strncpy(buf1,caseBuf+ii,l);ptr=buf1;buf1[l]=0;
		      break;
		    }
		}
	    }
	  dst+=sprintf(dst,"%s%s",i ? " " : "",ptr);
	}
    }
}


/* 
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
_/
_/ major functions
_/
_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
*/


static int titleGetCollectionSentences(char * stc,vMEG * dicWrd,int * cwrd)
{
  int wrds[maxWordsPerSentence],cntwrd=0,i,j,iFnd;
  vMEX * ref;
  vMEG dicRes;
  
  vSet0(dicRes);
  vtextDicSentenceConsume(stc,&dicRes,0,0);
  
  for(i=0;i<dicRes.cnt;i++)
    { /* list of words in current result */
      
      /* find this word in word list */
      if(!vmegFind(dicWrd,&iFnd,(char * )vmegId(&dicRes,i)))continue;
      
      /* list of all sentences of the current word */
      ref=vtextTreeWordRef(dicWrd,iFnd, vTEXTTREELNK); 
      
      if(i==0)
	{ /* copy the list for the first */
	  for(j=0;j<ref->cnt;j++)cwrd[j]=vmexArr(ref,j,int);
	  cntwrd=j;
	  continue;
	}
      else /* join this list with previous */
	cntwrd=titleJoinIntArraysByAnd(wrds,vDim(wrds),(int *)(ref->buf),ref->cnt,cwrd,cntwrd);
      
      /* copy the result to be compared at the next stage */
      for(j=0;j<cntwrd;j++)cwrd[j]=wrds[j];
    }
  
  vmegEmpty(&dicRes);
  return cntwrd;
  
}

static int titleComputeSentenceScore(vMEG* dicStc,vMEG * dicWrd,char * stc,int * pcnt)
{
  int cwrd[maxWordsPerSentence],cntwrd,j,score,allcnt;
  
  if(pcnt)*pcnt=0;
  
  cntwrd=titleGetCollectionSentences(stc,dicWrd,cwrd);
  for(allcnt=0,score=0,j=0;j<cntwrd;j++)
    {
      score+=vtextTreeWordScore(dicStc,cwrd[j])*vtextTreeWordCnt(dicStc,cwrd[j]);
      allcnt+=vtextTreeWordCnt(dicStc,cwrd[j]);
      if(pcnt)*pcnt+=vtextTreeWordCnt(dicStc,cwrd[j]);
    }
  if(0 && allcnt)score/=allcnt;
  
  return score;
}

static int titleThePhraseAndClassfier(char * excl, char * stc)
{
	int len=strlen(excl);
	int iDig, iAlf,iDot, iDas;
	char bfer[vTEXTMAXLINE],*ptr=0;

	strcpy(bfer,stc);
	if(	!(ptr=strstr(bfer,excl)) )return 0; /* skip if this sentence doesn't contain the phrase=res */
	if( (ptr[len] && !isspace((int)ptr[len])) || (ptr!=bfer && !isspace((int)(*(ptr-1)))) )return 0;
	memset(ptr,' ',len); /* clean our excl from the current sentence */

	/* scan until see non-classifier */
	for(ptr=bfer;*ptr && isspace(((int)(*ptr)));ptr++);
	for(;ptr;ptr=vstrSkipWords(ptr,1,0))
{
		titleCountDigAlf(ptr,&iDig,&iAlf,&iDot,&iDas);

		if( !titleIsClassif(ptr) &&		/* exclude classifier words */
			!(iAlf && iDig && !iDas)	/* exclude halfanumeric words */
			)break;	
	}
	if(!ptr)/* all were classfiers */
		return 1;
		
	return 0;
}


static int titleMeaningfulEnoughtToBeAlone(vMEG * dicStc,char * excl)
{
  int il;
  
  if(vmegFind(dicStc,0,excl))return -1; /* it occurs alone */
  for(il=0;il<dicStc->cnt;il++)
    { /* scan all sentences find the top sentence containing this kernel phrase */
      if(titleThePhraseAndClassfier(excl,(char *)vmegId(dicStc,il)))
	return il+1;
    }
  
  return 0;
}

static int titleCheckSubstringAndAdd(vMEG * dicFin,char * excl,vTEXTTREEWRD * lnk, int checkCaseSens)
{
  int i,iAdd,c1=0,c2=0,added,cnt2,cnt1;
  vMEG dicWrd;
  char * ptr=0,bfer[vTEXTMAXLINE];
  
  if(excl && strlen(excl)<4)return 0;
  
  vSet0(dicWrd);
  
  /* check if current sentence is subsentence of one already added */
  if(!vmegFind(dicFin,&iAdd,excl))
    {
      for(i=0;i<dicFin->cnt;i++)
	{
	  strcpy(bfer,(char * )vmegId(dicFin,i));
	  if(!checkCaseSens)vstrLowerCase(bfer);
	  vtextDicSentenceConsume(bfer,&dicWrd,&cnt1,0);
	  c1=dicWrd.cnt;
	  
	  strcpy(bfer,excl);
	  if(!checkCaseSens)vstrLowerCase(bfer);
	  added=vtextDicSentenceConsume(bfer,&dicWrd,&cnt2,0);
	  c2=dicWrd.cnt;
	  vmegEmpty(&dicWrd);
	  if(c2==c1)break;
	  /* if the second sentence contains the first one delete the first one  */
	  if(c1+added==c2 && lnk->score>=vtextTreeWordScore(dicFin,i) && cnt2>cnt1)
				/* vmegDel(dicFin,i,1);  */
	    *((char * )vmegId(dicFin,i))=0;
	}
      if(i<dicFin->cnt)return 0 ; /* part of sentence before */
    }
  
  
  /*
    for(ptr=notTheOnlyFinalWord; ptr && strcmp(excl,ptr);ptr=vstrNextList(ptr,1));
    if(ptr)return 0;
  */
  for(ptr=notTheOnlyFinalWord; 
      ptr && !titleThePhraseAndClassfier(ptr,excl);
      ptr=vstrNextList(ptr,1));
  if(ptr)return 0;
  
  
  for(ptr=donotStartTheFinalWord; ptr && strncmp(excl,ptr,strlen(ptr));ptr=vstrNextList(ptr,1));
  if(ptr)return 0;
  
  
  if(!vmegAdd(dicFin,&iAdd,excl,lnk,sizeof(vTEXTTREEWRD)))
    { /* add it to dictionary */
      /*		vtextTreeWordScore(dicFin,iAdd)+=lnk->score;*/
      vtextTreeWordCnt(dicFin,iAdd)+=lnk->cnt;
    }
  return 1;
}


static char *titleDoBestSuggestionANALYZER (char *hitInfoOriginal)
{
  vSTR out,dbg;
  int dodebug=1;
  static int firstTime=1;
  vMEG dicStc,dicWrd,dicRes,dicFin,dicSim,dicOrd;
  char excl[vTEXTMAXLINE],* res,* sim,* hitInfo,*ptr;
  int iSim,iRes,is,ir,scaleScore,consScore,base,cntwrd,isOK;
  int sortSim[maxWordsPerSentence],sortRes[maxWordsPerSentence],stcWrd[maxWordsPerSentence],wrds[maxWordsPerSentence];
  vTEXTTREEWRD mel;
  
  if(!hitInfoOriginal || !(hitInfo=strnew(hitInfoOriginal, 0)) )return 0;
  
  vSet0(out);
  vSet0(dicStc);vSet0(dicWrd);vSet0(dicRes);vSet0(dicFin);
  vSet0(dicSim);vSet0(dicOrd);
  vSet0(mel);
  
  if(dodebug)
    {
      vSet0(dbg);
      dbg.mode=vSTRDISCARD;
      dbg.outFile="/home/mieg/debug";
      if(firstTime)
	{firstTime=0; filremove(dbg.outFile, 0);}
      vstrPrintf(&dbg,"Original input\n%s\n\n\n",hitInfo);
    }
  /* clean junk from blast hits */
  titlePretreatJunkyBlast (hitInfo);
  if(dodebug)
    {
      vstrPrintf(&dbg,"Pretreated input\n%s\n\n\n",hitInfo);
    }
  
  
  /* build the text graph - tree */
  vstrFindReplaceSymbols(hitInfo,hitInfo,0,vSTRENDLINE,0,0,1,1); /* tokenize by double zeros also remove empty lines */
  vtextTreeCreate (hitInfo,&dicStc,&dicWrd,1,0); 
  vtextDicSortIndex(&dicStc,1,stcWrd,maxWordsPerSentence); /* sort the sentence dictionary */
  if(dodebug)
    {
      if(0)
	{
	  vtextTreePrint(&dbg,"Sentences",&dicStc,&dicWrd,vTEXTTREELNK,0);
	  vtextTreePrint(&dbg,"Words",&dicWrd,&dicStc,vTEXTTREELNK,0); 
	  vtextTreePrint(&dbg,"Next Neighbours",&dicWrd,&dicWrd,vTEXTTREELNKNXT,0);
	  vtextTreePrint(&dbg,"Previous Neighbours",&dicWrd,&dicWrd,vTEXTTREELNKPRV,0);
	}
      vstrPrintf(&dbg,"Dictionarized input\n");
      titlePrintDic(&dbg,&dicStc,stcWrd);
    }
  
  
  /* find the kernel phrases */
  if(dicStc.cnt && dicWrd.cnt)
    {
      vMEX * ref;
      /*		vMEX * ref=vtextTreeWordRef(&dicWrd,vTEXTTREEBEGINWRD, vTEXTTREELNK);
			titleFollowWord(&dicStc,&dicWrd,vTEXTTREEBEGINWRD,(int *)(ref->buf),ref->cnt,wrds,0,&dicRes,maxWordsPerSentence);*/
      for(iRes=0;iRes<dicWrd.cnt;iRes++)
	{ /* for every result */
	  ref=vtextTreeWordRef(&dicWrd,iRes, vTEXTTREELNK);
	  titleFollowWord(&dicStc,&dicWrd,iRes,(int *)(ref->buf),ref->cnt,wrds,0,&dicRes,maxWordsPerSentence);
	}
      vtextDicSortIndex(&dicRes,1,sortRes,maxWordsPerSentence); /* sort the kernel sentences */
    }
  if(dodebug)
    {
      vstrPrintf(&dbg,"\n kernel phrases \n\n");
      for(iRes=0;iRes<maxFinalSentences && iRes<dicRes.cnt;iRes++)
	{ /* for every result */
	  vstrPrintf(&dbg,":::::::: %d %d \"%s\"\n",vtextTreeWordScore(&dicRes,sortRes[iRes]),vtextTreeWordCnt(&dicRes,sortRes[iRes]),(char *)vmegId(&dicRes,sortRes[iRes]));
	}
    }
  
  
  /*  look for similarities between top score sentences. 
      First, count how many top sentences should be accounted 
      taking into account multiple repetitions of hits.
      Then Look for similarities and sort them */
  for(cntwrd=base=0;	cntwrd<dicStc.cnt && base<topSentencesForSimilarities ; cntwrd++)
    base+=vtextTreeWordCnt(&dicStc,stcWrd[cntwrd]);
  titleFindSimilarities(&dicStc,stcWrd,cntwrd,&dicSim);
  if(dicSim.cnt)vtextDicSortIndex(&dicSim,1,sortSim,maxWordsPerSentence);
  if(dodebug)
    {
      vstrPrintf(&dbg,"\n Sentences to look for similarities \n\n");
      for(iSim=0;iSim<cntwrd;iSim++)
	{ 
	  vstrPrintf(&dbg,":::::::: %d %d \"%s\"\n",vtextTreeWordScore(&dicStc,stcWrd[iSim]),vtextTreeWordCnt(&dicStc,stcWrd[iSim]),(char *)vmegId(&dicStc,stcWrd[iSim]));
	}
      vstrPrintf(&dbg,"\n Top Similarities \n\n");
      for(iSim=0;iSim<dicSim.cnt;iSim++)
	{
	  vstrPrintf(&dbg,":::::::: %d %d \"%s\"\n",vtextTreeWordScore(&dicSim,sortSim[iSim]),vtextTreeWordCnt(&dicSim,sortSim[iSim]),(char *)vmegId(&dicSim,sortSim[iSim]));
	}
							}
  
  
	/* try to enrich every kernel sentence */
  if(dodebug)
    {
      vstrPrintf(&dbg,"\n Kernel enrichments \n\n");
    }
  for(iRes=0;iRes<maxFinalSentences && iRes<dicRes.cnt;iRes++)
    { /* for every kernel phrase */
      ir=sortRes[iRes]; /* ir is the sorted index order */
      res=(char *)vmegId(&dicRes,ir);
      if(dodebug)
	{
	  vstrPrintf(&dbg,"==================================\n");
	  vstrPrintf(&dbg,"%d %d \"%s\"\n",vtextTreeWordScore(&dicRes,ir),vtextTreeWordCnt(&dicRes,ir),res);
	}
      isOK=0;
      excl[0]=0;																													 
      
      /* consider each similarity */
      for(scaleScore=0,iSim=0;iSim<dicSim.cnt;iSim++)
	{
	  is=sortSim[iSim];	/* is is the sorted index order */
	  sim=(char *)vmegId(&dicSim,is);

	  /* this similarity is subcontent of the kernel sentece ? 
	     skip it since it doesn't add any more information */
	  if( strstr(res,sim ) )
	    continue;
	  
	  /*  the kernel sentence is not cotained in this similarity ? 
	      Trying to extend it by this similarity can result in 
	      combining uncombinable subphrases. */
	  /* if(	!(ptr=strstr(sim,res)) || (ptr==sim || !isspace((int)(*(ptr-1)))) || !strcmp(sim,res) ) */
	  if(	!(ptr=strstr(sim,res)) || (ptr!=sim && !isspace((int)(*(ptr-1)))) || !strcmp(sim,res) ) 
	    continue;
	  /* if(	!(ptr=strstr(sim,res)) || !strcmp(sim,res) ) continue;*/
	  
	  /* join the words of the kernel phrase and the similarity 
	     in the natural order. The order - occuring in the all 
	     senctences of dictionary. */ 
	  vtextDicSentenceConsume(res,&dicOrd,0,sizeof(int));
	  vtextDicSentenceConsume(sim,&dicOrd,0,sizeof(int));
	  titleRefineOrder(&dicOrd,&dicStc,excl,stcWrd,cntwrd,hitInfoOriginal);
	  vmegEmpty(&dicOrd);
	  
	  {
				char lwrc[vTEXTMAXLINE];
				strcpy(lwrc,excl);vstrLowerCase(lwrc);
				consScore=titleComputeSentenceScore(&dicStc,&dicWrd,lwrc,0);
	  }
	  if(dodebug)
	    {
	      vstrPrintf(&dbg," considering : %5d %3d \"%s\"\n",vtextTreeWordScore(&dicSim,is),vtextTreeWordCnt(&dicSim,is),sim);
	      vstrPrintf(&dbg,"        score drop : from %5d to %5d %.2lf%%\n",vtextTreeWordScore(&dicRes,ir),consScore,consScore*100./vtextTreeWordScore(&dicRes,ir));
							}
	  
	  if(consScore<(percentSimilarityAllowedToReduceScore*vtextTreeWordScore(&dicRes,ir)/100) )
	    {
	      
	      if(dodebug)
		{
		  vstrPrintf(&dbg,"        discarded : the score of the mixture (%5d) is less than %d%% of the kernel phrase score\n",consScore,percentSimilarityAllowedToReduceScore);
		}
	      continue;
	    }
	  
	  
	  /*  compute and rescale its score so the best similarity 
	      updated sentence has the score of the current 
	      kernel phrase */
	  mel.score=vtextTreeWordScore(&dicSim,is);/*titleComputeSentenceScore(&dicStc,&dicWrd,excl,0,topSentencesForSimilarities);*/
	  mel.cnt=vtextTreeWordCnt(&dicSim,is);
	  if(!scaleScore)scaleScore=mel.score; /* the first is the the best similarity update */
	  if(scaleScore)mel.score=mel.score*vtextTreeWordScore(&dicRes,ir)/scaleScore;
	  
	  if(dodebug)
	    {
	      vstrPrintf(&dbg,"        mixing : %5d %3d \"%s\"\n",mel.score,mel.cnt,excl);
	    }
	  
	  /* add to final dictionary */
	  isOK=titleCheckSubstringAndAdd(&dicFin,excl,&mel,0);
	}
      
      
      /* if this kernel phrase didn't go through similarity 
	 selection process for some reason add it after careful 
	 attmepts to extend.
	 If the phrase occurs by itself as a hit, add it as it is.
	 If the Kernel phrase does'nt occur as a hit, add only if 
	 the top hit containing this kernel phrase has score drop 
	 less than a 50%.
      */
      if(!scaleScore)
	{
	  int presentAlone=0;
	  mel.score=vtextTreeWordScore(&dicRes,ir);
	  mel.cnt=vtextTreeWordCnt(&dicRes,ir);
	  
	  if(!(presentAlone=titleMeaningfulEnoughtToBeAlone(&dicStc,res)))
	    {
	      int il;
	      double maxScore=vtextTreeWordScore(&dicStc,0);
	      
	      for(il=0;il<dicSim.cnt;il++)
		{ /* find the top sentence containing this kernel phrase */
		  if((strstr((char *)vmegId(&dicStc,il),res)) && 
		     vtextTreeWordScore(&dicStc,il)>percentAllowedToDropKernel*maxScore/100 )
		    {
		      strcpy(excl,(char *)vmegId(&dicStc,il));
		      
		      vtextDicSentenceConsume(excl,&dicOrd,0,sizeof(int));
		      titleRefineOrder(&dicOrd,&dicStc,excl,stcWrd,cntwrd,hitInfoOriginal);
		      vmegEmpty(&dicOrd);
		      mel.score=titleComputeSentenceScore(&dicStc,&dicWrd,excl,0);
		      isOK=titleCheckSubstringAndAdd(&dicFin,excl,&mel,0);
		      break;
		    }
		}
	    }
	  else {
	    if(presentAlone>0)mel.score=vtextTreeWordScore(&dicStc,(presentAlone-1));
				/* strcpy(excl,res); */
	    vtextDicSentenceConsume(res,&dicOrd,0,sizeof(int));
	    titleRefineOrder(&dicOrd,&dicStc,excl,stcWrd,cntwrd,hitInfoOriginal);
	    vmegEmpty(&dicOrd);
	    isOK=titleCheckSubstringAndAdd(&dicFin,excl,&mel,0);
	  }
	}
      if(dodebug)
	{
	  if(!scaleScore)
	    vstrPrintf(&dbg," failed to extend by a similarity\n");
	  if(isOK)
	    vstrPrintf(&dbg,"        added with score %d as \"%s\"\n",mel.score,excl);
	  else 
	    vstrPrintf(&dbg,"        filtered out\n");
	}
      
		
    }
  
  /* sort the finale */	
  vtextDicSortIndex(&dicFin,1,stcWrd,maxWordsPerSentence);

  /* final tuning for the top hint */	
  if(dicFin.cnt)
    {
      int iDig,iAlf,iScor,score;
      
      /* if the best hit has digits inside, look for the one which 
	 doesn't. Maximum lost in score should be no more than percentDropInHalfanumeric */
      ptr=(char *)vmegId(&dicFin,stcWrd[0]);
      titleCountDigAlf(ptr,&iDig,&iAlf,0,0);
      
      if(	iDig && ( 
			 !(ptr=vstrSkipWords(ptr,1,0)) || /* has just one word */
			 (!vstrSkipWords(ptr,1,0) && vstrSearchSubstring((char *)vmegId(&dicFin,stcWrd[0]),notTheOnlyFinalWord,1,0,1) ) /* has two words and one is meaningless */
			 ) )
	{ /* and has digits */
	  
	  iScor=vtextTreeWordScore(&dicFin,stcWrd[0]); /* get first hint's score */
	  for(ir=1;ir<dicFin.cnt;ir++)
	    {
	      score=vtextTreeWordScore(&dicFin,stcWrd[ir]); /* get first hint's score */			
	      if(score<iScor*percentDropInHalfanumeric/100)
		{ir=dicFin.cnt;break;}
				/* look for digits on the next hit */
	      titleCountDigAlf((char *)vmegId(&dicFin,stcWrd[ir]),&iDig,&iAlf,0,0);
	      if(!iDig)break;
	    }
	  if(ir<dicFin.cnt)
	    { /* able to find one with good score without digits ? */
	      iDig=stcWrd[0];stcWrd[0]=stcWrd[ir];stcWrd[ir]=iDig;
	    }
	  if(dodebug)
	    {
	      vstrPrintf(&dbg,"==================================\n");
	      vstrPrintf(&dbg," The top hint is single halfanumeric word \"%s\" \n",(char *)vmegId(&dicFin,stcWrd[0]));
	      vstrPrintf(&dbg," It was exchanged with \"%s\" \n",(char *)vmegId(&dicFin,stcWrd[0]));
	    }
	}
    }
  
  
  /* final output */	
  titlePrintDic(&out,&dicFin,stcWrd);
  if (dodebug)
    {
      vstrPrintf(&dbg,"==================================\n");
      if (vstrPtr (&out))
	vstrPrintf(&dbg,"Final result\n%s\n\n", vstrPtr (&out)) ;
      else
	vstrPrintf(&dbg,"Final result\nEMPTY\n\n") ;
    }
 
  /* cleanup */
  vtextTreeEmpty(&dicSim);
  vtextTreeEmpty(&dicRes);
  vtextTreeEmpty(&dicWrd);
  vtextTreeEmpty(&dicStc);
  vtextTreeEmpty(&dicFin);
  
  messfree(hitInfo);
  
  vstrPrintf(&dbg,"-----------------------------------------------\n"
	     "-----------------------------------------------\n"
	     
	     "\n\n",hitInfo);
  return vstrPtr(&out);
}


/**********************************************************************/


static BOOL isHorribleTitle (char *title)
{
  char **cpp ;
  
  if (!title || !*title)
    return FALSE ;
   
  for (cpp = horribleTitlesToKill ; *cpp ; cpp++)
    if (pickMatch (title, *cpp))
      return TRUE ;
 
  return FALSE ;
} /* isHorribleTitle */

BOOL titleBestSuggestionANALYZER (vTXT blkp, char *blastHints)
{
  char *titleSuggestions = titleDoBestSuggestionANALYZER (blastHints) ; 
  char hint [vTEXTMAXFILESIZE] ;
  int i, k;
  BOOL isTitle = FALSE ;
  char *ptr, *cp;
  double score = 0 ;
  
  if (!titleSuggestions)
    return FALSE ;
      
  sscanf(titleSuggestions,"%lf",&score) ;
  ptr = vstrSkipWords(titleSuggestions,2,0) ;
  if (!ptr)
    ptr=titleSuggestions;

  cp = strstr (ptr, "\n") ;
  if (cp) *cp = 0 ;
  if (!isHorribleTitle (ptr))
    {
      char **cpp ;
      BOOL addLike = TRUE ;
      isTitle = TRUE ;
      
      if (score >= 200 && !pickMatch (ptr,"*like*"))
	{	  
	  addLike = FALSE ;
	  for (cpp = addLikeToTitle ; *cpp ; cpp++)
	    if (pickMatch (ptr, *cpp))
	      { addLike = TRUE ; break ; }
	}
      
      if (addLike && !pickMatch (ptr,"*like*"))
	vtxtPrintf (blkp, "Blastp_title \"%s like\"\n", ptr) ;
      else
	vtxtPrintf (blkp, "Blastp_title \"%s\"\n", ptr) ;
      
      if (1 && /* bug creates heavy chain chain chain chain */
	  (ptr = strstr (vtxtPtr (blkp), "yosin heavy"))
	  && strcmp (ptr+11, " chain"))
	{
	  vtxtReplaceString (blkp, "yosin heavy", "ZZosin heavy chain") ;
	  vtxtReplaceString (blkp, "ZZosin heavy chain chain", "yosin heavy chain") ;
	  vtxtReplaceString (blkp, "ZZosin heavy chain", "yosin heavy chain") ;
	}
    }
  if (cp) *cp = '\n' ;
  
  for (i = 0, k = 0;;i++)
    {
      hint[k] = titleSuggestions[i];
      if (hint[k] == '\n' || hint[k] == 0)
	{
	  if (k>0)
	    {		   
	      hint[k] = 0;
	      vtxtPrintf (blkp, "Title_hints %s\n", freeprotect (hint)) ;
	    }
	  if (titleSuggestions[i] == 0)
	    break;
	  k = -1;
	}
      k++;
    }
  messfree (titleSuggestions) ;

  return isTitle ;
}














