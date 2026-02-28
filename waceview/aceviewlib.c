#include <wac/ac.h>
#include <vtxt.h>
#include <waceview/cgi.h>

/* Main AceView program. Includes gif and html output (html includes
   lists, main graphics pages and tree displays.
*/
/* 2006:
the system works with frames with several call backs
this is slow
i tried to make it into a single table
but just now only the fiche text loads not any of the other 3 pieces

to activate that code set FR=FALSE (use table)
                             TRUE (use frames, the version that is operational)
june 2007: separate the page in 4 subpages: fgene,fmol, fexp, ffunc
*/
   /*
    * dbname as passed to the cgi 
    * database type ("h" for human, "w" for worm)
    * hostname that runs gifacemblyserver 
    * rpc number to reach server 
    * species name
    * title 
    *
    * '_' characters in species name and title are converted to space. 
    * These fields exist to be displayed to the user.
    */

#define VERSION "v80"
#define DATA_VERSION "v80a"
typedef struct dbConfigStruct {
  const char *dbName ; /* public name of the server, i.e. human */
  char dbSpecies ;
  const char *version ;
  const char *host ;
  int port  ;
  const char *speciesTitle ;
  const char *dbTitle ;
} DBCF ;

/* Never use a number as dbName, other than the taxon id */
static DBCF dbConfig[] =
{ /* dbName  species      version        host    port         speciesTitle            dbTitle */
  { "worm",  'w'     , VERSION /* "v47" */   , "deuce",   2000102, "Caenorhabditis elegans WS190", "Octobre 2008"} ,
  { "worm",  'w'     , VERSION , "ace",   2000102, "Caenorhabditis elegans WS190", "Octobre 2008"} ,
  { "ara",   'a'     , VERSION , "deuce",   2000300, "Arabidopsis thaliana",    "Ara"} ,
  { "ara",   'a'     , VERSION , "ace",   2000300, "Arabidopsis thaliana",    "Ara"} ,

  { "liverpool",   'h'     , "v47"   , "deuce", 2310435, "human",                  "August 2005" } ,

  { "human", 'h'     , VERSION , "deuce", 2310136, "human",                  "April 2007" } ,
  { "human", 'h'     , VERSION , "deuce", 2310236, "human",                  "April 2007" } ,
  { "36a",   'h'     , VERSION , "deuce", 2310136, "human",                  "April 2007" } ,
  { "36a",   'h'     , VERSION , "deuce", 2310236, "human",                  "April 2007" } ,

  /* redirected to 36a */
  { "35g", 'h'     , VERSION , "deuce", 2310136, "human",                  "April 2007" } ,
  { "35g", 'h'     , VERSION , "deuce", 2310236, "human",                  "April 2007" } ,

  { "mouse", 'h'     , VERSION , "deuce", 2320137, "mouse",                  "June 2007" } ,
  { "mouse", 'h'     , VERSION , "deuce", 2320237, "mouse",                  "June 2007" } ,
  { "mm_37", 'h'     , VERSION , "deuce", 2320137, "mouse",                  "June 2007" } ,
  { "mm_37", 'h'     , VERSION , "deuce", 2320237, "mouse",                  "June 2007" } ,

  { "rat", 'r'       , VERSION , "deuce", 2330018, "rat",                  "Sept 2008" } ,

  { "gold",  'h'     , VERSION , "deuce", 2310500, "human",                  "Gold project" } ,

  { "35g",   'h'     , "v47"   , "ace", 2310435, "human",                  "August 2005" } ,
  { "35g",   'h'     , "v47"   , "ace", 2310535, "human",                  "August 2005" } ,

  { "human", 'h'     , VERSION , "ace", 2310136, "human",                  "April 2007" } ,
  { "human", 'h'     , VERSION , "ace", 2310236, "human",                  "April 2007" } ,
  { "36a",   'h'     , VERSION , "ace", 2310136, "human",                  "April 2007" } ,
  { "36a",   'h'     , VERSION , "ace", 2310236, "human",                  "April 2007" } ,

  { "mouse", 'h'     , VERSION , "ace", 2320137, "mouse",                  "June 2007" } ,
  { "mouse", 'h'     , VERSION , "ace", 2320237, "mouse",                  "June 2007" } ,
  { "mm_37", 'h'     , VERSION , "ace", 2320137, "mouse",                  "June 2007" } ,
  { "mm_37", 'h'     , VERSION , "ace", 2320237, "mouse",                  "June 2007" } ,

  { "gold",  'h'     , VERSION , "ace", 2310500, "human",                  "Gold project" } ,

  { "debug", 'h'     , VERSION , "ace", 12345, "human",                  "April 2007" } ,
  { "test2", 'h'     , VERSION , "deuce", 54321, "human",                  "April 2007" } ,
  {  0, 0, 0, 0, 0, 0 }
} ;


typedef enum
{ ZERO_A = 0, TREE_A, DNA_A, PEP_A, PROBE_A, FICHE_A, FGENE_A, FMOL_A, FEXP_A, FFUNC_A, SECTION_A, GENE_A, MRNA_A,
  LIST_A, CLONES_A, TG_SUPPORT_A, GOLD_A, GOLDFILE_A, FASTA_A, ERROR_A, MEGA_A, LINK_A, INFO_A,
  LOG_A, GENE_P, VmRNA_P, Locator_P, LocatorBig_P, HmRNA_P, PROBE_P , LAST_A
} ACTION_TYPE ;

static char *actionName[] = {"ZERO", "TREE", "DNA", "PEP", "PROBE", "FICHE"
			     , "FGENE", "FMOL",  "FEXP", "FFUNC"
			     , "SECTION", "GENE", "MRNA"
			     , "LIST", "CLONES", "TG_SUPPORT", "GOLD", "GOLDFILE", "FASTA", "ERROR"
			     , "MEGA", "LINK", "INFO"
			     , "LOG", "GENE_P", "VmRNA_P", "Locator_P", "LocatorBig_P", "HmRNA_P", "PROBE_P"
} ;
 
typedef enum
  { ZERO_T = 0, GIF_T, BUBBLE_T, SWF_T, swf_T, SVG_T
} PICTURE_TYPE ;

typedef struct pStruct
{
  AC_DB  db ; 
  AC_KEYSET genes, newNames,
    mrnas; 
  AC_OBJ geneBox,
    clone,
    error, /* error message reported by acedb:webquery as a KeySet class name */
    tg,
    genefinder,
    mrna,
    product ;
  PICTURE_TYPE picture_type ;
  ACTION_TYPE action;
  char	*remote_addr, 
    *http_host, 
    *args ;
  char date[256] ;
  int host_offset ; /* is the offset in db_config of the host used last time */
  const char 
    *dbName, 
    /*
     * dbName is the database name as spoken by the web client
     */
    *version, 
    /*
     *  html pages and java scripts 
     */
    *host,
    /*
     * host is the name of the machine running the aceserver for
     * this database
     */
    *speciesTitle ,
    *dbTitle
  ;
  char 
    class[64], 
    org[64],
    query[2048], 
    callback[1024],	
    faction[64],
    actionName[64],
    anchor[128],
    dbSpecies,
		/*
		* dbSpecies is "h:human" or "w:worm".  use this when you need 
		* special code for a database that is too complex to be 
		* configurable.
		*/
    *picCommand,
    goldFile [1024] ; /* direct request for an existing html goldFile */
  unsigned char *pictureBuffer ;
  int pictureBufferSize ; /* size of the pictureBuffer received from the acedb server */
  int port ;
  char	*browser ;
  int width, height ; /* size of the image draing frame in pixcells */
  BOOL orderByNewName ;
  BOOL isGoogle, noRobot ; /* send the no robot tag */
  BOOL noJava ; /* not needed in the frameset pages */
  int  seqstart,seqlen,pagestart,details,newWindow;
  int nonMappedGenes ; /* number of genes without newName */
  int dna1, dna2, dnacolor ; 
  char new_context_name [1024], context_name[1024], context_error [1024] ; /* context */
  AC_KEYSET context_ks ; /* optional: recovered context from previous query */

} PP ;

enum { hide_mRNA = 0x1, show_cDNA=0x2} ;

#define DEFWIDTH        400
static BOOL showTreeObject (PP *ppp) ;
static int getOneKey (PP *ppp) ; /* used in TRRE_A case */
static void setJavaScripts (PP *ppp) ;
static char *fonts = "helvetica" ;

static char *error_message_not_found = 
" Some possible problems and run-arounds"
" are: </span></p>"
" "
" <ol start=1 type=1>"
"  <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;"
"      mso-list:l0 level1 lfo3;tab-stops:list 36.0pt'><b><span style='font-family:"
"      Tahoma;color:#007F7F'>Did you choose the correct species?</span></b><span"
"      style='font-family:Tahoma'><br/>"
"      The human, worm and arabidopsis genes are in separate databases: select"
"      from the second pull down menu, to the left of the main AceView page. </span></li>"
"  <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;"
"      mso-list:l0 level1 lfo3;tab-stops:list 36.0pt'><b><span style='font-family:"
"      Tahoma;color:#007F7F'>We did not recognize what you asked for, please try"
"      another identifier.</span></b><span style='font-family:Tahoma'><br/>"
"      There are many alias names, names may change, we may know this gene under"
"      another name, and despite our efforts to track gene names (see the <a"
"      href=\"HelpJan.htm#names\">FAQ</a>), a gene name from the previous build may"
"      not be recognized. But AceView should <span style='color:#3366FF'>recognize</span>"
"      a number of other NCBI identifiers that point to the gene: <span"
"      style='color:#3366FF'>LocusLink/GeneID, GenBank mRNA or EST accession,"
"      RefSeq accession, OMIM accession, or Unigene cluster name Hs<span"
"      class=GramE>.#</span>.</span> You might also <span style='color:#3366FF'>try"
"      a cDNA clone name, or a functional or descriptive name</span>. If you are"
"      interested by function, try the <span style='color:#3366FF'>Pfam query</span>."
"      But please do NOT look for a gi number, a protein or NP accession, an STS"
"      marker name (go to LocusLink then back to us), an older AceView name (e.g."
"      <i>Hs11_34081_29_1_661</i>, only archived on our ftp site), a WormPep"
"      accession, an NT contig or a genomic sequence accession since these"
"      queries are not supported in AceView. We do not yet support either queries"
"      by chromosome or by region. If you are interested in <i>IGH@, IGK@, IGL@,"
"      TRB@ or TRG@</i>, call for one gene in the region, click on the graphic"
"      &quot;Full page&quot; and zoom out from there (left icon) to see the whole"
"      area. </span></li>"
"  <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;"
"      mso-list:l0 level1 lfo3;tab-stops:list 36.0pt'><b style='mso-bidi-font-weight:"
"      normal'><span style='font-family:Tahoma;mso-fareast-font-family:\"Times New Roman\";"
"      mso-bidi-font-family:\"Courier New\";color:teal'>Please remain concise in"
"      your queries</span></b><span style='font-family:Tahoma;mso-fareast-font-family:"
"      \"Times New Roman\";mso-bidi-font-family:\"Courier New\";color:black'>: we do"
"      not search for partial matches and do not ignore spelling mistakes, like"
"      Google does. All you typed had to be found in the same gene entry. So,"
"      type the minimal number of words or word roots that translate your"
"      thought. </span><span style='font-family:Tahoma'></span></li>"
"  <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;"
"      mso-list:l0 level1 lfo3;tab-stops:list 36.0pt'><b><span style='font-family:"
"      Tahoma;color:#007F7F'>Are you sure you spelt your word(s) correctly?</span></b><span"
"      style='font-family:Tahoma'><br/>"
"      The character case does not matter (BRCA1 == brca1 == BrCa1) but we do not"
"      look for approximate matches: spelling mistakes are deadly. We extend what"
"      you typed (with any character) until we find a match (see the <a"
"      href=\"HelpQuery#search\" target=\"_top\">help</a> on query), so, for best"
"      results, avoid useless constraints: type only the roots of meaningful"
"      words that undoubtedly should be found in your gene. If you are unsure of"
"      the spelling, type only what you are sure of: prefer <i>lecit choleste</i>"
"      to <i>lecitin cholestelol</i> and <i>phosphor</i> (hardly ambiguous in"
"      fact) to <i>phosphorilation</i>. When in doubt, prefer a space or a"
"      question mark to an uncertain character. Prefer singular to plural (e.g."
"      use <i>claudin</i> rather than <i>claudins</i>) and avoid using dashes or"
"      dots in names unless strictly necessary. For example, <i>NF kappa b</i>, <i>NFkappaB</i>"
"      or <i>NF-KB</i> will not give the same answers: in human build 33, you get"
"      205 genes for the first query, 13 for the second and 11 for the third. </span><span"
"      style='font-family:Tahoma;mso-fareast-font-family:\"Times New Roman\";"
"      mso-bidi-font-family:\"Courier New\";color:black'>Keep spaces between words"
"      as usual: RNA pol III is more frequently used than RNApolIII, ubiquitin1"
"      will not be recognized.</span><span style='font-family:Tahoma'></span></li>"
"  <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;"
"      mso-list:l0 level1 lfo3;tab-stops:list 36.0pt'><b><span style='font-family:"
"      Tahoma;color:#007F7F'>If you come directly from LocusLink/Gene,</span></b><b"
"      style='mso-bidi-font-weight:normal'><span style='font-family:Tahoma;"
"      color:teal'> or were interested in a gene associated to a disease</span></b><span"
"      style='font-family:Tahoma'><br/>"
"      it is possible that we did not succeed yet in mapping the sequence(s) associated"
"      to the locus unambiguously, or that this locus has no associated sequence because"
"      it has not been cloned yet. In this case, OMIM or LocusLink would be more"
"      informative: AceView is only a complementary sequence oriented view. In"
"      the top Query box, search All databases at NCBI, or OMIM, or LocusLink/Gene."
"      </span></li>"
"  <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;"
"      mso-list:l0 level1 lfo3;tab-stops:list 36.0pt'><b><span style='font-family:"
"      Tahoma;color:#007F7F'>If you found nothing with a short query</span></b><span"
"      style='font-family:Tahoma'><br/>"
"      and no spelling mistake, we apologize, and suggest that you <span"
"      style='color:#3366FF'>try one of the two other ways to query AceView:"
"      through PFAM or through Blast</span>. Some words, such as actin (also"
"      found in acting and interacting) or transcription factor, are so frequent"
"      in the literature that our indexing system gets saturated, and you get no"
"      answer. </span></li>"
" </ol>"
" "
" <p class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;"
" margin-left:72.0pt'><b><span style='font-family:Tahoma;color:#007F7F'>If you"
" are interested by </span></b><b><span style='font-family:Tahoma;color:teal'>homologs</span></b><b"
" style='mso-bidi-font-weight:normal'><span style='font-family:Tahoma;"
" color:teal'> of a gene from another species</span></b><span style='font-family:"
" Tahoma'> (like for example <i>NF Kappa B</i>), you will get your answer,"
" complete and quantified, by entering the sequence from your favorite gene in"
" the &quot;query AceView via <span style='color:#3366FF'>Blast</span>&quot; that"
" we provide (thanks to <st1:place w:st=\"on\">S. Altschul</st1:place> and the"
" Blast team). </span></p>"
" "
" <p class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;"
" margin-left:72.0pt'><b><span style='font-family:Tahoma;color:#007F7F'>If you"
" asked for a DNA sequence, such as an oligonucleotide,</span></b><span"
" style='font-family:Tahoma'> directly in the query box, this will unfortunately"
" not bring an answer. Sorry, we do not yet provide a direct AceView align"
" service. Please try the <span style='color:#3366FF'>Blast</span> query in the"
" left <span class=GramE>margin,</span> it will align your DNA or protein against"
" all the AceView transcripts in the selected species</span></p>"
" "
" <p class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;"
" margin-left:72.0pt'><b><span style='font-family:Tahoma;color:#007F7F'>If you"
" are interested by genes belonging to a well known family, or that have a known"
" function </span></b><span style='font-family:Tahoma'><br/>"
" your best chances are to use the &quot;<span style='color:#3366FF'>query Pfam</span>&quot;"
" on the left. Here again you will benefit from the work of the Pfam and InterPro"
" teams, world experts in defining, recognizing, and describing conserved"
" families of proteins. We have just run their programs (slow and heavy) on all"
" the AceView transcripts and have compiled the data for you, so that you get"
" complete answers over the entire genome fast. </span></p>"
" "
" <ol start=7 type=1>"
"  <li class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;"
"      mso-list:l0 level1 lfo3;tab-stops:list 36.0pt'><b><span style='font-family:"
"      Tahoma;color:#007F7F'>If you requested a particular mRNA </span></b><b><span"
"      style='font-family:Tahoma;color:teal'>variant by name</span></b><b><span"
"      style='font-family:Tahoma;color:#007F7F'> and suffix (a, b or c<span"
"      class=GramE>)</span></span></b><span style='font-family:Tahoma'><br/>"
"      from a previous build and got no answer, this is expected. From the"
"      December 03 release on, mRNA variants have a unique identifier, to avoid"
"      confusion of having objects with different sequences named the same in"
"      successive versions of AceView.</span><span style='font-family:Tahoma'>And we"
"      do not have room to maintain more than the last two builds on our server. Archives"
"      may be available from Jim Kent, at UCSC.</span></li>"
" </ol>"
" "
" <p class=MsoNormal style='mso-margin-top-alt:auto;mso-margin-bottom-alt:auto;"
" margin-left:18.0pt'><span style='font-family:Tahoma'>You are welcome to <a"
" href=\"mailto:mieg@ncbi.nlm.nih.gov\">report"
" problems</a> or ask for help. We would be happy to try to fix any reproducible"
" bug as well as to add desirable features.</span></p>"
"" ;

static char *error_message_saturation = 
" If you got this message, you could try to "
"rephrase or extend the query. But it may "
"not always work, sorry. The problem is that some words, such as actin (also "
"found in acting and interacting) or transcription factor, are so "
"frequent in the literature that our indexing system gets saturated, and "
"you get no answer. <i>Actin filament</i> would give results, as would <i>actin "
"polymer</i> or <i>actin cytoskel</i> or <i>actin interacting protein</i>."
"But please remember to choose the best possible query among the three "
"types offered:"
"<ul>"
"<li><font color=#007f7f><b>  If you are interested by homologs</b></font><br/>\n of a gene from another species "
"(like for example NF Kappa B), you will get a much better answer, "
"complete and quantified, by entering the sequence from your favorite "
"gene in the \"query AceView via Blast\" we provide (thanks to S. "
"Altschul). Indeed, if you search in the query box directly, you only "
"search our annotation and titles, and if your gene homolog has an "
"interesting phenotype, we may have described it by its phenotype rather "
"than its homology."
"<li><font color=#007f7f><b>  If you are interested by genes belonging to a well known family,</b></font><br/>\n your "
"best chances are to use the \"query AceView via Pfam\" on the left. Here "
"again you will benefit from the work of the Pfam and Interpro teams, "
"world experts in defining, recognizing, and describing conserved "
"families of proteins. We have just run their programs (slow and heavy) "
"on all the AceView transcripts and have compiled the data for you, so "
"that you get complete answers over the entire genome fast."
"</ul>"
"<hr></hr>"
 ;

static char *pgetenv (char *a)
{
  char *s = getenv(a) ;
  if (!s)
    s = "(null)" ;
  return s ;
}

static char *logPrefix (PP *ppp, char *desc)
{
  char *s ;
  char *q, *rem, *p, *ref ;

  q = pgetenv("QUERY_STRING") ;
  rem = pgetenv("REMOTE_ADDR") ;
  p = pgetenv("HTTP_X_FWD_IP_ADDR") ;
  ref = pgetenv("HTTP_REFERER") ;
  s = messprintf("// aceview: %d %s %s %s %s %s %s %s", (int)time(0), ppp ? ppp->date : "(time)", desc,q,rem,p,ref, ppp->browser ? ppp->browser : "NULL-browser") ;
  return s ;
}

/******************************************************************************/

static char *url_sanitize (char *txt)
{
  char *cp0, *cp1 ;
  
  cp0 = strnew (txt, 0) ;
  url_decode_inplace(cp0) ;
  for (cp1 = cp0 ; *cp1 ; cp1++)
    if (*cp1 == '<' || *cp1 == '>') *cp1 = ' ' ;
  while ((cp1 = strstr (cp0, "script")))
    memset (cp1, ' ', 6) ;
  return cp0 ;
}

/******************************************************************************/
/***************************************************************************/

static void showMOI (PP * ppp )
{
  AC_TABLE gIntMap ;
  AC_OBJ oGene = ppp->geneBox, oGenefinder = 0, oTg = 0;
  const char *ttt ;

  printf ("\n<form  name=\"MOI\">\n") ;
  printf ("  <input type=hidden  name=\"currentQuery3\" value=\"%s\">\n", ppp->args) ;
  printf ("  <input type=hidden  name=\"currentDb\" value=\"%s\">\n", ppp->dbName) ;
  printf ("  <input type=hidden  name=\"currentVersion\" value=\"%s\">\n",  VERSION) ;
  printf ("  <input type=hidden  name=\"currentData\" value=\"%s\">\n",  DATA_VERSION) ;
  if (ppp->host_offset > 0)
    printf ("  <input type=hidden  name=\"currentHost\" value=%d>\n", ppp->host_offset) ;
  printf ("  <input type=hidden  name=\"currentSwf\" value=\"%s\">\n", "s") ; /* could be "B" */

  if (*ppp->new_context_name) 
    printf ("  <input type=hidden  name=\"ctx\" value=\"%s\">\n", 
	    ppp->new_context_name) ;

  if (ppp->mrna && !ppp->geneBox)
    {
      oTg = ac_tag_obj (ppp->mrna, "From_gene", 0) ;
      oGenefinder = ac_tag_obj (ppp->mrna, "From_prediction", 0) ;

      if (oTg)
        {
          oGene = ppp->geneBox =  ac_tag_obj (oTg, "Gene", 0) ;
          ppp->tg = oTg ;
        }
      if (!ppp->geneBox &&  oGenefinder)
        {
          oGene = ppp->geneBox = ac_tag_obj (oGenefinder, "Model_of_gene", 0) ;
          ppp->genefinder = oGenefinder ;
        }
      ppp->geneBox=ac_tag_obj (ppp->mrna, "From_gene", 0) ;
    }

  if (ppp->geneBox && !ppp->mrna)
    {
      if (!ppp->tg) ppp->tg = ac_tag_obj (ppp->geneBox, "Transcribed_gene", 0) ;
      if (ppp->tg)  ppp->mrna=ac_tag_obj (ppp->tg, "mRNA", 0) ;
      if (!ppp->mrna)
        {
          if (!ppp->genefinder) ppp->genefinder = ac_tag_obj (ppp->geneBox, "Genefinder", 0) ;
          if (ppp->genefinder)  ppp->mrna=ac_tag_obj (ppp->genefinder, "Predicted_mRNA", 0) ;
        }
    }

  if (ppp->geneBox)
    printf ("  <input type=hidden  name=\"currentGENE\" value=\"%s\">\n", ac_name (ppp->geneBox)) ;
  else
    printf ("  <input type=hidden  name=\"currentGENE\" value=0>\n") ;
  if (ppp->mrna)
    printf ("  <input type=hidden  name=\"currentMRNA\" value=\"%s\">\n", ac_name (ppp->mrna)) ;
  else
    printf ("  <input type=hidden  name=\"currentMRNA\" value=0>\n") ;
  if (ppp->product)
    printf ("  <input type=hidden  name=\"currentProduct\" value=\"%s\">\n", ac_name (ppp->product)) ;
  else
    printf ("  <input type=hidden  name=\"currentProduct\" value=0>\n") ;

  if (ppp->clone)
    printf ("  <input type=hidden  name=\"currentClone\" value=\"%s\">\n", ac_name (ppp->clone)) ;

  if (ppp->action == FICHE_A && ace_lower((int)ppp->class[0]) == 'g')
    printf ("  <input type=hidden  name=\"hasMap\" value=\"1\">\n") ;
  printf ("  <input type=hidden  name=\"currentAction\" value=\"%s\">\n"
	  ,  ppp->action >= ZERO_A && ppp->action < LAST_A ? actionName[ppp->action] : "0"
	  ) ;

  ttt = oGene ? ac_tag_printable (oGene, "LocusLink", 0) : 0 ;

  if (!ttt) ttt = ac_tag_printable (ppp->geneBox, "Genefinder", 0) ;
  if (ttt) printf ("  <input type=hidden  name=\"currentLOCUS\" value=\"%s\">\n", ttt) ;

  ttt = oGene ? ac_tag_printable (oGene, "LocusId", 0) : 0 ;
  if (!ttt && (oGenefinder = ac_tag_obj (ppp->geneBox, "Genefinder", 0)))
    ttt = ac_tag_printable (oGenefinder, "LocusId", 0) ;
  if (ttt)printf ("  <input type=hidden  name=\"currentLOCUSID\" value=\"%s\">\n", ttt) ;

  printf ("  <input type=hidden  name=\"currentCosmid\" value=\"%s\">\n",
          oGenefinder ? ac_name (ac_tag_obj (oGenefinder, "Source", 0)) : "NULL") ;
  gIntMap = ac_tag_table (ppp->geneBox, "IntMap", 0) ;
  if (gIntMap && gIntMap->cols >= 3)
    {
      const char * chrm ;

      chrm = ac_table_printable (gIntMap, 0, 0, "") ;
      if (strstr (chrm, "_"))chrm = strstr (chrm, "_")+1 ;
      printf ("  <input type=hidden  name=\"currentCHROM\" value=\"%s\">\n", chrm) ;
      printf ("  <input type=hidden  name=\"currentX1\" value=\"%d\">\n", ac_table_int (gIntMap, 0, 1, 0)) ;
      printf ("  <input type=hidden  name=\"currentX2\" value=\"%d\">\n", ac_table_int (gIntMap, 0, 2, 0)) ;
    }
  printf ("  <input type=hidden  name=\"aceviewVersion\" value=\"103\">\n") ;
  ac_free (gIntMap) ;

  printf ("</form >\n") ;
} /* showMOI */

/***************************************************************************/
/***************************************************************************/

static void setMetaDescription (PP *ppp) 
{
  printf ("  <META NAME=\"title\"\n"
	   " CONTENT=\"\n"
	  "AceView: %s:%s a comprehensive annotation"
	  " of human, mouse and worm genes with"
	  " mRNAs or EST\">\n\n"
	  , ppp->class, ppp->query
	  ) ;

  printf ("<META NAME=\"keywords\"\n"
	  " CONTENT=\"\n"
	  "AceView, genes, Acembly, AceDB, Homo sapiens, Human,\n"
	  " nematode, Worm, Caenorhabditis elegans , WormGenes, WormBase, mouse,\n"
	  " mammal, Arabidopsis, gene, alternative splicing variant, structure,\n"
	  " sequence, DNA, EST, mRNA, cDNA clone, transcript, transcription, genome,\n"
	  " transcriptome, proteome, peptide, GenBank accession, dbest, RefSeq,\n"
	  " LocusLink, non-coding, coding, exon, intron, boundary, exon-intron\n"
	  " junction, donor, acceptor, 3'UTR, 5'UTR, uORF, poly A, poly-A site,\n"
	  " molecular function, protein annotation, isoform, gene family, Pfam,\n"
	  " motif ,Blast, Psort, GO, taxonomy, homolog, cellular compartment,\n"
	  " disease, illness, phenotype, RNA interference, RNAi, knock out mutant\n"
	  " expression, regulation, protein interaction, genetic, map, antisense,\n"
	  " trans-splicing, operon, chromosome, domain, selenocysteine, Start, Met,\n"
	  " Stop, U12, RNA editing, bibliography"
 	  "\">\n" ) ;
  
  /* this text is duplicated in index.html */
  printf ("<META NAME=\"Description\" \n"
	  " CONTENT= \"\n"
	  "AceView offers a comprehensive annotation of human, mouse and nematode genes\n"
	  " reconstructed by co-alignment and clustering of all publicly available\n"
	  " mRNAs and ESTs on the genome sequence. Our goals are to offer a reliable\n"
	  " up-to-date resource on the genes, their functions, alternative variants,\n"
	  " expression, regulation and interactions, in the hope to stimulate\n"
	  " further validating experiments at the bench\n"
 	  "\">\n") ;
  printf ("<meta http-equiv=\"content-type\" content=\"text/html;charset=iso-8859-1\">\n") ;
  if (0) printf ("<meta http-equiv=\"expires\" content=\"Mon, 31 Dec 2007 23:59:59 UTC\">\n") ;
  {
    char tbuf[200] ; 
    timeShowFormat (timeNow() + 30*86400, "%a, %d %b %y %H:%M:%S UTC", tbuf, 199) ;
    printf ("<meta http-equiv=\"expires\" content=\"%s\">\n", tbuf) ;
  }
  printf ("<meta name=\"author\"\n"
	  " content=\"Danielle Thierry-Mieg and Jean Thierry-Mieg,\n"
	  " NCBI/NLM/NIH, mieg@ncbi.nlm.nih.gov"
	  "\">\n") ;
  /* content privacy policy needed to be able to set a cookie */
  printf ("<meta http-equiv=\"P3P\" content='CP=\"IDC DSP COR CURa ADM OUR IND PHY ONL COM STA\"'>") ;
 if (0) printf ("<meta name=\"Content-Language\" content=\"us\">\n") ;

 if (0) printf ("<meta name=\"robots\" content=\"index,nofollow\">\n") ;
} /* setMetaDescription */
  
/***************************************************************************/

static void setJavaScripts (PP *ppp)
{
  /* send start up vars needed by javascript/main.ppp->version.js and include the script */

  printf ("\n<link rel='stylesheet'  name='fixuplink' href='jscript/aceview.%s.css' type='text/css'>\n"
	  , ppp->version
	  ) ;
  
  printf ("\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
          "  <!--\n"
          "    var myurl=\"%s\" ;\n"
          "    var db=\"%s\" ;\n"
          "    var doSwf=\"s\" ;\n"
          "    var classe=\"%s\" ;\n"
	  , ppp->callback, ppp->dbName, ppp->class) ;

  /*
  if (ppp->gene)
    printf ("   var gene=\"%s\" ;\n" , ac_name (ppp->gene)) ;
  if (ppp->mrna)
    printf ("   var mrna=\"%s\" ;\n" , ac_name (ppp->mrna)) ;
  */
  printf ("  //-->\n"
	  "</script>\n"
	  ) ;

  printf ("\n<script language = \'JavaScript\' type=\'text/javascript\' "
	  "src=\'jscript/basic.%s.js\'> "
	  "\n</script>\n"
	  , ppp->version
	  ) ;

  printf ("\n<script language = \'JavaScript\' type=\'text/javascript\' "
	  "src=\'jscript/main.%s.js\'>"
	  "\n</script>\n"
	  , ppp->version
	  ) ;

  printf ("\n<script language = \'JavaScript\' type=\'text/javascript\' "
	  "src=\'jscript/utils.%s.js\'>"
	  "\n</script>\n"
	  , ppp->version
	  ) ;
  printf ("\n<script language = \'JavaScript\' type=\'text/javascript\' "
	  "src=\'jscript/toggle.%s.js\'>"
	  "\n</script>\n"
	  , ppp->version
	  ) ;

  if (0) printf ("\n<script language = \'JavaScript\' type=\'text/javascript\' "
	  "src=\'jscript/swfobject.%s.js\'>"
	  "\n</script>\n"
	  , ppp->version
	  ) ;

  printf ("\n<script language = \'JavaScript\' type=\'text/javascript\' "
	  "src=\'jscript/cssQuery-p.%s.js\'>"
	  "\n</script>\n"
	  , ppp->version
	  ) ;

  printf ("\n<script language = \'JavaScript\' type=\'text/javascript\' "
	  "src=\'jscript/notify.%s.js\'>"
	  "\n</script>\n"
	  , ppp->version
	  ) ;

  printf ("\n<script language = \'JavaScript\' type=\'text/javascript\' "
	  "src=\'jscript/fanfold.%s.js\'>"
	  "\n</script>\n"
	  , ppp->version
	  ) ;

  printf ("\n<script language = \'JavaScript\' type=\'text/javascript\' "
	  "src=\'jscript/egene.%s.js\'>"
	  "\n</script>\n"
	  , ppp->version
	  ) ;

  printf ("\n<script language = \'JavaScript\' type=\'text/javascript\' "
	  " src=\'jscript/wz_jsgraphics.%s.js\'>"
	  "\n</script>\n"
	  , ppp->version
	  ) ;

} /* setJavaScripts */

/***************************************************************************/

static void standardHtmlHeader (PP *ppp, BOOL hasBody)
{
  printf ("Content-type: text/html\n\n") ;
  printf ("<html>\n") ;

  printf ("<head>\n") ;
  if (ppp->noRobot) printf ("<META NAME=\"ROBOTS\" CONTENT=\"NOINDEX, NOFOLLOW\">\n") ;

  printf ("  <title>AceView: %s:%s, a comprehensive annotation"
	  " of human, mouse and worm genes with"
	  " mRNAs or ESTsAceView.</title>\n\n"
	  , ppp->class, ppp->query
	  ) ;

  setMetaDescription (ppp) ;
  if (! ppp->noJava) setJavaScripts (ppp) ;
  printf ("</head>\n\n") ;

  if (hasBody)
    {
      if (ppp->geneBox)
	showMOI (ppp) ;

      printf ("<body ") ;
      /*
	if (1 || ppp->action == FGENE_A || ppp->action == FMOL_A || ppp->action == FEXP_A || ppp->action == FFUNC_A)
	printf (" onLoad=aceViewOnLoad() ") ;
      */
      printf (" onunload=aceViewOnUnload() >\n\n") ;
      if (ppp->action  != FICHE_A && 
	  ppp->action  != FGENE_A && ppp->action  != FMOL_A && ppp->action  != FEXP_A && ppp->action  != FFUNC_A &&
	  ppp->action  != ERROR_A 
	  )
	{
	  if (strstr (ppp->browser, "Mozilla/") && !strstr (ppp->browser, "MSIE") && !strstr (ppp->browser, "Netscape6") )
	    
	    printf ("<layer id=\'javascriptRequired\' name=\'javascriptRequired\' STYLE=\'position: absolute ;\' >\n  <FONT size=-1 color=RED>\n<b>To see this page, please, enable javaScript and styleSheets on your browser: <br/>\nEdit -> Preferences -> Advanced -> Enable JavaScript.<br/>\nReload!</b><br/>\n  </FONT>\n</layer>\n") ;
	  else if (1)
	    printf ("<div id=\'javascriptRequired\' name=\'javascriptRequired\' STYLE=\'position: absolute ;\'>\n<FONT size=-1 color=RED><b>To see this page, please, enable javaScript and styleSheets on your browser: Tools -> Internet Options ... <br/>\nReload !<br/>\n</b>\n  </FONT>\n</div>\n" ) ;
	  printf ("\n\t<script language = \'JavaScript\' type=\'text/javascript\'>\n"
		  "\t\t if (document.layers)document.javascriptRequired.visibility=\"hide\" ;\n"
		  "\t\t if (!document.all && document.getElementById)document.getElementById (\"javascriptRequired\").style.visibility=\"hidden\" ;\n"
		  "\t\t if (document.all)document.all[\"javascriptRequired\"].style.visibility=\"hidden\" ;\n"
		  "\t</script>\n") ;
	} 
    }
} /* standardHtmlHeader */

/***************************************************************************/

#ifdef JUNK
/* This was a draft from vahan */
/* needed for Opening Closing paragraph in a Browser */
/*
  static void standardFicheHeader (PP * ppp)
  {
  printf ("<STYLE TYPE='text/css'>.hiddentext
  {display:none ; }</STYLE>") ;
  printf (      "<SCRIPT language=JavaScript>\n"
  "<!--\n"
  "    function ficheSectionToggler (snm)\n"
  "
  {\n"
  "        whichEl = gObjectStyle (snm) ;\n"
  "        whichEl.display = (whichEl.display == \"block\" ) ? \"none\" : \"block\" ;\n"
  "              isopen = (whichEl.display == \"block\") ? 1 : 0 ; \n"
  "              layerSet (\"button\"+snm, \"-\", \"-\", ficheSectionButton (snm, isopen), \"-\", \"-\") ; \n"
  "    }\n"
  "     function ficheSectionButton (snm, isopen)\n"
  "
  {\n"
  "         return \"<SPAN ID = 'button\"+snm+\"'>\"+ \n"
  "             \"<a href = 'javascript:ficheSectionToggler (\\\"\"+snm+\"\\\") ;'>\"+ \n"
  "             \"<small>< \"+ (isopen ? \"-\" : \"+\")+\" ></small></a>\"+ \n"
  "             \"</SPAN>\" ; \n"
  "     }\n"
  "    function ficheSectionStart (snm)\n"
  "
  {\n"
  "        if (isBrws (\"ns4\"))return ;\n"
  "        ttxt = ficheSectionButton (snm, 0)+ \n"
  "             \"<DIV name = '\"+snm+\"' ID='\"+snm+\"' CLASS='hiddentext' >\" ;\n"
  "        document.write (ttxt) ;\n"
                "    }\n"
                "    function ficheSectionClose ()\n"
                "
                {\n"
                "        if (isBrws (\"ns4\"))return ;\n"
                "        document.write (\"</DIV\") ;\n"
                "    }\n"
                "//-->\n"
                "</SCRIPT>\n" ) ;
                } */
#endif

/***************************************************************************/
/* display the various error messages
 *
 * Note that the html_encode and url_encode calls leak memory, but we are
 * about to exit anyway, so it is no problem.
 */

static void reportError (PP *ppp, char *s)
{
  char *t;
  t = url_sanitize(s) ;
  t = html_encode(t);

  printf ("<table border=0><tr><td><a href=\"index.html?%s\" target=\"_top\"><IMG src=\"images/aceview.%s.gif\" width=134 border=0></a></td></tr></table><br/>\n<br/>\n", ppp->dbName,ppp->version) ;

  if (0 && ! ppp->error)
    {
      printf("<h2>Unknown Error</h2>") ;
      printf("<p>\n\nSorry, the query --");
      printf("<b>%s</b> -- ",t);
      printf("returned no gene and no specific error message.\n") ;
      printf("This problem has been logged and we will try to fix it in the next few days.") ;
      printf ("Please try something else like a gene name or sequence accession of a phenotype.\n") ;  
    }
  else if ( ppp->error && strstr (ac_name(ppp->error), "too many hits"))
    {
      printf("<h2>Too Many Hits</h2>") ;
      printf("<p>\n\nSorry, the query --");
      printf("<b>%s</b> -- ",t);
      printf("is not specific enough. ");
      
      printf(error_message_saturation, ppp->version) ;
    }
  else if (! ppp->error || strstr (ac_name(ppp->error), "locus not found"))
    {
      printf("<h2>Locus Not Found</h2>") ;
      printf("<p>\n\nSorry, but we could not find any gene matching --");
      printf("<b>%s</b> -- ",t);
      printf("in the AceView system. ");
      
      printf(error_message_not_found, ppp->version) ;
    }
} /* reportError */

/***************************************************************************/
/***************************************************************************/
/*
 * query context is currently stored as a temp file on the database server.
 * This is the directory where we keep the context.  You need to set
 * something up to delete files from this directory pretty regularly
 */
#define CTXDIR "/home/mieg/ee/SERVER/QUERY_TMP"

static BOOL contextCreate (PP *ppp)
{
  char            host[100], day[50];
  time_t          t;
  struct tm      *tp;
  
  host[0] = '?';
  host[1] = 0;
  gethostname (host, sizeof(host));
  t = time (0);
  tp = gmtime (&t);
  strftime (day, 50, "%j%H", tp);
  
  sprintf (ppp->new_context_name, "ctx-%s-%s-%d", day, host, (int) getpid()) ;
  return TRUE ; /* cannot fail i think */
}

static BOOL contextSave (PP *ppp, AC_KEYSET aks)
{
  char fNam [2048] ;
  
  if (aks && ac_keyset_count (aks) &&
      contextCreate (ppp))
    {
      sprintf (fNam, "%s/%s", CTXDIR, ppp->new_context_name) ;
  
      if (ac_keyset_write (aks, fNam)) 
	return TRUE ;
      else
	*ppp->new_context_name = 0 ;
    }

  return FALSE ;
} /* contextSave */

static BOOL contextLoad (PP *ppp)
{
  char  *cp, buf[2056];
  
  if (!*ppp->context_name)
    return FALSE ;
  
  for (cp = ppp->context_name ; *cp ; cp++)
    if (!isalnum((int) *cp) && (*cp != '.') && (*cp != '_') && (*cp != '-'))
      return FALSE ;
  
  if (strchr(ppp->context_name, '/'))
    return FALSE ;
  sprintf (buf,"%s/%s", CTXDIR, ppp->context_name);
  
  ppp->context_ks = ac_read_keyset (ppp->db, buf, NULL) ;
  if (! ppp->context_ks)
    {
      *ppp->context_name = 0 ;
      sprintf (ppp->context_error,
	       "<h1>Query results expired </h1>\n"
	       "<p>\n\nQuery results are only stored for about an hour after you make "
	       "the query.  Click <a href='av.cgi/index.html?%s'>here</a> to start a new query\n"
	       , ppp->dbName);
      return FALSE ;
    }
  return TRUE ;
} /* contextLoad  */

/***************************************************************************/
/***************************************************************************/

static void executeAction (PP *ppp)
{
  unsigned char *fiche = 0 ;

  switch (ppp->action)
    {
    case FGENE_A:
    case FMOL_A:
    case FEXP_A:
    case FFUNC_A:
      ac_command (ppp->db, messprintf (" // %s start %s  %s aceviewlib.2019_02_17" 
				       , actionName[ppp->action]
				       , ppp->class, ppp->query), 0, 0) ;

      fiche = ac_command (ppp->db
			 , messprintf("view -v %s -c %s -n \"%s\" -xml -preview2"  
				      , actionName[ppp->action]
				      , ppp->class, ppp->query) 
			  , 0, 0) ;
      ac_command (ppp->db, messprintf (" // fiche done %s  %s" , ppp->class, ppp->query), 0, 0) ;
      break ;
    case FICHE_A:
      ac_command (ppp->db, messprintf (" // fiche start %s  %s" , ppp->class, ppp->query), 0, 0) ;
      
      fiche = ac_command (ppp->db
			  , messprintf("view -v %s -c %s -n \"%s\" -xml -preview2",  
				       (ppp->faction[0] ) ?  ppp->faction : "fiche" , 
				       ppp->class, ppp->query) 
			  , 0, 0) ;
      ac_command (ppp->db, messprintf (" // fiche done %s  %s" , ppp->class, ppp->query), 0, 0) ;
      break ;
    case SECTION_A:
      ac_command (ppp->db, messprintf (" // section start %s  %s  %s" , ppp->faction , ppp->class, ppp->query), 0, 0) ;
      
      fiche = ac_command (ppp->db
			  , messprintf("view -v section  -c %s -n \"%s\" -xml  -preview2  -s %s",  
				       ppp->class, ppp->query, ppp->faction) 
			  , 0, 0) ;
      ac_command (ppp->db, messprintf (" // section done %s  %s  %s" , ppp->faction , ppp->class, ppp->query), 0, 0) ;

      break ;
    case FASTA_A:
      fiche = ac_command (ppp->db
			 , messprintf ("view -v fasta -c %s -n \"%s\" -xml", ppp->class, ppp->query)
			 , 0, 0) ;
      break ;
    case CLONES_A:
      ac_command (ppp->db, messprintf (" // %s start %s  %s" 
				       , actionName[ppp->action]
				       , ppp->class, ppp->query), 0, 0) ;
      fiche = ac_command (ppp->db
			 , messprintf ("view -v CLONES -c %s -n \"%s\" -xml", ppp->class, ppp->query)
			 , 0, 0) ;
      ac_command (ppp->db, messprintf (" // fiche done %s  %s" , ppp->class, ppp->query), 0, 0) ;
      break ;
    case TG_SUPPORT_A:
      fiche = ac_command (ppp->db
			 , messprintf ("view -v SUPPORT -c %s -n \"%s\" -xml -x1 %d -x2 %d"
				       , ppp->class, ppp->query, ppp->dna1, ppp->dna2)
			 , 0, 0) ;
      break ;
    case GOLD_A:
      fiche = ac_command (ppp->db
			 , messprintf ("view -v Gold -c %s -n \"%s\" -xml", ppp->class, ppp->query)
			 , 0, 0) ;
      break ;
    case LOG_A:
      ac_command (ppp->db
		  , messprintf ("comment LOG %s %s %s", ppp->remote_addr, ppp->query, ppp->browser ? ppp->browser : "NULL-browser")
		  , 0, 0) ;
      fiche = (unsigned char *) strnew ("AceView", 0) ;
      break ;
    default:
      break;
    }
  if (fiche && strncmp ((char *)fiche, "// Error", 8)) /* ok case */
    {
      printf ("%s\n", fiche) ;  /* the actual fiche */
      messfree (fiche) ;
    }
  else /* error case, no fiche or fiche reports an error */
    {
      if (getOneKey (ppp))
	showTreeObject (ppp) ;
      printf ("</body></html>\n") ;
    }
}

/***************************************************************************/

static void showInfoFrameSet (PP *ppp)
{
  int mrnaMode = 0, directMode = 0 ;

  if (!strcasecmp (ppp->class, "mRNA"))
    mrnaMode = 1 ;
  if (!strcasecmp (ppp->class, "cDNA_clone"))
    directMode = 1 ;

  if (!strcmp (ppp->version,"v47") || directMode) /* old system */
    {
      printf ("<FORM NAME=\"vsim2\">\n") ;

      printf ("<FRAMESET COLS=\"80%%, *\" name=\"Info\">\n") ;
                                /* FICHE Column */
  
      if (directMode)
	{
	  printf (" <frame MARGINWIDTH=0 MARGINHEIGHT=0 SCROLLING=\"Auto\" FRAMEBORDER=\"0\" FRAMESPACING=0 name=\"Fiche\"") ;
	  printf (" src=\"av.cgi?db=%s&c=%s&a=%s&wh=%d&l=%s\">\n", ppp->dbName, ppp->class, ppp->faction, ppp->height,url_encode(ppp->query)) ;
	}
      else if (!mrnaMode)
	{
	  
	  {
	    printf (" <frame MARGINWIDTH=0 MARGINHEIGHT=0 SCROLLING=\"Auto\" FRAMEBORDER=\"0\" FRAMESPACING=0 name=\"Fiche\"") ;
	    printf (" src=\"av.cgi?db=%s&c=gene&a=v%s&wh=%d&l=%s\">\n", ppp->dbName, ppp->faction, ppp->height, url_encode(ac_name (ppp->geneBox))) ;
	  }
	}
      else 
	{
	  printf (" <frame MARGINWIDTH=0 MARGINHEIGHT=0 SCROLLING=\"Auto\" FRAMEBORDER=\"0\" FRAMESPACING=0 name=\"Fiche\"") ;
	  printf (" src=\"av.cgi?db=%s&c=%s&a=v%s&wh=%d&l=%s\">\n", ppp->dbName, ppp->class, ppp->faction, ppp->height, url_encode(ac_name (ppp->mrna))) ;
	}
  
      /* GENE OR MRNA Graphic */
      
      
      printf ("<frame MARGINWIDTH=0 MARGINHEIGHT=0 SCROLLING=\"Auto\" FRAMEBORDER=\"0\" FRAMESPACING=0 NAME=\"Graph\"") ;  
      if (!mrnaMode)
	printf (" src=\"av.cgi?db=%s&a=vgene&c=gene&wh=%d&l=%s\">\n", ppp->dbName, ppp->height, url_encode (ac_name (ppp->geneBox))) ;
      else 
	printf (" src=\"av.cgi?db=%s&a=vmrna&c=mrna&wh=%d&v=2&l=%s\">\n", ppp->dbName, ppp->height, url_encode (ac_name (ppp->mrna))) ;

      printf ("</FRAMESET>\n ") ;
    }
  else
    {
      char *SG = "S" ; /* flash by default */
      if (!strcmp (ppp->version,"v47")) SG = "G" ;

      printf ("\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
	      "  <!--\n" 
	      "  document.write (writeInfo (1, \"%s\",\"%s%s\",\"%s\",%d,\"%s\") ) ;\n"
	      "  //-->\n</script>\n"
	      ,	ppp->dbName
	      , mrnaMode ? "vmrna&" : "vgene&" 
	      , SG
	      , mrnaMode ? "mrna" : "gene" 
	      , ppp->details
	      , mrnaMode ? url_encode (ac_name (ppp->mrna)) : url_encode (ac_name (ppp->geneBox))
	      ) ;
    }
}

/***************************************************************************/

static void showFramesInfo (PP *ppp)
{
  int mrnaMode = 0 ;
  char *SG = "S" ; /* flash by default */
  if (!strcmp (ppp->version,"v47")) SG = "G" ;
  if (!strcasecmp (ppp->class, "mRNA"))
    mrnaMode = 1 ;

  standardHtmlHeader (ppp, FALSE) ;

  if (!ppp->newWindow) 
    showInfoFrameSet (ppp) ;
  else if (strcmp (ppp->version,"v47") || (ppp->action != GENE_A && ppp->action != MRNA_A)) /* old system */
    {
      printf ("<FORM NAME=\"vsim2\">\n") ;
      printf ("<frame MARGINWIDTH=0 MARGINHEIGHT=0 SCROLLING=\"Auto\" FRAMEBORDER=\"0\" FRAMESPACING=0 NAME=\"Graph\"") ;
      if (!mrnaMode)
	printf (" src=\"av.cgi?db=%s&a=%s%s&wh=%d&c=%s&l=%s&v=%x\">\n", 
		ppp->dbName, "vgene&",SG,  ppp->height, "Gene", url_encode (ac_name (ppp->geneBox)), ppp->details) ;
      else 
	printf (" src=\"av.cgi?db=%s&a=%s%s&wh=%d&c=mrna&v=2&l=%s&v=%x\">\n", 
		ppp->dbName, "vmrna&",SG, ppp->height, url_encode (ac_name (ppp->mrna)), ppp->details) ;
      printf ("</FORM>\n") ;
    }
  else /* in this case write info should just show the drawing but not the fiche */
    {
       printf ("\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
	      "  <!--\n" 
	      "  document.write (writeInfo (0, \"%s\",\"%s%s\",\"%s\",%d,\"%s\") ) ;\n"
	      "  //-->\n</script>\n"
	      ,	ppp->dbName
	      , mrnaMode ? "vmrna&" : "vgene&" 
	       , "G"
	      , mrnaMode ? "mrna" : "gene" 
	      , ppp->details
	      , mrnaMode ? url_encode (ac_name (ppp->mrna)) : url_encode (ac_name (ppp->geneBox))
	      ) ;
   }
  

  printf ("</html>\n") ;
}

/*************** general layout of the triple aceView graph ****************/
/*** Browser will later call back aceview to fill the different frames ****/
  /******************************************************************************/

static void showFramesMain (PP *ppp)
{
  int mrnaMode = 0 ;
  BOOL showTabFrame = TRUE ; /* needed for the bubbles */
  BOOL showInfoFrame = FALSE ;
  if (!strcmp (ppp->version,"v47")) showInfoFrame = TRUE ;
  if (!strcasecmp (ppp->class, "mRNA"))
    mrnaMode = 1 ;
  ppp->noJava = TRUE ;
  standardHtmlHeader  (ppp, FALSE) ;

  printf ("<FORM NAME=\"vsim\">\n") ;
  if (!ppp->newWindow)
    {
      printf ("<FRAMESET COLS=\"150, *\">\n") ;

      printf ("<frame  scrolling=\"yes\" MARGINWIDTH=0 MARGINHEIGHT=0 SCROLLING=\"Auto\" NORESIZE FRAMEBORDER=\"0\" FRAMESPACING=0 NAME=\"NCBI_banner\" src=\"ncbi_banner_%s.%s.html\">\n", ppp->dbName, ppp->version ) ;
    }
  
  if (showTabFrame)
    {
      printf ("<FRAMESET ROWS=\"50, *\">\n") ;
      
      printf ("<frame  MARGINWIDTH=0 MARGINHEIGHT=0 SCROLLING=\"No\" NORESIZE FRAMEBORDER=\"0\" FRAMESPACING=0 NAME=\"tabs\" src=\"tabs.%s.html?%d\" >\n", ppp->version, mrnaMode ? 2 : 1) ;
    }
  /* Info */
  if (showInfoFrame)
    { 
      printf ("<frame  MARGINWIDTH=0 MARGINHEIGHT=0 SCROLLING=\"Auto\" NORESIZE FRAMEBORDER=\"0\" FRAMESPACING=0 NAME=\"info\"") ;
      if (!mrnaMode)printf (" src=\"av.cgi?db=%s&c=gene&a=info&l=%s\">\n", ppp->dbName, url_encode (ac_name (ppp->geneBox))) ;
      else printf (" src=\"av.cgi?db=%s&c=%s&a=info&l=%s\">\n", ppp->dbName, ppp->class, url_encode (ac_name (ppp->mrna))) ;
    }
  else
    {
      printf ("<frame  MARGINWIDTH=0 MARGINHEIGHT=0 SCROLLING=\"Auto\" NORESIZE FRAMEBORDER=\"0\" FRAMESPACING=0 NAME=\"info\"") ;
      if (!mrnaMode)printf (" src=\"av.cgi?db=%s&c=gene&a=fiche&l=%s\">\n", ppp->dbName, url_encode (ac_name (ppp->geneBox))) ;
      else printf (" src=\"av.cgi?db=%s&c=%s&a=fiche&l=%s\">\n", ppp->dbName, ppp->class, url_encode (ac_name (ppp->mrna))) ;
    }
  if (showTabFrame) 
    printf ("</FRAMESET>\n") ;
  if (!ppp->newWindow)
    printf ("</FRAMESET>\n") ;

  printf ("<NOFRAMES>\n") ;
  printf ("<body>\n") ;
  printf ("AceView offers a comprehensive annotation of human, mouse and nematode genes\n"
	  " reconstructed by co-alignment and clustering of all publicly available\n"
	  " mRNAs and ESTs on the genome sequence. Our goals are to offer a reliable\n"
	  " up-to-date resource on the genes, their functions, alternative variants,\n"
	  " expression, regulation and interactions, in the hope to stimulate\n"
	  " further validating experiments at the bench\n"
 	  "<p>\n"
	  ) ;
  printf ("This site normally  requires a frame enabled viewer,\n"
	  ) ;
  if (!strcasecmp (ac_class (ppp->geneBox), "gene"))
    printf (" However, click here to access our <a href=\"av.cgi?db=%s&F=%s\">full description<a> of gene %s.</a>\n"
	    , ppp->dbName, url_encode (ac_name (ppp->geneBox))
	    , ac_name (ppp->geneBox)
	    ) ;
  printf (" , go to our <a href=\"help.html\">main page</a> or see our <href=\"AceViewNews.html\">latest news</a>\n"
	  ) ;
  if (ppp->isGoogle && ppp->geneBox)
    {
      unsigned char *fiche2 = ac_command (ppp->db
				 , messprintf("view -v FEXP -c %s -n \"%s\" -xml -preview2",  
					      ppp->class,  ac_name(ppp->geneBox))
				 , 0, 0) ;
      if (fiche2)
	{
	  printf ("%s\n", fiche2) ;  /* the actual fiche */
	  messfree (fiche2) ;
	}
    }
  printf ("</body>\n") ;
  printf ("</NOFRAMES>\n") ;
  printf ("</FRAMESET>\n ") ;

  printf ("</FORM>\n</html>\n") ;
}

/***************************************************************************/
/* Produces either a .swf file
   This is achieved using the gif commands available on the gifaceserver.
*/
/*******************/
static int exportLocatorPicture (PP *ppp)
{ 
  if (ppp->picCommand && !strcmp(ppp->picCommand, " swfdump -compile -"))
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " view "
						   " -c Gene -n %s " 
						   " -v DtGLOC -flash -preview2 "
						   , ac_protect (ac_name(ppp->geneBox), 0)
						   )
				     , &(ppp->pictureBufferSize) , 0) ;
  else if (ppp->picCommand && !strcmp(ppp->picCommand, " svg "))
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " view "
						   " -c Gene -n %s " 
						   " -v DtGLOC -svg -preview2 "
						   , ac_protect (ac_name(ppp->geneBox), 0)
						   )
				     , &(ppp->pictureBufferSize) , 0) ;
  else if (ppp->picCommand)
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " gif ;"
						   " dimensions %d %d ;"
						   " display -D DtGLOC Gene %s ;"
						   " %s "
						   , (int) (1.2 * ppp->width + .9)
						   , (int) (1.2 * ppp->height + .9) 
						   , ac_name(ppp->geneBox)
						   , ppp->picCommand)
				     ,&(ppp->pictureBufferSize) , 0) ;
  
   return 1 ;
}

/*******************/
static int exportLocatorBigPicture (PP *ppp)
{ 
  if (ppp->picCommand && !strcmp(ppp->picCommand, " swfdump -compile -"))
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " view "
						   " -c Gene -n %s " 
						   " -v DtGLOCBIG -flash -preview2 "
						   , ac_protect (ac_name(ppp->geneBox), 0)
						   )
				     , &(ppp->pictureBufferSize) , 0) ;
  else if (ppp->picCommand && !strcmp(ppp->picCommand, " svg "))
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " view "
						   " -c Gene -n %s " 
						   " -v DtGLOCBIG -svg -preview2 "
						   , ac_protect (ac_name(ppp->geneBox), 0)
						   )
				     , &(ppp->pictureBufferSize) , 0) ;
  else if (ppp->picCommand)
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " gif ;"
						   " dimensions %d %d ;"
						   " display -D DtGLOCBIG Gene %s ;"
						   " %s "
						   , (int) (1.2 * ppp->width + .9)
						   , (int) (1.2 * ppp->height + .9) 
						   , ac_name(ppp->geneBox)
						   , ppp->picCommand)
				     ,&(ppp->pictureBufferSize) , 0) ;
  
   return 1 ;
}

/*******************/

static int exportHmrnaPicture (PP *ppp)
{ 
  if (ppp->picCommand && !strcmp(ppp->picCommand, " swfdump -compile -"))
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " view "
						   " -c Gene -n %s " 
						   " -v DtHseq -flash -preview2 "
						   , ac_protect (ac_name(ppp->geneBox), 0)
						   )
				     , &(ppp->pictureBufferSize) , 0) ;
  else if (ppp->picCommand && !strcmp(ppp->picCommand, " svg"))
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " view "
						   " -c Gene -n %s " 
						   " -v DtHseq -svg -preview2 "
						   , ac_protect (ac_name(ppp->geneBox), 0)
						   )
				     , &(ppp->pictureBufferSize) , 0) ;
  else if (ppp->picCommand)
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " gif ;"
						   " dimensions %d %d ;"
						   " display -D DtHseq Gene %s ;"
						   " %s ",
						   (int) (1.2 * ppp->width + .9),
						   (int) (1.2 * ppp->height + .9), 
						   ac_name(ppp->geneBox),
						   ppp->picCommand)
				     ,&(ppp->pictureBufferSize) , 0) ;
  
   return 1 ;
}

/*******************/

static int exportVGeneFlashPicture (PP *ppp)
{
  if (ppp->picCommand && !strcmp(ppp->picCommand, " swfdump -compile -"))
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " view "
						   " -c Gene -n %s " 
						   " -v DtVGene -flash -preview2 "
						   , ac_protect (ac_name(ppp->geneBox), 0)
						   )
				     , &(ppp->pictureBufferSize) , 0) ;
  else if (ppp->picCommand && !strcmp(ppp->picCommand, " svg"))
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " view "
						   " -c Gene -n %s " 
						   " -v DtVGene -svg -preview2 "
						   , ac_protect (ac_name(ppp->geneBox), 0)
						   )
				     , &(ppp->pictureBufferSize) , 0) ;
 
   return 1 ;
}

/*******************/

static int exportVGeneSvgPicture (PP *ppp)
{
  if (ppp->picCommand && !strcmp(ppp->picCommand, " svg"))
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " view "
						   " -c Gene -n %s " 
						   " -v DtVGene -svg -preview2 "
						   , ac_protect (ac_name(ppp->geneBox), 0)
						   )
				     , &(ppp->pictureBufferSize) , 0) ;
 
   return 1 ;
}

/*******************/

static int exportVmrnaSvgPicture (PP *ppp)
{
  if (ppp->picCommand && !strcmp(ppp->picCommand, " svg"))
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " view "
						   " -c mRNA -n %s " 
						   " -svg -preview2 "
						   , ac_protect (ac_name(ppp->geneBox), 0)
						   )
				     , &(ppp->pictureBufferSize) , 0) ;
 
   return 1 ;
}

/*******************/

static int exportVGeneGifPicture (PP *ppp)
{
  int  x1 = 0, x2 = 0, dx, origin = 0, sign ;
  char  map[255] = "" ; /* char *detailOn = "on" ; */
  AC_TABLE tt = 0 ;
  
  if (ppp->details == 0 &&
      (
       ac_has_tag (ppp->geneBox, "cloud_gene")
       )
      )
    ppp->details = 2 ;

  if (ppp->seqlen>1000000) ppp->seqlen = 1000000 ;

  if (!ac_has_tag (ppp->geneBox, "SMAP"))
    {
      printf ("<br/>\n<br/>\nSorry, this gene has not yet been mapped.<br/>\n") ;
      printf ("If you know its molecular identification, ") ;
      printf ("<a href=\"mailto:mieg@ncbi.nlm.nih.gov\"><b>please</b><a>") ;
      printf (" let us know.") ;
      return 0 ;
    }
  /* get the map coordinates */
  tt = ac_tag_table (ppp->geneBox, "IntMap", 0) ;
  if (tt && tt->cols > 2)
    {
      strcpy (map, ac_table_printable (tt, 0, 0, "")) ;
      x1 = origin = ac_table_int (tt, 0, 1, 0) ;
      x2 = ac_table_int (tt, 0, 2, 0) ;
      sign = x1>x2 ? -1 : 1 ;
      /*   xmid = (x1 + x2)/2 ; offset = sign * (x2 - x1)/10 ; */
      dx = x1 < x2 ? x2 - x1 : x1 - x2 ;
      dx /= 10 ;
      if (dx < 500) dx = 500 ;
      if (dx > 8000) dx = 8000 ;
      if (x1 < x2) { x1 -= dx ; x2 += dx ; }
      else { x1 += dx ; x2 -= dx ; }

      if (ppp->seqstart) x1 += ppp->seqstart*sign ;
      if (ppp->seqlen) x2 = x1 + ppp->seqlen*sign ;
      ppp->seqlen = sign * (x2-x1) ;
      if (ppp->seqlen < 10)
	{
	  ppp->seqlen = 10 ;
	  x2 = x1 + ppp->seqlen ;
	}
    }
  ac_free (tt) ;
  /*
  if ( (ppp->dbSpecies == 'h' && ppp->seqlen >= 500000) ||
       (ppp->dbSpecies == 'w' && ppp->seqlen >= 100000) )
    detailOn = "off" ;
  
    width  += 10 * acCount (ppp->geneBox, "Follow Mrna") ;
    i = acCount (ppp->geneBox, "Follow Antisens_to ; Follow Mrna") ;
    if (i)
    width  += 10 * i + 50 ;
    if (!strcmp (cdnaOn, "on"))

    {
    int cnt = acCount (ppp->geneBox, "Follow cDNA_clone") ;
    width  += 13* (cnt < 60 ? cnt : 60) ;
    }
  */

  if (ppp->picCommand)
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   " gif ;"
						   " dimensions %d %d ;"
						   " seqget \"%s\" -coords %d %d -origin %d -view %s ;"
						   " seqdisplay -view %s ;"
						   " %s ",
						   ppp->details & show_cDNA ? 1100 : 700,
						   ppp->height, map, x1, x2, origin,
						   ppp->details & show_cDNA ? "av_tg_whole" : "av_tg" , /* seqget -view */
						   ppp->details & show_cDNA ? "av_tg_whole" : "av_tg" , /* seqdisplay -view */
						   ppp->picCommand)
				     ,&(ppp->pictureBufferSize) , 0) ;  
  return 1 ;
}

/***************************************************************************/
/* Same idea as exportGene but for mRNAs. Somewhat simpler at the moment as there is currently
   no option for zoom (view has no effect on this function).
*/

static int exportVmrnaPicture (PP *ppp)
{
  if (ppp->details == 0)
    ppp->details = 2 ;
	 
  /* width  += 10 * acCount (ppp->mrna, "Follow cDNA_clone") ; */
  if (ppp->picCommand)
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   "gif ;"
						   "dimensions %d %d ;"
						   " seqget -class mrna \"%s\" -view %s ;"
						   " seqdisplay -view %s;"
						   " %s"
						   , 20000  /* ppp->details&show_cDNA ? 1200 : 600 */
						   , ppp->height, ac_name (ppp->mrna)
						   , ppp->details&show_cDNA ? "av_mrna_whole" : "av_mrna"  /* seqget -view */
						   , ppp->details&show_cDNA ? "av_mrna_whole" : "av_mrna"  /* seqdisplay -view */
						   , ppp->picCommand)
				     ,&(ppp->pictureBufferSize) , 0) ;
  
  return 1 ;
}

/***************************************************************************/
/* Same idea as exportGene but for mRNAs. Somewhat simpler at the moment as there is currently
   no option for zoom (view has no effect on this function).
*/

static int exportProbePicture (PP *ppp)
{
  if (ppp->details == 0)
    ppp->details = 2 ;
  
  /* width  += 10 * acCount (ppp->mrna, "Follow cDNA_clone") ; */
  
  if (ppp->picCommand)
    ppp->pictureBuffer = ac_command (ppp->db
				     , messprintf (
						   "gif "
						   "dimensions %d %d ;"
						   " display -D DtHseq  %s %s"
						   " ; %s" 
						   , ppp->details&show_cDNA ? 1200 : 600
						   , ppp->height, ac_class (ppp->geneBox), ac_name (ppp->geneBox)
						   , ppp->picCommand
						   )
				     ,&(ppp->pictureBufferSize) , 0) ;
  
  return 1 ;
}

/***************************************************************************/
/***************************************************************************/
/* sets genes, mrnas lists and gene, mrna: the active top member
 * in mRNA case we do an exact search then go for a generic gene search
 * if there are mutiple genes, we shift to LIST_Action
 */

static int getOneKey (PP *ppp) /* used in TRRE_A case */
{
  int nkeys = 0 ;

  if (*ppp->class)
    {
      ppp->mrnas = ac_dbquery_keyset (ppp->db
				    , messprintf ("find %s IS \"%s\"", ppp->class, ppp->query)
				    , 0) ;

      nkeys = ac_keyset_count (ppp->mrnas) ;
      if (nkeys)
	{
	  AC_ITER ai = ac_keyset_iter (ppp->mrnas, 0, 0) ;
	  ppp->geneBox = ac_next_obj (ai) ;
	  ac_free (ai) ;
	}
    }
  return nkeys ;
}

/***************************************************************************/

static int getProducts (PP *ppp)
{
  int nproduct = 0, ngenes = 0 ;
  AC_KEYSET products = 0 ;

  if (!strcasecmp (ppp->class, "product"))
    {
      products = ac_dbquery_keyset (ppp->db
				    , messprintf ("find product IS \"%s\"", ppp->query) 
				    , 0) ;
      nproduct = ac_keyset_count (products) ;

      if (nproduct == 1)
        {
	  AC_ITER ai = ac_keyset_iter (products, 0, 0) ;
          ppp->product = ac_next_obj (ai) ;
	  ppp->mrnas = ac_objquery_keyset (ppp->product, ">mrna", 0) ;
          ppp->genes = ppp->mrnas ? ac_ksquery_keyset (ppp->mrnas
						       , "{> From_gene;>Gene} SETOR {>from_prediction;>Model_of_gene}"
						       , 0) : 0 ;
	  ngenes = ppp->genes ? ac_keyset_count (ppp->genes) : 0 ;
	  ac_free (ai) ;
          if (ppp->genes)
	    {
	      ai = ac_keyset_iter (ppp->genes, 0, 0) ;
	      ppp->geneBox = ac_next_obj (ai) ;
	      ac_free (ai) ;
	    }
          if (ppp->mrnas)
	    {
	      ai = ac_keyset_iter (ppp->mrnas, 0, 0) ;
	      ppp->mrna = ac_next_obj (ai) ;
	      ac_free (ai) ;
	      ac_free (ppp->mrnas) ;
	    }
	  ppp->mrnas = ppp->genes ? ac_ksquery_keyset (ppp->genes, ">mrna;>from_gene ;> mrna", 0) : 0 ;
	  if (ppp->product) /* fix the spelling */
	    strcpy (ppp->query, ac_name (ppp->product)) ;
        }
    }
  if (ngenes > 1)
    ppp->action = LIST_A ;
  ac_free (products) ;
  return ngenes ;
}

/***************************************************************************/

static int getMrnas (PP *ppp)
{
  int nmrna = 0, ngenes = 0 ;

  if (!strcasecmp (ppp->class, "mrna"))
    {
      ppp->mrnas = ac_dbquery_keyset (ppp->db
				    , messprintf ("find mrna IS \"%s\"", ppp->query) 
				    , 0) ;
      nmrna = ac_keyset_count (ppp->mrnas) ;

      if (nmrna)
        {
	  AC_ITER ai = ac_keyset_iter (ppp->mrnas, 0, 0) ;
          ppp->mrna = ac_next_obj (ai) ;
          ppp->genes = ac_objquery_keyset (ppp->mrna
					, "{> From_gene;>Gene} SETOR {>from_prediction;>Model_of_gene}"
					, 0) ;
	  ngenes = ac_keyset_count (ppp->genes) ;
	  ac_free (ai) ;
          if (ppp->genes)
	    {
	      ai = ac_keyset_iter (ppp->genes, 0, 0) ;
	      ppp->geneBox = ac_next_obj (ai) ;
	      ac_free (ai) ;
	    }
	  ai = ac_objquery_iter (ppp->mrna, ">product ; best_product", 0) ;
	  ppp->product = ai ? ac_next_obj (ai)  : 0 ;
	  ac_free (ai) ;
        }
    }
  if (ngenes > 1)
    ppp->action = LIST_A ; 
  return ngenes ;
}

/***************************************************************************/

static int getGenes (PP *ppp)
{
  int ngenes = 0 ;
  int pass = 0 ;
  char *ctx, cOld = 0 ;
 
  { /* spaces cannot be exported from netscape the are mapped to \+ */
    char *cp, *cq ;

    cp = cq = ppp->query ;
    while (*cp)
      {
        if (*cp == '\\' && * (cp+1) == '+')
          { *cq++ = ' ' ; cp  += 2 ; }
        else
          *cq++ = *cp++ ;
      }
    *cq = 0 ;
  }
  
  ctx = contextLoad (ppp) ? "-ctx" : "" ;
 lao:
  ppp->genes = ac_command_keyset (ppp->db
				, messprintf ("WebQuery %s %s \"%s\"",
					      ctx, *ppp->class ? ppp->class : "0", ppp->query)
				, ppp->context_ks, 0) ;
  if (pass == 5 && ppp->genes)
    { 
      strcpy (ppp->class, "Gene") ;
    }

  ngenes = ppp->genes ? ac_keyset_count (ppp->genes) : 0 ;
  if (ngenes)
    {
      AC_ITER ai = ac_keyset_iter (ppp->genes, 0, 0) ;
      ppp->geneBox = ac_next_obj (ai) ;
      if (! *ppp->class)
	strcpy (ppp->class, "Gene") ;
      ac_free (ai) ;
    }

  if ( strcasecmp (ppp->class, "geneid") &&
       strcasecmp (ppp->class, "unigene") &&
       (! ppp->geneBox ||
	(ppp->geneBox && !strcasecmp (ac_class (ppp->geneBox), "KeySet")) 
	)
       )
    { 
      switch (pass) /* try to find something less specific */
	{
	case 0:  /* worm connection from tatiana using wormbase locus-tag */
	  pass = 1 ; /* avoid recursions */
	  if (! ctx[0] && ppp->dbSpecies == 'w' && ! strcasecmp (ppp->class, "NewName"))
	    {
	      strcpy (ppp->class, "Predicted_gene") ;
	      goto lao ;
	    }
	  /* else fall thru */
	case 1:  /* worm connection from tatiana using wormbase locus-tag */
	  if (! ctx[0] && ppp->dbSpecies == 'w' && ! strcasecmp (ppp->class, "Predicted_gene"))
	    {
	      char *cp = ppp->query + strlen (ppp->query) - 1 ;
	      if (*cp >= 'a' && *cp <= 'z' &&
		  *(cp-1) >= '0' &&  *(cp-1) <= '9')
		{
		  cOld = *cp ;
		  *cp = 0 ; 
		  pass = 2 ;
		  goto lao ;
		}
	      else if (*cp >= '0' && *cp <= '9')
		{
		  *(cp+1) = 'a' ; 
		  pass = 3 ;
		  goto lao ;
		}
	    }
	  pass = 4 ; /* avoid recursions */
	  /* else fall thru */
	case 2:
	case 3:
	  if (pass == 2) /* may be we fell thru */
	    {
	      char *cp = ppp->query + strlen (ppp->query) ;
	      *cp = cOld ; /* reestablish */ 
	    }
	  else if (pass == 3) /* may be we fell thru */
	    {
	      char *cp = ppp->query + strlen (ppp->query) - 1 ;
	      *cp = 0 ; /* reestablish */ 
	    } 
	  pass = 4 ; /* avoid recursions */
	  /* else fall thru */
	case 4:
	  if (! ctx[0] && ! strncasecmp (ppp->query, "LOC", 3))
	    {
	      char *cp, *cq ;

	      cp = ppp->query ; cq = cp+3 ;
	      while ((*cp++ = *cq++)) ;
	      strcpy (ppp->class, "GeneId") ;
	      pass = 5 ;
	      goto lao ;
	    }
	  pass = 6 ; /* avoid recursions */ 
	  /* else fall thru */
	case 5:
	  if (pass == 5)  /* may be we fell thru */
	    {
	      char *cp, *cq ;

	      cq = ppp->query + strlen (ppp->query) ; cp = cq + 3 ;
	      while (cq >= ppp->query) *cp-- = *cq-- ;
	      cp = ppp->query ; cq = "LOC" ;
              while (*cq) *cp++ = *cq++ ; /* do not copy the zero */	      
	      pass = 6 ;
	    }
	  pass = 6 ; /* avoid recursions */  
	  /* else fall thru */
	case 6:  /* try a general query */
	  pass = 7 ; /* avoid recursions */
	  if (! ctx[0] && ppp->class[0])
	    {
	      ppp->class[0] = 0 ; 
	      if (ppp->action != CLONES_A) 
		ppp->action = ZERO_A ;
	      goto lao ;
	    }
	  /* else fall thru */
	case 7:
	  { /* spaces cannot be exported from netscape the are mapped to \+ */
	    char *cp = ppp->query + strlen (ppp->query) - 1 ;
	    while (*cp != '.' && cp > ppp->query + 2) cp-- ;
	    if (*cp == '.') /* try shortening complex names */
	      {  /* notice that we do not increase pass becuase since we shorten query
		    we cannot loop for ever */
		*cp = 0 ;
		ppp->class[0] = 0 ;
		if (ppp->action != CLONES_A) 
		  ppp->action = ZERO_A ;
		goto lao ;
	      }
	  }
	  pass = 8 ; /* avoid recursions */
	case 8:
	  { /* spaces cannot be exported from netscape the are mapped to \+ */
	    char *cp = strstr (ppp->query, "and") ;

	    if (cp && cp > ppp->query + 3) /* try shortening names GENE1andGENE2 */
	      { 
		*cp = 0 ;
		ppp->class[0] = 0 ;
		if (ppp->action != CLONES_A) 
		  ppp->action = ZERO_A ;
		goto lao ;
	      }
	  }
	  /* else fall thru */
	default: /* definitive failure */
	  ppp->error = ppp->geneBox ; 
	  ppp->genes = ppp->mrnas = 0 ; 
	  return 0 ;
	}
    }
  else if (ppp->geneBox && !strcasecmp (ac_class (ppp->geneBox), "cDNA_clone"))
    {
      if (1) ppp->action = ppp->action == GOLD_A ? GOLD_A : CLONES_A ;
      ppp->clone = ppp->geneBox ;
      ppp->geneBox = ac_tag_obj (ppp->clone, "From_gene", 0) ;
      strcpy (ppp->class, "cDNA_clone") ;
      strcpy (ppp->query, ac_name (ppp->clone)) ;
      ppp->mrna = ac_tag_obj (ppp->clone, "in_mrna", 0) ;
    }
  else if (ppp->geneBox && ppp->action == GOLD_A  && !strcasecmp (ac_class (ppp->geneBox), "sequence"))
    {
      ppp->clone = ppp->geneBox ;
      ppp->geneBox = ac_tag_obj (ppp->clone, "From_gene", 0) ;
      strcpy (ppp->class, "Sequence") ;
      strcpy (ppp->query, ac_name (ppp->clone)) ;
      ppp->mrna = 0 ;
    }
  else if (ngenes > 1)
    ppp->action = LIST_A ;
  else if (ppp->geneBox)
    {
      ppp->mrnas = 0 ;
      if (ppp->mrnas && ! strcasecmp (ppp->class, "mrna"))
	ppp->mrnas = ac_objquery_keyset (ppp->geneBox
					 , "{>Transcribed_gene ! shedded_from ;>mrna ; IS *.a* } SETOR {>genefinder;>predicted_mrna} SETELSE {>Transcribed_gene ;>mrna} "
					 , 0) ;
      if (! ppp->mrnas  || ! ac_keyset_count (ppp->mrnas))
	ppp->mrnas = ac_objquery_keyset (ppp->geneBox
					 , "{>Transcribed_gene ! shedded_from ;>mrna ; IS *.a* } SETOR {>genefinder;>predicted_mrna} SETELSE {>Transcribed_gene ;>mrna} "
					 , 0) ;
      if (ppp->mrnas)
	{
	  AC_ITER ai = ac_keyset_iter (ppp->mrnas, 0, 0) ;
	  ppp->mrna = ac_next_obj (ai) ;
	  ac_free (ai) ;
	  ai = ac_objquery_iter (ppp->mrna, ">product ; best_product", 0) ;
	  ppp->product = ai ? ac_next_obj (ai)  : 0 ;
	  ac_free (ai) ;
	}
      if (ppp->mrna && ! strcasecmp (ppp->class, "mrna"))
	strcpy (ppp->query, ac_name (ppp->mrna)) ;
      else if (ppp->action != TG_SUPPORT_A)
	strcpy (ppp->query, ac_name (ppp->geneBox)) ;
    }
  return ngenes ;
}
  
/***************************************************************************/
/* Prints tree view of provided object as an html table to stdout */

static BOOL showTreeObject (PP *ppp)
{
#ifdef JUNK
  int col, row, lasttcol = -1, i, j, k, l ;
  char * bgcolor[3], * myName, * class ;
  AC_OBJ obj,  bob ;
  char * skiplist[] =
  { "NCBI", "ABI", "NewName", "Type", "Similarity", "Hit", "Gene", "" } ;
  char * classlist[] =
  { "DNA", "Protein", "LocusLink", "Locus", "OMIM", "cDNA_clone", "" } ;
  char * classURLs[] =
  { "https://www.ncbi.nlm.nih.gov/entrez/viewer.cgi?db=Nucleotide&dopt=GenBank&val=%s",
    "https://www.ncbi.nlm.nih.gov/entrez/viewer.cgi?db=Protein&dopt=GenPept&val=%s",
    "",
    "https://www.ncbi.nlm.nih.gov/LocusLink/LocRpt.cgi?l=%s",
    /***********************************************/
    /*                           "https://www.ncbi.nlm.nih.gov/entrez/disp*/
    /*                           omim.cgi?id=%s",       */
    /***********************************************/
    "https://www.ncbi.nlm.nih.gov/entrez/viewer.cgi?db=Nucleotide&dopt=GenBank&val=%s"
    } ;
  char * classHints[]=
  { "Entrez nucleotide record for dna %s",
    "Entrez protein record for peptide %s",
    "Look-up locus %s in LocusLink",
    "LocusLink report for locus %s",
    "OMIM phenotype record %s associated with this sequence",
    "Danielle owes Vahan a box of cachous"} ;
  char * oldtagnames[]=
  { "UTR_3prime", "UTR_5prime", "Partial_Exon",
    "Alternative_Partial_Exon", "Length_3prime_UTR", "Length_5prime_UTR", "" } ;
  char * newtagnames[]=
  { "3\' UTR", "5\' UTR", "Exon", "Alternative Exon", "3\' UTR length", "5\' UTR length" } ;
  char tagname[255] ;

  if (! (obj = ppp->geneBox))
    return FALSE ;

  bgcolor[0] = "eeeeee" ;
  bgcolor[1] = "ccccff" ;
  bgcolor[2] = "a0a0ff" ;

  /*  classURLs[0] = vstrDynaPrintf ( "/AceView/av.cgi?db=%s&a=fiche&c=dna&q=%%s", ppp->dbName) ; */
  classURLs[0] = strnew (messprintf ( "/AceView/av.cgi?db=%s&a=fiche&c=dna&q=%%s", ppp->dbName), 0) ;
  classURLs[2] = "https://www.ncbi.nlm.nih.gov/LocusLink/list.cgi?ORG=Hs&V=1&Q=%s" ;

  srand (1235423) ;

  col = obj->cols ;
  row = obj->rows ;

  class = ac_class (obj) ;

  printf ("<table width = 100%%><tr><td>") ;
  printf ("<h3>%s %s</h3>", class, ac_name (obj)) ;
  if (!strcmp (class, "Gene") || !strcmp (class, "Gene") || !strcmp (class, "MRNA"))
    {
      printf ("<td><a href=\"") ;
      printf ("av.cgi?db=%s&q=%s", ppp->dbName, ac_name (obj)) ;
      printf ("&c=%s\" target=AceView onmouseover=\"window.status=\'", class) ;
      printf ("AceView of %s", ac_name (obj)) ;
      printf ("\'\" onmouseout=\"window.status=\'\'\">AceView of %s</a></td>\n", ac_name (obj)) ;
    }

  printf ("\n<script language = \'JavaScript\' type=\'text/javascript\'>\n") ;
  printf ("  function zz () { ; }\n") ;
  printf ("    if (window.name == \"morewin\") document.write (\"<td><a href=\'javascript:zz ()\' onClick=\\\"window.name=\'popup_\' + parseInt (1000*Math.random ()) ; return false ;\\\" onmouseover=\\\"window.status=\'Preserve this window\' ; return true ;\\\">Preserve window</a></td>\\\n\")\n") ;
  printf ("    if (history.length > 1) document.write (\"<td><a href=\'javascript:zz ()' onClick=\\\"window.history.back () ; return false ;\\\" onmouseover=\\\"window.status=\'Back to previous page\'\\\">Back</a></td>\") ;\n") ;
  printf ("</script>\n\n") ;
  printf ("<td width=15 align=right><a href=\"javascript:zz ()\" onClick=\"window.close () ; return false ;\" onmouseover=\"window.status=\'Close popup\' ; return true ;\"><img src=\"/AceView/images/close.%s.gif\" alt=\"Close popup\" border=0></a>\n</td>", ppp->version) ;
  printf ("</tr></table>\n") ;

  printf ("<table cellpadding=5 cellspacing=0>\n") ;
  for (i = 0 ; i < row ; i++)
    {
      printf ("<tr>") ;
      for (j = 0 ; j<col ; j++)

        {
          bob = acGetObj (obj, i, j, 0) ;

          if (!bob)
            continue ;

          for (l = 0 ; skiplist[l][0] && strcmp (skiplist[l], ac_name (acGetObj (obj, i, 0, 0))) ; l++) ;
          if (skiplist[l][0])  break ;

          if (acType (bob) == 'g')
            {
              myName = ac_name (bob) ;
              for (l = 0 ; oldtagnames[l][0] && strcmp (oldtagnames[l], ac_name (acGetObj (obj, i, j, 0))) ;l++)
                { ; }
              if (oldtagnames[l][0]) myName = newtagnames[l] ;
              printf ("<td not_a_bgcolor = #%s>", j ? bgcolor[1] : bgcolor[0]) ;
              for (k = lasttcol ; k  >= j ; k--) printf ("</td></tr>\n</table>\n") ;
              strcpy (tagname, myName) ;
              for (k = 0 ; tagname[k] ; k++) if (tagname[k] == '_') tagname[k] = ' ' ;
              if (j == 0)

                {
                  printf ("\n<table width=\"100%%\" cellpadding=4 cellspacing=0><tr><td colspan=%d not_a_bgcolor=#a0a0ff><b>%s</b></td></tr>\n<tr>",
                          col-j, tagname) ;
                  lasttcol=j ;
                }
              else if (j == 1)
                {
                  printf ("\n<table cellpadding=2 cellspacing=0>\n<tr><td valign=top><b>%s</b></td>",
                          tagname) ;
                  lasttcol=j ;
                }
              else
                {
                  printf ("<b>%s</b></td>", tagname) ;
                }
            }
          else
            {
              if (j) printf ("<td valign=top not_a_bgcolor=#%2x%2x%2x>", 254 - j*128 / col, 254 - j*128 / col, 192 + i*64 / row) ;
              if (acType (bob) == 't' ||
                  (
                   (bob->typ == 'K' || bob->typ == 'k') && 
		   (
		    !strcmp (ac_class (bob), "Text") ||
		    !strcmp (ac_class (bob), "Map") ||
		    !strcmp (ac_class (bob), "Database")
		    )))
                printf ("%s", ac_name (bob)) ;
              else if (bob->typ == 'k' || acType (bob) == 'K')

                {
                  myName=ac_name (bob) ;
                  class=ac_class (bob) ;
                  for (l = 0 ; classlist[l][0] && strcmp (ac_class (bob), classlist[l]) ; l++ )
                    { ; }
                  if (classlist[l][0])
                    {
                      printf ("<a href = \"") ;
                      printf (classURLs[l], myName) ;
                      printf ("\" target=\"info\" onmouseover=\"window.status=\'") ;
                      printf (classHints[l], myName) ;
                      printf ("\'\" onmouseout=\"window.status=\'\'\">%s</a>", myName) ;
                    }
                  else
                    {
                      printf ("<a href=\"javascript:openAceViewLink ('%s', '%s', 't') ;\">", ac_class (bob), ac_name (bob)) ;
                      printf (ac_name (bob)) ;
                      printf ("</a>\n") ;
                    }
                }
              else if (bob->typ == 'f')
                {
                  printf ("%f", acFloat (bob)) ;
                }
              else if (bob->typ == 'i' || bob->typ == 'd')
                {
                  printf ("%d", acInt (bob)) ;
                }
              else
                {
                }
              if (j) printf ("</td>") ;
            }
        }
      printf ("</tr>\n") ;
    }
  for (k=lasttcol ; k  >= 0 ; k--) printf ("</table></td></tr>\n") ;

  printf ("</table>\n") ;
#endif
  return TRUE ;
} /* showTreeObject */

/***************************************************************************/
/* Prints DNA sequence */

static BOOL showDnaObject (PP *ppp)
{
  unsigned char *answer ;

  if (ppp->dna1)
    answer = ac_command (ppp->db
			 , messprintf("view -v dna -xml -c %s -n %s -remote %s -x1 %d -x2 %d -color %d "
				      , ppp->class, ac_name(ppp->geneBox)
				      ,  ppp->remote_addr
				      ,ppp->dna1, ppp->dna2, ppp->dnacolor)
			 ,0 ,0) ;
  else
    answer = ac_command (ppp->db
			 , messprintf("view -v dna -xml -c %s -n %s -remote %s "
				      , ppp->class, ac_name(ppp->geneBox)
				      , ppp->remote_addr)
			 ,0 ,0) ;
  printf ("%s", answer) ;
  messfree (answer) ;
  
  return TRUE ;
} /* showDnaObject */

/***************************************************************************/
/* Prints goldFile */

static BOOL showGoldFile (PP *ppp)
{
  unsigned char *answer ;

  answer = ac_command (ppp->db
		       , messprintf("view -gf %s", ppp->goldFile)
		       ,0 ,0) ;
  if (ppp->db)
    {
      if (answer && ! strstr ((char *)answer, "ERROR file not found"))
	  ac_command (ppp->db
		      , messprintf (" // goldFile OK %s", ppp->goldFile)
		      , 0, 0) ;
      else
	{
	  messfree (answer) ;
	  ac_command (ppp->db
		      , messprintf (" // goldFile ERROR file not found %s", ppp->goldFile)
		      , 0, 0) ;
	}
    }
  if (! answer)
    answer = (unsigned char *) strnew ("<html><body>Sorry, i do not understand your request, please send us a mail if the bug is reproducible</body></html>", 0) ;
  printf ("Content-type: text/html\n\n") ;
  printf ("%s", answer) ;
  messfree (answer) ;
  
  return TRUE ;
} /* showDnaObject */

/***************************************************************************/
/* Prints Peptide sequence */

static BOOL showPepObject (PP *ppp)
{
  unsigned char *answer ;

  if (ppp->dna1)
    answer = ac_command (ppp->db
			 , messprintf("view -v pep -xml -c %s -n %s -x1 %d -x2 %d  "
				      , ppp->class, ac_name(ppp->geneBox), ppp->dna1, ppp->dna2)
			 ,0 ,0) ;
  else
    answer = ac_command (ppp->db
			 , messprintf("view -v pep -xml -c %s -n %s "
				      , ppp->class, ac_name(ppp->geneBox))
			 ,0 ,0) ;

  printf ("%s", answer) ;
  messfree (answer) ;
  
  return TRUE ;
} /* showDnaObject */

/***************************************************************************/

/* Prints a row of an html table for the given gene object - the <table> and </table> must be output separately */
static BOOL listOneGeneTT (PP *ppp, AC_TABLE tt, int index, int dt, int offset)
{
  const char * cyto = "", * fullname = "", * bob, *goodName = "" ;
  char chrom[256], goodNameBuf[1024] ;
  AC_OBJ gene = ppp->geneBox ;
  char * callback = ppp->callback ;
  int nn, nclo ;
      
  gene = ac_table_obj (tt, index, dt, 0) ;

  goodName = ac_name(gene) ;

  strcpy (goodNameBuf, ac_name (gene)) ;
  strcpy (chrom, ac_table_printable (tt, index, 7 + dt, "")) ;

  if (!strncmp (goodName, "G_t", 3))
    {
      const char *ccp = goodName + 4 ;
      if (*ccp == '_')
        { strncpy (chrom, goodName + 3, 1) ; goodName  += 5 ; }
      if (* (ccp+1) == '_')
        { strncpy (chrom, goodName + 3, 2) ; goodName  += 6 ; }
    }
 
  fullname = ac_table_printable (tt, index, 8 + dt, "") ; /* gene->title */
  if (index % 2) bob = "<td class = GRB>" ;
  else bob = "<td>" ;

  if (ppp->dbSpecies == 'h')
    { 
      if (!strncasecmp (fullname, "Homo sapiens ", 13))
	fullname += 13 ;
      if (!strncasecmp (fullname, "gene ", 5))
	fullname += 5 ;
      nn = strlen (ac_name (gene)) ;
      if (!strncasecmp (fullname, ac_name (gene), nn))
	fullname += nn+1 ;
      if (!strncasecmp (fullname, " ", 1))
	fullname += 1 ;
      if (!strncasecmp (fullname, "encoding ", 9))
	fullname += 9 ;
      if (!strncasecmp (fullname, "a ", 2))
	fullname += 2 ;
      if (!strncasecmp (fullname, "an ", 3))
	fullname += 3 ;
      cyto = ac_table_printable (tt, index, 3 + dt, "") ; /* gene->Cytogenetic */
      nclo = ac_table_int (tt, index, 11 + dt, 0) ;       /* gene->tg->cdna-clone.COUNT */
      if (!nclo)
	nclo = ac_table_int (tt, index, 12 + dt, 0) ;     /* gene->has_cdna_clone.COUNT */
      printf ("<tr valign=\"top\">\n%s%d</td>"
	      "%s<a href=\"%s&c=%s&l=%s\" target=\"_top\"><b>%s</b></a></td>"
	      "%s%s</td>%s%s</td>%s%d</td>%s%s</td></tr>\n",
	      bob, index+offset,
	      bob, callback, ac_class (gene), ac_name (gene), goodName,
	      bob, chrom,
	      bob, cyto,
	      bob, nclo, 
	      bob, fullname) ;
    }
  else /* worm, ara, fly */
    {
      float z = -99999 ;
      const char *map = 0, *ess = "" ;

      if (!strncasecmp (fullname, "Homo sapiens ", 13))
	fullname += 13 ;
      if (!strncasecmp (fullname, "Arabidopsis thaliana ", 13))
	fullname += 22 ;
      if (!strncasecmp (fullname, "essential ", 10))
	{ fullname += 10 ; ess = "Essential, " ; }
      if (!strncasecmp (fullname, "gene ", 5))
	fullname += 5 ;
      nn = strlen (ac_name (gene)) ;
      if (!strncasecmp (fullname, ac_name (gene), nn))
	fullname += nn+1 ;
      if (!strncasecmp (fullname, " ", 1))
	fullname += 1 ;
      if (!strncasecmp (fullname, "encoding ", 9))
	fullname += 9 ;
      if (!strncasecmp (fullname, "a ", 2))
	fullname += 2 ;
      if (!strncasecmp (fullname, "an ", 3))
	fullname += 3 ;

      if (z == -99999)
	{
	  map =  ac_table_printable (tt, index, 1 + dt, "") ; /*  gene->InterpolatedMap */
	  z = ac_table_float (tt, index, 2 + dt, -99999) ;
	}
      if (z == -99999)
	{
	  map =  ac_table_printable (tt, index, 5 + dt, "") ; /*  gene->Map  */
	  z = ac_table_float (tt, index, 6 + dt, -99999) ;
	}
      if (z == -99999)
	{
	  map =  ac_table_printable (tt, index, 7 + dt, "") ; /*  IntMap */
	  z = ac_table_float (tt, index, 6 + dt, -99999) ;
	}
      /* ess = ac_table_printable (tt, index, 9 + dt, "") ;   gene->gene_id == essential || phenotype */
      nclo = ac_table_int (tt, index, 11 + dt, 0) ;       /* gene->tg->cdna-clone.COUNT */
      if (!nclo)
	nclo = ac_table_int (tt, index, 12 + dt, 0) ;     /* gene->has_cdna_clone.COUNT */
      if (z != -99999)
	{
	  printf ("<tr valign=\"top\">\n%s%d</td>"
		  "%s<a href=\"%s&c=%s&l=%s\" target=\"_top\"><b>%s</b></a></td>"
		  "%s%s;%s%.2f</td>%s%d</td>%s%s%s</td></tr>\n",
		  bob, index+offset,
		  bob, callback, ac_class (gene), ac_name (gene), goodName,
		  bob, map, z > 0 ? "+" : "", z, 
		  bob, nclo,
		  bob, ess, fullname
		  ) ;
	}
      else
	{
	  printf ("<tr valign=\"top\">\n%s%d</td>"
		  "%s<a href=\"%s&c=%s&l=%s\" target=\"_top\"><b>%s</b></a></td>"
		  "%s%s</td>%s%d</td>%s%s%s</td></tr>\n",
		  bob, index+offset,
		  bob, callback, ac_class (gene), ac_name (gene), goodName,
		  bob, chrom,
		  bob, nclo,
		  bob, ess, fullname
		  ) ;
	}
    }
  return TRUE ;
} /* listOneGeneTT */

/***************************************************************************/
/* same as above for an mrna object */

static BOOL listOneMrna (AC_OBJ mrna, int index, char * callback)
{
  char * bob ;
  char chrom[256] ;

  strcpy (chrom, ac_tag_printable (mrna, "IntMap", "")) ;
  if (index % 2) bob="<td class=GRB>" ;
  else bob="<td>" ;

  printf ("<tr valign=\"top\">\n%s%d</td>%s%s</td>%s%s</td></tr>",
          bob, index,
          bob, ac_name (mrna), bob, chrom) ;
  return TRUE ;
} /* listOneMrna */

/***************************************************************************/
/***************************************************************************/

static void showGeneListPaging (PP *ppp, int page, int cnt2)
{
  printf ("<table width=600><tr><td align=left>\n") ;
  printf ("Page %d of %d ", page + 1, cnt2) ;

  if (cnt2 > 0) printf (" go to ") ;
  if (page > 0)
    {
      printf ("<a href=\"%s&q=%s&p=0&c=%s\" onmouseover=\"window.status=\'First page\' ; return true ;\" onmouseout=\"window.status=\'\' ; return true ;\">  <bold> &nbsp; first </bold></a> ",
              ppp->callback, url_encode(ppp->query), ppp->class) ;
      printf ("<a href=\"%s&q=%s&p=%d&c=%s\"  onmouseover=\"window.status=\'Previous page\' ; return true ;\" onmouseout=\"window.status=\'\' ; return true ;\">  <bold> &nbsp; previous </bold></a> ",
              ppp->callback, url_encode(ppp->query), page - 1, ppp->class) ;
    }
  if (page < cnt2 - 1 )
    {
      printf ("<a href=\"%s&q=%s&p=%d&c=%s\"  onmouseover=\"window.status=\'Next page\' ; return true ;\" onmouseout=\"window.status=\'\' ; return true ;\">  <bold> &nbsp; next </bold></a> ",
              ppp->callback, url_encode(ppp->query), page + 1, ppp->class) ;
      printf ("<a href=\"%s&q=%s&p=%d&c=%s\" onmouseover=\"window.status=\'Last page\' ; return true ;\" onmouseout=\"window.status=\'\' ; return true ;\">  <bold> &nbsp; last </bold></a>\n",
              ppp->callback, url_encode(ppp->query), cnt2 -1, ppp->class) ;
    } 
  if (cnt2 > 0) printf (" page. ") ;
  printf ("</td></tr></table>\n") ;
}

/***************************************************************************/

static void showGeneList (PP *ppp)
{
  int i, j, cnt, cnt2 = 0 ;
  BOOL insideTable = FALSE ;
  int page = 1, pagelength = 100 ;
  const char *ccp ;
  AC_OBJ pfam ;
  AC_ITER pfamList = 0 ;
  char linkBuf[1024] ;

  /* cnt = ppp->tblObj ? ppp->tblObj->rows : (ppp->geneBoxs ? ppp->geneBoxs->objCnt : 0 ) ; */
  cnt = ppp->genes ? ac_keyset_count (ppp->genes) : 0 ;
  page = ppp->pagestart ;
  
  if (ppp->dbTitle && ppp->speciesTitle)
    printf ("On the %s annotation of the %s genome<br/>\n",
            ppp->dbTitle, ppp->speciesTitle) ;

  if (!strcasecmp (ppp->class, "Pfam") && 
      (pfamList = ac_query_iter (ppp->db, 0, messprintf( "Find Pfam IS \"%s\"", ppp->query), 0, 0)) &&
      (pfam = ac_next_obj (pfamList)))
    {
      const char * txtAccession ;
      printf ("<b>%d genes contain the Pfam motif<br/>\n", cnt) ;

      if (! (ccp = ac_tag_printable (pfam, "Definition", 0)))
	ccp = ppp->query ;

      if ( (txtAccession = ac_tag_printable (pfam, "Accession", 0)))
        {
	  char *cq = strstr (txtAccession , ".") ; /* drop pfam version number */
	  if (cq) *cq = 0 ;
	  sprintf (linkBuf, "http://www.sanger.ac.uk/cgi-bin/Pfam/getacc?%s", txtAccession) ; 
          printf ("<a href='%s' target='_top' >%s (see InterPro)</a>", linkBuf, ccp) ;
        }
      else printf ("%s", ccp) ;

      /*  the comment is much too long
          if ( (cp=ac_tag_printable (acOrigin (pfamList), "Comment", 0)))
          E
          printf ("<br/>\n%s\n", cp) ;
      */
      printf ("</b>\n") ;
      ac_free (pfamList) ;
    }
  else if (!strcasecmp (ppp->class, "BlastHit"))
    printf ("<b>%d genes contain a BlastP hit to the Protein<br/>\n<a href=\"https://www.ncbi.nlm.nih.gov/entrez/viewer.cgi?db=Protein&dopt=GenPept&val=%s\" target=\"_top\"><b>%s</b></a>\n", cnt, ppp->query, ppp->query) ;
  else if (!strcasecmp (ppp->class, "Psort"))
    printf ("<b>At least %d genes contain the Psort motif<br/>\n%s</b>\n", cnt, ppp->query) ;
  else if (!strcasecmp (ppp->class, "LocusLink"))
    {
      insideTable = TRUE ;
      if (cnt  >= 2500)
        printf ("<table width=600><td><h3>At least %d LocusLink genes relate to %s </h3></td>\n",
                cnt, ppp->query) ;
      else if (cnt > 1)
        printf ("<table width=600><td><h3>%d LocusLink genes relate to %s </h3></td>\n",
                cnt, ppp->query) ;
      else
        printf ("<table width=600><td><h3>%d LocusLink gene relates to %s </h3></td>\n",
                cnt, ppp->query) ;
    }
  else
    {
      insideTable = TRUE ;
      if (cnt  >= 2500)
        printf ("<table width=600><td><h3>At least %d genes relate directly or indirectly to %s </h3></td>\n",
                cnt, ppp->query) ;
      else 
        printf ("<table width=600><td><h3>%d genes relate directly or indirectly to %s </h3></td>\n",
                cnt, ppp->query) ;
    }
  if (cnt > pagelength)
    {
      cnt2 = cnt ;
      cnt2 = (int) (1 + (cnt2 - 1)/pagelength) ;
      showGeneListPaging (ppp, page, cnt2) ;
    }
  if (insideTable)
    printf ("</td></tr></table>") ;
  if (ppp->orderByNewName)
    {
      ppp->newNames = ac_command_keyset (ppp->db, "webquery gene2newname", ppp->genes, 0) ;
      ppp->nonMappedGenes = ac_keyset_count (ppp->genes) - ac_keyset_count (ppp->newNames) ;
      if (ppp->nonMappedGenes > 0)
	printf ("There are also <b>%d</b> uncloned genes visible only in alphabetic ordering"
		, ppp->nonMappedGenes) ;
    }
  if (ppp->dbSpecies == 'h')
    printf ("<table width=600 border=0 cellpadding=2 cellspacing=0 bgcolor=#CCCCFF>"
	    "<tr><td>"
	    "<table width=100%% border=0 cellpadding=6 cellspacing=0 bgcolor=#ffffff>"
	    "  <tr valign=\"top\" bgcolor=#CCCCFF>\n"
	    "    <td class=\"table-head\"></td>\n"
	    "    <td class=\"table-head\">Gene Name</td>\n"
	    "    <td class=\"table-head\">Aligned<br>on chrom</td>\n"
	    "    <td class=\"table-head\">Cyto location</td>\n"  
	    "    <td class=\"table-head\">Supporting<br>cDNA<br>clones</td>\n"
	    "    <td class=\"table-head\">Description</td>\n"
	    "  </tr>\n") ;
  else /* if (ppp->dbSpecies == 'w') */
    printf ("<table width=600 border=0 cellpadding=2 cellspacing=0 bgcolor=#CCCCFF>"
	    "<tr><td>"
	    "<table width=100%% border=0 cellpadding=6 cellspacing=0 bgcolor=#ffffff>"
	    "  <tr valign=\"top\" bgcolor=#CCCCFF>\n"
	    "    <td class=\"table-head\"></td>\n"
	    "    <td class=\"table-head\">Gene Name</td>\n"
	    "    <td class=\"table-head\">Map</td>\n"
	    "    <td class=\"table-head\">Clones</td>\n"
	    "    <td class=\"table-head\">Description</td>\n"
	    "  </tr>\n") ;


  /*  if (ppp->tblObj)
      for (i=page * pagelength, j=0 ; j < pagelength && i< (ppp->tblObj->rows) ;i++, j++)
      listOneGene (acGetObj (ppp->tblObj, i, 2, 0), i+1, ppp->callback) ;
  else
  */
  if (ppp->genes)
    {
      AC_TABLE tt ;
      int x0, dt = 0 ;
      AC_KEYSET subset = 0, newNames = 0 ;
      /* table maker */

      /* int gMax = ac_keyset_count (ppp->genes) ; pas encore evalue */
      x0 =  page * pagelength + 1 ;
      
      if (ppp->orderByNewName)
	{
	  dt = 1 ; /* col=hidden new name */
	  subset = ac_subset_keyset (ppp->newNames, x0, pagelength, 0) ;

	  tt = ac_tablemaker_table (ppp->db, "web2genelist3", subset, ac_tablemaker_db, 0, 0, 0, 0) ;
	}
      else
	{
	  dt = 0 ; /* col=hidden new name */
	  subset = ac_subset_keyset (ppp->genes, x0, pagelength, 0) ;
	  tt = ac_tablemaker_table (ppp->db, "web2genelist4", subset, ac_tablemaker_db, 0, 0, 0, 0) ;
	}
	

      if (!tt)
	{
	  printf("\n\n\n\n\nCannot make table of genes - web2genelist3 table missing from db?\n\n\n\n\n");
	  exit(1);
	}
      if (0) printf ("<br/>\npage=%d  pagelength=%d x0=%d mxnewnames=%d tt->rows=%d<br/>\n",
	      page,        pagelength, x0, ac_keyset_count(newNames) , tt->rows) ;
      for (i = 0 ; i < pagelength && i < tt->rows; i++)
	{
	  listOneGeneTT (ppp, tt, i, dt, x0) ;
	}
      ac_command (ppp->db, messprintf (" // table done %d lines" , tt ? tt->rows : 0), 0, 0) ;

      ac_free (tt) ;  
      ac_free (newNames) ;
      ac_free (subset) ;
    }

  else if (ppp->mrnas)
    {
      AC_ITER ai = ac_keyset_iter (ppp->mrnas, 0, 0) ;
      AC_OBJ mrna ;
      
      mrna = ac_next_obj (ai) ;
      for (i = 0 ; mrna && i < page * pagelength ; i++)
	mrna = ac_next_obj (ai) ;
      for (i = page * pagelength, j = 0 ; j < pagelength && mrna; i++, j++)
	{
	  ppp->mrna = mrna ;
	  listOneMrna (ppp->mrna, i+1, ppp->callback) ;
	  mrna = ac_next_obj (ai) ;
	}
      ac_free (ai) ;  
    }
  printf ("</table></td></tr></table>\n") ;

  if (cnt2 > 0)
    showGeneListPaging (ppp, page, cnt2) ;

  /* create a special division that will be dynamically transfered to the tabs frame */
  printf ("<div id='999'>\n") ; 
  
  printf ("<table border=0><tr><td>Back to  <a href=\"index.html?%s\" target=\"_top\"><IMG src=\"images/aceview.%s.gif\" width=134 border=0></a> home page</td></tr></table>\n", ppp->dbName,ppp->version) ;

  printf ("</div>\n") ;
  printf ("\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
	  "  <!--\n"
	  "    aceViewOnLoad() ;\n "
	  "  //-->\n"
	  "</script>\n"
	  ) ;
  
  printf ("</body></html>\n") ;
} /* showGeneList */

/***************************************************************************/

static void showContextQuery (PP *ppp)
{
  if (0) printf ("<table border=0><tr><td><a href=\"index.html?%s\" target=\"_top\"><IMG src=\"images/aceview.%s.gif\" width=134 border=0> Back to home page</a></td></tr></table>\n", ppp->dbName,ppp->version) ;

  if (!contextSave (ppp, ppp->genes))
    return ;

  if (0 && ! strcasecmp (ppp->dbName, "debug"))
    printf ("<br/> db=%s context= %s<br/>\n", ppp->dbName, ppp->new_context_name) ;

  printf(
	 "<FORM method=get action='av.cgi'>\n"
	 );
  
 
  printf ("<INPUT type=hidden name=db value=%s>\n", ppp->dbName) ;
  printf ("<INPUT type=hidden name=ctx value=%s>\n", ppp->new_context_name) ;

  printf ("<a href='help.%s.html#context_search'>Refine</a>\n",  ppp->version) ;

  printf (" this search by adding criteria<br/>\n") ;

  printf(
	 /* "<br/>\nExample,  ada*, dehydrogen, nerv, sugar metabolism, NM_062321." */
	 /* 	 "<br/>\nExample,  ada*, sugar metabolism, NM_062321."  */
	 "<br/>\n\n"
	 );

  printf ("<font color=#007f0f><INPUT type=text name=q size=77 value=\"\" style=PDP></font>") ;
  printf ("<br/>") ;
  if (0)
    {
      printf ("Sort: <font color=RED> <SELECT name=\"N\" onchange=submit >\n"
	      "<OPTION value=0 %s>alphabetically</OPTION>\n"
	      "<OPTION value=1 %s >by chromosome position</OPTION>\n"
	      "</SELECT>\n</font><br/>\n\n"
	      , ppp->orderByNewName ? "" : "selected"
	      , ppp->orderByNewName ? "selected" : ""
	      ) ;
    }
  else
    {
      printf ("Sort genes: <b><font color=#003f0f><INPUT TYPE=RADIO NAME=N VALUE=0 %s>alphabetically\n"
	      , ppp->orderByNewName ? "" : "CHECKED") ;
      printf (" <INPUT TYPE=RADIO NAME=N VALUE=1 %s>by chromosome position</font></b>\n"
	      , ppp->orderByNewName ? "CHECKED" : "") ;
    }

  printf ("<INPUT type=submit value=Go style=PDP>") ; 
  printf ("<p>") ;
  printf("</FORM>\n");
} /* showContextQuery */

/***************************************************************************/

static void exportJavaMapCallback (void)
{
  printf ("\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
	  "  <!--\n"
          /*      "var d = document ;\n"*/
          /*      "if (classe == \"GENE\") writeGeneHeader (d, \"Map\") ;\n"*/
          "  document.write (writeMap (classe) ) ;\n"
	  "  aceViewOnLoad () ;\n"
          /*      "d.close () ;\n"               */
          "  //-->\n</script>\n") ;
} /* exportJavaMapCallback */

/***************************************************************************/

static void logClientInfo (PP *ppp)
{
  if (ppp->db)
    {
      ac_command (ppp->db
		 , messprintf (" // REMOTE_ADDR=%s\t%s\t%s", ppp->remote_addr, ppp->date, ppp->args)
		 , 0, 0) ;
      ac_command (ppp->db
		 , messprintf (" // HTTP_HOST=%s\t%s\tbrowser=%s", ppp->http_host, ppp->args,ppp->browser ? ppp->browser : "NULL-browser")
		 , 0, 0) ;
    }
}

/***************************************************************************/
/***************************************************************************/
/*
 * a: ACTION  :: gif tree fiche clones fasta genes mrna mega list link info dna pep
 * c: CLASS :: (optional) class_name (i.e. Gene, cDNA_clone, ...)
 * db: DATABASE :: human, worm, 30 ... or taxon number, mapped to human,worm... if registered
 * F: GENE_FICHE ::unc-32, implies class = gene & a = fiche
 * f: FACTION :: fiche clones blastp
 * l:   alias of q
 * n: NEWWINDOW    no-parameter-value
 * o: origin :: n (nextdb), u (ucsc), g (genecard), b (ncbi), e (ebi)
 * p: page_start :: used in gene lists
 * N: page_start :: orderByNewName used in gene lists
 * q: QUERY :: query to ask inside this class
 * G/V/B/S/s: picture_style :: G:gif file, V:SVG, b:bubble support of the gif file, S:flash/swfc kitchen, s:actual flash/swfc movie
 * t: SHIFT :: int = start_of_zone
 * v: details :: i.e. cdna_decorate
 * w: ZOOM  :: int = extent_of_zone
 * x: anchor :: #anchor hopefully present in the fiche
 * h: host :: int = previous offset in db_config
 * ww: width :: int = width of the drawing in pixcells (600)
 * wh: height :: int = height of the drawing in pixcells (1200)
 */

static BOOL getParams (PP *ppp)
{
  char buff[1024], *args = 0, *ctx ;
  int i = 0 ;
  BOOL isExtended = FALSE ;
  BOOL force36a = FALSE ;
  BOOL forcemm37 = FALSE ;
  char timeBuf[25] ;

  strcpy (ppp->date, timeShow(timeNow(), timeBuf, 25)) ;
  ppp->http_host = getenv("HTTP_HOST");
  if (!ppp->http_host)
    ppp->http_host = "unknown" ;
  else
    ppp->http_host = url_sanitize (ppp->http_host) ;

  if (! ppp->remote_addr)
    ppp->remote_addr = getenv("HTTP_X_FWD_IP_ADDR") ; /* ("REMOTE_ADDR") ; */
  if (! ppp->remote_addr)
    ppp->remote_addr = getenv("REMOTE_ADDR");
  if (! ppp->remote_addr)
    ppp->remote_addr = "REMOTE_ADDR=unknown" ;
  else
    ppp->remote_addr = url_sanitize (ppp->remote_addr) ;

  ppp->action = ZERO_A ;
  ppp->picture_type = ZERO_T ;

  ppp->browser = getenv ("HTTP_USER_AGENT") ;
  if (0 && !ppp->browser)
    ppp->browser = "Mozilla/4.0" ;
  if (!ppp->browser)
    ppp->browser = "";
  else
    ppp->browser =  url_sanitize (ppp->browser) ;

  if (strstr (ppp->browser, "Googlebot") || strstr (ppp->browser, "yahoo.com"))
    ppp->isGoogle = TRUE ;
  args = getenv("QUERY_STRING");
  if (!args)
    args="-args unknown: POST?-";
  ppp->args = strnew (args, 0) ;
  
  ctx = cgi_clean_arg("h") ;  
  if (ctx)
    {
      char *s = url_sanitize (ctx) ;
      sscanf (s, "%d", &(ppp->host_offset)) ;
    }
  ctx = cgi_clean_arg("ctx") ;
  if (ctx)
    {
      char *s = url_sanitize (ctx) ;
      strcpy (ppp->context_name, s) ;
    }

  ctx = cgi_clean_arg("ww") ; ppp->width = 600*1.2 ;
  if (ctx)
    {
      char *s = url_sanitize (ctx) ;
      sscanf (s, "%d", &(ppp->width)) ;
    }
  ctx = cgi_clean_arg("wh") ; ppp->height = 1000*1.2 ;
  if (ctx)
    {
      char *s = url_sanitize (ctx) ;
      sscanf (s, "%d", &(ppp->height)) ;
    }
  ppp->width /= 1.2 ;
  ppp->height /= 1.2 ;  /* allow for the margins */

  if (cgi_arg_exists ("G"))
    ppp->picture_type = GIF_T ;
  else if (cgi_arg_exists ("B"))
    ppp->picture_type = BUBBLE_T ;
  else if (cgi_arg_exists ("S"))
    ppp->picture_type = SWF_T ;
  else if (cgi_arg_exists ("V"))
    ppp->picture_type = SVG_T ;
  else if (cgi_arg_exists ("s"))
    ppp->picture_type = swf_T ;
  
  { /* class */
    char *s ;
    if ((s = cgi_clean_arg("term2")) && *s)
      {
	isExtended = TRUE ;
	strcpy(ppp->class,"extended");
	if (strlen(s) < sizeof(ppp->query))
	  strcpy(ppp->query, s);
	else
	  ppp->query[0] = 0;
      }
    else
      {
	s = cgi_clean_arg("c");
	if (!s)
	  s = "" ;
	if (strlen(s) < sizeof(ppp->class))
	  strcpy(ppp->class,s);
	else
	  ppp->class[0] = 0;
      }
    if (!strcasecmp (ppp->class, "pg"))
      strcpy (ppp->class, "predicted_gene") ;
    else if (!strcasecmp (ppp->class, "locusid"))
      strcpy (ppp->class, "geneid") ;
  }


  { /* ACTION */
    char *buff = cgi_clean_arg ("a") ;
    if (!buff) 
      buff = "" ;
    strncpy (ppp->actionName, buff, 63) ; ppp->actionName[63] = 0 ;
    if (!strcasecmp (buff, "vmrna"))
      ppp->action = VmRNA_P ;
    else if (!strcasecmp (buff, "mrna")) /* backward compatibility */
      {  ppp->action = VmRNA_P ; ppp->picture_type = BUBBLE_T ; }
    else if (!strcasecmp (buff, "locator"))
      ppp->action = Locator_P ;	
    else if (!strcasecmp (buff, "locatorbig"))
      ppp->action = LocatorBig_P ;	
    else if (!strcasecmp (buff, "hmrna"))
      ppp->action = HmRNA_P ;	
    else if (!strcasecmp (buff, "vgene"))
      ppp->action = GENE_P ;	
    else if (!strcasecmp (buff, "gene")) /* backward compatibility */
      { ppp->action = GENE_P ;	; ppp->picture_type = BUBBLE_T ; }
     else if (!strcasecmp (buff, "gifprobe"))
      ppp->action = PROBE_P ;
    else if (!strcasecmp (buff, "tree"))
      ppp->action = TREE_A ;
    else if (!strcasecmp (buff, "dna"))
      ppp->action = DNA_A ;
    else if (!strcasecmp (buff, "pep"))
      ppp->action = PEP_A ;
    else if (!strcasecmp (buff, "probe"))
      ppp->action = PROBE_A ;
    else if (!strcasecmp (buff, "fiche"))
      ppp->action = FICHE_A ;
    else if (!strcasecmp (buff, "vfiche")) /* backwards compatibility */
      ppp->action = FICHE_A ;
    else if (!strcasecmp (buff, "fgene"))
      ppp->action = FGENE_A ;
    else if (!strcasecmp (buff, "fmol"))
      ppp->action = FMOL_A ;
    else if (!strcasecmp (buff, "fexp"))
      ppp->action = FEXP_A ;
    else if (!strcasecmp (buff, "ffunc"))
      ppp->action = FFUNC_A ;
    else if (!strcasecmp (buff, "section"))
      ppp->action = SECTION_A ;
    else if (!strcasecmp (buff, "clones"))
      ppp->action = CLONES_A ;
    else if (!strcasecmp (buff, "fasta"))
      ppp->action = FASTA_A ;
    else if (!strcasecmp (buff, "gene"))
      ppp->action = GENE_A ;
    else if (!strcasecmp (buff, "mrna"))
      ppp->action = MRNA_A ;
    else if (!strcasecmp (buff, "mega"))
      ppp->action = MEGA_A ;
    else if (!strcasecmp (buff, "gold"))
      ppp->action = GOLD_A ;
    else if (!strcasecmp (buff, "l"))
      ppp->action = LOG_A ;
    else if (!strcasecmp (buff, "list"))
      ppp->action = LIST_A ;
    else if (!strcasecmp (buff, "link"))
      ppp->action = LINK_A ;
    else if (!strcasecmp (buff, "info"))
      ppp->action = INFO_A ;
    else if (!strcasecmp (buff, "Tg_support"))
      {
	char *cp = 0, *cq = 0 ;
	int x1, x2 ;
	ppp->action = TG_SUPPORT_A ; 
	cp = cgi_clean_arg("q");
	if (cp)
	  cq = strstr (cp, ":") ;
	if (cq) 
	  { 
	    *cq = 0 ; 
	    strncpy (ppp->query, cp, 999) ;
	    cp = cq + 1;  cq = strstr (cp, ":") ;
	    if (cq) 
	      { 
		*cq = 0 ; 
		if (sscanf (cp, "%d", &x1) != 1) x1 = 0 ;
		cp = cq + 1; 
		if (sscanf (cp, "%d", &x2) != 1) x2 = 0 ;
	      }
	  }		
	strcpy (ppp->class, "transcribed_gene") ; 
	if (x1 >0 && x2 > 0)
	  {
	    ppp->action = TG_SUPPORT_A ;
	    ppp->dna1 = x1 ; ppp->dna2 = x2 ; 
	  }
      }

  }

  { /* FICHECONTENT */
    char *buff = cgi_clean_arg("f") ;
    if (!buff) 
      buff = "" ;

    if (ppp->action == SECTION_A)
      strcpy (ppp->faction, buff) ;
    else if (!strcasecmp (buff, "clones"))
      strcpy (ppp->faction, "clones") ;
    else if (!strcasecmp (buff, "blastp"))
      strcpy (ppp->faction, "blastp") ;
    else
      {
        strcpy (ppp->faction, "fiche") ;
      }
  }

  { /* NEWWINDOW */
    char *buff = cgi_clean_arg("n") ;
    if (!buff) buff = "" ;
    if (!strcasecmp (buff, "new"))
      ppp->newWindow = 1 ;
  }

  { /* shift & zoom */
    char *buff = cgi_clean_arg("t") ;
    if (!buff)
      buff = "" ;
    if (buff[0])
      ppp->seqstart = atoi (buff) ;
    
    buff = cgi_clean_arg("w") ;
    if (!buff) buff = "" ;
    if (buff[0])
      ppp->seqlen = atoi (buff) ;
  }


  { /* details */
    char *buff = cgi_clean_arg("v") ;
    if (!buff) 
      buff = "" ;
    if (buff[0])
      ppp->details = 0 ;
    sscanf (buff, "%x", &ppp->details) ;
  }

  { /* anchor */
    char *buff = cgi_clean_arg("x") ;
    if (!buff)
      buff = "" ;
    if (buff[0])
      strncpy (ppp->anchor, buff, 127) ;
  }
  
  
  { /* PAGESTART */
    char *buff = cgi_clean_arg("p") ;
    if (!buff) 
      buff = "" ;

    if (buff[0])
      ppp->pagestart = atoi (buff) ;
  }

  { /* ORDERBYNEWNAME */
    char *buff = cgi_clean_arg("N") ;
    if (buff && buff[0] == '1') 
      ppp->orderByNewName = TRUE ;
  }

  { /* F geneFiche */
    char *buff = cgi_clean_arg("F") ;
    if (buff)
      {
	strcpy (ppp->class, "gene") ;
	strcpy (ppp->faction, "fiche") ;
	ppp->action = FICHE_A ;
	strcpy (ppp->query, buff) ;
      }
  }

  { /* google crawlers */
    if (ppp->isGoogle && !ppp->action)
      {
	if (!*ppp->class) strcpy (ppp->class, "gene") ;
	strcpy (ppp->faction, "fiche") ;
	ppp->action = FICHE_A ;
      }
  }

  { /* org */
    char *buff = cgi_clean_arg("org") ;

    *ppp->org = 0 ; 
    if (!buff) 
      buff = "" ;

    strcpy (ppp->org, buff) ;
  }


  if (ppp->action != TG_SUPPORT_A &&  ! isExtended && !ppp->query[0] )
  { /* query, q or l are equivalent */
    BOOL forceGene = FALSE, forceMrna = FALSE ;
    char *cp ;
    cp = cgi_clean_arg("term");
    if (!cp || !*cp)
      cp = cgi_clean_arg("q");
    if (!cp || ! *cp)
      cp = cgi_clean_arg("l");

    if (cp && (strlen(cp) < sizeof(ppp->query)))
      strcpy(ppp->query, cp);
    else
      ppp->query[0] = 0;

    /* remapping of unigene names */
    {
      int dummyInt = 0 ;
      char dummyc = 0 ;
      if (!strncasecmp (ppp->query, "Hs.", 3) &&
	  sscanf (ppp->query + 3,"%d%c", &dummyInt, &dummyc) == 1)
	strcpy (ppp->class, "Unigene") ;
    }

    /* remapping of names that have an explicit date */
    if ((i = strlen(ppp->query)) > 5)
      {
	cp = ppp->query + i - 5 ;
	if (!strcasecmp (cp, "Dec03"))
	  { force36a = TRUE ; forceGene = TRUE ; }
	if (!strcasecmp (cp, "Nov04"))
	  { force36a = TRUE ; forceGene = TRUE ; }
	if (!strcasecmp (cp, "Jun05"))
	  { force36a = TRUE ; forceGene = TRUE ; }
	if (!strcasecmp (cp, "Aug05"))
	  { force36a = TRUE ; forceMrna = TRUE ; }
	if (!strcasecmp (cp, "Apr07"))
	  { force36a = TRUE ; forceMrna = TRUE ; }
	if (!strcasecmp (cp, "Jun07"))
	  { forcemm37 = TRUE ; forceMrna = TRUE ; }
      }
    
    if (! ppp->class[0] &&  forceMrna)
      strcpy (ppp->class, "mrna") ;
    if (forceGene)
      {
	ppp->class[0] = 0 ; /* force a general query */
	i = strlen(ppp->query) ;
	cp = ppp->query + i - 1 ;	   
	while (*cp != '.' && cp > ppp->query + 2) cp-- ;
	if (*cp == '.') *cp = 0 ;
      }

    cp = ppp->query ;
    if (cp)
      {
        i = 0 ;
        cp-- ;
        while (*++cp)
	  switch (*cp)
            {
            case '0':case '2':case '4':case '6':case '8':
            case '1':case '3':case '5':case '7':case '9':
              i++ ;
              break ;
            case '*':
            case '?':
            case ' ':
            case '\t':
              break ;
            default:
              i++ ;
              break ;
            }
        if ( !i)
          { *ppp->class = 0 ; *ppp->query = 0 ; /* go down to default query */ }
      }
  }

  { /* goldFile kills all other actions */
    char *goldFile = cgi_clean_arg("gf");
    if (goldFile && strlen(goldFile) < 512)
      {
	ppp->action = GOLDFILE_A ;
	strcpy (ppp->goldFile, goldFile) ;
      }
  }

  { /* SELECT DATABASE */
    char cc ;
    char *s, *cp ;
    int j0, j1, jj, npp ;
    PP myPP[64] ;
    PP *mypp = 0 ;
    AC_DB db = 0 ;
    DBCF *dbcf ;

    memset (buff, 0, sizeof(buff)) ;
    /* db=9606 has priority over db=human etc */
    s = cgi_clean_arg("db");
    if (s)
      {
	int tx = atoi (s) ;
	char *txdb = 0 ;
	BOOL isNum = TRUE ;
	
	for (cp = s ; *cp ; cp++)
	  if (*cp < '0' || *cp > '9')
	    isNum = FALSE ;
	if (isNum)
	  {
	    switch (tx)
	      {
	      case 9606: txdb = "human" ; break ;
	      case 10090: txdb = "mouse" ; break ;
	      case 10116: txdb = "rat" ; break ;
	      case 6239: txdb = "worm" ; break ;
	      case 3702: txdb = "ara" ; break ;
	      case 0: break ;
	      default: 
		sprintf(ppp->query, "<p>Sorry, this taxon %d is not yet annotated in AceView." 
			" Please <a href=\"mailto:mieg@ncbi.nlm.nih.gov\">mail us</a> if you think it should.<br>\n" 
			, tx
			);
		fprintf(stderr,"%s unknown_taxon=%d\n", logPrefix(ppp, "db_open"), tx) ;
		return FALSE ;
	      }
	  }
	if (txdb) strcpy (buff, txdb) ;
	else  if (s && (strlen(s) < sizeof(buff)))
	  strcpy(buff,s);
      }
    
    if (!strcasecmp (buff, "hg16"))
      strcpy (buff, "human") ;
    else if (!strcasecmp (buff, "29"))
      strcpy (buff, "human") ;
    else if (!strcasecmp (buff, "30"))
      strcpy (buff, "human") ;
    else if (!strcasecmp (buff, "31"))
      strcpy (buff, "human") ;
    else if (!strcasecmp (buff, "32"))
      strcpy (buff, "human") ;
    else if (!strcasecmp (buff, "33"))
      strcpy (buff, "human") ;
    else if (!strcasecmp (buff, "34"))
      strcpy (buff, "human") ;
    else if (!strcasecmp (buff, "35g"))
      strcpy (buff, "35g") ;
    else if (!strncasecmp (buff, "35", 2))
      strcpy (buff, "human") ;
    else if (!strncasecmp (buff, "36", 2))
      strcpy (buff, "human") ;
    else if (!strcasecmp (buff, "hg17"))
      strcpy (buff, "human") ;
    else if (!strcasecmp (buff, "hg18"))
      strcpy (buff, "human") ;
    else if (!strncasecmp (buff, "hg", 2))
      strcpy (buff, "human") ;

    if (forcemm37 && !strcmp (buff, "mouse"))
      strcpy (buff, "mm_37") ;
    else if (buff[0] == '3' || !strcmp (buff, "human"))
      {
	if (force36a) strcpy (buff, "36a") ;
      }

    memset (myPP, 0, 64*sizeof(PP)) ;
    if (!*buff)
      {
	/*
	 * if there is no database specified, then this link probably
	 * came from locuslink - see if we can identify the organism
	 * and use the right database
	 */
	if (!strcmp (ppp->org, "6239"))
	  strcpy (buff, "worm") ;
	else if (!strcmp (ppp->org, "9606"))
	  strcpy (buff, "human") ;
	else if (!strcmp (ppp->org, "10090"))
	  strcpy (buff, "mouse") ;
	else
	  strcpy (buff, "human") ; /* default value */ 
      }
    /*
     * db = database name.
     *
     * oct 7 2005, i change db-config into a struct
     * defined at the top of this file
     */


    /* first we detect all the servers for this db and store their descriptors in myPP */
    npp = 0 ;
    
    for (dbcf = dbConfig ; dbcf->dbName ; dbcf++)
      {
	if (!strcmp (dbcf->dbName, buff))
	  {
	    /*
	     * ok, this line is for the database that we want - see
	     * what it says
	     */
	    mypp = myPP + npp ;
	    mypp->dbName = dbcf->dbName ;
	    mypp->dbSpecies = dbcf->dbSpecies ;
	    mypp->version = dbcf->version ;
	    mypp->host = dbcf->host ;
	    mypp->port = dbcf->port ;
	    mypp->speciesTitle = dbcf->speciesTitle ;
	    mypp->dbTitle = dbcf->dbTitle ;
	    
	    if (strcmp(mypp->host,"-") == 0)
	      {
		sprintf (ppp->query,"Sorry, this URL db=%s is discontinued, please try db=worm or db=human\n",buff);
		return FALSE;
	      }
	    npp++ ; /* success */
	  }
      }
    
    if (!npp)
      {
	sprintf(ppp->query,"Sorry, this URL is incorrect, you asked for db=%s, please try db=human or db=worm\n",buff);
	fprintf(stderr,"%s %s\n",logPrefix(ppp, "db_unknown"), ppp->dbName);
	return FALSE;
      }

    /* now we initialise our merry go round according to the machine ip of the client */
    if (ppp->host_offset) /* try same machine as before */
      j0 = (ppp->host_offset - 1) % npp ;
    else
      {
	cp = ppp->remote_addr ;
	cc = 0 ;
	while (*cp) cc ^= *cp++ ;
	j0 = cc % npp ;
      }
    /* now we try to open the corresponding database */
    /* we go round robin but start at a client dependant position */
    for (j1 = jj = j0 ; jj < j0 + npp ; jj++)
      {
	char buff[300] ;
	j1 = jj % npp ;
	mypp = myPP + j1 ;
		  
	sprintf(buff , "a:%s:%d -timeout 180 ", mypp->host, mypp->port) ;
	
	db = ac_open_db (buff, 0) ;
	if (db)
	  break ;
      }
    if (!db)
      {
	sprintf(ppp->query, "Sorry, the AceView %s server is currently unavailable. Please try again later"
		, !ppp->dbName ||  strstr(ppp->dbName, "null") ? "" : ppp->dbName
		);
	fprintf(stderr,"%s %s\n", logPrefix(ppp, "db_open"), ppp->dbName);
	return FALSE ;
      }

    /* success */
    ppp->host_offset = j1 + 1 ;
    ppp->db = db ;
    ppp->dbName = mypp->dbName ;
    ppp->dbSpecies = mypp->dbSpecies ;
    ppp->version = mypp->version ;
    ppp->host = mypp->host ;

    ppp->speciesTitle = mypp->speciesTitle ;
    ppp->dbTitle = mypp->dbTitle ;

    sprintf(ppp->callback, "av.cgi?db=%s", ppp->dbName);

    if (!strcmp (ppp->version,"v47") && 
	(! ppp->picture_type ||  ppp->picture_type == SWF_T)
	)
      {
	if (ppp->action == VmRNA_P || ppp->action == GENE_P)
	      ppp->picture_type = BUBBLE_T ;
      }
    return TRUE;
  }
} /* getParams  */

/***************************************************************************/
/*
* This creates the links page.  For each possible link, we look at
* whether it makes sense for this type of database.
*/

static char * colors[] = {"#ffffff", "#efefff"} ;
static int curCol = 0 ;

static void linkstart(char *url)
{
  curCol = ( curCol + 1 ) % (sizeof (colors)/sizeof (colors[0]));
  printf ("<tr><td bgcolor=%s><font face='%s'><a href='%s' target='_top'>",
	  colors[ curCol ], fonts, url);
}

static void linkmid()
{
  printf("</a></font></td><td bgcolor=%s>",colors[curCol]);
}

static void linkend()
{
  printf("</td></tr>\n");
}


static void oneLineLink(char *a, char *b, char *c) 
{
  linkstart(b);
  printf("%s",a);
  linkmid();
  printf("%s",c);
  linkend();
}

static void showLinkList (PP * ppp)
{
  char linkBuf[1024] ;
  const char *nnam ;
  AC_OBJ obj, obj2 ;
  
  printf ("<h3> About the gene %s.</h3>", ac_name(ppp->geneBox)) ;
  
  curCol = 0;
  printf ("<table border=1 width=100%%>\n") ;
  
  
  /*
   * dbSpecies is one of 
   *	human
   *	worm
   *
   * We construct a different set of links, depending on the species.
   * (We do not use ppp->species 
   */
  if (ppp->dbSpecies == 'h')
    { /****************** OMIM ************************/
      if ( (obj = ac_tag_obj (ppp->geneBox, "Gene", 0)) &&
	   (obj2 = ac_tag_obj (obj, "OMIM", 0)) &&
	   (nnam = ac_name (obj2)))
        {
	  sprintf (linkBuf, "https://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=%s", nnam) ;
	  oneLineLink ("OMIM", linkBuf, "provides a description of the gene and its associated diseases.");
        }
    }
  
  if ( (ppp->dbSpecies == 'h') || (ppp->dbSpecies == 'w'))
    { /****************** locus Link ************************/
      nnam = ac_tag_printable (ppp->geneBox, "LocusId", 0) ;
      if (!nnam && (obj = ac_tag_obj (ppp->geneBox, "Genefinder", 0)))
 	nnam = ac_tag_printable (obj, "LocusId", 0) ;
      if (nnam)
        {
	  sprintf (linkBuf, "https://www.ncbi.nlm.nih.gov/LocusLink/LocRpt.cgi?l=%s", nnam) ;
	  oneLineLink ("LocusLink", linkBuf, "at NCBI, provides a description of the gene as well as links.");
        }
    }
  if (0)
    printf ("<br><image src=\"https://www.ncbi.nlm.nih.gov/AceView/av.cgi?db=worm&q=unc-32&c=gene&a=swf\" width=400 height=400 border=2><br>\n") ;
  else
    printf ("<br><image src=\"https://www.ncbi.nlm.nih.gov/AceView/totowz.xml\" width=400 height=400 border=2><br>\n") ;
  
  
  if (ppp->dbSpecies == 'w')
    { /****************** refseq ************************/
      AC_TABLE tt = ac_tag_table (ppp->geneBox, "NM_Id", 0) ;
      int ii ;
      const char *nm, *seq ;

      for (ii = 0 ; ii < tt->rows ; ii++)
	{
	  nm = ac_table_printable (tt, ii, 0, 0) ;
	  seq = ac_table_printable (tt, ii, 1, 0) ;
	  if (nm && seq)
	    {
	      sprintf (linkBuf, "https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=nucleotide&term=%s", nm) ;
	      oneLineLink ("RefSeq", linkBuf, 
			   messprintf ("%s %s at NCBI, provides links and a short description derived from AceView.", nm, seq)) ;
	    }
	}
      ac_free (tt) ;
    }
  
  
  if (ppp->dbSpecies != 'w')
    { /****************** refseq ************************/
      AC_OBJ tg = ac_tag_obj (ppp->geneBox, "Transcribed_gene", 0) ;
      AC_TABLE tt = tg ? ac_tag_table (tg, "cDNA_clone", 0) : 0 ;
      int ii ;
      const char *nm ;

      for (ii = 0 ; ii < tt->rows ; ii++)
	{
	  nm = ac_table_printable (tt, ii, 0, 0) ;
	  if (strncasecmp (nm, "NM", 2))
	    continue ;
	  if (nm)
	    {
	      sprintf (linkBuf, "https://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=nucleotide&term=%s", nm) ;
	      oneLineLink ("RefSeq", linkBuf, 
			   messprintf ("%s RefSeq at NCBI, provides links and descriptions.", nm)) ;
	    }
	}
      ac_free (tt) ;
    }
  
  
  if (ppp->dbSpecies == 'h')
    { /************************** GeneCard *****************************/
      nnam = ac_tag_printable (ppp->geneBox, "LocusLink", 0) ;
      if (!nnam)
	nnam = ac_tag_printable (ppp->geneBox, "Genefinder", 0) ;
      if (nnam)
	{
	  sprintf (linkBuf,
		   "http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s&keywords=%s"
		   ,nnam, nnam) ;
	  oneLineLink ("GeneCards", linkBuf, " offers concise information about the functions of all human genes that have an approved symbol, as well as selected others.");
	}
    }
  
  if ( ppp->dbSpecies == 'w')
    {
      AC_TABLE tt = 0 ;
      AC_OBJ thing2;
      int x;
      
      /***************************************
       * two variants of wormbase
       */
      tt = ac_tag_table (ppp->geneBox,"genefinder", 0);
      if (tt)  /* wormbase graphic sequence */
	{
	  for (x=0; x < tt->rows; x++)
	    {
	      if (! strncmp (ac_table_printable (tt, x, 0, "x"), "hmm", 3))
		continue ;
	      thing2 = ac_table_obj (tt, x, 0, 0) ;
	      if (!thing2)
		continue;
	    
	      sprintf(linkBuf, "http://www.wormbase.org/db/get?class=Sequence;name=%s",
		      ac_name(thing2));
	      linkstart(linkBuf);
	      printf("WormBase");
	      linkmid();
	      printf("Sequence %s\n",ac_name(thing2));
	      linkend();
	    }
	}
      else
	{
	  tt = ac_tag_table (ppp->geneBox, "locus", 0);
	  if (tt)  /* wormbase text page on locus */
	    {
	      for (x=0; x<tt->rows; x++)
		{
		  thing2 = ac_table_obj (tt, x, 0, 0) ;
		  if (!thing2)
		    continue;
		
		  sprintf(linkBuf, "http://www.wormbase.org/db/get?class=Locus;name=%s",
			  ac_name(thing2));
		  linkstart(linkBuf);
		  printf("WormBase");
		  linkmid();
		  printf("Locus %s\n",ac_name(thing2));
		  linkend();
		}
	    }
	}
      ac_free (tt) ;
    }
  
  /*
    http://www.wormbase.org/db/get?class=Sequence;name=ZK637.8a
    http://www.wormbase.org/db/seq/gbrowse?source=wormbase;name=II:8000000..8001800
    You can also use cosmid coordinates:
    http://www.wormbase.org/db/seq/gbrowse?source=wormbase;name=DH11:1..1000
    
    or GenBank coordinates
    http://www.wormbase.org/db/seq/gbrowse?source=wormbase;name=Z49126:-1000..1000
    
    The last example shows that you can use negative coordinates.  However, 
    coordinates are relative to base 1 of the sequence.
    
    4) open the best page on a region of the genetic map
    http://www.wormbase.org/db/misc/epic?name=unc-9;class=Locus
    
    or say ask fro chromosome II from centimorgan 11 to 12.3
    http://www.wormbase.org/db/misc/epic?class=Map;name=II;map_start=11;map_stop=12.3
    
    help document with all the public links ?
    http://www.wormbase.org/about/linking.html
    
    This has got most of the information, but I just noticed that it is missing 
    information about linking to the genetic maps.  
  */
  
  
  if (ppp->dbSpecies == 'w')
    { 
      AC_TABLE rnai = 0, tt = 0 ;
      AC_OBJ thing, clone_group ;
      int number;
      const char *s;
      int x;
      int heidelberg = 0;
      
      /***************************************
       * worfdb link 
       */
      rnai = ac_tag_table (ppp->geneBox, "rnai", 0) ;
      if (rnai)
	{
	  const char *s;
	  int x;
	  for (x=0; x < rnai->rows; x++)
	    {
	      thing = ac_table_obj (rnai,x,0, 0);
	      if (! thing)
		continue;
	      s = ac_name(thing);
	      if ((s[0] == 'm') && (s[1] == 'v') && (s[2] == '_'))
		{
		  sprintf (linkBuf, "http://worfdb.dfci.harvard.edu/searchallwormorfs.pl?by=name&sid=%s",s+3);
		  linkstart(linkBuf);
		  printf("Worfdb");
		  linkmid();
		  printf("%s",s+3);
		  linkend();
		}
	      if (ac_has_tag (thing,"DIC_Movie_available"))
		heidelberg = 1;
	    }
	}
      
      /***************************************
       * nextdb
       */
      tt = ac_tag_table (ppp->geneBox, "clone_group", 0);
      if (tt)
	for (x=0; x < tt->rows; x++)
	  {
	    clone_group = ac_table_obj (tt, x, 0, 0) ;
	    if (! clone_group)
	      continue;
	    
	    s = ac_name (clone_group);
	    if ((s[0] == 'Y') && (s[1] == 'K') && (sscanf(s+2, "%d", &number) == 1))
	      {
		sprintf (linkBuf, "http://nematode.lab.nig.ac.jp/cgi-bin/db/ShowGeneInfo.sh?celk=CELK%05d",number);
		linkstart(linkBuf);
		printf("%s","NextDB");
		linkmid();
		printf("The Yuji Kohara Worm Transcriptome database, with in situ hybridisation data to all stages of development, information about yk clones, YK clusters, RNAi data and other useful tools") ;
		linkend();
	      }
	  }
      
      /***************************************
       * CBG movies
       */
      if (heidelberg)
	{
	  linkstart("http://worm-srv1.mpi-cbg.de/dbScreen");
	  printf("CBG");
	  linkmid();
	  printf("Movies of embryo development in RNA interference experiments are available from Cenix BioScience GmbH\n");
	  linkend();
	}
      ac_free (tt) ;
    }

  printf ("</table>\n") ;
  
  
  curCol = 0;
  printf ("<h3> About the genome area.</h3>") ;
  printf ("<table border=1 width=100%%>\n") ;
  
  
  if ( ppp->dbSpecies == 'w')
    {
      AC_TABLE oThing ;
      int ir ; float z ;

      /***************************************
       * two variants of wormbase
       */
      oThing = ac_tag_table (ppp->geneBox, "InterpolatedMap", 0);
      if (oThing)  /* wormbase graphic genetic map */
	{
	  for (ir=0; ir < oThing->rows; ir++)
	    {
	      if (!ac_table_obj (oThing, ir, 0, 0))
		continue;
	      z = ac_table_float (oThing, ir, 1, 0) ;
	      sprintf (linkBuf, "//www.wormbase.org/db/misc/epic?class=Map;name=%s;map_start=%g;map_stop=%g",
		       ac_table_printable (oThing, ir, 0, ""), z - .3, z + .3);
	      linkstart (linkBuf);
	      printf ("WormBase");
	      linkmid ();
	      printf ("Region of the genetic map") ;
	      linkend ();
	    }
	}
    }
  
  if ( ppp->dbSpecies == 'h')
    {/************************** UCSC *****************************/
      AC_TABLE tt = 0 ;  
      
      if ( (tt = ac_tag_table (ppp->geneBox, "IntMap", 0)) && tt->cols >= 3)
        {
	  char *ucsctarget = !strcmp (ppp->dbName, "34") ? "hg16" : "hg17" ;
	  nnam = ac_table_printable (tt, 0, 0, "") ;
	  if (strstr (nnam, "_"))
	    nnam = strstr (nnam, "_")+1 ;
	  sprintf (linkBuf,
		   "http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&acembly=full&position=chr%s:%d-%d"
		   , ucsctarget
		   , nnam, ac_table_int (tt, 0, 1, 0), ac_table_int (tt, 0, 2, 0)) ;
	  oneLineLink ("UCSC", linkBuf,"Human genome at Jim Kent's. This viewer offers a comparison of various types of gene constructions, including AceView. This is the fastest way to zoom out and get an overview of the AceView genes in the region.");
        }
      ac_free (tt) ;
    }
  
  
  if ( (ppp->dbSpecies == 'h') || (ppp->dbSpecies == 'w'))
    {/************************** MapView *****************************/
      AC_TABLE tt = 0 ;
      if (ppp->dbSpecies == 'h' && (tt =  ac_tag_table (ppp->geneBox, "IntMap", 0)) && tt->cols >= 3)
	{
	  nnam = ac_table_printable (tt, 0, 0, "") ;
	  if (!strncmp(nnam,"CHROMOSOME_",11)) nnam += 11 ;
	  sprintf (linkBuf,
		   "https://www.ncbi.nlm.nih.gov/mapview/maps.cgi?ORG=%s&CHR=%s&MAPS=cntg,scan,loc&BEG=%d&END=%d"
		   , "hum"
		   ,nnam
		   , ac_table_int (tt, 0, 1, 0), ac_table_int (tt, 0, 2, 0)
		   ) ;
	  oneLineLink ("MapView", linkBuf,"displays the NCBI view of the same genomic region.");
	  ac_free (tt) ;
	}
      if (ppp->dbSpecies == 'w' && (tt =  ac_tag_table (ppp->geneBox, "IntMap", 0)) && tt->cols >= 3)
	{
	  nnam = ac_table_printable (tt, 0, 0, "") ;
	  if (!strncmp(nnam,"CHROMOSOME_",11)) nnam += 11 ;
	  sprintf (linkBuf,
		   "https://www.ncbi.nlm.nih.gov/mapview/maps.cgi?ORG=%s&CHR=%s"
		   ,"celegans" 
		   ,nnam 
		   /*, ac_table_int (tt, 0, 1, 0), ac_table_int (tt, 0, 2, 0) */
		   ) ;
	  oneLineLink ("MapView", linkBuf,"displays the NCBI view of the same genomic region."); 
	  ac_free (tt) ;
	}
    }
  
  if (ppp->dbSpecies == 'h')
    {/************************** EBI *****************************/
      AC_TABLE tt = 0 ;
      if ( (tt = ac_tag_table (ppp->geneBox, "IntMap", 0)) && tt->cols >= 3)
	{
	  nnam = ac_table_printable (tt, 0, 0, "") ;
	  sprintf (linkBuf,
		   "http://www.ensembl.org/Homo_sapiens/contigview?chr=chr%s&vc_start=%d&vc_end=%d",
		   nnam, ac_table_int (tt, 0, 1, 0), ac_table_int (tt, 0, 2, 0)) ;
	  oneLineLink ("Ensembl", linkBuf,"displays the EBI view of the same genomic region. The AceView genes are displayed as an option of the DAS sources.");
	} 
      ac_free (tt) ;
    }
  if (ppp->dbSpecies == 'w')
    oneLineLink ("elegans net", "http://members.tripod.com/C.elegans/", "interesting site") ;
  printf ("</table><br/>\n<br/>\n") ;
}

/***************************************************************************/
/***************************************************************************/
/**Example: run '?db=human&&c=MRNA&a=fiche&l=G_t3_chr3_545_1867.a'  ********/
/***************************************************************************/

static void doAction (PP *ppp)
{
  char *target ;
  unsigned char *fiche = 0 ;
  char * line ;
  unsigned char *cp, *cq ;
  int i, size ;


  switch (ppp->picture_type)
    {
    case GIF_T:  /* on call back from javascript/main.ppp->version.js:writeMap () */
      printf ("Content-type: image/gif\n\n") ;
      ppp->picCommand = "gif -nobox -" ;
      sleep (1) ; /* to be sure the text is requested first */
      break ;
    case BUBBLE_T:
      standardHtmlHeader (ppp, TRUE) ;
      ppp->picCommand = "bubble -" ; 
      break ;
    case SWF_T:  /* complete SWF display */
      standardHtmlHeader (ppp, TRUE) ;
      if (ppp->action == GENE_P) { ppp->height = 700; ppp->width = 800 ; }
      target = messprintf ("av.cgi?db=%s&c=%s&q=%s&s&a=%s&ww=%d&wh=%d&v=%d&h=%d"
			   , ppp->dbName
			   , ppp->action == VmRNA_P ? "mrna" : "gene"
			   ,  ppp->action == VmRNA_P ? url_encode(ac_name (ppp->mrna)) : url_encode(ac_name (ppp->geneBox))
			   , ppp->actionName
			   , ppp->width
			   , ppp->height
			   , ppp->details 
			   , ppp->host
			   ) ;

      printf ("<EMBED SRC=\"%s\" bgcolor='#ffffff' "
	      "allowScriptAccess='sameDomain' "
	      "ALIGN=\"LEFT\" WIDTH=\"%d\" HEIGHT=\"%d\" QUALITY=\"high\"\n"
	      , target
	      , ppp->width
	      , ppp->height) ;
      printf ("TYPE=\"application/x-shockwave-flash\" \n") ;
      printf ("PLUGINSPAGE=\"http://www.macromedia.com/go/getflashplayer\">\n") ;
      printf ("</EMBED>\n") ;
      
  /* create a special division that will be dynamically transfered to the tabs frame */
      if (strcmp (ppp->version,"v47") && (ppp->action == GENE_P || ppp->action == VmRNA_P))
	{
	  printf ("<div id='999'>\n") ;
	  printf ("  <table width=\"98%%\" border=2 bordercolor='blue'>") ;
	  printf ("    <tr>\n") ;
	  printf ("      <td>?</td>\n") ;
	  printf ("      <td>\n") ;
	  printf (        "<a href='javascript:fireCommand (1)'>Gene Summary</a>") ;
	  printf ("      </td>\n") ;
	  printf ("      <td>\n") ;
	  printf ("This plot is in Flash, use right mouse button to zoom in, then drag. \n") ;
	  printf ("If you prefer, the same views are accessible in gif form the \'Gene Summary\' page\n") ;
	  printf ("      </td>\n") ;
	  printf ("    </tr>\n") ;
	  printf ("  </table>\n") ;
	  printf ("</div>\n") ;
	}
      printf ("\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
	      "  <!--\n"
	      "     aceViewOnLoad() ;\n"
	      "  //-->\n</script>\n"
	      ) ;
      printf ("<div id='executor' onclick='jsExecutor()'>\n") ;
      printf ("\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
	      "  <!--\n"
	      "     var cannedCommands = new Array (\n\"toto\", \"openAceViewAction ('gene', '%s', 'fgene')\");" 
	      "  //-->\n"
	      "</script>\n"
	      "</div>\n"
	      , url_encode (ppp->query)
	      ) ;
      printf ("</body></html>\n") ;
      ac_command (ppp->db, messprintf (" // SWF_T done %s %s gene %s" , ppp->date
				       , ppp->picCommand ? ppp->picCommand : ""
				       , ac_name(ppp->geneBox)), 0, 0) ;
      ppp->picCommand = 0 ;
      break ;
    case swf_T:  /* embeded subsgraph */
      if (1) printf ("Expires: Mon, 31 Dec 2007 23:59:59 UTC\n") ;
      if (strstr (ppp->browser, "MSIE")) printf ("Content-type: application/x-shockwave-flash\n\n") ; /* IE */
      else                               printf ("Content-type: application/x-shockwave-flash\n\n") ; /* Firefox */
      ppp->picCommand = " swfdump -compile -";
      break ;
    default:
      break ;
    }
 /* http header */
  switch (ppp->action) 
    {
    case ZERO_A:
      showFramesMain (ppp) ;
      break ;
    case LOG_A: 
      standardHtmlHeader (ppp, TRUE) ;
      executeAction (ppp) ;
      break ;
    case INFO_A:
      showFramesInfo (ppp) ;
      break ;
    case FGENE_A:
    case FMOL_A:
    case FEXP_A:
    case FFUNC_A:
    case FICHE_A:
    case GOLD_A:
    case FASTA_A:
    case TG_SUPPORT_A:
    case CLONES_A:
      standardHtmlHeader (ppp, TRUE) ;

      printf ("<script language=\'JavaScript\' type=\'text/javascript\'>\n") ;
      printf ("if (window.name.toLowerCase () == \"graph\")document.write (\"<a href='javascript:openFullPage (document.location)' >Full Page</a> &nbsp; | &nbsp; \") ;\n") ;
      printf ("</script>\n") ;

      executeAction (ppp) ; 
      if (0) system ("printenv >> /tmp/miegenv") ;
      printf ("</body></html>\n") ;
      break ;

    case SECTION_A:
      printf ("Content-type: text/html\n\n") ;
      executeAction (ppp) ; 
      break ;

    case MEGA_A:
      standardHtmlHeader (ppp, TRUE) ;
      fiche = ac_command (ppp->db
			 , messprintf ("view -v fiche -c %s -n \"%s\" -xml", ppp->class, ppp->query) 
			 ,0, 0) ;
      printf ("%s\n", fiche) ;
      messfree (fiche) ;
      printf ("</body></html>\n") ;
      break ;

    case  GENE_P:  /* on call back from javascript/main.ppp->version.js:writeMap () */
      if (ppp->picCommand && !strcmp(ppp->picCommand, " swfdump -compile -"))
	exportVGeneFlashPicture (ppp) ;
      else if (ppp->picCommand && !strcmp(ppp->picCommand, " svg"))
	exportVGeneSvgPicture (ppp) ;
      else if (ppp->picCommand)
	exportVGeneGifPicture (ppp) ;
      break ;
      
    case  VmRNA_P:  /* on call back from javascript/main.ppp->version.js:writeMap () */
     if (ppp->picCommand && !strcmp(ppp->picCommand, " svg"))
	exportVmrnaSvgPicture (ppp) ;
     else if (ppp->picCommand) 
       exportVmrnaPicture (ppp) ;
      break ;
      
    case HmRNA_P:  /* on call back from javascript/main.ppp->version.js:writeMap () */
     if (ppp->picCommand) 
       exportHmrnaPicture (ppp) ;
      break ;
      
    case Locator_P:  /* on call back from javascript/main.ppp->version.js:writeMap () */
     if (ppp->picCommand) 
       exportLocatorPicture (ppp) ;
      break ;
      
    case LocatorBig_P:  /* on call back from javascript/main.ppp->version.js:writeMap () */
     if (ppp->picCommand) 
       exportLocatorBigPicture (ppp) ;
      break ;
      
    case PROBE_P:  /* on call back from javascript/main.ppp->version.js:writeMap () */
     if (ppp->picCommand)
       exportProbePicture (ppp) ;
      break ;

    case LIST_A:
      standardHtmlHeader (ppp, TRUE) ;
      showContextQuery (ppp) ;
      showGeneList (ppp) ;
      printf ("</body></html>\n") ;
      break ;
      
    case LINK_A:
      standardHtmlHeader (ppp, TRUE) ;
      showLinkList (ppp) ;
      printf ("</body></html>\n") ;
      break ;

    case TREE_A:
      standardHtmlHeader (ppp, TRUE) ;
      if (getOneKey (ppp))
        showTreeObject (ppp) ;
      printf ("</body></html>\n") ;
      break ;

    case GOLDFILE_A:
      showGoldFile (ppp) ; /* no decoration needed for the moment */
      break ;

    case DNA_A:
      standardHtmlHeader (ppp, TRUE) ;
      if (getOneKey (ppp))
	{
	  if (ppp->geneBox)
	    showMOI (ppp) ;
	  showDnaObject (ppp) ;
	}
      printf ("</body></html>\n") ;
      break ;


    case PEP_A:
      standardHtmlHeader (ppp, TRUE) ;
      if (getOneKey (ppp))
        showPepObject (ppp) ;
      printf ("</body></html>\n") ;
      break ;


    case ERROR_A:
    default:   /* Error */
      ppp->noRobot = TRUE ;
      ppp->noJava = TRUE ;
      standardHtmlHeader (ppp, TRUE) ;
      reportError (ppp, ppp->query);
      printf ("</body></html>\n") ;
      break ;
    }

  if ((size = ppp->pictureBufferSize))
    switch (ppp->picture_type)
      {
      case swf_T:
	cp = (unsigned char *) strstr ((char *)ppp->pictureBuffer, "//") ;
	if (! cp) cp = ppp->pictureBuffer ;
	else { while (*cp++ != '\n') {} ; cp++ ; }
	cq = ppp->pictureBuffer + size - 1 ;
	while (size > 0 && *cq != '/') { cq-- ; size-- ; }
	cq -= 2 ; size -= 2 ;
	i = cp - ppp->pictureBuffer ;
	size = size - (cp - ppp->pictureBuffer) ;
	fwrite (cp, size, 1, stdout) ;
	exportJavaMapCallback () ;
	if (ppp->db) ac_command (ppp->db, "// swf file received by aceviewmain", 0, 0) ;
	break ;
      case SVG_T:
	cp = (unsigned char *) strstr ((char *)ppp->pictureBuffer, "//") ;
	if (! cp) cp = ppp->pictureBuffer ;
	else { while (*cp++ != '\n') {} ; cp++ ; }
	cq = ppp->pictureBuffer + size - 1 ;
	while (size > 0 && *cq != '/') { cq-- ; size-- ; }
	cq -= 2 ; size -= 2 ;
	i = cp - ppp->pictureBuffer ;
	size = size - (cp - ppp->pictureBuffer) ;
	fwrite (cp, size, 1, stdout) ;
	exportJavaMapCallback () ;
	if (ppp->db) ac_command (ppp->db, "// SVG file received by aceviewmain", 0, 0) ;
	break ;
      case GIF_T:
	if (ppp->db) ac_command (ppp->db, "// gif file received by aceviewmain", 0, 0) ;
	for (i = 0 ; i<size-5 ; i++)
	  if (ppp->pictureBuffer[i] == 'G' && !memcmp (ppp->pictureBuffer + i, "GIF89a", 5)) 
	    break ;
	if (ppp->pictureBuffer)
	  fwrite (ppp->pictureBuffer+i, size-i, 1, stdout) ;
	break ;
      case BUBBLE_T:   /* compress and dump bubble info */	
	printf ("\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
		"  <!--\n") ;
	printf ("  var viewDetails=%x ;\n", ppp->details) ;
	printf ("  var picture_type='G' ;\n") ;
	printf ("  var gene=\"%s\" ;\n", ac_name (ppp->geneBox)) ;
	if (ppp->mrna) printf ("  var mrna = \"%s\" ;\n", ac_name (ppp->mrna)) ;
	printf ("  var seqstart=%d ;\n", ppp->seqstart) ;
	printf ("  var seqlen=%d ;\n", ppp->seqlen) ;
	if (ppp->action == VmRNA_P)
	  {
	    int ir ;
	    AC_OBJ oTg ;
	    AC_TABLE gmrna = 0 ;
	    
	    oTg = ac_tag_obj (ppp->mrna, "From_gene", 0) ;
	    /* AC_OBJ	    oGenefinder = ac_tag_obj (ppp->mrna, "From_prediction", 0) ; */
	    
	    
	    if (
		(oTg && (gmrna = ac_tag_table (oTg, "mRNA", 0))) ||
		(ppp->geneBox && (gmrna = ac_tag_table (ppp->geneBox, "Product", 0)))
		)
	      {
		printf ("  var mrnaList = new Array (\n") ;
		for (ir = 0 ;ir < gmrna->rows ;ir++)
		  printf ("    \"%s\", ", ac_table_printable (gmrna, ir, 0, "")) ;
		printf ("    \"\") ;\n") ;
	      }
	    ac_free (gmrna) ;
	  }
	printf ("  var allFields=new Array (\n") ;
	/* compress and dump bubble info */
	for (line = strtok ((char *)ppp->pictureBuffer, "\n"), i = 0 ; line && line[0]  != '>' ; line = strtok (NULL, "\n"))
	  {
	    if (line[0] == '/') continue ;
	    printf ("    \"%s\", \n", line) ;
	  }
	printf ("    \"\" ) ;\n") ;
	printf ("  //-->\n</script>\n") ;
      
  /* create a special division that will be dynamically transfered to the tabs frame */
      if (strcmp (ppp->version,"v47") && (ppp->action == GENE_P || ppp->action == VmRNA_P))
	{
	  printf ("<div id='999'>\n") ;
	  printf ("  <table width=\"98%%\" border=2 bordercolor='blue'>") ;
	  printf ("    <tr>\n") ;
	  printf ("      <td>?</td>\n") ;
	  printf ("      <td>\n") ;
	  printf (        "<a href='javascript:fireCommand (1)'>Gene Summary</a>") ;
	  printf ("      </td>\n") ;
	  printf ("      <td>\n") ;
	  printf ("This plot is in Gif, use the icons on the left to scroll and zoom. \n") ;
	  printf ("If you prefer, a similar view is accessible in Flash format from the \'Gene Summary\' page\n") ;
	  printf ("      </td>\n") ;
	  printf ("    </tr>\n") ;
	  printf ("  </table>\n") ;
	  printf ("</div>\n") ;
	}

	exportJavaMapCallback () ;
	printf ("<div id='executor' onclick='jsExecutor()'>\n") ;
	printf ("\n<script language = \'JavaScript\' type=\'text/javascript\'>\n "
		"  <!--\n"
		"     var cannedCommands = new Array (\n\"toto\", \"openAceViewAction ('gene', '%s', 'fgene')\");" 
		"  //-->\n"
		"</script>\n"
		"</div>\n"
		, url_encode (ppp->query)
		) ;
	printf ("</body></html>\n") ;
	ac_command (ppp->db, messprintf (" // bubble done %s %s gene %s" , ppp->date, ppp->picCommand, ac_name(ppp->geneBox)), 0, 0) ;
	break ;
      default:
	break ;
      }
}

/***************************************************************************/
/** Example: run '?db=human&&c=MRNA&a=fiche&l=G_t3_chr3_545_1867.a' ********/
/***************************************************************************/

char *C_html_extra_message_generator (int nn)
{
  static ACTION_TYPE action = ZERO_A ;

  if (nn > 0 && nn < LAST_A) 
    action = (ACTION_TYPE) nn ;
  switch (action)
    {
    case ZERO_A: return "ZERO" ;
    case TREE_A: return "Tree display" ;
    case  DNA_A: return "DNA sequence" ;
    case  PEP_A: return "Peptide sequence" ;
    case  PROBE_A: return "Probe" ;
    case  FICHE_A: return "Fiche" ;
    case  FGENE_A: return "Gene page" ;
    case  FMOL_A: return "Molecular page" ;
    case  FEXP_A: return "Expression page" ;
    case  FFUNC_A: return "Function page" ;
    case  SECTION_A: return "Section" ;
    case  GENE_A: return "Gene graphics" ;
    case  MRNA_A: return "mRNA graphics" ;
    case  LIST_A: return "List of genes" ;
    case  CLONES_A: return "List of clones" ;
    case  TG_SUPPORT_A: return "Clones supporting a given gene" ;
    case  GOLD_A: return "Gold paper" ;
    case  GOLDFILE_A: return "Gold auxiliary data" ;
    case  FASTA_A: return "Fasta of all sequences related to a gene" ;
    case  ERROR_A: return "Error" ;
    case  MEGA_A: return "Mega" ;
    case  LINK_A: return "Link" ;
    case  INFO_A: return "Info container page" ;
    case  LOG_A: return "Log" ;
    case  GENE_P: return "Gene_p" ;
    case  VmRNA_P: return "mRNA vertical graphics" ;
    case  Locator_P: return "Locator graphics" ;
    case  LocatorBig_P: return "Locator_Big graphics" ;
    case  HmRNA_P: return "mRNA horizontal graphics" ;
    case  PROBE_P: return "Probe" ;
    case  LAST_A: return "LAST" ;
    }
  return "ZERO" ;
}

int C_html_page_generator (void)
{
  PP pp ;
  char *exdb ;
  char *arg = 0 ;
  memset (&pp, 0, sizeof (PP)) ;
  setbuf(stdout,NULL);

  /*
  * cgi debugging hack - you can pass in the query string as argv[1]
  * and the CGI system will believe it is in a real web server
  */
  if (arg)
    {
      char *s = messprintf("QUERY_STRING=%s",arg) ;
      s=strdup(s);		/* putenv needs the string forever */
      putenv(s) ;
      putenv("REQUEST_METHOD=GET");
    }

  /*
  * You can now use the functions in wrpc/cgi.[ch] in this program.
  * in particular cgi_clean_arg() 
  */
  cgi_init("aceview.log");

  /*
  * if an external database is specified, send a redirect.
  */
  exdb = cgi_clean_arg("exdb");
  if (exdb && (strcmp(exdb,"AceView") != 0))
    {
    char *term;
    term=cgi_clean_arg("term");
    if (!term)
      term="";
    fprintf(stderr,"%s %s %s\n",logPrefix(0, "redirect"),exdb,term);
    printf("Location: https://www.ncbi.nlm.nih.gov/genome/guide/gquery.cgi?db=%s&term=%s\n\n",
	url_encode(exdb), url_encode(term));
    return 0 ;
    }

  memset (&pp, 0, sizeof (PP)) ;
  if (! getParams (&pp))
    {
      /* 
       * missing cgi parameters, no such database, database 
       * valid but unreachable.
       *
       * no log message here - getParams did it
       */
      pp.action = ERROR_A;
      C_html_extra_message_generator (pp.action) ;
      pp.geneBox = NULL;
      standardHtmlHeader (&pp, FALSE) ;
      printf("<h1>cgi parameter error %s arg=%s</h1>\n",pp.query,arg);
      if (0) printf("<a href=\"index.html\" target='_top'>click here to return to query page</a>\n");
      return 1 ;
    }
  else
    {
      /*
       * initial parameters were good and we opened the database.
       * try to do something.
       */
      C_html_extra_message_generator (pp.action) ;
      if (pp.action == GOLDFILE_A) ; /* just return the file from the server */
      else if (pp.action == LOG_A) ; /* just log the comment on the server */
      else if (pp.action == DNA_A)
	{
	  char buf[1024], *cp, *cq, *cr ;
	  strcpy (buf, pp.class) ;

	  if (!strcasecmp (buf, "dna"))
	    strcpy (pp.class, "sequence") ;
	  else if (!strncasecmp (buf, "dna:", 4) &&
		   (cp = strstr (buf+4,":")) )
	    {
	      *cp = 0 ;
	      strcpy (pp.class,"sequence") ; 
	      cq = strstr (cp+1,":") ;
	      if (cq)
		{
		  *cq = 0 ; 
		  pp.dnacolor = atoi (cq+1) ;
		}
	      pp.dna1 = atoi(buf+4) ;
	      pp.dna2 = atoi (cp+1) ;
	    }
	  else if (!strncasecmp (buf, "dna_", 4) &&
		   (cp = strstr (buf+4,":")) )
	    {
	      *cp = 0 ;
	      strcpy (pp.class,buf+4) ;
	      cq = strstr (cp+1,":") ;
	      if (cq)
		{
		  *cq = 0 ; 
		  cr = strstr (cq+1,":") ; 
		  if (cr)
		    {
		      *cr = 0 ;
		      pp.dnacolor = atoi (cr+1) ;
		    }
		  pp.dna2 = atoi (cq+1) ;
		}
	      pp.dna1 = atoi(cp+1) ;
	    }
	}

      else if (pp.action == PEP_A)
	{
	  char buf[1024], *cp, *cq ;
	  strcpy (buf, pp.class) ;
	  
	  if (!strncasecmp (buf, "pep_", 4) &&
	      (cp = strstr (buf+4,":")) )
	    {
	      *cp = 0 ;
	      strcpy (pp.class,buf+4) ;
	      cq = strstr (cp+1,":") ;
	      *cq = 0 ;
	      pp.dna1 = atoi(cp+1) ;
	      pp.dna2 = atoi (cq+1) ;
	    }
	}

      else if ( pp.action != TREE_A 
	   && !getProducts (&pp) 
	   && !getMrnas (&pp) 
	   && !getGenes (&pp)  )
	{
	  /* pp.error is reported by acedb webquery system */
	  char *errp = 0 ;
	  if ( pp.error && strstr (ac_name(pp.error), "too many hits"))
	    errp = logPrefix (&pp, "saturated") ;
	  else
	    errp = logPrefix (&pp, "no_results") ;
	  fprintf(stderr,"%s\n", errp) ;
	  if (pp.db) ac_command (pp.db, errp, 0, 0) ;
	  pp.action = ERROR_A;
	  pp.geneBox = NULL;
	  pp.noRobot = TRUE ;
	  standardHtmlHeader (&pp, FALSE) ;
	  reportError (&pp, pp.query) ;
	  exit(0);
	}
    }
  
  logClientInfo (&pp) ;

  if (pp.browser && *pp.browser && pp.newWindow)
    showFramesMain (&pp) ;
  else
    doAction (&pp) ;

  messfree (pp.args) ;

  ac_free (pp.db) ;
  return 0 ;
}  /* main */

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
