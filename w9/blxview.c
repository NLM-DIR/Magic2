/*  Last edited: Dec 16 15:03 1998 (fw) */

/* $Id: blxview.c,v 1.9 2017/06/05 20:59:41 mieg Exp $ */

/* 

   BLIXEM - BLast matches In an X-windows Embedded Multiple alignment

-------------------------------------------------------------
|  File: blxview.c                                          |
|  Author: Erik Sonnhammer (esr@ncbi.nlm.nih.gov)           |
|  Copyright (C) E Sonnhammer 1993-1997                     |
-------------------------------------------------------------

  Date      Modification
--------  ---------------------------------------------------
93-02-20  Created
93-05-25  Dispstart/dispend fix by Richard for seq's < displen.
93-07-24  All boxes of a protein turn black when one is picked
          Sorting by protein name or score added
93-11-16  Added picking of Big Picture HSP's and Reverse Strand Big Picture
93-11-17  Added blastn support
94-01-29  Added Highlight sequences by names matching regexp
94-03-07  Fixed window limits. Added initial settings (BigPict, gotoNext)
94-03-08  Added 'always-both-strands' in blastn mode.
94-03-27  Added Tblastn support (works in seqbl mode).
94-12-01  [2.0] Added dotter calling.
95-02-02  Rewrote auxseq and padseq allocation to be fully dynamic and unlimited.
95-02-06  Added Tblastx support (works in seqbl mode).
95-06-01  [2.1] Added SEG display
95-06-23  [2.1] Added Settings window
95-07-21  2.1 announced--------------------------------------
95-08-01  [2.2] Initial Sorting mode on command line.
95-09-15  [2.2] Added Settings pull-down menu for faster manipulation.
95-09-29  [2.2] Improved WWW browser finding with findCommand() - doesn't get fooled so easily.
95-10-04  [2.2] Added "Print whole alignment"
95-10-27  [2.2] Added acedb-fetching at double clicking. 
          BLIXEM_FETCH_ACEDB makes this default.
          Reorganised Settings window to Toggles and Menus rows.
95-11-01  [2.3] Changed command line syntax to "-" for piping.
          Added X options capability (-acefont, -font).
96-01-05  [2.3] Force all tblastn HSP's to be qframe +1, to harmonize with MSPcrunch 1.4
          which gives sframe for tblastn (otherwise the output would be dead).
96-02-08  [2.3] Added option -S "Start display at position #" to stand-alone command line.
96-02-08  [2.3] Added checkmarks to pull-down settings menu.
96-03-05  2.3 Announced.
96-03-08  [2.4] Only look for WWW browser once.
96-05-09  [2.4] Proper grayscale print colours.
96-05-09  [2.4] Enabled piping of query sequence too, for Pepmap and WWW calls.
96-08-20  [2.4] Fixed minor bug in squashed mode and added restoring of previous sorting after squash.
97-05-28  [2.4] Fixed parsing to handle gapped matches.
                Added "Highlight lower case residues" for gapped alignments and
                "Show sequence descriptions" (for MSPcrunch 2.1).
		Added setting the color of matching residues in the Settings window.
97-06-17  [2.4] Fixed "Highlight differences" for gapped alignments ('.' -> '-').
                Changed "Highlight lower case residues" to "Highlight subject insertions" and 
		set this automatically for gapped alignments.  Works for both lower case and number 
		insert markers.
                Changed blviewRedraw to use strlen to accommodate reverse gapped alignments.
		Simplified (and thereby debugged) selection of Big Picture MSPs to be drawn.
		Made Big Picture ON/OFF control more logical and consistent.
		Added a calcID() step to fix sortById() at startup.
		Added "Hilight Upper/Lower case" - useful general purpose feature.
97-10-14  [2.4] Added Description box when picking sequences.

-------------------------------------------------------------

Known bugs:
-----------

revcomp broken when called from acedb.  Slen/send problem?


MSP score codes:
----------------
-1 exon                  -> Big picture + Alignment
-2 intron                -> Big picture + Alignment
-3 Any coloured segment  -> Big picture
-4 stringentSEGcolor     -> Big picture
-5 mediumSEGcolor        -> Big picture
-6 nonglobularSEGcolor   -> Big picture
-10 hidden by hand
*/

#include <math.h>
#include <ctype.h>
/*#include "acedb.h"*/
#include "graph.h"
#include "key.h"
#include "../wh/menu.h"
#include "dotter.h"
#include "utils.h"

#ifdef ACEDB
#include "display.h"
#endif

extern void externalCommand (char *command);
char *blixemVersion = "2.4";

#ifdef SUN
/*   extern void  pclose(FILE*); */ /* fw 940906 now defined in mystdlib.h */
  extern int   atoi(char*),
               re_exec(char*);
  extern char *re_comp(char*)/*,
               toupper(char),
               tolower(char)*/;
#endif

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define BPoffset 4
#define NAMESPACE 12
#define SEQ2BP(i) (float)plusmin*(i-BigPictStart-qoffset)*BPx/BigPictLen + BPoffset
#define MAXALIGNLEN 10000

typedef struct _BPMSP {
    char sname[FULLNAMESIZE+1];
    char *desc;
    int box;
    struct _BPMSP *next;
} BPMSP;

enum { UNSORTED, SORTBYSCORE, SORTBYID, SORTBYNAME, SORTBYPOS };

static BPMSP BPMSPlist, *bpmsp;
static int  fetchCount;
static void scrollRightBig(void);
static void scrollLeftBig(void);
static void scrollRight1(void);
static void scrollLeft1(void);
static void Goto(void);
static void gotoMatch(int direc);
static void prevMatch(void) ;
static void nextMatch(void);
static void keyboard (int key);
static void toggleColors (void);
static void blviewPick (int box);
static void blviewHelp(void);
static Graph blviewCreate (void);
static void blviewRedraw(void);
static void blviewDestroy(void);
static void blviewDestroy(void);
static char *getqseq(int start, int end, char *q);
static char *get3rd_base(int start, int end, char *q);
static void calcID(MSP *msp);
static void getsseq(MSP *msp),
            BigPictToggle(void),
            SEGtoggle(void),
            BigPictToggleRev(void),
            zoomIn(void),
            zoomOut(void),
            zoomWhole(void),
            MiddleDrag(double x, double y), 
            MiddleUp(double x, double y), 
            MiddleDown(double x, double y),
            setHighlight(void),
            clrHighlight(void),
            callDotter(void),
            callDotterHSPs(void),
            callDotterSelf(void),
/*            dotterPanel(void), */
            setDotterParams(void),
            autoDotterParams(void),
            allocAuxseqs(int len),
            blixemSettings(void),
            settingsRedraw(void),
            menuCheck(MENU menu, int mode, int thismode, char *str),
	    hidePicked(void),
            setMenuCheckmarks(void);

static Graph g,        /* The blxview window */
             settingsGraph;

static int   actstart,        /* Active region coordinates relative to query in BASES */
             actend,     
             dispstart,       /* Displayed region relative to query in BASES */
             displen=240;     /* Displayed sequence length in BASES */
static char  actframe[16]="(+1)";    /* Active frame */
static int   plusmin=1;       /* +1 for top strand, -1 */
static float queryy,          /* Screen coords of genome seq */
             separator_y,
             lastExonx,
             BPboxheight = 5.7,
             BPboxstart,
             BPboxwidth,
             BigPictZoom = 10,
             oldGraphWidth = 0, oldGraphHeight = 0;
static char *q;    /* The genomic sequence (query=the translated seq) */
static int   qlen;
static int   qoffset;         /* Offset to the 'real' DNA start */
static char *qname;
static MSP   DBhits,          /* List of MSP's */
            *pickedMSP;	      /* Last picked MSP */
static int   messagebox;
static char  message[1024], 
             HighlightSeq[NAMESIZE+4] = "", 
             searchSeq[NAMESIZE+4] = "", 
            *cp,
            *padseq = 0,
            *auxseq = 0,
            *auxseq2 = 0,
             dotterqname[NAMESIZE+1],
            *dottersseq=0,
             stringentSEGtx[10],
             mediumSEGtx[10],
             nonglobularSEGtx[10],
             sortModeStr[32] = "      ",
             fetchMode[32] = "efetch", /* Not done with enum to get menu-strings for free */
            *URL=0;
static int   lastbox = 0;
static int   colortoggle = 0;
static int   backgColor = LIGHTGRAY,
             IDcolor = CYAN,
             consColor = MIDBLUE,
             geneColor = BLUE,
             hiColor = YELLOW,
             oldcolor,
             BigPictON = 1,
             BigPictRev = 0,
             BigPictStart, BigPictLen, BPbox, BPx, 
             gridColor = YELLOW,
             nx, ny,
             blastx = 1,
             blastp = 0,
             blastn = 0,
             tblastn = 0,
             tblastx = 0,
             symbfact = 3,
             i,
             start_nextMatch = 0,
             dotter_first = 0,
             compN,
             IDdots = 0,
             squash = 0,
             verbose = 0,
             HiliteSins = 0,
             HiliteUpperOn = 0,
             HiliteLowerOn = 0,
             DESCon = 0,
             dotterZoom = 0,
             dotterStart = 0,
             dotterEnd = 0,
             dotterHSPs = 0,
             auxseqlen = 0,
             smartDotter = 1,
             SEGon = 0,
             stringentSEGcolor = LIGHTGREEN,
             stringentSEGbox,
             mediumSEGcolor = GREEN,
             mediumSEGbox,
             nonglobularSEGcolor = DARKGREEN,
             nonglobularSEGbox,
             alphabetsize,
             stringentSEGwin = 12,
             mediumSEGwin = 25,
             nonglobularSEGwin = 45,
             printColorsOn,
             wholePrintOn,
             settingsButton,
             sortMode = UNSORTED,
             HSPgaps = 0;

static double oldx, DNAstep;
static BOOL   dragFast ;
static Array stringentSEGarr, mediumSEGarr, nonglobularSEGarr;

static MENU blixemMenu ;	/* Main menu - contains function calls, etc. */
static MENU settingsMenu;	/* Contains toggles and 'states' */

#define SortByScoreStr      "Sort by score"
#define SortByIdStr         "Sort by identity"
#define SortByNameStr       "Sort by name"
#define SortByPosStr        "Sort by position"
#define BigPictToggleStr    "Big Picture"
#define BigPictToggleRevStr "Big Picture Other Strand"
#define toggleIDdotsStr     "Highlight differences"
#define squashMatchesStr    "Squash matches"
#define SEGtoggleStr        "Complexity/SEG panel"
#define printColorsStr      "B/W Print colours"
#define toggleColorsStr     "No colours"
#define toggleVerboseStr    "Verbose mode"
#define toggleHiliteSinsStr "Highlight subject insertions"
#define toggleHiliteUpperStr "Highlight upper case"
#define toggleHiliteLowerStr "Highlight lower case"
#define toggleDESCStr       "Show sequence descriptions"

static void toggleVerbose(void) {
    MSP *msp;

    verbose = (verbose ? 0 : 1);

    if (verbose) for (msp = DBhits.next; msp; msp = msp->next) {
#ifdef ACEDB
	printf("%d %s ", msp->key, name(msp->key));
#endif
	printf("%s %d %d %d %d %d %d :%s:\n", 
	       msp->sname, msp->score, msp->id, 
	       msp->qstart+qoffset, msp->qend+qoffset, 
	       msp->sstart, msp->send, (msp->sseq ? msp->sseq : ""));
    }
    blviewRedraw();
}
static void toggleHiliteSins(void) {
    HiliteSins = (HiliteSins ? 0 : 1);
    blviewRedraw();
}
static void toggleHiliteUpper(void) {
    HiliteUpperOn = (HiliteUpperOn ? 0 : 1);
    blviewRedraw();
}
static void toggleHiliteLower(void) {
    HiliteLowerOn = (HiliteLowerOn ? 0 : 1);
    blviewRedraw();
}
static void toggleDESC(void) {
    DESCon = (DESCon ? 0 : 1);
    blviewRedraw();
}
static void ToggleStrand(void) {
    dispstart += plusmin*(displen-1);
    plusmin = -plusmin;
    sprintf(actframe, "(%+d)", plusmin);
    blviewRedraw();
}
static void scrollRightBig(void)
{
    dispstart = dispstart + plusmin*displen*.5;
    blviewRedraw();
}
static void scrollLeftBig(void)
{
    dispstart = dispstart - plusmin*displen*.5;
    blviewRedraw();
}
static void scrollRight1(void)
{
    dispstart = dispstart + plusmin*symbfact;
    blviewRedraw();
}
static void scrollLeft1(void)
{
    dispstart = dispstart - plusmin*symbfact;
    blviewRedraw();
}

static void Goto(void)
{
    static char dfault[32];
    int i;

    /*sprintf(posstr, "%d", dispstart + qoffset);*/

    if (!graphPrompt ("Goto which position: ", 
		      dfault, 
		      "t"))
	return;

    freeint(&i);
    dispstart = i - qoffset;
    sprintf(dfault, "%d", i);

    blviewRedraw();
}

static void gotoMatch(int direc)
{
  MSP *msp ;
    int offset, closest_offset;
    char strand;

    if ( strchr(actframe, '+')) strand = '+';
    else strand = '-';

    if (direc != -1 && direc != 1) {
	messerror ( "gotoMatch must have -1 or 1 as argument\n" ) ;
	return;
    }

    closest_offset = 0;

    for (msp = DBhits.next; msp ; msp = msp->next)
	if ( strchr(msp->frame, strand) &&
	    (msp->qstart - dispstart)*plusmin*direc - 2 > 0)
	{
	    offset = (msp->qstart - dispstart)*plusmin*direc - 2;
	    if (!closest_offset) closest_offset = offset;
	    else if (offset < closest_offset) closest_offset = offset;
	}

    if (!closest_offset) {
	blviewRedraw();
	return;
    }

    if (direc < 0) dispstart -= 4*plusmin;

    dispstart = dispstart + direc*plusmin*(closest_offset );
    blviewRedraw();
}

static void prevMatch(void) { gotoMatch(-1); }
static void nextMatch(void) { gotoMatch(1); }

static void gotoBegin(void)
{
    dispstart = (plusmin > 0 ? 0 : qlen);
    blviewRedraw();
}
static void gotoEnd(void)
{
    dispstart = (plusmin > 0 ? qlen : 0);
    blviewRedraw();
}

static void keyboard (int key)
{
    switch (key)
	{
	case '<': case ',': scrollLeft1();  break;
	case '>': case '.': scrollRight1(); break;
	case UP_KEY:    blviewPick(lastbox - 1);  break;
	case DOWN_KEY:  blviewPick(lastbox + 1);  break;
	default: return;
	}
}


static void toggleColors (void)
{
    static int oldback, oldgrid, oldID, oldcons, oldgene, oldhi;

    graphActivate(settingsGraph);

    if (!colortoggle) {
	oldback = backgColor; backgColor = WHITE;
	oldgrid = gridColor; gridColor = BLACK;
	oldID = IDcolor; IDcolor = WHITE;
	oldcons = consColor; consColor = WHITE;
	oldgene = geneColor; geneColor = BLACK;
	oldhi = hiColor; hiColor = WHITE;
	colortoggle = 1;
    }
    else {
	backgColor = oldback;
	gridColor= oldgrid;
	IDcolor = oldID;
	consColor = oldcons;
	geneColor = oldgene;
	hiColor = oldhi;
	colortoggle = 0;
    }
    blviewRedraw();
}


static void printColors (void)
{
    static int oldback, oldgrid, oldID, oldcons, oldgene, oldhi;

    graphActivate(settingsGraph);

    if (!printColorsOn) {
	oldback = backgColor; backgColor = WHITE;
	oldgrid = gridColor; gridColor = LIGHTGRAY;
	oldID = IDcolor; IDcolor = GRAY;
	oldcons = consColor; consColor = PALEGRAY;
	oldgene = geneColor; geneColor = BLACK;
	oldhi = hiColor; hiColor = LIGHTGRAY;
	printColorsOn = 1;
    }
    else {
	backgColor = oldback;
	gridColor= oldgrid;
	IDcolor = oldID;
	consColor = oldcons;
	geneColor = oldgene;
	hiColor = oldhi;
	printColorsOn = 0;
    }
    blviewRedraw();
}


static void toggleIDdots (void)
{
    IDdots = !IDdots;
    blviewRedraw();
}


static void calcSEGarray(Array array, int win)
{
    int i, j, *rescount;
    float pi, sum;

    rescount = (int *)messalloc(24*sizeof(int));

    for (i=0; i < qlen; i++) {
	rescount[aa_atob[(int)q[i]]]++;
	if (i+1 >= win) {
	    for (sum = j = 0; j < 24; j++)
		if (rescount[j]) {
		    pi = (float)rescount[j]/win;
		    sum += pi*log(pi);
		}
	    arr(array, i-win/2, float) = -sum/log(2);
	    rescount[aa_atob[(int)q[i+1-win]]]--;
	}
    }

    messfree(rescount);

    /* TEST - delete later * /
    for (i=0; i < qlen; i++)
	printf ("%3d  %c  %f\n", i, q[i], arr(stringentSEGarr, i, float));
	*/
}

static void calcSEGarrays(BOOL force)
{
    /* force:
       FALSE - only calculate if necessary, i.e. first call.
       TRUE - force (re)calculation.
       */

    if (!stringentSEGarr) {
	calcSEGarray(stringentSEGarr = arrayCreate(qlen, float), stringentSEGwin);
	calcSEGarray(mediumSEGarr = arrayCreate(qlen, float), mediumSEGwin);
	calcSEGarray(nonglobularSEGarr = arrayCreate(qlen, float), nonglobularSEGwin);
    }
    else if (force) {
	calcSEGarray(stringentSEGarr = arrayCreate(qlen, float), stringentSEGwin);
	calcSEGarray(mediumSEGarr = arrayCreate(qlen, float), mediumSEGwin);
	calcSEGarray(nonglobularSEGarr = arrayCreate(qlen, float), nonglobularSEGwin);
    }
}


static void SEGtoggle (void)
{
    SEGon = !SEGon;
    if (SEGon) BigPictON = 1;
    calcSEGarrays(FALSE);
    blviewRedraw();
}


/* Find the expression 'query' in the string 'text'
 * Return 1 if found, 0 otherwise
 */
int strMatch(char *text, char *query)
{
    /* Non-ANSI bsd way: * /
       if (re_exec(text) == 1) return 1;
       else return 0;
    */

    /* ANSI way: */
    return pickMatch(text, query);
}

void highlightProteinboxes(void)
{
    MSP *msp;

    /* Highlight alignment boxes of current search string seq */
    if (*searchSeq) 
	for (msp = DBhits.next; msp ; msp = msp->next)
	    if (msp->box && strMatch(msp->sname, searchSeq))
		graphBoxDraw(msp->box, BLACK, RED);

    /* Highlight alignment boxes of currently selected seq */
    if (!squash)
	for (msp = DBhits.next; msp ; msp = msp->next)
	    if (msp->box && !strcmp(msp->sname, HighlightSeq)) 
		graphBoxDraw(msp->box, WHITE, BLACK);

    if (BigPictON)
    {
	/* Highlight Big Picture boxes of current search string seq */
	if (*searchSeq) 
	    for (bpmsp = BPMSPlist.next; bpmsp && *bpmsp->sname; bpmsp = bpmsp->next)
		if (strMatch(bpmsp->sname, searchSeq))
		    graphBoxDraw(bpmsp->box, RED, BLACK);

	/* Highlight Big Picture boxes of currently selected seq */
	for (bpmsp = BPMSPlist.next; bpmsp && *bpmsp->sname; bpmsp = bpmsp->next)
	    if (!strcmp(bpmsp->sname, HighlightSeq))
		graphBoxDraw(bpmsp->box, CYAN, BLACK);
    }
}

static void hidePicked (void)
{ 
    MSP *msp;

    for (msp = DBhits.next; msp ; msp = msp->next)
	if (!strcmp(msp->sname, HighlightSeq)) {
	    msp->id = msp->score;
	    msp->score = -999;
	}
    blviewRedraw () ;
}

static void blviewPick (int box)
{
    MSP *msp;
#ifndef ACEDB	       
    static char *browser=0;
#endif /* !ACEDB */

    extern int *findCommand (char *command, char **retp);

    if (!box) return ;

    if (box == lastbox)
    {  
	/* Second click - efetch this seq */
	for (msp = DBhits.next; msp ; msp = msp->next)
	    if (msp->box == box) break;

	if (msp && *msp->sname) 
	{
	    if (!strcmp(fetchMode, "efetch")) 
		externalCommand(messprintf("efetch '%s' &", msp->sname));
	    else if (!strcmp(fetchMode, "WWW-efetch"))
	    {
#ifdef ACEDB	       
	        /* use Netscape remote-controlled by X-atoms */
		graphWebBrowser (messprintf ("%s%s", URL, msp->sname));
#else
		if (!browser && !(browser = getenv("BLIXEM_WWW_BROWSER"))) {
		    printf("Looking for WWW browsers ...\n");
		    if (!findCommand("netscape", &browser) &&
			!findCommand("Netscape", &browser) &&
			!findCommand("Mosaic", &browser) &&
			!findCommand("mosaic", &browser) &&
			!findCommand("xmosaic", &browser)) {
			messout("Couldn't find any WWW browser.  Looked for "
				"netscape, Netscape, Mosaic, xmosaic & mosaic. "
				"System message: \"%s\"", browser);
			return;
		    }
		}
		printf("Using WWW browser %s\n", browser);
		fflush(stdout);
		system(messprintf("%s %s%s&", browser, URL, msp->sname));
#endif
	    }
#ifdef ACEDB	       
	    else if (!strcmp(fetchMode, "acedb")) {
		display(msp->key, 0, 0);
	    }
	    else if (!strcmp(fetchMode, "acedb text")) {
		display(msp->key, 0, TREE);
	    }
#endif
	    else 
		messout("Unknown fetchMode: %s", fetchMode);
	}
	return;
    }
    else
    {
        /* Reset all highlighted boxes */
	if (!squash)
	    for (msp = DBhits.next; msp ; msp = msp->next)
		if (msp->box && !strcmp(msp->sname, HighlightSeq)) 
		    graphBoxDraw(msp->box, BLACK, backgColor);
	
	if (BigPictON)
	    for (bpmsp = BPMSPlist.next; bpmsp && *bpmsp->sname; bpmsp = bpmsp->next)
		if (!strcmp(bpmsp->sname, HighlightSeq))
		    graphBoxDraw(bpmsp->box, BLACK, backgColor);


	/* Find clicked protein 
	 ***********************/

	for (msp = DBhits.next; msp ; msp = msp->next)
	    if (msp->box == box) break;

	if (msp) {
	    /* Picked box in alignment */
	    strcpy(HighlightSeq, msp->sname);
	    pickedMSP = msp;
	}
	else if (BigPictON)
        {
	    /* look for picked box in BigPicture */
	    for (bpmsp = BPMSPlist.next; bpmsp && *bpmsp->sname; bpmsp = bpmsp->next)
		if (bpmsp->box == box) break;
	    if (bpmsp && *bpmsp->sname) strcpy(HighlightSeq, bpmsp->sname);
	}
	
	if (msp || (bpmsp && *bpmsp->sname)) {
	    /* Put description in message box */
	    messagebox = graphBoxStart();
	    graphTextPtr (message, NAMESIZE+13, separator_y - .1, displen+17);
	    graphBoxEnd ();
	    if (msp) {
		strncpy(message, msp->sname, 1023);
		if (msp->desc) {
		    strcat(message, " ");
		    strncat(message, msp->desc, 1023-strlen(message));
		}
	    }
	    else {
		strncpy(message, bpmsp->sname, 1023);
		if (bpmsp->desc) {
		    strcat(message, " ");
		    strncat(message, bpmsp->desc, 1023-strlen(message));
		}
	    }
	    message[displen+17]=0;
	    graphBoxDraw(messagebox, BLACK, backgColor);

	    /* Highlight picked box */
	    highlightProteinboxes();
	    lastbox = box;
	}
    }
}

static void mspcpy(MSP *dest, MSP *src)
{
    dest->score  = src->score;
    dest->id     = src->id;
    strcpy(dest->frame, src->frame);
    dest->qstart = src->qstart; 
    dest->qend   = src->qend;
    strcpy(dest->sname, src->sname);
    dest->desc   = src->desc;
    dest->sstart = src->sstart; 
    dest->send   = src->send;
    dest->sseq   = src->sseq;
    dest->box    = src->box;
#ifdef ACEDB
    dest->key    = src->key;
#endif
}


static void sortMSPs(int (*func)(MSP*,MSP*))
{
    MSP tmpmsp, *msp1, *msp2;

    if (!DBhits.next) return;

    for (msp1 = DBhits.next ; msp1->next ; msp1 = msp1->next )
    {
	for (msp2 = msp1->next ; msp2 ; msp2 = msp2->next )
	{
	    if ( (*func)(msp1, msp2) > 0 )
	    {
		mspcpy(&tmpmsp, msp2);
		mspcpy(msp2, msp1);
		mspcpy(msp1, &tmpmsp);
	    }
	}
    }
    if (graphActivate(g)) blviewRedraw();
}


static int possort(MSP *msp1, MSP *msp2) {
    return ( (msp1->qstart > msp2->qstart) ? 1 : -1 );
}
static int namesort(MSP *msp1, MSP *msp2) {
    if (strcmp(msp1->sname, msp2->sname))
	return strcmp(msp1->sname, msp2->sname);
    else
	return possort(msp1, msp2);
}
static int scoresort(MSP *msp1, MSP *msp2) {
    return ( (msp1->score < msp2->score) ? 1 : -1 );
}
static int idsort(MSP *msp1, MSP *msp2) {
    return ( (msp1->id < msp2->id) ? 1 : -1 );
}


static void fetchByefetch(void) {
    strcpy(fetchMode, "efetch");
    if (graphActivate(g)) blviewRedraw();
}
static void fetchByWWWefetch(void) {
    strcpy(fetchMode, "WWW-efetch");
    if (graphActivate(g)) blviewRedraw();
}
#ifdef ACEDB
static void fetchByacedb(void) {
    strcpy(fetchMode, "acedb");
    if (graphActivate(g)) blviewRedraw();
}
static void fetchByacedbtext(void) {
    strcpy(fetchMode, "acedb text");
    if (graphActivate(g)) blviewRedraw();
}
#endif
static void sortByName(void) {
    sortMode = SORTBYNAME;
    strcpy(sortModeStr, "Name");
    sortMSPs(namesort);
}
static void sortByScore(void) {
    sortMode = SORTBYSCORE;
    strcpy(sortModeStr, "Score");
    sortMSPs(scoresort);
}
static void sortByPos(void) {
    sortMode = SORTBYPOS;
    strcpy(sortModeStr, "Position");
    sortMSPs(possort);
}
static void sortById(void) {
    sortMode = SORTBYID;
    strcpy(sortModeStr, "Identity");
    sortMSPs(idsort);
}

#ifdef CODE_NOT_USED
static void incBack(void) {
    backgColor++;
    if (!(backgColor % BLACK)) backgColor++;
    blviewRedraw();
}
static void decBack(void) {
    backgColor--;
    if (!(backgColor % BLACK)) backgColor--;
    blviewRedraw();
}
static void incGrid(void) {
    gridColor++;
    blviewRedraw();
}
static void decGrid(void) {
    gridColor--;
    blviewRedraw();
}
#endif /* CODE_NOT_USED */

static void squashMatches(void)
{
    static int oldSortMode;

    if (!squash) {
	oldSortMode = sortMode;
	sortByName();
	squash = 1;
    }
    else {
	switch (oldSortMode) {
	case SORTBYNAME : sortByName(); break;
	case SORTBYSCORE : sortByScore(); break;
	case SORTBYPOS : sortByPos(); break;
	case SORTBYID : sortById(); break;
	}
	squash = 0;
    }
    blviewRedraw();
}


static void blviewHelp(void)
{
    graphMessage (messprintf("\
\
BLIXEM - BLast matches\n\
         In an\n\
         X-windows\n\
         Embedded\n\
         Multiple alignment\n\
\n\
LEFT MOUSE BUTTON:\n\
  Pick on boxes and sequences.\n\
  Fetch annotation by double clicking on sequence (Requires 'efetch' to be installed.)\n\
\n\
MIDDLE MOUSE BUTTON:\n\
  Scroll horizontally.\n\
\n\
RIGHT MOUSE BUTTON:\n\
  Menu.  Note that the buttons Settings and Goto have their own menus.\n\
\n\
RESIDUE COLOURS:\n\
  Yellow = Query.\n\
  See Settings Panel for matching residues (click on Settings button).\n\
\n\
version %s\n\
(c) Erik Sonnhammer", blixemVersion));
}


static void wholePrint(void)
{
    int 
	tmp,
	dispstart_save = dispstart, 
	BigPictON_save = BigPictON;

    static int
	start=1, end=0;

    if (!end) end = qlen;

    /* Swap coords if necessary */
    if ((plusmin < 0 && start < end) || (plusmin > 0 && start > end )) {
	tmp = start;
	start = end;
	end = tmp;
    }

    /* Apply max limit MAXALIGNLEN */
    if ((abs(start-end)+1) > MAXALIGNLEN*symbfact) {
	start = dispstart - plusmin*MAXALIGNLEN*symbfact;
	if (start < 1) start = 1;
	if (start > qlen) start = qlen;

	end = start + plusmin*MAXALIGNLEN*symbfact;
	if (end > qlen) end = qlen;
	if (end < 1) end = 1;
    }

    if (!graphPrompt ("Please state the zone you wish to print", 
		      messprintf("%d %d", start, end),
		      "iiz")) 
	return;
    freeint(&start);
    freeint(&end);

    dispstart = start;
    displen = abs(end-start)+1;

    /* Validation */
    if (plusmin > 0 && start > end) {
	messout("Please give a range where from: is less than to:");
	return;
    }
    else if (plusmin < 0 && start < end) {
	messout("Please give a range where from: is more than to:");
	return;
    }
    if (displen/symbfact > MAXALIGNLEN) {
	messout("Sorry, can't print more than %d residues.  Anyway, think of the paper!", MAXALIGNLEN*symbfact);
	return;
    }

    wholePrintOn = 1;
    BigPictON = 0;

    blviewRedraw();
    graphPrint();

    /* Restore */
    wholePrintOn = 0;
    dispstart = dispstart_save;
    displen = dispstart_save;
    BigPictON = BigPictON_save;

    blviewRedraw();
}


/* BLXVIEW is called from some external program
 *
 * Interface: The string opts may contain a number of options that tell
 * blixem to start up with different defaults
 */
Graph blxview(char *seq, char *seqname, int start, int offset, MSP *msp, char *opts)
{
    char *opt;

    if (graphActivate(g)) blviewDestroy();
    graphPop();

    q = seq;
    qlen = actend = strlen(q);

    qname = seqname;

    dispstart = start;
    actstart=1;
    qoffset = offset;
    DBhits.next = msp;
    BPMSPlist.next = 0;
    BigPictZoom = 10;
    *HighlightSeq = 0;
    blastx = blastp =  blastn = tblastn = tblastx = 0;

    opt = opts;
    while (*opt) {

	/* Used options:

	   BGLMNPRSTXZ-+brsinp

	*/

	switch (*opt) {
	case 'G':		
	    /* Gapped HSPs */
	    HiliteSins = HSPgaps = 1;                break;
	case 'P':
	    blastp = 1;
	    blastx = blastn = tblastn = tblastx = 0;
	    alphabetsize = 24;
	    symbfact = 1;                           break;
	case 'N':
	    blastn = 1;
	    blastp = blastx = tblastn = tblastx = 0;
	    alphabetsize = 4;
	    symbfact = 1;                  
	    if (!strchr(opts, 'r')) BigPictRev = 1; 
                                                    break;
	case 'X':
	    blastx = 1;
	    blastp = blastn = tblastn = tblastx = 0;
	    alphabetsize = 4;
	    symbfact = 3;                           break;
	case 'T':
	    tblastn = 1;
	    blastp = blastx = blastn = tblastx = 0;
	    alphabetsize = 24;
	    symbfact = 1;                           break;
	case 'L':
	    tblastx = 1;
	    blastp = blastx = blastn = tblastn = 0;
	    alphabetsize = 24;
	    symbfact = 3;                           break;
	case '-':
	    strcpy(actframe, "(-1)"); 
	    plusmin = -1;                           break;
	case '+':
	    strcpy(actframe, "(+1)"); 
	    plusmin = 1;                            break;
	case 'B': 
	    BigPictON = 1;                          break;
	case 'b': 
	    BigPictON = 0;                          break;
	case 'd':
	    dotter_first = 1;                       break;
	case 'i':
	    for (msp = DBhits.next; msp; msp = msp->next)
		if (!msp->id) calcID(msp);
	    sortById();                             break;
	case 'M':
	    start_nextMatch = 1;                    break;
	case 'n':
	    sortByName();                           break;
	case 'p':
	    sortByPos();                            break;
	case 'R': 
	    BigPictRev = 1;                         break;
	case 'r': 
	    BigPictRev = 0;                         break;
	case 'S':
	    SEGon = 1;                              break;
	case 's':
	    sortByScore();                          break;
	case 'Z':
	    BigPictZoom = strlen(seq);              break;
	}

	opt++;
    }

    if (blastx + blastn + blastp + tblastn + tblastx == 0) {
	messout("Blast type not set - assuming blastx!");
	blastx = 1;
    }

#ifdef ACEDB
    if ((URL = getenv("BLIXEM_FETCH_WWW")))
	strcpy(fetchMode, "WWW-efetch");
    else if (getenv("BLIXEM_FETCH_EFETCH"))
	strcpy(fetchMode, "efetch");
    else 
	strcpy(fetchMode, "acedb");
#else
    if ((URL = getenv("BLIXEM_FETCH_WWW")))
	strcpy(fetchMode, "WWW-efetch");
    else if (getenv("BLIXEM_FETCH_EFETCH"))
	strcpy(fetchMode, "efetch");
    else {
	strcpy(fetchMode, "WWW-efetch");
    }
#endif

    if (!URL) {
	URL = messalloc(256);
	strcpy(URL, "http://www.sanger.ac.uk/cgi-bin/seq-query?");
    }

    return blviewCreate();
}    

#define autoDotterParamsStr "Automatic Dotter parameters"

/* BLVIEWCREATE initializes the display and the buttons
*/
static Graph blviewCreate (void)
{
  static MENUOPT mainMenu[] = 
  {   
    { graphDestroy,     "Quit" },
    { blviewHelp,       "Help" },
    { graphPrint,       "Print" },
    { wholePrint,       "Print whole alignment" },
    { blixemSettings,   "Change Settings" },
    { menuSpacer,       "" },
/*  { dotterPanel,       "Dotter panel" },*/
    { callDotter,       "Dotter" },
    { callDotterHSPs,   "Dotter HSPs only" },
    { callDotterSelf,   "Dotter query vs. itself" },
    { setDotterParams,  "Manual Dotter parameters" },
    { autoDotterParams, autoDotterParamsStr },
    { menuSpacer,       "" },
    { hidePicked,	  "Hide picked match" },
    { setHighlight,     "Highlight sequences by name" },
    { clrHighlight,     "Clear highlighted and unhide" },
    { 0, 0 }
  };

  static MENUOPT settingsMenuOpt[] = 
  {   
    { sortByScore,      SortByScoreStr },
    { sortById,         SortByIdStr },
    { sortByName,       SortByNameStr },
    { sortByPos,        SortByPosStr },
    { menuSpacer,       "" },
    { BigPictToggle,    BigPictToggleStr },
    { BigPictToggleRev, BigPictToggleRevStr },
    { SEGtoggle,        SEGtoggleStr },
    { menuSpacer,       "" },
    { toggleDESC,       toggleDESCStr },
    { toggleHiliteSins, toggleHiliteSinsStr },
    { squashMatches,    squashMatchesStr },
    { toggleIDdots,     toggleIDdotsStr },
    { toggleHiliteUpper,toggleHiliteUpperStr },
    { toggleHiliteLower,toggleHiliteLowerStr },
    { menuSpacer,       "" },
    { printColors,      printColorsStr },
    { toggleColors,     toggleColorsStr },
    { toggleVerbose,    toggleVerboseStr },
    { menuSpacer,       "" },
    { blixemSettings,   "Settings window" },
    { 0, 0 }
  };
    
    if (!graphActivate(g)) 
    {
	float w, h;
	if (!oldGraphWidth) {
	    graphScreenSize(&w, &h, 0, 0, 0, 0);
	    h *= 0.4;
	}
	else { 
	    w = oldGraphWidth;
	    h = oldGraphHeight;
	}
	g = graphCreate(TEXT_SCROLL, messprintf("Blixem %s", qname), 0, 0, w, h);

	graphRegister (DESTROY, blviewDestroy ); 
	graphRegister (RESIZE, blviewRedraw ); 
	graphRegister (PICK, blviewPick) ;
	graphRegister (KEYBOARD, keyboard) ;
	graphRegister (MIDDLE_DOWN, MiddleDown) ;
	graphRegister (KEYBOARD, keyboard) ;
#if defined(MACINTOSH)
	graphHelp ( "Blixem" ) ; /* Need Help for title on mac */
	graphMenu ( mainMenu ) ;
#else
	settingsMenu = menuInitialise ("Settings", (MENUSPEC*)settingsMenuOpt) ;
	if (blastp || tblastn) menuSetFlags(menuItem(settingsMenu, BigPictToggleRevStr), MENUFLAG_HIDE);
	graphNewMenu(settingsMenu);

	blixemMenu = menuInitialise ("blixem", (MENUSPEC*)mainMenu) ;
	if (!HSPgaps) menuSetFlags(menuItem(settingsMenu, toggleHiliteSinsStr), MENUFLAG_HIDE);
	menuSetFlags(menuItem(blixemMenu, autoDotterParamsStr), MENUFLAG_DISABLED);
	graphNewMenu(blixemMenu);

#endif
    }

    if ((cp = (char *)strrchr(qname, '/'))) qname = cp+1;
    if (strlen(qname) > NAMESPACE+4) qname[NAMESPACE+4] = 0;

    allocAuxseqs(INITDBSEQLEN);

    if (SEGon) calcSEGarrays(FALSE);

#ifdef ACEDB		
    /* Acedb passes revcomp exons incorrectly for blastn * /
    if (blastn) {
	messout ( "\nWarning: exons on the opposite strand may be invisible due "
	       "to a strange bug in blastn mode in acedb which inverts them, "

	       "although the identical routine works correctly in blastx mode.\n"
	       "Also, when Blixem is called from REVCOMP blastn mode, the HSPs "
	       "are not revcomped.\n");
    }*/
#endif

    if (dotter_first && DBhits.next && DBhits.next->sname[0]) {
	strcpy(HighlightSeq, DBhits.next->sname);
	callDotter();
    }

    if (start_nextMatch) 
	nextMatch();
    else 
	blviewRedraw();
    
    return FALSE /* g */ ;
}


static MENUOPT gotoMenu[] = 
{ 
  { gotoBegin, " <- Goto Beginning " },
  { gotoEnd,   "    Goto End ->    " },
  { 0, 0 }
};


/* BLVIEWDESTROY frees all the allocated memory
   N.B. This memory was allocated in the calling program (acedb)

   Could free auxseq, auxseq2 and padseq too, but they'd have to be remalloc'ed
   next time then.
*/
static void blviewDestroy(void)
{
    MSP *msp, *fmsp;
    BPMSP *tmsp;

    /* Free the allocated sequences */
    for (msp = DBhits.next; msp; msp = msp->next)
    {
	if (msp->sseq && msp->sseq != padseq) {
	    for (fmsp = msp->next; fmsp; fmsp = fmsp->next)
		if (fmsp->sseq == msp->sseq) fmsp->sseq = 0;

	    /* Bug in fmapfeatures.c causes introns to have stale sseq's */
	    if (msp->score >= 0) {
		messfree(msp->sseq);
		messfree(msp->desc);
	    }
	}
    }
    messfree(q);

    for (msp = DBhits.next; msp; )
    {
	fmsp = msp;
	msp = msp->next;
	messfree(fmsp);
    }

    for (bpmsp = BPMSPlist.next; bpmsp; )
    {
	tmsp = bpmsp;
	bpmsp = bpmsp->next;
	messfree(tmsp);
    }

    arrayDestroy(stringentSEGarr);
    arrayDestroy(mediumSEGarr);
    arrayDestroy(nonglobularSEGarr);
}


void drawBigPictMSP(MSP *msp, int BPx, char strand)
{
    float  msp_y, msp_sx, msp_ex, midx;
    float  oldLinew;

    msp_sx = max((float)plusmin*(msp->qstart-BigPictStart)*BPx/BigPictLen +BPoffset, 4);
    msp_ex = max((float)plusmin*(msp->qend-BigPictStart)*BPx/BigPictLen +BPoffset, 4);
    
    if (msp->score == -1) /* EXON */
    {
	oldcolor = graphColor(geneColor); oldLinew = graphLinewidth(.15);
	msp_y = 7.9 + queryy ;
	if (strand == 'R') msp_y += 1.5 ;
	graphRectangle(msp_sx, msp_y, msp_ex, msp_y + 0.7);
	graphColor(oldcolor); graphLinewidth(oldLinew);
    }
    else if (msp->score == -2) /* INTRON */
    {
	oldcolor = graphColor(geneColor); oldLinew = graphLinewidth(.15);
	msp_y = 7.9 + queryy ;
	if (strand == 'R') msp_y += 1.5 ;
	midx = 0.5 * (msp_sx + msp_ex) ;
	graphLine (msp_sx, msp_y+0.4, midx, msp_y) ;
	graphLine (msp_ex, msp_y+0.4, midx, msp_y) ;
	graphColor(oldcolor); graphLinewidth(oldLinew);
    }
    else if (msp->score == -3) /* COLOURED  FEATURE SEGMENT - Need to make this much more flexible */
    {
	/* Colour stored in ->id */
	oldcolor = graphColor(msp->id); oldLinew = graphLinewidth(.1);
	msp_y = 7.9 + queryy ;
	if (strand == 'R') msp_y += 1.5 ;
	graphFillRectangle(msp_sx, msp_y, msp_ex, msp_y + 0.7);
	graphColor(BLACK);
	graphRectangle(msp_sx, msp_y, msp_ex, msp_y + 0.7);
	graphColor(oldcolor); graphLinewidth(oldLinew);
    }
    else if (msp->score > 0) /* BLAST MATCH */
    {
        if (!msp->sseq) getsseq(msp);
	if (!msp->id) calcID(msp);
	msp_y = (float)(140 - msp->id)/20 + queryy ;
	if (strand == 'R') msp_y += 9 ;
	if (!bpmsp->next) bpmsp->next = (BPMSP *)messalloc(sizeof(BPMSP));
	bpmsp = bpmsp->next;
	strcpy(bpmsp->sname, msp->sname);
	bpmsp->desc = msp->desc;
	oldLinew = graphLinewidth(0.1);
	bpmsp->box = graphBoxStart();
	graphFillRectangle(msp_sx, msp_y, msp_ex, msp_y+.2);  
	/* graphLine doesn't want to be picked */
	graphBoxEnd();	graphBoxDraw(bpmsp->box, BLACK, TRANSPARENT);
	graphLinewidth(oldLinew);
    }
}

/* Checks if the msp is supposed to be drawn given its frame and position */
void selectBigPictMSP(MSP *msp, int BPx, int BigPictStart, int BigPictLen)
{
    if (actframe[1] == msp->frame[1])
	drawBigPictMSP(msp, BPx, 'M');
    if (BigPictRev) 
	if (actframe[1] != msp->frame[1])
	    drawBigPictMSP(msp, BPx, 'R');
    
    return;

#ifdef OLD_CODE_TOO_COMPLICATED_AND_SEEMS_BUGGED__JUST_DRAW_ALL
    if (( ((actframe[1] == '+' && msp->frame[1] == '+')) && 
	  (msp->qstart <  BigPictStart && msp->qend >  BigPictStart || /* left || span */
	   msp->qstart >= BigPictStart && msp->qend <= BigPictStart+BigPictLen || /* middle */
	   msp->qstart <  BigPictStart+BigPictLen && msp->qend > BigPictStart+BigPictLen)) /* right || span */
	||
	(actframe[1] == '-' && msp->frame[1] == '-' && 
	 (msp->qend > BigPictStart && msp->qstart < BigPictStart ||
	  msp->qend <= BigPictStart && msp->qstart >= BigPictStart-BigPictLen ||
	  msp->qend > BigPictStart-BigPictLen && msp->qstart < BigPictStart-BigPictLen)))
	drawBigPictMSP(msp, BPx, 'M');
    if (BigPictRev) 
	if ( (strchr(actframe, '+') && strchr(msp->frame, '-') && 
	      (msp->qend <  BigPictStart && msp->qstart < BigPictStart+BigPictLen || 
	       msp->qend >= BigPictStart && msp->qstart <= BigPictStart+BigPictLen ||
	       msp->qend <  BigPictStart+BigPictLen && msp->qstart > BigPictStart+BigPictLen))
	     ||
	     (strchr(actframe, '-') && strchr(msp->frame, '+') && 
	      (msp->qend >  BigPictStart && msp->qstart < BigPictStart ||
	       msp->qend <= BigPictStart && msp->qstart >= BigPictStart-BigPictLen ||
	       msp->qend > BigPictStart-BigPictLen && msp->qstart < BigPictStart-BigPictLen)))
	    drawBigPictMSP(msp, BPx, 'R');
#endif
} /* selectBigPictMSP */


void drawSEG(MSP *msp, int BPx, char strand)
{
    float  msp_y, msp_sx, msp_ex, y;

    msp_y = 0 ;
    msp_sx = max((float)plusmin*(msp->qstart-BigPictStart)*BPx/BigPictLen +BPoffset, 4);
    msp_ex = max((float)plusmin*(msp->qend-BigPictStart)*BPx/BigPictLen +BPoffset, 4);

    y = queryy + 8.5;
    
    if (msp->score == -4) {
	msp_y = y+1;
	oldcolor = graphColor(stringentSEGcolor); 
    }
    else if (msp->score == -5) {
	msp_y = y + .5;
	oldcolor = graphColor(mediumSEGcolor); 
    }
    else if (msp->score == -6) {
	msp_y = y;
	oldcolor = graphColor(nonglobularSEGcolor); 
    }

    if (strand == 'R') msp_y += 1.5 ;
    graphFillRectangle(msp_sx, msp_y, msp_ex, msp_y + 0.3);
    graphColor(oldcolor);
}


static void drawSEGcurve(int start, int end, Array array, int win)
{
    int i;
    float x, y, xold=0, yold=0, oldLinew = graphLinewidth(.3) /*, min=4, max=0 */;

    for (i = start; i < end; i++) {
	if (i > win/2 && i < qlen - win/2) {
	    x = SEQ2BP(i);
	    y = queryy + 9 - arr(array, i, float)*2;
	    if (xold) graphLine(xold, yold, x, y);
	    xold = x;
	    yold = y;

	    /*if (arr(array, i, float) < min) min = arr(array, i, float);
	    if (arr(array, i, float) > max) max = arr(array, i, float);*/
	}
	/**c = q[i];
	  graphText(c, SEQ2BP(i), queryy+5);*/
    }
    /* printf("min = %f , max = %f\n", min, max); */

    graphLinewidth(oldLinew);
}


/* BLVIEWREDRAW Redraws the view in the disp_region
*/
static void blviewRedraw(void)
{
    char    text[MAXALIGNLEN+1], *query=0, querysym, subjectsym, *lastname=0;
    int     i, 
            dbdisp,  /* number of msp's displayed */
            frame=0,
	    curframe, DNAline,
            TickStart, TickUnit, gridx,
            pos, msp_offset, exon_end, firstRes, slen=0;
    MSP    *msp;
    float   oldLinew, alignstart;

    settingsRedraw();

    graphActivate(g);
    graphClear();
    graphBoxDraw(0, backgColor, backgColor);
    graphTextFormat(FIXED_WIDTH);

    for (msp = DBhits.next; msp ; msp = msp->next) msp->box = 0;
    lastbox = messagebox = 0;
    queryy = separator_y = 0.2;


/* Calculate window sizes
**************************/

    graphFitBounds (&nx, &ny);
    graphWindowSize (0, 0, &oldGraphWidth, &oldGraphHeight);

    graphColor(backgColor); graphRectangle(0, 0, nx+100, ny+100);graphColor(BLACK);

    if (wholePrintOn) {
	nx = displen/symbfact + 43;
    }
    else {
	displen = (nx - 43)*symbfact ;        /* Keep spacious for printing */
    }

	

    /* Window boundary check */
    if (strchr(actframe, '+'))
    { 
	if (displen > qlen - (symbfact+1) )                 /* Window too big for sequence */
	{ 
	    dispstart = 1;
	    displen = qlen - (symbfact+1);
	    if (blastp || blastn) displen = qlen;
	}
	else if (dispstart < 1)                             /* Off the left edge  */
	    dispstart = 1 ;
	else if (dispstart + displen > qlen - (blastx||tblastx ? 3 : -1))     /* Off the right edge */
	    dispstart = qlen - (blastx||tblastx ? 3 : -1) - displen;
    }
    else if (strchr(actframe, '-'))
    { 
	if (displen > qlen - (symbfact+1))                  /* Window too big for sequence */
	{ 
	    dispstart = qlen;
	    displen = qlen - (symbfact+1);
	    if (blastn) displen = qlen;
	}
        else if (dispstart > qlen)                          /* Off the left edge  */
	    dispstart = qlen ;
	else if (dispstart - displen < (blastx||tblastx ? 4 : 0) )   /* Off the right edge */
	    dispstart =  (blastx||tblastx ? 4 : 0) + displen ;
    }


/* 0. Draw Big Picture
**************************/

    if (BigPictON) 
    {
	BigPictLen = displen * BigPictZoom;
	if (BigPictLen > qlen) {
	    BigPictLen = qlen;
	    BigPictZoom = (float)qlen/displen;
	}
	BigPictStart = dispstart + plusmin*displen/2 - plusmin*BigPictLen/2;
	if (plusmin > 0) {
	    if (BigPictStart + BigPictLen > qlen) BigPictStart = qlen - BigPictLen;
	    if (BigPictStart < 1) BigPictStart = 1;
	}
	else {
	    if (BigPictStart - BigPictLen < 1) BigPictStart = BigPictLen;
	    if (BigPictStart > qlen) BigPictStart = qlen;
	}
	BPx = nx - BPoffset;
	graphButton("Zoom In", zoomIn, .5, queryy +.1);
	graphButton("Zoom Out", zoomOut, 9, queryy +.1);
	graphButton("Whole", zoomWhole, 18.5, queryy +.1);

	/* Draw position lines (ticks) */
	TickUnit = (int)pow(10.0, floor(log10((float)BigPictLen)));                       /* Basic scaleunit */
	while (TickUnit*5 > BigPictLen) TickUnit /= 2;                                  /* Final scaleunit */
	TickStart = (int)TickUnit*( floor((float)(BigPictStart+qoffset)/TickUnit + .5));/* Round off       */
	for (i=TickStart; 
	     (plusmin == 1 ? i < BigPictStart+qoffset + BigPictLen : i > BigPictStart+qoffset - BigPictLen); 
	     i = i + plusmin*TickUnit) 
	{
	    gridx = (float)plusmin*(i-BigPictStart-qoffset)*BPx/BigPictLen + BPoffset;
	    if (gridx > 24 && gridx < nx-10)
	    {
		graphText(messprintf("%d", i), gridx, queryy);
		graphColor(gridColor);
		graphLine(gridx, queryy+1, gridx, queryy+7);
		if (BigPictRev) graphLine(gridx, queryy+11, gridx, queryy+16);
		graphColor(BLACK);
	    }
	}
	
	/* Draw Percentage lines */
	for (i=0; i<=100; i += 20) {
	    graphText(messprintf("%3i%%", i), 0, queryy + (float)(130-i)/20);
	    if (BigPictRev) graphText(messprintf("%3i%%", i), 0, queryy+9 + (float)(130-i)/20);
	    graphColor(gridColor);
	    graphLine(4, queryy + (float)(140-i)/20, nx+1, queryy + (float)(140-i)/20);
	    if (BigPictRev) graphLine(4, queryy+9 + (float)(140-i)/20, nx+1, queryy+9 + (float)(140-i)/20);
	    graphColor(BLACK);
	}


	/* Draw active window frame */
	oldLinew = graphLinewidth(.3);
	BPboxstart = (float)plusmin*(dispstart-BigPictStart)*BPx/BigPictLen + BPoffset;
	BPboxwidth = (float)displen*BPx/BigPictLen;
	BPbox = graphBoxStart();
	graphRectangle(BPboxstart, queryy +1.7, BPboxstart+BPboxwidth, queryy +1.7 + BPboxheight);
	graphBoxEnd(); graphBoxDraw(BPbox, geneColor, TRANSPARENT);
	graphLinewidth(.5);
        graphLinewidth(oldLinew);
	

	/* Draw Big Picture MSPs */
	lastExonx = 0;
	bpmsp = &BPMSPlist;
	for (msp = DBhits.next; msp; msp = msp->next)
        {
	    if (strcmp(msp->sname, HighlightSeq) && !strMatch(msp->sname, searchSeq))
		selectBigPictMSP(msp, BPx, BigPictStart, BigPictStart+BigPictLen);
	}
	/* Draw the highlighted later, so they're not blocked by black ones */
	for (msp = DBhits.next; msp; msp = msp->next)
	{
	    if (!strcmp(msp->sname, HighlightSeq) || strMatch(msp->sname, searchSeq))
		selectBigPictMSP(msp, BPx, BigPictStart, BigPictStart+BigPictLen);
	}
	if (bpmsp->next) *bpmsp->next->sname = 0;

	queryy += 10;
	if (BigPictRev) queryy += 8;

	/* Draw separator line with perforation * /
	oldLinew = graphLinewidth(1);
	if (squash) graphColor(RED);
	else graphColor(BLACK);
	graphLine(0, queryy-.5, nx+1, queryy-.5);
	graphLinewidth(.1);
	graphColor(WHITE);
	for (i = 0; i <= nx+1; i++) graphLine((float)i, queryy-.5, (float)i+.8, queryy-.5);
	graphColor(BLACK);
	graphLinewidth(oldLinew); */

	/* Draw separator line a la Mosaic <HR> */
#define separatorwidth 0.5
	oldLinew = graphLinewidth(separatorwidth);
	if (squash) 
	    graphColor(RED);
	else 
	    graphColor(GRAY);
	graphLine(0, queryy-1.0, nx+1, queryy-1.0);

	graphLinewidth(.2);
	graphColor(WHITE);
	graphLine(0, queryy-1.1+0.5*separatorwidth, nx+1, queryy-1.1+0.5*separatorwidth);
	graphColor(BLACK);
	graphLinewidth(oldLinew);

	
	/* 1. Draw SEG features
	 **************************/

	if (SEGon) {
	    int y;
		    
	    /* Draw scale */
	    for (i = 1; i <= 4; i++) {
		y = queryy+5 - (i-2)*2;
		graphText(messprintf("%d", i), 1, y-.5); 

		graphColor(DARKGRAY);
		graphLine (BPoffset, y, nx, y);
		graphColor(BLACK);
	    }

	    /* Draw entropy plots */
	    graphColor(stringentSEGcolor);
	    drawSEGcurve(BigPictStart, BigPictStart+BigPictLen, 
			 stringentSEGarr, stringentSEGwin);

	    graphColor(mediumSEGcolor);
	    drawSEGcurve(BigPictStart, BigPictStart+BigPictLen, 
			 mediumSEGarr, mediumSEGwin);

	    graphColor(nonglobularSEGcolor);
	    drawSEGcurve(BigPictStart, BigPictStart+BigPictLen, 
			 nonglobularSEGarr, nonglobularSEGwin);


	    /* Draw precalculated segments */
	    for (msp = DBhits.next; msp; msp = msp->next) {
		if (msp->score < -3 && msp->score > -7) drawSEG(msp, BPx, 'M');
	    }

	    queryy += 11;

	    oldLinew = graphLinewidth(.1);
	    if (squash) 
		graphColor(RED);
	    else 
		graphColor(DARKGRAY);
	    graphLine(0, queryy-1, nx+1, queryy-1);
	    graphColor(WHITE);
	    graphLine(0, queryy-.9, nx+1, queryy-.9);
	    graphColor(BLACK);
	    graphLinewidth(oldLinew);
	    

	}

	separator_y = queryy;
    }



/* 1. Draw buttons & boxes
**************************/


    if (BigPictON) {
	menuUnsetFlags(menuItem(settingsMenu, BigPictToggleRevStr), MENUFLAG_DISABLED);
    }
    else if (blastn || blastx || tblastx) {
	menuSetFlags(menuItem(settingsMenu, BigPictToggleRevStr), MENUFLAG_DISABLED);
    }

    settingsButton = graphButton("Settings", blixemSettings, 0.5, queryy);
    graphBoxMenu( graphButton("  Goto  ", Goto, 0.5, queryy + 1.6), gotoMenu);
    graphButton("<match", prevMatch,  10.5, queryy);
    graphButton("match>", nextMatch,  10.5, queryy + 1.6);
    graphButton("<<", scrollLeftBig,  18, queryy);
    graphButton(">>", scrollRightBig, 18, queryy + 1.6);
    graphButton("<", scrollLeft1,     21.5, queryy);
    graphButton(">", scrollRight1,    21.5, queryy + 1.6);
    if (blastx ||tblastx || blastn)
	graphButton("Strand^v", ToggleStrand, 24.5, queryy + 1.6);
    graphButton("Help", blviewHelp,   nx-5, queryy + 1.6);

    /* Draw frame lines */
    graphLine(NAMESIZE +  1.0, separator_y+3, NAMESIZE +  1.0, ny+100);
    graphLine(NAMESIZE +  7.5, separator_y+3, NAMESIZE +  7.5, ny+100);
    graphLine(NAMESIZE + 11.5, separator_y+3, NAMESIZE + 11.5, ny+100);
    graphLine(NAMESIZE + 21.5, separator_y+3, NAMESIZE + 21.5, ny+100);
    graphLine(NAMESIZE + 22.5 + displen/symbfact, separator_y+3,
	      NAMESIZE + 22.5 + displen/symbfact, ny+100);

    queryy += 3;
    DNAline = queryy-2;

    graphText( "Score %ID  Start", 14, queryy);
    graphText( "End", NAMESIZE + 23 + displen/symbfact, queryy);
    queryy++;

    for (curframe = 1; 
	 abs(curframe) <= symbfact; 
	 curframe += (blastn ? -2 : 1))
    {
	frame = plusmin * curframe;
	sprintf(actframe, "(%+d)", frame);

	if (blastn && curframe == -1) compN = 1;
	else compN = 0;

	alignstart = NAMESIZE + 22 + (blastx || tblastx ? (float)(curframe - 1)/3 : 0);

/* 2. Draw the query
********************/
    
	if (blastx || tblastx)
	{ 
	    /* Force dispstart to cohere to the frame */
	    if (strchr(actframe, '-')) 
		while ((qlen - dispstart +1 +frame) % symbfact) dispstart--; 
	    else
		while ((dispstart - frame) % symbfact) dispstart++;
	}

	if (!(query = getqseq(dispstart, dispstart + plusmin*(displen-1), q))) return;
	if (compN) compl(query);

	if (!colortoggle) graphColor(hiColor);
	else graphColor(WHITE);
	graphFillRectangle(0, queryy, nx, queryy+1);
	graphColor(BLACK);
	graphLine(0, queryy, nx, queryy);

	graphText(query, alignstart, queryy);

	graphText (messprintf ("%-8s%s", qname, actframe),  1, queryy);
	graphText (messprintf ("%9d", dispstart + qoffset),  NAMESPACE+12, queryy);
	sprintf (text, "%-9d", dispstart+plusmin*displen + qoffset -plusmin); 
	graphText (text,  alignstart + 1 + displen/symbfact, queryy);

	if (blastx || tblastx)
	{
	    /* Draw the DNA sequence */
	    graphText(get3rd_base(dispstart, dispstart + plusmin*(displen-1), q), 
		      alignstart, DNAline++);
	}


/* 3. Draw the dbhits
*********************/

	dbdisp = 0;
	fetchCount = 0;
	for (msp = DBhits.next; msp; msp = msp->next)
        {
	    /* Current frame MSP - display if in window */
	    if (!strcmp(msp->frame, actframe) || (blastn && msp->frame[1] == actframe[1]))
	    {
		if (   (!compN && frame>0 && msp->qstart <= dispstart + displen && msp->qend >= dispstart)
		    || (!compN && frame<0 && msp->qstart >= dispstart - displen && msp->qend <= dispstart)
		    || ( compN && frame<0 && msp->qend   <= dispstart + displen && msp->qstart >= dispstart)
		    || ( compN && frame>0 && msp->qend   >= dispstart - displen && msp->qstart <= dispstart))
		{
		    if (msp->score == -1) 
		    {
			/* EXON */
			if (blastx || tblastx ||  (blastn && msp->frame[1] == actframe[1]) )
			{

			    dbdisp++;
			    if (!compN) 
			    {
				pos = (msp->qstart - dispstart) / symbfact ;
				exon_end = (msp->qend - dispstart) / symbfact ;
			    }
			    else 
			    {
				pos = dispstart - msp->qend;
				exon_end = dispstart - msp->qstart;
			    }
			    if (frame < 0) {
				pos = -pos;
				exon_end = -exon_end;
			    }
			    exon_end++;
			    if (pos < 0) pos = 0 ;
			    if (exon_end > displen/symbfact) exon_end = displen/symbfact;
			    
			    graphColor(geneColor); oldLinew = graphLinewidth(.3);
			    graphRectangle(alignstart+pos, queryy, 
					   alignstart+exon_end, queryy+1);
			    
			    graphColor(hiColor);
			    graphFillRectangle(alignstart+pos, queryy+dbdisp,
					       alignstart+exon_end, queryy+dbdisp+1);
			    graphColor(BLACK); 
			    strncpy(text, msp->sname, NAMESIZE);
			    text[NAMESIZE] = 0;
			    graphText(text, 1, queryy+dbdisp);
			    sprintf(text, "%5d  %3d %9d", 
				    msp->score, msp->id, (compN? msp->send : msp->sstart));
			    graphText(text,  NAMESIZE+1, queryy+dbdisp);
			    sprintf(text, "%-6d", (compN? msp->sstart : msp->send));
			    graphText(text,  alignstart+1+displen/symbfact, queryy+dbdisp);
			    /*			  graphText(msp->sseq, alignstart+pos, queryy+dbdisp); */

			    graphLinewidth(oldLinew);
		        }
		    }
		    else if (msp->score >= 0)
		    {  
			if (!msp->sseq) getsseq(msp);
			if (!msp->id) calcID(msp);
			
			if (!squash || !lastname || strcmp(msp->sname, lastname)) dbdisp++;
			lastname = msp->sname;
			
			msp->box = graphBoxStart () ;
			if ((cp = (char *)strchr(msp->sname, ':'))) cp++;
			else cp = msp->sname;
			strncpy(text, cp, NAMESIZE);
			text[NAMESIZE] = 0;
			graphText(text, 1, queryy+dbdisp);
			sprintf(text, "%5d  %3d %9d", 
				msp->score, msp->id, (compN ? msp->send : msp->sstart));
			graphText(text,  NAMESIZE+1, queryy+dbdisp);
			sprintf(text, "%-6d", (compN ? msp->sstart : msp->send)); 
			graphText(text,  alignstart+1+displen/symbfact, queryy+dbdisp);
			
			firstRes = -1;

			if (compN) {
			    if (HSPgaps)
				slen = strlen(msp->sseq);
			    else
				slen = msp->send;
			}

			for (i=0; i< displen/symbfact; i++)
			{
			    text[i] = ' ';
			    if (   (!compN && frame>0 && msp->qstart <= dispstart + i*symbfact 
				    && msp->qend >= dispstart + i*symbfact)
				   || (!compN && frame<0 && msp->qstart >= dispstart - i*symbfact 
				       && msp->qend <= dispstart - i*symbfact)
				   
				   || (compN && frame<0 && msp->qend <= dispstart + i
				       && msp->qstart >= dispstart + i)
				   || (compN && frame>0 && msp->qend >= dispstart - i
				       && msp->qstart <= dispstart - i))
			    {
				if (firstRes == -1) firstRes = i;
				
				msp_offset = plusmin*(int)(dispstart - msp->qstart)/symbfact + i;
				
				if (tblastn || tblastx) {
				    /* Fudged for the moment */
				    msp_offset++;
				    /* msp_offset += (msp->sstart-1)/3 + 1; */
				}
				else 
				    msp_offset += msp->sstart;
				
				if (compN) {
				    msp_offset = slen - plusmin*(dispstart - msp->qend) - i;
				}
				subjectsym = msp->sseq[msp_offset-1];
				querysym   = query[i];
				
				if (ace_upper(querysym) == ace_upper(subjectsym))
				{
				    if (IDdots) subjectsym = '.';
				    if (!colortoggle && !IDdots) graphColor(IDcolor);
					else graphColor(backgColor);
				}
				else if (!blastn && PAM120[ aa_atob[(int)querysym]-1 ][ aa_atob[(int)subjectsym]-1 ] > 0)
				{
				    if (!colortoggle) graphColor(consColor);
				    else graphColor(backgColor);
				}
				else
				{ 
				    /* Mismatch */
				    if (IDdots && !colortoggle) {
					if (subjectsym == '.') subjectsym = '-';
					graphColor(IDcolor);
				    }
				    else 
					graphColor(backgColor);
				}
				
				if (HiliteSins && 
				    ((isdigit(subjectsym) || subjectsym == '(' || subjectsym == ')') ||
				     (isalpha(subjectsym) && subjectsym == ace_lower(subjectsym))))
				    graphColor(RED);
				
				if (HiliteUpperOn && isupper(subjectsym)) graphColor(RED);
				if (HiliteLowerOn && islower(subjectsym)) graphColor(RED);
				
				text[i] = subjectsym;
				graphFillRectangle(i+alignstart, queryy+dbdisp, i+alignstart+1, queryy+dbdisp+1);
			    }
			}
			text[i] = '\0';
			graphColor(BLACK);
			graphText(text, alignstart, queryy+dbdisp);
			
			graphBoxEnd () ;
			graphBoxDraw(msp->box, BLACK, TRANSPARENT);
			if (squash) {
			    graphColor(RED);
			    graphLine(alignstart+firstRes, queryy+dbdisp, 
				      alignstart+firstRes, queryy+dbdisp+1);
			    graphColor(BLACK);
			}
		    }
		    if (DESCon && msp->desc) {
			/* Draw DESC line */
			int box;
			
			dbdisp++;
			box = graphBoxStart();
			graphText(messprintf("%s %s", msp->sname, msp->desc), 1, queryy+dbdisp);
			graphBoxEnd();
			graphBoxDraw(box, BLACK, WHITE);
		    }

		}

	    } /* else not in active frame */
	    
	} /* End of HSPs */
	queryy += dbdisp + 2;

    } /* End of frames */
    
    if (blastn) {
	frame = -frame;
	sprintf(actframe, "(%+d)", frame);
    }
    
    if (blastx || tblastx) dispstart -= plusmin*2;
    
    /* Draw frame lines */

    /* Draw centre notch */
    graphLine(NAMESIZE + 22.5 + displen/symbfact/2, separator_y,
	      NAMESIZE + 22.5 + displen/symbfact/2, separator_y+1);
    

    graphTextBounds(nx, queryy);

    if (!wholePrintOn) graphRedraw();
    highlightProteinboxes();

    messfree(query);

    setMenuCheckmarks();
}


/* GET3RD_BASE returns every third base from string q
*/
static char *get3rd_base(int start, int end, char *q)
{
  static int  
      i, len;
  static char 
      *bases = 0,
      *aux,
      *aux2;

  if (start < 1 || end > qlen) {
      messerror ( "Genomic sequence coords out of range: %d - %d\n", 
		 start, end);
      return NULL;
  }

  if (!bases) {
      bases = messalloc(qlen+1);
      aux = messalloc(qlen+1);
      aux2 = messalloc(qlen+1);
  }

  len = abs(end-start) + 1;

  if (start < end)
    for (i=start; i < end; i += 3 ) 
      bases[(i-start)/3] = q[i-1];
  else if (start > end) /* Reverse and complement it first */
    {
      strncpy(aux2, q+end-1, start-end+1);
      aux2[start-end+1] = 0;
      if (!revcomp(aux, aux2))
	  messcrash ("Cannot reverse & complement") ;
	
      for (i=0; i < len; i += 3) bases[i/3] = aux[i];
    }
  else return NULL;

  bases[len/3] = '\0';

  return bases;
}


/* CALCID caculated percent identity of an MSP
*/
static void calcID(MSP *msp)
{
    static int id, i, len;
    char *qseq;
    
    if (msp->score >= 0 && msp->sseq != padseq) { 
	if (!(qseq = getqseq(msp->qstart, msp->qend, q))) { 
	    messdump ( "calcID failed: Don't have genomic sequence %d - %d\n", 
		    msp->qstart, msp->qend);
	    msp->id = 0;
	    return;
	}
	
	if (tblastn) {
	    len = msp->qend - msp->qstart + 1;

	    for (i=0, id=0; i < len; i++)
		if (ace_upper(msp->sseq[i]) == ace_upper(qseq[i])) id++;
	}
	else if (tblastx) {
	    len = abs(msp->qend - msp->qstart + 1)/3;

	    for (i=0, id=0; i < len; i++)
		if (ace_upper(msp->sseq[i]) == ace_upper(qseq[i])) id++;
	}
	else 
	{
	    len = msp->send - msp->sstart + 1;
	    
	    for (i=0, id=0; i< len; i++)
		if (ace_upper(msp->sseq[i + msp->sstart -1]) == ace_upper(qseq[i])) id++;
	}
	messfree(qseq);
	msp->id = (int)((float)100*id/len + .5);
    }
    else msp->id = 0 ;
}    


static void assignPadseq(MSP *msp)
{	
    static int padseqlen=0;
    char *oldpadseq;
    int len = max(msp->sstart, msp->send);
    MSP *hsp;
    
    if (!padseq) {
	padseq = messalloc(INITDBSEQLEN+1);
	memset(padseq, '-', INITDBSEQLEN); 
	padseqlen = INITDBSEQLEN;
    }

    if (len > padseqlen) {
	oldpadseq = padseq;
	messfree(padseq);

	padseq = messalloc(len+1);
	memset(padseq, '-', len);
	padseqlen = len;

	/* Change all old padseqs to new */
	for (hsp = DBhits.next; hsp ; hsp = hsp->next)
	    if (hsp->sseq == oldpadseq) hsp->sseq = padseq;
    }

    msp->sseq = padseq;
}



static void allocAuxseqs(int len)
{
    if (auxseq) {
	messfree(auxseq);
	messfree(auxseq2);
    }

    auxseq = messalloc(len+1);
    auxseq2 = messalloc(len+1);
    auxseqlen = len;
}
	    

/* GETQSEQ translates a segment of the query seq (with 'sequence' coords = 1...)
*/
static char *getqseq(int start, int end, char *q)
{
    char *query;

    if (start < 1 || end < 1 || start > strlen(q) || end > strlen(q) || end==start)
    {
	messdump ( "Requested query sequence %d - %d out of available range: 1 - %d\n", 
	start, end, strlen(q));
	return NULL;
    }

    if (blastp || tblastn)
    {
	query = messalloc(end-start+2);
	strncpy(query, q+start-1, end-start+1);
	query[end-start+1] = 0;
	return query;
    }

    if (abs(end-start)+1 > auxseqlen) allocAuxseqs(abs(end-start)+1);

    if (start < end)
    {
	strncpy(auxseq, q+start-1, end-start+1);
	auxseq[end-start+1] = 0;
    }
    else if (start > end) /* Reverse and complement it */
    {
	strncpy(auxseq2, q+end-1, start-end+1);
	auxseq2[start-end+1] = 0;
	if (!revcomp(auxseq, auxseq2))
	    messcrash ("Cannot reverse & complement") ;
    }
    else return NULL;

    if (blastn) 
    {
	query =  messalloc(strlen(auxseq)+1);
	strcpy(query, auxseq);
	return query;
    }

    /* Translate the DNA sequence */
    if (!(query = translate(auxseq, stdcode1)))
	messcrash ("Cannot translate the genomic sequence") ;

    return query;
}


/* GETSSEQ fetches the database sequence from an external database
*/
static void getsseq(MSP *msp)
{
#if !(defined(MACINTOSH) || defined(WIN32))
    static char 
	  fetchstr[256];
    char *cp;
    FILE *pipe;
    MSP  *auxmsp;
    int   len, c;

    if (!*msp->sname) {
	messout ( "Nameless HSP at %d-%d - skipping Efetch\n", 
		msp->qstart+qoffset, msp->qend+qoffset);
	assignPadseq(msp);
	return;
    }

    if (verbose)
	printf("Efetching %s\n", msp->sname);
    else {
	if (!fetchCount) {
	    fetchCount++;
	    printf("\nEfetching external sequences");
	}
	printf("."); fflush(stdout);
    }

    sprintf(fetchstr, "efetch -q '%s'", msp->sname);

    if (msp->score >= 0)
    { 
	len = 0;
	cp = auxseq;

	if (!(pipe = (FILE*)popen(fetchstr, "r")))
	    messcrash("Failed to open pipe %s\n", fetchstr);
	
	while (!feof(pipe))
	    if (isalpha(c=fgetc(pipe)) && len < auxseqlen) {
	        *cp++ = c;
		len++;
	    }
	pclose(pipe);
	
	/* Check if auxseq was long enough */
	if (len > auxseqlen) {
	    allocAuxseqs(len);
	    cp = auxseq;
	    
	    if (!(pipe = (FILE*)popen(fetchstr, "r")))
		messcrash("Failed to open pipe %s\n", fetchstr);

	    while (!feof(pipe))
	      {
		c=fgetc(pipe) ; /* need 2 lines. is a macro on mac_x */
		if (isalpha(c))
		  *cp++ = c;
	      }
	    pclose(pipe);
	}
	
	*cp = 0;

	if (len > 1) /* Otherwise failed */
	{
	    if ((cp = (char *)strchr(auxseq, '\n'))) *cp = 0;
	    msp->sseq = messalloc(strlen(auxseq)+1);
	    
	    /* Harmonize upper and lower case - (t)blastx always translate
	     * the query to upper case */
	    if (isupper(*q) || tblastx || blastx) {
		for (i=0; auxseq[i]; i++)
		    msp->sseq[i] = ace_upper(auxseq[i]);
		msp->sseq[i] = 0;
	    }
	    else  {
		for (i=0; auxseq[i]; i++)
		    msp->sseq[i] = ace_lower(auxseq[i]);
		msp->sseq[i] = 0;
	    }

	    /* Check illegal offsets */
	    len = strlen(msp->sseq);
	    for (auxmsp = DBhits.next; auxmsp ; auxmsp = auxmsp->next)
		if (!strcmp(auxmsp->sname, msp->sname) && auxmsp->send > len )
		{
		    printf("%s HSP with offset beyond sequence (%d > %d) - using pads\n", 
			   msp->sname, auxmsp->send, len);
		    assignPadseq(msp);
		    break;
		}
	}	    
	else 
	{
#if !defined(ACEDB)
	    messdump ( "Unable to efetch %s - using pads instead\n", msp->sname);
#endif
	    /* Sequence not in database - fill up with pads */
	    assignPadseq(msp);
	}

	/* Set sseq for all MSPs of this subject */
	for (auxmsp = DBhits.next; auxmsp ; auxmsp = auxmsp->next)
	    if (!strcmp(auxmsp->sname, msp->sname))
	    {
		if (msp->sseq == padseq)
		    assignPadseq(auxmsp);
		else
		    auxmsp->sseq = msp->sseq;
		calcID(auxmsp);
	    }
    }
#else
	messout ( "WIN32/MAC Warning: efetch external command not implemented for blxview::getsseq()\n"
			  "Cannot retrieve sequence \"%s\" HSP at %d-%d - using pads instead\n", 
				msp->sname, msp->qstart+qoffset, msp->qend+qoffset);
	assignPadseq(msp);
#endif /* !(defined(MACINTOSH) || defined(WIN32))*/
}


static char *fetchSeqRaw(char *seqname)
{
#if !(defined(MACINTOSH) || defined(WIN32))
    static char fetchstr[256];
    char *cp, *retval;
    FILE *pipe;
    int len;

    if (!*seqname) {
	messdump ( "Nameless sequence - skipping Efetch\n");
	return 0;
    }
    printf("Efetching %s...\n", seqname);
    
    sprintf(fetchstr, "efetch -q '%s'", seqname);

    len = 0;
    cp = auxseq;

    if (!(pipe = (FILE*)popen(fetchstr, "r")))
	messcrash("Failed to open pipe %s\n", fetchstr);
	
    while (!feof(pipe)) {
	if (len < auxseqlen) 
	    *cp++ = fgetc(pipe);
	else 
	    fgetc(pipe);
	
	len++;
    }
    pclose(pipe);
	
    /* Check if auxseq was long enough */
    if (len > auxseqlen) {
	allocAuxseqs(len);
	cp = auxseq;

	if (!(pipe = (FILE*)popen(fetchstr, "r")))
	    messcrash("Failed to open pipe %s\n", fetchstr);
	
	while (!feof(pipe)) *cp++ = fgetc(pipe); 
	pclose(pipe);
    }

    *cp = 0;

    if (len > 1) /* Otherwise failed */
    {
	if ((cp = (char *)strchr(auxseq, '\n'))) *cp = 0;
	retval = messalloc(strlen(auxseq)+1);
	strcpy(retval, auxseq);
	
    }	    
    else {
	messdump ( "Unable to efetch %s \n", seqname);
	retval = 0;
    }

    return retval;
#else
	messout ( "WIN32/MAC Warning: efetch external command not implemented for blxview::fetchSeqRaw()\n") ;
	return 0 ;
#endif /* !(defined(MACINTOSH) || defined(WIN32)) */
}


/********************************************************************************
**                            BIG PICTURE ROUTINES                            ***
********************************************************************************/

static void BigPictToggle (void) { 
    BigPictON = !BigPictON; 
    blviewRedraw();
}

static void BigPictToggleRev (void) { 
    BigPictRev = !BigPictRev; 
    blviewRedraw();
}


static void zoomOut(void) {
    BigPictZoom *= 2;
    blviewRedraw();
}
static void zoomIn(void) {
    BigPictZoom /= (float)2;
    if (BigPictZoom < 1) BigPictZoom = 1;
    blviewRedraw();
}
static void zoomWhole(void) {
    BigPictZoom = (float)qlen/displen;
    blviewRedraw();
}


/* If crosshair-coordinates are screwed up, change here!
 ********************************************************/
static int x_to_residue(float x)
{
    int retval;

    if (blastx || tblastx) 
	retval = dispstart + plusmin*(x - NAMESIZE - 22.3)*3;
    else
	retval = dispstart + plusmin*(x - NAMESIZE - 22);


    if (plusmin > 0) {
	if (retval < dispstart) retval = dispstart;
	if (retval > dispstart+displen-1) retval = dispstart+displen-1;
	return retval;
    }
    else {
	if (retval > dispstart-1) retval = dispstart-1;
	if (retval < dispstart-displen) retval = dispstart-displen;
	return retval +1;
    }
}


static void displayResidue(double x)
{
    int qpos, spos;
    static char queryname[NAMESIZE+1];

    if (!*queryname) {
	if (!*qname) strcpy(queryname, "Query");
	else strncpy(queryname, qname, NAMESIZE);
	queryname[NAMESIZE] = 0;
    }

    qpos = x_to_residue(x);

    if (pickedMSP) {
	if (blastx || tblastx) 
	{
	    if (plusmin > 0) spos = (qpos - pickedMSP->qstart)/3 + pickedMSP->sstart;
	    else spos = (pickedMSP->qstart - qpos)/3 + pickedMSP->sstart;
	}
	else if (blastn)
	{
	    if (pickedMSP->qstart < pickedMSP->qend) {
		if (plusmin > 0) spos = qpos - pickedMSP->qstart + pickedMSP->sstart;
		else spos = qpos - pickedMSP->qstart + pickedMSP->sstart;
	    }
	    else {
		if (plusmin > 0) spos = pickedMSP->qstart - qpos + pickedMSP->sstart;
		else spos = pickedMSP->qstart - qpos + pickedMSP->sstart;
	    }
	}
	else spos = qpos - pickedMSP->qstart + pickedMSP->sstart;
	
	if (spos < pickedMSP->sstart) spos = pickedMSP->sstart;
	if (spos > pickedMSP->send)   spos = pickedMSP->send;

	sprintf(message, "%s: %d   ", queryname, qpos + qoffset);
	if (HSPgaps) 
	    strcat(message, "Gapped HSP - no coords");
	else
	    strcat(message, messprintf("%s: %d", pickedMSP->sname, spos));
    }
    else sprintf(message, "%s: %d   No subject picked", queryname, qpos + qoffset);

    graphBoxDraw(messagebox, BLACK, backgColor);
}


static void markDNA(double y)
{
    if (y < 0) return;
    graphXorLine (NAMESIZE+22, separator_y + 1.35 + y, 
		  NAMESIZE+22+displen/3, separator_y + 1.35 +y);
}


static void MiddleDrag (double x, double y) 
{
    if (dragFast) /* Big Picture */
    {
	graphXorBox (BPbox, oldx - BPboxwidth/2, 1.85);
	graphXorBox (BPbox, x - BPboxwidth/2, 1.85);
    }
    else
    {
	graphXorLine (oldx, separator_y, oldx, ny) ;
	graphXorLine (x, separator_y, x, ny);

	if (blastx || tblastx) {
	    markDNA(DNAstep);
	    DNAstep =  abs ( (x_to_residue(x) - dispstart) % 3);
	    markDNA(DNAstep);
	}

	displayResidue(x);
    }

    oldx = x ;
}

static void MiddleUp (double x, double y) 
{ 
    if (dragFast)
	dispstart = BigPictStart + plusmin*((x-BPoffset)/(nx-BPoffset)*BigPictLen - displen/2);
    else
	dispstart = x_to_residue(x) - plusmin*displen/2;

    blviewRedraw();
}

static void MiddleDown (double x, double y) 
{ 
    float nyf;

    graphFitBounds (&nx, 0);
    graphWhere(0,0,0, &nyf);  ny = nyf-.5;

    dragFast = (y < separator_y) ;

    if (dragFast)
	graphXorBox (BPbox, x - BPboxwidth/2, 1.85);
    else {
	graphXorLine (x, separator_y, x, ny) ;

	/* Cleanse messagebox from pick */
	if (messagebox) {
	    *message = 0;
	    graphBoxDraw(messagebox, BLACK, backgColor);
	}
	    
	messagebox = graphBoxStart();
	graphTextPtr (message, NAMESIZE+22, separator_y - .1, 40);
	graphBoxEnd ();

	if (blastx || tblastx) {
	    DNAstep = (x_to_residue(x) - dispstart) % 3;
	    markDNA(DNAstep);
	}

	displayResidue(x);
    }
   
    oldx = x;
    graphRegister (MIDDLE_DRAG, MiddleDrag) ;	/* must redo */
    graphRegister (MIDDLE_UP, MiddleUp) ;
}


static void setHighlight(void)
{
    static char dfault[64];

    if (graphPrompt ("String: (wildcards: * ?)", dfault, "t"))
    {
	/* ANSI way */
	strncpy(searchSeq, freeword(), NAMESIZE+3);
	searchSeq[NAMESIZE+3] = 0;
	for (cp = searchSeq; *cp ; cp++) *cp = ace_upper(*cp);

	/* Non-ANSI bsd way :
	   if (!re_comp(searchSeq)) fprintf(stderr, "%s\n", re_comp(searchSeq));*/

	strncpy(dfault, searchSeq, 63);

	blviewRedraw();
    }
}


static void clrHighlight(void)
{
    MSP *msp;

    /* Clear highlighted */
    *searchSeq = *HighlightSeq = 0;
    pickedMSP = 0;

    /* Unhide hidden matches */
    for (msp = DBhits.next; msp ; msp = msp->next)
	if (msp->score == -999) {
	    msp->score = msp->id;
	    calcID(msp);
	}
	    
    blviewRedraw();
}


static int smartDotterRange(void)
{
    char strand;
    int qstart, qend, sstart, send, extend, len, mid;
    MSP *msp;
     qstart=qend=sstart=send=extend=len=mid=0;
    /* Estimate wanted query region from extent of HSP's */
    strand = ( plusmin > 0 ? '+' : '-');
    for (qstart = 0, msp = DBhits.next; msp; msp = msp->next)	
	if (!strcmp(msp->sname, HighlightSeq) && (msp->frame[1] == strand || blastn)) {
	    
	    if (!qstart) {
		qstart = msp->qstart;
		qend = msp->qend;
		sstart = msp->sstart;
		send = msp->send;
	    }
	    else {
		if (msp->frame[1] == '+') {
		    if (msp->qstart < qstart) {
			qstart = msp->qstart;
			sstart = msp->sstart;
		    }
		    if (msp->qend > qend) {
			qend = msp->qend;
			send = msp->send;
		    }
		}
		else {
		    if (msp->qstart > qstart) {
			qstart = msp->qstart;
			sstart = msp->sstart;
		    }
		    if (msp->qend < qend) {
			qend = msp->qend;
			send = msp->send;
		    }
		}
	    }
	}
    

    if (!qstart) {
	messout("No matches on this strand");
	return 0;
    }
    

    /* Invert qstart and qend if needed.
       This may happen when calling from reversed strand in blastn * /
    if ( (qend < qstart && plusmin == 1) || (qend > qstart && plusmin == -1)) {
	tmp = qstart;
	qstart = qend;
	qend = tmp;
    }*/
    
    /* Extrapolate start to start of vertical sequence */
    extend = sstart;
    if (blastx || tblastx) extend *= 3;
    qstart -= plusmin*extend;
    
    /* Extrapolate end to end of vertical sequence */
    if (!tblastn && (dottersseq = fetchSeqRaw(HighlightSeq))) {
	extend = strlen(dottersseq) - send;
    }
    else {
	extend = 200;
    }
    if (blastx || tblastx) extend *= 3;
    qend += plusmin*extend;

    /* Due to gaps, we might miss the ends - add some more */
    extend = 0.1 * plusmin * (qend-qstart);
    qstart -= extend;
    qend += extend;

    if (blastx || tblastx) {
	/* If sstart and send weren't in the end exons, we'll miss those - add some more */
	extend = 0.2 * plusmin * (qend-qstart);
	qstart -= extend;
	qend += extend;
    }

    /* Keep it within bounds */
    if (qstart < 1) qstart = 1;
    if (qend > qlen) qend = qlen;
    if (qstart > qlen) qstart = qlen;
    if (qend < 1) qend = 1;

    /* Apply min and max limits - min 500 residues, max 10 Mb dots */
    len = plusmin*(qend-qstart);
    mid = qstart + plusmin*len/2;

    if (len < 500) len = 500;
    if (len*(send-sstart) > 1e7) len = 1e7/(send-sstart);

    qstart = mid - plusmin*len/2;
    qend = mid + plusmin*len/2;

    /* Keep it within bounds */
    if (qstart < 1) qstart = 1;
    if (qend > qlen) qend = qlen;
    if (qstart > qlen) qstart = qlen;
    if (qend < 1) qend = 1;
    
    dotterStart = qstart;
    dotterEnd = qend;

    return 1;
}


static void callDotter(void)
{
    static char opts[] = "     ";
    char type=0, *queryseq, *sname;
    int offset;
    MSP *msp;

    if (!*HighlightSeq) {
	messout("Select a sequence first");
	return;
    }

    if (blastp || tblastn) type = 'P';
    else if (blastx) type = 'X';
    else if (blastn || tblastx) type = 'N';

    if (smartDotter) {
	if (!smartDotterRange()) return;
    }
    else 
	dottersseq = fetchSeqRaw(HighlightSeq);

    /* Try to get the subject sequence, in case we're in seqbl mode (only part of seq in MSP) */
    if (!dottersseq || tblastn) {

	/* Check if sequence is passed from acedb */ 
	if (!tblastx) {
	    printf("Looking for sequence stored internally ... ");
	    for (msp = DBhits.next; msp ; msp = msp->next) {
		if (!strcmp(msp->sname, HighlightSeq) && msp->sseq != padseq) {
		    dottersseq = messalloc(strlen(msp->sseq)+1);
		    strcpy(dottersseq, msp->sseq);
		    break;
		}
	    }
	    if (!dottersseq) printf("not ");
	    printf("found\n");
	}
	
	if (!dottersseq) {
	    printf("Can't fetch subject sequence for dotter - aborting\n");
	    messout("Can't fetch subject sequence for dotter - aborting\n");
	    return;
	}
    }

    if (strchr(dottersseq, '-') || tblastn ) messout("Note: the sequence passed to dotter is incomplete");

    if (!*dotterqname) {
	if (!*qname) strcpy(dotterqname, "Blixem-seq");
	else strncpy(dotterqname, qname, NAMESIZE);
	dotterqname[NAMESIZE] = 0;
    }

    /* Get query sequence */
    /* Avoid translating queryseq by pretending to be blastn - very sneaky and dangerous */
    if (blastx || tblastx) blastn = 1;
    if (!(queryseq = getqseq(dotterStart, dotterEnd, q))) return;
    if (blastx || tblastx) blastn = 0;

    if (plusmin > 0) {
	offset = dotterStart-1 + qoffset;
	opts[0] = ' ';
    }
    else {
	offset = dotterEnd-1 + qoffset;
	opts[0] = 'R';
    }

    if ((sname = strchr(HighlightSeq, ':')))
	sname++;
    else
	sname = HighlightSeq;

    opts[1] = (dotterHSPs ? 'H' : ' ');

    opts[2] = (HSPgaps ? 'G' : ' ');

    printf("Calling dotter with query sequence region: %d - %d\n", dotterStart, dotterEnd);

    dotter(type, opts, dotterqname, queryseq, offset, sname, dottersseq, 0, 
	   0, 0, NULL, NULL, NULL, NULL, 0.0, dotterZoom, DBhits.next, qoffset, 0, 0);
}


static void callDotterHSPs(void)
{
    dotterHSPs = 1;
    callDotter();
    dotterHSPs = 0;
}


static int smartDotterRangeSelf(void)
{
    int len, mid;
    
    len = 2000;
    mid = dispstart + plusmin*displen/2;

    dotterStart = mid - plusmin*len/2;
    dotterEnd = mid + plusmin*len/2;

    /* Keep it within bounds */
    if (dotterStart < 1) dotterStart = 1;
    if (dotterStart > qlen) dotterStart = qlen;
    if (dotterEnd > qlen) dotterEnd = qlen;
    if (dotterEnd < 1) dotterEnd = 1;
    
    return 1;
}


static void callDotterSelf(void)
{
    static char opts[] = "     ";
    char type=0, *queryseq;
    int  offset;

    if (blastp || tblastn || tblastx) type = 'P';
    else if (blastx || blastn) type = 'N';

    if (smartDotter)
	if (!smartDotterRangeSelf()) return;

    if (!*dotterqname) {
	if (!*qname) strcpy(dotterqname, "Blixem-seq");
	else strncpy(dotterqname, qname, NAMESIZE);
	dotterqname[NAMESIZE] = 0;
    }

    /* Get query sequence */
    /* Can't do reversed strand since Dotter can't reverse vertical scale */
    /* Avoid translating queryseq by pretending to be blastn - very sneaky and dangerous */
    if (blastx || tblastx) blastn = 1;
    if (!(queryseq = getqseq(min(dotterStart, dotterEnd), max(dotterStart, dotterEnd), q))) return;
    if (blastx || tblastx) blastn = 0;

    dottersseq = messalloc(strlen(queryseq)+1);
    strcpy(dottersseq, queryseq);

    offset = min(dotterStart, dotterEnd)-1 + qoffset;

    printf("Calling dotter with query sequence region: %d - %d\n", dotterStart, dotterEnd);

    dotter(type, opts, dotterqname, queryseq, offset, dotterqname, dottersseq, offset,
	   0, 0, NULL, NULL, NULL, NULL, 0.0, dotterZoom, DBhits.next, qoffset, 0, 0);
}


static void setDotterParams(void)
{
    if (!*dotterqname) {
	if (!*qname) strcpy(dotterqname, "Blixem-seq");
	else strncpy(dotterqname, qname, NAMESIZE);
	dotterqname[NAMESIZE] = 0;
    }

    if (!dotterStart) dotterStart = dispstart-100;
    if (!dotterEnd) dotterEnd = dispstart+100;

    if (!graphPrompt ("Dotter parameters: zoom (compression) factor, "
		     "start, end, Queryname", 
		     messprintf("%d %d %d %s", dotterZoom, 
				dotterStart+qoffset, 
				dotterEnd+qoffset,
				dotterqname),
		     "iiiw"))
	return;

    freeint(&dotterZoom);
    freeint(&dotterStart); dotterStart -= qoffset;
    freeint(&dotterEnd); dotterEnd -= qoffset;
    strncpy(dotterqname, freeword(), NAMESIZE);
    
    smartDotter = 0;
    menuUnsetFlags(menuItem(blixemMenu, autoDotterParamsStr), MENUFLAG_DISABLED);
    graphNewMenu(blixemMenu);

    blviewRedraw();
}


static void autoDotterParams(void)
{
    smartDotter = 1;
    menuSetFlags(menuItem(blixemMenu, autoDotterParamsStr), MENUFLAG_DISABLED);
    graphNewMenu(blixemMenu);
} 



/************************  BLIXEM SETTINGS  WINDOW  **************************/

static void blixemConfColourMenu(KEY key, int box)
{
    /* Taken from ? */
    int *colour;
    if (graphAssFind(assVoid(box+2000), &colour)) { 
	*colour = key;
	graphBoxDraw(box, BLACK, *colour);
	graphRedraw();
	blviewRedraw();
    }
}


static void blixemConfColour(int *colour, int init, float x, float *y, int len, char *text)
{ 
    int box;
    if (text)
	graphText(text, x+len+1, *y);

    box = graphBoxStart();
    graphRectangle(x, *y, x+len, *y+1);
    graphBoxEnd();
    graphBoxFreeMenu(box, blixemConfColourMenu, graphColors);
    *colour = init;
    graphBoxDraw(box, BLACK, init);
    graphAssociate(assVoid(box+2000), colour);

    *y += 1.5;
}


static void graphButtonCheck(char* text, void (*func)(void), float x, float *y, int On)
{
    graphButton(messprintf("%s %s", On ? "*" : " ", text), func, x, *y);

    /* Could be more fancy like this 
       if (On) graphFillRectangle(x-.5, *y, x-2, *y+1.2);
       else graphRectangle(x-.5, *y, x-2, *y+1.2);*/
    
    *y += 1.5;
}


static void graphButtonDisable(char* text, float x, float *y, int On)
{
    int box = graphBoxStart();
    graphText(messprintf("%s %s", On ? "*" : " ", text), x, *y);
    graphBoxEnd();
    graphBoxDraw (box, DARKGRAY, WHITE);

    *y += 1.5;
}


static void setStringentSEGwin(char *cp)
{
    stringentSEGwin = atoi(cp);
    calcSEGarrays(TRUE);
    blviewRedraw();
}
static void setMediumSEGwin(char *cp)
{
    mediumSEGwin = atoi(cp);
    calcSEGarrays(TRUE);
    blviewRedraw();
}
static void setNonglobularSEGwin(char *cp)
{
    nonglobularSEGwin = atoi(cp);
    calcSEGarrays(TRUE);
    blviewRedraw();
}



static void settingsPick(int box)
{
    if (box == stringentSEGbox) graphTextScrollEntry(stringentSEGtx,0,0,0,0,0);
    else if (box == mediumSEGbox) graphTextScrollEntry(mediumSEGtx,0,0,0,0,0);
    else if (box == nonglobularSEGbox) graphTextScrollEntry(nonglobularSEGtx,0,0,0,0,0);
}

static void settingsRedraw(void)
{
    float x1=1, x2=35, y;
    
    static MENUOPT sortMenu[] = { 	
      { sortByScore,      "Score" },
      { sortById,         "Identity" },
      { sortByName,       "Name" },
      { sortByPos,        "Position" },
      { 0, 0 }
    };

    static MENUOPT fetchMenu[] = { 	
#if !defined(WIN32)
	{ fetchByefetch,     "efetch" },
	{ fetchByWWWefetch,  "WWW-efetch" },
#endif
#ifdef ACEDB
	{ fetchByacedb,      "acedb" },
	{ fetchByacedbtext,  "acedb text" },
#endif
	{ 0, 0 }
    };

    if (!graphActivate(settingsGraph)) return;
    graphClear();

    /* Background */
    graphBoxDraw(0, backgColor, backgColor);
    graphFitBounds (&nx, &ny);
    graphColor(backgColor); graphRectangle(0, 0, nx+100, ny+100);graphColor(BLACK);

    y  = 1;
    graphText("Toggles:", x1, y);
    y += 1.5;

    graphButtonCheck(BigPictToggleStr, BigPictToggle, x1, &y, BigPictON);
    if (blastn || blastx || tblastx) {
	if (BigPictON) {
	    graphButtonCheck(BigPictToggleRevStr, BigPictToggleRev, x1, &y, BigPictRev);
	    menuUnsetFlags(menuItem(settingsMenu, BigPictToggleRevStr), MENUFLAG_DISABLED);
	}
	else {
	    graphButtonDisable(BigPictToggleRevStr, x1, &y, BigPictRev);
	    menuSetFlags(menuItem(settingsMenu, BigPictToggleRevStr), MENUFLAG_DISABLED);
	}
    }
    graphButtonCheck(SEGtoggleStr, SEGtoggle, x1, &y, SEGon);

    graphButtonCheck(toggleDESCStr, toggleDESC, x1, &y, DESCon);
    graphButtonCheck(squashMatchesStr, squashMatches, x1, &y, squash);
    graphButtonCheck(toggleIDdotsStr, toggleIDdots, x1, &y, IDdots);
    graphButtonCheck(printColorsStr, printColors, x1, &y, printColorsOn);

/* Old disused toggles ...? */
/*    graphButtonCheck("Display colours", toggleColors, x1, &y, );*/
/*    graphButtonCheck("Printerphilic colours", printColors, x1, &y, printColorsOn);*/


    y = 1;
    graphText("Menus:", x2, y);
    y += 1.5;

    blixemConfColour(&backgColor, backgColor, x2, &y, 3, "Background colour");
    blixemConfColour(&gridColor, gridColor, x2, &y, 3, "Grid colour");
    blixemConfColour(&IDcolor, IDcolor, x2, &y, 3, "Identical residues");
    blixemConfColour(&consColor, consColor, x2, &y, 3, "Conserved residues");


    graphText("Sort HSPs by ", x2, y);
    graphBoxMenu(graphButton(messprintf("%s", sortModeStr), settingsRedraw, x2+13, y), sortMenu);
    y += 1.5;

    graphText("Fetch by ", x2, y);
    graphBoxMenu(graphButton(messprintf("%s", fetchMode), settingsRedraw, x2+9, y), fetchMenu);
    y += 1.5;

    if (SEGon) {
	graphText("Complexity curves:", x2, y);
	y += 1.5;

	sprintf(stringentSEGtx, "%d", stringentSEGwin);
	stringentSEGbox = graphTextScrollEntry (stringentSEGtx, 6, 3, x2+19, y, setStringentSEGwin);
	blixemConfColour(&stringentSEGcolor, stringentSEGcolor, x2, &y, 3, "window length:");

	sprintf(mediumSEGtx, "%d", mediumSEGwin);
	mediumSEGbox = graphTextScrollEntry (mediumSEGtx, 6, 3, x2+19, y, setMediumSEGwin);
	blixemConfColour(&mediumSEGcolor, mediumSEGcolor, x2, &y, 3, "window length:");
	
	sprintf(nonglobularSEGtx, "%d", nonglobularSEGwin);
	nonglobularSEGbox = graphTextScrollEntry (nonglobularSEGtx, 6, 3, x2+19, y, setNonglobularSEGwin);
	blixemConfColour(&nonglobularSEGcolor, nonglobularSEGcolor, x2, &y, 3, "window length:");
    }

    graphRedraw();
}


static void blixemSettings(void)
{
    if (!graphActivate(settingsGraph)) {
	settingsGraph = graphCreate(TEXT_FIT, "Blixem Settings", 0, 0, .6, .3);
	graphRegister(PICK, settingsPick);
	settingsRedraw();
    }
    else graphPop();
}

/*
static void dotterPanel(void)
{
    graphCreate(TEXT_FIT, "Dotter Panel", 0, 0, .3, .1);
    graphButton("Dotter", callDotter, 1, 1);
    graphButton("Dotter HSPs only",callDotterHSPs, 1, 2.5);
    graphButton("Dotter query vs. itself",callDotterSelf, 1, 4);
    graphRedraw();
}
*/

/* Menu checkmarks */

static void menuCheck(MENU menu, int mode, int thismode, char *str)
{
    if (mode == thismode)
	menuSetFlags(menuItem(menu, str), MENUFLAG_TOGGLE_STATE);
    else
	menuUnsetFlags(menuItem(menu, str), MENUFLAG_TOGGLE_STATE);
}


static void setMenuCheckmarks(void)
{
    menuCheck(settingsMenu, sortMode, SORTBYSCORE, SortByScoreStr);
    menuCheck(settingsMenu, sortMode, SORTBYID, SortByIdStr);
    menuCheck(settingsMenu, sortMode, SORTBYNAME, SortByNameStr);
    menuCheck(settingsMenu, sortMode, SORTBYPOS, SortByPosStr);

    menuCheck(settingsMenu, 1, BigPictON, BigPictToggleStr);
    menuCheck(settingsMenu, 1, BigPictRev, BigPictToggleRevStr);
    menuCheck(settingsMenu, 1, IDdots, toggleIDdotsStr);
    menuCheck(settingsMenu, 1, squash, squashMatchesStr);
    menuCheck(settingsMenu, 1, SEGon, SEGtoggleStr);
    menuCheck(settingsMenu, 1, printColorsOn, printColorsStr);
    menuCheck(settingsMenu, 1, colortoggle, toggleColorsStr);
    menuCheck(settingsMenu, 1, verbose, toggleVerboseStr);
    menuCheck(settingsMenu, 1, HiliteSins, toggleHiliteSinsStr);
    menuCheck(settingsMenu, 1, HiliteUpperOn, toggleHiliteUpperStr);
    menuCheck(settingsMenu, 1, HiliteLowerOn, toggleHiliteLowerStr);
    menuCheck(settingsMenu, 1, DESCon, toggleDESCStr);
    graphNewBoxMenu(settingsButton, settingsMenu);
    graphNewMenu(blixemMenu);
}
 
