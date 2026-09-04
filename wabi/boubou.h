
#ifndef  BOU_DEF_00

#define  BOU_DEF_00
#define  IBC_0 1
#define  WORM_0 1
#define  WARREN_0 1
#define  SEQC_0 1

typedef struct pShowRunStruct { char *nam ; unsigned int flag ;} PSR ;

#define PSR_runs  0x1
#define PSR_groups  0x2
#define PSR_showOpenButton  0x4
#define PSR_showClosedButton  0x8


typedef struct pairGeneGroupStruct { char *nam ; unsigned int flag ;} PGG ;

/* all projects */
#define PGG_NOZOOM     0x1
#define PGG_SOLEXA     0x2
#define PGG_PA         0x4
#define PGG_NUC        0x8
#define PGG_TOTAL     0x10
#define PGG_CAP       0x20
#define PGG_rna       0x40
#define PGG_u         0x80

#define PGG_f      0x1000
#define PGG_r      0x2000
#define PGG_ns     0x4000
#define PGG_cover 0x10000
#define PGG_nu    0x20000
#define PGG_pp    0x40000
#define PGG_ends  0x80000


#define PGG_uf    0x11000
#define PGG_ur    0x12000
#define PGG_uns   0x14000
#define PGG_nuf   0x21000
#define PGG_nur   0x22000
#define PGG_ppf   0x41000
#define PGG_ppr   0x42000

#define PGG_ELF  0x181000
#define PGG_ERF  0x281000
#define PGG_ELR  0x182000
#define PGG_ERR  0x282000

#define PGG_endRatios   0x400000
#define PGG_endRatioLF  0x501000
#define PGG_endRatioRF  0x601000
#define PGG_endRatioLR  0x502000
#define PGG_endRatioRR  0x602000

#define PGG_isGroup   0x800000

/* Eric */
#define PGG_S          0x80
#define PGG_G1        0x100
#define PGG_ES        0x200
#define PGG_MS        0x400
#define PGG_RED       0x800
#define PGG_BM       0x1000
#define PGG_AFFY     0x2000
#define PGG_TF       0x4000
#define PGG_DYNA     0x8000

/* Warren */
#define PGG_STAT3    0x100
#define PGG_STAT5A   0x200
#define PGG_STAT5B   0x400
#define PGG_IRF4     0x800
#define PGG_IgG     0x1000
#define PGG_CD4     0x2000
#define PGG_BCELL   0x4000
#define PGG_UTR     0x8000

/* G26 */
#define PGG_AAA    0x100
#define PGG_KF     0x200
#define PGG_G26    0x400
#define PGG_small  0x800

/* SEQC */
#define PGG_ILM    0x100
#define PGG_LIF    0x200
#define PGG_HEL   0x400
#define PGG_R454  0x800
#define PGG_Brain  0x1000
#define PGG_UHR  0x2000
#define PGG_any  0x4000
/* IBC */
#define PGG_IBC  0x8000

#ifdef BOU_PNS
static PGG pggGroups[] =
  {
#ifdef BOU_DEF
    { "pA", PGG_PA},
    { "pA-", PGG_NUC},
    { "Total", PGG_TOTAL},
    { "S", PGG_S},
    { "G1", PGG_G1},
    { "ES", PGG_ES},
    { "MS", PGG_MS},
    { "RED", PGG_RED},
    { "BM", PGG_BM},
    { "DYNAMIC", PGG_DYNA}, 
    { "AFFY", PGG_AFFY},
    { "Solexa", PGG_SOLEXA},
    { "Transfrag", PGG_TF} ,
#endif

#ifdef WARREN
    { "STAT3", PGG_STAT3},
    { "STAT5A", PGG_STAT5A},
    { "STAT5B", PGG_STAT5B},
    { "IRF4", PGG_IRF4},
    { "IgG", PGG_IgG},
    { "CD4", PGG_CD4},
    { "Bcell", PGG_BCELL},
#endif

#ifdef WORM
    { "RNA-seq", PGG_rna},
    { "AAA", PGG_AAA},
    { "26G", PGG_G26},
    { "small", PGG_small},
    { "+", PGG_f},
    { "-", PGG_r},
    { "ns",PGG_ns},

#endif

#ifdef SEQC
    { "ILM", PGG_ILM},
    { "LIF", PGG_LIF},
    { "HEL", PGG_HEL},
    { "R454", PGG_R454},
    { "Brain", PGG_Brain},
    { "UHR", PGG_UHR},
    { "Any", PGG_rna},
    { "ns", PGG_ns},
    { "pA", PGG_PA},
    { "Total", PGG_TOTAL},
    { "Cap", PGG_CAP},
#endif

#ifdef IBC
    { "IBC", PGG_IBC},
    { "BodyMap", PGG_any},
    { "+", PGG_f},
    { "-", PGG_r},
    { "ns", PGG_ns},
#endif
    { "+", PGG_f},
    { "-", PGG_r},
    { "Ends", PGG_ends},
    { "Non-unique", PGG_nu},
    { "Partial", PGG_pp},
    { "Cover", PGG_cover},
    { "End-ratio", PGG_endRatios},

    {0, 0}
  } ;
#endif

#ifdef BOU_DEF
#define SOLEXA_STEP 1000
#else
#define SOLEXA_STEP 10
#endif

typedef struct pairNamStruct { const char *p; unsigned int flag ; float damper; char *nam, *nam2 ; float av, s99 ; int col ; BOOL open ; int isRna ; float maxRatio, maxIndex ; Array signal ; } PNS ;
typedef struct solexaNamStruct { const char *p; unsigned int flag ; float damper; char *nam ; float av, s99 ; int col ; BOOL open ; int noSmoothing ; float maxRatio, maxIndex ; Array signal ; } PNX ;


/* as of 2011, this table is contructed dynamically be exploring the class Run -> Wiggle */

/* up to jan 2011, these tables were hard coded */ 
#ifdef OLDJUNK
static PNX *solexaAll ;

static PNX solexaAllOld[] = 
  {
#ifdef BOU_DEF
    { "SX_S", PGG_SOLEXA, 0, "SX_SES", 0, 10, LIGHTGREEN, FALSE , FALSE, 0, 3, 0},
    { "SX_G1", PGG_SOLEXA , 0, "SX_G1ES",  0, 10, LIGHTBLUE, FALSE , FALSE, 0, 5, 0},
    { "SOLID_S", PGG_SOLEXA, 0, "SLD_SES", .0, 10, GREEN, FALSE , FALSE, 0, 40, 0},
    { "SOLID_G1", PGG_SOLEXA , 0, "SLD_G1ES",  .0, 10, BLUE, FALSE , FALSE, 0, 65, 0},
    { "S_G1_smoothed", PGG_SOLEXA, 0, "S_G1_smoothed", 200, 10, RED, FALSE , TRUE, 0, 2000, 0},
    { "S_G1_eric", PGG_SOLEXA , 0, "S_G1_eric",  200, 10, ORANGE, FALSE , TRUE, 0, 2000, 0},
#endif
    /* maxIndex is the total number of 'uniquely' aligned probes at 0/1/2 errors 
     * this number should eb used to rationalize the intensities, say to 3.000.000 
     */
#ifdef WARREN
{"STAT3pIL21IRF4KO", PGG_STAT3 , 0, "STAT3pIL21IRF4KO", 0, 0, VIOLET, TRUE, 0, 0, 4050573, 0},
{"STAT3mIL21IRF4KO", PGG_STAT3 , 0, "STAT3mIL21IRF4KO", 0, 0, PALEVIOLET, TRUE, 0, 0, 3718865, 0},
  /*
    {"STAT3mIL21ter", PGG_STAT3 , 0, "STAT3mIL21ter", 0, 0, RED2, TRUE, 0, 0, 3401892, 0},
    {"STAT3pIL21ter", PGG_STAT3 , 0, "STAT3pIL21ter", 0, 0, RED5, TRUE, 0, 0, 3423556, 0},
    {"Ig3pIL21ter", PGG_IgG , 0, "Ig3pIL21ter", 0, 0, LIGHTGRAY, TRUE, 0, 0, 2886982, 0},

    {"STAT5AmIL2WTa1", PGG_STAT5A , 0, "STAT5AmIL2WTa1", 0, 0, RED2, TRUE, 0, 0, 2716619, 0},
    {"STAT5AmIL2WTa2", PGG_STAT5A , 0, "STAT5AmIL2WTa2", 0, 0, RED3, TRUE, 0, 0, 2531415, 0},
    {"STAT5AmIL2WTb1", PGG_STAT5A , 0, "STAT5AmIL2WTb1", 0, 0, RED4, TRUE, 0, 0, 2656330, 0},
    {"STAT5AmIL2WTb2", PGG_STAT5A , 0, "STAT5AmIL2WTb2", 0, 0, RED5, TRUE, 0, 0, 4159330, 0},
    {"STAT5ApIL2WTa1", PGG_STAT5A , 0, "STAT5ApIL2WTa1", 0, 0, RED6, TRUE, 0, 0, 3226614, 0},
    {"STAT5ApIL2WTa2", PGG_STAT5A , 0, "STAT5ApIL2WTa2", 0, 0, RED7, TRUE, 0, 0, 2362832, 0},
    {"STAT5ApIL2WTb1", PGG_STAT5A , 0, "STAT5ApIL2WTb1", 0, 0, RED8, TRUE, 0, 0, 2480506, 0},
    {"STAT5ApIL2WTb2", PGG_STAT5A , 0, "STAT5ApIL2WTb2", 0, 0, DARKRED, TRUE, 0, 0, 4964830, 0},
    {"STAT5ApIL2WTc", PGG_STAT5A , 0, "STAT5ApIL2WTc", 0, 0, PURPLE, TRUE, 0, 0, 3297727, 0},
    {"STAT5AmIL2KIa1", PGG_STAT5A , 0, "STAT5AmIL2KIa1", 0, 0, BLUE2, TRUE, 0, 0, 3607158, 0},
    {"STAT5AmIL2KIa2", PGG_STAT5A , 0, "STAT5AmIL2KIa2", 0, 0, BLUE3, TRUE, 0, 0, 3527218, 0},
    {"STAT5AmIL2KIb1", PGG_STAT5A , 0, "STAT5AmIL2KIb1", 0, 0, BLUE4, TRUE, 0, 0, 3071480, 0},
    {"STAT5AmIL2KIb2", PGG_STAT5A , 0, "STAT5AmIL2KIb2", 0, 0, BLUE5, TRUE, 0, 0, 4096243, 0},
    {"STAT5ApIL2KIa1", PGG_STAT5A , 0, "STAT5ApIL2KIa1", 0, 0, BLUE6, TRUE, 0, 0, 2505149, 0},
    {"STAT5ApIL2KIa2", PGG_STAT5A , 0, "STAT5ApIL2KIa2", 0, 0, BLUE7, TRUE, 0, 0, 2432144, 0},
    {"STAT5ApIL2KIb1", PGG_STAT5A , 0, "STAT5ApIL2KIb1", 0, 0, BLUE8, TRUE, 0, 0, 2667682, 0},
    {"STAT5ApIL2KIb2", PGG_STAT5A , 0, "STAT5ApIL2KIb2", 0, 0, DARKBLUE, TRUE, 0, 0, 4989311, 0},
    {"STAT5BmIL2WTa1", PGG_STAT5B , 0, "STAT5BmIL2WTa1", 0, 0, GREEN2, TRUE, 0, 0, 2423314, 0},
    {"STAT5BmIL2WTa2", PGG_STAT5B , 0, "STAT5BmIL2WTa2", 0, 0, GREEN3, TRUE, 0, 0, 2819117, 0},
    {"STAT5BmIL2WTb1", PGG_STAT5B , 0, "STAT5BmIL2WTb1", 0, 0, GREEN4, TRUE, 0, 0, 2794641, 0},
    {"STAT5BmIL2WTb2", PGG_STAT5B , 0, "STAT5BmIL2WTb2", 0, 0, GREEN5, TRUE, 0, 0, 4565353, 0},
    {"STAT5BpIL2WTa1", PGG_STAT5B , 0, "STAT5BpIL2WTa1", 0, 0, GREEN6, TRUE, 0, 0, 2174197, 0},
    {"STAT5BpIL2WTa2", PGG_STAT5B , 0, "STAT5BpIL2WTa2", 0, 0, GREEN7, TRUE, 0, 0, 2885903, 0},
    {"STAT5BpIL2WTb1", PGG_STAT5B , 0, "STAT5BpIL2WTb1", 0, 0, GREEN8, TRUE, 0, 0, 2462385, 0},
    {"STAT5BpIL2WTb2", PGG_STAT5B , 0, "STAT5BpIL2WTb2", 0, 0, GREEN, TRUE, 0, 0, 4577985, 0},
    {"STAT5BmIL2KIa1", PGG_STAT5B , 0, "STAT5BmIL2KIa1", 0, 0, PALEVIOLET, TRUE, 0, 0, 2464933, 0},
    {"STAT5BmIL2KIa2", PGG_STAT5B , 0, "STAT5BmIL2KIa2", 0, 0, LIGHTVIOLET, TRUE, 0, 0, 3158397, 0},
    {"STAT5BmIL2KIb1", PGG_STAT5B , 0, "STAT5BmIL2KIb1", 0, 0, VIOLET, TRUE, 0, 0, 2984578, 0},
    {"STAT5BmIL2KIb2", PGG_STAT5B , 0, "STAT5BmIL2KIb2", 0, 0, DARKVIOLET, TRUE, 0, 0, 4417176, 0},
    {"STAT5BpIL2KIa1", PGG_STAT5B , 0, "STAT5BpIL2KIa1", 0, 0, PALEMAGENTA, TRUE, 0, 0, 2491909, 0},
    {"STAT5BpIL2KIa2", PGG_STAT5B , 0, "STAT5BpIL2KIa2", 0, 0, LIGHTMAGENTA, TRUE, 0, 0, 3154645, 0},
    {"STAT5BpIL2KIb1", PGG_STAT5B , 0, "STAT5BpIL2KIb1", 0, 0, MAGENTA, TRUE, 0, 0, 2641837, 0},
    {"STAT5BpIL2KIb2", PGG_STAT5B , 0, "STAT5BpIL2KIb2", 0, 0, PURPLE, TRUE, 0, 0, 5088322, 0},
   */

{"IRF4all_12", PGG_IRF4  , 0, "IRF4all_12",       0, 0, GREEN5, TRUE, 0, 0, 50000000, 0},
{"IRF4inIL21_8", PGG_IRF4  , 0, "IRF4inIL21_8" ,  0, 0, GREEN7, TRUE, 0, 0, 21521363, 0},


{"IRF4pIL21IRF4KO", PGG_IRF4  , 0, "IRF4pIL21IRF4KO",  0, 0, GREEN6, TRUE, 0, 0, 5088322, 0},
{"IRF4mIL21IRF4KO", PGG_IRF4  , 0, "IRF4mIL21IRF4KO",  0, 0, GREEN3, TRUE, 0, 0, 5088322, 0},

{"IRF4pIL21WT_4", PGG_IRF4  , 0, "IRF4pIL21WT_4", 0, 0, CYAN, TRUE, 0, 0, 10237360, 0},
{"IRF4mIL21WT_4", PGG_IRF4  , 0, "IRF4mIL21WT_4", 0, 0, PALECYAN, TRUE, 0, 0, 11284003, 0},


{"STAT3all_8", PGG_STAT3  , 0, "STAT3all_8", 0, 0, MAGENTA, TRUE, 0, 0, 24000000, 0},
{".STAT3all_8", PGG_STAT3  , 0, "STAT3all_8",0, 0, MAGENTA, TRUE, 0, 0, 24000000, 0},

{"STAT3pIL21IRF4KO_3", PGG_IRF4 | PGG_STAT3 , 0, "STAT3pIL21IRF4KO_3", 0, 0, VIOLET, TRUE    , 0, 0, 5000000, 0},
{"STAT3mIL21IRF4KO_2", PGG_IRF4 | PGG_STAT3 , 0, "STAT3mIL21IRF4KO_2", 0, 0, PALEVIOLET, TRUE, 0, 0, 5000000, 0},

{"STAT3pIL21WT_7", PGG_IRF4 | PGG_STAT3 , 0, "STAT3pIL21WT_7", 0, 0, RED6, TRUE, 0, 0, 18000000, 0},
{"STAT3mIL21WT_6", PGG_IRF4 | PGG_STAT3 , 0, "STAT3mIL21WT_6", 0, 0, RED2, TRUE, 0, 0, 17000000, 0},

{"Ig3all_5", PGG_IgG , 0, "Ig3all_5", 0, 0, GRAY, TRUE, 0, 0, 5000000, 0},
{"Ig4all_2", PGG_IgG , 0, "Ig4all_2", 0, 0, DARKGRAY, TRUE, 0, 0, 3633052, 0},
{"Ig3pIL21WT_3", PGG_IgG , 0, "Ig3pIL21WT_3", 0, 0, LIGHTGRAY, TRUE, 0, 0, 7000000, 0},

  /*
    {"Ig5pIL2WT", PGG_IgG , 0, "Ig5pIL2WT",   0, 0, BLACK, TRUE, 0, 0, 2708921, 0},
    {"Ig5pIL2WTb", PGG_IgG , 0, "Ig5pIL2WTb", 0, 0, BLACK, TRUE, 0, 0, 2593437, 0},
    {"Ig5mIL2WT", PGG_IgG , 0, "Ig5mIL2WT",   0, 0, DARKGRAY, TRUE, 0, 0, 2704182, 0},
    {"Ig5mIL2WTb", PGG_IgG , 0, "Ig5mIL2WTb", 0, 0, DARKGRAY, TRUE, 0, 0, 2845987, 0},
    {"Ig5pIL2KI", PGG_IgG , 0, "Ig5pIL2KI",   0, 0, GRAY, TRUE, 0, 0, 2246479, 0},
    {"Ig5pIL2KIb", PGG_IgG , 0, "Ig5pIL2KIb", 0, 0, LIGHTGRAY, TRUE, 0, 0, 2896308, 0},
    {"Ig5mIL2KI", PGG_IgG , 0, "Ig5mIL2KI",   0, 0, LIGHTGRAY, TRUE, 0, 0, 2950304, 0},
    {"Ig5mIL2KIb", PGG_IgG , 0, "Ig5mIL2KIb", 0, 0, GRAY, TRUE, 0, 0, 3202045, 0},
  */
#endif

#ifdef WORM
{"Embryo.f", PGG_rna | PGG_f , 0, "Emb+",  0, 0, GREEN3, TRUE, 0, 0, 2716619, 0},
{"L1.f",     PGG_rna | PGG_f , 0, "L1+",   0, 0, GREEN5, TRUE, 0, 0, 2716619, 0},
{"L2.f",     PGG_rna | PGG_f , 0, "L2+",   0, 0, GREEN7, TRUE, 0, 0, 2716619, 0},
{"L3.f",     PGG_rna | PGG_f , 0, "L3+",   0, 0, BLUE3, TRUE, 0, 0, 2716619, 0},
{"Dauer.f",  PGG_rna | PGG_f , 0, "Dauer+",0, 0, VIOLET, TRUE, 0, 0, 2716619, 0},
{"L4.f",     PGG_rna | PGG_f , 0, "L4+",   0, 0, BLUE5, TRUE, 0, 0, 2716619, 0},
{"Adult.f",  PGG_rna | PGG_f , 0, "Ad+",   0, 0, BLUE7, TRUE, 0, 0, 2716619, 0},
{"Mixed.f",  PGG_rna | PGG_f , 0, "Mix+",  0, 0, YELLOW, TRUE, 0, 0, 2716619, 0},
{"Male.f",   PGG_rna | PGG_f , 0, "Male+",  0, 0, CYAN, TRUE, 0, 0, 2716619, 0},

{"Embryo.r", PGG_rna | PGG_r , 0, "Emb-", 0, 0, GREEN3, TRUE, 0, 0, 2716619, 0},
{"L1.r",     PGG_rna | PGG_r , 0, "L1-",  0, 0, GREEN5, TRUE, 0, 0, 2716619, 0},
{"L2.r",     PGG_rna | PGG_r , 0, "L2-",  0, 0, GREEN7, TRUE, 0, 0, 2716619, 0},
{"L3.r",     PGG_rna | PGG_r , 0, "L3-",  0, 0, BLUE3, TRUE, 0, 0, 2716619, 0},
{"Dauer.r",  PGG_rna | PGG_r , 0, "Dauer-", 0, 0, VIOLET, TRUE, 0, 0, 2716619, 0},
{"L4.r",     PGG_rna | PGG_r , 0, "L4-",  0, 0, BLUE5, TRUE, 0, 0, 2716619, 0},
{"Adult.r",  PGG_rna | PGG_r , 0, "Ad-",  0, 0, BLUE7, TRUE, 0, 0, 2716619, 0},
{"Mixed.r",  PGG_rna | PGG_r , 0, "Mix-", 0, 0, RED, TRUE, 0, 0, 2716619, 0},
{"Male.r",   PGG_rna | PGG_r , 0, "Male-",  0, 0, CYAN, TRUE, 0, 0, 2716619, 0},

{"Embryo.ns", PGG_rna | PGG_ns , 0, "Emb", 0, 0, GREEN3, TRUE, 0, 0, 2716619, 0},
{"L1.ns",     PGG_rna | PGG_ns , 0, "L1",  0, 0, GREEN5, TRUE, 0, 0, 2716619, 0},
{"L2.ns",     PGG_rna | PGG_ns , 0, "L2",  0, 0, GREEN7, TRUE, 0, 0, 2716619, 0},
{"L3.ns",     PGG_rna | PGG_ns , 0, "L3",  0, 0, BLUE3, TRUE, 0, 0, 2716619, 0},
{"Dauer.ns",  PGG_rna | PGG_ns , 0, "Dauer", 0, 0, VIOLET, TRUE, 0, 0, 2716619, 0},
{"L4.ns",     PGG_rna | PGG_ns , 0, "L4",  0, 0, BLUE5, TRUE, 0, 0, 2716619, 0},
{"Adult.ns",  PGG_rna | PGG_ns , 0, "Ad",  0, 0, BLUE7, TRUE, 0, 0, 2716619, 0},
{"Mixed.ns",  PGG_rna | PGG_ns , 0, "Mix", 0, 0, GRAY, TRUE, 0, 0, 2716619, 0},
{"Male.ns",   PGG_rna | PGG_r , 0, "Male",  0, 0, CYAN, TRUE, 0, 0, 2716619, 0},

{"AAA.f", PGG_AAA | PGG_f , 0, "pA.f" , 0, 0, BLACK, FALSE, TRUE, 0, 2716619, 0},
{"AAA.r", PGG_AAA | PGG_r , 0, "pA.r",  0, 0, BLACK, FALSE, TRUE, 0, 2716619, 0},
{"Small.f", PGG_small | PGG_f, 0, "Small+", 0, 0, PALEMAGENTA, TRUE, 0, 0, 37000, 0},
{"Small.r", PGG_small | PGG_r, 0, "Small-", 0, 0, PALEORANGE, TRUE, 0, 0, 37000, 0},
{"G26.f",   PGG_small | PGG_G26 | PGG_f , 0, "26G+",  0, 0, MAGENTA, TRUE, 0, 0, 100000, 0},
{"G26.r",   PGG_small | PGG_G26 | PGG_r , 0, "26G-",  0, 0, ORANGE, TRUE, 0, 0, 100000, 0},
{"ATC26.f", PGG_small | PGG_G26 | PGG_f , 0, "ATC26+", 0, 0, PALEMAGENTA, TRUE, 0, 0, 37000, 0},
{"ATC26.r", PGG_small | PGG_G26 | PGG_r , 0, "ATC26-", 0, 0, PALEORANGE, TRUE, 0, 0, 25000, 0},
{"non26.f", PGG_small | PGG_G26 | PGG_f , 0, "small+", 0, 0, PALEMAGENTA, TRUE, 0, 0, 2000000, 0},
{"non26.r", PGG_small | PGG_G26  | PGG_r, 0, "small-", 0, 0, PALEORANGE, TRUE, 0, 0, 2000000, 0},
#endif

#ifdef SEQC_ZERO

{"Brain.f", PGG_rna | PGG_Brain | PGG_f , 0, "Brain+",  0, 0, GREEN4, TRUE, 0, 0, 2716619, 0},
{"UHR.f",   PGG_rna | PGG_UHR   | PGG_f , 0, "UHR+",    0, 0, GREEN2, TRUE, 0, 0, 2716619, 0},
{"any.f", PGG_rna | PGG_any | PGG_f , 0, "Body+",  0, 0, GREEN6, TRUE, 0, 0, 2716619, 0},

{"Brain.r", PGG_rna | PGG_Brain | PGG_r , 0, "Brain-",  0, 0, RED4, TRUE, 0, 2716619, 0},
{"UHR.r",   PGG_rna | PGG_UHR   | PGG_r , 0, "UHR-",    0, 0, RED2, TRUE, 0, 0, 2716619, 0},
{"any.r", PGG_rna | PGG_any | PGG_r , 0, "Body-",  0, 0, RED6, TRUE, 0, 2716619, 0},

{"Brain.ns", PGG_rna | PGG_Brain | PGG_ns , 0, "Brain_ns", 0, 0, BLUE4, TRUE, 0, 0, 2716619, 0},
{"UHR.ns",   PGG_rna | PGG_UHR   | PGG_ns , 0, "UHR_ns",   0, 0, BLUE2, TRUE, 0, 0, 2716619, 0},
{"any.ns", PGG_rna | PGG_any | PGG_ns , 0, "Body_ns", 0, 0, BLUE6, TRUE, 0, 0, 2716619, 0},

#endif

#ifdef IBC

{"s928.ns", PGG_IBC | PGG_ns , 0, "s928", 0, 0, GREEN4, TRUE, 0, 0, 2716619, 0},
{"s929.ns", PGG_IBC | PGG_ns , 0, "s929",   0, 0, BLUE4, TRUE, 0, 0, 2716619, 0},
{"s930.ns", PGG_IBC | PGG_ns , 0, "s930",   0, 0, BLUE5, TRUE, 0, 0, 2716619, 0},
{"s931.ns", PGG_IBC | PGG_ns , 0, "s931",   0, 0, BLUE6, TRUE, 0, 0, 2716619, 0},
{"s932.ns", PGG_IBC | PGG_ns , 0, "s932",   0, 0, BLUE7, TRUE, 0, 0, 2716619, 0},
{"s933.ns", PGG_IBC | PGG_ns , 0, "s933",   0, 0, BLUE8, TRUE, 0, 0, 2716619, 0},
{"s934.ns", PGG_IBC | PGG_ns , 0, "s934",   0, 0, RED4, TRUE, 0, 0, 2716619, 0},
{"s935.ns", PGG_IBC | PGG_ns , 0, "s935",   0, 0, RED5, TRUE, 0, 0, 2716619, 0},
{"s936.ns", PGG_IBC | PGG_ns , 0, "s936",   0, 0, RED6, TRUE, 0, 0, 2716619, 0},
{"s937.ns", PGG_IBC | PGG_ns , 0, "s937",   0, 0, RED7, TRUE, 0, 0, 2716619, 0},
{"s938.ns", PGG_IBC | PGG_ns , 0, "s938",   0, 0, RED8, TRUE, 0, 0, 2716619, 0},
{"any.f",   PGG_rna | PGG_any | PGG_f , 0, "BodyMap+",  0, 0, YELLOW, TRUE, 0, 0, 2716619, 0},
{"any.r",   PGG_rna | PGG_any | PGG_r , 0, "BodyMap-",  0, 0, LIGHTORANGE, TRUE, 0, 2716619, 0},
{"any.ns",  PGG_rna | PGG_any | PGG_ns , 0, "BodyMap_ns", 0, 0, GREEN4, TRUE, 0, 0, 2716619, 0},


#endif

#ifdef SEQC_ZERO
{"ILM_S.Brain.nu.f", PGG_nu | PGG_f | PGG_ILM | PGG_Brain  | PGG_rna , 0, "ILM_S.Brain.nu.f", 0, 0, GREEN5, TRUE, 0, 0, 3000000, 0},
{"LIF_S.Brain.nu.f", PGG_nu | PGG_f | PGG_LIF | PGG_Brain  | PGG_rna , 0, "LIF_S.Brain.nu.f", 0, 0, GREEN6, TRUE, 0, 0, 3000000, 0},
{"HELdge.Brain.nu.f", PGG_nu | PGG_f | PGG_HEL | PGG_Brain | PGG_rna , 0, "HELdge.Brain.nu.f", 0, 0, GREEN5, TRUE, 0, 0, 3000000, 0},

{"ILM_S.UHR.nu.f", PGG_nu | PGG_f | PGG_ILM | PGG_UHR | PGG_rna  , 0, "ILM_S.UHR.nu.f", 0, 0, RED4, TRUE, 0, 0, 3000000, 0},
{"LIF_S.UHR.nu.f", PGG_nu | PGG_f | PGG_LIF | PGG_UHR | PGG_rna  , 0, "LIF_S.UHR.nu.f", 0, 0, RED6, TRUE, 0, 0, 3000000, 0},
{"HELdge.UHR.nu.f", PGG_nu | PGG_f | PGG_HEL | PGG_UHR  | PGG_rna , 0, "HELdge.UHR.nu.f", 0, 0, RED5, TRUE, 0, 0, 3000000, 0},


{"ILM_S.Brain.nu.r", PGG_nu | PGG_r | PGG_ILM | PGG_Brain  | PGG_rna , 0, "ILM_S.Brain.nu.r", 0, 0, BLUE5, TRUE, 0, 0, 3000000, 0},
{"LIF_S.Brain.nu.r", PGG_nu | PGG_r | PGG_LIF | PGG_Brain  | PGG_rna , 0, "LIF_S.Brain.nu.r", 0, 0, BLUE7, TRUE, 0, 0, 3000000, 0},
{"HELdge.Brain.nu.r", PGG_nu | PGG_r | PGG_HEL | PGG_Brain | PGG_rna , 0, "HELdge.Brain.nu.r", 0, 0, BLUE6, TRUE, 0, 0, 3000000, 0},

{"ILM_S.UHR.nu.r", PGG_nu | PGG_r | PGG_ILM | PGG_UHR | PGG_rna  , 0, "ILM_S.UHR.nu.r", 0, 0, LIGHTVIOLET, TRUE, 0, 0, 3000000, 0},
{"LIF_S.UHR.nu.r", PGG_nu | PGG_r | PGG_LIF | PGG_UHR | PGG_rna  , 0, "LIF_S.UHR.nu.r", 0, 0, VIOLET, TRUE, 0, 0, 3000000, 0},
{"HELdge.UHR.nu.r", PGG_nu | PGG_r | PGG_HEL | PGG_UHR  | PGG_rna , 0, "HELdge.UHR.nu.r", 0, 0, LIGHTVIOLET, TRUE, 0, 0, 3000000, 0},





{"ILM_nS.Brain.ns", PGG_nu | PGG_ns | PGG_ILM | PGG_Brain  | PGG_rna , 0, "ILM_nS.Brain.ns", 0, 0, GREEN1, TRUE, 0, 0, 3000000, 0},
{"HEL.Brain.nu.f", PGG_nu | PGG_f | PGG_HEL | PGG_Brain  | PGG_rna , 0, "HEL.Brain.nu.f", 0, 0, GREEN2, TRUE, 0, 0, 3000000, 0},
{"HEL.Brain.nu.r", PGG_nu | PGG_r | PGG_HEL | PGG_Brain  | PGG_rna , 0, "HEL.Brain.nu.r", 0, 0, GREEN2, TRUE, 0, 0, 3000000, 0},
{"R454.Brain.ns", PGG_nu | PGG_ns | PGG_R454 | PGG_Brain  | PGG_rna , 0, "R454.Brain.ns", 0, 0, GREEN1, TRUE, 0, 0, 3000000, 0},


{"ILM_nS.UHR.ns", PGG_nu | PGG_ns | PGG_ILM | PGG_UHR | PGG_rna  , 0, "ILM_nS.UHR.ns", 0, 0, RED1, TRUE, 0, 0, 3000000, 0},
{"HEL.UHR.nu.f", PGG_nu | PGG_f | PGG_HEL | PGG_UHR | PGG_rna  , 0, "HEL.UHR.nu.f", 0, 0, RED2, TRUE, 0, 0, 3000000, 0},
{"HEL.UHR.nu.r", PGG_nu | PGG_r | PGG_HEL | PGG_UHR | PGG_rna  , 0, "HEL.UHR.nu.r", 0, 0, RED2, TRUE, 0, 0, 3000000, 0},
{"R454.UHR.ns", PGG_nu | PGG_ns | PGG_R454 | PGG_UHR | PGG_rna  , 0, "R454.UHR.ns", 0, 0, RED1, TRUE, 0, 0, 3000000, 0},



{"ILM_RS.Brain.nu.f", PGG_nu | PGG_f | PGG_ILM | PGG_Brain | PGG_TOTAL , 0, "ILM_RS.Brain.nu.f", 0, 0, LIGHTGREEN, TRUE, 0, 0, 3000000, 0},
{"LIF_R.Brain.nu.f", PGG_nu | PGG_f | PGG_LIF | PGG_Brain | PGG_TOTAL , 0, "LIF_R.Brain.nu.f", 0, 0, DARKGREEN, TRUE, 0, 0, 3000000, 0},
{"ILM_RS.UHR.nu.f", PGG_nu | PGG_f | PGG_ILM | PGG_UHR  | PGG_TOTAL , 0, "ILM_RS.UHR.nu.f", 0, 0, DARKRED, TRUE, 0, 0, 3000000, 0},


{"ILM_RS.Brain.nu.r", PGG_nu | PGG_r | PGG_ILM | PGG_Brain | PGG_TOTAL , 0, "ILM_RS.Brain.nu.r", 0, 0, CYAN, TRUE, 0, 0, 3000000, 0},
{"LIF_R.Brain.nu.r", PGG_nu | PGG_r | PGG_LIF | PGG_Brain | PGG_TOTAL , 0, "LIF_R.Brain.nu.r", 0, 0, DARKCYAN, TRUE, 0, 0, 3000000, 0},
{"ILM_RS.UHR.nu.r", PGG_nu | PGG_r | PGG_ILM | PGG_UHR  | PGG_TOTAL , 0, "ILM_RS.UHR.nu.r", 0, 0, DARKVIOLET, TRUE, 0, 0, 3000000, 0},



{"ILM_RnS.Brain.ns", PGG_nu | PGG_ns | PGG_ILM | PGG_Brain | PGG_TOTAL , 0, "ILM_RnS.Brain.ns", 0, 0, LAVANDER, TRUE, 0, 0, 3000000, 0},
{"ILM_RnS.UHR.ns", PGG_nu | PGG_ns | PGG_ILM | PGG_UHR  | PGG_TOTAL , 0, "ILM_RnS.UHR.ns", 0, 0, BROWN, TRUE, 0, 0, 3000000, 0},





  /* absent */
{"LIF_R.UHR.nu.r", PGG_nu | PGG_r | PGG_LIF | PGG_UHR  | PGG_TOTAL , 0, "LIF_R.UHR.nu.r", 0, 0, DARKVIOLET, TRUE, 0, 0, 3000000, 0},
{"LIF_R.UHR.nu.f", PGG_nu | PGG_f | PGG_LIF | PGG_UHR  | PGG_TOTAL , 0, "LIF_R.UHR.nu.f", 0, 0, GREEN6, TRUE, 0, 0, 3000000, 0},


  /* unique */
 
{"ILM_S.Brain.u.f", PGG_u | PGG_f | PGG_ILM | PGG_Brain  | PGG_rna , 0, "ILM_S.Brain.u.f", 0, 0, GREEN5, TRUE, 0, 0, 3000000, 0},
{"LIF_S.Brain.u.f", PGG_u | PGG_f | PGG_LIF | PGG_Brain  | PGG_rna , 0, "LIF_S.Brain.u.f", 0, 0, GREEN6, TRUE, 0, 0, 3000000, 0},
{"HELdge.Brain.u.f", PGG_u | PGG_f | PGG_HEL | PGG_Brain | PGG_rna , 0, "HELdge.Brain.u.f", 0, 0, GREEN5, TRUE, 0, 0, 3000000, 0},

{"ILM_S.UHR.u.f", PGG_u | PGG_f | PGG_ILM | PGG_UHR | PGG_rna  , 0, "ILM_S.UHR.u.f", 0, 0, RED4, TRUE, 0, 0, 3000000, 0},
{"LIF_S.UHR.u.f", PGG_u | PGG_f | PGG_LIF | PGG_UHR | PGG_rna  , 0, "LIF_S.UHR.u.f", 0, 0, RED6, TRUE, 0, 0, 3000000, 0},
{"HELdge.UHR.u.f", PGG_u | PGG_f | PGG_HEL | PGG_UHR  | PGG_rna , 0, "HELdge.UHR.u.f", 0, 0, RED5, TRUE, 0, 0, 3000000, 0},


{"ILM_S.Brain.u.r", PGG_u | PGG_r | PGG_ILM | PGG_Brain  | PGG_rna , 0, "ILM_S.Brain.u.r", 0, 0, BLUE5, TRUE, 0, 0, 3000000, 0},
{"LIF_S.Brain.u.r", PGG_u | PGG_r | PGG_LIF | PGG_Brain  | PGG_rna , 0, "LIF_S.Brain.u.r", 0, 0, BLUE7, TRUE, 0, 0, 3000000, 0},
{"HELdge.Brain.u.r", PGG_u | PGG_r | PGG_HEL | PGG_Brain | PGG_rna , 0, "HELdge.Brain.u.r", 0, 0, BLUE6, TRUE, 0, 0, 3000000, 0},

{"ILM_S.UHR.u.r", PGG_u | PGG_r | PGG_ILM | PGG_UHR | PGG_rna  , 0, "ILM_S.UHR.u.r", 0, 0, LIGHTVIOLET, TRUE, 0, 0, 3000000, 0},
{"LIF_S.UHR.u.r", PGG_u | PGG_r | PGG_LIF | PGG_UHR | PGG_rna  , 0, "LIF_S.UHR.u.r", 0, 0, VIOLET, TRUE, 0, 0, 3000000, 0},
{"HELdge.UHR.u.r", PGG_u | PGG_r | PGG_HEL | PGG_UHR  | PGG_rna , 0, "HELdge.UHR.u.r", 0, 0, LIGHTVIOLET, TRUE, 0, 0, 3000000, 0},



{"ILM_nS.Brain.u.ns", PGG_u | PGG_ns | PGG_ILM | PGG_Brain  | PGG_rna , 0, "ILM_nS.Brain.u.ns", 0, 0, GREEN1, TRUE, 0, 0, 3000000, 0},
{"HEL.Brain.u.f", PGG_u | PGG_f | PGG_HEL | PGG_Brain  | PGG_rna , 0, "HEL.Brain.u.f", 0, 0, GREEN2, TRUE, 0, 0, 3000000, 0},
{"HEL.Brain.u.r", PGG_u | PGG_r | PGG_HEL | PGG_Brain  | PGG_rna , 0, "HEL.Brain.u.r", 0, 0, GREEN2, TRUE, 0, 0, 3000000, 0},
{"R454.Brain.u.ns", PGG_u | PGG_ns | PGG_R454 | PGG_Brain  | PGG_rna , 0, "R454.Brain.u.ns", 0, 0, GREEN1, TRUE, 0, 0, 3000000, 0},


{"ILM_nS.UHR.u.ns", PGG_u | PGG_ns | PGG_ILM | PGG_UHR | PGG_rna  , 0, "ILM_nS.UHR.u.ns", 0, 0, RED1, TRUE, 0, 0, 3000000, 0},
{"HEL.UHR.u.f", PGG_u | PGG_f | PGG_HEL | PGG_UHR | PGG_rna  , 0, "HEL.UHR.u.f", 0, 0, RED2, TRUE, 0, 0, 3000000, 0},
{"HEL.UHR.u.r", PGG_u | PGG_r | PGG_HEL | PGG_UHR | PGG_rna  , 0, "HEL.UHR.u.r", 0, 0, RED2, TRUE, 0, 0, 3000000, 0},
{"R454.UHR.u.ns", PGG_u | PGG_ns | PGG_R454 | PGG_UHR | PGG_rna  , 0, "R454.UHR.u.ns", 0, 0, RED1, TRUE, 0, 0, 3000000, 0},



{"ILM_RS.Brain.u.f", PGG_u | PGG_f | PGG_ILM | PGG_Brain | PGG_TOTAL , 0, "ILM_RS.Brain.u.f", 0, 0, LIGHTGREEN, TRUE, 0, 0, 3000000, 0},
{"LIF_R.Brain.u.f", PGG_u | PGG_f | PGG_LIF | PGG_Brain | PGG_TOTAL , 0, "LIF_R.Brain.u.f", 0, 0, DARKGREEN, TRUE, 0, 0, 3000000, 0},
{"ILM_RS.UHR.u.f", PGG_u | PGG_f | PGG_ILM | PGG_UHR  | PGG_TOTAL , 0, "ILM_RS.UHR.u.f", 0, 0, DARKRED, TRUE, 0, 0, 3000000, 0},


{"ILM_RS.Brain.u.r", PGG_u | PGG_r | PGG_ILM | PGG_Brain | PGG_TOTAL , 0, "ILM_RS.Brain.u.r", 0, 0, CYAN, TRUE, 0, 0, 3000000, 0},
{"LIF_R.Brain.u.r", PGG_u | PGG_r | PGG_LIF | PGG_Brain | PGG_TOTAL , 0, "LIF_R.Brain.u.r", 0, 0, DARKCYAN, TRUE, 0, 0, 3000000, 0},
{"ILM_RS.UHR.u.r", PGG_u | PGG_r | PGG_ILM | PGG_UHR  | PGG_TOTAL , 0, "ILM_RS.UHR.u.r", 0, 0, DARKVIOLET, TRUE, 0, 0, 3000000, 0},



{"ILM_RnS.Brain.u.ns", PGG_u | PGG_ns | PGG_ILM | PGG_Brain | PGG_TOTAL , 0, "ILM_RnS.Brain.u.ns", 0, 0, LAVANDER, TRUE, 0, 0, 3000000, 0},
{"ILM_RnS.UHR.u.ns", PGG_u | PGG_ns | PGG_ILM | PGG_UHR  | PGG_TOTAL , 0, "ILM_RnS.UHR.u.ns",  0, 0, BROWN, TRUE, 0, 0, 3000000, 0},





#endif

    { 0, 0, 0, 0, 10.0, 10.0, BLACK, FALSE, FALSE, 0, 0, 0},
    { 0, 0, 0, 0, 10.0, 10.0, BLACK, FALSE, FALSE, 0, 0, 0}
  } ;
#endif

#ifdef BOU_PNS

static PNS pnsAll[] =
  {

#ifdef BOU_DEF
   /* total AFFY MS */
    { "AFX_3188", PGG_AFFY | PGG_MS,  32, "pA=XMS", "pA=XMS", 5.76553, 13.22, RED, FALSE, 2, 9.5, 11.0, 0} ,
    { ".AFX_3188", 0, 32, ".pA=XMS", "pA=XMS", 5.890, 13.22, RED, FALSE, FALSE, 9.5, 11.0, 0} ,

    { "AFX_3330", PGG_AFFY | PGG_MS, 32, "pA=XMSb", "pA=XMS", 5.6410, 13.22, RED, FALSE, 2, 9.5, 11.0, 0} ,
    { ".AFX_3330", 0, 32, ".pA=XMSb", "pA=XMS", 5.641, 13.22, RED, FALSE, FALSE, 9.5, 11.0, 0} ,

    { "AFX_3331", PGG_AFFY | PGG_MS, 53, "pA=XstarvedMS", "pA=XstarvedMS", 6.2000, 13.22, LIGHTRED, FALSE, 2, 9.5, 11.0, 0} ,
    { ".AFX_3331", 0, 53, ".pA=XstarvedMS", "pA=Xs", 6.21010, 13.22, LIGHTRED, FALSE, FALSE, 9.5, 11.0, 0} ,

     /* pA MS */
    { "86687_635", PGG_PA | PGG_MS, 38, "pA=12MS", "pA/G1=12MS", 8.0635, 13.27, RED, FALSE, 2, 11.15, 10.44, 0} ,
    { "86687_532", PGG_G1 | PGG_MS, 430, "G1=12MS",  "ratio", 10.6025, 12.70, BLUE3, FALSE, FALSE, 12.61, 11.79, 0} ,

    { "92677_635", PGG_PA | PGG_MS, 128, "pA=22MS", "pA/G1=22MS", 9.2126, 14.82, RED, FALSE, 2, 11.575, 11.48, 0} ,
    { "92677_532", PGG_G1 | PGG_MS, 90, "G1=22MS", "ratio", 10.1240, 12.87, BLUE4, FALSE, FALSE, 12.44, 11.64, 0} ,

    /* exchange 144047=22b really is 23b
     *          140787=23b really is 22b
     */
    {"140787_635", PGG_PA | PGG_MS, 128, "pA=22bMS", "pA/G1=22bMS", 9.1919, 13.35, RED, FALSE, 2, 12.04, 11.69, 0} ,
    {"140787_532", PGG_G1 | PGG_MS, 512, "G1=22bMS", "ratio", 11.0153, 13.31, BLUE5, FALSE, FALSE, 12.9, 12.17, 0} ,

    /* total/pA MS */
    { "84964_532", PGG_PA | PGG_MS, 2, "pA=11MS",  "pA/total=11MS",  7.577, 13.95, RED, FALSE, 2, 11.0, 10.38, 0} ,
    { "84964_635", PGG_TOTAL | PGG_MS, 152, "total=11MS", "ratio", 10.1606, 13.77, VIOLET, FALSE, 1, 12.21, 11.94, 0} ,

    { "103738_532", PGG_PA | PGG_MS, 38, "pA=21MS","pA/total=21MS",  8.9014, 14.72, RED , FALSE, 2, 11.645, 11.49, 0} ,
    { "103738_635", PGG_TOTAL | PGG_MS, 362, "total=21MS", "ratio",10.1819, 13.75, VIOLET, FALSE, 1, 12.20, 11.92, 0} ,

     /* total/G1 MS */
    { "86948_635", PGG_TOTAL | PGG_MS, 181, "total=14MS",  "total/G1=14MS", 9.644, 13.65, VIOLET, FALSE, 1, 11.85, 11.61, 0} ,
    { "86948_532", PGG_G1 | PGG_MS, 215, "G1=14MS",  "ratio", 10.682, 12.792, BLUE3, FALSE, FALSE, 12.85, 11.86, 0} ,

    { "94646_635", PGG_TOTAL | PGG_MS, 304, "total=24MS", "total/G1=24MS", 10.4397, 14.04, VIOLET, FALSE, 1, 12.52, 12.09, 0} ,
    { "94646_532", PGG_G1 | PGG_MS, 215, "G1=24MS",  "ratio", 10.5065, 13.45, BLUE4, FALSE, FALSE, 12.70, 12.10, 0} ,

    {"144049_635", PGG_TOTAL | PGG_MS, 304, "total=24bMS", "total/G1=24bMS", 9.8283, 13.45, DARKGREEN /*VIOLET*/, FALSE, 1, 11.83, 11.6, 0} ,    
    {"144049_532", PGG_G1 | PGG_MS, 27, "G1=24bMS",  "ratio", 10.5361, 13.28, BLUE5, FALSE, FALSE, 12.87, 12.13, 0} ,

    /* pA-/G1 MS */
    { "86827_635", PGG_NUC | PGG_MS, 152, "Nuc=13MS",  "Nuc/G1=13MS",  9.904, 12.98, BROWN, FALSE, 1, 12.02, 11.56, 0} ,
    { "86827_532", PGG_G1 | PGG_MS, 430, "G1=13MS", "ratio", 10.9145, 12.74, BLUE3, FALSE, FALSE, 12.91, 12.04, 0} ,

    { "92680_635", PGG_NUC | PGG_MS, 608, "Nuc=23MS",  "Nuc/G1=23MS", 10.9997, 14.23, BROWN, FALSE, 1, 12.79, 12.57, 0} ,
    { "92680_532", PGG_G1 | PGG_MS, 152, "G1=23MS",  "ratio", 10.5706, 13.67, BLUE4, FALSE, FALSE, 12.73, 12.16, 0} ,

    {"144047_635", PGG_NUC | PGG_MS, 304, "Nuc=23bMS",  "Nuc/G1=23bMS", 9.98, 14.54,  BROWN, FALSE, 1, 11.60, 11.71, 0} ,
    {"144047_532", PGG_G1 | PGG_MS, 64, "G1=23bMS",  "ratio", 10.5711, 14.36, BLUE5, FALSE, FALSE, 12.95, 12.86, 0} ,

    /* S/G1 MS */
    { "86696_635", PGG_S | PGG_MS, 304, "S=15MSr",   "S/G1=15MSr/v",  10.7428, 12.85,  LIGHTCYAN, FALSE, 0, 11.04, 12.07, 0} ,
    { "86696_532", PGG_G1 | PGG_MS, 430, "G1=15MSv",  "ratio", 10.8942, 12.798, BLUE3, FALSE, FALSE, 12.900, 12.00, 0} , 

    { "92676_635", PGG_S | PGG_MS, 215, "S=25MS",  "S/G1=25MS", 10.2622, 13.25, LIGHTCYAN, FALSE, FALSE, 10.7, 12.05, 0} ,
    { "92676_532", PGG_G1 | PGG_MS, 215, "G1=25MS",  "ratio", 10.5958, 13.04, BLUE4, FALSE, FALSE, 12.81, 12.01, 0} , 

    { "148457_635", PGG_S | PGG_MS, 430, "S=28MS",  "S/G1=28MS", 10.8563, 13.50, CYAN, FALSE , FALSE, 11.12, 12.37, 0},
    { "148457_532", PGG_G1 | PGG_MS, 256, "G1=28MS",  "ratio", 11.0832 , 13.51,  BLUE5, FALSE , FALSE, 13.20, 12.39, 0},

 
    { "178412_635", 0 && (PGG_S | PGG_MS), 430, "S=29MS",  "S/G1=29MS", 10.7124, 13.50, CYAN, FALSE , FALSE, 11.12, 12.37, 0},
    { "178412_532", 0 && (PGG_G1 | PGG_MS), 724, "G1=29MS",  "ratio",  11.8261, 13.51, BLUE5, FALSE , FALSE, 13.20, 12.39, 0},

    { "180576_635", PGG_S | PGG_MS, 1024, "S=30MS", "S/G1=30MS" ,  11.55, 13.63, DARKCYAN, FALSE, FALSE, 12.60, 12.38, 0} ,
    { "180576_532", PGG_G1 | PGG_MS, 304, "G1=30MS" , "ratio",  12.02, 13.58, BLUE7, FALSE, FALSE, 12.93, 12.37, 0} ,

    { "101320_635", PGG_S | PGG_MS, 304, "Samplified=26MS", "Samp/G1amp=26MS" ,  10.453, 13.63, LIGHTCYAN, FALSE, FALSE, 12.60, 12.38, 0} ,
    { "101320_532", PGG_G1 | PGG_MS, 304, "G1amplified=26MS" , "ratio",  10.7727, 13.58, BLUE6, FALSE, FALSE, 12.93, 12.37, 0} ,


   /* G1/G1 MS */
    { "86701_635", PGG_G1 | PGG_MS, 181, "G1=16MSr", "G1/G1=16MS", 10.5253, 12.59, BLUE3, FALSE, FALSE, 11.04, 11.75, 0} ,
    { "86701_532", PGG_G1 | PGG_MS, 256, "G1=16MSv",  "ratio", 10.5117, 12.5, BLUE3, FALSE, FALSE, 12.68, 11.59, 0} ,

    { "87300_635", PGG_G1 | PGG_MS, 724, "G1=16bMSr", "G1/G1=16bMS",  11.5314, 13.44, BLUE2, FALSE, FALSE, 11.04, 12.67, 0} ,
    { "87300_532", PGG_G1 | PGG_MS, 1024, "G1=16bMSv", "ratio", 11.5672, 13.36, BLUE2, FALSE, FALSE, 13.30, 12.58, 0} ,

    { "98731_635", PGG_G1 | PGG_MS, 256, "G1=27MSr", "G1/G1=27MS",  10.5431, 13.22, BLUE4, FALSE, FALSE, 12.60, 12.04, 0} ,
    { "98731_532", PGG_G1 | PGG_MS, 256, "G1=27MSv", "ratio",  10.6066, 13.18, BLUE4, FALSE, FALSE, 12.775, 12.04, 0} ,

    { "228322_532", PGG_G1 | PGG_RED, 181, "G1=50Crd", "G1/G1=50Crd/MStrvd" ,  11.988, 13.63, DARKGREEN, FALSE, FALSE, 12.60, 12.38, 0} ,
    { "228322_635", PGG_G1 | PGG_MS, 861, "G1=50MStrvd", "ratio" ,  11.68, 13.63,  BLUE3, FALSE, FALSE, 12.60, 12.38, 0} ,

   /***********************************/
    /* AFX ES */
    { "AFX_3329", PGG_AFFY | PGG_ES, 32, "pA=XES", "pA=XES", 5.7469, 13.22,  CERISE, FALSE, 2, 9.5, 11.0, 0} ,
    { ".AFX_3329", 0, 32, ".pA=XES", "pA=XES", 5.0, 13.22,  CERISE, FALSE, FALSE, 9.5, 11.0, 0} ,


    /* pA ES */
    { "197152_635", PGG_PA | PGG_ES, 304, "pA=33ESr", "pA/G1=33ES", 9.3952, 13.22, CERISE, FALSE, 2, 12.60, 12.04, 0} ,
    { "197152_532", PGG_G1 | PGG_ES, 304, "G1=33ESv", "ratio", 11.0103, 13.18, LAVANDER, FALSE, FALSE, 12.775, 12.04, 0} ,

    { "228324_532", PGG_PA | PGG_ES, 256, "pA=42ES", "pA/G1=42ES" ,  9.31, 13.13, DARKRED, FALSE, FALSE, 12.60, 12.38, 0} ,
    { "228324_635", PGG_G1 | PGG_ES, 861, "G1=42ES", "ratio" ,  11.3072, 13.63, LAVANDER, FALSE, FALSE, 12.60, 12.38, 0} ,

   /* total/G1 ES */
    { "198378_635", PGG_TOTAL | PGG_ES, 256, "total=35ESr", "total/G1=35ES",  10.4279, 13.22, LIGHTMAGENTA, FALSE, 1, 12.60, 12.04, 0} ,
    { "198378_532", PGG_G1 | PGG_ES, 304, "G1=35ESv", "ratio", 10.8712, 13.18, LAVANDER, FALSE, FALSE, 12.775, 12.04, 0} ,

    { "228323_532", PGG_TOTAL | PGG_ES, 256, "total=44ES", "total/G1=44ES" ,  9.75, 13.63, LIGHTMAGENTA, FALSE, FALSE, 12.60, 12.38, 0} ,
    { "228323_635", PGG_G1 | PGG_ES, 609, "G1=44ES", "ratio" ,  10.92, 13.63, LAVANDER, FALSE, FALSE, 12.60, 12.38, 0} ,

    /* pA-/G1 ES */
    /* 34 est degoutante */
    { "197157_635", 0, 1722, "Nuc=34ESr", "Nuc/G1=34ES", 12.05, 13.22, ORANGE, FALSE, 1, 12.60, 12.04, 0} ,
    { "197157_532", 0, 256, "G1=34ESv", "ratio",  11.2249, 13.18, LAVANDER, FALSE, FALSE, 12.775, 12.04, 0} ,

    { ".228330_532", PGG_NUC | PGG_ES, 256, "Nuc=43ES", "Nuc/G1=43ES" ,  11.0, 13.63,  ORANGE, FALSE, FALSE, 12.60, 12.38, 0} ,
    { ".228330_635", PGG_G1 | PGG_ES, 1000, "G1=43ES", "ratio" ,  11.773, 13.63, LAVANDER, FALSE, FALSE, 12.60, 12.38, 0} ,
    /* S/G1 ES */
    {"140789_635", PGG_S | PGG_ES, 512, "S=31ES",   "S/G1=31ES", 11.0025 , 13.61,  LIGHTGREEN, FALSE, FALSE, 11.06, 12.61, 0} ,
    {"140789_532", PGG_G1 | PGG_ES, 512, "G1=31ES",   "ratio", 11.1942, 13.54, LAVANDER, FALSE, FALSE, 12.98, 12.53, 0} ,
    
    { "226071_532", PGG_S | PGG_ES, 128, "S=41ES", "S/G1=41ES" ,  11.3273, 13.63, DARKGREEN, FALSE, FALSE, 12.60, 12.38, 0} ,
    { "226071_635", PGG_G1 | PGG_ES, 1024, "G1=41ES", "ratio" ,  11.6531, 13.63, LAVANDER, FALSE, FALSE, 12.60, 12.38, 0} ,

    /* G1/G1 ES */
    { "194537_635", PGG_G1 | PGG_ES, 608, "G1=32ESr", "G1/G1=32ES", 11.2665, 13.22, LAVANDER, FALSE, FALSE, 12.60, 12.04, 0} ,
    { "194537_532", PGG_G1 | PGG_ES, 304, "G1=32ESv", "ratio", 11.325, 13.18, LAVANDER, FALSE, FALSE, 12.775, 12.04, 0} ,

    /******************************************** ERY ******/
    { "AFX_3333", PGG_AFFY | PGG_RED, 32, "pA=XERY", "pA=XERY", 5.9402, 13.22, ORANGE, FALSE, 2, 9.5, 11.0, 0} ,
    { ".AFX_3333", 0, 32, ".pA=XERY", "pA=XERY", 5.0, 13.22, ORANGE, FALSE, FALSE, 9.5, 11.0, 0} ,

    /* PA/G1 RED */
    {"229118_532", PGG_PA | PGG_RED, 304, "pA=49Red",  "pA/G1=49Red", 10, 14.54,  BROWN, FALSE, 1, 11.60, 11.71, 0} ,
    {"229118_635", PGG_G1 | PGG_RED, 724, "G1=49Red",  "ratio", 11.0095, 14.36, BLUE5, FALSE, FALSE, 12.95, 12.86, 0} ,

    /* Total/G1 RED */
    {"231191_532", PGG_TOTAL | PGG_RED, 304, "total=48Red",  "total/G1=48Red", 10.8, 14.54,  BROWN, FALSE, 1, 11.60, 11.71, 0} ,
    {"231191_635", PGG_G1 | PGG_RED, 608, "G1=48Red",  "ratio", 10.9326, 14.36, BLUE5, FALSE, FALSE, 12.95, 12.86, 0} ,

    /* Nuc/G1 RED */
    {"233922_532", PGG_NUC | PGG_RED, 304, "Nuc=47Red",  "Nuc/G1=47Red", 10.15, 14.54,  BROWN, FALSE, 1, 11.60, 11.71, 0} ,
    {"233922_635", PGG_G1 | PGG_RED, 724, "G1=47Red",  "ratio", 11.1507, 14.36, BLUE5, FALSE, FALSE, 12.95, 12.86, 0} ,


    /* S/G1 RED */
    {"194302_635", PGG_S | PGG_RED, 512, "S=36Red",   "S/G1=36Red", 10.7854 , 13.61,  GREEN6, FALSE, FALSE, 11.06, 12.61, 0} ,
    {"194302_532", PGG_G1 | PGG_RED, 362, "G1=36Red",   "ratio", 11.7133, 13.54, MIDBLUE, FALSE, FALSE, 12.98, 12.53, 0} ,

    {"194272_635", PGG_S | PGG_RED, 512, "S=37Red",   "S/G1=37Red", 10.7822 , 13.61,  GREEN8, FALSE, FALSE, 11.06, 12.61, 0} ,
    {"194272_532", PGG_G1 | PGG_RED, 256, "G1=37Red",   "ratio", 11.545, 13.54, MIDBLUE, FALSE, FALSE, 12.98, 12.53, 0} ,
  
    { "232815_635", PGG_S | PGG_RED | PGG_DYNA, 608, "S1=45Red", "S1/G1=45Red" ,  10.758, 13.63, GREEN3, FALSE, FALSE, 12.60, 12.38, 0} ,
    { "232815_532", PGG_G1 | PGG_RED, 362, "G1=45Red", "ratio" ,  11.3148, 11.63,  MIDBLUE, FALSE, FALSE, 12.60, 12.38, 0} ,

    { "226089_532", PGG_S | PGG_RED | PGG_DYNA, 128, "S2=38Red", "S2/G1=38Red" ,  11.055, 13.63, GREEN5, FALSE, FALSE, 12.60, 12.38, 0} ,
    { "226089_635", PGG_G1 | PGG_RED, 861, "G1=38Red", "ratio" ,  11.3788, 11.63,  MIDBLUE, FALSE, FALSE, 12.60, 12.38, 0} ,
   
    { "233914_635", PGG_S | PGG_RED | PGG_DYNA, 724, "S3=46Red", "S3/G1=46Red" ,  10.882, 13.63, GREEN7, FALSE, FALSE, 12.60, 12.38, 0} ,
    { "233914_532", PGG_G1 | PGG_RED, 512, "G1=46Red", "ratio" ,  11.7514, 11.63,  MIDBLUE, FALSE, FALSE, 12.60, 12.38, 0} ,
   
    { "226091_532", PGG_S | PGG_RED , 27, "S123=39Red", "S123/G1=39Red" ,  11.0308, 13.63, GREEN, FALSE, FALSE, 12.60, 12.38, 0} ,
    { "226091_635", PGG_G1 | PGG_RED,430, "G1=39Red", "ratio" ,  10.81, 10.63,  MIDBLUE, FALSE, FALSE, 12.60, 12.38, 0} ,

    /* S/G1 VARIOUS */
    {"233372_635", PGG_S | PGG_RED, 512, "S=51FLiv",   "S/G1=51FLiv", 10.62 , 13.61,  GREEN6, FALSE, FALSE, 11.06, 12.61, 0} ,
    {"233372_532", PGG_G1 | PGG_RED, 512, "G1=51FLiv",   "ratio", 11.6232, 13.54, MIDBLUE, FALSE, FALSE, 12.98, 12.53, 0} ,

    {"232672_635", PGG_S | PGG_BM, 608, "S=52BM",   "S/G1=52BM", 10.745, 13.61,  GREEN6, FALSE, FALSE, 11.06, 12.61, 0} ,
    {"232672_532", PGG_G1 | PGG_BM, 152, "G1=52BM",   "ratio", 11.14, 13.54, MIDBLUE, FALSE, FALSE, 12.98, 12.53, 0} ,

    /* G1/G1 ES */
     /***************************************************/
#ifdef JUNK
    /* Nath Addition of one track for affy 3683 */

    { "AFX_3683", PGG_AFFY | PGG_ES, 24, "X_ESNL", "X_ESNL", 8.5, 13.22, BLACK, FALSE, 2, 9.5, 11.0, 0} ,
    { ".AFX_3683", 0, 24, ".X_ESNL", "X_ESNL", 8.5, 13.22, BLACK, FALSE, FALSE, 9.5, 11.0, 0} ,
    

    /* Nath Addition of  tracks for transfrag visualisation 144047 */ 

    /* Non exonic probes only for exp 144047=23bMS */

    { "N_transfrag144047-NonExonic-med10-0.25", 1, "Tranf23bMS-NonExonic-med10-0.25", "Tranf23bMS-NonExonic-med10-0.25", 0, 20, ORANGE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag144047-NonExonic-med10-0.25", 1, ".Tranf23bMS-NonExonic-med10-0.25", "Tranf23bMS-NonExonic-med10-0.25",0, 20, ORANGE, FALSE,  FALSE, 9.5, 11.0, 0},
    
    { "N_transfrag144047-NonExonic-med15-0.25", 0, "Tranf23bMS-NonExonic-med15-0.25", "Tranf23bMS-NonExonic-med15-0.25", 0, 20, GREEN, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag144047-NonExonic-med15-0.25", 0, ".Tranf23bMS-NonExonic-med15-0.25", "Tranf23bMS-NonExonic-med15-0.25",0, 20, GREEN, FALSE,  FALSE, 9.5, 11.0, 0},

    { "N_transfrag144047-NonExonic-med20-0.25", -1, "Tranf23bMS-NonExonic-med20-0.25", "Tranf23bMS-NonExonic-med20-0.25", 0, 20, MIDBLUE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag144047-NonExonic-med20-0.25", -1, ".Tranf23bMS-NonExonic-med20-0.25", "Tranf23bMS-NonExonic-med20-0.25",0, 20, MIDBLUE, FALSE,  FALSE, 9.5, 11.0, 0},


            /* End of Non exonic probes only for exp 144047=23bMS */

    { "N_transfrag144047-med15-0.25", 3, "Tranf23bMS-med15-0.25", "Tranf23bMS-med15-0.25", 0, 20, RED, FALSE, 2, 9.5, 11.0, 0} ,
    { ".N_transfrag144047-med15-0.25", 3, ".Tranf23bMS-med15-0.25", "Tranf23bMS-med15-0.25",0, 20, RED, FALSE, FALSE, 9.5, 11.0, 0} ,

    { "N_transfrag144047-med15-0.25-filtered-90", 2, "Tranf23bMS-med15-0.25-filtered-90", "Tranf23bMS-med15-0.25-filtered-90", 0, 20, BLUE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag144047-med15-0.25-filtered-90", 2, ".Tranf23bMS-med15-0.25-filtered-90", "Tranf23bMS-med15-0.25-filtered-90",0, 20, BLUE, FALSE,  FALSE, 9.5, 11.0, 0},
    { "N_transfrag144047-med15-0.25-filtered-80", 1, "Tranf23bMS-med15-0.25-filtered-80", "Tranf23bMS-med15-0.25-filtered-80", 0, 20, GREEN, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag144047-med15-0.25-filtered-80", 1, ".Tranf23bMS-med15-0.25-filtered-80", "Tranf23bMS-med15-0.25-filtered-80",0, 20, GREEN, FALSE,  FALSE, 9.5, 11.0, 0},
    { "N_transfrag144047-med15-0.25-filtered-70", 0, "Tranf23bMS-med15-0.25-filtered-70", "Tranf23bMS-med15-0.25-filtered-70", 0, 20, VIOLET, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag144047-med15-0.25-filtered-70", 0, ".Tranf23bMS-med15-0.25-filtered-70", "Tranf23bMS-med15-0.25-filtered-70",0, 20, VIOLET, FALSE,  FALSE, 9.5, 11.0, 0},
    { "N_transfrag144047-med15-0.25-filtered-60", -1, "Tranf23bMS-med15-0.25-filtered-60", "Tranf23bMS-med15-0.25-filtered-60", 0, 20, MAGENTA, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag144047-med15-0.25-filtered-60", -1, ".Tranf23bMS-med15-0.25-filtered-60", "Tranf23bMS-med15-0.25-filtered-60",0, 20, MAGENTA, FALSE,  FALSE, 9.5, 11.0, 0},
    { "N_transfrag144047-med15-0.25-filtered-50", -2, "Tranf23bMS-med15-0.25-filtered-50", "Tranf23bMS-med15-0.25-filtered-50", 0, 20, ORANGE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag144047-med15-0.25-filtered-50", -2, ".Tranf23bMS-med15-0.25-filtered-50", "Tranf23bMS-med15-0.25-filtered-50",0, 20, ORANGE, FALSE,  FALSE, 9.5, 11.0, 0},

    /*{ "N_transfrag144047-med20-0.25", 1, "Tranf23bMS-med20-0.25", "Tranf23bMS-med20-0.25", 0, 20, VIOLET, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag144047-med20-0.25", 1, ".Tranf23bMS-med20-0.25", "Tranf23bMS-med20-0.25",0, 20, VIOLET, FALSE,  FALSE, 9.5, 11.0, 0},

    { "N_transfrag92680-med15-0.35", 0, "Tranf23MS-med15-0.35", "Tranf23MS-med15-0.35", 0, 20, RED, FALSE, 2, 9.5, 11.0, 0} ,
    { ".N_transfrag92680-med15-0.35", 0, ".Tranf23MS-med15-0.35", "Tranf23MS-med15-0.35",0, 20, RED, FALSE,  FALSE, 9.5, 11.0, 0},*/
            

    /* End of Nath Addition of  tracks for transfrag visualisation 144047 */

    /* Nath Addition of  tracks for transfrag visualisation 86827 */ 

    /* Non exonic probes only for exp 86827=13MS */

    { "N_transfrag86827-NonExonic-med10-0.25", 1, "Tranf13MS-NonExonic-med10-0.25", "Tranf13MS-NonExonic-med10-0.25", 0, 20, ORANGE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag86827-NonExonic-med10-0.25", 1, ".Tranf13MS-NonExonic-med10-0.25", "Tranf13MS-NonExonic-med10-0.25",0, 20, ORANGE, FALSE,  FALSE, 9.5, 11.0, 0},
    
    { "N_transfrag86827-NonExonic-med15-0.25", 0, "Tranf13MS-NonExonic-med15-0.25", "Tranf13MS-NonExonic-med15-0.25", 0, 20, GREEN, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag86827-NonExonic-med15-0.25", 0, ".Tranf13MS-NonExonic-med15-0.25", "Tranf13MS-NonExonic-med15-0.25",0, 20, GREEN, FALSE,  FALSE, 9.5, 11.0, 0},

    { "N_transfrag86827-NonExonic-med20-0.25", -1, "Tranf13MS-NonExonic-med20-0.25", "Tranf13MS-NonExonic-med20-0.25", 0, 20, MIDBLUE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag86827-NonExonic-med20-0.25", -1, ".Tranf13MSS-NonExonic-med20-0.25", "Tranf13MS-NonExonic-med20-0.25",0, 20, MIDBLUE, FALSE,  FALSE, 9.5, 11.0, 0},


    /* End of Non exonic probes only for exp 86827=13MS */

    { "N_transfrag86827-med15-0.25", 1, "Tranf13MS-med15-0.25", "Tranf13MS-med15-0.25", 0, 20, CERISE, FALSE, 2, 9.5, 11.0, 0} ,
    { ".N_transfrag86827-med15-0.25", 1, ".Tranf13MS-med15-0.25", "Tranf13MS-med15-0.25",0, 20, CERISE, FALSE, FALSE, 9.5, 11.0, 0} ,

    { "N_transfrag86827-med10-0.25", 0, "Tranf13MS-med10-0.25", "Tranf13MS-med10-0.25", 0, 20, MIDBLUE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag86827-med10-0.25", 0, ".Tranf13MS-med10-0.25", "Tranf13MS-med10-0.25",0, 20, MIDBLUE, FALSE,  FALSE, 9.5, 11.0, 0},

    { "N_transfrag86827-med20-0.25", 2, "Tranf13MS-med20-0.25", "Tranf13MS-med20-0.25", 0, 20, LAVANDER, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag86827-med20-0.25", 2, ".Tranf13MS-med20-0.25", "Tranf13MS-med20-0.25",0, 20, LAVANDER, FALSE,  FALSE, 9.5, 11.0, 0},
            

    /* End of Nath Addition of  tracks for transfrag visualisation 86827 */


    /* Nath Addition of  tracks for transfrag visualisation 92680 */ 
    
    /* Non exonic probes only for exp 92680=23MS */

    { "N_transfrag92680-NonExonic-med10-0.25", 1, "Tranf23MS-NonExonic-med10-0.25", "Tranf23MS-NonExonic-med10-0.25", 0, 20, ORANGE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag92680-NonExonic-med10-0.25", 1, ".Tranf23MS-NonExonic-med10-0.25", "Tranf23MS-NonExonic-med10-0.25",0, 20, ORANGE, FALSE,  FALSE, 9.5, 11.0, 0},
    
    { "N_transfrag92680-NonExonic-med15-0.25", 0, "Tranf23MS-NonExonic-med15-0.25", "Tranf23MS-NonExonic-med15-0.25", 0, 20, GREEN, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag92680-NonExonic-med15-0.25", 0, ".Tranf23MS-NonExonic-med15-0.25", "Tranf23MS-NonExonic-med15-0.25",0, 20, GREEN, FALSE,  FALSE, 9.5, 11.0, 0},

    { "N_transfrag92680-NonExonic-med20-0.25", -1, "Tranf23MS-NonExonic-med20-0.25", "Tranf23MS-NonExonic-med20-0.25", 0, 20, MIDBLUE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag92680-NonExonic-med20-0.25", -1, ".Tranf23MSS-NonExonic-med20-0.25", "Tranf23MS-NonExonic-med20-0.25",0, 20, MIDBLUE, FALSE,  FALSE, 9.5, 11.0, 0},


    /* End of Non exonic probes only for exp 92680=23MS */

    { "N_transfrag92680-med15-0.25", 2, "Tranf23MS-med15-0.25", "Tranf23MS-med15-0.25", 0, 20, LIGHTRED, FALSE, 2, 9.5, 11.0, 0} ,
    { ".N_transfrag92680-med15-0.25", 2, ".Tranf23MS-med15-0.25", "Tranf23MS-med15-0.25",0, 20, LIGHTRED, FALSE, FALSE, 9.5, 11.0, 0} ,

    { "N_transfrag92680-med10-0.25", 1, "Tranf23MS-med10-0.25", "Tranf23MS-med10-0.25", 0, 20, CYAN, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag92680-med10-0.25", 1, ".Tranf23MS-med10-0.25", "Tranf23MS-med10-0.25",0, 20, CYAN, FALSE,  FALSE, 9.5, 11.0, 0},

    { "N_transfrag92680-med20-0.25", 0, "Tranf23MS-med20-0.25", "Tranf23MS-med20-0.25", 0, 20, PURPLE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag92680-med20-0.25", 0, ".Tranf23MS-med20-0.25", "Tranf23MS-med20-0.25",0, 20, PURPLE, FALSE,  FALSE, 9.5, 11.0, 0},
            

    /* End of Nath Addition of  tracks for transfrag visualisation 92680 */


    /* Nath addition of random buttons*/
    /*HBA 23bMS*/

    { "N_trfHBARandom2-med15-0.25", 2, "HBA-Rd2-med15-0.25", "HBA-Rd2-med15-0.25", 0, 20, LIGHTRED, FALSE, 2, 9.5, 11.0, 0} ,
    { ".N_trfHBARandom2-med15-0.25", 2, ".HBA-Rd2-med15-0.25", "HBA-Rd2-med15-0.25",0, 20, LIGHTRED, FALSE, FALSE, 9.5, 11.0, 0} ,

    { "N_trfHBARandom4-med15-0.25", 1, "HBA-Rd4-med15-0.25", "HBA-Rd4-med15-0.25", 0, 20, CYAN, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_trfHBARandom4-med15-0.25", 1, ".HBA-Rd4-med15-0.25", "HBA-Rd4-med15-0.25",0, 20, CYAN, FALSE,  FALSE, 9.5, 11.0, 0},

    { "N_trfHBARandom9-med15-0.25", 0, "HBA-Rd9-med15-0.25", "HBA-Rd9-med15-0.25", 0, 20, PURPLE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_trfHBARandom9-med15-0.25", 0, ".HBA-Rd9-med15-0.25", "HBA-Rd9-med15-0.25",0, 20, PURPLE, FALSE,  FALSE, 9.5, 11.0, 0},

    /*MYC 23bMS*/
    { "N_trfRandom3-med15-0.25", 2, "MYC-Rd3-med15-0.25", "MYC-Rd3-med15-0.25", 0, 20, LIGHTRED, FALSE, 2, 9.5, 11.0, 0} ,
    { ".N_trfRandom3-med15-0.25", 2, ".MYC-Rd3-med15-0.25", "MYC-Rd3-med15-0.25",0, 20, LIGHTRED, FALSE, FALSE, 9.5, 11.0, 0} ,

    { "N_trfRandom6-med15-0.25", 1, "MYC-Rd6-med15-0.25", "MYC-Rd6-med15-0.25", 0, 20, CYAN, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_trfRandom6-med15-0.25", 1, ".MYC-Rd6-med15-0.25", "MYC-Rd6-med15-0.25",0, 20, CYAN, FALSE,  FALSE, 9.5, 11.0, 0},

    { "N_trfRandom8-med15-0.25", 0, "MYC-Rd8-med15-0.25", "MYC-Rd8-med15-0.25", 0, 20, PURPLE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_trfRandom8-med15-0.25", 0, ".MYC-Rd8-med15-0.25", "MYC-Rd8-med15-0.25",0, 20, PURPLE, FALSE,  FALSE, 9.5, 11.0, 0},



    { "N_transfrag144047-gauss-250", 0, "Transf-gauss-250", "Transf-gauss-250", 0, 20, CYAN, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag144047-gauss-250", 0, ".Transf-gauss-250", "Transf-gauss-250",0, 20, CYAN, FALSE,  FALSE, 9.5, 11.0, 0},

    { "N_transfrag144047-gauss-500", 1, "Transf-gauss-500", "Transf-gauss-500", 0, 20, PURPLE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag144047-gauss-500", 1, ".Transf-gauss-500", "Transf-gauss-500",0, 20, PURPLE, FALSE,  FALSE, 9.5, 11.0, 0},

    { "N_transfrag144047-gauss-750", -1, "Transf-gauss-750", "Transf-gauss-750", 0, 20, RED, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".N_transfrag144047-gauss-500", -1, ".Transf-gauss-750", "Transf-gauss-750",0, 20, RED, FALSE,  FALSE, 9.5, 11.0, 0},


    { "144047-gauss-PERI-250", 0, "PERIgauss250", "PERIgauss250", 0, 20,GREEN, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".144047-gauss-PERI-250", 0, ".PERIgauss250", "PERIgauss250",0, 20, GREEN, FALSE,  FALSE, 9.5, 11.0, 0},

    { "144047-gauss-PERI-500", 0, "PERIgauss500", "PERIgauss500", 0, 20, BLUE, FALSE, 2 , 9.5, 11.0, 0} ,
    { ".144047-gauss-PERI-500", 0, ".PERIgauss500", "PERIgauss500",0, 20, BLUE, FALSE,  FALSE, 9.5, 11.0, 0},

#endif
    /* End of Nath addition of random buttons*/
#endif

    { 0, 0, 0, 0, 0, 10.0, 10.0, BLACK, FALSE, FALSE, 0.0, 16.0, 0} ,
    { 0, 0, 0, 0, 0, 10.0, 10.0, BLACK, FALSE, FALSE, 0.0, 16.0, 0}
  } ;

typedef struct pairGeneNormStruct { char *nam ; float smin, s20, s40, s60, s80, s99, smax ; } PGN ;
static PGN pgnAll[] =
{
    {"pA=12MS", 8.9428, 9.61443, 9.84471, 10.3435, 11.1876, 15.6626, 17.9373}
  , {"G1=12MS", 8.47816, 9.43077, 9.75773, 10.0785, 10.5233, 12.2351, 15.4069}
  , {"pA=22MS", 9.08917, 9.87545, 10.0448, 10.2638, 10.7597, 15.7578, 16.7902}
  , {"G1=22MS", 8.3965, 9.3723, 9.73466, 10.1249, 10.7085, 13.772, 15.878}
  , {"pA=22bMS", 8.95729, 9.74147, 9.90202, 10.1027, 10.5511, 15.1734, 16.8109}
  , {"G1=22bMS", 8.24157, 9.32356, 9.69436, 10.0989, 10.7035, 13.2992, 14.9959}
  , {"pA=11MS", 9.02781, 9.81532, 9.97421, 10.2311, 10.978, 16.8277, 18.423}
  , {"total=11MS", 9.02841, 9.83986, 10.0172, 10.2806, 10.8999, 13.9388, 15.8427}
  , {"pA=21MS", 9.1086, 9.79362, 9.9734, 10.2287, 10.8783, 15.7825, 17.0994}
  , {"total=21MS", 9.01012, 9.7585, 9.93111, 10.1655, 10.6932, 14.0589, 15.826}
  , {"total=14MS", 9.21709, 9.92099, 10.1304, 10.4114, 10.8757, 14.0896, 16.36}
  , {"G1=14MS", 8.24347, 9.4997, 9.8292, 10.1317, 10.4902, 12.2523, 15.3227}
  , {"total=24MS", 9.12042, 9.83371, 9.99667, 10.194, 10.6154, 13.9851, 15.567}
  , {"G1=24MS", 8.32183, 9.34317, 9.7191, 10.1298, 10.7481, 13.927, 15.4982}
  , {"total=24bMS", 9.11257, 9.77779, 9.94208, 10.1504, 10.5555, 13.4937, 16.1784}
  , {"G1=24bMS", 7.45258, 9.40413, 9.74699, 10.0983, 10.5885, 12.6388, 15.4645}
  , {"Nuc=13MS", 8.86087, 9.71378, 9.90692, 10.1858, 10.7003, 13.0596, 16.0993}
  , {"G1=13MS", 8.49536, 9.39762, 9.72354, 10.0246, 10.4383, 11.8329, 15.0949}
  , {"Nuc=23MS", 8.70863, 9.76692, 9.94111, 10.1624, 10.6399, 13.6515, 15.0136}
  , {"G1=23MS", 8.37027, 9.3129, 9.68771, 10.1077, 10.7986, 14.1098, 15.4327}
  , {"Nuc=23bMS", 9.01749, 9.72736, 9.88899, 10.0852, 10.4801, 13.2425, 16.0267}
  , {"G1=23bMS", 7.44448, 9.39195, 9.7404, 10.1004, 10.6061, 12.6296, 15.4303}
  , {"S=15MSr", 8.6015, 9.60125, 9.92333, 10.2245, 10.6041, 12.0904, 15.2639}
  , {"G1=15MSv", 8.44467, 9.46359, 9.78772, 10.0939, 10.4552, 11.9684, 15.1152}
  , {"S=25MS", 8.49492, 9.471, 9.86493, 10.2938, 10.9281, 13.7039, 15.7425}
  , {"G1=25MS", 8.51387, 9.37239, 9.72763, 10.1043, 10.6511, 13.1357, 15.4089}
  , {"S=28MS", 8.53004, 9.60393, 9.94953, 10.2956, 10.7506, 12.5765, 15.1531}
  , {"G1=28MS", 7.79605, 9.4609, 9.78947, 10.1113, 10.5321, 12.3581, 14.9224}
  , {"S=29MS", 8.19383, 9.16054, 9.69073, 10.3035, 11.1582, 13.6824, 15.297}
  , {"G1=29MS", 8.38051, 9.3692, 9.73921, 10.1208, 10.6257, 12.3953, 14.1897}
    , {"S=30MS", 9.06142, 9.96193, 10.3244, 10.7006, 11.2078 , 13.1766, 14.9723}
  , {"G1=30MS", 7.87415, 9.3747, 9.72567, 10.0821, 10.5678, 12.5692, 13.9867}
  , {"Samplified=26MS", 8.48593, 9.39823, 9.8257, 10.3097, 11.0375, 13.8658, 15.5537}
  , {"G1amplified=26MS", 8.46331, 9.26347, 9.67312, 10.1184, 10.7857, 13.5254, 15.234}
  , {"G1=16MSr", 8.36985, 9.56746, 9.8822, 10.1741, 10.5504, 12.1537, 15.4787}
  , {"G1=16MSv", 8.35317, 9.52767, 9.84683, 10.1481, 10.522, 12.0983, 15.4779}
  , {"G1=16bMSr", 8.50752, 9.43871, 9.77251, 10.0922, 10.5137, 11.9889, 14.4844}
  , {"G1=16bMSv", 8.63846, 9.45524, 9.77061, 10.0804, 10.4909, 11.9389, 14.4551}
  , {"G1=27MSr", 8.40982, 9.3426, 9.71899, 10.13, 10.7483, 13.633, 15.4625}
  , {"G1=27MSv", 8.52526, 9.40587, 9.76518, 10.1513, 10.7339, 13.5318, 15.399}
  , {"G1=50Crd", 8.05624, 9.29264, 9.67249, 10.0605, 10.5778, 12.6305, 14.016}
  , {"G1=50MStrvd", 7.73426, 8.75218, 9.23843, 9.74929, 10.4038, 12.6595, 13.9288}
  , {"pA=33ESr", 9.14524, 9.77445, 9.94144, 10.15, 10.6484, 14.7972, 16.6115}
  , {"G1=33ESv", 8.23948, 9.3728, 9.73219, 10.1141, 10.6688, 13.1539, 14.9964}
  , {"pA=42ES", 9.21399, 9.88967, 9.99022, 10.1013, 10.3466, 15.243, 16.6956}
  , {"G1=42ES", 8.52762, 9.36954, 9.72553, 10.1395, 10.6641, 12.7501, 14.7116}
  , {"total=35ESr", 8.10219, 9.28921, 9.70066, 10.1581, 10.8142, 13.3348, 15.5777}
  , {"G1=35ESv", 8.29732, 9.3532, 9.71418, 10.0974, 10.6452, 13.032, 15.1355}
  , {"total=44ES", 9.26959, 9.81965, 9.93456, 10.0746, 10.4234, 14.7852, 16.2556}
  , {"G1=44ES", 8.50253, 9.35418, 9.73245, 10.1229, 10.6982, 13.1938, 15.0933}
  , {"Nuc=34ESr", 8.90799, 9.30627, 9.53016, 9.8034, 10.2123, 11.6793, 13.9874}
  , {"G1=34ESv", 7.76537, 9.42475, 9.75722, 10.0904, 10.5497, 12.5444, 14.7807}
  , {"Nuc=43ES", 9.73396, 10.7531, 10.9325, 11.1504, 11.6625, 15.0696, 16.0056}
  , {"G1=43ES", 8.27519, 9.23278, 9.69058, 10.1526, 10.7744, 13.2924, 14.2488}
  , {"S=31ES", 8.5143, 9.52508, 9.88057, 10.2619, 10.7727, 12.5607, 15.0087}
  , {"G1=31ES", 8.16115, 9.43728, 9.77207, 10.1171, 10.573, 12.2897, 14.817}
  , {"S=41ES", 7.5434, 9.48683, 9.84811, 10.2323, 10.8123, 12.7647, 14.6755}
  , {"G1=41ES", 8.4088, 9.33574, 9.71777, 10.1194, 10.6585, 12.587, 14.3692}
  , {"G1=32ESr", 8.32741, 9.37194, 9.73811, 10.1154, 10.6313, 12.5251, 14.7468}
  , {"G1=32ESv", 7.84325, 9.40971, 9.753, 10.102, 10.5886, 12.6052, 14.6817}
  , {"pA=49Red", 9.24451, 9.88756, 10.0126, 10.1374, 10.3453, 13.8464, 16.0067}
  , {"G1=49Red", 8.58171, 9.38877, 9.74028, 10.1034, 10.6031, 12.8989, 15.0063}
  , {"total=48Red", 9.73376, 10.6088, 10.7643, 10.9141, 11.1419, 14.279, 16.0067}
  , {"G1=48Red", 8.43639, 9.26898, 9.69341, 10.1249, 10.6821, 12.7939, 15.0807}
  , {"Nuc=47Red", 9.33564, 9.95193, 10.1066, 10.2992, 10.7236, 13.7576, 16.0067}
  , {"G1=47Red", 8.46156, 9.37373, 9.7521, 10.1366, 10.6488, 12.891, 14.8651}
  , {"S=36Red", 8.58155, 9.66481, 10.0233, 10.4026, 10.9244, 12.8651, 15.3058}
  , {"G1=36Red", 8.11113, 9.39392, 9.74237, 10.0969, 10.5628, 12.1011, 14.2946}
  , {"S=37Red", 8.43221, 9.59154, 9.96837, 10.3597, 10.8803, 12.7616, 15.209}
  , {"G1=37Red", 7.70662, 9.39466, 9.75064, 10.1116, 10.5862, 12.2205, 14.4606}
  , {"S1=45Red", 8.65489, 9.46037, 9.86446, 10.294, 10.8714, 12.9685, 15.2975}
  , {"G1=45Red", 8.28656, 9.34856, 9.7098, 10.0794, 10.5787, 12.6031, 14.6931}
  , {"S2=38Red", 8.31108, 9.56756, 9.93692, 10.3283, 10.8836, 13.1321, 14.9478}
  , {"G1=38Red", 8.48108, 9.32096, 9.70771, 10.1148, 10.662, 12.7809, 14.64}
  , {"S3=46Red", 8.72065, 9.54531, 9.93768, 10.3413, 10.8783, 13.0062, 15.0838}
  , {"G1=46Red", 8.12426, 9.38445, 9.73337, 10.0785, 10.5442, 12.5156, 14.2598}
  , {"S123=39Red", 7.97249, 9.59362, 9.89295, 10.2156, 10.6998, 12.738, 14.8698}
  , {"G1=39Red", 8.04, 9.20123, 9.62254, 10.127, 10.7942, 12.9439, 15.1994}
  , {"S=51FLiv", 8.60506, 9.4949, 9.90029, 10.3333, 10.9217, 13.1542, 15.3912}
  , {"G1=51FLiv", 7.87421, 9.34634, 9.71399, 10.0906, 10.6007, 12.5961, 14.388}
  , {"S=52BM", 8.62545, 9.47685, 9.89996, 10.3467, 10.942, 13.1775, 15.2683}
  , {"G1=52BM", 8.07459, 9.35029, 9.71314, 10.0805, 10.5741, 12.4385, 14.8633}
  , {0,0,0,0,0,0,0,0}
  , {0,0,0,0,0,0,0,0}
} ;

static PGN pgnExonic[] =
{
    {"pA=12MS", 9.03105, 9.88095, 10.5269, 11.3625, 12.814, 17.4304, 17.9373}
  , {"pA=22MS", 9.32268, 10.1018, 10.4492, 11.1476, 12.988, 16.687, 16.7902}
  , {"pA=22bMS", 9.16222, 9.89232, 10.2245, 10.9184, 12.6589, 16.6686, 16.8109}
  , {"pA=11MS", 9.17789, 10.0295, 10.5425, 11.5614, 13.8524, 18.2128, 18.423}
  , {"total=11MS", 9.16158, 10.019, 10.3714, 10.9177, 11.8228, 15.124, 15.8427}
  , {"pA=21MS", 9.11663, 10.0075, 10.4208, 11.2302, 13.0766, 16.9277, 17.0994}
  , {"total=21MS", 9.2308, 9.87609, 10.1276, 10.524, 11.4261, 14.624, 15.826}
  , {"total=14MS", 9.26289, 10.1488, 10.4927, 10.9218, 11.7165, 15.7267, 16.36}
  , {"total=24MS", 9.21135, 9.93534, 10.1588, 10.4931, 11.3654, 15.0946, 15.567}
  , {"total=24bMS", 9.19437, 9.82734, 10.0623, 10.4461, 11.3068, 15.4234, 16.1784}
  , {"Nuc=13MS", 9.16028, 9.88408, 10.198, 10.6293, 11.2715, 14.5318, 16.0993}
  , {"Nuc=23MS", 9.22965, 9.87903, 10.1094, 10.4553, 11.1831, 14.3785, 15.0136}
  , {"Nuc=23bMS", 9.10924, 9.76465, 9.97629, 10.3142, 11.0742, 14.5183, 16.0267}
  , {"pA=33ESr", 9.14524, 9.92221, 10.2656, 10.9931, 12.5357, 16.3113, 16.6115}
  , {"total=35ESr", 8.16706, 9.34395, 9.82988, 10.3902, 11.1991, 14.2637, 15.5777}
  , {"Nuc=34ESr", 8.91643, 9.24768, 9.45112, 9.72611, 10.176, 12.0177, 13.9874}
  , {"pA=49Red", 9.30478, 9.93189, 10.1116, 10.4067, 11.4305, 15.7852, 16.0067}
  , {"total=48Red", 9.73376, 10.6418, 10.8456, 11.127, 12.0007, 15.7435, 16.0067}
  , {"Nuc=47Red", 9.36122, 9.98473, 10.2142, 10.6207, 11.5823, 15.2158, 16.0067}
, {0,0,0,0,0,0,0,0}
  , {0,0,0,0,0,0,0,0}
} ;

static PGN pgnRatio[] =
{
    {"pA=12MS", 8.9428, 9.61443, 9.84471, 10.3435, 11.1876, 15.6626, 17.9373}
  , {"pA=22MS", 9.08917, 9.87545, 10.0448, 10.2638, 10.7597, 15.7578, 16.7902}
  , {"pA=22bMS", 8.95729, 9.74147, 9.90202, 10.1027, 10.5511, 15.1734, 16.8109}
  , {"pA=11MS", 9.02781, 9.81532, 9.97421, 10.2311, 10.978, 16.8277, 18.423}
  , {"pA=21MS", 9.1086, 9.79362, 9.9734, 10.2287, 10.8783, 15.7825, 17.0994}
  , {"total=14MS", 9.21709, 9.84048, 10.0247, 10.2797, 10.7138, 13.7601, 16.332}
  , {"total=24MS", 9.12042, 9.83371, 9.99667, 10.194, 10.6154, 13.9851, 15.567}
  , {"total=24bMS", 9.11257, 9.77779, 9.94208, 10.1504, 10.5555, 13.4937, 16.1784}
  , {"Nuc=13MS", 8.86087, 9.71378, 9.90692, 10.1858, 10.7003, 13.0596, 16.0993}
  , {"Nuc=23MS", 8.70863, 9.76692, 9.94111, 10.1624, 10.6399, 13.6515, 15.0136}
  , {"Nuc=23bMS", 9.01749, 9.72736, 9.88899, 10.0852, 10.4801, 13.2425, 16.0267}
  , {"S=15MSr", 8.6015, 9.60125, 9.92333, 10.2245, 10.6041, 12.0904, 15.2639}
  , {"S=25MS", 8.49492, 9.471, 9.86493, 10.2938, 10.9281, 13.7039, 15.7425}
  , {"S=28MS", 8.53004, 9.60393, 9.94953, 10.2956, 10.7506, 12.5765, 15.1531}
  , {"S=29MS", 8.19383, 9.16054, 9.69073, 10.3035, 11.1582, 13.6824, 15.297}
  , {"S=30MS", 9.06142, 9.96193, 10.3244, 10.7006, 11.2078, 13.1766, 14.9723}
  , {"Samplified=26MS", 8.48593, 9.39823, 9.8257, 10.3097, 11.0375, 13.8658, 15.5537}
  , {"G1=16MSr", 8.36985, 9.56746, 9.8822, 10.1741, 10.5504, 12.1537, 15.4787}
  , {"G1=16bMSr", 8.50752, 9.43871, 9.77251, 10.0922, 10.5137, 11.9889, 14.4844}
  , {"G1=27MSr", 8.40982, 9.3426, 9.71899, 10.13, 10.7483, 13.633, 15.4625}
  , {"G1=50Crd", 8.05624, 9.29264, 9.67249, 10.0605, 10.5778, 12.6305, 14.016}
  , {"pA=33ESr", 9.14524, 9.77445, 9.94144, 10.15, 10.6484, 14.7972, 16.6115}
  , {"pA=42ES", 9.21399, 9.88967, 9.99022, 10.1013, 10.3466, 15.243, 16.6956}
  , {"total=35ESr", 8.10219, 9.28921, 9.70066, 10.1581, 10.8142, 13.3348, 15.5777}
  , {"total=44ES", 9.26959, 9.81965, 9.93456, 10.0746, 10.4234, 14.7852, 16.2556}
  , {"Nuc=34ESr", 8.90799, 9.30627, 9.53016, 9.8034, 10.2123, 11.6793, 13.9874}
  , {"Nuc=43ES", 9.73396, 10.7531, 10.9325, 11.1504, 11.6625, 15.0696, 16.0056}
  , {"S=31ES", 8.5143, 9.52508, 9.88057, 10.2619, 10.7727, 12.5607, 15.0087}
  , {"S=41ES", 7.5434, 9.48683, 9.84811, 10.2323, 10.8123, 12.7647, 14.6755}
  , {"G1=32ESr", 8.32741, 9.37194, 9.73811, 10.1154, 10.6313, 12.5251, 14.7468}
  , {"pA=49Red", 9.24451, 9.88756, 10.0126, 10.1374, 10.3453, 13.8464, 16.0067}
  , {"total=48Red", 9.73376, 10.6088, 10.7643, 10.9141, 11.1419, 14.279, 16.0067}
  , {"Nuc=47Red", 9.33564, 9.95193, 10.1066, 10.2992, 10.7236, 13.7576, 16.0067}
  , {"S=36Red", 8.58155, 9.66481, 10.0233, 10.4026, 10.9244, 12.8651, 15.3058}
  , {"S=37Red", 8.43221, 9.59154, 9.96837, 10.3597, 10.8803, 12.7616, 15.209}
  , {"S1=45Red", 8.65489, 9.46037, 9.86446, 10.294, 10.8714, 12.9685, 15.2975}
  , {"S2=38Red", 8.31108, 9.56756, 9.93692, 10.3283, 10.8836, 13.1321, 14.9478}
  , {"S3=46Red", 8.72065, 9.54531, 9.93768, 10.3413, 10.8783, 13.0062, 15.0838}
  , {"S123=39Red", 7.97249, 9.59362, 9.89295, 10.2156, 10.6998, 12.738, 14.8698}
  , {"S=51FLiv", 8.60506, 9.4949, 9.90029, 10.3333, 10.9217, 13.1542, 15.3912}
  , {"S=52BM", 8.62545, 9.47685, 9.89996, 10.3467, 10.942, 13.1775, 15.2683}
  , {0,0,0,0,0,0,0,0}
  , {0,0,0,0,0,0,0,0}
} ;

static PGN pgnRatioExonic[] =
{
    {"pA=12MS", 9.03105, 9.88095, 10.5269, 11.3625, 12.814, 17.4304, 17.9373}
  , {"pA=22MS", 9.32268, 10.1018, 10.4492, 11.1476, 12.988, 16.687, 16.7902}
  , {"pA=22bMS", 9.16222, 9.89232, 10.2245, 10.9184, 12.6589, 16.6686, 16.8109}
  , {"pA=11MS", 9.17789, 10.0295, 10.5425, 11.5614, 13.8524, 18.2128, 18.423}
  , {"pA=21MS", 9.11663, 10.0075, 10.4208, 11.2302, 13.0766, 16.9277, 17.0994}
  , {"total=14MS", 9.28739, 10.0466, 10.3464, 10.7408, 11.3913, 15.3574, 16.332}
  , {"total=24MS", 9.21135, 9.93534, 10.1588, 10.4931, 11.3654, 15.0946, 15.567}
  , {"total=24bMS", 9.19437, 9.82734, 10.0623, 10.4461, 11.3068, 15.4234, 16.1784}
  , {"Nuc=13MS", 9.16028, 9.88408, 10.198, 10.6293, 11.2715, 14.5318, 16.0993}
  , {"Nuc=23MS", 9.22965, 9.87903, 10.1094, 10.4553, 11.1831, 14.3785, 15.0136}
  , {"Nuc=23bMS", 9.10924, 9.76465, 9.97629, 10.3142, 11.0742, 14.5183, 16.0267}
  , {"S=15MSr", 8.62693, 9.53527, 9.86288, 10.1673, 10.5447, 11.9742, 15.2639}
  , {"S=25MS", 8.51038, 9.38626, 9.8012, 10.2664, 10.9508, 13.7077, 15.7425}
  , {"S=28MS", 8.65983, 9.5316, 9.87681, 10.2401, 10.7157, 12.4607, 15.1531}
  , {"S=29MS", 8.21372, 8.91493, 9.31166, 9.8316, 10.6403, 13.408, 15.297}
  , {"S=30MS", 9.082, 9.92266, 10.2942, 10.6943, 11.2274, 13.1309, 14.9723}
  , {"Samplified=26MS", 8.54511, 9.28641, 9.71706, 10.2501, 11.0362, 13.8792, 15.5537}
  , {"G1=16MSr", 8.46463, 9.52998, 9.84998, 10.1358, 10.5153, 12.0268, 15.4787}
  , {"G1=16bMSr", 8.50752, 9.34997, 9.68461, 10.0077, 10.427, 11.8638, 14.4844}
  , {"G1=27MSr", 8.46064, 9.26372, 9.66065, 10.1071, 10.7868, 13.6763, 15.4625}
  , {"G1=50Crd", 8.05624, 9.1845, 9.55261, 9.96385, 10.5104, 12.5233, 14.016}
  , {"pA=33ESr", 9.14524, 9.92221, 10.2656, 10.9931, 12.5357, 16.3113, 16.6115}
  , {"pA=42ES", 9.28287, 9.99378, 10.255, 10.9948, 12.9456, 16.6131, 16.6956}
  , {"total=35ESr", 8.16706, 9.34395, 9.82988, 10.3902, 11.1991, 14.2637, 15.5777}
  , {"total=44ES", 9.39323, 9.9177, 10.191, 10.8398, 12.4401, 16.0995, 16.2556}
  , {"Nuc=34ESr", 8.91643, 9.24768, 9.45112, 9.72611, 10.176, 12.0177, 13.9874}
  , {"Nuc=43ES", 9.89363, 10.8524, 11.1824, 11.7803, 12.9611, 15.9704, 16.0056}
  , {"S=31ES", 8.58038, 9.50992, 9.88167, 10.2884, 10.8259, 12.5911, 15.0087}
  , {"S=41ES", 7.60574, 9.46096, 9.85026, 10.2726, 10.8908, 12.8277, 14.6755}
  , {"G1=32ESr", 8.33692, 9.32555, 9.69864, 10.0917, 10.6311, 12.4883, 14.7468}
  , {"pA=49Red", 9.30478, 9.93189, 10.1116, 10.4067, 11.4305, 15.7852, 16.0067}
  , {"total=48Red", 9.73376, 10.6418, 10.8456, 11.127, 12.0007, 15.7435, 16.0067}
  , {"Nuc=47Red", 9.36122, 9.98473, 10.2142, 10.6207, 11.5823, 15.2158, 16.0067}
  , {"S=36Red", 8.64111, 9.6389, 10.024, 10.4368, 10.9977, 12.8969, 15.3058}
  , {"S=37Red", 8.58558, 9.53406, 9.9208, 10.3462, 10.899, 12.7285, 15.209}
  , {"S1=45Red", 8.67961, 9.42858, 9.86163, 10.3306, 10.9584, 12.9626, 15.2975}
  , {"S2=38Red", 8.39093, 9.53298, 9.93566, 10.3674, 10.9502, 13.0893, 14.9478}
  , {"S3=46Red", 8.74924, 9.46828, 9.86809, 10.2973, 10.8675, 12.8679, 15.0838}
  , {"S123=39Red", 8.11185, 9.53148, 9.82703, 10.1654, 10.6919, 12.6845, 14.8543}
  , {"S=51FLiv", 8.68022, 9.44878, 9.87542, 10.3447, 10.9737, 13.108, 15.3912}
  , {"S=52BM", 8.62545, 9.41136, 9.85046, 10.3315, 10.9688, 13.089, 15.2683}
  , {0,0,0,0,0,0,0,0}
  , {0,0,0,0,0,0,0,0}
} ;

#endif
#endif
