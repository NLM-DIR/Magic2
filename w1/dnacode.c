#include "acedna.h"
#include "peptide.h"

/* this is the mapping use to print out */
char dnaDecodeExtendedChar[] =
 { '-','a','t','w','g','r','k','d','c','m','y','h','s','v','b','n',
     0,  0,  0,  0, '.'} ;
char rnaDecodeExtendedChar[] =
 { '-','a','u','w','g','r','k','d','c','m','y','h','s','v','b','n',
     0,  0,  0,  0, '.'} ;
				/* RD added '.' to help alignment code */
char dnaDecodeChar[] =
 { '-','a','t','w','g','r','k','d','c','m','y','h','s','v','b','n','>','<'} ;
char rnaDecodeChar[] =
 { '-','a','u','w','g','r','k','d','c','m','y','h','s','v','b','n','>','<'} ;

/* this is the mapping use to parse in */
char dnaEncodeChar[] =
{  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   RS_,   0,   FS_,   0,

   0,  A_,  B_,  C_,  D_,   0,   0,  G_,  H_,   0,   0,  K_,   0,  M_,  N_,   0,
   0,   0,  R_,  S_,  T_,  U_,  V_,  W_,   0,  Y_,   0,   0,   0,   0,   0,   0,
   0,  A_,  B_,  C_,  D_,   0,   0,  G_,  H_,   0,   0,  K_,   0,  M_,  N_,   0,
   0,   0,  R_,  S_,  T_,  U_,  V_,  W_,  N_,  Y_,   0,   0,   0,   0,   0,   0,
} ;

char hardComplementBase[] =	
 { 0, T_,A_,W_,C_,Y_,M_,H_,G_,K_,R_,D_,S_,B_,V_,N_,RS_,FS_ } ;

char checkComplementBase (char cc)
{
  int n = cc ;
  if (n < 0 || n > 17)
    messcrash ("complementBase(%d) out of range", n) ;
  return hardComplementBase[n] ;
}

/********************************************************/
/*
The naming of proteins:

There are 4 names:
- a number: 0, 1, 2, etc.
- a single character: A, C, D, etc
- an abbreviated name, Ala, Cys, etc
- a full name, Alanine, Cysteine, etc

All but the number are standardized, though I have not actually seen
a standards document.  (ms - 7/2002)


The protein names:

0  A  Ala  Alanine
1  C  Cys  Cysteine
2  D  Asp  Aspartate (also Aspartic Acid)
3  E  Glu  Glutamate (also Glutamic Acid)
4  F  Phe  Phenylalanine
5  G  Gly  Glycine
6  H  His  Histidine
7  I  Ile  Isoleucine
8  K  Lys  Lysine
9  L  Leu  Leucine
10 M  Met  Methionine
11 N  Asn  Asparagine
12 P  Pro  Proline
13 Q  Gln  Glutamine
14 R  Arg  Arginine
15 S  Ser  Serine
16 T  Thr  Threonine
17 V  Val  Valine
18 W  Trp  Tryptophan
19 Y  Tyr  Tyrosine
20 X  Xxx  Unknown / Undetermined
21 .  ???  ??? ( nobody knows what this is doing here; only acedb has it. )
22 *  Ter  Stop / Terminate transcription
23 Z  Glx  Glu/Gln
24 B  Asx  Asp/Asn
25 U  Sec  Selenocysteine (present in human, encoded by UGA which is normally a stop, often missing from other lists )
26 O  Pyl  Pyrrolysine  (in archea and some bacteria, encoded by UAG which is normally a stop)

  notes for maintaining these tables:

  'A' is 0x41
  '*' is 0x2a
  '.' is 0x2e
  '0' is 0x30

  don't know why digits and backspace are -4 in pepEncodeChar

pepDecodeChar[] turns a protein number into the printable character.
pepEncodeChar[] turns a protein character into the correct number.
pepName[] turns a protein character into the long name.
pepShortName[] turns a protein character into an abbreviated name.

*/

char pepDecodeChar[] =

/* 0    1    2    3    4    5    6    7    8    9   10   11   12 */

{ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
  'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', '.', '*', 'Z', 'B', 'U', 'O', 0} ;

/* 13   14   15   16   17   18   19   20   21   22  23   24   25  26 */

signed char pepEncodeChar[] =
{ 
/* 00 */ -1, -1, -1, -1,  -1, -1, -1, -1,  -4, -1, -1, -1,  -1, -1, -1, -1,
/* 10 */ -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
/* 20 */ -4, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, 22, -1,  -1, -1, 21, -1,  /* '*' */
/* 30 */ -4, -4, -4, -4,  -4, -4, -4, -4,  -4, -4, -1, -1,  -1, -1, -1, -1,  /* digits */

/* 40 */ -1,  0, 24,  1,   2,  3,  4,  5,   6,  7, -1,  8,   9, 10, 11, 26,
/* 50 */ 12, 13, 14, 15,  16, 25, 17, 18,  20, 19, 23, -1,  -1, -1, -1, -1,
/* 60 */ -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,
/* 70 */ -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1

/* 80 */ -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1
/* 90 */ -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1
/* a0 */ -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1
/* b0 */ -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1

/* c0 */ -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1
/* d0 */ -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1
/* e0 */ -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1
/* f0 */ -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1,  -1, -1, -1, -1
} ;

char *pepName[] =
{  
/* 00 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* 10 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* 20 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  "Stop",
						       0,   0, 0,  "???",  0,
/* 30 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 

/* 40 */ 0, 		"Alanine", 	 "Asp/Asn", 	"Cysteine", 
/* 44 */"Aspartate", 	"Glutamate", 	 "Phenylalanine","Glycine", 
/* 48 */"Histidine", 	"Isoleucine", 	 0, 		"Lysine", 
/* 4c */"Leucine", 	"Methionine", 	 "Asparagine", 	"Pyrrolysine",

/* 50 */"Proline", 	"Glutamine", 	 "Arginine", 	"Serine", 
/* 54 */"Threonine", 	"Selenocysteine" ,"Valine", 	"Tryptophan", 
/* 58 */"Unknown", 	"Tyrosine", 	 "Glu/Gln",  	0,  
/* 5a */ 0,  		0,  		0,  		0, 

/* 60 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* 70 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 

/* 80 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* 90 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* a0 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* b0 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 

/* c0 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* d0 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* e0 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* f0 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
} ;		

char *pepShortName[] =
{  
/* 00 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* 10 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* 20 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  "stp",  0,   0, 0,  "???",  0,
/* 30 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 

/* 40 */      0, "Ala", "AsX", "Cys",    "Asp", "Glu", "Phe", "Gly", 
/* 48 */  "His", "Ile",     0, "Lys",    "Leu", "Met", "Asn", "Pyl",
/* 50 */  "Pro", "Gln", "Arg", "Ser",    "Thr", "Sec", "Val", "Trp", 
/* 58 */  "Xxx", "Tyr", "Glx",     0,        0,     0,     0,     0, 
 
/* 60 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* 70 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 

/* 80 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* 90 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* a0 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* b0 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* c0 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* d0 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* e0 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0, 
/* f0 */ 0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0
} ;		

/***************************************************************/
/* mhmp 23.09.98 */  
/* Antobodies A Laboratory Manual  Ed Harlow David Lane p. 661*/
/* Cold Spring Harbor Laboratory 1988 page 661 */
/* je mets le Stop a -1 et les Asp/Asn,  Glu/Gln  et Unknown a 0 */
int  molecularWeight[] =
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 

   0, 89, 133, 121, 133, 147, 165, 75, 155, 131, 0,146, 131, 149, 132, 255,
   115, 146, 174, 105, 119, 168, 117, 204, 120, 181, 147,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 

   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
} ;	
/**************************************************************/
/**************************************************************
* Protein Molecular Weight Section, after Jim Ostell *
*
*  Returns a protein molecular weight of an encoded peptide array
*    If it cannot calculate the value it returns -1.0
*
***************************************************************/


/* Values are in  pepDecodeChar order:
   B is really D or N, but they are close so is treated as D
   Z is really E or Q, but they are close so is treated as E
   X is hard to guess, so the calculation fails on X

  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y  X  .  *  Z  B  U O
*/

static char C_atoms[27] =
  { 3, 3, 4, 5, 9, 2, 6, 6, 6, 6, 5, 4, 5, 5, 6, 3, 4, 5,11, 9, 0, 0, 0, 5, 4, 3, 12};
static char H_atoms[27] =
  { 5, 5, 5, 7, 9, 3, 7,11,12,11, 9, 6, 7, 8,12, 5, 7, 9,10, 9, 0, 0, 0, 7, 5, 5, 21};
static char N_atoms[27] =
  { 1, 1, 1, 1, 1, 1, 3, 1, 2, 1, 1, 2, 1, 2, 4, 1, 1, 1, 2, 1, 0, 0, 0, 1, 1, 1, 3};
static char O_atoms[27] =
  { 1, 1, 3, 3, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 2, 1, 1, 2, 0, 0, 0, 3, 3, 1, 3};
static char S_atoms[27] =
  { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
static char Se_atoms[27] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0};

/**************************************************************
*
*  Returns a protein molecular weight of an encoded peptide array
*    If it cannot calculate the value it returns -1.0
*
***************************************************************/

int pepWeight (Array pep)  /* encoded array */
{
  char *cp ;
  int residue, w = -1 ;

  int Ccnt = 0;  /* initialize counters */
  int Hcnt = 2;  /* always start with water */
  int Ocnt = 1;  /* H20 */
  int Ncnt = 0;
  int Scnt = 0;
  int Secnt = 0;
  
  
  if (pep)
    {
      int i = arrayMax(pep) ;
      cp = arrp (pep, 0, char) - 1 ;
      while (cp++, i--)
	{
	  residue = *cp ;
	  if (residue == 22)
	    break ;
	  if (H_atoms[residue] == 0)  /* unsupported AA */
	    return w ;                /* bail out */
	  Ccnt += C_atoms[residue];
	  Hcnt += H_atoms[residue];
	  Ncnt += N_atoms[residue];
	  Ocnt += O_atoms[residue];
	  Scnt += S_atoms[residue];
	  Secnt += Se_atoms[residue];
	}
    }
  
  w = (12.01115 * Ccnt) + (1.0079 * Hcnt) +
    (14.0067 * Ncnt) + (15.9994 * Ocnt) +
    (32.064 * Scnt) + (78.96 * Secnt);
  
  return w ;
}

/****************************************/
/****************************************/
/* pI according to Expasy */
/* 
 * Data an algorithm received from EXPASY, sept 28, 2006 
 *    after i sent a request to the helpdesk at 
 *    http://ca.expasy.org/cgi-bin/pi_tool?protein=MPQR 
 * ----------------------------------------------------------------- 
 * Elisabeth Gasteiger 
 * Swiss Institute of Bioinformatics 
 * CMU - 1, rue Michel Servet                Tel. (+41 22) 379 58 79 
 * CH - 1211 Geneva 4 Switzerland            Fax  (+41 22) 379 58 58 
 * Elisabeth.Gasteiger@isb-sib.ch            http://www.expasy.org/ 
 * -----------------------------------------------------------------  
 *    
 * I hope this helps.
 * 
 * We would like to make sure that the scientific community
 * is aware of the services we offer on ExPASy. This is why
 * your email address has been added to the list of subscribers
 * to our Swiss-Flash service of electronic bulletins
 * (http://www.expasy.org/swiss-flash/).
 * If you do not wish to receive these bulletins by email, you
 * can go to the above URL address and unsubscribe.
 * 
 * Best regards
 * Elisabeth Gasteiger 
 *
 *
 * Table of pk values :
 *  Note: the current algorithm does not use the last two columns.
 * Each row corresponds to an amino acid starting with Ala. J, O and U are
 * inexistant, but here only in order to have the complete alphabet. 
 *
 * 
 *
 */

#define PH_MIN 0 /* minimum pH value */
#define PH_MAX 14 /* maximum pH value */
#define MAXLOOP 2000 /* maximum number of iterations */
#define EPSI 0.0001 /* desired precision */

float pepPI (Array pep) /* encoded array */
{
  int i, ii, n, max, aaCount[256] ;
  int cTermResidue, nTermResidue = 0 ;
  char *cp ;

  double cter, nter, carg, chis, clys, casp, cglu, ccys, ctyr, charge ;
  double phMin = PH_MIN ;
  double phMax = PH_MAX ;
  double phMid ;

  double pepPk[27][5] = {
    /* Ct   Nt  Sm     Sc     Sn  */
    { 3.55, 7.59, 0.   , 0.   , 0.    }, /* A*/
    { 3.55, 7.50, 0.   , 0.   , 0.    }, /* B*/
    { 3.55, 7.50, 9.00 , 9.00 , 9.00  }, /* C*/
    { 4.55, 7.50, 4.05 , 4.05 , 4.05  }, /* D*/
    { 4.75, 7.70, 4.45 , 4.45 , 4.45  }, /* E*/
    { 3.55, 7.50, 0.   , 0.   , 0.    }, /* F*/
    { 3.55, 7.50, 0.   , 0.   , 0.    }, /* G*/
    { 3.55, 7.50, 5.98 , 5.98 , 5.98  }, /* H*/
    { 3.55, 7.50, 0.   , 0.   , 0.    }, /* I*/
    { 0.00, 0.00, 0.   , 0.   , 0.    }, /* J*/
    { 3.55, 7.50, 10.00, 10.00, 10.00 }, /* K*/
    { 3.55, 7.50, 0.   , 0.   , 0.    }, /* L*/
    { 3.55, 7.00, 0.   , 0.   , 0.    }, /* M*/
    { 3.55, 7.50, 0.   , 0.   , 0.    }, /* N*/
    { 0.00, 0.00, 0.   , 0.   , 0.    }, /* O*/
    { 3.55, 8.36, 0.   , 0.   , 0.    }, /* P*/
    { 3.55, 7.50, 0.   , 0.   , 0.    }, /* Q*/
    { 3.55, 7.50, 12.0 , 12.0 , 12.0  }, /* R*/
    { 3.55, 6.93, 0.   , 0.   , 0.    }, /* S*/
    { 3.55, 6.82, 0.   , 0.   , 0.    }, /* T*/
    { 0.00, 0.00, 0.   , 0.   , 0.    }, /* U*/
    { 3.55, 7.44, 0.   , 0.   , 0.    }, /* V*/
    { 3.55, 7.50, 0.   , 0.   , 0.    }, /* W*/
    { 3.55, 7.50, 0.   , 0.   , 0.    }, /* X*/
    { 3.55, 7.50, 10.00, 10.00, 10.00 }, /* Y*/
    { 3.55, 7.50, 0.   , 0.   , 0.    }, /* Z*/
    { 3.55, 7.50, 0.   , 0.   , 0.    } /* O, i copy the values for Lysine */
  } ;
  
  if (!pep || ! (max = arrayMax(pep)))
    goto done ;
  
  /* Compute the amino-acid composition */
  /* Look up N-terminal and C-terminal residue. */
  
  memset (aaCount, 0, sizeof(aaCount)) ;
  for (i = n = 0, cp = arrp(pep,0,char) ; i < max ; i++, cp++)
    {
      n = pepDecodeChar[(int)*cp] - 'A' ;
      if (n < 0)  /* stop = '*'- 'A' = -23 */
	{
	  max = i ;
	  break ;
	}
      aaCount [n]++ ;
      if (!i)  nTermResidue = n ;
      cTermResidue = n ;
    }
  
  for (ii = 0, charge = 1.0 ; 
       ii < MAXLOOP && (phMax - phMin) > EPSI ;
       ii++)
    {
      phMid = phMin + (phMax - phMin) / 2;
      
      cter = exp (-pepPk[cTermResidue][0]) /
	(exp (-pepPk[cTermResidue][0]) + exp (-phMid));
      nter = exp (-phMid) /
	(exp (-pepPk[nTermResidue][1]) + exp (-phMid));
      
      carg = aaCount['R' - 'A'] * exp (-phMid) /
	(exp (-pepPk['R' - 'A'][2]) + exp (-phMid));
      chis = aaCount['H' - 'A'] * exp (-phMid) /
	(exp (-pepPk['H' - 'A'][2]) + exp (-phMid));
      clys = aaCount['K' - 'A'] * exp (-phMid) /
	(exp (-pepPk['K' - 'A'][2]) + exp (-phMid));
      
      casp = aaCount['D' - 'A'] * exp (-pepPk['D' - 'A'][2]) /
	(exp (-pepPk['D' - 'A'][2]) + exp (-phMid));
      cglu = aaCount['E' - 'A'] * exp (-pepPk['E' - 'A'][2]) /
	(exp (-pepPk['E' - 'A'][2]) + exp (-phMid));
      
      ccys = aaCount['C' - 'A'] * exp (-pepPk['C' - 'A'][2]) /
	(exp (-pepPk['C' - 'A'][2]) + exp (-phMid));
      ctyr = aaCount['Y' - 'A'] * exp (-pepPk['Y' - 'A'][2]) /
	(exp (-pepPk['Y' - 'A'][2]) + exp (-phMid)) ;
      
      charge = carg + clys + chis + nter -
	(casp + cglu + ctyr + ccys + cter) ;
      
      if (charge > 0.0)
	phMin = phMid ;
      else
	phMax = phMid ;
    }
  done :
    n = 100 *  (phMax + phMin) / 2.0 ;
    return  (float)(n/100.0) ;
} /* pepPiExpasy */

/**************************************************************/
/**************************************************************/
/*** next two functions were in dnacpt.c but should be here ***/

/*
* codon() is obsolete.  Use pepGetTranslationTable() and e_code()
* for new code.
*/

char codonMito (const char* cp)
{
  switch(*cp)
    { 
    case T_: switch(*(cp+1))
      {
      case T_: switch(*(cp+2))
	{ 
	case T_: case C_: case Y_: return 'F';
	case A_: case G_: case R_: return 'L';
	default: return 'X' ;
	}
      case C_: return 'S' ;
      case A_: switch(*(cp+2))
	{ 
	case T_: case C_: case Y_: return 'Y';
	case A_: case G_: case R_: return '*';
	default: return 'X' ;
	}
      case G_: switch(*(cp+2))
	{ 
	case T_: case C_: case Y_: return 'C';
	case A_: return 'W';
	case G_: return 'W';
	default: return 'X' ;
	}
      default: return 'X' ;
      }
    case C_: switch(*(cp+1))
      {
      case T_: return 'L' ;
      case C_: return 'P' ;
      case A_: switch(*(cp+2))
	{
	case T_: case C_: case Y_: return 'H';
	case A_: case G_: case R_: return 'Q';
	default: return 'X' ;
	}
      case G_: return 'R' ;
      default: return 'X' ;
      }
    case A_: switch(*(cp+1))
      {
      case T_: switch(*(cp+2))
	{
	case T_: case C_: 
	case Y_: case M_: case W_:
	case H_: return 'I';
	case G_: case A_: return 'M';
	default: return 'X' ;
	}
      case C_: return 'T' ;
      case A_: switch(*(cp+2))
	{
	case T_: case C_: case Y_: return 'N';
	case A_: case G_: case R_: return 'K';
	default: return 'X' ;
	}
      case G_: switch(*(cp+2))
	{
	case T_: case C_: case Y_: return 'S';
	case A_: case G_: case R_: return '*';
	default: return 'X' ;
	}
      default: return 'X' ;
      }
    case G_: switch(*(cp+1))
      {
      case T_: return 'V' ;
      case C_: return 'A' ;
      case A_: switch(*(cp+2))
	{
	case T_: case C_: case Y_: return 'D';
	case A_: case G_: case R_: return 'E';
	default: return 'X' ;
	}
      case G_: return 'G' ;
      default: return 'X' ;
      }
  /***************** Wobble on first Letter *******************/
    case Y_:  /* T or C */ switch(*(cp+1))
      {
      case U_: switch(*(cp+2))
	{
	case A_: case G_: case R_: return 'L';
	default: return 'X' ;
	}
      default: return 'X' ;
      }
    case M_:  /* A or C */ switch(*(cp+1))
      {
      case G_: switch(*(cp+2))
	{
	case A_: case G_: case R_: return 'R';
	default: return 'X' ;
	}
      default:return 'X' ;
      }
    case R_: switch(*(cp+1))
      {
      case A_: switch(*(cp+2))
	{
	case T_: case C_: case Y_: return 'B'; /*Asp, Asn */
	default: return 'X' ;
	}
      default: return 'X' ;
      }
    case S_: switch(*(cp+1))
      {
      case A_: switch(*(cp+2))
	{
	case A_: case G_: case R_: return 'Z'; /* Glu, Gln */
	default: return 'X' ;
	}
      default: return 'X' ;
      }
    default: return 'X' ;
    }
}

/****************/

char reverseCodon (const char* cp)
{
  char temp[3] ;

  temp[0] = complementBase(cp[2]) ;
  temp[1] = complementBase(cp[1]) ;
  temp[2] = complementBase(cp[0]) ;
  return codon (temp) ;
}
/****************/

char antiCodon (const char* cp)
{
  char temp[3] ;

  temp[0] = complementBase(cp[0]) ;
  temp[1] = complementBase(cp[-1]) ;
  temp[2] = complementBase(cp[-2]) ;
  return codon (temp) ;
}

/****************/
/****************/

char codon (const char* cp)
{
  switch(*cp)
    { 
    case T_: switch(*(cp+1))
      {
      case T_: switch(*(cp+2))
	{ 
	case T_: case C_: case Y_: return 'F';
	case A_: case G_: case R_: return 'L';
	default: return 'X' ;
	}
      case C_: return 'S' ;
      case A_: switch(*(cp+2))
	{ 
	case T_: case C_: case Y_: return 'Y';
	case A_: case G_: case R_: return '*';
	default: return 'X' ;
	}
      case G_: switch(*(cp+2))
	{ 
	case T_: case C_: case Y_: return 'C';
	case A_: return '*';
	case G_: return 'W';
	default: return 'X' ;
	}
      default: return 'X' ;
      }
    case C_: switch(*(cp+1))
      {
      case T_: return 'L' ;
      case C_: return 'P' ;
      case A_: switch(*(cp+2))
	{
	case T_: case C_: case Y_: return 'H';
	case A_: case G_: case R_: return 'Q';
	default: return 'X' ;
	}
      case G_: return 'R' ;
      default: return 'X' ;
      }
    case A_: switch(*(cp+1))
      {
      case T_: switch(*(cp+2))
	{
	case T_: case C_: case A_: 
	case Y_: case M_: case W_:
	case H_: return 'I';
	case G_: return 'M';
	default: return 'X' ;
	}
      case C_: return 'T' ;
      case A_: switch(*(cp+2))
	{
	case T_: case C_: case Y_: return 'N';
	case A_: case G_: case R_: return 'K';
	default: return 'X' ;
	}
      case G_: switch(*(cp+2))
	{
	case T_: case C_: case Y_: return 'S';
	case A_: case G_: case R_: return 'R';
	default: return 'X' ;
	}
      default: return 'X' ;
      }
    case G_: switch(*(cp+1))
      {
      case T_: return 'V' ;
      case C_: return 'A' ;
      case A_: switch(*(cp+2))
	{
	case T_: case C_: case Y_: return 'D';
	case A_: case G_: case R_: return 'E';
	default: return 'X' ;
	}
      case G_: return 'G' ;
      default: return 'X' ;
      }
  /***************** Wobble on first Letter *******************/
    case Y_:  /* T or C */ switch(*(cp+1))
      {
      case U_: switch(*(cp+2))
	{
	case A_: case G_: case R_: return 'L';
	default: return 'X' ;
	}
      default: return 'X' ;
      }
    case M_:  /* A or C */ switch(*(cp+1))
      {
      case G_: switch(*(cp+2))
	{
	case A_: case G_: case R_: return 'R';
	default: return 'X' ;
	}
      default:return 'X' ;
      }
    case R_: switch(*(cp+1))
      {
      case A_: switch(*(cp+2))
	{
	case T_: case C_: case Y_: return 'B'; /*Asp, Asn */
	default: return 'X' ;
	}
      default: return 'X' ;
      }
    case S_: switch(*(cp+1))
      {
      case A_: switch(*(cp+2))
	{
	case A_: case G_: case R_: return 'Z'; /* Glu, Gln */
	default: return 'X' ;
	}
      default: return 'X' ;
      }
    default: return 'X' ;
    }
}

/****************/
/****************/
/****************/

