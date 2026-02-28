/*  File: lexsubs4.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1991
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Mar 12 17:18 2002 (rd)
 * * Mar 11 00:26 2002 (rd): changed sessionUserKey() in readTimeStamps() 
 		to thisSession.userKey because a kernel operation.
 * * Mar 10 23:15 2002 (rd): added #include <ctype.h>.
 * * Sep  2 14:47 1998 (edgrif): Removed faulty declaration of lexAlphaMark,
 *              and lexAlphaSaveAll, which are defined in lex.h anyway.
 * * Aug 29 20:07 1994 (mieg): _Tags are now defined as KEY variables
 *              their value is computed dynamically, ensuring consistency
 *              accross databases
 * * Jan 10 00:54 1992 (mieg): added lexAlias, eliminated lastmodif
 * * Dec 11 12:22 1991 (mieg): lexcleanup to remove spaces in addkey
 * * Nov  5 21:28 1991 (mieg): introduced the hash package
 * Created: Tue Nov  5 21:28:40 1991 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: lexsubs4.c,v 1.37 2020/05/30 16:50:32 mieg Exp $ */

/***************************************************************/
/***************************************************************/
/**  File lexsubs4.c :                                         **/
/**  Handles as Arrays  the lexique of the ACeDB program.     **/
/***************************************************************/
/***************************************************************/
/*                                                             */
/*  Many    routines are public :                              */
/*     avail,show, make/dump,  2define,                        */
/*     iskey,lock/unlock, iskeylocked,                         */
/*     KEY2LEX,key2word/word2key, addkey, kill,                */
/*     randomkey,lexNext and lexstrcmp.                        */
/*                                                             */
/* NB lextstrcmp & lexstrCasecmp now moved to w1/utils.c       */
/*                                                             */
/*  There are 256 vocabularies, the first 128 are generic      */
/*  and common to all implementations of the base, the others  */
/*  can be handled freely by every user.                       */
/*                                                             */
/*  Voc[0] contains the flow control dialog, the others the    */
/*  names of the genes, alleles, bits of the map and other     */
/*  objects, and the keywords describing the phenotypes.       */
/*                                                             */
/*  Each word is coded by a key. A key is 4 bytes long, the    */
/*  first byte is the vocabulary number, thus, the maximum     */
/*  number of entry per voc. is 16M. Each key may then refer   */
/*  to a disk address (or 0 if the word is not the name of an  */
/*  object or a list), and to a cache control area if the      */
/*  object is loaded.                                          */
/*                                                             */
/*  Lexshow is for debugging.                                  */
/*  Lexmake is called on entering the session, it reads in the */
/*  lexiques from disk and return 0 if ok. LexSave writes      */
/*  them back to disk. Note that the lexique is not saved on   */
/*  every modification, therefore, if a crash occurs during    */
/*  a session, the lexiques must be reconstructed from the     */
/*  disk.                                                      */
/*                                                             */
/*  iskey returns 0 on unknown keys, 1 if the key is empty     */
/*  or 2 if an object is associated to this key.               */
/*  lock prevents the simultaneous updating by 2 processes     */
/*  it is invoked by BSlock_make                               */
/*                                                             */
/*  Name make a key into a word, and word2key a word           */
/*  into a key if the word is known.  Lexaddkey(w,k,t) adds    */
/*  a new word w and sets its key k in vocabulary t.           */
/*                                                             */
/*  Lexstrcmp is a case insensitive strcmp.                    */
/*                                                             */
/***************************************************************/
#define CHRONOxxxx
#include <ctype.h>

#if defined(ACEDB4) || defined(ACEDB5)

#define DEF_BP /* avoid typedefing it as void* in disk.h */
typedef struct block *BLOCKP ;
typedef BLOCKP BP ;

#include "acedb.h"
#include "lex_bl_.h"
#include "disk.h"    /*call to diskprepare in lexinit*/
#include "call.h"
#include "cache.h"
#include "chrono.h"
#include "client.h"

#define MAXLEX   ((1<<24) - 1)   /*number of entries per vocabulary <=FFFFFF*/

                                 /*number of char in a vocabulary<max(LEXOFFSET)*/
/*#define MAXVOCAB (1<<31 - 1)*/
#define MAXVOCAB (0x7fffffff)	/* compiler probably interpreted the above as 1<< (31-1) --> 1073741824
				   wheras it looks like this number is supposed to be 2147483647,
				   but the alpha-compiler complains about 1<<31 to be too large, 
				   so we make that number directly by its HEX representation */

#define MAXTABLE 256             /*number of vocabulary tables <= 256*/
                                 /* watchout 256 meaning MAXTABLE used in coresubs */
                                 /* and also in parse.c , maybe elsewhere */
                                 /* and corresponds to the first byte inside KEY */

/***********************/

/*
* Lexi1 is an array of LEXI1 arrays.  Lexi1[class number] is the array
* for that class.  Index that array by the key number to get the struct
* for a specific object.  In the LEXI1 struct, the nameoffset field is the
* offset into the Voc string where you can find the object name.
*/
static Array  Lexi1 [MAXTABLE];

static Array  Lexi2 = 0 ;
static int minFreeLexi2 = 0 ;
static KEYSET LexHashTable[MAXTABLE] ;
static int    nClass[MAXTABLE] ;  /* number of bits for hasher */

/*
* Voc is an array of characters for each class.  It is just a bunch of
* object names concatenated together, with just enough zeroes to
* respect word alignmnet.  Lexi1 contains information you
* can use to get at an object name.  Voc[class number] is the data
* for that class.
*/
static Stack Voc  [MAXTABLE];

static BOOL  lexIsRead[MAXTABLE] ;
static char* nakedName(KEY key) ;  /* Does not follow aliases */
static int vocmodif[MAXTABLE];
static int lexhmodif[MAXTABLE];
static int leximodif[MAXTABLE];
static int lexSessionStartSize[MAXTABLE] ;
KEYSET touchedByClient = 0 ;


static Array timeStamps[MAXTABLE];
static BOOL timeModif[MAXTABLE];
static void lexKeyUpdated(KEY key);

#define si1 ((U_Int)sizeof(LEXI1))
#define si2 ((U_Int)sizeof(LEXI2))

static void lexvocread(int table) ;
static void lexvocwrite(int table);
static void lexiread(int table) ;
static void lexiwrite(int table);
static void lexhread(int table) ;
static void lexhwrite(int table);
static void lexReadTimeStamps(int t);
static void lexWriteTimeStamps(int t);

extern void lexDefineSystemTags(void) ;   /* in whooks/sysclass.c */
extern void sysClassInit (void) ;

static void lexHashInsert(int classe, KEY key) ;
void lexReHashClass(int classe) ;

static LEXP1 KEY2LEX1(KEY kk) ;
static LEXP KEY2LEX2(KEY kk) ;

#define PRINT FALSE	/* set TRUE to show save progress */

/**************************************************************/

void lexavail(unsigned long *vocnum,
	      unsigned long *lex1num,
	      unsigned long *lex2num,
	      unsigned long *vocspace, 
	      unsigned long *hspace,
	      unsigned long *totalspace)
{
  register  int t = 256;
  unsigned long n = 0, v = 0, k1 = 0, k2 = 0, h = 0, tt = 0;
  while(t--)
    if(Lexi1[t])
      { n++;
        v += (unsigned long) stackMark(Voc[t]) ;
        k1 += (unsigned long) arrayMax(Lexi1[t]) ;
        h += LexHashTable[t] ? arrayMax(LexHashTable[t]) : 0 ;
	tt += Voc[t]->a->dim +( Lexi1[t]->dim)*(Lexi1[t]->size)
	  + (LexHashTable[t] ?  (LexHashTable[t]->dim)*(LexHashTable[t]->size) : 0 ) ;
      }
  k2 = (unsigned long) arrayMax(Lexi2) ;
  tt += Lexi2->dim * Lexi2->size ;
  *vocspace=v; *vocnum=n; *lex1num=k1; *lex2num=k2; 
  *hspace = h ; *totalspace = tt/1024;
}

/*************************************************/

void lexInit(void)
{ int n, i ;
  LEXP lxp ;
  register int t = MAXTABLE ;

  minFreeLexi2 = 1 ;
  n = 2000 ;
  Lexi2 = arrayCreate (n, LEXI2) ;
  array (Lexi2, n - 1, LEXI2).dlx.lx = 0 ; /* make room */
  for (i = 1 , lxp = arrp (Lexi2, i, LEXI2) ; i < n - 1 ; i++, lxp++)
   lxp->dlx.lx = i + 1 ; 

  while(t--)
    { vocmodif[t] = lexhmodif[t] = leximodif[t] = timeModif[t] = FALSE ;
      Voc[t] = 0 ; Lexi1[t] = LexHashTable[t] = timeStamps[t] = 0 ;
      nClass[t] = 0 ;
      lexIsRead[t] = FALSE ;
    }

 /*  diskprepare()  will be called implicitly on first disk access */

    /* Hard define the tags needed for sessionInit */
  lexDefineSystemTags() ;
    /* read the tag files, if they exist, for compatibility in case lex0
       has never been saved
     */
  
  lexIsRead[0] = lexIsRead[1] = TRUE ;
  lexReHashClass(0) ;
  lexReHashClass(1) ;
}

/********************************************/

static void lexi2clear(LEXP lxp)
{ int i ;
  if (lxp)
    { lxp->dlx.dk = 0 ;
      lxp->addr = 0 ;
      lxp->cache  = 0 ;
      i = lxp - arrp(Lexi2, 0, LEXI2) ;
      lxp->key = 0 ;
      lxp->dlx.lx = minFreeLexi2 ; /* overload addr */
      minFreeLexi2 = i ;
    }
}

/************************************************/

static void lexi1Clear (int t)
{ int i ; LEXI1* lx ;

  if (Lexi1[t])
    { i = arrayMax(Lexi1[t]) ;
      lx = arrp(Lexi1[t], 0, LEXI1) - 1 ;
      while (lx++, i--)
	if (lx->is2)
	  lexi2clear (arrp (Lexi2, lx->dlx.lx, LEXI2)) ;
    }
  arrayDestroy(Lexi1[t]) ;
  keySetDestroy(LexHashTable[t]) ;
  arrayDestroy(timeStamps[t]);
}

/*************************************************/
  /* To be used only from the session manager */
void lexClear(void)
{
  register int t ;
  
  for(t=3;t< MAXTABLE; t++)
    { vocmodif[t] = lexhmodif[t] = leximodif[t] = timeModif[t] = FALSE ;
      stackDestroy(Voc[t]) ;
      if (arrayExists(Lexi1[t]))
	lexSessionStartSize[t] = arrayMax(Lexi1[t]) ;
      lexi1Clear (t) ;
      keySetDestroy(LexHashTable[t]) ;
      nClass[t] = 0 ;
      lexIsRead[t] = FALSE ;
      /* lexClearClassStatus (t, TOUCHSTATUS) ; NON: this reads the lex back ! */
    }
}

BOOL lexTableRead(KEY key)
{
  return lexIsRead[class(key)];
}

/*************************************************/

static void lexReadTable(int t)
{ 
  char lexBuffer [4096] ;
  memset (lexBuffer, 0, sizeof (lexBuffer)) ;
  if (lexIsRead[t])
    return ;
  lexIsRead[t] = TRUE ;
  
  pickList[3].type = 'A' ;  /* needed during bootstrap */
  
  if (!t || pickList[t].name)
    { 
      lexvocread(t) ;
      if(Voc[t])  /* could be 0 if class is new */
	{
	  sprintf (lexBuffer, "Reading in class %s",  pickClass2Word(t)) ;
	  messStatus (lexBuffer) ;
	  lexiread(t) ;
	  if (t<5) lexhread(t) ;
	  lexSessionStartSize[t] = arrayMax(Lexi1[t]) ;
	  lexClearClassStatus (t, TOUCHSTATUS) ;
	  lexReadTimeStamps(t);
	}
      else  /* Reinitialise this lexique */
	{ 
	  KEY dummy ;
	  lexi1Clear (t) ;
	  nClass[t] = 0 ;
	  lexSessionStartSize[t] = 0 ;
	  lexaddkey(0,&dummy,t) ;
	}
     }
}

/*************************************************/

void lexRead (void)
{ 
  BOOL debug = FALSE ;
  char lexBuffer [64] ;
  KEY key = 0 ;
  /* table 0 is the tags, created in lexDefineTags() and dynamically
     table 1 is created by hand
     table 2 is the session chooser
     table 3 is the Voc lexique, it MUST be reinitialised
     table 4 is the BAT chooser
  */
  memset (lexBuffer, 0, sizeof (lexBuffer)) ;
  if (debug)
    fprintf (stderr, "lexRead1\n") ;
  if (lexword2key ("_voc3", &key, _VGlobal) && iskeyold(key))
    { lexIsRead[_VVoc] = FALSE ; 
      lexReadTable(_VVoc) ;      
    }

  if (debug)
    fprintf (stderr, "lexRead2\n") ;
  lexReadTable(_VBat) ;

  lexIsRead[_VDisplay] = FALSE ;
  lexIsRead[_VClass] = FALSE ;

/* We now recover the correct values of the tags */
  if (debug)
    fprintf (stderr, "lexRead3\n") ;
  if (lexword2key ("_voc0", &key, _VVoc) && iskeyold(key))
    { lexIsRead[0] = FALSE ; /* we recover the tags of the previous session */
      lexReadTable(0) ;      
    }
  if (debug)
    fprintf (stderr, "lexRead4\n") ;
  lexDefineSystemTags () ; /* Redefine the system tags, which have been overwritten */
  tagInit () ;  /* finally, initialises correctly the _tag variables */
  
/* We now recover the correct values of the classes */
  if (debug)
    fprintf (stderr, "lexRead5\n") ;
  sprintf (lexBuffer, "_voc%d", _VMainClasses) ;
  if (debug)
    fprintf (stderr, "lexRead6 %s %d %d\n", lexBuffer, key, _VVoc) ;
  if (lexword2key (lexBuffer, &key, _VVoc))
    {
      if (debug)
	fprintf (stderr, "lexRead7\n") ;
      if (iskeyold(key))
	{ 
	  if (debug)
	    fprintf (stderr, "lexRead8\n") ;
	  
	  lexIsRead[_VMainClasses] = FALSE ; /* we recover the tags of the previous session */
	  lexReadTable(_VMainClasses) ; 
	}
    }
  if (debug)
    fprintf (stderr, "sysclassinit\n") ;
  sysClassInit () ; /* We redefine the system classes which may have been lost */
  classInit () ;  /* finally, initialises correctly the _VClass variables */

  sprintf (lexBuffer, "_voc%d", _VModel) ;
  if (lexword2key (lexBuffer, &key, _VVoc) &&
      iskeyold(key))
    { lexIsRead[_VModel] = FALSE ; /* we recover the tags of the previous session */
      lexReadTable (_VModel) ;
    }
  /*  pickGetClassNames() ; */
  if (debug)
    fprintf (stderr, "sysclassoptions\n") ;
  sysClassOptions() ;
} /* lexRead */

/*************************************************/

void lexReadGlobalTables(void)
{ lexIsRead[1] = lexIsRead[0] = TRUE ;
  lexiread(1) ;
  lexReHashClass(1) ;

  diskPrepare() ;  /* Reads the Global Bat */
  lexReadTable(2) ;
}

/**************************************************************/
void lexmark(int tt)    /*enables blocksubs to modify leximodif
                        * when changing a disk address
                        * note that the memory address is never
                        * saved to disk, nor the locks
                        * and that only lexaddkey and  lexRename
			* touch lexh and voc
                        */
{
  leximodif[tt] = TRUE ; 
}

/******************************************************************/
      /*writes the modified lexiques to disk*/
      /*public ; called by saveAll() and dosomethingelse()*/

void  lexSave(void) 
{ int t = MAXTABLE ;
  
  lexAlphaSaveAll () ;
  while(t--)  /* DO save lex[0]. tags.wrm only read for compatibility if lex0 is missing */
    {
      if(t==1)
	lexiwrite(t);  /*t = 1 Global voc */
      else 
	if(Lexi1[t])
	  {
	    lexvocwrite(t);
	    lexiwrite(t);
	    lexhwrite(t);
	    lexWriteTimeStamps(t);
	      
	  }
    }
}

/********************************************/

static void lexvocread(int t)
{ KEY key ;
  char buffer[24] ;
  sprintf(buffer, "_voc%d", t) ;
  lexaddkey (buffer,&key, !t || t>_VVoc ? _VVoc : _VGlobal ) ;
  lexaddkey (buffer,&key, !t || t>_VVoc ? _VVoc : _VGlobal ) ;
  stackDestroy(Voc[t]) ;
  Voc[t] = stackGet(key) ;
  vocmodif[t] = FALSE ;
}

/**************************************************************/

static void lexvocwrite(int t)
{
  char lexBuffer [64] ;
  KEY key ;

  if (!isWriteAccess() || !vocmodif[t])
    return ;

  if(PRINT)
    printf("\n   writingVoc %d",t);
  sprintf (lexBuffer, "_voc%d",t) ;
  lexaddkey (lexBuffer,&key,!t || t>_VVoc ? _VVoc : _VGlobal ) ;
  lexaddkey(lexBuffer,&key,!t || t>_VVoc ? _VVoc : _VGlobal ) ;
  stackStore(key,Voc[t]) ;
  
  vocmodif [t] = FALSE ;
}

/**************************************************************/

static void lexhread(int t)
{
  char lexBuffer [64] ;
  KEY key ;

  sprintf (lexBuffer, "_lexh%d",t) ;
  lexaddkey(lexBuffer,&key,!t || t>_VVoc ? _VVoc : _VGlobal ) ;
  lexaddkey(lexBuffer,&key,!t || t>_VVoc ? _VVoc : _VGlobal ) ;
  keySetDestroy(LexHashTable[t]) ;
  LexHashTable[t] = arrayGet(key, KEY,"k") ;
  if (!LexHashTable[t])
    lexReHashClass(t) ;
  else
    { int n = arrayMax(LexHashTable[t]) , nBits = 0 ;
      while ((1 << ++nBits) < n) ;
      if(n != (1<<nBits))
	messcrash("Wrong size in lexhread %d", t) ;
      nClass[t] = nBits ;
    }

  lexhmodif [t] = FALSE ;
}

/**************************************************************/

static void lexhwrite(int t)
{
 /*  Commented out since i want to change the hashText 
     routine till it is really good
  char lexBuffer [64] ;
  KEY key ;
  
 
  if (!isWriteAccess())
    return ;

  if(!lexhmodif[t])
    return ;
  if(PRINT)
    printf("\n   writinglexh %d",t);
    sprintf (lexBuffer, "_lexh%d",t) ;
  lexaddkey(lexBuffer,&key,!t || t>_VVoc ? _VVoc : _VGlobal ) ;
  lexaddkey(lexBuffer,&key,!t || t>_VVoc ? _VVoc : _VGlobal ) ;
  arrayStore(key,LexHashTable[t]) ;
*/  
  lexhmodif[t] = FALSE ;
  return ;
}

/**************************************************************/

static void lexiread(int t)
{
  char lexBuffer [64] ;
  KEY key ;
  register LEXI1* lxi1 ;
  register int j ;
  Array dummy ;
  
  sprintf (lexBuffer, "_lexi%d",t) ;
  lexaddkey(lexBuffer,&key,!t || t>_VVoc ? _VVoc : _VGlobal ) ;
  dummy = arrayGet(key, LEXI1,lexiFormat) ;
  lexi1Clear (t) ;
  Lexi1[t] = dummy ;
  if (!Lexi1[t])
    messcrash("Lexi[%d] not found",t);
  
  lxi1 = arrp(Lexi1[t], 0, LEXI1) - 1;
  j = arrayMax(Lexi1[t]);
  while(lxi1++, j--)
    { lxi1->is2 = 0 ;
      lxi1->lock &= ~LOCKSTATUS ; /* NOT 3, to zero the last 2 bits */
      if (lxi1->lock & ALIASSTATUS)
	{
	  key = lxi1->dlx.alias ;
	  /* lx2 =  arrp(Lexi1[t], KEYKEY(key), LEXI1) ; unused */
	  lxi1->lock  |= ISALIASSTATUS ;
	}
    }
  leximodif[t] = FALSE ;
}

/**************************************************************/
    /* i must restore the disk address in the union */
static void lexiwrite(int t)
{ 
  char lexBuffer [64] ;
  Array lx ;
  KEY key ;
  DISK dsk ;
  LEXP1 lxp ;
  int i, j ;
  BOOL   needCopy = FALSE ;

  if (!isWriteAccess() || !leximodif[t]) 
    return ;
  sprintf (lexBuffer, "_lexi%d",t) ;
  lexaddkey(lexBuffer,&key,!t || t>_VVoc ? _VVoc : _VGlobal ) ;

  if(PRINT)
    printf("\n   writinglexi %3d,  %8d  entries.",
	   t, arrayMax(Lexi1[t]));  
  leximodif[t] = FALSE ;  /* must come before the write */

  lx = Lexi1[t] ;
  i = arrayMax(lx) ;
  lxp = arrp(lx, 0, LEXI1) - 1 ;

  while (lxp++, i--)
    if (lxp->is2)
      { needCopy = TRUE ; break ; }

  if (needCopy)  /* restore the disk address in dlx union */
    { lx = arrayCopy (Lexi1[t]) ;
      i = arrayMax(lx) ;
      lxp = arrp(lx, 0, LEXI1) - 1 ;
      while (lxp++, i--)
	if (lxp->is2)
	  { j = lxp->dlx.lx ;
	    if (j > arrayMax(Lexi2) || j <= 0)
	      messcrash("lexi1/2 inconsistency") ;
	    dsk = arrp(Lexi2, j, LEXI2)->dlx.dk ;
	    lxp->is2 = 0 ;
	    lxp->dlx.dk = dsk ;
	  }
    }     

  arrayStore(key,lx,lexiFormat) ;
  if (lx != Lexi1[t]) 
    arrayDestroy(lx) ;
}

/*************************************************************************/

static void lexReadTimeStamps(int t)
{ 
  char lexBuffer [64] ;
  KEY key;

  sprintf (lexBuffer, "_lext%d",t) ;
  lexaddkey (lexBuffer,&key,!t || t>_VVoc ? _VVoc : _VGlobal ) ;
  
  timeStamps[t] = arrayGet(key, TIMESTAMP, timeStampFormat);

  /* If this class has never had timestamps, make some now. */
  if (!timeStamps[t])
    timeStamps[t] = arrayCreate(arrayMax(Lexi1[t]), TIMESTAMP);
  
  /* below we deal with two cases, 
     1) adding timestamps to a class for the first time.
     2) objects having been added by non timestamp aware binaries. 
     In both cases we fill in the missing creation and update 
     timestamps with those of the current kernel session. */
  
  if (arrayMax(Lexi1[t]) > arrayMax(timeStamps[t]))
    { 
      int i;
      if (0 && t>3)
	messdump("Adding missing object-timestamps for class %s",
		 pickClass2Word(t));
      for (i = arrayMax(timeStamps[t]); i<arrayMax(Lexi1[t]); i++)
	{
	  array(timeStamps[t], i, TIMESTAMP).created = thisSession.userKey ;
	  arr(timeStamps[t], i, TIMESTAMP).updated = thisSession.userKey ;
	}
      timeModif[t] = TRUE;
    }
  else 
    timeModif[t] = FALSE;

  return ;
}

static void lexKeyUpdated(KEY key)
{
  int t = class(key);
  int i = KEYKEY(key);

  if (t<=1)
    return;

  if (!timeStamps[t])
    lexReadTimeStamps(t);

  if (i>arrayMax(timeStamps[t]))
    return;

  array(timeStamps[t], i, TIMESTAMP).updated = sessionUserKey();
}

static void lexWriteTimeStamps(int t)
{
  char lexBuffer [64] ;
  KEY key ;
 
  if (!isWriteAccess())
    return ;

  if(!timeModif[t])
    return ;
  sprintf(lexBuffer, "_lext%d",t) ;
  lexaddkey(lexBuffer,&key,!t || t>_VVoc ? _VVoc : _VGlobal ) ;
 
  arrayStore(key, timeStamps[t], timeStampFormat);

  timeModif[t] = FALSE;
}

void lexSetCreationUserSession(KEY key, KEY timeStamp)
{
  int t = class(key);
  int i = KEYKEY(key);
 
  if (t<=1)
    return;

  if (!timeStamps[t])
    lexReadTimeStamps(t);

  if (i>arrayMax(timeStamps[t]))
    return;

  arr(timeStamps[t], i, TIMESTAMP).created = timeStamp;
  timeModif[t] = TRUE;

}

KEY lexCreationUserSession(KEY key)
{
  int t = class(key);
  int i = KEYKEY(key);
 
  if (t<=1)
    return 0;

  if (!timeStamps[t])
    lexReadTimeStamps(t);

  if (i>arrayMax(timeStamps[t]))
    return 0;

  return arr(timeStamps[t], i, TIMESTAMP).created;
}

mytime_t lexCreationStamp(KEY key)
{
  return sessionUserSession2Time(lexCreationUserSession(key));
}

KEY lexUpdateUserSession(KEY key)
{
  int t = class(key);
  int i = KEYKEY(key);

  if (t<=1)
    return 0;

  if (!timeStamps[t])
    lexReadTimeStamps(t);

  if (i>arrayMax(timeStamps[t]))
    return 0;

  return arr(timeStamps[t], i, TIMESTAMP).updated;
}

mytime_t lexUpdateStamp(KEY key)        
{
  return sessionUserSession2Time(lexUpdateUserSession(key));
}
/******************************************************************/
  
void lexOverLoad(KEY key,DISK disk) /* used in session.c only, for bootstrap */
{ LEXP lx2 = KEY2LEX2(key) ;
  LEXP1 lx1 = KEY2LEX1(key) ;

   if (lx2)
     { lx2->dlx.dk = disk ;
       lx2->cache = 0 ;
       lx2->addr = 0 ;
       lx2->key = key ;
     }
  else
    lx1->dlx.dk = disk ;
}

/******************************************************************/

DISK lexDisk(KEY key)
{ LEXP lx2 = KEY2LEX2(key) ;
  LEXP1 lx1 = KEY2LEX1(key) ;

   return 
     lx2 ? lx2->dlx.dk : lx1->dlx.dk ;
}

/******************************************************************/
static LEXP1 myKey2Lex1 = 0 ;
static LEXP myKey2Lex2 = 0 ;

void lexSetDisk(KEY key, DISK dk)
{ DISK old = 0 ; 
  iskey(key) ;

  if (myKey2Lex2)  
    { 
      old = myKey2Lex2->dlx.dk ; myKey2Lex2->dlx.dk = dk ;
    }
  else if (myKey2Lex1) 
    { 
      old = myKey2Lex1->dlx.dk  ; myKey2Lex1->dlx.dk = dk ; 
    }
  else messcrash ("bad call to lexSetDisk(%s)", name(key)) ;

  if (dk != old) 
    {
      leximodif[class(key)] = TRUE ;
      if (class(key)>1)
	{ 
	  timeModif[class(key)] = TRUE;
	  array(timeStamps[class(key)], KEYKEY(key), TIMESTAMP).updated = 
	    sessionUserKey();
	}
    }
}

/******************************************************************/
   /* Remove leading and trailing spaces and non ascii */
   /* rd 971009: remove leading/trailing \t, change other \t to ' '
      return result in private buffer (i.e. do not change orginal)
   */

/* check that name has some non-blank chars */
BOOL lexIsGoodName(char *name)
{
  if (!name) 
    return FALSE;

  while (*name) 
    if (!isspace((int)*name++))
      return TRUE;
  
  return FALSE;
}

char *lexcleanup (const char *cp0, AC_HANDLE handle )
{
  char *cq, *cp ;
  
  while (*cp0 == ' ' || *cp0 == '\t')
    cp0++ ;
  
  cp = strnew (cp0, handle);
  
  for (cq = cp + strlen(cp) ; cq-- > cp ; )
    if (*cq == ' ' || *cq == '\t')
      *cq = 0 ;
    else
      break ;
  
  return cp ;
}

/******************************************************************/

void lexHardDefineKey (int table, KEY key,  char *cp) 
{ 
  char lexBuffer [64] ;
  int i ;
  KEY  key2 ;
  char *cq ;

  if(lexword2key(cp,&key2,table))
    {
      if (key == KEYKEY(key2))
	return ;
      else
	messcrash
	  ("The old %d:%d and new %d:%d values of key %s differ.\n%s\n%s",
	   table,KEYKEY(key2), table, key, cp,
	   "As of ace.3.4 files wspec/sys* and tags can be removed",
	   "before creating the database, but not later on") ;
    }

  i = lexMax(table) ;
  if(!i) i = 1 ; /* Because a zeroth key is created by lexaddkey */

  if (!key) 
    key = i ;
  if(i < key)
    { for(;i < key; i++)
      {
	sprintf (lexBuffer, "__sys%d",i) ;
	lexaddkey(lexBuffer,&key2,table) ;
      }
      lexaddkey(cp,&key,table);
    }
  else if(i == key)
    lexaddkey(cp,&key,table);
  else
    { 
      cq = name(KEYMAKE(table,key)) ;
      sprintf (lexBuffer, "__sys%d",key) ;
      if(key != _LongFloat && key != _LongInt /* introduced september 2019 */
	 && key != _zLong && key != _zLong_Unsigned /* introduced september 2019 */
	 && strcmp(cq,lexBuffer))
	messcrash("Tag %s = %d tries to overwrite tag %s",
		  cp, key, cq) ;
      /* Else I am overwriting a dummy tag */
      array(Lexi1[table], (int) key, LEXI1).nameoffset 
	= stackMark(Voc[table]) ;
      cq =  lexcleanup(cp, 0) ;
      pushText(Voc[table], cq) ;
      messfree (cq) ;
      
      leximodif[table] = lexhmodif[table] = vocmodif[table] = TRUE ;
      lexAlphaMark (table) ;
    }
}

/********************************************/
                      /* Returns 0 if kk is unknown,
                       * 1 if kk is just a vocabulary entry
                       * 2 if kk corresponds to a full object
		       *
		       * sets myKey2Lex
                       */
int isNakedKey(KEY kk)
{
  KEY k = KEYKEY(kk);
  int t = class(kk);
  if (!kk || (k >= lexMax(t)) )
   { myKey2Lex1 = 0 ;  myKey2Lex2 = 0 ; return 0 ; }

  myKey2Lex1 = arrp(Lexi1[t], k, LEXI1);
  myKey2Lex2 = myKey2Lex1->is2 ? arrp(Lexi2, myKey2Lex1->dlx.lx, LEXI2) : 0 ;

  if (myKey2Lex2 && myKey2Lex2->key != kk)
    messcrash ("inconsistency in isNakedkey") ;

  return
    myKey2Lex1->dlx.dk || myKey2Lex2  ?  2  :  1  ; /* disk or lex2 both mean full object */
}

/********************************************/
   /* prevent the formation of alias loops */
BOOL lexAliasLoop (KEY old, KEY new)
{
  KEY kk = new ;
  int result ;
  
  while (TRUE)
    { result = isNakedKey(kk) ;
      if (!result ||
	  ! (myKey2Lex1->lock & ALIASSTATUS))
	break ;
      kk = myKey2Lex1->dlx.alias ;
      if (kk == old)
	 /* to alias old to new would create a loop */
	return TRUE ;
    }
  return FALSE ;
}

#ifdef NOT_NEED_BACUSE_APPERENTLY_NOT_USED
/********************************************/
   /* Scans the alias list for occurence of newName */
static BOOL lexAliasExists(KEY key, char* newName)
{
  KEY kk = key ;
  char *name ;
  BOOL isCaseSensitive = pickList[class(key) & 255 ].isCaseSensitive ;
/*   int (*lexstrIsCasecmp)(const char *a,const char *b) = isCaseSensitive ? */
/*    lexstrCasecmp : lexstrcmp ; */
  int (*lexstrIsCasecmp)() = isCaseSensitive ?
    strcmp : strcasecmp ;
 
  while (TRUE)
    { name = nakedName(kk) ;  /* sets myKey2lex */
      if (!lexstrIsCasecmp(name, newName)) /* Equivalence  */
	return TRUE ;
      if (!(myKey2Lex1->lock & ALIASSTATUS))
      	return FALSE ;
      kk = myKey2Lex1->dlx.alias ;
    }
}
#endif /* NOT_NEED_BACUSE_APPERENTLY_NOT_USED */


/********************************************/
/* Called from parse.c, removes an aliased key and relinks the alias list */
BOOL lexDestroyAlias (KEY key)
{ KEY  under = 0 ;
  LEXI1* lx ;
  int t, i ;


      if (!isNakedKey(key) ||                     /* sets myKey2lex */
	  ! (myKey2Lex1->lock & ALIASSTATUS))
	return FALSE ;
      under = myKey2Lex1->dlx.alias ;
      myKey2Lex1->lock &= ~ALIASSTATUS ;  
      myKey2Lex1->lock |= EMPTYSTATUS ;
      myKey2Lex1->dlx.alias = 0 ;
      t = class(key) ;
      leximodif[t] = TRUE ;
      lexAlphaMark (t) ;
      lx = arrp(Lexi1[t], 0, LEXI1) - 1 ;
      i = arrayMax(Lexi1[t]) ;
      while(lx++, i--)
	if ((lx->lock & ALIASSTATUS) && lx->dlx.alias == key)
	  lx->dlx.alias = under ;
      return TRUE ;

}

/********************************************/
/* returns 0:no, 1:may_be, 2:yes, key is alias of another key 
*/

int lexIsAlias (KEY key)
{
  if (!lexIsRead [class(key)])
    return 1 ; /* may be */
  isNakedKey(key) ;
  if (myKey2Lex1->lock & ISALIASSTATUS)
    return 2 ;
  return 0 ;
}

/********************************************/
/* Scans the alias list for valid object, used by bsCreate etc  */
KEY lexAliasOf (KEY key)
{
  while (TRUE)
    {
      if (pickList[(class(key)) & 255 ].protected ||
	  !KEYKEY(key) ||	  
	  !isNakedKey(key) ||                     /* sets myKey2lex */
	  ! (myKey2Lex1->lock & ALIASSTATUS))
	return key ;
      key = myKey2Lex1->dlx.alias ;
    }
}

/***************************************************/

/* returns false if key does not exist, has aliasstatus or emptystatus, 
   is in class 0 (a tag), a private kernel class  or is a KEYKEY 0 model */
/* cannot use lexGetStatus for this as it follows the alias list */

BOOL lexIsKeyVisible (KEY key)
{
  if (!isNakedKey(key) ||
      (myKey2Lex1->lock & (ALIASSTATUS | EMPTYSTATUS)))
    return FALSE ;

  if (pickList[class(key)].private)		/* i.e. a tag */
    return FALSE ;

#ifndef NEW_MODELS
  if (KEYKEY(key) == 0 && 
      pickList[class(key)].type == 'B')	/* i.e. a model */
    return FALSE ;
#endif

  return TRUE ;
}

/********************************************/
   /* Iterates through the alias list */
int iskey(KEY kk)
{
  KEY start = kk ;
  int result ;
 
  while (TRUE)
    { result = isNakedKey(kk) ;
      if (!result ||
	  ! (myKey2Lex1->lock & ALIASSTATUS))
	break ;
      kk = myKey2Lex1->dlx.alias ;
      if (kk == start)
	messcrash ("Loop in alias list of Key %d : %s",
		   kk, nakedName(kk) ) ;
    }
  return result ;
}

/********************************************/

/* Returns TRUE
 * if kk has a primary cache or a disk address
 
 never use this function from higher layers of the code
 it is only imteresting for the kernel, when manipulating cache 1
*/

BOOL iskeyold(KEY kk)
{
  if (!iskey(kk))
    return FALSE ;
  if (!myKey2Lex2)
    return myKey2Lex1->dlx.dk ? TRUE : FALSE ;
  return
    myKey2Lex2->addr || myKey2Lex2->dlx.dk ? TRUE : FALSE ;
}

/********************************************/

static LEXP1 KEY2LEX1(KEY kk)
{
  return
    iskey(kk) ?  myKey2Lex1 : 0 ;
}

/********************************************/

static LEXP KEY2LEX2(KEY kk)
{
  return
    iskey(kk) ?  myKey2Lex2 : 0 ;
}

/********************************************/
/* for external code ; creates the lex2 entry */

static int newLexi2(void)
{ LEXP lx2 ;
  int i , j ;
  
  i = minFreeLexi2 ;
  if (!i)
    { i = arrayMax (Lexi2) ;
      array (Lexi2, i + 1000, LEXI2).dlx.lx = 0 ;
      for (j = i , lx2 = arrp (Lexi2, j, LEXI2) ; j < i + 1000 ; j++, lx2++)
	lx2->dlx.lx = j + 1 ;
    }
  lx2 = arrayp(Lexi2, i, LEXI2) ;
  minFreeLexi2 = lx2->dlx.lx ;
  lx2->dlx.dk = 0 ; /* un-overload dlx */
  return i ;  
}

LEXP KEY2LEX(KEY kk)
{ DISK disk ;
  int j2 ;
  int result;

  result = iskey(kk);		/* sets myKey2Lex1,2 */
  if (!result)
    messcrash("KEY2LEX called with invalid KEY %s", name(kk));
    
  if (myKey2Lex2)
    return myKey2Lex2 ;
  if (!myKey2Lex1)
    return 0 ;
  
  disk = myKey2Lex1->dlx.dk ;
  myKey2Lex1->dlx.lx = j2 = newLexi2() ;
  myKey2Lex1->is2 = 1 ;
  myKey2Lex2 = arrayp(Lexi2, j2, LEXI2) ;
  myKey2Lex2->dlx.dk = disk ;
  myKey2Lex2->cache = 0 ;
  myKey2Lex2->addr = 0 ;
  {		/* RD 980807 : important bug fix
		   must follow aliases!!
		   NB can't use lexAliasOf() because crash during bootstrap
		   but LEXI2->key only needed for checking
		   so could remove it completely - it was bugged!
		*/
    while (arrp(Lexi1[class(kk)], KEYKEY(kk), LEXI1)->lock & ALIASSTATUS)
      kk = arrp(Lexi1[class(kk)], KEYKEY(kk), LEXI1)->dlx.alias ;
  }
  myKey2Lex2->key = kk ;
  return
    myKey2Lex2 ;
}

/********************************************/

BOOL lexlock(KEY key)
{      /* code partially  duplicated in cacheLock*/
 LEXP1 q=KEY2LEX1(key);
  if(q)
    {
      if(q->lock & LOCKSTATUS)
	messerror("Double locking  of %s",
		  name(key));
      else
	{
	  q->lock |= LOCKSTATUS; return TRUE;
	}
    }
  return FALSE;
}

/********************************************/

void lexunlock(KEY key)
{
 LEXP1 q=KEY2LEX1(key);
  if(   q
     && (q->lock & LOCKSTATUS) )
    q->lock &= ~LOCKSTATUS; 
  else
    messcrash("Unbalanced unlock of %s",
	      name(key));
}

/********************************************/
   /* subClass mechanism */
void lexSetIsA (KEY key, unsigned char c) 
{     
  LEXP1 q=KEY2LEX1(key);
  if(q)
    q->isMask = c ;
}

/********************************************/
   /* subClass mechanism: test appartenance */
BOOL lexIsInClass(KEY key, KEY classe)
{ unsigned char mask ;   
  LEXP1 q=KEY2LEX1(key);

  if (classe < 256)
    return classe == class(key) ;

  if (!q)
    return FALSE ;

  pickIsA(&classe, &mask) ;
  return 
    classe == class(key)  &&
      (mask & q->isMask) == mask ;
}

/********************************************/
  /* alteration of LEXPRIVATESTATUS flags would 
     foul the lex system */
  /* code duplicated in cacheUpdate
     jan.97: i think this is no longer true
  */

/* do not mark as touched an object obtained from the server */
BOOL X_CLIENT_PARSING = FALSE ;

void lexSetStatus(KEY key, unsigned char c) 
{     
  LEXP1 q=KEY2LEX1(key);

  if (c == TOUCHSTATUS)
    { if (X_CLIENT_PARSING)
	return ;
      if (touchedByClient) 
	keySet(touchedByClient, keySetMax(touchedByClient)) = key;
    }

  if (q)
    q->lock |= (c & ~LEXPRIVATESTATUS) ;
  else
    messcrash("lexSetStatus called with bad key %s", name(key));
}

/********************************************/
  
  /* code duplicated in cacheUpdate
      jan.97: i think this is no longer true
  */
void lexUnsetStatus(KEY key, unsigned char c) 
{     
  LEXP1 q=KEY2LEX1(key);

  if (q)
    q->lock &= (~c) | LEXPRIVATESTATUS ;
  else
    messcrash("lexUnsetStatus called with bad key %s",
	      name(key));
}

/********************************************/
  
void lexClearClassStatus(KEY classe, unsigned char c) 
{ int i ;        
  LEXP1 q ;
 
  c = (~c) | LEXPRIVATESTATUS ;
  if ((i = lexMax(classe)))
    { q = arrp(Lexi1[classe], 0, LEXI1) - 1 ;
      while (q++, i--)
	q->lock &= c ;
    }
}

/********************************************/

unsigned char lexGetStatus(KEY key)
{ 
  LEXP1 q=KEY2LEX1(key);
  if(q)
    return q->lock ;
  else
    messcrash("lexGetStatus called with bad key %s",
	      name(key));
  return 0 ;  /* for compiler happiness */
}

/********************************************/
   /* Avoid the alias system, used by word2key and hasher */
static char * nakedName(KEY kk)  
{     
  if (!isNakedKey(kk))
    return  "\177\177(NULL KEY)" ;
  
  if (pickType (kk) == 'D')
    return (name(myKey2Lex1->nameoffset)) ;

  return
    stackText(Voc[class(kk)],
	      myKey2Lex1->nameoffset) ;
}

/********************************************/

BOOL lexKillDkey (KEY kk, KEY parent)
{
  if (!iskey(kk) ||
      myKey2Lex1->nameoffset != parent)
    return FALSE ;
  arrayKill (kk) ;
  myKey2Lex1->nameoffset = 0 ;
  return TRUE ;
}

/********************************************/

KEY lexCreateDkey (int t, KEY parent)
{
  KEY k ;
  LEXI1 *aip ;

  if (pickType (t<<24) == 'D')
    messcrash ("Attempt to create a D key, not in a D class") ;

  if (!lexIsRead[t])
    lexReadTable(t) ;
    
  k = (KEY) arrayMax(Lexi1[t]) ;

  aip = arrayp(Lexi1[t], (int) k, LEXI1) ;
  aip->nameoffset = parent ;
  aip->dlx.dk = 0 ;
  aip->lock = 0;
  aip->isMask = 0;
  aip->is2 = 0 ;

  leximodif[t] = TRUE ;
  /* notice that we do not need to hasha D class */
  return KEYMAKE (t, k) ;
}

/********************************************/

char *name(KEY kk)  
{   
  
  if (!iskey(kk))
    return  "\177\177(NULL KEY)" ;
  
  if (pickType (kk) == 'D')
    return (name(myKey2Lex1->nameoffset)) ;

  return
    stackText(Voc[class(kk)],
	      myKey2Lex1->nameoffset) ;
}

char* nameWithClass (KEY k)
{
  static char *buf = 0 ;
  static int buflen = 0 ;
  char *nam = name (k) ;
  char *classNam = className (k) ;

  if (strlen(nam) + strlen(classNam) + 2 > buflen)
    { if (buf) messfree (buf) ;
      buflen = strlen(nam) + strlen(classNam) + 2 ;
      buf = (char*) messalloc (buflen) ;
    }

  sprintf (buf, "%s:%s", classNam, nam) ;
  return buf ;
}

/********************************************/
  /* Iterates through the alias list */
BOOL nextName(KEY key, const char **cpp)  
{
  static LEXP1 lxp = 0 ;
  static KEY lastKey = 0 , kk ;
  
  if ( *cpp && key != lastKey)
    messcrash ("nextName %s called out of context", nakedName(key)) ;
 
  if (!*cpp)    /* initialise */
    kk = lastKey = key ;
  else         /* loop */
    { if (!lxp || !(lxp->lock & ALIASSTATUS)) /* shift once */
	return FALSE ;
      kk = lxp->dlx.alias ;
    }

  while (TRUE)
    { if (!isNakedKey(kk))
	return FALSE ;
      if (!(myKey2Lex1->lock & EMPTYSTATUS))   /* if empty */
	break ;
      if (!(myKey2Lex1->lock & ALIASSTATUS))  /* loop again */
	return FALSE ;
      kk = myKey2Lex1->dlx.alias ;
    }
  lxp = myKey2Lex1 ;                         /* used by shift once */
  *cpp = nakedName(kk) ;

  return TRUE ;
}

/********************************************/

BOOL lexReClass(KEY key,KEY *kp, int classe)
{
  if (classe > 0 && classe < MAXTABLE)
    {
      PICKLIST *p = pickList + classe ;
      if (p->name && p->visible && p->type == 'B' && ! p->protected && p->model)
	return 
	  lexword2key(name(key),kp,classe) ;
    }
  return FALSE ;    
}

/**********************************************/

static BOOL colonParse (char* note, char **classe, char **item)
{
  char *cp ;

  for (cp = note ; *cp != ':' ; ++cp) 
    if (!*cp)
      return FALSE ;

  *cp = 0 ;
  *classe = note ;
  *item = ++cp ;
  return TRUE ;
}

BOOL lexClassKey (char *text, KEY *kp)
{ 
  int classe ;
  char *cName, *kName ;

  *kp = 0 ;
  return
    colonParse (text, &cName, &kName) &&
    (classe = pickWord2Class (cName)) &&
    ((*(kName - 1) = ':')) &&     /* Restore text */
    lexword2key (kName, kp, classe) ;
}

/********************************************/
/********************************************/
/************ Hash Package ******************/

                     /* Hashing system */
#define HASH_RATE 1.6
#define  SIZEOFINT (8 * sizeof (int))

unsigned int hashText(char *cp, int n, BOOL isOdd)
{
  unsigned int 
    i, h = 0 , j = 0,
    rotate = (isOdd ? 21 : 13) ,
    leftover = SIZEOFINT - rotate ;

  /*  chrono("hashText") ;   */
  while(*cp)
    { h = ace_upper(*cp++) ^  /* XOR*/
	( ( h >> leftover ) | (h << rotate)) ; 
    }
  /* compress down to n bits */
  for(j = h, i = n ; i< SIZEOFINT; i+=n)
    j ^= (h>>i) ;
  /* chronoReturn() ; */

  if (isOdd) /* return odd number */
    j |= 0x1 ;
  return
     j &  ( (1<<n) - 1  )  ;
}

/***************************************************/

 void lexReHashClass(int classe)
{
  KEY key = 0 ; int n ;

  chrono("lexReHashClass") ;

  n = 7 ;  /* chose size */
  while((1 << ++n) <= HASH_RATE * lexMax(classe)) ;
  nClass[classe] = n ;
  
 
  LexHashTable[classe] = arrayReCreate(LexHashTable[classe], 1<<n, KEY) ;
  keySet(LexHashTable[classe], ( 1 << n ) - 1) = 0 ; /* make room */

  key = arrayMax(Lexi1[classe]) ;
  while (key--)
    lexHashInsert(classe, key) ;
  chronoReturn() ;
}

/****************************************************/
   /* should only be called if !lexword2key */
   /* key must be just the KEYKEY part */
int NII  , NIC  , NWC  , NWF  , NWG  ;
/* Simon Kelley sent this fix. so i have put it in (il) */
static void lexHashInsert(int classe, KEY key)
{
  static KEYSET hashTable ;
  unsigned int h, dh = 0 ;
  int max, n , nn ; 
  KEY *kp, *kpMax ;
  char *cp ;
  BOOL failedFirstPass = FALSE;

  /* chrono("lexHashInsert") ; */
  key = KEYKEY (key) ;
  hashTable  = LexHashTable[classe] ;
retry:
  if(failedFirstPass || !hashTable ||
     (arrayMax(Lexi1[classe]) * HASH_RATE > keySetMax(hashTable)))
    lexReHashClass(classe) ;

  cp = nakedName(KEYMAKE(classe,key)) ;
  if (!*cp)
    {
      /* chronoReturn () ; */
      return ;  /* may occur after a rename */
    }
  h = hashText(cp, nClass [classe], FALSE) ;

  nn = nClass[classe] ;
  n = max = 1 << nn ;

  kp = arrp(hashTable, h , KEY) ;
  kpMax = arrp(hashTable, n -1, KEY) ;
 
  while ( n--)  /* circulate thru table, starting at h */
    { if(!*kp)      /* found empty slot, do insert */
        { *kp = key ;
          NII++ ;
          /* chronoReturn() ; */
          return ;
        }
        /* circulate */
      if (!dh)
        dh = hashText(cp, nClass [classe], TRUE) ;

      kp += dh ;
      if (kp > kpMax)
            kp -= max ;
      NIC++ ;
    }
  if (failedFirstPass) 
    messcrash("Hash table full, NIC = %d, lexhashInsert is bugged", NIC)
;
  else
    { /* we may fail to find a slot because of the presence of stale
keys */
      /* after lexHardRename has run, try one more time, having purged
*/
      /* the stale keys by re-hashing - srk */
      failedFirstPass = TRUE;
      goto retry;
    }
}

/*
  I comment out this function 
  it needs plot so prevents acequery froom linking

  Note that i am not sure that it coincides with the 
  current strategy in hashInsert and word2key

#include "plot.h"
void lexHashTest(int classe)
{
  static KEYSET hashTable ;
  unsigned int h ;
  int n, nn , nn1 ;
  KEY *kp, *kpMin, *kpMax ;
  int empty = 0 , full = 0 , j ;
  Array histo = arrayCreate(200,int) ;
  char *cp ;

  hashTable  = LexHashTable[classe] ;
  if(!hashTable || !arrayMax(hashTable))
    {  messerror("empty table\n") ;
       return ;
     }
  
  kpMin = arrp(hashTable, 0 , KEY) ;
  kpMax = arrp(hashTable, n -1, KEY) ;
  nn = nClass[classe] ;
  nn1 = n = arrayMax(hashTable) ;
  for(h=0; h <nn; h++)
    { 
      kp = arrp(hashTable, h , KEY) ;
      j = 0 ;
      if (!*kp)
	empty++ ;
      else
	{ full++ ;
	  cp = nakedName(KEYMAKE(classe, *kp)) ;
	  n = nn ;
	  while ( n--)  
	    { if(!*kp)      
		{
		 array(histo,j,int)++ ;
		 break ;
		}

	     // circulate 
	      if (!dh)
		dh = hashText(cp, nClass [classe], TRUE) ;
	      kp += dh ;
	      if (kp > kpMax)
		kp -= max ;
	      j++ ;
	    }
	}
    }
  plotHisto
    (messprintf ("Hash correl %s", pickClass2Word(classe)), histo) ;
}
*/

/********************************************/

static void lexHardRename (KEY key, const char *newName)
{
  char 
    *cp, *oldName = nakedName(key) ;  /* sets myKey2Lex */
  int classe = class(key) , n = strlen(oldName) ;

  for (cp = oldName ; *cp ;)
    *cp++ = 0 ;  /* will be left out of hash on next rehashing */
  
  if (strlen(newName) > n )
    { myKey2Lex1->nameoffset = stackMark(Voc[classe]) ;
      leximodif[classe] = TRUE ;	/* since nameoffset is changed */
      pushText (Voc[classe], newName) ;
    }
  else
    strcpy (oldName, newName) ;
    
  if (myKey2Lex1->lock & EMPTYSTATUS)  /* because the new name should appear in lists */
    { myKey2Lex1->lock &= ~EMPTYSTATUS ;
      leximodif[class(key)] = TRUE ;
    }
  lexHashInsert(classe, key) ; /* re-insert in hash table */
  
  vocmodif[classe] = TRUE ;
}
  
/********************************************/

void lex2Check (void)
{return ;
/*
  LEXP lx = arrp(Lexi2, 0, LEXI2) - 1 ;
  int i = arrayMax(Lexi2) ;

  while (lx++, i--)
  if (lx->cache && (*(KEY*)lx->cache) != lx->key)
  messcrash("lex2Check error") ;
  */
}

void lex2clear (KEY key)
{ iskey (key) ;
  if (!myKey2Lex2 ||
      myKey2Lex2->cache ||
      myKey2Lex2->addr)
    return ;
/*   lex2Check() ; */
  myKey2Lex1->dlx.dk = myKey2Lex2->dlx.dk ;
  myKey2Lex1->is2 = 0 ;
  lexi2clear(myKey2Lex2) ;
}

/********************************************/

static void lexMakeAlias(KEY key, KEY alias, BOOL keepOldName)
{ 
  isNakedKey(key) ;
  
  myKey2Lex1->dlx.alias = alias ; /* overloading dlx union with alias key */
  myKey2Lex1->lock = ALIASSTATUS ;  
  if (!keepOldName)
    myKey2Lex1->lock |= EMPTYSTATUS ;  

      /* destroy the lex2 entry */
  myKey2Lex1->is2 = 0 ;
  if (myKey2Lex2)
    lexi2clear(myKey2Lex2) ;
      
  if (!keepOldName)              /* Clean up */
    { char *cp = nakedName(key) ;
      while (*cp) *cp++ = 0 ;
      myKey2Lex1->nameoffset = 0 ;
    }
  
  isNakedKey(alias) ;
  myKey2Lex1->lock &= ~EMPTYSTATUS ;
  myKey2Lex1->lock |= ISALIASSTATUS ;
  leximodif[class(key)] = TRUE ;
}
 
/********************************************/
  /* Returns TRUE if newName is accepted */
  /* if keepOldName, the 2 names will be used in the future,
   * so we need 2 keys anyway,
   * If not, we would rather keep a single key except if the old 
   * and new one allready exist because they may be referenced from
   * elsewhere.
   */

BOOL lexAlias(KEY *keyp, const char* newName, BOOL doAsk, BOOL keepOldName)
{
  char *oldName, *oldnamepreserve = 0 ;
  int  classe = *keyp ? class(*keyp) : 0 ;
  KEY  newKey, key ;
  extern BOOL bsFuseObjects(KEY key, KEY new) ;
  BOOL isCaseSensitive = pickList[classe & 255 ].isCaseSensitive ;
  int (*lexstrIsCasecmp)(const char *, const char *) = isCaseSensitive ?
    strcmp : strcasecmp ;
  AC_HANDLE h = 0 ;

     /********** Protections *************/ 

  if (!keyp) return FALSE ; 
  key = *keyp ;
  if (!isNakedKey(key))  /* key is not a vocabulary entry */
    return TRUE ; /* no error when parsing same ace file twice */

  if  (myKey2Lex1->lock & ALIASSTATUS) /* follow alias line */
    key = *keyp = lexAliasOf (*keyp) ;

  if (!isNakedKey(key) ||  /* key is not a vocabulary entry */
      (myKey2Lex1->lock & (ALIASSTATUS | EMPTYSTATUS)))
	/* already aliased or EMPTY, to forbid complex alias graphs */
    return TRUE ; /* no error when parsing same ace file twice */

  if (lexiskeylocked(key))
    { if (doAsk) 
       messout("Sorry, alias fails because %s is locked elsewhere",
	      name(key)) ;
      return FALSE ;
    }

      /************* Trivial cases *************/
  h = handleCreate () ;
  newName = lexcleanup(newName, h) ;
  oldName = name(key) ;

  if (!strcmp(oldName, newName)) /* Identity */
    { handleDestroy (h) ; return TRUE ; }

  if (externalAlias) 
    oldnamepreserve = strnew(oldName, h) ;
  if (!lexstrIsCasecmp(oldName, newName)) /* Equivalence and certainly same length */
    { 
      lexHardRename(key, newName) ;         /* but not exact (dealt with earlier */
      lexKeyUpdated(key);
      goto ok ;
    }

     /************* Hard Rename to unknown name case *************/

  if (!keepOldName
      && !lexword2key (newName, &newKey, classe)) 
    { lexHardRename(key, newName) ;
      lexKeyUpdated(key);
      goto ok ;
    }

     /************ Else I want the 2 keys defined ************/

  lexaddkey(newName, &newKey, classe) ;  /* attention, this may invalidate oldName pointer */

  if (lexiskeylocked(newKey))
    { if (doAsk)
	messout("Sorry, alias fails because %s is locked elsewhere",
	      name(newKey)) ;
 
      handleDestroy (h) ;
      return FALSE ;
    }
  
  if(lexAliasLoop(key, newKey))  /* No loops */
    {
      if (doAsk)
	messout("Sorry, aliasing %s to %s would create a loop.",
		oldName, newName) ;
      handleDestroy (h) ;
      return FALSE ;
    }

  /* should I fuse ? */
  if (isNakedKey(key) == 1 )                /* key is actually empty */
    lexMakeAlias(key, newKey, keepOldName) ;
  else     /*  isNakedkey(key) == 2 */
    {
      if (pickType(key) == 'A' &&
	  isNakedKey(newKey) == 2 )  
	arrayKill (newKey) ;
      if (isNakedKey(newKey) == 1 )             /* key is full but newkey is empty */
	{                                  /* Just shift things over */
	  LEXP1
	    newLex1 = KEY2LEX1(newKey) ,
	    lex1 = KEY2LEX1(key) ;
	  LEXOFFST offset = lex1->nameoffset ;
	  KEY kk ;

	  lex1->nameoffset = newLex1->nameoffset ;
	  newLex1->nameoffset = offset ;
	  
	  kk = key ; key = newKey ; newKey = kk ; /* swap before clean up */
	  lexMakeAlias(key, newKey, keepOldName) ;
	  lexKeyUpdated(newKey);
	}
      else  /* also  isNakedKey(newKey) == 2 */
	{
	  if (pickType(key) != 'B')
	    { if (doAsk)
	      messout (" I don t know how to Fuse non TREE objects, sorry") ;
             handleDestroy (h) ;
	    return FALSE ;
	    }

	  /* NOTE WELL, if interactive we ask the user if they really want   */
	  /* to merge one object with another. But even if err_msg_out is    */
	  /* non-NULL we don't return anything here, the fusing is just done.*/
	  if (doAsk &&
	      !messQuery
	      (messprintf ("%s and %s both exist.\n"
			   "In case of conflict i will rather keep the data of %s.\n\n"
			   "Do you want to proceed",
			   name(key), name(newKey), name(newKey)))
	      )
	    {  handleDestroy (h) ; return FALSE ; }
	  bsFuseObjects (key, newKey) ;
	  /* to late to return FALSE, old is already destroyed */
	  lexMakeAlias(key, newKey, keepOldName) ;
	  lexKeyUpdated(newKey);
	  lexKeyUpdated(key);
	}
    }
ok:
  lexHashInsert(classe, key) ; /* re-insert in hash table */

  if (pickType (key) == 'B')  /* set subclass lex bits */
    {
      OBJ obj = bsCreate (key) ;
      bIndexObject(key, obj) ;   /* needed always */
      bsDestroy (obj) ;
    }
  lexAlphaMark (class(key)) ;

  if (externalAlias) 
    externalAlias (key, oldnamepreserve, keepOldName) ;
  messfree (oldnamepreserve) ; 
  handleDestroy (h) ;
  return TRUE ;
}
  
int lexNakedAlphaOrder(void *a, void *b) /* called in lexalpha.c */
{ 
  KEY k1 = *(KEY*)a, k2 = *(KEY*)b ;
  int i = class(k1) - class(k2) ;

  if (i)   
    return i ;
  return lexstrcmp (nakedName(k1), nakedName(k2)) ;
}

/********************************************/
/********************************************/

/* following are private communication lexword2key -> lexaddkey */
static KEY *KP ;
static BOOL isSubclassFailInWord2key ;

/*
* lexword2key finds existing keys, fails on non existing keys
* lexaddkey creates non existing keys
*/

BOOL lexword2key (const char *cp, KEY *key, KEY classe)
                        /* given a word *cp, sets its key in *key */
                        /* returns TRUE if found, FALSE if not */
                        /* No completion performed */
			/* classe can be a table or key in class Class */
{
  unsigned int h, dh = 0 ;
  int max, n , nn ;
  KEY *kp, *kpMax ;
  KEYSET hashTable ;
  int t ; unsigned char mask ;
  LEXI1 *lexi ;
  char *voc, *cq ;
  BOOL isCaseSensitive ;
  int (*lexstrIsCasecmp)(const char *, const char *) ;
    
  chrono ("lexword2key") ;
  pickIsA(&classe, &mask) ; 
  t = classe ;
  isCaseSensitive = pickList[t & 255 ].isCaseSensitive ;
  lexstrIsCasecmp = isCaseSensitive ? strcmp : strcasecmp ;

  if (t < 0 || t >= MAXTABLE)
    messcrash("lexword2key called on impossible class %d", classe) ;

  if(!lexIsRead[classe])
    lexReadTable(classe) ;

  if (Lexi1[classe] && !LexHashTable [classe])
    lexhread (classe) ;
  hashTable = LexHashTable [classe] ;
  if ( !cp || !*cp || !nClass[classe] || !arrayMax(Lexi1[classe]) )
    { *key = 0;    
      chronoReturn() ;
      return FALSE ;
    }

  chrono("lexword2key") ;
  cq = lexcleanup(cp, 0) ;  /* remove spaces etc */

  voc = stackText(Voc[classe],0) ;
  lexi = arrp(Lexi1[classe], 0 ,LEXI1) ;  
  h = hashText(cq, nn = nClass [classe], 0) ;

  n = max = arrayMax(hashTable) ;
  kp = arrp(hashTable, h , KEY) ;
  kpMax = arrp(hashTable, n -1, KEY) ;
  while ( n--)  /* circulate thru table, starting at h */
    { 
      if(!*kp)      /* found empty slot, cp is unknown */
	{
	  *key = 0 ;
	  KP = kp ;  /* private to be used by lexaddkey */
	  chronoReturn() ;
	  NWG++ ;
	  messfree(cq);
	  return FALSE ;
	}
      if (lexstrIsCasecmp(cq, voc +  (lexi + *kp)->nameoffset) == 0)
	{ 
	  /*
	   * the name we are looking for matches exactly
	   */
	  *key = KEYMAKE(classe,*kp) ;
	  chronoReturn() ;
	  if (mask != (mask & lexi[*kp].isMask)) /* subset criteria fail */
	    { isSubclassFailInWord2key = TRUE ;
	      NWG++ ;
	      messfree(cq);
	      return FALSE ;
	    }
	  NWF++ ;
	  messfree(cq);
	  return TRUE ;   /* found */
	}
        /* circulate */
      if (!dh)
	dh = hashText(cq, nn, TRUE) ;

      kp += dh ;
      if (kp > kpMax)
	    kp -= max ;

      NWC++ ;
    }
  chronoReturn() ;
  KP = 0 ;
  messfree(cq);
  return FALSE ;   /* Whole table scanned without success */
} /* lexword2key */

/**************************************************************************/
                 /* To obtain the next key in KEY  order        */
                 /* or the first entry if *key=0.               */
                 /* Updtates the value of *key,                 */
                 /* Returns 0 if the vocabulary is empty        */
                 /* or if *key is its last entry, 1 otherwise.  */
 
KEY lexLastKey (int classe)
{
  KEY key = 0, k = 0 ;
  
  while(lexNext(classe,&k))
    key = k ;

  return key ;
} /* lexLastKey */
 
/********/
		/* RD 980713: removed static i to check sequential use on same
		   class.  Necessary for AQL, but discussion with jean concluded
		   this is also correct, not changing the semantics described above!
		*/

BOOL lexNext (KEY classe, KEY *kp)
{
  LEXP1 lxp1 ;
  int i ;
  unsigned char mask ;
  
  if (!pickIsA(&classe, &mask))
    { *kp = 0 ; return FALSE ; }

  if (!lexIsRead[classe])
    lexReadTable(classe) ;
  
  if (!*kp)       /* Find beginning of lexique */
    {
      if (!classe)
	{ *kp = 1 ; return TRUE ; }    /* prevents looping on class 0 */
      else
	{
#ifndef NEW_MODELS
	  i = pickList[classe].type == 'B' ? 0 : -1 ;
#else
	  i = -1 ;
#endif
	}
    }
  else				/* *kp is already in the classe (we hope!) */
    { if (class(*kp) != classe)
	{ messerror ("class mixup in lexNext: class(%s) is %s not %s",
		     name(*kp), className(*kp), className(KEYMAKE(classe,0))) ;
	  return FALSE ;
	}
      i = KEYKEY(*kp) ;
    }
  
  if (Lexi1[classe]                              /* lexique non-empty */
      && i+1 < arrayMax(Lexi1[classe]))
    { lxp1 = arrp(Lexi1[classe], i + 1, LEXI1) - 1 ;
      while (++lxp1, ++i < arrayMax(Lexi1[classe]))    /* not at end of lexique */
	  if (!(lxp1-> lock & EMPTYSTATUS) &&
	      ( mask == (mask & lxp1-> isMask)))
	    { *kp = KEYMAKE(classe, i) ;
	      return TRUE;
	    } 
    }
  *kp = 0 ;
  return FALSE ;
}

/********************************************/
void lexkill(KEY kk)      /*deallocates a key*/
{
  iskey(kk) ;
  
  if (!myKey2Lex1)
    return ;
  chrono("lexkill") ;
  
  /*  saveAll() ;  must be done before  to empty the cache 
      but not here because lexkill is used by sessionControl 
      */ 

  if (lexDestroyAlias(kk))
    return ;

  if(myKey2Lex1->dlx.dk)
   { myKey2Lex1->dlx.dk = 0 ;
     leximodif[class(kk)] = TRUE ;
     lexAlphaMark (class(kk)) ;
   }
  
  if ( ! (myKey2Lex1->lock & EMPTYSTATUS ))
    { myKey2Lex1->lock |= EMPTYSTATUS;
      leximodif[class(kk)] = TRUE ;
      lexAlphaMark (class(kk)) ;
    }
  myKey2Lex1->is2 = 0 ;
  if (myKey2Lex2)
    lexi2clear(myKey2Lex2) ;
  
  chronoReturn() ;
  return;
}

/**************************************************************/

KEY str2tag (const char *tagname)
{
  KEY k; 
 
  lexaddkey (tagname, &k, 0);
  return k;
}

/**************************************************************/
/* strnew */
BOOL lexaddkey (const char *cp, KEY *kptr, KEY t)
	/*
	* cp is an object name.  t is a class number.  We find
	* that name in that class.  *kptr is filled in with the
	* key of the object with that name.
	*
	* Returns TRUE if we needed to create the object, or
	* FALSE if the object already existed.
	*/
                   /*add to the t lexique the word *cp */
                   /*returns TRUE if added, FALSE if known */
{
 KEY k;
 LEXI1 ai;
 unsigned char mask ;
 char *cq ;
 BOOL debug = FALSE ;
/* extern BOOL READING_MODELS ;*/

 chrono("lexaddkey") ;
 if (debug) fprintf (stderr, "lexaddkey %s\n", cp ? cp : "NULL") ;
 pickIsA(&t, &mask) ;
 if (debug && cp && !strcmp (cp, "Greek"))
   exit (0) ;
#ifdef NEW_MODELS
 if (t == _VMainClasses && !READING_MODELS)
   {
     invokeDebugger() ;
     return cp ? lexword2key (cp, kptr, t) : 0 ;
   }
#endif


 KP = 0 ;
 isSubclassFailInWord2key = FALSE ;
 if (cp && (lexword2key (cp, kptr, t) || isSubclassFailInWord2key))
   { 
     /*
     * We found that the object name exists already.  Look to see if the object
     * wants the name updated when you edit it - if so, we write over the name
     * in the hash table so that it reflects the new case of the object name
     * that was entered by the user.
     */
     LEXI1 *lexi ;
     char *voc ;
     chronoReturn() ;
     t = class(*kptr);
     if (pickList[t & 255 ].updateNames)
       {
	 char *cp1;
         voc = stackText(Voc[t],0) ;
         lexi = arrp(Lexi1[t], 0 ,LEXI1) ;
	 cq = lexcleanup(cp, 0);
	 cp1 = voc +  (lexi + KEYKEY(*kptr))->nameoffset;
         if ( (strcmp(cp1, cq) != 0) && (strlen(cp1) == strlen(cq)) )
           {
	     vocmodif[t] = TRUE ; /* mark the buffer dirty */
	     strcpy(cp1, cq);     /* change the name en place */
           }
	 messfree (cq) ;
       }
     return FALSE ;
   }
				/* word already known */
 if (!lexIsRead[t])
   {
     lexReadTable(t) ;
     KP = 0 ;
   }
                                /* initialise the lexique */
  if (!Lexi1[t])
    {
      int nBits = 8 ;
      Voc[t] = stackCreate(1 << (nBits + 2)) ;
      Lexi1[t] = arrayCreate(1 << (nBits - 1), LEXI1);
      
      timeStamps[t] = arrayCreate(1 << (nBits - 1), TIMESTAMP);
      timeModif[t] = TRUE;
	
      LexHashTable[t] = arrayCreate(1<<nBits, KEY);
      keySet(LexHashTable[t], (1<<nBits) - 1) = 0 ; /* make room */
      nClass[t] = nBits ;
#ifndef NEW_MODELS
    if (t && (pickType(KEYMAKE(t,0)) == 'B'))
      {
	char buff2[1000] ; /* recursive call, we need a local buffer */
	sprintf (buff2, "?%s", pickClass2Word(t)) ;
	lexaddkey (buff2, &k, t) ;
      }
#endif
     if (lexaddkey("\177\176(NULL KEY)",&k,t) &&
	 KEY2LEX1(k) )
       KEY2LEX1(k)-> lock = EMPTYSTATUS ;  /* i.e. never display it */
     KP = 0 ;
   }
 
 k = (KEY) arrayMax(Lexi1[t]) ;
 if (!cp && k <= 2) /* call with cp=0 used for initialisation */
   {
     *kptr = KEYMAKE(t,0);
     chronoReturn() ;
     return TRUE ;
   }
 
 cq = cp ? lexcleanup (cp, 0) : 0 ;
 if (!cp || !*cp || !cq || !*cq)
   messcrash("Attempt to create a key with emtpy name in class%s",
	     className (KEYMAKE(t,0))) ;
 
 if ( (k>=MAXLEX) ||
      ( (stackMark(Voc[t]) + strlen(cq)) > MAXVOCAB ))
   messcrash("Lexique %d full: now %d keys, %d size stack",
	     t, k, stackMark(Voc[t])) ;
 
      /*create a lex record*/
 ai.nameoffset = stackMark(Voc[t]) ;
 ai.dlx.dk = 0 ;
 ai.lock = 0;
 ai.isMask = 0;
 ai.is2 = 0 ;

      /* add the word to the end of vocabulary */
 pushText(Voc[t], cq) ;
         /*add an entry to the end of lex */
 array(Lexi1[t], (int) k, LEXI1) = ai ; 
/* note that Voc and Lexi must be updated first 
   for lexHashInsert to work */

 if (KP)
   { chrono("Found KP") ; chronoReturn() ;
     *KP = k ;
     /* direct insertion , be careful if you touch the code */
     /* KP is inherited from the last call to lexword2key */

     /* SRK - restored this to original in lexsubs.c, we were
        re-hashing on every call to lexaddkey 
	if(arrayMax(Lexi1[t]) > .34 * keySetMax(LexHashTable[t]))
	*/
     if(arrayMax(Lexi1[t]) * HASH_RATE > keySetMax(LexHashTable[t]))
       lexReHashClass(t) ;
   }
 else
   lexHashInsert(t, k) ;
 
 leximodif[t] = lexhmodif[t] = vocmodif[t] = TRUE ;
 lexAlphaMark (t) ;
  *kptr = KEYMAKE(t,k) ;
  chronoReturn() ;
  messfree (cq) ;
  return TRUE;
}  /* lexaddkey */

/*************************************************/

int lexMax(int t)
{
  int nn = 0 ;
  if (t<0 || t>MAXTABLE)
    return 0 ;
  if (!lexIsRead[t])
    lexReadTable(t) ;
  if (Lexi1[t])
    nn = arrayMax(Lexi1[t]) ; /* nn introduced to please solaris */
  return nn ;
}

/*************************************************/

int lexHashMax(int t)
{
  return LexHashTable[t] ? arrayMax (LexHashTable[t]) : 0 ;
}

/*************************************************/

int vocMax(int t)
{
  return Voc[t] ? stackMark(Voc[t]) : 0 ;
}

/*************************************************/
/*************************************************/

void lexSessionStart (void)
{
  int t ;

  for (t = 0 ; t < MAXTABLE ; ++t)
    if (lexIsRead[t] && Lexi1[t])
      { lexSessionStartSize[t] = arrayMax(Lexi1[t]) ;
	lexClearClassStatus (t, TOUCHSTATUS) ;
      }
    else
      lexSessionStartSize[t] = 0 ; 
  keySetDestroy (touchedByClient) ;
  if (externalSaver) touchedByClient = keySetCreate () ;
}

void lexSessionEnd (void)
{
  char lexBuffer [ DIR_BUFFER_SIZE ] ;
  KEYSET newKeys = 0, touchedKeys = 0 ;
  int	 t, i ;
  LEXI1* q ;
  char *cp = 0 ;
  FILE *f ;
  static int lastSession = 0 ;
  AC_HANDLE h = 0 ;
#ifdef NEW_MODELS
#define KEYKEY0 1
#else
#define KEYKEY0 2
#endif

  if (lastSession == thisSession.session)
    return  ;			/* prevents loops */

  newKeys = keySetCreate () ;
  touchedKeys = keySetCreate () ;
  for (t = 0 ; t < MAXTABLE ; t++)
    if (!pickList[t].protected && lexIsRead[t] && Lexi1[t])
      { 
	q = 0 ;
	i = lexSessionStartSize[t]>KEYKEY0 ? lexSessionStartSize[t] : KEYKEY0 ;
	if (i < arrayMax(Lexi1[t]))
	  {  
	    q = arrp(Lexi1[t], i, LEXI1) ;
	    for(; i < arrayMax(Lexi1[t]) ; ++i, ++q)
	      if (!(q->lock & ALIASSTATUS))
		keySet(newKeys, arrayMax(newKeys)) = KEYMAKE(t,i) ;
	  }
	if (KEYKEY0 < arrayMax(Lexi1[t]))
	  for (i = KEYKEY0, q = arrp(Lexi1[t], i, LEXI1) ; 
	       i < arrayMax(Lexi1[t]) ; ++i, ++q)
	    if ((q->lock & TOUCHSTATUS) && !(q->lock & ALIASSTATUS))
	      keySet(touchedKeys, arrayMax(touchedKeys)) = KEYMAKE(t,i) ;
      }

  if (keySetMax (newKeys))
    { keySetSort (newKeys) ;
#ifdef NEW_INSIDE
    sprintf (lexBuffer, "new-%s", name(thisSession.userKey)) ;
    lexaddkey (lexBuffer, &key, _VKeySet) ;
    arrayStore (key, newKeys, "k") ;
#else
    sprintf (lexBuffer, "database/new/new-%s", name(thisSession.userKey)) ;
    cp = sessionFilName(lexBuffer, 0,0) ;
#if !( defined(MACINTOSH) || defined(WIN32) )
    seteuid(euid);
#endif
    f = cp ? fopen (cp, "w") : 0 ;
    if (!f) 
      { 
	sprintf (lexBuffer, "mkdir -p %s", sessionFilName("database/new", 0,0)) ;
	callSystem (lexBuffer) ;
	sprintf (lexBuffer, "database/new/new-%s", name(thisSession.userKey)) ;
	cp = sessionFilName(lexBuffer, 0, 0) ;
	f = cp ? fopen (cp, "w") : 0 ;
	}
      if (!f) 
	{ 
	  messerror("Please create a folder called 'database/new'") ;
	  sprintf (lexBuffer, "new-%s", name(thisSession.userKey)) ;
	  cp = sessionFilName(lexBuffer, 0, 0) ;
	  f = cp ? fopen (cp, "w") : 0 ;
	}
      if (f)
	{ 
	  keySetDump(f,0,newKeys) ;
	  fclose(f) ;
	} 
#if !( defined(MACINTOSH) || defined(WIN32) )
      seteuid(ruid);
#endif
#endif
    }
  keySetDestroy (newKeys) ;

  if (keySetMax (touchedKeys))
    {
      keySetSort (touchedKeys) ;
#ifdef NEW_INSIDE
      sprintf (lexBuffer, "touched-%s", name(thisSession.userKey))
      lexaddkey (lexBuffer, &key, _VKeySet) ;
      arrayStore (key, touchedKeys, "k") ;
#else  /* !NEW_INSIDE */
#if !( defined(MACINTOSH) || defined(WIN32) )
      seteuid(euid);
#endif
       f = cp ? fopen (cp, "w") : 0 ;
      if (!f) 
	{ 
	  if (!sessionFilName("database/touched", 0,"r"))
	    {
	      sprintf (lexBuffer, "mkdir -p %s", sessionFilName("database/touched", 0,0)) ;
	      callSystem (lexBuffer) ;
	    }
	  sprintf (lexBuffer, "database/touched/touched-%s", name(thisSession.userKey)) ;
	  cp = sessionFilName(lexBuffer, 0,0) ;
	   f = cp ? fopen (cp, "w") : 0 ;
	}
      if (!f) 
	{ 
	  messerror("Please create a folder called 'database/touched'") ;
	  sprintf (lexBuffer, "database/touched-%s", name(thisSession.userKey)) ;
	  cp = sessionFilName(lexBuffer, 0,0) ;
	   f = cp ? fopen (cp, "w") : 0 ;
	}
      if (f)
	{ 
	  keySetDump(f,0,touchedKeys) ;
	  fclose(f) ;
	}
#if !( defined(MACINTOSH) || defined(WIN32) )
      seteuid(ruid);
#endif /* !(MAC || WIN) */

#endif /* !NEW_INSIDE */
    }
  keySetDestroy (touchedKeys) ; 
  keySetDestroy (touchedByClient) ;
  handleDestroy(h);

  lastSession = thisSession.session ;

  return;
} /* lexSessionEnd */

void lexShutDown (void)
{
  int t = 256 ;
  while (t--)
    {
      arrayDestroy (Lexi1[t]) ;
      keySetDestroy (LexHashTable[t]) ;
      stackDestroy (Voc[t]) ;
      arrayDestroy (Lexi2) ;
      arrayDestroy(timeStamps[t]);
    }
}

/************** end of file **********************/

#endif /* ACEDB4 */



void lexSpecialFix (void)
{
  int t ;
  KEY key ;
  int nn, nMax ;
  char *cp, *cq ;
  Stack V ;
  int touched = 0, total = 0 ;
  LEXP1 lx ;

  t = 256 ;
  while (t--)
    {
      if (t <= _VMainClasses)
	continue ;  /* beware with kernel */

      if (t == _VText)
	continue ;
      if (pickList[t].protected)
	continue ;

      nMax = lexMax (t) ; /* will load the lex and the Voc */
      if (!nMax)
	continue ;
      key = t << 24 ;

      V = Voc[t] ;
      if (!stackExists (V))
	messcrash ("Null Voc[%d]", t) ;
      for (nn = touched = 0, lx = arrp(Lexi1[t], 0, LEXI1); nn < nMax ; nn++, lx++)
	{
	  cp = stackText(V, lx->nameoffset) ;
	  while (*cp++) /* cannot be the full name */
	    if (*cp == 'M')
	      {
		cq = cp ;
		if (*cq++ == 'M' &&
		    *cq++ == 'a' &&
		    *cq++ == 'r' &&
		    *cq++ == '0' &&
		    *cq++ == '4' &&
		    *cq == 0
		    )
		  {
		    cq = cp ;
		    *cq++ = 'A' ;
		    *cq++ = 'u' ;
		    *cq++ = 'g' ;
		    touched++ ;
		  }
	      }
	}
      if (touched)
	{
	  total += touched ;
	  vocmodif [t] = TRUE ;
	  printf ("Modified %10d entries in %s\n", touched, className (key)) ;
	}
    }
  printf ("total    %10d\n\n", total) ;
}
