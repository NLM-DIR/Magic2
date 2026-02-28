/*  File: adisk.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1996
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: @(#)adisk.c	1.6 3/21/97 
 * Description: stores arrays on disk
 		replacement for disknew.c, blocksubs.c
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 23 13:01 1998 (fw)
 * * Sep 15 09:20 1998 (edgrif): Add macro to define bitField from 
 *              bitset.h.
 * Created: Tue Dec  3 21:07:15 1996 (rd)
 *-------------------------------------------------------------------
 */


#define DEF_BP  /* avoid typedefing it as void* in disk.h */
typedef struct block *BLOCKP ;
typedef BLOCKP BP ;

#include "acedb.h"
#include "adisk.h"
#include "byteswap.h"
#include "bitset.h"
void lexSetDisk(KEY key, DISK dk) ; 
  /* each file starts with a superblock of fixed size
   * addressed in a particular way.
   * I choose 2048, which is rather large, but enhances, i hope, 
   * the proba of later exact page access 
   */

#define SUPERKEY __superKey
#define SUPERBLOCKSIZE  ((myoff_t)2048)
#define MIN_PARTITION_SIZE 64
#define SIZE_SHIFT_FACTOR 1 

static DISK globalAddress = 0 ;
static int 
  globalFileSize = 0, globalFileMax = 0, globalFileNumber = 0,
  globalBatFileSize = 0 ;
static BOOL debug = FALSE ;
static BOOL debug1 = FALSE ;
static BOOL debug2 = FALSE ;


typedef struct block /* sizeof(block) must be < BLOC_SIZE */
  {                     /* ATTENTION, used in read/write superblock */
    DISK gAddress ;
    int mainRelease ,  subDataRelease, subCodeRelease , subSubCodeRelease ;
    int session, byteSexIndicator;
    char dbName[32] ; /* added as of release 4.3 */
  } BLOCK;   /* super block info */


static int currentSession = 0 ;

typedef struct {
  int size ;            /* various files have various blocking size */
  int max ;             /* number of allocated blocs in that file */
  int number ;          /* 1+index of this file in array(partitions) */
  int position ;        /* position of target in block units */
  BitSet bat ;          /* bloc allocation table of this file */
  BitSet batPlus ;      /* blocs allocated during this session */
  BitSet batMinus ;     /* blocs liberated during this session */
  BOOL full ;           /* NOT YET TRATED */
  int wlastpos ;      /* to cluster writing */
  myoff_t rcurr, wcurr ;        /* current position of file pointer */
  int r, w ;            /* file descriptors, -1 if closed */
} PARTITION ;

typedef struct {	/* NB keep this a multiple of 8 bytes,*/
  int flags , dummy ;
  DISK disk ;
  KEY key, parent ;     /* edit swapHeader if you change this structure */
  KEY session ;
  int n1, s1, n2, s2 ;
} DISK_HEADER ;
 
int MAXKNOWNCACHE =  2000 ;    /* default number of cacheEntries in the cache, used by objcache.c */
/* should be configurable, that was in blocksubs */


static char *pName (PARTITION *p) ;
static PARTITION *newPartition (int targetsize) ;
static void aDiskFile (PARTITION *p,BOOL doCreate) ;
static PARTITION *aDiskLocate (DISK d) ;

static void aDiskBATwrite(void) ;
static BOOL batIsModified = FALSE ;
static KEY currentKey = 0 ;
static Array partitions = 0 ;
static int nRead = 0, nWrite = 0, nBytesRead = 0, nBytesWritten = 0 ;

static BOOL ps_debug = FALSE ;
static void  aDiskReadSuperBlock (BP bp) ;
static Array aDiskBatDecompress (Array caption, BitSet bat) ;

static BOOL aDiskBigArrayGet (KEY key, KEY *pp, Array *ap, BigArray *bp)  ;
static void aDiskBigArrayStore (KEY key, KEY parent, Array a, BigArray b) ;

/******* useless  missing stuff *******/
void blockInit(void) {;}
void blockavail(int *used,int *pinned,int *free,int *modif)
{ *used = *pinned = *free = *modif = 0 ; }
int blockMax (void) { return 0 ; }

void bsTreeKill  (KEY key)
{  aDiskKill (key) ; }

/************************************************************/
/**************** utilities  ********************************/

static void diskcrash (PARTITION *p, char *s)
{ 
  messcrash ("Disk access error - %s: file %s, key %s",
	     s, pName(p), name(currentKey)) ;
}

static char *pName (PARTITION *p)
{
  static char buffer[512] ;

  if (!p)
    return "No current partition" ;

  sprintf(buffer,"%s", sessionFilName
	  (messprintf("database/a_%d_%d",p->size,p->number + 1),"ace5",0)) ;
  return buffer ;
} /* pName  */

/********************************************/

static void aSwapHeader (DISK_HEADER *h)
{
  h->disk = swapKEY(h->disk) ;
  h->key = swapKEY(h->key) ;
  h->parent = swapKEY(h->parent) ;
  h->session = swapKEY(h->session) ;
  h->n1 = swapInt(h->n1) ;
  h->n2 = swapInt(h->n2) ;
  h->s1 = swapInt(h->s1) ;
  h->s2 = swapInt(h->s2) ;
  h->flags = swapInt(h->flags) ;
  h->dummy = swapInt(h->dummy) ;
} /* aSwapHeader */

/********************************************/

                  /*returns the number of used and free blocks*/
void diskavail(unsigned long int *free, unsigned long int *used, unsigned long int *plus, unsigned long int *minus, unsigned long int *dr, unsigned long int* dw)
{
  int i ;
 unsigned long int n ;
 long bused = 0, balloc = 0 ; 
 BOOL debug = *dw == 123 ? TRUE : FALSE ;
 PARTITION *p ;

  if (!partitions) 
    return ;

  *dr = nRead ;*dw = nWrite ;

  *used = *plus = *minus = *free = 0 ; 
  if (debug) freeOutf("//\n") ;
  for (i = 0 ; i <  arrayMax(partitions) ; i++)
    { 
      p = arrp(partitions, i, PARTITION) ;
      n = bitSetCount (p->bat) ;
      *used +=  n ;
      *plus += bitSetCount (p->batPlus) ;
      *minus += bitSetCount (p->batMinus) ;
      *free += p->max - n ;
      bused += n * (long)p->size ;
      balloc += p->max * (long)p->size ; 
      if (debug) freeOutf("// partition %4d  size %9d \tusing %9ld/%ld blocks\t  %9ld/%ld kb\n",
			i+1, p->size, n, p->max, 
			(n * p->size) / 1024, 
		        (p->max * p->size) / 1024) ;

    }
  if (debug) freeOutf("// Total diskspace  %ld/%ld kb \n", bused/1024, balloc/1024) ;
} /* diskavail */

/********************************************/
/* Count the number of Set bits */
static int aDiskBatSet (Array pp, int type)
{
  int i, n = 0 ; 
  PARTITION *p ;
  BitSet bb = 0 ;

  i = arrayMax(pp) ;
  while (i--)
    { 
      bb = 0 ;
      p = arrp(pp, i, PARTITION) ;
      switch(type)
	{
	case 1: bb = p->bat ; break ;
	case 2: bb = p->batPlus ; break ;
	case 3: bb = p->batMinus ; break ;
	}
      if (bb)
	n += bitSetCount (bb) ;
    }
  return n ;
} /* aDiskBatSet */

/************************************************************/
/***************** Initialisation ********* *****************/
/*************************************************************/
/* import cache size or used defaults */
static FREEOPT cacheSizeOptions[] =
{
  {3, "Cache"},
  {'1', "CACHE1"},
  {'2', "CACHE2"},
  {'3', "DISK"},
  } ;
int  MAXKNOWNACACHE = 1<<11 ; /* 2 megs */
static void getCacheSize()
{
  int n  ; KEY option ;
  char *cp = sessionFilName ("wspec/cachesize","wrm","r") ;
  FILE *fil = cp ? filopen(cp, 0, "r") : 0 ;
    
  if(!fil)
    return ;

  freespecial ("\n\t\"\\") ;
  while(freeread(fil))
    if (freekey(&option, cacheSizeOptions))
      switch (option) 
	{
	  case '1':
	  freenext() ; freestep ('='); freenext() ;
	  if (freeint(&n))
	    MAXKNOWNACACHE = n << 10 ;
	  else
	    messerror("In wspec/cachesize.wrm, I cannot read the number after CACHE1") ;
	  break ;
	  case '2':
	  freenext() ; freestep ('='); freenext() ;
	  if (freeint(&n))
	    MAXKNOWNCACHE = n ;
	  else
	    messerror("In wspec/cachesize, I cannot read the number after CACHE2") ;
	  break ;
	  case '3':
	    break ;
	}
  filclose(fil) ;
} /* getCacheSize */

/************************************************************/
/* verify global compatibility */
static void aDiskCheckInit (void)
{
  if (sizeof(DISK_HEADER) % 8)
    messcrash ("DISK_HEADER size not a multiple of 8") ;
  if (sizeof(DISK) != 4)
    messcrash ("DISK size differs from 4 bytes") ;
  getCacheSize() ;

#if defined(MACINTOSH)
  messcrash ("mac version not yet done") ;
#endif
} /* aDiskCheckInit */

/************************************************************/

BOOL dataBaseCreate (void)
{
  swapData = FALSE ;   /* The database file must be in our order */
  seteuid (euid);  /* switch these in the public functions */

  aDiskCheckInit() ;
  system(messprintf("\\rm -f %s/*.ace5",filName("database",0,0))) ;
  partitions = arrayCreate (8, PARTITION) ;
  newPartition (MIN_PARTITION_SIZE) ; /* create a trivial partition */
  messdump("Data Base Creation\n") ;
  seteuid (ruid);  
  return TRUE ;
} /* dataBaseCreate */

/************************************************************/

void aDiskAccessInit(BP bp)
{ 
  PARTITION *p = 0 ;

  aDiskCheckInit() ;
  if (debug1) messdump("Data Base Access Init\n") ; 
  partitions = arrayCreate (8, PARTITION) ; 
       /* fake partition to read the superblock */
  p = arrayp (partitions, 0, PARTITION) ;
  p->size =  MIN_PARTITION_SIZE ; 
  p->number = 0 ;
  p->max = 2 ;
  p->r = p->w = -1 ;
  aDiskReadSuperBlock (bp) ; /* do read super block */
  close (p->r) ;
     
       /* fake partition to read the global lex */
  lexSetDisk (__Global_Bat, globalAddress) ;
  p = arrayp (partitions, globalFileNumber, PARTITION) ;
  p->size =  globalFileSize ;
  p->number = globalFileNumber ;
  p->max = globalFileMax ;
  p->r = p->w = -1 ;
} /* aDiskAccessInit */

/************************************************************/

void aDiskClose (void)
{
  int i ; PARTITION *p ;
#if defined(MACINTOSH)
  messcrash ("mac version not yet done") ;
#endif
  if (!partitions) return ;

  seteuid (euid);  
  i = arrayMax(partitions) ;
  while (i--)
    { p = arrp(partitions, i, PARTITION) ;
      if (p->r > -1) close(p->r) ; 
      p->r = -1 ;
      if (p->w > -1) close(p->w) ; 
      p->w = -1 ;
    }
  seteuid (ruid);  
} /* aDiskClose */

/************************************************************/
/***************** Bloc allocation handling *****************/
/************************************************************/

static PARTITION *newPartition (int targetsize)
{
  PARTITION *p = arrayp (partitions, arrayMax(partitions), PARTITION) ;
  int max, size ;

  /* The size should be bigger than target size
     We construct a different file from MIN_PARTITION_SIZE bytes onwards 
      which seems good in the testing phase
     As a test, i allocate minimum 1 Mbyte or 4 times size
     */
  max = 1 << 20 ; /* one mega */
  size = MIN_PARTITION_SIZE ; 
  while (size < targetsize) size <<= SIZE_SHIFT_FACTOR ; /* duplicated in aDiskAssign */
  max /= size ;
  if (max < 8) max = 8 ;
  max = max/8 ; max = 8 *max ;  /* must be divisible by 8 for BAT to work */
  p->size =  size ;
  p->max =  max ;        
  p->full = FALSE ;
  p->wlastpos = 0 ; p->rcurr = p->wcurr = 0 ;
  p->number = arrayMax(partitions) - 1 ;

  batIsModified = TRUE ;
  p->bat = bitSetCreate (max, 0) ;   
  p->batPlus = bitSetCreate (max, 0) ;   
  p->batMinus = bitSetCreate (max, 0) ;   
  bitUnSet (p->bat, max - 1) ; /* make room */
  bitUnSet (p->batPlus, max - 1) ; /* make room */
  bitUnSet (p->batMinus, max - 1) ; /* make room */
  if (p->number == 0) bitSet (p->bat, 0) ; /* address 0:0 is not usable */

  aDiskFile (p, TRUE) ;
  messdump ("Created New partition %s\n", pName(p));
  return p ;
} /* newPartition */

/************************************************************/

static PARTITION *aDiskExtend (PARTITION *p) 
{
  if (p->max * p->size > (1 << 25)) /* with 1<<25, the largest data file is 64 Mb */
    return newPartition (p->size) ; 
  aDiskFile (p, FALSE) ;      
  batIsModified = TRUE ;
  bitExtend (p->bat, p->max) ;
  bitExtend (p->batPlus, p->max) ;
  bitExtend (p->batMinus, p->max) ;
  bitUnSet (p->bat, p->max - 1) ; /* make room */  
  bitUnSet (p->batPlus, p->max - 1) ; /* make room */
  bitUnSet (p->batMinus, p->max - 1) ; /* make room */
  return p ;
} /* aDiskExtend */

/************************************************************/
static DISK lastDiskLocate = 0 ;	/* info for two-phase diskRead*() */
static PARTITION *lastP = 0 ;
static BOOL isHeaderRead = FALSE ;
/*************************************************/
/* do free only if allocated during this session */
/*************************************************/
/* utility, could be in bitset.c if that existed */

static BOOL firstUnsetBit (PARTITION *p, int *ip)
{
  BitSet bb = p->bat ;
  register unsigned int uu, *uip ;
  register int i ; 
  unsigned long int ii, first, last ;

  if (! bitSetExists (bb) || ! bigArrayMax(bb) || bb->size != 32)
    return FALSE ;
  
  first = p->wlastpos >> 5 ; last = bigArrayMax(bb) ;
lao:
  uip = bigArrp(bb, first, unsigned int) - 1 ; ii = first - 1 ;
  i = last ;
  while (uip++, ii++, i--)
    if (~(*uip)) break ;
  if (i == -1)
    {
      if (first) 
	{ last = first ; first = 0 ; goto lao ; }
      else
	return FALSE ;
    }

  uu = *uip ;
  for (i = 0 ; i < 32  && (bitField[i] & uu) ; i++) ;
  /* note that 32 == 1 << 5 */
  if (i == 32)  messcrash ("bit manipulation error in firstUnsetBit") ;

  i = (ii << 5) | i ;
  if (i >= p->max) return FALSE ;  /* p->max may not be divisible by 32 */
  *ip = i ;
  return TRUE ;
} /* firstUnsetBit */

/*************************************************/
/* give a new empty adress adapted to size */
static PARTITION *aDiskAlloc (DISK *dp, int n1, int n2)
{
  int j1, j, jx, size, targetsize, pos = 0 ;
  PARTITION *p = 0 ;

  p = arrp (partitions, 0, PARTITION) - 1 ;

  /* choose the correct partition size */
  targetsize = sizeof(DISK_HEADER) + n1 + n2 ;
  size = MIN_PARTITION_SIZE ; 
  while (size < targetsize) size <<= SIZE_SHIFT_FACTOR ; /* duplicated in aDiskAssign */
  /* find a partition of the correct size with an empty block
   * or the smallest partion of the correct size
   */
  for (jx = j1 = j = 0, p = arrp (partitions, j, PARTITION); 
       j < arrayMax(partitions) ; p++, j++)
    {
      if(size != p->size)         /* wrong size */
	continue ;
      if (!jx || jx > p->max)     /* smallest of correct size */
	{ j1 = j ; jx = p->max ;} ;
      if (firstUnsetBit (p, &pos))
	goto ok ;
    }
  if (!jx)
    p = newPartition(size) ; 
  else
    { p = arrp(partitions, j1, PARTITION); p = aDiskExtend(p) ; }
  if (!firstUnsetBit (p, &pos))
    diskcrash (p, "Cant find free block in a new or extended partition") ;
ok:
  p->position = pos ;
  if (p->number >= (1 << 12) ||  p->position >= (1 << 20))
    diskcrash (p, "The kernel cannot memorize such  large disk address") ;
  *dp = ((p->number & 0xfff) << 20) | (p->position & 0xfffff) ;  /* see below the & 0xfffff clause */
  batIsModified = TRUE ;
  bitSet (p->bat, p->position) ;
  bitSet (p->batPlus, p->position) ;
  return p ;
} /* aDiskAlloc */

/************************************************************/
/*
static BOOL mybit(BitSet b, int pos)  // for debugging
{ return bit(b,pos) ; }
*/

/************************************************************/

static void aDiskFree(PARTITION *p)
{
  unsigned long int j = p->max, pos = p->position ; 
  if (pos >= j)
    messcrash("aDiskfree called with a wrong disk address");

  if (!bit (p->bat, pos))
    messcrash ("aDiskFree called on free block: %s", name(currentKey));

  /* on peut oublier les blocs crees pendant la session */
  if (bit (p->batPlus, pos))
    { bitUnSet (p->batPlus, pos) ;     /* free it */
      bitUnSet (p->bat, pos) ;     /* free it */
    }
  else
    bitSet (p->batMinus,pos) ;
  batIsModified = TRUE ;
} /*  aDiskFree */

/********************************************/
   /* The strategy is always to move blocks on disk when rewriting
      and never reuse a block during a given session.
      In this way, the former session certainly remains consistent
      
      There is however a very important exception to this rule,
      The session objects and the Bat arrays are not moved. The reason 
      is that during session 8 you may be fusing session 3 and 4,
      then you modify obj session4 and its bats. By doing it en place
      you insure that it is written in a block allocated during session 4.
      Otherwise, if you eventually destroy session 8 before session 4, the
      bat of session 8 would be lost.
      */

/************************************************************/

static PARTITION *aDiskDoAssign (DISK *dp, KEY key, int n1, int n2)
{
  DISK d = 0 ;
  PARTITION *p = 0 ;
  int size, targetsize ;
  
  if (key == SUPERKEY)
    { *dp = d = 1 ; return aDiskLocate(d) ; }
  
  d = lexDisk(key) ;
  targetsize = sizeof(DISK_HEADER) + n1 + n2 ;
  size = MIN_PARTITION_SIZE ; while (size < targetsize) size <<= SIZE_SHIFT_FACTOR ; /* duplicated in newPartition */
  if (d)
    p = aDiskLocate(d) ;
  if (p)
    { 
      *dp = d ;
      /* do not move these very special classes */
      if (key == __Global_Bat  && debug1) 
	freeOutf ("aDiskDoAssign _Global_BAT address %d\n", d) ;
      if (key == __voc3 && debug1)
	messdump("aDiskDoAssign voc3 address %d\n", *dp) ;
      if (class(key) == _VBat ||
	  class(key) == _VSession)
	{
	  if (size < 1000) size = 1000 ;
	  if(size == p->size) /* we use to crash on "Disk assign cannot reassign a Bat or Session object" */
	    return p ;
	}
      if (bit(p->batPlus, p->position) &&  /* allocated during this session */
	  size == p->size)                 /* still correct size */
	return p ;
      if (class(key) != _VGlobal) /* Freeing the global class in unconsistent */
	aDiskFree(p) ;
    }
  /* Else, we move the block on disk */
  p = aDiskAlloc(dp, n1, n2) ;
  lexSetDisk(key,*dp) ;
  if (key == __Global_Bat && debug1) 
    freeOutf ("aDiskDoAssign2 _Global_BAT address old=%d new=%d\n", d,*dp) ;
  if (key == __voc3 && debug1) 
    freeOutf ("aDiskDoAssign2 voc3 address %d\n", *dp) ;
  
  return p ;
} /* aDiskDoAssign */

/************************************************************/

static PARTITION *aDiskLocate (DISK d)
{ 
  PARTITION *p = 0 ;
  int i = d >> 20 ; /* same code, no check, in disk prepare */

  if (i >= 0 && i < arrayMax(partitions))
    p = arrp (partitions, i, PARTITION) ;
  else
    messcrash ("aDiskLocate received a bad disk address %x -> partition %d",
			d, i) ;
  p->position = d & 0xfffff ;           /* see above  the << 20 clause */
  if (p->position >= p->max)
    messcrash ("aDiskLocate received a bad disk address %x >= max(partition %d) = %d",
	       d, i, p->max) ;
  return p ;
} /* aDiskLocate */

/************************************************************/
/************************************************************/
/************************************************************/
/**************** low level read/write ****************/

static void aDiskFile (PARTITION *p, BOOL doCreate)
{
  FILE *f ;
  unsigned long int i = p->max ;
  unsigned int sbsize = SUPERBLOCKSIZE, umax = p->size, un = 1 ; /* pedantic prototyping */
  void *emptybp ; 
  char timeBuf[25] ;

  while (umax < (1 << 20) && (i > 1)) { umax <<= 1 ; i >>= 1 ;}
  emptybp = messalloc(umax) ;       /* zeroing implied */

  if (doCreate)
    {
      f =  fopen(pName(p), "wb") ;
      if (!f)
	diskcrash (p, "aDiskFileCreate cannot create new partition") ;
    }
  else
    {
      if (p->w >= 0) close (p->w) ; 
      if (p->r >= 0) close (p->r) ;
      f =  fopen(pName(p), "ab") ;
      if (!f)
	diskcrash (p, "aDiskFile cannot create new partition") ;
      p->max <<= 2 ; 
      i = 3*i ;  /*we multiply the size of the file by 4 */
    }
  

  /* first a super block */
  if (fwrite (emptybp, sbsize, un,  f) != un)
    diskcrash (p, "Cannot actually write the new disk file");

  while (i--)
    if (fwrite (emptybp, umax, un,  f) != un)
      diskcrash (p, "Cannot actually write the new disk file");
  
  fclose (f) ;
  messdump ("Partition %s: size %6d, now  %8d blocks,  %s\n",
	    timeShow(timeParse("now"), timeBuf, 25), p->size, p->max, pName(p)) ;
  messfree (emptybp) ;
  p->w = p->r = -1 ;
} /* aDiskFile */

/************************************************************/

static myoff_t aDiskSeek (KEY key, PARTITION *p, DISK d, BOOL isRead)
{  
  myoff_t pp = 0, pos, delta ;
  int ff = -1 ;
  lastDiskLocate = d ;
  isHeaderRead = FALSE ;

  if (isRead)
    { if (p->r < 0)
      { p->r = open(pName(p), O_RDONLY | O_BINARY,0); p->rcurr = 0 ; }
      ff = p->r ;
    }
  else
    { if (p->w < 0)
      { p->w = open(pName(p), O_WRONLY | O_BINARY,0); p->wcurr = 0 ; p->wlastpos = 0 ; }
      ff = p->w ;
    }
  if (ff < 0)
    diskcrash (p, "aDiskSeek cannot locate block") ;  
 
  pos = (key == SUPERKEY) ? 
    0 : SUPERBLOCKSIZE + (myoff_t)p->position * (myoff_t)p->size ;
  /*
  delta = pos  ; 
  if((pp = lseek(ff, delta, SEEK_SET)) != pos)
      diskcrash(p, "aDiskSeek error") ;
  */
  delta = pos  - (isRead ? p->rcurr : p->wcurr) ;
  if (delta)
    if((pp = lseek(ff, delta, SEEK_CUR)) != pos)
      messcrash("aDiskSeek error %s %s delta=%d, pos = %d, pp = %d p->cur = %d", pName(p), isRead ? "r": "w", delta, pos, pp, (isRead ? p->rcurr : p->wcurr)) ;
  if (isRead) p->rcurr = pos ; 
  else { p->wcurr = pos ; p->wlastpos = p->position ; }
  return pp ;
} /* aDiskSeek */

/************************************************************/

static void aDiskDoRead (PARTITION *p, char *vp, int size)
{ 
  int nb = -1 ;			/* number of bytes actually read */
  int nretry = 0 ; 

#if defined(MACINTOSH)
  if ((nb = FMPBread (p->r, 0, fsFromMark, vp, size)) != size)
    diskcrash (p, "Read failure") ;
#else
  p->rcurr += size ;
  while((nb = read (p->r,vp,size)) != size)
    { if (nretry++ == 5)
	diskcrash (p, "Read failure after 5 retries") ;
      if (nb > 0)
	{ vp += nb ;
	  size -= nb ;
	}
      if (size < 0) 
	diskcrash (p, "Read failure, negative size") ;
    }
#endif
} /* aDiskDoRead */

/************************************************************/

void aDiskDoWrite (PARTITION *p, char *vp, int size)
{ 
  int nb = -1 ;			/* number of bytes actually written */
  int nretry = 0 ; 

#if defined(MACINTOSH)
  if ((nb = FMPBwrite (p->w, 0, fsFromMark, vp, size)) != size)
    diskcrash ("Write failure") ;
#else
  seteuid (euid) ;
  p->wcurr += size ;
  while ((nb = write(p->w, vp, size)) != size)
    { if (nretry++ == 5)
	diskcrash (p, "Write failure after 5 retries") ;
      if (nb > 0)
	{ vp += nb ;
	  size -= nb ;
	} 
      if (size < 0) 
	diskcrash (p, "Write failure, negative size") ;
    } 
  seteuid (ruid) ;
#endif
}

/**************** low level read/write ends *****************/
/************************************************************/
/************************************************************/
/************************************************************/

/**************** public routines *********************/

DISK aDiskAssign (KEY key, 
		int n1, int s1,
		int n2, int s2)
{ 
  DISK d = 0 ;

  if (isWriteAccess())
    aDiskDoAssign (&d, key, n1*s1, n2*s2) ;
  return d ;
} /* aDiskAssign */

/**************************************************************/

DISK aDiskWrite (KEY key, KEY parent, 
		char* p1, int n1, int s1,
		char* p2, int n2, int s2)
{ 
  DISK d ;
  PARTITION *p ;
  DISK_HEADER h ;
  int n ; 

  if (!isWriteAccess())
    return 0 ;

  nWrite++ ;
  currentKey = key ;		/* for error messages */

  h.key = key ; h.parent = parent ;
  h.n1 = n1 ; h.s1 = s1 ;
  h.session = currentSession ;
  if (n2) 
    { h.n2 = n2 ; h.s2 = s2 ; }
  else
    { h.n2 = 0 ; h.s2 = 0 ; }

  if (swapData)
    aSwapHeader (&h) ;

  p = aDiskDoAssign (&d, key, n1*s1, n2*s2) ;
  /* myoff_t pp1 = */ aDiskSeek (key, p, d, FALSE) ;
  h.disk = d ;
  if (key == SUPERKEY) h.parent = 1 ;  /* unswap ! */
  aDiskDoWrite (p, (char*)(&h), sizeof(DISK_HEADER)) ;
  if (n1 && s1)
    aDiskDoWrite (p, p1, n1*s1) ;
  if (n2 && s2)
    aDiskDoWrite (p, p2, n2*s2) ;

  n = sizeof(DISK_HEADER) + h.n1 * h.s1 + h.n2 * h.s2 ;
  nBytesWritten += n ; 

  /* test */
  if (key == __Global_Bat) 
    {
    int u1, u2, u3, u4 ; KEY k2 ;

    if (debug) freeOutf("// aDiskWrite:_Global_Bata at %d:%d ", p->number, p->position) ;
    aDiskReadHeader (d,key,&k2,&u1,&u2,&u3,&u4) ;
    if (u1 != h.n1 || u2 != h.s1 || u3 != h.n2 || u4 != h.s2 || k2 != h.parent)
      messcrash("cannot read header back %s:%s", className(key), name(key)) ;

  }
  return d ;
} /* aDiskWrite */

/**************************************************************/
/****** read in two routines, 
        so calling function can allocate correct memory size
*******/

void aDiskReadHeader (DISK d, KEY key, 
		     KEY *parent, 
		     int *n1, int *s1, 
		     int *n2, int *s2)
{ 
  DISK_HEADER h ;
  PARTITION *p = 0 ;

  nRead++ ;
  currentKey = key ;		/* for error messages */

  lastP = p = aDiskLocate (d) ;
  if (0 && debug) 
    printf("aDiskReadHeader at %d:%d ", p->number, p->position) ;
  aDiskSeek (key, p, d, TRUE) ;
  aDiskDoRead (p, (char*)(&h), sizeof(DISK_HEADER)) ;

  if (key == SUPERKEY)
    { 
      KEY sex = h.parent ;
      swapData = FALSE ;
      if (sex != 1)
	{ swapData = TRUE ; sex = swapKEY(sex) ; }
      if (sex != 1)
	messcrash ("this disk has no sex appeal") ;
    }
  if (swapData)
    aSwapHeader (&h) ;

  if (h.disk != d)
    messerror ("aDiskReadHeade: Disk address mismatch on key: %s:%s", className(key), name(key)) ;
  if (h.key != key)
    messcrash ("aDiskReadHeade: key mismatch: searching %s:%s,  got %s:%s",  
	       className(key), name(key), className(h.key),name(h.key)) ;

  *parent = h.parent ;
  *n1 = h.n1 ; *s1 = h.s1 ;
  *n2 = h.n2 ; *s2 = h.s2 ;

  if (debug)  printf(" %d bytes, d=%d, p->r=%d, %s\n", 
		     d, (*n1) * (*s1) + (*n2) * (*s2), p->r, name(key)) ;
  isHeaderRead = TRUE ;
  nBytesRead += sizeof(DISK_HEADER) + h.n1 * h.s1 + h.n2 * h.s2 ;
} /* aDiskReadHeader */

/***************/
 
void aDiskReadData (DISK d,
		   char *p1, int n1, int s1,
		   char *p2, int n2, int s2)
{ 
  if (lastDiskLocate != d || !isHeaderRead)
    messcrash ("diskReadData() not preceeded by diskReadHeader()") ;

  if (n1 && s1)
    aDiskDoRead (lastP, p1, n1*s1) ;
  if (n2 && s2)
    aDiskDoRead (lastP, p2, n2*s2) ;
  isHeaderRead = FALSE ;
} /* aDiskReadData */

/************************************************************/
/************************************************************/
/************************************************************/
/************************************************************/
/*
The following code exists in disknew and should probably 
be duplicated in some form or another
*/
/************************************************************/
/************************************************************/
/************************************************************/
/**************** public functions **************************/

BOOL saveAll (void)
{
  int i = 0, j = 3 ;
  BOOL modif = FALSE ;

  /* ici, in regarde si une session est ouverte */
  if(!isWriteAccess())
    return FALSE ;

  while(j--)
    { int oldBlocks = nWrite ;
      if(batIsModified)
	{ j+=3 ;  /* Three additional loops for security ! */
	  modif = TRUE ;
	  aDiskBATwrite() ;
	}
      cacheSaveAll();  /* Secondary Cache */
      aCacheSaveAll();  /* Secondary Cache */
      
      bIndexSave () ;
      lexSave () ;                   /* Lexiques go to cache */
      if (modif)
	lexmark(1) ;
      if(modif)
	sessionRegister() ;
      if (ps_debug)
	printf("saveAll: pass %d: wrote %d blocks\n", i++, 
	       nWrite - oldBlocks);
    }
  return modif ;
} /* saveAll */

/************************************************************/
/**************************************************************/
/* Read the global BAT, after lexi(0)  */
void diskPrepare(void)
{ 
  BitSet bat = 0 ;
  Array caption = 0 ;
  DISK d = lexDisk(__Global_Bat) ;
  int position, number ;
  char *cp = 0 ;
  PARTITION *p ;

  if (!d) return ;  /* case of session 1 */
  /* since this partition does not yet exists at this stage, we must make it up */

  number = d >> 20 ; position = d & 0xfffff ;    
  p = arrayp (partitions, number, PARTITION) ;
  p->number = number ;
  p->max = position + 100 ; /* wild guess, temporary */
  p->r = p->w = -1 ;
  /* hard bit is to get p->size 
  size = MIN_PARTITION_SIZE ;   
  while (size < (1<<30))
    {
      if ((cp = sessionFilName
	  (messprintf("database/a_%d_%d",size,number + 1),"ace5","r")))
	break ;
      size <<= 2 ;
    }
*/
  if (debug1) messdump ("\nReading BAT, address %d, %s\n",  lexDisk(__Global_Bat), pName(p)) ;
  p->size =  globalBatFileSize ;
  cp = sessionFilName (pName(p), 0, 0) ;
  if (!cp)
    messcrash("Cannot a  the bat file of the form $ACEDB/database/a_*_%d\n",
	      number) ;
  aDiskBigArrayGet(__Global_Bat, 0, &caption, &bat) ;
  if (debug1) messdump ("aDiskArrayGet returned a bat of size %ld\n", bigArrayMax(bat)) ;
  if(!bat || !bigArrayMax(bat))
    {
      if(!isWriteAccess())
	messcrash("Cannot read the BAT");
    }
  if (partitions) arrayDestroy (partitions) ;
  partitions = aDiskBatDecompress (caption, bat) ;
  arrayDestroy (caption) ;  bitSetDestroy (bat) ;
} /* diskPrepare */

void diskShutDown (void)
{


} /* diskShutDown */

/**************************************************************/
    /* query extends the size of the main database file */

/**************************************************************/
/*************************************************************/

int aDiskGetGlobalAdress (void) 
{ 
  return globalAddress ;
} /* aDiskGetGlobalAdress */

/*************************************************************/

static void  aDiskReadSuperBlock(BP bp)
{
  int i ;
  Array a = 0 ;
  char *cp, *cq ;

  if (!aDiskArrayGet (SUPERKEY, 0, &a, 0))
    messcrash ("Cannot find the superblock") ;

  globalAddress = bp->gAddress  = array(a,0,KEY) ;
  bp->mainRelease = array(a,1,KEY) ;
  bp->subDataRelease = array(a,2,KEY) ;
  bp->subCodeRelease = array(a,3,KEY) ;
  bp->subSubCodeRelease = array(a,4,KEY) ;
  bp->session = currentSession =  array(a,5,KEY) ;
  globalFileNumber = array (a,6,KEY) ;
  globalFileSize = array (a,7,KEY) ;
  globalFileMax = array (a,8,KEY) ;
  globalBatFileSize = array (a,9,KEY) ;
  cp = (char*)(arrp(a,10,KEY)) ;
  cq = bp->dbName ;
  i = 32 ; while (i--) *cq++ = *cp++ ;
} /* aDiskReadSuperBlock */

/*************************************************************/
/* Called from session WriteSuperBlock */  
void diskWriteSuperBlock(BP bp)
{
  Array a = arrayCreate(10,KEY) ;
  PARTITION *p = aDiskLocate (bp->gAddress) ;
  char *cp, *cq ;
  int i ;

  array(a,0,KEY) = (KEY) bp->gAddress ;
  array(a,1,KEY) = (KEY) bp->mainRelease ;
  array(a,2,KEY) = (KEY) bp->subDataRelease ;
  array(a,3,KEY) = (KEY) bp->subCodeRelease ;
  array(a,4,KEY) = (KEY) bp->subSubCodeRelease ;
  array(a,5,KEY) = (KEY) bp->session ;
  array(a,6,KEY) = globalFileNumber =  p->number ;
  array(a,7,KEY) = globalFileSize = p->size ;
  array(a,8,KEY) = globalFileMax = (KEY) p->max ;

  p = aDiskLocate (lexDisk(__Global_Bat)) ;
  array(a,9,KEY) = globalBatFileSize = p->size ;

  array(a,10,KEY) = 0 ;
  array(a,11,KEY) = 0 ; 

  cp = (char*)(arrp(a,10,KEY)) ;
  cq = bp->dbName ;
  i = 32 ; while (i--) *cp++ = *cq++ ;
  
  if (swapData) 
    {
      i = 10 ;
      while (i--) arr(a,i,KEY) = swapKEY(arr(a,i,KEY)) ;
    }
  
  aDiskClose() ;
  aDiskArrayStore (SUPERKEY, 0, a, 0) ;
  aDiskClose() ;
  arrayDestroy (a) ;

  /*  bat plus/minus reinitialise */
  i = arrayMax(partitions) ;
  while (i--)
    { 
      unsigned long int max = p->max ;

      p = arrp(partitions, i, PARTITION) ;  
      p->batPlus = bitSetReCreate (p->batPlus, max) ;
      p->batMinus = bitSetReCreate (p->batMinus, max) ;
      bitUnSet(p->batPlus, max - 1) ; /* make room */
      bitUnSet(p->batMinus, max - 1) ; /* make room */ 
    }
} /* diskWriteSuperBlock */

/**************************************************************/
/*              lecture de la table des blocs                 */
/**************************************************************/

static Array aDiskBatDecompress (Array caption, BitSet bat)
{
  long int max, i1, ii, jj = 0, i, n1 = 0 ;
  PARTITION *p ;
  Array pp ;
  unsigned int *u1p, *u2p ;

  if (debug2) printf ("aDiskBatDecompress\n") ;
  u2p = bigArrp(bat, 0, unsigned int) ;
  ii = array(caption, jj++, int) ;
  pp = arrayCreate(ii, PARTITION) ;
  for (i1 = 0 ; i1 < ii ; i1++)
    {
      p = arrayp(pp, i1, PARTITION) ;
      p->number = i1 ;
      p->size = array(caption, jj++, int) ;
      p->max = max = array(caption, jj++, int) ;
      if (debug2) printf ("\t%lu", max) ;
      p->r = p->w = -1 ;  
      p->full = FALSE ;
      p->bat = bitSetCreate (max, 0) ;   
      p->batPlus = bitSetCreate (max, 0) ;   
      p->batMinus = bitSetCreate (max, 0) ;   
      bitUnSet(p->bat, max - 1) ; /* make room */
      bitUnSet(p->batPlus, max - 1) ; /* make room */
      bitUnSet(p->batMinus, max - 1) ; /* make room */

      u1p = bigArrp(p->bat, 0, unsigned int) ;
      i = bigArrayMax(p->bat) ; n1 += i ;
      /* i = (max + 31) >> 5; n1 += i ; */
      if (n1 > arrayMax(bat))
	messcrash ("inconsistency in the stored bat") ;
      while (i--) *u1p++ = *u2p++ ;
    }
  if (debug2) printf ("\nn1=%lu bigArrayMax(bat)=%ld\n", n1, bigArrayMax(bat)) ;
  if (n1 != bigArrayMax(bat)) /* should it be but this crash at init */
    messcrash ("inconsistency in the stored bat") ;
  if (debug1) messdump("DiskBatRead found %d partitions\n", ii) ;
  return pp ;
} /* aDiskBatDecompress  */

/********************************************/

static void aDiskBatCompress (Array pp, Array caption, BitSet bat, int state)
{
  PARTITION *p ;
  int ii, jj ;
  long int i, n1, n2 ;
  unsigned int *u1p, *u2p ;
  BitSet bb ;

  /* concatenate all the bats */  
  if (debug2) printf ("aDiskBatCompress\n") ;
  p = arrp(pp, 0, PARTITION)  - 1 ;
  ii = arrayMax(pp) ;
  jj = 0 ; array(caption, jj++, int) = ii ;
  n1 = n2 = 0 ;
  while (p++, ii--)
    { 
      array(caption, jj++, int) = p->size ;
      array(caption, jj++, int) = p->max ; 
      if (debug2) printf ("\t%d", p->max) ;
      n2 += bigArrayMax(p->bat) ; 
      if (debug2) printf ("aDiskBatCompress  p->size=%d p->max=%d n2=%lu\n", p->size,  p->max,  n2) ;
      bigArray(bat,n2 - 1,unsigned int) = 0 ; /* make room */
      u1p = bigArrp(bat, n1, unsigned int) ; /* may have been reallocated */
      bb = 0 ;
      switch(state)
	{
	  case 1: bb = p->bat ; break ;
	  case 2: bb = p->batPlus ; break ;
	  case 3: bb = p->batMinus ; break ;
	}
      if (bb)
	{
	  u2p = bigArrp(bb, 0, unsigned int) ;
	  i = bigArrayMax(bb) ; n1 += i ; 
	  if (debug2) printf (" bigArrayMax(bb) = %lu", i) ;
	  while (i--) *u1p++ = *u2p++ ;
	}
    }
  if (debug2) printf ("\nn1=%lu bigArrayMax(bat)=%ld\n", n1, bigArrayMax(bat)) ;
} /* aDiskBatCompress */

/********************************************/

static void aDiskBATwrite(void)
{  
  KEY kPlus, kMinus ;
  BitSet bat = bitSetCreate(16000*8, 0) ;
  Array caption = arrayCreate (12, int) ;

  if (!batIsModified) 
    return ;
  batIsModified = FALSE ; /* must come before actual write */ 
  if (!partitions) 
    return ;

  lexaddkey(messprintf("p-%d",thisSession.session),&kPlus,_VBat) ;
  lexaddkey(messprintf("m-%d",thisSession.session),&kMinus,_VBat) ;

  caption = arrayReCreate (caption, 48, int) ;
  bat = bitSetReCreate (bat, 16000*8) ;
  aDiskBatCompress (partitions, caption, bat, 1) ;
  if (debug1) messdump ("\naDiskBatCompress prepared a bat of size %ld\n", bigArrayMax(bat)) ;
  aDiskBigArrayStore(__Global_Bat, 0, caption, bat) ;
  if (debug1) messdump("Writing BAT, address %d, size %ld\n", lexDisk(__Global_Bat), bigArrayMax(bat)) ;

  caption = arrayReCreate (caption, 48, int) ;
  bat = bitSetReCreate (bat, 16000*8) ;
  aDiskBatCompress (partitions, caption, bat, 2) ;
  aDiskBigArrayStore(kPlus, 0, caption, bat) ;

  caption = arrayReCreate (caption, 48, int) ;
  bat = bitSetReCreate (bat, 16000*8) ;
  aDiskBatCompress (partitions, caption, bat, 3) ;
  aDiskBigArrayStore(kMinus, 0, caption, bat) ;
  arrayDestroy (caption) ;
  bitSetDestroy (bat) ;
} /* aDiskBATwrite */

/********************************************/

static BOOL aDiskFuseBats(BitSet old, BitSet new)
{ 
  register long int i;
  register unsigned int *o , *n ;
  BOOL result = FALSE ;

  if( !bigArrayExists(old) || !bigArrayExists(new) ) 
    return FALSE ;
  o = bigArrp(old,0,unsigned int) - 1 ;
  n = bigArrp(new,0,unsigned int) - 1 ;

  if(bigArrayMax(old) > bigArrayMax(new))
    messcrash("Inconsistency in fuseBats") ;
    
  i = bigArrayMax(old) ; 
  while(o++, n++, i--)
    if (*o)
      { if(*o & *n)
	  messcrash("Duplication in fuseBats") ;
	result = TRUE ;
	*n |= *o ;
      }
    
  return result ;
} /* aDiskFuseBats */

/********************************************/

static BOOL aDiskFusePartitions (Array ppf, Array pps)
{ 
  int i ;
  BOOL result = FALSE ;
  PARTITION *pp, *ss ;

  i = arrayMax (ppf) ;
  if (i > arrayMax(ppf))
    messcrash ("aDiskFusePartitions, father's larger than son's") ;
  while (i--)
    {
      pp = arrp (ppf, i, PARTITION) ;
      ss = arrp (pps, i, PARTITION) ;

      if (pp->number != ss->number || 
	  ss->size != pp->size)
	messcrash ("Size inconsistency in aDiskFusePartitions") ;
      result |= aDiskFuseBats (pp->bat, ss->bat) ;
    }
  return result ;
} /* aDiskFusePartitions */

/********************************************/
   /* compress block which have been allocated and freed
      during the same session, i.e. necessary after a batFusion
      */
static void  aDiskCoalesceBats (BigArray batArray, BigArray plusA, BigArray minusA)
{ 
  register long int i;
  register unsigned int *plus, *minus, *bat , c ;

  bat = bigArrp(batArray,0,unsigned int) ;
  plus = bigArrp(plusA,0,unsigned int) ;
  minus = bigArrp(minusA,0,unsigned int) ;

  i = bigArrayMax(plusA) ;
  while(i--)
    { c = *minus & *plus ;  /* reset c bits to zero */
      if(c)
	{ 
	  batIsModified = TRUE ;
	  if ( (*bat & c) != c )
	    messcrash("diskBatCompress releasing a non set block %c %c",
		      *bat , c ) ;
	  *bat &= ~c ;  /* tilda NOT c) */
	  *plus &= ~c ;   /* tilda NOT c) */
	  *minus &= ~c ;  /* tilda NOT c) */
	}
      bat++ ; plus++ ; minus++ ;
    }
} /* aDiskCoalesceBats */

/********************************************/

static void aDiskCoalescePartitions (Array ppp, Array ppm)
{ 
  int i ;
  PARTITION *rr, *pp, *mm ;

  i = arrayMax (ppp) ;
  if (i != arrayMax(ppm))
    messcrash ("aDiskCoalescePartitions inconsistency") ;
  if (i > arrayMax(partitions))
    messcrash ("aDiskCoalescePartitions imain nconsistency ") ;
  while (i--)
    {
      pp = arrp (ppp, i, PARTITION) ;
      mm = arrp (ppm, i, PARTITION) ;
      rr = arrp (partitions, i, PARTITION) ;

      if (pp->number != mm->number || 
	  pp->size != mm->size ||
	  pp->number != rr->number ||
	  pp->size != rr->size)
	messcrash ("Size inconsistency in aDiskCoalescePartitions") ;
       aDiskCoalesceBats (rr->bat, pp->bat, mm->bat) ;
    }
} /* aDiskCoalescePartitions */

/********************************************/

static void aDiskDestroyPartition (Array ppf)
{ 
  int i ;
  PARTITION *pp ;

  if (!ppf) return ;
  i = arrayMax (ppf) ;
  if (i > arrayMax(ppf))
    messcrash ("aDiskFusePartitions, father's larger than son's") ;
  while (i--)
    {
      pp = arrp (ppf, i, PARTITION) ;
      bitSetDestroy (pp->bat) ;
      bitSetDestroy (pp->batPlus) ;
      bitSetDestroy (pp->batMinus) ;
    }
  arrayDestroy (ppf) ;
} /* aDiskDestroyPartition */

/**************************************************************/

BOOL diskFuseBats(KEY fPlus, KEY fMinus, KEY sPlus, KEY sMinus,
		  int *nPlusp, int *nMinusp)
{
  BitSet 
    batfp = 0, batsp = 0,
    batfm = 0, batsm = 0;

  Array
    ppfp = 0, ppsp = 0, ppfm = 0, ppsm = 0,
    captfp = 0, captsp = 0,
    captfm = 0, captsm = 0 ;

  BOOL done = FALSE, modif = FALSE ; 
 
  *nPlusp = *nMinusp = 0 ; 
  aDiskBigArrayGet (fPlus, 0, &captfp, &batfp) ;
  aDiskBigArrayGet (sPlus, 0, &captsp, &batsp) ;
  aDiskBigArrayGet (fMinus, 0, &captfm, &batfm) ;
  aDiskBigArrayGet (sMinus, 0, &captsm, &batsm) ;
  
  if (!batfp || !batfm || !batsp || !batsm)
    goto abort ;

  ppfp = aDiskBatDecompress (captfp, batfp) ;
  ppfm = aDiskBatDecompress (captfm, batfm) ;
  ppsp = aDiskBatDecompress (captsp, batsp) ;
  ppsm = aDiskBatDecompress (captsm, batsm) ;

  modif |= aDiskFusePartitions (ppfp, ppsp) ;
  modif |= aDiskFusePartitions (ppfm, ppsm) ;
  if (modif)
    aDiskCoalescePartitions (ppsp, ppsm) ;
 
  *nPlusp = aDiskBatSet(ppsp, 1) ;
  *nMinusp = aDiskBatSet(ppsm, 1) ;

  captsp = arrayReCreate (captsp, 48, int) ;
  batsp = bitSetReCreate (batsp, 16000 * 8) ;
  aDiskBatCompress (ppsp, captsp, batsp, 1) ;
  aDiskBigArrayStore (sPlus, 0, captsp, batsp) ;

  captsm = arrayReCreate (captsm, 48, int) ;
  batsm = bitSetReCreate (batsm, 16000 * 8) ;
  aDiskBatCompress (ppsm, captsm, batsm, 1) ;
  aDiskBigArrayStore(sMinus, 0, captsm, batsm) ;

  arrayKill (fPlus) ;  /* not valid if there were session branches */
  arrayKill (fMinus) ; 


 done = TRUE ;

abort:
  aDiskDestroyPartition (ppfp) ;
  aDiskDestroyPartition (ppfm) ;
  aDiskDestroyPartition (ppsp) ;
  aDiskDestroyPartition (ppsm) ;

  arrayDestroy (captfp) ;  bitSetDestroy (batfp) ;
  arrayDestroy (captfm) ;  bitSetDestroy (batfm) ;
  arrayDestroy (captsp) ;  bitSetDestroy (batsp) ;
  arrayDestroy (captsm) ;  bitSetDestroy (batsm) ;

  return done ;
} /* diskFuseBats */

/**************************************************************/

void diskFuseOedipeBats(Array fPlusArray, Array fMinusArray) 
{ 
  /* BOOL modif = FALSE ;
    
    modif |= aFuseBats(fPlusArray, plusArray) ;
  modif |= aFuseBats(fMinusArray, minusArray) ;
  
  if(modif)
     diskBATcompress(plusArray, minusArray) ;
     */
} /* diskFuseOedipeBats */

/**************************************************************/
/* call CloseDataBase to flush any client side cached copy of the blocks,
   needed to close a concurrency hole in session.c */
void diskFlushLocalCache(void)
{
  aDiskClose () ;
} /* diskFlushLocalCache */

/************************************************************/
/************************************************************/
static void aDiskBigArrayStore (KEY key, KEY parent, Array a, BigArray b)
{
  messcrash ("aDiskBigArrayStore is not written") ;
}

void aDiskArrayStore (KEY key, KEY parent, Array a, Array b)
{
  if (a && b)
    aDiskWrite (key, 0, a->base, a->max, a->size, b->base, b->max, b->size) ;
  else if(a)
    aDiskWrite (key, 0, a->base, a->max, a->size, 0, 0, 0) ;
  else if (b)
    aDiskWrite (key, 0, 0, 0, 0, b->base, b->max, b->size) ;
  else
    messcrash("aDiskArrayStore called with a=b=0") ;
} /* aDiskArrayStore */

/************************************************************/

static BOOL aDiskBigArrayGet (KEY key, KEY *pp, Array *ap, BigArray *bp) 
{
  messcrash ("aDiskBigArrayGet  no idea why we have 2 arrays in this call, 2015_09_12") ;
  return FALSE ;
}


BOOL aDiskArrayGet (KEY key, KEY *pp, Array *ap, Array *bp) 
{
  int n1, n2, s1, s2 ;
  Array a = 0, b = 0 ;
  char *cp = 0, *cq = 0 ;
  DISK d = 0 ;
  KEY parent = 0 ;

  d = (key == SUPERKEY) ? 1 : lexDisk(key) ;
  if (!d) return FALSE ;
  aDiskReadHeader (d, key, &parent, &n1, &s1, &n2, &s2) ;

  if (n1 && s1)
    { 
      a = uArrayCreate (n1, s1, 0) ;
      arrayForceFeed (a, n1) ;
      cp = a->base ;
    }
  if (n2 && s2)
    { 
      b = uArrayCreate (n2, s2, 0) ;
      arrayForceFeed (b, n2) ;
      cq = b->base ;
    }
  aDiskReadData (d, cp, n1, s1, cq, n2, s2) ;

  if (swapData && s1==4) 
    { 
      int i = n1 ;
      while (i--) arr(a,i,KEY) = swapKEY(arr(a,i,KEY)) ;
    }

  if (ap) *ap = a ; else arrayDestroy (a) ;
  if (bp) *bp = b ; else arrayDestroy (b) ;
  if (pp) *pp = parent ;
  return TRUE ;
} /* aDiskArrayGet */

/************************************************************/

void aDiskKill (KEY key) 
{
  DISK d = 0 ;

  d = lexDisk(key) ;
  if (!d) return ;

  aDiskFree(aDiskLocate(d)) ;
  lexSetDisk(key, 0) ;
} /* aDiskKill */

/************************************************************/
/* to thest aDisk,c */
void testAdiskArrayStore (KEY key, Array a)
{
  Array b ;
  int n1, n2, s1, s2 ;
  char *cp, *cq ;
  DISK d = 0 ;
  KEY parent = 0 ;

  d = aDiskWrite (key, 0, a->base, a->max, a->size, 0, 0, 0) ;

  aDiskReadHeader (d, key, &parent, &n1, &s1, &n2, &s2) ;
  if (n1 != a->max || s1 != a->size || n2 || s2)
    messcrash ("error in testArrayStore") ;
  
  b = uArrayCreate (n1, s1, 0) ;
  arrayForceFeed (b, n1) ;
  cp = b->base ;
  aDiskReadData (d, cp, n1, s1, 0, 0, 0) ;
  cp = b->base ;
  cq = a->base ;
  n2 = n1*s1 ;
  while (n2--) if (*cp++ != *cq++) messcrash ("test error") ;
  arrayDestroy (b) ;
} /* testAdiskArrayStore */

/************************************************************/
/************************************************************/
/************************************************************/
