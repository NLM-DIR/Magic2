/*  File: disknew.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
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
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Apr  9 14:38 2003 (edgrif)
 * * Sep 21 12:04 1999 (edgrif): Add sync write mode code for disk writes.
 * * May 10 10:22 1999 (edgrif): Altered dataBaseReadDefinition for Jean.
 * * Apr 16 16:08 1999 (fw): fixed very subtle bug where totalNbpartions
 *              was wrong, because the global variable wasn't reset
 *              probably cause for SANgc02258 and other diskextend bugs
 * * Oct 15 11:33 1998 (fw): included the output from
 *              messSysErrorText in messcrash texts
 *               resulting from disk operations or other OS calls
 * *       *       *        *        *        *
 *	-	extern uid_t euid, ruid moved to session.h
 * * Jun  3 12:56 1996 (rd)
 * * Apr 11 - 94 : multi-partitions version. (pierre Simondon)
 * * Jan 22 10:33 1992 (mieg): Modify bat and session obj en place.
 * Created: Wed Jan 22 10:33:50 1992 (mieg)
 *-------------------------------------------------------------------
 */

#ifdef ACEDB5
void compilersDoNotLike (void *emptyFiles) {return ; }

#else			/* ACEDB 4 only */

/***************************************************************/
/***************************************************************/
/**  File disksubs.c :                                        **/
/**  Handles the disk of the ACeDB program.                   **/
/***************************************************************/
/***************************************************************/
/*                                                             */
/*  Six routines are public :                                  */
/*     avail,prepare, alloc/free, blockread/blockwrite.        */
/*                                                             */
/*  A single long sequential file large enough to hold         */
/*  DISKSIZE blocks is "prepared" once and for all.            */
/*                                                             */
/*  "diskalloc", using the static BAT table, returns a DISK    */
/* value which is the number of the first free block following */
/*   This number is then stored in the lexique as lxi->diskaddr*/
/* and used, after multiplication by BLOC_SIZE as an offset    */
/* in "diskblockread/write", which use fseek.                  */
/*  "diskfree" returns a bloc to the BAT                       */
/*   diskavail gives disk statistics.                          */
/*                                                             */
/*   The BAT table handling is private.                        */
/*                                                             */
/***************************************************************/


/*
Acedb databases are made up of a set of partition files.  A partition
file contains some number of blocks of data.  To read/write a block,
you open the appropriate partition file and seek to the right offset.

wspec/database.wrm contains a definition of the list of partition files
that a database may use.  database/database.map contains a list of
partition files that are actually being used, along with a count of
blocks stored in each one.  You can add more partition files to
database.wrm; the new files will be available next time you start the
database.

This file contains functions to:

- initialize the block input/output subsystem.  
	dataBaseAccessInit () for an existing database
	dataBaseCreate () for a new database

- disk block input/output
	diskblockread ()
	diskblockwrite ()
		The parameter is a BP, which is typedef struct block.
		 (there are several struct block definitions; I am not
		sure which it is.)  You can only read/write the specific
		database structure, not arbitrary data.

		These functions fix the data if the database was written
		on a machine with a different byte order.  They know
		quite a bit about the database.

- block allocation tables
	There are several functions that handle block allocation.

- closing/flushing files


This subsystem gives the appearance of a single very large data store
with data addressed by block number.  The blocks are handled by a "BAT"
 (block allocation table).  This is an array of bytes.  Each byte
represents 8 blocks.  The lsb represents the lowest numbered block in
the byte.  1 is allocated, 0 is free.

The database has these things called "sessions".  A session has a regular
BAT, a "plus BAT" which identifies blocks that have been allocated during
this session, and a "minus BAT" which identifies blocks that have been
deallocated during this session.  I don't really understand why free blocks
are tracked per-session instead of one table for the whole database, but
the "plus BAT" and "minus BAT" are somehow involved in tracking what blocks
are available.  See diskalloc () for how they are used.

For each partition file, we know the block number of the first block in
the file.  For a partition file with N blocks in it, the next N bits in
the BAT apply to blocks in that file.

Or, said another way,
	block N is in file X, 
		where	N >= start block of file X 
		and	N < start block of file X + 1
	block N is at file offset Y, where 
		Y = (N - start block of X ) * BLOC_SIZE

When all the blocks are used, the system tries to increase the number of
blocks available by extending the last partition file.  (You cannot
extend any partition file but the last, because of the above relations.)
If the partition file has not reached it's maximum size, it is extended
either to it's maximum size or by a growth increment, whichever is smaller.

When extending the database, the BAT is also made bigger and all the new
blocks are marked free.


There is a special case of copy-on-write databases.  In this case, the
database consists of a set of "new" files that are the partition files
for this database and a set of "old" files that are the partition files
for some older read-only database.  The system attempts to use the old
files whenever possible.  When the database wants to write to a block in
an old file, it copies the old file to make a new one.  From that point
on, only the new file is ever referenced.

 (The block allocation algorithm does not currently know about this
feature; it could try to preferentially allocate blocks from files that
have already been copied.)


Other notes:

seteuid () calls are around write () as well as open () because NFS checks access
on each read/write call, not just at open.  This code never had seteuid () calls
around read (), but presumably that is just because the database files are
typically 644 mode.

We do not use setuid programs here at NCBI anyway, so it does not matter to us.
For efficiency, I only make the seteuid () call if euid != ruid.

*/

#define DEF_BP /* avoid typedefing it as void* in disk.h */
typedef struct block *BLOCKP ;
typedef BLOCKP BP ;
#include "acedb.h"
#include "disk__.h"
#include "disk.h"
#include "block.h"
#include "lex_bl_.h"
#include "byteswap.h"
#include "chrono.h"

extern void swapABlock (BP bp) ;
extern void swapBBlock (BP bp) ;

static void diskBATread (void);
static void diskBATwrite (void);
static void diskBATcompress (Array plusA, Array minusA) ;
static void diskExtend (BOOL zeroing) ;
static BOOL isBat (DISK d) ;

/* relative address of a block in a partition (in bytes) */
static long partitionAddress (DISK d, int *partition);

/* full pathname of a partition : */
static char *partitionName (int numPartition);
static char *oldPartitionName (int numPartition);

/* close all open partitions of the base : */
static void closeDataBase ();

/* creation and update database map : */
static BOOL dataBaseReadDefinition (BOOL new, BOOL skip);
static BOOL dataBaseUpdateMap (void);
static void readPartitionMap (void);

int INITIALDISKSIZE = 64 ;
/*
 * INITIALDISKSIZE is the number of blocks that we try to add
 * when increasing the database size.  It must be a multiple
 * of 8 because a single block allocation table (batArray)
 * entry applies to 8 blocks.
 */

static Array batArray  = 0 ;      /* block allocation table */
static Array plusArray  = 0 ;     /* new blocks in this session */
static Array minusArray  = 0 ;    /* blocked freed by this sesion */

static BOOL batIsModified = FALSE ;
static int blocksRead = 0, blocksWritten = 0;

/* syncWrite - indicates whether acedb should open database files for        */
/*             writing using the O_SYNC flag. This should ensure that acedb  */
/*             will be able to detect timeouts for writes (applies to NFS    */
/*             mounted databases especially).                                */
/*  it slows the code ENORMOUSLY */

static int syncWrite = 0 ;

/* This could be an option that was turned on/off interactively, but since   */
/* sync mode writing is at least 5x slower, the option is currently set from */
/* session.c via an entry in passwd.wrm. It's the sort of thing a database   */
/* administrator should set when there are NFS problems.                     */
void diskSetSyncWrite(void)
{
#if !defined (DARWIN) && !defined (MAC_X)
	    syncWrite |= O_SYNC ;
#endif
} 
/* for debug purpose only : */
static int ps_debug  = 0 ;
static int ps_debug1 = 0 ; 

#define DEFAULT_DATABASE_SIZE (1024 * 1024 )
/*
 * if we do not have a wspec/database.wrm file, the database will be made of
 * one file named "blocks.wrm" containing DEFAULT_DATABASE_SIZE blocks.
 *
 * currently 1 GB (i.e. 1 meg of 1k blocks ) - this cannot be more than 2 GB
 * because file offsets are signed 32 bit quantities.
 */



/**************************************************************/

#define LG_NAME          128

/*
* This is all we need to know about a database partition file.  (actually,
* it is more than we need to know, but for compatibility with old databases,
* I have kept unused fields that are stored in database.map and database.wrm)
*
*/
typedef struct {
  char hostname[LG_NAME];
  /*
   * not used for anything meaningful
   */
  char fileSystem[LG_NAME];
  /*
   * not used for anything meaningful, but a lot of work seems to
   * go into making sure it is "" when in memory and "ACEDB" in
   * database.map file
   */
  char fileName[LG_NAME];
  /*
   * name of the file that this partition is in
   */
  int  maxBlocks;
  /*
   * the most blocks that will be allowed in this file
   */
  int  currentBlocks;
  /*
   * the number of blocks currently in this file
   */
  int fd;
  /*
   * file descriptor of the open file for this partition, or -1
   */
  int fd_writeable;
  /*
   * true if fd is open for writing.
   */
  int status_old;
  /*
   * We are using the old file - i.e. not yet copied for write.
   * used for status display only.
   */

} PARTITION;

/**************************************************************/

static int lastUsedPartition;
/*
 * lastUsedPartition is the number of the highest partition that
 * we have actually used.
 */

static int totalNbPartitions ;
/*
 * totalNbPartitions is the number of partition files we have defined
 */

static Array pTable = 0 ;
/*
 * pTable is an array of PARTITION - it contains a record for each
 * partition file that the database is permitted to use.  These files
 * might not all exist.
 */


static char *COW_partition_dir = NULL;
/*
 * COW_partition_dir is a directory where we can find read-only
* COW_partition_dir of this database.  We use the old database files as
* copy-on-write images.
 */

/**************************************************************/

/*
 * get_fd - get an open a file descriptor for a partition file; returns -1
 * on error so you can issue your own message
 *
 * partition is which partition file you want the file descriptor for
 *
 * mode is one of 
 *	O_WRONLY - I want a writable fd
 *	O_RDONLY - I want a readable fd
 *
 * (These values are just used because they are handy.  O_WRONLY really gets
 * you a read/write file descriptor. )
 *
 * This function remembers the open file descriptors.  It tries to optimize
 * the case where you switch from read-only to read/write access by getting
 * read/write when it can.  It also implements copy-on-write databases.
 *
 */
static int get_fd (int partition, int mode )
{
  PARTITION *pp;

  pp = arrayp (pTable, partition, PARTITION);

  /*
  * if the file is already open, see if we can use the existing fd.
  */
  if (pp->fd >= 0)
    {
      /*
       * If it is open for write, the fd is good for anything we need.
       */
      if (pp->fd_writeable)
	  return pp->fd;
      /*
       * If it is open for read and we don't need write, we can use
       * the existing fd.
       */
      if (mode != O_WRONLY)
	  return pp->fd;
      /*
       * ok, we want write but do not have it.  close it and re-open.
       */
      close (pp->fd);
    }

  /*
   * make uid correct for file access
   */
  if (euid != ruid)
    seteuid (euid);

  /*
   * In the event of a copy-on-write database, the file may be from the
   * old database.  We assert now that this is not the case, but change
   * it later if it is.
   */
  pp->status_old = 0;

  /*
   * If we want write access, open the partion name for read/write.
   */
  if (mode == O_WRONLY)
    {
      if (COW_partition_dir)
	{
	  /*
	   * If we have a copy-on-write database, we first try the new 
	   * partition file.  If it is not there, we copy the old partition 
	   * file, which results in a new partition file that we can use.
	   */ 
	  pp->fd = open (partitionName (partition), syncWrite | O_RDWR | O_BINARY , 0644 );
	  if (pp->fd < 0)
	    {
	      int ofd, nfd, n;
	      char b[65536];
	      nfd = open (partitionName (partition), syncWrite | O_RDWR | O_BINARY | O_CREAT , 0644 );
	      if (nfd < 0)
		messcrash ("cannot copy a copy-on-write block file");
	      ofd = open (oldPartitionName (partition), syncWrite | O_RDONLY | O_BINARY , 0644 );
	      while ((n = read (ofd, b, sizeof (b))) > 0 )
		if (write (nfd, b, n) != n)
		  messcrash ("write error copying copy-on-write block file");
	      close (ofd);
	      /*
 	      * do not close the newly copied file -- it is the one we want!
	      */
	      pp->fd = nfd;
	    }
	}
      else
	{
	  /*
	   * If we do not have a copy-on-write database, we just open the 
	   * partition file.
	   */
	  pp->fd = open (partitionName (partition),syncWrite | O_RDWR | O_BINARY | O_CREAT, 0644 );
	}
      pp->fd_writeable = 1;
    }
  else
    {
      if (COW_partition_dir)
	{
	  /*
	   * if this is a copy-on-write database, we try to:
	   *	if we have write access to the database,
	   *	   open the new partition file for writing, even though we
	   * 	   are only reading at the moment
	   *    if we do not have write access to the database,
	   *	   open the new partition file read-only
	   *    if there is no new partition file,
	   *	   open the old partition file read only
	   *
	   * The idea here is that we want the newer file if there is one, 
	   * and we want to avoid closing the file and re-opening it if we 
	   * later try to upgrade our access from read to read/write.
	   *
	   */
	  if (isWriteAccess ())
	    {
	      pp->fd = open (partitionName (partition), syncWrite | O_RDWR | O_BINARY , 0644 );
	      pp->fd_writeable = 1;
	    }
	  else
	    {
	      pp->fd = open (partitionName (partition), syncWrite | O_RDONLY | O_BINARY , 0644 );
	      pp->fd_writeable = 0;
	    }
	  if (pp->fd < 0)
	    {
	      pp->fd = open (oldPartitionName (partition), syncWrite | O_RDONLY | O_BINARY , 0644 );
	      pp->fd_writeable = 0;
	      pp->status_old=1;
	    }
	}
      else
	{
	  /*
	   * if this is not a copy-on-write databse, we just open the 
	   * partition file.
	   */
	  pp->fd = open (partitionName (partition), syncWrite | O_RDONLY | O_BINARY , 0644 );
	  pp->fd_writeable = 0;
	}
    }

  if (euid != ruid)
    seteuid (ruid);

  return pp->fd;
}


/**************************************************************/

static int diskBatSet (Array a) ;

extern VoidRoutine messcrashroutine ; 
extern void simpleCrashRoutine (void);
#define NOCRASH NOCRASH1 ; NOCRASH2 ;
#define NOCRASH1 VoidRoutine cacheCrash = messcrashroutine
#define NOCRASH2 messcrashroutine = simpleCrashRoutine /* don't allow crash recovery */
#define OKCRASH  messcrashroutine = cacheCrash

/*******************************************************************/
/* Partitions configuration : must be called on first session  :   */
/*******************************************************************/

/*
 * There are two functions here to initialize the information we need
 * to find the partition files in the databse.
 *
 * dataBaseAccessInit ()
 *	This one is called if we are opening an existing database.
 *
 * dataBaseCreate ()
 *	This one is called if the database is to be created
 *
 */

BOOL dataBaseCreate ()
{ 
  NOCRASH ;

  if (ps_debug)
    fprintf (stderr, "dataBaseCreate\n");

  pTable = arrayReCreate (pTable, 20, PARTITION) ;

  /*
   * first, we read the partition list from wspec/database.wrm
   */
  if (!dataBaseReadDefinition (TRUE, TRUE))
    messcrash (" dataBase creation : no description file found !!!\n");

  /*
   * then we write the partion list out to database/database.map
   */
  dataBaseUpdateMap () ;

  /*
   * and now that we have a valid database.map file, we read it to
   * get our list of partitions.
   */
  dataBaseAccessInit ();

  OKCRASH ;
  return TRUE;
}

/*******************************************************************/
/* must be called at session open.                                 */
/*******************************************************************/

void  dataBaseAccessInit ()
{ int part;
 PARTITION *pp ;

 /*
  * re-create the pTable because we may already have one if we are
  * creating the entire database.
  */

 pTable = arrayReCreate (pTable, 20, PARTITION) ;

 /*
  * read database/database.map into pTable
  */

 readPartitionMap ();

 /*
  * print various information if debugging
  */

 if (ps_debug)
   {
     fprintf (stderr, ">>dataBaseAccessInit : %d current partitions\n",
	     lastUsedPartition);
     pp = arrayp (pTable, 0, PARTITION) ;
     for (part = 0; part <= lastUsedPartition; part++, pp++) 
       {
	 if (ps_debug)
	   fprintf (stderr,"part %d f system %s,file %s Max %d, curr %d\n",
		   part, pp->fileSystem, pp->fileName, 
		   pp->maxBlocks, pp->currentBlocks);
       }
   }

 /*
  * note that we did not open any partition files - they will be opened
  * on demand.  
  */
 for (part = 0; part < totalNbPartitions; part++)
   arrayp (pTable, part, PARTITION)->fd = -1;

}

/**************************************************************/
/* On first call during a given session, 
   any modified item
   is moved by blocksave to a new disk position, 
   this alters the lex and BAT but the change is not 
   registered on disk, hence the need for a sequence of
   write */

BOOL saveAll (void)
{ extern void invokeDebugger (void) ;
 int i = 0, j = 3 ;
 BOOL modif = FALSE ;
 extern BOOL blockCheck (void), bIndexSave (void) ;
 NOCRASH1 ;

 /* ici, in regarde si une session est ouverte */
 /* (babelfish says: here, in looks at if a session is open ) */
 if (!isWriteAccess ())
   return FALSE ;
 NOCRASH2 ;

 if (!blockCheck ())
   messcrash ("Sorry - bad cache structure - in saveAll") ;

 while (j--)
   { int oldBlocks = blocksWritten;
   if (batIsModified)
     { j++ ;  /* Three additional loops for security ! */
     modif = TRUE ;
     }
   diskBATwrite () ;
      
   cacheSaveAll ();  /* Secondary Cache */
  
   bIndexSave () ;
   lexSave () ;                   /* Lexiques go to cache */
   blockSave (TRUE) ;                  /*  flush the cache */
   if (modif)
     sessionRegister () ;
   if (ps_debug)
     printf ("saveAll: pass %d: wrote %d blocks\n", i++, 
	    blocksWritten - oldBlocks);
   }
 OKCRASH ;
 return modif ;
}

/**************************************************************/
/* If no explicit call is made,
   when the first cache is full, blockalloc
   flushes a block, calling diskwrite calling
   BATread calling diskalloc with a resulting
   fatal knot in the global used blockend pointer
*/

void diskPrepare (void)
{ if (thisSession.session > 1 )   /* Except during session 0/1 */
  diskBATread () ;               /* la BAT existe deja : */
 else 
   { batArray = arrayCreate (100,unsigned char) ;
   plusArray = arrayCreate (100,unsigned char) ;
   minusArray = arrayCreate (100,unsigned char) ;
   /* select first current partition : */
   lastUsedPartition = 0;
   /* reset first partition file : */
   diskExtend (TRUE) ;
   }
}

void diskShutDown (void)
{
  arrayDestroy (batArray) ;
  arrayDestroy (plusArray) ;
  arrayDestroy (minusArray) ;
  arrayDestroy (pTable) ;
}

/**************************************************************/
/*
 * diskExtend increases the size of the database.  This includes
 * writing blocks of 0's to a database file (to force disk space
 * to be allocated) and increasing the size of the BAT (block
 * allocation table) to match the new size of the database.
 *
 * The new disk space is allocated in some reasonable increment,
 * without concern for how much new space is needed.
 */

static void diskExtend (BOOL zeroing)
{
  int i ;
  myoff_t pos,posr ; /* needed to check lseek success */
  unsigned long dd, newdd, block_increment;
  /*
   * dd is how many blocks in the database
   * block_increment is how many blocks we will add
   * newdd is how many blocks in the database after we extend it
   */
  PARTITION *pp =  arrayp (pTable, lastUsedPartition, PARTITION) ;
  /*
   * points at the partition that is eligible to be extended.
   * We can only extend in the last file currently in use.
   */
  void *zero_buffer;
  int zero_buffer_size;
  /*
   * a buffer of zeroes to write out
   */
  int fd;
  /*
   * file descriptor of the extending partition.
   */

  NOCRASH ;

  /*
   * you can't ever get here if you don't have write access
   * to the database.
   */

  if (!isWriteAccess ())
    messcrash ("Disk Extend called without Write Access") ;

  chrono ("diskExtend") ;

  if (!batArray)
    messcrash ("No batArray in diskextend (FALSE)") ;

  /*
   * dd is the current size of the database, in blocks.  (batArray
   * contains one bit for each block.
   */
  dd = 8 * arrayMax (batArray) ;
  if (ps_debug)
    fprintf (stderr, "diskExtend : current blocks number in bat = %ld\n",dd);

  if (pp->currentBlocks / 8 < pp->maxBlocks / 8 )
    {
      /*
       * If the current file does not have it's max number of blocks yet,
       * we will make it bigger.  
       */
    }
  else
    {
      /*
       * the current file is full, so we must start a new one.
       */
      lastUsedPartition++;
      if (ps_debug)
	fprintf (stderr, "diskExtend : creating block file lastUsedPartition = %d\n", lastUsedPartition) ;
      if (lastUsedPartition >= totalNbPartitions)
	messcrash ("out of partitions - add more in wspec/database.wrm\n");
      pp =  arrayp (pTable, lastUsedPartition, PARTITION) ;
    }
  
  /*
   * In either case, find out how many more blocks can fit
   * in this file, but limit it at some reasonably small size.
   */
  block_increment = (pp->maxBlocks / 8 - pp->currentBlocks / 8 ) * 8;
  if (block_increment > INITIALDISKSIZE)
    block_increment = INITIALDISKSIZE & ~7;
  
  /* 
   * get the file descriptor for the partition
   */
  fd = get_fd (lastUsedPartition, O_WRONLY );
  if (fd < 0)
    messcrash ("diskextend can't open blocks %s (%s)",
	       partitionName (lastUsedPartition), messSysErrorText ()); 
  
  /*
   * seek to the place where we believe the end of the file is.
   */
  pos = pp->currentBlocks * BLOC_SIZE ; 
  if ((posr = lseek (fd, pos, SEEK_SET)) != pos)
    messcrash (
	       "Diskextend failed to lseek the end of the partition to be extended %s, pos=%ld, returned=%ld (%s)\n",
	        partitionName (lastUsedPartition), 
	       (unsigned long) pos, (unsigned long)posr, messSysErrorText ()) ;
  
  /*
   * write a block of zeroes there.  This has some beneficial side-effects:
   * - the new space is filled with 0 even if we are reclaiming part of
   *   a data file that we have forgotten about
   * - we know for certain that the disk space has been allocated
   * - we can hope that we get a large number of contiguous blocks in the
   *   file system
   */
  zero_buffer_size = block_increment * BLOC_SIZE;
  zero_buffer = messalloc (zero_buffer_size );

  if (write (fd, zero_buffer, zero_buffer_size ) != zero_buffer_size)
    messcrash ("unable to extend partition file: %s\n",strerror (errno));

  /*
   * now that we have added blocks, we mark the partition with how many
   * more there are.
   */
  pp->currentBlocks += block_increment;

  /*
   * there used to be an fsync () call here, but I specifically do NOT
   * want to fsync the block of 0's.  We are extending the database, so
   * we know for certain that we will be re-writing at least some of those
   * blocks immediately on return from this function.  It is sufficient
   * that we did the write, which caused them to be allocated to actual
   * physical disk blocks.
   */

  /*
   * write out the list of files and how many blocks are in each
   */
  dataBaseUpdateMap () ;

  /* 
   * Now the database is bigger.  The block allocation tables must
   * be made bigger to match.  All the new blocks are free.
   */
  newdd = dd + block_increment;

  if (newdd & 7)
    messcrash ("in diskExtend (), newdd is not a multiple of 8: %d\n",newdd);

  for (i= dd/8; i<newdd/8 ;i++)
    {
      array (batArray,i,unsigned char) = 0;
      array (plusArray,i,unsigned char) = 0;
      array (minusArray,i,unsigned char) = 0;
    }

  /* first 3 blocks reserved : */
  array (batArray,0,unsigned char) |= 7 ;
  batIsModified = TRUE ;
  
  messfree (zero_buffer);

  /*
   *
   * I have removed the code that does this:
   *
   * We have already written blocks of 0's to the file to extend it, and closed
   * it.  Now we open it again and see if it's size matches what we expect.  If
   * it does not, we wait and try again.  Evidently, this is to work around some 
   * horribly buggy NFS implementation.
   */

  chronoReturn () ;
  OKCRASH ;
}

/**************************************************************/
/*************************************************************/
/* Called from session WriteSuperBlock */  
void diskWriteSuperBlock (BP bp)
{
  register unsigned char * p, *m ;
  if (ps_debug) fprintf (stderr, "diskWriteSuperBlock\n");

  closeDataBase ();

  /*  isWriteAccess () is  implied by the call */
  diskblockwrite (bp) ;     /* reopens file */

  closeDataBase ();

  /*  batReinitialise */
  p = arrp (plusArray,0,unsigned char) ;
  m = arrp (minusArray,0,unsigned char) ;
  memset (p, 0, arrayMax (plusArray));
  memset (m, 0, arrayMax (minusArray));
}

void diskFlushLocalCache (void)
     /* call closeDataBase to flush any client side cached copy of the blocks,
	needed to close a concurrency hole in session.c */
{
  closeDataBase ();
}

/**************************************************************/
/*              lecture de la table des blocs                 */
/* (babelfish says: reading of the table of the blocks )      */
/**************************************************************/

static void diskBATread (void)
{
  int i = INITIALDISKSIZE/8 ;
  int j,k ;
  NOCRASH1 ;

  /* le cas ou on a deja la bat */
  if (batArray) 
    return ;
  NOCRASH2 ;
  /* acces au bloc contenant la bat */
  /* access to the block containing the bat */
  batArray = arrayGet (__Global_Bat,unsigned char,"c") ;
  if (!batArray || !arrayMax (batArray))
    {
      if (!isWriteAccess ())
	messcrash ("Cannot read the BAT");
    }

  i = arrayMax (batArray) ;
  /* on cree les bat plus et moins , de meme taille */
  plusArray = arrayCreate (i,unsigned char) ;
  minusArray = arrayCreate (i,unsigned char) ;
  arrayMax (minusArray) = arrayMax (plusArray) = i;
  
  /* DISK address 0 1 2 are for the system */
  array (batArray,0,unsigned char) |= 7;     
  j = lexDisk (__Global_Bat) ;
  k = 1 << (j%8); j = j/8;
  if (array (batArray,j,unsigned char) & k)
    ; 
  /*    messdump ("\n Correct bat self referenced address %d, session %d",
	j, KEY2LEX (__Global_Bat)->addr ?
	KEY2LEX (__Global_Bat)->addr->p->h.session : 0);
  */
  else 
    messcrash ("no bat self reference for address %d",j);
  OKCRASH ;
}

/********************************************/

static void diskBATwrite (void)
{
  KEY kPlus, kMinus ;
  char buf[1000] ;

  if (!batArray || ! batIsModified) 
    return ;
  batIsModified = FALSE ; /* must come before actual write */ 
  if (ps_debug)
    printf ("\nWriting BAT, address %d", 
	   KEY2LEX (__Global_Bat)->dlx.dk) ;
  sprintf (buf, "p-%d",thisSession.session) ;
  lexaddkey (buf,&kPlus,_VBat) ;
  sprintf (buf, "m-%d",thisSession.session) ;
  lexaddkey (buf,&kMinus,_VBat) ;
  arrayStore (__Global_Bat, batArray,"c") ;
  arrayStore (kPlus, plusArray,"c") ;
  arrayStore (kMinus, minusArray,"c") ;
}

/********************************************/
/* compress block which have been allocated and freed
   during the same session, i.e. necessary after a batFusion
*/
static void diskBATcompress (Array plusA, Array minusA)
{ register int i;
 register unsigned char *plus, *minus, *bat , c ;
 NOCRASH ;

 bat = arrp (batArray,0,unsigned char) ;
 plus = arrp (plusA,0,unsigned char) ;
 minus = arrp (minusA,0,unsigned char) ;

 if (arrayMax (batArray) < arrayMax (plusA))
   messcrash ("bat size inconsistency in diskBATcompress") ;
 i = arrayMax (plusA) ;
 while (i--)
   { c = *minus & *plus ;  /* reset c bits to zero */
   if (c)
     { 
       batIsModified = TRUE ;
       if ((*bat & c) != c )
	 messcrash ("diskBatCompress releasing a non set block %c %c",
		   *bat , c ) ;
       *bat &= ~c ;  /* tilda NOT c) */
       *plus &= ~c ;   /* tilda NOT c) */
       *minus &= ~c ;  /* tilda NOT c) */
     }
   bat++ ; plus++ ; minus++ ;
   }
 OKCRASH ;
}

/********************************************/

static BOOL fuseBats (Array old, Array new)
{ register int i;
 register unsigned char *o , *n ;
 BOOL result = FALSE ;
 NOCRASH1 ;

 if (!arrayExists (old) || !arrayExists (new) ) 
   return FALSE ;
 NOCRASH2 ;
 o = arrp (old,0,unsigned char) ;
 n = arrp (new,0,unsigned char) ;

 if (arrayMax (old) > arrayMax (new))
   messcrash ("Inconsistency in fuseBats") ;

 i = arrayMax (old) ; o--; n-- ;
 while (o++, n++, i--)
   if (*o)
     { if (*o & *n)
       messcrash ("Duplication in fuseBats") ;
     result = TRUE ;
     *n |= *o ;
     }
 OKCRASH ;
 return result ;
}

/**************************************************************/

BOOL  diskFuseBats (KEY fPlus, KEY fMinus, KEY sPlus, KEY sMinus,
		   int *nPlusp, int *nMinusp)
{
  BOOL done = FALSE, modif = FALSE ;
  Array fPlusArray = 0, fMinusArray = 0, 
    sPlusArray = 0, sMinusArray = 0;

  fPlusArray = arrayGet (fPlus,unsigned char,"c") ;
  fMinusArray = arrayGet (fMinus,unsigned char,"c") ;
  sPlusArray = arrayGet (sPlus,unsigned char,"c") ;
  sMinusArray = arrayGet (sMinus,unsigned char,"c") ;
  
  if (!sPlusArray || !sMinusArray ||
      !fPlusArray || !fMinusArray)
    goto abort ;
      
  modif |= fuseBats (fPlusArray, sPlusArray) ;
  modif |= fuseBats (fMinusArray, sMinusArray) ;

  if (modif)
    diskBATcompress (sPlusArray, sMinusArray) ;
  *nPlusp = diskBatSet (sPlusArray) ;
  *nMinusp = diskBatSet (sMinusArray) ;

  arrayStore (sPlus, sPlusArray,"c") ;
  arrayStore (sMinus, sMinusArray,"c") ;

  arrayKill (fPlus) ;  /* not valid if there were session branches */
  arrayKill (fMinus) ; 

  done = TRUE ;
 abort:
  arrayDestroy (sPlusArray) ;
  arrayDestroy (sMinusArray) ;
  arrayDestroy (fPlusArray) ;
  arrayDestroy (fMinusArray) ;

  return done ;
}

/**************************************************************/

void diskFuseOedipeBats (Array fPlusArray, Array fMinusArray) 
{ 
  BOOL modif = FALSE ;
    
  modif |= fuseBats (fPlusArray, plusArray) ;
  modif |= fuseBats (fMinusArray, minusArray) ;
  
  if (modif)
    diskBATcompress (plusArray, minusArray) ;
}

/********************************************/

/* Returns the number of set bits */
static int diskBatSet (Array a)
{
  register unsigned char *cp ;
  register unsigned int i;
  register int j = arrayExists (a) ? arrayMax (a) : 0 , n = 0 ;

  if (!j)
    return 0 ;
  cp  = arrp (a,0,unsigned char) ;
  while (j--)
    { if (!*cp) n+=8 ;
    else if (*cp != 0xff)
      for (i=1;i<256;i<<=1) if (~ (*cp) & i) n++;
    cp++;
    }
  return 8*arrayMax (a) - n ;
}

/********************************************/

/*returns the number of used and free blocks*/
void diskavail(unsigned long int *free, unsigned long int *used, unsigned long int *plus, unsigned long int *minus, unsigned long int *dr, unsigned long int* dw)
{
  char timeBuf[25] ;

  if (!batArray)
    diskBATread ();

  *used = diskBatSet (batArray) ;
  *plus = diskBatSet (plusArray) ;
  *minus = diskBatSet (minusArray) ;
  *free = 8*bigArrayMax (batArray) - *used ;
  *dr = blocksRead ; *dw = blocksWritten ;
  if (getenv ("ACEDB_DEBUG"))
    messdump ("diskavail %s: used %d, plus %d, minus %d, free %d\n", 
	      timeShow (timeParse ("now"), timeBuf, 25), *used, *plus, *minus, *free) ;
}

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
void diskalloc (BP bp)
{
  register unsigned int i ;
  /*unsigned because we use a right shift operator on i*/
  register int j ;
  register unsigned  char *cp , *top, *end , minus ;
  static DISK d = 0;
  
  if (bp->h.disk && 
     (bp->h.session == thisSession.session ||
       class (bp->h.key) == _VBat ||
       class (bp->h.key) == _VSession 
       ) )
    return;
  
  if (!batArray)
    diskBATread () ;
  /* We move the block on disk */
  if (class (bp->h.key) != _VGlobal) /* Freeing the global class in unconsistent */
    diskfree (bp) ;
  
 lao:
  j = arrayMax (batArray) ;   /* must come after eventual BATread */
  top = arrp (batArray,0,unsigned char) ;
  end = arrp (batArray,arrayMax (batArray) - 1 ,unsigned char) + 1 ;
  cp = arrp (batArray, (int) (d/8),unsigned char) ;
  if (cp >= end)
    cp = top ;
                
  while (j--) /* explore the whole BAT starting in the middle */
    {
      if (*cp != 0xff)
	for (i=1;i<256;i<<=1)
	  if (~ (*cp) & i)  /* free block on the global BAT */
	    { d = (DISK) (8* (cp - top)) ;
	    minus = arr (minusArray, (int) (d/8),unsigned char)  ;
	    if (~minus & i)	   /* Not freed in this session */
	      { *cp |= (unsigned char)i ;
	      arr (plusArray, (int) (d/8),unsigned char) |= i ; 
	      while (i>>=1) d++;
	      bp->h.disk=d ;
	      bp->h.session = thisSession.session ;
	      batIsModified = TRUE ;
	      return ;
	      }
	    else
	      { NOCRASH2 ; messcrash ("disk alloc inconsistency") ; }
	    }
      cp++;
      if (cp >= end) 
	cp = top;
    }
  
  diskExtend (FALSE);
  goto lao ;
} /* diskalloc */

/*************************************************/
/* do free only if allocated during this session */
/*************************************************/

void diskfree (BP bp)
{
  register unsigned char i;  /*unsigned because we use a right shift operator on i*/
  register unsigned char *cp , *minus, *plus ;
  register int j ;
  DISK d = bp->h.disk ;

  bp->h.disk = 0 ; /* prevents recursion */

  if (!d)
    return ;

  chrono ("diskfree") ;
  
  if (!batArray)
    diskBATread () ;
  
  j = arrayMax (batArray) ;   /* must come after eventual BATread */
  if (d >= 8*j)
    { NOCRASH2 ; messcrash ("diskfree called with a wrong disk address"); }

  cp = arrp (batArray, (int) (d/8),unsigned char) ;
  minus = arrp (minusArray, (int) (d/8),unsigned char) ;
  plus = arrp (plusArray, (int) (d/8),unsigned char) ;
  
  i = 1 << (d%8);
  if (! (*cp & i))
    { NOCRASH2 ; messcrash ("diskfree called on free block "); }

  /* on peut oublier les blocs crees pendant la session */
  if (*plus & i)  /* allocated during this session */
    { (*cp) &= ~i;     /* free it */
    (*plus) &= ~i;     /* free it */
    }
  else
    (*minus) |= i ;
  batIsModified = TRUE ;
  chronoReturn () ;

  return;
} /* diskfree */

/********************************************/


static BOOL isBat (DISK d)
{
  register unsigned char i;  /*unsigned because we use a right shift operator on i*/
  register unsigned char *cp ;
  register int j ;
  
  if (!batArray)
    return TRUE ; /* during init, one does a few direct address reads */
  
  j = arrayMax (batArray) ;   /* must come after eventual BATread */
  if (d >= 8*j)
    /* hors limite de la bat pour ce bloc */
    { NOCRASH2 ; messcrash ("diskfree called with a wrong disk address"); }

  cp = arrp (batArray, (int) (d/8),unsigned char) ;
  i = 1 << (d % 8);
  if (*cp & i)
    return TRUE ;
  return FALSE ;
}

/**************************************************************/
/* to write a block to disk   : return 0 if success */
/*         called by blockunload  */
void diskblockwrite (BP bp)
{
  char *vp= (char *)bp;
  long adrPartition;
  int partition = 0 ;
  unsigned int nb,size=BLOC_SIZE;
  myoff_t pos,pp ;
  long d = bp->h.disk ; 
  int fd, nretry = 0 ;
  NOCRASH1 ;
  
  if (!isWriteAccess ())
    return ;
  NOCRASH2 ;
  chrono ("diskblockwrite") ;
  blocksWritten++;
  
  if (!batArray)
    diskBATread () ;
  
  if ((d < 0)|| (d >= 8 * arrayMax (batArray)))
    messcrash ("Diskblockwrite called with wrong DISK address d=%ld, max = %d",
	       d, 8 * arrayMax (batArray)); 
  
  /* select partition and block address on partition : */
  adrPartition = partitionAddress (d, &partition);
  if (adrPartition < 0) 
    messcrash (" Bad Block Address in partition "); 
  /* open partition if not yet */
  fd = get_fd (partition, O_WRONLY);
  if (fd < 0)
    messcrash ("diskblockwrite can't open blocks %s (%s)",
	       partitionName (partition), messSysErrorText ()); 
  
  pos     = (myoff_t) (adrPartition * (long)size ) ;
  
  if (ps_debug1) 
    fprintf (stderr, "diskblockwrite: d %ld, part %d adr part %d offset %ld\n",
	     (unsigned long)d, (int)partition, (int)adrPartition, (unsigned long)pos);
  
#ifdef ACEDB4
extern void swapSuperBlock(BLOCKP bp) ;

  { static BLOCK copyBlock;
  if (swapData)
    { memcpy (&copyBlock, bp, sizeof (BLOCK));
    vp = (char *)&copyBlock;
    
    switch (bp->h.type) /* type must be invariant internal->external */
      { 
      case 'S':
	swapSuperBlock (&copyBlock);
	break;
	
      case 'A': case 'C':
	swapABlock (&copyBlock);
	break;
	
      case 'B':
	swapBBlock (&copyBlock);
	break;
	
      default:
	messcrash ("unknown block type trying to swap data"); 
      }
    }
  }
#endif /* ACEDB4  */
  
  if ((pp = lseek (fd, pos, SEEK_SET))!= pos)
    messcrash (
	       "Diskblockwrite with wrong DISK address d=%ld pos=%ld, returned=%ld (%s)\n",
	       (unsigned long)d, (unsigned long) pos, (unsigned long)pp, messSysErrorText ());
  
  if (euid != ruid)
    seteuid (euid);
  
  while ((nb = write (fd, vp, size)) != size || 
	 ((unsigned int) (* (DISK*)bp) == 0)
	 )
    {
      if (nretry++ == 12)
	messcrash ("Diskblockwrite  cannot actually write the DISK,n=%d (%s)\n",
		   nb, messSysErrorText ());
      sleep (2*nretry) ;
    }
  
  if (euid != ruid)
    seteuid (ruid);
  
  chronoReturn () ;
  OKCRASH ;

  return;
} /* diskblockwrite */

/**************************************************************/
/* to read a block from disk   : return 0 if success */
/*     called by blockload */

void diskblockread (BLOCK *bp, DISK d)
{
  char *vp= (char *)bp;
  long adrPartition;
  int partition = 0 ;
  unsigned int size=BLOC_SIZE;
  myoff_t pos, pos2 ;
  int n, nretry = 0, fd ;

  NOCRASH ;

  chrono ("diskblockread") ;
  blocksRead++ ;

  /* select partition and block address on partition : */
  adrPartition = partitionAddress (d, &partition);
  if (adrPartition < 0) 
    messcrash (" Bad Block Address in partition ");

  /* open partition if not yet */
  fd = get_fd (partition, O_RDONLY);
  if (fd < 0)
    messcrash ("Diskblockread : can't open blocks %s (%s)",
	      partitionName (partition), messSysErrorText ());
  
  pos = (myoff_t) (adrPartition * (long)size) ; 
  
  if (ps_debug1)
    fprintf (stderr, "Diskblockread : d %ld, part %d adr part %d offset %ld\n",
	    (long)d, (int)partition, (int)adrPartition, (long)pos);
  
  if ((pos2 = lseek (fd, pos, SEEK_SET)) != pos)
    messcrash (
	      "Diskblockread : wrong DISK address d=%ld partAdr %d pos=%ld, returned=%d (%s)\n",
	      (long)d, (long)pos, adrPartition, (long)pos2, messSysErrorText ());
  
  while ((n=read (fd,vp,size))!=size || 
	 ((unsigned int) (* (DISK*)bp) == 0)
	 )
   {
     if (nretry++ == 12)
       messcrash ("Diskblockread : cannot actually read the DISK, n=%d (%s)\n",
		  n, messSysErrorText ());
     sleep (2*nretry) ;
   }
  
#ifdef ACEDB4
  if (swapData)
    switch (bp->h.type) /* type must be invariant internal->external */
      { 
	
      case 'S':
	break;   /* Don't swap superblock, we need to read it before swapData
		    is valid */
      case 'A': case 'C':
	swapABlock (bp);
	break;
	
      case 'B':
	swapBBlock (bp);
	break;
	
      default:
	messcrash ("Diskblockread : unknown block type trying to transform data");
      }

#endif /* ACEDB4  */
  if (bp->h.type != 'S' && d != * (DISK*)bp)
    {
      const char *cp = messSysErrorText () ;
      if (! cp || ! *cp) cp = "syserror message is null" ;
      messcrash ("Diskblockread : read a block (bp->type=%c:addr=%u) not matching its address %u : syserror=%s"
		 , bp->h.type ? bp->h.type : '0' , (unsigned int) (* (DISK*)bp), (unsigned int) d, cp);
    }

  if (!isBat (d))
    messerror ("Diskblockread : read a block d=%dnot registered in the BAT"
	       , (unsigned int) d) ;

  chronoReturn () ; 
  OKCRASH ;

  return ;
} /* diskblockread */

/****************************************************/
/*   donne l'adresse dans une partition a partir de */
/*             l'adresse dans la base               */
/****************************************************/

/*
 * (babelfish says: give the address in a partition 
 * from the address in the base )
 *
 * d is the block number, *partition is filled in with
 * the partition number that block is in
 */

static long partitionAddress (DISK d, int *partition)
{
  /* premiere version tres bete : on cherche sequentiellement */
  long adresseDisque = (long)d;
  long adrDebutPartition = (long)0;
  long adrFinPartition = (long)0;
  PARTITION *pp = arrayp (pTable, 0, PARTITION) ;
  int part;
  for (part = 0; part <= lastUsedPartition; part++, pp++) 
    {
      adrFinPartition = adrDebutPartition + pp->currentBlocks ;
      if (d <  adrFinPartition)
	{
	  if (ps_debug1)
	    fprintf (stderr, "pour adr disque %ld adresse part : (%d,%ld)\n",
		    (long)d, part, adresseDisque - adrDebutPartition);
	  *partition = part;
	  return adresseDisque - adrDebutPartition;
	}
      adrDebutPartition += pp->currentBlocks ;
    }
  { 
    NOCRASH2 ;
    messcrash ("impossible disk address %uld, currently %d partitions lastAdr=%uld", 
	       (unsigned long) d, lastUsedPartition+1, adrDebutPartition) ;
  }

  return 0 ; /* for compiler happiness */
} /* PartitionAddress */

/****************************************************/
/*  give full partition pathname from #number       */
/****************************************************/

static char *partitionName (int num)
{ 
  PARTITION *pp ;
  static char buf[4096] ;

  if (num < 0 || num >= arrayMax (pTable)) 
    { NOCRASH2 ; messcrash ("partitionName : invalid number %d\n", num) ; }

  pp = arrayp (pTable, num, PARTITION) ;
  if (*pp->fileSystem)
    sprintf (buf,"%s/%s", pp->fileSystem, pp->fileName) ;
  else
    sprintf (buf, "%s/%s",  sessionFilName ("database", 0, "rd"), 
	     pp->fileName) ;
  return buf ;
}

static char *oldPartitionName (int num )
{ 
  PARTITION *pp ;
  static char buf[4096] ;

  if (num < 0 || num >= arrayMax (pTable)) 
    { NOCRASH2 ; messcrash ("partitionName : invalid number %d\n", num) ; }

  pp = arrayp (pTable, num, PARTITION) ;
  sprintf (buf,"%s/%s", COW_partition_dir, pp->fileName) ;
  return buf ;
}

/****************************************************/
/*  Close all open partitions of the DataBase       */
/****************************************************/

static void closeDataBase ()
{ 
  int x, fd, n_parts ;

  NOCRASH ;

  n_parts =  arrayMax (pTable) ;
  for (x=0; x<n_parts; x++)
    {
      fd = arrp (pTable, 0, PARTITION )->fd;
      if (fd < 0)
	continue;
#ifdef HAVE_FSYNC
      fsync (fd);
#endif
      close (fd);
      arrp (pTable, 0, PARTITION )->fd = -1;
    }

  OKCRASH ;
}

/****************************************************/
/*        Database map creation and update          */
/****************************************************/

/* called with argument "new"
 * if True, it's a new DataBase. We need to read de definition file and
 * create the map.
 * if False, it's just an extend control. In This case, we can just extend
 * the current map with new partitions added at the end of the DataBase 
 * description file.
 * Others modifications are ignored.
 */
static BOOL dataBaseReadDefinition (BOOL new, BOOL skip)
{
  char *defName ;
  FILE *fpDef;
  char *cp ;
  int level ;
  PARTITION *pp = 0 ;
  int nbPartitions = 0;
  NOCRASH1 ;

  /* DataBase description file name : */
  defName = sessionFilName ("wspec/database", "wrm", "r") ;

  if (ps_debug) fprintf (stderr, "opening def file %s\n", defName);

  /* open dataBase description file. */

  if (!defName || ! (fpDef = fopen (defName, "r")))
    { 
      /*
       * if there is no wspec/database.wrm, we make one up with a single
       * partition.
       */
      if (!new) return FALSE ;
      messdump ("No dataBase Definition, default to standard unix file");
      pp =  arrayp (pTable, 0, PARTITION) ;
      strcpy (pp->hostname, "local");
      strcpy (pp->fileSystem, "ACEDB");
      strcpy (pp->fileName, "blocks.wrm");
      pp->maxBlocks = DEFAULT_DATABASE_SIZE;
      pp->currentBlocks = 0;
      lastUsedPartition = -1;
      totalNbPartitions  = 1;
      return TRUE;
    }

  level = freesetfile (fpDef,"") ;
  if (!new && skip) 
    { int nbCurrentPart = totalNbPartitions ;
    /* there is already a dataBase and a map. Just add new partitions : */
    /* if any at the end of the dataBase description file :             */
      
    /* skip existing partitions : */
    while (nbCurrentPart)  
      if (freecard (level)) /* mhmp 17.06.99*/
	{ if (! (cp = freeword ()))
	  continue ; /* skip empty lines and comments */
	if (!strcmp (cp, "PART"))
	  nbCurrentPart --;
        if (strcmp (cp, "COW") == 0) 
	  {
	  cp = freeword ();
	  if (cp)
	    COW_partition_dir = strdup (cp);
	  else
	    messcrash ("COW needs directory name in wspec/database.wrm");
	  continue;
	}
	}
      else break ;
    /* if (!nbCurrentPart) return TRUE ;  dont return, mieg may99 */
    /* set nb Partitions to current number : */
    nbPartitions = totalNbPartitions;
    }

  /* create or extend partitions map : */
  NOCRASH2 ;
  while (freecard (level))
    { 
      int junk_int;
      if (! (cp = freeword ()))
	continue ; /* skip empty lines and comments */
      if (strcmp (cp, "COW") == 0) 
	{
	  cp = freeword ();
	  if (cp)
	    COW_partition_dir = strdup (cp);
	  else
	    messcrash ("COW needs directory name in wspec/database.wrm");
	continue;
	}
      if (strcmp (cp, "PART"))	/* only recognise PART keyword */
	continue ;
      pp =  arrayp (pTable, nbPartitions++, PARTITION) ;
      freestep (':') ;
      if (! freeint (&junk_int) ||
	  ! (cp = freeword ()) ||
	  ! (strcpy (pp->hostname, cp), cp = freeword ()) ||
	  ! (strcpy (pp->fileSystem, cp), cp = freeword ()) ||
	  ! (strcpy (pp->fileName, cp)) ||
	  ! freeint (&pp->maxBlocks) )
	messcrash ("Wrong line %d in file %s", freestreamline (level), defName) ; 
      
      /* current, maximum partition size in blocks and offset : */
      
      if (!strcmp ("ACEDB", pp->fileSystem)) *pp->fileSystem = 0 ;
      if (!strcmp ("NONE", pp->fileName))
	*pp->fileName = 0 ;
      else
	{ /* this partition name cannot match an existing one !!!! */
	  int i = nbPartitions - 1 ;
	  PARTITION *pp1 = pp - 1 ;
	  while (i--)
	    if (!strcmp (pp1->fileName, pp->fileName))
	      messcrash ("\nDataBase : partition %s exists !!", pp->fileName);

	  /* delete partition if exists , but not old style blocks.wrm !! */
	  if (strcmp ("blocks.wrm", pp->fileName))
	    unlink (pp->fileName) ;
	}
      pp->currentBlocks = 0;
    }
  /* complete memory map initialization : */
  if (new)
    lastUsedPartition = -1;
  totalNbPartitions  = nbPartitions;
  /* in case of "old DataBase style" map, if more than this "blocks.wrm" */
  /* partition, adjust max blocks.wrm partition to its current size      */
  pp =  arrayp (pTable, 0, PARTITION) ;
  if (!strcmp (pp->fileName, "blocks.wrm") &&
      nbPartitions > 1 )
    pp->maxBlocks = pp->currentBlocks ;
  /*   donot fclose (fpDef); allready done, srk, may 98 */
  OKCRASH ;
  return TRUE;
}



/*
 * dataBaseUpdateMap - write our in-memory copy of the list of partition
 * files to database/database.map
 */

static BOOL dataBaseUpdateMap (void)
{
  char *mapName;
  FILE *fil = 0 ;
  PARTITION *pp = arrayp (pTable, 0, PARTITION) ;

  char fileSystem[256], fileName[256];
  int nbPartitions = totalNbPartitions;

  if (euid != ruid)
    seteuid (euid); 

  mapName = sessionFilName ("database/database", "map", 0);
  /* sessionFilname fails on file which need write access via euid */
  if (ps_debug1) fprintf (stderr, "opening map file %s\n", mapName);
  if (!mapName || ! (fil = fopen (mapName, "w")))
    return FALSE ;

  fprintf (fil, "%d %d\n", 
	  totalNbPartitions, lastUsedPartition) ;
  while (nbPartitions--) 
    {    /* system and file names : */
      strcpy (fileSystem, pp->fileSystem);
      strcpy (fileName, pp->fileName);
      if (!strlen (pp->fileSystem)) strcpy (fileSystem, "ACEDB");
      if (!strlen (pp->fileName))   strcpy (fileName, "NONE");

      /*
       * this 1 is historical - it doesn't seem to mean anything, but
       * we keep it to be compatible with previous database files
       */
      if (fprintf (fil, "1 %s %s %s %d %d %d\n",
		  pp->hostname, fileSystem, fileName,
		  pp->maxBlocks, pp->currentBlocks, 0) == EOF) {
	messdump ("\nUnable to update dataBase map!!!!");
	if (ps_debug) fprintf (stderr, "error writing database map\n");
	return FALSE;
      }
      pp++;
    }
  fclose (fil);

  if (euid != ruid)
    seteuid (ruid);

  return TRUE;
} /* dataBaseUpdateMap */

/*************************************************/
/* download in memory structure the database map */
/*************************************************/

/****************************************************/
/* try to build a dataBase map with first partition */
/* as old style dataBase file database/blocks.wrm   */
/****************************************************/


/*
 * readPartitionMap - read the list of partion files from database/database.map
 *
 */

static void readPartitionMap (void)
{
  char *mapName, *cp ;
  FILE *fpMap = 0;
  PARTITION *pp ;
  int i, level = 0, currentNbPartitions, block_count;
  NOCRASH ;

  mapName = sessionFilName ("database/database", "map", "r") ;

  /* if no map file, attempt to create one with old dataBase file blocks.wrm*/
  if (!mapName || ! (fpMap = fopen (mapName, "r")))
    { 
      messcrash ("no database.map");
    }

  level = freesetfile (fpMap, "") ;
  if (!freecard (level))
    messcrash ("\nUnable to jump the top line of database/database.map");
  
  block_count = 0;
  i = 0 ; 
  lastUsedPartition = -1 ; 
  totalNbPartitions = 0 ;
  while (freecard (level))
    {
      int junk_int;

      cp = freepos () ;
      if (!*cp)
	continue ;
      
      /*
       * find the pTable entry for this partition
       */
      pp = arrayp (pTable, i, PARTITION)  ;
      totalNbPartitions++ ;
  
      /*
       * read the next line of the database.map file into it.
       */
      if (!freeint (&junk_int) ||
	  ! (cp = freeword ()) ||
	  (strcpy (pp->hostname, cp), ! (cp = freeword ())) ||
	  (strcpy (pp->fileSystem, cp), ! (cp = freeword ())) ||
	  (strcpy (pp->fileName, cp), !freeint (&pp->maxBlocks)) ||
	  !freeint (&pp->currentBlocks))
	messcrash ("\nUnable to line %d of map %s", i, mapName) ;

      /*
       * no, we have not opened the partition yet.
       */
      pp->fd = -1;

      /* system and file names : */
      if (!strcmp (pp->fileSystem, "ACEDB")) strcpy (pp->fileSystem, "") ;
      if (!strcmp (pp->fileName,    "NONE")) strcpy (pp->fileName, "") ;

      /*
       * and the block count is increased by the number of blocks in
       * this file
       */
      block_count += pp->currentBlocks;

      if (pp->currentBlocks)
 	{
	  /*
	   * if there is any data in this file, then it is our current best
	   * guess for which is the last partition used.
	   */
	  lastUsedPartition = i ;
	}

      i++ ;
    }

  /* 
   * Now we have all the partition lines from database/database.map.  But the
   * user may have added more to wspec/database.wrm so we also attempt to
   * find new definitions there.
   */
  currentNbPartitions = totalNbPartitions;
  if (dataBaseReadDefinition (FALSE, TRUE))
    if (totalNbPartitions > currentNbPartitions)
      dataBaseUpdateMap ();
  OKCRASH ;
}

/*************************************************************/

static void myText (char *cp, int x, int y)
{ 
  char ww[1000] ;
  int i = 990 ;
  char *cq = ww ;

  if (x < 0) x = 0 ;
  while (x--) *cq++ = ' ' ;
  *cq++ = '/' ;   *cq++ = '/' ; *cq++ = ' ' ;
  while (i-- && *cp) *cq++ = *cp++ ;
  *cq++ = '\n' ; *cq++ = 0 ;
  freeOut (ww) ;
}

/**************************************************************/

void disknewStatus (void)
{
  PARTITION *pp;
  int n, block_num, end, used, debut, plus, minus;
  char *s;

  pp = arrayp (pTable, 0, PARTITION);

  
  s= messprintf ("%-20s %8s %8s %10s  %3s %c %c %8s %8s %8s %8s",
		"file name", "max", "current", "debut", "fd", ' ', ' ',
		"used", "plus","minus", "free"); 
  myText (s, 0,0);

  n = arrayMax (pTable);
  block_num = 0;
  while (n-- > 0)
    {
      used = 0;
      plus = 0;
      minus = 0;
      end = pp->currentBlocks + block_num;
      debut = block_num;
      while (block_num < end)
	{
	  int c;
	  c = array (batArray, block_num / 8 , unsigned char);
	  while (c)
	    {
	      if (c & 1)
		used++;
	      c = c >> 1;
	    }
	  c = array (plusArray, block_num / 8 , unsigned char);
	  while (c)
	    {
	      if (c & 1)
		plus++;
	      c = c >> 1;
	    }
	  c = array (minusArray, block_num / 8 , unsigned char);
	  while (c)
	    {
	      if (c & 1)
		minus++;
	      c = c >> 1;
	    }
	  block_num += 8;
	}

      s = messprintf ("%-20s %8d %8d %10d  %3d %c %c %8d %8d %8d %8d",
		    pp->fileName, pp->maxBlocks, pp->currentBlocks, debut, pp->fd,
		    pp->fd_writeable ? 'w': 'r', pp->status_old ? 'O':' ', used, plus, minus, 
		    pp->currentBlocks - used - plus + minus);
      myText (s, 0, 0);
      pp++;
    }

}

#endif /* ACEDB 4 only */

/**************************************************************/
/**************************************************************/
 
