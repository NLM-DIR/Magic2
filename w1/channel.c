/*  File:  channel.c
 *  Author:  D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) D and J Thierry-Mieg 2015
 * -------------------------------------------------------------------
 * AceView is free software; you can redistribute it and/or
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
 * This file is part of the AceView project developped by
 *	 D and J Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *
 *  This should be linked with the acedb libraries
 *  available from http://www.aceview.org
 *
 *  Please direct all question about this code to
 *  mieg@ncbi.nlm.nih.gov
 *
 * Purpose of this code
 *  Implement a channel class, inspired by the 'go' language channels
 *  to allow robust synchronisation of a multithreaded program.
 *  Under the condition that channels should be the only mean of communication between threads,
 *  no other type of synchronization should be needed.
 *  Read the header file wh/channel.h for a detailled documentation and a user guide.
 *
 *  The public functions are
 *    channelCreate (int size, type, h), where type, like in arrays, can be int, float, .., struct x..., HIT, Array,...
 *       or any other declared type.
 *       Size is the max number of filled slots, over this limit channelPut will wait until a channelGet call
 *       consumes a slot.
 *    channelDestroy, or equivalently ac_free(channel) or ac_free(h), releases memory
 *    BOOL channelPut (channel, void *vp, type) : put a record into the channel, wait if channel is full
 *    BOOL channelGet (channel, void *vp, type) : get a record from the channel, wait if channel is emtpty
 *   (TYPE) channelGive (channel, void *vp, type) : get a record from the channel, wait if channel is emtpty
 *
 * Implementation notes (2026 revision)
 *   The records are stored in a ring buffer of cMax slots (head, count), protected by one mutex.
 *   Two condition variables: notEmpty (waited on by getters), notFull (waited on by putters).
 *   getWaiting/putWaiting count the blocked threads, so we only signal when somebody is waiting.
 *   A call that moves k records wakes k waiters (broadcast if k > 1) and, if records or room are
 *   left, one more waiter, so that no thread sleeps while it could work.
 *   Closing broadcasts both condition variables: blocked getters drain the queue then return 0,
 *   blocked putters return 0.
 *   Timed waits use CLOCK_MONOTONIC, they are not disturbed by a change of the wall clock.
*/

#include "ac.h"
#include "channel.h"
#include <pthread.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <time.h>

/*******************************************************************************************************/
/* the address of these variables identifies a channel / a lock,
 * this is cleaner than the address of a function, which cannot portably be stored in a void* */
static char chanMagicTag, lockMagicTag ;
#define CHAN_MAGIC ((void *)&chanMagicTag)
#define LOCK_MAGIC ((void *)&lockMagicTag)

#define CNAME(_ch) ((_ch)->c.title[0] ? (_ch)->c.title : "no-name")

/*******************************************************************************************************/
/* private helpers, all of them are called with the channel mutex locked */

static void chanCheck (CHAN *c, int size, const char *caller)
{
  if (! c || c->c.magic != CHAN_MAGIC)
    messcrash ("Invalid channel passed to %s", caller) ;
  if (c->c.size != size)
    messcrash ("Invalid type passed to %s, received size=%d expected size=%d", caller, size, c->c.size) ;
}

/* wake n waiters, we know that nobody is waiting if waiting == 0 */
static void chanWake (pthread_cond_t *cond, int waiting, int n)
{
  if (waiting > 0)
    {
      if (n > 1)
	pthread_cond_broadcast (cond) ;
      else
	pthread_cond_signal (cond) ;
    }
}

/* copy n <= count records out of the ring, in at most 2 memcpy */
static void chanRead (CHAN *c, void *vp, int n)
{
  size_t size = c->c.size ;
  int first = c->c.cMax - c->c.head ;

  if (first > n) first = n ;
  memcpy (vp, c->vp + (size_t)c->c.head * size, (size_t)first * size) ;
  if (n > first)
    memcpy ((unsigned char *)vp + (size_t)first * size, c->vp, (size_t)(n - first) * size) ;
  c->c.head += n ;
  if (c->c.head >= c->c.cMax) c->c.head -= c->c.cMax ;
  c->c.count -= n ;
  c->c.nGet += n ;
}

/* copy n <= cMax - count records into the ring, in at most 2 memcpy */
static void chanWrite (CHAN *c, const void *vp, int n)
{
  size_t size = c->c.size ;
  int tail = c->c.head + c->c.count ;
  int first ;

  if (tail >= c->c.cMax) tail -= c->c.cMax ;
  first = c->c.cMax - tail ;
  if (first > n) first = n ;
  memcpy (c->vp + (size_t)tail * size, vp, (size_t)first * size) ;
  if (n > first)
    memcpy (c->vp, (const unsigned char *)vp + (size_t)first * size, (size_t)(n - first) * size) ;
  c->c.count += n ;
  c->c.nPut += n ;
}

static void chanAfterRead (CHAN *c, int n)
{
  chanWake (&c->c.notFull, c->c.putWaiting, n) ;
  if (c->c.count > 0)    /* records left, another getter may work */
    chanWake (&c->c.notEmpty, c->c.getWaiting, 1) ;
}

static void chanAfterWrite (CHAN *c, int n)
{
  chanWake (&c->c.notEmpty, c->c.getWaiting, n) ;
  if (c->c.count < c->c.cMax)   /* room left, another putter may work */
    chanWake (&c->c.notFull, c->c.putWaiting, 1) ;
}

static void chanTrace (CHAN *c, const char *what, int n)
{
  fprintf (stderr, " --- %s %s %p : %d record(s), now %d in cache\n"
	   , what, CNAME (c), (void *)c, n, c->c.count) ;
}

/*******************************************************************************************************/

void channelDebug (CHAN *c, int debug, const char *title)
{
  if (c && c->c.magic == CHAN_MAGIC)
    {
      pthread_mutex_lock (&(c->c.mutex)) ;
      c->c.debug = debug ;
      if (title)
	{
	  strncpy (c->c.title, title, sizeof (c->c.title) - 1) ;
	  c->c.title[sizeof (c->c.title) - 1] = 0 ;
	}
      pthread_mutex_unlock (&(c->c.mutex)) ;
    }
} /* channelDebug */

/*******************************************************************************************************/
/* number of packets received by channel c, negative is channel is still open */
int channelCount (CHAN *c)
{
  long n = 0 ;

  if (c && c->c.magic == CHAN_MAGIC)
    {
      pthread_mutex_lock (&(c->c.mutex)) ;
      n = c->c.nPut ;
      if (n > INT_MAX) n = INT_MAX ;
      if (! c->c.isClosed)
	n = -n ;
      if (c->c.debug)
	fprintf (stderr, " - channelCount %s : closed=%d count=%ld\n", CNAME (c), c->c.isClosed ? 1 : 0, n) ;
      pthread_mutex_unlock (&(c->c.mutex)) ;
    }

  return (int) n ;
} /* channelCount */

/*******************************************************************************************************/
/* close the channel, must be called with the mutex locked */
static void chanCloseLocked (CHAN *c)
{
  if (! c->c.isClosed)
    {
      c->c.isClosed = TRUE ;
      if (c->c.debug)
	fprintf (stderr, " - Closing channel %s\n", CNAME (c)) ;
      /* getters drain the remaining records then return 0, putters return 0 */
      pthread_cond_broadcast (&(c->c.notEmpty)) ;
      pthread_cond_broadcast (&(c->c.notFull)) ;
    }
} /* chanCloseLocked */

void channelClose (CHAN *c)
{
  if (c && c->c.magic == CHAN_MAGIC)
    {
      pthread_mutex_lock (&(c->c.mutex)) ;
      chanCloseLocked (c) ;
      pthread_mutex_unlock (&(c->c.mutex)) ;
    }
} /* channelClose */

/*******************************************************************************************************/

void channelAddSources (CHAN *c, int nSources)
{
  if (c && c->c.magic == CHAN_MAGIC && nSources >= 0)
    {
      pthread_mutex_lock (&(c->c.mutex)) ;
      if (c->c.debug)
	fprintf (stderr, "channel %s add %d sources, new total %d\n"
		 , CNAME (c), nSources, c->c.nSources + nSources) ;
      if (! c->c.isClosed)   /* a channel closed unconditionally is closed for ever */
	c->c.nSources += nSources ;
      pthread_mutex_unlock (&(c->c.mutex)) ;
    }
} /* channelAddSources */

/*******************************************************************************************************/

void channelCloseSource (CHAN *c)
{
  if (c && c->c.magic == CHAN_MAGIC)
    {
      pthread_mutex_lock (&(c->c.mutex)) ;
      if (c->c.isClosed)
	{
	  /* already closed, for example by channelClose on an error path: nothing to do */
	  if (c->c.debug)
	    fprintf (stderr, "channel %s received closeSource but is already closed\n", CNAME (c)) ;
	}
      else
	{
	  if (c->c.nSources <= 0)
	    messcrash ("channel %s closeSource called beyond the number of added sources", CNAME (c)) ;
	  if (c->c.debug)
	    fprintf (stderr, "channel %s received closeSource, %d remaining sources\n"
		     , CNAME (c), c->c.nSources - 1) ;
	  if (--c->c.nSources == 0)
	    chanCloseLocked (c) ;
	}
      pthread_mutex_unlock (&(c->c.mutex)) ;
    }
} /* channelCloseSource */

/*******************************************************************************************************/
/* The caller must close the channel and make sure that no thread uses it any more
 * before calling ac_free: destroying a mutex or a condition variable with waiters is undefined */

static void uChannelDestroy (void *vp)
{
  CHAN *c = (CHAN *) vp ;

  if (c && c->c.magic == CHAN_MAGIC)
    {
      c->c.magic = NULL ;

      pthread_mutex_destroy (&(c->c.mutex)) ;
      pthread_cond_destroy (&(c->c.notFull)) ;
      pthread_cond_destroy (&(c->c.notEmpty)) ;

      ac_free (c->c.h) ;
    }
  else if (c && c->c.magic == LOCK_MAGIC)
    {
      c->c.magic = NULL ;

      pthread_mutex_destroy (&(c->c.mutex)) ;

      ac_free (c->c.h) ;
    }
  return ;
} /* uChannelDestroy */

/*******************************************************************************************************/

CHAN *uChannelCreate (int cMax, int size, AC_HANDLE h0)
{
  AC_HANDLE h ;
  CHAN *c ;
  pthread_condattr_t ca ;

  if (cMax <= 0)
    messcrash ("cMax = %d <= 0 in uChannelCreate", cMax) ;
  if (size <= 0)
    messcrash ("record size = %d <= 0 in uChannelCreate", size) ;
  if ((double)cMax * (double)size > (double)INT_MAX)
    messcrash ("a channel of %d records of %d bytes is too large in uChannelCreate", cMax, size) ;

  h = handleCreate () ;
  c = (CHAN *) halloc (sizeof (CHAN), h0) ;
  memset (c, 0, sizeof (CHAN)) ;
  c->c.h = h ;
  c->c.magic = CHAN_MAGIC ;
  c->c.cMax = cMax ;
  c->c.size = size ;
  pthread_mutex_init (&(c->c.mutex), NULL) ;

  /* timed waits must not depend on the wall clock */
  pthread_condattr_init (&ca) ;
  pthread_condattr_setclock (&ca, CLOCK_MONOTONIC) ;
  pthread_cond_init (&(c->c.notFull), &ca) ;
  pthread_cond_init (&(c->c.notEmpty), &ca) ;
  pthread_condattr_destroy (&ca) ;

  c->vp = halloc (cMax * size, h) ;
  /* register destroy */
  blockSetFinalise (c, uChannelDestroy) ;
  return c ;
} /* uChannelCreate */

/*******************************************************************************************************/
/* very light blocking interface
 * providing a semaphore locking to non threadsafe functions like dictAdd
 */
CHAN *channelLockCreate (AC_HANDLE h0)
{
  CHAN *c = (CHAN *) halloc (sizeof (CHAN), h0) ;

  memset (c, 0, sizeof (CHAN)) ;
  c->c.magic = LOCK_MAGIC ;
  pthread_mutex_init (&(c->c.mutex), NULL) ;
  /* register destroy */
  blockSetFinalise (c, uChannelDestroy) ;
  return c ;
} /* channelLockCreate */

void channelLock (void *vp)
{
  CHAN *c = (CHAN *)vp ;
  if (! c || c->c.magic != LOCK_MAGIC)
    messcrash ("Bad call to channelLock") ;
  pthread_mutex_lock (&(c->c.mutex)) ;
}

void channelUnlock (void *vp)
{
  CHAN *c = (CHAN *)vp ;
  if (! c || c->c.magic != LOCK_MAGIC)
    messcrash ("Bad call to channelUnlock") ;
  pthread_mutex_unlock (&(c->c.mutex)) ;
}

/*******************************************************************************************************/
/* Get up to max records.
 * wait == TRUE : blocks until at least one record is available
 * Returns:  n > 0 : number of records read,
 *           0 : the channel is closed and empty (vp is zeroed)
 *          -1 : (wait == FALSE only) the channel is empty
 */
int uChannelMultiGet (CHAN *c, void *vp, int size, int max, BOOL wait)
{
  int nn = 0, rc ;

  chanCheck (c, size, "uChannelMultiGet") ;
  if (max <= 0)
    return 0 ;

  pthread_mutex_lock (&(c->c.mutex)) ;

  while (c->c.count == 0)   /* the channel is empty */
    {
      if (c->c.isClosed)
	{
	  memset (vp, 0, size) ; /* a close channel returns a zeroed buffer */
	  pthread_mutex_unlock (&(c->c.mutex)) ;
	  return 0 ;
	}
      if (! wait)
	{
	  pthread_mutex_unlock (&(c->c.mutex)) ;
	  return -1 ;
	}
      /* wait till there is data, releasing the mutex while waiting. relocking when data is present */
      c->c.getWaiting++ ;
      rc = pthread_cond_wait (&(c->c.notEmpty), &(c->c.mutex)) ;
      c->c.getWaiting-- ;
      if (rc)
	messcrash ("uChannelMultiGet error %d returned by pthread_cond_wait on channel %s", rc, CNAME (c)) ;
    }

  nn = c->c.count < max ? c->c.count : max ;
  chanRead (c, vp, nn) ;
  if (c->c.debug)
    chanTrace (c, "uChannelMultiGet", nn) ;
  chanAfterRead (c, nn) ;

  pthread_mutex_unlock (&(c->c.mutex)) ;
  return nn ;
} /*  uChannelMultiGet */

/*******************************************************************************************************/

void *uChannelGet (CHAN *c, void *vp, int size)
{
  uChannelMultiGet (c, vp, size, 1, TRUE) ; /* blocking call */
  return vp ;          /* pseudo return by value, the macro will cast back to *(type *) */
} /* uChannelGet */

/*******************************************************************************************************/
/* channelGetMaxWait : blocking get with maximum wait time             */
/* Returns:                                                            */
/*    1 : success, data copied                                         */
/*    0 : channel closed and empty                                     */
/*   -1 : timeout reached, no data available                           */
/* ------------------------------------------------------------------- */
int uChannelGetMaxWait (CHAN *c, void *vp, int size, double maxSeconds)
{
  struct timespec ts ;
  BOOL timedOut = FALSE ;
  int rc ;

  chanCheck (c, size, "uChannelGetMaxWait") ;

  if (! (maxSeconds > 0))   /* also catches NaN */
    maxSeconds = 0 ;
  if (maxSeconds > 1e9)
    maxSeconds = 1e9 ;
  clock_gettime (CLOCK_MONOTONIC, &ts) ;
  {
    time_t sec = (time_t) maxSeconds ;
    ts.tv_sec += sec ;
    ts.tv_nsec += (long) ((maxSeconds - (double) sec) * 1e9) ;
    if (ts.tv_nsec >= 1000000000L)
      {
	ts.tv_sec++ ;
	ts.tv_nsec -= 1000000000L ;
      }
  }

  pthread_mutex_lock (&(c->c.mutex)) ;
  for (;;)
    {
      /* always look at the state while holding the mutex, also after a timeout */
      if (c->c.count > 0)
	{
	  chanRead (c, vp, 1) ;
	  if (c->c.debug)
	    chanTrace (c, "uChannelGetMaxWait", 1) ;
	  chanAfterRead (c, 1) ;
	  pthread_mutex_unlock (&(c->c.mutex)) ;
	  return 1 ;
	}
      if (c->c.isClosed)
	{
	  memset (vp, 0, size) ; /* a close channel returns a zeroed buffer */
	  pthread_mutex_unlock (&(c->c.mutex)) ;
	  return 0 ;
	}
      if (timedOut)
	{
	  pthread_mutex_unlock (&(c->c.mutex)) ;
	  return -1 ;
	}

      c->c.getWaiting++ ;
      rc = pthread_cond_timedwait (&(c->c.notEmpty), &(c->c.mutex), &ts) ;
      c->c.getWaiting-- ;
      if (rc == ETIMEDOUT)
	timedOut = TRUE ;
      else if (rc)
	messcrash ("uChannelGetMaxWait error %d returned by pthread_cond_timedwait on channel %s", rc, CNAME (c)) ;
      /* spurious wakeup or real signal: recheck the condition */
    }
} /* uChannelGetMaxWait */

/*******************************************************************************************************/
/* Put max records.
 * wait == TRUE : blocks until all max records are written, even if max is larger than the capacity
 * Returns:  n > 0 : number of records written,
 *                   n == max, unless the channel was closed while we were blocked
 *           0 : the channel is closed, nothing was written
 *          -1 : (wait == FALSE only) the channel is full, nothing was written
 * with wait == FALSE, as many records as there is room for are written
 */
int uChannelMultiPut (CHAN *c, void *vp, int size, int max, BOOL wait)
{
  int nn = 0, rc ;
  const unsigned char *src = (const unsigned char *) vp ;

  chanCheck (c, size, "uChannelMultiPut") ;
  if (max <= 0)
    return 0 ;

  pthread_mutex_lock (&(c->c.mutex)) ;

  while (nn < max)
    {
      int room, k ;

      /* checked under the mutex, and again after each wait: one cannot write on a close channel */
      if (c->c.isClosed)
	break ;

      room = c->c.cMax - c->c.count ;
      if (room == 0)   /* the channel is full */
	{
	  if (! wait)
	    break ;
	  /* wait till there is room, releasing the mutex while waiting. relocking when there is room */
	  c->c.putWaiting++ ;
	  rc = pthread_cond_wait (&(c->c.notFull), &(c->c.mutex)) ;
	  c->c.putWaiting-- ;
	  if (rc)
	    messcrash ("uChannelMultiPut error %d returned by pthread_cond_wait on channel %s", rc, CNAME (c)) ;
	  continue ;
	}

      k = max - nn < room ? max - nn : room ;
      chanWrite (c, src + (size_t)nn * size, k) ;
      nn += k ;
      if (c->c.debug)
	chanTrace (c, "uChannelMultiPut", k) ;
      chanAfterWrite (c, k) ;
      if (! wait)
	break ;
    }

  if (nn == 0)
    nn = c->c.isClosed ? 0 : -1 ;   /* nothing written: closed, or (non blocking only) full */

  pthread_mutex_unlock (&(c->c.mutex)) ;
  return nn ;
} /* uChannelMultiPut */

/*******************************************************************************************************/

void channelTest (int nn)
{
  AC_HANDLE h = handleCreate () ;
  CHAN *c = channelCreate (10, int, h) ;
  int i, g ;

  for (i = nn ; i >= 0 ; i--)
    channelTryPut (c, &i, int) ; /* this will stop putting after 10 loops */
  channelClose (c) ;

  i = 0 ;
  while (channelMultiGet (c, &g, 1, int) > 0)
    printf ("%d\t%d\n", i++, g) ;

  printf ("// This test is expected to print at most 10 lines, numbered [0,9], because the channel is created with cMax=10\n") ;
  ac_free (h) ;
  return ;
} /* channelTest */

/*******************************************************************************************************/
/*******************************************************************************************************/
/*******************************************************************************************************/
