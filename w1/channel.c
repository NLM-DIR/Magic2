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
 *    channelDestroy, or equivalently ac_free(channel) or ac_free(h), closes the channel and releases memory
 *    BOOL channelPut (channel, void *vp, type, BOOL wait) : put a record into the channel, [wait] if channel is full 
 *    BOOL channelGet (channel, void *vp, type, BOOL wait) : get a record from the channel, [wait] if channel is emtpty
 *   (TYPE) channelGive (channel, void *vp, type, BOOL wait) : get a record from the channel, [wait] if channel is emtpty
*/

#include "ac.h"
#include "channel.h"
#include <pthread.h>

/*******************************************************************************************************/

void channelDebug (CHAN *c, int debug, char *title)
{
  if (c && ! c->c.isClosed)
    {
      pthread_mutex_lock (&(c->c.mutex)) ;
      c->c.debug = debug ;
      if (title)
	strncpy (c->c.title, title, 23) ;
      pthread_mutex_unlock (&(c->c.mutex)) ;
    }
} /* chanDebug */

/*******************************************************************************************************/
/* number of packets received by channel c, negative is channel is still open */
int channelCount (CHAN *c)
{
  int n = 0 ;

  if (c)
    {
      pthread_mutex_lock (&(c->c.mutex)) ;
      n = c->c.nPut1 ;
      if (! c->c.isClosed)
	n = -n ;
      if (c->c.debug)
	fprintf (stderr, " - channelCount %s : closed=%d count=%d\n", c->c.title[0] ? c->c.title : "no-name", c->c.isClosed ? 1 : 0, n) ;
      pthread_mutex_unlock (&(c->c.mutex)) ;
    }

  return n ;
} /* channelCount */

/*******************************************************************************************************/

void channelClose (CHAN *c)
{
  if (c)
    {
      pthread_mutex_lock (&(c->c.mutex)) ;
      if (c && ! c->c.isClosed /* && ! c->c.nLocks */ )
	{
	  c->c.isClosed = TRUE ;
	  if (c->c.debug)
	    fprintf (stderr, " - Closing channel %s\n", c->c.title[0] ? c->c.title : "no-name") ;
	  pthread_cond_signal(&(c->c.notEmpty)) ;
	}
      pthread_mutex_unlock (&(c->c.mutex)) ;
    }
} /* chanClose */

/*******************************************************************************************************/

void channelAddSources (CHAN *c, int nSources)
{
  if (c && nSources >= 0)
    {
      pthread_mutex_lock (&(c->c.mutex)) ;
      if (c->c.debug)
	fprintf (stderr, "channel %s add %d sources, new total %d\n"
		 , c->c.title[0] ? c->c.title : "no-name"
		 , nSources, c->c.nSources + nSources
		 ) ;
      if (! c->c.isClosed)
	c->c.nSources += nSources ;
      else
	messcrash ("channel %s Add %d Sources called on closed channel"
		 , c->c.title[0] ? c->c.title : "no-name"
		 , nSources
		 ) ;
      pthread_mutex_unlock (&(c->c.mutex)) ;
    }
} /* chanAddSources */

/*******************************************************************************************************/

void channelCloseSource (CHAN *c)
{
  if (c)
    {
      pthread_mutex_lock (&(c->c.mutex)) ;
      if (c->c.debug)
	fprintf (stderr, "channel %s received closeSources, %d remaining sources"
		 , c->c.title[0] ? c->c.title : "no-name"
		 , c->c.nSources - 1
		 ) ;
      if (! c->c.isClosed)
	c->c.nSources-- ;
      else
	messcrash ("channel %s closeSource called on closed channel"
		 , c->c.title[0] ? c->c.title : "no-name"
		 ) ;
      if (c->c.nSources < 0)
	messcrash ("channel %s closeSources called beyond the number of added sources"
		 , c->c.title[0] ? c->c.title : "no-name"
		 ) ;
      if (c->c.nSources == 0)
	{
	  c->c.isClosed = TRUE ;
	  if (c->c.debug)
	    fprintf (stderr, " - Closing channel %s\n"
		     , c->c.title[0] ? c->c.title : "no-name"
		     ) ;
	  pthread_cond_signal(&(c->c.notEmpty)) ;
	}
      pthread_mutex_unlock (&(c->c.mutex)) ;
    }
} /* chanCloseAt */

/*******************************************************************************************************/

static void uChannelDestroy (void *vp)
{
  CHAN *c = (CHAN *) vp ;
  if (c && c->c.magic == uChannelDestroy)
    {
      c->c.magic = NULL ;
 
      pthread_mutex_destroy (&(c->c.mutex)) ;
      pthread_cond_destroy (&(c->c.notFull)) ;
      pthread_cond_destroy (&(c->c.notEmpty)) ;

      ac_free (c->c.h) ;
    }
  if (c && c->c.magic == channelLockCreate)
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
  AC_HANDLE h = handleCreate () ;
  int nn = sizeof (CHAN1) + sizeof (CHANR) + (cMax+1) * size ;
  CHAN *c = (CHAN *) halloc (nn, h0) ;

  if (cMax <= 0)
    messcrash ("cMax = %d <= 0 in uChannelCreate", cMax) ;
  c->c.h = h ;
  c->c.cMax = cMax + 1 ;
  c->c.size = size ;
  c->c.magic = uChannelDestroy ;
  c->c.in = c->c.out = 0 ;
  pthread_mutex_init (&(c->c.mutex), NULL);
  pthread_cond_init (&(c->c.notFull), NULL);
  pthread_cond_init (&(c->c.notEmpty), NULL);
  c->vp = halloc ((cMax + 1) * size, h) ;
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
  int nn = sizeof (CHAN1) + sizeof (CHANR) ;
  CHAN *c = (CHAN *) halloc (nn, h0) ;

  c->c.magic = channelLockCreate ;
  pthread_mutex_init (&(c->c.mutex), NULL);
  /* register destroy */
  blockSetFinalise (c, uChannelDestroy) ;
  return c ;
}

void channelLock (void *vp)
{
  CHAN *c = (CHAN *)vp ;
  if (! c || c->c.magic != channelLockCreate)
    messcrash ("Bad call to channelLock") ;
  pthread_mutex_lock (&(c->c.mutex)) ;
}

void channelUnlock (void *vp)
{
  CHAN *c = (CHAN *)vp ;
  if (! c || c->c.magic != channelLockCreate)
    messcrash ("Bad call to channelUnlock") ;
  pthread_mutex_unlock (&(c->c.mutex)) ;
}

/*******************************************************************************************************/

int uChannelMultiGet (CHAN *c, void *vp, int size, int max, BOOL wait)
{
  int nn = 0, rc = 0 ;

  if (! c || c->c.magic != uChannelDestroy)
    messcrash ("Invalid channel passed to uChannelMultiGet") ;
  if (c->c.size != size)
    messcrash ("Invalid type passed to uChannelMutiGet, received size=%d expected size=%d", size, c->c.size) ;

  if (! wait) /* non blocking */
    {
      if (pthread_mutex_trylock (&(c->c.mutex)))
	return nn ;
    }
  else
    pthread_mutex_lock (&(c->c.mutex)) ;

  while (c->c.in == c->c.out && rc == 0)  /* the chammel is empty */
    {
      if (c->c.isClosed) /* return FALSE immediatly */
	{
	  memset (vp, 0, size) ; /* a close channel returns a zeroed buffer */
	  /* repeat the broadcast in case another client is waiting */
	  pthread_cond_signal(&(c->c.notEmpty)) ;
	  pthread_mutex_unlock (&(c->c.mutex)) ;
	  return 0 ;
	}
      else if (! wait)  /* return -1 immediatly */
	{
	  /* repeat the broadcast in case another client is waiting */
	  pthread_cond_signal(&(c->c.notEmpty)) ;
	  pthread_mutex_unlock (&(c->c.mutex)) ;
	  return -1 ;
	}
      else    /* wait till there is data, releasing the mutex while waiting. relocking when data is present */
	rc = pthread_cond_wait (&(c->c.notEmpty), &(c->c.mutex)) ;
    }
  if (rc)  /* error */
    {
      fprintf (stderr, " uChannelMultiGet error %d returned by  pthread_cond_wait ", rc) ;
      goto done ;
    }

	  /* actual work */
  while (max-- && c->c.in != c->c.out)
    {
      memcpy (vp + nn * c->c.size, c->vp + c->c.size * c->c.out, size) ;
      c->c.out = (c->c.out + 1) % c->c.cMax ;
      nn++ ;
    }
  c->c.nGet += nn ;  /* global counter */
  if (c->c.debug) 
    fprintf (stderr, " ---  uChannelMultiGet %s %ld  now %d in cache int:%d float:%f double:%g void:%p \n"
	     , c->c.title 
	     , (char *)c - (char *)0
	     , (c->c.in - c->c.out + 3* c->c.cMax) % c->c.cMax
	     , *(int*)vp, *(float*)vp, *(double*)vp, vp
	     ) ; 
 
  if (! 
      (((c->c.in + 1) % c->c.cMax) == c->c.out)  /* the channel is full */
      )
    {
      int i ;
      pthread_cond_signal (&(c->c.notFull)) ;
      for (i = 0 ; i < CHANNEL_SELECT_MAX ; i++)
      if (c->c.select[i])
	pthread_cond_broadcast (&((c->c.select[i])->c.notFull)) ;
    }
 done:
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
/*   >0 : success, data copied, number of records read (usually 1)    */
/*    0 : channel closed and empty                                     */
/*   -1 : timeout reached, no data available                           */
/* ------------------------------------------------------------------- */
int uChannelGetMaxWait(CHAN *c, void *vp, int size, double maxSeconds)
{
  if (! c || c->c.magic != uChannelDestroy)
    messcrash ("Invalid channel passed to uChannelGetMaxWait") ;

  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  time_t sec = (time_t)maxSeconds;
  long nsec = (long)((maxSeconds - sec) * 1e9) ;
  ts.tv_sec += sec;
  ts.tv_nsec += nsec;
  if (ts.tv_nsec >= 1000000000L)
    {
      ts.tv_sec++ ;
      ts.tv_nsec -= 1000000000L ;
    }
  
  pthread_mutex_lock (&c->c.mutex);
  
  while (TRUE)
    {
      // Check condition WHILE HOLDING the mutex (safe, no race)
      if (c->c.isClosed && c->c.in == c->c.out)
	{
	  memset (vp, 0, size) ; /* a close channel returns a zeroed buffer */
	  /* repeat the broadcast in case another client is waiting */
	  pthread_cond_signal(&(c->c.notEmpty)) ;
	  pthread_mutex_unlock (&(c->c.mutex)) ;
	  return 0 ;
	}
      
      if (c->c.in != c->c.out)
	{
	  /* data available, copy safely */
	  memcpy (vp, c->vp + c->c.size * c->c.out, size) ;
	  c->c.out = (c->c.out + 1) % c->c.cMax ;
	  c->c.nGet++ ;
	  if (! (((c->c.in + 1) % c->c.cMax) == c->c.out)  /* the channel is full */
	      )
	    {
	      int i ;
	      pthread_cond_signal (&(c->c.notFull)) ;
	      for (i = 0 ; i < CHANNEL_SELECT_MAX ; i++)
		if (c->c.select[i])
		  pthread_cond_broadcast (&((c->c.select[i])->c.notFull)) ;
	    }
	  pthread_mutex_unlock (&(c->c.mutex)) ;
	  return 1 ;
	}
      
      // empty and open → wait with timeout
      int ret = pthread_cond_timedwait(&(c->c.notEmpty), &(c->c.mutex), &ts);
      
      if (ret == ETIMEDOUT)
	{
	  pthread_mutex_unlock (&(c->c.mutex)) ;
	  return -1 ;
	}
      
        /*  spurious wakeup or real signal → recheck condition */
    }
  
  /* unreachable */
  pthread_mutex_unlock (&(c->c.mutex)) ;
  messcrash ("unreachable in uChannelGetMaxWait") ;
  return 0 ;
} /* channelGetMaxWait */

/*******************************************************************************************************/
/* return TRUE if one of the channel is probably ready
 * no garantee, since we relealse the mutex. the calling function must then use tryGet
 * on each channel
 */

BOOL channelSelect (CHAN **ccc)
{
  int ii, rc = 0 ;
  BOOL ok = FALSE , isOpen = FALSE ;
  CHAN *c, *c0 ;

  for (ii = 0 ; ccc[ii] && !ok ; ii++)
    {
      c = ccc[ii] ;
      pthread_mutex_lock (&(c->c.mutex)) ;
      if (c->c.in != c->c.out)  /* the chammel is empty */
	ok = TRUE ;
      if (! c->c.isClosed)
	isOpen = TRUE ;
      c0 = ccc[0] ;
      if (! ok && ii && c0 != c) /* ask for a braodcast to ccc[0] if ccc[iii] gets data */
	{
	  int i ;
	  for (i = 0 ; c0 && i < CHANNEL_SELECT_MAX ; i++)
	    {
	      if (c->c.select[i] == c0) /* already requested */
		{ c0 = 0 ; break ; }
	    }
	  for (i = 0 ; c0 && i < CHANNEL_SELECT_MAX ; i++)
	    {
	      if (c->c.select[i] == 0) /* already requested */
		{ c->c.select[i] = c0 ; c0 = 0 ; break ; }
	    }
	  if (c0)
	    messcrash ("More than CHANNEL_SELECT_MAX = %d channelSelec calls using different sets of channels all starting on channel %s\nPlease devise a less greedy algorithm or increase CHANNEL_SELECT_MAX in wh/channel_.h and recompile\n"
		       , CHANNEL_SELECT_MAX
		       , c0->c.title[0] ?  c0->c.title : "unknown"
		       ) ;
	}
      pthread_mutex_unlock (&(c->c.mutex)) ;
    }
  if (ok)   /* one of the channels has data */
    return TRUE ;
  if (! isOpen)
    return FALSE ; /* all channels are closed */
  
  c = ccc[0] ;
  pthread_mutex_lock (&(c->c.mutex)) ;
  rc = pthread_cond_wait (&(c->c.notEmpty), &(c->c.mutex)) ;
  pthread_mutex_unlock (&(c->c.mutex)) ;
  if (rc)  /* error */
    {
      fprintf (stderr, " channelSelect error %d returned by  pthread_cond_wait ", rc) ;
      return FALSE ;
    }

  return TRUE ;
} /* channelSelect */

/*******************************************************************************************************/

int uChannelMultiPut (CHAN *c, void *vp, int size, int max, BOOL wait)
{
  int nn = 0, rc = 0 ;

  if (! c || c->c.magic != uChannelDestroy)
    messcrash ("Invalid channel passed to uChannelGet") ;
  if (c->c.size != size)
    messcrash ("Invalid type passed to uChannelPut, received size=%d expected size=%d", size, c->c.size) ;
  if (c->c.isClosed)
    { 
      pthread_cond_signal(&(c->c.notEmpty)) ;
      return 0  ;  /* one cannot write on a close channel */
    }

  if (! wait) /* non blocking */
    {
      if (pthread_mutex_trylock (&(c->c.mutex)))
	  return -1 ;
    }
  else
    pthread_mutex_lock (&(c->c.mutex)) ;

  /* we always keep one slot empty in the circle otherwise at init time we would not be able to put */ 
  while (((c->c.in + 1) % c->c.cMax) == c->c.out && rc == 0)  /* the chammel is full */
    {
      if (! wait)  /* return FALSE immediatly */
	{
	  pthread_mutex_unlock (&(c->c.mutex)) ;
	  return FALSE ;
	}
      else    /* wait till there is data, releasing the mutex while waiting. relocking when data is present */
	rc = pthread_cond_wait (&(c->c.notFull), &(c->c.mutex)) ;
    }
  if (rc)  /* error */
    {
      fprintf (stderr, " uChannelMultiGet error %d returned by  pthread_cond_wait ", rc) ;
      goto done ;
    }

  	  /* actual work */
  c->c.nPut1 += max ;
  while (max-- && ((c->c.in + 1) % c->c.cMax) != c->c.out)
     {
       memcpy (c->vp + c->c.size * c->c.in, vp + nn * c->c.size, size) ;
       c->c.in = (c->c.in + 1) % c->c.cMax ;
       nn++ ;
     }
  c->c.nPut2 += nn ; /* global counter */
  if (c->c.debug) 
    fprintf (stderr, " ---  uChannelPut %s %ld now %d in cache int:%d float:%f double:%g  void*:%p\n"
	     , c->c.title 
	     , (char *)c - (char *)0
	     , (c->c.in - c->c.out + 3* c->c.cMax) % c->c.cMax 
	     , *(int*)vp, *(float*)vp, *(double*)vp, vp
	     ) ; 
   
  if (c->c.in != c->c.out)  
    pthread_cond_signal(&(c->c.notEmpty)) ;
 done:
  pthread_mutex_unlock (&(c->c.mutex)) ;

  return nn ;
} /* uChannelMultiPut */

/*******************************************************************************************************/

void channelTest (int nn)
{
  AC_HANDLE h = handleCreate () ;
  CHAN *c = channelCreate (10, int, h) ;
  int i ;
  BOOL g ;

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
