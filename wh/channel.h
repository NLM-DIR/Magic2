/*  File: channel.h
 *  channels interface : derived from the GO language channels
 *  Author: Danielle and Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov) Mark Sienkiewicz (sienkiew@stsci.edu)
 *  Copyright (C) J Thierry-Mieg, 2015
 * -------------------------------------------------------------------
 * Created: mars 2015, mieg
 *   May 2015: added channelClose
 *   2026: revised, channelSelect and the remote_channel extension were removed
 *-------------------------------------------------------------------
 */

#ifndef CHANNEL_H
#define CHANNEL_H

/*************************************************************************/
/*************************************************************************/
/* channels : C implementation of the GO language channels
 *
 * Channels allow communication and ensure robust synchronization of multi-threaded programs.
 * The channels are mutex protected, so that the different threads can freely read/write
 * in the channels and should not need any other means of synchronization.
 * They are buffered and strongly typed: they can only transmit records of the declared type.
 * The channel handles the records as opaque buffers: only the size of the type is checked.
 *
 * CREATION:
 * To create a channel, provide the number of records to be buffered, a type, and a handle.
 *    Example 1 :
 *       CHAN *c = channelCreate (24, double, h) ; // this channel may buffer 24 doubles.
 * The type can be any basic C type (int, char*, double ...) or any typedefined struct.
 * Warning, since execution is parallelized, great care must be exercised if passing pointers.
 * Either the routine which PUTs a pointer in a channel should forget it and delegate
 * the responsability of the pointed buffer to the routine which will get the pointer from that
 * channel, or it should be insured that neither routines will modify that buffer.
 * The parameter h, a memory handle, is optional. It facilitates the destruction phase.
 *    Example 2 :
 *       typedef struct x3Struct { int x,y,z ;} X3 ;
 *       CHAN *cx3 = channelCreate (8, X3, h) ;  // cx3 may buffer 8 structures of type X3.
 * USAGE:
 * Data is written to the channel using channelPut, and read using channelGive or channelGet.
 * The Get/Put can occur in any order in different threads, as channels are self-synchronizing
 *
 * channelPut : writes a record of the correct type into a channel.
 *              parameters: CHAN c: an existing channel of the correct type
 *                          &z : the address of a variable of the correct type
 *                               the address is needed to construct the macro, but
 *                               in reality the record is copied and passed by value
 *                          type : the type of the variable
 *           // the type of the variable should match the type of the channel
 *     Example 3 :
 *       double pi = 3.1415926535 ;
 *       channelPut (c, &pi, double) ;
 *           // Execution blocks if the channel buffer is full,
 *           // Execution resumes when a record is consumed by another thread calling channelGet.
 *           // Returns 1 if the record was written, 0 if the channel is closed
 *           // (a channel closed while the call was blocked also returns 0, the record is not written)
 *
 * channelGive : reads a record of the correct type from a channel, returns the value.
 *              parameters: CHAN c: an existing channel of the correct type
 *                          &z : the address of a variable of the correct type
 *                               the call, if successful, will edit the variable
 *                          type : the type of the variable
 *           // the type of the variable should match the type of the channel
 *
 *   If the channel is empty, channelGive waits until another thread puts a record on the channel,
 *   as soon as the channel has data, it returns the next record and removes it from the channel.
 *   The record is copied at the address provided, and returned by value. As a result,
 *   channelGive may be used as the (blocking) argument of a function call as shown here:
 *      Example 4 :
 *        In the main program : create the channel cx3, as in example 2, and start two threads.
 *        In thread 1 ;
 *               X3 x3 ; // allocate an instance of struct x3Struct
 *               printf("%f\n", sinus(channelGive(cx3, &x3, X3).y) ; // blocking
 *        In thread 2 : // may happen before or after the request from thread 1
 *               X3 aa = { 0, pi/3, 0 } ; // allocate and initialise an instance of X3
 *               channelPut (cx3, &aa, X3) ; // write a record on the channel.
 *          At that point, thread 1 will resume and print the value 1/2.
 *          If the channel is closed and empty, channelGive returns a zeroed record.
 *
 * channelGet : reads a record of the correct type from a channel, returns the channel state.
 *              parameters: CHAN c: an existing channel of the correct type
 *                          &z : the address of a variable of the correct type
 *                               the call, if successful, will edit the variable
  *                          type : the type of the variable
 *           // the type of the variable should match the type of the channel
 *   If the channel is empty, channelGet waits until another thread puts a record on the channel,
 *   as soon as the channel has data, it returns the next record and removes it from the channel.
 *   The record is copied at the address provided and the returned value is TRUE.
 *   But if the channel is closed, the returned value is FALSE. As a result,
 *   channelGet may be used as the (blocking) argument of a loop control:
 *      Example 4b: // In this example we want to loop
 *        In thread 1 ;
 *               X3 x3  { 0, pi/3, 0 } ; // allocate and initialise an instance of X3
 *               for (i = 0; i < 10 ; i++)  // export 10 instances of x3,
 *                 channelPut (cx3, &x3, X3) ; // write a record on the channel.
 *               channelClose (cx3) ;          // close the channel
 *        In thread 2 : // may happen before or after the request from thread 1
 *               X3 x3 ; // allocate an instance of struct x3Struct
 *               while (channelGet (cx3, &x3, X3))  // write a record on the channel.
 *                  {} // do something with x3
 *               // the loop exits after consuming the 10 records
 *
 * A multi get/put version of these functions is available for programs using channels in tight
 * loops. Using these calls will improve the load balance between readers and writers by
 * accelerating the slower side relative to the faster side.
 *     Example 5 :
 *       int max = 7 ;
 *       double vp[max] ;
 *       n = channelMultiGet (c,vp,max,double) ;
 *           // Execution blocks if the channel buffer is empty
 *           // Execution resumes when at least one record is available
 *           // Returns the number of records (between 1 and max) copied at address vp
 *           // Returns zero if the channel is closed
 *       n = channelMultiPut (c,vp,max,double) ;
 *           // Execution blocks until all max records have been written
 *           // (even if max is larger than the capacity of the channel)
 *           // Returns max, or the number of records written before the channel was closed
 *
 * A non-blocking version of these functions is available for more complex programs:
 * channelTryPut (respectively channelTryGet) return 1 if the record was written (read)
 * or return -1 if the channel was full (empty). or 0 if the channel is closed.
 * This allows the program to loop and do something else.
 *
 * Waiting with a time limit:
 *     Example 6 :
 *       n = channelGetMaxWait (c, &x, double, 2.5) ;   // wait at most 2.5 seconds
 *           // n == 1 : a record was read
 *           // n == 0 : the channel is closed and empty
 *           // n == -1 : timeout, no data available
 *
 *  Closing a channel is the generic way to break a channelGet or channelMultiGet loop
 *  Once a channel is closed, it will not accept any new channelPut, i.e. any new records
 *  and when all previously queued records have been consumed, channelMultiGet will return 0
 *  Closing also wakes up the threads blocked in channelPut, which then return 0.
       Example 5b:
         In thread 1:
           channelClose (c) ;
         In thread 2:
           while (channelMultiGet (c, vp, max, double))
             {
             }
          print "channel c is closed, all records have been analyzed" ;
 *
 *  Parallel agents: In case several agents perform the same task in parallel.
 *  If two or more agents a1, a2, ... consume in parallel the records exported by channel inChannel
 *  At some point one of them, say a2, receive a inChannel isClosed message, while a1 is still running
 *  If a2 closes its outChannel, a1 will not be able to export to outChannel and data are lost
 *
 *  Solution : set and count sources, to insure that outChannel will accept all packets fro all its sources
       Example 7:
         Create n agents feeding out channel, and signal it
	   channelAddSources (outChannel, n) ; // increase source count
         In each agent a1, a2, ...
           while (channelGet (inChannel, ...)
	     { channelPut (outChannel) ;}
	   channelCloseSource (outChannel) ;  // decrease source count
     // outChannel will actually close when the number of sources decreases to zero
     // This way, only the last running agent will actually close the  outChannel
     //   However channelClose (outChannel) ; would still close unconditionnally
     //   channelCloseSource and channelAddSources on an already closed channel do nothing.

 * DESTRUCTION
 * To destroy a channel and recover the associated memory, call ac_free on c or h
 *     Example 8 :
 *       ac_free (c) ; // free c or alternatively
 *       ac_free (h) ; // free h, which will automatically call ac_free (c) ;
 *        remarks: -calling ac_free (c) multiple times is licit
 *                 -using c after its destruction is a fatal error, the program will crash,
 *          hopefully with an explicit error message, but it may also just seg-fault.
 *                 -close the channel and make sure that no thread is still using it
 *          before destroying it

 * USING CHANNELS in conjunction with the wego interface developped by Mark Sienkiewicz
 * Define a communication structure type containing the channels and anything else you want
      Example 9.a :
         typedef struct test1Struct { CHAN *c, *done ; } TEST1 ;
 * Define some functions which will run in parallel and communicate via the channels :
      example 9.b :
         static void wegoTest1 (const void *vp)
	   {
             const TEST1 *tt = (const TEST1 *)vp ;
             int x ;
	     BOOL b = FALSE ;

             while (channelGet (tt->c, &x, int))  // wait for data, break when tt->c is closed
	       { b = TRUE ; // do something with x until the channel is closed }
             channelPut (tt->done, &b, BOOL) ; // happens when tt->c is closed and 'do something with last x' is over
             return ;
           }
 * In the main program :
       example 9.c :
         int main ()
           {
              // declarations
              int i ;
              BOOL b = TRUE ;
              AC_HANDLE h = ac_new_handle () ; // allocate a memory handle
              TEST1 t ;          // declare the communication structure

              // initialisations
              wego_max_threads (4) ;  // Initialize wego and allow at most 4
	                              // simultaneously executing threads
                                      // the number of started threads is not limited
              t.c = channelCreate (6, int, h) ;     // initialize the channels
              t.done = channelCreate (1, BOOL, h) ;
              wego_go (wegoTest1, &t, TEST1) ;      // launch the thread

              // actual work
              for (i = 0 ; i < 10 ; i++)
	        channelPut (t.c, &i, int) ;   // send 10 integer values along channel t.c
                // we can send more than 6, since wegoTest1 is consuming the data in parallel

              // termination
              channelClose (t.c) ;     // signal the end of the workload
              channelGet (t.done, &b, BOOL) ;  // wait on the done channel
              printf ("done\n") ;              // work is done

	      // exit
	      ac_free (h) ;                    // release the memory
              return 0 ;                       // success
           }
 *
 * A verbose toggle is available to help during development.
 *     example 10:
 *              CHAN *c = channelCreate (4, int, h) ;
 *              channelDebug (c, 1, "title") ;   // log all calls to stderr
 *              channelDebug (c, 0, 0) ;         // stop logging
 *
 * A simple test is implemented in channelTest ().
 * Several complete examples are provided in wego/wego_test.c
 *
 */
/*************************************************************************/
/*************************************************************************/

typedef struct channelStruct CHAN ;

/* Public interface: please use these prototypes
 * they implement the interface discribed above
 *
 * TYPE stands for a declared C-type : int, float, double, CHAN ...

 CHAN *channelCreate (int cMax, TYPE, AC_HANDLE h) ;

 TYPE channelGive(CHAN *chan, TYPE *vp, TYPE) ;  // blocking, returns the value, updates *vp
 BOOL channelGet(CHAN *chan, TYPE *vp, TYPE)  ; // blocking, updates vp and returns TRUE, or FALSE if the channel is closed
 int  channelPut(CHAN *chan, TYPE *vp, TYPE) ;  // blocking, returns 1  if success, 0 if channel is closed

 int channelMultiGet(CHAN *chan, TYPE *vp, int max, TYPE) ;  // blocking, returns number of records read (1 to max), 0 if channel is closed
 int channelMultiPut(CHAN *chan, TYPE *vp, int max, TYPE) ; // blocking, returns max, or fewer if the channel was closed meanwhile

 int channelTryGet(CHAN *chan,  TYPE *vp, TYPE) ; // non-blocking, returns number read, 0 if channel is closed, -1 if empty
 int channelTryPut(CHAN *chan,  TYPE *vp, TYPE) ; // non-blocking, returns number written, 0 if channel is closed, -1 if full

 int channelGetMaxWait(CHAN *chan, TYPE *vp, TYPE, double maxSeconds) ; // 1 read, 0 closed, -1 timeout

 * Notice that these pseudo-prototypes are commented out
 * They are actually implemented by macros defined in channel_.h
 * This is not very elegant, but this is the only way we now to
 * pass a TYPE to a C function
 *
 * All channelPut function take as argument the address of a variable of type TYPE
 * to allow dynamic typing, but the argument is effectivelly copied and passed by value.
*/

void channelClose (CHAN *c) ;
void channelDebug (CHAN *c, int debug, const char *title) ;

int  channelCount (CHAN *c) ;  /* number of records received by channel c, negative if channel is still open */

void channelAddSources (CHAN *c, int nSources) ; /* source_count += nSources */
void channelCloseSource (CHAN *c) ; /* source_count -= 1, channel closes at zero */


void channelTest (int nn) ; /* This function should print on stdout the numbers 0, 1,2, 3... up to min(9,nn) */

CHAN *channelLockCreate (AC_HANDLE h) ; /* light interface to lock/unlock a mutex */
void channelLock (void *lock) ;
void channelUnlock (void *lock) ;

/*************************************************************************/
/* Private implementation of the uChannel functions */

#include "channel_.h"
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

#endif
