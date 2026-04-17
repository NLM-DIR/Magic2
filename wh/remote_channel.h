/*  File: remote_channel.h
 *  remote channels interface : 
 *     C implementation of the GO concurrent programming paradigm.
 *     This interface extends the channel interface from  multi-threading to multi-tasking
 *     and allows concurrent programming accross a computer network.
 *  Author: Danielle and Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov) Mark Sienkiewicz (sienkiew@stsci.edu)
 *  Copyright (C) J Thierry-Mieg, 2015
 * -------------------------------------------------------------------
 * Created: mars 2015, mieg
 *-------------------------------------------------------------------
 */

#ifndef REMOTE_CHANNEL_H
#define REMOTE_CHANNEL_H


/* This interface generalizes the channel interface from
 * multithreading to multitasking.
 * The main program would follow the same concurrent programming
 * design advocated in the GO language manuals, but this time
 * the tasks will execute as distinct processes, either on the
 * same computer or elswhere on the network, as specified
 * in the configuration file 'scripts/submit'
 * Remote channels follow the restricted publish-subscribe paradigm.
 * The server starts subtasks, publishes a channel interface, and the 
 * substasks subscribe. Subscription is signature protected and limited
 * to those subtasks created by the server.

 * CONFIGURATION:
 *   Client tasks  are launched via a call to system()
          system ("submit task ....") ;
     The submit script, available in wchannels/submit is self documented, and should be
     configured to fit your hardware (multicore machine, Sun Grid Engine compute farm...) 
     
 * INITIALISATION:
 *   Initialisation is mandatory and must preceed any other remote_channel calls. 
 *     TASK *t0 = taskServerInit (argv[0], TASK_LOCAL, AC_HANDLE h) ; // server side
 *     TASK *t0 = taskClientInit (&argc, argv, AC_HANDLE h) ;         // client side
 *   where the arguments of taskClientInit are those received by
 *                   main (int argc, const char *argv[]).
 *
 *   On the server side, the program name argv[0] will be  used below to construct
 *   robust remote client calls
 *   On the client side, taskClientInit will consume the arguments added implicitely by 
 *   taskDispatch (see below) and used by the remote_channel library to communicate with
 *   the server, and leave untouched the parameters specific to the application 
 *   (second argument of taskDispatch described below) 

 * CREATION:
 * To dispatch a task which will execute somewhere on the network
 * call taskDispatch()
 *   Example 1a:
          AC_HANDLE h = ac_new_handle () ;
          TASK *task = taskDispatch (0, "-a 12 -b hello -clientId 1", h) ;    // self-call 
          TASK *task = taskDispatch ("my_code", "-a 12 -b hello", h) ;
          TASK *task = taskDispatch ("/bin/grep", "foo_bar  my_file", h) ;
 *
 * The first parameter is the name of a program. If the name starts with /, as in /bin/grep, 
 * the name is not modified, otherwise the executable is searched first in the same directory 
 * where the server executable resides, then elsewhere on the path. Finally, is the executable 
 * is not found, the server stops and reports a FATAL ERROR.
 * These conventions may appear unusual, but they allows robust portable programming. Indeed, 
 * it is very likely that the server and the client program were developped at the same time 
 * and reside in the same directory, or have a known relative path, but once the program is
 * put in production, the absolute location of the executables cannot be known in advance. 
 * A convenient way to maintain agreement between the server and client programs, is
 * to create a single executable, and to recursivelly call the server executable, 
 * taskDispatch (0, ...) with modified arguments. This guarantees that server and 
 * client always remain synchrone during the development phase and, in particular, 
 * agree on the naming of the remote communication channels. 

 * The parameters, passed as the second argument, will be transmitted as is to the program.
 * The program will start executing at some unkown time. 
 * 
 * COMMUNICATION CHANNELS:
 * Synchrone communication with that taks can then be
 * established via channels, provided the 2 programs create pairs
 * of task reader/writer channels using the same channelName
 *    Example 2a:
 *         // In the program which called taskDispatch()
           typedef struct tStructA { .... } TTa ;
           typedef struct tStructB { .... } TTb ;
           CHAN *writer = taskCreateWriterChannel (t, "red", TTa, 12, h) ;
           CHAN *reader = taskCreateReaderChannel (t, "blue", TTb, 12, h) ;

            // In the program launched by taskDispatch
           typedef struct tStructA { .... } TTa ;
           typedef struct tStructB { .... } TTb ;
	   TASK *task = taskClientCreate () ;
           CHAN *reader = taskCreateReaderChannel (task, "red", TTa, 12, h) ;
           CHAN *writer = taskCreateWriterChannel (task, "blue", TTb, 12, h) ;

      // task zero refers to the parent task
      // it is the responsibility of the programmer to create matching pairs of channels
      // such that the same type of structure is written and read in a given channel



  * USAGE:
  * The channels can now be used exactly as in the examples given in channel.h
  * The only differnce is that the communication are established not between concurrent
  * threads, sharing the same memory space, but between distinct tasks, possibly running
  * on distinct machines. Thus, it would be meaningless for a remote channel to transmit
  * pointer types.
  * Latency may be highly variable, and it is probably a good idea to transmit rather 
  * large structure. On the other hand, on a large hardware architecture
  * like a Grid Engine, it may be reasonable to expect hundreds of tasks
  * to run simultaneously without contention.
  *
  * CAVEAT:
  *  Exactly as with usual channels, you wish to loop on
  *    while (channelMultiGet()) {}
  *  since the loop will break when the channel is closed
  *  i.e. when the last value 'put' in the channel will have been processed
  *  However in remote channels, it is essential to test and break on
  *    if (channelPut ()) { exit in some clean way }
  *  since because of the large asynchrone nature of multitasking
  *  it is very possible that the calling server has terminated (using another client)
  *  and is no longer listening on the client writer channel.
  *  In that case channelPut returns 0, and the client should terminate
  *
  * DESTRUCTION:
  * 
  * When the server exits, all its children, created by taskDispatch
  * will exit in the following few seconds, if they are already running
  * or are soon are they are dispatched, if they were still waiting.
  *
  * SWITCHING CONTEXT
  * A more subbtle way to kill the running clients is to switch context
  *
  *  tS = taskServerCreate (...) ;
  *  t1 = taskDispatch (...) ;
  *  taskSwitchContext (tS)  ;   // close communications with and kills t1
  *                        // and all previously opened remote channels 
  *  t2 = taskDispatch (...) ;
  *
  * This allows to run a first phase of a program, spawning many children tasks,
  * and when this phase is over, kill all these children tasks and start new ones.
  * A favorable side effect will be to end all the threads of the server
  * program which were waiting on readind these remote channels.
  *
  * Because of the very variable delays in a multitasking environment
  * the server side listening channel cannot reliably counts its closed providers.
  * The best strategy is probably to close the listening channel when the correct
  * number of data packets have been received. This will automatically
  * kill the the client writers.
  *
  *   for (i = 0 ; i < nn ; i++)
  *     channelPut (writer, &z, TYPE) ; // write nn packets
  *   n = 0 ;
  *   while (channelGet (reader, &z, TYPE))
  *      {  // do something with z
  *         if (++n == nn) break ;
  *       }
  *   channelClose (reader) ;  // indirectly kills all clients writing on reader
  *   // or more explicitly switchContext() 							     *
  */

/*************************************************************************/
/*************************************************************************/
/* Public interface: */
#include "channel.h"


typedef enum { TASK_LOCAL=0, TASK_SUBMIT } TASK_CONFIG ;

TASK *taskServerInit (const char *program_name, TASK_CONFIG config, AC_HANDLE h) ;
TASK *taskClientInit (int *argcp, const char **argv, AC_HANDLE h) ;
TASK *taskDispatch (const char *programName, char *params, AC_HANDLE h) ;
void taskSwitchContext (TASK *task) ; 

/* pseudo prototypes: TYPE stands for any C type: int, double, ... 
 * they are implemented via macros in remote_channel_.h

CHAN*  taskCreateReaderChannel (TASK *task, char *name, int max, TYPE) ;
CHAN*  taskCreateWriterChannel (TASK *task, char *name, int max, TYPE) ;

*/
/*************************************************************************/
/* Private implementation of the remote uChannel functions */

#include "remote_channel_.h"

/*************************************************************************/
/*************************************************************************/

#endif
