/*
************************************************************************
*** wego - an easy library for simple SMP/multithreading in C
*** Author Mark Sienkiewicz, Jean Thierry-Mieg, NCBI
************************************************************************
* $Header: /home/mieg/aa/CVS_ACEVIEW/ace/wego/wego.c,v 1.37 2020/06/29 01:20:13 mieg Exp $
*/

#ifdef SUN
#include <malloc.h>
#endif
#include <pthread.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>

#include "wego_.h"
#include "ac.h"
#include "channel.h"

#define wego_debug 0

#if wego_debug
#define DEBUG(_x) { fprintf _x ; }
#else
#define DEBUG(_x) { }
#endif

/************************************************************************/

/*
* information about a single task
*/
struct wego_task
{
  void	(*function)();		/* function to call */
  void *	parameter_struct;	/* parameters from user */
  int	parameter_struct_size;
  void *	result_struct;		/* place to store return values */
  int	result_struct_size;
  
  
  int	flags;
  int lightInterface ;

#define WEGO_TASK_ALLOC_RESULT 	1
  
  struct wego_group *group;	/* what group are we in? */
  
  struct wego_task *q_next;	/* next in queue of all pending tasks */
  
  struct wego_task *g_next;	/* next in list of all tasks in group */

  /*
   * fields for applications waiting for our result
	*/
  pthread_mutex_t		finished_semaphore;
  pthread_cond_t		finished_condition;
  int			finished_flag;
  int			reportable_flag;
  int			ever_reported;
  int			running_flag;
  
};

/************************************************************************/

/*
* group of tasks
*/
struct wego_group
	{
	pthread_mutex_t		group_semaphore;
	pthread_cond_t		group_condition;
	struct wego_task *task_list;
	int	task_count;
	};


/************************************************************************/

/*
* Thread management variables
* You must have thread_take semaphore to access any of these variables:
*/

	/*
	* thread synchronization for these variables
	*/
	static pthread_cond_t	thread_take_condition = PTHREAD_COND_INITIALIZER;
	static pthread_mutex_t	thread_take = PTHREAD_MUTEX_INITIALIZER;

	struct threadrunner
		{
		int			main;		/* main program */
		int			running;	/* using a cpu */
		int			available;	/* available to run a task */
		int			waiting;	/* waiting to continue */
		int			thread_id;	/* small integer thread number */
		pthread_cond_t		condition;
		pthread_mutex_t		semaphore;
		struct threadrunner *	next;
		};

	static struct threadrunner 	main_thread = { 1, 1, 0, 0, 0 };
	static struct threadrunner *	all_threads = &main_thread;

	static pthread_cond_t		thread_startup_condition = PTHREAD_COND_INITIALIZER;
	static pthread_mutex_t		thread_startup_lock = PTHREAD_MUTEX_INITIALIZER;

	/*
	* queue of all tasks waiting to be run
	*/
	static struct wego_task *queued_tasks = NULL;

/************************************************************************/

/*
* forward declarations
*/
static void wego_want_to_sleep();
static void wego_want_to_run();
static void wego_schedule();
static void *wego_threadrunner(void *thread_parameter);



/************************************************************************/

pthread_key_t threadrunner_key;

#define TID ( ((struct threadrunner *)pthread_getspecific(threadrunner_key))->thread_id )

/************************************************************************/

/*
* This is the maximum number of threads that can be running concurrently.
* You can set this from your main program before you call any wego
* functions.  If you do not, we use the value from the environment variable 
* WEGO_THREADS.  Otherwise, the default is 8.
*
* If you change this value after starting some threads, the results are
* undefined.
*/
static int max_threads = 0 ;

/************************************************************************/
/************************************************************************/
/*
 * Create a new task.
 *
 ** There are no concurrency issues until we enter the task into the
 ** common data structures.
 */

struct wego_task *wego_run (struct wego_group *group, void (*function)(), 
			    void *parameter_struct, size_t parameter_struct_size, 
			    void *result_struct, size_t result_struct_size
			    )
{
  struct wego_task *r;
  
  DEBUG((stderr, "wego_run called from %d\n",TID))
    
    /*
     * create the struct
     */
    r = malloc(sizeof(struct wego_task));
  if (!r)
    return NULL;
  
  /*
   * all fields initialized to 0
   */
  memset(r,0,sizeof(*r));
  
  pthread_mutex_init(&r->finished_semaphore, NULL);
  pthread_cond_init(&r->finished_condition, NULL);
  
  r->function = function;
  
  if (parameter_struct)
    {
      r->parameter_struct = malloc(parameter_struct_size);
      if (! r->parameter_struct)
	{
	  free(r);
	  return NULL;
	}
      memcpy(r->parameter_struct,parameter_struct,parameter_struct_size);
    }
  else
    r->parameter_struct = NULL;
  r->parameter_struct_size = parameter_struct_size;
  
  if (result_struct)
    r->result_struct = result_struct;
  else
    {
      r->result_struct = malloc(result_struct_size);
      if (!r->result_struct)
	{
	  free(r->parameter_struct);
	  free(r);
	  return NULL;
	}
      r->flags |= WEGO_TASK_ALLOC_RESULT;
    }
  r->result_struct_size = result_struct_size;
  
  r->group = group;
  
  /*
   * add it to the list of everything to be done
   */
  pthread_mutex_lock(&thread_take);
  r->q_next = queued_tasks;
  queued_tasks = r;
  pthread_mutex_unlock(&thread_take);
  
  /*
   * if it is part of a group, add it to the group.
   */
  if (group)
    {
      pthread_mutex_lock(&group->group_semaphore);
      group->task_count++;
      r->g_next = group->task_list;
      group->task_list = r;
      pthread_mutex_unlock(&group->group_semaphore);
    }
  
  DEBUG((stderr, "wego_run: calling schedule -\n"))
    wego_schedule();
  return r;
} /* wego_run */

/************************************************************************/
/* Ligth interface expecting all communications to be handled via channels
* Create a new task.
*
** There are no concurrency issues until we enter the task into the
** common data structures.
*/

WEGO_TASK *wego_go_size (void(*f)(const void *vp), const void *vp, int size, struct wego_group *group) 
{
  struct wego_task *r;
  /* int parameter_struct_size = sizeof (void*) ;   */
  DEBUG((stderr, "wego_go_size: start in %d\n",TID))
    
    /*
     * create the struct
     */
    r = malloc(sizeof(struct wego_task));
  if (!r)
    return NULL;
  
  /*
   * all fields initialized to 0
   */
  memset(r,0,sizeof(*r));
  
  pthread_mutex_init(&r->finished_semaphore, NULL);
  pthread_cond_init(&r->finished_condition, NULL);
  
  r->function = f ;
  r->lightInterface  = 1 ;

  r->parameter_struct = malloc(size) ;
  r->parameter_struct_size = size ;
  memcpy(r->parameter_struct, vp, size) ;
  
  r->result_struct = 0 ;
  r->result_struct_size = 0 ;
  
  r->group = group ;
  
  /*
   * add it to the list of everything to be done
   */
  pthread_mutex_lock(&thread_take);
  r->q_next = queued_tasks;
  queued_tasks = r;
  pthread_mutex_unlock(&thread_take);
  
  /*
   * if it is part of a group, add it to the group.
   */
  if (group)
    {
      pthread_mutex_lock(&group->group_semaphore);
      group->task_count++;
      r->g_next = group->task_list;
      group->task_list = r;
      pthread_mutex_unlock(&group->group_semaphore);
    }
 
  DEBUG((stderr, "wego_go_size: - schedule -\n"))
    wego_schedule();
  DEBUG((stderr, "wego_go_size: - done -\n"))
  return r;
} /* wego_go_size */

void wego_go_ (void(*f)(const void *vp), const void *vp, int size)
{
  if (! max_threads)   /* a security if wego has not been initialized */
    wego_max_threads (8) ;
  wego_go_size(f, vp, size, 0) ;
} /* wego_go_ */

void wego_exit (int err) 
{ 
  pthread_exit(0) ;
} /* wego_exit */

/*************************************************************************************/
/*************************************************************************************/
static CHAN *wego_log_channel = 0 ;

void wego_log (const char *text)
{
  char *cp = strnew (text, 0) ;
  channelPut (wego_log_channel, &cp, const char**) ;
} /* wego_log  */

/************************************************************************/

static void wego_export_logs (const void *vp)
{
  const char *text = 0 ;
  while (channelMultiGet (wego_log_channel, &text, 1, const char*))
    if (text)
      {
	fprintf (stderr, "// %s: %s\n", timeShowNow (), text) ;
	ac_free (text) ;
      }
} /* wego_export_logs */

/************************************************************************/

void wego_flush (void)
{
  channelClose (wego_log_channel) ;
  sleep (1) ; /* to allow closure of wego_log */
} /* wego_flush */

/*************************************************************************************/
/*************************************************************************************/
/*
 * create a group
 *
 ** There are no concurrency issues until after we return the data structure.
 */

struct wego_group *wego_new_group (void)
{
  struct wego_group *g;
  g = malloc(sizeof(*g));
  if (!g)
    return NULL;
  memset (g,0,sizeof(*g));
  pthread_mutex_init (&g->group_semaphore, NULL);
  pthread_cond_init (&g->group_condition, NULL);
  DEBUG((stderr, "wego_new_group return %v\n",g)) ;
  return g;
} /* wego_new_group */

/************************************************************************/
/*
* extract the result of a task
*
* t->finished_semaphore controls access to the data about whether a task
* is finished.
*	t->finished_flag means it is finished
*	t->reportable means it is eligible to be reported by wego_fin()
*	t->ever_reported means the results have been reported at least once
*
*/

void *wego_result (struct wego_task *t)
{
  void *r;
  
  pthread_mutex_lock(&t->finished_semaphore);
  DEBUG((stderr, "wego_result called from thread %d\n",TID))
    while (! t->finished_flag)
      {
	wego_want_to_sleep();
	pthread_cond_wait( &t->finished_condition, &t->finished_semaphore );
	wego_want_to_run();
      }
  t->reportable_flag = 0;
  t->ever_reported = 1;
  r = t->result_struct;
  DEBUG((stderr, "wego_result done\n"))
    pthread_mutex_unlock(&t->finished_semaphore);
  return r;
} /* wego_result */

/************************************************************************/

struct wego_task *wego_fin (struct wego_group *group)
{
  struct wego_task *t;
  int reported_already;
  
  DEBUG((stderr, "wego_fin start\n")) ;
  
  reported_already=0;
  
  /*
   * We have to lock the group semaphore while we examine the list
   * of all tasks for one we have not reported on yet.
   *
   * We have to lock the finished semaphore while we examine an
   * individual task.
   *
   * There is no deadlock because this is the only place in the
   * library where we lock BOTH a group semaphore and a finished semaphore.
   */
  
  for (;;)
    {
      pthread_mutex_lock(&group->group_semaphore);
      /*
       * Look for any one task that is reportable.
       */
      for (t=group->task_list; t; t=t->g_next)
	{
	  pthread_mutex_lock(&t->finished_semaphore);
	  if (t->reportable_flag)
	    {
	      /*
	       * found one - report it.
	       */
	      t->reportable_flag = 0;
	      t->ever_reported = 1;
	      pthread_mutex_unlock(&t->finished_semaphore);
	      pthread_mutex_unlock(&group->group_semaphore);
	      DEBUG((stderr, "wego_fin return data\n"))
		return t;
	    }
	  if (t->ever_reported)
	    reported_already++;
	  pthread_mutex_unlock(&t->finished_semaphore);
	}
      
      /*
       * If all tasks on the list have been reported, return NULL.
       */
      if (group->task_count == reported_already)
	{
	  pthread_mutex_unlock(&group->group_semaphore);
	  DEBUG((stderr, "wego_fin return NULL\n"))
	    wego_want_to_run();
			return NULL;
	}
      pthread_mutex_unlock(&group->group_semaphore);
      
      
      /*
       * didn't find any reportable. sleep and try again.
       */
      wego_want_to_sleep();
      pthread_mutex_lock(&group->group_semaphore);
      pthread_cond_wait( &group->group_condition, &group->group_semaphore );
      pthread_mutex_unlock(&group->group_semaphore);
      wego_want_to_run();
      
    }
} /* wego_fin */

/************************************************************************/

static void wego_want_to_sleep (void)
{
  struct threadrunner *t;
  pthread_mutex_lock(&thread_take);
  DEBUG((stderr, "want_to_sleep called from thread %d\n",TID))
    t = pthread_getspecific(threadrunner_key);
  t->running = 0;
  t->waiting = 1;
  pthread_mutex_unlock(&thread_take);
  DEBUG((stderr, "want_to_sleep: calling schedule \n"))
    wego_schedule();
} /* wego_want_to_sleep */

/************************************************************************/

static void wego_want_to_run (void)
{
  struct threadrunner *t;
  
  t = pthread_getspecific(threadrunner_key);
  DEBUG((stderr, "wego_want_to_run from thread %d\n",TID))
    
    pthread_mutex_lock(&thread_take);
  while (! t->running)
    {
      DEBUG((stderr, "wego_want_to_run in %d waiting\n",TID))
	pthread_cond_wait(&thread_take_condition, &thread_take);
    }
  DEBUG((stderr, "wego_want_to_run in %d success running\n",TID))
    t->waiting = 0;
  pthread_mutex_unlock(&thread_take);
} /* wego_want_to_run */

/************************************************************************/
/*
 * This is the main program of each thread.
 */

static void *wego_threadrunner (void *thread_parameter)
{
  struct threadrunner *t;
  
  pthread_mutex_lock(&thread_startup_lock);
  t = thread_parameter;
  pthread_setspecific(threadrunner_key, thread_parameter);
  
  t->running = 0;
  t->available = 1;
  
  pthread_cond_signal(&thread_startup_condition);
  pthread_mutex_unlock(&thread_startup_lock);
  
  
  DEBUG((stderr, "threadrunner in %d: startup complete\n",TID)) ;
  
  /*
   * Wait for: we can run a task
   */
  for (;;)
    {
      struct wego_task *w;
      struct wego_group *g;
      
      DEBUG((stderr, "					threadrunner top of loop\n"))
	
	/*
	 * set our thread status to show that:
	 * - we are not using CPU
	 * - we are not waiting to run
	 * - we are available to run a task
	 */
	pthread_mutex_lock(&thread_take);
      t->running = 0;
      t->available = 1;
      t->waiting = 0;
      pthread_mutex_unlock(&thread_take);
      
      /*
       * make sure the scheduler thinks about our new state
       */
      wego_schedule();
      
      /*
       * sleep until the scheduler marks us as running.  the scheduler
       * set running to 1, available to 0.
       */
      pthread_mutex_lock(&t->semaphore);
      while (t->running == 0)
	pthread_cond_wait(&t->condition, &t->semaphore);
      pthread_mutex_unlock(&t->semaphore);
      
      /*
       * fetch a task from the queue.  it is possible that the
       * queue may be exhausted now and then.  If so, go back to sleep.
       */
      pthread_mutex_lock(&thread_take);
      w = queued_tasks;
      if (!w)
	{
	  DEBUG((stderr, "					threadrunner no more tasks\n"))
	    pthread_mutex_unlock(&thread_take);
	  continue;
	}
      queued_tasks = queued_tasks -> q_next;
      pthread_mutex_unlock(&thread_take);
      
      /*
       * execute the user's function
       */
      /*
	if (w->lightInterface)
	w->function(w->parameter_struct) ;
      else
	w->function(w, w->parameter_struct, w->result_struct);
      */
      
      if (w->lightInterface)
	((void(*)(const void*))w->function)(w->parameter_struct);
      else
	((void(*)(const void*, void*, void*))w->function)(w, w->parameter_struct, w->result_struct);
      /*
       * Flag the task as finished.
       * Wake up anybody who is waiting for the results.
       */
      pthread_mutex_lock(&w->finished_semaphore);
      w->finished_flag = 1;
      w->reportable_flag = 1;
      pthread_cond_signal(&w->finished_condition);
      pthread_mutex_unlock(&w->finished_semaphore);
      
      /*
       * If the task was part of a group, wake up anybody
       * who is waiting for group results.
       */
      g=w->group;
      if (g)
	{
	  pthread_mutex_lock(&g->group_semaphore);
	  pthread_cond_signal(&g->group_condition);
	  pthread_mutex_unlock(&g->group_semaphore);
	}
      
    }
  /* 
   * We never get here.  
   * Threads never go away.  If we don't need a thread, it
   * just sits idle.
   */
  return NULL ; 
} /* wego_threadrunner */

/************************************************************************/

static int wego_new_thread (void)
{
  pthread_t thread;
  struct threadrunner *t;
  static int thread_id = 1;
  int error = 0 ;
  int ok = 0 ;
  
  DEBUG((stderr, "wego_new_thread in %d\n",TID)) ;
  
  t = malloc(sizeof(*t));
  memset(t,0,sizeof(*t));
  t->running = 0;
  t->available = 1;
  t->main = 0;
  t->thread_id = thread_id++;
  t->next = all_threads;
  pthread_cond_init(&t->condition,NULL);
  pthread_mutex_init(&t->semaphore,NULL);
  
  pthread_mutex_lock(&thread_startup_lock);
  error = pthread_create(&thread, NULL, wego_threadrunner, (void *)t) ;
  if (! error)
    {
      pthread_cond_wait(&thread_startup_condition, &thread_startup_lock);
      all_threads = t;
      ok = 1 ;
    }
  else /* creation failed */
    {
      fprintf (stderr, "ERROR  wego_new_thread: pthread_create retuned error %dn", error) ;
      t->available = -1 ;  
      free (t) ;
    }
  pthread_mutex_unlock(&thread_startup_lock);
  return ok ;
} /* wego_new_thread */

/************************************************************************/
/*
 * this function is used by wego_schedule() to wake up a specific
 * thread that is waiting in threadrunner().
 */

void wego_wakeup_threadrunner (struct threadrunner *t)
{
  pthread_mutex_lock(&t->semaphore);
  pthread_cond_signal(&t->condition);
  pthread_mutex_unlock(&t->semaphore);
} /* wego_wakeup_threadrunner */

/************************************************************************/

static void wego_schedule (void)
{
  struct threadrunner *t;
  int running = 0;
  int num_available = 0;
  
  pthread_mutex_lock(&thread_take);
  DEBUG((stderr, "\n\nwego_schedule in %d, max_treads = %d\n",TID, max_threads)) ;    
  /*
   * See how many threads are running, how many are available.  I consider
   * counting them here to be more reliable than ++/-- all over the place.
   * Typically, I expect a small number of threads, so there is little
   * overhead here.
   */
  for (t = all_threads; t; t=t->next)
    {
      DEBUG((stderr, "wego_schedule: thread %d is in state run %d available %d waiting %d\n",t->thread_id, t->running, t->available, t->waiting))
	if (t->running)
	  running++;
      if (t->available)
	num_available++;
    }
  
  /*
   * Make sure we have enough threads to keep all the CPUs busy.  Sometimes
   * this creates more threads than we need, but that is ok.  They will just
   * sit idle.
   */
  while (num_available + running <= max_threads)
    {
      DEBUG((stderr, "wego_schedule: starting extra thread, available %d running %d < max %d\n", num_available , running , max_threads) )
	if (wego_new_thread())
	  num_available++;
	else
	  {
	    max_threads = num_available + running - 1 ;
	  }
    }
  
  DEBUG((stderr, "wego_schedule: has running %d available %d\n",running,num_available)) ;
  
  /*
   * First priority is to threads that are waiting for results.  These are
   * in the middle of their execution.
   */
  for (t=all_threads; t; t=t->next)
    {
      if (running >= max_threads)
	break;
      if (t->waiting)
	{
	  DEBUG((stderr, "wego_schedule: waking %d\n",t->thread_id))
	    t->waiting=0;
	  t->running=1;
	  /* should we try to wake anybody here? */
	  running++;
	}
    }
  
  /*
   * If we still do not have all the CPUs busy, look for new tasks.
   * If there are any, wake up threads to start running them.
   */
  if (queued_tasks)
    {
      DEBUG((stderr, "wego_schedule: look for thread starters\n"))
	for (t=all_threads; t; t=t->next)
	  {
	    if (t->available)
	      {
		DEBUG((stderr, "wego_schedule: waking available %d\n",t->thread_id))
		  t->available = 0 ;
		t->running = 1 ;
		wego_wakeup_threadrunner (t);
		running++ ;
	      }
	  }
    }
  
  DEBUG((stderr, "wego_schedule in %d done\n\n\n",TID))
    /*
     * this broadcast will wake up everybody who is sleeping
     * This is too many, but only the ones 
     * in  wego_want_to_run will really start
     */
    pthread_cond_broadcast(&thread_take_condition);
  pthread_mutex_unlock(&thread_take);
} /* wego_schedule */

/************************************************************************/
/************************************************************************/
#include "../wh/ac.h"  /* declares arrayReport */

static void wego_do_max_threads (int max)
{
  int n;

  if (max_threads)
    {
      if (max > max_threads)
	max_threads = max ;
      return ;
    }
  max_threads = max ;
  memSetIsMultiThreaded () ; /* set memsubs in thread-safe mode */

  arrayReport (-2) ;   /* blocks static array counts */
  bigArrayReport (-2) ;   /* blocks static bigArray counts */
  n = pthread_key_create (&threadrunner_key, NULL );
  if (n < 0)
    perror ("pthread_key_create");
  pthread_setspecific (threadrunner_key, &main_thread);
  aceInWaitAndRetry (1) ;
  wego_log_channel = channelCreate (256, char *, 0) ;
  wego_go (wego_export_logs, &n, int) ;
  return ;
} /* wego_do_max_threads */

/**************************/

void wego_max_threads (int max)
{
  return wego_do_max_threads (max) ;
} /* wego_max_threads */

/**************************/

void wego_init (void)
{
  return wego_do_max_threads (3) ;
} /* wego_init */

/************************************************************************/
/************************************************************************/

/* 
* You wish these two functions were implemented.  They will come later.
*
* Leaking this memory is undesirable, but not so bad for initial testing.
*/
void wego_task_destroy (WEGO_TASK *task)
{
}

void wego_group_destroy (struct wego_group *group )
{
}

/************************************************************************/
/*
 * This function is not called from anywhere.  You can use it in gdb
 * with "print wego_dump_threads()"
 */

int wego_dump_threads ()
{
  struct threadrunner *t;
  for (t = all_threads; t; t=t->next)
    {
      printf("%d\t",t->thread_id);
      printf("run %d avail %d wait %d main %d\n",t->running, t->available, t->waiting, t->main);
    }
  return 0;
} /* wego_dump_threads */

/*************************************************************************/
/*************************************************************************/

