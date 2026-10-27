/*
************************************************************************
*** wego - an easy library for simple SMP/multithreading in C
*** Author Mark Sienkiewicz, Jean Thierry-Mieg, NCBI
************************************************************************
* $Header: /home/mieg/aa/CVS_ACEVIEW/ace/wego/wego.c,v 1.37 2020/06/29 01:20:13 mieg Exp $
*
* Patched version:
*  - C23 clean: the task function is a tagged union of the two real
*    signatures (full / light), no more "void (*)()" prototypes.
*  - asserts check the tag before every call.
*  - FIFO task queue, per-group finished list (wego_fin is O(1) per call),
*    no lost wakeup in wego_fin, no thundering herd, wego_exit no longer
*    kills the pool thread, light tasks are freed, wego_flush no longer
*    uses sleep(1), thread-safe first initialisation.
*/

#ifdef SUN
#include <malloc.h>
#endif
#include <assert.h>
#include <errno.h>
#include <pthread.h>
#include <setjmp.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "wego_.h"
#include "ac.h"
#include "../wh/ac.h"  /* declares arrayReport */
#include "channel.h"

#define wego_debug 0

#if wego_debug
#define DEBUG(_x) do { fprintf _x ; } while (0)
#else
#define DEBUG(_x) do { } while (0)
#endif

/************************************************************************/

/*
* Which member of the function union is valid.  Set by the code that
* creates the task, checked (assert) by the thread that calls it.
* 0 means "never set" and is a bug.
*/
enum { WEGO_KIND_UNSET = 0, WEGO_KIND_FULL = 1, WEGO_KIND_LIGHT = 2 } ;

#define WEGO_TASK_ALLOC_RESULT 	1
#define WEGO_TASK_DETACHED 	2	/* nobody holds the task pointer: free it after the run */

/*
* information about a single task
*/
struct wego_task
{
  union
  {
    void (*full)(struct wego_task *t, void *p, void *r) ;   /* wego_run     */
    void (*light)(const void *vp) ;                         /* wego_go_size */
  } function ;
  int	kind ;			/* tag of the union: WEGO_KIND_xxx */

  void *	parameter_struct;	/* parameters from user */
  size_t	parameter_struct_size;
  void *	result_struct;		/* place to store return values */
  size_t	result_struct_size;

  int	flags;

  struct wego_group *group;	/* what group are we in? */

  struct wego_task *q_next;	/* next in queue of all pending tasks */
  struct wego_task *g_next;	/* next in list of all tasks in group */
  struct wego_task *done_next;	/* next in the group list of finished tasks */

  /*
   * fields for applications waiting for our result
   */
  pthread_mutex_t		finished_semaphore;
  pthread_cond_t		finished_condition;
  int			finished_flag;
  int			reportable_flag;
  int			ever_reported;
};

/************************************************************************/

/*
* group of tasks
*/
struct wego_group
{
  pthread_mutex_t		group_semaphore;
  pthread_cond_t		group_condition;
  struct wego_task *task_list;	/* all tasks of the group */
  int	task_count;
  struct wego_task *done_head;	/* finished tasks, in order of completion */
  struct wego_task *done_tail;
  int	unreported;		/* tasks not yet reported by wego_fin or wego_result */
};

/************************************************************************/

/*
* Thread management variables
* You must have thread_take semaphore to access any of these variables:
*/

	static pthread_mutex_t	thread_take = PTHREAD_MUTEX_INITIALIZER;

	struct threadrunner
		{
		int			main;		/* main program */
		int			running;	/* using a cpu */
		int			available;	/* available to run a task */
		int			waiting;	/* waiting to continue */
		int			thread_id;	/* small integer thread number */
		int			started;	/* set by the new thread at startup */
		int			in_task;	/* currently inside a user task */
		jmp_buf			exit_env;	/* target of wego_exit */
		pthread_cond_t		condition;	/* with semaphore: wake a idle runner */
		pthread_mutex_t		semaphore;
		pthread_cond_t		run_condition;	/* with thread_take: wait for a cpu slot */
		struct threadrunner *	next;
		};

	static struct threadrunner 	main_thread =
	  {
	    .main = 1,
	    .running = 1,
	    .condition = PTHREAD_COND_INITIALIZER,
	    .semaphore = PTHREAD_MUTEX_INITIALIZER,
	    .run_condition = PTHREAD_COND_INITIALIZER
	  } ;
	static struct threadrunner *	all_threads = &main_thread;

	static pthread_cond_t		thread_startup_condition = PTHREAD_COND_INITIALIZER;
	static pthread_mutex_t		thread_startup_lock = PTHREAD_MUTEX_INITIALIZER;

	/*
	* FIFO queue of all tasks waiting to be run
	*/
	static struct wego_task *queued_tasks = NULL;
	static struct wego_task *queued_tail = NULL;
	static int queued_count = 0;	/* tasks in the queue */
	static int starters = 0;	/* threads woken to take a task, not yet arrived */

/************************************************************************/

/*
* forward declarations
*/
static void wego_want_to_sleep (void);
static void wego_want_to_run (void);
static void wego_schedule (void);
static void *wego_threadrunner (void *thread_parameter);

/************************************************************************/

pthread_key_t threadrunner_key;

#define TID ( ((struct threadrunner *)pthread_getspecific(threadrunner_key))->thread_id )

/************************************************************************/

/*
* This is the maximum number of threads that can be running concurrently.
* You can set this from your main program before you call any wego
* functions.  If you do not, the default is 8.
*
* NOTE: new tasks are started as soon as an idle thread exists, even if
* that goes above max_threads, because a running task may well be blocked
* on a channel and we cannot know it.  Enforcing a hard cap could deadlock
* channel pipelines.
*
* If you change this value after starting some threads, the results are
* undefined.
*/
static int max_threads = 0 ;

/*
* first initialisation, protected against concurrent first calls
*/
static pthread_mutex_t wego_init_lock = PTHREAD_MUTEX_INITIALIZER ;
static int wego_is_initialized = 0 ;

static void wego_ensure_init (void)
{
  if (! __atomic_load_n (&wego_is_initialized, __ATOMIC_ACQUIRE))
    wego_max_threads (8) ;   /* a security if wego has not been initialized */
} /* wego_ensure_init */

/************************************************************************/
/************************************************************************/
/*
 * task allocation helpers
 */

static struct wego_task *wego_task_new (void)
{
  struct wego_task *r = calloc (1, sizeof (*r)) ;
  if (r)
    {
      pthread_mutex_init (&r->finished_semaphore, NULL) ;
      pthread_cond_init (&r->finished_condition, NULL) ;
    }
  return r ;
} /* wego_task_new */

static void wego_task_free (struct wego_task *r)
{
  free (r->parameter_struct) ;
  if (r->flags & WEGO_TASK_ALLOC_RESULT)
    free (r->result_struct) ;
  pthread_mutex_destroy (&r->finished_semaphore) ;
  pthread_cond_destroy (&r->finished_condition) ;
  free (r) ;
} /* wego_task_free */

/*
 * register the task in its group, then put it at the tail of the queue.
 * The group is done first so that a task can never finish before the
 * group knows about it.
 */
static void wego_enqueue (struct wego_task *r, struct wego_group *group)
{
  r->group = group ;

  if (group)
    {
      pthread_mutex_lock (&group->group_semaphore) ;
      group->task_count++ ;
      group->unreported++ ;
      r->g_next = group->task_list ;
      group->task_list = r ;
      pthread_mutex_unlock (&group->group_semaphore) ;
    }

  pthread_mutex_lock (&thread_take) ;
  r->q_next = NULL ;
  if (queued_tail)
    queued_tail->q_next = r ;
  else
    queued_tasks = r ;
  queued_tail = r ;
  queued_count++ ;
  pthread_mutex_unlock (&thread_take) ;
} /* wego_enqueue */

/************************************************************************/
/************************************************************************/
/*
 * Create a new task.
 *
 ** There are no concurrency issues until we enter the task into the
 ** common data structures.
 */

struct wego_task *wego_run (struct wego_group *group,
			    void (*function)(struct wego_task *, void *, void *),
			    void *parameter_struct, size_t parameter_struct_size,
			    void *result_struct, size_t result_struct_size
			    )
{
  struct wego_task *r;

  wego_ensure_init () ;
  DEBUG((stderr, "wego_run called\n"));

  r = wego_task_new () ;
  if (!r)
    return NULL;

  r->function.full = function ;
  r->kind = WEGO_KIND_FULL ;

  if (parameter_struct && parameter_struct_size)
    {
      r->parameter_struct = malloc (parameter_struct_size) ;
      if (! r->parameter_struct)
	{
	  wego_task_free (r) ;
	  return NULL ;
	}
      memcpy (r->parameter_struct, parameter_struct, parameter_struct_size) ;
    }
  r->parameter_struct_size = parameter_struct_size ;

  if (result_struct)
    r->result_struct = result_struct ;
  else if (result_struct_size)
    {
      r->result_struct = calloc (1, result_struct_size) ;
      if (! r->result_struct)
	{
	  wego_task_free (r) ;
	  return NULL ;
	}
      r->flags |= WEGO_TASK_ALLOC_RESULT ;
    }
  r->result_struct_size = result_struct_size ;

  wego_enqueue (r, group) ;
  wego_schedule () ;
  return r ;
} /* wego_run */

/************************************************************************/
/* Light interface expecting all communications to be handled via channels
* Create a new task.
*
** There are no concurrency issues until we enter the task into the
** common data structures.
*
* detached: nobody will ever see the task pointer (wego_go_), so the
* thread that runs the task frees everything when it is done.
*/

static WEGO_TASK *wego_go_internal (void (*f)(const void *vp), const void *vp, size_t size,
				    struct wego_group *group, int detached)
{
  struct wego_task *r ;

  assert (! (detached && group)) ;
  DEBUG((stderr, "wego_go_internal: start\n"));

  r = wego_task_new () ;
  if (!r)
    return NULL ;

  r->function.light = f ;
  r->kind = WEGO_KIND_LIGHT ;
  if (detached)
    r->flags |= WEGO_TASK_DETACHED ;

  if (vp && size)
    {
      r->parameter_struct = malloc (size) ;
      if (! r->parameter_struct)
	{
	  wego_task_free (r) ;
	  return NULL ;
	}
      memcpy (r->parameter_struct, vp, size) ;
    }
  r->parameter_struct_size = size ;

  wego_enqueue (r, group) ;
  wego_schedule () ;
  /* if detached, r may already be freed: do not touch it */
  return r ;
} /* wego_go_internal */

WEGO_TASK *wego_go_size (void(*f)(const void *vp), const void *vp, size_t size, struct wego_group *group)
{
  wego_ensure_init () ;
  return wego_go_internal (f, vp, size, group, 0) ;
} /* wego_go_size */

void wego_go_ (void(*f)(const void *vp), const void *vp, size_t size)
{
  wego_ensure_init () ;
  wego_go_internal (f, vp, size, 0, 1) ;
} /* wego_go_ */

/*
 * Terminate the current task early.  Inside a pool thread we jump back to
 * the runner loop, so the task is marked finished and the thread is
 * recycled.  (pthread_exit would kill the pool thread while it still
 * counts as running, and the task would never be reported.)
 */
void wego_exit (int err)
{
  struct threadrunner *t = NULL ;

  (void) err ;
  if (__atomic_load_n (&wego_is_initialized, __ATOMIC_ACQUIRE))
    t = pthread_getspecific (threadrunner_key) ;
  if (t && t->in_task)
    longjmp (t->exit_env, 1) ;
  pthread_exit (0) ;
} /* wego_exit */

/*************************************************************************************/
/*************************************************************************************/
static CHAN *wego_log_channel = 0 ;
static pthread_mutex_t wego_log_lock = PTHREAD_MUTEX_INITIALIZER ;
static pthread_cond_t wego_log_cond = PTHREAD_COND_INITIALIZER ;
static int wego_log_done = 0 ;

void wego_log (const char *text)
{
  char *cp ;

  if (! wego_log_channel)   /* wego not initialised: no serialisation possible */
    {
      fprintf (stderr, "// %s: %s\n", timeShowNow (), text) ;
      return ;
    }
  cp = strnew (text, 0) ;
  channelPut (wego_log_channel, &cp, const char**) ;
} /* wego_log  */

/************************************************************************/

static void wego_export_logs (const void *vp)
{
  const char *text = 0 ;

  (void) vp ;
  while (channelMultiGet (wego_log_channel, &text, 1, const char*))
    if (text)
      {
	fprintf (stderr, "// %s: %s\n", timeShowNow (), text) ;
	ac_free (text) ;
      }

  /* tell wego_flush that everything was written */
  pthread_mutex_lock (&wego_log_lock) ;
  wego_log_done = 1 ;
  pthread_cond_broadcast (&wego_log_cond) ;
  pthread_mutex_unlock (&wego_log_lock) ;
} /* wego_export_logs */

/************************************************************************/

void wego_flush (void)
{
  struct timespec ts ;

  if (! wego_log_channel)
    return ;
  channelClose (wego_log_channel) ;

  /* wait until wego_export_logs has drained the channel, but never hang
   * for more than 5 seconds if something unexpected happened */
  clock_gettime (CLOCK_REALTIME, &ts) ;
  ts.tv_sec += 5 ;
  pthread_mutex_lock (&wego_log_lock) ;
  while (! wego_log_done)
    if (pthread_cond_timedwait (&wego_log_cond, &wego_log_lock, &ts) == ETIMEDOUT)
      break ;
  pthread_mutex_unlock (&wego_log_lock) ;
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
  g = calloc (1, sizeof(*g));
  if (!g)
    return NULL;
  pthread_mutex_init (&g->group_semaphore, NULL);
  pthread_cond_init (&g->group_condition, NULL);
  DEBUG((stderr, "wego_new_group return %p\n", (void *) g)) ;
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
* Lock order: group_semaphore -> finished_semaphore (only in wego_fin),
* and group_semaphore -> thread_take (wego_want_to_sleep inside wego_fin).
* Nothing ever takes them in the reverse order.
*/

void *wego_result (struct wego_task *t)
{
  void *r;
  int first ;

  pthread_mutex_lock(&t->finished_semaphore);
  if (! t->finished_flag)
    {
      /* the lock is held from the test to the cond_wait: no lost wakeup */
      wego_want_to_sleep();
      while (! t->finished_flag)
	pthread_cond_wait( &t->finished_condition, &t->finished_semaphore );
      pthread_mutex_unlock(&t->finished_semaphore);
      wego_want_to_run();     /* never wait for a cpu slot holding a task lock */
      pthread_mutex_lock(&t->finished_semaphore);
    }
  t->reportable_flag = 0;
  first = ! t->ever_reported ;
  t->ever_reported = 1;
  r = t->result_struct;
  pthread_mutex_unlock(&t->finished_semaphore);

  /* the group must know that this task does not need to be reported any more */
  if (first && t->group)
    {
      pthread_mutex_lock (&t->group->group_semaphore) ;
      t->group->unreported-- ;
      pthread_cond_broadcast (&t->group->group_condition) ;
      pthread_mutex_unlock (&t->group->group_semaphore) ;
    }
  return r;
} /* wego_result */

/************************************************************************/

struct wego_task *wego_fin (struct wego_group *group)
{
  struct wego_task *t;
  int ok ;

  /*
   * The runner threads append each finished task to group->done list,
   * so we only look at tasks that are really finished.
   * The group semaphore is held from the examination of the list to the
   * cond_wait, so a task finishing in between cannot be missed.
   */
  pthread_mutex_lock (&group->group_semaphore) ;
  for (;;)
    {
      while ((t = group->done_head) != NULL)
	{
	  group->done_head = t->done_next ;
	  if (! group->done_head)
	    group->done_tail = NULL ;
	  t->done_next = NULL ;

	  pthread_mutex_lock (&t->finished_semaphore) ;
	  ok = t->reportable_flag ;     /* 0 if wego_result already consumed it */
	  if (ok)
	    {
	      t->reportable_flag = 0 ;
	      t->ever_reported = 1 ;
	    }
	  pthread_mutex_unlock (&t->finished_semaphore) ;

	  if (ok)
	    {
	      group->unreported-- ;
	      pthread_mutex_unlock (&group->group_semaphore) ;
	      DEBUG((stderr, "wego_fin return data\n")) ;
	      return t ;
	    }
	}

      /*
       * If all tasks of the group have been reported, return NULL.
       */
      if (group->unreported == 0)
	{
	  pthread_mutex_unlock (&group->group_semaphore) ;
	  DEBUG((stderr, "wego_fin return NULL\n")) ;
	  return NULL ;
	}

      /*
       * Nothing finished yet: give up our cpu slot and sleep.
       */
      wego_want_to_sleep () ;
      while (! group->done_head && group->unreported > 0)
	pthread_cond_wait (&group->group_condition, &group->group_semaphore) ;
      pthread_mutex_unlock (&group->group_semaphore) ;
      wego_want_to_run () ;
      pthread_mutex_lock (&group->group_semaphore) ;
    }
} /* wego_fin */

/************************************************************************/

static void wego_want_to_sleep (void)
{
  struct threadrunner *t;
  pthread_mutex_lock(&thread_take);
  t = pthread_getspecific(threadrunner_key);
  t->running = 0;
  t->waiting = 1;
  pthread_mutex_unlock(&thread_take);
  wego_schedule();
} /* wego_want_to_sleep */

/************************************************************************/

static void wego_want_to_run (void)
{
  struct threadrunner *t;

  t = pthread_getspecific(threadrunner_key);

  pthread_mutex_lock(&thread_take);
  while (! t->running)
    pthread_cond_wait(&t->run_condition, &thread_take);   /* signalled by wego_schedule only for us */
  t->waiting = 0;
  pthread_mutex_unlock(&thread_take);
} /* wego_want_to_run */

/************************************************************************/
/*
 * This is the main program of each thread.
 */

static void *wego_threadrunner (void *thread_parameter)
{
  struct threadrunner *t = thread_parameter ;

  pthread_mutex_lock(&thread_startup_lock);
  pthread_setspecific(threadrunner_key, t);

  /* wego_new_thread already created us with available = 1 */
  t->started = 1;

  pthread_cond_signal(&thread_startup_condition);
  pthread_mutex_unlock(&thread_startup_lock);

  /*
   * Wait for: we can run a task
   */
  for (int first = 1 ; ; first = 0)
    {
      struct wego_task *w;
      struct wego_group *g;

      /*
       * set our thread status to show that:
       * - we are not using CPU
       * - we are not waiting to run
       * - we are available to run a task
       *
       * Not on the very first pass: the scheduler may already have
       * woken this brand new thread (it is in all_threads and available
       * as soon as wego_new_thread returns), and resetting the state
       * here would erase that wakeup and strand the queued task.
       */
      if (! first)
	{
	  pthread_mutex_lock(&thread_take);
	  t->running = 0;
	  t->available = 1;
	  t->waiting = 0;
	  pthread_mutex_unlock(&thread_take);
	}

      /*
       * make sure the scheduler thinks about our new state
       */
      wego_schedule();

      /*
       * sleep until the scheduler marks us as running.  the scheduler
       * sets running to 1, available to 0.
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
      if (starters > 0)
	starters--;
      w = queued_tasks;
      if (!w)
	{
	  pthread_mutex_unlock(&thread_take);
	  continue;
	}
      queued_tasks = w->q_next;
      if (! queued_tasks)
	queued_tail = NULL;
      queued_count--;
      pthread_mutex_unlock(&thread_take);

      /*
       * execute the user's function, after checking the tag of the union.
       * wego_exit() jumps back here.
       */
      t->in_task = 1;
      if (setjmp (t->exit_env) == 0)
	{
	  assert (w->kind == WEGO_KIND_FULL || w->kind == WEGO_KIND_LIGHT) ;
	  if (w->kind == WEGO_KIND_LIGHT)
	    {
	      assert (w->function.light != NULL) ;
	      w->function.light (w->parameter_struct) ;
	    }
	  else
	    {
	      assert (w->function.full != NULL) ;
	      w->function.full (w, w->parameter_struct, w->result_struct) ;
	    }
	}
      t->in_task = 0;

      /*
       * Light tasks never look at their parameters again.
       * Detached tasks (wego_go_) have nobody waiting: free everything.
       */
      if (w->kind == WEGO_KIND_LIGHT)
	{
	  free (w->parameter_struct) ;
	  w->parameter_struct = NULL ;
	}
      if (w->flags & WEGO_TASK_DETACHED)
	{
	  wego_task_free (w) ;
	  continue ;
	}

      /*
       * Flag the task as finished.
       * Wake up anybody who is waiting for the results.
       */
      pthread_mutex_lock(&w->finished_semaphore);
      w->finished_flag = 1;
      w->reportable_flag = 1;
      pthread_cond_broadcast(&w->finished_condition);
      pthread_mutex_unlock(&w->finished_semaphore);

      /*
       * If the task was part of a group, put it on the list of finished
       * tasks and wake up anybody who is waiting for group results.
       */
      g=w->group;
      if (g)
	{
	  pthread_mutex_lock(&g->group_semaphore);
	  w->done_next = NULL ;
	  if (g->done_tail)
	    g->done_tail->done_next = w ;
	  else
	    g->done_head = w ;
	  g->done_tail = w ;
	  pthread_cond_broadcast(&g->group_condition);
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
/* called with thread_take held, from wego_schedule */

static int wego_new_thread (void)
{
  pthread_t thread;
  struct threadrunner *t;
  static int thread_id = 1;
  int error = 0 ;
  int ok = 0 ;

  t = calloc (1, sizeof(*t));
  if (! t)
    return 0 ;
  t->available = 1;
  t->thread_id = thread_id++;
  t->next = all_threads;
  pthread_cond_init(&t->condition,NULL);
  pthread_mutex_init(&t->semaphore,NULL);
  pthread_cond_init(&t->run_condition,NULL);

  pthread_mutex_lock(&thread_startup_lock);
  error = pthread_create(&thread, NULL, wego_threadrunner, (void *)t) ;
  if (! error)
    {
      pthread_detach (thread) ;
      while (! t->started)    /* predicate loop: spurious wakeups are harmless */
	pthread_cond_wait(&thread_startup_condition, &thread_startup_lock);
      all_threads = t;
      ok = 1 ;
    }
  else /* creation failed */
    {
      fprintf (stderr, "ERROR  wego_new_thread: pthread_create returned error %d\n", error) ;
      pthread_cond_destroy (&t->condition) ;
      pthread_mutex_destroy (&t->semaphore) ;
      pthread_cond_destroy (&t->run_condition) ;
      free (t) ;
    }
  pthread_mutex_unlock(&thread_startup_lock);
  return ok ;
} /* wego_new_thread */

/************************************************************************/
/*
 * this function is used by wego_schedule() to wake up a specific
 * thread that is waiting in threadrunner().  The running flag is
 * written under the thread semaphore, the same lock the runner uses
 * to read it.  Called with thread_take held.
 */

void wego_wakeup_threadrunner (struct threadrunner *t)
{
  pthread_mutex_lock(&t->semaphore);
  t->running = 1 ;
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
  DEBUG((stderr, "wego_schedule, max_threads = %d\n", max_threads)) ;
  /*
   * See how many threads are running, how many are available.  I consider
   * counting them here to be more reliable than ++/-- all over the place.
   * Typically, I expect a small number of threads, so there is little
   * overhead here.
   */
  for (t = all_threads; t; t=t->next)
    {
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
      if (wego_new_thread())
	num_available++;
      else
	{
	  max_threads = num_available + running - 1 ;
	  break ;
	}
    }

  /*
   * First priority is to threads that are waiting for results.  These are
   * in the middle of their execution.  Each one is signalled individually.
   */
  for (t=all_threads; t; t=t->next)
    {
      if (running >= max_threads)
	break;
      if (t->waiting)
	{
	  t->waiting=0;
	  t->running=1;
	  running++;
	  pthread_cond_signal (&t->run_condition) ;
	}
    }

  /*
   * Look for new tasks.  Wake up one idle thread per queued task which
   * is not already claimed by a thread that has been woken but has not
   * yet reached the queue.  (See the comment at max_threads for why this
   * does not stop at max_threads.)
   */
  for (t=all_threads; t && queued_count > starters; t=t->next)
    {
      if (t->available)
	{
	  t->available = 0 ;
	  starters++ ;
	  running++ ;
	  wego_wakeup_threadrunner (t);
	}
    }

  pthread_mutex_unlock(&thread_take);
} /* wego_schedule */

/************************************************************************/
/************************************************************************/

static void wego_do_max_threads (int max)
{
  int n;

  pthread_mutex_lock (&wego_init_lock) ;
  if (__atomic_load_n (&wego_is_initialized, __ATOMIC_ACQUIRE))
    {
      pthread_mutex_lock (&thread_take) ;
      if (max > max_threads)
	max_threads = max ;
      pthread_mutex_unlock (&thread_take) ;
      pthread_mutex_unlock (&wego_init_lock) ;
      return ;
    }
  max_threads = max ;
  memSetIsMultiThreaded () ; /* set memsubs in thread-safe mode */

  arrayReport (-2) ;   /* blocks static array counts */
  bigArrayReport (-2) ;   /* blocks static bigArray counts */
  n = pthread_key_create (&threadrunner_key, NULL );
  if (n)   /* pthread_key_create returns 0 on success, an errno value on failure */
    fprintf (stderr, "ERROR wego: pthread_key_create returned %d\n", n) ;
  pthread_setspecific (threadrunner_key, &main_thread);
  aceInWaitAndRetry (1) ;
  wego_log_channel = channelCreate (256, char *, 0) ;
  wego_go_internal (wego_export_logs, &n, sizeof (int), 0, 1) ;
  __atomic_store_n (&wego_is_initialized, 1, __ATOMIC_RELEASE) ;
  pthread_mutex_unlock (&wego_init_lock) ;
  return ;
} /* wego_do_max_threads */

/**************************/

void wego_max_threads (int max)
{
  wego_do_max_threads (max) ;
} /* wego_max_threads */

/**************************/

void wego_init (void)
{
  wego_do_max_threads (3) ;
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
  (void) task ;
}

void wego_group_destroy (struct wego_group *group )
{
  (void) group ;
}

/************************************************************************/
/*
 * This function is not called from anywhere.  You can use it in gdb
 * with "print wego_dump_threads()"
 */

int wego_dump_threads (void)
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
