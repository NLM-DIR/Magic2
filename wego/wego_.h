#ifndef WEGO__H
#define WEGO__H

/*
************************************************************************
*** wego - an easy library for simple SMP/multithreading in C
*** header to the full wego library as developped by mark Sienkiewicz
*** some headers are already defined in wego.h
************************************************************************
* $Header: /home/mieg/aa/CVS_ACEVIEW/ace/wh/wego_.h,v 1.1 2017/01/23 21:46:58 mieg Exp $
*/
#include <wego.h>

/*
* A "task" is a description of work that will be done in a thread.
* When you create a task, you specify a function that will perform
* the task, and you give parameters the function will use.
*
* Some time later, you can ask for the result from the function.
*/
typedef struct wego_task WEGO_TASK ;

/*
* A "group" is a way to collect several tasks into a group. You
* can ask for any task that is finished.  (That is, you can collect
* the results in the order that the tasks finish -- if you did
* not have a group, you would have to ask for results from
* specific tasks in a specific order.)
*/
typedef struct wego_group WEGO_GROUP ;

/*
* wego_new_group() creates a group with no tasks in it.
*/
extern WEGO_GROUP * wego_new_group();

/*
* wego_run() creates a new task. Parameters are:
*	group 	
*		is the wego_group that the task should be added to.
*		Use NULL if the task is not part of a group.
*	function 
*		is the function that will be called to perform the computation.
*	parameter_struct
*	parameter_struct_size
*		is a block of memory passed to the function as parameters.
*		It is copied by wego_run(); you can free or reuse the
*		memory immediately.
*	result_struct
*	result_struct_size
*		is a block of memory that the function may store returned
*		results in.  result_struct may be NULL; the library will
*		allocate memory for it.
*
* The return value is a descriptor of the task created.
*
* The task will be queued until a thread is available to run it.
*
*/
extern WEGO_TASK *wego_run (WEGO_GROUP *group, 
		void (*function)(WEGO_TASK *t, void *p, void *r),
		void *parameter_struct, size_t parameter_struct_size, 
		void *result_struct, size_t result_struct_size
	);

/* we copy the structure to respect the paradigm: aliased or mutable */
WEGO_TASK *wego_go_size (void(*f)(const void *vp), const void *vp, int size, WEGO_GROUP *group) ;

#define wego_group_go(_group,_f,_vp,type) wego_go_size((_f),(_vp),sizeof(type),(_group))

/*
* wego_result() returns the result pointer from the specified task.
* It blocks if the task is not finished.
*/
extern void * wego_result(WEGO_TASK *t);

/*
* wego_fin() searches the group for a wego_task that 
*	- is finished
*	- has not been returned by wego_fin() before
*	- has not been examined by wego_result()
* The idea is that you call wego_fin() several times to collect all
* the results.
*/
extern WEGO_TASK *wego_fin (WEGO_GROUP *group);

/*
* wego_task_destoy() destroys a wego_task.
*	
* wego_free_group() destroys a wego_group and all the wego_tasks
* that it contains.  
*
* You may free a task or group before it runs or while it is running.
* You may or may not save some CPU time by avoiding unneccesary work.
* The library does not abort running code, and it does not absolutely
* guarantee not to start a task that you have removed.
* 
* wego_task_destoy() and wego_group_destoy() are not implemented yet.
* The memory is leaked.
*/
extern void wego_task_destroy (WEGO_TASK *task);
extern void wego_group_destroy (WEGO_GROUP *task);

void wego_exit (int err) ;

/*
************************************************************************
*** end of file
************************************************************************
*/

#endif
