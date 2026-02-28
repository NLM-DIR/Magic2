/*  File: utils.c
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) J Thierry-Mieg and R Durbin, 1999
 * -------------------------------------------------------------------
  * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (Sanger Centre, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.crbm.cnrs-mop.fr
 *
 * Description: Contains utility functions that are "one-offs" and not
 *              really part of any other package.
 * Exported functions: getLogin
 * HISTORY:
 * Last edited: Sep 12 09:41 2001 (edgrif)
 * * Mar 22 14:42 1999 (edgrif): Added getSystemName().
 * * Jan 25 16:05 1999 (edgrif): Add getLogin from session.c
 * Created: Thu Jan 21 15:46:49 1999 (edgrif)
 * CVS info:   $Id: utils.c,v 1.53 2020/06/07 14:37:57 mieg Exp $
 *-------------------------------------------------------------------
 */

#include "ac.h"
#ifdef SGI
#include <math.h>					    /* work around SGI library bug */
#endif

#include <sys/resource.h>
#include <version.h>					    /* For UT_ macros. */

#include <ctype.h>                                          /* for toupper() */

/* moved here from freesubs.c, because freesubs.c is not thread safe
 * this is threadsafe provided nobody ever writes there
 * please access via ace_upper()/ace_lower
 */

const char UT_UPPER[128] =
{ 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
  16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
  32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,
  48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,
  64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,
  80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,
  96,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,
  80,81,82,83,84,85,86,87,88,89,90,123,124,125,126,127
} ;

const char UT_LOWER[128] =
{  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15, 
  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,  
  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  
  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  
  64,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122,  91,  92,  93,  94,  95,
  96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127
} ;

#ifdef IBM
/* I can't believe this is still required for IBM machines...sigh..          */

struct passwd {
  char    *pw_name;
  char    *pw_passwd;
  long   pw_uid;
  long   pw_gid;
  char    *pw_gecos;
  char    *pw_dir;
  char    *pw_shell;
} ;
extern struct passwd *getpwuid (long uid) ;
extern struct passwd *getpwnam (const char *name) ;

#else  /* !IBM */

#include <pwd.h>

#endif /* !IBM */


/* This is all pretty unsafe at the moment, it relies on the user having     */
/* called getLogin to initialise the below two globals before using them..   */
/* These globals are also separately set by session.c using the system       */
/* calls getuid() and geteuid(), this is not optimal.                        */
/*                                                                           */
uid_t ruid = -1, euid = -1;

#ifdef ED_G_NEVER_INCLUDE_THIS_CODE
/* Currently unused, there should be a unified interface to uids but this    */
/* requires more understanding of how they are used in the code before       */
/* implementing.                                                             */
uid_t utGetRuid(void)
{
  return ruid ;
}

void utSetRuid(uid_t newRuid)
{
  ruid = newRuid ;
}

uid_t utGetEuid(void)
{
  return euid ;
}

void utSetEuid(uid_t newEuid)
{
  euid = newEuid ;
}
#endif /* ED_G_NEVER_INCLUDE_THIS_CODE */



/*************************************************************/
/****************** getLogin() *******************************/

/**** getLogin() for real/effective user names 
***** original by 03.05.1992 kocab DKFZ  
*****/

char *getLogin (BOOL isReal)
     /* general public function */
{
  /* isReal for real or effective can't fail */

#if defined(MACINTOSH)

  return "acedb" ;
  
#elif defined (WIN32)

  char *name = NULL;
  if ((name = getenv("USERNAME")) != NULL )
    return name ;   /* Windows NT has usernames */
  else
    return "acedb" ;

#else  /* all UNIX */
  static char *rname = 0 ;
  static char *ename = 0 ;

/* RD 980417: changed so getlogin() was last resort.  It can return
   "root" inappropriately where getpwuid(ruid)->pw_name returns the
   right answer (ddts problem SANgc02099)
*/

  if (!ename)
    {
      if (ruid == -1)
	{ ruid = getuid() ;
	euid = geteuid () ;
	}
      
      if (!rname)
	if (ruid)
	  if (getpwuid(ruid))
	    rname = strnew(getpwuid(ruid)->pw_name, 0) ;
      if (!rname)
	if (getlogin())
	  rname = strnew(getlogin(), 0) ;
      if (!rname)
	rname = "getLoginFailed" ;
      
      if (!ename)
	if (getpwuid(euid))
	  ename = strnew(getpwuid(euid)->pw_name, 0) ;
      if (!ename)
	ename = rname ;
      /* mieg: jan5 99, new problem
	 on my dec alpha, if i run as a daemon, getLogin fails totally
	 printf ("rname=%s ename =%s",rname?rname:"null",ename ? ename:"null") ; 
	 */
    }

  return isReal ? rname : ename ;
#endif /* UNIX */

} /* getLogin */


/* Get the system name using the POSIX uname call.                           */
char *getSystemName(void)
{
  static char *sysname = NULL ;
  
  if (sysname == NULL)
    {
      struct utsname sys_details ;
      
      if (uname(&sys_details) != 0)
	{
	  /* POSIX failed... */
	  sysname = (char *) messalloc (101) ;
	  if (gethostname(sysname, 100) == -1)
	    messcrash("cannot retrieve system name.") ;
      }
      else
	sysname = strnew(sys_details.nodename, 0) ;
    }
  
  return sysname ;
}

char *getUserHomeDir (const char *user_name, AC_HANDLE h)
{
  struct passwd *pwd;
  char *homedir = NULL;

  if ((pwd = getpwnam (user_name)))
    homedir = strnew (pwd->pw_dir, h);
  else
    messerror ("Unknown user: %s", user_name);

  return homedir;
} /* getUserHomeDir */

/*****************************/

BOOL getCmdLineText (AC_HANDLE h, int *argcp, const char **argv,
		       const char *arg_name, const char **arg_val)
     /* call with (&argc, argv, "-option", &string)
      *   for options with value (e.g. -display <name>)
      * or with   (&argc, argv, "-option", NULL)
      *   for options that are simple switches (e.g. -help)
      *  --option and -option are treated as synonims 
      * RETURNS:
      *   TRUE if option was found (will set *arg_val if non-NULL)
      * SIDE-EFFECT:
      *   The option (and possibly its value) will be removed from
      *   the aguments (argc,argv). Call only once per option.
      */
{
  int i, num ;
  BOOL isFound = FALSE;

  for (i = 1; i < *argcp; i++)
    {
      if (strcmp(argv[i], arg_name) == 0 || (*argv[i] == '-' && strcmp(1 + argv[i], arg_name) == 0))
	{
	  if (arg_val)
	    {
	      if ((*argcp - i) < 2)
		return  FALSE ;

	      *arg_val = strnew(argv[i+1], h);
	      num = 2;
	    }
	  else
	    num = 1;

	  /* clear argument(s) from the list */
	  for (i += num; i < *argcp; i++)
	    argv[i - num] = argv[i];
	  argv[*argcp - num] = 0;
	  (*argcp) -= num;

	  isFound = TRUE;
	  break;
        }
    }
  return isFound;
} /* getCmdLineOption */

BOOL getCmdLineOption (int *argcp, const char **argv,
		       const char *arg_name, const char **arg_val)
{
  return getCmdLineText (0, argcp, argv, arg_name, arg_val) ;
}
/* same as above but returns a BOOL */
BOOL getCmdLineBool (int *argcp, const char **argv,
		     const char *arg_name)
{
  return getCmdLineOption (argcp, argv, arg_name, 0) ;
} /* getCmdLineBool */

BOOL getCmdLineInt (int *argcp, const char **argv,
		    const char *arg_name, int *val)
{
  const char *ccp ;
  if (getCmdLineOption (argcp, argv, arg_name, &ccp))
    {
      int x = 0 ;  
      if (sscanf (ccp, "%d", &x) == 1)
	{
	  *val = x ;
	  messfree (ccp) ;
	  return TRUE ;
	}
    }
  return FALSE ;
} /* getCmdLineInt */

BOOL getCmdLineLong (int *argcp, const char **argv,
		     const char *arg_name, long int *val)
{
  const char *ccp ;
  if (getCmdLineOption (argcp, argv, arg_name, &ccp))
    {
      long int x = 0 ;  
      if (sscanf (ccp, "%ld", &x) == 1)
	{
	  *val = x ;
	  messfree (ccp) ;
	  return TRUE ;
	}
    }
  return FALSE ;
} /* getCmdLineLong */

/* same as above but returns a float */
BOOL getCmdLineFloat (int *argcp, const char **argv,
		    const char *arg_name, float *val)
{
  const char *ccp ;
  if (getCmdLineOption (argcp, argv, arg_name, &ccp))
    {
      float x = 0 ;
      if (sscanf (ccp, "%f", &x) == 1)
	{
	  *val = x ;
	  messfree (ccp) ;
	  return TRUE ;
	}
    }
  return FALSE ;
} /* getCmdLineFloat */

/*****************************/

/* Check to see if user supplied the "-sleep secs" option, this is useful    */
/* for debugging when the acedb program is called from within a script. It   */
/* means that you can make the acedb program sleep until you can attach a    */
/* debugger to see what its doing.                                           */
/*                                                                           */
/* This routine is not perfect but its only intended for debugging use.      */
/*                                                                           */
void checkCmdLineForSleep(int *argcp, const char **argv)
{
  const char *secs_str = NULL ;

  if (getCmdLineOption(argcp, argv, "-sleep", &secs_str))
    {
      int secs = 0 ;

      if (secs_str != NULL)
	secs = atoi(secs_str) ;

      if (secs == 0)					    /* atoi returns 0 on error. */
	secs = 60 ;

      sleep(secs) ;
    }

  return ;
}

/*****************************/

#ifdef JUNK

/* This routine tries to set various process resource limits to their max    */
/* values. We especially need this for acedb which certainly requires lots   */
/* of memory and sometimes a large stack as well.                            */
/*                                                                           */
/* The routine reports on any resources that cannot be increased to a        */
/* reasonable level and gives the user the choice of whether to continue.    */
/*                                                                           */
/* The routine will be a noop for machines that don't support the required   */
/* get/setrlimit calls.                                                      */
/*                                                                           */
void utUnlimitResources(BOOL allow_user_abort)
{

#if defined(RLIMIT_DATA) && defined(RLIMIT_STACK)	    /* Usually in sys/resource.h */

  typedef struct ResourceLimit_ 
  {
    int resource ;					    /* resource id from sys/resource.h */
    char *resource_str ;				    /* Stringified resource id. */
    char *resname ;					    /* our string name for resource. */
    int min_limit ;					    /* our required minimum. */
  } ResourceLimit ;

  ResourceLimit our_limits[] = {{RLIMIT_DATA, UT_PUTSTRING(RLIMIT_DATA), "heap size", 52428800},
				{RLIMIT_STACK, UT_PUTSTRING(RLIMIT_STACK), "stack size", 4194304}} ;
  int num_resources = UtArraySize(our_limits) ;
  struct rlimit reslimits ;
  Stack errors = NULL ;
  int i ;
  int result = 0;


  /* Check the hard limits....                                               */
  for (i = 0 ; i < num_resources ; i++)
    {
      if (getrlimit(our_limits[i].resource, &reslimits) != 0)
	messcrash("Unable to retrieve current limits for %s (%s), reason was %s",
		  our_limits[i].resource_str, our_limits[i].resname, messSysErrorText()) ;
      else
	{
	  if (reslimits.rlim_max < our_limits[i].min_limit)
	    {
	      char *resource = NULL ;

	      if (errors == NULL)
		{
		  errors = stackCreate(20) ;
		  pushText(errors,
			   "Your environment has the following user process \"hard\" limits"
			   " which are sufficiently low that they may cause acedb to crash"
			   " (NOTE that you need root permission to raise them): ") ;
		}
	      else
		catText(errors, ", ") ;

	      resource = hprintf(0, "%s (%d bytes)", our_limits[i].resname, reslimits.rlim_max) ;
	      catText(errors, resource) ;
	      ac_free(resource) ;
	    }
	  else
	    {
	      if (reslimits.rlim_cur < reslimits.rlim_max)
		{
		  reslimits.rlim_cur = reslimits.rlim_max ;
		  if ((result = setrlimit (our_limits[i].resource, &reslimits)) != 0)
#ifdef __CYGWIN__
		    if (result != 22)   /* on Windows, ignore "can't increase stacksize" msg. */
#endif
		    messcrash("Unable to increase current limit for %s, reason was %s",
			      our_limits[i].resname, messSysErrorText()) ;
		}
	    }
	}
    }

  if (errors != NULL)
    {
      if (allow_user_abort)
	{
	  if (!messQuery("%s. Do you you want to continue ?", popText(errors)))
	    messExit("User initiated exit.") ;
	}
      else
	{
	  /* n.b. this may be a noop if the messdump context has not been    */
	  /* initialised yet.                                                */
	  messdump("%s", popText(errors)) ;
	}
    }

#endif /* defined(RLIMIT_DATA) && defined(RLIMIT_STACK) */

  return ;
}
#endif

/************************************************************/

/* put "break invokeDebugger" in your favourite debugger init file */
void invokeDebugger (void) 
{
  static BOOL reentrant = FALSE ;

  if (!reentrant)
    { reentrant = TRUE ;
      messalloccheck() ;
      reentrant = FALSE ;
    }
}

/* This SGI workaround make the comiler crash! I've disabled it for now,
   pending determination if it is still needed, or fixes a long-lost library
   bug. - srk */
#if 0

/* work around SGI library bug */
#ifdef  SGI 
double log10 (double x) { return log(x) / 2.3025851 ; }
#endif

#endif


/*************************************************************/

/* calling utilsInit will initialise all the statics makeing
   later calls to getLogin and getSystemName effectively threadsafe
   */

void utilsInit (void)
{
  static int n = 0 ;

  if (n)
    messcrash ("utilsInit is not thread safe, it should be called only once") ;
  n = 1 ;
  getSystemName () ;
  getLogin (TRUE) ; 
  oneByteInitialize (0) ;
} /* utilsInit */

/****************** little string routines ************************/

/* Translates newlines into escaped newlines, i.e. '\n' becomes '\''n'     */
char *utCleanNewlines(char *s, AC_HANDLE handle)
{ 
  int n = 0 ;
  char *cp, *cq ;
  char *copy = NULL ;

  for (cp = s ; *cp ; ++cp)
    if (*cp == '\n') ++n ;

  if (n == 0)
    copy = s ;
  else
    {
      copy = (char *) halloc (cp-s+n+1, handle) ;
      for (cp = s, cq = copy ; *cp ; ++cp)
	if (*cp == '\n') 
	  { *cq++ = '\\' ; *cq++ = 'n' ; }
	else
	  *cq++ = *cp ;
      *cq = 0 ;
    }

  return copy ;
}


/*************************************************************************************/
/* check if the counts on the 2 strands are compatible
 * if yes, return the type
 * else return -1 ;

 * chi^2 5/100=3.84  1/100=6.64   1/1000=10.83   1/10000=15.02
 */
/* incompatibility at risk 1/1000 */
int chi2 (int a1, int a2, int b1, int b2, float *c2p)
{
  double r1, r2, c1, c2, t, o, u, chi=0 ;
  int flag = 0 ;

  if (a1 >= 0 && a2 >= 0 && b1 >= 0 && b2 >= 0)
    {
      if (a1 + a2 >=6 && b1 + b2 >= 6)
	{
	  flag |= 0x1 ; /* both strands */
	}
      if (a1 + a2 + b1 + b2 >= 0)
	{
	  r1 = a1 + b1 ; r2 = a2 + b2 ;
	  c1 = a1 + a2 ; c2 = b1 + b2 ; t = c1 + c2 ;
	  if (!r1 || ! r2 || ! c1 || !c2) return flag ; /* they are clearly compatible */
	  o = r1 * c1 / t ; u = a1 - o ; if (u < 0) u = -u ; if (u > 0.5) u -= 0.5 ; chi += u*u/o ;
	  o = r1 * c2 / t ; u = b1 - o ; if (u < 0) u = -u ; if (u > 0.5) u -= 0.5 ; chi += u*u/o ;
	  o = r2 * c1 / t ; u = a2 - o ; if (u < 0) u = -u ; if (u > 0.5) u -= 0.5 ; chi += u*u/o ;
	  o = r2 * c2 / t ; u = b2 - o ; if (u < 0) u = -u ; if (u > 0.5) u -= 0.5 ; chi += u*u/o ;
	  
	  if (c2p) *c2p = chi ;
	  if (chi > 10.83)  flag = 0x8 ; /* risk  1/1000 */
	  else if (chi > 6.64)  flag = 0x4 ; /* risk 10/1000 == 1% */
	  else if (chi > 3.84)  flag = 0x2 ; /* risk 50/1000 == 5% */
	}
    }
  return flag ;
}

/****************** little arithmetic routines ************************/

int utArrondi (float x)
     /*   1.3 ->  1 ;  1.7 ->  2
      *  -1.3 -> -1 ; -1.7 -> -2 */
{
  if (x >= 0)
    return (int) (x + 0.5) ;
  else
    return (int) (x - 0.5) ;
}

/********************************************/
 /* returns 1, 2, 5, 10, 20, 50, 100 etc 
  * in absolute value the returned number is smaller than p
  *   19 -> 10  ; -19 -> -10  
  */
int utMainPart (float p)
{
  register int i = 1, sign = 1;
  
  if (!p)
    return 0 ;
  if (p < 0) 
    { sign  = -1 ; p = -p ; }
  if (p > 1000000000) p = 1000000000 ; /* avoid overflow */
  while(i <= p + .00001) 
    i = 10 * i;
  i /= 10 ;

  if (2 * i > p)
    return i * sign;
  if(5 * i > p)
    return 2 * i * sign;
  return
    5*i*sign;
} /* utMainPart */

/********************************************/
/* returns 1, 2, 5, 10, 20, 50, 100 etc 
 * but the returned number may be a bit bigger than p
 *   19 -> 20  ; -19 -> -20  
 */
int utMainRoundPart (float p)
{
  register int i = 1, sign = 1;

  if (!p)
    return 0 ;
  if (p < 0) 
    {sign = -1 ; p = -p ; } /* RD: ambiguous without spaces on SGI */
  if (p > 1000000000) p = 1000000000 ; /* avoid overflow */
  while (i < p + .1) 
    i = 10 * i;
  if (4 * i < 5 * p) 
    return i * sign;
  i /= 2;
  if (4 * i < 5 * p) 
    return i * sign;
  i =  (2 * i) / 5;
  if (4 * i < 5 * p)
    return i * sign;
  i = i/2; 
  return i * sign;
} /* utMainRoundPart */

/********************************************/
/* returns 1, 2, 5, 10, 20, 50, 100 etc 
 * but the returned number is guaranteed to >= p
 *   19 -> 20  ; -19 -> -20  
 */
int utUpperRoundPart (float p)
{
  register int i = 1, sign = 1;

  if (!p)
    return 0 ;
  if (p < 0) 
    {sign = -1 ; p = -p ; } /* RD: ambiguous without spaces on SGI */
  if (p > 1000000000) p = 1000000000 ; /* avoid overflow */
  if (p <= 1)
    return sign ;
  if (p <= 2)
    return 2*sign ;
  if (p <= 5)
    return 5*sign ;
  while (i < p + .1) 
    i = 10 * i;
  i /= 10 ;
   if (p <= i)
    return i * sign ;
  if (p <= 2 * i)
    return 2*i*sign ;
  if (p <= 5 * i)
    return 5*i*sign ;
  return 10*i*sign ;
} /* utUpperRoundPart */

/*********************************************************/
/* returns 1 , 2 , 5 ,10 , 20, 50 ,100 etc */
double utDoubleMainPart (double p)
{
  register double i = 1, sign = 1;

  if (!p)
    return 0.;
  if(p < 0) 
    {sign = -1; p = -p;}
  if (p > 1000000000) p = 1000000000 ; /* avoid overflow */
  i = (double) exp(log((double)10.0) * (double)(1 + (int)(log10((double)p)))) ;
  if (4 * i < 5 * p) 
    return i * sign;
  i /= 2;
  if (4 * i < 5 * p)
    return i * sign;
  i = (2 * i) / 5;
  if (4 * i < 5 * p)
    return i*sign;
  i = i/2;
  return i*sign;
}

/*************************************************************************************/

double utLogFac (int n)
{
  int i ;
  double z ;
  static double zz[1001] ;
  static BOOL firstPass = TRUE ;

  if (n <= 1) 
    return 1 ;
  if (firstPass)
    {
      firstPass = FALSE ;
      zz[0] = 0 ;
      for (i = 1, z = 1 ; i < 1000 ; i++)
	{
	  z += log (i) ;
	  zz[i] = z ;
	}
    }
  if (n < 1000)
    z = zz[n] ;
  else /* Stirling formula */
    {
      z = n ;  /* Stirling formula n! = sqrt(2 pi) n^(n+1/2) exp(-n) */
      z = (z + .5)* log(z) - z + log(2 * 3.141592)/2.0 ;
    }
  return z ;
}


static void  utSmoothHistoTestOne (int m, int n)
{
  Array hh = arrayCreate (101, double) ;
  double z = m/((double)n) ;
  int i, j, k, NN = 100 ;
  
  printf ("Test %d %d %.2f\n", m, n, z) ;

  utSmoothHisto (hh, m, n, 1) ;
  for (i = 0 ; i <= NN ; i++)
    {
      z = arr (hh, i, double) ;
      if (100 * z > 1)
	{
	  j = 100 * z ;
	  printf ("... %d\t%.1f\t", i, 100*z) ;
	  for (k = 0 ; k < j ; k++)
	    printf (".") ;
	  printf ("\n") ;
	}
    }
  arrayDestroy (hh) ;
} /* utSmoothHistoTestShow */
  
void utSmoothHistoTest (void)  /* call it from the debugger ! */
{
  int i, m, n ;
  
  for (n = 1 ; n < 100000 ; n *= 5)
    {
      for (i = 0 ; i <= 100 ; i += 20)
 	{
	  m = n * i /100.0 ;
	  utSmoothHistoTestOne (m, n) ;
	}
    }
} /* utSmoothHistoTest */

/*************************************************************************************/
/* accumulate in an Array of size 100 type double the smooth histo of surface 1 for m mutant among n cover 
 * for m, n large, delta (m/n) 
 * for m, n big, a Gaussian entered on m/n of width sqrt(n)
 * for m n small, Bayesian reciprocal binomial
 */
int utSmoothHisto (Array hh, int m, int n, int mult)
{
  int i, j, k, m1 ;
  static Array aaa = 0 ;
  double p, q, t, z, zmax ;
  int NN = 100 ;

  /* if n >= 100 and m <= 10 the Poisson is an excellent approximation
   * poisson(L) ~ gauss (mean=variance=L if L > 1000)
   * poisson(L) ~ gauss (mean=variance=L if L > 100) and P(X<x) replaced by P(X < x + 1/2)
   * if a variable x is Poisson(L), then sqrt(x) is Gauss(mean=sqtr(L), variance =1/4)
   * if number of arrival in time [0,t] follow Poisson(Lt), inter arrival time follow Exp(mean=1/L)
   */

  if (! hh)
    messcrash (" utSmoothHisto recieved a null Array") ;
  if (! arrayExists (hh))
    messcrash (" utSmoothHisto recieved a non exiting Array" ) ;
  if (hh->size != sizeof (double))
    messcrash (" utSmoothHisto should receive an Array of douubles") ;

  z = array (hh, NN, double) ;
  if (m > n) { i = m ; m = n ; n = i ; }
  if (m < 0 || n < 0)
    messcrash (" utSmoothHisto recieved negative integers cover = %d mutant = %d", n, m) ;

  if (n == 0) return 0 ;

  if (m == 0 || m == n)
    {
      double a[NN + 1] ;

      z = t = a[0] = 1 ;
      for (i = 1 ;  z > 1e-8 && i <= NN ; i++)
	{
	  p = i/((double)NN) ; q = 1 - p ;
	  for (z = 1, j = 0 ; z > 1e-6 &&  j < n ; j++)
	    z = q * z ;
	  t += z ;
	  a[i] = z ;
	}
     
      if (m == 0)
	{
	  for (j = 0 ; j < i ; j++)
	    arr (hh, j , double) += mult * a[j]/t ;
	}
      else /* m == n , symmetrize */
	{
	  for (j = 0 ; j < i ; j++)
	    arr (hh, NN - j , double) += mult * a[j]/t ;
	}
      return 1 ;
    }
  else if (n > 10000) /* keep a resolution of 1/NN */
    {
      z = NN *  m/((double)n) ; 
      i = z ; j = i + 1 ;
      if (j > NN) j = NN ;
      arr (hh, i, double) += mult * (1 - z + i) ;
      arr (hh, j, double) += mult * (z - i) ;

      return 1 ;
    }
  else if (n > 400) /* compute using Stirling formula */
    {
      double a[NN + 1] ;

      a[0] = a[NN] = t =  0 ;
      for (i = NN * m/((double)n), zmax = 0 ; i < NN ; i++)
	{
	  p = i/((double)NN) ; q = 1 - p ;
	  z = 0 ;
	  z += m * log (p) ;
	  z += (n-m) * log (q) ;
	    z += utLogFac(n) ;
	    z -= utLogFac(m) ;
	    z -= utLogFac(n-m) ;
	  a[i] = z ;
	  t += z ;
	  if (i == 1 || z > zmax) zmax = z ;
	  if (10000 * z < zmax) break ;
	}
      for (i =  NN * m/((double)n) - 1 ; i > 0 ; i--)
	{
	  p = i/((double)NN) ; q = 1 - p ;
	  z = 0 ;
	  z += m * log (p) ;
	  z += (n-m) * log (q) ;
	    z += utLogFac(n) ;
	    z -= utLogFac(m) ;
	    z -= utLogFac(n-m) ;
	    a[i] = exp (z) ;
	  t += z ;
	  if (i == 1 || z > zmax) zmax = z ;  
	  if (10000 * z < zmax) break ;
	}
      if (1000 * zmax < 1) /* fall back on large number estimation */
	{
	  z = NN * m/((double)n) ; 
	  i = z ; j = i + 1 ;
	  if (j > NN) j = NN ;
	  arr (hh, i, double) += mult * (1 - z + i) ;
	  arr (hh, j, double) +=  mult * (z - i) ;

	  return 1 ;
	}
      for (i = 1 ; i < NN ; i++)  	  
	arr (hh, i, double) += mult * a[i]/t ;
      return 1 ;
    }
  else /* precompute the distribs with small numbers in a lazy way */
    {
      double *a ;
      
      if (! aaa)
	{
	  aaa = arrayCreate (16000, double*) ;
	  a = array (aaa, 16000 - 1, double*) ; /* make room */
	}
      k = m + (n-2) * (n-2) ;
      a = array (aaa, k, double *) ;
      if (! a)  /* compute only once */
	{
	  a = array (aaa, k, double *) =  halloc ((NN+1)*sizeof (double), 0) ;
	  for (i = 0, t = 0 ; i <= NN ; i++)
	    {
	      /* Cmn = n!/m! (n-m)! p^m q^(n-m) */
	      p = i/((double)NN) ; q = 1 - p ;
	      if (m > n/2)
		{ z = p ; p = q ; q = z ; m1 = n - m ; }
	      else
		m1 = m ;
	      for (z = 1, j = 1 ; j <= m1 ;j++)
		z = z * p * (n + 1 - j)/j ; 
	      for (j = 1 ; j <= n - m1 ; j++)
		z = z * q ;
	      a[i] = z ; 
	      t += z ;
	    } 
	  for (i = 0 ; t> 0 && i <= NN ; i++)  
	    a[i] /= t ;
	}   
      for (i = 1 ; i < NN ; i++)  	  
	arr (hh, i, double) += mult * a[i] ;
      return 1 ;
    }
  
  return 1 ;
}  /* utSmoothHisto */


/*********************************************************/
#ifdef JUNK

/* Some of the unix number conversion routines are arcane to put it mildly,  */
/* here are some that (hopefully) work in a more straightforward way....     */

/* N.B. these should be integrated with the code in freesubs.c for freeint() */
/* etc. at some time.                                                        */
/* another version of these functionalities is implemented in acein.c        */

/* Given a string will attempt to convert it to an int, the whole string     */
/* must be a valid integer: [+-]?[0-9]*                                      */
/*                                                                           */
/* returns TRUE if the number could be converted and sets num_out to the     */
/* number, FALSE otherwise.                                                  */
       
BOOL utStr2LongInt (char *num_str, long int *num_out) ;      

BOOL utStr2Int (char *num_str, int *num_out)
{
  BOOL result = FALSE ;
  long int retval = 0 ;
  
  if ((result = utStr2LongInt(num_str, &retval)))
    {
      if (retval >= INT_MIN && retval <= INT_MAX)
	{
	  result = TRUE ;
	  *num_out = (int)retval ;
	}
    }

  return result ;
}


/* Given a string will attempt to convert it to a long int, the whole string */
/* must be a valid integer: [+-]?[0-9]*                                      */
/*                                                                           */
/* returns TRUE if the number could be converted and sets num_out to the     */
/* number, FALSE otherwise.                                                  */
/*                                                                           */
BOOL utStr2LongInt(char *num_str, long int *num_out)
{
  BOOL result = FALSE ;
  long int retval = 0 ;
  char *cp = NULL, *endptr = NULL ;
  int i, num_len, translation ;

  /* It is not possible to detect invalid numbers of the form "1232rubbish"  */
  /* because strtol() will convert as much as it can ("1232" in this case)   */
  /* and return as though everything is OK....sigh...                        */
  /* So here I check that the string is all digits...                        */
  num_len = strlen(num_str) ;
  cp = num_str ;
  if (*cp == '-' || *cp == '+')
    {
      cp++ ;
      num_len-- ;
    }
  for (i = 0, translation = 1 ; i < num_len && translation != 0 ; i++, cp++)
    {
      translation = isdigit((int)*cp) ;
    }

  if (translation != 0)
    {
      errno = 0;
      retval = strtol(num_str, &endptr, 10) ;
      if (retval == 0 && (errno != 0 || num_str == endptr))
	{
	  result = FALSE ;
	}
      else if (errno !=0 && (retval == LONG_MAX || retval == LONG_MIN))
	{
	  result = FALSE ;
	}
      else
	{
	  result = TRUE ;
	  *num_out = retval ;
	}
    }
  
  return result ;
}
#endif

/************************************************************/

char *strcasestr (const char *str1, const char *str2)
{
  char *s1 = strnew (str1, 0) ;
  char *s2 = strnew (str2, 0) ;

  bufferToUpper (s1) ;
  bufferToUpper (s2) ;

  const char *cp = strstr(s1, s2) ;
  if (cp)
    {
      int n = cp - s1 ;
      cp = str1 + n ;
    }
  ac_free (s1) ;
  ac_free (s2) ;
  return (char *) cp ; 
}


/*************************************************************************/
/*************************************************************************/

/* It's sometimes necessary to quote parts of the text so it doesn't get     */
/* interpreted too literally (i.e. you may want to keep \n as "\n"). The     */
/* aceInProtect() routine will do this. But then you will probably need to   */
/* unquote the text at some time.                                            */
/*                                                                           */
/* Quoting protects text by putting \" at the start and end of the text and  */
/* \\ in front of any special chars.                                         */
/*    - sometimes you want to completely reverse this - use aceInUnprotect() */
/*       to do this.                                                         */
/*    - sometimes you just want to remove the leading and trailing \", use   */
/*       aceInUnprotectQuote() to do this.                                   */
/*                                                                           */

static char *ac_doUnprotect (const char *text, BOOL just_quotes, AC_HANDLE h)
{
  char *result ;
  char *cp, *cq ; 
  const char *ccp = text ;

  if (! text) 
    return 0 ; /* empty string */ 

  while (*ccp == '\"' && *ccp == '\t' && *ccp == ' ') ccp++ ;
  if (!*ccp)
    return  halloc (1, h) ; /* empty string */

  result = halloc (1 + strlen(text), h) ; /* ensure long enough */
  
  /* remove leading white space, first quotes and any more white space.      */
  ccp = text ;
  while (*ccp == ' ' || *ccp == '\t') ccp++ ;
  if (*ccp == '"') ccp++ ;
  while (*ccp == ' ' || *ccp == '\t') ccp++ ;
  strcpy (result, ccp) ;

  /* remove trailing white space, last quotes and any more white space.     */
  cp = result ;
  cq = cp + strlen(cp) - 1 ;
  while (cq > cp && (*cq == ' ' || *cq == '\t')) *cq-- = 0 ;
  if (*cq == '"') /* remove one unprotected quote */
    {
      int i = 0 ; char *cr = cq - 1 ;
      while (cr > cp && *cr == '\\')
	{ i++ ; cr-- ; }
      if ( i%2 == 0)
	*cq-- = 0 ;  /* discard */
    }
  while (cq > cp && (*cq == ' ' || *cq == '\t')) *cq-- = 0 ;


  /* Optionally gobble the \ as well as the " */
  if (just_quotes) ;
  else
    {
      cq = cp = result ; cp-- ;
      while (*++cp)
	switch (*cp)
	  {
	  case '\\': 
	    if (*(cp+1) == '\\') { cp++ ; *cq++ = '\\' ; break ;}
	    if (*(cp+1) == '\n') { cp ++ ; break ; } /* skip backslash-newline */
	    if (*(cp+1) == 'n') { cp ++ ; *cq++ = '\n' ; break ; }
	    break ;
	  default: *cq++ = *cp ;
	  }
      *cq = 0 ;   /* terminate the string */
    }

  return result ; /* we must return result, so that ac_free(result) is a licit operation */
} /* ac_doUnprotect */

/****************************************/

char* ac_unprotect (const char *text, AC_HANDLE h)
{
  return  ac_doUnprotect(text, FALSE, h) ;
} /* ac_unprotectquote */

/****************************************/

char *ac_unprotectquote (const char *text, AC_HANDLE h)
{
  return  ac_doUnprotect(text, TRUE, h) ;
} /* ac_unprotect */

/****************************************/
/* acedb will read result back as text */
char* ac_protect (const char* text, AC_HANDLE h) 
{
  char *result ;
  const char *ccp = text ;
  const char *cp ; char *cq ;
  
  if (! text)
    return  0 ; /* null string */
  
  while (*ccp == '\t' && *ccp == ' ') ccp++ ;
  if (!*ccp)
    return  strnew("\"\"", h) ; /* emptyprotected  string */
  
  result = halloc ( 64 + 2*(1+strlen(text)), h) ; /* ensure long enough: 2020_05_31, was 2*() before, but failed an AQUILA test.escape */
  
  cq = result ;
  *cq++ = '"' ;
  for (cp = text ; *cp ; )
    { 
      if (*cp == '\\' || *cp == '"' || 		       /* protect these */
	  *cp == '/' || *cp == '%' || *cp == ';' ||
	  *cp == '\t'
	  /* || *cp == '\n' NO this is done on next line mieg:2017-11-21 */
	  )
	*cq++ = '\\' ;
      if (*cp == '\n') 
	{ /* do not acedump a \n, bad for other scripts */
	  *cq++ = '\\' ; *cq++ = 'n' ; cp++ ;
	} 
      else
	*cq++ = *cp++ ;
    }
  *cq++ = '"' ;
  *cq = 0 ;

  return result ;
}  /* ac_protect */

/***************************************************************/
/* Wrapper for lexstrcmp to allow it to be used with arraySort */
/* ie for sorting arrays of character pointers.                */
/* lexstrcmp does case-insensitive, intuitive comparison of    */
/* potentially mixed alpha/numeric strings, but requires char* */
/* parameters, while arraySort ultimately utilises qsort which */
/* requires void * parameters.                                 */
/***************************************************************/
int arrstrcmp(const void *s1, const void *s2)
{
  const char *a = *(const char **)s1;
  const char *b = *(const char **)s2;

  return lexstrcmp(a, b);
}

/**************************************************************/
/* Correctly sorts anything containing integers               */
/* lexRename relies on the fact that lexstrcmp must fail      */
/* if the length differ                                       */
/**************************************************************/

static BOOL lexPureNumber (const char *a, long int *x)
{
  int i, k, n = strlen (a) ;
  long int u = 0 ;
  const char *p = a ;
  switch ((int)*p)
    {
    case '+' : p++ ; n-- ; break ;
    case '-' : p++ ; n-- ; break ;
    }
  for (i = 0 ; i < n ; p++, i++)
    {
      k = (*p) - '0' ;
      if (k>=0 && k <=9)
	u = 10 * u + k ;
      else
	return FALSE ;
      if (u < 0)
	return FALSE ;
    }
  if (*a == '-')
    u = -u ;
  *x = u ;
  return TRUE ;
} /* lesPureNumber */

/**********/

int lexstrcmp(const char *a, const char *b)
{ 
  register char c,d ;
  register const char *p,*q ;
  register int  nbza, nbzb ; /* nb de zeros en tete */
  register int  nbzReturn = 0 ;
  long x ,y ;
  
  if (lexPureNumber (a, &x) && lexPureNumber (b, &y) )
    { /* mieg added 20250302 :    -7 should sort before +2 */
      if (a < b) return -1 ;
      if (a > b) return 1 ;
      if  (strlen(a) > strlen(b)) return -1 ;
      if  (strlen(a) < strlen(b)) return 1 ;
      return 0 ;
    }
  while (*a)
    {                /* Bond007< Bond07 < Bond7 < Bond68 */
      if (isdigit((int)*a) && isdigit((int)*b))
        { 
	  for (nbza = 0 ; *a == '0' ; ++a, nbza++) {} ;  /* saut des premiers zeros */
          for (nbzb = 0 ; *b == '0' ; ++b, nbzb++) {} ;
          for (p = a ; isdigit((int)*p) ; ++p) ;
          for (q = b ; isdigit((int)*q) ; ++q) ;
          if (p-a > q-b) { return 1 ; }  /* the longer number is the bigger */
          if (p-a < q-b) { return -1; }
          while (isdigit ((int)*a))
            { 
	      if (*a > *b)     { return 1 ; }
              if (*a++ < *b++) { return -1; }
            }
          if (!nbzReturn)
            { 
	      if (nbza < nbzb) nbzReturn = +1 ;
              if (nbza > nbzb) nbzReturn = -1 ;
            }
        }
      else
        { 
	  if((c=ace_upper(*a++))>(d=ace_upper(*b++)))
            { return 1 ; }
          if(c<d) { return -1 ; }
        }
    }
 if (!*b)
   { return nbzReturn ; }

 return -1 ;
}

/************************************************************/
/* Same as above but CASE SENSITIVE                         */
/* Correctly sorts anything containing integers             */
/* lexRename relies on the fact that lexstrcmp must fail    */
/*  if the length differ                                    */
/************************************************************/
int lexstrCasecmp(const char *a,const char *b)
{ 
  register char c,d ;
  register const char *p,*q;
  register int  nbza = 0, nbzb = 0 ; /* nb de zeros en tete */

  while (*a)
    {                /* Bond007< Bond07 < Bond7 < Bond68 */
      if (isdigit((int)*a) && isdigit((int)*b))
	{ 
	  nbza = nbza * 10 ; /* pour plusieurs series de 0 */
	  nbzb = nbzb * 10 ;
	  for (p = a ; *p == '0' ; ++p, nbza++) ;  /* saut des premiers zeros*/
	  a = p ;
	  for (q = b ; *q == '0' ; ++q, nbzb++) ;
	  b = q ;
	  for (p = a ; isdigit((int)*p) ; ++p) ;
	  for (q = b ; isdigit((int)*q) ; ++q) ;
	  if (p-a > q-b) { return 1 ; }  /* the longer number is the bigger */
	  if (p-a < q-b) { return -1; }
	  while (isdigit ((int)*a))
	    { 
	      if (*a > *b)     { return 1 ; }
	      if (*a++ < *b++) { return -1; }
	    }
	}
      else
        { 
	  if((c=(*a++))>(d=(*b++))) { return 1 ; }
	  if(c<d) { return -1 ; }
	}
    }
 if (!*b) 
   { 
     if (nbza < nbzb) { return +1 ; } 
     if (nbza > nbzb) { return -1 ; }
     return 0 ;
   }
 return -1 ;
}


/************************************************************/
/* mieg, 2003, measure the number of bytes for a vararg list */
int utPrintfSizeOfArgList (const char * formatDescription, va_list args)
{
  int len = 0, lenAdd, width, stopCountingWidth, isLong, nArg = 0 ;
  const char *fmt, *typeOf;
  void *pVal;
  /* we do not care if we slightly overestimate */
#define	textSIZEOFSLASH		1	 /* should be 1 logically */
#define	textSIZEOFINT   	128	 /* should be width or strlen(sprintf(intnumber)) */
#define	textSIZEOFREAL		128	 /* see above ... (example: %.2f 1.234e78 uses 81 bytes) */
#define	textSIZEOFADDRESS	18 /* should be (2 + sizeof(void*)*2) each byte is two hexadecimal characters */
#define textSIZEOFCHAR		1	/* 1 */
#define textSIZEOFDEFAULT	1	/* 1 */
  
  
  if (!formatDescription)
    return 0 ;
  
  for (fmt = formatDescription;*fmt;fmt++)
    {
      if(*fmt == '%')
	{ /* scan for format specification */
	  
	  /* see the width of the variable, example: %10.3lf */
	  width = 0; stopCountingWidth = 0;
	  for (typeOf = fmt + 1; *typeOf && *typeOf != '%' && !isalpha((int)(*typeOf)); typeOf++ )
	    {
	      
				/* stop counting if dot is enountered after some numbers */
	      if(*typeOf == '.')
		stopCountingWidth = 1; 
				/* accumulate the width */
	      if(isdigit((int)(*typeOf)) && (!stopCountingWidth))
		width = width*10+(*typeOf)-'0';
	    }
	  
	  if (*typeOf  ==  'l')
	    { isLong = 1; typeOf++;} /* long format */
	  else
	    isLong = 0 ;
	  switch (tolower((int)*typeOf))
	    {
	    case 'i':case 'u': case 'd': case 'x':
	      if( isLong) va_arg(args, long int);
	      else va_arg(args, int );
	      lenAdd = textSIZEOFINT;
	      break;
	    case 'f': case 'g': case 'e':
	      if (isLong) va_arg(args, double );
	      else va_arg(args, double );
	      lenAdd = textSIZEOFREAL;
	      break;	
	    case 's':
	      pVal = va_arg(args, char * );
	      lenAdd = strlen((char*)pVal);
	      break;
	    case 'p': 
	      pVal = va_arg(args, void * );
	      lenAdd = textSIZEOFADDRESS;
	      break;
	    case 'c': 
	      va_arg(args, int );
	      lenAdd = textSIZEOFCHAR;
	      break;
	    default: /* %% or alike */
	      lenAdd = textSIZEOFDEFAULT;
	      break;
	    }
	  /* never trust a shorter width because number are NOT trucated
	   * exampe %4.1f, 100.0 needs 5 bytes, not 4
	   */
	  if (width > lenAdd)
	    lenAdd = width;
	  fmt = typeOf;
	  len += lenAdd; 
          nArg++ ;
	}
      else 
	len++;
    }
  
  return len + 256 ;  /* overestimate, not very costly and safer */
} /* utPrintfSizeOfArgList */

int utPrintfSizeOf(char * formatDescription , ...)
{
  int len = 0 ;
  va_list args ;
  
  va_start (args, formatDescription) ;
  len =  utPrintfSizeOfArgList (formatDescription, args) ;
  va_end(args);
  
  return len;
}  /* utPrintfSizeOf */


/********************************************************************/
/********************************************************************/
/* utility to find and consume an argument on the unix command line */
/* get the existence of the arg */
BOOL getArg (int *argcp, const char **argv, const char *arg)
{
  int i = *argcp ;

  for (i = 1 ; i < *argcp ; i++)
    if (argv[i] && !strcmp (arg, argv[i]))
      {
	for (;i < *argcp - 1 ; i++)
	  argv[i] = argv[i+1] ;
	*argcp -= 1 ;
	return TRUE ;
      }
  return FALSE ;
}
/*************************************************************************************/
/* get the value of the arg */
const char* getArgV (int *argcp, const char **argv, const char *arg)
{
  int i = *argcp ;
  const char *cp = 0 ;

  for (i = 1 ; i < *argcp - 1 ; i++)
    if (argv[i] && !strcmp (arg, argv[i]))
      {
	cp = argv[i+1] ;
	for (;i < *argcp - 2 ; i++)
	  argv[i] = argv[i+2] ;
	*argcp -= 2 ;
	break ;	
      }
  return cp ;
}

/*******************************************************************/
/*******************************************************************/

/* match to template with wildcards.   Authorized wildchars are * ? #
     ? represents any single char
     * represents any set of chars
   case sensitive.   Example: *Nc*DE# fits abcaNchjDE23

   returns 0 if not found
           1 + pos of first sigificant match (i.e. not a *) if found
*/

static int doPickMatch (const char *cp,const  char *tp, BOOL caseSensitive)
{
  const char *c=cp, *t=tp;
  const char *ts = 0, *cs = 0, *s = 0 ;
  int star=0;

  while (TRUE)
    switch(*t)
      {
      case '\0':
 /*
        return (!*c ? ( s ? 1 + (s - cp) : 1) : 0) ;
*/
	if(!*c)
	  return  ( s ? 1 + (s - cp) : 1) ;
	if (!star)
	  return 0 ;
        /* else not success yet go back in template */
	t=ts; c=cs+1;
	if(ts == tp) s = 0 ;
	break ;
      case '?':
	if (!*c)
	  return 0 ;
	if(!s) s = c ;
        t++ ;  c++ ;
        break;
      case '*':
        ts=t;
        while( *t == '?' || *t == '*')
          t++;
        if (!*t)
          return s ? 1 + (s-cp) : 1 ;
	if (!caseSensitive)
	  {
	    while (ace_upper(*c) != ace_upper(*t))
	      if(*c)
		c++;
	      else
		return 0 ;
	  }
	else
	  {
	    while (*c != *t)
	      if(*c)
		c++;
	      else
		return 0 ;
	  }
        star=1;
        cs=c;
	if(!s) s = c ;
        break;
      default:
	if (!caseSensitive)
	  {      
	    if (ace_upper(*t++) != ace_upper(*c++))
	      { if(!star)
		return 0 ;
	      t=ts; c=cs+1;
	      if(ts == tp) s = 0 ;
	      }
	    else
	      if(!s) s = c - 1 ;
	  }
	else
	  {
	    if (*t++ != *c++)
	      { if(!star)
		return 0 ;
	      t=ts; c=cs+1;
	      if(ts == tp) s = 0 ;
	      }
	    else
	      if(!s) s = c - 1 ;
	  }
        break;
      }
} /* pickMatch */

int pickMatch (const char *cp, const char *tp)
{ return doPickMatch (cp,tp,FALSE) ; }
int pickMatchCaseSensitive (const  char *cp,const char *tp, BOOL caseSensitive)
{ return doPickMatch (cp,tp,caseSensitive) ; }

/*************************************************************************************/
/* exchangle lower/upper case */
void bufferSwitchCase (char *buf)
{
  register  char *cp ;
  for (cp = buf ; *cp ; cp++)
    {
      if (*cp == ace_upper(*cp))
	*cp= ace_lower(*cp) ;
      else
	*cp = ace_upper(*cp) ;
    }
}  /* bufferSwitchCase */

/*************************************************************************************/

void bufferToUpper (char *buf)
{
  register  char *cp ;
  for (cp = buf ; *cp ; cp++)
    *cp = ace_upper(*cp) ;
  return ;
}  /*  bufferToUpper */

/*************************************************************************************/

void bufferToLower (char *buf)
{
  register  char *cp ;
  for (cp = buf ; *cp ; cp++)
    *cp = ace_lower(*cp) ;
  return ;
}  /*  bufferToLower */

/*************************************************************************************/

int intOrder (const void *a, const void *b)
{
  int x1 = *(const int *) a, x2 = *(const int *) b ;
  return x1 < x2 ? -1 : (x1 > x2 ? 1 : 0) ;
}

/*************************************************************************************/

int unsignedIntOrder (const void *a, const void *b)
{
  int x1 = *(const unsigned int *) a, x2 = *(const unsigned int *) b ;
  return x1 < x2 ? -1 : (x1 > x2 ? 1 : 0) ;
}

/*************************************************************************************/

int unsignedIntReverseOrder (const void *a, const void *b)
{
  int x1 = *(const unsigned int *) a, x2 = *(const unsigned int *) b ;
  return x1 > x2 ? -1 : (x1 < x2 ? 1 : 0) ;
}

/*************************************************************************************/

int floatOrder (const void *a, const void *b)
{
  float x1 = *(const float *) a, x2 = *(const float *) b ;
  return x1 < x2 ? -1 : (x1 > x2 ? 1 : 0) ;
}

/*************************************************************************************/

int doubleOrder (const void *a, const void *b)
{
  double x1 = *(const double *) a, x2 = *(const double *) b ;
  return x1 < x2 ? -1 : (x1 > x2 ? 1 : 0) ;
}

/*************************************************************************************/
/* 2009_03_31
 * mieg@ncbi.nlm.nih.gov
 * The purpose of this library is to obtain a very compressed format for deep sequencing
 * tag counts, also called wiggles, attached to every position along the genome
 *
 * The problem of storing deep sequencing profiles is acute.
 * The formats proposed so far use lots of bytes. To be terse yet reasonnably precise
 * we propose to give the value at every base position, so we do not have
 * to indicate the position, saving 4bytes, and to be content with 1%
 * accuracy for each value, which should be sufficient for most applications.
 * Accepting these limitations, we can store all values up to 10^5 in a single
 * byte per position along the genome. Of course we could gain even more space
 * by storing a value only every n base (n=2,5,10 as wished).

 * The trick is to store the exact value for small numbers, then a kind of one byte
 * approximation of the log of the value with a base very close to 1. 
 * The library below shows the proposed encoder, fast decoder, and a test.

 * The test code shows that the average error is indeed 1% in the [1,10^5] range.
 * All comments welcome !

 */

static int oneByteValue[256] ;

/*********************************/

void oneByteInitialize (int showCode)
{
  int i, ix, a, old = 0 ;
  double x, z, logdeux = log(2.0) ;
  
  for (i = 0; i < 256; i++)
    {
      a = 19 ;
      z = exp (logdeux/a) ; 
      x = old * z ;
      ix = x + .4999 ;
      if (i && ix < old+1)
	ix = old + 1 ;
      old = ix ;
      oneByteValue[i] = ix ;
      if (showCode) 
	printf ("%d\t%d\n", i, ix) ;
    }

  return ;
} /* oneByteInialize */

/*********************************/
/* Encoding could be accelerated, but is probably used seldom */
unsigned char oneByteEncode (unsigned int x)
{
  unsigned char c = 255 ;
  int j ;
  
  for (j = 0 ; j < 256 ; j++)
    if(2*x < oneByteValue[j] + oneByteValue[j+1])
      { c = j ; break ; }
	  
  return c ;
} /* oneByteEncode */

/*********************************/
/* Decoding is very fast */
unsigned int oneByteDecode (unsigned char c)
{
  unsigned int x = oneByteValue[(int)c] ;
  return x ;
} /* oneByteEncode */

/*********************************/
/* This test shows that the relative error is between 0 and 2.4%,
 * with average 1.05% in the [0,10^5] range.
 */
int oneByteTest (void)
{
  int i, j, nn ; 
  unsigned int i1 ;
  double x, z, z2, z2max = 0 ;

  oneByteInitialize (1) ; 
  for (nn=0, x=0, i=1 ; i<100000; nn++, i++)
    {
      i1 = i ;
      j = oneByteDecode (oneByteEncode (i1)) ;
      z = (i - j) ; z/=i; z2 = z*z ;
      x += z2 ;
      if (z2 > z2max) z2max = z2 ;
    }
  printf ("Max and mean relative  errors: %g  %g\n", sqrt(z2max), sqrt(x/nn));
  
  return 0 ;
} /* oneByteTest */

/*************************************************************************************/
/* grep the multiplicity in a fastc style sequence name xxxx#421a200b121c100 */
int fastcMultiplicity (const char *ccp, int *mult, int multMax)
{
  int i, j, n = 0, n1 = 0 ;
  const char *cn = ccp + strlen (ccp) ;

  while (cn > ccp && *cn != '#') cn-- ;
  if (cn > ccp)
    {
      for (cn++ ; *cn && *cn >= '0' && *cn <= '9' ; cn++)
	{
	  n *= 10 ;
	  n += *cn - '0' ;
	}
    }
  if (mult && multMax)
    {
      while (cn > ccp && *cn != '#') cn-- ;
      while (cn > ccp && *cn != '!') cn-- ;
      if (*cn++ == '!')
	{
	  while (*cn != '#')
	    {
	      i = *cn++  - 'a' ;
	      j = 0 ;
	      while (*cn >= '0' && *cn <= '9')
		{ j = 10 * j + (*cn - '0') ; cn++ ; }
	      if (i >= 0 && i < multMax)
		mult[i] += j ;
	      n1 += j ;
	    }
	}
    }
  if (n1 > n) n = n1 ; /* a hack, if something went wrong in the name creation */
  return n > 0 ? n : 1 ;
} /* fastcMultiplicity */
 
/********************* eof ***********************************/
#ifdef JUNK
/* copied from  wikipedia: fetch-and-add 
 * this function is specific of the x86 pltforms and allows multithread synchro
*/
inline int fetch_and_add (int *variable, int value)
{
  asm volatile("lock; xaddl %%eax, %2;"
	       :"=a" (value)      // Output
	       :"a" (value), "m" (*variable)  //Input
	       :"memory") ;
  return value ;
} /* fetch_and_add */
#endif

/*************************************************************************************/
/*************************************************************************************/
/* pour construire la table en lazy
 * on choist n1, n2  < N  , w0 = n1 * N * N * N  + n1 * N * N ,
 * pour chaque config i, on trouve un score u1 et on incremente table[w0 + u1] et p2
 * a la fin on balaie pout trouver les cumul du score (N -1)^2 * jusqu'au score 0
 * et on divise par p2, ce qui nous donne la proba de la valeur >= U1
 * et on stoque p2 [w0 - 1], ce qui nous dit que la table a ete cree pour ce couple n1,n2
 * il suffit alors de consulter la table
 */

static void wilcoxonCompute (Array aa, int NN,  int n1, int n2)
{  /* we should tabulate and verify that the exact discrete table converges to the erfc function for large n1, n2 */
  int NN2 = NN * NN ;
  int NN3 = NN2 * NN ;
  int i, j, i1 ;
  int n = n1 + n2 ;
  int N0 = n1 * NN3 + n2 * NN  ;
  int iMax = 1 << n ;
  int p2 = 0 ; /* number of configs with score >= U1, total number of configs */
  int v1, v2 ; /* running number of 1 and of 0 we are beating */
  int u1 ; /* running u1 score (same as AUC) */
  double *up ;
  double cc ;
  
  if (n1 > n2) { i = n1 ; n1 = n2 ; n2 = i ; }  /* swap */
  /* compute all bit chains of size n = n1 + n2 */
  
  up = arrp (aa, N0, double) ;
  memset (up, 0, NN2 * sizeof(double)) ;
  
  p2 = 0 ;
  for (i = 0 ; i < iMax ; i++)
    {
      i1 = i ; 
      v1 = 0 ; v2 = n2 ; u1 = 0 ;
      for (j = 0 ; j < n ; j++, i1 >>= 1)
	{
	  if (i1 & 0x1) 
	    {
	      v1++ ; 
	      u1 += v2 ; /* running AUC */	
	      if (v1 > n1)
		continue ;  /* bad i config */
	    }
	  else
	    {
	      v2-- ;
	      if (v2 < 0)
		continue ; /* bad i config */
	    }
	}
      if (v1 == n1 && v2 == 0) /* found a config of type (n1,n2) */
	{
	  p2++ ;
	  up[u1] ++ ;
	}
    }
  /* compute the cumuls right to left, and divide by the total */
  for (i = NN2, cc = 0 ; i >= 0 ; i--)
    { cc += up[i] ; up[i] = cc/p2 ; }
  up[NN2 - 1] = p2 ; /* store the cu,ul here, since max(u1) = n1*n2 < NN*2 - 1 */
  return ;
}

/*************************************************************************************/
/* Wilcoxon rank sum test statitics
 * http://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
 * Mann–Whitney U test (also called the Mann–Whitney–Wilcoxon (MWW), Wilcoxon rank-sum test (WRS), or Wilcoxon–Mann–Whitney test) 
 * given 2 interdigited sets
 * compute the probability that the rank order of the members as separated as much or more
 *
 */

static void  testWilcoxon (void) ;
BOOL wilcoxon (int U1, int n1, int n2, double *pValuep, double *pGaussValuep) 
{
  int NN = 21 ;  /* the array will be of size NN^4 */
  int NN2 = NN * NN ;
  int NN3 = NN2 * NN ;
  int NN4 = NN3 * NN ;
  int N0 = n1 * NN3 + n2 * NN  ;
  static Array aa = 0 ; /* allow lzy calculation of the exact probabilities up to n1 < N and n2 < N */

  if (U1 == -99 && n1 == -99 && n2 == -99) /* secret convention */
    {
      testWilcoxon () ;
      return TRUE ;
    }

  if (2 * U1 < n1 * n2)
    U1 =  n1 * n2 - U1 ; 

  if (! pValuep)
    messcrash ("wilcoxon received a null pointer pValuep, please edit the source code of the calling function") ;
  *pValuep = -1 ;
  if (pGaussValuep)
    *pGaussValuep = -1 ;

  if (U1 < 0 || n1 <1 || n2 < 1)
    return FALSE ;

  if (!aa)
    {
      aa = arrayCreate (NN4, double) ;
      array (aa, NN4 - 1, double) = 0 ; /* make room */
    }
  /* wilcoxon is related to AUC

    For comparing two small sets of observations, a direct method is quick, and gives insight into the meaning of the U statistic, 
    it corresponds to the number of wins out of all pairwise contests (see the tortoise and hare example under Examples below). 
    For each observation in one set, count the number of times it wins over any observations in the other set. 
    Count 0.5 for any ties. The sum of wins and ties is U1 for the first set. 

    This is exactly the same as the AUC unnormalized surface, each time we go up the surface of the horizontal rectangle is
    we could start in the other direction get U2 and of course U1 + U2 = the surface of the rectangle = n1*n1

    one can also compute the sum of the rank R1 of the first set and use U1 = n1 * n2 + {n1 * (n1+1)/2.0} - R1 
    

    Normal approximation

    For large samples, U is approximately normally distributed, with mean m and  standard deviation sigma. 
    m = n1 * n2 / 2    (i.e., half the size of the rectangle)
       notice that (2 * (U - m)) is the un-normalized Ginni coefficient 
          The fact that the average value of U is half the surface of the rectangle follows from a geometric argument: 
          except the diagonal path (whose surface == half the rectangle)
          each path in the rectangle can be paired with he symmetric path around the center of the rectangle, 
          and the sum of the U of these 2 paths == the whole rectangle. QED.
    sigma = sqrt (n1 * n2 * (n1 + n2 + 1) / 12.0) ;
	  
    In that case, the standardized value

        z = (U - m)/sigma

    is approximately a standard normal deviate whose significance can be checked using the erfc function 
    i verified that I have no scaling error
       if z =  erfc(x/sqrt(2.0)))
          x = 1.96 => z = 0.0499958,  
       i.e. the well known 5% probablibility for x to be in the interval = average +- 1.96 sigma (two-tailed error)
  */
  if (n1 * n2 >= 1) /* always compute the Gauss value, this allows to verify that the tabulated values converge towards th Gauss value for large n */
    {
      double z, m, sigma, sq2 = sqrt((double)2.0) ;
      
      m = n1 * n2 / 2.0 ;
      sigma = sqrt (n1 * n2 * (n1 + n2 + 1) / 12.0) ;
      z = (U1 - m) / sigma ;
      *pValuep = erfc (z/sq2)/2.0 ;  /* divide by 2 to get the one sided tail area */
      if (pGaussValuep)
	*pGaussValuep = *pValuep ;
     }

  if (n1 < NN &&  n2 < NN) /* do not tabulate more than 1 Million cases */
    {
      if (arr (aa, N0 + NN - 1, double) == 0) /* this table has not been computed */
	wilcoxonCompute (aa, NN, n1, n2) ;
      *pValuep = arr (aa, N0 + U1, double) ;
      return TRUE ;  /* exact calculation */
    }
  else  /* numbers too large for a tabulation */
    if (n1 >= 5 && n2 >= 5)
      {
	return TRUE ;  /* trust the Gaussain estimation */
      }
    else /* n1 too small for Gauss approximation, n1 + n2 too large for direct tabulation */
      {
	*pValuep = -1 ;
	return FALSE ;
      }

  /* WARNING: i did not correct for ties 
     The formula for the standard deviation is more complicated in the presence of tied ranks; the full formula is given in the text books referenced below.[citation needed]
     However, if the number of ties is small (and especially if there are no large tie bands) ties can be ignored when doing calculations by hand. The computer statistical packages will use the correctly adjusted formula as a matter of routine.
     
     Note that since U1 + U2 = n1 n2, the mean m = n1 * n2/2 is both the mean of U1 and the mean of U2. Therefore, the absolute value of the z statistic can be computed with U1 or U2 indifferently
  */

  /* Generalization to c dimensions, say analyse at once stages 1,2,3,4,4s
     Because of its probabilistic form, the U statistic can be generalised to a measure of a classifier's separation power for more than two classes:[17]
     
     M = sum (AUC(i,j) )/ c (c-1)        summing over i=1..c  j=1...c   for 2 classes we recover the standard AUC
     the rest of the discussion is stupid, since AUK(k,k) is not zero but nk * nk /2 , so probably we should use Gini not AUC, and it not clear when or if we should divide by ni * nj, i.e. move from U(i,j)
     
     Where c is the number of classes, and the R_{k,l} term of AUC_{k,l} considers only the ranking of the items belonging to classes k and l (i.e., items belonging to all other classes are ignored) according to the classifier's estimates of the probability of those items belonging to class k. AUC_{k,k} will always be zero but, unlike in the two-class case, generally AUC_{k,l} \ne AUC_{l,k}, which is why the M measure sums over all (k, l) pairs, in effect using the average of AUC_{k,l} and AUC_{l,k}.

     Machine Learning, 45, 171–186, 2001
     2001 Kluwer Academic Publishers. Manufactured in The Netherlands.
     A Simple Generalisation of the Area Under the ROC
     Curve for Multiple Class Classification Problems
     DAVID J. HAND d.j.hand@ic.ac.uk
     ROBERT J. TILL
     r.till@ic.ac.uk
     Department of Mathematics, Imperial College, Huxley Building, 180 Queen’s Gate, London SW7 2BZ, UK  
  */
  
  /*
 
    In practice some of this information may already have been supplied and common sense should be used in deciding whether to repeat it. A typical report might run,

    "Median latencies in groups E and C were 153 and 247 ms; the distributions in the two groups differed significantly (Mann–Whitney U = 10.5, n1 = n2 = 8, P < 0.05 two-tailed)."
    
    A statement that does full justice to the statistical status of the test might run,
    
    "Outcomes of the two treatments were compared using the Wilcoxon–Mann–Whitney two-sample rank-sum test. The treatment effect (difference between treatments) was quantified using the Hodges–Lehmann (HL) estimator, which is consistent with the Wilcoxon test.[21] This estimator (HLΔ) is the median of all possible differences in outcomes between a subject in group B and a subject in group A. A non-parametric 0.95 confidence interval for HLΔ accompanies these estimates as does ρ, an estimate of the probability that a randomly chosen subject from population B has a higher weight than a randomly chosen subject from population A. The median [quartiles] weight for subjects on treatment A and B respectively are 147 [121, 177] and 151 [130, 180] kg. Treatment A decreased weight by HLΔ = 5 kg (0.95 CL [2, 9] kg, 2P = 0.02, ρ = 0.58)."

However it would be rare to find so extended a report in a document whose major topic was not statistical inference.
  */
  return FALSE ;
} /* wilcoxon */

/*************************************************************************************/
/* Test the Wilcoxon program, verify in particular that the tabulated values converge towards the Gaussian estimates */
static void  testWilcoxon (void)
{
  double p, pG ;
  int nn, n1, n2, U1 ;
  
  for (nn = 10 ; nn < 15 ; nn++)
    for (n1 = 10 ; n1 < 15 ; n1++)
      {
	printf ("\n%d\t%d", nn, n1) ;
	n2 = nn - n1 ;
	for (U1 = 0 ; U1 <= n1 * n2 ; U1++) 
	  {
	    wilcoxon (U1, n1, n2, &p, &pG) ;
	    printf ("\tU=%d:: %g : G=%g", U1, p, pG) ;
	  }
      }
  printf ("\n") ;
    
  return ;
} /* testWilcoxon  */



/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
#include <regex.h>

struct regExpStruct {
  int magic ;
  regex_t *regex;
  const char *pattern ;
} ;

static void regExpFinalise (void *vp)
{
  RegExp br = (RegExp)vp ;
  
  if (br && br->magic == 623562)
    {
      if (br->regex)
	regfree (br->regex) ;
      ac_free (br->pattern) ;
      br->magic = 0 ;
    }  
}

/*******************************************************************/
/* 0: bad bql, !NULL: use in regExpFind, then call regExpDoFree */
RegExp regExpCreate (const char *pattern, AC_HANDLE h)
{
  int nn ;
  RegExp br ;

  if (! pattern)
    messcrash ("regExpCreate called wiit NULL pattern") ;
  if (! *pattern)
    return 0 ;

  br = (RegExp) halloc (sizeof (struct regExpStruct), h) ;
  if (br)
    {
      br->magic = 623562 ;
      br->regex = halloc (sizeof (regex_t), h) ;
      br->pattern = strnew (pattern, 0) ;

      blockSetFinalise (br, regExpFinalise) ;
          
      /* man regcomp for details */
      if (0)
	nn = regcomp (br->regex, pattern, REG_EXTENDED | REG_ICASE) ;
      else
	nn = regcomp (br->regex, pattern, REG_EXTENDED | REG_ICASE | REG_NOSUB ) ;
      if (nn) /* error */
	ac_free (br) ;
    }
  return br ;
} /* regExpCreate */

/*******************************************************************/

int regExpFind (RegExp br, const char *data)
{
  int ok = -1 ; /* default value : null or empty 'data' */

  if (!br)
    messcrash ("regExpFind calle on NULL RegExp*") ;
  if (br->magic != 623562)
    messcrash ("regExpFind calle on invalid RegExp*") ;
    
  if (data && *data)
    {
      regmatch_t pmatch ;
      int nn = 0 ;

      nn = regexec (br->regex, data, 1, &pmatch, 0) ;
      if (nn) 
	ok = 0 ; /* pattern not found */
      else  /* match found */ 
	{ 	     
	  if (1)
	    ok = 1 + pmatch.rm_so ;
	  else
	    ok = 1 ;
	}
    }
  return ok ;
} /* regExpFind */

/*******************************************************************/

int regExpMatch (const char *data, const char *pattern)
{
  int ok = 0 ;
  
  if (pattern && data && *pattern && *data)
    {
      RegExp br = regExpCreate (pattern, 0) ;
      if (br)
	{
	  ok = regExpFind (br, data) ; /* 0: not found, >0 position */
	  ac_free (br) ;
	}
      else
	ok = -1 ; /* pattern is wrong */
    }
  else
    ok = -2 ;

  return ok ;
} /* regExpMatch  */

/*************************************************************************************/
/*  Check if executable exists in PATH */
BOOL isExecutableInPath (const char *name)
{
  const char *path = getenv("PATH") ;
  if (!path)
    return FALSE ;
  char *path_copy = strdup(path) ;
  if (path_copy)
    {
      char *dir = strtok(path_copy, ":") ;
      while (dir)
	{
	  char full_path[1024] ;
	  snprintf (full_path, sizeof(full_path), "%s/%s", dir, name) ;
	  if (access(full_path, X_OK) == 0)
	    {
	      free (path_copy) ;
	      return TRUE ;
	    }
	  dir = strtok (NULL, ":") ;
	}
      free (path_copy) ;
    }
  return FALSE ;
} /* isExecutableInPath */

/*******************************************************************/
/************************* end of file ****************************/
/*******************************************************************/


