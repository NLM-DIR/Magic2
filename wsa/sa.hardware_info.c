/* ================================================================== */
/*                     GPU  DETECTION                                 */
/* ================================================================== */
/*
 * saGetGpuInfo
 * ------------
 * Probes for CUDA-capable GPU devices via dlopen() and returns the
 * count.  Optionally fills *bestDevice with the index of the device
 * having the most free memory, suitable for passing to cudaSetDevice().
 *
 * The CUDA runtime is loaded at run time via dlopen() so this binary
 * has no link-time dependency on libcuda.  The function is therefore
 * safe to call on machines with no CUDA driver installed.
 *
 * Returns:
 *    > 0   number of CUDA devices found; *bestDevice set if non-NULL
 *      0   no CUDA-capable device, or driver / runtime not found
 *     -1   runtime found but device query failed (driver mismatch,
 *          device in exclusive-process mode, etc.); treat as 0 for
 *          respawn decisions, but consider logging a warning.
 *
 * Called exactly once, during the hardware-detection phase in
 * sa.main.c, before any threads are created.
 *
 * RESPAWN USAGE
 * -------------
 *   int bestDev = -1 ;
 *   int ngpu    = saGetGpuInfo (&bestDev) ;
 *   if (ngpu > 0)
 *     {
 *       // build "numactl ... magic2_gpu --numactl --gpu-device=N ..."
 *       // bestDev carries the selected device index
 *       return system (cmd) ;
 *     }
 *
 * Pass NULL for bestDevice if you only need presence/absence:
 *   if (saGetGpuInfo (NULL) > 0)
 *       use_gpu_binary () ;
 */

typedef cudaError_t (*pfn_cudaGetDeviceCount_t) (int *) ;
typedef cudaError_t (*pfn_cudaSetDevice_t)      (int) ;
typedef cudaError_t (*pfn_cudaMemGetInfo_t)     (size_t *, size_t *) ;

#if defined(__linux__) || defined(__APPLE__)
#  include <dlfcn.h>
#endif

int saGetGpuInfo (int *bestDevice)
{
    if (bestDevice) *bestDevice = -1 ;

#if !defined(__linux__) && !defined(__APPLE__)
    return 0 ;
#else

    /* ---- locate the CUDA runtime shared library ------------------- */
    void *libcuda = NULL ;
#ifdef __linux__
    libcuda = dlopen ("libcudart.so.12",  RTLD_LAZY | RTLD_LOCAL) ;
    if (!libcuda)
        libcuda = dlopen ("libcudart.so.11.0", RTLD_LAZY | RTLD_LOCAL) ;
    if (!libcuda)
        libcuda = dlopen ("libcudart.so",      RTLD_LAZY | RTLD_LOCAL) ;
#else  /* __APPLE__ */
    libcuda = dlopen ("libcudart.dylib",  RTLD_LAZY | RTLD_LOCAL) ;
#endif
    if (!libcuda)
        return 0 ;   /* no CUDA runtime installed — not an error */

    /* ---- resolve the three entry points we need ------------------- */
    pfn_cudaGetDeviceCount_t fn_count =
        (pfn_cudaGetDeviceCount_t) dlsym (libcuda, "cudaGetDeviceCount") ;
    pfn_cudaSetDevice_t fn_set =
        (pfn_cudaSetDevice_t)      dlsym (libcuda, "cudaSetDevice") ;
    pfn_cudaMemGetInfo_t fn_mem =
        (pfn_cudaMemGetInfo_t)     dlsym (libcuda, "cudaMemGetInfo") ;

    if (!fn_count || !fn_set || !fn_mem)
      { /* Library present but entry points missing — should not happen. */
        dlclose (libcuda) ;
        return -1 ;
      }

    /* ---- query device count --------------------------------------- */
    int ndev = 0 ;
    cudaError_t rc = fn_count (&ndev) ;
    if (rc != 0 || ndev <= 0)
      { dlclose (libcuda) ;
        return (rc == 0) ? 0 : -1 ;
      }

    /* ---- pick device with most free memory ------------------------ */
    int    best      = 0 ;
    size_t best_free = 0 ;

    for (int d = 0 ; d < ndev ; d++)
      {
        if (fn_set (d) != 0) continue ;   /* skip inaccessible device */
        size_t free_bytes = 0, total_bytes = 0 ;
        if (fn_mem (&free_bytes, &total_bytes) == 0)
            if (free_bytes > best_free)
              { best_free = free_bytes ; best = d ; }
      }

    /* Leave library handle open: process is about to exec() anyway.  */
    if (bestDevice) *bestDevice = best ;
    return ndev ;

#endif  /* __linux__ || __APPLE__ */
}
/*
 * sa.hardware_info.c
 *
 * Part of the sortalign package — RNA aligner
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg, Greg Boratyn  NCBI/NLM/NIH
 *
 * PURPOSE
 * -------
 * Interrogate the hardware once at program startup, before any threading,
 * to decide:
 *   (a) which NUMA node is least loaded  → passed to numactl on Linux
 *   (b) how many cpus, threads and RAM are available
 * See sa.hardware_info.c for the public API.
 *
 * CALL SEQUENCE (sa.main.c)
 * -------------------------
 *
 *  ── first invocation of the program ──────────────────────────────────
 *
 *  1.  saGetBestNumaNode()          [Linux only, guard with #ifdef __linux__]
 *        Measures CPU load across NUMA nodes (~100 ms).
 *        Returns the index of the least-loaded node, or -1 when there
 *        is only one node.  The caller re-execs via numactl only when
 *        the return value is >= 0.  If -1, continue without re-exec.
 *
 *  ── second invocation (or first on non-NUMA / non-Linux) ─────────────
 *
 *  2.  saGetNodeCpus()              [all platforms]
 *        Returns the number of logical CPUs available to this process.
 *        After numactl binding: CPUs of the bound node.
 *        Without numactl: total logical CPUs on the machine.
 *        No sleep; returns immediately.
 *
 *  3.  saGetMaxThreads()            [all platforms]
 *        Returns the OS ceiling on virtual thread registration for this
 *        process.  -1 means unlimited; use your own cap (e.g. 256).
 *        No sleep; returns immediately.
 *
 *  4.  saGetAvailableRamKb()        [all platforms, re-callable]
 *        Returns currently available RAM in KB, scoped to the bound NUMA
 *        node when possible.  Safe to call every few seconds.
 *        Do NOT call in a tight inner loop.
 *
 * PORTING
 * -------
 *   Linux  : /proc/stat, /sys/devices/system/node/, /proc/meminfo
 *   macOS  : sysctl(3), host_statistics64 (Mach)
 *   Other  : sysconf(3) — POSIX, works on AIX, Solaris, BSDs
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>      /* sysconf, usleep */

#ifdef __linux__
#  include <dirent.h>
#  include <sys/resource.h>   /* getrlimit / RLIMIT_NPROC */
#  include <sys/types.h>
#endif

#ifdef __APPLE__
#  include <sys/sysctl.h>
#  include <mach/mach.h>
#endif


/* ================================================================== */
/*                     LINUX  INTERNALS                               */
/* ================================================================== */
#ifdef __linux__

/* Maximum logical CPU index we will track.  Each cpu_times_t is 64 bytes,
 * so 512 entries → two heap buffers of 32 KB each during sampling.   */
#define SA_MAX_CPUS  512

/* Interval between the two /proc/stat samples.
 * 100 ms gives a clean signal without a perceptible startup delay.
 * Do not reduce below ~50 ms or the delta becomes noisy.             */
#define SA_SAMPLE_US  200000 /* 200 ms */

typedef struct {
    unsigned long long user, nice, sys, idle, iowait, irq, softirq, steal ;
} cpu_times_t ;

/* ------------------------------------------------------------------ */
/* Count nodeN entries under /sys/devices/system/node.
 * Returns 1 if the directory is absent (non-NUMA kernel, container).  */
static int count_numa_nodes (void)
{
    DIR *d = opendir ("/sys/devices/system/node") ;
    if (!d) return 1 ;
    int count = 0, dummy ;
    struct dirent *e ;
    while ((e = readdir (d)))
        if (sscanf (e->d_name, "node%d", &dummy) == 1)
            count++ ;
    closedir (d) ;
    return count > 0 ? count : 1 ;
}

/* ------------------------------------------------------------------ */
/* Parse /sys/devices/system/node/nodeN/cpulist  ("0-7,16-23" format).
 * Fills cpu_ids[], returns the count.  Returns 0 on any error.        */
static int node_cpu_list (int node, int *cpu_ids, int max_ids)
{
    char path[128] ;
    snprintf (path, sizeof (path),
              "/sys/devices/system/node/node%d/cpulist", node) ;
    FILE *f = fopen (path, "r") ;
    if (!f) return 0 ;

    char line[4096] ;
    int  n = 0 ;
    if (fgets (line, sizeof (line), f))
      {
        char *tok = line ;
        while (tok && n < max_ids)
          {
            char *comma = strchr (tok, ',') ;
            if (comma) *comma = '\0' ;
            char *dash  = strchr (tok, '-') ;
            if (dash)
              { int a = atoi (tok), b = atoi (dash + 1) ;
                for (int i = a ; i <= b && n < max_ids ; i++)
                    cpu_ids[n++] = i ;
              }
            else if (*tok >= '0' && *tok <= '9')
                cpu_ids[n++] = atoi (tok) ;
            tok = comma ? comma + 1 : NULL ;
          }
      }
    fclose (f) ;
    return n ;
}

/* ------------------------------------------------------------------ */
/* Read per-CPU cumulative jiffies from /proc/stat.
 * Returns the highest cpu index seen + 1, or -1 on error.            */
static int read_proc_stat (cpu_times_t *times, int max_cpus)
{
    FILE *f = fopen ("/proc/stat", "r") ;
    if (!f) return -1 ;

    char buf[512] ;
    int  num_cpus = 0 ;
    /* Skip the aggregate "cpu " (no digit) summary line. */
    if (!fgets (buf, sizeof (buf), f))
      { fclose (f) ; return -1 ; }

    while (fgets (buf, sizeof (buf), f))
      {
        int cpu_id = -1 ;
        cpu_times_t t = {0} ;
        int n = sscanf (buf,
            "cpu%d %llu %llu %llu %llu %llu %llu %llu %llu",
            &cpu_id,
            &t.user, &t.nice, &t.sys, &t.idle,
            &t.iowait, &t.irq, &t.softirq, &t.steal) ;
        if (n < 9 || cpu_id < 0 || cpu_id >= max_cpus) continue ;
        times[cpu_id] = t ;
        if (cpu_id + 1 > num_cpus) num_cpus = cpu_id + 1 ;
      }
    fclose (f) ;
    return num_cpus ;
}

/* ------------------------------------------------------------------ */
/* CPU utilisation [0..100] computed from two snapshots for one CPU.  */
static double cpu_usage (const cpu_times_t *a, const cpu_times_t *b)
{
    unsigned long long active =
          (b->user    - a->user)
        + (b->nice    - a->nice)
        + (b->sys     - a->sys)
        + (b->irq     - a->irq)
        + (b->softirq - a->softirq)
        + (b->steal   - a->steal) ;
    unsigned long long total = active
        + (b->idle   - a->idle)
        + (b->iowait - a->iowait) ;
    return total == 0 ? 0.0 : 100.0 * (double) active / (double) total ;
}

/* ------------------------------------------------------------------ */
/* Available RAM for one NUMA node, from /sys/.../nodeN/meminfo.
 * Returns MemFree + Cached + Inactive in KB; 0 on error.
 * Cached and Inactive pages are reclaimable on demand, so they count
 * as usable for planning purposes.                                    */
static long node_available_ram_kb (int node)
{
    char path[128] ;
    snprintf (path, sizeof (path),
              "/sys/devices/system/node/node%d/meminfo", node) ;
    FILE *f = fopen (path, "r") ;
    if (!f) return 0 ;

    char line[256] ;
    long free_kb = 0, cached_kb = 0, inactive_kb = 0 ;
    while (fgets (line, sizeof (line), f))
      {
        if (strstr (line, "MemFree:"))  sscanf (line, "%*s %ld", &free_kb) ;
        if (strstr (line, "Cached:"))   sscanf (line, "%*s %ld", &cached_kb) ;
        if (strstr (line, "Inactive:")) sscanf (line, "%*s %ld", &inactive_kb) ;
      }
    fclose (f) ;
    return free_kb + cached_kb + inactive_kb ;
}

/* ------------------------------------------------------------------ */
/* System-wide available RAM from /proc/meminfo (Linux fallback).
 * MemAvailable is the most accurate figure: it accounts for caches
 * the kernel will reclaim when a real allocation arrives.             */
static long linux_system_available_ram_kb (void)
{
    FILE *f = fopen ("/proc/meminfo", "r") ;
    if (!f) return 0 ;
    char line[256] ;
    long avail = 0 ;
    while (fgets (line, sizeof (line), f))
        if (strncmp (line, "MemAvailable:", 13) == 0)
          { sscanf (line + 13, "%ld", &avail) ; break ; }
    fclose (f) ;
    return avail ;
}

#endif  /* __linux__ */


/* ================================================================== */
/*                    macOS  INTERNALS                                */
/* ================================================================== */
#ifdef __APPLE__

/* Available RAM via Mach VM statistics.
 * free_count + inactive_count is the realistic "usable" figure:
 * inactive pages are reclaimed by the kernel on demand.              */
static long macos_available_ram_kb (void)
{
    vm_statistics64_data_t vm ;
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT ;
    if (host_statistics64 (mach_host_self (), HOST_VM_INFO64,
                           (host_info64_t) &vm, &count) == KERN_SUCCESS)
      {
        unsigned long long usable =
            ((unsigned long long) vm.free_count
           + (unsigned long long) vm.inactive_count)
           * (unsigned long long) vm_page_size ;
        return (long) (usable / 1024ULL) ;
      }
    /* Fallback: conservative 70% of total physical RAM via sysctl.   */
    int64_t total = 0 ;
    size_t  len   = sizeof (total) ;
    if (sysctlbyname ("hw.memsize", &total, &len, NULL, 0) == 0)
        return (long) ((unsigned long long) total * 70 / 100 / 1024) ;
    return 0 ;
}

#endif  /* __APPLE__ */


/* ================================================================== */
/*                     PUBLIC  FUNCTIONS                              */
/* ================================================================== */

/*
 * saGetBestNumaNode  [Linux only]
 * --------------------------------
 * Samples /proc/stat twice over ~100 ms, computes average CPU load per
 * NUMA node, and returns the index of the least-loaded node.
 *
 * Returns:
 *    >= 0   node index  → caller should re-exec under:
 *               numactl --cpunodebind=<N> --membind=<N> ./program ...
 *    -1     only one NUMA node, or node topology unreadable:
 *               caller should skip numactl and continue directly.
 *
 * This function exists only on Linux.  Wrap call sites with:
 *    #ifdef __linux__
 *    int node = saGetBestNumaNode () ;
 *    if (node >= 0) { ... numactl re-exec ... }
 *    #endif
 */
#ifdef __linux__
int saGetBestNumaNode (void)
{
    int num_nodes = count_numa_nodes () ;
    if (num_nodes <= 1)
        return -1 ;   /* nothing to choose between */

    /* Snapshot buffers on the heap — ~64 KB total, not on the stack. */
    cpu_times_t *before = calloc (SA_MAX_CPUS, sizeof (cpu_times_t)) ;
    cpu_times_t *after  = calloc (SA_MAX_CPUS, sizeof (cpu_times_t)) ;
    if (!before || !after)
      { free (before) ; free (after) ; return -1 ; }

    int num_cpus = read_proc_stat (before, SA_MAX_CPUS) ;
    if (num_cpus < 1)
      { free (before) ; free (after) ; return -1 ; }

    usleep (SA_SAMPLE_US) ;    /* 100 ms — only cost of this function */

    if (read_proc_stat (after, SA_MAX_CPUS) < 1)
      { free (before) ; free (after) ; return -1 ; }

    double usage[SA_MAX_CPUS] ;
    for (int c = 0 ; c < num_cpus ; c++)
        usage[c] = cpu_usage (&before[c], &after[c]) ;
    free (before) ;
    free (after) ;

    /* Find the node with the lowest mean CPU utilisation. */
    int    best_node = 0 ;
    double min_avg   = 101.0 ;   /* sentinel above any real value */
    int    cpu_ids[SA_MAX_CPUS] ;

    for (int node = 0 ; node < num_nodes ; node++)
      {
        int ncpus = node_cpu_list (node, cpu_ids, SA_MAX_CPUS) ;
        if (ncpus == 0) continue ;

        double sum = 0.0 ;
        int    counted = 0 ;
        for (int i = 0 ; i < ncpus ; i++)
          { int c = cpu_ids[i] ;
            if (c >= 0 && c < num_cpus)
              { sum += usage[c] ; counted++ ; }
          }
        if (counted == 0) continue ;

        double avg = sum / counted ;
        if (avg < min_avg)
          { min_avg = avg ; best_node = node ; }
      }

    return best_node ;
}
#endif  /* __linux__ */


/*
 * saGetNodeCpus
 * -------------
 * Returns the number of logical CPUs available to this process.
 * No sleep; returns immediately.
 *
 * On Linux:
 *   Reads Cpus_allowed_list from /proc/self/status.  After
 *   "numactl --cpunodebind=N" this reflects exactly the CPUs of node N.
 *   Without numactl it reflects all online CPUs.
 *
 * On macOS:
 *   Uses hw.logicalcpu via sysctl (counts hyperthreads).
 *
 * Fallback (other Unix):
 *   sysconf(_SC_NPROCESSORS_ONLN).
 *
 * Always returns at least 1.
 */
int saGetNodeCpus (void)
{
#ifdef __linux__
    int  n = 0 ;
    FILE *f = fopen ("/proc/self/status", "r") ;
    if (f)
      {
        char line[256] ;
        while (fgets (line, sizeof (line), f))
          {
            if (strncmp (line, "Cpus_allowed_list:", 18) != 0) continue ;
            /* Parse the range string "0-7" or "0-7,16-23" etc. */
            char *p = line + 18 ;
            while (*p == ' ' || *p == '\t') p++ ;
            while (p && *p >= '0' && *p <= '9')
              {
                char *comma = strchr (p, ',') ;
                if (comma) *comma = '\0' ;
                char *dash  = strchr (p, '-') ;
                if (dash)
                    n += atoi (dash + 1) - atoi (p) + 1 ;
                else
                    n++ ;
                p = comma ? comma + 1 : NULL ;
              }
            break ;
          }
        fclose (f) ;
      }
    if (n > 0) return n ;
    n = (int) sysconf (_SC_NPROCESSORS_ONLN) ;
    return n > 0 ? n : 1 ;

#elif defined(__APPLE__)
    int    n   = 0 ;
    size_t len = sizeof (n) ;
    if (sysctlbyname ("hw.logicalcpu", &n, &len, NULL, 0) == 0 && n > 0)
        return n ;
    int mib[2] = { CTL_HW, HW_NCPU } ;
    if (sysctl (mib, 2, &n, &len, NULL, 0) == 0 && n > 0)
        return n ;
    return 1 ;

#else
    int n = (int) sysconf (_SC_NPROCESSORS_ONLN) ;
    return n > 0 ? n : 1 ;
#endif
}


/*
 * saGetMaxThreads
 * ---------------
 * Returns the OS ceiling on virtual thread registrations for this process
 * (i.e. the number of pthreads you can create before pthread_create fails).
 *
 * This is NOT the recommended concurrency level.  Use saGetNodeCpus() and
 * saGetAvailableRamKb() to decide how many threads to actually run.
 *
 * Returns:
 *    > 0    soft RLIMIT_NPROC limit (Linux)
 *    -1     unlimited or not determinable; apply your own cap (e.g. 256)
 *
 * No sleep; returns immediately.
 */
int saGetMaxThreads (void)
{
#ifdef __linux__
    struct rlimit rl ;
    if (getrlimit (RLIMIT_NPROC, &rl) != 0) return -1 ;
    if (rl.rlim_cur == RLIM_INFINITY)       return -1 ;
    int n = rl.rlim_cur ;
    return n < 256 ? n : 256 ;
#else
    /* On macOS, RLIMIT_NPROC limits processes for the UID, not threads
     * per process — not a useful figure here.                          */
    return -1 ;
#endif
}


/*
 * saGetAvailableRamKb
 * -------------------
 * Returns currently available RAM in KB for this process's memory domain.
 *
 * On Linux after "numactl --membind=N":
 *   Reads /sys/devices/system/node/nodeN/meminfo (MemFree + Cached +
 *   Inactive).  Cached and Inactive pages are reclaimable on demand.
 *
 * On Linux without numactl:
 *   Reads MemAvailable from /proc/meminfo.
 *
 * On macOS:
 *   Uses Mach host_statistics64 (free + inactive pages × page size).
 *
 * Fallback:
 *   sysconf(_SC_AVPHYS_PAGES) × page size.
 *
 * TIME COST:
 *   Linux : one fopen + sequential scan of a small kernel file.
 *           Typically 0.02–0.10 ms.
 *   macOS : one Mach trap.  Typically < 0.05 ms.
 *   Other : two sysconf() calls.  Negligible.
 *
 * CALLING FREQUENCY:
 *   Safe every few seconds.  Do NOT call in a tight inner loop.
 *   The value is an estimate — treat it as a planning heuristic.
 */
long saGetAvailableRamKb (void)
{
#ifdef __linux__
    /* Detect which NUMA node we are memory-bound to (if any). */
    int   node = -1 ;
    FILE *f = fopen ("/proc/self/status", "r") ;
    if (f)
      { char line[256] ;
        while (fgets (line, sizeof (line), f))
            if (strncmp (line, "Mems_allowed_list:", 18) == 0)
              { sscanf (line + 18, " %d", &node) ; break ; }
        fclose (f) ;
      }
    if (node >= 0)
      { long kb = node_available_ram_kb (node) ;
        if (kb > 0) return kb ;
      }
    return linux_system_available_ram_kb () ;

#elif defined(__APPLE__)
    return macos_available_ram_kb () ;

#else
    long pages     = sysconf (_SC_AVPHYS_PAGES) ;
    long page_size = sysconf (_SC_PAGE_SIZE) ;
    if (pages > 0 && page_size > 0)
        return (long) ((unsigned long long) pages * page_size / 1024) ;
    return 0 ;
#endif
}
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/*
 * sa.hardware_info.c
 *
 * Part of the sortalign package — RNA aligner
 * Authors: Jean Thierry-Mieg, Danielle Thierry-Mieg, Greg Boratyn  NCBI/NLM/NIH
 *
 * PURPOSE
 * -------
 * Interrogate the hardware once at program startup, before any threading,
 * to decide:
 *   (a) which NUMA node is least loaded  → passed to numactl on Linux
 *   (b) how many cpus, threads and RAM are available
 * See sa.hardware_info.c for the public API.
 *
 * CALL SEQUENCE (sa.main.c)
 * -------------------------
 *
 *  ── first invocation of the program ──────────────────────────────────
 *
 *  1.  saGetBestNumaNode()          [Linux only, guard with #ifdef __linux__]
 *        Measures CPU load across NUMA nodes (~100 ms).
 *        Returns the index of the least-loaded node, or -1 when there
 *        is only one node.  The caller re-execs via numactl only when
 *        the return value is >= 0.  If -1, continue without re-exec.
 *
 *  ── second invocation (or first on non-NUMA / non-Linux) ─────────────
 *
 *  2.  saGetNodeCpus()              [all platforms]
 *        Returns the number of logical CPUs available to this process.
 *        After numactl binding: CPUs of the bound node.
 *        Without numactl: total logical CPUs on the machine.
 *        No sleep; returns immediately.
 *
 *  3.  saGetMaxThreads()            [all platforms]
 *        Returns the OS ceiling on virtual thread registration for this
 *        process.  -1 means unlimited; use your own cap (e.g. 256).
 *        No sleep; returns immediately.
 *
 *  4.  saGetAvailableRamKb()        [all platforms, re-callable]
 *        Returns currently available RAM in KB, scoped to the bound NUMA
 *        node when possible.  Safe to call every few seconds.
 *        Do NOT call in a tight inner loop.
 *
 * PORTING
 * -------
 *   Linux  : /proc/stat, /sys/devices/system/node/, /proc/meminfo
 *   macOS  : sysctl(3), host_statistics64 (Mach)
 *   Other  : sysconf(3) — POSIX, works on AIX, Solaris, BSDs
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>      /* sysconf, usleep */

#ifdef __linux__
#  include <dirent.h>
#  include <sys/resource.h>   /* getrlimit / RLIMIT_NPROC */
#  include <sys/types.h>
#endif

#ifdef __APPLE__
#  include <sys/sysctl.h>
#  include <mach/mach.h>
#endif


/* ================================================================== */
/*                     LINUX  INTERNALS                               */
/* ================================================================== */
#ifdef __linux__

/* Maximum logical CPU index we will track.  Each cpu_times_t is 64 bytes,
 * so 512 entries → two heap buffers of 32 KB each during sampling.   */
#define SA_MAX_CPUS  512

/* Interval between the two /proc/stat samples.
 * 100 ms gives a clean signal without a perceptible startup delay.
 * Do not reduce below ~50 ms or the delta becomes noisy.             */
#define SA_SAMPLE_US  200000 /* 200 ms */

typedef struct {
    unsigned long long user, nice, sys, idle, iowait, irq, softirq, steal ;
} cpu_times_t ;

/* ------------------------------------------------------------------ */
/* Count nodeN entries under /sys/devices/system/node.
 * Returns 1 if the directory is absent (non-NUMA kernel, container).  */
static int count_numa_nodes (void)
{
    DIR *d = opendir ("/sys/devices/system/node") ;
    if (!d) return 1 ;
    int count = 0, dummy ;
    struct dirent *e ;
    while ((e = readdir (d)))
        if (sscanf (e->d_name, "node%d", &dummy) == 1)
            count++ ;
    closedir (d) ;
    return count > 0 ? count : 1 ;
}

/* ------------------------------------------------------------------ */
/* Parse /sys/devices/system/node/nodeN/cpulist  ("0-7,16-23" format).
 * Fills cpu_ids[], returns the count.  Returns 0 on any error.        */
static int node_cpu_list (int node, int *cpu_ids, int max_ids)
{
    char path[128] ;
    snprintf (path, sizeof (path),
              "/sys/devices/system/node/node%d/cpulist", node) ;
    FILE *f = fopen (path, "r") ;
    if (!f) return 0 ;

    char line[4096] ;
    int  n = 0 ;
    if (fgets (line, sizeof (line), f))
      {
        char *tok = line ;
        while (tok && n < max_ids)
          {
            char *comma = strchr (tok, ',') ;
            if (comma) *comma = '\0' ;
            char *dash  = strchr (tok, '-') ;
            if (dash)
              { int a = atoi (tok), b = atoi (dash + 1) ;
                for (int i = a ; i <= b && n < max_ids ; i++)
                    cpu_ids[n++] = i ;
              }
            else if (*tok >= '0' && *tok <= '9')
                cpu_ids[n++] = atoi (tok) ;
            tok = comma ? comma + 1 : NULL ;
          }
      }
    fclose (f) ;
    return n ;
}

/* ------------------------------------------------------------------ */
/* Read per-CPU cumulative jiffies from /proc/stat.
 * Returns the highest cpu index seen + 1, or -1 on error.            */
static int read_proc_stat (cpu_times_t *times, int max_cpus)
{
    FILE *f = fopen ("/proc/stat", "r") ;
    if (!f) return -1 ;

    char buf[512] ;
    int  num_cpus = 0 ;
    /* Skip the aggregate "cpu " (no digit) summary line. */
    if (!fgets (buf, sizeof (buf), f))
      { fclose (f) ; return -1 ; }

    while (fgets (buf, sizeof (buf), f))
      {
        int cpu_id = -1 ;
        cpu_times_t t = {0} ;
        int n = sscanf (buf,
            "cpu%d %llu %llu %llu %llu %llu %llu %llu %llu",
            &cpu_id,
            &t.user, &t.nice, &t.sys, &t.idle,
            &t.iowait, &t.irq, &t.softirq, &t.steal) ;
        if (n < 9 || cpu_id < 0 || cpu_id >= max_cpus) continue ;
        times[cpu_id] = t ;
        if (cpu_id + 1 > num_cpus) num_cpus = cpu_id + 1 ;
      }
    fclose (f) ;
    return num_cpus ;
}

/* ------------------------------------------------------------------ */
/* CPU utilisation [0..100] computed from two snapshots for one CPU.  */
static double cpu_usage (const cpu_times_t *a, const cpu_times_t *b)
{
    unsigned long long active =
          (b->user    - a->user)
        + (b->nice    - a->nice)
        + (b->sys     - a->sys)
        + (b->irq     - a->irq)
        + (b->softirq - a->softirq)
        + (b->steal   - a->steal) ;
    unsigned long long total = active
        + (b->idle   - a->idle)
        + (b->iowait - a->iowait) ;
    return total == 0 ? 0.0 : 100.0 * (double) active / (double) total ;
}

/* ------------------------------------------------------------------ */
/* Available RAM for one NUMA node, from /sys/.../nodeN/meminfo.
 * Returns MemFree + Cached + Inactive in KB; 0 on error.
 * Cached and Inactive pages are reclaimable on demand, so they count
 * as usable for planning purposes.                                    */
static long node_available_ram_kb (int node)
{
    char path[128] ;
    snprintf (path, sizeof (path),
              "/sys/devices/system/node/node%d/meminfo", node) ;
    FILE *f = fopen (path, "r") ;
    if (!f) return 0 ;

    char line[256] ;
    long free_kb = 0, cached_kb = 0, inactive_kb = 0 ;
    while (fgets (line, sizeof (line), f))
      {
        if (strstr (line, "MemFree:"))  sscanf (line, "%*s %ld", &free_kb) ;
        if (strstr (line, "Cached:"))   sscanf (line, "%*s %ld", &cached_kb) ;
        if (strstr (line, "Inactive:")) sscanf (line, "%*s %ld", &inactive_kb) ;
      }
    fclose (f) ;
    return free_kb + cached_kb + inactive_kb ;
}

/* ------------------------------------------------------------------ */
/* System-wide available RAM from /proc/meminfo (Linux fallback).
 * MemAvailable is the most accurate figure: it accounts for caches
 * the kernel will reclaim when a real allocation arrives.             */
static long linux_system_available_ram_kb (void)
{
    FILE *f = fopen ("/proc/meminfo", "r") ;
    if (!f) return 0 ;
    char line[256] ;
    long avail = 0 ;
    while (fgets (line, sizeof (line), f))
        if (strncmp (line, "MemAvailable:", 13) == 0)
          { sscanf (line + 13, "%ld", &avail) ; break ; }
    fclose (f) ;
    return avail ;
}

#endif  /* __linux__ */


/* ================================================================== */
/*                    macOS  INTERNALS                                */
/* ================================================================== */
#ifdef __APPLE__

/* Available RAM via Mach VM statistics.
 * free_count + inactive_count is the realistic "usable" figure:
 * inactive pages are reclaimed by the kernel on demand.              */
static long macos_available_ram_kb (void)
{
    vm_statistics64_data_t vm ;
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT ;
    if (host_statistics64 (mach_host_self (), HOST_VM_INFO64,
                           (host_info64_t) &vm, &count) == KERN_SUCCESS)
      {
        unsigned long long usable =
            ((unsigned long long) vm.free_count
           + (unsigned long long) vm.inactive_count)
           * (unsigned long long) vm_page_size ;
        return (long) (usable / 1024ULL) ;
      }
    /* Fallback: conservative 70% of total physical RAM via sysctl.   */
    int64_t total = 0 ;
    size_t  len   = sizeof (total) ;
    if (sysctlbyname ("hw.memsize", &total, &len, NULL, 0) == 0)
        return (long) ((unsigned long long) total * 70 / 100 / 1024) ;
    return 0 ;
}

#endif  /* __APPLE__ */


/* ================================================================== */
/*                     PUBLIC  FUNCTIONS                              */
/* ================================================================== */

/*
 * saGetBestNumaNode  [Linux only]
 * --------------------------------
 * Samples /proc/stat twice over ~100 ms, computes average CPU load per
 * NUMA node, and returns the index of the least-loaded node.
 *
 * Returns:
 *    >= 0   node index  → caller should re-exec under:
 *               numactl --cpunodebind=<N> --membind=<N> ./program ...
 *    -1     only one NUMA node, or node topology unreadable:
 *               caller should skip numactl and continue directly.
 *
 * This function exists only on Linux.  Wrap call sites with:
 *    #ifdef __linux__
 *    int node = saGetBestNumaNode () ;
 *    if (node >= 0) { ... numactl re-exec ... }
 *    #endif
 */
#ifdef __linux__
int saGetBestNumaNode (void)
{
    int num_nodes = count_numa_nodes () ;
    if (num_nodes <= 1)
        return -1 ;   /* nothing to choose between */

    /* Snapshot buffers on the heap — ~64 KB total, not on the stack. */
    cpu_times_t *before = calloc (SA_MAX_CPUS, sizeof (cpu_times_t)) ;
    cpu_times_t *after  = calloc (SA_MAX_CPUS, sizeof (cpu_times_t)) ;
    if (!before || !after)
      { free (before) ; free (after) ; return -1 ; }

    int num_cpus = read_proc_stat (before, SA_MAX_CPUS) ;
    if (num_cpus < 1)
      { free (before) ; free (after) ; return -1 ; }

    usleep (SA_SAMPLE_US) ;    /* 100 ms — only cost of this function */

    if (read_proc_stat (after, SA_MAX_CPUS) < 1)
      { free (before) ; free (after) ; return -1 ; }

    double usage[SA_MAX_CPUS] ;
    for (int c = 0 ; c < num_cpus ; c++)
        usage[c] = cpu_usage (&before[c], &after[c]) ;
    free (before) ;
    free (after) ;

    /* Find the node with the lowest mean CPU utilisation. */
    int    best_node = 0 ;
    double min_avg   = 101.0 ;   /* sentinel above any real value */
    int    cpu_ids[SA_MAX_CPUS] ;

    for (int node = 0 ; node < num_nodes ; node++)
      {
        int ncpus = node_cpu_list (node, cpu_ids, SA_MAX_CPUS) ;
        if (ncpus == 0) continue ;

        double sum = 0.0 ;
        int    counted = 0 ;
        for (int i = 0 ; i < ncpus ; i++)
          { int c = cpu_ids[i] ;
            if (c >= 0 && c < num_cpus)
              { sum += usage[c] ; counted++ ; }
          }
        if (counted == 0) continue ;

        double avg = sum / counted ;
        if (avg < min_avg)
          { min_avg = avg ; best_node = node ; }
      }

    return best_node ;
}
#endif  /* __linux__ */


/*
 * saGetNodeCpus
 * -------------
 * Returns the number of logical CPUs available to this process.
 * No sleep; returns immediately.
 *
 * On Linux:
 *   Reads Cpus_allowed_list from /proc/self/status.  After
 *   "numactl --cpunodebind=N" this reflects exactly the CPUs of node N.
 *   Without numactl it reflects all online CPUs.
 *
 * On macOS:
 *   Uses hw.logicalcpu via sysctl (counts hyperthreads).
 *
 * Fallback (other Unix):
 *   sysconf(_SC_NPROCESSORS_ONLN).
 *
 * Always returns at least 1.
 */
int saGetNodeCpus (void)
{
#ifdef __linux__
    int  n = 0 ;
    FILE *f = fopen ("/proc/self/status", "r") ;
    if (f)
      {
        char line[256] ;
        while (fgets (line, sizeof (line), f))
          {
            if (strncmp (line, "Cpus_allowed_list:", 18) != 0) continue ;
            /* Parse the range string "0-7" or "0-7,16-23" etc. */
            char *p = line + 18 ;
            while (*p == ' ' || *p == '\t') p++ ;
            while (p && *p >= '0' && *p <= '9')
              {
                char *comma = strchr (p, ',') ;
                if (comma) *comma = '\0' ;
                char *dash  = strchr (p, '-') ;
                if (dash)
                    n += atoi (dash + 1) - atoi (p) + 1 ;
                else
                    n++ ;
                p = comma ? comma + 1 : NULL ;
              }
            break ;
          }
        fclose (f) ;
      }
    if (n > 0) return n ;
    n = (int) sysconf (_SC_NPROCESSORS_ONLN) ;
    return n > 0 ? n : 1 ;

#elif defined(__APPLE__)
    int    n   = 0 ;
    size_t len = sizeof (n) ;
    if (sysctlbyname ("hw.logicalcpu", &n, &len, NULL, 0) == 0 && n > 0)
        return n ;
    int mib[2] = { CTL_HW, HW_NCPU } ;
    if (sysctl (mib, 2, &n, &len, NULL, 0) == 0 && n > 0)
        return n ;
    return 1 ;

#else
    int n = (int) sysconf (_SC_NPROCESSORS_ONLN) ;
    return n > 0 ? n : 1 ;
#endif
}


/*
 * saGetMaxThreads
 * ---------------
 * Returns the OS ceiling on virtual thread registrations for this process
 * (i.e. the number of pthreads you can create before pthread_create fails).
 *
 * This is NOT the recommended concurrency level.  Use saGetNodeCpus() and
 * saGetAvailableRamKb() to decide how many threads to actually run.
 *
 * Returns:
 *    > 0    soft RLIMIT_NPROC limit (Linux)
 *    -1     unlimited or not determinable; apply your own cap (e.g. 256)
 *
 * No sleep; returns immediately.
 */
int saGetMaxThreads (void)
{
#ifdef __linux__
    struct rlimit rl ;
    if (getrlimit (RLIMIT_NPROC, &rl) != 0) return -1 ;
    if (rl.rlim_cur == RLIM_INFINITY)       return -1 ;
    int n = rl.rlim_cur ;
    return n < 256 ? n : 256 ;
#else
    /* On macOS, RLIMIT_NPROC limits processes for the UID, not threads
     * per process — not a useful figure here.                          */
    return -1 ;
#endif
}


/*
 * saGetAvailableRamKb
 * -------------------
 * Returns currently available RAM in KB for this process's memory domain.
 *
 * On Linux after "numactl --membind=N":
 *   Reads /sys/devices/system/node/nodeN/meminfo (MemFree + Cached +
 *   Inactive).  Cached and Inactive pages are reclaimable on demand.
 *
 * On Linux without numactl:
 *   Reads MemAvailable from /proc/meminfo.
 *
 * On macOS:
 *   Uses Mach host_statistics64 (free + inactive pages × page size).
 *
 * Fallback:
 *   sysconf(_SC_AVPHYS_PAGES) × page size.
 *
 * TIME COST:
 *   Linux : one fopen + sequential scan of a small kernel file.
 *           Typically 0.02–0.10 ms.
 *   macOS : one Mach trap.  Typically < 0.05 ms.
 *   Other : two sysconf() calls.  Negligible.
 *
 * CALLING FREQUENCY:
 *   Safe every few seconds.  Do NOT call in a tight inner loop.
 *   The value is an estimate — treat it as a planning heuristic.
 */
long saGetAvailableRamKb (void)
{
#ifdef __linux__
    /* Detect which NUMA node we are memory-bound to (if any). */
    int   node = -1 ;
    FILE *f = fopen ("/proc/self/status", "r") ;
    if (f)
      { char line[256] ;
        while (fgets (line, sizeof (line), f))
            if (strncmp (line, "Mems_allowed_list:", 18) == 0)
              { sscanf (line + 18, " %d", &node) ; break ; }
        fclose (f) ;
      }
    if (node >= 0)
      { long kb = node_available_ram_kb (node) ;
        if (kb > 0) return kb ;
      }
    return linux_system_available_ram_kb () ;

#elif defined(__APPLE__)
    return macos_available_ram_kb () ;

#else
    long pages     = sysconf (_SC_AVPHYS_PAGES) ;
    long page_size = sysconf (_SC_PAGE_SIZE) ;
    if (pages > 0 && page_size > 0)
        return (long) ((unsigned long long) pages * page_size / 1024) ;
    return 0 ;
#endif
}
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
