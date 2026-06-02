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
 *   (b) how many CPUs are on that node   → upper bound for thread spawning
 *   (c) how many virtual threads the OS allows this process to register
 *   (d) how much RAM is currently available  → guard before allocating large buffers
 *
 * PUBLIC API  (declare in sa.h)
 * ----------
 *   SA_HardwareInfo saGetHardwareInfo (void) ;
 *       Call once, at startup, before any threads are created.
 *       On Linux:  samples NUMA node load, picks the least-busy node,
 *                  reports CPUs and RAM for that node.
 *       On macOS / other Unix:  no NUMA, reports total CPUs and RAM.
 *
 *   long saGetAvailableRamKb (void) ;
 *       May be called periodically (not in a tight loop — see comment below).
 *       Returns currently available RAM in KB for the process's memory domain.
 *
 * PORTING NOTES
 * -------------
 *   Linux  : reads /proc/stat, /sys/devices/system/node/, /proc/meminfo
 *   macOS  : uses sysctl(3)
 *   Other  : falls back to sysconf(3)  — gives CPUs and physical pages
 *   numactl: only meaningful on Linux multi-node machines; best_node == 0
 *            is always a safe value on single-node or non-Linux builds.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>      /* sysconf, getpid, usleep */

#ifdef __linux__
#  include <dirent.h>
#  include <sys/resource.h>   /* getrlimit / RLIMIT_NPROC */
#  include <sys/types.h>
#  include <time.h>
#endif

#ifdef __APPLE__
#  include <sys/sysctl.h>
#  include <sys/types.h>
#  include <mach/mach.h>      /* host_statistics for free pages */
#endif

/* ------------------------------------------------------------------ */
/* Public struct — fill in saGetHardwareInfo()                         */
/* ------------------------------------------------------------------ */

typedef struct {
    int  best_node ;        /* NUMA node index to pass to numactl (Linux)
                             * Always 0 on macOS / single-node systems.    */
    int  node_cpus ;        /* Number of logical CPUs on best_node.
                             * On non-NUMA systems: total logical CPUs.     */
    int  max_threads ;      /* Max virtual threads this process may register
                             * (from RLIMIT_NPROC on Linux, or a safe
                             *  platform cap elsewhere).  -1 = unlimited.   */
    long available_ram_kb ; /* Available RAM in KB at the moment of the call,
                             * scoped to best_node when possible.           */
} SA_HardwareInfo ;


/* ================================================================== */
/*                     LINUX  IMPLEMENTATION                          */
/* ================================================================== */
#ifdef __linux__

/* ---------- constants -------------------------------------------- */

/* Maximum logical CPU index we will track.  If your machine has more,
 * increase this.  Each cpu_times_t is 80 bytes, so 512 → 80 KB on the
 * heap (we allocate dynamically below, not on the stack).            */
#define SA_MAX_CPUS  512

/* How long to wait between the two /proc/stat samples.
 * 100 ms is enough to distinguish busy from idle nodes.
 * Do not reduce below ~50 ms or the delta becomes noisy.            */
#define SA_SAMPLE_US  100000   /* 100 ms */

/* ---------- internal types --------------------------------------- */

typedef struct {
    unsigned long long user, nice, sys, idle, iowait, irq, softirq, steal ;
} cpu_times_t ;

/* ---------- count NUMA nodes ------------------------------------- */
/*
 * Returns the number of nodeN directories under /sys/devices/system/node.
 * Returns 1 if the directory is absent (non-NUMA kernel or container).
 */
static int count_numa_nodes (void)
{
    DIR *d = opendir ("/sys/devices/system/node") ;
    if (!d) return 1 ;

    int count = 0 ;
    struct dirent *e ;
    while ((e = readdir (d)))
      { int dummy ;
        if (sscanf (e->d_name, "node%d", &dummy) == 1)
            count++ ;
      }
    closedir (d) ;
    return count > 0 ? count : 1 ;
}

/* ---------- CPU list for one NUMA node --------------------------- */
/*
 * Parses /sys/devices/system/node/nodeN/cpulist  (format: "0-7,16-23")
 * Fills cpu_ids[] with the individual CPU indices, returns count.
 * Returns 0 on any error.
 */
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

            char *dash = strchr (tok, '-') ;
            if (dash)
              { int a = atoi (tok), b = atoi (dash + 1) ;
                for (int i = a ; i <= b && n < max_ids ; i++)
                    cpu_ids[n++] = i ;
              }
            else
              { cpu_ids[n++] = atoi (tok) ; }

            tok = comma ? comma + 1 : NULL ;
          }
      }
    fclose (f) ;
    return n ;
}

/* ---------- read /proc/stat -------------------------------------- */
/*
 * Fills times[cpu_id] for every "cpuN" line found.
 * Returns the highest cpu_id seen + 1  (i.e. the CPU count),
 * or -1 on error.
 */
static int read_proc_stat (cpu_times_t *times, int max_cpus)
{
    FILE *f = fopen ("/proc/stat", "r") ;
    if (!f) return -1 ;

    char buf[512] ;
    int  num_cpus = 0 ;

    /* Skip the first "cpu " aggregate line */
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
        if (cpu_id + 1 > num_cpus)
            num_cpus = cpu_id + 1 ;
      }

    fclose (f) ;
    return num_cpus ;
}

/* ---------- per-CPU utilisation from two snapshots --------------- */
/*
 * Returns utilisation in [0.0, 100.0] for a single CPU.
 * Returns 0.0 if the total delta is zero (CPU was halted / offline).
 */
static double cpu_usage (const cpu_times_t *before, const cpu_times_t *after)
{
    unsigned long long active =
          (after->user    - before->user)
        + (after->nice    - before->nice)
        + (after->sys     - before->sys)
        + (after->irq     - before->irq)
        + (after->softirq - before->softirq)
        + (after->steal   - before->steal) ;

    unsigned long long total  = active
        + (after->idle   - before->idle)
        + (after->iowait - before->iowait) ;

    if (total == 0) return 0.0 ;
    return 100.0 * (double) active / (double) total ;
}

/* ---------- max registerable virtual threads --------------------- */
/*
 * Returns the soft RLIMIT_NPROC limit for this process.
 * On Linux this limits total threads+processes for the real UID,
 * which is the relevant ceiling when registering pthreads.
 * Returns -1 if the limit is RLIM_INFINITY (treat as unlimited).
 * Returns  0 on error (caller should use a safe hard-coded fallback).
 */
static int linux_max_threads (void)
{
    struct rlimit rl ;
    if (getrlimit (RLIMIT_NPROC, &rl) != 0) return 0 ;
    if (rl.rlim_cur == RLIM_INFINITY)       return -1 ;
    return (int) rl.rlim_cur ;
}

/* ---------- available RAM for one NUMA node ---------------------- */
/*
 * Reads /sys/devices/system/node/nodeN/meminfo and returns
 * MemFree + Cached + Inactive in KB.
 * Returns 0 if the file is not available.
 *
 * TIME COST: one fopen + sequential line scan of a tiny kernel file,
 * typically < 0.1 ms.  Safe to call every few seconds; not in a tight loop.
 */
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

    /* MemFree alone is too conservative: the kernel will reclaim Cached
     * and Inactive pages on demand, so they count as usable.           */
    return free_kb + cached_kb + inactive_kb ;
}

/* ---------- system-wide available RAM (Linux fallback) ----------- */
/*
 * Used when no NUMA node is identified.
 * MemAvailable is the most accurate single figure for "how much can
 * a new allocation get" — it accounts for reclaimable caches.
 *
 * TIME COST: one fopen + scan of /proc/meminfo (~50 lines), < 0.1 ms.
 */
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

/* ---------- main Linux entry point ------------------------------- */

static SA_HardwareInfo linux_hardware_info (void)
{
    SA_HardwareInfo info = { 0, 1, 0, 0 } ;

    /* ---- 1. How many NUMA nodes?  If only one, nothing to measure. */
    int num_nodes = count_numa_nodes () ;

    if (num_nodes <= 1)
      {
        /* Single node or non-NUMA kernel: skip load sampling entirely. */
        int cpu_ids[SA_MAX_CPUS] ;
        info.best_node      = 0 ;
        info.node_cpus      = node_cpu_list (0, cpu_ids, SA_MAX_CPUS) ;
        if (info.node_cpus == 0)   /* cpulist unavailable, use sysconf */
            info.node_cpus  = (int) sysconf (_SC_NPROCESSORS_ONLN) ;
        info.max_threads    = linux_max_threads () ;
        info.available_ram_kb = node_available_ram_kb (0) ;
        if (info.available_ram_kb == 0)
            info.available_ram_kb = linux_system_available_ram_kb () ;
        return info ;
      }

    /* ---- 2. Allocate snapshot buffers on the heap (not the stack).
     *         Two snapshots × SA_MAX_CPUS × 64 bytes ≈ 64 KB total.  */
    cpu_times_t *before = calloc (SA_MAX_CPUS, sizeof (cpu_times_t)) ;
    cpu_times_t *after  = calloc (SA_MAX_CPUS, sizeof (cpu_times_t)) ;
    if (!before || !after)
      { free (before) ; free (after) ; return info ; }

    /* ---- 3. First snapshot, sleep, second snapshot. */
    int num_cpus = read_proc_stat (before, SA_MAX_CPUS) ;
    if (num_cpus < 1)
      { free (before) ; free (after) ; return info ; }

    usleep (SA_SAMPLE_US) ;

    if (read_proc_stat (after, SA_MAX_CPUS) != num_cpus)
      { free (before) ; free (after) ; return info ; }

    /* ---- 4. Compute per-CPU utilisation over the interval. */
    double usage[SA_MAX_CPUS] ;
    for (int c = 0 ; c < num_cpus ; c++)
        usage[c] = cpu_usage (&before[c], &after[c]) ;

    free (before) ;
    free (after) ;

    /* ---- 5. For each node, average the utilisation of its CPUs.
     *         The node with the lowest average wins.                  */
    int    best_node      = 0 ;
    double min_avg        = 101.0 ;   /* sentinel: above any real value */
    int    best_cpu_ids[SA_MAX_CPUS] ;
    int    best_cpu_count = 1 ;

    int cpu_ids[SA_MAX_CPUS] ;

    for (int node = 0 ; node < num_nodes ; node++)
      {
        int ncpus = node_cpu_list (node, cpu_ids, SA_MAX_CPUS) ;
        if (ncpus == 0) continue ;

        double sum = 0.0 ;
        int    counted = 0 ;
        for (int i = 0 ; i < ncpus ; i++)
          {
            int c = cpu_ids[i] ;
            if (c >= 0 && c < num_cpus)
              { sum += usage[c] ; counted++ ; }
          }
        if (counted == 0) continue ;

        double avg = sum / counted ;
        if (avg < min_avg)
          {
            min_avg = avg ;
            best_node = node ;
            best_cpu_count = ncpus ;
            memcpy (best_cpu_ids, cpu_ids, ncpus * sizeof (int)) ;
            (void) best_cpu_ids ;   /* stored for potential future use */
          }
      }

    /* ---- 6. Fill the result struct. */
    info.best_node        = best_node ;
    info.node_cpus        = best_cpu_count ;
    info.max_threads      = linux_max_threads () ;
    info.available_ram_kb = node_available_ram_kb (best_node) ;
    if (info.available_ram_kb == 0)
        info.available_ram_kb = linux_system_available_ram_kb () ;

    return info ;
}

#endif  /* __linux__ */


/* ================================================================== */
/*                     macOS  IMPLEMENTATION                          */
/* ================================================================== */
#ifdef __APPLE__

/*
 * macOS has no numactl and no per-node memory.
 * We report total logical CPUs and available RAM.
 * best_node is always 0.
 *
 * TIME COST of sysctl calls: < 0.05 ms each.
 */

static int macos_cpu_count (void)
{
    int n = 0 ;
    size_t len = sizeof (n) ;
    /* HW_LOGICALCPU counts hyperthreads; use HW_PHYSICALCPU for cores only.
     * For thread-count decisions, logical CPUs is the right figure.       */
    if (sysctlbyname ("hw.logicalcpu", &n, &len, NULL, 0) == 0 && n > 0)
        return n ;
    /* Fallback: HW_NCPU via MIB */
    int mib[2] = { CTL_HW, HW_NCPU } ;
    if (sysctl (mib, 2, &n, &len, NULL, 0) == 0 && n > 0)
        return n ;
    return 1 ;
}

/*
 * Available RAM on macOS via Mach VM statistics.
 * free_count + inactive_count gives a realistic "usable" figure;
 * inactive pages are reclaimable by the kernel on demand.
 *
 * TIME COST: one mach trap, < 0.1 ms.  Safe to call every few seconds.
 */
static long macos_available_ram_kb (void)
{
    vm_statistics64_data_t vm ;
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT ;
    if (host_statistics64 (mach_host_self (), HOST_VM_INFO64,
                           (host_info64_t) &vm, &count) == KERN_SUCCESS)
      {
        vm_size_t page = vm_page_size ;
        unsigned long long usable =
            ((unsigned long long) vm.free_count
           + (unsigned long long) vm.inactive_count) * page ;
        return (long) (usable / 1024ULL) ;
      }

    /* Fallback: conservative 70% of total physical RAM */
    int64_t total = 0 ;
    size_t  len   = sizeof (total) ;
    if (sysctlbyname ("hw.memsize", &total, &len, NULL, 0) == 0)
        return (long) ((unsigned long long) total * 70 / 100 / 1024) ;

    return 0 ;
}

static SA_HardwareInfo macos_hardware_info (void)
{
    SA_HardwareInfo info ;
    info.best_node        = 0 ;
    info.node_cpus        = macos_cpu_count () ;
    info.max_threads      = -1 ;   /* macOS imposes no useful RLIMIT_NPROC */
    info.available_ram_kb = macos_available_ram_kb () ;
    return info ;
}

#endif  /* __APPLE__ */


/* ================================================================== */
/*                  GENERIC UNIX  FALLBACK                            */
/* ================================================================== */
#if !defined(__linux__) && !defined(__APPLE__)

/*
 * sysconf is POSIX and works on most Unix variants (AIX, Solaris, BSDs).
 * We get logical CPU count and physical page count; no per-node info.
 *
 * TIME COST: negligible kernel call, < 0.01 ms.
 */
static SA_HardwareInfo generic_hardware_info (void)
{
    SA_HardwareInfo info ;
    info.best_node   = 0 ;
    info.node_cpus   = (int) sysconf (_SC_NPROCESSORS_ONLN) ;
    if (info.node_cpus < 1) info.node_cpus = 1 ;
    info.max_threads = -1 ;   /* unknown; caller should use own cap */

    long pages     = sysconf (_SC_AVPHYS_PAGES) ;   /* available pages */
    long page_size = sysconf (_SC_PAGE_SIZE) ;
    info.available_ram_kb =
        (pages > 0 && page_size > 0)
        ? (long) ((unsigned long long) pages * page_size / 1024)
        : 0 ;

    return info ;
}

#endif


/* ================================================================== */
/*                     PUBLIC  FUNCTIONS                              */
/* ================================================================== */

/*
 * saGetHardwareInfo
 * -----------------
 * Call ONCE at startup, before spawning any threads.
 *
 * On Linux with multiple NUMA nodes:
 *   - samples /proc/stat twice over SA_SAMPLE_US microseconds
 *   - picks the node with the lowest average CPU load
 *   - the caller should then re-exec the program under:
 *       numactl --cpunodebind=<best_node> --membind=<best_node> ...
 *
 * On Linux with a single node, macOS, or other Unix:
 *   - no load sampling, returns immediately
 *   - best_node is always 0 (numactl call is a no-op / not needed)
 *
 * Fields in the returned struct:
 *   best_node        — NUMA node index (0 on non-NUMA)
 *   node_cpus        — logical CPUs available for threading
 *   max_threads      — OS limit on virtual thread registration;
 *                      -1 means unlimited; 0 means query failed
 *                      (suggest capping at 256 if -1 or 0)
 *   available_ram_kb — RAM available right now, in KB
 */
SA_HardwareInfo saGetHardwareInfo (void)
{
#if defined(__linux__)
    return linux_hardware_info () ;
#elif defined(__APPLE__)
    return macos_hardware_info () ;
#else
    return generic_hardware_info () ;
#endif
}


/*
 * saGetAvailableRamKb
 * -------------------
 * Returns currently available RAM in KB.
 * May be called periodically to decide whether to launch another thread
 * or allocate a large buffer.
 *
 * TIME COST:
 *   Linux : opens and scans /proc/meminfo (~50 short lines) — typically
 *           0.02–0.10 ms.  Calling once per second is perfectly safe;
 *           once every 10–100 ms is acceptable under load.
 *           Do NOT call in a tight inner loop (millions of times).
 *   macOS : one Mach trap — < 0.05 ms.  Same guidance applies.
 *   Other : two sysconf() calls — negligible.
 *
 * Note: the returned value is an estimate.  The kernel may reclaim
 * cached pages at any time, so treat it as a planning heuristic,
 * not a guarantee.
 */
long saGetAvailableRamKb (void)
{
#if defined(__linux__)
    /* Try to stay consistent with the node chosen in saGetHardwareInfo.
     * We re-read /proc/self/status to find the current Mems_allowed node.
     * If the process has been bound to a single node via numactl, this
     * returns that node; otherwise it falls back to system-wide meminfo.  */
    int node = -1 ;
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
