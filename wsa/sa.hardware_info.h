
#ifndef HARWARE_INFO_H
#define HARWARE_INFO_H
/*
 * sa.harware_info.h
 *
 * Public API for sa.utils.c — hardware interrogation for sortalign.
 * Include this header wherever these functions are called.
 *
 * TYPICAL CALL SEQUENCE IN sa.main.c
 * ------------------------------------
 *
 *   ── first invocation ─────────────────────────────────────────────
 *
 *   #ifdef __linux__
 *   if (numactl_is_available && !already_under_numactl)
 *     {
 *       int node = saGetBestNumaNode () ;   // ~100 ms
 *       if (node >= 0)
 *         {
 *           // build "numactl --cpunodebind=N --membind=N ./prog ..." 
 *           // and re-exec via system() or execv()
 *           return system (cmd) ;
 *         }
 *       // node == -1: single-node machine, fall through
 *     }
 *   #endif
 *
 *   ── second invocation (or first on non-NUMA / non-Linux) ─────────
 *
 *   int   ncpus      = saGetNodeCpus () ;       // immediate
 *   int   maxthreads = saGetMaxThreads () ;      // immediate
 *   long  ramkb      = saGetAvailableRamKb () ;  // immediate
 *
 *   // ... later, before spawning a new thread or large allocation:
 *   if (saGetAvailableRamKb () > MIN_RAM_PER_THREAD_KB)
 *       spawn_thread () ;
 */


/*
 * saGetBestNumaNode  [Linux only]
 *
 * Measures CPU load across NUMA nodes (~100 ms) and returns the index
 * of the least-loaded node.
 *
 * Return value:
 *    >= 0   pass this to:  numactl --cpunodebind=<N> --membind=<N>
 *    -1     only one NUMA node present, or topology unreadable;
 *           skip the numactl re-exec and continue directly.
 *
 * Call only on Linux, only when numactl is available, and only before
 * any threads have been created.
 */
#ifdef __linux__
int saGetBestNumaNode (void) ;
#endif

/*
 * saGetNodeCpus  [all platforms]
 *
 * Returns the number of logical CPUs available to this process.
 *   - After numactl binding : CPUs of the bound node only.
 *   - Without numactl       : all online CPUs on the machine.
 *   - macOS / other Unix    : total logical CPUs.
 *
 * No sleep.  Always returns >= 1.
 * Use this value as the upper bound when deciding thread count.
 */
int saGetNodeCpus (void) ;

/*
 * saGetMaxThreads  [all platforms]
 *
 * Returns the OS ceiling on virtual thread registrations for this
 * process (the limit before pthread_create would fail).
 *
 * This is NOT the recommended concurrency level — use saGetNodeCpus()
 * and saGetAvailableRamKb() to decide that.
 *
 * Return value:
 *    > 0    soft RLIMIT_NPROC limit (Linux)
 *    -1     unlimited or not determinable; apply your own cap
 *           (256 is a safe default for most RNA-aligner workloads)
 *
 * No sleep.
 */
int saGetMaxThreads (void) ;

/*
 * saGetAvailableRamKb  [all platforms, re-callable]
 *
 * Returns currently available RAM in KB, scoped to the NUMA node the
 * process is bound to when possible.
 *
 * May be called periodically to decide whether to launch a new thread
 * or allocate a large buffer.  Safe every few seconds.
 * Do NOT call in a tight inner loop.
 *
 * The value is an estimate — the kernel may reclaim cached pages at
 * any time.  Treat it as a planning heuristic.
 *
 * Typical cost: 0.02–0.10 ms on Linux; < 0.05 ms on macOS.
 */
long saGetAvailableRamKb (void) ;


#endif  /* HARWARE_INFO_H */
