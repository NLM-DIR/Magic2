
#ifndef HARWARE_INFO_H
#define HARWARE_INFO_H

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

SA_HardwareInfo saGetHardwareInfo (void) ; /* hardware information */
long saGetAvailableRamKb (void) ; /* currently available RAM in KB */

#endif
