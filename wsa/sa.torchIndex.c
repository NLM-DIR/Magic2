/*
 * sa_torch_build.c
 *
 * Offline index builder for the LibTorch torcherator genome index.
 * Called by magic2 --createIndex after the genome CW arrays are ready.
 *
 * Single entry point:
 *
 *   saTorchBuildIndex (cws, cwSizes, NN, meta, IDX)
 *
 * Preconditions (caller's responsibility):
 *   - cws[p] is a flat array of CW records for partition p
 *   - records within each cws[p] are sorted by CW.seed (ascending)
 *   - cwSizes[p] is the number of records in cws[p]  (may be 0)
 *   - NN == TORCH_NN == 16
 *   - meta->maxTargetRepeats, meta->step, meta->seedLength,
 *     meta->genomeLength, meta->nChroms are already filled by caller
 *   - IDX directory exists and is writable
 *
 * What this function does for each partition p in a single linear pass:
 *   - Walk cws[p] identifying runs of identical seed values
 *   - Drop any run longer than meta->maxTargetRepeats (repeat filter)
 *   - Emit one KEYS entry per surviving unique seed
 *   - Emit one OFFSETS entry per surviving unique seed
 *   - Stream surviving CW records directly to the .cw file
 *   - Write sentinel OFFSETS[K_p] = M_p
 *   - Write .keys and .offsets binary files
 *   - Accumulate K[p] and M[p] into meta
 *
 * After all partitions:
 *   - Fill meta->K_total, meta->M_total, meta->version, meta->nPartitions
 *   - Call saTorchMetaWrite  → IDX/torch_index.meta
 *   - Call saTorchMetaValidate as final sanity check
 *
 * Memory strategy: allocate KEYS and OFFSETS scratch at src_size
 * (upper bound: at most one entry per input record).  Single pass,
 * no realloc needed.  The .cw data is streamed directly from the
 * input array with no intermediate copy.
 *
 * Peak extra memory per partition: 2 × src_size × 4 bytes
 * (~500 MB for a human genome partition of ~56M records).
 * Only one partition is in flight at a time.
 *
 * Returns 1 on success, 0 on any error.
 *
 * Authors: Magic2 team / NLM-NIH
 * Created: May 2026
 */

#include "sa.h"
// #include "sa_torch_index.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/* ------------------------------------------------------------------ */
/* Internal: write a flat binary array to a file.                     */
/* Returns 1 on success, 0 on error.                                  */
/* ------------------------------------------------------------------ */
static int torchWriteBinary (const char *path,
                              const void *data,
                              size_t      n_bytes)
{
    FILE *f = fopen (path, "wb") ;
    if (!f) {
        fprintf (stderr, "[sa_torch_build] cannot create %s: %s\n",
                 path, strerror(errno)) ;
        return 0 ;
    }
    if (n_bytes > 0) {
        size_t written = fwrite (data, 1, n_bytes, f) ;
        if (written != n_bytes) {
            fprintf (stderr,
                "[sa_torch_build] short write on %s: "
                "wrote %zu of %zu bytes\n",
                path, written, n_bytes) ;
            fclose (f) ;
            return 0 ;
        }
    }
    fclose (f) ;
    return 1 ;
}

/* ------------------------------------------------------------------ */
/* saTorchBuildIndex                                                   */
/* ------------------------------------------------------------------ */

int saTorchBuildIndex (CW    **cws,
                       long int     *cwSizes,
                       int           NN,
                       TORCH_META   *meta,
                       const char   *IDX)
{
    if (NN != TORCH_NN) {
        fprintf (stderr,
            "[sa_torch_build] NN=%d but TORCH_NN=%d — aborting\n",
            NN, TORCH_NN) ;
        return 0 ;
    }

    /* Fill the fields we own */
    strncpy (meta->version, TORCH_INDEX_VERSION,
             sizeof(meta->version) - 1) ;
    meta->nPartitions = TORCH_NN ;
    meta->K_total     = 0 ;
    meta->M_total     = 0 ;

    int p ;
    for (p = 0 ; p < NN ; p++) { meta->K[p] = 0 ; meta->M[p] = 0 ; }

    /* ----------------------------------------------------------------
     * Process each partition in turn
     * ---------------------------------------------------------------- */
    for (p = 0 ; p < NN ; p++) {

        const CW *src      = cws[p] ;
        long int        src_size = cwSizes[p] ;
        int             maxRep   = meta->maxTargetRepeats ;

        /* ----------------------------------------------------------
         * Allocate KEYS and OFFSETS at the upper bound (src_size).
         * At most one unique seed per input record.
         * One extra slot in offsets for the sentinel.
         * ---------------------------------------------------------- */
        uint32_t *keys    = NULL ;
        uint32_t *offsets = NULL ;

        if (src_size > 0) {
            keys = malloc ((size_t)src_size * sizeof(uint32_t)) ;
            if (!keys) {
                fprintf (stderr,
                    "[sa_torch_build] malloc keys p=%d n=%ld: %s\n",
                    p, src_size, strerror(errno)) ;
                return 0 ;
            }
            offsets = malloc ((size_t)(src_size + 1) * sizeof(uint32_t)) ;
            if (!offsets) {
                fprintf (stderr,
                    "[sa_torch_build] malloc offsets p=%d n=%ld: %s\n",
                    p, src_size, strerror(errno)) ;
                free (keys) ;
                return 0 ;
            }
        }

        /* ----------------------------------------------------------
         * Open .cw file for streaming output
         * ---------------------------------------------------------- */
        char path[4096] ;
        snprintf (path, sizeof(path), "%s/" CW_FMT, IDX, p) ;
        FILE *cw_file = fopen (path, "wb") ;
        if (!cw_file) {
            fprintf (stderr,
                "[sa_torch_build] cannot create %s: %s\n",
                path, strerror(errno)) ;
            free (keys) ; free (offsets) ;
            return 0 ;
        }

        /* ----------------------------------------------------------
         * Single linear pass over src
         * ---------------------------------------------------------- */
        long int k       = 0 ;   /* surviving unique seeds (KEYS index)  */
        long int cw_out  = 0 ;   /* surviving CW records written to file */
        long int dropped = 0 ;   /* unique seeds dropped by repeat filter*/
        long int i       = 0 ;

        while (i < src_size) {
            uint32_t seed  = src[i].seed ;
            long int start = i ;

            /* measure the run of identical seeds */
            while (i < src_size && src[i].seed == seed) i++ ;
            long int run = i - start ;

            if (run > maxRep) {
                dropped++ ;
                continue ;          /* skip this seed entirely */
            }

            /* surviving seed: record key and offset */
            keys[k]    = seed ;
            offsets[k] = (uint32_t) cw_out ;
            k++ ;

            /* stream CW records directly from input to file — no copy */
            size_t n_written = fwrite (src + start, sizeof(CW),
                                       (size_t)run, cw_file) ;
            if ((long int)n_written != run) {
                fprintf (stderr,
                    "[sa_torch_build] short write .cw p=%d at i=%ld: %s\n",
                    p, i, strerror(errno)) ;
                fclose (cw_file) ;
                free (keys) ; free (offsets) ;
                return 0 ;
            }
            cw_out += run ;
        }

        /* sentinel: offsets[k] = total records written */
        if (k > 0)
            offsets[k] = (uint32_t) cw_out ;

        fclose (cw_file) ;

        fprintf (stderr,
            "[sa_torch_build] partition %02d: "
            "input=%ld  unique_kept=%ld  records=%ld  dropped=%ld\n",
            p, src_size, k, cw_out, dropped) ;

        /* ----------------------------------------------------------
         * Write .keys  (only the k surviving entries)
         * ---------------------------------------------------------- */
        snprintf (path, sizeof(path), "%s/" TORCH_KEYS_FMT, IDX, p) ;
        if (!torchWriteBinary (path, keys,
                               (size_t)k * sizeof(uint32_t))) {
            free (keys) ; free (offsets) ; return 0 ;
        }

        /* ----------------------------------------------------------
         * Write .offsets  (k+1 entries: k offsets + sentinel)
         * ---------------------------------------------------------- */
        snprintf (path, sizeof(path), "%s/" TORCH_OFFSETS_FMT, IDX, p) ;
        if (!torchWriteBinary (path, offsets,
                               (size_t)(k + 1) * sizeof(uint32_t))) {
            free (keys) ; free (offsets) ; return 0 ;
        }

        free (keys) ;
        free (offsets) ;

        /* accumulate into meta */
        meta->K[p]     = k ;
        meta->M[p]     = cw_out ;
        meta->K_total += k ;
        meta->M_total += cw_out ;
    }

    /* ----------------------------------------------------------------
     * Write metadata and validate
     * ---------------------------------------------------------------- */
    if (!saTorchMetaWrite (IDX, meta)) {
        fprintf (stderr,
            "[sa_torch_build] failed to write %s/%s\n",
            IDX, TORCH_META_FILENAME) ;
        return 0 ;
    }

    int bad = saTorchMetaValidate (meta) ;
    if (bad > 0) {
        fprintf (stderr,
            "[sa_torch_build] %d inconsistencies in metadata — "
            "index may be corrupt\n", bad) ;
        return 0 ;
    }

    fprintf (stderr,
        "[sa_torch_build] index complete: "
        "K_total=%ld  M_total=%ld  dir=%s\n",
        meta->K_total, meta->M_total, IDX) ;

    return 1 ;
}

