/* ---------------------------------------------------------------
 * magic2  --help  output
 * Include this block inside your existing help dispatch, e.g.:
 *   if (strcmp(arg, "-h") == 0 || strcmp(arg, "--help") == 0)
 *     { magic2_help(); exit(0); }
 * ---------------------------------------------------------------*/
#include "stdio.h"

void saHelp (void)
{
  fprintf (stderr,
    "magic2  --  Deep RNA-sequencing aligner\n"
    "National Library of Medicine, NIH  |  https://github.com/NLM-DIR/Magic2\n"
    "\n"
    "Usage:\n"
    "  magic2 --createIndex <IDX> { -t <genome.fa> | -T <tConfig> } [options]\n"
    "  magic2 -x <IDX> { -i <reads> | --srr <SRR,...> } -o <outDir> [options]\n"
    "  magic2 --sraDownload <SRR,...> [--fastq] [--split_pairs] [--maxGB <INT>]\n"
    "  magic2 --help | --devhelp | --version\n"
    "\n"
  ) ;

  fprintf (stderr,
    "Quick examples:\n"
    "  # Build index once, then align\n"
    "  magic2 --createIndex IDX -t genome.fasta.gz\n"
    "  magic2 -x IDX -i reads.fastq.gz -o results/ --bam --full\n"
    "\n"
    "  # Annotated index (recommended for RNA-seq)\n"
    "  magic2 --createIndex IDX -T human.tConfig\n"
    "  magic2 -x IDX -i R1.fastq.gz+R2.fastq.gz -o results/ --bam --full\n"
    "\n"
    "  # Stream directly from SRA, no local storage of raw reads\n"
    "  magic2 -x IDX --srr SRR24511885 -o results/ --bam --full\n"
    "\n"
    "  # Quick QC survey across several SRA runs (no alignment files written)\n"
    "  magic2 -x IDX --srr SRR001,SRR002,SRR003 --no_ali -o qc/\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "INDEX CONSTRUCTION  (run once per reference genome)\n"
    "------------------------------------------------------------\n"
    "  --createIndex <DIR>      Directory to write the index into (required)\n"
    "  -t <FILE>                Single FASTA genome file (raw genome, no annotation)\n"
    "  -T <FILE>                Target configuration file (recommended; see below)\n"
    "\n"
    "  Target configuration file (-T):  two-column tab-separated text, one entry per line\n"
    "    G  genome.fasta.gz          Genomic sequence (required)\n"
    "    M  mito.fasta.gz            Mitochondrial sequence\n"
    "    C  chloro.fasta.gz          Chloroplast sequence\n"
    "    R  rRNA.fasta.gz            Ribosomal RNA (strongly recommended)\n"
    "    I  annotation.gtf.gz        Gene annotation: GTF, GFF/GFF3, or .introns\n"
    "    E  ERCC.fasta               Spike-in / control sequences\n"
    "    A  adaptors.fasta           Adaptor sequences (auto-detected if omitted)\n"
    "    B  bacteria.fasta           Bacterial contamination sequences\n"
    "    V  virus.fasta              Viral contamination sequences\n"
    "  Note: chromosome names in G and I files must match exactly.\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "INPUT\n"
    "------------------------------------------------------------\n"
    "  -x, --index <DIR>        Index directory created by --createIndex (required)\n"
    "  -i <FILE[s]>             Input reads.  Format detected from file extension.\n"
    "                             Single-end:    -i reads.fastq.gz\n"
    "                             Paired-end:    -i R1.fastq.gz+R2.fastq.gz\n"
    "                             Interleaved:   -i sample.sample_12.fastq.gz\n"
    "                             Multiple runs: -i run1.fa.gz,run2.fa.gz\n"
    "                             Pipe (fastq):  zcat *.fastq.gz | magic2 -x IDX -i - --fastq\n"
    "  -I <FILE>                Tab-separated run manifest (filename, run name, descriptors)\n"
    "  --srr <SRR,...>          Comma-separated SRA run identifiers; streamed from NCBI\n"
    "                             or read from local ./SRA/ cache if present\n"
    "  --sraDownload <SRR,...>  Download SRA runs to ./SRA/ then exit\n"
    "    --fastq                  Download in FASTQ format [default: FASTA]\n"
    "    --split_pairs            Write pairs as _R1 / _R2 files [default: interleaved]\n"
    "    --maxGB <INT>            Stop after downloading INT gigabases [default: unlimited]\n"
    "  --sraCaching             Cache SRA reads locally while aligning on the fly\n"
    "  -r, --run <STR>          Label for this run (used to name output subdirectories)\n"
    "\n"
    "  Explicit format overrides (use when format cannot be inferred from filename):\n"
    "    --fasta | --fastq | --fastc | --raw | --interleaved\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "OUTPUT\n"
    "------------------------------------------------------------\n"
    "  -o <DIR>                 Output directory; created if it does not exist (required)\n"
    "  --bam                    Export alignments in BAM format [default]\n"
    "  --hits                   Export alignments in Magic .hits format\n"
    "                             (explicit mismatch, intron, and overhang detail)\n"
    "  --tabular                Export alignments in tabular (Magic-BLAST) format\n"
    "  --no_ali                 Suppress alignment output (QC metrics and counts only)\n"
    "  --full                   Enable all outputs: --genes --wiggles --wiggle_ends\n"
    "  --genes                  Export gene expression count tables (TSF format)\n"
    "  --wiggles                Export coverage plots (UCSC wiggle / BF format)\n"
    "  --wiggle_step <INT>      Coverage plot resolution in bases [default: 10]\n"
    "  --wiggle_ends            Export transcription-start / end-site coverage tracks\n"
    "  --introns                Export intron support table (TSF format) [default: on]\n"
    "  --gzo                    Compress output files with gzip\n"
    "  --quality_factors        Echo base-quality scores into BAM/SAM output\n"
    "\n"
    "  Merging per-run TSF tables after multi-run analysis:\n"
    "    magic2 --mergeTSF run1/genes.tsf run2/genes.tsf -o all_genes.txt\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "ALIGNMENT BEHAVIOUR\n"
    "------------------------------------------------------------\n"
    "  --RNA                    Force RNA alignment mode (intron-aware, poly-A detection)\n"
    "  --DNA                    Force DNA alignment mode (continuous alignments only)\n"
    "                             [default: auto-detected from data]\n"
    "  --no_splice              Reject spliced alignments; require continuous alignment\n"
    "  --maxIntron <INT>        Maximum intron size [default: 1000000]\n"
    "  --strandedness <MODE>    Override library strandedness: auto | forward | reverse\n"
    "                             | unstranded  [default: auto-detected from GT-AG ratio]\n"
    "\n"
    "  Scoring (change only if the defaults produce poor results for your data):\n"
    "  --minAli <INT>           Minimum alignment length in bases [default: 30]\n"
    "  --minScore <INT>         Minimum alignment score [default: 30]\n"
    "                             Score = aligned_bases - errCost * error_count\n"
    "  --errCost <INT>          Penalty per substitution or short indel (<= 3 b)\n"
    "                             [default: 8]\n"
    "  --errMax <INT>           Hard cap on mismatches per alignment [default: none]\n"
    "  --errRateMax <INT>       Maximum mismatch percentage per alignment [default: 10]\n"
    "\n"
    "  Adaptor handling:\n"
    "  --adaptor1 <SEQ>         Read-1 3' exit adaptor sequence (auto-detected if omitted)\n"
    "  --adaptor2 <SEQ>         Read-2 5' exit adaptor sequence (auto-detected if omitted)\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "READ FILTERING\n"
    "------------------------------------------------------------\n"
    "  --minReadLength <INT>    Drop reads shorter than INT bases\n"
    "                             [default: 20, or max-read-length for short-RNA data]\n"
    "  --minEntropy <INT>       Drop low-complexity reads below this entropy threshold\n"
    "                             [default: 20, or max-read-length/2 for short-RNA data]\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "PERFORMANCE\n"
    "------------------------------------------------------------\n"
    "  --threads <INT>          Number of CPU threads to use [default: all available]\n"
    "                             magic2 parallelises across pipeline layers, not just\n"
    "                             individual steps; this is the primary tuning knob.\n"
    "                             Equivalent programs: -p in HISAT2, --runThreadN in STAR.\n"
    "\n"
    "  Note: on shared clusters, avoid requesting fewer than 16 threads; the internal\n"
    "  channel architecture requires a minimum concurrency to sustain throughput.\n"
    "  If no progress is logged within the first minute, try increasing --threads.\n"
    "\n"
    "  Run magic2 --devhelp for additional low-level pipeline parameters\n"
    "  (--nAgents, --nBlocks, --bMax, --NN, --max_threads) intended for developers.\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "MISCELLANEOUS\n"
    "------------------------------------------------------------\n"
    "  --do_not_align           Validate the index directory without aligning\n"
    "  --verbose                Emit detailed diagnostic messages\n"
    "  --version                Print version string and exit\n"
    "  --help, -h               Print this help message and exit\n"
    "  --devhelp                Print developer / pipeline-tuning parameters and exit\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "CITATION\n"
    "------------------------------------------------------------\n"
    "  If you use magic2 in published work, please cite:\n"
    "  Thierry-Mieg J, Thierry-Mieg D, Boratyn G.\n"
    "  Magic2: deep RNA-sequencing aligner.  https://github.com/NLM-DIR/Magic2\n"
    "\n"
    "  Related reference:\n"
    "  Thierry-Mieg D & Thierry-Mieg J.  AceView: a comprehensive cDNA-supported\n"
    "  gene and transcripts annotation.  Genome Biology 7(Suppl 1):S12, 2006.\n"
    "\n"
  ) ;
}


/* ---------------------------------------------------------------
 * magic2  --devhelp  output
 * Low-level pipeline parameters; not intended for end users.
 * ---------------------------------------------------------------*/

static void magic2_devhelp (void)
{
  fprintf (stderr,
    "magic2  --devhelp  (developer / pipeline-tuning parameters)\n"
    "These parameters control internal concurrency and index structure.\n"
    "The defaults are chosen to work well on hardware from a laptop to a\n"
    "large server.  Change them only if you understand the pipeline architecture.\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "PIPELINE CONCURRENCY\n"
    "------------------------------------------------------------\n"
    "  magic2 is parallelised using Go-style channels: all pipeline layers run\n"
    "  simultaneously and self-balance their workload.  Each agent and channel\n"
    "  runs in its own virtual thread.  The parameters below control the number\n"
    "  of agents and in-flight data blocks, not a simple thread pool.\n"
    "\n"
    "  --nAgents, --nA <INT>    Agents per pipeline layer [default: 10]\n"
    "  --nBlocks, --nB <INT>    Simultaneous data blocks in the pipeline [default: 20]\n"
    "  --bMax <INT>             Maximum megabases per data block [default: 20; range 1-1024]\n"
    "  --max_threads <INT>      Hard cap on UNIX threads [default: 128]\n"
    "                             Do not set below 64; the pipeline will stall.\n"
    "\n"
    "  Recommended configurations:\n"
    "    Small test:            --nBlocks 1\n"
    "    Default (< 32 GB RAM): --nAgents 10 --nBlocks 20 --max_threads 128\n"
    "    Large server:          --nAgents 20 --nBlocks 30 --max_threads 512\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "INDEX STRUCTURE\n"
    "------------------------------------------------------------\n"
    "  --NN <INT>               Split seed files into NN parts [default: 16]\n"
    "                             Allowed values: 1, 2, 4, 8, 16, 32, 64\n"
    "                             seedLength > 16 requires NN >= 4^(seedLength-16)\n"
    "  --noJump                 Do not insert jumpers in the index\n"
    "                             Reserved for future GPU version\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "INDEX CONSTRUCTION\n"
    "------------------------------------------------------------\n"
    "  --seedLength <INT>       Seed length in bases [default: 16 for genomes <1 Gb,\n"
    "                             18 for larger genomes (human / mouse)]\n"
    "  --step <INT>             During --createIndex: seed step [default: 1 for targets\n"
    "                             < 1 Mb, 6 for larger targets]\n"
    "                             Larger values (8-10) reduce RAM and run faster\n"
    "                             but reduce sensitivity\n"
    "	                          Smaller values (1-4) increase sensitivity, RAM and CPU cost\n"
    "  --maxTargetRepeats <INT> Do not index target words repeated more than this\n"
    "                             many times [default: 81]\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "ALIGNMENTS\n"
    "------------------------------------------------------------\n"
    "  --step <INT>             Step between query seeds [default: targetStep/2]\n"
    "                             --step 1 increases sensitivity at higher CPU cost\n"
    "  --ignoreIntronSeeds      Ignore known introns from the -T annotation during\n"
    "                             seed lookup [default: false]\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "INPUT / COMPRESSION\n"
    "------------------------------------------------------------\n"
    "  --gzi                    Force decompression of input files (useful when piping)\n"
    "                             Files named *.gz are decompressed automatically\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "DEBUGGING  (for C developers; requires understanding of the re-spawn model)\n"
    "------------------------------------------------------------\n"
    "  magic2 normally re-spawns itself via numactl to bind to idle cores before\n"
    "  declaring its agent/channel topology.  This means that attaching gdb to the\n"
    "  initial process and setting a breakpoint will not work as expected: the\n"
    "  breakpoint is set in the parent but execution continues in the child.\n"
    "\n"
    "  Two flags suppress re-spawning for debugging purposes:\n"
    "\n"
    "  --debug                  Suppress re-spawning and additionally forces --nAgents 1\n"
    "                             --nBlocks 1, giving a fully serial execution path\n"
    "                             suitable for step-by-step debugging.\n"
    "  --numactl                Suppress re-spawning.  The process runs in place,\n"
    "                             so gdb breakpoints and valgrind instrumentation\n"
    "                             work normally.  Parallelism is otherwise unchanged.\n"
    "                             Note: this flag is added automatically by the\n"
    "                             re-spawn call itself to prevent recursion; do not\n"
    "                             be surprised to see it in the child process argv.\n"
    "\n"
    "  Typical gdb sessions:\n"
    "    # Breakpoints work, fully serial (easiest to step through):\n"
    "    gdb --args magic2 --debug -x IDX -i test.fa -o /tmp/dbg/\n"
    "\n"
    "    # Breakpoints work, full parallelism:\n"
    "    gdb --args magic2 --numactl -x IDX -i test.fa -o /tmp/dbg/\n"
    "\n"
  ) ;

  fprintf (stderr,
    "  Run magic2 --help for the standard user-facing options.\n"
    "\n"
  ) ;
}
