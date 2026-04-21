/* ---------------------------------------------------------------
 * file sa.help.c
 * invoked by the command
 *      magic2  --help  [--devhelp]
 * ---------------------------------------------------------------*/
#include "sa.h"

void saHelp (void)
{
  fprintf (stderr,
    "------------------------------------------------------------\n"
    "MAGIC2 on line HELP\n"
    "------------------------------------------------------------\n"
    "magic2  --  A Deep RNA-sequencing analyzer\n"
    "National Library of Medicine, NIH  |  https://github.com/NLM-DIR/Magic2\n"
    "\n"
    "Usage:\n"
    "  magic2 --createIndex <IDX> { -t <genome.fa> | -T <tConfig> } [options]\n"
    "  magic2 -x <IDX> { -i <reads> | -I <rConfig> | --srr <SRR,...> } -o <outDir> [options]\n"
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
    "  # Align multiple samples described in a run manifest\n"
    "  magic2 -x IDX -I my_project.rConfig -o results/ --bam --full\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "INDEX CONFIGURATION  (needed once per reference genome)\n"
    "------------------------------------------------------------\n"
    "  --createIndex <DIR>      Directory to write the index into (required)\n"
    "  -t <FILE>                Single FASTA genome file (raw genome, no annotation)\n"
    "  -T <FILE>                Target configuration file (recommended; see below)\n"
    "  -t and -T are mutually exclusive.\n"
    "  -t accepts a single file; -T accepts a manifest describing multiple targets.\n"
    "\n"
    "  Target configuration file (-T):  two-column tab-separated text, one entry per line\n"
    "    G  genome.fasta.gz          Genomic sequence (required)\n"
    "    M  mito.fasta.gz            Mitochondrial sequence\n"
    "    C  chloro.fasta.gz          Chloroplast sequence\n"
    "    R  rRNA.fasta.gz            Ribosomal RNA (strongly recommended)\n"
    "    I  annotation.gtf.gz        Gene annotation: GTF, GFF/GFF3, or .introns\n"
    "    E  ERCC.fasta.gz            Spike-in / control sequences\n"
    "    B  bacteria.fasta.gz        Bacterial contamination sequences\n"
    "    V  virus.fasta.gz           Viral contamination sequences\n"
    "  Notes:\n"
    "    Chromosome names in G and I files must match exactly.\n"
    "    Multiple G lines are allowed if chromosomes are in separate files.\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "INDEX CONSTRUCTION\n"
    "------------------------------------------------------------\n"
    "  --seedLength <INT>       Seed length in bases [default: 16 for genomes <1 Gb,\n"
    "                             18 for larger genomes (human / mouse)]\n"
    "  --step <INT>             Seed step during index construction\n"
    "                             [default: 1 for targets < 1 Mb, 6 for larger genomes]\n"
    "                             Step is the maximum distance between consecutive seeds;\n"
    "                             average distance is step/2.\n"
    "                             Larger values (8-10) reduce RAM and index faster\n"
    "                             but reduce sensitivity.\n"
    "                             Smaller values (1-4) increase sensitivity at higher\n"
    "                             RAM and CPU cost.\n"
    "                             The alignment seed step is set automatically to step/2.\n"
    "  --maxTargetRepeats <INT> Do not index target words repeated more than this\n"
    "                             many times [default: 81]\n"
    "                             This effectively masks highly repeated regions of the genome.\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "SEQUENCES TO BE ANALYZED\n"
    "------------------------------------------------------------\n"
    "  -x, --index <DIR>        Index directory created by --createIndex (required)\n"
    "\n"
    "  -i and -I are mutually exclusive:\n"
    "    -i accepts a single file or paired-end pair; options come from the command line.\n"
    "    -I accepts a run manifest file; options are specified per sample within the file.\n"
    "\n"
    "  -i <FILE or R1+R2>       A single input file, or a paired-end pair joined with +.\n"
    "                             Format is detected from the file extension.\n"
    "                             Single-end:   -i reads.fastq.gz\n"
    "                             Paired-end:   -i R1.fastq.gz+R2.fastq.gz\n"
    "                             Interleaved:  -i sample.sample_12.fastq.gz\n"
    "                               (interleaved: consecutive records ending in /1 and /2;\n"
    "                                4 lines per pair in FASTA, 8 lines per pair in FASTQ)\n"
    "                             Pipe:         zcat *.fastq.gz | magic2 -x IDX -i - --fastq\n"
    "                             For multiple files or samples, use -I instead.\n"
    "\n"
    "  -I <FILE>                Run manifest file (rConfig): tab-separated, 3 columns.\n"
    "                             Column 1: filename (FASTA or FASTQ, plain or .gz).\n"
    "                                       For paired-end: R1.fastq.gz+R2.fastq.gz\n"
    "                             Column 2: run name.  Use only letters, digits, and\n"
    "                                       underscores: the run name becomes part of\n"
    "                                       output file and directory names.\n"
    "                                       Multiple files sharing a run name are merged\n"
    "                                       into a single output (BAM, counts, wiggles).\n"
    "                             Column 3: comma-separated properties (optional).\n"
    "                                       RNA or DNA          [default: auto-detected]\n"
    "                                       fasta, fastq, fastc [default: from suffix]\n"
    "                                       interleaved         paired reads in one file\n"
    "                                       Adaptor1=<seq>      read-1 exit adaptor\n"
    "                                       Adaptor2=<seq>      read-2 exit adaptor\n"
    "                             Example rConfig file:\n"
    "                               lib1_R1.fastq.gz+lib1_R2.fastq.gz  sample_A  RNA\n"
    "                               lib2.fasta.gz                      sample_A  RNA\n"
    "                               lib3_R1.fastq.gz+lib3_R2.fastq.gz  sample_B  RNA,Adaptor1=AGATCG\n"
    "                             (sample_A merges lib1 and lib2 into one output)\n"
    "\n"
    "  --srr <SRR,...>          Comma-separated SRA run identifiers; streamed from NCBI\n"
    "                             or read from local ./SRA/ cache if present.\n"
    "                             The SRR identifier is used as the run name.\n"
    "  --sraDownload <SRR,...>  Download SRA runs to ./SRA/ then exit\n"
    "    --fastq                  Download in FASTQ format [default: FASTA]\n"
    "    --split_pairs            Write pairs as _R1 / _R2 files [default: interleaved]\n"
    "    --maxGB <INT>            Stop after downloading INT gigabases [default: unlimited]\n"
    "  --sraCaching             Cache SRA reads locally while aligning on the fly\n"
    "  -r, --run <STR>          Run label for -i input (used in output file names)\n"
    "\n"
    "  Explicit format overrides for -i (use when format cannot be inferred from filename):\n"
    "    --fasta | --fastq | --fastc | --raw | --interleaved\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "OUTPUT\n"
    "------------------------------------------------------------\n"
    "  -o <DIR>                 Output directory (required).\n"
    "                             The full directory path is created if it does not exist\n"
    "                             (e.g. -o results/human/run1/ creates all intermediate\n"
    "                             directories). Relative paths are anchored to the current\n"
    "                             working directory; absolute paths (starting with /) are\n"
    "                             used as-is.\n"
    "  --hits                   Export alignments in Magic .hits format [default]\n"
    "                             (explicit mismatches, introns, and overhang details)\n"
    "  --tabular                Export alignments in tabular (Magic-BLAST) format\n"
    "  --bam                    Export alignments in BAM format\n"
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
    "ALIGNMENT FILTERS AND SCORING\n"
    "------------------------------------------------------------\n"
    "  Change these only if the defaults produce poor results for your data.\n"
    "  --minReadLength <INT>    Drop reads shorter than INT bases\n"
    "                             [default: 20, or max-read-length for short-RNA data]\n"
    "  --minEntropy <INT>       Drop low-complexity reads below this entropy threshold\n"
    "                             [default: 20, or max-read-length/2 for short-RNA data]\n"
    "  --minAli <INT>           Minimum alignment length in bases [default: 30]\n"
    "  --minScore <INT>         Minimum alignment score [default: 30]\n"
    "                             Score = aligned_bases - errCost * error_count\n"
    "  --errCost <INT>          Penalty per substitution or short indel (<= 3 b)\n"
    "                             [default: 4]\n"
    "  --errMax <INT>           Hard cap on mismatches per alignment [default: none,\n"
    "                             except 1 error maximum for reads shorter than 40 bases]\n"
    "  --errRateMax <INT>       Maximum mismatch percentage per alignment [default: 12]\n"
    "  --maxIntron <INT>        Maximum intron size [default: 1 Megabase]\n"
    "  --no_splice              Reject spliced alignments; require continuous alignment\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "OPTIONAL ALIGNMENT STRATEGY\n"
    "------------------------------------------------------------\n"
    "  These parameters are evaluated automatically, but can be forced if needed.\n"
    "  --RNA                    Force RNA alignment mode (intron-aware, poly-A detection)\n"
    "  --DNA                    Force DNA alignment mode (continuous alignments only)\n"
    "                             [default: auto-detected from data]\n"
    "  --strand <MODE>          Force library strandedness: forward | reverse | unstranded\n"
    "                             [default: auto-detected from GT-AG/CT-AC ratio]\n"
    "  --adaptor1 <SEQ>         Read-1 exit adaptor sequence (auto-detected if omitted)\n"
    "  --adaptor2 <SEQ>         Read-2 exit adaptor sequence (auto-detected if omitted)\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "READ ALIGNMENT TUNING\n"
    "------------------------------------------------------------\n"
    "  --step <INT>             Seed step during alignment [default: index step / 2]\n"
    "                             --step 1 maximises sensitivity at higher CPU cost\n"
    "  --ignoreIntronSeeds      Do not use known splice junctions from the -T annotation\n"
    "                             to guide alignment. Useful when evaluating novel organisms\n"
    "                             or benchmarking de-novo intron discovery. [default: off]\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "PARALLELISM\n"
    "------------------------------------------------------------\n"
    "  magic2 manages its own parallelism automatically.\n"
    "  It interrogates the hardware, binds to the least busy NUMA node,\n"
    "  and self-tunes the number of pipeline agents and data blocks based\n"
    "  on available RAM.  No thread count needs to be specified.\n"
    "  If no progress is logged within the first 60 seconds, run\n"
    "  magic2 --devhelp for low-level tuning options.\n"
    "\n"
    "  To process multiple runs in parallel, simply launch independent\n"
    "  magic2 processes: the index is memory-mapped and shared at no\n"
    "  extra RAM cost.\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "MISCELLANEOUS\n"
    "------------------------------------------------------------\n"
    "  --version                Print version string and exit\n"
    "  --gzi                    Force decompression of input files (useful when piping)\n"
    "                             Files named *.gz are decompressed automatically\n"
    "  --do_not_align           Validate the index directory without aligning\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "DOCUMENTATION\n"
    "------------------------------------------------------------\n"
    "  --help, -h               Print this help message and exit\n"
    "  --devhelp                Print developer / pipeline-tuning parameters and exit\n"
    "  Full documentation:      wsa_doc/ in the source tree\n"
    "                           https://github.com/NLM-DIR/Magic2\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "CITATION\n"
    "------------------------------------------------------------\n"
    "  If you use magic2 in published work, please cite:\n"
    "  Thierry-Mieg J, Thierry-Mieg D, Boratyn G.\n"
    "  Magic2: deep RNA-sequencing analyzer.  https://github.com/NLM-DIR/Magic2\n"
    "\n"
    "  Related reference:\n"
    "  Thierry-Mieg D & Thierry-Mieg J.  AceView: a comprehensive cDNA-supported\n"
    "  gene and transcripts annotation.  Genome Biology 7(Suppl 1):S12, 2006.\n"
    "\n"
  ) ;
  exit (0) ;
}


/* ---------------------------------------------------------------
 * magic2  --devhelp  output
 * Low-level pipeline parameters; not intended for end users.
 * ---------------------------------------------------------------*/

void saDevHelp (void)
{
  fprintf (stderr,
    "------------------------------------------------------------\n"
    "MAGIC2 HELP FOR DEVELOPERS AND ADVANCED USERS\n"
    "------------------------------------------------------------\n"
    "Try magic2 --help for the standard user information.\n"
    "\n"
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
    "  The program interrogates the hardware to identify the least busy NUMA node\n"
    "  and counts the number nCPU of CPUs attached to that node.\n"
    "\n"
    "  --nAgents, --nA <INT>    Agents per pipeline layer [default: 3/2 * nCPU]\n"
    "  --nBlocks, --nB <INT>    Simultaneous data blocks in the pipeline [default: 2 * nCPU]\n"
    "  --bMax <INT>             Maximum megabases per data block [default: 3; range 1-100]\n"
    "  --max_threads <INT>      Hard cap on UNIX threads [default: derived from hardware]\n"
    "                             Do not set below 64; the pipeline will stall.\n"
    "\n"
    "  Recommended configuration: do NOT alter these parameters.\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "READ ALIGNMENT TUNING\n"
    "------------------------------------------------------------\n"
    "  --rStep <INT>            Seed step during read alignment [default: index step / 2]\n"
    "                             Reducing to 1 maximises sensitivity at higher CPU cost.\n"
    "                             Normally set automatically; override only if benchmarking.\n"
    "\n"
  ) ;

  fprintf (stderr,
    "------------------------------------------------------------\n"
    "INDEX STRUCTURE\n"
    "------------------------------------------------------------\n"
    "  --NN <INT>               Split seed files into NN parts [default: 16]\n"
    "                             Allowed values: 1, 2, 4, 8, 16, 32, 64\n"
    "                             seedLength > 16 requires NN >= 4^(seedLength-16);\n"
    "                             this minimum is always enforced automatically.\n"
    "\n"
    "  Recommended configuration: do NOT alter this parameter.\n"
    "\n"
  ) ;

  fprintf (stderr,
	   "------------------------------------------------------------\n"
	   "DEBUGGING  (for C developers)\n"
	   "------------------------------------------------------------\n"
	   "  On machines with multiple NUMA nodes (multiple physical processor\n"
	   "  packages), magic2 re-spawns itself via numactl to bind to the least\n"
	   "  busy node, ensuring memory locality before declaring its agent/channel\n"
	   "  topology.  Re-spawning is skipped on single-node machines, non-NUMA\n"
	   "  kernels, and containers where NUMA topology is not exposed.\n"
	   "  When re-spawning does occur, attaching gdb to the initial process and\n"
	   "  setting a breakpoint will not work as expected: the breakpoint is set\n"
	   "  in the parent but execution continues in the child.\n"
	   "\n"
    "  Two flags suppress re-spawning for debugging purposes:\n"
    "\n"
    "  --debug                  Suppress re-spawning and force --nBlocks 1 --nAgents 1,\n"
    "                             giving a fully serial execution path suitable for\n"
    "                             step-by-step debugging. Verbose diagnostic messages\n"
    "                             are enabled automatically in --debug mode.\n"
    "  --numactl                Suppress re-spawning only. The process runs in place\n"
    "                             so gdb breakpoints and valgrind instrumentation work\n"
    "                             normally. Parallelism is otherwise unchanged.\n"
    "                             Note: this flag is also added automatically by the\n"
    "                             re-spawn call itself to prevent recursion; do not\n"
    "                             be surprised to see it in the child process argv.\n"
    "\n"
    "  Typical gdb sessions:\n"
    "    # Fully serial, easiest to step through:\n"
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
  exit (0) ;
}
