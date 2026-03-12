# Magic2

## Getting started

### Create a target index
```
sortalign --createIndex IDX -t target.fasta
```

### Map public SRA reads on this target
```
sortalign -x IDX -i SRR1234,SRR5678 -o results

### Map private reads on this target
```
sortalign -x IDX -i fasta.gz -o results

### Map paired-ends fastq files
```
sortalign -x IDX -i forward.fastq.gz+reverse.fastq.gz -o results
# SAM output
sortalign -x IDX -i forward.fastq.gz+reverse.fastq.gz -o results --sam
# Multiple runs
sortalign -x IDX -i forward.fastq.gz+reverse.fastq.gz,f1.fastq.gz+f2.fastq.gz -o results
# Report target coverage
sortalign -x IDX -i forward.fastq.gz+reverse.fastq.gz,f1.fastq.gz+f2.fastq.gz -o results --wiggle
# If you do not need alignments
sortalign -x IDX -i forward.fastq.gz+reverse.fastq.gz,f1.fastq.gz+f2.fastq.gz -o results --wiggle --noalign
```

## User's guide

Sortalign maps deep-sequencing runs and reports alignments as well as coverage plots and introns. In the first step, sortalign analyses the target genome (specified by -t or --target), or a set or targets and their annotations (specified by a configuration file -T) and creates an index. This is needed only once per target and takes about 10 minutes for the human genome. When the index is ready, sortalign can map private data (in FASTA or FASTQ format). or automatically download and map data from NCBI Sequence Read Archive. 

When mapping to a human genome, the program requires about 30GB of memory. Depending on the hardware, the it will launch n threads (between 6 and 32 on current computers) and can process about 1 Gigabase of human RNA-seq reads per minute.

In addition to the alignmnets, the program exports, at nearly no additional time cost, the coverage plots, the intron discovery and support, the gene expression indices, the location of the poly-A addition sites, the candidate start and ends of the transcripts, the details on the mismatches (substitution and indels with explict genoe coordinates) and many detailled metrics.

### Index generation

Below are additional options for generating a target index with the default values given on the command line.
```
sortalign --createIndex IDX -t target.fasta --seedLength 18 --maxTargetRepeats 81
```

### Mapping private sequence files
`sortalign` takes input files as a comma-separated list of file names or an NCBI SRA accession. File names for paired reads are separated by `+`. Use `-` for standard input. Files with names ending with ".gz" are automatically decompressed. Adding `--gzi` forces decompression of input files. Sortalign guesses file format from a file name, but one can also specify the format of input files with flags: `--fasta`, `--fastq`, `--fastc`, `--raw`, `-sra`.

```
sortalign -x IDX -i f1.fastq.gz+f2.fastq.gz -o results --sam
```

### NCBI SRA
`sortalign` can download and map reads from an NCBI Sequence Read Archive (SRA). Just specify SRA accession in `-i` option.
```
sortalign -x IDX -i SRR24511885 -o results --sam  [--sraCaching] [--maxGB <int>]
```
The --maxGB option will just download and analyze the first n gigabases of the SRA dataset.

The `--sraCaching` option will cache the reads in a local directory called ./SRA. This is only useful if you need the reads several times, because sortalign starts aligning immediatly as soon as some data is received, so on a slow network, all the data will be aligned as soon as it has been fully transmitted.

You can also use sortalign just to download SRA runs in your local cache, without analyzing them using the command
sortalign --sraDownload SRR24511885 [--maxGB <int>]


### Output

The `-o` option specifies a path to the output directory. By default `sortalign` reports alignments in its internal format. To request the SAM report, use `--sam` option. `--gzo` will gzip output files.

You can also generate target coverage wiggles in UCSC format with `--wiggles`.
```
sortalign -x IDX -i SRR24511885 -o results --wiggle
```
