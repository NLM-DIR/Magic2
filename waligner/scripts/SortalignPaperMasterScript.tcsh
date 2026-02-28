#!/bin/tcsh
setenv SV v4.aug13
setenv SV v5.aug17
setenv SV v5.aug17
setenv SV v6.aug18  # they all crash signal 9, memory leak
setenv SV v8.aug19
setenv SV v9.aug21
setenv SV v10.aug23
setenv SV v11.aug24
setenv SV v12.aug27
setenv SV v13.aug30
setenv SV v14.sep2
setenv SV v15.sep4      # before clustering
setenv SV v16.sep6      # 16mers, clustering on seed counts 1/4 mult^2, maxTargetRepeats 4
setenv SV v17.sep6      # 16mers, clustering on seed counts 1/4 mult^2
setenv SV v18.sep7      # 18mers, maxTargetRepeats 12
setenv SV v19.sep7      # 18mers, maxTargetRepeats 31
setenv SV v20.sep8
setenv SV v21.sep10     # new clustering counting seeds in read coordinates
setenv SV v22.sep12     # 2.5 min words
setenv SV v23.sep15     # min hash of all words up to step
setenv SV v24.sep22    # same code as 23, but new methods G5R5  G3R3 G3R1 with min hashing    
setenv SV v25.sep22        # fixed introns boundaries, prefer probe with less repeats, maxTargetRepeats 81, seedLength 18
setenv SV v26.sep23        # BAD on all metrics same as V25, but maxTargetRepeats 31, seedLength 18
setenv SV v27.sep25        # same as 25 maxTargetRepeats 81, seedLength 16
setenv SV v28.sep27        # maxTargetRepeats 31, seedLength 18, * for exon words
setenv SV v29.sep29        # maxTargetRepeats 81, seedLength 18, MAXJUMP 3
setenv SVlast v29.sep29    # maxTargetRepeats 81, seedLength 18, fixed topology bug, MAXJUMP 3
setenv SV v30.sep30        # maxTargetRepeats 31, seedLength 18,  MAXJUMP 3
setenv SVlast v30.sep30    # maxTargetRepeats 31, seedLength 18, fixed topology bug, MAXJUMP 3
setenv SV v31.oct6        # maxTargetRepeats 31, seedLength 18,  MAXJUMP 8, fixed snps and topology,  minScore 20
setenv SVlast v31.oct6
setenv SV v32.oct9        # maxTargetRepeats 31, seedLength 18,  MAXJUMP 8, fixed snps and topology,  minScore 20 # sam2bam works on Roche, iRefSeq38, ChipSeq1/2
setenv SVlast v32.oct9
setenv SV v33.oct10        # maxTargetRepeats 81, seedLength 18,  MAXJUMP 8,  same code as v32
setenv SVlast v33.oct10
setenv SV v34.oct10        # maxTargetRepeats 81, seedLength 18,  MAXJUMP 3,  same code as v32
setenv SVlast v34.oct10
setenv SV v35.oct13        # maxTargetRepeats 81, seedLength 18,  MAXJUMP 3,  align on image of transcript
setenv SVlast v35.oct13
setenv SV v36.oct16        # maxTargetRepeats 81, seedLength 18,  MAXJUMP 3, maxScore 0, restored step=1 on < 60bp
setenv SVlast v36.oct16
setenv SV v37.oct16        # maxTargetRepeats 81, seedLength 18,  MAXJUMP 3, maxScore 0, restored step=1 on < 60bp
setenv SVlast v37.oct16
setenv SV v42.nov2        # maxTargetRepeats 81, seedLength 18,  MAXJUMP 3, Yann, nA = nB = 40, numactl pin to best
setenv SVlast v42.nov2
setenv SV v43.nov2        # maxTargetRepeats 31, seedLength 18,  MAXJUMP 3, Yann, nA = nB = 40, numactl pin to best
setenv SVlast v43.nov2
setenv SV v44.nov12        # yes wiggle
setenv SVlast v44.nov12
setenv SV v45.nov21        # no wiggle, no sam,    fixed sam and pairs (but i do not use --sam ) , added new T2T index to Hisat, star and minimap
setenv SVlast v45.nov21
#setenv SV v46.dec1        # no wiggle, no sam, fixed pairs again => fixed multiAligned, modified sam2gold to not count introns
#setenv SVlast v46.dec1
setenv SV v47.dec9        # no wiggle, no sam, removed enhanced exon seeds, completed intron discovery and using seeds
#setenv SVlast v47.dec9
setenv SV v48.dec10        #no wiggle
#setenv SVlast v48.dec10
setenv SV v50.jan1         # no wiggle, use the gff to assess the introns while createIndex (do not parse the gff while aligning)
#setenv SVlast v50.jan1
#setenv SV v51.jan4         # no wiggle, same as 50 
#setenv SVlast v51.jan4
setenv SV v52.jan7         # no wiggle, aligned the BigArrays, complete compiler clean up
#setenv SVlast v52.jan7
setenv SV v53.jan20         # no wiggle, no introns (bug), fixing errors close to introns which created sam bugs this index has a NO INTRON seeds (gff error)
setenv SVlast v53.jan20
setenv SV v54.jan22         # no wiggle, reintroduced the introns, fixed  BAM lots of issues on exon boundaries: this index has a Sort ERROR  max31
setenv SVlast v54.jan22
setenv SV v55.81.jan23      # no wiggle, introns, fixing for BAM lots of issues on exon boundaries: this index has a Sorting ERROR  max81
setenv SVlast v55.81.jan23
setenv SV v56.81.jan23      # no wiggle, introns, fixed the error in the seed sorter:  max81
setenv SVlast v56.81.jan23
setenv SV v57.31.jan25      # no wiggle, introns, :  max31
setenv SVlast v57.31.jan25
setenv SV v58.31.16M.jan26      # no wiggle, introns, :  max31
setenv SVlast v58.31.16M.jan26
setenv SV v59.81.16M.jan26      # no wiggle, introns, :  max31
setenv SVlast v59.81.16M.jan26
setenv SV v60.81.18M.Wiggle.jan26      #    wiggle, introns, :  max31
setenv SVlast v60.81.18M.Wiggle.jan26
setenv SV     v61.81.18M.errCost4.jan27      # no wiggle, introns, :  max81
setenv SVlast v61.81.18M.errCost4.jan27zz
setenv SV     v62.81.18M.errCost4.feb4      # no wiggle, no sam , fix extendHits bug loosing errors in cdnaalign :  max81
setenv SVlast v62.81.18M.errCost4.feb4
setenv SV     v63.81.18M.errCost8.feb4      # no wiggle, no sam , fix extendHits bug loosing errors in cdnaalign :  max81
setenv SVlast v63.81.18M.errCost8.feb4
setenv SV     v64.31.18M.errCost4.feb6      # no wiggle, no sam , fix extendHits bug loosing errors in cdnaalign :  max81
setenv SVlast v64.31.18M.errCost4.feb6
#setenv SV     v65.31.18M.Wiggle.errCost4.feb7      # wiggle, no sam , new index for gene expressions
setenv SVlast v65.31.18M.Wiggle.errCost4.feb7
setenv SV     v66.31.18M.NoWiggle.errCost4.feb8    # adaptor clipping no wiggle, no sam , new index for gene expressions
setenv SVlast v66.31.18M.NoWiggle.errCost4.feb8
setenv SV     v67.81.18M.Wiggle.errCost4.feb9    # adaptor clipping wiggle, DNA/RNA automatic, no sam , new index for gene expressions
setenv SVlast v67.81.18M.Wiggle.errCost4.feb9
setenv SV     v68.81.18M.Wiggle.errCost4.feb15    # adaptor clipping wiggle, binary printf DNA/RNA automatic, no sam , new index for gene expressions
setenv SVlast v68.81.18M.Wiggle.errCost4.feb15
setenv SV     v69.81.18M.Wiggle.errCost4.feb16    # adaptor clipping wiggle, binary printf DNA/RNA automatic, no sam , new index for gene expressions
setenv SVlast v69.81.18M.Wiggle.errCost4.feb16
setenv SV     v70.81.18M.e4.W.feb21    # new load balancer for numactl
setenv SVlast v70.81.18M.e4.W.feb21
setenv SV     v71.81.18M.e4.W.feb22    # edited intron boundaries using the new randomTarget
setenv SVlast v71.81.18M.e4.W.feb22

if ($SV == $SVlast) then
  \cp  /home/mieg/ace/bin.LINUX_4_OPT/sortalign bin/sortalign.$SV
  if (-e bin.nativeZZZ/sortalign) \cp bin.native/sortalign bin/sortalign.$SV
  touch bin/wsa.$SV.toto
  \rm -rf bin/wsa.$SV*
  mkdir bin/wsa.$SV
  cp ~/ace/wsa/sa.*[ch] bin/wsa.$SV
endif
\cp bin/sortalign.$SV bin/sortalign

setenv NOINTRONSEEDS 0
setenv EXPORTSAM 0
setenv EXPORTWIGGLES 1
setenv EXPORTWIGGLEENDS 0
setenv seedLength 16
setenv seedLength 18
setenv maxTargetRepeats 31
setenv maxTargetRepeats 81



##       SortalignPaperMasterScript.tcsh
## Author, Greg Boratyn, Danielle Thierry-Mieg, Jean Thierry-Mieg
## email for this script:   mieg@ncbi.nlm.nih.gov

## This is a tcsh executable script, distributed as part of ($MAGIC_SRC)/waligner/scripts/*.tcsh
## To see the on line help, run it under the tcsh interpretor using the command
##       SortalignPaperMasterScript.tcsh

if ($# == 0) goto phase_Help
if ($1 == help) goto phase_Help
if ($1 == '-help') goto phase_Help
if ($1 == '--help') goto phase_Help

## This script is a nerwer version of the 2018 MagicBLAST_paper_master_script.tcsh

# git clone https://github.com/ncbi/magicblast.git
echo "#### SORTALIGN $SV"


#############################################################################
## Metadata

## Aligners
# List of aligners used in the analysis
# The number in front serves to order the tables in a systematic way
# one can insert new versions of each program by inserting new numbers
# but since the numbers are erased in the final tables, the number AND the names 
# must be unique
setenv methods "011_SortAlignG3R5"
setenv methods "012_SortAlignG5R1"
setenv methods "013_SortAlignG3R1"
setenv methods "02_ClipAlign"
setenv methods "10_MagicBLAST"
setenv methods "21_HISAT2"
setenv methods "31_STARlong"
setenv methods "50_Minimap2"


setenv allMethods "011_SortAlignG6R3 012_SortAlignG3R3 013_SortAlignG3R1 014_SortAlignG3R3.g 015_SortAlignG3R1.g 11_MagicBLAST_2018 12_MagicBLAST_2022 13_MagicBLAST_2024 21_HISAT2_4threads 22_HISAT2_8threads 23_HISAT2_16threads 31_STARlong 50_Minimap2 51_Minimap2_4threads 52_Minimap2_8threads 53_Minimap2_16threads 54_Minimap2_32threads"

setenv methods "$allMethods"
#setenv methods "012_SortAlignG3R3"

set createIndex=0
if ($createIndex == 1)  setenv methods "011_SortAlignG6R3 012_SortAlignG3R3"
if ($createIndex == 11)  setenv methods "011_SortAlignG6R3"
if ($createIndex == 12)  setenv methods "012_SortAlignG3R3"

#############################################################################
## Datasets
## Each dataset must be aligned on the reference genome carrying the relevant truth
# Experimental human datasets, to be aligned on the NCBI human genome
setenv main_runs "iRefSeq PacBio Roche Illumina"

# Baruzzo data
# sets, to be aligned on the Baruzzo human reference genome
setenv HG19_r1_runs "HG19t1r1 HG19t2r1 HG19t3r1"
setenv HG19_r2_runs "HG19t1r2 HG19t2r2 HG19t3r2"
setenv HG19_r3_runs "HG19t1r3 HG19t2r3 HG19t3r3"
setenv HG19_runs "$HG19_r1_runs $HG19_r2_runs $HG19_r3_runs"

# Baruzzo datasets, to be aligned on the Baruzzo malaria reference genome
setenv PFAL_r1_runs "PFALt1r1 PFALt2r1 PFALt3r1"
setenv PFAL_r2_runs "PFALt1r2 PFALt2r2 PFALt3r2"
setenv PFAL_r3_runs "PFALt1r3 PFALt2r3 PFALt3r3"
setenv PFAL_runs "$PFAL_r1_runs $PFAL_r2_runs $PFAL_r3_runs"



setenv bigIlluminaRuns "RNA_PolyA_A_1_50Gb RNA_PolyA_B_1_47Gb "

setenv ABRuns "A_100PE B_100PE A_polyA_151PE B_polyA_151PE A_Total_151PE B_Total_151PE B_Roche454_352 A_Nanopore.1234 B_Nanopore.1306 A_PacBio.1745 B_PacBio.1810 A_Total_PacBio.2822 HG02_Revio_PacBio.2044 "
setenv iRefSeqRuns "iRefSeq38 iRefSeq38.g"
setenv refseqRuns "iRefSeq38 iRefSeq38.g iRefSeqT2T iRefSeqT2T.g"
setenv dnaRuns "DNA_ChipSeq1_SE35 DNA_ChipSeq2_SE35 DNA_B_151PE"
setenv wormRuns "  Worm_RNA_150PE "
setenv monkeyRuns "Chimpanzee_0.92pc Macaque_2.69pc PTMacaque_2.71pc Baboon_2.92pc Marmoset_3.14 SquirrelMonkey_3.08pc MouseLemur_2.64pc"
setenv baruzzoRuns "HG19t1_0.543pc HG19t2_1.186pc HG19t3_6.024pc"

setenv allRuns "$ABRuns  $iRefSeqRuns $dnaRuns $wormRuns  $monkeyRuns  $baruzzoRuns "
setenv runs    "$ABRuns  $iRefSeqRuns $dnaRuns $wormRuns  $monkeyRuns  $baruzzoRuns "


# setenv runs "RNA_PolyA_AB_1_97G"
# setenv runs "Roche"
# setenv runs "iRefSeq38"

# to create all IDX use these runs
if ($createIndex > 0)  setenv runs "iRefSeqT2T iRefSeq38 HG19t1_0.543pc  Worm_RNA_150PE"

echo "### S.tcsh SV=$SV"
echo "$methods"
echo "$runs"


# This adapter is present in the PacBio SRR runs and gives a peak at 32 aligned bases = polyA + first 8 bp of adaptor
# AAAAAAAAAAAAAAAAAAAAAAAAAAAAGTACTCT  GCGTTGATACCACTGCTTAGATCGGAAGAG
#############################################################################
## Fasta and Fastq files
## All runs fasta or fastq files are in the directories Fasta/$run
##  If they are absent, the script will download them from NCBI
#
#  The script assumes that all files are gzipped, and called $run/$run*.fast[aq].gz
#  their logical name, i.e. PacBio, links to their SRA identificator, e.g. SRR5009494.
#    The iRefSeq and the Baruzzo files are given in fasta format
#    The Illumina, Pabio and Roche fils are given in .fastq format
#    Some runs are paired-ends:
#       -Illumina paired end run has 2 files called SRR534301_1 and SRR534301_2
#       -Baruzzo paired end runs are called .forward and .reverse
#    In the Roche file Fasta/Roche/SRR899420.fastq we removed all read_1 (all 4-bases long)
#      and kept only the 24577 read_2.
#
#  For convenience and completeness, we also copied in Fasta/iRefSeq the original gff file
#  The iRefSeq fasta file can be extracted from the gff file using the command option
#       MagicBLAST_paper_master_script.tcsh make_iRefSeq 
#  of the the present script

foreach run ($runs)
  if (! -d Fasta/$run) mkdir -p  Fasta/$run
end 

#############################################################################
## Reference genome
# The main genome is T2t release 2024, limited to the main chromosomes
setenv main_genome GRCh38.genome.fasta.gz
setenv main_genome T2T.genome.fasta.gz


# Baruzzo benchmark reference genome, available from
setenv HG19_genome HG19.Baruzzo.genome.fasta.gz
setenv PFAL_genome PFAL.Baruzzo.genome.fasta.gz

## Automatic download of the BenchMark fastq files from NCBI
set FTP="ftp://ftp.ncbi.nlm.nih.gov/blast/demo/magicblast_article/"
if (! -d Reference_genome) mkdir Reference_genome
pushd Reference_genome
  foreach ff (HG19.Baruzzo.genome.TM.txt  HG19.Baruzzo.genome.fasta.gz  PFAL.Baruzzo.genome.TM.txt  PFAL.Baruzzo.genome.fasta.gz)
    if (! -e $ff) then
      wget $FTP/Reference_genome/$ff
    endif
  end
popd

setenv T2T_runs "iRefSeq"

if (1) then
  foreach run ($T2T_runs)
    if (! -d Fasta/$run || -e Fasta/$run/target) continue
    pushd Fasta/$run
      ln -s ../../Reference_genome/T2T.genome.fasta.gz genome.gz
      echo T2T > target
    popd
  end
  foreach run ($HG19_runs)
    if (! -d Fasta/$run || -e Fasta/$run/target) continue
    pushd Fasta/$run
      ln -s ../../Reference_genome/HG19.Baruzzo.genome.fasta.gz genome.gz
      echo HG19 > target
    popd
  end
  foreach run ($PFAL_runs)
    if (! -d Fasta/$run || -e Fasta/$run/target) continue
    pushd Fasta/$run
      ln -s ../../Reference_genome/PFAL.Baruzzo.genome.fasta.gz genome.gz
      echo PFAL > target
    popd
  end
  touch Fasta/genomes_are_assigned

endif

#############################################################################
## Automatic download of the binaries from the NCBI ftp site

# /am/ftp-blast/demo/magicblast_article
set FTP="ftp://ftp.ncbi.nlm.nih.gov/blast/demo/magicblast_article"
if (! -d bin || ! -d scripts/HTSeq) then
  if (-e bin/binaries.linux64.tar.gz) then
    mv binaries.linux64.tar.gz .
  endif
  if (! -e binaries.linux64.tar.gz) then  
     wget $FTP/binaries.linux64.tar.gz
  endif
  if (! -e binaries.linux64.tar.gz) then  
    echo "FATAL ERROR: The automatic download of the binaries from $FTP/binaries.linux64.tar.gz failed"
    echo "May be the network connection did not work, please try manually to run the command"
    echo "    wget $FTP/binaries.linux64.tar.gz"
    echo "if it does not work please email mieg@ncbi.nlm.nih.gov"
  endif
  if (-e binaries.linux64.tar.gz) then
    echo "expanding binaries.linux64.tar.gz, please wait"
    gunzip -c binaries.linux64.tar.gz | tar xf -
    mv binaries.linux64.tar.gz bin
  endif
endif

#############################################################################
## BAM files
## The BAM files are named    $method/$run/$method.$run.bam
## All runs were aligned on the relevant appropriate genome by all aligners
## but it did not always work. Some files are missing, for example 30_STAR.iRefSeq.bam, 
## because STAR crashed on long reads

##############################################################################
## utilities
setenv TMPDIR /tmp
if (-d /export/home/TMP) setenv TMPDIR /export/home/TMP
if (! -d tmp) mkdir tmp
if (! -d RESULTS) mkdir RESULTS

##############################################################################
##############################################################################
## Executable and source code
## Our scripts are in the scripts directory, e.g. scripts/submit
## Our executable are compiled for generic LINUX 64 bits machine in the bin directory
##    e.g. magicblast, dna2dna, sam2gold 
## Our source code is available for analysis and recompilation in machine optimized mode
## in the source_code directory, together with instructions in the correspondng README file.

echo -n "## SortAlignPaperMasterScipt.tcsh $1 : "
date

echo "runs = $runs"
echo "methods = $methods"


if ($1 == init) goto phase_Init
if ($1 == download) goto phase_Download
if ($1 == align) goto phase_Align
if ($1 == Make_iRefSeq) goto phase_Make_iRefSeq
if ($1 == countFasta) goto phase_countFasta
if ($1 == samDownLoad) goto phase_SamDownLoadFromNCBI
if ($1 == timing) goto phase_Timing
if ($1 == accuracy) goto phase_Accuracy
if ($1 == introns) goto phase_introns
if ($1 == aliqc) goto phase_aliqc
if ($1 == perfect) goto phase_perfect
if ($1 == errors) goto phase_DirectErrorCount
if ($1 == export) goto phase_Export
if ($1 == aliLn) goto phase_aliLn
if ($1 == subs) goto phase_count_subtitutions_in_benchmark
echo "Unknown command : $1, please try --help"
goto phaseLoop

phase_Help:

echo "\nusage  scripts/SortAlignPaperMasterScipt.tcsh command, where command is one of:"
echo 'help   : This online help'
echo 'init   : in tcsh, "source README init" will set the variables $runs, $methods which may be convenient'
echo 'download  : Automatic download of the fastq files, please monitor carefully the results'
echo 'Make_iRefSeq : create the iRefSeq fasta file and intron file from the gff and the genome file'
echo 'countFasta  : Count the reads in each fasta/fastq file'
echo 'sam    : Automatically download the sam files from NCBI (rather that running the aligners locally)'
echo 'align  : run for all runs all aligners for which the script Aligners/$method/align.tcsh is defined'
echo 'aliqc : count reads and errors in all the sam files'
echo 'accuracy : measure intron and alignment accuracy relative to the gold standard'
echo 'errors : count the errors in the BAM files using the NM:i:x optional field'
echo 'subs   : count substitutions in the human and malaria Baruzzo benchmarks'
echo 'aliLn  : export histogram of aligned lengths'
echo 'export : export QC and ROC curve of intron discovery'

goto phaseLoop

phase_Init:
goto phaseLoop

##############################################################################
#############################################################################
## Automatic download of the BenchMark fastq files from NCBI
phase_Download:

## NCBI repository for the datafiles used in the paper
set FTP="ftp://ftp.ncbi.nlm.nih.gov/blast/demo/magicblast_article/"

echo "checking $HG19_runs  $PFAL_runs"
foreach run ( $HG19_runs  $PFAL_runs)
  if (! -e Fasta/$run/$run.reverse.fasta.gz) then
    pushd Fasta/$run
      echo "Trying to download $run from$FTP"
      wget $FTP/Fasta/$run/$run.cig.gz
      wget $FTP/Fasta/$run/$run.forward.fasta.gz
      wget $FTP/Fasta/$run/$run.reverse.fasta.gz
    popd
  endif
end
set run=HG19t1r1_50
if (! -d Fasta/$run) then
  echo "preparing the 50+50 clipped run"
  mkdir Fasta/$run
  pushd  Fasta/$run
    echo HG19 > target
    ln -s ../../Reference_genome/HG19.Baruzzo.genome.fasta.gz genome.gz
    ../../bin/dna2dna -i ../HG19t1r1/HG19t1r1.forward.fasta.gz -gzo -o $run.forward -rightClipAt 50
    ../../bin/dna2dna -i ../HG19t1r1/HG19t1r1.reverse.fasta.gz -gzo -o $run.reverse -rightClipAt 50
  popd
endif

set run=SRR5438850_120
if (! -d Fasta/$run) then
  echo "preparing the 120+120 clipped run"
  mkdir Fasta/$run
  pushd  Fasta/$run
    echo GRCh38 > target
    ln -s ../../Reference_genome/GRCh38.genome.fasta.gz genome..gz
    ln -s ../iRefSeq/iRefSeq.cig.gz $run.cig.gz
    ../../bin/dna2dna -i ../SRR5438850/SRR5438850_1.fasta.gz -gzo -o $run.forward -rightClipAt 120
    ../../bin/dna2dna -i ../SRR5438850/SRR5438850_2.fasta.gz -gzo -o $run.reverse -rightClipAt 120
  popd
endif

set run=SRR5438850_150
if (! -d Fasta/$run) then
  echo "preparing the 150+150 clipped run"
  mkdir Fasta/$run
  pushd  Fasta/$run
    echo GRCh38 > target
    ln -s ../../Reference_genome/GRCh38.genome.fasta.gz
    ln -s ../iRefSeq/iRefSeq.cig.gz $run.cig.gz
    ../../bin/dna2dna -i ../SRR5438850/SRR5438850_1.fasta.gz -gzo -o $run.forward -rightClipAt 150
    ../../bin/dna2dna -i ../SRR5438850/SRR5438850_2.fasta.gz -gzo -o $run.reverse -rightClipAt 150
  popd
endif

set run=SRR5437876_120
if (! -d Fasta/$run) then
  echo "preparing the 120+120 clipped run"
  mkdir Fasta/$run
  pushd  Fasta/$run
    echo GRCh38 > target
    ln -s ../../Reference_genome/GRCh38.genome.fasta.gz
    ln -s ../iRefSeq/iRefSeq.cig.gz $run.cig.gz
    ../../bin/dna2dna -i ../SRR5437876/SRR5437876_1.fasta.gz -gzo -o $run.forward -rightClipAt 120
    ../../bin/dna2dna -i ../SRR5437876/SRR5437876_2.fasta.gz -gzo -o $run.reverse -rightClipAt 120
  popd
endif

set run=SRR5437876_150
if (! -d Fasta/$run) then
  echo "preparing the 150+150 clipped run"
  mkdir Fasta/$run
  pushd  Fasta/$run
    echo GRCh38 > target
    ln -s ../../Reference_genome/GRCh38.genome.fasta.gz
    ln -s ../iRefSeq/iRefSeq.cig.gz $run.cig.gz
    ../../bin/dna2dna -i ../SRR5437876/SRR5437876_1.fasta.gz -gzo -o $run.forward -rightClipAt 150
    ../../bin/dna2dna -i ../SRR5437876/SRR5437876_2.fasta.gz -gzo -o $run.reverse -rightClipAt 150
  popd
endif

set run=PFALt1r1S
if (! -d Fasta/$run) then
  echo "preparing the subsampled run"
  mkdir Fasta/$run
  pushd  Fasta/$run
    echo PFAL > target
    ln -s ../../Reference_genome/PFAL.Baruzzo.genome.fasta.gz genome.gz
    ../../bin/dna2dna -i ../PFALt1r1/PFALt1r1.forward.fasta.gz -gzo -o $run.forward -subsample 100
    ../../bin/dna2dna -i ../PFALt1r1/PFALt1r1.reverse.fasta.gz -gzo -o $run.reverse -subsample 100
  popd
endif

echo "checking iRefSeq"
foreach run (iRefSeq)
  if (! -e Fasta/$run/$run.fasta.gz) then
    pushd Fasta/$run
      echo "Trying to download $run from$FTP"
      wget $FTP/Fasta/$run/$run.cig.gz
      wget $FTP/Fasta/$run/$run.fasta.gz
      wget $FTP/Fasta/$run/GRCh38_genome.gff.gz
      ln -s GRCh38_genome.gff.gz genome.gff.gz
      gunzip -c $run.fasta.gz | ../../bin/dna2dna -getTM > $run.TM
    popd
  endif
end

echo "checking Roche"
foreach run (Roche)
  if (! -e Fasta/$run/$run.fasta.gz) then
    pushd Fasta/$run
      echo "Trying to download $run from$FTP"
      wget $FTP/Fasta/$run/$run.fasta.gz
      gunzip -c $run.fasta.gz | ../../bin/dna2dna -getTM > $run.TM
    popd
  endif
end

#############################################################################
## Automatic download of the Illumina Roche pacBio from SRA

foreach run (PacBio Illumina)
  if (! -d Fasta/$run) continue
  if ($run == PacBio) set run2=SRR5009494
  if ($run == Roche)  set run2=SRR899420
  if ($run == Illumina)  set run2=SRR534301
  if (! -e Fasta/$run/$run2.fastq.gz && ! -e Fasta/$run/$run2'_1'.fastq.gz) then
     set n=`bin/fastq-dump --help | wc -l`
    if ($n < 10) then
      echo "Sorry, the executable bin/fastq-dump available from NCBI SRA and needed to download the $run run is not found"
      goto phaseLoop
    endif
    set sf=""
    if ($run == Illumina)  set sf="--split-files"
    echo "Trying to download $run2 from SRA"
    bin/fastq-dump $sf -O Fasta/$run $run2
    if (-e Fasta/$run/$run2.fastq || -e Fasta/$run/$run2'_1'.fastq) then
      gzip Fasta/$run/$run2*.fastq
      pushd Fasta/$run
        if (-e $run2.fastq.gz) ln -s $run2.fastq.gz $run.fastq.gz
        if (-e $run2'_1'.fastq.gz) ln -s $run2'_1'.fastq.gz $run'_1'.fastq.gz
        if (-e $run2'_2'.fastq.gz) ln -s $run2'_2'.fastq.gz $run'_2'.fastq.gz  
      popd
    endif
    if (-e ~/ncbi/public/sra/$run2) \rm -rf  ~/ncbi/public/sra/$run2
  endif
end

#############################################################################
## Automatic download of the fastq files from SRA

foreach run ($pacbio_runs)
  if (! -d Fasta/$run) continue
  if (! -e Fasta/$run/$run.fasta.gz && ! -e Fasta/$run/$run.fastq.gz) then
    set n=`bin/fastq-dump --help | wc -l`
    if ($n < 10) then
      echo "Sorry, the executable bin/fastq-dump available from NCBI SRA and needed to download the $run run is not found"
      goto phaseLoop
    endif
    echo "Trying to download $run from SRA"
    bin/fastq-dump -O Fasta/$run $run
    gzip Fasta/$run/$run.fastq
  endif
end

foreach run ($long_illumina_runs)
  if (! -d Fasta/$run) continue
  if (! -e Fasta/$run/$run.fasta'_1'.gz && ! -e Fasta/$run/$run'_1'.fastq.gz) then
    set n=`bin/fastq-dump --help | wc -l`
    if ($n < 10) then
      echo "Sorry, the executable bin/fastq-dump available from NCBI SRA and needed to download the $run run is not found"
      goto phaseLoop
    endif
    echo "Trying to download $run2 from SRA"
    bin/fastq-dump -O Fasta/$run --split-files $run 
    gzip Fasta/$run/$run*.fastq
  endif
end

goto phaseLoop

#############################################################################
## Automatic download of the AGLR2 fasta  files from SRA
# new worm
# Worm_on_Worm_75PE  DRR412604
# Worm_on_Worm_150PE   SRR26453245

# --maxGb did not work, i obtained 10Gb
bin/sortalign --sraDownload SRR24511885,SRR24511871,SRR24511867,SRR24511881,SRR24511877 --maxGB 2
bin/sortalign --sraDownload SRR26453245,DRR412604 --maxGB 2
bin/sortalign --sraDownload SRR24511885,SRR24511871,SRR24510080,SRR24510074,SRR24511773,SRR24511767,SRR24510080,SRR24510074,SRR899420,SRR24564586,SRR24564582,SRR29438430,SRR24511720,SRR24564588,SRR24564584,SRR24518509,SRR26453245,SRR1758917,SRR1758936,SRR1759007,SRR1758904,SRR1758979,SRR1759037,SRR1758993 --maxGB 2
 
A_capt_100PE:SRR24511885
B_capt_100PE:SRR24511871
A_polyA_151PE:SRR24510080
B_polyA_151PE:SRR24510074
A_Total_151PE:SRR24511773	
B_Total_151PE:SRR24511767	
A_polyA_151PE_50Gb:SRR24510080
B_polyA_151PE_47Gb:SRR24510074
B_Roche454_352:SRR899420

11_A_ROCR3_PacBio-F3	A_capt_PacBio.1745	SRR24564586	1745
12_B_ROCR3_PacBio-F3	B_capt_PacBio.1810	SRR24564582	1810
29_SRR29438430	HG02_Revio_PacBio.2044	SRR29438430	av 2044
10_A_WTS_PacBio	A_Total_PacBio.2822	SRR24511720	2822
13_A_ROCR3_Nanopore-F3	A_capt_Nanopore.1234	SRR24564588	1234
14_B_ROCR3_Nanopore-F3	B_capt_Nanopore.1306	SRR24564584	1306
1_iRefSeq38	iRefSeq_38		
2_iRefSeq38.g	iRefSeq_38_g		
5_ChipSeq1	DNA_ChipSeq1_SE35	 	
6_ChipSeq2	DNA_ChipSeq2_SE35	 	
7_B_ROCR2_Illumina_DNA	DNA_B_capt_151PE	SRR24518509	
 	 		
26.WormSRR548309	Worm_RNA_150PE	SRR26453245	
16_FrontalCortex_CHP_Chimpanzee	Chimpanzee_0.92pc	SRR1758917	
18_FrontalCortex_CMC_Macaque	Macaque_2.69pc	SRR1758936	
19_TemporalLobe_PTM_Macaque	PTMacaque_2.71pc	SRR1759007	
17_TemporalLobe_BAB_Baboon	Baboon_2.92pc	SRR1758904	
20_BrainRight_MST_Marmoset	Marmoset_3.14	SRR1758979	
22_TemporalLobe_SQM_SquirrelMonkey	SquirrelMnkey_3.08pc	SRR1759037	
21_TemporalLobe_MLM_MouseLemur	MouseLemur_2.64pc	SRR1758993	
23_HG19t1r1	23_HG19t1_0.543pc	synthetic truth (Baruzzo)	
24_HG19t2r1	24_HG19t2_1.186pc		
25_HG19t3r1	25_HG19t3_6.024pc		



foreach rSrr (A_100PE:SRR24511885 B_100PE:SRR24511871 A_polyA_151PE:SRR24510080 B_polyA_151PE:SRR24510074 A_Total_151PE:SRR24511773 B_Total_151PE:SRR24511767 A_polyA_151PE_50Gb:SRR24510080 B_polyA_151PE_47Gb:SRR24510074 B_Roche454_352:SRR899420 A_Nanopore.1234:SRR24564588 B_Nanopore.1306:SRR24564584 A_PacBio.1745:SRR24564586 B_PacBio.1810:SRR24564582 A_Total_PacBio.2822:SRR24511720 HG02_Revio_PacBio.2044:SRR29438430 iRefSeq_38:iRefSeq DNA_B_151PE:SRR24518509 Worm_RNA_150PE:SRR26453245 Chimpanzee_0.92pc:SRR1758917 Macaque_2.69pc:SRR1758936 PTMacaque_2.71pc:SRR1759007 Baboon_2.92pc:SRR1758904 Marmoset_3.14:SRR1758979 SquirrelMonkey_3.08pc:SRR1759037 MouseLemur_2.64pc:SRR1758993)

foreach rSrr (A_100PE:SRR24511885)
   set target=T2T
   set run=`echo $rSrr | gawk -F ':' '{print $1;}'`
   set srr=`echo $rSrr | gawk -F ':' '{print $2;}'`
   if ($run == Worm_RNA_150PE) set target=worm
   mkdir NewFasta/$run
   echo $srr >  NewFasta/$run/SRR
   echo $run SRA/$srr.sra.fasta.gz
   if (! -e NewFasta/$run/$srr.sra.fasta.gz) then
     zcat SRA/$srr.sra.fasta.gz | gawk 'BEGIN{n=0;}{if(n%4 < 2)print;n++;if(n>=40000000)exit;}' | gzip > NewFasta/$run/$srr.forward.fasta.gz &
     zcat SRA/$srr.sra.fasta.gz | gawk 'BEGIN{n=2;}{if(n%4 >= 2)print;n++;if(n>=40000002)exit;}' | gzip > NewFasta/$run/$srr.reverse.fasta.gz &
     zcat SRA/$srr.sra.fasta.gz | head -40000000 | gzip > NewFasta/$run/$srr.sra.fasta.gz &
   endif
   pushd NewFasta/$run
     ln -s $srr.sra.fasta.gz  $run.sra.fasta.gz
     echo $srr > SRR
     echo $target > target
     ln -s $srr.forward.fasta.gz $run.forward.fasta.gz
     ln -s $srr.reverse.fasta.gz $run.reverse.fasta.gz
   popd
end

goto phaseLoop

##############################################################################
##############################################################################
## Count the number of reads, the shortest, the longest read in every fasta/fastq file
## using the utility bin/dna2dna (compiled for Linux 64bits) 
## The source code is part of the aceview/magic distribution in the source_code directory

phase_countFasta:
echo 'counting the number of reads in each fasta/fastq file'
foreach run ($runs)
  echo $run
  if (-e Fasta/$run/$run.fasta.gz && ! -e Fasta/$run/$run.count) then
    echo "counting $run, please wait"
    echo "bin/dna2dna -i  Fasta/$run/$run.fasta.gz -I fasta -count -o Fasta/$run/$run"
          bin/dna2dna -i  Fasta/$run/$run.fasta.gz -I fasta -count -o Fasta/$run/$run
  endif
  if (-e Fasta/$run/$run.fastq.gz && ! -e Fasta/$run/$run.count) then
    echo "counting $run, please wait"
    echo "bin/dna2dna -i  Fasta/$run/$run.fastq.gz -I fastq -count -o Fasta/$run/$run"
          bin/dna2dna -i  Fasta/$run/$run.fastq.gz -I fastq -count -o Fasta/$run/$run
  endif
  if (-e Fasta/$run/$run'_1'.fastq.gz && ! -e Fasta/$run/$run'_1'.count) then
    echo "counting $run, please wait"
    bin/dna2dna -i  Fasta/$run/$run'_1'.fastq.gz -I fastq -count -o Fasta/$run/$run'_1'
    bin/dna2dna -i  Fasta/$run/$run'_2'.fastq.gz -I fastq -count -o Fasta/$run/$run'_2'
  endif
  if (-e Fasta/$run/$run'_1'.fasta.gz && ! -e Fasta/$run/$run'_1'.count) then
    echo "counting $run, please wait"
    bin/dna2dna -i  Fasta/$run/$run'_1'.fasta.gz -I fasta -count -o Fasta/$run/$run'_1'
    bin/dna2dna -i  Fasta/$run/$run'_2'.fasta.gz -I fasta -count -o Fasta/$run/$run'_2'
  endif
  if (-e Fasta/$run/$run'_R1'.fastq.gz && ! -e Fasta/$run/$run'_1'.count) then
    echo "counting $run, please wait"
    bin/dna2dna -i  Fasta/$run/$run'_R1'.fastq.gz -I fastq -count -o Fasta/$run/$run'_1'
    bin/dna2dna -i  Fasta/$run/$run'_R2'.fastq.gz -I fastq -count -o Fasta/$run/$run'_2'
  endif
  if (-e Fasta/$run/$run.forward.fasta.gz && ! -e Fasta/$run/$run.forward.count) then
    echo "counting $run, please wait"
    bin/dna2dna -i  Fasta/$run/$run.forward.fasta.gz -I fasta -count -o Fasta/$run/$run.forward 
    bin/dna2dna -i  Fasta/$run/$run.reverse.fasta.gz -I fasta -count -o Fasta/$run/$run.reverse 
  endif
  set nreads=`cat Fasta/$run/*.count  | gawk '/^Fragment_kept/{n+=$2}END{print n}'`
  echo "$run contains $nreads reads"
end

 set toto=Fasta/RunStat.txt
 echo -n "### toto : " > $toto
 date >> $toto
 echo "# Run\tnReads\tAverage read length\tMinimal read length\tMaximal read length\tBases" >> $toto
 foreach run ($runs)
   echo -n "$run" >> $toto
   echo ' ' > toto
   cat Fasta/$run/$run.count >> toto
   cat Fasta/$run/$run'_'[12].count >> toto
   cat Fasta/$run/$run.*.count >> toto
   cat toto | gawk '/Tags_processed/{nr+=$2;}/Bases_tags_processed/{bp+=$2;}/Min_probe_length/{x=$2;if(minx<1 || minx > x)minx=x;}/Max_probe_length/{x=$2;if(maxx+0<x)maxx=x;}END{printf("\t%d\t%d\t%d\t%d\t%d\n",nr,int(bp/(nr+.00001)),minx,maxx,bp) ;}' >> $toto

 end
echo $toto
\cp $toto COMPARE

goto phaseLoop

##############################################################################
##############################################################################
## Create the iRefSeq fasta file and intron file from the gff and the genome file

phase_Make_iRefSeq:

## In practice, the file Fasta/iRefSeq/iRefSeq.fasta.gz is downloaded from
## ftp://ftp.ncbi.nlm.nih.gov/blast/demo/magicblast_article/
## The scrip is given here for transparency and to allow the reconstruction
## of the iRefSeq in the future from a different gff file and reference genome

if (! -e Fasta/iRefSeq/iRefSeq.fasta.gz) then
  echo "Creating Fasta/iRefSeq/iRefSeq.fasta.gz using the genome and the gff3 annotation"
  if (! -e Fasta/iRefSeq/genome.gz) then
    echo "Missing file Fasta/iRefSeq/genome.gz, I cannot create the iRefSeq fasta file"
    goto phaseLoop
  endif
  if (! -e Fasta/iRefSeq/genome.gff.gz( then
     echo "Missing file Fasta/iRefSeq/genome.gff.gz, I cannot create the iRefSeq fasta file"
    goto phaseLoop
  endif

  echo "Found the genome and the gff file, constructing the fasta in Fasta/iRefSeq/tmp"
  if (! -d Fasta/iRefSeq/tmp) mkdir  Fasta/iRefSeq/tmp
  pushd  Fasta/iRefSeq/tmp
    # This script is surprisingly complex, sorry, because we are trying to identify the NMs which map as well
    # at different locus of the genome, but while doing so, we unfortunately discovered a number of
    # irregularities in the definition of the RefSeqs that we try to compensate

    # To simplify the matter, we directly provide the iRefSeq fasta and gff files on our ftp site. 

    # extract the NM_ from the gff file
    zcat ../genome.gtf.gz | grep NM_ | grep NC_ > NM.gtf
    # we could directly export the fasta file with the command 
    # ../../../bin/dna2dna -gff3 NM.gff -gtfRemap iRefSeq -gtfGenome ../genome.gz -o iRefSeq -O fasta
    # but some NM have a single indentifier and yet map on 2 chroms
    # by not providing the genome we only export the 6 columns sponge file
    ../../../bin/dna2dna -gtf NM.gtf -gtfRemap iRefSeq -o iRefSeq
    set nNM=`cat iRefSeq.[fr].sponge | cut -f 1 | sort -u | wc -l`
    echo  "Number of NM_ $nNM (is 45065)"
    set nNM_chr=`cat iRefSeq.[fr].sponge | cut -f 1,3 | sort -u | wc -l`
    echo  "Number of NM_chrom $nNM_chr (is 67046)"
    set nG=`cat iRefSeq.[fr].sponge | cut -f 6 | sort -u | wc -l`
    echo  "Number of genes with NM $nG (is 67046)
    echo "Evaluating the mapping multiplicity of the iRefSeq"
    # to fix the issue that the same NM may map on  several chromosomes (none in T2T, 43 cases on GRCh38)
    # we merge the chrom and the NM in column 1 to create a disambiguated sponge file
    cat  iRefSeq.[fr].sponge  | gawk -F '\t' '{printf("%s:%s\t%s\t%s\t%s\t%s\t%s\n",$1,$3,$2,$3,$4,$5,$6);}' > NM_chr.sponge
    # the sponge file has the NM the gene and the coordinates of all exons, hence the sequence
    # we now construct the fasta file
     ../../../bin/dna2dna -sponge NM_chr.sponge -i ../genome.gz -o iRefSeq -O fasta -maxLineLn 80 -gzo
    # measure the number of distinct NM with identical sponge (hence DNA) and mapping    
    cat iRefSeq.[fr].sponge | grep NM_ | sort > _t
    ../../../bin/dna2dna -i iRefSeq.fasta.gz  -count -o iRefSeq
    \mv iRefSeq.fasta.gz iRefSeq.count ..
    # count the distinct NM sequences : 67046 (no doublets in T2T)
    ../../../bin/dna2dna -i ../iRefSeq.fasta.gz -O raw | cut -f 1 | sort -u | wc
    ../../../bin/dna2dna -i ../iRefSeq.fasta.gz -getTM > ../iRefSeq.TM

    # map the NM on the NM to find who is identical or included in the other
    clipalign -i ../iRefSeq.fasta.gz -t ../iRefSeq.fasta.gz -errMax 0 -o nm2nm -maxHit 24 -minAli 140
    bestali -i  nm2nm.hits -exportBest -o nm2nm33
    sortalign -t ../iRefSeq.fasta.gz --createIndex RS
    sortalign -x RS -i  ../iRefSeq.fasta.gz --align -o nm2nm.sortali
    # now count NM mapping exactly in NM with a different geneid -> 143, we add 43+43 for the 43 NM which map on X and Y
    cat nm2nm33.hits | gawk '{if($2-$4==0 && index($1,"|"$9">")==0)print}' > nm2nm.2genes.hits
    wc nm2nm.2genes.hits
    cat nm2nm.2genes.hits | gawk  -F '\t' '{n[$1]++;}END{for(k in n)u[n[k]]++;for (v  in u) {if(v>0)k+=u[v];kk+=v*u[v];print v, u[v];}print k, kk}' | sort -k 1n
    # we now have 303 NM mapping on another NM with a different gene name, however some distinct genes have same gene coordinates
    # extract the extreme coords of the NM from the sponge file
    cat NM_chr.sponge | gawk -F '\t' '{nm=$1;a1=$4;if($5<a1)a1=$5;a2=$5;if($4>a2)a2=$4;if(aa2[nm]<a2)aa2[nm]=a2;if(aa1[nm]>-a1)aa1[nm]=-a1;}END{for(nm in aa1) printf("%s\t%d\t%d\n",nm,-aa1[nm],aa2[nm]);}' > NM_chr.segment
    # reanalize the nm 2 nm hits file and eliminate the lines with overlapping coordinates
    echo ZZZZZ > ZZZZZ
    cat NM_chr.segment ZZZZZ  nm2nm.2genes.hits | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{nm=$1;if(zz+0<1){aa1[nm]=$2;aa2[nm]=$3;split(nm,aa,":");chrom[nm]=aa[2];next;}}{split($1,aa,"|");nm1=aa[1];nm2=$11;ok=1;if (chrom[nm1]==chrom[nm2] && aa1[nm1]<aa2[nm2] && aa2[nm1] > aa1[nm2])ok=0;if (ok==1)print}' > nm2nm.2genes.hits.no_doublons
    cat NM_chr.segment ZZZZZ  nm2nm.2genes.hits | gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{nm=$1;if(zz+0<1){aa1[nm]=$2;aa2[nm]=$3;split(nm,aa,":");chrom[nm]=aa[2];next;}}{split($1,aa,"|");nm1=aa[1];nm2=$11;ok=1;if (chrom[nm1]==chrom[nm2] && aa1[nm1]<aa2[nm2] && aa2[nm1] > aa1[nm2])ok=0;if (ok==0)print}' > nm2nm.2genes.hits.doublons
    # final count of the repeated NM : 291 NM have several mappings
    cat nm2nm.2genes.hits.no_doublons | gawk  -F '\t' '{n[$1]++;}END{for(k in n)u[n[k]]++;for (v  in u) {if(v>0)k+=u[v];kk+=v*u[v];print v, u[v];}print kk}' | sort -k 1n
    # so finally we have 291 NM have multiple mapping just by looking at the annotated NM themselves + (43 + 43)  from the pseudo autosomal region with single NM and geneid total 291+86=379
    cat NM_chr.segment| gawk '{split($1,aa,":");n[aa[1]]++;chrom[aa[1]]=aa[2];}END{for (nm in n)if(n[nm]>1)print nm,n[nm],chrom[nm];}' > NM.pseudo_autosomal_region.mapping_twice 
    wc  NM.pseudo_autosomal_region.mapping_twice     
    cat nm2nm.2genes.hits.no_doublons | gawk '{split($1,aa,"|");print aa[3] "="$9}' | sed -e 's/>/</' | sort -u > gene_pairs
    ## construct a cig file for the refseq
    # Use NM_...:chrom... as NM identifiers because in the pseudo autosomal region, the same NM maps in 2 places: one NM_ 2 locus
    # whereas a usual palindromic exactly duplicated genes has 1 NM per locus, i.e. 2 NM 2 locus
    # this raises the number of NM supported introns from 210357 to 210509

    #############################################################################
## .cig TRUTH Files
## The Baruzzo benchmark is providing the original position of the simulated reads
## in their .cig format, which is analogous, but not identical, to the SAM format.
##     exon_length M intron_length N exon_length M .... in the orientation of the mRNA
## Since the fasta and the .cig files both come from Baruzzo, we located them in Fasta/$run.cig.gz
## For convenience, we reformatted the RefSeq gtf file into a similar Fasta/iRefSeq/iRefSeq.cig.gz
##
## To compare the BAM files produced by the different aligners to the .cig "truth"
## we developped a C code called sam2gold (see below)
   cat NM.gtf | head -100 |  gawk -F '\t' '{if ($3 != "exon") next;}{chrom=$1;a1=$4;a2=$5;strand=$7;z=$9;split(z,aa,"transcript_id");split(aa[2],bb,"\"") ;nm=bb[2];printf("%s\t%s\t%s\t%d\t%d\n",nm,chrom,strand,a1,a2);}' | sort -V | gawk -F '\t' '{nm=$1;chrom[nm]=$2;strand[nm]=$3;a1=$4;a2=$5;da=a2-a1+1;if(dda[nm]+0<1){aa1[nm]=a1;}else{intron=a1-aa2[nm]-1;cig[nm] = cig[nm] intron "N";}if(a2>aa2[nm])aa2[nm]=a2;dda[nm]+=da;cig[nm]=cig[nm] da "M";}END{for (nm in aa1) printf("%s:%s\t%d\t%d\t%s\t.\t%s\t.\n",nm,chrom[nm],aa1[nm],aa2[nm],cig[nm],strand[nm]);}' | gzip > ../iRefSeq.cig.gz
    


    endif

goto phaseLoop

##############################################################################
##############################################################################
## SAM 
## download the precomputed SAM files from NCBI
phase_SamDownLoadFromNCBI:
goto phaseLoop
set FTP="ftp://ftp.ncbi.nlm.nih.gov/blast/demo/magicblast_article/"
foreach run ($runs)
  foreach method ($methods)

    ## the preferred methos is to download the aligned files from NCBI

    # For HISAT and STAR we have a special version of the code to align long runs
    # so in these cases we do not atempt to align the long runs with the short code
    if ($method == 30_STAR || $method == 32_STAR.2.6c) then
       if (-e Fasta/$run/isLongRun) continue
    endif
    # and vice versa
    if ($method == 20_HISAT2_relaxed) then
       if (! -e Fasta/$run/isLongRun) continue
    endif
 
      if (! -e $method/$run/$method.$run.sam_sorted.gz) then
        mkdir -p $method/$run
        pushd  $method/$run
          wget $FTP/SAM/$method.$run.sam_sorted.gz
        popd
      endif

  end
end

goto phaseLoop

##############################################################################
##############################################################################
## ALIGN Run all aligners on all runs

phase_Align:

foreach run ($runs)
  foreach method ($methods)

    if (! -e Aligners/$method/align.tcsh) then
      echo "missing script Aligners/$method/align.tcsh"
      continue
    endif

    if (! -d RESULTS) mkdir RESULTS
    if (! -d RESULTS/$method) mkdir RESULTS/$method
    if (! -d RESULTS/$method/$run) mkdir RESULTS/$method/$run

    # For HISAT and STAR we have a special version of the code to align long runs
    # so in these cases we do not atempt to align the long runs with the short code
    if ($method == 30_STAR || $method == 32_STAR.2.6c) then
       if (-e Fasta/$run/isLongRun) continue
    endif
    # and vice versa
    if ($method == 20_HISAT2_relaxed) then
       if (! -e Fasta/$run/isLongRun) continue
    endif

    if (-e RESULTS/$method/$run/align.running) continue
    if (-e RESULTS/$method/$run/align.done) continue
    
    if (-e RESULTS/$method/$run/sam) continue
    if (-e RESULTS/$method/$run/sam.gz) continue
    if (-e RESULTS/$method/$run/sam_sorted) continue
    if (-e RESULTS/$method/$run/sam_sorted.gz) continue
		
    if (-e Aligners/$method/align.tcsh) then
      set read_1="x"
      set read_2=""
      if (-e Fasta/$run/$run.fasta.gz) set read_1=Fasta/$run/$run.fasta.gz
      if (-e Fasta/$run/$run.fastq.gz) set read_1=Fasta/$run/$run.fastq.gz
      if (-e Fasta/$run/$run.forward.fastq.gz) set read_1=Fasta/$run/$run.forward.fastq.gz
      if (-e Fasta/$run/$run.reverse.fastq.gz) set read_2=Fasta/$run/$run.reverse.fastq.gz
      if (-e Fasta/$run/$run.forward.fasta.gz) set read_1=Fasta/$run/$run.forward.fasta.gz
      if (-e Fasta/$run/$run.reverse.fasta.gz) set read_2=Fasta/$run/$run.reverse.fasta.gz
      if ($method == 011_SortAlignG3R5 && -e Fasta/$run/$run.sra.fasta.gz) set read_1=Fasta/$run/$run.sra.fasta.gz
      if ($method == 012_SortAlignG3R3 && -e Fasta/$run/$run.sra.fasta.gz) set read_1=Fasta/$run/$run.sra.fasta.gz
      if ($method == 013_SortAlignG3R1 && -e Fasta/$run/$run.sra.fasta.gz) set read_1=Fasta/$run/$run.sra.fasta.gz
      if ($method == 014_SortAlignG3R3.g && -e Fasta/$run/$run.sra.fasta.gz) set read_1=Fasta/$run/$run.sra.fasta.gz
      if ($method == 015_SortAlignG3R1.g && -e Fasta/$run/$run.sra.fasta.gz) set read_1=Fasta/$run/$run.sra.fasta.gz
      if ($read_1 == Fasta/$run/$run.sra.fasta.gz) set read_2=""
      
      if (! -e $read_1) then
        echo "Run $run Missing read file $read_1"
	ls -ls Fasta/$run/*fast*
	ls -ls Fasta/$run/$run'_R1.fasta.gz'
        continue
      endif
      set target=`cat Fasta/$run/target`
      if (! -d $method/$run) mkdir -p $method/$run
        echo "Aligners/$method/align.tcsh $method $run $target $read_1 $read_2"
        # scripts/submit $method/$run "Aligners/$method/align.tcsh $method $run $target $read_1 $read_2" local
	touch RESULTS/$method/$run/toto
	\rm RESULTS/$method/$run/*
        scripts/submit RESULTS/$method/$run/$run "Aligners/$method/align.tcsh $method $run $target $read_1 $read_2"  64G
        if (-e RESULTS/$method/$run/s2g.samSats) \rm RESULTS/$method/$run/s2g.samSats
        if (-e RESULTS/$method/$run/$run.s2g.samSats) \rm RESULTS/$method/$run/$run.s2g.samSats
        if (-e RESULTS/$method/$run/$run/s2g.samSats) \rm RESULTS/$method/$run/$run/s2g.samSats
      endif
  end
end

goto phaseLoop

##############################################################################
##############################################################################
## Intron, exon, insertion deletion comparison to the TRUTH
## Compare the alignment results, provided in BAM format 
## to the GOLD standard truth from iRefSeq and Baruzzo given in .cig format
## The source C-code is part of the aceview/magic distribution www.aceview.org/Software
## The executable for LINUX 64 bits is in bin
##
## sam2gold produces several output files
##    .qc a small self documented statistics table in tab delimited format
##    .aliqc the same statistics in a more computer friendly tag-values tab delimited format
##    .Intron a table giving the coordinates of all introns, with support in GOLD or BAM
##    .Deletion a table giving the coordinates of all deletions, with support in GOLD or BAM
##    .Insertion a table giving the coordinates of all insertions, with support in GOLD or BAM

phase_Timing:
echo -n "phase_Timing :"
date

# Create a tsf table
echo -n "### Timings of all methods version : $SV :" > COMPARE/timing.txt
#hs other
foreach species (any)

 set out=RESULTS/compare.$species
 \rm $out.*.tsf
  foreach type (cpu wallT maxMem parallel)
    echo "### $type" > $out.$type.tsf
  end

 setenv speciesRuns "$allRuns"
 if ($species == other) then
  setenv speciesRuns "$monkeyRuns"
 endif
 setenv speciesRuns "$allRuns $monkeyRuns"
 foreach mm ($allMethods)
  foreach run ($runs)

    if (! -e RESULTS/$mm/$run/align.done) continue

    set e=RESULTS/$mm/$run/$run.err
    if (! -e $e) continue

    set wallT=""
    set cpu=""
    set maxMem=""
    set multiP=""

    cat $e | grep TIMING | gawk -F '\t' '{n=split($5,aa,":");if(n==3)s=3600*aa[1]+60*aa[2]+aa[3];else s=60*aa[1]+aa[2];printf("%d\t",s);print;}' | sort -k 1nr | tail -1 > _t
    set n=`cat _t | gawk '/TIMING/{n++}END{print n+0}'`
    if ($n == 1) then
      set wallT=`cat _t | gawk -F '\t' '{print $6}'`
      set cpu=`cat _t | gawk -F '\t' '{print $8+$10}'`
      set maxMem=`cat _t | gawk -F '\t' '{print $12}'`
      set multiP=`cat _t | gawk -F '\t' '{print $14}'`
   endif

    set nRun=`echo "$run ZZZZZ $runs" | gawk '{r=$1;for (i=3;i<=NF;i++)if(r==$i)k=i-2;}END{printf ("%2d", k);}'`
   
    if (! -e RESULTS/$mm/$run/align.done) continue 
    echo "$nRun.$run\t$mm\tt\t$wallT" >> $out.wallT.tsf
    echo "$nRun.$run\t$mm\tt\t$cpu" >> $out.cpu.tsf
    echo "$nRun.$run\t$mm\tt\t$maxMem" >> $out.maxMem.tsf
    echo "$nRun.$run\t$mm\tt\t$multiP" >> $out.parallel.tsf
  end
end

  date >> COMPARE/timing.txt

  set type=wallT
  set title="Overall elapsed time in decimal minutes on the same hardware, as reported by the unix command /usr/bin/time %E"
  echo "### $title" >>  COMPARE/timing.txt
  cat $out.$type.tsf  | gawk -F '\t' '/^#/{next;}{n=split($4,aa,":");if(n==1)z=$4;if(n==2)z=60*aa[1]+aa[2];if(n==3)z=3600*aa[1]+60*aa[2]+aa[3];z=z/60.0;if(z+0>=0.002)printf("%s\t%s\tt\t%.2f\n",$1,$2,z);}' | bin/tsf -I tsf -O table --title $type.$SV >> COMPARE/timing.txt
  echo "" >> COMPARE/timing.txt

  set type=wallT
  set title="Overall elapsed time in hours:minutes:seconds on the same hardware, as reported by the unix command /usr/bin/time %E"
  echo "### $title" >>  COMPARE/timing.txt
  cat $out.$type.tsf  | gawk -F '\t' '/^#/{next;}{if (! $4)next;printf("%s\t%s\tt\t",$1,$2);n=split($4,aa,":");if(n==1)printf("%d",int($4));if (n==2)printf("%d:%02d",aa[1],int(aa[2]));if (n==3)printf("%d:%02d:%02d",aa[1],aa[2],int(aa[3]));printf("\n");}' | bin/tsf -I tsf -O table --title $type.$SV >> COMPARE/timing.txt
  echo "" >> COMPARE/timing.txt

  set type=maxMem
  set title="Peak RAM in Gigabytes, as reported by the unix command /usr/bin/time %M"
  echo "### $title" >>  COMPARE/timing.txt
  cat $out.$type.tsf  | gawk -F '\t' '/^#/{next;}{z=$4/1000000;if(z+0>=0.001)printf("%s\t%s\tt\t%.2f\n",$1,$2,z);}' | bin/tsf -I tsf -O table --title $type.$SV >> COMPARE/timing.txt
  echo "" >> COMPARE/timing.txt
  
  set type=parallel
  set title="Average number of running threads, as reported by the unix command /usr/bin/time %P"
  echo "### $title" >>  COMPARE/timing.txt
  cat $out.$type.tsf  | gawk -F '\t' '/^#/{next;}{z=$4;gsub("%","",z);z=int((z+50)/100);if(z+0>=1)printf("%s\t%s\tt\t%s\n",$1,$2,z);}' | bin/tsf -I tsf -O table --title $type.$SV >> COMPARE/timing.txt
  echo "" >> COMPARE/timing.txt
  	
  set type=cpu
  set title="CPU time in decimal minutes, sum of the system and user times, as reported by the unix command /usr/bin/time %U + %S"
  echo "### $title" >>  COMPARE/timing.txt
  cat $out.$type.tsf  | gawk -F '\t' '/^#/{next;}{z=$4/60.0;if($4+0>=10)printf("%s\t%s\tt\t%.2f\n",$1,$2,z);}' | bin/tsf -I tsf -O table --title $type.$SV >> COMPARE/timing.txt
  echo "" >> COMPARE/timing.txt
  echo "\n" >> COMPARE/timing.txt
end
  mv COMPARE/timing.txt COMPARE/timing.$SV.txt
    
goto phaseLoop

# list the unaligned RefSeq
  cat tatouI.*.hits | gawk '/^#/{next;}{split($1,aa,"/");print aa[2];}' | sort -u > _a
  zcat Fasta/iRefSeq/iRefSeq.fasta.gz | gawk '/^>/{print substr($1,2);}' | sort > _b
  cat _a _b | sort | gawk '{nn[$1]++;}END{for (k in nn)if(nn[k]==1)print k}' | sort | head


phase_aliqc:
echo -n "phase_aliqc :"
date

 foreach run ($runs)
  set nRawReads=`cat Fasta/$run/$run*.count | gawk '/^Sequence_kept/{n+=$2;}END{print n+0}'`
  set nRawBases=`cat Fasta/$run/$run*.count | gawk '/^Bases_seq_kept/{n+=$2;}END{print n+0}'`
  foreach mm ($allMethods)
    if (! -e RESULTS/$mm/$run/align.done) then
      if (-e RESULTS/$mm/$run/s2g.samStats) \rm RESULTS/$mm/$run/s2g.samStats
      continue
    endif
    if (-e RESULTS/$mm/$run/Aligned.out.sam && ! -e RESULTS/$mm/$run/sam) then
       mv RESULTS/$mm/$run/Aligned.out.sam RESULTS/$mm/$run/sam
    endif

    if (-e RESULTS/$mm/$run/$run/$run.sam && ! -e RESULTS/$mm/$run/s2g.samTools.txt) then
      samtools stats RESULTS/$mm/$run/$run/$run.sam | gawk -f scripts/sam_stats.awk > RESULTS/$mm/$run/s2g.samTools.txt  
    endif

    if (! -e RESULTS/$mm/$run/s2g.samStats) then
      set sam=RESULTS/$mm/$run/sam
      if (! -e $sam) set sam=RESULTS/$mm/$run/$run.sam
      if (! -e $sam) set sam=RESULTS/$mm/$run/Aligned.out.sam
      if (! -e $sam) set sam=RESULTS/$mm/$run/$mm.$run.sam
      if (! -e $sam) set sam=RESULTS/$mm/$run/sam.gz
      if (! -e $sam) set sam=RESULTS/$mm/$run/sam_sorted
      if (! -e $sam) set sam=RESULTS/$mm/$run/sam_sorted.gz
      if (! -e $sam) continue

      set exportitr=""
      if ($run == iRefSeq || $run == iRefSeq38) set exportitr="--exportIntronSupport"

      set exportitr="--exportIntronSupport"

      if (! -e RESULTS/$mm/$run/s2g.samTools.txt) then
        samtools stats $sam | gawk -f scripts/sam_stats.awk > RESULTS/$mm/$run/s2g.samTools.txt &
      endif
      touch RESULTS/$mm/$run/s2g.samStats
                                                echo "bin/sam2gold --method $mm --run $run --samStats --nRawBases $nRawBases --nRawReads $nRawReads -i $sam -o RESULTS/$mm/$run/s2g --addReadPairSuffixForce $exportitr"
                                                echo "bin/sam2gold --method $mm --run $run --samStats --nRawBases $nRawBases --nRawReads $nRawReads -i $sam -o RESULTS/$mm/$run/s2g --addReadPairSuffixForce $exportitr" > 	RESULTS/$mm/$run/s2g.cmd
            scripts/submit RESULTS/$mm/$run/s2g "bin/sam2gold --method $mm --run $run --samStats --nRawBases $nRawBases --nRawReads $nRawReads -i $sam -o RESULTS/$mm/$run/s2g --addReadPairSuffixForce $exportitr" 64G
	    # 64G
    
    endif
    # wc RESULTS/$mm/$run/s2g.samStats
  end
  foreach mm ($allMethods)
    if (! -e RESULTS/$mm/$run/align.done) continue
    if (-e RESULTS/$mm/$run/$run.s2g.samStats) then
      \cp RESULTS/$mm/$run/$run.s2g.samStats RESULTS/$mm/$run/s2g.samStats
    endif
    if (-e RESULTS/$mm/$run/$run/s2g.samStats) then
      \cp RESULTS/$mm/$run/$run/s2g.samStats RESULTS/$mm/$run/s2g.samStats
    endif
  end
end

echo -n "### Quality control for all methods and datasets : $SV : "  > COMPARE/samStats.$SV.txt
date  >> COMPARE/samStats.$SV.txt
echo "### True error rates in Baruzzo datasets:   t1=0.543,  t2=1.186, t3=6.024" >> COMPARE/samStats.$SV.txt
cat RESULTS/*/*/s2g.samStats | sed -e 's/nMultiAligned 0 times/nUnaligned/g' -e 's/nMultiAligned 1 times/nAlignedOnce/g' -e 's/nMultiAligned 2 times/nMultiAligned_2_sites/g' > RESULTS/allSamStats

\mv RESULTS/allSamStats RESULTS/allSamStats.old
cat RESULTS/allSamStats.old | sed -e 's/nAlignedReads/AlignedReads/' -e 's/nPerfectReads/Perfect_reads/' -e 's/nRawBases/RawBases/' -e 's/nRawReads/RawReads/' -e 's/nReads/Reads/' -e 's/nErrors/nMismatches_and_InDels/' -e 's/nMismaches_and_InDels/nMismatches_and_InDels/' > RESULTS/allSamStats

mv RESULTS/allSamStats RESULTS/allSamStats.1
foreach run ($runs)
    set nRun=`echo "$run ZZZZZ $runs" | gawk '{r=$1;for (i=3;i<=NF;i++)if(r==$i)k=i-2;}END{printf ("%2d", k);}'`
    cat RESULTS/allSamStats.1 | gawk -F '\t' '{if($1 == run){m=$2;if(substr(m,1,2)=="01")m=m"."SV;printf("%2d_%s\t%s",k,run,m);for(i=3;i<=NF;i++)printf("\t%s",$i);printf("\n");}}' k=$nRun run=$run SV=$SV>> RESULTS/allSamStats
end

setenv runsN `echo "$runs" | gawk '{for(k=1;k<=NF;k++)printf("%2d_%s ",k,$k);}'`

set tsfMethods=`echo $methods | gawk '{sep="";for(i=1;i<=NF;i++){printf("%s%s",sep,$i);if(substr($i,1,2)=="01")printf(".%s",SV);sep=",";}}END{printf("\n");}' SV=$SV`
echo $tsfMethods

foreach tag (AlignedReads nAlignedBases nMismatches_and_InDels Perfect_reads nUnaligned nAlignedOnce nMultiAligned)
  echo "\n$tag\t$SV" >> COMPARE/samStats.$SV.txt
  if (-e  toto.tag) \rm toto.tag
  cat RESULTS/allSamStats | gawk -F '\t' '{gsub (" ", "_",$3);if (length($5) >= 1 && $3 == tag) {printf("%s\t%s\tt\t%s\n", $1,$2,$5);}}' tag=$tag >> toto.tag
  foreach run ($runsN)
    foreach mm ($allMethods)
      echo "$run\t$mm\tt\t0" >> toto.tag
    end
  end
  cat toto.tag | bin/tsf --sampleSelect $tsfMethods  -I tsf -O table --title perCent.$tag >> COMPARE/samStats.$SV.txt
  echo "\n" >> COMPARE/samStats.$SV.txt
end

foreach tag (AlignedReads nAlignedBases nMismatches_and_InDels Perfect_reads RawBases RawReads nUnaligned nAlignedOnce nMultiAligned)
  echo $tag
  if (-e  toto.tag) \rm toto.tag
  echo "\n$tag\t$SV" >> COMPARE/samStats.$SV.txt
  foreach run ($runsN)
    foreach mm ($allMethods)
      echo "$run\t$mm\ti\t0" >> toto.tag
    end
  end
  cat RESULTS/allSamStats | gawk -F '\t' '{gsub (" ", "_",$3);if (length($4) >= 1 && $3 == tag) {printf("%s\t%s\ti\t%s\n", $1,$2,$4);}}' tag=$tag >> toto.tag  
  cat toto.tag | bin/tsf --sampleSelect $tsfMethods    -I tsf -O table --title $tag >> COMPARE/samStats.$SV.txt
  echo "\n" >> COMPARE/samStats.$SV.txt
end


foreach tag (nAlignments)
  echo $tag
  if (-e  toto.tag) \rm toto.tag
  echo "\n$tag\t$SV" >> COMPARE/samStats.$SV.txt
  foreach run ($runsN)
    foreach mm ($allMethods)
      echo "$run\t$mm\tf\t0" >> toto.tag
    end
  end
  echo "\nAverage_number_of_alignments per aligned read\t$SV" >> COMPARE/samStats.$SV.txt
  cat RESULTS/allSamStats | gawk -F '\t' '{gsub (" ", "_",$3);if (length($5) >= 1 && $3 == tag) {gsub("%","",$5);printf("%s\t%s\tf\t%s\n", $1,$2,$5);}}' tag=$tag >> toto.tag
   cat toto.tag | bin/tsf --sampleSelect $tsfMethods    -I tsf -O table --title "Average number of alignments" >> COMPARE/samStats.$SV.txt
  echo "\n" >> COMPARE/samStats.$SV.txt
end


foreach tag (nMultiAligned)
  echo $tag
  if (-e  toto.tag) \rm toto.tag
  echo "\n$tag\t$SV" >> COMPARE/samStats.$SV.txt
  foreach run ($runsN)
    foreach mm ($allMethods)
      echo "$run\t$mm\tf\t0" >> toto.tag
    end
  end
  echo "\nPercentaged of aligned read which are multi aligned\t$SV" >> COMPARE/samStats.$SV.txt
  cat RESULTS/allSamStats | gawk -F '\t' '{gsub (" ", "_",$3);if (length($5) >= 1 && $3 == tag) {gsub("%","",$5);printf("%s\t%s\tf\t%s\n", $1,$2,$5);}}' tag=$tag >> toto.tag
   cat toto.tag | bin/tsf --sampleSelect $tsfMethods    -I tsf -O table --title "Average number of alignments" >> COMPARE/samStats.$SV.txt
  echo "\n" >> COMPARE/samStats.$SV.txt
end

foreach tag (Supported_introns)
  echo $tag
  if (-e  toto.tag) \rm toto.tag
  echo "\n$tag\t$SV" >> COMPARE/samStats.$SV.txt
  foreach run ($runsN)
    foreach mm ($allMethods)
      echo "$run\t$mm\ti\t0" >> toto.tag
    end
  end
  echo "\nSupported introns\t$SV" >> COMPARE/samStats.$SV.txt
  cat RESULTS/allSamStats | gawk -F '\t' '{gsub (" ", "_",$3);if (length($4) >= 1 && $3 == tag) {gsub("%","",$5);printf("%s\t%s\ti\t%s\n", $1,$2,$4);}}' tag=$tag >> toto.tag
   cat toto.tag | bin/tsf --sampleSelect $tsfMethods    -I tsf -O table --title "Supported introns" >> COMPARE/samStats.$SV.txt
  echo "\n" >> COMPARE/samStats.$SV.txt
end


foreach tag (Intron_supports)
  echo $tag
  if (-e  toto.tag) \rm toto.tag
  echo "\n$tag\t$SV" >> COMPARE/samStats.$SV.txt
  foreach run ($runsN)
    foreach mm ($allMethods)
      echo "$run\t$mm\ti\t0" >> toto.tag
    end
  end
  echo "\nIntron supports\t$SV" >> COMPARE/samStats.$SV.txt
  cat RESULTS/allSamStats | gawk -F '\t' '{gsub (" ", "_",$3);if (length($4) >= 1 && $3 == tag) {gsub("%","",$5);printf("%s\t%s\ti\t%s\n", $1,$2,$4);}}' tag=$tag >> toto.tag
   cat toto.tag | bin/tsf --sampleSelect $tsfMethods    -I tsf -O table --title "Intron supports" >> COMPARE/samStats.$SV.txt
  echo "\n" >> COMPARE/samStats.$SV.txt
end
foreach tag (Intron_supports)
  echo $tag
  if (-e  toto.tag) \rm toto.tag
  echo "\n$tag\t$SV" >> COMPARE/samStats.$SV.txt
  foreach run ($runsN)
    foreach mm ($allMethods)
      echo "$run\t$mm\tf\t0" >> toto.tag
    end
  end
  echo "\nIntron stranding\t$SV" >> COMPARE/samStats.$SV.txt
  cat RESULTS/allSamStats | gawk -F '\t' '{gsub (" ", "_",$3);if (length($7) >= 1 && $3 == tag) {gsub("%","",$7);printf("%s\t%s\tf\t%s\n", $1,$2,$7);}}' tag=$tag >> toto.tag
   cat toto.tag | bin/tsf --sampleSelect $tsfMethods    -I tsf -O table --title "Intron strandings" >> COMPARE/samStats.$SV.txt
  echo "\n" >> COMPARE/samStats.$SV.txt
end

foreach tag (Non_compatible_pairs)
  echo $tag
  if (-e  toto.tag) \rm toto.tag
  echo "\n$tag\t$SV" >> COMPARE/samStats.$SV.txt
  foreach run ($runsN)
    foreach mm ($allMethods)
      echo "$run\t$mm\tf\t0" >> toto.tag
    end
  end
  echo "\nNon compatible pairs\t$SV" >> COMPARE/samStats.$SV.txt
  cat RESULTS/allSamStats | gawk -F '\t' '{gsub (" ", "_",$3);if (length($5) >= 1 && $3 == tag) {gsub("%","",$5);printf("%s\t%s\tf\t%s\n", $1,$2,$5);}}' tag=$tag >> toto.tag
   cat toto.tag | bin/tsf --sampleSelect $tsfMethods    -I tsf -O table --title "Compatible pairs" >> COMPARE/samStats.$SV.txt
  echo "\n" >> COMPARE/samStats.$SV.txt
end

foreach tag (Non_compatible_pairs)
  echo $tag
  if (-e  toto.tag) \rm toto.tag
  echo "\n$tag\t$SV" >> COMPARE/samStats.$SV.txt
  foreach run ($runsN)
    foreach mm ($allMethods)
      echo "$run\t$mm\tf\t0" >> toto.tag
    end
  end
  echo "\nIncompatible pairs\t$SV" >> COMPARE/samStats.$SV.txt
  cat RESULTS/allSamStats | gawk -F '\t' '{gsub (" ", "_",$3);if (length($5) >= 1 && $3 == tag) {gsub("%","",$5);printf("%s\t%s\tf\t%s\n", $1,$2,$5);}}' tag=$tag >> toto.tag
   cat toto.tag | bin/tsf --sampleSelect $tsfMethods    -I tsf -O table --title "Incompatible pairs" >> COMPARE/samStats.$SV.txt
  echo "\n" >> COMPARE/samStats.$SV.txt
end

\cp bin/sortalign COMPARE/sortalign.$SV.exe

goto phaseLoop

phase_Accuracy:
echo -n "phase_Accuracy :"
date

foreach run ($runs)
  if (-e Fasta/iRefSeq/iRefSeq.cig.gz && ! -e Fasta/$run/$run.cig.gz) then
    pushd Fasta/$run
      ln -s ../iRefSeq/iRefSeq.cig.gz   $run.cig.gz
    popd
  endif
end

# Illumina $HG19_runs  $PFAL_runs 
#$main_runs $pacbio_runs $long_illumina_runs  $HG19_runs  $PFAL_runs $methods
foreach run ($runs)
  foreach mm ($methods)
     if (-e Fasta/$run/$run.cig.gz && -e RESULTS/$mm/$run/sam && ! -e RESULTS/$mm/$run/sam_sorted.gz) then
       # sorting can be very slow
       echo -n "\n sort the sam file RESULTS/$mm/$run/sam : "
       date
       time sort RESULTS/$mm/$run/sam > RESULTS/$mm/$run/sam_sorted
       gzip RESULTS/$mm/$run/sam_sorted
       date
       ls -ls RESULTS/$mm/$run/sam_sorted.gz
       # \rm $out/sam
  endif

     endif
     if (-e Fasta/$run/$run.cig.gz && -e RESULTS/$mm/$run/sam_sorted.gz && ! -e RESULTS/$mm/$run/sam2gold.out) then
       echo RESULTS/$mm/$run/sam_sorted.gz
       set arp=""
       if (-e Fasta/$run/Paired_end) then
         set arp=`echo $mm | gawk 'BEGIN{arp="";}{if(index($1,"STAR")>0) arp="-addReadPairSuffix"; if(index($1,"60_")>0) arp="-addReadPairSuffix2"; if(index($1,"70_")>0) arp="-addReadPairSuffixForce";}END{printf("%s",arp);}'`
       endif
       \rm  RESULTS/$mm/$run/sam2gold.*
       set snpF=""
       if (-e Fasta/$run/$run.snps.gz) set snpF='-s '$run..$mm':'Fasta/$run/$run.snps.gz
       echo                                      "bin/sam2gold $arp $snpF -g $run..GOLD:Fasta/$run/$run.cig.gz -i  $run..$mm':'RESULTS/$mm/$run/sam_sorted.gz -r $run --method $mm -o  RESULTS/$mm/$run/$run > RESULTS/$mm/$run/sam2gold.out"
       scripts/submit $mm/$run/$mm.$run.sam2gold "bin/sam2gold $arp $snpF -g $run..GOLD:Fasta/$run/$run.cig.gz -i  $run..$mm':'RESULTS/$mm/$run/sam_sorted.gz -r $run --method $mm -o  RESULTS/$mm/$run/$run > RESULTS/$mm/$run/sam2gold.out" 64G
     endif
  end
end
 
goto phaseLoop

##############################################################################
##############################################################################
## Alignment quality control
## Evaluate in great details the intrinsic quality of the alignment results, provided in BAM format 
## This analysis does not refer to the gold truth
## This is a python.2.7 scripts given in scripts/AliQC.py
## It was developped in collaboration with Joe Meehan, FDA, for the MAQC/SEQC project
## There is a dependency, one must fisrt install HTSeq as explained in the previous section
##
## aliqc produces a computer friendly tag-values tab delimited format .aliqc.tsv
## In the following section aliqc is used again to merge these file into a single table

# OBSOLETE use phase_aliqc: sam2gold --samstats
phase_aliqc_py:

# create sam_sorted only once
set ok=1
  foreach mm ($methods)
    foreach run ($runs)
      if (-e RESULTS/$mm/$run/bam && ! -e RESULTS/$mm/$run/sam_sorted.gz) then
        set ok=0
        echo "transformng RESULTS/$mm/$run/bam into sam_sorted.gz"
        scripts/submit RESULTS/$mm/$run/samview "samtools view RESULTS/$mm/$run/bam | sort -T $TMPDIR | gzip >  RESULTS/$mm/$run/sam_sorted.gz" 
      endif
    end
  end
if ($ok == 0) goto phaseLoop

# use sam_sorted.gz rather than bam (aliqc --BAM also works, but it need to call sort which is very costly)
# PacBio Roche iRefSeq SRR5189652 SRR5189667 HG19t1r1 HG19t1r2 HG19t2r1 HG19t2r2 HG19t3r1 PFALt1r1  PFALt1r2  PFALt1r3  PFALt2r1 PFALt2r2 PFALt2r3 PFALt3r1 PFALt3r2 PFALt3r3 SRR5438850
foreach minAli (0 50 80)
  foreach run ($runs)
    foreach mm ($methods)
      if (-e RESULTS/$mm/$run/sam_sorted.gz && -e Fasta/$run/genome.gz && ! -e RESULTS/$mm/$run/minAliv2.$minAli.aliqc.tsv) then
        echo "Running AliQC.py on RESULTS/$mm/$run/sam_sorted.gz"
	set nreads=`cat Fasta/$run/*.count  | gawk '/^Tags_kept/{n+=$2}END{print n+0}'`
	set nbases=`cat Fasta/$run/*.count  | gawk '/^Bases_tags_kept/{n+=$2}END{print n+0}'`
	if ($nreads == 0) then
	  echo "missing file Fasta/$run/$run.count, please run phase 1"
        else
          echo "Running AliQC.py $mm $run --nreads $nreads --nbases $nbases"
	  if (-e RESULTS/$mm/$run/f2.20000.sam_sorted) then
            set ok=1
            foreach ii (`seq 1 1000`)
              if (! -e   RESULTS/$mm/$run/f2.$ii.sam_sorted) continue
              if (-e  RESULTS/$mm/$run/f2.$ii.min.$minAli.aliqc.tsv) continue
              scripts/submit RESULTS/$mm/$run/f2.$ii.minali$minAli "python3 scripts/AliQC.py --SAMSORTED -i RESULTS/$mm/$run/f2.$ii.sam_sorted -r $run..$mm.minAli$minAli -f Fasta/$run/genome.gz -o RESULTS/$mm/$run/f2.$ii.min.$minAli --minAli $minAli"
              set ok=0
            end
            if ($ok == 1) then
              cat  RESULTS/$mm/$run/f2.*.min.$minAli.aliqc.tsv >  RESULTS/$mm/$run/f2_all.aliqc.tsv
              python3 scripts/AliQC.py --merge -i RESULTS/$mm/$run/f2_all.aliqc.tsv -r $run..$mm.minAli$minAli -o RESULTS/$mm/$run/$mm.$run.minAliv2.$minAli --nreads $nreads --nbases $nbases --minAli $minAli"
            endif
          else
	    set ss=`echo $run | gawk '{r=$1;f="";if(mm=="70_Minimap2" && (r=="Illumina" || substr(r,1,3)=="SRR"))f="F";printf("SAMSORTEDGZ%s",f);}' mm=$mm`
            if (-e RESULTS/$mm/$run/aliqc.minali$minAli.err) \rm RESULTS/$mm/$run/aliqc.minali$minAli.err RESULTS/$mm/$run/aliqc.minali$minAli.out 
            scripts/submit RESULTS/$mm/$run/aliqc.minali$minAli "python3 scripts/AliQC.py --$ss -i RESULTS/$mm/$run/sam_sorted.gz -r $run..$mm.minAli$minAli -f Fasta/$run/genome.gz -o RESULTS/$mm/$run/minAliv2.$minAli --nreads $nreads --nbases $nbases --minAli $minAli"
          endif
         endif
      endif
    end
  end
end

goto phaseLoop

##############################################################################
##############################################################################
## Direct statistics of the error counts reported in the BAM files
## The sam2gold.c code, above, parses the bam file and the genome
## and computes its own evaluation of the number of error per aligned read
## In the present script, we rely on the presence in the BAM files of the NM:i:x
## optional field, collate the X values and report the statistics
## Hopefully, the 2 methods should be compatible, but they do not have
## to agree exactly, now aliqc counts a double or triple deletion as a 2 or 3 mismatches
## and some aligners may report it as 1, 2 or 3 errors
## and there is a big issue with deltions versus introns

phase_DirectErrorCount:
echo phase_DirectErrorCount
foreach mm ($methods)
  foreach run ($runs)
    echo "phase_DirectErrorCount $mm $run"
    if (-e RESULTS/$mm/$run/sam && ! -e RESULTS/$mm/$run/direct_eror_count) then
      scripts/submit RESULTS/$mm/$run/$mm.$run.nerrors "scripts/directErrorCount.tcsh $mm $run" 
     endif
  end
end

# export the results in a single human readable table
set toto=RESULTS/Error_counts_using_NM_optional_tag.txt
echo -n "## $toto :" > $toto
if (-e $toto.1) \rm $toto.1
date >> $toto
foreach mm ($methods)
  foreach run ($runs)
    if (-e  RESULTS/$mm/$run/$mm.$run.nerrors) then
      echo -n "$mm\t$run\t" >> $toto.1
      cat RESULTS/$mm/$run/$mm.$run.nerrors | gawk '/^#/{next;}{n[$1]=$2;if($1>max)max=$1;}END{for(i=-5;i<=max;i++)printf("\t%d",n[i]);printf("\n")}' >> $toto.1
    endif
  end
end
cat $toto.1 | gawk '{if(NF>max)max=NF;}END{printf("\t\t");max-=6;for(i=-5;i<=max;i++)printf("\t%d",i);printf("\n")}' > $toto.2
cat $toto.2 $toto.1 | gawk -F '\t' -f scripts/transpose.awk >> $toto
\rm $toto.1 $toto.2
echo "The table of errors using the optional NM tag of the BAM files is in"
echo $toto

goto phaseLoop

##############################################################################
##############################################################################
## Count the substitutions as declared in the Baruzzo Benchmark
## The statistics only measures the substitutions in the r3 reads
## fully and continuously aligned on the plus strand of chromosome 8
## This seems sufficient since it involves in each case at least 100,000 reads 

phase_count_subtitutions_in_benchmark:

# there are several phases in the calculation
# 1: select the full reads, characterized by a cigar string 100M
if (! -e SUBS) mkdir SUBS
#HG19 PFAL
foreach sp (HG19)
  foreach tt (t1 t2 t3)
    if (-e Fasta/$sp$tt'r1'/$sp$tt'r1'.cig.gz && ! -e  SUBS/subs.$sp$tt) then 
      zcat Fasta/$sp$tt'r1'/$sp$tt'r1'.cig.gz | grep chr8  | grep '+' |  grep 100M | cut -f 1,2,3,4,8 > SUBS/subs.$sp$tt
    endif
  end
end

# 2: construct a 6 columns tab delimited shadow file, summarizing the coordinate of the alignemnts
foreach sp (HG19)
  foreach tt (t1 t2 t3)
    if (-e  SUBS/subs.$sp$tt && ! -e SUBS/subs.$sp$tt.shadow ) then 
      cat SUBS/subs.$sp$tt | gawk -F '\t' '{printf("%s\t1\t100\t%s\t%d\t%d\n",$1,$2,$3,$4);}' > SUBS/subs.$sp$tt.shadow 
    endif
  end
end

# 3: isolate the genome of chromosome 8, using the dna2dna utility
foreach sp (HG19)
  if (-e Reference_genome/$sp.Baruzzo.genome.fasta.gz && ! -e Reference_genome/$sp.Baruzzo.chr8.fasta.gz) then
    bin/dna2dna -i Reference_genome/$sp.Baruzzo.genome.fasta.gz -I fasta -O fasta -keepName -o Reference_genome/$sp.Baruzzo.chr8 -gzo 
  endif
end

# 4: Export the corresponding genomic segment in raw format
# The raw format has just 2 tab delimited columns: atgcatgc  identifier
# Notice that dna2dna is very versatile, it can directly export messenger RNAs given a genome and a gff file.
# try bin/dna2dna --help for a full list of functioalities
foreach sp (HG19)
  foreach tt (t1 t2 t3)
    if (-e  SUBS/subs.$sp$tt.shadow && -e Reference_genome/$sp.Baruzzo.chr8.fasta.gz && ! -e SUBS/subs.$sp$tt.raw) then
      bin/dna2dna -i Reference_genome/$sp.Baruzzo.chr8.fasta.gz -shadow SUBS/subs.$sp$tt.shadow -O raw -keepName >  SUBS/subs.$sp$tt.raw 
    endif
  end
end

# 5: the first subs file contains in column 1 and 5 the name and sequence of each read
#    the raw file contains in column 2 and 1 the name and sequence of each corresponding genomic segment
#    Both sequences are exactly 100 bases long, hence a simple awk script is sufficient to count the mismatching bases
echo ZZZZZ > SUBS/ZZZZZ
foreach sp (HG19 PFAL)
  foreach tt (t1 t2 t3)
    cat SUBS/subs.$sp$tt SUBS/ZZZZZ  SUBS/subs.$sp$tt.raw | gawk -F '\t' '/^ZZZZZ/{zz=1;}{if(zz<1){seq[$1]=$5;next;}if (seq[$2]){s1=seq[$2];s2=toupper($1);n=0;for(i=1;i<=100;i++)if(substr(s1,i,1) != substr(s2,i,1))n++;print n}}' | tags | sort -V -k 1n > SUBS/subs.$sp$tt.txt &
  end
end

# 6: produce a final synthetic table
set toto=RESULTS/mm_stats.Baruzzo.txt
echo -n "### $toto : " > $toto
date >> $toto
foreach sp (HG19 PFAL)
  foreach tt (t1 t2 t3)
    cat  SUBS/subs.$sp$tt.txt | gawk '{n1+=$2;n2+=$1*$2;printf ("%s\t%d\t%d\t%d\t%d\n",t,$1,$2,n1,n2);}' t=$sp$tt >> $toto
  end
end
echo "The distribution of substitutions in chromosome 8 in tbe Baruzzo benchmark"
echo "have been exported in $toto"
head -12 $toto

goto phaseLoop

##############################################################################
##############################################################################
## Exportation of global, human readable, quality control tables
## These tables where used directly to prepare the plots and table of the
## Magic-BLAST paper

phase_Export:

if (! -d RESULTS) mkdir RESULTS


# collate the aliqc.tsv tables from all runs using again the AliQC.py scripts with --table option
# this will produce 3 output tables
foreach minAli (0 50 80)
  if (-e toto) \rm toto
  foreach method ($methods)
    foreach run ($runs)
      if (-e RESULTS/$method/$run/$method.$run.minAliv2.$minAli.aliqc.tsv) then
        cat RESULTS/$method/$run/$method.$run.minAliv2.$minAli.aliqc.tsv  >> toto
      else
        echo "Reference_File:Fasta/$run/genome.gz\t$run..$method.minAli$minAli\t0" >> toto
      endif
    end
  end
  cat toto | python3 scripts/AliQC.py --view table --split --minAli $minAli -o RESULTS/aliqc.minAli$minAli
end

# reformat the 3 output tables
foreach minAli (0 50 80)
  foreach type (mismatch_histo_and_types mismatch_per_cycle_per_kb_aligned aligned_reads_multiplicity_length_aligned)
    set toto=RESULTS/aliqc.minAli$minAli.$type
    cat $toto.tsv  | head -12 |  gawk -F '\t' '/^##/{next;}/^#/{print}' > $toto.txt
    cat $toto.tsv  | gawk -F '\t' '/^###/{next;}/^##/{print}' >> $toto.txt

    echo "Any\nR1\nR2\nZZZZZ" > $toto.tsv1
    foreach method ($methods)
      echo $method >>  $toto.tsv1
    end
    echo ZZZZZ >> $toto.tsv1
    foreach run ($runs)
      echo $run >>  $toto.tsv1
    end
    echo ZZZZZ >> $toto.tsv1
    cat $toto.tsv  | gawk -F '\t' '/^#/{next;}{print}' >> $toto.tsv1

    cat $toto.tsv1 | sed -e 's/32_STAR.2.6c/32_STAR_2/g' -e 's/32_STAR.2/32_STAR_2/g' -e 's/14_MagicBLAST_1.5.0/14_MagicBLAST_1_5/g'  -e 's/14_MagicBLAST_1.5/14_MagicBLAST_1_5/g' | gawk -F '\t' '/ZZZZZ/{zz++;next;}{if(zz+0==0){tt[na+0]=$1;na++;next;}}{if(zz+0==1){mm[nm+0]=$1;nm++;next;}}{if(zz+0==2){run[nr+0]=$1;nr++;next;}}{if(NF<5)next;uu=$5;for(i=6;i<=NF;i++)uu=uu "\t" $i; split($3,aa,".");split($4,bb,".");k= $2 bb[1] aa[1] ; z[k]=uu; if(0)printf( "TTT\t#%s#\t%s\n",k,z[k]);next;}END{if(0)for(k in z)print "XXX",k,z[k];for (ia = 0 ; ia < na ; ia++) {for (ir = 0 ; ir < nr ; ir++){for (im = 0 ; im < nm ; im++){k=tt[ia] mm[im] run[ir];printf("\t%s\t%s\t%s\t%s\n",tt[ia],mm[im],run[ir],z[k]);} printf("\n");}printf("\n\n");}if(0)for(k in z)printf("YYY\tAA%sBB\t%s\n",k,z[k]);}' > $toto.tsv2
    cat $toto.tsv2 | grep -v r2 | grep -v r3 | sed -e 's/r1//g' >> $toto.txt

    echo "\n\n\n\n\n\n" >>  $toto.txt
    cat $toto.tsv1 | sed -e 's/32_STAR.2.6c/32_STAR_2/g' -e 's/32_STAR.2/32_STAR_2/g' -e 's/14_MagicBLAST_1.5.0/14_MagicBLAST_1_5/g'  -e 's/14_MagicBLAST_1.5/14_MagicBLAST_1_5/g' | gawk -F '\t' '/ZZZZZ/{zz++;next;}{if(zz+0==0){tt[na+0]=$1;na++;next;}}{if(zz+0==1){mm[nm+0]=$1;nm++;next;}}{if(zz+0==2){run[nr+0]=$1;nr++;next;}}{if(NF<5)next;uu=$5;for(i=6;i<=NF;i++)uu=uu "\t" $i; split($3,aa,".");split($4,bb,".");k= $2 bb[1] aa[1] ; z[k]=uu; if(0)printf( "TTT\t#%s#\t%s\n",k,z[k]);next;}END{if(0)for(k in z)print "XXX",k,z[k];for (ia = 0 ; ia < na ; ia++) {for (im = 0 ; im < nm ; im++){for (ir = 0 ; ir < nr ; ir++){k=tt[ia] mm[im] run[ir];printf("\t%s\t%s\t%s\t%s\n",tt[ia],mm[im],run[ir],z[k]);} printf("\n");}printf("\n\n");}if(0)for(k in z)printf("YYY\tAA%sBB\t%s\n",k,z[k]);}' > $toto.tsv2
    cat $toto.tsv2 | grep -v r2 | grep -v r3 | sed -e 's/r1//g' >> $toto.txt
    \rm  $toto.tsv1 $toto.tsv2

    cat $toto.tsv  | gawk -F '\t' '/^###/{print}' >> $toto.txt
  end
end

set toto=RESULTS/Mapping_accuracy.txt
echo -n "## $toto :" > $toto
date >> $toto 
cat scripts/mapping_accuracy.header   >> $toto
 
if (-e $toto.1) \rm $toto.1
foreach run ($runs)
  foreach method ($methods)
      if (-e $method/$run/$method.$run.introns.tsv) then
        cat $method/$run/$method.$run.introns.tsv | sed  -e 's/\.\./\t/' | grep GoldMap | grep -v GOLD | sort -V -k 2,2 -k 3,3 > $toto.1
      else
        echo "Intron_1\t$run\tGOLD\t8\t210509\t210509\t210509\t0\t0\t1.0000\t1.0000\t1.0000"  > $toto.1
      endif
      cat $toto.1 |  gawk -F '\t' -f scripts/introns_precision_recall.awk  | sed  -e 's/t1r/T1R/g'  -e 's/t2r/T2R/g'  -e 's/t3r/T3R/g' | gawk -F '\t' '/^#/{next;}{print;}'  | sed -e 's/PFAL/Malaria/g' -e 's/STARlongzz/STAR long/g' >> $toto
  end
  echo >> $toto
end

echo "\n\n\n\n\n" >> $toto

if (-e $toto.1) \rm $toto.1
foreach method ($methods)
  foreach run ($runs)
      if (-e $method/$run/$method.$run.introns.tsv) then
        cat $method/$run/$method.$run.introns.tsv | sed  -e 's/\.\./\t/' | grep GoldMap | grep -v GOLD | sort -V -k 2,2 -k 3,3 > $toto.1
      else
        echo "Intron_1\t$run\tGOLD\t8\t210509\t210509\t210509\t0\t0\t1.0000\t1.0000\t1.0000"  > $toto.1
      endif
      cat $toto.1 |  gawk -F '\t' -f scripts/introns_precision_recall.awk  | sed  -e 's/t1r/T1R/g'  -e 's/t2r/T2R/g'  -e 's/t3r/T3R/g' | gawk -F '\t' '/^#/{next;}{print;}'  | sed -e 's/PFAL/Malaria/g' -e 's/STARlongzz/STAR long/g' >> $toto
  end
  echo >> $toto
end

# set toto2=RESULTS/Mapping_accuracy.light.txt
# echo -n "## $toto2 :" > $toto2
# date >> $toto2
# cat scripts/mapping_accuracy.header   >> $toto2

# cat  $toto.1 | grep -v R2 | grep -v R3 |  sed  -e 's/R1//g' | grep -v SRR | grep -v Simulated | sed  -e 's/\.\./\t/' | grep GoldMap | grep -v GOLD | sort -k 2,2 -k 3,3 | gawk -F '\t' -f scripts/introns_precision_recall.awk  | sort -V -u | sed  -e 's/r1//g'   | grep -v r2 | grep -v r3 | gawk -F '\t' '/^#/{print;next;}{if($2!=old)printf("\n");old=$2;print;}' >> $toto2

\rm $toto.1

# goto phaseLoop
##### analysis of the eRefSeq truth

set toto=RESULTS/Truth.snps_stats.txt
echo -n "### $toto : " > $toto
date >> $toto

echo "# Run\tMatch/kb\tMismatch/kb\tSubstitutions/kb\tTransition/kb\tTransversion/kb\tInsertion/kb\tDeletion/kb\tA>G\tT>C\tG>A\tC>T\tA>T\tT>A\tG>C\tC>G\tA>C\tT>G\tG>T\tC>A\tInsA\tInsT\tInsG\tInsC\tDelA\tDelT\tDelG\tDelC\tSimple substitutions\tDouble substitutions\tTriple substitutions\t\t\t\t\t\t\t\t\t\t\tSimple insertions\tDouble insertions\tTriple insertions\t\t\t\t\t\t\t\t\t\t\tSimple deletions\tDouble deletions\tTriple deletions\t\t\t\t\t\t\t\t\t\t" >> $toto
foreach run ($runs)
  if (-e Fasta/$run/$run.snps.gz && ! -e Fasta/$run/$run.spns_stats) then
    echo $run
    set nb=`cat Fasta/$run/$run.count | gawk '/^Bases_tags_kept/{print $2;}'`
    zcat Fasta/$run/$run.snps.gz | cut -f 2,3 | gawk -F '\t' '/>/{nSub++;nMM++;tv++;ntSub[$2]++;x=$1;if(x==x2+2)nkSub[3]++;else if(x==x1+1)nkSub[2]++;else nkSub[1]++;x2=x1;x1=x;}/A>G/{tr++;tv--;}/T>C/{tr++;tv--;}/G>A/{tr++;tv--;}/C>T/{tr++;tv--;}/Ins/{nIns++;k=length($2)-3;nMM+=k;nkIns[k]++;ntIns[substr($2,4)]++;}/Del/{nDel++;k=length($2)-3;nkDel[k]++;ntDel[substr($2,4)]++;}END{print run "\t" 1000.0*(nb-nMM)/nb "\t" 1000.0*nMM/nb "\t" 1000.0*nSub/nb "\t" 1000.0*tr/nb "\t" 1000.0*tv/nb  "\t" 1000.0*nIns/nb "\t" 1000.0*nDel/nb"\t" ntSub["A>G"] "\t" ntSub["T>C"] "\t" ntSub["G>A"] "\t" ntSub["C>T"] "\t" ntSub["A>T"] "\t" ntSub["T>A"] "\t" ntSub["G>C"] "\t" ntSub["C>G"] "\t" ntSub["A>C"] "\t" ntSub["T>G"] "\t" ntSub["G>T"]"\t" ntSub["C>A"] "\t" ntIns["A"]"\t" ntIns["T"]"\t" ntIns["G"]"\t" ntIns["C"]"\t" ntDel["A"]"\t" ntDel["T"]"\t" ntDel["G"]"\t" ntDel["C"]"\t" nkSub[1]+0"\t" nkSub[2]+0"\t" nkSub[3]+0"\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t" nkIns[1]+0"\t" nkIns[2]+0"\t" nkIns[3]+0"\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t" nkDel[1]+0"\t" nkDel[2]+0"\t" nkDel[3]+0"\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";}' nb=$nb run=$run >  Fasta/$run/$run.spns_stats
  endif
  if (-e Fasta/$run/$run.spns_stats) then
    cat  Fasta/$run/$run.spns_stats >> $toto
  endif
end

##### INTRONS report  Insertion Deletion

foreach type (Intron  Insertion Deletion)

  set toto1="RESULTS/$type"_per_coverage_stats.txt
  set toto1L="RESULTS/$type"_per_coverage_stats.light.txt
  set toto3="RESULTS/$type"_per_coverage_stats.best.txt
  set toto4="RESULTS/$type"_per_coverage_stats.1support.txt
  set toto3L="RESULTS/$type"_per_coverage_stats.best.light.txt
  set toto4L="RESULTS/$type"_per_coverage_stats.1support.light.txt
  echo -n "### $toto1 :" > $toto1
  echo -n "### $toto3 :" > $toto3
  echo -n "### $toto4 :" > $toto4
  echo -n "### $toto1L :" > $toto1L
  echo -n "### $toto3L :" > $toto3L
  echo -n "### $toto4L :" > $toto4L
  date >> $toto1
  date >> $toto3
  date >> $toto4
  date >> $toto1L
  date >> $toto3L
  date >> $toto4L

 if ($type == Intron) then
    echo "## An alignment supporting an intron is defined by a line in the SAM/BAM file where the CIGAR contains an N with minimal accepted intron length 50 bases" > $toto1.1
    echo "## When a read is aligned at multiple sites, each of its alignments supporting an intron is counted" >> $toto1.1
    echo "## Some spliced genes are truly repeated, some are very similar. If one rejects all multiply aligned reads, the introns of these genes cannot be detected," >> $toto1.1
    echo "## Therefore, we keep the introns detected by multiply aligned reads, and do not artificially overestimate the specificity of methods unable to select the true positions" >> $toto1.1
    echo "## Note that in the benchmark, all reads are uniquely aligned, yet some support multiple neighboring introns" >> $toto1.1
 else
    echo "## An alignment supporting an indel is defined by a line in the SAM/BAM file where the CIGAR contains an I or o D" > $toto1.1
    echo "## When a read is aligned at multiple sites, each of its alignments supporting an indel is counted" >> $toto1.1
  endif

  echo -n "# Species\tRun\tMethod\tMinimal $type support" >> $toto1.1
  echo "\t$type in benchmark\t$type discovered in method\tFP: False positive $type\tTP: True positive $type\tFN: False negative $type\t$type discovery precision p=TP/(TP+FP)\t$type discovery recall r=TP/(TP+FN)\t$type discovery F-score 2pr/(p+r)" >> $toto1.1

  cat $toto1.1 >> $toto1
  cat $toto1.1 >> $toto1L

  if (-e $toto1.2) \rm $toto1.2

  foreach run ($runs)
    foreach mm ($methods)
      if (-e $toto1.3) \rm $toto1.3
      touch  $toto1.3
      if (-e  RESULTS/$mm/$run/$mm.$run.delins.tsv) then
        cat RESULTS/$mm/$run/$mm.$run.delins.tsv | grep $type | sed -e 's/on_/on/' -e "s/$type//" | grep -v GOLD | sort -V -k 1,1n -k 2,2 |  sed  -e 's/\.\./\t/'  -e 's/g_//' | gawk '/^#/{next;}/GOLD/{next;}{support=$1;species="Human";run=$2;if(0 && support > 200 && run != "Illumina")next;if(substr(run,1,4)=="PFAL")species="Malaria";method=$3;printf("%s\t%s\t%s\t%d",species,run,method,support); printf("\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", $5, $6,$8,$7,$9,$10,$11,$12);}'  > $toto1.3
      endif
      set n=`cat $toto1.3 | wc -l`
      if ($n == 0) then
        echo "$run $mm" | gawk '{species="Human";run=$1;if(substr(run,1,4)=="PFAL")species="Malaria";method=$2;printf("%s\t%s\t%s\t%d\t0\t0\t0\t0\t0\t0\t0\t0\n",species,run,method,1);}' >> $toto1.2
      else
        cat  $toto1.3 >>  $toto1.2
      endif
    end
  end

  cat $toto1.2  >  $toto1.3
  cat $toto1.3 | gawk -F '\t' '/^#/{print;next;}{if($3!=old)printf("\n\n");old=$3;print;}'  >> $toto1
  cat $toto1.3 | grep -v r2 | grep -v r3 | grep -v SRR | grep -v Simulated | gawk -F '\t' '/^#/{print;next;}{if($3!=old)printf("\n\n");old=$3;print;}' >> $toto1L

  echo "## Limited to best support-depth" >> $toto3
  cat $toto1.1 >> $toto3
  cat $toto1.3 | gawk -F '\t' 'BEGIN{best=0;}{z=$1 "\t" $2 "\t" $3; if (z != old) {old=z;if(bestScore)print best;bestScore=$12;best=$0;}if($12>bestScore){bestScore=$12;best=$0;}}END{if(bestScore)print best;}' | gawk -F '\t' '/^#/{print;next;}{if($2!=old)printf("\n");old=$2;print;}' >> $toto3
  
  echo "## Limited to best support-depth" >> $toto3L
  cat $toto1.1 >> $toto3L
   cat $toto1.3 | grep -v r2 | grep -v r3 | grep -v SRR | grep -v Simulated | gawk -F '\t' 'BEGIN{best=0;}{z=$1 "\t" $2 "\t" $3; if (z != old) {old=z;if(bestScore)print best;bestScore=$12;best=$0;}if($12>bestScore){bestScore=$12;best=$0;}}END{if(bestScore)print best;}' | gawk -F '\t' '/^#/{print;next;}{if($2!=old)printf("\n");old=$2;print;}' >> $toto3L
  
  echo "## Limited to 1 support" >> $toto4
  cat $toto1.1 >> $toto4
  cat $toto1.3 | gawk -F '\t' '{if($4==1)print}' | gawk -F '\t' '/^#/{print;next;}{if($2!=old)printf("\n");old=$2;print;}' >> $toto4
  echo "\n\n\n\n\n\n" >> $toto4
  foreach mm ($methods)
    cat $toto1.3 | gawk -F '\t' '{if($4==1 && $3 == mm)print}' mm=$mm | gawk -F '\t' '/^#/{print;next;}{print;}' >> $toto4
    echo "\n" >> $toto4
  end

  echo "## Limited to 1 support" >> $toto4L
  cat $toto1.1 >> $toto4L
  cat $toto1.3 | grep -v r2 | grep -v r3 | grep -v SRR | grep -v Simulated | gawk -F '\t' '{if($4==1)print}' | gawk -F '\t' '/^#/{print;next;}{if($2!=old)printf("\n");old=$2;print;}' >> $toto4L


  \rm $toto1.[123]
end


set toto=RESULTS/AliQC.content.txt
cat RESULTS/aliqc.minAli0.aligned_reads_multiplicity_length_aligned.txt | gawk -F '\t' '{run=$3;m=$4;if(m=="")next;if(run!=oldrun){print z;z=run "\t" m;}else z=z "," m;oldrun=run;}END{print z}' >> $toto


  cat RESULTS/Intron_per_coverage_stats.txt | gawk -F '\t' '/^#/{next}{z=$1 "\t" $2 "\t" $3;k=0+n[z];if($4>k)n[z]=$4;}END{for(k in n)printf("%s\t%d\n",k,n[k]);printf("toto\n");}'  | gawk -F '\t' 'BEGIN{printf("# Species\tRun\tMethods\tMaximal intron support\n")}{z=$1 "\t" $2;if (z != oldz){if(length(oldz)>3)printf("%s\t%s\t%s\n",oldz,substr(m[oldz],2),substr(n[oldz],2));m[z]="";n[z]="";}oldz=z;m[z]=m[z]","$3;n[z]=n[z]","$4;}' > RESULTS/Intron_per_coverage_stats.title.txt
# \rm $toto.*

  cat RESULTS/Mapping_accuracy.txt | gawk -F '\t' '/^#/{next}{n=$5+0;if(n==0)next; z=$1 $2; if(z != oldz)k=0;oldz=z;if(k==0)printf("\n%s\t", z);else printf(",");k++;printf("%s",$3);}END{printf("\n");}'  > RESULTS/Mapping_accuracy.title.txt

  \rm RESULTS/*.tsv

set toto=RESULTS/Mismatch_calling_accuracy.txt
echo -n "### $toto : " > $toto
date >> $toto
echo "## Baratin" >> $toto
echo >> $toto

set ok=0
foreach run ($runs)
  foreach method ($methods)
    if ($ok == 0 && -e $method/$run/$method.$run.snp) then
      cat $method/$run/$method.$run.snp | gawk '/^#/{print}' >> $toto
      set ok=1
    endif
  end
end

foreach run ($runs)
  echo >> $toto
  foreach method ($methods)
    if (-e RESULTS/$method/$run/$method.$run.snp) then
      cat RESULTS/$method/$run/$method.$run.snp | gawk '/^#/{next;}/^$/{next;}{print}' >> $toto
    endif
  end
end

  echo >> $toto
  echo >> $toto
  echo >> $toto
  echo >> $toto

foreach method ($methods)
  echo >> $toto
    foreach run ($runs)
    if (-e RESULTS/$method/$run/$method.$run.snp) then
      cat RESULTS/$method/$run/$method.$run.snp | gawk '/^#/{next;}/^$/{next;}{print}' >> $toto
    endif
  end
end

cat $toto | transpose | cut -f 1-10

goto phaseLoop  # fall through to aliLn

##############################################################################
##############################################################################
####### aliLn

phase_aliLn:


set toto=RESULTS/Aligned_length.histo.txt
echo -n  "## $toto : " > $toto
date >> $toto

\rm $toto.*
  foreach run ($runs)

    set ok=`echo $run | gawk '{ok=1;z=substr($1,length($1-2)) ; if (z == "r2" || z == "r3")ok=0;print ok;}'`
    if ($ok == 0) continue
     
    if (! -e Fasta/$run/$run.TM) then
      if (-e Fasta/$run/$run.fastq.gz) then
        bin/dna2dna -i Fasta/$run/$run.fastq.gz -I fastq -getTM > Fasta/$run/$run.TM
      endif
    endif
    if (! -e Fasta/$run/$run.TM) then
      if (-e Fasta/$run/$run.fasta.gz) then
        bin/dna2dna -i Fasta/$run/$run.fasta.gz -I fasta -getTM > Fasta/$run/$run.TM
      endif
    endif
    if (! -e Fasta/$run/$run.TM) then
      if (-e Fasta/$run/$run'_1'.fasta.gz) then
        zcat  Fasta/$run/$run'_'?.fasta.gz | bin/dna2dna -I fasta -getTM > Fasta/$run/$run.TM
      endif
    endif
    if (! -e Fasta/$run/$run.TM) then
      if (-e Fasta/$run/$run'_R1'.fastq.gz) then
        zcat  Fasta/$run/$run'_'R?.fastq.gz | bin/dna2dna -I fastq -getTM > Fasta/$run/$run.TM
      endif
    endif

    set delta=1
    if ($run == Roche) set delta=10
    if ($run == PacBio) set delta=30
    if ($run == SRR5189652) set delta=30
    if ($run == SRR5189667) set delta=30
    if ($run == SRR5437876) set delta=30
    if ($run == SRR5438850) set delta=10
    if ($run == SRR5438850_120) set delta=1
    if ($run == SRR5438850_150) set delta=1
    if ($run == iRefSeq) set delta=100
    if ($run =~ "RS*") set delta=100
    if ($run == ONG) set delta=100
    if ($run == ONG2) set delta=100
    if ($run == ONG3) set delta=100
    if ($run == ONG4) set delta=100
    if ($run == ONG5) set delta=100

    foreach mm ($methods)
      if (-e  RESULTS/$mm/$run/$mm.$run.aliLn) then
        cat RESULTS/$mm/$run/$mm.$run.aliLn | gawk '{k=int(($1+delta - 1)/delta) ; if(k>900)k=900;nn[k] += $2;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=$delta mm=$mm rr=$run >> $toto.1
      endif  
    end
  end 

  echo "Illumina\ttruth\t101\t217498656" >>  $toto.1
  echo "SRR5437876\ttruth\t300\t32935604" >>  $toto.1
  # echo "SRR5437876_120\ttruth\t120\t32935604" >>  $toto.1
  # echo "SRR5437876_150\ttruth\t150\t32935604" >>  $toto.1
  echo "SRR5438850\ttruth\t250\t46617660" >>  $toto.1 
  echo "SRR5438850_120\ttruth\t120\t46617660" >>  $toto.1 
  echo "SRR5438850_150\ttruth\t150\t46617660" >>  $toto.1 

  foreach run ($runs)
    if ($run == iRefSeq || $run =~ "RS*") then
      cat Fasta/$run/$run.TM | gawk -F '\t' '/^#/{next;}{k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=100 mm="truth" rr=$run > $toto.1a
    endif
    if ($run == PacBio) then
      cat Fasta/PacBio/PacBio.TM | gawk -F '\t' '/^#/{next;}{k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=30 mm="truth" rr=PacBio >> $toto.1a
    endif
    if ($run == Roche) then
      cat Fasta/Roche/Roche.TM | gawk -F '\t' '/^#/{next;}{if($2<19)next;k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=10 mm="truth" rr=Roche >> $toto.1a
    endif
    if ($run ==  SRR5189652 ) then
      cat Fasta/SRR5189652/SRR5189652.TM | gawk -F '\t' '/^#/{next;}{if($2<19)next;k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=30 mm="truth" rr=SRR5189652n >> $toto.1a
    endif
    if ($run == SRR5189667) then
      cat Fasta/SRR5189667/SRR5189667.TM | gawk -F '\t' '/^#/{next;}{if($2<19)next;k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=30 mm="truth" rr=SRR5189667 >> $toto.1a
    endif
    if ($run == ONG) then
      cat Fasta/ONG/ONG.TM | gawk -F '\t' '/^#/{next;}{if($2<19)next;k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=100 mm="truth" rr=ONG >> $toto.1a
    endif
    if ($run == ONG2) then
      cat Fasta/ONG2/ONG2.TM | gawk -F '\t' '/^#/{next;}{if($2<19)next;k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=100 mm="truth" rr=ONG2 >> $toto.1a
    endif
    if ($run == ONG3) then
      cat Fasta/ONG3/ONG3.TM | gawk -F '\t' '/^#/{next;}{if($2<19)next;k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=100 mm="truth" rr=ONG3 >> $toto.1a
    endif
    if ($run == ONG4) then
      cat Fasta/ONG4/ONG4.TM | gawk -F '\t' '/^#/{next;}{if($2<19)next;k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=100 mm="truth" rr=ONG4 >> $toto.1a
    endif
    if ($run == ONG5) then
      cat Fasta/ONG5/ONG5.TM | gawk -F '\t' '/^#/{next;}{if($2<19)next;k=int(($2+delta - 1)/delta) ; if(k>900)k=900;nn[k]++;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=100 mm="truth" rr=ONG5 >> $toto.1a
    endif
  end

  echo "Human_T1\ttruth\t100\t20000000" > $toto.t0
  echo "Human_T2\ttruth\t100\t20000000" >>  $toto.t0
  echo "Human_T3\ttruth\t100\t20000000" >>  $toto.t0
  echo "Malaria_T1\ttruth\t100\t20000000" >>  $toto.t0
  echo "Malaria_T2\ttruth\t100\t20000000" >>  $toto.t0
  echo "Malaria_T3\ttruth\t100\t20000000" >>  $toto.t0
  
  set delta=1
  foreach tt (1 2 3)
    if (-e $toto.t$tt) \rm $toto.t$tt
    foreach mm ($methods)
      cat RESULTS/$mm/HG19t$tt'r1'/$mm.HG19t$tt'r1'.aliLn | gawk '{k=int(($1+delta - 1)/delta) ; if(k>900)k=900;nn[k] += $2;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=$delta mm=$mm rr=Human_T$tt >> $toto.t$tt
      cat RESULTS/$mm/PFALt$tt'r1'/$mm.PFALt$tt'r1'.aliLn | gawk '{k=int(($1+delta - 1)/delta) ; if(k>900)k=900;nn[k] += $2;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%s\t%d\t%d\n", rr, mm, delta*i, nn[i]+0);}' delta=$delta mm=$mm rr=Malaria_T$tt >> $toto.t$tt
    end
  end

  cat $toto.t[0123] >> $toto.1
  cat $toto.1 | sort -V -k 1,1 -k 3,3n -k 2,2 > $toto.2
  cat $toto.1a | sort -V -k 1,1 -k 3,3n -k 2,2 >> $toto.2
  if (-e $toto.3) \rm $toto.3
  foreach mm (truth $methods)
     cat $toto.2 | gawk -F '\t' '{if($2 == mm)print}' mm=$mm >> $toto.3
  end

  echo -n "## $toto :" > $toto
  if (-e $toto.5) \rm $toto.5
  date >> $toto
  echo "## Histogram of length to be aligned in truth dataset,and aligned by each program. Each read is counted only once, at the location of its BAM primary alignment (excluding the secondaries with flag 256)" >> $toto

if (0) then

  foreach rr ($runs)
    echo "# $rr\tTruth" > $toto.4
    cat $toto.3 | gawk -F '\t' '{if($1 != run)next;}{k=$3;n=$4;m=$2;if(k>kMax)kMax=k;kk[k]=1;nn[m,k]=n;}END{kk[0]=1;for (k=0;k<=kMax;k++)if(kk[k]>0)printf("%d\t%d\n", k,nn["truth",k]);}' run=$rr >> $toto.4
    echo "\n\n" >> $toto.4
    cat $toto.4 | scripts/transpose >> $toto.5
    echo "\n\n" >> $toto.5
  end
#

endif

  foreach rr ($runs)
    set ok=`echo $run | gawk '{ok=1;z=substr($1,length($1-1)) ; if (z == "r1" || z == "r2" || z == "r3")ok=0;print ok;}'`
    if ($ok == 0) continue
     
    echo -n "# $rr\tTruth" > $toto.4
    echo "truth" > $toto.4m
    foreach mm ($methods)
      echo -n "\t$mm" >> $toto.4
      echo $mm >> $toto.4m
    end

      cat $toto.4m ZZZZZ $toto.3 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){mm[0+imMax]=$1;imMax++;next;}}{if($1 != run)next;}{k=$3;n=$4;m=$2;if(k>kMax)kMax=k;kk[k]=1;nn[m,k]=n;}END{kk[0]=1;for (k=0;k<=kMax;k++)if(kk[k]>0){printf("%d",k); for(im=0;im<imMax;im++)printf("\t%d", nn[mm[im],k]) ;printf("\n");}}' run=$rr  >> $toto.4

    echo "\n\n" >> $toto.4
    cat $toto.4 | scripts/transpose >> $toto.5
    echo "\n\n" >> $toto.5
  end

  foreach rr (Human_T1 Human_T2 Human_T3 Malaria_T1 Malaria_T2 Malaria_T3)
    echo -n "# $rr\tTruth" > $toto.4
    foreach mm ($methods)
      echo -n "\t$mm" >> $toto.4
    end
    cat $toto.4m ZZZZZ $toto.3  |  gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){mm[0+imMax]=$1;imMax++;next;}}{if($1 != run)next;}{k=$3;n=$4;m=$2;if(k>kMax)kMax=k;kk[k]=1;nn[m,k]=n;}END{kk[0]=1;for (k=0;k<=kMax;k++)if(kk[k]>0){printf("%d",k); for(im=0;im<imMax;im++)printf("\t%d", nn[mm[im],k]) ;printf("\n");}}' run=$rr  >> $toto.4
    echo "\n\n" >> $toto.4
    cat $toto.4 | scripts/transpose >> $toto.5
    echo "\n\n" >> $toto.5
  end

  cat $toto.5 | scripts/transpose | sed -e 's/Truth/Actual reads/g' >> $toto

  \rm $toto.*
  \rm RESULTS/*.tsv
goto phaseLoop

##############################################################################

phase_perfect:

echo "Searching perfects\n\t$runs\n\t$methods"
foreach run ($runs)

  foreach mm ($methods)
    if (-e RESULTS/$mm/$run/$run.$run.2.hits && ! -e RESULTS/$mm/$run/s2g.perfect.list2) then
      echo -n "seraching perfect $mm $run "
      cat RESULTS/$mm/$run/$run.$run.*.hits | gawk -F '\t' '/^#/{next;}{if ($2==$4)print $1;}' | sed -e "s/$run\///" | sort -uV > RESULTS/$mm/$run/s2g.perfect.list
      wc RESULTS/$mm/$run/s2g.perfect.list
    endif
  end

  cat RESULTS/*/$run/*ect.list | sed  -e "s/$run\///" | sort -u > RESULTS/any.$run.perfect.list
 #    cat RESULTS/*/$run/*ect.list | sed -e 's/>//' -e 's/<//' -e "s/$run\///" | sort -u > RESULTS/any.$run.perfect.list
  wc RESULTS/any.$run.perfect.list
end

goto phaseLoop

##############################################################################

phase_introns2:
# locate the introns in baruzzo

if (! -e  HG19.INTRON_DB/known_introns.ace2) then
  echo "// known introns" >  HG19.INTRON_DB/known_introns.ace
    echo "#" >  HG19.INTRON_DB/known_introns.txt
  foreach type (t1 t2 t3)
    echo $type
    zcat Fasta/HG19$type*/HG19*.cig.gz  | gawk -F '\t' '{c=$2;s=$7;if(s=="-")next;n=split($6,aa,",");for(k=1;k<n;k++){split(aa[k],bb,"-");split(aa[k+1],cc,"-");x=bb[2]+1;y=cc[1]-1;printf("%s\t%d\t%d\t%d\t%s\n", c,x,y,y-x+1,type);}}' type=$type  >> HG19.INTRON_DB/known_introns.txt
        zcat Fasta/HG19$type*/HG19*.cig.gz  | gawk -F '\t' '{c=$2;s=$7;if(s=="+")next;n=split($6,aa,",");for(k=1;k<n;k++){split(aa[k],bb,"-");split(aa[k+1],cc,"-");x=bb[2]+1;y=cc[1]-1;printf("%s\t%d\t%d\t%d\t%s\n", c,y,x,y-x+1,type);}}' type=$type  >> HG19.INTRON_DB/known_introns.txt  
  end

  cat HG19.INTRON_DB/known_introns.txt | sort -u  >  HG19.INTRON_DB/known_introns.u.txt

  wc HG19.INTRON_DB/known_introns.txt HG19.INTRON_DB/known_introns.u.txt
  cat  HG19.INTRON_DB/known_introns.u.txt | gawk -F '\t' '/^#/{next}{c=$1;x=$2;y=$3;ln=$4;t=$5;printf("Intron %s__%d_%d\nIntMap %s %d %d\nLength %d\nIn_mRNA %s\nIntron\n\n",c,y,x,c,y,x,ln,t);}}'  > HG19.INTRON_DB/known_introns.ace
endif


setenv iRuns "$runs"
#setenv iRuns "iRefSeq iRefSeq38"
set iMethods="011_SortAlignG5R5 012_SortAlignG3R3 013_SortAlignG3R1 11_MagicBLAST_2018 12_MagicBLAST_2022 13_MagicBLAST_2024 21_HISAT2_4threads 23_HISAT2_16threads 31_STARlong 53_Minimap2_16threads"
set iMethods="$methods"

echo "Creating the INTRON_DB databases"
foreach target (T2T GRCh38 HG19)
  set IDB=$target.INTRON_DB
  if (! -d $IDB) then
    mkdir $IDB
    pushd $IDB
      mkdir database
      ln -s ~/ace/waligner/metaData/wspec.aceview_web_site wspec
      echo y | tace .
    popd
  endif

  if (! -e $target.GFF/mrna2intron.ace) then
     mkdir $target.GFF
     pushd $target.GFF
       ln -s ../Reference_genome/$target.genome.fasta.gz
       ln -s ../Reference_genome/$target.mito.fasta.gz
       if (-e ../Reference_genome/$target.gtf.gz && ! -e ../Reference_genome/$target.mrna2intron.ace) then
         ln -s ../Reference_genome/$target.gtf.gz
         dna2dna -gtf $target.gtf.gz -gtfRemap KT_RefSeq  -o gtf.$target
         touch 
         # construct the intron.ace file needed in INTRON_DB
         cat gtf.$target.introns | gawk -F '\t' '{g=$2;m=$3;c=$4;x1=$5;x2=$6;c=$6;i1=$7;i2=$8;printf("Intron %s__%d_%d\nIntMap %s %d %d\nGene %s\nIn_mRNA %s %d %d\n\n",c,i1,i2,c,i1,i2,g,m,x1,x2);}' > ../Reference_genome/$target.mrna2intron.ace
       endif
    popd
  endif
  if (! -e $IDB/genome.done) then
     tace $IDB <<EOF
       pparse $target.GFF/$target.genome.fasta.gz
       query find sequence DNA
       edit genomic
       pparse $target.GFF/$target.mito.fasta.gz
       pparse Reference_genome/$target.mrna2intron.ace
       save
       quit
EOF
  if ($target == HG19 && -e HG19.INTRON_DB/known_introns.ace) then
    echo "pparse HG19.INTRON_DB/known_introns.ace" | tace $IDB  -no_prompt
  endif



    touch  $IDB/genome.done
  endif
end

#### mrna2intron.ace should look like
# Intron chr1.T2T.NC_060925.1__354592_355554
# IntMap chr1.T2T.NC_060925.1 354592 355554
# Gene SAMD11
# In_mRNA NM_001385640.1 1027 0
#
###

echo "Creating the intron counts"
foreach target (T2T GRCh38 HG19)
  set IDB=$target.INTRON_DB
  if (-e $IDB/introns.ace) continue
    echo '#' > $IDB/introns.tsf
    foreach mm ($iMethods)
      foreach run ($iRuns)
        set target2=`cat Fasta/$run/target`
        echo "$run $target $target2"
        if ($target != $target2) continue
	# sortalign -> introns.trsf
        set ff=RESULTS/$mm/$run/$run/introns.tsf
        if (-e $ff) then
          wc $ff
          cat $ff | gawk '/^#/{next;}{printf("%s\t",mm);print;}' mm=$mm   >> $IDB/introns.tsf
	  continue
	endif
	# other methods -> s2g.Introns
        set ff=RESULTS/$mm/$run/s2g.Introns
        if (-e $ff) then
          wc $ff
          cat $ff | gawk '/^#/{next;}{split($1,aa,":");split(aa[3],bb,"-");printf("%s\t%s__%d_%d\t%s\tiit\t%d\n",mm,aa[1],bb[1],bb[2],run,$2);}' mm=$mm  run=$run >> $IDB/introns.tsf
	  continue
	endif
      end
    end
    cat $IDB/introns.tsf | sort | gawk -F '\t' '/^#/{next;}{if($2 != old){old=$2;split(old,aa,"__");split(aa[2],bb,"_");a1=bb[1];a2=bb[2];if(a1<100 || a2<100){old=0;next;}ln=a2-a1;if(ln<0)ln=-ln;ln+=1;printf("\nIntron %s\nIntron\nIntMap %s %d %d\nLength %d\n", old,aa[1],a1,a2, ln);}printf("de_duo %s__%s %d\n",$1,$3,$5);}END{printf("\n");}' > $IDB/introns.ace
    tace $IDB  <<EOF
      query find intron de_duo
      edit -D de_duo
      parse  $IDB/introns.ace
      bql -o  $IDB/introns.counts.txt select ii,mm,nn,mrna from ii in ?Intron where ii#de_duo, mm in ii->de_duo, nn in mm[1] , mrna in ii#in_mrna
      query find intron In_MRNA
      save
      quit
EOF
end

# in T2T we need the genome named as 1,2,3  to be able to read the intron feet 
if (! -e T2T.INTRON_DB/T2T.renamed.genome.fasta.gz) then
  zcat Reference_genome/T2T.genome.fasta.gz | gawk '/^>/{split($1,aa,".");gsub("chr","",aa[1]);print aa[1];next;}{print;}' | gzip > T2T.INTRON_DB/T2T.renamed.genome.fasta.gz
  echo "pparse  T2T.INTRON_DB/T2T.renamed.genome.fasta.gz" | tace T2T.INTRON_DB -no_prompt 
endif

echo "Creating the intron feet, only new introns will be analyzed"
foreach target (T2T GRCh38)
  set IDB=$target.INTRON_DB
  echo "altintrons --setFeet --db $IDB -p toto"
        altintrons --setFeet --db $IDB -p toto
endif
foreach target (T2T GRCh38)
  set IDB=$target.INTRON_DB
    tace $IDB  <<EOF	
      bql -o  $IDB/introns.counts.txt select ii,mm,nn,mrna from ii in ?Intron where ii#de_duo && ii#gt_ag , mm in ii->de_duo, nn in mm[1] , mrna in ii#in_mrna
      query find intron In_MRNA
      quit
EOF
end


 cat T2T.INTRON_DB/introns.counts.txt | gawk -F '\t' '{any[$2]++;anyC[$2]+=$3;}/NULL/{new[$2]++;newC[$2]+=$3;}/In_mRNA/{known[$2]++;knownC[$2]+=$3;}END{for (m in any)printf("T2T\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",m,any[m],known[m],new[m],anyC[m],knownC[m],newC[m]);}' | sort > _a
 cat GRCh38.INTRON_DB/introns.counts.txt | gawk -F '\t' '{any[$2]++;anyC[$2]+=$3;}/NULL/{new[$2]++;newC[$2]+=$3;}/In_mRNA/{known[$2]++;knownC[$2]+=$3;}END{for (m in any)printf("GRCh38\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",m,any[m],known[m],new[m],anyC[m],knownC[m],newC[m]);}' | sort > _b
 cat HG19.INTRON_DB/introns.counts.txt | gawk -F '\t' '{any[$2]++;anyC[$2]+=$3;}/NULL/{new[$2]++;newC[$2]+=$3;}/In_mRNA/{known[$2]++;knownC[$2]+=$3;}END{for (m in any)printf("HG19\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",m,any[m],known[m],new[m],anyC[m],knownC[m],newC[m]);}' | sort > _h

 phase_introns:
setenv iRuns "$runs"
#setenv iRuns "iRefSeq iRefSeq38"
set iMethods="011_SortAlignG5R5 012_SortAlignG3R3 013_SortAlignG3R1 11_MagicBLAST_2018 12_MagicBLAST_2022 13_MagicBLAST_2024 21_HISAT2_4threads 23_HISAT2_16threads 31_STARlong 53_Minimap2_16threads"
set iMethods="$methods"

 
echo -n "# Introns counts $SV : " > _c
date >> _c
echo -n "#Run\tMethod" >> _c

echo -n "# Introns supports $SV : " > _cC
date >> _cC
echo -n "#Run\tMethod" >> _cC

foreach mm ($iMethods)
    echo $mm | gawk '{printf("\t%s",$1);}' >> _c
    echo $mm | gawk '{printf("\t%s",$1);}' >> _cC
end
foreach run ($iRuns)
  echo -n "\n$run Known Intron\t" >> _c
  echo -n "\n$run Known Intron\t" >> _cC
  foreach mm ($iMethods)
    cat _a _b _h | gawk '{split($2,aa,"__");if(aa[1]==m && aa[2]==r)n=$4;}END{printf("\t%d",n+0);}' m=$mm r=$run >> _c
    cat _a _b _h | gawk '{split($2,aa,"__");if(aa[1]==m && aa[2]==r)n=$7;}END{printf("\t%d",n+0);}' m=$mm r=$run >> _cC    
  end
end
echo >> _c
echo >> _cC
foreach run ($iRuns)
  echo -n "\n$run New Intron\t" >> _c
  echo -n "\n$run New Intron\t" >> _cC
  foreach mm ($iMethods)
    cat _a _b _h | gawk '{split($2,aa,"__");if(aa[1]==m && aa[2]==r)n=$5;}END{printf("\t%d",n+0);}' m=$mm r=$run >> _c
    cat _a _b _h | gawk '{split($2,aa,"__");if(aa[1]==m && aa[2]==r)n=$8;}END{printf("\t%d",n+0);}' m=$mm r=$run >> _cC
  end
end
echo >> _c
echo >> _cC
\cp _c COMPARE/introns_counts_known_new.$SV.txt
\cp _cC COMPARE/introns_supports_known_new.$SV.txt
wc COMPARE/introns_*_known_new.$SV.txt

# ROC curve, truth = gt_ag  , not useful when we align the iRefSeq38
 cat GRCh38.INTRON_DB/introns.tsf | gawk '{if($1==m){print $7;}}' m=011_SortAlignG5R5 | tags | sort -k 2nr | head -8
 cat GRCh38.INTRON_DB/introns.tsf | gawk '{if($1==m){print $5, $7;}}' m=011_SortAlignG5R5 | sort -k 1nr | gawk '{if($2=="gt_ag")y++;else x++;if(y%1000 == 0)printf("%d\t%d\n",x,y);}END{printf("%d\t%d\n", x,y)}'

 
goto phaseLoop
# iRefSeq38 NR_024540.1 is aligned but has no produced intron
set mm=013_SortAlignG3R1
set run=iRefSeq
set ff=RESULTS/$mm/$run/$run/introns.tsf
set IDB=T2T.INTRON_DB
cat $ff | gawk '/^#/{next;}{printf("%s\t%s.NoI\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,$5,$6) ;}'  > $IDB/introns.tsf
    cat $IDB/introns.tsf | sort | gawk -F '\t' '/^#/{next;}{if($2 != old){old=$2;split(old,aa,"__");split(aa[2],bb,"_");a1=bb[1];a2=bb[2];if(a1<100 || a2<100){old=0;next;}ln=a2-a1;if(ln<0)ln=-ln;ln+=1;printf("\nIntron %s\nIntron\nIntMap %s %d %d\nLength %d\n", old,aa[1],a1,a2, ln);}printf("de_duo %s__%s %d\n",$1,$3,$5);}END{printf("\n");}' > $IDB/introns.ace

set run=iRefSeq38
set ff=RESULTS/$mm/$run/$run/s2g.Introns
set IDB=GRCh38.INTRON_DB
cat $ff | gawk '/^#/{next;}{printf("%s\t%s.NoI\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,$5,$6) ;}'  > $IDB/introns.tsf
    cat $IDB/introns.tsf | sort | gawk -F '\t' '/^#/{next;}{if($2 != old){old=$2;split(old,aa,"__");split(aa[2],bb,"_");a1=bb[1];a2=bb[2];if(a1<100 || a2<100){old=0;next;}ln=a2-a1;if(ln<0)ln=-ln;ln+=1;printf("\nIntron %s\nIntron\nIntMap %s %d %d\nLength %d\n", old,aa[1],a1,a2, ln);}printf("de_duo %s__%s %d\n",$1,$3,$5);}END{printf("\n");}' > $IDB/introns.ace

##############################################################################
##############################################################################
## Create a worm database to verify the wiggles

if (! -d WormDB) mkdir WormDB
if (! -e WormDB/database) then
  pushd WormDB
    ln -s ~/ace/waligner/metaData/wspec.aceview_web_site wspec
    mkdir database
    echo y | tace .
    if (! -d DATA) mkdir DATA
      pushd DATA
      rr
      tace .. < _r
    popd
  popd
endif

goto phaseLoop

cat T2T.chr7.fasta | gawk '/^>/{next}{for(i=1;i<=1000;i++)printf(">S%d\n%stgctagatcgaatc\n",i,substr($1,10000000+100*i+1,200));}' > test.chr7.fasta
cat T2T.chr7.fasta | gawk '/^>/{next}{for(i=1;i<=1000;i++)printf(">S%d\ntgctagatcgaatc%s\n",i,substr($1,10000000+100*i+1,200));}' > test.chr7.fasta

##############################################################################
##############################################################################
## Done

phase6:
phaseLoop:
  echo -n "$1 $SV done : "
  date

##############################################################################
##############################################################################
##############################################################################
##############################################################################
