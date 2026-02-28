#!bin/tcsh -f

set dd=$1
set ss=$2
 set srr=`echo $ss | gawk -F '___' '{print $1}'`
 set ff=`echo $ss | gawk -F '___' '{print $2}'`
 set p=`echo $ss | gawk -F '___' '{print $3}'`

echo "$srr $ff $p"
pushd $dd

  if ($p == Paired_end) then
    echo "fastq-dump --split-files --fasta 0 --gzip $srr"
    fastq-dump --split-files --fasta 0 --gzip $srr
  else
    echo "fastq-dump --fasta 0 --gzip $srr"
    fastq-dump --fasta 0 --gzip $srr
  endif
  if (-e $srr'_3.fasta.gz') then
  # HACK 2024_11_01 in illumina UMI the short UMI file may be in position 1 or in position 3, mv it to position 3
    set n1=`ls -ls $srr'_1.fasta.gz' | gawk '{print $6;}'`
    set n2=`ls -ls $srr'_2.fasta.gz' | gawk '{print $6;}'`
    set n3=`ls -ls $srr'_3.fasta.gz' | gawk '{print $6;}'`
    echo "ls -ls $srr"'*.fasta.gz'  > $srr.m
    echo $n1 $n2 $n3 | gawk '{if (2 * $1 < $2 && 2 * $1 < $3) printf("mv %s_1.fasta.gz %s_0.fasta.gz \n mv %s_2.fasta.gz %s_1.fasta.gz \n mv %s_3.fasta.gz %s_2.fasta.gz \n mv %s_0.fasta.gz %s_3.fasta.gz \n", r,r,r,r,r,r,r,r);}' r=$srr >> $srr.m
    echo "ls -ls $srr"'*.fasta.gz'  >> $srr.m
    source $srr.m
  endif
popd


 
