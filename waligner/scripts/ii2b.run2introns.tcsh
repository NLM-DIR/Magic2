#!bin/tcsh -f

set run=$1
set uu=$2

# centralize the introm support from lane to run
if (! -e tmp/INTRONRUNS/$run/$run.u.intronSupport.counts.gz) then
  zcat tmp/INTRONLANES/$run/f2.*.intronSupport.u.gz | gawk -F '\t' '/^Intron/{ii=$2"\t"$4;nn[ii]+=$6;}END{for(ii in nn)printf("%s\t%d\n",ii,nn[ii]);}' | gzip > tmp/INTRONRUNS/$run/$run.u.intronSupport.counts.gz
  zcat tmp/INTRONLANES/$run/f2.*.intronSupport.nu.gz | gawk -F '\t' '/^Intron/{ii=$2"\t"$4;nn[ii]+=$6;}END{for(ii in nn)printf("%s\t%d\n",ii,nn[ii]);}' | gzip > tmp/INTRONRUNS/$run/$run.nu.intronSupport.counts.gz
endif

if (! -e tmp/INTRONRUNS/$run/$run.u.intronSupport.ace.gz) then
  bin/bestali  -intronSupport2ace  -gzo -o tmp/INTRONRUNS/$run/$run.$uu -inFileList tmp/INTRONRUNS/$run/$run.$uu.list -run $run
endif


exit 0
