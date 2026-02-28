#!bin/tcsh -f

set group=$1
set zone=$2 

      set flist=tmp/SNP_BRS/$group/s2.$zone.file_list

echo $flist
cat $flist

gunzip -c `cat $flist` | bin/snp -BRS_merge -gzo -run $group -o tmp/SNP_BRS/$group/$zone.BRS.u 

touch tmp/SNP_BRS/$group/s2.$zone.done
