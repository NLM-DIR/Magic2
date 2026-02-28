#!bin/tcsh -f

set zone=$1
set nRuns=$2
set Strategy=$3

 echo -n "s21: select the characteristic SNPs in $zone start : "
 date
   #                                                          $ $24=measured_in $26=Reference $31= allelefreq in cohort
  cat tmp/SNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt | gawk -F '\t' '{if ($24  > .050 * n && $26 > 1 && $28 + $29 + $30 > 1 && ($31 >= 5 && $31 < 95)) print;}' n=$nRuns > tmp/SNP_DB/$zone/$MAGIC.characteristic_snp.txt

if ($Strategy == RNA_seq) then
  cat tmp/SNP_DB/$zone/$MAGIC.characteristic_snp.txt | gawk -F '\t' '/^#/{next;}/A:G/{next;}/a:g/{next;}{p=$5;if(index(p,":Sub:")==0)next;print $5;}' > tmp/SNP_DB/$zone/$MAGIC.characteristic.list
  
  cat tmp/SNP_DB/$zone/$MAGIC.characteristic.list ZZZZZ tmp/SNPH/$zone/$MAGIC.snp.sorted | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}{z=$1 ":" $2 ":" $3;if(ok[z]==1)print;}'  | gzip >  tmp/SNP_DB/$zone/$MAGIC.characteristic.snp.gz
else
  cat tmp/SNP_DB/$zone/$MAGIC.characteristic_snp.txt | gawk -F '\t' '/^#/{next;}{print  $5 ;}' > tmp/SNP_DB/$zone/$MAGIC.characteristic.list
  
  cat tmp/SNP_DB/$zone/$MAGIC.characteristic.list ZZZZZ tmp/SNPH/$zone/$MAGIC.snp.sorted | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}{z=$1 ":" $2 ":" $3;if(ok[z]==1)print;}' | gzip >  tmp/SNP_DB/$zone/$MAGIC.characteristic.snp.gz

endif


echo -n "s21 : count SNPs per type and per run"

if (0) then
  cat tmp/SNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt 

  cat $toto | gawk -F '\t' '{if (mx < 1){for(i=1;i<=NF;i++)if(tolower(substr($i,1,24)) == "maximal allele frequency")mx=i;for(i=mx+1;i<=NF;i++)run[i]=$i;nRuns=NF;}}{type=$11;k=0;for(i=mx+1;i<=NF;i++){if($i > -10) {x = 2 * $i/100; k+= x;zz[type,i]+= x;}} if(k>0) union[type]++;total[type]+=k;}END{for(t in union){printf("union\t%s\t%.2f\n",t,union[t]);printf("total\t%s\t%.2f\n",t,total[t]);for(i=mx+1;i<=nRuns;i++)if (zz[t,i]>0)printf("%s\t%s\t%.2f\n",run[i],t,zz[t,i]);}}'
endif


 echo -n "done : "
 date

exit 0
