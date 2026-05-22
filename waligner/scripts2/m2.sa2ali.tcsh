#!bin/tcsh
set run=$1
set dd=tmp/SA/$run 
set toto=$dd/sa2ali.ace

echo "Ali $run\nRun $run" > $toto
echo "-D Letter_profile" >> $toto
echo "-D stranding Introns" >> $toto
echo "-D Candidate_introns" >> $toto
echo "Counts\nStrandedness\nAli\nUnicity\nGene_expression\nSponge\nAlignments\nLetter_profile\nATGC_kb\nPair_fate\nErrors\nComputer_ressource" >> $toto
#RawCounts
set ff=$dd/runStats.tsf
cat $ff | gawk -F '\t' "/^$run/{print;}" > $ff.1
cat $ff.1 | gawk -F '\t'  '{if($2=="Reads"){nr=$4;next;}}{if($2=="Bases"){nb=$4;printf("Raw_data %s Id %s Accepted %g kb %.2f bp_per_tag\n", nr,nr, nb/1000,nb/nr);}}'  >> $toto
cat $ff.1 | gawk -F '\t'  '{if($2=="Reads"){nr=$4;next;}}{if($2=="Bases"){nb=$4;printf("Accepted %s Seq %s Tags %g kb %.2f bp_per_tag\n", nr,nr, nb/1000,nb/nr);}}'  >> $toto

set ff=$dd/$run.overhang.5prime.tsf
#Entry_adaptor_clipping
set ff=$dd/$run.overhang.3prime.tsf
#Exit_adaptor_clipping

set ff=$dd/runStats.tsf
cat $ff.1 | gawk -F '\t' '/Length_distribution_1_5_50_95_99_mode_av/{printf("Length_distribution_1_5_50_95_99_mode_av %s %s %s %s %s %s %s\n", $4,$5,$6,$7,$8,$9,$10);}' >> $toto
cat $ff.1 | gawk -F '\t'  '/Reads_aligned_in_class/{s=$10;p=$6;m=$8;cl=substr($2,24);if(p+m >100) printf("Stranding %s %s %s plus %s minus\n", cl,s,p,m);}' >> $toto
cat $ff.1 | gawk -F '\t'  '/Reads_aligned_in_class/{cl=substr($2,24);nr[cl]=$6+$8} /Bases_aligned_in_class/{cl=substr($2,24);nb[cl]=$6+$8;printf("nh_Ali %s %s seq %s tag %.3g kb %.2f bp\n",cl,nr[cl],nr[cl],nb[cl]/1000,nb[cl]/nr[cl]);}' >> $toto
cat $ff.1 | gawk -F '\t'  '/Reads_aligned_once/{nn[1]=$4;s[1]=$5;}/Reads_multi_aligned__/{cl=substr($2,22)+0;nn[cl]=$4;s[cl]=$5;}END{printf("Unicity any ");for(cl=1;cl<=10;cl++)printf(" %.3g", nn[cl]);printf("\n");}' >> $toto

set gg=$dd/geneCounts.tsf

set fw=$dd/wiggleCumuls.tsf
if (-e $fw.gz) gunzip -f $fw.gz
if (-e $fw) then
  cat $fw | gawk -F '\t' '/^Any/{CDS=$5;UTR=$6;intronic=$7;intergenic=$8;if (CDS+0>0)printf("S_1_CDS %.3f Mb aligned\n", CDS/1000000);if (UTR+0>0)printf("S_1_UTR %.3f Mb aligned\n", UTR/1000000);if (intronic+0>0)printf("S_1_intronic %.3f Mb aligned\n", intronic/1000000);if (intergenic+0>0)printf("S_1_intergenic %.3f Mb aligned\n", intergenic/1000000);}' >> $toto
endif

cat $ff | gawk -F '\t' '/Unaligned_reads/{nr=$4;}/Unaligned_bases/{nb=$4;printf("Unaligned %s Seq %d Tags %.3g kb\n", nr, nr, nb/1000) ;}' >> $toto

cat $ff | gawk -F '\t' '/Perfect_reads/{nr=$4;printf("Perfect_reads %.3g %d\n", nr, nr) ;}' >> $toto
cat $ff | gawk -F '\t' '/Complex_reads/{nr=$4;printf("Complex_reads %.3g %d\n", nr, nr) ;}' >> $toto
cat $ff | gawk -F '\t' '/Partial_reads/{nr=$4;printf("Partial_reads %.3g %d\n", nr, nr) ;}' >> $toto
cat $ff | gawk -F '\t' '/Aligned_pairs/{nr=$4;printf("Aligned_fragments %g\n", nr) ;}' >> $toto
cat $ff | gawk -F '\t' '/Compatible_pairs/{nr=$4;printf("Compatible_pairs %g\n", nr) ;}' >> $toto
cat $ff | gawk -F '\t' '/Circle_pairs/{nr=$4;printf("Circle_pairs %g\n", nr) ;}' >> $toto
cat $ff | gawk -F '\t' '/Non_compatible_pairs/{nr=$4;printf("Non_compatible_pairs %g\n", nr) ;}' >> $toto

set fe=$dd/salign.err
if (-e $fe) then
  cat $fe | gawk '/TIMING/{if ($4=="U") {t=$5;}}END{if(t+0>0)printf("CPU sortAlign %d seconds\n", t);}' >> $toto
  cat $fe | gawk '/TIMING/{if ($2=="E") {t=0;x=$3;n=split(x,aa,":");if(n==2)t=60*aa[1]+aa[2];if(n==3)t=3600*aa[1]+60*aa[2]+aa[3];}}END{if(t+0>0)printf("Elapsed sortAlign %d seconds\n", t);}' >> $toto
  cat $fe | gawk '/TIMING/{if ($8=="M") {m=$9;}}END{if(m+0>0)printf("Max_memory sortAlign %d Gb\n", m/1000000);}' >> $toto
  cat $fe | gawk '/TIMING/{if ($10=="P") {p=$11;}}END{if(p+0>0)printf("Multi_threading sortAlign %.2f average_running_threads\n", p/100.0);}' >> $toto
endif


set ff=$dd/$run.letterProfile.tsf
if (-e $ff) then
  echo "\nAli $run\nRun $run" >> $toto
  cat $ff | gawk -F '\t' '/^#/{next;}{fr=substr($1,length($1),1);if(fr=="r")fr="f2";else fr="f1"; if($2+0>0)printf("Letter_profile %s %d %d %d %d %d %d %d %d %d %d %d\n",fr,$2,$10,$11,$12,$13,$14,$4,$5,$6,$7,$8)}' >> $toto
endif

set ff=$dd/runErrors.tsf
if (-e $ff) then
  cat $ff | gawk -F '\t' '/^#/{next;}{if (run != $1) {printf("\nAli %s\n", $1);run=$1;}t=$2;if(t=="Any"){u=t;if(n>0)printf("Cumulated_mismatches %d\n",$4);}if(substr(t,2,1)==">")u=tolower(t);if(substr(t,1,3)=="Ins"){k=length(substr(t,4));u=substr("++++++++",1,k) tolower(substr(t,4));}if(substr(t,1,3)=="Del"){k=length(substr(t,4));u=substr("-------------",1,k) tolower(substr(t,4));}if ($4)printf("Error_profile f1 %s %d\n", u, $4);}END{print "\n";}' >> $toto
endif

set ff=$dd/polyA.tsf
if (-e $ff) then
  cat $ff | gawk -F '\t' '/^#/{next;}/^PolyA/{next;}{run=$2;nn[run]++;nnn[run]+=$5;}END{for(run in nn){if(nn[run]>0)printf("Ali %s\nSLs pA %d sites %d supports\n\n", run,nn[run],nnn[run]);}}' >> $toto
endif

set ff=$dd/SL.tsf
if (-e $ff) then
  cat $ff | gawk -F '\t' '/^#/{next;}/^SL/{next;}{k=split($1,aa,"___");if(k!=2)next;sl=0+substr(aa[2],3);if(sl<1||sl>20)next;run=$2;nn[run]++;ns[run,sl]++;nns[run,sl]+=$5;}END{for(run in nn){printf("\nAli %s\n",run);for(sl=1;sl<20;sl++)if(nns[run,sl]>0)printf("SLs SL%d %d sites %d supports\n", sl, ns[run,sl],nns[run,sl]);}}END{print "\n";}' >> $toto
endif


echo >> $toto

if (0) then
  set ff=$dd/runStats.tsf 
  cat $ff | gawk -F '\t'  'BEGIN{s2=100;n=0;}/Reads_aligned_in_class/{r=run;s=$10;p=$6;m=$8;cl=substr($2,24);if(p+m >100) printf("Run %s\nObserved_strandedness_in_ns_mapping %s %s %d plus %d minus\n\n", run,cl,s,p,m);}' run=$run > $toto.1
  cat $toto.1 >> $toto
endif

wc $toto
