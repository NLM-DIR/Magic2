#!bin/tcsh -f

# do not use -fe, some files may not exist in a perfect quality situation
set run=$1


######### C2 sec.tsv start
# export the alignment statistics in seqc.tsv format, 
# since seqc.tsv statitics acn also be derived from BAM files
# this allows to comprae MAGIC to other aligners
    
echo > tmp/$MAGIC_COUNT_DIR/$run/tsv2

set ok=1 

# collate the lane seqc.tsv stats
# if a lane is missing, (ok == 0), do not create the run report

foreach lane (`cat Fastc/$run/LaneList`)
  if (-e  tmp/$MAGIC_COUNT_DIR/$lane.seqc.tsv) then
    cat tmp/$MAGIC_COUNT_DIR/$lane.seqc.tsv >> tmp/$MAGIC_COUNT_DIR/$run/tsv2
  else
    set ok=0
  endif
end

# if all lanes are present prepare the run report
if ($ok == 1 && -e scripts/aliqc.pyZZZ) then

  # grab the number of reads to be aligned
  set nr=`cat MetaDB/$MAGIC/ali.ace | gawk '/^Ali/{ok=0;gsub(/\"/,"",$2);if($2==r)ok=1;next;}{if(ok<1)next;}/^Accepted/{nr=$4;next;}END{print 0+nr;}' r=$run`
  # synthetize the seqc.tsv stats of the whole run 
  cat  tmp/$MAGIC_COUNT_DIR/$run/tsv2 | scripts/aliqc.py -m -r $run..magic -o tmp/$MAGIC_COUNT_DIR/$run --nreads $nr

endif

# Clean up
\rm  tmp/$MAGIC_COUNT_DIR/$run/tsv2

if ($ok == 0) return 0

######### C2 sec.tsv end


# regroup all the lanes per run
# computer ressources
echo "report the CPU and RAM"


  set toto=tmp/$MAGIC_COUNT_DIR/$run/c2.toto.alistats
  
  echo "Run $run" >  $toto.ace
  echo "Ali $run" >> $toto.ace 

  echo >> $toto.ace
  echo "Ali $run" >> $toto.ace
  echo "-D Ali" >> $toto.ace
  echo "-D Computer_ressource" >> $toto.ace
  echo "-D Strandedness" >> $toto.ace
  echo "-D Unicity" >> $toto.ace
  echo "-D Rejected_alignments" >> $toto.ace
  echo "-D Rejected_too_many_hits" >> $toto.ace
  echo "-D Alignments" >> $toto.ace
  echo "-D Errors" >>  $toto.ace
  echo "-D Pair_fate" >>  $toto.ace

  if (-e  tmp/INTRONRUNS/$run/$run.u.intronSupport.ace.gz) then
    gunzip -c  tmp/INTRONRUNS/$run/$run.u.intronSupport.ace.gz | head -3 | gawk '/^Anti_run/{a=$6;next;}/^Run_U/{b=$6;next;}END{if(a+b>10000)printf("Stranding Exons_juntions %.2f  %d plus %d minus 0 ambiguous\n",run, 100*a/(a+b), b, a);}' run=$run >> $toto.ace
  endif

  cat Fastc/$run/LaneList | gawk '{n++}END{printf ("Number_of_lanes %d\n", n);}' >> $toto.ace

  cat tmp/$MAGIC_COUNT_DIR/$run/*.jobstats.preace | gawk '{if ($2 != r)next;z=$2 " " $4;run[z]=$2;tt[z]=$4;}/^CPU/{cpu[z]+=$5;}/^Memory/{s=$5;if(s>mem[z])mem[z]=s;}END{for(z in mem)printf("Max_Memory %s %d Mb\n",tt[z],mem[z]);for(z in cpu)printf("CPU %s %d seconds\n",tt[z],cpu[z]);}' r=$run | sort >> $toto.ace

#######
## h_Ali, nh_Ali,  stranding, multiplicity

# hits per prefix
set prefix=f2
if (-e $toto.hits1) \rm $toto.hits1

  # create a non empty $toto.txt
  echo 'toto' > $toto.txt
  foreach lane (`cat Fastc/$run/LaneList`)
    cat tmp/$MAGIC_COUNT_DIR/$lane.count | gawk  -F '\t' '/^#/{next}/^$/{next}{printf("%s\t",run);print}' run=$run >> $toto.txt
  end
  cat  $toto.txt | gawk -F '\t' '/HITS/{z= $1 "\t" $3 "." prefix  ; if (t[z]<1){t[z]=1;nt++;t2z[nt]=z;z2t[z]=nt;}it=z2t[z];if(imax<NF)imax=NF;for (i=5;i<=NF;i++)nn[it,i]+=$i;}END{for(it=1;it<=nt;it++){printf("\n%s",t2z[it]);for (i=5;i<=imax;i++)printf("\t%d",nn[it,i]);}printf("\n");}' prefix=$prefix  >> $toto.hits1

sort  $toto.hits1 >  $toto.hits
cat $toto.hits |  gawk -F '\t' '{z = $1  ; if($3<1)next;  n=$4+$5; namb = $6 ; if (n && $2 != "ZZZ0_SpikeIn") printf ("Stranding %s %.3f  %d plus %d minus %d ambiguous\n", $2, 100 * $4/n, $4, $5, namb); }'  >>  $toto.ace

\rm    $toto.hits1

# any hits
cat tmp/$MAGIC_COUNT_DIR/$run/*.count | gawk  -F '\t' '/^#/{next}/^$/{next}{printf("%s\t",run);print}' run=$run > $toto.txt
  cat  $toto.txt | gawk -F '\t' '/HITS/{z= $1 "\t" $3  ; if (t[z]<1){t[z]=1;nt++;t2z[nt]=z;z2t[z]=nt;}it=z2t[z];if(imax<NF)imax=NF;for (i=5;i<=NF;i++)nn[it,i]+=$i;}END{for(it=1;it<=nt;it++){printf("\n%s",t2z[it]);for (i=5;i<=imax;i++)printf("\t%d",nn[it,i]);}printf("\n");}'  > $toto.hits

cat tmp/$MAGIC_COUNT_DIR/$run/*.count > $toto.txt
  cat  $toto.txt | gawk -F '\t' '/MULT/{z= $1  "\t" $2 ; if (t[z]<1){t[z]=1;nt++;t2z[nt]=z;z2t[z]=nt;}it=z2t[z];if(imax<NF)imax=NF;for (i=3;i<=NF;i++)nn[it,i]+=$i;}END{for(it=1;it<=nt;it++){printf("\n%s",t2z[it]);for (i=3;i<=imax;i++)printf("\t%d",nn[it,i]);}printf("\n");}' | sort > $toto.mult

  cat $toto.hits |  gawk -F '\t' '{z = $1  ; if($3<1)next;printf ("nh_Ali  %s %.2f seq %.2f tags  %.3f kb_ali  %.2f bp_av_ali %.3f kb_clip %.2f bp_av_clip\n",$2, $3,$7,$8/1000,$8/$7,$9/1000,$9/$7);d=1;if ($14>0)d=$14; printf ("h_Ali %s %d seq %d tags %d kb_ali %.2f bp_av_ali %d kb_clip %.2f bp_av_clip\n", $2, $10,$14,$15/1000,$15/d,$16/1000,$16/d) ;  n=$4+$5; namb = $6 ; }'  >>  $toto.ace

  cat $toto.mult |  gawk -F '\t' '{z = $1 ; if($3<10)next;printf ("Unicity Number_of_targets 1 2 3 4 5 6 7 8 9 10 -2\n") ;printf ("Unicity %s",$2) ; for(i=3;i<=NF;i++)printf(" %s",$i);printf("\n");}' >>  $toto.ace

\rm $toto.txt $toto.hits $toto.mult

###### report the errors position, strand dependent, only works for COUNT

echo toto >  $toto.pos
echo toto >  $toto.amb
echo toto > $toto.FB
echo toto > $toto.LB
echo toto > $toto.err_spike

      gunzip -c  tmp/$MAGIC_COUNT_DIR/$run/*.errorProfile.gz   | grep RPOS >> $toto.err_spike
      gunzip -c tmp/$MAGIC_COUNT_DIR/$run/*.errorProfile.gz |  gawk -F '\t' '/^POS/{prefix=$2;x=$3;prefixes[prefix]=1;xx[x]=1;p[$2,x]+=$4;}END{for(prefix in prefixes)for(x in xx)if(0+p[prefix,x]>0 && x+0<=1000)printf("Error_position  %s %d %d \n", prefix, x, p[prefix,x]) ;}'  | sort -k 1,1 -k 2,2 -k 3,3n >> $toto.pos
      gunzip -c tmp/$MAGIC_COUNT_DIR/$run/*.errorProfile.gz |  gawk -F '\t' '/^FB/{prefix=$2;x=$3;prefixes[prefix]=1;xx[x]=1;p[$2,x]+=$4;}END{for(prefix in prefixes)for(x in xx)if(0+p[prefix,x]>0 && x+0<=1000)printf("FB %s %d %d \n", prefix, x, p[prefix,x]) ;}'  | sort -k 1,1 -k 2,2 -k 3,3n >> $toto.FB
      gunzip -c tmp/$MAGIC_COUNT_DIR/$run/*.errorProfile.gz |  gawk -F '\t' '/^LB/{prefix=$2;x=$3;prefixes[prefix]=1;xx[x]=1;p[$2,x]+=$4;}END{for(prefix in prefixes)for(x in xx)if(0+p[prefix,x]>0 && x+0<=1000)printf("LB %s %d %d \n", prefix, x, p[prefix,x]) ;}'  | sort -k 1,1 -k 2,2 -k 3,3n >> $toto.LB
      
gunzip -c tmp/$MAGIC_COUNT_DIR/$run/*.errorProfile.gz |  gawk -F '\t' '/^NNPOS/{prefix=$2;x=$3;prefixes[prefix]=1;xx[x]=1;p[$2,x]+=$4;}END{for(prefix in prefixes)for(x in xx)if(0+p[prefix,x]>0 && x+0<=1000)printf("Ambiguous_position %s  %s %d \n", prefix, x, p[prefix,x]) ;}' >>  $toto.amb

gunzip -c tmp/$MAGIC_COUNT_DIR/$run/*.errorProfile.gz | gawk -F '\t' '/^ERRALI/{if(NF==4){s=$2;nkb[s]+=$3;}if(NF==3){s="f1";nkb[s]+=$2;}next;}/^ERRS/{s=$2;ss[s]=1;t=$3;if(t>tmax)tmax=t;tt[t]=$4;ns[s,t]+=$5;nns[s]+=$5;next;}/^ERR/{s="f1";ss[s]=1;t=$2;if(t>tmax)tmax=t;tt[t]=$3;ns[s,t]+=$4;nns[s]+=$4;next;}END{tmax++;tt[tmax]="Any";for(s in ss){ns[s,tmax]=nns[s];for(t=1;t<=tmax;t++)if(nkb[s]>0)if(t==tmax || ns[s,t]>0)printf("Error_profile %s %s %d %d %.1f\n", s, tt[t],ns[s,t],nkb[s]/1000,1000* ns[s,t]/(nkb[s]));}}' >> $toto.ace


cat $toto.err_spike | gawk '/RPOS/{ z = $2 " " $3 " \"" $6 "\"" ; nt[z] += $5 ; ne[z] += $7 ;}END {for (z in nt)if(ne[z]>0 && 3*ne[z] > nt[z])printf ("Error_spike %s %d %d\n", z, ne[z], nt[z]);  }'  >>   $toto.ace

cat $toto.pos | sort  -k 1,1 -k 2,2 -k 3,3n | gawk '/^Error_position/{nf[$2]++;nn[$3,$2]=$4;if($3>mx)mx=$3;}END{for (x=1;x<=mx;x++){printf("Error_position %d", x); if(nf["f"]>0)printf(" %d",nn[x,"f"]); if(nf["f1"]>0)printf(" %d",nn[x,"f1"]); if(nf["f2"]>0)printf(" %d",nn[x,"f2"]); if(nf["f3"]>0)printf(" %d",nn[x,"f3"]);printf("\n");}}' >>   $toto.ace

cat $toto.FB | sort  -k 1,1 -k 2,2 -k 3,3n | gawk '/^FB/{nf[$2]++;nx[$3]++;nn[$3,$2]=$4;if($3>mx)mx=$3;}END{for (x=1;x<=mx;x++){if (nx[x]>0){printf("First_base_aligned %d", x); if(nf["f"]>0)printf(" %d",nn[x,"f"]); if(nf["f1"]>0)printf(" %d",nn[x,"f1"]); if(nf["f2"]>0)printf(" %d",nn[x,"f2"]); if(nf["f3"]>0)printf(" %d",nn[x,"f3"]);printf("\n");}}}' >>   $toto.ace

cat $toto.LB | sort  -k 1,1 -k 2,2 -k 3,3n | gawk '/^LB/{nf[$2]++;nx[$3]++;nn[$3,$2]=$4;if($3>mx)mx=$3;}END{for (x=1;x<=mx;x++){if (nx[x]>0){printf("Last_base_aligned %d", x); if(nf["f"]>0)printf(" %d",nn[x,"f"]); if(nf["f1"]>0)printf(" %d",nn[x,"f1"]); if(nf["f2"]>0)printf(" %d",nn[x,"f2"]); if(nf["f3"]>0)printf(" %d",nn[x,"f3"]);printf("\n");}}}' >>   $toto.ace

cat $toto.amb | gawk '/^Ambiguous_position/{nf[$2]++;nx[$3]++;nn[$2,$3]=$4;if($3>mx)mx=$3;}END{for (x=1;x<=mx;x++){if (nx[x]>0){printf("Ambiguous_position %d", x); if(nf["f"]>0)printf(" %d",nn["f",x]); if(nf["f1"]>0)printf(" %d",nn["f1",x]); if(nf["f2"]>0)printf(" %d",nn["f2",x]); if(nf["f3"]>0)printf(" %d",nn["f3",x]);printf("\n");}}}' >>   $toto.ace
 
     gunzip -c  tmp/$MAGIC_COUNT_DIR/$run/*.aliProfile.gz   | gawk '/^Perfect_reads/{n+=$4;}END{if(n > 0)printf ("Perfect_reads %d %d\n", n,n);}'  >> $toto.ace

###### Read fate

#this one is too slow and should be split, nic had an out of RAM problem 2013_03_25 

foreach lane (`cat Fastc/$run/LaneList`)
  if (! -e tmp/$MAGIC_COUNT_DIR/$lane.badcounts.txt) then
    gunzip -c tmp/$MAGIC_COUNT_DIR/$lane.hits.gz ZZZZZ.gz tmp/$MAGIC_COUNT_DIR/$lane.too_many_hits.gz |  gawk -F '\t' '/^ZZZZZ/{zz=1;next;}{if(zz<1)ok[$1]=1;}{if(ok[$1]==1)next;split($1,aa,"#");split(aa[2],bb,"/");s=1;if(bb[1]>0)t=bb[1];else t=1;}/BadScore/{bs++;bt+=t;next;}/Partial/{bs++;bt+=t;next;}/OK_Multi/{next}/Multi/{s10++;t10+=t;}END{if(s10>0)printf("At_least_10_sites %d seq %d tags \n", s10, t10) ;if(bt>0)printf("Low_quality_mapping %d seq %d tags \n", bs, bt) ;}' >  tmp/$MAGIC_COUNT_DIR/$lane.badcounts.txt 
  endif
end
cat tmp/$MAGIC_COUNT_DIR/$run/*.badcounts.txt | gawk '{a=$1;s[a]+=$2;t[a]+=$4;}END{for(a in s)printf("%s %d seq %d tags\n",a, s[a],t[a]);}' >> $toto.ace

gunzip -c tmp/PHITS_genome/$run/*.overRepresentedSeeds.gz | gawk '{split($1,aa,"#");split(aa[2],bb,"/");ns++;if(bb[1]>0)nt+=bb[1];else nt++;}END{printf("Diffuse_mapping  %d seq %d tags \n", ns,nt) ;}' >> $toto.ace

gunzip -c tmp/$MAGIC_COUNT_DIR/$run/*.noInsert.gz | gawk '/^#/{next}{if($5){ns++;nt+=$2;bp+=$2*($4-$3+1);}}END{printf("No_insert  %d seq %d tags %d kb \n", ns,nt,bp/1000) ;}' >> $toto.ace

gunzip -c tmp/$MAGIC_COUNT_DIR/$run/*.aliProfile.gz |  gawk -F '\t' '/^#/{next}/^Perfect_reads/{next;}/^Cumul/{next}{z=$1 "\t" $2; nn[z]=1;for(i=3;i<11;i++)n[z,i]+=0+$i;}END{for(k in nn){printf("%s",k);for(i=3;i<11 ;i++)printf("\t%d",n[k,i]);printf("\n");}}' | sort -k 1,1 -k 2,2n >> $toto.ace

gunzip -c tmp/$MAGIC_COUNT_DIR/$run/*.aliProfile.gz |  gawk -F '\t' '/^#/{next}/^Cumul/{z=$1; nn[z]=1;for(i=2;i<10;i++)n[z,i]+=0+$i;}END{for(k in nn){printf("%s",k);for(i=2;i<10 ;i++)printf("\t%d",n[k,i]);printf("\n");}}' | sort -k 1,1 -k 2,2n >> $toto.ace


## construct the pair length profile

cat tmp/$MAGIC_COUNT_DIR/$run/*.pairStats | gawk '/^#/{next;}/^$/{next}/^Orphan/{next;}{t=$1;i=t2i[t];if(i<1){imax++;i=imax;t2i[t]=i;i2t[i]=t}nam[i]=$1;}{nn[i]+=$2;}END{for(i=1;i<=imax;i++) {if(nn[i]>0) printf("%s\t%s\t%s\n", r,nam[i],nn[i]) ;}}' r=$run > tmp/$MAGIC_COUNT_DIR/$run/runPairStats.txt
cat tmp/$MAGIC_COUNT_DIR/$run/*.pairStats | gawk '/^#/{next;}/^$/{next}/^Orphan/{nn[$2]+=$3;next;}END{for (i in nn) if(nn[i]>=0)printf("%s\tOrphans %s\t%s\n", r,i,nn[i]) ;}' r=$run >> tmp/$MAGIC_COUNT_DIR/$run/runPairStats.txt
cat  tmp/$MAGIC_COUNT_DIR/$run/runPairStats.txt | gawk -F '\t' '{printf ("%s %s\n",$2,$3);}' >>  $toto.ace
# exclude the last point, which is the remnant of the distribution
cat tmp/$MAGIC_COUNT_DIR/$run/*.pairStats | gawk '/^#Number/{for(i=2;i<NF;i++)n[i]+=$i;if(NF > imax)imax=NF;}END{printf("Insert_length_in_pairs\t%s",r);for(i=4;i<=imax;i++)printf("\t%d",n[i]);printf("\n");}' r=$run >  tmp/$MAGIC_COUNT_DIR/$run/runPairHisto.txt

echo "\nAli $run" >> $toto.ace
if (-e  tmp/$MAGIC_COUNT_DIR/$run/runPairHisto.txt) then
  foreach run2 (`cat MetaDB/$MAGIC/RunPairedList`)
    if ($run == $run2) then
      echo "$run ZZZZZ ZZZZZ" > $toto.fl
      cat $toto.fl ZZZZZ tmp/$MAGIC_COUNT_DIR/$run/runPairHisto.txt |   gawk -f scripts/c5.pair_length_histo.awk | head -100 | gawk -F '\t' '/^Run/{run=$2;next;}/^Average/{av=$2;next;}/^Median/{median=$2;next;}/^Position/{mode=$2;next;}/^Low 1/{z1=$2;next;}/^Low 5/{z5=$2;next;}/^High 95/{z95=$2;next;}/^High 99/{z99=$2;next;}END{printf("\nAli %s\nFragment_length_1 %s\nFragment_length_5 %s\nFragment_length_mode %s\nFragment_length_median %s\nFragment_length_average %s\nFragment_length_95 %s\nFragment_length_99 %s\n\n",run,z1,z5,mode,median,av,z95,z99);}' >>  $toto.ace
      \rm $toto.fl
      touch  tmp/$MAGIC_COUNT_DIR/$run/Fragment_length.fixed
    endif
  end
endif

#######  collate the partials

 cat tmp/COUNT/$run/*.partial.p3p5 | gawk '{p5 += $2; p5b += $3; if ($5=="3p"){p3 +=$6;p3b+=$7;}else{p3+=$5;p3b+=$6;}}END{printf("\nAli %s\nPartial_5p %d %d\nPartial_3p %d %d\n\n",run,p5,p5b,p3,p3b);}' run=$run  >> $toto.ace

#######  finalize

echo ' ' >> $toto.ace


cat $toto.ace | gawk -f scripts/c2.left_blank.awk run=$run > $toto.left_blank.ace
cat  $toto.left_blank.ace >> $toto.ace

echo done
mv $toto.ace  tmp/$MAGIC_COUNT_DIR/$run/c2.alistats.ace
\rm $toto.*

exit 0


\rm toto.ace
foreach run (`cat MetaDB/$MAGIC/RunList`)
  cat tmp/COUNT/$run/*.partial.p3p5 | gawk '{p5 += $2; p5b += $3; if ($5=="3p"){p3 +=$6;p3b+=$7;}else{p3+=$5;p3b+=$6;}}END{printf("\nAli %s\nPartial_5p %d %d\nPartial_3p %d %d\n\n",run,p5,p5b,p3,p3b);}' run=$run  >> toto.ace
end
