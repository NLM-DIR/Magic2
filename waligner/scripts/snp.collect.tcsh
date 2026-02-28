#!bin/tcsh -f

set run=$1
set justMito=$2
set poolDummy=$3
set solidDummy=$4
set collect=$5
set zone=$6

set uu=u

set minSnpFrequency2=$minSnpFrequency
if (-e tmp/SNP_BRS/$MAGIC.minFrequency.txt) then
  # reset minSnpFrequency to at least twice the percent error rate, the error_profile is given in error per million 
  set z=`cat tmp/SNP_BRS/$MAGIC.minFrequency.txt | gawk '{if ($1 == run) z=($4 + 0)/5000 ; if (mf<z)mf=z;}END{print int(mf);}' mf=$minSnpFrequency run=$run`
  set minSnpFrequency2=$z
endif

echo "scripts/snp.collect.tcsh $1 $2 $3 $4 $5 $6 minSnpCover=$minSnpCover minSnpCount=$minSnpCount minFrequency=$minSnpFrequency2"

# avoid analysing MetaDB/*List in submitted scripts
set solid=""
if ($solidDummy == solid) set solid="-solid"

# echo the command, then parse the BRS table and export a list of snp counts in target coordinates detailling both strands
# the program collects the snps not wild type at risk 1/100, function snpNotLow() 
# and all the counts at the positions listed in -snp_list even if in the current run we see a wild type

echo "set vdb"
set vdb=""
if ($collect == count) then 
    set vdb="toto"
    if (-e  tmp/SNP_LIST/$zone/$MAGIC.Variant.list) set vdb="tmp/SNP_LIST/$zone/$MAGIC.Variant.list"
    if (-e  Global_SNP_LIST/$MAGIC.$zone.Variant.list) set vdb="Global_SNP_LIST/$MAGIC.$zone.Variant.list"
    if ($vdb == toto) exit (1)
    set vdb="-snp_list $vdb"
endif

####  Main work
echo -n "$collect main work $vdb : "

date

    set out=tmp/SNP/$run/$zone.$uu
    if ($collect == count) set out=tmp/SNP/$run/$MAGIC.$zone.$uu

    set select8kb=""
    set target=`echo $Etargets av | gawk '{print $1;}'`
    if ($collect == detect && -e tmp/METADATA/$target.selected8kbTranscriptList.txt) then
      set select8kb="-selected8kbList  tmp/METADATA/$target.selected8kbTranscriptList.txt"
    endif
    if (! -e $out.gz) then
  
      # in the new method, we add the BRS of the runs into the groups, so we can use the high thresholds at the detect stage
      set mins=" -minCover $minSnpCover -minMutant $minSnpCount -minFrequency $minSnpFrequency2"
      echo "bin/snp -BRS_$collect $solid $mins -run $run $vdb -strategy $Strategy  $select8kb -i  tmp/SNP_BRS/$run/$zone.BRS.$uu.gz  -o $out -gzo" 
            bin/snp -BRS_$collect $solid $mins -run $run $vdb -strategy $Strategy  $select8kb -i  tmp/SNP_BRS/$run/$zone.BRS.$uu.gz  -o $out -gzo 

    endif
    ls -ls $out.gz
######

    if ($collect == detect) then
      echo "Creating   tmp/SNP_LIST/$zone/Variant.$run.$uu.list"
      ls -ls    tmp/SNP/$run/$zone.$uu.snp_detect.tsf.gz
      gunzip -c tmp/SNP/$run/$zone.$uu.snp_detect.tsf.gz |  gawk '/^#/{next}{print $1}' | sort -u >   tmp/SNP_LIST/$zone/Variant.$run.$uu.list
     wc  tmp/SNP_LIST/$zone/Variant.$run.$uu.list
    endif

stats:

echo -n "$collect stats : "
date

set zone1="$zone"
if ($collect == count) set zone1="$MAGIC.$zone"

    gunzip -c  tmp/SNP/$run/$zone1.$collect.$uu.snp.gz  |  gawk -F '\t' '/Incompatible_strands/{next;}{t=$3;c=$9;m=$10;w=$11;if(c>minC && 100*m>c*minF){nn++;n[t]++;nnm+=m;mm[t]+=m;ww[t]+=w;nnw+=w;}}END{printf("%s\tTotal\t%d\t%d\t%d\t%d\t%d\n",run,minF,minC,nn,nnm,nnw);for(t in n)printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\n",run,t,minF,minC,n[t],mm[t],ww[t]);}' run=$run minC=$minSnpCover minF=$minSnpFrequency2 | sort -k 5n >   tmp/SNP/$run/$zone1.$collect.$uu.stats
    gunzip -c  tmp/SNP/$run/$zone1.$collect.$uu.snp.gz  |  gawk -F '\t' '/Incompatible_strands/{next;}{t=$3;c=$9;m=$10;w=$11;if(c>minC && 100*m>c*minF){nn++;n[t]++;nnm+=m;mm[t]+=m;ww[t]+=w;nnw+=w;}}END{printf("%s\tTotal\t%d\t%d\t%d\t%d\t%d\n",run,minF,minC,nn,nnm,nnw);for(t in n)printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\n",run,t,minF,minC,n[t],mm[t],ww[t]);}'  run=$run minC=$minSnpCover minF=95 | sort -k 5n >>   tmp/SNP/$run/$zone1.$collect.$uu.stats


    gunzip -c  tmp/SNP/$run/$zone1.$collect.$uu.snp.gz | gawk -F '\t' '/Incompatible_strands/{next;}{type=$3;if($10 + $11 > 0){split(type,aa,":");if(aa[1]=="Sub")type=aa[2] ">" aa[3];else{split(aa[1],bb,"_");k=0;if (bb[1]=="Ins")k=3;if (bb[1]=="Del")k=2;if(k==0)next;z=aa[k];type=bb[1];for(i=1;i<=length(z);i++){c=substr(z,i,1);if(c==toupper(c))type=type c;}}n[type]++;m[type]+=$10;w[type]+=$11;}}END{for(t in n)printf("%s\t%s\t%d\t%d\t%d\n",t,run,n[t],m[t],w[t]);}' run=$run | sort >  tmp/SNP/$run/$zone1.$collect.types.$uu.txt
  
  if ($collect == count) then
    gunzip -c  tmp/SNP/$run/$MAGIC.*.$collect.$uu.snp.gz | bin/histo -snp -o  tmp/SNP/$run/$zone1.snp_frequency
  endif

if ($justMito == 1) then
  touch tmp/SNP/$run/$zone1.snp.$collect.justmito.done
else
  touch tmp/SNP/$run/$zone1.snp.$collect.done
endif

echo -n "$collect done : "
date
echo done

exit 0

# HACK to back fix the new nomenclature in the .types.file
set uu=u
foreach collect (detect count)
  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    foreach zone (`cat tmp/SNP_ZONE/ZoneList`)
      set zone1=$zone
      if ($collect == count) set zone1="$MAGIC.$zone"
      gunzip -c  tmp/SNP/$run/$zone1.$collect.$uu.snp.gz | gawk -F '\t' '/Incompatible_strands/{next;}{type=$3;if($10 + $11 > 0){split(type,aa,":");if(aa[1]=="Sub")type=aa[2] ">" aa[3];else{split(aa[1],bb,"_");k=0;if (bb[1]=="Ins")k=3;if (bb[1]=="Del")k=2;if(k==0)next;z=aa[k];type=bb[1];for(i=1;i<=length(z);i++){c=substr(z,i,1);if(c==toupper(c))type=type c;}}n[type]++;m[type]+=$10;w[type]+=$11;}}END{for(t in n)printf("%s\t%s\t%d\t%d\t%d\n",t,run,n[t],m[t],w[t]);}' run=$run | sort >  tmp/SNP/$run/$zone1.$collect.types.$uu.txt
    end
  end
end

##########################################################
## verif scripts

foreach run ( `cat MetaDB/$MAGIC/RunList` )
  if (! -e tmp/SNP/$run/zoneg.1.detect.u.snp.gz) then
    echo tmp/SNP/$run
  endif
end


foreach run ( `cat MetaDB/$MAGIC/RunList` )
  if (! -e tmp/SNP/$run/zoneg.1.$collect.u.snp.gz) then
    echo tmp/SNP/$run
  else
    touch tmp/SNP/$run/snp.$collect.done
  endif
end

foreach run ( `cat MetaDB/$MAGIC/RunList` )
  if (! -e tmp/SNP/$run/s1.err) echo tmp/SNP/$run
  grep Exit tmp/SNP/$run/s1.err
end

# grep Exit tmp/SNP/*/s1.err | grep -v '= 0' | gawk '{split($1,aa,"/");printf("\\rm -rf tmp/SNP*/%s\n",aa[3]);}' > _rm

exit 0

# hack to look for saturation
\rm _toto
foreach ii (0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15)
  cat B$ii/Global.B$ii.Variant.list >> _toto
end
cat _toto | sort -u > Global.B0_15.Variant.list

# by hank i run the scripts 15 times stopping at say 15 .. 8 7 .. 0
\rm _toto
foreach ii (0 1   )
  cat B$ii/Global.B$ii.Variant.list >> _toto
end
cat _toto | sort -u > Global.B0_1.Variant.list

