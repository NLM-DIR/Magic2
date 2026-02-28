#!bin/tcsh -f
# attention  -ef kills the code in gunzip -c $mytmp/$run/*/$f1$zone.hits.$uu.gz  if a single lane lacks a single target
set run=$1
set justMito=$2
set geneSNPTarget=$3

setenv ici `pwd`

# we cannot ventilate in tmp space and parallelize the sections
set submitLocally="32G"
if (1) then
  set mytmp=$TMPDIR/aceview.s1.snp_accumulate.$$
  set submitLocally="local"
endif

# set mytmp=/export/home/TMP/aceview.s1.snp_accumulate.472
echo $mytmp

if (! -d $mytmp) mkdir $mytmp
if (! -d $mytmp/$run) mkdir $mytmp/$run

echo ZZZZZ | gzip > $mytmp/ZZZZZ.gz

        set qual=""
        if ($run == DAUER.FQ) set qual=--qual
	set qual2=""

        if (-e tmp/SNP_BRS/$run/Quality_profile.txt) then
          set n=`wc tmp/SNP_BRS/$run/Quality_profile.txt | gawk '{print 0+$1}'`
          if ($n > -10) then
	    set qual2=" --runQuality tmp/SNP_BRS/$run/Quality_profile.txt"
          else
            if ($justMito == 0) then
              echo "empty file  tmp/SNP_BRS/$run/Quality_profile.txt, I cannot compute the SNP BRS table"
              touch  tmp/SNP_BRS/$run/Quality_profile.txt
            endif
          endif
        else
          if ($justMito == 0) then
            echo "missing file  tmp/SNP_BRS/$run/Quality_profile.txt, I cannot compute the SNP BRS table"
            touch  tmp/SNP_BRS/$run/Quality_profile.txt
          endif
        endif

setenv ZONE ZONE

set dropMultiplicity=""
if ($?snpDropMultiplicity == 1) set dropMultiplicity="-dropMultiplicity"

if ($justMito == 1) goto laba
if ($NO_INTRON == 0 && $Strategy == RNA_seq) goto phaseRNA

  if (! -e tmp/SNP_ZONE/_allG ) then
    set strategyOk=1
    if (-d tmp/SNP_ZONERUN/$run) \rm  -rf tmp/SNP_ZONERUN/$run 
    if (! -d  tmp/SNP_ZONERUN/$run) mkdir tmp/SNP_ZONERUN/$run 
    setenv ZONE ZONERUN/$run
    if (! -e tmp/SNP_$ZONE/_allg) then
      set runS=$run
      if (-e MetaDB/$MAGIC/r2sublib) then
        set runS=`cat MetaDB/$MAGIC/r2sublib | gawk '{if($2==run)runS=$1;}END{if(runS)print runS;else print run;}' run=$run`
      endif
      foreach chrom ($chromSetAll)
        if (! -e tmp/WIGGLERUN/$runS/$chrom/coverome.$minCoveron.u.peaks) then
          echo "missing file  tmp/WIGGLERUN/$runS/$chrom/coverome.$minCoveron.u.peaks"
          continue
        endif
        cat tmp/WIGGLERUN/$runS/$chrom/coverome.$minCoveron.u.peaks | gawk '/^#/{next;}{print}' | cut -f 1,2,3 | sort -k 1,1 -k 2,2n | gawk -F '\t' '{a1=$2;a2=$3;if($1 != old)olda2=0;if(a1+0>10+olda2)a1-=10;a2+=10;olda2=a2;old=$1;printf("%s\t%09d\t%09d\n",$1,a1,a2);}' >> tmp/SNP_$ZONE/_allg
      end

      set n=`cat tmp/SNP_$ZONE/_allg | gawk '/^#/{next}{dx=$3-$2+1;n+=dx;}END{print n;}'`
      cat  tmp/SNP_$ZONE/_allg | gawk '/^#/{next}{bp+=$3-$2+1;if($1 != oldc || bp > 2500000) { bp=$3-$2+1 ; k++;} oldc=$1; if (k!=oldk) { out="tmp/SNP/" ZONE "/zoneg." k ".txt" ;} oldk=k ;printf("%s\t%s\t%s\t%d\n",$1,$2,$3,k) > out ; }'  nn=$n  ZONE=$ZONE

      touch  tmp/SNP_$ZONE/chom2zone.dummy
      \rm tmp/SNP_$ZONE/chom2zone.*
      foreach  ff (`ls tmp/SNP_$ZONE/zoneg.*.txt `)
        set n=`echo $ff | gawk '{i=index($1,"zoneg.")+6;z=substr($1,i);gsub(/\.txt/,"",z);print z+0}'`
        set chrom=`cat $ff | gawk '/^#/{next;}{print $1;exit;}'`
        ln -s $ici/TARGET/CHROMS/$species.chrom_$chrom.fasta.gz tmp/SNP_$ZONE/zoneg.$n.fasta.gz 
        echo $n >> tmp/SNP_$ZONE/chom2zone.$chrom
      end

    endif
  endif

if (! -e tmp/SNP_ZONE/_allG && -d tmp/SNP_ZONERUN) then
  foreach ff ( `ls tmp/SNP_ZONERUN/$run/zoneg.*.txt` )
    set n=`cat $ff | cut -f 1 | sort -u | wc -l`
    if ($n > 1) then
      echo "ERROR $n (> 1)  chromosomes in file  $ff"
      exit 1
    endif
  end
endif

echo -n "snp accumulate ventilate in directory $mytmp/$run "
date

# tmp/COUNT/$lane.hits.gz 

if (-e tmp/SNP_$ZONE/_allg && ! -e $mytmp/$run/all_zoneg) then
  cat tmp/SNP_$ZONE/zoneg.*.txt | sort -k 1,1 -k 2,2n > $mytmp/$run/all_zoneg
endif


if (-e tmp/SNP_ZONE/_allG) then
  set ff1=zoneG.
  mkdir $mytmp/$run
  foreach lane (`cat Fastc/$run/LaneList`)
    mkdir $mytmp/$lane
    echo "bin/snp -minAliPerCent 90 --ventilate --run $run -o $mytmp/$lane/zoneG -i tmp/COUNT/$lane.hits.gz  --select tmp/SNP_ZONE/_allG $dropMultiplicity"
          bin/snp -minAliPerCent 90 --ventilate --run $run -o $mytmp/$lane/zoneG -i tmp/COUNT/$lane.hits.gz  --select tmp/SNP_ZONE/_allG $dropMultiplicity
    ls -ls  $mytmp/$lane/zoneG*.hits*
  end
endif
if (-e tmp/SNP_ZONE/_allg) then
  set ff1=zoneg.
  foreach lane (`cat Fastc/$run/LaneList`)
    mkdir $mytmp/$lane
    echo "bin/snp -minAliPerCent 90 --ventilate --run $run -o $mytmp/$lane/zoneg -i tmp/COUNT/$lane.hits.gz  --select $mytmp/$run/all_zoneg "
          bin/snp -minAliPerCent 90 --ventilate --run $run -o $mytmp/$lane/zoneg -i tmp/COUNT/$lane.hits.gz  --select $mytmp/$run/all_zoneg 
  end
endif
goto laba

phaseRNA:

set target=$geneSNPTarget
source scripts/target2target_class.txt  

if (-e tmp/SNP_$ZONE/_allr && ! -e $mytmp/$run/all_zoner) then
  cat tmp/SNP_$ZONE/zoner.*.txt | sort -k 1,1 -k 2,2n > $mytmp/$run/all_zoner
endif

if (-e $mytmp/$run/all_zoner) then

  set ff1=zoner.
  foreach lane (`cat Fastc/$run/LaneList`)
    if (-e  $mytmp/$lane/zoner.1.hits.u.gz) continue 
    mkdir $mytmp/$lane
    echo "bin/snp --minAliPerCent 90 --ventilate --run $run -o $mytmp/$lane/zoner -i tmp/COUNT/$lane.hits.gz  --select $mytmp/$run/all_zoner "
          bin/snp --minAliPerCent 90 --ventilate --run $run -o $mytmp/$lane/zoner -i tmp/COUNT/$lane.hits.gz  --select $mytmp/$run/all_zoner 
  end

endif

###############################################################
laba:
echo -n "snp accumulate: ventilation done "
date
echo "start count"

if ($NO_INTRON == 0 && $Strategy == RNA_seq) then
  if (-e tmp/SNP_$ZONE/_allr) then
    set ff1=zoner.
    cat tmp/SNP_ZONE/ZoneList > $mytmp/$run/ZoneList
  endif
else
  if (-e tmp/SNP_$ZONE/_allG) then
    set ff1=zoneG.
    cat tmp/SNP_ZONE/ZoneList > $mytmp/$run/ZoneList
  endif
  if (-e tmp/SNP_$ZONE/_allg)  then
    set ff1=zoneg.
    ls $mytmp/$run/*/zoneg.*.hits.u | gawk '{n=split($1,aa,"/");split(aa[n],bb,".");printf ("%s.%s\n",bb[1],bb[2]);}' | sort -u  >> $mytmp/$run/ZoneList
  endif
endif

# SpikeIn rrna
foreach zone (`cat $mytmp/$run/ZoneList` )
  set target=xx
echo "hello from zone $zone strategy=$Strategy"


  if ($justMito == 1 && $zone != mito && $zone != SpikeIn) continue
#  if (($zone == mito || $zone == rrna || $zone == SpikeIn) && ! -d tmp/PHITS_$zone) continue
  if ($zone == mito || $zone == rrna || $zone == SpikeIn) set target=$zone
  if ($zone =~ zoneG.* && $Strategy != Exome && $Strategy != Genome && $NO_INTRON == 0) continue
  if ($zone =~ zoneg.* && $Strategy != Exome && $Strategy != Genome  && $NO_INTRON == 0) continue
  if ($zone =~ zoneG.*) set target=genome
  if ($zone =~ zoneg.*) set target=genome
  if ($zone =~ zoner.* && $Strategy != RNA_seq) continue
  if ($zone =~ zoner.* && $Strategy == RNA_seq) set target=$geneSNPTarget

  if ($target == xx) continue
  source scripts/target2target_class.txt  

       if  (! -e  tmp/SNP_BRS/$run/$zone.BRS.u.gz) then

        set slct="--fasta TARGET/Targets/$species.$target.fasta.gz"
        if ($target == mito) then
          set slct=" --mito "
        endif

        # if ($zone != $target)  set slct="--select tmp/SNP_$ZONE/$zone.txt --fasta tmp/SNP_$ZONE/$zone.fasta.gz"
        if ($zone != $target)  set slct=" --fasta tmp/SNP_$ZONE/$zone.fasta.gz"
        set f1=""
        if ($zone == $target) set f1="$ff1"
        if ($zone == $target) set slct=" --fasta TARGET/Targets/$species.$target.fasta.gz "
        set ff="--strategy $Strategy --minAliPerCent 90 --target_class $target_class "

        set rr="--run $run"

echo hello3 $target $zone $rr
           foreach uu (u)
             set unique=""
             if ($uu == u) set unique=" --unique "
             # set qual2=""
	     if ($justMito == 1) then
               set n=`ls $mytmp/ZZZZZ.gz tmp/COUNT/$run/*.mito.gz | wc | gawk '{print $1}'`
               if ($n < 2) continue
               if (-e $mytmp/$run/$zone.BRS.$uu.gz) continue
               echo "scripts/submit $mytmp/$run/$zone.BRS.$uu gunzip -c tmp/COUNT/$run/*.mito.gz  | bin/snp $slct $qual --hits2BRS $ff $rr  --gzo  -o $mytmp/$run/$zone.BRS.$uu"
               # if (-e tmp/SNP/$run/mito.BRS.$uu.txt.gz)    \rm tmp/SNP/$run/mito.BRS.$uu.txt.gz 
               # if (-e tmp/SNP/$run/mito.detect.$uu.snp.gz) \rm tmp/SNP/$run/mito.detect.$uu.snp.gz

               scripts/submit $mytmp/$run/$zone.BRS.$target.$uu "gunzip -c  tmp/COUNT/$run/*.mito.gz | bin/snp $slct $qual $qual2  $unique --hits2BRS $ff $rr  --gzo  -o $mytmp/$run/$zone.BRS.$uu" $submitLocally
             else
               echo "hello5 $f1$zone.hits.$uu"
               echo "ls -ls $mytmp/ZZZZZ.gz $mytmp/$run/*/$f1$zone.hits.$uu"
               set n=`ls $mytmp/ZZZZZ.gz $mytmp/$run/*/$f1$zone.hits.$uu | wc | gawk '{print $1}'`
               echo "n=$n"
               if ($n < 2) continue
               if (-e $mytmp/$run/$zone.BRS.$uu.gz) continue
               echo   "cat $mytmp/$run/*/$f1$zone.hits.$uu | bin/snp $slct $qual $qual2 $unique --hits2BRS $ff $rr  --gzo  -o $mytmp/$run/$zone.BRS.$uu"
               cat $mytmp/$run/*/$f1$zone.hits.$uu | bin/snp $slct $qual $qual2 $unique --hits2BRS $ff $rr  --gzo  -o $mytmp/$run/$zone.BRS.$uu
             endif
           end
        
        # mv $mytmp/$run/$zone.BRS.*.gz  tmp/SNP/$run/
      endif
end

# avoid double counting
 # $mytmp/$run/ZoneList
if (0 && -e tmp/SNP_$ZONE/_allG) then

  foreach zone (`cat $mytmp/$run/ZoneList`)

    set z=`echo $zone | sed -e 's/zoneG.//'`
    set uu=`cat tmp/SNP_ZONE/_allG | gawk -F '\t' '{if($4 == z){c=$1;a1=$2;a2=$3;ok=1;next;}if(ok==1 && c==$1 && $2<a2){ok=2;a11=$2;z2=$4;}}END{if(ok==2)printf("%s#%d#zoneG.%s\n",z,a11,z2);else print z;}' z=$z`
    echo "uu=$uu"
    set a11=`echo $uu | gawk '{split($1,aa,"#");print aa[2]+0;}'`
    if ($a11 == 0) continue 
    set zone2=`echo $uu | gawk '{split($1,aa,"#");print  aa[3];}'`
    echo "$zone $zone2 split at $a11"

    set uu=u
    gunzip -c $mytmp/$run/$zone.BRS.$uu.gz | gawk -F '\t' '/^#/{print > f2".1";next;}{if($2+0<a11 || $1 != "~")print > f2".1"; else print > f2".2"}' a11=$a11 f2=$mytmp/$run/$zone.BRS.x
    mv $f2.1  $mytmp/$run/$zone.BRS.$uu.new
    \rm  $mytmp/$run/$zone.BRS.$uu.gz 
    gzip $mytmp/$run/$zone.BRS.$uu.new

    set f3=$mytmp/$run/$zone2.BRS.x
    gunzip -c   $mytmp/$run/$zone2.BRS.$uu.gz > $f3.1
    cat $f3.1 | gawk '/^# Target/{print ;exit;}{print}' > $f3.2
    cat $f3.1 | gawk '/^# Target/{ok=1;next;}{if(ok==1)print}' > $f3.3
    cat $f2.2 $f3.3 |  sort -k 2,2n -k 4,4 -k 3,3r > $f3.4
   
    cat $f3.4 ZZZZZ | gawk -F '\t' '{if($2 != old[2] || $3 != old[3] || $4 != old[4]){if(old[2]>0){printf("%s",old[1]);for(i=2;i<=9;i++)printf("\t%s",old[i]);printf("\n");}for(i=1;i<=5;i++)old[i]=$i;for(i=6;i<=9;i++)old[i]=0;}for(i=6;i<=9;i++)old[i]+=$i;}' >> $f3.2
    mv $f3.2  $mytmp/$run/$zone2.BRS.$uu.new2
    gzip  $mytmp/$run/$zone2.BRS.$uu.new2
    \rm $f2.* $f3.*

  end

endif


#############################################################################
## cleanup

echo -n "snp accumulate: count done "
date
 

cleanup:
echo "clean up: $mytmp"

ls -ls  $mytmp/$run/*/*.gz
mv $mytmp/$run/s1.*  $mytmp/$run/*.gz tmp/SNP_BRS/$run
ls -ls  tmp/SNP_BRS/$run/*.gz

\rm -rf $mytmp

if ($justMito == 1) then
  touch tmp/SNP_BRS/$run/s1.justmito.done
else
  touch tmp/SNP_BRS/$run/s1.done
endif

echo -n "snp accumulate done "
date

##### wiggle a la base pres des petits targets
## gunzip -c tmp/SNP/Rhs434/SpikeIn.BRS.u.gz | gawk '/^ERR/{next}{if($1 != "~")s=$1;if(length($3)>1)next;printf("%s\t%d\t%s\t%s\t%d\n",s,$2,$3,$4,$6);}' > tata4

### grep Exit tmp/SNP/*/s1.err | grep -v '= 0' | gawk '{i=index($1,"s1");printf("\\rm -rf tmp/SNP*/%s\n",substr($1,9,i-2-8)); }' > _k
