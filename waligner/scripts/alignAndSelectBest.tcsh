#!bin/tcsh -f

set run=$1
set lane=$2
set minInsertLength=$3
set useMagicBlastTag=$4

date
echo "alignAndSelectBest.tcsh $1 $2 $3 $4"

set mIL=""
if ( $minInsertLength > 0) set mIL="-Remove_inserts_shorter_than $minInsertLength"
set COUNT=COUNT

set lane_clean=`echo $lane | sed -e 's/\//_/g'`
set mytmp=$TMPDIR/aceview.a123.$lane_clean.$$
set cleanUp=1
if (0) then
  set mytmp=/export/home/TMP/aceview.a123.1
  set mytmp=tmp/TMP/$MAGIC
  mkdir tmp/TMP 
  set cleanUp=0
else
  if (-d $mytmp) \rm -rf $mytmp
endif
echo "$mytmp"

  mkdir $mytmp
  mkdir $mytmp/MetaDB
  mkdir $mytmp/MetaDB/$MAGIC
  mkdir $mytmp/COUNT
  mkdir $mytmp/COUNT/$run

if (! $?targets) then
  setenv targets "$DNAtargets $RNAtargets"
endif
if ($NO_INTRON == 2) then
  setenv targets "$DNAtargets"
endif

echo "$mytmp\ntarget=$targets"
\cp MetaDB/$MAGIC/*  $mytmp/MetaDB/$MAGIC

  set v11=""
  set v22=""
  set adaptor1=""
  set adaptor2=""

    set plateforme=`cat $mytmp/MetaDB/$MAGIC/run2machine2adaptors.txt | gawk 'BEGIN{m="NULL";}{if($1 == run && $2 != machine) m=$2;}END{print m;}' run=$run`
    if ($plateforme == "Illumina" && $?exitAdaptorIllumina) then
      set v11="$exitAdaptorIllumina"
    endif
    if ($plateforme == "Roche_454" && $?exitAdaptorRoche_454) then
      set v11="$exitAdaptorRoche_454"
    endif
    if ($plateforme == "SOLiD" && $?exitAdaptorSOLiD) then 
      set v11="$exitAdaptorSOLiD"
    endif
    if ($plateforme == "PacBio" && $?exitAdaptorPacBio) then 
      set v11="$exitAdaptorPacBio"
    endif
    if ($plateforme == "PGM" && $?exitAdaptorPGM) then 
      set v11="$exitAdaptorPGM"
    endif
    if ($plateforme == "Helicos" && $?exitAdaptorHelicos) then
      set v11="$exitAdaptorHelicos"
    endif

set adaptor1=`cat $mytmp/MetaDB/$MAGIC/run2machine2adaptors.txt | gawk -F '\t' 'BEGIN{m="NULL";}{gsub(/\"/,"",$1);if($1 == run){if ($5 != "NULL") m=$5;if ($3 != "NULL") m=$3;}}END{m=substr(m,1,30);gsub(/\"/,"",m);gsub("Text:","",m);print m;}' run=$run`
set adaptor2=`cat $mytmp/MetaDB/$MAGIC/run2machine2adaptors.txt | gawk -F '\t' 'BEGIN{m="NULL";}{gsub(/\"/,"",$1);if($1 == run){if ($6 != "NULL") m=$6;if ($4 != "NULL") m=$4;}}END{m=substr(m,1,30);gsub(/\"/,"",m);gsub("Text:","",m);print m;}' run=$run`

set entryAdaptor1=`cat $mytmp/MetaDB/$MAGIC/run2machine2adaptors.txt | gawk -F '\t' 'BEGIN{m="NULL";}{gsub(/\"/,"",$1);if($1 == run){if ($7 != "NULL") m=$7;}}END{gsub(/\"/,"",m);gsub("Text:","",m);print m;}' run=$run`
set entryAdaptor2=`cat $mytmp/MetaDB/$MAGIC/run2machine2adaptors.txt | gawk -F '\t' 'BEGIN{m="NULL";}{gsub(/\"/,"",$1);if($1 == run){if ($7 != "NULL") m=$8;}}END{gsub(/\"/,"",m);gsub("Text:","",m);print m;}' run=$run`

  echo "# adaptor1=$adaptor1"
  echo "# adaptor2=$adaptor2"
  echo "# entryAdaptor1=$entryAdaptor1"
  echo "# entryAdaptor2=$entryAdaptor2"

set v1=""
set v2=""
if ("$adaptor1" != "" && "$adaptor1" != "NULL") then
  set v1="-exitAdaptor $adaptor1"
else
  if ("$v11" != "" ) then
    set v1="-exitAdaptor $v11"
  endif
endif
if ("$entryAdaptor1" != "" && "$entryAdaptor1" != "NULL") then
  set v1="$v1 -entryAdaptor $entryAdaptor1"
endif
if ("$entryAdaptor2" != "" && "$entryAdaptor2" != "NULL") then
  set v1="$v1 -entryAdaptor2 $entryAdaptor2"
endif

if ("$adaptor2" != "" && "$adaptor2" != "NULL") then
  set v2="-exitAdaptor2 $adaptor2"
else
  if ("$v22" != "" ) then
    set v2="-exitAdaptor2 $v22"
  endif
endif

echo "vectors:: $v1 $v2"


phasea1:
 echo -n "Phase a123: $run $lane  align on $targets and select best alignment "
 echo "Using tmp directory: $mytmp"
 date

set pair=""
foreach run2 (`cat $mytmp/MetaDB/$MAGIC/RunPairedList`)
  if ($run == $run2) then
    set pair="-pair 500"
  endif
end

# 2020_04_03: isMir kills -MRNAH whihc is annoying and small reads are aligned in a different way, nit using the present script
set isMir=0
foreach run2 (`cat $mytmp/MetaDB/$MAGIC/RunSmallRnaList`)
  if (0 && $run == $run2) then
    set isMir=1
  endif
end

set isSolid=" "
set iiSolid=0
foreach run2 (`cat $mytmp/MetaDB/$MAGIC/RunSolidList`)
  if ($run == $run2) set iiSolid=1
end
if ($iiSolid == 1) set isSolid="-solidC"

set isLong=""
foreach run2 (`cat $mytmp/MetaDB/$MAGIC/RunNanoporeList`)
  if ($run == $run2) set isLong="-nanopore"
end

foreach run2 (`cat $mytmp/MetaDB/$MAGIC/RunPacBioList`)
  if ($run == $run2) set isLong="-pacbio"
end

foreach run2 (`cat $mytmp/MetaDB/$MAGIC/RunUltimaList`)
  if ($run == $run2) set isLong="-ultima"
end

set dna2dnaFormat=fastc
if ($iiSolid == 1) set dna2dnaFormat=csfastc

# it is necessary to compute the genome before the transcriptome
# Yet after SpikeIn mito rrna chloro 
# so that we can apply the previous score limit and still find
# all the exons and hence all the introns

set hitfiles = ""

if (-e tmp/Unaligned/$lane.fastc.gz)  then
  \rm  tmp/Unaligned/$lane.*
endif

printf "dummy\t1\n" | gzip >  $mytmp/COUNT/$lane.best_score.gz
if (-e tmp/COUNT/$lane.hits.gz)  then
  gunzip -c tmp/COUNT/$lane.hits.gz | gawk -F '\t' '/^#/{next;}{p=$1;s=$2;if(p==oldp)next;printf("%s\t%s\n",p,s);}' | gzip > $mytmp/COUNT/$lane.best_score.gz
else
  if (-e tmp/COUNT/$lane.hits.1.gz)  then
    gunzip -c tmp/COUNT/$lane.hits.1.gz | gawk -F '\t' '/^#/{next;}{p=$1;s=$2;if(p==oldp)next;printf("%s\t%s\n",p,s);}' | gzip > $mytmp/COUNT/$lane.best_score.gz
  endif
endif

# impose genome last
set gTarget=0
foreach target (DNASpikeIn SpikeIn mito rrna chloro transposon $RNAtargets $DNAtargets)
  # after chloro we apply previous best score is available
  if ($gTarget == 1) set gTarget=2
  if ($target == chloro && $gTarget == 0) set gTarget=1
  # compute only once
  if (-e $mytmp/PHITS_$target/$lane.hits.gz) then
    continue  
  endif
  if (-e tmp/PHITS_$target/$lane.done && -e tmp/PHITS_$target/$lane.hits.gz) then
    if (! -d $mytmp/PHITS_$target/$run) mkdir $mytmp/PHITS_$target/$run
    \cp tmp/PHITS_$target/$lane.*  $mytmp/PHITS_$target/$run
    goto plusbas
    continue
  endif
  if (-e $mytmp/PHITS_$target/$lane.done || -e tmp/PHITS_$target/$lane.done || -e tmp/PHITS_$target/$lane.hits.gz) continue
  # only compute requested targets
  set ok=0
  foreach target2 ($targets)
    if ($target == $target2) set ok=1
  end
  if ($ok == 0) continue

  echo -n "$run $lane target=$target  "
  date

  source scripts/target2target_class.txt
  set targetFasta=$target
  if ($target == gdecoy) set targetFasta=genome
  if ($target == SpikeIn) then
    foreach run2 (`cat $mytmp/MetaDB/$MAGIC/RunSequinList`)
      if ($run == $run2) then
        set targetFasta=SpikeInSequin
      endif
    end
  endif

  if (! -e TARGET/Targets/$species.$targetFasta.fasta.gz) then
    echo "Missing file  TARGET/Targets/$species.$targetFasta.fasta.gz"
    exit 1 
  endif

  set Xintrons=""
  if ($target == introns) set Xintrons="-exactTargetStart 93 -exactTargetStop 108"

  set off3="1 -seedShift 5"
  if ($target == gdecoy) set off3="3 -seedShift 50 -decoy"
  if ($isMir == 1) then
    set off3="1 -seedShift 2"
    if ($target == gdecoy) set off3="3 -seedShift 5 -decoy"
  endif

  if (! -d $mytmp/PHITS_$target) mkdir $mytmp/PHITS_$target
  if (! -d $mytmp/PHITS_$target/$run) mkdir $mytmp/PHITS_$target/$run

  set isSplice=" "
  set slBonus=" "
  set maxHit2=$maxHit
  if ($target == virus) set maxHit2=30
  if ($target == bacteria) set maxHit2=30

# if the bonus is changed, synchronize in snp.confirmation.tcsh the line
# if($8=="A_mito" || $8 == "B_rrna")s--;
# and synchronize a123M.MagicBlast.tcsh
  set bonus=0
  if ($target == gdecoy) set bonus=2
  if ($target == genome) set bonus=0
  if ($target =~ chrom*) set bonus=0
  if ($target == cloud)  set bonus=-6
  if ($target == snp)    set bonus=-6
  if ($target == introns)    set bonus=-6
  if ($target == rnaGene) set bonus=0
  if ($target == smallRNA) set bonus=0
  if ($target == mito) set bonus=1
  if ($target == chloro) set bonus=1
  if ($target == rrna) set bonus=1
  if ($target == transposon) set bonus=1
  if ($target == av2008) set bonus=-1
  if ($target == RefSeqCurrent) set bonus=-2
  if ($target == Gaj) set bonus=-3
  if ($target == FBK) set bonus=-4
  if ($target == virus) set bonus=-30
  if ($target == bacteria) set bonus=-30
  set targetBonus="-targetBonus $bonus"

    set clipPolyA=" "
    set clipPolyT=" "
    if ($Strategy == RNA_seq) then
      if ($isMir == 0) then
        if ($target_class =~ [A-Z]T_*) then
          set clipPolyA="-SclipPolyA"
          set clipPolyT="-SclipPolyT"
        else
          set clipPolyA="-clipPolyA"
          set clipPolyT="-clipPolyT"
        endif
      endif
    endif
   if ($?NoPolyA) then
     if ($NoPolyA == 1) then
       set clipPolyA=" "
       set clipPolyT=" "
     endif
   endif

  set isStranded=" "
  if (0) then
    # 2014_10_20, we decided to align only in non stranded mode, but carefully count the read mapping as hit and anti hits
    foreach run2 (`cat $mytmp/MetaDB/$MAGIC/RunForwardList`)
      if ($run == $run2 && isMir == 0 && $target_class =~ [A-Z]T_*) set isStranded="-stranded $bonusStrand"
    end

    foreach run2 (`cat $mytmp/MetaDB/$MAGIC/RunReverseList`)
      if ($run == $run2 && $target_class =~ [A-Z]T_*) set isStranded="-stranded -$bonusStrand"
    end
  endif

  set isSplice=" "
  set isMRNAH=" "
  if ($Strategy == RNA_seq) then
    if ($isMir == 0 && $target_class =~ [A-Z]T_* ) then
      set isMRNAH="-MRNAH"
    endif
    if ($target_class =~ [09]_* || $target_class =~ B_* || $target_class =~ [A-Z]T_* || $intronMaxLength == 0 ) then
      set isSplice=" "
    else
      foreach run2 (`cat $mytmp/MetaDB/$MAGIC/RunRnaList`)
        if ($run == $run2) then
          set isSplice="-splice"
        endif
      end
    endif
  endif
  if ($species == worm) set slBonus="-slBonus 5"

  set A2G=" "
  if (0 && $MAGIC == ChinaG) set  A2G="-A2G"
  if (0 && $MAGIC == NBE) set  A2G="-A2G"
  if ($?MAGIC_A2G) then
    if ($MAGIC_A2G == 1 && $iiSolid == 0 && $target_class =~ [A-Z]T_*) set  A2G="-A2G"
  endif

  if ($Strategy != RNA_seq) then
    set isSplice=" "
    set intronMaxLength=1000000000
  endif
  if ($target == genome) then      
    # 2015_09_25, coli, we want to find all rearrangements
    set isSplice="-splice"
  endif
  echo "## $Strategy #$target $run $isSolid $clipPolyT  $isSplice $isMRNAH $slBonus $isStranded $targetBonus"

  set minEntropy4=`echo $minAli | gawk '{printf("%d",int((2.0*$1)/3));}'`

  set t2g=""
  if (-e  TARGET/MRNAS/$species.$target.transcript2gene.txt.gz) set t2g="-target2gene TARGET/MRNAS/$species.$target.transcript2gene.txt.gz"

set jump5=""
set jump5=`cat $mytmp/MetaDB/$MAGIC/runs.ace | gawk '{gsub(/\"/,"",$0);}/^Run/{ok=0;if($2 == run)ok=1;}/^Jump5/{if(ok==1){ if(0+$2>0 && 0+$2 < 6)printf(" -jump5 %d", 0+$2);if(0+$3>0 && 0+$3 < 6)printf(" -jump5_2 %d", 0+$3);}}' run=$run`
set forceRightClip=""
set forceRightClip=`cat $mytmp/MetaDB/$MAGIC/runs.ace | gawk '{gsub(/\"/,"",$0);}/^Run/{ok=0;if($2 == run)ok=1;}/^ForceRightClip/{if(ok==1){ if(0+$2>20)printf(" -forceRightClip %d", 0+$2);if(0+$3>0 && 0+$3 < 6)printf(" -forceRightClip_2 %d", 0+$3);}}' run=$run`


echo "prepare previous score "
date
  set prevScore=""
  set toto=""
  if  ($gTarget >= 2 && -e tmp/COUNT/$lane.best_score.gz) then
    # merge current best score with previous round best score
    if (-e $mytmp/COUNT/$lane.best_score.gz) then
      mv  $mytmp/COUNT/$lane.best_score.gz  $mytmp/COUNT/$lane.best_score.old.gz
      gunzip -c tmp/COUNT/$lane.best_score.gz $mytmp/COUNT/$lane.best_score.old.gz  |  gawk -F '\t' '/^#/{next;}{p=$1;s=$2;if (ss[p]<s)ss[p]=s;}END{for(p in ss)printf("%s\t%s\n",p,ss[p]);}' | gzip >  $mytmp/COUNT/$lane.best_score.gz 
      \rm  $mytmp/COUNT/$lane.best_score.old.gz 
    else
      cp tmp/COUNT/$lane.best_score.gz  $mytmp/COUNT/$lane.best_score.gz
    endif
  endif
  if (-e $mytmp/COUNT/$lane.best_score.gz) then
    set prevScore=" -previousScore  $mytmp/COUNT/$lane.best_score.gz "
    \cp  $mytmp/COUNT/$lane.best_score.gz tmp/PHITS_$target/$run
  endif

set avoidPseudo=""
if ($target == av) set avoidPseudo="-avoidPseudoGenes"
echo "align "
date
# if ($target == bacteria)set prevScore=""
 
set overh='-showOverhang -showTargetPrefix'
if ($target == bacteria || $target == virus)set overh=""
if ($VIRUS_PROJECT == 1 || $target == virus) set overh="-showSequence"

set SAM=""
if ($MAGIC_SAM == 1) set SAM="-sam"

set targetMask=""
if (-e TARGET/Targets/$species.$target.mask.txt) set targetMask="-targetMask TARGET/Targets/$species.$target.mask.txt"
echo "TARGET\t$target" >> $mytmp/PHITS_$target/$lane.err

         echo "bin/clipalign -best -i Fastc/$lane.$dna2dnaFormat.gz  $isSolid $isLong -t TARGET/Targets/$species.$targetFasta.fasta.gz  -maxHit $maxHit2 $clipPolyA  $clipPolyT -minEntropy $minEntropy4 -seedLength $seedLength -probeMinLength $minAli  -clipN $clipN -minAli $minAli $isSplice $A2G $isMRNAH $slBonus  $targetBonus  $Xintrons $isStranded -seedOffset $off3 -intronMaxLength  $intronMaxLength -target_class $target_class $t2g $jump5 $forceRightClip $overh -strategy $Strategy $targetMask $avoidPseudo $v1 $v2 $prevScore -gzo -o $mytmp/PHITS_$target/$lane $SAM"

  (bin/time -p bin/clipalign -best -i Fastc/$lane.$dna2dnaFormat.gz  $isSolid  $isLong -t TARGET/Targets/$species.$targetFasta.fasta.gz  -maxHit $maxHit2 $clipPolyA  $clipPolyT -minEntropy $minEntropy4 -seedLength $seedLength -probeMinLength $minAli  -clipN $clipN -minAli $minAli $isSplice $A2G $isMRNAH $slBonus  $targetBonus  $Xintrons $isStranded -seedOffset $off3 -intronMaxLength  $intronMaxLength  -target_class $target_class $t2g $jump5 $forceRightClip $overh -strategy $Strategy $targetMask $avoidPseudo $v1 $v2 $prevScore -gzo -o $mytmp/PHITS_$target/$lane $SAM) >>& $mytmp/PHITS_$target/$lane.err

 if ($status > 0) then
    set sss=$status
    date
    echo "FATAL ERROR bin/clipalign $target exited with non-zero status $sss"
    if (! -d tmp/PHITS_$target/$run) mkdir tmp/PHITS_$target/$run
    \mv $mytmp/PHITS_$target/$run/* tmp/PHITS_$target/$run
    ls -ls core.* | tail -1
    exit 1
 endif

plusbas:
  echo -n  "clipalign $target done : "
  ls -ls core.*
  date   
  ls -ls  $mytmp/PHITS_$target/$lane.hits.gz
  if (-e tmp/$COUNT/$run/c2.alistats.ace) \rm  tmp/$COUNT/$run/c2.*
  if (-e tmp/$COUNT/$run/Venn.$run.txt) \rm  tmp/$COUNT/$run/Venn.$run.txt

# do not export 'noise' detection of pure vectors
  if (-e $mytmp/PHITS_$target/$lane.noInsert.gz) then
    set n=`gunzip -c $mytmp/PHITS_$target/$lane.noInsert.gz | wc -l`
    if ($n < 200) \rm $mytmp/PHITS_$target/$lane.noInsert.gz
  endif

  if (-e $mytmp/PHITS_$target/$lane.hits.gz) then
    set hitfiles="$hitfiles $mytmp/PHITS_$target/$lane.hits.gz"
    mv  $mytmp/COUNT/$lane.best_score.gz  $mytmp/COUNT/$lane.best_score.old.gz
    gunzip -c $mytmp/COUNT/$lane.best_score.old.gz  $mytmp/PHITS_$target/$lane.hits.gz |  gawk -F '\t' '/^#/{next;}{p=$1;s=$2;if (ss[p]<s)ss[p]=s;}END{for(p in ss)printf("%s\t%s\n",p,ss[p]);}' | gzip >  $mytmp/COUNT/$lane.best_score.gz 
    \rm  $mytmp/COUNT/$lane.best_score.old.gz 
  endif

  if (! -e tmp/COUNT/$lane.noInsert.gz && -e $mytmp/PHITS_$target/$lane.noInsert.gz) then
    cp $mytmp/PHITS_$target/$lane.noInsert.gz tmp/COUNT/$lane.noInsert.gz
  endif

  touch $mytmp/PHITS_$target/$lane.done
end # target

# -doubleSeed

#######################################################################################
## Select and export the best alignments
## speyrubo  is the archetypal mito pseudogene present in the genome

phasea3:
echo "#############################"
echo "#############################"
echo -n "## a3: select best $run $lane "
date

if ($COUNT == COUNT) then
  set toto=""
  if (-e tmp/COUNT/$lane.hits.gz) then
    foreach ii (10 9 8 7 6 5 4 3 2 1)
      if (-e tmp/COUNT/$lane.hits.$ii.gz) then
        set jj=$ii
        @ jj = $ii + 1
        mv tmp/COUNT/$lane.hits.$ii.gz tmp/COUNT/$lane.hits.$jj.gz
      endif
    end
    mv tmp/COUNT/$lane.hits.gz tmp/COUNT/$lane.hits.1.gz
    set hitfiles="tmp/COUNT/$lane.hits.1.gz $hitfiles "
  else
    if (-e tmp/COUNT/$lane.hits.1.gz) then
      set hitfiles="tmp/COUNT/$lane.hits.1.gz $hitfiles "
    endif
  endif
endif

  set geneRemap=" -geneRemap tmp/METADATA/mrnaRemap.gz "

  set filter=$species
  if ($minAli < 24) set filter=$minAli
  set sig=""
  if (-e TARGET/Targets/$species.av.signature_transcripts_list.txt) set sig="-sigTargets  TARGET/Targets/$species.av.signature_transcripts_list.txt"
  echo "hitfiles= # $hitfiles #"
  echo "gunzip -c $hitfiles | bin/bestali  -filter $filter  $mIL $geneRemap -maxHit $maxHit -countBest -seqc -strategy $Strategy -exportBest  -exportSuffix -exportVenn -errorProfile $pair -aliProfile -exportMito $sig -gzo -o $mytmp/COUNT/$lane"
  (bin/time -p gunzip -c $hitfiles | sort -k 1,1 -k 2,2nr | bin/bestali  -filter  $filter   $mIL $geneRemap -maxHit $maxHit -countBest -seqc -strategy $Strategy -exportBest  -exportVenn -exportSuffix -errorProfile $pair -aliProfile -exportMito $sig -gzo -o $mytmp/COUNT/$lane) >&  $mytmp/COUNT/$lane.err
  touch $mytmp/COUNT/$lane.mrnaOrder.done
  ls -ls $mytmp/COUNT/$lane.hits.gz
  echo -n "bestali done : "
  date
  if ($cleanUp == 1 && -e tmp/COUNT/$lane.hits.1.gz) \rm  tmp/COUNT/$lane.hits.1.gz
  
    set ff=$mytmp/COUNT/$lane.hits.gz
    if (-e $ff) then
      mv $ff $ff.old
      gunzip -c $ff.old | scripts/tab_sort -k 1,1 -k 2,2nr -k 8,8 -k 11,11 -k 26,26  -k 6,6n -k 7,7n  | gzip > $ff 
      \rm $ff.old
    endif
    touch  tmp/COUNT/$run/snp1_sorting.done

# recover the original identifiers
# since it is really costly i blank out this option
set needAlias=0
if (0) then
  foreach run2 (`cat $mytmp/MetaDB/$MAGIC/RunPairedList`)
    if ($RecoverId == 1 && $run2 == $run && -e Fastc/$lane.id.gz) set needAlias=1
  end
  if ($needAlias == 1) then
    if (! -d $mytmp/ALIa) mkdir $mytmp/ALIa
    if (! -d $mytmp/ALIa/$run) mkdir $mytmp/ALIa/$run

      echo -n "a123: ventilate to recover the original identifiers $lane "
      date
      (bin/time -p gunzip -c $mytmp/COUNT/$lane.hits.gz | bin/bestali -ventilate 64 -id_file Fastc/$lane.id.gz -gzo -o $mytmp/ALIa/$lane) >&  $mytmp/ALIa/$lane.err
      echo -n "a3: recover the original identifiers done "
      date
 
  endif
endif

if (-e $mytmp/COUNT/$lane.hits.gz ) then
  echo ZZZZZ | gzip > $mytmp/ZZZZZ.gz
  touch $mytmp/COUNT/$lane.too_many_hits.gz
  gunzip -c $mytmp/COUNT/$lane.hits.gz $mytmp/ZZZZZ.gz $mytmp/COUNT/$lane.too_many_hits.gz |  gawk -F '\t' '/^ZZZZZ/{zz=1;next;}{if(zz<1)ok[$1]=1;}{if(ok[$1]==1)next;split($1,aa,"#");split(aa[2],bb,"/");s=1;if(bb[1]>0)t=bb[1];else t=1;}/BadScore/{bs++;bt+=t;next;}/Partial/{bs++;bt+=t;next;}/OK_Multi/{next}/Multi/{s10++;t10+=t;}END{if(s10>0)printf("At_least_10_sites %d seq %d tags \n", s10, t10) ;if(bt>0)printf("Low_quality_mapping %d seq %d tags \n", bs, bt) ;}' > $mytmp/COUNT/$lane.badcounts.txt
endif

date

# extract the fasta file of unaligned sequences

if (! -e tmp/Unaligned/$lane.$dna2dnaFormat.gz) then
  echo "extract the fasta file of unaligned sequences"
  echo "//" > $mytmp/COUNT/$lane.ok
  if (-e $mytmp/COUNT/$lane.hits.gz ) then
    gunzip -c $mytmp/COUNT/$lane.hits.gz | gawk '{gsub (/>$/,"",$1);gsub (/<$/,"",$1);if($1 != old) print $1;old=$1;}' >> $mytmp/COUNT/$lane.ok
  endif
  if (-e tmp/COUNT/$lane.noInsert.gz) then
    gunzip -c tmp/COUNT/$lane.noInsert.gz | cut -f 1  | gawk '{gsub (/>$/,"",$1);gsub (/<$/,"",$1);if($1 != old) print $1;old=$1;}' >>  $mytmp/COUNT/$lane.ok
  endif
  echo "bin/dna2dna -I $dna2dnaFormat -i Fastc/$lane.$dna2dnaFormat.gz -gzi -gzo  -minEntropy $minEntropy -reject   $mytmp/COUNT/$lane.ok -O $dna2dnaFormat -count -o tmp/Unaligned/$lane "
        bin/dna2dna -I $dna2dnaFormat -i Fastc/$lane.$dna2dnaFormat.gz -gzi -gzo  -minEntropy $minEntropy -reject   $mytmp/COUNT/$lane.ok -O $dna2dnaFormat -count -o tmp/Unaligned/$lane 
  \rm  $mytmp/COUNT/$lane.ok
endif

date

#######################################################################################
## cleanup
echo cleanup

  if (! -d tmp/$COUNT) mkdir tmp/$COUNT
  if (! -d tmp/$COUNT/$run) source scripts/mkDir $COUNT $run
  if ($needAlias == 1 && ! -d tmp/ALIa)   mkdir tmp/ALIa
  if ($needAlias == 1 && ! -d tmp/ALIa/$run) scripts/mkDir ALIa $run

  touch tmp/$COUNT/$lane.toto

  # the genome lane hits are needed for double_collect
  if (-e $mytmp/COUNT/$lane.err) mv $mytmp/COUNT/$run/* tmp/$COUNT/$run
  if (-e $mytmp/ALIa/$lane.err) mv $mytmp/ALIa/$run/* tmp/ALIa/$run


foreach target ($targets)

  if (! -d tmp/PHITS_$target) mkdir tmp/PHITS_$target
  if (! -d tmp/PHITS_$target/$run) source scripts/mkDir PHITS_$target $run


  if (-e $mytmp/PHITS_$target/$lane.sam.gz) then
    mv $mytmp/PHITS_$target/$lane.sam.gz tmp/PHITS_$target/$lane.sam.gz
  endif
  if (-e $mytmp/PHITS_$target/$lane.hits.gz) then
    if (0 && $Strategy == RNA_seq && $target == genome ) then 
       echo -n "a3: ventilate the genome hits in BG format $run $lane "
       date
       (bin/time -p bin/wiggle -ventilate -I BHIT -O BHIT -i $mytmp/PHITS_genome/$lane.hits.gz -o tmp/PHITS_genome/$lane -gzo) >& tmp/PHITS_genome/$lane.BG.err
       echo  "a3: ventilate done "
    endif
    # mv $mytmp/PHITS_$target/$lane.hits.gz tmp/PHITS_$target/$lane.hits.gz
    if ($cleanUp == 1) \rm $mytmp/PHITS_$target/$lane.hits.gz
  endif
  if ($target == gdecoy) \rm  $mytmp/PHITS_$target/$run/*overhang*
  if (-e $mytmp/PHITS_$target/$lane.err) mv $mytmp/PHITS_$target/$run/* tmp/PHITS_$target/$run

end

echo  "a3: jobstats "
scripts/jobstats.tcsh $run $lane 1 0
touch tmp/$COUNT/$lane.align.done

echo  "a3: cleanup "

if ($cleanUp == 1 && -d $mytmp) \rm -rf $mytmp
echo -n "done "
date

  if (-e tmp/$COUNT/$MAGIC.job_statistics.txt) \rm tmp/$COUNT/$MAGIC.job_statistics.txt
  if (-d tmp/GENERUNS/$run) \rm -rf  tmp/GENERUNS/$run
  if (-d tmp/GENELANES/$run) \rm tmp/GENELANES/$lane.*
  if (-e tmp/$COUNT/$run/c2.alistats.ace) \rm tmp/$COUNT/$run/c2.alistats.ace

ls -ls TARGET/Targets/$species.virus.fasta.gz TARGET/Targets/$species.bacteria.fasta.gz > tmp/$COUNT/$run/virus_bacteria.date

exit 0


##############################################################

foreach run (`cat MetaDB/$MAGIC/RunList`)
  set n=`cat Fastc/$run/LaneList | wc -l`
  if ($n < 150) continue
  foreach ii (`seq 0 1 9`)
    mkdir Fastc/$run.$ii
    mv Fastc/$run/f2.*$ii.* Fastc/$run.$ii 
  end
end

foreach run (`cat MetaDB/$MAGIC/RunList`)
  set n=`cat Fastc/$run/LaneList | wc -l`
  if ($n < 150) continue
  foreach ii (`seq 0 1 9`)
    \cp Fastc/$run/*.ace Fastc/$run/Max_probe_length Fastc/$run/first_line.txt Fastc/$run/original.count Fastc/$run.$ii 
    pushd Fastc/$run.$ii 
      ls *.fastc.gz | gawk '{gsub(/.fastc.gz/,"",$1); printf("%s.%d/%s\n",run,ii,$1);}' run=$run ii=$ii > LaneList
    popd
  end
end

echo ' ' > subsublib.ace
foreach run (`cat MetaDB/$MAGIC/RunList`)
  set n=`cat Fastc/$run/LaneList | wc -l`
  if ($n < 150) continue
  foreach ii (`seq 0 1 9`)
    echo "$run $ii" | gawk '{printf("Run %s.%d\nRNA\nIllumina\nSublibrary_of %s\nPaired_end\nProject %s\n\n", run,ii,run,magic);}' run=$run ii=$ii magic=$MAGIC >>  subsublib.ace
  end
end

foreach run (`cat MetaDB/$MAGIC/RunList`)
  set n=`cat Fastc/$run/LaneList | wc -l`
  if ($n < 150) continue
  \rm Fastc/$run/LaneList
end

foreach run (`cat tutu`)
  foreach ii (0 1 2 3 4 5 6 7 8 9)
  cat Fastc/$run/original_counts.ace | gawk '/^Raw_data/{printf ("%s %d %s %d %s %d %s\n",$1,$2/10,$3,$4/10,$5,$6/10,$7);next;}' > Fastc/$run.$ii/original_counts.ace
  end
end

foreach ff (`ls Fastc/SRR1000*.?/original_counts.ace`)
  set run=`echo $ff | gawk '{split($1,aa,"/");print aa[2];}'`
  cat $ff | gawk '/^Raw_data/{printf ("Ali %s\n",run);print;printf("\n");}' run=$run  >_f
  mv _f $ff
end

# for some reason i doubled the LaneList files
foreach run (`cat tutu`)
  foreach ii (0 1 2 3 4 5 6 7 8 9)
    cat Fastc/$run.$ii/LaneList | sort -u | sort -V > _f
    mv _f Fastc/$run.$ii/LaneList
  end
end
