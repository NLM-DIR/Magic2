#!bin/tcsh -f

set run=$1
set lane=$2
set s=$3
set paired=$4

# without lane the code would do all the lanes of the run
# but the output name should be changed   either as $lane/K    or as $run/K

set out_step="-out_step 10"
if ($?wiggle_step) then
   set out_step="-out_step $wiggle_step"
endif

set param="-O BV $out_step"
set stranded=""
if ($NO_INTRON == 2) then
  set stranded=""
  if ($s == strand) set stranded="-strand"
  if ($s == antistrand) set stranded="-antistrand"
endif

echo "$run strand=$stranded $param"

set ici=`pwd`

set mytmp=$TMPDIR/aceview.wigglerun.$$
# set mytmp=$TMPDIR/aceview.wigglerun.11702
# set mytmp=tmp/WIGTMP

echo  $mytmp
mkdir $mytmp
mkdir $mytmp/$run
mkdir $mytmp/$lane
if (! -d tmp/WIGGLELANE/$run) mkdir  tmp/WIGGLELANE/$run
if (! -d tmp/WIGGLELANE/$lane) mkdir  tmp/WIGGLELANE/$lane

cp bin/wiggle $mytmp
cp scripts/target2target_class.txt $mytmp

set remapon=""
echo $Strategy
echo RNAtargets="$RNAtargets"
  if ($Strategy == RNA_seq && $NO_INTRON != 2 && "$RNAtargets" != "") then
echo hello1a
    if (-e tmp/METADATA/mrnaRemap.gz) then
      set remapon="  -mapon  tmp/METADATA/mrnaRemap.gz "
    else
      echo "missing files tmp/METADATA/mrnaRemap.gz"
      echo "phase wg1 strategy=$Strategy cannot proceed"
      exit 1
    endif
  endif

# goto laba 2024
# goto laba 

echo -n "splitting the hits file per chromosomes "
date
echo $mytmp/$run 

    if (! -e  tmp/COUNT/$lane.hits.gz) exit 1
    if (! -d  $mytmp/$lane) mkdir $mytmp/$lane
 ls -ls tmp/COUNT/$lane.hits.gz

 set pair=""
 if ($paired == pairs)  set pair="-pair 500"

set isIlm=`cat  MetaDB/$MAGIC/runs.ace | gawk '/^Run /{gsub(/\"/,"",$2);ok=0;if($2==run)ok=1;}/^Illumina/{if(ok==1)okk=1;}END{print 0+okk;}' run=$run`
set isNanopore=`cat  MetaDB/$MAGIC/runs.ace | gawk '/^Run /{gsub(/\"/,"",$2);ok=0;if($2==run)ok=1;}/anopore/{if(ok==1)okk=1;}/PacBio/{if(ok==1)okk=1;}END{print 0+okk;}' run=$run`

# BAD: no partial blocks the export of the read-ends  2024-03-12
# set noPartial=`cat  MetaDB/$MAGIC/runs.ace | gawk '/^Run /{gsub(/\"/,"",$2);ok=0;if($2==run)ok=1;}/Entry_adaptor/{if(ok==1)okk=1;}END{print 0+okk;}' run=$run`
set noPartial=""

if ($isIlm == 1) then
  set maxWigErr=3
  set maxWigErrRate=3
else
  set maxWigErr=0
  set maxWigErrRate=0
endif



 # u nu pp   2024
 foreach uu (u nu pp)
   if (-e tmp/WIGGLERUN/$run/wg2a.$uu.done) continue
   set filter="" 
   if ($uu == u) set filter="-unique"
   if ($uu == nu) set filter="-non_unique"
   if ($uu == pp) set filter="-partial"
   if ($isNanopore == 1 && $uu == pp) continue
   if ($isNanopore == 1 || $noPartial == 1) set filter="$filter -noPartial"
  ## split recursivelly the hits so we scan once the big file and split it per chromosome or per target
  ## ATTENTION do nor -gzo: there would be to many gzip pipes generated on the computer
    # if we set -minErrRate 7, we obtain in addition a second set of K7 wiggle using only error rich reads
    echo "$mytmp/wiggle -strategy $Strategy -ventilate $filter $stranded $remapon -i tmp/COUNT/$lane.hits.gz -o $mytmp/$lane -minErrRate 0  -minAliRate 70 -minAliLength 70 -maxErr $maxWigErr -maxErrRate $maxWigErrRate $pair -I BHIT -O BG"
          $mytmp/wiggle -strategy $Strategy -ventilate $filter $stranded $remapon -i tmp/COUNT/$lane.hits.gz -o $mytmp/$lane -minErrRate 0  -minAliRate 70 -minAliLength 70 -maxErr $maxWigErr -maxErrRate $maxWigErrRate $pair -I BHIT -O BG
  end

echo -n "splitting done tralala"
date

#exit 0
laba:
# RefSeq takes at lest 14Gb of memory
# $chromSetAll
# setenv chromSetAll CHROMOSOME_III
foreach chrom ($chromSetAll)
  cd $ici
  set tag="null"
  setenv target $chrom
  foreach chrom2 ($chromSetAll)
    if ($chrom == $chrom2) setenv target genome
  end
  if ($chrom == genome) continue
  if ($chrom == gdecoy) continue

  source $mytmp/target2target_class.txt
echo ## --- $chrom $target
  set myf1=$mytmp/$run
  set myf2=$target_class
  set doremapon=0
  if ($target == genome) then
    set myf2=$chrom
    set remapon=""
    if ($Strategy == RNA_seq && -e  tmp/METADATA/mrnaRemap.gz) then
      set doremapon=1
      gunzip -c tmp/METADATA/mrnaRemap.gz | gawk -F '\t' '{if ($5 == chrom && $6<$7)print;}' chrom=$chrom > $mytmp/$lane.remapF
      gunzip -c tmp/METADATA/mrnaRemap.gz | gawk -F '\t' '{if ($5 == chrom && $6>$7)print;}' chrom=$chrom > $mytmp/$lane.remapR
      set remaponF="  -mapon $mytmp/$lane.remapF"
      set remaponR="  -mapon $mytmp/$lane.remapR"
    endif
  endif

  #################### 
set frs="f r"
if ($Strategy == RNA_seq) set frs="f r ELF ELR ERF ERR"

  foreach map (mapped remapped)
    # if ($map == remapped && $target != genome) continue 
    if ($map == remapped && $Strategy != RNA_seq) continue 
    echo "### $target $map $chrom"
    # (u nu pp 2024)
    foreach uu (u nu pp)
      if (-e tmp/WIGGLERUN/$run/wg2a.$uu.done) continue
      set frs="f r"
      if ($Strategy == RNA_seq && $uu == u) set frs="f r ELF ELR ERF ERR"
      set out="$mytmp/$lane/K.$map.$chrom.$uu"
          # in mapped and in remapped case, merge for each category the individual BG hits into a wiggle
       foreach fr ($frs)
            if (-e  $mytmp/$lane.$map.$myf2.$fr.minerr0.$uu.BG) then
              echo "cat $mytmp/$lane.$map.$myf2.$fr.minerr0.$uu.BG | sort | $mytmp/wiggle -I BG $param -strand -o $out.$fr "
              date
                    cat $mytmp/$lane.$map.$myf2.$fr.minerr0.$uu.BG | sort | $mytmp/wiggle -I BG $param -strand -o $out.$fr 
              date
            else
              echo "missing file $mytmp/$lane.$map.$myf2.$fr.minerr0.$uu.BG"
            endif
       end
       # in mapped case, we are done
       # in remap case, we must remap the genes to the chromosomes
       # if a gene is covered 1000 times, it is 1000 times faster to remap the wiggle than to remap the individual BG hits
       if ($map == remapped) then

          if ($doremapon == 1) then 
             # while remapping data may change strand, for each class we separate the 2 cases
             foreach fr ($frs)
                if (-e   $out.$fr.BV) then

		  echo "###### cat  $out.$fr.BV | $mytmp/wiggle -I BV $param $remaponF -o $mytmp/$lane/F.$fr.$map.$chrom.$uu"
                               cat  $out.$fr.BV | $mytmp/wiggle -I BV $param $remaponF -o $mytmp/$lane/F.$fr.$map.$chrom.$uu
		  echo "###### cat  $out.$fr.BV | $mytmp/wiggle -I BV $param $remaponR -o $mytmp/$lane/R.$fr.$map.$chrom.$uu"
                               cat  $out.$fr.BV | $mytmp/wiggle -I BV $param $remaponR -o $mytmp/$lane/R.$fr.$map.$chrom.$uu
		endif
             end
             # then we rename the flipped files and merge them
             foreach fff ( F.f_R.r_f F.r_R.f_r  F.EL_R.ER_EL F.ER_R.EL_ER F.ELF_R.ERR_ELF  F.ERF_R.ELR_ERF   F.ELR_R.ERF_ELR   F.ERR_R.ELF_ERR )
               set x1=`echo $fff | gawk -F _ '{print $1}'`
               set x2=`echo $fff | gawk -F _ '{print $2}'`
               set x3=`echo $fff | gawk -F _ '{print $3}'`
               if (-e $mytmp/$lane/$x1.$map.$chrom.$uu.BV || -e $mytmp/$lane/$x2.$map.$chrom.$uu.BV || -e $mytmp/$lane/K.$map.$chrom.$uu.$x3.BV) then
                  echo "cat $mytmp/$lane/$x1.$map.$chrom.$uu.BV  $mytmp/$lane/$x2.$map.$chrom.$uu.BV  |  $mytmp/wiggle -I BV -O BV $out_step -gzo -o $mytmp/$lane/KK.$map.$chrom.$uu.$x3"
                        cat $mytmp/$lane/$x1.$map.$chrom.$uu.BV  $mytmp/$lane/$x2.$map.$chrom.$uu.BV  |  $mytmp/wiggle -I BV -O BV $out_step -gzo -o $mytmp/$lane/KK.$map.$chrom.$uu.$x3
                        mv  $mytmp/$lane/KK.$map.$chrom.$uu.$x3.BV.gz  $mytmp/$lane/K.$map.$chrom.$uu.$x3.BV.gz 
                endif
             end
          endif

       endif

    end  # end uu
  end   # end map
  ####################
gzip  $mytmp/$lane/K.mapped.*.BV
mv $mytmp/$lane/K.*.gz tmp/WIGGLELANE/$lane
# optionally
# mv $mytmp/$lane/* tmp/WIGGLELANE/$lane

end  # chrom/target

#############################################################################
## cleanup

cleanup:
\rm -rf $mytmp

touch tmp/WIGGLELANE/$lane/wg1.done
