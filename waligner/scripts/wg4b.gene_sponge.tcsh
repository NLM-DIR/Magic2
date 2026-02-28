#! bin/tcsh -f

## scan the wiggle files and count the genebox support 

set phase=$1

if ($phase == Cumul) goto phaseCumul
if ($phase == Parse) goto phaseParse

set WG=$2
set run=$3
set target=$4
if (! -d tmp/$WG/$run) exit 0
    
echo "wg4b.gene_sponge.tcsh $phase $WG $run $target"
if ($phase == noCapture) goto phaseNoCapture 
goto phaseCapture

#####################################################

phaseNoCapture:

set tt=tmp/$WG/$run/$target.$phase.tsf
echo "# Gene\tRun\tf\tchrom\ta1\ta2\tkb" > $tt
foreach chrom ($chromSetAll)
  foreach fr (f r)
    set w=tmp/$WG/$run/$chrom/R.chrom.u.$fr.BF.gz
    if (! -e $w) continue
    set spongeF=tmp/METADATA/gtf.$target.$fr.sponge.gz
    if (! -e $spongeF) continue
    geneelements -wiggle $w -sponge 1 -spongeFile $spongeF  -sxxChromosome $chrom -run $run | gawk -F '\t' '/^level_1/{r=$2;g=$3;c=$5;a1=$6;a2=$7;z=$11/1000.0;printf("%s\t%s\ttiif\t%s\t%d\t%d\t%.1f\n",g,r,c,a1,a2,z);}' >> $tt
  end
end

exit 0

#####################################################

phaseCapture:


set tt=tmp/$WG/$run/$target.$phase.tsf
echo "# Gene\tRun\tf\tchrom\ta1\ta2\tkb" > $tt
foreach chrom ($chromSetAll)
  foreach fr (f r)
    set w=tmp/$WG/$run/$chrom/R.chrom.u.$fr.BF.gz
    if (0 && -e tmp/$WG/$run/$chrom/R.chrom.frns.pp.BF.gz)  w="tmp/$WG/$run/$chrom/R.chrom.frns.pp.BF.gz"
    if (! -e $w) continue
    set spongeF=tmp/METADATA/gtf.$target.$fr.sponge.gz
    set spongeF=tmp/METADATA/$target.$fr.any_genebox.$capture.sponge
    if (! -e $spongeF) continue
    geneelements -wiggle $w -sponge 1 -spongeFile $spongeF  -sxxChromosome $chrom -run $run | gawk -F '\t' '/^level_1/{r=$2;g=$3;c=$5;a1=$6;a2=$7;z=$11/1000.0;printf("%s\t%s\ttiif\t%s\t%d\t%d\t%.1f\n",g,r,c,a1,a2,z);}' >> $tt
  end
end

exit 0

#####################################################
## CAPTURE draft scripts, to be moved in their own wscript when ready
# capture wiggle stats extracted from SPONGE
           tmp/WIGGLERUN/Nanopore_ROCR3_B-F1+2/wg4b.av.ns.any_genebox.R1.sponge.count
set tt=tmp/$WG/$run/$target.$phase.tsf
phaseCumul:

set target=`echo $Etargets | gawk '{print $1; last}'`
set type=genebox

foreach capture ($CAPTURES)
  set toto=tmp/WIGGLERUN/$MAGIC.$capture.capture.count
  echo "//" >  $toto.ace
  foreach fr (f r ns)
    echo '#' > $toto
    foreach run (`cat MetaDB/$MAGIC/RunsList`)
      set WG=toto
      if (-d tmp/WIGGLERUN/$run) set WG=WIGGLERUN
      if (-d tmp/WIGGLEGROUP/$run) set WG=WIGGLEGROUP
      if (-e  tmp/$WG/$run/wg4b.$target.$fr.any_$type.$capture.sponge.count ) then
        cat tmp/$WG/$run/wg4b.$target.$fr.any_$type.$capture.sponge.count | gawk '/^level_10/{next;}/^level_1/{n++;nn+=$11/1000.0;}END{printf("%s\t%g\t%d\n",run,nn,n);}' run=$run >> $toto
      endif
    end

    set groupLevelMax=`cat  MetaDB/$MAGIC/g2r | cut -f 3 | sort -k 1n | tail -1`
    foreach level (`seq 1 1 $groupLevelMax`)
      foreach group (`cat MetaDB/$MAGIC/g2r | gawk '{if($3 == level) print $1;}' level=$level | sort -u `)
        cat MetaDB/$MAGIC/g2r ZZZZZ $toto | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){if($1 == g)rr[$2]=1;next;}}{if(rr[$1]==1){bp+=$2;nn+=$3;}}END{printf("%s\t%g\t%d\n", g,bp,nn);}' g=$group > $toto.1
        cat $toto >> $toto.1
        cat $toto.1 | sort -u > $toto
        \rm $toto.1
      end
    end
    cat $toto | gawk -F'\t' '{if ($3>0)printf("Ali %s\nCapture %s.%s.%s %s kb %d blocks\n\n",$1,c,type,fr,$2,$3);}' c=$capture type=$type fr=$fr >> $toto.ace
    \rm $toto
  end
  wc $toto.ace

end





####
## capture ali stats extracted from GENEINDEX


foreach capture ($CAPTURES)
  set toto=tmp/WIGGLERUN/$MAGIC.$capture.capture.ali_count
  echo "//" >  $toto.ace

  cat  TARGET/GENES/$capture.capture.av.gene_list ZZZZZ  tmp/GENEINDEX/$MAGIC.av.GENE.u.ace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$1]=1;next;}}/^Transcript/{ok=0;next;}/^Gene/{ok=0;gsub(/\"/,"",$0);g=$2;if(gg[g]==1)ok=1;next;}/^Run_U/{gsub(/\"/,"",$2);r=$2;if(ok==1){nr[r]+=$6; mb[r]+=$8;ng[r]++;}}END{for (r in nr)printf("%s\t%g\tkb\t%d\treads\t%d\tgenes\n", r, mb[r],nr[r],ng[r]);}' > $toto

  set groupLevelMax=`cat  MetaDB/$MAGIC/g2r | cut -f 3 | sort -k 1n | tail -1`
  foreach level (`seq 1 1 $groupLevelMax`)
    set okk=0
    foreach group (`cat MetaDB/$MAGIC/GroupList`)
      cat MetaDB/$MAGIC/g2r ZZZZZ $toto | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){if($1 == g)rr[$2]=1;next;}}{if(rr[$1]==1){bp+=$2;nr+=$4;}}END{printf("%s\t%.2f\tkb\t%d\treads\n", g,bp,nr);}' g=$group > $toto.1
      cat $toto >> $toto.1
      cat $toto.1 | sort -u > $toto
      \rm $toto.1
    end
  end
  cat $toto | gawk -F'\t' '{printf("Ali %s\nCapture %s_ali %s kb %d reads\n\n",$1,c,$2,$4);}' c=$capture  >> $toto.ace
  \rm $toto
  wc $toto.ace
end

exit 0

#####################################################

phaseParse:

cat tmp/WIGGLERUN/$MAGIC.*capture*.ace > tmp/WIGGLERUN/$MAGIC.Capture.ace
cat tmp/WIGGLERUN/$MAGIC.Capture.ace | gawk '/^Ali/{a=$2;next;}/^Capture/{printf("%s\t",a);print;}' | sort | gawk -F '\t' '{if(a != $1)printf("\nAli %s\n",$1);a=$1;printf("%s\n",$2);}' > tmp/WIGGLERUN/$MAGIC.Capture.sorted.ace
\rm tmp/WIGGLERUN/$MAGIC.Capture.ace

echo "pparse  tmp/WIGGLERUN/$MAGIC.Capture.sorted.ace" | bin/tacembly MetaDB -noprompt

exit 0

#####################################################
