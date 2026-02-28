#!bin/tcsh -f

set run=$1
set uu=$2
set phase=$3

set out_step="-out_step 10"
if ($?wiggle_step) then
   set out_step="-out_step $wiggle_step"
endif

if ($phase == lane2run) goto lane2run
if ($phase == restrand) goto restrand

exit 0

lane2run:
foreach chrom (mito SpikeIn $chromSetAll)
  echo "... wg2a $run $chrom"
  if (! -d tmp/WIGGLERUN/$run/$chrom) mkdir tmp/WIGGLERUN/$run/$chrom

      foreach fr (f r ELF ELR ERF ERR)

          # 2013_03_25  the genes cases is not programmed in WIGGLELANE, the -ventilate option does not support it
          #   tmp/WIGGLERUN/$run/$chrom/R.genes.$uu.$fr.BF.gz

          if (! -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr.BF.gz) then
            set n=`ls -ls tmp/WIGGLELANE/$run/*/K.*mapped.$chrom.$uu.$fr.BV.gz | gawk '{if($6==20)next;n++;}END{print 0+n;}'`
            if ($n > 0) then
              echo "Constructing tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr.BF.gz"
              echo     "tmp/WIGGLELANE/$run/*/K.*mapped.$chrom.$uu.$fr.BV.gz | bin/wiggle -I BV -O BF $out_step -gzo -o tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr"
              gunzip -c tmp/WIGGLELANE/$run/*/K.*mapped.$chrom.$uu.$fr.BV.gz | bin/wiggle -I BV -O BF $out_step -gzo -o tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr
	    endif
          endif
      end
      if ((-e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.f.BF.gz || -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.r.BF.gz) && ! -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.frns.$uu.BF.gz) then
        gunzip -c tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.[fr].BF.gz | bin/wiggle -I BF -O BF $out_step -gzo -o  tmp/WIGGLERUN/$run/$chrom/R.chrom.frns.$uu -cumul
     endif
      if (( -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.ELF.BF.gz || -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.ELR.BF.gz) && ! -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.EL.$uu.BF.gz) then
        gunzip -c tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.EL[FR].BF.gz | bin/wiggle -I BF -O BF $out_step -gzo -o  tmp/WIGGLERUN/$run/$chrom/R.chrom.EL.$uu -cumul
     endif
      if (( -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.ERF.BF.gz || -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.ERR.BF.gz) && ! -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.ER.$uu.BF.gz) then
        gunzip -c tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.ER[FR].BF.gz | bin/wiggle -I BF -O BF $out_step -gzo -o  tmp/WIGGLERUN/$run/$chrom/R.chrom.ER.$uu -cumul
     endif

  # 2015_06_12  Nettoyage des fuites de strand, est ce necessaire ?


end


touch tmp/WIGGLERUN/$run/wg2a.$uu.done
#\rm -rf tmp/WIGGLELANE/$run/*/*.$uu.*.BV.gz
goto done

restrand:
set WG=WIGGLERUN
# reevaluate the strand of the run
    set tt=tmp/$WG/$run/projected.stranding.tsf
    echo '#' > $tt
    echo '#' > $tt.1
      set ok=0
      foreach chrom ($chromSetAll)
        set sBV=TARGET/CHROMS/$chrom.stranding.BV
	if (! -e $sBV) continue
        set wf=tmp/$WG/$run/$chrom/R.chrom.u.f.BF.gz
        set wr=tmp/$WG/$run/$chrom/R.chrom.u.r.BF.gz
	if (! -e $wf || ! -e $wr) continue
        bin/wiggle -i $wf  -I BF -O BV | gawk '{n=$1+0;if(n>0)printf("%d\tf\tii\t%d\t0\n",$1,$2);}' > $tt.a
        bin/wiggle -i $wr  -I BF -O BV | gawk '{n=$1+0;if(n>0)printf("%d\tr\tii\t0\t%d\n",$1,$2);}' > $tt.b
        cat $tt.a $tt.b | bin/tsf -I tsf -O tsf -s ns >$tt.c

        if (0) then   # do this once using run==the stranded group
          cat $tt.c | gawk -F '\t' '/^#/{next;}{x=$1;p=($4+20)/($4+$5+20);q=int(100*p);printf("%d\t%d\n",x,q);}' > $tt.d
	  cat $tt.d | gawk '{p=$2;if(p>90 || p < 10) print;}' > $sBV
        endif
	
        cat $sBV $tt.c | gawk -F '\t' '/^#/{next}{if(NF==2){pp[$1]=2*($2-50);next;}}{p=pp[$1];if(p>97){sp+=$4;sm+=$5;};if(p<-97){sp+=$5;sm+=$4;}}END{printf("%s\t%s\tiif\t%d\t%d\t%.2f\n",run,chrom,sp,sm,100*sp/(sp+sm+.001));}' run=$run chrom=$chrom >> $tt.1
	set ok=1
	\rm $tt.a $tt.b $tt.c
      end
    if ($ok == 1) then
      if (-e tmp/WIGGLELANE/$run/wg1.antistranded) then
        bin/tsf -I tsf -O tsf  -i $tt.1 -s any | gawk -F '\t' '/^#/{next;}{sp=$5;sm=$4;printf("%s\tP_genome\tiif\t%d\t%d\t%.3f\n",run,sp,sm,100*sp/(sp+sm+.001));}' run=$run >> $tt
      else
        bin/tsf -I tsf -O tsf  -i $tt.1 -s any | gawk -F '\t' '/^#/{next;}{sp=$4;sm=$5;printf("%s\tP_genome\tiif\t%d\t%d\t%.3f\n",run,sp,sm,100*sp/(sp+sm+.001));}' run=$run >> $tt
      endif
      \rm $tt.1
    endif
 
# if the stranding is below 40, we flip the stand of all the files
if ($Strategy == RNA_seq && ! -e tmp/$WG/$run/wg2a.restrand.done) then
    set ss=`cat $tt | gawk '/^#/{s=50}{s=$6}END{print int(s);}'`
    echo "restrand ss=$ss"
    if ($ss < 40 && ! -e tmp/WIGGLELANE/$run/wg1.antistranded) then
       echo $ss >  tmp/$WG/$run/wg2a.restrand.yes
      echo "restanding needed ss=$ss"

      foreach chrom ($chromSetAll)
         foreach uu (u nu pp)
           set wf=tmp/$WG/$run/$chrom/R.chrom.$uu.f.BF.gz
           set wr=tmp/$WG/$run/$chrom/R.chrom.$uu.r.BF.gz
	   echo "exchanging $wf $wr"
	   cp $wf $wf.ok
	   cp $wr $wr.ok
	   mv $wf $wf.z
	   mv $wr $wf
	   mv $wf.z $wr
        end
        set wf=tmp/$WG/$run/$chrom/R.chrom.u.ELF.BF.gz
        set wr=tmp/$WG/$run/$chrom/R.chrom.u.ELR.BF.gz
        mv $wf $wf.z
        mv $wr $wf
        mv $wf.z $wr
        set wf=tmp/$WG/$run/$chrom/R.chrom.u.ERF.BF.gz
        set wr=tmp/$WG/$run/$chrom/R.chrom.u.ERR.BF.gz
        mv $wf $wf.z
        mv $wr $wf
        mv $wf.z $wr
      end
   else
       echo $ss >  tmp/$WG/$run/wg2a.restrand.no
    endif
endif
touch tmp/$WG/$run/wg2a.restrand.done

  foreach chrom ($chromSetAll)
    echo "Construct the transcriptsEnds  $WG/$run"
    echo "  bin/wiggle  -transcriptsEnds tmp/$WG/$run/$chrom/R.chrom.u -gzi -I BF -O COUNT -o tmp/$WG/$run/$chrom/wg2a -minCover 30 -wiggleRatioDamper 5"
            bin/wiggle  -transcriptsEnds tmp/$WG/$run/$chrom/R.chrom.u -gzi -I BF -O COUNT -o tmp/$WG/$run/$chrom/wg2a -minCover 30 -wiggleRatioDamper 5

    foreach cover (5 10 20 50)
        gunzip -c  tmp/$WG/$run/$chrom/R.chrom.frns.$uu.BF.gz | bin/wiggle -I BF -O BV  -gauss 20 -minCover $cover -peaks -o tmp/$WG/$run/$chrom/coverome.$cover.$uu
    end

    if ($uu == u) then
      cat  tmp/$WG/$run/$chrom/coverome.$minCoveron.u.peaks | gawk '/^#/{next;}{nn++;chrom=$1;a1=$2;a2=$3;if(a1>a2){a0=a1;a1=a2;a2=a0;}printf("%s\t%09d\t%09d\tCOVERON\tCoveron.%s.%s_%d\t%d\t%s\t%s\t%s\n",chrom,a1,a2,run,chrom,nn,$4,$5,$6,$7); }' run=$run >> toto.sp1
      cat toto.sp1 | sort -k 1,1 -k2,2n -k3,3n |  gawk -F '\t' 'BEGIN{printf("#Target\ta1\ta2\tCoveron\tGene\tLength bp\tmax cover\taverage cover\taligned bp\n");}/^#/{next;}{chrom=$1;a1=$2;a2=$3;type=$4;if(type=="GENE"){g++;gnam[g]=$5;gch[g]=chrom;g1[g]=a1;g2[g]=a2;gz=0;if(a2>amax){amax=a2;top=g;}next;}if(type=="COVERON"){gz=0;for(i=g;i>=top;i--){if(gch[i]==chrom && a1<g2[i] && a2>g1[i]){gz=gnam[i];break;}}if(gz==0){ngz++;gz="New_" chrom "_"  ngz;}printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$5,gz,$6,$7,$8,$9);}}'  > tmp/WIGGLEGROUP/$run/$chrom/coverome2gene.$minCoveron.u.peaks

    \rm toto.sp[01]
    endif
end


done:
 exit 0



       
