#!bin/tcsh -ef

set group=$1
set chrom=$2
set justMito=$3

# construct the combined wiggle of runs with low genomic contamination

# prepare the list of BF runs

set out_step="-out_step 10"
if ($?wiggle_step) then
   set out_step="-out_step $wiggle_step"
endif

set out=tmp/WIGGLEGROUP/$group/$chrom
if ($justMito == 1) set out=tmp/WIGGLEGROUP/$group/$chrom.step1
if (! -d  tmp/WIGGLEGROUP/$group) mkdir tmp/WIGGLEGROUP/$group
if (! -d  tmp/WIGGLEGROUP/$group/$chrom) mkdir tmp/WIGGLEGROUP/$group/$chrom
if ($justMito == 1 && ! -d tmp/WIGGLEGROUP/$group/$chrom.step1) mkdir tmp/WIGGLEGROUP/$group/$chrom.step1

if (! -e $out/wg2b.done) then
# check that all runs are available

  foreach uu (u nu pp)
    foreach fr (f r  ELF ELR ERF ERR)
      if (-e $out/$fr.$uu.chrom.list) \rm  $out/$fr.$uu.chrom.list
      if (-e $out/$fr.$uu.genes.list) \rm $out/$fr.$uu.genes.list
      touch $out/$fr.$uu.chrom.list
      set ok=0
      foreach run (`cat tmp/WIGGLEGROUP/$group/RunList`)
        if (-e MetaDB/$MAGIC/WiggleDropEndList) then
          set drop=0
          if ($fr ==  ELF || $fr ==  ELR || $fr == ERF || $fr == ERR) then
            foreach run2 (`cat MetaDB/$MAGIC/WiggleDropEndList`)
              if ($run2 == $run) set drop=1
            end
            if ($drop == 1) continue 
          endif
	endif
        if (-e      tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$uu.$fr.BF.gz) then
          echo -n " tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$uu.$fr.BF.gz " >>   $out/$fr.$uu.chrom.list
          set ok=1
        else
	  if (-e      tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr.BF.gz) then
            echo -n " tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr.BF.gz " >>   $out/$fr.$uu.chrom.list
            set ok=1
	  endif
        endif

        if ($Strategy != Exome && $justMito == 1) then
          if ($fr == f && -e tmp/SNP/$run/SpikeIn.plus.u.BV.gz) then
            echo -n " tmp/SNP/$run/SpikeIn.plus.u.BV.gz " >>   $out/$fr.$uu.chrom.list
            set ok=1
          endif
          if ($fr == r && -e tmp/SNP/$run/SpikeIn.minus.u.BV.gz) then
            echo -n " tmp/SNP/$run/SpikeIn.minus.u.BV.gz " >>   $out/$fr.$uu.chrom.list
            set ok=1
          endif
        endif
      end
      echo "ok=$ok uu=$uu"
      if ($ok == 1) then
          # contruct the combined wiggles
        if (-e  $out/$fr.$uu.genes.list) then
          if ($Strategy != Exome && $justMito == 0 && ! -e tmp/WIGGLEGROUP/$group/$chrom/R.genes.$uu.$fr.BF.gz) then
            echo  $out/$fr.$uu.genes
            gunzip -c `cat $out/$fr.$uu.genes.list` | bin/wiggle  -I BF -gzo -O BF $out_step  -o tmp/WIGGLEGROUP/$group/$chrom/R.genes.$uu.$fr -cumul >&   tmp/WIGGLEGROUP/$group/$chrom/R.genes.$uu.$fr.cumul
          endif
        endif
        cat $out/$fr.$uu.chrom.list
        echo
        if ( ! -e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.$fr.BF.gz) then 
          echo  $out/$fr.$uu.chrom
          cat $out/$fr.$uu.chrom.list
          echo aa
          gunzip -c `cat $out/$fr.$uu.chrom.list` | bin/wiggle  -I BV -gzo -O BF $out_step  -o  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.$fr -cumul >&   tmp/WIGGLEGROUP/$group/$chrom/R.genes.$uu.$fr.cumul 
        endif
      endif
    end
   
    if (-e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.f.BF.gz || -e  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.r.BF.gz) then
      echo "bb "
      ls -ls tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.[fr].BF.gz
      gunzip -c tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.[fr].BF.gz | bin/wiggle -I BF -O BF  $out_step -gzo -o tmp/WIGGLEGROUP/$group/$chrom/R.chrom.frns.$uu -cumul
    endif
    
    if (-e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ELF.BF.gz || -e  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ELR.BF.gz) then
       gunzip -c  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.EL*.BF.gz  | bin/wiggle -I BF -O BF   $out_step -o tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.EL
       if (-e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.EL.BF.gz) \rm tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.EL.BF.gz
       gzip  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.EL.BF
    endif
    if (-e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ERF.BF.gz || -e  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ERR.BF.gz) then
       gunzip -c  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ER*.BF.gz  | bin/wiggle -I BF -O BF   $out_step -o tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ER
       if (-e tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ER.BF.gz) \rm tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ER.BF.gz
       gzip  tmp/WIGGLEGROUP/$group/$chrom/R.chrom.$uu.ER.BF
    endif
  
  end


    set WG=WIGGLEGROUP
    set tt=tmp/$WG/$group/projected.stranding.tsf
    echo '#' > $tt
    echo '#' > $tt.1
      set ok=0
      foreach chrom ($chromSetAll)
        set sBV=TARGET/CHROMS/$chrom.stranding.BV
	if (! -e $sBV) continue
        set wf=tmp/$WG/$group/$chrom/R.chrom.u.f.BF.gz
        set wr=tmp/$WG/$group/$chrom/R.chrom.u.r.BF.gz
	if (! -e $wf || ! -e $wr) continue
        bin/wiggle -i $wf -I BF -O BV | gawk '{n=$1+0;if(n>0)printf("%d\tf\tii\t%d\t0\n",$1,$2);}' > $tt.a
        bin/wiggle -i $wr -I BF -O BV | gawk '{n=$1+0;if(n>0)printf("%d\tr\tii\t0\t%d\n",$1,$2);}' > $tt.b
        cat $tt.a $tt.b | bin/tsf -I tsf -O tsf -s ns >$tt.c

        if (0) then   # do this once using run==the stranded group
          cat $tt.c | gawk -F '\t' '/^#/{next;}{x=$1;p=($4+20)/($4+$5+20);q=int(100*p);printf("%d\t%d\n",x,q);}' > $tt.d
	  cat $tt.d | gawk '{p=$2;if(p>90 || p < 10) print;}' > $sBV
        endif
	
        cat $sBV $tt.c | gawk -F '\t' '/^#/{next}{if(NF==2){pp[$1]=2*($2-50);next;}}{p=pp[$1];if(p>97){sp+=$4;sm+=$5;};if(p<-97){sp+=$5;sm+=$4;}}END{printf("%s\t%s\tiif\t%d\t%d\t%.2f\n",run,chrom,sp,sm,100*sp/(sp+sm+.001));}' run=$group chrom=$chrom >> $tt.1
	set ok=1
	\rm $tt.a $tt.b $tt.c
      end
    if ($ok == 1) then
      bin/tsf -I tsf -O tsf  -i $tt.1 -s any | gawk -F '\t' '/^#/{next;}{sp=$4;sm=$5;printf("%s\tP_genome\tiif\t%d\t%d\t%.3f\n",run,sp,sm,100*sp/(sp+sm+.001));}' run=$group >> $tt
      \rm $tt.1
    endif


    echo "Construct the transcriptsEnds  $WG/$group"
    echo "  bin/wiggle  -transcriptsEnds tmp/$WG/$group/$chrom/R.chrom.u -gzi -I BF -O COUNT -o tmp/$WG/$group/$chrom/wg2a  -minCover 300 -wiggleRatioDamper 5"
            bin/wiggle  -transcriptsEnds tmp/$WG/$group/$chrom/R.chrom.u -gzi -I BF -O COUNT -o tmp/$WG/$group/$chrom/wg2a -minCover 300 -wiggleRatioDamper 5

endif

touch $out/wg2b.done
exit 0

