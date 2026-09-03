#!bin/tcsh -ef

set WRG=$1
set group=$2
set chrom=$3
set out=tmp/$WRG/$group/$chrom
if (! -d  $out) mkdir $out

# construct the combined run wiggle from sublibs, or group wiggle form runs

set out_step="-out_step 10"
if ($?wiggle_step) then
   set out_step="-out_step $wiggle_step"
endif


# check that all runs are available   
  foreach uu (u nu pp)
    foreach fr (f r ELF ELR ERF ERR)
      if (-e $out/$fr.$uu.chrom.list) \rm  $out/$fr.$uu.chrom.list
      touch $out/$fr.$uu.chrom.list
      set ok=0
      foreach run (`cat MetaDB/$MAGIC/r2sublib MetaDB/$MAGIC/g2r | gawk -F '\t' '{if ($1==g)print $2;}' g=$group | sort -u`)
        if (-e MetaDB/$MAGIC/WiggleDropEndList) then
          set drop=0
          if ($fr ==  ELF || $fr ==  ELR || $fr == ERF || $fr == ERR) then
            foreach run2 (`cat MetaDB/$MAGIC/WiggleDropEndList`)
              if ($run2 == $run) set drop=1
            end
            if ($drop == 1) continue 
          endif
	endif
        set out2=tmp/$WRG/$group/$group.$chrom.$uu.$fr
	set BBB=tmp/SA/$run/wiggles/$run.$chrom.$uu.$fr.az
        if (-e $BBB) then
          echo $BBB >>  $out2.chrom.list
            set ok=1
	  endif
        endif
	set BBB=tmp/WIGGLERUN/$run.$chrom.$uu.$fr.az
        if (-e $BBB) then
          echo $BBB >>   $out2.chrom.list
            set ok=1
	  endif
        endif
      end
      echo "ok=$ok uu=$uu $fr"
      if ($ok == 1) then
          # contruct the combined wiggles
	bin/wiggle -f $out2.chrom.list -I AZ -O AZ $out_step  -o $out2 -cumul >&   $out/R.genes.$uu.$fr.cumul
      endif
    end
  end

touch $out/wg2b.done
exit 0

