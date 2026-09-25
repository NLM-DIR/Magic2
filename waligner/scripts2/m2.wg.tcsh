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

set totoSS=tmp/Peaks/$MAGIC.strand_shifts.max.txt
# check that all runs are available   
  foreach uu (u nu pp)
    foreach fr (f r ELF ELR ERF ERR)
      set out2=tmp/$WRG/$group/$group.$chrom.$uu.$fr
      set wigList=$out2.$fr.wigList
      if (-e $out2.az) continue
      if (-e $wigList ) \rm  $wigList
      touch $wigList
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
	set BBB=tmp/SA/$run/wiggles/$run.$chrom.$uu.$fr.az
	if (! -e $BBB) set BBB=tmp/WIGGLERUN/$run/$run.$chrom.$uu.$fr.az
	if (! -e $BBB) set BBB=tmp/WIGGLEGROUP/$run/$run.$chrom.$uu.$fr.az	
        set dx=0
	if (-e $totoSS) then
          set dx=`cat $totoSS | gawk '{if($1==run)dx=int($2/2);}END{dx2=0;if(fr=="f")dx2=dx;if (fr=="r")dx2=-dx;print dx2;}' fr=$fr run=$run`
        endif
	if (-e $BBB) then
          echo "$BBB\t$dx" >>  $wigList
            set ok=1
	  endif
        endif
      end
      echo "ok=$ok uu=$uu $fr"            # contruct the combined wiggles
      if ($ok == 1) then
	bin/wiggle -f $wigList -I AZ -O AZ $out_step  -o $out2 -cumul
      endif
    end
  end

touch $out/wg2b.done
exit 0

