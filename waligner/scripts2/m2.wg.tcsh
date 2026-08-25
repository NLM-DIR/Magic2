#!bin/tcsh -ef

set WRG=$1
set group=$2
set chrom=$3

# construct the combined wiggle of runs with low genomic contamination

# prepare the list of BF runs

set out_step="-out_step 10"
if ($?wiggle_step) then
   set out_step="-out_step $wiggle_step"
endif

if (! -d  tmp/$WRG/$group) mkdir tmp/$WRG/$group
set out=tmp/$WRG/$group/$chrom
if (! -d  $out) mkdir $out


if (! -e $out/m2.wg.doneZZZ) then
# check that all runs are available
# 
  foreach uu (u nu pp)
    foreach fr (f r  ELF ELR ERF ERR)
      if (-e $out/$fr.$uu.chrom.list) \rm  $out/$fr.$uu.chrom.list
      touch $out/$fr.$uu.chrom.list
      set ok=0
      foreach run (`cat MetaDB/$MAGIC/RunWiggleList MetaDB/$MAGIC/g2r | gawk -F '\t' '{if (NF==1 && $1==g)print $1; if(NF>=2 && $1==g)print $2;}' g=$group | sort -u`)
        if (-e MetaDB/$MAGIC/WiggleDropEndList) then
          set drop=0
          if ($fr ==  ELF || $fr ==  ELR || $fr == ERF || $fr == ERR) then
            foreach run2 (`cat MetaDB/$MAGIC/WiggleDropEndList`)
              if ($run2 == $run) set drop=1
            end
            if ($drop == 1) continue 
          endif
	endif
	echo BBB tmp/SA/$run/wiggles/$run.$chrom.$uu.$fr.BF.gz
	ls -ls  tmp/SA/$run/wiggles/$run.$chrom.$uu.$fr.az
        if (-e  tmp/SA/$run/wiggles/$run.$chrom.$uu.$fr.az) then
          echo "tmp/SA/$run/wiggles/$run.$chrom.$uu.$fr.az " >>   $out/$fr.$uu.chrom.list
            set ok=1
	  endif
        endif

      end
      echo "ok=$ok uu=$uu $fr"
      if ($ok == 1) then
          # contruct the combined wiggles
        if ( ! -e $out/R.chrom.$uu.$fr.az) then 
          echo  $out/$fr.$uu.chrom
          cat $out/$fr.$uu.chrom.list
	  set nr=`wc -l $out/$fr.$uu.chrom.list | gawk '{print $1}'`
	  echo "CCCC nr=$nr"
	  if ($nr > 1) then
            cat `cat $out/$fr.$uu.chrom.list` | bin/wiggle  -I AZ -O AZ $out_step  -o $out/R.chrom.$uu.$fr -cumul >&   $out/R.genes.$uu.$fr.cumul
	  endif
	  if ($nr == 1) then
	    pushd $out
	      pwd
	      ln -s ../../../SA/$run/wiggles/$run.$chrom.$uu.$fr.az R.chrom.$uu.$fr.az
              ln -s ../../../SA/$run/wiggles/$run.$chrom.$uu.$fr    R.chrom.$uu.$fr
	    popd
	  endif
        endif
      endif
    end


# evaluate the coverage per threshold

  if ($uu == u && ($fr == f || $fr == r) && -e  $out/R.chrom.$uu.$fr.az) then

    foreach cover (1 2 5 10 20 50 100 200 500 1000)

      if (! -e  $out/$chrom.pseudoexon.$cover.$uu.txt) then
        bin/wiggle -i  $out/$chrom.frns.$uu.az  -I BF -O BV -minCover $cover | gawk -F '\t' '/^#/{next;}/^track/{next;}/^variableStep/{if(oldx2)printf("%s\t%d\t%d\t%f\n",chrom,oldx1,oldx2,oldy);nn+=oldx2 - oldx1 + 1;oldx1=0;oldx2=0;oldy=0;if(x2==0)x2=1;printf("#Chromosome %s Length %.3f Mb, Usable %.3f Mb Covered over %d : %.3f%%\n",chrom, x2/1000000,nn/1000000,minCover,100*nn/x2);nn=0;x2=0;chrom=$1;oldy=0;gsub(/variableStep chrom=/,"",chrom);next;} {n++;if(n<4)next;x1=$1-10;x2=$1+10;if($2>oldy)oldy=$2;if(x1>oldx2+50){if(oldx2)printf("%s\t%d\t%d\t%f\n",chrom,oldx1,oldx2,oldy);nn+=oldx2 - oldx1 + 1;oldx1=x1;oldy=0;}oldx2=x2;if($2>oldy)oldy=$2;}END{if(oldx2)printf("%s\t%d\t%d\t%f\n",chrom,oldx1,oldx2,oldy);nn+=oldx2 - oldx1 + 1;if(x2==0)x2=1;printf("#Chromosome %s Length %.3f Mb, Usable %.3f Mb Covered over %d : %.3f%%\n",chrom, x2/1000000,nn/1000000,minCover,100*nn/x2);}'  minCover=$cover > $out/$chrom.pseudoexon.$cover.$uu.txt
      endif
    end

  endif

    if (-e $out/R.chrom.$uu.ELF.az || -e  $out/R.chrom.$uu.ELR.az) then
       cat  $out/R.chrom.$uu.EL*.az  | bin/wiggle -I AZ -O AZ   $out_step -o $out/R.chrom.$uu.EL
    endif
    if (-e $out/R.chrom.$uu.ERF.az || -e  $out/R.chrom.$uu.ERR.az) then
       cat  $out/R.chrom.$uu.ER*.az  | bin/wiggle -I AZ -O AZ  $out_step -o $out/R.chrom.$uu.ER
    endif
  
  end

    if (-e tmp/$WRG/$group/$chrom/R.chrom.u.ERF.az) then
      echo "Construct the transcriptsEnds  $WRG/$group"
      echo "  bin/wiggle  -transcriptsEnds tmp/$WRG/$group/$chrom/R.chrom.u -I AZ -O COUNT -o tmp/$WRG/$group/$chrom/wg2a  -minCover 300 -wiggleRatioDamper 5"
              bin/wiggle  -transcriptsEnds tmp/$WRG/$group/$chrom/R.chrom.u -I AZ -O COUNT -o tmp/$WRG/$group/$chrom/wg2a -minCover 300 -wiggleRatioDamper 5
    endif
endif

touch $out/wg2b.done
exit 0

