#!bin/tcsh -ef

set group=$1
set chrom=$2

# construct the combined wiggle of runs with low genomic contamination

# prepare the list of BF runs

set out_step="-out_step 10"
if ($?wiggle_step) then
   set out_step="-out_step $wiggle_step"
endif

if (! -d  tmp/WIGGLEGROUP/$group) mkdir tmp/WIGGLEGROUP/$group
set out=tmp/WIGGLEGROUP/$group/$chrom
if (! -d  $out) mkdir $out

if (! -e $out/m2.wg.done) then
# check that all runs are available

  foreach uu (u nu pp)
    foreach fr (f r  ELF ELR ERF ERR)
      if (-e $out/$fr.$uu.chrom.list) \rm  $out/$fr.$uu.chrom.list
      touch $out/$fr.$uu.chrom.list
      set ok=0
      foreach run (`cat MetaDB/$MAGIC/g2r | gawk -F '\t' '{if($1==g)print $2;}' g=$group | sort`)
        if (-e MetaDB/$MAGIC/WiggleDropEndList) then
          set drop=0
          if ($fr ==  ELF || $fr ==  ELR || $fr == ERF || $fr == ERR) then
            foreach run2 (`cat MetaDB/$MAGIC/WiggleDropEndList`)
              if ($run2 == $run) set drop=1
            end
            if ($drop == 1) continue 
          endif
	endif
        if (-e  tmp/SA/$run/wiggles/$run.$chrom.$uu.$fr.BF.gz) then
          echo -n "tmp/SA/$run/wiggles/$run.$chrom.$uu.$fr.BF.gz " >>   $out/$fr.$uu.chrom.list
            set ok=1
	  endif
        endif

      end
      echo "ok=$ok uu=$uu $fr"
      if ($ok == 1) then
          # contruct the combined wiggles
        if ( ! -e $out/R.chrom.$uu.$fr.BF.gz) then 
          echo  $out/$fr.$uu.chrom
          cat $out/$fr.$uu.chrom.list
          gunzip -c `cat $out/$fr.$uu.chrom.list` | bin/wiggle  -I BF -gzo -O BF $out_step  -o $out/R.chrom.$uu.$fr -cumul >&   $out/R.genes.$uu.$fr.cumul 
        endif
      endif
    end
    
    if (-e $out/R.chrom.$uu.ELF.BF.gz || -e  $out/R.chrom.$uu.ELR.BF.gz) then
       gunzip -c  $out/R.chrom.$uu.EL*.BF.gz  | bin/wiggle -I BF -O BF   $out_step -o $out/R.chrom.$uu.EL
       if (-e $out/R.chrom.$uu.EL.BF.gz) \rm $out/R.chrom.$uu.EL.BF.gz
       gzip  $out/R.chrom.$uu.EL.BF
    endif
    if (-e $out/R.chrom.$uu.ERF.BF.gz || -e  $out/R.chrom.$uu.ERR.BF.gz) then
       gunzip -c  $out/R.chrom.$uu.ER*.BF.gz  | bin/wiggle -I BF -O BF   $out_step -o $out/R.chrom.$uu.ER
       if (-e $out/R.chrom.$uu.ER.BF.gz) \rm $out/R.chrom.$uu.ER.BF.gz
       gzip  $out/R.chrom.$uu.ER.BF
    endif
  
  end


    echo "Construct the transcriptsEnds  $WG/$group"
    echo "  bin/wiggle  -transcriptsEnds tmp/$WG/$group/$chrom/R.chrom.u -gzi -I BF -O COUNT -o tmp/$WG/$group/$chrom/wg2a  -minCover 300 -wiggleRatioDamper 5"
            bin/wiggle  -transcriptsEnds tmp/$WG/$group/$chrom/R.chrom.u -gzi -I BF -O COUNT -o tmp/$WG/$group/$chrom/wg2a -minCover 300 -wiggleRatioDamper 5

endif

touch $out/wg2b.done
exit 0

