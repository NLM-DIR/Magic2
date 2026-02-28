#!bin/tcsh -ef

set group=$1
set chrom=$2
set unique=$3
set chrom1=$4

set out_step="-out_step 10"
if ($?wiggle_step) then
   set out_step="-out_step $wiggle_step"
endif

pushd tmp
foreach ss (f r ELF ELR ERF ERR)
   if (-e         TABIX/$group/R.$chrom.$unique.$ss.BF.gz) continue

   set iiiG=WIGGLEGROUP/$group/$chrom/R.chrom.$unique.$ss.BF.gz 
   set iiiR=WIGGLERUN/$group/$chrom/R.chrom.$unique.$ss.BF.gz 
   set iii=$iiiG
   if (! -e $iiiG && -e $iiiR) set iii=$iiiR
   if (! -e $iii) continue

   pushd  TABIX/$group
     ln -s ../../$iii R.$chrom.$unique.$ss.BF.gz 
   popd

   if (0) then
     if (-e         TABIX/$group/$chrom.$unique.$ss.tabix.gz) continue
     echo "gunzip -c $iii | bin/wiggle -I BF -O BV $out_step  |  gawk "'/^#/{next;}/^track/{next;}/^variableStep/{split($2,aa,"chrom=");zz=0;if(aa[2]==chrom)zz=1;next;}{if(zz==1)printf("%s\t%d\t%d\n",chrom,$1,$2);}'" chrom=$chrom1 | bgzip -c > TABIX/$group/$chrom.$unique.$ss.tabix.gz"
      gunzip -c $iii | bin/wiggle -I BF -O BV $out_step |  gawk '/^#/{next;}/^track/{next;}/^variableStep/{split($2,aa,"chrom=");zz=0;if(aa[2]==chrom)zz=1;next;}{if(zz==1)printf("%s\t%d\t%d\n",chrom,$1,$2);}' chrom=$chrom1 | bgzip -c > TABIX/$group/$chrom.$unique.$ss.tabix.gz
      tabix -s 1 -b 2 -e 2  tmp/TABIX/$group/$chrom.$unique.$ss.tabix.gz
    endif
end
popd

touch  tmp/TABIX/$group/$chrom.$unique.tabix.done

# usage tabix tata.gz $chom:40668-40700
# tabix TABX/Mixed.gonly.f.gz X:900668-1095700 | wc
# tabix TABX/Mixed.gonly.f.gz II:90668-195700 | wc
