#!bin/tcsh

set tConfig=$1
set run=$2


set rConfig=tmp/SortAlign/$run/rConfig

foreach lane (`cat Fastc/$run/LaneList`)
  echo "Fastc/$lane.fastc.gz $run" >> $rConfig
end

echo "numactl --cpunodebind=0 --membind=0 bin/sortalign -x tmp/SortAlign/targetIndex -I $rConfig --align --nB 50 --nA 25 -o tmp/SortAlign/$run/sali --max_threads 256 "
      numactl --cpunodebind=0 --membind=0 bin/sortalign -x tmp/SortAlign/targetIndex -I $rConfig --align --nB 50 --nA 25 -o tmp/SortAlign/$run/sali --max_threads 256 >  tmp/SortAlign/$run/sali.out
touch  tmp/SortAlign/$run/sali.done


  
