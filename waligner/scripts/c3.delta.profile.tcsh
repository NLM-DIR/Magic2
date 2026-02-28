#!bin/tcsh -f
set run=$1

if (! -e  tmp/Profiles/$run/readsBeforeAlignment.delta.txt) then 
  echo "preparing tmp/Profiles/$run/readsBeforeAlignment.delta.txt"
  gunzip -c Fastc/$run/*.*.fastc.gz | bin/dna2dna  -I fastc -minEntropy $minEntropy  -plot 1000000 -o tmp/Profiles/$run/readsBeforeAlignment.1M 

  set n=0
  if (-e Fastc/$run/Max_probe_length) set n=`cat Fastc/$run/Max_probe_length`
  set delta=1
  if ($n >= 1000) set delta=10

  cat tmp/Profiles/$run/readsBeforeAlignment.1M.profile.txt | gawk -F '\t' '/^#/{next}{k=int(($1+delta - 1)/delta) ; if(k>900)k=900;nn[k] += $2;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%d\t%d\n", run, delta*i, nn[i]+0);}' delta=$delta run=$run >  tmp/Profiles/$run/readsBeforeAlignment.delta.txt
  \rm tmp/Profiles/$run/readsBeforeAlignment.1M.*
endif


if (0 && ! -e  tmp/Profiles/$run/readsAfterAlignment.delta.txt) then 
  echo "preparing  tmp/Profiles/$run/readsAfterAlignment.delta.txt"
  foreach lane (`cat Fastc/$run/LaneList`)
    bin/bestali -aliProfile -i tmp/COUNT/$lane.hits.gz -gzo -o tmp/COUNT/$lane 
  end
  gunzip -c tmp/COUNT/$run/*.aliProfile.gz | gawk '/^Aligned_length/{n[$2]+=$3;}END{for(k in n)printf("%d\t%d\n",k,n[k]);}' | sort -k 1n >  tmp/Profiles/$run/readsAfterAlignment.1M.profile.txt 

  set n=0
  if (-e Fastc/$run/Max_probe_length) set n=`cat Fastc/$run/Max_probe_length`
  set delta=1
  if ($n > 1000) set delta=10


  cat tmp/Profiles/$run/readsAfterAlignment.1M.profile.txt | gawk -F '\t' '/^#/{next}{k=int(($1+delta - 1)/delta) ; if(k>900)k=900;nn[k] += $2;if(k>kMax)kMax =k;}END{for (i = 0 ; i <= kMax ; i++) printf ("%s\t%d\t%d\n", run, delta*i, nn[i]+0);}' delta=$delta run=$run >  tmp/Profiles/$run/readsAfterAlignment.delta.txt
  \rm tmp/Profiles/$run/readsAfterAlignment.1M.profile.txt
endif
