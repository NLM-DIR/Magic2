#!bin/tcsh -f

set isRL=$1
set RL=$2
set H=$3

if ($H != H) then
  set H=""
  set HH=""
else
  set HH="-hierarchic"
endif
set step=10

set target=av
source scripts/target2target_class.txt
    
set mw=TARGET/MRNAS/mRNA2Wiggle.txt
if (-e mRNA2Wiggle.txt) set mw=mRNA2Wiggle.txt

if ($isRL == lane) then
  set lane=$RL
  if (-e tmp/COUNT/$lane.hits.gz) then
    set mw=TARGET/MRNAS/mRNA2Wiggle.txt
    if (-e mRNA2Wiggle.txt) set mw=mRNA2Wiggle.txt
    gunzip -c tmp/COUNT/$lane.hits.gz |  bin/wiggle -I BHIT -gzo -o tmp/MRNAWIGGLELANE$H/$lane -mapon $mw -ventilate $HH
    touch tmp/MRNAWIGGLELANE$H/$lane.done
  endif
endif


if ($isRL == run) then
  set run=$RL
  echo "# $run\tmRNA\tLength\tBegin\t0\t$step.." > tmp/MRNAWIGGLERUN$H/$run/mrna.f.wiggles.txt
  echo "# $run\tmRNA\tLength\tBegin\t0\t$step.." > tmp/MRNAWIGGLERUN$H/$run/mrna.r.wiggles.txt
  foreach mrna (`cat $mw | cut -f 2`)
    echo $mrna

    gunzip -c  tmp/MRNAWIGGLELANE$H/$run/*:$mrna.f.minerr0.hits.gz |  bin/wiggle -I BHIT -O BF -gzo -o tmp/MRNAWIGGLERUN$H/$run/$mrna.f -out_step $step
    gunzip -c  tmp/MRNAWIGGLELANE$H/$run/*:$mrna.r.minerr0.hits.gz |  bin/wiggle -I BHIT -O BF -gzo -o tmp/MRNAWIGGLERUN$H/$run/$mrna.r -out_step $step

    gunzip -c  tmp/MRNAWIGGLELANE$H/$run/*:$mrna.f.minerr0.hits.gz |  bin/wiggle -I BHIT -O BF -gzo -o tmp/MRNAWIGGLERUN$H/$run/$mrna.lengthCoverage.f -out_step $step -lengthCoverage
    gunzip -c  tmp/MRNAWIGGLELANE$H/$run/*:$mrna.r.minerr0.hits.gz |  bin/wiggle -I BHIT -O BF -gzo -o tmp/MRNAWIGGLERUN$H/$run/$mrna.lengthCoverage.r -out_step $step -lengthCoverage

    echo -n "Backward_plot\t$run\t$mrna" >> tmp/MRNAWIGGLERUN$H/$run/mrna.f.wiggles.txt
    echo -n "Backward_plot\t$run\t$mrna" >> tmp/MRNAWIGGLERUN$H/$run/mrna.r.wiggles.txt
    gunzip -c tmp/MRNAWIGGLERUN$H/$run/$mrna.f.BF.gz | gawk '/^fixedStep/{n1=step;i=index($3,"start=");if(i>0)n1=substr($3,7)+0;n2=n1-step;next;}{n2+=step;z[n2]=$1;nt+=step*$1;}END{printf("\t%d\t%d\t%d",n2,nt,n1);for(i=n2;i>=n1;i-=step)printf("\t%d",z[i]);printf("\n");}' step=$step >> tmp/MRNAWIGGLERUN$H/$run/mrna.f.wiggles.txt
    gunzip -c tmp/MRNAWIGGLERUN$H/$run/$mrna.r.BF.gz | gawk '/^fixedStep/{n1=step;i=index($3,"start=");if(i>0)n1=substr($3,7)+0;n2=n1-step;next;}{n2+=step;z[n2]=$1;nt+=step*$1;}END{printf("\t%d\t%d\t%d",n2,nt,n1);for(i=n2;i>=n1;i-=step)printf("\t%d",z[i]);printf("\n");}' step=$step >> tmp/MRNAWIGGLERUN$H/$run/mrna.r.wiggles.txt
 
    echo -n "Backward_lengthCoverage_plot\t$run\t$mrna" >> tmp/MRNAWIGGLERUN$H/$run/mrna.lengthCoverage.f.wiggles.txt
    echo -n "Backward_lengthCoverage_plot\t$run\t$mrna" >> tmp/MRNAWIGGLERUN$H/$run/mrna.lengthCoverage.r.wiggles.txt
    gunzip -c tmp/MRNAWIGGLERUN$H/$run/$mrna.lengthCoverage.f.BF.gz | gawk '/^fixedStep/{n1=step;i=index($3,"start=");if(i>0)n1=substr($3,7)+0;n2=n1-step;next;}{n2+=step;z[n2]=$1;nt+=step*$1;}END{printf("\t%d\t%d\t%d",n2,nt,n1);for(i=n2;i>=n1;i-=step)printf("\t%d",z[i]);printf("\n");}' step=$step >> tmp/MRNAWIGGLERUN$H/$run/mrna.lengthCoverage.f.wiggles.txt
    gunzip -c tmp/MRNAWIGGLERUN$H/$run/$mrna.lengthCoverage.r.BF.gz | gawk '/^fixedStep/{n1=step;i=index($3,"start=");if(i>0)n1=substr($3,7)+0;n2=n1-step;next;}{n2+=step;z[n2]=$1;nt+=step*$1;}END{printf("\t%d\t%d\t%d",n2,nt,n1);for(i=n2;i>=n1;i-=step)printf("\t%d",z[i]);printf("\n");}' step=$step >> tmp/MRNAWIGGLERUN$H/$run/mrna.lengthCoverage.r.wiggles.txt
 
   \rm  tmp/MRNAWIGGLERUN$H/$run/$mrna.[fr].BF.gz
  end

  touch tmp/MRNAWIGGLERUN$H/$run/run.done
endif

# ATTENTION the Run wiggles are drawn Backward, we want the group wiggles Forward

if ($isRL == group) then
  set group=$RL 
  foreach fr (f r)
    if (-e tmp/MRNAWIGGLEGROUP$H/$group/toto.$fr) \rm tmp/MRNAWIGGLEGROUP$H/$group/toto.$fr
    if (-e tmp/MRNAWIGGLEGROUP$H/$group/totoLC.$fr) \rm tmp/MRNAWIGGLEGROUP$H/$group/totoLC.$fr
    foreach run (`cat MetaDB/$MAGIC/_q2 | gawk '/Is_run/{if ($1 == g || $4 == g)print $1;}' g=$group | sort -u`)
      cat  tmp/MRNAWIGGLERUN$H/$run/mrna.$fr.wiggles.txt | sed -e 's/^Backward plot/Backward_plot/' |  gawk '/^#/{next}{print}' >> tmp/MRNAWIGGLEGROUP$H/$group/toto.$fr
      cat  tmp/MRNAWIGGLERUN$H/$run/mrna.lengthCoverage.$fr.wiggles.txt | sed -e 's/^Backward plot/Backward_plot/' |  gawk '/^#/{next}{print}' >> tmp/MRNAWIGGLEGROUP$H/$group/totoLC.$fr
    end
    sort -k 3,3 -k 2,2 tmp/MRNAWIGGLEGROUP$H/$group/toto.$fr > tmp/MRNAWIGGLEGROUP$H/$group/totos.$fr
    sort -k 3,3 -k 2,2 tmp/MRNAWIGGLEGROUP$H/$group/totoLC.$fr > tmp/MRNAWIGGLEGROUP$H/$group/totoLCs.$fr

    cat  tmp/MRNAWIGGLEGROUP$H/$group/totos.$fr ZZZZZ | gawk -F '\t' '{m=$3;if(old && m != old){printf("Forward_plot\t%d\t%d\t%s:%s\t%s\t%d",mx,mt,group,fr,old,0);for(x=step;x<=mx;x+=step){printf("\t%d",mm[x]);mm[x]=0;}if(mx>mxx)mxx=mx;printf("\n");mt=0;mx=0;}old=m;mt+=$5;x=$4+step;for(j=7;j<=NF;j++){x-=step;mm[x]+=$j;if (0) print $2,x,mm[x];}if($4>mx)mx=$4;}END{printf("0 Forward plot\tLength\tTotal\tGroup:strand\tmRNA\tBegin");for(x=0;x<=mxx;x+=10)printf("\t%d",x);printf("\n");}' group=$group fr=$fr step=$step | sort > tmp/MRNAWIGGLEGROUP$H/$group/mrna.$fr.wiggles.txt
    cat  tmp/MRNAWIGGLEGROUP$H/$group/totoLCs.$fr ZZZZZ | gawk -F '\t' '{m=$3;if(old && m != old){printf("Forward_plot\t%d\t%d\t%s:%s\t%s\t%d",mx,mt,group,fr,old,0);for(x=step;x<=mx;x+=step){printf("\t%d",mm[x]);mm[x]=0;}if(mx>mxx)mxx=mx;printf("\n");mt=0;mx=0;}old=m;mt+=$5;x=$4+step;for(j=7;j<=NF;j++){x-=step;mm[x]+=$j;if (0) print $2,x,mm[x];}if($4>mx)mx=$4;}END{printf("0 Forward plot\tLength\tTotal\tGroup:strand\tmRNA\tBegin");for(x=0;x<=mxx;x+=10)printf("\t%d",x);printf("\n");}' group=$group fr=$fr step=$step | sort > tmp/MRNAWIGGLEGROUP$H/$group/mrna.LengthCoverage.$fr.wiggles.txt

    #\rm tmp/MRNAWIGGLEGROUP$H/$group/toto*

  end
  touch tmp/MRNAWIGGLEGROUP$H/$group/group.done
endif

exit 0
# example
cat tmp/MRNAWIGGLEGROUPH/CH20a/mrna.f.wiggles.txt | grep LOC399942 | scripts/transpose | gawk '{n+=10;print n-60,$1;}' | more
