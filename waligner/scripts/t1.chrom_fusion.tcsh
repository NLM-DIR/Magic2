#!bin/tcsh -f

set phase=$1
set run=$2
setenv target genome

if ($phase == t1) goto phaseT1
exit (0)

##### Phase 1 collect the chrom_fusions
phaseT1:
  # collate the chrom fusions
  source scripts/target2target_class.txt
        
  set minAli=0

  if (! -e  tmp/ChromFusion/$run/t1.chrom_fusion.$target.tsf.gz1) then
    zcat tmp/COUNT/$run/*.hits.gz | gawk -F '\t' '/^#/{next;}{if($8 != "Z_genome")next;if($7-$6<minAli)next;r=$1;gsub(/>/,"",r);gsub(/</,"",r);ch=$11;if(r==oldr && ch != oldch){c1=ch;c2=oldch;if(c1>c2){c0=c1;c1=c2;c2=c0;}nn[c1 "___" c2]+=$3;} oldr=r;oldch=ch;}END{for(c in nn)printf("%s\t%s\ti\t%d\n",c,run,nn[c]);}' minAli=$minAli run=$run > tmp/ChromFusion/$run/t1.chrom_fusion.$target.tsf
    gzip tmp/ChromFusion/$run/t1.chrom_fusion.$target.tsf
  endif

goto phaseLoop

phaseLoop:
 echo done

