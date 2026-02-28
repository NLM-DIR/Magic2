#!/bin/tcsh -f

set phase=$1
set run=$2
set lane=$3

echo hello from virus.tcsh $phase

if ($phase == vir2) goto phaseVir2

##########################################################################################################
## phase vir1

if ($phase == vir1) then

# VIRUS
  foreach tt (v b t)
    if ($tt == v) then
      set VV=VIRUS
      set vv=virus
      set tc=v_virus
    endif
    if ($tt == b) then
      set VV=BACTERIA
      set vv=bacteria
      set tc=b_bacteria
    endif
    if ($tt == t) then
      set VV=BACTERIA
      set vv=transposon
      set tc=D_transposon
    endif

    if (-e tmp/$VV/$run/vir1.$vv.count) continue 

    # extract the $vv hits from the COUNT/best hits tables
    echo ' ' >  tmp/$VV/$run/_f.$vv
    set ok=0
    foreach lane (`cat Fastc/$run/LaneList`)
      if (-e tmp/COUNT/$lane.hits.gz && ! -e tmp/$VV/$lane.$tc.hits.gz ) then
        zcat  tmp/COUNT/$lane.hits.gz | gawk -F '\t' '{if($8 == tc)print}' tc=$tc | gzip > tmp/$VV/$lane.$tc.hits.gz
      endif
      if (-e tmp/$VV/$lane.$tc.hits.gz ) then
        echo tmp/$VV/$lane.$tc.hits.gz  >>  tmp/$VV/$run/_f.$vv
        if ($vv == bacteria && -e tmp/$VV/$lane.clean_hits.gz) echo tmp/$VV/$lane.clean_hits.gz >>  tmp/$VV/$run/_f.$vv
        set ok=1
      endif
    end

    if ($ok == 1) then
      bin/bestali -autoH -target_class $tc -inFileList  tmp/$VV/$run/_f.$vv -run $run > tmp/$VV/$run/vir1.$vv.count 
      bin/bestali -inFileList  tmp/$VV/$run/_f.$vv -target_class $tc -run $run -countBest | gawk '/^HITS/{print}/^MULT/{print}' | grep $vv  > tmp/$VV/$run/vir1.$vv.HitMult.1
      cat tmp/$VV/$run/vir1.$vv.HitMult.1 | gawk -F '\t' '/^HITS/{z=$1"\t"$2"\t"$3 ;i=z2i[z]+0;if(i==0){iMax++;i=iMax;z2i[z]=i;i2z[i]=z;}for(j=4;j<=NF;j++)n[i,j]+=$j;if(nf<NF)nf=NF;}END{for(i=1;i<=iMax;i++){printf("%s",i2z[i]);for(j=4;j<=nf;j++)printf("\t%d",0+n[i,j]);printf("\n");}}' > tmp/$VV/$run/vir1.$vv.HitMult
      cat tmp/$VV/$run/vir1.$vv.HitMult.1 | gawk -F '\t' '/^MULT/{z=$1"\t"$2;i=z2i[z]+0;if(i==0){iMax++;i=iMax;z2i[z]=i;i2z[i]=z;}for(j=3;j<=NF;j++)n[i,j]+=$j;if(nf<NF)nf=NF;}END{for(i=1;i<=iMax;i++){printf("%s",i2z[i]);for(j=3;j<=nf;j++)printf("\t%d",0+n[i,j]);printf("\n");}}' >> tmp/$VV/$run/vir1.$vv.HitMult
      \rm tmp/$VV/$run/vir1.$vv.HitMult.1
    endif
  end

  goto phaseLoop
endif

##########################################################################################################
# HACK : hand realignment

if ($phase == bact1) then
  bin/clipalign -best -i Fastc/$lane.fastc.gz -t TARGET/Targets/$species.bacteria.fasta.gz  -maxHit 30 -minEntropy 30 -seedLength 16 -probeMinLength 100  -clipN 3 -minAli 60 -targetBonus -30 -target_class b_bacteria -strategy Genome  -previousScore  tmp/COUNT/$lane.best_score.gz -gzo  -o tmp/BACTERIA/$lane

  goto phaseLoop
endif
  
######

# HACK realign because we changed the bacteria list 2020-04-20
if (0) then
  foreach run (`cat MetaDB/$MAGIC/RunList `)
    if (! -d tmp/BACTERIA/$run) mkdir tmp/BACTERIA/$run
    foreach lane (`cat Fastc/$run/LaneList`)
      if (! -e tmp/BACTERIA/$lane.hits.gz) then
        scripts/submit tmp/BACTERIA/$lane.bact "scripts/virus.tcsh bact1 $run $lane"
      endif
    end
  end

  foreach run (`cat MetaDB/$MAGIC/RunList `)
    if (-d tmp/BACTERIA/$run) then
      ls  tmp/BACTERIA/$run/*.hits.gz >  tmp/BACTERIA/$run/_f
      set n=`cat tmp/BACTERIA/$run/_f | wc -l` 
      if ($n > 0 && ! -e   tmp/BACTERIA/$run/bact1.virus.count) then
        bin/bestali -autoH -target_class b_bacteria -inFileList  tmp/BACTERIA/$run/_f -run $run > tmp/BACTERIA/$run/bact1.virus.count 
      endif
    endif
  end

# HACK rattrapage des scores (because by error 20 2020-04-20, i had not used the COUNT/beast_score

  foreach lane ( `cat MetaDB/MM/LaneList ` )
    if (-e  tmp/BACTERIA/$lane.hits.gz  && ! -e   tmp/BACTERIA/$lane.clean_hits.gz) then
      scripts/submit tmp/BACTERIA/$lane.hack2 "scripts/hack.tcsh $lane"
    endif
  end

endif

########################################################################################################################
## vir2 report

phaseVir2:

if (-d VirusDB && ! -e  VirusDB/virus.metadata.txt) then
  bin/tace VirusDB <<EOF
    query find sequence Virus
    select -o VirusDB/virus.metadata.txt  type,s,ln,t,a from type="Virus", s in @, ln in s->length, t in s->title, a in s->accession
    query find sequence microbes
    select -o VirusDB/bacteria.metadata.txt  type,s,ln,t,a from type="Microbes", s in @, ln in s->length, t in s->title, a in s->accession
    query find sequence Transposon
    select -o VirusDB/transposon.metadata.txt  type,s,ln,t,a from type="Transposon", s in @, ln in s->length, t in s->title, a in s->accession
    query find sequence sp_species && Virus
    select -o VirusDB/virus_groups.metadata.txt  type,g,ln,t,s from type="Virus", g in @, ln in g->length, t in g->title, s in g->sp_species
    query find sequence sp_species && Microbes
    select -o VirusDB/bacteria_groups.metadata.txt  type,g,ln,t,s from type="Microbes", g in @, ln in g->length, t in g->title, s in g->sp_species
EOF

endif

  set totov=tmp/VIRUS/$MAGIC
  set totob=tmp/BACTERIA/$MAGIC
  set totot=tmp/BACTERIA/$MAGIC.tr
## associate the Accession to its full title using  VirusDB/virus.metadata.txt which is exported in metadata.tcsh
  cat MetaDB/$MAGIC/runs.ace MetaDB/$MAGIC/ali.ace > $totov.info.ace
  cat VirusDB/virus.metadata.txt | gawk -F '\t' '{printf("Gene \"%s\"\nTitle \"%s\"\nLength %d\nTargeted\nGeneId %s\n\n",$2,$4,$3,$5);}' >> $totov.info.ace
  cat MetaDB/$MAGIC/runs.ace MetaDB/$MAGIC/ali.ace > $totob.info.ace
  cat VirusDB/bacteria.metadata.txt | gawk -F '\t' '{printf("Gene \"%s\"\nTitle \"%s\"\nLength %d\nTargeted\nGeneId %s\n\n",$2,$4,$3,$5);}' >> $totob.info.ace
  cat MetaDB/$MAGIC/runs.ace MetaDB/$MAGIC/ali.ace > $totot.info.ace
  cat VirusDB/transposon.metadata.txt | gawk -F '\t' '{printf("Gene \"%s\"\nTitle \"%s\"\nLength %d\nTargeted\nGeneId %s\n\n",$2,$4,$3,$5);}' >> $totot.info.ace
  cat VirusDB/virus_groups.metadata.txt     | gawk -F '\t' '{printf("Gene_group \"%s\"\nTitle \"%s\"\nLength %d\nGenes %s\n%s\nAccession %s\n\n",$2,$4,$3,$5,$1,$5);}' > tmp/VIRUS/groups.ace
  cat VirusDB/bacteria_groups.metadata.txt | gawk -F '\t' '{printf("Gene_group \"%s\"\nTitle \"%s\"\nLength %d\nGenes %s\n%s\nAccession %s\n\n",$2,$4,$3,$5,$1,$5);}' > tmp/BACTERIA/groups.ace

## create a count file in the style expected by geneindex.c and collate the relevant runs in $toto.777

  echo ' ' > $totov.counts.ace
  echo ' ' > $totob.counts.ace
  echo ' ' > $totot.counts.ace
  foreach run (`cat MetaDB/$MAGIC/RunList `)
    foreach target (virus)
      if (-e  tmp/VIRUS/$run/vir1.$target.count) then
        cat  tmp/VIRUS/$run/vir1.$target.count | sed -e 's/NC_045512.2\t/NC_045512\t/g' | gawk  '/^#/{next;}{if($7<0)$7=int(($7+2**32)/$4);g=$2;n=split(g,aa,"|");if(substr(g,1,3)=="gi|"  && n>=5){split(aa[4],bb,".");g=bb[1];}printf("Gene \"%s\"\nRun_U \"%s\" 0.00 %d seqs %d tags %f kb %d raw_tags\n\n",g,$1,$3,$3,$7*$3/1000.0,$4);}' >> $totov.counts.ace
      endif
    end
    foreach target (bacteria)
      if (-e  tmp/BACTERIA/$run/vir1.$target.count) then
        cat  tmp/BACTERIA/$run/vir1.$target.count | gawk  '/^#/{next;}{if($7<0)$7=int(($7+2**32)/$4);g=$2;n=split(g,aa,"|");if(substr(g,1,3)=="gi|" && n>=5){split(aa[4],bb,".");g=aa[4];}printf("Gene \"%s\"\nRun_U \"%s\" 0.00 %d seqs %d tags %f kb %d raw_tags\n\n",g,$1,$3,$3,$7*$3/1000,0,$4);}' >> $totob.counts.ace
      endif
    end
    foreach target (transposon)
      if (-e  tmp/BACTERIA/$run/vir1.$target.count) then
        cat  tmp/BACTERIA/$run/vir1.$target.count | gawk  '/^#/{next;}{if($7<0)$7=int(($7+2**32)/$4);g=$2;n=split(g,aa,"|");if(substr(g,1,3)=="gi|" && n>=5){split(aa[4],bb,".");g=aa[4];}printf("Gene \"%s\"\nRun_U \"%s\" 0.00 %d seqs %d tags %f kb %d raw_tags\n\n",g,$1,$3,$3,$7*$3/1000.0,$4);}' >> $totot.counts.ace
      endif
    end
  end

  tace MetaDB <<EOF
    query find project IS $MAGIC
    select -o MetaDB/$MAGIC/Mb_in_genes_with_GeneId_minus_high_genes.txt r,m from p in @, r in p->run, a in r->ali, m in a->Mb_in_genes_with_GeneId_minus_high_genes[2] where m
EOF
  set MG=MetaDB/$MAGIC/Mb_in_genes_with_GeneId_minus_high_genes.txt

# we wish to sort the runs, then the group, but we add the sublibs so that their count is parsed, but anyway geneindex.c considers sublibs as private
  cat  MetaDB/$MAGIC/RunsListSorted   > tmp/VIRUS/$MAGIC.GroupsRunsListSorted
  bin/geneindex -export btm -deepGene $totov.counts.ace -u -runList tmp/VIRUS/$MAGIC.GroupsRunsListSorted -runAce $totov.info.ace -o $totov  -geneGroup  tmp/VIRUS/groups.ace
  cat  MetaDB/$MAGIC/RunsListSorted    > tmp/BACTERIA/$MAGIC.GroupsRunsListSorted
  bin/geneindex -export btm -deepGene $totob.counts.ace -u -runList tmp/BACTERIA/$MAGIC.GroupsRunsListSorted -runAce $totob.info.ace -o $totob  -geneGroup  tmp/BACTERIA/groups.ace
  bin/geneindex -export btm -deepGene $totot.counts.ace -u -runList tmp/BACTERIA/$MAGIC.GroupsRunsListSorted -runAce $totot.info.ace -o $totot

  if (-d RESULTS/VIRUS && ! -d RESULTS/Microbiome) mv RESULTS/VIRUS RESULTS/Microbiome
  if (! -d RESULTS/Microbiome) mkdir RESULTS/Microbiome

  set toto2a=RESULTS/Microbiome/$MAGIC.virus_bacteria.reads_aligned_per_run.txt
  echo -n "## $toto2a " > $toto2a
  cat $totov.gene_group.reads_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2a
  echo >> $toto2a
  cat $toto2a | gawk '/#Run/{print;last;}' > $totov.999
 
 cat $totov.gene_group.reads_aligned_per_gene.txt  | gawk '/^#/{next;}{print}' | scripts/tab_sort -k 8nr >> $totov.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ  ZZZZZ VirusDB/virus.metadata.txt ZZZZZ $totov.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Virus | scripts/tab_sort -k 8,8nr >> $toto2a
  echo "\n\n" >> $toto2a
  cat $totov.reads_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2a
  echo >> $toto2a

  cat $toto2a | gawk '/#Run/{print;last;}' > $totov.999
  cat $totov.reads_aligned_per_gene.txt  | gawk '/^#/{next;}{print}' | scripts/tab_sort -k 8nr >> $totov.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/virus.metadata.txt ZZZZZ $totov.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Virus | scripts/tab_sort -k 8,8nr >> $toto2a
  echo "\n\n" >> $toto2a

  cat $totob.gene_group.reads_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2a
  cat $toto2a | gawk '/#Run/{print;last;}' > $totob.999
  cat $totob.gene_group.reads_aligned_per_gene.txt  | gawk '/^#/{next;}{print}' | scripts/tab_sort -k 8nr >> $totob.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/bacteria.metadata.txt ZZZZZ $totob.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Microbes | scripts/tab_sort -k 8,8nr >> $toto2a
  echo "\n\n" >> $toto2a

  cat $totob.reads_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2a
  cat $toto2a | gawk '/#Run/{print;last;}' > $totov.999
  cat $totob.reads_aligned_per_gene.txt  | gawk '/^#/{next;}{print}' | scripts/tab_sort -k 8nr >> $totob.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/bacteria.metadata.txt ZZZZZ $totob.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Microbes | scripts/tab_sort -k 8,8nr >> $toto2a
  echo "\n\n" >> $toto2a

  cat $totot.reads_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2a
  cat $toto2a | gawk '/#Run/{print;last;}' > $totov.999
  cat $totot.reads_aligned_per_gene.txt  | gawk '/^#/{next;}{print}' | scripts/tab_sort -k 8nr >> $totot.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/transposon.metadata.txt ZZZZZ $totot.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Transposon | scripts/tab_sort -k 8,8nr >> $toto2a
  echo "\n\n" >> $toto2a

  set toto2b=RESULTS/Microbiome/$MAGIC.virus_bacteria.reads_aligned_per_target_per_million_raw_reads.txt 
  echo -n "## $toto2b : " > $toto2b
  cat $totov.gene_group.reads_aligned_per_gene_per_million_raw_reads.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2b
  echo >> $toto2b
  cat $toto2b | gawk '/#Run/{print;last;}' > $totov.999
  cat $totov.gene_group.reads_aligned_per_gene_per_million_raw_reads.txt  | gawk '/^#/{next;}{print}' | scripts/tab_sort -k 8nr >> $totov.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/virus.metadata.txt ZZZZZ $totov.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Virus | scripts/tab_sort -k 8,8nr >> $toto2b
  echo "\n\n" >> $toto2b

  cat $totov.reads_aligned_per_gene_per_million_raw_reads.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2b
  echo >> $toto2b
  cat $toto2b | gawk '/#Run/{print;last;}' > $totov.999
  cat $totov.reads_aligned_per_gene_per_million_raw_reads.txt  | gawk '/^#/{next;}{print}' | scripts/tab_sort -k 8nr >> $totov.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/virus.metadata.txt ZZZZZ $totov.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Virus | scripts/tab_sort -k 8,8nr >> $toto2b
  echo "\n\n" >> $toto2b

  cat $totob.gene_group.reads_aligned_per_gene_per_million_raw_reads.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2b
  cat $toto2b | gawk '/#Run/{print;last;}' > $totov.999
  cat $totob.gene_group.reads_aligned_per_gene_per_million_raw_reads.txt  | gawk '/^#/{next;}{print}' | scripts/tab_sort -k 8nr >> $totob.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/bacteria.metadata.txt ZZZZZ $totob.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Microbes | scripts/tab_sort -k 8,8nr >> $toto2b
  echo "\n\n" >> $toto2b
 
 cat $totob.reads_aligned_per_gene_per_million_raw_reads.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2b
  cat $toto2b | gawk '/#Run/{print;last;}' > $totov.999
  cat $totob.reads_aligned_per_gene_per_million_raw_reads.txt  | gawk '/^#/{next;}{print}' | scripts/tab_sort -k 8nr >> $totob.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/bacteria.metadata.txt ZZZZZ $totob.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Microbes | scripts/tab_sort -k 8,8nr >> $toto2b
  echo "\n\n" >> $toto2b

  cat $totot.reads_aligned_per_gene_per_million_raw_reads.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2b
  cat $toto2b | gawk '/#Run/{print;last;}' > $totov.999
  cat $totot.reads_aligned_per_gene_per_million_raw_reads.txt  | gawk '/^#/{next;}{print}' | scripts/tab_sort -k 8nr >> $totot.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/transposon.metadata.txt ZZZZZ $totot.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Transposon | scripts/tab_sort -k 8,8nr >> $toto2b
  echo "\n\n" >> $toto2b

  set toto2c=RESULTS/Microbiome/$MAGIC.virus_bacteria.length_aligned_per_run.txt
  echo -n "## $toto2c " > $toto2c
  cat $totov.gene_group.reads_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2c
  echo >> $toto2c
  cat $totov.gene_group.reads_aligned_per_gene.txt ZZZZZ $totov.gene_group.kb_aligned_per_gene.txt  | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){for(i=5;i<=NF;i++)u[$1,i]=$i;next;} printf("%s",$1);for(i=2;i<=NF;i++){z=u[$1,i]+0;if(z>0)printf("\t%d",int(1000*$i/z));else printf("\t%s",$i);}printf("\n");}' | scripts/tab_sort -k 8nr > $totov.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/virus.metadata.txt ZZZZZ $totov.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Virus | scripts/tab_sort -k 8,8nr >> $toto2c
  echo "\n\n" >> $toto2c
  cat $totov.reads_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2c
  echo >> $toto2c
  cat $totov.reads_aligned_per_gene.txt ZZZZZ $totov.kb_aligned_per_gene.txt  | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){for(i=5;i<=NF;i++)u[$1,i]=$i;next;} printf("%s",$1);for(i=2;i<=NF;i++){z=u[$1,i]+0;if(z>0)printf("\t%d",int(1000*$i/z));else printf("\t%s",$i);}printf("\n");}' | scripts/tab_sort -k 8nr > $totov.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/virus.metadata.txt ZZZZZ $totov.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Virus | scripts/tab_sort -k 8,8nr >> $toto2c
  echo "\n\n" >> $toto2c
  cat $totob.gene_group.reads_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2c
  cat $totob.gene_group.reads_aligned_per_gene.txt ZZZZZ $totob.gene_group.kb_aligned_per_gene.txt  | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){for(i=5;i<=NF;i++)u[$1,i]=$i;next;} printf("%s",$1);for(i=2;i<=NF;i++){z=u[$1,i]+0;if(z>0)printf("\t%d",int(1000*$i/z));else printf("\t%s",$i);}printf("\n");}' | scripts/tab_sort -k 8nr > $totob.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/bacteria.metadata.txt ZZZZZ $totob.999    | gawk  -F '\t'  -f scripts/vir2.k2.awk t=Microbes | scripts/tab_sort -k 8,8nr >> $toto2c
  echo "\n\n" >> $toto2c
  cat $totob.reads_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2c
  cat $totob.reads_aligned_per_gene.txt ZZZZZ $totob.kb_aligned_per_gene.txt  | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){for(i=5;i<=NF;i++)u[$1,i]=$i;next;} printf("%s",$1);for(i=2;i<=NF;i++){z=u[$1,i]+0;if(z>0)printf("\t%d",int(1000*$i/z));else printf("\t%s",$i);}printf("\n");}' | scripts/tab_sort -k 8nr > $totob.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/bacteria.metadata.txt ZZZZZ $totob.999    | gawk  -F '\t'  -f scripts/vir2.k2.awk t=Microbes | scripts/tab_sort -k 8,8nr >> $toto2c
  echo "\n\n" >> $toto2c
  cat $totot.reads_aligned_per_gene.txt  | gawk '/^#High/{next;}/^#/{print}'  >> $toto2c
  cat $totot.reads_aligned_per_gene.txt ZZZZZ $totot.kb_aligned_per_gene.txt  | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){for(i=5;i<=NF;i++)u[$1,i]=$i;next;} printf("%s",$1);for(i=2;i<=NF;i++){z=u[$1,i]+0;if(z>0)printf("\t%d",int(1000*$i/z));else printf("\t%s",$i);}printf("\n");}' | scripts/tab_sort -k 8nr > $totot.999
  cat MetaDB/$MAGIC/RunsList ZZZZZ ZZZZZ VirusDB/transposon.metadata.txt ZZZZZ $totot.999    | gawk  -F '\t' -f scripts/vir2.k2.awk t=Transposon | scripts/tab_sort -k 8,8nr >> $toto2c
  echo "\n\n" >> $toto2c



  scripts/vir2.k3.tcsh high 132 299 0
  scripts/vir2.k3.tcsh mid  110 132 0
  scripts/vir2.k3.tcsh low   60 110 0 

  scripts/vir2.k3.tcsh high 132 299 1
  scripts/vir2.k3.tcsh mid  110 132 1
  scripts/vir2.k3.tcsh low   60 110 1

 # scripts/vir2.k3.tcsh high 132 299 2
 # scripts/vir2.k3.tcsh mid  110 132 2
 # scripts/vir2.k3.tcsh low   60 110 2

# RESULTS/Microbiome/$MAGIC.virus_bacteria.kb_aligned_per_target.$type'_similarity.txt' 
# RESULTS/Microbiome/$MAGIC.virus_bacteria.coverage_per_target.$type'_similarity.txt' 
cat tmp/vir2.high.0 tmp/vir2.mid.0 tmp/vir2.low.0 > RESULTS/Microbiome/$MAGIC.virus_bacteria.kb_aligned_per_target.HighMidLowSimilarity.txt
cat tmp/vir2.high.1 tmp/vir2.mid.1 tmp/vir2.low.1 > RESULTS/Microbiome/$MAGIC.virus_bacteria.BasesPerHostMappedMb.HighMidLowSimilarity.txt
# NOT ADDITIVE does not work fro Sp(ecies)_group cat tmp/vir2.high.2 tmp/vir2.mid.2 tmp/vir2.low.2 > RESULTS/Microbiome/$MAGIC.virus_bacteria.coverage_per_target.HighMidLowSimilarity.txt

\rm tmp/vir2.high.[0-2] tmp/vir2.mid.[0-2] tmp/vir2.low.[0-2]
\rm tmp/vir2.high.[0-2].* tmp/vir2.mid.[0-2].* tmp/vir2.low.[0-2].*

######## create a .ace file : Ali->High_genes Virus/bacteria counting the 5 bests
  set toto4=tmp/VIRUS/$MAGIC.virus_bacteria.high_genes.ace
  cat $toto2a | gawk -F '\t' '/^#Run/{for(i=2;i<=NF;i++) {r[i]=$i;if($i=="#Run")iMin=i+2;}next;}{if($6=="Virus"   ){for(i=iMin;i<=NF;i++)printf("%s\t%s\t%s\t%s\n",r[i],$1,$i,$5);}}' | sort -k 1,1 -k 3,3nr > toto
  cat toto |  gawk -F '\t' '{if(length($1)<2)next;}{if($1!=r){n=0;r=$1;printf("\nAli %s\n-D High_genes virus\n",$1);}if(n>5)next;n++;printf("High_genes virus %s %s \"%s\"\n",$2,$3,$4);}END{printf("\n");}' > $toto4
  cat $toto2a | gawk -F '\t' '/^#Run/{for(i=2;i<=NF;i++) {r[i]=$i;if($i=="#Run")iMin=i+2;}next;}/^#/{next;}{if($6=="Microbe"){for(i=iMin;i<=NF;i++)printf("%s\t%s\t%s\t%s\n",r[i],$1,$i,$5);}}' | sort -k 1,1 -k 3,3nr > toto
  cat toto |  gawk -F '\t' '{if(length($1)<2)next;}{if($1!=r){n=0;r=$1;printf("\nAli %s\n-D High_genes microbe\n-D High_genes bacteria\n",$1);}if(n>5)next;n++;printf("High_genes microbe %s %s \"%s\"\n",$2,$3,$4);}END{printf("\n");}' >> $toto4
  cat $toto2a | gawk -F '\t' '/^#Run/{for(i=2;i<=NF;i++) {r[i]=$i;if($i=="#Run")iMin=i+2;}next;}/^#/{next;}{if($6=="Transposon"){for(i=iMin;i<=NF;i++)printf("%s\t%s\t%s\t%s\n",r[i],$1,$i,$5);}}' | sort -k 1,1 -k 3,3nr > toto
  cat toto |  gawk -F '\t' '{if(length($1)<2)next;}{if($1!=r){n=0;r=$1;printf("\nAli %s\n-D High_genes transposon\n",$1);}if(n>5)next;n++;printf("High_genes transposon %s %s \"%s\"\n",$2,$3,$4);}END{printf("\n");}' >> $toto4

  set toto=tmp/BACTERIA/$MAGIC.counts.ace
  cat $toto | gawk '/^Gene/{g=$2;next;}/^Run_U/{r=$2;kb=$8;if(kb>1)printf("%s\t%s\t%s\n",r,g,kb);}'| sort -k 1,1 -k3,3nr | gawk -F '\t' '{r=$1;g=$2;kb=$3;if(r!=old)n=0;n++;if(n<5)printf("Ali %s\nHigh_genes BACTERIA %s %.3f Mb\n\n",r,g,kb/1000.0);}' >> $toto4


# collate the counts
  foreach run (`cat MetaDB/$MAGIC/RunList `)
    echo "Ali $run"  >> $toto4
    cat tmp/VIRUS/$run/vir1.virus.HitMult tmp/BACTERIA/$run/vir1.bacteria.HitMult tmp/BACTERIA/$run/vir1.transposon.HitMult |  gawk -F '\t' '{if(substr($2,1,3)=="any")next;}/HITS/{z= $2; if (t[z]<1){t[z]=1;nt++;t2z[nt]=z;z2t[z]=nt;}it=z2t[z];if(imax<NF)imax=NF;for (i=4;i<=NF;i++)nn[it,i]+=$i;}END{for(it=1;it<=nt;it++){printf("HITS\t%s.new",t2z[it]);for (i=4;i<=imax;i++)printf("\t%d",nn[it,i]);printf("\n");}}'  > $toto4.hits
    cat $toto4.hits | gawk -F '\t' '{z=$1; if($3<1)next;printf ("nh_Ali %s %.2f seq %.2f tags  %.3f kb_ali  %.2f bp_av_ali %.3f kb_clip %.2f bp_av_clip\n",$2, $3,$7,$8/1000,$8/$7,$9/1000,$9/$7);d=1;if ($14>0)d=$14; printf ("h_Ali %s %d seq %d tags %d kb_ali %.2f bp_av_ali %d kb_clip %.2f bp_av_clip\n", $2, $10,$14,$15/1000,$15/d,$16/1000,$16/d) ;  n=$4+$5; namb = $6 ; }' >> $toto4
    cat tmp/VIRUS/$run/vir1.virus.HitMult tmp/BACTERIA/$run/vir1.bacteria.HitMult tmp/BACTERIA/$run/vir1.transposon.HitMult | gawk -F '\t' '/any/{next;}/MULT/{z= $1  "\t" $2 ; if (t[z]<1){t[z]=1;nt++;t2z[nt]=z;z2t[z]=nt;}it=z2t[z];if(imax<NF)imax=NF;for (i=3;i<=NF;i++)nn[it,i]+=$i;}END{for(it=1;it<=nt;it++){printf("\n%s",t2z[it]);for (i=3;i<=imax;i++)printf("\t%d",nn[it,i]);}printf("\n");}' | sort > $toto4.mult 
    cat $toto4.mult |  gawk -F '\t' '{z = $1 ; if($3<10)next;printf ("Unicity %s",$2) ; for(i=3;i<=NF;i++)printf(" %s",$i);printf("\n");}' >> $toto4
    echo >> $toto4
  end

  echo " pparse  $toto4" | bin/tacembly MetaDB -noprompt

goto phaseLoop
exit 0

cat <<EOF > toto
NC_001806.2 gi|820945227|ref|NC_001806.2| Human herpesvirus 1 strain 17, complete genome
NC_001798.2 gi|820945149|ref|NC_001798.2| Human herpesvirus 2 strain HG52, complete genome
NC_001348.1 gi|9625875|ref|NC_001348.1| Human herpesvirus 3, complete genome
NC_009334.1 gi|139424470|ref|NC_009334.1| Human herpesvirus 4, complete genome
NC_007605.1 gi|82503188|ref|NC_007605.1| Human herpesvirus 4 complete wild type genome
NC_006273.2 gi|155573622|ref|NC_006273.2| Human herpesvirus 5 strain Merlin, complete genome
NC_001664.2 gi|224020395|ref|NC_001664.2| Human herpesvirus 6A, complete genome
NC_000898.1 gi|9633069|ref|NC_000898.1| Human herpesvirus 6B, complete genome
NC_001716.2 gi|51874225|ref|NC_001716.2| Human herpesvirus 7, complete genome
NC_009333.1 gi|139472801|ref|NC_009333.1| Human herpesvirus 8, complete genome
EOF
cat toto | gawk '{print $1}' >  tmp/VIRUS/herpes.list

# export the intron support
if (-e tmp/VIRUS/$MAGIC.intron.pre) \rm tmp/VIRUS/$MAGIC.intron.pre
foreach run (`cat MetaDB/$MAGIC/RunList`)
  if (-d tmp/PHITS_virus/$run ) then
    gunzip -c tmp/PHITS_virus/$run/*.introns.gz | gawk -F '\t' '/^#/{next}{n[$3 "\t" $5 "\t" $7 "\t" $9]+=$11;}END{for(k in n) printf ("%s\t%d\n",k,n[k])}' | sort -k 1 -k 2n | gawk -F '\t' '{if($2<$3){dx=$3-$2+1;s="+";}else{dx=$2-$3+1;s="-";}if($5>0)printf("%s\t%d\t%d\t%s\t%d\t%s\t%d\n",$1,$2+0,$3+0,s,dx,$4,$5);}' >> tmp/VIRUS/$MAGIC.intron.pre
  endif
end
cat  tmp/VIRUS/$MAGIC.intron.pre | gawk -F '\t' '{n[$1 "\t" $2  "\t" $3 "\t" $4 "\t" $5 "\t" $6] += $7;}END{for (k in n)printf("%s\t%d\n",k,n[k]);}' | sort -k 1,1 -k 7nr > tmp/VIRUS/$MAGIC.intron.pre2
 
if (-e tmp/VIRUS/$MAGIC.intron.preB) \rm tmp/VIRUS/$MAGIC.intron.preB
foreach run (`cat MetaDB/$MAGIC/RunList`)
  if (-d tmp/PHITS_bacteria/$run ) then
    gunzip -c tmp/PHITS_bacteria/$run/*.introns.gz | gawk -F '\t' '/^#/{next}{n[$3 "\t" $5 "\t" $7 "\t" $9]+=$11;}END{for(k in n) printf ("%s\t%d\n",k,n[k])}' | sort -k 1 -k 2n | gawk -F '\t' '{if($2<$3){dx=$3-$2+1;s="+";}else{dx=$2-$3+1;s="-";}if($5>0)printf("%s\t%d\t%d\t%s\t%d\t%s\t%d\n",$1,$2+0,$3+0,s,dx,$4,$5);}' >> tmp/VIRUS/$MAGIC.intron.preB
  endif
end
cat  tmp/VIRUS/$MAGIC.intron.preB | gawk -F '\t' '{n[$1 "\t" $2  "\t" $3 "\t" $4 "\t" $5 "\t" $6] += $7;}END{for (k in n)printf("%s\t%d\n",k,n[k]);}' | sort -k 1,1 -k 7nr > tmp/VIRUS/$MAGIC.intron.pre2B
 
echo -n "## Support for introns or small deletions in viruses in project $MAGIC : " > RESULTS/Microbiome/$MAGIC.virus_intron_support.txt
date >> RESULTS/Microbiome/$MAGIC.virus_intron_support.txt
echo "## Introns are discovered by mapping discontinuously using Magic on a number of pre-selected viruses" >> RESULTS/Microbiome/$MAGIC.virus_intron_support.txt
cat  tmp/VIRUS/virus_names.txt ZZZZZ  tmp/VIRUS/$MAGIC.intron.pre2 |  gawk  -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){split($2,aa," ");ok[aa[1]]=$2;next}}BEGIN{printf("# Accession\tTitle\tCoordinate of intron first base on accession\tCoordinate of intron last base\tIntron boundary\tIntron length (nt)\tStrand\tNumber of reads supporting the intron\n");}{if($1 != old)printf("\n");old=$1;printf("%s\t%s\t%d\t%d\t%s\t%d\t%s\t%d\n",$1,ok[$1],$2,$3,$6,$5,$4,$7);}' >> RESULTS/Microbiome/$MAGIC.virus_intron_support.txt


# extract the herpesvirus names and export their global coverage
cat  tmp/VIRUS/virus_names.txt | grep herpes | cut -f 1 >  tmp/VIRUS/herpes.list
\cp  tmp/VIRUS/herpes.list  tmp/VIRUS/star.list
# corona wuhan # herpes 7
echo NC_045512 > tmp/VIRUS/star.list
echo NC_001716 >> tmp/VIRUS/star.list    

foreach v (`cat tmp/VIRUS/star.list`)
  echo ".... virus.tcsh: cumulating the wiggles : $v"
  echo >  tmp/VIRUS/_h
  echo ' ' > tmp/VIRUS/$MAGIC.$v.hits
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    zcat tmp/VIRUS/$run/*.v_virus.hits.gz | gawk -F '\t' '{if ($11 == v) print;}' v=$v >>  tmp/VIRUS/$run/$v.hits
    cat tmp/VIRUS/$run/$v.hits >> tmp/VIRUS/$MAGIC.$v.hits
  end
  cat tmp/VIRUS/$MAGIC.$v.hits | wiggle -I BHIT -O BV -o tmp/VIRUS/$MAGIC.$v -force_unique 
end

set v=NC_045512
cat tmp/VIRUS/MM.$v.hits  | gawk -F '\t' '{if(index($1,">")>0)ir=1;else ir=-1;da=ir*($13-$12);if(da<0)print;}' > tmp/VIRUS/MM.$v.f.hits
cat tmp/VIRUS/MM.$v.hits  | gawk -F '\t' '{if(index($1,">")>0)ir=1;else ir=-1;da=ir*($13-$12);if(da>0)print;}' > tmp/VIRUS/MM.$v.r.hits


cat tmp/VIRUS/MM.$v.f.hits | wiggle -I BHIT -O BV -o tmp/VIRUS/$MAGIC.$v.f -force_unique
cat tmp/VIRUS/MM.$v.r.hits | wiggle -I BHIT -O BV -o tmp/VIRUS/$MAGIC.$v.r -force_unique
cat tmp/VIRUS/MM.$v.f.BV | wiggle -I BV -O count -cumul
cat tmp/VIRUS/MM.$v.r.BV | wiggle -I BV -O count -cumul


foreach v (`cat tmp/VIRUS/star.list`)
  echo ".... virus.tcsh: detecting the SNPs : $v"
  echo >  tmp/VIRUS/_h
  echo ' ' > tmp/VIRUS/$MAGIC.$v.hits
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    bin/snp -i tmp/VIRUS/$run/$v.hits -minAliPerCent 60 --fasta tmp/VIRUS/$v.fasta  --hits2BRS --run $run -o tmp/VIRUS/$run/$v.BRS
    bin/snp -i tmp/VIRUS/$run/$v.BRS -db MetaDB -minFrequency 1 -minCover 1000 -BRS_detect > tmp/VIRUS/$run/$v.detect
  end
end

set v=NC_045512
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    bin/snp -i tmp/VIRUS/$run/$v.hits -minAliPerCent 60 --fasta tmp/VIRUS/$v.fasta  --hits2BRS --run $run -o tmp/VIRUS/$run/$v.BRS
  end
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    bin/snp -i tmp/VIRUS/$run/$v.BRS -db MetaDB -minFrequency 10 -minCover 1000 -BRS_detect > tmp/VIRUS/$run/$v.detect
  end
cat tmp/VIRUS/*/$v.detect | gawk -F '\t' '/^#/{next;}{print $1;}' | sort -u > tmp/VIRUS/$v.snp_list
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    bin/snp -i tmp/VIRUS/$run/$v.BRS -db MetaDB -minFrequency 10 -minCover 1000 -BRS_count -snp_list tmp/VIRUS/$v.snp_list > tmp/VIRUS/$run/$v.count 
  end
cat tmp/VIRUS/*/$v.count | gawk -F '\t' '{if($9 >= m)print}' m=1000 | sort -k 1,1 -k 2,2n -k 3,3 >>  tmp/VIRUS/$v.count.sorted



# analyse the SNPs

# creta a special virome TARGET, then run
# MAGIC a0C wait a0D wait a123 wait c1 wait c2 wait c3 wait c4 wait c5 wait c6 wait c7








 set toto=RESULTS/Microbiome/$MAGIC.virus_wiggles.BF
 echo -n '## ' > $toto
 date >> $toto
 echo "## Coverage plots of selected viruses in all samples of project $MAGIC" >> $toto
 cat $MAGIC.toto  | scripts/transpose | gawk -F '\t' '{ln++; if (ln<3 || ln == 4)next;if(ln == 3) {printf ("Coordinate");for(i=1;i<=NF;i++){split($i,aa," ") ;printf("\t%s",substr(aa[2],7));}printf ("\n");next;}printf("%d\t", 10 * (ln-4)) ; print}' >> $toto

 cat $toto | head -2 > $toto.1
 cat $toto | head -3 | tail -1 | transpose > $toto.2
 cat tmp/VIRUS/virus_names.txt ZZZZZ $toto.2 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$2]=1;next}}{k0=$1;for(k in ok){if(index(k,$1)>0)k0=k;}print k0}' | transpose  >> $toto.1
 cat $toto | tail -n +4 >> $toto.1
 mv $toto.1 $toto
 \rm $toto.2


\rm  $MAGIC.toto.Rh $MAGIC.toto



# extract the herpesvirus names and export their global coverage

foreach v (NZ_CP007181.1 AL111168 HE978252.1 CP000814.1)
  echo ".... virus.tcsh: cumulating the wiggles : $v"
  echo >  tmp/VIRUS/_h
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    cat tmp/VIRUS/$run/*.vir1.hits | gawk -F '\t' '{split ($11,aa,"|");if($11 == v || aa[4] == v) print;}' v=$v >>  tmp/VIRUS/_h
  end
  set n=`wc -l  tmp/VIRUS/_h | gawk '{print $1}' `
  if ($n > 1) then
    wiggle -i  tmp/VIRUS/_h -I BHIT -O BV -o tmp/VIRUS/$MAGIC.$v -force_unique 
  endif
end

set run=COV-20200314-P9-B05-P 
set lane=COV-20200314-P9-B05-P/f2.2 
bin/variant_caller -target_class Z_genome -fastcFile Fastc/$lane.fastc.gz -hitFile tmp/COUNT/$lane.hits.gz -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 1000 -minSnpCount 100 -dx 8 
bin/variant_caller -target_class Z_genome -run $run -target_fasta TARGET/Targets/corona.genome.fasta.gz -t NC_045512 -method VIRUS -minSnpCover 1000 -minSnpCount 100 -dx 12 -o tatou12



######
## extract the profile of the homozygous snps of each run

setenv MAGIC MMM
set zone=zoneG.NC_045512_a
set target=NC_045512
set ff=tmp/SNPH/$zone/$MAGIC.snp.sorted.homozygous

cat $ff | gawk -F '\t' /^$target/'{t=$1;r=$6;runs[r]=1;if($8>90)w[r]=w[r]","$2":"$3;}END{for (r in runs)printf("%s\t%s\t%s\n",t,r,substr(w[r],2));}' | head -5

attaaaggtttataccttcccaggtaacaaaccaaccaactttcgatctcttgtagatctgttctctaaacgaactttaaaatctgtgtggctgtcactcggctgcatgcttagtgcactcacgcagtataattaataactaattactgtcgttgacaggacacgagtaactcgtctatcttctgcaggctgcttacggtttcgtccgtgttgcagccgatcatcagcacatctaggtttcgtccgggtgt
          tataccttcccaggtaacaaaccaaccaactttcgatctcttgtagatctgttctctaaacgaactttaaaatctgtgtggctgtcactcggctgcatgcttagtgcactcacgcagtataattaataactaattactgtcgttgacaggacacgagtaactcgtctatcttctgcaggctgcttacggtttcgtccgtgttgcagccgatcatcagcacatctaggttttgtccgggtgtgaccgaaagg


# export all wiggles as a single file

set target=NC_045512
set type=u.f
set toto=tmp/WIGGLERUN/$MAGIC.all_wiggles_unique_forward.txt
echo -n "### $toto : " > $toto
echo " " > $toto.1
date >> $toto
echo "### Coverage plots on the forward strand by uniquely aligned pairs" >> $toto


set type=pp.f
set toto=tmp/WIGGLERUN/$MAGIC.all_wiggles_partial_forward.txt
echo -n "### $toto : " > $toto
echo " " > $toto.1
date >> $toto
echo "### Coverage plots on the forward strand by uniquely partially aligned pairs" >> $toto

foreach run (`cat MetaDB/MMM/GroupListSorted`)
  set f1=tmp/WIGGLEGROUP/$run/$target/R.chrom.$type.BF.gz
  if (! -e $f1) continue
  echo $run
  wiggle -I BF -O BV -i $f1 | gawk '{if($1+0 >0) printf("%s\t%d\t%d\n",run,$1,$2);}' run=$run >> $toto.1
end
foreach run (`cat MetaDB/MMM/RunListSorted`)
  set f1=tmp/WIGGLERUN/$run/$target/R.chrom.$type.BF.gz
  if (! -e $f1) continue
  wiggle -I BF -O BV -i $f1 | gawk '{if($1+0 >0) printf("%s\t%d\t%d\n",run,$1,$2);}' run=$run >> $toto.1
end
cat $toto.1 | gawk -F '\t' '{r=$1;i=0+r2i[r];if(i==0){iMax++;i=iMax;r2i[r]=i;i2r[i]=r;}x=$2;z[i,x]=$3;}END {printf("# Run");for (i=1;i<=iMax;i++)printf("\t%s",i2r[i]);for(x=0;x<=30000;x+=10){printf("\n%d",x);for(i=1;i<=iMax;i++)printf("\t%d",z[i,x]);}printf("\n");}' >> $toto
\cp $toto RESULTS/Coverage_and_exons


phaseLoop:
 echo done


