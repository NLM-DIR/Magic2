#!bin/tcsh -f

set phase=$1
set run=$2
setenv target $3

if ($phase == t1) goto phaseT1
if ($phase == t2) goto phaseT2
if ($phase == t3) goto phaseT3
if ($phase == t4) goto phaseT4
if ($phase == OTHERS) goto phaseOTHERS
exit (0)

##### Phase 1 collect the gene_fusions
phaseT1:

  if ($target == magic) continue
  # collate the gene fusions
  source scripts/target2target_class.txt
        
  set spm=""
  if (-e tmp/METADATA/av.split_mrnas.gz) set spm="-splitMrnas tmp/METADATA/av.split_mrnas.gz"
  set ma="-minAli 100"
  set ma="-minAli 30"

  if (! -e  tmp/GeneFusion/$run/t1.gene_fusion.$target.txt.gz1 && -e tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt) then
    echo ' ' >  tmp/GeneFusion/$run/t1.gene_fusion.$target.txt
    if (-e tmp/GeneFusion/$run/t1.gene_fusion.LaneList) \rm tmp/GeneFusion/$run/t1.gene_fusion.LaneList
    foreach lib ($run `cat MetaDB/$MAGIC/r2sublib | gawk '{if(r==$1)print " "$2;}' r=$run`)
        mkdir  tmp/GeneFusion/$run/$lib
        if (-e Fastc/$lib/LaneList) cat Fastc/$lib/LaneList >> tmp/GeneFusion/$run/t1.gene_fusion.LaneList
    end
    foreach lane (`cat tmp/GeneFusion/$run/t1.gene_fusion.LaneList`)
        echo "t1 $lane"
        if (-e tmp/COUNT/$lane.hits.gz) then
          echo "bin/tricoteur  -hitFile tmp/COUNT/$lane.hits.gz -run $run -lane $lane -target_class  $target_class -geneFusion -o tmp/GeneFusion/$run/$lane -geneMap tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt -mrnaRemap tmp/METADATA/mrnaRemap.gz $spm $ma"
                bin/tricoteur  -hitFile tmp/COUNT/$lane.hits.gz -run $run -lane $lane -target_class  $target_class -geneFusion -o tmp/GeneFusion/$run/$lane -geneMap tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt -mrnaRemap tmp/METADATA/mrnaRemap.gz $spm $ma
          sleep 1
          if (-e tmp/GeneFusion/$run/$lane.geneFusion.txt) then
            ls -ls tmp/GeneFusion/$run/$lane.geneFusion.txt
            cat tmp/GeneFusion/$run/$lane.geneFusion.txt >>  tmp/GeneFusion/$run/t1.gene_fusion.$target.txt
            \rm  tmp/GeneFusion/$run/$lane.geneFusion.txt
          endif
        endif
    end
    gzip tmp/GeneFusion/$run/t1.gene_fusion.$target.txt
  endif
# count
  zcat  tmp/GeneFusion/$run/t1.gene_fusion.$target.txt.gz | gawk -F '\t' '/^#/{next;}{z = $3 ;} /PAIR/{nn[z]++;np[z]++;next;}/READ/{nn[z]++;nr[z]++;}END{for(k in nn)if(nn[k]>0)printf("%s\t%d\t%d\n",k,nr[k],np[k]);}' | sort -k 2nr > tmp/GeneFusion/$run/t1.gene_fusion.$target.count
  if (-e tmp/GeneFusion/$run/t1.gene_fusion.$target.txt.gz) then
    ls -ls tmp/GeneFusion/$run/t1.gene_fusion.$target.*
    echo "bin/tricoteur -run $run -t $target -geneFusionFile tmp/GeneFusion/$run/t1.gene_fusion.$target.txt.gz -o   tmp/GeneFusion/$run/t1.gene_fusion.$target"
          bin/tricoteur -run $run -t $target -geneFusionFile tmp/GeneFusion/$run/t1.gene_fusion.$target.txt.gz -o   tmp/GeneFusion/$run/t1.gene_fusion.$target
    ls -ls tmp/GeneFusion/$run/t1.gene_fusion.$target.*
  endif

  
  touch tmp/GeneFusion/$run/t1.gene_fusion.$target.done

goto phaseLoop

zcat tmp/GeneFusion/PacBio.ROCR3.CL2-F1*/t1.gene_fusion.av.txt.gz | grep BCAS4__BCAS3 | grep f2.26 | cut -f 2 | sort -u > toto.list
dna2dna -i Fastc/PacBio.ROCR3.CL2-F1-2_cell9/f2.26.fastc.gz -select toto.list -I fastc -O fastc -o toto
\rm toto.hits
foreach r (`cat toto.list`)
  zcat tmp/COUNT/PacBio.ROCR3.CL2-F1-2_cell9/f2.26.hits.gz | grep $r | grep ET_av >> toto.hits
end

cat toto.hits | grep n.2346430#1 > toto1.hits

bin/tricoteur -hitFile toto1.hits -run PacBio.ROCR3.CL2-F1-2_cell9 -lane PacBio.ROCR3.CL2-F1-2_cell9/f2.26 -target_class ET_av -geneFusion  -geneMap tmp/METADATA/av.mrna_map_ln_gc_gene_geneid.txt -mrnaRemap tmp/METADATA/mrnaRemap.gz 

#######################################################################################
## phase t2:   cumul the gene_fusions per group
phaseT2:

goto phaseLoop

## Create a project table au gene pres
phaseT3:

set minFusion=50
set target=av
if (! -e tmp/GeneFusion/t3.$MAGIC.a) then
  echo " " >  tmp/GeneFusion/t3.$MAGIC.a
  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    set ff=tmp/GeneFusion/$run/t1.gene_fusion.$target.txt.gz
    if (! -e $ff) continue
#   zcat $ff  | grep BCAS4__BCAS3 >> tmp/GeneFusion/t3.$MAGIC.a
    zcat $ff  >> tmp/GeneFusion/t3.$MAGIC.a
  end
  cat tmp/GeneFusion/t3.$MAGIC.a | gawk -F '\t' '/^#/{next;}{f=$3;nnDS[f]+=$20;nnS[f]+=$21;nnG1[f]+=$22;nnG2[f]+=$23;next;}END{for(f in nnS) if(nnS[f]>=min)printf("%s\t%d\t%d\t%d\t%d\n",f,nnDS[f],nnS[f],nnG1[f],nnG2[f]);}' min=$minFusion | sort -k 2nr > tmp/GeneFusion/t3.$MAGIC.b
endif


if (! -e tmp/GeneFusion/t3.$MAGIC.c1) then
  echo $target > tmp/GeneFusion/t3.$MAGIC.c1
  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    set ff=tmp/GeneFusion/$run/t1.gene_fusion.$target.count
    if (! -e $ff) continue
    cat tmp/GeneFusion/t3.$MAGIC.b ZZZZZ $ff | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz+0<1){ok[$1]=1;next;}}{if(ok[$1]==1){printf("%s\t",run);print;}}' run=$run >> tmp/GeneFusion/t3.$MAGIC.c1
  end
endif

if (! -e tmp/GeneFusion/t3.$MAGIC.c2) then
  echo $target > tmp/GeneFusion/t3.$MAGIC.c2
  set groupLevelMax=`cat  MetaDB/$MAGIC/g2r | cut -f 3 | sort -k 1n | tail -1`
  foreach level (`seq 1 1 $groupLevelMax`)
    echo level=$level
    set okk=0
    foreach group (`cat MetaDB/$MAGIC/GroupListSorted`)
      set gLevel=`cat MetaDB/$MAGIC/g2r |  gawk -F '\t' '{if($1==g)level=$3}END{print 0+level}' g=$group` 
      if ($gLevel != $level) continue
      if (! -d tmp/GeneFusion/$group) mkdir tmp/GeneFusion/$group
      set ff=tmp/GeneFusion/$group/t1.gene_fusion.$target.count
      ls -ls $ff
      echo " " >  $ff.1
      foreach run2 (`cat MetaDB/$MAGIC/g2r | gawk '{if($1==group) print $2;}' group=$group| sort -u`)
        set ff2=tmp/GeneFusion/$run2/t1.gene_fusion.$target.count
        if (! -e $ff) continue
        cat $ff2 >> $ff.1
      end
      cat $ff.1 | gawk -F '\t' '{g=$1;n1[g]+=$2;n2[g]+=$3;}END{for (g in n1)printf("%s\t%d\t%d\n",g,n1[g],n2[g]);}' > $ff
      cat tmp/GeneFusion/t3.$MAGIC.b ZZZZZ $ff | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz+0<1){ok[$1]=1;next;}}{if(ok[$1]==1){printf("%s\t",group);print;}}' group=$group >> tmp/GeneFusion/t3.$MAGIC.c2
    end
  end
endif

phaseT3a:
  set target=av
  set toto=tmp/GeneFusion/t3.$MAGIC.$target.txt
  cat  MetaDB/$MAGIC/GroupListSorted2 ZZZZZ MetaDB/$MAGIC/RunListSorted2  ZZZZZ MetaDB/$MAGIC/Run2Title.txt ZZZZZ tmp/GeneFusion/t3.$MAGIC.c[12]    | gawk -F '\t' '/^ZZZZZ/{zz++;if(zz==1)lastG=iMax;next;}{gsub(/\"/,"",$0);r=$1;i=r2i[r]+0;if(i==0){iMax++;i=iMax;r2i[r]=i;i2r[i]=r;title[r]=r;}  if(zz+0==2)title[r]=$3;if(zz+0<=2)next;k=$2;if(i>lastG){nn[k]+=$3+$4;nn1[k]+=$3;nn2[k]+=$4;nnn1+=$3;nnn2+=$4;}nni[i]+=$3+$4;n1[k,i]+=$3;n2[k,i]+=$4;}END{printf("# Run\tAny\tSupporting reads");for(i=1;i<=iMax;i++)if (nni[i]>0)printf("\t%s",i2r[i]);if(nnn2>0){printf("\t\t# Run\tAny\tAny2");for(i=1;i<=iMax;i++)if (nni[i]>0)printf("\t%s",i2r[i]);}printf("\n# Title\tAny\tAny1");for(i=1;i<=iMax;i++)if (nni[i]>0){r=i2r[i];t=title[r];if(length(t)==0)t=r;printf("\t%s",t);}if(nnn2>0){printf("\t\t# Title\tAny\tSupporting read-pairs");for(i=1;i<=iMax;i++)if (nni[i]>0){r=i2r[i];t=title[r];if(length(t)==0)t=r;printf("\t%s",t);}}for (k in nn){ printf("\n%s\t%d\t%d",k,nn[k],nn1[k]);for(i=1;i<=iMax;i++)if (nni[i]>0)printf("\t%d",n1[k,i]);if(nnn2>0){printf("\t\t%s\t%d\t%d",k,nn[k],nn2[k],nn2[k]);for(i=1;i<=iMax;i++)if (nni[i]>0)printf("\t%d",n2[k,i]);}}}END{printf("\n");}'  > tmp/GeneFusion/t3.$MAGIC.d
  echo -n "### File $toto : " > $toto
  date >> $toto
  echo "### Table of candidate gene fusions in project $MAGIC. Left table: single read  support, right table: pair support" >> $toto
  cat tmp/GeneFusion/t3.$MAGIC.b | gawk '/^#/{printf("#\t\t");print}' >> $toto

  if (! -e tmp/GeneFusion/$target.gene2intmap.txt) then
    cat tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt | gawk -F '\t' '{gene=$5;split($2,aa,":");chrom[gene]=aa[1];split(aa[2],bb,"-");a1=bb[1];a2=bb[2];if(a1>a2){a0=a1;a1=a2;a2=a0;}if(0+aa1[gene]==0){aa1[gene]=a1;aa2[gene]=a2;}if(a1<aa1[gene])aa1[gene]=a1;if(a2>aa2[gene])aa2[gene]=a2;}END{for (g in chrom)printf("%s\t%s\t%d\t%d\n",g,chrom[g],aa1[g],aa2[g]);}' | sort > tmp/GeneFusion/$target.gene2intmap.txt
  endif
  cat  tmp/GeneFusion/t3.$MAGIC.d | gawk -F '\t' '/^#/{printf("#\t\t");print;}' >> $toto
  cat  tmp/GeneFusion/$target.gene2intmap.txt ZZZZZ tmp/GeneFusion/t3.$MAGIC.d | gawk -F '\t' '/^ZZZZZ/{zz++;next;}/^#/{next;}{if(zz+0<1){g=$1;gc[g]=$2;g1[g]=$3;g2[g]=$4;next;}}{split($1,aa,"__");ga=aa[1];gb=substr(aa[2],1,length(aa[2])-2);printf("%s(%s:%d-%d__%s:%d-%d)",$1,gc[ga],g1[ga],g2[ga],gc[gb],g1[gb],g2[gb]);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' | sort -k 2nr > tmp/GeneFusion/t3.$MAGIC.e
  cat tmp/METADATA/$MAGIC.$target.captured_genes.g2c  ZZZZZ tmp/GeneFusion/t3.$MAGIC.e | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){cc[$1]=$2;next;}}{split($1,aa,"__");g1=aa[1];split(aa[2],bb,"(");k=length(bb[1]);g2=substr(bb[1],1,k-2);printf("%s\t%s\t",cc[g1],cc[g2]);print;}' >> $toto

  \cp $toto RESULTS/GeneFusion


goto phaseLoop




## Create a project table a la base pres: 
phaseT4:

set minFusion=50
set target=av
  echo " " >  tmp/GeneFusion/t3.$MAGIC.a
  foreach run (`cat MetaDB/$MAGIC/RunsList `)
    set ff=tmp/GeneFusion/$run/t1.gene_fusion.av.txt.gz 
    if (! -e $ff) continue
    zcat  $ff | gawk -F '\t' '/^#/{next;}/++/{if(($10 < $11 && $11 >= $15-1 && $11 < $15+12)||($10 > $11 && $11 <= $15+1 && $11 > $15-12))print ;}'  >> tmp/GeneFusion/t4.$MAGIC.a
  end
  cat tmp/GeneFusion/t4.$MAGIC.a | gawk -F '\t' '/^#/{next;}{f=$3"\t"$12"\t"$14"\t"$17"\t"$18;nn[f]++;next;}END{for(f in nn) if(nn[f]>=min)printf("%s\t%d\n",f,nn[f]);}' min=$minFusion | sort -k 6nr > tmp/GeneFusion/t4.$MAGIC.b

#  cat tmp/GeneFusion/t3.$MAGIC.a | gawk -F '\t' '/^#/{next;}{f=$1"\t"$12"\t"$17;nn[f]+=$2+$3;next;}END{for(f in nn) if(nn[f]>=min)printf("%s\t%d\n",f,nn[f]);}' min=$minFusion | sort -k 2nr > tmp/GeneFusion/t3.$MAGIC.b

  echo $target > tmp/GeneFusion/t4.$MAGIC.c
  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    set ff=tmp/GeneFusion/$run/t1.gene_fusion.av.txt.gz 
    if (! -e $ff) continue
    zcat tmp/GeneFusion/t4.$MAGIC.b ZZZZZ $ff | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz+0<1){ok[$1]=1;next;}}{if(ok[$1]==1){printf("%s\t",run);print;}}' run=$run >> tmp/GeneFusion/t4.$MAGIC.c
  end
  foreach run (`cat MetaDB/$MAGIC/GroupListSorted`)
    if (! -d tmp/GeneFusion/$run) mkdir tmp/GeneFusion/$run
    set ff=tmp/GeneFusion/$run/t1.gene_fusion.av.txt.gz 
    echo " " >  tmp/GeneFusion/$run/t1.gene_fusion.$target.count
    foreach run2 (`cat MetaDB/$MAGIC/g2r | gawk '{if($1==group) print $2;}' group=$run | sort -u`)
      set ff2=tmp/GeneFusion/$run2/t1.gene_fusion.av.txt.gz 
      if (! -e $ff) continue
      zcat $ff2 >> $ff
    end
    cat tmp/GeneFusion/t4.$MAGIC.b ZZZZZ $ff | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz+0<1){ok[$1]=1;next;}}{if(ok[$1]==1)n1[$1]+=$2;n2[$1]+=$3;}END{for (k in n1) {printf("%s\t%s\t%d\t%d\n",run,k,n1[k],n2[k]);}}'  run=$run  >> tmp/GeneFusion/t4.$MAGIC.c
  end

phaseT4a:
  set target=av
  set toto=tmp/GeneFusion/t4.$MAGIC.$target.txt
  cat  MetaDB/$MAGIC/GroupListSorted2 ZZZZZ MetaDB/$MAGIC/RunListSorted2  ZZZZZ MetaDB/$MAGIC/Run2Title.txt ZZZZZ tmp/GeneFusion/t4.$MAGIC.c    | gawk -F '\t' '/^ZZZZZ/{zz++;if(zz==1)lastG=iMax;next;}{gsub(/\"/,"",$0);r=$1;i=r2i[r]+0;if(i==0){iMax++;i=iMax;r2i[r]=i;i2r[i]=r;title[r]=r;}  if(zz+0==2)title[r]=$3;if(zz+0<=2)next;k=$2;if(i>lastG){nn[k]+=$3+$4;nn1[k]+=$3;nn2[k]+=$4;nnn1+=$3;nnn2+=$4;}nni[i]+=$3+$4;n1[k,i]+=$3;n2[k,i]+=$4;}END{printf("# Run\tAny\tSupporting reads");for(i=1;i<=iMax;i++)if (nni[i]>0)printf("\t%s",i2r[i]);if(nnn2>0){printf("\t\t# Run\tAny\tAny2");for(i=1;i<=iMax;i++)if (nni[i]>0)printf("\t%s",i2r[i]);}printf("\n# Title\tAny\tAny1");for(i=1;i<=iMax;i++)if (nni[i]>0){r=i2r[i];t=title[r];if(length(t)==0)t=r;printf("\t%s",t);}if(nnn2>0){printf("\t\t# Title\tAny\tSupporting read-pairs");for(i=1;i<=iMax;i++)if (nni[i]>0){r=i2r[i];t=title[r];if(length(t)==0)t=r;printf("\t%s",t);}}for (k in nn){ printf("\n%s\t%d\t%d",k,nn[k],nn1[k]);for(i=1;i<=iMax;i++)if (nni[i]>0)printf("\t%d",n1[k,i]);if(nnn2>0){printf("\t\t%s\t%d\t%d",k,nn[k],nn2[k],nn2[k]);for(i=1;i<=iMax;i++)if (nni[i]>0)printf("\t%d",n2[k,i]);}}}END{printf("\n");}'  > tmp/GeneFusion/t4.$MAGIC.d
  echo -n "### File $toto : " > $toto
  date >> $toto
  echo "### Table of candidate gene fusions in project $MAGIC. Left table: single read  support, right table: pair support" >> $toto
  cat tmp/GeneFusion/t4.$MAGIC.b | gawk '/^#/{print}' >> $toto

  if (! -e tmp/GeneFusion/$target.gene2intmap.txt) then
    cat tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt | gawk -F '\t' '{gene=$5;split($2,aa,":");chrom[gene]=aa[1];split(aa[2],bb,"-");a1=bb[1];a2=bb[2];if(a1>a2){a0=a1;a1=a2;a2=a0;}if(0+aa1[gene]==0){aa1[gene]=a1;aa2[gene]=a2;}if(a1<aa1[gene])aa1[gene]=a1;if(a2>aa2[gene])aa2[gene]=a2;}END{for (g in chrom)printf("%s\t%s\t%d\t%d\n",g,chrom[g],aa1[g],aa2[g]);}' | sort > tmp/GeneFusion/$target.gene2intmap.txt
  endif
  cat  tmp/GeneFusion/t4.$MAGIC.d | gawk -F '\t' '/^#/{print;}' >> $toto
  cat  tmp/GeneFusion/$target.gene2intmap.txt ZZZZZ tmp/GeneFusion/t4.$MAGIC.d | gawk -F '\t' '/^ZZZZZ/{zz++;next;}/^#/{next;}{if(zz+0<1){g=$1;gc[g]=$2;g1[g]=$3;g2[g]=$4;next;}}{split($1,aa,"__");ga=aa[1];gb=substr(aa[2],1,length(aa[2])-2);printf("%s(%s:%d-%d__%s:%d-%d)",$1,gc[ga],g1[ga],g2[ga],gc[gb],g1[gb],g2[gb]);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' | sort -k 2nr >> $toto


  \cp $toto RESULTS/GeneFusion


goto phaseLoop

#######################################################################################################################################
## Topology of the Fusion network
## 2021_03_26
# count objects at different support scales, then count numbers of pairs, of cliques, then the scaling dimension 
# count object per scale
## analyse des gene_fusion seqc2
tace MetaDB <<EOF
  query find project IS Transloc ; >run
  s -o tmp/GeneFusion/AB_CL.groups @
  query find project IS TranslocShort ; >run
  s -o tmp/GeneFusion/GeneFusionShort.runs @
EOF
foreach run (`cat  tmp/GeneFusion/AB_CL.groups`)
  set sample=`echo $run | gawk -F _ '/ABCCL/{print $0;next;}{print $2"_"$3;}'`
  echo $sample
end

if (-e _t) \rm _t _s
foreach run (`cat  tmp/GeneFusion/AB_CL.groups`)
  set sample=`echo $run | gawk -F _ '/ABCCL/{print $0;next;}{print $2"_"$3;}'`
  echo $sample >> _s
  wc  tmp/GeneFusion/$run/t1.gene_fusion.av.count
  cat tmp/GeneFusion/$run/t1.gene_fusion.av.count | gawk -F '\t' 'BEGIN{p["+"]="";p["-"]="anti.";}/^#/{next;}/^__/{next}{split ($1,aa,"__");g1=aa[1];g2=aa[2];k=length(g2);s1=substr(g2,k-1,1);s2=substr(g2,k,1);g2=substr(g2,1,k-2);printf("%s%s__%s%s\t%s\t%d\n",p[s1],g1,p[s2],g2,s,$2+$3);}' s=$sample >> _t
end
cat _t | gawk -F '\t' '{n[$1"\t"$2]+=$3;}END{for(g in n){k=n[g];if(k>=1)printf("%s\t%d\n",g,k);}}' | sort > tmp/GeneFusion/AB_CL.count

if (-e _ts) \rm _ts _ss
foreach run (`cat  tmp/GeneFusion/GeneFusionShort.runs`)
  set sample=`echo $run | sed -e 's/^RNA_//'`
  echo $sample >> _ss
  wc  tmp/GeneFusion/$run/t1.gene_fusion.av.count
  cat tmp/GeneFusion/$run/t1.gene_fusion.av.count | gawk -F '\t' 'BEGIN{p["+"]="";p["-"]="anti.";}/^#/{next;}/^__/{next}{split ($1,aa,"__");g1=aa[1];g2=aa[2];k=length(g2);s1=substr(g2,k-1,1);s2=substr(g2,k,1);g2=substr(g2,1,k-2);printf("%s%s__%s%s\t%s\t%d\n",p[s1],g1,p[s2],g2,s,$2+$3);}' s=$sample >> _ts
end
cat _ts | gawk -F '\t' '{n[$1"\t"$2]+=$3;}END{for(g in n){k=n[g];if(k>=1)printf("%s\t%d\n",g,k);}}' | sort > tmp/GeneFusion/GeneFusionShort.count

# get the Magic fusion events
set toto=RESULTS/GeneFusion/Fusion_stats.txt
echo -n "### $toto :" > $toto
date >> $toto
echo "#N\tGenes fused to N other genes\tGenes with N fusion supports\tPairs with N supports" >> $toto
cat tmp/GeneFusion/AB_CL.count  | gawk -F '\t' '{k=$3;if(k+0<5)next;split($1,aa,"__");g3=$1;g1=aa[1];g2=aa[2];gs[g1]+=k;gs[g2]+=k;gg[g1]++;gg[g2]++;if(k>1000)k=1000;for (i = 1; i<k ;i++) np[i]++;}END{for(g in gs){for (i=1; i<1000 && i <= gs[g] ; i++) ng[i]++;if(gs[g]>=1000)ng[1000]++;}for(g in gg)for (i = 1; i <= gg[g] ; i++) ngg[i]++; for (i=1;i<1000 ;i++)printf("%d\t%d\t%d\t%d\n",i,ngg[i],ng[i],np[i]);}' >> $toto


# grab the number of reads in all genes
if (! -e  tmp/GeneFusion/gene2reads.txt) then
   zcat RESULTS/Expression/AceFiles/SEQC2.AceView.GENE.u.ace.gz | gawk '/^Gene/{g=$2;gsub(/\"/,"",g);split(g,aa,"(");g=aa[1];next;}/_SumOfAllReadsInProject/{printf("%s\t%s\t%s\t%s\n",g,$3,$6,$8);}' >  tmp/GeneFusion/gene2reads.txt
endif

# grab the number of neighbours relative to the number of supporting reads
set totog2r=RESULTS/GeneFusion/Fusion_expression_comparison.txt
echo -n "### $totog2r :" > $totog2r
date >> $totog2r
echo "# cGene\tFusions reads\tNeighbours(k>=100)\tExpression index\tReads supporting the gene\tkb supporting the gene" >> $totog2r
cat  tmp/GeneFusion/gene2reads.txt ZZZZZ  tmp/GeneFusion/AB_CL.count  | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){g2r[$1]=$2 "\t" $3 "\t" $4;next;}}{k=$3;if(k+0<100)next;split($1,aa,"__");g3=$1;g1=aa[1];g2=aa[2];gs[g1]+=k;gs[g2]+=k;gg[g1]++;gg[g2]++;}END{for(g in gg)printf("%s\t%d\t%d\t%s\n",g,gs[g],gg[g],g2r[g]);}' | sort -k 3nr | head -1000  >> $totog2r


################################################################################################################
################## OTHERS
## Get the fusion events from other teams participating in SEQC2
# data extraction in NCI_STAR_fusion
if (! -e RESULTS/GeneFusion/OTHERS/STARfusion_final/NCI_STARfusion.txt) then
  pushd RESULTS/GeneFusion/OTHERS/STARfusion_final
    foreach ff (*.abridged.tsv)
      set run=`echo $ff | gawk '{split($1,aa,"_star");print aa[1];}'`
      cat $ff | gawk -F '\t' '{printf("%s\t%s\tii\t%d\t%d\n",$1,run,$2,$3);}' run=$run > $run.tsf
    end
    cat *.tsf | grep FusionName | head -1 > xx
    cat *.tsf | grep -v FusionName | sort  >> xx
    mv xx NCI_STARfusion.txt
  popd
endif
pushd RESULTS/GeneFusion/OTHERS/Simon_AGLR2_CupCake_LongGF_ACB_PacBio
  if (-e xx) \rm xx
  foreach ff ( *.tsv )
    set run=`echo $ff|sed -e 's/\.tsv//'`
    cat $ff | gawk -F '\t' '{gg=$1 "--" $2;pp=$4;printf("%s\t%s__%s\t%d\n",gg,run,pp,$3);}' run=$run >> xx
  end
  mv xx Simon.counts.txt

  set tati=Simon_AGLR2_titration.txt
  echo -n "## $tati : ">> $tati
  date >> $tati
  cat Simon.counts.txt  | gawk -F '\t' '{gg=$1;m=$2;k=$3;split(m,aa,"_");s=index("ACB",aa[2]);m=aa[3]"_"aa[5];ng[gg]+=k;mm[m]+=k;nz[gg,m]+=k;nzs[gg,m,s]+=k;}END{printf("# Fusion\tMethod\tA\tC\tB\tTitrating\n");for (g in ng)for(m in mm){if(nz[g,m]>0){a=nzs[g,m,1];c=nzs[g,m,2];b=nzs[g,m,3];if(a>b){u=a;v=b;}else{u=b;v=a;}d=u-v;t="Flat";if (u<10)t="Not_enough_counts";else if(c>=v+d/10 && c <= v+9*d/10)t="Titrating";else if (c>2*(a+b)|| c<u/10)t="Incompatible";printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%s\n",g,ng[g],m,nz[gg,m],a,c,b,t);}}}' > $tati.1
  head -1 $tati.1 >> $tati
  cat $tati.1 | sort -k 2nr >> $tati
  \rm $tati.1
 cat $tati | gawk '/^#/{next;}{z=$3"\t"$8;nn[z]++;}END{for (k in nn)printf("%s\t%d\n",k,nn[k]);}' | sort
  
popd

## multicomptage des paires chez Anne
phaseOTHERS:
cat RESULTS/GeneFusion/OTHERS/2020Dec_Anne_STARFusionResultsAGLR_20200923.txt | sort -k 1,1 -k 2,2 -k 3,3 -k 8,8 | gawk -F '\t' '/^Panel/{next;}{gg=$8;split (gg,aa,"--");z= $1 $2 $3 $8;if(z==old)$10=0;old=z;k=$9+$10;g1=aa[1];g2=aa[2];if(g1=="" || g2=="")next;if(0 && g1>g2){g0=g1;g1=g2;g2=g0;};gg=g1"__"g2;ng[gg]+=k;i=index("0123456789XXYYYRRRAECDBLLL",$2);nz[gg,i]+=k;}END{for (g in ng){nz[g,13]=nz[g,19];nz[g,14]=nz[g,21];nz[g,15]=nz[g,23];printf("%s\tStarFusion_Anne\t%d",g,ng[g]);for(i=1;i<=12;i++)printf("\t");for(i=13;i<=15;i++)printf("\t%d",nz[g,i]);printf("\t\t\t");for(i=19;i<=23;i++)printf("\t%d",nz[g,i]);printf("\t\t\t\n");}}' > tmp/GeneFusion/gene2gene.support.other_groups.0

echo "### " >  tmp/GeneFusion/others.short.0
date >>  tmp/GeneFusion/others.short.0
echo toto | gawk '{printf("#Fusion\tMethod\tCumul");N=split("All_shorts\tAGLR1\tAGLR2\tROCR1\tROCR2\tILMR1\tILMR2\tILMR3\tGARR1\tTotal\tpolyA",aa,"\t");for(i=1;i<=N;i++){m=aa[i];printf("\t%s A\t%s E\t%s C\t%s D\t%s B\t%s TNA",m,m,m,m,m,m) ;}printf("\n");}' > tmp/GeneFusion/others.short.0

cat RESULTS/GeneFusion/OTHERS/2020Dec_Anne_STARFusionResultsAGLR_20200923.txt | sort -k 1,1 -k 2,2 -k 3,3 -k 8,8 | gawk -F '\t' '/^Panel/{next;}{t1=index("XXXX,AGLR1,AGLR2",$1);t1=int(t1/6);if(t1<1)next;t2=index("AECDB",$2);if(t2<1)next;t3=5*(t1)+t2-1;gg=$8;split (gg,aa,"--");z= $1 $2 $3 $8;if(z==old)$10=0;old=z;k=$9+$10;g1=aa[1];g2=aa[2];if(g1=="" || g2=="")next;gg=g1"__"g2;ng[gg]+=k;nz[gg,t2-1]+=k;nz[gg,t3]+=k;}END{for(g in ng){printf("%s\tStarFusion_Anne\t%d",g,ng[g]);for(i=0;i<15;i+=5){z=0;u=1;dj=1;if(i>=45)dj=2;for(j=dj;j<5;j+=dj){z1=nz[g,i+j]-nz[g,i+j-dj];if(z*z1<0)u=0;if(z1!=0)z=z1;}if(z==0)u=2;kk=0;for(j=0;j<5;j++){k=nz[g,i+j];kk+=k;printf("\t%d",k);}if(kk==0)u=3;printf("\t%s",substr("NTAA",u+1,1));}for(i=11;i<=55;i++)printf("\t");printf("\n");}}' >> tmp/GeneFusion/others.short.0

cat RESULTS/GeneFusion/OTHERS/20201218_Liz_PacBio_Fusion_Call_v3_sheet4_CL1-10_AceView_59fusionsOnly.txt | gawk -F '\t' '/^#/{next;}{g3=$6;for(i=26;i<=35;i++)if($i>=0)printf("%d\t%s\t%d\n",i-26,g3,$i);}' > _liz
cat RESULTS/GeneFusion/OTHERS/20201218_Liz_PacBio_Fusion_Call_v3_Sheet2_Pass_titration_1100AceView.txt  | cut -f 2,6,16  >> _liz
cat _liz | gawk -F '\t' '{gg=$2;gsub("--","__",gg);cl=index("0123456789XXACBRRRSSSSSLLL",$1);k=$3;ng[gg]+=k;nz[gg,cl]+=k;}END{for(g in ng){nz[g,24]=nz[g,13];nz[g,25]=nz[g,14];nz[g,26]=nz[g,15];nz[g,27]=nz[g,24];nz[g,28]=nz[g,25];printf("%s\tLiz_CupCake_PacBio\t%d",g,ng[g]);for(i=1;i<=15;i++)printf("\t%d",nz[g,i]);printf("\t\t\t\t\t\t\t\t");for(i=24;i<=29;i++)printf("\t%d",nz[g,i]);printf("\n");}}' >>  tmp/GeneFusion/gene2gene.support.other_groups.0
\rm _liz
cat RESULTS/GeneFusion/OTHERS/20210126_Liz_PacBio_Fusion_Call_GMAP_AGLR-2-only-v4_5206Gencode35.txt | cut -f 2,6,16 > _liz2
cat _liz2 | gawk -F '\t' '{gg=$2;gsub("--","__",gg);cl=index("0123456789XXACBRRRSSSSSLLL",$1);k=$3;ng[gg]+=k;nz[gg,cl]+=k;}END{for(g in ng){nz[g,24]=nz[g,13];nz[g,25]=nz[g,14];nz[g,26]=nz[g,15];nz[g,27]=nz[g,24];nz[g,28]=nz[g,25];printf("%s\tLiz_GMap_A2\t%d",g,ng[g]);for(i=1;i<=12;i++)printf("\t");for(i=13;i<=15;i++)printf("\t%d",nz[g,i]);printf("\t\t\t\t\t\t\t\t");for(i=24;i<=29;i++)printf("\t%d",nz[g,i]);printf("\n");}}' >>  tmp/GeneFusion/gene2gene.support.other_groups.0

cat RESULTS/GeneFusion/OTHERS/Pizzly_Pawel-3-AGLR2_subNA_20057.txt  | cut -f 2,6,15,16 | gawk -F '\t' '{gg=$2;gsub("--","__",gg);cl=index("0123456789XXYYYRRRAECDBLLL",$1);k=$3+$4;ng[gg]+=k;nz[gg,cl]+=k;}END{for (g in ng){nz[g,13]=nz[g,19];nz[g,14]=nz[g,21];nz[g,15]=nz[g,23];}for(g in ng){if(ng[g]>=1){printf("%s\tPizzly_Pawel_A2\t%d",g,ng[g]);for(i=1;i<=12;i++)printf("\t");for(i=13;i<=15;i++)printf("\t%d",nz[g,i]);printf("\t\t\t");for(i=19;i<=23;i++)printf("\t%d",nz[g,i]);printf("\t\t\t\n");}}}' >>  tmp/GeneFusion/gene2gene.support.other_groups.0
cat RESULTS/GeneFusion/OTHERS/Pizzly_Pawel-5-ROCR2_subNA_17106.txt | cut -f 2,6,15,16 | gawk -F '\t' '{gg=$2;gsub("--","__",gg);cl=index("0123456789XXYYYRRRAECDBLLL",$1);k=$3+$4;ng[gg]+=k;nz[gg,cl]+=k;}END{for (g in ng){nz[g,13]=nz[g,19];nz[g,14]=nz[g,21];nz[g,15]=nz[g,23];}for(g in ng){if(ng[g]>=1){printf("%s\tPizzly_Pawel_R2\t%d",g,ng[g]);for(i=1;i<=12;i++)printf("\t");for(i=13;i<=15;i++)printf("\t%d",nz[g,i]);printf("\t\t\t");for(i=19;i<=23;i++)printf("\t%d",nz[g,i]);printf("\t\t\t\n");}}}' >>  tmp/GeneFusion/gene2gene.support.other_groups.0


cat RESULTS/GeneFusion/OTHERS/Pizzly_Pawel-3-AGLR2_subNA_20057.txt  RESULTS/GeneFusion/OTHERS/Pizzly_Pawel-5-ROCR2_subNA_17106.txt | cut -f 1,2,6,15,16 |  gawk -F '\t' '/^#/{next;}{split($1,aa,"-");t1=index("XXXX,AGLR1,AGLR2,ROCR1,ROCR2,ILMR1,ILMR2,ILMR3,GARR1",aa[1]);t1=int(t1/6);if(t1<1)next;t2=index("AECDB",$2);if(t2<1)next;t3=5*(t1)+t2-1;gg=$3;split (gg,aa,"--");k=$4+$5;g1=aa[1];g2=aa[2];if(g1=="" || g2=="")next;gg=g1"__"g2;ng[gg]+=k;nz[gg,t2-1]+=k;nz[gg,t3]+=k;}END{for(g in ng){printf("%s\tPizzly_Pawel\t%d",g,ng[g]);for(i=0;i<45;i+=5){z=0;u=1;dj=1;if(i>=45)dj=2;for(j=dj;j<5;j+=dj){z1=nz[g,i+j]-nz[g,i+j-dj];if(z*z1<0)u=0;if(z1!=0)z=z1;}if(z==0)u=2;kk=0;for(j=0;j<5;j++){k=nz[g,i+j];kk+=k;printf("\t%d",k);}if(kk==0)u=3;printf("\t%s",substr("NTAA",u+1,1));}for(i=46;i<=55;i++)printf("\t");printf("\n");}}' >> tmp/GeneFusion/others.short.0

cat RESULTS/GeneFusion/OTHERS/20201207_MichaelSalmans_ILMN_FusionCall_sheet1hg38fusion.txt | cut -f 2,6,15,16 | gawk -F '\t' '{gg=$2;gsub("--","__",gg);cl=index("AECDB",$1);k=$3+$4;ng[gg]+=k;nz[gg,18+cl]+=k;}END{for(g in ng){if(ng[g]>=1){nz[g,13]=nz[g,19];nz[g,14]=nz[g,21];nz[g,15]=nz[g,23];printf("%s\tSalmans_Illumina\t%d",g,ng[g]);for(i=1;i<=12;i++)printf("\t");for(i=13;i<=15;i++)printf("\t%d",nz[g,i]);printf("\t\t\t");for(i=19;i<=23;i++)printf("\t%d",nz[g,i]);printf("\t\t\t\n");}}}' >>  tmp/GeneFusion/gene2gene.support.other_groups.0
cat RESULTS/GeneFusion/OTHERS/20201207_MichaelSalmans_ILMN_FusionCall_sheet1hg38fusion.txt | cut -f 1,2,6,15,16 | gawk -F '\t' '/^#/{next;}{split($2,aa,"_");t1=index("XXXX,AGLR1,AGLR2,ROCR1,ROCR2,ILMR1,ILMR2,ILMR3,GARR1",$1);t1=int(t1/6);if(t1<1)next;t2=index("AECDB",$2);if(t2<1)next;t3=5*(t1)+t2-1;gg=$3;split (gg,aa,"--");k=$4+$5;g1=aa[1];g2=aa[2];if(g1=="" || g2=="")next;gg=g1"__"g2;ng[gg]+=k;nz[gg,t2-1]+=k;nz[gg,t3]+=k;}END{for(g in ng){printf("%s\tSalmans_Illumina\t%d",g,ng[g]);for(i=0;i<45;i+=5){z=0;u=1;dj=1;if(i>=45)dj=2;for(j=dj;j<5;j+=dj){z1=nz[g,i+j]-nz[g,i+j-dj];if(z*z1<0)u=0;if(z1!=0)z=z1;}if(z==0)u=2;kk=0;for(j=0;j<5;j++){k=nz[g,i+j];kk+=k;printf("\t%d",k);}if(kk==0)u=3;printf("\t%s",substr("NTAA",u+1,1));}for(i=46;i<=55;i++)printf("\t");printf("\n");}}' >> tmp/GeneFusion/others.short.0

## multicomptage des paires chez Moha
cat RESULTS/GeneFusion/OTHERS/20201202_Mohammad_ILMN_FusionCatcher_FusionCall_v2.txt | cut -f 1,2,3,6,16,23 | sed -e 's/\r//' | sort | gawk -F '\t' '{z=$2 $3 $4;if(1       )$6=0;old = z; k=$5+$6;split ($4,aa,"--");g1=aa[1];g2=aa[2];if(g1=="" || g2=="")next;gg=g1"__"g2;s=index("AECDB",$2);s1=index("ACB",$2);ng[gg]+=k;nz[gg,s1+12]+=k;if(index($1,"PolyA_enriched")>0||index($1,"Total_riboDepleted")>0)nz[gg,s1+15]+=k;else nz[gg,s+18]+=k;}END{for(g in ng){printf("%s\tFusionCatcher_Sahraeian\t%d",g,ng[g]);for(i=1;i<=12;i++)printf("\t");for(i=13;i<=23;i++)printf("\t%d",nz[g,i]);printf("\t\t\t\n");}}' >>  tmp/GeneFusion/gene2gene.support.other_groups.0
cat RESULTS/GeneFusion/OTHERS/20201202_Mohammad_ILMN_FusionCatcher_FusionCall_v2.txt | cut -f 1,2,3,6,16 | gawk -F '\t' '/^#/{next;}{split($2,aa,"_");split($1,bb,"_");t1=index("XXXX,AGLR1,AGLR2,ROCR1,ROCR2,ILMR1,ILMR2,ILMR3,GARR1,Total,PolyA",bb[1]);t1=int(t1/6);if(t1<1)next;t2=index("AECDB",$2);if(t2<1)next;t3=5*(t1)+t2-1;gg=$4;split (gg,aa,"--");k=$5;g1=aa[1];g2=aa[2];if(g1=="" || g2=="")next;gg=g1"__"g2;ng[gg]+=k;nz[gg,t2-1]+=k;nz[gg,t3]+=k;}END{for(g in ng){printf("%s\tFusionCatcher_Sahraeian\t%d",g,ng[g]);for(i=0;i<55;i+=5){z=0;u=1;dj=1;if(i>=45)dj=2;for(j=dj;j<5;j+=dj){z1=nz[g,i+j]-nz[g,i+j-dj];if(z*z1<0)u=0;if(z1!=0)z=z1;}if(z==0)u=2;kk=0;for(j=0;j<5;j++){k=nz[g,i+j];kk+=k;printf("\t%d",k);}if(kk==0)u=3;printf("\t%s",substr("NTAA",u+1,1));}for(i=56;i<=55;i++)printf("\t");printf("\n");}}' >> tmp/GeneFusion/others.short.0

cat RESULTS/GeneFusion/OTHERS/2021March_Michael_StarFusionResults_PacBio_Actually_ILMN_reads.tsv | cut -f 2,3,4,19 | gawk -F '\t' '{split($4,aa,"_");gg=$1;gsub("--","__",gg);k=$2+$3;n=0;s=index("AECDB",aa[3]);nz[gg,18+s]+=k;ng[gg]+=k;}END{for(g in ng){nz[g,13]=nz[g,19];nz[g,14]=nz[g,21];nz[g,15]=nz[g,23];printf("%s\tMichael_StarFusion\t%d",g,ng[g]);for(i=1;i<=12;i++)printf("\t");for(i=13;i<=15;i++)printf("\t%d",nz[g,i]);printf("\t\t\t");for(i=19;i<=23;i++)printf("\t%d",nz[g,i]);printf("\t\t\t\n");}}'  >>  tmp/GeneFusion/gene2gene.support.other_groups.0
cat RESULTS/GeneFusion/OTHERS/20210413_Michael_STARFusion_Illumina_trimmed_reads.tsv | cut -f 2,3,4,19 | gawk -F '\t' '{split($4,aa,"_");gg=$1;gsub("--","__",gg);k=$2+$3;n=0;s=index("AECDB",aa[3]);nz[gg,18+s]+=k;ng[gg]+=k;}END{for(g in ng){nz[g,13]=nz[g,19];nz[g,14]=nz[g,21];nz[g,15]=nz[g,23];printf("%s\tMichael2_StarFusion\t%d",g,ng[g]);for(i=1;i<=12;i++)printf("\t");for(i=13;i<=15;i++)printf("\t%d",nz[g,i]);printf("\t\t\t");for(i=19;i<=23;i++)printf("\t%d",nz[g,i]);printf("\t\t\t\n");}}'  >>  tmp/GeneFusion/gene2gene.support.other_groups.0
cat RESULTS/GeneFusion/OTHERS/20210413_Michael_STARFusion_Illumina_trimmed_reads.tsv | cut -f 2,3,4,19 | sort -k 1,1 -k 4,4 | gawk -F '\t' '/^#/{next;}{split($4,aa,"_");t1=index("XXXX,AGLR1,AGLR2.ROCR1,ROCR2,ILMR1,ILMR2,ILMR3,GARR1",aa[2]);t1=int(t1/6);if(t1<1)next;t2=index("AECDB",aa[3]);if(t2<1)next;t3=5*(t1)+t2-1;gg=$1;split (gg,aa,"--");k=$2+$3;g1=aa[1];g2=aa[2];if(g1=="" || g2=="")next;gg=g1"__"g2;ng[gg]+=k;nz[gg,t2-1]+=k;nz[gg,t3]+=k;}END{for(g in ng){printf("%s\tStarFusion_Michael2\t%d",g,ng[g]);for(i=0;i<45;i+=5){z=0;u=1;dj=1;if(i>=45)dj=2;for(j=dj;j<5;j+=dj){z1=nz[g,i+j]-nz[g,i+j-dj];if(z*z1<0)u=0;if(z1!=0)z=z1;}if(z==0)u=2;kk=0;for(j=0;j<5;j++){k=nz[g,i+j];kk+=k;printf("\t%d",k);}if(kk==0)u=3;printf("\t%s",substr("NTAA",u+1,1));}for(i=46;i<=55;i++)printf("\t");printf("\n");}}' >> tmp/GeneFusion/others.short.0

cat RESULTS/GeneFusion/OTHERS/STARfusion_final/NCI_STARfusion.txt  | gawk -F '\t' '{gg=$1;gsub("--","__",gg);k=$4+$5;split($2,aa,"_");s=index("AECDB",aa[3]);s1=index("ACB",aa[3]);ng[gg]+=k;nz[gg,s1+12]+=k;if(index($2,"Total")>0 ||index($2,"PolyA")>0)nz[gg,s1+15]+=k;else nz[gg,s+18]+=k;}END{for(g in ng){printf("%s\tNCI_STARfusion\t%d",g,ng[g]);for(i=1;i<=12;i++)printf("\t");for(i=13;i<=23;i++)printf("\t%d",nz[g,i]);printf("\t\t\t\t\t\t\n");}}' >>  tmp/GeneFusion/gene2gene.support.other_groups.0
cat RESULTS/GeneFusion/OTHERS/STARfusion_final/NCI_STARfusion.txt  | sort -k 1,1 -k 2,2 | gawk -F '\t' '/^#/{next;}{split($2,aa,"_");t1=index("XXXX,AGLR1,AGLR2.ROCR1,ROCR2,ILMR1,ILMR2,ILMR3,GARR1,Total,PolyA",aa[2]);t1=int(t1/6);if(t1<1)next;t2=index("AECDB",aa[3]);if(t2<1)next;t3=5*(t1)+t2-1;gg=$1;split (gg,aa,"--");k=$4+$5;g1=aa[1];g2=aa[2];if(g1=="" || g2=="")next;gg=g1"__"g2;ng[gg]+=k;nz[gg,t2-1]+=k;nz[gg,t3]+=k;}END{for(g in ng){printf("%s\tStarFusion_NCI\t%d",g,ng[g]);for(i=0;i<45;i+=5){z=0;u=1;dj=1;if(i>=45)dj=2;for(j=dj;j<5;j+=dj){z1=nz[g,i+j]-nz[g,i+j-dj];if(z*z1<0)u=0;if(z1!=0)z=z1;}if(z==0)u=2;kk=0;for(j=0;j<5;j++){k=nz[g,i+j];kk+=k;printf("\t%d",k);}if(kk==0)u=3;printf("\t%s",substr("NTAA",u+1,1));}for(i=46;i<=55;i++)printf("\t");printf("\n");}}' >> tmp/GeneFusion/others.short.0

cat RESULTS/GeneFusion/OTHERS/Simon_AGLR2_CupCake_LongGF_ACB_PacBio/Simon.counts.txt | gawk -F '\t' '/^#/{next;}/_LongGF_/{gg=$1;gsub("--","__",gg);k=$3;split($2,aa,"_");pp[gg]=aa[5];s=index("ACB",aa[2]);ng[gg]+=k;nz[gg,12+s]+=k;s=index("AECDB",aa[2]);nz[gg,18+s]+=k;nz[gg,23+s]+=k;nz[gg,26+s]+=k;}END{for(g in ng){printf("%s\tSimon_LongGF_%s\t%d",g,pp[g],ng[g]);for(i=1;i<=12;i++)printf("\t");for(i=13;i<=15;i++)printf("\t%d",nz[g,i]);printf("\t\t\t");for(i=19;i<=23;i++)printf("\t%d",nz[g,i]);printf("\t\t\t\t\t\t\n");}}'  >>  tmp/GeneFusion/gene2gene.support.other_groups.0
cat RESULTS/GeneFusion/OTHERS/Simon_AGLR2_CupCake_LongGF_ACB_PacBio/Simon.counts.txt | gawk -F '\t' '/^#/{next;}/_cupcake_/{gg=$1;gsub("--","__",gg);k=$3;split($2,aa,"_");pp[gg]=aa[5];s=index("ACB",aa[2]);ng[gg]+=k;nz[gg,12+s]+=k;s=index("AECDB",aa[2]);nz[gg,18+s]+=k;nz[gg,23+s]+=k;nz[gg,26+s]+=k;}END{for(g in ng){printf("%s\tSimon_cupcake_%s\t%d",g,pp[g],ng[g]);for(i=1;i<=12;i++)printf("\t");for(i=13;i<=15;i++)printf("\t%d",nz[g,i]);printf("\t\t\t");for(i=19;i<=23;i++)printf("\t%d",nz[g,i]);printf("\t\t\t\t\t\t\n");}}'  >>  tmp/GeneFusion/gene2gene.support.other_groups.0
if (0) then
 cat RESULTS/GeneFusion/OTHERS/Simon_AGLR2_CupCake_LongGF_ACB_PacBio/Simon.counts.txt | gawk -F '\t' '/^#/{next;}/_cupcake_/{split($2,aa,"_");t1=index("XXXX,AGLR1,AGLR2.ROCR1,ROCR2,ILMR1,ILMR2,ILMR3,GARR1",aa[1]);t1=int(t1/6);if(t1<1)next;t2=index("AECDB",aa[2]);if(t2<1)next;t3=5*(t1)+t2-1;gg=$1;split (gg,aa,"--");k=$3;g1=aa[1];g2=aa[2];if(g1=="" || g2=="")next;gg=g1"__"g2;ng[gg]+=k;nz[gg,t2-1]+=k;nz[gg,t3]+=k;}END{for(g in ng){printf("%s\tSimon_cupcake\t%d",g,ng[g]);for(i=0;i<45;i+=5){z=0;u=1;dj=1;if(i>=45)dj=2;for(j=dj;j<5;j+=dj){z1=nz[g,i+j]-nz[g,i+j-dj];if(z*z1<0)u=0;if(z1!=0)z=z1;}if(z==0)u=2;kk=0;for(j=0;j<5;j++){k=nz[g,i+j];kk+=k;printf("\t%d",k);}if(kk==0)u=3;printf("\t%s",substr("NTAA",u+1,1));}for(i=46;i<=55;i++)printf("\t");printf("\n");}}' >> tmp/GeneFusion/others.short.0
  cat RESULTS/GeneFusion/OTHERS/Simon_AGLR2_CupCake_LongGF_ACB_PacBio/Simon.counts.txt | gawk -F '\t' '/^#/{next;}/_LongGF_/{split($2,aa,"_");t1=index("XXXX,AGLR1,AGLR2.ROCR1,ROCR2,ILMR1,ILMR2,ILMR3,GARR1",aa[1]);t1=int(t1/6);if(t1<1)next;t2=index("AECDB",aa[2]);if(t2<1)next;t3=5*(t1)+t2-1;gg=$1;split (gg,aa,"--");k=$3;g1=aa[1];g2=aa[2];if(g1=="" || g2=="")next;gg=g1"__"g2;ng[gg]+=k;nz[gg,t2-1]+=k;nz[gg,t3]+=k;}END{for(g in ng){printf("%s\tSimon_LongGF\t%d",g,ng[g]);for(i=0;i<45;i+=5){z=0;u=1;dj=1;if(i>=45)dj=2;for(j=dj;j<5;j+=dj){z1=nz[g,i+j]-nz[g,i+j-dj];if(z*z1<0)u=0;if(z1!=0)z=z1;}if(z==0)u=2;kk=0;for(j=0;j<5;j++){k=nz[g,i+j];kk+=k;printf("\t%d",k);}if(kk==0)u=3;printf("\t%s",substr("NTAA",u+1,1));}for(i=46;i<=55;i++)printf("\t");printf("\n");}}' >> tmp/GeneFusion/others.short.0
endif

cat tmp/GeneFusion/GeneFusionShort.count  | gawk -F '\t' '/^#/{next;}{split($2,bb,"_");t1=index("XXXX,AGLR1,AGLR2,ROCR1,ROCR2,ILMR1,ILMR2,ILMR3,GARR1,Total,PolyA",bb[1]);t1=int(t1/6);if(t1<1)next;t2=index("AECDB",bb[2]);if(t2<1)next;t3=5*(t1)+t2-1;gg=$1;k=$3;ng[gg]+=k;if(t1<=8)nz[gg,t2-1]+=k;nz[gg,t3]+=k;if(t1<=88&&k>=10)ok[gg]=1;}END{for(g in ng)if(ok[g]>0){printf("%s\tDanJean_Magic\t%d",g,ng[g]);for(i=0;i<55;i+=5){z=0;u=1;dj=1;if(i>=45)dj=2;for(j=dj;j<5;j+=dj){z1=nz[g,i+j]-nz[g,i+j-dj];if(z*z1<0)u=0;if(z1!=0)z=z1;}if(z==0)u=2;kk=0;for(j=0;j<5;j++)if(1){k=nz[g,i+j];kk+=k;printf("\t%d",k);}if(kk==0)u=3;printf("\t%s",substr("NTAA",u+1,1));}for(i=56;i<=55;i++)printf("\t");printf("\n");}}' | sort -k 3nr >> tmp/GeneFusion/others.short.0

set tutu8=RESULTS/GeneFusion/Fusion.short.stats.txt
echo -n "### $tutu8 :" > $tutu8
date >> $tutu8
cat tmp/GeneFusion/others.short.0 | gawk -F '\t' '{m=$2;mm[m]=1;for(i=1;i<=11;i++){j=3+6*i;z[m,i,$j]++;}}END{printf("#Method\tCumul T\tN\tAGLR1 T\tN\tAGLR2 T\tN\tROCR1 T\tN\tROCR2 T\tN\tILMR1 T\tN\tILMR2 T\tN\tILMR3 T\tN\tGARR1 T\tN\tTotal T\tN\tPolyA T\tPolyA N\t%% Titrating\tCumul\tAGLR1\tAGLR2\tROCR1\tROCR2\tILMR1\tILMR2\tILMR3\tGARR1\tTotal\tPolyA");for(m in mm){if(m!="Method"){printf("\n%s",m);for(i=1;i<=11;i++){printf("\t%d\t%d",z[m,i,"T"],z[m,i,"N"]);}printf("\t");for(i=1;i<=11;i++){printf("\t%.2f",100*z[m,i,"T"]/(.000001+z[m,i,"T"]+z[m,i,"N"]));}printf("\t");}}printf("\n");}' | sort >> $tutu8

set tutu9=RESULTS/GeneFusion/Fusion.shorts.counts.txt
echo -n "### $tutu9 :" > $tutu9
date >> $tutu9 
echo -n "# Line\tNumber of methods\tTotal support\t" >> $tutu9
head -1  tmp/GeneFusion/others.short.0  >> $tutu9
cat PROBES/hs.av.split_mrnas.txt  PROBES/REMAP/gencode.tr2ensg2refseqGene ZZZZZ  tmp/GeneFusion/others.short.0 ZZZZZ  tmp/GeneFusion/others.short.0  | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){if($2=="*"){ok[$1]=1;g2g[$1]=$5;}next;}}{split($1,aa,"__");if(ok[aa[1]]==1)aa[1]=g2g[aa[1]];if(ok[aa[2]]==1)aa[2]=g2g[aa[2]];$1=aa[1]"__"aa[2];if(zz==1){nm[$1]++;ng[$1]+=$3;next;}else printf("%d\t%d\t",nm[$1],ng[$1]);print;}' | sed -e 's/ /\t/g' | sort -k 1nr -k 2,2nr -k 3,3nr -k 4,4 | gawk '/^#/{next;}{line++;printf("%d\t",line);print;}'  | sed -e 's/\t0\t0\t0\t0\t0\tA/\t\t\t\t\t\t/g'  >> $tutu9



cat tmp/GeneFusion/gene2gene.support.other_groups.0 | grep __ | sort > tmp/GeneFusion/gene2gene.support.other_groups
cat tmp/GeneFusion/gene2gene.support.other_groups | grep BCAS3__ATXN7 > RESULTS/GeneFusion/toto.txt
cat tmp/GeneFusion/gene2gene.support.other_groups | grep BCAS4__BCAS3 > RESULTS/GeneFusion/toto.txt
cat RESULTS/GeneFusion/toto.txt

#### rename from ensembl to refseq using the split_mrna format, we wnat to move from a weird transcript name to a ensembl gene name then toa refseq gene name 
if (! -e PROBES/REMAP/gencode.tr2ensg2refseqGene) then
  zcat  PROBES/REMAP/gencode.v36.transcripts.fa.gz | gawk -F '|' '{split ($2,aa,".");printf ("%s\t%s\n",$6,aa[1]);]' | head
  zcat PROBES/REMAP/gencode.v36.transcripts.fa.gz | gawk -F '|' '{if(NF>6) {split ($2,aa,".");printf ("%s\t%s\n",$6,aa[1]);}}' | sort -u | gzip > PROBES/REMAP/gencode.tr2ensg.gz
  zcat PROBES/REMAP/gencode.tr2ensg.gz ZZZZZ.gz /home/mieg/AW/Human_DATA/20210427_gene2ensembl.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$2]=$1;next;}}{t=gg[$3];if(length(t)>1)printf("%s\t%s\t%s\n",t,$3,$4);}'  > PROBES/REMAP/gencode.tr2ensg2refseq
  cat PROBES/REMAP/gencode.tr2ensg2refseq ZZZZZ  PROBES/REMAP/RefSeq_109.title_line   | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){nm2t[$3]=$1;nz[$3]=$1"\t"$2"\t"$3;next;}}{split($1,aa,"|");g2=aa[3];t2=substr(aa[1],2);t1=nm2t[t2];if(length(t1)>1)printf("%s\t%s\n",g2,nz[t2]);}' | gawk '{if($1 == $2)next;}{printf("%s\t*\t0\t0\t%s\n",$2,$1);}'  > PROBES/REMAP/gencode.tr2ensg2refseqGene
endif
goto done

#############################################################################################################
####
## select the fusions with high coverage seen in only a few clones 
# if we do not remove k>10 in col N-2 (PacBio.raw) we add 20k bad fusions, but 200 if we rm col N-1
set toto=RESULTS/GeneFusion/Fusion_filtered.txt
cat tmp/GeneFusion/gene2gene.support.other_groups ZZZZZ tmp/GeneFusion/AB_CL.count | gawk -F '\t' 'BEGIN{cls="CL1-Brain_6,CL2-Breast_6,CL3-Cervix_6,CL4-Liver_6,CL5-Lipo_6,CL6-Blym_6,CL7-Tlym_6,CL8-Macr_6,CL9-Skin_6,CL10-Testis_6,B-Male_10,A-UHR_13,A_68,C_47,B_55,A_8NoCapture,C_8NoCapture,B_8NoCapture,A_33sC8,E_31sC8,C_31sC8,D_31sC8,B_31sC8,A_8l,C_8l,B_8l,A_PacBioccs,C_PacBioccs,B_PacBioccs,A_Nanop,C_Nanop,B_Nanop,PacBio.flnc3_39_ABCCL,PacBio.ccs3_39_ABCCL,PacBio.raw_27_ABCCL,Nanopore_48_ABCCL";N=split(cls,aa,",");for(i=1;i<=N;i++)cl2i[aa[i]]=i;}/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}{k=$3;if(k+0<1)next;cl=cl2i[$2];g=$1;if(cl<1)next;if(cl==20 || cl==22 || (cl<=15 && cl!=11 && cl != 12))gt[g]+=k;if(cl<=10){if(k>=10)n20[g]++;if(0 && k>=9)n9[g]++;}if(cl>=13 && k>10 && cl!=N-1)n50[g]++;z[g,cl]=k;}END{printf("# Fusion\tMethod\tCumul");for(i=1;i<=N;i++)printf("\t%s",aa[i]);for(g in ok)n20[g]++;for(g in n50){a=z[g,13];b=z[g,15];c=z[g,14];if(a>b){u=a;a=b;b=u}if(1 || b>2*a && c>a+d/4&&c<b+3*d/4)n20[g]++;}for(g in n20)if(ok[g]>0 || n9[g]<4)if(gt[g]>0){printf("\n%s\tDanJean_MagicExtended\t%d",g,gt[g]);for(i=1;i<=N;i++)printf("\t%d",0+z[g,i]);}printf("\n");}'  > $toto.a
cat ZZZZZ tmp/GeneFusion/AB_CL.count | gawk -F '\t' 'BEGIN{cls="CL1-Brain_6,CL2-Breast_6,CL3-Cervix_6,CL4-Liver_6,CL5-Lipo_6,CL6-Blym_6,CL7-Tlym_6,CL8-Macr_6,CL9-Skin_6,CL10-Testis_6,B-Male_10,A-UHR_13,A_68,C_47,B_55,A_8NoCapture,C_8NoCapture,B_8NoCapture,A_33sC8,E_31sC8,C_31sC8,D_31sC8,B_31sC8,A_8l,C_8l,B_8l,A_PacBioccs,C_PacBioccs,B_PacBioccs,A_Nanop,C_Nanop,B_Nanop,PacBio.flnc3_39_ABCCL,PacBio.ccs3_39_ABCCL,PacBio.raw_27_ABCCL,Nanopore_48_ABCCL";N=split(cls,aa,",");for(i=1;i<=N;i++)cl2i[aa[i]]=i;}/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}{k=$3;if(k+0<1)next;cl=cl2i[$2];g=$1;if(cl<1)next;if(cl==20 || cl==22 || (cl<=15 && cl!=11 && cl!=12))gt[g]+=k;if(cl<=10){if(k>=10)n20[g]++;if(0&&k>=9)n9[g]++;}if(cl>=13 && k>10 && cl!=N-1)n50[g]++;z[g,cl]=k;}END{printf("# Fusion\tMethod\tCumul");for(i=1;i<=N;i++)printf("\t%s",aa[i]);for(g in ok)n20[g]++;for(g in n50){a=z[g,13];b=z[g,15];c=z[g,14];if(a>b){u=a;a=b;b=u}if(1 || b>2*a && c>a+d/4&&c<b+3*d/4)n20[g]++;}for(g in n20)if(ok[g]>0 || n9[g]<4)if(gt[g]>0){printf("\n%s\tDanJean_Magic\t%d",g,gt[g]);for(i=1;i<=N;i++)printf("\t%d",0+z[g,i]);}printf("\n");}'  > $toto.b

cat ZZZZZ tmp/GeneFusion/AB_CL.count | gawk -F '\t' 'BEGIN{cls="CL1-Brain_6,CL2-Breast_6,CL3-Cervix_6,CL4-Liver_6,CL5-Lipo_6,CL6-Blym_6,CL7-Tlym_6,CL8-Macr_6,CL9-Skin_6,CL10-Testis_6,B-Male_10,A-UHR_13,A_68,C_47,B_55,A_8NoCapture,C_8NoCapture,B_8NoCapture,A_33sC8,E_31sC8,C_31sC8,D_31sC8,B_31sC8,A_8l,C_8l,B_8l,A_PacBioccs,C_PacBioccs,B_PacBioccs,A_Nanop,C_Nanop,B_Nanop,PacBio.flnc3_39_ABCCL,PacBio.ccs3_39_ABCCL,PacBio.raw_27_ABCCL,Nanopore_48_ABCCL";N=split(cls,aa,",");for(i=1;i<=N;i++)cl2i[aa[i]]=i;}/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}{k=$3;if(k+0<1)next;cl=cl2i[$2];g=$1;if(cl<1)next;if(cl==20 || cl==22 || (cl<=15 && cl!=11 && cl!=12))gt[g]+=k;if(cl<=12){if(cl!=11){if(k>=10)n20[g]++;if(k>=4)n9[g]++;}}if(0&&cl>=13 && k>10 && cl != N-1)n50[g]++;z[g,cl]=k;}END{printf("# Fusion\tMethod\tCumul");for(i=1;i<=N;i++)printf("\t%s",aa[i]);for(g in ok)n20[g]++;for(g in n50){a=z[g,13];b=z[g,15];c=z[g,14];if(a>b){u=a;a=b;b=u}if(1 || b>2*a && c>a+d/4&&c<b+3*d/4)n20[g]++;}for(g in n20)if(ok[g]>0 || n9[g]<4)if(gt[g]>0){printf("\n%s\tDanJean_Magic\t%d",g,gt[g]);for(i=1;i<=N;i++)printf("\t%d",0+z[g,i]);}printf("\n");}'  > $toto.b3
echo -n "### $toto.b4 :" > $toto.b4
date >> $toto.b4
echo "## Fusions with at least 10 support in a single cell line sample, but not more than 3 supports in more than 3 samples" >> $toto 
\cat $toto.b3 | head -12 | gawk '/^#/{printf("# Ordering\t Number of methods\t");print}' >> $toto.b4
cat $toto.b3  ZZZZZ $toto.b3   | gawk '/^#/{next}/^ZZZZZ/{zz++;next;}/novel/{next;}{if(zz<1){g=$1;nm[g]++;next;}}{g=$1;if(nm[g]>0)printf("%d\t",nm[g]);print;}'   > $toto.b5
cat PROBES/hs.av.split_mrnas.txt  PROBES/REMAP/gencode.tr2ensg2refseqGene ZZZZZ $toto.b5 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){if($2=="*"){ok[$1]=1;g2g[$1]=$5;}next;}}{split($2,aa,"__");if(ok[aa[1]]==1)aa[1]=g2g[aa[1]];if(ok[aa[2]]==1)aa[2]=g2g[aa[2]];$2=aa[1]"__"aa[2];next;}{print}' | sed -e 's/ /\t/g' | sort -k 1nr -k 5,5nr -k 2,2 -k 3,3 | gawk '{line++;printf("%d\t",line);print;}'   >> $toto.b4

################################################################################################################

echo -n "### $toto :" > $toto
date >> $toto
echo "## Fusions in different methods (for Magic: with at least 10 support in a single cell line sample, but not more than x supports in more than y samples)" >> $toto 
cat $toto.a | head -12 | gawk '/^#/{printf("# Ordering\t Number of methods\tFusion\tMethod\tSupporting fragments\tTitrating");for(i=5;i<=NF;i++)printf("\t%s",$i);printf("\n");}' >> $toto
cat PROBES/hs.av.split_mrnas.txt  PROBES/REMAP/gencode.tr2ensg2refseqGene ZZZZZ $toto.a $toto.b  tmp/GeneFusion/gene2gene.support.other_groups  | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){if($2=="*"){ok[$1]=1;g2g[$1]=$5;}next;}}{split($1,aa,"__");if(ok[aa[1]]==1)aa[1]=g2g[aa[1]];if(ok[aa[2]]==1)aa[2]=g2g[aa[2]];$1=aa[1]"__"aa[2];}{print}' | sed -e 's/ /\t/g' > $toto.c
cat $toto.c ZZZZZ $toto.c | gawk -F '\t' '/^#/{next}/^ZZZZZ/{zz++;next;}/novel/{next;}{m=$2;if (m=="DanJean_Magic"){if($3<10)next;}}{if(zz<1){g=$1;nm[g]++;next;}}{g=$1;if(nm[g]>0)printf("%d",nm[g]);for(i=1;i<=3;i++)printf("\t%s",$i);a=$19+$22+$27;c=$20+$24+$28;b=$21+$26+$29;d=$23;e=$25;de=d+e;;if(a>b){u=a;v=b;}else{u=b;v=a;}duv=u-v;du = 3 * sqrt(u) ; t="Flat";if (de > 10 && de > 2*(u+du)) t="Incompatible"; else if (u< 1.3 *v)t="Flat"; else if (u<10)t="Not_enough_counts"; else if (duv > du){if (c > v + duv/10 && c < u - duv/10) t = "Titrating"; else t="Incompatible";};printf("\t%s:%d:%d:%d",t,a,c,b);for(i=4;i<=NF;i++)printf("\t%s",$i);printf("\n");}'  | sort -k 1nr -k 2,2 -k 3,3 | gawk '/DanJean_MagicExtended/{if($1<2)next;}{line++;printf("%d\t",line);print;}'   >> $toto

################################################################################################################
set toto2=RESULTS/GeneFusion/Fusion_filtered.titrating_rates.txt
echo -n "### $toto2 :" > $toto2
date >> $toto2
echo "## In each proposed fusion, if the support in either A, B, or C sample is at least 10, the titration can be evaluated" >> $toto2
echo "## if A > B and A > 10 and d = A - B and B + d/10 <= C <= B + 9d/10, or vice-versa, the fusion is reported as titrating" >> $toto2
echo "## if C < A/10 or C < B/10 or C > 2 * (A+B), the counts are estimated incompatible with titration " >> $toto2
cat $toto | gawk -F '\t' '/^#/{next;}{m=$4;mm[m]++;split($6,aa,":");t=aa[1];if (t!="Not_enough_counts")mm10[m]++;tt[t]++;z[m,t]++;}END{printf("# Method\tNumber of fused gene pairs\tFusions with 10 or more supports in at least one sample\tPercent");for(t in tt)printf("\t%s",t);for(t in tt)printf("\t%s",t);for (m in mm){printf("\n%s\t%d\t%d\t%.2f",m,mm[m],mm10[m],100*mm10[m]/mm[m]);for(t in tt)printf("\t%d",z[m,t]);for(t in tt)printf("\t%.2f",100*z[m,t]/mm[m]);}printf("\n");}' > $toto2.a
head -1 $toto2.a >> $toto2
cat $toto2.a | gawk '/^#/{next;}{print}' | sort -k 2nr >> $toto2

echo "\n\n## Compare Moha and Magic"  >> $toto2
cat $toto | gawk -F '\t' '/^#/{next;}{g=$3;m=$4;split($6,aa,":");t=aa[1];if (m != "DanJean_MagicExtended" && m!= "Moha_NoPairs")next;tt[t]++;mm[m]++;gg[g]++;z[g,t,m]++;}END{for (g in gg)if(gg[g]==2){for (t1 in tt)for(t2 in tt)for (m1 in mm) for (m2 in mm) if(z[g,t1,m1]==1&&z[g,t2,m2]==1)nn[t1,t2,m1,m2]++;ggg++;}for (m1 in mm)for (m2 in mm)if (m1<m2){printf("%s\\%s",m1,m2);for (t1 in tt)printf("\t%s",t1);for(t2 in tt){printf("\n%s",t2);for(t1 in tt)printf("\t%d",nn[t2,t1,m2,m1]);}}printf("\nggg=%d\n",ggg);}' >> $toto2

echo "\n\n## Compare Moha and Magic"  >> $toto2
cat $toto | gawk -F '\t' '/^#/{next;}{g=$3;m=$4;split($6,aa,":");t=aa[1];if (m != "DanJean_MagicExtended" && m!= "Moha_FusionCatcher")next;tt[t]++;mm[m]++;gg[g]++;z[g,t,m]++;}END{for (g in gg)if(gg[g]==2){for (t1 in tt)for(t2 in tt)for (m1 in mm) for (m2 in mm) if(z[g,t1,m1]==1&&z[g,t2,m2]==1)nn[t1,t2,m1,m2]++;ggg++;}for (m1 in mm)for (m2 in mm)if (m1<m2){printf("%s\\%s",m1,m2);for (t1 in tt)printf("\t%s",t1);for(t2 in tt){printf("\n%s",t2);for(t1 in tt)printf("\t%d",nn[t2,t1,m2,m1]);}}printf("\nggg=%d\n",ggg);}' >> $toto2

echo "\n\n## Compare Moha and Moha"  >> $toto2
cat $toto | gawk -F '\t' '/^#/{next;}{g=$3;m=$4;split($6,aa,":");t=aa[1];if (m != "Moha_NoPairs" && m!= "Moha_FusionCatcher")next;tt[t]++;mm[m]++;gg[g]++;z[g,t,m]++;}END{for (g in gg)if(gg[g]==2){for (t1 in tt)for(t2 in tt)for (m1 in mm) for (m2 in mm) if(z[g,t1,m1]==1&&z[g,t2,m2]==1)nn[t1,t2,m1,m2]++;ggg++;}for (m1 in mm)for (m2 in mm)if (m1<m2){printf("%s\\%s",m1,m2);for (t1 in tt)printf("\t%s",t1);for(t2 in tt){printf("\n%s",t2);for(t1 in tt)printf("\t%d",nn[t2,t1,m2,m1]);}}printf("\nggg=%d\n",ggg);}' >> $toto2


## histo of compatibility
set toto3=RESULTS/GeneFusion/Fusion_compatibility.txt
echo -n "### $toto3 :" > $toto3
date >> $toto3
echo "## For each fusion method, count how many of its fusions are detected by n methods"  >> $toto3
cat $toto | gawk -F '\t' '/^#/{next;}{m=$4;n=$2;if(n>N)N=n;mm[m]++;z[m,n]++;}END{printf("2000000");for(m in mm)printf("\t%s",m);printf("\n1000000");for(m in mm)printf("\t%d",mm[m]);for(n=1;n<=N;n++){printf("\n%d",n);for(m in mm)printf("\t%d",z[m,n]);}printf("\n");}' | transpose | sort -k 2nr | transpose | sort -k 1nr >> $toto3

#################################################################################################################################
## enter the fusions in the GeneIndexDB database

## parse the captures
if (! -e GeneIndexDB/captured_genes.ace) then
  if (-e GeneIndexDB/captured_genes.txt) \rm   GeneIndexDB/captured_genes.txt
  foreach cap ($CAPTURES)
    set ff=TARGET/GENES/$cap.capture.av.gene_list
    if (! -e $ff) continue
    cat $ff | gawk '{printf("%s\t%s\n",$1,cap);}' cap=$cap >>  GeneIndexDB/captured_genes.txt
  end
  cat GeneIndexDB/captured_genes.txt | sort | gawk -F '\t' '{k=split($1,aa,"(");g=aa[1];if(g!= old){printf("\nGene %s\n",g);if(k==2){split(aa[2],bb,")");g2=bb[1];printf("AvGene %s\n",g2);}}old=g;printf("Capture %s\n",$2);}END{printf("\n\n");}' > GeneIndexDB/captured_genes.ace
  echo "pparse  GeneIndexDB/captured_genes.ace" | tace GeneIndexDB -noprompt
endif

# parse the support
if (1) then
  cat RESULTS/GeneFusion/Fusion_filtered.txt | gawk -F '\t' '/^#/{next;}{f=$3;m=$4;n=$5;k=split($6,aa,":");split(f,bb,"__");g1=bb[1];g2=bb[2];if(k>3)printf("Fusion %s\nGene1 %s\nGene2 %s\nMethod %s %d %s %d %d %d\n\n", f,g1,g2,m,n,aa[1],aa[2],aa[3],aa[4]);}' > tmp/GeneFusion/fusions.counts.ace 
  echo "pparse  tmp/GeneFusion/fusions.counts.ace" | tace GeneIndexDB -noprompt
endif

tace GeneIndexDB <<EOF
  select -o GeneIndexDB/map2centromere.txt  m,c1,c2 from m in ?map, c1 in m->Centromere_telomeres[2],c2 in c1[1] where c2
  find fusion
  select -o GeneIndexDB/fusion2map.txt f,g1,c1,a1,a2,g2,c2,b1,b2 from f in @,m in f->method where m=="StarFusion_Anne", g1 in f->gene1, c1 in g1->intmap, a1 in c1[1],a2 in c1[2],g2 in f->gene2, c2 in g2->intmap, b1 in c2[1],b2 in c2[2]
  select -o GeneIndexDB/fusion2map.txt f,g1,c1,a1,a2,g2,c2,b1,b2 from f in @,m in f->method where m=="DanJean_Magic", k in m[1] where k>30, ti in m[2] where ti=="Titrating", g1 in f->gene1, c1 in g1->intmap, a1 in c1[1],a2 in c1[2],g2 in f->gene2, c2 in g2->intmap, b1 in c2[1],b2 in c2[2]
  find fusion
  select -o GeneIndexDB/fusion2map.txt f,g1,c1,a1,a2,g2,c2,b1,b2 from f in @, g1 in f->gene1, c1 in g1->intmap, a1 in c1[1],a2 in c1[2],g2 in f->gene2, c2 in g2->intmap, b1 in c2[1],b2 in c2[2]
EOF

cat  GeneIndexDB/map2centromere.txt ZZZZZ GeneIndexDB/fusion2map.txt | gawk -F '\t' '/^ZZZZZ/{zz++;split("CF,CP",FF,",");next;}{if(zz<1){chr=$1c;cc[chr]=($2+$3)/2;next;}}/NULL/{next;}{f=$1;chr1=$3;a1=$4;a2=$5;chr2=$7;b1=$8;b2=$9;c1=cc[chr1];c2=cc[chr2];if(c1*c2==0)next;if ((a1-c1)*(a2-a1)>0)s1=1;else s1=2;if ((b1-c2)*(b2-b1)>0)s2=1;else s2=2;printf("%s_%s\t%s\t%d\t%s\t%d\t",FF[s1],FF[s2],chr1,c1,chr2,c2);print;}END{print pp,nn;}' > RESULTS/GeneFusion/Fusion_filtered.centric.txt
cat  RESULTS/GeneFusion/Fusion_filtered.centric.txt | gawk -F '\t' '{printf("Fusion %s\n%s\n\n",$6,$1);}' >  GeneIndexDB/fusion.centric.ace
cat  GeneIndexDB/map2centromere.txt ZZZZZ GeneIndexDB/fusion2map.txt | gawk -F '\t' '/^ZZZZZ/{zz++;split("CF,CP",FF,",");next;}{if(zz<1){chr=$1c;cc[chr]=($2+$3)/2;next;}}/NULL/{next;}{f=$1;chr1=$3;a1=$4;a2=$5;chr2=$7;b1=$8;b2=$9;c1=cc[chr1];c2=cc[chr2];if(c1*c2==0)next;if ((a1-c1)*(a2-a1)>0)s1=1;else s1=2;if ((b1-c2)*(b2-b1)>0)s2=1;else s2=2;printf("%s_%s\t%s\t%d\t%s\t%d\t",FF[s1],FF[s2],chr1,c1,chr2,c2);print;}END{print pp,nn;}' | tags




#################################################################################################################################
## 2021_04_22
## Testing d'ou viennent les fusions: les chromosomes, la transcription dans la cellule, le protocole moleculaire, le mapping
## si chromosome, on regarde si le chromosome rearrange a un ou zero ou deux centromeres, seul le cas 1 est viable, et les fusions seront introniques
## si machine a transcription de la cellule: les fusions seront introniques, et orientees AB plutot que BA mais sans respecter les chromosomes
## si le bruit vient du protocole, les molecules arrivent completes et epicees, donc la cassure n'a pas de raison de coincider avec les introns
## De plus on attend autant de AB que de BA car on n'est pas en cours de synthese, meme si la proba depend de la concntration totale en A et en B
## On calcule donc l'hitogramme lisse, au sens de la these 3e cycle de Danielle, AB/(AB+BA)
## On veut aussi mesurer la taille du glissement: centre sur 0 avec exponentielle descendante, ou avec un max en troispeut etre
## Enfin, pour une paire de mRNA fixee, a t'on un ou plusieurs points de cassure ?
##
## On veut aussi la frequence allelique dans A no-cap, dans B no-cap et dans le meilleur des cell
## il nous faux l'index du gene 1 et du gene 2, qui les capture, la couverture de l'exon donneur, de l'exon suivant local ou des introns locaux et de la fusion

\rm _bibi
foreach group (Sample_A_Nanop_A2R2R3_18lib Sample_A_PacBioccs_A2R3_7lib Sample_ABC_PacBio.ccs3_19lC2 Sample_ABC_PacBio.flnc3_19lC2 Sample_ABC_PacBio.raw_7lC1 Sample_ABC_Nanopore_28lC3 PacBio.fl_39_ABCCL PacBio.ccs3_39_ABCCL PacBio.raw_27_ABCCL Nanopore_48_ABCCL)
  cat tmp/GeneFusion/$group/t1.gene_fusion.av.tsf | gawk -F '\t' '{n=4;d=$14;if(d<=500 && d > -500)nn[d]+=n;}END{printf("## Group %s, duplication or gap in the alignemnet of the fusion read\n# Overlap (<0) or gap (>0)\tFusion reads\n",g);for(d=-500;d<=500;d++)printf("%s\t%d\t%d\n",g,d,nn[d]);}' g=$group >> _bibi
end
set toto4=RESULTS/GeneFusion/Duplications_gaps_in_fusion_reads.txt
echo -n "### $toto4 :" > $toto4
date >> $toto4
echo "## # Histogram of duplications (<0) or gaps (>0) in the alignemnet of the fusion reads"  >> $toto4
cat _bibi | gawk '/^#/{next;}{m=$1;mm[m]=1;z[m,$2]=$3;}END{printf("#");for (m in mm)printf ("\t%s",m);for(d=-500;d<=500;d++){printf("\n%d",d);for(m in mm)printf("\t%d",z[m,d]);}printf("\n");}' >> $toto4

# there is a peak around 84, where does it come from 
cat tmp/GeneFusion/$group/t1.gene_fusion.av.tsf | gawk -F '\t' '{n=4;d=$14;if (d>=81 && d <=87){nn[$7]+=n;nn[11]+=n;}}END{for(g in nn)printf("%s\t%d\n",g,nn[g]);}' | sort -k 2nr | head -12
# measure the ratio AB/(AB+BA) limited to the ++ case
# raw measure
set group=Sample_A_Nanop_A2R2R3_18lib
cat tmp/GeneFusion/$group/t1.gene_fusion.av.tsf | gawk -F '\t' '{if(index($1,"++")<1)next;n=4;g1=$6;g2=$10;if(g1<g2){g3 = g1 g2; gg1[g3]+=n;}else {g3 = g2 g1 ; gg2[g3]+=n;}gg[g3]+=n;}END{for(g in gg){p=gg1[g];m=gg2[g];if(p<m){q=p;p=m;m=q;}if(p+m>5)printf("%d\t%d\n",p,p+m);}}' | histo -smooth -title '# Ratio of fusion to reciprocal 100*AB/(AB+BA)' > RESULTS/GeneFusion/Fusion_orientation_ratio.txt

# distance between exact junction and closest exon boundary ? may be to be done in tricoteur.c and entered in a database


#################################################################################################################################
## grab the corresponding lists of reads
zcat  tmp/GeneFusion/gene2gene.count.filtered.gz  ZZZZZ.gz tmp/GeneFusion/*CL7*/t1.gene_fusion.av.txt.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ggg[$1]=1;next;}}{g1=$4;g2=$5;if(0 && g1>g2){g0=g1;g1=g2;g2=g0;}gg=g1"__"g2;if(ggg[gg]==1){;print;}}' | sort >  tmp/GeneFusion/gene2gene.support

cat tmp/GeneFusion/gene2gene.support.other_goups ZZZZZ tmp/GeneFusion/gene2gene.count | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ggg[$1]=1;ggg0[$1]=$0;next;}}{if(ggg[$1]==1){printf("%s\t%s\n",$0,ggg0[$1]);ok[$1]=1;}}END{for(gg in ggg)if(ok[gg]<1){for(i=1;i<=16;i++)printf("\t");printf("%s\n",ggg0[gg]);}}'

## try to grab the exact coordinates of the interesting cases
## on the case of CAS4 CAS3, we observe that the donor or CAS4 goes to any kind of acceptor os CAS3 downstream of the suppose break point
echo '#' > _t7
foreach run (`cat  tmp/GeneFusion/AB_CL.groups`)
  cat tmp/GeneFusion/$run/t1.gene_fusion.av.tsf  >> _t7
end
bin/tsf --merge --setSample  CL_AECDB -i _t7 > tmp/GeneFusion/t3.CL_AECDB.av.tsf
## get the exact donot acceptors
cat $toto.a ZZZZZ tmp/GeneFusion/t3.CL_AECDB.av.tsf | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){if($2=="DanJean_Magic")ok[$1]=1;next;}}{f=substr($1,1,length($1)-2);if(ok[f]==1)print;}' | head tmp tmp/GeneFusion/t3.CL_AECDB.av.filtered.tsf

##############################################################################
## 2021_04_02
## analyze detaillee de BCAS4__BCAS3++
## amazing, we can localize the genomic breakpoints by walking in nanopore and pacbio runs around coords seen in the wiggle
## we conclude that the fused chrom contains in successesion
# BCAS3 first exon , 20:49411431-49411710
# first intron up to 20:49411710-49430979  (then 2bp overlap)
# not exon 2 :       20:49434748-49434819  (then 3bp overlap)
# a reverse exon     20:58845848-58845400  
# then BCAS3 intron  17:59287702-592890---

zcat tmp/COUNT/PacBio.ccs3_ROCR3_CL02-F3/f2.1.hits.gz | gawk -F '\t' '{c=$11;a1=$12;a2=$13; if(a1>a2){a0=a1;a1=a2;a2=a0}if(c=="20" && x>a1&&x<a2)print;}' x=49430300 | cut -f 1

zcat tmp/COUNT/PacBio.ccs3_ROCR3_CL02-F3/f2.*.hits.gz | gawk -F '\t' '{c=$11;a1=$12;a2=$13; if(a1>a2){a0=a1;a1=a2;a2=a0}if(c=="20" && x>a1&&x<a2)print;}' x=49430300 | cut -f 1
zcat tmp/COUNT/PacBio.ccs3_ROCR3_CL02-F3/f2.*.hits.gz | grep n.18693#1

##############################################################################
## 2021_04_02
## Rationalize the RNA_seq breakpoints on donor acceptors if sliding or within 2bp
## deduce the genome coordinate and transfer to the earliest (gene.a .b .c) possible transcript)
## assume that the acceptor side breakpoint is downstream of the previous intron
## this gives us a localisation for the breakpoint in any given strain

bin/tricoteur -run CL_AECDB -t $target -geneFusionFile tmp/GeneFusion/t3.CL_AECDB.av.tsf -mrnaMap tmp/METADATA/av.MRNA.splicing.ace >  tmp/GeneFusion/t3.CL_AECDB.av.genome.tsf

# where splicedFusion has counts and the coords are the rationalized chrm/gene/mrna coords (to maintain compatibility) followed by chrom intervals 

## compute the ratios relative to the standard introns
## we expect the ratio on the donor side to be the ratio of the proba of the acceptors
## but the ration on the acceptor side to be the ratio of the strangth of the promotors of the 2 genes


##############################################################################
## 2021_04_02
## we would like to know is the proba of a donor acceptor pair is the prodiuct of its proba
## count for each intron, each donor, each acceptor the total number of supporting reads in human
## then for each intron contruct its donor and acceptor 8 base concensus, exctract the 4 base concensus xGyz (Gy = GT or Gc)
## export the counts for the donors, the acceptors and the intron
## check if the proba of the intron is the product of the proba of the donor and acceptor

if (! -e GeneIndexDB/parse_genome.done) then
  tace GeneIndexDB <<EOF
    pparse TARGET/Targets/$species.genome.fasta.gz
    save
EOF
touch GeneIndexDB/parse_genome.done

# document left and the right  motif, select 50 bases in exon then 50 in intron
if (! -e GeneIndexDB/donor_motifs.done) then
  foreach chrom ($chromSetAll)
    echo $chrom
    tace GeneIndexDB <<EOF
      date
      query find intron $chrom'_'*
      // select -o toto.motifs.$chrom.f ii,c,a1,a2,s,m1,m2 from ii in @, c in ii->IntMap, a1 in c[1], a2 in c[2] where a1 < a2, s in ?Sequence where s == c.name, m1 in DNA(s,a1-50,a1+49), m2 in DNA(s,a2-49,a2+50), 
      date
      query find intron $chrom*
      select -o toto.motifs.$chrom.r ii,c,a1,a2,s,m1,m2 from ii in @, c in ii->IntMap, a1 in c[1], a2 in c[2] where a1 > a2, s in ?Sequence where s == c.name, m1 in DNA(s,a1+50,a1-49), m2 in DNA(s,a2+49,a2-50), 
      date
      quit
EOF
    cat toto.motifs.$chrom.? | gawk '{printf("Intron %s\nDonor %s\nAcceptor %s\n\n",$1,$6,$7);}' > toto.motifs.$chrom.ace
    echo "pparse toto.motifs.$chrom.ace" | tace GeneIndexDB -noprompt
  end
  touch GeneIndexDB/donor_motifs.done
endif

cat toto.motifs.*.[fr] | gawk '{f=substr($6,51,2) "_" substr ($7,49,2) ; o=""; if (index("gt_ag,gc_ag,ct_ac,at_ag",f)==0)o="Other ";printf("Intron %s\n%s %s\n\n", $1,o, f);}' > toto.intron_feet.ace

# export the total support of donor/acceptor
    tace GeneIndexDB <<EOF
      date
      pparse  toto.intron_feet.ace
      query find intron 
      select -o toto.motifs_support.txt i,D,A,n,m1,m2 from i in @, D in i->D, A in i->A,n in i->RNA_seq, m1 in i->donor, m2 in i->acceptor
      date
      save
      quit
EOF

## Verify that most reads support introns of type gt_ag: it works, i get 98.2 gt_ag, 0.75 ct_ac, 0.75 gc_ag, 0.07 at_ac, 1 other
set toto=RESULTS/Introns_exons_polyA/donor_acceptor_motif.supports.txt
echo -n "### $toto : " > $toto
date >> $toto
echo "## Verify that most reads support introns of type gt_ag" >> $toto
echo "## Total number of reads supporting different combination of donor and acceptor sites" >> $toto
cat toto.motifs_support.txt | gawk -F '\t' '{n=$4;d=substr($5,51,2);a=substr($6,49,2);if(a=="nn" || d=="nn"||a==""||d=="")next;dd[d]+=n;aa[a]+=n;ii[d,a]+=n;}END{printf("-3\t-2");for (a in aa)printf("\t%s",a);printf("\n-2\t-1");for (a in aa)printf("\t%d",aa[a]);for(d in dd){printf("\n%s\t%d",d,dd[d]);for (a in aa)printf("\t%d",ii[d,a]);}printf("\n");}' | tab_sort -k 2n | transpose | tab_sort -k 2n | transpose  >> $toto

echo "\n\n\n## Percentage of reads supporting different combination of donor and acceptor sites" >> $toto
cat toto.motifs_support.txt | gawk -F '\t' '{n=$4;d=substr($5,51,2);a=substr($6,49,2);if(a=="nn" || d=="nn"||a==""||d=="")next;dd[d]+=n;aa[a]+=n;ii[d,a]+=n;iii+=n;}END{printf("-3\t-2");for (a in aa)printf("\t%s",a);printf("\n-2\t-1");for (a in aa)printf("\t%.4f",100*aa[a]/iii);for(d in dd){printf("\n%s\t%.4f",d,100*dd[d]/iii);for (a in aa)printf("\t%.4f",100*(ii[d,a]/iii -0*dd[d]*aa[a]/(iii*iii)));}printf("\n");}' | tab_sort -k 2n | transpose | tab_sort -k 2n | transpose >> $toto

## Go in the global wiggle and grab 4 numbers:  <a1 >a1 <a2 >a2, in each zone i take the median of 3 successive values, notice that i am not stranded
## to do that we need to create the shadow wiggle and then to populate it


set group=AB_CL
tace GeneIndexDB <<EOF
  select -o tmp/OR/d6.donors.list d from d in ?donor where d#IntMap
  select -o tmp/OR/d6.acceptors.list a from a in ?acceptor where a#IntMap
  status
EOF


set chrom=22

foreach chrom ($chromSetAll)
if (! -d tmp/OR/$group) mkdir tmp/OR/$group
  cat tmp/OR/d6.donors.list | gawk -F '\t' '{i=$1;split($1,aa,"_");c=aa[1];a=aa[3];fr=aa[4];if(c==chrom && fr=="f")printf("%s.DX\t1\t%s\t%d\t%d\t%s\n%s.DI\t2\t%s\t%d\t%d\t%s\n",i,c,a-30,a-1,i,i,c,a+1,a+30,i);}' chrom=$chrom >  tmp/OR/$group/d6.introns.$chrom.sponge
  cat tmp/OR/d6.donors.list | gawk -F '\t' '{i=$1;split($1,aa,"_");c=aa[1];a=aa[3];fr=aa[4];if(c==chrom && fr=="r")printf("%s.DI\t1\t%s\t%d\t%d\t%s\n%s.DX\t2\t%s\t%d\t%d\t%s\n",i,c,a-30,a-1,i,i,c,a+1,a+30,i);}' chrom=$chrom >>  tmp/OR/$group/d6.introns.$chrom.sponge
  cat tmp/OR/d6.acceptors.list | gawk -F '\t' '{i=$1;split($1,aa,"_");c=aa[1];a=aa[3];fr=aa[4];if(c==chrom && fr=="r")printf("%s.AX\t1\t%s\t%d\t%d\t%s\n%s.AI\t2\t%s\t%d\t%d\t%s\n",i,c,a-30,a-1,i,i,c,a+1,a+30,i);}' chrom=$chrom >>  tmp/OR/$group/d6.introns.$chrom.sponge
  cat tmp/OR/d6.acceptors.list | gawk -F '\t' '{i=$1;split($1,aa,"_");c=aa[1];a=aa[3];fr=aa[4];if(c==chrom && fr=="f")printf("%s.AI\t1\t%s\t%d\t%d\t%s\n%s.AX\t2\t%s\t%d\t%d\t%s\n",i,c,a-30,a-1,i,i,c,a+1,a+30,i);}' chrom=$chrom >>  tmp/OR/$group/d6.introns.$chrom.sponge

end

set group=Sample_SEQC2_AECDBCL_276_189sh87lg
set group=AGLR1_A_4libs

  set ww1=tmp/WIGGLEGROUP/$group/$chrom/R.chrom.frns.u.BF.gz
  set ww2=tmp/WIGGLEGROUP/$group/$chrom/R.chrom.frns.pp.BF.gz

  if (! -e $ww1  && ! -e $ww2) then
    echo "$ww1 not found"
    continue
  else if (! -e $ww1) then
    set ww=$ww2
  else if (! -e $ww2) then
    set ww=$ww1
  else
    set ww="$ww1,$ww2"
  endif

  set limit=1
  set mask=introns:tmp/OR/AB_CL/d6.introns.$chrom.sponge
  bin/geneelements -sponge $limit -spongeFile $mask  -sxxChromosome $chrom -wiggle $ww >  tmp/SPONGE/$group/introns.$chrom.u.ns.1
  cat tmp/SPONGE/$group/introns.$chrom.u.ns.1 | head -8 | gawk -F '\t' '/^introns/{split($3,aa,".");ii=aa[1];j=int(index("uDXDIAIAX",aa[2])/2);if(j==1 || j==2)dd[ii]=1;if(j>2)aa[ii]=1;z[ii,j]=int($12);}END{for(ii in dd)printf("Donor %s\nSponge %s %d %d\n\n",ii,run,z[ii,1]+0,z[ii,2]+0);for(ii in aa)printf("Acceptor %s\nSponge %s %d %d\n\n",ii,run,z[ii,3]+0,z[ii,4]+0);}' 
end


##############################################################################
## Intron fusions
    tace GeneIndexDB <<EOF
      date
      find intron 
      select -o intron_fusions.txt select ii,g1,g2,ln,t,r,n from ii in @,a in ii->A,d in ii->D,g1 in d->gene where g1,g2 in a->gene where g2 != g1,t in ii->type, ln in ii->length,r in ii->de_uno,n in r[1]
      date
      quit
EOF
# count gene fusions
cat intron_fusions.txt | grep gt_ag | cut -f 2,3 | sort -u | wc

##############################################################################
## confirmation des gene_fusion nanopore par illumina 
## cas plus ancien

# le concensus des coordonnees est analyse par tricoteur

foreach run (`cat MetaDB/$MAGIC/RunsList`)
  if (-e tmp/GeneFusion/$run/t1.gene_fusion.$target.txt.gz) then
    echo $run
    bin/tricoteur -run $run -t av -geneFusionFile tmp/GeneFusion/$run/t1.gene_fusion.$target.txt.gz -o   tmp/GeneFusion/$run/t1.gene_fusion.$target
  endif
end



zcat  tmp/GeneFusion/TotR5-100_S8/t1.gene_fusion.av.txt.gz  | head
#### Gene fusion candidates	2020-02-21_12:18:30	file=tmp/GeneFusion/TotR5-100_S8/f2.1.geneFusion.txt
#Run	Read	Fusion	Gene_A	Gene_B	Chrom_A	Chrom_B	Distance	Type	x1 A	x2 A	mRNA_A	from	to	x1 B	x2 B	mRNA_B	from	to	Distinct_supports	Support	Gene_A supports	Gene_B supports	score A	score B	Ali A	Ali B	c1 A	c2 A	c1 B	c2 B
TotR5-100_S8/f2.1	n.223121#1	ABCC4__skerdorbo++	ABCC4	skerdorbo	13	3	0	READ	1	74	ABCC4.aAug10	2595	2668	73	131	skerdorbo.cAug10	266	324	1	1	132	2	0	0	74	59	1	74	73	131
TotR5-100_S8/f2.1	n.223121#1	ABCC4__skerdorbo++	ABCC4	skerdorbo	13	3	0	READ	151	123	ABCC4.aAug10	2640	2668	124	66	skerdorbo.cAug10	266	324	1	1	132	2	0	0	29	59	123	151	66	124

# le fichier tmp/GeneFusion/t3.crn.av.txt donne la liste des gene fusions interessantes (couver > $minC dans au moins 1 run)
head tmp/GeneFusion/t3.crn.av.txt
### File tmp/GeneFusion/t3.crn.av.txt : Thu Feb 20 17:36:37 EST 2020
### Table of candidate gene fusions in project crn. Left table: single read  support, right table: pair support
# Run	Any	Supporting reads	coco	C1_S1	C2_S2	D1_S3	E1_S4	E2_S5	F1_S6	F2_S7	S8	S9	BL10_200_B2_F1_S7	QuartetD5_500_B3_C1_S9	BL10-200-2-BL10-50-1-TotR5-200_S5	D5-25-2-TotR10-50_S7	TotR-100-TotR8-100-BL10-100-2_S8		# Run	Any	Any2	coco	C1_S1	C2_S2	D1_S3	E1_S4	E2_S5	F1_S6	F2_S7	S8	S9	BL10_200_B2_F1_S7	QuartetD5_500_B3_C1_S9	BL10-200-2-BL10-50-1-TotR5-200_S5	D5-25-2-TotR10-50_S7	TotR-100-TotR8-100-BL10-100-2_S8
# Title	Any	Any1	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL		# Title	Any	Supporting read-pairsNULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL
RN7SL2andRN7SL3__LOC723809++(14:50320341-50329627__7:104535070-104567092)	111	111	151	0	0	0	35	38	14	13	0	0	11	0	RN7SL2andRN7SL3__LOC723809++	111	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0




phaseLoop:
 echo done

