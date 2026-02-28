#!/bin/tcsh -f
goto parse
foreach zone (`cat tmp/SNP_ZONE/ZoneList `)
  if (-e tmp/TSNP_DB/$zone/DanLi.parse.done) goto parse
end

# convert the dan li files to hg37 coordinates
foreach pp (AGLR1 AGLR2 ROCR1 ROCR2)
  cat 20220525_DanLi_$pp'_targetedRNA-seq_variant_list.txt' | gawk -F '\t' '/^Variant/{next}{w=$1;split($1,aa,"-");printf("%s\t%d\t%d\t%s-%s-%s-%s\n",aa[1],aa[2],aa[2]+1,w,$2,$3,$4);}'  | sort -u   > DanLi.$pp.hg38.bed
  ~/bin/liftOver DanLi.$pp.hg38.bed ~/VV/CODE/LIFTOVER/T2T/hg38ToHg19.over.chain.gz DanLi.$pp.hg37.bed DanLi.$pp.hg37.dead
end

# compare the lists of Dan Li and Wendell
 cat DanLi.ROCR2.*.hg37.bed | gawk '{if (length($3)<12)printf("%s-%s\n",$1,$3);}' | sort -u  > _u1  
zcat ../WENDELL_TRUTH/KnownPositives_hg19.vcf.wendellJones_GenBiol2021_40k.gz | gawk -F '\t' '{b=$4;bb=$5;if (length(b)<=4 && length(bb)<=4)printf("%s-%d-%s>%s\n",$1,$2,$4,$5);}' | sort  -u > _u2
wc _u[12]
cat _u[12] | sort -u | wc

# collate the VCF of all BRS snps
# run snp3.tcsh phase bed

foreach zone (`cat ../tmp/SNP_ZONE/ZoneList `)
  if (! -e  ../tmp/TSNP_DB/$zone/vcf.hg37.bed) then
    echo "missind file ../tmp/TSNP_DB/$zone/vcf.hg37.bed
  endif
end


# identify wendel true positives in vcf format + frequency
zcat ../WENDELL_TRUTH/KnownPositives_hg19.vcf.wendellJones_GenBiol2021_40k.gz | gawk -F '\t' '/^#/{next;}{b=$4;bb=$5;split($8,aa,";");split(aa[1],ff,"=");if (length(b)<=4 && length(bb)<=4)printf("%s\t%d\t%s\t%s\t%.2f\n",$1,$2,$4,$5,100*ff[2]);}' > wendell.true_positives.txt

foreach zone (`cat ../tmp/SNP_ZONE/ZoneList `)
  cat ../tmp/TSNP_DB/$zone/vcf.hg37.txt ../ZZZZZ wendell.true_positives.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){z="chr" $2 "-" toupper($3) "-" toupper($4) "-" toupper($5); m[z]++;mm[z,m[z]]=$1;next;}f=$5;z=$1"-"$2"-"$3"-"$4;if (m[z]>=1)for(i=1;i<=m[z] && i<2;i++){printf("%s\t%s\n",mm[z,i],f);}}'  >     wendell.true_positives.$zone.txt
  cat wendell.true_positives.$zone.txt | gawk -F '\t' '{printf("Variant \"%s\"\nWtrue %.2f\n\n",$1,$2);}' > wendell.true_positives.$zone.ace
  cat ../tmp/TSNP_DB/$zone/vcf.hg37.txt ../ZZZZZ wendell.true_positives.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){z="chr" $2 "-" toupper($3) "-" toupper($4) "-" toupper($5); m[z]++;mm[z,m[z]]=$1;next;}f=$5;z=$1"-"$2"-"$3"-"$4;if (m[z]>=1)for(i=1;i<=m[z] && i<200;i++){printf("%s\t%s\n",mm[z,i],f);}}'  >     wendell.true_positives.$zone.txt2
  cat wendell.true_positives.$zone.txt2 | gawk -F '\t' '{printf("Variant \"%s\"\nWtrue2 %.2f\n\n",$1,$2);}' > wendell.true_positives.$zone.ace2
end

# identify the danli to known varriants
foreach zone (`cat ../tmp/SNP_ZONE/ZoneList `)
  foreach pp (AGLR1 AGLR2 ROCR1 ROCR2)
    cat ../tmp/TSNP_DB/$zone/vcf.hg37.txt ../ZZZZZ DanLi.$pp.hg37.bed | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){z="chr" $2 "-" toupper($3) "-" toupper($4) "-" toupper($5); m[z]++;mm[z,m[z]]=$1;next;}split($4,aa,"-");z=$1"-"$2"-"aa[3]"-"aa[4];if (m[z]>=1)for(i=1;i<=m[z] && i<2;i++){printf("%s\t",mm[z,i],f);print;}}' > DanLi.$pp.hg37.$zone.txt
  end
end

# verify relative to the already registered danli
foreach zone (`cat ../tmp/SNP_ZONE/ZoneList `)
  tace ../tmp/TSNP_DB/$zone <<EOF
    query find Variant dan_li
    select -o DanLi/DanLi.known.$zone.list v from v in @
    quit 
EOF
end

# compare
cat DanLi.*.hg37.zo*.txt | cut -f 1 | sort -u > _v1
cat DanLi.known.*.list | sort -u > _v2
cat _v[12] | sort -u > _v3
wc _v[123]

# create an ace file
foreach zone (`cat ../tmp/SNP_ZONE/ZoneList `)
  echo ' ' > DanLi.hg37.$zone.ace
  foreach pp (AGLR1 AGLR2 ROCR1 ROCR2)
    cat DanLi.$pp.hg37.$zone.txt | gawk -F '\t' '{split($5,aa,"-");printf("Variant \"%s\"\nDan_Li %s-%s-%s-%s\nDanLi_counts %s %d %d %d Frequency %.2f\n\n", $1, aa[1],aa[2],aa[3],aa[4], pp,aa[7],aa[6]-aa[7],aa[6],100*aa[5] );}' pp=$pp >> DanLi.hg37.$zone.ace
  end
end

# false poitives BRS snps
zcat ../WENDELL_TRUTH/KnownNegatives_hg19.bed.WendellJones_GenBiol_2021_1.10M.gz  | gawk -F '\t' '/^#/{next;}{printf("%s-%s\n",$1,$2);}' > wendell.false_positives.txt
foreach zone (`cat ../tmp/SNP_ZONE/ZoneList `)
  cat ../tmp/TSNP_DB/$zone/vcf.hg37.txt ../ZZZZZ wendell.false_positives.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){z="chr" $2 "-" $3; m[z]++;mm[z,m[z]]=$1;next;}z=$1; if (m[z]>=1)for(i=1;i<=m[z] && i<2;i++){printf("Variant %s\n",mm[z,i]);}}'    >  brs.false_positive.$zone.list
  cat ../tmp/TSNP_DB/$zone/vcf.hg37.txt ../ZZZZZ wendell.false_positives.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){z="chr" $2 "-" $3; m[z]++;mm[z,m[z]]=$1;next;}z=$1; if (m[z]>=1)for(i=1;i<=m[z] && i<200;i++){printf("Variant %s\n",mm[z,i]);}}'  >  brs.false_positive.$zone.list2
end

# Fatigue monomodal
foreach zone (`cat ../tmp/SNP_ZONE/ZoneList `)
  if (-e Fatigue.monomodal.$zone.list) continue
  tace ~/Fatigue123/tmp/TSNP_DB/$zone <<EOF
    query find variant monomodal
    list -a -f Fatigue.monomodal.$zone.list
EOF
end

# Check which DanLi snps have not been remapped
# list the mapped ids
cat DanLi.zone*.ace | grep Dan_Li | gawk '{gsub("\"","",$2);print $2}' | sort -u > _a
# list the DanLi ids
cat DanLi.*.hg37.txt | gawk '{print $7}' | sort -u  > _b

cat _a _b | gawk '{n[$1]++;}END{for(k in n) if (n[k]==1)print k;}' > _c
wc _[abc]

cat _c ../ZZZZZ DanLi.*.hg37.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if (ok[$7]==1)print;}'


# parse the DanLi snps using snp3.tcsh


# Dan Li histograms

set DanLi_A1=20220525_DanLi_AGLR1_targetedRNA-seq_variant_list.txt  
set DanLi_A2=20220525_DanLi_AGLR2_targetedRNA-seq_variant_list.txt  
set DanLi_R1=20220525_DanLi_ROCR1_targetedRNA-seq_variant_list.txt  
set DanLi_R2=20220525_DanLi_ROCR2_targetedRNA-seq_variant_list.txt


cat $DanLi_A2 ZZZZZ $DanLi_R2 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;nextt;}if(ok[$1]<1)print;}' > DanLi_R2notA2
cat $DanLi_R2 ZZZZZ $DanLi_A2 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;nextt;}if(ok[$1]<1)print;}' > DanLi_A2notR2
cat $DanLi_A2 ZZZZZ $DanLi_R2 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;nextt;}if(ok[$1]==1)print;}' > DanLi_R2inA2
cat $DanLi_R2 ZZZZZ $DanLi_A2 | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;nextt;}if(ok[$1]==1)print;}' > DanLi_A2inR2

set toto=DanLi.AlleleFrequency.histos.txt
echo -n  "## Histogram of Dan Li frequencies" > $toto
date >> $toto
echo tutu | gawk '{printf("#Type\tTotal");for(i=0;i<=100;i++)printf("\t%d",i);printf("\n");}' >> $toto
foreach tt (A1 R1    A2 A2inR2 A2notR2  R2 R2inA2 R2notA2)
  echo -n $tt >> $toto
  if ($tt == A1) set ff=$DanLi_A1
  if ($tt == A2) set ff=$DanLi_A2
  if ($tt == R1) set ff=$DanLi_R1
  if ($tt == R2) set ff=$DanLi_R2
  if ($tt == A2inR2) set ff=DanLi_A2inR2
  if ($tt == A2notR2) set ff=DanLi_A2notR2
  if ($tt == R2inA2) set ff=DanLi_R2inA2
  if ($tt == R2notA2) set ff=DanLi_R2notA2
  cat $ff | gawk -F '\t' '/^chr/{nn++;f=int(100*$2+.5);hh[f]++;}END{printf("\t%d",nn);for(i=0;i<=100;i++)printf("\t%d",hh[i]);printf("\n");}' >> $toto
end


######################   june 23 2023: we presented our snp results, dan li sends his AGLR2 mask to find why he has only 2000 SNPs
set ff38=DanLi.AGLR2_liftover_all_merge_s_19to38.bed
# count overlapping region
cat $ff38 | gawk -F '\t' '{if($4=="+")print;}' | sort -k 1,1 -k 2,2n | gawk -F '\t' '{if ($3<x)print;x=$4}'  | wc
set ff37=DanLi.AGLR2.pulledBack.37
  ~/bin/liftOver $ff38 ~/VV/CODE/LIFTOVER/T2T/hg38ToHg19.over.chain.gz $ff37.bed $ff37.dead


 cat ../tmp/TSNP_DB/zoner.*/SNP_A.wendell.3.zoner.*.groups.SNP_summary.txt | gawk -F '\t' '/^#/{next;}{if (index($29,"ROCR2_A_4libs")>0)printf( "chr%s\t%d\tSNP\n",$3,$4)}' | sort -u | sort -k 1,1 -k 2,2n > _R2
 cat ../tmp/TSNP_DB/zoner.*/SNP_A.wendell.3.zoner.*.groups.SNP_summary.txt | gawk -F '\t' '/^#/{next;}{if (index($29,"AGLR2_A_4libs")>0)printf( "chr%s\t%d\tSNP\n",$3,$4)}' | sort -u | sort -k 1,1 -k 2,2n > _A2
cat $ff37.bed _R2 | sort -k 1,1 -k 2,2n > _R2bed
cat $ff37.bed _A2 | sort -k 1,1 -k 2,2n > _A2bed
cat _R2bed | gawk '{if($3=="SNP"){ns++;if ($1==c && $2>=a1 &&$2 <= a2)nn++;next;}c=$1;a1=$2;a2=$3;}END{print ns,nn;}'
cat _A2bed | gawk '{if($3=="SNP"){ns++;if ($1==c && $2>=a1 &&$2 <= a2)nn++;next;}c=$1;a1=$2;a2=$3;}END{print ns,nn;}'
wc _A2 _R2

###
# Liste des snps detectes 2 libs, comparaison
# new colonne detected
 cat ../tmp/TSNP_DB/zoner.*/SNP_A.wendell.3.zoner.*.groups.SNP_summary.txt | gawk -F '\t' '/^#/{next;}{if (index($29,"AGLR2_A_4libs")>0)printf( "%s\n",$8)}' | sort -u | sort -k 1,1 -k 2,2n > _a
 cat ../tmp/TSNP_DB/zoner.*/SNP_A.wendell.3.zoner.*.groups.detected_snps.tsf | grep GLR2_A_4libs | cut -f 1 | sort -u > _b
wc _a _b
# ls liste _A2 (3528) differe de a liste _a2 (3554) car on a garde 2 lignes pur le emme snp genomique affectant  2 genes disticnts, on a seuelement indentifies les snps affetcanbts plusieurs transcrits
# 3554 contre 3553 dans le rapport a case du vrai positif qui est monomodal

set dlA1=20220525_DanLi_AGLR1_targetedRNA-seq_variant_list.txt  
set dlA2=20220525_DanLi_AGLR2_targetedRNA-seq_variant_list.txt  
set dlR1=20220525_DanLi_ROCR1_targetedRNA-seq_variant_list.txt  
set dlR2=20220525_DanLi_ROCR2_targetedRNA-seq_variant_list.txt

cat $dlA2 $dlR2 $dlR2 | gawk '{n[$1]++;}END{for (s in n)nn[n[s]]++; for (k in nn) print k,nn[k];}'
cat $dlA2 $dlR2 $dlR2 | gawk '{n[$1]++;}END{for (s in n) if (n[s]==1)print s;}' > _bbA2
cat $dlA2 $dlR2 $dlR2 | gawk '{n[$1]++;}END{for (s in n) if (n[s]==2)print s;}' > _bbR2
cat $dlA2 $dlR2 $dlR2 | gawk '{n[$1]++;}END{for (s in n) if (n[s]==3)print s;}' > _bbA2R2
cat $dlA2 $dlR2 | cut -f 1 | sort -u > _bbAR
cat $dlA1 $dlA2 $dlR1 $dlR2 | cut -f 1 | sort -u > _bbAny
wc _bb*

# intersect with DanLi A2 mask
cat _bbA2 | gawk '/^ch*/{split($1,aa,"-");printf ("%s\t%d\tSNP\n",aa[1],aa[2]);}'   > _bbbA2
cat _bbR2 | gawk '/^ch*/{split($1,aa,"-");printf ("%s\t%d\tSNP\n",aa[1],aa[2]);}'   > _bbbR2
cat _bbA2R2 | gawk '/^ch*/{split($1,aa,"-");printf ("%s\t%d\tSNP\n",aa[1],aa[2]);}' > _bbbA2R2
cat _bbAR | gawk '/^ch*/{split($1,aa,"-");printf ("%s\t%d\tSNP\n",aa[1],aa[2]);}' > _bbbAR
cat _bbAny | gawk '/^ch*/{split($1,aa,"-");printf ("%s\t%d\tSNP\n",aa[1],aa[2]);}' > _bbbAny

cat $ff38 _bbbA2   | sort -k 1,1 -k 2,2n  > _bbbA2bed
cat $ff38 _bbbR2   | sort -k 1,1 -k 2,2n  > _bbbR2bed
cat $ff38 _bbbA2R2 | sort -k 1,1 -k 2,2n  > _bbbA2R2bed
cat $ff38 _bbbAR | sort -k 1,1 -k 2,2n  > _bbbARbed
cat $ff38 _bbbAny | sort -k 1,1 -k 2,2n  > _bbbAnybed

cat _bbbA2bed | gawk '{if($3=="SNP"){ns++;if ($1==c && $2>=a1 &&$2 <= a2)nn++;next;}c=$1;a1=$2;a2=$3;}END{print ns,nn;}'
cat _bbbA2R2bed | gawk '{if($3=="SNP"){ns++;if ($1==c && $2>=a1 &&$2 <= a2)nn++;next;}c=$1;a1=$2;a2=$3;}END{print ns,nn;}'
cat _bbbR2bed | gawk '{if($3=="SNP"){ns++;if ($1==c && $2>=a1 &&$2 <= a2)nn++;next;}c=$1;a1=$2;a2=$3;}END{print ns,nn;}'
cat _bbbARbed | gawk '{if($3=="SNP"){ns++;if ($1==c && $2>=a1 &&$2 <= a2)nn++;next;}c=$1;a1=$2;a2=$3;}END{print ns,nn;}'
cat _bbbAnybed | gawk '{if($3=="SNP"){ns++;if ($1==c && $2>=a1 &&$2 <= a2)nn++;next;}c=$1;a1=$2;a2=$3;}END{print ns,nn;}'

# Conclusion. Of 1064 R2 only SNPs in R2 only 680 belong to the A2 mask

cat mrnaR2only ZZZZZ ../tmp/TSNP_DB/zoner.*/SNP_A.wendell.3.zoner.*.groups.SNP_summary.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){mm[$1]=1;next;}split($8,aa,":");m=aa[1];if(mm[m]==1)print;}' > R2only.groups.SNP_summary.txt


# exported these files to Dan Li june 24 2023
 cat RESULTS/SNV/SnpA2_A.*.6.groups.SNP_summary.june16.txt | gawk '/^#/{print}' > _y
 cat RESULTS/SNV/SnpA2_A.*.3.groups.SNP_summary.june16.txt | gawk '/^#/{next;}/AGLR2_A_4libs/{print}' >> _y
 cat RESULTS/SNV/SnpA2_A.*.6.groups.SNP_summary.june16.txt | gawk '/^#/{next;}/AGLR2_A_4libs/{print}' | grep -v onomo >> _y
 mv _y RESULTS/SNV/SnpA2_A.for_Dan_Li.SNP_summary.jun24.txt 

 cat RESULTS/SNV/SnpR2_A.*.6.groups.SNP_summary.june16.txt | gawk '/^#/{print}' > _y
 cat RESULTS/SNV/SnpR2_A.*.3.groups.SNP_summary.june16.txt | gawk '/^#/{next;}/ROCR2_A_4libs/{print}' >> _y
 cat RESULTS/SNV/SnpR2_A.*.6.groups.SNP_summary.june16.txt | gawk '/^#/{next;}/ROCR2_A_4libs/{print}' | grep -v onomo >> _y
 wc _y
 mv _y RESULTS/SNV/SnpR2_A.for_Dan_Li.SNP_summary.jun24.txt
 

## compare the DanLi mask A2 to our mapping of the probes
pushd ../PROBES/A2
  zcat A2.genome.GRCh37.hits.gz  | gawk -F '\t' '/^#/{next;}{printf("%s\t%d\t%d\n",$11,$12,$13);}' > ~/SEQC2_2020/DanLi/Probes.A2.hg37
  zcat A2.genome.GRCh38.hits.gz  | gawk -F '\t' '/^#/{next;}{printf("%s\t%d\t%d\n",$11,$12,$13);}' > ~/SEQC2_2020/DanLi/Probes.A2.hg38
popd

cat ~/SEQC2_2020/DanLi/Probes.A2.hg38 | gawk '{nn++;x=$3-$2;if(x<0)x=-x;xx+=x;}END{print nn,xx,xx/nn;}'
cat ~/SEQC2_2020/DanLi/Probes.A2.hg38 | gawk '{if($2<$3){a1=$2;a2=$3;}else{a1=$3;a2=$2;}printf("%s\t%d\t%d\tPPP\n",$1,a1,a2);}' > A2.mapping.hg38
cat $ff38     A2.mapping.hg38 | sort -k 1,1 -k 2,2n >  A2.mapping.hg38.mix
cat  A2.mapping.hg38.mix | gawk -F '\t' '{if ($4=="PPP"){pc=$1;p1=$2;p2=$3;np++;}else {ac=$1;a1=$2;a2=$3;na++;}if(a1>p1)x1=a1;else x1=p1;if(a2<p2)x2=a2;else x2=p2;if(ac==pc && x2>x1){nn++;xx+=x2-x1+1;}}END{print np,na,nn,xx,xx/nn;}'


cat ~/SEQC2_2020/DanLi/Probes.A2.hg37 | gawk '{nn++;x=$3-$2;if(x<0)x=-x;xx+=x;}END{print nn,xx,xx/nn;}'
cat ~/SEQC2_2020/DanLi/Probes.A2.hg37 | gawk '{if($2<$3){a1=$2;a2=$3;}else{a1=$3;a2=$2;}printf("chr%s\t%d\t%d\tPPP\n",$1,a1,a2);}' > A2.mapping.hg37
cat $ff37.bed     A2.mapping.hg37 | sort -k 1,1 -k 2,2n >  A2.mapping.hg37.mix
cat  A2.mapping.hg37.mix | gawk -F '\t' '{if ($4=="PPP"){pc=$1;p1=$2;p2=$3;np++;}else {ac=$1;a1=$2;a2=$3;na++;}if(a1>p1)x1=a1;else x1=p1;if(a2<p2)x2=a2;else x2=p2;if(ac==pc && x2>x1){nn++;xx+=x2-x1+1;}}END{print np,na,nn,xx,xx/nn;}'

#####
# create  a bed file of the aligned probes
foreach gg  (GRCh37 GRCh38 T2T-CHM13v2)
  foreach pp (AGLR1 AGLR2 ILMR1 ILMR2 ILMR3 ROCR1 ROCR2 ROCR3)
    set ff=Probes_mapping/$pp/$pp.genome.$gg
    zcat $ff.hits.gz ZZZZZ.gz | gawk -F '\t' '/^#/{next;}{a1=$12;a2=$13;if(a1>a2){a0=a1;a1=a2;a2=a0;}printf("%s\t%d\t%d\n",$11,a1-dx,a2+dx);}' dx=0 | sort -k 1,1 -k 2,2n -k 3,3n | gawk '{if($1==old && a2>=$2-1){a2=$3;next} if(old)printf("%s\t%d\t%d\n",old,a1,a2);old=$1;a1=$2;a2=$3}' | gzip > $ff.bed.gz
    zcat $ff.bed.gz | gawk -F '\t' '{n++;nn+=$3-$2+1;}END{print gg,pp,"zone=",n,"\tprojection on genome=",nn}' gg=$gg pp=$pp
  end
end

foreach gg  (GRCh37 GRCh38 T2T-CHM13v2)
  foreach pp (AGLR1 AGLR2 ROCR1 ROCR2 ROCR3 ILMR1 ILMR2 ILMR3 )
    set ff=Probes_mapping/$pp/$pp.genome.$gg.hits.gz
    zcat $ff | gawk -F '\t' '/^#/{next;}{printf("%s\t%s\n",$1,$4);}' | sort -u | gawk -F '\t' '{n++;nn+=$2;}END{print gg,pp, "mapped probes=",n, "cumulated length=",nn}' gg=$gg pp=$pp
  end
end

foreach gg  (GRCh37 GRCh38 T2T-CHM13v2)
  foreach pp (AGLR1 AGLR2 ROCR1 ROCR2 ROCR3 ILMR1 ILMR2 ILMR3 )
    set ff=Probes_mapping/$pp/$pp.genome.$gg.hits.gz
    echo -n $gg,$pp, "nb mappings probes=" 
    zcat $ff | gawk -F '\t' '/^#/{next;}{print $1;}' | wc
  end
end
foreach gg  (GRCh37 GRCh38 T2T-CHM13v2)
  foreach pp (AGLR1 AGLR2 ROCR1 ROCR2 ROCR3 ILMR1 ILMR2 ILMR3 )
    set ff=Probes_mapping/$pp/$pp.genome.$gg.hits.gz
    echo -n $gg,$pp, "projected length="
    zcat $ff | gawk -F '\t' '/^#/{next;}{print $1;}' | sort -u | wc
  end
end

foreach pp (AGLR1 AGLR2 ILMR1 ILMR2 ILMR3 ROCR1 ROCR2 ROCR3)
  set qq=`echo $pp | gawk '{printf("%s%s",substr($1,1,1),substr($1,5,1));}'`
  cat tmp/METADATA/AGLR2.av.captured_genes.ace | gawk '/^Capture_touch/{if($2==qq)t++;next;}/^Capture/{if($2==qq)n++;next;}END{print pp,qq,n,t}' pp=$pp qq=$qq
end
foreach pp (A1 A2 I1 I2 I3 R1 R2 R3)
  set tt=AceView_2010.GRCh37
  set ff=~/ftp-SEQC/SEQC2/Capture_and_Gene_annotation/Targeted_Genes/$tt/av.$pp.well_covered.gene_list.txt 
  cat $ff | gawk '/^#/{next;}{nn++;}END{print pp,nn}' pp=$pp
  set ff=~/ftp-SEQC/SEQC2/Capture_and_Gene_annotation/Targeted_Genes/$tt/av.$pp.touched.gene_list.txt 
  cat $ff | gawk '/^#/{next;}{nn++;}END{print pp,nn}' pp=$pp
end
foreach pp (A1 A2 I1 I2 I3 R1 R2 R3)
  set tt=RefSeq_105.GRCh37
  set ff=~/ftp-SEQC/SEQC2/Capture_and_Gene_annotation/Targeted_Genes/$tt/RefSeq_105.$pp.well_covered.gene_list.txt 
  cat $ff | gawk '/^#/{next;}{nn++;}END{print pp,nn}' pp=$pp
  set ff=~/ftp-SEQC/SEQC2/Capture_and_Gene_annotation/Targeted_Genes/$tt/RefSeq_105.$pp.touched.gene_list.txt 
  cat $ff | gawk '/^#/{next;}{nn++;}END{print pp,nn}' pp=$pp
end

foreach pp (A1 A2 I1 I2 I3 R1 R2 R3)
  set tt=RefSeq_109.GRCh38
  set ff=~/ftp-SEQC/SEQC2/Capture_and_Gene_annotation/Targeted_Genes/$tt/RefSeq_109.$pp.well_covered.gene_list.txt 
  cat $ff | gawk '/^#/{next;}{nn++;}END{print pp,nn}' pp=$pp
  set ff=~/ftp-SEQC/SEQC2/Capture_and_Gene_annotation/Targeted_Genes/$tt/RefSeq_109.$pp.touched.gene_list.txt 
  cat $ff | gawk '/^#/{next;}{nn++;}END{print pp,nn}' pp=$pp
end
foreach pp (A1 A2 I1 I2 I3 R1 R2 R3)
  set tt=Gencode_36.GRCh38
  set ff=~/ftp-SEQC/SEQC2/Capture_and_Gene_annotation/Targeted_Genes/$tt/Gencode_36.$pp.well_covered.gene_list.txt 
  cat $ff | gawk '/^#/{next;}{nn++;}END{print pp,nn}' pp=$pp
  set ff=~/ftp-SEQC/SEQC2/Capture_and_Gene_annotation/Targeted_Genes/$tt/Gencode_36.$pp.touched.gene_list.txt 
  cat $ff | gawk '/^#/{next;}{nn++;}END{print pp,nn}' pp=$pp
end
foreach gg (_105 _109 .T2T-CHM13v2)
  foreach pp (A1 A2 R1 R2 R3 I1 I2 I3)
    echo -n $gg,$pp
    zcat $pp/$pp.RefSeq$gg.hits.gz | gawk '/^#/{next;}{if(length($9)>1)print $9}' | sort -u | wc
  end
end

###### comaprative table with Binshengs 2023_10_17

query find run Ali && SEQC2_title && IS RNA_Total_A_1	
select   t,s,r,srr, ns,nr, gs, gr1, ms, m1,m5,mode,median,av,m95,m99 from r in @, t in r->sorting_title, srr in r->SEQC2_SRR,  a in r->ali, s in r->SEQC2_title, nr in a->accepted[3], ns in s[1], gr in nr[2], gs in ns[2], m1 in a->Length_distribution_1_5_50_95_99_mode_av,m5 in m1[1],mode in m1[5],median in m1[2],av in m1[6],m95 in m1[3],m99 in m1[4], gr1=gr/1000000, ms in gs[3]


query find run Ali && SEQC2_title && IS Nanopore.titration_AGLR2_A-F3
query find run Ali && SEQC2_title && ! paired_end
select -o toto.r2s.txt t,s,r,srr, ns,nr, gs, gr1, ms, m1,m5,mode,median,av,m95,m99 from r in @, t in r->sorting_title, srr in r->SEQC2_SRR,  a in r->ali, s in r->SEQC2_title, nr in a->accepted[3], ns in s[1], gr in nr[2], gs in ns[2], m1 in a->Length_distribution_1_5_50_95_99_mode_av,m5 in m1[1],mode in m1[5],median in m1[2],av in m1[6],m95 in m1[3],m99 in m1[4], gr1=gr/1000000, ms in gs[3]

query find run Ali && SEQC2_title && paired_end
select -o toto.r2s.txt2 t,s,r,srr, ns,nr, gs, gr1, ms, m1,m5,mode,median,av,m95,m99 from r in @, t in r->sorting_title, srr in r->SEQC2_SRR,  a in r->ali, s in r->SEQC2_title, nr in a->accepted[3], ns in s[1], gr in nr[2], gs in ns[2], m1 in a->Fragment_length_1, m5 in a->Fragment_length_5, m95 in a->Fragment_length_95, m99 in a->Fragment_length_99, mode in a->Fragment_length_mode, median in a->Fragment_length_median, av in a->Fragment_length_average, gr1=gr/1000000, ms in gs[3]

echo "# Sorting title\tSRR id\tMagic\tSRR\tReads in SRR\tReads in Magic\tGb in SRR\tGb in Magic\tMedian in SRR\t1%\t5%\tMode\tMedian\tAverage\t95%\t99%" > RESULTS/Binsheng2Magic.txt
cat toto.r2s.txt2 >> RESULTS/Binsheng2Magic.txt                                                         
cat toto.r2s.txt >> RESULTS/Binsheng2Magic.txt                            

