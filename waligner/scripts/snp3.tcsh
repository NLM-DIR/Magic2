#!/bin/tcsh
# may 2023   analyse the tsmps and the BRS snps using snpsummary for the SEQC2_capture project

set phase=$1
set  zone=$2

echo "snp3 $1 $2"

if ($phase == 1) then 
  bin/snpsummary -db tmp/TSNP_DB/$zone -p $MAGIC -o RESULTS/SNV/$MAGIC.snp3a.$zone --snpType 0  # -minSnpFrequency 5 -minSnpCover 20
  goto done
endif

#kill snp not A2R2
if ($phase == kill) then
  set pp=SnpA2R2 
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF 
    find variant
    spush
    find gene
    select v from g in @ where g->capture == "A2" and g->capture == "R2", v in g->variant ;
    sminus
    spop
    kill
    save
    quit
EOF
  goto done
endif


if ($phase == G) then
# export all well covered counts c>=20
  set pp=$MAGIC
if (1) then
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF 
    read-models
    parse MetaDB/$MAGIC/runs.ace
    save
    select -o tmp/TSNP_DB/$zone/snp3.$pp.cmw20.txt v,r,c,m,w from v in ?Variant, r in v->BRS_counts where r->project == "$MAGIC" , c in r[1] where c >=20 , m in r[2], w in r[3]
    quit
EOF
endif
# select variants, well contrated both in runs S* and C*
if ($MAGIC == NASAZZZ) then
  cat tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.txt | gawk '{v=$1;if(v!=old){if (fMaxC>80 && fMinC<10 && fMaxS>80 && fMinS<10)print old;fMaxC=0;fMinC=100;fMaxS=0;fMinS=100;old=v;}c=$3;m=$4;w=$5;if (c<20)next;r=substr($2,1,1);f=100*m/c;if (r=="C"){if(f>fMaxC)fMaxC=f;if(f<fMinC)fMinC=f;}if (r=="S"){if(f>fMaxS)fMaxS=f;if(f<fMinS)fMinS=f;}}' >  tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.good_list
# get all their data
  cat tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.good_list ZZZZZ tmp/TSNP_DB/$zone/snp3.direct.cmw20.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if(ok[$1]==1){v=$1;r=$2;c=$3;m=$4;w=$5;f=100*m/c;printf("%s\t%s\tfi3\t%.2f\t%d\t%d\t%d\n",v,r,f,c,m,w);}}' >  tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.good.tsf
  cat   tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.good.tsf | gawk -F '\t' '{v=$1;r=$2;f=$4;if (v!= old){printf("\nVariant %s\n",v);old=v;} if(f>80)printf("mm %s\n",r);else if(f>40)printf("wm %s\n",r);if(f<10)printf("ww %s\n",r);}END{printf("\n");}' > tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.good_type.ace

  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF 
    pparse  tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.good_type.ace
    save
    quit
EOF

else
  cat tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.txt | gawk '{v=$1;if(v!=old){if (fMax>45 && fMin<5)print old;fMax=0;fMin=100;old=v;}c=$3;m=$4;w=$5;if (c<20)next;f=100*m/c;if(f>fMax)fMax=f;if(f<fMin)fMin=f;}' >  tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.good_list
  #  cat tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.txt | gawk '{v=$1;if(v!=old){if (fMax>80 && fMin<10)print old;fMax=0;fMin=100;old=v;}c=$3;m=$4;w=$5;if (c<20)next;f=100*m/c;if(f>fMax)fMax=f;if(f<fMin)fMin=f;}' >  tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.good_list
if ($MAGIC == NASA && -e tmp/TSNP_DB/$zone/snp3.Nano.cmw20.txt) then
   cat tmp/TSNP_DB/$zone/snp3.Nano.cmw20.txt | gawk '{v=$1;if(v!=old){if (fMax>40 && fMin<10)print old;fMax=0;fMin=100;old=v;}c=$3;m=$4;w=$5;if (c<20)next;f=100*m/c;if(f>fMax)fMax=f;if(f<fMin)fMin=f;}' >>  tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.good_list
endif

endif
# get all their data
  cat tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.good_list ZZZZZ tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}if(ok[$1]==1){v=$1;r=$2;c=$3;m=$4;w=$5;f=100*m/c;printf("%s\t%s\tfi3\t%.2f\t%d\t%d\t%d\n",v,r,f,c,m,w);}}' >  tmp/TSNP_DB/$zone/snp3.$MAGIC.cmw20.good.tsf

endif




if ($phase == p) then
  set pp=SnpA2R2 
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF 
    read-models
    find project $pp
    kill
    find groups
    kill
    find compare
    kill
    pparse MetaDB/$pp/runs.ace
    pparse MetaDB/$pp/groups.ace
    pparse MetaDB/$pp/samples.ace
    pparse MetaDB/$pp/compares.ace
    query find run "NULL"
    kill
    parse toto.project.ace   
    parse capture.geneAandB.ace
    find gene toto
    query find gene capture && ! IntMap && ! variant 
    kill
    save
    quit
EOF
  goto done
endif

if ($phase == M) then 
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF 
    query find variant monomodal == snpa2r2
    edit -D monomodal snpa2r2
    query find variant monomodal == Mole
    edit -D monomodal Mole
    select s,r,c,m from s in @, r in s->brs_counts where r like "SRR*", c in r[1], m in r[2] where m > 1 and 100*m < 20*c
    edit monomodal Mole
    save
    quit
EOF

  goto done
endif

if ($phase == G) then
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF 
    query find run union_of && project == $MAGIC
    edit -D project $MAGIC
    edit -D union_of
    pparse MetaDB/$MAGIC/runs.ace
    pparse MetaDB/$MAGIC/groups.ace
    save
    quit
EOF
  bin/tsnp --db_group_count  --db tmp/TSNP_DB/$zone --project $MAGIC 
  goto done
endif

if ($phase == R) then 
    set remap2g=remap2genome
# if -o filename is not provided, tsnp directly edits the database
    # bin/tsnp --db_$remap2g  tmp/METADATA/mrnaRemap.gz  --db tmp/TSNP_DB/$zone --force 
    bin/tace  tmp/TSNP_DB/$zone <<EOF
      read-models
      pparse tmp/METADATA/RvY.GENE.info.ace
      pparse tmp/METADATA/gtf.av.goodProduct.ace
      save
      quit
EOF
    echo -n "parse "
    ls -ls tmp/METADATA/gtf.av.goodProduct.ace
    echo -n "translate start: "
    date
    bin/tsnp --db_translate --db tmp/TSNP_DB/$zone  -p $MAGIC
    echo "translate done :"
    date
  goto done
endif

if ($phase == bed) then 
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF
    query find variant VCF
    select -o tmp/TSNP_DB/$zone/vcf.hg37.txt   s,m,a1,b,bb from s in @, m in s->VCF, a1 in m[1], b in m[2], bb in m[3]
EOF
 cat tmp/TSNP_DB/$zone/vcf.hg37.txt  | gawk -F '\t' '{printf("chr%s\t%d\t%d\t%s\t%s\t%s\n", $2, $3, $3+1,$1,$4,$5);}' | sort -k 1,1 -k 2,2n > tmp/TSNP_DB/$zone/vcf.hg37.bed
  /home/mieg/bin/liftOver tmp/TSNP_DB/$zone/vcf.hg37.bed /home/mieg/VV/CODE/LIFTOVER/T2T/hg19ToHg38.over.chain.gz  tmp/TSNP_DB/$zone/vcf.hg38.bed  tmp/TSNP_DB/$zone/vcf.hg38.dead
  cat  tmp/TSNP_DB/$zone/vcf.hg38.bed  | gawk -F '\t' '/:Ins:/{next}{printf("Variant %s\nVCF_hg38 %s-%s-%s-%s\n\n", $4,$1,$2,$5,$6);}' > tmp/TSNP_DB/$zone/vcf.hg38.ace
  cat  tmp/TSNP_DB/$zone/vcf.hg38.bed  | gawk -F '\t' '/:Ins:/{split($4,aa,":");printf("Variant %s\nVCF_hg38 %s-%s-%s-%s\n\n", $4,$1,$2,aa[4],aa[5]);}' >> tmp/TSNP_DB/$zone/vcf.hg38.ace
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF
    query find variant VCF_hg38
    edit -D VCF_hg38
    pparse tmp/TSNP_DB/$zone/vcf.hg38.ace
    save
    quit
EOF

  goto done
endif


if ($phase == D ) then 
  pushd tmp/TSNP_DB/$zone
    if (! -d dumpDir) mkdir dumpDir
    touch dumpDir/toto
    \rm dumpDir/*
    tace . <<EOF 
      dump -s dumpDir
      quit
EOF
  popd
  goto done
endif

if ($phase == E) then 
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF 
    s -o Global_SNP_LIST/brs.$zone.list v from v in ?Variant where v#VCF 
    quit
EOF
  goto done
endif

if ($phase == W) then 
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF 
    query find run PacBio.ROCR3.*
    kill
    query find run PacBio.*flnc*
    kill
    query find variant BRS_counts
    // select -o tmp/TSNP_DB/$zone/brs.fix v,r from v in @, r in v->brs_counts where ! r#project
    // pparse tmp/TSNP_DB/$zone/brs.fix.ace
    save
    quit
EOF
  goto done
endif

if ($phase == w) then 
  bin/tace tmp/TSNP_DB/$zone -no_prompt <<EOF 
    read-models
y
    pparse DanLi/DanLi.hg37.$zone.ace

    query find variant  WTrue
    spush
    pparse DanLi/wendell.true_positives.$zone.ace
    sminus
    spop
    edit -D Wtrue

    query find variant  WTrue2
    spush
    pparse DanLi/wendell.true_positives.$zone.ace2
    sminus
    spop
    edit -D Wtrue2

    query find variant  Wfalse
    spush
    key  DanLi/brs.false_positive.$zone.list
    query ! Wfalse
    edit Wfalse
    key  DanLi/brs.false_positive.$zone.list
    sminus
    spop
    edit -D Wfalse

    query find variant  Wfalse2
    spush
    key  DanLi/brs.false_positive.$zone.list2
    query ! Wfalse2
    edit Wfalse2
    key  DanLi/brs.false_positive.$zone.list2
    sminus
    spop
    edit -D Wfalse2

    key DanLi/Fatigue.monomodal.$zone.list
    edit monomodal Fatigue
    save
    quit
EOF
  goto done
endif

# SnpA2_A  SnpI2_goodSeq_A  SnpI2_lowQ_A  SnpR2_A SnpA2_B  SnpI2_goodSeq_B  SnpI2_lowQ_B  SnpR2_B SnpA2_C  SnpI2_goodSeq_C  SnpI2_lowQ_C  SnpR2_C


# SnpA2_A  SnpR2_A SnpA2_B  SnpR2_B (1 5 6 7 8) :: count Wtrue Wtrue2 Wfalse Wfalse2   (1 3   5 6 7 8     20 21 22 23) SNP_B
if ($phase == dc) then 

  \rm  tmp/TSNP_DB/$zone.out tmp/TSNP_DB/$zone.err
 foreach stype (3 4 5)
  if ($stype == 1) set titre=any
  if ($stype == 3) set titre=Wtrue
  if ($stype == 4) set titre=Wfalse
  if ($stype == 5) set titre=DanLiNonWendell

  if ($stype == 6) set titre=NonWendellNonDanLi
  if ($stype == 20) set titre=coding
  if ($stype == 21) set titre=UTR_3prime
  if ($stype == 22) set titre=A2G
  if ($stype == 23) set titre=G2A
  foreach pp (SnpA2_A  SnpR2_A)
#    bin/snpsummary -db tmp/TSNP_DB/$zone -o tmp/TSNP_DB/$zone/$pp.$titre.$stype.$zone --snpType $stype -e VIQSgd --doubleDetect -p $pp --histos --countLibs --unique
  end

  set pp=SnpA2R2
  set pp=$MAGIC  
  \rm  tmp/TSNP_DB/$zone/$pp.$titre.$stype.$zone.*
#   bin/snpsummary -db tmp/TSNP_DB/$zone -o tmp/TSNP_DB/$zone/$pp.$titre.$stype.$zone.groups     --snpType $stype  -e VIQTSgdDGR2 -p $pp --histos --countLibs --doubleDetect --titration  --unique --tsf --justDetected
   bin/snpsummary -db tmp/TSNP_DB/$zone -o tmp/TSNP_DB/$zone/$pp.$titre.$stype.$zone.groups.SNV.AF_Counts --snpType $stype  -e VIQTSgdDG2 -p $pp --histos --countLibs --doubleDetect --titration  --unique --tsf
   bin/snpsummary -db tmp/TSNP_DB/$zone -o tmp/TSNP_DB/$zone/$pp.$titre.$stype.$zone.libs.SNV.Counts  --snpType $stype  -e VIQTSR     -p $pp  --unique
   bin/snpsummary -db tmp/TSNP_DB/$zone -o tmp/TSNP_DB/$zone/$pp.$titre.$stype.$zone.libs.SNV.AF       --snpType $stype  -e VIQTSr     -p $pp  --unique
 end
 goto done
endif

if ($phase == cap) then 
  \rm  tmp/TSNP_DB/$zone.out tmp/TSNP_DB/$zone.err
 foreach stype (3)
  set capture=A1A2I2I3R1R2
  if ($stype == 1) set titre=any
  if ($stype == 3) set titre=$capture
  set pp=$MAGIC
  \rm  tmp/TSNP_DB/$zone/$pp.$titre.$stype.$zone.*
  bin/snpsummary -db tmp/TSNP_DB/$zone -o tmp/TSNP_DB/$zone/$pp.$titre.$stype.$zone.groups.SNV.AF_Counts     --snpType $stype  -e VIQTSgdDG2 -p $pp --histos --countLibs --doubleDetect --titration --capture $capture  --unique
  bin/snpsummary -db tmp/TSNP_DB/$zone -o tmp/TSNP_DB/$zone/$pp.$titre.$stype.$zone.libs.SNV.Counts  --snpType $stype  -e VIQTSR     -p $pp  --unique
  bin/snpsummary -db tmp/TSNP_DB/$zone -o tmp/TSNP_DB/$zone/$pp.$titre.$stype.$zone.libs.SNV.AF       --snpType $stype  -e VIQTSr     -p $pp --unique
 end
 goto done
endif



done:
  echo done
  exit 0

  
  foreach zone (`cat tmp/SNP_ZONE/ZoneList `)
    scripts/submit tmp/TSNP_DB/$zone "scripts/snp3.tcsh G $zone"
  end
 qusage 1

set pp=$MAGIC
#  NonWendellNonDanLi.6
foreach titre (Wtrue.3  Wfalse.4 DanLiNonWendell.5)
    foreach GR (libs.SNV.AF libs.SNV.Counts groups.SNV.AF_Counts)
      set toto=RESULTS/SNV/$pp.$titre.$GR.aug24
      echo -n "## SEQC2 Genes captured by A2 and R2 in project $pp limited to  $titre : " > $toto.txt
      echo -n "## $pp : SNV list counts and properties in project $pp restricted to genes captured by A2 and R2  limited to  $titre : " > $toto.txt
      date >> $toto.txt 
      cat tmp/TSNP_DB/zoner.*/$pp.$titre.zoner.*.$GR.SNP_summary.txt  | head -12 | tail -11 | gawk '/^#/{printf ("#\t"); print ;}'  >> $toto.txt
      cat tmp/TSNP_DB/zoner.*/$pp.$titre.zoner.*.$GR.SNP_summary.txt  | gawk '/^#/{next;}{print}' | tab_sort -k 3,3n -k 4,4n  | gawk '{n++;printf("%d\t", n); print;}' >> $toto.txt
    end
end


qusage 1

cat tmp/TSNP_DB/zoner.*/SnpA2_A.any22.5.zone*double_counts.txt | gawk '/libs/{next;}/^RNA/{nn[$1 "\t" $2]+=$3;}END{for(k in nn)printf("%s\t%d\n",k,nn[k]);}' | sort
cat  tmp/TSNP_DB/zoner.*/SnpA2_A.wendell.3.zoner.*.nDetectingLibsPerSnp.tsf | gawk -F '\t' '/^#/{next;}{for(i=4;i>=$4;i--)n[i]++;}END{print n[1],n[2],n[3],n[4];}'


# histo of snp covered at least 10000 times
cat RESULTS/SNV/SNP_A.any.1.groups.SNP_summary.june16.txt | head -200 | gawk '/^#/{print}' | transpose | grep -n Wende
cat RESULTS/SNV/SNP_A.any.1.groups.SNP_summary.june16.txt | head -200 | gawk '/^#/{print}' | transpose | grep -n 70
cat RESULTS/SNV/SNP_A.any.1.groups.SNP_summary.june16.txt | gawk -F '\t' '/^#/{next;}{x=$42 + 0;if (x==-10)n10++;if(x>=0 && $85>=10000){hh[int(100*x+.49)]++;n2++;}nn++;}END{print n10,n2,nn;for(i=0;i<=100;i++)print i,hh[i];}'
cat RESULTS/SNV/SNP_B.any.1.groups.SNP_summary.june16.txt | head -200 | gawk '/^#/{print}' | transpose | grep -n 55
cat RESULTS/SNV/SNP_B.any.1.groups.SNP_summary.june16.txt | gawk -F '\t' '/^#/{next;}{x=$42 + 0;if (x==-10)n10++;if(x>=0 && 80>=10000){hh[int(100*x+.49)]++;n2++;}nn++;}END{print n10,n2,nn;for(i=0;i<=100;i++)print i,hh[i];}'

#venn diagrams

cat tmp/TSNP_DB/zoner.*/SNP_A.wendell.3.zoner.*.groups.detected_snps.tsf | grep Total_A_4libs | grep Wtrue | cut -f 1 | sort -u > _tot
cat tmp/TSNP_DB/zoner.*/SNP_A.wendell.3.zoner.*.groups.detected_snps.tsf | grep PolyA_A_4libs | grep Wtrue | cut -f 1 | sort -u > _pA
 cat tmp/TSNP_DB/zoner.*/SNP_A.wendell.3.zoner.*.groups.detected_snps.tsf | grep AGLR2_A_4libs | grep Wtrue | cut -f 1 | sort -u > _a2
 cat tmp/TSNP_DB/zoner.*/SNP_A.wendell.3.zoner.*.groups.detected_snps.tsf | grep ROCR2_A_4libs | grep Wtrue | cut -f 1 | sort -u > r2
cat _tot  _a2 _a2  _r2 _r2 _r2 _r2 _r2 | gawk '{n[$1]++;}END{for (k in  n) nn[n[k]]++;for (k in nn) print k, nn[k];}'



set pp=$MAGIC

#  coding.20 UTR_3prime.21 A2G.22 G2A.23)
  # titration

set pp=$MAGIC
foreach titre (any.1 Wtrue.3  Wfalse.4 A1A2I2I3R1R2.3 A1A2I2I3R1R2_Wfalse.3)
  set toto=RESULTS/SNV/$pp.$titre.SNP_summary.june16
  echo -n "## $pp :  List of titratimg SNV by platform : " > $toto.titration.txt
  date >> $toto.titration.txt 
  cat tmp/TSNP_DB/zoner.*/$pp.$titre.zoner.*.titration.tsf | bin/tsf -I tsf --skip1 1  >>  $toto.titration.txt
end

foreach titre (any.1 wendell.3  Wfalse.6 A1A2I2I3R1R2.3 A1A2I2I3R1R2_Wfalse.3)
end


  # histo of our frequencies
foreach titre (Wtrue.3  Wfalse.4)
  echo -n "## $pp :  Histogram of SNV allele frequency project $pp restricted to genes captured by A2 and R2  limited to  $titre : " > $toto.histos.txt
  date >> $toto.histos.txt 
  cat tmp/TSNP_DB/zoner.*/$pp.$titre.zoner.*.groups.group_histos.tsf | bin/tsf -I tsf -O tabular | grep -v '##' | transpose > $toto.histos.txt2
  cat $toto.histos.txt2 | head -2 | tail -1 | gawk -F '\t' '{printf("\tSorting title"); for (i=1;i<=NF ; i++) {n=split ($i, aa, "___") ; if (n>=2) printf ("\t%s", aa[1], aa[2]) ;} printf ("\n") ;}' >>  $toto.histos.txt
  cat $toto.histos.txt2 | head -2 | tail -1 | gawk -F '\t' '{printf("\tAllele frequency percent") ; for (i=1;i<=NF ; i++) {n=split ($i, aa, "___") ; if (n>=2) printf ("\t%s", aa[2]) ;} printf ("\n") ;}' >>  $toto.histos.txt
  cat $toto.histos.txt2 | tail -n +3  >>  $toto.histos.txt
  \rm  $toto.histos.txt2


  # detectLibs
  echo -n "## $pp :  For each group, how many SNPs are seen as compatible or contradictory in n libs depending on the global allele frequency : ref low mid high pure (0,5,20,20,90) : " > $toto.detect.tsf
  date >> $toto.detect.tsf
  cat tmp/TSNP_DB/zoner.*/$pp.$titre.zoner.*.groups.detect.tsf | bin/tsf -I tsf -O tabular | grep -v '##' | transpose > $toto.detectLibs.txt2
  cat $toto.detectLibs.txt2 | head -2 | tail -1 | gawk -F '\t' '{printf("\tSorting title"); for (i=1;i<=NF ; i++) {n=split ($i, aa, "___") ; if (n>=2) printf ("\t%s", aa[1], aa[2]) ;} printf ("\n") ;}' >>  $toto.detect.tsf
  cat $toto.detectLibs.txt2 | head -2 | tail -1 | gawk -F '\t' '{printf("\tClass") ; for (i=1;i<=NF ; i++) {n=split ($i, aa, "___") ; if (n>=2) printf ("\t%s", aa[2]) ;} printf ("\n") ;}' >>  $toto.detect.tsf
  cat $toto.detectLibs.txt2 | tail -n +3  >>  $toto.detect.tsf
  \rm  $toto.detectLibs.txt2

  # calledLibs
  echo -n "## $pp :  For each group, how many SNPs are seen as compatible or contradictory in n libs depending on the global allele frequency : ref low mid high pure (0,5,20,20,90) : " > $toto.calledLibs.txt
  date >> $toto.calledLibs.txt 
  cat tmp/TSNP_DB/zoner.*/$pp.$titre.zoner.*.groups.calledLibs.tsf | bin/tsf -I tsf -O tabular | grep -v '##' | transpose > $toto.calledLibs.txt2
  cat $toto.calledLibs.txt2 | head -2 | tail -1 | gawk -F '\t' '{printf("\tSorting title"); for (i=1;i<=NF ; i++) {n=split ($i, aa, "___") ; if (n>=2) printf ("\t%s", aa[1], aa[2]) ;} printf ("\n") ;}' >>  $toto.calledLibs.txt
  cat $toto.calledLibs.txt2 | head -2 | tail -1 | gawk -F '\t' '{printf("\tClass") ; for (i=1;i<=NF ; i++) {n=split ($i, aa, "___") ; if (n>=2) printf ("\t%s", aa[2]) ;} printf ("\n") ;}' >>  $toto.calledLibs.txt
  cat $toto.calledLibs.txt2 | tail -n +3 | gawk -F '___' '{print $2;}'  >>  $toto.calledLibs.txt
  \rm  $toto.calledLibs.txt2

end


# histo of dan li frequencies
cat DanLi/DanLi.zoner.*.ace | gawk '/^DanLi_counts AGLR2/{x=int($7/10);print x}' | tags | sort -k 1n


# compare the coverage of magic/danli
cat tmp/TSNP_DB/zoner.*/SnpA2_A.any.3.zoner.*.SNP_summary.txt | gawk -F '\t' '{split($14,c1,":"); split($23,c2,":");n1=0+c1[2];n2=0+c2[2];if (n1*n2>0)printf("%d\t%d\t%d\n",n1,n2,2*n1/n2);}' | cut -f 3 | tags | sort -k 1n
cat tmp/TSNP_DB/zoner.*/SnpR2_A.any.3.zoner.*.SNP_summary.txt | gawk -F '\t' '{split($14,c1,":"); split($25,c2,":");n1=0+c1[2];n2=0+c2[2];if (n1*n2>0)printf("%d\t%d\t%d\n",n1,n2,2*n1/n2);}' | cut -f 3 | tags | sort -k 1n


  foreach zone (`cat tmp/SNP_ZONE/ZoneList `)
    cat tmp/TSNP_DB/$zone/brs.fix | gawk '{printf ("Variant %s\n-D BRS_counts %s\n\n", $1,$2);}' > tmp/TSNP_DB/$zone/brs.fix.ace
  end

#########

if (! -d tmp/TSNP_DB/ANY) then
  mkdir tmp/TSNP_DB/ANY
  pushd tmp/TSNP_DB/ANY
  mkdir database
  ln -s ../../../metaData/wspec
  tace . <<EOF
y
EOF
  popd
endif

set toto=tmp/TSNP_DB/ANY/_r
echo ' ' >  $toto
foreach zone (`cat tmp/SNP_ZONE/ZoneList `)
  foreach ff (`ls tmp/TSNP_DB/$zone/dumpDir/*.ace`)
    echo "pparse $ff" >> $toto
  end
end
echo save >> $toto
echo quit >> $toto

tace  tmp/TSNP_DB/ANY <  $toto &


######
# Venn diagrams of MAGIC SNPs
set ff=RESULTS/SNV/SnpA2R2.Wtrue.3.groups.SNP_summary.aug16.txt

cat $ff | gawk -F '\t' '/^#/{print}' | transpose | grep -n etected
cat $ff | gawk -F '\t' '/^#/{next;}{print $32}' | tags | grep A2




gene_name "HPD";
 annotation "4-hydroxyphenylpyruvate dioxygenase";
 Pfam "[Glyoxalase] Glyoxalase/Bleomycin resistance protein/Dioxygenase superfamily|[Glyoxalase_4] Glyoxalase/Bleomycin resistance protein/Dioxygenase superfamily|[Glyoxalase_3] Glyoxalase-like domain|[Glyoxalase_2] Glyoxalase-like domain";
 CDD "4-hydroxyphenylpyruvate dioxygenase. This protein oxidizes 4-hydroxyphenylpyruvate, a tyrosine and phenylalanine catabolite, to homogentisate. Homogentisate can undergo a further non-enzymatic oxidation and polymerization into brown pigments that protect some bacterial species from light. A similar process occurs spontaneously in blood and is hemolytic (see PMID:8000039). In some bacterial species, this enzyme has been studied as a hemolysin.";
 EC "EC1.13.11.27";
 KO "K00457";
 similar_genes "HPPD_DANRE(e-value:1.9e-90;
bit-score:330.0)";

zcat kumamushi_gtf_RvY_cds_scaf199.gtf.gz | cut -f 9 | grep gene_name | sed -e 's/\;/\n/g' -e 's/\"\"//g' | gawk '/^ transcript_id/{printf("\nSequence %s\n",$2);next;}/^ gene_name/{gsub("gene_name","LocusLink",$0);print;next;}/^ annotation/{gsub("annotation","Remark",$0);print;}/^ Pfam/{gsub("Pfam","Remark",$0);print;}/^ CDD/{gsub(/^ CDD/,"Remark",$0);print;next;}' > toto.ace
zcat kumamushi_gtf_RvY_cds_scaf199.gtf.gz | cut -f 9 | grep annotation | grep -v gene_name | sed -e 's/\;/\n/g' -e 's/\"\"//g' | gawk '/^ transcript_id/{printf("\nSequence %s\n",$2);next;}/^ annotation/{gsub("annotation","Remark",$0);print;}/^ Pfam/{gsub("Pfam","Remark",$0);print;}/^ CDD/{gsub(/^ CDD/,"Remark",$0);print;next;}' >> toto.ace
zcat kumamushi_gtf_RvY_cds_scaf199.gtf.gz | cut -f 9 | grep Pfam| grep -v annotation | grep -v gene_name | sed -e 's/\;/\n/g' -e 's/\"\"//g' | gawk '/^ transcript_id/{printf("\nSequence %s\n",$2);next;}/^ Pfam/{gsub("Pfam","Remark",$0);print;}/^ CDD/{gsub(/^ CDD/,"Remark",$0);print;}' >> toto.ace
zcat kumamushi_gtf_RvY_cds_scaf199.gtf.gz | cut -f 9 | grep CDD | grep -v Pfam| grep -v annotation | grep -v gene_name | sed -e 's/\;/\n/g' -e 's/\"\"//g' | gawk '/^ transcript_id/{printf("\nSequence %s\n",$2);next;}/^ CDD/{gsub(/^ CDD/,"Remark",$0);print;}' >> toto.ace
