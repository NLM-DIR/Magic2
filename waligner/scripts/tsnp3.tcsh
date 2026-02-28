#!bin/tcsh -f

set Strategy=$1
set zone=$2
set snp_type=".homozygous"
set snp_type=".filtered"
set snp_type=".differential"
set snp_type=""

# if ($Strategy == RNA_seq) set TTT=SNP_DB

#########################################################
## HACK used in CORONA to create the remap the snps

if (0 && $MAGIC == JUNK && ! -e tmp/METADATA/mrnaRemap.gz) then
  bin/tace tmp/TSNP_DB/$zone <<EOF
    select  -o tmp/METADATA/tsnp3.mrnaRemap t,m,x1,x2,map,a1,a2,g from t = "KT_RefSeq",m in ?mrna, map in m->intMap, a1 in map[1],a2 in map[2],x1=1,x2=a2-a1+1,g in m->gene
    quit
EOF
endif

#########################################################
## Flag the selected SNPs

if (-e tmp/TSNP_DB/$MAGIC.snp.selected.list) then 
  bin/tacembly tmp/TSNP_DB/$zone <<EOF
    read-models
    query find snp selected
    spush
    key  tmp/TSNP_DB/$MAGIC.snp.selected.list
    sminus
    edit Selected
    spop
    edit -D selected
    save
    quit
EOF
endif

#########################################################
## update the Runs, their groups, titles and other metadata

bin/tacembly tmp/TSNP_DB/$zone <<EOF
       read-models
       query find run project == $MAGIC
       edit -D project $MAGIC
       query union_of
       edit -D union_of
       pparse  MetaDB/$MAGIC/runs.ace
       pparse  MetaDB/$MAGIC/samples.ace
       save
       quit
EOF

if ($Strategy == RNA_seq) goto RNA_seq
#########################################################
## Export a nice table


  setenv ici `pwd`
  set toto=tmp/TSNP_DB/$zone/$MAGIC.tsnp3._r
  echo ' ' > $toto
  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    if (-e       tmp/TSNP/$run/$zone/$MAGIC.tsnp2c.val.tsf.gz) then
      echo `pwd`/tmp/TSNP/$run/$zone/$MAGIC.tsnp2c.val.tsf.gz >> $toto
    endif
  end

bin/tsnp --merge   -o tmp/TSNP_DB/$zone/$MAGIC.tsnp3.merged -f $toto
echo     "bin/tsnp -i tmp/TSNP_DB/$zone/$MAGIC.tsnp3.merged.snp_frequency.tsf --db tmp/TSNP_DB/$zone -p $MAGIC --db_report  -o tmp/TSNP_DB/$zone/tsnp3.$MAGIC --select NYC_JOR,Treagen_TX --minSnpFrequency $minSnpFrequency --minSnpCount $minSnpCount --wiggleDir tmp/WIGGLERUN"
# min zero, to parse the tsf table
time      bin/tsnp -i tmp/TSNP_DB/$zone/$MAGIC.tsnp3.merged.snp_frequency.tsf --db tmp/TSNP_DB/$zone -p $MAGIC --db_report  -o tmp/TSNP_DB/$zone/tsnp3.$MAGIC --select NYC_JOR,Treagen_TX --minSnpFrequency $minSnpFrequency --minSnpCount $minSnpCount --wiggleDir tmp/WIGGLERUN
# min 5 to export the nice table
# time    bin/tsnp -i tmp/TSNP_DB/$zone/$MAGIC.tsnp3.merged.snp_frequency.tsf --db tmp/TSNP_DB/$zone -p $MAGIC --db_report  -o tmp/TSNP_DB/$zone/tsnp3.$MAGIC --select NYC_JOR,Treagen_TX --minSnpFrequency $minSnpFrequency --minSnpCount $minSnpCount --wiggleDir tmp/WIGGLERUN

mv  tmp/TSNP_DB/$zone/$MAGIC.tsnp3.snp_frequency_table.txt _u$$

cat _u$$ | head -20 | gawk '/^#/{print}'                             >  tmp/TSNP_DB/$zone/$MAGIC.tsnp3.snp_frequency_table.txt
cat _u$$ | gawk '/^#/{next;}{print}' | scripts/tab_sort -k 8n -k 6n >>  tmp/TSNP_DB/$zone/$MAGIC.tsnp3.snp_frequency_table.txt

\cp  tmp/TSNP_DB/$zone/$MAGIC.tsnp3.snp_frequency_table.txt $MAGIC.$zone.snp_frequency_table.txt
\rm _u$$

exit 0

##########################################################
##########################################################
## It is not so clear to understand the rest of the code below


RNA_seq:
exit 0
#########################################################
## Export a nice table

if (! -e tmp/TSNP_DB/$zone/$MAGIC.mRNA.$zone.snp_frequency_table.txt2) then
  setenv ici `pwd`
  set toto=tmp/TSNP_DB/$zone/_r.snp3
  echo ' ' > $toto
  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    if (-e tmp/TSNP/$run/tsnp2.$MAGIC.$zone.val.txt.gz) then
      echo `pwd`/tmp/TSNP/$run/tsnp2.$MAGIC.$zone.val.txt.gz >> $toto
    else if (-e tmp/TSNP/$run/tsnp2.$MAGIC.$zone.val.txt) then
      echo `pwd`/tmp/TSNP/$run/tsnp2.$MAGIC.$zone.val.txt >> $toto
    else if (-e tmp/TSNP/$run/tsnp2.val.txt) then
      echo `pwd`/tmp/TSNP/$run/tsnp2.val.txt >> $toto
    endif
  end

  echo "bin/tsnp -f tmp/TSNP_DB/$zone/_r.snp3 --db tmp/TSNP_DB/$zone -p $MAGIC --db_report                -o tmp/TSNP_DB/$zone/$MAGIC.mRNA.$zone       --select NYC_JOR --wiggleDir tmp/WIGGLERUN --minSnpFrequency $minSnpFrequency --minSnpCount $minSnpCount --Reference_genome  $Reference_genome"
  time  bin/tsnp -f tmp/TSNP_DB/$zone/_r.snp3 --db tmp/TSNP_DB/$zone -p $MAGIC --db_report                -o tmp/TSNP_DB/$zone/$MAGIC.mRNA.$zone --select NYC_JOR  --minSnpFrequency $minSnpFrequency --minSnpCount $minSnpCount 

  \cp tmp/TSNP_DB/$zone/$MAGIC.mRNA.$zone.snp_frequency_table.txt RESULTS/SNV
  \cp tmp/SNP_DB/$zone/$MAGIC.mRNA.$zone.snp_wiggle_delta.txt  RESULTS/SNV

endif



## It is not so clear to undertstand the rest of the code below

exit 0


################ Find the 51 mer repeats
## i do not uinderstand the output format of this program, it only give the stats
# bin/jumpalign -t TARGET/Targets/$species.genome.fasta.gz -exportRepeats  -wMax 51 -o tmp/TSNP_DB/$zone/genome.51 -gzo


#########################################################
## Flag the genomic repeats
## i do not understand the formts and what this code does, we used it for the FDA project in 2018

if (0 && ! -e tmp/TSNP_DB/$zone/snp.Genome_51mer_repeats.ace && -e tmp/TSNP_DB/$zone/genome.51.repeats.gz) then

   bin/tacembly tmp/TSNP_DB/$zone <<EOF
     select -o tmp/TSNP_DB/$zone/snp_positions.txt  snp,chrom,pos from snp in ?variant, chrom in snp->IntMap, pos in chrom[1] where pos order_by +2+3
EOF
   gzip tmp/TSNP_DB/$zone/snp_positions.txt
   date

   gunzip -c tmp/TSNP_DB/$zone/snp_positions.txt.gz ZZZZZ.gz TARGET/Genome_repeats/Genome_51mer_repeats.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){z=$2 "\t" $3;chrom[$2]=1;snp[z]=$1;next;}}{if(chrom[$1]<1)next;for (i=$3;i<=$4;i++){z=$1 "\t" i;s=snp[z];if(s)printf("Variant %s\nGenome_51mer_repeats %d\n\n",s,$5);}}' >  tmp/TSNP_DB/$zone/snp.Genome_51mer_repeats.ace
   gunzip -c tmp/TSNP_DB/$zone/snp_positions.txt.gz ZZZZZ.gz TARGET/Genome_repeats/Genome_101mer_repeats.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){z=$2 "\t" $3;chrom[$2]=1;snp[z]=$1;next;}}{if(chrom[$1]<1)next;for (i=$3;i<=$4;i++){z=$1 "\t" i;s=snp[z];if(s)printf("Variant %s\nGenome_101mer_repeats %d\n\n",s,$5);}}' >  tmp/TSNP_DB/$zone/snp.Genome_101mer_repeats.ace
 

  gunzip -c ALL_SNP/snp_positions.txt ZZZZZ.gz TARGET/Genome_repeats/Genome_51mer_repeats.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){z=$2 "\t" $3;chrom[$2]=1;snp[z]=$1;next;}}{if(chrom[$1]<1)next;for (i=$3;i<=$4;i++){z=$1 "\t" i;s=snp[z];if(s)printf("Variant %s\nGenome_51mer_repeats %d\n\n",s,$5);}}' > ALL_SNP/snp.Genome_51mer_repeats.ace
  gunzip   -c ALL_SNP/snp_positions.txt ZZZZZ.gz TARGET/Genome_repeats/Genome_101mer_repeats.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){z=$2 "\t" $3;chrom[$2]=1;snp[z]=$1;next;}}{if(chrom[$1]<1)next;for (i=$3;i<=$4;i++){z=$1 "\t" i;s=snp[z];if(s)printf("Variant %s\nGenome_101mer_repeats %d\n\n",s,$5);}}' > ALL_SNP/snp.Genome_101mer_repeats.ace
endif

   date
   bin/tacembly tmp/TSNP_DB/$zone <<EOF
     read-models
     pparse tmp/TSNP_DB/$zone/snp.Genome_51mer_repeats.ace
     save
     quit
EOF
  date
  
endif

#########################################################
## Flag the selected SNPs

if (0 && ! -e   tmp/TSNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace) then
  echo -n "snp.prepare file  tmp/TSNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace "
  date

   # split the file because gawk chokes on very large files and acedb can reconstitute the pieces
    cat  tmp/TSNPH/$zone/$MAGIC.snp$snp_type.sorted   | sort -u > tmp/TSNPH/toto.$$.s
    split -l 20000  -a 6 tmp/TSNPH/toto.$$.s tmp/TSNPH/toto.$$.s.p.
    echo ' ' >  tmp/TSNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace

    set targetType=IntMap
    if ($Strategy == RNA_seq) set targetType=mRNA
    foreach ff (`ls  tmp/TSNPH/toto.$$.s.p.*`)
      echo $ff
      cat $ff | gawk -F '\t' -f scripts/snp2ace.awk format=txt targetType=$targetType >>   tmp/TSNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace
    end

    \rm  tmp/TSNPH/toto.$$  tmp/TSNPH/toto.$$.*  
endif

if (0 && ! -e tmp/TSNP_DB/$zone/snp.parseSnp.done) then
  if (-e   tmp/TSNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace) then
   echo -n "snp.parse file  tmp/TSNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace "
   date

    bin/tacembly tmp/TSNP_DB/$zone <<EOF
       pparse  tmp/TSNP_DB/$zone/$MAGIC.$zone.snp.sorted$snp_type.ace
       save
       quit
EOF

    touch  tmp/TSNP_DB/$zone/snp.parseSnp.done
  endif
endif




if (0 && ! -e  tmp/TSNP_DB/$zone/snp.differential.done2 && -e tmp/TSNPH/$zone/$MAGIC.differential_snps ) then
  
  cat  tmp/TSNPH/$zone/$MAGIC.differential_snps | gawk -F '\t' '/^#/{next}{split($1,aa,":"); f1=$5;f2=$12;if (f1<f2){if(f1<1)typ="Gained_in"; else typ="Up_in";}else{if(f2<1)typ="Lost_in";else typ="Down_in";}printf ("Variant \"%s:%s:%s\"\nDifferential %s \"%s\" \"%s\" %s %s\n\n",aa[1],aa[2],aa[3], typ, $9, $2, f1, f2);}' | grep -v "D_HR_surv_noMYCNA_NB399_exome_Normal" | grep -v Ghs3987 | grep -v Ghs4007 >  tmp/TSNP_DB/$zone/snp.differential.ace
  bin/tacembly tmp/TSNP_DB/$zone <<EOF
   query find variant Differential
   edit -D Differential
   read-models
   pparse  tmp/TSNP_DB/$zone/snp.differential.ace
   query find variant differential ; >Gene
   list -a -f tmp/TSNP_DB/$zone/$MAGIC.differential_snp2gene.list
   query find variant Gained_in ; >Gene
   list -a -f tmp/TSNP_DB/$zone/$MAGIC.differential_snp2gene.Gained.list
   query find variant Up_in  ; >Gene
   list -a -f tmp/TSNP_DB/$zone/$MAGIC.differential_snp2gene.Up.list
   query find variant Down_in ; >Gene
   list -a -f tmp/TSNP_DB/$zone/$MAGIC.differential_snp2gene.Down.list
   query find variant Lost_in  ; >Gene
   list -a -f tmp/TSNP_DB/$zone/$MAGIC.differential_snp2gene.Lost.list
   save
   quit
EOF
  touch  tmp/TSNP_DB/$zone/snp.differential.done
endif

laba:
if (0 && -e RESULTS/StLouisSNP/zone.stl.ace) then
  echo "pparse  RESULTS/StLouisSNP/zone.stl.ace"  | bin/tacembly tmp/TSNP_DB/$zone -no_prompt 
endif
echo "bin/snp -i tmp/TSNPH/$zone/$MAGIC.snp.sorted$snp_type -db tmp/TSNP_DB/$zone -project $MAGIC -db_frequencyTable -db_frequencyHisto -dropMonomodal -o  `pwd`/tmp/TSNP_DB/$zone/$MAGIC -Reference_genome  $Reference_genome"
      bin/snp -i tmp/TSNPH/$zone/$MAGIC.snp.sorted$snp_type -db tmp/TSNP_DB/$zone -project $MAGIC -db_frequencyTable -db_frequencyHisto -dropMonomodal -o  `pwd`/tmp/TSNP_DB/$zone/$MAGIC -Reference_genome  "$Reference_genome" 

if (-e tmp/TSNP_DB/$zone/$MAGIC.characteristic_snp.txt) \rm tmp/TSNP_DB/$zone/$MAGIC.characteristic_snp.txt

echo 'bin/snp -i tmp/TSNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt              -db tmp/TSNP_DB/$zone -project $MAGIC -db_prevalenceTable -dropMonomodal -o  `pwd`/tmp/TSNP_DB/$zone/$MAGIC -Reference_genome   "$Reference_genome "' 
      bin/snp -i tmp/TSNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt              -db tmp/TSNP_DB/$zone -project $MAGIC -db_prevalenceTable -dropMonomodal -o  `pwd`/tmp/TSNP_DB/$zone/$MAGIC -Reference_genome   "$Reference_genome " 

if ( -e tmp/TSNPH/$zone/$MAGIC.differential_snps ) then
  bin/snp -i tmp/TSNPH/$zone/$MAGIC.snp.sorted$snp_type -db tmp/TSNP_DB/$zone -project $MAGIC -db_frequencyTable -db_frequencyHisto -dropMonomodal -differential -o  `pwd`/tmp/TSNP_DB/$zone/$MAGIC.differential -Reference_genome  "$Reference_genome "
  bin/snp -i tmp/TSNP_DB/$zone/$MAGIC.differential.snp_list_with_allele_frequency_per_sample.txt -db tmp/TSNP_DB/$zone -project $MAGIC -db_prevalenceTable -dropMonomodal -o  `pwd`/tmp/TSNP_DB/$zone/$MAGIC.differential -Reference_genome   "$Reference_genome" 
endif

mv  tmp/TSNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt     tmp/TSNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt.1
cat tmp/TSNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt.1 | gawk '/^#/{print;next;}' > tmp/TSNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt
cat tmp/TSNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt.1 | gawk '/^#/{next;}{print;next;}' >> tmp/TSNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt
\rm tmp/TSNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt.1

# extract the list of genomic SNPs so i can reinject them in the exome study
cat   tmp/TSNP_DB/$zone/$MAGIC.snp_list_with_allele_frequency_per_sample.txt | gawk -F '\t' '/^#/{next}{chrom=$1;pos=$2;typ=$4; gsub(">","2",typ);printf("%s:%d_%s\n",chrom,pos,typ);}' | sort -u >  tmp/TSNP_DB/$zone/$MAGIC.genomic.snp.list

#### snp per gene should we do this here or after filtering
## by the way do we filter here or after gathering all zones 
#echo "## For each gene, we count all SNPs affecting this gene and report the number of observations as wild-type, low, mid, high, pure and the number of samples with at least 0 to 10 doses of any of these SNPs" >> $toto
# cat tmp/TSNP_DB/zone*/$MAGIC.snp_prevalence_per_gene.txt 

exit 0

cat tmp/METADATA/RefSeq.MRNA.ln.ace | gawk '/^Transcript/{gsub(/\"/,"",$2);g=$2;next;}/^IntMap/{m=$2;a1=$3;a2=$4;if(a1<a2){ln=a2-a1+1;}else{ln=a1-a2+1;};printf("Product %s\nBest_product\n\nmRNA %s\nGene %s\nProduct %s 1 %d\nSplicing 1 %d 1 %d Exon 1\nDNA mRNA:%s %d\n",g,g,g,g,ln,ln,ln,g,ln);print;printf("\n");}' >  tmp/TSNP_DB/RefSeq.ace
gunzip -c  TARGET/Targets/ec.RefSeq.fasta.gz  | gawk '/^>/{split($1,aa,"|");m=substr(aa[1],2);printf(">mRNA:%s\n",m);next;}{print}' >  tmp/TSNP_DB/RefSeq.fasta

# analyse the concensus around a particular error : T>G is abundant behing GGt, the worst beeing GGTGG then GGTGn in E.coli or human exome or human RNA-seq, but disappears in LifeTech SOLiD
cat RESULTS/SNV/*.snp_list_with_allele_frequency_per_sample.txt | gawk -F '\t' '{if($28+$29+$30>20 && $11 == motif) {i=index($6,"t");print substr($6,i-2,5);}}' motif="T>G" | tags | sort -k 2n | tail -10


#########################################
#########################################
## identify the rsnumbers
## this script recovers a 3 column table : chrom  coord rsNumber for each chromosome in  tmp/TSNP_RS
# info is in /am/ftp-snp/00README
if (! -d  tmp/TSNP_RS) mkdir tmp/TSNP_RS
set rsdir=/am/ftp-snp/organisms/human_9606_b141_GRCh37p13/chr_rpts/
set ref=GRCh37.p13
if (! -d  tmp/TSNP_RS/$ref) mkdir tmp/TSNP_RS/$ref
foreach chrom ($chromSetAll)
  if (-e tmp/TSNP_RS/pos2rs.$chrom.txt2) continue
  set ff=$rsdir/chr_$chrom.txt.gz
  gunzip -c $ff | gawk -F '\t' '{rs=$1;c=$7;x=$12;withdrawn=$3;suspect=$23;clinvar=$24; if($22 == ref && c==chrom && x != " " && withdrawn == 0) printf ("%s\t%s\t%s\t%s\t%s\n",c,x,rs,suspect,clinvar);}' ref=$ref chrom=$chrom | sort -u | sort -k 1,1 -k2,2n > tmp/TSNP_RS/$ref/pos2rs.$chrom.txt &
end

#########################################
#########################################

if (0) bin/snp -pheno -db tmp/TSNP_DB/$zone

laba2:
bin/tacembly tmp/TSNP_DB/$zone <<EOF
  query find variant
  show -a -f  tmp/TSNP_DB/$zone/snp.translated.ace
  list -a -f tmp/TSNP_DB/$zone/v.list
  query find variant mm
  table -active -o tmp/TSNP_DB/$zone/s20.mm_snps.txt -f tables/s20.mm_snps.def
  query find variant FQ
  list -a -f tmp/TSNP_DB/$zone/vfq.list
  query find variant mrna 
  list -a -f tmp/TSNP_DB/$zone/vm.list
  query find variant translate
  list -a -f tmp/TSNP_DB/$zone/vt.list
  query find Variant coding
  table -active -o tmp/TSNP_DB/$zone/snp2coding.txt -f tables/snp2coding.def
  // table -o tmp/TSNP_DB/$zone/snp.genotype.$zone -f tables/snp.genotype.run.def
  // table -o tmp/TSNP_DB/$zone/snp.pos2pheno.$zone -f tables/snp.pos2pheno.def
  query find variant phenotype
  // table -active -o PHENO/chgMaxAggrColl2.chrom.$zone.txt  -f ../PHENO/variant2map2pheno.def CATEG chgMaxAggrColl2
  quit
EOF
cat  tmp/TSNP_DB/$zone/snp2coding.txt  | gawk -F '\t' '{gsub(/\"/,"",$0);snp=$1;if(snp==old)next;old=snp;typ=$2;mrna=$7;if(substr(mrna,1,1)=="_");mrna=substr(mrna,2);x=1+(($8-$9)/3);a=$3;ap=$4;b=$5;bp=$6;printf("%s\t%s:%s:%s:%s:%d:%s:%s\n",snp,mrna,typ,a,ap,x,b,bp);}' > tmp/TSNP_DB/$zone/snp2coding.nickname.txt

touch tmp/TSNP_DB/$zone/snp.parse.done


echo -n "snp.done"
date

exit 0

set toto=RESULTS/SNV/$MAGIC.SNP_list_with_translation.txt
echo -n "# File $toto : " > $toto
date >> $toto
echo "# Only read pairs aligning uniquely, over at least 90% of their length and as compatible pairs (in close proximity and facing each other) are considered." >> $toto
echo "# To call a SNP or a 1 to 3 bases insertion or deletion, the variant has to be seen in at least one sample, where the variant allele frequency is above $minSnpFrequency %, the coverage above $minSnpCover fold and the variant has at least $minSnpCount supporting reads." >> $toto


cat MetaDB/$MAGIC/RunListSorted  MetaDB/$MAGIC/GroupSnpAdditiveList ZZZZZ MetaDB/$MAGIC/gtitle.txt ZZZZZ  tmp/TSNP_DB/*/snp.translated.ace | gawk -f scripts/snp.translated.awk minSnpCover=$minSnpCover  >  $toto.1
cat $toto.1 | head -20 | gawk '/^#/{print}' >> $toto
cat  $toto.1 | gawk '/^#/{next;}{print}' | sort -k 14,14nr -k 15,15nr -k 16,16nr | sed -e 's/UUUU//g'  >> $toto
\rm $toto.1

###################
## SNP_WIGGLE tracks for UCSC

# 7 Synonymous
# 8 AA_substitution
# 9 Length_variation
# 10 Frameshift
# 11 UTR_5prime
# 12 UTR_3prime
# 13 Non_coding_transcript
#   ......   Kent has splicing i.e. near the intron boundary

if (! -d tmp/TSNPW) mkdir tmp/TSNPW
foreach group (`cat MetaDB/$MAGIC/GroupSnpList | gawk '{print substr($1,1,3);}' | sort -u`)
  if (! -e tmp/TSNPW/$group.bed ) then
    cat  tmp/TSNP_DB/*/s20.mm_snps.txt | gawk -F '\t' '{gsub(/\"/,"",$0);snp=$1;g=$2;z=$3;c=$4;a1=$5;a2=$6;if(a1>a2){a0=a1;a1=a2;a2=a0};a1--;a2--;if (g==group)printf("chr%s\t%d\t%d\t%s\n",c,a1,a2,z);}' group=$group | sort -k 1,1 -k 2,2n >  tmp/TSNPW/$group.bed &
  endif
end
foreach group (`cat MetaDB/$MAGIC/GroupSnpList | gawk '{print substr($1,1,3);}' | sort -u`)
  if (-e  tmp/TSNPW/$group.bed     && ! -e  tmp/TSNPW/$group.bb    ) then
    bedToBigBed  -type=bed4 -tab tmp/TSNPW/$group.bed  tmp/WIGGLE/$species.chrom.sizes  tmp/TSNPW/$group.bb &
  endif
end

\rm toto
foreach group (`cat MetaDB/$MAGIC/GroupSnpList | gawk '{print substr($1,1,3);}' | sort -u`)
  cat <<EOF >> toto
track $group'_SNP'
type bigBed
visibility dense
parent Primates_SNPS
maxHeightPixels 500:150:1
bigDataUrl ftp://ftp.ncbi.nlm.nih.gov/repository/acedb/human/bigwig.hg19/Tissue/$group.snp.bb
priority 4

EOF

end


###################
#########################################
#########################################

