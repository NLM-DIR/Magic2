#!bin/tcsh -f

set Strategy=$1
set zone=$2

  if (! -d tmp/TSNP_DB/$zone/database) then
    pushd tmp/TSNP_DB/$zone
        mkdir database
        ln -s ../../../MetaDB/wspec wspec
        ln -s ../../../metaData/tables.VariantDB tables
        whoami >> wspec/passwd.wrm
        echo y | ../../../bin/tacembly .
    popd

    if (-e TARGET/CHROMS/Centromere_telomeres.ace) then
      echo "pparse TARGET/CHROMS/Centromere_telomeres.ace" | bin/tacembly tmp/TSNP_DB/$zone -no_prompt
    endif

    if (! -e TARGET/CHROMS/Centromere_telomeres.ace && -e TARGET/CHROMS/GCF_report.txt) then
       cat GCF_000001405.22_GRCh37.p10_assembly_report.txt | grep assembled-molecule | cut -f 1,7 | sed -e 's/MT/mito/' | gawk '{printf("Map %s\nTitle %s\n\n",$1,$2);}' >  Chromosome_official_name.ace
       # find Map ; list -a -f titi.list
       cat titi.list | grep NT_ | gawk '/NT_/{print ;split($2,aa,"|");printf("Title \"%s\n\n",aa[2]);}' > titi2
       cat titi2 >> Chromosome_official_name.ace

    endif

    if (-e TARGET/CHROMS/Chromosome_official_name.ace) then
      # expecting  Map chr11 \n  Title NC_0012345.11 \n\n
      # this title will be used in the VCF and the HGVS Human Genome Variation Society full name
      echo "pparse TARGET/CHROMS/Chromosome_official_name.ace" | bin/tacembly tmp/TSNP_DB/$zone -no_prompt
    endif
  endif
  echo "read-models" > tmp/TSNP_DB/$zone/_rs3.$zone

   
  foreach target (mito SpikeIn $Etargets $Ttargets)
    if (-e  tmp/METADATA/$target.MRNA.info.ace) then
      # gene gene_type geneId mrna mrna-length GC
      echo "pparse  tmp/METADATA/$target.GENE.info.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
      echo "pparse  tmp/METADATA/$target.MRNA.info.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
      echo "pparse  tmp/METADATA/$target.MRNA.splicing.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
    endif 
    if (-e tmp/METADATA/gtf.$target.goodProduct.ace) then
      # mrna -> product -> good/best product 
      echo "pparse   tmp/METADATA/gtf.$target.goodProduct.ace"  >> tmp/TSNP_DB/$zone/_rs3.$zone
    endif
    if (-e TARGET/GTF/$species.magic.good_product.ace.gz) then
      # mrna -> product -> good/best product 
      echo "pparse TARGET/GTF/$species.magic.good_product.ace.gz"  >> tmp/TSNP_DB/$zone/_rs3.$zone
    endif
    if (-e  tmp/TSNP_DB/$species.$target.mrna.fasta.gz) then
      echo "pparse tmp/TSNP_DB/$species.$target.mrna.fasta.gz"  >> tmp/TSNP_DB/$zone/_rs3.$zone
    endif 
   end

  # gene mrna product mrna-splicing intmap
  echo "parse TARGET/MRNAS/mrnaStructure.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
  echo "pparse tmp/TSNP_DB/RefSeq.cds.ace.gz" >> tmp/TSNP_DB/$zone/_rs3.$zone
  # product quality should be in the gff file
  echo "pparse  TARGET/MRNAS/very_good_product.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
  echo "pparse  TARGET/MRNAS/good_product.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
  # mRNA-dna extracted from TARGET/Targets
  echo "pparse tmp/TSNP_DB/$species.av.mrna.fasta.gz" >> tmp/TSNP_DB/$zone/_rs3.$zone
  if (-e  tmp/TSNP_ZONE/$zone.fasta.gz && $Strategy != RNA_seq)   echo "pparse tmp/TSNP_ZONE/$zone.fasta.gz"  >> tmp/TSNP_DB/$zone/_rs3.$zone
  # Gene titles and intmap
  echo "pparse  tmp/METADATA/$target.GENE.info.ace" >> tmp/TSNP_DB/$zone/_rs3.$zone
  echo "pparse tmp/SNP_ZONE/$zone.fasta.gz"  >> tmp/TSNP_DB/$zone/_rs3.$zone
  echo save >> tmp/TSNP_DB/$zone/_rs3.$zone
  echo quit >> tmp/TSNP_DB/$zone/_rs3.$zone

  echo -n "snp.parseChrom "
  date
  bin/tacembly tmp/TSNP_DB/$zone < tmp/TSNP_DB/$zone/_rs3.$zone

  touch  tmp/TSNP_DB/$zone/snp.parseChrom.done
  touch  tmp/TSNP_DB/$zone/tsnp0.done

endif

exit 0


if (! -e   tmp/SNP_DB/$zone/snp.Genome_51mer_repeats.ace && -e TARGET/Genome_repeats/Genome_51mer_repeats.gz) then

   bin/tacembly tmp/SNP_DB/$zone <<EOF
     q -o tmp/SNP_DB/$zone/snp_positions.txt select snp,chrom,pos from snp in ?variant, chrom in snp->IntMap, pos in chrom[1] where pos order_by +2+3
EOF
   gzip tmp/SNP_DB/$zone/snp_positions.txt
   date

   gunzip -c tmp/SNP_DB/$zone/snp_positions.txt.gz ZZZZZ.gz TARGET/Genome_repeats/Genome_51mer_repeats.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){z=$2 "\t" $3;chrom[$2]=1;snp[z]=$1;next;}}{if(chrom[$1]<1)next;for (i=$3;i<=$4;i++){z=$1 "\t" i;s=snp[z];if(s)printf("Variant %s\nGenome_51mer_repeats %d\n\n",s,$5);}}' >  tmp/SNP_DB/$zone/snp.Genome_51mer_repeats.ace
   gunzip -c tmp/SNP_DB/$zone/snp_positions.txt.gz ZZZZZ.gz TARGET/Genome_repeats/Genome_101mer_repeats.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){z=$2 "\t" $3;chrom[$2]=1;snp[z]=$1;next;}}{if(chrom[$1]<1)next;for (i=$3;i<=$4;i++){z=$1 "\t" i;s=snp[z];if(s)printf("Variant %s\nGenome_101mer_repeats %d\n\n",s,$5);}}' >  tmp/SNP_DB/$zone/snp.Genome_101mer_repeats.ace
 
if (0) then
  gunzip -c ALL_SNP/snp_positions.txt ZZZZZ.gz TARGET/Genome_repeats/Genome_51mer_repeats.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){z=$2 "\t" $3;chrom[$2]=1;snp[z]=$1;next;}}{if(chrom[$1]<1)next;for (i=$3;i<=$4;i++){z=$1 "\t" i;s=snp[z];if(s)printf("Variant %s\nGenome_51mer_repeats %d\n\n",s,$5);}}' > ALL_SNP/snp.Genome_51mer_repeats.ace
  gunzip -c ALL_SNP/snp_positions.txt ZZZZZ.gz TARGET/Genome_repeats/Genome_101mer_repeats.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){z=$2 "\t" $3;chrom[$2]=1;snp[z]=$1;next;}}{if(chrom[$1]<1)next;for (i=$3;i<=$4;i++){z=$1 "\t" i;s=snp[z];if(s)printf("Variant %s\nGenome_101mer_repeats %d\n\n",s,$5);}}' > ALL_SNP/snp.Genome_101mer_repeats.ace
endif

   date
   bin/tacembly tmp/SNP_DB/$zone <<EOF
     read-models
     pparse tmp/SNP_DB/$zone/snp.Genome_51mer_repeats.ace
     save
     quit
EOF
  date
  
endif

if (! -e  tmp/SNP_DB/$zone/snp.differential.done2 && -e tmp/SNPH/$zone/$MAGIC.differential_snps ) then
  
  cat  tmp/SNPH/$zone/$MAGIC.differential_snps | gawk -F '\t' '/^#/{next}{split($1,aa,":"); f1=$5;f2=$12;if (f1<f2){if(f1<1)typ="Gained_in"; else typ="Up_in";}else{if(f2<1)typ="Lost_in";else typ="Down_in";}printf ("Variant \"%s:%s:%s\"\nDifferential %s \"%s\" \"%s\" %s %s\n\n",aa[1],aa[2],aa[3], typ, $9, $2, f1, f2);}' | grep -v "D_HR_surv_noMYCNA_NB399_exome_Normal" | grep -v Ghs3987 | grep -v Ghs4007 >  tmp/SNP_DB/$zone/snp.differential.ace
  bin/tacembly tmp/SNP_DB/$zone <<EOF
   query find variant Differential
   edit -D Differential
   read-models
   pparse  tmp/SNP_DB/$zone/snp.differential.ace
   query find variant differential ; >Gene
   list -a -f tmp/SNP_DB/$zone/$MAGIC.differential_snp2gene.list
   query find variant Gained_in ; >Gene
   list -a -f tmp/SNP_DB/$zone/$MAGIC.differential_snp2gene.Gained.list
   query find variant Up_in  ; >Gene
   list -a -f tmp/SNP_DB/$zone/$MAGIC.differential_snp2gene.Up.list
   query find variant Down_in ; >Gene
   list -a -f tmp/SNP_DB/$zone/$MAGIC.differential_snp2gene.Down.list
   query find variant Lost_in  ; >Gene
   list -a -f tmp/SNP_DB/$zone/$MAGIC.differential_snp2gene.Lost.list
   save
   quit
EOF
  touch  tmp/SNP_DB/$zone/snp.differential.done
endif

if (-e RESULTS/StLouisSNP/zone.stl.ace) then
  echo "pparse  RESULTS/StLouisSNP/zone.stl.ace"  | bin/tacembly tmp/SNP_DB/$zone -no_prompt 
endif

#########################################
## identify the rsnumbers
## this script recovers a 3 column table : chrom  coord rsNumber for each chromosome in  tmp/SNP_RS
# info is in /am/ftp-snp/00README
if (! -d  tmp/SNP_RS) mkdir tmp/SNP_RS
set rsdir=/am/ftp-snp/organisms/human_9606_b141_GRCh37p13/chr_rpts/
set ref=GRCh37.p13
if (! -d  tmp/SNP_RS/$ref) mkdir tmp/SNP_RS/$ref
foreach chrom ($chromSetAll)
  if (-e tmp/SNP_RS/pos2rs.$chrom.txt2) continue
  set ff=$rsdir/chr_$chrom.txt.gz
  gunzip -c $ff | gawk -F '\t' '{rs=$1;c=$7;x=$12;withdrawn=$3;suspect=$23;clinvar=$24; if($22 == ref && c==chrom && x != " " && withdrawn == 0) printf ("%s\t%s\t%s\t%s\t%s\n",c,x,rs,suspect,clinvar);}' ref=$ref chrom=$chrom | sort -u | sort -k 1,1 -k2,2n > tmp/SNP_RS/$ref/pos2rs.$chrom.txt &
end


###################
#########################################
#########################################
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

if (! -d tmp/SNPW) mkdir tmp/SNPW
foreach group (`cat MetaDB/$MAGIC/GroupSnpList | gawk '{print substr($1,1,3);}' | sort -u`)
  if (! -e tmp/SNPW/$group.bed ) then
    cat  tmp/SNP_DB/*/s20.mm_snps.txt | gawk -F '\t' '{gsub(/\"/,"",$0);snp=$1;g=$2;z=$3;c=$4;a1=$5;a2=$6;if(a1>a2){a0=a1;a1=a2;a2=a0};a1--;a2--;if (g==group)printf("chr%s\t%d\t%d\t%s\n",c,a1,a2,z);}' group=$group | sort -k 1,1 -k 2,2n >  tmp/SNPW/$group.bed &
  endif
end
foreach group (`cat MetaDB/$MAGIC/GroupSnpList | gawk '{print substr($1,1,3);}' | sort -u`)
  if (-e  tmp/SNPW/$group.bed     && ! -e  tmp/SNPW/$group.bb    ) then
    bedToBigBed  -type=bed4 -tab tmp/SNPW/$group.bed  tmp/WIGGLE/$species.chrom.sizes  tmp/SNPW/$group.bb &
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
