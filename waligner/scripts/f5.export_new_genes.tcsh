#!bin/tcsh -f

# bql to export mrna shadows
# select -o /home/mieg/toto.unc116.mrna.preshadow m,type,c,a1,a2,g from m in @, c in m->intmap, m1 in c[1], m2 in c[2], s1 in m->splicing, s2 in s1[1], type in s1[4], g in m->from_gene, a1=m1-s1+1,a2=m1-s2+1
# cat  /home/mieg/toto.unc116.mrna.preshadow | gawk -F '\t' '/Exon/{if($1!=old)n=0;n++;old=$1;printf("%s\t%d\t%s\t%s\t%s\t%s\n",$1,n,$3,$4,$5,$6);}' > /home/mieg/toto.unc116.mrna.shadow
# cat  /home/mieg/toto.unc116.mrna.preshadow | gawk -F '\t' '/Intron/{printf("%s\t%d\t%d\n",$3,$4,$5);}' > /home/mieg/toto.unc116.introns
# select -o /home/mieg/toto.unc32.mrna.preshadow m,type,c,a1,a2,g from m in @, c in m->intmap, m1 in c[1], m2 in c[2], s1 in m->splicing, s2 in s1[1], type in s1[4], g in m->from_gene, a1=m1+s1-1,a2=m1+s2-1
# cat  /home/mieg/toto.unc32.mrna.preshadow | gawk -F '\t' '/Exon/{if($1!=old)n=0;n++;old=$1;printf("%s\t%d\t%s\t%s\t%s\t%s\n",$1,n,$3,$4,$5,$6);}' > /home/mieg/toto.unc32.mrna.shadow
# cat  /home/mieg/toto.unc32.mrna.preshadow | gawk -F '\t' '/Intron/{printf("%s\t%d\t%d\n",$3,$4,$5);}' > /home/mieg/toto.unc32.introns

set chrom=$1

set target=magic
source scripts/target2target_class.txt
set ici=`pwd`
if (! -d  tmp/X.$MAGIC/$chrom) exit 1



echo -n 'Start of phase f5 create the gene boxes '
date

printf "Clone Main_clone\nAceview_release_suffix March13\nname_by_gene\n\n" >  tmp/X.$MAGIC/$chrom/f5.main_clone.ace

bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
      parse tmp/X.$MAGIC/$chrom/f5.main_clone.ace 
      find tg
      acem 
        gene_intron
        quit
      save
      quit
EOF


bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
     query find tg to_be_fused_with
      edit -D to_be_fused_with
      query find tg Antisens_to
      edit -D Antisens_to
      query find tg Overlaps
      edit -D Overlaps
      query find tg Shedded_from
      edit -D Shedded_from
      query find sequence Matching_transcribed_gene
      edit -D Matching_transcribed_gene
      query find sequence Matching_mRNA  
      edit -D  Matching_mRNA  
      save
      quit
EOF

bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
  acembly
    cdna_21  // it seems we did not create identical_to tags
    cdna_50
    quit
  save
  quit
EOF


# find the relevant geneid and export the coords of the future genebox

bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
      table -o tmp/X.$MAGIC/$chrom/f5.tg2genebox.txt -f tables/f5.tg2genebox.def // exclude shedded tg
      save
      quit
EOF

# construct the geneboxes
  cat tmp/X.$MAGIC/$chrom/f5.tg2genebox.txt | gawk -F '\t' '/^\"/{ printf ("Sequence %s\nGenes %s %s %s\n\nGene %s\nTranscribed_gene %s\n\n", $4, $1, $2, $3, $1, $1) ;}'  >  tmp/X.$MAGIC/$chrom/f5.tg2genebox.ace

# parse the genebox and the geneid objects
# flag the shedded genes 
# and reexport the product <->genebox assiociation
bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
      read-models
      pparse tmp/X.$MAGIC/$chrom/f5.tg2genebox.ace
      table -o tmp/X.$MAGIC/$chrom/f5.tgshedded2tg2gg.txt  -f tables/f5.tgshedded2tg2gg.def
      save
      quit
EOF

# hook the shedded tg to the parent genebox
touch  tmp/X.$MAGIC/$chrom/f5.tgshedded2tg2gg.txt
cat tmp/X.$MAGIC/$chrom/f5.tgshedded2tg2gg.txt | gawk -F '\t' '/^\"/{printf("Gene %s\nTranscribed_gene %s\n\n", $3, $1);}' d8.tgshedded2tg2gg.txt >!  tmp/X.$MAGIC/$chrom/f5.tgshedded2tg2gg.ace

# table export the read/tg/mrna/product tags that should appear in the genebox
bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
      query find tg intron ; > intron
      list -a -f  tmp/X.$MAGIC/$chrom/f5.tg2intron.list
      pparse tmp/X.$MAGIC/$chrom/f5.tgshedded2tg2gg.ace // hook shedded tg to parent genebox
      table -o tmp/X.$MAGIC/$chrom/f5.gene2shedtg2coords.txt -f tables/f5.gene2shedtg2coords.def
      table -o tmp/X.$MAGIC/$chrom/f5.geneid2seq2tg2gg.txt   -f tables/f5.geneid2seq2tg2gg.def 
      table -o tmp/X.$MAGIC/$chrom/f5.tg2genebox2product.txt -f tables/f5.tg2genebox2product.def
      table -o tmp/X.$MAGIC/$chrom/f5.pg2tg2gg.txt           -f tables/f5.pg2tg2gg.def
      save
      quit
EOF
 
cat  tmp/X.$MAGIC/$chrom/f5.tg2intron.list | gawk '/^Intron/{gsub(/\"/,"",$2);split($2,aa,"_");printf("%s\t%09d\t%09d\n",aa[1],aa[3],aa[4]);}' >  tmp/X.$MAGIC/$chrom/f5.tg2intron
cat  tmp/X.$MAGIC/$chrom/f5.gene2shedtg2coords.txt |   gawk -F '\t' '/^\"/{g=$1;c=$3;a1=$4;a2=$5;if(g!=old){if(old) printf (" %d %d\n\n", b1,b2);printf("Sequence %s\nGenes %s ",c,g);old=g;b1=a1;b2=a2;}if(a1<a2){if(a1<b1)b1=a1;if(a2>b2)b2=a2;}else{if(a1>b1)b1=a1;if(a2<b2)b2=a2;}}END{if(g) printf (" %d %d\n\n", b1,b2);}'  >! tmp/X.$MAGIC/$chrom/f5.gene2shedtg2coords.ace

echo "pparse  tmp/X.$MAGIC/$chrom/f5.gene2shedtg2coords.ace" | bin/tacembly tmp/X.$MAGIC/$chrom -no_prompt

# Usual coordinate assignments
    echo "gene2chrom"
    bin/gene2chrom2 -any -gg -i  tmp/X.$MAGIC/$chrom  >!  tmp/X.$MAGIC/$chrom/f5.g2c.ggi.ace

# create the geneId and product and genefinder tags tag in the relevant geneboxes

cat  tmp/X.$MAGIC/$chrom/f5.geneid2seq2tg2gg.txt | gawk -F '\t' '/\"/{printf ("GeneId %s\nGene %s\nOther_gene %s\n\n", $1, $2, $2);}' >! tmp/X.$MAGIC/$chrom/f5.geneid2seq2tg2gg.ace
cat tmp/X.$MAGIC/$chrom/f5.tg2genebox2product.txt | gawk  -F '\t' '/\"/ {printf("Gene %s\nProduct %s\n\n", $2, $3); }' >  tmp/X.$MAGIC/$chrom/f5.tg2genebox2product.ace
cat  tmp/X.$MAGIC/$chrom/f5.pg2tg2gg.txt | gawk  -F '\t' '/\"/ {printf("Gene %s\nGeneId %s\nGeneFinder %s\n\n", $4, $3, $1); }' >  tmp/X.$MAGIC/$chrom/f5.pg2tg2gg.ace

bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
      read-models
      query find gene geneid
      edit -D geneid
      pparse tmp/X.$MAGIC/$chrom/f5.pg2tg2gg.ace // parse the contact geneid 
      query find gene geneid
      kstore contacts
      pparse  tmp/X.$MAGIC/$chrom/f5.tg2genebox2product.ace
      pparse  tmp/X.$MAGIC/$chrom/f5.g2c.ggi.ace
      pparse  tmp/X.$MAGIC/$chrom/f5.geneid2seq2tg2gg.ace
      kget contacts
      edit -D geneid // remove the gid from aligned sequence, they survive in Other_geneId
      pparse tmp/X.$MAGIC/$chrom/f5.pg2tg2gg.ace // restore the contact geneid 
      query find geneid IntMap  // should not have been created
      edit -D IntMap
      save
      quit
EOF

  pushd tmp/X.$MAGIC/$chrom 
  ../../../bin/tacembly . <<EOF 
    bql -o f5.g2pg2ll select g,pg,ll from g in ?gene, pg in g->genefinder, ll in pg->LocusLink where ll
    quit
EOF
  cat f5.g2pg2ll | gawk -F '\t' '{printf("Gene %s\nLocusLink %s\n\n",$1,$3);}'  > f5.g2pg2ll.ace
  ../../../bin/tacembly . <<EOF 
    pparse f5.g2pg2ll.ace
    save
    quit
EOF
popd  

## try to rename using the geneboxes 
## recompute since we now have the geneboxes
bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
  query find clone strategy
  edit -D NoKantorInfo
  edit -D name_by_section
  edit name_by_gene
  find tg
  acembly 
    // cdna_71 -split_cloud
    quit
  save
  query find clone strategy
  list -a -f  tmp/X.$MAGIC/$chrom/f5.rename.done
   quit
EOF

## rename
pushd tmp/X.$MAGIC/$chrom
  ../../../bin/tacembly  . <<EOF
    find mrna
    bql -o f5.names select @
    quit
EOF
  cat f5.names | gawk -F '\t' '{m=$1;m0=m;n=length(mm);if (substr(m,1,n)==mm)next;}/^G_t_/{m=substr(m,5);}{printf("-R mRNA %s %s.%s\n\n-R DNA mRNA:%s mRNA:%s.%s\n\n", m0,mm,m,m0,mm,m);}' mm=magic1 > f5.rename.ace
  ../../../bin/tacembly  . <<EOF
    pparse f5.rename.ace
    save
    quit
EOF
popd


## dump : export all including the cloud

if (-e tmp/X.$MAGIC/$chrom/f5.rename.done) then
  bin/gffdump tmp/X.$MAGIC/$chrom >  tmp/X.$MAGIC/$chrom/f5.genes.gtf
  bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
     find mrna
     dna  -noClassName tmp/X.$MAGIC/$chrom/f5.mrnas.fasta
     quit
EOF
  gzip tmp/X.$MAGIC/$chrom/f5.mrnas.fasta tmp/X.$MAGIC/$chrom/f5.genes.gtf

endif

##################################################
## Export the good mRNA for an iteration

if (! -e tmp/X.$MAGIC/$chrom/f5.goodMrnas.fasta.gz) then
  bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
     query find product good_product ; > mrna ; gt_ag
     spush
     query find mrna gt_ag > 1
     sor
     spop
     dna  -noClassName tmp/X.$MAGIC/$chrom/f5.goodMrnas.fasta
     quit
EOF
  cat tmp/X.$MAGIC/$chrom/f5.goodMrnas.fasta | sed -e 's/G_t_/Magic_/g' | gzip >  tmp/X.$MAGIC/$chrom/f5.goodMrnas.fasta.gz
  \rm  tmp/X.$MAGIC/$chrom/f5.goodMrnas.fasta
endif

##################################################


if (0 && ! -e  tmp/X.$MAGIC/$chrom/f5.dump.done) then
  pushd  tmp/X.$MAGIC/$chrom
  if (! -d f5.dumpdir) mkdir  f5.dumpdir
  \rm -rf f5.dumpdir
  mkdir  f5.dumpdir f5.dumpdir/wspec
  cp wspec/* f5.dumpdir/wspec
  $ici/bin/tacembly . <<EOF
      dump -s f5.dumpdir
EOF

  popd
endif 

# export the good products, they are used to beautify the SNP file

set toto=tmp/X.$MAGIC/$chrom
bin/tacembly $toto <<EOF
     bql -o $toto/f5.good_product.txt select m,p,x1,x2,vg from m in ?mRNA where m#from_gene, p in m->product where (p#good_product AND p#best_product) OR p#very_good_product, x1 in p[1], x2 in p[2], vg in p#very_good_product
EOF
wc  $toto/f5.good_product.txt
cat  $toto/f5.good_product.txt  | gawk -F '\t' '{printf("mRNA \"%s\"\nProduct \"%s\" %d %d\n\n", $1, $2, $3, $4);}'  > $toto/f5.good_product.ace
cat  $toto/f5.good_product.txt  | gawk -F '\t' '{printf("Product \"%s\"\nGood_product\n",$2);}/ery_good/{printf("Very_good_product\n");}{printf("\n");}' >> $toto/f5.good_product.ace

gzip  $toto/f5.good_product.ace
\rm  $toto/f5.good_product.txt

##################################################
## done

done:

touch tmp/X.$MAGIC/$chrom/f5.dump.done

echo f5.done
exit 0

bin/clipalign -i Fastc/SRR077419/f.17.fastc.gz -gzo    -t TARGET/Targets/Dmelanogaster.DNASpikeIn.fasta.gz  -maxHit 10 -clipPolyA  -clipPolyT -minEntropy 16 -seedLength 16 -probeMinLength 24  -clipN 2 -minAli 24 -splice         -targetBonus 0     -seedOffset 1 -seedShift 5 -intronMaxLength  100000 -o toto  -showTargetPrefix -target_class 1_DNASpikeIn  -exitAdaptor ACACGCAAACTTCTCAACAGGCCGTACCAATATCCGCAGCTGGTCTCCAAGGTGA,AGGGCAGAGGATGGATGCAAGGATAAGT,AGGTTTGGTCCTAGCCTTTGTATTAGCT,AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT     -showOverhang -strategy RNA_seq

foreach chrom ($chromSetAll)
  end
