#!bin/tcsh -f

set phase=$1

if ($phase == subsampling) goto subsampling
set SNPCHROM=""
if ($phase == ii4) set SNPCHROM=$2

# plase avoid setting mySeaLevel
# This parameter is a way to loop in g4 in debugging mode
# if set, you must run several part of the script by hand
# and there is no guarantee whatsoever on the results
# cat scripts/geneindex.tcsh | gawk -f scripts/csh.beautifier.awk >! scripts/geneindex.beau

if ($phase == ii4) then
  set mySeaLevel=""
  set TGx2=""
else
  set mySeaLevel=$2
  set TGx2=$2
endif

set TGx3=$3
set myRemoveLimit=$3
set mySeaWall=$4
set mySeaUU=$5

setenv GeneIndexDB `pwd`/GeneIndexDB

set refG=""
if ($?Reference_genome) set refG="-referenceGenome $Reference_genome"

if ($phase == TGx) then
  echo phaseTGx.$TGx2.$TGx3
  goto phaseTGx.$TGx2.$TGx3
endif


# if set to N in the directory GENERUNS, the parsing of the lane hits with randomly keep 1/N of the fragments
set SUBSAMPLING=0
if (-e tmp/GENERUNS/SUBSAMPLING) then
  set SUBSAMPLING=`cat tmp/GENERUNS/SUBSAMPLING | gawk '/^SUBSAMPLING/{if ($2>1)n=$2;}END{print n+0;}'`
endif

if ($phase == g1a) goto phaseg1a
if ($phase == g1b) goto phaseg1b
if ($phase == g1c) goto phaseg1c

if ($phase == b0) goto phaseb0

if ($phase == g2a) goto phaseg2a
if ($phase == m2a) goto phasem2a
if ($phase == m2b) goto phasem2b
if ($phase == m2aH) goto phasem2aH
if ($phase == m2bH) goto phasem2bH
if ($phase == g2b) goto phaseg2b
if ($phase == g3a) goto phaseg3a
if ($phase == g3b) goto phaseg3b
if ($phase == m3a) goto phasem3a
if ($phase == m3aH) goto phasem3aH
if ($phase == m3bH) goto phasem3bH
if ($phase == g4) goto phaseg4
if ($phase == g4sp) goto phaseg4sp
if ($phase == g4spx) goto phaseg4spx
if ($phase == gsnp4) goto phasegsnp4
if ($phase == m4) goto phasem4
if ($phase == m4H) goto phasem4H
if ($phase == ma4) goto phasema4
if ($phase == snp4) goto phasesnp4
if ($phase == klst4) goto phaseklst4
if ($phase == r2g4) goto phaser2g4
if ($phase == g5) goto phaseg5  # obsolete
if ($phase == g6) goto phaseg6 # obsolete
if ($phase == g7) goto phaseg7 # obsolete
if ($phase == g8) goto phaseg8 # obsolete

if ($phase == ii1) goto phaseii1
if ($phase == ii2a) goto phaseii2a
if ($phase == ii2b) goto phaseii2b
if ($phase == ii2c) goto phaseii2c
if ($phase == ii2d) goto phaseii2d
if ($phase == ii3a) goto phaseii3a
if ($phase == ii3b) goto phaseii3b
if ($phase == ii4) goto phaseii4
if ($phase == ii5) goto phaseii5
if ($phase == ii6) goto phaseii6
if ($phase == ii99) goto phaseii99

if ($phase =~ SF*) goto phase$phase

echo "g1a:   Create the GeneIndexDB database"
echo "g1b:   Import the metadata of the known genes to the GeneIndexDB database"
echo "g1c:   Exportation des metadata ERCC et MAQC"
echo "b0:    Search for empty lanes, which should not exist"
echo "g2a:   Gene index collect the unique gene expression data for each lane"
echo "m2a:   mRNA index collect the unique gene expression data for each lane"
echo "m2aH:  Hierarchic mRNA index collect the unique gene expression data for each lane"
echo "g2b:   Gene consolidate per run in .u.ace.gz"
echo "m2b:   mRNA consolidate per run in .u.ace.gz"
echo "m2bH:  Hierarchic mRNA consolidate per run in .u.ace.gz"
echo "g3a:   Gene index collect the quasi-unique gene expression data for each lane"
echo "m3a:   mRNA index collect the quasi-unique gene expression data for each lane"
echo "m3aH:  Hierarchic  mRNA index collect the quasi-unique gene expression data for each lane"
echo "g3b:   Gene consolidate per run in .nu.ace.gz"
echo "m3b:   mRNA Consolidate per run in .nu.ace.gz"
echo "m3bH:  Hierarchic  mRNA Consolidate per run in .nu.ace.gz"
echo "g3c:   Not used : Gene consolidate per library in .u. and .nu.ace.gz"
echo "m3c:   Not used : mRNA Consolidate per library in .u. and  .nu.ace.gz"
echo "m3cH:  Not used : Hierarchic  mRNA Consolidate per library in .u. and  .nu.ace.gz"
echo "g4:    Export gene expression tables  in RESULTS/Expression and analyse differential expression according to the Compare class"
echo "g4sp:    Export gene expression tables based on the genebox sponge file"
echo "g4spx:    Export gene expression tables based on the exon/mrna sponge file"
echo "m4:    Export mRNA expression tables  in RESULTS/Expression and analyse differential expression according to the Compare class"
echo "m4H:   Export hierarchic mRNA expression tables  in RESULTS/Expression and analyse differential expression according to the Compare class"
echo "ma4:    Analyse the micro-array data declared in OTHER_PIPELINES as another RESULTS/Expression"
echo "snp4:   Analyse differential SNPs according to the Compare class"
echo "gsnp4:  Export genes differing by the number of protein coding SNPs they contain"
echo "klst4:  Analyse differential mrna using Kallisto.tpm according to the Compare class"
echo "r2g4:   Analyse differential mrna using Kallisto.tpmRNA aditing according to the Compare class"
echo "g5 :    obsolete:  Recalcul des index, par run et par groupe"
echo "g6 :     obsolete : Exportation a partir de GeneIndexDB des tables d'index par genes annotes"
echo "g7 :     obsolete : Export best 1000 genes in each run->Gene_element_index and run->compare_to"
echo "g8 :    Gene expression compare_to runs"

echo "ii2a:   Known_intron counts per lane"
echo "ii2b:   Known_intron counts per run"
echo "ii2c:   Known_intron counts per run with sublibraries"
echo "ii2d:   Histogram of the Cumulated support of annotated exon junctions"
echo "ii3a:   Known_intron counts per lane in non unique alignments"
echo "ii3b:   Known_intron counts per run in non unique alignments"
echo "ii4:    Export introns expression tables  in RESULTS/Expression and analyse differential expression according to the Compare class"
echo "ii5 :   Centralization des elements mrna/gene + gene origin"
echo "ii99 :  Special export to feed the public aceview server"


# 2009_01_28
# CF README.introns for the construction of the intron summary

goto phaseLoop
#######################################################################################
#######################################################################################
# Create the GeneIndexDB database
phaseg1a:

    if (! -d GeneIndexDB) then
      mkdir GeneIndexDB
      pushd GeneIndexDB
        mkdir database
        mkdir wspec
        cd wspec
        foreach ff (`ls ../../metaData/wspec.aceview_web_site/*.wrm`)
          ln -s $ff
        end
        cd ..
        \rm wspec/passwd.wrm
        \rm wspec/layout.wrm
        cp ../metaData/wspec/passwd.wrm wspec
        set mynam=`whoami`
        set n=`cat wspec/passwd.wrm | gawk '{if($1==mynam)n++}END{printf("%d", 0+n)}' mynam=$mynam`
        if ($n == 0) whoami >> wspec/passwd.wrm
        \cp wspec/layout.GeneIndexDB.wrm wspec/layout.wrm
        ln -s ../metaData/tables.GeneIndexDB tables
        echo "initializing GeneIndexDB"
        echo y | ../bin/tacembly .
      popd
    endif

goto phaseLoop

#######################################################################################
#######################################################################################
# Importation of metadata of known genes
phaseg1b:

if (! -d  TARGET || ! -d TARGET/GENES || -e GeneIndexDB/g1b.done) goto phaseLoop
set aceview_reference_db="/net/lmem01/export/home/mieg/37"
set aceview_reference_db="../../../rat_18_1_from_ace11_web_reannotated_dec_2011"
set aceview_reference_db="../../ZZ"

if (! -e  TARGET/GENES/av.tg2clo.ace && -d $aceview_reference_db ) then
  pushd  TARGET/GENES

  $bin/tacembly  $aceview_reference_db  <<EOF
    query find gene ! cloud
    show -a -f av.gene.preace

    query find gene NOT cloud ; >transcribed_gene gene 
    kstore tg1
    query mrna
    bql -a  -o av.tg2g2mrna.txt select tg,g,m from tg in @ , g in tg->gene, m in tg->mrna
    query find predicted_gene Model_of AND NM_id
    bql -a  -o RefSeq.gene_model.txt select RefSeq,model,nm,gid,t,chrom,a1,a2 from RefSeq in @, model in RefSeq->Model_of, nm in RefSeq->NM_id, gid in RefSeq->GeneId_pg, t in gid->title, chrom in RefSeq->intmap, a1 in chrom[1], a2 in a1[1] 
   
    kget tg1
    show -a -f av.tg2clo.preace cdna_clone

    quit
EOF

  cat av.tg2g2mrna.txt | gawk -F '\t' '{printf("Gene %s\nTranscribed_gene %s\nmRNA %s\n\n",$2,$1,$3);}' > av.tg2g2mrna.ace
  cat av.tg2g2mrna.txt | gawk -F '\t' '{printf("Transcribed_gene %s\nmRNA %s\n\n",$1,$3);}' >> av.tg2g2mrna.ace


  cat RefSeq.gene_model.txt | gawk -F '\t' '{printf("Gene %s\nRefSeq\nGene_model %s\n",$2,$1);if($3!="NULL")printf("NM_id %s\n",$3);if($4!="NULL")printf("GeneId %s\n",$4);if($5!="NULL")printf("Title %s\n",$5);if($6!="NULL")printf("\nGenefinder %s\nIntMap %s %d %d\n",$1,$6,$7,$8);printf("\n");}' > RefSeq.gene_model.ace

  cat << EOF >! geneClean.awk
    /^Gene /{printf("Gene %s\nav\n",\$2);next;}
    /^GeneId/{print;next;}
    /^Genefinder/{print;next;}
    /^Transcribed/{print;next;}
    /^Transcribed_gene/{print;next;}
    /^Title/{print;next;}
    /^IntMap/{print;next;}
    /^Main/{print;next;}
    /^Spliced_gene/{print "Main"; next; }
    /^Putative/{print "Main"; next;}
    /^Cloud_gene/{printf("Cloud\n");next;}
    /^Cloud/{print;next;}
    /^\$/{print;next;}
EOF
  gawk -f geneClean.awk av.gene.preace > av.gene.ace
  \rm av.gene.preace
  cat <<EOF >! tg2clo.awk
    /^Transcribed_gene/{if(nnm+nxm+nest>0){gg[gene]=1;gnest[gene]=nest;gnm[gene]=nnm;}tg=\$2;gsub(/\"/,"",tg);gene=tg;j = index(gene, "Apr07") ;  if (j>1) { for (j1 = j-2 ; j1>j-5 && j1>1 ; j1--) if(substr(gene,j1,1)=="."){ gene = substr(gene,1,j1-1); break;}}nnm=0;nxm=0;nest=0;next;}
    /^cDNA_clone/{clo=\$2;if(clo~/XM_/)nxm++;else {if(clo~/NM_/)nnm++;else nest++;}next;}
   END {for (gene in gg){printf("Gene \"%s\"\n",gene);if(gnm[gene]>0)printf("NM\n");if(gnest[gene]>0)printf("cDNA  %d\n",gnest[gene]);printf("\n");}}
EOF
  gawk -f tg2clo.awk av.tg2clo.preace > av.tg2clo.ace
  \rm av.tg2clo.preace


  popd 


endif 


# zcat /home/mieg/37lm5/all.gene.ace.gz | gawk '/^Gene /{printf("\n%s\n",$0);next;}/^[Pp]astille/{print}' > gene2pastille.ace

cat tmp/METADATA/av.f.mrna.sponge | gawk -F '\t' '{a2=$5;if($2==1)m1=$4;if($2>1){m=$1;c=$3;a1=$4;a2=$5;g=$6;if(a1>b2-20)printf("Intron %s__%d_%d\nIntMap %s %d %d\nAV\nFrom_gene %s\nIn_mRNA %s %d %d\nLength %d\n\n",c,b2+1,a1-1,c,b2+1,a1-1,g,m,b2-m1+1,a1-m1-1,a1-b2-1);}b2=a2;}'  >  tmp/METADATA/av.f.mrna2intron.ace
cat tmp/METADATA/av.r.mrna.sponge | gawk -F '\t' '{a2=$5;if($2==1)m1=$4;if($2>1){m=$1;c=$3;a1=$4;a2=$5;g=$6;if(b2>a1+20)printf("Intron %s__%d_%d\nIntMap %s %d %d\nAV\nFrom_gene %s\nIn_mRNA %s %d %d\nLength %d\n\n",c,b2-1,a1+1,c,b2-1,a1+1,g,m,m1-b2+1,m1-a1-1,b2-a1-1);}b2=a2;}'  >  tmp/METADATA/av.r.mrna2intron.ace



if (-e GeneIndexDB/database && ! -e  GeneIndexDB/g1b.done) then
  # Parse all the ace files in TARGET/GENES
  setenv ici `pwd`
  pushd  TARGET/GENES
    $ici/scripts/rr $ici/_r
    cat $ici/_r | sed -e 's/pparse/parse/' >  $ici/_r2
    $ici/bin/tacembly $ici/GeneIndexDB < $ici/_r2
  popd
  # parse all the gene data in METADATA
  set toto=GeneIndexDB/_r.g1b
  echo 'read-models' > $toto
  foreach target ($Etargets)
      echo "parse tmp/METADATA/$target.counts.ace" >> $toto
      echo "parse tmp/METADATA/$target.GENE.info.ace" >> $toto
      echo "parse tmp/METADATA/$target.GENE.ln.ace" >> $toto
      echo "parse tmp/METADATA/$target.MRNA.info.ace" >> $toto
      echo "parse tmp/METADATA/$target.MRNA.ln.ace" >> $toto
      echo "parse tmp/METADATA/gtf.$target.goodProduct.ace" >> $toto
      echo "parse tmp/METADATA/gtf.$target.f.intron.ace" >> $toto
      echo "parse tmp/METADATA/gtf.$target.r.intron.ace" >> $toto
      if (-e tmp/METADATA/$target.split_mrnas.gene2length.ace) then
        echo "parse tmp/METADATA/$target.split_mrnas.gene2length.ace" >> $toto
      endif
  end 
  # echo "parse tmp/METADATA/av.f.mrna2intron.ace" >> $toto
  # echo "parse tmp/METADATA/av.r.mrna2intron.ace" >> $toto
  echo "parse TARGET/MRNAS/good_product.ace" >> $toto
  echo "parse TARGET/MRNAS/very_good_product.ace" >> $toto
  echo "parse TARGET/MRNAS/$species.introns.ace.gz" >> $toto
  echo "query find intron from_gene && ! AV" >> $toto
  echo "edit AV" >> $toto
  echo save >> $toto
  echo quit >> $toto

  bin/tacembly GeneIndexDB < $toto 

  touch GeneIndexDB/g1b.done
endif

goto phaseLoop

########################################################################################
## Importation des donnees ERCC and MAQC dans GeneIndexDB

phaseg1c:

echo "phase geneIndex.a2  Formatage des donnees ERCC"

if (-e GeneIndexDB/g1c.done) goto phaseLoop

if (-e TARGET/MRNAS/Ambion_ERCC_Sequence_cms_095046.txt && ! -e TARGET/MRNAS/Ambion_ERCC_Sequence_cms_095046.ace) then
  pushd TARGET/MRNAS
  cat Ambion_ERCC_Sequence_cms_095046.txt | gawk '/subgroup/{next}{printf("Gene \"%s\"\nMix1 %s %.2f\nMix2 %s\t%.2f\n\n",$2,$4,log($4)/log(2),$5,log($5)/log(2));}' > Ambion_ERCC_Sequence_cms_095046.ace
  bin/dna2dna -getTm -i Ambion_ERCC_Concentration_cms_095047_plus_phiX.fasta.gz |  gawk '/^#/{next}{printf("Gene \"%s\"\nLength %s\nGC_percent %d\n\n",$1,$2,int($5));}' >> Ambion_ERCC_Sequence_cms_095046.ace
  popd
endif

if (-e TARGET/MRNAS/Ambion_ERCC_Sequence_cms_095046.ace) then
  echo "pparse TARGET/MRNAS/Ambion_ERCC_Sequence_cms_095046.ace" | $tacembly GeneIndexDB -noprompt
endif

echo "phase geneIndex.g1c  Exportation des resultats MAQC1"
if ($species != hs) goto phaseLoop

setenv ici `pwd`
if (! -d  TARGET || ! -d TARGET/GENES || -e GeneIndexDB/g1c.done) goto phaseLoop

if (! -d TARGET/MAQC_DATA) mkdir TARGET/MAQC_DATA
pushd  TARGET/MAQC_DATA

if (! -e TARGET/MAQC_DATA/maqc.ace) then
 pushd  TARGET/MAQC_DATA
  $ici/bin/tacembly ~/35f/MaqcInfo/ZZ <<EOF
    query find gene MAQC_acdb_plus
    show -a -f MAQC_acdb_plus.preace MAQC_acdb_plus
    query find gene MAQC_acdb_minus
    show -a -f MAQC_acdb_minus.preace MAQC_acdb_minus

    query find gene  MAQC_avg_signal  UNIQUE Float
    show -a -f MAQC_avg_signal.preace MAQC_signal  
    query find gene MAQC_fold_change UNIQUE Float
    show -a -f MAQC_fold_change.preace MAQC_signal  
          
    query find probe taq*
    query Single_exact_gene  // 1004
    follow Single_exact_gene // 998
    list -a -f geneWithTaq.list
    bql -a -o gene35f2taq2signal.txt select g,p,a,b from g in @, p in g->MAQC_probe where p like "TAQ*", a in p->Signal_A, b  in p->Signal_B
    query find gene MAQC_signal
    query find gene MAQC_avg_signal AND MAQC_fold_change 
    bql -a -o gene35f2maqcSignal2foldChange.txt select g,x,y from g in @, x in g->MAQC_avg_signal, y in g->MAQC_fold_change
    quit
EOF

# use the whole lot, or just the titrating
# cat TAQ2ELEMENT.ace  maqc.A_B.ace  >! MAQC.ace
# cat TAQ2ELEMENT.ace  maqc.titrating.ace  >! MAQC.ace

# rename as elements
  cat MAQC*.preace | gawk '{print}' > maqc.ace
  cat geneWithTaq.list  | gawk '/^KeySet/{next}/^Gene /{printf("Gene  %s\nTAQ\n\n", $2);next;}{print}' > maqc.list.ace
  cat gene35f2maqcSignal2foldChange.txt |  gawk '/^\"/{printf("Gene %s\nMAQC_A %f\nMAQC_B %f\n\n", $1,$2 + $3/2, $2 - $3/2);next;}{print}' > maqc.A_B.ace
  cat gene35f2taq2signal.txt |  gawk '/^\"/{printf("Gene %s\nTAQ_A %f\nTAQ_B %f\n\n", $1,$3, $4);next;}{print}' > TAQ.A_B.ace

# parse and check for elements with tag MAQC and no tag g
  $ici/bin/tacembly $ici/GeneIndexDB << EOF
    // pparse  maqc.ace
    // pparse  maqc.A_B.ace
    // pparse  TAQ.A_B.ace
    query find gene MAQC_acdb_minus OR MAQC_acdb_plus OR TAQ_A OR TAQ_B
    show -a -f maqc.titrating.ace MAQC
    save 
    quit
EOF
 popd
endif # TARGET/MAQC_DATA/maqc.ace

touch GeneIndexDB/g1c.done

goto phaseLoop

#########################################################################################
#########################################################################################
s#########################################################################################
#########################################################################################
## This code is obsolete, because the files were magic SNPs are stored has been modified
## the objective of that codee was to identify the corresponding dbSNP number

# Parsing the SNPdb ASN.1 flat file or the EBI gvf file
# /am/ftp-snp/organisms/human_9606/ASN1_flat/ds_flat_ch18.flat.gz

phasedbsnp2ace:

# read in the tmp/VariantDB so we can translate
if (! -e  tmp/VariantDB/RefSeq2dna.ace) then 
  gunzip -c TARGET/MRNAS/$species.RefSeq.fasta.gz | gawk '/^>/{t=substr($1,2);split(t,aa,"|");printf("\nGene %s\nmRNA %s\nGeneID %s\n\nDNA \"mRNA:%s\"\n",aa[3],aa[1],aa[5],aa[1]);next;}{print}' > tmp/VariantDB/RefSeq2dna.ace
  cat TARGET/MRNAS/RefSeqStructure.txt ZZZZZ | gawk -F '\t' '{m=$1;if(m != oldm){if(oldm && c2>c1)printf("mRNA %s\nCDS %d %d\n\n",oldm,c1,c2);oldm=m;c1=0;c2=0;dx=0;}x1=$5;x2=$6;t=$7;if(t=="\"Exon\""){if(c1==0)c1=dx+1;c2=dx+1+x2-x1;}dx+=x2-x1+1;}' >  tmp/VariantDB/RefSeq2cds.ace
endif

# repris de EBI Homo_sapiens.gvf www.ensembl.org/info/data/ftp/index.html -> easy interface
# 29801527 snp
# 11 dbSNP  SNV 70624  70624  . + . ID=30;Variant_seq=G;Dbxref=dbSNP:rs28602246;Reference_seq=A
# 11 dbSNP  deletion 70625  70625 . +  . ID=31;Variant_seq=-;Dbxref=dbSNP:rs60571627;Reference_seq=G
# gvf release73 ensembl downloaded 2013_10_09

cat ~/RESULTS/EBI/Homo_sapiens.gvf | grep dbSNP | gawk -F '\t' '/dbSNP/{if($3=="SNV")printf ("%d %d %s\n",$4,$5,$7);}' | sort -u | wc

# export from tmp/VariantDB
$tacembly . <<EOF
  find variant
  bql -a -o snp2dbsnp.txt select v,m,x,y from v in class "variant", m in v->intmap, x in m[1], y in m[2]
  query find variant high 
  bql -a -o snp2dbsnp2.txt select v,m,x,y from v in @, m in v->intmap, x in m[1], y in m[2]
  query find variant not_low AND substitution AND IntMap
  bql -a -o snp2dbsnp1.txt select v,m,x,y from v in @, m in v->intmap, x in m[1], y in m[2]
  query find variant high AND substitution AND IntMap
  bql -a -o snp2dbsnp2.txt select v,m,x,y from v in @, m in v->intmap, x in m[1], y in m[2]
  query find variant high AND substitution AND IntMap && NOT low
  bql -a -o snp2dbsnp3.txt select v,m,x,y from v in @, m in v->intmap, x in m[1], y in m[2]
EOF

# rationalize the 2 files
cat ~/RESULTS/EBI/Homo_sapiens.gvf | gawk -F '\t' '/dbSNP/{if($3=="SNV")printf ("%s %09d %s\n",$1,$4,$7);}' | sort -u > snv.pos.ebi.txt
# forget the strand
cat  snv.pos.ebi.txt | gawk '{printf ("%s %09d\n",$1,$2,$3);}' | sort -u  >  snv.pos.ebi.txt2 

# check if the AceView substitution have a dbsnp at the same location
# yes for half of the not_low
cat snp2dbsnp1.txt | grep -v NULL | gawk '/\"/{gsub(/\"/,"",$0);x=$3;y=$4;printf("%s %09d\n",$2,$3);}'  | sort -u > snv.pos.av.high.txt1
wc  snv.pos.ebi.txt2  snv.pos.av.high.txt1
cat  snv.pos.ebi.txt2  snv.pos.av.high.txt1 | wc
cat  snv.pos.ebi.txt2  snv.pos.av.high.txt1 | sort -u | wc

cat snp2dbsnp2.txt | grep -v NULL | gawk '/\"/{gsub(/\"/,"",$0);x=$3;y=$4;printf("%s %09d\n",$2,$3);}'  | sort -u > snv.pos.av.high.txt2
wc  snv.pos.ebi.txt2  snv.pos.av.high.txt2
cat  snv.pos.ebi.txt2  snv.pos.av.high.txt2 | wc
cat  snv.pos.ebi.txt2  snv.pos.av.high.txt2 | sort -u | wc

# yes for 90% of the high AND NOT low
cat snp2dbsnp1.txt | grep -v NULL | gawk '/\"/{gsub(/\"/,"",$0);x=$3;y=$4;printf("%s %09d\n",$2,$3);}'  | sort -u > snv.pos.av.high.txt3
wc  snv.pos.ebi.txt2  snv.pos.av.high.txt3
cat  snv.pos.ebi.txt2  snv.pos.av.high.txt3 | wc
cat  snv.pos.ebi.txt2  snv.pos.av.high.txt3 | sort -u | wc

# Capture the relevant dbsnp
echo ZZZZZ > ZZZZZ 
cat  snp2dbsnp1.txt  ZZZZZ ~/RESULTS/EBI/Homo_sapiens.gvf | gawk '/ZZZZZ/{zz=1;next;}{gsub(/\"/,"",$2);}{if(zz<1){vv[$2 "__" $3]=$1;next;}}/dbSNP/{v=vv[$1 "__" $4]; if(v){i=index($0,"dbSNP:");t=substr($0,i+6);j=index(t,";");printf("dbSNP %s\nVariant %s\nDetails \"%s\"\n\n",substr(t,1,j-1),v,$0);}}' | head

# count in Richa her most homozygote snps
cat /home/agarwala/one_ind.data | cut -f 2,3 | gawk '{if($1+$2>30 && 4*$1 > $1+$2)printf("%d\t%d\n",$1,$1+$2);}' | bin/histo -smooth -plot -o richa.filt5
cat /home/agarwala/one_ind.data | cut -f 2,3 | gawk '{if ($1+$2> 30 &&  100*$1 > 98*($1+$2))printf("%d\t%d\n",$1,$1+$2);}' | wc

# try to rationalize with richa, nothing matches even remotely
cat snp2dbsnp.txt ZZZZZ /home/agarwala/one_ind.data | gawk -F '\t' '/ZZZZZ/{zz=1;next;}{if(zz<1){gsub(/\"/,"",$2);vv[$2 "_" $3]=$1;next;}if ($2+$3> 10 &&  100*$2 > 98*($2+$3)){split($1,aa,":");v=vv[aa[1] "_" aa[2]];printf("%s\t%s\t%s\tv=%s\n",$1,$2,$3,v);}}'
 grep 10042460 snp2dbsnp.txt
 grep 1004246 snp2dbsnp.txt

goto phaseLoop

#########################################################################################
#########################################################################################
# Search for empty lanes, which should not exist" 

phaseb0:

echo "Search for empty lanes in project $MAGIC, which should probably not exist" 
foreach target ($Etargets)
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    ls -ls tmp/GENELANES/$run.*.$target.gene.ace | gawk '{if ($6 == 0) print }'
  end
end
popd

goto phaseLoop

#########################################################################################
#########################################################################################
## centralize the data per lane

phaseg2a:
phaseg3a:
set GM=GENE
set gm=gene
goto phasegm1

# mrna count  
phasem2a:
phasem3a:
set GM=MRNA
set gm=mrna
goto phasegm1

phasem2aH:
phasem3aH:
set GM=MRNAH
set gm=mrna
goto phasegm1

## actual work
phasegm1:

if (! -d tmp/GENERUNS) mkdir tmp/GENERUNS
if (! -d tmp/GENELANES) mkdir tmp/GENELANES

# gene count
#($Etargets mito SpikeIn)
foreach target ($Etargets )
  foreach run (`cat MetaDB/$MAGIC/RunList`)

    if ($target == introns && $GM == MRNA ) continue
    if ($target == mito && $GM == MRNA ) continue
    if ($target == SpikeIn && $GM == MRNA ) continue

    if ($target == introns && $GM == MRNAH ) continue
    if ($target == mito && $GM == MRNAH ) continue
    if ($target == SpikeIn && $GM == MRNAH ) continue

    set  set stranded=0
    foreach run2 (`cat MetaDB/$MAGIC/RunForwardList`)
      if ($run == $run2)  set stranded=1
    end
    foreach run2 (`cat MetaDB/$MAGIC/RunReverseList`)
      if ($run == $run2)  set stranded=-1
    end
    # the experimental value takes precedence */
    set s=`cat MetaDB/$MAGIC/runs.ace | gawk 'BEGIN{x=-1;}/^Run[ \t]/{gsub(/\"/,"",$2);ok=0;if($2==run)ok=1;next;}/^Observed_strandedness_in_ns_mapping/{if(ok==1)x=$3;next;}END{print int(100*x)}' run=$run`
    if ($s >= 0 && $s < 500) set stranded=-1
    if ($s >= 9500 && $s <= 10000) set stranded=1

    set uGeneSupport=""
    set uu=u
    set runs=$run
    if (-e MetaDB/$MAGIC/r2sublib) then
      set runs=`gawk 'BEGIN{r=0;}{if($2==run)r=$1;}END{if(r==0)r=run;print r;}' run=$run MetaDB/$MAGIC/r2sublib`
    endif
 #   echo "SUBLIB $run $runs $target"
    if ($phase == g3a || $phase == m3a || $phase == m3aH) then
      set uGeneSupport=none
      if (-e tmp/GENERUNS/$run/$run.$target.$GM.u."$gm"Support.ace.gz) then
        set uGeneSupport="tmp/GENERUNS/$run/$run.$target.$GM.u."$gm"Support.ace.gz "
      else
        if ($phase == g3a && -e tmp/GENERUNS/$runs/$runs.$target.$GM.u.geneSupport.ace.gz) then
          set uGeneSupport="tmp/GENERUNS/$runs/$runs.$target.$GM.u.geneSupport.ace.gz "
        endif
      endif
     if ($uGeneSupport == none) continue
      set uu=nu
    endif
#echo "AA $run"
    if ($GM == GENE && -e tmp/GENERUNS/$runs/$runs.$target.$uu.geneSupport.ace.gz) continue
#echo "BB $run tmp/GENERUNS/$runs/$runs.$target.$GM.$uu.$gm.Support.ace.gz"
    if (-e tmp/GENERUNS/$runs/$runs.$target.$GM.$uu."$gm"Support.ace.gz) then
      # ls -ls tmp/GENERUNS/$runs/$runs.$target.$GM.$uu."$gm"Support.ace.gz
      continue
    endif
#echo "CC $run"
# ls -ls  tmp/GENERUNS/$run/$run.$target.$GM.$uu."$gm"Support.ace.gz
    #echo "*** $GM $run stranded=$stranded"

    set COUNT=COUNT
    foreach lane (`cat Fastc/$run/LaneList`)
      if (-e tmp/COUNT/$lane.hits.gz) then
        echo "$GM $target $lane"
        if (! -d tmp/GENELANES/$run) source scripts/mkDir GENELANES $run
        if (! -e  tmp/GENELANES/$lane.$target.$GM."$gm"Support.$uu.gz) then
          if (-e  tmp/GENELANES/$lane.$target.$GM.$phase.err) \rm  tmp/GENELANES/$lane.$target.$GM.$phase.*
	  #echo  "scripts/elements.hits2genes.tcsh $GM $target $run  COUNT $lane $stranded $SUBSAMPLING $uGeneSupport "
          scripts/submit tmp/GENELANES/$lane.$target.$GM.$phase "scripts/elements.hits2genes.tcsh $GM $target $run  COUNT $lane $stranded $SUBSAMPLING $uGeneSupport "
        endif
      endif
    end
  end
end

echo "end of loop, no more looping"
goto phaseLoop

###########################################################################
## phase ii1: Import the introns annotated in RefSeq and aceview

phaseii1:
 echo -n "Phase ii1: Import the introns annotated in RefSeq and aceview"

if (-e GeneIndexDB/ii1.done) goto phaseLoop
bin/tacembly GeneIndexDB << EOF
   parse TARGET/MRNAS/$species.introns.ace.gz
   parse TARGET/MRNAS/$species.introns_RNA_seq.ace.gz
   query find transcribed_gene gene AND Intron
   bql -a -o tmp/GENEINDEX/tg2g2ii.txt select tg,g,ii from tg in @, g in tg->gene, ii in tg->intron
   save
   quit
EOF

if (! -e  tmp/GENEINDEX/gene2intron.ace) then
cat  tmp/GENEINDEX/tg2g2ii.txt | cut -f 2,3 | sort -u | gawk  -F '\t' '{g=$1;if(g != old)printf("\n\nGene %s", g);old=g;printf("\nIntron %s", $2);}END{printf("\n\n");}' >  tmp/GENEINDEX/tg2intron.ace
endif

bin/tacembly GeneIndexDB << EOF
   parse tmp/GENEINDEX/tg2intron.ace
   save
   quit
EOF
touch GeneIndexDB/ii1.done
goto phaseLoop

###########################################################################
## phase ii2: Intron support read in the best alignment to the transcriptome

phaseii3a:
  goto phaseLoop
phaseii2a:
 echo -n "Phase ii2a: Intron support read in the best alignemnt to the transcriptome"
 date

if (! -d tmp/INTRON_INDEX) mkdir tmp/INTRON_INDEX
if (! -d tmp/INTRONLANES) mkdir tmp/INTRONLANES

foreach run (`cat MetaDB/$MAGIC/RunList`)
  if ( -e tmp/INTRONRUNS/$run/$run.u.intronSupport.ace.gz) continue
  if (! -d tmp/INTRONLANES/$run) source scripts/mkDir INTRONLANES $run
 
  set pair=0
  foreach run2 (`cat MetaDB/$MAGIC/RunPairedList`)
    if ($run == $run2) then
      set pair="500"
    endif
  end

  foreach lane (`cat Fastc/$run/LaneList`)

    if (-e tmp/METADATA/mrnaRemap.gz && ! -e tmp/INTRONLANES/$lane.intronSupport.u.gz && -e tmp/COUNT/$lane.hits.gz) then
        if (-e tmp/INTRON_INDEX/Exons_juntions.stranding.ace) \rm tmp/INTRON_INDEX/Exons_juntions.stranding.ace
        if (-e tmp/INTRONLANES/$lane.$phase.err) \rm tmp/INTRONLANES/$lane.$phase.*
	echo "scripts/$phase.intron_support.tcsh $run $lane $pair"
        scripts/submit tmp/INTRONLANES/$lane "scripts/ii3a.intron_support.tcsh $run $lane $pair"
    endif
  end  
end

goto phaseLoop

##############################################################################
## centralize the gene unique support data per run as ace files

phaseg2b:
  set uu=u
  set GM=GENE
  set gm=gene
  set support2ace='-geneSupport2ace'
goto phasegm2

phaseg3b:
  set uu=nu
  set GM=GENE
  set gm=gene
  set support2ace='-geneSupport2ace'
goto phasegm2

# mrna count  
phasem2b:
  set uu=u
  set GM=MRNA
  set gm=mrna
  set support2ace='-mrnaSupport2ace'
goto phasegm2

phasem2bH:
  set uu=u
  set GM=MRNAH
  set gm=mrna
  set support2ace='-mrnaSupport2ace'
goto phasegm2

phasem3b:
  set uu=nu
  set GM=MRNA
  set gm=mrna
  set support2ace='-mrnaSupport2ace'
goto phasegm2

phasem3bH:
  set uu=nu
  set GM=MRNAH
  set gm=mrna
  set support2ace='-mrnaSupport2ace'
goto phasegm2


phasegm2:

if (! -d tmp/GENEINDEX) mkdir tmp/GENEINDEX
#  ($Etargets mito SpikeIn)
foreach target ($Etargets)
  echo "target=$target"
  foreach lib (`cat MetaDB/$MAGIC/RunsList`)
    # back compatiblity
    if ($GM == GENE && -e tmp/GENERUNS/$lib/$lib.$target.$uu.geneSupport.ace.gz) then
      mv  tmp/GENERUNS/$lib/$lib.$target.$uu.geneSupport.ace.gz tmp/GENERUNS/$lib/$lib.$target.$GM.$uu.geneSupport.ace.gz
    endif

    if (-e tmp/GENERUNS/$lib/$lib.$target.$GM.$uu."$gm"Support.ace.gz) continue
    if (! -d tmp/GENERUNS/$lib) source scripts/mkDir GENERUNS $lib
    set ok=1

    if (-e  tmp/GENERUNS/$lib/$lib.$target.$GM.$uu.list) mv  tmp/GENERUNS/$lib/$lib.$target.$GM.$uu.list   tmp/GENERUNS/$lib/$lib.$target.$GM.$uu.list.old
    if (-e   tmp/GENERUNS/$lib/$lib.$target.$GM.3pHisto.list) mv   tmp/GENERUNS/$lib/$lib.$target.$GM.3pHisto.list   tmp/GENERUNS/$lib/$lib.$target.$GM.3pHisto.list.old
                                                                  
    foreach run (`cat MetaDB/$MAGIC/r2sublib | gawk -F '\t' '{if($1 == lib)print $2;}END{print lib}' lib=$lib | cut -f 1 | sort -u`)
      if (! -d Fastc/$run) continue
        foreach lane (`cat Fastc/$run/LaneList`)
	  if (-e tmp/GENELANES/$lane.$target.$GM."$gm"Support.$uu.gz) then
            echo  tmp/GENELANES/$lane.$target.$GM."$gm"Support.$uu.gz >> tmp/GENERUNS/$lib/$lib.$target.$GM.$uu.list
          else

            if ($ok == 1) echo "missing file  tmp/GENELANES/$lane.$target.$GM."$gm"Support.$uu.gz"
            set ok=0
          endif

          foreach kb (8kb 5kb)	  
            if (-e tmp/GENELANES/$lane.$target.$GM.3pHisto.$kb.txt) then 
              echo  tmp/GENELANES/$lane.$target.$GM.3pHisto.$kb.txt >> tmp/GENERUNS/$lib/$lib.$target.$GM.3pHisto.$kb.list
            endif
          end
        end
     end
     if ($ok == 0) continue
     
     scripts/submit tmp/GENERUNS/$lib/$lib.$target.$GM.$uu "scripts/g2m.tcsh $support2ace $lib $target $GM $gm $uu $phase"

  end
end

goto phaseLoop

##############################################################################
## centralize the introns unique support data per run as ace files
phaseii2b:
phaseii3b:
  set uu=u
  if ($phase == ii3b) set uu=nu

if (! -d tmp/INTRONRUNS) mkdir tmp/INTRONRUNS
foreach target ($Etargets)
  break
end
  echo "... $phase $target"
  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    if (-e tmp/INTRONRUNS/$run/$run.$uu.intronSupport.ace.gz && -e  tmp/INTRONRUNS/$run/known_introns.$target.ace) continue
    if (! -d tmp/INTRONRUNS/$run) source scripts/mkDir INTRONRUNS $run
    set ok=1
    echo "... $phase $target ok=1"
    if (-e  tmp/INTRONRUNS/$run/$run.$uu.list) \rm  tmp/INTRONRUNS/$run/$run.$uu.list
    foreach lib ($run `cat MetaDB/$MAGIC/r2sublib | gawk '{if(r==$1)print " "$2;}' r=$run`)
      if (! -e Fastc/$lib/LaneList) continue
        foreach lane (`cat Fastc/$lib/LaneList`)
	  if (-e tmp/INTRONLANES/$lane.intronSupport.$uu.gz) then
            echo  tmp/INTRONLANES/$lane.intronSupport.$uu.gz >> tmp/INTRONRUNS/$run/$run.$uu.list
          else
            echo "ERROR phase ii2b MISSING FILE tmp/INTRONLANES/$lane.intronSupport.$uu.gz"
            echo "ERROR phase ii2b MISSING FILE tmp/INTRONLANES/$lane.intronSupport.$uu.gz" >> tmp/INTRONRUNS/$run/$run.$uu.list
            set ok=0
          endif
        end
     end
     echo "... $phase $target ok2=$ok"
     if ($ok == 0) continue
     if (-e tmp/introns/$MAGIC.cumulated_support.txt) \rm  tmp/introns/$MAGIC.cumulated_support.txt
     if (-e tmp/introns/$MAGIC.introns.stranding.ace) \rm  tmp/introns/$MAGIC.introns.stranding.ace 
     if (-e tmp/introns/$MAGIC.intronSupportPerRun.txt) \rm  tmp/introns/$MAGIC.intronSupportPerRun.txt
     scripts/submit tmp/INTRONRUNS/$run/ii2b.run2introns  "scripts/ii2b.run2introns.tcsh $run $uu"
  end


goto phaseLoop

##############################################################################
## centralize the gene support data per project as ace files
phaseg3b_doublon:

foreach target ($Etargets mito SpikeIn)
  echo "target=$target"
    foreach uu (u nu)
      if (! -e tmp/GENEINDEX/$MAGIC.$target.$uu.geneSupport.ace.gz) then
        if (-e GeneIndexDB/parse.index.done) \rm GeneIndexDB/parse.index.done
        if (-e  tmp/GENEINDEX/$MAGIC.$target.$uu.list) \rm  tmp/GENEINDEX/$MAGIC.$target.$uu.list
        foreach lane (`cat MetaDB/$MAGIC/LaneList`)
	  if (-e tmp/GENELANES/$lane.$target.GENE.geneSupport.$uu.gz) then
            echo  tmp/GENELANES/$lane.$target.GENE.geneSupport.$uu.gz >>  tmp/GENEINDEX/$MAGIC.$target.$uu.list
          endif
        end
        scripts/submit tmp/GENEINDEX/$MAGIC.$target.$uu "scripts/g3b.gene_support.tcsh tmp/GENEINDEX/$MAGIC.$target.$uu  tmp/GENEINDEX/$MAGIC.$target.$uu.list"
        if (-e tmp/GENEINDEX/$MAGIC.$target.$uu.ace) then \rm tmp/GENEINDEX/$MAGIC.$target.$uu.ace
      endif
    end
end

goto phaseLoop

# tmp/GENELANES/*/*.$target.GENE.geneSupport.$uu.gz 

#######################################################################################
phaseii2c:

 echo "... ii2c Known_intron counts in Runs with sublibraries"
 foreach run (`cat MetaDB/$MAGIC/r2sublib | cut -f 1 | sort -u`)
    set ok=0

    foreach target ($Etargets)    
       if  (-e tmp/INTRONRUNS/$run/known_introns.$target.ace) continue
       set ok=1
    end
    if ($ok == 1) then
      if (! -d  tmp/INTRONRUNS) mkdir tmp/INTRONRUNS
      if (! -d  tmp/INTRONRUNS/$run) source scripts/mkDir INTRONRUNS $run

      scripts/submit  tmp/INTRONRUNS/$run/ii2c "scripts/ii2c.run2introns.tcsh $run 0"
   endif
end

goto phaseLoop

#######################################################################################
phaseii2d:

 set level=1

 echo "... ii2d Known_intron counts in groups level $level"
 foreach run (`cat MetaDB/$MAGIC/g2r | gawk -F '\t' '{if($3 == level) print $1;}' level=$level | sort -u`)
    set ok=0

    foreach target ($Etargets)    
       if  (-e tmp/INTRONRUNS/$run/known_introns.$target.ace) continue
       set ok=1
    end
    if ($ok == 1) then
      if (! -d  tmp/INTRONRUNS) mkdir tmp/INTRONRUNS
      if (! -d  tmp/INTRONRUNS/$run) source scripts/mkDir INTRONRUNS $run

      scripts/submit  tmp/INTRONRUNS/$run/ii2c "scripts/ii2c.run2introns.tcsh $run $level"
   endif
end

goto phaseLoop

#######################################################################################
phaseg3c:
  set GM=GENE
  set GMS=geneSupport
  goto phasegm3c
phasem3c:
  set GM=MRNA
  set GMS=mrnaSupport
  goto phasegm3c
phasem3cH:
  set GM=MRNAH
  goto phasegm3c

phasegm3c:
 echo "... $phase $GM counts in Runs with sublibraries"
 foreach lib (`cat MetaDB/$MAGIC/r2sublib | cut -f 1 | sort -u`)
   foreach target ($Etargets)    
     foreach uu (u nu)
       set out=tmp/GENERUNS/$lib/$target.$GM.$uu.$GMS
       if (-e $out.ace.gz) continue
       if (-e GeneIndexDB/parse.index.done) \rm GeneIndexDB/parse.index.done
       set inFileList tmp/GENERUNS/$lib/$target.$uu.list 
     if (-e $inFileList) \rm $inFileList

       if (! -d  tmp/GENERUNS/$lib) source scripts/mkDir GENERUNS $lib
       set ok=1
       set ok2=0
       foreach run (`cat MetaDB/$MAGIC/r2sublib | gawk -F '\t' '{if($1 == lib)print $2;}' lib=$lib | cut -f 1 | sort -u`)
         if  (-e tmp/GENERUNS/$run/$target.$GM.$uu.$GMSgeneSupport.ace.gz) then
           set ok2=1
           echo tmp/GENERUNS/$run/$target.$GM.$uu.$GMSgeneSupport.ace.gz >> $inFileList 
         else
           set ok=0
         endif
       end
       if ($ok == 0 || $ok2 == 0) continue

       scripts/submit tmp/GENERUNS/$lib/$phase "scripts/gm3c.runs2libs.tcsh $inFileList $out"
     end
   end
end

goto phaseLoop

#######################################################################################
#######################################################################################
# centralization des elements mrna/gene + gene origin

# AUC parse FDR counts Diff_

phaseii4:


  set chrom=$SNPCHROM
  if ($chrom == "") then
    echo "Missing parameter chrom in call to : genetindex.tcsh ii4 chrom"
    exit 1
  endif

# goto phaseLoop

# g4sp expression based on sponge file is not yet written as of 2019_10_14
phaseg4sp:
phaseg4spx:
phasesnp4:
touch tmp/GENEINDEX/toto.ace
\rm tmp/GENEINDEX/$MAGIC.*.ace

if ($SNPCHROM == 000) set mySeaLevel=""
echo "phase $phase $mySeaLevel"
phaseklst4:
phaser2g4:
phaseg4:
phasem4:
phasem4H:
phasema4:
phaseSF24A:
phaseSF24B:
phaseSF31A:
phaseSF31B:

if ($mySeaLevel != "") goto justRunG4

echo LALALA0========A

bin/bestali -checkGroupHierarchy -project $MAGIC -db MetaDB

echo LALALA0========B

if (! -d RESULTS/Expression)  mkdir RESULTS/Expression
if (! -d RESULTS/Expression/unique) mkdir RESULTS/Expression/unique
if (! -d RESULTS/Expression/quasi_unique) mkdir RESULTS/Expression/quasi_unique
if (! -d RESULTS/Expression/ERCC) mkdir RESULTS/Expression/ERCC
if (! -d RESULTS/Expression/unique/Mitochondria) mkdir RESULTS/Expression/unique/Mitochondria
if (! -d RESULTS/Expression/quasi_unique/Mitochondria) mkdir RESULTS/Expression/quasi_unique/Mitochondria
if (! -d tmp/GENEINDEX/Results) mkdir tmp/GENEINDEX/Results

echo LALALA0========C

# one may wih to use -minIndex 10
# av RefSeq seqc SpikeIn
justRunG4:
echo "justRunG4 $Etargets"
#  av RefSeq seqc SpikeIn magic
if (! $?other_targets) set other_targets=""
if (! $?Etargets) then
  set Etargets="$targets"
endif
# $Etargets introns snp $other_targets
# ($Etargets introns snp $other_targets)
foreach target ($Etargets introns snp $other_targets) 

  set myace=""
  set targetBeau=$target
  if ($target == av) set targetBeau=AceView

  if ($phase == ma4) then
    if (! -e OTHER_PIPELINES/$target.ace) then
      continue
    endif
  endif
  # if ($target != av && $target != RefSeq) continue 

  if ($target == gdecoy) continue
  if ($target == DNASpikeIn ) continue
  if ($target == rrna ) continue
  if ($target == rnaGene ) continue

  if ($phase == ii4 && $target != introns) continue
  if ($phase != ii4 && $target == introns) continue

  if ($phase == snp4 && $target != snp) continue
  if ($phase != snp4 && $target == snp) continue

  if (! -d RESULTS/Expression/unique/$target) mkdir RESULTS/Expression/unique/$target
  if (! -d RESULTS/Expression/quasi_unique/$target) mkdir RESULTS/Expression/quasi_unique/$target

  source scripts/target2target_class.txt

  set ok=0
  foreach target2 ($Etargets toto)
    if ($phase == snp4 || $phase == ii4 || $phase == ma4 || $phase == klst4 || $target == $target2) set ok=1
  end
  if ($ok == 0) continue

  echo -n "phase $phase : geneindex export table and correlation analysis: $target "
  date

echo LALALA0==================

# if ($mySeaLevel == "") 

set ok=1
if (1) then

  set GM=GENE
  if ($phase == ma4 ) set GM=MA
  if ($phase == klst4 ) set GM=Kallisto
  if ($phase == r2g4 ) set GM=r2g
  if ($phase == snp4) set GM=SNP
  if ($phase == g4sp) set GM=GENESP
  if ($phase == g4spx) set GM=GENESPX
  if ($phase == g4spx) set phase=g4sp
  if ($phase == gsnp4) set GM=PROTEIN_CHANGING_SNP_PER_GENE
  if ($phase == m4)   set GM=MRNA
  if ($phase == m4H)  set GM=MRNAH
  if ($phase == ii4)  set GM=INTRON
  if ($phase =~ SF* )  set GM=$phase

   set CAPT=toto
   if (1 && $MAGIC == CL1) set CAPT=A1R3
   if ($MAGIC == CL2) set CAPT=A1
   if ($MAGIC == CL3) set CAPT=A1
   if ($phase == ii4zzzz) set CAPT=A1
   if (0) then
     set CAPT=A1A2I3R1R2   # A2R2 ... see TARGET/GENES/$capture.av.gene_list  
     set CAPT=A1A2I2I3R1R2   # A2R2 ... see TARGET/GENES/$capture.av.gene_list  
     if (! -e TARGET/GENES/$CAPT.capture.$target.gene_list) then
       echo "ERROR: missing file TARGET/GENES/$CAPT.capture.$target.gene_list"
       continue
     endif
    endif
    set capt=""
    if ($CAPT == toto) then
      set CAPT=""
    else
      set capt="-captured $CAPT"
      set CAPT=".$CAPT" 
    endif

echo LALALA0==================---
  if ($?chrom) then
    if ($phase == ii4) then
      bin/tace tmp/INTRON_DB/$chrom <<EOF
        find run
	edit -D project
        pparse MetaDB/$MAGIC/runs.ace
	save
        query find intron de_duo
	bql -o  tmp/INTRON_DB/$chrom/ii4.de_duo.txt select ii,r,n from ii in @, r in ii->de_duo where r->project == "$MAGIC", n in r[1]
EOF
      cat tmp/INTRON_DB/$chrom/ii4.de_duo.txt | gawk -F '\t' '{if($1 != ii){printf("\nIntron %s\n", $1);ii=$1;}printf ("Run_U %s 0 %d seqs %d tags %d kb\n", $2,$3,$3,$3);}' >  tmp/INTRON_DB/$chrom/ii4.index.ace
    endif
  endif

echo LALALA0==================+++
  echo  "phase $phase : geneindex export table and correlation analysis: $target GM=$GM"

# export the runs and groups
  cat MetaDB/$MAGIC/runs.ace  > tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  cat MetaDB/$MAGIC/groups.ace  >> tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  cat MetaDB/$MAGIC/compares.ace  >> tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace

# export the sea level
  set chrom=""
  if ($phase == ii4) set chrom=$SNPCHROM
  if ($phase != ii4 && -e MetaDB/$MAGIC/ali.ace) cat MetaDB/$MAGIC/ali.ace | gawk '/^Ali /{print}/h_Ali/{print}/^Intergenic/{print}/^Accessible_length/{print}/^$/{print}' >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ($phase == ii4 && -e MetaDB/$MAGIC/ali.ace) cat MetaDB/$MAGIC/ali.ace | gawk '/^Ali /{print}/h_Ali/{print}/^Candidate_introns any/{print}/^$/{print}' >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
echo LALALA1====================++99
  if ($phase == ii4 && -e tmp/METADATA/$MAGIC.av.captured_genes.ace) cat tmp/METADATA/$MAGIC.av.captured_genes.ace >> tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ($phase == ii4 && -e tmp/INTRON_DB/$chrom/d5.$MAGIC.info$CAPT.ace) cat tmp/INTRON_DB/$chrom/d5.$MAGIC.info$CAPT.ace >> tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace

echo LALALA1====================

# export the general METADATA.ace
  if (($GM == GENE || $GM == GENESP || $GM == GENESPX || $GM =~ SF* ) && -e tmp/METADATA/$target.GENE.info.ace) cat  tmp/METADATA/$target.GENE.info.ace >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if (($phase == m4 || $phase == m4H) && -e tmp/METADATA/$target.MRNA.info.ace) cat tmp/METADATA/$target.MRNA.info.ace >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace

  if ( -e TARGET/GENES/$target.gene.ace) cat  TARGET/GENES/$target.gene.ace >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace

  if ( $phase == g4 && -e TARGET/GENES/$target.gene2geneid.ace) cat TARGET/GENES/$target.gene2geneid.ace  >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ( $phase == g4 && -e tmp/METADATA/$MAGIC.$target.captured_genes.ace) cat  tmp/METADATA/$MAGIC.$target.captured_genes.ace  >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ( $phase == g4sp && -e TARGET/GENES/$target.gene2geneid.ace) cat TARGET/GENES/$target.gene2geneid.ace  >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ( $phase == g4sp && -e tmp/METADATA/$MAGIC.$target.captured_genes.ace) cat  tmp/METADATA/$MAGIC.$target.captured_genes.ace  >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ( $phase == g4 && -e tmp/METADATA/$target.split_mrnas.gene2length.ace) cat tmp/METADATA/$target.split_mrnas.gene2length.ace >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ( $phase == g4sp && -e tmp/METADATA/$target.split_mrnas.gene2length.ace) cat tmp/METADATA/$target.split_mrnas.gene2length.ace >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace


  # if ( $target == RefSeq && -e TARGET/GENES/RefSeq.gene_model.ace) cat  TARGET/GENES/RefSeq.gene_model.ace >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ($phase == g4 && -e TARGET/GENES/$target.gene2nm_id.ace) cat  TARGET/GENES/$target.gene2nm_id.ace >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ($phase == g4sp && -e TARGET/GENES/$target.gene2nm_id.ace) cat  TARGET/GENES/$target.gene2nm_id.ace >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ($phase == m4H && -e TARGET/MRNAS/$target.transcript2nm_id.ace) cat TARGET/MRNAS/$target.transcript2nm_id.ace >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace


  if ($target == intronsZZZ && -e TARGET/GENES/gene2intron.ace) cat TARGET/GENES/gene2intron.ace >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ($target == intronsZZZ && -e TARGET/MRNAS/$species.introns_RNA_seq.ace.gz) gunzip -c TARGET/MRNAS/$species.introns_RNA_seq.ace.gz >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ($target == intronsZZZ && -e TARGET/MRNAS/$species.introns.ace.gz ) gunzip -c TARGET/MRNAS/$species.introns.ace.gz >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  foreach target2 ($Etargets)
    if ($target == intronsZZZ && -e tmp/METADATA/$target2.introns.info.ace ) cat  tmp/METADATA/$target2.introns.info.ace >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  end

  # if ( -e TARGET/GENES/$species.a5.geneid2nm_id.ace) cat TARGET/GENES/$species.a5.geneid2nm_id.ace >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace

  if ( $target == RefSeq && -e TARGET/Targets/pg2intmap.RefSeq2008.ace) cat TARGET/Targets/pg2intmap.RefSeq2008.ace  >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ($species == rn && -e MicroArray/Hits/Affy.Rat230_2.$target.probeset2gene.quasi_unique.ace) cat  MicroArray/Hits/Affy.Rat230_2.$target.probeset2gene.quasi_unique.ace >> tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ($phase == g4 && ($MAGIC == Kidney || $MAGIC == NB || $MAGIC == NB_DEG || $MAGIC == NB_export) && -e MicroArray/Hits/AGLuK.$target.probe2gene.unique.ace) cat MicroArray/Hits/AGLuK.$target.probe2gene.unique.ace >> tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace

   if ($phase == m4H && $MAGIC == NB && -e MicroArray/Hits/AGLuK.$target.probe2mrna.unique.txt) then
      cat MicroArray/Hits/AGLuK.$target.probe2mrna.unique.txt | gawk -F '\t' '{printf("Transcript %s\nMicroArray %s\n\n",$2,$1);}'  >> tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
    endif

  if ($phase == ma4 && -e MicroArray/Hits/$target.av.probe2mrna.unique.txt) cat MicroArray/Hits/$target.av.probe2mrna.unique.txt | gawk -F '\t' '{printf("Gene %s\nFrom_transcript %s\n\n",$1,$2);}'  >> tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
  if ($phase == ma4 && -e MicroArray/Hits/$target.av.probe2mrna.unique.txt) cat MicroArray/Hits/$target.av.probe2gene.unique.txt | gawk -F '\t' '{printf("Gene %s\nFrom_gene %s\n\n",$1,$2);}'  >> tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace

  if ($phase == ma4 && -e MicroArray/Hits/$target.probe2chrom.unique.txt) cat  MicroArray/Hits/$target.probe2chrom.unique.txt | gawk -F '\t' '{split($2,aa,":");gsub(/chr/,"",aa[1]);split(aa[2],bb,"-");printf("Gene %s\nIntMap %s %s %s\n\n",$1,aa[1],bb[1],bb[2]);}' >> tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace




  if (($MAGIC == SEQC_main || $MAGIC == SEQC) && -e MicroArray/Hits/Affy_U133_plus_2.$target.probeset2gene.unique.ace) cat  MicroArray/Hits/Affy_U133_plus_2.$target.probeset2gene.unique.ace >> tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace

  if (-e RESULTS/Coverage_and_exons/$MAGIC.intergenicKb.ace) cat RESULTS/Coverage_and_exons/$MAGIC.intergenicKb.ace >>  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace

  if (0 && $phase == m4H && -e TARGET/MRNAS/$species.m4H.NMNR2avMRNA.hack) then
    cat  tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace | grep -v NM_id > toto.$$
    cat toto.$$ TARGET/MRNAS/$species.m4H.NMNR2avMRNA.hack >> tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace
    \rm toto.$$
  endif


endif

echo LALALA2

  set sg=""
  if (-e TARGET/Targets/$species.$target.stable_genes.txt) set sg="-stableGenes TARGET/Targets/$species.$target.stable_genes.txt"


if ($SNPCHROM == 000) exit 0
if ($ok == 0) continue

  set method=""
  if ($target == SpikeIn) set method="-method ERCC -malus 6 "
  if ($target == RefSeq) set method="-method Gene_RefSeq"
  if ($target == av) set method="-method Gene_AceView"
  if ($target == av && species == hs) set method="-method Gene_AceView2010"
  if ($target == av && $species == rat) set method="-method Gene_AceView08"

  if ($target == av && $phase == m4) set method="-method AceView_Transcript"
  if ($target == RefSeq && $phase == m4) set method="-method RefSeq_Transcript"
  if ($target == av && $phase == m4H) set method="-method AceView_Hierarchic_Transcript"
  if ($target == RefSeq && $phase == m4H) set method="-method RefSeq_Hierarchic_Transcript"

  if ($phase == ma4) set method="-method $target -MA 60 -keepIndex -minFoldChange .5" # detrimental for NB_AGLuK 2013_11_13
  if ($phase == ma4) set method="-method $target -MA 60 -keepIndex "
  if ($phase == snp4) set method="-method $target  "
  if ($phase == klst4) set method="-method kallisto.$target -keepIndex -threshold 4"
  if ($phase == r2g4) set method="-method r2g.$target -threshold 4"

  if ($phase =~ SF*) set method="-method sailfish$phase "

  set chromAlias=""
  if (-e TARGET/Targets/$species.chromAlias.txt) set chromAlias="-chromAlias TARGET/Targets/$species.chromAlias.txt"
  if (-e tmp/METADATA/$target.split_mrnas.gz) set chromAlias="$chromAlias -split_mRNAs tmp/METADATA/$target.split_mrnas.gz"



# -export abti $correl 
  set seedGene=""
  if ($target == av && $MAGIC == NB) set seedGene="-seedGene KDM5D"
  if ($target == av && $MAGIC == NB) set seedGene="-seedGene MYCN"
  if ($target == RefSeq && $MAGIC == NB) set seedGene="-seedGene X__MYCN"
  if ($target == av && $MAGIC == Klein) set seedGene="-seedGene Aanat"
  set correl="" 
  set correl="-exportDiffGenes "  
  if ($target == av && $MAGIC == Liver) set correl="-TGx "

  set compare="-compare -correlation"
  if ($MAGIC == SEQC2ZZZZZ) set compare=""
  set shA=" "
  if ($MAGIC == Klein) set shA="-showAnyRun"

   set uu=u 
   if ($?myUU) then
     if ( $myUU != "") set uu=$myUU      
   endif
   set seaWall=""
   if ($mySeaWall != "") set seaWall="-wall $mySeaWall"

   if ($MAGIC == NBZZZZ) set uu=u 
   if ($MAGIC == HR) set uu=u 
    
   if ($phase == ma4) set uu=u
   if ($phase == klst4) set uu=u
   if ($phase == r2g4) set uu=u
   # if ($phase == m4) set uu=u
   if ($phase == snp4) set uu=u
   # if($phase == m4H)set uu=u
   if ($target == introns) set uu=u

   set seaLevel=0
   set seaLevel=3
   if ($mySeaLevel != "") set seaLevel=$mySeaLevel
   set removeLimit=50
   if ($myRemoveLimit != "") set removeLimit=$myRemoveLimit
# echo "hello2 $removeLimit"

   if ($target == SpikeIn) then
      set seaWall="-wall 0"
      set seaLevel=0
      set removeLimit=0
   endif

   set mNam=""
   if ($mySeaLevel != "") set mNam=".seaK$seaLevel.rm$removeLimit"
   if ($mySeaWall != "") set mNam=".seaW$mySeaWall.rm$removeLimit"
# echo "hello3 seaWall=$seaWall mNam=$mNam"
   if ($mNam != "") set compare=""


   touch tmp/GENEINDEX/$MAGIC.$target.$uu.mask
   set mask="-mask tmp/GENEINDEX/$MAGIC.$target.$uu.mask"
   if (-e TARGET/Targets/$species.$target.mask.txt) cat  TARGET/Targets/$species.$target.mask.txt >>  tmp/GENEINDEX/$MAGIC.$target.$uu.mask
   if (-e TARGET/Targets/$species.$target.g4mask.txt) cat  TARGET/Targets/$species.$target.g4mask.txt >>  tmp/GENEINDEX/$MAGIC.$target.$uu.mask
   if (-e TARGET/Targets/$species.SpikeIn.mask.txt) cat  TARGET/Targets/$species.SpikeIn.mask.txt >>  tmp/GENEINDEX/$MAGIC.$target.$uu.mask

  if ($phase == r2g4) then
    set toto=tmp/RNA_editing/$MAGIC.run2gene.r2g    
    if (! -e $toto) then
      echo -n "## /RNA_editing/run2gene.r2g : " > $toto
      date >> $toto
      foreach run (`cat MetaDB/$MAGIC/RunList`)
        cat tmp/RNA_editing/$run/f2.*.r1.av.5_A2G.genes  | gawk '{n[$1]+=$2;}END{for (g in n) printf("%s\t%s\t%d\n",g,run,n[g]);}' run=$run >> $toto
      end
      # cat $toto | gawk '/^#/{print;next;}{n[$1]+=$3;}END{for (g in n) printf("%s\tany\t%d\n",g,n[g]);}' | sort -k 3n | tail -60
    endif
    if(-e $toto) then
      cat $toto | gawk '/^#/{print;next;}{printf ("Gene %s\nRun_U %s %f %f seq %f tag %d bp ok\n\n", $1, $2,0,$3,$3,100*$3);}' > tmp/GENEINDEX/$MAGIC.$target.r2g.$uu.ace
      set myace=tmp/GENEINDEX/$MAGIC.$target.r2g.$uu.ace
    endif
  endif

   if (1) then
     cat MetaDB/$MAGIC/RunsList > tmp/GENEINDEX/$MAGIC.$target.$GM.list2
     touch  tmp/GENEINDEX/$MAGIC.$target.$GM.list
     if (-e tmp/GENEINDEX/$MAGIC.$target.$GM.$uu.ace) then
       set n=`diff tmp/GENEINDEX/$MAGIC.$target.$GM.list tmp/GENEINDEX/$MAGIC.$target.$GM.list2 | wc -l`
       if ($n > 0) \rm   tmp/GENEINDEX/$MAGIC.$target.$GM.$uu.ace
     endif
     \mv  tmp/GENEINDEX/$MAGIC.$target.$GM.list2 tmp/GENEINDEX/$MAGIC.$target.$GM.list
   endif

   if ($phase == g4 && ! -e tmp/GENEINDEX/$MAGIC.$target.GENE.$uu.ace) then 
             # if case 2 exist do not retry by using the OTHER_RUNS symbolic link below
     foreach run (`cat MetaDB/$MAGIC/RunsList`)
         if (-e  tmp/GENERUNS/$run/$run.SpikeIn.GENE.$uu.geneSupport.ace.gz) then
           gunzip -c  tmp/GENERUNS/$run/$run.SpikeIn.GENE.$uu.geneSupport.ace.gz >> tmp/GENEINDEX/$MAGIC.$target.GENE.$uu.ace
         endif
         if (-e  tmp/GENERUNS/$run/$run.$target.GENE.$uu.geneSupport.ace.gz) then
           gunzip -c  tmp/GENERUNS/$run/$run.$target.GENE.$uu.geneSupport.ace.gz >> tmp/GENEINDEX/$MAGIC.$target.GENE.$uu.ace
           continue  
         endif
         if (-e  tmp/GENEINDEX/OTHER_RUNS/$run/$run.$target.GENE.$uu.geneSupport.ace.gz) then
           gunzip -c  tmp/GENEINDEX/OTHER_RUNS/$run/$run.$target.GENE.$uu.geneSupport.ace.gz >> tmp/GENEINDEX/$MAGIC.$target.GENE.$uu.ace
         endif
     end
   endif
   
   set deepTsf=""
   if ($phase == g4sp && ! -e tmp/GENEINDEX/$MAGIC.$target.$GM.$uu.ace) then
     set deepTsf="-deepTsf"
     foreach run (`cat MetaDB/$MAGIC/RunsList `)
       foreach WG (WIGGLERUN WIGGLEGROUP)
         if (! -d tmp/$WG/$run) continue
         set tt=tmp/$WG/$run/$target.noCapture.tsf
         if (-e $tt) cat $tt >> tmp/GENEINDEX/$MAGIC.$target.$GM.$uu.ace
       end
     end
   endif
   if ($phase == g4sp && ! -e tmp/GENEINDEX/$MAGIC.$target.$GM.$uu.ace) then 
     if (-e                   tmp/GENEINDEX/$MAGIC.$target.$GM.$uu.ace) \rm tmp/GENEINDEX/$MAGIC.$target.$GM.$uu.ace
             # if case 2 exist do not retry by using the OTHER_RUNS symbolic link below
     foreach run (`cat MetaDB/$MAGIC/RunsList `)
         set long=`cat MetaDB/$MAGIC/RunNanoporeList  MetaDB/$MAGIC/RunPacBioList |gawk '{if($1==run)ok=1;}END{print ok+0;}' run=$run`
         if ($long == 1) then
           # ls -ls tmp/SPONGE/$run/$target.mrna.v4.3.$uu.ns.1
           set myGM=gene
           if ($GM == GENESPX) set myGM=mrna
           set n=`ls  tmp/SPONGE/$run/$target.$myGM.v4.*.$uu.[fr].1 | gawk '{n++;}END{print n+0;}'`
           if ($n  > 0) then
             cat tmp/SPONGE/$run/$target.$myGM.v4.*.$uu.[fr].1 | gawk -F '\t' '/^#/{next;}{g=$3;nn[g]+=$11;}END{for(g in nn){z=nn[g]/100;printf("Gene \"%s\"\nRun_U %s 0.00 %.2f seqs %.2f tags %.2f kb\n\n",g,run,z,z,z/10);}}' run=$run >> tmp/GENEINDEX/$MAGIC.$target.$GM.$uu.ace
           endif
         else
           if (-e  tmp/GENERUNS/$run/$run.$target.GENE.$uu.geneSupport.ace.gz) then
             gunzip -c  tmp/GENERUNS/$run/$run.$target.GENE.$uu.geneSupport.ace.gz >> tmp/GENEINDEX/$MAGIC.$target.$GM.$uu.ace
           endif
         endif
     end
   endif

   if ($phase =~ SF* && ! -e tmp/GENEINDEX/$MAGIC.$target.$phase.$uu.ace) then 
     foreach run (`cat MetaDB/$MAGIC/RunList`)
         if (-e  tmp/Sailfish/$run/$run.$target.$phase.$uu.geneSupport.ace.gz) then
           gunzip -c tmp/Sailfish/$run/$run.$target.$phase.$uu.geneSupport.ace.gz >> tmp/GENEINDEX/$MAGIC.$target.$phase.$uu.ace
         endif
     end
   endif

# ls -ls tmp/GENEINDEX/$MAGIC.$target.$uu.ace
   if ($phase == m4 && ! -e  tmp/GENEINDEX/$MAGIC.$target.MRNA.$uu.ace) then 
     foreach run (`cat MetaDB/$MAGIC/RunsList`)
         if (-e  tmp/GENERUNS/$run/$run.$target.MRNA.$uu.mrnaSupport.ace.gz) then
           gunzip -c  tmp/GENERUNS/$run/$run.$target.MRNA.$uu.mrnaSupport.ace.gz >> tmp/GENEINDEX/$MAGIC.$target.MRNA.$uu.ace
         endif
     end
   endif
   if ($phase == m4H && ! -e  tmp/GENEINDEX/$MAGIC.$target.MRNAH.$uu.ace) then 
     foreach run (`cat MetaDB/$MAGIC/RunsList`)
         if (-e  tmp/GENERUNS/$run/$run.$target.MRNAH.$uu.mrnaSupport.ace.gz) then
           gunzip -c  tmp/GENERUNS/$run/$run.$target.MRNAH.$uu.mrnaSupport.ace.gz >> tmp/GENEINDEX/$MAGIC.$target.MRNAH.$uu.ace
         endif
     end
   endif


# -seaLevel $seaLevel $seaWall -removeLimit $removeLimit
set out=junkyard

set rjm=""
if (-e RESULTS/baddy_batchies.txt) set rjm="-rejectMarker RESULTS/baddy_batchies.txt" 

if ($target == introns) then

    set out=$MAGIC$mNam.$target.INTRON.u.$SNPCHROM$CAPT
    if (-e tmp/INTRON_DB/$chrom/d5.$MAGIC.de_uno$CAPT.ace) then
      echo "bin/geneindex -deepIntron tmp/INTRON_DB/$chrom/d5.$MAGIC.de_uno$CAPT.ace -u $mask $chromAlias -runList MetaDB/$MAGIC/GroupsRunsListSorted -runAce tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace -compare -o  tmp/GENEINDEX/Results/$out -gzo $method $seedGene $sg $capt $rjm $refG -export ait  $shA "

# $seedGene  $compare $correl
      \rm tmp/GENEINDEX/$out.*
      if (1) then
            bin/geneindex -deepIntron tmp/INTRON_DB/$chrom/d5.$MAGIC.de_uno$CAPT.ace -u $mask $chromAlias -runList MetaDB/$MAGIC/GroupsRunsListSorted -runAce tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace  -compare -o  tmp/GENEINDEX/Results/$out -gzo $method $seedGene $sg $capt $rjm $refG -export ait $shA
         if (! -e tmp/GENEINDEX/Results/$out.done) then
           echo "FATAL ERROR inside bin/geneindex : failed to create  tmp/GENEINDEX/Results/$out.done"
           exit 1
         endif
      endif
    endif
endif 

echo "geneindex.tcsh: phase=$phase GM=$GM target=$target"
#cat tmp/GENEINDEX/Results/testA2.introns.INTRON.u.testA2.score.genes.profiles.txt  | grep 9__33988970_33986836


if ($phase == g4 || $phase == g4sp || $phase == gsnp4 ||  $phase == ma4  || $phase == snp4  || $phase == r2g4 || $phase =~ SF*) then

  set myace=JUNK1209983485    
  if ($phase == r2g4 && -e tmp/GENEINDEX/$MAGIC.$target.r2g.$uu.ace) then
    set myace=tmp/GENEINDEX/$MAGIC.$target.r2g.$uu.ace
  endif
  if (-e   tmp/GENEINDEX/$MAGIC.$target.$GM.$uu.ace) set myace=tmp/GENEINDEX/$MAGIC.$target.$GM.$uu.ace
  if (-e   tmp/GENEINDEX/$MAGIC.$target.$phase.$uu.ace) set myace=tmp/GENEINDEX/$MAGIC.$target.$phase.$uu.ace
  if ($phase == ma4 && -e OTHER_PIPELINES/$target.ace) set myace=OTHER_PIPELINES/$target.ace
  set targetSNP=$targetBeau
  set selectSNP=""
  if ($phase == snp4) then
    if (! -e tmp/GENEINDEX/$MAGIC.characteristic_snps.tsf.gz) then
      echo "missing file tmp/GENEINDEX/$MAGIC.characteristic_snps.tsf.gz, please =run please run phase SNV"
      exit 1
    endif
    set myace=tmp/GENEINDEX/$MAGIC.characteristic_snps.tsf.gz
    set targetSNP=$target.$SNPCHROM
    if ($SNPCHROM == "" ) set targetSNP=$target
    set mNam=""
    set tgcl=""
    set  selectSNP="-digital 1 -ratio_bound 1.5"
    set  selectSNP=""
  else 
    if ($phase == gsnp4) then
       set myace=tmp/GENEINDEX/$MAGIC.snp_per_gene.protein.ace
       set  tgcl="-target_class $target_class"    
    else
      if ($phase != ma4) then
        set  tgcl="-target_class $target_class"    
      else
        set tgcl=""
      endif
    endif
  endif

  if (! -e $myace) then 
    echo "FATAL ERROR: Phase $phase cannot access the expression data $myace"
    exit 1
  endif
# ls -ls $myace
   \rm  tmp/GENEINDEX/Results/$MAGIC$mNam.$targetSNP.$GM.$uu.*

     set out=$MAGIC$mNam.$targetSNP.$GM$CAPT.$uu
     set dg=deepGene
     if ($phase == snp4) set dg=deepSNP

echo "geneindex.tcsh8: phase=$phase GM=$GM target=$target out=$out MAGIC=$MAGIC"

     set testtest=""
     \rm tmp/GENEINDEX/Results/$out.*
     
 # gene_group
set GeneGroup=""
if ($MAGIC_GENE_GROUP == 1 &&   -e TARGET/GENES/Gene_groups.ace) then
  set GeneGroup="-geneGroup TARGET/GENES/Gene_groups.ace"
endif

     echo "aaa9 phase $phase  myace = $myace"
     set splitM=""
     if ($phase == snp4) then
       set params="-$dg $myace -$uu $mask $chromAlias -runList MetaDB/$MAGIC/GroupsRunsListSorted -runAce tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace  -o  tmp/GENEINDEX/Results/$out -gzo $method $selectSNP  $shA $sg $capt $refG  $tgcl  $compare  $rjm -htmlSpecies $species "
     else
       if (0 && $phase == g4 && -e TARGET/GENES/$species.$target.split_mrnas.txt) then  
         # we cannot do this here, it must be done in phase g2a: collect hits per gene in pgm bestali
         set splitM="-split_mRNAs TARGET/GENES/$species.$target.split_mrnas.txt"
       endif
       set params="-$dg $myace -$uu $mask $chromAlias -runList MetaDB/$MAGIC/GroupsRunsListSorted -runAce tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace $deepTsf -o  tmp/GENEINDEX/Results/$out -gzo  -pA $method $selectSNP  $shA $sg $capt $refG  $tgcl $GeneGroup $correl $compare  $rjm  -htmlSpecies $species $splitM  -export aitvz "
     endif

     if (1) then
         echo "AAAA bin/geneindex $params"
               bin/geneindex $params
        if (! -e tmp/GENEINDEX/Results/$out.done) then
          echo "FATAL ERROR inside bin/geneindex : failed to create  tmp/GENEINDEX/Results/$out.done"
          exit 1
       endif
    endif

endif

#### 
# hack for sept 9

if (0) then 
  cat MetaDB/$MAGIC/gtitle.txt ZZZZZ RESULTS/SNPclassif.txt | gawk -F '\t' '/^ZZZZZ/{zz++;nextt;}{if(zz<1){gsub(/\"/,"",$0);r=$1;st[r]=$11;ot[r]=$12;sample[r]=$3;next;}}/^# Run/{print;pass++;if(pass!=2)next;printf("# Run\tSorting_title");for(i=3;i<=NF;i++)printf("\t%s",st[$i]);printf("\n");printf("# Run\tOther_title");for(i=3;i<=NF;i++)printf("\t%s",ot[$i]);printf("\n");printf("# Run\tSample");for(i=3;i<=NF;i++)printf("\t%s",sample[$i]);printf("\n");next;}' > RESULTS/SNPclassif.with_titles.txt

endif

####

if ($phase == m4 || $phase == m4H || $phase == klst4 || $phase == ii4) then

  set  tgcl="-target_class $target_class"    
  set GM=MRNA
  if ($phase == m4H) set GM=MRNAH
  set isMRNAH=''
  if ($phase == m4H)  set isMRNAH='-MRNAH'
  if ($phase == klst4 ) set GM=Kallisto

  if ($phase == klst4) then
    set myace=tmp/GENEINDEX/$MAGIC.$target.Kallisto.$uu.ace
    if (! -e $myace) then 
      echo "... constructing $myace"
      echo ' ' > $myace
      foreach run (`cat MetaDB/$MAGIC/RunList`)
        if (-e  tmp/Kallisto/$target/$run/abundance.tsv) then
          cat tmp/Kallisto/$target/$run/abundance.tsv | gawk -F '\t' '{x=$5+0;z=0;if(x>0)z=log(1000*x)/log(2);if(z<4)z=4;printf("Transcript %s\nRun_U %s %f %f seq %f tag %d bp ok\n\n", $1, run, z,$4,$4,100*$4);}' run=$run >> $myace
          # Transcript UKv4_A_23_P314216
          # Run_U Rhs1044 8.471  65 seq 65 tag 6500 bp  ok
        endif
      end
    else
      echo "... found $myace"
    endif
  endif


# ls -ls tmp/GENEINDEX/$MAGIC.$target.$GM.$uu.ace
  set myace=JUNK1209983485
  if ($phase == ii4 && -e tmp/INTRON_DB/$chrom/ii4.index.ace) set myace=tmp/INTRON_DB/$chrom/ii4.index.ace

echo "... $phase $target myace=$myace"
# gene_group

set GeneGroup=""
if ($MAGIC_GENE_GROUP == 1 && -e TARGET/GENES/Gene_groups.ace) then
  set GeneGroup="-geneGroup TARGET/GENES/Gene_groups.ace"
endif

  if (-e $myace) then 
echo "... $phase $target $myace MAGIC=$MAGIC GeneGroup=$GeneGroup"
     set out=$MAGIC$mNam.$targetBeau.$GM$CAPT.$uu
     \rm tmp/GENEINDEX/Results/$out.*
     \rm RESULTS/Expression/AceFiles/$out.*
     touch tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace 
     echo "bin/geneindex -deepTranscript $myace -$uu $mask $chromAlias -runList MetaDB/$MAGIC/GroupsRunsListSorted -runAce tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace  -o tmp/GENEINDEX/Results/$out -gzo -pA $method $isMRNAH $tgcl  $shA $sg $capt  $refG $rjm $GeneGroup -htmlSpecies $species   -export aitvz  $correl $compare" 
           bin/geneindex -deepTranscript $myace -$uu $mask $chromAlias -runList MetaDB/$MAGIC/GroupsRunsListSorted -runAce tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace  -o tmp/GENEINDEX/Results/$out -gzo -pA $method $isMRNAH $tgcl  $shA $sg $capt  $refG $rjm  $GeneGroup -htmlSpecies $species   -export aitvz  $correl $compare
     if (! -e tmp/GENEINDEX/Results/$out.done) then
       echo "FATAL ERROR inside bin/geneindex : failed to create  tmp/GENEINDEX/Results/$out.done"
       exit 1
     endif
  endif

endif

set nruns=`wc MetaDB/$MAGIC/RunsList | gawk '{print $1;}'`
echo -n "uuu0  geneindex.c $target phase=$phase done succesfully nruns=$nruns : "
# cat tmp/GENEINDEX/Results/testA2.introns.INTRON.u.testA2.score.genes.profiles.txt  | grep 9__33988970_33986836

date
   if ($nruns < 24 && $phase == g4 && -d GeneIndexDB && ! -e GeneIndexDB/database/lock.wrm &&  -e tmp/GENEINDEX/Results/$MAGIC$mNam.$target.$uu.ace.gz) then
     bin/tacembly GeneIndexDB <<EOF
       read-models
       pparse tmp/GENEINDEX/Results/$MAGIC$mNam.$target.$uu.ace.gz
       save
       quit
EOF
   endif
 
  foreach type (expression_index reads_aligned_per_gene sFPKM)
    if (-e tmp/GENEINDEX/Results/$out.$type.txt.gz) then
      set ff=tmp/GENEINDEX/Results/$out.$type
      gunzip -c  $ff.txt.gz  | head -300 > $ff.topOfFile.txt
      gunzip -c  $ff.txt.gz  | sort  | gawk 'BEGIN{print "\n\n\n\n\n\n\n\n"}/^#/{next;}{n++;if(n<=300)print;}'  >> $ff.topOfFile.txt
    endif
  end

## sex ratio
  if ($phase == g4 && $uu == u) then
    set ok=0 
    if (-e TARGET/GENES/Gene_groups.ace) then
      set ok=`cat TARGET/GENES/Gene_groups.ace | grep Magic_Female | wc -l`
    endif
    if ($ok == 1) then
      set ff=tmp/GENEINDEX/Results/$out.ace.gz
      if (-e $ff) then
        gunzip -c $ff | gawk '/^Gene/{ok=0;if(okok==3)last;}/^Gene \"Magic_Male/{ok=1;okok+=1;next;}/^Gene \"Magic_Female/{ok=2;okok+=2;next;}{if(ok == 0)next;}/^(Run|Group)_U/{run=$2;runs[run]+=ok;z[run,ok]=$3;next;}/^Group_U/{run=$2;runs[run]+=ok;z[run,ok]=$3;next;}END{for(run in runs){ if(runs[run]==3){z1 = z[run,1]+0;;z2 = z[run,2]+0;z3=(z1 - z2 + 7);printf("Ali %s\nSex_ratio %s %.2f %.1f Male %.1f Female\n\n",run, target,z3, z1,z2);}}}' target=$target > tmp/GENEINDEX/Results/$out.sex_ratio.ace
        echo "pparse  tmp/GENEINDEX/Results/$out.sex_ratio.ace" | bin/tacembly MetaDB -noprompt 
      endif
    endif
  endif


  set qu="unique"
  if ($uu == nu) set qu="quasi_unique"

  if (! -d RESULTS/Expression) mkdir RESULTS/Expression 
  if (! -d RESULTS/Expression/$qu) mkdir RESULTS/Expression/$qu
  if (! -d RESULTS/Expression/$qu/$target) mkdir RESULTS/Expression/$qu/$target
  if (! -d RESULTS/Expression/AceFiles) mkdir RESULTS/Expression/AceFiles
  if (! -d RESULTS/Expression/$qu/$target/Diff_genes) mkdir RESULTS/Expression/$qu/$target/Diff_genes
  if (-d RESULTS/Expression/$qu/$target/Differential_runs  && ! -d RESULTS/Expression/$qu/$target/Pairs) mv RESULTS/Expression/$qu/$target/Differential_runs  RESULTS/Expression/$qu/$target/Pairs
  if (! -d RESULTS/Expression/$qu/$target/Pairs) mkdir RESULTS/Expression/$qu/$target/Pairs

  \rm RESULTS/Expression/*/$out.* RESULTS/Expression/$qu/$out.*
     echo -n "mv tmp/GENEINDEX/Results/$out.*  RESULTS/Expression/$qu/$target : "
# cat tmp/GENEINDEX/Results/testA2.introns.INTRON.u.testA2.score.genes.profiles.txt  | grep 9__33988970_33986836

     date
  gunzip tmp/GENEINDEX/Results/$out.sFPKM.txt.gz
  mv tmp/GENEINDEX/Results/$out.*.ace tmp/GENEINDEX/Results/$out.*.ace.gz  tmp/GENEINDEX/Results/$out.ace.gz  RESULTS/Expression/AceFiles
  gzip -f RESULTS/Expression/AceFiles/*.ace
# cat RESULTS/Expression/unique/introns/testA2.introns.INTRON.u.testA2.score.genes.profiles.txt  | grep 9__33988970_33986836

# Reorganize the compare files 
  if (1) then
    foreach ff (`ls tmp/GENEINDEX/Results/$out.*.compare.txt`)
      mv $ff $ff.1
      cat $ff.1 | gawk -F '\t' '/^#/{if (NF < 20){print;next;}}/^$/{zz=1;line=1;next;}{zz=zz+0;line++;if(nf<NF)nf=NF;for(i=1;i<=NF;i++)aa[i,line,zz]=$i;if(line > lineMax)lineMax = line;}END{for (line=1;line<=lineMax;line++){printf("%s",aa[1,line,0]); for(i=2;i<=nf;i++)printf("\t%s",aa[i,line,0]);for(i=nf+1;i<100;i++)printf("\t");u=1;if(line==1)u=0;for(i=1;i<=nf;i++)printf("\t%s",aa[i,line,u]);printf("\n");}}' > $ff
     \rm $ff.1
    end
  endif

  \rm  tmp/GENEINDEX/Results/$out.*done
  mv tmp/GENEINDEX/Results/$out.*  RESULTS/Expression/$qu/$target
  mv RESULTS/Expression/$qu/$target/$out.*.DEG_FDR_selection.txt tmp/GENEINDEX/Results
  mv RESULTS/Expression/$qu/$target/$out.*.histo.*.txt* RESULTS/Expression/$qu/$target/Pairs
  mv RESULTS/Expression/$qu/$target/$out.*.alpha.*.txt* RESULTS/Expression/$qu/$target/Pairs
  mv RESULTS/Expression/$qu/$target/$out.*.beta.*.txt* RESULTS/Expression/$qu/$target/Pairs
  mv RESULTS/Expression/$qu/$target/$out.*.diffGenes.*.txt* RESULTS/Expression/$qu/$target/Diff_genes 

  # gene type profile per prevalence in the population

  foreach ff (`ls RESULTS/Expression/$qu/$target/Diff_genes/$MAGIC.$targetBeau.$GM.$uu.*diffGenes.[0].txt`)
    cat $ff | gawk -F '\t' '{n++;k=NF;score=$5;if(n>200)k=14;if (n<20 || score > 180){printf("%s",$1);for(i=2;i<=k;i++)printf("\t%s",$i);printf("\n");}}' > $ff.light.txt
  end

    foreach ff (`ls RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.*.$uu.*.*diffGenes.0.txt`)
      touch  tmp/METADATA/$target.metadata.txt
 
     cat  tmp/METADATA/$target.metadata.txt ZZZZZ $ff | gawk -F '\t' '/log2/{for(i=1;i<=NF;i++)if(index($i,"log2")>0)nFC=i;}/^#/{ndiese++;if (ndiese==1)printf("%s\t%s.with_metadata.txt\n",$1,$2);else print; next;}/^ZZZZZ/{zz++;next;}{if(zz<1){g = $1 ; gMap[g]=$2 ; gid[g]=$3 ;gType[g]=$4 ;gTitle[g]=$5;next;}}/^Genes/{print;next;}/\tScore/{printf("Species\tComparison\tOrdering number\tGene or transcript\tDifferential score, range [0,200]\tlog2 of expression fold change\tAverage expression index of low group (log2(1000 sFPKM))\tAverage expression index of high group (log2(1000 sFPKM))\tNCBI GeneId\tMap on chromosome:from-to on build 37\tType\tGene description\tLink (in Excel please click twice and change line to activate the link)\t\n");next;}{if($4==old4 && ((nn<20 && $5 >= 100) || $5 >= minScore)){g=$4;nn++;printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\twww.ncbi.nlm.nih.gov/AceView/av.cgi?db=%s&q=%s\n",old1,old2,old3,$4,$5,- $nFC,$(nFC + 1),z,gid[g],gMap[g],gType[g],gTitle[g],db,g,g);}z=$(nFC+1);old1=$1;old2=$2;old3=$3;old4=$4;next;}' db=$species minScore=1  > $ff.with_metadata.txt
      set n=`wc -l $ff.with_metadata.txt | gawk '{print $1}'`
      if ($n > 2000) then
        cat $ff.with_metadata.txt | head -1000 > $ff.with_metadata.1000.txt
       endif
    end
  # end

echo "phase $phase fuse the FDR selections files into a single file per project"
# fuse the FDR selections files into a single file per project and collect the thresholds at the top of the file
set toto1=tmp/GENEINDEX/Results/$MAGIC.$targetBeau.$GM.$uu
set toto=RESULTS/Expression/$qu/$target/$MAGIC.$targetBeau.$GM.$uu.DEG_FDR_stats.txt

echo -n "##\t" > $toto
cat  $toto1.*.DEG_FDR_selection.txt | head -1 | cut -f 2 | gawk '{printf("%s\t",$1);}' >> $toto
echo " file=$toto" >> $toto
echo >> $toto
echo "# Compare\tGenes overexpressed in\trelative to\tNumber of runs in first group\tNumber of runs in second group\tMagic DEG score threshold (0 to 200)\tNumber of DEGs after noise substraction\t% FDR" >> $toto

echo "cat   $toto1.*.DEG_FDR_selection.txt "
cat   $toto1.*.DEG_FDR_selection.txt | gawk -F '\t' '{if ($6 == "Selected threshold")print}' | sort | cut -f 1,2,3,4,5,7,9,11 >> $toto

set toto=RESULTS/Expression/$qu/$target/$MAGIC.$targetBeau.$GM.$uu.DEG_FDR_selection.txt
echo -n "## $toto : " > $toto
date >> $toto
echo >> $toto
cat  $toto1.*.DEG_FDR_selection.txt >> $toto


# hack for the primates project
if (0) then

 cat  tmp/METADATA/$target.metadata.txt ZZZZZ  RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.av.*.$uu.*[Nn]o*.*diffGenes.0.txt | gawk -F '\t' '/log2/{if(nlog>0)next;nlog=1;for(i=1;i<=NF;i++)if(index($i,"log2")>0)nFC=i;}/^#.*file/{if(uuFile<1)print;uuFile=1; next;}/^ZZZZZ/{zz++;next;}/Testis/{next;}/Ovary/{next;}{if(zz<1){g = $1 ; gMap[g]=$2 ; gid[g]=$3 ;gType[g]=$4 ;gTitle[g]=$5;next;}}/^Genes/{if(uuG<1)print;uuG=1;next;}/\tScore/{printf("Species\tComparison\tOrdering number\tGene or transcript\tDifferential score, range [0,200]\tlog2 of expression fold change\tAverage expression index of low group (log2(1000 sFPKM))\tAverage expression index of high group (log2(1000 sFPKM))\tNCBI GeneId\tMap on chromosome:from-to on build 37\tType\tGene description\tLink (in Excel please click twice and change line to activate the link)\t\t\t\t\t\t\t\t\t\t\t\t\t");for(i=0;i<=260;i+=5)printf("\t%.1f",i/10);printf("\n");next;}{if(ok[$4]<1 && $4==old4 && (oldzz[12] == 200 || $5 >= minScore)){g=$4;ok[g]=1;nn++;for(ii=0;ii<2;ii++){printf("%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\twww.ncbi.nlm.nih.gov/AceView/av.cgi?db=%s&q=%s",old1,old2,old3,$4,$5,- $nFC,$(nFC + 1),z,gid[g],gMap[g],gType[g],gTitle[g],db,g,g);for(i=6;i<=NF;i++){if(ii==0)printf("\t%s",oldzz[i]);else printf("\t%s",$i);}printf("\n");}}z=$(nFC+1);old1=$1;old2=$2;old3=$3;old4=$4;for(i=5;i<=NF;i++)oldzz[i]=$i;next;}' db=$species minScore=200 > RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.av.GM.u.score.AUC200.with_metadata_and_histo.txt

cat  RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.av.GM.u.score.AUC200.with_metadata_and_histo.txt | gawk '/Testis/{next;}{print}END{print "ADCY8";}' | cut -f 4 | gzip > toto.glist.gz
gunzip -c toto.glist.gz Dan.score200.list.gz | sort -u | gzip > toto.glist2.gz
gunzip -c  toto.glist2.gz ZZZZZ.gz RESULTS/Expression/unique/av/Tissue.av.*.u.expression_index.txt.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}/^#/{print;next;}{if(ok[$1]==1)print;}' >  RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.av.GM.u.score.AUC200.expression_index.txt
wc  RESULTS/Expression/unique/$target/Diff_genes/$MAGIC.av.GM.u.score.AUC200.expression_index.txt

gunzip -c RESULTS/Expression/unique/av/Tissue.av.MRNAH.u.expression_index.txt.gz | gawk -F '\t' '/^#/{print}/^ZFAND6/{print}' >  RESULTS/Expression/unique/$target/ZAAND6.expressio.MRNAH.txt
mv RESULTS/Expression/unique/$target/ZFAND6.expressio.MRNAH.txt RESULTS/Expression/unique/$target/ZFAND6.expression.MRNAH.txt

endif


  if (-e RESULTS/Expression/AceFiles/$out.withIndex.ace.gz && $uu == u && ($phase == g4 || $phase == m4H || $phase == g4sp)) then 
       echo "pparse RESULTS/Expression/AceFiles/$out.withIndex.ace.gz" | bin/tacembly MetaDB -no_prompt
  endif

  if (0) then
    set tutu=RESULTS/Expression/$qu/$target/$MAGIC.$target.$GM.$uu.selected_index_variance.txt
    echo " which are the least variable genes $tutu"
    echo -n "# $tutu : " > $tutu
    date >> $tutu
  
    gunzip -cf MetaDB/$MAGIC/GroupVarianceList ZZZZZ RESULTS/Expression/AceFiles/$out.ace.gz  | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){nok++;g2n[$1]=nok;n2g[nok]=$1;next;}}/^Gene/{g=$2;gsub(/\"/,"",g);z0=-1;z1=-1;z2=-1;s0=-1;s1=-1;s2=-1;next;}/^Trancript/{g=$2;gsub(/\"/,"",g);z0=-1;z1=-1;z2=-1;s0=-1;s1=-1;s2=-1;next;}/^Transcript/{g=$2;gsub(/\"/,"",g);z0=-1;z1=-1;z2=-1;s0=-1;s1=-1;s2=-1;next;}/^Group_/{gr=$2;n=g2n[gr];if(n>0){z[n]=$3;if($29!="ok")ne[n]=$29 ":";var[n]=$30;}next;}/^$/{if (g != ""){ngg++;if(ngg==1){printf("Gene\tType");for(n=1;n<=nok;n++)printf("\t%s index\t%s variance",n2g[n],n2g[n]);printf("\n");}printf("%s\t%s",g,type);for(n=1;n<=nok;n++){printf("\t%s%.2f\t%s%.2f",ne[n],0+z[n],ne[n],0+var[n]);z[n]=0;ne[n]="";var[n]="";var[n]=-999;}}g="";printf("\n") ;}' type="$target.$GM.$uu" >> $tutu &

  endif


  if ($phase == gsnp4) then
    if (! -d RESULTS/SNV_per_gene) mkdir  RESULTS/SNV_per_gene
    if (! -d RESULTS/SNV_per_gene/$target) mkdir  RESULTS/SNV_per_gene/$target
    mv RESULTS/Expression/$qu/$target/$MAGIC.$targetBeau.$GM.* RESULTS/Expression/$qu/$target/Di*/$MAGIC.$targetBeau.$GM.* RESULTS/SNV_per_gene/$target
  endif
  if ($phase == snp4) then
    if (! -d RESULTS/SNV) mkdir  RESULTS/SNV
    \mv RESULTS/Expression/$qu/$target/$MAGIC.* RESULTS/SNV
  endif

end 



goto phaseLoop

foreach ff (`ls RESULTS/Expression/unique/*/NB.*.u.expression_index.txt.gz`)
  gunzip -c $ff | head -100 > $ff.test.txt
  gunzip -c $ff | head -1000 | tail -100 >> $ff.test.txt
  gunzip -c $ff | tail -100 >> $ff.test.txt
end

# cat tutu.txt | gawk '{gsub(/\"/,"",$0);}/_A_/{ra[$1]=$2;next;}/_C_/{gsub(/_C_/,"_A_",$1);rc[$1]=$2;next;}/_C_/{gsub(/_C_/,"_A_",$1);rc[$1]=$2;next;}/_C_/{gsub(/_D_/,"_A_",$1);rd[$1]=$2;next;}/_B_/{gsub(/_B_/,"_A_",$1);rb[$1]=$2;next;}END{for (a in ra)printf("Run %s\nTitration %s %s %s %s\n\n",ra[a],ra[a],rc[a],rd[a],rb[a]);}' > tutu.ace


# special script to launch in parallel several version for optimizing the sea level
if (0) then
# sea: 0 1 2 3 7 10   rm: 0 50 100 200 500 1000
foreach sea (3)
  foreach rm (50)
    foreach seaWall (2 )
      foreach seaUU (nu)
        if (! -e tmp/GENEINDEX/Results/SEQC_main.seaW$seaWall.rm$rm.av.$seaUU.TitratingGenes.txt2) then
           \rm  toto.$MAGIC.seaW$seaWall.rm$rm.$seaUU
           scripts/submit toto.$MAGIC.seaW$seaWall.rm$rm.av.$seaUU.R.$seaR "scripts/geneindex.tcsh g4 $sea $rm $seaWall $seaUU" 32G
        endif
      end
    end
  end
end
endif

# special script to analyse the output of the above
set uu=u
set toto=RESULTS/$MAGIC.seaLevelSummary.$uu.txt

\rm $toto
date > $toto
echo "Number of genes with index" >> $toto 
cat tmp/GENEINDEX/Results/$MAGIC.sea3.rm50.av.$uu.NumberOfGenesPerRunAboveGivenIndex.txt | head -6 | gawk '{printf("\t");print} ' >> $toto
foreach sea ( 0 1 2 3 5 7 10 )
  foreach lm (0 50 100 200 500 1000)
    if (-e tmp/GENEINDEX/Results/titration.sea$sea.rm$lm.av.$uu.NumberOfGenesPerRunAboveGivenIndex.txt) then
      echo -n "sea$sea.rm$lm\t" >> $toto
      cat tmp/GENEINDEX/Results/titration.sea$sea.rm$lm.av.$uu.NumberOfGenesPerRunAboveGivenIndex.txt | head -7 | tail -1 >> $toto
    endif
  end
  echo >> $toto
end
echo "\n\nTitrating genes at threshold .1 (1.071  fold)" >> $toto
cat tmp/GENEINDEX/Results/$MAGIC.sea3.rm50.av.$uu.NumberOfTitratingGeneAboveGivenIndex.txt | head -20 | gawk '{printf("\t");print} ' >> $toto
foreach sea ( 0 1 2 3 5 7 10 )
  foreach lm (0 50 100 200 500 1000)
     echo -n "sea$sea.rm$lm\t" >> $toto
     cat tmp/GENEINDEX/Results/$MAGIC.sea$sea.rm$lm.av.$uu.NumberOfTitratingGeneAboveGivenIndex.txt | head -22 | tail -1 >> $toto
  end
  echo >> $toto
end
echo "\n\nTitrating genes at threshold 1.0 (2 fold)" >> $toto
cat tmp/GENEINDEX/Results/$MAGIC.sea3.rm50.av.$uu.NumberOfTitratingGeneAboveGivenIndex.txt | head -20 | gawk '{printf("\t");print} ' >> $toto
foreach sea ( 0 1 2 3 5 7 10 )
  foreach lm (0 50 100 200 500 1000)
     echo -n "sea$sea.rm$lm\t" >> $toto
     cat tmp/GENEINDEX/Results/$MAGIC.sea$sea.rm$lm.av.$uu.NumberOfTitratingGeneAboveGivenIndex.txt | head -32 | tail -1 >> $toto
  end
  echo >> $toto
end

echo "\n\nTitrating genes at threshold 2.0 (4 fold)" >> $toto
cat tmp/GENEINDEX/Results/$MAGIC.sea3.rm50.av.$uu.NumberOfTitratingGeneAboveGivenIndex.txt | head -20 | gawk '{printf("\t");print} ' >> $toto
foreach sea ( 0 1 2 3 5 7 10 )
  foreach lm (0 50 100 200 500 1000)
     echo -n "sea$sea.rm$lm\t" >> $toto
     cat tmp/GENEINDEX/Results/$MAGIC.sea$sea.rm$lm.av.$uu.NumberOfTitratingGeneAboveGivenIndex.txt | head -42 | tail -1 >> $toto
  end
  echo >> $toto
end

## same script but more restricted top some conditions
\rm $toto
date > $toto
echo "\n\nTitrating genes at threshold .1 (1.071  fold)" >> $toto
cat tmp/GENEINDEX/Results/$MAGIC.sea3.rm50.av.$uu.NumberOfTitratingGeneAboveGivenIndex.txt | head -20 | gawk '{printf("\t");print} ' >> $toto
foreach sea (1 2 3 5 7)
  foreach lm (50 100 200)
     echo -n "sea$sea.rm$lm\t" >> $toto
     cat tmp/GENEINDEX/Results/$MAGIC.sea$sea.rm$lm.av.$uu.NumberOfTitratingGeneAboveGivenIndex.txt | head -22 | tail -1 >> $toto
  end
  echo >> $toto
end


# generate the macro groups
# A1, B1 unmistakable: seens 14 ILM and 6 LIF among the groups 
# A2, B2 harder 4 ILM 2 LIF
# A3, B3 dubious > 4 total


if ($MAGIC == SEQC_main) then
  # get the list of candidate titrating genes in column 6 to 29 even->A>B   odd->B>A

foreach sea (3)
  foreach rm (1.5 2 2.5 3.0 4)
    foreach seaWall (4)
      foreach seaUU (nu)
          if (-e tmp/GENEINDEX/Results/$MAGIC.seaW$seaWall.rm$rm.av.$seaUU.TitratingGenes.txt) then
            if (-e tmp/GENEINDEX/Results/$MAGIC.seaW$seaWall.rm$rm.av.$seaUU.TitratingGenes.txt.gz) then
              \rm tmp/GENEINDEX/Results/$MAGIC.seaW$seaWall.rm$rm.av.$seaUU.TitratingGenes.txt.gz
            endif
            gzip tmp/GENEINDEX/Results/$MAGIC.seaW$seaWall.rm$rm.av.$seaUU.TitratingGenes.txt
          endif
          if (-e tmp/GENEINDEX/Results/$MAGIC.seaW$seaWall.rm$rm.av.$seaUU.TitratingGenes.txt.gz) then
            \rm  toto.compyatibility.seaW$seaWall.rm$rm.out  toto.compatibility.seaW$seaWall.rm$rm.err
            scripts/submit toto.compatibility.seaW$seaWall.rm$rm "gunzip -c tmp/GENEINDEX/Results/$MAGIC.seaW$seaWall.R$rm.av.$seaUU.ace.gz ZZZZZ.gz tmp/GENEINDEX/Results/$MAGIC.seaW$seaWall.rm$rm.av.$seaUU.TitratingGenes.txt.gz | gawk -F '\t' -f scripts/g4.titration_compatibility.awk title=sea$seaWall.R$rm.$seaUU OGG=tmp/GENEINDEX/Results/$MAGIC.seaW$seaWall.R$rm.$seaUU.titrating_gene_list.txt >  tmp/GENEINDEX/Results/$MAGIC.seaW$seaWall.R$rm.$seaUU.titrating_profile.txt"
          else
            \rm  tmp/GENEINDEX/Results/$MAGIC.sea$sea.rm$rm.titrating_profile.txt
          endif
      end
    end
  end
end 

gzip tmp/GENEINDEX/Results/SEQC_main.av.GENE.$uu.TitratingGenes.txt
gunzip -c tmp/GENEINDEX/Results/SEQC_main.av.GENE.u.TitratingGenes.txt.gz | scripts/transpose | head -119 | cut -f 1,2,3,4,5
gunzip -c tmp/GENEINDEX/Results/SEQC_main.av.GENE.u.TitratingGenes.txt.gz | scripts/transpose | gawk -F '\t' '{if(substr($2,1,3)=="Rhs")next;print}' | scripts/transpose | gzip >  tmp/GENEINDEX/Results/SEQC_main.av.GENE.u.TitratingGenes.groups.txt.gz

gunzip -c tmp/GENEINDEX/Results/$MAGIC.av.GENE.$uu.ace.gz ZZZZZ.gz tmp/GENEINDEX/Results/$MAGIC.av.GENE.$uu.TitratingGenes.groups.txt.gz | gawk -F '\t' -f scripts/g4.titration_compatibility.awk title=av.nu OGG=tmp/GENEINDEX/Results/$MAGIC.$uu.titrating_gene_list.txt >  tmp/GENEINDEX/Results/$MAGIC.GENE.$uu.titrating_profile.txt

echo -n "# " > RESULTS/SEQC_main.GENE.$uu.titrating_gene_list.txt
date >>  RESULTS/SEQC_main.GENE.$uu.titrating_gene_list.txt
head -1 tmp/GENEINDEX/Results/$MAGIC.$uu.titrating_gene_list.txt | gawk -F '\t' '{printf("Gene\tNCBI GeneId\tChromosome\tStrand\tFrom\tTo\tTitle");for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' >> RESULTS/SEQC_main.GENE.$uu.titrating_gene_list.txt
cat  TARGET/GENES/av.gene.ace ZZZZZ  tmp/GENEINDEX/Results/$MAGIC.$uu.titrating_gene_list.txt | gawk -F '\t' '/^IntMap/{gsub(/\"/,"",$2);split($2,aa," "); a1=aa[2];a2=aa[3];s="+";if(a1>a2){s="-";a0=a1;a1=a2;a2=a0;} mm[g]= aa[1] "\t" s "\t" a1 "\t" a2 ; next;}/^Title/{gsub(/\"/,"",$2);gsub(/^Title /,"",$2);title[g]=$2;next;}/^GeneId/{gsub(/\"/,"",$2);gid[g]=$2;next;}/^Gene /{gsub(/\"/,"",$1);gsub(/^Gene /,"",$1);g=$1;}/^ZZZZZ/{zz=1;next;}{if(zz<1)next;}{g=$1;printf("%s\t%s\t%s\t%s",g,gid[g],mm[g],title[g]);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}'  >> RESULTS/SEQC_main.GENE.$uu.titrating_gene_list.txt

# read the titrating flags in GeneIndexDB
cat RESULTS/SEQC_main.GENE.$uu.titrating_gene_list.txt | gawk -F '\t' '/^#/{next;}{if($1 == "Gene")next;printf("Gene %s\n",$1); if($11>0)printf("SEQC_acdb_minus\n");else printf("SEQC_acdb_plus\n");printf("SEQC_A %s\nSEQC_B %s\n\n",$8,$9);}' > tmp/GENEINDEX/Results/$MAGIC.titrating_gene_list.GENE.ace
tbly GeneIndexDB <<EOF
  read-models
  parse  tmp/GENEINDEX/Results/$MAGIC.titrating_gene_list.GENE.ace
  save
  quit
EOF

# construct the POG list among the A class
cat RESULTS/SEQC_main.$uu.titrating_gene_list.txt | cut -f 1,11,20,21,22,23,24,25 | gawk '/^#/{next}{if($2>1 || $2<-1)next;zi=0;zl=0;if($4=="AI" || $4=="BI")zi+=1;if($5=="AL" || $5=="BL")zl+=1;if($7=="AI" || $7=="BI")zi+=2;if($8=="AL" || $8=="BL")zl+=2;if(zi > 0)N++;if(zl>0)N++;if(zi > 0 && zi == zl)ok+=2 ;}END{printf("N=%d ok=%d POG=%.2f\n",N, ok, 100.0*ok/N);}'


# construct the histograms
cat RESULTS/SEQC_main.$uu.titrating_gene_list.txt | gawk -F '\t' '/^ERCC/{next}/^#/{next}{ii=$20;jj=$23;if(ii)t=ii;else t=jj;if (t+0>0){printf("t=%s\t",t);print;}if(t){d=$11;if(d<-20)d=-20;if(d>20)d=20;u=1;if(d<0){d=-d;u=-1;}d=u*int(10*d);tt[t]=1;fc[t,d]++;ffc[t]++;}}END{printf("Delta\tA\tB\tC\tD\tE\tF\tG\nTotal");for(i=1;i<=7;i++)printf("\t%d",ffc[substr("ABCDEFG",i,1)]);for(d=-200;d<=200;d++){printf("\n%.1f",d/10.0);for(i=1;i<=7;i++)printf("\t%d",fc[substr("ABCDEFG",i,1),d]);}printf("\n");}' >  RESULTS/SEQC_main.$uu.titrating_gene.volcano.txt


# cat titi1 | gawk -F '\t' '{if($6 == "A")print int(100 * $5)}' | bin/histo -plain -o toto88
# cat titi1 | gawk -F '\t' '{if($11>0 && $12>0 && $4 > 10 )print $5}' | bin/histo -plot -plain -o toto.c


foreach col (pc pp gc gp tc tp bad)
  if ($col == pc) set tt=perfect_count
  if ($col == pp) set tt=perfect_ratio
  if ($col == gc) set tt=good_count
  if ($col == gp) set tt=good_ratio
  if ($col == tc) set tt=titrating5_count
  if ($col == tp) set tt=titrating5_ratio

  if ($col == bad) set tt=baddies
  set toto=RESULTS/$MAGIC.titration_compatibility.$uu.$tt.txt
  date > $toto

   if ($col == pc) echo "Counts of genes called titrating in each experiment belonging to the consistant-and-easy-to-measure titrating lists (category A1 or B1, seen in at least 14 ILM and 6 LIF in the 49 macro groups" >> $toto
   if ($col == pp) echo "Percentage of genes called titrating in each experiment belonging to the consistant-and-easy-to-measure titrating lists (category A1 or B1, seen in at least 14 ILM and 6 LIF in the 49 macro groups" >> $toto


   if ($col == gc) echo "Counts of genes called titrating in each experiment confirmed (category A1/A2 or B1/B2, seen in at least 4 ILM and 2 LIF in the 49 macro groups and not contradicted in any macro group" >> $toto
   if ($col == gp) echo "Percentage of genes called titrating in each experiment confirmed (category A1/A2 or B1/B2, seen in at least 4 ILM and 2 LIF in the 49 macro groups and not contradicted in any macro group" >> $toto

   if ($col == tc) echo "Count genes called titrating in each experiment belonging to the confirmed titrating lists (category A1/A2/A3 or B1/B2/B3, seen in at least 5 ILM or LIF in the 54 macro groups" >> $toto
   if ($col == tp) echo "Percentage of genes called titrating in each experiment belonging to the confirmed titrating lists (category A1/A2/A3 or B1/B2/B3, seen in at least 5 ILM or LIF in the 54 macro groups" >> $toto

   if ($col == bad) echo "In each run, count the titrating gene contradictory to the union of the confirmed set A1 A2 A3 B1 B2 B3" >> $toto
    \rm toto
foreach sea (3)
  foreach rm (50)
    foreach seaWall (2 4)
      foreach seaUU (u nu)
        if (-e  tmp/GENEINDEX/Results/$MAGIC.seaW$seaWall.rm$rm.$seaUU.titrating_profile.txt) then
          echo sea$seaWall.rm$rm.$seaUU >> toto
        endif
      end
    end
  end
end

    cat toto ZZZZZ MetaDB/$MAGIC/GroupList MetaDB/$MAGIC/RunsList ZZZZZ   tmp/GENEINDEX/Results/$MAGIC.seaW*.rm*.titrating_profile.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){nc++;i2c[nc]=$1;c2i[$1]=nc;next;}if(zz==1){nr++;i2r[nr]=$1;r2i[$1]=nr;next;}r=$2;c=$3;if(0+$4<1)$4=1;ic=c2i[c];ir=r2i[r];irOk[ir]=1;z=0;if(col == "bad")z=0+$7;if(col== "pc")z=$9+$12-$7;if(col== "pp")z=100*($9+$12-$7)/$4;if(col== "gc")z=$9+$10+$12+$13-$7;if(col== "gp")z=100*($9+$10+$12+$13-$7)/$4;if(col== "tc")z=$5+$6-$7;if(col== "tp")z=100*($5+$6-$7)/$4;ncr[ic,ir]=z;}END{for(ic=1;ic<=nc;ic++)printf("\t%s",i2c[ic]);for(ir=1;ir<=nr;ir++){if(irOk[ir]){printf("\n%s",i2r[ir]);for(ic=1;ic<=nc;ic++)printf("\t%.1f",ncr[ic,ir]);}}printf("\n");}' col=$col >> $toto
end

# enif of MAGIC == SEQC_Main
endif

foreach sea (2 7)
  foreach lm (50 500)
     cat tmp/GENEINDEX/Results/$MAGIC.sea$sea.rm$lm.av.nu.TitratingGenes.txt | gawk -F '\t' '/^0/{zz=1;}{if(zz<1)next;for(i=4;i<=11;i++){gg[i,$i]=1;ggg[$i]++;}}END{for(g in ggg){if(ggg[g]==4)good[g]=1;a=0;b=0;for(i=4;i<=10;i+=2){a+=gg[i,g];b+=gg[i+1,g];}if(a>0 && b>0)bad[g]=1;}for(g in bad)print g;}' > bad.$sea.$lm.list
     cat tmp/GENEINDEX/Results/$MAGIC.sea$sea.rm$lm.av.nu.TitratingGenes.txt | gawk -F '\t' '/^0/{zz=1;}{if(zz<1)next;for(i=4;i<=11;i++){gg[i,$i]=1;ggg[$i]++;}}END{for(g in ggg){if(ggg[g]==4)good[g]=1;a=0;b=0;for(i=4;i<=10;i+=2){a+=gg[i,g];b+=gg[i+1,g];}if(a>0 && b>0)bad[g]=1;}for(g in good)print g;}' > good.$sea.$lm.list
  end
end
gzip bad.*.list

gunzip -c bad.2.50.list.gz ZZZZZ.gz tmp/GENEINDEX/Results/$MAGIC.sea2.rm50.av.nu.ace.gz | gawk '/^ZZZZZ/{zz=1;next;}{if(zz<1){bad[$1]=1;next;}}/^Genes/{next;}/^Gene/{g=$2;gsub(/\"/,"",g);a=-9999;ok=0;if(bad[g])ok=1;next;}{if(ok==0)next;}/^Group_nU/{z=$2;gsub(/\"/,"",z);if(z=="D_I_BGI_A_5libs_80runs"){a=$3;aa[int($3+.5)]++;}if(z=="D_I_BGI_B_5libs_80runs"){b=$3;bb[int($3+.5)]++;if(a>-9999)dd[int(1000+ 10*(b-a))]++;}}END{for(x in dd)x1=0+x;x2=x1;for(x in dd){if(x+0<x1)x1=x+0;if(x+0>x2)x2=x+0;}for(x=x1;x<=x2;x++)printf("diff\t%.1f\t%d\n",(x-1000)/10,dd[x]);for(x in aa)x1=x+0;x2=x1;for(x in aa){if(x+0<x1)x1=x+0;if(x+0>x2)x2=x+0;}for(x=x1;x<=x2;x++)printf("A\t%d\t%d\n",x,aa[x]);for(x in bb)x1=x+0;x2=x1;for(x in bb){if(x+0<x1)x1=x.0;if(x+0>x2)x2=x+0;}for(x=x1;x<=x2;x++)printf("B\t%d\t%d\n",x,bb[x]);}' > RESULTS/baddy.histo.txt


gunzip -c good.2.50.list.gz ZZZZZ.gz tmp/GENEINDEX/Results/$MAGIC.sea2.rm50.av.nu.ace.gz | gawk '/^ZZZZZ/{zz=1;next;}{if(zz<1){bad[$1]=1;next;}}/^Genes/{next;}/^Gene/{g=$2;gsub(/\"/,"",g);a=-9999;ok=0;if(bad[g])ok=1;next;}{if(ok==0)next;}/^Group_nU/{z=$2;gsub(/\"/,"",z);if(z=="D_I_BGI_A_5libs_80runs"){a=$3;if(a<2)a=2;aa[int(a+.5)]++;}if(z=="D_I_BGI_B_5libs_80runs"){b=$3;if(b<2)b=2;bb[int(b+.5)]++;if(a>-9999)dd[int(1000+10*(b-a))]++;}}END{for(x in dd)x1=0+x;x2=x1;for(x in dd){if(x+0<x1)x1=x+0;if(x+0>x2)x2=x+0;}for(x=x1;x<=x2;x++)printf("diff\t%.1f\t%d\n",(x-1000)/10,dd[x]);for(x in aa)x1=x+0;x2=x1;for(x in aa){if(x+0<x1)x1=x+0;if(x+0>x2)x2=x+0;}for(x=x1;x<=x2;x++)printf("A\t%d\t%d\n",x,aa[x]);for(x in bb)x1=x+0;x2=x1;for(x in bb){if(x+0<x1)x1=x+0;if(x+0>x2)x2=x+0;}for(x=x1;x<=x2;x++)printf("B\t%d\t%d\n",x,bb[x]);}'> RESULTS/gooddy.histo.txt


cat toto88.sampleClassificationBySignature.txt | gawk -F '\t' '/^#Sample/{for(i=3;i<=NF;i++)s[i]=$i;next;}/^#Run\t/{for(i=3;i<=NF;i++)r[i]=$i;next;}/Ratio NB_female\/NB_male/{for(i=3;i<=NF;i++)ff[i]=$i;next;}/NB_Female_measured\/NB_Male_measured/{for(i=3;i<=NF;i++)ff2[i]=$i;nf=NF;next;}END{for(i=3;i<=nf;i++){uu="";if(index(s[i],"_F_")>1)uu="###";if(index(s[i],"_M_")>1)uu="---";printf("%s\t%s\t%.2f\t%.2f\t%s\n",r[i],s[i],ff[i],ff2[i],uu);}}' | sort -k 3n

cat toto88.sampleClassificationBySignature.txt | gawk -F '\t' '/^#Sample/{for(i=3;i<=NF;i++)s[i]=$i;next;}/^#Run\t/{for(i=3;i<=NF;i++)r[i]=$i;next;}/f_NB_191/{for(i=3;i<=NF;i++)ff[i]=$i;nf=NF;next;}/MYCN_low/{for(i=3;i<=NF;i++)ff2[i]=$i;nf=NF;next;}/MYCN2_low/{for(i=3;i<=NF;i++)ff3[i]=$i;nf=NF;next;}END{for(i=3;i<=nf;i++){uu="";if(substr(s[i],1,1)=="f")uu="###";if(substr(s[i],1,1)=="u")uu="---";printf("%s\t%s\t%.2f\t%.2f\t%.2f\t%s\n",r[i],s[i],ff[i],ff2[i],ff3[i],uu);}}' | sort -k 3n


cat tmp/GENEINDEX/Results/NB.av.u.sampleClassificationBySignature.txt | gawk -F '\t' '/^#Sample/{for(i=3;i<=NF;i++)s[i]=$i;next;}/^#Run\t/{for(i=3;i<=NF;i++)r[i]=$i;next;}/Ratio NB_female\/NB_male/{for(i=3;i<=NF;i++)ff[i]=$i;next;}/NB_Female_measured\/NB_Male_measured/{for(i=3;i<=NF;i++)ff2[i]=$i;nf=NF;next;}END{for(i=3;i<=nf;i++){uu="";if(index(s[i],"_F_")>1)uu="###";if(index(s[i],"_M_")>1)uu="---";printf("%s\t%s\t%.2f\t%.2f\t%s\n",r[i],s[i],ff[i],ff2[i],uu);}}' | sort -k 3n | gawk '{if($5=="---")f++;if($5=="###"){u++;s+=f}n++;if (n<200)a[$5]++;if(n>330)b[$5]++;aa[$5]=1;}END{beta=s/(u*f);printf("s=%d uf=%d u=%d f=%d AUC=%.2f\t",s,u*f,u,f,100*beta);for(k in aa)print k,a[k],b[k]}'

# quantify the different constructions 

# add the 3 outcome
set target=HHong
set target=introns
set target=RefSeq
set target=av

set parity=2
cat runs_MYCN_over_16.5.list ZZZZZ tmp/GENEINDEX/Results/$MAGIC.$target.$uu.sampleClassificationBySignature.txt |  gawk -F '\t' -f scripts/NB.roc.awk justLowMYCN=0 parity=$parity >  tmp/GENEINDEX/Results/$MAGIC.$target.u.outcome.$parity.txt
\cp  tmp/GENEINDEX/Results/$MAGIC.$target.u.outcome.$parity.txt RESULTS/Expression

# RESULTS/Expression.Mix4.nu.min10.max50.NA1/$MAGIC.$target.u.outcome.$parity.txt
# construct a ROC curve
scripts/NB.roc.tcsh $target $parity .

cat _roc.$parity | grep uf
head -4 _km.$parity
\cp _roc.1           _roc.$target.D13a.nu.min10.max50.NA1.txt
\cp _roc.1   RESULTS/_roc.$target.D13a.nu.min10.max50.NA1.txt

# a hack
# construct a classifier for Devan
set tutu=RESULTS/HISTOS/Devan.txt
echo "# " > $tutu
date >> $tutu
echo -n "# HR classifier derived from Devan's list of 19 genes :\n#" >> $tutu
cat RESULTS/HISTOS/DevanList19genesHR.txt | sort |  gawk '{if($1)printf (" %s",$1);}' >> $tutu
echo "\n# NB\tHigh-risk\tClassifier\tBinary_classifier\tSurvived/Died\tFavorable/Unfavorable\tEvent-free" >> $tutu
# cat  RESULTS/Expression.Mix4.nu.min10.max50.NA1/NB.av.u.outcome.txt
set ii=8
cat tmp/GENEINDEX/Results/$MAGIC.$target.u.outcome.txt | sort -k $ii'n' | gawk -F '\t' '/^#/{next}{h="";s=$3;if(index(s,"_HR")>1)h="HR";x=$ii-z;b=1;if(x<0)b=-1;s=$3;split(s,aa,"_");ph=substr(aa[1],length(aa[1]));ef=substr(aa[1],length(aa[1]-1),1);k=substr($2,length($2));f="";if(k%2 ==1){ff=substr(s,1,1);f="";if(ff=="f")f=+.7;if(ff=="u")f=-.7;}printf("%s\t%s\t%.3f\t%s\t%.1f\t%s\t%.2f\t%s\n",$2,h,-x,-b,0.5-ph,f,-(ef-.5)/2,s);}' ii=$ii z=2.560 | grep NB >> $tutu

# 1.425 for Devan #  2.560 for acevHR # -5.925  OS ii=7
# cat ~/ACEVIEWHELP/exonJunction.txt | gawk -F '\t' '/^#/{print;next;}{if(length($1)==0){print;next;}}{a=$5-$4-$3-$2;b=0;if($6+$7 > 0)b=a*100/($6+$7);printf("%s\t%d",$0,a) ;if(a<10000)printf("\t");else printf("\t%.2f",b);printf("\n");next;}' >  ~/ACEVIEWHELP/exonJunction.fixed.txt

# bin/geneindex -deepGene tmp/GENEINDEX/NB.av.u.ace -u -mask tmp/GENEINDEX/NB.av.nu.mask  $compare $refG -runList tmp/GENEINDEX/NB.GroupsRunsListSorted -runAce tmp/GENEINDEX/NB.info.ace  -o tmp/GENEINDEX/Results/NB.av.u -gzo -pA -method Gene_AceView08 -seedGene MYCN -seaLevel 3  -removeLimit 50    -export aitsvz -stableGenes TARGET/Targets/$species.av.stable_genes.txt

# Index1 > 8 no rule on ok
# RefSeq  Outcome u/f     Classifier: All f/u     s=139 uf=17460 u=97 f=180 AUC=0.80
# RefSeq  Outcome u/f     Classifier: Training f/u        s=523 uf=17460 u=97 f=180 AUC=3.00
# RefSeq  Outcome u/f     Classifier: Training Event free s=10307 uf=17280 u=96 f=180 AUC=59.65
# RefSeq  Outcome u/f     Classifier: Training Overall survival   s=201 uf=17460 u=97 f=180 AUC=1.15
# RefSeq  Outcome u/f     Classifier: MYCN april  s=14542 uf=17280 u=96 f=180 AUC=84.16
# RefSeq  Outcome u/f     Classifier: MYCN april iterated s=16674 uf=17280 u=96 f=180 AUC=96.49
# RefSeq  Outcome u/f     Classifier: MYCN april 2 iter   s=16244 uf=17280 u=96 f=180 AUC=94.00
# RefSeq  Outcome u/f     Classifier: MYCN seeded s=15331 uf=17280 u=96 f=180 AUC=88.72
# RefSeq  Outcome u/f     Classifier: MYCN seeded iter    s=16490 uf=17280 u=96 f=180 AUC=95.43
# RefSeq  Outcome u/f     Classifier: MYCN seeded 2 iter  s=8007 uf=17460 u=97 f=180 AUC=45.86
# RefSeq  Outcome EF      Classifier: All f/u     s=11189 uf=59148 u=186 f=318 AUC=18.92
# RefSeq  Outcome EF      Classifier: Training f/u        s=12974 uf=59148 u=186 f=318 AUC=21.93
# RefSeq  Outcome EF      Classifier: Training Event free s=29656 uf=58830 u=185 f=318 AUC=50.41
# RefSeq  Outcome EF      Classifier: Training Overall survival   s=12221 uf=59148 u=186 f=318 AUC=20.66
# RefSeq  Outcome EF      Classifier: MYCN april  s=40329 uf=58830 u=185 f=318 AUC=68.55
# RefSeq  Outcome EF      Classifier: MYCN april iterated s=45753 uf=58830 u=185 f=318 AUC=77.77
# RefSeq  Outcome EF      Classifier: MYCN april 2 iter   s=44355 uf=58830 u=185 f=318 AUC=75.40
# RefSeq  Outcome EF      Classifier: MYCN seeded s=42953 uf=59148 u=186 f=318 AUC=72.62
# RefSeq  Outcome EF      Classifier: MYCN seeded iter    s=45843 uf=58830 u=185 f=318 AUC=77.92
# RefSeq  Outcome EF      Classifier: MYCN seeded 2 iter  s=27227 uf=59148 u=186 f=318 AUC=46.03
# RefSeq  Outcome OS      Classifier: All f/u     s=4926 uf=44067 u=111 f=397 AUC=11.18
# RefSeq  Outcome OS      Classifier: Training f/u        s=5849 uf=44067 u=111 f=397 AUC=13.27
# RefSeq  Outcome OS      Classifier: Training Event free s=26477 uf=43670 u=110 f=397 AUC=60.63
# RefSeq  Outcome OS      Classifier: Training Overall survival   s=5176 uf=44067 u=111 f=397 AUC=11.75
# RefSeq  Outcome OS      Classifier: MYCN april  s=33391 uf=43670 u=110 f=397 AUC=76.46
# RefSeq  Outcome OS      Classifier: MYCN april iterated s=37571 uf=43670 u=110 f=397 AUC=86.03
# RefSeq  Outcome OS      Classifier: MYCN april 2 iter   s=36734 uf=43670 u=110 f=397 AUC=84.12
# RefSeq  Outcome OS      Classifier: MYCN seeded s=34884 uf=44067 u=111 f=397 AUC=79.16
# RefSeq  Outcome OS      Classifier: MYCN seeded iter    s=37295 uf=43670 u=110 f=397 AUC=85.40
# RefSeq  Outcome OS      Classifier: MYCN seeded 2 iter  s=19451 uf=44067 u=111 f=397 AUC=44.14
# av      Outcome u/f     Classifier: All f/u     s=149 uf=17460 u=97 f=180 AUC=0.85
# av      Outcome u/f     Classifier: Training f/u        s=170 uf=17460 u=97 f=180 AUC=0.97
# av      Outcome u/f     Classifier: Training Event free s=4498 uf=17280 u=96 f=180 AUC=26.03
# av      Outcome u/f     Classifier: Training Overall survival   s=142 uf=17460 u=97 f=180 AUC=0.81
# av      Outcome u/f     Classifier: MYCN april  s=16428 uf=17280 u=96 f=180 AUC=95.07
# av      Outcome u/f     Classifier: MYCN april iterated s=16884 uf=17280 u=96 f=180 AUC=97.71
# av      Outcome u/f     Classifier: MYCN april 2 iter   s=16742 uf=17280 u=96 f=180 AUC=96.89
# av      Outcome u/f     Classifier: MYCN seeded s=16399 uf=17280 u=96 f=180 AUC=94.90
# av      Outcome u/f     Classifier: MYCN seeded iter    s=16041 uf=17280 u=96 f=180 AUC=92.83
# av      Outcome u/f     Classifier: MYCN seeded 2 iter  s=8007 uf=17460 u=97 f=180 AUC=45.86
# av      Outcome EF      Classifier: All f/u     s=12198 uf=59148 u=186 f=318 AUC=20.62
# av      Outcome EF      Classifier: Training f/u        s=11707 uf=59148 u=186 f=318 AUC=19.79
# av      Outcome EF      Classifier: Training Event free s=20702 uf=59148 u=186 f=318 AUC=35.00
# av      Outcome EF      Classifier: Training Overall survival   s=11963 uf=59148 u=186 f=318 AUC=20.23
# av      Outcome EF      Classifier: MYCN april  s=44794 uf=58830 u=185 f=318 AUC=76.14
# av      Outcome EF      Classifier: MYCN april iterated s=46428 uf=58830 u=185 f=318 AUC=78.92
# av      Outcome EF      Classifier: MYCN april 2 iter   s=45856 uf=58830 u=185 f=318 AUC=77.95
# av      Outcome EF      Classifier: MYCN seeded s=44861 uf=58830 u=185 f=318 AUC=76.26
# av      Outcome EF      Classifier: MYCN seeded iter    s=43723 uf=58830 u=185 f=318 AUC=74.32
# av      Outcome EF      Classifier: MYCN seeded 2 iter  s=27227 uf=59148 u=186 f=318 AUC=46.03
# av      Outcome OS      Classifier: All f/u     s=5036 uf=44067 u=111 f=397 AUC=11.43
# av      Outcome OS      Classifier: Training f/u        s=4893 uf=44067 u=111 f=397 AUC=11.10
# av      Outcome OS      Classifier: Training Event free s=14874 uf=44067 u=111 f=397 AUC=33.75
# av      Outcome OS      Classifier: Training Overall survival   s=4696 uf=44067 u=111 f=397 AUC=10.66
# av      Outcome OS      Classifier: MYCN april  s=37007 uf=43670 u=110 f=397 AUC=84.74
# av      Outcome OS      Classifier: MYCN april iterated s=37987 uf=43670 u=110 f=397 AUC=86.99
# av      Outcome OS      Classifier: MYCN april 2 iter   s=37756 uf=43670 u=110 f=397 AUC=86.46
# av      Outcome OS      Classifier: MYCN seeded s=37103 uf=43670 u=110 f=397 AUC=84.96
# av      Outcome OS      Classifier: MYCN seeded iter    s=36588 uf=43670 u=110 f=397 AUC=83.78
# av      Outcome OS      Classifier: MYCN seeded 2 iter  s=19451 uf=44067 u=111 f=397 AUC=44.14

# sig itere a .7
goto plusBas55
av      Outcome u/f     Classifier: All f/u     s=187 uf=17472 u=96 f=182 AUC=1.07
av      Outcome u/f     Classifier: Training f/u        s=1202 uf=17472 u=96 f=182 AUC=6.88
av      Outcome u/f     Classifier: Training Event free s=276 uf=17472 u=96 f=182 AUC=1.58
av      Outcome u/f     Classifier: Training Overall survival   s=17148 uf=17290 u=95 f=182 AUC=99.18
av      Outcome u/f     Classifier: High Risk   s=6107 uf=17472 u=96 f=182 AUC=34.95
av      Outcome u/f     Classifier: MYCN april  s=14939 uf=17290 u=95 f=182 AUC=86.40
av      Outcome u/f     Classifier: MYCN april iterated s=16806 uf=17290 u=95 f=182 AUC=97.20
av      Outcome u/f     Classifier: MYCN april 2 iter   s=16665 uf=17290 u=95 f=182 AUC=96.39
av      Outcome u/f     Classifier: MYCN seeded s=16271 uf=17290 u=95 f=182 AUC=94.11
av      Outcome u/f     Classifier: MYCN seeded iter    s=16620 uf=17290 u=95 f=182 AUC=96.12
av      Outcome u/f     Classifier: MYCN seeded 2 iter  s=8088 uf=17472 u=96 f=182 AUC=46.29
av      Outcome EF      Classifier: All f/u     s=12118 uf=59706 u=186 f=321 AUC=20.30
av      Outcome EF      Classifier: Training f/u        s=16421 uf=59706 u=186 f=321 AUC=27.50
av      Outcome EF      Classifier: Training Event free s=11760 uf=59385 u=185 f=321 AUC=19.80
av      Outcome EF      Classifier: Training Overall survival   s=47375 uf=59385 u=185 f=321 AUC=79.78
av      Outcome EF      Classifier: High Risk   s=26334 uf=59706 u=186 f=321 AUC=44.11
av      Outcome EF      Classifier: MYCN april  s=41257 uf=59385 u=185 f=321 AUC=69.47
av      Outcome EF      Classifier: MYCN april iterated s=46536 uf=59385 u=185 f=321 AUC=78.36
av      Outcome EF      Classifier: MYCN april 2 iter   s=46045 uf=59385 u=185 f=321 AUC=77.54
av      Outcome EF      Classifier: MYCN seeded s=45072 uf=59385 u=185 f=321 AUC=75.90
av      Outcome EF      Classifier: MYCN seeded iter    s=46931 uf=59385 u=185 f=321 AUC=79.03
av      Outcome EF      Classifier: MYCN seeded 2 iter  s=27626 uf=59706 u=186 f=321 AUC=46.27
av      Outcome OS      Classifier: All f/u     s=4862 uf=44289 u=111 f=399 AUC=10.98
av      Outcome OS      Classifier: Training f/u        s=7789 uf=44289 u=111 f=399 AUC=17.59
av      Outcome OS      Classifier: Training Event free s=5023 uf=44289 u=111 f=399 AUC=11.34
av      Outcome OS      Classifier: Training Overall survival   s=39209 uf=43890 u=110 f=399 AUC=89.33
av      Outcome OS      Classifier: High Risk   s=20645 uf=44289 u=111 f=399 AUC=46.61
av      Outcome OS      Classifier: MYCN april  s=34283 uf=44289 u=111 f=399 AUC=77.41
av      Outcome OS      Classifier: MYCN april iterated s=38102 uf=44289 u=111 f=399 AUC=86.03
av      Outcome OS      Classifier: MYCN april 2 iter   s=37847 uf=43890 u=110 f=399 AUC=86.23
av      Outcome OS      Classifier: MYCN seeded s=37169 uf=43890 u=110 f=399 AUC=84.69
av      Outcome OS      Classifier: MYCN seeded iter    s=37797 uf=43890 u=110 f=399 AUC=86.12
av      Outcome OS      Classifier: MYCN seeded 2 iter  s=19587 uf=44289 u=111 f=399 AUC=44.23

plusBas55

# construct the PCA clustering

gunzip -c RESULTS/Expression/unique/av/$MAGIC.av.u.expression_index.txt.gz |  gawk -F '\t' '{ok=0;}/^#Run\t/{ok=1;nr++;if(nr>1)next;}/^#/{if(ok==0)next;}/^$/{next;}{if(ok==0)printf("%s",$1);}{for(i=12;i<=121+12-1;i++)printf("\t%f",0+$i);printf("\n");next;}' > RESULTS/Expression/unique/av/$MAGIC.av.u.expression_index.R.txt

R << EOF
genes = read.table("RESULTS/Expression/unique/av/$MAGIC.av.u.expression_index.R.txt")
gg=as.matrix(genes)
cc=cor(gg)
# heatmap takes 5 minutes to draw
pdf("RESULTS/Expression/unique/av/$MAGIC.av.u.expression_index.PCA.correl.pdf") ;
hh=heatmap(cc)
write(colnames(gg)[hh$rowInd],"RESULTS/Expression/unique/av/$MAGIC.av.u.expression_index.PCA.correl.Legend.txt")
# pp=prcomp(genes,center=TRUE,scale=FALSE)
quit()
EOF
gv RESULTS/Expression/unique/av/$MAGIC.av.u.expression_index.PCA.correl.pdf

### special list for rat TGx

\rm titi
foreach ff (`ls  tmp/GENEINDEX/Results/Liver.av.u.diff.*.txt`)
  set n=`echo $ff | gawk '/4Heal/{n=1;}/AnyHealthy/{n=1;}/lb27/{n=1;}/13rats/{n=1;}/Atypic/{k=1;}/Patho/{k=1;}/Control/{i=index($1,"Control");j=index(substr($1,i+5),"Control");if(j>0)k=1;}{if(k>0)n=0;print 0+n;}'`
  if ($n == 0) continue
  cat $ff | gawk -F '\t' '/^# Genes/{title=substr($1,46);next;}/^#/{next;}{if(0 && $2=="")next;ix1=$11;dIndex=$13;sig=$17;z=ix1+dIndex+2*sig;printf("%d\t%s\t%s\t%s\n",int(100*z),$1,$9,title);}' >> titi
end
cat titi | sort -k 1nr | head -5000 | gawk -F '\t' '{g=$2;t[g]=$3;s=$4;ss[s]=1;z[g,s]=1;next;}END{printf("#Gene\tTitle");for(s in ss)printf("\t%s",s);for(g in t){printf("\n%s\t%s",g,t[g]);for(s in ss)printf("\t%s",z[g,s]);}printf("\n");}' > titi2

head -1 titi2 >  RESULTS/Expression/unique/av/$MAGIC.top_differential_genes.txt
tail -n +2 titi2 | sort >>  RESULTS/Expression/unique/av/$MAGIC.top_differential_genes.txt


head -1 titi2 >  RESULTS/Expression/unique/av/$MAGIC.top_differential_genes_with_geneid.txt
tail -n +2 titi2 | sort >>  RESULTS/Expression/unique/av/$MAGIC.top_differential_genes_with_geneid.txt


###

foreach target (av)
  set targetBeau=$target
  if ($target == av) set targetBeau=AceView

  foreach uu (u)
    cat  MetaDB/$MAGIC/RunsListSorted | sed -e 's/\"//g'  | gawk '/^Rrn/{printf("%s\n",$1);}'  > toto.rna_affy.runList
    cat  MetaDB/$MAGIC/RunsListSorted | sed -e 's/\"//g'  | gawk '/^Rrn/{printf("%sPB\n",$1);}'  >> toto.rna_affy.runList
    cat  MetaDB/$MAGIC/RunsListSorted | sed -e 's/\"//g'  | gawk '/^Rrn/{next;}{print}'  >> toto.rna_affy.runList
    cat MicroArray/Hits/Affy.Rat230_2.$target.probeset2gene.unique.ace tmp/GENEINDEX/$MAGIC.$target.$GM.info.ace  > toto.rna_affy.info.ace  
    gunzip -c tmp/GENEINDEX/$target.$uu.geneSupport.ace.gz   > toto.rna_affy.genesupport.ace
    cat MicroArray/PierreBushel/PB.training.ace2 >> toto.rna_affy.genesupport.ace

    bin/geneindex -deepGene toto.rna_affy.genesupport.ace -runList toto.rna_affy.runList -runAce  toto.rna_affy.info.ace  $refG -o tmp/GENEINDEX/Results/$MAGIC.$target.RNA_seq_Affy -gzo -MA -export i -method "Affymetrix_Pierre_Bushel"
    \cp tmp/GENEINDEX/Results/Liver.av.RNA_seq_Affy.expression_index.txt.gz RESULTS/Expression/unique/av/
    \cp tmp/GENEINDEX/Results/Liver.av.RNA_seq_Affy.sampleClassificationByCovariance.txt  RESULTS/Expression/unique/av/

    gunzip -c MicroArray/ProbeCount/Affy.Rat230_2.probeset.ace.gz  > toto.affyPB_PS.genesupport.ace
    cat MicroArray/PierreBushel/PB.training.ace2 >> toto.affyPB_PS.genesupport.ace
    bin/geneindex -deepGene toto.affyPB_PS.genesupport.ace -runList toto.rna_affy.runList -runAce  toto.rna_affy.info.ace  $refG -o tmp/GENEINDEX/Results/$MAGIC.$target.AffyPB_PS -gzo -MA 60 -export i 
    \cp tmp/GENEINDEX/Results/Liver.av.AffyPB_PS.expression_index.txt.gz RESULTS/Expression/unique/av/
    \cp tmp/GENEINDEX/Results/Liver.av.AffyPB_PS.sampleClassificationByCovariance.txt  RESULTS/Expression/unique/av/


    cat  MetaDB/$MAGIC/RunsListSorted | sed -e 's/\"//g' | gawk '/^Rrn/{printf("%s\n",$1);}'  > toto.rna_affy.runList
    cat  MetaDB/$MAGIC/RunsListSorted | sed -e 's/\"//g'  | gawk '/^Rrn/{printf("%sPS\n",$1);}'  >> toto.rna_affy.runList
    cat  MetaDB/$MAGIC/RunsListSorted | sed -e 's/\"//g'  | gawk '/^Rrn/{next;}{print}'  >> toto.rna_affy.runList
    gunzip -c tmp/GENEINDEX/$target.$uu.geneSupport.ace.gz   > toto.rna_gene_probeset.genesupport.ace
    gunzip -c MicroArray/ProbeCount/Affy.Rat230_2.probeset.ace.gz | gawk '/^Run_u/{printf("Run_u %sPS ", $2);for(i=3;i<=NF;i++)printf(" %s",$i);printf("\n");next;}{print}' >> toto.rna_gene_probeset.genesupport.ace
     bin/geneindex -deepGene  toto.rna_gene_probeset.genesupport.ace  -runList toto.rna_affy.runList -runAce  toto.rna_affy.info.ace  $refG -o tmp/GENEINDEX/Results/$MAGIC.$target.RNA_gene_probeset -gzo -MA -export i 
    \cp tmp/GENEINDEX/Results/Liver.av.RNA_gene_probeset.expression_index.txt.gz RESULTS/Expression/unique/av/
    \cp tmp/GENEINDEX/Results/Liver.av.RNA_gene_probeset.sampleClassificationByCovariance.txt RESULTS/Expression/unique/av/

  end
end

goto phaseLoop
######################################################
## March 18, 2012, export to the TGX sharepoint of Pierre Bushel

########## Individual Probes
set toto1=RESULTS.$MAGIC/MicroArray/Affy.Rat230_2.probe.expression_index.txt
set toto=RESULTS.$MAGIC/MicroArray/TGx.NCBI.120419.File1.RNAseq_expression_index_at_Affymetrix_individual_probe_positions.txt
echo -n "# " > $toto
date >> $toto
echo "#RNA-seq expression index measured at the location of Affymetrix Rat230.2 probes uniquely mapped on the main chromosomes." >> $toto
echo "#The RNA-seq index of expression is directly comparable to a microarray logarithmic luminosity." >> $toto
echo "#The Affymetrix probe sequences were mapped to the genome allowing at most one mismatch." >> $toto
echo "#The log base 2 of the RNA-seq coverage of the genome at the location of each uniquely mapped probe is measured." >> $toto
echo "#As in microarray experiments, each run (i.e. a complete column) could be shifted to fix the median or average of the column." >> $toto
echo "#Absent values were considered too low to be reliably reported, once taking into account the noise in the particular RNA-seq experiment." >> $toto
gunzip -c  $toto1.gz | head -50 | gawk -F '\t' '/^#RunId/{n++;if(n==1)print;next;}' >> $toto
echo "#Probe identifier" >> $toto
gunzip -c  $toto1.gz | gawk -F '\t' '/^Outside/{next}/^#/{next;}/^$/{next;}{i=index($2,"Rn");if(i>=1)next;print}' | sort >> $toto
\rm $toto.gz
gzip $toto

 gunzip -c $toto.gz | head -100 > RESULTS.$MAGIC/test1.txt             
gunzip -c $toto.gz | head -10000 | tail -100 >> RESULTS.$MAGIC/test1.txt   
gunzip -c $toto.gz |  tail -100 >> RESULTS.$MAGIC/test1.txt   

# nb probes in file
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^Rat/{print $1;next;}' | sort -u | wc
# from  probeset in file
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^Rat/{j=index($1,"_at");print substr($1,1,j);next;}' | sort -u | wc
# nb probes with some expressiom
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}{for (i=7;i<=NF;i++){if($i>3){print $1;next;}}}' | sort -u | wc
# nb probes with some expressiom
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}{for (i=7;i<=NF;i++){if($i>3){j=index($1,"_at");print substr($1,1,j+2);next;}}}' | sort -u | wc

########## ProbeSets
set toto1=RESULTS.$MAGIC/MicroArray/Affy.Rat230_2.probeset.expression_index.txt
set toto=RESULTS.$MAGIC/MicroArray/TGx.NCBI.120330.File2.RNAseq_expression_index_at_Affymetrix_probeset_positions.txt

echo -n "# " > $toto
date >> $toto
echo "# Affymetrix level of expression of 31099 probesets normalized by Pierre Bushel compared to" >> $toto
echo "# RNA-seq index of expression of 28976 microarray probe sets mapped uniquely on the main chromosomes, excluding the floating contigs." >> $toto
echo "# The RNA-seq index of expression is directly comparable to a microarray logarithmic luminosity." >> $toto
echo "# The Affymetrix probe sequences were remapped to the genome allowing at most one mismatch." >> $toto
echo "# Probes mapping at several locations were removed, in effect reducing the number of probes from some probesets." >> $toto
echo "# The RNA-seq coverage of the genome at the location of each uniquely mapped probe is measured," >> $toto
echo "# and the average of these values is attributed to each probeset. An interesting refinement would be to discard" >> $toto
echo "# the individual probes whose differential expression across all runs does not correlate well with the rest of" >> $toto
echo "# the probe set. This would discard dead probes and cross hybridizing or ill designed probes." >> $toto
echo "# As in microarray experiments, each run (i.e. a complete column) could be shifted to fix the median or average of the column." >> $toto
echo "# Absent values were considered too low to be reliably reported, once taking into account the noise in the particular RNA-seq experiment." >> $toto



gunzip -c $toto1.gz | head -30 | gawk '/^#RunId/{n++;if(n==1)print;}' | sed -e 's/Measured object/Method/' > $toto.3

cat $toto.3 | gawk -F '\t' '{printf("%s",$1);for (i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' >> $toto
gunzip -c  MicroArray/PierreBushel/PB.training.renamed.expression_index.txt.gz | head -15 | gawk -F '\t' '/^#Sample/{n++;if(n==1)print}' | sed -e 's/Sample/Affymetrix array/' -e 's/Measured object/Method/' >>  $toto
echo "Probeset\tMethod" >> $toto


gunzip -c  MicroArray/PierreBushel/PB.training.renamed.expression_index.txt.gz | gawk -F '\t' '/^Rat/{printf("%s\tAffymetrix_Pierre_Bushel",$1);for (i=3;i<=NF;i++)printf("\t%s",$i);printf("\n");}'  > $toto.1
gunzip -c  $toto1.gz   | gawk -F '\t' '/^Rat/{printf("%s\tProbeset_RNAseq",$1);for (i=3;i<=NF;i++)printf("\t%s",$i);printf("\n");}' >> $toto.1

cat $toto.1  | sort >> $toto
\rm $toto.gz $toto1 $toto.3
gzip $toto

 gunzip -c $toto.gz | head -100 > RESULTS.$MAGIC/test2.txt             
gunzip -c $toto.gz | head -10000 | tail -100 >> RESULTS.$MAGIC/test2.txt   
gunzip -c $toto.gz |  tail -100 >> RESULTS.$MAGIC/test2.txt   

# nb probeset in file
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/RNA/{print $1;next;}' | sort -u | wc
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/Pierre/{print $1;next;}' | sort -u | wc
# nb probes with some expressiom
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/RNA/{for (i=7;i<=NF;i++){if($i>3){print $1;next;}}}' | sort -u | wc

## compute the correlation coefficient
set toto=/home/mieg/ACEVIEWHELP/RESULTS.Rat_liver_apr4/MicroArray/TGx.NCBI.120330.File2.RNAseq_expression_index_at_Affymetrix_probeset_positions

set median=0
set type1=Affymetrix
set type2=RNAseq
set plage=0

if ($plage==1) then
  set toto2=$toto.correl.$type1.$type2.median$median.plageMinMax.txt
  set toto3=$toto.correl.$type1.$type2.median$median.plageMinMax.txt3
else
  set toto2=$toto.correl.$type1.$type2.median$median.plageMinAny.txt
  set toto3=$toto.correl.$type1.$type2.median$median.plageMinAny.txt3
endif

echo  -n "# " > $toto2
date >> $toto2
gunzip -c $toto.txt.gz | gawk -F '\t' '/^#RunId/{printf("\t\t");print;}' >> $toto2
echo  -n "# " > $toto3

foreach min (0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)
  foreach max (1000)
    if ($plage == 1) then
      @ max = $min + 1
    endif
    set toto1=$toto.correl.$type1.$type2.median$median.min$min.max$max.txt
    echo  -n "# " > $toto1
    date >> $toto1
    echo "# Correlation between Affymetrix MASS5 measures and RNA-seq logarithmic index of coverage at the location of the uniquely mapped probe sets" >> $toto1
    echo "# Remove index below $min, renormalize each gene to 10 yes/no=$median" >> $toto1
    gunzip -c $toto.txt.gz | gawk -F '\t' '/^#RunId/{print;for(i=12;i<=NF;i++)rid[i]=$i;nf=NF;next;}/^#/{next;}{split($1,aa,"(");$1=aa[1];}{n=0;m=0;for(i=12;i<=NF;i++){if($i+0>0){m+=0+$i;n++}}if(n==0)next;m=m/n;if(median==1){for(i=12;i<=NF;i++)if($i>0)$i=$i-m+10;}}/'$type1'/{p=$1;for(i=12;i<=NF;i++)px[p,i]=0+$i;pp[p]=1;next;}/'$type2'/{if(m<0+qxMin || m>=0+qxMax)next;q=$1;for(i=12;i<=NF;i++)qx[q,i]=0+$i;qq[q]=1;next;}END{np=0;nq=0;npq=0;for(i=12;i<=nf;i++)npqi[i]=0;for(p in pp)np++;for(p in qq){nq++;if(pp[p]){for(i=12;i<=nf;i++){if(( px[p,i]>= 0 && qx[p,i]>= 0)){npq++;npqi[i]++;ZX+=px[p,i];X[i]+=px[p,i];ZX2+=px[p,i]*px[p,i];X2[i]+=px[p,i]*px[p,i];ZY+=qx[p,i];Y[i]+=qx[p,i];ZY2+=qx[p,i]*qx[p,i];Y2[i]+=qx[p,i]*qx[p,i];ZXY+=px[p,i]*qx[p,i];XY[i]+=px[p,i]*qx[p,i];}}}}printf("\nN\t%d",npq);for(i=12;i<=nf;i++){printf("\t%d",npqi[i]);if(npq==0)npq=1;if(npqi[i]==0)npqi[i]=1}printf("\nX\t%.3f",ZX);for(i=12;i<=nf;i++)printf("\t%.3f",X[i]);printf("\nY\t%.3f",ZY);for(i=12;i<=nf;i++)printf("\t%.3f",Y[i]);printf("\nX2\t%.3f",ZX2);for(i=12;i<=nf;i++)printf("\t%.3f",X2[i]);printf("\nY2\t%.3f",ZY2);for(i=12;i<=nf;i++)printf("\t%.3f",Y2[i]);printf("\nXY\t%.3f",ZXY);for(i=12;i<=nf;i++)printf("\t%.3f",XY[i]);x=ZX2/npq-ZX*ZX/(npq*npq);y=ZY2/npq-ZY*ZY/(npq*npq);w=ZXY/npq-ZX*ZY/(npq*npq);if(x*y>0)w=w/sqrt(x*y);printf("\nCorrel\t%.3f",w);for(i=12;i<=nf;i++){x=X2[i]/npqi[i]-X[i]*X[i]/(npqi[i]*npqi[i]);y=Y2[i]/npqi[i]-Y[i]*Y[i]/(npqi[i]*npqi[i]);w=XY[i]/npqi[i]-X[i]*Y[i]/(npqi[i]*npqi[i]);w=w/sqrt(.1+x*y);printf("\t%.3f",w);}printf("\n");}' qxMin=$min qxMax=$max median=$median  >> $toto1

    cat $toto1 | gawk '/^Correl/{printf("%s\t%s\t",mi,mx);print}' mi=$min mx=$max  >> $toto2
    cat $toto1 | gawk '/^N/{printf("%s\t%s\t",mi,mx);print}' mi=$min mx=$max  >> $toto3
  end
end

cat $toto3 >> $toto2

########## AceView

set exportType=expression_index

phaseTGx:
set exportType=expression_index
set exportType=reads_aligned_per_gene

set toto1=RESULTS/Expression/quasi_unique/av/$MAGIC.av.nu.$exportType.txt
set toto=RESULTS/Expression/unique/av/TGx.NCBI.120330.File3.RNAseq_expression_index_per_gene.txt
set toto=RESULTS/Expression/unique/av/TGx.NCBI.120830.File3.RNAseq_expression_index_per_gene.txt
set toto=TGx.NCBI.120830.File3.RNAseq_expression_index_per_gene.txt.gz 
set toto=RESULTS/Expression/quasi_unique/av/TGxSEQC_noNA_GeneExpression_quasi_unique_Magic_20120831.txt
set toto=RESULTS/Expression/quasi_unique/av/TGxSEQC_GeneExpression_ReadCounts_quasi_unique_Magic_20120831.txt

# exportation SEQC_main 2012_06_23 towards SAS SEQC sharepoint 
# aug26 I save all previous results in new directory 
# mv RESULTS/Expression RESULTS/Expression.Exported.2012_06_26
# and I reconstruct the index table with a new code, no NA : interpolate below 4 reads and add a malus using useBonus:ZZ=6
# on sept 4 we reposted the gene index 
# on nov 17 we export the raw counts for Rick Jensent to the share point
set GM=GENE
set exportType=expression_index
set exportType=reads_aligned_per_gene
set toto1=RESULTS/Expression/quasi_unique/av/$MAGIC.av.nu.$exportType.txt
set toto=RESULTS/Expression/quasi_unique/av/$MAGIC.AceViewGenes2010_ExpressionIndex.2012_06_23.txt
set toto=RESULTS/Expression/quasi_unique/av/$MAGIC.AceViewGenes2010_ExpressionIndex.2012_08_27.txt
set toto=RESULTS/Expression/quasi_unique/av/$MAGIC.AceViewGenes2010_ExpressionIndex.2012_09_04.txt
set toto=RESULTS/Expression/quasi_unique/av/$MAGIC.AceViewGenes2010_ReadCounts.2012_09_04.txt
set toto=RESULTS/Expression/quasi_unique/av/$MAGIC.AceViewGenes2010_ReadCounts.NWU_EF2_fixed.MoreGroups.2012_11_23.txt
set toto=RESULTS/Expression/quasi_unique/av/$MAGIC.AceViewGenes2010_ExpressionIndex.NWU_EF2_fixed.MoreGroups.2012_11_23.txt
set toto=RESULTS/Expression/quasi_unique/av/$MAGIC.AceViewGenes2010_ReadCounts.EF_fixed.MoreGroups.2012_12_01.txt
set toto=RESULTS/Expression/quasi_unique/av/$MAGIC.AceViewGenes2010_ExpressionIndex.EF_fixed.MoreGroups.2012_12_01.txt

# NB AceView Genes
phaseTGx.av.GENE:
set exportType=expression_index
set GM=GENE
set target=av
set toto1=RESULTS/Expression/unique/av/$MAGIC.av.GENE.u.$exportType.txt
set toto=RESULTS/Expression/unique/av/NB.NCBI.120429.RNAseq_expression_index_per_gene.txt
set toto=RESULTS/Expression/unique/av/NB.NCBI.120821.RNAseq_expression_index_per_gene.txt
set toto=RESULTS/Expression/quasi_unique/av/NB.NCBI.20120828.GeneExpression.AceView.quasiUnique.txt
set toto=RESULTS/Expression/quasi_unique/av/NB.NCBI.20121127.GeneExpression.AceView.quasiUnique.txt
set toto=RESULTS/Expression/unique/introns/SEQC_NB_MAV_G_log2.20121127.txt
set toto=RESULTS/Expression/unique/av/NB.NCBI.131207.AceView_expression_index_per_gene.unique_pairs.txt
goto TGX2

# NB AceView transcripts
phaseTGx.av.MRNAH:
set exportType=reads_aligned_per_gene
set exportType=expression_index
set GM=MRNAH
set target=av
set toto1=tmp/GENEINDEX/Results/NB.av.$GM.u.$exportType.txt
set toto1=RESULTS/Expression/unique/av/$MAGIC.av.MRNAH.u.$exportType.txt
set toto=RESULTS/Expression/unique/introns/SEQC_NB_MAV_T_log2.20121127.txt
set toto=RESULTS/Expression/unique/introns/SEQC_NB_MAV_Th_log2.20121210.txt
set toto=RESULTS/Expression/unique/introns/SEQC_NB_MAV_Th_log2.20130211.txt
set toto=RESULTS/Expression/unique/introns/SEQC_NB_MAV_Th_$exportType.20130211.txt
set toto=RESULTS/Expression/unique/av/NB.NCBI.131207.AceView_expression_index_per_transcript.hierarchic.unique_pairs.txt
goto TGX2

# NB RefSeq GENE
phaseTGx.RefSeq.GENE:
set exportType=expression_index
set GM=GENE
set target=RefSeq
set toto1=RESULTS/Expression/quasi_unique/RefSeq/$MAGIC.RefSeq.nu.$exportType.txt
set toto=RESULTS/Expression/quasi_unique/RefSeq/NB.NCBI.20120911.GeneExpression.RefSeq.quasiUnique.txt
set toto1=RESULTS/Expression/unique/RefSeq/$MAGIC.RefSeq.GENE.u.$exportType.txt
set toto=RESULTS/Expression/unique/RefSeq/NB.NCBI.131207.RefSeq_expression_index_per_gene.unique_pairs.txt
goto TGX2


# NB RefSeq Transcript
phaseTGx.RefSeq.MRNAH:
set exportType=expression_index
set GM=MRNAH
set target=RefSeq
set toto1=RESULTS/Expression/quasi_unique/RefSeq/$MAGIC.RefSeq.nu.$exportType.txt
set toto=RESULTS/Expression/quasi_unique/RefSeq/NB.NCBI.20120911.GeneExpression.RefSeq.quasiUnique.txt
set toto1=RESULTS/Expression/unique/RefSeq/$MAGIC.RefSeq.MRNAH.u.$exportType.txt
set toto=RESULTS/Expression/unique/RefSeq/NB.NCBI.131207.RefSeq_expression_index_per_transcript.hierarchic.unique_pairs.txt
goto TGX2

# NB EBI GENE
phaseTGx.EBI.GENE:
set exportType=expression_index
set GM=GENE
set target=EBI
set toto1=RESULTS/Expression/unique/EBI/$MAGIC.EBI.GENE.u.$exportType.txt
set toto=RESULTS/Expression/unique/EBI/NB.NCBI.131207.Ensembl_expression_index_per_gene.unique_pairs.txt
goto TGX2

# NB EBI Transcript
phaseTGx.EBI.MRNAH:
set exportType=expression_index
set GM=MRNAH
set target=EBI
set toto1=RESULTS/Expression/unique/$target/$MAGIC.$target.$GM.u.$exportType.txt
set toto=RESULTS/Expression/unique/$target/NB.NCBI.131207.Ensembl_expression_index_per_transcript.hierarchic.unique_pairs.txt
goto TGX2

# NB AGLuK
phaseTGx.AGLuK.MA:
set exportType=expression_index
set GM=MA
set GM_at=MA
set target=AGLuK
set toto1=RESULTS/Expression/unique/$target/$MAGIC.$target.$GM.u.$exportType.txt
set toto=RESULTS/Expression/unique/$target/NB.NCBI.131207.Agilent_expression_index.txt
goto TGX2

# NB RNA_at_AGLuK
phaseTGx.RNA_at_AGLuK.MA:
set exportType=expression_index
set GM=MA
set GM_at=MA_at
set target=RNA_at_AGLuK
set toto1=RESULTS/Expression/unique/$target/$MAGIC.$target.$GM.u.$exportType.txt
set toto=RESULTS/Expression/unique/$target/NB.NCBI.131207.RNA_at_Agilent_expression_index.txt
goto TGX2

# Primates, just simplify the header
 gunzip -c RESULTS/Expression/unique/av/Tissue_export.av.GENE.u.expression_index.txt.gz | gawk -F '\t' '{n++;if(n<=1 || n>=28)print;}' | gzip  > ~/ftp-private/Primates/Tissue.av.GENE.u.expression_index.2014_09_17.txt.gz



########################################
# NB Introns
set exportType=expression_index
set target=Introns
set toto1=RESULTS/Expression/unique/introns/$MAGIC.introns.u.$exportType.txt
set toto=RESULTS/Expression/unique/introns/NB.NCBI.20120911.GeneExpression.Introns.unique.txt
set toto=RESULTS/Expression/unique/introns/SEQC_NB_MAV_J_log2.20121127.txt

gunzip -c $toto1.gz | head -70 | gawk -F '\t' 'function pp() {printf("%s",$1);for(i=2;i<=NF;i++)if(bad[i]<1)printf("\t%s",$i);printf("\n");}/^#Title/{ntitle++;if(ntitle>1)next; bad[2]=1;bad[8]=1;bad[9]=1;for(i=10;i<=NF;i++){if(index($i,"b")>0)bad[i]=1;gsub(/a/,"",$i);gsub(/Title/,"Junction",$1);$i="SEQC_" $i;}pp();}/^#/{next;}{if(NF>10)pp();}'| gawk -F '\t' '/^#/{print}' > $toto
gunzip -c $toto1.gz | gawk -F '\t' 'function pp() {printf("%s",$1);for(i=2;i<=NF;i++){if(0+bad[i]<1)printf("\t%s",$i);}printf("\n");}/^#Title/{ntitle++;if(ntitle>1)next;bad[2]=1;bad[8]=1;bad[9]=1;for(i=10;i<=NF;i++){if(index($i,"b")>0)bad[i]=1;gsub(/a/,"",$i);$i="SEQC_" $i;}pp();}/^#/{next;}{if(NF>10)pp();}' |  gawk -F '\t' '/^#/{next;}{if(! $3)next;}{split($1,aa,"_");printf("%2s:%09d:%09d\tchr%s:%s-%s",aa[1],aa[3],aa[4],aa[1],aa[3],aa[4]);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' |  gawk -F '\t' '{printf("%s",$2);for(i=3;i<=NF;i++)printf("\t%s",$i);printf("\n");}' >> $toto

head -100 $toto > ../RESULTS/introns.test3.txt
head -100000 $toto | tail -100 >> ../RESULTS/introns.test3.txt

gunzip -c $toto1.gz | head -100 | gawk '/^#Sample/{next;}/^#RunId/{nrid++;if(nrid==1)print;next;}/^#Run/{nr++;if(nr==1)print;next;}/^#/{print;}' | sed -e 's/Gene/Junction/g'  > $toto
gunzip -c $toto1.gz |  gawk -F '\t' '/^#/{next;}{printf("%02d_%s__%09d_%09d\t",0+$3,$3,$5,$6);print}' | sort | gawk -F '\t' '{printf ("%s",$2);for(i=3;i<=NF;i++)printf("\t%s",$i);printf("\n");}' >> $toto

cat $toto | head -40 | cut -f 1,2,3,4,5,6,7
head -30 $toto > RESULTS/$MAGIC.introns.test3.txt
head -30000 $toto | tail -1000  >> RESULTS/$MAGIC.introns.test3.txt

########################################
### General analysis

TGX2:

echo -n "# " > $toto
date >> $toto
# cat scripts/ReadMe_for_NB_AceView_nuIndex_NAhint_August28_2012.txt >> $toto
# cat scripts/ReadMe_for_NB_AceView_Index_unique_paired_Dec9_2013.txt >> $toto
# cat scripts/ReadMe_for_TGx_AceView_nuIndex_NAhint_August30_2012.txt >> $toto
# cat scripts/ReadMe_for_SEQC_main_nuIndex_NAhint_nov23_2012.txt >> $toto
# cat scripts/ReadMe_for_NB_RefSeq_nuIndex_NAhint_August28_2012.txt >> $toto
 
echo '#' | gzip >   toto.complexlocus.txt.gz
if ($GM == GENE && $target == av) then
  echo "phaseTGx.$target.$GM :: complexlocus"
  gunzip -c $toto1.gz | gawk -F '\t' '/^#/{next}{g1[$1]=1;n=split($3,aa,";");if(n>1){for(i=1;i<=n;i++)g2[aa[i]]=$1;}}END{for(k in g2)if(g1[k]<1)printf("%s\t%s\n", g2[k],k);}' | sort -u | gzip >  toto.complexlocus.txt.gz


# we verify at least in human that we can as well throw away all the AandB they are sufficiently represented in the RefSeq part
  gunzip -c $toto1.gz | gawk -F '\t' '/^#/{next}{isAnd=index($1,"and");if(isAnd<2)next;g1[$1]=1;n=split($3,aa,";");if(n>1){for(i=1;i<=n;i++)g2[aa[i]]=$1;}}END{for(k in g2)if(g1[k]<1)printf("%s\t%s\n", g2[k],k);}' | sort -u  >  toto.complexAandBlocus.txt

  gunzip -c $toto1.gz | gawk -F '\t' '/^#/{next}{g1[$1]=1;n=split($3,aa,";");if(n>0){for(i=1;i<=n;i++)g2[aa[i]]=$1;}}END{for(k in g2)if(g1[k]<1)printf("%s\n", k);}'  > toto.knownGeneId
  echo 302729  >> toto.knownGeneId
  cat toto.knownGeneId | sort -u | gzip > _toto.gz
  mv _toto.gz toto.knownGeneId.gz
endif

echo "phaseTGx.$target.$GM :: header"
if ($MAGIC == Liver || $MAGIC == Liver_export) then
  gunzip -c $toto1.gz | gawk -F '\t' '/^#RunId/{gsub(/Micro array/,"Microarray",$0);printf("%s",$1); for(i=2;i<=11;i++)printf("\t%s",$i);for(i=12;i<=NF;i++)printf("\tCOH_WANG_%s",$i);printf("\n");exit(0);}' > $toto.3
  gunzip -c $toto1.gz | gawk -F '\t' '/^#Also/{print;exit(0);}' >> $toto.3
  gunzip -c $toto1.gz | gawk -F '\t' '/^#Low /{print;exit(0);}' >> $toto.3
  gunzip -c $toto1.gz | gawk -F '\t' '/^#RunId/{gsub(/Micro array/,"Microarray",$0);printf("%s",$1); for(i=2;i<=11;i++)printf("\t%s",$i);for(i=12;i<=NF;i++)printf("\tCOH_WANG_%s",$i);printf("\n");exit(0);}' >> $toto.3

else if ($MAGIC == NB || $MAGIC == NB_export) then
  gunzip -c $toto1.gz | head -30 | gawk -F '\t' '/^#Sample/{next;}/^#Title/{ntitle++;gsub(/, misannotated_or_technical problem_ignore/,"",$0);if(ntitle==1)print;next;}/^#/{nn[$1]++;if(nn[$1]==1)print;next;}/^#Run/{nr++;if(nr==1)print;next;}/^#/{print;}' > $toto.3

else
  # AnyRun is in column 11 (case gene) and column 8 (case MRNA)
  gunzip -c $toto1.gz | head -20 | gawk -F '\t' '/^#Sample/{ns++;if(ns==1)print;next;}/^#Low/{nl++;if(nl==1)print;next;}/^#Title/{nt++;if(nt==1)print;next;}/^#RunId/{nrid++;if(nrid==1)print;next;}/^#Run/{nr++;if(nr==1)print}' | sed -e 's/RhsFake3977/L_LIV_E_1_s1_fastq_not_color/'  -e 's/detectable/measurable/' > $toto.3
endif

# locate the group columns and accept the first group
set colMax=1000000
if ($MAGIC == NB) then
  set colMax=`head -100 $toto.3 | gawk -F '\t' '/^#Run/{nr++;if(nr>1)next;for(i=80;i<=NF;i++)if(substr($i,1,3)!= "Rhs"){z=1;if(target=="AGLuK" && GM="MA")z=0;print i+z;next;}}'` target=$target GM=$GM
endif
echo "colMax=$colMax"
cat $toto.3 | gawk -F '\t' '{printf("%s",$1);mx=0;for (i=2;i<=NF && i<colMax;i++){z=$i;if(z=="#Measured object")z=$1;printf("\t%s",z);}printf("\n");}'  colMax=$colMax  >> $toto


if (-e $toto.8) then
  cat $toto.8 | gawk -F '\t' '/^#/{print ; next;}{nl++;z1=0;z2=0;n1=0;n2=0;mx2=0;for(i=200;i<1800;i++){z1+=$i;n1++;}for(i=1850;i<2650;i++){if($i>1){z2+=$i;n2++;if($i>mx2)mx2=$i;}} if (z1>10 * n1 && z1>10*z2 && (1+n2)*z1 > 5 * (1 + n1)*z2 && mx2 < 5)print }' > RESULTS/MostDiff.txt
endif

# export av genes with geneid  and the RefSeq to complement the aceview genes with multiple geneid 
if (-e $toto.2) \rm  $toto.2
touch $toto.2


if (0) then
  # put the aceview gene name then the RefSeq name in paranthesis 
  gunzip -c toto.complexlocus.txt.gz ZZZZZ.gz RESULTS/Expression/unique/RefSeq/$MAGIC.RefSeq.GENE.u.$exportType.txt.gz | gawk -F '\t' '/^ZZZZZ/{zz=1;next;}/^#/{next;}{if(zz<1){gid2g[$2]=$1;next;}z=$3;if(gid2g[z]){printf("%s",gid2g[z]);gsub(/X__/,"",$1);if($1!=gid2g[z])printf("(%s)",$1);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}}' | sed -e 's/X__//g' -e 's/\\//g' > $toto.2
endif

if ($GM == GENE && $target == av) then
  # just the refseq name
  echo "complete av.GENE with RefSeq"
  gunzip -c toto.complexlocus.txt.gz ZZZZZ.gz RESULTS/Expression/unique/RefSeq/$MAGIC.RefSeq.GENE.u.$exportType.txt.gz | gawk -F '\t' '/^ZZZZZ/{zz=1;next;}/^#/{next;}{if(zz<1){gid2g[$2]=$1;next;}z=$3;if(gid2g[z]){gsub(/X__/,"",$1);printf("%s",$1);for(i=2;i<=3;i++)printf("\t%s",$i);if (0)printf("\t\t\t");for(i=4;i<=NF && i <= NF;i++)printf("\t%s",$i);printf("\n");}}'  | sed -e 's/X__//g' -e 's/\\//g' > $toto.2
endif

touch $toto1.mask
if (-e TARGET/Targets/$species.$target.gene_to_mask.list) cat  TARGET/Targets/$species.$target.gene_to_mask.list | gawk '{gsub(/\"/,"",$0);print $2;}' > $toto1.mask
\rm $toto1.mask.gz
if (0 && $GM == GENE) cat toto.complexAandBlocus.txt >>  $toto1.mask
gzip  $toto1.mask


#duplicate the Title line and locate critical columns
\cp $toto $toto.7

if ($GM == MA) then
  if ($GM_at == MA) then
    cat $toto | grep -v Average  | grep -v 'Mb aligned' | grep -v '#denominator'  >!  $toto.7
  endif
  if ($GM_at == MA_at) then
    cat $toto | grep -v 'Mb aligned' | grep -v '#denominator' | grep -v intergenic >!  $toto.7
  endif
endif

cat $toto | gawk '/^#Title/{nt++;if(nt==1)print}' >>  $toto.7
set iSumAll=`gunzip -c $toto.gz | head -30 | gawk -F '\t' '/^#Title/{for (i=2;i<=NF;i++)if($i=="SumOfAllReadsInProject")sumall=i;}END{print 0+sumall}'`
set iRST=0
if ($target != RefSeq && ($GM != GENE || $target != av)) then
  set iRST=`cat $toto.7 | head -20 | gawk -F '\t' '/^#/{for(i=1;n<1 && i<=NF;i++)if($i=="#RefSeq transcript Id"){n++;if(n==1)print i;last;}}'`
endif

# export gene not bad AND geneid
if ($GM == GENE) then
  # export all not bad genes 
  gunzip -c $toto1.mask.gz  ZZZZZ.gz $toto1.gz | gawk -F '\t' '/^ZZZZZ/{zz=1;next;}/^#/{next}{if(zz<1){bad[$1]=1;next;}if(bad[$1])next;if (1 || $3>1 || $6){print; next;}}' GM=$GM >> $toto.2

  # complete with RefSeq for the AandB genes or the missing RefSeq genes 
  if ($target == av) then
    gunzip -c toto.knownGeneId.gz ZZZZZ.gz RESULTS/Expression/unique/RefSeq/$MAGIC.RefSeq.GENE.u.$exportType.txt.gz | gawk -F '\t' '/^ZZZZZ/{zz=1;next;}/^#/{next;}{if (zz < 1){gid[$1]=1;next;}if(length($3)<1 || gid[$3]==1)next;printf("%s",$1);for(i=2;i<=3;i++)printf("\t%s",$i);if(0)printf("\t\t\t");for(i=4;i <= NF;i++)printf("\t%s",$i);printf("\n");}' | sed -e 's/X__//g' -e 's/\\//g' >> $toto.2
  endif


else
  # export gene not bad AND not geneid and one value not NA

  if (-e $toto.gz) \rm $toto.gz

  gunzip -c $toto1.mask.gz ZZZZZ.gz $toto1.gz | gawk -F '\t' '/^ZZZZZ/{zz=1;next;}/^#/{next}{if(zz<1){bad[$1]=1;next;}if(bad[$1])next;for (i=iSumAll+1;i<=NF && i<colMax-1;i++){if(substr($i,1,3)=="NA/"||$i+0<1)continue;print;next;}}'  colMax=$colMax  GM=$GM iSumAll=$iSumAll>> $toto.2

endif


# drop the group columns and reorder by highest gene expression index
set NN=1
if (exportType == expression_index) set NN=1000
cat $toto.7  | gawk -F '\t' '{printf("%s",$1);mx=0;for (i=2;i<=NF && i<=colMax;i++){if(i == iRST)continue;printf("\t%s",$i);}printf("\n");}' iRST=$iRST colMax=$colMax > $toto
cat $toto.2 | gawk -F '\t' 'BEGIN{sumall=12;}/^ZZZZZ/{zz=1;next;}/^#Title/{nTitleLine=1;titleLine=$0;for (i=2;i<=NF;i++)if($i=="SumOfAllReadsInProject")sumall=i;next;}/^ERCC-00/{next}/^ERCC_vector/{next}{allNA=1;mx=0;for (i=2;i<=NF && i<colMax;i++){if(i>sumall){if ($i+0 > mx)mx=$i+0;}}printf("%09d\t",999999999-int(NN*mx));print;}' NN=$NN colMax=$colMax | sort | gawk -F '\t' '{printf("%s",$2);mx=0;for (i=3;i<=NF && i<=colMax;i++){if(i == iRST+1)continue;printf("\t%s",$i);}printf("\n");}' iRST=$iRST colMax=$colMax >> $toto

########## add the SpikeIn
## add this file, behind the genes, the 92 ERCC, and the ERCC vector in alphabetic order
if ($GM == GENE) then

  gunzip -c $toto1.gz | gawk -F '\t' '/^ERCC-00/{print}/^ERCC_vector/{print}' | sed -e 's/Gene_AceView/ERCC/g' | sort |  gawk -F '\t' '{ok=0;for (i=iSumAll+1;i<=NF && i<colMax-1;i++){if(substr($i,1,3)=="NA/"||$i+0<1)continue;ok=1;}if(ok==1)print}' iSumAll=$iSumAll colMax=$colMax | gawk -F '\t' '{printf("%s",$1);for (i=2;i<=NF && i<colMax;i++){if(i == iRST)continue;printf("\t%s",$i);}printf("\n");}' iRST=$iRST colMax=$colMax >> $toto
endif

########

if ($GM == GENE) then
  #eliminate column 2 and the X__
  mv $toto $toto.8 
  if ($MAGIC == Liver) then
    cat $toto.8 | gawk -F '\t' '{printf("%s",$1); gsub(/Gene_RefSeq/,"Gene_RefSeq08",$11);for(i=2;i<=10;i++)printf("\t%s",$i);for(i=11;i<=NF;i++)printf("\t%s",$i);printf("\n");}' > $toto
  else
    cat $toto.8 | sed -e 's/Gene_AceView08/Gene_AceView2010/g' -e 's/X__//g' > $toto
  endif
endif

if ($GM == MA) then
  if ($GM_at == MA) then
    # eliminate iSumAll column
    \mv $toto $toto.3  
    cat $toto.3 | gawk -F '\t' '{printf("%s",$1);for(i=2;i<=NF;i++){if(i==iSumAll)continue;printf("\t%s",$i);}printf("\n");}'  iSumAll=$iSumAll > $toto
  endif
endif
\rm  $toto.[1-9]

if ($toto == RESULTS/Expression/unique/introns/SEQC_NB_MAV_G_log2.20121127.txt) then

  mv $toto $toto.3
  cat $toto.3 | head -70 | gawk -F '\t' 'function pp() {printf("%s",$1);for(i=2;i<=NF;i++)if(bad[i]<1)printf("\t%s",$i);printf("\n");}/^#Title/{ntitle++;if(ntitle>1)next;bad[2]=1;bad[11]=1;for(i=12;i<=NF;i++){if(index($i,"b")>0 || substr($i,1,2)!="NB")bad[i]=1;gsub(/a/,"",$i);gsub(/Title/,"Gene",$1);$i="SEQC_" $i;}pp();}/^#/{next;}{if(NF>10)pp();}'| gawk -F '\t' '/^#/{print}' > $toto
  cat $toto.3 | gawk -F '\t' 'function pp() {printf("%s",$1);for(i=2;i<=NF;i++){if(0+bad[i]<1)printf("\t%s",$i);}printf("\n");}/^#Title/{ntitle++;if(ntitle>1)next;bad[2]=1;bad[11]=1;for(i=12;i<=NF;i++){if(index($i,"b")>0)bad[i]=1;gsub(/a/,"",$i);$i="SEQC_" $i;}pp();}/^#/{next;}{if(NF>10)pp();}' |  gawk -F '\t' '/^#/{next;}{printf("%s",$1);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' | sed -e 's/NA\///g' >> $toto

endif

if (0) then
  # 2013_01_12, we exported using Liver to reuse the files froam august, so we need to clip the groups
  mv $toto $toto.7
  head -1 $toto.7 > $toto
  tail -n +3 $toto.7 | gawk -F '\t' '{printf("%s",$1);for(i=2;i<=135;i++)printf("\t%s",$i);printf("\n");}' >> $toto
endif

echo "construct test3 : RESULTS/$MAGIC.$exportType.$target.$GM.test3.txt"
if (-e $toto.gz) \rm $toto.gz
gzip $toto
 gunzip -c $toto.gz | head -300 > RESULTS/$MAGIC.$exportType.$target.$GM.test3.txt             
gunzip -c $toto.gz | grep -v ERCC | head -10000 | tail -100 >> RESULTS/$MAGIC.$exportType.$target.$GM.test3.txt   
gunzip -c $toto.gz | grep -v ERCC |  tail -100 >> RESULTS/$MAGIC.$exportType.$target.$GM.test3.txt   
gunzip -c $toto.gz | grep 'ERCC-'  >> RESULTS/$MAGIC.$exportType.$target.$GM.test3.txt   

# nb of genes in file
echo "Data Lines in file"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^$/{next}{n++}END{print n}'
echo "AceView genes"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^$/{next}/RefSeq/{next}/AceView/{print $1;next;}' | sort -u | wc
echo "AceView genes nogid"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^$/{next}/RefSeq/{next}/AceView/{if ($3=="")print $1;next;}' | sort -u | wc
echo "AceView genes gid"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^$/{next}/RefSeq/{next}/AceView/{if ($3=="")next;print $1;next;}' | sort -u | wc
echo "AceView multi gid"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^$/{next}/RefSeq/{next}/AceView/{n=0;n=split($3,aa,";");if (n>1)ng++;}END{print ng;}' 
echo "RefSeq genes"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/RefSeq/{print $1;next;}' | sort -u | wc
echo "ERCC"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^ERCC-/{print $1;next;}' | sort -u | wc

echo "RefSeq expressed"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/RefSeq/{for (i=11;i<=NF;i++){if(0+$i>3){print $1;next;}}}' | sort -u | wc
echo "RefSeq not expressed"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/RefSeq/{ok=0;for (i=11;i<=NF;i++){if(0+$i>3)ok=1;}if (ok<1){print $1;next;}}' | sort -u | wc

echo "EBI expressed"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^ENS/{{print $1;next;}}' | sort -u | wc
echo "EBI expressed"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^ENS/{for (i=11;i<=NF;i++){if(0+$i>3){print $1;next;}}}' | sort -u | wc
echo "EBI not expressed"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^ENS/{ok=0;for (i=11;i<=NF;i++){if(0+$i>3)ok=1;}if (ok<1){print $1;next;}}' | sort -u | wc

echo "AceView expressed"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/AceView/{for (i=11;i<=NF;i++){if(0+$i>3){print $1;next;}}}' | sort -u | wc

echo "AceView expressed in at least 5 runs not NA"
gunzip -c $toto.gz | gawk -F '\t' '/^#Run/{nnr++;if(nnr==1){for(i=14;i<NF;i++)if(substr($i,1,3)!="Rhs"){imax=i;next}imax=NF;}}/^#/{next}/RefSeq/{next;}/AceView/{g=$1;n=0;for (i=14;i<=NF && i<imax;i++){if(substr($i,1,3)=="NA/")continue;if($i>5)n++;}nn[n]++}END{for(n in nn)printf("%d\t%d\n",n,nn[n]);}' | sort -k 1n > titi
echo -n '# ' > $toto.shared_in_population.histo.txt
date >> $toto.shared_in_population.histo.txt
echo "Number of runs\tNumber of genes not NA/\tCumul" >>  $toto.shared_in_population.histo.txt
cat titi | sort -k 1nr | gawk '{n+=$2;printf("%d\t%d\t%d\n",$1,$2,n);}' | sort -k 1n  >>  $toto.shared_in_population.histo.txt

if (0) then

endif



echo "AceView not expressed"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/AceView/{ok=0;for (i=11;i<=NF;i++){if(0+$i>3)ok=1;}if (ok<1){print $1;next;}}' | sort -u | wc
echo "AceView not expressed no gid"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/AceView/{ok=0;if($3!="")next;for (i=11;i<=NF;i++){if(0+$i>3)ok=1;}if (ok<1){print $1;next;}}' | sort -u | wc
echo "AceView expressed"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/AceView/{k=0;for (i=11;i<=NF;i++){if(0+$i>3){k++;n1++;}else n2++;}if(k>0)ng++;}END{printf("%d AceView genes contain %d values in addition to %d NA/ values\n",ng,n1,n2);}'
echo "AceView expressed"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/AceView/{k=0;if($3=="")next;for (i=11;i<=NF;i++){if(0+$i>3){k++;n1++;}else n2++;}if(k>0)ng++;}END{printf("%d AceView genes with geneid contain %d values in addition to %d NA/ values\n",ng,n1,n2);}'
echo "AceView with Affy"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^$/{next}/RefSeq/{next}/AceView/{if($6=="")next;n=0;n1++;n=split($6,aa,";");if(n>1)ng++;np+=n;}END{printf("%d genes aceview contain %d probes, %d genes contain several\n", n1,np,ng);}' 
echo "RefSeq with Affy"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^$/{next}/RefSeq/{if($6=="")next;n=0;n1++;n=split($6,aa,";");if(n>1)ng++;np+=n;}END{printf("%d RefSeq genes contain %d probes, %d genes contain several\n", n1,np,ng);}' 
echo "Gene with NM"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^$/{next}/RefSeq/{if($4=="")next;n=0;n1++;n=split($4,aa,";");if(n>1)ng++;np+=n;}END{printf("%d RefSeqd contain %d NMs, %d genes contain several\n", n1,np,ng);}' 
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^$/{next}/RefSeq/{n=0;n=split($3,aa,";");if (n>1)ng++;}END{print ng;}' 
echo "RefSeq expressed"
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/RefSeq/{k=0;for (i=11;i<=NF;i++){if(0+$i>3){k++;n1++;}else n2++;}if(k>0)ng++;}END{printf("%d RefSeq genes contain %d values in addition to %d NA/ values\n",ng,n1,n2);}'
echo 'Number of distinct geneid'
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^$/{next}{n=0;n=split($3,aa,";");for(i=1;i<=n;i++)print aa[i];}' | sort -u | wc
echo 'Number of distinct affy'
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^$/{next}{n=0;n=split($5,aa,";");for(i=1;i<=n;i++)print aa[i];}' | sort -u | wc

goto phaseLoop

# 664 affy are not seen as  expressed
cat $toto.8 | gawk -F '\t' '/^#/{next}{for (i=12;i<=NF;i++){if($i>1){print ;next;}}}' | gawk -F '\t' '/^#/{next}/^$/{next}{if (length($4)>1)print ;next;}' | gawk -F '\t' '/^#/{next}{if (length($4)<2)next;n=split($4,aa,";");for(i=1;i<=n;i++)print aa[i];}' | sort -u > tata1
 cat $toto.8 | gawk -F '\t' '/^#/{next}/^$/{next}{if (length($4)>1)print ;next;}' | gawk -F '\t' '/^#/{next}{if (length($4)<2)next;n=split($4,aa,";");for(i=1;i<=n;i++)print aa[i];}' | sort -u > tata2
diff tata1 tata2 | gawk '/^>/{i=index($2,"(");print substr($2,1,i-1);}' | sort > affy.not_expressed.list
wc  affy.not_expressed.list

cat affy.not_expressed.list ZZZZZ MicroArray/PierreBushel/SEQC_TGMX_Training_set-MAS5_normalized.PierreBushel.2012_03_01.txt | gawk -F '\t' '/ZZZZZ/{zz=1;next;}{if(zz<1){ok[substr($1,12)]=1;next;}if(ok[$1])print}' > affy.not_expressed.pierre.expression.txt

cat  affy.not_expressed.pierre.expression.txt | gawk '{for (i=2;i<= NF;i++)print int(12*$i)}' | bin/histo -plain -plot -o  affy.not
cat MicroArray/PierreBushel/SEQC_TGMX_Training_set-MAS5_normalized.PierreBushel.2012_03_01.txt | gawk '{n++;if(n<12)next;for (i=2;i<= NF;i++)print int(12*$i)}' | bin/histo -plain -plot -o  affy.any

# variable
cat $toto.8 | gawk -F '\t' '/^#/{next}{if (length($5)>1)next;print}' | gawk '/AceView/{print;next;}/RefSeq/{print}' | gawk -F '\t' '/^#/{next}{ok=0;for (i=13;i<=NF;i++){if($i>3)ok=1;}if (ok>=0){print ;next;}}' | gawk -F '\t' '/^#/{next}{if (length($3)>1)next;print}' | cut -f 1 | wc

cat SEQC_main.stable_genes.txt ZZZZZ $toto.8 | gawk -F '\t' '/ZZZZZ/{zz=1;next;}/^#/{print}{if(zz<1){gsub(/\"/,"",$1);ok[$1]=1;next;}if(ok[$1]==1)print}' > $toto.stable.txt



#######################################################
##### RefSeq
set toto1=RESULTS.$MAGIC/Expression/unique/RefSeq/Liver.RefSeq.u.expression_index.txt
set toto=RESULTS.$MAGIC/Expression/unique/RefSeq/TGx.NCBI.120330.File5.RefSeq_genes_expression_index.txt


#\cp  $toto.gz $toto.1.gz
echo -n "# " > $toto
date >> $toto
cat <<EOF >> $toto
# RNA-seq expression index of genes annotated in RefSeq genes
# The RNA-seq read-pairs were aligned to the RefSeq transcripts and the number of bases
# aligned to all the alternative RefSeq transcripts of the same genes were added up.
# The index was then computed as explained in the documentation for File3 AceView gene table
EOF


gunzip -c $toto1.gz | head -30 | gawk '/^#RunId/{n++; if(n==1)print;}' > $toto.3

set nn=`cat $toto.3 | gawk -F '\t' '/^#RunId/{for(i=10;i<=NF;i++)if($i)n=i;}END{print n}'`
cat $toto.3 | gawk -F '\t' '{printf("%s",$1);for (i=2;i<=nn;i++)printf("\t%s",$i);printf("\n");}' nn=$nn >> $toto

gunzip -c $toto1.gz | gawk -F '\t' '/^#/{next;}{printf("%s",$1);for (i=2;i<=nn;i++)printf("\t%s",$i);printf("\n");}' nn=$nn  | sed -e 's/X__//g' | sort >> $toto

\rm $toto.gz
gzip $toto


 cat $toto.8 | head -100 > RESULTS.$MAGIC/test4.txt             
cat $toto.8 | head -10000 | tail -100 >> RESULTS.$MAGIC/test4.txt   
cat $toto.8 |  tail -100 >> RESULTS.$MAGIC/test4.txt   

##### RefSeq in read-counts
set toto1=RESULTS.$MAGIC/Expression/unique/RefSeq/Liver.RefSeq.u.reads_aligned_per_gene.txt
set toto=RESULTS.$MAGIC/Expression/unique/RefSeq/TGx.NCBI.120330.File6.RefSeq_genes_reads_counts.txt


#\cp  $toto.gz $toto.1.gz
echo -n "# " > $toto
date >> $toto
cat <<EOF >> $toto
# RNA-seq read counts
# The RNA-seq read-pairs were aligned to the RefSeq transcripts, only unique mapping were kept
# The table reports the number of reads aligning in the union of all XM/NM transcripts of a given gene
EOF


gunzip -c $toto1.gz | head -30 | gawk '/^#RunId/{n++;if(n==1)print;}' > $toto.3

set nn=`cat $toto.3 | gawk -F '\t' '/^#RunId/{for(i=10;i<=NF;i++)if($i)n=i;}END{print n}'`
cat $toto.3 | gawk -F '\t' '{printf("%s",$1);for (i=2;i<=nn;i++)printf("\t%s",$i);printf("\n");}' nn=$nn >> $toto

gunzip -c $toto1.gz | gawk -F '\t' '/^#/{next;}{printf("%s",$1);for (i=2;i<=nn;i++)printf("\t%s",$i);printf("\n");}' nn=$nn  | sed -e 's/X__//g' | sort >> $toto

\rm $toto.gz
gzip $toto


gunzip -c $toto.gz | head -100 > RESULTS.$MAGIC/test6.txt             
gunzip -c $toto.gz | head -10000 | tail -100 >> RESULTS.$MAGIC/test6.txt   
gunzip -c $toto.gz |  tail -100 >> RESULTS.$MAGIC/test6.txt   

########## from previous NB export


gunzip -c $toto1.gz | gawk -F '\t' '/^#/{next}{if ($3>1){print; next;}}' >> $toto.2
gunzip -c $toto1.gz | gawk -F '\t' '/^#/{next}{if ($3>1){next;}for (i=12;i<=nn;i++)if($i>3){print;next;}}' nn=$nn >> $toto.2
gunzip -c toto.knownGeneId.gz ZZZZZ.gz RESULTS/Expression/unique/RefSeq/$MAGIC.RefSeq.u.expression_index.txt.gz | gawk -F '\t' '/^ZZZZZ/{zz=1;next;}/^#/{next;}{if (zz < 1){gid[$1]=1;next;}if(gid[$3]==1)next;print}' | sed -e 's/X__//g' -e 's/\\//g' >> $toto.2


#### comparative map

set toto=RESULTS.$MAGIC/Expression/unique/av/$MAGIC.expression_index.comparison.mieg.2012_03_19.txt
set toto=RESULTS.$MAGIC/Expression/unique/av/TGx.NCBI.120330.File4.Comparison_Affymetrix_RNAseq_per_gene.txt
echo -n "# " > $toto
date >> $toto
cat <<EOF >> $toto
# Comparison of expression index measured by RNAseq and by Affymetrix arrays in 63 rats from the FDA toxicogenomics project.
# This table is limited to genes with uniquely mapped Affymetrix probesets. A complete table of RNAseq 
# expression index of uniquely mapped Affymetrix probes is given separately. 
# Reciprocally, this table is also limited to genes containing a uniquely mapped Affymetrix probeset. 
# Tables giving the RNAseq expression index for all genes defined by AceView 2008 transcripts,  and separately 
# for all genes defined by RefSeq transcripts are provided separately. 
#
# This table allows to compare levels of expression of genes measured in different ways: by microarray 
# analysis, by RNAseq limited to the region of the probeset, or by RNAseq of the entire gene defined 
# either by its AceView 2008 or RefSeq transcripts.
# 1) Affymetrix array: The expression of each probeset was measured on an Affymetrix microarray and
#  normalized by Pierre Bushel.
# This value is denoted Affymetrix_Pierre_Bushel in the second column of the table.
# 2) RNA-Seq at probeset location: Each Affymetrix probe was mapped to the genome and the genes 
# (AceView or RefSeq transcripts), allowing at most one mismatch. 
# The best mappings were recorded and all non-uniquely mapped probes were discarded. RNA-seq coverage 
# of the genome at the location of each uniquely mapped probe per probeset was measured,
# and the average of these values was attributed to each probeset. The index was computed as 
# log base 2 of this number, divided by the sum of coverage in all uniquely mapped probes in the microarray.  
# These values are denoted Probeset_RNAseq in column 2.
# The list of well mapped probes and their RNAseq index are provided as another table.
# 3) Each read-pair was mapped to the AceView genes, pairs mapping to a unique gene (although possibly
# to several alternative transcripts of that gene) contribute to the expression index of the gene. This
# value is denoted AceView gene.
# 4) Each read-pair was also mapped to the RefSeqs, pairs mapping to a unique RefSeq gene (although possibly
# to several alternative NM transcripts of that gene) inherited the expression index of the gene. This
# value is denoted RefSeq gene.
# The table was then sorted by genes. Missing values were considered too small to be called reliably.
# ATTENTION: the different measures have not been normalized in the same way. The values are shifted
# but the important question is the covariance of all measures in each gene across all experiments.
# Some genes called A_and_B are merged in AceView and usually contain several probesets. The split
# values are given in the other tables.
#
# EXPRESSION INDEX
EOF

# recover the list is the different files
gunzip -c   tmp/GENEINDEX/Results/$MAGIC.av.u.expression_index.txt.gz  | gawk -F '\t' '/^#/{next;}{if(length($5)>1)print}' > $toto.1

# recover the affy list
cat $toto.1 | gawk -F '\t' '{n=split($5,aa,";");for(j=1;j<=n;j++){i=index(aa[j],":");printf("%s\t%s\t%s\t%s\n", $1,aa[j],substr(aa[j],i+1),$3);}}' | sort -u | gzip > toto.affy.list.gz
cat $toto.1 | gawk -F '\t' '{n=split($3,aa,";");for(j=1;j<=n;j++){printf("%s\tzorglub\tzorglub\t%s\n", $1,aa[j]);}}' | sort -u | gzip >> toto.affy.list.gz
# recover the RefSeq
gunzip -c  toto.affy.list.gz ZZZZZ.gz  tmp/GENEINDEX/Results/$MAGIC.RefSeq.u.expression_index.txt.gz  | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){aa[$2]=1;a2g[$2]=$1;a2b[$2]=$3;gid2g[$4]=$1;gg[$1]=1;next;}z=$1;gsub(/X__/,"",z);if (gg[z]==1){print;next;}if(gid2g[$3]){printf("%s(%s)",gid2g[$3],$1);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");next;}n=split($5,bb,";");for(j=1;j<=n;j++){if(aa[bb[j]]==1){printf("%s(%s)",a2g[bb[j]],$1);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");next;}}}' | sed -e 's/X__//g' >> $toto.1


# recover the affy counts in ProbeSet
gunzip -c toto.affy.list.gz ZZZZZ.gz  RESULTS.$MAGIC/MicroArray/Affy.Rat230_2.probeset.expression_index.txt.gz | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){aa[$2]=1;a2g[$2]=$1;a2b[$2]=$3;gg[$1]=1;next;};if(aa[$1]==1){printf("%s\tProbeset_RNAseq\t\t\t%s\t\t\t\t\t\t%s:%s:RNAseq",a2g[$1],$1,a2g[$1],a2b[$1]);for(i=3;i<=NF;i++)printf("\t%s",$i);printf("\n");}}' >> $toto.1
# recover the affy counts in Pierre Bushel
gunzip -c toto.affy.list.gz ZZZZZ.gz MicroArray/PierreBushel/PB.training.renamed.expression_index.txt.gz | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){aa[$2]=1;a2g[$2]=$1;a2b[$2]=$3;gg[$1]=1;next;};if(aa[$1]==1){printf("%s\tAffymetrix_Pierre_Bushel\t\t\t%s\t\t\t\t\t\t%s:%s:Affy",a2g[$1],$1,a2g[$1],a2b[$1]);for(i=3;i<=NF;i++)printf("\t%s",$i);printf("\n");}}' >> $toto.1

gunzip -c tmp/GENEINDEX/Results/$MAGIC.av.u.expression_index.txt.gz | head -20 | gawk -F '\t' '/#RunId/{nnr1++;if(nnr1==1)print;next;}/#Run/{next;nnr2++;if(nnr2==1)print;next;}/^#Title/{nnt++;if(nnt==1)print;next;}' > $toto.2
set n=`cat  $toto.2 | gawk -F '\t' '/^#RunId/{mx=0;for(i=12;i<NF;i++)if(length($i)>2)mx=i;}END{print mx;}'`


cat  $toto.2 | gawk -F '\t' '{printf("%s",$1);for (i=2;i<=n;i++)printf("\t%s",$i);printf("\n");}' n=$n >> $toto
gunzip -c tmp/GENEINDEX/Results/$MAGIC.av.u.expression_index.txt.gz | gawk '/^RunId/{next}/^#Run/{print ; exit ;}' > $toto.3
gunzip -c  MicroArray/PierreBushel/PB.training.renamed.expression_index.txt.gz | head -15 | gawk -F '\t' '/^#Sample/{n++;if(n==1)print}' | sed -e 's/Sample/Affymetrix array/' -e 's/Measured object/Method\t\t\t\t\t\t\t\t\t/' >>  $toto
gunzip -c tmp/GENEINDEX/Results/$MAGIC.av.u.expression_index.txt.gz | gawk '/^#Gene/{print ; exit ;}' >> $toto

cat $toto.1 | sort | gawk -F '\t' '{printf("%s",$1);for (i=2;i<=n;i++)printf("\t%s",$i);printf("\n");}' n=$n >> $toto
\rm $toto.[12]
\rm $toto.gz
gzip $toto
 gunzip -c $toto.gz | head -100 > RESULTS.$MAGIC/test5.txt             
gunzip -c $toto.gz | head -10000 | tail -100 >> RESULTS.$MAGIC/test5.txt   
gunzip -c $toto.gz |  tail -100 >> RESULTS.$MAGIC/test5.txt   



# nb of genes in file
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^$/{next}/RefSeq/{next}/av/{print $1;next;}' | sort -u | wc
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/^$/{next}/RefSeq/{next}/av/{if ($3=="")print $1;next;}' | sort -u | wc
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/RefSeq/{print $1;next;}' | sort -u | wc
# nb genes with some expressiom
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/RefSeq/{ok=0;for (i=12;i<=NF;i++){if($i>3)ok=1;}if (ok<1){print $1;next;}}' | sort -u | wc
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/RefSeq/{for (i=12;i<=NF;i++){if($i>3){print $1;next;}}}' | sort -u | wc
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}/av/{for (i=12;i<=NF;i++){if($i>3){print $1;next;}}}' | sort -u | wc
gunzip -c $toto.gz | gawk -F '\t' '/^#/{next}{for (i=12;i<=NF;i++){if($i>3){print $1;next;}}}' | sort -u | wc

## compute the correlation coefficient
set toto=/home/mieg/ACEVIEWHELP/RESULTS.Rat_liver_apr4/Expression/unique/av/TGx.NCBI.120330.File4.Comparison_Affymetrix_RNAseq_per_gene

## This tcsh scripts compute the correlation coefficient from the public File4
set toto=TGx.NCBI.120330.File4.Comparison_Affymetrix_RNAseq_per_gene

set min=2
set median=1

set toto2=$toto.correl.any.median$median.min$min.txt
set toto3=$toto.correl.any.median$median.min$min.txt2
echo  -n "# " > $toto2
date >> $toto2
echo "# Correlation between all measures in File4" >> $toto2
gunzip -c $toto.txt.gz | gawk -F '\t' '/^#RunId/{printf("\t\t");print;}' >> $toto2
echo  -n "# " > $toto3
set n1=0
foreach type1 (Affymetrix_Pierre_Bushel Gene_AceView2010 Gene_RefSeq Probeset_RNAseq)
  @ n1 = $n1 + 1
  set n2=0
  foreach type2 (Affymetrix_Pierre_Bushel Gene_AceView2010 Gene_RefSeq Probeset_RNAseq)
    @ n2 = $n2 + 1
    if ($n1 >= $n2) continue

    set toto1=$toto.correl.$type1.$type2.median$median.min$min.txt
    echo  -n "# " > $toto1
    date >> $toto1
    echo "# Correlation between $type1 and $type2 measures in File4" >> $toto1
    echo "# Remove index below $min, renormalize each gene to 10 yes/no=$median" >> $toto1
    gunzip -c $toto.txt.gz | gawk -F '\t' '/^#RunId/{printf("\t\tRunId\tAll");for(i=12;i<=NF;i++){rid[i]=$i;printf("\t%s",$i);}nf=NF;next;}/^#/{next;}{split($1,aa,"(");$1=aa[1];}{n=0;m=0;for(i=12;i<=NF;i++){if($i+0>0){m+=0+$i;n++}}if(n==0)next;m=m/n;if(median==1){if(m<0+qxMin)next;for(i=12;i<=NF;i++)if($i>0)$i=$i-m+10;}}/'$type1'/{p=$1;for(i=12;i<=NF;i++)px[p,i]=0+$i;pp[p]=1;next;}/'$type2'/{q=$1;for(i=12;i<=NF;i++)qx[q,i]=0+$i;qq[q]=1;next;}END{np=0;nq=0;npq=0;for(i=12;i<=nf;i++)npqi[i]=0;for(p in pp)np++;for(p in qq){nq++;if(pp[p]){for(i=12;i<=nf;i++){if(( px[p,i]>=qxMin && qx[p,i]>=qxMin)){npq++;npqi[i]++;ZX+=px[p,i];X[i]+=px[p,i];ZX2+=px[p,i]*px[p,i];X2[i]+=px[p,i]*px[p,i];ZY+=qx[p,i];Y[i]+=qx[p,i];ZY2+=qx[p,i]*qx[p,i];Y2[i]+=qx[p,i]*qx[p,i];ZXY+=px[p,i]*qx[p,i];XY[i]+=px[p,i]*qx[p,i];}}}}printf("\nN\t%d",npq);for(i=12;i<=nf;i++){printf("\t%d",npqi[i]);if(npq==0)npq=1;if(npqi[i]==0)npqi[i]=1}printf("\nX\t%.3f",ZX);for(i=12;i<=nf;i++)printf("\t%.3f",X[i]);printf("\nY\t%.3f",ZY);for(i=12;i<=nf;i++)printf("\t%.3f",Y[i]);printf("\nX2\t%.3f",ZX2);for(i=12;i<=nf;i++)printf("\t%.3f",X2[i]);printf("\nY2\t%.3f",ZY2);for(i=12;i<=nf;i++)printf("\t%.3f",Y2[i]);printf("\nXY\t%.3f",ZXY);for(i=12;i<=nf;i++)printf("\t%.3f",XY[i]);x=ZX2/npq-ZX*ZX/(npq*npq);y=ZY2/npq-ZY*ZY/(npq*npq);w=ZXY/npq-ZX*ZY/(npq*npq);if(x*y>0)w=w/sqrt(x*y);printf("\nCorrel\t%.3f",w);for(i=12;i<=nf;i++){x=X2[i]/npqi[i]-X[i]*X[i]/(npqi[i]*npqi[i]);y=Y2[i]/npqi[i]-Y[i]*Y[i]/(npqi[i]*npqi[i]);w=XY[i]/npqi[i]-X[i]*Y[i]/(npqi[i]*npqi[i]);w=w/sqrt(.1+x*y);printf("\t%.3f",w);}printf("\n");}' qxMin=$min median=$median  >> $toto1

    echo -n "$type1\t$type2\t" >>  $toto2
    cat $toto1 | gawk '/^N/{print}'  | gawk '{printf("%s\t%s\t",t1,t2);print;}' t1=$type1 t2=$type2 >> $toto3
    tail -1 $toto1 >> $toto2
  end
end

cat $toto3  >> $toto2
\rm $toto3

#############################
##  84 genes tith 80 fold change

cat <<EOF > 80fold.gene.list
Abcb1b
Acot1
ADAM_CR.2
Adam8
Agpat9
Akr1b10
Akr1b10
Akr1b7
Akr1b8
Apoa4
Aqp3
Aqp7
Atf3
beyloy
Car3
Ccl2
Chrna2
chugar
Cpt1b
Cyp1a1
Cyp1b1
Cyp2b1
Cyp2c
Cyp2c24
Cyp2c7
Defb1
Dhrs7
Fabp3
Fam25a
Fosl1
gawdar
Gdf15
glawsa
Glul
Gpx2
Gstp1
Hba-a2
Hbb
Hdc
Hspa1b
Inmt
korflo
Krt23
Lcn2
Lilrb4
LOC259246
LOC287167
LOC292543
LOC298109
LOC298109
LOC360504
LOC688514
LOC690226
LOC690768
Maff
Marco
MGC72973
Mlc1
nyto
Obp3
Olr1
Otop1
p450.1
p450.3
Pla2g7
pleejor
RGD1305928
RGD1561849
RGD1562284
RGD1562373
RGD1562844
Rgs13
SAA.1
Serpine1
Sfn
shywey
skarchy
smuflu
speejar
speypey
Stac3
sterbor
Ucp3
Ugt2b
EOF

gzip  80fold.gene.list
## extract the restricted comparison table
set toto=RESULTS.$MAGIC/Expression/unique/av/$MAGIC.expression_index.comparison.mieg.2012_03_19.txt
set toto1=RESULTS.$MAGIC/Expression/unique/av/$MAGIC.expression_index.comparison.84_differential_genes.mieg.2012_03_19.txt
echo -n "# " > $toto1
date >> $toto1
gunzip -c $toto.gz | head -100 | gawk '/^#/{print}' >> $toto1
gunzip -c  80fold.gene.list.gz ZZZZZ.gz $toto.gz | gawk -F '\t' '/^#/{next}/^ZZZZZ/{zz=1;}{if(zz<1){gg[$1]=1;next;}}{split($1,bb,"(");if(gg[bb[1]]==1){print;next;}}'   >> $toto1


################
# extract the SU list of 16 genes

# construct the list
cat <<EOF > su.gene.list
Adamtsl4  idem pierr, plat sauf dip pour les 2 lef delta av-pierr = 6
Ahrr      tres faible presque absent en RNA seq, variable mais plat chez pierre delta av-pierr = 6
Aldh3a1   OUI sureexpression av/refseq/pierre delta av-pierr = 6
Cdk2      ***super plat av/refseq/peirre delta av-pierr = 6
Ceacam10  surexprime un peu dans LEF, variable dna les controles  
Ceacam20  surexprime un peu dans LEF, pas mesure chez pierre
Cyp1a1    OUI  plus fluctuant chez pierre
Cyp1a2    OUI  pas surexprime chez pierre
Dcun1d1   ***super plat, chez pierre et av
Deaf1     ***super plat, chez pierre et av   ## les groupes chutent
Dmgdh     ***plat, chez pierre et av legere surexpression  2X mais aussi dans les controles
Il11      plat et sous le radar
Mark3     ***super plat, chez pierre et av   ## les groupes chutent
Stmn1     un peu sous exprime, mais seulement dans LEF et tres fluctuant dans le debut
Tmem86b   2x seuelemnt dans av, fluctue dans les controles
Vav2      ***super plat, chez pierre et av
EOF

# good list
cat <<EOF > su.gene.list
Adam8
Atf3
SAA.1
Lvn2
Hspa1b
Scn4b
Pnpla3
Pnpla5
Acot1 
EOF

if (-e  su.gene.list.gz) \rm  su.gene.list.gz
gzip su.gene.list

# recover the list is the different files
gunzip -c  su.gene.list.gz ZZZZZ.gz tmp/GENEINDEX/Results/$MAGIC.av.u.expression_index.txt.gz  | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){split($1,aa," ");gg[aa[1]]=1;next;}if(gg[$1]){printf("%s\tGene_AceView",$1);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}}' > tmp/GENEINDEX/g4.su.txt
gunzip -c  su.gene.list.gz ZZZZZ.gz tmp/GENEINDEX/Results/$MAGIC.RefSeq.u.expression_index.txt.gz  | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){split($1,aa," ");gg[aa[1]]=1;next;}z=$1;gsub(/X__/,"",z);gsub(/_predicted/,"",z);if(gg[z]){printf("%s\tGene_RefSeq",z);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}}' | sed -e 's/X__//g' >> tmp/GENEINDEX/g4.su.txt
# recover the affy list
cat tmp/GENEINDEX/g4.su.txt | gawk -F '\t' '/Gene_AceView/{n=split($5,aa,";");for(j=1;j<=n;j++){i=index(aa[j],":");printf("%s\t%s\n", $1,substr(aa[j],i+1));}}' | sort -u | gzip > su.affy.list.gz
# recover the affy counts in ProbeSet
gunzip -c su.affy.list.gz ZZZZZ.gz  RESULTS.$MAGIC/MicroArray/Affy.Rat230_2.probeset.expression_index.txt.gz | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){gene[$2]=$1;gg[$2]=1;next;};split($1,aa,":");if(gg[aa[2]]){printf("%s\tProbeset_RNAseq\t\t\t%s\t\t\t\t\t\t%s",gene[aa[2]],$1,gene[aa[2]]);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}}' >> tmp/GENEINDEX/g4.su.txt
# recover the affy counts in Pierre Bushel
gunzip -c su.affy.list.gz ZZZZZ.gz MicroArray/PierreBushel/PB.training.renamed.expression_index.txt.gz | gawk -F '\t' '/^#/{next;}/^ZZZZZ/{zz++;next;}{if(zz<1){gene[$2]=$1;gg[$2]=1;next;};split($1,aa,":");if(gg[aa[2]]){printf("%s\tPierre_Affy\t\t\t%s\t\t\t\t\t\t%s",gene[aa[2]],$1,gene[aa[2]]);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}}' >> tmp/GENEINDEX/g4.su.txt
mv  tmp/GENEINDEX/g4.su.txt tutu.$$

gunzip -c tmp/GENEINDEX/Results/$MAGIC.av.u.expression_index.txt.gz | head -4 | gawk -F '\t' '/#RunId/{printf("#EXPRESSION INDEX\t#Method");for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");next;}/#Run/{printf("Magic_id\t#Method");for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' > toto.$$
cat  toto.$$ | gawk '/^#EXPRESSION INDEX/{print}' > tmp/GENEINDEX/g4.su.txt
set n=`cat  tmp/GENEINDEX/g4.su.txt | gawk -F '\t' '/^#EXP/{mx=0;for(i=12;i<NF;i++)if(length($i)>2)mx=i;print mx;}'`
cat MicroArray/PierreBushel/run2cel.txt toto.$$ | gawk -F '\t' '/^Magic_id/{NF=n;printf("#Affymetrix id");for(i=3;i<=NF;i++)printf("\t%s",affy[$i]);printf("\n");next;}{gsub(/PB/,"",$1);affy[$1]=$2;}' n=$n >> tmp/GENEINDEX/g4.su.txt

cat tutu.$$ | sort | gawk -F '\t' '{printf("%s",$1);for (i=2;i<=n;i++)printf("\t%s",$i);printf("\n");}' n=$n >> tmp/GENEINDEX/g4.su.txt
\rm tutu.$$
\cp  tmp/GENEINDEX/g4.su.txt RESULTS.$MAGIC/g4.su2.txt


if (-e GeneIndexDB/parse.index.done)  goto phaseLoop


  touch  GeneIndexDB/parse.index.done

goto phaseLoop

# AUC parse FDR ounts Diff_
#######################################################################################
#######################################################################################
# centralization des elements mrna/gene + gene origin
phaseg5:

if (-e  GeneIndexDB/parse.genes.done) goto phaseLoop
  date

# destroy previous annot
if (0) then
  echo -n 'destroy old intron support counts in GeneIndexDB'
  $bin/tacembly GeneIndexDB << EOF 
    query find element Gene && Deep
    edit -D Deep
    save
    quit
EOF
endif

# parse new one

  echo >! _r
  echo 'read-models' >> _r

foreach target ($Etargets)
  foreach uu (u nu)
      if (-e RESULTS/Expression/AceFiles/$MAGIC.$target.$uu.ace.gz) then
        echo "pparse RESULTS/Expression/AceFiles//$MAGIC.$target.$uu.ace.gz" >> _r
      endif
  end
end
  echo save >> _r
  echo quit >> _r

  echo -n 'parse the new gene support counts in GeneIndexDB'
  bin/tacembly GeneIndexDB < _r

  touch  GeneIndexDB/parse.genes.done


goto phaseLoop

#######################################################################################
#######################################################################################
# centralization des elements mrna/gene + gene origin
phaseii5:

# goto phaseLoop

if (-e  GeneIndexDB/parse.introns.done) goto phaseLoop
  date

# destroy previous annot
if (0) then
  echo -n 'destroy old intron support counts in GeneIndexDB'
  $bin/tacembly GeneIndexDB << EOF 
    query find element Intron && Deep
    edit -D Deep
    save
    quit
EOF
endif

# parse new one

  echo >! _r
  echo 'read-models' >> _r

  foreach uu (u nu)
    foreach target (intron)
      if (-e RESULTS/Expression/AceFiles/$MAGIC.introns.$uu.ace.gz) then
        echo "pparse  RESULTS/Expression/AceFiles/$MAGIC.introns.$uu.ace.gz" >> _r
      endif
    end
  end
  echo save >> _r
  echo quit >> _r

  echo -n 'parse the new intron support counts in GeneIndexDB'
  bin/tacembly GeneIndexDB < _r

  touch  GeneIndexDB/parse.introns.done

goto phaseLoop

#######################################################################################
#######################################################################################
#  Exportation a partir de GeneIndexDB des tables d'index par genes annotes
phaseii6:

# goto phaseLoop

if (-e  GeneIndexDB/parse.introns.done) goto phaseLoop
  date

if (! -e  tmp/introns/ii6.$MAGIC.intronsDeUno2runs.txt) then
  bin/tacembly GeneIndexDB <<EOF
    query find intron de_uno
    select -o tmp/introns/ii6.$MAGIC.intronsDeUno2runs.txt select ii,r,n from ii in @,r in ii->de_uno,n in r[1]
EOF
endif
New_minS
bin/tacembly MetaDB <<EOF
  query find project IS $MAGIC ; >run ; intron
  select -o MetaDB/$MAGIC/ii6.run_deuno_info.txt r,t,kb from r in @,a in r->ali, nh in a->nh_ali where nh == "any", t in nh[3], kb in nh[5]
  quit
EOF

set toto=RESULTS/Introns_exons_polyA/$MAGIC.intronsDeUno2runs.txt 
echo -n "### $toto " > $toto
date >> $toto
echo "## Read counts per intron in project $MAGIC" >> $toto
cat  MetaDB/$MAGIC/ii6.run_deuno_info.txt ZZZZZ tmp/introns/ii6.Transloc.intronsDeUno2runs.txt ZZZZZ| head -5000 | gawk -F '\t' '/^ZZZZZ/{if(zz<1){zz=1;printf("##\t\t\t\t\t");for(r=1;r<=rMax;r++)printf("\t%s",i2r[r]) ;printf("##\t\t\t\t\tReads aligned\t%d",allT);for(r=1;r<=rMax;r++)printf("\t%d",i2t[r]) ;printf("\n## Intron\t\t\t\t\tkb aligned\t%d",allK);for(r=1;r<=rMax;r++)printf("\t%d",i2k[r]) ;next;}}{if(zz<1){i++;r=$1;i2r[i]=r;r2i[i]=r;i2t[i]=$2;i2k[i]=$3;allR+=$2;allK+=$3;rMax=i;if(0)print ("YYYY %d %s\n",i,r);next;}}{ii=$1;if(ii!=old){old=ii;if(nn>=0){printf("\n%s\t%s\t%s\t%d",ii,gene,mrna,nn);nn=0;for(r=1;r<=rMax;r++){printf("\t%d",ii2r[r]);ii2r[r]=0;}}ii=$1;r=$2;i=r2i[r]+0;}i=r2i[$2];if(i<1)next;ii2r[i]=$3;nn+=$3;}END{printf("\n");}'  | head




goto phaseLoop

#######################################################################################
#######################################################################################
# special export to feed the public aceview server
phaseii99:
goto phaseLoop

# pick in the different project and feed a dedicated tag
# to add a tag, please edit gall/acedata/wspec/models.wrm
# and grep for nTagsP in wabi/hseqdisp.c

cd ~/SEQC_2013/TOTOA
cat  tmp/INTRON_INDEX/Seqc_A.introns.u.ace | gawk '/^Intron/{ii=$0;next;}/^Run_U/{n[ii]+=$6;next;}END{for(k in n)if(n[k]>0)printf("%s\tUHR\t%d\n",k,n[k]);}' > /tmp/toto345.UHR
cd ~/SEQC_2013/TOTOC
cat  tmp/INTRON_INDEX/Seqc_C.introns.u.ace | gawk '/^Intron/{ii=$0;next;}/^Run_U/{n[ii]+=$6;next;}END{for(k in n)if(n[k]>0)printf("%s\tC\t%d\n",k,n[k]);}' > /tmp/toto345.C
cd ~/SEQC_2013/TOTOD
cat  tmp/INTRON_INDEX/Seqc_D.introns.u.ace | gawk '/^Intron/{ii=$0;next;}/^Run_U/{n[ii]+=$6;next;}END{for(k in n)if(n[k]>0)printf("%s\tD\t%d\n",k,n[k]);}' > /tmp/toto345.D
cd ~/SEQC_2013/TOTOA
cat  tmp/INTRON_INDEX/Seqc_B.introns.u.ace | gawk '/^Intron/{ii=$0;next;}/^Run_U/{n[ii]+=$6;next;}END{for(k in n)if(n[k]>0)printf("%s\tBrain\t%d\n",k,n[k]);}' > /tmp/toto345.Brain


cd ~/Fatigue
cat tmp/INTRON_INDEX/Fatigue.introns.u.ace | gawk '/^Intron/{ii=$0;next;}/^Run_U/{n[ii]+=$6;next;}END{for(k in n)if(n[k]>0)printf("%s\tOther_sample\t%d\n",k,n[k]);}' > /tmp/toto345.fatigue
cd ~/NB
cat tmp/INTRON_INDEX/NB.introns.u.ace | gawk '/^Intron/{ii=$0;next;}/^Run_U/{n[ii]+=$6;next;}END{for(k in n)if(n[k]>0)printf("%s\tNB\t%d\n",k,n[k]);}' > /tmp/toto345.NB
cd ~/NB
cat tmp/INTRON_INDEX/Pain.introns.u.ace | gawk '/^Intron/{ii=$0;next;}/^Run_U/{n[ii]+=$6;next;}END{for(k in n)if(n[k]>0)printf("%s\tOther_sample\t%d\n",k,n[k]);}' > /tmp/toto345.iada
cd ~/Primates_Tissue
cat  tmp/INTRON_INDEX/Tissue.introns.u.ace | gawk '/^Intron/{ii=$0;next;}/^Run_U/{n[ii]+=$6;next;}END{for(k in n)if(n[k]>0)printf("%s\tPrimates\t%d\n",k,n[k]);}' > /tmp/toto345.Primates

cat ~/ace//toto345.* | gawk -F '\t' '{if($3>0)print;}' | sort > ~/ace/toto346

cat ~/ace/toto346 | gawk '{i=$2;if(i!=old)printf("\nIntron %s\n-D Other_sample\n",$2);old=$2;printf("%s %d\n",$3,$4);}' > ~/ace/toto346.ace





goto phaseLoop

#########################################################################################
#########################################################################################

$tacembly GeneIndexDB <<EOF
  table -o GeneUnicityIndex.txt -f $metaData/tables.GeneIndexDB/GeneUnicityIndex.def
EOF
cat  GeneIndexDB/GeneUnicityIndex.txt | gawk -F '\t' '/\"/{n = $5 ; if(n == 0)n = 1 ;printf("Gene %s\nUnicity_index %.1f\n\n", $1, 100*$3/n);}' > TARGET/GENES//GeneUnicityIndex.ace

foreach target (av)
  cat  GeneIndexDB/gene2length.$target.ace > GeneIndexDB/toto.$target.g2d
  cat  GeneIndexDB/gene2deep.$target.ace | gawk '/^$/{print}/^Gene/{print}/^Group_nU/{next}/sample_/{print}' >> GeneIndexDB/toto.$target.g2d
  date ; bin/bestali -gene2deep2index -db GeneIndexDB -i toto.$target.g2d -target_class $target > toto.g2d8.$target.txt ; date
  date ; bin/bestali -gene2deep2variance -db GeneIndexDB -i toto.$target.g2d -target_class $target > toto.g2variance.$target.txt ; date
  cat toto.g2d8.$target.txt | bin/histo -columns 2+ -o toto16.8.$target -plain -plot -min 5 -title $target.3kb_limit -w 50
  cat toto.g2d8.$target.txt | bin/histo -columns 2,3,4,5 -o toto16.8.ABCD -plain -plot -min 5 -title $target.3kb_limit -w 50
  cp toto16.8.$target.txt RESULTS/Expression/$MAGIC.geneindex.$target.3kb_limit.histo.txt
  cp toto16.8.ACDB.txt RESULTS/Expression/$MAGIC.geneindex.ACDB.histo.txt

  cat toto.g2d8.$target.txt | gawk -F '\t' 'BEGIN{if(0){dd[2]=.01;dd[3]= +.0;dd[4]=.01;}}{line++;if(line==1){ok[1]=1;for(j=1;j<=4;j++)z[j]=-1;for(i=2;i<=NF;i++)if(index($i,"sample_")==1){j++;ok[i]=j;}}for(j=1;j<=4;j++)z[j]=0;printf("%s",$1);j=0;for(i=2;i<=NF;i++)if($1 && ok[i]){j++;z[j]=$i+dd[j];if(0+z[j]<0){z[j]=0;$i="NA";}printf("\t%s",$i+dd[j]);}if(z[1]<z[2] && z[2]<z[3]&&z[3]<z[4] && z[3]>0 && z[4]>=7){if(z[3]>0)dz=z[3];if(z[2]>0)dz=z[2];if(z[1]>0)dz=z[1];dz=z[4]-dz;if(dz>delta)printf("\tTitrating_woman\t%.2f",dz);}if(z[1]>=7 && z[1]>z[2] && z[2]>0 && z[2]>z[3]&& z[3]>z[4]){if(z[2]>0)dz=z[2];if(z[3]>0)dz=z[3];if(z[4]>0)dz=z[4];dz=z[1]-dz;if(dz>delta)printf("\tTitrating_man\t%.2f",dz);}printf("\n");}' delta=.3 | gawk '{n++;if(n==1){print;next;}}/Titrating/{print}' > _t1

  head -1 _t1 > titrating_genes.man.txt 
  cat _t1 | tail -n +2 | sort -k 7nr | grep Titrating_man >>  titrating_genes.man.txt 
  head -1 _t1 > titrating_genes.woman.txt 
  cat _t1 | tail -n +2 | sort -k 7nr | grep Titrating_woman >>  titrating_genes.woman.txt 
wc titrating_genes.*.txt
\cp titrating_genes.*.txt RESULTS.Stan

# method by comparing the variances
cat toto.g2d8.$target.txt ZZZZZ toto.g2variance.$target.txt | gawk -F '\t' '/ZZZZZ/{zz=1;next;}{if(zz<1){g=$1;z[g]=$0;a=$2;b=$5;if(a<=7)a=7;if(b<=7)b=7;d[g]=b-a;next;}g=$1;a=$2;b=$5;ab=0;if(a>0 && b>0)ab=a+b;if(a<=0 && b>0)ab=2*b;ab=2*ab;if(a<=0 && b>0)ab=2*b;d2=d[g];if(d2<0)d2=-d2;if(d2>delta && ab>0 && d2>ab){printf("%s\t%s\t%s\t%.2f",z[g],$2,$5,d[g]);if(d[g]<0)printf("\tTitrating_man\n");else printf("\tTitrating_woman\n");}}' delta=.3  > _t3
  echo "#Titrating by variance" > titrating_genes_by_variance.man.txt 
  cat _t3 | tail -n +2 | sort -k 8nr | grep Titrating_man >>  titrating_genes_by_variance.man.txt 
  echo "#Titrating by variance" > titrating_genes_by_variance.woman.txt 
  cat _t3 | tail -n +2 | sort -k 8nr | grep Titrating_woman >>  titrating_genes_by_variance.woman.txt 
wc titrating_genes_by_variance.*.txt
\cp titrating_genes_by_variance.*.txt RESULTS.Stan


# count genes in the 1000 best woman by variance not seen in ACDB method
head -1000 titrating_genes_by_variance.woman.txt | cut -f 1 | sort > tata1
cat  titrating_genes.woman.txt  | cut -f 1 | sort > tata2
cat tata2 ZZZZZ tata1 | gawk '/ZZZZZ/{zz=1;next;}{if(zz<1)g[$1]=1;else {if(g[$1]<1)print;}}' | wc

# how mane woman genes are higher in 75 versus 100 => an alarming 1000 relative to 1700 ok
cat titrating_genes_by_variance.woman.txt | gawk '{n=$5 - $4;if(n<0)print}' | wc
cat titrating_genes_by_variance.woman.txt | gawk '{n=$5 - $4;if(n>0)print}' | wc



  echo "delta\tMan\tWoman"
  cat titrating_genes.*.txt | gawk '{x=int(2*$7+.49)/2;}/woman/{nf[x]++;nn[x]++;next;}/man/{nh[x]++;nn[x]++;}END{for(k in nn)print k,0+nh[k],0+nf[k];}' | sort -k 1n

end

# liste des read counts pour les samples ACDB de tous les genes, pour soustraire les quantiles
  cat  GeneIndexDB/gene2deep.$target.ace | gawk '/^Gene/{ok=1;next;}/G_Any/{ok=0;}/^Group_nU/{next}/sample_/{if(ok==1)printf ("%s\t%s\n",$2,$8);}' > _u
  cat _u | grep -v 100F | grep 0F | gawk '{printf("%s\t%d\n",$1,0+$2);}' | sort -k 2n > _u.0
  cat _u | grep 25F |  gawk '{printf("%s\t%d\n",$1,0+$2);}' | sort -k 2n > _u.25
  cat _u | grep 75F |  gawk '{printf("%s\t%d\n",$1,0+$2);}' | sort -k 2n > _u.75
  cat _u | grep 100F | gawk '{printf("%s\t%d\n",$1,0+$2);}' | sort -k 2n > _u.100

\rm _t6
foreach z (0 25 75 100)
  set tot=`cat _u.$z | gawk '{t+=$2;}END{print t}'`
  echo -n "$z"  >> _t6
   set tot1=`cat _u.$z | gawk '{if(1000*$2 < tot)t+=$2;}END{print t}' tot=$tot`
  echo -n "\t$tot1"  >> _t6
   set tot1=`cat _u.$z | gawk '{if(100*$2 < tot)t+=$2;}END{print t}' tot=$tot`
  echo -n "\t$tot1"  >> _t6
   set tot1=`cat _u.$z | gawk '{if(20*$2 < tot)t+=$2;}END{print t}' tot=$tot`
  echo -n "\t$tot1"  >> _t6
   set tot1=`cat _u.$z | gawk '{if(10*$2 < tot)t+=$2;}END{print t}' tot=$tot`
  echo -n "\t$tot1"  >> _t6
   set tot1=`cat _u.$z | gawk '{if(5*$2 < tot)t+=$2;}END{print t}' tot=$tot`
  echo  -n "\t$tot1"  >> _t6
  echo "\t$tot"  >> _t6
end
echo "%F\tExclude 1/1000\t1/100\t1/20\t1/10\t1/5\t1/1\t1/1000\t1/100\t1/20\t1/10\t1/5\t1/1" > _t7
cat _t6 | gawk '{printf("%s",$1);for(i=2;i<=7;i++)printf("\t%d",$i);for(i=2;i<=7;i++)printf("\t%.3f",log($i/$7)/log(2));printf("\n");}' >> _t7




foreach z (0 25 75 100)
  set tot=`cat _u.$z | gawk '{t+=$2;}END{print t}'`
  echo -n "$z\t$tot" 
  cat _u.$z | gawk '{if(5*$2 >= tot)print}' tot=$tot
end



#########################################################################################
#########################################################################################
# 

phaseii5:

goto phaseLoop

#########################################################################################
#########################################################################################
# export the tables, done in g4

phaseg6:


goto phaseLoop

#########################################################################################
#########################################################################################
# export the tables, somehow the output will be named GeneIndexDB/geneIndexTable.$target.gs.txt

phaseg7:

echo "g7: gene index statitics"

goto phaseLoop


#### variablity : genes with highest product index * variance
set toto=tmp/GENEINDEX/g7.gene_av_with_support.any.ace
set toto3=RESULTS/Expression/$MAGIC.most_variable_genes.txt
echo -n "# " > $toto3
date >> $toto3
echo "# Number of highly expressed AceView genes with most variable index" >> $toto3
cat  MetaDB/$MAGIC/GroupW_new_exonList MetaDB/$MAGIC/GroupList MetaDB/$MAGIC/RunsList ZZZZZ  MetaDB/$MAGIC/gtitle.txt ZZZZZ $toto  | gawk -F '\t' -f scripts/g7.most_variable_genes.awk  > $toto3.tmp

cat $toto3.tmp | head -15 | gawk '/^ZZZZZ/{zz=1;next;}{if(zz < 1)print}' >> $toto3
cat $toto3.tmp | gawk  '/^ZZZZZ/{zz=1;next;}{if(zz < 1)next;print}' | sort -k 2,2nr | head -1000 >> $toto3


#### 

set toto1=RESULTS/Expression/unique/$MAGIC.GeneIndex.av.u.txt
set toto2=RESULTS/Expression/unique/$MAGIC.GeneIndex.AceViewNovel.u.txt

cat $toto1 | gawk '/^Gene/{n++;if(n==2)ok=1;next;}{if(ok==1)print}' > toto3.tmp
cat $toto2 | gawk '/^Gene/{n++;if(n==2)ok=1;next;}{if(ok==1)print}' >> toto3.tmp
foreach group (`cat MetaDB/$MAGIC/GroupGene_element_indexList`)
  echo $group
  set toto3=RESULTS/Expression/1000genes_with_best_support/$MAGIC.$group.1000_av_genes_with_best_support.u.txt
  echo $group > $toto3
  cat $toto1 | head -30 | gawk '{print}/^Gene/{n++;if(n==2)exit}' >> $toto3
  
  set n=`cat $toto3 | gawk -F '\t' '{for(i=1;i<=NF ; i++)if($i == group " Index Unique")n=i;}END{print n+0 ;}' group=$group`

  if ($n < 2) continue 
  cat toto3.tmp | gawk -F '\t' '/^#/{next;}/^G_Any/{next}{x=0+$n;printf("%09d\t",int(100000000 - 1000*x));print}' n=$n | sort | gawk -F '\t' '{printf("%s",$2);for(i=3;i<=NF;i++){z=$i;if(z=="NA")z="";printf("\t%s",z);}printf("\n");}'  | head -1000 >> $toto3
  ls -ls $toto3
end

\rm toto3.tmp

goto phaseLoop

## R autocorrelation
# create a table

cat  MetaDB/$MAGIC/GroupList ZZZZZ tmp/GENEINDEX/g7.gene_av_with_support.any.ace | gawk '/^ZZZZZ/{zz=1;next;}{gsub(/\"/,"",$0);}{if(zz<1){n++;ok[$1]=n;printf("\t%s",$1);next;}}/^GeneId/{next;}/^Genefinder/{next;}/^Gene/{if(n8>0){printf("\n%s",gene);for(i=1;i<=n;i++){printf("\t%s",z[i]);z[i]=NN;}}gene=$2;n8=0;next;}/^Group_U/{k=ok[$2];if(k<1)next;x=$3;if(x<NN)x=NN;else n8++;z[k]=x;}' NN=9 > tmp/GENEINDEX/g7.group.R.txt

cat  MetaDB/$MAGIC/RunsList ZZZZZ tmp/GENEINDEX/g7.gene_av_with_support.any.ace | gawk '/^ZZZZZ/{zz=1;next;}{gsub(/\"/,"",$0);}{if(zz<1){n++;ok[$1]=n;printf("\t%s",$1);next;}}/^GeneId/{next}/^Genefinder/{next;}/^Gene/{if(n8>0){printf("\n%s",gene);for(i=1;i<=n;i++){printf("\t%s",z[i]);z[i]=NN;}}gene=$2;n8=0;next;}/^Run_U/{k=ok[$2];if(k<1)next;x=$3;if(x<NN)x=NN;else n8++;z[k]=x;}' NN=11 > tmp/GENEINDEX/g7.run.R.txt

R
genes = read.table("tmp/GENEINDEX/g7.group.R.txt")
gg=as.matrix(genes)
cc=cor(gg)
# heatmap takes 5 minutes to draw
pdf("tmp/GENEINDEX/g7.group.correl.pdf") ;
hh=heatmap(cc)
write(colnames(gg)[hh$rowInd],"tmp/GENEINDEX/g7.group.Legend.txt")
# pp=prcomp(genes,center=TRUE,scale=FALSE)
quit()

R
genes = read.table("tmp/GENEINDEX/g7.run.R.txt")
gg=as.matrix(genes)
cc=cor(gg)
# heatmap takes 5 minutes to draw
pdf("tmp/GENEINDEX/g7.run.correl.pdf") ;
hh=heatmap(cc)
write(colnames(gg)[hh$rowInd],"tmp/GENEINDEX/g7.run.Legend.txt")
# pp=prcomp(genes,center=TRUE,scale=FALSE)
quit()


#########################################################################################
## affy probesets
# gunzip -c ~/ftp/rat/MicroArrayProbeMapping_Affy.Rat230_2.2008_09_13.tar.gz | tar xf -
cd MicroArrayProbeMapping_Affy.Rat230_2.2008_09_13/  
cat Affy.Rat230_2.aceviewmrna.strand.txt | gawk -F '\t' '/^MRNA/{i=index($4,"_at");printf("AFFY\t%s\t%s\n",substr($1,6),substr($4,1,i+2));}' | sort -u > mrna2probeset
../../bin/tacembly ../../GeneIndexDB <<EOF
  query find mrna gene
  show -a -f mrna2gene.ace gene
  quit
EOF
cat mrna2gene.ace mrna2probeset | gawk '{gsub(/\"/,"",$0);}/mRNA/{m=$2;next;}/^Gene/{g[m]=$2;next;}/^AFFY/{m=$2;p=$3;if(g[m]){printf("Gene \"%s\"\nAffy \"%s\"\n\n",g[m],p);}}' > affy2gene.ace
cat mrna2gene.ace mrna2probeset | gawk '{gsub(/\"/,"",$0);}/mRNA/{m=$2;next;}/^Gene/{g[m]=$2;next;}/^AFFY/{m=$2;p=$3;if(g[m])next;{printf("Gene \"%s\"\nAffy \"%s\"\n\n",m,p);}}' >> affy2gene.ace

#########################################################################################
#########################################################################################
# compare different runs

phaseg8:

echo "g8: gene index statitics"
goto phaseLoop

bin/tacembly MetaDB <<EOF
  query find run compare_to AND project == $MAGIC 
  show -a -f GeneIndexDB/g7.r2c.ace compare_to  
  query find run project == $MAGIC ; > Ali ; Low_index
  bql -a -o GeneIndexDB/g7.r2c.low_index.txt  select  r,i from r in @, i in r->Low_index
EOF

bin/tacembly GeneIndexDB <<EOF
  query find gene deep && affy
  bql -a -o GeneIndexDB/g7.r2c.affy.txt select  g,a from g in @, a in g->Affy
EOF
  
cat GeneIndexDB/g7.r2c.ace ZZZZZ GeneIndexDB/g7.r2c.ace   | gawk '/^ZZZZZ/{zz=1;next;}{gsub(/\"/,"",$0);}/^Run/{ z=$2 ; r = r2i[z] ; if (zz < 1) { if(r<1){nr++;r=nr;i2r[nr]=z;r2i[z]=nr;} ; next; }}/^Compare_to/{z=$2 ; r2 = r2i[z] ;if (r<r2)printf ("%s/%s\n",i2r[r],i2r[r2]);}' > GeneIndexDB/g7.r2c.pairs


set toto="GeneIndexDB/$MAGIC.geneTable.u.av.*eneId.Index.txt"
foreach r2c (`cat GeneIndexDB/g7.r2c.pairs`)
  set nam1=`echo $r2c | gawk '{split($1,aa,"/");print aa[1];}'`
  set nam2=`echo $r2c | gawk '{split($1,aa,"/");print aa[2];}'`

  set toto3=GeneIndexDB/g7.genes_differentially_expressed_between.$nam1.$nam2.txt
  set NN=`cat GeneIndexDB/g7.r2c.low_index.txt | gawk '{gsub(/\"/,"",$0);if($1==nam1 || $2==nam2){if(x<1 || $2<x)x=$2;}}END{print x}'  nam1=$nam1 nam2=$nam2`

  echo "$nam1 <--> $nam2  NN=$NN"

    echo -n "# " > $toto3
    date >> $toto3
    echo "# Genes with significantly different support between $nam1 and $nam2" >> $toto3
    echo "Only genes with one expression index above 9 and at least 4 fold differences are listed" >> $toto3
    echo toto | gawk -F '\t' '/^Gene/{printf("\n%s",$1);for (i=2;i<=10;i++)printf("\t%s",$i); printf("\t%s\t%s\tIndex variation\n",nam1,nam2);}' nam1=$nam1 nam2=$nam2  >> $toto3
    cat  $toto  | gawk -F '\t' '/^kb/{next}/GC in run/{next}/^Gene/{ok++;gsub(/ Index Unique/,"",$0);for(i=1;i<=NF;i++){if($i==nam1)u1=i;if($i==nam2)u2=i;}next;}{if(ok<2)next;}{a=0+$u1;b=0+$u2;if(a<NN)a=NN;if(b<NN)b=NN;{da=a-b;if(da<0){da=-da;}if(da>=2){printf("%.1f\t%s",a-b,$1);for (i=2;i<=10;i++)printf("\t%s",$i);printf("\t%s\t%s\t%.1f\n",a,b,a-b);}}}' nam1=$nam1 nam2=$nam2 NN=$NN | sort -k 1,1nr | gawk -F '\t' '{x="";for(i=2;i<=NF;i++){printf("%s%s",x,$i);x="\t";}printf("\n");}' >> $toto3

  \cp $toto3  RESULTS/Expression
  cat $toto3 | gawk '/#/{next;}{n[int(10*(.5+$13))/10]++;}END{for(k in n)if(int(k) != 0)printf("%d\t%d\n", k,n[k])}' | sort -k 1nr >  RESULTS/Expression/g7.genes_differentially_expressed_between.$nam1.$nam2.stats.txt
end


goto phaseLoop


#########################################################################################
#########################################################################################
## Global intron support statistics

phaseii2d:
# collect each type
if (! -d tmp/INTRON_INDEX) goto phaseLoop
if (! -d tmp/introns) mkdir tmp/introns
if (! -d RESULTS/Introns_exons_polyA) mkdir  RESULTS/Introns_exons_polyA

if (! -e tmp/introns/$MAGIC.introns.stranding.ace) then
  echo "--- ii2d: introns.stranding"
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    if (-e tmp/INTRONRUNS/$run/$run.u.intronSupport.ace.gz) then
      gunzip -c  tmp/INTRONRUNS/$run/$run.u.intronSupport.ace.gz | head -3 | gawk '/^Anti_run/{a=$6;next;}/^Run_U/{b=$6;next;}END{if(a+b>10000)printf("Ali %s\nStranding Exons_juntions %.2f  %d plus %d minus 0 ambiguous\n\n",run, 100*a/(a+b), b, a);}' run=$run >> tmp/introns/$MAGIC.introns.stranding.ace
    endif
  end
  echo "pparse tmp/introns/$MAGIC.introns.stranding.ace tmp/INTRON_INDEX/Exons_juntions.stranding.ace" | bin/tacembly MetaDB -noprompt 
  bin/bestali -groupLetterProfile -db MetaDB -project $MAGIC
endif

goto phaseLoop


# histo in true log2 scale
set toto1=tmp/introns/introns.de_uno.count
set toto=RESULTS/Introns_exons_polyA/$MAGIC.de_uno_support_of_exons_junctions_discovered_in_the_genome.histo.txt
set toto1=tmp/introns/$MAGIC.cumulated_support.txt
set toto=RESULTS/Expression/unique/$MAGIC.cumulated_support_of_annotated_exons_junctions.histo.txt

if (! -e tmp/introns/$MAGIC.cumulated_support.txt) then
  echo "--- ii2d: cumulated_support per run"
  echo ' ' > tmp/introns/$MAGIC.intronSupportPerRun.txt
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    if (-e tmp/INTRONRUNS/$run/$run.u.intronSupport.ace.gz) then
      gunzip -c tmp/INTRONRUNS/$run/$run.u.intronSupport.ace.gz | gawk '/^Intron/{g=$2;next;}/^Run_U/{pgg[g]+=$6;gg[g]+=$6;next;}/^Anti_run/{agg[g]+=$6;gg[g]+=$6;next;}END{for(g in gg)printf ("%s\t%s\t%d\t%d\t%d\n",run, g,gg[g],pgg[g],agg[g]);}' run=$run >>  tmp/introns/$MAGIC.intronSupportPerRun.txt
    endif
  end

 cat  tmp/introns/$MAGIC.intronSupportPerRun.txt | gawk -F '\t' '{gg[$2]+=$3;}END{for(g in gg)if(gg[g]>0)printf ("%s\t%d\n", g,gg[g]);}' | sort -k 2nr > tmp/introns/$MAGIC.cumulated_support.txt
 \rm tmp/introns/$MAGIC.intronSupportPerRun.txt 
 if (-e $toto) \rm $toto 
endif


if (! -e $toto) then

  if (-e tmp/introns/$MAGIC.introns.any.txt) \rm tmp/introns/$MAGIC.introns.any.txt
  if (-e tmp/introns/$MAGIC.introns.annot.txt) \rm tmp/introns/$MAGIC.introns.annot.txt
  foreach target ($Etargets introns)
    if (! -e  tmp/introns/introns.$target.txt) then
      cat TARGET/MRNAS/introns_$target.txt | gawk '{printf("%s__%d_%d\n", $1,$2,$3);}' > tmp/introns/introns.$target.txt
    endif
    if ($target != introns) cat tmp/introns/introns.$target.txt >> tmp/introns/$MAGIC.introns.annot.txt
    cat tmp/introns/introns.$target.txt >> tmp/introns/$MAGIC.introns.any.txt
  end

  echo '#' > $toto.1
  foreach target ($Etargets introns)

     echo "--- ii2d: cumulated_support per target : $target"
     cat  tmp/introns/introns.$target.txt ZZZZZ $toto1  | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$1]++;if(gg[$1]==1)n1++;next;}}{ok=0;g=$1;gsub(/\"/,"",g);ok=gg[g];if(ok<1)next;n=$2 ;if(n==0)next;if(ggn[g]<1)n2++;ggn[g]+=n; k=n;j=0;m=1;while(k>0){j++;k=int(k/2);m=2*m;}j--;nn[j,1]++;if(j>jmax)jmax=j;}END{printf("%d\t%d\n",0,n1-n2);m=1;for(j=0;j<=jmax;j++){printf("%d\t%d\n",m,nn[j,1]);m=2*m;}}' | sort -k 1,1nr | gawk '{if(zz<1)printf("# Level\t%s\t%s Cumul\n",tt,tt);zz=1;}{n+=$2;printf("Intron\t%s\t%d\t%d\t%d\n",tt,$1,$2,n);}' tt=$target >> $toto.1
  end
  set annot=`echo "$Etargets" | gawk '{gsub(/ /,"_",$0);print}'`
  cat tmp/introns/$MAGIC.introns.annot.txt ZZZZZ $toto1 | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$1]++;if(gg[$1]==1)n1++;next;}}{ok=0;g=$1;gsub(/\"/,"",g);ok=gg[g];if(ok<1)next;n=$2 ;if(ggn[g]<1)n2++;ggn[g]+=n; if(n==0)next;k=n;j=0;m=1;while(k>0){j++;k=int(k/2);m=2*m;}j--;nn[j,1]++;if(j>jmax)jmax=j;}END{printf("%d\t%d\n",0,n1-n2);m=1;for(j=0;j<=jmax;j++){printf("%d\t%d\n",m,nn[j,1]);m=2*m;}}' | sort -k 1,1nr | gawk '{if(zz<1)printf("# Level\t%s\t%s Cumul\n",tt,tt);zz=1;}{n+=$2;printf("Intron\t%s\t%d\t%d\t%d\n",tt,$1,$2,n);}' tt=$annot >> $toto.1
  cat tmp/introns/$MAGIC.introns.any.txt ZZZZZ $toto1 | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$1]++;if(gg[$1]==1)n1++;next;}}{ok=0;g=$1;gsub(/\"/,"",g);ok=gg[g];if(ok<1)next;n=$2 ;if(ggn[g]<1)n2++;ggn[g]+=n; if(n==0)next;k=n;j=0;m=1;while(k>0){j++;k=int(k/2);m=2*m;}j--;nn[j,1]++;if(j>jmax)jmax=j;}END{printf("%d\t%d\n",0,n1-n2);m=1;for(j=0;j<=jmax;j++){printf("%d\t%d\n",m,nn[j,1]);m=2*m;}}' | sort -k 1,1nr | gawk '{if(zz<1)printf("# Level\t%s\t%s Cumul\n",tt,tt);zz=1;}{n+=$2;printf("Intron\t%s\t%d\t%d\t%d\n",tt,$1,$2,n);}' tt=Union >> $toto.1


  cat  $toto.1 | gawk -F '\t'  '/^#/{next;}{z=$3;x=$1 " " $2;zz[z]=1;xx[x]=1;nn1[z,x]=$4;nn2[z,x]=$5;}END{for(x in xx)printf("\tcount %s\tcumul %s",x,x);for(z in zz){printf("\n%s",z);for(x in xx)printf("\t%s\t%s",nn1[z,x],nn2[z,x]);}printf("\n");}' > $toto.2
head -1 $toto.2 > $toto.3
tail -n +2 $toto.2 | sort -k 1n >> $toto.3

  echo -n '# ' > $toto
  date >> $toto
  echo "# Cumulated support of annotated exon junctions in project $MAGIC" >> $toto

  cat $toto.3 | scripts/transpose | head -1 > $toto.4
  cat $toto.3 | scripts/transpose | sort | grep count | grep -v Union >> $toto.4
  cat $toto.3 | scripts/transpose | sort | grep count | grep Union >> $toto.4
  cat $toto.3 | scripts/transpose | sort | grep cumul | grep -v Union >> $toto.4
  cat $toto.3 | scripts/transpose | sort | grep cumul | grep Union >> $toto.4
  cat $toto.4 | scripts/transpose >> $toto
 
  echo -n " Cumulated support of annotated exon junctions : "
  echo $toto
  \rm $toto.[1234]
endif

goto phaseLoop

####
## table of counts limited to ebi/RefSeq/av

cat  ZZZZZ totox.RefSeq.txt   ZZZZZ  totox.EBI.txt  ZZZZZ   totox.av.txt  ZZZZZ | gzip > _toto.gz
echo -n "# " > RESULTS/test.introns.stats.txt
date >>  RESULTS/test.introns.stats.txt
foreach mm (A B C D)
   gunzip -c RESULTS/Expression/unique/introns/SEQC_$mm.introns.INTRON.u.reads_aligned_per_gene.txt.gz | head -60 | gawk '/^# /{print}/^#xxxxExon/{next;}/^#/{n[$1]++;if(n[$1]==1)print}' >> RESULTS/test.introns.stats.txt
  gunzip -c RESULTS/Expression/unique/introns/SEQC_$mm.introns.INTRON.u.reads_aligned_per_gene.txt.gz | head -60 | gawk '/^# /{print}/^#Exon/{next;}/^#/{n[$1]++;if(n[$1]==1)print}' > RESULTS/Expression/unique/introns/SEQC_$mm.introns_RefSeq_EBI_AceView.INTRON.u.reads_aligned_per_gene.txt
echo "#The refseq transcript ID column actually shows in which annotation the junction is present"  >> RESULTS/Expression/unique/introns/SEQC_$mm.introns_RefSeq_EBI_AceView.INTRON.u.reads_aligned_per_gene.txt
echo "#The _SumOfAllReadsInProject shows the total number of reads supporting the junction in sample_$mm"  >> RESULTS/Expression/unique/introns/SEQC_$mm.introns_RefSeq_EBI_AceView.INTRON.u.reads_aligned_per_gene.txt
  gunzip -c  _toto.gz RESULTS/Expression/unique/introns/SEQC_$mm.introns.INTRON.u.reads_aligned_per_gene.txt.gz | gawk -F '\t' '/^ZZZZZ/{if(zz<1)zz=1;else zz=2*zz;next;}{if(zz<8){gg[$1]+=zz;next;}}/^#/{next;}{ok=gg[$1];if(ok>0){ printf("%s\t",$1);if(ok%2==1)printf("RefSeq,");if(int(ok/2)%2==1)printf("Encode,");if(int(ok/4)%2==1)printf("AceView,");for(i=3;i<=NF;i++)printf("\t%s",$i);printf("\n");}}' | sort -k 8nr  >> RESULTS/Expression/unique/introns/SEQC_$mm.introns_RefSeq_EBI_AceView.INTRON.u.reads_aligned_per_gene.txt
  \rm  RESULTS/Expression/unique/introns/SEQC_$mm.introns_RefSeq_EBI_AceView.INTRON.u.reads_aligned_per_gene.txt.gz
  gzip RESULTS/Expression/unique/introns/SEQC_$mm.introns_RefSeq_EBI_AceView.INTRON.u.reads_aligned_per_gene.txt
end

   cat  totox.av.txt  totox.seqc.txt ZZZZZ  toto.introns.support.ace  | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$1]=1;next;}}/^Intron/{ok=0;z=$2;gsub(/\"/,"",z);ok=gg[z];}/^Group_U/{if(ok<1)next;n=$6 ; if(n==0)next;k=n;j=0;m=1;while(k>0){j++;k=int(k/10);m=10*m;}x=int(10*n/m);if(x>=5)x=5;else{ if(x>=2)x=2;else x=1;}j--;nn[j,x]++;if(j>jmax)jmax=j;}END{m=1;for(j=0;j<=jmax;j++){for(y=0;y<3;y++){x=1;if(y==1)x=2;if(y==2)x=5;printf("%d\t%d\n",x*m,nn[j,x]);}m=10*m;}}' | sort -k 1,1nr | gawk 'BEGIN{printf("Level\t%s\t%s Cumul\n",tt,tt);}{n+=$2;printf("%d\t%d\t%d\n",$1,$2,n);}' tt=$target > tatax.av_seqc

# venn of introns
 cat toto.introns.support.ace ZZZZZ  totox.av.txt ZZZZZ totox.encode.txt ZZZZZ totox.RefSeq.txt | gawk '/ZZZZZ/{if(zz<1)zz=1;else zz=2*zz;next;}/^Intron/{if(zz<1){g=$2;gsub(/\"/,"",g);ok[g]=1;next;}}/^Group/{next;}{if(zz>0){gg[$1]+=zz;if(ok[$1]>0)ggok[$1]+=zz;}}END{for(k in gg){nn[0]++;nn[gg[k]]++;nnok[ggok[k]]++;}for (n in nn)printf("%d\t%d\t%d\t%.1f\n", n,nn[n],nnok[n],100*nnok[n]/nn[n]);}'| sort -k 1n

cat totox.NB22.txt  ZZZZZ  toto.introns.support.ace  | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$1]=1;next;}}/^Intron/{ok=0;z=$2;gsub(/\"/,"",z);ok=gg[z];}/^Group_U/{if(ok<1)next;n=$6 ; if(n==0)next;if(n==10)print z;}' | head

cat  totox.av.txt ZZZZZ  totox.seqc.txt ZZZZZ  toto.introns.support.ace  | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){av[$1]=1;next;}}{if(zz<2){seqc[$1]=1;next;}}/Intron/{ok=0;z=$2;gsub(/\"/,"",z);ok=av[z];ok2=seqc[z];}/^Group_U/{if(ok==1 || ok2==0)next;if($6>10000)print z,$6;}'

goto phaseLoop

#########################################################################################

phaseLoop:
echo "Normal exit"
exit 0

#########################################################################################
#########################################################################################
#########################################################################################

## Table limited to a few genes
mkdir RESULTS/TRPV
foreach GM (GENE MRNAH)
  foreach type (expression_index reads_aligned_per_gene sFPKM)
    gunzip -c RESULTS/Expression/unique/av/$MAGIC.AceView.$GM.u.$type.txt.gz | gawk '/^#/{print;next;}/^TRPV/{print;next;}' > RESULTS/TRPV/$MAGIC.AceView.$GM.u.$type.txt
  end
end



pushd RESULTS/Expression/unique/RefSeq
  mkdir MATT
  foreach ff (`ls DRG_SRA_AxotomyIada.RefSeq.GENE.u.*`)
    zcat -f $ff | gawk '/^#/{print;next;}{ok=0;}/Pdyn/{ok=1;}/Penk/{ok=1;}/Pomc/{ok=1;}/Pnoc/{ok=1;}/Gal/{for(i=1;i<NF;i++)if($i=="Gal")ok=1;}{if(ok==1)print}' > MATT/$ff
  end
popd  

  

cat tmp/GENEINDEX/DD.EBI.GENE.u.ace | gawk '{gsub(/\"/,"",$2);}/^Gene/{g=$2;ok=0;if(g=="FBgn0002938")ok=1;next;}/^Run_U/{r=$2;if(ok==1)z1[r]=$4;next;}/^Anti/{r=$2;if(ok==1)z2[r]=$4;next;}END{for(r in z2);if(z2[r]>0)printf("%s\t%d\t%d\n",r,z1[r],z2[r]);}'

## counting SNPs in 2 runs
gunzip -c tmp/SNP_DB/Who.characteristic.snp.gz | grep Rhs4005 > toto80                                  
gunzip -c tmp/SNP_DB/Who.characteristic.snp.gz | grep F2_MCS-015 >> toto80                                  


cat toto81 | gawk -F '\t' '/A>G/{next;}{k++;if(k%1)next;}/>/{x=$8;if($7<.4)x=0;if($7>1.8)x=100;if($9<10)next;if(x<=4)x=0;else if(x>=98)x=200;else x=100;}{if($6 != old6)f2++;old6=$6;}{s=$1 $2 $3;ss[s]++;z[s,f2]=1+x/100;}END{for(s in ss) {x = 0+z[s,1] ; y = 0+z[s,2] ; if(x*y==0 || x*y == 100 || x == 20 || y == 20 )continue;uu[x-1,y-1]++ ; X += x ; X2 += x*x ; Y += y ; Y2 += y * y ; XY += x*y ;N++; if(x*y==60 || x*y==9)print s,x,y; }printf("%d\t%d\t%d\n%d\t%d\t%d\n%d\t%d\t%d\n\n",uu[0,0],uu[0,1],uu[0,2],uu[1,0],uu[1,1],uu[1,2],uu[2,0],uu[2,1],uu[2,2]);x = X2/N - X*X/(N*N);y=Y2/N - Y*Y/(N*N);w=XY/N - X*Y/(N*N);printf("\nCorrel %.2f\n",100*w/sqrt(x*y));print X,Y,X2,Y2,x,y,w;}' | tail -20




exit 0

# hack 2022_06_12    pour rattrapper le fait que l'histo des lg des fragments dans bestali avait une erreur ,introduite en juillet 2021 (covid)
foreach run (`cat toto.bad.runs`)
  foreach lane (`cat Fastc/$run/LaneList`)
    scripts/submit tmp/COUNT/$lane.pair_stats_fix "bin/bestali -target_class ET_av -i tmp/COUNT/$lane.hits.gz -run $run -pair 500 -seqc  -strategy RNA_seq -o tmp/COUNT/$lane.fix"
  end
end

foreach run (`cat toto.bad.runs`)
  foreach lane (`cat Fastc/$run/LaneList`)
    scripts/submit tmp/COUNT/$lane.pair_stats_fix "bin/bestali -target_class ET_av -i tmp/COUNT/$lane.hits.gz -run $run -pair 500 -seqc  -strategy RNA_seq -o tmp/COUNT/$lane.fix"
  end
end
foreach run (`cat toto.bad.runs`)
  foreach lane (`cat Fastc/$run/LaneList`)
    cat tmp/COUNT/$lane.pairStats | head -25 > tmp/COUNT/$lane.fix2
    cat tmp/COUNT/$lane.fix.pairStats | tail -2 >> tmp/COUNT/$lane.fix2
  end
end
foreach run (`cat toto.bad.runs`)
  foreach lane (`cat Fastc/$run/LaneList`)
    \mv  tmp/COUNT/$lane.fix2  tmp/COUNT/$lane.pairStats 
  end
end




exit 0

# analyse des ali de type -2, pour voir la distance sur le genome
# dans pacbio on pense que c'est peut etre une fausse interpretation de la strand dans la boucle de sequencage

# en fait l'immense majorite represente deux reads aligne a la meme position x1,x2, donc deux annotations antisense avec un sequencage non-strande

set run=PacBio2.ccs3_AGLR2_A-F1_run1
zcat tmp/COUNT/$run/f2.*.hits.gz | gawk -F '\t' '/ET_av/{if ( $10 == -2 ) { gsub("<",">",$1) ; }' | cut -f 1,6,7,9 > toto1
cat toto1 | gawk -F '\t' '{if($1!=r){aa=$0;r=$1;x1=$2;x2=$3;g1=$4;next;}if(x1==$2 && x2==$3)next;print aa;print}' > toto2
cat toto1 | gawk -F '\t' '{if($1!=r){aa=$0;r=$1;x1=$2;x2=$3;g1=$4;next;}if(x1==$2 && x2==$3)next;print aa;print}' > toto2
cat toto1 | gawk -F '\t' '{if($1!=r){aa=$0;r=$1;x1=$2;x2=$3;g1=$4;next;}if(x1==$2 && x2==$3)next;if(g1=="MALAT1" || $4=="MALAT1")print aa;print}' > toto3
wc toto[123]
cat toto2 | cut -f 4 | tags | sort -k 2nr | head -12


