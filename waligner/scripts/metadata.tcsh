#!bin/tcsh -f
# Create in a lazy way all needed metadafiles
# in order to minimize the number of files needed in TARGET

if (! -d tmp/METADATA) mkdir tmp/METADATA
if (-e tmp/METADATA/metadata.ok) \rm tmp/METADATA/metadata.ok

if (! -e tmp/METADATA/v2024_10_10) \rm  tmp/METADATA/*
touch  tmp/METADATA/v2024_10_10

echo "=== check the TARGET data"

set allRNAtargets=`echo $Etargets $Ttargets $RNAtargets | gawk '{for (i=1 ; i<= NF ; i++){if (z[$i]<1){z[$i]=1;zz = zz " " $i;}}}END{print zz}'`

######################################################
# verify the existence of the target fasta files

echo "---- verify the existence of the target fasta files"
foreach target ($DNAtargets)
  if ($target == gdecoy) continue
  # echo "------- verify $target"
  if (! -e TARGET/Targets/$species.$target.fasta.gz) then
    echo "\nFATAL ERROR : TARGET/Targets/$species.$target.fasta.gz missing"
    echo "Please remove this target from setenv DNAtargets in ./LIMITS or add the target fasta file"
    goto done
  endif
end

foreach target ($allRNAtargets)
  # echo "------- verify $target"
  if (! -e TARGET/Targets/$species.$target.fasta.gz && ! -e TARGET/GTF/$species.$target.gff.gz ) then
    echo "\nFATAL ERROR : TARGET/Targets/$species.$target.fasta.gz missing"
    echo "Please remove this target from setenv RNAtargets in ./LIMITS or add the target fasta file"
    goto done
  endif

  if ($target == rrna) continue
  if ($target == gdecoy) continue
  if ($target == gdecoy) continue

  if (0 && ! -e TARGET/GTF/$species.$target.gtf.gz) then
    echo "---- WARNING: Missing file TARGET/GTF/$species.$target.gtf.gz"
  endif
end

######################################################
######################################################
# compute the length of the DNA references based on the fasta file

echo "---- compute the length of the DNA references based on the fasta file"
foreach target ($DNAtargets)
  if ($target == gdecoy) continue
  if (! -e  tmp/METADATA/$species.$target.TM.txt.gz) then
    echo "---- measuring the length of each sequence in DNA target : $target"
    (bin/dna2dna -i  TARGET/Targets/$species.$target.fasta.gz  -getTm -o  tmp/METADATA/$species.$target -gzo) >& tmp/METADATA/$target.TM.err
    gunzip -c tmp/METADATA/$species.$target.TM.txt.gz | wc
  endif
end

if (-e tmp/METADATA/$species.genome.TM.txt.gz && ! -e tmp/METADATA/genome.ns.sponge1) then
  gunzip -c  tmp/METADATA/$species.genome.TM.txt.gz | gawk -F '\t' '/^#/{next;}{split($1,aa," ");printf("%s\t1\t%s\t1\t%d\t%s\n",aa[1],aa[1],$2,aa[1]);}' >  tmp/METADATA/genome.ns.sponge.1
  gunzip -c  tmp/METADATA/$species.genome.TM.txt.gz | gawk -F '\t' '/^#/{next;}{split($1,aa," ");printf("%sR\t1\t%s\t%d\t1\t%sR\n",aa[1],aa[1],$2,aa[1]);}' >>  tmp/METADATA/genome.ns.sponge.1
  if (-e tmp/METADATA/genome.ns.sponge) \rm tmp/METADATA/genome.ns.sponge
  foreach chrom ($chromSetAll)
    cat tmp/METADATA/genome.ns.sponge.1  | gawk  '{if($1 == chrom || $1 == chrom "R")print}' chrom=$chrom >> tmp/METADATA/genome.ns.sponge
  end
  \rm  tmp/METADATA/genome.ns.sponge.1 
endif

######################################################
# compute the length of the transcripts based on the fasta file

echo "---- compute the length of the transcripts based on the fasta file"
set new_fasta=0
foreach target ($allRNAtargets)
  if (-e  tmp/METADATA/$target.mrna_ln_gc_gene_geneid.txt) continue
  
  echo "---- measuring the length of each sequence in RNA target : $target"

  set new_fasta=1
  if (! -e  tmp/METADATA/$species.$target.TM.txt.gz) then
    (bin/dna2dna -i  TARGET/Targets/$species.$target.fasta.gz  -getTm -gzo -o tmp/METADATA/$species.$target) >>& tmp/METADATA/$target.TM.err
  endif
  gunzip -c tmp/METADATA/$species.$target.TM.txt.gz | gawk -F '\t' 'BEGIN{printf("#Transcript\tLength\tGC_percent\tGene\tGeneId\n");}/^#/{next;}{gsub(/>/,"",$1);split($1,aa,"|");t=aa[1];t=aa[1];g="";gid="";if(aa[2]=="Gene" || aa[2]=="GENE")g=aa[3];if(aa[4]=="GeneId")gid=aa[5];ln=$2;gc=$5;printf("%s\t%s\t%d\t%s\t%s\n",t,ln,gc,g,gid);}' | sort > tmp/METADATA/$target.mrna_ln_gc_gene_geneid.txt

  cat  tmp/METADATA/$target.mrna_ln_gc_gene_geneid.txt |   gawk -F '\t' '/^#/{next;}{m=$1;ln=$2;gc=$3;g=$4;gid=$5;if(length(m)<1)next;printf("mRNA \"%s\"\nTargeted\nLength %d\n",m,ln);if(gc+0>0)printf("GC_percent %d\n",gc);if(length(g)>1)printf("Gene \"%s\"\n",g);if(length(gid)>1)printf("GeneId \"%s\"\n",gid);printf("\n");}' >  tmp/METADATA/$target.MRNA.ln.ace
end

# identification des longueurs des genes   2014_06_07   from april to this date we used the length of the .a
# we switch back  to length of longest as used in february 2014 for NB export
# using the first variant of each gene, for AceView the .a
    
if (-e tmp/METADATA/$MAGIC.mrna_ln_gc_gene_geneid.txt) \rm tmp/METADATA/$MAGIC.mrna_ln_gc_gene_geneid.txt
foreach target ($Etargets)
    cat  tmp/METADATA/$target.mrna_ln_gc_gene_geneid.txt  >> tmp/METADATA/$MAGIC.mrna_ln_gc_gene_geneid.txt
    cat  tmp/METADATA/$target.mrna_ln_gc_gene_geneid.txt |  gawk -F '\t' '/^#/{next;}{g=$4;ln=$2;gc=$3;gid=$5;if(length(g)<1)next;if(length(gid)>0)g2gig[g]=gid;if(ln>g2ln[g] && g2n[g]<5){g2ln[g]=ln;g2gc[g]=gc;}g2n[g]++;}END{for(g in g2ln) { printf("Gene \"%s\"\n%s\nTargeted\nLength %d\nGC_percent %d\n",g,tag,g2ln[g],g2gc[g]);if(g2gid[g])printf("GeneId \"%s\"\n",g2gid[g]);printf("\n");}}' tag=$target >  tmp/METADATA/$target.GENE.ln.ace
end

######################################################
# Initialize the info.ace file with the data derived from the .fasta file

if (1 || $new_fasta == 1) then
  foreach target ($allRNAtargets)
    cat  tmp/METADATA/$target.GENE.ln.ace >  tmp/METADATA/$target.GENE.info.ace
    cat  tmp/METADATA/$target.MRNA.ln.ace >  tmp/METADATA/$target.MRNA.info.ace
  end
endif

######################################################
######################################################
# Process the GTF files (needed for Exome interpretation of CNV and SNP)

echo "---- process the gtf files"
set gtf_active=0
foreach target ($allRNAtargets)
  source scripts/target2target_class.txt

  if (-e TARGET/GTF/$species.$target.gtf.gz || -e TARGET/GTF/$species.$target.gff.gz || -e TARGET/GTF/$species.$target.gff.Bacteria.gz) then
    if (-e tmp/METADATA/$target.gtf.done) continue

    echo "---- processing the gtf file of RNA target : $target"

    if (-e TARGET/GTF/$species.$target.gtf.gz) then
      echo " bin/dna2dna -gtf TARGET/GTF/$species.$target.gtf.gz  -gtfRemap $target_class  -o tmp/METADATA/gtf.$target "
            (bin/dna2dna -gtf TARGET/GTF/$species.$target.gtf.gz  -gtfRemap $target_class  -o tmp/METADATA/gtf.$target ) >& tmp/METADATA/gtf.$target.err
      \rm tmp/METADATA/gtf.$target.*gz
      gzip  -f tmp/METADATA/gtf.$target.*
    else if (-e TARGET/GTF/$species.$target.gff.Bacteria.gz) then
      (bin/dna2dna -gff3 TARGET/GTF/$species.$target.gff.Bacteria.gz  -gtfRemap $target_class  -o tmp/METADATA/gtf.$target ) >& tmp/METADATA/gtf.$target.err
      \rm tmp/METADATA/gtf.$target.*gz
      gzip  -f tmp/METADATA/gtf.$target.*
    else if (0 && -e TARGET/GTF/$species.$target.gff.gz) then
      (bin/dna2dna -gff3 TARGET/GTF/$species.$target.gff.gz  -gtfRemap $target_class  -o tmp/METADATA/gtf.$target ) >& tmp/METADATA/gtf.$target.err
      \rm tmp/METADATA/gtf.$target.*gz
      gzip  -f tmp/METADATA/gtf.$target.*
      # remove the parasite .fasta file that should not have been created 
      if (-e  tmp/METADATA/gtf.$target.fasta.gz) \rm  tmp/METADATA/gtf.$target.fasta.gz 

      if (! -e TARGET/Targets/$species.$target.fasta.gz && -e  TARGET/Targets/$species.genome.fasta.gz) then
        # we need to order the mrna file, this can only done inside acedb
        mkdir tmp/Sort.$target
        pushd tmp/Sort.$target
          mkdir database
          ln -s ../../metaData/wspec.aceview_web_site wspec
          tace  . <<EOF
y
            parse ../../TARGET/Targets/$species.genome.fasta.gz
            parse ../../TARGET/Targets/$species.mito.fasta.gz
            parse ../../TARGET/Targets/$species.chloro.fasta.gz
            parse ../METADATA/gtf.$target.transcripts.ace.gz
            query find sequence Gene_model
            select -o mrnas.txt g,m,d from g in @, m in g->gene_model, d in DNA(m)
            save
            quit
EOF
          cat  mrnas.txt | gawk -F '\t' '{gsub(/^X__/,"",$1);printf(">%s|Gene|%s\n%s\n",$2,$1,$3);}' | gzip  > $species.$target.fasta.gz
          mv $species.$target.fasta.gz ../../TARGET/Targets
        popd
        #\rm -rf tmp.Sort.$target
      endif
    endif
    # ATTENTION -gzo crashes on mouse (although it works on rat), we must gzip as a postpocessing, probaly the files are too large
    if (-e  tmp/METADATA/gtf.$target.r.sponge.gz) then
      gunzip -c tmp/METADATA/gtf.$target.f.sponge.gz | gawk -F '\t' '{if($1==old)printf("Intron %s__%d_%d\n%s\nSupports %s\n\n",$3,a2+1,$4-1,target,$1);old=$1;a2=$5;}' target=$target > tmp/METADATA/gtf.$target.f.intron.ace
      gunzip -c tmp/METADATA/gtf.$target.r.sponge.gz | gawk -F '\t' '{if($1==old)printf("Intron %s__%d_%d\n%s\nSupports %s\n\n",$3,a2-1,$4+1,target,$1);old=$1;a2=$5;}' target=$target > tmp/METADATA/gtf.$target.r.intron.ace
    endif



    gunzip -c  tmp/METADATA/gtf.$target.mrnaRemap.gz | gawk -F '\t' '{m=$2;c[m]=$5;a1=$6;a2=$7;if(a1>a2){a0=a1;a1=a2;a2=a0;isUp[m]=1;}if(m1[m]+0==0){m1[m]=a1;m2[m]=a2;}if(a1<m1[m])m1[m]=a1;if(a2>m2[m])m2[m]=a2;}END{for(m in c){printf("%s\t%s",m,c[m]);if(isUp[m]==1)printf("\t%d\t%d\n",m2[m],m1[m]);else printf("\t%d\t%d\n",m1[m],m2[m]);}}' | sort >  tmp/METADATA/gtf.$target.mrna2intMap.txt 
    gunzip -c  tmp/METADATA/gtf.$target.mrnaRemap.gz | gawk -F '\t' '{g=$8;c[g]=$5;a1=$6;a2=$7;if(a1>a2){a0=a1;a1=a2;a2=a0;isUp[g]=1;}if(g1[g]+0==0){g1[g]=a1;g2[g]=a2;}if(a1<g1[g])g1[g]=a1;if(a2>g2[g])g2[g]=a2;}END{for(g in c){printf("%s\t%s\t",g,c[g]);if(isUp[g]==1)printf("%d\t%d\n",g2[g],g1[g]);else printf("%d\t%d\n",g1[g],g2[g]);}}' | sort >  tmp/METADATA/gtf.$target.gene2intMap.txt 

# reject geneboxes longer than geneBoxMaxLn
set geneBoxMaxLn=1000000
#    cat tmp/METADATA/gtf.$target.gene2intMap.txt | gawk -F '\t' '{a1=$3;a2=$4;if(a1>a2)next;n=a2-a1;if(n<0){n=-n;} if (n<nMax)printf("%s\t1\t%s\t%d\t%d\t%s\n",$1,$2,a1,a2,$1);}' nMax=$geneBoxMaxLn > tmp/METADATA/$target.f.gene.sponge
#    cat tmp/METADATA/gtf.$target.gene2intMap.txt | gawk -F '\t' '{a1=$3;a2=$4;if(a1<a2)next;n=a2-a1;if(n<0){n=-n;} if (n<nMax)printf("%s\t1\t%s\t%d\t%d\t%s\n",$1,$2,a1,a2,$1);}' nMax=$geneBoxMaxLn > tmp/METADATA/$target.r.gene.sponge
    cat tmp/METADATA/gtf.$target.gene2intMap.txt | gawk -F '\t' '{a1=$3;a2=$4;n=a2-a1;if(n<0){n=-n;} if (n<nMax)printf("%s\t1\t%s\t%d\t%d\t%s\n",$1,$2,a1,a2,$1);}' nMax=$geneBoxMaxLn > tmp/METADATA/$target.ns.gene.sponge

    if (-e  tmp/METADATA/gtf.$target.goodProduct.ace.gz) then
      \rm  tmp/METADATA/gtf.$target.goodProduct.ace.gz
    endif
    if (-e tmp/METADATA/gtf.$target.transcripts.ace.gz && ! -e  tmp/METADATA/gtf.$target.goodProduct.aceZZ) then
      gunzip -c tmp/METADATA/gtf.$target.transcripts.ace.gz | gawk '/^Sequence/{if(p1>0)printf("%s\t%d\t%d\n",m,p1,p2);m=$2;dx=1;p1=0;next;}/^Source_exons/{a1=$2;a2=$3;if($4=="CDS"){if(p1==0)p1=a1;p2=a2;}}END{if(p1>0)printf("%s\t%d\t%d\n",m,p1,p2);}' | gawk -F '\t' '{printf ("mRNA %s\nProduct %s %d %d\n\nProduct %s\nBest_product\nGood_product\n\n", $1,$1,$2,$3,$1);}' >  tmp/METADATA/gtf.$target.goodProduct.ace
 \cp tmp/METADATA/gtf.$target.goodProduct.ace  tmp/METADATA/gtf.$target.goodProduct.ace42
    endif

   gunzip -c   tmp/METADATA/gtf.$target.introns.gz | gawk -F '\t' '{if ($9 == "+") printf("%s\t%09d\t%09d\t+\n",$6,$7,$8);}'  | sort -u > tmp/METADATA/$target.f.introns
   gunzip -c   tmp/METADATA/gtf.$target.introns.gz | gawk -F '\t' '{if ($9 == "-") printf("%s\t%09d\t%09d\t-\n",$6,$7,$8);}'  | sort -u > tmp/METADATA/$target.r.introns

   cat tmp/METADATA/gtf.$target.[fr].cds_sponge.gz >  tmp/METADATA/gtf.$target.ns.cds_sponge.gz
   gunzip -c  tmp/METADATA/gtf.$target.ns.cds_sponge.gz > tmp/METADATA/$target.ns.cds.sponge
   gunzip -c  tmp/METADATA/gtf.$target.[fr].sponge.gz > tmp/METADATA/$target.ns.mrna.sponge
   touch tmp/METADATA/$target.gtf.done
   if (-e tmp/METADATA/gtf.$target.mrnaRemap.done) \rm tmp/METADATA/gtf.$target.mrnaRemap.done
   set gtf_active=1

  endif
end

#######################################################################
## create a non redundant genebox sponge file
## Stranded
## construct a projected stranded exon box for av
if (1) then
foreach target ($Etargets)
  # restrict the genebox sponge to the capture
  if (! $?CAPTURES) continue
  foreach capture ($CAPTURES)
    if (! -e TARGET/GENES/$capture.capture.$target.gene_list) continue
    wc  TARGET/GENES/$capture.capture.$target.gene_list
    if (-e tmp/METADATA/$target.ns.any_genebox.$capture.spongeZZ) continue
      foreach fr (f r ns)
        if (! -e  tmp/METADATA/$target.$fr.gene.sponge) continue
        cat  TARGET/GENES/$capture.capture.$target.gene_list ZZZZZ  tmp/METADATA/$target.$fr.gene.sponge | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$1]=1;next;}}{if(gg[$1])print;}' > tmp/METADATA/$target.$fr.gene.$capture.sponge
        # ls -ls tmp/METADATA/$target.$fr.gene.$capture.sponge
        cat tmp/METADATA/$target.$fr.gene.$capture.sponge | gawk -F '\t' '{if($4>$5){u=$4;$4=$5;$5=u;}printf("%s\t%09d\t%09d\n",$3,$4,$5);}' | sort > tmp/METADATA/$target.$fr.genebox.$capture.sponge.1
        cat tmp/METADATA/$target.$fr.genebox.$capture.sponge.1 ZZZZZ  | gawk -F '\t' '{if (chr != $1 || $2 + 0 > a2 + 0){nn++;if(a2+0>1){if(fr == "r"){a0=a1;a1=a2;a2=a0;} printf("Any_%s_%s\t%d\t%s\t%09d\t%09d\t%s\n",gene,fr,nn,chr,a1,a2,gene);}chr=$1;a1=$2;a2=$3;gene=$4;}if(a2<$3)a2=$3;}' t=$target fr=$fr >  tmp/METADATA/$target.$fr.any_genebox.$capture.sponge
      end
    if (-e tmp/METADATA/$target.$fr.genebox.$capture.sponge.1) \rm tmp/METADATA/$target.$fr.genebox.$capture.sponge.1  
    if (-e tmp/METADATA/$target.$fr.exon.sponge.1) \rm tmp/METADATA/$target.$fr.exon.sponge.1
  end
end 
endif
########################################################################
## mrnaRemap, used in GENE_EXPRESSION
## coordinates of all exons, derived from the gtf file

# check that the list of RNAtargets did not change
echo "----  check that the list of RNAtargets did not change"
echo "$allRNAtargets" >  tmp/METADATA/RNAtargets.new
touch tmp/METADATA/RNAtargets
set n=`diff  tmp/METADATA/RNAtargets tmp/METADATA/RNAtargets.new | wc -l`

if ($n > 0) then
  if (-e  tmp/METADATA/mrnaRemap.gz) \rm  tmp/METADATA/mrnaRemap.gz 
  mv   tmp/METADATA/RNAtargets.new  tmp/METADATA/RNAtargets
endif

# construct  tmp/METADATA/mrnaRemap.gz 
if (! -e tmp/METADATA/mrnaRemap.gz && -e TARGET/WIGGLEREMAP/mrnaRemap.av_RefSeq_seqc.txt.gz) then
  cp  TARGET/WIGGLEREMAP/mrnaRemap.av_RefSeq_seqc.txt.gz  tmp/METADATA/mrnaRemap.gz 
endif
# ls -ls  tmp/METADATA/mrnaRemap.gz 

if ($gtf_active == 1 || ! -e tmp/METADATA/mrnaRemap.gz) then

  echo "---- processing mrnaRemap coordinates informations"

  foreach target ($allRNAtargets)
    echo ' ' > tmp/METADATA/$target.MRNA.info.ace
    if (-e tmp/METADATA/gtf.$target.gene2intMap.txt) then 
       # derived from the gtf file
       if (-e tmp/METADATA/gtf.$target.transcripts.ace.gz) then
	 zcat tmp/METADATA/gtf.$target.transcripts.ace.gz | gawk '/^Sequence/{s=$2;next;}/^Title/{printf("Sequence %s\nTitle %s\n\n",s, substr($0,7));}/^LocusLink/{printf("Sequence %s\n",s);print;printf("\n");}' >>  tmp/METADATA/$target.MRNA.info.ace
         echo ' ' | gzip >! tmp/METADATA/toto.gz 
	 zcat tmp/METADATA/gtf.$target.transcripts.ace.gz | gawk '/^$/{out=1;intron=0;next;}/^Intron/{if(out==1){out=0;intron=$2;}next;}/^Gene/{if(intron!=0){printf("Intron %s\nGene %s\n\n",intron, substr($0,5));}}' | gzip >   tmp/METADATA/toto.gz 
         zcat tmp/METADATA/gtf.$target.gene_title.ace.gz ZZZZZ.gz  tmp/METADATA/toto.gz | gawk '/^ZZZZZ/{zz=1;next;}/^$/{intron=0;gene=0;title=0;next;}/^Gene /{gene=$2;gsub(/\"/,"",gene);if(intron==1){print;if(length(g2t[gene])>0)printf("Title %s\n",g2t[gene]);}}/^Title/{title=substr($0,6);if(zz<1 && gene!=0){g2t[gene]=title;next;}}/^Intron/{printf("\n");print;intron=1;gsub(/\"/,"",$2);n1=split($2,aa,"__");n2=split(aa[2],bb,"_");if(n1+n2==4)printf("IntMap %s %s %s\n",aa[1],bb[1],bb[2]);}'  > tmp/METADATA/$target.introns.info.ace
         \rm   tmp/METADATA/toto.gz 
         zcat tmp/METADATA/gtf.$target.transcripts.ace.gz | gawk '/^Sequence/{s=$2;next;}/^Title/{printf("mRNA %s\nTitle %s\n\n",s, substr($0,7));}/^Model_ofZZZ/{printf("mRNA %s\nModel_of %s\n\n",s, substr($0,10));}' >>  tmp/METADATA/$target.MRNA.info.ace
         # cat  tmp/METADATA/$target.GENE.info.ace | gawk '/Gene /{print; printf("Targeted\n\n");}' > tmp/METADATA/_x
         # cat _x >> tmp/METADATA/$target.GENE.info.ace
       endif
       gunzip -c tmp/METADATA/gtf.$target.geneTitle.gz >> tmp/METADATA/$target.GENE.info.ace
       echo ' ' >> tmp/METADATA/$target.GENE.info.ace
       cat tmp/METADATA/gtf.$target.gene2intMap.txt | gawk -F '\t' '/^#/{next;}{if(length($1)>0)printf("Gene \"%s\"\nIntMap \"%s\" %d %d\n\n",$1,$2,$3,$4);}' >>  tmp/METADATA/$target.GENE.info.ace
       cat tmp/METADATA/gtf.$target.mrna2intMap.txt | gawk -F '\t' '/^#/{next;}{if(length($1)>0)printf("mRNA \"%s\"\nIntMap \"%s\" %d %d\n\n",$1,$2,$3,$4);}' >>  tmp/METADATA/$target.MRNA.info.ace
        cat tmp/METADATA/$target.MRNA.info.ace | gawk '/^Sequence /{s=$2;ss[s]=1;next;}/^LocusLink/{LL[s]=$2;next;}/^Title/{tt[s]=substr($0,8);next;}END{for(s in ss){printf("%s\t%s\t%s\n",s,LL[s],tt[s])}}' | sed -e 's/\"//g' > tmp/METADATA/$target.symbol.title.txt
    else
       # backwards compatibility with hand prepared files
       if (-e  TARGET/GENES/$target.gene2intmap.ace) then
          cat  TARGET/GENES/$target.gene2intmap.ace  >>  tmp/METADATA/$target.GENE.info.ace
       endif
       if (-e TARGET/MRNAS/$target.mrna2intmap.ace.gz) then
          gunzip -c TARGET/MRNAS/$target.mrna2intmap.ace.gz  >>  tmp/METADATA/$target.MRNA.info.ace
       endif
    endif
  end

  if (-e tmp/METADATA/mrnaRemap.1) \rm tmp/METADATA/mrnaRemap.1
  touch  tmp/METADATA/mrnaRemap.1 tmp/METADATA/mrnaRemap.2
  if (-e  TARGET/WIGGLEREMAP/mrnaRemap.av_RefSeq_seqc.txt.gz) then
    gunzip -c TARGET/WIGGLEREMAP/mrnaRemap.av_RefSeq_seqc.txt.gz > tmp/METADATA/mrnaRemap.2
  endif

  set mrnaRemap_needed=0
  foreach target ($allRNAtargets)
    if (! -e tmp/METADATA/gtf.$target.mrnaRemap.gz) continue 
 

    source scripts/target2target_class.txt
    mv tmp/METADATA/mrnaRemap.2  tmp/METADATA/mrnaRemap.3
    cat tmp/METADATA/mrnaRemap.3 | gawk -F '\t' "/^$target_class/{next;}{print}" >  tmp/METADATA/mrnaRemap.2
    gunzip -c  tmp/METADATA/gtf.$target.mrnaRemap.gz >> tmp/METADATA/mrnaRemap.1

    touch tmp/METADATA/gtf.$target.mrnaRemap.done
    set mrnaRemap_needed=1
  end

  if ($mrnaRemap_needed == 1) then
    cat tmp/METADATA/mrnaRemap.1 >  tmp/METADATA/mrnaRemap.3
    cat tmp/METADATA/mrnaRemap.2 | sort >>  tmp/METADATA/mrnaRemap.3

    cat  tmp/METADATA/mrnaRemap.3 | gzip > tmp/METADATA/mrnaRemap.gz
  endif

  \rm  tmp/METADATA/mrnaRemap.[123]
endif
echo "----  check done"
# loop on RNAtargets because it is licit to have an RNA_seq project with no previously known annotations
foreach target ($allRNAtargets)
  if ($target == rrna) continue
  if (0 && $Strategy == RNA_seq && ! -e tmp/METADATA/mrnaRemap.gz) then
    echo "---- FATAL ERROR missing file  tmp/METADATA/mrnaRemap.gz"
    goto done
  endif
end

touch  tmp/METADATA/mrnaRemap.gz 
foreach target ($allRNAtargets)
    source scripts/target2target_class.txt
    if (-e tmp/METADATA/mrnaRemap.gz && ! -e   tmp/METADATA/$target.MRNA.splicing.ace) then 
      gunzip -c  tmp/METADATA/mrnaRemap.gz | gawk -F '\t' '/^#/{next;}{if($1 == tc){m=$2;x1=$3;x2=$4;chr=$5;a1=$6;a2=$7;if(x1 == 1) { a0=a1 ;nx=0; printf("\nmRNA %s\n-D Splicing\n",m);}if(a1<a2){u1=a1-a0+1;u2=a2-a0+1;}else{u1=a0-a1+1;u2=a0-a2+1;}nx++;if(nx>1 && u1>v2+20){ln=u1-v2-1;printf ("Splicing %d %d %d %d Intron %d Length %d bp\n",v2+1, u1-1, y2,x1,nx-1,ln);}  }ln=x2-x1+1; if(u1>0) printf ("Splicing %d %d %d %d Exon %d Length %d bp\n",u1, u2, x1,x2,nx,ln);v1=u1;v2=u2;y1=x1;y2=x2;b1=a1;b2=a2;}END{printf("\n");}' tc=$target_class >   tmp/METADATA/$target.MRNA.splicing.ace
    endif
end


foreach target ($allRNAtargets)
  if ($target == introns) continue
  if ($target == rrna) continue
  if ($target == smallRNA) continue

  if (! -e tmp/METADATA/gtf.$target.mrna2intMap.txt && -e TARGET/WIGGLEREMAP/mrnaRemap.$target.txt) then
     printf "ZZZZZ\tZZZZZ\n" > ZZZZZ2
     cat  TARGET/WIGGLEREMAP/mrnaRemap.$target.txt ZZZZZ | gawk -F '\t' '{gsub(/\"/,"",$0);}{m=$2;c[m]=$5;a1=$6;a2=$7;if(a1>a2){a0=a1;a1=a2;a2=a0;isUp[m]=1;}if(m1[m]+0==0){m1[m]=a1;m2[m]=a2;}if(a1<m1[m])m1[m]=a1;if(a2>m2[m])m2[m]=a2;}END{for(m in c){printf("%s\t%s",m,c[m]);if(isUp[m]==1)printf("\t%d\t%d\n",m2[m],m1[m]);else printf("\t%d\t%d\n",m1[m],m2[m]);}}' | sort >  tmp/METADATA/gtf.$target.mrna2intMap.txt 
  endif

  if (-e  tmp/METADATA/$target.mrna_ln_gc_gene_geneid.txt || -e tmp/METADATA/gtf.$target.mrna2intMap.txt) then 
    cat   tmp/METADATA/gtf.$target.mrna2intMap.txt ZZZZZ tmp/METADATA/$target.mrna_ln_gc_gene_geneid.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){m=$1;m2map[m]=$2 ":" $3 "-" $4;next;}printf("%s\t%s",$1,m2map[$1]);for(i=2;i<=NF;i++)printf("\t%s",$i);printf("\n");}' > tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt
# reformat in .ace format tmp/METADATA/av.mrna_map_ln_gc_gene_geneid.txt
  

    cat tmp/METADATA/$target.mrna_map_ln_gc_gene_geneid.txt |   gawk -F '\t' '/^#/{next;}{m=$1;chr=$2;ln=$3;gc=$4;g=$5;gid=$6;if(length(m)<1)next;printf("mRNA \"%s\"\nLength %d\n",m,ln);if(gc+0>0)printf("GC_percent %d\n",gc);k1=split(chr,aa,":");k2=split(aa[2],bb,"-");if(k1==2 && k2==2)printf("IntMap %s %s %s\n",aa[1],bb[1],bb[2]);if(length(g)>1)printf("Gene \"%s\"\n",g);if(length(gid)>1)printf("GeneId \"%s\"\n",gid);printf("\n");}' >  tmp/METADATA/$target.MRNA.ln.ace

  else
    echo "---- FATAL ERROR XXX $allRNAtargets YYY missing file tmp/METADATA/$target.mrna_ln_gc_gene_geneid.txt OR tmp/METADATA/gtf.$target.mrna2intMap.txt"
    goto done
  endif

end

set ok=0
foreach target (av $Etargets)
  if ($ok == 0 && -e tmp/METADATA/$target.mrna_ln_gc_gene_geneid.txt && ! -e tmp/METADATA/$target.selected5kbTranscriptList.txt) then
    # MALAT is extremelly expressed over its firts kb, it is not really an 8kb gene
    cat  tmp/METADATA/$target.mrna_ln_gc_gene_geneid.txt | gawk -F '\t' '/MALAT/{next;}{gene=$4;if(gene == old)next;old=gene;ln=$2;if(ln>=8000)print $1;}' >  tmp/METADATA/$target.selected8kbTranscriptList.txt
    cat  tmp/METADATA/$target.mrna_ln_gc_gene_geneid.txt | gawk -F '\t' '/MALAT/{next;}{gene=$4;if(gene == old)next;old=gene;ln=$2;if(ln>=5000)print $1;}' >  tmp/METADATA/$target.selected5kbTranscriptList.txt
    set ok=1
  endif
end

# tmp/METADATA/mrnaRemap.gz  TARGET/WIGGLEREMAP/mrnaRemap.av_RefSeq_seqc.txt.gz

######################################################
######################################################
# sponge files


echo "---- sponge files"
echo $allRNAtargets

# WIGGLEREMAP, obsolete, backward compatibility if we do not have the gtf file
if (-d TARGET/WIGGLEREMAP) then

  foreach target ($allRNAtargets)
    foreach fr (f r ns)
      if (-e TARGET/WIGGLEREMAP/$target.GENE.$fr.sponge && ! -e TARGET/WIGGLEREMAP/$target.$fr.gene.sponge) then
        mv  TARGET/WIGGLEREMAP/$target.GENE.$fr.sponge  TARGET/WIGGLEREMAP/$target.$fr.gene.sponge
      endif
      if (-e TARGET/WIGGLEREMAP/$target'_cds'.$fr.sponge && ! -e TARGET/WIGGLEREMAP/$target.$fr.cds.sponge) then
        mv  TARGET/WIGGLEREMAP/$target'_cds'.$fr.sponge  TARGET/WIGGLEREMAP/$target.$fr.cds.sponge
      endif
      if (-e TARGET/WIGGLEREMAP/$target.$fr.sponge && ! -e TARGET/WIGGLEREMAP/$target.$fr.mrna.sponge) then
        mv  TARGET/WIGGLEREMAP/$target.$fr.sponge  TARGET/WIGGLEREMAP/$target.$fr.mrna.sponge
      endif
    end

   foreach fr (f r ns)
     if (-e TARGET/WIGGLEREMAP/$target.$fr.mrna.sponge && ! -e TARGET/WIGGLEREMAP/$target.$fr.gene.sponge) then
        cat TARGET/WIGGLEREMAP/$target.$fr.mrna.sponge | gawk -F '\t' '{m=$1;n=$2;c=$3;a1=$4;a2=$5;g=$6;if(ga1[g]+0<1){ga1[g]=a1;gaa2[g]=a2;gc[g]=c;}if(a1<ga1[g])ga1[g]=a1;if(a2>ga2[g])ga2[g]=a2;}END{for(g in gc)printf("%s\t1\t%s\t%d\t%d\t%s\n",g,gc[a],ga1[g],ga2[g],g);}' | sort -V > TARGET/WIGGLEREMAP/$target.$fr.gene.sponge
     endif
   end
   foreach fr (f r)
     if (-e TARGET/WIGGLEREMAP/$target.$fr.mrna.sponge && -e TARGET/WIGGLEREMAP/$target.ns.cds.sponge && ! -e TARGET/WIGGLEREMAP/$target.$fr.cds.sponge) then
        cat TARGET/WIGGLEREMAP/$target.$fr.mrna.sponge ZZZZZ  TARGET/WIGGLEREMAP/$target.ns.cds.sponge | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz+0<1){gg[$1 "_CDS"]=1;next;}}{if(gg[$1]==1)print;}'  | sort -V > TARGET/WIGGLEREMAP/$target.$fr.cds.sponge
     endif
   end

    if (-e TARGET/WIGGLEREMAP/mrnaRemap.$target.txt && ! -e TARGET/WIGGLEREMAP/$target.ns.mrna.sponge) then
      cat TARGET/WIGGLEREMAP/mrnaRemap.$target.txt | gawk -F '\t' '{m=$2;g=m;if($8)g=$8;if(m!=oldm)n=0;n++;oldm=m;a1=$6;a2=$7; if(a1<a2)printf("%s\t%d\t%s\t%d\t%d\t%s\n",m,n,$5,a1,a2,g);}' | sort -k 1,1 -k 2,2n  >  TARGET/WIGGLEREMAP/$target.f.mrna.sponge
      cat TARGET/WIGGLEREMAP/mrnaRemap.$target.txt | gawk -F '\t' '{m=$2;g=m;if($8)g=$8;if(m!=oldm)n=0;n++;oldm=m;a1=$6;a2=$7; if(a1>a2)printf("%s\t%d\t%s\t%d\t%d\t%s\n",m,n,$5,a1,a2,g);}' | sort -k 1,1 -k 2,2n  >  TARGET/WIGGLEREMAP/$target.r.mrna.sponge
      cat TARGET/WIGGLEREMAP/mrnaRemap.$target.txt | gawk -F '\t' '{m=$2;g=m;if($8)g=$8;if(m!=oldm)n=0;n++;oldm=m;a1=$6;a2=$7; printf("%s\t%d\t%s\t%d\t%d\t%s\n",m,n,$5,a1,a2,g);}' | sort -k 1,1 -k 2,2n  >  TARGET/WIGGLEREMAP/$target.ns.mrna.sponge
    endif

    foreach fr (f r ns)
      if (-e TARGET/WIGGLEREMAP/$target.$fr.mrna.sponge && ! -e tmp/METADATA/$target.$fr.mrna.sponge) then
        cp TARGET/WIGGLEREMAP/$target.$fr.mrna.sponge  tmp/METADATA
      endif
      if (-e TARGET/WIGGLEREMAP/$target.$fr.gene.sponge && ! -e tmp/METADATA/$target.$fr.gene.sponge) then
        cp TARGET/WIGGLEREMAP/$target.$fr.gene.sponge  tmp/METADATA
      endif
      if (-e TARGET/WIGGLEREMAP/$target.$fr.cds.sponge && ! -e tmp/METADATA/$target.$fr.cds.sponge) then
        cp TARGET/WIGGLEREMAP/$target.$fr.cds.sponge  tmp/METADATA
      endif
    end
  end
endif


# hack if no gtf
foreach target ($allRNAtargets)
  if (-e tmp/METADATA/$target.ns.sponge && ! -e tmp/METADATA/$target.ns.mrna.sponge) then
    mv  tmp/METADATA/$target.ns.sponge tmp/METADATA/$target.ns.mrna.sponge 
  endif
  if (-e tmp/METADATA/$target'_cds.ns.sponge' && ! -e  tmp/METADATA/$target.ns.cds.sponge) then
    mv tmp/METADATA/$target'_cds'.ns.sponge  tmp/METADATA/$target.ns.cds.sponge
  endif
end

######################################################
## known introns

echo "---- known introns"
foreach target ($allRNAtargets)

  if ($target == rrna) continue
  if (0 && -e TARGET/MRNAS/introns_$target.txt && ! -e tmp/METADATA/$target.f.manual_introns) then
    cat TARGET/MRNAS/introns_$target.txt  | gawk -F '\t' '{if ($2<$3)  printf("%s\t%09d\t%09d\t+\n",$1,$2, $3);}' >  tmp/METADATA/$target.f.manual_introns
    cat TARGET/MRNAS/introns_$target.txt  | gawk -F '\t' '{if ($2>$3)  printf("%s\t%09d\t%09d\t-\n",$1,$2, $3);}' >  tmp/METADATA/$target.r.manual_introns
    \cp tmp/METADATA/$target.f.manual_introns tmp/METADATA/$target.f.introns
    \cp tmp/METADATA/$target.r.manual_introns tmp/METADATA/$target.r.introns
  endif

  foreach fr (f r)
    if (0 && -e tmp/METADATA/$target.$fr.sponge && ! -e tmp/METADATA/$target.$fr.introns) then
      cat  tmp/METADATA/$target.$fr.sponge | gawk -F '\t' '{if (old==$1 && chrom==$3 && num == $2 - 1){a1=$4;dx=1;ss="+";if(b2>a1){dx=-1;ss="-";}printf("%s\t%09d\t%09d\t%s\n",chrom,b2+dx,a1-dx,ss);} old=$1;num=$2;chrom=$3;b1=$4;b2=$5;}' | sort -u > tmp/METADATA/$target.$fr.introns
    endif
  end
  foreach fr (f r)
    if (-e  tmp/METADATA/gtf.$target.introns.gz && ! -e tmp/METADATA/$target.$fr.introns) then
      zcat  tmp/METADATA/gtf.$target.introns.gz | gawk -F '\t' '{chrom=$6;a1=$7;a2=$8;if (a1<a2){ss="+";if ($fr=="r")next;}else{a0=a1;a1=a2;a2=a0;ss="-";if ($fr=="f")next;}printf("%s\t%09d\t%09d\t%s\n",chrom,b2+dx,a1-dx,ss);}' fr=$fr  | sort -u > tmp/METADATA/$target.$fr.introns
    endif
  end
  if (! -e tmp/METADATA/$target.fr.introns.gz) then
      cat  tmp/METADATA/$target.[fr].introns | gzip > tmp/METADATA/$target.fr.introns.gz
  endif

end 


######################################################
######################################################
# supplementary information not from .fasta and not from .gtf

if (-e TARGET/GENES/av.split_mrnas.gz && ! -e tmp/METADATA/av.split_mrnas.gene2length.aceZZZ) then
  cp  TARGET/GENES/av.split_mrnas.gz  tmp/METADATA/av.split_mrnas.gz
# find their length
  cat TARGET/GENES/av.genes2length.ace | gawk '/^Gene/{gsub(/\"/,"",$2);g=$2;}/^Length/{x=$2+0;if(x>0)printf("%s\t%d\n",g,x);}' > _tln
  echo ZZZZZ >> _tln
  gunzip -c TARGET/Targets/$species.av.TM.txt.gz ZZZZZ.gz tmp/METADATA/av.split_mrnas.gz >> _tln
  cat _tln | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){g2l[$1]=$2;next;}}{if(zz==1){split($1,aa,"|");m2l[aa[1]]=$2;next;}}{g=$7;x=0;if($2=="*")x=g2l[$1];if(x==0){a=$3;b=$4;x=b-a;}if(x==0)x=m2l[$2]+0;if(x>0 && x>ln[g]+0)ln[g]=x;}END{for(g in ln)printf("Gene %s\nLength %d\n\n",g,ln[g]);}' > tmp/METADATA/av.split_mrnas.gene2length.ace
  \rm _tln
endif

######################################################
######################################################
#### analyse the telomeric motifs in the targets

scripts/a0B.Bloom.tcsh METADATA

######################################################
#### construct the metadata from an AceView database using table make

if (0 && ! -e TARGET/GENES/AceView_table.metadata.txt) then
  bin/tacembly 37a_1 << EOF
    table -o $MagicRootDir/TARGET/GENES/AceView_table.metadata.txt -f $MagicRootDir/metaData/tables/AceView_table.metadata.def
    quit
EOF
endif
foreach target ($allRNAtargets)
if ($target == av && -e TARGET/GENES/AceView_table.metadata.txt && ! -e TARGET/GENES/av.metadata.txt) then

  echo "---- processing AceView metadata"

  # set toto=RESULTS/AceView_2010.metadata.txt
  set toto=TARGET/GENES/av.metadata.txt
  echo -n "# " > $toto
  date >> $toto
  echo "# Metadada for the 72695 AceView-2010." >> $toto
  printf "# Gene\tChromosome:from base - to base (genome GRCh37)\tNCBI GeneId (associated to RefSeq)\tGene type\tGene descriptor\tTranscripts\tUniquely mapped Agilent probes\n" >> $toto

  cat  TARGET/GENES/AceView_table.metadata.txt | gawk -F '\t' '/^#/{next;}{gsub(/\"/,"",$0);gsub(/\\; /,";",$0);}{g=$1;map[g]=$2 ":" $3 "-" $4;if($5 != "NULL")gid[g]=$5;if($6 != "NULL") coding[g]=2;if($7 != "NULL") coding[g]=1;if($8!="NULL")conserved[g]=1;if($9 != "NULL")spl_non_coding[g]=1;if($10 != "NULL")spliced[g]=2;if($11 != "NULL")spliced[g]=1;if($12 != "NULL")as[g]=1;if($13 != "NULL")title[g]=$13;if($14 != "NULL")mrna[g]=mrna[g] ";" $14;if($15 != "NULL")microarray[g]= $15;}END{for(g in map){if(length(map[g])<5)continue;if(coding[g]<1)nam="Non_coding";if(coding[g]==1)nam="Marginally_coding";if(coding[g]==2)nam="Coding";if(coding[g]>0){if(conserved[g]==1) nam=nam "_ancient" ; else nam=nam "_recent";}if(spliced[g]==2)nam = nam "_spliced" ;if(spliced[g]==1)nam = nam "_single_exon";if(spliced[g]==2){if(as[g]==1) nam=nam "_antisense" ; else nam = nam "_free";}printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",g,map[g],gid[g],nam,title[g],substr(mrna[g],2),microarray[g]);}}' | sort >> $toto

  #### special hack for av 2010
  #    cat  tmp/METADATA/$target.mrna_ln_gc_gene_geneid.txt ZZZZZ TARGET/GENES/av.metadata.txt | gawk -F '\t' '/^#/{printf("#\t");print;next;}/^ZZZZZ/{zz++;next;}{if(zz<1){g=$4;ok[g]=1;next;}}{g=$1;if(ok[g]==1)printf("kept\t");else printf("rejected\t") ;print}'
  ####
endif

if (! -e tmp/METADATA/$target.metadata.txt) then
  # to insure consistency give priority to the gene list in the fasta file
  cat  tmp/METADATA/$target.mrna_ln_gc_gene_geneid.txt > tmp/METADATA/toto
  echo ZZZZZ >>  tmp/METADATA/toto
  if (-e TARGET/GENES/av.metadata.txt) cat  TARGET/GENES/av.metadata.txt >>  tmp/METADATA/toto
  cat  tmp/METADATA/toto | gawk -F '\t' '/^#/{printf("#\t");print;next;}/^ZZZZZ/{zz++;next;}{if(zz<1){m=$1;g=$4;ok[g]=1;g2m[g]=g2m[g]";"m;next;}}{g=$1;if(ok[g]){printf("%s",g);for(i=2;i<=5;i++)printf("\t%s",$i);printf("\t%s\t%s\n",substr(g2m[g],2),$7);}}' > tmp/METADATA/$target.metadata.txt
  \rm  tmp/METADATA/toto

endif

    if (-e tmp/METADATA/$target.metadata.txt) then
      cat  tmp/METADATA/$target.metadata.txt | gawk -F '\t' '/^#/{next;}{printf("Gene \"%s\"\n", $1) ; if($3)printf("GeneId \"%s\"\n",$3);if($4)printf("Gene_type \"%s\"\n",$4);split($2,aa,":");if($5)printf("Title \"%s\"\n",$5);split($2,aa,":");n=split(aa[2],bb,"-");if(n==2)printf("IntMap \"%s\" %d %d\n",aa[1],bb[1],bb[2]); printf("\n");}'  >>  tmp/METADATA/$target.GENE.info.ace
     cat  tmp/METADATA/$target.metadata.txt | gawk -F '\t' '/^#/{next;}{n=split($6,aa,";");for(i=1;i<=n;i++){if($5)printf("mRNA \"%s\"\nTitle \"%s\"\n\n",aa[i],$5);}}'  >>  tmp/METADATA/$target.MRNA.info.ace
    endif

    if (-e TARGET/GENES/$target.gene2title.ace)  then
      cat TARGET/GENES/$target.gene2title.ace  >>  tmp/METADATA/$target.GENE.info.ace
      cat TARGET/GENES/$target.gene2title.ace  >>  tmp/METADATA/$target.MRNA.info.ace
    endif

  if (-e TARGET/GENES/$target.SEQC_sigma.ace) then
    cat TARGET/GENES/$target.SEQC_sigma.ace >>  tmp/METADATA/$target.GENE.info.ace
  endif 
echo "---- processing AceView metadata done"


end

########################################################

foreach target ($allRNAtargets) 
  if (-e tmp/METADATA/$target.Gene_type.done) continue
  touch tmp/METADATA/$target.Gene_type.done

 if (-e  TARGET/GENES/$target.gene2title.ace) then
     cat TARGET/GENES/$target.gene2title.ace >>  tmp/METADATA/$target.GENE.info.ace
  endif
  if (-e  TARGET/GENES/$target.gene2WbId.ace) then
     cat TARGET/GENES/$target.gene2WbId.ace >>  tmp/METADATA/$target.GENE.info.ace
  endif
  if (-e  TARGET/GENES/$target.introns.type.ace) then
     cat TARGET/GENES/$target.introns.type.ace >>  tmp/METADATA/$target.GENE.info.ace
  endif

  if (-e  tmp/METADATA/$species.$target.TM.txt.gz) then
    gunzip -c  tmp/METADATA/$species.$target.TM.txt.gz | gawk '{n=split($1,aa,"|");if(n>=5 && (aa[2]=="Gene" || aa[2]=="GENE") && aa[4]=="GeneId")printf("Gene \"%s\"\nGeneId \"%s\"\n\n",aa[3],aa[5]);}' >>  tmp/METADATA/$target.GENE.info.ace
  endif

  if ($target == RefSeq && -e  tmp/METADATA/$species.$target.TM.txt.gz) then
    gunzip -c  tmp/METADATA/$species.$target.TM.txt.gz | gawk '{n=split($1,aa,"|");if(n>=3 && (aa[2]=="Gene" || aa[2]=="GENE") && substr($1,2,2) == "M_")printf("Gene \"%s\"\nPastille_coding\n\n",aa[3]);}' >>  tmp/METADATA/$target.GENE.info.ace
  endif

end

######################################################
######################################################
# Manipulate the known transcript ends

    if (-e TARGET/GENES/Ecocyc.ace  &&  -e tmp/METADATA/Eco.knownEnds.done) then
      foreach chrom ($chromSetAll)
          set toto=tmp/METADATA/Eco.knownStarts.$chrom
          cat TARGET/GENES/Ecocyc.ace | gawk '/^IntMap/{if($2 == chrom){printf("%d\t%d\n",$3,$4);}}' chrom=$chrom > $toto.a
          cat $toto.a |  gawk '{if($1<$2)print $1; }' | sort -u | sort -k 1n > $toto.f
          cat $toto.a |  gawk '{if($1>$2)print $1; }'  | sort -u | sort -k 1n > $toto.r
          set toto1=tmp/METADATA/Eco.knownEnds.$chrom
          cat $toto.a |  gawk '{if($1<$2)print $2; }' | sort -u | sort -k 1n > $toto1.f
          cat $toto.a |  gawk '{if($1>$2)print $2;}'  | sort -u | sort -k 1n > $toto1.r
	  \rm $toto.a
      end
      touch tmp/METADATA/Eco.knownEnds.done
    endif

######################################################
######################################################
# Manipulate the masking files

foreach target ($allRNAtargets $DNAtargets)
  if (-e TARGET/Targets/$species.genome.mask.txt && ! -e  TARGET/Targets/$species.$target.mask.txt) then
    if (-e TARGET/WIGGLEREMAP/mrnaRemap.$target.txt) then

      echo "---- masking the genes falling inside the maked genomic regions : $target"

      source scripts/target2target_class.txt
      cat TARGET/Targets/$species.genome.mask.txt | gawk '/^\//{next}/^#/{next}{gsub(/,/,"",$0);if(0+$2>1)printf("%s\t%d\t%d\n",$1,0+$2,0+$3);}' | sort -k 1,1 -k 2,2n > tmp/METADATA/toto.mask 
      cat ZZZZZ >> tmp/METADATA/toto.mask 
      gunzip -c tmp/METADATA/mrnaRemap.gz | gawk '{if($1==tcl)print;}' tcl=$target_class >> tmp/METADATA/toto.mask 
      cat tmp/METADATA/toto.mask  | gawk '/^#/{next}/^\//{next;}/ZZZZZ/{zz=1;next;}{if(zz<1){nm++;c1=$2;c2=$3;if(c1>c2){c0=c1;c1=c2;c2=c0;}mc[nm]=$1;mc1[nm]=c1;mc2[nm]=c2;if (chrom1[$1]<1)chrom1[$1]=nm;chrom2[$1]=nm ;next;}gene=$2;x1=$3;x2=$4;c=$5; gsub(/^c_/,"",c);if(chrom1[c]<1)next;a1=$6;a2=$7;if(a1>a2){b1=a2;b2=a1;isDown=0;}else{b1=a1;b2=a2;isDown=1;}for(i=chrom1[c];i<=chrom2[c];i++){c1=mc1[i];c2=mc2[i];if(c1>b1)u1=c1;else u1=b1;if(c2<b2)u2=c2;else u2=b2;if(0)print "i=" i "u1=" u1 "u2=" u2 ;if(u1>u2)continue;if(isDown==1){w1=x1+u1-b1;w2=x1+u2-b1;}else{w1=x1+b2-u2;w2=x1+b2-u1;}printf("%s\t%d\t%d\n",gene,w1,w2);}}'  >  TARGET/Targets/$species.$target.mask.txt
      \rm tmp/METADATA/toto.mask 
    
    endif
  endif
end

######################################################
######################################################
# Document the targets in MetaDB

echo "--- Document the targets in MetaDB"
if (! -e  MetaDB/$MAGIC/target_counts.ace) then

  set toto = MetaDB/$MAGIC/target_counts.ace
  echo "// $MAGIC" > $toto

  foreach target ($allRNAtargets)
    # the .gz file is used in ii2b
    if (! -e tmp/METADATA/$target.fr.introns.gz) then
      cat  tmp/METADATA/$target.[fr].introns | gzip > tmp/METADATA/$target.fr.introns.gz
    endif

    if (! -e tmp/METADATA/$target.counts.ace) then
      cat tmp/METADATA/$target.[fr].introns | gawk '{n++;}END{printf("Target %s\nNumber_of_introns %d\n\n",t,n);}' t=$target >  tmp/METADATA/$target.counts.ace
    endif
    cat  tmp/METADATA/$target.counts.ace >> $toto
  end

  echo "pparse $toto" | bin/tacembly MetaDB -noprompt
endif

######################################################
# captures

if ($?CAPTURES) then
  foreach target ($Etargets)
    if (-e tmp/METADATA/$MAGIC.$target.captured_genes.txt) \rm tmp/METADATA/$MAGIC.$target.captured_genes.txt
    foreach capture ($CAPTURES)

      if (-e TARGET/GENES/$capture.capture.$target.gene_list) then
        cat  TARGET/GENES/$capture.capture.$target.gene_list | gawk '/^#/{next;}{if(length($1)>0)printf("%s\t%s\t%s\n",$1,cap,"Capture");}' cap=$capture >> tmp/METADATA/$MAGIC.$target.captured_genes.txt
      endif
      if (-e TARGET/GENES/$capture.captureTouch.$target.gene_list) then
        cat  TARGET/GENES/$capture.captureTouch.$target.gene_list | gawk '/^#/{next;}{if(length($1)>0)printf("%s\t%s\t%s\n",$1,cap,"Capture_touch");}' cap=$capture >> tmp/METADATA/$MAGIC.$target.captured_genes.txt
      endif
    end
    cat   tmp/METADATA/$MAGIC.$target.captured_genes.txt | sort -V | gawk -F '\t' '{if($1 != old)printf ("\nGene %s\n",$1);old=$1;printf("%s %s\n",$3, $2);}END{printf("\n");}' > tmp/METADATA/$MAGIC.$target.captured_genes.ace
    if (-e tmp/METADATA/$MAGIC.$target.captured_genes.txt) \rm tmp/METADATA/$MAGIC.$target.captured_genes.txt
    if (-e tmp/METADATA/$MAGIC.$target.captured_genes.ace) then
      cat  tmp/METADATA/$MAGIC.$target.captured_genes.ace | gawk '/^Gene/{g=$2;gsub(/\"/,"",g);split(g,aa,"(");g=aa[1];}/touch/{next;}/^Capture/{c=$2;if(length(c)==2)print("%s\t%s\n",g,c);}' | sort -u | gawk '{g=$1;c[g]=c[g] $2;}END{for (g in c)printf("%s\t%s\n",g,c[g]);}' > tmp/METADATA/$MAGIC.$target.captured_genes.g2c
    endif
  end
endif

goto ok

######################################################
# Load GeneIndexDB

bin/tacembly GeneIndexDB << EOF
   parse TARGET/MRNAS/$species.introns.ace.gz
   parse TARGET/MRNAS/$species.introns_RNA_seq.ace.gz
   query find transcribed_gene gene AND Intron
   bql -a -o TARGET/GENES/tg2g2ii.txt select tg,g,ii from tg in @, g in tg->gene, ii in tg->intron
   save
   quit

if (! -e  TARGET/GENES/gene2intron.ace) then
  
  cat  TARGET/GENES/tg2g2ii.txt | cut -f 2,3 | sort -u | gawk  -F '\t' '/^#/{next;}{g=$1;if(g != old)printf("\n\nGene %s", g);old=g;printf("\nIntron %s", $2);}END{printf("\n\n");}' >  TARGET/GENES/gene2intron.ace
endif


ok:
  echo "=== All Target metadata are available"
  touch tmp/METADATA/metadata.ok
  exit 0
done:
  exit 1
