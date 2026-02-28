#!bin/tcsh -f
set chrom=$1
setenv ici `pwd`

if (! -d tmp/X.$MAGIC/$chrom) mkdir  tmp/X.$MAGIC/$chrom
if (! -d tmp/X.$MAGIC/$chrom/database) then
   pushd tmp/X.$MAGIC/$chrom
      mkdir database
      ln -s ../../../metaData/tables
      if (-e TABIX) \rm TABIX
      ln -s ../../TABIX
      mkdir wspec
        \cp  $ici/metaData/wspec.aceview_web_site/*.wrm wspec
        pushd wspec
          \rm models.wrm
          ln -s ../../../../metaData/wspec.aceview_web_site/models.wrm
        popd
        set mynam=`whoami`
        set n=`cat wspec/passwd.wrm | gawk '{if($1==mynam)n++}END{printf("%d", 0+n)}' mynam=$mynam`
     if ($n == 0) then
           whoami >> wspec/passwd.wrm
        endif
      echo y |  $ici/bin/tacembly .
   popd
endif

echo "Clone R\nMain_clone\nMainTitle $MAGIC\nSpecies $species\n" > tmp/X.$MAGIC/$chrom/Main_clone.ace

set target=`echo $Etargets | gawk '{print $1}'`

if (! -e tmp/X.$MAGIC/$chrom/f3.genome.done2) then

    if ($target == EBI && $species == Dmelanogaster && -e  TARGET/GTF/$species.$target.gtf.gz && -d TARGET/GENES && ! -e TARGET/GENES/f3.gtf.$target.gene2FB_symbol.ace) then
       gunzip -c TARGET/GTF/$species.$target.gtf.gz | gawk -F '\t' '{if($3 == "gene"){z=$9;split(z,aa,";");split(aa[1],bb," ");split(aa[3],cc," ");if(bb[1]=="gene_id" && cc[1]=="gene_name")printf("Sequence X__%s\nFB_symbol %s\n\n",bb[2],cc[2]);}}' > TARGET/GENES/f3.gtf.$target.gene2FB_symbol.ace
    endif
    if ($target == FlyBase && -e  TARGET/GTF/$species.$target.gtf.gz && -d TARGET/GENES && ! -e TARGET/GENES/f3.gtf.$target.gene2FB_symbol.ace) then
       gunzip -c TARGET/GTF/$species.$target.gtf.gz | gawk -F '\t' '{if($3 == "gene"){z=$9;gsub(/\"/,"",z);split(z,aa,";");split(aa[1],bb," ");split(aa[2],cc," ");if(bb[1]=="gene_id" && cc[1]=="gene_symbol")printf("Sequence X__%s\nFB_symbol %s\n\n",bb[2],cc[2]);}}' > TARGET/GENES/f3.gtf.$target.gene2FB_symbol.ace
    endif

printf "-R Sequence c_$chrom $chrom\n\n" >  tmp/X.$MAGIC/$chrom/f3.rename_chrom.ace

    bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
      read-models
      find est XH*
      kill
      pparse TARGET/CHROMS/$species.chrom_$chrom.fasta.gz
      pparse  tmp/f1.strategy.ace
      // pparse  tmp/pA/$chrom/$MAGIC.pA.$chrom.feature.ace 
      pparse  tmp/TABIX/$MAGIC/$chrom.tabix.ace
      pparse MetaDB/$MAGIC/runs.ace
      parse tmp/METADATA/gtf.$target.f.intron.ace
      parse tmp/METADATA/gtf.$target.r.intron.ace
      parse tmp/METADATA/gtf.RefSeq.transcripts.ace.gz
      parse tmp/METADATA/$MAGIC.$target.captured_genes.ace
      parse tmp/METADATA/av.GENE.info.ace
      query find sequence locuslink
      spush
      query intmap == $chrom
      sminus
      spop
      kill
      query find sequence $chrom
      kstore ss
      acem
        make_subseq -dna c t 1000000 10000 // this breaks my 6 contigs into tiles
        quit                    // oct 15 2001, i changed from 400 kb to 600 kb
      kget ss
      Follow DNA
      kill 
      undo
      edit -D DNA
      parse tmp/X.$MAGIC/$chrom/f3.rename_chrom.ace
      query find sequence genomic
      follow source
      edit  YBR_contig  // NT_sequences 
 
      save
      quit
EOF

# in Xcds the score is independant of the coverage, it just reflects the length of the CDS, so that all long CDS are shown with all their exons
foreach target ($Etargets)
  if ($target == magic) continue
  set ffcds=tmp/METADATA/gtf.$target.ns.cds.spongeZZZZZ
  if (! -e $ffcds) then
    set ffcds=tmp/METADATA/$target.ns.cds.sponge
  endif
  echo AAA $ffcds
  if (-e $ffcds) then
    cat $ffcds ZZZZZ $ffcds | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if (chrom != $3)next;}{a1=$4;a2=$5;ln=a2-a1;if(ln<0)ln = -ln;ln++;}{if(zz<1){t2ln[$1]+=ln;next;}}{nnn=t2ln[$1];if(nnn >= 270)printf("Xcds_%s__%d_%d\t1\t%d\t%s\t%d\t%d\t%d\n",chrom,a1,a2,ln,chrom,a1,a2,nnn);}' chrom=$chrom | sort -u | sort -V >  tmp/X.$MAGIC/$chrom/f3.$target.cds.txt2
      bin/dna2dna -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow  tmp/X.$MAGIC/$chrom/f3.$target.cds.txt2 >  tmp/X.$MAGIC/$chrom/f3.$target.cds.fasta
      cat  tmp/X.$MAGIC/$chrom/f3.$target.cds.txt2 | gawk -F '\t' '{z=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6;if($7>nn[z]+0)nn[z]=$7+0;}END{for(z in nn)printf("%s\t%d\n",z,nn[z]);}' > tmp/X.$MAGIC/$chrom/f3.$target.cds.txt3
      cat  tmp/X.$MAGIC/$chrom/f3.$target.cds.txt3 | gawk -F '\t' '{printf("Sequence %s\ncdna_clone %s\nForward\nColour LIGHTVIOLET\nIntMap %s %d %d\nIs_read\nComposite %d\n\n",$1,$1,$4,$5,$6,$7);}' >  tmp/X.$MAGIC/$chrom/f3.$target.cds.ace
      bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
        pparse  tmp/X.$MAGIC/$chrom/f3.$target.cds.ace
        pparse  tmp/X.$MAGIC/$chrom/f3.$target.cds.fasta
        save
        quit
EOF
  endif
end

    bin/gene2chrom2 -any -gs -i tmp/X.$MAGIC/$chrom  >! tmp/X.$MAGIC/$chrom/f3.g2c.gsi.ace
    bin/gene2chrom2 -any -pg -i tmp/X.$MAGIC/$chrom  >! tmp/X.$MAGIC/$chrom/f3.g2c.pgi.ace

    bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
      pparse  tmp/X.$MAGIC/$chrom/f3.g2c.gsi.ace
      pparse  tmp/X.$MAGIC/$chrom/f3.g2c.pgi.ace
      save
      quit
EOF

    touch  tmp/X.$MAGIC/$chrom/f3.genome.done
endif

if (! -e   tmp/X.$MAGIC/$chrom/f3.cosmid2map.$chrom.txt) then
    echo "get cosmid2map $chrom"
    bin/tacembly tmp/X.$MAGIC/$chrom <<EOF > /dev/null
      query find sequence genomic  &&  NOT Is_gene_tile 
      bql -a -o tmp/X.$MAGIC/$chrom/f3.cosmid2map.$chrom.txt  select s,m,a1,a2 from s in @, m in s->intmap, a1 in m[1], a2 in a1[1]
      quit
EOF
endif

echo "XXXXXXX f3.parse"
if (! -e tmp/X.$MAGIC/$chrom/f3.parse.done) then

  touch tmp/X.$MAGIC/$chrom/f3.genes.ace
  if (1 && -e TARGET/GENES/av.gene.ace) then
    # cree des gene box
    cat TARGET/GENES/av.gene.ace | gawk '/^IntMap/{m=$2;gsub(/\"/,"",m);if(m!=chrom)next;{ok=1;gg = gg "IntMap " m " " $3 " "  $4 "\n"; cc=cc "Genes " g " " $3 " " $4 "\n";}next;}/^av/{next;}/^Transcribed_gene/{next;}/^Gene /{if(ok==1)print gg "\n";ok=0;gg=$0 "\n"; g = $2;next;}{gg = gg $0 "\n" ;next;}END{if(ok)print gg "\n";print "Sequence " chrom "\n" cc "\n"}' chrom=$chrom > tmp/X.$MAGIC/$chrom/f3.genes.ace
  else 
    # cree des gene models de type le premier de Etargets (may be magic, on preferererait av ou RefSeq)
    set target=`echo $Etargets | gawk '{print $1}'`
    set target=av
    if (1 && -e tmp/METADATA/gtf.$target.transcripts.ace.gz) then
      gunzip -c tmp/METADATA/gtf.$target.transcripts.ace.gz > tmp/X.$MAGIC/$chrom/f3.genes.ace
    endif
  endif

  set ch=""
  if ($species == worm) set ch=CHROMOSOME_

  pushd tmp/X.$MAGIC/$chrom
    $ici/scripts/rrf
    $ici/bin/tacembly  $ici/tmp/X.$MAGIC/$chrom < _r
  popd



  bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
    read-models
    parse metaData/methods.ace
    parse tmp/TABIX/$MAGIC/$chrom.tabix.ace
    pparse MetaDB/$MAGIC/runs.ace
    pparse MetaDB/$MAGIC/samples.ace
    pparse MetaDB/$MAGIC/ali.ace
    // pparse my.chrom.alias.ace
    pparse tmp/X.$MAGIC/$chrom/f1.XI.ace.gz
    pparse tmp/X.$MAGIC/$chrom/f1.XI.fasta.gz
    pparse tmp/X.$MAGIC/$chrom/f2.XY.ace
    pparse tmp/X.$MAGIC/$chrom/f2.XY.fasta
    pparse tmp/X.$MAGIC/$chrom/f2.XW.ace
    pparse tmp/X.$MAGIC/$chrom/f2.XW.fasta
    pparse tmp/X.$MAGIC/$chrom/f2.XA.ace
    pparse tmp/X.$MAGIC/$chrom/f2.XA.fasta
    pparse tmp/X.$MAGIC/$chrom/f2.XSL.ace
    pparse tmp/X.$MAGIC/$chrom/f2.XSL.fasta
    pparse  tmp/X.$MAGIC/$chrom/f3.genes.ace
    save
    query find intron
    list -a -f tmp/X.$MAGIC/$chrom/f3.intron.list
    quit
EOF

if (-d tmp/INTRON_DB/$chrom) then
  bin/tacembly tmp/INTRON_DB/$chrom << EOF
    key tmp/X.$MAGIC/$chrom/f3.intron.list
    show -a -f tmp/X.$MAGIC/$chrom/f3.introns.preace
    quit
EOF

  if (-e tmp/X.$MAGIC/$chrom/f3.introns.preace) then
    cat tmp/X.$MAGIC/$chrom/f3.introns.preace  | gawk '/^Magic/{next;}{print;}' >  tmp/X.$MAGIC/$chrom/f3.introns.ace
    \rm   tmp/X.$MAGIC/$chrom/f3.introns.preace 
  endif
endif

hello:

set ggNS="toto"
if (-e MetaDB/$MAGIC/GroupW_new_exonList) then
  set ggNS1=` cat  MetaDB/$MAGIC/GroupW_new_exonList | gawk '{printf("%s",$1);exit;}'`  
  if ($ggNS1 != "") set ggNS=$ggNS1
endif
set ggS="toto"
ls -ls  MetaDB/$MAGIC/GroupW_strandedList
if (-e MetaDB/$MAGIC/GroupW_strandedList) then
  set ggS1=` cat  MetaDB/$MAGIC/GroupW_strandedList | gawk '{printf("%s",$1);exit;}'`
  if ($ggS1 != "") set ggS=$ggS1
endif

echo " " > tmp/X.$MAGIC/$chrom.tabix.ace
if ($ggS != "") then
  echo "Map $chrom\nWiggle $ggS $ggS/$chrom\n\n" >> tmp/X.$MAGIC/$chrom.tabix.ace
endif
if ($ggNS != "") then
  echo "Map $chrom\nWiggle $ggNS $ggNS/$chrom\n\n" >> tmp/X.$MAGIC/$chrom.tabix.ace
endif
echo "parse the tabix file"
echo "pparse tmp/X.$MAGIC/$chrom.tabix.ace" | bin/tacembly tmp/X.$MAGIC/$chrom -no_prompt
set WGR="toto"
if (-d  tmp/WIGGLEGROUP/$ggNS) set WGR=WIGGLEGROUP

# restore the original data
echo "restore the original data"
if (-e tmp/$WGR/$ggNS/$chrom/R.chrom.u0.f.BF.gz) then
  \mv  tmp/$WGR/$ggNS/$chrom/R.chrom.u0.f.BF.gz  tmp/$WGR/$ggNS/$chrom/R.chrom.u.f.BF.gz
  \mv  tmp/$WGR/$ggNS/$chrom/R.chrom.u0.r.BF.gz  tmp/$WGR/$ggNS/$chrom/R.chrom.u.r.BF.gz
endif

###### grab the coverons and the gene-ends using multiPeaks (XG stranded, XH non stranded)
\rm  tmp/X.$MAGIC/$chrom/f3.Xends.*
set WGR=toto
echo "grab the XG/XH: $ggS/$ggNS"
foreach XGH (XG XH)
  if ($XGH == XG) set stranding=98
  if ($XGH == XG && $ggS == toto) continue 
  if ($XGH == XH && $ggNS == toto) continue 
  set ggs=$ggS
  if ($XGH == XH) set ggs=$ggNS
  if ($XGH == XH) set stranding=50
  if ($XGH == XH && $ggNS == $ggS) continue

  if (-d  tmp/WIGGLERUN/$ggs) set WGR=WIGGLERUN
  if (-d  tmp/WIGGLEGROUP/$ggs) set WGR=WIGGLEGROUP
  if (! -d  tmp/$WGR/$ggs) continue
  set ok=0

echo "# exon multipeaks"
  echo "Construct the multipeaks  $WGR/$ggs"
  echo "  bin/wiggle  -multiPeaks 4 -minCover $minExonCover -I BF -O COUNT -wiggle1 tmp/$WGR/$ggs/$chrom/R.chrom.u.f.BF.gz  -wiggle2 tmp/$WGR/$ggs/$chrom/R.chrom.u.r.BF.gz -stranding 98 -o  tmp/X.$MAGIC/$chrom/f3.$XGH.f"
echo "# grab the exons from the wiggle "
  bin/wiggle  -multiPeaks 4 -minCover $minExonCover -I BF -O COUNT -wiggle1 tmp/$WGR/$ggs/$chrom/R.chrom.u.f.BF.gz  -wiggle2 tmp/$WGR/$ggs/$chrom/R.chrom.u.r.BF.gz -stranding 98 -o  tmp/X.$MAGIC/$chrom/f3.$XGH.f
  bin/wiggle  -multiPeaks 4 -minCover $minExonCover -I BF -O COUNT -wiggle1 tmp/$WGR/$ggs/$chrom/R.chrom.u.r.BF.gz  -wiggle2 tmp/$WGR/$ggs/$chrom/R.chrom.u.f.BF.gz -stranding 98 -o  tmp/X.$MAGIC/$chrom/f3.$XGH.r

  cat tmp/X.$MAGIC/$chrom/f3.$XGH.f.multiPeaks ZZZZZ | scripts/tab_sort -k 1,1 -k 2n |  gawk -F '\t' '/^#/{next;}{if($6+0<0)next;if(a10<1)a10=$2;if($1 != chrom || $2 >= a2+10){chrom=$1;if(a2-a1 > 0){color=1;s=score;while (color<7 && s>100){s/=5;color++;}if(score>=1)printf("Sequence %s_%s__%d_%d\nForward\ncDNA_clone %s_%s__%d_%d\nIntMap %s %d %d\nIs_read\nTags %d\nColour Green%d%s\n\n",XGH,c,a10,a2,XGH,c,a10,a2,c,a10,a2,score,color,cov);cov="";a10=$2;score=0;}}c=$1; a1=$2;a2=$3;if(a1<1)a1=1;s=int($6);if (s>0)cov= cov "\nComposite " a1 - a10 + 1 " " a2 - a10 " " s; if(s>score)score=s;} ' XGH=$XGH >  tmp/X.$MAGIC/$chrom/f3.$XGH.f.ace
  wc  tmp/X.$MAGIC/$chrom/f3.$XGH.f.ace
  cat tmp/X.$MAGIC/$chrom/f3.$XGH.f.ace  | gawk '/^Sequence/{split($2,aa,"__");chrom=substr(aa[1],4);split(aa[2],bb,"_");a1=bb[1];a2=bb[2];ln=a2-a1;if(ln<0)ln=-ln;ln++;printf("%s\t1\t%d\t%s\t%d\t%d\n",$2,ln,chrom,a1,a2);}' >  tmp/X.$MAGIC/$chrom/f3.$XGH.f.shadow
  bin/dna2dna  -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow   tmp/X.$MAGIC/$chrom/f3.$XGH.f.shadow | gawk '/^>/{print;next;}{printf("%s\n",$1);}' >  tmp/X.$MAGIC/$chrom/f3.$XGH.f.fasta

  cat tmp/X.$MAGIC/$chrom/f3.$XGH.r.multiPeaks | scripts/tab_sort -k 1,1 -k2nr > tmp/X.$MAGIC/$chrom/f3.$XGH.r.multiPeaks.r
  cat tmp/X.$MAGIC/$chrom/f3.$XGH.r.multiPeaks.r ZZZZZ | gawk -F '\t' '/^#/{next;}{if($6+0<0)next;if(a10<1){a1=$3;a10=$3;}if($1!=chrom || $3 <= a2-10){chrom=$1;if(a2-a1 < 0){color=1;s=score;while (color<7 && s>100){s/=5;color++;}if(score>=1)printf("Sequence %s_%s__%d_%d\nForward\ncDNA_clone %s_%s__%d_%d\nIntMap %s %d %d\nIs_read\nTags %d\nColour Green%d%s\n\n",XGH,c,a10,a2,XGH,c,a10,a2,c,a10,a2,int(score),color,cov);cov="";a10=$3;score=0;}}c=$1; a1=$3;a2=$2;if(a2<1)a2=1;s=int($6);if(s>0)cov= cov "\nComposite " a10 -a1  + 1 " " a10 -a2 " " s; if(s>score)score=s; }' XGH=XG >  tmp/X.$MAGIC/$chrom/f3.$XGH.r.ace
  wc  tmp/X.$MAGIC/$chrom/f3.$XGH.r.ace
  cat tmp/X.$MAGIC/$chrom/f3.$XGH.r.ace  | gawk '/^Sequence/{split($2,aa,"__");chrom=substr(aa[1],4);split(aa[2],bb,"_");a1=bb[1];a2=bb[2];ln=a2-a1;if(ln<0)ln=-ln;ln++;printf("%s\t1\t%d\t%s\t%d\t%d\n",$2,ln,chrom,a1,a2);}' >  tmp/X.$MAGIC/$chrom/f3.$XGH.r.shadow
  bin/dna2dna  -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow   tmp/X.$MAGIC/$chrom/f3.$XGH.r.shadow | gawk '/^>/{print;next;}{printf("%s\n",$1);}' >  tmp/X.$MAGIC/$chrom/f3.$XGH.r.fasta

# transcriptsEnds :  do not use -stranding option, this creates end echos rather than damping them

    echo "Construct the transcriptsEnds  $WGR/$ggs"
    echo "  bin/wiggle  -transcriptsEnds tmp/$WGR/$ggs/$chrom/R.chrom.u -gzi -I BF -O COUNT -o tmp/X.$MAGIC/$chrom/f3.Xends.$ggs -minCover $minExonCover -wiggleRatioDamper 5"
            bin/wiggle  -transcriptsEnds tmp/$WGR/$ggs/$chrom/R.chrom.u -gzi -I BF -O COUNT -o tmp/X.$MAGIC/$chrom/f3.Xends.$ggs -minCover $minExonCover -wiggleRatioDamper 5

end

# cumulate X and XG ?  
foreach fr (ELF ELR ERF ERR)
    cat tmp/X.$MAGIC/$chrom/f3.Xends.*.$fr.transcriptsEnds >> tmp/X.$MAGIC/$chrom/f3.Xends.$fr.transcriptsEnds
end


foreach fr (ELF ELR ERF ERR)
  set ff=tmp/X.$MAGIC/$chrom/f3.Xends.$fr.transcriptsEnds
  if (-e  $ff) then

     if ($fr == ELF) then
       cat $ff | gawk -F '\t' '/^#/{next;}{if($2<1)$2=1;printf("Xends_%s.%s__%d_%d\t1\t%d\t%s\t%d\t%d\t%d\n", fr,$1,$2,$3,$4,$1,$2,$3,$7/100);}' fr=$fr > $ff.shadow
       cat $ff.shadow | gawk -F '\t' '{printf("Sequence %s\ncDNA_clone %s\nIs_read\nIntMap %s %d %d\nColour PALECYAN\n-D Composite\nComposite %d\nForward\nmForward\n\n", $1, $1,$4,$5,$6,$7);}' > $ff.ace
     endif
     if ($fr == ERF) then
       cat $ff | gawk -F '\t' '/^#/{next;}{if($2<1)$2=1;printf("Xends_%s.%s__%d_%d\t1\t%d\t%s\t%d\t%d\t%d\n", fr,$1,$3,$2,$4,$1,$3,$2,$7/100);}' fr=$fr > $ff.shadow
       cat $ff.shadow | gawk -F '\t' '{printf("Sequence %s\ncDNA_clone %s\nIs_read\nIntMap %s %d %d\nColour PALEYELLOW\n-D Composite\nComposite %d\nReverse\nmReverse\n\n", $1, $1,$4,$5,$6,$7);}' > $ff.ace
     endif
     if ($fr == ERR) then
       cat $ff | gawk -F '\t' '/^#/{next;}{if($2<2)$1=1;printf("Xends_%s.%s__%d_%d\t1\t%d\t%s\t%d\t%d\t%d\n", fr,$1,$3,$2,$4,$1,$3,$2,$7/100);}' fr=$fr > $ff.shadow
       cat $ff.shadow | gawk -F '\t' '{printf("Sequence %s\ncDNA_clone %s\nIs_read\nIntMap %s %d %d\nColour PALECYAN\n-D Composite\nComposite %d\nForward\nmForward\n\n", $1, $1,$4,$5,$6,$7);}' > $ff.ace
     endif
     if ($fr == ELR) then
       cat $ff | gawk -F '\t' '/^#/{next;}{if($2<1)$2=1;printf("Xends_%s.%s__%d_%d\t1\t%d\t%s\t%d\t%d\t%d\n", fr,$1,$2,$3,$4,$1,$2,$3,$7/100);}' fr=$fr > $ff.shadow
       cat $ff.shadow | gawk -F '\t' '{printf("Sequence %s\ncDNA_clone %s\nIs_read\nIntMap %s %d %d\nColour PALEYELLOW\n-D Composite\nComposite %d\nReverse\nmReverse\n\n", $1, $1,$4,$5,$6,$7);}' > $ff.ace
     endif

     bin/dna2dna  -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow $ff.shadow > $ff.fasta
  endif
end


    cat   tmp/X.$MAGIC/$chrom/f3.X*.ace >    tmp/X.$MAGIC/$chrom/f3.XFC2.any.ace 
    cat   tmp/X.$MAGIC/$chrom/f3.X*.shadow >    tmp/X.$MAGIC/$chrom/f3.XFC2.any.shadow
    bin/dna2dna  -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow   tmp/X.$MAGIC/$chrom/f3.XFC2.any.shadow >    tmp/X.$MAGIC/$chrom/f3.XFC2.any.fasta



laba:

  bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
    pparse tmp/X.$MAGIC/$chrom/f3.introns.ace
    pparse tmp/X.$MAGIC/$chrom/f3.XG.f.ace 
    pparse tmp/X.$MAGIC/$chrom/f3.XG.f.fasta
    pparse tmp/X.$MAGIC/$chrom/f3.XG.r.ace 
    pparse tmp/X.$MAGIC/$chrom/f3.XG.r.fasta
    pparse tmp/X.$MAGIC/$chrom/f3.XH.f.ace 
    pparse tmp/X.$MAGIC/$chrom/f3.XH.f.fasta
    pparse tmp/X.$MAGIC/$chrom/f3.XH.r.ace 
    pparse tmp/X.$MAGIC/$chrom/f3.XH.r.fasta
    pparse tmp/X.$MAGIC/$chrom/f3.XFC2.any.ace 
    pparse tmp/X.$MAGIC/$chrom/f3.XFC2.any.fasta
    pparse tmp/X.$MAGIC/$chrom/f3.introns.ace
    pparse tmp/X.$MAGIC/$chrom/f3.transcriptsIntronSupport.ace
    query find Sequence IS Xends_ELF.* && Forward
    edit Real_5prime
    edit Colour PALECYAN
    query find Sequence IS Xends_ERR.* && Forward
    edit Real_5prime
    edit Colour PALECYAN
    query find Sequence IS Xends_ERF.* && Reverse
    edit Real_3prime
    edit Colour PALEYELLOW
    query find Sequence IS Xends_ELR.* && Reverse
    edit Real_3prime
    edit Colour PALEYELLOW
    query find Run IS $ggs && ! W_colour_plus
    edit W_colour_plus BLUE
    edit W_colour_minus GREEN
    query find Runs Wiggle && ! W_stranded == $MAGIC &&  ! W_new_exon == $MAGIC
    edit -D Wiggle
    save
    quit
EOF



hello:

# these lines add the RefSeq mapping to this chromosome as additional pseudo composite reads
if (0) then
# this is FATAL ERROR if the NM are also entered as gene models, becuase it contradicts the genome because of the introns
  cat tmp/METADATA/RefSeq.mrna_map_ln_gc_gene_geneid.txt | gawk -F '\t' '{gsub(/\"/,"",$0); split($2,aa,":");if(aa[1] != chrom)next;split(aa[2],bb,"-");m=$1;c[m]=aa[1];a1[m]=bb[1];a2[m]=bb[2];}END{for(m in c)printf("Sequence %s\nIntMap %s %d %d\nForward\nComposite 10\nIs_read\n\n",m,c[m],a1[m],a2[m]);}' chrom=$chrom | gzip > tmp/X.$MAGIC/$chrom/f3.RefSeq.intmap.ace.gz

  gunzip -c tmp/X.$MAGIC/$chrom/f3.RefSeq.intmap.ace.gz ZZZZZ.gz $ici/TARGET/Targets/$species.RefSeq.fasta.gz |   gawk '/^ZZZZZ/{zz++;next;}/^Sequence/{okk[$2]=1;next;}{if(zz<1)next;}/^>/{s=substr($1,2);i=index(s,"|");if(i>0)s=substr(s,1,i-1);ok=0;if(okk[s]==1){ok=1;print ">" s;}next;}{if(ok==1)print}' | gzip > tmp/X.$MAGIC/$chrom/f3.RefSeq.fasta.gz

  bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
    pparse tmp/X.$MAGIC/$chrom/f3.RefSeq.intmap.ace.gz
    // pparse tmp/X.$MAGIC/$chrom/f3.RefSeq.fasta.gz  
EOF

endif

  bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
    query find sequence is_read && ! IS X?_* && ! intmap==$chrom
    edit -D is_read
    find method encode
    edit Show_up_strand
    query find sequence IS XY_* && (Other  || ct_ac)
    spush
    query find sequence IS XY_* && ! gt_ag && ! gc_ag && COUNT {> intron ; gt_ag || gc_ag} == 0
    sor
    query find sequence IS XW_* && (Other  || ct_ac)
    sor
    query find sequence IS XW_* && ! gt_ag && ! gc_ag && COUNT {> intron ; gt_ag || gc_ag} == 0
    sor
    spop
    spush
    follow cdna_clone
    sor
    spop
    kill
    save
    quit
EOF

  bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
    find sequence Is_read && ! cdna_clone
    acem
      cdna_80
      quit
    query find run w_stranded
    edit W_colour
    query find run w_new_exon
    edit W_colour
    save
    quit
EOF

  scripts/f3.kill_ct_ac_introns.tcsh $chrom 1

# transfer the RNA_seq support read in XI->composite back into the Intron class, which is visible in the tg
bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
  select -o tmp/X.$MAGIC/$chrom/f3.intron_support.txt  xi, ii, n from xi in ?Sequence where xi ~ "XI_*", n in xi->composite, ii in xi->intron where n && ii
  date
EOF
cat tmp/X.$MAGIC/$chrom/f3.intron_support.txt | gawk -F '\t' '{nn[$2]+=$3;}END{for (ii in nn)printf("%s\t%d\n",ii,nn[ii]);}' | sort -V | gawk -F '\t' '{printf ("Intron %s\nRNA_seq %d\n\n",$1,$2);}' >  tmp/X.$MAGIC/$chrom/f3.intron_support.ace
bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
  parse tmp/X.$MAGIC/$chrom/f3.intron_support.ace
  save
  quit
EOF


##### parse the RefSeq
if (-e TARGET/GENES/Ecocyc.ace) then
  echo "pparse TARGET/GENES/Ecocyc.ace" | bin/tacembly tmp/X.$MAGIC/$chrom -no_prompt
endif

foreach target ($RNAtargets)
  if ($target == av) continue
  if ($target == magic) continue
  if ($target == rrna) continue
  if (-e tmp/METADATA/gtf.$target.transcripts.ace.gz) then
    bin/tacembly tmp/X.$MAGIC/$chrom <<EOF  
      query find predicted_gene method == $target
      kill
      pparse TARGET/GENES/RvY.LocusLink.ace
      pparse tmp/METADATA/gtf.$target.transcripts.ace.gz 
      query is_predicted_gene
      edit -D method
      edit method $target
      spush
      query intmap == $chrom
      sminus
      spop
      kill
      parse metaData/methods.ace
      save
      quit
EOF
  endif
end


# only keep the group wiggles
bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
  query find run wiggle  && ! union_of
  edit -D wiggle
  save
  quit
EOF



# les introns mangent 1% de leur soutien sur les zones inclues des XH sur les 2 brins

bin/tacembly tmp/X.$MAGIC/$chrom <<EOF
    pparse yk2sl.xsl.$chrom.ace
    pparse yk2sl.xsl.$chrom.fasta
    save
  quit
EOF


  echo "parse f4.genes2intmap.ace when available"
  pushd tmp/X.$MAGIC/$chrom
    if (-e TABIX) \rm TABIX
    ln -s ../../TABIX
    if (! -e tables) ln -s ../../../metaData/tables
  popd
  if ($species == worm && ! -e tmp/X.$MAGIC/$chrom/genes.ace) then
    tbly ~/yknew <<EOF
      query find gene transcribed_gene && IntMap == $chrom
      show -a -f tmp/X.$MAGIC/$chrom/f4.genes2intmap.ace IntMap
EOF
    if (-e tmp/X.$MAGIC/$chrom/f4.genes2intmap.ace) then
      cat  tmp/X.$MAGIC/$chrom/f4.genes2intmap.ace | gawk '/^Gene/{g=$2;next;}/^IntMap/{printf("Sequence %s\nGenes %s %s %s\n\n", $2, g, $3,$4);}' >   tmp/X.$MAGIC/$chrom/f4.geneParents.ace
      bin/tacembly tmp/X.$MAGIC/$chrom << EOF
        read-models
        parse tmp/X.$MAGIC/$chrom/f4.genes2intmap.ace
        parse tmp/X.$MAGIC/$chrom/f4.geneParents.ace
        save
      quit
EOF
    endif
  endif



touch tmp/X.$MAGIC/$chrom/f3.parse.done
endif

exit 0

echo 284702 > titi
echo ZZZZZ >> titi

gunzip -c PUBLIC/SEQC_NB_MAV_G_log2.20121127.txt.gz | grep  284702 >> titi


cat titi | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$1]=1;next;}}{n=split($2,aa,";");for(i=1;i<=n;i++)if(gg[aa[i]]==1)print;}' 


gunzip -c  titi.gz  ZZZZZ.gz PUBLIC/SEQC_NB_MAV_G_log2.20121127.txt.gz


  gunzip -c RESULTS/MAV_lostGeneID.txt.gz ZZZZZ.gz PUBLIC/SEQC_NB_MAV_G_log2.20121127.txt.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$1]=1;next;}}{n=split($2,aa,";");for(i=1;i<=n;i++)if(gg[aa[i]]==1)print;}' >> RESULTS/MAV_lostGeneID.SEQC_NB_MAV_G_log2.20121127.txt



cat tatou.17.noncoding.list ZZZZZ tatou.201.coding.list ZZZZZ tatou.29.spliced.list ZZZZZ tatou.6.pseudogenes.list | gawk '/^ZZZZZ/{zz++;next;}/^Sequence/{if(zz==0)g[$2]=g[$2] " non_coding";if(zz==1)g[$2]=g[$2] " coding";if(zz==2)g[$2]=g[$2] " spliced";if(zz==3)g[$2]=g[$2] " pseudogene";}END{for(gg in g)printf("%s\t%s\n",gg,g[gg]);}' | sed -e 's/\"//g' -e 's/X__//' | gawk '{printf("%s\t",$1);if(index($0,"coding")>0)nam="Coding"; else nam="Non_coding" ; if(index($0,"spliced")>0)nam= nam "_spliced"; else nam = nam "_single_exon" ;printf("%s\n",nam);}' | sort


  
gunzip -c RESULTS/MAV_lostGeneID.txt.gz ZZZZZ.gz > tutu
cat tutu tmp/METADATA/av.metadata.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$1]=1;next;}}{n=split($3,aa,";");for(i=1;i<=n;i++)if(gg[aa[i]]==1)print;}'


# table transcripts2spliced_RefSeq_coding.def used on ace11
scp transcripts2spliced_RefSeq_coding.def ace11:SERVER

cat transcripts2spliced_RefSeq_coding.txt | gawk -F '\t' '{gsub(/\"/,"",$0);g=$1;gg[g]=1;if($2 != "NULL")spl[g]=1;if($3 != "NULL" && index(rs[g],$3)==0)rs[g]=rs[g] "_" $3;if($4 != "NULL")coding[g]=1;if($5 != "NULL")gene[g]=$5;}END{for (g in gg){if(coding[g]==1) nam="Coding"; else nam="Non_coding" ; if(spl[g]==1)nam= nam "_spliced"; else nam = nam "_single_exon" ; if(length(rs[g])>1)nam= nam "_with_RefSeq" ; printf("%s\t%s\t%s\t%s\n",gene[g],g,nam,substr(rs[g],2));}}' | sort | gzip > transcripts2spliced_RefSeq_coding.txt.1.gz

gunzip -c TARGET/Targets/hs.av.fasta.gz  ZZZZZ.gz transcripts2spliced_RefSeq_coding.txt.1.gz | gawk '/^ZZZZZ/{zz++;next;}/^#/{next;}/^>/{split(substr($1,2),aa,"|");s=aa[1];ss[s]=1;next;}{if(zz<1)next;}{k=0;s=$2;if(ss[s]==1)k="kept";else k="rejected";printf("%s\t",k);print;}' > tutu
cat tutu | sort > RESULTS/transcripts2spliced_RefSeq_coding.kept_rejected.txt

##### expression of a few genes
pushd RESULTS/Expression/unique/RefSeq

echo ZZZZZ > ZZZZZ
foreach ff (`ls $MAGIC.*.txt  $MAGIC.*.txt.gz `)
  gunzip -c -f SELECTED/gene.list ZZZZZ $ff | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz+0<1){gg[$1]=1;next;}}/^#/{print;next;}{if(gg[$1] + gg[$2]> 0)print;next;}' > SELECTED/$ff
end
  
done:
echo 'done'

exit 0

set ggNS=zozo
foreach chrom ($chromSetAll)
  mv tmp/TABIX/$ggNS/$chrom.u0.f.tabix.gz    tmp/TABIX/$ggNS/$chrom.u.f.tabix.gz
  mv tmp/TABIX/$ggNS/$chrom.u0.r.tabix.gz     tmp/TABIX/$ggNS/$chrom.u.r.tabix.gz
  mv tmp/TABIX/$ggNS/$chrom.u0.f.tabix.gz.tbi tmp/TABIX/$ggNS/$chrom.u.f.tabix.gz.tbi
  mv tmp/TABIX/$ggNS/$chrom.u0.r.tabix.gz.tbi tmp/TABIX/$ggNS/$chrom.u.r.tabix.gz.tbi
end



