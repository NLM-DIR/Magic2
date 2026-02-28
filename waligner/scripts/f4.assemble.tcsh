#!bin/tcsh -f

set chrom=$1
setenv ici `pwd`

echo -n "f4.assemble.tcsh start "
date
goto laba0
laba0:
if (! -e tmp/X.$MAGIC/$chrom/f3.parse.done) then
  echo "missing file tmp/X.$MAGIC/$chrom/f3.parse.done, please run phase f3"
  exit 1
endif

  bin/tacembly tmp/X.$MAGIC/$chrom << EOF
    read-models
    query find sequence Is_read && ! cdna_clone
    comment "cdna_80 create cDNA_clone info"
    acembly
      cdna_80
      quit
    comment "cdna_77 Kill premRNA echoes on other strand"
    acembly
      cdna_77
      quit
    save
    comment "kill all old tg and old mrnas"
    find tg
    kill
    find mrna
    kill
    comment "cdna_1 align all est"    
    acembly
      cdna_1 $chrom
      quit
    save
    comment "cdna_71 realign all tg -split_cloud"
    acembly
      cdna_71 -split_cloud
      quit
    save

    comment "test XG"
    query find sequence XG_chr3__1976795_1970805 ; from_gene
    list
    comment "export suspect polyA as tables"
    table -o tmp/X.$MAGIC/$chrom/f4.10.polyAsuspect1.txt  -f tables/10.polyAsuspect1.def
    table -o tmp/X.$MAGIC/$chrom/f4.10.polyAsuspect2.txt  -f tables/10.polyAsuspect2.def
    table -o tmp/X.$MAGIC/$chrom/f4.10.polyAsuspect4.txt  -f tables/10.polyAsuspect4.def
    quit
EOF

laba1:
## remove echo introns
## we need the intron coordinates, some are missing
  echo "find intron NOT IntMap"
  bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
     query find intron ! IntMap
     list -a -f  tmp/X.$MAGIC/$chrom/f4.intron_nomap.list
     quit
EOF
  touch  tmp/X.$MAGIC/$chrom/f4.intron_nomap.list
  ls -ls  tmp/X.$MAGIC/$chrom/f4.intron_nomap.list

  cat  tmp/X.$MAGIC/$chrom/f4.intron_nomap.list | gawk '/^Intron/{z=$2;gsub(/\"/,"",z);n1=split (z,aa,"__");n2=split(aa[2],bb,"_") ;printf("Intron %s\nIntMap %s %s %s\n\n",$2,aa[1],bb[1],bb[2]);}' > tmp/X.$MAGIC/$chrom/f4.intron_nomap.ace
  bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
     pparse  tmp/X.$MAGIC/$chrom/f4.intron_nomap.ace
     save
     query find intron ! IntMap
     quit
EOF


if (-e tmp/X.$MAGIC/$chrom/database/lock.wrm) exit 1
echo "kill echo ct_ac introns"
scripts/f3.kill_ct_ac_introns.tcsh $chrom 2
# find all antisens gene, with non classic or ct_ac intron supported 5 times less than an approximately antisense g[tc]_ag intron
  bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
     table -o tmp/X.$MAGIC/$chrom/f4.killEchoIntron.out -f tables/f4.killEchoIntron.def
     quit
EOF
echo " export list of reads and gene"
cat tmp/X.$MAGIC/$chrom/f4.killEchoIntron.out | gawk -F '\t' '/^"/{printf("Sequence %s\n",$3);printf("Transcribed_gene %s\n",$1);}' > tmp/X.$MAGIC/$chrom/f4.killEchoIntron.list
# kill the bad reads, recompute the genes, they will usually vanish
echo "kill tmp/X.$MAGIC/$chrom/f4.killEchoIntron.list"
if (-e tmp/X.$MAGIC/$chrom/database/lock.wrm) exit 1

  bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
     key  tmp/X.$MAGIC/$chrom/f4.killEchoIntron.list
     spush
     query find sequence IS XY_* && (Other  || ct_ac)
     sor
     query find sequence IS XW_* && (Other  || ct_ac)
     sor
     follow from_gene
     sor
     spop
     spush
     query CLASS Sequence
     edit -D Is_read
     spop
     query CLASS Transcribed_gene
     list -a -f  tmp/X.$MAGIC/$chrom/f4.killEchoIntron.tg.list
     save
     quit
EOF

laba2:

if (-e tmp/X.$MAGIC/$chrom/database/lock.wrm) exit 1
# split this code out of previous one to monitor eventual errors
  bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
     key  tmp/X.$MAGIC/$chrom/f4.killEchoIntron.tg.list
    list
     acem
        cdna_73 -split_cloud
        quit
    comment "test XG after split73"
    query find sequence XG_chr3__1976795_1970805 ; from_gene
    list
     query find mrna DNA:2 > 10000 ; > from_gene
     acem
        cdna_73 -split_cloud // 2017_03_11 a hack to rm stupid extra long empty exons
        quit
     save
     find clone
     list -a -f tmp/X.$MAGIC/$chrom/f4.killEchoIntron.done
     query find transcribed_gene Intron
     list -a -f f4.spliced_tg.list
     quit
EOF

touch  tmp/X.$MAGIC/$chrom/f4.assemble.done
exit 0  # Composite case, exit now, no need for the complications below
if (-e tmp/X.$MAGIC/$chrom/database/lock.wrm) exit 1
echo -n "f4.killEchoIntron done"
date
scripts/f3.kill_ct_ac_introns.tcsh $chrom 3


if (-e tmp/X.$MAGIC/$chrom/database/lock.wrm) exit 1


foreach chrom ($chromSetAll)
  bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
     query find transcribed_gene Intron
     list -a -f tmp/X.$MAGIC/$chrom/f4.spliced_tg.list
     query find predicted_gene Intron ; >model_of
     list -a -f tmp/X.$MAGIC/$chrom/f4.spliced_pg.list
     quit
EOF
end

wc tmp/X.$MAGIC/*/f4.spliced_pg.list
wc tmp/X.$MAGIC/*/f4.spliced_tg.list

# suspect 1 finds reverse read starting inside the ORF
# suspect 2 finds forward polyA read endding inside the ORF
# suspect 4 removes polyA of forward read if a reverse read assembles further down

echo "reverse polyA suspect reads "
pushd tmp/X.$MAGIC/$chrom
  gawk  -F '\t' '/\"/ {printf ("Sequence %s\nColour PALEGREEN\n\ncDNA_clone %s\nInternal_priming\n\n",$8, $9);}'  f4.10.polyAsuspect1.txt >!  f4.10.polyAsuspect1.ace
  gawk  -F '\t' '/\"/ {printf ("Sequence %s\nColour PALEGREEN\n\ncDNA_clone %s\nInternal_priming\n\n",$8, $9);}'  f4.10.polyAsuspect2.txt >!  f4.10.polyAsuspect2.ace
  gawk  -F '\t' '/\"/ {printf ("Sequence %s\nColour PALEGREEN\n-D polya_after_base\n-D Number_of_terminal_A\n\n",$4);}'  f4.10.polyAsuspect4.txt >!  f4.10.polyAsuspect4.ace

  wc f4.10.polyAsuspect*.txt

  $ici/bin/tacembly . << EOF >! f4.10.reverse.log
    pparse f4.10.polyAsuspect1.ace
    spush 
    pparse f4.10.polyAsuspect2.ace
    sor
    pparse f4.10.polyAsuspect4.ace
    sor
    query find est Number_of_terminal_A && !PolyA_after_base
    edit -D  Number_of_terminal_A // since they were removed by the above scripts
    sor
    query find cdna_clone double_fuzzy
    query follow read ; NOT ref_mrna && NOT ref_seq
    edit -D Is_read
    sor

    comment "test XG in split_cloud "
    query find sequence XG_chr3__1976795_1970805 ; from_gene
    list
    
    comment "restore missing flags gt_ag in est"
    query find est Flipped
    edit -D Flipped
    query find tg ct_ac
    acem
      cdna_73  -split_cloud // restore missing flags gt_ag in est
      quit
    comment "FlipAllGenes"
    acembly
      cDNA_FlipAllGenes
      quit
    query find est Flipped
    sor
    spop
    query follow from_gene
    comment "realign, will kill all the non_best_mrna and resize"
    acem
      cdna_73 -split_cloud // will kill all the non_best_mrna and resize and fuse a first time
      quit 
    comment "test XG in after split_cloud1 "
    query find sequence XG_chr3__1976795_1970805 ; from_gene
    list
    find tg
    comment cDNA_Flag_suspected_internal_deletion
    acem
      cDNA_Flag_suspected_internal_deletion -ignore // ignore non informative clones
      quit
    comment "test XG in after split_cloud2 "
    query find sequence XG_chr3__1976795_1970805 ; from_gene
    list
    query find tg to_be_fused_with
    comment "fuse_genes"
    acem
      cdna_73 -split_cloud // -locally incompatible avec to_be_fused_with
      quit
    comment "test XG in after split_cloud3 "
    query find sequence XG_chr3__1976795_1970805 ; from_gene
    list
    query find mrna NOT from_gene
    spush
    follow dna
    sor
    spop
    kill
    save
    acembly
      cdna_50
      gene_intron
      quit
    save
    quit
EOF

  scripts/f3.kill_ct_ac_introns.tcsh $chrom 4
  scripts/f3.kill_ct_ac_introns.tcsh $chrom 5
  scripts/f3.kill_ct_ac_introns.tcsh $chrom 6

#  /home/mieg/bin/gene2gene.2007_12_17 . // obsolete, replaced by cdna_21 cdna_50

  $ici/bin/tacembly . << EOF
    acembly
      cdna_21
      cdna_50
      quit
    save
    find intron
    list -a -f f4.introns.list
    find clone
    list -a -f f4.assemble.done
    quit
EOF

  echo ' ' >  f4.introns.preace
  echo ' ' >  f4.introns.ace
  if (-d $ici/GeneIndexDB) then
    $ici/bin/tacembly $ici/GeneIndexDB << EOF
      key f4.introns.list
      show -a -f f4.introns.preace
      quit
EOF

    cat  f4.introns.preace | gawk '/^$/{print}{gsub(/\"/,"",$0);}/^Intron/{printf("Intron %s\n",$2);next;}/^Group/{if($2 == nx)print;}/Validated_u/{if($2 == "any")print}/de_[du][nu]o/{if($2 == "any")print}' nx=Rat  > f4.introns.ace
  endif

  $ici/bin/tacembly . << EOF
    pparse f4.introns.ace
    find transcribed_gene
    list -a -f f4.tg.list
    query find transcribed_gene Intron_boundaries
    list -a -f f4.spliced_tg.list
    save
    quit
EOF

if (0) then
  $ici/bin/tacembly . << EOF
    find transcribed_gene
    acem
      cdna_71 -split_cloud
      quit
    save
    quit
EOF
endif

popd

# edit away the ct_ac XY and XI
scripts/f4.fix.tcsh $chrom




echo -n "f4.assemble.tcsh end "
date

exit 0

bin/tace ~/yk <<EOF
  date
  table -f $ici/tables/f5.tg2genebox2product.def -o  t1
  date
  table -f $ici/tables/f5.tg2genebox2product.def -bqlo  t2
  date
  table -f $ici/tables/f5.tg2genebox2product.def -o  t3
  date
  table -f $ici/tables/f5.tg2genebox2product.def -bqlo  t4
  date
EOF
wc t[1234]
ls -ls t?
