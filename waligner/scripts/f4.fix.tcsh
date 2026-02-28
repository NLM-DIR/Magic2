#!bin/tcsh -f

set chrom=$1

bin/tacembly tmp/X.$MAGIC/$chrom << EOF
  query find predicted_gene IS magic.*
  kill
  read-models
  parse metaData/methods.ace
  query find gene capture == A2
  edit Colour YELLOW
  query find gene capture == R3
  edit Colour RED3
  query capture == A2
  edit Colour ORANGE
  query find est XI_*
  edit COLOUR BLUE1
  query Composite > 1
  edit COLOUR BLUE2
  query Composite > 4
  edit COLOUR BLUE3
  query Composite > 16
  edit COLOUR BLUE4
  query Composite > 64
  edit COLOUR BLUE5
  query Composite > 256
  edit COLOUR BLUE6
  query Composite > 1024
  edit COLOUR BLUE7
  query Composite > 4048
  edit COLOUR BLUE8
  save
  quit
EOF






bin/tacembly tmp/X.$MAGIC/$chrom << EOF
  query find read XY* ; ct_ac || other ; ! gt_ag
  edit -D Is_read
  follow from_gene
  acembly
    cdna_73 -clean_killed_mRNA
    quit
  query find tg (other || ct_ac) && ! gt_ag
  follow read
  query ct_ac || other
  edit -D Is_read
  follow from_gene
  acembly
    cdna_73 -clean_killed_mRNA
    quit
  save
  quit
EOF

# clean XH echo of intron
pushd tmp/X.$MAGIC/$chrom
  ../../../bin/tacembly .  << EOF
    query find read XH_* 
    bql -o f4.fix.xh.txt select chr,a2,a1,s,x,n from x in @,chr in x->intmap, a1 in chr[1],a2 in chr[2],s in x->strand,x1 in x->composite,x2 in x1[1],n in x2[1]
    query find read XI_* 
    bql -o f4.fix.xi.txt select chr,a1,a2,s,x,n from x in @,chr in x->intmap, a1 in chr[1],a2 in chr[2],s in x->strand,n in x->composite
    quit
EOF
  cat f4.fix.xi.txt f4.fix.xh.txt | gawk -F '\t' '{chr=$1;a1=$2;a2=$3;if(a1<a2)print;}' | sort -k 1,1 -k 2,2n > f4.fix.xhi.f.txt
  cat f4.fix.xi.txt f4.fix.xh.txt | gawk -F '\t' '{chr=$1;a1=$2;a2=$3;if(a1>a2)print;}' | sort -k 1,1 -k 3,3n > f4.fix.xhi.r.txt

  cat f4.fix.xhi.f.txt | gawk -F '\t' '{c=$1;a1=$2;a2=$3;n=$6;if(substr($5,1,3)== "XI_"){if(c != iC || n>iN || a1>i2-35){iN=n;iC=c;i1=a1;i2=a2}next;}if(substr($5,1,3)== "XH_"){if(c == iC && 20*n<iN&& a1>i1 && a2<i2){printf("Sequence %s\n",$5);}}}' > f4.xh.kill
  cat f4.fix.xhi.r.txt | gawk -F '\t' '{c=$1;a1=$3;a2=$2;n=$6;if(substr($5,1,3)== "XI_"){if(c != iC || n>iN || a1>i2-35){iN=n;iC=c;i1=a1;i2=a2;if(0)print n,i1,i2;}next;}if(substr($5,1,3)== "XH_"){if(0)print $5,n,a1,a2,"hh",iN,i1,i2; if(c == iC && 20*n<iN&& a1>i1 && a2<i2){printf("Sequence %s\n",$5);}}}'   >> f4.xh.kill


  ../../../bin/tacembly .  << EOF
    key f4.xh.kill
    spush
    follow from_gene
    kstore gg
    spop
    kill
    kget gg
    acem
      cdna_73 -locally
      quit
    save
    quit
EOF

# restrand the XH using reconstructed mrna with introns 
  ../../../bin/tacembly .  << EOF
    query find mrna gt_ag > 2 
    bql -o f4.fix.in_gtag.txt select chr,a1,a2,t,g from x in @, t="good",chr in x->intmap, a1 in chr[1],a2 in chr[2],g in x->from_gene
    query find sequence IS XH_* || IS Xends_* ; COUNT {>in_mrna; gt_ag}=0
    bql -o f4.fix.no_gtag.txt select chr,a2,a1,t,s,x,n from x in @, t="bad",chr in x->intmap, a1 in chr[1],a2 in chr[2],s in x->strand,x1 in x->composite,x2 in x1[1],n in x2[1]

EOF
    cat f4.fix.in_gtag.txt f4.fix.no_gtag.txt | gawk -F '\t' '{chr=$1;a1=$2;a2=$3;if(a1<a2){if($4=="good"){$2-=30;$3+=30};print;}}' | sort -k 1,1 -k 2,2n > f4.fix.gt_ag.f.txt
    cat f4.fix.in_gtag.txt f4.fix.no_gtag.txt | gawk -F '\t' '{chr=$1;a1=$2;a2=$3;if(a1>a2){if($4=="good"){$2+=30;$3-=30};print;}}' | sort -k 1,1 -k 3,3n > f4.fix.gt_ag.r.txt

  cat f4.fix.gt_ag.f.txt | gawk  '{c=$1;a1=$2;a2=$3;n=$6;if($4=="good"){if(c== iC && a1<i2){if (i2<a2)a2=i2;g=$5;next}iC=c;i1=a1;i2=a2;g=$5;next;}if(c == iC && a1>i1-30 && a2<i2+30){printf("Sequence %s\nFrom_gene %s\n\n",$6,g);}}' > f4.gt_ag.ace
  cat f4.fix.gt_ag.r.txt | gawk  '{c=$1;a1=$3;a2=$2;n=$6;if($4=="good"){if(c== iC && a1<i2){if (i2<a2)a2=i2;g=$5;next}iC=c;i1=a1;i2=a2;g=$5;next;}if(c == iC && a1>i1-30 && a2<i2+30){printf("Sequence %s\nFrom_gene  %s\n\n",$6,g);}}' >> f4.gt_ag.ace

  ../../../bin/tacembly .  << EOF
    parse f4.gt_ag.ace
    kstore rr
    spush
    query forward
    edit reverse
    sminus
    spop
    edit forward
    kget rr
    follow from_gene
    sor
    spop
    acem
      cdna_73 -locally
      cdna_31
      quit
    query find est IS Xends_* &&  Flipped
    spush
    follow cdna_clone
    sor
    spop 
    kill
    save
    quit
EOF

  

popd
