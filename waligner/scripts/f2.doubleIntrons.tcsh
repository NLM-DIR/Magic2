#!bin/tcsh -f

set chrom=$1
set dd=tmp/X.$MAGIC
set ddc=tmp/X.$MAGIC/$chrom
echo -n "f2 $dd/f2.allDoubleIntrons.txt $chrom start : "
date


  zcat $dd/f2.allDoubleIntrons.txt.gz | gawk '{split($1,aa,"__");if(aa[1]==chrom)print;}' chrom=$chrom > $ddc/f2.XW.txt
  cat $ddc/f2.XW.txt | gawk -F '\t' '{if($3<1)next;split($1,dd,"___");split(dd[1],aa,"__");c=aa[1];split(aa[2],bb,"_");a1=bb[1];a2=bb[2];split(dd[2],bb,"_");b1=bb[1];b2=bb[2];s=1;if(a1>a2)s=-1;printf("Sequence XW_%s__%s_%s_%s_%s\nColour ORANGE\nForward\nComposite %d\nIntron %s__%d_%d\nIntron %s__%d_%d\nIs_read\ncDNA_clone XW__%s_%s_%s_%s_%s\nIntMap %s %d %d\n\n",c,a1,a2,b1,b2,$3,c,a1,a2,c,b1,b2,c,a1,a2,b1,b2,c,a1-30*s,b2+30*s);}' > $ddc/f2.XW.ace
  cat $ddc/f2.XW.txt | gawk -F '\t' '{if($3<1)next;split($1,dd,"___");split(dd[1],aa,"__");c=aa[1];split(aa[2],bb,"_");a1=bb[1];a2=bb[2];split(dd[2],bb,"_");b1=bb[1];b2=bb[2];s=1;if(a1>a2)s=-1;if(a1>0 && b2 + 30*s >0 && b1>0 && b2 > 0 && s*(b1-a2)+1+30 > 0 && s*(b1-a2)+60 > 0){printf("XW_%s__%s_%s_%s_%s\t1\t30\t%s\t%d\t%d\n",c,a1,a2,b1,b2,c,a1-30*s,a1-s);printf("XW_%s__%s_%s_%s_%s\t%d\t%d\t%s\t%d\t%d\n",c,a1,a2,b1,b2,31,29+s*(b1-a2),c,a2+s,b1-s);printf("XW_%s__%s_%s_%s_%s\t%d\t%d\t%s\t%d\t%d\n",c,a1,a2,b1,b2,s*(b1-a2)+1+30,s*(b1-a2)+60,c,b2+s,b2+30*s);}}' > $ddc/f2.XW.shadow
 
    bin/dna2dna -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow  $ddc/f2.XW.shadow > $ddc/f2.XW.fasta




# DOUBLEINTRON de_uno found on the fly by the aligner
echo AAAAA
# export the ace file
set ff=$dd/f2.allDoubleIntronsGenomic.txt
set ffc=$ddc/f2.XY
zcat $ff.gz |  gawk  -F '\t' '/^DOUBLEINTRON/{if($2!=chrom)next;if($9<2)next;print;}' chrom=$chrom > $ffc.txt
cat $ffc.txt |  gawk  -F '\t' '{c=$2;s=1;if($4>$5)next;s=1;dx1=$4-$3+1;dx2=$6-$5+1;dx3=$8-$7+1;di1=$5-$4-1;di2=$7-$6-1;if(di1<0 || di2<0 || dx1<0 || dx2<0 || dx3<0)next;printf ("Sequence XY_%s__%d_%d_%d_%d\n", $2,$4+s,$5-s,$6+s,$7-s); printf ("cDNA_clone XY_%s__%d_%d_%d_%d\n", $2,$4+s,$5-s,$6+s,$7-s); printf("COLOUR CERISE\nForward\nComposite %d\nIs_read\n", $9);printf("IntMap %s %d %d\n",$2,$3,$8);printf("Intron %s__%s_%s\nIntron %s__%d_%d\n\n",$2,$4+s,$5-s,$2,$6+s,$7-s);}'  chrom=$chrom > $ffc.ace
cat $ffc.txt |  gawk  -F '\t' '{c=$2;s=1;if($4<$5)next;s=-1;dx3=$7-$8+1;dx2=$5-$6+1;dx1=$3-$4+1;di2=$6-$7-1;di1=$4-$5-1;if(di1<0 || di2<0 || dx1<0 || dx2<0 || dx3<0)next;printf ("Sequence XY_%s__%d_%d_%d_%d\n", $2,$4+s,$5-s,$6+s,$7-s); printf ("cDNA_clone XY_%s__%d_%d_%d_%d\n",$2,$4+s,$5-s,$6+s,$7-s); printf("COLOUR CERISE\nForward\nComposite %d\nIs_read\n", $9);printf("IntMap %s %d %d\n",$2,$3,$8);printf("Intron %s__%s_%s\nIntron %s__%d_%d\n\n",$2,$4+s,$5-s,$2,$6+s,$7-s);}' >>  $ffc.ace

  # create the shadows
  cat $ffc.txt |  gawk  -F '\t' '{c=$2;s=$1;if($4>$5)next;nam="XY_" c "__" $4+1 "_" $5-1 "_" $6+1 "_" $7-1; dx1=$4-$3+1;dx2=$6-$5+1;dx3=$8-$7+1;di1=$5-$4-1;di2=$7-$6-1;if(di1<0 || di2<0 || dx1<0 || dx2<0 || dx3<0)next;printf("%s\t%d\t%d\t%s\t%d\t%d\n",nam,1,dx1,c,$3,$4);printf("%s\t%d\t%d\t%s\t%d\t%d\n",nam,dx1+1,dx1+dx2,c,$5,$6);printf("%s\t%d\t%d\t%s\t%d\t%d\n",nam,dx1+dx2+1,dx1+dx2+dx3,c,$7,$8);}' chrom=$chrom >  $ffc.shadow
cat $ffc.txt |  gawk  -F '\t' '{c=$2;s=1;if($4<$5)next;nam="XY_" c "__" $4-1 "_" $5+1 "_" $6-1 "_" $7+1; ;s=-1;dx3=$7-$8+1;dx2=$5-$6+1;dx1=$3-$4+1;di2=$6-$7-1;di1=$4-$5-1;if(di1<0 || di2<0 || dx1<0 || dx2<0 || dx3<0)next;printf("%s\t%d\t%d\t%s\t%d\t%d\n",nam,1,dx1,c,$3,$4);printf("%s\t%d\t%d\t%s\t%d\t%d\n",nam,dx1+1,dx1+dx2,c,$5,$6);printf("%s\t%d\t%d\t%s\t%d\t%d\n",nam,dx1+dx2+1,dx1+dx2+dx3,c,$7,$8);}' chrom=$chrom >>  $ffc.shadow

  # export the fasta
  dna2dna -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow $ffc.shadow > $ffc.fasta

#  bug dans dna2dna -shadow      2024/12/08 XY.shadow avec donnees redondantes => fausse insertion dans le fasta du XY
#                                                                                                 gatccaatatgctttatgcagtgtg
#genome                                                                                       tcaggatccaatatgctttatgcagtgtg                            cacatcgcatcagagataactattga
#w    ttgcacaaaatgcacatgtacaagtgatggg gtgaaatgtgctttcctgtctcctgactgctcatctcacacaaaagaggcactggctgaggatccaatatgctttatgcagtgtg gatccaatatgctttatgcagtgtg  cacatc
#w3   gatccaatatgctttatgcagtgtgcacatc gtgaaatgtgctttcctgtctcctgactgctcatctcacacaaaagaggcactggctgag

#########################################
##  XA XSL
# grab the pA and the SL
  echo "XXXXXXXXXXXXXXXXXXXXXX Grab the XA_"


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

ls -ls   tmp/SLpA/$ggNS.SLpA.gz tmp/SLpA/$ggS.SLpA.gz
touch $ddc/f2.XA.ace
touch $ddc/f2.XA.fasta
touch $ddc/f2.XSL.ace
touch $ddc/f2.XSL.fasta
if (-e tmp/SLpA/$ggNS.SLpA.gz || -e tmp/SLpA/$ggS.SLpA.gz) then
  echo "Grab the XA_"
  gunzip -c tmp/SLpA/$ggNS.SLpA.gz tmp/SLpA/$ggS.SLpA.gz | gawk -F '\t' '/^pA/{if($2 == chrom && $5 > 0) {a1=$3+0;if($4=="Forward")a2=a1-30;else a2=a1+30;support=$5;printf("Sequence XA_%s__%s_%d\ncDNA_clone XA_%s__%s_%d\nIntMap %s %d %d\nIs_read\nColour LIGHTORANGE\nComposite %d\nReverse\nmReverse\nPolyA_after_base 9\n\n",chrom,a1,a2,chrom,a1,a2,chrom,a1,a2,support);}}' chrom=$chrom >  $ddc/f2.XA.ace
  cat $ddc/f2.XA.ace | gawk '/^Sequence/{split($2,aa,"__");chrom=substr(aa[1],4);split(aa[2],bb,"_");a1=bb[1];a2=bb[2];ln=a2-a1;if(ln<0)ln=-ln;ln++;printf("%s\t1\t%d\t%s\t%d\t%d\n",$2,ln,chrom,a1,a2);}' >  $ddc/f2.XA.shadow
  bin/dna2dna  -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow $ddc/f2.XA.shadow | gawk '/^>/{print;next;}{printf("AAAAAAAA%s\n",$1);}' >  $ddc/f2.XA.fasta
ls -ls $ddc/f2.XA.*

  echo "Grab the XSL_"
  gunzip -c tmp/SLpA/$ggNS.SLpA.gz tmp/SLpA/$ggS.SLpA.gz | gawk -F '\t' '/^SL/{if($2 == chrom && $5 > 0) {a1=$3+0;if($4=="Forward")a2=a1+30;else a2=a1-30;support=$5;printf("Sequence X%s_%s__%s_%d\ncDNA_clone X%s_%s__%s_%d\nTranspliced_to %s\nIntMap %s %d %d\nIs_read\nColour BLACK\nComposite %d\nForward\nmForward\n\n",$1,chrom,a1,a2,$1,chrom,a1,a2,$1,chrom,a1,a2,support);}}' chrom=$chrom >  $ddc/f2.XSL.ace
  cat $ddc/f2.XSL.ace | gawk '/^Sequence/{split($2,aa,"__");i=index(aa[1],"_");chrom=substr(aa[1],i+1);split(aa[2],bb,"_");a1=bb[1];a2=bb[2];ln=a2-a1;if(ln<0)ln=-ln;ln++;printf("%s\t1\t%d\t%s\t%d\t%d\n",$2,ln,chrom,a1,a2);}' > $ddc/f2.XSL.shadow
  bin/dna2dna  -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow $ddc/f2.XSL.shadow | gawk '/^>/{print;next;}{printf("%s\n",$1);}' >  $ddc/f2.XSL.fasta
ls -ls $ddc/f2.XSL.*
endif



touch  $ddc/f2.done
sleep 10

exit 0
