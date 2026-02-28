#!bin/tcsh -f

set MAGIC=$1
set chrom=$2


if (! -d tmp/X.$MAGIC) mkdir tmp/X.$MAGIC
if (! -d tmp/X.$MAGIC) mkdir tmp/X.$MAGIC
if (! -d tmp/X.$MAGIC/$chrom) mkdir tmp/X.$MAGIC/$chrom

set sMin=3
set sUno=""
set iDuo="tmp/X.$MAGIC/$chrom/f1.XI"
set sDuo="-sxxNewIntronsFileName  $iDuo "

# export intron with at least 1 percent forking at donor and acceptor sites, and the donor and acceptor are at least 10 percent of the reatained intron 
if (! -e $iDuo) then
  bin/tace tmp/INTRON_DB/$chrom <<EOF
    find intron 
    date
    select -o $iDuo ii, g, s, av, t, t2, c, c1, c2, dnaD, dnaA from ii in @ where ! ii#is_echo, g in ii->group where g#Intron && g->Intron[1] == $MAGIC, s in g[1] where s >= $sMin, t in ii->type, t2 in ii->Other, av in ii->stype[1], d in ii->D, a in ii->A, dg in d->group where dg == g, ag in a->group where ag == g, d1 in dg[1], d2 in dg[2], d3 in dg[3], a1 in ag[1], a2 in ag[2], a3 in ag[3] where 10*d3 >= d2 && 10*a3 >= a2 && 100 * s > d3 && 100 *s > a3, dnaD in d->motifs, dnaA in a->motifs, c in ii->IntMap, c1 in c[1], c2 in c[2]
    date
EOF
endif


if (0) then
  bin/tace tmp/INTRON_DB/chrIII <<EOF
    find intron *_834966_*
    select  ii, g, s, av, t, t2, c, c1, c2, from ii in @ where ! ii#is_echo, g in ii->group where g#Intron && g->Intron[1] == "Good_Long279V", s in g[1] where s >= 3, t in ii->type, t2 in ii->Other, av in ii->stype[1], d in ii->D, a in ii->A, dg in d->group where dg == g, ag in a->group where ag == g, d1 in dg[1], d2 in dg[2], d3 in dg[3], a1 in ag[1], a2 in ag[2], a3 in ag[3] where 10*d3 >= d2 && 10*a3 >= a2 , c in ii->IntMap, c1 in c[1], c2 in c[2]
EOF
endif




cat $iDuo |  gawk -F '\t' '{ii=$1;g=$2;s=$3;av=$4;t=$5;t2=$6;c=$7;c1=$8+0;c2=$9+0;if(c1<c2){c1-=35;c2+=35;}else{c1+=35;c2-=35}printf("Sequence XI_%s_%s\n",g,ii);printf("cDNA_clone XI_%s_%s\nIs_read\nForward\nComposite %d\n",g,ii,s);if(av!="NULL")print av;if(t2!="NULL")printf("Other %s\n",t2);else print t;printf("Intron %s\nIntMap %s %d %d\n\n", ii, c,c1,c2);}'  | gzip > $iDuo.ace.gz

cat $iDuo |  gawk -F '\t' '{ii=$1;g=$2;dD=$10;dA=$11;i=index(dD,"--");sD=substr(dD,i-35,35);i=index(dA,"--");sA=substr(dA,i+2,35);printf(">XI_%s_%s\n%s%s\n",g,ii,sD,sA);}' |  gzip > $iDuo.fasta.gz
touch tmp/X.$MAGIC/$chrom/f1.done

exit 0
 # sponge_exon sponge_intron splice splice_percent
  select intron chrIII__12040_12350
  select ii, g, s, av, t, d1, d2, d3, a1, a2, a3 from ii in @ where ! ii#is_echo, g in ii->group where g#Intron && g->Intron[1] == Good_Long279V, s in g[1] where s >= 3, t in ii->type, av in ii->stype[1], d in ii->D, a in ii->A, dg in d->group where dg == g, ag in a->group where ag == g, d1 in dg[1], d2 in dg[2], d3 in dg[3], a1 in ag[1], a2 in ag[2], a3 in ag[3] where 10*d3 >= d2 && 10*a3 >= a2

##########################################################################
##########################################################################


###
if (0) then
set minIntronSupport=`cat MetaDB/$MAGIC/RunList | gawk '{n++}END{printf("%d",int(n/20)+1);}'`
if (-e  tmp/introns/d4.$MAGIC.candidate_introns.ace) then
  set n=`cat tmp/introns/d4.$MAGIC.candidate_introns.ace | gawk '/^Ali/{ok=1;}/^Candidate_introns/{if(ok==1)n=$19+0;last;}END{print n+0}' `
  if ($n >0) set minIntronSupport=$n
endif

set sMin="-minimalIntronSupport $minIntronSupport"
echo "Minimal Intron support sMin=$sMin"
if (! -d  tmp/X.$MAGIC/$chrom) mkdir  tmp/X.$MAGIC/$chrom

set sDuo=" -newDoubleIntrons "
  if (! -e tmp/X.$MAGIC/$chrom/introns.de_duo) then
  echo "collating de_duo"

    # touch  tmp/X.$MAGIC/any.$chrom.limit1
    foreach run (`cat MetaDB/$MAGIC/RunList`)
        # was OR pA pT SL1 SL2 SL3 SL4 SL5 SL6 SL7 SL8 SL9 SL10 SL11 SL12, 2016_04_01 we now prefer to do the XA in phase f3
       foreach dd (OR)
         echo $run $chrom $dd
         set ok=0
         foreach n (1 2 5 10 20 50 100 200 1000)
           if ($ok == 1) continue
           if (-e tmp/$dd/$run/$chrom/$run.$n.0.txt) then
             set ok=1
             cat  tmp/$dd/$run/$chrom/$run.$n.*.txt | gawk -F '\t' '{if($3 != chrom)next;}/^pA/{print}/^pT/{print}/^SL/{print}/^INTRON/{print}'  chrom=$chrom |  sort   >>  tmp/X.$MAGIC/$chrom/introns.de_duo
           endif
         end
       end
    end

endif

set  iDuo=tmp/X.$MAGIC/$chrom/introns.de_duo
if (-e $iDuo) then
  set sDuo="-sxxNewIntronsFileName  $iDuo "
  ls -ls $iDuo
endif


deUno:

if (! -e tmp/X.$MAGIC/$chrom/f1.introns.de_uno.gz) then
  echo "collating de_uno"
  if (-e tmp/X.$MAGIC/$chrom/f1.introns.de_uno.1) \rm tmp/X.$MAGIC/$chrom/f1.introns.de_uno.*
  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    gunzip -c tmp/OR/$run/d4.de_uno.txt.gz | gawk -F '\t' '{if($1==chrom)print;}' chrom=$chrom >> tmp/X.$MAGIC/$chrom/f1.introns.de_uno.1
  end

  cat tmp/X.$MAGIC/$chrom/f1.introns.de_uno.1 | gawk -F '\t' '{z=$1 "\t" $2 "\t" $3 ; n[z] += $4;if (length($5)>1) t[z]=$5;}END{for (k in n) printf("%s\t%d\t%s\n",k,n[k],t[z]);}' >  tmp/X.$MAGIC/$chrom/f1.introns.de_uno.2
  cat  tmp/X.$MAGIC/$chrom/f1.introns.de_uno.2 | gawk  -F '\t' '{chrom=$1;a1=$2;a2=$3;nn=$4;type=$5;ln=a2-a1;if(ln<0)ln=-ln;ln+=1;if(a1<a2){s="Forward";b1=a1-35;b2=a1-1;c1=a2+1;c2=a2+35;}else{s="Reverse";b1=a1+35;b2=a1+1;c1=a2-1;c2=a2-35;}printf("INTRON\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t\t%d\t%d\n",s,chrom,b1,b2,chrom,c1,c2,ln,nn);}' >  tmp/X.$MAGIC/$chrom/f1.introns.de_uno

  \rm  tmp/X.$MAGIC/$chrom/f1.introns.de_uno.[12]
  gzip tmp/X.$MAGIC/$chrom/f1.introns.de_uno
endif

set iUno=tmp/X.$MAGIC/$chrom/f1.introns.de_uno.gz

set sUno=""

if (-e $iUno) then
  set sUno="-sxxDeUnoIntronsFileName  $iUno"
  ls -ls $iUno
endif

# if ($species == worm)   set chr="CHROMOSOME_"

endif
###

f1_txts:
if (! -e tmp/X.$MAGIC/$chrom/f1.txts) then

    echo "running bin/geneelements"
    echo " bin/geneelements -newDoubleIntrons -stranded $sUno $sDuo $sMin -sxxChromosome $chrom -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz"
           bin/geneelements -newDoubleIntrons -stranded $sUno $sDuo $sMin -sxxChromosome $chrom -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz  | sort -u  > tmp/X.$MAGIC/$chrom/f1.txts
    ls -ls  tmp/X.$MAGIC/$chrom/f1.txts

endif


echo "collate ace       "
# intron name $8 $9 was u1 u2   2024_02 20
  cat   tmp/X.$MAGIC/$chrom/f1.txts | gawk -F '\t' '/^#/{next}/^\/\//{next}{chrom=$1;if (0) gsub(/CHROMOSOME_/,"",chrom);u1=$6;u2=$11;dnaa=$12;dnab=$13;dnac=$14;support=$15;feet=$16;if(length(feet)==5){if(feet == "gt_ag" || feet == "gc_ag" || feet == "at_ac" || feet == "ct_ac"){feet=feet  "\n";}else {feet= "Other "  feet  "\n" ;}}else feet="Julot\n" ;typea=$3;typeb=$4;typec=$5; if(typeb == "Exon")next;xx="XI_" ; col="PALEYELLOW";nam=xx group "_" chrom "__" $8 "_" $9 ;names[nam]++; n = names[nam] ;nam = nam "." n ; printf("Sequence %s\ncDNA_clone %s\nIntMap %s %d %d\nIs_Read\nForward\nComposite %d\nColour %s\n%s",nam, nam, $1,u1,u2,support, col,feet) ;sl=0;if(substr(typea,1,2)=="SL"){x1=length(dnaa);sl=substr(typea,3,1)+0;}if(substr(typeb,1,2)=="SL"){x1=length(dnaa);sl=substr(typeb,3,1)+0;};if(sl>0)printf("Transpliced_to SL%d %d\n",sl,x1);pA=0;if(typeb=="pA"){x1=length(dnaa);}if(typec=="pA"){x1=length(dnaa)+length(dnab);pA=length(dnac);};if(pA>0)printf("PolyA_after_base %d\n",x1);if(typeb=="Intron")printf("Intron %s__%d_%d\n",chrom,$8,$9);if($2=="Reverse")dx=-1;else dx=1;if(typea=="Intron")printf("Intron %s__%d_%d\n",chrom,$7+dx,$8-dx);if(typec=="Intron")printf("Intron %s__%d_%d\n",chrom,$9+dx,$10-dx);printf("\n");}' group=$MAGIC  > tmp/X.$MAGIC/$chrom/f1.ace
 
echo "collate fasta"

cat   tmp/X.$MAGIC/$chrom/f1.txts  | gawk -F '\t' '/^#/{next}/^\/\//{next}{chrom=$1;gsub(/CHROMOSOME_/,"",chrom);u1=$8;u2=$9;typeb=$4;if(typeb == "Exon") next; xx="XI_" ; nam=xx group "_" chrom "__" u1 "_" u2 ; names[nam]++; n = names[nam] ;nam = nam "." n ; printf(">%s\n%s%s%s\n",nam, $12,$13,$14) ;}' group=$MAGIC  > tmp/X.$MAGIC/$chrom/f1.fasta

gzip  tmp/X.$MAGIC/$chrom/f1.ace
gzip  tmp/X.$MAGIC/$chrom/f1.fasta
gzip  tmp/X.$MAGIC/$chrom/f1.txts
touch tmp/X.$MAGIC/$chrom/f1.done

exit 0

