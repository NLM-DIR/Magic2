#!bin/tcsh -f

set chrom=$1

if (! -d tmp/X.$MAGIC) mkdir tmp/X.$MAGIC
if (! -d tmp/X.$MAGIC) mkdir tmp/X.$MAGIC
if (! -d tmp/X.$MAGIC/$chrom) mkdir tmp/X.$MAGIC/$chrom

set sMin=3
set sUno=""
set iDuo="tmp/X.$MAGIC/$chrom/f1.XI"
set sDuo="-sxxNewIntronsFileName  $iDuo "

# export intron with at least 1 percent forking at donor and acceptor sites, and the donor and acceptor are at least 10 percent of the reatained intron 
if (! -e $iDuo) then
  echo "tace INTRON_DB/$chrom "
  bin/tace tmp/INTRON_DB/$chrom <<EOF
    parse MetaDB/$MAGIC/runs.ace
    parse MetaDB/$MAGIC/groups.ace
    save
    date
    find intron
    select -o $iDuo ii, g, s, av, t, t2, c, c1, c2, dnaD, dnaA from ii in @ where ! ii#is_echo, g in ii->group where g#Intron && g->Intron[1] == $MAGIC, s in g[1] where s >= $sMin, t in ii->type, t2 in ii->Other, av in ii->stype[1], d in ii->D, a in ii->A, dg in d->group where dg == g, ag in a->group where ag == g, d1 in dg[1], d2 in dg[2], d3 in dg[3], a1 in ag[1], a2 in ag[2], a3 in ag[3] where 10*d3 >= d2 && 10*a3 >= a2 && 100 * s > d3 && 100 *s > a3, dnaD in d->motifs, dnaA in a->motifs, c in ii->IntMap, c1 in c[1], c2 in c[2]
    date
EOF
endif

cat $iDuo |  gawk -F '\t' '{ii=$1;g=$2;s=$3;av=$4;t=$5;t2=$6;c=$7;c1=$8+0;c2=$9+0;if(c1<c2){c1-=35;c2+=35;}else{c1+=35;c2-=35}printf("Sequence XI_%s_%s\n",g,ii);printf("cDNA_clone XI_%s_%s\nIs_read\nForward\nComposite %d\n",g,ii,s);if(av!="NULL")print av;if(t2!="NULL")printf("Other %s\n",t2);else print t;printf("Intron %s\nIntMap %s %d %d\n\n", ii, c,c1,c2);}'  | gzip > $iDuo.ace.gz

cat $iDuo |  gawk -F '\t' '{ii=$1;g=$2;dD=$10;dA=$11;i=index(dD,"--");sD=substr(dD,i-35,35);i=index(dA,"--");sA=substr(dA,i+2,35);printf(">XI_%s_%s\n%s%s\n",g,ii,sD,sA);}' |  gzip > $iDuo.fasta.gz
touch tmp/X.$MAGIC/$chrom/f1.done

exit 0

##########################################################################
##########################################################################



