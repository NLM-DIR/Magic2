#!bin/tcsh -f

set chrom=$1
set phases="$2"

# This loop ends in phaseLoop
foreach phase ($phases)

echo "phase=$phase  project=$MAGIC"

if (0 && $chrom != 20) exit 0
if ($phase == chromDB) goto chromDB
if ($phase == parseGenome) goto parseGenome
if ($phase == parseTargetIntrons) goto parseTargetIntrons
if ($phase == setFeet) goto setFeet
if ($phase == deUno) goto deUno
if ($phase == deMrna) goto deMrna
if ($phase == setDA) goto setDA
if ($phase == setDAsupport) goto setDAsupport
if ($phase == setSponge) goto setSponge
if ($phase == setGroups) goto setGroups
if ($phase == intronCounts) goto intronCounts
if ($phase == altIntrons) goto altIntrons

if ($phase == deDuo) goto deDuo
if ($phase == collate) goto collate
if ($phase == donorAcceptor) goto donorAcceptor
if ($phase == diffAcceptor) goto diffAcceptor
if ($phase == capture) goto capture
if ($phase == intronComfirmation) goto  intronComfirmation

echo "unknown phase $phase"
exit 1
goto phaseLoop

###############################################################
chromDB:
echo -n "I1 phase $phase start :"
date

if (! -d tmp/INTRON_DB)   mkdir tmp/INTRON_DB
if (! -d tmp/INTRON_DB/$chrom)   mkdir tmp/INTRON_DB/$chrom
if (! -d tmp/INTRON_DB/$chrom/database) then
    pushd tmp/INTRON_DB/$chrom
    mkdir database
    ln -s ../../../metaData/wspec.aceview_web_site wspec
    ../../../bin/tace . <<EOF
y
EOF
    popd
endif
goto phaseLoop

###############################################################
### parse the genome 
parseGenome:
echo -n "I1 phase $phase start :"
date

if (! -e tmp/INTRON_DB/$chrom/I1.parse_genome.done) then
  pushd tmp/INTRON_DB/$chrom
    zcat ../../../TARGET/CHROMS/$species.chrom_$chrom.fasta.gz | gawk '/^>/{split($1,aa,"|");print aa[1] ;next;}{print}' | gzip > chrom_$chrom.fasta.gz
    ../../../bin/tace . <<EOF
      parse chrom_$chrom.fasta.gz
      find sequence $chrom
      edit genomic
      parse ../../../TARGET/Targets/$species.mito.fasta.gz
      save
      quit
EOF
    if (-e I1.intron_feet.done) \rm I1.intron_feet.done
    touch I1.parse_genome.done
  popd
endif

goto phaseLoop

###############################################################
### parse the gene, mrna, intron relations  : intron->supports

parseTargetIntrons:
echo -n "I1 phase $phase start :"
date

# parse in INTRON_DB/$chrom the introns/genes/mRNA/chromosomes

if (! -e tmp/INTRON_DB/$chrom/I1.parse_genes.done) then  
  zcat tmp/METADATA/gtf.*.introns.gz | gawk -F '\t' '{type=substr($1,4);c=$6;split (c,cc,"|") ;c=cc[1];if(c != chrom)next;i1=$7+0;i2=$8+0;ln=i2-i1;;if(ln<0)ln=-ln; ln=ln+1;printf("Intron %s__%d_%d\nIntMap %s %d %d\nLength %d\n%s\nGene \"%s\"\nIn_mRNA \"%s\" %d %d\n\n",c,i1,i2,c,i1,i2,ln,type,$2,$3,$4,$5);}' chrom=$chrom | gzip >  tmp/INTRON_DB/$chrom/I1.TargetIntrons.ace.gz
  pushd tmp/INTRON_DB/$chrom
    ../../../bin/tace . <<EOF
      read-models
      parse I1.TargetIntrons.ace.gz
      query find intron av
      edit AceView
      save
      quit
EOF
    touch I1.parseTargetIntrons.done
    if (-e I1.setFeet.done) \rm I1.setFeet.done
  popd
endif

goto phaseLoop

########################################################################
### parse the genomic de_uno counts

deUno:
echo -n "I2 phase $phase start :"
date

set ok=0
set ff=tmp/INTRON_DB/$chrom/$MAGIC.I2.deUno
echo ' ' > $ff

foreach run1 (`cat MetaDB/$MAGIC/RunList`) 
  set run2=$run1
  set run2=`cat MetaDB/$MAGIC/r2sublib | gawk -F '\t' '{if($2==run1){r=$1;last;}}END{if(r)print r;else print run1;}' run1=$run1`
  if (-e tmp/INTRON_DB/$chrom/I2.$run2.deUno.gz) continue
  set ok=1
  # echo "... deUno $run1 $run2"
  if ($USEMAGICBLAST == 1) then
    zcat tmp/MAGICBLAST/$run1/f2.*.introns.tsf.gz | bin/tsf -s $run2 --merge | gawk -F '\t' '{split($1,aa,"__");c=aa[1];if (c != chrom) next;plit(aa[2],bb,"_");a2=bb[1];b1=bb[2];if(a2<b1){strand="Forward";a2--;a1=a2-12;b1++;b2=b1+10;ln=b1-a2-1;}else{strand="Reverse";a2++;a1=a2+10;b1--;b2=b1-10;ln=a2-b1-1;}nr=$4;feet=$7;if(nr>0)printf("INTRON\t%s\t%s\t%09d\t%09d\t%s\t%09d\t%09d\t%s\t%d\t%d\t%d\n",strand,chrom,a1,a2,chrom,b1,b2,feet,ln,nr,nperfect);}' chrom=$chom  | sort | gzip > tmp/INTRON_DB/$chrom/I2.$run2.deUno.gz
  else
    gunzip -c tmp/PHITS_genome/$run1/*.introns.gz  tmp/PHITS_mito/$run1/*.introns.gz  | gawk -F '\t' '/^INTRON/{split($3,cc,"|");c=cc[1];split($6,cc,"|");c2=cc[1];if (c != chrom) next;if (c2 != chrom) next;i1=$5+0;i2=$7+0;if ($2=="Forward"){i1++;i2--;}else {i1--;i2++;}nn=$11; printf("%s__%d_%d\t%d\t%s\n",c,i1,i2,nn,$9);}' chrom=$chrom | gawk -F '\t' '{i=$1;n[i]+=$2;t[i]=$3;}END{for (i in n)printf("%s\t%s\t%d\t%s\n",i,run,n[i],t[i]);}' run=$run2  | sort | gzip > tmp/INTRON_DB/$chrom/I2.$run2.deUno.gz
  endif

  zcat tmp/INTRON_DB/$chrom/I2.$run2.deUno.gz >> $ff
end

if ($ok == 1) then
  date
  cat $ff | sort  > $ff.sorted
  date
  echo "bin/altintrons --deMrna $ff.sorted --db tmp/INTRON_DB/$chrom -p $MAGIC"
        bin/altintrons --deMrna $ff.sorted --db tmp/INTRON_DB/$chrom -p $MAGIC
  date
  if (-e  tmp/INTRON_DB/$chrom/I1.setFeet.done) \rm  tmp/INTRON_DB/$chrom/I1.setFeet.done
  \rm $ff $ff.ace
endif

goto phaseLoop

########################################################################
### parse the intron counts from mRNA alignments

deMrna:
echo -n "I2 phase $phase start :"
date

set ok=0
set f0=tmp/INTRON_DB/$chrom/I2.deMrna
if (-e $f0) \rm $f0
foreach run (`cat MetaDB/$MAGIC/RunsList`) 
  set f1=tmp/INTRON_DB/$chrom/I2.$run.deMrna
  if (-e $f1.gz) then
    set ok=`ls -ls $f1.gz | gawk '{print $6}'`
    if ($ok > 100) continue
  endif
  set ff=tmp/INTRONRUNS/$run/$run.u.intronSupport.counts.gz
  if (! -e $ff) continue 
  set ok=1
  zcat $ff | gawk -F '\t' '{split($1,cc,"__"); c = cc[1] ; if(c != chrom) next; print;}' chrom=$chrom > $f1
  cat $f1 >> $f0
  gzip -f $f1
end
date

if ($ok == 1) then
  cat $f0 | sort > $f0.sorted
  date
  echo "bin/altintrons --deMrna $f0.sorted --db tmp/INTRON_DB/$chrom -p $MAGIC"
        bin/altintrons --deMrna $f0.sorted --db tmp/INTRON_DB/$chrom -p $MAGIC
  date
  if (-e  tmp/INTRON_DB/$chrom/I1.setFeet.done) \rm  tmp/INTRON_DB/$chrom/I1.setFeet.done
  if (-e  tmp/INTRON_DB/$chrom/$MAGIC.I2.setDA.done)   \rm  tmp/INTRON_DB/$chrom/$MAGIC.I2.setDA.done
endif

goto phaseLoop

########################################################################
### clean up
cleanUp:

\rm _killIntrons
foreach run (`cat MetaDB/$MAGIC/RunsList`) 
  set f1=tmp/INTRON_DB/$chrom/I2.$run.deMrna
  if (-e $f1.gz) then
    set ok=`ls -ls $f1.gz | gawk '{print $6}'`
    if ($ok < 100) then
      echo "\\rm $f1.gz" >> _killIntrons
    endif
  endif
  set f1=tmp/INTRONRUNS/$run/$run.u.intronSupport.counts
  if (-e $f1.gz) then
    set ok=`ls -ls $f1.gz | gawk '{print $6}'`
    if ($ok < 100) then
      echo "\\rm -rf tmp/INTRONLANES/$run tmp/INTRONRUNS/$run"  >> _killIntrons
    endif
  endif
end
wc  _killIntrons

goto phaseLoop

########################################################################
### check we have all intron length feet intmap

setFeet:
echo -n "I1 phase $phase start :"
date

if (-d tmp/INTRON_DB/$chrom/database && ! -e tmp/INTRON_DB/$chrom/I1.setFeet.done) then
  echo "bin/altintrons --setFeet --db tmp/INTRON_DB/$chrom -p $MAGIC"
        bin/altintrons --setFeet --db tmp/INTRON_DB/$chrom -p $MAGIC
  touch tmp/INTRON_DB/$chrom/I1.setFeet.done
endif

goto phaseLoop

########################################################################
### create the donor/acceptors of all introns

setDA:
echo -n "I1 phase $phase start :"
date

if (-d tmp/INTRON_DB/$chrom/database && ! -e tmp/INTRON_DB/$chrom/$MAGIC.I2.setDA.done) then
  echo "bin/altintrons --setDA --db tmp/INTRON_DB/$chrom -p $MAGIC"
        bin/altintrons --setDA --db tmp/INTRON_DB/$chrom -p $MAGIC
  touch tmp/INTRON_DB/$chrom/$MAGIC.I2.setDA.done
endif

goto phaseLoop

########################################################################
### create the donor/acceptors of all introns

setDAsupport:
echo -n "I3 phase $phase start :"
date

if (-d tmp/INTRON_DB/$chrom/database && ! -e tmp/INTRON_DB/$chrom/$MAGIC.I3.setDAsupport.done) then
    echo "bin/altintrons --setDAsupport --db tmp/INTRON_DB/$chrom -p $MAGIC"
          bin/altintrons --setDAsupport --db tmp/INTRON_DB/$chrom -p $MAGIC
  touch tmp/INTRON_DB/$chrom/$MAGIC.I3.setDAsupport.done
endif

goto phaseLoop

########################################################################
### create the donor/acceptors sponge 

setSponge:
echo -n "I4 phase $phase start :"
date

if (-d tmp/INTRON_DB/$chrom/database && -e tmp/INTRON_DB/$chrom/$MAGIC.I3.setDAsupport.done) then
   bin/tacembly tmp/INTRON_DB/$chrom  <<EOF
     read-models
     pparse MetaDB/$MAGIC/runs.ace
     save
     quit
EOF
   echo "bin/altintrons --setSponge --db tmp/INTRON_DB/$chrom --setFeet --chrom $chrom -p $MAGIC"
         bin/altintrons --setSponge --db tmp/INTRON_DB/$chrom --setFeet --chrom $chrom -p $MAGIC
  touch tmp/INTRON_DB/$chrom/$MAGIC.I4.setSponge.done
endif

goto phaseLoop

########################################################################
### create the intron/donor/acceptors group counts

setGroups:
echo -n "I5 phase $phase start :"
date

if (-d tmp/INTRON_DB/$chrom/database &&  -e tmp/INTRON_DB/$chrom/$MAGIC.I4.setSponge.done) then
   bin/tacembly tmp/INTRON_DB/$chrom  <<EOF
     read-models
     pparse MetaDB/$MAGIC/runs.ace
     pparse MetaDB/$MAGIC/groups.ace
     save
     quit
EOF
   echo "bin/altintrons --setGroups --db tmp/INTRON_DB/$chrom  --chrom $chrom -p $MAGIC"
         bin/altintrons --setGroups --db tmp/INTRON_DB/$chrom  --chrom $chrom -p $MAGIC
  touch tmp/INTRON_DB/$chrom/$MAGIC.I5.setGroups.done
endif

goto phaseLoop

########################################################################
## 
# count introns export tsf

intronCounts:
echo -n "I6 phase $phase start :"
date

if (-d tmp/INTRON_DB/$chrom/database &&  -e tmp/INTRON_DB/$chrom/$MAGIC.I5.setGroups.done) then
   echo "bin/altintrons -db tmp/INTRON_DB/$chrom --counts -p $MAGIC -o tmp/INTRON_DB/$chrom/$MAGIC.I6"
         bin/altintrons -db tmp/INTRON_DB/$chrom --counts -p $MAGIC -o tmp/INTRON_DB/$chrom/$MAGIC.I6
  touch tmp/INTRON_DB/$chrom/$MAGIC.I6.counts.done
endif

goto phaseLoop

########################################################################
########################################################################
## 
# altIntrons

altIntrons:
echo -n "I7 phase $phase start :"
date

   echo "pparse MetaDB/$MAGIC/runs.ace" | bin/tace tmp/INTRON_DB/$chrom -no_prompt
   echo "bin/altintrons -db tmp/INTRON_DB/$chrom -p $MAGIC -o tmp/INTRON_DB/$chrom/$MAGIC.I7.altIntrons"
         bin/altintrons -db tmp/INTRON_DB/$chrom -p $MAGIC -o tmp/INTRON_DB/$chrom/$MAGIC.I7.altIntrons
goto phaseLoop

########################################################################
## deDuoCumul

deDuo:
echo -n "I1 phase $phase start :"
date
set toto = tmp/introns/$MAGIC.I1.de_duo.$chrom

if (! -e tmp/INTRON_DB/$chrom/$MAGIC.I1.deDuo.done) then

  # hack, because of the schema, the supports are set in intron->de_duo rather than de_uno in the final ace file
  set toto = tmp/INTRON_DB/$chrom/$MAGIC.I1.de_duo
  echo ' ' > $toto.1
  echo ' ' > $toto.txt
  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    set ff=tmp/OR/$run/d4.de_uno.txt.gz
    if (! -e $ff) set ff=tmp/OR/$run/I1.de_uno.txt.gz
    if (-e  $ff) then
      set minS=`cat  tmp/introns/$MAGIC.I1.minS | gawk -F '\t' 'BEGIN{n=1;}{if($1==run)n=$2+0;}END{print n}' run=$run`
      zcat $ff | gawk -F '\t' "/^$chrom\t/"'{a1=$2+0 ; a2=$3+0 ; ii=$1 "__" a1 "_" a2 ; if ($4>=minS) printf ( "%s\t%s\t%s\t%s\t%d\n",ii,$5,"any",run,$4 ) ; }' run=$run minS=$minS chrom=$chrom >> $toto.txt
    endif
  end

  cat $toto.txt | sort -V | gzip > $toto.txt.gz
  \rm  $toto.txt  $toto.1 $toto.2.gz
  # hack, because of the schema, the supports are set in intron->de_duo rather than de_uno in the final ace file
  zcat $toto.txt.gz  | gawk -F '\t' '{ii=$1; if($5+0<1)next;if(ii!=old){if (n>0)printf("RNA_seq %d\n",n);n=0;printf("\nIntron \"%s\"\n",ii);}old=ii;printf("de_duo %s %d\n",$4,$5);n+=$5;}END{if (n>0)printf("RNA_seq %d\n",n);printf("\n");}' | gzip > tmp/INTRON_DB/$chrom/$MAGIC.I1.de_duo.ace.gz

  if (-e tmp/INTRON_DB/$chrom/$MAGIC.I1.collate.done) \rm tmp/INTRON_DB/$chrom/$MAGIC.I1.collate.done
  touch tmp/INTRON_DB/$chrom/$MAGIC.I1.deDuo.done
endif
ls -ls tmp/INTRON_DB/$chrom/$MAGIC.I1.de_duo.ace.gz
goto phaseLoop

#############################################################################
## 
# check we have all intron length feet intmap
donorAcceptor:
echo -n "I1 phase $phase start :"
date

if (! -e tmp/INTRON_DB/$chrom/$MAGIC.I1.collate.done) goto phaseLoop



goto phaseLoop

###############################################################
## 
goto phaseLoop

###############################################################
### CAPTURE

capture:
  if ($?CAPTURES && -e tmp/METADATA/$MAGIC.av.captured_genes.ace && ! -e tmp/INTRON_DB/$chrom/$MAGIC.I1.capture.done) then 
    pushd tmp/INTRON_DB/$chrom
      ../../../bin/tace . << EOF
        read-models
        find gene
        spush
        parse  ../../../tmp/METADATA/$MAGIC.av.captured_genes.ace
        find gene
        sxor
        spop
        kill
        save
        bql -o  $MAGIC.I1.captured_introns.txt select ii,c from g in ?Gene, c in g->capture where c, ii in g->Intron  
        quit
EOF

     cat $MAGIC.I1.captured_introns.txt  | gawk -F '\t' '{if($1 != old)printf ("\nIntron %s\n",$1);old=$1;printf("Capture %s\n",$2);}END{printf("\n");}' > $MAGIC.I1.captured_introns.ace
     ../../../bin/tace .  << EOF
       parse $MAGIC.I1.captured_introns.ace
       save
       quit
EOF

      touch tmp/INTRON_DB/$chrom/$MAGIC.I1.capture.done
    popd  
   endif
  touch tmp/INTRON_DB/$chrom/$MAGIC.I1.capture.done
goto phaseLoop


###############################################################
###############################################################
## phaseLoop

phaseLoop:
  echo -n "I1 phase $phase done :"
  date
end
  echo I1.intronDB.tcsh $phases  done
  exit 0


###############################################################
###############################################################
