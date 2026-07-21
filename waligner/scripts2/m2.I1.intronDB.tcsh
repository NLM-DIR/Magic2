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

if (! -e tmp/INTRON_DB/$chrom/I1.parse_introns.done) then
  tace GeneDB <<EOF
    query find intron IntMap == $chrom
    show -a -f  tmp/INTRON_DB/$chrom/I1.TargetIntrons.ace
    quit
EOF
  gzip tmp/INTRON_DB/$chrom/I1.TargetIntrons.ace
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

if (-e tmp/INTRON_DB/$chrom/I1.deUno.done)  goto phaseLoop

set ok=1
set ff=tmp/INTRON_DB/$chrom/$MAGIC.I2.deUno
echo ' ' > $ff

foreach run (`cat MetaDB/$MAGIC/RunList`) 
  cat tmp/SA/$run/introns.tsf | gawk -F '\t' '/^#/{next;}{split($1,aa,"_");if(aa[1]!=chrom)next;printf("%s\t%s\t%d\t%s\n",$1,$2,$5+$5,$6);}' chrom=$chrom >> $ff
end

if ($ok == 1) then
  date
  cat $ff | sort  > $ff.sorted
  date
  echo "bin/altintrons --deMrna $ff.sorted --db tmp/INTRON_DB/$chrom -p $MAGIC"
        bin/altintrons --deMrna $ff.sorted --db tmp/INTRON_DB/$chrom -p $MAGIC
  date
  if (-e  tmp/INTRON_DB/$chrom/I1.setFeet.done) \rm  tmp/INTRON_DB/$chrom/I1.setFeet.done
  if (-e  tmp/INTRON_DB/$chrom/I1.setFeet.done) \rm  tmp/INTRON_DB/$chrom/I1.setDA.done
  \rm $ff
  touch tmp/INTRON_DB/$chrom/I1.deUno.done
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
     pparse MetaDB/$MAGIC/groups.ace
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
