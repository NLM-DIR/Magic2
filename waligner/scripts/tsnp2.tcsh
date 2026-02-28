#!bin/tcsh -f

set phase=$1
set zone=$2
set run=$3

if ($phase == snp2a) goto snp2a
if ($phase == tsnp2a) goto tsnp2a
if ($phase == tsnp2b) goto tsnp2b
if ($phase == tsnp2c) goto tsnp2c
if ($phase == tsnp2g) goto tsnp2g

echo "bad phase in tsnp2.tcsh $1 $2 $3"
date
exit 1


snp2a:
tsnp2a:

if (1 && ! -e tmp/TSNP_DB/$zone/database) then
   echo "scripts/tsnp_DB.tcsh $Strategy $zone"
   scripts/tsnp_DB.tcsh $Strategy $zone
endif


echo -n "tsnp2a: start "
date


  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    if (-e tmp/TSNP/$run/$zone/tsnp1.2.0.deUno.tsf && ! -e tmp/TSNP/$run/$zone/tsnp2.2.deUno.tsf) then
      echo "cat  tmp/TSNP/$run/$zone/tsnp1.2.*.deUno.tsf | bin/tsf --merge --setSample $run > tmp/TSNP/$run/$zone/tsnp2.2.deUno.tsf"
            cat  tmp/TSNP/$run/$zone/tsnp1.2.*.deUno.tsf | bin/tsf --merge --setSample $run > tmp/TSNP/$run/$zone/tsnp2.2.deUno.tsf
    endif
    if (-e tmp/TSNP/$run/$zone/tsnp1.MB.0.deUno.tsf && ! -e tmp/TSNP/$run/$zone/tsnp2.MB.deUno.tsf) then
      cat  tmp/TSNP/$run/$zone/tsnp1.MB.*.deUno.tsf | bin/tsf --merge --setSample $run > tmp/TSNP/$run/$zone/tsnp2.MB.deUno.tsf
    endif    
  end

  echo "-R Map NC_045512.2  NC_045512\n" >  tmp/TSNP_DB/map.rename.ace

  set toto=tmp/TSNP_DB/$zone/$MAGIC.tsnp2a._r
    echo ' ' > $toto.ace
  echo "read-models" > $toto
  echo "query find variant BRS_counts" >> $toto
  echo "edit -D BRS_counts" >> $toto
  echo "pparse $toto.ace" >> $toto
  if (-e DanLi/DanLi.$zone.ace)  echo "pparse DanLi/DanLi.$zone.ace" >> $toto 
  echo "pparse MetaDB/$MAGIC/runs.ace" >> $toto

  echo "pparse MetaDB/$MAGIC/groups.ace" >> $toto
  echo "pparse MetaDB/$MAGIC/samples.ace" >> $toto

  set foundIn=Found_in_genome
  set mapIn=IntMap
  if ($NO_INTRON == 0 && $Strategy == RNA_seq)  then
    set foundIn=Found_in_mRNA
    set mapIn=mRNA
  endif


  echo '#' > $toto.M.tsf  
  echo '#' > $toto.MB.tsf  
  echo '#' > $toto.BRS.tsf  

  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    set minSnpFrequency2=$minSnpFrequency
    if (-e tmp/TSNP/$MAGIC.minFrequency.txt) then
      # reset minSnpFrequency to at least twice the percent error rate, the error_profile is given in error per million 
      set z=`cat tmp/TSNP/$MAGIC.minFrequency.txt | gawk '{if ($1 == run) z=($4 + 0)/5000 ; if (mf<z)mf=z;}END{print int(mf);}' mf=$minSnpFrequency run=$run`
      set minSnpFrequency2=$z
    endif

    if (-e tmp/TSNP/$run/$zone/tsnp2.2.deUno.tsf) then
      cat tmp/TSNP/$run/$zone/tsnp2.2.deUno.tsf | gawk -F '\t' '/^#/{next;}{print}' >> $toto.M.tsf  
    endif
    if (-e tmp/TSNP/$run/$zone/tsnp2.MB.deUno.tsf) then
      cat tmp/TSNP/$run/$zone/tsnp2.MB.deUno.tsf |  gawk -F '\t' '/^#/{next;}{print;}'  >> $toto.MB.tsf 
    endif    
    if (-e  tmp/SNP/$run/$MAGIC.$zone.u.snp_count.tsf.gz) then
      gunzip -c tmp/SNP/$run/$MAGIC.$zone.u.snp_count.tsf  >> $toto.BRS.tsf 
    endif

  end

# throw in the groups
foreach gr (`cat MetaDB/$MAGIC/g2r | cut -f 1`)
  continue
  echo > _x.$$
  foreach run (`cat MetaDB/$MAGIC/g2r | gawk -F '\t' '{if($1==g)print $2;}' g=$gr`)
    echo $2 >> _x.$$
  end
  cat _x.$$ ZZZZZ $toto.BRS.tsf | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]==1;next;}if (ok[$2]==1)print;}' > toto.$gr
end 

# extract a list of good snps:  MB is method MagicBlast
echo 'toto'  > $toto.list
cat $toto.M.tsf |  gawk -F '\t' '/^#/{next;}{v=$1;m=$7;c=$9;if(m>c)c=m;r=c-m;if(c>minC){f=100*m/c;if(f>=minF) print v;}}' minF=$minSnpFrequency minC=$minSnpCover >> $toto.list
cat $toto.MB.tsf |  gawk -F '\t' '/^#/{next;}{v=$1;m=$7;c=$9;if(m>c)c=m;r=c-m;if(c>minC){f=100*m/c;if(f>=minF) print v;}}' minF=$minSnpFrequency minC=$minSnpCover >> $toto.list

cat $toto.list | sort -u | wc
cat $toto.BRS.tsf |  gawk -F '\t' '/^#/{next;}{v=$1;c=$7;m=$8;if(c>minC){f=100*m/c;if(f>=minF) print v;}}' minF=$minSnpFrequency minC=$minSnpCover >> $toto.list
cat $toto.list | sort -u | wc

cat $toto.list | sort -u > $toto.sorted_list
\mv $toto.sorted_list $toto.list

# export an ace file
cat $toto.M.tsf | sort -k 1,1 -k 2,2n -k 4,4 > $toto.M_sorted.tsf 
cat $toto.MB.tsf | sort -k 1,1 -k 2,2n -k 4,4 > $toto.MB_sorted.tsf 
cat $toto.BRS.tsf | gawk '/^#/{next;}{print;}' |  sort -k 1,1 -k 2,2n  > $toto.BRS_sorted.tsf 
\rm $toto.M.tsf  $toto.MB.tsf  $toto.BRS.tsf 

cat $toto.list ZZZZZ $toto.M_sorted.tsf |  gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}{v=$1;if(ok[v]<1)next;run=$2;a1=$4;a2=$5;m=$7;c=$9;tag=$12;if(substr($12,1,6)=="Multi_")tag=$12" " $13" "$14;if(m>c)c=m;r=c-m;if(c>minC)f=100*m/c;else f=-10;if(f>=minF){if(v!=oldV){split(v,aa,":");seq=aa[1];printf("\nVariant %s\n%s\nParent_sequence %s\n%s %s %d %d\n%s\n",v,foundIn,seq,mapIn,seq,a1,a2,tag);oldV=v;}printf("MCounts %s %d %d %d Frequency %.2f\n",run,m,r,c,f);}}' minF=$minSnpFrequency minC=$minSnpCover foundIn=$foundIn mapIn=$mapIn > $toto.ace
cat $toto.list ZZZZZ $toto.MB_sorted.tsf |  gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){ok[$1]=1;next;}}{v=$1;if(ok[v]<1)next;run=$2;a1=$4;a2=$5;m=$7;c=$9;tag=$12;if(substr($12,1,6)=="Multi_")tag=$12" " $13" "$14;if(m>c)c=m;r=c-m;if(c>minC)f=100*m/c;else f=-10;if(f>=minF){if(v!=oldV){split(v,aa,":");seq=aa[1];printf("\nVariant %s\n%s\nParent_sequence %s\n%s %s %d %d\n%s\n",v,foundIn,seq,mapIn,seq,a1,a2,tag);oldV=v;}printf("MBCounts %s %d %d %d Frequency %.2f\n",run,m,r,c,f);}}' minF=$minSnpFrequency minC=$minSnpCover foundIn=$foundIn mapIn=$mapIn >> $toto.ace
cat $toto.BRS_sorted.tsf | gawk -F '\t' '{v=$1;run=$2;split($1,aa,":");seq=aa[1];a1=aa[2]+0;type=aa[3];cp=$4;mp=$5;wp=$6;cm=$7;mm=$8;wm=$9;c=cp+cm;m=mp+mm;w=wp+wm;if(type=="Del" || type == "Ins"){wp=cp;wm=cm;} if (c==0)c=1;f=100.0*m/c ;if (v!=oldV){printf("\nVariant %s\n%s\n",v,foundIn);oldV=v;};printf("BRS_counts %s %d %d %d %d %d %d %d Frequency %.2f\n",run,c,m,w,mp,wp,mm,wm,f);}END{printf("\n");}'  foundIn=$foundIn   >> $toto.ace
 

if (-e  tmp/TSNP_DB/$zone/av.parse.done) then
# rattrapage of new runs, say nanopore, in an existing variant database  
   if (! -e tmp/TSNP_DB/$zone/$MAGIC.$phase.done) then 
     echo "pparse $toto.ace" |  bin/tacembly tmp/TSNP_DB/$zone -noprompt
     touch tmp/TSNP_DB/$zone/$MAGIC.$phase.done
  endif
  goto done1
  exit 0
endif


  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    if (-e tmp/TSNP/$run/$zone/tsnp1.2.0.bridges.ace) then
      cat  tmp/TSNP/$run/$zone/tsnp1.2.*.bridges.ace | gawk '/^Observed_genomic_sequence/{next;}/^Reference_genomic_sequence/{next;}/^[fr]Counts/{next;}/n[su]Counts/{next;}{print}' >> $toto.ace
      cat  tmp/TSNP/$run/$zone/tsnp1.2.*.bridges.ace | gawk '/^Variant/{v=$2;}/^nsCounts/{m[v]+=$3;r[v]+=$4;c[v]+=$5;}END{for(v in c){f=100*m[v]/c[v];if (f>=-95)printf("Variant %s\n%s %s %d %d %d Frequency %.2f\n\n",v,cc,run,m[v],r[v],c[v],100.0*m[v]/c[v]);}}' run=$run cc=MCounts  >> $toto.ace
    endif
    if (-e tmp/TSNP/$run/$zone/tsnp1.MB.0.bridges.ace) then
      cat  tmp/TSNP/$run/$zone/tsnp1.MB.*.bridges.ace | gawk '/^Observed_genomic_sequence/{next;}/^Reference_genomic_sequence/{next;}/^[fr]Counts/{next;}/n[su]Counts/{next;}{print}' >> $toto.ace
      cat  tmp/TSNP/$run/$zone/tsnp1.MB.*.bridges.ace | gawk '/^Variant/{v=$2;}/^nsCounts/{m[v]+=$3;r[v]+=$4;c[v]+=$5;}END{for(v in c){f=100*m[v]/c[v];if (f>=-95)printf("Variant %s\n%s %s %d %d %d Frequency %.2f\n\n",v,cc,run,m[v],r[v],c[v],100.0*m[v]/c[v]);}}' run=$run cc=MBCounts  >> $toto.ace
    endif
  end

  if (-e tmp/METADATA/gtf.av.goodProduct.ace) then
    echo "pparse tmp/METADATA/gtf.av.goodProduct.ace" >> $toto
  endif
  if (-e tmp/METADATA/$MAGIC.av.captured_genes.ace) then
    echo "pparse tmp/METADATA/$MAGIC.av.captured_genes.ace" >> $toto
  endif
  echo "pparse MetaDB/$MAGIC/runs.ace" >> $toto
  echo "pparse MetaDB/$MAGIC/groups.ace" >> $toto
  echo "pparse MetaDB/$MAGIC/samples.ace" >> $toto
  echo "pparse tmp/SNP_ZONE/$zone.fasta.gz" >> $toto

  echo "pparse tmp/TSNP_DB/map.rename.ace" >> $toto
  echo 'save' >> $toto
  echo 'quit' >> $toto

  echo -n " tsnp2a: parse the tsnp "
date
  mv $toto.ace $toto.preace
  cat $toto.preace | gawk '/^-/{next}{print}' > $toto.ace
  \rm $toto.preace
  bin/tace tmp/TSNP_DB/$zone <  $toto

  if (! -e tmp/TSNP_DB/$zone/av.parse.doneZZZ) then
    tace tmp/TSNP_DB/$zone <<EOF
      find mrna
      spush
      follow gene
      sor
      undo
      follow product 
      sor
      spop
      list -a -f tmp/TSNP_DB/$zone/mrna.list
      quit
EOF
    tace tmp/TSNP_DB/Genes <<EOF
      key tmp/TSNP_DB/$zone/mrna.list
      spush
      follow dna
      sor
      spop
      show -a -f tmp/TSNP_DB/$zone/mrna.ace
      quit
EOF

    tace tmp/TSNP_DB/$zone <<EOF
      pparse tmp/TSNP_DB/$zone/mrna.ace
      find gene toto
      query find mrna ; ! DNA ; > Variant
      save
      quit
EOF

    touch tmp/TSNP_DB/$zone/av.parse.done
  endif




  echo -n " remap and translate "
date

done1:
  set remap2g=remap2genes
  if ($NO_INTRON == 0 && $Strategy == RNA_seq) then
    set remap2g=remap2genome
  endif
# if -o filename is not provided, tsnp directly edits the database
  echo "  bin/tsnp --db_$remap2g  tmp/METADATA/mrnaRemap.gz  --db tmp/TSNP_DB/$zone "
          bin/tsnp --db_$remap2g  tmp/METADATA/mrnaRemap.gz  --db tmp/TSNP_DB/$zone 
  echo "  bin/tsnp --db_translate --db tmp/TSNP_DB/$zone  -p $MAGIC"
          bin/tsnp --db_translate --db tmp/TSNP_DB/$zone  -p $MAGIC

#          bin/tsnp --db_report --db tmp/TSNP_DB/$zone  -p $MAGIC
#    bin/tsnp --db_snp_profile --db tmp/TSNP_DB/$zone  -p $MAGIC
#    bin/tsnp --db_group_count --db tmp/TSNP_DB/$zone  -p $MAGIC

# phase may be snp2a (BRS) or tsnp2a) tricoteur
touch tmp/TSNP_DB/$zone/$MAGIC.$phase.done
echo -n "tsnp2a: done "
date

goto phaseLoop

########################################################################################
## tsnp2b make words

tsnp2b:
echo -n "tsnp2b: start makewords  "
date

# bin/tsnp --db tmp/TSNP_DB/$zone -p $MAGIC --makeWords --zone $zone --gzo -o tmp/TSNP_DB/$zone/tsnp2b.$MAGIC --filter "MBcounts && substitution"
bin/tsnp --db tmp/TSNP_DB/$zone -p $MAGIC --makeWords --zone $zone -gzo -o tmp/TSNP_DB/$zone/tsnp2b.$MAGIC 
touch tmp/TSNP_DB/$zone/tsnp2b.done
echo -n "tsnp2b done "
date

goto phaseLoop

########################################################################################
## tsnp2c count words in run

tsnp2c:
echo -n "tsnp2c: count words in runs start "
date

  bin/tricoteur -count -wLn 31 -wordFile tmp/TSNP_DB/$zone/tsnp2b.$MAGIC.w31.gz -run $run -gzo -o tmp/TSNP/$run/$zone/tsnp2c.$MAGIC

echo -n "tsnp2c: count words in runs done "
date
goto phaseLoop

########################################################################################
## tsnp2g study the GGG

tsnp2g:
echo -n "tsnp2g: study the GGG"
date

    bin/tace tmp/TSNP_DB/$zone <<EOF
      read-models
      save
      quit     
EOF
  bin/tsnp --db_GGG --db tmp/TSNP_DB/$zone -o tmp/TSNP_DB/$zone/tsnp2g
echo -n "tsnp2g: study the GGGG"
date
goto phaseLoop

########################################################################################
## junk
if (0) then



    bin/tace tmp/TSNP_DB/$zone <<EOF
      pparse MetaDB/$MAGIC/runs.ace
      pparse tmp/SNP_ZONE/$zone.fasta.gz
      save
      quit     
EOF
    if (-d tmp/SNP_DB/$zone/database && ! -e tmp/TSNP_DB/$zone/tsnp2.$MAGIC.w31.gz) then
      bin/tsnp --db tmp/TSNP_DB/$zone -p $MAGIC --makeWords --zone $zone --gzo -o tmp/TSNP_DB/$zone/tsnp2.$MAGIC --filter "MBcounts && substitution"
    endif
 

tace tmp/TSNP_DB/$zone <<EOF
  Find Variant
  edit -D Reference_genomic_sequence
  edit -D Observed__genomic_sequence 
  parse tmp/TSNP_DB/zoneG.mutated_cov_May7_a/tsnp2._r.ace
  save
  quit
EOF

 ~/ace/bin.LINUX_4/tsnp -db_translate -db tmp/TSNP_DB/$zone -p $MAGIC  > tmp/TSNP_DB/$zone/tsnp2.translate.ace
 date
echo "pparse tmp/TSNP_DB/$zone/tsnp2.translate.ace" |  bin/tacembly tmp/TSNP_DB/$zone -noprompt
date
tace tmp/TSNP_DB/$zone <<EOF
  find variant
  select --title "Twist SNPS in Dan Modified genome"  -o RESULTS/SNV/Twist.tsf.txt v,r,s,w,f,ref,ob from v in @, r in v->MCounts, s in r[1], w in r[3], f in r[5],ref in v->Reference_genomic_sequence,ob in v->Observed__genomic_sequence TITLE v:SNP, r:Run,s:Support,w:Wiggle,f:Frequency,ref:Reference sequence,ob:Observed sequence
  quit
EOF
cat RESULTS/SNV/Twist.tsf.txt | grep Ins
  

endif

########################################################################################

phaseLoop:
  echo "phase $phase done"

########################################################################################
########################################################################################
