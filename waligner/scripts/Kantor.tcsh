#!bin/tcsh -f

set phase=$1
set chrom=$2
set db=$3

echo "scripts/Kantor.tcsh $*"

if ($phase == k0) goto phasek0
if ($phase == k1) goto phasek1
if ($phase == k2) goto phasek2
if ($phase == k9) goto phasek9

# if ($species == hs && $chrom != chr4) goto phaseLoop

############################################################
## Phase k0 : create the Kantor/$chrom database

phasek0:

# Create the Kantor database
  set KDB=Kantor/$chrom/KantorDB
  if (! -d Kantor/$chrom/done) mkdir Kantor/$chrom/done
  if (! -d Kantor/$chrom/done/creation) mkdir Kantor/$chrom/done/creation
echo bb $species $species_kantor
  if (! -d Kantor/$chrom/tmp/$species_kantor.data) mkdir Kantor/$chrom/tmp/$species_kantor.data
  if (! -d Kantor/$chrom/tmp/$species_kantor.data/$chrom) mkdir Kantor/$chrom/tmp/$species_kantor.data/$chrom
echo cc
  if (! -d $KDB) then
     echo "// hello" > touch Kantor/$chrom/done/creation/ace.creation
     mv Kantor/$chrom/done/* Kantor/$chrom/tmp
     mkdir $KDB
     pushd $KDB
       mkdir database
       ln -s ../../../metaData/wspec.aceview_web_site wspec
       mkdir database
       tacembly . <<EOF
y
        save
        quit
EOF
     popd
  endif
  touch $db/k0.done
  \rm $db/k0.start

goto phaseLoop

############################################################
## Phase k1 : recover the available kantor and parse them in the kantor database

# Hack, recover in the trace DB edited interactivelly the kantor precomputed in magic
set myZH=ZHabi_traces
tace $myZH <<EOF
  query find kantor product
  list -a -f $myZH/kantor.list
  quit
EOF

foreach chrom ($chromSetAll)
   tace Kantor/$chrom/KantorDB << EOF
      key $myZH/kantor.list
      show -a -f $myZH/kantor.$chrom.ace
      quit
EOF
end

cat $myZH/kantor.*.ace > $myZH/all.kantor.ace
tacembly $myZH <<EOF
  parse $myZH/all.kantor.ace
    query  > product ; >mrna 
    acembly
      cdna_kantor // kantorizes active mrna set
      quit
      // dna mmMrna.$chrom.dna
    save
  quit
EOF



phasek1:
echo -n "Start of phase k1 chrom=$chrom "
date

set KDB=Kantor/$chrom/KantorDB

# construct the list of existing kantor files, and parse them
  echo 'read-models' >  Kantor/$chrom/_r
 
  set ok=0
  echo "" >  Kantor/$chrom/k1._r
  if (-d Kantor/$chrom/tmp) then
     foreach ff (`ls Kantor/$chrom/tmp/*/ace.*`)
        echo "pparse $ff" >>  Kantor/$chrom/k1._r
        set ok=1
        # \rm -rf $dd
      end
  endif
  echo "save" >> Kantor/$chrom/k1._r
  echo "quit" >> Kantor/$chrom/k1._r


  # a hack to clean a current PSort bug 
  if (0) then
    foreach ff (`ls  Kantor/$chrom/psortZZZ/ace.*`)
      cat $ff | gawk '/^Psort Localization \"\"/{next;}{print}' > $ff.bis
      mv  $ff.bis $ff
    end
  endif

# parse the data in the kantor database

  if ($ok == 0) touch $db/k1.done
  if ($ok == 1  && -e $KDB/database/log.wrm && ! -e $KDB/database/lock.wrm) then
    if (-e $db/k2.done) \rm $db/k2.done
    if (-e $db/k9.done) \rm $db/k9.done


    tacembly $KDB <  Kantor/$chrom/k1._r
    if (! -d Kantor/$chrom/done) mkdir Kantor/$chrom/done
    mv Kantor/$chrom/tmp/*  Kantor/$chrom/done
  endif
  \rm $db/$phase.start
  touch $db/k1.done

goto phaseLoop

############################################################
############################################################
## Phase k2 : parse the data in the active database

#to change the kantor because of the model change

phasek2:

set KDB=Kantor/$chrom/KantorDB

  tacembly $db << EOF
    query find kantor product
    list -a -f $db/k2.list
    quit
EOF

  tacembly $KDB << EOF
    query find pfam IS Pox_polyA_pol
    kill
    key $db/k2.list
    show -a -f $db/k2.ace
    quit
EOF

echo "\n\n-R Method Pfam.v32 Pfam\n\n" >> $db/k2.ace

  tacembly $db << EOF
    query find pfam IS Pox_polyA_pol
    kill
    pparse $db/k2.ace
    query  > product ; >mrna 
    acembly
      cdna_kantor // kantorizes active mrna set
      quit
      // dna mmMrna.$chrom.dna
    save
    find model
    list -a -f $db/k2.done
    quit
EOF


 if (-e $db/k9.done) \rm $db/k9.done
 if (-e $db/$phase.start) \rm $db/$phase.start
 touch $db/k2.done

goto phaseLoop

############################################################
############################################################
## Phase k9 : megaRun this chromosome

#to change the kantor because of the model change

phasek9:

echo -n "Start of phase k9 megaRun chrom=$chrom "
date

setenv ici `pwd`
set KDB=Kantor/$chrom/KantorDB
echo k9aaa

echo k9bbb
    setenv megaRun ~/MEGA3/scripts/megaRun
    echo -n 'starting MEGA3/scripts/megaRun :'
    date

    setenv chrom $chrom
    pushd Kantor/$chrom/tmp
      # if (! -d tmp/$species_kantor.data) mkdir tmp/$species_kantor.data
       $megaRun $ici/$db  psort $species_kantor
      # $megaRun $ici/$db  acekog $species_kantor  
       $megaRun $ici/$db  pfam $species_kantor  
       $megaRun $ici/$db blastp $species_kantor  
      # $megaRun $ici/$db oligo $species_kantor   
      # $megaRun $ici/$db acekog_n $species_kantor  
      date
    popd

    echo -n "removing $db/k1.done,  pwd="
    pwd
    ls -ls $db/k1.done
    if (-e $db/k1.done) \rm $db/k1.done
    if (-e $db/k2.done) \rm $db/k2.done
    if (-e $db/$phase.start) \rm $db/$phase.start
    touch  $db/k9.done


echo -n 'End of phase k9'
date

goto phaseLoop

############################################################
############################################################
## Phase k9_server : megaRun this chromosome

#to change the kantor because of the model change

phasek9_server:

echo -n "Start of phase k9 megaRun chrom=$chrom "
date

  # setenv chrom is used by megaRun
  setenv chrom $chrom
  setenv ici `pwd`
  pushd $db
  set k=`cat wspec/acetcp_access.wrm | grep '0.0'0.255' | wc -l`
  if ($k < 1) echo "\n0.0.0.255 w" >>  wspec/acetcp_access.wrm

  setenv port 1235901
  set chrom2=`echo $chrom | sed -e 's/^chr//'`
    if ($chrom2 == "1") setenv port 1235901
    if ($chrom2 == "2") setenv port 1235002
    if ($chrom2 == "2L") setenv port 2235002
    if ($chrom2 == "2R") setenv port 3235002
    if ($chrom2 == "3") setenv port 1235003
    if ($chrom2 == "3L") setenv port 2235003
    if ($chrom2 == "3R") setenv port 3235003
    if ($chrom2 == "4") setenv port 1235004
    if ($chrom2 == "5") setenv port 1235005
    if ($chrom2 == "6") setenv port 1235006
    if ($chrom2 == "7") setenv port 1235007
    if ($chrom2 == "8") setenv port 1235008
    if ($chrom2 == "9") setenv port 1235009
    if ($chrom2 == "10") setenv port 1235010
    if ($chrom2 == "11") setenv port 1235011
    if ($chrom2 == "12") setenv port 1235012
    if ($chrom2 == "13") setenv port 1235013
    if ($chrom2 == "14") setenv port 1235014
    if ($chrom2 == "15") setenv port 1235015
    if ($chrom2 == "16") setenv port 1235016
    if ($chrom2 == "17") setenv port 1235017
    if ($chrom2 == "18") setenv port 1235018
    if ($chrom2 == "19") setenv port 1235019
    if ($chrom2 == "20") setenv port 1235020
    if ($chrom2 == "21") setenv port 1235021
    if ($chrom2 == "22") setenv port 1235022
    if ($chrom2 == "X") setenv port 1235023
    if ($chrom2 == "Y") setenv port 1235024
    if ($chrom2 == "Un") setenv port 1235025
    if ($chrom2 == "MT") setenv port 1235026
    if ($chrom2 == "Pltd") setenv port 1235027

    set iiStart=0

  iciStart: 
    echo "warning: at Try 0 we expect the error message connect: Connection refused"

    echo -n "Try $iiStart "  
    taceclient a:localhost:$port  <<EOF >! k9.test
      date
      quit
EOF
    @ iiStart = $iiStart + 1
    set startOk = `grep -c 'Active Objects' k9.test `
    echo " ok=$startOk"  
    if ($startOk > 0) goto iciStartDone 
    if ($iiStart > 1 && $iiStart < 60) then
      sleep 10
      goto iciStart
    endif


    echo -n "launching  taceserver $db $port : "
    date
    #    cat /etc/hosts | grep $HOST | gawk '{printf("%s w\n", $1)}' >>  wspec/acetcp_access.wrm

    taceserver . $port 3600:3600:0 &
    echo 'sleep 10'
    sleep 10
    goto iciStart

  iciStartDone:
    popd

    setenv megaRun ~/MEGA3/scripts/megaRun
    echo -n 'starting MEGA3/scripts/megaRun :'
    date
    touch  $db/k9.megaRun.start
    if (! -d Kantor) mkdir Kantor
    if (! -d Kantor/$chrom) mkdir Kantor/$chrom
    if (! -d Kantor/$chrom/tmp) mkdir Kantor/$chrom/tmp
    pushd Kantor/$chrom/tmp
      # if (! -d tmp/$species_kantor.data) mkdir tmp/$species_kantor.data
      $megaRun a:localhost:$port  psort $species_kantor
      # $megaRun a:localhost:$port acekog $species_kantor  
      # $megaRun a:localhost:$port   pfam $species_kantor  
      $megaRun a:localhost:$port blastp $species_kantor  
      # $megaRun a:localhost:$port oligo $species_kantor   
      # $megaRun a:localhost:$port acekog_n $species_kantor  
      date
    popd

    echo -n "shutting down taceserver $db port=$port : "
    date
    taceclient a:localhost:$port <<EOF > /dev/null
      shutdown now
      quit
EOF
    if (-e $db/k1.kantor_parse.done) \rm $db/k1.kantor_parse.done
    touch  $db/k9.megaRun.done
    \rm $db/k9.megaRun.start
  endif
  popd
end

echo -n 'End of phase k9_server'
date

goto phaseLoop

############################################################
phaseLoop:
 echo "phase $phase done"
 exit 0

