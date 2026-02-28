#!bin/tcsh -ef

set phase=$1
set run=$2
set chrom=$3
set uu=$4
set fr=$5
set limit=$6
set mask="$7"
set mainTarget=$8

# genome gene_and_new_exon new_and_known_exon known_exon RefSeq
echo "sponge_chrom.tcsh $phase $run $chrom $uu $fr $limit $mask"

if ($phase == stats) goto phaseStats
if ($phase == pressGene) goto phasePressGene

goto phaseLoop

phaseStats:

# prepare the frns wiggles if not available, now always done in wg2a for runs and wg2b for groups 
set ww1=toto$$
set ww2=toto$$
set ok=0
if (-e  tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$fr.$uu.BF.gz || -e  tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$fr.pp.BF.gz) then
  set ww1=tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$fr.$uu.BF.gz 
  set ww2=tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$fr.pp.BF.gz
  set ok=1
else if (-e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.$uu.BF.gz || -e  tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.pp.BF.gz) then
  set ww1=tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.$uu.BF.gz
  set ww2=tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.pp.BF.gz
  set ok=1
endif

if ($ok == 0) then
  echo "tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.[u|pp].BF.gz not found"
  exit 0
endif
if (! -e $ww1) then
  set ww=$ww2
else if (! -e $ww2) then
  set ww=$ww1
else
  set ww="$ww1,$ww2"
endif

echo "##Phase 4444444444444444444"
if ($fr == frns && $mainTarget != introns) then        
  set toto=tmp/SPONGE/$run/Total.$chrom.$uu.$fr.$limit
  if (! -e $toto.txt) then
    if (! -e $toto) then
      echo Total > $toto
      echo "... bin/geneelements -sponge $limit -spongeFile $mask  -sxxChromosome $chrom -wiggle $ww >  $toto"
                bin/geneelements -sponge $limit -spongeFile $mask  -sxxChromosome $chrom -wiggle $ww >  $toto
    endif

    grep Total $toto |   sed -e "s/^Total/$chrom\t$fr/" | gawk -F '\t' '{r=run;}/^#/{r="#Run\tChrom"}{printf("%s\t",r);print;}' run=$run > $toto.txt
    \rm $toto
  endif
endif


goto phaseLoop

phasePressGene:

echo "##Phase pressGene v4   limit=$limit fr=$fr"

set ww=NULL
#case f,r
if ($fr == f || $fr == r) then
  echo tt
if (-e tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr.BF.gz) set ww=tmp/WIGGLERUN/$run/$chrom/R.chrom.$uu.$fr.BF.gz
  if (-e tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$uu.$fr.BF.gz) set ww=tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$uu.$fr.BF.gz
endif
#case frns
if ($fr == frns) then
  if (-e tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.$uu.BF.gz) set tmp/WIGGLERUN/$run/$chrom/R.chrom.$fr.$uu.BF.gz 
  if (-e tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$fr.$uu.BF.gz) set tmp/WIGGLEGROUP/$run/$chrom/R.chrom.$fr.$uu.BF.gz 
  set fr=ns
endif
echo $ww
if (! -e $ww) goto phaseLoop



set long=`cat MetaDB/$MAGIC/RunNanoporeList  MetaDB/$MAGIC/RunPacBioList |gawk '{if($1==run)ok=1;}END{print ok+0;}' run=$run`
set long=1
if ($long == 1) then
  foreach mainTarget ($Etargets)
    foreach GM (gene mrna)
      set geneMask=tmp/METADATA/$mainTarget.$fr.$GM.sponge
      ls -ls $geneMask
      if (-e $geneMask && $limit == 1) then
        set toto=tmp/SPONGE/$run/$mainTarget.$chrom.$uu.$fr.$limit
        if (-e $toto)  \rm $toto
        set toto=tmp/SPONGE/$run/$mainTarget.$GM.$chrom.$uu.$fr.$limit
        if (-e $toto)  \rm $toto
        set toto=tmp/SPONGE/$run/$mainTarget.$GM.v4.$chrom.$uu.$fr.$limit
        set split=""
        if (-e tmp/METADATA/$mainTarget.split_mrnas.gz)  then
          set split="-split_mRNAs  tmp/METADATA/$mainTarget.split_mrnas.gz"
        endif
        ls -ls tmp/METADATA/$mainTarget.split_mrnas.gz
        if (! -e $toto) then
          echo "bin/geneelements -sponge $limit $split -spongeFile $geneMask  -sxxChromosome $chrom -wiggle $ww >  $toto"
                bin/geneelements -sponge $limit $split -spongeFile $geneMask  -sxxChromosome $chrom -wiggle $ww >  $toto
        endif
      endif
    end
  end
endif
goto phaseLoop

phaseLoop:
 echo done


