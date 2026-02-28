#!/bin/tcsh -f

# 2015_10_02
# Start from SraRunInfo.csv
# Using the SRA interface, select relevant runs
# use menu save-as file -> format RunInfo
#    the data should be set in excell to the format 12-Mar-15 (for 2015)
#    open the file in Excell on the PC, save as txt tab delimited, 
#    this transforms the format of the dates (which in csv contain comas)
#                  file -> format XML

# please edit the configuration
set phaseSet="$*"
if ("$phaseSet" == "" || "$phaseSet" == "help"  || "$phaseSet" == "-h" || "$phaseSet" == "--h" || "$phaseSet" == "-help" || "$phaseSet" == "--help") goto usage

set phase=`echo $phaseSet | gawk '{print $1}'`

if ($species == corona) then
  set date=2020_05_02
  set ff=~/aaa6/zoo/human/CRN_Mason2/SRX_CORONA/SRA_covid-19_may2_SraRunInfo.txt 
endif

if ($species == Ecoli) then
  set date=2015_10_03
  set ff=/home/mieg/ACEVIEWHELP/Ecoli_DATA/2015/1734E.coliK12_Oct7_2015_SRApublic_SraRunInfo.txt
  set date=2017_10_05
  set ff=/home/mieg/ACEVIEWHELP/Ecoli_DATA/2017/Ecoli_RNA_SRA_3449runs_8.2Tb.RunInfo.limit_to_taxon_511145.748_runs.2017_10_05.txt
  set date=2025_01_25
  set ff=/home/mieg/ACEVIEWHELP/Ecoli_DATA/2025/Ecoli_RNA_SRA_5852runs_12.8Tb.RunInfo.limit_to_taxon_511145.748_runs.2025_01_25.txt
endif

if ($species == droso) then
  set date=2016_04_19
  set ff=/home/mieg/ACEVIEWHELP/Droso/DATA/Droso_allRNA.2016_04_18.SraRunInfo.txt
  set date=2016_10_12 
  set ff=/home/mieg/ACEVIEWHELP/Droso/DATA/SraRunInfo_projectsMissingFromSRXdb_oct12_2016_andInFBfile.txt
  set date=2016_10_14
  set ff=/home/mieg/ACEVIEWHELP/Droso/DATA/2016Oct15_Drosophila_melanogaster_Organism_NOT_biomol_dna_Properties_13863_SraRunInfo.txt
  set ff=/home/mieg/ACEVIEWHELP/Droso_DATA/DATA/DrosophilaMelano_2017Mar20_4640RNA_since2016Apr10_PreviousWasApr19_SraRunInfo.txt
  set date=2025_01_21
  set ff=/home/mieg/ACEVIEWHELP/Droso_DATA/20250121_DrosoMelano_SraRunInfo.txt
endif

if ($species == hs) then
  set date=2016_10_14
  set ff='/home/mieg/ACEVIEWHELP/Human_DATA/*SraRunInfo.txt' 
  set date=2017_06_21  
  set ff='/home/mieg/ACEVIEWHELP/Human_DATA/194552_Human_RNA_exp_June22_2017_SraRunInfo.txt'
  set date=2017_08_31
  set ff='/home/mieg/ACEVIEWHELP/Human_DATA/194552_Human_RNA_exp_June22_2017_SraRunInfo.txt'
  set ff='/home/mieg/ACEVIEWHELP/Human_DATA/NA12878_GIAB_inSRA_Sept1_2017_102genome_9exomes_SraRunInfo.txt'
  set date=2018_08_13
  set ff='/home/mieg/ACEVIEWHELP/Human_DATA/GIAB_NA12878_runForMagicBlastSraRunInfo.txt' 
  set date=2018_10_12
  set ff='/home/mieg/ACEVIEWHELP/Human_DATA/20180916_RNAseq_public_human_16sept2018_SraRunInfo.txt'
  set ff='/home/mieg/ACEVIEWHELP/Human_DATA/20180916_1881_runs_longerThan596bases_.txt'
  set date=2018_10_29
  set ff='/home/mieg/ACEVIEWHELP/Human_DATA/2018_ANTE_SraRunInfo_165runs.txt'
  set date=2019_01_25
  set ff='/home/mieg/ACEVIEWHELP/Human_DATA/93Nanopore_Jan15_2019_HumanRNApublic_SraRunInfo.txt'
  set ff='/home/mieg/ACEVIEWHELP/Human_DATA/3Nanopore_Jan15_2019_HumanRNApublic_SraRunInfo.txt'
  set date=2019_11_27
  set ff='/home/mieg/AW/Human_DATA/128_exRNA_Sheng_Nov2019_SraRunInfo.txt'
  set date=2019_12_10
  set ff='/home/mieg/AW/Human_DATA/20191128_461_exRNA_fromSRA_11SRP_homosaporgn_exRNA_public_RNA_SraRunInfo.txt'
  set date=2019_12_14
  set ff='/home/mieg/AW/Human_DATA/Nature_Biotech_paper_Galas_148_SraRunInfo.txt'
  set date=2019_12_16
  set ff='/home/mieg/AW/Human_DATA/MAQC_SEQC2_data_from_SRA_3projects_SraRunInfo.txt'
  set date=2021_02_16
  set ff='/home/mieg/AW/Human_DATA/ToReferee//2021Feb_Ghimire_review_4ScReports_SraRunInfo.txt'
  set date=2021_09_20
  set ff='/home/mieg/AW/Human_DATA/BBS9_forKans_SraRunInfo.txt'
  set date=2023_04_17
  set ff='/home/mieg/AW/Human_DATA/CircularRNA_in_Frailty_SraRunInfo.txt'  
  set date=2023_06_18
  set ff='/home/mieg/AW/Human_DATA/20230618_RNA_Hydatidiform_mole_SraRunInfo.txt'
  set date=2023_08_13
  set ff='/home/mieg/AW/Human_DATA/20230813_453_Human_mole_RNAseq_SraRunInfo.txt'
  set date=2023_08_13
  #set ff='/home/mieg/AW/Human_DATA/20230813_840_Human_mole_DNA-seq_SraRunInfo.txt'
  set date=2025_05_13
  set ff='/home/mieg/AW/2025_Chris_Stacy/Chris_HUH7_SRA_RNA_5687_SraRunInfo.txt'
  set ff='/home/mieg/AW/Human_DATA/SEQC2/PacBio_A_excellenteQuality_SraRunInfo.txt'
  set date=2025_09_13
  set ff='/home/mieg/AW/Human_DATA/20250904_Roberto-Landsman_HepG2-CTCF_ChIP/SRA_2025_HEPG2_CTCF_ChIP_and_controls_Philadelphia_runinfo.txt'
  set date=2026_02_01
  set ff='/home/mieg/AW/Human_DATA/GIAB_eskSonHG002_RNA/20260202_HG002_SraRunInfo.txt'
endif
if ($species == rn) then
  set date=2016_11_04
  set ff='/home/mieg/ACEVIEWHELP/Rat_DATA/2016Nov05_6365_Rattus_norvegicus_Organism_NOT_biomol_dna_PropertiesSraRunInfo.txt'
  set date=2019_02_27
  set ff='/home/mieg/ACEVIEWHELP/Rat_DATA/20190227_Rattus_norvegicus_RNA_12094_SraRunInfo.txt'
  set date=2021_11_07
  set ff='/home/mieg/ACEVIEWHELP/Rat_DATA/20211107_Rat_miRNA_SRA_forLemingReview_SraRunInfo.txt'
  set date=2021_11_08
  set ff='/home/mieg/ACEVIEWHELP/Rat_DATA/20211107_3809miRNA_Rat_inSra_RunInfo.txt'
endif
if ($species == worm) then
  set date=2019_08_24
  set ff='/home/mieg/ACEVIEWHELP/Worm_DATA/20190828_SRA_RNA_runinfo.txt'
  set date=2024_03_21
  set ff='/home/mieg/ACEVIEWHELP/Worm_DATA/20230321_26255_RNAinSRA_Celegans_SraRunInfo.txt'
  set date=2024_10_01
endif
if ($species == mm) then
  set date=2016_11_04
  set ff='/home/mieg/ACEVIEWHELP/Mouse_DATA/2016Nov5_131768_mus_musculus_Organism_NOT_biomol_dnaPropertiesSraRunInfo.txt'
endif
if ($species == dog) then
  set date=2016_11_04
  set ff='/home/mieg/ACEVIEWHELP/Dog_DATA/Dog_3strandedRunsForTestTarget.txt'
  set date=2017_06_21
  set ff='/home/mieg/ACEVIEWHELP/Dog_DATA/Dog_RNA_729_June22_2017_SraRunInfo.txt'
endif
if ($species == Campylobacter || $species == jejuni) then
  set date=2017_03_07
  set ff='/home/mieg/ACEVIEWHELP/Viruses_microbes_DATA/Campylobacter_Skesa_OXA_SraRunInfo.txt'
  set ff='/home/mieg/ACEVIEWHELP/Viruses_microbes_DATA/CampylobacterProject_SRP063302_1952entries_SraRunInfo.txt'
  set ff='/home/mieg/ACEVIEWHELP/Viruses_microbes_DATA/Campylobacter_jejuni_genome_Truseq_SraRunInfo.txt'

  set ff='/home/mieg/ACEVIEWHELP/Viruses_microbes_DATA/Campylobacter_jejuni.Selected_Nextera_CDC_notNano.txt'
  set ff='/home/mieg/ACEVIEWHELP/Viruses_microbes_DATA/20170522_TruSeqNano_CDC.txt'
endif
if ($species == Lysobacter) then
  set date=2019_01_28
  set ff='/home/mieg/ACEVIEWHELP/Viruses_microbes_DATA/2019_Jan28Michael_Galperin_Lysobacter_enzymogenes_SraRunInfo.txt'
endif
if ($species == TD) then
  set date=2023_12_11
  set ff1='/home/mieg/AW/tardigrade_DATA/20231208_16_DNA_Ramazzottius_varieornatus_SraRunInfo.txt'
  set ff2='/home/mieg/AW/tardigrade_DATA/20231208_161_RNA_Ramazzottius_varieornatus_SraRunInfo.txt'
  cat $ff1 $ff2 > SRX_DB/TD.runInfo.txt
  set ff=`pwd`/SRX_DB/TD.runInfo.txt
endif

echo "species=$species date=$date"
set today=`echo $date | sed -e 's/_/-/g'`
set dd=DATA/$date
setenv ici `pwd` 
set ff1=$dd/runInfo.tsv


# Create a work directory
if (! -d DATA) mkdir DATA
if (! -d DATA/$date) mkdir DATA/$date
if (! -d $dd) mkdir $dd

# copy the RunInfo file in Unix format and transform to tab delimited
cat $ff | sed -e 's/\r//' | grep -v "Synthetic-Long-Read" > $dd/runInfo.tsv
set ff1=$dd/runInfo.tsv

echo "..... requested phases : $phaseSet"
foreach  phase ($phaseSet)

echo -n "=== phase $phase start "
date


if ($phase == SRR) goto phaseSRR
if ($phase == SRP) goto phaseSRP
if ($phase == PRJ) goto phasePRJ
if ($phase == GEO) goto phaseGEO
if ($phase == Sample) goto phaseSample
if ($phase == SRX) goto phaseSRX
if ($phase == Files) goto phaseFiles
if ($phase == Sublibs) goto phaseSublibs
if ($phase == Titles) goto phaseTitles
if ($phase == srr2run) goto phase_srr2run
if ($phase == srr2srx) goto phase_srr2srx
if ($phase == Papers) goto phasePapers

if ($phase == sraDownload) goto phase_sraDownload
if ($phase == sraDownloadTest) goto phase_sraDownload

usage:
echo "usage: SRX_import.tcsh SRR PRJ SRP GEO Sample SRX Files Papers  Sublibs Titles  srr2srx srr2run |  sraDownload sraDownloadTest"
goto phaseLoop

############
phaseSRR:

if (! -d SRX_DB) then
  ln -s ~/ace/waligner/metaData
  ln -s ~/ace/waligner/scripts
  ln -s ~/ace/bin.ICC_centos7 bin
  mkdir SRX_DB
  pushd  SRX_DB
    ln -s ../metaData/wspec.SRX wspec
    mkdir database
    echo y | ../bin/tacembly .
  popd
endif

############

# parse the SRR number, use it as intermediary Run name
cat $ff1 | gawk -F '\t' -f scripts/SRX_import.1.awk today=$today | grep -v Submission_dateOK > $dd/run_info.ace

bin/tacembly SRX_DB <<EOF
  read-models
  parse $dd/run_info.ace
  edit project $MAGIC
  save
  quit
EOF
goto phaseLoop

####### Parse the SRP objects i.e. title and abstracts for the runs
phaseSRP:
if (! -d $dd/SRP) mkdir $dd/SRP
bin/tacembly SRX_DB <<EOF
  query find srp srr && ! title && is srp*
  select -o $dd/srp.list select srp from srp in class srp where srp#srr and not srp#title
  quit
EOF

if (-e  $dd/SRP/_wget) \rm  $dd/SRP/_wget

ls  $dd/SRP | gawk '/html/{gsub(".html","",$1);print $1;}' >   $dd/srp.list2
cat $dd/srp.list $dd/srp.list2 $dd/srp.list2 | gawk '{n[$1]++;}END{for (k in n) if (n[k]==1)print k}' >  $dd/srp.list3

cat $dd/srp.list3 | gawk '{printf("wget -O %s/SRP/%s.html \"https://trace.ncbi.nlm.nih.gov/Traces/sra?acc=%s\"\n",dd,$1,$1);}' dd=$dd >  $dd/SRP/_wget

if (-e  $dd/SRP/_wget) then
  wc  $dd/SRP/_wget
  source  $dd/SRP/_wget
endif

if (-e  $dd/srp.ace) \rm  $dd/srp.ace
foreach srp (`cat $dd/srp.list`)
   if (! -e  $dd/SRP/$srp.html) continue 
   cat  $dd/SRP/$srp.html | gawk -f scripts/SRX_import.SRP.awk srp=$srp | sed -e 's/\\\"//g' -e s'/\\//g'  >> $dd/srp.ace
end

bin/tacembly SRX_DB <<EOF
  read-models
  parse $dd/srp.ace
  save
EOF
goto phaseLoop

####### Parse the PRJ objects i.e. title and abstracts for the runs
phasePRJ:
if (! -d $dd/PRJ) mkdir $dd/PRJ
bin/tacembly SRX_DB <<EOF
  query find srp srr && ! title && IS prj*
  select -o $dd/prj.list @
  quit
EOF

if (-e  $dd/PRJ/_wget) \rm  $dd/PRJ/_wget

ls  $dd/PRJ | gawk '/html/{gsub(".html","",$1);print $1;}' >   $dd/prj.list2
cat $dd/prj.list $dd/prj.list2 $dd/prj.list2 | gawk '{n[$1]++;}END{for (k in n) if (n[k]==1)print k}' >  $dd/prj.list3

cat $dd/prj.list3 | gawk '{printf("wget -O %s/PRJ/%s.html \"https://ncbi.nlm.nih.gov/bioproject/%s\"\n",dd,$1,$1);}' dd=$dd >  $dd/PRJ/_wget

if (-e  $dd/PRJ/_wget) then
  wc  $dd/PRJ/_wget
  source  $dd/PRJ/_wget
endif

if (-e  $dd/prj.ace) \rm  $dd/prj.ace
foreach prj (`cat $dd/prj.list`)
   if (! -e  $dd/PRJ/$prj.html) continue 
   cat  $dd/PRJ/$prj.html | gawk '{printf(" %s",$0);}' | gawk -f scripts/SRX_import.PRJ.awk prj=$prj | sed -e 's/\\\"//g' -e s'/\\//g'  >> $dd/prj.ace
end

bin/tacembly SRX_DB <<EOF
  read-models
  parse $dd/prj.ace
  save
EOF
goto phaseLoop

####### Parse the GEO to find the geo->contibutors and affiliations
phaseGEO:
if (! -d $dd/GEO) mkdir $dd/GEO

bin/tacembly SRX_DB <<EOF
  select -o $dd/geo.list g from g in ?geo 
EOF

if (-e $dd/GEO/_wget) \rm $dd/GEO/_wget
foreach geo (`cat $dd/geo.list`)
   if ( -e  $dd/GEO/$geo.html) continue 
   echo "wget -O "$dd"/GEO/"$geo".html  "'"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='$geo'"' | grep -v '(' | grep -v ')' >>  $dd/GEO/_wget
end

if (! -e $dd/GEO/_wget) goto phaseLoop

wc $dd/GEO/_wget
source  $dd/GEO/_wget

if (-e  $dd/geo.ace) \rm  $dd/geo.ace
foreach geo (`cat $dd/geo.list`)
   if (! -e  $dd/GEO/$geo.html) continue 
   cat  $dd/GEO/$geo.html | gawk -f scripts/SRX_import.3.awk geo=$geo  >> $dd/geo.ace
end

bin/tacembly SRX_DB <<EOF
  read-models
  parse $dd/geo.ace
  save 
EOF

goto phaseLoop

####### Parse the biosamples to get the attributes (organism, tissues etc)
phaseSample:
if (! -e $dd/BIOSAMPLE) mkdir $dd/BIOSAMPLE

bin/tacembly SRX_DB <<EOF
  query find biosample ! Title
  // query find project IS $MAGIC ; > srr ; >biosample
  select -o $dd/biosample.list @
EOF


ls  $dd/BIOSAMPLE | gawk '/html/{gsub(".html","",$1);print $1;}' >   $dd/biosample.list2
cat $dd/biosample.list $dd/biosample.list2 $dd/biosample.list2 | gawk '{n[$1]++;}END{for (k in n) if (n[k]==1)print k}' >  $dd/biosample.list3

cat $dd/biosample.list3 | gawk '{printf("wget -O %s/BIOSAMPLE/%s.html \"https://www.ncbi.nlm.nih.gov/biosample/?term=%s\"\n",dd,$1,$1);}' dd=$dd >  $dd/BIOSAMPLE/_wget

wc $dd/BIOSAMPLE/_wget
source  $dd/BIOSAMPLE/_wget

if (-e  $dd/biosample.ace) \rm  $dd/biosample.ace
foreach biosample (`cat $dd/biosample.list`)
   if (! -e  $dd/BIOSAMPLE/$biosample.html) continue 
   cat  $dd/BIOSAMPLE/$biosample.html | gawk -f scripts/SRX_import.4.awk biosample=$biosample  >> $dd/biosample.ace
end

bin/tacembly SRX_DB <<EOF
  read-models
  parse $dd/biosample.ace
  save 
EOF

goto phaseLoop

####### Parse the SRX get details on the preparation of the library
phaseSRX:
if (! -d $dd/SRX) mkdir $dd/SRX

bin/tacembly SRX_DB <<EOF
  query find project IS $MAGIC ; >SRR ; >SRX ; ! title
  select -o  $dd/srx.list @
EOF

ls  $dd/SRX | gawk '/html/{gsub(".html","",$1);print $1;}' >   $dd/srx.list2
cat $dd/srx.list $dd/srx.list2 $dd/srx.list2 | gawk '{n[$1]++;}END{for (k in n) if (n[k]==1)print k}' >  $dd/srx.list3

cat $dd/srx.list3 | gawk '{printf("wget -O %s/SRX/%s.html \"https://www.ncbi.nlm.nih.gov/sra/?term=%s\"\n",dd,$1,$1);}' dd=$dd >  $dd/SRX/_wget

wc  $dd/SRX/_wget
source  $dd/SRX/_wget

if (-e  $dd/srx.ace) \rm  $dd/srx.ace
foreach srx (`cat $dd/srx.list`)
   if (! -e  $dd/SRX/$srx.html) continue 
   cat  $dd/SRX/$srx.html | gawk -f scripts/SRX_import.5.awk srx=$srx  >> $dd/srx.ace
end

bin/tacembly SRX_DB <<EOF
  read-models
  parse $dd/srx.ace
  save 
  // query find SRX
  // list -a -f $dd/SRX.list
EOF

goto phaseLoop

######################
## add the file description

phaseFiles:
tbly SRX_DB <<EOF
  query find SRR !File && !Solid && !Sublibraries && Paired_end && project == $MAGIC
  select -o  $dd/r2f.p @
  query find SRR !File && !Solid && !Sublibraries && !Paired_end && project == $MAGIC
  select -o  $dd/r2f.s @
  query find SRR !File && Solid && !Sublibraries && Paired_end && project == $MAGIC
  select -o  $dd/r2f.p.s @
  query find SRR !File && Solid && !Sublibraries && !Paired_end && project == $MAGIC
  select -o  $dd/r2f.s.s @
EOF

cat $dd/r2f.p  | gawk  '{gsub(/\"/,"",$0);printf("SRR %s\n-D File\nFile fasta/1 %s/SRA/%s_1.fasta.gz\nFile fasta/2 %s/SRA/%s_2.fasta.gz\n\n",$1,dd,$1,dd,$1);}' dd=$dd > $dd/r2f.ace
cat $dd/r2f.s  | gawk  '{gsub(/\"/,"",$0);printf("SRR %s\n-D File\nFile fasta %s/SRA/%s.fasta.gz\n\n",$1,dd,$1,dd,$1);}' dd=$dd >> $dd/r2f.ace
cat $dd/r2f.s.s  | gawk  '{gsub(/\"/,"",$0);printf("SRR %s\n-D File\nFile csfasta %s/SRA/%s.fasta.gz\n\n",$1,dd,$1,dd,$1);}' dd=$dd >> $dd/r2f.ace
cat $dd/r2f.p.s  | gawk  '{gsub(/\"/,"",$0);printf("SRR %s\n-D File\nFile csfasta/1 %s/SRA/%s_1.fasta.gz\nFile csfasta/2 %s/SRA/%s_2.fasta.gz\n\n",$1,dd,$1,dd,$1);}' dd=$dd >> $dd/r2f.ace

echo "pparse  $dd/r2f.ace" | tbly SRX_DB -no_prompt

goto phaseLoop

#############
## add the super titles
phaseTitles:

echo "create the titles "
bin/sra_metadata -db SRX_DB -a -dbEdit -p $MAGIC 
bin/sra_metadata -db SRX_DB -s -dbEdit -p $MAGIC 

tbly SRX_DB <<EOF
  date
  query find project IS $MAGIC ; > srr
  bql -a -o r2b2t.txt select r,b,m,T from r in @, b in r->biosample, m in r->magic_sample2 where m,T in r->Manual_title
  date
EOF

echo "list the attributes"
tbly SRX_DB <<EOF
  date
  query find project IS $MAGIC ; > srr ; >biosample
  bql -a -o b2t.txt select b,t from b in @, t in b->title
  query find project IS $MAGIC ; > srr ; biosample
  bql -a -o rb2p.txt select r,b,t,t2 from r in @, b in r->biosample, t in b->biosample_attribute, t2 in t[1]
  quit
EOF

echo "establish the collection dates"
cat rb2p.txt | gawk -F '\t' '/collection date/{z=$4;gsub("collection date","",z) ;gsub("\"","",z);z=substr(z,1,10);i=index(z,"-");if (i>1)printf("SRR %s\nCollection_date %s\n\n",$1,z);}' > r2collection_date.ace
cat rb2p.txt | gawk -F '\t' '/ENA-LAST-UPDATE/{print;}/last update/{printf("%s\t%s\n",$1,$4);}' | sort -V | gawk -F '\t' '{z=$2;gsub("\"","",z);z=substr(z,1,10);i=index(z,"-");if (i>1)printf("SRR %s\nSubmission_date %s\n\n",$1,z);}' > r2submission_date.ace

cat rb2p.txt | gawk -F '\t' '/collection date/{print $4;}' | grep \- | grep -v "not collected" | sed -e 's/collection date//' | sort -u

cat r2b2t.txt | gawk -F '\t' '{T=$3;if($4 != "NULL")t=$4;printf("SRR %s\nTitle %s\n\n", $1,$3,$4);}' > r2b2t.ace
echo "parse file r2b2t.ace "
tags r2b2t.ace
tbly SRX_DB <<EOF
  pparse r2b2t.ace
  save
  quit
EOF

#############
## count the SRX spots cumulating the spots of the SRR sublibraries

echo "count the spots in sublibs "
if (! -e ZZZZZ) echo ZZZZZ > ZZZZZ

tace SRX_DB <<EOF
  query find SRR sublibraries
  select -a -o $dd/SRX_count_spots.txt srx,r,s1,s2,s3,s4,s5,s6,s7,s8,s9 from srx in @, r in srx->sublibraries where r,s in r#spots where s,s1 in s[1],s2 in s[2],s3 in s[3],s4 in s[4],s5 in s[5],s6 in s[6],s7 in s[7],s8 in s[8],s9 in s[9]
EOF
cat  $dd/SRX_count_spots.txt ZZZZZ | gawk -F '\t' '{gsub(/\"/,"",$0);}{if($1 != srx) {if(ns>0){printf("SRX %s\nSpots %d ",srx,ns); if(nb>0) { printf(" bases_in_SRA %d ", nb); if(nA >0){printf(" Average_length %d ",nA); if (nSize > 0){printf(" Insert_size %s ", nSize); if( nMates > 0) printf(" spots_with_mate %d ", nMates);}}} printf("\n\n");} ns = 0 ; nb = 0 ; nA = 0 ; nSize = 0 ; nMates = 0 ;} srx =$1;if ($3 > 0) { ns += $3;nb += $5; nMates += $11; nA = (nA * (ns - $3) + $7 * $3)/(ns); nSize = (nSize * (ns - $3) + $9 * $3)/(ns);}}' >  $dd/SRX_count_spots.ace

echo "parse the spots counts"
tbly SRX_DB <<EOF
  parse  $dd/SRX_count_spots.ace
  save
  quit
EOF

goto phaseLoop

tbly SRX_DB << EOF
    query find srr biosample AND project == $MAGIC
    bql -a -o $dd/title.txt  select srr,t from srr in @, b in srr->biosample , t in b->magic_sample2 
    save
    quit
EOF
cat  $dd/title.txt | gawk -F '\t' '/^#/{next;}{printf ("SRR %s\nTitle %s\n\n", $1,$2);}' > $dd/titles.ace
tbly SRX_DB << EOF
    pparse  $dd/titles.ace
    save
    quit
EOF

goto phaseLoop

#############
## get the papers
phasePapers:
setenv CVSROOT /home/mieg/VV/CODE/CVSROOT
if (! -d $dd/PAPERS) then
  mkdir $dd/PAPERS
  pushd  $dd/PAPERS
  cvs checkout IMPORT_DATA/BIBLIO/PmImport
  mv IMPORT_DATA/BIBLIO .
  popd
endif

if (-e $dd/PAPERS/a5.pmnocit.list) \rm $dd/PAPERS/a5.pmnocit.list
tbly SRX_DB <<EOF
  query find paper IS pm* AND NOT citation
  list -a -f $dd/PAPERS/a5.pmnocit.list
EOF

if (! -e $dd/PAPERS/a5.pmnocit.list) goto phaseLoop
set n=`wc -l  $dd/PAPERS/a5.pmnocit.list | gawk '{print 0+$1;}'`
if ($n == 0) goto phaseLoop

wc -l  $dd/PAPERS/a5.pmnocit.list 
pushd  $dd/PAPERS
  perl  BIBLIO/PmImport/medlineGet.pl < a5.pmnocit.list >! a5.pmnocit.gb
  perl  BIBLIO/PmImport/medline2ace.pl < a5.pmnocit.gb >! a5.pmnocit.preace
  cat a5.pmnocit.preace | grep -v EC_symbol > a5.pmnocit.ace
popd


tbly SRX_DB << EOF
    pparse   $dd/PAPERS/a5.pmnocit.ace
    save
    quit
EOF

goto phaseLoop

#######  create the sublibraries
phaseSublibs:

echo ' ' > $dd/multisrx.txt
tbly SRX_DB << EOF
    query find project IS $MAGIC ; >srr ; >srx ; COUNT { >SRR ; File}  > 1
    bql -a -o $dd/multisrx.txt select srx, srr from srx in @, srr in srx->srr where srr#File
    save
    quit
EOF

echo ' ' >  $dd/multisrx.ace
cat $dd/multisrx.txt | gawk -F '\t' '/^#/{next;} {printf("SRR %s\nSublibraries %s\n\n", $1,$2);}' > $dd/multisrx.ace
tbly SRX_DB << EOF
    pparse  $dd/multisrx.ace
    save
    quit
EOF
echo ' ' >  $dd/sublib.ace
tbly SRX_DB << EOF
    query  find project $MAGIC ; >srr ; sublibrary_of
    show -a -f $dd/sublib.ace
    save
    quit
EOF

cat $dd/multisrx.txt ZZZZZ $dd/sublib.ace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){srx[$2] = $1;next;}}/^$/{if(ok==1)print;ok=0;next;}/^SRR_download/{next;}/^SRR/{gsub(/\"/,"",$2);ok=0;z = srx[$2];if(length(z)<3)next;ok=1;printf("\nSRR %s\n",z);next;}/^File/{next;}/^Sublibrary_of/{next;}{if(ok==1) print}' >  $dd/sublib.ace2
tbly SRX_DB << EOF
    // pparse  $dd/sublib.ace2
    save
    quit
EOF

goto phaseLoop

#######  done

bly SRX_DB &

#############
## download the actual fasta files

phase_sraDownload:

if (! -d  $dd/SRA ) mkdir $dd/SRA
if (! -d  Fastc ) mkdir Fastc
# check which we already know   // key toto.a2.list

tbly SRX_DB <<EOF
  select  -o $dd/SRA/srrP.txt1  r from p in ?project where p == $MAGIC, r in p->SRR where r#file
  bql -o $dd/SRA/srrP.txt2 select r,p,d from r in @, p in r#paired_end, d in r->SRR_download
EOF
wc $dd/SRA/srrP.txt[12] 
 cat  $dd/SRA/srrP.txt2 | gawk -F '\t' '{srr=$1;pair=$2;file=$3;if( index(file,"https")==1)printf ("%s___%s___%s\n",srr,file,pair);}' > $dd/SRA/srr.todo
wc $dd/SRA/srr.todo

if (-e totosraNeeded) \rm totosraNeeded
touch totosraNeeded
foreach ss (`cat  $dd/SRA/srr.todo`)
  set srr=`echo $ss | gawk -F '___' '{print $1}'`
  if (-e Fastc/$srr/f2.2.fastc.gz) then
    if (-e $dd/SRA/$srr.fasta.gz || -e $dd/SRA/$srr'_1'.fasta.gz || -e $dd/SRA/$srr.csfasta.gz || -e $dd/SRA/$srr'_1'.csfasta.gz) then
      echo "removing $dd/SRA/$srr...fasta.gz"
      \rm $dd/SRA/$srr.*fasta.gz  $dd/SRA/$srr'_'*fasta.gz
    endif
  endif
  if (! -d Fastc/$srr &&  ! -e $dd/SRA/$srr && ! -e $dd/SRA/$srr.fasta.gz &&  ! -e $dd/SRA/$srr'_1'.fasta.gz && ! -e $dd/SRA/$srr.csfasta.gz &&  ! -e $dd/SRA/$srr'_1'.csfasta.gz) then
    echo $ss >> totosraNeeded
   if ($phase == sraDownload) then
     echo "scripts/submit $dd/SRA/$srr scripts/SRX_import_fasta.tcsh $dd/SRA $ss"
     scripts/submit $dd/SRA/$srr "scripts/SRX_import_fasta.tcsh $dd/SRA $ss"
     sleep 1
   endif
  endif
end
wc totosraNeeded

goto phaseLoop
ls DATA/2016_01_08/SRA/*fasta* | gawk '{split($1,aa,"/");split(aa[4],bb,"_");split(bb[1],cc,".");printf("Run %s\n", cc[1]);}' | sort -u > sraReady.list
ls DATA/2016_01_08/SRA/*fasta* | gawk '{split($1,aa,"/");split(aa[4],bb,"_");split(bb[1],cc,".");printf("SRR %s\n", cc[1]);}' | sort -u > sraReady.srr.list
######################
## clean up
######################
## Transfer to RunMasterDB

phase_srr2run:

tbly SRX_DB <<EOF
  query find project IS ?MAGIC ; > SRR ; > SRP ;  > geo ; reference && srp
  bql -a -o $dd/geo2srp2ref.txt  select geo,srp,ref from geo in @, srp in geo->srp, ref in geo->reference
  query find longtext
  show -a -f  $dd/longtext.ace
  quit
EOF
 cat $dd/geo2srp2ref.txt | gawk '/^"/{printf ("SRP %s\nReference \"%s\"\n\n",$2,$3);}' > $dd/geo2srp2ref.ace
 echo "pparse $dd/geo2srp2ref.ace " | tbly  SRX_DB -no_prompt

tbly SRX_DB <<EOF
  query find project $MAGIC ; >srr
  // key sraReady.srr.list
  spush
  follow sublibraries
  sor
  spop
  kstore kk
  show -a -f $dd/srr.dump.ace
  query follow srp
  show -a -f $dd/srp.preace
  kget kk
  query follow srx
  show -a -f $dd/srx.preace
  kget kk
  query Nucleic_Acid_Extraction
  show -a -f $dd/srr.nuc.ace Nucleic_Acid_Extraction
  kget kk
  query sraNucleic_Acid_Extraction
  show -a -f $dd/srr.sranuc.ace sraNucleic_Acid_Extraction
  find paper
  show -a -f $dd/srr.paper.preace
  follow abstract
  show -a -f $dd/srr.paper_abstract.ace
  query find project $MAGIC ; >srr
  follow biosample
  show -a -f $dd/srr.biosample.preace
  kget kk
  query  Sample_builder
  show -a -f $dd/srr.sample_builder.preace Sample_builder
  kget kk
  follow compare
   show -a -f $dd/compare.preace
  quit
EOF

echo ZZZZZ > ZZZZZ
# cat $dd/srr.dump.ace   ZZZZZ $dd/srr.dump.ace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){if($1 == "SRR"){ srr=$2;srr2run[srr];}if($1=="Run")srr2run[srr]=srr;next;}}/^Run/{next}/^SRR /{printf("Run %s\nSRR %s\n",srr2run[$2],$2);next;}{print}' > $dd/srr2run.preace
cat $dd/srr.dump.ace > $dd/srr2run.preace

cat <<EOF > $dd/srr2run.awk
/^Adult/{print;next;}
/^Age/{print;next;}
/^Annotation_problem/{print;next;}
/^Biosample/{next;}
/^Body_site/{next;}
/^Cap_CAGE_RACE/{print;next;}
/^Carcass/{print;next;}
/^Cell_line/{print;next;}
/^Center_name/{next;}
/^Connective/{print;next;}
/^Control/{print;next;}
/^Date_received/{print;next;}
/^Dauer/{print;next;}
/^Digestive/{print;next;}
/^Driver_gene/{print;next;}
/^EE/{print;next;}
/^ERROR/{print;next;}
/^Embryo/{print;next;}
/^Epithelium/{print;next;}
/^Endocrine_gland/{print;next;}
/^Exocrine_gland/{print;next;}
/^Endocrine/{print;next;}
/^Female/{print;next;}
/^Fetal/{print;next;}
/^File/{print;next;}
/^Forward/{next;}
/^Genotype/{print;next;}
/^Gene_selection/{print;next;}
/^Gene_ectopic/{print;next;}
/^Gene_down/{print;next;}
/^Gene_up/{print;next;}
/^Germline_and_development/{print;next;}
/^Group/{print ; next;}
/^Head/{print;next;}
/^Helicos/{print;next;}
/^Hemic_and_immune/{print;next;}
/^Hermaphrodite/{print;next;}
/^Illumina/{print;next;}
/^Ion_Torrent/{print;next;}
/^Juvenile/{print;next;}
/^L1/{print;next;}
/^L2/{print;next;}
/^L3/{print;next;}
/^L4/{print;next;}
/^Any_larva/{next;}
/^Larva/{print;next;}
/^Library_name/{next;}
/^Lymphoid/{print;next;}
/^Magic_author2/{gsub("Magic_author2" ,"Author", \$0);print; next;}
/^Male/{print;next;}
/^Males_and_females/{print;next;}
/^Marker_gene/{print;next;}
/^Membrane/{print;next;}
/^Microbiome/{print;next;}
/^Mixed_embryo/{print;next;}
/^Mixed_larvae/{print;next;}
/^Mixed_sex/{next;}
/^Mixed_stages/{print;next;}
/^Molt/{print;next;}
/^Muscle/{print;next;}
/^Musculoskeletal_and_cardiovascular/{print;next;}
/^Nanopore/{print;next;}
/^Nascent_RNA/{print;next;}
/^Nerve/{print;next;}
/^Nervous_and_sense_organ/{print;next;}
/^No_tissue/{print;next;}
/^OA/{print;next;}
/^Other/{print;next;}
/^Oxford_nanopore/{print "Nanopore";next;}
/^Paired_end/{print;next;}
/^PacBio/{print;next;}
/^Phenotype/{print;next;}
/^Project/{print;next;}
/^Pupa/{print;next;}
/^RIP_CLIP/{print;next;}
/^RNA/{print;next;}
/^Reference/{print;next;}
/^Release_date/{print;next;}
/^Replicate/{print;next;}
/^Respiratory/{print;next;}
/^Ribosome_profiling/{print;next;}
/^Roche_454/{print;next;}
/^Run/{next;printf("\n");print;next;}
/^SNP/{print;next;}
/^SOLiD/{print;next;}
/^SRP/{next;}
/^SRR_download/{print; next;}
/^SRR/{printf ("\n\nRun %s\n-D Nucleic_Acid_Extraction\n-D Origin\n-D group\n-D Title\n-D SRR\n-D Sample\n-D SRP\n",\$2); x=substr(\$2,2,3); if(x=="DRR" || x=="SRR" || x=="ERR")printf("SRR %s\n",\$2); next;}
/^SRX/{next;}
/^Sample_name/{next;}
/^Small_RNA/{print;next;}
/^SOLiD/{print;next;}
/^Sorting_title/{print;next;}
/^Species/{print;next;}
/^Spots/{print;next;}
/^Strain/{print;next;}
/^Stranded/{next;print "Strand";next;}
/^Sublibraries/{print;next;}
/^Sublibrary_of/{print;next;}
/^Submission_date/{print;next;}
/^Title/{print;next;}
/^Treatment/{print;next;}
/^Total_RNA/{print;next;}
/^Total/{print "Total_RNA" ; next;}
/^Union_of/{print;next;}
/^Unspecified_RNA/{print;next;}
/^Urogenital/{print;next;}
/^Warning/{next;}
/^Whole_genome/{print;next;}
/^Whole_organism/{print;next;}
/^YA/{print;next;}
/^nanopore/{print "Nanopore";next;}
/^nonStranded/{print;next;}
/^polyA/{print;next;}
/^sraUnspecified_RNA/{next;}
/^sraGenesraGene_selection/{next;}
/^sraGenesraGenesraPolyA/{next;}
/^sraPolyA/{next;}
/^sraUnspecified_RNA/{next;}
/^sraRIP_CLIP/{next;}
/^sraSmall_RNA/{next;}
/^sraUnspecified_RNA/{next;}

END {printf("\n");}
EOF

cat  $dd/srr2run.preace | gawk -f $dd/srr2run.awk >  $dd/srr2run.ace

cat <<EOF > $dd/srr2srr.awk
/^SRR_download/{next;}
/^SRR/{printf("\n");print;next;}
/^Species/{print;next;}
/^Stranded/{print;next;}
END {printf("\n");}
EOF
cat  $dd/srr2run.preace | gawk -f $dd/srr2srr.awk >  $dd/srr2srr.ace

cat <<EOF > $dd/srp2run.awk
/^Abstract/{print;next;}
/^Description/{print;next;}
/^Comment/{print;next;}
/^GEO/{next;}
/^Identifier/{print;next;}
/^Reference/{print;next;}
/^SRP/{printf("\n");print;print "-D T" ; next;}
/^SRR_download/{next;}
/^SRR/{print;printf ("Run %s\n",\$2);next;}
/^Species/{print;next;}
/^Title/{print;next;}
END {printf("\n");}
EOF

cat  $dd/srp.preace | gawk -f $dd/srp2run.awk >  $dd/srp.ace

cat <<EOF > $dd/srx2run.awk
/^SRX/{printf("\n");print; print "-D T" ; next;}
/^Construction_protocol/{print;next;}
/^Design/{print;next;}
/^SRR/{print;next;}
/^Submitted_by/{print;next;}
/^Spots/{print;next;}
/^Title/{print;next;}
END {printf("\n");}
EOF

cat  $dd/srx.preace | gawk -f $dd/srx2run.awk >  $dd/srx.ace
cat  $dd/compare.preace |  sed -e 's/^SRR[ \t]/Runs /' > $dd/compare.ace
tbly MetaDB -no_prompt <<EOF
  read-models
  query find project IS $MAGIC ; > run
  edit -D group
  query find project IS $MAGIC 
  kill
  pparse  $dd/srr2run.ace
  pparse  $dd/srr2srr.ace
  pparse  $dd/srp.ace
  pparse  $dd/srx.ace
  pparse  $dd/srr.nuc.ace
  pparse  $dd/srr.sranuc.ace
  pparse  $dd/compare.ace
  kill
  pparse  $dd/compare.ace
  save
  pparse  $dd/longtext.ace
  save
  quit
EOF


cat <<EOF > $dd/srx2paper.awk
/^GEO/{next;}
{print}
END {printf("\n");}
EOF

cat  $dd/srr.paper.preace | gawk -f $dd/srx2paper.awk >  $dd/srr.paper.ace
echo "pparse  $dd/srr.paper.ace" | tbly MetaDB -no_prompt
echo "pparse  $dd/srr.paper_abstract.ace" | tbly MetaDB -no_prompt

cat <<EOF > $dd/srx2paper.awk
/^GEO/{next;}
{print}
END {printf("\n");}
EOF

cat $dd/srr.biosample.preace ZZZZZ  $dd/srr.sample_builder.preace | gawk '{line++;}/^ZZZZZ/{zz++;inside=0;next;}/^$/{inside=0;print;next;}/^Biosample_/{next;}/^Sample_name/{next;}/^Biosample/{bio=$2;if(length(bio)>0){printf("Sample %s\n-D T \n",$2);inside=1;}next;}/^Magic_sample2/{gsub("Magic_sample2","Title", $0);if(inside==1)print;next;}/^SRR/{if(zz<1){print; srr2bio[$2]=bio;}next;}{if(zz==1 && inside==1)print;next;}' > $dd/srr.sample.ace

echo "pparse  $dd/srr.sample.ace" | tbly MetaDB -no_prompt

tbly MetaDB <<EOF
  query find project IS $MAGIC ; >run ; sublibrary_of && Group
  bql -a -o rSubGr.txt  select r,sub,g from r in @ , sub in r->sublibrary_of, g in r->group
  undo
  bql -a -o rSubGr.txt1  select r,sub,g from r in @ , sub in r->sublibrary_of, g in r->group
EOF
cat  rSubGr.txt | gawk -F '\t' '{printf ( "Run %s\n-D Union_of %s\nUnion_of %s\n\n",$3,$1,$2 ) ; }' >  rSubGr.fix.ace

echo "pparse  rSubGr.fix.ace" | tbly MetaDB -no_prompt
\rm  rSubGr.fix.ace  rSubGr.txt

tbly SRX_DB  <<EOF
  query find project IS $MAGIC ; > srr ; > biosample
  show -a -f $MAGIC.biosample.preace T
  quit
EOF

cat  $MAGIC.biosample.preace | sed -e 's/^Biosample/Sample/' -e 's/^Sample_attribute/Biosample_attribute/' -e 's/^Submission/\!Submission/' >  $MAGIC.biosample.ace

tbly MetaDB <<EOF
  select -o srx2srr.txt srx,srr from srx in ?srx, srr in srx->srr
  quit
EOF
cat  srx2srr.txt | gawk -F '\t' '{printf ("Run %s\nSRX %s\n\n", $2,$1);}' >  srx2run.ace
tbly MetaDB <<EOF
  pparse srx2run.ace
  pparse $MAGIC.biosample.ace 
  save
  quit
EOF

tace SRX_DB <<EOF
  s -o subIsDate.txt a,b from a in ?author where a ~ "????-??-??", b in a->biosample
  quit
EOF

cat subIsDate.txt | gawk -F '\t' '{printf("Biosample %s\n-D Submission %s\nSubmission_date %s\n\n", $2,$1,$1);}' > subIsDate.ace1
cat subIsDate.txt | gawk -F '\t' '{printf("Sample %s\n-D Submission %s\nSubmission_date %s\n\n", $2,$1,$1);}' > subIsDate.ace

tbly MetaDB <<EOF
  select -o r2srr2s.txt run,srr,sample from run in ?run, srr in run->srr, sample in srr->sample where sample
  select -o r2srr2s2.txt run,srr,sample from run in ?run, run2 in run->sublibraries,  srr in run2->srr, sample in srr->sample where sample
EOF

cat r2srr2s.txt r2srr2s2.txt | gawk -F '\t' '{printf ("Run %s\nSample %s\n\n", $1, $3) ;} ' > run2sample.ace
tace SRX_DB <<EOF
  find keyset
  show -a -f ks.ace
EOF
tbly MetaDB <<EOF
  parse run2sample.ace
  parse ks.ace
  save
  quit
EOF
goto phaseLoop

#######################################################
## Inside SRX_db, transfer the manual annotations from SRR to SRX
phase_srr2srx:
tbly SRX_DB <<EOF
  query find project IS $MAGIC ; >srr ;
  query file && sublibrary_of && ! union_of
  select -a -o srr2srx.txt  srr,srx from srr in @, srx in srr->sublibrary_of
  show -a -f $MAGIC.srr.ace
EOF

echo ZZZZZ > ZZZZZ
cat  srr2srx.txt  ZZZZZ  $MAGIC.srr.ace | gawk '/^ZZZZZ/{zz++;next;}{if(zz<1){r2x[$1]=$2;next;}}/^File/{next;}/^SRR /{printf("SRR %s\n",r2x[$2]);next;}/^Sublib/{next;}/^Spots/{next;}{print}' >  $MAGIC.srr2srx.ace

tbly SRX_DB <<EOF
  pparse $MAGIC.srr2srx.ace
  save
  quit
EOF

goto phaseLoop

#######################################################


##########
# danielle score for droso runs
cat RESULTS/scoring.txt | sed -e 's/\r//' | gawk -F '\t' '{run=$2;k=0;x=$4;if(x>10000)k+=5;else if (x>8000)k+4;else if (x>5000)k+=2;else k++;x=$5;if(x<5)k+=3;else if (x<10)k+=2;else if (x<20)k+=.5;else if (x>50)k-=1;x=$6;if(x>4500)k+=4;else if (x>3000)k+=4;else if (x >2000)k+=1.5;else if (x>1000)k+=.5; x=$7;z="";if(x>99)z="++";else if (x >94)z="+";else if (x<1)z="--";else if (x<6)z="-";printf("Run %s\t%f\t%s\n",run,k,z);}' >  RESULTS/scoring.2.txt

cat RESULTS/scoring.2.txt | gawk '/^Rdm/{if($3=="++" && $2 >=11)print "Run " $1;}' > toto.Aplus.list
cat RESULTS/scoring.2.txt | gawk '/^Rdm/{if($3=="--" && $2 >=11)print "Run " $1;}' > toto.Aminus.list
cat RESULTS/scoring.2.txt | gawk -F '\t' '/^#/{next;}{printf("Run %s\nScore %s \"%s\"\n\n",$1,$2,$3);}' >  RESULTS/scoring.2.ace

#############


phaseLoop:
 date
end
