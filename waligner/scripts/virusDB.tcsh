#!/bin/tcsh -f

if (! -d VirusDB) then
  mkdir VirusDB
  pushd VirusDB
  ln -s ../metaData/wspec
  whoami >> wspec/passwd.wrm
  mkdir database
  ../bin/tace . << EOF
y
  quit
EOF
  popd
endif


pushd VirusDB
set ff=/home/mieg/AW/Viruses_microbes_DATA/ChrisMason_MissingVirusesApr2020/New_Virus_list_556_virus_HandEditedDanielle_April2020.txt
cat $ff | gawk -F '\t' '{for(i=1;i<=NF;i++)print $i;}' | gawk -F ',' '{for(i=1;i<=NF;i++)print $i;}' | gawk  '/^N[CR]/{print}/^KP/{print}' | gawk '{split($1,aa,".");print aa[1];}' | sort -u > ~/AW/Danou.list
# marche pas: on obtient des melanges
#  gawk '{printf("seqfetch -t 5 -dn -s %s > %s.fasta1\n", $1,$1);}' > _m ; source _m > titi ; cat titi | dna2dna -I fasta -o titi2 >$ /dev/null
# web NCBI all-ressource; cherchez ncbi/sites/batchentrez (choose file Danou.list, db=nucleotide; puis retrieve => download fasta at full-GB)
 cat /home/mieg/AW/Viruses_microbes_DATA/ChrisMason_MissingVirusesApr2020/Danou_list.fasta | gawk '/^>/{s=substr($1,2);t=substr($0,length($1)+1); printf("Sequence %s\nTitle \"%s\"\n\n",s, t) ;}' > titi.ace
cat /home/mieg/AW/Viruses_microbes_DATA/ChrisMason_MissingVirusesApr2020/Danou_list.fasta  >> titi.ace

cat $ff | gawk -F '\t' '{n = split($1,aa,",");if(n==1)printf("Sequence %s\nAccession %s\nVirus\n\n",$1,$1);}' > v2acc.1.ace
cat $ff | gawk -F '\t' '{n = split($1,aa,",");if(n>1){printf("Sequence %s__\nVirus\n",aa[1]);for(i=1;i<=n;i++)printf("Subsequence %s\nAccession %s\n",aa[i],aa[i]);printf("\n");}}' > v2acc.n1.ace
cat $ff | gawk -F '\t' '{if ($3) {printf("Sequence %s__\nVirus\n",$3);for(i=3;i<=NF;i++){k=split($i,aa,".");if(k==2)printf("Subsequence %s\nAccession %s\n",aa[1],aa[1]);}printf("\n");}}' > v2acc.n2.ace

cat $ff | gawk -F '\t' '{n=split($1,aa,",");if(n>1)printf("Sequence %s__\nTitle \"%s\"\n\n",aa[1],$2);}' > title.ace
cat $ff | gawk -F '\t' '{if(length($3)>1)printf("Sequence %s__\nTitle \"%s\"\n\n",$3,$2);}' >> title.ace

# eliminate the .version number from all sequennces
tace . << EOF
  s -o toto s from s in ?sequence where s ~ "*.*"
EOF
cat toto | gawk -F '.' '{if (NF == 2) {printf("-R Sequence %s %s\n\n",$0,$1);}}' | grep -v SRR > v.rename.ace
cat toto | gawk -F '.' '{if (NF == 2) {printf("-R DNA %s %s\n\n",$0,$1);}}' | grep -v SRR >> v.rename.ace
tace . << EOF
  parse title.ace
  pparse v.rename.ace
  save
EOF
tace . << EOF
  find dna
  s -o toto @
EOF
cat toto | gawk -F '.' '{if (NF==2) {printf("-R DNA %s %s\n\n",$0,$1);}}' > dna.rename.ace
tace . << EOF
  pparse v.rename.ace
  save
EOF
tace . << EOF
  s ?sequence subsequence ;
  s -o toto s,s2,n from s in @, s2 in s->subsequence, n in s2->dna[2]
EOF
cat toto | gawk '{if($1 != old) x=1 ; old=$1; printf("Sequence %s\nSubsequence %s %d %d\n\n",$1,$2,x,x+$3-1);x=x+$3+100;}' > s2subseq.ace


# some sequencs were missing : find sequence virus ! subsequence && ! dna ; s -o danou2.list @ 
# then export from web

cat /home/mieg/AW/Viruses_microbes_DATA/ChrisMason_MissingVirusesApr2020/Danou2_78_Genbankfull_sequence.fasta | gawk '/^>/{s=substr($1,2);t=substr($0,length(s)+3);split(s,aa,".");a=aa[1];printf(">%s\n",a);next;}{print}' > toto.fasta
cat /home/mieg/AW/Viruses_microbes_DATA/ChrisMason_MissingVirusesApr2020/Danou2_78_Genbankfull_sequence.fasta | gawk '/^>/{s=substr($1,2);t=substr($0,length(s)+3);split(s,aa,".");a=aa[1];printf("Sequence %s\nAccession %s\nTitle \"%s\"\n\n",a,s,t);}' > toto.ace


tace . << EOF
  query find sequence virus
  dna -f danou.final.virus.fasta
  select -o toto.ln  s,ss,a,b from s in ?sequence, ss in s->subsequence, a in ss[1], b in ss[2] where b > 0
EOF
cat toto.ln | gawk '{s=$1;n[s] += $4-$3+1;}END{for (s in n) printf("Sequence %s\nLength %d\n\n", s,n[s]);}' > toto.ln.ace
tace . << EOF
  parse toto.ln.ace
  save
  quit
EOF
\rm toto.ln toto.ln.ace

mv danou.final.virus.fasta virus_manual_selection.2020_04_14.fasta
gzip virus_manual_selection.2020_04_14.fasta
\cp virus_manual_selection.2020_04_14.fasta.gz ~/CRN/TARGET/Targets


###############################################################################
## Bacteria

cd ~/CRN_Mason/VirusDB
pushd  /home/mieg/AW/Viruses_microbes_DATA/ChrisMason_missingImportantBacteria/
  foreach ff ( `ls *.fna.gz` )
    set s=`echo $ff | gawk -F _GCF '{print $1}'`
    zcat $ff | gawk '/^>/{a=substr($1,2);t=substr($0,length(s)+2);printf("Sequence %s\nTitle %s\nSubsequence %s\nBacteria\n\nSequence %s\nTitle \"%s\"\n\n",s,s,a,a,t);}' s=$s >> /tmp/tutu
  end
  zcat *.gz > /tmp/tutu1.fasta
popd


mv  /tmp/tutu1.fasta missingImportantBacteria.fasta
gzip missingImportantBacteria.fasta

dna2dna -i missingImportantBacteria.fasta,gz -getTM > missingImportantBacteria.TM
cat missingImportantBacteria.TM | gawk -F '\t' '/^#/{next;}{split($1,aa," ");s=aa[1];ln=$2;printf("Sequence %s\nLength %d\n\n",s,ln);}' > tutu1.ace

tace . << EOF
  select s from s in ?sequence where s#bacteria && s#subsequence ;
  s -o toto s,s2,n from s in @, s2 in s->subsequence, n in s2->length
EOF
cat toto | gawk '{if($1 != old) x=1 ; old=$1; printf("Sequence %s\nSubsequence %s %d %d\n\n",$1,$2,x,x+$3-1);x=x+$3+100;}' > s2subseq.ace
tace . << EOF
  parse s2subseq.ace
  select s from s in ?sequence where s#bacteria and s#subsequence ;

  s -o toto4 s,s2,x1,x2 from s in @, s2 in s->subsequence, x1 in s2[1], x2 in s2[2]
EOF
# we have a problem of ordering, which is fixed by exporting in bql and reimporting
cat toto4  | gawk '{if($1 != old) x=1 ; old=$1; printf("Sequence %s\nSubsequence %s %d %d\n\n",$1,$2,x,x+$3-1);x=x+$3+100;}' > s2subseq2.ace
tace . << EOF
  parse s2subseq.4.ace
  edit -D Subsequence
  parse s2subseq.4.ace
  save
  quit
EOF
# we can now reexport as single concatenated fasta files and edit the lengths
zcat missingImportantBacteria.fasta.gz | gawk '/^>/{s=substr($1,2);t=substr($0,length(s)+3);printf("Sequence %s\nTitle \"%s\"\n\n",s,t);}' > toto7.ace
tace . << EOF
  parse missingImportantBacteria.fasta.gz 
  parse toto7.ace
  query find sequence bacteria
  dna -f newbacteria.fasta
  query find sequence bacteria
  select -o toto5 s,ln,t,s2,ln2,t2 from s in @,ln in s->length , t in s->title, s2 in s->subsequence, ln2 in s2->length, t2 in s2->title
  select -o toto.ln  s,ss,a,b from s in ?sequence, ss in s->subsequence, a in ss[1], b in ss[2] where b > 0
  save
EOF
cat toto.ln | gawk '{s=$1;n[s] += $4-$3+1;}END{for (s in n) printf("Sequence %s\nLength %d\n\n", s,n[s]);}' > toto.ln.ace
tace . << EOF
  parse toto.ln.ace
  save
  quit
EOF
\rm toto.ln toto.ln.ace

###############################################################################
## Transposon

cd ~/CRN_Mason/VirusDB
pushd  /home/mieg/AW/Viruses_microbes_DATA/ChrisMason_missingImportantBacteria/
# there is a problem with the names in the fasta file, they come up as 3 columns in TM, iedit by hand
zcat TARGET/Targets/hs.transposon.TM.txt.gz | gawk -F '\t' '/^#/{next;}{split($1,aa," ");printf("Sequence %s\nTransposon\nLength %d\n\n",aa[1],$2);}' > toto.ln.ace
tace VirusDB << EOF
  parse toto.ln.ace
  save
  quit
EOF
\rm toto.ln.ace





 



