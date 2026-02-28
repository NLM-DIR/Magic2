#!bin/tcsh

set chrom=$1
set pass=$2

# export all introns support and type
  bin/tacembly  tmp/X.$MAGIC/$chrom << EOF
    s -o tmp/X.$MAGIC/$chrom/XI_ct_ac.txt t,s,c from s in ?sequence where s like "xi_*", c in s->composite, t in s->Intron_boundaries
    quit
EOF

# locate  the ct_ac/Other  XI_ echo introns, absorb them in gt_ag
  cat tmp/X.$MAGIC/$chrom/XI_ct_ac.txt  | sed -e 's/^Other/ct_ac/' | sort | gawk -F '\t' '{s=$2;split(s,aa,"__");split(aa[2],bb,".");split(bb[1],cc,"_");x1=cc[1];x2=cc[2];}/^ct_ac/{k++;ss[k]=s;z[k]=$3;for(i=-8;i<=8;i++){u=x2+i "_" x1+i;uu[u]=k;}next;}/gt_ag/{k=uu[bb[1]];if (k<1)next;if(z[k]*10 < $3)printf("Sequence %s\n",ss[k]);}' > tmp/X.$MAGIC/$chrom/XI_ct_ac.$pass.list

# purge the ct_ac XI_ echo introns
  bin/tacembly  tmp/X.$MAGIC/$chrom << EOF > tmp/X.$MAGIC/$chrom/f3.kill_ct_ac.pass$pass.out
    key tmp/X.$MAGIC/$chrom/XI_ct_ac.$pass.list
    spush
    find Est X*_-*  // parasites
    sor
    query find est XY_* ct_ac > 1
    sor
    spop
    spush
    follow from_gene
    kstore tg
    spop
    spush
    follow dna
    sor
    undo
    follow cdna_clone
    sor
    spop
    kill
    kget tg
    acem
      cdna_73
      quit
    save
    quit
EOF
