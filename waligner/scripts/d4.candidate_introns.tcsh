#!bin/tcsh -f

set group=$1

####### 2019_04 Start from d1.de_uno, reformat to match the known introns, merge the 2 files using the max number, report up to chosen min support
# combine the de_uno
# zcat tmp/OR/WG2_AGLR1_A_100_1_HS25_BHVYWVBCX2_L02_AACGCTTA/d1.WG2_AGLR1_A_100_1_HS25_BHVYWVBCX2_L02_AACGCTTA.de_uno.txt.gz | head
##Type	Strand	From target	a1	a2	To target	b1	b2	Boundary	Length	Number of supporting reads	Number of quasi perfect reads
#INTRON	Forward	16	106553	106613	16	107086	107124	gt_ag	472	1	1
#INTRON	Forward	Un|NT_167214.1	108944	109018	Un|NT_167214.1	109060	109084	gc_tc	41	1	1
#INTRON	Forward	Un|NT_167214.1	109360	109403	Un|NT_167214.1	109460	109515	gt_cc	56	1	1

# with the known introns from the intron analysis
# zcat tmp/INTRONRUNS/WG2_AGLR1_A_100_1_HS25_BHVYWVBCX2_L02_AACGCTTA/WG2_AGLR1_A_100_1_HS25_BHVYWVBCX2_L02_AACGCTTA.u.intronSupport.counts.gz | head
# G	000000000	000000000	5845727
# 1	000017232	000017056	1
# 1	000017605	000017369	1


###### Reformat the d1.de_uno in 4 columns -> d4.de_uno
## group case
if(! -e  tmp/OR/$group/d4.de_uno.txt.gz) then
  echo ' ' > tmp/OR/$group/d4.de_uno.txt.1
  set ok=0
  foreach run (`cat MetaDB/$MAGIC/g2r MetaDB/$MAGIC/r2sublib |  gawk -F '\t' '{if($1==g)print $2;}END{print g}' g=$group  | sort -u`) 
    if (-e  tmp/OR/$run/d1.$run.de_uno.txt.gz) then
      gunzip -c tmp/OR/$run/d1.$run.de_uno.txt.gz | gawk -F '\t' '/^INTRON/{a1=$5;a2=$7;if(a1<a2){a1++;a2--;}else{a1--;a2++;}z=$1 "\t" $2 "\t" $3 "\t0\t" a1 "\t" $6 "\t" a2 "\t0\t" $9"\t" $10 ; n11[z] += 0+$11 ; n12[z] += 0+$12;}END{for (z in n11)printf ("%s\t%d\t%d\n",z,n11[z],n12[z]) ;}' | gawk -F '\t' '{printf("%s\t%09d\t%09d\t%d\t%s\t%d\n",$3,$5,$7,$11,$9, $10)}' >>  tmp/OR/$group/d4.de_uno.txt.1
      set ok=1
    else if (-e  tmp/OR/$run/d4.de_uno.txt.gz) then
      zcat tmp/OR/$run/d4.de_uno.txt.gz >>  tmp/OR/$group/d4.de_uno.txt.1
      set ok=1
    endif
  end
  if ($ok == 1) then
    cat  tmp/OR/$group/d4.de_uno.txt.1  | gawk -F '\t' '{z = $1 "\t" $2 "\t" $3 ; n[z] += $4 ; t[z]= $5 "\t" $6 ; } END {for (z in n)printf ("%s\t%d\t%s\n", z, n[z], t[z]) ;}' | sort | gzip > tmp/OR/$group/d4.de_uno.txt.gz
  endif
  \rm tmp/OR/$group/d4.de_uno.txt.1
endif
  
## run case
if (0 && -e  tmp/OR/$group/d1.$group.de_uno.txt.gz && ! -e  tmp/OR/$group/d4.de_uno.txt.gz) then
    gunzip -c tmp/OR/$group/d1.$group.de_uno.txt.gz | gawk -F '\t' '/^INTRON/{a1=$5;a2=$7;if(a1<a2){a1++;a2--;} else{a1--;a2++;} z=$1 "\t" $2 "\t" $3 "\t0\t" a1 "\t" $6 "\t" a2 "\t0\t" $9"\t" $10 ; n11[z] += 0+$11 ; n12[z] += 0+$12;}END{for (z in n11)printf ("%s\t%d\t%d\n",z,n11[z],n12[z]) ;}' | gawk -F '\t' '{printf("%s\t%09d\t%09d\t%d\t%s\t%d\n",$3,$5,$7,$11,$9, $10)}' | gzip >  tmp/OR/$group/d4.de_uno.txt.gz
endif

##### Throw in the known introns

if(! -e  tmp/OR/$group/d4.known.txt.gz) then
  echo ' ' > tmp/OR/$group/d4.known.txt.1
  set ok=0
  if ($USEMAGICBLAST == 1) set ok=1  
  foreach run (`cat MetaDB/$MAGIC/g2r MetaDB/$MAGIC/r2sublib |  gawk -F '\t' '{if($1==g)print $2;}END{print g}' g=$group  | sort -u`) 
    if (-e  tmp/INTRONRUNS/$run/$run.u.intronSupport.counts.gz) then
      gunzip -c tmp/INTRONRUNS/$run/$run.u.intronSupport.counts.gz >>  tmp/OR/$group/d4.known.txt.1
      set ok=1
    endif
  end
  if ($ok == 1) then
    cat  tmp/OR/$group/d4.known.txt.1  | gawk -F '\t' '/G_Any_Intron/{next}{if($4+0<1)next;}{z = $1 "\t" $2 "\t" $3 ; n[z] += $4 ; } END {for (z in n)printf ("%s\t%d\n", z, n[z]) ;}' | sort | gzip > tmp/OR/$group/d4.known.txt.gz
  endif
  \rm tmp/OR/$group/d4.known.txt.1
endif
  
if (-e  tmp/INTRONRUNS/$group/$group.u.intronSupport.counts.gz && ! -e tmp/OR/$group/d4.known.txt.gz) then
  cp   tmp/INTRONRUNS/$group/$group.u.intronSupport.counts.gz tmp/OR/$group/d4.known.txt.gz
endif

###### merge de_uno and known
if (-e  tmp/OR/$group/d4.de_uno.txt.gz && -e tmp/OR/$group/d4.known.txt.gz) then
   gunzip -c tmp/OR/$group/d4.de_uno.txt.gz ZZZZZ.gz tmp/OR/$group/d4.known.txt.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}/G_Any_Intron/{next}{if($4+0<1)next;}{z = $1 "\t" $2 "\t" $3 ; nn[z] = 1 ; n[z,zz+1] += $4 ; if(NF > 5)  t[z]= $5 "\t" $6 ; if(zz == 1) known[z]="known" ; else known[z] = "new" ; } END { for (z in nn) {k = 0 + n[z,1]; if (0 + n[z,2] > k) k = 0 + n[z,2]; if (! t[z]) t[z]= "-\t-" ; printf ("%s\t%d\t%s\t%s\n", z, k, t[z],known[z]) ;}}' | sort > tmp/OR/$group/d4.de_uno.txt
  \rm tmp/OR/$group/d4.de_uno.txt.gz
endif
if (-e  tmp/OR/$group/d4.de_uno.txt.gz) then
  gunzip  tmp/OR/$group/d4.de_uno.txt.gz
endif

# select a min support so that at least 95% of the supporting reads are counted
if (! -e tmp/OR/SMAGIC.any.introns) then
  foreach target ($Etargets)
    cat  tmp/METADATA/$target.[fr].introns >> tmp/OR/$MAGIC.any.introns
  end
  cat tmp/OR/$MAGIC.any.introns | sort -u | sort -V > tmp/OR/$MAGIC.any.introns.1
  mv tmp/OR/$MAGIC.any.introns.1  tmp/OR/$MAGIC.any.introns
endif

# flag the de.uno introns donor acceptor as doubtfull
if (-e  tmp/OR/$group/d4.de_uno.txt) then
  set toto=tmp/OR/$group/d4.de_uno.txt
  # count the donors flag non gt_ag non ct_ac
  cat $toto ZZZZZ $toto | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{c=$1;a1=$2+0;a2=$3+0;n=$4+0;f=$5;if (zz<1){if (f=="ct_ac" || f=="gt_ag"){nn[c,a1]+=n;nn[c,a2]+=n;}next;}if (f=="ct_ac" || f=="gt_ag"){print ; next;}for (i=-2;i<=2;i++){n1=nn[c,a1+i];n2=nn[c,a2+i];if (0) print i,a1,n,n1,n2;if(n1>k*n || n2 > k*n){printf("%s\t*\n",$0);next;}}print;}' k=5 > $toto.1
  \mv $toto.1 $toto
endif

if (-e  tmp/OR/$group/d4.de_uno.txt) then
  set minS=`cat tmp/OR/$group/d4.de_uno.txt | gawk -F '\t' '{nn+=$4; n[$4]+=$4;}END{for(k in n)printf("%d\t%d\t%s\n",k,n[k],nn);}' | sort -k 1,1nr | gawk -F '\t' '{k=$1;n+=$2;nn=$3;if(n>.99*nn){if (k > 10) k = 10 ; if (k < 3) k = 3 ; print k;exit;}}'`
  echo -n '// ' > tmp/OR/$group/d4.candidate_introns.ace
  date >>  tmp/OR/$group/d4.candidate_introns.ace
  foreach target ($Etargets)
    cat  tmp/METADATA/$target.[fr].introns ZZZZZ tmp/OR/$group/d4.de_uno.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){a1=$2+0;a2=$3+0;av[$1,a1+0,a2+0] =1 ; nAv++; next; }}{if ($7=="*")next;if($2+0==0)next;min1=factor*minS;if($5 == "gt_ag" || $5 == "gt_ag")min1=minS;if($4<min1)next;nn++;if(av[$1,0+$2,0+$3]==1){ if($5 == "gt_ag" || $5 == "gt_ag") { inGtAg++;}inAv++;n4in+=$4;}else { if($5 == "gt_ag" || $5 == "gt_ag") { outGtAg++;} outAv++;n4Out+=$4;}}END{if(nn==0)nn=1;if(nAv==0)nAv=1;sp=inAv/nn;sn=inAv/nAv;printf("Ali %s\n-D Candidate_introns uno\n-D Candidate_introns %s\n",group,target);if (0){printf("Candidate_introns %s %d:%d:%d ",target, nn,inGtAg+outGtAg,nn-inGtAg-outGtAg); printf(" In_%s %d:%d:%d ",target, inAv, inGtAg,inAv-inGtAg); printf ("Not_in_%s %d:%d:%d ", target, outAv, outGtAg, outAv - outGtAg); printf (" Specificity %.2f Sensitivity %.2f Known_Support %d New_support %d New_minS %d\n\n", 100*sp, 100*sn, n4in, n4Out, minS) ;}else printf("Candidate_introns %s %d In_%s %d Not_in_%s %d  new_gt_ag %d Specificity %.2f Sensitivity %.2f Known_Support %d New_support %d New_minS %d\n\n", target, inAv+outAv, target, inAv, target, outAv, outGtAg, 100*sp, 100*sn, n4in, n4Out, minS);}' target=$target group=$group minS=$minS factor=10 >> tmp/OR/$group/d4.candidate_introns.ace
  end
  cat tmp/OR/$MAGIC.any.introns ZZZZZ tmp/OR/$group/d4.de_uno.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){a1=$2+0;a2=$3+0;av[$1,a1+0,a2+0] =1 ; nAv++; next; }}{if ($7=="*")next;if($2+0==0)next;min1=factor*minS;if($5 == "gt_ag" || $5 == "gt_ag")min1=minS;if($4<min1)next;nn++;if(av[$1,0+$2,0+$3]==1){ if($5 == "gt_ag" || $5 == "gt_ag") { inGtAg++;}inAv++;n4in+=$4;}else { if($5 == "gt_ag" || $5 == "gt_ag") { outGtAg++;} outAv++;n4Out+=$4;}}END{if(nn==0)nn=1;if(nAv==0)nAv=1;sp=inAv/nn;sn=inAv/nAv;printf("Ali %s\n-D Candidate_introns uno\n-D Candidate_introns %s\n",group,target);if (0){printf("Candidate_introns %s %d:%d:%d ",target, nn,inGtAg+outGtAg,nn-inGtAg-outGtAg); printf(" In_%s %d:%d:%d ",target, inAv, inGtAg,inAv-inGtAg); printf ("Not_in_%s %d:%d:%d ", target, outAv, outGtAg, outAv - outGtAg); printf (" Specificity %.2f Sensitivity %.2f Known_Support %d New_support %d New_minS %d\n\n", 100*sp, 100*sn, n4in, n4Out, minS) ;}else printf("Candidate_introns %s %d In_%s %d Not_in_%s %d  new_gt_ag %d Specificity %.2f Sensitivity %.2f Known_Support %d New_support %d New_minS %d\n\n", target, inAv+outAv, target, inAv, target, outAv, outGtAg, 100*sp, 100*sn, n4in, n4Out, minS);}' target=Any group=$group minS=$minS factor=10 >> tmp/OR/$group/d4.candidate_introns.ace
 
endif


if (-e tmp/OR/$group/d4.de_uno.txt && ! -e tmp/OR/$group/d4.de_uno.txt.gz) then
  gzip tmp/OR/$group/d4.de_uno.txt
endif

exit 0
