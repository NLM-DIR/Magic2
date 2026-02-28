BEGIN { 
  nt = 0 ; nr=0;
  m = "any" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["any"] = it ; i2nam[nt] = "aligned" ;  
  m = "B_rrna" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["rrna"] = it ; i2nam[nt] = "rRNA" ; 
  m = "A_mito" ; 
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["mito"] = it ; i2nam[nt] = "mitochondria" ; TT1 = it ;
  m = "C_chloro" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["chloro"] = it ; i2nam[nt] = "chloroplast" ;
  m = "DT_magic" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["seqc"] = it ; i2nam[nt] = "Magic_reannotation" ; ok[it] = 1 ;
  m = "ET_av" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["av"] = it ; i2nam[nt] = "AceView_2011" ;
  m = "FT_Eco" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["av"] = it ; i2nam[nt] = "Ecocyc" ;
  m = "FT_Tjaden" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["av"] = it ; i2nam[nt] = "Tjaden" ;
  m = "FT_av2008" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["av"] = it ; i2nam[nt] = "AceView_2008" ;
  m = "FT_extra" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["av"] = it ; i2nam[nt] = "extra" ;
  m = "FT_simul" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["av"] = it ; i2nam[nt] = "Simul" ;
  m = "FT_cloud" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["av"] = it ; i2nam[nt] = "Cloud" ;
  m = "KT_RefSeq" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["RefSeq"] = it  ; i2nam[nt] = "RefSeq_37.104_models" ;
  m = "KT_RefSeqCurrent" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["RefSeq"] = it  ; i2nam[nt] = "RefSeq_37.104_models" ;
  m = "LT_seqc" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["seqc"] = it ; i2nam[nt] = "SEQC_Transcriptome_2013" ;
  m = "LT_magic" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["seqc"] = it ; i2nam[nt] = "magic reconstruction" ;
  m = "LT_UCSC" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ;  tt2i["UCSC"] = it  ; i2nam[nt] = "ucscGene_Nov_2011" ;
  m = "MT_EBI" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ;  tt2i["EBI"] = it  ; i2nam[nt] = "Encode.37.70_Jan2013" ;
  m = "NT_MiT" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "MichiganFeb2015" ;
  m = "NT_miRNA" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "MichiganFeb2015" ;
  m = "NT_HINV" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "HINV" ;
  m = "NT_FBK" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "FBK" ;
  m = "NT_FlyBase" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; i2nam[nt] = "FlyBase" ;
  m = "OT_rnaGene" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i[ "rnaGene"] = it ; i2nam[nt] = "RNA_genes" ;
  m = "PT_tRNA" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["tRNA"] = it  ; i2nam[nt] = "tRNA" ; TT2 = it ;
  m = "QT_smallRNA" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["smallRNA"] = it  ; i2nam[nt] = "smallRNA" ; TT2 = it ;
  m = "Any_transcriptome" ;  
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["Any_transcriptome"] = it  ; i2nam[nt] = "Any previously annotated transcript" ;
  m = "T_GeneBox_RefSeq" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ;  tt2i["GeneBox_RefSeq"] = it  ; i2nam[nt] = "RefSeq_37.104_genes" ;
  m = "T_GeneBox_av" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ;  tt2i["GeneBox_av"] = it  ; i2nam[nt] = "AceView_2011_genes" ;
  m = "T_GeneBox_EBI" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ;  tt2i["GeneBox_EBI"] = it  ; i2nam[nt] = "Encode.37.70_Jan2013_genes" ;
  m = "T_GeneBox_MiT" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ;  tt2i["GeneBox_MiT"] = it  ; i2nam[nt] = "MichiganFeb2015" ;
  m = "U_introns" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ;  tt2i["introns"] = it  ; i2nam[nt] = "introns" ;
  m = "0_SpikeIn" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["SpikeIn"] = it ; i2nam[nt] = "ERCC_SpikeIn" ;
  m = "1_DNASpikeIn" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["DNASpikeIn"] = it ; i2nam[nt] = "DNA_SpikeIn" ;
  m = "X_Bamy" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["Bamy"] = it ; i2nam[nt] = "Bacillus amyloliquefaciens" ; ok[it] = 1 ;
  m = "W_Line" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["Line"] = it ; i2nam[nt] = "LINE repeats" ; ok[it] = 1 ;
  m = "Y_Pfluo" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["Pfluo"] = it ; i2nam[nt] = "Pseudomonas fluorescens" ; ok[it] = 1 ;
  m = "Z_genome" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["genome"] = it  ; i2nam[nt] = "Genome" ;
  m = "z_gdecoy" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["gdecoy"] = it ; i2nam[nt] = "Imaginary_genome_specificity_control" ;
  m = "b_bacteria" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["bacteria"] = it ; i2nam[nt] = "bacterial collection"; ok[it] = 1 ;
  m = "v_virus" ;
  t = m ;  nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ; tt2i["virus"] = it ; i2nam[nt] = "virus collection"; ok[it] = 1 ;
}
/^ZZZZZ/{zz++ ; if(zz==2){  RR = nr ; } next;}
{
  gsub(/\\\//,"/",$0);
}
{  # select the good targets
  if (0+zz < 1) 
    {
      it = tt2i[$1] ; ok[it] = 1 ; 
      next ;
    }  
}
{
  if (zz < 3) 
    {
      z = $1 ; ir = r2i[z] ; 
      if (ir < 1) 
	{ 
	  nr++ ; ir = nr ; i2r[nr] = z ; r2i[z] = nr ;
	}
      if (zz < 2) 
	{ isgroup[ir] = 1 ; ngroups++ ; }
      else 
	isgroup[ir] = 0 ;
      good[ir] = 1 ;
      next ;
    }
}
/^\"/{ gsub("\"","",$0);split($0,aa,"\t");rg=aa[1];machine[rg]=aa[2];sample[rg]=aa[3];systm[rg]=aa[4];tissue[rg]=aa[5];titre[rg]=aa[6];stitre[rg]=aa[11];stitre2[rg]=aa[12];otitre[rg]=aa[13];systm[rg] = systm[rg] " " aa[7];runid[rg]=aa[8];sytitre[rg]=aa[10];next;}
{
  gsub("\"","",$0);
}
/^Run/{ z=$2; ir = r2i[z]; next ; }
/^Ali/{ z=$2; ir = r2i[z]; next ; }
/^Accepted/ { seq[ir] = $2 ; tag[ir] = $4 ; kb[ir] = $6 ; next; }
/^Exit_adaptor_clipping/ { vecx[ir] = $2 ; vecxkb[ir] += $4 ; next; }
/^Entry_adaptor_clipping/ { vece[ir] = $2 ; vecekb[ir] += $4 ; next; }

/^h_Ali/{ t = $2; it = 0 + t2i[t] ; if(it == 0) { nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ;} tt[it]=1; tseq[ir,it] = $3 ; ttag[ir,it] = $5 ; tkb[ir,it] = $7 ; tbp[ir,it] = $9 ; next ; }
/^nh_Ali/{ t = $2; it = 0 + t2i[t] ; if(it == 0) { nt++; it = nt ;t2i[t] = nt ; i2t[nt]=t ;} tt[it]=1; tseq[ir,it] = $3 ; ttag[ir,it] = $5 ; tkb[ir,it] = $7 ; tbp[ir,it] = $9 ; next ; }
/^ATGC_kb/ {gc[ir]=100*($9+$10)/($7+$8+$9+$10+1);}
/^CPU/{cpu[ir]+=$3;}
/^Max_memory/{if($3 > maxmem[ir])maxmem[ir]=$3;}
END {

  if (title == "Hierarchic") 
    {
      it= tt2i["genome"] ; i2nam[it] = "New_exonic_intergenic_or_Intronic" ;
      it2= tt2i["Any_transcriptome"] ;   ok[it2] = 1 ;
     
      for (r = 1 ; r <= nr ; r++) 
	if (good[r]==1)
	  {
	    for (it = TT1 ; it <= TT2 ; it++)
	      {
		ttag[r,it2] +=  ttag[r,it] ;
		tseq[r,it2] +=  tseq[r,it] ;
		tkb[r,it2]  +=  tkb[r,it] ;
	      }
	  }
    }
  printf ("\n") ;
 
  printf ("10000\tGroup or Run\tRunId\tTitle\tSorting_title\tSorting_title2\tOther_title\tSy_title\tReads\tAligned reads\t%% aligned\tread/sequence ratio\taligned read/sequence ratio\t%% GC\tMillion reads aligned  on all targets per CPU-hour\tMax RAM Mbytes\t\tGroup or Run\tTitle\tSorting_title\tSorting_title2\tOther_title") ;
  if (0) printf("\tEntry adaptor\tExit adaptor\t%% Entry adaptor\t%% Exit adaptor") ;
  for (it = 2 ; it <= nt ; it++) { if (ok[it]) printf("\t%s reads", i2nam[it]);}
  printf ("\t\tGroup or Run\tRunId\tTitle\tSorting_title\tSorting_title2\tOther_title\tSy_title") ;
  for (it = 2 ; it <= nt ; it++) { if (ok[it]) printf("\t%% %s reads", i2nam[it]);}
  printf ("\tPool\tTissue\tSample") ;
  ng = 0 ;
  for (r = 1 ; r <= nr ; r++) 
    if (good[r]==1)
      { 
	if (0 && nr >= 24)
	  {
	    split (titre[r], gnam, "_") ;
	    if (ng == 0) 
	      printf("\n10000\t.") ;
	    else
	      {
		if ((gnam[1] != oldgnam1 || (ngroups > 24 && gnam[2] != oldgnam2)) && ! (gnam[3] == "any" && oldgnam3 == "any"))
		  printf("\n10000\t.") ;
	      }
	    oldgnam1 = gnam[1] ; oldgnam2 = gnam[2] ; oldgnam3 = gnam[3] ;
	    ng++ ;
	  }
	printf ("\n%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d", 10000+r, i2r[r], runid[i2r[r]], titre[i2r[r]], stitre[i2r[r]], stitre2[i2r[r]], otitre[i2r[r]], sytitre[i2r[r]], tag[r]) ;
	printf ("\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t\t%s", ttag[r,1], 100*ttag[r,1]/(1+tag[r]), tag[r]/(1+seq[r]), ttag[r,1]/(1+tseq[r,1]), gc[r], (3600*tag[r])/(1000000*cpu[r]+1), maxmem[r], i2r[r]) ;
	printf ("\t%s\t%s\t%s\t%s", titre[i2r[r]], stitre[i2r[r]], stitre2[i2r[r]], otitre[i2r[r]]) ;
	if (0) printf ("\n%d\t%d\t%.2f\t%.2f", vece[r], vecx[r], 100*vece[r]/(1+tag[r]), 100*vecx[r]/(1+tag[r])) ;
	for (it = 2 ; it <= nt ; it++) 
	  if (ok[it])printf ("\t%d", ttag[r,it]) ;
	printf ("\t\t%s\t%s", i2r[r], runid[i2r[r]]) ;
	printf ("\t%s\t%s\t%s\t%s\t%s", titre[i2r[r]], stitre[i2r[r]], stitre2[i2r[r]], otitre[i2r[r]], sytitre[i2r[r]]) ;
	for (it = 2 ; it <= nt ; it++) 
	  {
	    if (0 + ok[it] < 1) continue ;
	    z = tag[r] ;
	    t = "any" ;      jt = t2i[t] ;       any = ttag[r,jt] ;
	    t = "z_gdecoy" ; jt = t2i[t] ;    gdecoy = ttag[r,jt] ;
	    t = "0_SpikeIn" ;   jt = t2i[t] ;      SpikeIn = ttag[r,jt] ;
	    t = "1_DNASpikeIn" ;   jt = t2i[t] ;      DNASpikeIn = ttag[r,jt] ;
	    t = "B_rrna" ;   jt = t2i[t] ;      rrna = ttag[r,jt] ;
	    
	    z = any - gdecoy - SpikeIn - DNASpikeIn - rrna ;

	    if (i2t[it] == "B_rrna") z = any - gdecoy - SpikeIn -DNASpikeIn ;
	    if (i2t[it] == "any") z = tag[r] ;
	    if (i2t[it] == "0_SpikeIn") z = any - gdecoy ;
	    if (i2t[it] == "1_DNASpikeIn") z = any - gdecoy ;
	    if (i2t[it] == "z_gdecoy") z = any ;
	   
            z = ttag[r,1] ; 
	    if (z <= 0) z = 1 ;
	    
	    
	    if (ttag[r,it] * 10000 < any || ttag[r,it] > 2 * z  ||  (2*ttag[r,it] > z &&  ttag[r,it] * 50 < any ))
	      printf ("\t") ;
	    else	    
	      printf ("\t%.2f", 100 * ttag[r,it]/z) ;
	  }
	printf ("\t%s", pool[i2r[r]]) ;
	printf ("\t%s", tissue[i2r[r]]) ;
	printf ("\t%s", sample[i2r[r]]) ;
      }
  
  t = "any" ; it = t2i[t] ; i2nam[it] = "Any_target_kb" ;
    printf ("\n") ;

  printf ("30000\tGroup or Run\tRunId\tTitle\tSorting_title\tSorting_title2\tOther_title\tSy_title\tkb\tkb aligned\t%% aligned\tread/sequence ratio\taligned read/sequence ratio\t%% GC\tGigabases aligned  on all targets per CPU-hour\tMax RAM Mbytes\t\tGroup or Run\tTitle\tSorting_title\tSorting_title2\tOther_title") ;
  if (0) printf("\tEntry adaptor\tExit adaptor\t%% Entry adaptor\t%% Exit adaptor") ;
  for (it = 2 ; it <= nt ; it++) { if (ok[it]) printf("\t%s kb", i2nam[it]);}
  printf ("\t\tGroup or Run\tRunId\tTitle\tSorting_title\tSorting_title2\tOther_title\tSy_title") ;
  for (it = 2 ; it <= nt ; it++) { if (ok[it]) printf("\t%% %s kb", i2nam[it]);}
  printf ("\tPool\tTissue\tSample") ;
  ng = 0 ;
  for (r = 1 ; r <= nr ; r++) 
    if (good[r]==1)
      { 
	if (0 && nr >= 24)
	  {
	    split (titre[r], gnam, "_") ;
	    if (ng == 0) 
	      printf("\n30000") ;
	    else
	      {
		if ((gnam[1] != oldgnam1 || (ngroups > 24 && gnam[2] != oldgnam2)) && ! (gnam[3] == "any" && oldgnam3 == "any"))
		  printf("\n30000") ;
	      }
	    oldgnam1 = gnam[1] ; oldgnam2 = gnam[2] ; oldgnam3 = gnam[3] ;
	    ng++ ;
	  }

	printf ("\n%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d", 30000+r, i2r[r], runid[i2r[r]], titre[i2r[r]], stitre[i2r[r]], stitre2[i2r[r]], otitre[i2r[r]], sytitre[i2r[r]], kb[r]) ;
	printf ("\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t\t%s", tkb[r,1], 100*tkb[r,1]/(1+kb[r]), tag[r]/(1+seq[r]), ttag[r,1]/(1+tseq[r,1]), gc[r],(3600*kb[r])/(1000000*cpu[r]+1), maxmem[r] , i2r[r]) ;
	printf ("\t%s\t%s\t%s\t%s", titre[i2r[r]], stitre[i2r[r]], stitre2[i2r[r]], otitre[i2r[r]]) ;
	if (0) printf ("\t%d\t%d\t%.2f\t%.2f", vecekb[r], vecxkb[r], 100*vecekb[r]/(1+kb[r]), 100*vecxkb[r]/(1+kb[r])) ;
	for (it = 2 ; it <= nt ; it++) 
	  if (ok[it]) printf ("\t%d", tkb[r,it]) ;
	printf ("\t\t%s\t%s", i2r[r], runid[i2r[r]]) ;
	printf ("\t%s\t%s\t%s\t%s\t%s", titre[i2r[r]], stitre[i2r[r]], stitre2[i2r[r]], otitre[i2r[r]], sytitre[i2r[r]]) ;
	for (it = 2 ; it <= nt ; it++) 
	  {
	    if (0 + ok[it] < 1) continue ;
	    z = kb[r] ; 
	    
	    t = "any" ;      jt = t2i[t] ;       any = tkb[r,jt] ;
	    t = "z_gdecoy" ; jt = t2i[t] ;    gdecoy = tkb[r,jt] ;
	    t = "0_SpikeIn" ;   jt = t2i[t] ;      SpikeIn = tkb[r,jt] ;
	    t = "1_DNASpikeIn" ;   jt = t2i[t] ;      DNASpikeIn = tkb[r,jt] ;
	    t = "B_rrna" ;   jt = t2i[t] ;      rrna = tkb[r,jt] ;
	    
	    z = any - gdecoy - SpikeIn - DNASpikeIn - rrna ;
	    
	    if (i2t[it] == "B_rrna") z = any - gdecoy - SpikeIn - DNASpikeIn ;
	    if (i2t[it] == "any") z = kb[r] ;
	    if (i2t[it] == "0_SpikeIn") z = any - z_gdecoy ;
	    if (i2t[it] == "1_DNASpikeIn") z = any - z_gdecoy ;
	    if (i2t[it] == "z_gdecoy") z = any ;
	    
            z = tkb[r,1] ; 
	    if (z <= 0) z = 1 ;
	    
	    
	    if (tkb[r,it] * 10000 < any || tkb[r,it] > 2* z  ||  (2*tkb[r,it] > z &&  tkb[r,it] * 50 < any ))
	      printf ("\t") ;
	    else	    
	      printf ("\t%.2f", 100 * tkb[r,it]/z) ;
	  }
	printf ("\t%s", pool[i2r[r]]) ;
	printf ("\t%s", tissue[i2r[r]]) ;
	printf ("\t%s", sample[i2r[r]]) ;
      }
  
  printf ("\n") ;
  t = "any" ; it = t2i[t] ; i2nam[it] = "Any_target_sequence" ;

  for (it = 1 ; it <= nt ; it++) 
    if (! i2nam[it]) 
      i2nam[it] = i2t[it] ;
  
  for (r in rr) { for(im=1;im<=nManip[r];im++){m=manip[r,im];nm2r[m]++;}}

  printf ("50000\tGroup or Run\tRunId\tTitle\tSorting_title\tOther_title\tSorting_title2\tSy_title\tSequences\tAligned sequences\t%% aligned\tread/sequence ratio\taligned read/sequence ratio\t%% GC\tMillion distinct sequences aligned on all targets per CPU-hour\tMax RAM Mbytes\t\tGroup or Run\tTitle\tSorting_title\tSorting_title2\tOther_title") ;

  for (it = 2 ; it <= nt ; it++) { if (ok[it]) printf("\t%s sequences", i2nam[it]);}
  printf ("\t\tGroup or Run\tRunId\tTitle\tSorting_title\tSorting_title2\tOther_title\tSy_title") ;
  for (it = 2 ; it <= nt ; it++) { if (ok[it]) printf("\t%% %s sequences", i2nam[it]);}
  printf ("\tPool\tTissue\tSample") ;
  ng = 0 ;
  for (r = 1 ; r <= nr ; r++) 
    if (good[r]==1)
      { 
	if (0 && nr >= 24)
	  {
	    split (titre[r], gnam, "_") ;
	    if (ng == 0) 
	      printf("\n50000") ;
	    else
	      {
		if ((gnam[1] != oldgnam1 || (ngroups > 24 && gnam[2] != oldgnam2)) && ! (gnam[3] == "any" && oldgnam3 == "any"))
		  printf("\n50000") ;
	      }
	    oldgnam1 = gnam[1] ; oldgnam2 = gnam[2] ; oldgnam3 = gnam[3] ;
	    ng++ ;
	  }

	printf ("\n%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d", 50000+r, i2r[r], runid[i2r[r]], titre[i2r[r]], stitre[i2r[r]], stitre2[i2r[r]], otitre[i2r[r]], sytitre[i2r[r]], seq[r]) ;
	printf ("\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t\t%s", tseq[r,1], 100*tseq[r,1]/(1+seq[r]), tag[r]/(1+seq[r]), ttag[r,1]/(1+tseq[r,1]), gc[r], (3600*seq[r])/(1000000*cpu[r]+1), maxmem[r], i2r[r]) ;
	printf ("\t%s\t%s\t%s\t%s", titre[i2r[r]], stitre[i2r[r]], stitre2[i2r[r]], otitre[i2r[r]]) ;
	for (it = 2 ; it <= nt ; it++) 
	  if (ok[it]) printf ("\t%d", tseq[r,it]) ;
	printf ("\t\t%s\t%s", i2r[r], runid[i2r[r]]) ;
	printf ("\t%s\t%s\t%s\t%s\t%s", titre[i2r[r]], stitre[i2r[r]], stitre2[i2r[r]], otitre[i2r[r]], sytitre[i2r[r]]) ;
	for (it = 2 ; it <= nt ; it++) 
	  {
	    if (0 + ok[it] < 1) continue ;
	    z = seq[r] ;
	    t = "any" ;      jt = t2i[t] ;       any = tseq[r,jt] ;
	    t = "z_gdecoy" ; jt = t2i[t] ;    gdecoy = tseq[r,jt] ;
	    t = "0_SpikeIn" ;   jt = t2i[t] ;      SpikeIn = tseq[r,jt] ;
	    t = "1_DNASpikeIn" ;   jt = t2i[t] ;      DNASpikeIn = tseq[r,jt] ;
	    t = "B_rrna" ;   jt = t2i[t] ;      rrna = tseq[r,jt] ;
	    
	    z = any - gdecoy - SpikeIn - rrna ;
	    
	    if (i2t[it] == "B_rrna") z = any - gdecoy - SpikeIn - DNASpikeIn ;
	    if (i2t[it] == "any") z = seq[r] ;
	    if (i2t[it] == "0_SpikeIn") z = any - z_gdecoy ;
	    if (i2t[it] == "1_DNASpikeIn") z = any - z_gdecoy ;
	    if (i2t[it] == "z_gdecoy") z = any ;
	    
            z = tseq[r,1] ; 
	    if (z <= 0) z = 1 ;

	    if (tseq[r,it] * 100000000 < any || tseq[r,it] > 2* z ||  (2*tseq[r,it] > z &&  tseq[r,it] * 50 < any ))
	      printf ("\t") ;
	    else if (tseq[r,it] * 1000000 < any)	    
	      printf ("\t%.6f", 100 * tseq[r,it]/z) ;
	    else if (tseq[r,it] * 10000 < any)	    
	      printf ("\t%.4f", 100 * tseq[r,it]/z) ;
	    else	    
	      printf ("\t%.2f", 100 * tseq[r,it]/z) ;
	  }
	printf ("\t%s", pool[i2r[r]]) ;
	printf ("\t%s", tissue[i2r[r]]) ;
	printf ("\t%s", sample[i2r[r]]) ;
      }
  printf ("\n") ;
}

