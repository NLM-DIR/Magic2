BEGIN {
    nam["A1"] = "AGLR1" ; 
    nam["A2"] = "AGLR2" ; 
    nam["I1"] = "ILMR1" ; 
    nam["I2"] = "ILMR2" ; 
    nam["I3"] = "ILMR3" ; 
    nam["R1"] = "ROCR1" ; 
    nam["R2"] = "ROCR2" ; 
    nam["R3"] = "ROCR3" ; 

    xpu = 1 ; xpnu = 2 ; xbu = 3 ; xbnu = 4 ; xcu = 5 ; xcnu = 6 ;
}
{   # discard the comments
    if (substr($2,1,1) == "#") next ; 
}
/^ZZZZZ/ {zz++ ; next ; }
{
    if (zz < 1)
    { # parse the metadata
	g2nam[$1] = $2  ; g2ln[$1] = $3 ; g2chrom[$1] = $4 ; g2a1[$1] = $5 ; g2a2[$1] = $6 ; 
	next ; 
    }

    # parse the data
    pp = $1 ; ppp[pp]++ ;    # probe list
    g = $2 ;  gg[g]++ ;      # gene list
    if ($3+$4 > 0) gHasP[g,pp] = 1 ; 
    anyU[g]+= $3 ; anyNU[g]+= $4 ; anyBU[g]+= $5 ; anyBNU[g]+= $6 ; 

    ppg[pp,g,xpu] = $3 ;     # unique probes
    ppg[pp,g,xpnu] = $4 ;    # unique bases
    ppg[pp,g,xbu] = $5 ;     # non unique probes 
    ppg[pp,g,xbnu] = $6 ;    # non unique bases

    ln = 0+g2ln[g] ; 
    split(g,aa,"(") ; 
    if (ln == 0)
    {
	gsub(")","",aa[2]) ; 
	ln = g2ln[aa[2]]+0 ; 
    }
    if (ln == 0) ln = g2ln[aa[1]]+0 ; 
    if (0) print "GGG",g,aa[1],aa[2],ln ; 
    if (ln>0)
    {
	ppg[pp,g,xcu] = int(.5+100*$5/ln)/100 ; 
	ppg[pp,g,xcnu] = int(.5+100*$6/ln)/100 ; 
    }
}
END {
    split(".8,1,1,.1,.4,.2,.5,.2",mmm,",") ; 
    split(".4,.8,.2,1,1,.5,.2,.1",mmm,",") ; 
    np = split("A2|R2|R3|I3|A1|R1|I1|I2",aa,"|") ; 
    np = split("A1|A2|R1|R2|R3|I1|I2|I3",aa,"|") ; 
    split("probes mapped in a single gene|probes mapped in several genes|probe bases mapped in a single gene|probe bases mapped in several genes|coverage mapped in a single gene|coverage mapped in several genes",tit,"|") ; 
    split("Probes mapped in a single gene|Probe bases mapped in a single gene|Fold coverage by probes mapped only to this gene|Probes mapped in several genes|Probe bases mapped in several genes|Fold coverage by probes mapped to several genes",tit,"|") ; 
    
    printf("# Gene\tGene ID\tChromosome\tfrom\tto\tLength of the longest transcript (nt)\tGene lightly touched\tGene well covered by capture probes in A:Agilent, R:Roche.I:Illumina\tProbes mapped only in this gene (any platform)\tProbes mapped in several genes (any platform)\tPercentage of probes mapped only in this gene\tProbe bases mapped only in this gene\tProbe bases mapped in several genes") ; 
    
    for (p = 1 ; p <= np ; p++)
    { 
	printf("\t%s",nam[aa[p]]) ;
	for(j = 1 ; j <= 6 ; j++)
	    printf("\t%s %s",nam[aa[p]],tit[j]) ; 
    }
    for(g in gg)
    {
	ok = 0 ;
	for (p = 1 ; p <= np ; p++)
	{
	    pp = aa[p] ; 
	    if (ppg[aa[p],g,xcu]+ppg[aa[p],g,xcnu] >= mmm[p])
		ok = 1 ;
	} 
	if (ok == 0)
	    continue ;
	g2 = g ; 
	k = split(g,kk,"(") ; 
	if (k>1)
	{
	    split(kk[2],k2k,")") ; g2 = k2k[1] ; 
	}
	g3 = g ; ln = g2ln[g]+0 ; 
	if (ln == 0)g3 = g2 ; 
	
	printf("\n%s\t%s\t%s\t%d\t%d\t%d\t",g,g2nam[g2],g2chrom[g2],g2a1[g3],g2a2[g3],g2ln[g3]) ; 
	
	for (p = 1 ; p <= np ; p++)
	{
	    pp = aa[p] ; 
	    if (gHasP[g,pp] == 1 && (ppg[aa[p],g,xcu]+ppg[aa[p],g,xcnu] < mmm[p]))
		printf("%s",pp) ;
	} 
	printf("\t") ; 
	for (p = 1 ; p <= np ; p++)
	{
	    pp = aa[p] ; 
	    if (ppg[aa[p],g,xcu]+ppg[aa[p],g,xcnu] >= mmm[p])
		printf("%s",pp) ;
	} 
	printf("\t%d\t%d\t%.2f\t%d\t%d",anyU[g],anyNU[g],100.0*anyU[g]/(0.0001+anyU[g]+anyNU[g]),anyBU[g],anyBNU[g]) ;  
	for (p = 1 ; p <= np ; p++)
	{
	    printf ("\t") ;
	    printf("\t%.2f",ppg[aa[p],g,xpu]) ; 
	    printf("\t%.2f",ppg[aa[p],g,xbu]) ; 
	    printf("\t%.2f",ppg[aa[p],g,xcu]) ; 
	    printf("\t%.2f",ppg[aa[p],g,xpnu])		    
	    printf("\t%.2f",ppg[aa[p],g,xbnu]) ; 
	    printf("\t%.2f",ppg[aa[p],g,xcnu]) ; 
	}

    }
}
