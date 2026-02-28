function parseDate(d) {
    n = split(d, aa, "/") ;
    if (n == 3)
    {
	m= aa[1] + 0 ; d = aa[2] + 0 ; y = aa[3] + 0 ;
	if (y > 1000 && m < 13 && d < 32)
	    return  y "-" m "-" d ;
    }
    else
    {
	split (d, aa, "-") ; an=aa[3]+0;if(an <50)an += 2000 ; else an += 1900 ;
	m = index("xxJanFebMarAprMayJunJulAugSepOctNovDec",aa[2])/3;if (m<10) m = "0" m ;
	j=aa[1];;if (j+0<10)j="0" j ;
	return an "-" m "-" j ;
    }
}

{    gsub (/\"/,"",$0) ; }

{
    if ($13 == "CTS")  # concatenated tag sequencing
	next ;
    if ($13 == "WCS")  # whole chromosome sequencing
	next ;
    if ($13 == "MeDIP-Seq") # methylation immuno precipitation
	next ;
}

/^Run/{next;}
{ srr = $1 ;
    if (! srr) next ; 
  # use the SRR as the temporary Run name, it will be aliased later
    print "SRR \"" srr "\"" ;
    print "-D Nucleic_Acid_Extraction" ;
    print "-D sraNucleic_Acid_Extraction" ;
    # The SRA file gives 2 different dates, reales and Load,  Load is meaningless
    doOk = 0 ;
    if ($3) { d1 = parseDate($3) ; printf ("Submission_date %s\n", d1) ; }
    if ($2) { d1 = parseDate($2) ; printf ("Release_date %s\n", d1) ; }

    # spots if the number of fragments (to be verified after download in Ali)
    printf ("Spots %d bases_in_SRA %s Average_length %d", $4, $5, $7) ;
    if ($6 + $17 + $18 + 0 > 0) printf (" Insert_size %d ", 0+$17) ;
    # if ($6 + $18 + 0 > 0) printf ("\"+-\" %d", 0+$18) ;
    if ($6 > 0) printf (" spots_with_mate %d", $6) ;
    printf ("\n") ;

    if ($6 > 0) printf ("Paired_end\n") ;

    if ($6 > 0 && $4 != $6) printf ("ERROR \"spots  and paired end differ %d %s\"\n", $4,$6) ;
   
    # download path, report  ERROR if missing 
    if ($10) 
    {  
	z = $10 ;
	gsub (/http:\/\//,"http:\\/\\/",z) ;
	gsub (/https:\/\//,"https:\\/\\/",z) ;
	printf ("SRR_download \"%s\"\n",  z)
    }
    else
	printf ("ERROR \"missing download address\n") ;
     
    ####### column 11 SRX identifier
    srx = $11 ;
    if (srx)  print "SRX " srx ;
     
    ####### column 12 Library_name
    lib = $12 ;
    if (lib == "unspecified") lib = "" ;
    if (length(lib) > 1) print "Library_name \"" lib "\"" ;
    
    ####### column 13, 14, 15 Sequencing strategy , 16 drop, 17, 18 spots, 19, 20 platform
    print "-D RNA" ;
    if ($13 == "WGS")
    {
	if (1) # $15 == "GENOMIC" || $15 == "METAGENOMIC" || $15 == " TRANSCRIPTOMIC" || $15 == "METATRANSCRIPTOMIC" || $15 == "VIRAL"
	{
	    if ($14 == "Hybrid Selection" )
		printf ("sraExome") ;
	    else
	    {
		z = $14 ;
		if (z == "other" || z == "unspecified") $14 = "" ;
		if (z == "MDA") $14 = "MDA often for single cell sequencing" ;
		
		printf ("Whole_genome") ;
	    }

	    z = "" ;
	    if (length($13) > 0 && length($14) > 0) z = $13 ", " $14 ;
	    else if (length($13) > 0) z = $13 ;
	    else if (length($14) > 0) z = $14 ;

	    if (length (z) > 0)
		printf (" \"%s\"", z) ;
	    printf ("\n") ;

	    if ($15 == "METAGENOMIC") 
		printf ("Microbiome\n") ;
	}
	else
	    printf ("ERROR \"This WGS is neither GENOMIC nor METAGENOMIC\n") ;
    }

    else if ($13 == "WXS") #  && ! $14 == "cDNA" && ! $14 == "RT-PCR" && ! $14 == "polyA" && ! $14 == "Oligo-dT")
    {
	printf ("sraExome") ;
	z = $14 ;
	if (z == "other" || z == "unspecified") $14 = "" ;
	if (z == "MDA") $14 = "MDA often for single cell sequencing" ;

	    z = "" ;
	    if (length($13) > 0 && length($14) > 0) z = $13 ", " $14 ;
	    else if (length($13) > 0) z = $13 ;
	    else if (length($14) > 0) z = $14 ;

	    if (length (z) > 0)
		printf (" \"%s\"", z) ;
	    printf ("\n") ;

	if ($15 == "METAGENOMIC")  
	    printf ("sraMicrobiome\n", z) ;
    }
    
    else  # RNA-seq
    {
	if ($15 == "TRANSCRIPTOMIC" || $15 == "METATRANSCRIPTOMIC" || $15 == "VIRAL")
	{
	    if ($15 == "METATRANSCRIPTOMIC" || $15 == "VIRAL")
		printf ("sraMicrobiome\n") ;
	    else  # CLONE CLONEEND EST FL-cDNA other RNA_Seq ncRNA_Seq
		printf ("sraRNA\n") ;
	    
	    if ($14 == "CAGE" || $14 == "RACE")
		printf ("sraRNA sraCap_CAGE_RACE") ;
	    else if ($14 == "cDNA")
		printf ("sraRNA sraUnspecified_RNA") ;
	    else if ($14 == "ChIP")
		printf ("sraRNA sraUnspecified_RNA") ;
	    else if ($14 == "DNase")
		printf ("sraRNA sraUnspecified_RNA") ;
	    else if ($14 == "Inverse rRNA")
		printf ("sraRNA sraUnspecified_RNA") ;
	    else if ($14 == "MDA")
		printf ("sraRNA sraUnspecified_RNA") ;
	    else if ($14 == "PolyA" || $14 == "Oligo-dT")
	    { printf ("sraRNA sraPolyA") ; $14 = "" ; }
	    else if ($14 == "other")
		printf ("sraRNA sraUnspecified_RNA") ;
	    else if ($14 == "PCR")
		printf ("sraRNA sraUnspecified_RNA") ;
	    else if ($14 == "RANDOM")
		printf ("sraRNA sraUnspecified_RNA") ;
	    else if ($14 == "RANDOM PCR")
		printf ("sraRNA sraUnspecified_RNA") ;
	    else if ($14 == "Reduced Representation")
		printf ("sraRNA sraUnspecified_RNA") ;
	    else if ($7 < 100 && $14 == "size fractionation")
		printf ("sraRNA sraSmall_RNA") ;
	    else if ($7 >= 100 && $14 == "size fractionation")
		printf ("sraRNA sraUnspecified_RNA") ;
	    else if ($13 == "RNA-Seq")
		printf ("sraRNA sraUnspecified_RNA") ;
	    else if ($13 == "miRNA-Seq")
		printf ("sraRNA sraSmall_RNA") ;
	    else if ($14 == "Hybrid Selection" || $14 == "RT-PCR" )
		printf ("sraRNA sraGene_selection") ;
	    else if ($14 == "unspecified")
		printf ("sraRNA sraUnspecified_RNA") ;
	    else
		printf ("sraRNA sraUnspecified_RNA") ;
	    
	    if ($13 == "OTHER") $13 = "" ;
	    if ($14 == "other" || $14 == "unspecified") $14 = "" ;
			if ($14 == "MDA") $14 = "MDA often for single cell sequencing" ;
	    z = "" ;
	    if (length($13) > 0 && length($14) > 0) z = $13 ", " $14 ;
	    else if (length($13) > 0) z = $13 ;
	    else if (length($14) > 0) z = $14 ;

	    if (length (z) > 0)
		printf (" \"%s\"", z) ;
	    printf ("\n") ;

	    if ($14 == "Reduced Representation")
		printf ("sraNormalized \"Reduced sepresentation\"\n") ;

             ####### column 13 special tags to be added 

	    if ($13 == "AMPLICON")
		printf ("-D sraUnspecified_RNA\nsraGene_selection\n") ;
	    else if ($13 == "RIP-Seq")
		printf ("-D sraUnspecified_RNA\nsraRIP_CLIP \"%s\"\n", z) ;
	}
    }

    ####### column 16 drop
    ####### column 17 18 insert size, incorporated in Spots
    
    ####### column 19-20 platform-model i.e. Illumina Hiseq
    z = "" ;
    if ($19 == "ABI_SOLID")
	z = "SOLiD" ;
    else if ($19 == "HELICOS")
	z = "Helicos" ;
    else if ($19 == "ILLUMINA")
	z = "Illumina" ;
    else if ($19 == "ULTIMA")
	z = "Ultima" ;
    else if ($19 == "DNBSEQ")
	z = "DNBSEQ" ;
    else if ($19 == "BGISEQ")
	z = "Illumina" ;
    else if ($19 == "ION_TORRENT")
	z = "Ion_Torrent" ;
    else if ($19 == "LS454")
	z = "Roche_454" ;
    else if ($19 == "COMPLETE_GENOMICS")
	z = "Complete_genomics" ;
    else if ($19 == "CAPILLARY")
	z = "Capillary" ;
    else if ($19 == "PACBIO_SMRT")
	z = "PacBio" ;
    else if ($19 == "OXFORD_NANOPORE")
	z = "Nanopore" ;
    else
	z = "ERROR \"unkown platform " $19 "\"" ;
    z1 = $20 ;
    if (z1 == "unspecified") $z1 = "" ;
    if (z && z1) gsub (z, "", z1) ;
    while (substr (z1,1,1) == " ") z1 = substr (z1, 2) ; # remove leading blanks 
    if (length(z) > 0 && length(z1) > 0)
	z = z " \"" z1 "\"" ;
    if (length(z) > 0)
	printf ("%s\n", z) ;
    
    ####### column 21 Project (obsolete as of 2024_03 but needed as a safegard)
    srp = $21 ;  # SRP number
    if ($22=="") printf ("-D Biosample\nSRP %s\n", srp) ;

    ####### column 22 the newer PRJNA number
    prj = $22 ;  # PRJNA number
    if (prj) printf ("-D Biosample\n-D SRP\nSRP %s\n", prj) ;

    ####### column 26 biosample
    biosample = $26 ;
    if (biosample) printf ("Biosample %s\n", biosample) ;

    ####### column 28/29 taxon
    taxName = $29 ;
    taxid = $28 ;
    if (taxid && taxid != 32644) # unidentified
    {
	printf ("Species \"") ;
	if (taxName)
	    printf ("%s (taxid:%s)", taxName, taxid) ;
	else
	    printf ("taxid:%s", taxid) ;
	printf ("\"\n") ;
    }

    ####### column 30 sample_name
    print "-D Sample_name" ;
    z = $30 ;
    if (index (z, "E-MTAB-") == 1)
    {
	j = index (z, ":") ;
	if (j > 1)
	    z = substr (z, j+1) ;
    }
    if (( substr(z,1,3) == "DRS" ||  substr(z,1,3) == "SRS" ||  substr(z,1,3) == "GSM")  &&  0 + substr(z,4) > 1)
	z = "" ;
    if (( substr(z,1,4) == "SAMD" ||  substr(z,1,4) == "SAME" ||  substr(z,1,4) == "SAMN")  &&  0 + substr(z,5) > 1)
	z = "" ;
			
				if (z == "na") z = "" ;
				if (z == "NA") z = "" ;
				if (z == "N/A") z = "" ;
				if (z == "n/a") z = "" ;
				if (z == "missing") z = "" ;
				if (z == "not applicable") z = "" ;
				if (z == "non provided") z = "" ;
				if (z == "no description") z = "" ;
				if (z == "unspecified") z = "" ;

    if (length (z) > 1) print "Sample_name \"" z "\"" ;

    ####### column 35 sex
    z = $35 ;
    if (z == "XX" || z == "F" || z == "f" || z == "Female" || z == "female" )
	print "Female" ;
    if (z == "XY" || z == "M" || z == "m" || z == "Male" || z == "male" )
	print "Male" ;

    ####### column 37 tumor
    z = $37 ;
    if (z == "yes") 
	print "Phenotype Tumor" ;
    # if (z == "no") 
    #    print "Phenotype No_tumor" ;
    
    ####### column 41 Body_site
    z = $41 ;
     if (length (z) > 1) print "Body_site \"" z "\"" ;

    ####### column 42 Center
     z = $42 ;
     if (z == "GEO") z = "" ;
     if (length (z) > 1) print "Center_name \"" z "\"" ;

    ####### column 44 dbgap
    z = $44 ;
     if (length (z) > 1) print "dbGAP \"" z "\"" ;

    ####### done
    printf ("\n") ;
}


