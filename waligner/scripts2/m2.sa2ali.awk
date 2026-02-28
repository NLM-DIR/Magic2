BEGIN {
    FS = "\t" ;
}
{
    run = $1 ;
    tag =$2 ;
    nn=$4 ;
    #print ;
    tags1 = "Compatible_pairs__Non_compatible_pairs__Circle_pairs__Perfect_reads" ;
    n = split(tags1,tags,"__") ;
    for (i = 1 ; i <= n ; i++)
	if (tag == tags[i])
	{
	    if (nn +0 > 0) printf ("%s %d\n", tag, nn) ;
	    next;
	}
        if (tag == "Reads_too_short")
    if (tag == "Aligned_pairs")
    {
	printf ("%s %d\n", "Aligned_fragments", nn) ;
	next ;
    }
    if (tag == "Low_entropy_reads")
    {
	lowEntropy = nn ;
	next ;
    }
    if (tag == "RawReads")
    {
	rawReads = nn ;
	next ;
    }
    if (tag == "Reads")
    {
	rawReads = nn ;
	reads = nn ;
	next ;
    }
    if (tag == "Aligned_reads")
    {
	alignedReads = nn ;
	next ;
    }
    if (tag == "Reads_Aligned_once")
    {
	multiAli[1] = nn ;
	allMultiAli += nn ;
	next ;
    }
    if (substr(tag, 1, 21) == "Reads_multi_aligned__")
    {
	k = 0 + substr (tag, 22) ;
	if (k > 1 && k <= 10)
	    multiAli[k] = nn ;
	allMultiAli += nn ;
	next ;
    }
    if (tag == "Reads_aligned_in_class_G")
    {
	alignedReads_G = nn ;
	next ;
    }
    if (tag == "Reads_aligned_in_class_M")
    {
	alignedReads_M = nn ;
	next ;
    }
    if (tag == "Reads_aligned_in_class_R")
    {
	alignedReads_R = nn ;
	next ;
    }
    if (tag == "Reads_aligned_in_class_E")
    {
	alignedReads_E = nn ;
	next ;
    }
    if (tag == "Bases_aligned_in_class_G")
    {
	alignedBases_G = nn ;
	next ;
    }
    if (tag == "Bases_aligned_in_class_M")
    {
	alignedBases_M = nn ;
	next ;
    }
    if (tag == "Bases_aligned_in_class_R")
    {
	alignedBases_R = nn ;
	next ;
    }
    if (tag == "Bases_aligned_in_class_E")
    {
	alignedBases_E = nn ;
	next ;
    }

    if (tag == "Aligned_bases")
    {
	alignedBases = nn ;
	next ;
    }
    if (tag == "Min_read_length")
    {
	printf ("Min_probe_length %d\n", nn) ;
	next ;
    }
    if (tag == "Max_read_length")
    {
	printf ("Max_probe_length %d\n", nn) ;
	next ;
    }
    if (tag == "xxx")
    {
	printf ("%s %d\n", tag, nn) ;
	next ;
    }
    if (tag == "Bases")
    {
	bases = nn ;
	rawBases = nn ;
	next ;
    }
    if (tag == "RawBases")
    {
	rawBases = nn ;
	next ;
    }
    if (tag == "Missmatches")
    {
	printf ("Cumulated_mismatches %d\n", nn) ;
	next ;
    }
    if (tag == "Stranding_in_class_E")
    {
	printf ("Stranding 0_SpikeIn.f2 %s %d plus %d minus\n", $4, $5, $6) 
    }
    if (tag == "Stranding_in_class_R")
    {
	printf ("Stranding B_rrna.f2 %s %d plus %d minus\n", $4, $5, $6) 
    }
    if (tag == "Stranding_in_class_M")
    {
	printf ("Stranding A_mito.f2 %s %d plus %d minus\n", $4, $5, $6) 
    }
    if (tag == "Stranding_in_class_G")
    {
	printf ("Stranding Z_genome.f2 %s %d plus %d minus\n", $4, $5, $6) 
    }
    if (tag == "Intron_supports")
    {
	intronSupport = nn ; iiPlus = $5 ; iiMinus = $6 ;
	printf ("Stranding I_introns %.5f %d plus %d minus\n", iiPlus/(nn+0.00000001), iiPlus, iiMinus) ; 
	next ;
    }
    if (tag == "Supported_introns")
    {
	printf ("Candidate_introns any %d In_any %d Not_in_any 0 New_gt_ag 0 Sensitivity 0 Specificity 0 Known_support %d New_support 0\n", nn, nn, intronSupport) ;
    }

    if (tag == "ATGCN")
    {
	A = $4 ; T = $5 ; G = $6 ; C = $7 ; N= $8 ; t = A + T + G + C + N ;
	if (t > 0) printf ("ATGC_kb %d %d %d %d %d    %d %d %d %d %d\n", 1000*A/t, 1000*T/t, 1000*G/t, 1000*C/t, 1000*N/t,A/1000, T/1000, G/1000, C/1000, N/1000) ;
    }
    if (tag == "Length_distribution_1_5_50_95_99_mode_av")
    {
	printf ("%s %s %s %s %s %s %s %s \n", tag, $4,$5,$6,$7,$8,$9,$10) ;
    }
    if (tag == "Exonic_intronic_intergenic_Bases")
    {
	x = int (($8+0)/1000) ;
	if (x > 0)
	{
	    printf ("Intergenic %d kb\n", x) ;
	    intergenic = $8 ;
	}
	next ;
    }
    if (tag == "genomeLength")
    {
	genomeLength = nn ;
	next ;
    }
    if (tag == "Clipped_3prime_Adaptor_read1")
    {
	adpt = $5 ;
	if (length(adpt)> 2)
	    printf ("Adaptor1 %s %d\n", adpt, $4) ;
    }
    if (tag == "Clipped_3prime_Adaptor_read2")
    {
	adpt = $5 ;
	if (length(adpt)> 2)
	    printf ("Adaptor2 %s %d\n", adpt, $4) ;
    }
    if (tag == "CDS_utr_intronic_intergenic_Bases")
    {
	if ($4 + 0 > 0) printf ("S_1_CDS %.3f \"Mb aligned\"\n", $4) ;
	if ($4 + 0 > 0) printf ("S_1_UTR %.3f \"Mb aligned\"\n", $6) ;
	if ($4 + 0 > 0) printf ("S_1_intronic %.3f \"Mb aligned\"\n", $8) ;
	if ($4 + 0 > 0) printf ("S_1_intergenic %.3f \"Mb aligned\"\n", $10) ;
    }

}
END {
    if (reads + 0 > 0)
    {
	bb = int((bases + reads/2)/reads) ;
	kb = bases/1000 ;
	printf ("Raw_data %d Id %d Accepted %d kb\n", rawReads, reads, kb) ;
	printf ("Accepted %d Seq %d Tags %d kb %d bp\n", rawReads - lowEntropy - tooShort, reads  - lowEntropy - tooShort, kb - lowEntropyBases - tooShortBases, bb) ;
	if (alignedReads_E + 0 > 0)
	    printf ("nh_ali E_ercc %d seq %d tags %d kb_ali %d bp_av_ali %d kb_clip %d bp_av_clip\n", alignedReads_E, alignedReads_E, alignedBases_E/1000, alignedBases_E/alignedReads , alignedBases_E/1000, alignedBases_E/alignedReads_E ) ;
	if (alignedReads_M + 0 > 0)
	    printf ("nh_ali A_mito %d seq %d tags %d kb_ali %d bp_av_ali %d kb_clip %d bp_av_clip\n", alignedReads_M, alignedReads_M, alignedBases_M/1000, alignedBases_M/alignedReads , alignedBases_M/1000, alignedBases_M/alignedReads_M ) ;
	if (alignedReads_R + 0 > 0)
	    printf ("nh_ali B_rrna %d seq %d tags %d kb_ali %d bp_av_ali %d kb_clip %d bp_av_clip\n", alignedReads_R, alignedReads_R, alignedBases_R/1000, alignedBases_R/alignedReads , alignedBases_R/1000, alignedBases_R/alignedReads_R ) ;
	if (alignedReads_G + 0 > 0)
	    printf ("nh_ali Z_genome %d seq %d tags %d kb_ali %d bp_av_ali %d kb_clip %d bp_av_clip\n", alignedReads_G, alignedReads_G, alignedBases_G/1000, alignedBases_G/alignedReads , alignedBases_G/1000, alignedBases_G/alignedReads_G ) ;
	if (alignedReads + 0 > 0)
	    printf ("nh_ali any %d seq %d tags %d kb_ali %d bp_av_ali %d kb_clip %d bp_av_clip\n", alignedReads, alignedReads, alignedBases/1000, alignedBases/alignedReads , alignedBases/1000, alignedBases/alignedReads ) ;

	if (allMultiAli > 0)
	{
	    printf ("Unicity Number_of_targets 1 2 3 4 5 6 7 8 9 10\n") ;
	    printf ("Unicity any ") ;
	    for (k = 1 ; k <= 10 ; k++)
		printf (" %d ", 0 + multiAli[k]) ;
	    printf ("\n") ;
	}
	if (genomeLength > 0 && intergenic > 0)
	{
	    printf ("Intergenic_density %.3f\n", intergenic/genomeLength) ;
	}
    }
    printf ("\n") ;
}
