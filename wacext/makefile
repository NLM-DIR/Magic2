# @(#)makefile	1.58 10/1/97

SHELL = /bin/sh
#  MAKE = /usr/local/bin/gmake
#  MAKE = /usr/ccs/bin/make on solaris 

#.SILENT:

# Suppress auto SCCS extraction.
.SCCS_GET:

####  Program list 

PGM = gene2gene gene2chrom2 gene2chrom3 go2gene bly2hosac kscount skscount asndump sasndump gene_nm_remap sgene_nm_remap biblio_as_dna sbiblio_as_dna resizegenebox sresizegenebox abstractchi2 sabstractchi2 e8 su21 su31 irrep sgenecardexport genecardexport gintegra alibaba alibabaclient blystats laure timehisto geturl dnachip dnachipclient akg sakg buildstats sbuildstats worf resequence  maggie omim2ace maqc maqccl av clipalign clipalign_chan paraclipalign estscan UCSCtrackCompare ts mute genexml2ace synonyms genetiling geneelements nico_worf david boubou taxtree UMI  variant_caller tricoteur tsnp jumpalign dnamotif wiggle_smoother samclustering  dna2dna mir wiggle histo snp bestali mega3psort2ace calireport googletable geneindex pg2pg seqc dnastitch kmercount g2x qcsummary dnabloom chanclient chanserver chanscript chantest snpbloom cnv sv sra_metadata demazure dnawalk pi_channel_test introns table_histo magicx sam2gold tfali postMagicBlast tsf altintrons snpsummary 

#### EXTERNAL acedb programs linked against acelib

all $(PGM) : bin-dir copy-makefile.acc
	echo 'making' $@ 'for $(ACEDB_MACHINE)' ; echo
	cd ../bin.$(ACEDB_MACHINE) ; \
	$(MAKE) -f makefile.acc $@

copy-makefile.acc:
	cp makefile.acc ../bin.$(ACEDB_MACHINE) ; \
	cp ../wmake/$(ACEDB_MACHINE)_DEF ../bin.$(ACEDB_MACHINE)/deffile ; \
	chmod u+w ../bin.$(ACEDB_MACHINE)/makefile.acc

#### binary directory bin.$(ACEDB_MACHINE)

bin-dir:
	if test $(ACEDB_MACHINE) ; then exit ; \
	else \
	echo  'To compile, please setenv ACEDB_MACHINE (now not set) to one of: ' ; \
	cd ../wmake ; echo ; ls *DEF | sed 's/_DEF//' ; \
	echo     '      e.g.  add in your .login file the line '; \
	echo  ' setenv ACEDB_MACHINE  SUN' ; \
	echo '  or some other machine name defined in the directory wmake.' ;\
	echo '     Thank you' ; echo ; echo ; exit 1 ; fi
	if test -f ../wmake/$(ACEDB_MACHINE)_DEF ; then exit ; \
	else echo  'To compile, please setenv ACEDB_MACHINE (now = $(ACEDB_MACHINE) ) to one of: ' ; \
	cd ../wmake ; echo ; ls *DEF | sed 's/_DEF//' ; \
	echo     '      i.e.  add in your .login file the line '; \
	echo  ' setenv ACEDB_MACHINE  SUN' ; \
	echo '  or some other machine name defined in the directory wmake.' ;\
	echo '     Thank you' ; echo ; echo ; exit 1 ; fi
	if test -d ../bin.$(ACEDB_MACHINE)  ; then exit ; \
	else mkdir ../bin.$(ACEDB_MACHINE)  ; fi

#### end of makefile

 
 
 
 
