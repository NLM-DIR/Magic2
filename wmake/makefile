#####################################################################
#### Main makefile for Acedb, contig tools and Image
#####################################################################

# $Id: makefile,v 1.84 2019/05/24 21:47:50 mieg Exp $


# don't change this, using other shells can lead to portability problems.
SHELL = /bin/sh

# for debugging makefile, comment out .SILENT
# .SILENT:

# Suppress auto SCCS extraction.
#.SCCS_GET:

#####################################################################
#### ACEDB programs


#####################################################################
all: libs ace aceview acembly acext acenn other tools misc

# Sub categories of code, please add new stuff to one of these groups if possible

LIBS= libfree.a libtsfree.a libgraph.a libgif.a libace.a libacedna.a libstaden.a  libmd5.a libvsim.a libaccl.a  libacs.a libtcpfn.a libace5.a libwego.a libchannels.a libxace.a



libs: $(LIBS)

ACE= xace tace giface taceserver taceclient
ace: $(ACE)

ACEMBLY= makeSCF xacembly tacembly gifacembly \
	pubid kantormegaparse \
	locuslinkgo cl2acc \
	acekog2ace tgifacemblyserver \
	kantorcount ckantorcount \
	asnbeautifier cbeautifier expasy scf2ctf ctf2scf 

acembly: $(ACEMBLY)

OTHER= blixem dotter efetch

other: $(OTHER)

TOOLS= acediff acediffsorted pmapace diskdump diskfix homonym tagcount 
tools : $(TOOLS)

MISC= acctest acstest \
		aceview gffdump   \
		acesplit pseudoword_test \
		av.cgi pfam.cgi \
		blastd blast.cgi blastd_convert_table allblast blastct \
		wego_test


misc: $(MISC)

# RPCACE is no longer used, it is replaced by taceclient taceserver
# although taceserver does not suypport some fancy authentication developped at the sanger around 2003
RPCACE = aceserver aceclient jade2sybase xaceclient gifaceserver gifacemblyserver gffread 


rpcace: $(RPCACE)

#
# BROKEN is a list of things that do not compile, but that I do not
# know that we really need

BROKEN=xtest coltest giftest boxtest jadeserver tqp testvahan xseq  testprog 

#
# ACE5 support is incomplete

ACE5= tace5 xace5 xacembly5 tacembly5

#### ACEDB programs

$(LIBS) $(ACE) $(ACEMBLY) $(OTHER) $(TOOLS) $(MISC) $(ACE5) $(BROKEN) dist : bin-dir
	echo '...Making' $@ 'for $(ACEDB_MACHINE)'
	cd bin.$(ACEDB_MACHINE) ; \
	cp ../wmake/$(ACEDB_MACHINE)_DEF ./deffile; \
	cp ../wmake/truemake ./makefile ; \
	chmod u+w ./makefile ./deffile; touch ./makefile ./deffile; \
	$(MAKE) $@

# Really this should be merged with the above, all that's different is the message...
libfree_clean libgraph_clean libace_clean libxace_clean \
		xace_clean tace_clean \
		aceclient_clean aceserver_clean \
		sace_clean \
		xacembly_clean tacembly_clean \
		gifacemblyserver_clean : bin-dir
	@echo "removing `echo $@ | cut -d_ -f1` .o files for $(ACEDB_MACHINE)"
	@cd bin.$(ACEDB_MACHINE) ; \
	cp ../wmake/truemake ./makefile ; \
	cp ../wmake/$(ACEDB_MACHINE)_DEF ./deffile; \
	chmod u+w ./makefile ./deffile; touch ./makefile ./deffile; \
	$(MAKE) $@

#### CONTIG programs

contig_all:
	@echo ""
	@echo " Compiling ContigC programs"
	@echo "-------------------------------------------------"
	$(MAKE) contig pmap fphinx fpabi openmapdb mapsub extend \
		getascdb lastclones contasp logdb rancl randb \
		model anytoany checkdb rmvec
	@echo ""
	@echo "These binaries are available :"
	@echo ""
	@ls -l bin.$(ACEDB_MACHINE)  | grep '^-..x'
	@echo ""

acext:
	echo "Making all AceC applications in wacext"
	cd wacext; make  all

acenn:
	echo "Making all NeuralNet applications in wnn"
	cd wnn ; make  all

magic:
	echo "Making magic code"
	cd wsa ; make  all

#####################################################################
#### binary directory bin.$(ACEDB_MACHINE)
####
#### Dependency used by all of the above targets to check that user
#### has set up the ACEDB_MACHINE environment variable correctly.
####
bin-dir:
	@if test $(ACEDB_MACHINE) ; then exit ; \
	else \
	echo  'To compile, please setenv ACEDB_MACHINE (now not set) to one of: ' ; \
	cd wmake ; echo ; ls *DEF ; \
	echo	 '      without the _DEF ending,' ; \
	echo     '      e.g.  add in your .login file the line '; \
	echo  ' setenv ACEDB_MACHINE  SUN' ; \
	echo '  or some other machine name defined in the directory wmake.' ;\
	echo '     Thank you' ; echo ; echo ; exit 1 ; fi
	@if test -f wmake/$(ACEDB_MACHINE)_DEF ; then exit ; \
	else echo  'To compile, please setenv ACEDB_MACHINE (now = $(ACEDB_MACHINE) ) to one of: ' ; \
	cd wmake ; echo ; ls *DEF ; \
	echo	 '      without the _DEF ending,' ; \
	echo     '      i.e.  add in your .login file the line '; \
	echo  ' setenv ACEDB_MACHINE  SUN' ; \
	echo '  or some other machine name defined in the directory wmake.' ;\
	echo '     Thank you' ; echo ; echo ; exit 1 ; fi
	@if test -d bin.$(ACEDB_MACHINE)  ; then exit ; \
	else mkdir bin.$(ACEDB_MACHINE)  ; fi

#####################################################################
#### end of makefile
#####################################################################

