# @(#)makefile	1.58 10/1/97

SHELL = /bin/sh
#  MAKE = /usr/local/bin/gmake
#  MAKE = /usr/ccs/bin/make on solaris 

#.SILENT:

# Suppress auto SCCS extraction.
.SCCS_GET:

####  Program list 

PGM = sortalign

#### EXTERNAL acedb programs linked against acelib

all $(PGM) : bin-dir copy-makefile.sa
	echo 'making' $@ 'for $(ACEDB_MACHINE)' ; echo
	cd ../bin.$(ACEDB_MACHINE) ; \
	$(MAKE) -f makefile.sa $@

copy-makefile.sa:
	cp makefile.sa ../bin.$(ACEDB_MACHINE) ; \
	cp ../wmake/$(ACEDB_MACHINE)_DEF ../bin.$(ACEDB_MACHINE)/deffile ; \
	chmod u+w ../bin.$(ACEDB_MACHINE)/makefile.sa

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

 
 
 
 
