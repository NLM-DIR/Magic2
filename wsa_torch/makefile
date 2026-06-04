# wsa_torch/makefile
# Mirrors the wsa/makefile pattern exactly.
# Copies makefile.sa_torch into bin.$(ACEDB_MACHINE) and builds there.
# The wsa_torch/ source directory is never written to.

SHELL = /bin/sh

.SCCS_GET:

PGM = sortalign_torch

all $(PGM) : bin-dir copy-makefile.sa_torch
	echo 'making' $@ 'for $(ACEDB_MACHINE)' ; echo
	cd ../bin.$(ACEDB_MACHINE) ; \
	$(MAKE) -f makefile.sa_torch $@

copy-makefile.sa_torch:
	cp makefile.sa_torch ../bin.$(ACEDB_MACHINE) ; \
	cp ../wmake/$(ACEDB_MACHINE)_DEF ../bin.$(ACEDB_MACHINE)/deffile ; \
	chmod u+w ../bin.$(ACEDB_MACHINE)/makefile.sa_torch

bin-dir:
	if test $(ACEDB_MACHINE) ; then exit ; \
	else \
	echo 'To compile, please setenv ACEDB_MACHINE (now not set) to one of: ' ; \
	cd ../wmake ; echo ; ls *DEF | sed 's/_DEF//' ; \
	exit 1 ; fi
	if test -f ../wmake/$(ACEDB_MACHINE)_DEF ; then exit ; \
	else echo 'Unknown ACEDB_MACHINE=$(ACEDB_MACHINE)' ; exit 1 ; fi
	if test -d ../bin.$(ACEDB_MACHINE) ; then exit ; \
	else mkdir ../bin.$(ACEDB_MACHINE) ; fi

# end of wsa_torch/makefile
