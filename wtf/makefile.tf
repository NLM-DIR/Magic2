#  %W% %G%
#################################################################
############### acedb: R.Durbin and J.Thierry-Mieg ##############
############### wext/makefile.ext  Feb 1998        ##############
#################################################################
#
# ATTENTION:
#   Define in this file machine independant directives
#   to compile external codes linked against acelib
#
#   The makefile is decomposed into 2 parts
#     makefile: which selects the machine and creates ./bin.*
#     then copies the present file into ../bin.*
#     and the present file which gives details on the linking
#  
#   Machine specific choices, like linker and compiler
#     are specified in ../wmake/*_DEF
#
#   Object modules and exe are created in ../bin.*
#     to keep the present directory clean and allow
#     multiple architectures
#
#################################################################
# Before running the make command
# setenv ACEDB_MACHINE to indicate one of the predefined
# Machine dependant option files, or write your own.
# These files are called  $(ACEDB_SRC)/wmake/$(ACEDB_MACHINE)_DEF
#
# See the explanation at the top of ../wmake/truemake
#################################################################

# LINK_ACC is the set of object files needed if you use
# the client version of the new Ace C library
LINK_ACC= libacedna.a libtcpfn.a $(LINK_ACC_RPC) libaccl.a  libfree.a  libtsfree.a
LINK_LIGHT=libacedna.a libwego.a libfree.a  libtsfree.a -lpthread 

# LINK_ACS is the set of object files needed if you use
# the standalone acinside version of the new Ace C library
LINK_ACS=libace.a libfree.a libtsfree.a libwego.a

.KEEP_STATE: 

# suppress auto SCCS extraction
.SCCS_GET:

#define defaults  overridable in $(ACEDB_MACHINE)_DEF
RANLIB_NEEDED = true	
AR_OPTIONS = rlu

LINKER = cc 
COMPILER = cc
TEST = test

#################################################################
# most versions of lex require you to use -ll
# but not Linux - reset LEX_LIBS in $(ACEDB_MACHINE)_DEF file

LEX_LIBS = -ll

#################################################################
# include your choice from a machine dependant file
# do not edit the present makefile for this purpos
# this simplifies multiple architecture maintainance

include deffile

# Note that you can keep different DEF files for the same machine
# setting  various compiler  options

###########################################################
##  Compiler and library options
## CC, LIBS, NAME are defined in $(ACEDB_MACHINE)_DEF
##
 
IDIR = -I. -I.. -I../wh -I/netopt/sge/include -I ../wtf/TENSORFLOW/clib/include

# Do not use -I/usr/include
# it prevents gcc from picking up its own includes
# (cc goes to /usr/include anyway)

## to undefine any rubbish
CCFLAGS =
GCFLAGS =

## Different platforms use CC or COMPILE.c
#  (USEROPTS - see comments at top of file)
#
CC =        $(COMPILER) $(USEROPTS) $(IDIR) -DAC_TEST -D$(NAME) -c
COMPILE.c = $(COMPILER) $(USEROPTS) $(IDIR) -DAC_TEST -D$(NAME) -c


#################################################################
#################################################################
# Add here your definitions
#################################################################
#################################################################acc
 
ALL_TF_OBJS = tf1.o tf2c.o tfnnali8b.o
ALL_TF_SOURCES =  $(ALL_TF_OBJS:.o=.c) 
LINK_TF = tf2c.o ../wtf/TENSORFLOW/clib/lib/libtensorflow.so ../wtf/TENSORFLOW/clib/lib/libtensorflow_framework.so -L ../wtf/TENSORFLOW/clib/lib

$(ALL_TF_SOURCES)  :
	$(TEST) -L $@ || ln -s ../wtf/$@ . 

# "all" should always be the first target so that it is the default make action.
all: tf1 tfnnali8b

tf1: tf1.o $(LINK_ACS) tf2c.o
	$(LINKER)  -o $@  $@.o  $(LINK_ACS) $(LIBS) $(LINK_TF)

tfnnali8b: tfnnali8b.o $(LINK_ACS) tf2c.o
	$(LINKER)  -o $@  $@.o  $(LINK_ACS) $(LIBS) $(LINK_TF)

###########################################################
########### end of     acedb makefile.nn   ################
###########################################################
# Requires the TensorFlow C library library and headers, so download those.

# if (! -d "clib") then
#  echo "Downloading TensorFlow C library into TENSORFLOW/clib"
#  mkdir TENSORFLOW
#  mkdir TENSORFLOW/clib
#  pushd  TENSORFLOW/clib
#    wget "https://storage.googleapis.com/tensorflow/libtensorflow/libtensorflow-cpu-linux-x86_64-1.4.0.tar.gz" 
#    # | tar -C clib -xz
#    tar xf libtensorflow-cpu-linux-x86_64-1.4.0.tar.gz
#  popd
# endif

# gcc -std=c99 -I clib/include -L clib/lib train.c -ltensorflow -ltensorflow_framework
# LD_LIBRARY_PATH=clib/lib ./a.out