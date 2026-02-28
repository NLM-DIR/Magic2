# Makefile

ifeq ($(shell uname),Darwin)
	STLINC = /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1

	ifeq ($(shell uname -m),aarch64)
		ARCH = arm
	else
		ARCH = x86
	endif
	OS = mac-$(ARCH)
else
	GCC_VERSION = "$(shell g++ --version | awk '/^g++/ {split($NF, a, "."); print a[1];}')"
	STLINC = /usr/include/c++/$(GCC_VERSION)
	OS = linux
endif


# NGS-SDK version
VERSION = 3.2.1

NGSSDK = ngs-sdk.$(VERSION)-$(OS)
URL = https://ftp-trace.ncbi.nlm.nih.gov/sra/ngs/$(VERSION)/$(NGSSDK).tar.gz

LIBS = -L$(NGSSDK)/lib64 -lncbi-ngs-c++-static -lncbi-ngs-static -lncbi-vdb-static -lngs-c++-static -lpthread -ldl



all: $(NGSSDK)/README
	$(MAKE) clean
	$(MAKE) libsraread.a


$(NGSSDK)/README:
	curl -O $(URL)
	tar xfz $(NGSSDK).tar.gz


%.o: %.c
	cc -c -o $@ $< -I ../wh

%.o: %.cpp
	cc -c -I$(NGSSDK)/include -I$(NGSSDK)/include/ncbi-vdb -I$(STLINC) -o $@ $<


libsraread.a: sra_read.o
	ar rcs $@ $^


demo: demo.o libsraread.a
	g++ -o $@ $^ $(LIBS)

clean:
	rm -f *.o
	rm -f *.a
	rm -f demo
