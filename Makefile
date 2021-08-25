#Compiler
CC=g++
C=gcc

#compiler flags
# -g adds debug info to the executable file
# -Wall turns on most warnings from compiler
CCFLAGS = -g -std=c++0x -DNANOSEQ_VERSION='"$(NANOSEQ_VERSION)"'  -Wall
CFLAGS  = -g -Wall -O2 -DNANOSEQ_VERSION='"$(NANOSEQ_VERSION)"' 

PREFIX?=/usr/local/

#Define locations of header files
ifdef HTSLIB
HTSINC=$(HTSLIB)
else
HTSLIB=$(PREFIX)/lib
HTSINC=$(PREFIX)/include
endif

INCLUDES= -I$(HTSINC)
LFLAGS= -L$(HTSLIB)
INSTALL=$(PREFIX)/bin

# define any libraries to link into executable:
LIBS =-l:libhts.a -l:libdeflate.a -lgzstream -lz -lpthread -lcurl -ldl -llzma -lbz2 -lm

# define the C source files
SRCS=./src/bamaddreadbundles.cc ./src/statistics.cc ./src/dsa.cc ./src/read_bundler.cc ./src/randomreadinbundle.cc ./gridsearch/gridsearch.cc ./src/bed_reader.cc ./src/pileup.cc ./src/writeout.cc ./src/variantcaller.cc ./src/bamcov.c ./src/sam_opts.c ./src/sam_utils.c

MD := mkdir -p

#Build target executable
BAMFILTER_TARGET=./src/bamaddreadbundles
GRIDSEARCH_TARGET=./src/gridsearch
DSA_TARGET=./src/dsa
STATISTICS_TARGET=./src/statistics
VARIANTCALLER_TARGET=./src/variantcaller
RANDOMREADINBUNDLE_TARGET=./src/randomreadinbundle
BAMCOV_TARGET=./src/bamcov

TARGETS=$(BAMFILTER_TARGET) $(GRIDSEARCH_TARGET) $(DSA_TARGET) $(STATISTICS_TARGET) $(VARIANTCALLER_TARGET) $(RANDOMREADINBUNDLE_TARGET) $(BAMCOV_TARGET)

.PHONY: all install clean

all: $(BAMFILTER_TARGET) $(GRIDSEARCH_TARGET) $(DSA_TARGET) $(STATISTICS_TARGET) $(VARIANTCALLER_TARGET) $(RANDOMREADINBUNDLE_TARGET)
	@echo  Binaries have been compiled.

$(HTSLIB)/libhts.a:
	@echo Need HTSLIB to point to htslib installation
	@exit 1

$(BAMFILTER_TARGET):$(HTSLIB)/libhts.a ./src/bamaddreadbundles.cc ./src/bamaddreadbundles.h
	$(CC) $(INCLUDES) $(CCFLAGS) -o $(BAMFILTER_TARGET) $(LFLAGS) $(CAT_LFLAGS) ./src/bamaddreadbundles.cc $(LIBS)

$(GRIDSEARCH_TARGET):$(HTSLIB)/libhts.a ./src/gridsearch.h
	$(CC) $(INCLUDES) $(CCFLAGS) -o $(GRIDSEARCH_TARGET) $(LFLAGS) ./src/gridsearch.cc $(LIBS)

$(DSA_TARGET):$(HTSLIB)/libhts.a ./src/dsa.cc ./src/pileup.cc ./src/bed_reader.cc ./src/read_bundler.cc ./src/writeout.cc ./src/pileup.h ./src/bed_reader.h ./src/read_bundler.h ./src/writeout.h
	$(CC) $(INCLUDES) $(CCFLAGS) -o $(DSA_TARGET) $(LFLAGS) ./src/dsa.cc ./src/pileup.cc ./src/bed_reader.cc ./src/read_bundler.cc ./src/writeout.cc $(LIBS)

$(STATISTICS_TARGET):$(HTSLIB)/libhts.a ./src/statistics.cc ./src/statistics.h
	$(CC) $(INCLUDES) $(CCFLAGS) -o $(STATISTICS_TARGET) $(LFLAGS) ./src/statistics.cc $(LIBS)

$(VARIANTCALLER_TARGET):$(HTSLIB)/libhts.a ./src/variantcaller.cc ./src/variantcaller.h
	$(CC) $(INCLUDES) $(CCFLAGS) -o $(VARIANTCALLER_TARGET) $(LFLAGS) ./src/variantcaller.cc $(LIBS)

$(RANDOMREADINBUNDLE_TARGET):$(HTSLIB)/libhts.a ./src/randomreadinbundle.cc ./src/randomreadinbundle.h
	$(CC) $(INCLUDES) $(CCFLAGS) -o $(RANDOMREADINBUNDLE_TARGET) $(LFLAGS) ./src/randomreadinbundle.cc $(LIBS)

$(BAMCOV_TARGET):$(HTSLIB)/libhts.a ./src/bamcov.c ./src/sam_opts.c ./src/sam_utils.c
	$(C) $(INCLUDES) $(CFLAGS) -o $(BAMCOV_TARGET) $(LFLAGS) ./src/bamcov.c ./src/sam_opts.c ./src/sam_utils.c $(LIBS)


install:$(TARGETS)
	$(MD) $(INSTALL)
	install $(TARGETS) $(INSTALL)

test :$(TARGETS)
	@echo
	@echo Running tests
	@cd ./test/;./run1;
	@cd ./test/;./run2;
	@cd ./test/;./run3;
	@cd ./test/;./run4;


clean:
	$(RM) $(TARGETS)
	$(RM) ./test/table.test.bed.gz ./test/variant.test ./test/test.neat.bam ./test/test.post.bam
# DO NOT DELETE THIS LINE -- make depend needs it
