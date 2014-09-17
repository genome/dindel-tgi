PREFIX=/usr
BINDIR=$(PREFIX)/bin
SHAREDIR=$(PREFIX)/share
SAMTOOLSDIR=samtools
SEQANDIR=seqan_library
DINDEL=dindel$(SUFFIX)

CPPFLAGS=-DNDEBUG -D_IOLIB=2 -DMINREADS=2 -DDINDEL
CXXFLAGS=-I$(SAMTOOLSDIR) -I$(SEQANDIR) -I./ -Wno-deprecated  -O3 
LDFLAGS=$(SAMTOOLSDIR)/libbam.a -lz -lboost_program_options -lpthread

SRCSDINDEL=DInDel.cpp HapBlock.cpp HaplotypeDistribution.cpp ObservationModelFB.cpp GetCandidates.cpp Faster.cpp
OBJSDINDEL=$(SRCSDINDEL:%.cpp=%.o)

all: $(DINDEL)

$(DINDEL):$(OBJSDINDEL) Read.hpp DInDel.hpp HapBlock.hpp Haplotype.hpp HaplotypeDistribution.hpp MyBam.hpp GetCandidates.hpp Variant.hpp Fasta.hpp OutputData.hpp MLAlignment.hpp ObservationModelSeqAn.hpp VariantFile.hpp ReadIndelErrorModel.hpp Library.hpp Faster.hpp
	$(CXX) -o $@ $(CXXFLAGS) $(DINDELFLAGS) $(OBJSDINDEL) $(LDFLAGS) 

clean:
	rm -f $(OBJSDINDEL) $(DINDEL)

install: all
	install $(DINDEL) $(DESTDIR)$(BINDIR)
	cp -dr --no-preserve=ownership python/* $(DESTDIR)$(SHAREDIR)/$(DINDEL)/
