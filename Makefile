SAMTOOLDIR=samtools/
SEQANDIR=seqan_library/

CPPFLAGS= -DNDEBUG -D_IOLIB=2 -DMINREADS=2 -DDINDEL
CXXFLAGS= -I$(SAMTOOLDIR) -I$(SEQANDIR) -I./  -Wno-deprecated  -O3 
LDFLAGS= -L$(SAMTOOLDIR)  -lbam -lz -lboost_program_options -static 

SRCSDINDEL=DInDel.cpp HapBlock.cpp HaplotypeDistribution.cpp ObservationModelFB.cpp GetCandidates.cpp Faster.cpp
OBJSDINDEL=$(SRCSDINDEL:%.cpp=%.o)  

dindel:$(OBJSDINDEL) Read.hpp DInDel.hpp HapBlock.hpp Haplotype.hpp HaplotypeDistribution.hpp MyBam.hpp GetCandidates.hpp Variant.hpp Fasta.hpp OutputData.hpp MLAlignment.hpp ObservationModelSeqAn.hpp VariantFile.hpp ReadIndelErrorModel.hpp Library.hpp Faster.hpp
	$(CXX) -o $@ $(CXXFLAGS) $(DINDELFLAGS) $(OBJSDINDEL) $(LDFLAGS) 

clean:
	rm -f $(OBJSDINDEL) $(OBJSCOMPAREVARIANTS)  $(OBJSMAKEGLF)
