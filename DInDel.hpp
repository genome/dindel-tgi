/*    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef DINDEL_HPP_
#define DINDEL_HPP_
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <boost/tuple/tuple.hpp>
#include <ext/hash_map>

#include "MyBam.hpp"
#include "faidx.h"
#include "Haplotype.hpp"
#include "ObservationModel.hpp"
#include "HaplotypeDistribution.hpp"
#include "ObservationModelFB.hpp"
//#include "Fast.hpp"
#include "MLAlignment.hpp"
#include "Read.hpp"
#include "StringHash.hpp"

#include "OutputData.hpp"
#include "Library.hpp"
#include "VariantFile.hpp"

const int SHIFTSTRAND = 1000000; // used to keep track of forward and reverse matches in ::filterHaplotypes

using namespace std;
using namespace boost;
using __gnu_cxx::hash;
typedef struct
{
	double pOff, pOn;
} HapReadLik;





class VariantCoverage {
public:
	VariantCoverage()
	{
		nf=0;
		nr=0;
	}
	VariantCoverage(int _nf, int _nr) 
	{
		nf = _nf;
		nr = _nr;
	}
	int nf, nr; // forward and reverse
};

class DetInDel
{
public:
	//DetInDel(const string & bfName, const string & tid, const string &outputFileName, const string & modelType) : params(tid, outputFileName, modelType) { fai=NULL; };
	static int fetchFuncFindInDel(const bam1_t *b, void *data);
	void findInDels(uint32_t start, uint32_t end, bool report);
	void detectIndels(const string & variantsFileName);
	void callVariants(const string & variantsFile);
	void findInDelsPositionsFile(const string & fileName);
	string getRefSeq(uint32_t lpos, uint32_t rpos);
	void empiricalDistributionMethod(int index, const vector<Read> & reads, uint32_t pos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants, OutputData & oData, OutputData & glfData);
	void fastMethod(const vector<Read> & reads, uint32_t pos, uint32_t leftPos, uint32_t rightPos, ostream & output);
	bool getHaplotypes(vector<Haplotype> & haps, const vector<Read> & reads, uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, const AlignedCandidates & candidateVariants);
	const vector<MyBam *> & getMyBams() const { return myBams; }
	class HapPairLik {
	public:
		double ll;
		int h1, h2;
		int numIndFirst;
		int numIndSecond;
		int numFirst, numSecond; // number of reads mapped to first and second
		int numOffBoth; // number of reads that do not map to either haplotype
		double numOffBothError;
		map<int, VariantCoverage> hapIndelCoverage1, hapSNPCoverage1,hapIndelCoverage2, hapSNPCoverage2; // indels and snps in the _haplotype_ covered by the read

		operator double() const { return ll;};
	};

	class HapEstResult {
		public:
			HapEstResult();
			HapEstResult(const AlignedVariant & _av, int _pos, double _prob, double _freq, int _nrf, int _nrr) {
				av=_av;
				pos=_pos;
				prob=_prob;
				freq=_freq;
				nrf=_nrf;
				nrr=_nrr;
			};
			AlignedVariant av;
			int pos;
			double prob;
			double freq;
			int nrf; // number of reads on reverse strand
			int nrr; // number of reads on forward strand
	};

	void addLibrary ( const string & name, const Library & lib)
	{
		libraries[name.c_str()]=lib;
	}
	void addLibrary ( const string & fileName)
	{
		libraries.addFromFile(fileName);
	}


protected:
	void outputHapsAndFreqs(ostream *output, const string & prefix, const vector<Haplotype> & haps, const vector<double> & freqs, uint32_t leftPos);
	//void getReads(uint32_t leftPos, uint32_t rightPos, vector<Read> & reads);
	void getReads(uint32_t leftPos, uint32_t rightPos, vector<Read> & reads, uint32_t & oldLeftPos, uint32_t  & oldRightFetchReadPos, vector<Read *> & readBuffer, bool reset);

	double getMaxHap(Haplotype & h1, Haplotype &h2, HapPairLik & hpl, const vector<Haplotype> & haps, vector<HapPairLik> & likPairs);
	void outputMaxHap(ostream *output, const string & prefix, const vector<Haplotype> & haps, vector<HapPairLik> & likPairs);
	void outputTopHaps(ostream *output, const string & prefix, const vector<Haplotype> & haps, vector<HapPairLik> & likPairs, int n);
	bool alignHaplotypes(vector<Haplotype> & haps, uint32_t pos, uint32_t & leftPos, uint32_t & rightPos,  map<int, set<AlignedVariant> > & variants);
	bool generateHaplotypes(vector<Haplotype> & haps, uint32_t pos, uint32_t & leftPos, uint32_t & rightPos,  const map<int, set<Variant> > & variants);
	double getHaplotypePrior(const Haplotype & h1, const Haplotype & h2, int leftPos, const AlignedCandidates & candidateVariants);
	void computeLikelihoods(const vector<Haplotype> &haps, const vector<Read> & reads, vector<vector<MLAlignment> > & liks, uint32_t leftPos, uint32_t rightPos, vector<int> & onHap);

	void computeHapPosition(const Haplotype & hap, const Read & read, vector<int> & alPos, int leftPos);
	void computeLikelihoodsFaster(const vector<Haplotype> &haps, const vector<Read> & reads, vector<vector<MLAlignment> > & liks, uint32_t leftPos, uint32_t rightPos, vector<int> & onHap);

	void computePairLikelihoods(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<HapPairLik> & likPairs, bool usePrior, const AlignedCandidates & candidateVariants, int leftPos);
	void statisticsHaplotypePair(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, HapPairLik & hpl,OutputData::Line & line);

	void estimateHaplotypeFrequencies(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<double> & hapFreqs);
	void estimateHaplotypeFrequenciesPosterior(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<double> & hapFreqs, map <int, vector<tuple<AlignedVariant, double,double> > > & posteriors, uint32_t pos, uint32_t leftPos, ostream & glfOutput);
	void estimateHaplotypeFrequenciesBayesEM(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<double> & hapFreqs, vector <HapEstResult > & posteriors,  uint32_t candPos, uint32_t leftPos,   uint32_t rightPos, OutputData & glfData, int index, const AlignedCandidates & candidateVariants,string program);
	void diploidGLF(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<double> & hapFreqs, vector <HapEstResult > & posteriors,  uint32_t candPos, uint32_t leftPos,  uint32_t rightPos, OutputData & glfData, int index, const AlignedCandidates & candidateVariants, string program);


	void debug(const pair<Haplotype, Haplotype> & hp, const vector<Read> & reads,  uint32_t leftPos, uint32_t rightPos);
	void debug(const pair<Haplotype, Haplotype> & hp1, const pair<Haplotype, Haplotype> & hp2, const vector<Read> & reads,  uint32_t leftPos, uint32_t rightPos);
	void analyzeDifference(const pair<Haplotype, Haplotype> & hp1, const vector<Read> & reads,  uint32_t leftPos, uint32_t rightPos);
	void showAlignments(const pair<Haplotype, Haplotype> & hp1, const vector<Read> & reads,  uint32_t leftPos, uint32_t rightPos);
	void showAlignmentsPerHaplotype(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, uint32_t candPos, uint32_t leftPos);

	double getPairPrior(const AlignedVariant & av1, const AlignedVariant & av2, int leftPos,const AlignedCandidates & candidateVariants);

	void filterHaplotypes(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks,  vector<int> & filtered,  map<pair<int, AlignedVariant>, VariantCoverage> & varCoverage, bool doFilter);

	//MyBam myBam;

	//MyBam myBam;
	vector<MyBam *> myBams;
	vector<string> myBamsFileNames;
	LibraryCollection libraries;


	class CIGAR : public vector<pair<int,int> >
	{
	public:
		typedef pair<int,int> CIGOp;
		int refPos;
	};
	CIGAR getCIGAR(const Haplotype & hap, const Read & read, const MLAlignment & ml, int refSeqStart);
	void writeRealignedBAMFile(const string & fileName, const vector<CIGAR> & cigars, const vector<Read> & reads, const vector<int> & onHap, const bam_header_t *bh);
	void writeUnalignedBAMFile(const string & fileName, const vector<Read> & reads, const vector<int> & onHap, const bam_header_t *bh);
	class InDel {
	public:
		InDel()
		{
			count[0]=0;
			count[1]=0;
		}
		typedef enum { In, Del} Type;
		Type type;
		size_t count[2];
	};

public:
	class Parameters {
	public:
		Parameters(const string & _tid, string _fileName, const string & modelType) : obsParams(modelType)
		{
			tid=_tid;
			fileName=_fileName;
			setDefaultValues();
		}
		void setDefaultValues()
		{
			bayesa0=0.001;
			width=30;
			maxHap=100;
			skipMaxHap=1000;
			maxReads=500;
			mapQualThreshold=0.995;
			glfNumHap=5;
			inferenceMethod="empirical";
			minReadOverlap=5;
			minCount=2;
			maxReadLength=40;
			numOutputTopHap=5;
			checkAllCIGARs=1;
			bayesType="all";

			fastWidth=4;
			analyzeLowFreq=false;
			analyzeLowFreqDiffThreshold=1.0;
			showHapDist=true;
			showHapAlignments=false;
			showCandHap=false;
			showReads=false;
			fastWidthOverlap=4;
			noIndelWindow=-1;
			mapUnmappedReads=false;
			priorIndel=1.0/10000;
			priorSNP=1.0/1000.0;
			filterReadAux=string("");
			quiet=true;
			computeML=false;
			computeMAP=false;
			doDiploid=false;
			slower=true;
			estimateHapFreqs=false;
			printCallsOnly=true;
			outputPooledLikelihoods=false;
			filterHaplotypes = false;

			outputGLF=true;
			outputRealignedBAM=false;
			processRealignedBAM="no";
			changeINStoN = false;


			EMtol=1e-4;
		}
		OutputData makeOutputData(ostream & out)
		{
			OutputData oData(out);
			oData("msg")("index");
			oData("analysis_type");
			oData("tid")("lpos")("rpos")("center_position")("realigned_position");
			oData("ref_all")("num_reads")("num_hqreads");
			oData("post_prob_variant")("est_freq")("was_candidate_in_window");

			oData("num_mapped_to_first")("num_mapped_to_second");
			oData("num_off_hap")("loglik_hap_pair")("loglik_next_hap_pair");
			oData("first_var_cover_forward")("first_var_cover_reverse")("second_var_cover_forward")("second_var_cover_reverse");
			oData("first_called_all")("second_called_all")("loglik_called_genotype")("loglik_ref_ref")("alt_genotypes");
			return oData;
		}

		OutputData makeGLFOutputData(ostream & out)
		{
			OutputData oData(out);
			oData("msg")("index");
			oData("analysis_type");
			oData("tid")("lpos")("rpos")("center_position")("realigned_position")("was_candidate_in_window");
			oData("ref_all")("nref_all")("num_reads");
			oData("post_prob_variant")("qual")("est_freq")("logZ")("hapfreqs");

			oData("indidx")("msq")("numOffAll")("num_indel")("num_cover_forward")("num_cover_reverse")("num_unmapped_realigned");
			oData("var_coverage_forward")("var_coverage_reverse");
			oData("nBQT")("nmmBQT")("mLogBQ")("nMMLeft")("nMMRight");
			oData("glf");
			return oData;
		}

		OutputData makeGLFv2OutputData(ostream & out)
		{
			OutputData oData(out);
			oData("msg")("index");
			oData("analysis_type");
			oData("tid")("candidate_position")("realigned_position");
			oData("ref_all")("nref_all")("num_reads");
			oData("post_prob_variant")("est_freq");

			oData("indidx")("msq")("num_cover_forward")("num_cover_reverse");
			oData("glf");
			return oData;
		}


		void print()
		{
			cout << "DetInDel parameters: " << endl;
			cout << "\ttid: " << tid << " width: " << width << " maxHap: " << maxHap << " maxReads: " << maxReads << " skipMaxHap: " << skipMaxHap << endl;
			cout << "\toutputFilename: " << fileName << endl;
			cout << "\tmapQualThreshold: " << mapQualThreshold << endl;
			//cout << "\tscaleError: " << scaleErr << endl;
			cout << "\tinferenceMethod: " << inferenceMethod << endl;
			//cout << "\tglfNumHap: " << glfNumHap << endl;
			cout << "\tanalyzeLowFreq: " << analyzeLowFreq << endl;
			cout << "\tanalyzeLowFreqDiffThreshold: " << analyzeLowFreqDiffThreshold << endl;
			cout << "\tshowHapDist: " << showHapDist << endl;
			cout << "\tminReadOverlap: " << minReadOverlap << endl;
			cout << "\tmaxReadLength: " << maxReadLength << endl;
			//cout << "\tminCount: " << minCount << endl;
			cout << "\tmaxHapReadProd: " << maxHapReadProd << endl;
			//cout << "\tfastWidth: " << fastWidth << endl;
			//cout << "\tfastWidthOverlap: " << fastWidthOverlap << endl;
			cout << "\tshowCandHap: " << showCandHap << endl;
			cout << "\tshowReads: " << showReads << endl;
			cout << "\tfilterHaplotypes: " << filterHaplotypes << endl;
			cout << "\tnoIndelWindow: " << noIndelWindow << endl;
			cout << "\tmapUnmappedReads: " << mapUnmappedReads << endl;

			cout << "\tnumOutputTopHap: " << numOutputTopHap << endl;

			cout << "\tcheckAllCIGARs: " << checkAllCIGARs << endl;
			cout << "\tchangeINStoN: " << changeINStoN << endl;
			


			
			cout << endl;
			cout << "\tquiet: " << quiet << endl;
			cout << "\tprintCallsOnly: " << printCallsOnly << endl;
			cout << "\tfaster: " << !slower << endl;
			cout << "\tdoDiploid: " << doDiploid << endl;
			cout << "\tdoEM: " << estimateHapFreqs << endl;

			cout << "\toutputPooledLikelihoods: " << outputPooledLikelihoods << endl;
			cout << "\toutputRealignedBAM: " << outputRealignedBAM << endl;
			cout << "\tprocessRealignedBAM: " << processRealignedBAM << endl;
			cout << "\tshowHapAlignments: " << showHapAlignments << endl;

			cout << "\tEM tol: " << EMtol << endl;
			cout << "\tbayesEM a0: " << bayesa0 << endl;
			cout << "\tbayesType: " << bayesType << endl;


			cout << "\tpriorIndel: " << priorIndel << endl;
			cout << "\tpriorSNP: " << priorSNP << endl;

			//cout << "\tmeanInsert: " << meanInsert << endl;
			//cout << "\tstdInsert: " << stdInsert << endl;

			cout << "\tfilterReadAux: " << filterReadAux << endl;

			cout << "Observation model parameters: " << endl;
			obsParams.print();
		}
		int noIndelWindow, numOutputTopHap, checkAllCIGARs, minReadOverlap, maxHapReadProd;
		uint32_t width, maxHap, maxReads, skipMaxHap, glfNumHap,  maxReadLength, minCount, fastWidth, fastWidthOverlap;
		double checkBaseQualThreshold;
		double mapQualThreshold, scaleErr, priorIndel, priorSNP, EMtol, bayesa0;
		string fileName, inferenceMethod, refFileName, tid, filterReadAux, bayesType, processRealignedBAM;
		bool analyzeLowFreq, showHapDist, showCandHap, showReads, showHapAlignments, alignAgainstReference, mapUnmappedReads, quiet, estimateHapFreqs, doDiploid, computeML, computeMAP, slower,printCallsOnly, outputPooledLikelihoods, filterHaplotypes;
		bool outputRealignedBAM, outputGLF, varFileIsOneBased, changeINStoN;
		double analyzeLowFreqDiffThreshold;
		double meanInsert, stdInsert;
		ObservationModelParameters obsParams;
	};
	Parameters params;

	DetInDel(const string & bfName, const Parameters & _params, int multipleFiles);
	~DetInDel();


	map<uint32_t, InDel> indels;
	class ScanStats
	{
		public:
		ScanStats()
		{
			numUnmappedMate=0;
		}
		int numUnmappedMate;
	};
	ScanStats scanStats;
protected:
	faidx_t *fai;

};



class FFData
{
public:
	uint32_t start, end;
	DetInDel *det;
	map<string, int> unmappedMate;
	map<int, int > insHisto, delHisto;
};

#endif /*DINDEL_HPP_*/
