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
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <boost/foreach.hpp>
#include "bam.h"
#include "DInDel.hpp"
#include "Haplotype.hpp"
#include "HaplotypeDistribution.hpp"
#include "ObservationModelFB.hpp"
#include "Utils.hpp"
#include "faidx.h"
#include "GetCandidates.hpp"
#include "ObservationModelSeqAn.hpp"
#include "VariantFile.hpp"
#include "Faster.hpp"
#include <ext/hash_map>
#include <exception>

const int USECALLWINDOW=0;
//using namespace seqan;
using namespace seqan;
namespace po = boost::program_options;

using namespace std;
//using namespace fasta;


DetInDel::DetInDel(const string & bfName, const Parameters & _params, int multipleFiles) : params(_params)

{
	fai=NULL;
	if (params.alignAgainstReference) {
		fai = fai_load(params.refFileName.c_str());
		if (!fai) {
			cerr << "Cannot open reference sequence file." << endl;
			params.alignAgainstReference=false;
			exit(1);
		}
	}

	if (multipleFiles==0) {
		myBams.push_back(new MyBam(bfName));
	} else {

		ifstream file(bfName.c_str());
		if (!file.is_open()) {
			   cout << "Cannot open file with BAM files:  " << bfName << endl;
			   throw string("File open error.");
		}
		while (!file.eof()) {
			string line;
			getline(file, line);
			if (!line.empty()) {
				istringstream is(line);
				string fname;
				is >> fname;
				if (!fname.empty()) {
					cout << "Reading BAM file " << fname << endl;
					myBams.push_back(new MyBam(fname));
					myBamsFileNames.push_back(fname);
				}
			}
		}
		file.close();
	}
}

DetInDel::~DetInDel()
{
	if (params.alignAgainstReference && fai) {
		fai_destroy(fai);
	}
	for (size_t b=0;b<myBams.size();b++) delete myBams[b];
}

void DetInDel::analyzeDifference(const pair<Haplotype, Haplotype> & hp1, const vector<Read> & reads,  uint32_t leftPos, uint32_t rightPos)
{
	cout << "Inference results" << endl;

	if (params.analyzeLowFreqDiffThreshold<-100.0) {
		size_t offset=50;
		cout << "Haplotype pair: " << endl;
		cout << hp1.first << endl << hp1.second << endl;


		cout << "h.1 alignment: " << endl;
		cout << string(offset,' ') << hp1.first.seq << endl;
		for (size_t r=0;r<reads.size();r++) {
			ObservationModelFBMax om(hp1.first, reads[r], leftPos,params.obsParams);
			MLAlignment ml=om.calcLikelihood();

			double lm=ml.ll;
			ObservationModelFBMax op(hp1.second, reads[r], leftPos, params.obsParams);
			MLAlignment ml2=op.calcLikelihood();
			double lp=ml2.ll;

			if (lm>=lp) om.printAlignment(offset);
		}

		cout << "h.2 alignment: " << endl;
		cout << string(offset,' ') << hp1.second.seq << endl;
		for (size_t r=0;r<reads.size();r++) {
			ObservationModelFBMax om(hp1.first, reads[r], leftPos,params.obsParams);
			MLAlignment ml=om.calcLikelihood();

			double lm=ml.ll;
			ObservationModelFBMax op(hp1.second, reads[r], leftPos, params.obsParams);
			MLAlignment ml2=op.calcLikelihood();
			double lp=ml2.ll;

			if (lp>=lm) op.printAlignment(offset);


		}



	} else {

		double ll=0.0,l1=0.0, l2=0.0;
		vector<size_t > show;
		for (size_t r=0;r<reads.size();r++) {
			ObservationModelFBMax om(hp1.first, reads[r], leftPos,params.obsParams);
			MLAlignment ml=om.calcLikelihood();

			double lm=ml.ll;
			ObservationModelFBMax op(hp1.second, reads[r], leftPos, params.obsParams);
			MLAlignment ml2=op.calcLikelihood();
			double lp=ml2.ll;
			double dll=log(exp(lp)+exp(lm))+log(.5);
			l1+=(addLogs(lm,lm)+log(.5));
			l2+=(addLogs(lp,lp)+log(.5));
			ll+=dll;

			if (lp-lm>params.analyzeLowFreqDiffThreshold) {
				show.push_back(r);
			}
	//		cout << "read[" << r <<"]: 1-mq: " << 1.0-reads[r].mapQual << " first hap lik: " << lm << " second hap lik: " << lp << " combined: " << dll << " lp+lm: " << ll << " lm+lm: " << l1 << " lp+lp: " << l2 << endl;

		}

		size_t offset=50;
		if (show.size()) {
			cout << "Haplotype pair: " << endl;
			cout << hp1.first << endl << hp1.second << endl;


			cout << "h.1 alignment: " << endl;
			cout << string(offset,' ') << hp1.first.seq << endl;

			for (size_t i=0;i<show.size();i++) {
				Read rr=reads[show[i]];
				ObservationModelFBMax om2(hp1.first, rr, leftPos, params.obsParams);
				om2.calcLikelihood();

		//		om2.computeMarginals();
				om2.printAlignment(offset);
			}

			cout << endl << endl;
			cout << "h.2 alignment: " << endl;
			cout << string(offset,' ') << hp1.second.seq << endl;

			for (size_t i=0;i<show.size();i++) {
				Read rr=reads[show[i]];
				ObservationModelFBMax om2(hp1.second, rr, leftPos, params.obsParams);
				om2.calcLikelihood();
		//		om2.computeMarginals();
				om2.printAlignment(offset);
			}
		}
		else { cout << "No differences in log-likelihoods over threshold." << endl; };
	}

}

void DetInDel::showAlignments(const pair<Haplotype, Haplotype> & hp1, const vector<Read> & reads,  uint32_t leftPos, uint32_t rightPos)
{
	cout << "Inference results" << endl;

	double ll=0.0;
	int offset=50;
	vector <double> lf(reads.size(),0.0), ls(reads.size(),0);
	cout << "h.1 alignment: " << endl;
	cout << string(offset,' ') << hp1.first.seq << endl;
	for (size_t r=0;r<reads.size();r++) {
		ObservationModelFBMax om(hp1.first, reads[r], leftPos,params.obsParams);
		double lm=om.getLogLikelihood();
		lf[r]=lm;
		if (lm<params.analyzeLowFreqDiffThreshold) {
			om.printAlignment(offset);
		}
	}

	cout << "h.2 alignment: " << endl;
	cout << string(offset,' ') << hp1.second.seq << endl;
	for (size_t r=0;r<reads.size();r++) {
		ObservationModelFBMax om(hp1.second, reads[r], leftPos,params.obsParams);
		double lm=om.getLogLikelihood();
		ls[r]=lm;
		if (lm<params.analyzeLowFreqDiffThreshold) {
			om.printAlignment(offset);
		}
		ll+=addLogs(lf[r],ls[r])+log(.5);
	}
	cout << "Total loglikelihood: " << ll << endl;


}

void DetInDel::showAlignmentsPerHaplotype(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, uint32_t candPos, uint32_t leftPos)
{
	cout << "ALIGNMENTS" << endl;

	vector<std::set<size_t> > maxHap(haps.size());
	for (size_t r=0;r<reads.size();r++) {
		size_t idx=0;
		double ml=-HUGE_VAL;
		for (size_t h=0;h<haps.size();h++) {
			if (liks[h][r].ll>ml) {
				ml=liks[h][r].ll;
				idx=h;
			}
		}
		maxHap[idx].insert(r);
	}

	int offset=50;
	for (size_t h=0;h<haps.size();h++) {
		cout << "*******************************************" << endl;
		cout << endl << "HAPLOTYPE " << h << endl << endl;
		cout << string(offset,' ') << haps[h].seq << endl;
		BOOST_FOREACH(size_t r, maxHap[h]) {
			ObservationModelFBMax om(haps[h], reads[r], leftPos,params.obsParams);
			om.calcLikelihood();
			om.printAlignment(offset);
		}
	}


}




string DetInDel::getRefSeq(uint32_t lpos, uint32_t rpos)
{
	if (!fai) throw string("FAI error.");

	char *str;
	char *ref;

	str = (char*)calloc(strlen(params.tid.c_str()) + 30, 1);
	sprintf(str, "%s:%d-%d", params.tid.c_str(), lpos, rpos);
	int len;
	ref = fai_fetch(fai, str, &len);
	if (len==0) throw string("faidx error: len==0");
	free(str);
	string res(ref);
	free(ref);

	transform(res.begin(), res.end(), res.begin(), ::toupper);
	return res;
}


double DetInDel::getMaxHap(Haplotype & h1, Haplotype &h2, HapPairLik & hpl, const vector<Haplotype> & haps, vector<HapPairLik> & likPairs)
{

	size_t idx=0, midx;
	double maxll=-HUGE_VAL;
	for (idx=0;idx<likPairs.size();idx++) {
			double ll=likPairs[idx].ll;
			if (ll>maxll) {
				maxll=ll;
				midx=idx;
			}
	}
	h1=haps[likPairs[midx].h1];
	h2=haps[likPairs[midx].h2];

	/*
	cout << "getMaxHap: " << midx <<  " h1: " << likPairs[midx].h1 << " h2: " << likPairs[midx].h2 << endl;
	cout << "indelcoverage h1: ";
	for (map<int, VariantCoverage>::const_iterator it=likPairs[midx].hapIndelCoverage1.begin();it!=likPairs[midx].hapIndelCoverage1.end();it++) {
		cout << "[" << it->second.nf << "," << it->second.nr << "]";
	}
	cout << endl;
	cout << "indelcoverage h2: ";
	for (map<int, VariantCoverage>::const_iterator it=likPairs[midx].hapIndelCoverage2.begin();it!=likPairs[midx].hapIndelCoverage2.end();it++) {
		cout << "[" << it->second.nf << "," << it->second.nr << "]";
	}
	cout << endl;
	*/

	hpl=likPairs[midx];
	return maxll;
}

void DetInDel::outputMaxHap(ostream *output, const string & prefix, const vector<Haplotype> & haps, vector<HapPairLik> & likPairs)
{

	Haplotype h1, h2;
	HapPairLik hpl;
	getMaxHap(h1,h2, hpl, haps, likPairs);
	*output << prefix << " " << hpl.ll << " " << hpl.numFirst << " " << hpl.numSecond << " " << hpl.numIndFirst << " " << hpl.numIndSecond << " " << hpl.numOffBoth << " " << h1.seq << " " << h2.seq << " ";
	for (map<int, AlignedVariant>::const_iterator it=h1.indels.begin();it!=h1.indels.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
	*output << "!";
	for (map<int, AlignedVariant>::const_iterator it=h2.indels.begin();it!=h2.indels.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
	*output << "!";
	for (map<int, AlignedVariant>::const_iterator it=h1.snps.begin();it!=h1.snps.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
	*output << "!";
	for (map<int, AlignedVariant>::const_iterator it=h2.snps.begin();it!=h2.snps.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
	*output << endl;


}

void DetInDel::outputTopHaps(ostream *output, const string & prefix, const vector<Haplotype> & haps, vector<HapPairLik> & likPairs, int n)
{
	// output n most likely haplotype pairs

	for (int ns=0;ns<n && ns<int(likPairs.size());ns++) {
		const Haplotype & h1 = haps[likPairs[ns].h1];
		const Haplotype & h2 = haps[likPairs[ns].h2];
		const HapPairLik & hpl = likPairs[ns];
		*output << prefix << " " << ns+1 << " " << hpl.ll << " " << hpl.numFirst << " " << hpl.numSecond << " " << hpl.numIndFirst << " " << hpl.numIndSecond << " " << hpl.numOffBoth << " " << h1.seq << " " << h2.seq << " ";
		for (map<int, AlignedVariant>::const_iterator it=h1.indels.begin();it!=h1.indels.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
		*output << "!";
		for (map<int, AlignedVariant>::const_iterator it=h2.indels.begin();it!=h2.indels.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
		*output << "!";
		for (map<int, AlignedVariant>::const_iterator it=h1.snps.begin();it!=h1.snps.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
		*output << "!";
		for (map<int, AlignedVariant>::const_iterator it=h2.snps.begin();it!=h2.snps.end();it++) if (it->second.getString()!="*REF") *output << "[" << it->first << "," << it->second.getStartRead() << "," << it->second.getString() << "]";
		*output << endl;
	}

}

void DetInDel::outputHapsAndFreqs(ostream *output, const string & prefix, const vector<Haplotype> & haps, const vector<double> & freqs, uint32_t leftPos)
{
	// output n most likely haplotype pairs

	for (size_t h=0;h<haps.size();h++) {
		const Haplotype & h1 = haps[h];
		*output << prefix << " " << h+1 << " " << freqs[h] << " ";
		for (map<int, AlignedVariant>::const_iterator it=h1.indels.begin();it!=h1.indels.end();it++) if (it->second.getString()!="*REF") *output << leftPos+it->first << "," << it->second.getString() << "|";
		//*output << "!";
		//for (map<int, AlignedVariant>::const_iterator it=h1.snps.begin();it!=h1.snps.end();it++) if (it->second.getString()!="*REF") *output << leftPos+it->first << "," << it->second.getString() << "|";
		*output << endl;
	}

}



void DetInDel::empiricalDistributionMethod(int index, const vector<Read> & reads, uint32_t pos, uint32_t leftPos, uint32_t rightPos, const AlignedCandidates & candidateVariants,  OutputData & oData, OutputData & glfData)
{
	vector<Haplotype> haps;

	vector<vector<MLAlignment> > liks;
	vector<HapPairLik> likPairs;



	// get the haplotypes
	// changes leftPos and rightPos to haplotype blocks in HDIterator

	// NOTE leftPos will be the left position of the reference sequence each haplotype will be aligned to
	bool skip=getHaplotypes(haps, reads, pos, leftPos, rightPos, candidateVariants);

	if (int(reads.size()*haps.size())>params.maxHapReadProd) {
		stringstream os;
		os << "skipped_numhap_times_numread>" << params.maxHapReadProd;
		throw os.str();
	}

	int refSeqPos=leftPos;

	if (skip) {
		cerr << "tid: " << params.tid << " pos:  " << pos << " SKIPPING!" << endl;
	}
	else {

		if (!params.quiet) cout << "[empiricalDistributionMethod] Number of haplotypes: " << haps.size() << endl;
		// compute likelihood of every read given every haplotype

		if (params.estimateHapFreqs) {
			vector<double> hapFreqs;
			map <int, vector<tuple<AlignedVariant, double, double> > > posteriors;
			vector<HapEstResult> her;

			OutputData::Line prefilledLine(oData);
			prefilledLine.set("index", index);
			prefilledLine.set("tid", params.tid);
			prefilledLine.set("center_position",pos);
			prefilledLine.set("num_reads", reads.size());
			prefilledLine.set("msg","ok");
			// string rseq = getRefSeq(leftPos+1, rightPos+1);
			prefilledLine.set("lpos",leftPos);
			prefilledLine.set("rpos",rightPos);
			// prefilledLine.set("refseq", rseq);


			vector<int> onHap(reads.size(),1); // which reads were mapped inside the haplotype window given an artificially high mapping quality

			if (params.slower) {
				computeLikelihoods(haps, reads, liks, leftPos, rightPos, onHap);
			} else {
				computeLikelihoodsFaster(haps, reads, liks, leftPos, rightPos, onHap);
				// int nrOffAll=0; // number of reads mapped outside all haplotypes
				// for (size_t x=0;x<onHap.size();x++) if (!onHap[x]) nrOffAll++;
			}


			int numReadOffAllHaps=0;
			int numHQReads=0;
			for (size_t r=0;r<reads.size();r++) if (reads[r].mapQual > (1.0-1e-6) ) {
				numHQReads++;
				int offall=1;
				for (size_t h=0;h<haps.size();h++) if (!liks[h][r].offHap) offall=0;
				if (offall) numReadOffAllHaps++;
			}
			prefilledLine.set("num_off_hap", numReadOffAllHaps);
 			prefilledLine.set("num_hqreads", numHQReads);


			//estimateHaplotypeFrequenciesPosterior(haps, reads,liks, hapFreqs, posteriors, pos, leftPos, glfOutput);
			/*
			estimateHaplotypeFrequenciesBayesEM(haps, reads,liks, hapFreqs, her, pos, leftPos, glfOutput, "all");
			BOOST_FOREACH(HapEstResult hr, her) {
				//cout << "EMA " << params.tid << " "  << pos << " " << leftPos+hr.pos << " " <<reads.size() << " " << hr.av.getString() << " " << hr.prob << " " << hr.freq << " " << hr.freq*double(reads.size()) << " " << hr.nrf << " " << hr.nrr << endl;
				OutputData::Line EMALine = prefilledLine;
				EMALine.set("analysis_type", "EMA");
				EMALine.set("realigned_position", leftPos+hr.pos);
				EMALine.set("first_called_all", hr.av.getString());
				EMALine.set("post_prob_variant", hr.prob);
				EMALine.set("est_freq", hr.freq);
				EMALine.set("first_var_cover_forward", hr.nrf);
				EMALine.set("first_var_cover_reverse", hr.nrr);
				oData.output(EMALine);
			}
			*/

			hapFreqs.clear();
			her.clear();
			estimateHaplotypeFrequenciesBayesEM(haps, reads,liks, hapFreqs, her, pos, leftPos, rightPos, glfData, index, candidateVariants, params.bayesType);

			for (size_t x=0;x<her.size();x++) {
				HapEstResult & hr = her[x];
				//cout << "EMSV " << params.tid << " "  << pos << " " << leftPos+hr.pos << " " <<reads.size() << " " << hr.av.getString() << " " << hr.prob << " " << hr.freq << " " << hr.freq*double(reads.size()) << " " << hr.nrf << " " << hr.nrr << endl;

				int var_in_window=0;
				const AlignedVariant & avar = hr.av;
				const AlignedVariant *av = candidateVariants.findVariant(hr.pos+leftPos, avar.getType(), avar.getString());
				if (av!=NULL) {
					var_in_window=1;
				}

				OutputData::Line EMSVLine = prefilledLine;
				EMSVLine.set("analysis_type", params.bayesType);
				EMSVLine.set("realigned_position", leftPos+hr.pos);
				EMSVLine.set("first_called_all", hr.av.getString());
				EMSVLine.set("post_prob_variant", hr.prob);
				EMSVLine.set("est_freq", hr.freq);
				EMSVLine.set("first_var_cover_forward", hr.nrf);
				EMSVLine.set("first_var_cover_reverse", hr.nrr);
				EMSVLine.set("was_candidate_in_window",var_in_window);
			//	oData.output(EMSVLine);
			}




			if (params.outputRealignedBAM && params.slower) {
				stringstream os;
				os << index << "_" << params.tid << "_" << leftPos+params.minReadOverlap << "_" << rightPos-params.minReadOverlap  << ".bam";

				vector < CIGAR > cigars(reads.size());
				for (size_t r=0;r<reads.size();r++) {
					if (onHap[r]) {
						double llmax=-HUGE_VAL;
						int hidx=0;
						for (size_t h=0;h<haps.size();h++) if (liks[h][r].ll>llmax) {
							llmax=liks[h][r].ll;
							hidx=h;
						}
						cigars[r]=getCIGAR(haps[hidx], reads[r], liks[hidx][r], refSeqPos);
					}
				}
				int leftOk = leftPos + params.minReadOverlap;
				int rightOk = rightPos - params.minReadOverlap;
	
				string newBAMFileName=params.fileName;
				newBAMFileName.append(".ra.").append(os.str());
				writeRealignedBAMFile(newBAMFileName, cigars, reads, onHap, myBams[0]->bh);

				if (params.processRealignedBAM != "no")  {
					stringstream cmd;
					cmd << params.processRealignedBAM << " " << newBAMFileName << " " << params.fileName << "_realigned" << " " << params.tid << " " << leftOk << " " << rightOk;
					cout << "Executing: " << cmd.str() << endl;
					system(cmd.str().c_str());
				}


				/*
				newBAMFileName=params.fileName;
				newBAMFileName.append(".ua.").append(os.str());
				writeUnalignedBAMFile(newBAMFileName, reads, onHap, myBams[0].bh);
				*/
			}

		}
		if (params.showHapAlignments) {
			showAlignmentsPerHaplotype(haps, reads, liks, pos, leftPos);
		}

		if (params.doDiploid) {
		//	cout << "A" << endl;
			vector<double> hapFreqs;
			map <int, vector<tuple<AlignedVariant, double, double> > > posteriors;
			vector<HapEstResult> her;

			vector<int> onHap(reads.size(),1); // which reads were mapped inside the haplotype window given an artificially high mapping quality

			if (params.slower) {
				computeLikelihoods(haps, reads, liks, leftPos, rightPos, onHap);
			} else {
				computeLikelihoodsFaster(haps, reads, liks, leftPos, rightPos, onHap);
			}

			diploidGLF(haps, reads, liks, hapFreqs, her, pos, leftPos, rightPos, glfData, index,  candidateVariants,"dip");


			//statisticsHaplotypePair(haps,reads,liks, hpl, prefilledLine);

			/*
			string prefix("FMAP");


			callIndel(haps, reads, liks, likPairs, pos, leftPos, rightPos, prefix, prefilledLine, oData);
			callSNP(haps, reads, liks, likPairs, pos, leftPos, rightPos, prefix, prefilledLine, oData);
			*/
			/*
			if (!(nci==0 && ncs==0 && params.printCallsOnly)) {
				stringstream os; os << " " << params.tid << " " << pos << " " << leftPos; prefix.append(os.str());
				outputTopHaps(&output, prefix, haps, likPairs, params.numOutputTopHap);
			}
			*/

			/*
			if (params.analyzeLowFreq) {
				pair<Haplotype, Haplotype> oh;
				oh.first=h2;
				oh.second=h1;
				analyzeDifference(oh, reads, leftPos, rightPos);

				oh.first=h1;
				oh.second=h2;
				analyzeDifference(oh, reads, leftPos, rightPos);
				//oh.first.printHaps();
				//oh.second.printHaps();
			}
			*/

			if (params.outputRealignedBAM && params.slower) {
				// computes pair likelihoods using priors
				computePairLikelihoods(haps, reads, liks, likPairs, true,candidateVariants, leftPos);
				Haplotype h1, h2; HapPairLik hpl;
				getMaxHap(h1, h2, hpl, haps, likPairs);

				vector < CIGAR > cigars(reads.size());
				for (size_t r=0;r<reads.size();r++) {
						int hmax = hpl.h1;
						if (fabs(liks[hpl.h1][r].ll-liks[hpl.h2][r].ll)<1e-8) {
							if (haps[hpl.h1].countIndels()<haps[hpl.h2].countIndels()) hmax = hpl.h1; else hmax = hpl.h2;
						} else {
							if (liks[hpl.h1][r].ll>liks[hpl.h2][r].ll) {
								hmax = hpl.h1;
							} else {
								hmax = hpl.h2;
							}
						}
						const Haplotype & hx = haps[hmax];
						const Read & rd = reads[r];
						const MLAlignment & ml = liks[hmax][r];
						cigars[r]=getCIGAR(haps[hmax], reads[r], liks[hmax][r], refSeqPos);
				}

				stringstream os;
				int leftOk = leftPos + params.minReadOverlap;
				int rightOk = rightPos - params.minReadOverlap;
				os << index << "_" << params.tid << "_" << leftPos+params.minReadOverlap << "_" << rightPos-params.minReadOverlap  << ".bam";
				string newBAMFileName=params.fileName;
				newBAMFileName.append(".ra.").append(os.str());
				writeRealignedBAMFile(newBAMFileName, cigars, reads, onHap, myBams[0]->bh);
				/*
				newBAMFileName=params.fileName;
				newBAMFileName.append(".ua.").append(os.str());
				writeUnalignedBAMFile(newBAMFileName, reads, onHap, myBams[0]->bh);
				*/
				if (params.processRealignedBAM != "no")  {
					stringstream cmd;
					cmd << params.processRealignedBAM << " " << newBAMFileName << " " << params.fileName << "_realigned" << " " << params.tid << " " << leftOk << " " << rightOk;
					cout << "Executing: " << cmd.str() << endl;
					system(cmd.str().c_str());
				}


			}


		}
//		glf.writeToFile(string(""), output);
	}

}



void DetInDel::writeUnalignedBAMFile(const string & fileName, const vector<Read> & reads, const vector<int> & onHap, const bam_header_t *bh=NULL)
{

	if (onHap.size()!=reads.size()) return;
	bool hasUnaligned=false;
	for (size_t x=0;x<onHap.size();x++) if (!onHap[x]) {
		hasUnaligned=true;
		break;
	}
	if (!hasUnaligned) return;

	bamFile bf = bam_open(fileName.c_str(),"wb");
	if (bf==NULL) throw string("Cannot open bamfile ").append(fileName).append(" for writing!");

	if (bh!=NULL) {
		bam_header_write(bf, bh);
	}

	for (size_t r=0;r<reads.size();r++) if (!onHap[r]) {
		bam1_t *b=reads[r].bam;
		if (bam_write1(bf,b)<=0) throw string("Error writing to unalignable read to bamfile.");
	}

	bam_close(bf);
}

void DetInDel::writeRealignedBAMFile(const string & fileName, const vector<CIGAR> & cigars, const vector<Read> & reads, const vector<int> & onHap, const bam_header_t *bh=NULL)
{
		if (cigars.size()!=reads.size()) throw string("Problem with the cigars.");

		bamFile bf = bam_open(fileName.c_str(),"wb");
		if (bf==NULL) throw string("Cannot open bamfile ").append(fileName).append(" for writing!");

		if (bh!=NULL) {
				bam_header_write(bf, bh);
		}

		for (size_t r=0;r<reads.size();r++) {

				bam1_t *b=reads[r].bam;

				if (onHap[r]) {
						bam1_t *nb=bam_init1();

						uint32_t old_ncig=b->core.n_cigar;
						uint32_t new_ncig=cigars[r].size();
						int old_data_len=b->data_len;
						int new_data_len=old_data_len - old_ncig*4 + new_ncig*4;

						*nb=*b;
						nb->data = (uint8_t*)calloc(new_data_len, 1);
						nb->data_len=new_data_len;
						nb->m_data=nb->data_len;
						// copy cigar
						for (uint32_t n=0;n<new_ncig;n++) {
								bam1_cigar(nb)[n]=(  (((uint32_t) cigars[r][n].second) << BAM_CIGAR_SHIFT) | ( ( (uint32_t) cigars[r][n].first ) )  );
						}
						nb->core.n_cigar=(unsigned int) new_ncig;

						for (size_t n=0;n<b->core.l_qname;n++) {
								nb->data[n]=b->data[n];
						}

						int y=b->core.l_qname+4*new_ncig;
						for (int n=b->core.l_qname+4*old_ncig;n<b->data_len;n++,y++) {
								nb->data[y]=b->data[n];
						}

						// update position of read

						nb->core.pos=cigars[r].refPos;
						// update insert size if mapped
						nb->core.isize=cigars[r].refPos-nb->core.mpos;
						if (bam_write1(bf,nb)<=0) throw string("Error writing alignment to realigned BAM file.");
						bam_destroy1(nb);
				} else {
						if (bam_write1(bf,b)<=0) throw string("Error writing alignment to realigned BAM file.");
				}
		}

		bam_close(bf);
}


DetInDel::CIGAR DetInDel::getCIGAR(const Haplotype & hap, const Read & read, const MLAlignment & ml, int refSeqStart)
{
	if (hap.ml.hpos.size()!=hap.size()) throw string("Haplotype has not been aligned!");
	if (ml.hpos.size()!=read.size()) throw string("Read is not properly aligned!");
	const MLAlignment & hml=hap.ml; // alignment of haplotype to reference

	//string qname = bam1_qname(read.getBam());
	const int debug = 0;
	/*
	if (qname == "IL8_4337:8:102:11530:1494") {
		cout << "YES" << endl;
		cout << qname << endl;
		debug = 1;
	}
	*/

	vector<int> npos(read.size()); // npos records position of read base on the reference sequence
	for (int b=0;b<int(read.size());b++) {
		if (ml.hpos[b]>=0) npos[b]=hml.hpos[ml.hpos[b]]; else npos[b]=ml.hpos[b];
	}

	if (debug) {
		for (size_t h=0;h<npos.size();h++) {
			cout << "[" << h << "," << npos[h] << "]";
		}
		cout << endl;
		cout << endl;

		for (size_t h=0;h<ml.hpos.size();h++) {
			cout << "[" << h << "," << ml.hpos[h] << "]";
		}
		cout << endl;
		cout << endl;
		for (size_t h=0;h<hml.hpos.size();h++) {
				cout << "[" << h << "," << hml.hpos[h] << "]";
			}
			cout << endl;
	}
	CIGAR cig;

	int b=0, prevponr=0; // position on reference previous base aligned on the reference (ie no deletion/insertion/LO/RO)

	// determine last base in read aligned to the haplotype
	b=read.size()-1;
	while (npos[b]<0) b--;
	int lastbonh=b;

	if (lastbonh<0) {
		// clip the whole read, read is de facto off haplotype
		cig.push_back(CIGAR::CIGOp(BAM_CSOFT_CLIP, read.size()));
		return cig;
	}


	if (debug) {
		cout << "lastbonh: " << lastbonh << endl;
	}
	// find first base in read that is aligned to haplotype and to reference
	// all sequence before that is considered 'soft clipped', ie not aligned. This may include sequence that matches perfectly to the reference
	b=0;
	while (npos[b]<0) b++;
 	if (b>0) cig.push_back(CIGAR::CIGOp(BAM_CSOFT_CLIP, b));
 	prevponr=npos[b];
 	cig.refPos=refSeqStart+prevponr;

 	int curr_cop=BAM_CMATCH;
 	int len_curr_cop=1;



 	while (b<lastbonh) {

 		int chp=npos[b]; // position on reference of current base in read
 		int nhp=npos[b+1];

 		if (nhp==MLAlignment::INS) {
 			if (chp==MLAlignment::INS) {
 				// stay on inserted sequence
 				if (curr_cop!=BAM_CINS) throw string("Error(1)!");
 				len_curr_cop++;
 			} else {
 				if (chp>=0) {
 					// going from on reference to insertions
 					if (curr_cop!=BAM_CMATCH) throw string("Error(2)!");
 					// write CIGAR
					cig.push_back(CIGAR::CIGOp(BAM_CMATCH,len_curr_cop));

					// update current CIGAR operation
					len_curr_cop=1;
					curr_cop=BAM_CINS;

					prevponr=chp;
 				} else throw string("How is this possible? (1)");
 			}

 		} else if (chp>=0 && nhp>=0 && nhp-chp==1) {
 	 	  // MATCH to MATCH
 			if (curr_cop!=BAM_CMATCH) {
 				cout << "b: " << b << " chp: " << chp << " nhp: " << nhp << endl;
 				throw string("Error(3)!");
 			}
 			len_curr_cop++;
 			prevponr=nhp;
 		} else if (chp>=0 && nhp>=0 && nhp-chp>1) {
 			// deletion
 			if (curr_cop!=BAM_CMATCH) throw string("Error(4)!");
 			// write CIGAR
 			cig.push_back(CIGAR::CIGOp(BAM_CMATCH,len_curr_cop));

 			// write deletion CIGAR
 			cig.push_back(CIGAR::CIGOp(BAM_CDEL,nhp-chp-1));

 			curr_cop=BAM_CMATCH;
 			len_curr_cop=1;

 			prevponr=nhp;
 		} else if (chp==MLAlignment::INS && nhp-prevponr==1) {
 			// from inserted bases to bases matched to the reference
 			cig.push_back(CIGAR::CIGOp(BAM_CINS,len_curr_cop));

 			//
 			curr_cop=BAM_CMATCH;
 			len_curr_cop=1;

 			prevponr=nhp;
 		} else if (chp==MLAlignment::INS && nhp-prevponr>1) {
 			// next base is again on reference but more than 1 reference base from the last read base aligned to the haplotype
 			cig.push_back(CIGAR::CIGOp(BAM_CINS,len_curr_cop));
 			cig.push_back(CIGAR::CIGOp(BAM_CDEL,nhp-prevponr-1));

 			curr_cop=BAM_CMATCH;
 			len_curr_cop=1;

 			prevponr=nhp;
 		}
 		b++;
 	}

 	// write last cigar
 	cig.push_back(CIGAR::CIGOp(curr_cop,len_curr_cop));

 	// write soft_clip at the end
 	if (read.size()-1 - lastbonh>0) {
 		cig.push_back(CIGAR::CIGOp(BAM_CSOFT_CLIP,read.size()-1 - lastbonh));
 	}

 	/*
 	cout << "cig: ";
 	BOOST_FOREACH(CIGAR::CIGOp cop, cig) {
 		cout << "(" << cop.first << "," << cop.second << ")" ;
 	}
 	cout << endl;
	*/
 	return cig;
}


void DetInDel::getReads(uint32_t leftPos, uint32_t rightPos, vector<Read> & reads, uint32_t & oldLeftPos, uint32_t & oldRightFetchReadPos, vector<Read *> & readBuffer, bool reset)
{
	// filter using map quality
	class SortFunc {
	public:
		static bool sortFunc(const Read & r1, const Read & r2)
		{
			// sort in decreasing order
			if (r1.mapQual>r2.mapQual) return true; else return false;
		}
	};

	if (leftPos<oldLeftPos) {
		cerr << "Windows are not sorted!" << endl;
		exit(3);
	}

	reads.clear();

	if (int(rightPos-leftPos)<3*params.minReadOverlap) throw string("Choose a larger width or a smaller minReadOverlap.");


	int maxDev = int ( libraries.getMaxInsertSize());

	//maxDev = 100;
	//cerr << "CHANGE THIS CHANGE THIS" << endl;

	string_hash< list<int> > mapped_name_to_idx, unmapped_name_to_idx; // query name to read idx
	string_hash< list<int> >::const_iterator hash_it;

	int numUnknownLib = 0;
	string_hash <int> unknownLib;
	const int LEFTPAD = 200;
	// note the idea is to get only consider reads starting at position leftMostReadPos (and not ones merely overlapping)
	// LEFTPAD should take care of overlap effects (note that leftMostReadPos is already generous, based on library insert size)


	uint32_t rightFetchReadPos = rightPos+maxDev;
	uint32_t rightMostReadPos = rightPos+maxDev;

	uint32_t leftFetchReadPos = leftPos-maxDev-LEFTPAD;
	uint32_t leftMostReadPos = leftPos-maxDev-LEFTPAD; // left most position of reads we want to seriously consider

	// reset indicates whether we want to remake the readBuffer

	vector<Read*> newReadBuffer;
	bool leftOverlapsPrevious = false;
	if (reset) {
		for (size_t r=0;r<readBuffer.size();r++) {
			if (readBuffer[r]!=NULL) delete readBuffer[r];
		}
		readBuffer.clear();
		oldRightFetchReadPos = rightFetchReadPos;
	} else {
		// clear reads that do not overlap with new window [leftMostReadPos, rightMostReadPos]

		for (size_t r=0;r<readBuffer.size();r++) {
			uint32_t rend = readBuffer[r]->getEndPos();
			uint32_t rbeg = readBuffer[r]->getBam()->core.pos;
			if (rbeg<leftMostReadPos) {
				delete readBuffer[r];
				readBuffer[r]=NULL;
			} else {
				newReadBuffer.push_back(readBuffer[r]);
			}
		}

		// note that it is required that the new leftPos of the window >= the old leftPos
		// therefore if leftFetchReadPos<=oldRightFetchReadPos the new w
		if (leftMostReadPos<oldRightFetchReadPos) {
			leftFetchReadPos = oldRightFetchReadPos;
			leftOverlapsPrevious = true;
		}
	}




	// cout << "leftFetchReadPos: " << leftFetchReadPos << " rightFetchReadPos: " << rightFetchReadPos << " oldRightFetchReadPos: " << oldRightFetchReadPos << endl;
	// cout << "leftMostReadPos: " << leftMostReadPos << " rightMostReadPos: " << rightMostReadPos << " leftOverlapsPrevious: " << int(leftOverlapsPrevious) << endl;
	// store updated readBuffer
	readBuffer.swap(newReadBuffer);

//	cout << "leftPos : " << leftPos << " rightPos: " << rightPos << " maxDev: " << maxDev << endl;


	// first clean readbuffer


	int numReads = readBuffer.size();

	vector<Read> newReads;
	if (leftFetchReadPos<=rightFetchReadPos) {
		cout << "Fetching reads...." << endl;
		for (size_t b=0;b<myBams.size();b++) {
			//bam_fetch(myBams[b].bf, myBams[b].idx, myBams[b].getTID(params.tid), leftPos+params.minReadOverlap, rightPos-params.minReadOverlap, &reads, &Read::fetchFuncVector);
			Read::FetchReadData data(&newReads, int(b), &(this->libraries), &myBams, numReads, params.maxReads*100);
			bam_fetch(myBams[b]->bf, myBams[b]->idx, myBams[b]->getTID(params.tid), leftFetchReadPos , rightFetchReadPos,&data , &Read::fetchFuncVectorPooled);
			numUnknownLib += data.numUnknownLib;
			numReads = data.numReads;
		}
		oldRightFetchReadPos = rightFetchReadPos;
	}

	// add new reads to readBuffer

	for (size_t r=0;r<newReads.size();r++) {
		if (newReads[r].getBam()->core.pos>=leftFetchReadPos) {
			// only store reads that do not overlap with the boundary;
			// reads overlapping with boundary will have been picked up before.
			readBuffer.push_back(new Read(newReads[r]));
		}
	}

	if (0) {
		// check with regular fetch
		vector<Read> tmpReads;
		for (size_t b=0;b<myBams.size();b++) {
			//bam_fetch(myBams[b].bf, myBams[b].idx, myBams[b].getTID(params.tid), leftPos+params.minReadOverlap, rightPos-params.minReadOverlap, &reads, &Read::fetchFuncVector);
			Read::FetchReadData data(&tmpReads, int(b), &(this->libraries), &myBams, numReads, params.maxReads*100);
			bam_fetch(myBams[b]->bf, myBams[b]->idx, myBams[b]->getTID(params.tid), leftMostReadPos , rightMostReadPos,&data , &Read::fetchFuncVectorPooled);
		}
		for (size_t r=0;r<tmpReads.size();r++) {
			if (tmpReads[r].getBam()->core.pos>=leftMostReadPos) {
				string qname = string(bam1_qname(tmpReads[r].getBam()));
				cout << "glp: "  << leftPos << " qname: " << qname << " pos: " << tmpReads[r].pos << " end: " << tmpReads[r].getEndPos() << endl;
			}
		}
	}



	// check readbuffer for duplicates (debugging)
	if (1) {
		string_hash <int> qnameCount;
		for (size_t r=0;r<readBuffer.size();r++) {
			string qname = string(bam1_qname(readBuffer[r]->getBam()));
			//cout << "lp: "  << leftPos << " qname: " << qname << " pos: " << readBuffer[r]->pos << " end: " << readBuffer[r]->getEndPos() << endl;
			string_hash<int>::iterator it = qnameCount.find(qname);
			if (it == qnameCount.end()) {
				qnameCount[qname]=1;
			} else {
				qnameCount[qname]++;
				if (qnameCount[qname]>2) {
					cerr << "Duplicate reads: readbuffer problem!" << endl;
					throw string("duplicate reads!");
				}
			}
		}
	}



	newReads.clear();

	size_t oldNumReads=readBuffer.size();


	// copy readBuffer to reads

	for (size_t r=0;r<readBuffer.size();r++) {
		reads.push_back(Read(*readBuffer[r]));
	}


	// get query names

	vector<int> unmapped;
	for (size_t r=0; r<reads.size();r++) {
		if (reads[r].isUnmapped())  {
			unmapped.push_back(r);
			unmapped_name_to_idx[ string(bam1_qname(reads[r].getBam())) ].push_back(r);
		//	 cout << " __reads[" << r  << "]: " << bam1_qname(reads[r].getBam()) << " UNMAPPED" << endl;
		} else {
			mapped_name_to_idx[ string(bam1_qname(reads[r].getBam())) ].push_back(r);
		//	cout << " __reads[" << r  << "]: " << bam1_qname(reads[r].getBam()) << " pos: " << reads[r].pos << " " << reads[r].getBam()->core.pos << " mpos: " << reads[r].getBam()->core.mpos << " mu: " << reads[r].mateIsUnmapped() << endl;
		}
	}


	// filter reads based on haplotype window, the minimum read overlap for mapped reads, minimum mapping quality and maximum read length

	/*
	cout << "name_to_idx.size: " << mapped_name_to_idx.size() << endl;
	for (hash_it = mapped_name_to_idx.begin();hash_it!=mapped_name_to_idx.end();hash_it++) {
		cout << "hit: " << hash_it->first;
		BOOST_FOREACH(int x, hash_it->second) {
			cout << " " << x;
		}
		cout << endl;
	}
	*/
	int numTIDmismatch = 0, numOrphan =0, numOrphanUnmapped = 0, numInRegion = 0;

	// reads are filtered by setting mapping quality to -1
	vector<Read> filteredReads;
	double minMapQual = params.mapQualThreshold;
	if (minMapQual<0.0) minMapQual=0.0;
	for (int r=0;r<int(reads.size());r++) {
		//cout << "***" << endl;
		//cout << "reads.mapQual " << reads[r].mapQual <<  " bq: " << reads[r].getBam()->core.qual << endl;
		bool filter = false;
		int tf = 0;
		if (reads[r].size()>params.maxReadLength) filter=true;

		if (reads[r].getEndPos()<leftMostReadPos || reads[r].pos>rightMostReadPos) filter=true;

		if (!reads[r].isUnmapped()) {
	//		cout << "mapped" << endl;
			if (reads[r].pos+int(reads[r].size())<int(leftPos)+params.minReadOverlap || reads[r].pos>int(rightPos)-params.minReadOverlap) {
		//		cout << " { " << reads[r].pos+reads[r].size() << " " << leftPos+params.minReadOverlap << " " << reads[r].pos << " " << rightPos-params.minReadOverlap << " } " << endl;
				filter=true;
				tf = 1;
			} else if (reads[r].mateIsUnmapped() == false ){
				if (reads[r].getBam()->core.mtid != reads[r].getBam()->core.tid) {
					// filter = true;
					// cout << "TIDERR: reads[" << r << "]: " << bam1_qname(reads[r].getBam()) << " matePos: " << reads[r].matePos << " mateLen: " << reads[r].mateLen << endl;
					numTIDmismatch++;
				} else {
					// find mate of mapped read
					// filter if we cannot find it (mapped to another chromosome, those are a bit suspicious)

					// lookup mapped read
					hash_it = mapped_name_to_idx.find(string(bam1_qname(reads[r].getBam())));
					if (hash_it == mapped_name_to_idx.end()) { numOrphan++; filter=true; } else {
						if (hash_it->second.size()>2) cerr << "HUH? DUPLICATE READ LABELS???" << endl;
						if (reads[r].mateIsUnmapped() == false) {
							filter = true;
						}

						BOOST_FOREACH(int idx, hash_it->second) {
							if (idx != r) {
								reads[r].mateLen = reads[idx].size();
								reads[r].matePos = reads[idx].pos;
								filter = false;
								if (reads[r].matePos != reads[r].getBAMMatePos()) {
									cerr << "matepos inconsistency!" << endl;
									cerr << reads[r].matePos << " " << reads[r].getBAMMatePos() << endl;
									exit(1);
								}
							}
						}

						if (filter == true) {
							numOrphan++;
							tf = 2;
						}
						//cout << "mapped read: " << r << " " << qname << " pos: " << reads[r].pos << " " << reads[r].getBam()->core.mtid << " " <<  reads[r].getBam()->core.mpos << " mateunmap: " << reads[r].mateIsUnmapped() << endl;
					}
				}
			} else if (reads[r].mateIsUnmapped() == true) {
				reads[r].matePos=reads[r].pos;
				hash_it = unmapped_name_to_idx.find(string(bam1_qname(reads[r].getBam())));
				if (hash_it == unmapped_name_to_idx.end()) { filter=true; } else {
					filter = true;
					if (hash_it->second.size()>2) cerr << "HUH? DUPLICATE READ LABELS???" << endl;
					BOOST_FOREACH(int idx, hash_it->second) {
						if (idx != r) {
							reads[r].mateLen = reads[idx].size();
							filter = false;
						}
					}

				}
				if (filter==true) {
					numOrphan++;
					tf = 3;
				}
			}
			if (filter == false) numInRegion++;

		} else {
	//		cout << " unmapped" << endl;
			// read is unmapped
			if (params.mapUnmappedReads) {
	//			cout << " unmapped " << qname << " ; ";

				// lookup mapped read
				hash_it = mapped_name_to_idx.find(string(bam1_qname(reads[r].getBam())));
				int idx;
				if (hash_it == mapped_name_to_idx.end()) { numOrphanUnmapped++; filter=true; tf = 4;} else {
					if (hash_it->second.size()!=1) {
						cerr << "UNMAPPED READ HAS MORE THAN ONE MATE!" << endl;
						exit(1);
					}
					idx = *(hash_it->second.begin());
					//cout << " FOUND " << idx << endl ;
					// check if mate will overlap with haplotype
					uint32_t range_l, range_r; // range of mate

					int maxInsert = (int) reads[idx].getLibrary().getMaxInsertSize();
					int minInsert = 0;
					uint32_t rpos = reads[idx].pos;

	//				cout << "idx: " << idx << " unmapped: " << reads[idx].isUnmapped() << " rpos: " << rpos << " isreverse: " << reads[idx].isReverse() << endl;
					if (reads[idx].isReverse()) {
						range_l = rpos-maxInsert;
						range_r = rpos-minInsert;
					} else {
						range_l = rpos+minInsert;
						range_r = rpos+maxInsert;
					}

	//				cout << "insert: " << insert << " std: " << std << " range_l : " << range_l << " range_r: " << range_r << " leftPos: " << leftPos << "  rightPos: " << rightPos << endl;

					if (range_r>leftPos && range_l<rightPos) {
						numInRegion++;
						filter=false;
						reads[r].mapQual = reads[idx].mapQual;
						reads[r].matePos = reads[idx].pos;
						reads[r].mateLen = reads[idx].size();
						if (reads[r].isReverse() == reads[idx].isReverse()) {
							reads[r].reverse();
							reads[r].complement();
						}
					} else {
						filter=true;
						tf = 5;
					}
				}
			} else {
				filter = true;
			}

		}
		if (filter == true) reads[r].mapQual = -1.0;

//		cout << "reads[" << r << "]: " << bam1_qname(reads[r].getBam()) << " matePos: " << reads[r].matePos << " mateLen: " << reads[r].mateLen << " Filter: " << tf << " filter: " << filter <<  " mq: " << reads[r].mapQual << endl;
	}


	int nUnmapped = 0;
	int nMateposError = 0;
	sort(reads.begin(), reads.end(), SortFunc::sortFunc);
	size_t max; for (max=0;max<params.maxReads && max<reads.size();max++) if (!(reads[max].mapQual<minMapQual)) {
		if (reads[max].matePos==-1 && reads[max].isPaired() && !reads[max].mateIsUnmapped() ) {
			nMateposError++;
			reads[max].matePos = reads[max].pos;
		};
		filteredReads.push_back(Read(reads[max]));
		if (reads[max].isUnmapped()) nUnmapped++;
	} else break;


	filteredReads.swap(reads);
	filteredReads.clear();
	//Read::filterReads(reads, params.maxReads, params.mapQualThreshold, params.maxReadLength, params.minReadOverlap, leftPos, rightPos);

	if (params.filterReadAux!="") {
		if (params.filterReadAux.size()>1) {
			size_t before=reads.size();
			int exclude=1;
			if (params.filterReadAux[0]=='+') exclude=0;
			string match=params.filterReadAux.substr(1,params.filterReadAux.size());
			Read::filterReads(reads, exclude, match);
			size_t after=reads.size();
			if (!params.quiet) cout << "filterAux: " << before-after << " reads were filtered based on match string " << params.filterReadAux << endl;
		}
	}

	if (!params.quiet) cout << "Number of reads: " << reads.size() << " out of " << oldNumReads << " # unmapped reads: " << nUnmapped << " numReadsUnknownLib: " << numUnknownLib << " numChrMismatch: " << numTIDmismatch << " numMappedWithoutMate: " << numOrphan << " numUnmappedWithoutMate: " << numOrphanUnmapped << endl;
	if (nMateposError) {
		cerr << "The mate position of " << nMateposError << " reads was recorded as -1 in the BAM file" << endl;
	}

	if (params.showReads) {
		for (size_t r=0;r<reads.size();r++) {
			cout << "read[" << r << "]: " << reads[r] << endl;
		}
	}

	if (reads.size()<2) {
		throw string("too_few_reads");
	} else if (reads.size()>=params.maxReads) {
		throw string("above_read_count_threshold");
	}

}


void DetInDel::detectIndels(const string & variantsFileName)
{

	ofstream output;
	ofstream glfOutput;

	string callsFile=params.fileName; callsFile.append(".calls.txt");
	string glfFile=params.fileName; glfFile.append(".glf.txt");

	OutputData oData=params.makeOutputData(output);

	/*
	output.open(callsFile.c_str());
	if (!output.is_open()) {
		throw(string("Cannot open file ").append(callsFile).append(" for writing."));
	}

	oData.outputLine(oData.headerString());
	*/


	glfOutput.open(glfFile.c_str());
	if (!glfOutput.is_open()) {
		throw(string("Cannot open file ").append(glfFile).append(" for writing."));
	}

	OutputData glfData=params.makeGLFOutputData(glfOutput);
	glfData.outputLine(glfData.headerString());

	VariantFile vf(variantsFileName);

	int index=0;
	//for (map<uint32_t,InDel>::const_iterator it=indels.begin();it!=indels.end();it++, cnt++) {

	vector<Read *> readBuffer;
	uint32_t oldLeftPos=0, oldRightFetchReadPos=0;


	string oldTid("-1");

	// NOTE ReadBuffer should be reset on first usage or on chromosome change!
	bool resetReadBuffer = true;



	while (!vf.eof()) {
		AlignedCandidates candidateVariants;
		candidateVariants=vf.getLineVector(params.varFileIsOneBased);
		if (candidateVariants.variants.size()==0) continue;

		vector<Read> reads;

		uint32_t pos, leftPos, rightPos;
		// get lowest and highest position

		//leftPos=(candidateVariants.leftPos>int(params.width))?(candidateVariants.leftPos-params.width):0;
		//rightPos=(candidateVariants.rightPos+params.width-1);
		leftPos = candidateVariants.leftPos;
		rightPos = candidateVariants.rightPos;


		pos = candidateVariants.centerPos;
		params.tid=candidateVariants.tid;

		if (params.tid!=oldTid) {
			// reinit
			resetReadBuffer = true;
			oldTid = params.tid;
			oldLeftPos = 0;
		}

		if (leftPos < oldLeftPos) {
			cerr << "leftPos: " << leftPos << " oldLeftPos: " << oldLeftPos << endl;
			cerr << "Candidate variant files must be sorted on left position of window!" << endl;
			exit(1);
		}


// TODO either add tid to AlignedVariant or infer it from the vector of aligned variants
// change alige
		index++;
		bool skipped = false;

		if (!params.quiet) cout << "****" << endl << " tid: " << params.tid << " pos: " << pos << " leftPos: " << leftPos << " " << " rightPos: " << rightPos << endl;

		string message="ok";
		/*
		if (!(indels[pos].count[0]>=params.minCount || indels[pos].count[1]>=params.minCount)) {
			message="below_indel_count_threshold";
			skipped=true;
			goto _end;
		}
		*/



		try {
			getReads(leftPos, rightPos, reads, oldLeftPos, oldRightFetchReadPos, readBuffer, resetReadBuffer);


			if (params.inferenceMethod=="empirical") empiricalDistributionMethod(index, reads, pos, leftPos, rightPos,  candidateVariants, oData, glfData);
			else throw string("Unknown inference method");

		}
		catch (string s) {
			for (size_t x=0;x<s.size();x++) if (s[x]==' ') s[x]='_';
			message=string("error_").append(s);
			skipped=true;
			goto _end;
		}
		catch (std::bad_alloc) {
			message = string("error_bad_alloc");
			skipped = true;
			goto _end;
		}
		catch (std::exception& e) {
			message = string("error_exception_").append(e.what());
			skipped = true;
			goto _end;
		}


		_end:

		if (skipped) {
			cerr << "skipped " << params.tid << " " << pos << " reason: " << message << endl;
			//OutputData::Line line(oData);
			//line.set("msg", message);
			//line.set("index", index);
			//oData.output(line);

			OutputData::Line gline(glfData);
			gline.set("msg", message);
			gline.set("index", index);
			gline.set("tid", params.tid);
			gline.set("lpos", leftPos);
			gline.set("rpos", rightPos);
			glfData.output(gline);

			// reset read buffer: all reads will be fetched again
			resetReadBuffer = true;
		} else {
			resetReadBuffer = false;
		}

		oldLeftPos = leftPos;
	}

	//output.close();
	glfOutput.close();

	// clean up read buffer
	for (size_t r=0;r<readBuffer.size();r++) {
		if (readBuffer[r]!=NULL) delete readBuffer[r];
	}


}




bool DetInDel::alignHaplotypes(vector<Haplotype> & haps,  uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, map<int, std::set<AlignedVariant> > & variants)
{
	uint32_t start=leftPos;
	uint32_t end=rightPos+1;

	variants.clear();

	int print=0;

	seqan::Score<int> score(-1, -460, -100,-960);

	Read rh1;
	rh1.pos=0;
	rh1.posStat.first=0;
	rh1.mapQual=1.0-1e-32;
	ObservationModelParameters alignParams("probabilistic");
	alignParams.pError=0.0001;
	alignParams.pMut=0.01;
	alignParams.maxLengthDel=50;
	alignParams.forceReadOnHaplotype=true;
	alignParams.bMid=0;
	//alignParams.maxLengthIndel=12;
	//alignParams.numIndels=2;
	//alignParams.indelDist="uniform";

	vector<Haplotype> tmp_haps;
	for (size_t h=0;h<haps.size();h++) {
		rh1.seq.seq=haps[h].seq;
		rh1.setAllQual(1.0-1e-16);

		Haplotype hRef;
		uint32_t start=leftPos;
		uint32_t end=rightPos;

		string refSeq=getRefSeq(start+1, end+1);




		hRef.append(refSeq);
		/*
		char lc = (haps[h].seq[haps[h].seq.size()-1]);
		char lcl;
		if (lc == 'T') lcl = 'A'; else if (lc == 'A') lcl = 'T'; else if (lc=='G') lcl = 'C'; else if (lc=='C') lcl = 'G';

		hRef.seq+= lcl;
		*/
		/*
		ObservationModelFBMax om(hRef, rh1, 0, alignParams);
		*/
		ObservationModelSeqAn om(hRef, rh1, 0, alignParams, score);
		haps[h].indels.clear();
		haps[h].snps.clear();
		//om.reportVariants(haps[h].indels, haps[h].snps, haps[h].align);
		//om.calcLikelihood();
		om.align();
		const MLAlignment & ml=om.getMLAlignment();
		haps[h].indels=ml.indels;
		haps[h].snps=ml.snps;
		haps[h].align=ml.align;
		haps[h].ml=ml;
		bool hasStartEndIndel = false;
		if (ml.hpos[0] == MLAlignment::LO) hasStartEndIndel = true;
		int hs = ml.hpos.size()-1;
		if (hs>0 && ml.hpos[hs] == MLAlignment::RO) hasStartEndIndel = true;
		//if (params.showCandHap) {
	//			cout << "hap " << h << endl;om.printAlignment(20);
	//			cout << string(20,' ') << haps[h].align << endl;
	//	}

		for (map<int, AlignedVariant>::const_iterator it=haps[h].indels.begin(); it!=haps[h].indels.end();it++) variants[it->first].insert(it->second);
		for (map<int, AlignedVariant>::const_iterator it=haps[h].snps.begin(); it!=haps[h].snps.end();it++) variants[it->first].insert(it->second);
		if (!hasStartEndIndel) {
			tmp_haps.push_back(haps[h]);
		}

	}

	haps.swap(tmp_haps);

	// add REF allele as variant to each haplotype in order to compute coverage statistics
	for (map<int, std::set<AlignedVariant> >::const_iterator it=variants.begin();it!=variants.end();it++) {
		for (size_t h=0;h<haps.size();h++) haps[h].addRefVariant(it->first);
	}

	if (!params.quiet) {
		for (map<int, std::set<AlignedVariant> >::const_iterator it=variants.begin();it!=variants.end();it++) {
			cout << "aligned_var@pos " << pos << " " << leftPos+it->first;
			BOOST_FOREACH(AlignedVariant av, it->second) {
				cout << " " << av.getString();
			}
			cout << endl;
		}
	}


	return true;
}

bool DetInDel::getHaplotypes(vector<Haplotype> &haps, const vector<Read> & reads,uint32_t pos, uint32_t & leftPos, uint32_t & rightPos, const AlignedCandidates & candidateVariants)
{


	uint32_t rs=(int(leftPos)>params.minReadOverlap)?(leftPos-params.minReadOverlap):0;
	uint32_t re=rightPos+params.minReadOverlap;
	string refSeq=getRefSeq(rs+1, re+1);

	HaplotypeDistribution hd(pos, refSeq, rs);

	// infer empirical haplotype distribution from WS alignment
	for (size_t r=0;r<reads.size();r++) {
		hd.insertRead(reads[r].getBam());
	}
	hd.setFrequencies();

	haps.clear();
	vector<Variant> indelVariants;

	/*
	if (!(candidateVariants.size()>0 && params.checkAllCIGARs==0)) {
		indelVariants=hd.getIndelVariantsAtMidPos();
	}


	// add any prespecified variants
	//indelVariants.insert(indelVariants.begin(), variants.begin(), variants.end());
	*/
	if (!params.quiet) {
		cout << "candidate_var@pos: " << pos ;
		BOOST_FOREACH(AlignedVariant v, candidateVariants.variants) {
			cout << " " << v.getStartHap() << "," << v.getString();
		}
		cout << endl;
	}


	// get haplotypes from empirical distribution

	try {
		HDIterator2 hdi(hd, params.maxHap, pos, leftPos, rightPos, params.noIndelWindow);

		double logNumHaps=hdi.getLogNumHaps();
		if (logNumHaps>log(params.skipMaxHap)) {
			cerr << "tid: " << params.tid << " pos: " << pos << " too many haplotypes [>exp(" << logNumHaps << ")]" << endl;
			return true;
		}

		//hdi.generateHapsWithIndels(haps, indels);
		vector<Haplotype> tmp_haps;
		hdi.generateHapsWithAlignedVariants(haps, candidateVariants, 0, params.changeINStoN);


	

		if (haps.size()>params.skipMaxHap || haps.size()*reads.size()>params.maxHapReadProd) {
			cerr << "tid: " << params.tid << " pos: " << pos << " too many haplotypes [>" << haps.size() << "]" << " numreads: " << reads.size() << endl;
			return true;
		}

		if (params.showHapDist) {
			cout << endl << "Empirical distribution: " << endl;
			cout << hdi << endl;
		}

		leftPos=hdi.start();
		rightPos=hdi.end();


		

		map<int, std::set<AlignedVariant> > variants;
		alignHaplotypes(haps,pos, leftPos, rightPos, variants);

		// remove duplicate reference-haplotypes of different length
		bool foundRef = false;
		for (size_t th=0;th<haps.size();th++) {
			const Haplotype & hap=haps[th];
			int num_indels =  hap.countIndels();
			int num_snps = hap.countSNPs();
			if (num_indels == 0 && num_snps == 0) {
				

				if (!foundRef) {
					tmp_haps.push_back(Haplotype(haps[th]));
					foundRef = true;
				}
			} else {
				tmp_haps.push_back(Haplotype(haps[th]));
			}
		}
		/*
		if (params.showCandHap) {
			for (size_t i=0;i<haps.size();i++) {
				cout << "PRE FILTER hdi[" << i << "]:" << haps[i] << endl;
			}
		}
		*/

		typedef map<int, AlignedVariant>::const_iterator It;
		haps.swap(tmp_haps);

		int nh=0;
		if (params.showCandHap) {
			for (size_t i=0;i<haps.size();i++) {
				cout << "POSTFILTER hdi[" << nh++ << "]:" << haps[i] << endl;
			}
		}
	}
	catch (string s) {
		if (s=="Blocks are not consecutive.") {
			cerr << "tid: " << params.tid <<  "pos: " << pos << s << endl;
			//return true;
			throw string("hapblock");
		} else {
			throw string(s);
		}
	}
	return false;
}




/*
// HaplotypeDistribution method
void DetInDel::computeLikelihoods(const vector<Haplotype> &haps, const vector<Read> & reads, vector<vector<MLAlignment> > & liks, uint32_t leftPos, uint32_t rightPos, vector<int> & onHap)
{
//	cout << "Computing likelihoods for all reads and haplotypes.\n";
	map<size_t, vector<size_t> > hapSizeToIdx;
	onHap = vector<int>(reads.size(),0); // records whether a read was aligned onto at least one haplotype

	typedef map<size_t, vector<size_t> >::const_iterator hapsCIt;

	for (size_t x=0;x<haps.size();x++) {
		hapSizeToIdx[haps[x].size()].push_back(x);
	}

	liks=vector<vector<MLAlignment> >(haps.size(),vector<MLAlignment>(reads.size()));

	for (hapsCIt it=hapSizeToIdx.begin();it!=hapSizeToIdx.end();it++) {
		const vector<size_t> hapIdx=it->second;
		// setup observation models for all the reads
		vector<ObservationModelFBMaxErr> oms; oms.reserve(reads.size());
		for (size_t r=0;r<reads.size();r++) {
			const Haplotype & hap=haps[hapIdx[0]];
			oms.push_back(ObservationModelFBMaxErr(hap, reads[r], leftPos, params.obsParams));
		}

		for (size_t h=0;h<hapIdx.size();h++) {
			size_t hidx=hapIdx[h];
			const Haplotype & hap=haps[hidx];

//			cout << "Haplotype[" << hidx << "]: " << endl;


			for (size_t r=0;r<reads.size();r++) {
				oms[r].changeHaplotype(hap);
				liks[hidx][r]=oms[r].calcLikelihood();
				if (!liks[hidx][r].offHapHMQ) onHap[r]=1;
//				cout << "[" << r << ": " << liks[hidx][r].ll << "] ";
				if (liks[hidx][r].ll>0.1) {
				  ObservationModelFBMaxErr om(hap, reads[r], leftPos, params.obsParams);
				 liks[hidx][r]=om.calcLikelihood();
				 cout << string(25,' ') << hap.seq << endl;
				  om.printAlignment(25);
				  cout << "h: " << h << " r: " << r << endl;
				  cout << bam1_qname(reads[r].getBam()) << endl;
				  cerr << "Likelihood>0" << endl;
				  exit(1);
				}
				if (isnan(liks[hidx][r]) || isinf(liks[hidx][r])) {
					cout << "NAN/Inf error" << endl;
					throw string("Nan detected");
				}
			}
	//		cout << endl;
		}
	}
}
*/
void DetInDel::computeLikelihoods(const vector<Haplotype> &haps, const vector<Read> & reads, vector<vector<MLAlignment> > & liks, uint32_t leftPos, uint32_t rightPos, vector<int> & onHap)
{
//	cout << "Computing likelihoods for all reads and haplotypes.\n";
	onHap = vector<int>(reads.size(),0); // records whether a read was aligned onto at least one haplotype

	typedef map<size_t, vector<size_t> >::const_iterator hapsCIt;

	liks=vector<vector<MLAlignment> >(haps.size(),vector<MLAlignment>(reads.size()));
	for (size_t hidx=0;hidx<haps.size();hidx++) {
		const Haplotype & hap=haps[hidx];
		for (size_t r=0;r<reads.size();r++) {
			ObservationModelFBMaxErr oms(hap, reads[r], leftPos, params.obsParams);
			liks[hidx][r]=oms.calcLikelihood();
			if (!liks[hidx][r].offHapHMQ) onHap[r]=1;
//				cout << "[" << r << ": " << liks[hidx][r].ll << "] ";
			if (liks[hidx][r].ll>0.1) {
				ObservationModelFBMaxErr om(hap, reads[r], leftPos, params.obsParams);
				liks[hidx][r]=om.calcLikelihood();
				cout << string(25,' ') << hap.seq << endl;
				om.printAlignment(25);
				cout << "hidx: " << hidx << " r: " << r << endl;
				cout << bam1_qname(reads[r].getBam()) << endl;
				cerr << "Likelihood>0" << endl;
				exit(1);
			}
			if (isnan(liks[hidx][r]) || isinf(liks[hidx][r])) {
				cout << "NAN/Inf error" << endl;
				throw string("Nan detected");
			}
		}
//		cout << endl;
	}
}


void DetInDel::computeHapPosition(const Haplotype & hap, const Read & read, vector<int> & alPos, int leftPos)
{
	// get position on haplotype of read alignment to reference from the aligned position of first and last base in the read

	const bam1_t *b=read.getBam();
	const bam1_core_t *c=&b->core;
	uint32_t* cigar=bam1_cigar(b);
	int k, end, start;
	end = c->pos;

	int offs=0, l=0, lend; // offset due to SOFT_SKIP at beginning
							// lend is base for which end is computed (there might be a SOFT_CLIP at the end of the read)

	bool al=false;
	for (k = 0; k < (int) c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;

		if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CREF_SKIP) al=true;

		if (!al && op == BAM_CSOFT_CLIP) offs += cigar[k] >> BAM_CIGAR_SHIFT;

		if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP)
					l += cigar[k] >> BAM_CIGAR_SHIFT;

		if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP) {
			end += cigar[k] >> BAM_CIGAR_SHIFT;
			lend=l;
		}
	}

	start=c->pos-leftPos; // make relative to alignment of haplotypes to reference
	end-=leftPos;

	// lookup start and end in alignment of haplotype to reference

	for (int x=0;x<int(hap.ml.hpos.size());x++) if (hap.ml.hpos[x]==start) {
		alPos.push_back(hap.ml.hpos[x]-offs);
		break;
	}

	for (int x=int(hap.ml.hpos.size())-1;x>=0;x--) if (hap.ml.hpos[x]==end) {
		alPos.push_back(hap.ml.hpos[x]-lend);
		break;
	}


}

void DetInDel::computeLikelihoodsFaster(const vector<Haplotype> &haps, const vector<Read> & reads, vector<vector<MLAlignment> > & liks, uint32_t leftPos, uint32_t rightPos, vector<int> & onHap)
{
	liks.clear();
	liks=vector<vector<MLAlignment> >(haps.size(),vector<MLAlignment>(reads.size()));
	onHap = vector<int>(reads.size(),0); // records whether a read was aligned onto at least one haplotype

	const unsigned int kmer=4;

	for (size_t h=0;h<haps.size();h++) {
		const Haplotype & hap = haps[h];
		//cout << "Haplotype[" << h << "]: " << endl;
		HapHash hash(kmer, hap);
		for (size_t r=0;r<reads.size();r++) {
			// given BWA alignment of read to reference, estimate a number of likely alignments to the haplotype
			vector<int> alPos;
			computeHapPosition(hap, reads[r], alPos, leftPos);
		//	cout << "[" << r << ": " << alPos.size() ;


			ObservationModelS om(hap, reads[r], leftPos, params.obsParams);

			// align using guessed alignments and the haplotype hash
			liks[h][r]=om.align(hash);
		//	cout << "," << liks[h][r].ll << "] ";
			if (!liks[h][r].offHapHMQ) onHap[r]=1; // if on-haplotype with artificial high mapping quality

			/*
			seqan::Score<int> score(-1, -460, -100,-960);
			ObservationModelSeqAn om2(hap, reads[r], leftPos, params.obsParams, score);
			om2.align();
			*/

 		}
		//cout << "done" << endl;
	}

	// check HMQ off-haplotype reads

	// realign a couple of high-mapping quality reads to obtain new candidate indels
	// propose new set of candidate haplotypes by realigning all reads to the new set of candidate haplotypes


	// --bamFiles /Users/caa/Documents/workspace/DInDelFastProb/bamfiles.txt --output test --region 12036338-12036340 --maxReadLength 60 --tid 17 --maxHap 8 --showEmpirical  --minReadOverlap 20 --width 60 --maxLengthIndel 10 --ref /Users/caa/data/human_b36_female.Y.fa --pError 0.0005  --maxRead 2000 --computeMAP
}

double DetInDel::getPairPrior(const AlignedVariant & av1, const AlignedVariant & av2, int leftPos, const AlignedCandidates & candidateVariants)
{
	std::set<AlignedVariant> vars;  vars.insert(av1); vars.insert(av2);
	double ll = 0.0;
	BOOST_FOREACH(AlignedVariant avar, vars) {
		double lnf = 0.0;
		int type = avar.getType();
		const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());

		if (type==Variant::SNP) lnf = log(params.priorSNP); else if (type==Variant::DEL || type==Variant::INS) lnf = log(params.priorIndel);
		if (av==NULL) {
			ll += lnf;
		} else {
			double prior = av->getFreq();
			if (prior<0.0) ll += lnf; else ll+=log(prior);
		}
	}

	return ll;

}

double DetInDel::getHaplotypePrior(const Haplotype & h1, const Haplotype & h2, int leftPos, const AlignedCandidates & candidateVariants)
{
	// returns log prior probability of the pair of haplotypes, based on settings in params
	// one day maybe change prior for known SNPs
	double ll=0.0;

	// count indels
	typedef  map <int, AlignedVariant>::const_iterator AVIt;
	//hash_map <int, int> indels, snps;
	std::set <AlignedVariant> indels, snps;
	for (AVIt it=h1.indels.begin();it!=h1.indels.end();it++) if (it->second.getString().find("*REF")==string::npos && it->second.getString().find("=>")==string::npos ) {
			//indels[it->first]=1;
			indels.insert(it->second);
	}
	for (AVIt it=h2.indels.begin();it!=h2.indels.end();it++) if (it->second.getString().find("*REF")==string::npos && it->second.getString().find("=>")==string::npos ) {
		//indels[it->first]=1;
		indels.insert(it->second);
	}

	for (AVIt it=h1.snps.begin();it!=h1.snps.end();it++) if (it->second.getString().find("*REF")==string::npos && it->second.getString().find("=>D")==string::npos) {
		snps.insert(it->second);
		//snps[it->first]=1;
	}
	for (AVIt it=h2.snps.begin();it!=h2.snps.end();it++) if (it->second.getString().find("*REF")==string::npos && it->second.getString().find("=>D")==string::npos) {
		snps.insert(it->second);
		//snps[it->first]=1;
	}

	BOOST_FOREACH(AlignedVariant indel, indels) {
	//	cout << "indel: " << indel.getString() << " " << indel.getStartHap();
		const AlignedVariant *av = candidateVariants.findVariant(indel.getStartHap()+leftPos, indel.getType(), indel.getString());
		if (av==NULL) {
		//	cout << " not found. " << endl;
			ll += log(params.priorIndel);
		} else {
			double prior = av->getFreq();
		//	cout << " found: " << prior << " " << av->getStartHap() << " " << av->getFreq() << endl;
			if (prior<0.0) ll += log(params.priorIndel); else ll+=log(prior);
		}
	}
	BOOST_FOREACH(AlignedVariant indel, snps) {
	//	cout << "snp: " << indel.getString() << " " << indel.getStartHap();
		const AlignedVariant *av = candidateVariants.findVariant(indel.getStartHap()+leftPos, indel.getType(), indel.getString());
		if (av==NULL) {
		//	cout << " not found. " << endl;
			ll += log(params.priorIndel);
		} else {
			double prior = av->getFreq();
		//	cout << " found: " << prior << " " << av->getStartHap() << " " << av->getFreq() << endl;
			if (prior<0.0) ll += log(params.priorIndel); else ll+=log(prior);
		}
	}
	/*
	BOOST_FOREACH(AlignedVariant snp, snps) {
		const AlignedVariant *av = candidateVariants.findVariant(snp.getStartHap()+leftPos, snp.getString());
		if (av==NULL) ll += log(params.priorSNP); else {
			double prior = av->getFreq();
			if (prior<0.0) ll += log(params.priorSNP); else ll+=log(prior);
		}
	}
	*/
	/*
	int numIndels=int(indels.size());
	int numSNPs=int(snps.size());

	ll+=double(numIndels)*log(params.priorIndel);
	ll+=double(numSNPs)*log(params.priorSNP);
	*/
//	cout << "ll: " << ll << endl;
	return ll;
}

#define FILTERHAPS
#ifdef FILTERHAPS

void DetInDel::filterHaplotypes(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks,  vector<int> & filtered, map<pair<int, AlignedVariant>, VariantCoverage> & varCoverage, bool doFilter)
{

	const int debugfh = 0;
	int numFiltered = 0;
	int numHaps = int(haps.size());
	filtered = vector<int>(haps.size(),0);

	varCoverage.clear();


	typedef pair<int, AlignedVariant> PAV;
	map<PAV, vector< std::set<int> > > hVarCoverage; // coverage of each variant per haplotype, so that eventually only coverage of reads that were not filtered out are reported

	for (int h=0;h<int(haps.size());h++) {
		// check all reads aligned to this haplotype and select the ones that are not off-haplotype and aligned without indels and at most two high-quality mismatches
		std::set<int> selReads;
		for (size_t r=0;r<reads.size();r++) {
			int sel = 0;
			if (!liks[h][r].offHapHMQ && liks[h][r].numIndels == 0) { // && liks[h][r].numMismatch<=3) {
				selReads.insert(int(r));
				sel = 1;
			}
			if (debugfh) cout << "sel: " << "h: " << h << " " << bam1_qname(reads[r].getBam()) << " mpos: " << reads[r].matePos << " selected: " << sel << endl;

		}

		// check each variant in the haplotype

		bool allCovered = true; // all variants in haplotype should be covered by at least one read.
		for (map<int, AlignedVariant>::const_iterator it = haps[h].indels.begin();it!= haps[h].indels.end();it++) {
				const AlignedVariant & av = it->second;

				PAV pav(it->first, av);

				map<PAV,vector<std::set<int> > >::iterator pit  = hVarCoverage.find(pav);
				if (pit == hVarCoverage.end()) {
					hVarCoverage[pav] = vector< std::set<int> >(haps.size()*2);
				}

				if (av.getType() == Variant::INS || av.getType() == Variant::DEL) {
					int left = av.getLeftFlankRead() - params.obsParams.padCover;     // readFlankLeft is the first left unique base flanking the indel
					int right = av.getRightFlankRead() + params.obsParams.padCover;
					int leftV = av.getLeftFlankRead();
					int rightV = av.getRightFlankRead();


					int len = right-left+1;
					bool covered = false;
					int numdelcovered = 0;
					//cout << "left: " << left << " right: " << right << " len: " << len << endl;
					if (av.getType() == Variant::DEL) {
						// see if there is at least one read spanning the interval with at most one mismatch
						BOOST_FOREACH(int r, selReads) {
							std::set <int> c;
							int strand = 0;
							if (reads[r].isUnmapped()) {
								if (!reads[r].mateIsReverse()) strand = 1;
							} else {
								if (reads[r].isReverse()) strand = 1;
							}

							int nmm = 0;
							for (int b=0;b<=int(liks[h][r].hpos.size());b++) {
								int hb = liks[h][r].hpos[b];
								if (hb>=left && hb<=right) {
									c.insert(hb);
									if ( haps[h].seq[hb]!='N' && haps[h].seq[hb]!=reads[r].seq.seq[b]) nmm++;
								}

							}
							int cov = 0;
							if (int(c.size())>= len && nmm<=params.obsParams.maxMismatch) {
								cov = 1;
								numdelcovered++;
								hVarCoverage[pav][h+strand*numHaps].insert(r);
							}
							//cout << "RC" << bam1_qname(reads[r].getBam()) << " cov: " << cov << endl;
						}

						if (numdelcovered>=1) {
							covered = true;
						}
					} else if (av.getType() == Variant::INS) {
						// see if all bases in the haplotype from left to right are covered by at least one read that matches the haplotype without indels and at most 2 mismatches
						vector<int> hapBaseCovered(len, 0);
						vector<int> thisReadCovered(len, 0);
						int lenins = av.getSeq().size();

						BOOST_FOREACH(int r, selReads) {
							for (int x=0;x<len;x++) thisReadCovered[x]=0;
							int nmm = 0;
							std::set <int> c;
							
							// determine read strand
							int strand = 0;
							if (reads[r].isUnmapped()) {
								if (!reads[r].mateIsReverse()) strand = 1;
							} else {
								if (reads[r].isReverse()) strand = 1;
							}

							for (int b=0;b<=int(liks[h][r].hpos.size());b++) {
								int hb = liks[h][r].hpos[b];
								if (hb>=left && hb<=right) {
									// covered even if there is a mismatch
									thisReadCovered[hb-left]+=1;
									c.insert(hb);
									// count number of mismatches
									if (haps[h].seq[hb]!=reads[r].seq.seq[b]) nmm++;
								}
							}

						        bool thisread_covered = false;
							// for insertion <= 10 bp, whole insertion+padCover must be covered with at most one error by at least one read (ie not just covered by multiple reads together)
							if ( (lenins>10 && nmm<=params.obsParams.maxMismatch) || (lenins<=10 && int(c.size())>=len && nmm<=params.obsParams.maxMismatch)) {
								thisread_covered = true;
								for (size_t x=0;x<thisReadCovered.size();x++) {
									hapBaseCovered[x] += thisReadCovered[x];
									if (thisReadCovered[x]==0) {
										thisread_covered = false;
									}
									if (debugfh) cout << " " << hapBaseCovered[x];
								}
								if (thisread_covered) {
									hVarCoverage[pav][h+strand*numHaps].insert(r);
								}
							}
							

							if (0) cout << " hap " << h << " var: " << av.getString() << " len: " << len << " " << bam1_qname(reads[r].getBam()) << " nmm: " << nmm << " c.size(): " << c.size() << " mpos: " << reads[r].matePos << " covered: " << thisread_covered << endl;
							if (thisread_covered) covered=true;


						}
					}

					if (!covered) {
						allCovered = false;
						break;
					}
					if (debugfh) cout << "hap" << h << " var: " << av.getString() << " COVERED:" << covered << endl;

				}

		}
		if (doFilter) {
			if (!allCovered) {
				numFiltered++;
				filtered[h]=1;
			}
		}
		if (debugfh) cout << "Haplotype[" << h << "]: filtered " << filtered[h] << endl;

	}
	cout << "Filtered " << numFiltered << " haplotypes." << endl;
	// determine coverage of each variant
	for (map<PAV, vector <std::set<int> > >::const_iterator it = hVarCoverage.begin();it != hVarCoverage.end(); it++) {
		const PAV & pav = it->first;
		std::set<int> rf, rr; // forward and reverse strand reads
		for (int h=0;h<numHaps;h++) if (filtered[h]!=1) {
			rf.insert(hVarCoverage[pav][h].begin(), hVarCoverage[pav][h].end());
			rr.insert(hVarCoverage[pav][h+numHaps].begin(), hVarCoverage[pav][h+numHaps].end());
		}
		varCoverage[pav]=VariantCoverage(int(rf.size()), int(rr.size()));
	}


}
#endif

void DetInDel::estimateHaplotypeFrequenciesBayesEM(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<double> & hapFreqs, vector <HapEstResult > & posteriors,  uint32_t candPos, uint32_t leftPos,   uint32_t rightPos, OutputData & glfData, int index, const AlignedCandidates & candidateVariants, string program="all")

{

	// estimate haplotype frequencies using EM
	hapFreqs.clear();

	size_t nh=haps.size();
	size_t nr=reads.size();

	vector<double> rl(nh*nr,0.0); // read given haplotype likelihoods

	vector<double> z(nh*nr,0.0); // expectations of read-haplotype indicator variables
	vector<double> pi(nh); // haplotype frequencies
	vector<double> lpi(nh); // expectation of log frequencies
	vector<double> nk(nh,0.0), ak(nh,0.0); // counts for each haplotype

	hapFreqs=nk;

	int numUnmappedRealigned=0;
	int idx=0;
	int numReadOffAllHaps=0;
	for (size_t r=0;r<nr;r++) {
		int offallhap=1;
		for (size_t h=0;h<nh;h++) {
			// initialize read-haplotype likelihoods
			rl[idx]=liks[h][r].ll;
			if (!liks[h][r].offHap) offallhap=0;
			idx++;
		}
		if (offallhap) {
			numReadOffAllHaps++;

		}else {
			if (reads[r].isUnmapped()) numUnmappedRealigned++;
		}
	}


	// filter reads

	vector<int> filtered(nh, 0);
	map<pair<int, AlignedVariant>, VariantCoverage> varCoverage;

	//if (params.filterHaplotypes) {
	cout << "ALWAYS CALLING ::filterHaplotypes" << endl;
	filterHaplotypes(haps, reads,liks, filtered, varCoverage, params.filterHaplotypes);
	//}


	typedef pair<int, AlignedVariant> PAV;

	std::set< PAV > allVariants;
	map<int, std::set<PAV> > allVariantsByPos; //

	typedef map<int, AlignedVariant>::const_iterator It;
	typedef map<int, std::set<PAV> >::const_iterator PIt;

	for (size_t th=0;th<nh;th++) {
		const Haplotype & hap=haps[th];
		for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
			if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
				allVariants.insert(PAV(it->first,it->second));
				allVariantsByPos[it->first].insert(PAV(it->first,it->second));
			}
		}
	}

	// set active variants, and divide into snps and indels
	vector< std::set< PAV > >  activeVariants, activeSNPs, activeIndels;

	int nav=0;
	int PRID=-1;
	if (program=="all") {
		std::set<PAV> snps, indels;
		BOOST_FOREACH(PAV pav, allVariants) {
			if (pav.second.isSNP()) {
				snps.insert(pav);
			} else if (pav.second.isIndel()) {
				indels.insert(pav);
			}
		}

		// both (double prior)
		activeVariants.push_back(allVariants);
		activeIndels.push_back(indels);
		activeSNPs.push_back(snps);
		nav++;
		PRID=1;
	} else if (program=="singlevariant") {
		std::set < std::set<PAV> > ssPAV;
		for (size_t h=0;h<haps.size();h++) if (filtered[h]==0) {
			const Haplotype & hap=haps[h];


			//cout << "hap[" << h << "].seq: " << hap.seq << endl;

			//cout << "vars:";
			std::set<PAV> act;
			for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
				if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
					act.insert(PAV(it->first, it->second));
			//		cout << "[" << it->first << "," << it->second.getString() << "]";
				}
			}
			//cout << endl;
			ssPAV.insert(act);
		}

		nav=0;
		BOOST_FOREACH(std::set<PAV> s, ssPAV) {
			std::set<PAV> snps, indels;
			BOOST_FOREACH(PAV pav, s) {
				if (pav.second.isSNP()) {
					snps.insert(pav);
				} else if (pav.second.isIndel()) {
					indels.insert(pav);
				}
			}

			// both (double prior)
			activeVariants.push_back(s);
			activeIndels.push_back(indels);
			activeSNPs.push_back(snps);
			nav++;
		}

		PRID=2;
	} else if (program == "priorpersite") {
		nav = 0;

		// add reference haplotype
		activeVariants.push_back(std::set<PAV>());
		activeIndels.push_back(std::set<PAV>());
		activeSNPs.push_back(std::set<PAV>());



		for (map<int, std::set<PAV> >::const_iterator site_it = allVariantsByPos.begin();site_it!=allVariantsByPos.end();site_it++) {
			std::set<PAV> snps, indels;
			BOOST_FOREACH(PAV pav, site_it->second) {
				if (pav.second.isSNP()) snps.insert(pav);
				else if (pav.second.isIndel()) indels.insert(pav);
			}

			int maxStateSnp = (snps.size()==0)?1:2;
			int maxStateIndel = (indels.size()==0)?1:2;

			int prevNumActive = activeVariants.size();
			int sSnp = 1, sIndel = 1;
			//for (int sSnp = 0;sSnp<maxStateSnp;sSnp++) {
			//	for (int sIndel = 0; sIndel < maxStateIndel; sIndel++) {
			//
			//		if (sSnp == 1 || sIndel == 1) {
						// extend previous activeVariants

						for (int pna = 0;pna<prevNumActive;pna++) {
							std::set<PAV> av = activeVariants[pna];
							std::set<PAV> aIndels = activeIndels[pna];
							std::set<PAV> aSNPs = activeSNPs[pna];
							if (sSnp == 1) {
								av.insert(snps.begin(),snps.end());
								aSNPs.insert(snps.begin(),snps.end());
							}
							if (sIndel == 1) {
								av.insert(indels.begin(),indels.end());
								aIndels.insert(indels.begin(),indels.end());
							}

							activeVariants.push_back(av);
							activeIndels.push_back(aIndels);
							activeSNPs.push_back(aSNPs);
						}

			//		}
		//		}
		//	}

		}
		nav = activeVariants.size();

		PRID = 3;

	} else {
		cerr << "Unknown EM option" << endl;
		exit(1);
	}

	vector<int> compatible(nh,0);
	vector<double> logliks(nav,0.0);
	vector<double> logpriors(nav, 0.0);
	vector<double> post(nav,0.0);
	vector<double> freqs(nav*nh,0.0);

	// create matrix of which variant is active in which set
	idx=0;
	int nv=int(allVariants.size());
	vector<int> active(nav*nv,0), hapHasVar(nh*nv,0);

	//cout << "active: " << endl;
	BOOST_FOREACH(PAV pav, allVariants) {
	//	 cout << pav.first << " " << pav.second.getString() << " ";
		for (int a=0;a<nav;a++) {
			if (activeVariants[a].find(pav)!=activeVariants[a].end()) active[a*nv+idx]=1;
	//		cout << " " << active[a*nv+idx];
		}
	//	 cout << endl;
		for (size_t h=0;h<nh;h++) {
			It it=haps[h].indels.find(pav.first);
			if (it!=haps[h].indels.end() && it->second.getString()==pav.second.getString()) hapHasVar[h*nv+idx]=1;
//			cout << "[" <<  active[h*nv+idx] << " " << hapHasVar[h*nv+idx]<< " ]";

		}
//		cout << endl;
		idx++;
	}

	/*
	cout << "allVariants: ";
	BOOST_FOREACH(PAV pav, allVariants) {
		cout << " [" << pav.first << " " << pav.second.getString() << "]";
	}
	cout << endl;
	*/





	double logz=-HUGE_VAL;

	double a0=params.bayesa0;

	for (int th=0;th<nav;th++) {

		// set active variants

		double logprior=0.0;

		map <int, int> sites;
		BOOST_FOREACH(PAV pav, activeSNPs[th]) {
			sites[pav.first]=1;

			const AlignedVariant & avar = pav.second;
			const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());
			int type = pav.second.getType();

			double lnf = 0.0;

			if (type==Variant::SNP) lnf = log(params.priorSNP); else if (type==Variant::DEL || type==Variant::INS) lnf = log(params.priorIndel);
			if (av==NULL) {
				logprior += lnf;
			} else {
				double prior = av->getFreq();
				if (prior<0.0) logprior += lnf; else logprior+=log(prior);
			}

		}
		BOOST_FOREACH(PAV pav, activeIndels[th]) {

			const AlignedVariant & avar = pav.second;
			const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());
			int type = pav.second.getType();

			double lnf = 0.0;

			if (type==Variant::SNP) lnf = log(params.priorSNP); else if (type==Variant::DEL || type==Variant::INS) lnf = log(params.priorIndel);
			if (av==NULL) {
				logprior += lnf;
			} else {
				double prior = av->getFreq();
				if (prior<0.0) logprior += lnf; else logprior+=log(prior);
			}

			sites[pav.first]=2;
		}


		/*
		for (map<int,int>::const_iterator it=sites.begin();it!=sites.end();it++) {
			if (it->second==2) logprior+=log(params.priorIndel); else if (it->second==1) logprior+=log(params.priorSNP);
		}
		*/

		logpriors[th]=logprior;

//		cout << "Number of indels: " << ni << " number of SNPs: " << ns << endl;

		// check haplotypes

		int numah=0; // number of haplotypes for which frequencies will be estimated
		for (size_t h=0;h<nh;h++) {
			compatible[h]=1;
			if (filtered[h]!=0) {
				compatible[h]=0;
			} else {
				for (It it=haps[h].indels.begin();it!=haps[h].indels.end();it++) {
					if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D') && activeVariants[th].find(PAV(it->first,it->second))==activeVariants[th].end()) {
						// haplotype h has a non-ref variant that is not one of the active variants
						compatible[h]=0;
						break;
					}
				}
			}
			if (compatible[h]) numah++;
		}


		// run EM for this set of active variants
		bool converged=false;
		double tol=params.EMtol;

		double eOld=-HUGE_VAL, eNew;

		// initialize frequencies
		for (size_t h=0;h<nh;h++) if (compatible[h]) lpi[h]=log(1.0/double(numah)); else lpi[h]=-100;

		/*
		for (size_t r=0;r<nr;r++) {
				cout << "rl[" << r << "]:";
				for (size_t h=0;h<nh;h++) {
					cout << " " << rl[r*nh+h];
				}
				cout << endl;
			}
		*/
		double loglik, llNew, llOld=-HUGE_VAL;
		int iter=0;
		while (!converged) {
			//cout << endl << "EM[" << iter << "]:" << endl;
			// compute expectation of indicator variables
			for (size_t h=0;h<nh;h++) nk[h]=0.0;

			loglik=0.0;

			int idx=0;
			for (size_t r=0;r<nr;r++) {
				double lognorm=-HUGE_VAL;
				// compute responsibilities
				for (size_t h=0;h<nh;h++) {
					z[idx]=lpi[h]+(rl[idx]);
					lognorm=addLogs(lognorm, z[idx]);
					idx++;
				}
				idx-=nh;
				// normalize and compute counts
				for (size_t h=0;h<nh;h++) {
					z[idx]-=lognorm;
					z[idx]=exp(z[idx]);
					nk[h]+=z[idx];
					idx++;
				}
				loglik+=lognorm;

			}
			// compute frequencies
			//cout << "pi: ";

			double ahat=0.0;
			for (size_t h=0;h<nh;h++) if (compatible[h]) {
					ak[h]=nk[h]+a0; // a0 is Dirichlet prior parameter
					ahat+=ak[h];
			}
			double dahat=boost::math::digamma(ahat);

			for (size_t h=0;h<nh;h++) {

				// do variational bayes
				if (compatible[h]) {
					lpi[h]=boost::math::digamma(ak[h])-dahat;
					pi[h]=log((a0+nk[h])/(double(numah)*a0+double(nr)));
				} else {
					lpi[h]=-100;
					pi[h]=-100;
				}
			//	cout << " " << pi[h];
			//	zp+=exp(pi[h]);
			}
			//cout << " zp: " << zp << endl;


			idx=0;
			eNew=0.0;
			for (size_t r=0;r<nr;r++) {
				for (size_t h=0;h<nh;h++) {
				// compute responsibilities
					eNew+=z[idx]*( pi[h]+rl[idx]);
					idx++;
				}
			}
			//cout << "eOld: " << eOld << " " << " eNew: " << eNew; for (size_t h=0;h<nh;h++) { cout << " " << pi[h]; } cout << endl;
			/*
			for (size_t r=0;r<nr;r++) {
				cout << "z[" << r << "]:";
				for (size_t h=0;h<nh;h++) {
					cout << " " << z[r*nh+h];
				}
				cout << endl;
			}
			*/
			llNew=loglik;
			//cout << "loglik: " << loglik << endl;
			if (0 && llOld>llNew+1e-5)  {
				cerr << "ERROR: nr: " << nr << " eOld: " << eOld << " eNew: " << eNew << " diff: " << eOld-eNew << endl;
				cout << "ERROR: nr: " << nr << " eOld: " << eOld << " eNew: " << eNew << " diff: " << eOld-eNew << endl;
				cerr << "ERROR: nr: " << nr << " llOld: " << llOld << " eNew: " << llNew << " diff: " << llOld-llNew << endl;
				cout << "ERROR: nr: " << nr << " llOld: " << llOld << " eNew: " << llNew << " diff: " << llOld-llNew << endl;

				//throw string("EM Error in estimateHapFreqs");
				//iter=100;

			}
			converged=(fabs(eOld-eNew))<tol || iter>25;
			//cout << "iter: " << iter << " eOld: " << eOld << " eNew: " << eNew << endl;

			eOld=eNew;
			llOld=llNew;
			iter++;


		}


		// check sum

		double z=0.0;
		for (size_t y=0;y<nh;y++) {
			z+=exp(pi[y]);
		}

		if (0) {
		cout << "th: " << th << endl;
		for (size_t y=0;y<nh;y++) {
			cout << "[" << y << "," << exp(pi[y]) << "]";
		}
		cout << endl << endl;
		}

		logliks[th]=loglik;
		logz=addLogs(logz, logliks[th]+logprior);
		for (size_t h=0;h<nh;h++) { freqs[th*nh+h]=exp(pi[h])/z; }
		//for (size_t h=0;h<nh;h++) { cout << " " << freqs[th*nh+h]; } cout << endl;

		//cout << "loglik: " << loglik << " " << logliks[th] << " logprior: " << logprior << endl << endl;

	}

	for (int a=0;a<nav;a++) {
		post[a]=exp(logliks[a]+logpriors[a]-logz);
		//cout << "post[" << a << "]: " << post[a] << endl;
	}

	for (size_t h=0;h<nh;h++) {
		hapFreqs[h]=0.0;
	}

	for (int th=0;th<nav;th++) for (size_t h=0;h<nh;h++) {

		hapFreqs[h]+=exp(logliks[th]+logpriors[th]-logz)*freqs[th*nh+h];
	}

	/*
	cout << "hapFreqs:\n ";
	for (size_t th=0;th<nh;th++) {
		cout << "hapFreq[" << th << "]: " << hapFreqs[th] << endl;
		cout << "H" << th << ": " << logliks[th] << " " << logpriors[th] << " ";
		for (size_t h=0;h<nh;h++) {
			cout << " " << freqs[th*nh+h];
		}
		cout << endl;
	}
	*/

	cout << endl;

	// compute marginal posteriors for the individual variants
	vector< std::set<int> > readidx(myBams.size());
	for (int r=0;r<int(nr);r++) readidx[reads[r].poolID].insert(r);

	vector<double> prior(nh*nh,0.0);

	idx=-1;
	BOOST_FOREACH(PAV pav, allVariants) {
		idx++;

		double logp=-HUGE_VAL;
		double freq=0.0;
		for (int th=0;th<nav;th++) {
			if (active[th*nv+idx]) {
				logp=addLogs(logp, logliks[th]+logpriors[th]);
			}
		}

		for (size_t h=0;h<nh;h++) {
			if (hapHasVar[h*nv+idx]) {
				freq+=hapFreqs[h];
			}
		}

		logp-=logz;

		// compute marginalized haplotype frequencies
		vector<double> prior(nh*nh,0.0);

		bool doGLF=false; //(candPos==leftPos+pav.first)?true:false;

		const AlignedVariant & avar = pav.second;
		const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());
		doGLF = (av==NULL)?false:true;

		//change this
if (params.outputGLF && doGLF) {
		map<int,int> otn; // old to new haplotype index
		// vector for mapping old haplotype index to new haplotype index
		vector<int> marsum(nv, 0);
		int s=1;
		for (int y=0;y<nv;y++) {
			if (y!=idx) {
				marsum[y]=s;
				s*=2;
			} else marsum[y]=0;
		//	cout << "marsum[" << y << "]: " << marsum[y] << endl;
		}


		// make a map of new marginalized states to the corresponding new index in the haplotype arrays.
		map<int, int> mar_states;
		for (int h=0;h<int(nh);h++) {
			int nidx=0;
			for (int v=0;v<nv;v++) nidx+=marsum[v]*hapHasVar[h*nv+v];
			map<int,int>::iterator it=mar_states.find(nidx);
			if (it!=mar_states.end()) otn[h]=it->second; else {
				int ns=int(mar_states.size());
				mar_states[nidx]=ns;
				otn[h]=ns;
			}
		//	cout << "oth["<<h<<"]: " << otn[h] << endl;
		}

		int nmarhap=mar_states.size();

		vector<double> marFreqs(nmarhap,0.0); //marginal frequencies
		for (size_t h=0;h<nh;h++) {
			int newh=otn[h];
			marFreqs[newh]+=hapFreqs[h];
		}
		// convert back to conditional frequencies/probabilities
		for (int h=0;h<nmarhap;h++) {
			if (marFreqs[h]<1e-16) marFreqs[h]=-50; else marFreqs[h]=log(marFreqs[h]);
		}
		//	cout << "marFreq[]: "; for (size_t x=0;x<nmarhap;x++) cout << " " << marFreqs[x]; cout << endl;

		// compute for every pair of haplotypes the likelihood of drawing it given priors on the variants and the estimated frequencies
		for (size_t h1=0;h1<nh;h1++) {
			for (size_t h2=h1;h2<nh;h2++) {
				prior[h1*nh+h2]=marFreqs[otn[h1]]+marFreqs[otn[h2]];
		//		cout << "prior: " << h1 << " " << h2 << ": " << prior[h1*nh+h2] << endl;
			}
		}
}
		// check which how many reads to a haplotype with this variant

		std::set<int> support;
		int totnf=0, totnr=0;
		double log5=log(0.5);
		if (1) for (size_t b=0;b<myBams.size();b++) {
			double msq=0.0;

			int nf=0;
			int nr=0;
			vector<double> lik(3,0.0);

			if (readidx[b].size()) {
				// compute RMS of mapping qualities

				int pos=pav.first;
				msq=0.0;
				int n=0;
				if (params.outputGLF && doGLF) {
				for (int i=0;i<3;i++) lik[i]=-HUGE_VAL;
				for (size_t h1=0;h1<nh;h1++) for (size_t h2=h1;h2<nh;h2++) {
					int genotype=hapHasVar[h1*nv+idx]+hapHasVar[h2*nv+idx];
					double pr=prior[h1*nh+h2];
				//	cout << "prior: " << h1 << " " << h2 << ": " << prior[h1*nh+h2];
					double ll=pr;
					BOOST_FOREACH(int r, readidx[b]) {
						ll+=log5+addLogs(rl[r*nh+h1],rl[r*nh+h2]);
					}

					lik[genotype]=addLogs(lik[genotype],ll);
				//	cout << " ll: " << ll << " " << lik[genotype] << endl;
				}
				}

				BOOST_FOREACH(int r, readidx[b]) {
					//cout << " " << 10*log10(1-reads[r].mapQual);

					size_t mlidx; double ml=-HUGE_VAL;
					std::set<size_t> mlis;
					for (size_t hi=0;hi<nh;hi++) {
						if (liks[hi][r].ll>=ml) {
							mlidx=hi;
							ml=liks[hi][r].ll;
						}
					}
					for (size_t hi=0;hi<nh;hi++) {
						if (liks[hi][r].ll>=ml-1e-7) {
							mlis.insert(hi);
						}
					}
					bool nrt=false, nft=false;

					map<int,bool>::const_iterator it;
					BOOST_FOREACH(size_t h, mlis) {
						bool covered=false;
						if (pav.second.isIndel()) {
							it=liks[h][r].hapIndelCovered.find(pos);
							if (it!=liks[h][r].hapIndelCovered.end() && it->second) covered=true;
						} else if (pav.second.isSNP()) {
							it=liks[h][r].hapSNPCovered.find(pos);
							if (it!=liks[h][r].hapSNPCovered.end() && it->second) covered=true;
						}
						if (covered) { // hapHasVar[] is to check whether haplotype truely has variant or only the reference
							if (hapHasVar[h*nv+idx]) {
								//if (pav.first+leftPos==43017596) {
								//	cout << "h: " << h << " r: " << r << " read: " << reads[r] << endl;
								//}


								if  (reads[r].onReverseStrand) nrt=true; else nft=true;
							}

						}
					}
					double mq=-10*log10(1.0-reads[r].mapQual);
					msq+=mq*mq;
					n++;
					if (nft) nf++;
					if (nrt) nr++;
				} // foreach read r
			//	cout << endl;
				if (n!=0) msq=sqrt(msq/double(n)); else msq=0.0;
				if (nf+nr>0) support.insert(b);
				totnf+=nf;
				totnr+=nr;
			}  // readidx[b].size()

			if (params.outputGLF && doGLF) {
				OutputData::Line line(glfData);
				line.set("msg", "ok");
				line.set("index", index);
				line.set("tid", params.tid);
				line.set("analysis_type",program);
				line.set("indidx",b);
				line.set("was_candidate_in_window",1);
				line.set("lpos",leftPos);
				line.set("rpos",rightPos);
				line.set("center_position",candPos);
				line.set("realigned_position",pav.first+leftPos);
				line.set("post_prob_variant", exp(logp));
				line.set("est_freq", freq);
				line.set("logZ", logz);
				line.set("nref_all", pav.second.getString());
				line.set("num_reads", readidx[b].size());
				line.set("msq",msq);
				line.set("num_cover_forward", nf);
				line.set("num_cover_reverse", nr);
				line.set("num_unmapped_realigned", numUnmappedRealigned);
				line.set("var_coverage_forward", varCoverage[pav].nf);
				line.set("var_coverage_reverse", varCoverage[pav].nr);

				if (b==0) {
					// output haplotypes and frequencies

					ostringstream os;
					bool ifh = true;
					for (size_t h=0;h<haps.size();h++) {
						if (hapFreqs[h]>1.0/double(2*reads.size())) {
							bool isfirst = true;
							if (!ifh) os << ";";
							ifh = false;
							int nvars = 0;
							for (map<int, AlignedVariant>::const_iterator it=haps[h].indels.begin();it!=haps[h].indels.end();it++) if (it->second.getString()!="*REF") {
								nvars++;
								if (!isfirst) os << ",";
								isfirst = false;
								os << leftPos+it->first << "," << it->second.getString();
							}
							if (nvars==0) os << "REF";
							os << ":" << hapFreqs[h];
						}
					}
					line.set("hapfreqs", os.str());

				}

				string likstring;

				for (int n=0;n<3;n++) {
					ostringstream o;
					o << lik[n];
					if (n==0) likstring.append("0/0:"); else if (n==1) likstring.append("0/1:"); else likstring.append("1/1:");
					likstring.append(o.str());
					if (n<2) likstring.append(";");
				}
				line.set("glf",likstring);
				glfData.output(line);
			}
			// glfOutput <<  PRID << " "  << b << " " << params.tid << " " << candPos << " " << pav.first+leftPos << " " << pav.second.getString() << " " << reads.size() <<  " " << msq << " "  << nf << " " << nr << " " << lik[0] << " " << lik[1] << " " << lik[2] << endl;


		} // foreach b
		posteriors.push_back(HapEstResult(pav.second, pav.first,exp(logp),freq, totnf, totnr));
	}


	cout << "candPos: " << candPos << " numReadOffAllHaps: " << numReadOffAllHaps << " logz: " << logz << endl;



    if (params.outputPooledLikelihoods) {

            // output which variants are active in which haplotype


            stringstream os; os << params.fileName << "." << params.tid << "." << candPos;
	    string oprefix = os.str();

	    string fname = oprefix;
	    fname.append(".hapvars");
            ofstream of(fname.c_str());
            if (!of.is_open()) {
                    throw string("Cannot open file ").append(fname).append(" for writing .hapvars file");
            }

            idx=0;
            BOOST_FOREACH(PAV pav, allVariants) {
                    stringstream o;
                    o << params.tid << " " << leftPos+pav.first << " " << pav.second.getString();
                    of << o.str() << string(50-o.str().length(),' ');
                    for (size_t h=0;h<nh;h++) {
                            of << " " << hapHasVar[h*nv+idx];
                    }

                    of << endl;
                    idx++;

            }

            of.close();


	    string prefix;
            stringstream os5; os5 << "EM " << params.tid << " " << candPos << " " << reads.size(); prefix.append(os5.str());

	    fname = oprefix;
	    fname.append(".hapfreqs");
	    of.open(fname.c_str());
	    outputHapsAndFreqs(&of,prefix,haps,hapFreqs, leftPos);
	    of.close();

            fname.clear();
            fname=oprefix;
	    oprefix.append(".liks");

	    cout << "fname: " << fname << endl;
            of.open(fname.c_str());
            if (!of.is_open()) {
                    throw string("Cannot open file ").append(fname).append(" for writing .liks file");
            }



            // output all likelihoods

            for (size_t r=0;r<nr;r++) {

                    of << r << " " << bam1_qname(reads[r].getBam()) << " " << log(1.0-reads[r].mapQual) << " " << reads[r].poolID;


                    for (size_t h=0;h<nh;h++) {
                            of << " " << liks[h][r].ll;
                    }


                    for (size_t h=0;h<nh;h++) {
                            of << " " << liks[h][r].offHap;
                    }
                    of << endl;
            }
	    of.close();

	    // output alignments
	    fname = oprefix;
	    fname.append(".alignments");
            cout << "fname: " << fname << endl;
            of.open(fname.c_str());
            if (!of.is_open()) {
                    throw string("Cannot open file ").append(fname).append(" for writing .liks file");
            }



	    for (size_t r=0;r<nr;r++) {
		cout << "###" << endl;
		cout << "read: " << bam1_qname(reads[r].getBam()) << " mpos: " << reads[r].matePos << endl;
		cout << "isUnmapped: " << reads[r].isUnmapped() << endl;
		// compute maximum alignment likelihood
		double lq = 0.0;
		for (size_t b=0;b<reads[r].size();b++) lq += log(reads[r].qual[b]);

		cout << "Max alignment loglik: " << lq << endl;

		double maxll = -HUGE_VAL;
		std::set <int> mlhaps;
		for (int h=nh-1;h>=0;h--) if (liks[h][r]>maxll) { maxll = liks[h][r]; }
		for (int h=nh-1;h>=0;h--) mlhaps.insert(h); //if (fabs(liks[h][r]-maxll)<0.01) mlhaps.insert(h);
		BOOST_FOREACH(int hidx, mlhaps) {
			cout << "r: " << r << " hidx: " << hidx << " maxll:" << maxll << endl;
			ObservationModelFBMaxErr obs(haps[hidx], reads[r], leftPos, params.obsParams);
			cout << string(50,' ') << haps[hidx].seq << endl;
			obs.printAlignment(50);
		}
	    }

            of.close();
    }
}

#ifndef DIPLOIDGLF
void DetInDel::diploidGLF(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<double> & hapFreqs, vector <HapEstResult > & posteriors,  uint32_t candPos, uint32_t leftPos, uint32_t rightPos,  OutputData & glfData, int index, const AlignedCandidates & candidateVariants, string program="all")

{
	size_t nh=haps.size();
	size_t nr=reads.size();

	vector<int> filtered(nh, 0);
	map<pair<int, AlignedVariant>, VariantCoverage> varCoverage;
	filterHaplotypes(haps, reads, liks, filtered, varCoverage, params.filterHaplotypes);

	vector<double> rl(nh*nr,0.0); // read given haplotype likelihoods

	// get read-hap likelihoods and get number of reads that are off all-haplotypes.

	int idx=0;
	int numReadOffAllHaps=0;
	for (size_t r=0;r<nr;r++) {
		int no=1;
		for (size_t h=0;h<nh;h++) {
			// initialize read-haplotype likelihoods
			rl[idx]=liks[h][r].ll;
			if (!liks[h][r].offHap) no=0;
			idx++;
		}
		if (no) {
			numReadOffAllHaps++;
		}
	}

	// get all variants

	typedef pair<int, AlignedVariant> PAV;
	const int VARSNP=1;
	const int VARINDEL=2;


	std::set< PAV > allVariants;
	map<int, std::set<PAV> > allVariantsByPos; //


	typedef map<int, AlignedVariant>::const_iterator It;
	typedef map<int, std::set<PAV> >::const_iterator PIt;

	vector<int> hap_num_indels(nh, 0);
	vector<int> hap_num_candidate_indels(nh, 0);
	vector<int> hap_num_snps(nh, 0);

	int ref_hap_idx = -1;
	for (size_t th=0;th<nh;th++) {
		const Haplotype & hap=haps[th];
		//cout << "hap[" << th << "]: ";
		hap_num_indels[th] = hap.countIndels();
		hap_num_snps[th] = hap.countSNPs();

		if (hap_num_indels[th] == 0 && hap_num_snps[th] == 0) {
			//if (ref_hap_idx != -1) cout << string("Already have ref-hap!") << " " << ref_hap_idx << endl;
			//if (ref_hap_idx!=-1) cout << "RH: " << haps[ref_hap_idx].seq << endl;
			//cout << "TH: " << haps[th].seq << endl;
			ref_hap_idx = th;
			//int h1 = th;
			//cout << "RRRR IN: " << hap_num_indels[h1] << " SNP: " << hap_num_snps[h1]; cout << " H1:"; for (It it=haps[h1].indels.begin();it!=haps[h1].indels.end();it++) cout << it->first << "," << it->second.getString() << ";" << endl;

		}
		if (hap_num_indels[th] != 0) {
			//check how many were candidates in the input file
			int nc = 0;
			for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
				const AlignedVariant & avar = it->second;
				const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());
				if (av!=NULL) nc+=1;
			}
			hap_num_candidate_indels[th] = nc;
		}
	


		for (It it=hap.indels.begin();it!=hap.indels.end();it++) {
			if (!it->second.isRef() && !(it->second.isSNP() && it->second.getString()[3]=='D')) {
				//cout << " " << it->first << "," << it->second.getString();
				allVariants.insert(PAV(it->first,it->second));
				allVariantsByPos[it->first].insert(PAV(it->first,it->second));
			}
		}
		//cout << endl;
	}

	idx=0;
	map<int,int> posToPosIdx;
	vector<int> varPositions;
	for (PIt pit=allVariantsByPos.begin();pit!=allVariantsByPos.end(); pit++) {
		posToPosIdx[pit->first]=idx;
		varPositions.push_back(pit->first);
		idx++;
	}

	int numVarPos = allVariantsByPos.size();
	int nv=int(allVariants.size());

	vector<int> hapVar(nh*numVarPos,0);
	vector<int> varType(nv+1);
	vector<PAV> variants(nv+1);
	idx=1;


	BOOST_FOREACH(PAV pav, allVariants) {
		int type=VARSNP;
		if (pav.second.isIndel()) type=VARINDEL;
		varType[idx]=type;
		int posIdx=posToPosIdx[pav.first];
		for (size_t h=0;h<nh;h++) {
			It it=haps[h].indels.find(pav.first);
			if (it!=haps[h].indels.end() && it->second.getString()==pav.second.getString()) hapVar[h*numVarPos+posIdx]=idx;
		}
		variants[idx]=pav;



		idx++;
	}



	/*
	cout << "allVariants: ";
	BOOST_FOREACH(PAV pav, allVariants) {
		cout << " [" << pav.first << " " << pav.second.getString() << "]";
	}
	cout << endl;
	*/

	// check which how many reads to a haplotype with this variant

	std::set<int> readidx;
	for (size_t r=0;r<nr;r++) readidx.insert(r);



	// compute marginal posteriors for the individual variants
	vector<double> prior(nh*nh,0.0), pairs_posterior(nh*nh, 0);
	// compute for every pair of haplotypes the likelihood of drawing it given priors on the variants and the estimated frequencies
	for (size_t h1=0;h1<nh;h1++) {
		for (size_t h2=h1;h2<nh;h2++) {
			prior[h1*nh+h2]=getHaplotypePrior(haps[h1],haps[h2], leftPos, candidateVariants);
		}
	}


	vector<int> max_indel_pair(2,-1);
	vector<int> max_noindel_pair(2,-1);

	double max_ll_indel = -HUGE_VAL;
	double max_ll_noindel = -HUGE_VAL;
	for (size_t h1=0;h1<nh;h1++) if (filtered[h1]==0) for (size_t h2=h1;h2<nh;h2++) if (filtered[h2]==0) {
		double ll=0.0;
		for (size_t r=0;r<reads.size();r++){
			ll+=log(0.5)+addLogs(rl[r*nh+h1],rl[r*nh+h2]);
		}
		// now we have the log likelihood, store posterior
		pairs_posterior[h1*nh+h2] = ll + prior[h1*nh+h2];

		// prinhaplotype pair
		/*		
		cout << "POST: " << h1 << " " << h2 << " " <<  pairs_posterior[h1*nh+h2] << " PR: " << prior[h1*nh+h2];
		
		cout << " IN: " << hap_num_indels[h1] << " SNP: " << hap_num_snps[h1]; cout << " H1:"; for (It it=haps[h1].indels.begin();it!=haps[h1].indels.end();it++) cout << it->first << "," << it->second.getString() << ";";
		cout << " IN: " << hap_num_indels[h2] << " SNP: " << hap_num_snps[h2]; cout << " H2:"; for (It it=haps[h2].indels.begin();it!=haps[h2].indels.end();it++) cout << it->first << "," << it->second.getString() << ";";
		cout << endl;
		*/
		
		if (pairs_posterior[h1*nh+h2]>max_ll_indel && (hap_num_candidate_indels[h1]>0 || hap_num_candidate_indels[h2]>0)) {
			max_ll_indel = pairs_posterior[h1*nh+h2];
			max_indel_pair[0] = h1;
			max_indel_pair[1] = h2;
		}
		if (pairs_posterior[h1*nh+h2]>max_ll_noindel && (hap_num_candidate_indels[h1]==0 && hap_num_candidate_indels[h2]==0)) {
			max_ll_noindel = pairs_posterior[h1*nh+h2];
			max_noindel_pair[0] = h1;
			max_noindel_pair[1] = h2;
		}

	}

	// output map-based variant calls
	
	double qual = 0.0;
	double ll_ref = max_ll_noindel;
	qual = - 10.0*( ll_ref - addLogs(max_ll_indel, ll_ref) )/log(10.0);
	cout << "ll_ref: " << ll_ref << " max_ll_indel: " << max_ll_indel << " qual: " << qual << endl;
	if (max_indel_pair[0]==-1 || max_indel_pair[1]==-1) throw string("Could not find indel allele");	
	if (1) {
		int numUnmappedRealigned = 0;
		int hx1 = max_indel_pair[0];
		int hx2 = max_indel_pair[1];
		for (size_t r=0;r<reads.size();r++) {
			if (reads[r].isUnmapped()) {
				if (liks[hx1][r].offHap == false || liks[hx2][r].offHap == false) {
					numUnmappedRealigned++;
				}
			}
		}

		map<int, std::set<AlignedVariant> > indels;
		for (int i=0;i<2;i++) {
			
			const Haplotype & hap = haps[ max_indel_pair[i] ];

			for (It it=hap.indels.begin();it!=hap.indels.end();it++) if (!it->second.isRef() || (it->second.isSNP() && it->second.getString()[3]=='D')){
				indels[it->first].insert(it->second);
			}
		}

		for (map<int, std::set<AlignedVariant> >::const_iterator it = indels.begin();it!=indels.end();it++) {
			
			double msq = 0;
			int numf=0, numr=0, n=0;
			int m =2;
			if (max_indel_pair[0]==max_indel_pair[1]) m = 1;
			for (int i=0;i<m;i++) {
				int h=max_indel_pair[i];
				It iter = haps[h].indels.find(it->first);
				if (iter != haps[h].indels.end() && iter->second.isIndel()) {
				
					for (int r=0;r<nr;r++) {
						bool covered=false, nft=false, nrt=false;

					map<int, bool>::const_iterator it2=liks[h][r].hapIndelCovered.find(it->first);
						if (it2!=liks[h][r].hapIndelCovered.end() && it2->second) covered=true;
					
						if (covered) { // hapHasVar[] is to check whether haplotype truely has variant or only the reference
							if  (reads[r].onReverseStrand) nrt=true; else nft=true;
							double mq=-10*log10(1.0-reads[r].mapQual);
							msq+=mq*mq;
							n++;
						}
						if (nft) numf++;
						if (nrt) numr++;
					}
				}
			}

			if (n!=0) msq=sqrt(msq/double(n)); else msq=0.0;
	
			int was_candidate = 0;
			
			// determine genotype
			const std::set< AlignedVariant> & alleles = it->second;
			string genotype;
			std::set<string> all_genotype;
			string nref_all;

			int vc_f = 0;
			int vc_r = 0;
			if (1) {
				const AlignedVariant & avar = *alleles.begin();
				const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());
				if (av!=NULL) was_candidate=1;
				vc_f += varCoverage[PAV(it->first, avar)].nf;
				vc_r += varCoverage[PAV(it->first, avar)].nr;
			}
			
			string a1="*REF", a2="*REF";
			bool a1_ref=true, a2_ref=true;
			It ita1 = haps[hx1].indels.find(it->first);
			It ita2 = haps[hx2].indels.find(it->first);
					
			if (ita1 != haps[hx1].indels.end() && !ita1->second.isRef()) {
				a1 = ita1->second.getString();
				a1_ref=false;
			}
			if (ita2 != haps[hx2].indels.end() && !ita2->second.isRef()) {
				a2 = ita2->second.getString();
				a2_ref=false;
			}
			all_genotype.insert(a1);
			all_genotype.insert(a2);

			//cout << "a1: " << a1 << " a2: " << a2 << " a1ref " << a1_ref << " a2ref " << a2_ref << endl;

			if (a1_ref && a2_ref) throw string("genotyping error");

			if (a1==a2) {
				genotype = string("1/1");
				nref_all = a1;
			} else {
				
				if (a1_ref) { 
					genotype = string("0/1");
					nref_all = a2;
				} else if (a2_ref) {
					genotype = string("0/1");
					nref_all = a1;
				} else  {
					nref_all = a1+','+a2;
					genotype = string("1/2");
					
					if (1) {
						const AlignedVariant & avar = *alleles.rbegin();
						const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());
						if (av!=NULL) was_candidate=1;
						vc_f += varCoverage[PAV(it->first, avar)].nf;
						vc_r += varCoverage[PAV(it->first, avar)].nr;
					}
				}
			}

			// determine genotype quality
			//cout << "POS: " << it->first << endl;
			double max_ll_altgeno = -HUGE_VAL;
			for (size_t h1=0;h1<nh;h1++) if (filtered[h1]==0) for (size_t h2=h1;h2<nh;h2++) if (filtered[h2]==0) {
				if (!( (h1==hx1 && h2 == hx2) || (h2==hx1 && h1 == hx2))) {
					std::set<string> _alt_genotype;
					 It it2=haps[h1].indels.find(it->first);
					 if (it2 == haps[h1].indels.end() || it2->second.isRef()) {
						 _alt_genotype.insert(string("*REF"));
					 } else {
						 _alt_genotype.insert(it2->second.getString());
					}
					it2=haps[h2].indels.find(it->first);
					 if (it2 == haps[h2].indels.end() || it2->second.isRef()) {
						 _alt_genotype.insert(string("*REF"));
					 } else {
						 _alt_genotype.insert(it2->second.getString());
					}
			//		 cout << "CALLED: " << h1 << " " << h2 << " geno: " <<  *all_genotype.begin() << " " << *all_genotype.rbegin() << " ALT: " << *_alt_genotype.begin() << " " << *_alt_genotype.rbegin() <<  " " << pairs_posterior[h1*nh+h2] << " " << max_ll_altgeno << endl;

				 	if (_alt_genotype != all_genotype && max_ll_altgeno<pairs_posterior[h1*nh+h2]) {
						 max_ll_altgeno = pairs_posterior[h1*nh+h2];
							
					}

				}
			}
			double genoqual = 0.0;
			genoqual = - 10.0*( max_ll_altgeno - addLogs(max_ll_indel, max_ll_altgeno) )/log(10.0);




	

			ostringstream glfs; glfs << genotype << ":" << genoqual;


	 	
			OutputData::Line line(glfData);
			line.set("msg", "ok");
			line.set("index", index);
			line.set("tid", params.tid);
			line.set("analysis_type",string("dip.map"));
			line.set("indidx",0);
			line.set("lpos",leftPos);
			line.set("rpos",rightPos);
			line.set("center_position",candPos);
			line.set("realigned_position",it->first+leftPos);
			line.set("was_candidate_in_window",was_candidate);
			line.set("qual", qual );
			//line.set("post_prob_variant", exp(logp));
			//line.set("est_freq", freq);
			line.set("nref_all", nref_all);
			line.set("num_reads", readidx.size());
			line.set("msq",msq);
			//line.set("numOffAll",numOffBoth);
			//line.set("num_indel",numMappedIndels);
			line.set("num_cover_forward", numf);
			line.set("num_cover_reverse", numr);
			line.set("var_coverage_forward", vc_f);
			line.set("var_coverage_reverse", vc_r);
			line.set("num_unmapped_realigned", numUnmappedRealigned);
			line.set("glf", glfs.str());
			//line.set("likrr", lik[0]);
			//line.set("likrn", lik[1]);
			//line.set("liknn", lik[2]);
			glfData.output(line);
		}
	}
	

	for (map<int, std::set<PAV> >::const_iterator it= allVariantsByPos.begin();it!=allVariantsByPos.end();it++)
	{
		bool has_variants_in_window = 0;
		BOOST_FOREACH(PAV pav, it->second) {
			const AlignedVariant & avar = pav.second;
			const AlignedVariant *av = candidateVariants.findVariant(avar.getStartHap()+leftPos, avar.getType(), avar.getString());
			//cout << "VAR "  << avar.getString() << endl;
			if (av!=NULL) {
				has_variants_in_window=1;
				break;
			}
		}



		std::set<int> support;
		int totnf=0, totnr=0;
		double log5=log(0.5);
		int nf=0;
		int nr=0;

		// compute RMS of mapping qualities

		int pos=it->first;
		int posIdx = posToPosIdx[pos];
		double msq=0.0;
		int n=0;

		//cout << endl;
		//cout << "pos: " << pos << " posIdx " << posIdx << endl;

		typedef std::set<int> IntGenotype;
		map<IntGenotype, double> genLiks;

		typedef map<IntGenotype,double>::iterator IGIt;

		double maxll=-HUGE_VAL;
		int hx1, hx2;

		for (size_t h1=0;h1<nh;h1++) if (filtered[h1]==0) for (size_t h2=h1;h2<nh;h2++) if (filtered[h2]==0) {
			IntGenotype genotype;
			int v1=hapVar[h1*numVarPos+posIdx];
			int v2=hapVar[h2*numVarPos+posIdx];

			genotype.insert(hapVar[h1*numVarPos+posIdx]);
			genotype.insert(hapVar[h2*numVarPos+posIdx]);

			double logPriorPos=0.0;
			//cerr << "FIX THIS!\n" << endl;
//			if ( (v1>0 && varType[v1]==VARSNP) || (v2>0 && varType[v2]==VARSNP)) { logPriorPos += log(params.priorSNP);  } else if ( (v1>0 && varType[v1]==VARINDEL) || (v2>0 && varType[v2]==VARINDEL) ) { logPriorPos+=log(params.priorIndel); };

			AlignedVariant av1, av2;
			if (v1) av1=variants[v1].second; else av1=AlignedVariant("*REF",-1);
			if (v2) av2=variants[v2].second; else av2=AlignedVariant("*REF",-1);

			logPriorPos = getPairPrior(av1,av2,leftPos, candidateVariants);
			double pr=prior[h1*nh+h2]-logPriorPos; // substract prior for this site to obtain likelihood

			//cout << "prior: " << h1 << " " << h2 << ": " << prior[h1*nh+h2] << " logPriorPos " << logPriorPos << endl;
			double ll=pr;
			for (size_t r=0;r<reads.size();r++){
				ll+=log5+addLogs(rl[r*nh+h1],rl[r*nh+h2]);
			}

			//cout << "genotype: " << *(genotype.begin()) << " " << *(genotype.rbegin()) << " lik: " << ll <<endl;


			IGIt igit = genLiks.find(genotype);
			if (igit==genLiks.end()) {
				genLiks[genotype]=ll;
			} else {
				genLiks[genotype]=addLogs(genLiks[genotype],ll);
			}

			if (ll>maxll) {
				maxll=ll;
				hx1=h1;
				hx2=h2;
			}

		}
		//cout << "hx1: " << hx1 << " hx2: " << hx2 << " postprob: " << maxll << endl;

		// see how many unmapped reads were realigned to the MAP haplotypes

		int numUnmappedRealigned = 0;
		for (size_t r=0;r<reads.size();r++) {
			if (reads[r].isUnmapped()) {
				if (liks[hx1][r].offHap == false || liks[hx2][r].offHap == false) {
					numUnmappedRealigned++;
				}
			}
		}

		if (params.outputPooledLikelihoods) {
			string ofLiksFile =  params.fileName;
			ofLiksFile.append(".check.txt");
			ofstream ofLiks(ofLiksFile.c_str());
			if (!ofLiks.is_open()) {
				throw string("Cannot open file for writing");
			}

			// output haplotypes
			//
			ofLiks << "HAPLOTYPES" << endl;
			for (size_t h=0;h<haps.size();h++) {
				ofLiks << h;
				stringstream varss;

				for(map<int, AlignedVariant>::const_iterator it = haps[h].indels.begin(); it != haps[h].indels.end(); it++) {
					varss << leftPos+it->first << "," << it->second.getString() << ";";
				}
				ofLiks << "\t" << varss.str() << endl;
			}

			ofLiks << "READS" << endl;


			for (size_t r=0;r<reads.size();r++) {
				int offBoth = 0;
				if (liks[hx1][r].offHap == true && liks[hx2][r].offHap == true) {
					offBoth =1;
				}
				ofLiks << r << "\t" << bam1_qname(reads[r].getBam()) << "\t" << reads[r].pos << "\t" << reads[r].mapQual;
				for (size_t h=0;h<haps.size();h++) {
					ofLiks << "\t" << liks[h][r].ll;
				}
				for (size_t h=0;h<haps.size();h++) {
					ofLiks << "\t" << int(liks[h][r].offHap);
				}
				ofLiks << endl;

			}
			ofLiks.close();
		}
#define DEBUGDIPLOIDGLF
#ifdef DEBUGDIPLOIDGLF
		if (params.outputPooledLikelihoods) {
			for (size_t r=0;r<reads.size();r++) {
				int mhx1 = 0;
				int mhx2 = 1;
				
				cout  << endl << "**READ** " << r << " " << bam1_qname(reads[r].getBam()) << " mapQual: " << reads[r].mapQual << " liks: " << liks[mhx1][r].ll << " " << liks[mhx2][r].ll << " unMapped: " << reads[r].isUnmapped() << endl;

				if (1) {
					/*
					cout << string(50,' '); cout << reads[r].seq.seq << endl;
					cout << haps[hx1].seq << endl;
					cout << haps[hx2].seq << endl;
					*/
					Read newread(reads[r]);

					/*
					newread.mapQual = 1.0 - 1e-20;
					newread.complement();
					newread.reverse();
					*/

					cout << "first: " << endl;
					ObservationModelFBMaxErr obs(haps[mhx1], newread, leftPos, params.obsParams);
					cout << string(50,' ') << haps[mhx1].seq << endl;
					obs.printAlignment(50);


					cout << "second: " << endl;
					ObservationModelFBMaxErr obs2(haps[mhx2], newread, leftPos, params.obsParams);
					cout << string(50,' ') << haps[mhx2].seq << endl;
					obs2.printAlignment(50);


				} else {
					int hidx = hx1; if (liks[hx2][r].ll>liks[hx1][r].ll) hidx = hx2;
					cout << "hidx: " << hidx << " hx1: " << hx1 << " hx2: " << hx2 << endl;
					ObservationModelFBMaxErr obs(haps[hidx], reads[r], leftPos, params.obsParams);
					cout << string(50,' ') << haps[hidx].seq << endl;
					obs.printAlignment(50);
				}
			}
		}

#endif




		double allmsq=0.0;
		int numMappedIndels=0;

		int nBQT=0, nmmBQT=0;
		double mLogBQ=0.0;
		int nMMLeft=0;
		int nMMRight=0;
		int numOffBoth =0;

		BOOST_FOREACH(int r, readidx) {
			double mq=-10*log10(1.0-reads[r].mapQual);
			allmsq+=(mq*mq);
			//cout << " " << 10*log10(1-reads[r].mapQual);

			int mlidx; double ml=-HUGE_VAL;
			std::set<int> mlis;

			if (liks[hx1][r].offHap && liks[hx2][r].offHap) numOffBoth++;

			if (liks[hx1][r].ll>=liks[hx2][r]) {
				mlidx=hx1;
				ml=liks[hx1][r].ll;
			} else {
				mlidx=hx2;
				ml=liks[hx2][r].ll;
			}

			mlis.insert(mlidx);

			bool nrt=false, nft=false;

			map<int,bool>::const_iterator it;
			BOOST_FOREACH(int h, mlis) {
				bool covered=false;
				numMappedIndels += int(liks[h][r].indels.size());
				nBQT+=liks[h][r].nBQT;
				nmmBQT+=liks[h][r].nmmBQT;
				mLogBQ+=liks[h][r].mLogBQ;
				if (liks[h][r].nMMLeft>=2) nMMLeft++;
				if (liks[h][r].nMMRight>=2) nMMRight++;


				map<int, AlignedVariant>::const_iterator hit=haps[h].indels.find(pos);
				if (hit->second.isIndel()) {
					it=liks[h][r].hapIndelCovered.find(pos);
					if (it!=liks[h][r].hapIndelCovered.end() && it->second) covered=true;
				} else if (hit->second.isSNP()) {
					it=liks[h][r].hapSNPCovered.find(pos);
					if (it!=liks[h][r].hapSNPCovered.end() && it->second) covered=true;
				}
				if (covered) { // hapHasVar[] is to check whether haplotype truely has variant or only the reference
					if  (reads[r].onReverseStrand) nrt=true; else nft=true;
					double mq=-10*log10(1.0-reads[r].mapQual);
					msq+=mq*mq;
					n++;
				}
			}
			if (nft) nf++;
			if (nrt) nr++;
		} // foreach read r
	//	cout << endl;

		if (n!=0) msq=sqrt(msq/double(n)); else msq=0.0;
		totnf+=nf;
		totnr+=nr;

		allmsq=(readidx.size()!=0)?sqrt(allmsq/double(readidx.size())):0;


		// recode variant indexes to glf/vcf indexes
		int nidx=1;
		map <int, int> toVCFidx;
		vector<string> alleles;
		alleles.push_back("R");
		toVCFidx[0]=0;

		ostringstream oAlleles;
		ostringstream oCovForward;
		ostringstream oCovReverse;
		int first=1;
		for (size_t h=0;h<nh;h++) {
			int v = hapVar[h*numVarPos+posIdx];
			if (v!=0) {
				map<int,int>::iterator tit = toVCFidx.find(v);
				if (tit==toVCFidx.end()) {
					toVCFidx[v]=nidx++;
					alleles.push_back(variants[v].second.getString());
					string str=(first==1) ? string(""):string(",");

					oAlleles << str << variants[v].second.getString();
					oCovForward << str << varCoverage[variants[v]].nf;
					oCovReverse << str << varCoverage[variants[v]].nr;
					first = 0;
				}
			}
		}


		// compute genotype posterior

		ostringstream o;

		first=1;
		for (map<IntGenotype, double>::iterator git=genLiks.begin();git!=genLiks.end();git++) {
			int v1 = *(git->first.begin());
			int v2 = *(git->first.rbegin());
			int a1 = toVCFidx[ *(git->first.begin()) ];
			int a2 = toVCFidx[ *(git->first.rbegin()) ];

			string str=(first==1) ? string(""):string(",");
			o << str << a1 << "/" << a2 << ":" << git->second;

			double logPrior=0.0;
			if ( (v1>0 && varType[v1]==VARSNP) || (v2>0 && varType[v2]==VARSNP)) { logPrior += log(params.priorSNP);  } else if ( (v1>0 && varType[v1]==VARINDEL) || (v2>0 && varType[v2]==VARINDEL) ) { logPrior+=log(params.priorIndel); };

			git->second -= logPrior;
			first = 0;
		}





		if (params.outputGLF) {
			OutputData::Line line(glfData);
			line.set("msg", "ok");
			line.set("index", index);
			line.set("tid", params.tid);
			line.set("analysis_type",program);
			line.set("indidx",0);
			line.set("lpos",leftPos);
			line.set("rpos",rightPos);
			line.set("center_position",candPos);
			line.set("realigned_position",pos+leftPos);
			line.set("was_candidate_in_window",has_variants_in_window);
			line.set("logZ", maxll );
			//line.set("post_prob_variant", exp(logp));
			//line.set("est_freq", freq);
			line.set("nBQT", nBQT);
			line.set("nmmBQT", nmmBQT);
			line.set("mLogBQ", mLogBQ/double(nBQT));
			line.set("nMMLeft", nMMLeft);
			line.set("nMMRight", nMMRight);
			line.set("nref_all", oAlleles.str());
			line.set("num_reads", readidx.size());
			line.set("msq",allmsq);
			line.set("numOffAll",numOffBoth);
			line.set("num_indel",numMappedIndels);
			line.set("num_cover_forward", nf);
			line.set("num_cover_reverse", nr);
			
			line.set("var_coverage_forward", oCovForward.str());
			line.set("var_coverage_reverse", oCovReverse.str());


			line.set("glf", o.str());
			line.set("num_unmapped_realigned", numUnmappedRealigned);
			//line.set("likrr", lik[0]);
			//line.set("likrn", lik[1]);
			//line.set("liknn", lik[2]);
			glfData.output(line);
		}
	// glfOutput <<  PRID << " "  << b << " " << params.tid << " " << candPos << " " << pav.first+leftPos << " " << pav.second.getString() << " " << reads.size() <<  " " << msq << " "  << nf << " " << nr << " " << lik[0] << " " << lik[1] << " " << lik[2] << endl;


  }

}
#endif

void DetInDel::estimateHaplotypeFrequencies(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<double> & hapFreqs)
{

	// estimate haplotype frequencies using EM
	hapFreqs.clear();

	size_t nh=haps.size();
	size_t nr=reads.size();




	vector<double> rl(nh*nr,0.0); // read given haplotype likelihoods

	vector<double> z(nh*nr,0.0); // expectations of read-haplotype indicator variables
	vector<double> pi(nh); // log of haplotype frequencies
	vector<double> nk(nh,0.0); // counts for each haplotype

	hapFreqs=nk;

	// initialize frequencies
	for (size_t h=0;h<nh;h++) pi[h]=1.0/double(nh);

	for (size_t h=0;h<nh;h++) for (size_t r=0;r<nr;r++) {
		// initialize read-haplotype likelihoods
		rl[h*nr+r]=liks[h][r].ll;

		// initialize expectations of indicator variables
		z[h*nr+r]=0.5;
	}


	bool converged=false;
	double tol=params.EMtol;

	double eOld=-HUGE_VAL, eNew;

	cout << "EM HapFreqs:";

	int iter=0;
	while (!converged) {

		// compute expectation of indicator variables
		for (size_t h=0;h<nh;h++) nk[h]=0.0;

		int idx=0;
		for (size_t r=0;r<nr;r++) {
			double lognorm=-HUGE_VAL;
			// compute responsibilities
			for (size_t h=0;h<nh;h++) {
				z[h*nr+r]=pi[h]+(rl[h*nr+r]);
				lognorm=addLogs(lognorm, z[h*nr+r]);
			}
			// normalize and compute counts
			for (size_t h=0;h<nh;h++) {
				z[nr*h+r]-=lognorm;
				z[nr*h+r]=exp(z[nr*h+r]);

				nk[h]+=z[nr*h+r];
			}
		}

		// compute frequencies

		for (size_t h=0;h<nh;h++) {
			double nph=nk[h]/nr;
			pi[h]=log(nph);
		}


		idx=0;
		eNew=0.0;
		for (size_t h=0;h<nh;h++) {

		for (size_t r=0;r<nr;r++) {
			// compute responsibilities
				eNew+=z[idx]*( pi[h]+rl[idx]);
				idx++;
			}
		}
		//cout << "[" << eNew << "]" << endl;
		//
		if (eOld>eNew) throw string("EM Error in estimateHapFreqs");
		converged=(fabs(eOld-eNew))<tol || iter>25;

		eOld=eNew;


		iter++;
	}

	for (size_t h=0;h<nh;h++) { cout << " " << exp(pi[h]); }
	cout << endl;

	// output haplotypes and estimated frequencies

	for (size_t h=0;h<nh;h++) hapFreqs[h]=exp(pi[h]);
}


void DetInDel::computePairLikelihoods(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, vector<HapPairLik> & likPairs, bool usePrior, const AlignedCandidates & candidateVariants, int leftPos)
{
//	cout << "Computing pair likelihoods...\n";
	likPairs.clear();
	size_t lh=haps.size();
	likPairs.reserve(lh*(lh/2));
	//const double log10=log(10);
	double maxLL=-HUGE_VAL; int hpm, hmm;
	size_t midx;
	for (size_t hp=0;hp<lh;hp++) for (size_t hm=hp;hm<lh;hm++) {
		double ll=0.0;
		HapPairLik hpl;
		hpl.numIndFirst=0;
		hpl.numIndSecond=0;
		hpl.numOffBoth=0;
		hpl.numOffBothError=0.0;
		hpl.numFirst=0;
		hpl.numSecond=0;
		hpl.h1=hp;
		hpl.h2=hm;
		for (size_t r=0;r<reads.size();r++) {
			ll+=(addLogs(liks[hp][r].ll,liks[hm][r].ll)+log(.5));
			const MLAlignment & ml1=liks[hp][r];
			const MLAlignment & ml2=liks[hm][r];

			// record which indel was on which haplotype
			if (!ml1.offHap && ml1.ll>ml2.ll && ml1.indels.size()!=0) hpl.numIndFirst++;
			else if (!ml2.offHap && ml2.ll>ml1.ll && ml2.indels.size()!=0) hpl.numIndSecond++;
			if (ml1.offHapHMQ && ml2.offHapHMQ) { hpl.numOffBoth++; hpl.numOffBothError+=reads[r].mapQual; };
			if (ml1.ll>=ml2.ll) hpl.numFirst++;
			if (ml2.ll>=ml1.ll) hpl.numSecond++;

	//		determine read coverage of all the variants
	//		TODO this part is really slow!
			for (map<int, bool>::const_iterator it=ml1.hapIndelCovered.begin();it!=ml1.hapIndelCovered.end();it++) if (it->second) if (ml1.ll>=ml2.ll) { if (reads[r].onReverseStrand) hpl.hapIndelCoverage1[it->first].nr++; else hpl.hapIndelCoverage1[it->first].nf++; }
			for (map<int, bool>::const_iterator it=ml2.hapIndelCovered.begin();it!=ml2.hapIndelCovered.end();it++) if (it->second) if (ml2.ll>=ml1.ll) { if (reads[r].onReverseStrand) hpl.hapIndelCoverage2[it->first].nr++; else hpl.hapIndelCoverage2[it->first].nf++; }
			for (map<int, bool>::const_iterator it=ml1.hapSNPCovered.begin();it!=ml1.hapSNPCovered.end();it++) if (it->second) if (ml1.ll>=ml2.ll) { if (reads[r].onReverseStrand) hpl.hapSNPCoverage1[it->first].nr++; else hpl.hapSNPCoverage1[it->first].nf++; }
			for (map<int, bool>::const_iterator it=ml2.hapSNPCovered.begin();it!=ml2.hapSNPCovered.end();it++) if (it->second) if (ml2.ll>=ml1.ll) { if (reads[r].onReverseStrand) hpl.hapSNPCoverage2[it->first].nr++; else hpl.hapSNPCoverage2[it->first].nf++; }


		}
		if (usePrior) ll+=getHaplotypePrior(haps[hp], haps[hm], leftPos, candidateVariants);

		hpl.ll=ll;
		if (ll>maxLL) {
			maxLL=ll;
			hpm=hp;
			hmm=hm;
			midx=likPairs.size();
		}
	//	cout << "hp: " << hp << " hm: " << hm << " ll: " << ll << endl;
		likPairs.push_back(hpl);
	}


//	cout << "ML hap: " << hpm << " " << hmm << " midx: " << midx << endl;
	/*
	cout << haps[hpm] << endl;
	cout << "indels: "; for (map<int, AlignedVariant>::const_iterator it=haps[hpm].indels.begin();it!=haps[hpm].indels.end();it++) {
		cout << "[" << it->first << "," << it->second.getString() << "]";
	}
	cout << endl;
	cout << "coverage: ";
	for (map<int, VariantCoverage>::const_iterator it=likPairs[midx].hapIndelCoverage1.begin();it!=likPairs[midx].hapIndelCoverage1.end();it++) {
		cout << "[" << it->first << "," << it->second.nf << "," << it->second.nr << "]";
	}
	cout << endl;
	cout << haps[hmm] << endl;
	cout << "indels: "; for (map<int, AlignedVariant>::const_iterator it=haps[hmm].indels.begin();it!=haps[hmm].indels.end();it++) {
		cout << "[" << it->first << "," << it->second.getString() << "]";
	}
	cout << endl;
	for (map<int, VariantCoverage>::const_iterator it=likPairs[midx].hapIndelCoverage2.begin();it!=likPairs[midx].hapIndelCoverage2.end();it++) {
		cout << "[" << it->first << "," << it->second.nf << "," << it->second.nr << "]";
	}

	cout << endl;
	*/

	class SortFunc {
	public:
		static bool sortFunc(const HapPairLik & hpl1, const HapPairLik & hpl2)
		{
			// sort in decreasing order
			return hpl1.ll>hpl2.ll;
		}
	};
	sort(likPairs.begin(), likPairs.end(),SortFunc::sortFunc);
}

void DetInDel::statisticsHaplotypePair(const vector<Haplotype> & haps, const vector<Read> & reads, const vector<vector<MLAlignment> > & liks, HapPairLik & hpl, OutputData::Line & line)
{
	hpl.numIndFirst=0;
	hpl.numIndSecond=0;
	hpl.numOffBoth=0;
	hpl.numOffBothError=0.0;
	hpl.numFirst=0;
	hpl.numSecond=0;
	int hp=hpl.h1;
	int hm=hpl.h2;
	for (size_t r=0;r<reads.size();r++) {
		const MLAlignment & ml1=liks[hp][r];
		const MLAlignment & ml2=liks[hm][r];

		if (!ml1.offHap && ml1.ll>ml2.ll && ml1.indels.size()!=0) hpl.numIndFirst++;
		else if (!ml2.offHap && ml2.ll>ml1.ll && ml2.indels.size()!=0) hpl.numIndSecond++;
		if (ml1.offHapHMQ && ml2.offHapHMQ) { hpl.numOffBoth++; hpl.numOffBothError+=reads[r].mapQual; };
		if (ml1.ll>=ml2.ll) hpl.numFirst++;
		if (ml2.ll>=ml1.ll) hpl.numSecond++;

//		determine read coverage of all the variants
//		TODO this part is really slow!
		for (map<int, bool>::const_iterator it=ml1.hapIndelCovered.begin();it!=ml1.hapIndelCovered.end();it++) if (it->second) if (ml1.ll>=ml2.ll) { if (reads[r].onReverseStrand) hpl.hapIndelCoverage1[it->first].nr++; else hpl.hapIndelCoverage1[it->first].nf++; }
		for (map<int, bool>::const_iterator it=ml2.hapIndelCovered.begin();it!=ml2.hapIndelCovered.end();it++) if (it->second) if (ml2.ll>=ml1.ll) { if (reads[r].onReverseStrand) hpl.hapIndelCoverage2[it->first].nr++; else hpl.hapIndelCoverage2[it->first].nf++; }
		for (map<int, bool>::const_iterator it=ml1.hapSNPCovered.begin();it!=ml1.hapSNPCovered.end();it++) if (it->second) if (ml1.ll>=ml2.ll) { if (reads[r].onReverseStrand) hpl.hapSNPCoverage1[it->first].nr++; else hpl.hapSNPCoverage1[it->first].nf++; }
		for (map<int, bool>::const_iterator it=ml2.hapSNPCovered.begin();it!=ml2.hapSNPCovered.end();it++) if (it->second) if (ml2.ll>=ml1.ll) { if (reads[r].onReverseStrand) hpl.hapSNPCoverage2[it->first].nr++; else hpl.hapSNPCoverage2[it->first].nf++; }

		// record which indel was on which haplotype
	}

	line.set("num_off_hap", hpl.numOffBoth);
	line.set("num_mapped_to_first",hpl.numFirst);
	line.set("num_mapped_to_second",hpl.numSecond);
}



void parseRegionString(const string & region, int & start, int & end)
{
	string filtered;
	for(size_t x=0;x<region.size();x++) {
		char c=region[x];
		if (c=='-') filtered+=' ';
		else if (c!=',') filtered+=c;
	}
	istringstream is(filtered);
	string e; is >> e;
	if (!from_string(start,e,std::dec)) throw string("Cannot parse region for start!");
	is >> e;
	if (!from_string(end,e,std::dec)) throw string("Cannot parse region end!");
}

void getParameters(po::variables_map & vm, DetInDel::Parameters & params)
{
	params.maxHap=vm["maxHap"].as<uint32_t> ();
	params.maxReads=vm["maxRead"].as<uint32_t> ();
	params.width=vm["width"].as<uint32_t> ();
	params.mapQualThreshold=vm["mapQualThreshold"].as<double>();
	params.skipMaxHap=vm["skipMaxHap"].as<uint32_t>();
	//params.glfNumHap=vm["glfNumHap"].as<uint32_t>();
	params.inferenceMethod=vm["inferenceMethod"].as<string>();
	params.minReadOverlap=vm["minReadOverlap"].as<uint32_t>();
	params.maxReadLength=vm["maxReadLength"].as<uint32_t>();
	//params.scaleErr=vm["mapScaleError"].as<double>();
	//params.minCount=vm["minCount"].as<uint32_t>();
	params.maxHapReadProd=vm["maxHapReadProd"].as<uint32_t>();

	params.priorSNP=vm["priorSNP"].as<double>();
	params.priorIndel=vm["priorIndel"].as<double>();
	params.bayesa0=vm["bayesa0"].as<double>();
	params.bayesType=vm["bayesType"].as<string>();

	


	if (vm.count("ref")) {
		params.alignAgainstReference=true;
		params.refFileName=vm["ref"].as<string>();
	} else {
		params.alignAgainstReference=false;
	}

	params.obsParams.pError=vm["pError"].as<double>();
	params.obsParams.pMut=vm["pMut"].as<double>();

	//params.obsParams.baseQualThreshold=vm["baseQualThreshold"].as<double>();
	//	params.obsParams.fixedBaseQual=vm["fixedBaseQual"].as<double>();
	params.obsParams.maxLengthIndel=vm["maxLengthIndel"].as<int>();
	params.obsParams.maxLengthDel=params.obsParams.maxLengthIndel;
	params.obsParams.mapQualThreshold=vm["capMapQualThreshold"].as<double>();
	params.obsParams.capMapQualFast=vm["capMapQualFast"].as<double>();
	//params.obsParams.scaleErr=vm["obsScaleError"].as<double>();
	//params.obsParams.numE=vm["numE"].as<int>();
	params.obsParams.padCover = vm["flankRefSeq"].as<int>();
	params.obsParams.maxMismatch = vm["flankMaxMismatch"].as<int>();
	params.checkAllCIGARs=vm["checkAllCIGARs"].as<int>();

	params.varFileIsOneBased=vm.count("varFileIsOneBased")?true:false;
	params.outputRealignedBAM=vm.count("outputRealignedBAM")?true:false;
	params.analyzeLowFreq=vm.count("compareReadHap")?true:false;
	params.analyzeLowFreqDiffThreshold=vm["compareReadHapThreshold"].as<double>();
	params.showHapDist=vm.count("showEmpirical")?true:false;
	params.showCandHap=vm.count("showCandHap")?true:false;
	params.showReads=vm.count("showReads")?true:false;
	params.quiet=vm.count("quiet")?true:false;
	params.computeML=vm.count("computeML")?true:false;
	params.computeMAP=vm.count("computeMAP")?true:false;
	params.doDiploid=vm.count("doDiploid")?true:false;

	params.filterHaplotypes=vm.count("filterHaplotypes")?true:false;

	params.printCallsOnly=vm.count("printCallsOnly")?true:false;
	params.estimateHapFreqs=vm.count("doPooled")?true:false;
	params.outputPooledLikelihoods=vm.count("opl")?true:false;
	params.showHapAlignments=vm.count("showHapAlignments")?true:false;
	if (vm.count("filterReadAux")) params.filterReadAux=vm["filterReadAux"].as<string>();
	if (vm.count("processRealignedBAM")) params.processRealignedBAM=vm["processRealignedBAM"].as<string>();

	params.slower=vm.count("faster")?false:true;
	params.changeINStoN=vm.count("changeINStoN")?true:false;
	params.outputGLF=true;


	// removed options
/*
	params.outputRealignedBAM=vm.count("outputRealignedBAM")?true:false;
	params.obsParams.modelType=vm["modelType"].as<string>();
	params.mapUnmappedReads=vm.count("mapUnmapped")?true:false;
	params.obsParams.mapUnmappedReads=vm.count("mapUnmapped")?true:false;
	params.obsParams.pFirstgLO=vm["pFirstgLO"].as<double>();
	//params.numOutputTopHap=vm["numOutputTopHap"].as<int>();

*/

}



int smain(int argc, char *argv[])
{
	if (1) {
		Haplotype hap;
		Read read;

		//hap.seq=          "ATCGTGTAGCTCTCTGGCTGGCTAGCTGATTGGCTCTTGCC";
		//read.seq.seq=              "CTCTCTGGCTGGCTAGCGAT";
		Haplotype ref;
		//                 012345678901234567890123456789
		hap.seq=          "ATCGATTCGTGATATATATATTCAATGTAGTCGCTAG";
		read.seq.seq=     "ATCGATTCGTGATAATATTCAATGTAGTCGCTAG";


		//hap.seq=          "ATCGATTCGTGATATATATATTCAATGTAGTCGCTAG";
		//read.seq.seq=     "ATCGATTCGTGATATATATATAATTCAATGTAGTCGCTAG";

		//                 012345678901234567890123456789
		//hap.seq=          "ATCGATTCGTGTTTTTTCAATGTAGTCGCTAG";
		//read.seq.seq=     "ATCGATTCGTGTTTTTCAATGTAGTCGCTAG";

		read.mapQual=1-1e-16;

		ObservationModelParameters obsParams;
		read.setAllQual(0.99);

		ObservationModelFBMaxErr omfbe(hap, read, 0, obsParams);
		/*
		ObservationModelS oms(hap, read, 0, obsParams);
		HapHash hash(4, hap);

		oms.align(hash);
		*/

		double ll= omfbe.calcLikelihood();

		cout << "ll: " << ll << endl;
		cout << string(50,' ') << hap.seq << endl;
		omfbe.printAlignment(50);
	}
	if (0) {
		Haplotype hap;
		Read read;

		//hap.seq=          "ATCGTGTAGCTCTCTGGCTGGCTAGCTGATTGGCTCTTGCC";
		//read.seq.seq=              "CTCTCTGGCTGGCTAGCGAT";

		hap.seq=          "AAAATCACCAACACTTCATAATCTATTTTTTCCCCTGAGGAACTTCCTAAAATGAATAAAAAAAAACCCCAGCCACATCTGCATTTGCAAACAGGAAACTCTGCAAGCCATACTAAGACCAAAGCTTAGTT";
		read.seq.seq=     "CAAACAGGAAACTCTGCAAGCCATACTAAGACCAAAGCTTAGTTA";


		read.mapQual=1-1e-16;

		ObservationModelParameters obsParams;
		read.setAllQual(0.99);

		ObservationModelFBMaxErr omfbe(hap, read, 0, obsParams);
		/*
		ObservationModelS oms(hap, read, 0, obsParams);
		HapHash hash(4, hap);

		oms.align(hash);
		*/

		double ll= omfbe.calcLikelihood();

		cout << "ll: " << ll << endl;
		cout << string(50,' ') << hap.seq << endl;
		omfbe.printAlignment(50);
	}

	return 0;

}






#ifdef DINDEL
int main(int argc, char *argv[])
{
	po::options_description which("[Required] Program option");
	which.add_options()
	("analysis", po::value<string>()->default_value("indels"),"Analysis type:\n"
													          "getCIGARindels:  Extract indels from CIGARs of mapped reads, and infer libary insert size distributions\n"
															  "indels: infer indels\n"
															  "realignCandidates: Realign/reposition candidates in candidate file\n");

	po::options_description required("[Required] ");
	required.add_options()
	("ref", po::value<string>(),"fasta reference sequence (should be indexed with .fai file)")
	("outputFile", po::value<string>(),"file-prefix for output results");

	po::options_description baminput("[Required] BAM input. Choose one of the following");
	baminput.add_options()
	("bamFile",po::value<string>(), "read alignment file (should be indexed)")
	("bamFiles",po::value<string>(), "file containing filepaths for BAMs to be jointly analysed (not possible for --analysis==indels");


	po::options_description regioninput("[Required for analysis == getCIGARindels]: \nRegion to be considered for extraction of candidate indels.");
	regioninput.add_options()
	("region", po::value<string>(),"region to be analysed in format start-end, eg. 1000-2000")
	("tid", po::value<string>(),"target sequence (eg 'X') ");

	po::options_description varfileinput("[Required for analysis == indels]");
	varfileinput.add_options()
	("varFile", po::value<string>(), "file with candidate variants to be tested.")
	("varFileIsOneBased", "coordinates in varFile are one-based");

	po::options_description output_options("Output options");
	output_options.add_options()
	("outputRealignedBAM", "output BAM file with realigned reads")
	("processRealignedBAM", po::value<string>(),"ABSOLUTE path to script to process realigned BAM file")
	//("outputGLF", "outputGLF for individuals in each bam file")
	("quiet", "quiet output");
	//("printCallsOnly", "print only genotypes where call_lik_ref>0.0001 (only affects --single)");

	po::options_description single_analysis("parameters for analysis==indels option");
	single_analysis.add_options()
	("doDiploid", "analyze data assuming a diploid sequence")
	("doPooled", "estimate haplotype frequencies using Bayesian EM algorithm.\nMay be applied to single individual and pools.");

	po::options_description analysis_opt("General algorithm parameters");
	analysis_opt.add_options()
	//("mapUnmapped", "remap unmapped reads for which mate is mapped")
	("faster","use faster but less accurate ungapped read-haplotype alignment model")
	("filterHaplotypes","prefilter haplotypes based on coverage")
	("flankRefSeq",po::value<int>()->default_value(2),"#bases of reference sequence of indel region")
	("flankMaxMismatch",po::value<int>()->default_value(2),"max number of mismatches in indel region")
	("priorSNP", po::value<double>()->default_value(1.0/1000.0), "prior probability of a SNP site")
	("priorIndel", po::value<double>()->default_value(1.0/10000.0), "prior probability of a detected indel not being a sequencing error")
	("width", po::value<uint32_t>()->default_value(60), "number of bases to left and right of indel")
	("maxHap", po::value<uint32_t>()->default_value(8), "maximum number of haplotypes in likelihood computation")
	("maxRead", po::value<uint32_t>()->default_value(10000), "maximum number of reads in likelihood computation")
	("mapQualThreshold", po::value<double>()->default_value(0.99), "lower limit for read mapping quality")
	("capMapQualThreshold", po::value<double>()->default_value(100.0), "upper limit for read mapping quality in observationmodel_old (phred units)")
	("capMapQualFast", po::value<double>()->default_value(45.0), "cap mapping quality in alignment using fast ungapped method\n (WARNING: setting it too high (>50) might result in significant overcalling!)")
	("skipMaxHap", po::value<uint32_t>()->default_value(200), "skip computation if number of haplotypes exceeds this number")
	//("glfNumHap", po::value<uint32_t>()->default_value(5), "number of haplotypes per glf-class")
	//("numOutputTopHap", po::value<int>()->default_value(5), "number of haplotype pairs output to haplotype file")
	("minReadOverlap", po::value<uint32_t>()->default_value(20),"minimum overlap between read and haplotype")
	("maxReadLength", po::value<uint32_t>()->default_value(500),"maximum length of reads")
	("minCount", po::value<uint32_t>()->default_value(1), "minimum number of WS observations of indel")
	("maxHapReadProd",po::value<uint32_t>()->default_value(10000000), "skip if product of number of reads and haplotypes exceeds this value")
	("changeINStoN", "change sequence of inserted sequence to 'N', so that no penalty is incurred if a read mismatches the inserted sequence");
	po::options_description pooled_analysis("parameters for --pooled option");
	pooled_analysis.add_options()
	("bayesa0", po::value<double>()->default_value(0.001), "Dirichlet a0 parameter haplotype frequency prior")
	("bayesType",po::value<string>()->default_value("singlevariant"), "Bayesian EM program type (all or singlevariant or priorpersite)");


	po::options_description option_filter("General algorithm filtering options");
	option_filter.add_options()
	("checkAllCIGARs",po::value<int>()->default_value(1),"include all indels at the position of the call site")
	("filterReadAux", po::value<string>(), "match string for exclusion of reads based on auxilary information");


	po::options_description obsModel("Observation model parameters");
	obsModel.add_options()
	("pError", po::value<double>()->default_value(5e-4), "probability of a read indel")
	//("modelType", po::value<string>()->default_value("probabilistic"), "probabilistic/threshold")
	("pMut", po::value<double>()->default_value(1e-5), "probability of a mutation in the read")
	("maxLengthIndel", po::value<int>()->default_value(5), "maximum length of a _sequencing error_ indel in read [not for --faster option]");
	//("pFirstgLO",po::value<double>()->default_value(0.01),"probability of transition from off the haplotype to on the haplotype");

	po::options_description libParams("Library options");
	libParams.add_options()
	("libFile", po::value<string>(), "file with library insert histograms (as generated by --analysis getCIGARindels)");


	po::options_description miscAnalysis("Misc results analysis options");
	miscAnalysis.add_options()
	("compareReadHap",  "compare likelihood differences in reads against haplotypes")
	("compareReadHapThreshold", po::value<double>()->default_value(0.5), "difference threshold for viewing")
	("showEmpirical", "show empirical distribution over nucleotides")
	("showCandHap", "show candidate haplotypes for fast method")
	("showHapAlignments","show for each haplotype which reads map to it")
	("showReads","show reads")
	("inferenceMethod",po::value<string>()->default_value("empirical"), "inference method")
	("opl","output likelihoods for every read and haplotype");

	required.add(which).add(baminput).add(regioninput).add(varfileinput).add(output_options).add(single_analysis).add(analysis_opt).add(pooled_analysis).add(option_filter).add(obsModel).add(libParams).add(miscAnalysis);

	po::variables_map vm;

	try {
			po::store(po::parse_command_line(argc, argv, required), vm);
	} catch (boost::program_options::error) {
			cout << "Error parsing input options. Usage: \n\n" << required <<"\n";
			exit(1);
	}
	po::notify(vm);

	// analysis
	if (!(vm.count("analysis"))) {
		cerr << "Error: Specify which analysis (--analysis) is required." << endl;
		exit(1);
	}

	 // required
	if (!(vm.count("ref") && vm.count("outputFile"))) {
		cerr << "Error: One of the following options was not specified:  --ref --tid or --outputFile" << endl;
		exit(1);
	}

	if (vm.count("getCIGARindels") && vm.count("region") && !vm.count("tid")) {
		cerr << "--tid must be specified if analysis==getCIGARindels and --region option is used. " << endl;
		exit(1);
	}
//#define DEBUGGING
#ifndef DEBUGGING
	try {
#endif
		// extract required parameters
		string file;
		int multipleFiles=0;
		string analysis=vm["analysis"].as<string>();

		// baminput
		if (analysis=="indels" || analysis=="getCIGARindels") {
			if (!(vm.count("bamFile") || vm.count("bamFiles"))) {
					cerr << "Error: Specify either --bamFile or --bamFiles." << endl;
					exit(1);
			}

			if (vm.count("bamFile")) {
			  file=vm["bamFile"].as< string >();
			  cout << "Reading BAM file: " << file << endl;
			} else if (vm.count("bamFiles")) {
			  file=vm["bamFiles"].as<string>();
			  multipleFiles=1;
			}
		}

		string faFile=vm["ref"].as<string>();
		string outputFile=vm["outputFile"].as< string >();

		string modelType="probabilistic"; //vm["modelType"].as< string >();
		DetInDel::Parameters params(string("1"), outputFile, modelType);
		getParameters(vm, params);

		if (analysis=="getCIGARindels") {
			GetCandidatesFromCIGAR gcfc;
			string outputFile=vm["outputFile"].as< string >();// outputFile.append(".variants.txt");
			fasta::Fasta fa(faFile);
			if (vm.count("region")) {
				string tid=vm["tid"].as<string>();
				string region=vm["region"].as<string>();
				int start, end;
				parseRegionString(region, start, end);
				DetInDel detInDel(file, params, multipleFiles);
				const vector<MyBam *> &  bams = detInDel.getMyBams();

				cout << "Getting indels from CIGARs in mapped reads from region " << tid << ":" << start << "-" << end << endl;
				gcfc.getIndelFromCIGARRegion(bams,tid, start, end, outputFile, fa);

			} else {
				if (multipleFiles) {
					cerr << "Can extract the full set of indels from only BAM file at a time." << endl;
					exit(1);
				}
				gcfc.get(file, outputFile, faFile);
			}
		} else if (analysis=="indels") {
			if (!vm.count("varFile")) {
				cerr << "Please specify the file with the candidate variants." << endl;
				exit(1);
			}

			string varFile = vm["varFile"].as<string>();

			DetInDel detInDel(file, params, multipleFiles);

			if (vm.count("libFile")) {
				cout << "Detected library file..." << endl;
				detInDel.params.mapUnmappedReads=true;
				detInDel.params.obsParams.mapUnmappedReads=true;
				detInDel.addLibrary(vm["libFile"].as<string>());
			}
			detInDel.params.print();


			detInDel.detectIndels(varFile);
		} else if (analysis == "realignCandidates") {
			GetCandidatesFromCIGAR gcfc;
			string outputFile=vm["outputFile"].as< string >(); outputFile.append(".variants.txt");

			if (!vm.count("varFile")) {
				cerr << "Please specify the file with the candidate variants." << endl;
				exit(1);
			}

			string varFile = vm["varFile"].as<string>();

			if (varFile == outputFile) {
				cerr << "outputFile is same as variant file used for input!" << endl;
				exit(1);
			}

			gcfc.realignCandidateFile(varFile, params.varFileIsOneBased,outputFile, faFile);
		} else {
			cerr << "Unrecognized --analysis option." << endl;
			exit(1);
		}
#ifndef DEBUGGING
   }
   catch (string s) {
	cout << "Exception: " << s << endl;
    exit(1);
   }
#endif
   return 0;
}
#endif

