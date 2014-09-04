 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: graph_align_tcoffee_msa.h 1911 2008-05-02 09:28:04Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_MSA_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_MSA_H

namespace SEQAN_NAMESPACE_MAIN
{

template<typename TAlphabet, typename TScore>
struct MsaOptions {
public:
	// Rescore segment matches after refinement
	bool rescore;

	// Output format
	// 0: Fasta
	// 1: Msf
	unsigned int outputFormat;

	// Scoring object
	TScore sc;

	// All methods to compute a guide tree
	// 0: Neighbor-joining
	// 1: UPGMA single linkage
	// 2: UPGMA complete linkage
	// 3: UPGMA average linkage
	// 4: UPGMA weighted average linkage
	unsigned int build;

	// All methods to compute segment matches
	// 0: global alignment
	// 1: local alignments
	// 2: overlap alignments
	// 3: longest common subsequence
	String<unsigned int> method;

	// Various input and output file names
	String<std::string> alnfiles;		// External alignment files
	String<std::string> libfiles;		// T-Coffee library files
	String<std::string> blastfiles;		// Blast match files
	String<std::string> mummerfiles;	// MUMmer files
	std::string outfile;				// Output file name
	std::string seqfile;				// Sequence file name
	std::string infile;					// Alignment file for alignment evaluation
	std::string treefile;				// Guide tree file
	
	// Initialization
	MsaOptions() : rescore(true), outputFormat(0), build(0) {}
};

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore>
inline void
evaluateAlignment(MsaOptions<TAlphabet, TScore> const& msaOpt) {
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TAlphabet> TSequence;
	typedef typename Size<TSequence>::Type TSize;
	StringSet<TSequence, Owner<> > origStrSet;
	StringSet<String<char> > names;

	// Read the sequences
	std::fstream strm;
	strm.open(msaOpt.infile.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
	read(strm,origStrSet,names,FastaAlign());	
	strm.close();

	// Make a dependent StringSet
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	TDepSequenceSet strSet(origStrSet);
	
	// Read the alignment
	typedef String<Fragment<> > TFragmentString;
	String<TScoreValue> scores;
	TFragmentString matches;
	std::fstream strm_lib;
	strm_lib.open(msaOpt.infile.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
	read(strm_lib,matches, scores, names, FastaAlign());	
	strm_lib.close();

	// Build the alignment graph
	typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
	TGraph g(strSet);
	buildAlignmentGraph(matches, g, FrequencyCounting() );

	// Print the scoring information
	TScoreValue gop = msaOpt.sc.data_gap_open;
	TScoreValue gex = msaOpt.sc.data_gap_extend;
	std::cout << "Scoring parameters:" << std::endl;
	std::cout << "*Gap opening: " << gop << std::endl;
	std::cout << "*Gap extension: " << gex << std::endl;
	std::cout << "*Scoring matrix: " << std::endl;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	std::cout << "   ";
	for(TSize col = 0; col<alphSize; ++col) std::cout << TAlphabet(col) << ',';
	std::cout << std::endl;
	for(TSize row = 0; row<alphSize; ++row) {
		for(TSize col = 0; col<alphSize; ++col) {
			if (col == 0) std::cout << TAlphabet(row) << ": ";
			std::cout << score(msaOpt.sc, TAlphabet(row), TAlphabet(col));
			if (col < alphSize - 1) std::cout << ',';
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	// Print the alignment information
	TSize numGapEx = 0;
	TSize numGap = 0;
	TSize numPairs = 0;
	TSize alignLen = 0;
	String<TSize> pairCount;
	String<char> mat;
	if (convertAlignment(g, mat)) {
		TScoreValue alignScore = alignmentEvaluation(g, msaOpt.sc, numGapEx, numGap, numPairs, pairCount, alignLen);
		std::cout << "Alignment Score: " << alignScore << std::endl;
		std::cout << "Alignment Length: " << alignLen << std::endl;
		std::cout << "#Match-Mismatch pairs: " << numPairs << std::endl;
		std::cout << "Score contribution by match-mismatch pairs: " << (alignScore - (((TScoreValue) numGap * gop) + ((TScoreValue) numGapEx * gex))) << std::endl;
		std::cout << "#Gap extensions: " << numGapEx << std::endl;
		std::cout << "Score contribution by gap extensions: " << ((TScoreValue) numGapEx * gex) << std::endl;
		std::cout << "#Gap openings: " << numGap << std::endl;
		std::cout << "Score contribution by gap openings: " << ((TScoreValue) numGap * gop) << std::endl;
		std::cout << std::endl;
		std::cout << "#Pairs: " << std::endl;
		std::cout << "   ";
		for(TSize col = 0; col<alphSize; ++col) std::cout << TAlphabet(col) << ',';
		std::cout << std::endl;
		for(TSize row = 0; row<alphSize; ++row) {
			for(TSize col = 0; col<alphSize; ++col) {
				if (col == 0) std::cout << TAlphabet(row) << ": ";
				std::cout << value(pairCount, row * alphSize + col);
				if (col < alphSize - 1) std::cout << ',';
			}
			std::cout << std::endl;
		}
	} else {
		std::cout << "No valid alignment!" << std::endl;
	}
}


//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TStringSet1, typename TNames, typename TAlphabet, typename TScore>
inline void
globalMsaAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gAlign, 
				   TStringSet1& sequenceSet,
				   TNames& sequenceNames,
				   MsaOptions<TAlphabet, TScore> const& msaOpt)
{
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	typedef Graph<Alignment<TStringSet, TSize> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	typedef double TDistanceValue;
	
	
	// Initialize alignment object
	clear(gAlign);
	assignStringSet(gAlign, sequenceSet);

#ifdef SEQAN_PROFILE
	std::cout << "Scoring parameters:" << std::endl;
	std::cout << "*Gap opening: " << msaOpt.sc.data_gap_open << std::endl;
	std::cout << "*Gap extension: " << msaOpt.sc.data_gap_extend << std::endl;
	//std::cout << "*Scoring matrix: " << std::endl;
	//TSize alphSize = ValueSize<TAlphabet>::VALUE;
	//std::cout << "   ";
	//for(TSize col = 0; col<alphSize; ++col) std::cout << TAlphabet(col) << ',';
	//std::cout << std::endl;
	//for(TSize row = 0; row<alphSize; ++row) {
	//	for(TSize col = 0; col<alphSize; ++col) {
	//		if (col == 0) std::cout << TAlphabet(row) << ": ";
	//		std::cout << score(msaOpt.sc, TAlphabet(row), TAlphabet(col));
	//		if (col < alphSize - 1) std::cout << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
#endif

	// Some alignment constants
	TStringSet& seqSet = stringSet(gAlign);
	TSize nSeq = length(seqSet);
	TSize threshold = 30;

	// Select all possible pairs for global and local alignments
	String<TSize> pList;
	selectPairs(seqSet, pList);

	// Set-up a distance matrix
	typedef String<TDistanceValue> TDistanceMatrix;
	TDistanceMatrix distanceMatrix;

	// Containers for segment matches and corresponding scores 
	typedef String<Fragment<> > TFragmentString;
	TFragmentString matches;
	typedef String<TScoreValue> TScoreValues;
	TScoreValues scores;

	// Include segment matches from subalignments
	if (!empty(msaOpt.alnfiles)) {
#ifdef SEQAN_PROFILE
		std::cout << "Parsing external alignment files:" << std::endl;
#endif
		typedef typename Iterator<String<std::string>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.alnfiles, Standard() );
		TIter begItEnd = end(msaOpt.alnfiles, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
#ifdef SEQAN_PROFILE
			std::cout << "*External file " << (*begIt) << std::endl;
#endif
			std::ifstream strm_lib;
			strm_lib.open((*begIt).c_str(), ::std::ios_base::in | ::std::ios_base::binary);
			read(strm_lib, matches, scores, sequenceNames, FastaAlign());
			strm_lib.close();
		}
#ifdef SEQAN_PROFILE
		std::cout << "External segment matches done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}

	// Include computed segment matches
	if (!empty(msaOpt.method)) {
#ifdef SEQAN_PROFILE
		std::cout << "Computing segment matches:" << std::endl;
#endif
		typedef typename Iterator<String<unsigned int>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.method, Standard() );
		TIter begItEnd = end(msaOpt.method, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
#ifdef SEQAN_PROFILE
			if (*begIt == 0) std::cout << "*Method: global" << std::endl;
			else if (*begIt == 1) std::cout << "*Method: local" << std::endl;
			else if (*begIt == 2) std::cout << "*Method: overlap" << std::endl;
			else if (*begIt == 3) std::cout << "*Method: lcs" << std::endl;
#endif
			if (*begIt == 0) appendSegmentMatches(seqSet, pList, msaOpt.sc, matches, scores, distanceMatrix, GlobalPairwise_Library() );
			else if (*begIt == 1) appendSegmentMatches(seqSet, pList, msaOpt.sc, matches, scores, LocalPairwise_Library() );
			else if (*begIt == 2) {
				Nothing noth;
				AlignConfig<true,true,true, true> ac;
				appendSegmentMatches(seqSet, pList, msaOpt.sc, matches, scores, noth, ac, GlobalPairwise_Library() );
			} 
			else if (*begIt == 3) appendSegmentMatches(seqSet, pList, matches, scores, Lcs_Library() );
			else {
#ifdef SEQAN_PROFILE
				std::cout << "*Unknown method!!!" << std::endl;
#endif
			}
#ifdef SEQAN_PROFILE
			std::cout << "*Done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
		}	
	}

	// Include a T-Coffee library
	if (!empty(msaOpt.libfiles)) {
#ifdef SEQAN_PROFILE
		std::cout << "Parsing a T-Coffee Library:" << std::endl;
#endif
		typedef typename Iterator<String<std::string>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.libfiles, Standard() );
		TIter begItEnd = end(msaOpt.libfiles, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
#ifdef SEQAN_PROFILE
			std::cout << "*T-Coffee library: " << (*begIt) << std::endl;
#endif
			std::ifstream strm_lib;
			strm_lib.open((*begIt).c_str(), std::ios_base::in | std::ios_base::binary);
			read(strm_lib, matches, scores, sequenceNames, TCoffeeLib());
			strm_lib.close();
		}
#ifdef SEQAN_PROFILE
		std::cout << "Parsing done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}

	// Include MUMmer segment matches
	if (!empty(msaOpt.mummerfiles)) {
#ifdef SEQAN_PROFILE
		std::cout << "Parsing MUMmer segment matches:" << std::endl;
#endif
		typedef typename Iterator<String<std::string>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.mummerfiles, Standard() );
		TIter begItEnd = end(msaOpt.mummerfiles, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
#ifdef SEQAN_PROFILE
			std::cout << "*MUMmer file: " << (*begIt) << std::endl;
#endif
			std::ifstream strm_lib;
			strm_lib.open((*begIt).c_str(), std::ios_base::in | std::ios_base::binary);
			read(strm_lib, matches, scores, seqSet, sequenceNames, MummerLib());		
			strm_lib.close();
		}
#ifdef SEQAN_PROFILE
		std::cout << "Parsing done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}

	// Include BLAST segment matches
	if (!empty(msaOpt.blastfiles)) {
#ifdef SEQAN_PROFILE
		std::cout << "Parsing BLAST segment matches:" << std::endl;
#endif
		typedef typename Iterator<String<std::string>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.blastfiles, Standard() );
		TIter begItEnd = end(msaOpt.blastfiles, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
#ifdef SEQAN_PROFILE
			std::cout << "*BLAST file: " << (*begIt) << std::endl;
#endif
			std::ifstream strm_lib;
			strm_lib.open((*begIt).c_str(), std::ios_base::in | std::ios_base::binary);
			read(strm_lib, matches, scores, sequenceNames, BlastLib());
			strm_lib.close();
		}
#ifdef SEQAN_PROFILE
		std::cout << "Parsing done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}

#ifdef SEQAN_PROFILE
	std::cout << "Total number of segment matches: " << length(matches) << std::endl;
#endif

	// Use these segment matches for the initial alignment graph
	TGraph g(seqSet);
	if (!msaOpt.rescore) buildAlignmentGraph(matches, scores, g, FractionalScore() );
	else buildAlignmentGraph(matches, scores, g, msaOpt.sc, ReScore() );
	clear(matches);
	clear(scores);
#ifdef SEQAN_PROFILE
	std::cout << "Alignment graph construction done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Guide tree
	Graph<Tree<TDistanceValue> > guideTree;
	if (!empty(msaOpt.treefile)) {
#ifdef SEQAN_PROFILE
		std::cout << "Guide Tree: " << msaOpt.treefile << std::endl;
#endif
		std::fstream strm_tree;
		strm_tree.open(msaOpt.treefile.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
		read(strm_tree, guideTree, sequenceNames, NewickFormat());	// Read newick tree
		strm_tree.close();
	} else {
#ifdef SEQAN_PROFILE
		if (msaOpt.build == 0) std::cout << "Guide Tree: Neighbor Joining" << std::endl;
		else if (msaOpt.build == 1) std::cout << "Guide Tree: UPGMA single linkage" << std::endl;
		else if (msaOpt.build == 2) std::cout << "Guide Tree: UPGMA complete linkage" << std::endl;
		else if (msaOpt.build == 3) std::cout << "Guide Tree: UPGMA average linkage" << std::endl;
		else if (msaOpt.build == 4) std::cout << "Guide Tree: UPGMA weighted average linkage" << std::endl;
#endif
		// Check if we have a valid distance matrix
		if (empty(distanceMatrix)) getDistanceMatrix(g, distanceMatrix, KmerDistance());
		if (msaOpt.build == 0) njTree(distanceMatrix, guideTree);
		else if (msaOpt.build == 1) upgmaTree(distanceMatrix, guideTree, UpgmaMin());
		else if (msaOpt.build == 2) upgmaTree(distanceMatrix, guideTree, UpgmaMax());
		else if (msaOpt.build == 3) upgmaTree(distanceMatrix, guideTree, UpgmaAvg());
		else if (msaOpt.build == 4) upgmaTree(distanceMatrix, guideTree, UpgmaWeightAvg());
	}
	clear(distanceMatrix);
#ifdef SEQAN_PROFILE
	std::cout << "Guide tree done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
		
	// Triplet extension
	if (nSeq < threshold) tripletLibraryExtension(g);
	else tripletLibraryExtension(g, guideTree, threshold / 2);
#ifdef SEQAN_PROFILE
	std::cout << "Triplet extension done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Progressive Alignment
	progressiveAlignment(g, guideTree, gAlign);
#ifdef SEQAN_PROFILE
	std::cout << "Progressive alignment done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	clear(guideTree);
	clear(g);
#ifdef SEQAN_PROFILE
	std::cout << "Clean-up done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
}


//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline void
globalMsaAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gAlign,
				   TScore const& scoreObject)
{
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Value<TString>::Type TAlphabet;
	TStringSet sequenceSet = stringSet(gAlign);
	String<String<char> > sequenceNames;
	fill(sequenceNames, length(sequenceSet), String<char>("tmpName"));
	MsaOptions<AminoAcid, TScore> msaOpt;
	msaOpt.sc = scoreObject;
	appendValue(msaOpt.method, 0);  // Global pairwise
	appendValue(msaOpt.method, 1);	// Local pairwise
	globalMsaAlignment(gAlign, sequenceSet, sequenceNames, msaOpt);
}

//////////////////////////////////////////////////////////////////////////////


template<typename TSource, typename TSpec, typename TScore>
inline void
globalMsaAlignment(Align<TSource, TSpec>& align,
				   TScore const& scoreObject)
{
	typedef StringSet<TSource, Dependent<> > TStringSet;
	TStringSet sequenceSet = stringSet(align);
	Graph<Alignment<TStringSet, void, WithoutEdgeId> > gAlign(sequenceSet);
	globalMsaAlignment(gAlign, scoreObject);
	
	// Pipe into Align data structure
	String<char> mat;
	convertAlignment(gAlign, mat);
	typedef Align<TSource, TSpec> TAlign;
	typedef typename Size<TAlign>::Type TSize;
	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow>::Type TRowIterator;
	clearGaps(align);
	TSize nseq = length(sequenceSet);
	String<TRowIterator> rowIter;
	resize(rowIter, nseq);
	for(TSize i = 0; i<nseq; ++i) value(rowIter, i) = begin(row(align, i));
	TSize lenMat = length(mat);
	TSize colLen = lenMat / nseq;
	TSize gapCount = 0;
	char gapChar = gapValue<char>();
	for(TSize alignRow = 0; alignRow < nseq; ++alignRow) {
		for(TSize pos = alignRow * colLen; pos < (alignRow + 1) * colLen; ++pos) {
			if (value(mat, pos) != gapChar) {
				if (gapCount) {
					insertGaps(value(rowIter, alignRow), gapCount);
					goFurther(value(rowIter, alignRow), gapCount);
					gapCount = 0;
				}
				goNext(value(rowIter,alignRow));
			} else ++gapCount;
		}
		if (gapCount) {
			insertGaps(value(rowIter, alignRow), gapCount);
			goFurther(value(rowIter, alignRow), gapCount);
			gapCount = 0;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// Just two testing functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TMatches>
void
_debugMatches(TStringSet& str, 
			  TMatches& matches)
{
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Size<TStringSet>::Type TSize;

	// Print all the matches
	std::cout << "The sequences:" << std::endl;
	for(TSize i = 0;i<length(str);++i) {
		std::cout << positionToId(str,i) << ':' << str[i] << std::endl;
	}
	std::cout << "The matches:" << std::endl;
	for(TSize i = 0;i<length(matches);++i) {
		TId tmp_id1 = sequenceId(matches[i],0);
		std::cout << tmp_id1 << ',' << fragmentBegin(matches[i],tmp_id1) << ',';
		for(TSize j = fragmentBegin(matches[i],tmp_id1); j < fragmentBegin(matches[i],tmp_id1) + fragmentLength(matches[i],tmp_id1); ++j) {
			std::cout << str[idToPosition(str, tmp_id1)][j];
		}
		TId tmp_id2 = sequenceId(matches[i],1);
		std::cout << ',' <<	tmp_id2 << ',' << fragmentBegin(matches[i],tmp_id2) << ',';
		for(TSize j = fragmentBegin(matches[i],tmp_id2); j < fragmentBegin(matches[i],tmp_id2) + fragmentLength(matches[i],tmp_id2); ++j) {
			std::cout << str[idToPosition(str, tmp_id2)][j];
		}
		std::cout << std::endl;

		SEQAN_TASSERT(sequenceId(matches[i],0) != sequenceId(matches[i],1))
		SEQAN_TASSERT(fragmentBegin(matches[i],tmp_id1) < length(str[idToPosition(str, tmp_id1)]))
		SEQAN_TASSERT(fragmentBegin(matches[i],tmp_id1) + fragmentLength(matches[i],tmp_id1) <= length(str[idToPosition(str, tmp_id1)]))
		SEQAN_TASSERT(fragmentBegin(matches[i],tmp_id2) < length(str[idToPosition(str, tmp_id2)]))
		SEQAN_TASSERT(fragmentBegin(matches[i],tmp_id2) + fragmentLength(matches[i],tmp_id2) <= length(str[idToPosition(str, tmp_id2)]))
		SEQAN_TASSERT(fragmentLength(matches[i],tmp_id2) == fragmentLength(matches[i],tmp_id1))
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph>
void
_debugRefinedMatches(TGraph& g)
{
	typedef typename Id<TGraph>::Type TId;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	
	std::cout << "Refined matches" << std::endl;
	TEdgeIterator it_tmp(g);
	for(;!atEnd(it_tmp);++it_tmp) {
		TId id1 = sequenceId(g,sourceVertex(it_tmp));
		TId id2 = sequenceId(g,targetVertex(it_tmp));
		std::cout << id1 << ',' << fragmentBegin(g,sourceVertex(it_tmp)) << ',';
		std::cout << label(g,sourceVertex(it_tmp));
		std::cout << ',' <<	id2 << ',' << fragmentBegin(g,targetVertex(it_tmp)) << ',';
		std::cout << label(g,targetVertex(it_tmp));
		std::cout << " (" << cargo(*it_tmp) << ")";
		std::cout << std::endl;	

		SEQAN_TASSERT(sequenceId(g,sourceVertex(it_tmp)) != sequenceId(g,targetVertex(it_tmp)))
		SEQAN_TASSERT(fragmentBegin(g,sourceVertex(it_tmp)) < length((stringSet(g))[idToPosition((stringSet(g)), id1)]))
		SEQAN_TASSERT(fragmentBegin(g,sourceVertex(it_tmp)) + fragmentLength(g,sourceVertex(it_tmp)) <= length((stringSet(g))[idToPosition((stringSet(g)), id1)]))
		SEQAN_TASSERT(fragmentBegin(g,targetVertex(it_tmp)) < length((stringSet(g))[idToPosition((stringSet(g)), id2)]))
		SEQAN_TASSERT(fragmentBegin(g,targetVertex(it_tmp)) + fragmentLength(g,targetVertex(it_tmp)) <= length((stringSet(g))[idToPosition((stringSet(g)), id2)]))
		SEQAN_TASSERT(fragmentLength(g,sourceVertex(it_tmp)) == fragmentLength(g,targetVertex(it_tmp)))

	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

