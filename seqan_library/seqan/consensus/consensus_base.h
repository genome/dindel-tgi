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
  $Id: graph_consensus_base.h 2103 2008-05-23 07:57:13Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_CONSENSUS_BASE_H
#define SEQAN_HEADER_CONSENSUS_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// Segment Match Generation tag
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Segment Match Generation.value.Overlap_Library:
	Segment matches from overlap alignments.
*/

struct Overlap_Library_;
typedef Tag<Overlap_Library_> const Overlap_Library;



//////////////////////////////////////////////////////////////////////////////
// Consensus tag
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Consensus Calling:
..summary:A tag that specifies how to call the consensus.
*/


/**
.Tag.Consensus Calling.value.Majority_Vote:
	A consensus based on the most common character.
*/

struct Majority_Vote_;
typedef Tag<Majority_Vote_> const Majority_Vote;

/**
.Tag.Consensus Calling.value.Bayesian:
	A consensus based on bayesian probability.
*/

struct Bayesian_;
typedef Tag<Bayesian_> const Bayesian;



//////////////////////////////////////////////////////////////////////////////
// Read alignment and Consensus Generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////

struct ConsensusOptions {
public:
	// Method
	// 0: graph-based multiple sequence alignment
	// 1: realign
	int method;

	// ReAlign Method
	// 0: Needleman-Wunsch
	// 1: Gotoh
	int rmethod;

	// Bandwidth of overlap alignment
	int bandwidth;

	// Number of computed overlaps per read (at the beginning and end of a read)
	int overlaps;

	// Minimum match length of a computed overlap
	int matchlength;

	// Minimum quality (in percent identity) of a computed overlap
	int quality;

	// Window size, only relevant for insert sequencing
	// If window == 0, no insert sequencing is assumed
	int window;
	
	// Consensus calling
	// 0: majority
	// 1: bayesian
	int consensus;

	// Output
	// 0: seqan style
	// 1: afg output format
	int output;

	// Multi-read alignment
	bool noalign;

	// Offset all reads, so the first read starts at position 0
	bool moveToFront;

	// Include reference genome
	bool include;

	// Scoring object for overlap alignments
	Score<int> sc;

	// Various input and output files
	std::string readsfile;				// File of reads in FASTA format
	std::string afgfile;				// AMOS afg file input
	std::string outfile;				// Output file name
	
	// Initialization
	ConsensusOptions() 
	{
		sc = Score<int>(2,-6,-4,-9);
	}
};


//////////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TStrSpec, typename TPosPair, typename TStringSpec, typename TSpec, typename TConfig, typename TId>
inline void 
_loadContigReads(StringSet<TValue, Owner<TStrSpec> >& strSet,
				 String<TPosPair, TStringSpec>& startEndPos,
				 FragmentStore<TSpec, TConfig> const& fragStore,
				 TId const contigId)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadPos TReadPos;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

	// Sort aligned reads according to contig id
	sortAlignedReads(fragStore.alignedReadStore, SortContigId());
	resize(strSet, length(fragStore.alignedReadStore));

	// Retrieve all reads, limit them to the clear range and if required reverse complement them
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TSize numRead = 0;
	TReadPos begClr = 0;
	TReadPos endClr = 0;
	TSize lenRead = 0;
	TSize offset = 0;
	for(;alignIt != alignItEnd; ++alignIt) {
		offset = _min(alignIt->beginPos, alignIt->endPos);
		getClrRange(fragStore, *alignIt, begClr, endClr);
		strSet[numRead] = infix(fragStore.readSeqStore[alignIt->readId], begClr, endClr);
		lenRead = endClr - begClr;
		if (alignIt->beginPos < alignIt->endPos) appendValue(startEndPos, TPosPair(offset, offset + lenRead), Generous());
		else {
			reverseComplementInPlace(strSet[numRead]);
			appendValue(startEndPos, TPosPair(offset + lenRead, offset), Generous());
		}
		++numRead;
	}
	resize(strSet, numRead, Exact());
}

//////////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TSize, typename TConfigOptions>
inline void
consensusAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gOut,
				   String<Pair<TSize, TSize> >& begEndPos,
				   TConfigOptions const& consOpt) 
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TOutGraph;
	typedef typename Id<TOutGraph>::Type TId;

	// Initialization
	TStringSet& seqSet = stringSet(gOut);

	// Select all overlapping reads and record the diagonals of the band
	String<Pair<TId, TId> > pList;
	String<Pair<int, int> > diagList;
	if (consOpt.window == 0) selectPairs(seqSet, begEndPos, consOpt.bandwidth, pList, diagList);
	else selectPairsIndel(seqSet, begEndPos, consOpt.window, pList, diagList);

	// Estimate the number of overlaps we want to compute
#ifdef SEQAN_PROFILE
	std::cout << "Pair selection done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Set-up a sparse distance matrix
	Graph<Undirected<double> > pairGraph;
	
	// Containers for segment matches and corresponding scores 
	typedef String<Fragment<> > TFragmentString;
	TFragmentString matches;
	typedef String<int> TScoreValues;
	TScoreValues scores;

	// Compute segment matches from global pairwise alignments
	appendSegmentMatches(seqSet, pList, diagList, begEndPos, consOpt.sc, consOpt.matchlength, consOpt.quality, consOpt.overlaps, matches, scores, pairGraph, Overlap_Library() );
	clear(pList);
	clear(diagList);
#ifdef SEQAN_PROFILE
	std::cout << "Overlap done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Use these segment matches for the initial alignment graph
	typedef Graph<Alignment<TStringSet, TSize> > TGraph;
	TGraph g(seqSet);
	buildAlignmentGraph(matches, scores, g, consOpt.sc, ReScore() );
	clear(matches);
	clear(scores);
#ifdef SEQAN_PROFILE
	std::cout << "Construction of Alignment Graph done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Guide Tree
	Graph<Tree<double> > guideTree;
	upgmaTree(pairGraph, guideTree);
#ifdef SEQAN_PROFILE
	std::cout << "Guide tree done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	clear(pairGraph);

	// Triplet library extension
	if ( ((2 * numEdges(g)) / numVertices(g) ) < 50 ) graphBasedTripletLibraryExtension(g);
	else reducedTripletLibraryExtension(g);
#ifdef SEQAN_PROFILE
	std::cout << "Triplet done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Perform a progressive alignment
	progressiveAlignment(g, guideTree, gOut);
	clear(g);
	clear(guideTree);
#ifdef SEQAN_PROFILE
	std::cout << "Progressive alignment done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
}

//////////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TSize>
inline void
consensusAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gOut,
				   String<Pair<TSize, TSize> >& begEndPos) 
{
	ConsensusOptions consOpt;
	consensusAlignment(gOut, begEndPos, consOpt);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TFragSpec, typename TConfig, typename TStringSet, typename TCargo, typename TSpec, typename TContigId>
inline void
updateContig(FragmentStore<TFragSpec, TConfig>& fragStore,
			 Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
			 TContigId contigId)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef std::map<TSize, TSize> TComponentLength;
	typedef char TValue;

	// Initialization
	TStringSet& strSet = stringSet(g);
	TSize nseq = length(strSet);
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';
	TSize maxCoverage = 0;
	TSize len = 0;
	String<TValue> mat;

	// Store for each read the begin position, the end position and the row in the alignment matrix
	String<TSize> readBegEndRowPos;
	resize(readBegEndRowPos, 3*nseq);

	// Strongly Connected Components, topological sort, and length of each component
	String<TSize> component;
	String<TSize> order;
	TComponentLength compLength;
	if (convertAlignment(g, component, order, compLength)) {
		TSize numOfComponents = length(order);
		
		// Assign to each sequence the start and end (in terms of component ranks)
		typedef String<std::pair<TSize, TSize> > TComponentToRank;
		TComponentToRank compToRank;
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) 
			appendValue(compToRank, std::make_pair(order[compIndex], compIndex), Generous());
		::std::sort(begin(compToRank, Standard()), end(compToRank, Standard()));

		typedef Pair<TSize, TSize> TRankPair;
		typedef String<TRankPair> TSequenceToRanks;
		TSequenceToRanks seqToRank;
		resize(seqToRank, nseq);
		typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
		TVertexIterator itVertex(g);
		for(;!atEnd(itVertex);++itVertex) {
			TVertexDescriptor vert = value(itVertex);
			TSize seq = idToPosition(strSet, sequenceId(g, vert));
			if (fragmentBegin(g, vert) == 0) 
				seqToRank[seq].i1 = ::std::lower_bound(begin(compToRank, Standard()), end(compToRank, Standard()), ::std::make_pair((TSize) component[vert], (TSize) 0))->second;
			if (fragmentBegin(g, vert) + fragmentLength(g, vert) == length(strSet[seq]))
				seqToRank[seq].i2 = ::std::lower_bound(begin(compToRank, Standard()), end(compToRank, Standard()), ::std::make_pair((TSize) component[vert], (TSize) 0))->second;
		}
		clear(compToRank);

		// Assign the sequences to rows
		String<TSize> seqToRow;
		resize(seqToRow, nseq);
		maxCoverage = 0;
		typedef String<bool> TLeftOver;
		typedef typename Iterator<TLeftOver, Standard>::Type TLeftOverIter;
		TLeftOver leftOver;
		fill(leftOver, nseq, true);
		typedef String<std::pair<TSize, TSize> > TSeqToBegin;
		typedef typename Iterator<TSeqToBegin, Standard>::Type TSeqToBeginIter;
		TSeqToBegin seqToBegin;
		TSize finishedSeq = 0;
		while(finishedSeq < nseq) {
			TLeftOverIter itL = begin(leftOver, Standard());
			TLeftOverIter itLEnd = end(leftOver, Standard());
			for(TSize pos = 0; itL != itLEnd; ++itL, ++pos) 
				if (*itL) appendValue(seqToBegin, std::make_pair((seqToRank[pos]).i1, pos), Generous());
			::std::sort(begin(seqToBegin, Standard()), end(seqToBegin, Standard()));
			
			TSize endPos = 0;
			TSeqToBeginIter itSB = begin(seqToBegin, Standard());
			TSeqToBeginIter itSBEnd = end(seqToBegin, Standard());
			for(;itSB != itSBEnd;++itSB) {
				if (endPos <= (*itSB).first) {
					TSize currentSeq = (*itSB).second;
					seqToRow[currentSeq] = maxCoverage;
					endPos = (seqToRank[currentSeq]).i2 + 2;
					leftOver[currentSeq] = false;
					++finishedSeq;
				}	
			}
			clear(seqToBegin);
			++maxCoverage;
		}
		clear(leftOver);

		// Create the matrix
		len = 0;
		String<TSize> compOffset;
		resize(compOffset, numOfComponents);
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			compOffset[order[compIndex]] = len;
			len+=compLength[order[compIndex]];
		}
		fill(mat, len * maxCoverage, gapChar);

		// Fill in the segments
		typedef typename Infix<TString>::Type TInfix;
		typedef typename Iterator<TInfix, Standard>::Type TInfixIter;
		typedef typename TGraph::TPosToVertexMap TPosToVertexMap;
		for(typename TPosToVertexMap::const_iterator it = g.data_pvMap.begin();it != g.data_pvMap.end(); ++it) {
			TInfix str = label(g,it->second);
			TSize c = property(component, it->second);
			TSize row = seqToRow[idToPosition(strSet, it->first.first)];
			//if (row == 0) {
			//	std::cout << sequenceId(g, it->second) << ':' << str << ',' << strSet[sequenceId(g, it->second)] << std::endl;
			//	std::cout << getProperty(component, it->second) << ',' << order[compIndex] << std::endl;
			//	std::cout << (seqToRank[sequenceId(g, it->second)]).i1 << ',' << (seqToRank[sequenceId(g, it->second)]).i2 << std::endl;
			//}
			TInfixIter sIt = begin(str, Standard());
			TInfixIter sItEnd = end(str, Standard());
			TSize i = compOffset[c];
			for(TSize pCol = i;sIt!=sItEnd;++sIt, ++pCol, ++i) 
				mat[row * len + pCol] = *sIt;
		}
		String<bool> active;
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			TSize offset = compOffset[order[compIndex]];
			TSize currentCompLength = compLength[order[compIndex]];

			clear(active);
			fill(active, maxCoverage, false);

			// Find the empty rows
			for(TSize i=0;i<nseq; ++i) {
				if (((seqToRank[i]).i1 <= compIndex) && ((seqToRank[i]).i2 >= compIndex)) 
					active[(seqToRow[i])] = true;
			}
			
			// Substitute false gaps with special gap character
			for(TSize i = 0; i < maxCoverage; ++i) {
				if (!(active[i])) {
					for(TSize pCol = offset;pCol < offset + currentCompLength;++pCol) 
						mat[i * len + pCol] = specialGap;
				}
			}
		}

		// Get the new begin and end positions
		for(TSize i=0;i<nseq; ++i) {
			TVertexDescriptor lastVertex = findVertex(const_cast<TGraph&>(g), positionToId(strSet, i), length(strSet[i]) - 1);
			TSize readBegin = compOffset[getProperty(component, findVertex(const_cast<TGraph&>(g), positionToId(strSet, i), 0))];
			TSize readEnd = compOffset[getProperty(component, lastVertex)] + fragmentLength(const_cast<TGraph&>(g), lastVertex);
			readBegEndRowPos[3*i] = readBegin;
			readBegEndRowPos[3*i+1] = readEnd;
			readBegEndRowPos[3*i+2] = seqToRow[i];
		}

	}
	clear(component);
	clear(order);
	compLength.clear();

	
	//// Debug code
	//for(TSize row = 0; row<maxCoverage; ++row) {
	//	for(TSize col = 0; col<len; ++col) {
	//		std::cout << mat[row * len + col];			
	//	}
	//	std::cout << std::endl;
	//}

	// Update the contig
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigStoreElement;
	TContigStoreElement& contigEl = fragStore.contigStore[contigId];
	clear(contigEl.gaps);
	clear(contigEl.seq);


	typedef typename Value<TString>::Type TAlphabet;
	String<TAlphabet> consensus;
	String<TValue> gappedCons;
	String<TSize> coverage;
	consensusCalling(mat, consensus, gappedCons, coverage, maxCoverage, Majority_Vote());

	// Create the gap anchors
	typedef typename Iterator<String<TValue>, Standard>::Type TStringIter;
	TStringIter seqIt = begin(gappedCons, Standard());
	TStringIter seqItEnd = end(gappedCons, Standard());
	typedef typename TFragmentStore::TReadPos TReadPos;
	typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
	TReadPos ungappedPos = 0;
	TReadPos gappedPos = 0;
	bool gapOpen = false;
	for(;seqIt != seqItEnd; goNext(seqIt), ++gappedPos) {
		if (value(seqIt) == gapChar) gapOpen = true;				
		else {
			if (gapOpen) {
				appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());
				gapOpen = false;
			}
			Dna5Q letter = value(seqIt);
			assignQualityValue(letter, 'D');
			appendValue(contigEl.seq, letter);
			++ungappedPos;
		}
	}
	if (gapOpen) 
		appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());


	// Update all aligned reads
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	sortAlignedReads(fragStore.alignedReadStore, SortContigId());
	TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	for(TSize i = 0;alignIt != alignItEnd; ++alignIt, ++i) {
		TSize lenRead = length(fragStore.readSeqStore[alignIt->readId]);
		TReadPos begClr = 0;
		TReadPos endClr = 0;
		getClrRange(fragStore, *alignIt, begClr, endClr);
		clear(alignIt->gaps);
		ungappedPos = begClr;
		if (alignIt->beginPos > alignIt->endPos) ungappedPos = lenRead - endClr;
		if (ungappedPos != 0) appendValue(alignIt->gaps, TContigGapAnchor(ungappedPos, 0));
		gappedPos = 0;
		gapOpen = false;
		for(TSize column = readBegEndRowPos[3*i]; column<readBegEndRowPos[3*i + 1]; ++column, ++gappedPos) {
			if (mat[readBegEndRowPos[3*i + 2] * len + column] == gapChar) gapOpen = true;				
			else {
				if (gapOpen) {
					appendValue(alignIt->gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());
					gapOpen = false;
				}
				++ungappedPos;
			}
		}
		if (gapOpen) appendValue(alignIt->gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous());
		if (alignIt->beginPos < alignIt->endPos) {
			if (endClr != lenRead) 
				appendValue(alignIt->gaps, TContigGapAnchor(lenRead, lenRead + (gappedPos - ungappedPos) - (lenRead - endClr)), Generous());
		} else {
			if (begClr != 0) 
				appendValue(alignIt->gaps, TContigGapAnchor(lenRead, lenRead + (gappedPos - ungappedPos) - begClr), Generous());
		}

		// Set new begin and end position
		if (alignIt->beginPos < alignIt->endPos) {
			alignIt->beginPos = readBegEndRowPos[3*i];
			alignIt->endPos = readBegEndRowPos[3*i+1];
		} else {
			alignIt->beginPos = readBegEndRowPos[3*i+1];
			alignIt->endPos = readBegEndRowPos[3*i];
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TCounters, typename TCoverage, typename TSize, typename TAlphabet>
inline void
__countLetters(String<TValue, TSpec> const& mat,
			   TCounters& counterValues,
			   TCoverage& coverage,
			   TSize alignDepth,
			   TAlphabet)
{
	SEQAN_CHECKPOINT
	typedef String<TValue, TSpec> TMatrix;
	typedef typename Iterator<TMatrix>::Type TMatIter;
	typedef typename Iterator<TCoverage>::Type TCovIter;

	// Initialization
	TSize len = length(mat) / alignDepth;
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';
	TSize alphabetSize = ValueSize<TAlphabet>::VALUE;
	clear(coverage);
	fill(coverage, len, 0);


	// Set-up counter values
	typedef typename Value<TCounters>::Type TCounter;
	typedef typename Iterator<TCounters>::Type TCounterIt;
	clear(counterValues);
	resize(counterValues, len);
	for(TSize i=0;i<len; ++i) {
		TCounter counter;
		fill(counter, alphabetSize + 1, 0);
		value(counterValues, i) = counter;
	}

	// Count all 
	TMatIter matIt = begin(mat);
	TMatIter matItEnd = end(mat);
	TCounterIt countIt = begin(counterValues);
	TCovIter covIt = begin(coverage);
	TSize pos = 0;
	for(; matIt != matItEnd; goNext(matIt), goNext(countIt), goNext(covIt), ++pos) {
		if (pos % len == 0) {
			countIt = begin(counterValues);
			covIt = begin(coverage);
		}
		TValue c = value(matIt);
		if (c == specialGap) continue;
		else {
			++value(covIt);
			if (c == gapChar) ++value(value(countIt), alphabetSize);
			else ++value(value(countIt), ordValue(TAlphabet(c)));
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TAlphabet, typename TSpec2, typename TGappedConsensus, typename TCoverage, typename TSize>
inline void
consensusCalling(String<TValue, TSpec> const& mat,
				 String<TAlphabet, TSpec2>& consensus,
				 TGappedConsensus& gappedConsensus,
				 TCoverage& coverage,
				 TSize alignDepth,
				 Bayesian)
{
	SEQAN_CHECKPOINT
	typedef double TProbability;
	typedef String<TProbability> TProbabilityDistribution;
	typedef String<TProbabilityDistribution> TPositionalPrDist;
	typedef typename Iterator<TPositionalPrDist>::Type TPosPrDistIter;

	// Initialization
	TSize len = length(mat) / alignDepth;
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';
	TSize alphabetSize = ValueSize<TAlphabet>::VALUE;
	TProbabilityDistribution backroundDist;
	fill(backroundDist, alphabetSize + 1, (1.0 / (TProbability) (alphabetSize + 1)));
	
	// Set-up the counters
	typedef String<TSize> TCounter;
	typedef String<TCounter> TCounters;
	TCounters counterValues;
	__countLetters(mat, counterValues, coverage, alignDepth, TAlphabet() );

	// Get an initial consensus
	typedef typename Iterator<TCounters>::Type TCounterIt;
	TCounterIt countIt = begin(counterValues);
	TCounterIt countItEnd = end(counterValues);
	TPositionalPrDist posPrDist;
	for(;countIt != countItEnd; goNext(countIt)) {
		TSize max = 0;
		TValue c = TValue();
		typedef typename Iterator<TCounter>::Type TCIt;
		TCIt cIt = begin(value(countIt));
		TCIt cItEnd = end(value(countIt));
		TSize pos = 0;
		for(;cIt != cItEnd; goNext(cIt), ++pos) {
			if (value(cIt) > max) {
				max = value(cIt);
				if (pos == alphabetSize) c = gapChar;
				else c = TAlphabet(pos);
			}
		}
		TProbabilityDistribution prDist;
		fill(prDist, alphabetSize + 1, 0);
		if (c == gapChar) value(prDist, alphabetSize) = 1;
		else value(prDist, ordValue(TAlphabet(c))) = 1;
		appendValue(posPrDist, prDist);
	}

	bool run = false;
	TProbabilityDistribution pI;
	TProbabilityDistribution pIJ;
	TProbabilityDistribution pIOld;
	TProbabilityDistribution pIJOld;
	std::cout << "Bayesian Consensus";
	while((run) || (empty(pIOld))) {
		// Store the values from the last iteration
		pIOld = pI;
		pIJOld = pIJ;

		// Count all letters in the consensus
		TProbabilityDistribution nI;
		fill(nI, alphabetSize + 1, 0);
		TPosPrDistIter itPosPrDist = begin(posPrDist);
		TPosPrDistIter itPosPrDistEnd = end(posPrDist);
		for(;itPosPrDist!=itPosPrDistEnd; ++itPosPrDist) {
			for(TSize i = 0; i<(alphabetSize + 1); ++i) {
				value(nI, i) += value(value(itPosPrDist), i);
			}
		}
	
		// Composition probabilities
		clear(pI);
		resize(pI, alphabetSize + 1);
		for(TSize i = 0; i<length(pI); ++i) {
			value(pI, i) = (TProbability) value(nI, i) / (TProbability) length(posPrDist);
		}

		// Count all letters that agree / disagree with the consensus
		TProbabilityDistribution nIJ;
		fill(nIJ, (alphabetSize + 1) * (alphabetSize + 1), 0);
		typedef String<TValue, TSpec> TMatrix;
		typedef typename Iterator<TMatrix>::Type TMatIter;
		TMatIter matIt = begin(mat);
		TMatIter matItEnd = end(mat);
		itPosPrDist = begin(posPrDist);
		TSize pos = 0;
		for(; matIt != matItEnd; goNext(matIt), goNext(itPosPrDist), ++pos) {
			if (pos % len == 0) {
				itPosPrDist = begin(posPrDist);
			}
			TValue c = value(matIt);
			if (c == specialGap) continue;
			else {
				TSize fragJ = alphabetSize;
				if (c != gapChar) fragJ = ordValue(TAlphabet(c));
				for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) {
					value(nIJ, consI * (alphabetSize + 1) + fragJ) += 1.0 * value(value(itPosPrDist), consI);
				}
			}
		}

		// Sequencing error probabilities
		clear(pIJ);
		resize(pIJ, (alphabetSize + 1) * (alphabetSize + 1));
		TProbability sumIJ = 0;
		for(TSize diag = 0; diag<(alphabetSize + 1); ++diag) sumIJ += value(nIJ, diag * (alphabetSize + 1) + diag);
		for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) {
			for(TSize fragJ = 0; fragJ<(alphabetSize + 1); ++fragJ) {
				value(pIJ, consI * (alphabetSize + 1) + fragJ) = value(nIJ, consI * (alphabetSize + 1) + fragJ) / sumIJ;
			}
		}
	
		//// Debug Code
		//std::cout << "A " << value(pI, 0) << std::endl;
		//std::cout << "C " << value(pI, 1) << std::endl;
		//std::cout << "G " << value(pI, 2) << std::endl;
		//std::cout << "T " << value(pI, 3) << std::endl;
		//std::cout << "- " << value(pI, 4) << std::endl;
		//std::cout << "AA " << value(pIJ, 0) << std::endl;
		//std::cout << "AC " << value(pIJ, 1) << std::endl;
		//std::cout << "AG " << value(pIJ, 2) << std::endl;
		//std::cout << "AT " << value(pIJ, 3) << std::endl;
		//std::cout << "A- " << value(pIJ, 4) << std::endl;
		//std::cout << "CA " << value(pIJ, 5) << std::endl;
		//std::cout << "CC " << value(pIJ, 6) << std::endl;
		//std::cout << "CG " << value(pIJ, 7) << std::endl;
		//std::cout << "CT " << value(pIJ, 8) << std::endl;
		//std::cout << "C- " << value(pIJ, 9) << std::endl;
		//std::cout << "GA " << value(pIJ, 10) << std::endl;
		//std::cout << "GC " << value(pIJ, 11) << std::endl;
		//std::cout << "GG " << value(pIJ, 12) << std::endl;
		//std::cout << "GT " << value(pIJ, 13) << std::endl;
		//std::cout << "G- " << value(pIJ, 14) << std::endl;
		//std::cout << "TA " << value(pIJ, 15) << std::endl;
		//std::cout << "TC " << value(pIJ, 16) << std::endl;
		//std::cout << "TG " << value(pIJ, 17) << std::endl;
		//std::cout << "TT " << value(pIJ, 18) << std::endl;
		//std::cout << "T- " << value(pIJ, 19) << std::endl;
		//std::cout << "-A " << value(pIJ, 20) << std::endl;
		//std::cout << "-C " << value(pIJ, 21) << std::endl;
		//std::cout << "-G " << value(pIJ, 22) << std::endl;
		//std::cout << "-T " << value(pIJ, 23) << std::endl;
		//std::cout << "-- " << value(pIJ, 24) << std::endl;

		// Recompute positional probability distribution
		itPosPrDist = begin(posPrDist);
		TSize col = 0;
		for(;itPosPrDist!=itPosPrDistEnd; goNext(itPosPrDist), ++col) {
			TProbabilityDistribution prDist;
			resize(prDist, alphabetSize + 1);
			for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) {
				TProbability numerator = value(pI, consI);
				TProbability denominator = 0;
				for(TSize allI = 0; allI<(alphabetSize + 1); ++allI) {
					TProbability denominatorSub = value(pI, allI);
					for(TSize row = 0; row < alignDepth; ++row) {
						TValue c = value(mat, row * len + col);
						if (c == specialGap) continue;
						TSize fragJ = alphabetSize;
						if (c != gapChar) fragJ = ordValue(TAlphabet(c));
						if (allI == consI) {
							numerator *= value(pIJ, allI * (alphabetSize + 1) + fragJ); 
						}
						denominatorSub *= value(pIJ, allI * (alphabetSize + 1) + fragJ); 
					}
					denominator += denominatorSub;
				}
				value(prDist, consI) = numerator / denominator;
			}
			value(itPosPrDist) = prDist;
		}	

		// Check termination criterion
		TProbability eps = 0.00001;
		typedef typename Iterator<TProbabilityDistribution>::Type TProbIter;
		TProbIter pIter = begin(pIOld);
		TProbIter pIterCompare = begin(pI);
		TProbIter pIterEnd = end(pIOld);
		run = false;
		for(;pIter != pIterEnd; goNext(pIter), goNext(pIterCompare)) {
			if (value(pIter) > value(pIterCompare)) {
				if (value(pIter) - value(pIterCompare) > eps) {
					run = true;
					break;
				}
			} else {
				if (value(pIterCompare) - value(pIter) > eps) {
					run = true;
					break;
				}
			}
		}
		if (!run) {
			pIter = begin(pIJOld);
			pIterCompare = begin(pIJ);
			pIterEnd = end(pIJOld);
			for(;pIter != pIterEnd; goNext(pIter), goNext(pIterCompare)) {
				if (value(pIter) > value(pIterCompare)) {
					if (value(pIter) - value(pIterCompare) > eps) {
						run = true;
						break;
					}
				} else {
					if (value(pIterCompare) - value(pIter) > eps) {
						run = true;
						break;
					}
				}
			}
		}
		std::cout << '.';
	}
	std::cout << std::endl;
	
	// Compute the most likely consensus
	TPosPrDistIter itPosPrDist = begin(posPrDist);
	TPosPrDistIter itPosPrDistEnd = end(posPrDist);
	clear(consensus);
	clear(gappedConsensus);
	for(;itPosPrDist!=itPosPrDistEnd; goNext(itPosPrDist)) {
		TProbability max = 0;
		TSize ind = 0;
		for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) {
			if (value(value(itPosPrDist), consI) > max) {
				max = value(value(itPosPrDist), consI);
				ind = consI;
			}
		}
		if (ind == alphabetSize) appendValue(gappedConsensus, gapChar);
		else {
			appendValue(consensus, TAlphabet(ind));
			appendValue(gappedConsensus, TAlphabet(ind));
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TStringSet, typename TScore, typename TBegEndRow, typename TSize>
inline TSize
fixDisruptedReads(String<TValue, TSpec>& mat,
				  TStringSet const& str,
				  TScore& scType,
				  TBegEndRow& begER,
				  TSize alignDepth)
{
	typedef typename Value< typename Value<TStringSet>::Type >::Type TAlphabet;

	// The limit expand threshold (the lower the more reads will be realigned)
	TSize limitExpand = 60;
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';
	TSize len = length(mat) / alignDepth;

	// Find all bad reads
	typedef typename Iterator<TBegEndRow>::Type TBEIter;
	typedef typename Iterator<TStringSet const>::Type TStrIter;
	TStrIter strIt = begin(str);
	TStrIter strItEnd = end(str);
	TBEIter itBE = begin(begER);
	TBEIter itBEEnd = end(begER);
	typedef String<TSize> TBadReads;
	TBadReads badReads;
	TSize counter = 0;
	for(;strIt != strItEnd; goNext(strIt), goNext(itBE), ++counter) {
		TSize b = value(itBE).i1;
		TSize e = value(itBE).i2;
		TSize diff = e - b;
		if (b > e) diff = b - e;
		if (length(value(strIt)) + limitExpand <= diff) {
			appendValue(badReads, counter);
		}
	}

	// Call the consensus
	String<unsigned int> coverage;
	String<char> gappedConsensus;
	String<Dna> consensusSequence;
	consensusCalling(mat, consensusSequence, gappedConsensus, coverage, alignDepth, Majority_Vote() );
	
	// Realign to the consensus
	typedef typename Iterator<TBadReads>::Type TBadIter;
	TBadIter badIt = begin(badReads);
	TBadIter badItEnd = end(badReads);
	for(;badIt != badItEnd; goNext(badIt)) {
		TSize b = value(begER, value(badIt)).i1;
		TSize e = value(begER, value(badIt)).i2;
		if (b > e) { TSize tmp = b; b = e; e = tmp; }
		String<char> cons = infix(gappedConsensus, b, e);
		String<TAlphabet> noGapCons;
		typedef typename Iterator<String<char> >::Type TCharIter;
		TCharIter charIt = begin(cons);
		TCharIter charItEnd = end(cons);
		for(;charIt != charItEnd; goNext(charIt)) {
			if (value(charIt) != gapChar) appendValue(noGapCons, value(charIt));
		}
		
		// Make a pairwise string-set
		TStringSet pairSet;
		//std::cout << noGapCons << std::endl;
		//std::cout << value(str, value(badIt)) << std::endl;
		appendValue(pairSet, noGapCons);
		appendValue(pairSet, value(str, value(badIt)));
		
		// Re-align
		// Maybe LCS is better???
		Graph<Alignment<TStringSet, TSize> > tmp(pairSet);
		globalAlignment(tmp, pairSet, scType, AlignConfig<true,true,true,true>(), Gotoh() );
		//std::cout << tmp << std::endl;		

		// Walk through all 3 sequences in parallel
		String<char> localAlign;
		convertAlignment(tmp, localAlign);
		String<char> consAlign = infix(localAlign, 0, length(localAlign) / 2);
		String<char> seqAlign = infix(localAlign, length(localAlign) / 2, length(localAlign));
		
		//std::cout << cons << std::endl;
		//std::cout << consAlign << std::endl;
		//std::cout << seqAlign << std::endl;
		TCharIter consIt = begin(cons);
		TCharIter consItEnd = end(cons);
		TCharIter consAlignIt = begin(consAlign);
		TCharIter seqAlignIt = begin(seqAlign);
		TSize row = value(begER, value(badIt)).i3;
		TSize col = b;
		while(consIt != consItEnd) {
			//std::cout << value(consIt) << ',' << value(consAlignIt) << ',' << value(seqAlignIt) << std::endl;
			if (value(consIt) == value(consAlignIt)) {
				// Both consensi have a gap
				value(mat, row * len + col) = value(seqAlignIt);
				goNext(consIt);
				goNext(consAlignIt);
				goNext(seqAlignIt);
				++col;
			} else {
				if (value(consIt) == gapChar) {
					while ((consIt != consItEnd) && (value(consIt) != value(consAlignIt))) {
						value(mat, row * len + col) = gapChar;
						goNext(consIt);
						++col;
					}
				} else {
					// This kind of read should be deleted
					while (value(consIt) != value(consAlignIt)) {
						goNext(consAlignIt);
					}
				}
			}
		}	
	}

	// Fix begin and end gaps
	itBE = begin(begER);
	itBEEnd = end(begER);
	for(;itBE!=itBEEnd; goNext(itBE)) {
		TSize b = value(itBE).i1;
		TSize e = value(itBE).i2;
		TSize row = value(itBE).i3;
		while (value(mat, row * len + b) == gapChar) {
			value(mat, row * len + b) = specialGap;
			++b;
		}
		while (value(mat, row * len + (e-1)) == gapChar) {
			value(mat, row * len + (e-1)) = specialGap;
			--e;
		}
		value(itBE).i1 = b;
		value(itBE).i2 = e;
	}


	return (length(badReads));
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TAlphabet, typename TSpec2, typename TGappedConsensus, typename TCoverage, typename TSize>
inline void
consensusCalling(String<TValue, TSpec> const& mat,
				 String<TAlphabet, TSpec2>& consensus,
				 TGappedConsensus& gappedConsensus,
				 TCoverage& coverage,
				 TSize alignDepth,
				 Majority_Vote)
{
	// Initialization
	TSize alphabetSize = ValueSize<TAlphabet>::VALUE;
	TValue gapChar = gapValue<TValue>();

	// Set-up the counters
	typedef String<TSize> TCounter;
	typedef String<TCounter> TCounters;
	TCounters counterValues;
	__countLetters(mat, counterValues, coverage, alignDepth, TAlphabet() );

	// Get the consensus
	typedef typename Iterator<TCounters>::Type TCounterIt;
	TCounterIt countIt = begin(counterValues);
	TCounterIt countItEnd = end(counterValues);
	clear(consensus);
	clear(gappedConsensus);
	for(;countIt != countItEnd; goNext(countIt)) {
		TSize max = 0;
		TValue c = TValue();
		typedef typename Iterator<TCounter>::Type TCIt;
		TCIt cIt = begin(value(countIt));
		TCIt cItEnd = end(value(countIt));
		TSize pos = 0;
		for(;cIt != cItEnd; goNext(cIt), ++pos) {
			if (value(cIt) > max) {
				max = value(cIt);
				if (pos == alphabetSize) c = gapChar;
				else c = TAlphabet(pos);
			}
		}
		if (c != gapChar) appendValue(consensus, c);
		appendValue(gappedConsensus, c);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TCargo, typename TSpec, typename TBegEndRowPos, typename TNewLibraryGraph>
inline unsigned int
realignLowQualityReads(Graph<Alignment<TStringSet, TCargo, TSpec> > const& gIn,
					   TBegEndRowPos const& readBegEndRowPos,
					   TNewLibraryGraph& gOut)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TBegEndRowPos>::Type TBegEndIter;

	// Initialization
	TStringSet& str = stringSet(gIn);
	clearVertices(gOut);
	
	// Find disrupted reads
	std::set<TId> unalignedRead;
	TBegEndIter beIt = begin(readBegEndRowPos);
	TBegEndIter beItEnd = end(readBegEndRowPos);
	TSize pos = 0;
	for(;beIt != beItEnd; ++beIt, ++pos) {
		TSize lenStr = length(value(str, pos));
		if (((value(beIt)).i2 - (value(beIt)).i1) > (lenStr + lenStr / 3 + 10)) unalignedRead.insert(positionToId(str, pos));
	}

	// Any disrupted reads
	if (unalignedRead.empty()) return 0;
	else exit(0);

	
	//// String of fragments to combine all pairwise alignments into a multiple alignment
	//typedef Fragment<> TFragment;
	//typedef String<TFragment> TFragmentString;
	//typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	//TFragmentString matches;
	//reserve(matches, numEdges(gIn));

	//// Insert all overlaps from the previous alignment
	//TEdgeIterator it_tmp(gIn);
	//for(;!atEnd(it_tmp);++it_tmp) appendValue(matches, TFragment(sequenceId(gIn, sourceVertex(it_tmp)),fragmentBegin(gIn, sourceVertex(it_tmp)), sequenceId(gIn, targetVertex(it_tmp)), fragmentBegin(gIn, targetVertex(it_tmp)), fragmentLength(gIn, sourceVertex(it_tmp))));
	//
	//// Recompute the overlap of interesting pairs
	//Score<int> score_type = Score<int>(2,-1,-4,-6);
	//TPairIter pairIt = begin(pList);
	//TPairIter pairItEnd = end(pList);
	//for(;pairIt != pairItEnd; ++pairIt) {
	//	TId id1 = (value(pairIt)).i1;
	//	TId id2 = (value(pairIt)).i2;

	//	if ((unalignedRead.find(id1) != unalignedRead.end()) || (unalignedRead.find(id2) != unalignedRead.end())) {
	//		// Make a pairwise string-set
	//		TStringSet pairSet;
	//		assignValueById(pairSet, str, id1);
	//		assignValueById(pairSet, str, id2);

	//		// Overlap alignment with a small mismatch score
	//		TSize from = length(matches);
	//		globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );

	//		// Determine a sequence weight
	//		TSize matchLen = 0;
	//		TSize overlapLen = 0;
	//		TSize alignLen = 0;
	//		getAlignmentStatistics(matches, pairSet, from, matchLen, overlapLen, alignLen);
	//		double quality = (double) matchLen / (double) overlapLen;

	//		// Take all overlaps of good quality
	//		if ((quality < 0.75) || (matchLen < 8)) {
	//			resize(matches, from);
	//		}
	//	}
	//}

	//// Refine all matches
	//matchRefinement(matches,stringSet(gOut),score_type, gOut);

	//return unalignedRead.size();
}


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
