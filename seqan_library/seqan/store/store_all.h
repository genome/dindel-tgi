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
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_STORE_ALL_H
#define SEQAN_HEADER_STORE_ALL_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Contig Store Configuration
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec = void>
struct FragmentStoreConfig 
{
	typedef String<Dna5Q>	TReadSeq;
	typedef String<Dna5Q>	TContigSeq;
	
	typedef double			TMean;
	typedef double			TStd;
	typedef unsigned char	TMappingQuality;
		
	typedef void					TReadStoreElementSpec;
	typedef Owner<ConcatDirect<> >	TReadSeqStoreSpec;
	typedef void					TMatePairStoreElementSpec;
	typedef void					TLibraryStoreElementSpec;
	typedef void					TContigStoreElementSpec;
	typedef void					TAlignedReadStoreElementSpec;
	typedef void					TAnnotationStoreElementSpec;
	typedef Owner<ConcatDirect<> >	TAlignedReadTagStoreSpec;
};

//////////////////////////////////////////////////////////////////////////////
// Fragment Store
//////////////////////////////////////////////////////////////////////////////

template <typename TSpec = void, typename TConfig = FragmentStoreConfig<TSpec> >
struct FragmentStore
{
private:
	typedef typename TConfig::TReadStoreElementSpec			TReadStoreElementSpec;
	typedef typename TConfig::TReadSeqStoreSpec				TReadSeqStoreSpec;
	typedef typename TConfig::TMatePairStoreElementSpec		TMatePairStoreElementSpec;
	typedef typename TConfig::TLibraryStoreElementSpec		TLibraryStoreElementSpec;
	typedef typename TConfig::TContigStoreElementSpec		TContigStoreElementSpec;
	typedef typename TConfig::TAlignedReadStoreElementSpec	TAlignedReadStoreElementSpec;
	typedef typename TConfig::TAlignedReadTagStoreSpec		TAlignedReadTagStoreSpec;
	typedef typename TConfig::TAnnotationStoreElementSpec	TAnnotationStoreElementSpec;

public:
	typedef typename TConfig::TMean					TMean;
	typedef typename TConfig::TStd					TStd;
	typedef typename TConfig::TMappingQuality		TMappingQuality;
	
	typedef typename TConfig::TReadSeq				TReadSeq;
	typedef typename TConfig::TContigSeq			TContigSeq;
	typedef typename Position<TReadSeq>::Type		TReadPos;
	typedef typename Position<TContigSeq>::Type		TContigPos;
	
	typedef GapAnchor<TReadPos>						TReadGapAnchor;
	typedef GapAnchor<TContigPos>					TContigGapAnchor;

	typedef String< ReadStoreElement< TReadSeq, TReadPos, TReadStoreElementSpec > >							TReadStore;
	typedef String< MatePairStoreElement< TMatePairStoreElementSpec > >										TMatePairStore;
	typedef String< LibraryStoreElement< TMean, TStd, TLibraryStoreElementSpec > >							TLibraryStore;
	typedef String< ContigStoreElement< TContigSeq, TContigGapAnchor, TContigStoreElementSpec > >			TContigStore;
	typedef String< AlignedReadStoreElement< TContigPos, TReadGapAnchor, TAlignedReadStoreElementSpec > >	TAlignedReadStore;
	typedef String< AnnotationStoreElement< TContigPos, TAnnotationStoreElementSpec > >						TAnnotationStore;
	
	// main containers
	TReadStore			readStore;			// readId     -> matePairId
	TMatePairStore		matePairStore;		// matePairId -> readId0, readId1, libraryId
	TLibraryStore		libraryStore;		// libraryId  -> libSizeMean, libSizeStd
	TContigStore		contigStore;		// contigId   -> contigSeq, contigGaps
	TAlignedReadStore	alignedReadStore;	//            -> id, readId, contigId, pairMatchId (not matePairId!), beginPos, endPos, gaps
	TAnnotationStore	annotationStore;	// annoId     -> parentId, contigId, beginPos, endPos

											// REMARKS: 
											// 1)
											//    beginPos <= endPos     forward strand
											//    beginPos >  endPos     backward strand (reverse complement)
											// 2) 
											//    The alignedReadStore can arbitrarily be resorted. The unique identifier id should
											//    be used to address additional information for each alignedRead in additional tables.

	// we store the read sequences in a seperate stringset to reduce the memory overhead 
	StringSet<TReadSeq, TReadSeqStoreSpec>	readSeqStore;

	// extra SAM fields
	String<TMappingQuality>	alignedReadQualityStore;
	StringSet<CharString, TAlignedReadTagStoreSpec> alignedReadTagStore;

	// retrieve the names of reads, mate-pairs, libraries, contigs, annotations by their ids
	StringSet<CharString>	readNameStore;
	StringSet<CharString>	matePairNameStore;
	StringSet<CharString>	libraryNameStore;
	StringSet<CharString>	contigNameStore;
	StringSet<CharString>	annotationNameStore;
};

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Read Store Accessors
//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TConfig>
inline void
clearReads(FragmentStore<TSpec, TConfig> &me)
{
	clear(me.readStore);
	clear(me.readSeqStore);
}

template <typename TSpec, typename TConfig, typename TRead, typename TId>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read, 
	TId matePairId)
{
	SEQAN_ASSERT(length(me.readStore) == length(me.readSeqStore))

	typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
	typename Value<TReadStore>::Type r;
	r.matePairId = matePairId;

	appendValue(me.readStore, r, Generous());
	appendValue(me.readSeqStore, read, Generous());
	return length(me.readStore);
}

template <typename TSpec, typename TConfig, typename TRead>
inline typename Size<typename FragmentStore<TSpec, TConfig>::TReadStore>::Type
appendRead(
	FragmentStore<TSpec, TConfig> &me, 
	TRead const &read)
{
	typedef typename FragmentStore<TSpec, TConfig>::TReadStore TReadStore;
    typedef typename Value<TReadStore>::Type TReadStoreElement;
    
	return appendRead(me, read, TReadStoreElement::INVALID_ID);
}

template <typename TSpec, typename TConfig, typename TId>
inline typename Value<typename FragmentStore<TSpec, TConfig>::TReadSeqStore>::Type
getRead(
	FragmentStore<TSpec, TConfig> &me, 
	TId id)
{
	return value(me.readSeqStore, id);
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
