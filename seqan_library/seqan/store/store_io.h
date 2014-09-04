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

#ifndef SEQAN_HEADER_STORE_IO_H
#define SEQAN_HEADER_STORE_IO_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Amos message file:
	Amos message file.
*/
struct TagAmos_;
typedef Tag<TagAmos_> const Amos;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Alignment Graph Format.value.FastaReadFormat:
	Fasta read format to write a multi-read alignment.
*/

struct FastaReadFormat_;
typedef Tag<FastaReadFormat_> const FastaReadFormat;

//////////////////////////////////////////////////////////////////////////////
// Auxillary functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////



template <typename TSpec, typename TConfig, typename TPos, typename TGapAnchor, typename TSpecAlign, typename TBeginClr, typename TEndClr>
inline void
getClrRange(FragmentStore<TSpec, TConfig> const& fragStore,
			AlignedReadStoreElement<TPos, TGapAnchor, TSpecAlign> const& alignEl,
			TBeginClr& begClr,		// Out-parameter: left / begin position of the clear range
			TEndClr& endClr)		// Out-parameter: right / end position of the clear range
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Iterator<String<TGapAnchor>, Standard>::Type TGapIter;
	
	TSize lenRead = length(fragStore.readSeqStore[alignEl.readId]);
	TGapIter itGap = begin(alignEl.gaps, Standard() );
	TGapIter itGapEnd = end(alignEl.gaps, Standard() );
	
	// Any gaps or clipped characters?
	if (itGap == itGapEnd) {
		begClr = 0;
		endClr = lenRead;
	} else {
		// Begin clear range
		begClr = (itGap->gapPos == 0) ? itGap->seqPos : 0;
		// End clear range
		--itGapEnd;
		if (itGapEnd->seqPos != lenRead) endClr = lenRead;
		else {
			int diff = (itGap != itGapEnd) ? (*(itGapEnd - 1)).gapPos - (*(itGapEnd-1)).seqPos : 0;
			int newDiff = itGapEnd->gapPos - itGapEnd->seqPos;
			endClr = (newDiff < diff) ? lenRead - (diff - newDiff) : lenRead;	
		}
	}

	// For reverse reads adapt clear ranges
	if (alignEl.beginPos > alignEl.endPos) {
		TBeginClr tmp = begClr;
		begClr = lenRead - endClr;
		endClr = lenRead - tmp;
	}
}



//////////////////////////////////////////////////////////////////////////////
// Read / Write of AMOS message files (*.afg)
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TSpec, typename TConfig>
inline void 
read(TFile & file,
	 FragmentStore<TSpec, TConfig>& fragStore,
	 Amos) 
{
	SEQAN_CHECKPOINT
	// Basic types
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Id<TFragmentStore>::Type TId;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef typename TFragmentStore::TReadSeq TReadSeq;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
	typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

	// All maps to mirror file ids to our ids
	typedef std::map<TId, TId> TIdMap;
	TIdMap libIdMap;
	TIdMap frgIdMap;
	TIdMap readIdMap;

	// Parse the file and convert the internal ids
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	while (!_streamEOF(file)) {
		if (_streamEOF(file)) break;

		// New block?
		if (c == '{') {
			c = _streamGet(file);
			String<char> blockIdentifier;
			_parse_readIdentifier(file, blockIdentifier, c);
			_parse_skipLine(file, c);

			// Library block
			if (blockIdentifier == "LIB") {
				TLibraryStoreElement libEl;
				TId id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				while (c != '}') {
					clear(fieldIdentifier);
					_parse_readIdentifier(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c, Generous());
							c = _streamGet(file);
						}
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "mea") {
						c = _streamGet(file);
						libEl.mean = _parse_readDouble(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "std") {
						c = _streamGet(file);
						libEl.std = _parse_readDouble(file, c);
						_parse_skipLine(file, c);
					} else {
						_parse_skipLine(file, c);
					}
				}
				libIdMap.insert(std::make_pair(id, length(fragStore.libraryStore)));
				appendValue(fragStore.libraryStore, libEl, Generous() );
				appendValue(fragStore.libraryNameStore, eid, Generous() );
			} else if (blockIdentifier == "FRG") {  // Fragment block
				TMatePairElement matePairEl;
				TId id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				bool foundRds = false;
				while (c != '}') {
					clear(fieldIdentifier);
					_parse_readIdentifier(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c, Generous() );
							c = _streamGet(file);
						}
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "lib") {
						c = _streamGet(file);
						matePairEl.libId = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "rds") {
						foundRds = true;
						c = _streamGet(file);
						matePairEl.readId[0] = _parse_readNumber(file, c);
						c = _streamGet(file);
						matePairEl.readId[1] = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else {
						_parse_skipLine(file, c);
					}
				}
				// Only insert valid mate pairs
				if (foundRds) {
					frgIdMap.insert(std::make_pair(id, length(fragStore.matePairStore)));
					appendValue(fragStore.matePairStore, matePairEl, Generous() );
					appendValue(fragStore.matePairNameStore, eid, Generous() );
				}
			} else if (blockIdentifier == "RED") {   // Read block
				TId id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				String<char> qual;
				TId matePairId = 0;
				TReadSeq seq;
				while (c != '}') {
					clear(fieldIdentifier);
					_parse_readIdentifier(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c, Generous() );
							c = _streamGet(file);
						}
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "frg") {
						c = _streamGet(file);
						matePairId = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "seq") {
						c = _streamGet(file);
						_parse_skipWhitespace(file, c);
						while (c != '.') {
							_parse_readSequenceData(file,c,seq);
							_parse_skipWhitespace(file, c);
						}
					} else if (fieldIdentifier == "qlt") {
						clear(qual);
						c = _streamGet(file);
						_parse_skipWhitespace(file, c);
						while (c != '.') {
							if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) appendValue(qual, c, Generous() );
							c = _streamGet(file);
						}
					} else {
						_parse_skipLine(file, c);
					}
				}
				// Set quality
				typedef typename Iterator<TReadSeq, Standard>::Type TReadIter;
				typedef typename Iterator<String<char> >::Type TQualIter;
				TReadIter begIt = begin(seq, Standard() );
				TQualIter qualIt = begin(qual);
				TQualIter qualItEnd = end(qual);
				for(;qualIt != qualItEnd; goNext(qualIt), goNext(begIt)) assignQualityValue(value(begIt), value(qualIt));

				// Insert the read
				readIdMap.insert(std::make_pair(id, length(fragStore.readStore)));
				appendRead(fragStore, seq, matePairId);
				appendValue(fragStore.readNameStore, eid, Generous() );
			} else if (blockIdentifier == "CTG") {   // Contig block
				TContigElement contigEl;
				TSize fromAligned = length(fragStore.alignedReadStore);
				TId id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				String<char> contigSeq;
				String<char> contigQual;
				while (c != '}') {
					// Are we entering a TLE block
					if (c == '{') {
						TAlignedElement alignEl;
						String<char> fdIdentifier;
						typedef typename TFragmentStore::TContigPos TContigPos;
						TContigPos offsetPos = 0;
						TContigPos clr1 = 0;
						TContigPos clr2 = 0;
						String<TContigPos> gaps;
						while (c != '}') {
							clear(fdIdentifier);
							_parse_readIdentifier(file, fdIdentifier, c);
							if (fdIdentifier == "src") {
								c = _streamGet(file);
								alignEl.readId = _parse_readNumber(file, c);
								_parse_skipLine(file, c);
							} else if (fdIdentifier == "off") {
								c = _streamGet(file);
								if (c != '-') offsetPos = _parse_readNumber(file, c);
								else offsetPos = 0;
								_parse_skipLine(file, c);
							} else if (fdIdentifier == "clr") {
								c = _streamGet(file);
								clr1 = _parse_readNumber(file, c);
								c = _streamGet(file);
								clr2 = _parse_readNumber(file, c);
								_parse_skipLine(file, c);
							} else if (fdIdentifier == "gap") {
								c = _streamGet(file);
								_parse_skipWhitespace(file, c);
								while (c != '.') {
									if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) {
										TSize nextGap = _parse_readNumber(file, c);
										appendValue(gaps, nextGap, Generous() );
									}
									c = _streamGet(file);
								}
							} else {
								_parse_skipLine(file, c);
							}
						}
						_parse_skipLine(file, c);

						// Get the length of the read
						TId readId = (readIdMap.find(alignEl.readId))->second;
						TSize lenRead = length(value(fragStore.readSeqStore, readId));

						// Create the gap anchors
						typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
						int offset = 0;
						if ((clr1 < clr2) && (clr1>0)) offset = clr1;
						else if ((clr1 > clr2) && (clr1 < lenRead)) offset = lenRead - clr1;
						int diff = -1 * (int) (offset);
						// Clipped begin
						if (offset != 0) appendValue(alignEl.gaps, TContigGapAnchor(offset, 0), Generous() );
						// Internal gaps
						typedef typename Iterator<String<TContigPos>, Standard>::Type TPosIter;
						TPosIter posIt = begin(gaps, Standard() ); 
						TPosIter posItEnd = end(gaps, Standard() );
						TContigPos lastGap = 0;
						TSize gapLen = 0;
						TSize totalGapLen = 0;
						for(;posIt!=posItEnd; goNext(posIt)) {
							if (gapLen == 0) {
								++gapLen; ++totalGapLen;
								++diff;
								lastGap = value(posIt);
							} 
							else if (lastGap == value(posIt)) {
								++gapLen; ++totalGapLen;
								++diff;
							}
							else {
								appendValue(alignEl.gaps, TContigGapAnchor(offset + lastGap, offset + lastGap + diff), Generous() );
								gapLen = 1; ++totalGapLen;
								lastGap = value(posIt);
								++diff;
							}
						}
						if (gapLen > 0) appendValue(alignEl.gaps, TContigGapAnchor(offset + lastGap, offset + lastGap + diff), Generous() );
						// Clipped end
						if ((clr1 < clr2) && (clr2 < lenRead)) {
							diff -= (lenRead - clr2);				
							appendValue(alignEl.gaps, TContigGapAnchor(lenRead, lenRead + diff), Generous() );
						} else if ((clr1 > clr2) && (clr2 > 0)) {
							diff -= clr2;
							appendValue(alignEl.gaps, TContigGapAnchor(lenRead, lenRead + diff), Generous() );
						}
						
						// Set begin and end position
						if (clr1 < clr2) {
							alignEl.beginPos = offsetPos;
							alignEl.endPos = offsetPos + totalGapLen + (clr2 - clr1);
						} else {
							alignEl.beginPos = offsetPos + totalGapLen + (clr1 - clr2);
							alignEl.endPos = offsetPos;
						}
		
						// Append new align fragment, note: contigId must still be set
						alignEl.id = length(fragStore.alignedReadStore);
						appendValue(fragStore.alignedReadStore, alignEl, Generous() );
					} else {
						clear(fieldIdentifier);
						_parse_readIdentifier(file, fieldIdentifier, c);
						if (fieldIdentifier == "iid") {
							c = _streamGet(file);
							id = _parse_readNumber(file, c);
							_parse_skipLine(file, c);
						} else if (fieldIdentifier == "eid") {
							c = _streamGet(file);
							while ((c != '\n') && (c != '\r')) {
								appendValue(eid, c, Generous() );
								c = _streamGet(file);
							}
							_parse_skipLine(file, c);
						} else if (fieldIdentifier == "seq") {
							c = _streamGet(file);
							_parse_skipWhitespace(file, c);
							while (c != '.') {
								do {
									_parse_readSequenceData(file,c,contigSeq);
								} while (c == '-');
								_parse_skipWhitespace(file, c);
							}
						} else if (fieldIdentifier == "qlt") {
							c = _streamGet(file);
							_parse_skipWhitespace(file, c);
							while (c != '.') {
								if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) {
									appendValue(contigQual, c, Generous() );
								}
								c = _streamGet(file);
							}
						} else {
							_parse_skipLine(file, c);
						}
					}
				}

				// Create the gap anchors
				char gapChar = gapValue<char>();
				typedef typename Iterator<String<char> >::Type TStringIter;
				TStringIter seqIt = begin(contigSeq);
				TStringIter seqItEnd = end(contigSeq);
				TStringIter qualIt = begin(contigQual);
				typedef typename TFragmentStore::TReadPos TPos;
				typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
				TPos ungappedPos = 0;
				TPos gappedPos = 0;
				bool gapOpen = false;
				for(;seqIt != seqItEnd; goNext(seqIt), goNext(qualIt), ++gappedPos) {
					if (value(seqIt) == gapChar) gapOpen = true;				
					else {
						if (gapOpen) {
							appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous() );
							gapOpen = false;
						}
						Dna5Q letter = value(seqIt);
						assignQualityValue(letter, value(qualIt));
						appendValue(contigEl.seq, letter, Generous() );
						++ungappedPos;
					}
				}
				if (gapOpen) appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos), Generous() );

				// Set the contigId in all aligned reads
				TSize toAligned = length(fragStore.alignedReadStore);
				TId newContigId = length(fragStore.contigStore);
				for(; fromAligned < toAligned; ++fromAligned) {
					(value(fragStore.alignedReadStore, fromAligned)).contigId = newContigId;
				}

				// Insert the contig
				appendValue(fragStore.contigStore, contigEl, Generous() );
				appendValue(fragStore.contigNameStore, eid, Generous() );
			} else {
				_parse_skipLine(file, c);
			}	
		} else {
			_parse_skipLine(file, c);
		}
	}

	// Renumber all ids
	typedef typename TIdMap::const_iterator TIdMapIter;
	typedef typename Iterator<typename TFragmentStore::TMatePairStore>::Type TMateIter;
	TMateIter mateIt = begin(fragStore.matePairStore);
	TMateIter mateItEnd = end(fragStore.matePairStore);
	for(;mateIt != mateItEnd; goNext(mateIt)) {
		if (mateIt->libId != TMatePairElement::INVALID_ID) {
			TIdMapIter libIdPos = libIdMap.find(mateIt->libId);
			if (libIdPos != libIdMap.end()) mateIt->libId = libIdPos->second;
			else mateIt->libId = TMatePairElement::INVALID_ID;
		}
		if (mateIt->readId[0] != TMatePairElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(mateIt->readId[0]);
			if (readIdPos != readIdMap.end()) mateIt->readId[0] = readIdPos->second;
			else mateIt->readId[0] = TMatePairElement::INVALID_ID;
		}
		if (mateIt->readId[1]!= TMatePairElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(mateIt->readId[1]);
			if (readIdPos != readIdMap.end()) mateIt->readId[1] = readIdPos->second;
			else mateIt->readId[0] = TMatePairElement::INVALID_ID;
		}
	}
	typedef typename Iterator<typename TFragmentStore::TReadStore>::Type TReadIter;
	TReadIter readIt = begin(fragStore.readStore);
	TReadIter readItEnd = end(fragStore.readStore);
	for(;readIt != readItEnd; goNext(readIt)) {
		if (readIt->matePairId != TReadStoreElement::INVALID_ID) {
			TIdMapIter mateIdPos = frgIdMap.find(readIt->matePairId);
			if (mateIdPos != frgIdMap.end()) readIt->matePairId = mateIdPos->second;
			else readIt->matePairId = TReadStoreElement::INVALID_ID;
		}
	}
	TId myPairMatchId = 0;  // Dummy variable to count the matches
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		if (alignIt->readId != TAlignedElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(alignIt->readId);
			if (readIdPos != readIdMap.end()) alignIt->readId = readIdPos->second;
			else alignIt->readId = TAlignedElement::INVALID_ID;
		}
		alignIt->pairMatchId = myPairMatchId++;
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec, typename TConfig>
inline void 
write(TFile & target,
	  FragmentStore<TSpec, TConfig>& fragStore,
	  Amos) 
{
	SEQAN_CHECKPOINT
	// Basic types
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Id<TFragmentStore>::Type TId;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
	typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

	// Write Header
	_streamWrite(target,"{UNV\niid:1\neid:seqan\ncom:\nafg file created with SeqAn\n.\n}\n");
	
	// Write Libraries
	typedef typename Iterator<typename TFragmentStore::TLibraryStore, Standard>::Type TLibIter;
	TLibIter libIt = begin(fragStore.libraryStore, Standard() );
	TLibIter libItEnd = end(fragStore.libraryStore, Standard() );
	bool noNamesPresent = (length(fragStore.libraryNameStore) == 0);
	for(TSize idCount = 0;libIt != libItEnd; goNext(libIt), ++idCount) {
		_streamWrite(target,"{LIB\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.libraryNameStore, idCount));
			_streamPut(target, '\n');
		}
		_streamWrite(target,"{DST\n");
		_streamWrite(target,"mea:");
		_streamPutFloat(target, libIt->mean);
		_streamPut(target, '\n');
		_streamWrite(target,"std:");
		_streamPutFloat(target, libIt->std);
		_streamPut(target, '\n');
		_streamWrite(target,"}\n");	
		_streamWrite(target,"}\n");
	}

	// Write Fragments / mate pairs
	typedef typename Iterator<typename TFragmentStore::TMatePairStore, Standard>::Type TMateIter;
	TMateIter mateIt = begin(fragStore.matePairStore, Standard() );
	TMateIter mateItEnd = end(fragStore.matePairStore, Standard() );
	noNamesPresent = (length(fragStore.matePairNameStore) == 0);
	for(TSize idCount = 0;mateIt != mateItEnd; goNext(mateIt), ++idCount) {
		_streamWrite(target,"{FRG\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.matePairNameStore, idCount));
			_streamPut(target, '\n');
		}
		_streamWrite(target,"lib:");
		_streamPutInt(target, mateIt->libId + 1);
		_streamPut(target, '\n');
		if ((mateIt->readId[0] != TMatePairElement::INVALID_ID) && (mateIt->readId[1] != TMatePairElement::INVALID_ID)) {
			_streamWrite(target,"rds:");
			_streamPutInt(target, mateIt->readId[0] + 1);
			_streamPut(target, ',');
			_streamPutInt(target, mateIt->readId[1] + 1);
			_streamPut(target, '\n');
		}
		_streamWrite(target,"}\n");
	}

	// Get clear ranges
	typedef Pair<typename TFragmentStore::TReadPos, typename TFragmentStore::TReadPos> TClrRange;
	String<TClrRange> clrRange;
	fill(clrRange, length(fragStore.readStore), TClrRange(0,0));
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore, Standard>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore, Standard() );
	TAlignIter alignItEnd = end(fragStore.alignedReadStore, Standard() );
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		typename TFragmentStore::TReadPos begClr = 0;
		typename TFragmentStore::TReadPos endClr = 0;
		getClrRange(fragStore, value(alignIt), begClr, endClr);
		value(clrRange, alignIt->readId) = TClrRange(begClr, endClr);
	}

	// Write reads
	typedef typename Iterator<typename TFragmentStore::TReadStore, Standard>::Type TReadIter;
	TReadIter readIt = begin(fragStore.readStore, Standard() );
	TReadIter readItEnd = end(fragStore.readStore, Standard() );
	noNamesPresent = (length(fragStore.readNameStore) == 0);
	for(TSize idCount = 0;readIt != readItEnd; ++readIt, ++idCount) {
		_streamWrite(target,"{RED\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.readNameStore, idCount));
			_streamPut(target, '\n');
		}
		_streamWrite(target,"seq:\n");
		typedef typename Iterator<typename TFragmentStore::TReadSeq>::Type TSeqIter;
		typedef typename Value<typename TFragmentStore::TReadSeq>::Type TAlphabet;
		TSeqIter seqIt = begin(value(fragStore.readSeqStore, idCount));
		TSeqIter seqItEnd = end(value(fragStore.readSeqStore, idCount));
		for(TSize k = 0;seqIt!=seqItEnd;goNext(seqIt), ++k) {
			if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
			_streamPut(target, getValue(seqIt));
		}
		_streamWrite(target, "\n.\n");
		_streamWrite(target,"qlt:\n");
		seqIt = begin(value(fragStore.readSeqStore, idCount));
		for(TSize k = 0;seqIt!=seqItEnd;goNext(seqIt), ++k) {
			if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
			Ascii c = ' ';
			convertQuality(c, getQualityValue(value(seqIt)));
			_streamPut(target, c);
		}
		_streamWrite(target, "\n.\n");
		if (readIt->matePairId != TReadStoreElement::INVALID_ID) {
			_streamWrite(target,"frg:");
			_streamPutInt(target, readIt->matePairId + 1);
			_streamPut(target, '\n');
		}
		if ((value(clrRange, idCount)).i1 != (value(clrRange, idCount)).i2) {
			_streamWrite(target,"clr:");
			_streamPutInt(target, (value(clrRange, idCount)).i1);
			_streamPut(target, ',');
			_streamPutInt(target, (value(clrRange, idCount)).i2);
			_streamPut(target, '\n');
		}
		_streamWrite(target,"}\n");
	}

	// Sort aligned reads according to contigId
	sortAlignedReads(fragStore.alignedReadStore, SortContigId());

	// Write Contigs
	typedef typename Iterator<typename TFragmentStore::TContigStore, Standard>::Type TContigIter;
	TContigIter contigIt = begin(fragStore.contigStore, Standard() );
	TContigIter contigItEnd = end(fragStore.contigStore, Standard() );
	alignIt = begin(fragStore.alignedReadStore);
	alignItEnd = end(fragStore.alignedReadStore);
	noNamesPresent = (length(fragStore.contigNameStore) == 0);
	for(TSize idCount = 0;contigIt != contigItEnd; goNext(contigIt), ++idCount) {
		_streamWrite(target,"{CTG\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.contigNameStore, idCount));
			_streamPut(target, '\n');
		}
		String<char> qlt;
		_streamWrite(target,"seq:\n");
		typedef typename Iterator<typename TFragmentStore::TContigSeq>::Type TContigIter;
		TContigIter seqContigIt = begin(contigIt->seq);
		TContigIter seqContigItEnd = end(contigIt->seq);
		typedef typename Iterator<String<typename TFragmentStore::TContigGapAnchor> >::Type TGapsIter;
		TGapsIter itGaps = begin(contigIt->gaps);
		TGapsIter itGapsEnd = end(contigIt->gaps);
		int diff = 0;
		char gapChar = gapValue<char>();
		typename TFragmentStore::TContigPos mySeqPos = 0;
		TSize k = 0;
		for(;itGaps != itGapsEnd; goNext(itGaps)) {
			while (mySeqPos < itGaps->seqPos) {
				if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
				++k;
				_streamPut(target, value(seqContigIt));
				Ascii c = ' ';
				convertQuality(c, getQualityValue(value(seqContigIt)));
				appendValue(qlt, c, Generous() );
				goNext(seqContigIt);++mySeqPos;
			}
			for(int i = 0; i < ((int) itGaps->gapPos - (int) itGaps->seqPos) - diff; ++i) {
				if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
				++k;
				_streamPut(target, gapChar);
				appendValue(qlt, '0', Generous() );
			}
			diff = (itGaps->gapPos - itGaps->seqPos);
		}
		for(;seqContigIt != seqContigItEnd; goNext(seqContigIt)) {
			if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
			++k;
			_streamPut(target, value(seqContigIt));
			Ascii c = ' ';
			convertQuality(c, getQualityValue(value(seqContigIt)));
			appendValue(qlt, c, Generous() );
		}
		_streamWrite(target, "\n.\n");
		_streamWrite(target,"qlt:\n");
		for(TSize k = 0;k<length(qlt); k+=60) {
			TSize endK = k + 60;
			if (endK > length(qlt)) endK = length(qlt);
			_streamWrite(target, infix(qlt, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		
		while ((alignIt != alignItEnd) && (idCount < alignIt->contigId)) goNext(alignIt);
		for(;(alignIt != alignItEnd) && (idCount == alignIt->contigId); goNext(alignIt)) {
			_streamWrite(target,"{TLE\n");
			_streamWrite(target,"src:");
			_streamPutInt(target, alignIt->readId + 1);
			_streamPut(target, '\n');
			typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor> >::Type TReadGapsIter;
			TReadGapsIter itGaps = begin(alignIt->gaps);
			TReadGapsIter itGapsEnd = end(alignIt->gaps);

			// Create the gaps string and the clear ranges
			typename TFragmentStore::TReadPos lenRead = length(value(fragStore.readSeqStore, alignIt->readId));
			TSize clr1 = 0;
			TSize clr2 = lenRead;
			// Create first clear range
			if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) clr1 = itGaps->seqPos;
			int diff = clr1;
			String<unsigned int> gaps;
			for(;itGaps != itGapsEnd; goNext(itGaps)) {
				for(int i = 0; i< diff - ((int) itGaps->seqPos - (int) itGaps->gapPos); ++i) {
					appendValue(gaps, itGaps->seqPos - clr1, Generous() );
				}
				// Clipped sequence
				if (diff - ((int) itGaps->seqPos - (int) itGaps->gapPos) < 0) {
					clr2 = lenRead + diff - ((int) itGaps->seqPos - (int) itGaps->gapPos);
				}
				diff = ((int) itGaps->seqPos - (int) itGaps->gapPos);
			}
			if (alignIt->beginPos > alignIt->endPos) {
				clr1 = lenRead - clr1;
				clr2 = lenRead - clr2;
			}
			_streamWrite(target,"off:");
			if (alignIt->beginPos < alignIt->endPos) _streamPutInt(target, alignIt->beginPos);
			else _streamPutInt(target, alignIt->endPos);
			_streamPut(target, '\n');
			_streamWrite(target,"clr:");
			_streamPutInt(target, clr1);
			_streamPut(target, ',');
			_streamPutInt(target, clr2);
			_streamPut(target, '\n');
			if (length(gaps)) {
				_streamWrite(target,"gap:\n");
				for(TSize z = 0;z<length(gaps); ++z) {
					_streamPutInt(target, value(gaps, z));
					_streamPut(target, '\n');
				}
				_streamWrite(target, ".\n");
			}
			_streamWrite(target,"}\n");
		}
		_streamWrite(target,"}\n");
	}
}


//////////////////////////////////////////////////////////////////////////////
// Read simulator format: Simple fasta read file with positions
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TSpec, typename TConfig, typename TFilePath>
inline bool 
_convertSimpleReadFile(TFile& file,
					   FragmentStore<TSpec, TConfig>& fragStore,
					   TFilePath& filePath, 
					   bool moveToFront)
{
	SEQAN_CHECKPOINT
	// Basic types
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Id<TFragmentStore>::Type TId;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef typename TFragmentStore::TContigPos TPos;
	typedef typename TFragmentStore::TReadSeq TReadSeq;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
	typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
	

	// All maps to mirror file ids to our internal ids
	typedef std::map<TId, TId> TIdMap;
	TIdMap libIdMap;
	TIdMap frgIdMap;
	TIdMap readIdMap;


	// Parse the file and convert the internal ids
	TPos maxPos = 0;
	TPos minPos = SupremumValue<TPos>::VALUE;
	TId count = 0;
	TValue c;
	if ((!file) || (_streamEOF(file))) return false;
	else c = _streamGet(file);
	while (!_streamEOF(file)) {
		if (_streamEOF(file)) break;

		// New read?
		if (c == '>') {
			TAlignedElement alignEl;
			TId id = count;
			TId fragId = count;
			TId repeatId = 0;
			
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);

			// Get the layout positions
			alignEl.beginPos = _parse_readNumber(file, c);
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
			alignEl.endPos = _parse_readNumber(file, c);
			
			// Any attributes?
			String<char> eid;
			String<char> qlt;
			TReadSeq seq;
			if (c == '[') {
				String<char> fdIdentifier;
				while (c != ']') {
					c = _streamGet(file);
					_parse_skipWhitespace(file, c);
					clear(fdIdentifier);
					_parse_readIdentifier(file, fdIdentifier, c);
					if (fdIdentifier == "id") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
					} else if (fdIdentifier == "fragId") {
						c = _streamGet(file);
						fragId = _parse_readNumber(file, c);
					} else if (fdIdentifier == "repeatId") {
						c = _streamGet(file);
						repeatId = _parse_readNumber(file, c);
					} else if (fdIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != ',') && (c != ']')) {
							appendValue(eid, c, Generous());
							c = _streamGet(file);
						}
					} else if (fdIdentifier == "qlt") {
						c = _streamGet(file);
						while ((c != ',') && (c != ']')) {
							appendValue(qlt, c, Generous());
							c = _streamGet(file);
						}
					} else {
						// Jump to next attribute
						while ((c != ',') && (c != ']')) {
							c = _streamGet(file);
						}
					}
				}
			}
			_parse_skipLine(file, c);
			_parse_skipWhitespace(file, c);
			while ((!_streamEOF(file)) && (c != '>')) {
				_parse_readSequenceData(file,c, seq);
				_parse_skipWhitespace(file, c);
			}
			
			// Set quality
			typedef typename Iterator<TReadSeq, Standard>::Type TReadIter;
			typedef typename Iterator<String<char> >::Type TQualIter;
			TReadIter begIt = begin(seq, Standard() );
			TReadIter begItEnd = begin(seq, Standard() );
			if (length(qlt)) {
				TQualIter qualIt = begin(qlt);
				TQualIter qualItEnd = end(qlt);
				for(;qualIt != qualItEnd; goNext(qualIt), goNext(begIt)) assignQualityValue(value(begIt), value(qualIt));
			} else {
				for(;begIt != begItEnd; goNext(begIt)) assignQualityValue(value(begIt), 'D');
			}

			// Set eid if not given
			if (empty(eid)) {
				std::stringstream input;
				input << "R" << id;
				input << "-" << repeatId;
				eid = input.str().c_str();
			}

			// Insert the read
			readIdMap.insert(std::make_pair(id, length(fragStore.readStore)));
			appendRead(fragStore, seq, fragId);
			appendValue(fragStore.readNameStore, eid, Generous());

			// Insert an aligned read
			TSize readLen = length(seq);
			if (alignEl.beginPos < alignEl.endPos) {
				if (readLen != alignEl.endPos - alignEl.beginPos) {
					alignEl.endPos = alignEl.beginPos + readLen;
				}
				if (alignEl.beginPos < minPos) minPos = alignEl.beginPos;
				if (alignEl.endPos > maxPos) maxPos = alignEl.endPos;
			} else {
				if (readLen != alignEl.beginPos - alignEl.endPos) {
					alignEl.beginPos = alignEl.endPos + readLen;
				}
				if (alignEl.endPos < minPos) minPos = alignEl.endPos;
				if (alignEl.beginPos > maxPos) maxPos = alignEl.beginPos;
			}
			alignEl.readId = id;
			alignEl.pairMatchId =  fragId;
			alignEl.contigId = 0;
			alignEl.id = length(fragStore.alignedReadStore);
			appendValue(fragStore.alignedReadStore, alignEl, Generous());
			++count;
		} else {
			_parse_skipLine(file, c);
		}
	}

	// Read contig or reference sequence
	TContigElement contigEl;
	std::string fileName = filePath + 'S';
	FILE* strmRef = fopen(fileName.c_str(), "rb");
	String<char> contigEid = "C0";
	if ((strmRef) && (!_streamEOF(strmRef))) {
		c = _streamGet(strmRef);
		while (!_streamEOF(strmRef)) {
			if (_streamEOF(strmRef)) break;
			if (c == '>') {
				clear(contigEid);
				c = _streamGet(strmRef);
				while ((c != '\r') && (c != '\n')) {
					appendValue(contigEid, c, Generous());
					c = _streamGet(strmRef);
				}
				_parse_skipLine(strmRef, c);
				_parse_skipWhitespace(strmRef, c);
				while ((!_streamEOF(strmRef)) && (c != '>')) {
					_parse_readSequenceData(strmRef,c,contigEl.seq);
					_parse_skipWhitespace(strmRef, c);
				}
			} else {
				_parse_skipLine(strmRef, c);
			}
		}
		fclose(strmRef);
	}
	if (empty(contigEl.seq)) {
		typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
		if (moveToFront) appendValue(contigEl.gaps, TContigGapAnchor(0, maxPos - minPos), Generous());
		else appendValue(contigEl.gaps, TContigGapAnchor(0, maxPos), Generous());
	}
	appendValue(fragStore.contigStore, contigEl, Generous());
	appendValue(fragStore.contigNameStore, contigEid, Generous());


	// Read fragments
	fileName = filePath + 'F';
	FILE* strmFrag = fopen(fileName.c_str(), "rb");
	if ((strmFrag) && (!_streamEOF(strmFrag))) {
		c = _streamGet(strmFrag);
		while (!_streamEOF(strmFrag)) {
			if (_streamEOF(strmFrag)) break;
			if (c == '>') {
				TMatePairElement matePairEl;
				c = _streamGet(strmFrag);
				_parse_skipWhitespace(strmFrag, c);

				// Get the fragment id
				TId id = _parse_readNumber(strmFrag, c);
			
				// Any attributes?
				std::stringstream input;
				input << "F" << id;
				String<char> eid(input.str().c_str());
				if (c == '[') {
					String<char> fdIdentifier;
					while (c != ']') {
						c = _streamGet(strmFrag);
						_parse_skipWhitespace(strmFrag, c);
						clear(fdIdentifier);
						_parse_readIdentifier(strmFrag, fdIdentifier, c);
						if (fdIdentifier == "libId") {
							c = _streamGet(strmFrag);
							matePairEl.libId = _parse_readNumber(strmFrag, c);
						} else if (fdIdentifier == "eid") {
							clear(eid);
							c = _streamGet(strmFrag);
							while ((c != ',') && (c != ']')) {
								appendValue(eid, c, Generous());
								c = _streamGet(strmFrag);
							}
						} else {
							// Jump to next attribute
							while ((c != ',') && (c != ']')) {
								c = _streamGet(strmFrag);
							}
						}
					}
				}
				_parse_skipLine(strmFrag, c);
				_parse_skipWhitespace(strmFrag, c);

				// Read the two reads belonging to this mate pair
				matePairEl.readId[0] = _parse_readNumber(strmFrag, c);
				c = _streamGet(strmFrag);
				_parse_skipWhitespace(strmFrag, c);
				matePairEl.readId[1] = _parse_readNumber(strmFrag, c);
				_parse_skipLine(strmFrag, c);

				// Insert mate pair
				if (matePairEl.readId[0] != matePairEl.readId[1]) {
					frgIdMap.insert(std::make_pair(id, length(fragStore.matePairStore)));
					appendValue(fragStore.matePairStore, matePairEl, Generous());
					appendValue(fragStore.matePairNameStore, eid, Generous());
				}
			} else {
				_parse_skipLine(strmFrag, c);
			}
		}
		fclose(strmFrag);
	}
	

	// Read libraries
	fileName = filePath + 'L';
	FILE* strmLib = fopen(fileName.c_str(), "rb");
	if ((strmLib) && (!_streamEOF(strmLib))) {
		c = _streamGet(strmLib);
		while (!_streamEOF(strmLib)) {
			if (_streamEOF(strmLib)) break;
			if (c == '>') {

				TLibraryStoreElement libEl;
				c = _streamGet(strmLib);
				_parse_skipWhitespace(strmLib, c);

				// Get the fragment id
				TId id = _parse_readNumber(strmLib, c);
			
				// Any attributes?
				std::stringstream input;
				input << "L" << id;
				String<char> eid(input.str().c_str());
				if (c == '[') {
					String<char> fdIdentifier;
					while (c != ']') {
						c = _streamGet(strmLib);
						_parse_skipWhitespace(strmLib, c);
						clear(fdIdentifier);
						_parse_readIdentifier(strmLib, fdIdentifier, c);
						if (fdIdentifier == "eid") {
							clear(eid);
							c = _streamGet(strmLib);
							while ((c != ',') && (c != ']')) {
								appendValue(eid, c, Generous());
								c = _streamGet(strmLib);
							}
						} else {
							// Jump to next attribute
							while ((c != ',') && (c != ']')) {
								c = _streamGet(strmLib);
							}
						}
					}
				}
				_parse_skipLine(strmLib, c);
				_parse_skipWhitespace(strmLib, c);

				// Read the mean and standard deviation
				libEl.mean = _parse_readNumber(strmLib, c);
				c = _streamGet(strmLib);
				_parse_skipWhitespace(strmLib, c);
				libEl.std = _parse_readNumber(strmLib, c);
				_parse_skipLine(strmLib, c);

				// Insert mate pair
				libIdMap.insert(std::make_pair(id, length(fragStore.libraryStore)));
				appendValue(fragStore.libraryStore, libEl, Generous());
				appendValue(fragStore.libraryNameStore, eid, Generous());
			} else {
				_parse_skipLine(strmLib, c);
			}
		}
		fclose(strmLib);
	}
	
	
	// Renumber all ids
	typedef typename TIdMap::const_iterator TIdMapIter;
	typedef typename Iterator<typename TFragmentStore::TMatePairStore>::Type TMateIter;
	TMateIter mateIt = begin(fragStore.matePairStore);
	TMateIter mateItEnd = end(fragStore.matePairStore);
	for(;mateIt != mateItEnd; goNext(mateIt)) {
		if (mateIt->libId != TMatePairElement::INVALID_ID) {
			TIdMapIter libIdPos = libIdMap.find(mateIt->libId);
			if (libIdPos != libIdMap.end()) mateIt->libId = libIdPos->second;
			else mateIt->libId = TMatePairElement::INVALID_ID;
		}
		if (mateIt->readId[0] != TMatePairElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(mateIt->readId[0]);
			if (readIdPos != readIdMap.end()) mateIt->readId[0] = readIdPos->second;
			else mateIt->readId[0] = TMatePairElement::INVALID_ID;
		}
		if (mateIt->readId[1]!= TMatePairElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(mateIt->readId[1]);
			if (readIdPos != readIdMap.end()) mateIt->readId[1] = readIdPos->second;
			else mateIt->readId[0] = TMatePairElement::INVALID_ID;
		}
	}
	typedef typename Iterator<typename TFragmentStore::TReadStore>::Type TReadIter;
	TReadIter readIt = begin(fragStore.readStore);
	TReadIter readItEnd = end(fragStore.readStore);
	for(;readIt != readItEnd; goNext(readIt)) {
		if (readIt->matePairId != TReadStoreElement::INVALID_ID) {
			TIdMapIter mateIdPos = frgIdMap.find(readIt->matePairId);
			if (mateIdPos != frgIdMap.end()) readIt->matePairId = mateIdPos->second;
			else readIt->matePairId = TReadStoreElement::INVALID_ID;
		}
	}
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		if (alignIt->readId != TAlignedElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(alignIt->readId);
			if (readIdPos != readIdMap.end()) alignIt->readId = readIdPos->second;
			else alignIt->readId = TAlignedElement::INVALID_ID;
		}
		if (moveToFront) {
			alignIt->beginPos -= minPos;
			alignIt->endPos -= minPos;
		}
	}
	return true;
}


//////////////////////////////////////////////////////////////////////////////
// Old proprietary FastaReadFormat
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TMatrix, typename TSize2, typename TSize, typename TReadSlot> 
inline bool
convertAlignment(FragmentStore<TSpec, TConfig>& fragStore,
				 TMatrix& mat,
				 TSize2 contigId,
				 TSize& coverage,
				 TReadSlot& slot)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<TMatrix>::Type TValue;

	// Gap char
	TValue gapChar = gapValue<TValue>();

	// Sort according to contigId
	sortAlignedReads(fragStore.alignedReadStore, SortContigId());
	
	// Find range of the given contig
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore, Standard>::Type TAlignIter;
	TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());

	// Sort the reads according to the begin position
	sortAlignedReads(infix(fragStore.alignedReadStore, alignIt - begin(fragStore.alignedReadStore, Standard()), alignItEnd - begin(fragStore.alignedReadStore, Standard())), SortBeginPos());
	TAlignIter alignItBegin = alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());
	alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, contigId, SortContigId());

	// Get the maximum coverage and the slot for each read
	typedef String<TSize> TFirstFreePos;
	typedef typename Iterator<TFirstFreePos, Standard>::Type TPosIter;
	TFirstFreePos freePos;
	TSize pos = 0;
	TSize maxTmp = 0;
	TSize numCol = 0;
	reserve(slot, alignItEnd - alignIt);
	for(;alignIt != alignItEnd; ++alignIt) {
		TPosIter itPos = begin(freePos, Standard());
		TPosIter itPosEnd = end(freePos, Standard());
		pos = 0;
		for(;itPos != itPosEnd; ++itPos, ++pos) 
			if (*itPos < _min(alignIt->beginPos, alignIt->endPos)) break;
		if (pos + 1 > length(freePos)) resize(freePos, pos+1, Generous());
		maxTmp = _max(alignIt->beginPos, alignIt->endPos);
		freePos[pos] = maxTmp;
		if (maxTmp > numCol) numCol = maxTmp;
		appendValue(slot, pos);
	}
	coverage = length(freePos);
	clear(freePos);

	// Fill the matrix
	typedef typename Iterator<TMatrix, Standard>::Type TMatIter;
	fill(mat, coverage * numCol, '.');
	alignIt = alignItBegin;
	TSize readPos = 0;
	TMatIter matIt = begin(mat, Standard());
	typename TFragmentStore::TReadSeq myRead;
	for(;alignIt != alignItEnd; ++alignIt, ++readPos) {
		typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor>, Standard>::Type TReadGapsIter;
		TReadGapsIter itGaps = begin(alignIt->gaps, Standard());
		TReadGapsIter itGapsEnd = end(alignIt->gaps, Standard());
		
		// Place each read inside the matrix
		myRead = fragStore.readSeqStore[alignIt->readId];
		TSize lenRead = length(myRead);
		TSize offset = alignIt->beginPos;
		if (alignIt->beginPos > alignIt->endPos) {
			reverseComplementInPlace(myRead);
			offset = alignIt->endPos;
		}
		matIt = begin(mat, Standard());
		matIt += (slot[readPos] * numCol + offset);

		typedef typename Iterator<typename TFragmentStore::TReadSeq, Standard>::Type TReadIter;
		TReadIter seqReadIt = begin(myRead, Standard());

		// First clear range
		TSize mySeqPos = 0;
		int diff = 0;
		if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) {
			mySeqPos = itGaps->seqPos;
			diff = -1 * mySeqPos;
			seqReadIt += mySeqPos;
		}
		TSize clr2 = lenRead;
		TSize stop = 0;
		for(;itGaps != itGapsEnd; ++itGaps) {
			// Any clipped sequence at the end
			stop =  itGaps->seqPos;
			if (diff - ((int) itGaps->gapPos - (int) itGaps->seqPos) > 0) 
				clr2 = stop = lenRead - (diff - ((int) itGaps->gapPos - (int) itGaps->seqPos));
			
			for(;mySeqPos < stop; ++matIt, ++seqReadIt, ++mySeqPos) 
				*matIt = *seqReadIt;

			for(int i = 0; i < ((int) itGaps->gapPos - (int) itGaps->seqPos) - diff; ++i, ++matIt) 
				*matIt = gapChar;
	
			diff = (itGaps->gapPos - itGaps->seqPos);
		}
		for(;mySeqPos < clr2; ++mySeqPos, ++seqReadIt, ++matIt) 
			*matIt = *seqReadIt;
	}
	//for(TSize row = 0; row < coverage; ++row) {
	//	for(TSize col = 0; col<numCol; ++col) {
	//		std::cout << mat[row * numCol + col];
	//	}
	//	std::cout << std::endl;
	//}
	return true;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TMatrix, typename TSize2, typename TSize> 
inline bool
convertAlignment(FragmentStore<TSpec, TConfig>& fragStore,
				 TMatrix& mat,
				 TSize2 contigId,
				 TSize& coverage)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	String<TSize> slot;
	return convertAlignment(fragStore, mat, contigId, coverage, slot);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TMatrix> 
inline bool
convertAlignment(FragmentStore<TSpec, TConfig>& fragStore,
				 TMatrix& mat)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	TSize coverage;
	return convertAlignment(fragStore, mat, 0, coverage);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TConfig, typename TGappedConsensus, typename TSize> 
inline void
_getGappedConsensusSeq(FragmentStore<TSpec, TConfig>& fragStore,
					   TGappedConsensus& gappedConsensus,
					   TSize contigId)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Value<TGappedConsensus>::Type TValue;
	
	TValue gapChar = gapValue<TValue>();
	typedef typename Iterator<typename TFragmentStore::TContigSeq, Standard>::Type TContigIter;
	TContigIter seqContigIt = begin(fragStore.contigStore[contigId].seq, Standard());
	TContigIter seqContigItEnd = end(fragStore.contigStore[contigId].seq, Standard());
	typedef typename Iterator<String<typename TFragmentStore::TContigGapAnchor>, Standard>::Type TGapsIter;
	TGapsIter itGaps = begin(fragStore.contigStore[contigId].gaps, Standard());
	TGapsIter itGapsEnd = end(fragStore.contigStore[contigId].gaps, Standard());
	int diff = 0;
	TSize mySeqPos = 0;
	for(;itGaps != itGapsEnd; goNext(itGaps)) {
		for(;mySeqPos < itGaps->seqPos; ++seqContigIt, ++mySeqPos) 
			appendValue(gappedConsensus, *seqContigIt, Generous());
			
		for(int i = 0; i < ((int) itGaps->gapPos - (int) itGaps->seqPos) - diff; ++i) 
			appendValue(gappedConsensus, gapChar, Generous());
			diff = (itGaps->gapPos - itGaps->seqPos);
	}
	for(;seqContigIt != seqContigItEnd; ++seqContigIt) 
		appendValue(gappedConsensus, *seqContigIt, Generous());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec, typename TConfig>
inline void 
write(TFile & file,
	  FragmentStore<TSpec, TConfig>& fragStore,
	  FastaReadFormat) 
{
	SEQAN_CHECKPOINT
	// Basic types
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef char TMultiReadChar;
	TMultiReadChar gapChar = gapValue<TMultiReadChar>();

	typedef typename Iterator<typename TFragmentStore::TContigStore, Standard>::Type TContigIter;
	TContigIter contigIt = begin(fragStore.contigStore, Standard() );
	TContigIter contigItEnd = end(fragStore.contigStore, Standard() );
	for(TSize idCount = 0;contigIt != contigItEnd; ++contigIt, ++idCount) {
		// Alignment matrix
		typedef String<TMultiReadChar> TAlignMat;
		TAlignMat mat;
		TSize maxCoverage;
		String<TSize> readSlot;
		convertAlignment(fragStore, mat, idCount, maxCoverage, readSlot);
		TSize len = length(mat) / maxCoverage;
		
		// Gapped consensus sequence
		typedef String<TMultiReadChar> TGappedConsensus;
		TGappedConsensus gappedConsensus;
		_getGappedConsensusSeq(fragStore, gappedConsensus, idCount);

		// Print the alignment matrix
		String<TSize> coverage;
		fill(coverage, len, 0);
		typedef typename Iterator<TGappedConsensus, Standard>::Type TConsIter;
		TConsIter itCons = begin(gappedConsensus, Standard());
		TSize winSize = 10000000;
		int offset = 2;
		TSize column = 0;
		while (column<len) {
			TSize window_end = column + winSize;
			if (window_end >= len) window_end = len;
			// Position
			for(int i = 0; i<offset - 2; ++i) _streamPut(file,' ');
			_streamWrite(file,"Pos: ");
			_streamPutInt(file, column);
			_streamPut(file,'\n');
			// Ruler
			for(int i = 0; i<offset + 3; ++i) _streamPut(file,' ');
			for(TSize local_col = 1; local_col<window_end - column + 1; ++local_col) {
				if ((local_col % 10)==0) _streamPut(file, ':');
				else if ((local_col % 5)==0) _streamPut(file, '.');
				else _streamPut(file, ' ');
			}
			_streamPut(file,'\n');
			// Matrix
			for(TSize row = 0; row<maxCoverage; ++row) {
				TSize tmp = row;
				int off = 0;
				while (tmp / 10 != 0) {
					tmp /= 10;
					++off;
				}
				for(int i = 0; i<offset - off; ++i) _streamPut(file,' ');
				_streamPutInt(file, row);
				_streamPut(file,':');
				_streamPut(file,' ');
				for(TSize local_col = column; local_col<window_end; ++local_col) {
					_streamPut(file, mat[row * len + local_col]);
					if (mat[row * len + local_col] != '.') ++coverage[local_col];
				}
				_streamPut(file,'\n');
			}
			_streamPut(file,'\n');
	
			// Consensus
			for(int i = 0; i<offset; ++i) _streamPut(file,' ');
			_streamWrite(file,"C: ");
			for(unsigned int local_col = column; local_col<window_end; ++local_col, ++itCons) 
				_streamPut(file, *itCons);
			_streamPut(file,'\n');
			for(int i = 0; i<offset-1; ++i) _streamPut(file,' ');
			_streamWrite(file,">2: ");
			for(unsigned int local_col = column; local_col<window_end; ++local_col) {
				if (coverage[local_col] > 2) _streamPut(file, gappedConsensus[local_col]);
				else _streamPut(file, gapChar);
			}
			_streamPut(file,'\n');
			_streamPut(file,'\n');
			column+=winSize;
		}
		_streamPut(file,'\n');
		_streamPut(file,'\n');


		// Print all aligned reads belonging to this contig

		// Sort according to contigId
		sortAlignedReads(fragStore.alignedReadStore, SortContigId());
	
		// Find range of the given contig
		typedef typename Iterator<typename TFragmentStore::TAlignedReadStore, Standard>::Type TAlignIter;
		TAlignIter alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigId());
		TAlignIter alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigId());

		// Sort the reads according to the begin position
		sortAlignedReads(infix(fragStore.alignedReadStore, alignIt - begin(fragStore.alignedReadStore, Standard()), alignItEnd - begin(fragStore.alignedReadStore, Standard())), SortBeginPos());
		TAlignIter alignItTmp = lowerBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigId());
		TAlignIter alignItTmpEnd = upperBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigId());
		String<std::pair<TSize, TSize> > idToPos;
		reserve(idToPos, alignItTmpEnd - alignItTmp);
		for(TSize iCount = 0; alignItTmp!=alignItTmpEnd; ++iCount, ++alignItTmp) 
			appendValue(idToPos, std::make_pair(alignItTmp->id, readSlot[iCount]));
		::std::sort(begin(idToPos, Standard()), end(idToPos, Standard()));

		// Sort the reads according to the id
		sortAlignedReads(infix(fragStore.alignedReadStore, alignIt - begin(fragStore.alignedReadStore, Standard()), alignItEnd - begin(fragStore.alignedReadStore, Standard())), SortId());
		alignIt = lowerBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigId());
		alignItEnd = upperBoundAlignedReads(fragStore.alignedReadStore, idCount, SortContigId());

		bool noNamesPresent = (length(fragStore.readNameStore) == 0);
		for(TSize iCount = 0;alignIt != alignItEnd; ++alignIt, ++iCount) {

			// Print all reads
			_streamWrite(file,"typ:");
			if (!noNamesPresent) {
				_streamPut(file,'R');
				_streamPutInt(file, iCount);
			} else _streamWrite(file, fragStore.readNameStore[alignIt->readId]);
			_streamPut(file,'\n');
			_streamWrite(file,"seq:");
			_streamWrite(file, fragStore.readSeqStore[alignIt->readId]);
			_streamPut(file,'\n');
			_streamWrite(file,"Pos:");
			_streamPutInt(file, alignIt->beginPos);
			_streamPut(file,',');
			_streamPutInt(file, alignIt->endPos);
			_streamPut(file,'\n');

			std::stringstream gapCoords;
			TSize letterCount = 0;
			TSize gapCount = 0;
			for(TSize column = _min(alignIt->beginPos, alignIt->endPos); column < _max(alignIt->beginPos, alignIt->endPos); ++column) {
				if (mat[idToPos[iCount].second * len + column] == gapChar) {
					++gapCount;
					gapCoords << letterCount << ' ';
				} else ++letterCount;
			}
			_streamWrite(file,"dln:");
			_streamPutInt(file, gapCount);
			_streamPut(file,'\n');
			_streamWrite(file,"del:");
			_streamWrite(file, gapCoords.str().c_str());
			_streamPut(file,'\n');
			_streamPut(file,'\n');
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// Rudimentary write functions for CeleraFrg and Celera Cgb
//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec, typename TConfig>
inline void 
_writeCeleraFrg(TFile& target,
				FragmentStore<TSpec, TConfig>& fragStore) 
{

	SEQAN_CHECKPOINT
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename TFragmentStore::TReadPos TReadPos;

	// Iterate over all aligned reads to get the clear ranges
	typedef Pair<TReadPos, TReadPos> TClrRange;
	String<TClrRange> clearStr;
	resize(clearStr, length(fragStore.readStore));
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		TReadPos begClr = 0;
		TReadPos endClr = 0;
		getClrRange(fragStore, value(alignIt), begClr, endClr);
		value(clearStr, alignIt->readId) = TClrRange(begClr, endClr);
	}

	// Write Reads
	typedef typename Iterator<typename TFragmentStore::TReadStore>::Type TReadIter;
	TReadIter readIt = begin(fragStore.readStore);
	TReadIter readItEnd = end(fragStore.readStore);
	bool noNamesPresent = (length(fragStore.readNameStore) == 0);
	for(TSize idCount = 0;readIt != readItEnd; goNext(readIt), ++idCount) {
		_streamWrite(target,"{FRG\n");
		_streamWrite(target,"act:");
		_streamPut(target, 'A');
		_streamPut(target, '\n');
		_streamWrite(target,"acc:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		_streamWrite(target,"typ:");
		_streamPut(target, 'R');
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"src:\n");
			_streamWrite(target, value(fragStore.readNameStore, idCount));
			_streamWrite(target, "\n.\n");
		}
		_streamWrite(target,"etm:");
		_streamPut(target, '0');
		_streamPut(target, '\n');
		_streamWrite(target,"seq:\n");
		typedef typename Iterator<typename TFragmentStore::TReadSeq>::Type TSeqIter;
		typedef typename Value<typename TFragmentStore::TReadSeq>::Type TAlphabet;
		TSeqIter seqIt = begin(value(fragStore.readSeqStore, idCount));
		TSeqIter seqItEnd = end(value(fragStore.readSeqStore, idCount));
		for(TSize k = 0;seqIt!=seqItEnd;goNext(seqIt), ++k) {
			if ((k % 70 == 0) && (k != 0)) _streamPut(target, '\n');
			_streamPut(target, getValue(seqIt));
		}
		_streamWrite(target, "\n.\n");
		_streamWrite(target,"qlt:\n");
		seqIt = begin(value(fragStore.readSeqStore, idCount));
		for(TSize k = 0;seqIt!=seqItEnd;goNext(seqIt), ++k) {
			if ((k % 70 == 0) && (k != 0)) _streamPut(target, '\n');
			Ascii c = ' ';
			convertQuality(c, getQualityValue(value(seqIt)));
			_streamPut(target, c);
		}
		_streamWrite(target, "\n.\n");
		// Note: Clear range does not have to be ordered, e.g. no indication for reverse complemented reads, this is happening in cgb records
		_streamWrite(target,"clr:");
		_streamPutInt(target, (value(clearStr, idCount)).i1);
		_streamPut(target, ',');
		_streamPutInt(target, (value(clearStr, idCount)).i2);
		_streamPut(target, '\n');
		_streamWrite(target,"}\n");
	}
}


//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TSpec, typename TConfig>
inline void 
_writeCeleraCgb(TFile& target,
				FragmentStore<TSpec, TConfig>& fragStore) 
{
	SEQAN_CHECKPOINT
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Id<TFragmentStore>::Type TId;
	typedef typename TFragmentStore::TReadPos TReadPos;


	// Write the first contig
	TId contigId = 0;

	// Sort the reads according to position
	sortAlignedReads(fragStore.alignedReadStore, SortBeginPos());

	// Write Header
	_streamWrite(target,"{IUM\nacc:0\nsrc:\ngen> @@ [0,0]\n.\ncov:0.000\nsta:X\nfur:X\nabp:0\nbbp:0\n");
	_streamWrite(target,"len:");
	_streamPutInt(target, length((value(fragStore.contigStore, contigId)).seq));
	_streamPut(target, '\n');
	_streamWrite(target,"cns:\n.\nqlt:\n.\nfor:0\n");
	_streamWrite(target,"nfr:");
	_streamPutInt(target, length(fragStore.readStore));
	_streamPut(target, '\n');

	// Write reads
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	TSize offsetLeft = _min(alignIt->beginPos, alignIt->endPos);
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		if (contigId != alignIt->contigId) continue;
		_streamWrite(target,"{IMP\n");
		_streamWrite(target,"typ:");
		_streamPut(target, 'R');
		_streamPut(target, '\n');
		_streamWrite(target,"mid:");
		_streamPutInt(target, alignIt->readId + 1);
		_streamPut(target, '\n');
		_streamWrite(target,"con:");
		_streamPut(target, '0');
		_streamPut(target, '\n');
		_streamWrite(target,"pos:");
		_streamPutInt(target, alignIt->beginPos - offsetLeft);
		_streamPut(target, ',');
		_streamPutInt(target, alignIt->endPos - offsetLeft);
		_streamPut(target, '\n');
		_streamWrite(target,"dln:0\n");
		_streamWrite(target,"del:\n");
		_streamWrite(target,"}\n");
	}
	_streamWrite(target,"}\n");
}






}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
