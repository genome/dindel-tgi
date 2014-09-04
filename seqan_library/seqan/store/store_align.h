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

#ifndef SEQAN_HEADER_STORE_ALIGN_H
#define SEQAN_HEADER_STORE_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Aligned Read Store
//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec = void>
struct AlignedReadStoreElement
{
	typedef typename Id<AlignedReadStoreElement>::Type TId;

	static const TId INVALID_ID;
	
	TId					id;
	TId					readId;
	TId					contigId;
	TId					pairMatchId;	// unique id. for multiple mate-pair matches (not matePairId)
	TPos				beginPos;		// begin position of the gapped sequence in gapped contig sequence
	TPos				endPos;			// end position of ..., for reverse aligned reads holds end < begin
	String<TGapAnchor>	gaps;

	AlignedReadStoreElement() : id(INVALID_ID), readId(INVALID_ID), contigId(INVALID_ID), pairMatchId(INVALID_ID), beginPos(0), endPos(0) {}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec>
const typename Id<AlignedReadStoreElement<TPos, TGapAnchor, TSpec> >::Type
AlignedReadStoreElement<TPos, TGapAnchor, TSpec>::INVALID_ID = SupremumValue<typename Id<AlignedReadStoreElement<TPos, TGapAnchor, TSpec> >::Type>::VALUE;


//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Sorting tags
//////////////////////////////////////////////////////////////////////////////


struct _SortContigId;
typedef Tag<_SortContigId> const SortContigId;


struct _SortId;
typedef Tag<_SortId> const SortId;

struct _SortBeginPos;
typedef Tag<_SortBeginPos> const SortBeginPos;

struct _SortEndPos;
typedef Tag<_SortEndPos> const SortEndPos;

struct _SortPairMatchId;
typedef Tag<_SortPairMatchId> const SortPairMatchId;

struct _SortReadId;
typedef Tag<_SortReadId> const SortReadId;


//////////////////////////////////////////////////////////////////////////////
// Sorting functors
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template <typename TAlignedRead, typename TTag>
struct _LessAlignedRead;

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return (a1.id) < (a2.id);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortContigId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.contigId < a2.contigId;
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortBeginPos> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return _min(a1.beginPos, a1.endPos) < _min(a2.beginPos, a2.endPos);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortEndPos> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return _max(a1.beginPos, a1.endPos) < _max(a2.beginPos, a2.endPos);
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortPairMatchId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.pairMatchId < a2.pairMatchId;
	}
};

template <typename TAlignedRead>
struct _LessAlignedRead<TAlignedRead, SortReadId> :
	public ::std::binary_function<TAlignedRead, TAlignedRead, bool>
{
	inline bool 
	operator() (TAlignedRead const& a1, TAlignedRead const& a2) const {
		return a1.readId < a2.readId;
	}
};

//////////////////////////////////////////////////////////////////////////////
// Sorting function
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign& alignStore, Tag<TSortSpec>) 
{
	std::stable_sort(
		begin(alignStore, Standard() ), 
		end(alignStore, Standard() ), 
		_LessAlignedRead<typename Value<TAlign>::Type, Tag<TSortSpec> const>() );
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSortSpec>
inline void
sortAlignedReads(TAlign const& alignStore, Tag<TSortSpec>) 
{
	std::stable_sort(
		begin(const_cast<TAlign&>(alignStore), Standard() ), 
		end(const_cast<TAlign&>(alignStore), Standard() ), 
		_LessAlignedRead<typename Value<TAlign>::Type, Tag<TSortSpec> const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.id = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortContigId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.contigId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortContigId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortBeginPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.beginPos = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortBeginPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.endPos = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortEndPos) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.endPos = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortEndPos const>() );
}


//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortPairMatchId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.pairMatchId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortPairMatchId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
lowerBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::lower_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadId const>() );
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSearchValue>
inline typename Iterator<TAlign, Standard>::Type
upperBoundAlignedReads(TAlign const& alignStore, 
					   TSearchValue const val,
					   SortReadId) 
{
	typedef typename Value<TAlign>::Type TAlignElement;
	TAlignElement el;
	el.readId = val;
	return ::std::upper_bound(
		begin(alignStore, Standard()), 
		end(alignStore, Standard()), 
		el,
		_LessAlignedRead<typename Value<TAlign>::Type, SortReadId const>() );
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
