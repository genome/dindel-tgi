#ifndef SEQAN_HEADER_BLAST_ITERATOR_H
#define SEQAN_HEADER_BLAST_ITERATOR_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Blast Iterators
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TSpec>
struct SimpleBlastIterator;	// for BlastReport<TBlastHsp,StoreBlast,TSpec>


// Hit Iterator
struct HitIterator;

// Hsp Iterator
struct HspIterator;


//////////////////////////////////////////////////////////////////////////////



template <typename TSpec>
struct StreamBlastIterator;	// for BlastReport<TBlastHsp,StreamBlast,TSpec>



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
