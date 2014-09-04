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

#include <iostream>

#ifndef SEQAN_HEADER_STORE_IO_SAM_OUT_H
#define SEQAN_HEADER_STORE_IO_SAM_OUT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// write
    
    template<typename TFile, typename TFragmentStore>
    inline void write(TFile & target,
                      TFragmentStore const & store,
                      SAM)
    {
        // write header
        
        // write aligments
        _writeAlignments(target, store, SAM());
    }
    
//////////////////////////////////////////////////////////////////////////////
// _writeAlignments

    template<typename TFile, typename TFragmentStore>
    inline void _writeAlignments(TFile & target,
                                 TFragmentStore const & store,
                                 SAM)
    {
        typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;
        typedef typename Iterator<String<TAlignedElement> >::Type TAlignIter;
        typedef typename Id<TFragmentStore>::Type TId;

        TAlignIter it = begin(store.alignedReadStore);

        for(; it != end(store.alignedReadStore); ++it){
            TId alignedId = value(it).id;
            _streamWrite(target, value(store.readNameStore, alignedId));
            _streamPut(target, '\t');
            
            //flag
            _streamPut(target, '\t');
            
            _streamWrite(target, value(store.contigNameStore, value(it).contigId));
            _streamPut(target, '\t');
            
            _streamPutInt(target, min(value(it).beginPos, value(it).endPos));
            _streamPut(target, '\t');
            
            _streamPutInt(target, value(store.alignedReadQualityStore, alignedId));
            _streamPut(target, '\t');
            
            // cigar
            _streamPut(target, '\t');
            
            // mrnm
            _streamPut(target, '\t');
            
            // mpos
            _streamPut(target, '\t');
            
            // isize
            _streamPut(target, '\t');
            
            _streamWrite(target, value(store.readSeqStore, value(it).readId));
            _streamPut(target, '\t');
            
            // qual
            _streamPut(target, '\t');
            
            _streamWrite(target, value(store.alignedReadTagStore, alignedId));
            
            _streamPut(target, '\n');
        }
        
    }
    
    
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
