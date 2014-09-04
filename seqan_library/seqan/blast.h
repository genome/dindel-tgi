#ifndef SEQAN_HEADER_BLAST_H
#define SEQAN_HEADER_BLAST_H

//____________________________________________________________________________
// prerequisites

#include <iostream>
#include <fstream>

#include <seqan/align.h>
#include <seqan/graph_types.h>
#include <seqan/misc/misc_parsing.h>
#include <seqan/graph_algorithms.h>


#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/blast/blast_generated_forwards.h>
#endif


#include <seqan/blast/blast_base.h>
#include <seqan/blast/blast_hsp.h>
#include <seqan/blast/blast_hit.h>
#include <seqan/blast/blast_report.h>
#include <seqan/blast/blast_parsing.h>
#include <seqan/blast/blast_iterator.h>
#include <seqan/blast/blast_hit_iterator.h>
#include <seqan/blast/blast_hsp_iterator.h>


#include <seqan/blast/blast_stream_report.h>
#include <seqan/blast/blast_stream_hit.h>
#include <seqan/blast/blast_stream_hit_iterator.h>
#include <seqan/blast/blast_stream_hsp_iterator.h>


#include <seqan/blast/blast_run.h>


#endif //#ifndef SEQAN_HEADER_...
