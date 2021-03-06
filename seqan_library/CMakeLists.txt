cmake_minimum_required (VERSION 2.6) 
project (SeqAn)

set (CMAKE_BUILD_TYPE Debug)
#set (CMAKE_BUILD_TYPE Release)

################################################################################
# Set path variables
################################################################################

set (SEQAN_BASE .)

################################################################################
# Library sources
################################################################################

file(GLOB SEQAN_SOURCES_ALIGN ${SEQAN_BASE}/seqan/align/*.h)
#file(GLOB_RECURSE SEQAN_SOURCES_BASIC ${SEQAN_BASE}/seqan FOLLOW_SYMLINKS ${SEQAN_BASE}/seqan/basic/*.h)
#file(GLOB_RECURSE SEQAN_SOURCES_GRAPH ${SEQAN_BASE}/seqan FOLLOW_SYMLINKS ${SEQAN_BASE}/seqan/graph/*.h)
#file(GLOB_RECURSE SEQAN_SOURCES_INDEX ${SEQAN_BASE}/seqan FOLLOW_SYMLINKS ${SEQAN_BASE}/seqan/index/*.h)
#file(GLOB_RECURSE SEQAN_SOURCES_SEQUENCE ${SEQAN_BASE}/seqan FOLLOW_SYMLINKS ${SEQAN_BASE}/seqan/sequence/*.h)
#add_executable (razers2 razers.cpp  ${SEQAN_SOURCES})

################################################################################
# Demos
################################################################################

aux_source_directory (${SEQAN_BASE}/demos SEQAN_DEMOS)
add_subdirectory (demos)
foreach (DEMOFILE ${SEQAN_DEMOS})
	get_filename_component (DEMO ${DEMOFILE} NAME_WE)
	add_executable (demos/${DEMO} ${DEMOFILE})
endforeach (DEMOFILE)

source_group("Header Files\\align" FILES ${SEQAN_SOURCES_ALIGN})

# Make sure the compiler can find include files from the SeqAn and TBB library. 
include_directories (${SEQAN_BASE})

# Make sure the compiler can find object files from the TBB library. 
link_directories (${SEQAN_BASE}/lib)

set (RAZERS_PREFIX apps/razers/)
set (HEADERS ${RAZERS_PREFIX}recognitionRateDP.h ${RAZERS_PREFIX}readSimulator.h ${RAZERS_PREFIX}razers.h ${RAZERS_PREFIX}razers_matepairs.h ${RAZERS_PREFIX}paramChooser.h)
set (RAZERS_HEADERS ${RAZERS_PREFIX}razers_parallel.h ${RAZERS_PREFIX}outputFormat.h ${RAZERS_PREFIX}mmap_fasta.h ${HEADERS})
list (SORT HEADERS)
list (SORT RAZERS_HEADERS)

add_executable (dfi apps/dfi/dfi.cpp)
add_executable (paramChooser apps/razers/paramChooser.cpp ${HEADERS})
add_executable (razers apps/razers/razers.cpp ${RAZERS_HEADERS})
add_executable (seqan_tcoffee apps/seqan_tcoffee/seqan_tcoffee.cpp)
add_executable (seqcons apps/seqcons/seqcons.cpp)

#target_link_libraries (paramChooser rt)
#target_link_libraries (razers tbb_debug pthread)
#target_link_libraries (razers tbb_debug)
target_link_libraries (razers)

################################################################################
# Set compiler flags
################################################################################

if (CMAKE_COMPILER_IS_GNUCXX)
	set (CMAKE_CXX_WARNING_LEVEL 4)
endif (CMAKE_COMPILER_IS_GNUCXX)

if (MSVC)
	add_definitions (/W2 /wd4996)
endif (MSVC)
