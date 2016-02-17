bin_PROGRAMS = pacbiokanga
pacbiokanga_SOURCES= SQLiteSummaries.cpp SQLiteSummaries.h SSW.cpp SSW.h SWAlign.cpp SWAlign.h PBAssemb.cpp PBAssemb.h PBECContigs.cpp PBECContigs.h \
                     SeqStore.cpp SeqStore.h PBFilter.cpp PBFilter.h pacbiocommon.h PacBioUtility.cpp PacBioUtility.h pacbiokanga.cpp pacbiokanga.h \
                     PBSim.cpp PBSim.h PBErrCorrect.cpp PBErrCorrect.h MAConsensus.cpp MAConsensus.h AssembGraph.cpp AssembGraph.h \
                     MAFKMerDist.cpp MAFKMerDist.h PBSWService.cpp PBSWService.h BKSProvider.cpp BKSProvider.h BKSRequester.cpp BKSRequester.h BKScommon.h 

# set the include path found by configure
INCLUDES= $(all_includes)

# our library
LDADD = ../libbiokanga/libbiokanga.a ../libbiokanga/zlib/libz.a ../libBKPLPlot/libBKPLPlot.a

