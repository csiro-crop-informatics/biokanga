bin_PROGRAMS = biokanga
biokanga_SOURCES= biokanga.cpp biokanga.h csv2sqlite.cpp SimReads.cpp Markers.cpp Markers.h SQLiteSummaries.cpp SQLiteSummaries.h SQLiteMarkers.cpp SQLiteMarkers.h \
                  SQLiteDE.cpp SQLiteDE.h psl2sqlite.cpp SQLitePSL.cpp SQLitePSL.h kanga.cpp kanga.h Aligner.cpp Aligner.h kangade.cpp Kangadna.cpp Kangadna.h \
                  FastaNxx.cpp FastaNxx.h kangax.cpp kangax.h genmarkerseq.cpp MarkerSeq.cpp MarkerSeq.h genDESeq.cpp genpseudogenome.cpp \
                  maploci2features.cpp MapLoci2Feat.cpp MapLoci2Feat.h \
		  mergeoverlaps.cpp MergeReadPairs.cpp MergeReadPairs.h fastaextract.cpp Assemble.cpp \
                  LocKMers.cpp LocKMers.h genkmarkers.cpp gensnpmarkers.cpp kmermarkers.cpp MarkerKMers.cpp MarkerKMers.h \
		  PEScaffold.cpp PEScaffold.h SSRdiscovery.cpp SSRdiscovery.h FilterSAMAlignments.cpp FilterSAMAlignments.h \
                  deNovoAssemb.cpp deNovoAssemb.h ArtefactReduce.cpp ArtefactReduce.h Scaffolder.cpp Scaffolder.h \
                  ReadStats.cpp ReadStats.h Blitz.cpp Blitz.h RemapLoci.cpp RemapLoci.h LocateROI.cpp LocateROI.h 

# set the include path found by configure
INCLUDES= $(all_includes)

# our library
LDADD = ../libbiokanga/libbiokanga.a ../libbiokanga/zlib/libz.a ../libBKPLPlot/libBKPLPlot.a

