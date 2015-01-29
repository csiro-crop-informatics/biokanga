noinst_LIBRARIES = libbiokanga.a
libbiokanga_a_SOURCES = AlignValidate.cpp argtable2.cpp BEDfile.cpp BioSeqFile.cpp \
	Centroid.cpp Conformation.cpp ConfSW.cpp CSVFile.cpp CVS2BED.cpp DataPoints.cpp \
	Diagnostics.cpp Endian.cpp ErrorCodes.cpp Fasta.cpp FeatLoci.cpp \
	FilterLoci.cpp FilterRefIDs.cpp GOAssocs.cpp GOTerms.cpp \
	HashFile.cpp HyperEls.cpp GFFFile.cpp GTFFile.cpp GOAssocs.cpp GOTerms.cpp Contaminants.cpp \
	MAlignFile.cpp Random.cpp SimpleRNG.cpp RsltsFile.cpp sais.cpp SAMfile.cpp SeqTrans.cpp SfxArray.cpp SfxArrayV2.cpp Shuffle.cpp \
	SSW.cpp SSW.h SmithWaterman.cpp NeedlemanWunsch.cpp Stats.cpp StopWatch.cpp Twister.cpp Utility.cpp ProcRawReads.cpp MTqsort.cpp \
        bgzf.cpp sqlite3.c

# set the include path found by configure
INCLUDES= -I../hdrs -I@top_srcdir@/ -I$(all_includes)

# using zlib
LDADD = ./zlib/libz.a
