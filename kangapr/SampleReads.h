#pragma once
const int cAllocRawSeqLen = 50000; // raw reads are expected to be less than this length
typedef enum TAG_ePMode {
	ePMDefault = 0,				// default processing mode is for single ended reads
	ePMPairedEnd,				// paired end processing with paired end reads in separate P1 and P2 files
	ePMplaceholder				// used to set the enumeration range
} etPMode;

class CSampleReads
{
	bool m_bInitialised;		// true after instantiation initialisation completed
	etPMode m_PMode;			// process mode, single or paired end
	bool m_bSOLiD;				// if true then input reads are expected to be in SOLiD colour space
	bool m_bIsFastq;			// if basespace then input reads are fastq
	bool m_bAppendOut;			// if true then append if output files exist, otherwise trunctate existing output files 

	int m_SampleOfs;				// start at this read instance (1 = first read)
	int m_SampleIncr;				// sample every Nth from SampleOfs
	int m_SampleMax;				// sample at most this number of reads

	CFasta m_P1InFasta;			// used for parsing in P1 reads
	CFasta m_P2InFasta;			// used for parsing in optional P2 reads
	int m_hOut5Reads;			// file handle for output 5' sampled reads
	int m_hOut3Reads;			// file handle for output 3' sampled reads

	char m_szIn5ReadsFile[_MAX_PATH];	// 5' reads are in this file
	char m_szIn3ReadsFile[_MAX_PATH];	// 3' reads are in this file
	char m_szOut5File[_MAX_PATH];	// write 5' sampled reads to this file
	char m_szOut3File[_MAX_PATH];	// write 3' sampled reads to this file

public:
	CSampleReads(void);
	~CSampleReads(void);
	void Reset(bool bSync = true);   // if bDync true then opened output files will be sync'd (commited) to disk before closing
	int
		GenSamples(etPMode PMode,		// process mode, single or paired end
				bool bAppendOut,			// if true then append if output files exist, otherwise trunctate existing output files 
			  int SampleOfs,				// start at this read instance (1 = first read)
			  int SampleIncr,				// sample every Nth from SampleOfs
			  int SampleMax,				// sample at most this number of reads	  
			  char *pszIn5ReadsFile,		// 5' reads are in this file
			  char *pszIn3ReadsFile,		// 3' reads are in this file
			  char *pszOut5ReadsFile,		// write sampled 5' reads to this file
			  char *pszOut3ReadsFile);		// write sampled 5' reads to this file


	int
		ProcOverlapPairs(void);				// overlap paired end processing

	int 
		OpenFiles(void);					// open files (as set by MergeOverlaps) for processing
};

