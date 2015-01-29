#pragma once

const unsigned int cRRMinSeqLen = 10;			// sequences must be at least this length otherwise user is warned and sequence sloughed

const int cBSFRRRdsVersion = 6;				// current file header version
const int cBSFRRRdsVersionBack= 5;			// backward compatible to this version

const int cRRMaxInFiles = 100;				// allow for at most this many input files
const int cRRMaxInFileSpecs = 100;			// allow for at most this many input file specs which may be wildcarded
const int cRRWrtBuffSize = 0x01ffffff;		// use a 32 MB write buffer size (must be at least cMaxFastQSeqLen + 1000)
const int cRRMaxDescrLen = 128;				// allow for descriptors of upto this length

const int cRRDataBuffAlloc = 0x0fffffff;		// alloc to hold reads in this byte sized increments, initial allocation will be 2x
const int cRRRdsBuffAlloc = 0x0fffff;			// alloc to hold preprocessed reads (for stats) in this byte sized allocation

// processing modes
typedef enum TAG_ePRRMode {		
	ePMRRNewSingle,				// generate new processed single end reads file
	ePMRRNewPaired,				// generate new processed paired end reads file
	ePMRRStats,					// generate statistics on processed reads
	ePMRRDumpFasta,				// dump rds file as fasta
	ePMRRplaceholder			// used to set the enumeration range
	} etPRRMode;

// fastq quality scoring method
typedef enum TAG_eFQMethod {
	eFQSanger,					// Sanger or Illumina 1.8+ Phred
	eFQIllumia,					// Illumina 1.3+
	eFQSolexa,					// Solexa/Illumia pre-1.3
	eFQIgnore,					// ignore fastq quality scores
	eFQplaceholder				// used to set the enumeration range
	} etFQMethod;

class CProcRawReads  : public CErrorCodes
{

	UINT32 m_NumDescrReads;					// number of tsDescrRead's
	size_t m_DataBuffAllocMem;				// total memory allocated 
	size_t m_DataBuffOfs;					// offset at which to next write
	UINT8 *m_pDataBuff;							// NOTE - memory for holding reads is allocated and reallocated using C's alloc() and realloc() not c++'s new()!!!!!
												// reason for this is that it reduces memory requirements when resizing this buffer

	tsBSFRdsHdr m_FileHdr;						// processed reads file header

	UINT32 *m_pDimerCnts;						// counts of all dimers
	UINT32 *m_pTrimerCnts;						// counts of all trimers
	UINT32 *m_pTetramerCnts;					// counts of all tetramers
	UINT32 *m_pDistQualScores;					// quality score distributions
	UINT8 *m_pWrtBuff;							// used to buffer output writes to file
	UINT8 *m_pRdsBuff;							// used to buffer reads

	size_t m_ReadsIdxBuffAllocMem;				// size of memory allocated for holding the reads index
	tsRawReadV6 **m_ppReadsIdx;					// reads index

	bool m_bIsSOLiD;							// true if any SOLiD or colorspace reads file processed

	bool m_bActivity;							// true if processing activity to be reported to user


	int m_hInFile;								// input file handle
	int m_hOutFile;								// output results file

	teBSFrsltCodes Hdr2Disk(char *pszRdsFile);  // write header to disk
	teBSFrsltCodes Disk2Hdr(char *pszRdsFile);  // load header from disk

	int
		AddEntry(bool bIsPairRead,		// true if this is the paired read
				 UINT32 PairReadID,		// identifies partner of this read if paired read processing
				 UINT8 FileID,			// identifies file from which this read was parsed
				 int DescrLen,			// length of following descriptor
				 char *pszReadDescr,	// copy of descriptor upto first space char, used to pair reads with matching descriptors
				 int ReadLen,			// length of following read
				 UINT8 *pszReadBuff);	// packed read + phred score

	int CompareRead(tsRawReadV6 *pRead1,tsRawReadV6 *pRead2);

	CMTqsort m_MTqsort;					// multi-threaded qsort

	static int SortReads(const void *arg1, const void *arg2);
	static int SortReadsPair(const void *arg1, const void *arg2);


public:
	CProcRawReads(void);
	~CProcRawReads(void);

	void Reset(void);

	void ReportActivity(bool bActivity);

	teBSFrsltCodes
	LoadAndProcessReads(etPRRMode PMode,		// processing mode
		UINT32 NumReadsLimit,				// limit processing to this many reads
		bool bKeepDups,						// true if duplicate reads not to be filtered out
		etFQMethod Quality,					// fastq quality value method
		int Trim5,							// trim this many bases off leading sequence 5' end
		int Trim3,							// trim this many bases off trailing sequence 3' end
		int NumInputFileSpecs,				// number of input file specs 
		char *pszInfileSpecs[],				// names of inputs file (wildcards allowed unless in dump mode) containing raw reads
		char *pszInPairFile,				// if paired reads processing then file containing paired ends
		char *pszOutFile);					// output into this file, NULL or '\0' if reads to be parsed etc but not written to file

	teBSFrsltCodes
		LoadAndProcessReadsDE(etPRRMode PMode,						// processing mode
		UINT32 NumReadsLimit,				// limit processing to this many reads
		int Trim5,							// trim this many bases off leading sequence 5' end
		int Trim3,							// trim this many bases off trailing sequence 3' end
		int	MinSampleCnts,					// minimum sample counts
		int	MinTotalCnts,					// minimum total samples counts
		int NumInputFileSpecs,				// number of input file specs 
		char *pszInfileSpecs[],				// names of inputs file (wildcards allowed unless in dump mode) containing raw reads
		char *pszOutFile);					// output into this file

	teBSFrsltCodes
		WriteToFile(etPRRMode PMode,			// processing mode
			etFQMethod Quality,				// fastq quality value method
			bool bKeepDups,					// true if duplicate reads were not filtered out
			int Trim5,						// trimmed this many bases off leading sequence 5' end
			int Trim3,						// trimmed this many bases off trailing sequence 3' end
			char *pszOutFile);				// output into this file

	teBSFrsltCodes
	GenStats(etPRRMode PMode,					// processing mode
		char *pszInfile,					// input file
		char *pszOutfile);					// output into this file

	teBSFrsltCodes
	GenDump(etPRRMode PMode,						// processing mode
   		int NumReadsLimit,					// limit output to this many reads
		char *pszInfile,					// input file
		char *pszOutfile);					// output into this file




	teBSFrsltCodes LoadReads(bool bIsPairRead,	// true if this file to process contains the paired reads
						UINT32 PairReadID,	// if non-zero then start paired reads identifiers from this value and increment after each read processed
						 UINT32 NumReadsLimit,				// limit to at most this number of reads, 0 if no limit
						 etFQMethod Quality,				// fastq quality value method
				 		 int Trim5,							// trim this many bases off leading sequence 5' end
						 int Trim3,							// trim this many bases off trailing sequence 3' end
    			 		 int FileID,						// uniquely identifies source file
						 char *pszFile);					// process from this file


	char *LocateFileName(int FileID);			// locates and returns filename corresponding to FileID
};

