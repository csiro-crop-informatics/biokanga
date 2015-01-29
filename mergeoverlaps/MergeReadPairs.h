#pragma once
// Merge read pairs where there is an expected 3' read overlap onto the 5' read

// output merged file format
typedef enum TAG_ePMode {
	ePMdefault = 0,				// Standard processing
	ePMcombined,				// combine overlaps and orphan reads into same file
	ePMseparate,					// overlaps and orphans into separate files
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// output merged file format
typedef enum TAG_eOFormat {
	eOFauto = 0,				// default is to merge and output in same format as inputs
	eOFfasta,					// merge and output as fasta
	eOFfastq,					// merge and output as fastq
	eOFplaceholder				// used to set the enumeration range
	} etOFormat;


const int cMinOverlapLen = 1;		// minimium allowed used specified number of overlap bases
const int cMaxOverlapPercSubs = 10;	// maximum allowed substitutions as a percentage of the overlap

const int cAllocOutBuffLen = 0x07ffffff; // 128M output buffering when writing to file

const char cPhredLowScore = 'A';	// generated Phred score used when there was a merge induced substitution in overlay and input sequences were non-fastq
const char cPhredOrphanScore = 'H';	// generated Phred score used if orphan reads and input sequences were non-fastq
const char cPhredHiScore = 'J';		// generated Phred score used when there was no merge substitution in overlay and input sequences were non-fastq

typedef enum TAG_eProcPhase {
	ePPUninit,					// uninitialised
	ePPReset,					// reset ready for file processing
	ePPRdyOpen,					// ready for input files to be opened 
	ePRdyProc					// files opened, ready for processing
} etProcPhase;

class CMergeReadPairs
{
	etProcPhase m_ProcPhase;	// current processing phase
	etPMode m_PMode;			// processing mode 
	etOFormat m_OFormat;		// output file format
	bool m_bIsFastq;			// if basespace then input reads are fastq
	bool m_bAppendOut;			// if true then append if output files exist, otherwise trunctate existing output files 
	int m_MinOverlap;			// reads must overlap by at least this many bases
    int m_MaxOverlapPropSubs;	// and overlap can have at most this proportion of substitutions
	int m_hIn5Reads;			// file handle for opened 5' reads
	int m_hIn3Reads;			// file handle for opened 3' reads
	int m_hOutMerged;			// file handle for output merged read microcontigs
	int m_hOut5Unmerged;		// file handle for output 5' unmerged reads
	int m_hOut3Unmerged;		// file handle for output 3' unmerged reads
	
	char m_szIn5ReadsFile[_MAX_PATH];	// 5' reads are in this file
	char m_szIn3ReadsFile[_MAX_PATH];	// 3' reads are in this file
	char m_szMergeOutFile[_MAX_PATH];	// write merged overlaping reads to this file
	char m_sz5UnmergedOutFile[_MAX_PATH];	// write 5' unmerged reads to this file
	char m_sz3UnmergedOutFile[_MAX_PATH];   // write 3' unmerged reads to this file

	int m_CurMSeqLen;			// m_pszMSeqs currently holds this many chars
	int m_AllocdMSeqs;			// m_pszMSeqs can hold at most this many chars
	char *m_pszMSeqs;			// will be allocated for holding merged sequences ready for writing to file
	int m_CurUnmergedP1Seqs;	// m_pszUnmergedP1Seqs currently holds this many chars
	int m_AllocdUnmergedP1Seqs;	// m_pszUnmergedP1Seqs can hold at most this many chars
	char *m_pszUnmergedP1Seqs;  // will be allocated for holding unmerged P1 sequences ready for writing to file
	int m_CurUnmergedP2Seqs;	// m_pszUnmergedP2Seqs currently holds this many chars
	int m_AllocdUnmergedP2Seqs;	// m_pszUnmergedP2Seqs can hold at most this many chars
	char *m_pszUnmergedP2Seqs;  // will be allocated for holding unmerged P2 sequences ready for writing to file

public:
	CMergeReadPairs(void);
	~CMergeReadPairs(void);
	void Reset(bool bSync = true);   // if bDync true then opened output files will be sync'd (commited) to disk before closing
	int									// returns number of sequences merged and written to output file
		MergeOverlaps(etPMode PMode,	// processing mode 
			  etOFormat OFormat,				// output file format
				bool bAppendOut,			// if true then append if output files exist, otherwise trunctate existing output files
			  int StartNum,					// initial starting sequence identifier
			  int MinOverlap,				// reads must overlap by at least this many bases
              int MaxOverlapPropSubs,		// and overlap can have at most this proportion of substitutions, if > 0 then floor of 1 sub allowed
			  char *pszIn5ReadsFile,		// 5' reads are in this file
			  char *pszIn3ReadsFile,		// 3' reads are in this file
			  char *pszMergeOutFile);		// write merged overlaping reads to this file

	int													// returns number of sequences merged and written to output file
		ProcOverlapPairs(int StartNum);		// initial starting sequence identifier

	int 
		OpenFiles(void);					// open files (as set by MergeOverlaps) for processing
};

