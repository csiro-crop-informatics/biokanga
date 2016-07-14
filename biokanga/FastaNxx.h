#pragma once

const int cDfltMinLengthRange = 4;		// accept sequence lengths from this min length up
const int cWrtBuffSize = 0x0fffff;		// use a 1 MB write buffer size (must be at least cMaxFastQSeqLen + 100)
const int cMaxTruncLength = 2000;	    // max sequence length accepted, overlength sequences are truncated
const int cDfltTruncLength = 200;	    // default truncate length

const int cMaxN50Length = 100000000;	// if N50 processing then can handle at most this length contigs
const int cDfltN50Length = 100000000;	// if N50 processing then this is the default max length contig
const UINT32 cAllocCtgLenDist = 10000000;	// allocate to hold this number of contig lengths for calc of N50 etc

const int cMaxInFileSpecs = 50;			// allow at most this many input file specs

typedef enum {
	ePMdefault = 0,						// default is for Nxx calculations
	ePMKMerDist							// K-mer distributions
	} etNxxPMode;

class CFastaNxx
{
	size_t m_AllocMemSeq;				// allocation size of m_pSeq
	etSeqBase *m_pSeq;					// allocated to fasta sequences
	int m_NumContigLens;				// current number of contig lengths in m_pContigLengths
	size_t m_AllocdCtgLens;				// memory size currently alloc'd to m_pContigLengths 
	int *m_pContigLengths;				// used to hold contig lengths for N50 and other stats generation 
	INT64 m_TotLenCovered;				// sum of all contig lengths
	UINT32 m_BaseCnts[6];				// total counts for each base including 'N's
	UINT32 m_DistBaseCnts[6][cMaxFastQSeqLen]; // distribution of base counts along reads
	UINT32 m_SeqNsCnts[cMaxFastQSeqLen]; // distribution of number of 'N's over all reads 
	UINT32 *m_pDimerCnts;					// counts of all dimers
	UINT32 *m_pTrimerCnts;				// counts of all trimers
	UINT32 *m_pTetramerCnts;			// counts of all tetramers
	UINT32 *m_pBins;					// bin counts
	int m_MaxLengthRead;
	CMTqsort m_MTqsort;					// use this instance for sorts
	static int SortByContigLen(const void *arg1, const void *arg2); // compare length function used by m_MTqsort

public:
	CFastaNxx();
	~CFastaNxx();

	void Reset(void);

	int 
	Process(etNxxPMode Mode,				// processing mode - 0  N50 distributions, 1 k-mer count distributions
			int MinLength,					// core elements must be of at least this length
			int MaxLength,					// will be truncated to this length
			int NumBins,				// when generating length distributions then use this many bins - 0 defaults to using 1000
			int BinDelta,				// when generating length distributions then each bin holds this length delta - 0 defaults to auto determine from NunBins and longest sequence length
			int NumInputFiles,				// number of input sequence files
			char **pszInFastaFile,			// names of input sequence files (wildcards allowed)
			char *pszRsltsFile);			// file to write fasta into

	int ProcessFile(etNxxPMode Mode,	// processing mode - 0  N50 distributions
			int MinLength,				// core elements must be of at least this length
			int MaxLength,				// process sequences of upto this length
			char *pszFastaFile);
};

