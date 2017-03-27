#pragma once

const int cMinKMerDistLen = 16;		// minimum allowed length KMer length processed when analysing Kmer distributions in read sequences	
const int cMaxKMerDistLen = 256;	// max allowed length KMer length 

const size_t cInitAllocKMerSeqInsts = 10000000;						// initially allocate to hold this many KMer sequence instances
const size_t cReallocKMerSeqInsts = (cInitAllocKMerSeqInsts/2);		// then extend as may be required by this many instances
const size_t cMaxKMerSeqInsts = 0xfffffff0;							// subject to available memory, can process at most this many unique KMer sequence instances

const UINT32 cKMerSeqHashMask =    0x0fffff;						// hashing mask (cMaxKMerSeqHashArrayEntries - 1)
const UINT32 cMaxKMerSeqHashArrayEntries = (cKMerSeqHashMask + 1);  // alloc hash array to hold this many entries, must be at least 1 + maximal sized sized hash

const size_t cWorkThreadStackSize = (1024*1024*2);					// working threads (can be multiple) stack size

const int cMinAceptSeqLen = 50;			// user can specify that reads must be of at least this min length after any end trimming to be further processed
const int cDfltAceptSeqLen = 80;		// default is that reads must be of at least this min length after any end trimming
const int cMaxAceptSeqLen = 500;		// user can specify that reads must be of at least this min length after any end trimming

const int cMaxTrimSeqLen = 1000;        // user can specify that reads are to be trimmed down to this length

const int cMinOverlapbp = (cMinAceptSeqLen+1)/2;			// any putative overlap must be at least this many bp before being further explored 
const int cMinOverlappc = 50;			// user can specify down to this required overlap as a percentage of read length
const int cDfltOverlappc = 70;			// default overlap as a percentage of read length
const int cMaxOverlappc = 95;			// user can specify at most this required overlap as a percentage of read length


typedef enum TAG_eARPMode {
	eAR2Fasta = 0,		// artefact reduce reads to multifasta
	eAR2Packed,			// artefact reduce reads to packed file for subsequent assembly
	eARPacked2fasta		// load packed and output as multifasta files
} etARPMode;

#pragma pack(1)

typedef struct TAG_sThreadIdentDuplicatePars {
	int ThreadIdx;					// index of this thread (1..m_NumThreads)
	void *pThis;					// will be initialised to pt to CKangadna instance
	etPMode PMode;					// processing mode
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// returned result code
	bool bStrand;					// if true then duplicate checking on original read orientation, if false then also check duplicates on antisense
	bool bPEdups;					// sequences are from paired end reads
	bool bDedupeIndependent;		// if paired end preprocessing then treat as if single ended when deuping
	tSeqID StartingSeqID;			// process starting with this sequence identifier (0 to start with 1st)
	tSeqID EndingSeqID;				// process finishing with this sequence identifier (0 to finish with last)
	UINT32 AllocMemProbeSubSeq;		// memory allocated to pProbeSubSeq
	void *pProbeSubSeq;				// to hold subsequence of packed probe sequence as used when exploring overlaps with other sequences
	UINT32 AllocMemPE1SeqWrds;		// memory allocated to pPE1SeqWrds
	void *pPE1SeqWrds;				// used to hold PE1 packed SeqWrds when reverse complementing
	UINT32 AllocMemPE2SeqWrds;		// memory allocated to pPE2SeqWrds
	void *pPE2SeqWrds;				// used to hold PE1 packed SeqWrds when reverse complementing
	UINT32 NumProcessed;			// number processed
	UINT32 NumDuplicates;			// of which this number are duplicates
	UINT32 MaxDuplicates;			// highest number of duplicates encountered by this thread
	UINT32 NumPE1Overlapping;		// number of PE1 sequences which overlapped other sequences
	UINT32 NumPE2Overlapping;		// number of PE2 sequences which overlapped other sequences
	UINT32 NumDupInstances[cMaxDupInstances+1]; // to hold duplicate instances counts

} tsThreadIdentDuplicatePars;


typedef struct TAG_sKMerSeqInst {
	UINT32 NxtSeq;				// offset (multiply by m_KMerSeqInstSize ) into m_pKMerSeqs[] at which next KMer sequence with same hash starts or 0 if no same hashed KMer 
	UINT32 NumInstances;		// number of instances of this KMer
	UINT8 Flags:8;				// to hold sundry flags
	UINT8 PackedSeqs[1];		// to contain the packed m_KMerSeqLen sequence, packed at 4 bases per byte 
	} tsKMerSeqInst;

typedef struct TAG_sThreadKmerDistPars {
	int ThreadIdx;					// index of this thread (1..m_NumThreads)
	void *pThis;					// will be initialised to pt to CReadStats instance
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// returned result code
	int KMerLen;					// KMer length
	INT64 TotNumReads;				// total number of reads processed by this thread
} tsThreadKmerDistPars;


typedef struct TAG_sThreadIdentOverlapPars {
	int ThreadIdx;					// index of this thread (1..m_NumThreads)
	void *pThis;					// will be initialised to pt to CKangadna instance
	etPMode PMode;					// processing mode
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// returned result code
	bool bStrand;					// if true then strand specific processing
	bool bRevCpl;					// if true then all reads have been reverse complemented
	etOvlFlankPhase OvlFlankPhase;	// current flank overlap processing phase
	int ProbeMinOverlap;			// minimum required probe overlap in % of probe length
	int TargMinOverlap;             // minimum required overlap of probe onto target in % of target length
	int MinFlankLen;				// minimum required non-overlap flank
	int OverlapSense;				// 0 sense overlaps sense, 1 antisense overlaps sense, 2 sense overlaps antisense
	tSeqID StartingSeqID;			// process starting with this sequence identifier (0 to start with 1st)
	tSeqID EndingSeqID;				// process finishing with this sequence identifier (0 to finish with last)
	UINT32 AllocMemOverlapSeq;		// memory allocated to pOverlapSeq and pOverlapFlankSeq
	void *pOverlapSeq;				// to hold a copy of packed probe sequence as used (may have been revcpl) when exploring overlaps with other sequences
	void *pOverlapFlankSeq;			// to hold subsequence (flank sequence) of pOverlapSeq
	UINT32 NumProcessed;			// number processed
	UINT32 NumOverlapping;			// number of sequences which overlapped other sequences
	UINT32 NumOverlapped;			// number of sequences determined as being overlapped
	UINT16 FlgOverlapping;			// use this flag as marker for sequences which overlap at least one other sequence 
	UINT16 FlgOverlapped;			// use this flag as marker for sequences which are overlapped by at least one other sequence 
} tsThreadIdentOverlapPars;

#pragma pack()

class CArtefactReduce : public CKangadna
{

	int m_MinAcceptSeqLen;           // filter out input sequences (after any trimming) which are less than this length

	int m_LoadedMeanSeqLen;			// reads loaded post contaminate filtering and flank trimming have this mean sequence length
	int m_LoadedMinSeqLen;			// min length of any read loaded post contaminate filtering and flank trimming
	int m_LoadedMaxSeqLen;			// max length of any read loaded post contaminate filtering and flank trimming 

	int m_KMerSeqLen;				// current KMer sequence length
	int m_KMerSeqInstSize;			// size of a complete tsKMerSeqInst containing m_KMerSeqLen bases
	int m_NumKMerSeqHashes;			// number of hashes currently used in m_pKMerSeqHashes
	UINT32 *m_pKMerSeqHashes;		// array, indexed by hashes over KMerSeq sequences, holding offsets (mult by m_KMerSeqInstSize) into packed m_pKMerSeqs at which KMer sequence instance starts 

	UINT32 m_UsedKMerSeqInsts;		// this many KMer tsKMerSeqInst instances have been used 
	UINT32 m_AllocdKMerSeqInsts;	// allocation was for this many KMer tsKMerSeqInst instances
	size_t m_AllocdKMerSeqInstsMem;	// allocation was for this size
	tsKMerSeqInst *m_pKMerSeqs;		// allocated to hold tsKMerSeqInst instances

	int
		RemoveDuplicates(bool bPEdups,			// can optionally request that duplicates are for both PE1 and PE2 being duplicates
									bool bStrand,			// if true then strand specific duplicates
									char *pszDupDist);		// optionally specify a file to which duplicate distributions are to be written
	int
		RemoveNonOverlaps(int MinOverlap,			// minimum required overlap (in percentage of actual read length)
						int MinFlankLen,            // minimum required non-overlap flank (in bp)
						int NumIterations = 1); 	// because of artefact errors tending to be at end of reads (both 5' and 3') then by default 1 iterations of passes are utilised 

	int									// returns number of KMers of length m_KMerSeqLen accepted from pRead, 0 if none, < 0 if errors
		AddReadKMers(int ReadLen,		// number of bases in read
				etSeqBase *pRead);			// read sequence

	int									 // returns < 0 if errors, 1 if this is the first instance of the KMer sequence, 2..N if multiple instances previously added
		AddReadKMer(etSeqBase *pKMerSeq);	// KMer bases, expected to be of at least length m_KMerSeqLen

public:
	CArtefactReduce(void);
	~CArtefactReduce(void);

	void ARReset(void);
	void ARInit(void);
	
	int
		Process(etARPMode PMode,			// processing mode, currently eAR2Fasta,  eAR2Packed
			char *pszCheckpointFile,		// if file of this name exists and is a checkpoint then resume processing from this checkpoint, otherwise create a checkpoint file 
			etSfxSparsity SfxSparsity,		// suffix sparsity
			int IterativePasses,			// iterative passes of overlap processing
			int MinPhredScore,				// only accept reads for filtering if their minimum Phred score is at least this threshold
			bool bNoDedupe,					// if true then do not remove all duplicate reads as per bDedupeIndependent parameter
			bool bStrand,					// true if read strand specific filtering
			int MaxNs,						// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 0, range 0..5)
			int Trim5,						// trim this number of 5' bases from input sequences (default is 0, range 0..20)
			int Trim3,						// trim this number of 3' bases from input sequences (default is 0, range 0..20)
			int MinSeqLen,		            // filter out input sequences (after any trimming) which are less than this length (default is 80bp, range 50..500)
			int TrimSeqLen,					// trim sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...10000)
			int MinOverlap,					// minimum required overlap (in % of read length) or <= 0 if no overlap processing
			int MinFlankLen,				// non-overlapping flank must be at least this length (defults to 15%% of mean read length if 0, else range 1bp to 25bp, only applies if overlap processing)
			int SampleNth,					// process every Nth reads
			int Zreads,						// maximum number of reads to accept for processing from any file
			bool bDedupeIndependent,		// if paired end preprocessing then treat as if single ended when deuping
			int NumThreads,					// number of worker threads to use
			bool bAffinity,					// thread to core affinity
			int NumPE1InputFiles,			// number of PE1 input files
			char *pszInPE1files[],			// input PE1 5' read files
			int NumPE2InputFiles,			// number of PE2 input files
			char *pszInPE2files[],		    // input PE2 3' read files
			char *pszContaminantFile,		// contaminants fasta file
			char *pszOutFile,				// where to write filtered sequences
			char *pszDupDistFile);			// write duplicate sequence distributions to this file

	int
		IdentifyDuplicates(bool bPEdups,			// request that duplicates are for both PE1 and PE2 being duplicates
									bool bStrand,			// if true then strand specific duplicates
									char *pszDupDist);		// optionally specify a file to which duplicate distributions are to be written



	int
		IdentifyOverlaps(etOvlFlankPhase OvlFlankPhase,		// overlap flank processing phase
						int MinOverlap = 70,					// sequences must flank overlap by at least this percentage of read length
						int MinFlankLen = 1,				// minimum required non-overlap flank (in bp)
						bool bRevCpl = false);				// if true then all sequences have been reverse complemented and sfx index is over these sequences

	int ProcIdentDuplicates(tsThreadIdentDuplicatePars *pPars);		// potentially called by multiple threads!

	int ProcIdentOverlaps(tsThreadIdentOverlapPars *pPars);


};

