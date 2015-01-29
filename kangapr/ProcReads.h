#pragma once


const int cAllocRawSeqLen = 50000;	// allow for P1/P2 read sequences loaded from file of up to this length
const UINT32 cMaxSfxBlkEls = 4000000000;	// construct suffix block arrays with 4byte array elements if no more than this, otherwise use 5 byte array elements

const INT64 cMaxConcatSeqLen = (INT64)cMaxSfxBlkEls * 50; // limit concatenated read sequences length to total no more than this length   
const size_t cMinConcatSeqLen = (size_t)100000;	// always allocate at least this for holding concatenated read sequences

const int cReallocNumReads = 2000000;  // incrementally allocate/realloc for this number of reads
const int cReallocConcatSeqs = (cReallocNumReads * 100); // incrementally allocate for concatenated sequences in this sized memory increments

const UINT8 cCSeqBOS = 0xf0;     // marks beginning of concatenated sequences, 1st sequence immediately follows
const UINT8 cCSeqSep = 0x80;	 // separators between concatenated sequences will always have bit 7 set
const UINT8 cCSeqEOS = 0xff;	 // marks end of concatenated sequences

#pragma pack(1)
typedef struct TAG_sSeqStarts {
	UINT32 PairID;					// uniquely identifies this pair
	UINT16 P1SeqLen;				// P1 sequence length
	UINT16 P2SeqLen;				// P2 sequence length
	size_t P1SeqOfs;				// P1 sequence starts at this offset + 1 in m_pP1Seqs2Assemb
	size_t P2SeqOfs;				// P2 sequence starts at this offset + 1 in m_pP2Seqs2Assemb
} tsSeqStarts;

typedef struct TAG_sRdsSrcFile {
	UINT32  NumReads;					// number of reads loaded from this source file
	UINT32	SrcFileID;					// uniquely identifies source file (1..N)
	UINT8   SrcFileName[_MAX_PATH];		// reads source file name
} tsRdsSrcFile;

typedef struct TAG_sRdsSfxHdr {
	UINT32 NumSrcFiles;			// number of source files from which reads were loaded
	UINT64 SumReadLens;		    // sum total of all read lengths
	UINT64 ConcatSeqLen;		// length of all concatenated sequences excluding initial cCSeqBOS and final cCSeqEOS
	tsRdsSrcFile RdsSrcFiles[cMaxSrcFiles]; // reads were loaded from these source files
	} tsRdsSfxHdr;


#pragma pack()

class CProcReads
{
	bool m_bIsPairedEndProc;	// if true then processing is for paired end P1 plus P2, if false then single ended P1 processing
	int m_NumThreads;			// max number of processing threads to use
	CMTqsort m_MTqsort;			// multithreaded sorting

	UINT32 m_TotSeqsParsed;		// total number of sequences parsed
    UINT32 m_TotSeqsUnderLen;	// total number of sequences filtered out because underlength
	UINT32 m_TotSeqsExcessNs;	// total number of sequences filtered out because too many Ns

	UINT32 m_TotP1Seqs2Assemb;	// original number of sequences to be assembled after length filtering
	UINT32 m_NumP1Seqs2Assemb;	// number of sequences in m_pP1Seqs2Assemb
	size_t m_P1Seqs2AssembLen;	// current length of sequences in m_pP1Seqs2Assemb
	UINT32 m_TotP2Seqs2Assemb;	// original number of sequences to be assembled after length filtering
	UINT32 m_NumP2Seqs2Assemb;	// number of sequences in m_pP2Seqs2Assemb
	size_t m_P2Seqs2AssembLen;	// current length of sequences in m_pP2Seqs2Assemb
	size_t m_AllocMemP1Seqs2Assemb;	// memory currently allocated for P1 m_pP1Seqs2Assemb 
	size_t m_AllocMemP2Seqs2Assemb;	// memory currently allocated for P2 m_pP2Seqs2Assemb 
	UINT8 *m_pP1Seqs2Assemb;		// holds P1 sequences (always in basespace) to assemble into contigs, concatenated with UINT32 lengths prepended to each sequence
	UINT8 *m_pP2Seqs2Assemb;		// holds P2 sequences (always in basespace) concatenated with UINT32 lengths prepended to each sequence

	UINT32 m_AllocSeqStarts;			// number of allocated read sequence starts
	size_t	m_AllocMemSeqStarts;		// memory allocated for holding the sequence starts
	UINT32 m_NumSeqStarts;				// number of sequence starts
	tsSeqStarts *m_pSeqStarts;			// pts to array of sequence starts


	int m_MeanReadLen;			// mean length of all reads

	int m_SfxElSize;			// suffix element size - either 4, <= cMaxSfxBlkEls, or 5 if > cMaxSfxBlkEls 
	INT64 m_NumSuffixEls;		// number of elements in suffix array
	INT64 m_AllocMemSfx;	    // allocated memory size for suffix array
	UINT32 *m_pSuffixArray;     // to hold suffix array for concatenated read/contig sequences

	size_t m_CurMaxMemWorkSetBytes;     // currently set max working set in bytes, need to convert to pages when setting working set
	UINT32 m_WinPageSize;				// windows memory page size in bytes (0 if process not on Windows)
	size_t m_BaseWinMinMem;				// windows base min working set memory in bytes when process initially started
	size_t m_BaseWinMaxMem;				// windows base max working set memory in bytes when process initially started

	tsRdsSfxHdr m_P1RdsSfxHdr;			// concatenated P1 reads header
	tsRdsSfxHdr m_P2RdsSfxHdr;			// concatenated P2 reads header

	CFasta m_P1InFasta;							// used to parse P1 reads
	CFasta m_P2InFasta;							// used to parse P2 reads

	char m_szP1OutFile[_MAX_PATH];				// accepted P1 reads output file name
	char m_szP2OutFile[_MAX_PATH];				// accepted P2 reads output file name
	char m_szP1RjtFile[_MAX_PATH];				// rejected P1 reads output file name
	char m_szP2RjtFile[_MAX_PATH];				// rejected P2 reads output file name

	char *PrfxFname(char *pszPfx,			// prefix to use 
		char *pszFile);						// file path + name to be prefixed

	int m_NumRawFiles;						// number of raw reads files processed

	CStopWatch m_StopWatch;

	// serialisations and locks
	bool m_bMutexesCreated;						// set true if mutexes and rwlocks created/initialised
	void AcquireSerialise(void);
	void ReleaseSerialise(void);
	void AcquireLock(bool bExclusive);				// defaults as read only lock
	void ReleaseLock(bool bExclusive);
#ifdef _WIN32
	HANDLE m_hMtxIterReads;
	HANDLE m_hMtxMHReads;
	SRWLOCK m_hRwLock;
	static unsigned __stdcall ThreadedReadsAssemb(void * pThreadPars);
#else
	pthread_mutex_t m_hMtxIterReads;
	pthread_mutex_t m_hMtxMHReads;
	pthread_rwlock_t m_hRwLock;
	static void *ThreadedReadsAssemb(void * pThreadPars);
#endif

	int CreateMutexes(void);
	void DeleteMutexes(void);

	bool SetMaxMemWorkSetSize(size_t Bytes);

	teBSFrsltCodes AllocReadsMemory(size_t P1ReqAllocSize,size_t P2ReqAllocSize); // allocate memory as may be required


public:
	CProcReads(void);
	~CProcReads(void);
		// Transforms 31bit identifier into a 40bit identifier such that the MSB of each byte is set to 1
	// Objective is to allow the identifier to be easily distinguished from ordinary sequence data which has the MSB reset
	static inline UINT64 IDtoXForm(UINT32 ID)
		{
		UINT64 XFormID;
		if(ID > 0x07fffffff)
			return( 0x0ffffffff);
		XFormID = 0x8080808080 | ((ID & 0x07f) | (ID & 0x3f80) << 1 | (ID & 0x01fc000) << 2 | (ID & 0x0fe00000) << 3 | (ID & 0x07f0000000) << 4);  
		return(XFormID);
		}

	// XFormToID
	// Transforms a transformed 40bit identifier (generated by IDtoXForm()) back into it's original 31bit identifier
	static inline UINT32 XFormToID(UINT64 XFormID)
		{
		UINT32 ID;
		ID = (UINT32)((XFormID & 0x07f) | (XFormID & 0x07f00) >> 1 | (XFormID & 0x07f0000) >> 2 | (XFormID & 0x07f000000) >> 3 | (XFormID & 0x07f00000000) >> 4);  
		return(ID);
		}

	// SetNumThreads
	// Set allowed number of threads
	void SetNumThreads(int maxThreads);
	void Init(void);	// class instance initiation
	void Reset(bool bFlush = false);	// if true then flush output to file before closing file handles, reset back to state immediately following instantiation of this class instance

	void SetPairedEndProc(bool bIsPairedEndProc);

	// LoadRawReads
	// Load reads from fasta or fastq formated raw reads file
	// basic assumption is that the paired ends are in the same order in their respective paired end files
	// so when reads are loaded from file P1 and discarded for whatever reason, then the corresponding read can be discarded from P2
	
	teBSFrsltCodes
		LoadRawReads(int MaxNs,			// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
					int Trim5,				// trim this fixed number of 5' bases from input sequences (default is 0, range 0..1000)
					int Trim3,				// trim this fixed number of 3' bases from input sequences (default is 0, range 0..1000)
					int MinSeqLen,		    // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)
					int FileID,				// uniquely identifies source file P1 and P2
					char *pszP1File,		// process from this P1 file 
					char *pszP2File,		// optionally process from this P2 file
					char *pszAcceptPfx,		// output accepted reads into a file with same name as corresponding pszP1File/pszP2File but prefixed with this  
					char *pszRejectPfx);		// if not null then output rejected reads into a file with same name as corresponding pszP1File/pszP2File but prefixed with this

	// AddSeq
	// Add sequences which are to be subsequently suffix indexed and assembled
	// basespace only sequences expected
	teBSFrsltCodes								// returned SeqID or if <= 0 then error result code
			AddSeq( int P1SeqLen,			// P1 sequence length
						UINT8 *pP1Seq,		// ptr to P1 sequence
						int P2SeqLen,		// P2 sequence length
						UINT8 *pP2Seq);		// ptr to P1 sequence
};

