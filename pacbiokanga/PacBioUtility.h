#pragma once

const int cMaxQuerySeqIdentLen = 80;		// allow for fasta sequence identifer of upto this length
const int cAllocQuerySeqLen = 0x030000;	// initially allocate to hold a query sequence of up to this length, will be realloc'd if needed for longer query sequences
const int cMaxQuerySeqLen = (cAllocQuerySeqLen * 32);    // can handle query sequences of up to this maximal length, longer sequences will be truncated to this length and the user warned
const int cMaxReadAheadQuerySeqs = 4000;	// read ahead and enqueue up to at most this many query sequences


const etSeqBase cSmartBellAdaptorSeq[] = // "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT" PacBio SMRTBell hairpin sequence
			{ eBaseA,eBaseT,eBaseC,eBaseT,eBaseC,eBaseT,eBaseC,eBaseT,eBaseC,eBaseT,eBaseT,eBaseT,eBaseT
			,eBaseC,eBaseC,eBaseT,eBaseC,eBaseC,eBaseT,eBaseC,eBaseC,eBaseT,eBaseC,eBaseC,eBaseG,eBaseT
			,eBaseT,eBaseG,eBaseT,eBaseT,eBaseG,eBaseT,eBaseT,eBaseG,eBaseT,eBaseT,eBaseG,eBaseA,eBaseG
			,eBaseA,eBaseG,eBaseA,eBaseG,eBaseA,eBaseT}; // PacBio SMRTBell hairpin sequence
const int cSmartBellAdaptorSeqLen = sizeof(cSmartBellAdaptorSeq); // number of bases in PacBio SMRTBell hairpin sequence

const int cHomopolymerTruncLen = 4;           // default is to truncate homopolymers of longer than this length to this length, min 1
const int cMinSmartBellTetramers = 7;		  // default is to require at least this many SMRTBell hairpin sequence tetramers to be discovered and in expected order				


#pragma pack(1)
typedef struct TAG_sSMRTBellHit {
	int NumTetramers;		// number of SMRTBell sequence tetramers identified
	int LocOfs;				// offset + 1 within the pacbio read at which the initial tetramer was identified
} tsSMRTBellHit;

typedef struct TAG_sQuerySeq {
    int SeqID;						// monotonically increasing unique sequence identifier
	char szQueryIdent[cMaxQuerySeqIdentLen+1];	// fasta identifier
	int QuerySeqLen;				// query sequence length
	UINT8 *pQuerySeq;				// allocated to hold sequence 
} tsQuerySeq;

typedef struct TAG_sLoadQuerySeqsThreadPars {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to CBlitz instance

#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int *pRslt;						// write intermediate result codes to this location
	int Rslt;						// returned result code
} tsLoadQuerySeqsThreadPars;


#pragma pack()

class CPacBioUtility
{
	int m_NumInputFiles;			// number of input files
	char **m_ppszInputFiles;		// names (wildcards allowed) of input files containing reads to be filtered
	char m_szAsyncReadsFile[_MAX_PATH];	// reads curently being loaded from this file
	bool m_bAsyncLoading;				// true if async loading still being processed and not yrt completed
	bool m_bMutexesCreated;				// set true if synchronisation mutexes etc created

	UINT8 m_TermBackgoundThreads;		// if non-zero then all background threads are to immediately terminate processing

	int m_TotSeqIDs;				// total number of query sequences which have been parsed and enqueued
	int m_NumQuerySeqs;				// number of query sequences currently in m_pQuerySeqs
	int m_NxtQuerySeqIdx;			// index into m_pQuerySeqs[] at which to dequeue the next query sequence
	int m_AllocdQuerySeqs;			// number of query sequences allocated
	tsQuerySeq *m_pQuerySeqs;		// allocated to hold array of query sequences
	
	bool m_bAllQuerySeqsLoaded;			// set true when all query sequences have been parsed and loaded
	teBSFrsltCodes m_LoadQuerySeqsRslt;	// set with exit code from background query sequences load thread, read after checking if m_bAllQuerySeqsLoaded has been set
	int m_ThreadLoadQuerySeqsRslt;		// returned by query sequence loading thread
#ifdef _WIN32
	unsigned int m_ThreadLoadQuerySeqsID;	// query sequences loading thread identifier
#else
	pthread_t m_ThreadLoadQuerySeqsID;
#endif


	int InitLoadQuerySeqs(void);		// query sequences are loaded asynchronously

	int										// returned enqueued query identifier
		EnqueueQuerySeq(char *pszQueryIdent,    // query identifier
			int QuerySeqLen,				// query sequence length
			UINT8 *pQuerySeq);				// query sequence

	int CreateMutexes(void);
	void DeleteMutexes(void);
	void AcquireSerialise(void);
    void ReleaseSerialise(void);
	void AcquireLock(bool bExclusive);
	void ReleaseLock(bool bExclusive);

#ifdef _WIN32
	HANDLE m_hMtxIterReads;
	SRWLOCK m_hRwLock;
	HANDLE m_hThreadLoadQuerySeqs;
#else
	pthread_mutex_t m_hMtxIterReads;
	pthread_rwlock_t m_hRwLock;
#endif

public:
	CPacBioUtility();
	~CPacBioUtility();

	void Reset(void);
	
	int ProcLoadReadsFile(tsLoadQuerySeqsThreadPars *pPars);

	teBSFrsltCodes StartAsyncLoadSeqs(int NumInputFiles,			// number of input files
									char **ppszInputFiles,		// names (wildcards allowed) of input files containing reads to be filtered
									int MaxReadahead = cMaxReadAheadQuerySeqs);	 // read ahead for at most this number of sequences

	UINT8 *										// returned dequeued sequence, caller is responsible for deleting memory allocated to hold the returned sequence (delete pRetSeq;)
		DequeueQuerySeq(int WaitSecs,		// if no sequences available to be dequeued then wait at most this many seconds for a sequence to become available
			int MaxLenQueryIdent,			// maximum length query identifier
			int *pSeqID,					// returned sequence identifier
			char *pszQueryIdent,			// where to return query identifier
			int *pQuerySeqLen);				// where to return query sequence length

	int										// identified SMRTBell at this offset + 1, 0 if non-detected
		DetectSMRTBell(int StartOfs,		// search from this offset
						int *pNumTetramers,	// returned number of tetramers from SmartBell sequence detected and which were in expected order 
						int InSeqLen,		// number of bases in sequence to search for SmartBell
						etSeqBase *pInSeq,		// sequence to search for SmartBell
						int MinTetramers = cMinSmartBellTetramers);	 // only report SmartBell if at least this many tetramers in expected order detected

	int													// number of putatve SMRTBell hairpins identified
		IdentifySMRTBells(int MaxSMRTBells,	// identify at most this many SMRTBells to return in pSMRTBells
						int SeqLen,			// length of sequence to search for SMRTBell hairpins
						etSeqBase *pSeq,	// identify all putative SMRTBell hairpins in this sequence
						tsSMRTBellHit *pSMRTBells, // returned identified SMRTBells
						int MinTetramers = cMinSmartBellTetramers);	// only report SmartBell if at least this many tetramers in expected order detected

	int		// length after hompolymer reduction
		ReduceHomopolymers(int InSeqLen,	// number of bases in sequence
						   etSeqBase *pInSeq,		// sequence which may contain homopolymers
						   etSeqBase *pOutSeq,		// homopolymer reduced sequence copied into this sequence buffer
						   int TruncLen);			// homopolymers longer than this length (min 1) to be truncated at this length

};

