#pragma once

const int cMinBootstraps = 1;			// min allowed query bootstrap iterations
const int cDfltBootstraps = 1000;		// default bootstrap iterations
const int cMaxBootstraps = 10000;		// max allowed query bootstrap iterations

const int cDfltIterTargs = 1000;		// default allowed target bootstrap iterations
const int cMaxIterTargs = 10000;		// max allowed target bootstrap iterations

const UINT32 cDfltBootstrappingAttempts = 1000;	// default is to allow this many bootstrapping attempts at generating a sample set before returning an error 
const UINT32 cDfltSamplingAttempts = 10000;		// default is to allow this many sampling attempts before starting a new bootstrapping attempt

const int cMaxWorkerThreads = 128;			// limiting max number of threads to this many

const int cMaxAssembSeqLen = 0x5fffffff;	// accepting individual assembly sequences if no longer than this length

const int cMaxQuerySeqLen = 0x0ffff;		// accepting individual query sequences if no longer than this length
const int cMaxTargSeqLen = 0x0fffff;		// accepting individual target sequences if no longer than this length

const int cNumSeqDescrs = 1000;				// initially allocate for this number of block descriptors
const int cAvgSeqDescrLen = 40;				// assuming descriptors will average this length, not incl trailing '\0' terminator
const int cTruncSeqDescrLen = 80;			// descriptors will be truncated to be no longer than this length, not incl trailing '\0' terminator
const int cDfltSeqBlocks = cNumSeqDescrs;	// initially allocating for this many sequence blocks in any sequence source, will realloc in increments of this many
const int cDfltSeqAlloc = 0x0ffffff;		// initially allocating for concatenated sequences totaling this size from any sequence source, will be realloc in increments of at least this size

typedef enum {
	ePMBSAdefault = 0, // default processing mode 
	ePMSAreportseqs,   // report on actual sequences used when bootstrapping
	ePMBSAPlaceholder
} ePMBSAlign;

#pragma pack(1)

typedef enum {
	ePMBSSQuerySeqs = 0, // sequence loaded from query sequence
	ePMBSSTargSeqs,			// sequence loaded from target sequence
	ePMBSSQueryAssemb,  // sequence loaded from query assembly
	ePMBSSTargAssemb,  // sequence loaded from target assembly
	ePMBSSrcPlaceholder
} ePMBSSeqSrc;

typedef struct TAG_sSeqDescr {
	ePMBSSeqSrc SeqSrc;     // descriptor source - 0: query seqs, 1: target sequences, 2: query assembly, 3: target assembly
	UINT32 SeqBlockID;		// descriptor is for this sequence block
	UINT8 DescrLen;			// descriptor has been truncated to this length including terminator
	UINT8 Descr[1];			// \0 terminated descriptor
	} tsSeqDescr;

typedef struct TAG_sSeqBlock {
	ePMBSSeqSrc SeqSrc;     // sequence source - 0: query seqs, 1: target sequences, 2: query assembly, 3: target assembly
	UINT32 SeqBlockID;		// uniquely identifies this sequence block
	UINT64 DescrOfs;		// offset in m_pszDescr of descriptor for this sequence block
	UINT32 SeqLen;			// sequence is this length
	UINT32 NumQueryHits;	// number of query hits
	INT64 SmplSeqOfs;		// offset in  pSeqs at which sequence starts 
	INT64 PopSeqOfs;		// sampled sequence starts at this offset in population pSeqs, -1 if unknown
	} tsSeqBlock; 

typedef struct TAG_sSeqAllocs {
	ePMBSSeqSrc SeqSrc;					// sequence source - 0: query seqs, 1: target sequences, 2: query assembly, 3: target assembly
	UINT32 NumSeqs;						// number of sequences loaded into m_pSeqs
	size_t UsedSeqsSize;				// m_pSeqs used size
	size_t AllocdSeqsSize;				// m_pSeqs allocation size
	etSeqBase *pSeqs;					// concatenated sequences from which to bootstrap sample
	UINT32 UsedSeqBlocks;				// number of used blocks
	size_t UsedSeqBlocksSize;			// number of sequence blocks currently allocated
	size_t SeqBlocksAllocSize;			// m_pSeqBlocks allocation size
	tsSeqBlock *pSeqBlocks;				// sequence blocks containing sequence lengths, descriptors, and where each sequence starts in m_pSeqs
	UINT32 NumSeqDescrs;				// number of sequence descriptors in m_pSeqDescrs
	size_t UsedSeqDescrsSize;			// sequence descriptors used size
	size_t AllocdSeqDescrsSize;		// m_pSeqs allocation size	
	tsSeqDescr *pSeqDescrs;			// allocated to hold all sequence descriptors as parsed from multifasta file
	} tsSeqAllocs;

typedef struct TAG_sQueryHit {
	UINT8 flgHit:1;					// set if hit discovered
	UINT8 flgAntisense:1;			// set if reported hit was with antisense query sequence
	UINT8 MMCnt;					// hit was with this many mismatches
	UINT32 TargIdx;					// query hit was to this target
	UINT32 TargOfs;					// query hit starts at this target offset
	} tsQueryHit;

typedef struct TAG_sWorkerInstance {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to class instance
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	UINT32 threadID;				// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	UINT32 AlignReqID;				// alignment request last processed by this thread; alignments only performed when AlignReqID != m_AlignReqID
	int Rslt;						// processing result
	UINT32 StartQuerySeqIdx;		// thread instance to process alignments starting with this query sequence
	UINT32 EndQuerySeqIdx;			// through to this query sequence inclusive
} tsWorkerInstance;

#pragma pack()

class CAlignsBootstrap
{
	ePMBSAlign m_PMode;		// bootstrap processing mode
	int m_RandSeed;			// random generator seed
	bool m_bSenseOnly;		// true if to align sense only, false if to align both sense and antisense
	int m_MaxSubs;			// allowing at most this many subs as percentage of query length before accepting alignment
	bool m_bWORreplacement;	// sample without replacement, default is to sample with replacement
	bool m_bNoOverlaps;		// sample without overlap, default is to sample allowing overlaps

	int m_NumBootstraps;	// number of bootstrap iterations, excludes initial original query sequences aligned onto initial target sequences	

	int m_CurBootstrap;		// current bootstrap iteration
	bool m_bUseTargBS;		// true if to align against target bootstraps, false if aligning against original target sequences 
	bool m_bUseQueryBS;		// true if to align with query bootstraps, false if align with original query sequences

	int m_hCSVQRslts;		// CSV query results file handle
	int m_hCSVTRslts;		// CSV target results file handle

	int m_NumThreads;		// max number of threads in worker thread pool

	tsWorkerInstance m_WorkerInstances[cMaxWorkerThreads];	// to hold all worker instance thread parameters

#ifdef WIN32
	alignas(4) volatile UINT32  m_NumWorkerInsts;				// number of worker instance threads actually started
	alignas(4) volatile UINT32 m_AlignReqID;					// bootstrap sampleset identifier, incremented if new sampleset available to be aligned
	alignas(4) volatile UINT32 m_CompletedWorkerInsts;			// number of worker instance threads completed current bootstrap alignments
	alignas(4) volatile UINT32 m_TermAllThreads;                // will be set to 1 if all worker threads are to terminate
#else
	__attribute__((aligned(4))) volatile UINT32  m_NumWorkerInsts;				// number of worker instance threads actually started
	__attribute__((aligned(4))) volatile UINT32 m_AlignReqID;			// bootstrap sampleset identifier, incremented if new sampleset available to be aligned
	__attribute__((aligned(4))) volatile UINT32 m_CompletedWorkerInsts;			// number of worker instance threads completed current bootstrap alignments
	__attribute__((aligned(4))) volatile UINT32 m_TermAllThreads;                  // will be set to 1 if all worker threads are to terminate
#endif

	CRandomMersenne *m_pRandomMersenne;		// used for generating random numbers larger than RAND_MAX (32767)

	tsSeqAllocs m_Seqs[ePMBSSrcPlaceholder]; // sequences loaded from each source - 0: query seqs, 1: target sequences, 2: query assembly, 3: target assembly
	tsQueryHit *m_pQueryHits;			// allocated to hold all query hits to targets

	int m_AllocRsltsBuff;				// summary results buffers allocated to hold at most this many chars
	int m_CurQRsltsOfs;					// offset into m_pszQRsltsBuff at which to write next query hits counts
	int m_CurTRsltsOfs;					// offset into m_pszTRsltsBuff at which to write next target hits counts
	char *m_pszQRsltsBuff;				// allocated for query summary results buffering
	char *m_pszTRsltsBuff;				// allocated for target summary results buffering

	int
		AlignBootstrap(int NumRepeats);  // number of times current set of counts are to be reported

	int
		LoadFastaSeqs(ePMBSSeqSrc SeqSrc,   // descriptor source - 0: query seqs, 1: target sequences, 2: query assembly, 3: target assembly
				int MinSeqLen,			// only accepting sequences which are at least this length
				char *pszFastaFile);	// fasta file to load from

	int
		ReportBootstrapSeqs(bool bTargSeqs,		// true if target sequences to be reported 
					int Iteration,		// which bootstrap iteration (1..n)
					char *pszSeqsFile);  // write bootstraps into this file, will have bootstrap iteration specific suffix appended

	int
	AddSeq(ePMBSSeqSrc SeqSrc,     // descriptor source - 0: query seqs, 1: target sequences, 2: query assembly, 3: target assembly
			 char *pszDescr,		// sequence descriptor
			 UINT32 SeqLen,			// sequence is this length
			 UINT8 *pSeqBuff);		// sequence

	INT64		// returned random number will be at most 60bits (2^60)
		GenRand60(INT64 Limit);	// generate random number between 0 and Limit inclusive where Limit is <= 2^60

	int GenBootstrap(UINT32 BootstrapAttempts = cDfltBootstrappingAttempts, // allow at most this many attempts at bootstrapping a set of samples before returning error
					UINT32 SampleAttempts = cDfltSamplingAttempts, // allow at most this many attempts at randomly locating a sample before restarting the bootstrap
					bool bTargs = false,   // false: generate bootstrap sampling from query assembly sequences, true: bootstrap sampling from target assembling sequences
					bool bWithoutReplacement = false, // false: sampling with replacement, true: sampling without replacement (currently not implemented)
					bool bNonOverlapping = false);	   // false: samples may be overlapping, true: samples must be non-overlapping (currently not implemented)

	int	// index, -1 if no matches, at which Query matched onto target with at most MaxSubs
		AlignQueriesToTargs(bool bSenseOnly,			// true if to align sense only, default is to align both sense and antisense
						UINT32 StartQuerySeqIdx,		// alignments starting with this query sequence
						UINT32 EndQuerySeqIdx,			// through to this query sequence inclusive
						int MaxSubs);				// accepting at most this percentage of bases of query length to be mismatches

	// initialise and start pool of worker threads
	int		StartWorkerThreads(UINT32 NumThreads,		// there are this many threads in pool
					UINT32 NumQuerySeqs);	// which will be processing a total of this many query sequences

	bool	// true if any worker threads in pool to start alignments, false if no worker threads 
			StartAlignments(void);			// signal worker pool of threads that there is a new bootstrap set to be aligned

	bool	// true if pool of worker threads completed current bootstrap set within WaitSecs, false if at least thread still processing	
			WaitAlignments(int WaitSecs=60);	// allow at most this many seconds for pool of worker threads to complete aligning current bootstrap set

	int		TerminateWorkerThreads(int WaitSecs = 120);				// alow at most this many seconds before force terminating threads



public:
	CAlignsBootstrap();
	~CAlignsBootstrap();

	void Reset(void);
	int Init(int RandSeed);				// if > 0 then random generator seed, otherwise time() used as the seed

	int Process(ePMBSAlign PMode,			// bootstrap processing mode
				int RandSeed,				// if > 0 then random generator seed, , otherwise time() used as the seed
				bool bSenseOnly,			// true if to align sense only, false if to align both sense and antisense
				int MaxSubs,				// allowing at most this many subs as percentage of query length before accepting alignment
				int NumBootstraps,			// number of bootstrap iterations, excludes initial original query sequences aligned onto initial target sequences
				bool bWORreplacement,		// sample without replacement, default is to sample with replacement
				bool bNoOverlaps,			// sample without overlap, default is to sample allowing overlaps
				char *pszQuerySeqsFile,		// fasta file containing initial query sequences from which to derive query length distributions
				char *pszTargSeqsFile,		// fasta file containing initial target sequences from which to derive target length distributions
				char *pszQueryAssembFile,	// file containing fasta assembly to be bootstrap sampled for query sequences with same length distributions as sequences in pszQuerySeqsFile 
				char *pszTargAssembFile,	// file containing fasta assembly to be bootstrap sampled for target sequences with same length distributions as sequences in pszTargSeqsFile
				char *pszQRsltsFile,		// summary number of query hits onto at least one target bootstrap results to this file 
				char *pszTRsltsFile,		// summary number of targets hit by at least one query bootstrap results to this file
				int NumThreads);			// number of worker threads to use 

	int ProcWorkerThread(tsWorkerInstance *pThreadPar);	// worker thread parameters

};

