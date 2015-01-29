#pragma once

// processing modes
typedef enum TAG_ePMode {
	ePMdefault,					// default processing mode
	ePMMoreSens,				// more sensitive - slower
	ePMUltraSens,				// ultra sensitive - much slower
	ePMLessSens,				// less sensitive - quicker
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// processing result codes returned by LocateOverlaidTarg()
typedef enum etLOTRslt {
	eLOTnone = 0,				// unable to overlap any target sequences
	eLOTprobe,					// probe sequence has been overlapped by some another probe
	eLOTtarg,				    // putative target has been overlapped by some other probe
	eLOThit 				    // accepted overlap onto a target sequence
} tLOTRslt;

const int cMaxSubs = 15;		// max allowed substitution rate per 100
const int cDfltSubs = 7;		// default allowed substitution rate per 100

const int cDfltPMCBaseRate = 7; // default polymorthic base correction rate
const int cMaxPMCBaseRate = 15;  // max allowed polymorthic base correction rate

const int cDfltMinOverlap = 30; // default min overlap of probe sequence onto annother sequence must be at least this number of bases
const int cAbsMinOverlap = 20;	// lower limit on user specified required sequence overlap
const int cAbsMaxOverlap = 200;	// upper limit on user specified required sequence overlap

const int cAbsMinCoreLen = 12;	// core length is 1/3 of minimum overlap but is constrained to be in the range cAbsMinCoreLen to cAbsMaxCoreLen
const int cAbsMaxCoreLen = 50;	// always limit cores to be <= maximum core length
const int cChkIterDepth = 100;  // if more than this many matches of core then check if there are too many matching cores

const int cMaxMergeIters = 10000;    // max allowed merge processing iterations - if more than this number then very likely a processing logic error

const int cAllocMergerSeqLen = 100000;    // allocate to hold merged sequences in these increments - extended as may be required
const int cAllocPMCSeq = 1000;			  // allocate to hold consensus (tsPMCBase) sequences in these increments - extended as may be required

const UINT32 cAllocPathIDs = 100000;		// allocate for path identifiers in this many increments

const UINT32 cAllocCtgLenDist = 10000000;	// allocate to hold this number of contig lengths for calc of N50 etc

const int cReallocNumContigs = 10000000;  // incrementally allocate/realloc (m_pContigs) for this number of sClustRead's
const int cReallocConcatSeqs = cReallocNumContigs; // incrementally allocate for concatenated sequences (m_pConcatSeqs) in this sized memory increments

const int cAllocContigSeqCnts = 10000;  // alloc for contigs in this sized incr
const int cAllocSeqIDsInContig = cAllocContigSeqCnts * 2;	// alloc for this number of reads in contigs

const int cAllocRawSeqLen = 500000000;	// allow for sequences loaded from file of upto this length

const int cAllocLineBuffLen = 1000000;	// sets buffer output length

const UINT32 cMaxSfxBlkEls = 4000000000;	// construct suffix block arrays with 4byte array elements if no more than this, otherwise use 5 byte array elements

const UINT64 cMaxConcatSeqLen = (UINT64)cMaxSfxBlkEls * 50; // limit concatenated read sequences length   

const UINT8 cCSeqBOS = 0xf0;     // marks beginning of concatenated sequences, 1st sequence immediately follows
const UINT8 cCSeqSep = 0x80;	 // separators between concatenated sequences will always have bit 7 set
const UINT8 cCSeqEOS = 0xff;	 // marks end of concatenated sequences



const int cMaxOverlayDist = 101;			// overlaid reads distribution counts

const int cHashAEntries = 0x07fff;		// TargSeqID start loci are hashed into this many entries
//const int cHashAEntries = 0x0fffff;		// TargSeqID start loci are hashed into this many entries
const int cMaxNumAIdentNodes = (cHashAEntries * 128);	// allow at most this many TargSeqIDs to be hash linked per thread

const int cMaxWorkerThreads = 64;			// can handle at most 64 threads
const int cMaxProbesPerBlock = 0x03fff;     // each thread can process at most this many probe reads per block

#pragma pack(1)

typedef struct TAG_sPMCBase {
	etSeqBase PMCBase;					// consensus base which is simply the base with highest count at this loci or if equal counts then the probe base
	UINT32 BaseCnts[5];					// counts for eBaseA..eBaseT and eBaseN
	} tsPMCBase;

typedef struct TAG_sAIdentNode {						// TargIDs are hash linked into tsIdentNodes
	UINT64 TargID;										// used to identify if a subsequence in the target has already been processed
	struct TAG_sAIdentNode *pNxt;						// will be NULL if last linked in current hash chain
} tsAIdentNode;

// concatenated sequences start with tsASeqSep followed by read sequences separated by tsASeqSep with last tsASeqSep set to 0xFFFFFFFFFF
// currently size of this sequence separator is expected to be 6 bytes
typedef struct TAG_sASeqSep {
	UINT8  SeqSep;		// cCSeqBOS if starting 1st sequence, cCSeqSep if starting an intermediate sequence, cCSeqEOS if previous sequence was the final sequence
	UINT32 SeqIDlo;		// 40bit unique sequence identifier with bit 7 set on each byte so these can be distinguished from sequence base bytes
	UINT8  SeqIDhi;		// unique sequence identifier requires 5 bytes
	} tsASeqSep;

// suffixed sequence (initially reads but could be the merge of multiple reads as contigs are assembled)
typedef struct TAG_tsSfxdSeq {
	UINT32 SeqID;		// unique sequence identifier
	UINT64 ConcatSeqOfs; // offset in concatenated sequences (m_pConcatSeqs) at which this sequence starts
	UINT32 ReadLen;		// length of this sequence
	UINT8 Flags;		// holds various combinations of processing state flags for this sequence CAUTION: if changed larger than 1byte then need to synchronise access
} tsSfxdSeq;

// tsSfxdSeq Flags
const UINT8 cFlagPalindrome = 0x01;	// this is a palindromic read
const UINT8 cFlagNA = 0x02;			// this read is not to be processed for assignment to any contig
const UINT8 cFlagMerged = 0x04;	    // this sequence has been merged with another sequence
const UINT8 cFlagDup = 0x08;		// sequence identified as being a duplicate of another

typedef struct TAG_sxRdsSfxDirEl {
	UINT32 BlockID;				// identifies (1..n) this suffix block within this file
    UINT32 SeqID;				// SeqID for initial read sfx'd in this block
	UINT64 NumSuffixEls;	    // number of suffix elements
	UINT32 ElSize;				// each element in this block is this many bytes long - currently will be either 4 or 5
	UINT64 SfxOfs;				// file offset at which suffix array starts
	UINT64 SizeOfSfx;			// total memory (and disk space) required for holding suffix array
	UINT32 NumContigs;			// number of reads sfx'd in this block
	} tsxRdsSfxDirEl;

typedef struct TAG_sRdsSrcFile {
	UINT32  NumContigs;					// number of reads loaded from this source file
	UINT32	SrcFileID;					// uniquely identifies source file (1..N)
	UINT8   SrcFileName[_MAX_PATH];		// reads source file name

} tsRdsSrcFile;

typedef struct TAG_sRdsSfxHdr {
	UINT32 NumSrcFiles;			// number of source files from which reads were loaded
	UINT64 SumReadLens;		    // sum total of all read lengths
	UINT64 ConcatSeqLen;		// length of all concatenated sequences excluding initial cCSeqBOS and final cCSeqEOS
	tsRdsSrcFile RdsSrcFiles[cMaxSrcFiles]; // reads loaded from these source files
	} tsRdsSfxHdr;


typedef struct TAG_sContigSeqID {
	UINT32 ContigID;			 // read sequence is assembled into this contig
	UINT32 RelOfs;				 // read sequence starts at this relative offset in contig
	UINT32 SeqID;				 // read sequence identifier
	} tsContigSeqID;

typedef struct TAG_sAssembThreadPars {
	int ThreadIdx;				// index of this thread (1..m_NumThreads)
	void *pThis;				// will be initialised to pt to CHomozyReduce instance
	etPMode PMode;				// processing mode
	UINT8 *pProbeSeq;			// allocated to hold probe sequence currently being processed by this thread
	UINT32 AllocProbeSeq;		// current allocation for pProbeSeq
	UINT32 ElSize;					// size of each suffix array element - currently either 4 or 5

#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int CurBlockID;					// current suffix block identifier
	int Rslt;						// returned result code
	UINT32 NumContigsProc;			// returned number of contigs processed by this thread instance

	bool bStrand;					// if true then assemble with original read orientation, if false then assemble as strand independent
	int CoreLen;					// core length to use

	int MinCtgLen;					// filter out homozygotic region reduced contigs of less than this length
	int MaxHomozySubs;				// characterise as homozygotic if substitution rate between regions <= this rate per 100bp
	int MinHomozyLen;				// homozygotic regions to be at least this length

	int CurMaxIter;					// max allowed iterations per subsegmented sequence when matching that subsegment
} tsAssembThreadPars;

typedef struct TAG_sProbesBlock {
	int NumContigs;			// number of reads for processing in this block
	int MaxContigs;			// block can hold at most this number of reads
	UINT32 ReadSeqIDs[cMaxProbesPerBlock]; // sequence identifiers for all probe reads to be processed in this block
} tsProbesBlock;
#pragma pack()

class CHomozyReduce
{
	CMTqsort m_MTqsort;		// multithreaded sorting


	UINT32 m_AllocdPathIDs;	// m_pPathIDs currently allocated to hold this many path identifiers
	UINT32 m_NumPathIDs;	// number of path identifiers currently in m_pPathIDs
	UINT32 *m_pPathIDs;		// allocated to hold path identifiers

	int m_NumRawFiles;		// number of raw reads files processed
	int m_NumRdsFiles;		// number of preprocessed (kangar) reads files processed

	int m_hInFile;				// input file handle
	char m_szInFile[_MAX_PATH];	// input file

	int m_hOutFile;				// output file handle
	char m_szOutFile[_MAX_PATH];// output file

	char m_szCtgDescr[80];		// contig descriptor prefix

	int m_NumThreads;			// max number of processing threads to use

	tsRdsSfxHdr m_RdsSfxHdr;	// concatenated reads and suffix array file header

	UINT32 m_TotSeqsParsed;		// total number of sequences parsed
    UINT32 m_TotSeqsUnderLen;	// total number of sequences filtered out because underlength
	UINT32 m_TotSeqsExcessNs;	// total number of sequences filtered out because too many Ns

	UINT32 m_TotSeqs2Assemb;	// original number of sequences to be assembled after length filtering
	UINT64 m_TotSeqs2AssembBases; // original number of nucleotides in sequences to be assembled after flank trimming and length filtering	

	UINT32 m_NumSeqs2Assemb;	// number of sequences in m_pSeqs2Assemb
	size_t m_Seqs2AssembLen;	// current length of sequences in m_pSeqs2Assemb
	size_t m_AllocMemSeqs2Assemb;		// memory currently allocated for m_pSeqs2Assemb 
	UINT8 *m_pSeqs2Assemb;				// holds sequences (always in basespace) concatenated with UINT32 lengths prepended to each sequence which are to be subsequently suffix indexed and assembled

	UINT32 m_NumSfxdSeqs;		// actual number of tsSfxdSeqs in m_pSfxdSeqs
	UINT32 m_AllocdNumSfxdSeqs;	// alloc'd number of tsSfxdSeq
	INT64 m_AllocdMemSfxdSeqs;	// memory alloc'd for holding tsSfxdSeqs
	tsSfxdSeq *m_pSfxdSeqs;		// pts to suffixed sequences or partially assembled contigs

	int m_MeanReadLen;			// mean length of all reads

	UINT64 m_AllocMemConcat;	// allocated memory size for concatenated read sequences
	UINT8 *m_pConcatSeqs;		// to hold all concatenated read/contig sequences
	
	int m_SfxElSize;			// suffix element size - either 4, <= cMaxSfxBlkEls, or 5 if > cMaxSfxBlkEls 
	INT64 m_NumSuffixEls;		// number of elements in suffix array
	INT64 m_AllocMemSfx;	    // allocated memory size for suffix array
	UINT32 *m_pSuffixArray;     // to hold suffix array for concatenated read/contig sequences

	UINT32 m_BuffNumOverlaidContigs;	// buffered number of overlaid reads in m_pFwdOvlAdjacencyArray ready for writing to disk
		
	UINT8 *m_pContigSeq;			    // consensus contig sequence
	size_t m_AllocContigSeq;			// how many UINT8 have been alloc'd to m_pContigSeq 
	size_t m_AllocSeqIDsinContig;		// how many  read identifiers have been alloc'd to m_pSeqIDsInContig
	size_t m_AllocSeqIDsinContigMem;	// memory allocated to m_pSeqIDsInContig
	UINT32 m_NumContigsinContig;			// used number of read identifiers in m_pSeqIDsInContig

	int m_TotNumReduced;				// total number, before any filtering, of contigs length reduced
	int m_TotNumNonReduced;				// total number of contigs not reduced

	int m_NumGenContigs;				// number of contigs written to file and total number of contigs in m_pContigLengths
	size_t m_AllocdCtgLens;				// memory size currently alloc'd to m_pContigLengths 
	int *m_pContigLengths;				// used to hold contig lengths for N50 and other stats generation 
	UINT64 m_TotLenCovered;				// sum of all contig lengths

	char *m_pszLineBuff;				// allocd for buffering of output assembled contigs
	int m_LineBuffLen;					// current number of chars buffered in m_pszLineBuff
	int m_AllocLineBuff;				// m_pszLineBuff allocated to hold at most this number of chars

	size_t m_CurMaxMemWorkSetBytes;     // currently set max working set in bytes, need to convert to pages when setting working set
	UINT32 m_WinPageSize;				// windows memory page size in bytes (0 if process not on Windows)
	size_t m_BaseWinMinMem;				// windows base min working set memory in bytes when process initially started
	size_t m_BaseWinMaxMem;				// windows base max working set memory in bytes when process initially started

	INT64 m_InvalidRefs;				// used to count detected invalid references which currently are simply sloughed, not reported

	int m_MinCtgLen;					// filter out homozygotic region reduced contigs of less than this length
	int m_MaxHomozySubs;				// characterise as homozygotic if substitution rate between regions <= this rate per 100bp
	int m_MinHomozyLen;					// homozygotic regions to be at least this length 

	teBSFrsltCodes GenRdsSfx(void);		// generate suffix array as either 4 or 5byte sfx els
	
	int CmpProbeTarg(etSeqBase *pEl1,etSeqBase *pEl2,int Len); // compare probe to target accounting for any sequence concatenator markers or XFormID etc
	int GetConcatSeqOfs(UINT8 *pSeq); // get offset (0..n) of base ptd at by pSeq within a concatenated sequence
	UINT8 *IterConcatSeqs(UINT8 *pCurSeq); // iterate over sequences in m_pConcatSeqs
	UINT8 *GetpConcatSeqMarker(UINT8 *pSeq); // returns ptr to marker byte terminating concatenated sequence ptd at by pSeq
	int    GetConcatSeqLen(UINT8 *pSeq); // returns length of concatenated sequence ptd at by pSeq
	int GetConcatSeqID(UINT8 *pSeq);	// returns sequence identifier for concatenated sequence ptd at by pSeq
	etSeqBase *GetConcatSeqStart(UINT8 *pSeq);			// returns ptr to first base of concatenated sequence

	teBSFrsltCodes SortContigsByLen(void);  // sort contigs by length

	bool m_bSorted;									// set TRUE after overlaids have been both sorted
	
 	static int SfxSortFunc(const void *arg1, const void *arg2);
	static int Sfx5SortFunc(const void *arg1, const void *arg2);
	static int ContigsSortFunc(const void *arg1, const void *arg2);
	static int SortByContigLen(const void *arg1, const void *arg2);

	int SortOverlaids(bool bResetMarkers=false,bool bForce = false); // bResetMarkers if any overlaids marked for removal are to be removed, bForce if overlays are to be force sorted

	teBSFrsltCodes ChunkedWrite(INT64 WrtOfs,UINT8 *pData,INT64 WrtLen);
	teBSFrsltCodes ChunkedWrite(int hFile,char *pszFile,INT64 WrtOfs,UINT8 *pData,INT64 WrtLen);
	teBSFrsltCodes ChunkedRead(int hFile,char *pszFile,INT64 RdOfs,UINT8 *pData,INT64 RdLen);
	teBSFrsltCodes ChunkedRead(INT64 RdOfs,UINT8 *pData,INT64 RdLen);

	int AddRead(int ReadID,		// read identifier - must be >= 1 and unique
			    int NumDups,	// number of other reads known to have same sequence
				int ReadLen,	// number of bases in *pSeq following
				etSeqBase *pSeq); // read sequence
				
	teBSFrsltCodes Disk2Hdr(tsBSFRdsHdr *pRdsHeader,	// load header into this 
						char *pszRdsFile);			// loading is from this file

	UINT32 ApproxNumContigsAligned(void);
	void ResetThreadedIterContigs(void); // must be called by master thread prior to worker threads calling ThreadedIterContigs()
	void AcquireSerialise(void);
	void ReleaseSerialise(void);
	void AcquireLock(bool bExclusive);				// defaults as read only lock
	void ReleaseLock(bool bExclusive);

	int m_NumSeqsProc;
	int m_NxtReadProcSeqID;
	unsigned long m_ProcessingStartSecs;

	bool	// returns false if no more reads availing for processing by calling thread
			ThreadedIterProbes(tsProbesBlock *pRetBlock);	// iterate and return blocks of read probes to be processed by each thread

	CStopWatch m_StopWatch;

	bool m_bMutexesCreated;						// set true if mutexes and rwlocks created/initialised

#ifdef _WIN32
	HANDLE m_hMtxIterContigs;
	HANDLE m_hMtxMHContigs;
	SRWLOCK m_hRwLock;
	static unsigned __stdcall ThreadedContigsAssemb(void * pThreadPars);
#else
	pthread_mutex_t m_hMtxIterContigs;
	pthread_mutex_t m_hMtxMHContigs;
	pthread_rwlock_t m_hRwLock;
	static void *ThreadedContigsAssemb(void * pThreadPars);
#endif

	int CreateMutexes(void);
	void DeleteMutexes(void);

	static inline UINT64		// unpacks the 5 bytes ptd at by p5Bytes and returns as UINT64 
		Unpack5(UINT8 *p5Bytes)
		{
		// ensures that only 5 bytes are actually accessed, can't just access as a UINT64 and mask retain bottom 40 bits...
		return((UINT64)((UINT64)p5Bytes[4] << 32 | *(UINT32 *)p5Bytes));
		}

	static inline UINT8 *Pack5(UINT64 Value, UINT8 *p5Bytes)
		{
		*(UINT32 *)p5Bytes = (UINT32)Value;
		p5Bytes[4] = (UINT8)(Value >> 32);
		return(p5Bytes + 5);
		}

	teBSFrsltCodes GenSfxdSeqs(void);		// generate suffixed sequences
	teBSFrsltCodes AddSeq(int SeqLen,		// sequence length
						UINT8 *pSeq);		// ptr to read sequence

	teBSFrsltCodes							// load from raw fasta or fastq file
		LoadRawContigs(int MaxNs,			// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
					int Trim5,				// trim this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int MinSeqLen,		    // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)
					int FileID,				// uniquely identifies source file
					char *pszFile);

	void DumpContigSeqs(char *pszDescr);				// dump all contigs to screen
	void DumpContigSeq(int IdentID,char *pszDescr,etSeqBase *pSeq,int SeqLen); // dump individual contig to screen 

	int WriteContigSeq(etSeqBase *pSeq,int SeqLen);


	void ValidateSeq(char *pszDescr,int SeqLen, etSeqBase *pSeq);		// validates that pSeq is valid containing only bases in the range eBaseA..eBaseN

	bool SetMaxMemWorkSetSize(size_t Bytes);

public:
	CHomozyReduce(void);
	~CHomozyReduce(void);

	int Reset(bool bSync = true);

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

	void SetNumThreads(int maxThreads);

	void SetCtgDescr(char *pszCtgDescr);	// set contig descriptor prefix

	teBSFrsltCodes 
		LoadContigs(int MaxNs,				// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
					int Trim5,				// trim this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int MinSeqLen,		    // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)
					 int FileID,			// uniquely identifies source file
					 char *pszFile);		// file containing input contigs

	int GenSuffixArray(void);	// concatenate sequences and generate suffix array over these concatenations

	int ReduceHomozygosity(etPMode PMode,			// processing mode
					char *pszRsltFile,				// write reduced homozygosity contigs into this file
					bool bStrand,					// non-strand specific homozygous region reduction - homozygous regions between any two contigs can be in different orientation
					int MaxHomozySubs,				// characterise as homozygotic if substitution rate between regions <= this rate per 100bp
					int MinHomozyLen,				// homozygotic regions to be at least this length
					int MinHetrozyLen,				// island (marked as homozygotic either flank) hetrozygotic regions must be at least this length otherwise treat as homozygotic
					int MinCtgLen);					// filter out homozygotic region reduced contigs of less than this length
	
	int GenContigSeqs(void);	// process overlapped reads and generate contig sequences
				

		UINT64			// index+1 in pSfxArray of first exactly matching probe or 0 if no match				
		LocateFirstExact(int ElSize,		// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
						etSeqBase *pProbe,  // pts to probe sequence
					  int ProbeLen,					// probe length to exactly match over
					  etSeqBase *pTarg,				// target sequence
					  UINT8 *pSfxArray,				// target sequence suffix array
					  UINT64 TargStart,				// position in pTarg (0..n) corresponding to start of suffix array
					  UINT64 SfxLo,					// low index in pSfxArray
					  UINT64 SfxHi);					// high index in pSfxArray

			UINT64			// index+1 in pSfxArray of last exactly matching probe or 0 if no match					
			LocateLastExact(int ElSize,		// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
					  etSeqBase *pProbe, // pts to probe sequence
					  int ProbeLen,					// probe length to exactly match over
					  etSeqBase *pTarg,				// target sequence
					  UINT8  *pSfxArray,			// target sequence suffix array
					  UINT64 TargStart,				// position in pTarg (0..n) corresponding to start of suffix array
					  UINT64 SfxLo,					// low index in pSfxArray
					  UINT64 SfxHi,					// high index in pSfxArray
					  UINT32 Limit = 0);				// if non-zero then need only iterate towards last exactly matching this for this Limit iterations


			tLOTRslt									// < 0 if errors, eLOTnone if no homozygotic regions identified in this probe, eLOThit if at least one homozygotic region in probe
				MarkHomozygoticRegions(UINT32 ProbeSeqID,// identifies probe sequence
						 etSeqBase *pProbeSeq,			// probe sequence 
						 int ProbeLen,					// probe length 
						 tsAssembThreadPars *pPars);	    // calling thread parameters

};


