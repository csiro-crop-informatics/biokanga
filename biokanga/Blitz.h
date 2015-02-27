// Copyright 2013 CSIRO  ( http://www.csiro.au/ ) 
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License
//   Please contact stuart.stephen@csiro.au for support or 
//   to submit modifications to this source

#pragma once

const int cMinCoreLen = 5;			// minimum allowed core or seed length
const int cDfltCoreLen = 10;		// default core or seed length
const int cMaxCoreLen = 16;			// max allowed core or seed length

const int cDfltMinQueryLenAlignedPct = 25;  // to be accepted a query sequence must align over at least this percentage (1..100) of it's length onto target

const int cMinPathScore = 50;		// user specified min allowed path score before that path will be reported
const int cDfltPathScore = 75;		// default minimum path score before that path will be reported
const int cMaxPathScore = 500;		// user specified max minimum allowed path score before that path will be reported

const int cDfltMaxExtnScoreThres = 12; // default extension score threshold, core overlap extensions with extension score above this threshold are terminated; extension score += 2 if mismatch, extension score -= 1 if match and score > 0
const int cMaxMaxExtnScoreThres = 30;  // maximum accepted extension score threshold, core overlap extensions with extension score above this threshold are terminated; extension score += 2 if mismatch, extension score -= 1 if match and score > 0

const int cMaxOverlapFloat = 8;		// allowing for overlap float of at most this many bases, needed because the cores with extensions are independent and may be overextended

const int cDfltMaxPathsToReport = 10;	// by default report at most this many scoring paths for any query sequence

// using affine gap scoring but limiting the gap extension cost to just the first 100bp
const int cGapOpenCost = 5;			// cost for opening path gap when scoring path
const int cGapExtendCost = 1;		// cost for extending gap per 10bp extension when scoring path
const int cGapExtendCostLimit = 10;		// clamp gap extension cost to be no more than this
const int cGapMaxLength = 100000;    // treat any gaps longer than this length as being not on same path

const int cMinCoreDelta = 1;		// minimum allowed core shift delta in bp
const int cMaxCoreDelta = 50;		// max allowed core shift delta in bp

const int cMaxWorkerThreads = 128;			// limiting max number of threads to this many

const int cDfltSensCoreIters  = 1500;	// default sensitivity core explore depth 
const int cMoreSensCoreIters  = 2000;	// more sensitivity core explore depth
const int cUltraSensCoreIters = 3000;	// ultra sensitivity core explore depth
const int cMinSensCoreIters   = 750;	// min sensitivity core explore depth

const int cMinOccKMerDepth = 100;       // user can override the core exploration search depth from this minimum 
const int cMaxOccKMerDepth = 20000;		// up to this maximum 

const int cMaxQuerySeqIdentLen = 80;		// allow for fasta sequence identifer of upto this length
const int cAllocQuerySeqLen = 0x0200000;	// initially allocate to hold a query sequence of up to this length (2Mbp), will be realloc'd if needed for longer query sequences
const int cMaxQuerySeqLen = (cAllocQuerySeqLen * 8);    // can handle query sequences of up to this maximal length (16Mbp), longer sequences will be truncated to this length and the user warned
const int cMaxReadAheadQuerySeqs = 4000;	// read ahead and enqueue up to at most this many query sequences

const int cNumAllocdAlignNodes = 200000;  // allow each query sequence to have up to this many aligned subsequences

const int cAlignRprtBufferSize = 500000; // buffer for buffering alignment results ready to write to file

#pragma pack(1)

typedef enum TAG_eBLZPMode {
	eBLZPMdefault,						// default processing mode
	eBLZPMplaceholder					// used as a placeholder and sets the range of these enumerations
}etBLZPMode;


typedef enum TAG_eBLZSensitivity {
	eBLZSdefault = 0,					// default processing sensitivity
	eBLZSMoreSens,						// more sensitive - slower
	eBLZSUltraSens,					// ultra sensitive - much slower
	eBLZSLessSens,						// less sensitive - quicker
	eBLZSplaceholder					// used as a placeholder and sets the range of these enumerations
}etBLZSensitivity;


typedef enum TAG_eBLZRsltsFomat {
	eBLZRsltsPSL = 0,	// default results format is PSL
	eBLZRsltsPSLX,		// results format is PSLX
	eBLZRsltsMAF,		// default results format is MAF
	eBLZRsltsBED,		// results as BED
	eBLZRsltsplaceholder   // used as a placeholder and flags the range of these enumerations
}etBLZRsltsFomat;


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

typedef struct TAG_sThreadQuerySeqsPars {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to CBlitz instance

#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	UINT32 NumAllocdAlignNodes;				// number of allocated alignment nodes
	tsQueryAlignNodes *pAllocdAlignNodes;	// allocated to hold aligned subsequences
	tsQueryAlignNodes **ppFirst2Rpts;		// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
	int *pRslt;						// write intermediate result codes to this location
	int Rslt;						// returned result code
} tsThreadQuerySeqsPars;

#pragma pack()


class CBlitz
{
	etBLZPMode m_ProcMode;			// processing mode
	etBLZSensitivity m_Sensitivity;  // alignment sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity
	eALStrand m_AlignStrand;		// align on to watson, crick or both strands of target
	int m_ExtnScoreThres;			// seed extensions terminate if the mismatches score is increased above this threshold; when extending seed matches are scored with -1 if current score > 0, and mismatches scored with 2
	int m_CoreLen;					// use this core length as the exactly matching seed length to be 5' and 3' extended whilst no more than m_MaxSubRate
	int m_CoreDelta;				// offset cores by this many bp
	int m_QueryLenAlignedPct;    	// only report alignment paths if the percentage of total aligned bases to the query sequence length is at least this percentage (1..100)
	int m_MaxIter;					// max allowed iterations (depth) per subsegmented sequence (core) when matching that subsegment
	int m_MinExtdCoreLen;			// cores with flank extensions must be at least this length
	int m_MaxOccKMerDepth;			// maximum depth to explore over-occurring core K-mers
	int m_MinPathScore;				// to be reported alignment paths on any target sequence must score at least this value
	int m_MaxPathsToReport;			// report at most this many alignment paths for any query
	int m_AlignPathID;				// alignment path identifier

	etBLZRsltsFomat m_RsltsFormat;	// output results format
	char *m_pszInputFile;			// name of input file containting query sequences
	char *m_pszSfxFile;				// target as suffix array
	char *m_pszOutFile;				// where to write alignments
	UINT32 m_ReportedPaths;			// total number of aligned paths reported
	UINT32 m_QueriesPaths;			// this many query sequences had at least one reported path
	UINT32 m_NumQueriesProc;		// total number of query sequences processed
	int m_hInFile;					// input file handle

	int m_hOutFile;					// results output file handle
	int m_szLineBuffIdx;			// offset into m_pszLineBuff at which to next write
	char *m_pszLineBuff;			// allocated to hold output line buffering
	CSfxArrayV3 *m_pSfxArray;		// suffix array holds genome of interest
	char m_szTargSpecies[cMaxDatasetSpeciesChrom+1]; // suffix array was generated over this targeted species

	int m_TotSeqIDs;				// total number of query sequences which have been parsed and enqueued
	int m_NumQuerySeqs;				// number of query sequences currently in m_pQuerySeqs
	int m_NxtQuerySeqIdx;			// index into m_pQuerySeqs[] at which to dequeue the next query sequence
	int m_AllocdQuerySeqs;			// number of query sequences allocated
	tsQuerySeq *m_pQuerySeqs;		// allocated to hold array of query sequences
	
	UINT8 m_TermBackgoundThreads; // if non-zero then all background threads are to immediately terminate processing


	int m_NumThreads;				// number of worker threads to use

	bool m_bAllQuerySeqsLoaded;			// set true when all query sequences have been parsed and loaded
	teBSFrsltCodes m_LoadQuerySeqsRslt;	// set with exit code from background query sequences load thread, read after checking if m_bAllQuerySeqsLoaded has been set
	int m_ThreadLoadQuerySeqsRslt;		// returned by query sequence loading thread
#ifdef _WIN32
	unsigned int m_ThreadLoadQuerySeqsID;	// query sequences loading thread identifier
#else
	pthread_t m_ThreadLoadQuerySeqsID;
#endif


	CMTqsort m_mtqsort;				// muti-threaded qsort

	void Init(void);			// initialise state to that immediately following construction
	void Reset(bool bSync);		// reset state, if bSync true then fsync before closing output file handles

	int InitQuerySeqThreads(int NumThreads,			// use this many threads
							int AlignNodes);			// each thread is allocatd this many subsequence alignment nodes


	int InitLoadQuerySeqs(void);		// query sequences are loaded asynchronously to the alignments
	teBSFrsltCodes LoadQuerySeqs(char *pszSeqsFile);	// processing query Sequences from this file


	bool m_bMutexesCreated;			// will be set true if synchronisation mutexes have been created
	int CreateMutexes(void);
	void DeleteMutexes(void);
	void AcquireSerialise(void);
    void ReleaseSerialise(void);
    void AcquireSerialiseMH(void);
	void ReleaseSerialiseMH(void);
	void AcquireLock(bool bExclusive);
	void ReleaseLock(bool bExclusive);


#ifdef _WIN32
	HANDLE m_hMtxIterReads;
	HANDLE m_hMtxMHReads;
	SRWLOCK m_hRwLock;
	HANDLE m_hThreadLoadQuerySeqs;
#else
	pthread_mutex_t m_hMtxIterReads;
	pthread_mutex_t m_hMtxMHReads;
	pthread_rwlock_t m_hRwLock;
#endif

	static int SortQueryAlignNodes(const void *arg1, const void *arg2); // Sort alignment nodes by TargSeqID.QueryID.FlgStrand.QueryStartOfs.TargStartOfs ascending
	static int SortHighScoreDescend(const void *arg1, const void *arg2); // Sort alignment nodes which are the first in path by score descending

public:
	CBlitz();
	~CBlitz();
	int
	Process(etBLZPMode PMode,			// processing mode
			etBLZSensitivity Sensitivity,	// sensitivity 0 - standard, 1 - high, 2 - very high, 3 - low sensitivity
			eALStrand AlignStrand,			// align on to watson, crick or both strands of target
			int  CoreLen,					// use this core length as the exactly matching seed length to be 5' and 3' extended whilst no more than m_MaxSubRate
			int  CoreDelta,					// offset cores by this many bp
			int MaxExtnScoreThres,			// terminate overlap extension if curent extension score more than this; if mismatch then extension score += 2, if match and score > 0 then score -= 1 
			int MaxOccKMerDepth,			// maximum depth to explore over-occurring core K-mers
			int  MinPathScore,				// only report alignment paths on any target sequence if the path score is >= this minimum score
			int QueryLenAlignedPct,				// only report alignment paths if the percentage of total aligned bases to the query sequence length is at least this percentage (1..100)
			int  MaxPathsToReport,			// report at most this many alignment paths for any query
			etBLZRsltsFomat RsltsFormat,	// output results format
			char *pszInputFile,				// name of input file containting query sequences
			char *pszSfxFile,				// target as suffix array
			char *pszOutFile,				// where to write alignments
			int NumThreads);				// number of worker threads to use

	int ProcLoadQuerySeqsFile(tsLoadQuerySeqsThreadPars *pPars);
	int ProcAlignQuerySeqs(tsThreadQuerySeqsPars *pPars);

	teBSFrsltCodes LoadRawQuerySeqs(char *pszSeqsFile);	// process query Sequences from this file
	int										// returned enqueued query identifier
		EnqueueQuerySeq(char *pszQueryIdent,    // query identifier
			int QuerySeqLen,				// query sequence length
			UINT8 *pQuerySeq);				// query sequence
	UINT8 *										// returned dequeued sequence, caller is responsible for deleting memory allocated to hold the returned sequence (delete pRetSeq;)
		DequeueQuerySeq(int WaitSecs,		// if no sequences available to be dequeued then wait at most this many seconds for a sequence to become available
			int MaxLenQueryIdent,			// maximum length query identifier
			int *pSeqID,					// returned sequence identifier
			char *pszQueryIdent,			// where to return query identifier
			int *pQuerySeqLen);				// where to return query sequence length

	int AlignSeqs(void);

	int
		BlocksAlignStats(UINT32 *pMatches,	// returned number of bases that match that aren't repeats
					UINT32 *pmisMatches,	// returned number of bases that don't match
					UINT32 *prepMatches,	// returned number of bases that match but are part of repeats
					UINT32 *pnCount,		// returned number of 'N' bases
				char  Strand,				// query sequence strand, '+' or '-')
				UINT8 *pQuerySeq,			// the query sequence
				UINT32 qSize,				// Query sequence size
				UINT32 TargSeqID,			// CSfxArray sequence identifier
				UINT32 tSize,				// Target sequence size 
				UINT32 TargPathStartOfs,	// at this starting offset
				UINT32 TargPathEndOfs,		// ending at this offset (inclusive)
				UINT32 NumPathNodes,		// number of alignment nodes in alignment path
				int SortedPathIdx,
				tsQueryAlignNodes *pAlignNodes,		// alignment nodes
				tsQueryAlignNodes **ppFirst2Rpts); // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt

	int Report( int MinPathScore,			// only report paths having at least this minimum score
				int  MaxPathsToReport,		// report at most this many alignment paths for any query
				char *pszQuerySeqIdent,		// query sequence name 
				UINT32 QueryLen,			// query length
				UINT8 *pQuerySeq,			// the query sequence
				UINT32 NumNodes,			// number of alignment nodes
				tsQueryAlignNodes *pAlignNodes, // alignment nodes
				tsQueryAlignNodes **ppFirst2Rpts);	// allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt

	int		// reporting alignment as PSL format
			ReportAsPSL(UINT32 Matches,			// Number of bases that match that aren't repeats
					UINT32 misMatches,			// Number of bases that don't match
					UINT32 repMatches,			// Number of bases that match but are part of repeats
					UINT32 nCount,				// Number of 'N' bases
					UINT32	qNumInsert,			// Number of inserts in query
					UINT32 qBaseInsert,			// Number of bases inserted in query
					UINT32 tNumInsert,			// Number of inserts in target
					UINT32 tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					UINT32 qSize,				// Query sequence size
					UINT32 qStart,				// Alignment start position in query
					UINT32 qEnd,				// Alignment end position in query
					char *pszTargName,			// aligning to this target
					UINT32 tSize,				// Target sequence size 
					UINT32 TargPathStartOfs,	// at this starting offset
					UINT32 TargPathEndOfs,		// ending at this offset (inclusive)
					UINT32 NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts);  // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt

	int		// reporting alignment as PSLX format
			ReportAsPSLX(UINT32 Matches,			// Number of bases that match that aren't repeats
					UINT32 misMatches,			// Number of bases that don't match
					UINT32 repMatches,			// Number of bases that match but are part of repeats
					UINT32 nCount,				// Number of 'N' bases
					UINT32	qNumInsert,			// Number of inserts in query
					UINT32 qBaseInsert,			// Number of bases inserted in query
					UINT32 tNumInsert,			// Number of inserts in target
					UINT32 tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					UINT32 qSize,				// Query sequence size
					UINT8 *pQuerySeq,			// the query sequence
					UINT32 qStart,				// Alignment start position in query
					UINT32 qEnd,				// Alignment end position in query
					UINT32 TargSeqID,			// CSfxArray sequence identifier
					char *pszTargName,			// aligning to this target
					UINT32 tSize,				// Target sequence size 
					UINT32 TargPathStartOfs,	// at this starting offset
					UINT32 TargPathEndOfs,		// ending at this offset (inclusive)
					UINT32 NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts);  // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt


	int		// reporting alignment as MAF format
			ReportAsMAF(int PathScore,			// score for this path
					UINT32 Matches,				// Number of bases that match that aren't repeats
					UINT32 misMatches,			// Number of bases that don't match
					UINT32 repMatches,			// Number of bases that match but are part of repeats
					UINT32 nCount,				// Number of 'N' bases
					UINT32	qNumInsert,			// Number of inserts in query
					UINT32 qBaseInsert,			// Number of bases inserted in query
					UINT32 tNumInsert,			// Number of inserts in target
					UINT32 tBaseInsert,			// Number of bases inserted in target
					char  Strand,				// query sequence strand, '+' or '-'
					char *pszQuerySeqIdent,     // this query sequence
					UINT32 qSize,				// Query sequence size
					UINT8 *pQuerySeq,			// the query sequence
					UINT32 qStart,				// Alignment start position in query
					UINT32 qEnd,				// Alignment end position in query
					UINT32 TargSeqID,			// CSfxArray sequence identifier
					char *pszTargName,			// aligning to this target
					UINT32 tSize,				// Target sequence size 
					UINT32 TargPathStartOfs,	// at this starting offset
					UINT32 TargPathEndOfs,		// ending at this offset (inclusive)
					UINT32 NumPathNodes,		// number of alignment nodes in alignment path
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts);  // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt
	
	int			// reporting alignment as BED format
			ReportAsBED(char *pszQuerySeqIdent,     // this query sequence
					char  Strand,				// query sequence strand, '+' or '-'
					UINT32 AlignScore,			// alignment has this score
					char *pszTargName,			// aligning to this target
					UINT32 NumPathNodes,		// number of alignment nodes in alignment path
					UINT32 TargPathStartOfs,	// at this starting offset
					UINT32 TargPathEndOfs,		// ending at this offset (inclusive)
					int SortedPathIdx,
					tsQueryAlignNodes *pAlignNodes,		// alignment nodes
					tsQueryAlignNodes **ppFirst2Rpts); // allocated to hold ptrs to alignment nodes which are marked as being FlgFirst2tRpt

	int	IdentifyHighScorePaths(UINT32 QueryLen,			// query length
			UINT32 TargSeqLen,					// targeted sequence length
			bool bStrand,						// reporting for paths on this strand - false if sense, true if antisense
			UINT32 NumNodes,					// reporting on this number of nodes starting from StartNodeIdx
			UINT32 StartNodeIdx,				// report for nodes starting at this node index (1..NumNodes) which is expected to be the first alignment node of a new target sequence
			tsQueryAlignNodes *pAlignNodes,		// alignment nodes
			int MinPathScore,					// only report those series having at least this score
			int  MaxPathsToReport);				// report at most this many alignment paths for any query

	int											// returned best score for paths starting at pAlignNodes[ExploreNodeIdx]
		HighScoreSW(UINT32 QueryLen,			// query length
			UINT32 TargSeqLen,					// targeted sequence length
 			bool bStrand,						// scoring for series on this strand - false if sense, true if antisense
			UINT32 ExploreNodeIdx,				// node to be explored for maximally scored path
			UINT32 NumNodes,					// total number of alignment nodes 
			tsQueryAlignNodes *pAlignNodes);	// alignment nodes
			
};

