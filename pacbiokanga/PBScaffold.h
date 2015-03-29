#pragma once
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

#include "./pacbiocommon.h"

const int cMaxInFileSpecs = 100;			// user can specify upto this many input files
const int cMinSeqLen = 50;					// minimum sequence length accepted for processing
const int cMaxDupEntries = 10;				// report 1st 10 duplicate entry names
const int cMaxAllocBuffChunk = 0x00ffffff;	// buffer for fasta sequences is realloc'd in this sized chunks

typedef enum TAG_ePBPMode {
	ePBPMOverlaps,										// overlap discovery mode, generates output overlap loci detail csv file to use as input in subsequent ePBPMScaffold processing mode
	ePBPMScaffold										// scaffolding mode, uses previously generated overlap loci detail csv file and scaffolds
	} etPBPMode;

#pragma pack(1)

typedef struct TAG_sAntisenseKMerOfs {					// no KMer in antisense strand if both MinOfs and MaxOfs are 0
	UINT32 MinOfs;										// at least one antisense KMer is located between this minimum and
	UINT32 MaxOfs;                                      // and this maximum offset (sub 1) in the antisense sequence
} tsAntisenseKMerOfs;

const int cMaxClusters = 10;	// at most this many clusters will be processed for any probe core hits on to a target
typedef	struct TAG_sCluster {
		UINT32 TargNodeID;		// cluster is with hits to this target
		UINT32 TargSeqLen;		// target is this length
		UINT32 ClustProbeOfs;	// first hit in cluster starts at this probe offset
		UINT32 ClustTargOfs;    // first hit in cluster starts at this target offset
		UINT32 ClustLen;		// cluster length
		UINT32 NumClustHits;	// number of consistency checked aligned hits accepted within this cluster
		UINT32 SumClustHitLens; // sum of all core hit lengths in this cluster  
		double ClustScore;		// score for this cluster
		} tsCluster; 

typedef struct TAG_sCoreHitsClusters {
	UINT32 ProbeID;			// clusters of hits between this probe
	UINT32 ProbeSeqLen;     // probe is of this length bp
	int NumClusters;        // number of clusters in Cluster[]
	tsCluster Clusters[cMaxClusters]; 
} tsCoreHitsClusters;

// seed core hits 
typedef struct TAG_sCoreHit {
	UINT32 ProbeNodeID;				// core hit was from this probe  node 
	UINT32 TargNodeID;				// hit was onto this target  node
	UINT32 ProbeOfs;                // hit was from this probe offset
	UINT32 TargOfs;					// onto this target offset
	UINT32 HitLen;					// hit was of this length
	UINT8 flgRevCpl:1;				// 1 if core sequence was revcpl'd before matching
	UINT8 flgMulti:1;				// 1 if core sequence was target multiloci and this instance to not be further processed
	} tsCoreHit;

// identified overlap between probe and target sequence
typedef struct TAG_sPBOverlaps {
	UINT8 flgAntisense:1;           // probe sequence was reverse complemented
	UINT32 ProbeEntryID;            // probe sequence suffix array identifier
	UINT32 TargEntryID;				// overlap from probe was onto this target suffix array identifier
	UINT32 ProbeStartOfs;           // overlap starts at this probe offset
	UINT32 TargStartOfs;            // overlap starts at this target offset
	UINT32 ProbeOverlapLen;         // probe overlap is of this length
	UINT32 TargOverlapLen;			// target overlap is of this length
} sPBOverlaps;

typedef struct TAG_sPBScaffNode {
	UINT32 NodeID;					// uniquely identifies this node
	UINT32 VertexID;				// assembly graph vertex identifier
	UINT32 EntryID;					// suffix array entry identifier for indexed sequence
	UINT32 SeqLen;					// length in bp of this scaffolding node sequence
	UINT32 flgCurProc:1;			// sequence is currently being processed
	UINT32 flgContained:1;			// sequence is fully contained within another sequence
	UINT32 flgUnderlength:1;        // sequence is under length
} tsPBScaffNode;

typedef struct TAG_sPBCoreHitCnts {
	UINT32 TargNodeID;				// node identifier for hit sequence	
	UINT32 NumSHits;				// number of hits onto target sequence from sense probe
	UINT32 NumAHits;				// number of hits onto target sequence from antisense probe
} sPBCoreHitCnts;

typedef struct TAG_sThreadPBScaffold {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to class instance
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	UINT32 threadID;				// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// processing result

	CMTqsort *pmtqsort;				// muti-threaded qsort

	CSSW *pSW;						// Smith-Waterman class

	bool bSelfHits;					// if true then processing for self hits with objective of detecting retained hairpins

	UINT32 MinOverlapLen;			// the putative overlap would be of at least this length
	bool bRevCpl;					// true if probe sequence is to be revcpl when looking for overlaps
	UINT32 MaxSeedCoreDepth;		// only further extend a seed core if there are no more than this number of matching cores in all targeted sequences
	UINT32 DeltaCoreOfs;			// offset core windows of coreSeqLen along the probe sequence when checking for overlaps 
	UINT32 CoreSeqLen;				// putative overlaps are explored if there are cores of at least this length in any putative overlap
	UINT32 MinNumCores;				// and if the putative overlap contains at least this many cores
	UINT32 MaxAcceptHitsPerSeedCore; // limit accepted hits per seed core to no more this many
	UINT32 MinScaffSeqLen;			// only process scaffold sequences which are at least this length

	UINT32 NumTargCoreHitCnts;		// current number of summary target core hit counts in TargCoreHitCnts
	sPBCoreHitCnts TargCoreHitCnts[cSummaryTargCoreHitCnts]; // top targets by core hit counts
	UINT32 NumCoreHits;				// currently this many core hits in m_pCoreHits
	UINT32 AllocdCoreHits;				// m_pCoreHits currently allocated to hold at most this many core hits
	size_t AllocdCoreHitsSize;		// m_pCoreHits current allocation size
	tsCoreHit *pCoreHits;			// allocated to hold all core hits	

	UINT32 AllocdProbeSeqSize;		// current allocation size for buffered probe sequence in pProbeSeq 	
	etSeqBase *pProbeSeq;			// allocated to hold the current probe sequence

	UINT32 AllocdTargSeqSize;		// current allocation size for buffered target sequence in pTargSeq 	
	etSeqBase *pTargSeq;			// allocated to hold the current probe sequence

	UINT32  AllocdAntisenseKmersSize;  // current allocation size for antisense KMer offset ranges in pAntisenseKmers
	tsAntisenseKMerOfs *pAntisenseKmers;	// used to hold antisense Kmer offset ranges when filtering

	UINT32 AlignErrMem;				// number of times alignments failed because of memory allocation errors
	UINT32 AlignExcessLen;			// number of times alignments failed because length of probe * target was excessive

} tsThreadPBScaffold;


#pragma pack()


class CPBScaffold
{
	etPBPMode m_PMode;						// processing mode

	UINT32 m_DeltaCoreOfs;					// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
	UINT32 m_MaxSeedCoreDepth;				// only further process a seed core if there are no more than this number of matching cores in all targeted sequences

	UINT32 m_MinSeedCoreLen;				// use seed cores of this length when identifying putative overlapping scaffold sequences
	UINT32 m_MinNumSeedCores;				// require at least this many seed cores between overlapping scaffold sequences

	int m_SWMatchScore;						// SW score for matching bases (0..100)
	int m_SWMismatchPenalty;				// SW mismatch penalty (-100..0)
	int m_SWGapOpenPenalty;					// SW gap opening penalty (-100..0)
	int m_SWGapExtnPenalty;					// SW gap extension penalty (-100..0)
	int m_SWProgExtnPenaltyLen;				// SW progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
	
	UINT32 m_NumOverlapProcessed;			// number of PacBio reads processed for overlapping other PacBio reads
	UINT32 m_ProvOverlapping;               // number of PacBio reads overlapping at least one other PacBio read
	UINT32 m_ProvOverlapped;				// number of PacBio reads provisionally overlapped, could be containing, another PacBio read
	UINT32 m_ProvContained;					// number of PacBio readsprovisionally contained within another PacBio read

	UINT32 m_MinScaffSeqLen;				// individual target scaffold sequences must be of at least this length (defaults to 5Kbp)
	UINT32 m_MinScaffOverlap;				// pairs of targeted scaffold sequences must overlap by at least this many bp to be considered for merging into a longer scaffold sequence (defaults to 5Kbp) 
	int m_NumPacBioFiles;					// number of input pacbio file specs
	char m_szPacBioFiles[cMaxInFileSpecs][_MAX_PATH];		// input pacbio files
	int m_NumHiConfFiles;					// number of input hiconfidence file specs
	char m_szHiConfFiles[cMaxInFileSpecs][_MAX_PATH];		// input hiconfidence files		
	char m_szPacBioSfxFile[_MAX_PATH];		// name of input file containing suffix indexed pacbio sequences
	char m_szOutFile[_MAX_PATH];			// where to write merged scaffolded sequences


	char m_szOutScaffFile[_MAX_PATH];		// scaffolding overlap distribution file name
	int m_hScaffDist;						// handle for file into which scaffolding overlap distributions are to be written
	char m_szScaffLineBuff[0x07fff];		// buffering for overlap distributions
	int m_ScaffLineBuffIdx;					// where to write next scaffold overlap

	int m_NumThreads;							// maximum number of worker threads to use

	UINT32 m_NumPBScaffNodes;					// m_pPBScaffNodes currently holds many scaffolding nodes
	UINT32 m_AllocdPBScaffNodes;				// m_pPBScaffNodes allocated to hold this many scaffolding nodes
	tsPBScaffNode *m_pPBScaffNodes;				// allocated to hold scaffolding nodes
	UINT32 *m_pMapEntryID2NodeIDs;				// used to map from suffix array entry identifiers to the corresponding scaffolding node identifier

	CSfxArrayV3 *m_pSfxArray;					// suffix array file (m_szTargFile) is loaded into this

	CAssembGraph *m_pAssembGraph;				// overlapping sequences are assembled into scaffolds using this class

	void Init(void);							// initialise state to that immediately following construction
	void Reset(bool bSync);						// reset state, if bSync true then fsync before closing output file handles
	int LoadTargetSeqs(char *pszTargFile);		// load sequences in this file into in memory suffix array; file expected to contain preindexed sequences 

	int LoadTargetSeqs(int MinSeqLen,int NumTargFiles,char **pszTargFiles);		// parse, and index sequences in this file into in memory suffix array; file expected to contain either fasta or fastq sequences

	int ProcessBioseqFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
				 char *pszFile);				// file containing sequences

	int ProcessFastaFile(int MinSeqLen,			// only accept for indexing sequences of at least this length
				char *pszFile);					// file containing sequences

	int IdentifyScaffoldOverlaps(int MaxSeqLen,		// max length sequence to be overlapped
							int NumOvlpThreads);		// identify all scaffold overlaps using this many threads

	int IdentifyCoreHits(UINT32 ProbeNodeID,	// identify all overlaps of this probe sequence PBScaffNodeID onto target sequences
				tsThreadPBScaffold *pPars);		// thread specific

	int					// returns index 1..N of just added core hit or -1 if errors
		AddCoreHit(UINT32 ProbeNodeID,			// core hit was from this probe scaffold node 
			   bool bRevCpl,					// true if core sequence was revcpl'd before matching
			   UINT32 ProbeOfs,                 // hit started at this probe offset
			   UINT32 TargNodeID,               // probe core matched onto this target scaffold node
			   UINT32 TargOfs,                  // probe core matched starting at this target loci
			   UINT32 HitLen,					// hit was of this length
               tsThreadPBScaffold *pPars);		// thread specific

	int	ClusterSpatial(tsThreadPBScaffold *pThreadPar, 
			   UINT32 ProbeLen,					// probe from which cores were used to align against targets was this length
			   tsCoreHitsClusters *pClusters,	// returned clusters
			   UINT32 MinClusterHits,				// clusters must contain at least this number of consistency checked hits
			   UINT32 MinClusterLen);				// clusters must be at least this length
	
	UINT32										// returned tsPBScaffNode node identifier
		MapEntryID2NodeID(UINT32 EntryID);		// suffix array entry identifier

	CMTqsort m_mtqsort;				// muti-threaded qsort

static int SortLenDescending(const void *arg1, const void *arg2);
static int SortCoreHitsByTargProbeOfs(const void *arg1, const void *arg2);
static int SortCoreHitsByProbeTargOfs(const void *arg1, const void *arg2);


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

public:
	CPBScaffold();
	~CPBScaffold();

	int ThreadIdentScaffOverlaps(tsThreadPBScaffold *pThreadPar);

	int
	Process(etPBPMode PMode,		// processing mode
		UINT32 DeltaCoreOfs,		// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		UINT32 MaxSeedCoreDepth,    // only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		UINT32 MinSeedCoreLen,		// use seed cores of this length when identifying putative overlapping scaffold sequences
		UINT32 MinNumSeedCores,     // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..100)
		int SWMismatchPenalty,		// mismatch penalty (-100..0)
		int SWGapOpenPenalty,		// gap opening penalty (-100..0)
		int SWGapExtnPenalty,		// gap extension penalty (-100..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		UINT32 MinScaffSeqLen,		// individual target scaffold sequences must be of at least this length (defaults to 5Kbp)
		UINT32 MinScaffOverlap,		// pairs of targeted scaffold sequences must overlap by at least this many bp to be considered for merging into a longer scaffold sequence (defaults to 5Kbp) 
		char *pszPacBioOvlps,		// pregenerated PacBio sequence overlap loci details
		char *pszPacBioSfxFile,		// pre-indexed PacBio sequences
		int NumPacBioFiles,			// number of input pacbio file specs
		char *pszPacBioFiles[],		// input pacbio files
		int NumHiConfFiles,			// number of input hiconfidence file specs
		char *pszHiConfFiles[],		// input hiconfidence files		
		char *pszOutFile,			// where to write merged scaffolded sequences
		int NumThreads);			// maximum number of worker threads to use

	int LoadPacBioOvlps(char *pszPacBioOvlps, bool bValidateOnly = false);	// load pregenerated PacBio sequence overlap loci CSV file
};


