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
#include "./PacBioUtility.h"

const int cMinFiltSeedCoreLen = 8;					// filtering seed cores can be specified down to this bp
const int cDfltFiltSeedCoreLen = 10;				// default filtering seed core length
const int cMaxFiltSeedCoreLen = 12;					// filtering seed cores can be specified up to this bp

const int cDfltTrim5 = 250;							// default is to after any filtering then to trim 5' end by this many bp
const int cDfltTrim3 = 250;							// default is to after any filtering then to trim 5' end by this many bp
const int cDfltMinReadLen = 5000;					// default is to only report reads of at least this minimum read length after any end trimming

const int cDfltMaxPacBioSeqLen = 0x1ffff;			// default is to allow for PacBio reads of <= 128Kbp
const int cAllocOutBuffSize = 0x0fffff;				// m_pOutBuff allocate output buffer size 


const int cDfltMinClustLen = 1000;                  // default min cluster length
const int cDfltMinCoresCluster = 3;                 // require at least this many cores per cluster
const int cDfltCoreSeparation = 700;                // default max core separation between cores in the same cluster


typedef enum TAG_ePBPMode {
	ePBPMFilter										// filter and optionally trim PacBio reads
	} etPBPMode;

#pragma pack(1)

typedef struct TAG_sAntisenseKMerOfs {					// no KMer in antisense strand if both MinOfs and MaxOfs are 0
	UINT32 MinOfs;										// at least one antisense KMer is located between this minimum and
	UINT32 MaxOfs;                                      // and this maximum offset (1+) in the antisense sequence
} tsAntisenseKMerOfs;


typedef struct TAG_sFiltCoreHitsCluster {
	char Status;			// character code representing cluster status - 'U' probe sequence is under length, 'A' sequence accepted, 'R' sequence rejected as score is above threshold 
	UINT32 ProbeID;			// clusters of hits between this probe and it's antisense
	UINT32 ProbeSeqLen;     // probe is of this length bp
	UINT32 ClustProbeOfs;	// first hit in cluster starts at this probe sense offset
	UINT32 ClustTargOfs;    // first hit in cluster starts at this probe antisense offset
	UINT32 ClustLen;		// cluster length
	UINT32 NumClustHits;	// number of consistency checked aligned hits accepted within this cluster
	UINT32 SumClustHitLens; // sum of all core hit lengths in this cluster  
	double ClustScore;		// score for this cluster
	} tsFiltCoreHitsCluster;

// seed core hits 
typedef struct TAG_sFiltCoreHit {
	UINT32 ProbeNodeID;				// core hit was from this probe  node 
	UINT32 ProbeOfs;                // hit was from this probe offset
	UINT32 TargOfs;					// onto this antisense target offset
	UINT32 HitLen;					// hit was of this length
	} tsFiltCoreHit;


typedef struct TAG_sPBCoreHitCnts {
	UINT32 TargNodeID;				// node identifier for hit sequence	
	UINT32 NumSHits;				// number of hits onto target sequence from sense probe
	UINT32 NumAHits;				// number of hits onto target sequence from antisense probe
} sPBCoreHitCnts;

typedef struct TAG_sThreadPBFilter {
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
	UINT32 NumTargCoreHitCnts;		// current number of summary target core hit counts in TargCoreHitCnts
	sPBCoreHitCnts TargCoreHitCnts[cSummaryTargCoreHitCnts]; // top targets by core hit counts
	UINT32 NumCoreHits;				// currently this many core hits in m_pCoreHits
	UINT32 AllocdCoreHits;				// m_pCoreHits currently allocated to hold at most this many core hits
	size_t AllocdCoreHitsSize;		// m_pCoreHits current allocation size
	tsFiltCoreHit *pCoreHits;			// allocated to hold all core hits	

	UINT32 AllocdTargSeqSize;		// current allocation size for buffered target sequence in pTargSeq 	
	etSeqBase *pTargSeq;			// allocated to hold the current probe sequence

	UINT32  AllocdAntisenseKmersSize;  // current allocation size for antisense KMer offset ranges in pAntisenseKmers
	tsAntisenseKMerOfs *pAntisenseKmers;	// used to hold antisense Kmer offset ranges when filtering

	UINT32 AlignErrMem;				// number of times alignments failed because of memory allocation errors
	UINT32 AlignExcessLen;			// number of times alignments failed because length of probe * target was excessive

} tsThreadPBFilter;


#pragma pack()


class CPBFilter
{
	etPBPMode m_PMode;						// processing mode
	double m_MaxScore;						// only accept reads with scores <= this score
	int m_MinSeedCoreLen;					// use seed cores of this length when identifying putative antisense subsequences indicative of SMRTBell read past onto antisense
	int m_MaxCoreSeparation;				// cores part of the same cluster must not be separated by more than this many bp
	int m_MinClustLen;						// any putative sense overlap with antisense subsequence cluster must be of at least this length
	int m_MinCoresCluster;                 // require at least this many cores per cluster
	int m_Trim5;							// 5' trim accepted reads by this many bp
	int m_Trim3;							// 3' trim accepted reads by this many bp
	int m_MinReadLen;						// read sequences must be at least this length after any end timming
	char m_szInputFile[_MAX_PATH];			// name of input file containing reads to be filtered
	
	char m_szOutFile[_MAX_PATH];			// name of file in which to write filter accepted and trimmed sequences
	int m_hOutFile;							// handle for file into which to write filter accepted and trimmed sequences
	char *m_pOutBuff;						// used to buffer sequences passing filtering ready to be written to file
	int m_AllocOutBuffSize;					// m_pOutBuff allocated to hold at most this many chars
	int m_OutBuffIdx;						// current number of chars in m_pOutbuff

	char m_szOutFilt[_MAX_PATH];			// write filtered out read details to this file
	int m_hOutFilt;							// handle for file into which filtered out read details are to be written
	char m_szFiltLineBuff[0x07fff];			// buffering for filtered out reads details
	int m_FiltLineBuffIdx;					// where to write next scaffold overlap

	int m_NumThreads;							// maximum number of worker threads to use
	int m_TotProcessed;						// total reads processed
	int	m_TotAccepted;						// after filtering accepted this number of reads
	int	m_TotRejected;                      // rejected this number of reads because of antisense cluster of core hits detected meeting score threshold etc
	int m_TotUnderLen;						// this number reads not accepted because they were underlength

	void Init(void);							// initialise state to that immediately following construction
	void Reset(void);							// reset state
	int LoadTargetSeqs(char *pszTargFile);		// loading sequences form this file 

	int ProcessFiltering(int MaxSeqLen,			// max length sequence expected
							int NumOvlpThreads);	// filtering using at most this many threads

	int IdentifyCoreHits(UINT32 SeqID,			// sequence identifier 
				int SeqLen,						// sequence is this long bp
				etSeqBase *pProbeSeq,				// sequence to process for core hits from sense on to the antisense
				tsThreadPBFilter *pPars);		// thread specific

	int												// number of putatve SmartBell hairpins identified
		IdentifySmartBells(UINT32 ProbeNodeID,	// identify all overlaps of this probe sequence ProbeNodeID onto target sequences
					tsThreadPBFilter *pPars);	// thread specific


	int					// returns index 1..N of just added core hit or -1 if errors
		AddCoreHit(UINT32 ProbeNodeID,			// core hit was from this probe read 
			   UINT32 ProbeOfs,                 // hit started at this probe offset
			   UINT32 TargOfs,                  // probe core matched starting at this antisense offset
			   UINT32 HitLen,					// hit was of this length
               tsThreadPBFilter *pPars);		// thread specific

	int	ClusterSpatial(tsThreadPBFilter *pThreadPar, 
			   UINT32 ProbeLen,					// sequence from which antisense cores were identified was this length
			   tsFiltCoreHitsCluster *pCluster,		// returned cluster
			   UINT32 MinClusterHits,			// cluster must contain at least this number of consistency checked hits
			   UINT32 MinClusterLen);			// cluster must be at least this length
	

	// Sort core hits by ProbeOfs.TargOfs ascending, note no need to sort on ProbeNodeID as always same
	CMTqsort m_mtqsort;				// muti-threaded qsort for the sorting
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

	CPacBioUtility m_PacBioUtility;

public:
	CPBFilter();
	~CPBFilter();

	int PBFilterReads(tsThreadPBFilter *pThreadPar);

	int			// number of targets identified having core hits in concordance with probe cores
		IdentConcordTargs(tsThreadPBFilter *pThreadPar);

	int
	Process(etPBPMode PMode,	// processing mode
		double MaxScore,			// only accept reads with scores <= this score
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative antisense subsequences indicative of SMRTBell read past onto antisense
		int MaxCoreSeparation,      // max core separation between cores in the same cluster
		int MinClustLen,			// any putative sense overlap with antisense subsequence cluster must be of at least this length
		int Trim5,					// 5' trim accepted reads by this many bp
		int Trim3,					// 3' trim accepted reads by this many bp
		int MinReadLen,				// read sequences must be at least this length after any end timming
		char *pszInputFile,			// name of input file containing reads to be filtered
		char *pszOutFile,			// name of file in which to write filter accepted and trimmed sequences
		int NumThreads);			// maximum number of worker threads to use

	int
	GenTargIndex(UINT32 MinScaffoldLen,	// individual indexed scaffold sequences must be at least this length
				 char *pszTargFile);	// load sequences from this file

};



