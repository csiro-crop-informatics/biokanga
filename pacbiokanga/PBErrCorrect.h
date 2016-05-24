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

const int cMaxInFileSpecs = 200;			// user can specify upto this many input files
const int cMinSeqLen = 500;					// minimum sequence length accepted for processing
const int cMaxDupEntries = 10;				// report 1st 10 duplicate entry names
const int cMaxAllocBuffChunk = 0x00ffffff;	// buffer for fasta sequences is realloc'd in this sized chunks

const int cFlgLCSeq = 0x01;					// sequence is a low confidence read (PacBio) containing many sequence errors
const int cFlgHCSeq = 0x02;					// sequence a high confidence sequence with few sequence errors

const int cDfltMinHCSeqLen = 1000;			// default min high confidence sequence length
const int cDfltMinHCSeqOverlap = 500;		// default min high confidence sequence overlap length onto PacBio sequence
const int cDfltHCRelWeighting = 3;              // default high confidence sequence overlap relative score weighting (standard overlaps have a weighting of 1)

const UINT8 cLCWeightingFactor = (0x01);   // if a low confidence sequence used to generate a consensus then use this relative weighting (1x) on bases from this sequence, note bit 7 reset as this is a low confidence sequence
const UINT8 cHCDfltWeightingFactor = (0x83);   // if a high confidence sequence used to generate a consensus then use this relative weighting (4x) on bases from this sequence, note bit 7 set as this is a high confidence sequence

const int cMinMaxArtefactDev = 1;			// user can specify down to this minimum or 0 to disable
const int cDfltMaxArtefactDev = 50;			// default percentage deviation from the mean allowed when classing overlaps as being artefactual
const int cDfltMaxConsolidateArtefactDev = 20;	// default percentage deviation from the mean allowed when classing overlaps as being artefactual if consolidating
const int cDfltScaffMaxArtefactDev = 20;		// but when scaffolding with error corrected reads there should a much smaller percentage deviation from the mean
const int cMaxMaxArtefactDev = 70;			// user can specify up to this maximum 

const int cReqConsensusCoverage = 20;   // targeting this mean probe sequence coverage to have confidence in consensus error correction

const int cMaxPacBioErrCorLen = cMaxSWQuerySeqLen;			// allowing for error corrected read sequences of up to this length
const int cMaxPacBioMAFLen    = cMaxSWMAFBuffSize;			// allowing for multialignment format buffering of up to this length

const int cDfltRMIReqDataSize = cMaxSWReqPayloadSize;		// RMI: each worker thread default allocates to process up to this much request data
const int cDfltRMIReqParamSize = cMaxSWParamLen;			// RMI: each worker thread default allocates to process up to this much parameterisation data
const int cDfltRMIRespDataSize = cMaxSWRespPayloadSize;		// RMI: each worker thread default allocates to return up to this much response data
const int cDfltRMIBufferSize =   cMaxSWMAFBuffSize;			// each worker thread default allocates to hold at most this sized MAlignCols2fasta/MAlignCols2MFA alignments

const UINT32 cRMI_SecsTimeout = 180;				// allowing for most RMI SW requests to take at most this many seconds to complete (request plus response)
const UINT32 cRMI_AlignSecsTimeout = 600;			// allowing for a RMI SW alignment request to take at most this many seconds to complete (request plus response)
const UINT32 cRMIThreadsPerCore = 8;				// current guesstimate is that 1 server core can support this many RMI SW threads ( 1 core per Non-RMI SW thread)
                                                    // predicated on assuming that the qualifying of read pairs for SW requires around 20% of per core time, the other 80% is spent on SW

typedef enum TAG_ePBPMode {								// processing mode
	ePBPMErrCorrect,									// error correct
	ePBPMConsensus,										// generate consensus from previously generated multiple alignments
	ePBPMOverlapDetail,									// generate overlap detail from previously generated consensus sequences
	ePBMConsolidate										// consolidate error corrected sequences into representative sequences (usually used to generate representative transcripts) 
	} etPBPMode;


#pragma pack(1)

// seed core hits 
typedef struct TAG_sPBECoreHit {
	UINT32 ProbeNodeID;				// core hit was from this probe  node 
	UINT32 TargNodeID;				// hit was onto this target  node
	UINT32 ProbeOfs;                // hit was from this probe offset
	UINT32 TargOfs;					// onto this target offset
	UINT32 HitLen;					// hit was of this length
	UINT32 WinHits;					// number of core hits onto target relative to this core which are within a window of probelen
	UINT8 flgRevCpl:1;				// 1 if core sequence was revcpl'd before matching
	UINT8 flgMulti:1;				// 1 if core sequence was target multiloci and this instance to not be further processed
	UINT8 flgClustered:1;			// 1 if core hit identified as part of a cluster of hits
	} tsPBECoreHit;

typedef struct TAG_sPBEScaffNode {
	UINT32 NodeID;					// uniquely identifies this node
	UINT32 VertexID;				// assembly graph vertex identifier
	UINT32 EntryID;					// suffix array entry identifier for indexed sequence
	UINT32 SeqLen;					// length in bp of this scaffolding node sequence
	UINT8 flgCurProc:1;				// sequence is currently being processed
	UINT8 flgCpltdProc:1;			// sequence processing completed
	UINT8 flgContained:1;			// sequence is fully contained within at least one other sequence
	UINT8 flgContains:1;			// sequence fully contains at least one other sequence
	UINT8 flgUnderlength:1;         // sequence is under length
	UINT8 flgHCseq:1;				// loaded as a high confidence (non-PacBio) sequence
} tsPBEScaffNode;

typedef struct TAG_sECChkPt {
	UINT32 NodeID;					// uniquely identifies this node
	UINT32 EntryID;					// suffix array entry identifier for indexed sequence
	UINT32 SeqLen;					// length in bp of this scaffolding node sequence
	INT64  ECFileOfs;				// error corrected reads file offset post write of error corrected sequences for this node
	UINT8 flgCpltdProc:1;			// sequence processing completed
	UINT8 flgContained:1;			// sequence is fully contained within at least one other sequence
	UINT8 flgContains:1;			// sequence fully contains at least one other sequence
	UINT8 flgUnderlength:1;         // sequence is under length
	UINT8 flgHCseq:1;				// loaded as a high confidence (non-PacBio) sequence	
	} tsECChkPt;


typedef struct TAG_sPBECoreHitCnts {
	UINT32 TargNodeID;				// node identifier for hit target sequence	
	UINT32	STargStartOfs;			// lowest target offset for any sense hit from probe
	UINT32	STargEndOfs;			// highest target offset for any sense hit from probe
	UINT32	ATargStartOfs;			// lowest target offset for any antisense hit from probe
	UINT32	ATargEndOfs;			// highest target offset for any antisense hit from probe
	UINT32	SProbeStartOfs;			// lowest probe offset for any sense hit onto target
	UINT32	SProbeEndOfs;			// highest probe offset for any sense hit onto target
	UINT32	AProbeStartOfs;			// lowest probe offset for any antisense hit onto target
	UINT32	AProbeEndOfs;			// highest probe offset for any antisense hit onto target
	UINT32 NumSHits;				// number of hits onto target sequence from sense probe
	UINT32 NumAHits;				// number of hits onto target sequence from antisense probe
	UINT8 flgProbeHCseq:1;          // set if probe was loaded as a high confidence (non-PacBio) sequence
	UINT8 flgTargHCseq:1;           // set if target was loaded as a high confidence (non-PacBio) sequence
} sPBECoreHitCnts;

typedef struct TAG_sThreadPBErrCorrect {
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
	UINT32 MinPBSeqLen;				// only process PacBio sequences which are at least this length

	UINT32 NumTargCoreHitCnts;		// current number of summary target core hit counts in TargCoreHitCnts
	sPBECoreHitCnts TargCoreHitCnts[cSummaryTargCoreHitCnts+1]; // top targets by core hit counts
	UINT32 NumCoreHits;				// currently this many core hits in m_pCoreHits
	UINT32 AllocdCoreHits;				// m_pCoreHits currently allocated to hold at most this many core hits
	size_t AllocdCoreHitsSize;		// m_pCoreHits current allocation size
	tsPBECoreHit *pCoreHits;			// allocated to hold all core hits	

	UINT32 AllocdProbeSeqSize;		// current allocation size for buffered probe sequence in pProbeSeq 	
	etSeqBase *pProbeSeq;			// allocated to hold the current probe sequence

	UINT32 AllocdTargSeqSize;		// current allocation size for buffered target sequence in pTargSeq 	
	etSeqBase *pTargSeq;			// allocated to hold the current target sequence


	UINT32 AlignErrMem;				// number of times alignments failed because of memory allocation errors
	UINT32 AlignExcessLen;			// number of times alignments failed because length of probe * target was excessive

	UINT32 ErrCorBuffIdx;			// index into m_szErrCorLineBuff at which to next copy a corrected sequence
	UINT32 AllocdErrCorLineBuff;	// allocation size for m_pszErrCorLineBuff
	char *pszErrCorLineBuff;		// allocated buffering for error corrected and scored sequences
	UINT32 MultiAlignBuffIdx;		// index into m_pszMultiAlignLineBuff at which to write next multialignment
	UINT32 AllocdMultiAlignLineBuff;// allocation size for m_pszMultiAlignLineBuff
	char *pszMultiAlignLineBuff;	// allocated buffering for multialignments

// RMI support
	bool bRMI;						// set true if using RMI SW methods
	teBKSPType ServiceType;			// requesting this service type - currently only eBKSPTSmithWaterman supported
	CBKSRequester *pRequester;		// requests for service are made through this class

	tsSSWCell RMIHighScoreCell;		// highest scoring cell in SW alignment
	UINT32 RMIReqDataSize;			// pRMIReqData allocation size
	UINT8 *pRMIReqData;				// allocated to hold request data
	UINT32 RMIParamDataSize;		// pRMIParamData allocation size
	UINT8 *pRMIParamData;			// allocated to hold parameter data
	UINT32 RMIRespDataSize;			// pRMIRespData allocation size
	UINT8 *pRMIRespData;			// allocated to hold response data
	UINT32 RMIBufferSize;			// RMIBuffer allocation size
	UINT8 *pRMIBuffer;				// allocated to hold buffered data needing a life time extending past the RMI_ function

} tsThreadPBErrCorrect;


#pragma pack()


class CPBErrCorrect
{
	etPBPMode m_PMode;						// processing mode

	bool m_bRMI;							// set true if using RMI SW methods

#ifdef WIN32
	FILETIME m_PrevIdleTime;				// used to load balance between RMI and non-RMI threads
	FILETIME m_PrevKernelTime;				// by checking overall process times every 60 secs and adjusting number of non-RMI threads to minimise idle time yet ensure there is at least 1 sec idle per 60 secs elapsed 
	FILETIME m_PrevUserTime;
#endif

	char m_szRMIHostName[cMaxHostNameLen];	// listening on this host name if RMI SW methods utilised
	char m_szRMIServiceName[cMaxServiceNameLen];	// listening on this service/port name for RMI SW service providers
	CBKSRequester *m_pRequester;			// if using RMI SW methods then will be initialised to pt to an instance of CBKSRequester

	UINT32 m_DeltaCoreOfs;					// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
	UINT32 m_MaxSeedCoreDepth;				// only further process a seed core if there are no more than this number of matching cores in all targeted sequences

	UINT32 m_MinSeedCoreLen;				// use seed cores of this length when identifying putative overlapping scaffold sequences
	UINT32 m_MinNumSeedCores;				// require at least this many seed cores between overlapping scaffold sequences

	int m_SWMatchScore;						// SW score for matching bases (0..100)
	int m_SWMismatchPenalty;				// SW mismatch penalty (-100..0)
	int m_SWGapOpenPenalty;					// SW gap opening penalty (-100..0)
	int m_SWGapExtnPenalty;					// SW gap extension penalty (-100..0)
	int m_SWProgExtnPenaltyLen;				// only apply gap extension penalty if gap at least this length (1..63) - use if aligning PacBio
	
	UINT32 m_NumOverlapProcessed;			// number of PacBio reads processed for overlapping other PacBio reads
	UINT32 m_ProvOverlapping;               // number of PacBio reads overlapping at least one other PacBio read
	UINT32 m_ProvOverlapped;				// number of PacBio reads provisionally overlapped, could be containing, another PacBio read
	UINT32 m_ProvContained;					// number of PacBio reads provisionally contained within another PacBio read
	UINT32 m_ProvArtefact;					// number of PacBio reads provisionally only partially, likely an alignment artefact, contained within another PacBio read
	UINT32 m_ProvSWchecked;					// number of times SW used to identify overlaps

	UINT32 m_ExactKmerDists[100];			// accumulates exactly matching K-mer length distributions
	UINT64 m_TotAlignSeqLen;				// were over this total alignment length
	UINT32 m_TotAlignSeqs;					// between this number of sequence pairs

	UINT32 m_OverlapFloat;					// allow up to this much float on overlaps to account for the PacBio error profile

	UINT32 m_TranscriptomeLens;				// 0 if disabled, processing transcript reads, putatively overlapping reads must have length differential no more than this percentage and overlaps to be nearly full length
	UINT32 m_MinPBSeqLen;					// individual target PacBio sequences must be of at least this length
	UINT32 m_MaxPBRdSeqLen;					// and no longer than this length
	UINT32 m_MinPBSeqOverlap;				// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
	UINT32 m_MaxArtefactDev;				// classify overlaps as artefactual if sliding window of 500bp over any overlap deviates by more than this percentage from the overlap mean
	UINT32 m_MinHCSeqLen;					// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
	UINT32 m_MinHCSeqOverlap;				// any overlap of a hiconfidence read onto a target PacBio read must be of at least this many bp to be considered for contributing towards error correction (defaults to 1Kbp) 
	UINT8 m_HCRelWeighting;					// hiconfidence read overlaps are usually weighted higher than normal lesser confidence read overlaps when calling consensus bases  

	UINT32 m_MinErrCorrectLen;				// error corrected sequences must be at least this minimum length
	UINT32 m_MinConcScore;					// error corrected sequences trimmed until mean 100bp concensus score is at least this threshold

	UINT32 m_SampleRate;					// sample input sequences at this rate (1..100)
	bool m_bAntisenseOvlps;					// true if to process for both sense and antisense overlaps


	int m_NumPacBioFiles;					// number of input pacbio file specs
	char m_szPacBioFiles[cMaxInFileSpecs][_MAX_PATH];		// input pacbio files
	int m_NumHiConfFiles;					// number of input hiconfidence file specs
	char m_szHiConfFiles[cMaxInFileSpecs][_MAX_PATH];		// input hiconfidence files		

	char m_szErrCorFile[_MAX_PATH];			// name of file into which write error corrected and scored sequences
	int m_hErrCorFile;						// file handle for writing error corrected and scored sequences
	UINT32 m_ErrCorFileUnsyncedSize;		// count of chars which have been written to file handle m_hErrCorFile but may not be on disk

	char m_szChkPtsFile[_MAX_PATH];			// name of file used for checkpointing in case resume processing is required
	int m_hChkPtsFile;						// opened file handle for checkpointing


	int ErrCorBuffIdx;						// index into m_szErrCorLineBuff at which to next copy a corrected sequence
	int AllocdErrCorLineBuff;				// allocation size for m_pszErrCorLineBuff
	char *pszErrCorLineBuff;				// allocated buffering for error corrected and scored sequences
	int MultiAlignBuffIdx;					// index into m_pszMultiAlignLineBuff at which to write next multialignment
	int AllocdMultiAlignLineBuff;			// allocation size for m_pszMultiAlignLineBuff
	char *pszMultiAlignLineBuff;			// allocated buffering for multialignments


	char m_szMultiAlignFile[_MAX_PATH];		// name of file into which write multialignments
	int m_hMultiAlignFile;					// file handle for writing multialignments
	UINT32 m_MultiAlignFileUnsyncedSize;	// count of chars which have been written to file handle m_hMultiAlignFile but may not be on disk

	int m_ScaffLineBuffIdx;					// offset in m_szScaffLineBuff to write next overlap detail
	char m_szScaffLineBuff[0x07fff];		// buffering for overlap distributions

	UINT32 m_NumPBScaffNodes;					// m_pPBScaffNodes currently holds many scaffolding nodes
    UINT32 m_MaxPBSeqLen;						// max length of any scaffolding node sequence
	UINT32 m_AllocdPBScaffNodes;				// m_pPBScaffNodes allocated to hold this many scaffolding nodes
	tsPBEScaffNode *m_pPBScaffNodes;			// allocated to hold scaffolding nodes
	UINT32 *m_pMapEntryID2NodeIDs;				// used to map from suffix array entry identifiers to the corresponding scaffolding node identifier

	CSfxArrayV3 *m_pSfxArray;					// suffix array file (m_szTargFile) is loaded into this

	void Init(void);							// initialise state to that immediately following construction
	void Reset(void);						// reset state

	INT64 EstSumSeqLens(int NumTargFiles,char **pszTargFiles);		// guestimate and return total sequence length by simply summing the lengths of each file - likely to grossly over estimate

	int LoadSeqs(int MinSeqLen,					 // only accept for indexing sequences of at least this length
				 int MaxSeqLen,                  // and no longer than this length (bp)
				int NumTargFiles,char **pszTargFiles,	// parse, and index sequences in this file into memory resident suffix array; file expected to contain either fasta or fastq sequences
				int Flags = cFlgLCSeq);			// which by default are low confidence PacBio read sequences

	int ProcessBioseqFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
						  int MaxSeqLen,        // and no longer than this length (bp)
				 char *pszFile,					// file containing sequences
				int Flags = cFlgLCSeq);			// which by default are low confidence PacBio read sequences

	int ProcessFastaFile(int MinSeqLen,			// only accept for indexing sequences of at least this length (bp)
                int MaxSeqLen,                  // and no longer than this length (bp)
				char *pszFile,					// file containing sequences
				int Flags = cFlgLCSeq);			// which by default are low confidence PacBio read sequences

	int IdentifySequenceOverlaps(int MaxSeqLen,		// max length sequence to be overlapped
							int NumOvlpCores,		// targeting number of running threads to maximise usage of this many cores
							int NumOvlpThreads,		// identify all read overlaps using at most this this many threads ( max of m_MaxAllowedThreads and m_MaxRMIInstances)
							teBKSPType BKSPType  = eBKSPTSmithWaterman,		// workers are requesting this service type
							UINT32 RMIBufferSize =	cDfltRMIBufferSize,	// each worker thread default allocates to process up to this much buffered data
							UINT32 RMIParamDataSize = cDfltRMIReqParamSize, // each worker thread default allocates to process up to this much parameter data
							UINT32 RMIReqDataSize = cDfltRMIReqDataSize, // each worker thread default allocates to process up to this much request data
							UINT32 RMIRespDataSize =cDfltRMIRespDataSize); // each worker thread default allocates to process up to this much response data


	int IdentifyCoreHits(UINT32 ProbeNodeID,	// identify all overlaps of this probe sequence PBScaffNodeID onto target sequences
				tsThreadPBErrCorrect *pPars);		// thread specific

	int					// returns 0 if core overlapped (uses a non-exhaustive search) a previously added core, index 1..N of just added core hit or -1 if errors
		AddCoreHit(UINT32 ProbeNodeID,			// core hit was from this probe scaffold node 
			   bool bRevCpl,					// true if core sequence was revcpl'd before matching
			   UINT32 ProbeOfs,                 // hit started at this probe offset
			   UINT32 TargNodeID,               // probe core matched onto this target scaffold node
			   UINT32 TargOfs,                  // probe core matched starting at this target loci
			   UINT32 HitLen,					// hit was of this length
               tsThreadPBErrCorrect *pPars);	// thread specific

	UINT32										// returned tsPBScaffNode node identifier
		MapEntryID2NodeID(UINT32 EntryID);		// suffix array entry identifier

	CMTqsort m_mtqsort;				// muti-threaded qsort

static int SortLenDescending(const void *arg1, const void *arg2); // Sort scaffolding nodes by length descending
static int SortCoreHitsByTargProbeOfs(const void *arg1, const void *arg2);
static int SortCoreHitsByProbeTargOfs(const void *arg1, const void *arg2);
static int SortCoreHitsDescending(const void *arg1, const void *arg2);

	bool m_bMutexesCreated;			// will be set true if synchronisation mutexes have been created
	int CreateMutexes(void);
	void DeleteMutexes(void);

#ifdef WIN32
	alignas(4)	volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access - replaces AcquireSerialise() as much more efficient
	alignas(4)	volatile unsigned int m_CASLock; // used with synchronous compare and swap (CAS) for serialising access -  - replaces AcquireLock(True) as much more efficient
	alignas(4)	volatile unsigned int m_CASThreadPBErrCorrect; // used with synchronous compare and swap (CAS) for serialising access - replaces AcquireSerialise() as much more efficient
	alignas(4)	volatile UINT32 m_RMINumCommitedClasses;			// when RMI service provider classes utilised then current number of class instances instantiated over all provider sessions
	alignas(4)	volatile UINT32 m_RMINumUncommitedClasses;		// when RMI service provider classes utilised then current number of class instances available to be instantiated over all provider sessions
	alignas(4)	volatile UINT32 m_LowestCpltdProcNodeID;		// lowest processing completed sequence node
	alignas(4)	volatile UINT32 m_ReduceNonRMIThreads;			// reduce number of actively processing non-RMI SW processing threads by this delta
	alignas(4)	volatile UINT32 m_CurAllowedActiveOvlpThreads;	// currently allowing a total of this many SW (total of both RMI and non-RMI) threads to be actively processing, always <= m_MaxActiveOvlpThreads
	alignas(4)	volatile UINT32 m_CurActiveRMIThreads;			// currently there are this many actively processing RMI SW processing threads
	alignas(4)	volatile UINT32 m_CurActiveNonRMIThreads;		// currently there are this many actively processing non-RMI SW processing threads
#else
	__attribute__((aligned(4))) volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access - replaces AcquireSerialise() as much more efficient
	__attribute__((aligned(4))) volatile unsigned int m_CASLock; // used with synchronous compare and swap (CAS) for serialising access -  - replaces AcquireLock(True) as much more efficient
	__attribute__((aligned(4))) volatile unsigned int m_CASThreadPBErrCorrect; // used with synchronous compare and swap (CAS) for serialising access - replaces AcquireSerialise() as much more efficient
	__attribute__((aligned(4))) volatile UINT32 m_RMINumCommitedClasses;			// when RMI service provider classes utilised then current number of class instances instantiated over all provider sessions
	__attribute__((aligned(4))) volatile UINT32 m_RMINumUncommitedClasses;		// when RMI service provider classes utilised then current number of class instances available to be instantiated over all provider sessions
	__attribute__((aligned(4))) volatile UINT32 m_LowestCpltdProcNodeID;		// lowest processing completed sequence node
	__attribute__((aligned(4))) volatile UINT32 m_ReduceNonRMIThreads;			// reduce number of actively processing non-RMI SW processing threads by this delta
	__attribute__((aligned(4))) volatile UINT32 m_CurAllowedActiveOvlpThreads;	// currently allowing a total of this many SW (total of both RMI and non-RMI) threads to be actively processing, always <= m_MaxActiveOvlpThreads
	__attribute__((aligned(4))) volatile UINT32 m_CurActiveRMIThreads;			// currently there are this many actively processing RMI SW processing threads
	__attribute__((aligned(4))) volatile UINT32 m_CurActiveNonRMIThreads;		// currently there are this many actively processing non-RMI SW processing threads
#endif

	void AcquireCASSerialise(void);
	void ReleaseCASSerialise(void);
	void AcquireCASLock(void);
	void ReleaseCASLock(void);

	void AcquireCASThreadPBErrCorrect(void);
	void ReleaseCASThreadPBErrCorrect(void);


	time_t m_ProcessStatsThen;						// used to determine if sufficient time has elapsed since last reporting of process stats
	UINT32 m_NumCPUCores;                           // total number of CPU cores 
	UINT32 m_NumOvlpCores;							// processing is targeting this number of cores - will always be <= m_NumCPUCores
	UINT32 m_MaxRMIInstances;						// max number of RMI service provider instances allowed
	UINT32 m_MaxNonRMIThreads;						// max number of non-RMI SW threads allowed
	UINT32 m_MaxActiveOvlpThreads;					// allowing a maximum of this total number of SW (total of both RMI and non-RMI) threads to be actively processing


	int UpdateProcessStats(void);					// determine  numbers of commited and uncommited service provider classses

	int												// marshaled parameter required this many bytes
			MarshalReq(UINT8 *pInto,				// marshal into this list
					teRMIParamType Type,			// parameter type
					void *pValue,					// parameter value
					UINT32 ValLen);					// length of parameter ptd to by pValue, only used if parameter type is pUint8

	int
		UnmarshalResp(UINT32 DataLen,
					UINT8 *pFrom,		// unmarshal from this marshalled parameter list
					void *pValue);


	UINT64 RMI_new(tsThreadPBErrCorrect *pThreadPar,UINT32 Timeout,						// class instances constructor, returns the ClassInstanceID
											UINT32 SessionID = 0);						// requesting class instance on this specific session, 0 if on least loaded session
	
	void RMI_delete(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,UINT64 ClassInstanceID);	// delete the class instance

	bool RMI_SetScores(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
				int MatchScore= cSSWDfltMatchScore,			// score for match
				int MismatchPenalty  = cSSWDfltMismatchPenalty,	// penalty for mismatch
				int GapOpenPenalty  = cSSWDfltGapOpenPenalty,	// penalty for opening a gap
				int GapExtnPenalty  = cSSWDfltGapOpenPenalty,	// penalty if extending already opened gap
				int DlyGapExtn = cSSWDfltDlyGapExtn,			// delayed gap penalties, only apply gap extension penalty if gap at least this length
				int ProgPenaliseGapExtn = cSSWDfltProgPenaliseGapExtn,	// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles
				int AnchorLen = cSSWDfltAnchorLen);				// identified first and last anchors in alignment to be of at least this length


	bool RMI_SetCPScores(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
				int MatchScore= cSSWDfltMatchScore,		// ClassifyPath() score for match
				int MismatchPenalty  = cSSWDfltMismatchPenalty,	// ClassifyPath() penalty for mismatch
				int GapOpenPenalty  = cSSWDfltGapOpenPenalty,	// ClassifyPath() penalty for opening a gap
				int GapExtnPenalty  = cSSWDfltGapOpenPenalty);	// ClassifyPath() penalty if extending already opened gap

	bool RMI_SetMaxInitiatePathOfs(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
				int MaxInitiatePathOfs = cMaxInitiatePathOfs);	// require SW paths to have started within this many bp (0 to disable) on either the probe or target - effectively anchoring the SW 

	bool RMI_PreAllocMaxTargLen( tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
							 UINT32 MaxTargLen,					// preallocate to process targets of this maximal length
							 UINT32 MaxOverlapLen = 0);			// allocating tracebacks for this maximal expected overlap, 0 if no tracebacks required

	int
		RMI_StartMultiAlignments(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
					int SeqLen,						// probe sequence is this length
					etSeqBase *pProbeSeq,			// probe sequence 
					int Alignments,					// number of pairwise alignments to allocate for
					UINT8 Flags);					// bit 0 set true if probe sequence loaded as a high confidence sequence

	bool RMI_SetProbe(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID, 
					UINT32 Len,etSeqBase *pSeq);					// set probe sequence to use in subsequent alignments

	bool RMI_SetTarg(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,  
					UINT32 Len,etSeqBase *pSeq);					// set target sequence to use in subsequent alignments

	int RMI_SetAlignRange(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,  
						UINT32 m_ProbeStartRelOfs,	// when aligning then start SW from this probe sequence relative offset
					  UINT32 m_TargStartRelOfs, 	// and SW starting from this target sequence relative offset
						UINT32 m_ProbeRelLen = 0,	// and SW with this probe relative length starting from m_ProbeStartRelOfs - if 0 then until end of probe sequence
						UINT32 m_TargRelLen = 0);	// and SW with this target relative length starting from m_TargStartRelOfs - if 0 then until end of target sequence

		tsSSWCell *									// smith-waterman style local alignment, returns highest accumulated exact matches scoring cell
				RMI_Align(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID, 
						tsSSWCell *pPeakScoreCell = NULL,	// optionally also return conventional peak scoring cell
						UINT32 MaxOverlapLen = 0);		// process tracebacks for this maximal expected overlap, 0 if no tracebacks required

      // methods which combines the functionality of SetTarg, SetAlignRange, Align, ClassifyPath, TracebacksToAlignOps, and AddMultiAlignment into a single method 
	int					// -3: timeout waiting for job to complete, -2: parameter errors, -1: class instance no longer exists, 0: currently no available service instance 1: if job accepted and was processed
		RMI_CombinedTargAlign(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
						tsCombinedTargAlignPars *pAlignPars, // input alignment parameters
						tsCombinedTargAlignRet *pAlignRet);		// returned alignment results

	int					//  -3: timeout waiting for job to complete, -2: parameter errors, -1: class instance no longer exists, 0: currently no available service instance 1: if job accepted and processed
		RMI_CombinedTargAlign(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
						UINT8 PMode,              // processing mode: 0 error correct , 1 generate consensus from previously generated multiple alignments, 2  generate overlap detail from previously generated consensus sequences
						UINT32 NumTargSeqs,			// current probe putatively overlaying this many targets
						UINT32 ProbeSeqLen,         // use this probe sequence length 
						UINT8 TargFlags,		    // bit 7 set if target loaded as a high confidence sequence, bits 0..3 is weighting factor to apply when generating consensus bases
						UINT32 TargSeqLen,          // target sequence length
						etSeqBase *pTargSeq,        // target sequence
						UINT32 ProbeStartRelOfs,	// when aligning then start SW from this probe sequence relative offset
						UINT32 TargStartRelOfs, 	// and SW starting from this target sequence relative offset
						UINT32 ProbeRelLen,		    // and SW with this probe relative length starting from m_ProbeStartRelOfs - if 0 then until end of probe sequence
						UINT32 TargRelLen,		    // and SW with this target relative length starting from m_TargStartRelOfs - if 0 then until end of target sequence
						UINT32 OverlapFloat,		// allowing up to this much float on overlaps to account for the PacBio error profile
						UINT32 MaxArtefactDev,		// classify overlaps as artefactual if sliding window of 500bp over any overlap deviates by more than this percentage from the overlap mean
						UINT32 MinOverlapLen,       // minimum accepted overlap length
						UINT32 MaxOverlapLen,      // max expected overlap length
						UINT8 *pRetClass,			// returned overlap classification
						tsSSWCell *pRetPeakMatchesCell, // returned peak matches cell
						UINT32 *pRetProbeAlignLength, // probe alignment length
						UINT32 *pRetTargAlignLength, // target alignment length
						bool *pRetbProvOverlapping,  // probe overlapping target
						bool *pRetbProvArtefact,	// set true if overlap classified as artefact
						bool *pRetbProvContained,	// probe was contained
						bool *pRetbAddedMultiAlignment); // added as a multialignment


int												// attempting to determine if path is artfact resulting from aligning to a paralogous fragment
		RMI_ClassifyPath(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
					int MaxArtefactDev,			// classify path as artefactual if sliding window of 500bp over any overlap deviates by more than this percentage from the overlap mean
					UINT32 ProbeStartOfs,			// alignment starts at this probe sequence offset (1..n)
					UINT32 ProbeEndOfs,				// alignment ends at this probe sequence offset
					UINT32 TargStartOfs,			// alignment starts at this target sequence offset (1..n)
					UINT32 TargEndOfs);				// alignment ends at this target sequence offset

int												// number of alignment ops generated
		RMI_TracebacksToAlignOps(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
					UINT32 ProbeStartOfs,	// alignment starts at this probe sequence offset (1..n)
					UINT32 ProbeEndOfs,				// alignment ends at this probe sequence offset
					UINT32 TargStartOfs,			// alignment starts at this target sequence offset (1..n)
					UINT32 TargEndOfs,				// alignment ends at this target sequence offset
					tMAOp **ppAlignOps = NULL);     // optionally return ptr to alignment operations

int
		RMI_AddMultiAlignment(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
					  UINT32 ProbeStartOfs,			// alignment starts at this probe sequence offset (1..n)
					  UINT32 ProbeEndOfs,			// alignment ends at this probe sequence offset inclusive
					  UINT32 TargStartOfs,			// alignment starts at this target sequence offset (1..n)
					  UINT32 TargEndOfs,			// alignment ends at this target sequence offset inclusive
					  UINT32 TargSeqLen,			// target sequence length
					  etSeqBase *pTargSeq,			// alignment target sequence
    					UINT8 Flags);				// bit 7 set if target loaded as a high confidence sequence, bits 0..3 is weighting factor to apply when generating consensus bases


int
		RMI_GenMultialignConcensus(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID);

int      // total number of returned chars in pszBuffer for the textual representation of the error corrected consensus sequence (could be multiple consensus sequences)
		RMI_MAlignCols2fasta(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
					UINT32 ProbeID,	// identifies sequence which was used as the probe when determining the multialignments
				  int MinConf,				// sequence bases averaged over 100bp must be of at least this confidence (0..9)
				  int MinLen,				// and sequence lengths must be of at least this length 
				  UINT32 BuffSize,			// buffer allocated to hold at most this many chars
				  char *pszBuffer);			// output error corrected sequences to this buffer

int      // total number of returned chars in pszBuffer for the textual representation of the multialignment 
		RMI_MAlignCols2MFA(tsThreadPBErrCorrect *pThreadPar, UINT32 Timeout,  UINT64 ClassInstanceID,
						 UINT32 ProbeID,		// identifies sequence which was used as the probe when determining the multialignments
					    UINT32 BuffSize,	// buffer allocated to hold at most this many chars
					    char *pszBuffer);	// output multialignment textual representation to this buffer

public:
	CPBErrCorrect();
	~CPBErrCorrect();

	int ThreadPBErrCorrect(tsThreadPBErrCorrect *pThreadPar);

	int
	Process(etPBPMode PMode,		// processing mode
		char *pszHostName,			// listening on this host name or IPv4/IPv5 address for connections by service providers 
		char *pszServiceName,			// Listen on this service name or port for for connections by service providers
		int MaxRMI,					// max number of RMI SW service provider instances supported
		int MaxNonRMI,				// max number of non-RMI SW threads supported
		int SampleRate,				// sample input sequences at this rate (1..100)
		bool bAntisenseOvlps,		// true if to process for both sense and antisense overlaps
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int TranscriptomeLens,		// 0 if disabled, processing transcript reads, putatively overlapping reads must have length differential no more than this percentage and overlaps to be nearly full length
		int MinPBSeqLen,			// only accepting PacBio reads of at least this length (defaults to 10Kbp) and if
		int MaxPBSeqLen,			// no more than this length (defaults to 30Kbp)
		int MinPBSeqOverlap,		// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
		int MaxArtefactDev,			// classify overlaps as artefactual if sliding window of 500bp over any overlap deviates by more than this percentage from the overlap mean
		int MinHCSeqLen,			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
		int MinHCSeqOverlap,		// any overlap of a hiconfidence read onto a target PacBio read must be of at least this many bp to be considered for contributing towards error correction (defaults to 1Kbp) 
		int HCRelWeighting,             // hiconfidence read overlaps are usually weighted higher than normal lesser confidence read overlaps when calling consensus bases 
	    int MinErrCorrectLen,		// error corrected and trimmed sequences must be at least this minimum length
		int MinConcScore,			// error corrected sequences trimmed until mean 100bp concensus score is at least this threshold
		int NumPacBioFiles,			// number of input pacbio file specs
		char *pszPacBioFiles[],		// input pacbio files
		int NumHiConfFiles,			// number of input hiconfidence file specs
		char *pszHiConfFiles[],		// input hiconfidence files		
	    char *pszErrCorFile,		// name of file into which write error corrected sequences
		char *pszMultiAlignFile,	// name of file into which write multiple alignments
		char *pszChkPtsFile,        // name of file used for checkpointing in case resume processing is required
		int NumThreads);			// maximum number of worker threads to use

		int
		GenConsensusFromMAF(int MinErrCorrectLen,		// error corrected sequences must be at least this minimum length
					 int MinConcScore,			// error corrected sequences trimmed until mean 100bp concensus score is at least this threshold
					char *pszErrCorFile,		// name of file into which write error corrected sequences
					char *pszMultiAlignFile);	// name of file containing multiple alignments to process

};




