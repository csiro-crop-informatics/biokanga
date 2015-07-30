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
const int cMinSeqLen = 500;					// minimum sequence length accepted for processing
const int cMaxDupEntries = 10;				// report 1st 10 duplicate entry names
const int cMaxAllocBuffChunk = 0x00ffffff;	// buffer for fasta sequences is realloc'd in this sized chunks

const int cFlgLCSeq = 0x01;					// sequence is a low confidence read (PacBio) containing many sequence errors
const int cFlgHCSeq = 0x02;					// sequence a high confidence sequence with few sequence errors

const int cDfltMinHCSeqLen = 1000;			// default min high confidence sequence length
const int cDfltMinHCSeqOverlap = 1000;		// default min high confidence sequence overlap length onto PacBio sequence

const int cMinMaxArtefactDev = 1;			// user can specify down to this minimum or 0 to disable
const int cDfltMaxArtefactDev = 50;			// default percentage deviation from the mean allowed when classing overlaps as being artefactual
const int cDfltScaffMaxArtefactDev = 15;		// but when scaffolding with error corrected reads there should a much smaller percentage deviation from the mean
const int cMaxMaxArtefactDev = 70;			// user can specify up to this maximum 


const int cMaxPacBioErrCorLen = 250000;		// allowing for error corrected read sequences of up to this length
const int cMaxPacBioMAFLen = (cMaxPacBioErrCorLen * 100);	// allowing for multialignment format buffering of up to this length


typedef enum TAG_ePBPMode {								// processing mode
	ePBPMErrCorrect,									// error correct
	ePBPMConsensus,										// generate consensus from previously generated multiple alignments
	ePBPMOverlapDetail									// generate overlap detail from previously generated consensus sequences
	} etPBPMode;


#pragma pack(1)

// seed core hits 
typedef struct TAG_sPBECoreHit {
	UINT32 ProbeNodeID;				// core hit was from this probe  node 
	UINT32 TargNodeID;				// hit was onto this target  node
	UINT32 ProbeOfs;                // hit was from this probe offset
	UINT32 TargOfs;					// onto this target offset
	UINT32 HitLen;					// hit was of this length
	UINT8 flgRevCpl:1;				// 1 if core sequence was revcpl'd before matching
	UINT8 flgMulti:1;				// 1 if core sequence was target multiloci and this instance to not be further processed
	} tsPBECoreHit;

typedef struct TAG_sPBEScaffNode {
	UINT32 NodeID;					// uniquely identifies this node
	UINT32 VertexID;				// assembly graph vertex identifier
	UINT32 EntryID;					// suffix array entry identifier for indexed sequence
	UINT32 SeqLen;					// length in bp of this scaffolding node sequence
	UINT8 flgCurProc:1;			// sequence is currently being processed
	UINT8 flgContained:1;			// sequence is fully contained within at least one other sequence
	UINT8 flgContains:1;			// sequence fully contains at least one other sequence
	UINT8 flgUnderlength:1;        // sequence is under length
	UINT8 flgHCseq:1;				// loaded as a high confidence (non-PacBio) sequence
} tsPBEScaffNode;

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
	UINT32 MinPropBinned;			// and if the putative overlap contains at least this proportion (1..100) of 250bp bins binned cores
	UINT32 MaxAcceptHitsPerSeedCore; // limit accepted hits per seed core to no more this many
	UINT32 MinPBSeqLen;				// only process PacBio sequences which are at least this length

	UINT32 NumTargCoreHitCnts;		// current number of summary target core hit counts in TargCoreHitCnts
	sPBECoreHitCnts TargCoreHitCnts[cSummaryTargCoreHitCnts]; // top targets by core hit counts
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
} tsThreadPBErrCorrect;


#pragma pack()


class CPBErrCorrect
{
	etPBPMode m_PMode;						// processing mode

	UINT32 m_DeltaCoreOfs;					// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
	UINT32 m_MaxSeedCoreDepth;				// only further process a seed core if there are no more than this number of matching cores in all targeted sequences

	UINT32 m_MinSeedCoreLen;				// use seed cores of this length when identifying putative overlapping scaffold sequences
	UINT32 m_MinNumSeedCores;				// require at least this many seed cores between overlapping scaffold sequences
	UINT32 m_BinClusterSize;				// clustering seed cores into this sized bins when determing if too few bins with at least 1 core; these few bins likely to result in SW artefacts 
	UINT32 m_MinPropBinned;					// require that the putative overlap contains at least this proportion (1..100) of m_BinClusterSize clustered binned cores


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
	UINT32 m_MinPBSeqLen;					// individual target PacBio sequences must be of at least this length
	UINT32 m_MinPBSeqOverlap;				// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
	UINT32 m_MaxArtefactDev;				// classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage from the overlap mean
	UINT32 m_MinHCSeqLen;					// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
	UINT32 m_MinHCSeqOverlap;				// any overlap of a hiconfidence read onto a target PacBio read must be of at least this many bp to be considered for contributing towards error correction (defaults to 1Kbp) 

	UINT32 m_MinErrCorrectLen;				// error corrected sequences must be at least this minimum length
	UINT32 m_MinConcScore;					// error corrected sequences trimmed until mean 100bp concensus score is at least this threshold

	UINT32 m_SampleRate;					// sample input sequences at this rate (1..100)

	int m_NumPacBioFiles;					// number of input pacbio file specs
	char m_szPacBioFiles[cMaxInFileSpecs][_MAX_PATH];		// input pacbio files
	int m_NumHiConfFiles;					// number of input hiconfidence file specs
	char m_szHiConfFiles[cMaxInFileSpecs][_MAX_PATH];		// input hiconfidence files		

	char m_szErrCorFile[_MAX_PATH];			// name of file into which write error corrected and scored sequences
	int m_hErrCorFile;						// file handle for writing error corrected and scored sequences

	int ErrCorBuffIdx;						// index into m_szErrCorLineBuff at which to next copy a corrected sequence
	int AllocdErrCorLineBuff;				// allocation size for m_pszErrCorLineBuff
	char *pszErrCorLineBuff;				// allocated buffering for error corrected and scored sequences
	int MultiAlignBuffIdx;					// index into m_pszMultiAlignLineBuff at which to write next multialignment
	int AllocdMultiAlignLineBuff;			// allocation size for m_pszMultiAlignLineBuff
	char *pszMultiAlignLineBuff;			// allocated buffering for multialignments


	char m_szMultiAlignFile[_MAX_PATH];		// name of file into which write multialignments
	int m_hMultiAlignFile;					// file handle for writing multialignments

	int m_ScaffLineBuffIdx;					// offset in m_szScaffLineBuff to write next overlap detail
	char m_szScaffLineBuff[0x07fff];		// buffering for overlap distributions

	int m_NumThreads;							// maximum number of worker threads to use

	UINT32 m_NumPBScaffNodes;					// m_pPBScaffNodes currently holds many scaffolding nodes
    UINT32 m_MaxPBSeqLen;						// max length of any scaffolding node sequence
	UINT32 m_AllocdPBScaffNodes;				// m_pPBScaffNodes allocated to hold this many scaffolding nodes
	tsPBEScaffNode *m_pPBScaffNodes;				// allocated to hold scaffolding nodes
	UINT32 *m_pMapEntryID2NodeIDs;				// used to map from suffix array entry identifiers to the corresponding scaffolding node identifier

	CSfxArrayV3 *m_pSfxArray;					// suffix array file (m_szTargFile) is loaded into this

	void Init(void);							// initialise state to that immediately following construction
	void Reset(bool bSync);						// reset state, if bSync true then fsync before closing output file handles

	INT64 EstSumSeqLens(int NumTargFiles,char **pszTargFiles);		// guestimate and return total sequence length by simply summing the lengths of each file - likely to grossly over estimate

	int LoadSeqs(int MinSeqLen,int NumTargFiles,char **pszTargFiles,	// parse, and index sequences in this file into memory resident suffix array; file expected to contain either fasta or fastq sequences
				int Flags = cFlgLCSeq);			// which by default are low confidence PacBio read sequences

	int ProcessBioseqFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
				 char *pszFile,					// file containing sequences
				int Flags = cFlgLCSeq);			// which by default are low confidence PacBio read sequences

	int ProcessFastaFile(int MinSeqLen,			// only accept for indexing sequences of at least this length
				char *pszFile,					// file containing sequences
				int Flags = cFlgLCSeq);			// which by default are low confidence PacBio read sequences

	int IdentifySequenceOverlaps(int MaxSeqLen,		// max length sequence to be overlapped
							int NumOvlpThreads);	// identify all overlaps using this many threads

	int IdentifyCoreHits(UINT32 ProbeNodeID,	// identify all overlaps of this probe sequence PBScaffNodeID onto target sequences
				tsThreadPBErrCorrect *pPars);		// thread specific

	int					// returns index 1..N of just added core hit or -1 if errors
		AddCoreHit(UINT32 ProbeNodeID,			// core hit was from this probe scaffold node 
			   bool bRevCpl,					// true if core sequence was revcpl'd before matching
			   UINT32 ProbeOfs,                 // hit started at this probe offset
			   UINT32 TargNodeID,               // probe core matched onto this target scaffold node
			   UINT32 TargOfs,                  // probe core matched starting at this target loci
			   UINT32 HitLen,					// hit was of this length
               tsThreadPBErrCorrect *pPars);		// thread specific

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

	volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access - replaces AcquireSerialise() as much more efficient
	void AcquireCASSerialise(void);
	void ReleaseCASSerialise(void);
	volatile unsigned int m_CASLock; // used with synchronous compare and swap (CAS) for serialising access -  - replaces AcquireLock(True) as much more efficient
	void AcquireCASLock(void);
	void ReleaseCASLock(void);


public:
	CPBErrCorrect();
	~CPBErrCorrect();

	int ThreadPBErrCorrect(tsThreadPBErrCorrect *pThreadPar);

	int
	Process(etPBPMode PMode,		// processing mode
		int SampleRate,				// sample input sequences at this rate (1..100)
		int DeltaCoreOfs,			// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
		int MaxSeedCoreDepth,		// only further process a seed core if there are no more than this number of matching cores in all targeted sequences
		int MinSeedCoreLen,			// use seed cores of this length when identifying putative overlapping scaffold sequences
		int MinNumSeedCores,        // require at least this many seed cores between overlapping scaffold sequences
		int SWMatchScore,			// score for matching bases (0..50)
		int SWMismatchPenalty,		// mismatch penalty (-50..0)
		int SWGapOpenPenalty,		// gap opening penalty (-50..0)
		int SWGapExtnPenalty,		// gap extension penalty (-50..0)
		int SWProgExtnPenaltyLen,	// progressive gap scoring then only apply gap extension score if gap at least this length (0..63) - use if aligning PacBio
		int MinPBSeqLen,			// only accepting PacBio reads of at least this length (defaults to 5Kbp)
		int MinPBSeqOverlap,		// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
		int MaxArtefactDev,			// classify overlaps as artefactual if sliding window of 500bp over any overlap deviates by more than this percentage from the overlap mean
		int MinHCSeqLen,			// only accepting hiconfidence reads of at least this length (defaults to 1Kbp)
		int MinHCSeqOverlap,		// any overlap of a hiconfidence read onto a target PacBio read must be of at least this many bp to be considered for contributing towards error correction (defaults to 1Kbp) 
		int MinErrCorrectLen,		// error corrected and trimmed sequences must be at least this minimum length
		int MinConcScore,			// error corrected sequences trimmed until mean 100bp concensus score is at least this threshold
		int NumPacBioFiles,			// number of input pacbio file specs
		char *pszPacBioFiles[],		// input pacbio files
		int NumHiConfFiles,			// number of input hiconfidence file specs
		char *pszHiConfFiles[],		// input hiconfidence files		
	    char *pszErrCorFile,		// name of file into which write error corrected sequences
		char *pszMultiAlignFile,	// name of file into which write multiple alignments
		int NumThreads);			// maximum number of worker threads to use

		int
		GenConsensusFromMAF(int MinErrCorrectLen,		// error corrected sequences must be at least this minimum length
					 int MinConcScore,			// error corrected sequences trimmed until mean 100bp concensus score is at least this threshold
					char *pszErrCorFile,		// name of file into which write error corrected sequences
					char *pszMultiAlignFile);	// name of file containing multiple alignments to process

};



