#pragma once
/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "../libbiokanga/commdefs.h"
#include "SSW.h"

const int cAllocTargetSeqSize = 100000;				// allocation buffer size to hold target sequences of at least this length
const int cAllocTProbetSeqSize = 100000;			// allocation buffer size to hold probe sequences of at least this length
const int cDfltMinTargetSeqLen    = 1000;			// minimum target sequence length accepted  
const int cDfltMinTargetOvlpLen   = 500;			// requiring target sequences to end overlap by at least this many bp or else be fully contained or containing  
const int cDfltMaxArtefactDev = 75;					// default percentage deviation from the mean over a 1Kbp alignment window allowed when classing overlaps as being artefactual
const int cDfltSSWDInitiatePathOfs = 250;			// default is to require SW paths to have started within this many bp on either the probe or target - effectively anchoring the SW
const int cDfltSSWDOvlpFloat = 100;					// allowed float in bp on SW overlaps

const int cMaxProbePBSSWs = 100;						// explore with SW at most this many probe alignments against target sequences
const int cPBSSWSummaryTargCoreHitCnts = 100;		// summary core hit counts on at most this many targets

#pragma pack(1)

typedef enum TAG_ePBSSWOverlapClass {
	ePBSSWOLCOverlapping = 0,	// probe classified as overlapping target, either 5' or 3'
	ePBSSWOLCcontains,			// probe completely contains the target
	ePBSSWOLCcontained,			// probe is completely contained within the target
	ePBSSWOLCartefact			// probe contains a subsequence of target, classifying as an artefact overlap and not further processed
} ePBSSWOverlapClass;

// seed core hits 
typedef struct TAG_sPBSSWCoreHit {
	UINT32 TargSeqID;				// hit was onto this target  node
	UINT32 ProbeOfs;                // hit was from this probe offset
	UINT32 TargOfs;					// onto this target offset
	UINT32 HitLen;					// hit was of this length
	UINT32 WinHits;					// number of core hits onto target relative to this core which are within a window of probelen
	UINT8 flgRevCpl:1;				// 1 if core sequence was revcpl'd before matching
	UINT8 flgMulti:1;				// 1 if core sequence was target multiloci and this instance to not be further processed
	UINT8 flgClustered:1;			// 1 if core hit identified as part of a cluster of hits
	} tsPBSSWACoreHit;

typedef struct TAG_sPBSSWCoreHitCnts {
	UINT32 TargSeqID;				// node identifier for hit target sequence	
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
} sPBSSWCoreHitCnts;

typedef struct TAG_sPBSSWInstance {
	UINT32 InstanceID;					// each instance is uniquely identified
	UINT8  FlgActive:1;                 // 1 if this instance is currently actively performing an alignment, new alignments can only be initiated if FlgActive == 0
	UINT32 NumCoreHits;					// currently this many core hits in m_pCoreHits
	UINT32 AllocdCoreHits;				// m_pCoreHits currently allocated to hold at most this many core hits
	size_t AllocdCoreHitsSize;			// m_pCoreHits current allocation size
	tsPBSSWACoreHit *pCoreHits;			// allocated to hold all core hits	
	UINT32 NumTargCoreHitCnts;			// current number of summary target core hit counts in TargCoreHitCnts
	sPBSSWCoreHitCnts TargCoreHitCnts[cPBSSWSummaryTargCoreHitCnts]; // top targets by core hit counts

	UINT32 ProbeSeqLen;					// current probe sequence length
	UINT32 AllocdProbeSeqSize;			// current allocation size for buffered probe sequence in pProbeSeq 	
	etSeqBase *pProbeSeq;				// allocated to hold the current probe sequence

	UINT32 TargSeqLen;					// current target sequence length
	UINT32 AllocdTargSeqSize;			// current allocation size for buffered target sequence in pTargSeq 	
	etSeqBase *pTargSeq;				// allocated to hold the current target sequence

	CMTqsort *pmtqsort;					// muti-threaded qsort
	CSSW *pSW;							// Smith-waterman class instance used by this alignment instance
} tsPBSSWInstance;

#pragma pack()

class CSWAlign
{

	int m_SWMatchScore;						// SW score for matching bases (0..100)
	int m_SWMismatchPenalty;				// SW mismatch penalty (-100..0)
	int m_SWGapOpenPenalty;					// SW gap opening penalty (-100..0)
	int m_SWGapExtnPenalty;					// SW gap extension penalty (-100..0)
	int m_SWProgExtnPenaltyLen;				// only apply gap extension penalty if gap at least this length (1..63) - use if aligning PacBio
	int m_SWDlyGapExtn;						// delayed gap penalties, only apply gap extension penalty if gap at least this length
	int m_SWProgPenaliseGapExtn;			// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles
	int m_SWAnchorLen;						// identified first and last anchors in alignment to be of at least this length

	int m_CPMatchScore;						// path class, ClassifyPath(), score for matching bases (0..100)
	int m_CPMismatchPenalty;				// path class mismatch penalty (-100..0)
	int m_CPGapOpenPenalty;					// path class gap opening penalty (-100..0)
	int m_CPGapExtnPenalty;					// path class gap extension penalty (-100..0)
	int m_MaxInitiatePathOfs;				// if non-zero then only allow new paths to start if within that offset (0 to disable) on either probe or target - effectively an anchored SW

	UINT32 m_MinOverlapLen;					// the putative overlap would be of at least this length
	UINT32 m_MaxSeedCoreDepth;				// only further extend a seed core if there are no more than this number of matching cores in all targeted sequences
	UINT32 m_DeltaCoreOfs;					// offset core windows of coreSeqLen along the probe sequence when checking for overlaps 
	UINT32 m_SeedCoreLen;					// putative overlaps are explored if there are cores of at least this length in any putative overlap
	UINT32 m_MinNumCores;					// and if the putative overlap contains at least this many cores
	UINT32 m_MaxAcceptHitsPerSeedCore;		// limit accepted hits per seed core to no more this many
	UINT32 m_MinNumSeedCores;				// require at least this many seed cores between overlapping scaffold sequences
	UINT32 m_OverlapFloat;					// allow up to this much float on overlaps to account for the PacBio error profile
	UINT32 m_MinPBSeqLen;					// individual target PacBio sequences must be of at least this length
	UINT32 m_MinPBSeqOverlap;				// any overlap of a PacBio onto a target PacBio must be of at least this many bp to be considered for contributing towards error correction (defaults to 5Kbp) 
	UINT32 m_MaxArtefactDev;				// classify overlaps as artefactual if sliding window of 1Kbp over any overlap deviates by more than this percentage from the overlap mean



	UINT32 m_DfltMaxProbeSeqLen;	// initially allocate for this length probe sequence to be aligned, will be realloc'd as may be required

	UINT32 m_NumSWInstances;			// number of SW instances currently used
	UINT32 m_AllocdSWInstances;			// allocated to hold at most this many instances
	tsPBSSWInstance *m_pSWInstances;	// allocated to hold SW instances

	UINT32 m_MaxTargSeqLen;				// longest target sequence length loaded
	CSfxArrayV3 *m_pSfxArray;			// targeted sequences are loaded/indexed into this suffix array

	int					// returns 0 if core overlapped (uses a non-exhaustive search) a previously added core, index 1..N of just added core hit or -1 if errors
	AddCoreHit(tsPBSSWInstance *pInstance,		// using this instance
				bool bRevCpl,					// true if core sequence was revcpl'd before matching
			   UINT32 ProbeOfs,                 // hit started at this probe offset
			   UINT32 TargSeqID,               // probe core matched onto this target sequence
			   UINT32 TargOfs,                  // probe core matched starting at this target loci
			   UINT32 HitLen);					// hit was of this length

	int									// returns index 1..N of core hits remaining or -1 if errors
	RemoveAddedCoreHits(tsPBSSWInstance *pInstance,		// using this instance
				int NumToRemove);                   // removing the last NumToRemove AddCoreHit() added

	int
	IdentifyCoreHits(tsPBSSWInstance *pInstance,// using this instance
					 bool bRevCpl,			// true if probe sequence to be reverse complemented
					 UINT32 MinTargLen = 1,		// minimum accepted target length
					 UINT32 MaxTargLen = 0);		// maximum accepted target length

	int LoadTargFastaFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
				char *pszFile,					// file containing sequences
				int Flags = 0);					// flags are user defined

	bool m_bLocksCreated;					// set true when locks created
	int CreateLocks(void);
	void AcquireLock(bool bExclusive);
	void ReleaseLock(bool bExclusive);
	void DeleteLocks(void);

#ifdef _WIN32
	SRWLOCK m_hRwLock;
#else
	pthread_rwlock_t m_hRwLock;
#endif

static int SortCoreHitsByTargProbeOfs(const void *arg1, const void *arg2);
static int SortCoreHitsByProbeTargOfs(const void *arg1, const void *arg2);
static int SortCoreHitsDescending(const void *arg1, const void *arg2);

public:
	CSWAlign();
	~CSWAlign();

	void Reset(void);
	int Initialise(UINT32 MaxSWAInstances,	// initialise for this many instances
			UINT32 MinOverlapLen,			// the putative overlap would be of at least this length
			UINT32 MaxSeedCoreDepth,		// only further extend a seed core if there are no more than this number of matching cores in all targeted sequences
			UINT32 DeltaCoreOfs,			// offset core windows of coreSeqLen along the probe sequence when checking for overlaps 
			UINT32 CoreSeqLen,				// putative overlaps are explored if there are cores of at least this length in any putative overlap
			UINT32 MinNumCores,				// and if the putative overlap contains at least this many cores
			UINT32 MaxAcceptHitsPerSeedCore, // limit accepted hits per seed core to no more this many
			UINT32 DfltMaxProbeSeqLen);		// initially allocate for this length probe sequence to be aligned, will be realloc'd as may be required


	int									// returns number of target sequences loaded and indexed, or error result if < 0
		LoadTargetSeqs(int MinSeqLen,	// only accepting target sequences of at least this length
			char *pszTargSeqsFile,		// load target sequences from this file
			int NumThreads = 4);		// use at most this number of threads when indexing target sequences

	UINT32			// returned alignment instance identifier or 0 if errors
			InitInstance(void);        // initialise instance - uses scores etc., as previously set with SetScores() and/or SetCPScores() and/or SetMaxInitiatePathOfs()

	bool SetScores(int MatchScore = cSSWDfltMatchScore,			// score for match
				   int MismatchPenalty = cSSWDfltMismatchPenalty,	// penalty for mismatch
				   int GapOpenPenalty = cSSWDfltGapOpenPenalty,	// penalty for opening a gap
				   int GapExtnPenalty = cSSWDfltGapExtnPenalty,	// penalty if extending already opened gap
				   int DlyGapExtn = cSSWDfltDlyGapExtn,			// delayed gap penalties, only apply gap extension penalty if gap at least this length
				   int ProgPenaliseGapExtn = cSSWDfltProgPenaliseGapExtn,	// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles
				   int AnchorLen = cSSWDfltAnchorLen);				// identified first and last anchors in alignment to be of at least this length

	bool SetCPScores(int MatchScore = cSSWDfltMatchScore,		// ClassifyPath() score for match
					 int MismatchPenalty = cSSWDfltMismatchPenalty,	// ClassifyPath() penalty for mismatch
					 int GapOpenPenalty = cSSWDfltGapOpenPenalty,	// ClassifyPath() penalty for opening a gap
					 int GapExtnPenalty = cSSWDfltGapExtnPenalty);	// ClassifyPath() penalty if extending already opened gap

	bool SetMaxInitiatePathOfs(int MaxInitiatePathOfs = cDfltSSWDInitiatePathOfs,	// require SW paths to have started within this many bp (0 to disable) on either the probe or target - effectively anchoring the SW 
								int OverlapFloat = cDfltSSWDOvlpFloat);		// with this overlap float 

	int
			AlignProbeSeq(UINT32 SWAInstance,			// using this alignment instance
							UINT32 ProbeSeqLen,			// sequence to align is this length
							etSeqBase *pProbeSeq,      // probe sequence to align
							UINT32 MinTargLen = 1,          // aligned to targets must be at least this long
							UINT32 MaxTargLen = 0,			// and if > 0 then target no longer than this many bp
							bool bSenseOnly = false,   // true if to align probe sense only, false to align both sense and antisense
	    					tsSSWCell *pRetMatched = NULL);    // optional (if not NULL) returned match detail
		
};

