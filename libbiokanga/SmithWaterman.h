#pragma once
#include "./commdefs.h"

const UINT32 cSWMinProbeOrTargLen = 5;			// require probe and target lengths to be at least this number of bp
const UINT32 cSWMinCells = (cSWMinProbeOrTargLen * cSWMinProbeOrTargLen);	// not worth the effort if < this number of cells! minimum query and target length is 5bp
const UINT32 cSWMaxProbeOrTargLen = 1000000;  // require either probe or target length to be no longer than this limit
const UINT64 cSWMaxCells = ((UINT64)cSWMaxProbeOrTargLen * cSWMaxProbeOrTargLen/10);	// 10^10 is getting big, big! lets hope this limit is never reached

// default smith-waterman scores
const int cSWDfltMatchScore = 1;		// score for matching bases
const int cSWDfltMismatchPenalty = -1;	// mismatch penalty
const int cSWDfltGapOpenPenalty = -3;	// gap opening penalty
const int cSWDfltGapExtnPenalty = -1;	// gap extension penalty
const int cSWDfltDlyGapExtn = 2;        // apply gap extension penalty for gaps of at least 2bp (must be in range 1..63)
const int cSWDfltProgPenaliseGapExtn = 0; // if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles


typedef UINT32 tSWTrcBckCell;		// traceback cells are 32 bits
const  UINT32 cSWScoreMsk   =     0x001fffff; // score is clamped to 2^21-1  max
const UINT32  cSWInDelLenMsk   =     0x07e00000; // these 6 bits used to hold the current gap opened length (clamped to be no more than 63)
const int cSWInDelLenShf = 21;						// shift factor to move InDel length into LSBs or from LSBs into cSWInDelLenMsk
const  UINT32 cSWGapOpnFlg  =     0x08000000; // set in cell if gap has been opened
const  UINT32 cSWTrcBckMatchFlg = 0x10000000; // set if probe and target base were exactly matching
const  UINT32 cSWTrcBckMsk  =     0xE0000000; // 3 bits to hold traceback direction (0=none,001=diag,010=down,100=left)
const  UINT32 cSWTrcBckDiagFlg =  0x20000000; // traceback diagonally, bases align, could be exact match (cSWTrcBckMatchFlg set) or mismatched bases
const  UINT32 cSWTrcBckDownFlg =  0x40000000; // traceback down, base insert into target (or deletion from probe)
const  UINT32 cSWTrcBckLeftFlg =  0x80000000; // traceback left, base insert into probe (or deletion from target)

#pragma pack(1)
typedef struct TAG_sSWColBand {
	UINT32 TrcBckCellPsn;				// cells in this band have been allocated from cells starting at m_pTrcBckCells[TrcBckCellPsn-1], 0 if none allocated
	UINT32 StartTargBasePsn;			// cells in this probe column start at this target base position 1..m_TargLen
	UINT32 EndTargBasePsn;			    // cells in this probe column end at this target  1..m_TargLen 
} tsSWColBand;
#pragma pack()

class CSmithWaterman
{
	bool m_bAligned;					// set true following successful Needleman-Wunsch scoring alignment
	bool m_bBanded;						// set true if banded processing
	tSWTrcBckCell *m_pTrcBckCells;		// array used for holding traceback cells
	size_t m_TrcBckCellsAllocdSize;		// allocation size
	size_t m_TrcBckCellsAllocd;			// actual cells in alloc'd m_pTrcBckCells (in tSWTrcBckCell's)
	size_t m_TrcBckCellsUsed;			// current cells used

	UINT32 m_ColBandsUsed;				// number of column active bands used - should be same as m_ProbeLen
	UINT32 m_ColBandsAllocd;			// number of column active bands allocated
	tsSWColBand *m_pColBands;			// each column, or probe base, in the matrix will contain a band of active target base cells

	etSeqBase *m_pProbe;				// alloc'd probe sequence memory
	UINT32 m_ProbeAllocd;				// actual m_pProbe alloc'd size (in teSWSeqBases's)
	UINT32 m_ProbeLen;					// current probe length

	etSeqBase *m_pTarg;					// alloc'd target sequence memory
	UINT32 m_TargAllocd;					// actual m_pTarg alloc'd size (in teSWSeqBases's)
	UINT32 m_TargLen;						// current targ length

	int m_MatchScore;					// score for matching bases (0..100)
	int m_MismatchPenalty;				// mismatch penalty (-100..0)
	int m_GapOpenPenalty;				// gap opening penalty (-100..0)
	int m_GapExtnPenalty;				// gap extension penalty (-100..0)
	int m_DlyGapExtn;					// delayed gap penalties, only apply gap extension penalty if gap at least this length
    int m_ProgPenaliseGapExtn;			// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles

	UINT32 m_SWBandInitial;				// expecting initial overlap for max scoring alignment path to have non-overlap flank of at most this many bases
	double m_SWPathLenDiff;				// restrain path length differential between query and probe for max score be at most this proportion of query

	int m_PeakScore;					// peak score in any cell
	UINT32 m_NumBasesAligned;				// this many bases (exact and subs) were aligned between probe and target
	UINT32 m_NumBasesExact;				// of the m_NumBasesAligned, this many were exact matches, remainder were subs
	UINT32 m_NumProbeInserts;				// number of inserted bases into the probe relative to the target
	UINT32 m_NumTargInserts;				// number of inserted bases into the target relative to the probe

	UINT32 m_ProbeAlignStartOfs;			// offset in the probe at which the alignment starts
	UINT32 m_TargAlignStartOfs;			    // offset in the target at which the alignment starts

	UINT32 m_PeakProbeIdx;					// highest scoring cell is at this matrix column or probe ofs
	UINT32 m_PeakTargIdx;					// highest scoring cell is at this matrix row or target ofs 

	inline tSWTrcBckCell *						// returns ptr to newly allocated cell  or NULL if errors
			AllocBandCell(UINT32 ProbeBasePsn,			// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn);		// current target base position 1..m_TargLen


	inline tSWTrcBckCell *						// returns ptr to cell 
			DerefBandCell(UINT32 ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn);		// current target base position 1..m_TargLen

	inline tSWTrcBckCell *					// returns ptr to cell  or NULL if errors
			DerefBandCellLeft(UINT32 ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn);		// current target base position 1..m_TargLen


	inline tSWTrcBckCell *					// returns ptr to cell  or NULL if errors
			DerefBandCellDiag(UINT32 ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn);		// current target base position 1..m_TargLen

	inline tSWTrcBckCell *									// returns ptr to cell  or NULL if errors
			DerefBandCellDown(UINT32 ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn);		// current target base position 1..m_TargLen

	


public:
	CSmithWaterman(void);
	~CSmithWaterman(void);

	void Reset(void);

	bool SetScores(int MatchScore= cSWDfltMatchScore,			// score for match
				int MismatchPenalty  = cSWDfltMismatchPenalty,	// penalty for mismatch
				int GapOpenPenalty  = cSWDfltGapOpenPenalty,	// penalty for opening a gap
				int GapExtnPenalty  = cSWDfltGapOpenPenalty,	// penalty if extending already opened gap
				int DlyGapExtn = cSWDfltDlyGapExtn,				// delayed gap penalties, only apply gap extension penalty if gap at least this length
				int ProgPenaliseGapExtn = cSWDfltProgPenaliseGapExtn);	// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles


	bool SetProbe( UINT32 Len,etSeqBase *pSeq);					// set probe sequence to use in subsequent alignments
	bool SetTarg( UINT32 Len,etSeqBase *pSeq);					// set target sequence to use in subsequent alignments

	int // smith-waterman style local alignment, returns highest score
		Align(bool bBanded = false,					// true (currently experimental) to use banded or constrained SW 
				 UINT32 MaxStartNonOverlap = 20,	// if banded then initial non-overlapping path expected to be at most this many bp, 0 for no limits
				 double MaxPathLenDiff = 0.10);		// if banded then path length differential between query and probe for max score be at most this proportion of query length

	int											// returned total alignment length between probe and target including InDels
		GetAlignStats(UINT32 *pNumAlignedBases=NULL,// returned number of bases aligning between probe and target
				 UINT32 *pNumExactBases=NULL,          // of the aligning bases there were this many exact matches, remainder were substitutions
				 UINT32 *pNumProbeInsertBases=NULL,    // this many bases were inserted into the probe relative to the target
				 UINT32 *pNumTargInsertBases=NULL,		// this many bases were inserted into the target relative to the probe
				 UINT32 *pProbeStartOfs=NULL,			// alignment starts at this probe offset (1 based)
				 UINT32 *pTargStartOfs=NULL);			// alignment starts at this target offset (1 based)
	
	int GetProbeStartOfs(void);	// get offset (1..n) in probe at which alignment starts
	int GetTargStartOfs(void);  // get offset (1..n) in target at which alignment starts
	int GetNumAlignedBases(void);	    // get number of bases which align (exactly plus subs), excluding InDels 

	int GetProbeAlign( UINT32 Len, etSeqBase *pBuff); // get probe alignment
	int GetTargAlign( UINT32 Len, etSeqBase *pBuff);	 // get target alignment

	int						// < 0 errors, 0 - no anchors, 1 - anchors returned (note that retuned anchor offsets are 1 based inclusive)
		GetAnchors(UINT32 MinAnchorLen,	// anchors must be of at least this length
							   UINT32 *pProbeOfs5,	// returned probe 5' anchor starts at this probe sequence offset 
							   UINT32 *pTargOfs5,	// returned target 5' anchor starts at this target sequence offset 
							   UINT32 *pProbeOfs3,	// returned probe 3' anchor ends at this probe sequence offset
							   UINT32 *pTargOfs3);	// returned target 3' anchor ends at this target sequence offset

	int DumpScores(char *pszFile,		// dump Smith-Waterman matrix to this csv file
					char Down = '<',	// use this char to represent cell down link representing base inserted into target relative to probe
					char Left = '^',	// use this char to represent cell left link representing base inserted into probe relative to target
				    char Diag = '\\');	// use this char to represent cell diagonal  representing matching base either exact or mismatch

	};
