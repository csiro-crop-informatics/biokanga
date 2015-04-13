#pragma once
#include "./commdefs.h"

const UINT32 cNWMinProbeOrTargLen = 5;			// require probe and target lengths to be at least this number of bp
const UINT32 cNWMinCells = (cNWMinProbeOrTargLen * cNWMinProbeOrTargLen);	// not worth the effort if > this number of cells! minimum query and target length is 5bp
const UINT32 cNWMaxCells = 2000000000;		// 2 billion cells is getting big, big! restriction is that must be < 31bits, probelen * targlen limit
const UINT32 cNWMaxProbeOrTargLen = (cNWMaxCells/cNWMinProbeOrTargLen);  // require either probe or target length to be no longer than this limit

// default Needleman-Wunsch scores
const int cNWDfltMatchScore = 1;			// score for matching bases
const int cNWDfltMismatchScore = -1;		// mismatch penalty
const int cNWDfltGapOpenScore = -3;			// gap opening penalty
const int cNWDfltGapExtScore = -1;			// gap extension penalty

typedef unsigned int tNWTrcBckCell;		  // traceback cells are 32 bits
const UINT32 cNWScoreMsk     =   0x001fffff; // score is clamped to +/- 2^20-1 max with sign in cNWScoreNegFlg
const UINT32 cNWScoreNegFlg  =   0x00100000; // treat score as negative if this bit is set
const UINT32 cNWSpareBits=       0x07e00000; // these 6 bits currently not used
const UINT32 cNWGapOpnFlg  =     0x08000000; // set in cell if gap has been opened
const UINT32 cNWTrcBckMatchFlg = 0x10000000; // set if probe and target base were exactly matching
const UINT32 cNWTrcBckMsk  =     0xE0000000; // 3 bits to hold traceback direction (0=none,001=diag,010=up,100=left)
const UINT32 cNWTrcBckDiagFlg =  0x20000000; // traceback diagonally, bases align, could be exact match (cSWTrcBckMatchFlg set) or mismatched bases
const UINT32 cNWTrcBckDownFlg =  0x40000000; // traceback down, bases inserted into target (or deletion from probe)
const UINT32 cNWTrcBckLeftFlg =  0x80000000; // traceback left, bases inserted into probe (or deletion from target)


#pragma pack(1)
typedef struct TAG_sNWColBand {
	UINT32 TrcBckCellPsn;				// cells in this band have been allocated from cells starting at m_pTrcBckCells[TrcBckCellPsn-1], 0 if none allocated
	UINT32 StartTargBasePsn;			// cells in this probe column start at this target base position 1..m_TargLen
	UINT32 EndTargBasePsn;			    // cells in this probe column end at this target  1..m_TargLen 
} tsNWColBand;
#pragma pack()


class CNeedlemanWunsch
{
	bool m_bAligned;					// set true following successful Needleman-Wunsch scoring alignment
	bool m_bBanded;						// set true if banded processing

	tNWTrcBckCell *m_pTrcBckCells;		// array used for holding traceback cells
	size_t m_TrcBckCellsAllocdSize;		// allocation size
	UINT32 m_TrcBckCellsAllocd;			// actual cells in alloc'd m_pTrcBckCells (in tNWTrcBckCell's)
	UINT32 m_TrcBckCellsUsed;			// current cells used

	UINT32 m_ColBandsUsed;				// number of column active bands used - should be same as m_ProbeLen
	UINT32 m_ColBandsAllocd;			// number of column active bands allocated
	tsNWColBand *m_pColBands;			// each column, or probe base, in the matrix will contain a band of active target base cells

	etSeqBase *m_pProbe;				// alloc'd probe sequence memory
	UINT32 m_ProbeAllocd;				// actual m_pProbe alloc'd size (in teSWSeqBases's)
	UINT32 m_ProbeLen;					// current probe length

	etSeqBase *m_pTarg;					// alloc'd target sequence memory
	UINT32 m_TargAllocd;				// actual m_pTarg alloc'd size (in teSWSeqBases's)
	UINT32 m_TargLen;					// current targ length

	int m_MatchScore;					// score for matching bases
	int m_MismatchScore;				// mismatch penalty
	int m_GapOpenScore;					// gap opening penalty
	int m_GapExtScore;					// gap extension penalty

	UINT32 m_NWBandInitial;				// expecting initial overlap for max scoring alignment path to have non-overlap flank of at most this many bases
	double m_NWPathLenDiff;				// restrain path length differential between query and probe for max score be at most this proportion of query

	int m_PeakScore;					// peak score in any cell
	UINT32 m_NumBasesAligned;			// this many bases (exact and subs) were aligned between probe and target
	UINT32 m_NumBasesExact;				// of the m_NumBasesAligned, this many were exact matches, remainder were subs
	UINT32 m_NumProbeInserts;			// number of inserted bases into the probe relative to the target
	UINT32 m_NumTargInserts;			// number of inserted bases into the target relative to the probe

	UINT32 m_ProbeAlignStartOfs;		// offset in the probe at which the alignment starts
	UINT32 m_TargAlignStartOfs;			// offset in the target at which the alignment starts

	UINT32 m_PeakProbeIdx;				// highest scoring cell is at this matrix column or probe ofs
	UINT32 m_PeakTargIdx;				// highest scoring cell 

	tNWTrcBckCell *						// returns ptr to newly allocated cell  or NULL if errors
			AllocBandCell(UINT32 ProbeBasePsn,			// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn);		// current target base position 1..m_TargLen


	tNWTrcBckCell *						// returns ptr to cell 
			DerefBandCell(UINT32 ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn);		// current target base position 1..m_TargLen

	tNWTrcBckCell *					// returns ptr to cell  or NULL if errors
			DerefBandCellLeft(UINT32 ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn);		// current target base position 1..m_TargLen


	tNWTrcBckCell *					// returns ptr to cell  or NULL if errors
			DerefBandCellDiag(UINT32 ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn);		// current target base position 1..m_TargLen

	tNWTrcBckCell *									// returns ptr to cell  or NULL if errors
			DerefBandCellDown(UINT32 ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn);		// current target base position 1..m_TargLen


public:
	CNeedlemanWunsch(void);
	~CNeedlemanWunsch(void);

	void Reset(void);

    bool SetScores(int MatchScore = cNWDfltMatchScore,		// set scoring
					int MismatchScore = cNWDfltMismatchScore,
					int GapOpenScore = cNWDfltGapOpenScore,
					int GapExtScore = cNWDfltGapExtScore);

	bool SetProbe(UINT32 Len,etSeqBase *pSeq);
	bool SetTarg(UINT32 Len,etSeqBase *pSeq);

	int Align(void);			// Needleman-Wunsch style global alignment, returns highest score

	int											// returned total alignment length between probe and target including InDels
		GetAlignStats(UINT32 *pNumAlignedBases=NULL,// returned number of bases aligning between probe and target
				 UINT32 *pNumExactBases=NULL,          // of the aligning bases there were this many exact matches, remainder were substitutions
				 UINT32 *pNumProbeInsertBases=NULL,    // this many bases were inserted into the probe relative to the target
				 UINT32 *pNumTargInsertBases=NULL);	// this many bases were inserted into the target relative to the probe
	
	int GetNumAlignedBases(void);	    // get number of bases which align (exactly plus subs), excluding InDels 
	int GetProbeAlign(UINT32 Len, etSeqBase *pBuff); // get probe alignment
	int GetTargAlign(UINT32 Len, etSeqBase *pBuff);	 // get target alignment

	int DumpScores(char *pszFile,		// dump Needleman-Wunsch matrix to this csv file
					char Down = '<',	// use this char to represent cell down link representing base inserted into target relative to probe
					char Left = '^',	// use this char to represent cell left link representing base inserted into probe relative to target
				    char Diag = '\\');	// use this char to represent cell diagonal  representing matching base either exact or mismatch

    };


