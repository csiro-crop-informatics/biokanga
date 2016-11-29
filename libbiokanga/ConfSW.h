#pragma once
#include "./commdefs.h"

typedef unsigned int tConfTrcBckCell;			  // traceback cells are 32 bits, currently top 4 bits are unused
const unsigned int cConfScoreMsk       = 0x00ffffff;	// score is clamped to (2^24)-1 (16777215) max
const unsigned int cConfTrcBckMatchFlg = 0x01000000;    // set if match
const unsigned int cConfTrcBckMsk      = 0x06000000;	// 2 bits to hold traceback direction (0=none,1=diag,2=up,3=left)
const unsigned int cConfTrcBckDiag     = 0x02000000;	// traceback diagonally
const unsigned int cConfTrcBckUp       = 0x04000000;	// traceback up
const unsigned int cConfTrcBckLeft     = 0x06000000;	// traceback left
const unsigned int cConfGapOpnFlg      = 0x08000000;	// set in cell if gap has been opened


class CConfSW
{
	tConfTrcBckCell *m_pTrcBckCells;		// array used for holding traceback cells
	int m_TrcBckCellsAllocd;			// actual cells in alloc'd m_pTrcBckCells (in tConfTrcBckCell's)
	int m_TrcBckCellsLen;				// current cells used

	int *m_pProbe;						// alloc'd probe conformation memory
	int m_ProbeAllocd;					// actual m_pProbe alloc'd size (in int's)
	int m_ProbeLen;						// current probe length

	int *m_pTarg;						// alloc'd target conformation memory
	int m_TargAllocd;					// actual m_pTarg alloc'd size (in int's)
	int m_TargLen;						// current targ length

	teOctStructStats m_ConfParam;			// which conformational parameter will be used in alignments
    int m_MatchDiff;					// differential match threshold

	int m_MatchScore;					// score for matching bases
	int m_MismatchPenalty;				// mismatch penalty
	int m_GapOpenPenalty;				// gap opening penalty
	int m_GapExtnPenalty;				// gap extension penalty

	unsigned int m_PeakScore;
	int m_AlignLen;

	unsigned int m_ProbeAlignStartIdx;
	unsigned int m_TargAlignStartIdx;

	unsigned int m_PeakProbeIdx;
	unsigned int m_PeakTargIdx;
	tConfTrcBckCell *m_pPeakCell;

public:
	CConfSW(void);
	~CConfSW(void);

	void Reset(void);

    bool SetScores(teOctStructStats ConfParam, // which conformational parameter will be used in alignments
			    int MatchDiff,		// differential match threshold
				int MatchScore= cSWDfltMatchScore,			// score for match
				int MismatchPenalty  = cSWDfltMismatchPenalty,	// penalty for mismatch
				int GapOpenPenalty  = cSWDfltGapOpenPenalty,	// penalty for opening a gap
				int GapExtnPenalty  = cSWDfltGapOpenPenalty);	// penalty if extending already opened gap

	bool SetProbe(int ProbeLen,			// probe length
					  int *pConf);		// where to return ints
	
	bool SetTarg( int TargLen,			// target length
				  int *pConf);			// where to return ints

	int Align(void);			// smith-waterman style local alignment
	
	int GetProbeAlignPsn(void); // get psn (0..n) in probe at which alignment starts
	int GetTargAlignPsn(void);  // get psn (0..n) in target at which alignment starts
	int GetAlignLength(void);	 // get alignment length
	int GetAlignLength(tConfTrcBckCell *pPeakCell);

	// get probe alignment
	int GetProbeAlign(int InDelValue,		// value to use if InDel
					   int Len,	// max number of conf values to return
   					  int *pBuff);			// where to return ints

	int GetTargAlign(int InDelValue,		// value to use if InDel
					   int Len,	// max number of conf values to return
					  int *pBuff);			// where to return ints
	int ScoreAlignments(int MinAlignScore);
	int ProcessAlignment(int AlignScore,tConfTrcBckCell *pCell);
	
    };
