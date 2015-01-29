#pragma once
// Striped Smith-Waterman
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

#include "./commdefs.h"

const UINT32 cSSWMinProbeOrTargLen = 5;			// require probe and target lengths to be at least this number of bp
const UINT32 cSSWMaxProbeOrTargLen = 100000000;  // require either probe or target length to be no longer than this limit
const UINT32 cSSWMaxPenaltyGap = 1000000;		 // can specify gaps of upto this length before gap extension penalties are applied

const int cSSWDfltMatchScore = 1;		// score for matching bases
const int cSSWDfltMismatchPenalty = -1;	// mismatch penalty
const int cSSWDfltGapOpenPenalty = -3;	// gap opening penalty
const int cSSWDfltGapExtnPenalty = -1;	// gap extension penalty
const int cSSWDfltDlyGapExtn = 2;       // delay applying gap extension penalty until 2nd base in gap (must be range 1..cSSWMaxPenaltyGap)
const int cSSWDfltProgPenaliseGapExtn = 0; // if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles

const int cSSWMinAnchorLen = 3;			// can specify exactly matching anchors down to this length 
const int cSSWDfltAnchorLen = 10;		// default exactly matching anchors of at least this length to be identified 
const int cSSWMaxAnchorLen = 1000;      // can specify exactly matching anchors of at most this length  

const int cMaxInitiatePathOfs = 500;    // default is to require SW paths to have started within this many bp on either the probe or target - effectively anchoring the SW
const int cMinNumExactMatches = 100;    // default is only consider a path as being a peak path if that path contains at least this many exactly matching bases

const int cMaxTopNPeakMatches = 100;     // can process for at most this many peak matches in any probe vs target SW alignment

#pragma pack(1)
// each cell in stripe or column
typedef struct TAG_sSSWCell {
	UINT32 StartPOfs;			// path started at this probe sequence base offset + 1, 0 if no path
	UINT32 StartTOfs;			// path started at this target sequence base offset + 1, 0 if no path
	UINT32 EndPOfs;				// path ended at this probe sequence base offset + 1, 0 if no path
	UINT32 EndTOfs;				// path ended at this target sequence base offset + 1, 0 if no path

	UINT32 PFirstAnchorStartOfs; // first anchor of at least m_AnchorLen exactly matching bases started at this probe sequence base offset + 1, 0 if no anchor
	UINT32 TFirstAnchorStartOfs; // first anchor of at least m_AnchorLen exactly matching bases started at this target sequence base offset + 1, 0 if no anchor
	UINT32 PLastAnchorEndOfs;	 // last anchor of at least m_AnchorLen exactly matching bases ended at this probe sequence base offset + 1, 0 if no anchor
	UINT32 TLastAnchorEndOfs;	 // last anchor of at least m_AnchorLen exactly matching bases ended at this probe sequence base offset + 1, 0 if no anchor

	UINT32 NumMatches;			// total number of bases accepted along the path as matching, either exactly or as mismatches
	UINT32 NumExacts;			// total number of bases accepted along the path as exactly matching
	UINT32 NumGapsIns;			// total number of insertion gaps opened along the path
	UINT32 NumGapsDel;		    // total number of deletion gaps along the path
	UINT32 NumBasesIns;		    // total number of bases, over all gaps, inserted into probe
	UINT32 NumBasesDel;		    // total number of bases, over all gaps, deleted from probe

	INT32 PeakScore;			// peak score along path
	INT32 CurScore;				// current score
	UINT32 CurGapLen;			// current gap length
	UINT32 CurExactLen;			// current exactly matching sequence length
	
} tsSSWCell;


#pragma pack()

class CSSW
{
	etSeqBase *m_pProbe;				// alloc'd probe sequence memory
	UINT32 m_ProbeAllocd;				// actual m_pProbe alloc'd size (in teSWSeqBases's)
	UINT32 m_ProbeLen;					// current probe length

	etSeqBase *m_pTarg;					// alloc'd target sequence memory
	UINT32 m_TargAllocd;				// actual m_pTarg alloc'd size (in teSWSeqBases's)
	UINT32 m_TargLen;					// current targ length

	int m_MatchScore;					// score for matching bases (0..100)
	int m_MismatchPenalty;				// mismatch penalty (-100..0)
	int m_GapOpenPenalty;				// gap opening penalty (-100..0)
	int m_GapExtnPenalty;				// gap extension penalty (-100..0)
	int m_DlyGapExtn;					// delayed gap penalties, only apply gap extension penalty if gap at least this length
    int m_ProgPenaliseGapExtn;			// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles

	int m_MaxInitiatePathOfs;			// if non-zero then only allow new paths to start if within that offset (0 to disable) on either probe or target - effectively an anchored SW

	UINT32 m_AnchorLen;				// identified anchors between aligned probe and target must be at least this length
	UINT32 m_UsedCells;				// number of currently allocated cells used
	UINT32 m_AllocdCells;			// number of currently allocated cells
	size_t m_AllocdCellSize;        // total current allocation size for m_pAllocdCells 
	tsSSWCell *m_pAllocdCells;		// allocated to hold cells	

	int m_MinNumExactMatches;		// peak matches must contain at least this many exact matches to qualify as a peak match
	tsSSWCell m_PeakMatchesCell;	// cell identified as path containing highest number of matches 
	tsSSWCell m_PeakScoreCell;		// cell identified as usual SW peak scoring as conventional in most SW scoring schemes

	int m_MaxTopNPeakMatches;		// can identify at most this many top peak matches for current probe vs target
	int m_NumTopNPeakMatches;		// currently identified this many top peak matches for current probe vs target
	tsSSWCell m_TopPeakMatches[cMaxTopNPeakMatches];  // used to hold the tp[

public:
	CSSW();
	~CSSW();
	void Reset(void);			    // reset state back to that immediately following instantiation

	bool SetScores(int MatchScore= cSSWDfltMatchScore,			// score for match
				int MismatchPenalty  = cSSWDfltMismatchPenalty,	// penalty for mismatch
				int GapOpenPenalty  = cSSWDfltGapOpenPenalty,	// penalty for opening a gap
				int GapExtnPenalty  = cSSWDfltGapOpenPenalty,	// penalty if extending already opened gap
				int DlyGapExtn = cSSWDfltDlyGapExtn,			// delayed gap penalties, only apply gap extension penalty if gap at least this length
				int ProgPenaliseGapExtn = cSSWDfltProgPenaliseGapExtn,	// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles
				int AnchorLen = cSSWDfltAnchorLen);				// identified first and last anchors in alignment to be of at least this length

	bool SetMaxInitiatePathOfs(int MaxInitiatePathOfs = cMaxInitiatePathOfs);	// require SW paths to have started within this many bp (0 to disable) on either the probe or target - effectively anchoring the SW 

	bool SetMinNumExactMatches(int MinNumExactMatches = cMinNumExactMatches);		// require at least this many exactly matching in path to further process that path
	bool SetTopNPeakMatches(int MaxTopNPeakMatches = cMaxTopNPeakMatches/2);		// can process for at most this many peak matches in any probe vs target SW alignment

	int GetTopNPeakMatches(tsSSWCell **pPeakMatches = NULL);		// returns number of peak matches in any probe vs target SW alignment


	bool SetProbe( UINT32 Len,etSeqBase *pSeq);					// set probe sequence to use in subsequent alignments
	bool SetTarg( UINT32 Len,etSeqBase *pSeq);					// set target sequence to use in subsequent alignments

	bool PreAllocMaxTargLen( UINT32 MaxTargLen);				// preallocate to process targets of this maximal length

	tsSSWCell *								// smith-waterman style local alignment, returns highest accumulated exact matches scoring cell
				Align(tsSSWCell *pPeakScoreCell = NULL);		// optionally also return conventional peak scoring cell
};

