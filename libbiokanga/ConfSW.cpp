/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

CConfSW::CConfSW(void)
{
m_pTrcBckCells = NULL;
m_pProbe = NULL;
m_pTarg = NULL;
Reset();
}

CConfSW::~CConfSW(void)
{
Reset();
}

// Reset
// Resets instance to that immediately following that of instance construction
// Note that all buffers will be deallocated..
void
CConfSW::Reset(void)
{
if(m_pTrcBckCells != NULL)		// array used for holding traceback cells
	{
	delete m_pTrcBckCells;
	m_pTrcBckCells = NULL;
	}
m_TrcBckCellsAllocd = 0;	// actual cells in alloc'd m_pTrcBckCells (in tConfTrcBckCell's)
m_TrcBckCellsLen = 0;		// current cells used

if(m_pProbe != NULL)
	{
	delete m_pProbe;
	m_pProbe = NULL;
	}

m_ProbeAllocd = 0;	// actual m_pProbe alloc'd size (in teSWSeqBases's)
m_ProbeLen = 0;			// current probe length

if(m_pTarg != NULL)				// alloc'd target sequence memory
	{
	delete m_pTarg;
	m_pTarg = NULL;
	}
m_TargAllocd = 0;			// actual m_pTarg alloc'd size (in teSWSeqBases's)
m_TargLen = 0;				// current targ length

m_MatchScore = cSWDfltMatchScore;				// score for matching bases
m_MismatchPenalty = cSWDfltMismatchPenalty;		// mismatch penalty
m_GapOpenPenalty = cSWDfltGapOpenPenalty;			// gap opening penalty
m_GapExtnPenalty = cSWDfltGapExtnPenalty;			// gap extension penalty

m_PeakScore = 0;
m_PeakProbeIdx = 0;
m_PeakTargIdx = 0;
m_pPeakCell = NULL;

m_ProbeAlignStartIdx = 0;
m_TargAlignStartIdx = 0;

m_AlignLen = 0;
}

// SetScores
// Set match, mismatch, gap opening and gap extension scores
// Note that a check is made to ensure that MatchScore is > 0 and
// that all other scores are < 0
bool 
CConfSW::SetScores(teOctStructStats ConfParam, // which conformational parameter will be used in alignments
				   int MatchDiff,		// differential match threshold
				int MatchScore,			// score for match
				int MismatchPenalty,	// penalty for mismatch
				int GapOpenPenalty,	// penalty for opening a gap
				int GapExtnPenalty)	// penalty if extending already opened gap
{
if(MatchScore <= 0 || MismatchPenalty > 0 || MatchDiff < 0 || 
   GapOpenPenalty > 0 || GapExtnPenalty > 0)
	return(false);
m_ConfParam = ConfParam;
m_MatchDiff = MatchDiff;			// differential threshold at which matches are accepted
m_MatchScore = MatchScore;				// score for matching bases
m_MismatchPenalty = MismatchPenalty;		// mismatch penalty
m_GapOpenPenalty = GapOpenPenalty;			// gap opening penalty
m_GapExtnPenalty = GapExtnPenalty;			// gap extension penalty
m_AlignLen = 0;
m_pPeakCell = NULL;
return(true);
}

// SetProbe
// Set probe to use in subsequent alignments
// Probe is used as the row increment
// Only restriction is that the length must be >= 5
bool 
CConfSW::SetProbe(int ProbeLen,			   // probe length
				  int *pConf)			
{
if(ProbeLen < 5 || pConf == NULL) 	       // can't be bothered with very short probes!
	return(false);
if(m_pProbe == NULL || m_ProbeAllocd < ProbeLen)
	{
	if(m_pProbe != NULL)
		delete m_pProbe;
	m_pProbe = new int [ProbeLen + 1000];
	if(m_pProbe == NULL)
		return(false);
	m_ProbeAllocd = ProbeLen + 1000;
	}
int *pVal = m_pProbe;
for(int Idx=0;Idx < ProbeLen; Idx++,pVal++)
	*pVal = pConf[Idx];
m_ProbeLen = ProbeLen;
m_AlignLen = 0;
m_ProbeAlignStartIdx = 0;
m_TargAlignStartIdx = 0;
m_pPeakCell = NULL;
return(true);
}

// SetTarg
// Set target sequence to use in subsequent alignments
// Target sequence is used as the column increment
// Only restriction is that the length must be >= 5
bool 
CConfSW::SetTarg( int TargLen,			// target length starting at StartOfs
				  int *pConf)			
{
if(TargLen < 5 || pConf == NULL)		// can't be bothered with very short targets!
	return(false);
if(m_pTarg == NULL || m_TargAllocd < TargLen)
	{
	if(m_pTarg != NULL)
		delete m_pTarg;
	m_pTarg = new int [TargLen + 1000];
	if(m_pTarg == NULL)
		return(false);
	m_TargAllocd = TargLen + 1000;
	}
int *pVal = m_pTarg;
for(int Idx=0;Idx < TargLen; Idx++,pVal++)
	*pVal = pConf[Idx];
m_TargLen = TargLen;
m_AlignLen = 0;
m_ProbeAlignStartIdx = 0;
m_TargAlignStartIdx = 0;
m_pPeakCell = NULL;
return(true);
}

// Align
// Align probe SetProbe() against target SetTarg() using scores specified via SetScores() 
// using smith-waterman dynamic programing with optimisations
// Returns the peak score of all aligned subsequences
int										// peak score of all subsequence alignments
CConfSW::Align(void)
{
int IdxP;						// current index into m_Probe[]
int IdxT;						// current index into m_Targ[]
int NumCells;					// m_ProbeLen * m_TargLen - total number of cells
int NewScore;
int DiagScore;							// putative diagonal score
int LeftScore;							// putative left score
int UpScore;							// putative up score
int PrevScore;							// score in back referenced cell
unsigned int  DiagDir;					// to hold back refernce direction as diagonal  

int *pProbe;						// current &m_Probe[IdxP]
int ProbeBase;					// current m_Probe[IdxP]

int *pTarg;						// current &m_Targ[IdxT]
int TargBase;						// current m_Targ[IdxT]

tConfTrcBckCell *pCell;
tConfTrcBckCell *pPrevCell;
bool bMatch;

// realloc m_pTrcBckCells as may be required to hold all cells
NumCells = m_ProbeLen * m_TargLen;
if(m_pTrcBckCells == NULL || m_TrcBckCellsAllocd < NumCells)
	{
	if(m_pTrcBckCells != NULL)
		{
		delete m_pTrcBckCells;
		m_pTrcBckCells = NULL;
		}
	m_pTrcBckCells = new tConfTrcBckCell [NumCells + 10000]; // alloc more, may save a new alloc later..
	if(m_pTrcBckCells == NULL)
		return(0);
	m_TrcBckCellsAllocd = NumCells + 10000;
	}
m_TrcBckCellsLen = NumCells;

// cell defaults are score = 0, no backreferences
memset(m_pTrcBckCells,0,NumCells * sizeof(tConfTrcBckCell));

m_PeakScore = 0;
m_AlignLen = 0;
m_ProbeAlignStartIdx = 0;
m_TargAlignStartIdx = 0;
m_PeakProbeIdx = 0;
m_PeakTargIdx = 0;
m_pPeakCell = NULL;
pProbe = m_pProbe;
pCell = m_pTrcBckCells;
for(IdxP = 0; IdxP < m_ProbeLen; IdxP++)
	{
	ProbeBase = *pProbe++;
	pTarg = m_pTarg;
	for(IdxT = 0; IdxT < m_TargLen; IdxT++, pCell++)
		{
		TargBase = *pTarg++;
		if(ProbeBase > (TargBase - m_MatchDiff) &&
		   ProbeBase < (TargBase + m_MatchDiff))
			bMatch = true;
		else
			bMatch = false;

		// calc the 3 scores
		// diagonal is either MatchScore or MismatchScore added to prev diagonal score
		if(IdxT > 0 && IdxP > 0)
			{
			pPrevCell = pCell - m_TargLen - 1;
			PrevScore = (*pPrevCell & cConfScoreMsk);
			DiagScore = PrevScore + (bMatch ? m_MatchScore : m_MismatchPenalty);
			DiagDir = (DiagScore > 0) ? cConfTrcBckDiag : 0;
			if(DiagScore > (int)cConfScoreMsk)
				DiagScore = cConfScoreMsk;
			}
		else
			{
			DiagScore = bMatch ? m_MatchScore : m_MismatchPenalty;
			DiagDir = 0;
			}

		// leftscore is either GapExtScore (if gap already opened) or GapOpenScore added to prev left score
		if(IdxT > 0)
			{
			pPrevCell = pCell - 1;
			PrevScore = (*pPrevCell & cConfScoreMsk);
			LeftScore = PrevScore + ((*pPrevCell & cConfGapOpnFlg) ? m_GapExtnPenalty : m_GapOpenPenalty);
			}
		else
			LeftScore = m_GapOpenPenalty;

		// upscore is either GapExtnPenalty (if gap already opened) or GapOpenPenalty added to prev up score
		if(IdxP > 0)
			{
			pPrevCell = pCell - m_TargLen;
			PrevScore = (*pPrevCell & cConfScoreMsk);
			UpScore = PrevScore + ((*pPrevCell & cConfGapOpnFlg) ? m_GapExtnPenalty : m_GapOpenPenalty);
			}
		else
			UpScore = m_GapOpenPenalty;
		
		// if no score was > 0 then cell score stays as 0
		if(DiagScore <= 0 && UpScore <= 0 && LeftScore <= 0)
			{
			// *pCell = 0;	// array cells were initialised to be 0
			continue;		// leave cell as score = 0, direction = none
			}

		// select highest score into cell together with traceback and gap opened flag..
		if(DiagScore >= UpScore)
			{
			if(DiagScore >= LeftScore)
				{
				*pCell = DiagScore | DiagDir | (bMatch ? 0 : cConfTrcBckMatchFlg);
				NewScore = DiagScore;
				}
			else
				{			
				*pCell = LeftScore | cConfTrcBckLeft | (bMatch ? 0 : cConfGapOpnFlg);
				NewScore = LeftScore;
				}
			}
		else		// DiagScore < LeftScore
			if(UpScore >= LeftScore)
				{
				*pCell = UpScore | cConfTrcBckUp | (bMatch ? 0 : cConfGapOpnFlg);
				NewScore = UpScore;
				}
			else
				{			
				*pCell = LeftScore | cConfTrcBckLeft | (bMatch ? 0 : cConfGapOpnFlg);
				NewScore = LeftScore;
				}

		if(NewScore > (int)m_PeakScore)
			{
			m_PeakScore = NewScore;
			m_PeakProbeIdx = IdxP;
			m_PeakTargIdx = IdxT;
			m_pPeakCell = pCell;
			}
		}
	}
return(m_PeakScore);
}


// finds all the alignments above a specified alignment score
int
CConfSW::ScoreAlignments(int MinAlignScore) 
{
int NumAlignments;
int AlignScore = 0;
int IdxP;
int IdxT;
int TrcBckDir;
bool bInGap = false;
tConfTrcBckCell *pCurCell;
tConfTrcBckCell *pRefCell;
tConfTrcBckCell *pAlignStart;

pCurCell = &m_pTrcBckCells[m_TrcBckCellsLen-1];
NumAlignments = 0;
for(IdxT = m_TargLen -1 ; IdxT >= 0; IdxT--)
	{
	for(IdxP = m_ProbeLen - 1; IdxP >= 0; IdxP--, pCurCell--)
		{
			// skip if no score
		if(!(*pCurCell & cConfScoreMsk))
			continue;

		// skip if any other cell has a backlink to this cell
		if(IdxT < m_TargLen - 1)
			{
			pRefCell = pCurCell + 1;
			if((*pRefCell & cConfTrcBckMsk) == cConfTrcBckLeft)
				continue;
			}

		if(IdxP < m_ProbeLen - 1)
			{
			pRefCell = pCurCell + m_TargLen;
			if((*pRefCell & cConfTrcBckMsk) == cConfTrcBckUp)
				continue;
			}
		
		if(IdxP < m_ProbeLen - 1 && IdxT < m_TargLen - 1)
			{
			pRefCell = pCurCell + m_TargLen + 1;
			if((*pRefCell & cConfTrcBckMsk) == cConfTrcBckDiag)
				continue;
			}

		pRefCell = pCurCell;
		AlignScore = 0;
		bInGap = false;
		pAlignStart = NULL;
		do {
			switch(TrcBckDir = (*pRefCell & cConfTrcBckMsk)) {
				case cConfTrcBckDiag:		// back on the diagonal
					if(pAlignStart == NULL &&		// only start scoring at first known match
						*pRefCell & cConfTrcBckMatchFlg)
						pAlignStart = pRefCell;
					
					if(bInGap)			// if was in a gap then gap is now closed
						{
						AlignScore -= 2;		// gap penalty
						bInGap = false;
						}
					if(*pRefCell & cConfTrcBckMatchFlg) // if probe matched target
						AlignScore += 1;		// match score
					else
						AlignScore -= 1;		// mismatch score
		
					pRefCell -= m_TargLen + 1;
					break;

				case cConfTrcBckUp:				  // up .. must be in gap	
					if(pAlignStart!=NULL)
						bInGap = true;
					pRefCell -= m_TargLen;
					break;

				case cConfTrcBckLeft:		      // left .. must be in gap
					if(pAlignStart!=NULL)
						bInGap = true;
					pRefCell -= 1;
					break;

				default:
					break;
				}
			}
		while(TrcBckDir);
		if(bInGap)
			AlignScore -= 1;

		if(AlignScore >= MinAlignScore)
			{
			NumAlignments++;
			ProcessAlignment(AlignScore,pAlignStart);
			}
		}
	}
return(NumAlignments);
}



// GetAlignLength
// Returns the alignment length
// Also internally updates m_ProbeAlignStart and m_TargAlignStart
int
CConfSW::GetAlignLength(void)
{
unsigned int TrcBckDir;
tConfTrcBckCell *pPeakCell;

if(m_pPeakCell == NULL)
	return(-1);

if(m_AlignLen != 0)
	return(m_AlignLen);

m_ProbeAlignStartIdx = m_PeakProbeIdx;
m_TargAlignStartIdx = m_PeakTargIdx;

pPeakCell = m_pPeakCell;
do {
	m_AlignLen+=1;
	switch(TrcBckDir = (*pPeakCell & cConfTrcBckMsk)) {
		case cConfTrcBckDiag:		// back on the diagonal
			pPeakCell -= m_TargLen + 1;
			m_ProbeAlignStartIdx -= 1;
			m_TargAlignStartIdx -= 1;
			break;

		case cConfTrcBckUp:			// up
			pPeakCell -= m_TargLen; 
			m_ProbeAlignStartIdx -= 1;
			break;

		case cConfTrcBckLeft:		// left
			pPeakCell -= 1;
			m_TargAlignStartIdx -= 1;
			break;

		default:
			break;
		}
	}
while(TrcBckDir);
return(m_AlignLen);
}

// GetAlignLength
// Returns the alignment length
int
CConfSW::GetAlignLength(tConfTrcBckCell *pPeakCell)
{
unsigned int TrcBckDir;
unsigned int AlignLen;

if(pPeakCell == NULL)
	return(-1);
AlignLen = 0;
do {
	AlignLen+=1;
	switch(TrcBckDir = (*pPeakCell & cConfTrcBckMsk)) {
		case cConfTrcBckDiag:		// back on the diagonal
			pPeakCell -= m_TargLen + 1;
			break;

		case cConfTrcBckUp:			// up
			pPeakCell -= m_TargLen;
			break;

		case cConfTrcBckLeft:		// left
			pPeakCell -= 1;
			break;

		default:
			break;
		}
	}
while(TrcBckDir);
return(AlignLen);
}

int 
CConfSW::GetProbeAlignPsn(void) // get psn (0..n) in probe at which alignment starts
{
if(GetAlignLength() < 1)
	return(-1);
return(m_ProbeAlignStartIdx);
}

int 
CConfSW::GetTargAlignPsn(void)  // get psn (0..n) in target at which alignment starts
{
if(GetAlignLength() < 1)
	return(-1);
return(m_TargAlignStartIdx);
}

// get probe alignment
int 
CConfSW::GetProbeAlign(int InDelValue,		// value to use if InDel
					   int Len,				// max number of conf values to return 
					  int *pBuff)			// where to return ints
{
int BuffIdx;
unsigned int TrcBckDir;
tConfTrcBckCell *pPeakCell;
int *pProbe;

if(GetAlignLength() < 1)
	return(-1);
if(Len < m_AlignLen)
	return(-1);

pPeakCell = m_pPeakCell;
pProbe = m_pProbe + m_AlignLen + m_ProbeAlignStartIdx - 1;
BuffIdx = m_AlignLen - 1;

// highest score will always be at an exact match cell
TrcBckDir = cConfTrcBckDiag;
while(TrcBckDir  && BuffIdx >= 0)
	{
	switch(TrcBckDir) {
		case cConfTrcBckDiag:		// back on the diagonal
			pBuff[BuffIdx--] = *pProbe--; 
			pPeakCell -= m_TargLen + 1;
			break;

		case cConfTrcBckUp:			// up
			pBuff[BuffIdx--] = *pProbe--; 
			pPeakCell -= m_TargLen;
			break;

		case cConfTrcBckLeft:		// left
			pBuff[BuffIdx--] = eBaseInDel; 
			pPeakCell -= 1;
			break;
		}
	TrcBckDir = *pPeakCell--  & cConfTrcBckMsk;
	}
return(m_AlignLen);
}

// get target alignment
int 
CConfSW::GetTargAlign(int InDelValue,		// value to use if InDel
					 int Len,				// max number of conf values to return 
					  int *pBuff)			// where to return ints
{
int BuffIdx;
unsigned int TrcBckDir;
tConfTrcBckCell *pPeakCell;
int *pTarg;

if(GetAlignLength() < 1)
	return(-1);
if(Len < m_AlignLen)
	return(-1);

pPeakCell = m_pPeakCell;
pTarg = m_pTarg + m_AlignLen + m_ProbeAlignStartIdx - 1;
BuffIdx = m_AlignLen-1;

// highest score will always be at an exact match cell
TrcBckDir = cConfTrcBckDiag;
while(TrcBckDir && BuffIdx >= 0)
	{
	switch(TrcBckDir) {
		case cConfTrcBckDiag:		// back on the diagonal
			pBuff[BuffIdx--] = *pTarg--; 
			pPeakCell -= m_TargLen + 1;
			break;

		case cConfTrcBckUp:			// up
			pBuff[BuffIdx--] = *pTarg--; 
			pPeakCell -= m_TargLen;
			break;

		case cConfTrcBckLeft:		// left
			pBuff[BuffIdx--] = InDelValue; 
			pPeakCell -= 1;
			break;
		}
	TrcBckDir = *pPeakCell--  & cConfTrcBckMsk;
	}
return(m_AlignLen);
}

int 
CConfSW::ProcessAlignment(int AlignScore,tConfTrcBckCell *pCell)
{
int AlignLen = GetAlignLength(pCell);
return(AlignLen);
}



