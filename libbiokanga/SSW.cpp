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

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

#include "SSW.h"

CSSW::CSSW()
{
m_pAllocdCells = NULL;
m_pProbe = NULL;
m_pTarg = NULL;
Reset();
}


CSSW::~CSSW()
{
if(m_pAllocdCells != NULL)
	{
#ifdef _WIN32
	free(m_pAllocdCells);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocdCells != MAP_FAILED)
		munmap(m_pAllocdCells,m_AllocdCellSize);
#endif
	}
if(m_pProbe != NULL)
	delete m_pProbe;

if(m_pTarg != NULL)				// alloc'd target sequence memory
	delete m_pTarg;

}

void 
CSSW::Reset(void)					// reset state back to that immediately following instantiation
{
if(m_pAllocdCells != NULL)
	{
#ifdef _WIN32
	free(m_pAllocdCells);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocdCells != MAP_FAILED)
		munmap(m_pAllocdCells,m_AllocdCellSize);
#endif
	m_pAllocdCells = NULL;
	}
m_UsedCells = 0;
m_AllocdCells = 0;
m_AllocdCellSize = 0;


if(m_pProbe != NULL)
	{
	delete m_pProbe;
	m_pProbe = NULL;
	}

m_ProbeAllocd = 0;		// actual m_pProbe alloc'd size (in teSWSeqBases's)
m_ProbeLen = 0;			// current probe length

if(m_pTarg != NULL)				// alloc'd target sequence memory
	{
	delete m_pTarg;
	m_pTarg = NULL;
	}
m_TargAllocd = 0;			// actual m_pTarg alloc'd size (in teSWSeqBases's)
m_TargLen = 0;				// current targ length

m_MatchScore = cSSWDfltMatchScore;			
m_MismatchPenalty = cSSWDfltMismatchPenalty;		
m_GapOpenPenalty = cSSWDfltGapOpenPenalty;		
m_GapExtnPenalty = cSSWDfltGapExtnPenalty;			
m_DlyGapExtn = cSSWDfltGapOpenPenalty;				 
m_ProgPenaliseGapExtn = cSSWDfltProgPenaliseGapExtn;
m_AnchorLen = cSSWDfltAnchorLen;
m_MaxInitiatePathOfs = cMaxInitiatePathOfs;
m_MinNumExactMatches = cMinNumExactMatches;
m_MaxTopNPeakMatches = 0;
m_NumTopNPeakMatches = 0;
memset(m_TopPeakMatches,0,sizeof(m_TopPeakMatches));
}

bool 
CSSW::SetMaxInitiatePathOfs(int MaxInitiatePathOfs)	// require SW paths to have started within this many bp (0 to disable) on either the probe or target - effectively anchoring the SW 
{
if(MaxInitiatePathOfs < 0 || MaxInitiatePathOfs > 10000)
	return(false);
m_MaxInitiatePathOfs = MaxInitiatePathOfs;
return(true);
}


bool 
CSSW::SetMinNumExactMatches(int MinNumExactMatches)		// can process for at most this many peak matches in any probe vs target SW alignment
{
if(MinNumExactMatches < 10)								// must be a reasonable lower limit!
	return(false);
m_MinNumExactMatches = MinNumExactMatches;
return(true);
}

bool 
CSSW::SetTopNPeakMatches(int MaxTopNPeakMatches)		// can process for at most this many peak matches in any probe vs target SW alignment
{
if(MaxTopNPeakMatches < 0 || MaxTopNPeakMatches > cMaxTopNPeakMatches)		
	return(false);
m_MaxTopNPeakMatches = MaxTopNPeakMatches;
return(true);
}

int 
CSSW::GetTopNPeakMatches(tsSSWCell **pPeakMatches)		// returns number of peak matches in any probe vs target SW alignment
{
if(m_MaxTopNPeakMatches == 0 || m_NumTopNPeakMatches == 0)
	{
	if(pPeakMatches != NULL)
		*pPeakMatches = NULL;
	return(0);
	}
if(pPeakMatches != NULL)
	*pPeakMatches = m_TopPeakMatches;
return(m_NumTopNPeakMatches);
}

// SetScores
// Set match, mismatch, gap opening and gap extension scores
bool 
CSSW::SetScores(int MatchScore,			// score for match
				int MismatchPenalty,	// penalty for mismatch
				int GapOpenPenalty,		// penalty for opening a gap
				int GapExtnPenalty,		// penalty if extending already opened gap
				int DlyGapExtn,			// delayed gap penalties, only apply gap extension penalty if gap at least this length
				int ProgPenaliseGapExtn, // if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles
				int AnchorLen)			 // identified first and last anchors in alignment to be of at least this length
{
if(MatchScore <= 0 || MatchScore > 100 || 
	DlyGapExtn < 1 || DlyGapExtn > cSSWMaxPenaltyGap || 
	ProgPenaliseGapExtn < 0 || ProgPenaliseGapExtn > cSSWMaxPenaltyGap || 
	MismatchPenalty < -100 || MismatchPenalty > 0 || 
	GapOpenPenalty < -100 || GapOpenPenalty > 0 || 
	GapExtnPenalty < -100 ||  GapExtnPenalty > 0 ||
	AnchorLen < cSSWMinAnchorLen || AnchorLen > cSSWMaxAnchorLen)
	return(false);
m_MatchScore = MatchScore;				
m_MismatchPenalty = MismatchPenalty;		
m_GapOpenPenalty = GapOpenPenalty;			
m_GapExtnPenalty = GapExtnPenalty;			
m_DlyGapExtn = DlyGapExtn;				
if(ProgPenaliseGapExtn > 0 && ProgPenaliseGapExtn < DlyGapExtn)
	ProgPenaliseGapExtn = DlyGapExtn;	
m_ProgPenaliseGapExtn = ProgPenaliseGapExtn;	
m_AnchorLen = AnchorLen;
return(true);
}

// SetProbe
// Set probe sequence to use in subsequent alignments
bool 
CSSW::SetProbe(UINT32 Len,etSeqBase *pSeq)
{
if(Len < cSSWMinProbeOrTargLen || Len > cSSWMaxProbeOrTargLen || pSeq == NULL || *pSeq > eBaseN) 	// can't be bothered with very short or very long probes!
	return(false);

if(m_pProbe == NULL || m_ProbeAllocd < Len)
	{
	if(m_pProbe != NULL)
		delete m_pProbe;
	m_ProbeAllocd = Len + 100;
	m_pProbe = new etSeqBase [m_ProbeAllocd];
	if(m_pProbe == NULL)
		return(false);
	}
memmove(m_pProbe,pSeq,Len);
m_ProbeLen = Len;
return(true);
}

// SetTarg
// Set target sequence to use in subsequent alignments
bool 
CSSW::SetTarg( UINT32 Len,etSeqBase *pSeq)
{
if(Len < cSSWMinProbeOrTargLen || Len > cSSWMaxProbeOrTargLen || pSeq == NULL || *pSeq > eBaseN)	// can't be bothered with very short or very long targets!
	return(false);
if(m_pTarg == NULL || m_TargAllocd < Len)
	{
	if(m_pTarg != NULL)
		delete m_pTarg;
	m_TargAllocd = Len + 100;
	m_pTarg = new etSeqBase [m_TargAllocd];
	if(m_pTarg == NULL)
		return(false);
	}
memmove(m_pTarg,pSeq,Len);
m_TargLen = Len;
return(true);
}

bool 
CSSW::PreAllocMaxTargLen( UINT32 MaxTargLen)					// preallocate to process targets of this maximal length
{
if(m_pAllocdCells != NULL && (m_AllocdCells + 1000 >= MaxTargLen && m_AllocdCells <= MaxTargLen * 2))
	return(true);

if(m_pAllocdCells != NULL)
	{
#ifdef _WIN32
	free(m_pAllocdCells);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocdCells != MAP_FAILED)
		munmap(m_pAllocdCells,m_AllocdCellSize);
#endif
	m_pAllocdCells = NULL;
	}
m_AllocdCells = (UINT32)( ((UINT64)MaxTargLen * 120) / 100);		// a little extra safety margin and saves on reallocations
m_AllocdCellSize = sizeof(tsSSWCell) * m_AllocdCells;
#ifdef _WIN32
m_pAllocdCells = (tsSSWCell *) malloc(m_AllocdCellSize);
if(m_pAllocdCells == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for traceback cells",(INT64)m_AllocdCellSize);
	m_AllocdCellSize = 0;
	m_AllocdCells = 0;
	return(false);
	}
#else
// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
m_pAllocdCells = (tsSSWCell *)mmap(NULL,m_AllocdCellSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pAllocdCells == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for traceback cells",(INT64)m_AllocdCellSize);
	m_AllocdCellSize = 0;
	m_AllocdCells = 0;
	return(false);
	}
#endif
return(true);
}

tsSSWCell *								// smith-waterman style local alignment, returns highest accumulated exact matches cell
CSSW::Align(tsSSWCell *pPeakScoreCell)		// optionally also return conventional peak scoring cell
{
UINT32 IdxP;							// current index into m_Probe[]
UINT32 IdxT;							// current index into m_Targ[]
bool bMatch;
UINT32 StartIdxT;
UINT32 EndIdxT;

int DiagScore;							// putative diagonal score
int DiagPeakScore;
int LeftScore;							// putative left score
int LeftPeakScore;
int DownScore;							// putative down score
int DownPeakScore;
int MismatchPenalty;
int PrevScore;							// score in back referenced cell
etSeqBase *pProbe;
etSeqBase ProbeBase;
etSeqBase *pTarg;
etSeqBase TargBase;
int LeftInDelLen;
int DownInDelLen;

tsSSWCell *pCell;
tsSSWCell *pPrevCell;
tsSSWCell LeftCell;
tsSSWCell DiagCell;

if(m_ProbeLen < cSSWMinProbeOrTargLen || m_ProbeLen > cSSWMaxProbeOrTargLen ||  m_TargLen < cSSWMinProbeOrTargLen || m_TargLen > cSSWMaxProbeOrTargLen)
	return(NULL);	
memset(&m_PeakMatchesCell,0,sizeof(m_PeakMatchesCell));
memset(&m_PeakScoreCell,0,sizeof(m_PeakScoreCell));

if(m_AllocdCells + 10 < m_TargLen && !PreAllocMaxTargLen(m_TargLen))
	return(NULL);

if(m_MaxTopNPeakMatches)
	{
	m_NumTopNPeakMatches = 0;
	memset(m_TopPeakMatches,0,sizeof(m_TopPeakMatches));
	}

m_UsedCells = m_TargLen;

// cell defaults are score = 0, no gap extensions ...
memset(m_pAllocdCells,0,m_TargLen * sizeof(tsSSWCell));

pProbe = m_pProbe;
memset(&LeftCell,0,sizeof(tsSSWCell));
memset(&DiagCell,0,sizeof(tsSSWCell));
UINT32 NumCellsSkipped = 0;
UINT32 NumCellsChecked = 0;
UINT32 LastCheckedIdxT;
LastCheckedIdxT = m_TargLen;
for(IdxP = 0; IdxP < m_ProbeLen; IdxP++)
	{
	ProbeBase = *pProbe++ & ~cRptMskFlg;
	StartIdxT = 0;
	EndIdxT = min(m_TargLen,LastCheckedIdxT+2);
	pTarg = &m_pTarg[StartIdxT];
	pCell = m_pAllocdCells;

	for(IdxT = StartIdxT; IdxT < EndIdxT; IdxT++,pCell++)
		{
		// if m_MaxOverlapStartOfs > 0 then only starting new paths if within that max offset
		if(m_MaxInitiatePathOfs && IdxT >= (UINT32)m_MaxInitiatePathOfs && IdxP >= (UINT32)m_MaxInitiatePathOfs)
			{
			if(pCell->PeakScore == 0 && LeftCell.PeakScore == 0 && pCell[-1].PeakScore == 0)
				{
				DiagCell = LeftCell;
				LeftCell = *pCell;
				pTarg += 1;
				NumCellsSkipped += 1;
				continue;
				}
			NumCellsChecked += 1;
			}

		LastCheckedIdxT = IdxT;

		TargBase = *pTarg++ & ~cRptMskFlg;
		bMatch = ProbeBase == TargBase;
		LeftInDelLen = 0;
		DownInDelLen = 0;
		// calc the 3 scores (DiagScore, LeftScore, DownScore) so can determine which highest scoring path direction to take
		// DiagScore is either MatchScore or MismatchScore added to prev DiagScore score
		if(IdxT > 0 && IdxP > 0)
			{
			DiagCell = LeftCell;
			LeftCell = *pCell;
			PrevScore = DiagCell.CurScore;
			DiagPeakScore =  DiagCell.PeakScore;
			MismatchPenalty = m_MismatchPenalty;

			// if a mismatch and part of an exactly matching sequence, then look ahead and reduce penalty to same as for opening a gap if at least next three bases exactly match
			if(!bMatch && m_MismatchPenalty < m_GapOpenPenalty && (IdxP < m_ProbeLen-4) && (IdxT < m_TargLen - 4) && DiagCell.CurExactLen >= 3)
				{
				if((*pProbe & ~cRptMskFlg) == (*pTarg & ~cRptMskFlg) && 
					(pProbe[1] & ~cRptMskFlg) == (pTarg[1] & ~cRptMskFlg) &&
					(pProbe[2] & ~cRptMskFlg) == (pTarg[2] & ~cRptMskFlg) &&
					(pProbe[3] & ~cRptMskFlg) == (pTarg[3] & ~cRptMskFlg))
					MismatchPenalty = max(m_MismatchPenalty,m_GapOpenPenalty);
				}

			DiagScore = PrevScore + (bMatch ? m_MatchScore : MismatchPenalty);

			if(DiagScore > DiagPeakScore)
				DiagPeakScore = DiagScore;
			}
		else // else either IdxT or IdxP was zero
			{
			memset(pCell,0,sizeof(tsSSWCell));
			memset(&LeftCell,0,sizeof(tsSSWCell));
			memset(&DiagCell,0,sizeof(tsSSWCell));
			if(bMatch)
				{
				pCell->StartPOfs = IdxP+1;
				pCell->StartTOfs = IdxT+1;
				pCell->CurScore = m_MatchScore;
				pCell->PeakScore = m_MatchScore;
				pCell->NumMatches = pCell->NumExacts = pCell->CurExactLen = 1;
				}
			continue;
			}

		// leftscore is either GapExtnPenalty (if gap already opened and gap at least m_DlyGapExtn) or GapOpenScore added to prev left score
        // leftscore relative to diagscore and downscore determines if an insertion into probe is called
		int GapExtnPenalty;
		pPrevCell = &LeftCell;
		PrevScore = pPrevCell->CurScore;
		LeftPeakScore =  pPrevCell->PeakScore;
		GapExtnPenalty = 0;
		if(pPrevCell->LeftInDelLen)  // was there a gap previously opened?
			{
			LeftInDelLen =  1 + pPrevCell->LeftInDelLen; // extending existing gap
			if(m_GapExtnPenalty < 0 && LeftInDelLen >= m_DlyGapExtn)
				{
				if(m_ProgPenaliseGapExtn == 0 || LeftInDelLen < m_ProgPenaliseGapExtn)   
					GapExtnPenalty = m_GapExtnPenalty;	
				else
					{
					if(ProbeBase == (*pProbe & ~cRptMskFlg) && ProbeBase == (pProbe[1] & ~cRptMskFlg)) // homopolymer where k >= 3 runs are usually PacBio artefacts
						GapExtnPenalty = m_GapExtnPenalty;			// reduced penalty if homopolymer
					else
						GapExtnPenalty = m_GapExtnPenalty  * 2;
					}
				}
			LeftScore = PrevScore + GapExtnPenalty; 
			}
		else
			{
			LeftInDelLen = 1;			// opening a gap
			LeftScore = PrevScore + m_GapOpenPenalty;
			}
		if(LeftScore > LeftPeakScore)
			LeftPeakScore = LeftScore;

		// down score is either GapExtnPenalty (if gap already opened) or GapOpenScore added to prev down score
        // downscore relative to diagscore and leftscore determines if a deletion from probe is called
		pPrevCell = pCell - 1;
		PrevScore = pPrevCell->CurScore;
		DownPeakScore =  pPrevCell->PeakScore;
		GapExtnPenalty = 0;
		if(pPrevCell->DownInDelLen)  // was there a gap previously opened?
			{
			DownInDelLen =  1 + pPrevCell->DownInDelLen; // extending gap
			if(m_GapExtnPenalty < 0 && DownInDelLen >= m_DlyGapExtn)
				{
				if(m_ProgPenaliseGapExtn == 0 || DownInDelLen < m_ProgPenaliseGapExtn)   
					GapExtnPenalty = m_GapExtnPenalty;	
				else
					{
					if(TargBase == (*pTarg & ~cRptMskFlg) && TargBase == (pTarg[1] & ~cRptMskFlg)) // homopolymer where k >= 3 runs are usually PacBio artefacts
						GapExtnPenalty = m_GapExtnPenalty;				// reduced penalty if homopolymer
					else
						GapExtnPenalty = m_GapExtnPenalty * 2; 
					} 
				}
			DownScore = PrevScore + GapExtnPenalty; 
			}
		else
			{
			DownInDelLen = 1;
			DownScore = PrevScore + m_GapOpenPenalty;
			}
		if(DownScore > DownPeakScore)
			DownPeakScore = DownScore;
		
		// if no score was > 0 then reset cell, path can't be extended
		if(DiagScore <= 0 && DownScore <= 0 && LeftScore <= 0)
			{
#ifdef USETHISCODe
			if(pCell->PeakScore != 0)
#endif
				memset(pCell,0,sizeof(tsSSWCell));	
			continue;		
			}

		// select highest score into cell together with traceback and gap opened flag..
		if(DiagScore >= DownScore && DiagScore >= LeftScore) // if diag score at least equal highest then preference matches
			{
			*pCell = DiagCell;
			pCell->CurScore = DiagScore;
			pCell->PeakScore = DiagPeakScore;
			pCell->NumMatches += 1;
			pCell->DownInDelLen = 0;
			pCell->LeftInDelLen = 0;
			if(bMatch)
				{
				if(pCell->StartPOfs == 0)
					{
					pCell->StartPOfs = IdxP + 1;
					pCell->StartTOfs = IdxT + 1;
					}
				pCell->NumExacts+=1;
				pCell->CurExactLen += 1;
				if(pCell->CurExactLen >= m_AnchorLen)  // succession of exact matches qualify as anchor?
					{
					pCell->PLastAnchorEndOfs = IdxP + 1;
					pCell->TLastAnchorEndOfs = IdxT + 1;
					if(pCell->PFirstAnchorStartOfs == 0)
						{
						pCell->PFirstAnchorStartOfs = IdxP + 2 - m_AnchorLen;
						pCell->TFirstAnchorStartOfs = IdxT + 2 - m_AnchorLen;
						}
					}
				}
			else
				pCell->CurExactLen = 0;
			}
		else
			if(DownScore >= LeftScore) // down score at least as high as left score, note this means that insertions into the target are being preferenced
				{			
				*pCell = pCell[-1];
				pCell->CurScore = DownScore;
				pCell->PeakScore = DownPeakScore;
				pCell->DownInDelLen = DownInDelLen;
				pCell->LeftInDelLen = 0;

				pCell->NumBasesDel += 1;
				if(DownInDelLen == 1)
					pCell->NumGapsDel += 1;
				pCell->CurExactLen = 0;
				}
			else
				{			
				*pCell = LeftCell;							// left score the highest
				pCell->CurScore = LeftScore;
				pCell->PeakScore = LeftPeakScore;
				pCell->LeftInDelLen = LeftInDelLen;
				pCell->DownInDelLen = 0;
				pCell->NumBasesIns += 1;
				if(LeftInDelLen == 1)
					pCell->NumGapsIns += 1;
				pCell->CurExactLen = 0; 
				}
 
		if(pCell->NumExacts >= (UINT32)m_MinNumExactMatches)
			{
			if(pCell->PeakScore >= m_PeakMatchesCell.PeakScore && pCell->CurScore > ((m_PeakMatchesCell.PeakScore * 95) / 100))
				{
				m_PeakMatchesCell = *pCell;
				m_PeakMatchesCell.EndPOfs = IdxP + 1;
				m_PeakMatchesCell.EndTOfs = IdxT + 1;
				}

			if(m_MaxTopNPeakMatches > 0)
				{
				int TopNIdx;
				tsSSWCell *pCurMinCell;
				tsSSWCell *pTop = m_TopPeakMatches;
				if(m_NumTopNPeakMatches > 0)
					{
					pCurMinCell = NULL;
					for(TopNIdx = 0; TopNIdx < m_NumTopNPeakMatches; TopNIdx++,pTop++)
						{
						if(pCell->StartPOfs == pTop->StartPOfs && pCell->StartTOfs == pTop->StartTOfs) // updating an existing path?
							{
							if(pCell->NumExacts > pTop->NumExacts)
								{
								*pTop = *pCell;
								pTop->EndPOfs = IdxP + 1;
								pTop->EndTOfs = IdxT + 1;
								}
							break;
							}
						if(pCurMinCell == NULL || pCurMinCell->NumExacts < pTop->NumExacts)
							pCurMinCell = pTop;
						}
					if(TopNIdx == m_NumTopNPeakMatches)
						{
						if(m_NumTopNPeakMatches < m_MaxTopNPeakMatches)
							{
							*pTop = *pCell;
							pTop->EndPOfs = IdxP + 1;
							pTop->EndTOfs = IdxT + 1;
							m_NumTopNPeakMatches += 1;
							}
						else
							if(pCurMinCell != NULL && pCurMinCell->NumExacts < pCell->NumExacts)
								{
								*pCurMinCell = *pCell;
								pCurMinCell->EndPOfs = IdxP + 1;
								pCurMinCell->EndTOfs = IdxT + 1;
								}
						}
					}
				else
					{
					*pTop = *pCell;
					pTop->EndPOfs = IdxP + 1;
					pTop->EndTOfs = IdxT + 1;
					m_NumTopNPeakMatches = 1;
					}
				}
			}

#ifdef _PEAKSCOREACCEPT
		if(pPeakScoreCell != NULL && pCell->PeakScore >= m_PeakScoreCell.PeakScore)
			{
			m_PeakScoreCell = *pCell;
			m_PeakScoreCell.EndPOfs = IdxP + 1;
			m_PeakScoreCell.EndTOfs = IdxT + 1;
			}
#endif
		}
	}
#ifdef _PEAKSCOREACCEPT
if(pPeakScoreCell != NULL)
	*pPeakScoreCell = m_PeakScoreCell;
#endif
return(&m_PeakMatchesCell);
} 
