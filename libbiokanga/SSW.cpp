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
m_pAllocdTracebacks = NULL;
m_pMACols = NULL;  
m_pMAAlignOps = NULL;
m_pAllWinScores = NULL;
m_pParsimoniousBuff = NULL;
m_pszMAFAlignBuff = NULL;
m_pszConsensusBuff = NULL;
m_pConsConfSeq = NULL;

m_gzFile = NULL;
m_hMAFFile = -1;
m_hConsSeqFile = -1;
Reset();
}


CSSW::~CSSW()
{

if(m_gzFile != NULL)
	gzclose(m_gzFile);

if(m_hMAFFile != -1)
	close(m_hMAFFile);

if(m_hConsSeqFile != -1)
	close(m_hConsSeqFile);

if(m_pAllocdCells != NULL)
	{
#ifdef _WIN32
	free(m_pAllocdCells);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocdCells != MAP_FAILED)
		munmap(m_pAllocdCells,m_AllocdCellSize);
#endif
	}

if(m_pAllocdTracebacks != NULL)
	{
#ifdef _WIN32
	free(m_pAllocdTracebacks);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocdTracebacks != MAP_FAILED)
		munmap(m_pAllocdTracebacks,m_AllocdTracebacksSize);
#endif
	}

if(m_pMACols != NULL)
	{
#ifdef _WIN32
	free(m_pMACols);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pMACols != MAP_FAILED)
		munmap(m_pMACols,m_AllocMAColsSize);
#endif
	}

if(m_pMAAlignOps != NULL)
	{
#ifdef _WIN32
	free(m_pMAAlignOps);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pMAAlignOps != MAP_FAILED)
		munmap(m_pMAAlignOps,m_AllocdMAAlignOpsSize);
#endif
	}

if(m_pAllWinScores != NULL)
	delete m_pAllWinScores;

if(m_pProbe != NULL)
	delete m_pProbe;

if(m_pTarg != NULL)				// alloc'd target sequence memory
	delete m_pTarg;

if(m_pParsimoniousBuff != NULL)
	delete m_pParsimoniousBuff;

if(m_pszMAFAlignBuff != NULL)
	delete m_pszMAFAlignBuff;

if(m_pszConsensusBuff != NULL)
	delete m_pszConsensusBuff;

if(m_pConsConfSeq != NULL)
	delete m_pConsConfSeq;
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

if(m_pAllocdTracebacks != NULL)
	{
#ifdef _WIN32
	free(m_pAllocdTracebacks);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocdTracebacks != MAP_FAILED)
		munmap(m_pAllocdTracebacks,m_AllocdTracebacksSize);
#endif
	m_pAllocdTracebacks = NULL;
	}

if(m_pMACols != NULL)
	{
#ifdef _WIN32
	free(m_pMACols);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pMACols != MAP_FAILED)
		munmap(m_pMACols,m_AllocMAColsSize);
#endif
	m_pMACols = NULL;
	}

if(m_pMAAlignOps != NULL)
	{
#ifdef _WIN32
	free(m_pMAAlignOps);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pMAAlignOps != MAP_FAILED)
		munmap(m_pMAAlignOps,m_AllocdMAAlignOpsSize);
#endif
	m_pMAAlignOps = NULL;
	}

if(m_pszConsensusBuff != NULL)
	{
	delete m_pszConsensusBuff;
	m_pszConsensusBuff = NULL;
	}

if(m_pConsConfSeq != NULL)
	{
	delete m_pConsConfSeq;
	m_pConsConfSeq = NULL;
	}
if(m_pParsimoniousBuff != NULL)
	{
	delete m_pParsimoniousBuff;
	m_pParsimoniousBuff = NULL;
	}

if(m_pAllWinScores != NULL)
	{
	delete m_pAllWinScores;
	m_pAllWinScores = NULL;
	}

if(m_pProbe != NULL)
	{
	delete m_pProbe;
	m_pProbe = NULL;
	}

if(m_pTarg != NULL)				// alloc'd target sequence memory
	{
	delete m_pTarg;
	m_pTarg = NULL;
	}

if(m_pszMAFAlignBuff != NULL)
	{
	delete 	m_pszMAFAlignBuff;
	m_pszMAFAlignBuff = NULL;
	}

if(m_hConsSeqFile != -1)
	{
	close(m_hConsSeqFile);
	m_hConsSeqFile = -1;
	}

if(m_hMAFFile >= 0)
	{
	close(m_hMAFFile);
	m_hMAFFile = -1;
	}
if(m_gzFile != NULL)
	{
	gzclose(m_gzFile);
	m_gzFile = NULL;
	}
m_bIsGZ = false;
m_MAFFileOfs = 0;				
m_MAFAlignBuffIdx = 0;
m_MAFAlignBuffered = 0;		
m_MAFFileLineNum = 0;	 
m_AllocMAFAlignBuffSize = 0;	
m_AllocParsimoniousSize = 0;
m_UsedCells = 0;
m_AllocdCells = 0;
m_AllocdCellSize = 0;
m_UsedTracebacks = 0;
m_AllocdTracebacks = 0;
m_AllocdTracebacksSize = 0;
m_MAAlignOps = 0;
m_AllocdMAAlignOpsSize = 0;  
m_ProbeAllocd = 0;		
m_ProbeLen = 0;			
m_TargAllocd = 0;			
m_TargLen = 0;				
m_MAProbeSeqLen = 0;
m_MACols = 0;				
m_MADepth = 0;				
m_MAColSize = 0;
m_ProbeStartRelOfs = 0;			
m_TargStartRelOfs = 0;			
m_ProbeRelLen = 0;
m_TargRelLen = 0;
m_AllocMACols = 0;			
m_AllocMAColsSize = 0;	
m_MACoverage = 0;

m_ErrCorSeqID = 0;

m_ConfWinSize = cDfltConfWind;
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
m_bStartedMultiAlignments = false;
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
CSSW::PreAllocMaxTargLen( UINT32 MaxTargLen,		// preallocate to process targets of this maximal length
							bool bNoTracebacks)		// if true then don't preallocate for tracebacks			
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
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for cells",(INT64)m_AllocdCellSize);
	m_AllocdCellSize = 0;
	m_AllocdCells = 0;
	return(false);
	}
#else
// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
m_pAllocdCells = (tsSSWCell *)mmap(NULL,m_AllocdCellSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pAllocdCells == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for cells",(INT64)m_AllocdCellSize);
	m_pAllocdCells = NULL;
	m_AllocdCellSize = 0;
	m_AllocdCells = 0;
	return(false);
	}
#endif

// now prealloc for the tracebacks
if(m_pAllocdTracebacks != NULL)
	{
#ifdef _WIN32
	free(m_pAllocdTracebacks);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocdTracebacks != MAP_FAILED)
		munmap(m_pAllocdTracebacks,m_AllocdTracebacksSize);
#endif
	m_pAllocdTracebacks = NULL;
	}
m_AllocdTracebacks = 0;
m_AllocdTracebacksSize = 0;

if(!bNoTracebacks)
	{
	m_AllocdTracebacks = (UINT32)((UINT64)MaxTargLen * 500);		
	m_AllocdTracebacksSize = sizeof(tsSSWTraceback) * m_AllocdTracebacks;
#ifdef _WIN32
	m_pAllocdTracebacks = (tsSSWTraceback *) malloc(m_AllocdTracebacksSize);
	if(m_pAllocdTracebacks == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for traceback cells",(INT64)m_AllocdTracebacksSize);
		m_AllocdTracebacksSize = 0;
		m_AllocdTracebacks = 0;
		return(false);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pAllocdTracebacks = (tsSSWTraceback *)mmap(NULL,m_AllocdTracebacksSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pAllocdTracebacks == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for traceback cells",(INT64)m_AllocdTracebacksSize);
		m_pAllocdTracebacks = NULL;
		m_AllocdTracebacksSize = 0;
		m_AllocdTracebacks = 0;
		return(false);
		}
#endif
	}

// and lastly prealloc for the multialignment operators
if(m_pMAAlignOps != NULL)
	{
#ifdef _WIN32
	free(m_pMAAlignOps);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pMAAlignOps != MAP_FAILED)
		munmap(m_pMAAlignOps,m_AllocdMAAlignOpsSize);
#endif
	m_pMAAlignOps = NULL;
	}
m_pMAAlignOps = NULL;
m_AllocdMAAlignOpsSize = 0;

if(!bNoTracebacks)
	{		
	m_AllocdMAAlignOpsSize = MaxTargLen*3;
#ifdef _WIN32
	m_pMAAlignOps = (UINT8 *) malloc(m_AllocdMAAlignOpsSize);
	if(m_pMAAlignOps == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for multialignment operators",(INT64)m_AllocdMAAlignOpsSize);
		m_AllocdMAAlignOpsSize = 0;
		return(false);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pMAAlignOps = (UINT8 *)mmap(NULL,m_AllocdMAAlignOpsSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pMAAlignOps == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for multialignment operators",(INT64)m_AllocdMAAlignOpsSize);
		m_pMAAlignOps = NULL;
		m_AllocdMAAlignOpsSize = 0;
		return(false);
		}
#endif
	}
return(true);
}



int
CSSW::ValidateTracebacks(tsSSWCell *pCell)		// valdate that the path can be traced back for this cell
{
tsSSWTraceback *pCurTraceback;
tsSSWTraceback *pPrevTraceback;
UINT32 ExpIdxP = pCell->EndPOfs;
UINT32 ExpIdxT = pCell->EndTOfs;
int PathLen;

if((pCurTraceback = InitiateTraceback(pCell->EndPOfs,pCell->EndTOfs))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: InitiateTraceback for EndPOfs: %u EndTOfs: %u",pCell->EndPOfs,pCell->EndTOfs);
	return(0);
	}	

PathLen = 0;
do {
	if(ExpIdxP != (pCurTraceback->IdxP & cTrBkIdxMsk))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: ValidateTracebacks expected %u got %u",ExpIdxP,pCurTraceback->IdxP & cTrBkIdxMsk);
		return(0);
		}
	if(ExpIdxT != (pCurTraceback->IdxT & cTrBkIdxMsk))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: ValidateTracebacks expected %u got %u",ExpIdxT,pCurTraceback->IdxT & cTrBkIdxMsk);
		return(0);
		}

	switch(pCurTraceback->IdxP & cTrBkFlgsMsk) {
		case cTrBkFlgMatch:		// traceback direction: matching base, trace back to IdxT-1, IdxP-1
			ExpIdxP -= 1;
			ExpIdxT -= 1;
			break;
		case cTrBkFlgIns:		// traceback direction: base inserted into probe relative to target, trace back to IdxT, IdxP-1
			ExpIdxP -= 1;
			break;
		case cTrBkFlgDel:		// traceback direction: base deleted from probe relative to target, trace back to IdxT-1, IdxP
			ExpIdxT -= 1;
			break;
		case cTrBkFlgStart:		// traceback direction: 5' start of alignment, no further traceback
			break;
		}
	PathLen += 1;
	pPrevTraceback = pCurTraceback;
	}
while((pCurTraceback = NxtTraceback(pPrevTraceback))!=NULL);

if(pCell->StartPOfs != (pPrevTraceback->IdxP & cTrBkIdxMsk))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: ValidateTracebacks start path expected %u got %u",pCell->StartPOfs,pPrevTraceback->IdxP & cTrBkIdxMsk);
		return(0);
		}	
if(pCell->StartTOfs != (pPrevTraceback->IdxT & cTrBkIdxMsk))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: ValidateTracebacks start path expected %u got %u",pCell->StartTOfs,pPrevTraceback->IdxT & cTrBkIdxMsk);
		return(0);
		}	
return(PathLen);
}


tsSSWTraceback *
CSSW::InitiateTraceback(UINT32 IdxP, UINT32 IdxT)
{
UINT32 UsedTracebacks;
tsSSWTraceback *pCur;
if(IdxP == 0 || IdxT == 0 || (UsedTracebacks = m_UsedTracebacks) == 0 || m_AllocdTracebacks == 0 || m_pAllocdTracebacks == NULL)
	return(NULL);
pCur = &m_pAllocdTracebacks[UsedTracebacks-1];

while((pCur->IdxT & cTrBkIdxMsk) != IdxT || (pCur->IdxP & cTrBkIdxMsk) != IdxP)
	{
	if(--UsedTracebacks == 0 || (pCur->IdxP & cTrBkIdxMsk) < IdxP)
		return(NULL);
	pCur -= 1;
	}
return(pCur);
}


tsSSWTraceback *						 // next traceback along path; NULL if at 5' start of path
CSSW::NxtTraceback(tsSSWTraceback *pCur) // use to iterate over tracebacks
{
UINT32 IdxP;
UINT32 IdxT;

if(pCur == NULL || (pCur->IdxP & cTrBkIdxMsk) == 0 || (pCur->IdxT & cTrBkIdxMsk) == 0 || m_pAllocdTracebacks == NULL)
	return(NULL);

switch(pCur->IdxP & cTrBkFlgsMsk) {
	case cTrBkFlgMatch:		// traceback direction: matching base, trace back to IdxT-1, IdxP-1
		IdxP = (pCur->IdxP & cTrBkIdxMsk) - 1; 
		IdxT = (pCur->IdxT & cTrBkIdxMsk) - 1;
		break;
	case cTrBkFlgIns:		// traceback direction: base inserted into probe relative to target, trace back to IdxT, IdxP-1
		IdxP = (pCur->IdxP & cTrBkIdxMsk) - 1; 
		IdxT = pCur->IdxT & cTrBkIdxMsk;
		break;
	case cTrBkFlgDel:		// traceback direction: base deleted from probe relative to target, trace back to IdxT-1, IdxP
		IdxP = pCur->IdxP & cTrBkIdxMsk; 
		IdxT = (pCur->IdxT & cTrBkIdxMsk) - 1;
		break;
	case cTrBkFlgStart:		// traceback direction: 5' start of alignment, no further traceback
		return(NULL);
	}
if(IdxP == 0 || IdxT == 0)
	return(NULL);

while((pCur->IdxT & cTrBkIdxMsk) != IdxT || (pCur->IdxP & cTrBkIdxMsk) != IdxP)
	{
	if((pCur->IdxP & cTrBkIdxMsk) < IdxP)
		return(NULL);
	pCur -= 1;
	}
return(pCur);
}

UINT32											// number of tracebacks marked, can be less than actual path length if some tracebacks already marked
CSSW::MarkTracebackPath(UINT32 MarkFlag,		// mark the traceback path which ends at 3' IdxP and IdxT with this marker flag(s) into the tsSSWTraceback.IdxT
			UINT32 IdxP, UINT32 IdxT)
{
UINT32 NumMarked;
tsSSWTraceback *pCur;
MarkFlag &= cTrBkFlgsMsk;		
if((pCur = InitiateTraceback(IdxP,IdxT))==NULL)
	return(0);
NumMarked = 0;
do {
	if((pCur->IdxT & MarkFlag) == MarkFlag)		// if already marked then assume remainder of traceback path has been marked so no need to traceback further
		break;
	pCur->IdxT |= MarkFlag;
	NumMarked += 1;
	}
while((pCur=NxtTraceback(pCur))!=NULL);
return(NumMarked);
}


UINT32                                      // after reduction there are this many tracebacks retained
CSSW::ReduceTracebacks(UINT32 RetainFlag,	// reduce tracebacks by removing tracebacks which have NOT been marked with this flag in IdxT
				 UINT32 ResetFlags)		    // and reset these flags in the retained traceback in IdxT
{
UINT32 Idx;
UINT32 NumRetained;
tsSSWTraceback *pCur;
tsSSWTraceback *pNxt;

if(!m_UsedTracebacks || m_AllocdTracebacks == 0 || m_pAllocdTracebacks == NULL)
	return(0);

RetainFlag &= cTrBkFlgsMsk;
ResetFlags &= cTrBkFlgsMsk;

pCur = m_pAllocdTracebacks;
pNxt = pCur;
NumRetained = 0;
for(Idx = 0; Idx < m_UsedTracebacks; Idx++, pCur++)
	{
	if(pCur->IdxT & RetainFlag)
		{
		pCur->IdxT &= ~ResetFlags;
		*pNxt++ = *pCur;
		NumRetained += 1;
		}
	}
m_UsedTracebacks = NumRetained;
return(NumRetained);
} 

int													// NOTE: both probe and target sequences must have been set (SetProbe and SetTarg) before setting the alignment range
CSSW::SetAlignRange(UINT32 ProbeStartRelOfs,		// when aligning then start SW from this probe sequence relative offset
					  UINT32 TargStartRelOfs, 	// and SW starting from this target sequence relative offset
						UINT32 ProbeRelLen,		// and SW with this probe relative length starting from m_ProbeStartRelOfs - if 0 then until end of probe sequence
						UINT32 TargRelLen)		// and SW with this target relative length starting from m_TargStartRelOfs - if 0 then until end of target sequence
{
if(m_ProbeLen < cSSWMinProbeOrTargLen || m_ProbeLen > cSSWMaxProbeOrTargLen ||  m_TargLen < cSSWMinProbeOrTargLen || m_TargLen > cSSWMaxProbeOrTargLen)
	return(eBSFerrParams);	

if((ProbeStartRelOfs + cSSWMinProbeOrTargLen) > m_ProbeLen || (TargStartRelOfs + cSSWMinProbeOrTargLen) > m_TargLen)
	return(eBSFerrParams);	

if(ProbeRelLen == 0)
	ProbeRelLen = m_ProbeLen - ProbeStartRelOfs;	
if(TargRelLen == 0)
	TargRelLen = m_TargLen - ProbeStartRelOfs;

if((ProbeStartRelOfs + ProbeRelLen) > m_ProbeLen || (TargStartRelOfs + TargRelLen) > m_TargLen)
	return(eBSFerrParams);	

m_ProbeStartRelOfs = ProbeStartRelOfs;
m_TargStartRelOfs = TargStartRelOfs;
m_ProbeRelLen = ProbeRelLen;
m_TargRelLen = TargRelLen;
return(eBSFSuccess);
}

tsSSWCell *								// smith-waterman style local alignment, returns highest accumulated exact matches cell
CSSW::Align(tsSSWCell *pPeakScoreCell,	// optionally also return conventional peak scoring cell
				bool bNoTracebacks)		// if false then not using tracebacks
{
UINT32 IdxP;							// current index into m_Probe[]
UINT32 IdxT;							// current index into m_Targ[]
bool bMatch;
bool bMatchNxt3;
UINT32 StartIdxT;
UINT32 EndIdxT;

UINT32 TargRelLen;
UINT32 ProbeRelLen;

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

tsSSWTraceback *pTraceback;

if(m_ProbeLen < cSSWMinProbeOrTargLen || m_ProbeLen > cSSWMaxProbeOrTargLen ||  m_TargLen < cSSWMinProbeOrTargLen || m_TargLen > cSSWMaxProbeOrTargLen)
	return(NULL);	

if(m_ProbeRelLen == 0)
	ProbeRelLen = m_ProbeLen - m_ProbeStartRelOfs;	
else
	ProbeRelLen = m_ProbeRelLen;
if(m_TargRelLen == 0)
	TargRelLen = m_TargLen - m_ProbeStartRelOfs;
else
	TargRelLen = m_TargRelLen;

memset(&m_PeakMatchesCell,0,sizeof(m_PeakMatchesCell));
memset(&m_PeakScoreCell,0,sizeof(m_PeakScoreCell));

// allocating to hold full length even if relative length a lot shorter to reduce number of reallocations which may be subsequently required
if(((m_AllocdCells + 10 < TargRelLen) || (!bNoTracebacks && m_pAllocdTracebacks == NULL)) &&  !PreAllocMaxTargLen(m_TargLen,bNoTracebacks))
	return(NULL);

if(!bNoTracebacks && m_pAllocdTracebacks != NULL)
	memset(m_pAllocdTracebacks,0,sizeof(tsSSWTraceback));

if(m_MaxTopNPeakMatches)
	{
	m_NumTopNPeakMatches = 0;
	memset(m_TopPeakMatches,0,sizeof(m_TopPeakMatches));
	}

m_UsedCells = TargRelLen;

// cell defaults are score = 0, no gap extensions ...
memset(m_pAllocdCells,0,TargRelLen * sizeof(tsSSWCell));


memset(&LeftCell,0,sizeof(tsSSWCell));
memset(&DiagCell,0,sizeof(tsSSWCell));
UINT32 NumCellsSkipped = 0;
UINT32 NumCellsChecked = 0;
UINT32 LastCheckedIdxT;
LastCheckedIdxT = TargRelLen;
m_UsedTracebacks = 0;
pTraceback = m_pAllocdTracebacks;						// NOTE: NULL if tracebacks not required
pProbe = &m_pProbe[m_ProbeStartRelOfs];
for(IdxP = 0; IdxP < ProbeRelLen; IdxP++)
	{
	if(m_pAllocdTracebacks != NULL && (m_UsedTracebacks + 10) > (m_AllocdTracebacks - TargRelLen)) // ensure that sufficent memory has been allocated to hold any tracebacks in next sweep over the target
		{
		// try to reduce the number of tracebacks
		if(m_PeakMatchesCell.PeakScore > 0)  // peak scoring cell's path may have already terminated so mark that independently of those still in m_pAllocdCells[]
			MarkTracebackPath(cTrBkFlgRetain,m_PeakMatchesCell.EndPOfs,m_PeakMatchesCell.EndTOfs);
		pCell = m_pAllocdCells;
		for(IdxT = 0; IdxT < TargRelLen; IdxT++,pCell++)
			{
			if(pCell->PeakScore > 0)
				MarkTracebackPath(cTrBkFlgRetain,pCell->EndPOfs,pCell->EndTOfs);
			}
		ReduceTracebacks(cTrBkFlgRetain,cTrBkFlgRetain);
		if((m_UsedTracebacks + 10) > (m_AllocdTracebacks - TargRelLen))
			{
			UINT32 trbsreq;
			size_t memreq;
			void *pAllocd;
			trbsreq = (m_AllocdTracebacks + (TargRelLen * 500));
			memreq = (size_t)trbsreq * sizeof(tsSSWTraceback);

#ifdef _WIN32
			pAllocd = realloc(m_pAllocdTracebacks,memreq);
#else
			pAllocd = mremap(m_pAllocdTracebacks,m_AllocdTracebacksSize,memreq,MREMAP_MAYMOVE);
			if(pAllocd == MAP_FAILED)
				pAllocd = NULL;
#endif
			if(pAllocd == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Align: Memory re-allocation to %lld bytes - %s",memreq,strerror(errno));
				return(NULL);
				}

			m_pAllocdTracebacks = (tsSSWTraceback *)pAllocd;
			m_AllocdTracebacksSize = memreq;
			m_AllocdTracebacks = trbsreq; 
			}
		pTraceback = &m_pAllocdTracebacks[m_UsedTracebacks];
		}

	ProbeBase = *pProbe++ & ~cRptMskFlg;
	StartIdxT = 0;
	EndIdxT = min(TargRelLen,LastCheckedIdxT+2);
	pTarg = &m_pTarg[m_TargStartRelOfs];
	pCell = m_pAllocdCells;
	for(IdxT = 0; IdxT < EndIdxT; IdxT++,pCell++)
		{
		// if m_MaxOverlapStartOfs > 0 then only starting new paths if within that max offset
		// and early terminate low confidence paths; these are paths of longer than 1Kbp with coverage outside of 15% and with mean length of matching subseqs of less than 3.0bp
		if(m_MaxInitiatePathOfs && IdxT >= (UINT32)m_MaxInitiatePathOfs && IdxP >= (UINT32)m_MaxInitiatePathOfs)
			{
			if(pCell->PeakScore > 0 && (pCell->NumGapsIns + pCell->NumGapsDel) > 10) // only check for path termination if a score and there are some InDel events in path
				{
				int ProbeCovLen;    
				int TargCovLen;
				int CovRatio;
				int ScoreRatio;
				int MeanMatchLen;
				if((ProbeCovLen = (m_ProbeStartRelOfs + IdxP) - pCell->StartPOfs) > 1000 && (TargCovLen = (m_TargStartRelOfs + IdxT) - pCell->StartTOfs) > 1000)
					{
					ScoreRatio=(pCell->CurScore * 1000)/pCell->PeakScore; // expecting the score to not have dropped by more than 0.9 from the peak
					CovRatio = (ProbeCovLen*1000)/TargCovLen;   // expecting the probe and targ coverage lengths to be within 15% of each other
					MeanMatchLen = (pCell->NumMatches * 1000)/(pCell->NumGapsIns + pCell->NumGapsDel); // expecting the mean lengths of matching subsequences between probe and targ to be at least 3bp
					if(ScoreRatio < 900 || CovRatio < 850 || CovRatio > 1150 || MeanMatchLen < 3000)
						{
						memset(pCell,0,sizeof(tsSSWCell));
						memset(&LeftCell,0,sizeof(tsSSWCell));
						memset(&pCell[-1],0,sizeof(tsSSWCell));
						}
					}
				}
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
		bMatchNxt3 = (IdxP < ProbeRelLen-4) && (IdxT < TargRelLen - 4) && 
					((*pProbe & ~cRptMskFlg) == (*pTarg & ~cRptMskFlg) && 
					(pProbe[1] & ~cRptMskFlg) == (pTarg[1] & ~cRptMskFlg) &&
					(pProbe[2] & ~cRptMskFlg) == (pTarg[2] & ~cRptMskFlg));

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
			if(!bMatch && bMatchNxt3 && m_MismatchPenalty < m_GapOpenPenalty && DiagCell.CurExactLen >= 3)
				MismatchPenalty = max(m_MismatchPenalty,m_GapOpenPenalty);

			if(PrevScore || bMatchNxt3)
				DiagScore = PrevScore + (bMatch ? m_MatchScore : MismatchPenalty);
			else
				DiagScore = PrevScore + (bMatch ? 0 : MismatchPenalty);

			if(DiagScore > DiagPeakScore)
				DiagPeakScore = DiagScore;
			}
		else // else either IdxT or IdxP was zero
			{
			if(IdxP > 0)
				LeftCell = *pCell;
			else
				memset(&LeftCell,0,sizeof(tsSSWCell));
			memset(pCell,0,sizeof(tsSSWCell));
			memset(&DiagCell,0,sizeof(tsSSWCell));
			if(bMatch && bMatchNxt3)
				{
				pCell->StartPOfs = m_ProbeStartRelOfs + IdxP + 1;
				pCell->StartTOfs = m_TargStartRelOfs + IdxT + 1;
				pCell->EndPOfs = pCell->StartPOfs;
				pCell->EndTOfs = pCell->StartTOfs;
				pCell->CurScore = m_MatchScore;
				pCell->PeakScore = m_MatchScore;
				pCell->NumMatches = pCell->NumExacts = pCell->CurExactLen = 1;
				 
				if(pTraceback != NULL)
					{
					pTraceback->IdxP = pCell->StartPOfs | cTrBkFlgStart;
					pTraceback->IdxT = pCell->StartTOfs;
					pTraceback += 1;
					m_UsedTracebacks += 1;
					}
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
			memset(pCell,0,sizeof(tsSSWCell));	
			continue;		
			}

		// select highest score into cell together with traceback and gap opened flag..
		if(DiagScore >= DownScore && DiagScore >= LeftScore) // if diag score at least equal highest then preference matches
			{
			if(DiagCell.StartPOfs == 0 && !(bMatch && bMatchNxt3))           // if starting a new path then check if this path starts with at least 4 exacts
				{
				memset(pCell,0,sizeof(tsSSWCell));	
				continue;		
				}
			*pCell = DiagCell;
			pCell->CurScore = DiagScore;
			pCell->NumMatches += 1;
			pCell->DownInDelLen = 0;
			pCell->LeftInDelLen = 0;
			pCell->EndPOfs = m_ProbeStartRelOfs + IdxP + 1;
			pCell->EndTOfs = m_TargStartRelOfs + IdxT + 1;

			if(pTraceback != NULL)
				{
				pTraceback->IdxP = pCell->EndPOfs;
				pTraceback->IdxT = pCell->EndTOfs;
				if(!bMatch)
					pTraceback->IdxT |= cTrBkFlgSub; // flags that the match was not exact
				if(pCell->PeakScore > 0)
					pTraceback->IdxP |= cTrBkFlgMatch;     // if not starting path then diagonal traceback 
				pTraceback += 1;
				m_UsedTracebacks += 1;
				}
			pCell->PeakScore = DiagPeakScore;

			if(bMatch)
				{
				if(pCell->StartPOfs == 0)           // starting a new path?
					{
					pCell->StartPOfs = m_ProbeStartRelOfs + IdxP + 1;
					pCell->StartTOfs = m_TargStartRelOfs + IdxT + 1;
					}
				pCell->NumExacts+=1;
				pCell->CurExactLen += 1;
				if(pCell->CurExactLen >= m_AnchorLen)  // succession of exact matches qualify as anchor?
					{
					pCell->PLastAnchorEndOfs = m_ProbeStartRelOfs + IdxP + 1;
					pCell->TLastAnchorEndOfs = m_TargStartRelOfs + IdxT + 1;
					if(pCell->PFirstAnchorStartOfs == 0)
						{
						pCell->PFirstAnchorStartOfs = m_ProbeStartRelOfs + IdxP + 2 - m_AnchorLen;
						pCell->TFirstAnchorStartOfs = m_TargStartRelOfs + IdxT + 2 - m_AnchorLen;
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
				pCell->EndPOfs = m_ProbeStartRelOfs + IdxP + 1;
				pCell->EndTOfs = m_TargStartRelOfs + IdxT + 1;
				if(pTraceback != NULL)
					{
					pTraceback->IdxP = pCell->EndPOfs | cTrBkFlgDel;			// down traceback;
					pTraceback->IdxT = pCell->EndTOfs;
					pTraceback += 1;
					m_UsedTracebacks += 1;
					}
				pCell->NumBasesDel += 1;
				if(DownInDelLen == 1)
					pCell->NumGapsDel += 1;
				pCell->CurExactLen = 0;
				}
			else
				{			
				*pCell = LeftCell;												// left score the highest
				pCell->CurScore = LeftScore;
				pCell->PeakScore = LeftPeakScore;
				pCell->LeftInDelLen = LeftInDelLen;
				pCell->DownInDelLen = 0;
				pCell->EndPOfs = m_ProbeStartRelOfs + IdxP + 1;
				pCell->EndTOfs = m_TargStartRelOfs + IdxT + 1;
				if(pTraceback != NULL)
					{
					pTraceback->IdxP = pCell->EndPOfs | cTrBkFlgIns;            // left traceback;
					pTraceback->IdxT = pCell->EndTOfs;
					pTraceback += 1;
					m_UsedTracebacks += 1;
					}
				pCell->NumBasesIns += 1;
				if(LeftInDelLen == 1)
					pCell->NumGapsIns += 1;
				pCell->CurExactLen = 0; 
				}
 
		if(pCell->NumExacts >= (UINT32)m_MinNumExactMatches && pCell->CurExactLen >= 4) // only interested in putative paths which are terminating with at least 4 exact matches at terminal end - paths needed at least 4 exacts to start
			{
			if(pCell->PeakScore > m_PeakMatchesCell.PeakScore)
				m_PeakMatchesCell = *pCell;
			else
				{
				if(pCell->PeakScore == m_PeakMatchesCell.PeakScore)   // if this putative path has same peak score as previous best path then may be an extension of previous but has a localised lower score 
					{
					if((double)pCell->CurScore >= (m_PeakMatchesCell.PeakScore * 95.0)/100.0)	// allow current score to have dropped by at most ~ 95% of peak score
						{
						if((double)pCell->NumExacts >= (m_PeakMatchesCell.NumExacts * 101.0) / 100.0) // provided that the number of exactly matching bases has increased ~ 1.0%
							m_PeakMatchesCell = *pCell;
						}
					}
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
								pTop->EndPOfs = m_ProbeStartRelOfs + IdxP + 1;
								pTop->EndTOfs = m_TargStartRelOfs + IdxT + 1;
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
							pTop->EndPOfs =  m_ProbeStartRelOfs + IdxP + 1;
							pTop->EndTOfs = m_TargStartRelOfs + IdxT + 1;
							m_NumTopNPeakMatches += 1;
							}
						else
							if(pCurMinCell != NULL && pCurMinCell->NumExacts < pCell->NumExacts)
								{
								*pCurMinCell = *pCell;
								pCurMinCell->EndPOfs =  m_ProbeStartRelOfs + IdxP + 1;
								pCurMinCell->EndTOfs = m_TargStartRelOfs + IdxT + 1;
								}
						}
					}
				else
					{
					*pTop = *pCell;
					pTop->EndPOfs = m_ProbeStartRelOfs + IdxP + 1;
					pTop->EndTOfs = m_TargStartRelOfs + IdxT + 1;
					m_NumTopNPeakMatches = 1;
					}
				}
			}

#ifdef _PEAKSCOREACCEPT
		if(pPeakScoreCell != NULL && pCell->PeakScore >= m_PeakScoreCell.PeakScore)
			{
			m_PeakScoreCell = *pCell;
			m_PeakScoreCell.EndPOfs = m_ProbeStartRelOfs + IdxP + 1;
			m_PeakScoreCell.EndTOfs = m_TargStartRelOfs + IdxT + 1;
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


int    // total number of returned chars in pszBuffer for the textual representation of error corrected consensus sequence (could be multiple consensus sequences)
CSSW::MAlignCols2fasta(UINT32 ProbeID,	// identifies sequence which was used as the probe when determining the multialignments
					int MinConf,		// sequence bases averaged over a 50bp window must be of at least this confidence (0..9) with the initial and final bases having at least this confidence
				  int MinLen,			// and sequence lengths must be of at least this length 
				  UINT32 BuffSize,		// buffer allocated to hold at most this many chars
				  char *pszBuffer)		// output error corrected sequences to this buffer
{
int CurSeqLen;
int MaxCurSeqLen;
int LineLen;
int MeanConf;
int ConfWinTot;
UINT32 BuffOfs;
UINT32 StartBuffOfs;
UINT32 BelowMinConfBuffOfs;
char *pBuff;
int NewSeqStartOfs;
char szDescrLine[80];
int DescrLen;
etSeqBase Base;
char ChrBase;
int ConsIdx;
tsConfBase ConfWin[cMaxConfWindSize+1];

tsMAlignCol *pCol;

if(m_MACoverage < 1 || pszBuffer == NULL || BuffSize < m_MACols)
	return(-1);

ConsIdx = 0;
memset(ConfWin,0,sizeof(ConfWin));
CurSeqLen = 0;
MaxCurSeqLen = 0;
LineLen = 0;
BuffOfs = 0;
StartBuffOfs = 0;
BelowMinConfBuffOfs = 0;
pBuff = pszBuffer;
m_ErrCorSeqID = 0;

pCol = (tsMAlignCol *)m_pMACols;
while(pCol != NULL)
	{
	Base = pCol->ConsBase;
	if(Base <= eBaseN && !(CurSeqLen == 0 && MinConf > (int)pCol->ConsConf))
		{
		// base would be accepted into a sequence
		if(CurSeqLen == 0)									// potentially starting a new sequence?
			{
			ConsIdx = 0;
			CurSeqLen = 0;
			LineLen = 0;
			ConfWinTot = 0;
			BelowMinConfBuffOfs = 0;
			StartBuffOfs = BuffOfs;							// note where in buffer this potential sequence started in case sequence later needs to be retracted because it is not at least MinLen long
			NewSeqStartOfs = BuffOfs;
			pBuff = &pszBuffer[BuffOfs];
			}

		switch(Base) {                                      
			case eBaseA:
				ChrBase = 'a';
				break;
			case eBaseC:
				ChrBase = 'c';
				break;
			case eBaseG:
				ChrBase = 'g';
				break;
			case eBaseT:
				ChrBase = 't';
				break;
			case eBaseN:
				ChrBase = 'n';
				break;
			}

		if(LineLen >= 80)
			{
			*pBuff++ = '\n';
			BuffOfs += 1;
			LineLen = 0;
			}
		*pBuff++ = ChrBase;
		LineLen += 1;
		BuffOfs += 1;
		CurSeqLen += 1;
		if(CurSeqLen > MaxCurSeqLen)
			MaxCurSeqLen = CurSeqLen;

		if(MinConf > 0)
			{
			if(MinConf > (int)pCol->ConsConf)					// note offset of last observed low confidence in case mean confidence drops below minimum 
				BelowMinConfBuffOfs = BuffOfs - 1;              // if mean below minimum then truncate the sequence at this last observed

			if(CurSeqLen > m_ConfWinSize)
				ConfWinTot -= ConfWin[ConsIdx].Conf;
			ConfWin[ConsIdx].Base = Base;
			ConfWin[ConsIdx++].Conf = (int)pCol->ConsConf;
			if(ConsIdx == m_ConfWinSize)
				ConsIdx = 0;	
			ConfWinTot += pCol->ConsConf;

			if(CurSeqLen >= m_ConfWinSize)									// if sequence is at least m_ConfWin then can calculate the mean confidence over the previous 50bp
				{
				MeanConf = (ConfWinTot+(m_ConfWinSize/2))/m_ConfWinSize;
				if(BelowMinConfBuffOfs > 0 && MinConf > MeanConf)	// check if mean has dropped below the minimum required
					{
					CurSeqLen -= BuffOfs - BelowMinConfBuffOfs; 
					BuffOfs = BelowMinConfBuffOfs;	
					pBuff = &pszBuffer[BuffOfs];

					if(CurSeqLen < MinLen)							// if sequence is now not of an acceptable length then slough the sequence 
						{
						BuffOfs = StartBuffOfs;
						pBuff = &pszBuffer[BuffOfs];
						}
					else    // else accepting sequence
						{
						m_ErrCorSeqID += 1;
						if(m_ErrCorSeqID == 1)
							DescrLen = sprintf(szDescrLine,">ecseq%u_%d %d|%d\n",ProbeID,m_ErrCorSeqID,CurSeqLen,MinConf);
						else
							DescrLen = sprintf(szDescrLine,"\n>ecseq%u_%d %d|%d\n",ProbeID,m_ErrCorSeqID,CurSeqLen,MinConf);	
						pBuff = &pszBuffer[NewSeqStartOfs];
						memmove(pBuff + DescrLen,pBuff,BuffOfs - NewSeqStartOfs);	// make room for an inserted fasta descriptor line
						memcpy(pBuff,szDescrLine,DescrLen);
						BuffOfs += DescrLen;
						NewSeqStartOfs = BuffOfs;
						pBuff = &pszBuffer[BuffOfs];
						}
					CurSeqLen = 0;
					ConsIdx = 0;
					ConfWinTot = 0;
					BelowMinConfBuffOfs = 0;
					LineLen = 0;
					}
				}
			}
		}
	if(pCol->NxtColIdx == 0)
		break;
	pCol = (tsMAlignCol *)&m_pMACols[(pCol->NxtColIdx - 1) * (size_t)m_MAColSize];
	}

if(CurSeqLen > 0)
	{
	if(CurSeqLen < MinLen)	  // if sequence is not of an acceptable length then slough the sequence								
		BuffOfs = StartBuffOfs;
	else    // else accepting sequence
		{
		m_ErrCorSeqID += 1;
		if(m_ErrCorSeqID == 1)
			DescrLen = sprintf(szDescrLine,">ecseq%u_%d %d|%d\n",ProbeID,m_ErrCorSeqID,CurSeqLen,MinConf);
		else
			DescrLen = sprintf(szDescrLine,"\n>ecseq%u_%d %d|%d\n",ProbeID,m_ErrCorSeqID,CurSeqLen,MinConf);
		pBuff = &pszBuffer[NewSeqStartOfs];
		memmove(pBuff + DescrLen,pBuff,BuffOfs - NewSeqStartOfs);	// make room for an inserted fasta descriptor line
		memcpy(pBuff,szDescrLine,DescrLen);
		BuffOfs += DescrLen;
		}
	}
if(m_ErrCorSeqID)						// ensure that if sequences were accepted then the last sequence was terminated with a NL
	{
	if(pszBuffer[BuffOfs-1] != '\n')
		{
		pszBuffer[BuffOfs] = '\n';
		BuffOfs += 1;
		}
	}
return(BuffOfs);
}



int    // total number of returned chars in pszBuffer for the textual representation of error corrected consensus sequence (could be multiple consensus sequences)
CSSW::ConsConfSeq2Seqs(UINT32 ProbeID,	// identifies sequence which was used as the probe when determining the multialignments
					int MinConf,		// sequence bases averaged over a m_ConfWin window must be of at least this confidence (0..9) with the initial and final bases having at least this confidence
				  int MinLen)			// and sequence lengths must be of at least this length 
{
int CurSeqLen;
int MaxCurSeqLen;
int LineLen;
int MeanConf;
int ConfWinTot;
UINT32 BuffOfs;
UINT32 StartBuffOfs;
UINT32 BelowMinConfBuffOfs;
char *pBuff;
int NewSeqStartOfs;
char szDescrLine[80];
int DescrLen;
etSeqBase Base;
char ChrBase;
int ConsIdx;
tsConfBase ConfWin[cMaxConfWindSize+1];

tsConfBase *pCol;

ConsIdx = 0;
memset(ConfWin,0,sizeof(ConfWin));
CurSeqLen = 0;
MaxCurSeqLen = 0;
LineLen = 0;
BuffOfs = 0;
StartBuffOfs = 0;
BelowMinConfBuffOfs = 0;
pBuff = m_pszConsensusBuff;
m_ErrCorSeqID = 0;

pCol = (tsConfBase *)m_pConsConfSeq;
while((Base = pCol->Base) != eBaseEOS)
	{
	if(Base <= eBaseN && !(CurSeqLen == 0 && MinConf > (int)pCol->Conf))
		{
		// base would be accepted into a sequence
		if(CurSeqLen == 0)									// potentially starting a new sequence?
			{
			ConsIdx = 0;
			CurSeqLen = 0;
			LineLen = 0;
			ConfWinTot = 0;
			BelowMinConfBuffOfs = 0;
			StartBuffOfs = BuffOfs;							// note where in buffer this potential sequence started in case sequence later needs to be retracted because it is not at least MinLen long
			NewSeqStartOfs = BuffOfs;
			pBuff = &m_pszConsensusBuff[BuffOfs];
			}

		switch(Base) {                                      
			case eBaseA:
				ChrBase = 'a';
				break;
			case eBaseC:
				ChrBase = 'c';
				break;
			case eBaseG:
				ChrBase = 'g';
				break;
			case eBaseT:
				ChrBase = 't';
				break;
			case eBaseN:
				ChrBase = 'n';
				break;
			}

		if(LineLen >= 80)
			{
			*pBuff++ = '\n';
			BuffOfs += 1;
			LineLen = 0;
			}
		*pBuff++ = ChrBase;
		LineLen += 1;
		BuffOfs += 1;
		CurSeqLen += 1;
		if(CurSeqLen > MaxCurSeqLen)
			MaxCurSeqLen = CurSeqLen;

		if(MinConf > 0)
			{
			if(MinConf > (int)pCol->Conf)					// note offset of last observed low confidence in case mean confidence drops below minimum 
				BelowMinConfBuffOfs = BuffOfs - 1;              // if mean below minimum then truncate the sequence at this last observed

			if(CurSeqLen > m_ConfWinSize)
				ConfWinTot -= ConfWin[ConsIdx].Conf;
			ConfWin[ConsIdx].Base = Base;
			ConfWin[ConsIdx++].Conf = (int)pCol->Conf;
			if(ConsIdx == m_ConfWinSize)
				ConsIdx = 0;	
			ConfWinTot += (int)pCol->Conf;

			if(CurSeqLen >= m_ConfWinSize)									// if sequence is at least m_ConfWinSize long then can calculate the mean confidence over the previous m_ConfWinSize bp
				{
				MeanConf = (ConfWinTot+(m_ConfWinSize/2))/m_ConfWinSize;
				if(BelowMinConfBuffOfs > 0 && MinConf > MeanConf)	// check if mean has dropped below the minimum required
					{
					CurSeqLen -= BuffOfs - BelowMinConfBuffOfs; 
					BuffOfs = BelowMinConfBuffOfs;	
					pBuff = &m_pszConsensusBuff[BuffOfs];

					if(CurSeqLen < MinLen)							// if sequence is now not of an acceptable length then slough the sequence 
						{
						BuffOfs = StartBuffOfs;
						pBuff = &m_pszConsensusBuff[BuffOfs];
						}
					else    // else accepting sequence
						{
						m_ErrCorSeqID += 1;
						if(m_ErrCorSeqID == 1)
							DescrLen = sprintf(szDescrLine,">ecseq%u_%d %d|%d\n",ProbeID,m_ErrCorSeqID,CurSeqLen,MinConf);
						else
							DescrLen = sprintf(szDescrLine,"\n>ecseq%u_%d %d|%d\n",ProbeID,m_ErrCorSeqID,CurSeqLen,MinConf);	
						pBuff = &m_pszConsensusBuff[NewSeqStartOfs];
						memmove(pBuff + DescrLen,pBuff,BuffOfs - NewSeqStartOfs);	// make room for an inserted fasta descriptor line
						memcpy(pBuff,szDescrLine,DescrLen);
						BuffOfs += DescrLen;
						NewSeqStartOfs = BuffOfs;
						pBuff = &m_pszConsensusBuff[BuffOfs];
						}
					CurSeqLen = 0;
					ConsIdx = 0;
					ConfWinTot = 0;
					BelowMinConfBuffOfs = 0;
					LineLen = 0;
					}
				}
			}
		}
	pCol += 1;
	}

if(CurSeqLen > 0)
	{
	if(CurSeqLen < MinLen)	  // if sequence is not of an acceptable length then slough the sequence								
		BuffOfs = StartBuffOfs;
	else    // else accepting sequence
		{
		m_ErrCorSeqID += 1;
		if(m_ErrCorSeqID == 1)
			DescrLen = sprintf(szDescrLine,">ecseq%u_%d %d|%d\n",ProbeID,m_ErrCorSeqID,CurSeqLen,MinConf);
		else
			DescrLen = sprintf(szDescrLine,"\n>ecseq%u_%d %d|%d\n",ProbeID,m_ErrCorSeqID,CurSeqLen,MinConf);
		pBuff = &m_pszConsensusBuff[NewSeqStartOfs];
		memmove(pBuff + DescrLen,pBuff,BuffOfs - NewSeqStartOfs);	// make room for an inserted fasta descriptor line
		memcpy(pBuff,szDescrLine,DescrLen);
		BuffOfs += DescrLen;
		}
	}
if(m_ErrCorSeqID)						// ensure that if any sequences were accepted then the last sequence was terminated with a NL
	{
	if(m_pszConsensusBuff[BuffOfs-1] != '\n')
		{
		m_pszConsensusBuff[BuffOfs] = '\n';
		BuffOfs += 1;
		}
	m_pszConsensusBuff[BuffOfs] = '\0';
	}
else
	BuffOfs = 0;
return(BuffOfs);
}

// States
// 0 parsing starting from m_pszMAFAlignBuff[m_MAFAlignBuffIdx] until line starting with 'CBnnnnnn' where <nnnnnn> is the probe sequence identifier
// 1 parsing remainder of 'CBnnnnnn' line for consensus bases
// 2 parsing for line starting with 'CCnnnnnn' where <nnnnnn> is the probe sequence identifier
// 3 parsing remainder of 'CCnnnnnn' line for consensus confidence scores
// MAF assembly block processing is terminated on either errors or on the next occurance of  'CBnnnnnn'
//  
int    // parse out the next multialignment block consensus bases and consensus confidence scores into m_pConsConfSeq
CSSW::ParseConsConfSeq(bool bCpltdReadMAF,			// true if m_pConsConfSeq contains all remaining multialignment blocks loaded from file 
						 int MinErrCorrectLen,		// error corrected sequences must be at least this minimum length
						 UINT32 *pProbeID)			// returned probe identifier
{
int State;
int Cnt;
int Ofs;
int NumConBases;
int NumConConfs;
int ExpNumConConfs;
int SeqLen;
bool bSloughLine;
UINT32 ProbeID;
char Chr;
char *pBuff;
etSeqBase Base;
tsConfBase *pConsConfSeq;

State = 0;
NumConBases = 0;
NumConConfs = 0;
SeqLen = 0;
*pProbeID = 0;
bSloughLine = false;
pBuff = &m_pszMAFAlignBuff[m_MAFAlignBuffIdx];
pConsConfSeq = m_pConsConfSeq;
pConsConfSeq->Base = eBaseEOS;
pConsConfSeq->Conf = 0;
for(; m_MAFAlignBuffIdx < m_MAFAlignBuffered; m_MAFAlignBuffIdx += 1,pBuff+=1)
	{
	if(!bCpltdReadMAF && State == 0 && (m_MAFAlignBuffered - m_MAFAlignBuffIdx) < (cMaxMAFBlockErrCorLen * 2))
		{
		m_pConsConfSeq->Base = eBaseEOS;
		m_pConsConfSeq->Conf = 0;
		return(0);					// force load of more multialignments
		}

	if(bCpltdReadMAF && (m_MAFAlignBuffered - m_MAFAlignBuffIdx) < 100) // insufficent remaining for completion of an assembly block?
		{
		m_pConsConfSeq->Base = eBaseEOS;
		m_pConsConfSeq->Conf = 0;
		return(0);					// force load of more multialignments
		}
				
	Chr = toupper(*pBuff);
	if(Chr == ' ' || Chr == '\t') // slough any space or tabs
		continue;
	if(bSloughLine && Chr != '\n')	    // sloughing until start of next line?
		continue;
	if(State == 0 && Chr == '\n')
		{
		m_MAFFileLineNum += 1;
		bSloughLine = false;
		continue;
		}
	if(State == 0 && !(Chr == 'C' && pBuff[1] == 'B'))		// when expecting consensus lines to start then slough all lines not starting with 'CB', treating those sloughed as comment lines
		{
		bSloughLine = true;
		continue;
		}
	bSloughLine = false;

	switch(State) {
		case 0:						// expecting consensus line to start
			if((Cnt=sscanf(pBuff,"CB_%u%n",&ProbeID,&Ofs))==1)
				{
				if(ProbeID == 0)
					{
					m_pConsConfSeq->Base = eBaseEOS;
					m_pConsConfSeq->Conf = 0;
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadMAFBlock: Probe identifier is 0 at line %lld",m_MAFFileLineNum);
					return(eBSFerrAlignBlk);
					}
				State = 1;			// next expecting to parse out the consensus bases
				NumConBases = 0;
				NumConConfs = 0;	
				SeqLen = 0;
				pConsConfSeq = m_pConsConfSeq;
				pConsConfSeq->Base = eBaseEOS;
				pConsConfSeq->Conf = 0;
				m_MAFAlignBuffIdx += Ofs;
				pBuff += Ofs;
				continue;
				}
			bSloughLine = true;		// not a consensus sequence line so treating as a comment line or perhaps sequence lines from the previous multialignment block
			continue;

		case 1:						// parsing out consensus bases
			switch(Chr) {
				case '\r':			// treating as whitespace
					continue;
				case '\n':			// expecting consensus confidence scores line as next line
					State = 2;
					continue;
				case 'A': 
					Base = eBaseA;
					break;
				case 'C': 
					Base = eBaseC;
					break;
				case 'G': 
					Base = eBaseG;
					break;
				case 'T': case 'U':
					Base = eBaseT;
					break;
				case 'N': 
					Base = eBaseN;
					break;
				case '.':
					Base = eBaseUndef;
					break;
				default:
					m_pConsConfSeq->Base = eBaseEOS;
					m_pConsConfSeq->Conf = 0;
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadMAFBlock: Unexpected chars in consensus sequence at line %lld",m_MAFFileLineNum);
					return(eBSFerrAlignBlk);
				}
			if(NumConBases + 10 >= cMaxMAFBlockErrCorLen)
				{
				m_pConsConfSeq->Base = eBaseEOS;
				m_pConsConfSeq->Conf = 0;
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadMAFBlock: Excessively long consensus sequence at line %lld",m_MAFFileLineNum);
				return(eBSFerrAlignBlk);
				}
			NumConBases += 1;
			if(Base != eBaseUndef)
				SeqLen += 1;
			pConsConfSeq->Base = Base;
			pConsConfSeq->Conf = 0;
			pConsConfSeq += 1;
			pConsConfSeq->Base = eBaseEOS;
			pConsConfSeq->Conf = 0;
			continue;

		case 2:						// expecting consensus confidence scores line to start
			if(NumConBases < MinErrCorrectLen)  // if previous consensus sequence less than the minimum required then no point in continuing with the current multialignment
				{
				bSloughLine = true;
				State = 0;
				continue;
				}
			// earlier releases did not include the number of confidence scores
			if(Chr == 'C' && pBuff[1] == 'C' && pBuff[2] == ' ')
				{
				ExpNumConConfs = NumConBases;
				Ofs = 2;
				}
			else
				{
				if((Cnt=sscanf(pBuff,"CC_%u%n",&ExpNumConConfs,&Ofs))==1)
					{
					if(ExpNumConConfs != NumConBases)
						{
						m_pConsConfSeq->Base = eBaseEOS;
						m_pConsConfSeq->Conf = 0;
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadMAFBlock: Inconsistency in number of consensus confidence scores at line %lld",m_MAFFileLineNum);
						return(eBSFerrAlignBlk);
						}
					}
				else
					{
					m_pConsConfSeq->Base = eBaseEOS;
					m_pConsConfSeq->Conf = 0;
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadMAFBlock: Expected consensus confidence scores at line %lld",m_MAFFileLineNum);
					return(eBSFerrAlignBlk);
					}
				}

			NumConConfs = 0;
			State = 3;
			m_MAFAlignBuffIdx += Ofs;
			pBuff += Ofs;
			pConsConfSeq = m_pConsConfSeq;
			continue;

		case 3:					// parsing out consensus confidence scores
			if(Chr == '\r')		// treating as whitespace
				continue;
			if(Chr >= '0' && Chr <='9')
				{
				if(NumConConfs >= ExpNumConConfs)
					{
					m_pConsConfSeq->Base = eBaseEOS;
					m_pConsConfSeq->Conf = 0;
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadMAFBlock: Inconsistency in numbers of consensus confidence scores and consensus bases at line %lld",m_MAFFileLineNum);
					return(eBSFerrAlignBlk);
					}
				NumConConfs += 1;
				pConsConfSeq->Conf = (UINT8)(Chr - '0');	
				pConsConfSeq += 1;	
				continue;
				}
			if(Chr == '\n')
				{
				if(NumConConfs == ExpNumConConfs)	// must be same
					{
					m_pConsConfSeq[NumConConfs].Base = eBaseEOS;
					m_pConsConfSeq[NumConConfs].Conf = 0;
					*pProbeID = ProbeID;
					return(NumConConfs);
					}
				m_pConsConfSeq->Base = eBaseEOS;
				m_pConsConfSeq->Conf = 0;
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadMAFBlock: Inconsistency in numbers of consensus confidence scores and expected scores at line %lld",m_MAFFileLineNum);
				return(eBSFerrAlignBlk);
				}
			m_pConsConfSeq->Base = eBaseEOS;
			m_pConsConfSeq->Conf = 0;
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadMAFBlock: Unexpected char '%c' in line %lld",Chr,m_MAFFileLineNum);
			return(eBSFerrAlignBlk);
		}
	}

if(State == 3 &&  NumConConfs == ExpNumConConfs)
	{
	m_pConsConfSeq[NumConConfs].Base = eBaseEOS;
	m_pConsConfSeq[NumConConfs].Conf = 0;
	*pProbeID = ProbeID;
	return(NumConConfs);
	}

m_pConsConfSeq->Base = eBaseEOS;
m_pConsConfSeq->Conf = 0;
if(State != 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadMAFBlock: File truncated at line %lld",m_MAFFileLineNum);
	return(eBSFerrAlignBlk);
	}
return(0);
}

int
CSSW::GenConsensusFromMAF(int MinErrCorrectLen,		// error corrected sequences must be at least this minimum length
			 int MinConcScore,			// error corrected sequences trimmed until mean m_ConfWin concensus score is at least this threshold
			char *pszErrCorFile,		// name of file into which write error corrected sequences
			char *pszMultiAlignFile)	// name of file containing multiple alignments to process
{
int Rslt;
UINT32 CurProbeID;
int BuffCnt;
int BuffTopUp;
bool bCpltdReadMAF;
UINT32 NumParsedBlocks;
int ConsSeqLen;

m_AllocMAFAlignBuffSize = cMaxMAFBlockLen * 10;
if(m_pszMAFAlignBuff == NULL)
	{
	if((m_pszMAFAlignBuff = new char [m_AllocMAFAlignBuffSize])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for MAF block buffering");
		Reset();
		return(eBSFerrMem);
		}
	}
m_MAFAlignBuffIdx = 0;
m_MAFFileOfs = 0;
bCpltdReadMAF = false;

if(m_pszConsensusBuff == NULL)
	{
	if((m_pszConsensusBuff = new char [cMaxMAFBlockErrCorLen])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for consensus sequence buffering");
		Reset();
		return(eBSFerrMem);
		}
	}

if(m_pConsConfSeq == NULL)
	{
	if((m_pConsConfSeq = new tsConfBase [cMaxMAFBlockErrCorLen + 10])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for consensus confidence + bases buffering");
		Reset();
		return(eBSFerrMem);
		}
	}

m_MAFAlignBuffIdx = 0;
m_MAFFileOfs = 0;
m_MAFAlignBuffered = 0;

#ifdef _WIN32
m_hConsSeqFile = open(pszErrCorFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
if((m_hConsSeqFile = open(pszErrCorFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	if(ftruncate(m_hConsSeqFile,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszErrCorFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
#endif
if(m_hConsSeqFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszErrCorFile);
	Reset();
	return(eBSFerrCreateFile);
	}

// now open the pszMultiAlignFile
// if file has extension of ".gz' then assume that this file has been compressed and needs processing with gzopen/gzread/gzclose
int NameLen = (int)strlen(pszMultiAlignFile);
if(NameLen >= 4 && !stricmp(".gz",&pszMultiAlignFile[NameLen-3]))
	{
	if((m_gzFile = gzopen(pszMultiAlignFile,"r"))==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open %s as a gzip'd file - %s",pszMultiAlignFile,strerror(errno));
		Rslt = eBSFerrOpnFile;
		Reset();
		return(Rslt);
		}

	if(gzbuffer(m_gzFile,cgzAllocInBuffer)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to set gzbuffer size to %d",cgzAllocInBuffer);
		Rslt = eBSFerrMem;
		Reset();
		return(Rslt);
		}

	if(!gzdirect(m_gzFile))			// false if file has been compressed
		m_bIsGZ = true;
	else
		m_bIsGZ = false;
	}
else
	{
	m_bIsGZ = false;
#ifdef _WIN32
	m_hMAFFile = open(pszMultiAlignFile, O_READSEQ );		// file access is normally sequential..
#else
	m_hMAFFile = open64(pszMultiAlignFile, O_READSEQ );		// file access is normally sequential..
#endif
	if(m_hMAFFile == -1)							// check if file open succeeded
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open %s - %s",pszMultiAlignFile,strerror(errno));
		Rslt = eBSFerrOpnFile;
		Reset();
		return(eBSFerrOpnFile);
		}
	}

if(m_gzFile != NULL)
	{
	m_MAFFileOfs = gztell(m_gzFile);
	m_MAFAlignBuffered = gzread(m_gzFile, m_pszMAFAlignBuff, m_AllocMAFAlignBuffSize);
	}
else
	{
	m_MAFFileOfs = _lseeki64(m_hMAFFile,0,SEEK_CUR);
	m_MAFAlignBuffered = read(m_hMAFFile, m_pszMAFAlignBuff, m_AllocMAFAlignBuffSize);
	}

if(m_MAFAlignBuffered < 100)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open %s - %s",pszMultiAlignFile,strerror(errno));
	Rslt = eBSFerrOpnFile;
	Reset();
	return(eBSFerrOpnFile);
	}

NumParsedBlocks = 0;
do {
	if(!bCpltdReadMAF)
		{
		BuffTopUp = m_MAFAlignBuffered - m_MAFAlignBuffIdx;
		if(BuffTopUp < (m_AllocMAFAlignBuffSize / 2))
			{
			if(m_MAFAlignBuffIdx > 0 && BuffTopUp > 0)
				memmove(m_pszMAFAlignBuff,&m_pszMAFAlignBuff[m_MAFAlignBuffIdx],BuffTopUp);
			m_MAFAlignBuffered = BuffTopUp;
			m_MAFAlignBuffIdx = 0;
			if(m_gzFile != NULL)
				{
				m_MAFFileOfs = gztell(m_gzFile);
				BuffCnt = gzread(m_gzFile, &m_pszMAFAlignBuff[m_MAFAlignBuffered], m_AllocMAFAlignBuffSize - m_MAFAlignBuffered);
				}
			else
				{
				m_MAFFileOfs = _lseeki64(m_hMAFFile,0,SEEK_CUR);
				BuffCnt = read(m_hMAFFile, &m_pszMAFAlignBuff[m_MAFAlignBuffered], m_AllocMAFAlignBuffSize - m_MAFAlignBuffered);
				}
			if(BuffCnt <= 0)
				bCpltdReadMAF = true;
			else
				m_MAFAlignBuffered += BuffCnt;
			}
		}

	Rslt = ParseConsConfSeq(bCpltdReadMAF,MinErrCorrectLen,&CurProbeID);		// parses next complete multialignment block consensus confidence score and consensus bases into m_pConsConfSeq
	if(Rslt < 0 || (Rslt == 0 && bCpltdReadMAF))
		break;
	if(Rslt == 0)
		continue;

	// have consensus confidence scores and consensus bases - can now apply new MinConcScore threshold
	ConsSeqLen =			// total number of returned chars in m_pszConsensusBuff for the textual representation of error corrected consensus sequence (could be multiple consensus sequences)
		ConsConfSeq2Seqs(CurProbeID,				// identifies sequence which was used as the probe when determining the multialignments
					     MinConcScore,				// sequence bases averaged over a m_ConfWin window must be of at least this confidence (0..9) with the initial and final bases having at least this confidence
						 MinErrCorrectLen);			// and sequence lengths must be of at least this length 

	if(ConsSeqLen > 0)
		CUtility::SafeWrite(m_hConsSeqFile,m_pszConsensusBuff,ConsSeqLen);
	NumParsedBlocks += 1;
	}
while(Rslt >= 0);

if(m_hConsSeqFile != -1)
	{
#ifdef _WIN32
	_commit(m_hConsSeqFile);
#else
	fsync(m_hConsSeqFile);
#endif
	close(m_hConsSeqFile);
	m_hConsSeqFile = -1;
	}
if(m_hMAFFile >= 0)
	{
	close(m_hMAFFile);
	m_hMAFFile = -1;
	}
if(m_gzFile != NULL)
	{
	gzclose(m_gzFile);
	m_gzFile = NULL;
	}

if(NumParsedBlocks == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to parse any multiple alignment blocks from file '%s', is this a multiple alignment file?",pszMultiAlignFile);
	Rslt = eBSFerrAlignBlk;
	}
return(Rslt);
}


int      // total number of returned chars in pszBuffer for the textual representation of the multialignment 
CSSW::MAlignCols2MFA(UINT32 ProbeID,		// identifies sequence which was used as the probe when determining the multialignments
					UINT32 BuffSize,		// buffer allocated to hold at most this many chars
					char *pszBuffer)		// output multialignment textual representation to this buffer
{
UINT32 BuffOfs;
char *pBuff;
UINT32 DepthIdx;
etSeqBase Base;
int ConsConf;
char ChrBase;

tsMAlignCol *pCol;

if(m_MACoverage < 2 || pszBuffer == NULL || BuffSize < (m_MACoverage * (m_MACols + 20)))
	return(-1);

BuffOfs = 0;
pBuff = pszBuffer;
BuffOfs += sprintf(pBuff,"CB_%07u ",ProbeID);
pBuff = &pszBuffer[BuffOfs];
pCol = (tsMAlignCol *)m_pMACols;

while(pCol != NULL)
	{
	Base = pCol->ConsBase;
	switch(Base) {
		case eBaseA:
			ChrBase = 'a';
			break;
		case eBaseC:
			ChrBase = 'c';
			break;
		case eBaseG:
			ChrBase = 'g';
			break;
		case eBaseT:
			ChrBase = 't';
			break;
		case eBaseN:
			ChrBase = 'n';
			break;
		case eBaseUndef:
			ChrBase = '.';
			break;
		case eBaseInDel:
			ChrBase = '-';
			break;
		default:
			ChrBase = '?';
			break;
		}
	*pBuff++ = ChrBase;
	BuffOfs += 1;
	if(pCol->NxtColIdx == 0)
		break;
	pCol = (tsMAlignCol *)&m_pMACols[(pCol->NxtColIdx - 1) * (size_t)m_MAColSize];
	}
BuffOfs += sprintf(pBuff,"\n");
pBuff = &pszBuffer[BuffOfs];

BuffOfs += sprintf(pBuff,"CC_%07u ",m_MACols);
pBuff = &pszBuffer[BuffOfs];
pCol = (tsMAlignCol *)m_pMACols;
while(pCol != NULL)
	{
	ConsConf = (int)pCol->ConsConf;
	if(ConsConf > 9)
		ConsConf = 9;
	BuffOfs += sprintf(pBuff,"%d",ConsConf);
	pBuff = &pszBuffer[BuffOfs];
	if(pCol->NxtColIdx == 0)
		break;
	pCol = (tsMAlignCol *)&m_pMACols[(pCol->NxtColIdx - 1) * (size_t)m_MAColSize];
	}
BuffOfs += sprintf(pBuff,"\n");
pBuff = &pszBuffer[BuffOfs];

for(DepthIdx = 0; DepthIdx < m_MACoverage; DepthIdx++)
	{
	pBuff = &pszBuffer[BuffOfs];
	if((BuffSize - BuffOfs) < m_MACols + 100)
		{
		BuffOfs += sprintf(pBuff,"Truncated, insufficent buffering for any additional rows\n");
		pBuff = &pszBuffer[BuffOfs];
		*pBuff = '\0';
		return(BuffOfs);
		}

	BuffOfs += sprintf(pBuff,"CS_%03d     ",DepthIdx+1);
	pBuff = &pszBuffer[BuffOfs];
	pCol = (tsMAlignCol *)m_pMACols;
	while(pCol != NULL)
		{
		Base = pCol->Bases[DepthIdx];
		switch(Base) {
			case eBaseA:
				ChrBase = 'A';
				break;
			case eBaseC:
				ChrBase = 'C';
				break;
			case eBaseG:
				ChrBase = 'G';
				break;
			case eBaseT:
				ChrBase = 'T';
				break;
			case eBaseN:
				ChrBase = 'N';
				break;
			case eBaseUndef:
				ChrBase = '.';
				break;
			case eBaseInDel:
				ChrBase = '-';
				break;
			default:
				ChrBase = '?';
				break;
			}
		*pBuff++ = ChrBase;
		BuffOfs += 1;
		if(pCol->NxtColIdx == 0)
			break;
		pCol = (tsMAlignCol *)&m_pMACols[(pCol->NxtColIdx - 1) * (size_t)m_MAColSize];
		}
	BuffOfs += sprintf(pBuff,"\n");
	pBuff = &pszBuffer[BuffOfs];
	}
return(BuffOfs);
}

int
CSSW::StartMultiAlignments(int SeqLen,		// probe sequence is this length
					etSeqBase *pProbeSeq,	// probe sequence 
					int Alignments)			// number of pairwise alignments to allocate for
{
int Idx;
UINT8 *pBase;
tsMAlignCol *pCol;
size_t memreq;
m_MACols = 0;
m_MADepth = 0;
m_MACoverage = 0;

if(Alignments < 2 || Alignments > 200 || SeqLen < 50 || pProbeSeq == NULL)
	return(0);

m_MAProbeSeqLen = SeqLen;
m_MAColSize = (sizeof(tsMAlignCol) + Alignments);
memreq = (size_t)m_MAColSize * SeqLen * 2; // allocate to hold at least SeqLen columns with Depth bases and a 2x overallocation to reduce chances of a reallocation later on as sequence insertions are discovered
if(m_pMACols == NULL || memreq > m_AllocMAColsSize)
	{
	if(m_pMACols != NULL)
		{
#ifdef _WIN32
		free(m_pMACols);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(m_pMACols != MAP_FAILED)
			munmap(m_pMACols,m_AllocMAColsSize);
#endif		
		m_pMACols = NULL;
		m_AllocMAColsSize = 0;
		}
	
#ifdef _WIN32
	m_pMACols = (UINT8 *) malloc(memreq);
	if(m_pMACols == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for multialignment",(INT64)memreq);
		return(eBSFerrMem);
		}
#else
		// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pMACols = (UINT8 *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pMACols == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for multialignment",(INT64)memreq);
		m_pMACols = NULL;
		return(eBSFerrMem);
		}
#endif
	m_AllocMAColsSize = memreq;	
	}

m_AllocMACols = (UINT32)(m_AllocMAColsSize/m_MAColSize); 
memset(m_pMACols,0,m_AllocMAColsSize);
m_MACols = SeqLen;
m_MADepth = Alignments+1;

// initialise with the reference sequence bases and column depth as 1
pBase = pProbeSeq; 
pCol = (tsMAlignCol *)m_pMACols;
for(Idx = 1; Idx <= SeqLen; Idx++, pBase++)
	{
	pCol->Ofs = Idx;
	pCol->ColIdx = Idx;
	pCol->PrevColIdx = Idx-1;
	pCol->NxtColIdx = Idx == SeqLen ? 0 : Idx + 1;
	pCol->Depth = 1;
	pCol->Extn = 0;
	pCol->ConsBase = *pBase;
	pCol->Bases[0] = *pBase;
	pCol = (tsMAlignCol *) ((UINT8 *)pCol + m_MAColSize);
	}
m_MACoverage = 1;
m_bStartedMultiAlignments = true;
return(eBSFSuccess);
}

int												// attempting to determine if path is artfact resulting from aligning to a paralogous fragment
CSSW::ClassifyPath(int MaxArtefactDev,			// classify path as artefactual if sliding window (currently 1Kbp) over any overlap deviates by more than this percentage from the overlap mean
				    UINT32 ProbeStartOfs,		// alignment starts at this probe sequence offset (1..n)
					UINT32 ProbeEndOfs,			// alignment ends at this probe sequence offset
					UINT32 TargStartOfs,		// alignment starts at this target sequence offset (1..n)
					UINT32 TargEndOfs)			// alignment ends at this target sequence offset
{
UINT32 TrcBkOp;
int Score;
int ProbeOfs;
bool bGapOpened;
int TraceBackIdx;
int m_NumWinScores;
tsTraceBackScore TraceBackScores[cTraceBackWin+1];

int WindowAlignScore;

tsSSWTraceback *pTraceBack;

if(MaxArtefactDev < 1)							// <= 0 to disable path classification
	return(0);
if((ProbeEndOfs - ProbeStartOfs) < cTraceBackWin) // need sufficent overlap to check for window deviations
	return(0);

if(m_pAllWinScores == NULL)
	{
	if((m_pAllWinScores = new int [cMaxMAFBlockErrCorLen+1])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: ClassifyPath unable to allocate for windowed traceback scores");
		return(eBSFerrMem);
		}
	}
m_NumWinScores = 0;
if(MaxArtefactDev > 25)			// clamp to be in range 1 to 25
	MaxArtefactDev = 25;

// calculate a score over the path using same scoring as when SW'ing
// now check for localised significant deviations from the overall rate along the target
// using a window of 1Kbp and expecting the peak deviation from overall mean to be less than MaxArtefactDev
Score = 0;
bGapOpened = false;
ProbeOfs = 0;
TraceBackIdx = 0;
pTraceBack = InitiateTraceback(ProbeEndOfs,TargEndOfs);
do {
	switch(TrcBkOp = (pTraceBack->IdxP & cTrBkFlgsMsk)) {
		case cTrBkFlgStart:							// start of alignment path; process as if cTrBkFlgMatch
		case cTrBkFlgMatch:							// match - may not be an exact match
			if(pTraceBack->IdxT & cTrBkFlgSub)
				Score += m_MismatchPenalty;
			else
				Score += m_MatchScore;
			bGapOpened = false;
			ProbeOfs += 1;
			break;	
		case cTrBkFlgIns:							// base inserted into probe relative to target
			ProbeOfs += 1;
		case cTrBkFlgDel:							// base deleted from probe relative to target
			if(!bGapOpened)
				Score += m_GapOpenPenalty;
			else
				Score += m_GapExtnPenalty;
			bGapOpened = true;
		}
	if(Score < 0)
		Score = 0;
	int LatestScoreIdx;
	int OldestScoreIdx;
	LatestScoreIdx = TraceBackIdx % cTraceBackWin;
	TraceBackScores[LatestScoreIdx].Score = Score;
	TraceBackScores[LatestScoreIdx].ProbeOfs = ProbeOfs;
	if(TraceBackIdx >= (cTraceBackWin - 1))
	    {
		OldestScoreIdx = (LatestScoreIdx + 1) % cTraceBackWin;
		WindowAlignScore = (1000 * (TraceBackScores[LatestScoreIdx].Score - TraceBackScores[OldestScoreIdx].Score)) / (TraceBackScores[LatestScoreIdx].ProbeOfs - TraceBackScores[OldestScoreIdx].ProbeOfs);
		if(WindowAlignScore < 0)
			WindowAlignScore = 0;
		if(TraceBackIdx < cMaxMAFBlockErrCorLen)
			m_pAllWinScores[m_NumWinScores++] = WindowAlignScore;
		}
	TraceBackIdx += 1;
	if(TrcBkOp == cTrBkFlgStart || ((pTraceBack->IdxP & cTrBkIdxMsk) == ProbeStartOfs && (pTraceBack->IdxT & cTrBkIdxMsk) == TargStartOfs))
		break;
	}
while(pTraceBack = NxtTraceback(pTraceBack));

if(m_NumWinScores < (cTraceBackWin * 3) / 2)				
	return(0);

// determine mean and min max scores so can then check for significant variances
int WinScoreIdx;
int *pWinScore;
int CurWinScore;
int MaxWinScore;
int MinWinScore;
int MeanWinScore;
INT64 SumWinScores;
 
SumWinScores = 0;
pWinScore = m_pAllWinScores;
SumWinScores = MaxWinScore = MinWinScore = *pWinScore++;
for(WinScoreIdx = 1; WinScoreIdx < m_NumWinScores; WinScoreIdx+=1,pWinScore+=1 )
	{
	CurWinScore = *pWinScore;
	SumWinScores += (INT64)CurWinScore;
	if(CurWinScore > MaxWinScore)
		MaxWinScore = CurWinScore;
	if(CurWinScore < MinWinScore)
		MinWinScore = CurWinScore;
	}
MeanWinScore = (int)(SumWinScores / m_NumWinScores);
if(MaxWinScore > (MeanWinScore * (MaxArtefactDev + 100))/100 || (MinWinScore * (MaxArtefactDev + 100))/100 < MeanWinScore)	
	return(1);			// classify as being artefact
return(0);
}


int													// number of alignment ops generated
CSSW::TracebacksToAlignOps(UINT32 ProbeStartOfs,	// alignment starts at this probe sequence offset (1..n)
					UINT32 ProbeEndOfs,				// alignment ends at this probe sequence offset
					UINT32 TargStartOfs,			// alignment starts at this target sequence offset (1..n)
					UINT32 TargEndOfs)				// alignment ends at this target sequence offset
{
tMAOp AOp;
tMAOp *pOps;
tMAOp *pOpsXchg;
UINT32 TrcBkOp;
int Ofs;
int NumOps;
tsSSWTraceback *pTraceBack;

m_MAAlignOps = 0;
NumOps = 0;
pOps = m_pMAAlignOps;
pTraceBack = InitiateTraceback(ProbeEndOfs,TargEndOfs);
do {
	if((m_MAAlignOps + 10) > m_AllocdMAAlignOpsSize)	// realloc as may be required
		{
		size_t memreq;
		void *pAllocd;
		memreq = (m_AllocdMAAlignOpsSize * 5)/4;        // increase by 25%
 
#ifdef _WIN32
		pAllocd = realloc(m_pMAAlignOps,memreq);
#else
		pAllocd = mremap(m_pMAAlignOps,m_AllocdMAAlignOpsSize,memreq,MREMAP_MAYMOVE);
		if(pAllocd == MAP_FAILED)
			pAllocd = NULL;
#endif
		if(pAllocd == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"TracebacksToAlignOps: Memory re-allocation to %lld bytes - %s",memreq,strerror(errno));
			return(eBSFerrMem);
			}
	
		m_pMAAlignOps = (UINT8 *)pAllocd;
		memset(&m_pMAAlignOps[m_MAAlignOps],0,memreq - m_AllocdMAAlignOpsSize);
		m_AllocdMAAlignOpsSize = memreq;
		}
	switch(TrcBkOp = (pTraceBack->IdxP & cTrBkFlgsMsk)) {
		case cTrBkFlgMatch:							// match - may not be an exact match
			AOp = cMAMatch;
			break;	
		case cTrBkFlgIns:							// base inserted into probe relative to target
			AOp = cMAInsert;
			break;	
		case cTrBkFlgDel:							// base deleted from probe relative to target
			AOp = cMADelete;
			break;
		case cTrBkFlgStart:							// start of alignment path; process as if cTrBkFlgMatch
			AOp = cMAMatch;
			break;
		}
	NumOps += 1;
	*pOps++ = AOp;
	if(TrcBkOp == cTrBkFlgStart || ((pTraceBack->IdxP & cTrBkIdxMsk) == ProbeStartOfs && (pTraceBack->IdxT & cTrBkIdxMsk) == TargStartOfs))
		break;
	}
while(pTraceBack = NxtTraceback(pTraceBack));

// now reverse the alignment ops
pOps = m_pMAAlignOps;
pOpsXchg = &m_pMAAlignOps[NumOps-1];
for(Ofs = 0; Ofs < NumOps/2; Ofs++,pOps++,pOpsXchg--)
	{
	AOp = *pOps;
	*pOps = *pOpsXchg;
	*pOpsXchg = AOp;
	}
m_pMAAlignOps[NumOps] = cMACompleted;
m_MAAlignOps = NumOps;
return(NumOps);
}

int
CSSW::AddMultiAlignment(UINT32 ProbeStartOfs,		// alignment starts at this probe sequence offset (1..n)
					  UINT32 ProbeEndOfs,			// alignment ends at this probe sequence offset inclusive
					  UINT32 TargStartOfs,			// alignment starts at this target sequence offset (1..n)
					  UINT32 TargEndOfs,			// alignment ends at this target sequence offset inclusive
					  etSeqBase *pTargSeq)			// alignment target sequence
{
tMAOp Op;
tMAOp *pAlignOps;
UINT32 NumAlignOps;
int DepthIdx;
tsMAlignCol *pNewCol;
tsMAlignCol *pCol;
tsMAlignCol *pPrevCol;
etSeqBase *pBase;

if(m_MACoverage == m_MADepth)
	return(-1);

pAlignOps = m_pMAAlignOps;
NumAlignOps = m_MAAlignOps;

// iterate over probe columns, seting column base as unaligned, until the ProbeStartOfs column at which alignment starts
pCol = (tsMAlignCol *)m_pMACols;
while(pCol->Ofs != 0 && pCol->Ofs < ProbeStartOfs)
	{
	pCol->Bases[pCol->Depth++] = eBaseUndef;
	if(pCol->NxtColIdx == 0)         // alignment should have started!
		return(-1);
	pCol = (tsMAlignCol *)&m_pMACols[(pCol->NxtColIdx - 1) * (size_t)m_MAColSize];
	}

pBase = &pTargSeq[TargStartOfs-1];

while(NumAlignOps && pCol != NULL)
	{
	Op=*pAlignOps++;
	NumAlignOps -= 1;
	if(Op < cMADelete)				 // if not an InDel then skip over any probe InDel using the probe base as the targets base until a probe match
		{
		while(pCol->Extn != 0)
			{
			pCol->Bases[pCol->Depth++] = pCol->Bases[0];
			if(pCol->NxtColIdx == 0)  // alignment should not be finishing on an InDel!
				return(-1);
			pCol = (tsMAlignCol *)&m_pMACols[(pCol->NxtColIdx-1) * (size_t)m_MAColSize];
			}	
		}

	switch(Op) {
		case cMAMatch:		 // base match between probe and target; note that may not be an exact match
			pCol->Bases[pCol->Depth++] = *pBase++;
			break;

		case cMAInsert:      // base inserted into probe relative to target - or could be base deleted from target relative to probe
			pCol->Bases[pCol->Depth++] = eBaseInDel;
			break;

		case cMADelete:      // base deleted from probe relative to target- or could be base inserted into target relative to probe
			if(pCol->Extn == 0)		// if not positioned on a probe deletion column then will need to create one
				{
				if((pNewCol = InsertCol(pCol->PrevColIdx)) == NULL)
					return(-1);
				pPrevCol = (tsMAlignCol *)&m_pMACols[(pNewCol->PrevColIdx-1) * m_MAColSize];
				pNewCol->Ofs = pPrevCol->Ofs;
				pNewCol->Extn = pPrevCol->Extn + 1;
				pNewCol->ConsBase = eBaseInDel;
				pNewCol->Depth = pPrevCol->Depth - 1;
				for(DepthIdx = 0; DepthIdx < pNewCol->Depth; DepthIdx++)
					pNewCol->Bases[DepthIdx] = pPrevCol->Bases[DepthIdx] == eBaseUndef ? eBaseUndef : eBaseInDel;
				pCol = pNewCol;
				}
			pCol->Bases[pCol->Depth++] = *pBase++;
			break;

	   }

	if(pCol->NxtColIdx)   // if not last
		pCol = (tsMAlignCol *)&m_pMACols[(pCol->NxtColIdx-1) * (size_t)m_MAColSize];
	else
		pCol = NULL;
	}

while(pCol != NULL)
	{
	pCol->Bases[pCol->Depth++] = eBaseUndef;
	if(pCol->NxtColIdx == 0)
		break;
	pCol = (tsMAlignCol *)&m_pMACols[(pCol->NxtColIdx - 1) * (size_t)m_MAColSize];
	}
m_MACoverage += 1;
return(1);
}

// generate multiple alignment consensus from multialignment columns at m_pMACols
// writes consensus base and score back into m_pMACols
// score = TotalBases * (MostAbundant - NextAbundant) / (MostAbundant + NextAbundant)
int
CSSW::GenMultialignConcensus(void)
{
int BaseCnts[eBaseEOS+1];
int AbundIdxs[eBaseEOS+1];
int NumIndelCols;
int *pCnts;
int *pNxtCnts;
int XAbundIdxs;
int TotBaseCnts;
UINT8 Base;
UINT8 *pBase;
int BaseIdx;
tsMAlignCol *pCol;
tsMAlignCol *pInDelStartCol;
tsMAlignCol *pInDelEndCol;

if(!m_bStartedMultiAlignments)
	return(0);

// firstly iterate over columns and where there are deletions of more than 1bp in the probe then try to find the most parsimonious alignment of bases in these insertion columns
pInDelStartCol = NULL;
pInDelEndCol = NULL;
pCol = (tsMAlignCol *)m_pMACols;  
do
	{
	if(pCol->Extn != 0)
		{
		if(pInDelStartCol == NULL)
			pInDelStartCol = pCol;
		pInDelEndCol = pCol;
		}
	else
		{
		if(pInDelStartCol != NULL)
			{
			NumIndelCols = 1 + pInDelEndCol->Extn -pInDelStartCol->Extn; 
			if(NumIndelCols >= 2)			// obtain parsimonious alignments only if InDel is at least 2bp long
				{
				ParsimoniousMultialign(pInDelEndCol->Depth,	NumIndelCols, pInDelStartCol);
				}
			pInDelStartCol = NULL;
			pInDelEndCol = NULL;
			}
		}
	
	if(pCol->NxtColIdx == 0)
		pCol = NULL;
	else
		pCol = (tsMAlignCol *)&m_pMACols[(pCol->NxtColIdx - 1) * (size_t)m_MAColSize];
	}
while(pCol != NULL);

pCol = (tsMAlignCol *)m_pMACols;
do
	{
	memset(BaseCnts,0,sizeof(BaseCnts));
	TotBaseCnts = 0;
	pBase = pCol->Bases;
	for(BaseIdx = 0; BaseIdx < pCol->Depth; BaseIdx++,pBase++)
		{
		if((Base = (*pBase & eBaseEOS)) != eBaseUndef) // don't bother counting the undefined
			{
			BaseCnts[*pBase & eBaseEOS] += 1;
			TotBaseCnts += 1;
			}
		}

	// order base counts by counts descending
	for(BaseIdx = 0; BaseIdx <= eBaseInDel; BaseIdx++)
		AbundIdxs[BaseIdx] = BaseIdx;
	do {
		XAbundIdxs = -1;
		for(BaseIdx = 0; BaseIdx < eBaseInDel; BaseIdx++)
			{
			pCnts = &BaseCnts[AbundIdxs[BaseIdx]];
			pNxtCnts = &BaseCnts[AbundIdxs[BaseIdx+1]];
			if(*pCnts < *pNxtCnts)
				{
				XAbundIdxs = AbundIdxs[BaseIdx];
				AbundIdxs[BaseIdx] = AbundIdxs[BaseIdx+1];
				AbundIdxs[BaseIdx+1] = XAbundIdxs; 
				}
			}
		}
	while(XAbundIdxs != -1);

	int ConsConf;
	int ConsMostAbundCnt;
	int ConsMNextAbundCnt;
	// simply choosing the most abundant base (could be equally abundant!) as being the consensus and then scoring according to relative abundance to next most abundant base
	pCol->ConsBase = AbundIdxs[0] < eBaseInDel ? AbundIdxs[0] : eBaseUndef;
	ConsMostAbundCnt = BaseCnts[AbundIdxs[0]];
	ConsMNextAbundCnt = BaseCnts[AbundIdxs[1]];
	ConsConf = TotBaseCnts * (ConsMostAbundCnt - ConsMNextAbundCnt) / (ConsMostAbundCnt + ConsMNextAbundCnt);
	pCol->ConsConf = ConsConf > 9 ? 9 : ConsConf;
	if(pCol->NxtColIdx == 0)
		pCol = NULL;
	else
		pCol = (tsMAlignCol *)&m_pMACols[(pCol->NxtColIdx - 1) * (size_t)m_MAColSize];
	}
while(pCol != NULL);
return(eBSFSuccess);
}

tsMAlignCol *								// inserted column or NULL if errors inserting
CSSW::InsertCol(UINT32 PrevCol)				// allocate and insert new column after Prev
{
tsMAlignCol *pCol;
tsMAlignCol *pPrevCol;
tsMAlignCol *pNextCol;

if(m_MACols + 10 >= m_AllocMACols)			// need to extend previously allocated columns to hold this new column?
	{
	size_t memreq;
	void *pAllocd;
	int AllocMACols;
	AllocMACols = (m_AllocMACols * 5) / 4;   // increase allocated cols by 25%
	memreq = (size_t)AllocMACols * m_MAColSize;

#ifdef _WIN32
	pAllocd = realloc(m_pMACols,memreq);
#else
	pAllocd = mremap(m_pMACols,m_AllocMAColsSize,memreq,MREMAP_MAYMOVE);
	if(pAllocd == MAP_FAILED)
		pAllocd = NULL;
#endif
	if(pAllocd == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"InsertCol: Memory re-allocation to %lld bytes - %s",memreq,strerror(errno));
		return(NULL);
		}
	m_pMACols = (UINT8 *)pAllocd;
	memset(&m_pMACols[m_AllocMAColsSize],0,memreq - m_AllocMAColsSize);
	m_AllocMACols = AllocMACols;
	m_AllocMAColsSize = memreq;
	}
pCol = (tsMAlignCol *)&m_pMACols[m_MACols++ * m_MAColSize];
pCol->ColIdx = m_MACols;
pPrevCol = (tsMAlignCol *)&m_pMACols[(PrevCol-1) * m_MAColSize];
if(pPrevCol->NxtColIdx != 0)
	{
	pNextCol = (tsMAlignCol *)&m_pMACols[(pPrevCol->NxtColIdx-1) * m_MAColSize];
	pNextCol->PrevColIdx = pCol->ColIdx;
	}
pCol->NxtColIdx = pPrevCol->NxtColIdx;
pCol->PrevColIdx = pPrevCol->ColIdx;
pPrevCol->NxtColIdx = pCol->ColIdx;
return(pCol);
}




UINT64		// returned mask
CSSW::MSBBitMsk(UINT64 BitsSet) // get bit mask of most significant bit set in BitsSet
{
UINT64 Msk;
if(BitsSet == 0)
	return(0);

if(BitsSet & 0x7fff000000000000)
	Msk =    0x4000000000000000;
else
	if(BitsSet & 0x0000ffff00000000)
		Msk =    0x0000800000000000;
	else            
		if(BitsSet & 0x00000000ffff0000)
				Msk =0x0000000080000000;
			else
		if(BitsSet & 0x000000000000ff00)
				Msk =0x0000000000008000;
			else
				Msk =0x0000000000000080;
while(!(BitsSet & Msk))
	Msk >>= 1;
return(Msk);
}


UINT64		// returned mask
CSSW::NSBBitMsk(int Nth,UINT64 BitsSet) // get bit mask of Nth (1 if MSB) significant bit set in BitsSet
{
UINT64 Msk;
if(BitsSet == 0)
	return(0);
while(BitsSet != 0)
	{
	if(BitsSet & 0x7fff000000000000)
		Msk =    0x4000000000000000;
	else
		if(BitsSet & 0x0000ffff00000000)
			Msk =    0x0000800000000000;
		else            
			if(BitsSet & 0x00000000ffff0000)
					Msk =0x0000000080000000;
				else
			if(BitsSet & 0x000000000000ff00)
					Msk =0x0000000000008000;
				else
					Msk =0x0000000000000080;
	while(!(BitsSet & Msk))
		Msk >>= 1;
	if(!--Nth)
		return(Msk);
	// Msk now holds the MSB; set it to 0 and iterate
	BitsSet &= ~Msk;
	}
return(0);
}


int		// number of bits set
CSSW::NumBitsSet(UINT64 BitsSet) // get total number of bits set in BitsSet
{
int BitCnt;
for(BitCnt=0; BitsSet != 0; BitsSet &= BitsSet-1)
	BitCnt += 1;
return(BitCnt);
}

static UINT64 _nCrMsks[] = { 0x0000000000000001, 0x0000000000000003, 0x0000000000000007, 0x000000000000000f,	// masks for {n or r} = 1..4
				      0x000000000000001f, 0x000000000000003f, 0x000000000000007f, 0x00000000000000ff,           // masks for {n or r} = 5..8
				      0x00000000000001ff, 0x00000000000003ff, 0x00000000000007ff, 0x0000000000000fff,			// masks for {n or r} = 9..12
                      0x0000000000001fff, 0x0000000000003fff, 0x0000000000007fff, 0x000000000000ffff,			// masks for {n or r} = 13..16
                      0x000000000001ffff, 0x000000000003ffff, 0x000000000007ffff, 0x00000000000fffff,			// masks for {n or r} = 17..20
                      0x00000000001fffff, 0x00000000003fffff, 0x00000000007fffff, 0x0000000000ffffff,			// masks for {n or r} = 21..24
                      0x0000000001ffffff, 0x0000000003ffffff, 0x0000000007ffffff, 0x000000000fffffff,			// masks for {n or r} = 25..28
                      0x000000001fffffff, 0x000000003fffffff, 0x000000007fffffff, 0x00000000ffffffff,			// masks for {n or r} = 29..32
					  0x00000001ffffffff, 0x00000003ffffffff, 0x00000007ffffffff, 0x0000000fffffffff,			// masks for {n or r} = 32..36
				      0x0000001fffffffff, 0x0000003fffffffff, 0x0000007fffffffff, 0x000000ffffffffff,           // masks for {n or r} = 37..40
				      0x000001ffffffffff, 0x000003ffffffffff, 0x000007ffffffffff, 0x00000fffffffffff,			// masks for {n or r} = 41..44
                      0x00001fffffffffff, 0x00003fffffffffff, 0x00007fffffffffff, 0x0000ffffffffffff,			// masks for {n or r} = 45..48
                      0x0001ffffffffffff, 0x0003ffffffffffff, 0x0007ffffffffffff, 0x000fffffffffffff,			// masks for {n or r} = 49..52
                      0x001fffffffffffff, 0x003fffffffffffff, 0x007fffffffffffff, 0x00ffffffffffffff,			// masks for {n or r} = 53..56
                      0x01ffffffffffffff, 0x03ffffffffffffff, 0x07ffffffffffffff, 0x0fffffffffffffff,			// masks for {n or r} = 57..60
                      0x1fffffffffffffff, 0x3fffffffffffffff, 0x7fffffffffffffff, 0xffffffffffffffff,			// masks for {n or r} = 61..64
					};


UINT64
CSSW::Nxt_nCr(int nPossibilities,		 // from n possibilities
				int rOutComes,			 // choose r outcomes 
			   UINT64 PrevCombination)	 // bit mask of previous combination generated
{
UINT64 nCrMsk1;
UINT64 nCrMsk2;
UINT64 nCrMsk3;

// check if combination already at it's most significant; if so then use the initial least significant combination as the next combination
nCrMsk1 = _nCrMsks[rOutComes-1] << (nPossibilities - rOutComes);
if((nCrMsk2 = (nCrMsk1 & PrevCombination)) == PrevCombination)
	return(_nCrMsks[rOutComes-1]);

nCrMsk3 = MSBBitMsk(PrevCombination);
PrevCombination &= ~nCrMsk3;
nCrMsk3 <<= 1;
if(nCrMsk3 & _nCrMsks[nPossibilities-1])
	{
	PrevCombination |= nCrMsk3;
	return(PrevCombination);
	}

// determine just how many are set in the most significant
nCrMsk1 = (UINT64)0x01 << (nPossibilities-2);
nCrMsk3 >>= 1;
while(nCrMsk1 & nCrMsk2)
	{
	nCrMsk3 |= nCrMsk1;
	nCrMsk1 >>= 1;
	}

// reset the bits
PrevCombination &= ~nCrMsk3;

// get mask of what is now the most significant bit set
nCrMsk1 = MSBBitMsk(PrevCombination);
PrevCombination &= ~nCrMsk1;
nCrMsk1 <<= 1;
while(!(nCrMsk3 & nCrMsk1))
	nCrMsk3 >>= 1;

PrevCombination |= (nCrMsk3 << 1) | nCrMsk1;
return(PrevCombination);
}


bool											// true: sequence has been completely permutated, true, false: intermediate permutation 
CSSW::PermuteInDels(tsPermInDels *pPerm)
{
UINT64 BitMsk;
int BaseIdx;
etSeqBase *pBase;
etSeqBase *pXBase;
etSeqBase TmpBase;
etSeqBase TmpBaseX;

if(pPerm->SeqLen <= 1 ||					// can't permute a single base sequence or 
	pPerm->NumInDels < 1 ||					// if sequence contains no InDels or
	pPerm->NumInDels == pPerm->SeqLen)      // sequence consists of only InDels
	return(true);

pPerm->CurInDelPsns = Nxt_nCr(pPerm->SeqLen,pPerm->NumInDels,pPerm->CurInDelPsns);

// pPerm->CurInDelPsns contains the bitmap of required permutated InDel positions 
pBase = pPerm->InPermSeq;
BitMsk = 0x01;
for(BaseIdx = 0; BaseIdx < pPerm->SeqLen; BaseIdx++, pBase++, BitMsk <<= 1) 
	{
	if(!(BitMsk & pPerm->CurInDelPsns))		// if non-InDel required at this base position
		{
		if(*pBase == eBaseInDel)			// if there is an indel then need to exchange with 1st 3' non-Indel
			{			
			pXBase = pBase+1;		
			while(*pXBase == eBaseInDel)
				pXBase += 1;
			*pBase = *pXBase;
			*pXBase = eBaseInDel;
			}
		}
	else									// else InDel required at this base position
		{
		if(*pBase != eBaseInDel)			// if there is a non-InDel then need to insert
			{
			TmpBaseX = *pBase;
			*pBase = eBaseInDel;
			pXBase = pBase+1;	
			do {	
				TmpBase = *pXBase;
				*pXBase++ = TmpBaseX;
				TmpBaseX = TmpBase;
				}
			while(TmpBaseX != eBaseInDel);
			}
		}
	}

return(pPerm->InitialInDelPsns == pPerm->CurInDelPsns); // if same permutation as initial then completely permutated
}


// The probe sequence has a gap when aligned to other sequences
// This could be because bases have been deleted from the probe sequence or bases have been inserted into one or more aligned to sequences, or a combination of both
// For alignments meeting a few rules, each aligned to sequence will have it's InDel bases moved around relative to all other aligned to sequences with the 
// peak parsimony (proportion of matching bases over all columns) used to determine which combination of InDels to report as being the most parsimonious multialignment.
// The potential combinations of InDels can explode so there is a heuristic rule used to determine if a multalignment is worth the cost of locating the most parsimonious.
// Multiple alignments are only evaluated for maximal parsimony if:
// Rule: mulialigned only if the sequence length is no longer than 10bp;   PacBio deletions are less common than insertions, est. 4 deletions for every 10 insertions. Deletion lengths average around 1-3bp per deletion.
//
typedef struct Tag_sExpCombs {
		int SeqLen;							// for sequences of this length
		int NumInDels;						// an individual sequence has this number of InDels
		UINT64 ExpTotCombs;					// number of InDel combinations for this sequence instance to discover the combination resulting in maximal parsimony  
	} tsExpCombs;

static tsExpCombs _ExpCombinations[] = {
							{2,1,2},	
							{3,1,3},{3,2,3},
							{4,1,4},{4,2,6},{4,3,4},
							{5,1,5},{5,2,10},{5,3,10},{5,4,5},
							{6,1,6},{6,2,15},{6,3,20},{6,4,15},{6,5,6},
							{7,1,7},{7,2,21},{7,3,35},{7,4,35},{7,5,21},{7,6,7},
							{8,1,8},{8,2,28},{8,3,56},{8,4,70},{8,5,56},{8,6,28},{8,7,8},
							{9,1,9},{9,2,36},{9,3,84},{9,4,126},{9,5,126},{9,6,84},{9,7,36},{9,8,9},
							{10,1,10},{10,2,45},{10,3,120},{10,4,210},{10,5,252},{10,6,210},{10,7,120},{10,8,45},{10,9,10},
							{11,1,11},{11,2,55},{11,3,165},{11,4,330},{11,5,462},{11,6,462},{11,7,330},{11,8,165},{11,9,55},{11,10,11},
							{12,1,12},{12,2,66},{12,3,220},{12,4,495},{12,5,792},{12,6,924},{12,7,792},{12,8,495},{12,9,220},{12,10,66},{12,11,12},
							{13,1,13},{13,2,78},{13,3,286},{13,4,715},{13,5,1287},{13,6,1716},{13,7,1716},{13,8,1287},{13,9,715},{13,10,286},{13,11,78},{13,12,13},
							{14,1,14},{14,2,91},{14,3,364},{14,4,1001},{14,5,2002},{14,6,3003},{14,7,3432},{14,8,3003},{14,9,2002},{14,10,1001},{14,11,364},{14,12,91},{14,13,14},
							{15,1,15},{15,2,105},{15,3,455},{15,4,1365},{15,5,3003},{15,6,5005},{15,7,6435},{15,8,6435},{15,9,5005},{15,10,3003},{15,11,1365},{15,12,455},{15,13,105},{15,14,15},
							{16,1,16},{16,2,120},{16,3,560},{16,4,1820},{16,5,4368},{16,6,8008},{16,7,11440},{16,8,12870},{16,9,11440},{16,10,8008},{16,11,4368},{16,12,1820},{16,13,560},{16,14,120},{16,15,16},
							{17,1,17},{17,2,136},{17,3,680},{17,4,2380}, {17,5,6188}, {17,6,12376},{17,7,19448},{17,8,24310}, {17,9,24310}, {17,10,19448}, {17,11,12376},{17,12,6188},{17,13,2380},{17,14,680},{17,15,136},{17,16,17},
							{18,1,18},{18,2,153},{18,3,816},{18,4,3060}, {18,5,8568}, {18,6,18564},{18,7,31824},{18,8,43758}, {18,9,48620}, {18,10,43758}, {18,11,31824},{18,12,18564},{18,13,8568},{18,14,3060},{18,15,816},{18,16,153},{18,17,18},
							{19,1,19},{19,2,171},{19,3,969},{19,4,3876}, {19,5,11628},{19,6,27132},{19,7,50388},{19,8,75582}, {19,9,92378}, {19,10,92378}, {19,11,75582},{19,12,50388},{19,13,27132},{19,14,11628},{19,15,3876},{19,16,969},{19,17,171},{19,18,19},
							{20,1,20},{20,2,190},{20,3,1140},{20,4,4845},{20,5,15504},{20,6,38760},{20,7,77520},{20,8,125970},{20,9,167960},{20,10,184756},{20,11,167960},{20,12,125970},{20,13,77520},{20,14,38760},{20,15,15504},{20,16,4845},{20,17,1140},{20,18,190},{20,19,20},
							{0,0,0} // flags end of _ExpCombinations 					
							};

 
double										// parsimony of returned (pSequences) multiple alignment, 0.0 (min) to 1.0 (max) or -1 if no sequences could be processed 
CSSW::ParsimoniousBasesMultialign(int NumSeqs,	// number of sequences (alignment depth)
		   int SeqLen,						// each sequence is this length 
		   etSeqBase *pSequences)			// input: ptr to each sequence concatenated together, sequences may contain eBaseA, eBaseC,eBaseG,eBaseT,eBaseN, eBaseUndef, and eBaseInDel
											//        NOTE: expectation is that if a sequence contains eBaseUndefs then the complete sequence is eBaseUndefs and if sequence contains any eBaseInDels then these will be at the 3' end of sequence
											// output: pSequences updated with most parsimonious  				
{
int TrimSeqLen;
double PeakParsimonyFactor;
double CurParsimonyFactor;
UINT32 NumAllInDels;
UINT32 NumNoInDels;
UINT32 NumPermutable;
UINT64 ExpTotCombs;
int NumSeqsNoUndefs;
int NxtIdx;
int HiCnts;
tsPermInDels *pPermInDels;
tsPermInDels *pChkPermInDels;
tsPermInDels *pPrevPermInDels;
int PeakSeqIdx;
int BaseIdx;
int SeqIdx;
int CntsIdx;
UINT64 InDelMsk;
UINT64 SeqInDelMsk;
etSeqBase Base;
etSeqBase *pBase;
etSeqBase *pSrcBase;
etSeqBase *pParsePermSeq;
tsExpCombs *pSeqLenCombs;
int BaseCnts[eBaseEOS+1];
int NumBases;
int LongestSeqLen;
int NxtLongestSeqLen;

if(SeqLen < 1 || NumSeqs < 2 || NumSeqs > cMaxMultiAlignSeqs)
	return(-1.0);

if(SeqLen >= 6)
	{
	LongestSeqLen = 0;
	NxtLongestSeqLen = 0;
	for(SeqIdx = 0; SeqIdx < NumSeqs; SeqIdx++)
		{
		NumBases = 0;
		pSrcBase = &pSequences[SeqIdx*SeqLen];
		for(BaseIdx = 0; BaseIdx < SeqLen; BaseIdx++,pSrcBase++)
			if((*pSrcBase & 0x07) <= eBaseN)
				NumBases += 1;
		if(NumBases > 0)
			{
			if(NumBases > LongestSeqLen)
				{
				NxtLongestSeqLen = LongestSeqLen;
				LongestSeqLen = NumBases;
				}
			else
				if(NumBases == LongestSeqLen || NumBases > NxtLongestSeqLen)
					NxtLongestSeqLen = NumBases;
			}
		}
	if(NxtLongestSeqLen == 0 || NxtLongestSeqLen >= cMaxParsimoniousAlignLen)
		return(-1.0);

	TrimSeqLen = min(cMaxParsimoniousAlignLen,LongestSeqLen);
	if(NxtLongestSeqLen < TrimSeqLen / 2)
		TrimSeqLen = NxtLongestSeqLen;
	}
else
	TrimSeqLen = SeqLen;

if(TrimSeqLen == 1)								// if only a single base then can't permutate so simply calculate parsimony
	{ 
	NumSeqsNoUndefs = 0;
	memset(BaseCnts,0,sizeof(BaseCnts));
	pBase = pSequences;
	for(SeqIdx = 0; SeqIdx < NumSeqs; SeqIdx++,pBase+=SeqLen)
		if((Base = *pBase & 0x07) != eBaseUndef)
			{
			BaseCnts[Base] += 1;
			NumSeqsNoUndefs += 1;
			}
	HiCnts = 0;
	for(CntsIdx = 0; CntsIdx <= eBaseInDel; CntsIdx++)
		{
		if(CntsIdx != eBaseUndef)  // undefined bases do not count towards parsimony
			{
			if(BaseCnts[CntsIdx] > HiCnts)
				HiCnts = BaseCnts[CntsIdx];
			}
		}
	if(NumSeqsNoUndefs == 0)        // if all bases were undefined then 
		PeakParsimonyFactor = -1.0; // flag as an error
	else
		if(NumSeqsNoUndefs == 1)		// if only 1 undefined then report as no parsimony
			PeakParsimonyFactor = 0.0;
		else
			PeakParsimonyFactor = HiCnts/(double)NumSeqsNoUndefs;
	return(PeakParsimonyFactor);
	}

// have at least 2 sequences to be permuted
// iterate each sequence to initialise it's corresponding tsPermInDels
SeqInDelMsk = ~(0xffffffffffffffff << (UINT64)TrimSeqLen);
pPermInDels = m_PermIndels;
NumSeqsNoUndefs = 0;
for(SeqIdx = 0; SeqIdx < NumSeqs; SeqIdx++,pPermInDels++)
	{
	pPermInDels->SeqID = SeqIdx+1;
	pPermInDels->NxtSeqID = 0;
	pPermInDels->SeqLen = TrimSeqLen;
	pPermInDels->flgUndef = 0;
	pPermInDels->flgNoPermutate = 0;
	pPermInDels->flgAllInDels = 0;
	pPermInDels->flgRepInstance = 0;
	pPermInDels->flgCharacterised = 0;
	pPermInDels->flgRepInstance = 0;
	pPermInDels->flgNoInDels = 0;
	pPermInDels->NumCopies = 0;
	pPermInDels->NumInDels = 0;
	pPermInDels->InitialInDelPsns = 0;
	InDelMsk = (UINT64)0x01 << (TrimSeqLen - 1);
	pSrcBase = &pSequences[SeqIdx*SeqLen];
	pBase = pPermInDels->InPermSeq;
	pParsePermSeq = pPermInDels->ParsePermSeq;
	NumSeqsNoUndefs += 1;
	for(BaseIdx = 0; BaseIdx < TrimSeqLen; BaseIdx++,pSrcBase++,pParsePermSeq++)
		{
		Base = *pSrcBase & 0x07;
		*pParsePermSeq = Base;
		*pBase++ = Base;
		if(!pPermInDels->flgUndef)
			{
			if(Base == eBaseInDel)
				{
				pPermInDels->NumInDels += 1;
				pPermInDels->InitialInDelPsns |= InDelMsk;
				InDelMsk >>= 1;
				}
			else
				if(Base >= eBaseUndef)	// note: treating eBaseUndef and higher valued as being undefined, eBaseInDel already checked
					{
					pPermInDels->flgUndef = 1;
					pPermInDels->flgNoPermutate = 1;
					pPermInDels->NumInDels = 0;
					pPermInDels->InitialInDelPsns = 0;
					pPermInDels->flgCharacterised = 1;
					NumSeqsNoUndefs -= 1;
					}
			}
		}
	pPermInDels->SeqInDelMsk = SeqInDelMsk;
	pPermInDels->CurInDelPsns = pPermInDels->InitialInDelPsns;
	}

if(NumSeqsNoUndefs == 0)  // if all bases were undefined then flag as being error
	return(-1.0);
if(NumSeqsNoUndefs <= 2)  // if no more than 2 sequences marked as having at least one InDel then don't bother permuting
	return(0.0);

// need at least 3 sequences to characterise each sequence according to if it contains no InDels, containing all InDels and the number of exactly matching sequences
NumAllInDels = 0;
NumNoInDels = 0;
NumPermutable = 0;
pPermInDels = m_PermIndels;
for(SeqIdx = 0; SeqIdx < NumSeqs; SeqIdx++,pPermInDels++)
	{
	if(pPermInDels->flgCharacterised == 1)
		continue;

	if(pPermInDels->NumInDels == TrimSeqLen)
		{
		NumAllInDels += 1;
		pPermInDels->flgNoPermutate = 1;
		pPermInDels->flgAllInDels = 1;
		}
	else
		if(pPermInDels->NumInDels == 0)
			{
			NumNoInDels += 1;
			pPermInDels->flgNoPermutate = 1;
			pPermInDels->flgNoInDels = 1;
			}
	pPermInDels->NumCopies = 1;
	pChkPermInDels = pPermInDels + 1;
	for(NxtIdx = SeqIdx+1; NxtIdx < NumSeqs; NxtIdx++,pChkPermInDels++)
		{
		if(!pChkPermInDels->flgCharacterised && pPermInDels->InPermSeq[0] == pChkPermInDels->InPermSeq[0] && !memcmp(pPermInDels->InPermSeq,pChkPermInDels->InPermSeq,TrimSeqLen))
			{
			if(pPermInDels->NxtSeqID == 0)
				pPermInDels->NxtSeqID = pChkPermInDels->SeqID;
			else
				pPrevPermInDels->NxtSeqID = pChkPermInDels->SeqID;	
			pChkPermInDels->flgNoPermutate = 1;
			pChkPermInDels->flgAllInDels = pPermInDels->flgAllInDels;
			pChkPermInDels->flgNoInDels = pPermInDels->flgNoInDels;
			pChkPermInDels->flgCharacterised = 1;
			pPrevPermInDels = pChkPermInDels;
			pPermInDels->NumCopies += 1;
			}
		}
	pPermInDels->flgCharacterised = 1;
	pPermInDels->flgRepInstance = 1;
	if(pPermInDels->flgNoPermutate == 0)
		NumPermutable += 1;
	}
CurParsimonyFactor = 1.0 - (double)NumPermutable/NumSeqsNoUndefs;
if(CurParsimonyFactor >= 0.75)		// if less than 25% of sequences are permutable then accept
	return(CurParsimonyFactor);

// check if it would take too long to explore all possible conbinations of InDels in the unique sequences
if(SeqLen >= 4)
	{
	pSeqLenCombs = _ExpCombinations;
	while(pSeqLenCombs->SeqLen != TrimSeqLen && pSeqLenCombs->SeqLen != 0)
		pSeqLenCombs += 1;
	if(pSeqLenCombs->SeqLen == 0)
		return(CurParsimonyFactor);
	ExpTotCombs = 1;
	pPermInDels = m_PermIndels;
	for(SeqIdx = 0; SeqIdx < NumSeqs; SeqIdx++,pPermInDels++)
		{
		if(pPermInDels->flgAllInDels || !pPermInDels->flgRepInstance || (pPermInDels->flgRepInstance && pPermInDels->flgNoInDels))
			continue;
		ExpTotCombs *= pSeqLenCombs[pPermInDels->NumInDels-1].ExpTotCombs;
		if(ExpTotCombs > cMaxAcceptExpCombs)
			return(CurParsimonyFactor);
		}
	}

// start permutating the relative positions of the InDels
PeakSeqIdx = 0;
PeakParsimonyFactor = 0.0;
do {
	// Check the permutations here and if the current permutation is the most parsimonious then copy the permutated sequence into OutPermSeq
	// The most parsimonious is that permutation which maximises the number of loci for which there is a single majority base
	CurParsimonyFactor = 1.0;
	for(BaseIdx = 0; BaseIdx < TrimSeqLen; BaseIdx++)
		{
		memset(BaseCnts,0,sizeof(BaseCnts));
		pPermInDels = m_PermIndels;
		for(SeqIdx = 0; SeqIdx < NumSeqs; SeqIdx++,pPermInDels++)
			{
			if(!pPermInDels->flgUndef && pPermInDels->flgRepInstance)
				BaseCnts[pPermInDels->InPermSeq[BaseIdx]] += pPermInDels->NumCopies;
			}
		HiCnts = 0;
		for(CntsIdx = 0; CntsIdx <= eBaseInDel; CntsIdx++)
			{
			if(BaseCnts[CntsIdx] > HiCnts)
				HiCnts = BaseCnts[CntsIdx];
			}
		CurParsimonyFactor *= HiCnts/(double)NumSeqsNoUndefs;
		}
	if(CurParsimonyFactor > PeakParsimonyFactor)
		{
		pPermInDels = m_PermIndels;
		for(SeqIdx = 0; SeqIdx < NumSeqs; SeqIdx++,pPermInDels++)
			{
			if(pPermInDels->flgRepInstance && !(pPermInDels->flgAllInDels || pPermInDels->flgNoInDels))
				{
				pChkPermInDels = pPermInDels;
				do {
					memcpy(pChkPermInDels->ParsePermSeq,pPermInDels->InPermSeq,TrimSeqLen);
					if(pChkPermInDels->NxtSeqID > 0)
						pChkPermInDels = &m_PermIndels[pChkPermInDels->NxtSeqID-1];
					else
						pChkPermInDels = NULL;
					}
				while(pChkPermInDels != NULL);
				}
			}
		PeakParsimonyFactor = CurParsimonyFactor;
		}

	pPermInDels = m_PermIndels;
	for(SeqIdx = 0; SeqIdx < NumSeqs; SeqIdx++,pPermInDels++)
		{
		if(!pPermInDels->flgRepInstance || pPermInDels->flgNoPermutate)
			continue;
		if(SeqIdx > PeakSeqIdx)
			PeakSeqIdx = SeqIdx;
		if(!PermuteInDels(pPermInDels))
			break;
		}
	}
while(SeqIdx < NumSeqs);

pPermInDels = m_PermIndels;
for(SeqIdx = 0; SeqIdx < NumSeqs; SeqIdx++,pPermInDels++)
	{
	memcpy(pSequences,pPermInDels->ParsePermSeq,TrimSeqLen);
	pSequences+=SeqLen;
	}
return(PeakParsimonyFactor);
}

double												// parsimony of multiple alignment, 0 (min) to 1.0 (max) 
CSSW::ParsimoniousMultialign(int Depth,				// parsimony for this number of bases in each column
							int NumCols,			// alignment parsimony over this many columns
							tsMAlignCol *pInDelStartCol) // starting at this column
{
double Score;
UINT8 *pColBase;
int ColIdx;
int BaseIdx;
tsMAlignCol *pCol;
int ConcatSeqLen;
etSeqBase *pParsimoniousBase;

ConcatSeqLen = NumCols * Depth;
if(m_pParsimoniousBuff == NULL || (ConcatSeqLen + 10) > m_AllocParsimoniousSize)
	{
	if(m_pParsimoniousBuff != NULL)
		delete m_pParsimoniousBuff;
	m_AllocParsimoniousSize = ConcatSeqLen + 1000;
	m_pParsimoniousBuff = new etSeqBase [m_AllocParsimoniousSize];
	}

ColIdx = 0;
pCol = pInDelStartCol;
do {
	pColBase = &pCol->Bases[0];
	pParsimoniousBase = &m_pParsimoniousBuff[ColIdx++];
	for(BaseIdx = 0; BaseIdx < Depth; BaseIdx++,pColBase++,pParsimoniousBase+=NumCols)
		*pParsimoniousBase = *pColBase;
	if(ColIdx == NumCols || pCol->NxtColIdx == 0)
		pCol = NULL;
	else
		pCol = (tsMAlignCol *)&m_pMACols[(pCol->NxtColIdx - 1) * (size_t)m_MAColSize];
	}
while(pCol != NULL);

Score = ParsimoniousBasesMultialign(Depth,NumCols,m_pParsimoniousBuff);

// copy back sequences back into columns
ColIdx = 0;
pCol = pInDelStartCol;
do {
	pColBase = &pCol->Bases[0];
	pParsimoniousBase = &m_pParsimoniousBuff[ColIdx++];
	for(BaseIdx = 0; BaseIdx < Depth; BaseIdx++,pColBase++,pParsimoniousBase+=NumCols)
		*pColBase = *pParsimoniousBase;
	if(ColIdx == NumCols || pCol->NxtColIdx == 0)
		pCol = NULL;
	else
		pCol = (tsMAlignCol *)&m_pMACols[(pCol->NxtColIdx - 1) * (size_t)m_MAColSize];
	}
while(pCol != NULL);
return(Score);
}


const UINT32 cQSWPGap =     0x80000000;      // if set then gap has been opened in probe
const UINT32 cQSWTGap =     0x40000000;      // if set then gap has been opened in targ
const UINT32 cQSWScoreMsk = 0x00ffffff; 

int								// highest score for path anchored at the pProbe[0] and pTarg[0]  
CSSW::ChkStartPath(int Scorethres, // looking for a path of at least this SW score
		  int SeqLen,			// both pProbe and pTarg are at least this length, max of 100 will be processed
		  etSeqBase *pProbe,	// SW this sequence against
		  etSeqBase *pTarg,		// this sequence
		  int SeedScore,		// seed the path with this initial score >= 1
		  int MatchScore,		// exact match score
          int MismatchPenalty,	// mismatch penalty
		  int GapOpenPenalty,	// gap open penalty
		  int GapExtnPenalty)   // gap extension penalty
{
int IdxP;
int IdxT;
etSeqBase *pTBase;
etSeqBase PBase;
etSeqBase TBase;

bool bMatch;
int DiagScore;
int PrevScore;
int LeftScore;
int DownScore;
int HighestScore;

UINT32 LeftCell;
UINT32 DiagCell;
UINT32 PrevCell;
UINT32 *pCell;
UINT32 QSWCells[100+1];

memset(QSWCells,0,sizeof(QSWCells));
LeftCell = 0;
DiagCell = 0;
HighestScore = 0;
for(IdxP = 0; IdxP < SeqLen; IdxP++)
	{
	PBase = *pProbe++ & ~cRptMskFlg;
	pTBase = pTarg;
	pCell = QSWCells;
	for(IdxT = 0; IdxT < SeqLen; IdxT++)
		{
		TBase = *pTBase++ & ~cRptMskFlg;
		bMatch = PBase == TBase;
		if(IdxT > 0 && IdxP > 0)
			{
			DiagCell = LeftCell;
			LeftCell = *pCell;
			PrevScore = DiagCell & cQSWScoreMsk;
			DiagScore = PrevScore + (bMatch ? MatchScore : MismatchPenalty);
			}
		else
			{
			LeftCell = 0;
			DiagCell = 0;
			if(IdxT == 0 && IdxP == 0)
				{
				if(bMatch)
					*pCell = SeedScore + m_MatchScore;
				else
					*pCell = SeedScore + MismatchPenalty;
				}
			else
				*pCell = 0;
			continue;
			}


		// leftscore is either GapExtnPenalty (if gap already opened and gap at least m_DlyGapExtn) or GapOpenScore added to prev left score
        // leftscore relative to diagscore and downscore determines if an insertion into probe is called
		int GapExtnPenalty;
		PrevCell = LeftCell;
		PrevScore = PrevCell & cQSWScoreMsk;
		GapExtnPenalty = 0;
		if(PrevCell & cQSWPGap)  // was there a gap previously opened?
			LeftScore = PrevScore + GapExtnPenalty; 
		else
			LeftScore = PrevScore + m_GapOpenPenalty;

		// down score is either GapExtnPenalty (if gap already opened) or GapOpenScore added to prev down score
        // downscore relative to diagscore and leftscore determines if a deletion from probe is called
		PrevCell = pCell[-1];
		PrevScore = PrevCell & cQSWScoreMsk;
		GapExtnPenalty = 0;
		if(PrevCell & cQSWTGap)  // was there a gap previously opened?
			DownScore = PrevScore + GapExtnPenalty; 
		else
			DownScore = PrevScore + m_GapOpenPenalty;
		
		// if no score was > 0 then reset cell, path can't be extended
		if(DiagScore <= 0 && DownScore <= 0 && LeftScore <= 0)
			{
			*pCell=0;	
			continue;		
			}

		if(DiagScore >= DownScore && DiagScore >= LeftScore)
			*pCell = (UINT32)DiagScore;
		else
			if(DownScore >= LeftScore)
				*pCell = (UINT32)DownScore | cQSWTGap;
			else
				*pCell = (UINT32)LeftScore | cQSWPGap;
		}

	}
pCell = QSWCells;
for(IdxP = 0; IdxP < SeqLen; IdxP++,pCell++)
	{
	if((DiagScore = (*pCell & cQSWScoreMsk)) > HighestScore)
		HighestScore = DiagScore;
	}
return(HighestScore);
}
