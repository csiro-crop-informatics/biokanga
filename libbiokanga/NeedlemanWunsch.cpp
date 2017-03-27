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

CNeedlemanWunsch::CNeedlemanWunsch(void)
{
m_pTrcBckCells = NULL;
m_pProbe = NULL;
m_pTarg = NULL;
m_pColBands = NULL;
Reset();
}

CNeedlemanWunsch::~CNeedlemanWunsch(void)
{
Reset();
}

// Reset
// Resets instance to that immediately following that of instance construction
// Note that all buffers will be deallocated..
void
CNeedlemanWunsch::Reset(void)
{
if(m_pTrcBckCells != NULL)
	{
#ifdef _WIN32
	free(m_pTrcBckCells);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pTrcBckCells != MAP_FAILED)
		munmap(m_pTrcBckCells,m_TrcBckCellsAllocdSize);
#endif
	}
m_TrcBckCellsAllocdSize = 0;
m_TrcBckCellsAllocd = 0;	
m_TrcBckCellsUsed = 0;		

if(m_pProbe != NULL)
	{
	delete m_pProbe;
	m_pProbe = NULL;
	}

m_ProbeAllocd = 0;	// actual m_pProbe alloc'd size (in teNWSeqBases's)
m_ProbeLen = 0;			// current probe length

if(m_pColBands != NULL)			// each column in the matrix will contain a band of active cells
	{
	delete m_pColBands;
	m_pColBands = NULL;
	}
m_ColBandsAllocd = 0;
m_ColBandsUsed = 0;

if(m_pTarg != NULL)				// alloc'd target sequence memory
	{
	delete m_pTarg;
	m_pTarg = NULL;
	}
m_TargAllocd = 0;			// actual m_pTarg alloc'd size (in teNWSeqBases's)
m_TargLen = 0;				// current targ length

m_MatchScore = cNWDfltMatchScore;			// score for matching bases
m_MismatchScore = cNWDfltMismatchScore;		// mismatch penalty
m_GapOpenScore = cNWDfltGapOpenScore;		// gap opening penalty
m_GapExtScore = cNWDfltGapExtScore;			// gap extension penalty

m_PeakScore = 0;

m_NumBasesAligned = 0;
m_NumBasesExact = 0;
m_NumProbeInserts = 0;
m_NumTargInserts = 0;
m_ProbeAlignStartOfs = 0;
m_TargAlignStartOfs = 0;
m_bBanded = false;
m_bAligned = false;
}

// SetScores
// Set match, mismatch, gap opening and gap extension scores
// Note that a check is made to ensure that MatchScore is > 0 and
// that all other scores are <= 0
bool 
CNeedlemanWunsch::SetScores(int MatchScore,int MismatchScore,int GapOpenScore,int GapExtScore)
{
if(MatchScore <= 0 || MismatchScore > 0 || GapOpenScore > 0 || GapExtScore > 0)
	return(false);
m_MatchScore = MatchScore;				// score for matching bases
m_MismatchScore = MismatchScore;		// mismatch penalty
m_GapOpenScore = GapOpenScore;			// gap opening penalty
m_GapExtScore = GapExtScore;			// gap extension penalty
m_ColBandsUsed = 0;
m_TrcBckCellsUsed = 0;
m_NumBasesAligned = 0;
m_NumBasesExact = 0;
m_NumProbeInserts = 0;
m_NumTargInserts = 0;
m_ProbeAlignStartOfs = 0;
m_TargAlignStartOfs = 0;
m_bAligned = false;
m_bAligned = false;
return(true);
}

// SetProbe
// Set probe sequence to use in subsequent alignments
// Probe sequence is used as the row increment
// Only restriction is that the length must be in the range of cNWMinProbeTargLen to cNWMaxProbeTargLen
bool 
CNeedlemanWunsch::SetProbe(UINT32 Len,etSeqBase *pSeq)
{
if(Len < cNWMinProbeOrTargLen || Len > cNWMaxProbeOrTargLen || pSeq == NULL) 	// can't be bothered with very short probes!
	return(false);

if(m_pColBands == NULL || m_ColBandsAllocd < Len || m_ColBandsAllocd > Len * 2)
	{
	if(m_pColBands != NULL)
		delete m_pColBands;
	m_ColBandsAllocd = Len + 100;
	m_pColBands = new tsNWColBand [m_ColBandsAllocd];
	if(m_pColBands == NULL)
		return(false);
	}

if(m_pProbe == NULL || m_ProbeAllocd < Len)
	{
	if(m_pProbe != NULL)
		delete m_pProbe;
	m_ProbeAllocd = Len + 100;
	m_pProbe = new etSeqBase [Len + 100];
	if(m_pProbe == NULL)
		return(false);
	}
memmove(m_pProbe,pSeq,Len);
m_ProbeLen = Len;
m_ColBandsUsed = 0;
m_TrcBckCellsUsed = 0;
m_NumBasesAligned = 0;
m_NumBasesExact = 0;
m_NumProbeInserts = 0;
m_NumTargInserts = 0;
m_ProbeAlignStartOfs = 0;
m_TargAlignStartOfs = 0;
m_bAligned = false;
return(true);
}

// SetTarg
// Set target sequence to use in subsequent alignments
// Target sequence is used as the column increment
// Only restriction is that the length must be in the range of cNWMinProbeTargLen to cNWMaxProbeTargLen
bool 
CNeedlemanWunsch::SetTarg(UINT32 Len,etSeqBase *pSeq)
{
if(Len < cNWMinProbeOrTargLen || Len > cNWMaxProbeOrTargLen || pSeq == NULL)	// can't be bothered with very short targets!
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
m_ColBandsUsed = 0;
m_TrcBckCellsUsed = 0;
m_NumBasesAligned = 0;
m_NumBasesExact = 0;
m_NumProbeInserts = 0;
m_NumTargInserts = 0;
m_ProbeAlignStartOfs = 0;
m_TargAlignStartOfs = 0;
m_bAligned = false;
return(true);
}

// Align
// Align probe SetProbe() against target SetTarg() using scores specified via SetScores() 
// Needleman-Wunsch style global alignment dynamic programing with optimisations
// Returns the peak score of all aligned subsequences
int									// peak score of all subsequence alignments
CNeedlemanWunsch::Align(void)
{
UINT32 IdxP;						// current index into m_Probe[]
UINT32 IdxT;						// current index into m_Targ[]
UINT32 NumCells;					// m_ProbeLen * m_TargLen - total number of cells
int NewScore;
int DiagScore;							// putative diagonal score
int LeftScore;							// putative left score
int DownScore;							// putative down score
int PrevScore;							// score in back referenced cell
tNWTrcBckCell DiagDir;					// to hold back reference direction as diagonal  

etSeqBase *pProbe;						// current &m_Probe[IdxP]
etSeqBase ProbeBase;					// current m_Probe[IdxP]

etSeqBase *pTarg;						// current &m_Targ[IdxT]
etSeqBase TargBase;						// current m_Targ[IdxT]

tNWTrcBckCell *pCell;
tNWTrcBckCell *pPrevCell;
bool bMatch;

m_bAligned = false;
if(m_ProbeLen < cNWMinProbeOrTargLen || m_TargLen < cNWMinProbeOrTargLen || ((INT64)m_ProbeLen * (INT64)m_TargLen > (INT64)cNWMaxCells))
	return(eBSFerrMaxEntries);	

NumCells = m_ProbeLen * m_TargLen;

if(m_pTrcBckCells == NULL || m_TrcBckCellsAllocd < NumCells || ((UINT64)m_TrcBckCellsAllocd > (UINT64)NumCells * 2))
	{
	NumCells += 100;						// small overallocation as a saftety margin
	if(m_pTrcBckCells != NULL)
		{
#ifdef _WIN32
		free(m_pTrcBckCells);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(m_pTrcBckCells != MAP_FAILED)
			munmap(m_pTrcBckCells,m_TrcBckCellsAllocdSize);
#endif
		}
	m_TrcBckCellsAllocd = NumCells;
	m_TrcBckCellsAllocdSize = sizeof(tNWTrcBckCell) * m_TrcBckCellsAllocd;
#ifdef _WIN32
	m_pTrcBckCells = (tNWTrcBckCell *) malloc(m_TrcBckCellsAllocdSize);
	if(m_pTrcBckCells == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for traceback cells",(INT64)m_TrcBckCellsAllocdSize);
		m_TrcBckCellsAllocdSize = 0;
		m_TrcBckCellsAllocd = 0;
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pTrcBckCells = (tNWTrcBckCell *)mmap(NULL,m_TrcBckCellsAllocdSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pTrcBckCells == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for traceback cells",(INT64)m_TrcBckCellsAllocdSize);
		m_TrcBckCellsAllocdSize = 0;
		m_TrcBckCellsAllocd = 0;
		return(eBSFerrMem);
		}
#endif
	}
m_TrcBckCellsUsed = NumCells;


m_PeakScore = 0;
m_NumBasesAligned = 0;
m_NumBasesExact = 0;
m_ProbeAlignStartOfs = 0;
m_TargAlignStartOfs = 0;
pProbe = m_pProbe;
pCell = m_pTrcBckCells;
for(IdxP = 0; IdxP < m_ProbeLen; IdxP++)
	{
	ProbeBase = *pProbe++ & ~cRptMskFlg;
	pCell = &m_pTrcBckCells[(IdxP * m_TargLen)];
	pTarg = &m_pTarg[0];
	for(IdxT = 0; IdxT < m_TargLen; IdxT++, pCell++)
		{
		TargBase = *pTarg++ & ~cRptMskFlg;
		bMatch = ProbeBase == TargBase;

		// calc the 3 scores
		// diagonal is either MatchScore or MismatchScore added to prev diagonal score
		if(IdxT > 0 && IdxP > 0)
			{
			pPrevCell = pCell - m_TargLen - 1;
			PrevScore = (*pPrevCell & cNWScoreMsk);
			if(PrevScore & cNWScoreNegFlg)			// if sign of score is negative then sign extend to left making score into negative int
				PrevScore |= ~cNWScoreMsk;
			DiagScore = PrevScore + (bMatch ? m_MatchScore : m_MismatchScore);
			}
		else   // either first column or row, these use implied pre-existing scores related to the cell distance from the matrix origin 
			{
			DiagScore = bMatch ? m_MatchScore : m_MismatchScore;
			
			if(IdxT > 0)
				{
				pPrevCell = pCell - 1;
				PrevScore = (*pPrevCell & cNWScoreMsk);
				if(PrevScore & cNWScoreNegFlg)			// if sign of score is negative then sign extend to left making score into negative int
					PrevScore |= ~cNWScoreMsk;
				DiagScore += m_GapOpenScore;
				DiagScore += PrevScore;
				DiagDir = cNWTrcBckDownFlg;
				}
			else
				if(IdxP > 0)
					{
					pPrevCell = pCell - m_TargLen;
					PrevScore = (*pPrevCell & cNWScoreMsk);
					if(PrevScore & cNWScoreNegFlg)			// if sign of score is negative then sign extend to left making score into negative int
						PrevScore |= ~cNWScoreMsk;
					DiagScore += m_GapOpenScore;
					DiagScore += PrevScore;
					DiagDir = cNWTrcBckLeftFlg;
					}
				else
					DiagDir = 0;
			*pCell = (DiagScore & cNWScoreMsk) | DiagDir | (bMatch ? cNWTrcBckMatchFlg : 0);
			if((IdxT == 0 && IdxP == 0) || DiagScore > m_PeakScore)
				m_PeakScore = DiagScore;
			continue;
			}

		// leftscore is either GapExtScore (if gap already opened) or GapOpenScore added to prev left score
		pPrevCell = pCell - m_TargLen;
		PrevScore = (*pPrevCell & cNWScoreMsk);
		if(PrevScore & cNWScoreNegFlg)			// if sign of score is negative then sign extend to left making score into negative int
			PrevScore |= ~cNWScoreMsk;
		LeftScore = PrevScore + ((*pPrevCell & cNWGapOpnFlg) ? m_GapExtScore : m_GapOpenScore);

		// down score is either GapExtScore (if gap already opened) or GapOpenScore added to prev down score
		pPrevCell = pCell - 1;
		PrevScore = (*pPrevCell & cNWScoreMsk);
		if(PrevScore & cNWScoreNegFlg)			// if sign of score is negative then sign extend to left making score into negative int
			PrevScore |= ~cNWScoreMsk;
		DownScore = PrevScore + ((*pPrevCell & cNWGapOpnFlg) ? m_GapExtScore : m_GapOpenScore);
		
		// select highest score into cell together with traceback and gap opened flag..
		if(DiagScore >= DownScore && DiagScore >= LeftScore)
			{
			*pCell = (DiagScore & cNWScoreMsk) | cNWTrcBckDiagFlg | (bMatch ? cNWTrcBckMatchFlg : 0);
			NewScore = DiagScore;
			}
		else
			if(DownScore >= LeftScore)
				{			
				*pCell = (DownScore & cNWScoreMsk) | cNWTrcBckDownFlg | (bMatch ? 0 : cNWGapOpnFlg);
				NewScore = DownScore;
				}
			else
				{			
				*pCell = (LeftScore & cNWScoreMsk) | cNWTrcBckLeftFlg | (bMatch ? 0 : cNWGapOpnFlg);
				NewScore = LeftScore;
				}


		if(NewScore > m_PeakScore)
			m_PeakScore = NewScore;
		}
	}
m_bAligned = true;
return(m_PeakScore);
}


// GetNumAlignedBases
// get number of bases which align (exactly plus subs), excluding InDels
// Also internally updates m_ProbeAlignStart and m_TargAlignStart
int
CNeedlemanWunsch::GetNumAlignedBases(void)	 // get number of bases which align (exactly plus subs), excluding InDels
{
tNWTrcBckCell TrcBckDir;
tNWTrcBckCell *pPeakCell;
UINT32 ProbeIdx;
UINT32 TargIdx;

if(!m_bAligned)
	return(0);

if(m_NumBasesAligned > 0)
	return(m_NumBasesAligned);

m_ProbeAlignStartOfs = m_ProbeLen - 1;
m_TargAlignStartOfs = m_TargLen - 1;
m_NumBasesAligned = 0;
m_NumBasesExact = 0;
m_NumProbeInserts = 0;
m_NumTargInserts = 0;
if(m_bBanded)
	{
	ProbeIdx = m_ProbeLen;
	TargIdx =  m_TargLen;
	pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
	}
else
	pPeakCell = &m_pTrcBckCells[(m_ProbeLen * m_TargLen) - 1];
do {
	switch(TrcBckDir = (*pPeakCell & cNWTrcBckMsk)) {
		case cNWTrcBckDiagFlg:		// back on the diagonal, either a match or mismatch
			m_NumBasesAligned += 1;
			if(*pPeakCell & cNWTrcBckMatchFlg)
				m_NumBasesExact += 1;
			if(m_bBanded)
				{
				ProbeIdx -= 1;
				TargIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= m_TargLen + 1;
			m_ProbeAlignStartOfs -= 1;
			m_TargAlignStartOfs -= 1;
			break;

		case cNWTrcBckLeftFlg:			// left, insertion into probe or deletion from target
			if(m_bBanded)
				{
				ProbeIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= m_TargLen;		// treating as insertion into probe
			m_ProbeAlignStartOfs -= 1;
			m_NumProbeInserts += 1;
			break;

		case cNWTrcBckDownFlg:			// down, insertion into target or deletion 
			if(m_bBanded)
				{
				TargIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else 
				pPeakCell -= 1;				// treating as insertion into target
			m_TargAlignStartOfs -= 1;
			m_NumTargInserts += 1;
			break;

		default:					// cell has no trace back so must be final cell in path
			m_NumBasesAligned += 1;
			if(*pPeakCell & cNWTrcBckMatchFlg)
				m_NumBasesExact += 1;
			break;		
		}
	}
while(TrcBckDir);
return(m_NumBasesAligned);
}


int 
CNeedlemanWunsch::GetProbeAlign(UINT32 Len, etSeqBase *pBuff) // get probe alignment
{
tNWTrcBckCell TrcBckDir;
tNWTrcBckCell *pPeakCell;
etSeqBase *pProbe;
etSeqBase *pProbeStart;
UINT32 ProbeIdx;
UINT32 TargIdx;

if(!GetNumAlignedBases())
	return(0);
if(Len < (m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts))
	return(0);

if(m_bBanded)
	{
	ProbeIdx = m_ProbeLen;
	TargIdx =  m_TargLen;
	pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
	}
else
	pPeakCell = &m_pTrcBckCells[(m_ProbeLen * m_TargLen) - 1];

pProbeStart = m_pProbe + m_ProbeAlignStartOfs;
pProbe = pProbeStart + m_NumBasesAligned + m_NumProbeInserts - 1;
pBuff += m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts - 1;

do
	{
	TrcBckDir = *pPeakCell  & cNWTrcBckMsk;
	switch(TrcBckDir) {
		case cNWTrcBckDiagFlg:			// back on the diagonal, report base
			*pBuff-- = *pProbe--; 
			if(m_bBanded)
				{
				ProbeIdx -= 1;
				TargIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= m_TargLen + 1;
			break;

		case cNWTrcBckLeftFlg:			// left, treating as insertion into probe so report base
			*pBuff-- = *pProbe--; 
			if(m_bBanded)
				{
				ProbeIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= m_TargLen;
			break;

		case cNWTrcBckDownFlg:		     // down, treating as insertion into target so report as InDel
			*pBuff-- = eBaseInDel; 
			if(m_bBanded)
				{
				TargIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= 1;
			break;

		default:					// no direction, must be final cell so report probe base
			*pBuff = *pProbe; 
			break;	
		}
	}
while(TrcBckDir);
return(m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts);
}


int 
CNeedlemanWunsch::GetTargAlign(UINT32 Len, etSeqBase *pBuff) // get target alignment
{
tNWTrcBckCell TrcBckDir;
tNWTrcBckCell *pPeakCell;
etSeqBase *pTarg;
etSeqBase *pTargStart;
UINT32 ProbeIdx;
UINT32 TargIdx;

if(!GetNumAlignedBases())
	return(0);
if(Len < (m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts))
	return(0);
if(m_bBanded)
	{
	ProbeIdx = m_ProbeLen;
	TargIdx =  m_TargLen;
	pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
	}
else
	pPeakCell = &m_pTrcBckCells[(m_ProbeLen * m_TargLen) - 1];
pTargStart = m_pTarg + m_TargAlignStartOfs;
pTarg = pTargStart +  + m_NumBasesAligned + m_NumTargInserts - 1;
pBuff += m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts - 1;

do
	{
	TrcBckDir = *pPeakCell  & cNWTrcBckMsk;
	switch(TrcBckDir) {
		case cNWTrcBckDiagFlg:		// back on the diagonal, report base
			*pBuff-- = *pTarg--; 
			if(m_bBanded)
				{
				ProbeIdx -= 1;
				TargIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= m_TargLen + 1;
			break;

		case cNWTrcBckLeftFlg:			// left, treating as insertion into probe so report InDel
			*pBuff-- = eBaseInDel; 
			if(m_bBanded)
				{
				ProbeIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else			
				pPeakCell -= m_TargLen;
			break;

		case cNWTrcBckDownFlg:		     // down, treating as insertion into target so report as target base
			*pBuff-- = *pTarg--; 
			if(m_bBanded)
				{
				TargIdx -= 1;
				pPeakCell = DerefBandCell(ProbeIdx,TargIdx);
				}
			else
				pPeakCell -= 1;			
			break;

		default:				// no direction, must be final cell so report target base
			*pBuff = *pTarg; 
			break;	
		}
	}
while(TrcBckDir);
return(m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts);
}


int												// returned total alignment length between probe and target including InDels
CNeedlemanWunsch::GetAlignStats(UINT32 *pNumAlignedBases,			// returned number of bases aligning between probe and target
				 UINT32 *pNumExactBases,           // of the aligning bases there were this many exact matches, remainder were substitutions
				 UINT32 *pNumProbeInsertBases,     // this many bases were inserted into the probe relative to the target
				 UINT32 *pNumTargInsertBases)		// this many bases were inserted into the target relative to the probe
{
if(!GetNumAlignedBases())
	return(-1);
if(pNumAlignedBases != NULL)
	*pNumAlignedBases = m_NumBasesAligned;
if(pNumExactBases != NULL)
	*pNumExactBases = m_NumBasesExact;
if(pNumProbeInsertBases != NULL)
	*pNumProbeInsertBases = m_NumProbeInserts;
if(pNumTargInsertBases != NULL)
	*pNumTargInsertBases = m_NumTargInserts;
return(m_NumBasesAligned + m_NumProbeInserts + m_NumTargInserts);
}

//
// Helper functions for Band cell dereferencing

tNWTrcBckCell *						// returns ptr to newly allocated cell  or NULL if errors
CNeedlemanWunsch::AllocBandCell(UINT32 ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn)		// current target base position 1..m_TargLen 
{
tsNWColBand *pNWColBand;
tNWTrcBckCell *pTrcBckCell;

if(ProbeBasePsn == 0 || ProbeBasePsn > m_ProbeLen || TargBasePsn == 0 || TargBasePsn > m_TargLen)
	return(NULL);
if((ProbeBasePsn - 1) > m_ColBandsUsed || ProbeBasePsn < m_ColBandsUsed)
	return(NULL);

if(m_TrcBckCellsUsed == m_TrcBckCellsAllocd)			// need to allocate more?
	{
	size_t ReallocSize;
	UINT32 EstReqCells;
	tNWTrcBckCell *pRealloc;
	EstReqCells = 10000 + (UINT32)(((double)m_ProbeLen / ProbeBasePsn) * m_TrcBckCellsAllocd); 
	ReallocSize = sizeof(tNWTrcBckCell) * EstReqCells;
#ifdef _WIN32
	pRealloc = (tNWTrcBckCell *)realloc(m_pTrcBckCells,ReallocSize);
#else
	pRealloc = (tNWTrcBckCell *)mremap(m_pTrcBckCells,m_TrcBckCellsAllocdSize,ReallocSize,MREMAP_MAYMOVE);
	if(pRealloc == MAP_FAILED)
		pRealloc = NULL;
#endif
	if(pRealloc == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocBandCell: traceback cell memory re-allocation to %lld bytes - %s",(INT64)ReallocSize,strerror(errno));
		return(NULL);
		}
	m_TrcBckCellsAllocdSize = ReallocSize;
	m_TrcBckCellsAllocd = EstReqCells;
	m_pTrcBckCells = pRealloc;
	}

pNWColBand = &m_pColBands[ProbeBasePsn - 1];
pTrcBckCell = &m_pTrcBckCells[m_TrcBckCellsUsed++];
*pTrcBckCell = 0;
if(ProbeBasePsn > m_ColBandsUsed)
	{
	m_ColBandsUsed += 1;
	pNWColBand->StartTargBasePsn = TargBasePsn;
	pNWColBand->EndTargBasePsn = TargBasePsn;
	pNWColBand->TrcBckCellPsn = m_TrcBckCellsUsed;
	return(pTrcBckCell);
	}
if(pNWColBand->EndTargBasePsn != TargBasePsn - 1)
	{
	m_TrcBckCellsUsed -= 1;
	return(NULL);
	}
pNWColBand->EndTargBasePsn = TargBasePsn;
return(pTrcBckCell);
}

tNWTrcBckCell *					// returns ptr to cell or NULL if errors
CNeedlemanWunsch::DerefBandCell(UINT32 ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn)		// current target base position 1..m_TargLen
{
UINT32 BandIdx;
tsNWColBand *pNWColBand;
if(m_ColBandsUsed == 0 || ProbeBasePsn == 0 || ProbeBasePsn > m_ColBandsUsed || TargBasePsn == 0 || TargBasePsn > m_TargLen)
	return(NULL);
pNWColBand = &m_pColBands[ProbeBasePsn - 1];
if(TargBasePsn < pNWColBand->StartTargBasePsn || TargBasePsn > pNWColBand->EndTargBasePsn || pNWColBand->TrcBckCellPsn == 0)
	return(NULL);
BandIdx = pNWColBand->TrcBckCellPsn - 1 + (TargBasePsn - pNWColBand->StartTargBasePsn);
return(&m_pTrcBckCells[BandIdx]);
}

tNWTrcBckCell *					// returns ptr to cell  or NULL if errors
CNeedlemanWunsch::DerefBandCellLeft(UINT32 ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn)		// current target base position 1..m_TargLen
{
UINT32 BandIdx;
tsNWColBand *pNWColBand;
if(m_ColBandsUsed == 0 || ProbeBasePsn < 2 || ProbeBasePsn > m_ColBandsUsed || TargBasePsn == 0 || TargBasePsn > m_TargLen)
	return(NULL);
pNWColBand = &m_pColBands[ProbeBasePsn - 2];
if(TargBasePsn < pNWColBand->StartTargBasePsn || TargBasePsn > pNWColBand->EndTargBasePsn || pNWColBand->TrcBckCellPsn == 0)
	return(NULL);
BandIdx = pNWColBand->TrcBckCellPsn - 1 + (TargBasePsn - pNWColBand->StartTargBasePsn);
return(&m_pTrcBckCells[BandIdx]);
}

tNWTrcBckCell *					// returns ptr to cell  or NULL if errors
CNeedlemanWunsch::DerefBandCellDiag(UINT32 ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn)		// current target base position 1..m_TargLen
{
UINT32 BandIdx;
tsNWColBand *pNWColBand;
if(m_ColBandsUsed == 0 || ProbeBasePsn < 2 || ProbeBasePsn > m_ColBandsUsed || TargBasePsn < 2 || TargBasePsn > m_TargLen)
	return(NULL);
TargBasePsn -= 1;
pNWColBand = &m_pColBands[ProbeBasePsn - 2];
if(TargBasePsn < pNWColBand->StartTargBasePsn || TargBasePsn > pNWColBand->EndTargBasePsn || pNWColBand->TrcBckCellPsn == 0)
	return(NULL);
BandIdx = pNWColBand->TrcBckCellPsn - 1 + (TargBasePsn - pNWColBand->StartTargBasePsn);
return(&m_pTrcBckCells[BandIdx]);
}

tNWTrcBckCell *					// returns ptr to cell  or NULL if errors
CNeedlemanWunsch::DerefBandCellDown(UINT32 ProbeBasePsn,	// current probe base position 1..m_ProbeLen
							  UINT32 TargBasePsn)		// current target base position 1..m_TargLen
{
UINT32 BandIdx;
tsNWColBand *pNWColBand;
if(m_ColBandsUsed == 0 || ProbeBasePsn == 0 || ProbeBasePsn > m_ColBandsUsed || TargBasePsn == 0 || TargBasePsn > m_TargLen)
	return(NULL);
pNWColBand = &m_pColBands[ProbeBasePsn - 1];
if(TargBasePsn < pNWColBand->StartTargBasePsn || TargBasePsn > pNWColBand->EndTargBasePsn || pNWColBand->TrcBckCellPsn == 0)
	return(NULL);
BandIdx = pNWColBand->TrcBckCellPsn - 1 +  (TargBasePsn - pNWColBand->StartTargBasePsn);
return(&m_pTrcBckCells[BandIdx - 1]);
}


int
CNeedlemanWunsch::DumpScores(char *pszFile,		// dump Smith-Waterman matrix to this csv file
					char Down,	// use this char to represent cell down link representing base inserted into target relative to probe
					char Left,	// use this char to represent cell left link representing base inserted into probe relative to target
				    char Diag)	// use this char to represent cell diagonal  representing matching base either exact or mismatch
{
char *pszBuff;
char *pRow;
UINT32 BuffIdx;
UINT32 TargIdx;
UINT32 ProbeIdx;
tNWTrcBckCell TrcBckDir;
int TrcBckScore;
int hDumpFile;
tNWTrcBckCell *pCell;

UINT32 EstDumpLen;

if(m_bBanded)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Banded DumpScores not currently supported");
	return(eBSFerrParams);
	}

if(!m_bAligned || pszFile == NULL || pszFile[0] == '\0')
	return(eBSFerrParams);

EstDumpLen = (UINT32)min((UINT64)m_TargLen * m_ProbeLen * 10,(UINT64)0x03ffff);
if((pszBuff = new char [EstDumpLen])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate dump buffer memory %d bytes",EstDumpLen);
	return(eBSFerrMem);
	}

#ifdef _WIN32
if((hDumpFile = open(pszFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hDumpFile = open(pszFile,O_RDWR | O_CREAT | O_TRUNC,S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszFile,strerror(errno));
	return(eBSFerrCreateFile);
	}

// target bases along top as columns, probe bases vertically as rows
pRow = pszBuff;
for(BuffIdx = TargIdx = 0; TargIdx < m_TargLen; TargIdx++)
	{
	BuffIdx += sprintf(&pszBuff[BuffIdx],",\" \",\"%c\"", CSeqTrans::MapBase2Ascii(m_pTarg[TargIdx]));
	if(BuffIdx + 1000 > EstDumpLen)
		{
		CUtility::SafeWrite(hDumpFile,pszBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
pCell = m_pTrcBckCells;
for(ProbeIdx = 0; ProbeIdx < m_ProbeLen; ProbeIdx++)
	{
	BuffIdx += sprintf(&pszBuff[BuffIdx],"\n%c",CSeqTrans::MapBase2Ascii(m_pProbe[ProbeIdx]));
	for(TargIdx = 0; TargIdx < m_TargLen; TargIdx++,pCell++)
		{
		TrcBckDir = *pCell & cNWTrcBckMsk;
		TrcBckScore = *pCell & cNWScoreMsk;
		if(TrcBckScore & cNWScoreNegFlg)			// if sign of score is negative then sign extend to left making score into negative int
				TrcBckScore |= ~cNWScoreMsk;

		if(TrcBckDir == cNWTrcBckDiagFlg)
			BuffIdx += sprintf(&pszBuff[BuffIdx],",\"%c\",%d",Diag,TrcBckScore);
		else 
			{
			if(TrcBckDir == cNWTrcBckDownFlg)
				BuffIdx += sprintf(&pszBuff[BuffIdx],",\"%c\",%d",Down,TrcBckScore);
			else
				{
				if(TrcBckDir == cNWTrcBckLeftFlg)
					BuffIdx += sprintf(&pszBuff[BuffIdx],",\"%c\",%d",Left,TrcBckScore);
				else
					BuffIdx += sprintf(&pszBuff[BuffIdx],",\" \",%d",TrcBckScore);
				}
			}
		if(BuffIdx + 1000 > EstDumpLen)
			{
			CUtility::SafeWrite(hDumpFile,pszBuff,BuffIdx);
			BuffIdx = 0;
			}
		}
	}
if(BuffIdx)
	CUtility::SafeWrite(hDumpFile,pszBuff,BuffIdx);
#ifdef _WIN32
_commit(hDumpFile);
#else
fsync(hDumpFile);
#endif
close(hDumpFile);
return(eBSFSuccess);
}



