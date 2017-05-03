/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */
#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include <process.h>
#include "../libbiokanga/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libbiokanga/commhdrs.h"
#endif

#include "pacbiokanga.h"
#include "SSW.h"
#include "MAConsensus.h"


CMAConsensus::CMAConsensus()
{
m_pMACols = NULL;  
m_pRefSeqs = NULL;
Reset();
}


CMAConsensus::~CMAConsensus()
{
if(m_pMACols != NULL)
	{
#ifdef _WIN32
	free(m_pMACols);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pMACols != MAP_FAILED)
		munmap(m_pMACols,m_AllocMAColsSize);
#endif
	}

if(m_pRefSeqs != NULL)
	delete m_pRefSeqs;

}

void 
CMAConsensus::Reset(void)					// reset state back to that immediately following instantiation
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
	}

if(m_pRefSeqs != NULL)
	{
	delete m_pRefSeqs;
	m_pRefSeqs = NULL;
	}

m_NumRefSeqs = 0;
m_AllocdRefSeqs = 0;
m_MATotRefSeqLen = 0;
m_MACurTotRefSeqLen = 0;
m_MACurCols = 0;				
m_AllocMACols = 0;			
m_AllocMAColsSize = 0;	
m_bStartedMultiAlignments = false;
}


int													// eBSFSuccess or error code 
CMAConsensus::Init(UINT32 NumRefSeqs,				// max number of reference sequences which will be added
				UINT64 TotRefSeqLen)				// max total bases of all reference sequences which will be added
{
size_t memreq;
if(NumRefSeqs < 1 || NumRefSeqs > cMaxNumRefSeqs || TotRefSeqLen < (NumRefSeqs * cMinRefSeqLen) || TotRefSeqLen > cMaxTotRefSeqLens)
	return(eBSFerrParams);
Reset();
m_MATotRefSeqLen = TotRefSeqLen;
memreq = (sizeof(tsMAlignConCol) * TotRefSeqLen * 110)/100; // allocate to hold at least SeqLen columns with Depth bases and a 10% overallocation to reduce chances of a reallocation later on as sequence insertions are discovered
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
	m_pMACols = (tsMAlignConCol *) malloc(memreq);
	if(m_pMACols == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for multialignment",(INT64)memreq);
		return(eBSFerrMem);
		}
#else
		// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pMACols = (tsMAlignConCol *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pMACols == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for multialignment",(INT64)memreq);
		m_pMACols = NULL;
		return(eBSFerrMem);
		}
#endif
	m_AllocMAColsSize = memreq;	
	}

m_AllocMACols = (UINT64)(m_AllocMAColsSize/sizeof(tsMAlignConCol)); 
memset(m_pMACols,0,m_AllocMAColsSize);
m_MACurCols = 0;

m_AllocdRefSeqs = NumRefSeqs;
if((m_pRefSeqs = new tsMARefSeq [NumRefSeqs]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %d bytes memory for %d reference sequences",(int)(sizeof(tsMARefSeq) * m_AllocdRefSeqs),NumRefSeqs);
	return(eBSFerrMem);
	}

m_NumRefSeqs = 0;
memset(m_pRefSeqs,0,sizeof(tsMARefSeq) * m_AllocdRefSeqs);
return(eBSFSuccess);
}

UINT32												// returned reference sequence indentifier to be used when adding alignments to this reference sequence
CMAConsensus::AddRefSeq(UINT32 SeqLen,				// reference sequence to which other sequences are aligned is this length
					etSeqBase *pRefSeq)			// reference sequence
{
UINT32 Idx;
UINT32 Idy;
tsMARefSeq *psRefSeq;
UINT8 *pBase;
tsMAlignConCol *pCol;

if(pRefSeq == NULL || SeqLen < cMinRefSeqLen || SeqLen > cMaxRefSeqLen || (SeqLen +  m_MACurTotRefSeqLen) > m_MATotRefSeqLen)
	return(0);

if(m_NumRefSeqs == m_AllocdRefSeqs)
	return(0);

psRefSeq = &m_pRefSeqs[m_NumRefSeqs++];
psRefSeq->RefID = m_NumRefSeqs;
psRefSeq->StartColIdx = m_MACurCols+1;
m_MACurCols += SeqLen;
psRefSeq->EndColIdx = m_MACurCols;
psRefSeq->SeqLen = SeqLen;

// initialise with the reference sequence bases and column depth as 1
pBase = pRefSeq; 
pCol = (tsMAlignConCol *)&m_pMACols[psRefSeq->StartColIdx-1];
for(Idx = 1; Idx <= SeqLen; Idx++, pBase++,pCol++)
	{
	pCol->RefID = m_NumRefSeqs;
	pCol->Ofs = Idx;
	pCol->ColIdx = psRefSeq->StartColIdx + Idx - 1;
	pCol->PrevColIdx = Idx == 1 ? 0 : pCol->ColIdx - 1;
	pCol->NxtColIdx = Idx == SeqLen ? 0 : pCol->ColIdx + 1;
	pCol->Depth = 1;
	pCol->Extn = 0;
	pCol->ConsBase = *pBase & 0x07;
	for(Idy = 0; Idy <= eBaseInDel; Idy++)
		pCol->BaseCnts[Idy] = 0;
	pCol->BaseCnts[pCol->ConsBase] = 1;
	}

m_bStartedMultiAlignments = true;
return(m_NumRefSeqs);
}



int
CMAConsensus::AddMultiAlignment(UINT32 RefSeqID,	// alignment is against this sequence
					  UINT32 RefStartOfs,			// alignment starts at this reference sequence offset (1..n)
					  UINT32 RefEndOfs,				// alignment ends at this reference sequence offset inclusive
					  UINT32 ProbeStartOfs,			// alignment starts at this probe sequence offset (1..n)
					  UINT32 ProbeEndOfs,			// alignment ends at this probe sequence offset inclusive
					  etSeqBase *pProbeSeq,			// alignment probe sequence
					  UINT32 NumMAAlignOps,			// number of alignment operators
					   tMAOp *pMAAlignOps)			// alignment operators
{
tsMARefSeq *pRefSeq;
tMAOp Op;
tMAOp *pAlignOps;
UINT32 NumAlignOps;
int Idy;
tsMAlignConCol *pNewCol;
tsMAlignConCol *pCol;
tsMAlignConCol *pPrevCol;
etSeqBase *pBase;

if(RefSeqID == 0 || RefSeqID > m_NumRefSeqs)
	return(-1);
pRefSeq = &m_pRefSeqs[RefSeqID-1];
pAlignOps = pMAAlignOps;
NumAlignOps = NumMAAlignOps;

// pt to reference column at which alignment starts and probe starting base
pCol = &m_pMACols[pRefSeq->StartColIdx + RefStartOfs - 2];
pBase = &pProbeSeq[ProbeStartOfs-1];

while(NumAlignOps && pCol != NULL)
	{
	Op=*pAlignOps++;
	NumAlignOps -= 1;
	if(Op < cMADelete)				 // if not an InDel then skip over any reference InDel using the reference base as the probes base until a reference match
		{
		while(pCol->Extn != 0)
			{
			pCol->BaseCnts[pCol->ConsBase] += 1;
			pCol->Depth += 1;
			if(pCol->NxtColIdx == 0)  
				return(-1);
			pCol = (tsMAlignConCol *)&m_pMACols[pCol->NxtColIdx - 1];
			}	
		}

	switch(Op) {
		case cMAMatch:		 // base match between probe and reference; note that may not be an exact match
			pCol->BaseCnts[0x07 & *pBase++] += 1;
			break;

		case cMAInsert:				// base inserted into probe relative to reference - or could be base deleted from reference relative to probe
			if(pCol->Extn == 0)		// if not positioned on a reference deletion column then will need to create one
				{
				if((pNewCol = InsertCol(pCol->PrevColIdx)) == NULL)
					return(-1);
				pPrevCol = (tsMAlignConCol *)&m_pMACols[pNewCol->PrevColIdx-1];
				pNewCol->RefID = RefSeqID;
				pNewCol->Depth = pPrevCol->Depth - 1;
				pNewCol->Ofs = pPrevCol->Ofs;
				pNewCol->Extn = pPrevCol->Extn + 1;
				pNewCol->ConsBase = eBaseInDel;
				for(Idy = 0; Idy <= eBaseInDel; Idy++)
					pNewCol->BaseCnts[Idy] = 0;
				pNewCol->BaseCnts[eBaseInDel] = 1;
				pCol = pNewCol;
				}
			pCol->BaseCnts[*pBase++] += 1;
			break;

		case cMADelete:      // base deleted from probe relative to reference - or could be base inserted into reference relative to probe
			pCol->BaseCnts[eBaseInDel] += 1;
			break;
	   }

	pCol->Depth += 1;
	if(pCol->NxtColIdx)   // if not last
		pCol = (tsMAlignConCol *)&m_pMACols[pCol->NxtColIdx-1];
	else
		pCol = NULL;
	}

return(eBSFSuccess);
}

// generate multiple alignment consensus from multialignment columns at m_pMACols
// writes consensus base back into m_pMACols
int
CMAConsensus::GenMultialignConcensus(void)
{
int BaseCnts[eBaseInDel+1];
int AbundIdxs[eBaseInDel+1];
int *pCnts;
int *pNxtCnts;
int XAbundIdxs;
int TotBaseCnts;
UINT32 *pBaseCnts;
int BaseIdx;
tsMAlignConCol *pCol;

if(!m_bStartedMultiAlignments)
	return(0);

pCol = (tsMAlignConCol *)m_pMACols;
do
	{
	memset(BaseCnts,0,sizeof(BaseCnts));
	TotBaseCnts = 0;
	pBaseCnts = pCol->BaseCnts;
	for(BaseIdx = 0; BaseIdx <= eBaseInDel; BaseIdx++,pBaseCnts++)
		{
		BaseCnts[BaseIdx] = *pBaseCnts;
		TotBaseCnts += *pBaseCnts;
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

	// simply choosing the most abundant base (could be equally abundant!) as being the consensus
	pCol->ConsBase = AbundIdxs[0] < eBaseInDel ? AbundIdxs[0] : eBaseUndef;
	if(pCol->NxtColIdx == 0)
		pCol = NULL;
	else
		pCol = &m_pMACols[pCol->NxtColIdx - 1];
	}
while(pCol != NULL);
return(eBSFSuccess);
}

tsMAlignConCol *										// inserted column or NULL if errors inserting
CMAConsensus::InsertCol(UINT64 PrevCol)				// allocate and insert new column after Prev
{
tsMAlignConCol *pCol;
tsMAlignConCol *pPrevCol;
tsMAlignConCol *pNextCol;

if(m_MACurCols + 10 >= m_AllocMACols)			// need to extend previously allocated columns to hold this new column?
	{
	size_t memreq;
	void *pAllocd;
	UINT64 AllocMACols;
	AllocMACols = (m_AllocMACols * 11) / 10;   // increase allocated cols by 10%
	memreq = (size_t)AllocMACols * sizeof(tsMAlignConCol);

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
	m_pMACols = (tsMAlignConCol *)pAllocd;
	memset(&m_pMACols[m_AllocMAColsSize],0,memreq - m_AllocMAColsSize);
	m_AllocMACols = AllocMACols;
	m_AllocMAColsSize = memreq;
	}
pCol = &m_pMACols[m_MACurCols++];
pCol->ColIdx = m_MACurCols;
pPrevCol = &m_pMACols[PrevCol-1];
if(pPrevCol->NxtColIdx != 0)
	{
	pNextCol = (tsMAlignConCol *)&m_pMACols[pPrevCol->NxtColIdx-1];
	pNextCol->PrevColIdx = pCol->ColIdx;
	}
pCol->NxtColIdx = pPrevCol->NxtColIdx;
pCol->PrevColIdx = pPrevCol->ColIdx;
pPrevCol->NxtColIdx = pCol->ColIdx;
return(pCol);
}


int												// number of bases in consensus 
CMAConsensus::GetConsensus(UINT32 RefSeqID)	// for this reference sequence identifier
{
int ConsensusLen;
tsMARefSeq *pRefSeq;
tsMAlignConCol *pCol;
if(RefSeqID == 0 || RefSeqID > m_NumRefSeqs)
	return(-1);
pRefSeq = &m_pRefSeqs[RefSeqID-1];

// pt to reference starting column
pCol = &m_pMACols[pRefSeq->StartColIdx - 1];
ConsensusLen = 0;
while(pCol != NULL)
	{
	if(pCol->ConsBase < eBaseInDel)
		ConsensusLen += 1;
	if(pCol->NxtColIdx == 0)
		break;
	pCol = &m_pMACols[pCol->NxtColIdx - 1];
	}

return(ConsensusLen);
}

int												// number of bases returned in pRetBases 
CMAConsensus::GetConsensus(UINT32 RefSeqID,		// consensus sequence for this reference identifier
					    UINT32 RefStartOfs,		// starting offset 1..N
						UINT32 RefLen,			// return at most this many bases
						etSeqBase *pRetBases)	// where to return bases
{
int ConsensusLen;
tsMARefSeq *pRefSeq;
tsMAlignConCol *pCol;
if(RefSeqID == 0 || RefSeqID > m_NumRefSeqs)
	return(-1);
pRefSeq = &m_pRefSeqs[RefSeqID-1];

// pt to reference starting column
pCol = &m_pMACols[pRefSeq->StartColIdx - 1];
ConsensusLen = 0;
while(pCol != NULL)
	{
	if(pCol->ConsBase <= eBaseN)
		{
		if(RefStartOfs == 1)
			{
			ConsensusLen += 1;
			*pRetBases++ = pCol->ConsBase;
			}
		else
			RefStartOfs -= 1;
		}
	if(pCol->NxtColIdx == 0)
		break;
	pCol = &m_pMACols[pCol->NxtColIdx - 1];
	}
return(ConsensusLen);
}

