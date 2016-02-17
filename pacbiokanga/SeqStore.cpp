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
#include "SeqStore.h"


CSeqStore::CSeqStore()
{
m_pDescrSeqs = NULL;
m_pSeqHdrs = NULL;
m_pSeqDescrIdx = NULL;
Reset();
}


CSeqStore::~CSeqStore()
{
if(m_pSeqHdrs != NULL)
	{
#ifdef _WIN32
	free(m_pSeqHdrs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqHdrs != MAP_FAILED)
		munmap(m_pSeqHdrs,m_AllocdSeqHdrSize);
#endif
	}

if(m_pDescrSeqs)
	{
#ifdef _WIN32
	free(m_pDescrSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pDescrSeqs != MAP_FAILED)
		munmap(m_pDescrSeqs,m_AllocdDescrSeqsSize);
#endif
	}

if(m_pSeqDescrIdx)
	{
#ifdef _WIN32
	free(m_pSeqDescrIdx);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqDescrIdx != MAP_FAILED)
		munmap(m_pSeqDescrIdx,m_AllocdSeqDescrIdxSize);
#endif
	}

}

int
CSeqStore::Reset(void)
{
if(m_pSeqHdrs != NULL)
	{
#ifdef _WIN32
	free(m_pSeqHdrs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqHdrs != MAP_FAILED)
		munmap(m_pSeqHdrs,m_AllocdSeqHdrSize);
#endif
	m_pSeqHdrs = NULL;
	}

if(m_pDescrSeqs)
	{
#ifdef _WIN32
	free(m_pDescrSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pDescrSeqs != MAP_FAILED)
		munmap(m_pDescrSeqs,m_AllocdDescrSeqsSize);
#endif
	m_pDescrSeqs = NULL;
	}


if(m_pSeqDescrIdx)
	{
#ifdef _WIN32
	free(m_pSeqDescrIdx);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqDescrIdx != MAP_FAILED)
		munmap(m_pSeqDescrIdx,m_AllocdSeqDescrIdxSize);
#endif
	m_pSeqDescrIdx = NULL;
	}

m_AllocdSeqHdrSize = 0;
m_AllocdDescrSeqsSize = 0;
m_AllocdSeqDescrIdxSize = 0;
m_NumSeqDescrIdxd = 0;
m_NumStoredSeqs = 0;
m_TotStoredSeqsLen = 0;

m_MinSeqLen = 0;
m_MaxSeqLen = 0;
m_UsedSeqHdrMem = 0; 
m_AllocdSeqHdrSize = 0;
m_UsedDescrSeqsMem = 0; 
m_AllocdDescrSeqsSize = 0;

return(eBSFSuccess);
}


tSeqID									// identifier by which this sequence can later be retreived (0 if unable to add sequence)
CSeqStore::AddSeq(UINT32 Flags,			// any flags associated with this sequence
			   char *pszDescr,			// sequence descriptor
			   UINT32 SeqLen,			// sequence is this length
			   etSeqBase *pSeq)			// sequence to add
{
int DescrLen;
size_t memreq;
void *pAllocd;
tsSeqHdr *pSeqHdr;
UINT8 *pDescrSeq;
char szDescr[cMaxDescrIDLen+1];

if(pSeq == NULL || SeqLen < cMinSeqStoreLen || SeqLen > cMaxSeqStoreLen || m_NumStoredSeqs == cMaxStoredSeqs)
	return(0);


if(m_pSeqDescrIdx)						// force index to be regenerated as adding any new sequences will invalidate index!
	{
#ifdef _WIN32
	free(m_pSeqDescrIdx);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqDescrIdx != MAP_FAILED)
		munmap(m_pSeqDescrIdx,m_AllocdSeqDescrIdxSize);
#endif
	m_pSeqDescrIdx = NULL;
	m_AllocdSeqDescrIdxSize = 0;
	m_NumSeqDescrIdxd = 0;
	}

if(pszDescr == NULL || pszDescr[0] == '\0')
	sprintf(szDescr,"Seq%u",m_NumStoredSeqs+1);
else
	{
	strncpy(szDescr,pszDescr,sizeof(szDescr)-1);
	szDescr[sizeof(szDescr)-1] = '\0';
	}
DescrLen = (int)strlen(szDescr);

memreq = sizeof(tsSeqHdr);
if(m_pSeqHdrs == NULL)
	{
	m_UsedSeqHdrMem = 0;
	m_NumStoredSeqs = 0;
	m_TotStoredSeqsLen = 0;
	m_AllocdSeqHdrSize = cAllocNumSeqHdrs * sizeof(tsSeqHdr);
#ifdef _WIN32
	m_pSeqHdrs = (tsSeqHdr *) malloc(m_AllocdSeqHdrSize);
	if(m_pSeqHdrs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for sequence headers",(INT64)m_AllocdSeqHdrSize);
		return(0);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pSeqHdrs = (tsSeqHdr *)mmap(NULL,m_AllocdSeqHdrSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSeqHdrs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for sequence headers",(INT64)m_AllocdSeqHdrSize);
		m_pSeqHdrs = NULL;
		m_AllocdSeqHdrSize = 0;
		return(0);
		}
#endif
	}
else
	{
	if(m_AllocdSeqHdrSize < m_UsedSeqHdrMem + (memreq * 10))
		{
		memreq = m_AllocdSeqHdrSize + sizeof(tsSeqHdr) * (size_t)cAllocNumSeqHdrs;
#ifdef _WIN32
		pAllocd = realloc(m_pSeqHdrs,memreq);
#else
		pAllocd = mremap(m_pSeqHdrs,m_AllocdSeqHdrSize,memreq,MREMAP_MAYMOVE);
		if(pAllocd == MAP_FAILED)
			pAllocd = NULL;
#endif
		if(pAllocd == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeq: Memory re-allocation to %lld bytes - %s",memreq,strerror(errno));
			return(0);
			}
		m_AllocdSeqHdrSize = memreq;
		m_pSeqHdrs = (tsSeqHdr *)pAllocd;
		}
	}

memreq = SeqLen + DescrLen + 1;
if(m_pDescrSeqs == NULL)
	{
	m_TotStoredSeqsLen = 0;
	m_UsedDescrSeqsMem = 0;
	m_AllocdDescrSeqsSize = max(memreq+1,cAllocDescrSeqMem);
#ifdef _WIN32
	m_pDescrSeqs = (UINT8 *) malloc(m_AllocdDescrSeqsSize);
	if(m_pDescrSeqs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for sequences",(INT64)m_AllocdDescrSeqsSize);
		return(0);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pDescrSeqs = (UINT8 *)mmap(NULL,m_AllocdDescrSeqsSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pDescrSeqs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for sequence headers",(INT64)m_AllocdDescrSeqsSize);
		m_pDescrSeqs = NULL;
		m_AllocdDescrSeqsSize = 0;
		return(0);
		}
#endif
	}
else
	{
	if(m_AllocdDescrSeqsSize <= m_UsedDescrSeqsMem + memreq + 10)
		{
		memreq = m_AllocdDescrSeqsSize + max(memreq+1,cAllocDescrSeqMem);
#ifdef _WIN32
		pAllocd = realloc(m_pDescrSeqs,memreq);
#else
		pAllocd = mremap(m_pDescrSeqs,m_UsedDescrSeqsMem,memreq,MREMAP_MAYMOVE);
		if(pAllocd == MAP_FAILED)
			pAllocd = NULL;
#endif
		if(pAllocd == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeq: Memory re-allocation to %lld bytes - %s",memreq,strerror(errno));
			return(0);
			}
		m_AllocdDescrSeqsSize = memreq;
		m_pDescrSeqs = (UINT8 *)pAllocd;
		}
	}

m_UsedSeqHdrMem += sizeof(tsSeqHdr);
pSeqHdr = &m_pSeqHdrs[m_NumStoredSeqs++];
pSeqHdr->SeqID = m_NumStoredSeqs;
pSeqHdr->Flags = Flags;
pSeqHdr->SeqLen = SeqLen;
pSeqHdr->DescrLen = (UINT8)DescrLen;
pSeqHdr->DescrSeqOfs = m_UsedDescrSeqsMem;
pDescrSeq = &m_pDescrSeqs[m_UsedDescrSeqsMem]; 
memcpy(pDescrSeq,pszDescr,DescrLen);
pDescrSeq[DescrLen] = '\0';
pDescrSeq += DescrLen + 1;
m_UsedDescrSeqsMem += DescrLen + 1;
memcpy(pDescrSeq,pSeq,SeqLen);
m_UsedDescrSeqsMem += SeqLen;
m_TotStoredSeqsLen += SeqLen;
if(SeqLen > m_MaxSeqLen)
	m_MaxSeqLen = SeqLen;
if(m_MinSeqLen == 0 || m_MinSeqLen > SeqLen)
	m_MinSeqLen = SeqLen;		
return(m_NumStoredSeqs);	
}

static CSeqStore *pThis;        // qsort references this instance
int
CSeqStore::GenSeqDescrIdx(void) // generate index over all sequence descriptors, index used when retrieving sequence identifier for given descriptor by GetSeqID()
{
UINT32 *pSeqDescrIdx;
UINT32 Idx;

if(m_NumStoredSeqs == 0)		// can't index if no sequence descriptors!
	return(eBSFerrFastaDescr);

m_NumSeqDescrIdxd = 0;
m_AllocdSeqDescrIdxSize = sizeof(UINT32) * (m_NumStoredSeqs + 1000);	// allocating more than actually required to reduce potential for future realloc's 
#ifdef _WIN32
m_pSeqDescrIdx = (UINT32 *) malloc(m_AllocdSeqDescrIdxSize);
if(m_pSeqDescrIdx == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for sequences",(INT64)m_AllocdSeqDescrIdxSize);
	m_AllocdSeqDescrIdxSize = 0;
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
m_pSeqDescrIdx = (UINT32 *)mmap(NULL,m_AllocdSeqDescrIdxSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pSeqDescrIdx == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for sequence headers",(INT64)m_AllocdSeqDescrIdxSize);
	m_pSeqDescrIdx = NULL;
	m_AllocdSeqDescrIdxSize = 0;
	return(eBSFerrMem);
	}
#endif
m_NumSeqDescrIdxd = m_NumStoredSeqs;
pSeqDescrIdx = m_pSeqDescrIdx;
for(Idx = 1; Idx <= m_NumSeqDescrIdxd; Idx++, pSeqDescrIdx += 1)
	*pSeqDescrIdx = Idx;
m_MTqsort.SetMaxThreads(cDfltSortThreads);
pThis = this;
m_MTqsort.qsort(m_pSeqDescrIdx,(INT64)m_NumSeqDescrIdxd,sizeof(UINT32),SortSeqDescr);
return(eBSFSuccess);
}

int
CSeqStore::SortSeqDescr(const void *arg1, const void *arg2)
{
UINT32 Idx1 = *(UINT32 *)arg1;
UINT32 Idx2 = *(UINT32 *)arg2;

tsSeqHdr *pHd1;
tsSeqHdr *pHd2;

pHd1 = &pThis->m_pSeqHdrs[Idx1-1];
pHd2 = &pThis->m_pSeqHdrs[Idx2-1];

char *pDescr1 = (char *)&pThis->m_pDescrSeqs[pHd1->DescrSeqOfs];
char *pDescr2 = (char *)&pThis->m_pDescrSeqs[pHd2->DescrSeqOfs];

return(stricmp(pDescr1,pDescr2));
}

UINT32									// previous flags
CSeqStore::SetFlags(tSeqID SeqID,		// set flags associated with sequence identified by SeqID
					UINT32 Flags)		// flags to be associated with this sequence
{
UINT32 PrevFlags;
tsSeqHdr *pSeqHdr;
if(m_pSeqHdrs == NULL || m_pDescrSeqs == NULL || SeqID < 1 || SeqID > m_NumStoredSeqs)
	return(0);
pSeqHdr = &m_pSeqHdrs[SeqID - 1];
PrevFlags = pSeqHdr->Flags;
pSeqHdr->Flags = Flags;
return(PrevFlags);
}

tSeqID								// sequence identifier of sequence matching on same descriptor
CSeqStore::GetSeqID(char *pszDescr)	// must match on this descriptor
{
if(m_pSeqHdrs == NULL || m_pDescrSeqs == NULL || m_NumStoredSeqs == 0)
	return(0);

if(m_pSeqDescrIdx == NULL || m_NumSeqDescrIdxd == 0)	// generate index if not already generated
	if(GenSeqDescrIdx() != eBSFSuccess)
		return(0);

return(LocateSeqID(pszDescr));
}


tSeqID								// sequence identifier of sequence matching on same descriptor
CSeqStore::LocateSeqID(char *pszDescr)	// must match on this descriptor
{
tsSeqHdr *pSeqHdr;
char *pszSeqDescr;
int Rslt;
if(m_pSeqDescrIdx == NULL || m_NumSeqDescrIdxd == 0 || pszDescr == NULL || *pszDescr == '\0')
	return(0);
int Lo,Mid,Hi;	// search limits
Lo = 0; Hi = m_NumSeqDescrIdxd-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pSeqHdr = &m_pSeqHdrs[m_pSeqDescrIdx[Mid]-1];
	pszSeqDescr = (char *)&m_pDescrSeqs[pSeqHdr->DescrSeqOfs];
	Rslt = stricmp(pszDescr,pszSeqDescr);
	if(Rslt < 0)	
		{
		Hi = Mid - 1;
		continue;
		}
	if(Rslt > 0)	
		{
		Lo = Mid + 1;
		continue;
		}
	return(pSeqHdr->SeqID);
	}
return(0);
}

UINT32								// returned number of currently stored sequences
CSeqStore::GetNumSeqs(void)			// get number of currently stored sequences
{
if(m_pSeqHdrs == NULL || m_pDescrSeqs == NULL)
	return(0);
return(m_NumStoredSeqs);
}


UINT32								// returned max length of any currently stored sequence
CSeqStore::GetMaxSeqLen(void)		// get max length of any currently stored sequence
{
if(m_pSeqHdrs == NULL || m_pDescrSeqs == NULL)
	return(0);
return(m_MaxSeqLen);
}

UINT32								// returned min length of any currently stored sequence
CSeqStore::GetMinSeqLen(void)		// get min length of any currently stored sequence
{
if(m_pSeqHdrs == NULL || m_pDescrSeqs == NULL)
	return(0);
return(m_MinSeqLen);
}

size_t								// returned total length of all currently stored sequences
CSeqStore::GetTotalLen(void)			// get total length of all currently stored sequences
{
if(m_pSeqHdrs == NULL || m_pDescrSeqs == NULL || m_NumStoredSeqs == 0)
	return(0);
return(m_TotStoredSeqsLen);
}

UINT32									// returned sequence length
CSeqStore::GetLen(tSeqID SeqID)		// get sequence length identified by SeqID
{
if(m_pSeqHdrs == NULL || m_pDescrSeqs == NULL || SeqID < 1 || SeqID > m_NumStoredSeqs)
	return(0);
return(m_pSeqHdrs[SeqID - 1].SeqLen);
}

UINT32									// returned flags
CSeqStore::GetFlags(tSeqID SeqID)		// get flags associated with sequence identified by SeqID
{
if(m_pSeqHdrs == NULL || m_pDescrSeqs == NULL || SeqID < 1 || SeqID > m_NumStoredSeqs)
	return(0);
return(m_pSeqHdrs[SeqID - 1].Flags);
}

UINT32									// returned flags
CSeqStore::GetDescr(tSeqID SeqID,		// get descriptor associated with sequence identified by SeqID
					int MaxDescrLen,	// return at most this many descriptor chars copied into
					char *pszDescr)		// this returned sequence descriptor
{
tsSeqHdr *pSeqHdr;
UINT8 *pDescrSeq;

if(pszDescr == NULL || MaxDescrLen == 0 || m_pSeqHdrs == NULL || m_pDescrSeqs == NULL || SeqID < 1 || SeqID > m_NumStoredSeqs)
	return(0);

pSeqHdr = &m_pSeqHdrs[SeqID - 1];
pDescrSeq = &m_pDescrSeqs[pSeqHdr->DescrSeqOfs];
strncpy(pszDescr,(const char *)pDescrSeq,MaxDescrLen - 1);
pszDescr[MaxDescrLen - 1] = '\0';
return(pSeqHdr->Flags);
}


UINT32									// returned sequence copied into pRetSeq is this length; 0 if errors
CSeqStore::GetSeq(tSeqID SeqID,			// get sequence identified by SeqID
				UINT32 StartOfs,		// starting from this sequence offset
				UINT32 MaxSeqSize,		// return at most this many sequence bases  
				etSeqBase *pRetSeq)		// copy sequence into this caller allocated buffer 
{
tsSeqHdr *pSeqHdr;
UINT8 *pDescrSeq;
UINT32 SeqLen;

if(m_pSeqHdrs == NULL || m_pDescrSeqs == NULL || SeqID < 1 || SeqID > m_NumStoredSeqs)
	return(0);

pSeqHdr = &m_pSeqHdrs[SeqID - 1];
if(StartOfs >= pSeqHdr->SeqLen)
	return(0);
SeqLen = min(MaxSeqSize,pSeqHdr->SeqLen - StartOfs);

pDescrSeq = &m_pDescrSeqs[pSeqHdr->DescrSeqOfs + (UINT64)(pSeqHdr->DescrLen + 1 + StartOfs)];
memcpy(pRetSeq,pDescrSeq,SeqLen);
return(SeqLen);
}


