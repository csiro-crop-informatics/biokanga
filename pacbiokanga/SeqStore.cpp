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

m_NumStoredSeqs = 0;
m_TotStoredSeqsLen = 0;
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
		return(false);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pSeqHdrs = (tsSeqHdr *)mmap(NULL,m_AllocdSeqHdrSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSeqHdrs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for sequence headers",(INT64)m_AllocdSeqHdrSize);
		m_pSeqHdrs = NULL;
		m_AllocdSeqHdrSize = 0;
		return(false);
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
		return(false);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pDescrSeqs = (UINT8 *)mmap(NULL,m_AllocdDescrSeqsSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pDescrSeqs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to allocate %lld bytes contiguous memory for sequence headers",(INT64)m_AllocdDescrSeqsSize);
		m_pDescrSeqs = NULL;
		m_AllocdDescrSeqsSize = 0;
		return(false);
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
return(m_NumStoredSeqs);	
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
UINT32 Idx;
tsSeqHdr *pSeqHdr;
UINT8 *pDescrSeq;
UINT8 DescrLen = (UINT8)strlen(pszDescr);

if(m_pSeqHdrs == NULL || m_pDescrSeqs == NULL || m_NumStoredSeqs == 0)
	return(0);
pSeqHdr = m_pSeqHdrs;
for(Idx = 0; Idx < m_NumStoredSeqs; Idx++,pSeqHdr++)
	{
	if(pSeqHdr->DescrLen != DescrLen)
		continue;
	pDescrSeq = &m_pDescrSeqs[pSeqHdr->DescrSeqOfs];
	if(!stricmp((const char *)pDescrSeq,pszDescr))
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
strncpy(pszDescr,(const char *)pDescrSeq,MaxDescrLen);
pszDescr[MaxDescrLen] = '\0';
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


