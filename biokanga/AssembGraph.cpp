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

#include "AssembGraph.h"

static tsGraphOutEdge *pStaticGraphOutEdges = NULL;		// static for sorting, will be initialised with m_pGraphOutEdges immediately prior to sorting
static  tEdgeID *pStaticGraphInEdges = NULL;			// static for sorting, will be initialised with m_pGraphInEdges immediately prior to sorting
static tsComponent *pStaticComponents = NULL;			// static for sorting, will be initialised with m_pComponents immediately prior to sorting
static tsGraphVertex *pStaticGraphVertices = NULL;		// static for sorting, will be initialised with m_pGraphVertices immediately prior to sorting

CAssembGraph::CAssembGraph(void)
{
m_pGraphVertices = NULL;
m_pGraphOutEdges = NULL;
m_pGraphInEdges = NULL;
m_pTransitStack = NULL;
m_pComponents = NULL;
m_bMutexesCreated = false;
Reset();
}

CAssembGraph::~CAssembGraph(void)
{
Reset();
}

int
CAssembGraph::CreateMutexes(void)
{
if(m_bMutexesCreated)
	return(eBSFSuccess);

#ifdef _WIN32
InitializeSRWLock(&m_hRwLock);
#else
if(pthread_rwlock_init (&m_hRwLock,NULL)!=0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create rwlock");
	return(eBSFerrInternal);
	}
#endif

#ifdef _WIN32
if(!InitializeCriticalSectionAndSpinCount(&m_hSCritSect,1000))
	{
#else
if(pthread_spin_init(&m_hSpinLock,PTHREAD_PROCESS_PRIVATE)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

#ifdef _WIN32
if((m_hMtxMHReads = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxMHReads,NULL)!=0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
#ifndef _WIN32
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	return(eBSFerrInternal);
	}

m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CAssembGraph::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
#ifdef _WIN32
CloseHandle(m_hMtxMHReads);

#else
pthread_mutex_destroy(&m_hMtxMHReads);
pthread_rwlock_destroy(&m_hRwLock);

#endif
m_bMutexesCreated = false;
}

void
CAssembGraph::AcquireSerialise(void)
{
int SpinCnt = 1000;
#ifdef _WIN32
while(!TryEnterCriticalSection(&m_hSCritSect))
	{
	if(SpinCnt -= 1)
		continue;
	SwitchToThread();
	SpinCnt = 100;
	}
#else
while(pthread_spin_trylock(&m_hSpinLock)==EBUSY)
	{
	if(SpinCnt -= 1)
		continue;
	pthread_yield();
	SpinCnt = 100;
	}
#endif
}

void
CAssembGraph::ReleaseSerialise(void)
{
#ifdef _WIN32
LeaveCriticalSection(&m_hSCritSect);
#else
pthread_spin_unlock(&m_hSpinLock);
#endif
}

void 
inline CAssembGraph::AcquireLock(bool bExclusive)
{
#ifdef _WIN32
if(bExclusive)
	AcquireSRWLockExclusive(&m_hRwLock);
else
	AcquireSRWLockShared(&m_hRwLock);
#else
if(bExclusive)
	pthread_rwlock_wrlock(&m_hRwLock);
else
	pthread_rwlock_rdlock(&m_hRwLock);
#endif
}

void
inline CAssembGraph::ReleaseLock(bool bExclusive)
{

#ifdef _WIN32
if(bExclusive)
	ReleaseSRWLockExclusive(&m_hRwLock);
else
	ReleaseSRWLockShared(&m_hRwLock);
#else
pthread_rwlock_unlock(&m_hRwLock);
#endif
}

void 
CAssembGraph::Reset(void)	
{
#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
if(m_pGraphVertices != NULL)
	{
#ifdef _WIN32
	free(m_pGraphVertices);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pGraphVertices != MAP_FAILED)
		munmap(m_pGraphVertices,m_AllocGraphVertices * sizeof(tsGraphVertex));
#endif	
	m_pGraphVertices = NULL;
	}

if(m_pGraphOutEdges != NULL)
	{
#ifdef _WIN32
	free(m_pGraphOutEdges);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pGraphOutEdges != MAP_FAILED)
		munmap(m_pGraphOutEdges,m_AllocGraphOutEdges  * sizeof(tsGraphOutEdge));
#endif	
	m_pGraphOutEdges = NULL;
	}

if(m_pGraphInEdges != NULL)
	{
#ifdef _WIN32
	free(m_pGraphInEdges);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pGraphInEdges != MAP_FAILED)
		munmap(m_pGraphInEdges,m_AllocGraphInEdges  * sizeof(tEdgeID));
#endif	
	m_pGraphInEdges = NULL;
	}

if(m_pTransitStack != NULL)
	{
#ifdef _WIN32
	free(m_pTransitStack);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pTransitStack != MAP_FAILED)
		munmap(m_pTransitStack,m_AllocTransitStack  * sizeof(tVertID));
#endif	
	m_pTransitStack = NULL;
	}

if(m_pComponents != NULL)
	{
#ifdef _WIN32
	free(m_pComponents);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pComponents != MAP_FAILED)
		munmap(m_pComponents,m_AllocComponents  * sizeof(tsComponent));
#endif	
	m_pComponents = NULL;
	}


m_AllocGraphVertices = 0;
m_AllocGraphOutEdges = 0;
m_AllocGraphInEdges = 0;

m_AllocTransitStack = 0;
m_UsedGraphVertices = 0;
m_UsedGraphOutEdges = 0;
m_UsedGraphInEdges = 0;

m_CurTransitDepth = 0;
m_MaxTransitDepth = 0;

m_VerticesSortOrder = eVSOUnsorted;
m_bOutEdgeSorted = false;
m_bVertexEdgeSet = false;
m_bInEdgeSorted = false;

m_bReduceEdges = true;
m_NumReducts = 0;

m_UsedComponents = 0;	
m_AllocComponents= 0;

m_NumDiscRemaps = 0;
m_bTerminate = false;
DeleteMutexes();
}

teBSFrsltCodes 
CAssembGraph::Init(void)	
{
Reset();
CreateMutexes();
return(eBSFSuccess);
}


teBSFrsltCodes
CAssembGraph::SetNumThreads(int maxThreads)
{
if(maxThreads < 0 || maxThreads > cMaxWorkerThreads)
		return(eBSFerrParams);
m_NumThreads = maxThreads;
m_MTqsort.SetMaxThreads(maxThreads);
CreateMutexes();
return(eBSFSuccess);
}

teBSFrsltCodes	 
CAssembGraph::AddEdges(int NumSeqs,				// number of overlapped sequences
		tsOverlappedSeq *pOverlappedSeqs) // overlapped sequences
{
int Idx;
tsGraphOutEdge *pOutEdge;
tsOverlappedSeq *pOverlap;
size_t AllocMem;
tSeqID CurOverlapSeqID;
tSeqID CurOverlappedSeqID;
tsGraphVertex *pVertex;
tVertID OverlappingVertexID;
tVertID OverlappedVertexID;

if(NumSeqs == 0 || pOverlappedSeqs == NULL)
	return(eBSFerrParams);
AcquireSerialise();
if(m_bTerminate)		// check if should early exit
	{
	ReleaseSerialise();
	return(eBSFerrInternal);
	}
if((m_UsedGraphOutEdges + (UINT32)NumSeqs) >= cMaxNumEdges)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: too many edges, can only accept at most %u",cMaxNumEdges);
	m_bTerminate = true;
	ReleaseSerialise();
	return(eBSFerrMaxEntries);
	}

if(m_VerticesSortOrder != eVSOSeqID)		// vertices need to have been sorted by SeqID ascending
	{
	if(m_UsedGraphVertices > 1)
		{
		pStaticGraphVertices = m_pGraphVertices;
		m_MTqsort.qsort(pStaticGraphVertices,m_UsedGraphVertices,sizeof(tsGraphVertex),SortVerticesSeqID);
		}
	m_VerticesSortOrder = eVSOSeqID;
	}

// ensure m_pGraphOutEdges allocated to hold an additional NumSeqs
if(m_pGraphOutEdges == NULL)				// initialisation may be required
	{
	AllocMem = cInitialAllocVertices * sizeof(tsGraphOutEdge);
#ifdef _WIN32
	m_pGraphOutEdges = (tsGraphOutEdge *) malloc(AllocMem);	
	if(m_pGraphOutEdges == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: graph forward edge (%d bytes per edge) allocation of %d edges failed - %s",
								(int)sizeof(tsGraphOutEdge),cInitialAllocVertices,strerror(errno));
		m_bTerminate = true;
		ReleaseSerialise();
		return(eBSFerrMem);
		}
#else
	m_pGraphOutEdges = (tsGraphOutEdge *)mmap(NULL,AllocMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pGraphOutEdges == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: graph forward edge (%d bytes per edge) allocation of %d edges failed - %s",
								(int)sizeof(tsGraphOutEdge),cInitialAllocVertices,strerror(errno));
		m_bTerminate = true;
		ReleaseSerialise();
		return(eBSFerrMem);
		}
#endif
	m_AllocGraphOutEdges = cInitialAllocVertices;
	m_UsedGraphOutEdges = 0;
	}
else
	{
	if((m_UsedGraphOutEdges + NumSeqs + 4) >= m_AllocGraphOutEdges) // need to add additional forward graph edges, 4 is simply to allow a small margin of error
		{
		tsGraphOutEdge *pTmp;
		UINT64 AllocEdges;
		UINT64 ReallocEdges;

		ReallocEdges =  (UINT64)((double)m_AllocGraphOutEdges * cReallocEdges);
		AllocEdges = (UINT64)m_AllocGraphOutEdges + ReallocEdges;
		if(AllocEdges > 0x0ffffffff)
			AllocEdges = 0x0ffffffff;
		AllocMem = (size_t)(AllocEdges * sizeof(tsGraphOutEdge));
#ifdef _WIN32
		pTmp = (tsGraphOutEdge *)realloc(m_pGraphOutEdges,AllocMem);
#else
		pTmp = (tsGraphOutEdge *)mremap(m_pGraphOutEdges,m_AllocGraphOutEdges * sizeof(tsGraphOutEdge),AllocMem,MREMAP_MAYMOVE);
		if(pTmp == MAP_FAILED)
			pTmp = NULL;
#endif
		if(pTmp == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: graph forward edge (%d bytes per edge) re-allocation to %lld edges from %lld failed - %s",
														(int)sizeof(tsGraphOutEdge),AllocMem,m_AllocGraphOutEdges,AllocEdges,strerror(errno));
			m_bTerminate = true;
			ReleaseSerialise();
			return(eBSFerrMem);
			}
		m_AllocGraphOutEdges += (UINT32)ReallocEdges;
		m_pGraphOutEdges = pTmp;
		}
	}

CurOverlapSeqID = 0;
CurOverlappedSeqID = 0;
pOverlap = pOverlappedSeqs;
pOutEdge = &m_pGraphOutEdges[m_UsedGraphOutEdges];
memset(pOutEdge,0,NumSeqs* sizeof(tsGraphOutEdge));
for(Idx = 0; Idx < NumSeqs; Idx++,pOverlap++,pOutEdge++)
	{
	// ensure up front that the sequence identifiers are known to have been associated with a vertex
	if(pOverlap->OverlappingSeqID != CurOverlapSeqID)
		{
		CurOverlapSeqID = pOverlap->OverlappingSeqID;
		pVertex = LocateVertexSeqID(CurOverlapSeqID);
		if(pVertex == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: unable to locate vertex for SeqID %u", CurOverlapSeqID);
			m_bTerminate = true;
			ReleaseSerialise();
			return(eBSFerrParams);
			}
		OverlappingVertexID = pVertex->VertexID;
		}
	if(pOverlap->OverlappedSeqID != CurOverlappedSeqID)
		{
		CurOverlappedSeqID = pOverlap->OverlappedSeqID;
		pVertex = LocateVertexSeqID(CurOverlappedSeqID);
		if(pVertex == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: unable to locate vertex for SeqID %u", CurOverlappedSeqID);
			m_bTerminate = true;
			ReleaseSerialise();
			return(eBSFerrParams);
			}
		OverlappedVertexID = pVertex->VertexID;
		}
	pOutEdge->FromVertexID = OverlappingVertexID;
	pOutEdge->ToVertexID = OverlappedVertexID;
	pOutEdge->SeqOfs = pOverlap->SeqOfs;
	pOutEdge->Contains = pOverlap->bContains ? 1 : 0;
	pOutEdge->OverlapSense = pOverlap->OverlapSense;
	}
m_UsedGraphOutEdges += NumSeqs;
m_bReduceEdges = true;
m_bOutEdgeSorted = false;
m_bInEdgeSorted = false;
ReleaseSerialise();
return(eBSFSuccess);
}

teBSFrsltCodes	 
CAssembGraph::AddEdge(tSeqID OverlappingSeqID, // identifies the overlapping sequence
				tSeqID OverlappedSeqID,			// identifies overlapped sequence or a completely contained sequence
				int SeqOfs,						// overlap (0..n) starts at this base relative to OverlappingSeqID 5'
				bool bContains,					// true if SeqID1 completely contains DnSeqID
				int OverlapSense)				// 0 sense overlaps sense, 1 antisense overlaps sense, 2 sense overlaps antisense 
{
size_t AllocMem;
tsGraphOutEdge *pOutEdge;
tsGraphVertex *pVertex;
tVertID OverlappingVertexID;
tVertID OverlappedVertexID;

if(OverlappingSeqID < 1 || OverlappedSeqID < 1 ||
	SeqOfs < 0)
	return(eBSFerrParams);

AcquireSerialise();
if(m_bTerminate)		// check if should early exit
	{
	ReleaseSerialise();
	return(eBSFerrInternal);
	}

if(m_UsedGraphOutEdges >= cMaxNumEdges)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge: too many edges, can only accept at most %u",cMaxNumEdges);
	m_bTerminate = true;
	ReleaseSerialise();
	return(eBSFerrMaxEntries);
	}

if(m_VerticesSortOrder != eVSOSeqID)		// vertices need to have been sorted by SeqID ascending
	{
	if(m_UsedGraphVertices > 1)
		{
		pStaticGraphVertices = m_pGraphVertices;
		m_MTqsort.qsort(pStaticGraphVertices,m_UsedGraphVertices,sizeof(tsGraphVertex),SortVerticesSeqID);
		}
	m_VerticesSortOrder = eVSOSeqID;
	}

// ensure up front that the sequence identifiers are known to have been associated with a vertex
pVertex = LocateVertexSeqID(OverlappingSeqID);
if(pVertex == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge: unable to locate vertex for SeqID %u", OverlappingSeqID);
	m_bTerminate = true;
	ReleaseSerialise();
	return(eBSFerrParams);
	}
OverlappingVertexID = pVertex->VertexID;

pVertex = LocateVertexSeqID(OverlappedSeqID);
if(pVertex == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge: unable to locate vertex for SeqID %u", OverlappedSeqID);
	m_bTerminate = true;
	ReleaseSerialise();
	return(eBSFerrParams);
	}
OverlappedVertexID = pVertex->VertexID;

if(m_pGraphOutEdges == NULL)				// initialisation may be required
	{
	AllocMem = cInitialAllocVertices * sizeof(tsGraphOutEdge);
#ifdef _WIN32
	m_pGraphOutEdges = (tsGraphOutEdge *) malloc(AllocMem);	
	if(m_pGraphOutEdges == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge: graph forward edge (%d bytes per edge) allocation of %d edges failed - %s",
								(int)sizeof(tsGraphOutEdge),cInitialAllocVertices,strerror(errno));
		m_bTerminate = true;
		ReleaseSerialise();
		return(eBSFerrMem);
		}
#else
	m_pGraphOutEdges = (tsGraphOutEdge *)mmap(NULL,AllocMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pGraphOutEdges == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge: graph forward edge (%d bytes per edge) allocation of %d edges failed - %s",
								(int)sizeof(tsGraphOutEdge),cInitialAllocVertices,strerror(errno));
		m_bTerminate = true;
		ReleaseSerialise();
		return(eBSFerrMem);
		}
#endif
	m_AllocGraphOutEdges = cInitialAllocVertices;
	m_UsedGraphOutEdges = 0;
	}
else
	{
	if((m_UsedGraphOutEdges + 4) > m_AllocGraphOutEdges) // need to add additional forward graph edges, 4 is simply to allow a small margin of error
		{
		tsGraphOutEdge *pTmp;
		UINT64 AllocEdges;
		UINT64 ReallocEdges =  (UINT64)((double)m_AllocGraphOutEdges * cReallocEdges);
		AllocEdges = (UINT64)m_AllocGraphOutEdges + ReallocEdges;
		if(AllocEdges > 0x0ffffffff)
			AllocEdges = 0x0ffffffff;
		AllocMem = (size_t)(AllocEdges * sizeof(tsGraphOutEdge));
#ifdef _WIN32
		pTmp = (tsGraphOutEdge *)realloc(m_pGraphOutEdges,AllocMem);
#else
		pTmp = (tsGraphOutEdge *)mremap(m_pGraphOutEdges,m_AllocGraphOutEdges * sizeof(tsGraphOutEdge),AllocMem,MREMAP_MAYMOVE);
		if(pTmp == MAP_FAILED)
			pTmp = NULL;
#endif
		if(pTmp == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge: graph forward edge (%d bytes per edge) re-allocation to %lld edges from %lld failed - %s",
													(int)sizeof(tsGraphOutEdge),AllocMem,m_AllocGraphOutEdges,AllocEdges,strerror(errno));
			m_bTerminate = true;
			ReleaseSerialise();
			return(eBSFerrMem);
			}
		m_AllocGraphOutEdges += (UINT32)ReallocEdges;
		m_pGraphOutEdges = pTmp;
		}
	}

// at least one unused outgoing edge available
pOutEdge = &m_pGraphOutEdges[m_UsedGraphOutEdges++];
memset(pOutEdge,0,sizeof(tsGraphOutEdge));
m_bOutEdgeSorted = false;
m_bInEdgeSorted = false;
pOutEdge->FromVertexID = OverlappingVertexID;
pOutEdge->ToVertexID = OverlappedVertexID;
pOutEdge->SeqOfs = SeqOfs;
pOutEdge->Contains = bContains ? 1 : 0;
pOutEdge->OverlapSense = OverlapSense;
m_bReduceEdges = true;
ReleaseSerialise();
return(eBSFSuccess);
}

UINT32										// returned vertex identifier
CAssembGraph::AddVertex(UINT32 SeqLen,		// sequence length
						tSeqID SeqID,		// allocates and initialises a new graph vertex which references this sequence
						tSeqID PESeqID)		// if paired end assembly then the paired (partner) sequence identifier, 0 if no paired end

{
size_t AllocMem;
tVertID VertexID;
tsGraphVertex *pVertex;
AcquireSerialise();

if(m_UsedGraphVertices >= cMaxNumVertices)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeqOverlap: too many vertices, can only accept at most %u",cMaxNumVertices);
	ReleaseSerialise();
	Reset();
	return(eBSFerrMaxEntries);
	}

if(m_pGraphVertices == NULL)		// initial allocation may be required
	{
	AllocMem = cInitialAllocVertices * sizeof(tsGraphVertex);
#ifdef _WIN32
	m_pGraphVertices = (tsGraphVertex *) malloc(AllocMem);	
	if(m_pGraphVertices == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddVertex: graph vertices (%d bytes per vertex) allocation of %u vertices failed - %s",
													(int)sizeof(tsGraphVertex),cInitialAllocVertices,strerror(errno));
		ReleaseSerialise();
		Reset();
		return(eBSFerrMem);
		}
#else
	m_pGraphVertices = (tsGraphVertex *)mmap(NULL,AllocMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pGraphVertices == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddVertex: graph vertices (%d bytes per vertex) allocation of %u vertices failed - %s",
													(int)sizeof(tsGraphVertex),cInitialAllocVertices,strerror(errno));
		m_pGraphVertices = NULL;
		ReleaseSerialise();
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocGraphVertices = cInitialAllocVertices;
	m_UsedGraphVertices = 0;
	}
else
	{
	if(m_UsedGraphVertices >= m_AllocGraphVertices) // need to alloc additional graph vertices?
		{
		tsGraphVertex *pTmp;
		UINT32 ReallocVertices;
		size_t memreq;
		ReallocVertices = (UINT32)(m_AllocGraphVertices * cReallocVertices);

		memreq = (m_AllocGraphVertices + ReallocVertices) * sizeof(tsGraphVertex);
#ifdef _WIN32
		pTmp = (tsGraphVertex *)realloc(m_pGraphVertices,memreq);
#else
		pTmp = (tsGraphVertex *)mremap(m_pGraphVertices,m_AllocGraphVertices * sizeof(tsGraphVertex),memreq,MREMAP_MAYMOVE);
		if(pTmp == MAP_FAILED)
			pTmp = NULL;
#endif
		if(pTmp == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddVertex: graph vertices (%d bytes per vertex) re-allocation from %u to %u failed - %s",
														(int)sizeof(tsGraphVertex),m_UsedGraphVertices,m_UsedGraphVertices + ReallocVertices,strerror(errno));
			ReleaseSerialise();
			return(eBSFerrMem);
			}
		memset(&pTmp[m_AllocGraphVertices],0,memreq - (m_AllocGraphVertices * sizeof(tsGraphVertex)));
		m_AllocGraphVertices += ReallocVertices;
		m_pGraphVertices = pTmp;
		}
	}
m_UsedGraphVertices+=1;
m_VerticesSortOrder = eVSOUnsorted;
VertexID = m_UsedGraphVertices;
pVertex = &m_pGraphVertices[VertexID-1];
memset(pVertex,0,sizeof(tsGraphVertex));
pVertex->VertexID = VertexID;
pVertex->PEVertexID = (tVertID)PESeqID;
pVertex->SeqID = SeqID;
pVertex->SeqLen = SeqLen;
ReleaseSerialise();
return(VertexID);
}

tsGraphVertex *			// ptr to vertex corresponding to SeqID , or NULL if unable to locate				
CAssembGraph::LocateVertexSeqID(tSeqID SeqID)		// match this SeqID
{
tsGraphVertex *pEl2;
int CmpRslt;
UINT32 TargPsn;
UINT32 NodeLo;				
UINT32 NodeHi;				
if( SeqID < 1,								// must be valid
	m_pGraphVertices == NULL,				// must be allocated
  	m_UsedGraphVertices < 1 ||				// must have at least 1 vertex
	m_VerticesSortOrder != eVSOSeqID)		// and vertices must be sorted by SeqID ascending
	return(NULL);

NodeLo = 0;
NodeHi = m_UsedGraphVertices - 1;
do {
	TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;
	pEl2 = &m_pGraphVertices[TargPsn];

	if(SeqID == pEl2->SeqID)
		return(pEl2);

	if(SeqID > pEl2->SeqID)
		CmpRslt = 1;
	else
		CmpRslt = -1;

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		NodeHi = TargPsn - 1;
		}
	else
		NodeLo = TargPsn+1;
	}
while(NodeHi >= NodeLo);

return(0);	// unable to locate any instance of SeqID
}

tsGraphVertex *			// ptr to vertex corresponding to VertexID , or NULL if unable to locate				
CAssembGraph::LocateVertex(tVertID VertexID)		// match this VertexID
{
tsGraphVertex *pEl2;
int CmpRslt;
UINT32 TargPsn;
UINT32 NodeLo;				
UINT32 NodeHi;				
if( VertexID < 1,								// must be valid
	m_pGraphVertices == NULL,				// must be allocated
  	m_UsedGraphVertices < 1 ||				// must have at least 1 vertex
	m_VerticesSortOrder != eVSOVertexID)		// and vertices must be sorted by VertexID ascending
	return(NULL);

NodeLo = 0;
NodeHi = m_UsedGraphVertices - 1;
do {
	TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;
	pEl2 = &m_pGraphVertices[TargPsn];

	if(VertexID == pEl2->VertexID)
		return(pEl2);

	if(VertexID > pEl2->VertexID)
		CmpRslt = 1;
	else
		CmpRslt = -1;

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		NodeHi = TargPsn - 1;
		}
	else
		NodeLo = TargPsn+1;
	}
while(NodeHi >= NodeLo);

return(0);	// unable to locate any instance of VertexID
}


// when all sequences have been associated to vertices then FinaliseVertices must be called which
// will firstly sort the vertices in SeqID ascending order and then remap the PEVertexIDs ready for edges between adjacent vertexes to be added with AddEdge()
UINT32									// 0 if errors else number of vertices finalised
CAssembGraph::FinaliseVertices(void)
{
UINT32 Idx;
tsGraphVertex *pVertex;
tsGraphVertex *pPEVertex;
if(m_pGraphVertices == NULL ||				// must be allocated
  	m_UsedGraphVertices < 1)				// must have at least 1 vertex
	return(0);
if(m_VerticesSortOrder != eVSOSeqID)
	{
	if(m_UsedGraphVertices > 1)
		{
		pStaticGraphVertices = m_pGraphVertices;
		m_MTqsort.qsort(pStaticGraphVertices,m_UsedGraphVertices,sizeof(tsGraphVertex),SortVerticesSeqID);
		}
	m_VerticesSortOrder = eVSOSeqID;
	}
// now remap the PEVertexIDs from the AddVertex() initialised PESeqID to actual PEVertexIDs
pVertex = m_pGraphVertices;
for(Idx = 0; Idx < m_UsedGraphVertices; Idx++,pVertex++)
	{
	if(pVertex->PEVertexID != 0)
		{
		if((pPEVertex = LocateVertexSeqID((tSeqID)pVertex->PEVertexID))!=NULL)
			pVertex->PEVertexID = pPEVertex->VertexID;
		}
	}

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
return(m_UsedGraphVertices);
}



// when all edges have been added then FinaliseEdges must be called
UINT32									// 0 if errors else number of edges finalised
CAssembGraph::FinaliseEdges(void)
{
tsGraphOutEdge *pOutEdge;
tEdgeID *pInEdge;
UINT32 EdgeIdx;
tVertID CurVertexID;
tsGraphVertex *pVertex;
#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif

if(m_pGraphVertices == NULL || m_UsedGraphVertices < 2 ||
   m_pGraphOutEdges == NULL || m_UsedGraphOutEdges < 1)
	return(0);
if(m_bReduceEdges)
	ReduceEdges();
m_bReduceEdges = false;

// need vertices sorted by vertex identifier ascending
if(m_VerticesSortOrder != eVSOVertexID)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %u vertices ...",m_UsedGraphVertices);
	if(m_UsedGraphVertices > 1)
		{
		pStaticGraphVertices = m_pGraphVertices;
		m_MTqsort.qsort(pStaticGraphVertices,m_UsedGraphVertices,sizeof(tsGraphVertex),SortVerticesVertexID);
		}
	m_VerticesSortOrder = eVSOVertexID;
	m_bOutEdgeSorted = false;
	m_bInEdgeSorted = false;
	m_bVertexEdgeSet = false;
	}

// ensure outgoing edges are sorted FromVertexID.ToVertexID ascending order
if(!m_bOutEdgeSorted)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %u outgoing edges ...",m_UsedGraphOutEdges);
	if(m_UsedGraphOutEdges >= 2)
		{
		pStaticGraphOutEdges = m_pGraphOutEdges;
		m_MTqsort.qsort(pStaticGraphOutEdges,m_UsedGraphOutEdges,sizeof(tsGraphOutEdge),SortOutEdgeFromVertexID);
		}
	m_bOutEdgeSorted = true;
	m_bInEdgeSorted = false;
	m_bVertexEdgeSet = false;
	}

if(!m_bInEdgeSorted)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assigning %u incoming edges ...",m_UsedGraphOutEdges);
	// firstly allocate to hold the incoming edges
	size_t AllocMem;
	if(m_pGraphInEdges == NULL)		// initial allocation may be required
		{
		m_UsedGraphInEdges = 0; 
		m_AllocGraphInEdges = m_UsedGraphOutEdges;	
		AllocMem = m_AllocGraphInEdges * sizeof(tEdgeID);
#ifdef _WIN32
		m_pGraphInEdges = (tEdgeID *) malloc(AllocMem);	
		if(m_pGraphInEdges == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"FinaliseEdges: graph vertex input edges (%d bytes per edge) allocation of %u edges failed - %s",
														(int)sizeof(tEdgeID),m_AllocGraphInEdges,strerror(errno));
			Reset();
			return(eBSFerrMem);
			}
#else
		m_pGraphInEdges = (tEdgeID *)mmap(NULL,AllocMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
		if(m_pGraphInEdges == MAP_FAILED)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"FinaliseEdges: graph vertex input edges (%d bytes per edge) allocation of %u edges failed - %s",
														(int)sizeof(tEdgeID),m_AllocGraphInEdges,strerror(errno));
			m_pGraphInEdges = NULL;
			Reset();
			return(eBSFerrMem);
			}
#endif
		}
	else
		{
		if(m_AllocGraphInEdges < m_UsedGraphOutEdges)	// need to alloc additional graph vertices?
			{
			tEdgeID *pTmp;
			size_t memreq;
			memreq = m_UsedGraphOutEdges * sizeof(tEdgeID);
#ifdef _WIN32
			pTmp = (tEdgeID *)realloc(m_pGraphInEdges,memreq);
#else
			pTmp = (tEdgeID *)mremap(m_pGraphInEdges,m_AllocGraphInEdges * sizeof(tEdgeID),memreq,MREMAP_MAYMOVE);
			if(pTmp == MAP_FAILED)
				pTmp = NULL;
#endif
			if(pTmp == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"FinaliseEdges: graph input edges (%d bytes per edge) re-allocation from %u to %u failed - %s",
															(int)sizeof(tEdgeID),m_AllocGraphInEdges,m_UsedGraphOutEdges,strerror(errno));
				return(eBSFerrMem);
				}
			m_AllocGraphInEdges = m_UsedGraphOutEdges;
			m_pGraphInEdges = pTmp;
			}
		}
	m_UsedGraphInEdges = m_UsedGraphOutEdges;
	pInEdge = m_pGraphInEdges;
	tEdgeID Idx;
	for(Idx = 1; Idx <= m_UsedGraphInEdges; Idx++, pInEdge++)
		*pInEdge = Idx;
	pStaticGraphInEdges = m_pGraphInEdges;
	pStaticGraphOutEdges = m_pGraphOutEdges;
	m_MTqsort.qsort(pStaticGraphInEdges,m_UsedGraphInEdges,sizeof(tEdgeID),SortInEdgesToVertexID);
	m_bInEdgeSorted = true;
	}

if(!m_bVertexEdgeSet)
	{
	// iterate all outgoing edges and update vertices with starting outgoing edges
	CurVertexID = 0;
	pOutEdge = m_pGraphOutEdges;
	for(EdgeIdx = 1; EdgeIdx <= m_UsedGraphOutEdges; EdgeIdx++,pOutEdge++)
		{
		if(CurVertexID == pOutEdge->FromVertexID)
			continue;
		CurVertexID = pOutEdge->FromVertexID;
		pVertex = &m_pGraphVertices[CurVertexID-1];
		pVertex->OutEdgeID = EdgeIdx; 
		}

	// iterate all incoming edges and update vertices with starting incoming edges
	CurVertexID = 0;
	pInEdge = m_pGraphInEdges;
	for(EdgeIdx = 1; EdgeIdx <= m_UsedGraphInEdges; EdgeIdx++,pInEdge++)
		{
		pOutEdge = &m_pGraphOutEdges[*pInEdge - 1];
		if(CurVertexID == pOutEdge->ToVertexID)
			continue;
		CurVertexID = pOutEdge->ToVertexID;
		pVertex = &m_pGraphVertices[CurVertexID-1];
		pVertex->InEdgeID = EdgeIdx; 
		}

	m_bVertexEdgeSet = true;
	}

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
return(m_UsedGraphOutEdges);
}

UINT32 
CAssembGraph::GetNumGraphVertices(void)		// returns current number of graph vertices
{
UINT32 Num;
AcquireSerialise();
Num = m_UsedGraphVertices;
ReleaseSerialise();
return(Num);
}

UINT32 
CAssembGraph::GetNumGraphOutEdges(void)		// returns current number of graph forward edges
{
UINT32 Num;
AcquireSerialise();
Num = m_UsedGraphOutEdges;
ReleaseSerialise();
return(Num);
}

UINT32 
CAssembGraph::GetNumReducts(void)		// returns current number of edge reductions
{
UINT32 Num;
AcquireSerialise();
Num = m_NumReducts;
ReleaseSerialise();
return(Num);
}


UINT32								
CAssembGraph::ClearEdgeTravFwdRevs(void)
{
tsGraphOutEdge *pEdge;
UINT32 EdgeIdx;
pEdge = m_pGraphOutEdges;
for(EdgeIdx = 0; EdgeIdx < m_UsedGraphOutEdges; EdgeIdx++,pEdge++)
	{
	pEdge->TravFwd = 0;
	pEdge->TravRev = 0;
	}
return(m_UsedGraphOutEdges);
}

UINT32								
CAssembGraph::ClearDiscCompIDs(void)
{
tsGraphVertex *pVertex;
UINT32 VertexIdx;
pVertex = m_pGraphVertices;
for(VertexIdx = 0; VertexIdx < m_UsedGraphVertices; VertexIdx++,pVertex++)
	pVertex->DiscGraphID = 0;
return(m_UsedGraphVertices);
}

UINT32									// total number of vertices with both inbound and outbound edges
CAssembGraph::VertexConnections(void)	// identify and mark vertices which have multiple in/out bound edges
{
UINT32 NumOutEdges;
UINT32 NumInEdges;
UINT32 NumMultiInVertices;
UINT32 NumMultiOutVertices;
UINT32 NumMultiInOutVertices;
UINT32 NumIsolatedVertices;
UINT32 NumStartVertices;
UINT32 NumEndVertices;
UINT32  NumInternVertices;
tEdgeID EdgeID;
tVertID VertexID;
tsGraphVertex *pVertex;
tsGraphOutEdge *pEdge;

NumMultiInVertices = 0;
NumMultiOutVertices = 0;
NumMultiInOutVertices = 0;
NumIsolatedVertices = 0;
NumStartVertices = 0;
NumInternVertices = 0;
NumEndVertices = 0;
pVertex = m_pGraphVertices;
for(VertexID = 1; VertexID <= m_UsedGraphVertices; VertexID++, pVertex++)
	{
	NumOutEdges = 0;
	NumInEdges = 0;
	if(pVertex->OutEdgeID != 0)				// check for multiple outgoing edges
		{
		EdgeID = pVertex->OutEdgeID;
		while(EdgeID <= m_UsedGraphOutEdges)
			{
			pEdge = &m_pGraphOutEdges[EdgeID-1];
			if(pEdge->FromVertexID != VertexID)
				break;
			NumOutEdges += 1;
			EdgeID += 1;
			}
		}
	pVertex->DegreeOut = min(15,NumOutEdges);

	if(pVertex->InEdgeID != 0)				// check for multiple incoming edges
		{
		EdgeID = pVertex->InEdgeID;
		while(EdgeID <= m_UsedGraphInEdges)
			{
			pEdge = &m_pGraphOutEdges[m_pGraphInEdges[EdgeID-1]-1];
			if(pEdge->ToVertexID != VertexID)
				break;
			NumInEdges += 1;
			EdgeID += 1;
			}
		}
	pVertex->DegreeIn = min(15,NumInEdges);

	if(NumOutEdges == 0 && NumInEdges == 0)
		{
		NumIsolatedVertices += 1;
		continue;
		}
	if(NumOutEdges > 1)
		NumMultiOutVertices += 1;
	if(NumInEdges > 1)
		NumMultiInVertices += 1;
	if(NumOutEdges == 0)
		NumEndVertices += 1;
	if(NumInEdges == 0)
		NumStartVertices += 1;
	if(NumInEdges > 0 && NumOutEdges > 0)
		NumInternVertices += 1;
	if(NumInEdges > 1 && NumOutEdges > 1)
		NumMultiInOutVertices += 1;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Vertices degree of connectivity: %u Isolated, %u Start, %u End, %u Internal, %u MultiDegreeIn, %u MultiDegreeOut, %u MultiInOut",
						NumIsolatedVertices,NumStartVertices,NumEndVertices,NumInternVertices,NumMultiInVertices,NumMultiOutVertices,NumMultiInOutVertices);
return(NumInternVertices);
}


static UINT32 m_NumReplacedDiscGraphIDs;
// IdentifyDisconnectedSubGraphs
// Within the graph there are likely to be many (could be millions) of completely disconnected subgraphs
// These disconnected subgraphs have no sequences which overlay, or are overlaid by, sequences in any other subgraph 
// This function iterates all vertices of the graph and locates all other vertices which are connected to the orginal vertice and marks these
// as belonging to an disconnected subgraph
// Currently this function is single threaded, could be worthwhile to multithread
// 
UINT32						// returned number of subgraphs identified
CAssembGraph::IdentifyDisconnectedSubGraphs(void)
{
UINT32 VertexIdx;
UINT32 NumVertices;
UINT32 MaxVertices;
UINT32 Cntr;
tDiscGraphID CurDiscGraphID;
tsGraphVertex *pVertex;

// need to have at least 4 vertices and 2 edges ( min for 2 disconnected graph components)
if(m_pGraphVertices == NULL || m_UsedGraphVertices < 4 ||
   m_pGraphOutEdges == NULL || m_UsedGraphOutEdges < 2)
	return(0);

if(!FinaliseEdges())
	return(0);

ClearEdgeTravFwdRevs();
ClearDiscCompIDs();

// initial alloc for the identified components as may be required
if(m_pComponents == NULL)
	{
	size_t AllocMem = (size_t)cInitalComponentsAlloc * sizeof(tsComponent);
#ifdef _WIN32
	m_pComponents = (tsComponent *) malloc(AllocMem);	
	if(m_pComponents == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifyDisconnectedSubGraphs: components (%d bytes per entry) allocation of %d entries failed - %s",
								(int)sizeof(tsComponent),cInitalComponentsAlloc,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	m_pComponents = (tsComponent *)mmap(NULL,AllocMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pComponents == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifyDisconnectedSubGraphs: components (%d bytes per entry) allocation of %d entries failed - %s",
								(int)sizeof(tsComponent),cInitalComponentsAlloc,strerror(errno));
		m_pComponents = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocComponents = cInitalComponentsAlloc;
	m_UsedComponents = 0;
	}
// determine, flag and report, on vertix degree of connectivity
VertexConnections();

// iterate vertices
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identifying disconnected graph components ...");
m_NumDiscRemaps = 0;
m_CurTransitDepth = 0;
CurDiscGraphID = 0;
MaxVertices = 0;
Cntr = 0;
pVertex = m_pGraphVertices;
time_t Started = time(0);
for(VertexIdx = 0; VertexIdx < m_UsedGraphVertices; VertexIdx++, pVertex++)
	{
	if(!(VertexIdx % 100))
		{
		time_t Now = time(0);
		unsigned long ElapsedSecs = (unsigned long) (Now - Started);
		if(ElapsedSecs >= 60)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processed %d vertices (%1.2f%%), identified components: %u, max vertices: %u",
								VertexIdx,(100.0 * (double)VertexIdx)/m_UsedGraphVertices,CurDiscGraphID,MaxVertices);
			Started = Now;
			}
		}
	if(pVertex->DiscGraphID > 0 && pVertex->DiscGraphID != (tDiscGraphID)-1)	// if vertex already associated to a disconnected graph then try next vertex
		continue;

	NumVertices = IdentifyDiscComponent(pVertex->VertexID,CurDiscGraphID + 1);
	if(NumVertices > 0)
		{
		// realloc for the identified components as may be required
		if((m_UsedComponents + 16) >= m_AllocComponents)
			{
			size_t AllocMem;
			tsComponent *pTmp;
			UINT32 ReallocComponents;
			ReallocComponents = (UINT32)(m_AllocComponents * cReallocComponents);
			AllocMem = (size_t)((m_AllocComponents + ReallocComponents) * sizeof(tsComponent));
		#ifdef _WIN32
			pTmp = (tsComponent *)realloc(m_pComponents,AllocMem);
		#else
			pTmp = (tsComponent *)mremap(m_pComponents,m_AllocComponents * sizeof(tsComponent),AllocMem,MREMAP_MAYMOVE);
			if(pTmp == MAP_FAILED)
				pTmp = NULL;
		#endif
			if(pTmp == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"PushTransitStack: graph Transit stack (%d bytes per entry) re-allocation to %lld from %lld failed - %s",
																	(int)sizeof(tsComponent),m_AllocComponents  + (UINT64)ReallocComponents,m_AllocComponents,strerror(errno));
				return(eBSFerrMem);
				}
			m_AllocComponents += ReallocComponents;
			m_pComponents = pTmp;
			}

		m_pComponents[m_UsedComponents].ComponentID = CurDiscGraphID + 1;
		m_pComponents[m_UsedComponents].VertexID = pVertex->VertexID;
		m_pComponents[m_UsedComponents++].NumVertices = NumVertices;

		CurDiscGraphID += 1;
		if(NumVertices > MaxVertices)
			MaxVertices = NumVertices;
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Number of disconnected graph components: %u, max vertices in any graph: %u",CurDiscGraphID,MaxVertices);


// sort the components by NumVertices and report the top 100
if(m_UsedComponents)
	{
	if(m_UsedComponents > 1)
		{
		pStaticComponents = m_pComponents;
		m_MTqsort.qsort(pStaticComponents,m_UsedComponents,sizeof(tsComponent),SortComponentNumVertices);
		}
	for(UINT32 Idx = 0; Idx < min(100,m_UsedComponents); Idx++)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"ComponentID %u, Vertices %u",m_pComponents[Idx].ComponentID,m_pComponents[Idx].NumVertices);
		}
	}
return(CurDiscGraphID);
}

int			// stack depth or < 1 if errors
CAssembGraph::PushTransitStack(tVertID VertexID)
{
// initial alloc for the transitative stack as may be required
if(m_pTransitStack == NULL)
	{
	size_t AllocMem = (size_t)cTransitStackAlloc * sizeof(tVertID);
#ifdef _WIN32
	m_pTransitStack = (tVertID *) malloc(AllocMem);	
	if(m_pTransitStack == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"PushTransitStack: Transit stack (%d bytes per transition) allocation of %d entries failed - %s",
								(int)sizeof(tVertID),cTransitStackAlloc,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	m_pTransitStack = (tVertID *)mmap(NULL,AllocMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pTransitStack == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifyDisconnectedSubGraphs: Transit stack (%d bytes per transition) allocation of %d entries failed - %s",
								(int)sizeof(tVertID),cTransitStackAlloc,strerror(errno));
		m_pTransitStack = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocTransitStack = cTransitStackAlloc;
	}
else
	{
				// check if need to deepen transition stack
	if((m_CurTransitDepth + 16) >= m_AllocTransitStack)
		{
		size_t AllocMem;
		tVertID *pTmp;
		UINT32 ReallocTransitStack;
		ReallocTransitStack = (UINT32)(m_AllocTransitStack * cReallocTransitStack); 
		AllocMem = (size_t)((m_AllocTransitStack + ReallocTransitStack) * sizeof(tVertID));
	#ifdef _WIN32
		pTmp = (tVertID *)realloc(m_pTransitStack,AllocMem);
	#else
		pTmp = (tVertID *)mremap(m_pTransitStack,m_AllocTransitStack * sizeof(tVertID),AllocMem,MREMAP_MAYMOVE);
		if(pTmp == MAP_FAILED)
			pTmp = NULL;
	#endif
		if(pTmp == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"PushTransitStack: graph Transit stack (%d bytes per entry) re-allocation to %lld from %lld failed - %s",
																(int)sizeof(tVertID),m_AllocTransitStack  + (UINT64)ReallocTransitStack,m_AllocTransitStack,strerror(errno));
			return(eBSFerrMem);
			}
		m_AllocTransitStack += ReallocTransitStack;
		m_pTransitStack = pTmp;
		}
	}
m_pTransitStack[m_CurTransitDepth++] = VertexID;
return((int)m_CurTransitDepth);
}

tVertID
CAssembGraph::PopTransitStack(void)
{
if(m_CurTransitDepth)
	return(m_pTransitStack[--m_CurTransitDepth]);
return(0);
}

tVertID
CAssembGraph::PeekTransitStack(void)
{
if(m_CurTransitDepth)
	return(m_pTransitStack[m_CurTransitDepth]);
return(0);
}

void
CAssembGraph::ClearTransitStack(void)
{
m_CurTransitDepth = 0;
}


UINT32												 // number of vertices marked as members of this component
CAssembGraph::IdentifyDiscComponent(tVertID VertexID, // start component traversal from this vertex
				tDiscGraphID DiscGraphID)			 // mark all traversed vertices as members of this component
{
bool bTraverse;
UINT32 NumMembers;
tEdgeID EdgeID;
tVertID CurVertexID;
tsGraphVertex *pVertex;
tsGraphOutEdge *pEdge;

if(VertexID == 0 || VertexID > m_UsedGraphVertices)
	return(0);
pVertex = &m_pGraphVertices[VertexID-1];

// if vertex already member of component then simply return without further processing
if(pVertex->DiscGraphID != 0)
	return(0);

pVertex->DiscGraphID = DiscGraphID;


// possible that the initial vertex does not have any incomimg or outgoing edges to any other another vertex
if(pVertex->InEdgeID == 0 && pVertex->OutEdgeID == 0)
	{
	pVertex->DiscGraphID = DiscGraphID;
	return(1);
	}

// vertex has at least one, outgoing or incoming, edge to or from some other vertex
ClearTransitStack();
PushTransitStack(VertexID);
NumMembers = 0;

while((CurVertexID = PopTransitStack()) != 0)
	{
	bTraverse = false;
	pVertex = &m_pGraphVertices[CurVertexID-1];

	// for reasons unknown (Ok, a bug) occasionally (perhaps 1 in 100 million) a vertex has been classified into a different component and
	// it is worth noting that this different component identifier is that of a component not previously classified so I think that somehow there was
	// a memory overwrite ... 
	// currently just report this bug and continue...
	if(pVertex->DiscGraphID != 0 && pVertex->DiscGraphID != DiscGraphID)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"IdentifyDiscComponent: Found vertex %u classified as %u when processing for component %u",CurVertexID, pVertex->DiscGraphID,DiscGraphID);
	
	if(pVertex->DiscGraphID == 0 || pVertex->DiscGraphID > DiscGraphID)
		{
		pVertex->DiscGraphID = DiscGraphID;
		NumMembers += 1;
		}
	if(pVertex->OutEdgeID != 0)				// try outgoing not yet traversed ...
		{
		EdgeID = pVertex->OutEdgeID;
		while(EdgeID <= m_UsedGraphOutEdges)
			{
			pEdge = &m_pGraphOutEdges[EdgeID-1];
			if(pEdge->FromVertexID != CurVertexID)
				break;
			if(!pEdge->TravFwd)
				{
				pEdge->TravFwd = 1;
				bTraverse = true;
				break;
				}
			EdgeID += 1;
			}
		if(bTraverse)
			{
			PushTransitStack(CurVertexID);
			PushTransitStack(pEdge->ToVertexID);
			continue;
			}
		}

	if(pVertex->InEdgeID != 0)				// if incoming not yet traversed then ...
		{
		EdgeID = pVertex->InEdgeID;
		while(EdgeID <= m_UsedGraphInEdges)
			{
			pEdge = &m_pGraphOutEdges[m_pGraphInEdges[EdgeID-1]-1];
			if(pEdge->ToVertexID != CurVertexID)
				break;
			if(!pEdge->TravRev)
				{
				pEdge->TravRev = 1;
				bTraverse = true;
				break;
				}
			EdgeID += 1;
			}
		if(bTraverse)
			{
			PushTransitStack(CurVertexID);
			PushTransitStack(pEdge->FromVertexID);
			continue;
			}
		}
	}
return(NumMembers);
}


// treat initial vertex as sense with all other vertices in component having a relative sense
// Set DirContext = 0
//
// If DirContext == 0
//	If overlapped read is sense then that read is downstream, DirContext stays same*-
//	If overlapped read was antisense then that read is downstream, DirContext = 1
// If DirContext == 1
//	If overlapped read is sense then that read is upstream, DirContext stays same 
//  If overlapped read is antisense then that read is downstream, DirContext = 0
 


// 
// Process all vectices in specified component and generate output fragments
// Generated fragments contain vertice sequences whereby the 5' vertex has 0 or more than 1 inbound edge, intermediate vertexs 1 inbound and 1 outbound edge, and
// the 5' vertex sequence had 0 or more than 1 outbound edge
UINT32
CAssembGraph::GenSeqFragment(tsGraphVertex *pVertex)		// initial seed vertex
{
UINT32 FragLen;
tsGraphOutEdge *pEdge;
bool bSense = true;

if(pVertex->flgEmitted)			// can't emit internal vertices (0 or 1 inbound and 0 or 1 outbound edges) multiple times
	{
	if(pVertex->DegreeIn <= 1 && pVertex->DegreeOut <= 1)
		return(0);
	}

FragLen = 0;

// follow back links to the vertex representing the 5' end of the fragment
while(pVertex->DegreeIn == 1)
	{
	pEdge = &m_pGraphOutEdges[m_pGraphInEdges[pVertex->InEdgeID-1]-1];
	pVertex = &m_pGraphVertices[pEdge->FromVertexID-1];
	}

// iterate towards 3' of fragment emitting sequence
// stop emiting sequence when vertex has 0 or more than 1 outgoing edges
do
	{
	// emit sequence for vertex here ... .... ...

	pVertex->flgEmitted = 1;
	if(pVertex->DegreeOut != 1)
		break;
	pEdge = &m_pGraphOutEdges[m_pGraphInEdges[pVertex->OutEdgeID-1]-1];
	pVertex = &m_pGraphVertices[pEdge->ToVertexID-1];
	}
while(1);
return(FragLen);
}


UINT32					// number of replacements
CAssembGraph::ReplaceDiscGraphIDs(void) 
{
UINT32 VertexIdx;
tsGraphVertex *pVertex;
tsRemapDiscGraphID *pRemap;
bool bReplaced;
UINT32 RemapIdx;
UINT32 NumReplaced = 0;

if(!m_NumDiscRemaps)
	return(0);

pVertex = m_pGraphVertices;
for(VertexIdx = 0; VertexIdx < m_UsedGraphVertices; VertexIdx++, pVertex++)
	{
	pRemap = m_DiscRemaps;
	bReplaced = false;
	for(RemapIdx = 0; RemapIdx < m_NumDiscRemaps; RemapIdx++,pRemap++)
		{
		if(pRemap->From == pVertex->DiscGraphID)
			{
			pVertex->DiscGraphID = pRemap->To;
			bReplaced = true;
			}
		}
	if(bReplaced)
		NumReplaced += 1;
	}
m_NumDiscRemaps = 0;
return(NumReplaced);
}

UINT32					// number of replacements
CAssembGraph::ReplaceDiscGraphID(tDiscGraphID ToReplaceID,		// replace existing disconnected graph identifiers 
		tDiscGraphID ReplaceID)					// with this identifier
{
UINT32 VertexIdx;
tsGraphVertex *pVertex;
UINT32 NumReplaced = 0;
pVertex = m_pGraphVertices;
for(VertexIdx = 0; VertexIdx < m_UsedGraphVertices; VertexIdx++, pVertex++)
	{
	if(pVertex->DiscGraphID == ToReplaceID)
		{
		pVertex->DiscGraphID = ReplaceID;
		NumReplaced += 1;
		}
	}
return(NumReplaced);
}


UINT32			// index+1 in m_pGraphOutEdges of first matching FromVertexID, or 0 if non matching				
CAssembGraph::LocateFirstVertID(tVertID FromVertexID)	// find first matching  
{
tsGraphOutEdge *pEl2;
int CmpRslt;
UINT32 Mark;
UINT32 TargPsn;
UINT32 NodeLo = 0;
UINT32 NodeHi = m_UsedGraphOutEdges-1;

if(m_bOutEdgeSorted == false)		// out edges must be sorted by FromVertexID ascending
	return(0);

do {
	TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;
	pEl2 = &m_pGraphOutEdges[TargPsn];

	if(FromVertexID > pEl2->FromVertexID)
		CmpRslt = 1;
	else
		if(FromVertexID < pEl2->FromVertexID)
			CmpRslt = -1;
		else
			CmpRslt = 0;

	if(!CmpRslt)	// if a match then may not be the lowest indexed match
		{
		if(TargPsn == 0 || NodeLo == TargPsn) // check if already lowest
			return(TargPsn + 1);
		// iterate until lowest located
		while(1) {
			if(CmpRslt == 0)
				{
				Mark = TargPsn;
				if(Mark == 0)
					return(Mark+1);
				NodeHi = TargPsn - 1;
				}	
			TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;

			pEl2 = &m_pGraphOutEdges[TargPsn];
			if(FromVertexID > pEl2->FromVertexID)
				CmpRslt = 1;
			else
				if(FromVertexID < pEl2->FromVertexID)
					CmpRslt = -1;
				else
					{
					CmpRslt = 0;
					continue;
					}
			NodeLo = TargPsn + 1;
			if(NodeLo == Mark)
				return(Mark+1);
			}
		}

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		NodeHi = TargPsn - 1;
		}
	else
		NodeLo = TargPsn+1;
	}
while(NodeHi >= NodeLo);

return(0);	// unable to locate any instance of FromVertexID
}



UINT64 
CAssembGraph::WriteContigs(char *pszOutFile)  // write assembled contigs to this output file
{
// given an initial starting node then objective is to locate an upstream head node 
//
// OverlayState = SenseToSense
// 
// If OverlayState == SenseToSense
//		Accept either maximal overlaid by Sense or AntiSense
//			If maximal overlaid was Sense-->Sense then
//               OverlayState = SenseToSense
//               CurReadID = overlaying read identifier
//			If maximal overlaid was AntiSense-->Sense then
//               OverlayState = AntiSenseToSense
//               CurReadID = overlaying read identifier
//

// objective to find longest traversal path through related nodes from a head node through to tail node; this is a contig
// assume a head node sequence has been identified - this is overlaid by no other sequences
// CurReadID = Head sequence identifier
// OverlayState = SenseToSense
// Start the contig with the head sequence
// 
// if OverlayState == SenseToSense
//     Accept the read sequence which is maximally overlaid either as either Sense-->Sense or Sense-->AntiSense
//			If maximal overlap is Sense-->Sense then
//				merge contig and overlaid sequence
//				OverlayState = SenseToSense
//				CurReadID = overlaid read identifier
//				continue
//			If maximal overlap is Sense-->AntiSense then
//				merge contig and overlaid revcpl sequence
//				OverlayState = SenseToAntiSense
//				CurReadID = overlaid read identifier
//				continue
//
// if OverlayState == SenseToAntiSense
//     Accept the read sequence which is maximally overlaid either as either AntiSense-->AntiSense or AntiSense-->Sense
//			If maximal overlap is AntiSense-->AntiSense
//				merge contig and revcpl overlaid sequence
//				OverlayState = SenseToAntiSense
//				CurReadID = overlaid read identifier
//				continue
//			If maximal overlap is AntiSense-->Sense then
//				merge contig and overlaid sequence
//				OverlayState = SenseToSense
//				CurReadID = overlaid read sequence identifier
//				continue

return(0);
}



// ReduceEdges
// As edges are added independently of all other edges then graph may contain many extraneous edges
// This function firstly identifies these extraneous edges and finally removes them
UINT32								// number of extraneous edges removed	 
CAssembGraph::ReduceEdges(void)		// reduce graph by detecting and removing extraneous edges
{
tsGraphOutEdge *pEdge;
tsGraphOutEdge *pDstEdge;
tVertID ToVertexID;
UINT16 SeqOfs;
UINT32 Idx;

UINT32 Num2Remove;
UINT32 NumRemoved;


if(m_pGraphOutEdges == NULL || m_UsedGraphOutEdges < 2)	// any edges to be reduced?
	return(0);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reduce %u edges starting ...",m_UsedGraphOutEdges);
m_NumReducts += 1;

// firstly, edges are sorted by ToVertexID.OverlapOfs ascending..
pStaticGraphOutEdges = m_pGraphOutEdges;
m_MTqsort.qsort(pStaticGraphOutEdges,m_UsedGraphOutEdges,sizeof(tsGraphOutEdge),SortOutEdgeToVertexIDSeqOfs);
m_bOutEdgeSorted = false;

// next identify those extraneous nodes to be removed
// retained edges are simply for all ToVertexIDs those with the shortest SeqOfs
Num2Remove = 0;
pEdge = m_pGraphOutEdges; 
ToVertexID = 0;
for(Idx = 0; Idx < m_UsedGraphOutEdges; Idx++,pEdge++)
	{
	if(pEdge->OverlapSense != 0)	
		continue;
	if(ToVertexID == pEdge->ToVertexID)
		{
		if(pEdge->SeqOfs > SeqOfs)
			{
			pEdge->bRemove = 1;
			Num2Remove += 1;
			continue;
			}
		}
	SeqOfs = pEdge->SeqOfs;
	ToVertexID = pEdge->ToVertexID;
	pEdge->bRemove = 0;
	}


if(!Num2Remove)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"ReduceEdges() completed, removed 0 edges");
	return(0);
	}

// extraneous edges have been identified and marked for removal, remove these marked edges
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reduce edges, %u edges identified for removal",Num2Remove);
pDstEdge = m_pGraphOutEdges;
pEdge = pDstEdge;
NumRemoved = 0;
for(Idx = 0; Idx < m_UsedGraphOutEdges; Idx++,pEdge++)
	{
	if(!pEdge->bRemove)
		{
		if(NumRemoved)
			*pDstEdge = *pEdge;
		pDstEdge += 1;
		}
	else
		NumRemoved +=1;
	}
m_UsedGraphOutEdges -= NumRemoved;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reduce edges completed, removed %d edges, %u edges retained",NumRemoved,m_UsedGraphOutEdges);
return(NumRemoved);
}

#ifdef USETHISCODE

tsGraphOutEdge *			// ptr to matching SeqID1.DnSeqID, or NULL if unable to locate				
CAssembGraph::LocateSeqIDDnSeqID(tSeqID SeqID,			// match this SeqID1
				  tSeqID DnSeqID)				// with this DnSeqID
{
tsGraphOutEdge *pEl2;
int CmpRslt;
UINT32 TargPsn;
UINT32 NodeLo;				
UINT32 NodeHi;				
if(m_NodeSortMode != eNSMSeqIDDnSeqID)		// nodes must be sorted by SeqID1.DnSeqID ascending
	return(NULL);

NodeLo = 0;
NodeHi = m_UsedGraphNodes - 1;
do {
	TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;
	pEl2 = &m_pGraphOutEdges[TargPsn];

	if(SeqID == pEl2->SeqID && DnSeqID == pEl2->DnSeqID)
		return(pEl2);

	if(SeqID > pEl2->SeqID)
		CmpRslt = 1;
	else
		{
		if(SeqID < pEl2->SeqID)
			CmpRslt = -1;
		else
			{
			if(DnSeqID > pEl2->DnSeqID)
				CmpRslt = 1;
			else
				CmpRslt = -1;
			}
		}

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		NodeHi = TargPsn - 1;
		}
	else
		NodeLo = TargPsn+1;
	}
while(NodeHi >= NodeLo);

return(0);	// unable to locate any instance of SeqID1.DnSeqID
}

tsGraphOutEdge *			// ptr to matching DnSeqID.SeqID1, or NULL if unable to locate				
CAssembGraph::LocateDnSeqIDSeq(tSeqID DnSeqID,	// match this DnSeqID
				  tSeqID SeqID)				// with this SeqID
{
tsGraphOutEdge *pEl2;
int CmpRslt;
UINT32 TargPsn;
UINT32 NodeLo;				
UINT32 NodeHi;				
if(m_NodeSortMode != eNSMDnSeqIDSeqID)		// nodes must be sorted by DnSeqID.SeqID1 ascending
	return(NULL);

NodeLo = 0;
NodeHi = m_UsedGraphNodes - 1;
do {
	TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;
	pEl2 = &m_pGraphOutEdges[TargPsn];

	if(DnSeqID == pEl2->DnSeqID && SeqID == pEl2->SeqID)
		return(pEl2);

	if(DnSeqID > pEl2->DnSeqID)
		CmpRslt = 1;
	else
		{
		if(DnSeqID < pEl2->DnSeqID)
			CmpRslt = -1;
		else
			{
			if(DnSeqID > pEl2->DnSeqID)
				CmpRslt = 1;
			else
				CmpRslt = -1;
			}
		}

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		NodeHi = TargPsn - 1;
		}
	else
		NodeLo = TargPsn+1;
	}
while(NodeHi >= NodeLo);

return(0);	// unable to locate any instance of DnSeqID.SeqID1
}

UINT32			// index+1 in m_pGraphOutEdges of first matching SeqID1, or 0 if non matching				
CAssembGraph::LocateFirstSeqID(tSeqID SeqID)			// find first matching 
{
tsGraphOutEdge *pEl2;
int CmpRslt;
UINT32 Mark;
UINT32 TargPsn;
UINT32 NodeLo = 0;
UINT32 NodeHi = m_UsedGraphNodes-1;

if(m_NodeSortMode != eNSMSeqIDDnSeqID)		// nodes must be sorted by SeqID1 ascending
	return(0);

do {
	TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;
	pEl2 = &m_pGraphOutEdges[TargPsn];

	if(SeqID > pEl2->SeqID)
		CmpRslt = 1;
	else
		if(SeqID < pEl2->SeqID)
			CmpRslt = -1;
		else
			CmpRslt = 0;

	if(!CmpRslt)	// if a match then may not be the lowest indexed match
		{
		if(TargPsn == 0 || NodeLo == TargPsn) // check if already lowest
			return(TargPsn + 1);
		// iterate until lowest located
		while(1) {
			if(CmpRslt == 0)
				{
				Mark = TargPsn;
				if(Mark == 0)
					return(Mark+1);
				NodeHi = TargPsn - 1;
				}	
			TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;

			pEl2 = &m_pGraphOutEdges[TargPsn];
			if(SeqID > pEl2->SeqID)
				CmpRslt = 1;
			else
				if(SeqID < pEl2->SeqID)
					CmpRslt = -1;
				else
					{
					CmpRslt = 0;
					continue;
					}
			NodeLo = TargPsn + 1;
			if(NodeLo == Mark)
				return(Mark+1);
			}
		}

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		NodeHi = TargPsn - 1;
		}
	else
		NodeLo = TargPsn+1;
	}
while(NodeHi >= NodeLo);

return(0);	// unable to locate any instance of SeqID
}


UINT32			// index+1 in m_pGraphOutEdges of first matching DnSeqID, or 0 if non matching				
CAssembGraph::LocateFirstDnSeqID(tSeqID DnSeqID)			// find first matching 
{
tsGraphOutEdge *pEl2;
int CmpRslt;
UINT32 Mark;
UINT32 TargPsn;
UINT32 NodeLo = 0;
UINT32 NodeHi = m_UsedGraphNodes-1;

if(!(m_NodeSortMode == eNSMDnSeqIDSeqID || m_NodeSortMode == eNSMDnSeqIDOfs))		// nodes must be sorted by DnSeqID ascending
	return(0);

do {
	TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;
	pEl2 = &m_pGraphOutEdges[TargPsn];

	if(DnSeqID > pEl2->DnSeqID)
		CmpRslt = 1;
	else
		if(DnSeqID < pEl2->DnSeqID)
			CmpRslt = -1;
		else
			CmpRslt = 0;

	if(!CmpRslt)	// if a match then may not be the lowest indexed match
		{
		if(TargPsn == 0 || NodeLo == TargPsn) // check if already lowest
			return(TargPsn + 1);
		// iterate until lowest located
		while(1) {
			if(CmpRslt == 0)
				{
				Mark = TargPsn;
				if(Mark == 0)
					return(Mark+1);
				NodeHi = TargPsn - 1;
				}	
			TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;

			pEl2 = &m_pGraphOutEdges[TargPsn];
			if(DnSeqID > pEl2->DnSeqID)
				CmpRslt = 1;
			else
				if(DnSeqID < pEl2->DnSeqID)
					CmpRslt = -1;
				else
					{
					CmpRslt = 0;
					continue;
					}
			NodeLo = TargPsn + 1;
			if(NodeLo == Mark)
				return(Mark+1);
			}
		}

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		NodeHi = TargPsn - 1;
		}
	else
		NodeLo = TargPsn+1;
	}
while(NodeHi >= NodeLo);

return(0);	// unable to locate any instance of SeqID
}
#endif


int  // function for sorting graph vertices on SeqID.VertexID ascending 
CAssembGraph::SortVerticesSeqID(const void *arg1, const void *arg2)
{
tsGraphVertex *pVertex1 = (tsGraphVertex *)arg1;
tsGraphVertex *pVertex2 = (tsGraphVertex *)arg2;

if(pVertex1->SeqID < pVertex2->SeqID)
	return(-1);
else
	if(pVertex1->SeqID > pVertex2->SeqID)
		return(1);
if(pVertex1->VertexID < pVertex2->VertexID)
	return(-1);
else
	if(pVertex1->VertexID > pVertex2->VertexID)
		return(1);
return(0);
}

int  // function for sorting graph vertices on DiscGraphID.VertexID ascending 
CAssembGraph::SortVerticesDiscGraphID(const void *arg1, const void *arg2)
{
tsGraphVertex *pVertex1 = (tsGraphVertex *)arg1;
tsGraphVertex *pVertex2 = (tsGraphVertex *)arg2;

if(pVertex1->DiscGraphID < pVertex2->DiscGraphID)
	return(-1);
else
	if(pVertex1->DiscGraphID > pVertex2->DiscGraphID)
		return(1);
if(pVertex1->VertexID < pVertex2->VertexID)
	return(-1);
else
	if(pVertex1->VertexID > pVertex2->VertexID)
		return(1);
return(0);
}

int  // function for sorting graph vertices on VertexID.SeqID ascending 
CAssembGraph::SortVerticesVertexID(const void *arg1, const void *arg2)
{
tsGraphVertex *pVertex1 = (tsGraphVertex *)arg1;
tsGraphVertex *pVertex2 = (tsGraphVertex *)arg2;

if(pVertex1->VertexID < pVertex2->VertexID)
	return(-1);
else
	if(pVertex1->VertexID > pVertex2->VertexID)
		return(1);
if(pVertex1->SeqID < pVertex2->SeqID)
	return(-1);
else
	if(pVertex1->SeqID > pVertex2->SeqID)
		return(1);
return(0);
}

int  // function for sorting outgoing edges on FromVertexID.ToVertexID ascending 
CAssembGraph::SortOutEdgeFromVertexID(const void *arg1, const void *arg2)
{
tsGraphOutEdge *pEdge1 = (tsGraphOutEdge *)arg1;
tsGraphOutEdge *pEdge2 = (tsGraphOutEdge *)arg2;

if(pEdge1->FromVertexID < pEdge2->FromVertexID)
	return(-1);
else
	if(pEdge1->FromVertexID > pEdge2->FromVertexID)
		return(1);
if(pEdge1->ToVertexID < pEdge2->ToVertexID)
	return(-1);
else
	if(pEdge1->ToVertexID > pEdge2->ToVertexID)
		return(1);
return(0);
}

int  // function for sorting outgoing edges on FromVertexID.SeqOfs ascending 
CAssembGraph::SortOutEdgeFromVertexIDSeqOfs(const void *arg1, const void *arg2)
{
tsGraphOutEdge *pEdge1 = (tsGraphOutEdge *)arg1;
tsGraphOutEdge *pEdge2 = (tsGraphOutEdge *)arg2;

if(pEdge1->FromVertexID < pEdge2->FromVertexID)
	return(-1);
else
	if(pEdge1->FromVertexID > pEdge2->FromVertexID)
		return(1);

if(pEdge1->SeqOfs < pEdge2->SeqOfs)
	return(-1);
else
	if(pEdge1->SeqOfs > pEdge2->SeqOfs)
		return(1);
return(0);
}

int  // function for sorting outgoing edges on ToVertexID.SeqOfs ascending 
CAssembGraph::SortOutEdgeToVertexIDSeqOfs(const void *arg1, const void *arg2)
{
tsGraphOutEdge *pEdge1 = (tsGraphOutEdge *)arg1;
tsGraphOutEdge *pEdge2 = (tsGraphOutEdge *)arg2;

if(pEdge1->ToVertexID < pEdge2->ToVertexID)
	return(-1);
else
	if(pEdge1->ToVertexID > pEdge2->ToVertexID)
		return(1);

if(pEdge1->SeqOfs < pEdge2->SeqOfs)
	return(-1);
else
	if(pEdge1->SeqOfs > pEdge2->SeqOfs)
		return(1);
return(0);
}

int  // function for sorting incomimg edges on ToVertexID.FromVertexID ascending 
CAssembGraph::SortInEdgesToVertexID(const void *arg1, const void *arg2)
{
tsGraphOutEdge *pEdge1 = &pStaticGraphOutEdges[(*(tEdgeID *)arg1)-1];
tsGraphOutEdge *pEdge2 = &pStaticGraphOutEdges[(*(tEdgeID *)arg2)-1];

if(pEdge1->ToVertexID < pEdge2->ToVertexID)
	return(-1);
else
	if(pEdge1->ToVertexID > pEdge2->ToVertexID)
		return(1);
if(pEdge1->FromVertexID < pEdge2->FromVertexID)
	return(-1);
else
	if(pEdge1->FromVertexID > pEdge2->FromVertexID)
		return(1);
return(0);
}


int  // function for sorting components on NumVertices decending 
CAssembGraph::SortComponentNumVertices(const void *arg1, const void *arg2)
{
tsComponent *pComp1 = (tsComponent *)arg1;
tsComponent *pComp2 = (tsComponent *)arg2;

if(pComp1->NumVertices > pComp2->NumVertices)
	return(-1);
else
	if(pComp1->NumVertices < pComp2->NumVertices)
		return(1);
if(pComp1->ComponentID < pComp2->ComponentID)
	return(-1);
else
	if(pComp1->ComponentID > pComp2->ComponentID)
		return(1);
return(0);
}
