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

#include "pacbiocommon.h"
#include "SeqStore.h"
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
m_pPathTraceBacks = NULL;
m_bMutexesCreated = false;
m_CASSerialise = 0;
m_CASLock = 0;
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
m_CASSerialise = 0;
m_CASLock = 0;

m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CAssembGraph::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
m_bMutexesCreated = false;
}


void
CAssembGraph::AcquireCASSerialise(void)
{
int SpinCnt = 5000;
int BackoffMS = 5;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_CASSerialise,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 1000;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#else
while(__sync_val_compare_and_swap(&m_CASSerialise,0,1)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 1000;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#endif
}

void
CAssembGraph::ReleaseCASSerialise(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASSerialise,0,1);
#else
__sync_val_compare_and_swap(&m_CASSerialise,1,0);
#endif
}


void
CAssembGraph::AcquireCASLock(void)
{
int SpinCnt = 5000;
int BackoffMS = 5;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_CASLock,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 1000;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#else
while(__sync_val_compare_and_swap(&m_CASLock,0,1)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 1000;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#endif
}

void
CAssembGraph::ReleaseCASLock(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASLock,0,1);
#else
__sync_val_compare_and_swap(&m_CASLock,1,0);
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


if(m_pPathTraceBacks != NULL)
	{
#ifdef _WIN32
	free(m_pPathTraceBacks);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pPathTraceBacks != MAP_FAILED)
		munmap(m_pPathTraceBacks,m_AllocdTraceBacks  * sizeof(tsPathTraceBack));
#endif	
	m_pPathTraceBacks = NULL;
	}

m_AllocGraphVertices = 0;
m_AllocGraphOutEdges = 0;
m_AllocGraphInEdges = 0;

m_AllocTransitStack = 0;
m_UsedGraphVertices = 0;
m_UsedGraphOutEdges = 0;
m_UsedGraphInEdges = 0;

m_AllocdTraceBacks = 0;
m_UsedTraceBacks = 0;

m_CurTransitDepth = 0;
m_MaxTransitDepth = 0;

m_VerticesSortOrder = eVSOUnsorted;
m_bOutEdgeSorted = false;
m_bVertexEdgeSet = false;
m_bInEdgeSorted = false;

m_bReduceEdges = true;
m_NumReducts = 0;

m_NumComponents = 0;	
m_AllocComponents= 0;

m_MinScaffScoreThres = cDfltMinAcceptScoreThres;

m_NumDiscRemaps = 0;
m_bTerminate = false;
DeleteMutexes();
}

teBSFrsltCodes 
CAssembGraph::Init(int ScaffScoreThres,		// accepted edges must be of at least this overlap score
					int MaxThreads)		// initialise with maxThreads threads, if maxThreads == 0 then number of threads not set	
{
if(MaxThreads < 0 || MaxThreads > cMaxWorkerThreads ||
		ScaffScoreThres < 0)
		return(eBSFerrParams);
m_MinScaffScoreThres = ScaffScoreThres;
Reset();
if(MaxThreads > 0)
	{
	m_NumThreads = MaxThreads;
	m_MTqsort.SetMaxThreads(MaxThreads);
	}
CreateMutexes();
return(eBSFSuccess);
}


UINT32									// returns total number of edges, including these edges if accepted, thus far accepted; if > cMaxValidID used as processing error indicator, cast to teBSFrsltCodes for actual error
CAssembGraph::AddEdges(UINT32 NumSeqs,			// number of overlapped sequences
		tsOverlappedSeq *pOverlappedSeqs)		// overlapped sequences
{
UINT32 Idx;
tsGraphOutEdge *pOutEdge;
tsGraphOutEdge *pPrevEdge;
tsOverlappedSeq *pOverlap;
size_t AllocMem;
tSeqID CurOverlapSeqID;
tSeqID CurOverlappedSeqID;
tsGraphVertex *pFromVertex;
tsGraphVertex *pToVertex;
tVertID OverlappingVertexID;
tVertID OverlappedVertexID;

if(NumSeqs == 0 || pOverlappedSeqs == NULL)
	return((UINT32)eBSFerrParams);
AcquireCASSerialise();
if(m_bTerminate)		// check if should early exit
	{
	ReleaseCASSerialise();
	return((UINT32)eBSFerrInternal);
	}
if((m_UsedGraphOutEdges + NumSeqs) > cMaxValidID)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: too many edges, can only accept at most %u",cMaxValidID);
	m_bTerminate = true;
	ReleaseCASSerialise();
	return((UINT32)eBSFerrMaxEntries);
	}

if(m_VerticesSortOrder != eVSOSeqID)		// vertices need to have been sorted by SeqID ascending so can check that edge referenced sequences are known to this class
	{
	if(m_UsedGraphVertices > 1)
		{
		pStaticGraphVertices = m_pGraphVertices;
		m_MTqsort.qsort(pStaticGraphVertices,m_UsedGraphVertices,sizeof(tsGraphVertex),SortVerticesSeqID);
		}
	m_VerticesSortOrder = eVSOSeqID;
	}

// ensure m_pGraphOutEdges allocated to hold additional NumSeqs
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
		ReleaseCASSerialise();
		return((UINT32)eBSFerrMem);
		}
#else
	m_pGraphOutEdges = (tsGraphOutEdge *)mmap(NULL,AllocMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pGraphOutEdges == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: graph forward edge (%d bytes per edge) allocation of %d edges failed - %s",
								(int)sizeof(tsGraphOutEdge),cInitialAllocVertices,strerror(errno));
		m_bTerminate = true;
		ReleaseCASSerialise();
		return(eBSFerrMem);
		}
#endif
	m_AllocGraphOutEdges = cInitialAllocVertices;
	m_UsedGraphOutEdges = 0;
	}
else
	{
	if((m_UsedGraphOutEdges + (2 * NumSeqs) + 10) >= m_AllocGraphOutEdges) // need to add additional forward graph edges, 10 is simply to allow a small margin of error
		{
		tsGraphOutEdge *pTmp;
		UINT64 AllocEdges;
		UINT64 ReallocEdges;

		ReallocEdges =  (UINT64)((double)m_AllocGraphOutEdges * cReallocEdges);
		AllocEdges = (UINT64)m_AllocGraphOutEdges + ReallocEdges;
		if(AllocEdges > cMaxGraphEdges)
			AllocEdges = cMaxGraphEdges;
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
			ReleaseCASSerialise();
			return((UINT32)eBSFerrMem);
			}
		m_AllocGraphOutEdges += (UINT32)ReallocEdges;
		m_pGraphOutEdges = pTmp;
		}
	}

CurOverlapSeqID = 0;
CurOverlappedSeqID = 0;
pOverlap = pOverlappedSeqs;
pOutEdge = &m_pGraphOutEdges[m_UsedGraphOutEdges];
memset(pOutEdge,0,NumSeqs * 2 * sizeof(tsGraphOutEdge));
for(Idx = 0; Idx < NumSeqs; Idx++,pOverlap++)
	{
	if(pOverlap->Score < m_MinScaffScoreThres || pOverlap->OverlapClass == eOLCOverlapping)   // only accepting if previously classified as being overlapping and score at least threshold
		continue;

	// ensure up front that the sequence identifiers can be have been associated to the correct vertex
	if(pOverlap->FromSeqID != CurOverlapSeqID)
		{
		CurOverlapSeqID = pOverlap->FromSeqID;
		pFromVertex = LocateVertexSeqID(CurOverlapSeqID);
		if(pFromVertex == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: unable to locate vertex for SeqID %u", CurOverlapSeqID);
			m_bTerminate = true;
			ReleaseCASSerialise();
			return((UINT32)eBSFerrParams);
			}
		OverlappingVertexID = pFromVertex->VertexID;
		}

	// check that from start/end loci are within 'from' sequence length; note that start/end loci are 0 based and min alignment length of 100bp 
	if((int)pOverlap->FromSeq5Ofs < 0 || (pOverlap->FromSeq5Ofs + 100) > pFromVertex->SeqLen ||
		pOverlap->FromSeq3Ofs < (pOverlap->FromSeq5Ofs + 99) || pOverlap->FromSeq3Ofs >= pFromVertex->SeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: Inconsistencies in overlap offsets for SeqID %u", CurOverlapSeqID);
		m_bTerminate = true;
		ReleaseCASSerialise();
		return((UINT32)eBSFerrParams);
		}

	if(pOverlap->ToSeqID != CurOverlappedSeqID)
		{
		CurOverlappedSeqID = pOverlap->ToSeqID;
		pToVertex = LocateVertexSeqID(CurOverlappedSeqID);
		if(pToVertex == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: unable to locate vertex for SeqID %u", CurOverlappedSeqID);
			m_bTerminate = true;
			ReleaseCASSerialise();
			return((UINT32)eBSFerrParams);
			}
		OverlappedVertexID = pToVertex->VertexID;
		}

	// check that start/end loci are within 'to' sequence length; note that start/end loci are 0 based and min alignment length of 100bp
	if((int)pOverlap->ToSeq5Ofs < 0 || (pOverlap->ToSeq5Ofs + 100) > pToVertex->SeqLen ||
		pOverlap->ToSeq3Ofs < (pOverlap->ToSeq5Ofs + 99) || pOverlap->ToSeq3Ofs >= pToVertex->SeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: Inconsistencies in overlap offsets for SeqID %u", CurOverlappedSeqID);
		m_bTerminate = true;
		ReleaseCASSerialise();
		return((UINT32)eBSFerrParams);
		}

	pOutEdge->FromVertexID = OverlappingVertexID;
	pOutEdge->ToVertexID = OverlappedVertexID;
	pOutEdge->FromSeqLen = pOverlap->FromSeqLen;
	pOutEdge->ToSeqLen = pOverlap->ToSeqLen;
	pOutEdge->FromSeq5Ofs = pOverlap->FromSeq5Ofs;
	pOutEdge->FromSeq3Ofs = pOverlap->FromSeq3Ofs;
	pOutEdge->ToSeq5Ofs = pOverlap->ToSeq5Ofs;
	pOutEdge->ToSeq3Ofs = pOverlap->ToSeq3Ofs;
	pOutEdge->OverlapClass = pOverlap->OverlapClass;
	pOutEdge->Score = pOverlap->Score;
	pOutEdge->ScoreAlignLen = pOverlap->ScoreAlignLen;
	pOutEdge->flgContains = pOverlap->OverlapClass == eOLCcontains ? 1 : 0;
	pOutEdge->flgContained = pOverlap->OverlapClass == eOLCcontained ? 1 : 0;
	pOutEdge->flgArtefact =  pOverlap->OverlapClass == eOLCartefact ? 1 : 0;
	pOutEdge->flgsOvlSense = pOverlap->bAntisense ? 2 : 0;;
	pPrevEdge = pOutEdge;
	pOutEdge += 1;
	m_UsedGraphOutEdges += 1;
	// add back edge from the OverlappedVertexID back to the OverlappingVertexID in case caller doesn't provide one
	// will mark as being an inferred edge so later can be removed if caller does provide it
	memcpy(pOutEdge,pPrevEdge,sizeof(tsGraphOutEdge));
	pOutEdge->FromVertexID = OverlappedVertexID;
	pOutEdge->ToVertexID = OverlappingVertexID;
	pOutEdge->FromSeq5Ofs = pOverlap->ToSeq5Ofs;
	pOutEdge->FromSeq3Ofs = pOverlap->ToSeq3Ofs;
	pOutEdge->ToSeq5Ofs = pOverlap->FromSeq5Ofs;
	pOutEdge->ToSeq3Ofs = pOverlap->FromSeq3Ofs;
	pOutEdge->flgsOvlSense = pOverlap->bAntisense ? 1 : 0;

	pOutEdge->flgInfBackEdge = 1;	// mark as inferenced so can be later removed if user supplies actual aligned edge
	pOutEdge += 1;
	m_UsedGraphOutEdges += 1;
	}
m_bReduceEdges = true;
m_bOutEdgeSorted = false;
m_bInEdgeSorted = false;
ReleaseCASSerialise();
return(m_UsedGraphOutEdges);
}


UINT32										// returns total number of edges, including this edge if accepted, thus far accepted; if > cMaxValidID used as processing error indicator, cast to teBSFrsltCodes for actual error
CAssembGraph::AddEdge(tSeqID FromSeqID,		// identifies the 'From' or overlapping sequence
				tSeqID ToSeqID,				// identifies the 'To' overlapped sequence
				UINT32 FromSeqLen,			// 'From' sequence length is this length
				UINT32 ToSeqLen,			// 'To' sequence length is this length
				UINT32 Score,				// score associated with this overlap, higher scores represent higher confidence in the overlap
				UINT32 ScoreAlignLen,		// scored for this alignment length (ProbeAlignLen + TargAlignLen) / 2
				UINT32 FromSeq5Ofs,			// overlap of FromSeqID onto ToSeqID starts at this base relative to FromSeqID 5' start
				UINT32 FromSeq3Ofs,			// overlap of FromSeqID onto ToSeqID ends at this base relative to FromSeqID 5' start
				UINT32 ToSeq5Ofs,			// overlap onto ToSeqID from FromSeqID starts at this base relative to ToSeqID 5' start
				UINT32 ToSeq3Ofs,			// overlap onto ToSeqID from FromSeqID ends at this base relative to ToSeqID 5' start
				eOverlapClass OverlapClass,	// classification of overlap from FromSeqID onto ToSeqID, note that classification must be eOLCOverlapping
				bool bAntisense)			// false: 'From' sense overlaps 'To' sense, true: 'From' sense overlaps 'To' antisense 
{
size_t AllocMem;
tsGraphOutEdge *pOutEdge;

tsGraphOutEdge *pPrevEdge;
tsGraphVertex *pFromVertex;
tsGraphVertex *pToVertex;
tVertID OverlappingVertexID;
tVertID OverlappedVertexID;

if(Score < m_MinScaffScoreThres || OverlapClass != eOLCOverlapping)			// only accepting edges which have been classified by caller as overlapping and score at least the minimum threshold
	return(m_UsedGraphOutEdges);			

if(FromSeqID < 1 || ToSeqID < 1 || FromSeqID > cMaxValidID || ToSeqID > cMaxValidID)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge: Sequence identifiers must be in range 1..%u",cMaxValidID);
	m_bTerminate = true;
	return((UINT32)eBSFerrMaxEntries);
	}

AcquireCASSerialise();
if(m_bTerminate)		// check if should early exit
	{
	ReleaseCASSerialise();
	return((UINT32)eBSErrSession);
	}

if(m_UsedGraphOutEdges >= cMaxValidID)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge: too many edges, can only accept at most %u",cMaxValidID);
	m_bTerminate = true;
	ReleaseCASSerialise();
	return((UINT32)eBSFerrMaxEntries);
	}

if(m_VerticesSortOrder != eVSOSeqID)		// vertices need to have been sorted by SeqID ascending so can check that edge referenced sequences are known to this class
	{                                     
	if(m_UsedGraphVertices > 1)
		{
		pStaticGraphVertices = m_pGraphVertices;
		m_MTqsort.qsort(pStaticGraphVertices,m_UsedGraphVertices,sizeof(tsGraphVertex),SortVerticesSeqID);
		}
	m_VerticesSortOrder = eVSOSeqID;
	}

// ensure up front that the sequence identifiers are known to have been associated with a vertex and that the lengths match
pFromVertex = LocateVertexSeqID(FromSeqID);
if(pFromVertex == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge: unable to locate vertex for SeqID %u", FromSeqID);
	m_bTerminate = true;
	ReleaseCASSerialise();
	return((UINT32)eBSFerrParams);
	}

if(pFromVertex->SeqLen != FromSeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge: FromSeqLen: %u not equal to vertex sequence length: %u for FromSeqID: %u",FromSeqLen, pFromVertex->SeqLen,FromSeqID);
	m_bTerminate = true;
	ReleaseCASSerialise();
	return((UINT32)eBSFerrParams);
	}

// check that from start/end loci are within 'from' sequence length; note that start/end loci are 0 based and min alignment length of 100bp 
if((int)FromSeq5Ofs < 0 || (FromSeq5Ofs + 100) > pFromVertex->SeqLen ||
	FromSeq3Ofs < (FromSeq5Ofs + 99) || FromSeq3Ofs >= pFromVertex->SeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: Inconsistencies in overlap offsets for SeqID %u", FromSeqID);
	m_bTerminate = true;
	ReleaseCASSerialise();
	return((UINT32)eBSFerrParams);
	}
OverlappingVertexID = pFromVertex->VertexID;

pToVertex = LocateVertexSeqID(ToSeqID);
if(pToVertex == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge: unable to locate vertex for SeqID %u", ToSeqID);
	m_bTerminate = true;
	ReleaseCASSerialise();
	return((UINT32)eBSFerrParams);
	}

if(pToVertex->SeqLen != ToSeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge: ToSeqLen: %u not equal to vertex sequence length: %u for ToSeqID: %u",ToSeqLen, pToVertex->SeqLen,ToSeqID);
	m_bTerminate = true;
	ReleaseCASSerialise();
	return((UINT32)eBSFerrParams);
	}

// check that start/end loci are within 'to' sequence length; note that start/end loci are 0 based and min alignment length of 100bp
if((int)ToSeq5Ofs < 0 || (ToSeq5Ofs + 100) > pToVertex->SeqLen ||
	ToSeq3Ofs < (ToSeq5Ofs + 99) || ToSeq3Ofs >= pToVertex->SeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdges: Inconsistencies in overlap offsets for SeqID %u", ToSeqID);
	m_bTerminate = true;
	ReleaseCASSerialise();
	return((UINT32)eBSFerrParams);
	}
OverlappedVertexID = pToVertex->VertexID;

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
		ReleaseCASSerialise();
		return((UINT32)eBSFerrMem);
		}
#else
	m_pGraphOutEdges = (tsGraphOutEdge *)mmap(NULL,AllocMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pGraphOutEdges == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEdge: graph forward edge (%d bytes per edge) allocation of %d edges failed - %s",
								(int)sizeof(tsGraphOutEdge),cInitialAllocVertices,strerror(errno));
		m_bTerminate = true;
		ReleaseCASSerialise();
		return((UINT32)eBSFerrMem);
		}
#endif
	m_AllocGraphOutEdges = cInitialAllocVertices;
	m_UsedGraphOutEdges = 0;
	}
else
	{
	if((m_UsedGraphOutEdges + 10) > m_AllocGraphOutEdges) // need to add additional forward graph edges, 10 is simply to allow a small margin of error
		{
		tsGraphOutEdge *pTmp;
		UINT64 AllocEdges;
		UINT64 ReallocEdges =  (UINT64)((double)m_AllocGraphOutEdges * cReallocEdges);
		AllocEdges = (UINT64)m_AllocGraphOutEdges + ReallocEdges;
		if(AllocEdges > cMaxGraphEdges)
			AllocEdges = cMaxGraphEdges;
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
			ReleaseCASSerialise();
			return((UINT32)eBSFerrMem);
			}
		m_AllocGraphOutEdges += (UINT32)ReallocEdges;
		m_pGraphOutEdges = pTmp;
		}
	}

// at least four unused outgoing edges available
pOutEdge = &m_pGraphOutEdges[m_UsedGraphOutEdges++];
memset(pOutEdge,0,sizeof(tsGraphOutEdge));
m_bOutEdgeSorted = false;
m_bInEdgeSorted = false;
pOutEdge->FromVertexID = OverlappingVertexID;
pOutEdge->ToVertexID = OverlappedVertexID;
pOutEdge->FromSeqLen = FromSeqLen;
pOutEdge->ToSeqLen = ToSeqLen;
pOutEdge->Score = Score;
pOutEdge->ScoreAlignLen = ScoreAlignLen;
pOutEdge->FromSeq5Ofs = FromSeq5Ofs;
pOutEdge->FromSeq3Ofs = FromSeq3Ofs;
pOutEdge->ToSeq5Ofs = ToSeq5Ofs;
pOutEdge->ToSeq3Ofs = ToSeq3Ofs;
pOutEdge->OverlapClass = OverlapClass;
pOutEdge->flgContains = OverlapClass == eOLCcontains ? 1 : 0;  // redundant as only overlapping, non-contained, edges accepted
pOutEdge->flgContained = OverlapClass == eOLCcontained ? 1 : 0;
pOutEdge->flgArtefact = OverlapClass == eOLCartefact ? 1 : 0;
pOutEdge->flgsOvlSense = bAntisense ? 2 : 0;
pOutEdge->flgInfBackEdge = 0;

pPrevEdge = pOutEdge;

// add back edge from the OverlappedVertexID back to the OverlappingVertexID in case caller doesn't provide one
// will mark as being an inferred edge so later can be removed if caller does provide it
pOutEdge = &m_pGraphOutEdges[m_UsedGraphOutEdges++];
memcpy(pOutEdge,pPrevEdge,sizeof(tsGraphOutEdge));
pOutEdge->FromVertexID = OverlappedVertexID;
pOutEdge->ToVertexID = OverlappingVertexID;
pOutEdge->FromSeqLen = ToSeqLen;
pOutEdge->ToSeqLen = FromSeqLen;
pOutEdge->FromSeq5Ofs = ToSeq5Ofs;
pOutEdge->FromSeq3Ofs = ToSeq3Ofs;
pOutEdge->ToSeq5Ofs = FromSeq5Ofs;
pOutEdge->ToSeq3Ofs = FromSeq3Ofs;
pOutEdge->flgsOvlSense = bAntisense ? 1 : 0;
pOutEdge->flgInfBackEdge = 1;	// mark as inferenced so can be later removed if user supplies actual aligned edge


// add edge as antisense overlapping antisense if was originally sense overlapping sense
// will mark as being an inferred edges so later can be removed if caller does provide it
if(!bAntisense)
	{
	tsGraphOutEdge *pOutEdgeA;

	pOutEdgeA = &m_pGraphOutEdges[m_UsedGraphOutEdges++];
	memcpy(pOutEdgeA,pPrevEdge,sizeof(tsGraphOutEdge));
	pOutEdgeA->FromSeq5Ofs = FromSeqLen - (FromSeq3Ofs + 1);
	pOutEdgeA->FromSeq3Ofs = FromSeqLen - (FromSeq5Ofs + 1);
	pOutEdgeA->ToSeq5Ofs = ToSeqLen - (ToSeq3Ofs + 1);
	pOutEdgeA->ToSeq3Ofs = ToSeqLen - (ToSeq5Ofs + 1);
	pOutEdgeA->flgsOvlSense = 3;
	pOutEdgeA->flgInfBackEdge = 1;	// mark as inferenced so can be later removed if user supplies actual aligned edge


	pOutEdgeA = &m_pGraphOutEdges[m_UsedGraphOutEdges++];
	memcpy(pOutEdgeA,pOutEdge,sizeof(tsGraphOutEdge));
	pOutEdgeA->FromSeq5Ofs = pOutEdge->FromSeqLen - (pOutEdge->FromSeq3Ofs + 1);
	pOutEdgeA->FromSeq3Ofs = pOutEdge->FromSeqLen - (pOutEdge->FromSeq5Ofs + 1);
	pOutEdgeA->ToSeq5Ofs = pOutEdge->ToSeqLen - (pOutEdge->ToSeq3Ofs + 1);
	pOutEdgeA->ToSeq3Ofs = pOutEdge->ToSeqLen - (pOutEdge->ToSeq5Ofs + 1);
	pOutEdgeA->flgsOvlSense = 3;
	pOutEdgeA->flgInfBackEdge = 1;	// mark as inferenced so can be later removed if user supplies actual aligned edge
	}

m_bReduceEdges = true;
ReleaseCASSerialise();
return(m_UsedGraphOutEdges);
}

UINT32										// returned vertex identifier
CAssembGraph::AddVertex(UINT32 SeqLen,		// sequence length
						tSeqID SeqID)		// allocates and initialises a new graph vertex which references this sequence
{
size_t AllocMem;
tVertID VertexID;
tsGraphVertex *pVertex;

AcquireCASSerialise();
if(m_UsedGraphVertices >= cMaxValidID)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeqOverlap: too many vertices, can only accept at most %u",cMaxValidID);
	ReleaseCASSerialise();
	Reset();
	return((UINT32)eBSFerrMaxEntries);
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
		ReleaseCASSerialise();
		Reset();
		return((UINT32)eBSFerrMem);
		}
#else
	m_pGraphVertices = (tsGraphVertex *)mmap(NULL,AllocMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pGraphVertices == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddVertex: graph vertices (%d bytes per vertex) allocation of %u vertices failed - %s",
													(int)sizeof(tsGraphVertex),cInitialAllocVertices,strerror(errno));
		m_pGraphVertices = NULL;
		ReleaseCASSerialise();
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
			ReleaseCASSerialise();
			return((UINT32)eBSFerrMem);
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
pVertex->SeqID = SeqID;
pVertex->SeqLen = SeqLen;
ReleaseCASSerialise();
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

return(NULL);	// unable to locate any instance of VertexID
}


// when all sequences have been associated to vertices then FinaliseVertices must be called which
// will firstly sort the vertices in SeqID ascending order and then remap the PEVertexIDs ready for edges between adjacent vertexes to be added with AddEdge()
UINT32									// 0 if errors else number of vertices finalised
CAssembGraph::FinaliseVertices(void)
{
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
UINT32 NumRemoved;
UINT32 Num2Remove;
tsGraphOutEdge *pOutEdge;
tsGraphOutEdge *pEdge;
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

// ensure vertices sorted by vertex identifier ascending
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

// mark for removal all back edges which were originally inferred but for which an actual alignment edge exists
pEdge = m_pGraphOutEdges;
pOutEdge = pEdge+1;
Num2Remove = 0;
for(EdgeIdx = 0; EdgeIdx < m_UsedGraphOutEdges - 1; EdgeIdx++,pOutEdge++,pEdge++)
	{
	if(pEdge->FromVertexID == pOutEdge->FromVertexID &&
	   pEdge->ToVertexID == pOutEdge->ToVertexID &&
	   pEdge->flgsOvlSense == pOutEdge->flgsOvlSense &&
		pOutEdge->flgInfBackEdge == 1)
			{
			pOutEdge->flgRemove = 1;
			Num2Remove += 1;
			}
	}
if(Num2Remove)
	{
	pOutEdge = m_pGraphOutEdges;
	pEdge = pOutEdge;
	NumRemoved = 0;
	for(EdgeIdx = 0; EdgeIdx < m_UsedGraphOutEdges; EdgeIdx++,pEdge++)
		{
		if(!pEdge->flgRemove)
			{
			if(NumRemoved)
				*pOutEdge = *pEdge;
			pOutEdge += 1;
			}
		else
			NumRemoved +=1;
		}
	m_UsedGraphOutEdges -= NumRemoved;
	}

// ensure incoming edges are sorted ToVertexID.FromVertexID ascending order
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


if(m_bReduceEdges)
	{
	ReduceEdges();  // identify and remove any redundant edges, these could be parallel edges between same pairs of vertices
	m_bReduceEdges = false;
	if(m_UsedGraphOutEdges == 0)	// check if all edges were removed - very unlikely but this is a world of unlikeliness ....
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"FinaliseEdges: All edges were removed when reducing, nothing to assemble!");
		return(m_UsedGraphOutEdges);
		}
	}

if(!m_bVertexEdgeSet)
	{
	// iterate all outgoing edges and update vertices with the initial outgoing edge
	CurVertexID = 0;
	pOutEdge = m_pGraphOutEdges;
	for(EdgeIdx = 1; EdgeIdx <= m_UsedGraphOutEdges; EdgeIdx++,pOutEdge++)
		{
		if(CurVertexID == pOutEdge->FromVertexID)
			{
			if(pVertex->DegreeOut < cMaxEdges)
				pVertex->DegreeOut += 1;
			continue;
			}
		CurVertexID = pOutEdge->FromVertexID;
		pVertex = &m_pGraphVertices[CurVertexID-1];
		pVertex->OutEdgeID = EdgeIdx; 
		pVertex->DegreeOut = 1;
		}

	// iterate all incoming edges and update vertices with the initial incoming edge
	CurVertexID = 0;
	pInEdge = m_pGraphInEdges;
	for(EdgeIdx = 1; EdgeIdx <= m_UsedGraphInEdges; EdgeIdx++,pInEdge++)
		{
		pOutEdge = &m_pGraphOutEdges[*pInEdge - 1];
		if(CurVertexID == pOutEdge->ToVertexID)
			{
			if(pVertex->DegreeIn < cMaxEdges)
				pVertex->DegreeIn += 1;
			continue;
			}
		CurVertexID = pOutEdge->ToVertexID;
		pVertex = &m_pGraphVertices[CurVertexID-1];
		pVertex->InEdgeID = EdgeIdx; 
		pVertex->DegreeIn = 1;
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
AcquireCASSerialise();
Num = m_UsedGraphVertices;
ReleaseCASSerialise();
return(Num);
}

UINT32 
CAssembGraph::GetNumGraphOutEdges(void)		// returns current number of graph forward edges
{
UINT32 Num;
AcquireCASSerialise();
Num = m_UsedGraphOutEdges;
ReleaseCASSerialise();
return(Num);
}

UINT32 
CAssembGraph::GetNumReducts(void)		// returns current number of edge reductions
{
UINT32 Num;
AcquireCASSerialise();
Num = m_NumReducts;
ReleaseCASSerialise();
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
	pEdge->flgTravFwd = 0;
	pEdge->flgTravRev = 0;
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
	pVertex->ComponentID = 0;
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
	pVertex->DegreeOut = min(cMaxEdges,NumOutEdges);

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
	pVertex->DegreeIn = min(cMaxEdges,NumInEdges);

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


static UINT32 m_NumReplacedComponentIDs;
// IdentifyDisconnectedSubGraphs
// Within the graph there are likely to be many (could be millions) of completely disconnected subgraphs (components)
// These disconnected subgraphs have no sequences which overlay, or are overlaid by, sequences in any other subgraph 
// This function iterates all vertices of the graph and locates all other vertices which are connected to the original vertex and marks these
// as belonging to an disconnected subgraph
// all subgraphs or components are uniquely identified
// Currently this function is single threaded, could be worthwhile to multithread
// 
UINT32						// returned number of subgraphs identified
CAssembGraph::IdentifyDiscComponents(void)
{
UINT32 VertexIdx;
UINT32 NumVertices;
UINT32 MaxVertices;
UINT32 Cntr;
tComponentID CurComponentID;
tsGraphVertex *pVertex;
tsComponent *pComponent;

if(m_pGraphVertices == NULL || m_UsedGraphVertices < 1)
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
	m_NumComponents = 0;
	}
// determine, flag and report, on vertex degree of connectivity
VertexConnections();

// iterate vertices
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identifying disconnected graph components ...");
m_NumDiscRemaps = 0;
m_CurTransitDepth = 0;
CurComponentID = 0;
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
								VertexIdx,(100.0 * (double)VertexIdx)/m_UsedGraphVertices,CurComponentID,MaxVertices);
			Started = Now;
			}
		}
	if(pVertex->ComponentID > 0 && pVertex->ComponentID != (tComponentID)-1)	// if vertex already associated to a disconnected graph then try next vertex
		continue;

	NumVertices = IdentifyDiscComponent(pVertex->VertexID,CurComponentID + 1);
	if(NumVertices > 0)
		{
		// realloc for the identified components as may be required
		if((m_NumComponents + 16) >= m_AllocComponents)
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

		pComponent = &m_pComponents[m_NumComponents];
		memset(pComponent,0,sizeof(tsComponent));
		pComponent->ComponentID = CurComponentID + 1;
		pComponent->VertexID = pVertex->VertexID;
		pComponent->NumVertices = NumVertices;
		m_NumComponents += 1;
		CurComponentID += 1;
		if(NumVertices > MaxVertices)
			MaxVertices = NumVertices;
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Number of disconnected graph components: %u, max vertices in any graph: %u",CurComponentID,MaxVertices);

// sort the subgraphs or components by NumVertices and report the top 50 with 2 or more vertices in that component
if(m_NumComponents)
	{
	if(m_NumComponents > 1)
		{
		pStaticComponents = m_pComponents;
		m_MTqsort.qsort(pStaticComponents,m_NumComponents,sizeof(tsComponent),SortComponentNumVertices);
		}
	for(UINT32 Idx = 0; Idx < min(50,m_NumComponents); Idx++)
		{
		if(m_pComponents[Idx].NumVertices >= 2)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"ComponentID %u, Vertices %u",m_pComponents[Idx].ComponentID,m_pComponents[Idx].NumVertices);
		}
	// back to 
	pStaticComponents = m_pComponents;
	m_MTqsort.qsort(pStaticComponents,m_NumComponents,sizeof(tsComponent),SortComponentID);
	}
ClearEdgeTravFwdRevs();
return(m_NumComponents);
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
	m_CurTransitDepth = 0;
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


int
CAssembGraph::DumpVertex(tsGraphVertex *pVertex)
{
int Idx;
int BuffIdx;
char szBuff[2000];

BuffIdx = sprintf(szBuff,"Vertex - VertexID:%u,SeqID:%u,SeqLen:%u,ComponentID:%u",
							pVertex->VertexID,pVertex->SeqID,pVertex->SeqLen,pVertex->ComponentID);
for(Idx = 0; Idx < 2; Idx++)
	BuffIdx += sprintf(&szBuff[BuffIdx],"\n RelSense: %d flgPathScored:%u,flgPathTerm:%u, RecurseDepth:%u,PathScore:%llu,PathScoreEdgeID:%u",
							Idx, pVertex->RelSensePaths[Idx].flgPathScored,pVertex->RelSensePaths[Idx].flgPathTerm,pVertex->RelSensePaths[Idx].RecurseDepth,pVertex->RelSensePaths[Idx].PathScore,pVertex->RelSensePaths[Idx].PathScoreEdgeID);		
printf("\n%s",szBuff);
return(BuffIdx);
}

int
CAssembGraph::DumpEdge(tsGraphOutEdge *pEdge)
{
int BuffIdx;
char szBuff[2000];
BuffIdx = sprintf(szBuff,"Edge - FromVertexID:%u,ToVertexID:%u,FromSeqLen:%u,ToSeqLen:%u,FromSeq5Ofs:%u,FromSeq3Ofs:%u,ToSeq5Ofs:%u,ToSeq3Ofs:%u,Score:%u,ScoreAlignLen:%u,",
			pEdge->FromVertexID,pEdge->ToVertexID,pEdge->FromSeqLen,pEdge->ToSeqLen,pEdge->FromSeq5Ofs,pEdge->FromSeq3Ofs,pEdge->ToSeq5Ofs,pEdge->ToSeq3Ofs,pEdge->Score,pEdge->ScoreAlignLen);
BuffIdx += sprintf(&szBuff[BuffIdx]," flgOvlAntisense:%d,flgInfBackEdge:%d",
			pEdge->flgsOvlSense,pEdge->flgInfBackEdge);
printf("\n%s",szBuff);
return(BuffIdx);
}

// components contain those vertices which are linked to at least one other vertex in that component and to no other vertex in any another component
UINT32												 // number of vertices marked as members of this component
CAssembGraph::IdentifyDiscComponent(tVertID VertexID, // start component traversal from this vertex
				tComponentID ComponentID)			 // mark all traversed vertices as members of this component
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
if(pVertex->ComponentID != 0)
	return(0);

pVertex->ComponentID = ComponentID;


// possible that the initial vertex does not have any incoming or outgoing edges to any other another vertex
if(pVertex->InEdgeID == 0 && pVertex->OutEdgeID == 0)
	{
	pVertex->ComponentID = ComponentID;
	return(1);
	}

// vertex has at least one, outgoing or incoming, edge to or from some other vertex
ClearTransitStack();
PushTransitStack(VertexID);
NumMembers = 1;

while((CurVertexID = PopTransitStack()) != 0)
	{
	bTraverse = false;
	pVertex = &m_pGraphVertices[CurVertexID-1];

	if(pVertex->ComponentID != 0 && pVertex->ComponentID != ComponentID)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"IdentifyDiscComponent: Found vertex %u classified as %u when processing for component %u",CurVertexID, pVertex->ComponentID,ComponentID);
		return(-1);
		}

	if(pVertex->ComponentID == 0)
		{
		pVertex->ComponentID = ComponentID;
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
			if(!pEdge->flgTravFwd)
				{
				pEdge->flgTravFwd = 1;
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
			if(!pEdge->flgTravRev)
				{
				pEdge->flgTravRev = 1;
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

// OverlapAcceptable
// Determines if the overlap from the 'From' vertex onto the 'To' vertex would extend the 'From' vertex in the 3' direction 
INT32					// returned From sequence extension; -1 if no sequence extension
CAssembGraph::OverlapAcceptable(tsGraphOutEdge *pEdge,			// overlap edge
			  					bool bFromIsAntisense,			// 'From' vertex is this sense relative to initial starting vertex
								bool *pbToAntisense)			// accepted overlap would result in 'To' vertex having this sense relative to initial starting vertex						

{
bool bToAntisense;
INT32 ExtdFromSeqLen;

switch(pEdge->flgsOvlSense) {
	case 0:				// sense overlapping sense; accept overlaps
		if(bFromIsAntisense)
			return(-1);
		bToAntisense = false;
		break;
 
	case 2:				// sense overlapping antisense
		if(bFromIsAntisense)
			return(-1);
		bToAntisense = true;
		break;

	case 1:				// antisense overlapping sense; accept 3' overlap extension only if bIsAntisense true
		if(!bFromIsAntisense)
			return(-1);
		bToAntisense = false; // 'To' would be relative sense
		break;

	case 3:				// antisense overlapping antisense; accept 3' overlap only if bIsAntisense true
		if(!bFromIsAntisense)
			return(-1);
		bToAntisense = true; // 'To' would be relative antisense
		break;
	}

		// will overlap 3' extend the 'From' sequence?
if((pEdge->FromSeqLen - pEdge->FromSeq3Ofs) < (pEdge->ToSeqLen - pEdge->ToSeq3Ofs))
	{
	ExtdFromSeqLen = pEdge->FromSeq3Ofs + pEdge->ToSeqLen - pEdge->ToSeq3Ofs;
	*pbToAntisense = bToAntisense;
	}
else
	ExtdFromSeqLen = -1;
return(ExtdFromSeqLen);
}

static int CurStackDepth = 0;
static int MaxStackDepth = 0;

UINT64
CAssembGraph::ScorePaths(tVertID VertexID,	// score all paths starting with outgoing edges from this 'From' vertex
						bool bFromIsAntisense)	// initial starting vertex relative sense
{
bool bToIsAntisense;
bool bHighestToIsAntisense;
UINT64 EdgeScore;
UINT64 HighestScore;
UINT32 HighestScoreEdgeID;
UINT32 Idx;
int FromExtdLen;
tsGraphVertex *pVertex;
tsGraphVertex *pToNxtVertex;
tsGraphOutEdge *pEdge;
tsRelSense *pRelSense;

pVertex = &m_pGraphVertices[VertexID-1];

// vertex may have already been committed as part of a previously accepted highest scoring path
if(pVertex->flgPathAccepted)
	return(0);
pRelSense = &pVertex->RelSensePaths[bFromIsAntisense==false ? 0 : 1];

// if score already generated for paths originating at this vertex then no need to regenerate ...
if(pRelSense->flgPathScored)
	return(pRelSense->PathScore);

pVertex->RecurseDepth = 1;
pRelSense->RecurseDepth = 1;
pEdge = &m_pGraphOutEdges[pVertex->OutEdgeID-1];
for(Idx = 0; Idx < pVertex->DegreeOut; Idx++,pEdge++)
	{
	pToNxtVertex = &m_pGraphVertices[pEdge->ToVertexID-1];
	pToNxtVertex->RecurseDepth = 2;
	}

// iterate each of the outgoing edges scoring paths
HighestScore = 0;
HighestScoreEdgeID = 0;
bHighestToIsAntisense = false;
pEdge = &m_pGraphOutEdges[pVertex->OutEdgeID-1];
for(Idx = 0; Idx < pVertex->DegreeOut; Idx++,pEdge++)
	{
	if((FromExtdLen=OverlapAcceptable(pEdge,bFromIsAntisense,&bToIsAntisense)) > 0)	// explore if edge represents a 'From' extension assuming initial sense
		{
		CurStackDepth += 1;
		EdgeScore = ScorePath(pRelSense->RecurseDepth+1,pEdge,bToIsAntisense);		// score paths starting with this edge
		CurStackDepth -= 1;
		if(EdgeScore > HighestScore || (!bToIsAntisense && EdgeScore == HighestScore))	// if path score higher than any previously seen for VertexID then note the edge
			{
			bHighestToIsAntisense = bToIsAntisense;
			HighestScore = EdgeScore;
			HighestScoreEdgeID =  pVertex->OutEdgeID + Idx;
			}
		}
	}

// iterated all paths starting with vertex outgoing edges and noted the highest scoring
pRelSense->flgPathScored = 1;
pRelSense->ToIsAntisense = bHighestToIsAntisense ? 1 : 0;
if(HighestScore == 0)	// 0 if a terminating vertex with no accepted outgoing edges
	{
	pRelSense->PathScore = (UINT64)pVertex->SeqLen;
	pRelSense->flgPathTerm = 1;
	}
else
	pRelSense->PathScore = HighestScore + pEdge->ScoreAlignLen;
pRelSense->PathScoreEdgeID =  HighestScoreEdgeID;
return(pRelSense->PathScore);
}

void
CAssembGraph::ReportEdgeOverlap(const char *pszDescr,		// description used to tag the reported overlap
			 UINT32 Depth,								// current recursive depth - used to detect circular paths
			 tsGraphOutEdge *pEdge)						// edge overlap
{
const char *pszOvlSense;

switch(pEdge->flgsOvlSense) {
	case 0:
		pszOvlSense = "Sense/Sense";
		break;
	case 1:
		pszOvlSense = "Antisense/Sense";
		break;
	case 2:
		pszOvlSense = "Sense/Antisense";
		break;
	case 3:
		pszOvlSense = "Antisense/Antisense";
		break;
	}
gDiagnostics.DiagOut(eDLWarn,gszProcName,"'%s' from: %d To: %u with overlap sense '%s'",pszDescr,pEdge->FromVertexID,pEdge->ToVertexID, pszOvlSense);
}

UINT64											// highest scoring of any path from pEdge
CAssembGraph::ScorePath(UINT32 Depth,			// current recursive depth - used to detect circular paths
				  tsGraphOutEdge *pEdge,		// score paths starting with this edge
				  bool bFromIsAntisense)		// path from vertex has this path starting vertex relative sense

{
bool bToIsAntisense;
UINT32 Idx;
int FromSeqExtn;
tsGraphVertex *pToVertex;
tsGraphVertex *pToNxtVertex;
tsGraphOutEdge *pToEdge;
tsRelSense *pRelSense;
tVertID HighestScoreToVertexID;
UINT64 EdgeScore;
UINT64 HighestScore;
UINT32 HighestScoreEdgeID;
UINT32 HighScoreAlignLen;
bool bHighestToIsAntisense;
if(CurStackDepth > MaxStackDepth)
	{
	MaxStackDepth = CurStackDepth;
	if(MaxStackDepth > 1 && (MaxStackDepth % 5000) == 0)
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"ScorePath: Recursion depth %u",MaxStackDepth);
	}

// if score already generated for paths originating at this vertex then no need to regenerate ...
pToVertex = &m_pGraphVertices[pEdge->ToVertexID-1];
pRelSense = &pToVertex->RelSensePaths[bFromIsAntisense==false ? 0 : 1];
if(pRelSense->flgPathScored)
	return(pRelSense->PathScore);
pRelSense->RecurseDepth = Depth;
pRelSense->FromIsAntisense = bFromIsAntisense ? 1 : 0;

// need to mark all ToVertices with current Depth so can check circular references
pToEdge = &m_pGraphOutEdges[pToVertex->OutEdgeID-1];
for(Idx = 0; Idx < pToVertex->DegreeOut; Idx++,pToEdge++)
	{
	pToNxtVertex = &m_pGraphVertices[pToEdge->ToVertexID-1];
	if(pToNxtVertex->RecurseDepth == 0)
		pToNxtVertex->RecurseDepth = Depth;
	}

	// recursively iterate each of the outgoing edges scoring paths
HighScoreAlignLen = 0;
HighestScore = 0;
HighestScoreEdgeID = 0;
HighestScoreToVertexID = 0;
bHighestToIsAntisense = false;
bToIsAntisense = false;
pToEdge = &m_pGraphOutEdges[pToVertex->OutEdgeID-1];
for(Idx = 0; Idx < pToVertex->DegreeOut; Idx++,pToEdge++)
	{
	if(pToEdge->ToVertexID == pEdge->FromVertexID)
		continue;
	pToNxtVertex = &m_pGraphVertices[pToEdge->ToVertexID-1];
	if(pToNxtVertex->RecurseDepth > 0 &&  (pToNxtVertex->RecurseDepth + 50) < Depth)
		continue;
	if((FromSeqExtn = OverlapAcceptable(pToEdge,bFromIsAntisense,&bToIsAntisense)) > 0)
		{
		if(!pToNxtVertex->RelSensePaths[bToIsAntisense==false ? 0 : 1].flgPathScored)
			{
			if(pToNxtVertex->RelSensePaths[bToIsAntisense==false ? 0 : 1].RecurseDepth > 0)  // circular reference?
				continue;
			CurStackDepth += 1;
			EdgeScore = ScorePath(Depth+1,pToEdge,bToIsAntisense);
			CurStackDepth -= 1;
			}
		else
			EdgeScore = pToNxtVertex->RelSensePaths[bToIsAntisense==false ? 0 : 1].PathScore;
		if(EdgeScore > HighestScore || (EdgeScore == HighestScore && pToEdge->ScoreAlignLen > HighScoreAlignLen))
			{
			HighScoreAlignLen = pToEdge->ScoreAlignLen;
			HighestScore = EdgeScore;
			HighestScoreEdgeID = pToVertex->OutEdgeID + Idx;
			HighestScoreToVertexID = pToEdge->ToVertexID;
			bHighestToIsAntisense = bToIsAntisense;
			}
		}
	}

pRelSense->flgPathScored = 1;
if(HighestScore == 0)	// 0 if classing as terminating vertex with no accepted outgoing edges
	{
	pRelSense->flgPathTerm = 1;
	pRelSense->PathScore = pToVertex->SeqLen;
	}
else
	pRelSense->PathScore = HighestScore + HighScoreAlignLen;
pRelSense->ToIsAntisense =  bHighestToIsAntisense ? 1 : 0;
pRelSense->PathScoreEdgeID =  HighestScoreEdgeID;
pRelSense->RecurseDepth = 0;
return(pRelSense->PathScore);
}

int
CAssembGraph::AddTraceBackPath(bool bFirst,		// true if path starting
				 tComponentID ComponentID,      // path for this component
				 tVertID VertexID,				// path includes this vertex
				 UINT32 SeqLen,					// vertex sequence length
				 UINT32 Off5,					// sequence from this 5' offset
				 UINT32 Off3,					// to this 3' offset inclusive
				 bool bIsAntisense)            // sequence relative strand
{
tsComponent *pComponent;
tsPathTraceBack *pPathTraceback;

if((m_UsedTraceBacks + 4) >= m_AllocdTraceBacks) // 4 is simply to allow a small margin of error
	{
	size_t AllocMem;
	tsPathTraceBack *pTmp;
	size_t AllocTraceBacks;
	size_t ReallocTraceBacks;

	ReallocTraceBacks =  (size_t)((double)m_AllocdTraceBacks * cReallocTraceBacks);
	AllocTraceBacks = (size_t)m_AllocdTraceBacks + ReallocTraceBacks;
	AllocMem = (size_t)(AllocTraceBacks * sizeof(tsPathTraceBack));
#ifdef _WIN32
	pTmp = (tsPathTraceBack *)realloc(m_pPathTraceBacks,AllocMem);
#else
	pTmp = (tsPathTraceBack *)mremap(m_pPathTraceBacks,m_AllocdTraceBacks * sizeof(tsPathTraceBack),AllocMem,MREMAP_MAYMOVE);
	if(pTmp == MAP_FAILED)
		pTmp = NULL;
#endif
	if(pTmp == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"FindHighestScoringPaths: tracebacks re-allocation failed - %s",strerror(errno));
		return(eBSFerrMem);
		}
	m_AllocdTraceBacks =(UINT32) AllocTraceBacks;
	m_pPathTraceBacks = pTmp;
	}


pComponent = &m_pComponents[ComponentID-1];
if(bFirst)
	{
	if(pComponent->NumTraceBacks == 0)
		pComponent->StartTraceBackID = m_UsedTraceBacks + 1;
	pPathTraceback = &m_pPathTraceBacks[pComponent->StartTraceBackID-1];
	pComponent->NumTraceBacks = 0;
	}
m_UsedTraceBacks += 1;
pPathTraceback = &m_pPathTraceBacks[pComponent->StartTraceBackID + pComponent->NumTraceBacks - 1];
pComponent->NumTraceBacks += 1;
pPathTraceback->ComponentID = ComponentID;
pPathTraceback->Off5 = Off5;
pPathTraceback->Off3 = Off3;
pPathTraceback->IsAntisense = bIsAntisense ? 1 : 0;
pPathTraceback->VertexID = VertexID;
pPathTraceback->SeqLen = SeqLen;
return(pComponent->NumTraceBacks);
}

int														// returns path length
CAssembGraph::GenTraceBackPath(tsComponent *pComponent, // generate traceback path for this component
								bool bIsAntisense)		// initial starting vertex relative sense
{
UINT32 Vertices;
UINT32 PathLength;
tsGraphVertex *pCurVertex;
tsGraphOutEdge *pCurEdge;
bool bTermVertex;
bool bFirst;
bool bIsNxtAntisense;

pCurVertex = &m_pGraphVertices[pComponent->PathStartVertexID-1];
bTermVertex = false;
PathLength = 0;
Vertices = 0;
bFirst = true;
bIsAntisense = false;
pCurEdge = NULL;
do {
	Vertices += 1;

	if(pCurVertex->RelSensePaths[bIsAntisense==false ? 0 : 1].flgPathTerm)
		break;
	pCurEdge = &m_pGraphOutEdges[pCurVertex->RelSensePaths[bIsAntisense==false ? 0 : 1].PathScoreEdgeID-1];
	bIsNxtAntisense = pCurVertex->RelSensePaths[bIsAntisense==false ? 0 : 1].ToIsAntisense == 0 ? false : true;
	PathLength += pCurEdge->FromSeq5Ofs;

	AddTraceBackPath(bFirst,pCurVertex->ComponentID,pCurVertex->VertexID,pCurVertex->SeqLen,0,pCurEdge->FromSeq5Ofs-1,bIsAntisense);
	bFirst = false;
	bIsAntisense = bIsNxtAntisense;
	pCurVertex = &m_pGraphVertices[pCurEdge->ToVertexID-1];	
	}
while(!bTermVertex);
if(pCurEdge != NULL)
	{
	PathLength += pCurEdge->ToSeqLen;
	AddTraceBackPath(bFirst,pCurVertex->ComponentID,pCurVertex->VertexID,pCurVertex->SeqLen,0,pCurEdge->ToSeqLen-1,bIsAntisense);
	}
else
	{
	PathLength = pCurVertex->SeqLen;
	AddTraceBackPath(bFirst,pCurVertex->ComponentID,pCurVertex->VertexID,pCurVertex->SeqLen,0,pCurVertex->SeqLen-1,bIsAntisense);
	}	
return(PathLength);
}

// components contain those vertices which are linked to at least one other vertex in that component and no other vertices in another component
// identify the maximal scoring paths for each component
int												 // eBSFSuccess or otherwise
CAssembGraph::FindHighestScoringPaths(void)		 // score all possible paths and record highest scoring path for each component
{
size_t AllocMem;
tVertID Idx;
tComponentID CurCompID;
tVertID CurVertexID;
tsGraphVertex *pVertexA;
tsGraphVertex *pVertex;
tsComponent *pComponent;
UINT64 PathScore;

// graph processing
gDiagnostics.DiagOut(eDLInfo,gszProcName,"FindHighestScoringPaths: Starting ...");


// ensure m_pPathTraceBacks allocated to hold tracebacks
if(m_pPathTraceBacks == NULL)				// initialisation may be required
	{
	AllocMem = cInitialAllocTraceBacks * sizeof(tsPathTraceBack);
#ifdef _WIN32
	m_pPathTraceBacks = (tsPathTraceBack *) malloc(AllocMem);	
	if(m_pPathTraceBacks == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"FindHighestScoringPaths: memory allocation for %u tracebacks failed ",cInitialAllocTraceBacks);
		return(eBSFerrMem);
		}
#else
	m_pPathTraceBacks = (tsPathTraceBack *)mmap(NULL,AllocMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pPathTraceBacks == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"FindHighestScoringPaths: memory allocation for %u tracebacks failed ",cInitialAllocTraceBacks);
		return(eBSFerrMem);
		}
#endif
	m_AllocdTraceBacks = cInitialAllocTraceBacks;
	}
m_UsedTraceBacks = 0;

pComponent = m_pComponents;
for(CurCompID = 1; CurCompID <= m_NumComponents; CurCompID++,pComponent++)
	{
	pComponent->NumTraceBacks = 0;
	pComponent->StartTraceBackID = 0;
	pVertex = m_pGraphVertices;	
	for(CurVertexID = 1; CurVertexID <= m_UsedGraphVertices; CurVertexID++,pVertex++)
		{
		if(pVertex->ComponentID != CurCompID)
			continue;
		pVertexA = m_pGraphVertices;
		for(Idx = 0; Idx < m_UsedGraphVertices; Idx++,pVertexA++)
			{
			if(pVertexA->ComponentID != CurCompID)
				continue;
			pVertexA->RecurseDepth = 0;
			memset(pVertexA->RelSensePaths,0,sizeof(pVertexA->RelSensePaths));
			}
		PathScore = ScorePaths(CurVertexID);
		if(PathScore > 0)
			{
			if(pComponent->PathScore < PathScore)
				{
				pComponent->PathScore = PathScore;
				pComponent->PathStartVertexID = CurVertexID;
				pComponent->PathLength = GenTraceBackPath(pComponent);
				}
			}
		}
	}

UINT32 ComponentIdx;
pComponent = m_pComponents;
for(ComponentIdx = 0; ComponentIdx < min(50,m_NumComponents); ComponentIdx++,pComponent++)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Component: %d Path length: %u Path Vertices: %u Path Score: %llu", pComponent->ComponentID, pComponent->PathLength, pComponent->NumTraceBacks, pComponent->PathScore);

	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"FindHighestScoringPaths: Completed");

return(0);
}


UINT32					// number of replacements
CAssembGraph::ReplaceComponentIDs(void) 
{
UINT32 VertexIdx;
tsGraphVertex *pVertex;
tsRemapComponentID *pRemap;
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
		if(pRemap->From == pVertex->ComponentID)
			{
			pVertex->ComponentID = pRemap->To;
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
CAssembGraph::ReplaceComponentID(tComponentID ToReplaceID,		// replace existing disconnected graph identifiers 
						tComponentID ReplaceID)					// with this identifier
{
UINT32 VertexIdx;
tsGraphVertex *pVertex;
UINT32 NumReplaced = 0;
pVertex = m_pGraphVertices;
for(VertexIdx = 0; VertexIdx < m_UsedGraphVertices; VertexIdx++, pVertex++)
	{
	if(pVertex->ComponentID == ToReplaceID)
		{
		pVertex->ComponentID = ReplaceID;
		NumReplaced += 1;
		}
	}
return(NumReplaced);
}


UINT32			// index+1 in m_pGraphOutEdges of first edge with matching FromVertexID, or 0 if non matching				
CAssembGraph::LocateFirstFwdEdgeID(tVertID FromVertexID)	// find first matching  
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

tsGraphOutEdge *					// ptr to edge with matching FromVertexID and ToVertexID, or NULL if unable to locate				
CAssembGraph::LocateFromToEdge(tVertID FromVertexID,	// match edge with this FromVertexID which is
				  tVertID ToVertexID)			// to this ToVertexID
{
tsGraphOutEdge *pEl2;
int CmpRslt;
UINT32 TargPsn;
UINT32 NodeLo;				
UINT32 NodeHi;				
if(!m_bOutEdgeSorted)		// edges must be sorted by FromVertexID.ToVertexID ascending
	return(NULL);

NodeLo = 0;
NodeHi = m_UsedGraphOutEdges - 1;
do {
	TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;
	pEl2 = &m_pGraphOutEdges[TargPsn];

	if(FromVertexID == pEl2->FromVertexID && ToVertexID == pEl2->ToVertexID)
		return(pEl2);

	if(FromVertexID > pEl2->FromVertexID)
		CmpRslt = 1;
	else
		{
		if(FromVertexID < pEl2->FromVertexID)
			CmpRslt = -1;
		else
			{
			if(ToVertexID > pEl2->ToVertexID)
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

return(NULL);	// unable to locate any edge instance matching both FromVertexID and ToVertexID
}

tsGraphOutEdge *			// ptr to edge with matching FromVertexID and ToVertexID, or NULL if unable to locate				
CAssembGraph::LocateToFromEdge(tVertID ToVertexID,	// match edge with this ToVertexID which is 
				  tVertID FromVertexID)				// from this FromVertexID
{
tsGraphOutEdge *pEl2;
int CmpRslt;
UINT32 TargPsn;
UINT32 NodeLo;				
UINT32 NodeHi;				
if(!m_bInEdgeSorted)		// in edges must be sorted by ToVertexID.FromVertexID ascending
	return(NULL);

NodeLo = 0;
NodeHi = m_UsedGraphInEdges - 1;
do {
	TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;
	pEl2 = &m_pGraphOutEdges[m_pGraphInEdges[TargPsn]-1];

	if(ToVertexID == pEl2->ToVertexID && FromVertexID == pEl2->FromVertexID)
		return(pEl2);

	if(ToVertexID > pEl2->ToVertexID)
		CmpRslt = 1;
	else
		{
		if(ToVertexID < pEl2->ToVertexID)
			CmpRslt = -1;
		else   // else now matching on pEl2->ToVertexID, look for match on the pEl2->FromVertexID 
			{
			if(FromVertexID > pEl2->FromVertexID)
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

return(NULL);	
}


UINT32			// index+1 in m_pGraphInEdges of first matching ToVertexID, or 0 if non matching				
CAssembGraph::LocateFirstDnSeqID(tVertID ToVertexID)			// find first matching 
{
tsGraphOutEdge *pEl2;
int CmpRslt;
UINT32 Mark;
UINT32 TargPsn;
UINT32 NodeLo = 0;
UINT32 NodeHi = m_UsedGraphVertices-1;

if(!m_bOutEdgeSorted || !m_bInEdgeSorted)		// out edges must have been sorted by FromVertexID.ToVertexID and in edges must be sorted by ToVertexID.FromVertexID ascending
	return(0);

do {
	TargPsn = ((UINT64)NodeLo + NodeHi) / 2L;
	pEl2 = &m_pGraphOutEdges[m_pGraphInEdges[TargPsn]];

	if(ToVertexID > pEl2->ToVertexID)
		CmpRslt = 1;
	else
		if(ToVertexID < pEl2->ToVertexID)
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
			if(ToVertexID > pEl2->ToVertexID)
				CmpRslt = 1;
			else
				if(ToVertexID < pEl2->ToVertexID)
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

UINT64 
CAssembGraph::WriteContigSeqs(char *pszOutFile,  // write assembled contigs to this output file
								CSeqStore *pSeqStore)  // holds sequences used to assemble contig
{
int hOutFile;
etSeqBase *pBase;
UINT32 TraceBackIdx;
tVertID VertexID;
int ContigLen;
int SubSeqLen;
UINT8 *pContigSeq;
int ChrIdx;
int BasesLine;
char *pszFastaBuff;
tsGraphVertex *pVertex;
uint32 ComponentIdx;
tsComponent *pCurComponent;
tsPathTraceBack *pPathTraceback;
UINT32 AllocPathLength;

pCurComponent = m_pComponents;
for(AllocPathLength = ComponentIdx = 0; ComponentIdx < m_NumComponents; ComponentIdx+=1,pCurComponent+=1)
	{
	if(pCurComponent->PathLength > AllocPathLength)
		AllocPathLength = pCurComponent->PathLength;
	}
AllocPathLength = ((AllocPathLength + 100000) * 11) / 10; // note: additional allocation is to allow for fasta nl's every 80 chars plus descriptor lines

#ifdef _WIN32
hOutFile = open(pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
if((hOutFile = open(pszOutFile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	if(ftruncate(hOutFile,0)!=0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutFile,strerror(errno));
		return(eBSFerrCreateFile);
		}
#endif
if(hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutFile);
	return(eBSFerrCreateFile);
	}

if((pContigSeq = new UINT8 [AllocPathLength]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: WriteContigSeqs memory allocation failed");
	close(hOutFile);
	return(eBSFerrMem);
	}

if((pszFastaBuff = new char [AllocPathLength]) == NULL)   
	{
	delete pContigSeq;
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: WriteContigSeqs memory allocation failed");
	close(hOutFile);
	return(eBSFerrMem);
	}

pCurComponent = m_pComponents;

for(ComponentIdx = 0; ComponentIdx < m_NumComponents; ComponentIdx+=1,pCurComponent+=1)
	{
	ContigLen = 0;
	pPathTraceback = &m_pPathTraceBacks[pCurComponent->StartTraceBackID-1];
	for(TraceBackIdx = 1; TraceBackIdx <= pCurComponent->NumTraceBacks; TraceBackIdx++,pPathTraceback++)
		{
		VertexID = pPathTraceback->VertexID;
		pVertex = &m_pGraphVertices[VertexID-1];
		SubSeqLen = 1 + pPathTraceback->Off3 - pPathTraceback->Off5;
		if((SubSeqLen + ContigLen) > (int)AllocPathLength)
			{
			delete pContigSeq;
			delete pszFastaBuff;
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: Excess length (%d) contig > %d expected allowed",SubSeqLen + ContigLen, AllocPathLength);
			close(hOutFile);
			return(eBSFerrMem);
			}

		if(TraceBackIdx != pCurComponent->NumTraceBacks && pPathTraceback->IsAntisense == 1)
			{
			pSeqStore->GetSeq(pVertex->SeqID,pVertex->SeqLen - SubSeqLen,SubSeqLen,&pContigSeq[ContigLen]);
			CSeqTrans::ReverseComplement(SubSeqLen,&pContigSeq[ContigLen]);
			}
		else
			pSeqStore->GetSeq(pVertex->SeqID,pPathTraceback->Off5,SubSeqLen,&pContigSeq[ContigLen]);
		ContigLen += SubSeqLen;
		}
		
	ChrIdx = sprintf(pszFastaBuff,">Contig%d %d\n",pCurComponent->ComponentID,ContigLen);
	pBase = pContigSeq;
	int BaseIdx;

	BaseIdx = 0;
	do {
		BasesLine = min(80,ContigLen - BaseIdx);
		CSeqTrans::MapSeq2Ascii(&pContigSeq[BaseIdx],BasesLine,&pszFastaBuff[ChrIdx]);
		BaseIdx += BasesLine;
		ChrIdx += BasesLine;
		pszFastaBuff[ChrIdx] = '\n';
		ChrIdx += 1;
		}
	while(BaseIdx < ContigLen);

	CUtility::SafeWrite(hOutFile,pszFastaBuff,ChrIdx);
	}

#ifdef _WIN32
_commit(hOutFile);
#else
fsync(hOutFile);
#endif
close(hOutFile);
delete pContigSeq;
delete pszFastaBuff;
return(0);
}



// ReduceEdges
// As edges are added independently of all other edges then graph may contain many extraneous edges
// This function firstly identifies these extraneous edges and then finally removes them
// Types of extraneous edges are -
// Parallel  instances of an edges from one vertex to another vertex - may be because of retained SMRTbell hairpins
//
UINT32								// number of extraneous edges removed	 
CAssembGraph::ReduceEdges(void)		// reduce graph by detecting and removing extraneous edges
{
tsGraphOutEdge *pEdge;
tsGraphOutEdge *pPrevEdge;
tsGraphOutEdge *pDstEdge;
tsGraphVertex *pFromVertex;
tsGraphVertex *pToVertex;
tVertID ToVertexID;
tVertID FromVertexID;
int flgsOvlSense;
UINT32 NumFlgVertices;
UINT32 NumFlgEdges;
UINT32 Idx;

UINT32 Num2Remove;
UINT32 NumRemoved;


if(m_UsedGraphVertices < 2 || m_pGraphOutEdges == NULL || m_UsedGraphOutEdges < 2)	// any edges to be reduced?
	return(0);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reduce %u edges starting ...",m_UsedGraphOutEdges);
m_NumReducts = 0;

// firstly, ensure out edges are sorted by FromVertexID.ToVertexID ascending..
pStaticGraphOutEdges = m_pGraphOutEdges;
if(!m_bOutEdgeSorted || !m_bInEdgeSorted)
	return(0);

Num2Remove = 0;
NumFlgVertices = 0;
NumFlgEdges = 0;

// next identify those extraneous edges to be removed
// any read with parallel overlaps onto the same other read is assumed to be a read containing SMRTbell retained hairpins
// have no confidence in that sequence so all edges into and out of that read marked for removal
pEdge = m_pGraphOutEdges;
pPrevEdge = NULL;
FromVertexID = 0; 
ToVertexID = 0;
flgsOvlSense = 0;
for(Idx = 0; Idx < m_UsedGraphOutEdges; Idx++,pEdge++)
	{
	if(pPrevEdge != NULL && pPrevEdge->FromVertexID == pEdge->FromVertexID && pPrevEdge->ToVertexID == pEdge->ToVertexID && pPrevEdge->flgsOvlSense == pEdge->flgsOvlSense)
		{
		pFromVertex = &m_pGraphVertices[FromVertexID-1];
		pFromVertex->flgRmvEdges = 1;
		NumFlgVertices += 1;
		continue;
		}
	pPrevEdge = pEdge;
	}

if(NumFlgVertices == 0 && NumFlgEdges == 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"ReduceEdges() completed, no edges removed");
	return(0);
	}

// iterate over edges and if either the FromVertexID or ToVertexID vertices the flgRmvEdges set then mark the edge for removal
pDstEdge = m_pGraphOutEdges;
pEdge = pDstEdge;
for(Idx = 0; Idx < m_UsedGraphOutEdges; Idx++,pEdge++)
	{
	pFromVertex = &m_pGraphVertices[pEdge->FromVertexID-1];
	pToVertex = &m_pGraphVertices[pEdge->ToVertexID-1];
	if(pFromVertex->flgRmvEdges || pToVertex->flgRmvEdges)
		{
		if(pEdge->flgRemove == 0)
			{
			NumFlgEdges += 1;
			pEdge->flgRemove = 1;
			}
		}
	}

// extraneous edges have been identified and marked for removal, remove these marked edges
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reduce edges, %u edges identified for removal",NumFlgEdges);
pDstEdge = m_pGraphOutEdges;
pEdge = pDstEdge;
NumRemoved = 0;
for(Idx = 0; Idx < m_UsedGraphOutEdges; Idx++,pEdge++)
	{
	if(!pEdge->flgRemove)
		{
		if(NumRemoved)
			*pDstEdge = *pEdge;
		pDstEdge += 1;
		}
	else
		NumRemoved +=1;
	}
m_UsedGraphOutEdges -= NumRemoved;
AcquireCASSerialise();
m_NumReducts = NumRemoved;
ReleaseCASSerialise();
m_UsedGraphInEdges = m_UsedGraphOutEdges;
pStaticGraphInEdges = m_pGraphInEdges;
pStaticGraphOutEdges = m_pGraphOutEdges;
if(m_UsedGraphOutEdges >= 2)
	m_MTqsort.qsort(pStaticGraphOutEdges,m_UsedGraphOutEdges,sizeof(tsGraphOutEdge),SortOutEdgeFromVertexID);
tEdgeID *pInEdge = m_pGraphInEdges;
for(Idx = 1; Idx <= m_UsedGraphInEdges; Idx++, pInEdge++)
	*pInEdge = (tEdgeID)Idx;
if(m_UsedGraphInEdges >= 2)
	m_MTqsort.qsort(pStaticGraphInEdges,m_UsedGraphInEdges,sizeof(tEdgeID),SortInEdgesToVertexID);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reduce edges completed, removed %d edges, %u edges retained",NumRemoved,m_UsedGraphOutEdges);
return(NumRemoved);
}



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

int  // function for sorting graph vertices on ComponentID.VertexID ascending 
CAssembGraph::SortVerticesComponentID(const void *arg1, const void *arg2)
{
tsGraphVertex *pVertex1 = (tsGraphVertex *)arg1;
tsGraphVertex *pVertex2 = (tsGraphVertex *)arg2;

if(pVertex1->ComponentID < pVertex2->ComponentID)
	return(-1);
else
	if(pVertex1->ComponentID > pVertex2->ComponentID)
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
if(pEdge1->flgsOvlSense < pEdge2->flgsOvlSense)
	return(-1);
else
	if(pEdge1->flgsOvlSense > pEdge2->flgsOvlSense)
		return(1);
if(pEdge1->flgInfBackEdge < pEdge2->flgInfBackEdge)
	return(-1);
else
	if(pEdge1->flgInfBackEdge > pEdge2->flgInfBackEdge)
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

if(pEdge1->FromSeq5Ofs < pEdge2->FromSeq5Ofs)
	return(-1);
else
	if(pEdge1->FromSeq5Ofs > pEdge2->FromSeq5Ofs)
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

if(pEdge1->FromSeq5Ofs < pEdge2->FromSeq5Ofs)
	return(-1);
else
	if(pEdge1->FromSeq5Ofs > pEdge2->FromSeq5Ofs)
		return(1);
return(0);
}

int  // function for sorting incoming edges on ToVertexID.FromVertexID ascending 
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

if(pEdge1->flgsOvlSense < pEdge2->flgsOvlSense)
	return(-1);
else
	if(pEdge1->flgsOvlSense > pEdge2->flgsOvlSense)
		return(1);
if(pEdge1->flgInfBackEdge < pEdge2->flgInfBackEdge)
	return(-1);
else
	if(pEdge1->flgInfBackEdge > pEdge2->flgInfBackEdge)
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

int  // function for sorting components on NumVertices decending 
CAssembGraph::SortComponentID(const void *arg1, const void *arg2)
{
tsComponent *pComp1 = (tsComponent *)arg1;
tsComponent *pComp2 = (tsComponent *)arg2;

if(pComp1->ComponentID < pComp2->ComponentID)
	return(-1);
else
	if(pComp1->ComponentID > pComp2->ComponentID)
		return(1);
return(0);
}

