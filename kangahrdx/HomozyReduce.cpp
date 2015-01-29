// de novo read assembly into contigs
// 



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


// #define _xDEBUG 1

#include "HomozyReduce.h"

// Following DiagSequences and LogReadHit are purely for use during low level debugging
#ifndef USEWHILEDBUG1

char *
DiagCmpSequences(int SeqLen,etSeqBase *pProbeSeq,etSeqBase *pTargSeq)
{
int Idx;

static char szBuff[(cMaxFastQSeqLen + 10) * 3];
int BuffIdx;
BuffIdx=sprintf(szBuff,"P:\"%s\"\n   ",CSeqTrans::MapSeq2Ascii(pProbeSeq,SeqLen));
// compare with probe sequence and highlight the differences
for(Idx = 0; Idx < SeqLen; Idx++)
	{
	if((pTargSeq[Idx] & 0x07) != (pProbeSeq[Idx]  & 0x07))
		szBuff[BuffIdx++] = '|';
	else
		szBuff[BuffIdx++] = '.';
	}

BuffIdx+=sprintf(&szBuff[BuffIdx],"\nT:\"%s\"\n",CSeqTrans::MapSeq2Ascii(pTargSeq,SeqLen));
return(szBuff);
}

#endif

CHomozyReduce::CHomozyReduce(void)
{
m_pSfxdSeqs = NULL;
m_pConcatSeqs = NULL;
m_pSuffixArray = NULL;
m_pSeqs2Assemb = NULL;
m_pContigSeq = NULL;
m_pszLineBuff = NULL;
m_pPathIDs = NULL;
m_pContigLengths = NULL;

m_hInFile = -1;
m_hOutFile = -1;

m_bMutexesCreated = false;
m_NumThreads = 0;
memset(&m_RdsSfxHdr,0,sizeof(m_RdsSfxHdr));
m_StopWatch.Start();
Reset(false);

// on Windows need to get the baseline memory working set as will almost certainly be increasing these limits later!  
#ifdef _WIN32
SYSTEM_INFO SystemInfo;
HANDLE hProcess;
hProcess = GetCurrentProcess();
GetSystemInfo(&SystemInfo);
m_WinPageSize = SystemInfo.dwPageSize;
BOOL bRslt = GetProcessWorkingSetSize(hProcess,(PSIZE_T)&m_BaseWinMinMem,(PSIZE_T)&m_BaseWinMaxMem);
if(bRslt == false)
	{
	m_BaseWinMinMem = 0;
	m_BaseWinMaxMem = 0;
	}

#else
//m_BaseWinMinMem = 0;
//m_BaseWinCurMaxMem = 0;	
#endif
}

CHomozyReduce::~CHomozyReduce(void)
{
Reset(false);
}

int
CHomozyReduce::Reset(bool bSync)
{

if(m_hOutFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutFile);
#else
		fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}


if(m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
	}


if(m_pSuffixArray != NULL)
	{
#ifdef _WIN32
	free(m_pSuffixArray);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSuffixArray != MAP_FAILED)
		munmap(m_pSuffixArray,m_AllocMemSfx);
#endif	
	m_pSuffixArray = NULL;
	}

if(m_pConcatSeqs != NULL)
	{
#ifdef _WIN32
	free(m_pConcatSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pConcatSeqs != MAP_FAILED)
		munmap(m_pConcatSeqs,m_AllocMemConcat);
#endif	
	m_pConcatSeqs = NULL;
	}

if(m_pSeqs2Assemb != NULL)
	{
#ifdef _WIN32
	free(m_pSeqs2Assemb);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqs2Assemb != MAP_FAILED)
		munmap(m_pSeqs2Assemb,m_AllocMemSeqs2Assemb);
#endif	
	m_pSeqs2Assemb = NULL;
	}

if(m_pPathIDs != NULL)
	{
	free(m_pPathIDs);
	m_pPathIDs = NULL;
	}

if(m_pContigLengths != NULL)
	{
#ifdef _WIN32
	free(m_pContigLengths);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pContigLengths != MAP_FAILED)
		munmap(m_pContigLengths,m_AllocdCtgLens);
#endif	
	m_pContigLengths = NULL;
	}

if(m_pSfxdSeqs != NULL)
	{
#ifdef _WIN32
	free(m_pSfxdSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSfxdSeqs != MAP_FAILED)
		munmap(m_pSfxdSeqs,m_AllocdMemSfxdSeqs);
#endif	
	m_pSfxdSeqs = NULL;
	}

if(m_pContigSeq != NULL)
	{
#ifdef _WIN32
	free(m_pContigSeq);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pContigSeq != MAP_FAILED)
		munmap(m_pContigSeq,m_AllocContigSeq);
#endif	
	m_pContigSeq = NULL;
	}

if(m_pszLineBuff != NULL)
	{
	delete m_pszLineBuff;
	m_pszLineBuff = NULL;
	}

if(m_bMutexesCreated)
	{
	DeleteMutexes();
	m_bMutexesCreated = false;
	}

m_AllocMemConcat = 0;
m_AllocMemSfx = 0;
m_NumSfxdSeqs = 0;
m_AllocdNumSfxdSeqs = 0;
m_AllocdMemSfxdSeqs = 0;
m_MeanReadLen = 0;
m_CurMaxMemWorkSetBytes = 0;
m_BuffNumOverlaidContigs = 0;

m_TotSeqsParsed = 0;
m_TotSeqsUnderLen = 0;
m_TotSeqsExcessNs = 0;
m_TotSeqs2Assemb = 0;
m_NumSeqs2Assemb = 0;
m_Seqs2AssembLen = 0;
m_AllocMemSeqs2Assemb = 0; 
m_TotSeqs2AssembBases = 0;

m_InvalidRefs = 0;

m_NumSuffixEls = 0;
m_AllocContigSeq = 0;
m_AllocSeqIDsinContig = 0;
m_AllocSeqIDsinContigMem = 0;
m_NumContigsinContig = 0;
m_szInFile[0] = '\0';
m_szOutFile[0] = '\0';

m_LineBuffLen = 0;

m_TotLenCovered = 0;

m_AllocdPathIDs = 0;
m_NumPathIDs = 0;

m_NumRawFiles = 0;
m_NumRdsFiles = 0;
m_AllocdCtgLens = 0; 
m_NumGenContigs = 0;
m_TotNumReduced = 0;
m_TotNumNonReduced = 0;
m_MinCtgLen = 0;

memset(&m_RdsSfxHdr,0,sizeof(m_RdsSfxHdr));
return(0);
}

void
CHomozyReduce::SetNumThreads(int maxThreads)
{
m_NumThreads = maxThreads;
m_MTqsort.SetMaxThreads(maxThreads);
}

int
CHomozyReduce::CreateMutexes(void)
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
if((m_hMtxIterContigs = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxIterContigs,NULL)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

#ifdef _WIN32
if((m_hMtxMHContigs = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxMHContigs,NULL)!=0)
	{
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
#ifdef _WIN32
	CloseHandle(m_hMtxIterContigs);
#else
	pthread_rwlock_destroy(&m_hRwLock);
	pthread_mutex_destroy(&m_hMtxIterContigs);
#endif
	return(eBSFerrInternal);
	}

m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CHomozyReduce::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
#ifdef _WIN32
CloseHandle(m_hMtxIterContigs);
CloseHandle(m_hMtxMHContigs);

#else
pthread_mutex_destroy(&m_hMtxIterContigs);
pthread_mutex_destroy(&m_hMtxMHContigs);
pthread_rwlock_destroy(&m_hRwLock);

#endif
m_bMutexesCreated = false;
}


void
CHomozyReduce::ValidateSeq(char *pszDescr,int Len, etSeqBase *pSeq)
{
etSeqBase *pBase;
int Idx;
if(pszDescr == NULL || Len == 0 || pSeq == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ValidateSeq '%s' invalid parameters: Len: %d pSeq: %s",pszDescr,Len,pSeq == NULL ? "NULL" : "non-null");
	return;
	}

pBase = pSeq;
for(Idx = 0; Idx < Len; Idx++, pBase++)
	{
	if(*pBase > eBaseN)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ValidateSeq '%s' have a illegal base in sequence (length %d) problem at offset %d",pszDescr,Len,Idx);

		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Sequence'%s'", CSeqTrans::MapSeq2Ascii(pSeq,Len));
		return;
		}
	}
}

#ifdef _WIN32
unsigned __stdcall CHomozyReduce::ThreadedContigsAssemb(void * pThreadPars)
#else
void * CHomozyReduce::ThreadedContigsAssemb(void * pThreadPars)
#endif
{
int Rslt;

CHomozyReduce *pThis;
int Idx;
tsAssembThreadPars *pPars = (tsAssembThreadPars *)pThreadPars; // makes it easier not having to deal with casts!
tsProbesBlock ProbesBlock;
UINT32 ProbeSeqID;
tsSfxdSeq *pSfxdSeq;

UINT16 ProbeLen;
UINT8 *pProbeSeq;
int NumNonReduced;
int NumReduced;

NumNonReduced = 0;
NumReduced = 0;

pThis = (CHomozyReduce *)pPars->pThis;
memset(&ProbesBlock,0,sizeof(ProbesBlock));
ProbesBlock.MaxContigs = cMaxProbesPerBlock;

while(pThis->ThreadedIterProbes(&ProbesBlock))
	{
	for(Idx = 0; Idx < ProbesBlock.NumContigs; Idx++)
		{
		pPars->NumContigsProc += 1;

		ProbeSeqID = ProbesBlock.ReadSeqIDs[Idx];
		pSfxdSeq = &pThis->m_pSfxdSeqs[ProbeSeqID-1];
				
		pProbeSeq = &pThis->m_pConcatSeqs[pSfxdSeq->ConcatSeqOfs+5];	// returns ptr to first base of concatenated sequence
		ProbeLen = pSfxdSeq->ReadLen;					

		UINT32 TargSeqID = 0;
		UINT32 MergeLen = 0;
		etSeqBase *pMergeSeq = NULL;

		Rslt =  pThis->MarkHomozygoticRegions(ProbeSeqID,	// identifies probe sequence
						 pProbeSeq,					// probe
						 ProbeLen,					// probe length
						 pPars);					// calling thread parameters


		switch(Rslt) {
			case eLOTnone:							// no homozygotic regions in this sequence
				NumNonReduced += 1;
				break;

			case eLOThit:							// able to overlay onto another sequence (or could be that probe has been corrected for polymorphisms)
				NumReduced += 1;
				break;

			default:
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ThreadedContigsAssemb: LocateOverlaidTarg internal error - returned %d",Rslt);
				pPars->Rslt = Rslt;
#ifdef _WIN32
				_endthreadex(0);
				return(eBSFSuccess);
#else
				pthread_exit(NULL);
#endif
			}
		}
	}

pThis->AcquireSerialise();

pThis->m_TotNumReduced += NumReduced;	
pThis->m_TotNumNonReduced += NumNonReduced;

pThis->ReleaseSerialise();

pPars->Rslt = 1;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

void
CHomozyReduce::AcquireSerialise(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxIterContigs,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterContigs);
#endif
}

void
CHomozyReduce::ReleaseSerialise(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxIterContigs);
#else
pthread_mutex_unlock(&m_hMtxIterContigs);
#endif
}

void 
inline CHomozyReduce::AcquireLock(bool bExclusive)
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
inline CHomozyReduce::ReleaseLock(bool bExclusive)
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
CHomozyReduce::ResetThreadedIterContigs(void) // must be called by master thread prior to worker threads calling ThreadedIterContigs()
{
m_NumSeqsProc = 0;
m_NxtReadProcSeqID = 0;
m_ProcessingStartSecs = m_StopWatch.ReadUSecs();
}


// Returns the number of reads thus far aligned 
UINT32
CHomozyReduce::ApproxNumContigsAligned(void)
{
UINT32 NumAligned;
AcquireLock(false);
NumAligned = m_NumSeqsProc;
ReleaseLock(false);
return(NumAligned);
}

// ThreadedIterContigs
// Iterates over all reads and returns a block of reads for processing to the calling thread
// The body of this function is serialised  
bool	// returns false if no more reads available for processing by calling thread
CHomozyReduce::ThreadedIterProbes(tsProbesBlock *pRetBlock)	// iterate and return blocks of read probes to be processed by each thread
{
UINT32 NumContigsLeft;
UINT32 MaxContigs2Proc;
pRetBlock->NumContigs = 0;

AcquireSerialise();

if(m_pSfxdSeqs == NULL || 
	m_AllocdNumSfxdSeqs == 0 || (UINT32)m_NumSeqsProc >= m_NumSfxdSeqs) // if all sequences processed then time to move onto next processing phase
	{
	pRetBlock->NumContigs = 0;
	ReleaseSerialise();
	return(false);
	}

// adjust pRetBlock->MaxContigs according to the number of reads remaining and threads still processing these reads
// idea is to maximise the number of threads still processing when most reads have been processed so that
// the last thread processing doesn't end up with a large block of reads needing lengthly processing
NumContigsLeft = m_NumSfxdSeqs - m_NumSeqsProc;
if(NumContigsLeft < cMaxProbesPerBlock/8)	// if < cMaxProbesPerBlock/8 yet to be processed then give it all to the one thread
	MaxContigs2Proc = NumContigsLeft;
else
	MaxContigs2Proc = min((UINT32)pRetBlock->MaxContigs,10 + (NumContigsLeft / (UINT32)m_NumThreads));
MaxContigs2Proc = min(MaxContigs2Proc,NumContigsLeft);
if(!m_NumSeqsProc)			   // 0 if first
	m_NxtReadProcSeqID = 1;

while(MaxContigs2Proc && (UINT32)m_NumSeqsProc < m_NumSfxdSeqs)
	{
	// check if this read sequence is still available to be aligned
	if(!(m_pSfxdSeqs[m_NxtReadProcSeqID-1].Flags & (cFlagNA | cFlagMerged)))	// no lock required as only single byte thus always consistent
		{
		pRetBlock->ReadSeqIDs[pRetBlock->NumContigs++] = m_NxtReadProcSeqID;
		MaxContigs2Proc -= 1;
		}
	m_NxtReadProcSeqID += 1;
	m_NumSeqsProc += 1;
	}

if((UINT32)m_NumSeqsProc > m_NumSfxdSeqs)
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ThreadedIterProbes: unexpected iteration, please report to stuart.stephen@csiro.au");

ReleaseSerialise();
return(true);
}

bool
CHomozyReduce::SetMaxMemWorkSetSize(size_t Bytes)
{
// need to determine available physical memory
#ifdef _WIN32
static int m_WorkingSetSizeRejected = 0;	// count of number of times SetProcessWorkingSetSize() fails
BOOL bRslt;
HANDLE hProcess;
size_t ReqMaxSize;
size_t ReqMinSize;

if(Bytes < 10000000)		// ensure there will be a reasonable workingset size
	Bytes = 10000000;

Bytes += (m_WinPageSize-1);		// round up to a windows memory page size
Bytes /= m_WinPageSize;
Bytes *= m_WinPageSize;

if(m_WinPageSize == 0 || (m_CurMaxMemWorkSetBytes > 0 && 
   ((double)Bytes  >  (m_CurMaxMemWorkSetBytes * 0.9) && (double)Bytes < m_CurMaxMemWorkSetBytes * 1.1)))
	return(true);

hProcess = GetCurrentProcess();

ReqMaxSize = Bytes;
ReqMinSize = (ReqMaxSize+1)/2;

bRslt = SetProcessWorkingSetSize(hProcess,ReqMinSize,ReqMaxSize);	// can only but try...

if(bRslt == false && (m_WorkingSetSizeRejected++ > 5))
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"SetMaxMemWorkSetSize: unable to SetProcessWorkingSetSize for min %lld max %lld bytes (page size is %d)",
					ReqMinSize,ReqMaxSize,m_WinPageSize);
	}
if(bRslt)
	m_CurMaxMemWorkSetBytes = (UINT64)ReqMaxSize;
return(bRslt > 0 || m_WorkingSetSizeRejected < 5 ? true : false);
#else
return(true);	
#endif
}


int
CHomozyReduce::ReduceHomozygosity(etPMode PMode,			// processing mode
					char *pszRsltFile,				// write reduced homozygosity contigs into this file
					bool bStrand,					// strand specific homozygous region reduction - homozygous regions between any two contigs must be in same orientation
					int MaxHomozySubs,				// characterise as homozygotic if substitution rate between regions <= this rate per 100bp
					int MinHomozyLen,				// homozygotic regions to be at least this length
					int MinHetrozyLen,				// island (marked as homozygotic either flank) hetrozygotic regions must be at least this length otherwise treat as homozygotic
					int MinCtgLen)					// filter out homozygotic region reduced contigs of less than this length
{
int Rslt;
int ThreadIdx;
tsAssembThreadPars WorkerThreads[cMaxWorkerThreads];

UINT8 *pProbe;

int PrevNumContigsAligned;
int CurNumContigsAligned;

UINT32 NumAcceptedAsMerged;

int CurMaxIter;

if(m_TotSeqs2Assemb == 0)	// have to be assembling at least one read!
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No sequences to assemble");
	return(-1);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %u sequences totalling %lld bases for homozygotic regions...",m_TotSeqs2Assemb,m_TotSeqs2AssembBases);


// processing sensitivity determines CurMaxIter
// In general the higher the sensitivity then more iterations are allowed
// Higher sensitivities increase resource requirements and require much longer run times
switch(PMode) {
	case ePMUltraSens:			// ultra sensitive - much slower
		CurMaxIter = 25000;
		break;

	case ePMMoreSens:			// more sensitive - slower
		CurMaxIter = 15000;
		break;

	case ePMLessSens:			// less sensitive - quicker
		CurMaxIter = 2000;
		break;

	case ePMdefault:			// default processing mode
		CurMaxIter = 10000;
		break;
	}

m_MinCtgLen = MinCtgLen;
m_MaxHomozySubs = MaxHomozySubs;
m_MinHomozyLen = MinHomozyLen;

strcpy(m_szOutFile,pszRsltFile);

// where to output assembled contigs
#ifdef _WIN32
m_hOutFile = open(pszRsltFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(pszRsltFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszRsltFile,strerror(errno));
			Reset(false);
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszRsltFile);
	Reset(false);
	return(eBSFerrCreateFile);
	}

// allocate for output buffering
if(m_pszLineBuff == NULL)
	{
	if((m_pszLineBuff = new char [cAllocLineBuffLen]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: output buffering allocation of %lld bytes - %s",(INT64)cAllocLineBuffLen,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}
	m_AllocLineBuff = cAllocLineBuffLen;
	}


m_NumGenContigs = 0;

// allocate for contig length distributions
if(m_pContigLengths == NULL)
	{
	size_t MemReq = cAllocCtgLenDist * sizeof(int);
#ifdef _WIN32
	m_pContigLengths = (int *) malloc(MemReq);	// initial and perhaps the only allocation
	if(m_pContigLengths == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation of %d bytes for contig length distributions - %s",MemReq,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pContigLengths = (int *)mmap(NULL,MemReq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pContigLengths == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation of %d bytes through mmap()  failed - %s",MemReq,strerror(errno));
		m_pContigLengths = NULL;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	m_AllocdCtgLens = MemReq;
	m_NumGenContigs = 0;
	m_TotLenCovered = 0;
	}

pProbe = NULL;
int ProbeSeqID = 0;
int ExpectedID = 1;
int NoOverlayCnt = 0;
int HasSingleOverlay = 0;
int HasMultiOverlays = 0;
int TotNodeInstances = 0;
int OverDist[cMaxOverlayDist+2];
memset(OverDist,0,sizeof(OverDist));

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting contig assembly processing...");

CreateMutexes();
ResetThreadedIterContigs();

NumAcceptedAsMerged = 0;
m_NumSeqs2Assemb = m_TotSeqs2Assemb;

if((Rslt = GenSfxdSeqs()) < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to generate suffixed sequences");
	Reset(false);
	return(Rslt);
	}

	// now generate suffix array over sequences
if((Rslt=GenRdsSfx()) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Generating suffix array failed");
	return(Rslt);
	}


// startup the processing threads
memset(WorkerThreads,0,sizeof(WorkerThreads));

int TmpNumThreads = m_NumThreads;
#ifdef _xDEBUG
printf("\nWARNING: Currently Restricting number of threads....\n");
m_NumThreads = 1;
#endif
for(ThreadIdx = 0; ThreadIdx < m_NumThreads; ThreadIdx++)
	{
	WorkerThreads[ThreadIdx].ThreadIdx = ThreadIdx + 1;
	WorkerThreads[ThreadIdx].pThis = this;
	WorkerThreads[ThreadIdx].PMode = PMode;
	WorkerThreads[ThreadIdx].ElSize = m_SfxElSize;
	WorkerThreads[ThreadIdx].CoreLen = cAbsMinCoreLen;
	WorkerThreads[ThreadIdx].MinCtgLen = m_MinCtgLen;
	WorkerThreads[ThreadIdx].MaxHomozySubs = m_MaxHomozySubs;
	WorkerThreads[ThreadIdx].MinHomozyLen = m_MinHomozyLen;

	if((WorkerThreads[ThreadIdx].pProbeSeq = new UINT8 [cAllocMergerSeqLen+10])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for probe sequences...");
		exit(1);
		}
	WorkerThreads[ThreadIdx].AllocProbeSeq = cAllocMergerSeqLen;

	WorkerThreads[ThreadIdx].bStrand = bStrand;
	WorkerThreads[ThreadIdx].CurMaxIter = CurMaxIter;

#ifdef _WIN32
	WorkerThreads[ThreadIdx].threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,ThreadedContigsAssemb,&WorkerThreads[ThreadIdx],0,&WorkerThreads[ThreadIdx].threadID);
#else
	WorkerThreads[ThreadIdx].threadRslt =	pthread_create (&WorkerThreads[ThreadIdx].threadID , NULL , ThreadedContigsAssemb , &WorkerThreads[ThreadIdx] );
#endif
	}

	// let user know that we are working hard...
PrevNumContigsAligned = ApproxNumContigsAligned();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: %u contigs processed",PrevNumContigsAligned);

	// wait for all threads to have completed
for(ThreadIdx = 0; ThreadIdx < m_NumThreads; ThreadIdx++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject( WorkerThreads[ThreadIdx].threadHandle, 60000))
		{
		CurNumContigsAligned = ApproxNumContigsAligned();
		if(CurNumContigsAligned > PrevNumContigsAligned)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: %u contigs processed",CurNumContigsAligned);
			PrevNumContigsAligned = CurNumContigsAligned;
		}
	CloseHandle( WorkerThreads[ThreadIdx].threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while((JoinRlt = pthread_timedjoin_np(WorkerThreads[ThreadIdx].threadID, NULL, &ts)) != 0)
		{
		CurNumContigsAligned = ApproxNumContigsAligned();
		if(CurNumContigsAligned > PrevNumContigsAligned)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: %u contigs processed",CurNumContigsAligned);
		PrevNumContigsAligned = CurNumContigsAligned;
		ts.tv_sec += 60;
		}
#endif
	if(WorkerThreads[ThreadIdx].pProbeSeq != NULL)
		{
		delete WorkerThreads[ThreadIdx].pProbeSeq;
		WorkerThreads[ThreadIdx].pProbeSeq = NULL;
		}
	WorkerThreads[ThreadIdx].AllocProbeSeq = 0;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: completed processing %u contigs",m_NumSeqsProc);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: %u (%1.2f%%) identified as containing homozygotic regions",m_TotNumReduced,(100.0 * m_TotNumReduced)/m_NumSeqsProc);
	
if(m_pSuffixArray != NULL)
	{
#ifdef _WIN32
	free(m_pSuffixArray);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSuffixArray != MAP_FAILED)
		munmap(m_pSuffixArray,m_AllocMemSfx);
#endif	
	m_pSuffixArray = NULL;
	}

UINT32 SumNumContigsProc = 0;				// returned number of reads processed

tsAssembThreadPars *pThread;
pThread = WorkerThreads;

// reduction ended
// iterate over sequences and output to file as contigs those sequence regions not marked as being homozygotic
// if a region is less than 100bp then extend flanks of that region until at least 100bp
UINT32 Idx;
UINT32 SeqIdx;
UINT32 SeqLen;
UINT32 MarkedLen;
UINT32 UnmarkedLen;
UINT32 UnmarkedStartIdx;
UINT32 FlankLen;
tsSfxdSeq *pSfxdSeq;
UINT8 *pSeq;
UINT8 *pUnmarkedSeq;
for(Idx = 0; Idx < m_NumSfxdSeqs; Idx++)
	{
	pSfxdSeq = &m_pSfxdSeqs[Idx];
	pSeq = &m_pConcatSeqs[pSfxdSeq->ConcatSeqOfs+5];
	SeqLen = pSfxdSeq->ReadLen;
	if(SeqLen < (UINT32)MinCtgLen)
		continue;

	MarkedLen = 0;
	UnmarkedLen = 0;
	UnmarkedStartIdx = 0;
	pUnmarkedSeq = NULL;
	for(SeqIdx = 0; SeqIdx < pSfxdSeq->ReadLen; SeqIdx++,pSeq++)
		{
		if((*pSeq & 0xf0) == cMarkMskFlg)		// has been marked?
			{
			if(UnmarkedLen >= (UINT32)MinHetrozyLen)		// has to be worth the effort!
				{
				if(UnmarkedLen < (UINT32)MinCtgLen)
					{
					FlankLen = ((UINT32)MinCtgLen - UnmarkedLen)/2;
					if(UnmarkedStartIdx >= FlankLen)
						UnmarkedStartIdx -= FlankLen;
					else
						UnmarkedStartIdx = 0;
					if(UnmarkedStartIdx + MinCtgLen > pSfxdSeq->ReadLen)
						UnmarkedStartIdx = pSfxdSeq->ReadLen - MinCtgLen;
					pUnmarkedSeq = &m_pConcatSeqs[pSfxdSeq->ConcatSeqOfs+ 5 + UnmarkedStartIdx];
					UnmarkedLen = MinCtgLen;
					}
				WriteContigSeq(pUnmarkedSeq,UnmarkedLen);
				}
			MarkedLen += 1;
			UnmarkedLen = 0;
			pUnmarkedSeq = NULL;
			}
		else							// else unmarked
			{
			if(pUnmarkedSeq == NULL)
				{
				pUnmarkedSeq = pSeq;
				UnmarkedStartIdx = SeqIdx;
				MarkedLen = 0;
				}
			UnmarkedLen += 1;
			}
		}
	if(UnmarkedLen >= (UINT32)MinHetrozyLen)
		{
		if(UnmarkedLen < (UINT32)MinCtgLen)
			{
			FlankLen = ((UINT32)MinCtgLen - UnmarkedLen)/2;
			if(UnmarkedStartIdx >= FlankLen)
				UnmarkedStartIdx -= FlankLen;
			else
				UnmarkedStartIdx = 0;
			if(UnmarkedStartIdx + MinCtgLen > pSfxdSeq->ReadLen)
				UnmarkedStartIdx = pSfxdSeq->ReadLen - MinCtgLen;
			pUnmarkedSeq = &m_pConcatSeqs[pSfxdSeq->ConcatSeqOfs+ 5 + UnmarkedStartIdx];
			UnmarkedLen = MinCtgLen;
			}
		WriteContigSeq(pUnmarkedSeq,UnmarkedLen);
		}
	}

if(m_LineBuffLen > 0)
	{
	if(write(m_hOutFile,m_pszLineBuff,m_LineBuffLen)!=m_LineBuffLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write to disk on file %s - error %s",m_szOutFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	}


#ifdef _WIN32
_commit(m_hOutFile);
#else
fsync(m_hOutFile);
#endif
close(m_hOutFile);
m_hOutFile = -1;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d homozygotic region reduced contigs were generated covering %lld bases",m_NumGenContigs,m_TotLenCovered);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Contig Nx lengths:");

INT64 NxLens[11];	// to hold contig lengths for all Nx where Nx varies from 10 to 100 in increments of N10
INT64 SumCtgLens;
INT64 NxSum;
INT64 *pNxLens;
int CtgLenIdx;
int CtgNxIdx;
int *pCnts;
int NumCtgs;

SumCtgLens = 0;
pCnts = m_pContigLengths;
for(CtgLenIdx = 0; CtgLenIdx < m_NumGenContigs; CtgLenIdx++,pCnts++)
	SumCtgLens += *pCnts;

m_MTqsort.qsort(m_pContigLengths,(UINT64)m_NumGenContigs,sizeof(int),SortByContigLen);

pNxLens = NxLens;
for(CtgNxIdx = 1; CtgNxIdx <= 10; CtgNxIdx++,pNxLens++)
	*pNxLens = (SumCtgLens * CtgNxIdx)/10;

NxSum = 0;
CtgNxIdx = 0;
NumCtgs = 0;
pCnts = m_pContigLengths;
for(int Idx = 0; Idx < m_NumGenContigs; Idx++,pCnts++)
	{
	if(*pCnts)
		{
		NxSum += *pCnts;
		NumCtgs += 1;
		}

	if(NxSum >= NxLens[CtgNxIdx])
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"  N%d contig length is %d, Contigs: %d (%1.2f%%)",(CtgNxIdx+1)*10,*pCnts,NumCtgs, ((double)NumCtgs * 100)/(double)m_NumGenContigs);
		NxLens[CtgNxIdx++] = *pCnts;
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed");
Reset(true);
return(eBSFSuccess);
}




// IterConcatSeqs
// use to iterate over concatenated sequences returning ptr to next sequence following the current read
// To start from first read then pass in NULL as pCurSeq
// pCurSeq must pointing to the XFormID for the current sequence, this sequence is iterated until the next cCSeqSep and ptr to the XFormID is returned
// Returns NULL if all read hits have been iterated (cCSeqEOS encountered starting next sequence)
UINT8 *
CHomozyReduce::IterConcatSeqs(UINT8 *pCurSeq)
{
if(pCurSeq == NULL)
	return(&m_pConcatSeqs[1]); // skip over initial cCSeqBOS, returns ptr to the XFormID for the first sequence
pCurSeq += 5;				   // skip over XFormID
while(*pCurSeq < cCSeqSep)     // skip to end of current sequence
	pCurSeq += 1;
if(*pCurSeq == cCSeqEOS)	  // check if last has been already iterated
	return(NULL);
return(pCurSeq+1);			  // returns ptr to the XFormID for the returned sequence
}

// GetpSeqmarker
// Returns ptr to marker terminating current sequence ptd at by pSeq
// *pSeq must be pointing to the sequence and not to the separator or XFormID
// Returns ptr to separator or NULL if errors
UINT8 *
CHomozyReduce::GetpConcatSeqMarker(UINT8 *pSeq)
{
if(pSeq == NULL || *pSeq >= cCSeqSep)
	return(NULL);
while(*pSeq++ < cCSeqSep);
return(&pSeq[-1]);
}

// GetConcatSeqID
// Get sequence identifier for concatenated sequence ptd to by pSeq
// *pSeq must be pointing to the sequence or XFormID
// Returns sequence identifier or -1 if errors
INT32
CHomozyReduce::GetConcatSeqID(UINT8 *pSeq)
{
if(pSeq == NULL)
	return(-1);
while(*pSeq < cCSeqSep)			// backup pSeq until pointing into the XForm'd sequence separator immediately preceding the sequence
	pSeq -= 1;					
while(pSeq[1] >= cCSeqSep)      // forward until last byte of XForm'd sequence separator  
	pSeq++;
return(XFormToID(*(UINT64 *)(pSeq-4) & 0x0fffffffff));
}

// GetConcatSeqLen
// returns length of concatenated sequence ptd at by pSeq
// *pSeq must be pointing to the XFormID or into the sequence
// Returns sequence length or -1 if errors
int    
CHomozyReduce::GetConcatSeqLen(UINT8 *pSeq)
{
int SeqLen;
if(pSeq == NULL)
	return(-1);

while(pSeq[-1] < cCSeqSep)
	pSeq -= 1;
while(*pSeq >= cCSeqSep)  // if not pointing into sequence then skip into start of sequence
	pSeq++;
SeqLen = 0;							// start counting..
while(*pSeq++ < cCSeqSep)
	SeqLen += 1;
return(SeqLen);
}

etSeqBase *
CHomozyReduce::GetConcatSeqStart(UINT8 *pSeq)			// returns ptr to first base of concatenated sequence
{
if(pSeq == NULL)
	return(NULL);
while(pSeq[-1] < cCSeqSep)
	pSeq -= 1;
while(*pSeq >= cCSeqSep)  // if not pointing into sequence then skip into start of sequence
	pSeq++;
return(pSeq);
}

// GetConcatSeqOfs
// returns offset (0..n) at which pSeq is ptd at within concatenated sequence
// *pSeq must be pointing to the sequence and not to the separator or XFormID
// Returns offset or -1 if errors
int    
CHomozyReduce::GetConcatSeqOfs(UINT8 *pSeq)
{
int SeqOfs;
if(pSeq == NULL || *pSeq >= cCSeqSep)
	return(-1);
SeqOfs = 0;
while(pSeq[-1] < cCSeqSep)			// backup to start of sequence
	{
	SeqOfs++;						// whilst counting..
	pSeq -= 1;
	}
return(SeqOfs);
}

// CmpProbeTarg
// Compares probe against target taking into account any concatenated sequence marker or XFormID etc
int
CHomozyReduce::CmpProbeTarg(etSeqBase *pEl1,etSeqBase *pEl2,int Len)
{
int Psn;
for(Psn=0; Psn < Len; Psn++,pEl1++,pEl2++)
	{
	if(*pEl2 >= cCSeqSep)
		return(-1);
	if(*pEl1 > *pEl2)
		return(1);
	if(*pEl1 < *pEl2)
		return(-1);
	}
return(0);
}

UINT64			// index+1 in pSfxArray of first exactly matching probe or 0 if no match				
CHomozyReduce::LocateFirstExact(int ElSize,		// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
				  etSeqBase *pProbe,			// pts to probe sequence
				  int ProbeLen,					// probe length to exactly match over
				  etSeqBase *pTarg,				// target sequence
				  UINT8 *pSfxArray,				// target sequence suffix array
				  UINT64 TargStart,				// position in pTarg (0..n) corresponding to start of suffix array
				  UINT64 SfxLo,					// low index in pSfxArray
				  UINT64 SfxHi)					// high index in pSfxArray
{
etSeqBase *pEl1;
etSeqBase *pEl2;
UINT8 El1;
UINT8 El2;
UINT8 Byte;

int CmpRslt;
int Ofs;
UINT64 Mark;
UINT64 TargPsn;
UINT64 TargEl;

do {
	pEl1 = pProbe;
	TargPsn = ((UINT64)SfxLo + SfxHi) / 2L;
	TargEl = (TargPsn + TargStart) * ElSize;
	if(ElSize == 4)
		pEl2 = &pTarg[*(UINT32 *)&pSfxArray[TargEl]];
	else
		pEl2 = &pTarg[Unpack5(&pSfxArray[TargEl])];

	CmpRslt = 0;
	for(Ofs=0; Ofs < ProbeLen; Ofs++)
		{
		El2 = (Byte = *pEl2++) & 0x0f;
		if(Byte >= cCSeqSep || El2 == eBaseEOS)
			{
			CmpRslt = -1;
			break;
			}
		El1 = (Byte = *pEl1++) & 0x0f;
		if(El1 > El2)
			{
			CmpRslt = 1;
			break;
			}
		if(El1 < El2)
			{
			CmpRslt = -1;
			break;
			}
		}

	if(!CmpRslt)	// if a match then may not be the lowest indexed match
		{
		if(TargPsn == 0 || SfxLo == TargPsn) // check if already lowest
			return(TargPsn + 1);
		// iterate until lowest located
		while(1) {
			if(CmpRslt == 0)
				{
				Mark = TargPsn;
				if(Mark == 0)
					return(Mark+1);
				SfxHi = TargPsn - 1;
				}	
			TargPsn = ((UINT64)SfxLo + SfxHi) / 2L;

			TargEl = (TargPsn + TargStart) * ElSize;
			if(ElSize == 4)
				pEl2 = &pTarg[*(UINT32 *)&pSfxArray[TargEl]];
			else
				pEl2 = &pTarg[Unpack5(&pSfxArray[TargEl])];

		
			pEl1 = pProbe;
			CmpRslt = 0;
			for(Ofs=0; Ofs < ProbeLen; Ofs++)
				{
				El2 = (Byte = *pEl2++) & 0x0f;
				if(Byte >= cCSeqSep || El2 == eBaseEOS)
					{
					CmpRslt = -1;
					break;
					}
				El1 = (Byte = *pEl1++) & 0x0f;
				if((Byte & 0x80) > El2 || El1 > El2)
					{
					CmpRslt = 1;
					break;
					}
				if(El1 < El2)
					{
					CmpRslt = -1;
					break;
					}
				}
			if(CmpRslt == 0)				// 0 if still matching
				continue;
			SfxLo = TargPsn + 1;
			if(SfxLo == Mark)
				return(Mark+1);
			}
		}

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		SfxHi = TargPsn - 1;
		}
	else
		SfxLo = TargPsn+1;
	}
while(SfxHi >= SfxLo);

return(0);	// unable to locate any instance of pProbe
}


UINT64			// index+1 in pSfxArray of last exactly matching probe or 0 if no match					
CHomozyReduce::LocateLastExact(int ElSize,		// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
				  etSeqBase *pProbe,			// pts to probe sequence
				  int ProbeLen,					// probe length to exactly match over
				  etSeqBase *pTarg,				// target sequence
				  UINT8 *pSfxArray,				// target sequence suffix array
				  UINT64 TargStart,				// position in pTarg (0..n) corresponding to start of suffix array
				  UINT64 SfxLo,					// low index in pSfxArray
				  UINT64 SfxHi,					// high index in pSfxArray
				  UINT32 Limit)					// if non-zero then need only iterate towards last exactly matching this for this Limit iterations
{
etSeqBase *pEl1;
etSeqBase *pEl2;
UINT8 El1;
UINT8 El2;
UINT8 Byte;

int CmpRslt;
int Ofs;

UINT64 TargEl;
UINT64 Mark;
UINT64 TargPsn;
UINT64 LoTargPsn = TargStart;
UINT64 SfxHiMax = SfxHi;
UINT64 SfxLoMax = SfxLo;

if(Limit > 0)					// Limit is a soft limit so add a few more on 
	Limit += 10;

do {
	pEl1 = pProbe;
	TargPsn = (SfxLo + SfxHi) / 2L;
	TargEl = (TargPsn + TargStart) * ElSize;

	if(ElSize == 4)
		pEl2 = &pTarg[*(UINT32 *)&pSfxArray[TargEl]];
	else
		pEl2 = &pTarg[Unpack5(&pSfxArray[TargEl])];

	CmpRslt = 0;
	for(Ofs=0; Ofs < ProbeLen; Ofs++)
		{
		El2 = (Byte = *pEl2++) & 0x0f;
		if(Byte >= cCSeqSep || El2 == eBaseEOS)
			{
			CmpRslt = -1;
			break;
			}
		El1 = (Byte = *pEl1++) & 0x0f;
		if((Byte & 0x80) > El2 || El1 > El2)
			{
			CmpRslt = 1;
			break;
			}
		if(El1 < El2)
			{
			CmpRslt = -1;
			break;
			}
		}

	if(!CmpRslt)	// if a match then may not be the highest indexed match
		{
		LoTargPsn = TargPsn;
		if(TargPsn == SfxHiMax || SfxHi == TargPsn) // check if already highest
			return(TargPsn + 1);
		// iterate until highest located
		while(1) {
			if(CmpRslt == 0)
				{
				Mark = TargPsn;
				if(Mark == SfxHi || (Limit > 0 && ((TargPsn - LoTargPsn) > Limit)))
					return(Mark+1);
				SfxLo = TargPsn + 1;
				}	
			TargPsn = ((UINT64)SfxLo + SfxHi) / 2L;

			TargEl = (TargPsn + TargStart) * ElSize;
			if(ElSize == 4)
				pEl2 = &pTarg[*(UINT32 *)&pSfxArray[TargEl]];
			else
				pEl2 = &pTarg[Unpack5(&pSfxArray[TargEl])];

			pEl1 = pProbe;
			CmpRslt = 0;
			for(Ofs=0; Ofs < ProbeLen; Ofs++)
				{
				El2 = (Byte = *pEl2++) & 0x0f;
				if(Byte >= cCSeqSep || El2 == eBaseEOS)
					{
					CmpRslt = -1;
					break;
					}
				El1 = (Byte = *pEl1++) & 0x0f;
				if((Byte & 0x080) > El2 || El1 > El2)
					{
					CmpRslt = 1;
					break;
					}
				if(El1 < El2)
					{
					CmpRslt = -1;
					break;
					}
				}
				
			if(CmpRslt == 0)				// 0 if still matching
				continue;
			SfxHi = TargPsn - 1;
			if(SfxHi == Mark)
				return(Mark+1);
			}
		}

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		SfxHi = TargPsn - 1;
		}
	else
		SfxLo = TargPsn+1;
	}
while(SfxHi >= SfxLo);

return(0);	// unable to locate any instance of pProbe
}


// MarkHomozygoticRegions
// Locate and mark homozygotic regions shared between contigs except for one instance of each region
// These marked regions will be subsequently removed resulting in contigs which are more consensus representative
tLOTRslt						// < 0 if errors, eLOTnone or eLOTprobe or eLOTtarg if no matches, eLOThit if overlap onto target accepted 
CHomozyReduce::MarkHomozygoticRegions(UINT32 ProbeSeqID,// identifies probe sequence
						 etSeqBase *pProbeSeq,			// probe sequence 
						 int ProbeLen,					// probe length 
						 tsAssembThreadPars *pPars)	    // calling thread parameters
{
int ElSize;
int CurCoreSegOfs;				// current core segment relative start
int IterCnt;					// count iterator for current segment target matches

int CurOverlayLen;				// current overlay length

int NumMarkedRegions;
int CurMMCnt;					// current number of mismatches for current target sequence being processed

etSeqBase *pProbeBase;
etSeqBase *pTargBase;
UINT64 TargIdx;

int CoreDelta;					// core window offset increment (1..n)
int CoreLen;					// core window length 

int Ofs;
UINT8 Base1;
UINT8 Base2;
etSeqBase *pEl1;
etSeqBase *pEl2;

UINT32 NumTargSeqProc;

UINT32 TargSeqID;

bool bFirstIter;				// set false after the first subsequence core returned by LocateFirstExact has been processed

UINT64 LastTargIdx;
UINT32 NumCopies;

int TargMatchLen;

UINT64 SfxElOfs;
UINT64 SfxElVal;

UINT8 ProbeBase;
UINT8 TargBase;

int CurCoreDelta;

etSeqBase *pTarg;			// target sequence
UINT8 *pSfxArray;			// target sequence suffix array
UINT64 SfxLen;				// number of suffixs in pSfxArray

tsSfxdSeq *pTargSfxdSeq;
tsSfxdSeq *pProbeSfxdSeq;

bool bSelfHit;
int AllowedMismatches;

int Passes;
bool bRevCpl;

UINT8 SubWin[100];			// used to hold aligner induced substutions over last 100 aligned bases
int SubWinIdx;		    	// index into SubWin of next psn to write probe to target match
int SubsInWin;				// current number of substitutions within SubWin

etSeqBase *pTargLeftStart;

int ProbeMatchStart;
int ProbeMatchEnd;
int ProbeMatchSubs;

// ensure suffix array block loaded for iteration!
if(m_pSuffixArray == NULL || m_NumSuffixEls == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"MarkHomozygoticRegions: Invalid input - m_pSuffixArray == NULL or m_NumSuffixEls == 0");
	return((tLOTRslt)eBSFerrInternal);
	}

pTarg = (etSeqBase *)&m_pConcatSeqs[1];
pSfxArray = (UINT8 *)m_pSuffixArray;
SfxLen = (UINT64)m_NumSuffixEls;

CoreLen = min(cAbsMaxCoreLen,max(cAbsMinCoreLen, pPars->CoreLen));
CoreDelta = CoreLen;

NumTargSeqProc = 0;

pProbeSfxdSeq = &m_pSfxdSeqs[ProbeSeqID-1];

if(pPars->bStrand)			
	Passes = 1;				
else
	Passes = 2;				

if(pPars->AllocProbeSeq < (UINT32)ProbeLen + 100)	// safety margin!
	{
	delete pPars->pProbeSeq;
	pPars->AllocProbeSeq = ProbeLen + 10000; 
	if((pPars->pProbeSeq = new UINT8 [pPars->AllocProbeSeq]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"MarkHomozygoticRegions: Memory allocation for pPars->pProbeSeq of %d bytes through new() failed, ProbeLen = %d",pPars->AllocProbeSeq,ProbeLen);
		return((tLOTRslt)eBSFerrMem);
		}
	}
memcpy(pPars->pProbeSeq,pProbeSeq,ProbeLen);
pPars->pProbeSeq[ProbeLen] = 0x80;
NumMarkedRegions = 0;

while(Passes--)
	{
	if(pProbeSfxdSeq->Flags & (cFlagNA | cFlagMerged))
		return(eLOTprobe);

	if(Passes == 0 && !pPars->bStrand)				// if 1 then strand independent so need to revcpl the probe sequence
		{
		CSeqTrans::ReverseComplement(ProbeLen,pPars->pProbeSeq);
		bRevCpl = true;
		}
	else
		bRevCpl = false;

	bSelfHit = false;

	ElSize = pPars->ElSize;
	CurCoreDelta = CoreDelta;

	for(CurCoreSegOfs = 0; 
		CurCoreSegOfs <= (ProbeLen - CoreLen) && 
		CurCoreDelta > 0;	
		CurCoreSegOfs += CurCoreDelta)
		{
		if((CurCoreSegOfs + CoreLen + CurCoreDelta) > ProbeLen)
			CurCoreDelta = ProbeLen - (CurCoreSegOfs + CoreLen);

		// try and find at least one exact match of the core against a subsequence of same length in any sequence - will later check for self hits 
		TargIdx = LocateFirstExact(ElSize,&pPars->pProbeSeq[CurCoreSegOfs],CoreLen,pTarg,pSfxArray,0,0,SfxLen-1);

		if(!TargIdx)        // 0 if no core segment matches
			continue;
		
		// have at least one core match
		TargIdx -= 1;
		IterCnt = 0;
		NumCopies = 0;
		bFirstIter = true;		// set false after the first subsequence core returned by LocateFirstExact has been processed
		while(!pPars->CurMaxIter || IterCnt < pPars->CurMaxIter)	// only check a limited number of putative core extensions	
			{
			if(!bFirstIter)	   // if not first core match (subsequent matches for same core are only putative)
				{				
				// ensure not about to iterate past end of suffix array!
				if((TargIdx + 1) >= SfxLen)
					break;

				SfxElOfs = (TargIdx+1) * ElSize;
				if(ElSize == 4)
					SfxElVal = (UINT64)*(UINT32 *)&pSfxArray[SfxElOfs];
				else
					SfxElVal = Unpack5(&pSfxArray[SfxElOfs]);


				if((SfxElVal + CoreLen) >  SfxLen) 
					break;

				if(IterCnt == cChkIterDepth && !NumCopies)
					{
					// check how many more exact copies there are of the current probe subsequence, if too many then don't bother exploring these
					LastTargIdx = LocateLastExact(ElSize,&pPars->pProbeSeq[CurCoreSegOfs],CoreLen,pTarg,pSfxArray,0,TargIdx-1,SfxLen-1,(UINT32)pPars->CurMaxIter+1);
					NumCopies = (UINT32)(LastTargIdx > 0 ? 1 + LastTargIdx - TargIdx : 0);
					if(pPars->CurMaxIter && NumCopies > (UINT32)pPars->CurMaxIter)		// only checking at the cChkIterDepth iteration is not too cpu resource intensive
						break;
					}

				// check that this new putative core is still matching
				pTargBase = &pTarg[SfxElVal]; 
				pProbeBase = &pPars->pProbeSeq[CurCoreSegOfs];

				pEl1= pProbeBase;
				pEl2 = pTargBase;
				for(Ofs=0; Ofs < CoreLen; Ofs++,pEl1++,pEl2++)
					{
					Base1 = *pEl1;
					Base2 = *pEl2;
					if(Base1 >= 0x80)
						Base1 = eBaseEOS;
					else
						Base1 &= 0x07;
					if(Base2 >= 0x80)
						Base2 = eBaseEOS;
					else
						Base2 &= 0x07;
					if((Base1 != Base2) || Base1 == eBaseEOS)
						break;
					}
				if(Ofs != CoreLen)			// will be not equal if target no longer matches
					break;					// try next core segment

				// confirmed that core still exactly matches target core - no longer putative!
				TargIdx += 1;
				}

			bFirstIter = false;
			SfxElOfs = TargIdx * ElSize;
			if(ElSize == 4)
				SfxElVal = (UINT64)*(UINT32 *)&pSfxArray[SfxElOfs];
			else
				SfxElVal = Unpack5(&pSfxArray[SfxElOfs]);

			if(SfxElVal < 1)			// shouldn't occur but who knows...
				{
				AcquireLock(true);
				m_InvalidRefs += 1;
				ReleaseLock(true);
				continue;
				}

			// not interested in self-hits so check if target identifier same as probe identifier
			// if same then is a self-hit...
			pTargBase = &pTarg[SfxElVal]; 
			TargSeqID = GetConcatSeqID(pTargBase);	
			if(TargSeqID < 1 || TargSeqID > m_NumSfxdSeqs)			// shouldn't occur but who knows...
				{
				AcquireLock(true);
				m_InvalidRefs += 1;
				ReleaseLock(true);
				continue;
				}

			if(TargSeqID == ProbeSeqID)
				{
				bSelfHit = true;	
				continue;
				}

			// confirmed as not being a self hit
			// check target and slough target hits which are shorter than probe or if same length have identifier less than probe
			pTargSfxdSeq = &m_pSfxdSeqs[TargSeqID - 1];
			if(pTargSfxdSeq->Flags & cFlagNA || (pTargSfxdSeq->ReadLen < (UINT32)ProbeLen || (pTargSfxdSeq->ReadLen == (UINT32)ProbeLen && TargSeqID < ProbeSeqID)))
				continue;

			TargMatchLen = ProbeLen; 

			// ensure comparisons are still within start/end range of target sequence/assembly
			if((SfxElVal + CoreLen + 1) > SfxLen)		// added 1 purely for safety until fully debugged!
				continue;
 
			NumTargSeqProc += 1;
			IterCnt += 1;

			// now do the matching allowing for mismatches
			// ProbeSeq pts to first base in probe
			// CurCoreSegOfs contains the current core relative offset from the probe start
			// pSfxArray[TargIdx] pts to base in TargSeq corresponding to probe core start
			
			// determine if probe is overlapping onto target, and set pProbeBase to probe sequence overlapping target start
			// checks if 5' probe maps on or before target 5' sequence
			// does not allow for for 5' probe to start after 5' sequence
			// what is needed is to allow probe to be extended left and right, counting number of aligner induced substitutions, until
			// either probe or target boundaries encountered. Then can determine if min overlap length has been exceeded and number
			// of substitutions is acceptable

			int TargLeftOfs;				// relative offset in target at which probe starts
			int ProbeLeftOfs;				// relative offset in probe at which target starts

			pTargBase = &pTarg[SfxElVal];	// pTargBase set to 5' base of matching target core
			pTargLeftStart = pTargBase;
			pProbeBase = &pPars->pProbeSeq[CurCoreSegOfs];

			// CurCoreSegOfs is the probe core segment start
			ProbeLeftOfs = CurCoreSegOfs;
			ProbeMatchStart = CurCoreSegOfs;
			ProbeMatchSubs = 0;
			TargLeftOfs = 0;

			CurMMCnt = 0;
			CurOverlayLen = CoreLen;
			pTargBase -= 1;
			pProbeBase -= 1;
			
			memset(SubWin,0,sizeof(SubWin));
			SubWinIdx = CoreLen;
			SubsInWin = 0;
			AllowedMismatches = pPars->MaxHomozySubs;
			
			do
				{
				TargBase = *pTargBase--;
				if(TargBase >= 0x80)
					break;
				else
					TargBase &= 0x07;
				
				if(ProbeLeftOfs > 0)
					{
					if(SubWinIdx >= 100 && SubWin[(SubWinIdx - 100) % 100] == 1)
						SubsInWin -= 1;
					ProbeBase = 0x07 & *pProbeBase--;
					if(TargBase != ProbeBase)
						{
						CurMMCnt += 1;
						SubsInWin += 1;
						SubWin[SubWinIdx % 100] = 1;
						}
					else
						{
						SubWin[SubWinIdx % 100] = 0;
						ProbeMatchStart = ProbeLeftOfs;
						ProbeMatchSubs = SubsInWin;
						}
					SubWinIdx += 1;
					ProbeLeftOfs -= 1;
					pTargLeftStart -= 1;
					}
				else
					TargLeftOfs += 1;
				}
			while(SubsInWin <= AllowedMismatches);

			// now extend right
			SubsInWin = ProbeMatchSubs;
			pTargBase = &pTarg[SfxElVal];					// pTargBase set to 5' base of matching target core
			pProbeBase = &pPars->pProbeSeq[CurCoreSegOfs];
			ProbeMatchEnd = CurCoreSegOfs + CoreLen;
			pTargBase += CoreLen;							
			pProbeBase += CoreLen;
			do
				{
				TargBase = *pTargBase++;
				if(TargBase >= 0x80)
					break;
				else
					TargBase &= 0x07;
				
				ProbeBase = *pProbeBase++;
				if(ProbeBase >= 0x80)
					break;
				if(SubWinIdx >= 100 && SubWin[(SubWinIdx - 100) % 100] == 1)
					SubsInWin -= 1;
				ProbeBase &= 0x07;
				if(TargBase != ProbeBase)
					{
					CurMMCnt += 1;
					SubsInWin += 1;
					SubWin[SubWinIdx % 100] = 1;
					}
				else
					{
					SubWin[SubWinIdx % 100] = 0;
					ProbeMatchEnd += 1;
					}
				SubWinIdx += 1;
				}
			while(SubsInWin <= AllowedMismatches);

			CurOverlayLen = ProbeMatchEnd - ProbeMatchStart;
			if(CurOverlayLen < pPars->MinHomozyLen) 
				continue;

			// if overlay is less than 100 then need to prorate the assembler induced substitutions to be proportionally the same as 100bp overlay 
			if(CurOverlayLen < 100 && SubsInWin > 0)
				SubsInWin = (100*SubsInWin) / CurOverlayLen;
			if(SubsInWin > AllowedMismatches)
				continue;

			// mark the probe sequence with this identified homozygotic region
			// as only can be one writer, and write/reads at byte level then no need for serialisation
			if(!bRevCpl)
				{
				pTargBase = &pProbeSeq[ProbeMatchStart];
				while(ProbeMatchStart++ < ProbeMatchEnd)
					*pTargBase++ |= cMarkMskFlg;
				}
			else
				{
				pTargBase = &pProbeSeq[ProbeLen- 1 - ProbeMatchStart];
				while(ProbeMatchStart++ < ProbeMatchEnd)
					*pTargBase-- |= cMarkMskFlg;
				}
			NumMarkedRegions += 1;
			}
		}
	}

return(NumMarkedRegions > 0 ? eLOThit : eLOTnone);				    // at least one region
}




void
CHomozyReduce::DumpContigSeqs(char *pszDescr)
{
// iterate and output..
UINT32 Idx;
tsSfxdSeq *pContigSeq;
etSeqBase *pSeq;
UINT32 ReadLen;
for(Idx = 0; Idx < m_NumSfxdSeqs; Idx++)
	{
	pContigSeq = &m_pSfxdSeqs[Idx];
	ReadLen = pContigSeq->ReadLen;
	pSeq = &m_pConcatSeqs[pContigSeq->ConcatSeqOfs + 5]; 
	DumpContigSeq(Idx+1,pszDescr,pSeq,ReadLen);
	}
return;
}

void
CHomozyReduce::DumpContigSeq(int IdentID,char *pszDescr,etSeqBase *pSeq,int SeqLen)
{
int CurSubSeqLen;
char *pContigSeqTxt;
printf("\n>%s%d %s\n",m_szCtgDescr,IdentID,pszDescr);
while(SeqLen)
	{
	CurSubSeqLen = min(79,SeqLen);
	pContigSeqTxt = CSeqTrans::MapSeq2Ascii(pSeq,CurSubSeqLen);
	printf("%s\n",pContigSeqTxt);
	SeqLen -= CurSubSeqLen;
	pSeq += CurSubSeqLen;
	}
}

void
CHomozyReduce::SetCtgDescr(char *pszCtgDescr)
{
strncpy(m_szCtgDescr,pszCtgDescr,sizeof(m_szCtgDescr)-1);
m_szCtgDescr[sizeof(m_szCtgDescr)-1] = '\0';
}

int 
CHomozyReduce::WriteContigSeq(etSeqBase *pSeq,int SeqLen)
{
int CurSubSeqLen;

if(pSeq ==  NULL || SeqLen == 0)
	return(eBSFerrParams);

AcquireSerialise();

// is more memory required to hold this new contig length?
if(((m_NumGenContigs + 10) * sizeof(int)) > m_AllocdCtgLens)	// allow a little safety margin
	{
	int *pTmp;
	size_t AllocNeeded = m_AllocdCtgLens + (cAllocCtgLenDist * sizeof(int));
#ifdef _WIN32
	pTmp = (int *) realloc(m_pContigLengths,AllocNeeded);
#else
	pTmp = (int *)mremap(m_pContigLengths,m_AllocdCtgLens,AllocNeeded,MREMAP_MAYMOVE);
	if(pTmp == MAP_FAILED)
		pTmp = NULL;
#endif
	if(pTmp == NULL)
		{
		ReleaseSerialise();
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenContig memory reallocation for contig lengths distribution of %d bytes failed..",AllocNeeded);
		return(eBSFerrMem);
		}
	m_pContigLengths = pTmp;
	m_AllocdCtgLens = AllocNeeded;
	}
m_pContigLengths[m_NumGenContigs++] = SeqLen;
m_TotLenCovered += (size_t)SeqLen;

if((m_LineBuffLen + 500) > m_AllocLineBuff)
	{
	if(write(m_hOutFile,m_pszLineBuff,m_LineBuffLen)!=m_LineBuffLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write to disk on file %s - error %s",m_szOutFile,strerror(errno));
		ReleaseSerialise();
		return(eBSFerrFileAccess);
		}
	m_LineBuffLen = 0;
	}

m_LineBuffLen += sprintf(&m_pszLineBuff[m_LineBuffLen],">%s%d Len:%d\n",m_szCtgDescr,m_NumGenContigs,SeqLen);
while(SeqLen)
	{
	CurSubSeqLen = min(79,SeqLen);
	m_LineBuffLen += sprintf(&m_pszLineBuff[m_LineBuffLen],"%s\n",CSeqTrans::MapSeq2Ascii(pSeq,CurSubSeqLen));
	
	if((m_LineBuffLen + 500) > m_AllocLineBuff)
		{
		if(write(m_hOutFile,m_pszLineBuff,m_LineBuffLen)!=m_LineBuffLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write to disk on file %s - error %s",m_szOutFile,strerror(errno));
			ReleaseSerialise();
			return(eBSFerrFileAccess);
			}
		m_LineBuffLen = 0;
		}

	SeqLen -= CurSubSeqLen;
	pSeq += CurSubSeqLen;
	}

ReleaseSerialise();
return(eBSFSuccess);
}


// SortByContigLen
// Sorts contig lengths descending

int
CHomozyReduce::SortByContigLen(const void *arg1, const void *arg2)
{
int *pEl1 = (int *)arg1;
int *pEl2 = (int *)arg2;

if(*pEl1 < *pEl2)
	return(1);
if(*pEl1 > *pEl2)
	return(-1);
return(0);
}


