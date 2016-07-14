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

#include "./biokanga.h"
#include "./Kangadna.h"

#ifdef _WIN32
#pragma intrinsic(_InterlockedCompareExchange)
#endif

static UINT8 *m_xpConcatSeqs;	// to hold all concatenated packed sequences
static int SeqWrdBytes;			// sequence packing used in m_xpConcatSeqs

CKangadna::CKangadna(void)
{
m_pContaminants = NULL;
m_pPartialSeqs2Assemb = NULL;
m_pAcceptLevDist = NULL; 
m_pBlockNsLoci = NULL;
memset(&m_Sequences,0,sizeof(m_Sequences));
m_pszLineBuff = NULL;
m_hInFile = -1;
m_hOutFile = -1;
m_hOrphansFile = -1;
m_hInSeqTypesFile = -1;
m_hOutSeqTypesFile = -1;
m_hOutFastaSE = -1;
m_hOutFastaR1 = -1;
m_hOutFastaR2 = -1;
m_bMutexesCreated = false;
m_NumThreads = 0;
m_PMode = 0;
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

CKangadna::~CKangadna(void)
{
Reset(false);
}

int
CKangadna::Reset(bool bSync)
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

if(m_hOutSeqTypesFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutSeqTypesFile);
#else
		fsync(m_hOutSeqTypesFile);
#endif
	close(m_hOutSeqTypesFile);
	m_hOutSeqTypesFile = -1;
	}


if(m_hOutFastaSE != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutFastaSE);
#else
		fsync(m_hOutFastaSE);
#endif
	close(m_hOutFastaSE);
	m_hOutFastaSE = -1;
	}
if(m_hOutFastaR1 != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutFastaR1);
#else
		fsync(m_hOutFastaR1);
#endif
	close(m_hOutFastaR1);
	m_hOutFastaR1 = -1;
	}
if(m_hOutFastaR2 != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutFastaR2);
#else
		fsync(m_hOutFastaR2);
#endif
	close(m_hOutFastaR2);
	m_hOutFastaR2 = -1;
	}


if(m_hOrphansFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOrphansFile);
#else
		fsync(m_hOrphansFile);
#endif
	close(m_hOrphansFile);
	m_hOrphansFile = -1;
	}

if(m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
	}

if(m_hInSeqTypesFile != -1)
	{
	close(m_hInSeqTypesFile);
	m_hInSeqTypesFile = -1;
	}

if(m_pContaminants != NULL)
	{
	delete m_pContaminants;
	m_pContaminants = NULL;
	}

ResetTypeSeqs();

if(m_pBlockNsLoci != NULL)
	{
#ifdef _WIN32
	free(m_pBlockNsLoci);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pBlockNsLoci != MAP_FAILED)
		munmap(m_pBlockNsLoci,m_AllocdBlockNsLociSize);
#endif	
	m_pBlockNsLoci = NULL;
	}

if(m_pPartialSeqs2Assemb != NULL)
	{
#ifdef _WIN32
	free(m_pPartialSeqs2Assemb);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pPartialSeqs2Assemb != MAP_FAILED)
		munmap(m_pPartialSeqs2Assemb,m_AllocdPartialSeqs2Assemb);
#endif	
	m_pPartialSeqs2Assemb = NULL;
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

if(m_pAcceptLevDist != NULL)
	{
	delete m_pAcceptLevDist;
	m_pAcceptLevDist = NULL;
	}
m_LevDistKMerLen = 0;
m_MaxLev = 0;
m_PenaltyMatch = 0;
m_PenaltyMisMatch = cLevenshteinDefault;
m_PenaltyInsert = cLevenshteinDefault;
m_PenaltyDelete = cLevenshteinDefault;
m_AcceptLevWidth = 0;
m_CurMaxMemWorkSetBytes = 0;
m_StartProcSeqID = 1;
m_PMode = 0;
m_pszLineBuff = NULL;		
m_AllocLineBuff = 0;
memset(&m_SeqEsts,0,sizeof(m_SeqEsts));
m_szInFile[0] = '\0';
m_szOutFile[0] = '\0';
m_szOrphansFile[0] = '\0';
m_szProcSeqTypesFile[0] = '\0';
m_szOutSeqTypesFile[0] = '\0';
m_LineBuffLen = 0;
m_NumRawFiles = 0;
m_NumRdsFiles = 0;
m_Sequences.SfxSparsity = eSSparsity15;
m_Sequences.SeqWrdBytes = 4;
SeqWrdBytes = 4;
m_Sequences.SeqHdrLen = 3;
m_bAffinity = false;
m_bDedupeIndependent = false;
m_bRawDedupe = false;

m_AllocdPartialSeqs2Assemb = 0;
m_PartialSeqs2AssembOfs = 0;
m_NumPartialSeqs2Assemb = 0;
m_LenPartialSeqs2Assemb = 0;

m_NumBlockNsLoci = 0;
m_AllocdBlockNsLoci = 0;
m_AllocdBlockNsLociSize = 0;

m_MinOverlapReductFact = 0;	
m_InitialReqPEPrimOverlap = 0;
m_InitialReqPESecOverlap = 0;
m_InitialReqPESumOverlap = 0;
m_InitialReqSEPrimOverlap = 0;
m_MinReqPEPrimOverlap = 0;
m_MinReqPESecOverlap = 0;
m_MinReqPESumOverlap = 0;
m_MinReqSEPrimOverlap = 0;
m_MinPE2SEOvlp = 0;
m_MinReqPESepDist = 0;
m_MaxReqPESepDist = 0;

memset(&m_SrcFiles,0,sizeof(m_SrcFiles));
return(0);
}


void
CKangadna::ResetTypeSeqs(void)
{
if(m_Sequences.pSuffixArray != NULL)
	{
#ifdef _WIN32
	free(m_Sequences.pSuffixArray);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_Sequences.pSuffixArray != MAP_FAILED)
		munmap(m_Sequences.pSuffixArray,m_Sequences.AllocMemSfx);
#endif	
	m_Sequences.pSuffixArray = NULL;
	}

if(m_Sequences.pSeqs2Assemb != NULL)
	{
#ifdef _WIN32
	free(m_Sequences.pSeqs2Assemb);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_Sequences.pSeqs2Assemb != MAP_FAILED)
		munmap(m_Sequences.pSeqs2Assemb,m_Sequences.AllocMemSeqs2Assemb);
#endif	
	m_Sequences.pSeqs2Assemb = NULL;
	}

if(m_Sequences.pSeqStarts != NULL)
	{
#ifdef _WIN32
	free(m_Sequences.pSeqStarts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_Sequences.pSeqStarts != MAP_FAILED)
		munmap(m_Sequences.pSeqStarts,m_Sequences.AllocMemSeqStarts);
#endif	
	m_Sequences.pSeqStarts = NULL;
	}


if(m_Sequences.pSeqFlags != NULL)
	{
#ifdef _WIN32
	free(m_Sequences.pSeqFlags);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_Sequences.pSeqFlags != MAP_FAILED)
		munmap(m_Sequences.pSeqFlags,m_Sequences.AllocMemSeqFlags);
#endif	
	m_Sequences.pSeqFlags = NULL;
	}

if(m_Sequences.pTmpRevCplSeqs != NULL)
	delete (UINT8 *)m_Sequences.pTmpRevCplSeqs;

memset(&m_Sequences,0,sizeof(tsSequences));
}

teBSFrsltCodes
CKangadna::SetNumThreads(int maxThreads,bool bAffinity)
{
if(maxThreads < 0 || maxThreads > cMaxWorkerThreads)
		return(eBSFerrParams);
m_NumThreads = maxThreads;
m_bAffinity = bAffinity;
m_MTqsort.SetMaxThreads(maxThreads);
CreateMutexes();
return(eBSFSuccess);
}

void
CKangadna::SetDedupePE(bool bDedupeIndependent)		// dedupe policy on paired ends
{
m_bDedupeIndependent = bDedupeIndependent;
}

void 
CKangadna::SetPMode(int PMode)								// set processing mode
{
m_PMode = PMode;
m_Sequences.SeqWrdBytes = 4;
}

int
CKangadna::FreeSfx(void)
{
if(m_Sequences.pSuffixArray != NULL)
	{
#ifdef _WIN32
	free(m_Sequences.pSuffixArray);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_Sequences.pSuffixArray != MAP_FAILED)
		munmap(m_Sequences.pSuffixArray,m_Sequences.AllocMemSfx);
#endif	
	m_Sequences.pSuffixArray = NULL;
	}
m_Sequences.NumSuffixEls = 0;			// number of elements in suffix array
m_Sequences.AllocMemSfx = 0;				// allocated memory size for suffix array
return(eBSFSuccess);
}

int
CKangadna::FreeSeqStarts(bool bFreeFlags)	// optionally also free flags array
{
if(m_Sequences.pSeqStarts != NULL)
	{
#ifdef _WIN32
	free(m_Sequences.pSeqStarts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_Sequences.pSeqStarts != MAP_FAILED)
		munmap(m_Sequences.pSeqStarts,m_Sequences.AllocMemSeqStarts);
#endif	
	m_Sequences.pSeqStarts = NULL;
	}
m_Sequences.NumSeqStarts = 0;			// number of sequence descriptors used in m_pSeqStarts;
m_Sequences.AllocMemSeqStarts = 0;		// memory (bytes) currently allocated for m_pSeqStarts 
if(bFreeFlags)
	{
	if(m_Sequences.pSeqFlags != NULL)
		{
#ifdef _WIN32
		free(m_Sequences.pSeqFlags);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(m_Sequences.pSeqFlags != MAP_FAILED)
			munmap(m_Sequences.pSeqFlags,m_Sequences.AllocMemSeqStarts);
#endif	
		m_Sequences.pSeqFlags = NULL;
		}
	m_Sequences.NumSeqFlags = 0;
	m_Sequences.AllocMemSeqFlags = 0;
	}
return(eBSFSuccess);
}


int 
CKangadna::GetSeqWrdBytes(void)			// returns size of tSeqWd 
{
return(m_Sequences.SeqWrdBytes);
}

// set suffix sparsity
teBSFrsltCodes 
CKangadna::SetSfxSparsity(etSfxSparsity SfxSparsity) 
{
if(!(SfxSparsity == eSSparsity1 || SfxSparsity == eSSparsity15))
	return(eBSFerrParams);

m_Sequences.SfxSparsity = SfxSparsity;
m_Sequences.SeqWrdBytes = SfxSparsity == eSSparsity1 ? 1 : 4;
m_Sequences.SeqHdrLen = 3;
return(eBSFSuccess);
}


int
CKangadna::CreateMutexes(void)
{
if(m_bMutexesCreated)
	return(eBSFSuccess);

m_CASSeqFlags = 0;
m_CASNxtProcRead = 0;
m_CASReadsCtrl = 0;

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
if((m_hMtxIterReads = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxIterReads,NULL)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

#ifdef _WIN32
if((m_hMtxIterNxtProcRead = CreateMutex(NULL,false,NULL))==NULL)
	{
#else
if(pthread_mutex_init (&m_hMtxIterNxtProcRead,NULL)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
	pthread_mutex_destroy(&m_hMtxIterReads);
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

#ifdef _WIN32
if(!InitializeCriticalSectionAndSpinCount(&m_hSCritSectSeqHdrs,1000))
	{
#else
if(pthread_spin_init(&m_hSpinLockSeqHdrs,PTHREAD_PROCESS_PRIVATE)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
	pthread_mutex_destroy(&m_hMtxIterReads);
	pthread_mutex_destroy(&m_hMtxIterNxtProcRead);
#endif
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create mutex");
	return(eBSFerrInternal);
	}

#ifdef _WIN32
if(!InitializeCriticalSectionAndSpinCount(&m_hSCritSectSeqFlags,1000))
	{
#else
if(pthread_spin_init(&m_hSpinLockSeqFlags,PTHREAD_PROCESS_PRIVATE)!=0)
	{
	pthread_rwlock_destroy(&m_hRwLock);
	pthread_mutex_destroy(&m_hMtxIterReads);
	pthread_mutex_destroy(&m_hMtxIterNxtProcRead);
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
#ifdef _WIN32
	CloseHandle(m_hMtxIterReads);
#else
	pthread_rwlock_destroy(&m_hRwLock);
	pthread_mutex_destroy(&m_hMtxIterReads);
	pthread_mutex_destroy(&m_hMtxIterNxtProcRead);
	pthread_spin_destroy(&m_hSpinLockSeqHdrs);
#endif
	return(eBSFerrInternal);
	}

m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CKangadna::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
#ifdef _WIN32
CloseHandle(m_hMtxIterReads);
CloseHandle(m_hMtxIterNxtProcRead);
CloseHandle(m_hMtxMHReads);
DeleteCriticalSection(&m_hSCritSectSeqHdrs);

#else
pthread_mutex_destroy(&m_hMtxIterReads);
pthread_mutex_destroy(&m_hMtxIterNxtProcRead);
pthread_mutex_destroy(&m_hMtxMHReads);
pthread_spin_destroy(&m_hSpinLockSeqHdrs);
pthread_rwlock_destroy(&m_hRwLock);

#endif
m_bMutexesCreated = false;
}

void
CKangadna::AcquireSerialise(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
}

void
CKangadna::ReleaseSerialise(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
}

void
CKangadna::AcquireSerialiseNxtProcRead(void)
{
int SpinCnt = 1000;
int BackoffMS = 5;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_CASNxtProcRead,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#else
while(__sync_val_compare_and_swap(&m_CASNxtProcRead,0,1)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#endif
}

void
CKangadna::ReleaseSerialiseNxtProcRead(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASNxtProcRead,0,1);
#else
__sync_val_compare_and_swap(&m_CASNxtProcRead,1,0);
#endif
}

void
CKangadna::AcquireSerialiseReadsCtrl(void)
{
int SpinCnt = 1000;
int BackoffMS = 5;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_CASReadsCtrl,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#else
while(__sync_val_compare_and_swap(&m_CASReadsCtrl,0,1)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#endif
}

void
CKangadna::ReleaseSerialiseReadsCtrl(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASReadsCtrl,0,1);
#else
__sync_val_compare_and_swap(&m_CASReadsCtrl,1,0);
#endif
}

void
CKangadna::AcquireSerialiseSeqHdr(void)
{
int SpinCnt = 500;
int BackoffMS = 5;

#ifdef _WIN32
while(!TryEnterCriticalSection(&m_hSCritSectSeqHdrs))
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
	if(BackoffMS < 500)
		BackoffMS += 5;
	}
#else
while(pthread_spin_trylock(&m_hSpinLockSeqHdrs)==EBUSY)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
	if(BackoffMS < 500)
		BackoffMS += 5;
	}
#endif
}

void
CKangadna::ReleaseSerialiseSeqHdr(void)
{
#ifdef _WIN32
LeaveCriticalSection(&m_hSCritSectSeqHdrs);
#else
pthread_spin_unlock(&m_hSpinLockSeqHdrs);
#endif
}


void
CKangadna::AcquireSerialiseSeqFlags(void)
{
int SpinCnt = 500;
int BackoffMS = 5;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_CASSeqFlags,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#else
while(__sync_val_compare_and_swap(&m_CASSeqFlags,0,1)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
	if(BackoffMS < 500)
		BackoffMS += 2;
	}
#endif
}

void
CKangadna::ReleaseSerialiseSeqFlags(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASSeqFlags,0,1);
#else
__sync_val_compare_and_swap(&m_CASSeqFlags,1,0);
#endif
}

void 
CKangadna::AcquireLock(bool bExclusive)
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
CKangadna::ReleaseLock(bool bExclusive)
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

bool
CKangadna::SetMaxMemWorkSetSize(size_t Bytes)
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

void
CKangadna::SetCtgDescr(char *pszCtgDescr)
{
strncpy(m_szCtgDescr,pszCtgDescr,sizeof(m_szCtgDescr)-1);
m_szCtgDescr[sizeof(m_szCtgDescr)-1] = '\0';
}

UINT32											// total number of reads accepted for processing into next phase
CKangadna::GetNumReads(UINT32 *pNumPE1Reads,	// returned total number of single ended or 5' paired end read sequences accepted parsed
			UINT32 *pNumPE2Reads,				// returned total number of 3' paired end read sequences accepted parsed
			UINT32 *pTotPE1Seqs,				// number of PE1 sequences remaining at end of each phase completion
			UINT32 *pTotPE2Seqs,				// number of PE2 sequences remaining at end of each phase completion
			UINT32 *pTotSeqsParsed,				// returned total number of sequences parsed for 3' and 5' combined
			UINT32 *pTotSeqsUnderLen,			// returned total number of under length sequences parsed for 3' and 5' combined
			UINT32 *pTotSeqsExcessNs,			// returned total number of under length sequences parsed for 3' and 5' combined
			UINT32 *pMeanSeqLen,				// returned mean length (rounded up) of all sequences
			UINT32 *pMinSeqLen,					// returned minimum sequence length
			UINT32 *pMaxSeqLen)					// returned maximum sequence length
{
if(pNumPE1Reads != NULL)
	*pNumPE1Reads = m_Sequences.NumPE1Seqs2Assemb;

if(pNumPE2Reads != NULL)
	*pNumPE2Reads = m_Sequences.NumPE2Seqs2Assemb;

if(pTotPE1Seqs != NULL)
	{
	if(m_Sequences.NumPE2Seqs2Assemb)
		*pTotPE1Seqs = m_Sequences.NumSeqs2Assemb/2;
	else
		*pTotPE1Seqs = m_Sequences.NumSeqs2Assemb;
	}

if(pTotPE2Seqs != NULL)
	{
	if(m_Sequences.NumPE2Seqs2Assemb)
		*pTotPE2Seqs = m_Sequences.NumSeqs2Assemb/2;
	else
		*pTotPE2Seqs = 0;
	}

if(pTotSeqsParsed != NULL)
	*pTotSeqsParsed = m_Sequences.TotSeqsParsed;

if(pTotSeqsUnderLen != NULL)
	*pTotSeqsUnderLen = m_Sequences.TotSeqsUnderLen;

if(pTotSeqsExcessNs != NULL)
	*pTotSeqsExcessNs = m_Sequences.TotSeqsExcessNs;


if(pMeanSeqLen != NULL)
	*pMeanSeqLen = (UINT32)(ceil(m_Sequences.MeanSeqLen));

if(pMinSeqLen != NULL)
	*pMinSeqLen = m_Sequences.MinSeqLen;

if(pMaxSeqLen != NULL)
	*pMaxSeqLen = m_Sequences.MaxSeqLen;

return(m_Sequences.NumSeqs2Assemb);
}

teBSFrsltCodes
CKangadna::DumpHeader(char *pszTypeSeqFile)
{
char szHeader[0x07fff];
int BuffIdx;
UINT32 Idx;
tsPPCRdsFileHdr PPCRdsHdr;
int hInSeqTypesFile;

hInSeqTypesFile = open(pszTypeSeqFile,O_READSEQ);
if(hInSeqTypesFile == -1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"DumpHeader: Unable to open %s - %s",pszTypeSeqFile,strerror(errno));
	Reset(false);
	return(eBSFerrOpnFile);
	}

if(_lseeki64(hInSeqTypesFile,0,SEEK_SET)!=0)			
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"DumpHeader: Seek failed to offset 0 - %s",pszTypeSeqFile,strerror(errno));
	close(hInSeqTypesFile);
	return(eBSFerrFileAccess);
	}

// read in header..
if(cSizeofPPCRdsFileHdr != read(hInSeqTypesFile,&PPCRdsHdr,cSizeofPPCRdsFileHdr))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"DumpHeader: Failed reading file header - %s",pszTypeSeqFile);
	close(hInSeqTypesFile);
	return(eBSFerrNotBioseq);
	}

// header read, validate it as being sparse suffix file header
if(tolower(PPCRdsHdr.Magic[0]) != 'p' ||
	tolower(PPCRdsHdr.Magic[1]) != 'r' ||
	tolower(PPCRdsHdr.Magic[2]) != 'd' ||
	tolower(PPCRdsHdr.Magic[3]) != 's')
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"DumpHeader: File header does not have magic signiture of 'prds' - %s",pszTypeSeqFile);
	close(hInSeqTypesFile);
	return(eBSFerrNotBioseq);
	}

	// can we handle this version?
if(PPCRdsHdr.Version < cPPCRdsFileVersionBack || PPCRdsHdr.Version > cPPCRdsFileVersion)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"%s opened as a pre-processed concatenated reads file - expected between version %d and %d, file version %d",
					pszTypeSeqFile,cPPCRdsFileVersionBack,cPPCRdsFileVersion,PPCRdsHdr.Version);
	close(hInSeqTypesFile);			// closes opened file..
	return(eBSFerrFileVer);
	}

// quick check that file is of expected size
if(_lseeki64(hInSeqTypesFile,PPCRdsHdr.FileSize,SEEK_SET) != PPCRdsHdr.FileSize)			
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"DumpHeader: Expected file to be of size %lld but seek failed '%s' - %s",PPCRdsHdr.FileSize,pszTypeSeqFile,strerror(errno));
	close(hInSeqTypesFile);			// closes opened file..
	return(eBSFerrFileAccess);
	}

int DiagLevel = gDiagnostics.GetFileDiagLevel();
if(DiagLevel >= eDLDiag)
	{
	BuffIdx = sprintf(szHeader,"\tHeader length: %d\n",cSizeofPPCRdsFileHdr);
	BuffIdx += sprintf(&szHeader[BuffIdx],"\tMagic [0] %c [1] %c [2] %c [3] %c\n", (char)PPCRdsHdr.Magic[0], (char)PPCRdsHdr.Magic[1],(char)PPCRdsHdr.Magic[2],(char)PPCRdsHdr.Magic[3]);
	BuffIdx += sprintf(&szHeader[BuffIdx],"\tVersion: %u, FileSize: %llu, NumRawFiles: %d\n",PPCRdsHdr.Version,PPCRdsHdr.FileSize,PPCRdsHdr.NumRawFiles);
	for(Idx = 0; Idx < PPCRdsHdr.NumRawFiles; Idx++)
		{
		BuffIdx += sprintf(&szHeader[BuffIdx],"\tFileID: %d, PEFileID: %d, SeqsParsed: %u, SeqsUnderLen: %u, SeqsExcessNs: %u, SeqsAccepted: %u\n",
											PPCRdsHdr.SrcFiles[Idx].FileID,PPCRdsHdr.SrcFiles[Idx].PEFileID,PPCRdsHdr.SrcFiles[Idx].SeqsParsed,PPCRdsHdr.SrcFiles[Idx].SeqsUnderLen,PPCRdsHdr.SrcFiles[Idx].SeqsExcessNs,PPCRdsHdr.SrcFiles[Idx].SeqsAccepted);

		BuffIdx += sprintf(&szHeader[BuffIdx],"\tFile: '%s'\n",	PPCRdsHdr.SrcFiles[Idx].szFile);									
		}
	BuffIdx += sprintf(&szHeader[BuffIdx],"\tSfxSparsity: %u, SeqWrdBytes: %u, SeqHdrLen: %d, TotSeqsParsed: %u, TotSeqsUnderLen: %u, TotSeqsExcessNs: %u\n",
										PPCRdsHdr.Sequences.SfxSparsity,PPCRdsHdr.Sequences.SeqWrdBytes,PPCRdsHdr.Sequences.SeqHdrLen,PPCRdsHdr.Sequences.TotSeqsParsed,PPCRdsHdr.Sequences.TotSeqsUnderLen,PPCRdsHdr.Sequences.TotSeqsExcessNs);
	printf("\nFile Header Dump\n%s",szHeader);

	printf("\n\tSeqWrdBytes: %d, SeqHdrLen: %d\n\tTotSeqs2Assemb: %u, NumSeqs2Assemb: %u, Seqs2AssembLen: %lld\n\tMeanSeqLen: %1.1f, MinSeqLen: %d, MaxSeqLen: %d\n\n",
			PPCRdsHdr.Sequences.SeqWrdBytes,	PPCRdsHdr.Sequences.SeqHdrLen,
			PPCRdsHdr.Sequences.TotSeqs2Assemb,PPCRdsHdr.Sequences.NumSeqs2Assemb, PPCRdsHdr.Sequences.Seqs2AssembLen,
			PPCRdsHdr.Sequences.MeanSeqLen,PPCRdsHdr.Sequences.MinSeqLen,PPCRdsHdr.Sequences.MaxSeqLen);
	}
UINT32 NumSeqs;
NumSeqs = PPCRdsHdr.Sequences.NumSeqs2Assemb;
if(PPCRdsHdr.Sequences.bPESeqs)
	NumSeqs /= 2;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"DumpHeader: Sequences: %s %u, MeanSeqLen: %1.1f, MinSeqLen: %d, MaxSeqLen: %d",PPCRdsHdr.Sequences.bPESeqs ? "Paired end" : "Single ended",
									NumSeqs,PPCRdsHdr.Sequences.MeanSeqLen,PPCRdsHdr.Sequences.MinSeqLen,PPCRdsHdr.Sequences.MaxSeqLen);
return(eBSFSuccess);
}

teBSFrsltCodes			// load saved concatenated and packed sequences from file 
CKangadna::LoadPackedSeqsFromFile(char *pszTypeSeqFile)	// loading is from this file
{
teBSFrsltCodes Rslt;
tsPPCRdsFileHdr PPCRdsHdr;

ResetTypeSeqs();			// ensure any memory previously allocated will be freed

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadPackedSeqsFromFile: Loading artefact reduced packed reads from file: '%s'",pszTypeSeqFile);
if((Rslt = DumpHeader(pszTypeSeqFile))!= eBSFSuccess)
	return(Rslt);

m_hInSeqTypesFile = open(pszTypeSeqFile,O_READSEQ);
if(m_hInSeqTypesFile == -1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadPackedSeqsFromFile: Unable to open %s - %s",pszTypeSeqFile,strerror(errno));
	Reset(false);
	return(eBSFerrOpnFile);
	}

if(_lseeki64(m_hInSeqTypesFile,0,SEEK_SET)!=0)			
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadPackedSeqsFromFile: Seek failed to offset 0 - %s",pszTypeSeqFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

// read in header..
if(sizeof(tsPPCRdsFileHdr) != read(m_hInSeqTypesFile,&PPCRdsHdr,sizeof(tsPPCRdsFileHdr)))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadPackedSeqsFromFile: Failed reading file header - %s",pszTypeSeqFile);
	Reset(false);
	return(eBSFerrNotBioseq);
	}

// header read, validate it as being sparse suffix file header
if(tolower(PPCRdsHdr.Magic[0]) != 'p' ||
	tolower(PPCRdsHdr.Magic[1]) != 'r' ||
	tolower(PPCRdsHdr.Magic[2]) != 'd' ||
	tolower(PPCRdsHdr.Magic[3]) != 's')
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadPackedSeqsFromFile: File header does not have magic signature of 'prds' - %s",pszTypeSeqFile);
	Reset(false);
	return(eBSFerrNotBioseq);
	}

	// can we handle this version?
if(PPCRdsHdr.Version < cPPCRdsFileVersionBack || PPCRdsHdr.Version > cPPCRdsFileVersion)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"%s opened as a packed sequence file - expected between version %d and %d, file version %d",
					pszTypeSeqFile,cPPCRdsFileVersionBack,cPPCRdsFileVersion,PPCRdsHdr.Version);
	Reset(false);			// closes opened file..
	return(eBSFerrFileVer);
	}

// quick check that file is of expected size
if(_lseeki64(m_hInSeqTypesFile,PPCRdsHdr.FileSize,SEEK_SET) != PPCRdsHdr.FileSize)			
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadPackedSeqsFromFile: Expected file to be of size %lld but seek failed '%s' - %s",PPCRdsHdr.FileSize,pszTypeSeqFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

if(PPCRdsHdr.Sequences.NumSeqs2Assemb == 0)			// surely must be at least one sequence to be assembled!!
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"LoadPackedSeqsFromFile: Nothing to do, file '%s' contains no sequences",pszTypeSeqFile);
	Reset(false);			// closes opened file..
	return(eBSFerrNoEntries);
	}

memmove(&m_Sequences,&PPCRdsHdr.Sequences,sizeof(m_Sequences));
SeqWrdBytes = m_Sequences.SeqWrdBytes;
m_NumRawFiles = PPCRdsHdr.NumRawFiles;
memmove(&m_SrcFiles,&PPCRdsHdr.SrcFiles,sizeof(m_SrcFiles));


UINT64 CurWorkSetSize = 0;
if(m_Sequences.OfsSeqs2Assemb && m_Sequences.AllocMemSeqs2Assemb)
	CurWorkSetSize = m_Sequences.AllocMemSeqs2Assemb;

if(CurWorkSetSize != m_CurMaxMemWorkSetBytes)
	{
	if(!SetMaxMemWorkSetSize((size_t)CurWorkSetSize))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadPackedSeqsFromFile: Unable to set maximum working set size (bytes) to %lld bytes",m_CurMaxMemWorkSetBytes);
		Reset(false);
		return(eBSFerrMaxDirEls);
		}
	}



m_Sequences.AllocMemSeqStarts = 0;
m_Sequences.AllocMemSeqs2Assemb = 0;
m_Sequences.AllocMemSfx = 0;
m_Sequences.pSeqStarts = NULL;
m_Sequences.pSeqs2Assemb = NULL;
m_Sequences.pSuffixArray = NULL;
m_Sequences.pSeqFlags = NULL;
m_Sequences.AllocMemSeqFlags = 0;
Rslt = eBSFSuccess;
if(PPCRdsHdr.Sequences.OfsSeqs2Assemb && PPCRdsHdr.Sequences.AllocMemSeqs2Assemb)
	{
	if((Rslt = AllocLoadBlock(pszTypeSeqFile,PPCRdsHdr.Sequences.OfsSeqs2Assemb,PPCRdsHdr.Sequences.AllocMemSeqs2Assemb,&m_Sequences.pSeqs2Assemb,&m_Sequences.AllocMemSeqs2Assemb))!=eBSFSuccess)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadPackedSeqsFromFile: Loading file %s failed",pszTypeSeqFile);
	else
		m_Sequences.bSeqs2AssembDirty = true;			// force suffix indexes to be re-generated
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadPackedSeqsFromFile: Completed loading artefact reduced packed reads from file: '%s'",pszTypeSeqFile);
return(Rslt);
}


teBSFrsltCodes 
CKangadna::AllocLoadBlock(char *pszInFile, // loading from this file
					 UINT64 FileOfs,				// load from this file offset
					 UINT64  AllocBlockSize,		// load, and allocate for, this block size from disk
					 void **ppLoadedBlock,		// returned ptr to allocated memory
					 UINT64 *pAllocBlockSize)	// size of allocated memory
{
teBSFrsltCodes Rslt;

if(FileOfs == 0 || AllocBlockSize == 0 || ppLoadedBlock == NULL || pAllocBlockSize == NULL || m_hInSeqTypesFile == -1)
	return(eBSFerrParams);

*ppLoadedBlock = NULL;
*pAllocBlockSize = 0;

	// now try and allocate memory
#ifdef _WIN32
*ppLoadedBlock = (void *) malloc((size_t)AllocBlockSize);	
if(*ppLoadedBlock == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocLoadBlock: array memory allocation of %llu bytes - %s",AllocBlockSize,strerror(errno));
	Reset(false);
	return(eBSFerrMem);
	}
#else
if((*ppLoadedBlock = (void *)mmap(NULL,AllocBlockSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocLoadBlock: array memory allocation of %llu bytes through mmap()  failed - %s",AllocBlockSize,strerror(errno));
	*ppLoadedBlock = NULL;
	*pAllocBlockSize = 0;
	Reset(false);
	return(eBSFerrMem);
	}
#endif
*pAllocBlockSize = AllocBlockSize;
memset(*ppLoadedBlock,0,(size_t)AllocBlockSize); // commits the memory!
	// memory allocated, now initialise from file
if((Rslt = ChunkedRead(m_hInSeqTypesFile,pszInFile,FileOfs,(UINT8 *)*ppLoadedBlock,AllocBlockSize)) != eBSFSuccess)
	Reset(false);
return(Rslt);
}

teBSFrsltCodes			
CKangadna::SavePackedSeqsToFile(char *pszTypeSeqFile)	// save to this file
{
teBSFrsltCodes Rslt;
tsPPCRdsFileHdr PPCRdsHdr;

if(m_PMode == 0)
	return(SaveAsFasta(pszTypeSeqFile));

#ifdef _WIN32
if((m_hOutSeqTypesFile = open(pszTypeSeqFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutSeqTypesFile = open(pszTypeSeqFile, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",pszTypeSeqFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output preprocessed reads packed sequences file created/truncated: '%s'",pszTypeSeqFile);

if(_lseeki64(m_hOutSeqTypesFile,0,SEEK_SET)!=0)			
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Seek failed to offset 0 - %s",pszTypeSeqFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

// initialise header..
memset(&PPCRdsHdr,0,sizeof(PPCRdsHdr));
PPCRdsHdr.Magic[0] = 'p';
PPCRdsHdr.Magic[1] = 'r';
PPCRdsHdr.Magic[2] = 'd';
PPCRdsHdr.Magic[3] = 's';
PPCRdsHdr.Version = cPPCRdsFileVersion;
memmove(&PPCRdsHdr.Sequences,&m_Sequences,sizeof(PPCRdsHdr.Sequences));
PPCRdsHdr.NumRawFiles = m_NumRawFiles;
memmove(&PPCRdsHdr.SrcFiles,&m_SrcFiles,sizeof(m_SrcFiles));

PPCRdsHdr.FileSize = sizeof(PPCRdsHdr);

if(m_Sequences.pSeqs2Assemb != NULL && m_Sequences.Seqs2AssembOfs > 0 && m_Sequences.pSeqs2Assemb != NULL)
	{
	UpdateAllSeqHeaderFlags(0,~(cFlgSeqPE | cFlgSeqPE2),false);
	PPCRdsHdr.Sequences.AllocMemSeqs2Assemb = (m_Sequences.Seqs2AssembOfs + 16) * m_Sequences.SeqWrdBytes; // could have overallocated...
	if(PPCRdsHdr.Sequences.AllocMemSeqs2Assemb > m_Sequences.AllocMemSeqs2Assemb)
		PPCRdsHdr.Sequences.AllocMemSeqs2Assemb = m_Sequences.AllocMemSeqs2Assemb;
	PPCRdsHdr.Sequences.OfsSeqs2Assemb = PPCRdsHdr.FileSize;
	PPCRdsHdr.FileSize += PPCRdsHdr.Sequences.AllocMemSeqs2Assemb;
	}

if((Rslt = ChunkedWrite(m_hOutSeqTypesFile,pszTypeSeqFile,0,(UINT8 *)&PPCRdsHdr,sizeof(PPCRdsHdr)))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"SaveTypeSeqsToFile: Write header to file %s failed",pszTypeSeqFile);
	return(Rslt);
	}

if(m_Sequences.pSeqs2Assemb != NULL && m_Sequences.Seqs2AssembOfs > 0 && m_Sequences.pSeqs2Assemb != NULL)
	if((Rslt = ChunkedWrite(m_hOutSeqTypesFile,pszTypeSeqFile,PPCRdsHdr.Sequences.OfsSeqs2Assemb,(UINT8 *)m_Sequences.pSeqs2Assemb,PPCRdsHdr.Sequences.AllocMemSeqs2Assemb))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"SaveTypeSeqsToFile: Write to file %s failed",pszTypeSeqFile);
		return(Rslt);
		}

#ifdef _WIN32
_commit(m_hOutSeqTypesFile);
#else
fsync(m_hOutSeqTypesFile);
#endif
close(m_hOutSeqTypesFile);
m_hOutSeqTypesFile = -1;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"SaveTypeSeqsToFile: Completed with out errror");
return(eBSFSuccess);
}

teBSFrsltCodes
CKangadna::AllocSeqs2AssembMem(UINT64 ReqAllocSize)		// allocate to hold at least this many bytes
{
size_t SizeT = (size_t)ReqAllocSize;
if(ReqAllocSize < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqs2AssembMem: ReqAllocSize < 1");
	Reset(false);
	return(eBSFerrParams);
	}
if((UINT64)SizeT != ReqAllocSize)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqs2AssembMem: Requested allocation of %lld too large for 32bit application",ReqAllocSize);
	Reset(false);
	return(eBSFerrParams);
	}
ReqAllocSize = (ReqAllocSize + 3) & (UINT64)0x0ffffffffffffc;						// ensure allocation is for integral number of tSeqWrd4's

if(m_Sequences.pSeqs2Assemb != NULL && (m_Sequences.AllocMemSeqs2Assemb < ReqAllocSize || ((m_Sequences.AllocMemSeqs2Assemb * 10) > (ReqAllocSize * 12))))
	{
#ifdef _WIN32
	free(m_Sequences.pSeqs2Assemb);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_Sequences.pSeqs2Assemb != MAP_FAILED)
		munmap(m_Sequences.pSeqs2Assemb,m_Sequences.AllocMemSeqs2Assemb);
#endif	
	m_Sequences.pSeqs2Assemb = NULL;
	m_Sequences.AllocMemSeqs2Assemb = 0;
	}

if(m_Sequences.pSeqs2Assemb == NULL)
	{
	m_Sequences.AllocMemSeqs2Assemb = (size_t)ReqAllocSize;
#ifdef _WIN32
	m_Sequences.pSeqs2Assemb = malloc((size_t)m_Sequences.AllocMemSeqs2Assemb);	
	if(m_Sequences.pSeqs2Assemb == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqs2AssembMem: Concatenated packed sequences memory allocation of %llu bytes - %s",m_Sequences.AllocMemSeqs2Assemb,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}
#else
	if((m_Sequences.pSeqs2Assemb = (void *)mmap(NULL,m_Sequences.AllocMemSeqs2Assemb, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqs2AssembMem: Concatenated packed sequences memory of %llu bytes through mmap()  failed - %s",m_Sequences.AllocMemSeqs2Assemb,strerror(errno));
		m_Sequences.pSeqs2Assemb = NULL;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	}

memset(m_Sequences.pSeqs2Assemb,0,(size_t)m_Sequences.AllocMemSeqs2Assemb);	// commits the memory!
m_Sequences.Seqs2AssembLen = 0;
m_Sequences.Seqs2AssembOfs = 0;
m_Sequences.NumSeqs2Assemb = 0;
if(m_Sequences.SfxSparsity == 0 || m_Sequences.SeqWrdBytes == 0)
	{
	m_Sequences.SfxSparsity = 15;
	m_Sequences.SeqWrdBytes = 4;
	SeqWrdBytes = m_Sequences.SeqWrdBytes;
	m_Sequences.SeqHdrLen = 3;
	}
return(eBSFSuccess);
}


teBSFrsltCodes
CKangadna::AllocBlockNsLoci(UINT32 ReqAllocBlocks)		// alloc/realloc to at least ReqAllocSize (tsBlockNsLoci)
{
size_t SizeT = (size_t)ReqAllocBlocks * sizeof(tsBlockNsLoci);
if(ReqAllocBlocks < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocBlockNsLoci: ReqAllocBlocks < 1");
	Reset(false);
	return(eBSFerrParams);
	}

if(m_pBlockNsLoci != NULL && (m_AllocdBlockNsLoci < ReqAllocBlocks || ((m_AllocdBlockNsLoci * 10) > (ReqAllocBlocks * 12))))
	{
#ifdef _WIN32
	free(m_pBlockNsLoci);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pBlockNsLoci != MAP_FAILED)
		munmap(m_pBlockNsLoci,m_AllocdBlockNsLociSize);
#endif	
	m_pBlockNsLoci = NULL;
	m_AllocdBlockNsLociSize = 0;
	m_AllocdBlockNsLoci = 0;
	}

if(m_pBlockNsLoci == NULL)
	{
	m_AllocdBlockNsLociSize = SizeT;
	m_AllocdBlockNsLoci = ReqAllocBlocks;
#ifdef _WIN32
	m_pBlockNsLoci = (tsBlockNsLoci *)malloc(m_AllocdBlockNsLociSize);	
	if(m_pBlockNsLoci == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocBlockNsLoci: memory allocation of %llu bytes - %s",m_AllocdBlockNsLociSize,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}
#else
	if((m_pBlockNsLoci = (tsBlockNsLoci *)mmap(NULL,m_AllocdBlockNsLociSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocBlockNsLoci: memory of %llu bytes through mmap()  failed - %s",m_AllocdBlockNsLociSize,strerror(errno));
		m_pBlockNsLoci = NULL;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	}

memset(m_pBlockNsLoci,0,m_AllocdBlockNsLociSize);	// commits the memory!
m_NumBlockNsLoci = 0;
return(eBSFSuccess);
}

UINT32
CKangadna::GetEstSeqsToProc(void)  
{
if(m_SeqEsts.bCalcMeanSeqLen)
	return(0);
return(m_SeqEsts.NumSeqsToProc);
}

tsEstSeqs * 
CKangadna::GetEstMemReq(int Trim5,				// will be trimming this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// will be trimming trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int TrimSeqLen,			// will be trimming sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...10000)
					int SampleNth,			// will only be processing every Nth reads
					int Zreads)				// will only accept this number of reads for processing from any file
{
UINT64 MemReq;

MemReq = EstMemReq(Trim5,Trim3,TrimSeqLen,SampleNth,Zreads);
return(MemReq == 0 ? NULL : &m_SeqEsts);
}

teBSFrsltCodes 
CKangadna::EstMemReq(char *pszInFile)      // accumulate estimates from this file
{
CFasta Fasta;
UINT32 NumSeqs;
INT32 DescrLen;
INT32 SeqLen;

if(pszInFile == NULL || pszInFile[0] == '\0')
	return(eBSFerrParams);

// get estimate of number of sequences, mean sequence length, and mean descriptor length
if((NumSeqs = Fasta.FastaEstSizes(pszInFile,NULL,NULL,&DescrLen,NULL,&SeqLen)) == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory requirements for file '%s'",pszInFile);
	return(eBSFerrOpnFile);
	}

m_SeqEsts.bCalcMeanSeqLen = true;
m_SeqEsts.NumFiles += 1;
m_SeqEsts.TotSeqLens += (UINT64)SeqLen * (NumSeqs+10);	// 10 is just as a slight safety margin
m_SeqEsts.TotNumSeqs += NumSeqs;
m_SeqEsts.NumSeqsToProc = m_SeqEsts.TotNumSeqs;
return(eBSFSuccess);
}


// estimates (with a 5% safety additional)  how much memory needs to be allocated to hold packed sequences including their headers
UINT64 
CKangadna::EstMemReq(int Trim5,				// will be trimming this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// will be trimming trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int TrimSeqLen,			// will be trimming sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...10000)
					int SampleNth,			// will only be processing every Nth reads
					int Zreads)				// will only accept this number of reads for processing from any file
	
{
UINT32 NumSeqWords;
UINT64 MemReq;
UINT32 MemPerSeq;
UINT32 MeanLenTrimmed;
UINT64 NumSeqsToProc;

if(m_SeqEsts.NumFiles == 0 || m_SeqEsts.TotSeqLens == 0 || m_SeqEsts.TotNumSeqs == 0)
	return(0);

if(m_SeqEsts.bCalcMeanSeqLen == true)
	{
	m_SeqEsts.MeanSeqLen = MeanLenTrimmed = (UINT32)((m_SeqEsts.TotSeqLens +  m_SeqEsts.TotNumSeqs - 1)  / (UINT64)m_SeqEsts.TotNumSeqs);
	m_SeqEsts.bCalcMeanSeqLen = false;
	}
else
	MeanLenTrimmed = m_SeqEsts.MeanSeqLen;

if(MeanLenTrimmed < (UINT32)(Trim5 + Trim3 + 10))	// ensure will have a mean length of at least 10bp after trimming
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Est. mean read length too short after end trimming...");
	return(0);
	}

// adjust sequence mean length for any trimming required
MeanLenTrimmed -= (Trim5 - Trim3);
if(TrimSeqLen && MeanLenTrimmed > (UINT32)TrimSeqLen)
	MeanLenTrimmed = (UINT32)TrimSeqLen;

// calc number bytes required per sequence
// packed 15 bases per 32bit tSeqWrd4 with 12 byte header
NumSeqWords = (MeanLenTrimmed+14)/15;	// number of 4byte words
MemPerSeq = NumSeqWords * 4;			// back to bytes
MemPerSeq += 12;						// plus header

if(Zreads)
	{
	NumSeqsToProc = Zreads * m_SeqEsts.NumFiles;
	if(NumSeqsToProc > m_SeqEsts.TotNumSeqs)
		NumSeqsToProc = m_SeqEsts.TotNumSeqs;
	}
else
	if(SampleNth > 1)
		NumSeqsToProc = (m_SeqEsts.TotNumSeqs + SampleNth - 1) / SampleNth;
	else
		NumSeqsToProc = m_SeqEsts.TotNumSeqs;

m_SeqEsts.NumSeqsToProc = (UINT32)NumSeqsToProc;
// add 5% for a little safety margin...
NumSeqsToProc = (NumSeqsToProc * 105) / 100;
MemReq = NumSeqsToProc * MemPerSeq;
MemReq = (MemReq + 3) & (UINT64)0x0ffffffffffffc;  // roundup to 32 bit word boundary
return(MemReq);
}

bool				// false if mean Phred score is below minimum threshold or any individual base Phred is < 10, true if Phred accepted
CKangadna::MeetsMinPhredScoreThres(int QSSchema,	// quality scoring schema - guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger,
			int MinMeanPhredScore,			// minimum allowed mean (over all read bases) Phred score
			char *pszPhredScores)			// read ascii Phred scores
{
char Phred;
char Phred0;
int Score;

if(QSSchema == 0 || MinMeanPhredScore < 20 || pszPhredScores == NULL || pszPhredScores[0] == '\0')
	return(true);

switch(QSSchema) {						// quality scoring schema - guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or Sanger 
	case 1:	// Solexa  ';' (-5) to 'h' (40)
	case 2: // Illumina 1.3+ '@' (0) to 'h' (40)
	case 3: // Illumina 1.5+ 'B' (0) to 'h' (40)
		Phred0 = '@';
		break;

	case 4: // Illumina 1.8+ or Sanger. Sanger is '!' (0) to 'J' (41) and Illumina is '#' (2) to 'J' (41)
		Phred0 = '!';
		break;
	}

int SumScores;
int NumScores;

NumScores = 0;
SumScores = 1;
while(Phred = *pszPhredScores++)
	{
	if((Score = (int)(Phred - Phred0)) < 10)	// if less than Phred 10 then this read is of dubious quality ...
		return(false);
	SumScores += Score;
	NumScores += 1;
	}
return( SumScores/NumScores >= MinMeanPhredScore ? true : false);
}


#ifdef _WIN32
unsigned __stdcall ThreadedFiltReads(void * pThreadPars)
#else
void * ThreadedFiltReads(void * pThreadPars)
#endif
{
int Rslt = 0;
tsThreadFiltReadsPars *pPars = (tsThreadFiltReadsPars *)pThreadPars; // makes it easier not having to deal with casts!
CKangadna *pThis = (CKangadna *)pPars->pThis;
Rslt = pThis->ProcReadsThread(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

teBSFrsltCodes										// returns < eBSFSuccess if any errors, eBSFSuccess if all reads already processed, 1 if SE or 2 if PE reads being returned for processing
CKangadna::GetNxtProcRead(tsThreadFiltReadsPars *pPars)		// uniquely identifies calling thread, used to manage buffers for quality scores and sequences
{
teBSFrsltCodes Rslt;

tsProcReadsCtrl *pCtrl;
int PE1DescrLen;
UINT8 szPE1DescrBuff[4096];
int PE2DescrLen;
UINT8 szPE2DescrBuff[4096];

pPars->PE1QSLen = 0;
pPars->szPE1QScoresBuff[0] = '\0';
pPars->PE1ReadLen = 0;
pPars->PE1RawReadsBuff[0] = 0;
pPars->PE2QSLen = 0;
pPars->szPE2QScoresBuff[0] = '\0';
pPars->PE2ReadLen = 0;
pPars->PE2RawReadsBuff[0] = 0;

if(pPars == NULL || (pCtrl =  pPars->pProcReadsCtrl) == NULL)
	return(eBSFerrParams);

AcquireSerialiseReadsCtrl();
if((Rslt = pCtrl->Rslt) < eBSFSuccess || pCtrl->bPE1PE2AllReadsProc || pCtrl->pFastaPE1 == NULL)
	{
	ReleaseSerialiseReadsCtrl();     
	return(eBSFSuccess);		
	}

// at least one read, or two if PE, can be dequeued
if((Rslt = (teBSFrsltCodes)(pPars->PE1ReadLen = pCtrl->pFastaPE1->ReadSequence(pPars->PE1RawReadsBuff,sizeof(pPars->PE1RawReadsBuff)-1,true,false))) > eBSFSuccess)
	{
	pPars->NumPE1ParsedReads += 1;
	pCtrl->NumPE1ParsedReads += 1;
	if(pPars->PE1ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		pCtrl->NumPE1DescrReads += 1;
		PE1DescrLen = pCtrl->pFastaPE1->ReadDescriptor((char *)szPE1DescrBuff,sizeof(szPE1DescrBuff)-1);
		szPE1DescrBuff[sizeof(szPE1DescrBuff)-1] = '\0';

		if(pCtrl->bIsPE1Fastq && pCtrl->MinPhredScore > 0)
			{
			pPars->PE1QSLen = pCtrl->pFastaPE1->ReadQValues((char *)pPars->szPE1QScoresBuff, sizeof(pPars->szPE1QScoresBuff) - 1);
			pPars->szPE1QScoresBuff[pPars->PE1QSLen] = '\0';
			}
		else
			{
			pPars->PE1QSLen = 0;
			pPars->szPE1QScoresBuff[0] = '\0';
			}


		pPars->PE1ReadLen = pCtrl->pFastaPE1->ReadSequence(pPars->PE1RawReadsBuff,sizeof(pPars->PE1RawReadsBuff)-1);
		if(pPars->PE1ReadLen < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence after %d sequences parsed from file '%s'",pCtrl->NumPE1DescrReads,pCtrl->szPE1File);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE1DescrBuff);
			pCtrl->pFastaPE1->Close();
			if(pCtrl->bProcPE)
				pCtrl->pFastaPE2->Close();
			pCtrl->Rslt = eBSFerrParse;
			ReleaseSerialiseReadsCtrl();
			return(eBSFerrParse);
			}

		if (pPars->PE1QSLen > 0 && pPars->PE1ReadLen != pPars->PE1QSLen)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Problem parsing sequence after %d sequences parsed from file '%s'", pCtrl->NumPE1DescrReads, pCtrl->szPE1File);
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected sequence length to be same as quality length. Last descriptor parsed: %s", szPE1DescrBuff);
			pCtrl->pFastaPE1->Close();
			if(pCtrl->bProcPE)
				pCtrl->pFastaPE2->Close();
			pCtrl->Rslt = eBSFerrParse;
			ReleaseSerialiseReadsCtrl();
			return(eBSFerrParse);
			}

		if(!pCtrl->bProcPE)
			{
			pCtrl->CurTotPE1ReadsParsed += 1;
			ReleaseSerialiseReadsCtrl();
			return((teBSFrsltCodes)1);			// returning 1 SE read
			}

		Rslt = (teBSFrsltCodes)(pPars->PE2ReadLen = pCtrl->pFastaPE2->ReadSequence(pPars->PE2RawReadsBuff,sizeof(pPars->PE2RawReadsBuff)-1,true,false));
		if(Rslt <= eBSFSuccess)
			{
			if(Rslt < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors processing file: %s ",pCtrl->szPE2File);
				while(pCtrl->pFastaPE2->NumErrMsgs())
					gDiagnostics.DiagOut(eDLFatal,gszProcName,pCtrl->pFastaPE2->GetErrMsg());
				}
			else
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Insufficent reads in file: %s ",pCtrl->szPE2File);
			pCtrl->pFastaPE1->Close();
			pCtrl->pFastaPE2->Close();
			pCtrl->Rslt = Rslt;
			ReleaseSerialiseReadsCtrl();
			return(Rslt);
			}
		pPars->NumPE2ParsedReads += 1;
		pCtrl->NumPE2ParsedReads += 1;
		if(pPars->PE2ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
			{
			pCtrl->NumPE2DescrReads += 1;
			PE2DescrLen = pCtrl->pFastaPE2->ReadDescriptor((char *)szPE2DescrBuff,sizeof(szPE2DescrBuff)-1);
			szPE2DescrBuff[sizeof(szPE2DescrBuff)-1] = '\0'; 

			if(pCtrl->bIsPE2Fastq && pCtrl->MinPhredScore > 0)
				{
				pPars->PE2QSLen = pCtrl->pFastaPE2->ReadQValues((char *)pPars->szPE2QScoresBuff, sizeof(pPars->szPE2QScoresBuff)-1);
				pPars->szPE2QScoresBuff[pPars->PE2QSLen] = '\0';
				}
			else
				{
				pPars->PE2QSLen = 0;
				pPars->szPE2QScoresBuff[0] = '\0';
				}
			pPars->PE2ReadLen = pCtrl->pFastaPE2->ReadSequence(pPars->PE2RawReadsBuff,sizeof(pPars->PE2RawReadsBuff) - 1);
			if(pPars->PE2ReadLen < 1)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence after %d sequences parsed from file '%s'",pCtrl->NumPE2DescrReads,pCtrl->szPE2File);
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE2DescrBuff);
				pCtrl->pFastaPE1->Close();
				pCtrl->pFastaPE2->Close();
				pCtrl->Rslt = eBSFerrParse;
				ReleaseSerialiseReadsCtrl();
				return(eBSFerrParse);
				}

			if (pPars->PE2QSLen > 0 && pPars->PE2ReadLen != pPars->PE2QSLen)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Problem parsing sequence after %d sequences parsed from file '%s'", pCtrl->NumPE2DescrReads, pCtrl->szPE2File);
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected sequence length to be same as quality length. Last descriptor parsed: %s", szPE2DescrBuff);
				pCtrl->pFastaPE1->Close();
				pCtrl->pFastaPE2->Close();
				pCtrl->Rslt = eBSFerrParse;
				ReleaseSerialiseReadsCtrl();
				return(eBSFerrParse);
				}
			}
		else
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"raw sequence file '%s' processing error: %s ",pCtrl->szPE2File, 
						Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
			pCtrl->pFastaPE1->Close();
			pCtrl->pFastaPE2->Close();
			pCtrl->Rslt = eBSFerrParse;
			ReleaseSerialiseReadsCtrl();
			return(eBSFerrParse);
			}
		}
	pCtrl->CurTotPE1ReadsParsed += 1;
	ReleaseSerialiseReadsCtrl();
	return((teBSFrsltCodes)2);		// returning 2 (PE1 and PE2) reads
	}
pCtrl->bPE1PE2AllReadsProc = true;
pCtrl->Rslt = Rslt;
ReleaseSerialiseReadsCtrl();
return(Rslt);
}

 
int
CKangadna::ProcReadsThread(tsThreadFiltReadsPars *pPars)
{
teBSFrsltCodes Rslt;

int ContamLen5PE1;
int ContamLen3PE1;
int ContamLen5PE2;
int ContamLen3PE2;

etSeqBase *pPE1Seq;
etSeqBase *pPE2Seq;

bool bIsPE1Fastq;
bool bIsPE2Fastq;

bool bFailsPhred;

pPars->NumContamFlankTrim = 0;
pPars->NumContamVector = 0;
pPars->NumPE1AcceptedReads = 0;
pPars->NumPE1ExcessNs = 0;
pPars->NumPE1ParsedReads = 0;
pPars->NumPE1Underlen = 0;
pPars->NumPE1UnderPhredScore = 0;
pPars->NumPE2AcceptedReads = 0;
pPars->NumPE2ExcessNs = 0;
pPars->NumPE2ParsedReads = 0;
pPars->NumPE2Underlen = 0;
pPars->NumPE2UnderPhredScore = 0;
pPars->PE1QSLen = 0;
pPars->szPE1QScoresBuff[0] = '\0';
pPars->PE1ReadLen = 0;
pPars->PE1RawReadsBuff[0] = 0;
pPars->PE2QSLen = 0;
pPars->szPE2QScoresBuff[0] = '\0';
pPars->PE2ReadLen = 0;
pPars->PE2RawReadsBuff[0] = 0;

pPars->AcceptedTotSeqLen = 0;
pPars->MaxAcceptedReadLen = 0;	 
pPars->MinAcceptedReadLen = 0;	 

// get next SE or PE to process 
while((Rslt = GetNxtProcRead(pPars)) > eBSFSuccess)				// < eBSFSuccess if errors, eBSFSuccess if all reads processed, 1 if SE, 2 if PE reads	
	{
	if(pPars->PE1QSLen > 0 && pPars->szPE1QScoresBuff[0] != '\0')
		bIsPE1Fastq = true;						// need to check on quality scores
	else
		bIsPE1Fastq = false;

	if(pPars->PE2QSLen > 0 && pPars->szPE2QScoresBuff[0] != '\0')
		bIsPE2Fastq = true;						// need to check on quality scores
	else
		bIsPE2Fastq = false;

				// if trimming any contaminant sequences then check for contaminant overlaps onto this read(s)
	if(m_pContaminants != NULL)
		{
		bool bContamVector = false;

		// currently treating any contaminant match errors as if simply there was no overlap - should really report these!!!
		if((ContamLen5PE1 = m_pContaminants->MatchContaminants(eAOF5PE1Targ,1,pPars->Trim5+1,pPars->PE1ReadLen,pPars->PE1RawReadsBuff)) <= pPars->Trim5)
			ContamLen5PE1 = 0;
		else
			{
			if(ContamLen5PE1 == pPars->PE1ReadLen)	// will be equal if totally contained
				bContamVector = true;
			else
				if(ContamLen5PE1 > pPars->Trim5)
					ContamLen5PE1 -= pPars->Trim5;
			}
		if(!bContamVector)
			{
			if((ContamLen3PE1 = m_pContaminants->MatchContaminants(eAOF3PE1Targ,1,pPars->Trim3+1,pPars->PE1ReadLen,pPars->PE1RawReadsBuff)) <= pPars->Trim3)
				ContamLen3PE1 = 0;
			else
				{
				if(ContamLen3PE1 == pPars->PE1ReadLen)	// will be equal if totally contained
					bContamVector = true;
				else
					if(ContamLen3PE1 > pPars->Trim3)
						ContamLen3PE1 -=pPars-> Trim3;
				}
			}

		if(m_Sequences.bPESeqs && !bContamVector)
			{
			if((ContamLen5PE2 = m_pContaminants->MatchContaminants(eAOF5PE2Targ,1,pPars->Trim5+1,pPars->PE2ReadLen,(etSeqBase *)pPars->PE2RawReadsBuff)) <= pPars->Trim5)
				ContamLen5PE2 = 0;
			else
				{
				if(ContamLen5PE2 == pPars->PE2ReadLen)	// will be equal if totally contained
					bContamVector = true;
				else
					if(ContamLen5PE2 > pPars->Trim5)
						ContamLen5PE2 -= pPars->Trim5;
				}
			if(!bContamVector)
				{
				if((ContamLen3PE2 = m_pContaminants->MatchContaminants(eAOF3PE2Targ,1,pPars->Trim3+1,pPars->PE2ReadLen,(etSeqBase *)pPars->PE2RawReadsBuff)) <= pPars->Trim3)
					ContamLen3PE2 = 0;
				else
					{
					if(ContamLen3PE2 == pPars->PE2ReadLen)	// will be equal if totally contained
						bContamVector = true;
					else
					if(ContamLen3PE2 > pPars->Trim3)
						ContamLen3PE2 -= pPars->Trim3;
					}
				}
			}
		else
			{
			ContamLen5PE2 = 0;
			ContamLen3PE2 = 0;
			}

		if(bContamVector)
			{
			pPars->NumContamVector += 1;
			continue;
			}
		if(ContamLen5PE1 || ContamLen5PE2)
			pPars->NumContamFlankTrim += 1;		
		if(ContamLen3PE1 || ContamLen3PE2)
			pPars->NumContamFlankTrim += 1;
		}
	else
		{
		ContamLen5PE1 = 0;
		ContamLen3PE1 = 0;
		ContamLen5PE2 = 0;
		ContamLen3PE2 = 0;
		}

	// ensure would still have a sequence of at least cMinSeqLen after any end trims were applied
	if((pPars->Trim5 + pPars->Trim3 + ContamLen5PE1 + ContamLen3PE1 + (int)pPars->MinSeqLen) > (int)pPars->PE1ReadLen)
		{
		pPars->NumPE1Underlen += 1;
		continue;
		}

	if(m_Sequences.bPESeqs)
		{
		if((pPars->Trim5 + pPars->Trim3 + ContamLen5PE2 + ContamLen3PE2  + (int)pPars->MinSeqLen) > (int)pPars->PE2ReadLen)
			{
			pPars->NumPE2Underlen += 1;
			continue;
			}
		}

	bFailsPhred = false;
	if(bIsPE1Fastq && pPars->PE1QCSchema > 0 && pPars->MinPhredScore > 0)
		{
		pPars->szPE1QScoresBuff[pPars->PE1QSLen - pPars->Trim3 - ContamLen3PE1] = '\0';

		if(!MeetsMinPhredScoreThres(pPars->PE1QCSchema,pPars->MinPhredScore,(char *)&pPars->szPE1QScoresBuff[pPars->Trim5 + ContamLen5PE1]))
			{
			pPars->NumPE1UnderPhredScore += 1;
			bFailsPhred = true;
			}
		}

	if(bIsPE2Fastq && pPars->PE2QCSchema > 0 && pPars->MinPhredScore > 0)
		{
		pPars->szPE2QScoresBuff[pPars->PE2QSLen - pPars->Trim3 - ContamLen3PE2] = '\0';
		if(!MeetsMinPhredScoreThres(pPars->PE2QCSchema,pPars->MinPhredScore,(char *)&pPars->szPE2QScoresBuff[pPars->Trim5 + ContamLen5PE2]))
			{
			pPars->NumPE2UnderPhredScore += 1;
			bFailsPhred = true;
			}
		}
	if(bFailsPhred)
		continue;

	// trim 5' and 3' as needed
	pPE1Seq = (etSeqBase *)pPars->PE1RawReadsBuff;
	if(pPars->Trim5 || ContamLen5PE1)
		{
		pPE1Seq += (pPars->Trim5 + ContamLen5PE1);
		pPars->PE1ReadLen -= (pPars->Trim5 + ContamLen5PE1);
		}
	pPars->PE1ReadLen -= (pPars->Trim3 + ContamLen3PE1);
	if(pPars->TrimSeqLen > 0 && pPars->PE1ReadLen > (UINT32)pPars->TrimSeqLen)
		pPars->PE1ReadLen = pPars->TrimSeqLen;

	if(pPars->PE2ReadLen)
		{
		pPE2Seq = (etSeqBase *)pPars->PE2RawReadsBuff;
		if(pPars->Trim5 || ContamLen5PE2)
			{
			pPE2Seq += (pPars->Trim5 + ContamLen5PE2);
			pPars->PE2ReadLen -= (pPars->Trim5 + ContamLen5PE2);
			}
		pPars->PE2ReadLen -= (pPars->Trim3 + ContamLen3PE2);
		if(pPars->TrimSeqLen > 0 && pPars->PE2ReadLen > (UINT32)pPars->TrimSeqLen)
			pPars->PE2ReadLen = (UINT32)pPars->TrimSeqLen;
		}

	// check for excessive number of Ns 
	UINT32 Idx;
	UINT32 NumNs = 0;		// number of indeterminate bases in last 100bp window
	etSeqBase *pBase = pPE1Seq;
	for(Idx = 0; Idx < pPars->PE1ReadLen; Idx++,pBase++)
		{
		if(Idx >= 100 && (pBase[-100] & 0x07) == eBaseN)
				NumNs -= 1;
		if((*pBase & 0x07) == eBaseN)
			NumNs += 1;
		if(NumNs > (UINT32)pPars->MaxNs)
			break;
		}

	if(NumNs > (UINT32)pPars->MaxNs ||
			(pPars->PE1ReadLen <= 100 && NumNs > ((UINT32)(pPars->MaxNs * pPars->PE1ReadLen)/100)))
		{
		pPars->NumPE1ExcessNs += 1;
		continue;
		}
	if(pPars->PE2ReadLen)
		{
		NumNs = 0;		// number of indeterminate bases in last 100bp window
		pBase = pPE2Seq;
		for(Idx = 0; Idx < (UINT32)pPars->PE2ReadLen; Idx++,pBase++)
			{
			if(Idx >= 100 && (pBase[-100] & 0x07) == eBaseN)
					NumNs -= 1;
			if((*pBase & 0x07) == eBaseN)
				NumNs += 1;
			if(NumNs > (UINT32)pPars->MaxNs)
				break;
			}

		if(NumNs > (UINT32)pPars->MaxNs ||
				(pPars->PE2ReadLen <= 100 && NumNs > (UINT32)((pPars->MaxNs * pPars->PE2ReadLen)/100)))
			{
			pPars->NumPE2ExcessNs += 1;
			continue;
			}
		}
		
	if(pPars->SampleNth > 1 && (pPars->NumPE1ParsedReads % pPars->SampleNth))
		continue;

	if(pPars->MinAcceptedReadLen == 0 || pPars->MinAcceptedReadLen > (int)pPars->PE1ReadLen)
		pPars->MinAcceptedReadLen = pPars->PE1ReadLen;
	if(pPars->MaxAcceptedReadLen == 0 || pPars->MaxAcceptedReadLen <  (int)pPars->PE1ReadLen)
		pPars->MaxAcceptedReadLen = pPars->PE1ReadLen;

	AcquireSerialiseNxtProcRead();
	if(m_Sequences.MinSeqLen == 0 || m_Sequences.MinSeqLen > pPars->PE1ReadLen)
		m_Sequences.MinSeqLen = pPars->PE1ReadLen;
	if(m_Sequences.MaxSeqLen == 0 || m_Sequences.MaxSeqLen <  pPars->PE1ReadLen)
		m_Sequences.MaxSeqLen = pPars->PE1ReadLen;
	if((Rslt=AddSeq(pPars->PE1FileID,m_Sequences.bPESeqs ? cFlgSeqPE : 0,pPars->PE1ReadLen,pPE1Seq)) !=eBSFSuccess)
		{
		ReleaseSerialiseNxtProcRead();
		break;
		}
	pPars->NumPE1AcceptedReads += 1;
	pPars->AcceptedTotSeqLen += pPars->PE1ReadLen;

	if(pPars->PE2ReadLen)
		{
		if(pPars->MinAcceptedReadLen == 0 || pPars->MinAcceptedReadLen > (int)pPars->PE2ReadLen)
			pPars->MinAcceptedReadLen = pPars->PE1ReadLen;
		if(pPars->MaxAcceptedReadLen == 0 || pPars->MaxAcceptedReadLen <  (int)pPars->PE2ReadLen)
			pPars->MaxAcceptedReadLen = pPars->PE2ReadLen;

		if(m_Sequences.MinSeqLen == 0 || m_Sequences.MinSeqLen > pPars->PE2ReadLen)
			m_Sequences.MinSeqLen = pPars->PE2ReadLen;
		if(m_Sequences.MaxSeqLen == 0 || m_Sequences.MaxSeqLen <  pPars->PE2ReadLen)
			m_Sequences.MaxSeqLen = pPars->PE2ReadLen;
		CSeqTrans::ReverseComplement(pPars->PE2ReadLen,pPE2Seq);
		if((Rslt=AddSeq(pPars->PE2FileID,cFlgSeqPE2 | cFlgSeqPE,pPars->PE2ReadLen,pPE2Seq)) !=eBSFSuccess)
			{
			ReleaseSerialiseNxtProcRead();
			break;
			}
		pPars->NumPE2AcceptedReads += 1;
		pPars->AcceptedTotSeqLen += pPars->PE2ReadLen;
		}
	ReleaseSerialiseNxtProcRead();

	AcquireSerialiseReadsCtrl();
	pPars->pProcReadsCtrl->CurTotPE1ReadsAccepted += 1;
	ReleaseSerialiseReadsCtrl();

	if(pPars->Zreads > 0 && pPars->NumPE1AcceptedReads >= pPars->Zreads)
		{
		Rslt = eBSFSuccess;
		break;
		}
	}
pPars->Rslt = Rslt;
return(Rslt);
}


// LoadReadsThreaded
// Load reads from fasta, fastq or csfasta formated raw reads file
teBSFrsltCodes
CKangadna::LoadReadsThreaded(int MaxNs,				// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
					int MinPhredScore,		// filter out input sequences with mean Phred score lower than this threshold
					int Trim5,				// trim this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int MinSeqLen,		    // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 20..10000)
					int TrimSeqLen,			// trim sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...10000)
					int SampleNth,			// process every Nth reads
					int Zreads,				// maximum number of reads to accept for processing from any file
					char *pszPE1File,		// file containing reads (kangar or raw fasta/fastq)
					char *pszPE2File,		// if paired end processing then PE2 3' file containing reads
					int MaxNumThreads)		// use at most this any threads, 0 if no limit
{
teBSFrsltCodes Rslt;

tsProcReadsCtrl ProcReadsCtrl;              // reads ctrl shared over all worker threads
tsThreadFiltReadsPars *pFiltReadsPars;      // allocated to hold parameterisation and buffers for all worker threads
tsThreadFiltReadsPars *pCurThread;
int ThreadIdx;

bool bIsPE1Fastq;
int NumPE1DescrReads;
int NumPE1Underlen;
int NumPE1UnsupportedBases;
int NumPE1ExcessNs;
int NumPE1AcceptedReads;
int NumPE1ParsedReads;
int NumPE1UnderPhredScore;

int CurAcceptedMaxSeqLen;
int CurAcceptedMinSeqLen;
INT64 CurAcceptedTotSeqLen;
double CurAcceptedMeanSeqLen;

bool bIsPE2Fastq;
int NumPE2DescrReads;
int NumPE2Underlen;
int NumPE2UnsupportedBases;
int NumPE2ExcessNs;
int NumPE2AcceptedReads;
int NumPE2ParsedReads;
int NumPE2UnderPhredScore;

UINT32 NumContamVector;
UINT32 NumContamFlankTrim;


CFasta FastaPE1;
CFasta FastaPE2;

tsReadFile *pPE1ReadFile;
tsReadFile *pPE2ReadFile;

UINT32 PE1FileID;
UINT32 PE2FileID;

int MaxAllowedFiles;

if(MaxNumThreads == 0 || MaxNumThreads > min(m_NumThreads, 15))		// limiting to 15 threads as a) ~600MB required per thread and b) threads are serialised when loading from file so threads likely to be frequently blocked
	MaxNumThreads = min(m_NumThreads, 15);

if(MaxNumThreads == 1)
	return(LoadReads(MaxNs,MinPhredScore,Trim5,Trim3,MinSeqLen,TrimSeqLen,SampleNth,Zreads,pszPE1File,pszPE2File));

if((pFiltReadsPars = new tsThreadFiltReadsPars[MaxNumThreads])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for %d worker threads...",MaxNumThreads);
	Reset(false);
	return(eBSFerrMem);
	}

memset(&ProcReadsCtrl,0,sizeof(ProcReadsCtrl));

m_Sequences.bPESeqs = (pszPE2File == NULL || pszPE2File[0] == '\0') ? false : true;

// error if memory not preallocated to hold sequences
if(m_Sequences.pSeqs2Assemb == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected memory to have been preallocated for holding compacted sequences");
	delete pFiltReadsPars;
	return(eBSFerrInternal);
	}

PE1FileID = m_NumRawFiles + 1;

pPE1ReadFile = &m_SrcFiles[m_NumRawFiles];
if(m_Sequences.bPESeqs)
	{
	pPE2ReadFile = &m_SrcFiles[m_NumRawFiles+1];
	PE2FileID = m_NumRawFiles + 2;
	}
m_NumRawFiles += m_Sequences.bPESeqs ? 2 : 1;

MaxAllowedFiles = (m_Sequences.bPESeqs ? 2 : 1) * cKDNAMaxInFileSpecs;
if(m_NumRawFiles > MaxAllowedFiles)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Too many files to process - max files: %d files...",MaxAllowedFiles);
	delete pFiltReadsPars;
	Reset(false);
	return(eBSFerrParams);
	}

memset(pPE1ReadFile,0,sizeof(tsReadFile));
pPE1ReadFile->FileID = PE1FileID;
strcpy((char *)pPE1ReadFile->szFile,pszPE1File);
if(m_Sequences.bPESeqs)
	{
	pPE1ReadFile->PEFileID = PE2FileID;
	memset(pPE2ReadFile,0,sizeof(tsReadFile));
	pPE2ReadFile->FileID = PE2FileID;
	pPE2ReadFile->PEFileID = PE1FileID;
	strcpy((char *)pPE2ReadFile->szFile,pszPE2File);
	}

int PE1QCSchema;
int PE2QCSchema;

if((Rslt=(teBSFrsltCodes)FastaPE1.FastaEstSizes(pszPE1File,NULL,NULL,NULL,NULL,NULL,&PE1QCSchema))<=1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",pszPE1File,FastaPE1.ErrText((teBSFrsltCodes)Rslt),FastaPE1.GetErrMsg());
	delete pFiltReadsPars;
	Reset(false);
	return(Rslt);
	}

if((Rslt=(teBSFrsltCodes)FastaPE1.Open(pszPE1File,true,cMaxStageBuffSize))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",pszPE1File,FastaPE1.ErrText((teBSFrsltCodes)Rslt),FastaPE1.GetErrMsg());
	delete pFiltReadsPars;
	Reset(false);
	return(Rslt);
	}

if(FastaPE1.IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to load '%s', SOLiD colorspace not supported",pszPE1File);
	delete pFiltReadsPars;
	FastaPE1.Close();
	Reset(false);
	return(eBSFerrFileType);
	}
bIsPE1Fastq = FastaPE1.IsFastq();

if(m_Sequences.bPESeqs)
	{
	if((Rslt=(teBSFrsltCodes)FastaPE2.FastaEstSizes(pszPE2File,NULL,NULL,NULL,NULL,NULL,&PE2QCSchema))<=1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",pszPE2File,FastaPE2.ErrText((teBSFrsltCodes)Rslt),FastaPE2.GetErrMsg());
		delete pFiltReadsPars;
		FastaPE1.Close();
		Reset(false);
		return(Rslt);
		}

	if((Rslt=(teBSFrsltCodes)FastaPE2.Open(pszPE2File,true,cMaxStageBuffSize))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",pszPE2File,FastaPE2.ErrText((teBSFrsltCodes)Rslt),FastaPE2.GetErrMsg());
		delete pFiltReadsPars;
		FastaPE1.Close();
		Reset(false);
		return(Rslt);
		}

	if(FastaPE2.IsSOLiD())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to load '%s', SOLiD colorspace not supported",pszPE2File);
		delete pFiltReadsPars;
		FastaPE1.Close();
		FastaPE2.Close();
		Reset(false);
		return(eBSFerrFileType);
		}
	bIsPE2Fastq = FastaPE2.IsFastq();
	if(bIsPE1Fastq != 	bIsPE2Fastq)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Paired end file formats not of same type: '%s' is %s, '%s' is %s", 
										bIsPE1Fastq ? "Fastq" : "Fasta",pszPE2File,bIsPE2Fastq ? "Fastq" : "Fasta", pszPE2File);
		delete pFiltReadsPars;
		FastaPE1.Close();
		FastaPE2.Close();
		Reset(false);
		return(eBSFerrFileType);		
		}
	}
else
	bIsPE2Fastq = false;

UINT64 CurWorkSetSize;
CurWorkSetSize = m_Sequences.AllocMemSeqs2Assemb + m_Sequences.AllocMemSeqStarts + m_Sequences.AllocMemSfx + m_Sequences.AllocMemSeqFlags + (sizeof(tsThreadFiltReadsPars) * MaxNumThreads);

if(CurWorkSetSize != m_CurMaxMemWorkSetBytes)
	{
	if(!SetMaxMemWorkSetSize((size_t)CurWorkSetSize))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to set maximum working set size (bytes) to %lld bytes",CurWorkSetSize);
		FastaPE1.Close();
		if(m_Sequences.bPESeqs)
			FastaPE2.Close();
		Reset(false);
		return(eBSFerrMaxDirEls);
		}
	}

NumPE1DescrReads = 0;
NumPE1Underlen = 0;
NumPE1UnderPhredScore = 0;
NumPE1UnsupportedBases = 0;
NumPE1ExcessNs = 0;
NumPE1AcceptedReads = 0;
NumPE1ParsedReads = 0;

NumPE2DescrReads = 0;
NumPE2Underlen = 0;
NumPE2UnderPhredScore = 0;
NumPE2UnsupportedBases = 0;
NumPE2ExcessNs = 0;
NumPE2AcceptedReads = 0;
NumPE2ParsedReads = 0;

CurAcceptedMaxSeqLen = 0;
CurAcceptedMinSeqLen = 0;
CurAcceptedTotSeqLen = 0;
CurAcceptedMeanSeqLen = 0.0;

NumContamVector = 0;
NumContamFlankTrim = 0;

// initialise filter parameters and start up worker threads

#ifndef _WIN32
// increase the default stack of just 2MB
size_t defaultStackSize;
pthread_attr_t threadattr; 
pthread_attr_init(&threadattr);
pthread_attr_getstacksize(&threadattr, &defaultStackSize);
if(defaultStackSize != cReadsThreadStackSize)
	pthread_attr_setstacksize(&threadattr, cReadsThreadStackSize);
#endif

ProcReadsCtrl.bProcPE = m_Sequences.bPESeqs == 0 ? false : true;
ProcReadsCtrl.bIsPE1Fastq = bIsPE1Fastq;
ProcReadsCtrl.bIsPE2Fastq = bIsPE2Fastq;
ProcReadsCtrl.bPE1PE2AllReadsProc = false;
ProcReadsCtrl.MinPhredScore = MinPhredScore;
ProcReadsCtrl.NumPE1DescrReads = 0;
ProcReadsCtrl.NumPE1ParsedReads = 0;
ProcReadsCtrl.NumPE2DescrReads = 0;
ProcReadsCtrl.NumPE2ParsedReads = 0;
ProcReadsCtrl.CurTotPE1ReadsParsed = 0;
ProcReadsCtrl.CurTotPE1ReadsAccepted = 0;
ProcReadsCtrl.pFastaPE1 = &FastaPE1;
ProcReadsCtrl.pFastaPE2 = &FastaPE2;
ProcReadsCtrl.PE1FileID = PE1FileID;
ProcReadsCtrl.PE2FileID = PE2FileID;
ProcReadsCtrl.Rslt = eBSFSuccess;
strcpy(ProcReadsCtrl.szPE1File,pszPE1File);
if(bIsPE2Fastq)
	strcpy(ProcReadsCtrl.szPE2File,pszPE2File);
else
	ProcReadsCtrl.szPE2File[0] = '\0';

pCurThread = pFiltReadsPars;
for(ThreadIdx = 0; ThreadIdx < MaxNumThreads; ThreadIdx++,pCurThread++)
	{
	memset(pCurThread,0,sizeof(tsThreadFiltReadsPars));
	pCurThread->ThreadIdx = ThreadIdx + 1;
	pCurThread->pThis = this;
	pCurThread->MaxNs = MaxNs;
	pCurThread->MinPhredScore = MinPhredScore;
	pCurThread->PE1QCSchema = PE1QCSchema;
	pCurThread->PE2QCSchema = PE2QCSchema;
	pCurThread->Trim5 = Trim5;
	pCurThread->Trim3 = Trim3;
	pCurThread->MinSeqLen = MinSeqLen;
	pCurThread->TrimSeqLen = TrimSeqLen;
	pCurThread->SampleNth = SampleNth;
	pCurThread->Zreads = Zreads;
	pCurThread->PE1FileID = PE1FileID;
	pCurThread->PE2FileID = PE2FileID;
	pCurThread->pProcReadsCtrl = &ProcReadsCtrl;

#ifdef _WIN32
	pCurThread->threadHandle = (HANDLE)_beginthreadex(NULL,cReadsThreadStackSize,ThreadedFiltReads,pCurThread,0,&pCurThread->threadID);
#else
	pCurThread->threadRslt = pthread_create (&pCurThread->threadID , &threadattr , ThreadedFiltReads , pCurThread);
#endif
	}
#ifndef _WIN32
pthread_attr_destroy(&threadattr);		// no longer required
#endif

// allow threads a few seconds to startup
#ifdef _WIN32
	Sleep(3000);
#else
	sleep(3);
#endif
UINT32 CurNumParsed;
UINT32 CurNumAccepted;

pCurThread = pFiltReadsPars;
for(ThreadIdx = 0; ThreadIdx < MaxNumThreads; ThreadIdx++,pCurThread++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject(pCurThread->threadHandle, 60000))
		{
		AcquireSerialiseReadsCtrl();
		CurNumParsed =  ProcReadsCtrl.CurTotPE1ReadsParsed;
		CurNumAccepted = ProcReadsCtrl.CurTotPE1ReadsAccepted;
		ReleaseSerialiseReadsCtrl(); 
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u %s sequences parsed, %d accepted",CurNumParsed,m_Sequences.bPESeqs ? "PE" : "SE",CurNumAccepted);
		}
	CloseHandle( pCurThread->threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while((JoinRlt = pthread_timedjoin_np(pCurThread->threadID, NULL, &ts)) != 0)
		{
		AcquireSerialiseReadsCtrl();
		CurNumParsed =  ProcReadsCtrl.CurTotPE1ReadsParsed;
		CurNumAccepted = ProcReadsCtrl.CurTotPE1ReadsAccepted;
		ReleaseSerialiseReadsCtrl(); 
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: %u %s sequences parsed, %d accepted",CurNumParsed,m_Sequences.bPESeqs ? "PE" : "SE",CurNumAccepted);
		ts.tv_sec += 60;
		}
#endif
	NumPE1ParsedReads += pCurThread->NumPE1ParsedReads;
	NumPE1AcceptedReads += pCurThread->NumPE1AcceptedReads;

	if(CurAcceptedMinSeqLen == 0 || CurAcceptedMinSeqLen > (int)pCurThread->MinAcceptedReadLen)
		CurAcceptedMinSeqLen = pCurThread->MinAcceptedReadLen;
	if(CurAcceptedMaxSeqLen == 0 || CurAcceptedMaxSeqLen <  (int)pCurThread->MaxAcceptedReadLen)
		CurAcceptedMaxSeqLen = pCurThread->MaxAcceptedReadLen;

	CurAcceptedTotSeqLen += pCurThread->AcceptedTotSeqLen;

	NumPE1Underlen += pCurThread->NumPE1Underlen;
	NumPE2Underlen += pCurThread->NumPE2Underlen; 

	NumPE1ExcessNs += pCurThread->NumPE1ExcessNs;
	NumPE2ExcessNs += pCurThread->NumPE2ExcessNs;

	NumPE1UnderPhredScore += pCurThread->NumPE1UnderPhredScore;
	NumPE2UnderPhredScore += pCurThread->NumPE2UnderPhredScore;

	NumPE2ParsedReads += pCurThread->NumPE2ParsedReads;
	NumPE2AcceptedReads += pCurThread->NumPE2AcceptedReads;

	NumContamVector += pCurThread->NumContamVector;
	NumContamFlankTrim += pCurThread->NumContamFlankTrim;
	}

// threads have all terminated
delete pFiltReadsPars;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Parsed %d, Accepted %d %s",NumPE1ParsedReads,NumPE1AcceptedReads,m_Sequences.bPESeqs ? "paired sequences":"sequences");

if(Rslt != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors processing file: %s ",pszPE1File);
	while(FastaPE1.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,FastaPE1.GetErrMsg());
	FastaPE1.Close();
	if(m_Sequences.bPESeqs)
		FastaPE2.Close();
	Reset(false);
	return(Rslt);
	}
FastaPE1.Close();
if(m_Sequences.bPESeqs)
	FastaPE2.Close();

pPE1ReadFile->SeqsParsed = NumPE1ParsedReads;
pPE1ReadFile->SeqsUnderLen = NumPE1Underlen;
pPE1ReadFile->SeqsExcessNs = NumPE1ExcessNs;
pPE1ReadFile->SeqsAccepted = NumPE1AcceptedReads;

if(m_Sequences.bPESeqs)
	{
	pPE2ReadFile->SeqsParsed = NumPE2ParsedReads;
	pPE2ReadFile->SeqsUnderLen = NumPE2Underlen;
	pPE2ReadFile->SeqsExcessNs = NumPE2ExcessNs;
	pPE2ReadFile->SeqsAccepted = NumPE2AcceptedReads;
	}

m_Sequences.TotSeqsParsed += NumPE1ParsedReads + NumPE2ParsedReads;
m_Sequences.TotSeqsUnderLen += NumPE1Underlen + NumPE2Underlen;
m_Sequences.TotSeqsExcessNs += NumPE1ExcessNs + NumPE2ExcessNs; 
m_Sequences.TotSeqs2Assemb += NumPE1AcceptedReads + NumPE2AcceptedReads;
m_Sequences.NumPE1Seqs2Assemb += NumPE1AcceptedReads;
m_Sequences.NumPE2Seqs2Assemb += NumPE2AcceptedReads;
if(m_Sequences.TotSeqs2Assemb > 0)
	m_Sequences.MeanSeqLen = (double)m_Sequences.Seqs2AssembLen / m_Sequences.TotSeqs2Assemb;
else 
	m_Sequences.MeanSeqLen = 0.0;
m_Sequences.bSeqs2AssembDirty = true;	

if(NumPE1AcceptedReads + NumPE2AcceptedReads)
	CurAcceptedMeanSeqLen = (double)CurAcceptedTotSeqLen / ((double)NumPE1AcceptedReads + NumPE2AcceptedReads);
else
	CurAcceptedMeanSeqLen = 0.0;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Completed parsing %d, accepted %d %s, accepted sequences min length: %d max length: %d mean length: %1.1f",
								NumPE1ParsedReads,NumPE1AcceptedReads,m_Sequences.bPESeqs ? "paired sequences":"sequences",CurAcceptedMinSeqLen,CurAcceptedMaxSeqLen,CurAcceptedMeanSeqLen);

if(m_pContaminants != NULL)
	{
	if(m_Sequences.bPESeqs)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Removed %u paired sequences because at least one end was totally contained in a vector contaminate sequence", NumContamVector);
	else
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Removed %u sequences because these were totally contained in a vector contaminate sequence", NumContamVector);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Trimmed %u sequences because these flank overlapped onto a contaminate sequence", NumContamFlankTrim);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: total of %d sequences were under length %d and %d had excessive indeterminate (%d Ns)",NumPE1Underlen + NumPE2Underlen,MinSeqLen,NumPE1ExcessNs+NumPE2ExcessNs,MaxNs);

if((PE1QCSchema > 0 || PE2QCSchema > 0) && MinPhredScore > 0)
	{
	if(!m_Sequences.bPESeqs)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: total of %d sequences had Phred scores under %d",NumPE1UnderPhredScore,MinPhredScore);
	else
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: total of %d (%d PE1, %d PE2) sequences had Phred scores under %d",NumPE1UnderPhredScore + NumPE2UnderPhredScore,NumPE1UnderPhredScore,NumPE2UnderPhredScore,MinPhredScore);
	}
else
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: No Phred score filtering");

return((teBSFrsltCodes)(NumPE1AcceptedReads + NumPE2AcceptedReads));
}


// LoadReads
// Load reads from fasta, fastq or csfasta formated raw reads file
teBSFrsltCodes
CKangadna::LoadReads(int MaxNs,				// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
					int MinPhredScore,		// filter out input sequences with mean Phred score lower than this threshold
					int Trim5,				// trim this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int MinSeqLen,		    // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 20..10000)
					int TrimSeqLen,			// trim sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...10000)
					int SampleNth,			// process every Nth reads
					int Zreads,				// maximum number of reads to accept for processing from any file
					char *pszPE1File,		// file containing reads (kangar or raw fasta/fastq)
					char *pszPE2File)		// if paired end processing then PE2 3' file containing reads
{
teBSFrsltCodes Rslt;

bool bIsPE1Fastq;
int NumPE1DescrReads;
int NumPE1Underlen;
int NumPE1UnsupportedBases;
int NumPE1ExcessNs;
int NumPE1AcceptedReads;
int NumPE1ParsedReads;
int CurAcceptedMaxSeqLen;
int CurAcceptedMinSeqLen;
INT64 CurAcceptedTotSeqLen;
double CurAcceptedMeanSeqLen;
UINT32 PE1ReadLen;
int PE1DescrLen;
UINT8 szPE1DescrBuff[4096];
UINT32 PE1QSLen;
UINT8 *pszPE1QScoresBuff;
int NumPE1UnderPhredScore;

bool bIsPE2Fastq;
int NumPE2DescrReads;
int NumPE2Underlen;
int NumPE2UnsupportedBases;
int NumPE2ExcessNs;
int NumPE2AcceptedReads;
int NumPE2ParsedReads;
UINT32 PE2ReadLen;
int PE2DescrLen;
UINT8 szPE2DescrBuff[1024];

UINT32 PE2QSLen;
UINT8 *pszPE2QScoresBuff;
int NumPE2UnderPhredScore;

int ContamLen5PE1;
int ContamLen3PE1;
int ContamLen5PE2;
int ContamLen3PE2;

UINT32 NumContamVector;
UINT32 NumContamFlankTrim;


etSeqBase *pPE1Seq;
etSeqBase *pPE2Seq;
UINT8 *pPE1RawReadsBuff;
UINT8 *pPE2RawReadsBuff;
CFasta FastaPE1;
CFasta FastaPE2;

tsReadFile *pPE1ReadFile;
tsReadFile *pPE2ReadFile;

UINT32 PE1FileID;
UINT32 PE2FileID;

int MaxAllowedFiles;

if(m_NumThreads > 1)
	return(LoadReadsThreaded(MaxNs,MinPhredScore,Trim5,Trim3,MinSeqLen,TrimSeqLen,SampleNth,Zreads,pszPE1File,pszPE2File,0));

m_Sequences.bPESeqs = (pszPE2File == NULL || pszPE2File[0] == '\0') ? false : true;

// error if memory not preallocated to hold sequences
if(m_Sequences.pSeqs2Assemb == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected memory to have been preallocated for holding compacted sequences");
	return(eBSFerrInternal);
	}

PE1FileID = m_NumRawFiles + 1;

pPE1ReadFile = &m_SrcFiles[m_NumRawFiles];
if(m_Sequences.bPESeqs)
	{
	pPE2ReadFile = &m_SrcFiles[m_NumRawFiles+1];
	PE2FileID = m_NumRawFiles + 2;
	}
m_NumRawFiles += m_Sequences.bPESeqs ? 2 : 1;

MaxAllowedFiles = (m_Sequences.bPESeqs ? 2 : 1) * cKDNAMaxInFileSpecs;
if(m_NumRawFiles > MaxAllowedFiles)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Too many files to process - max files: %d files...",MaxAllowedFiles);
	Reset(false);
	return(eBSFerrParams);
	}

memset(pPE1ReadFile,0,sizeof(tsReadFile));
pPE1ReadFile->FileID = PE1FileID;
strcpy((char *)pPE1ReadFile->szFile,pszPE1File);
if(m_Sequences.bPESeqs)
	{
	pPE1ReadFile->PEFileID = PE2FileID;
	memset(pPE2ReadFile,0,sizeof(tsReadFile));
	pPE2ReadFile->FileID = PE2FileID;
	pPE2ReadFile->PEFileID = PE1FileID;
	strcpy((char *)pPE2ReadFile->szFile,pszPE2File);
	}

if((pPE1RawReadsBuff = new UINT8 [cMaxRawSeqLen+16]) == NULL)					// 16 as a small safety margin!
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for PE1 raw reads buffering...");
	Reset(false);
	return(eBSFerrMem);
	}
if((pszPE1QScoresBuff = new UINT8 [cMaxQScoreLen+16]) == NULL)					// 16 as a small safety margin!
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for PE1 raw reads quality score buffering...");
	delete pPE1RawReadsBuff;
	Reset(false);
	return(eBSFerrMem);
	}

if(m_Sequences.bPESeqs)
	{
	if((pPE2RawReadsBuff = new UINT8 [cMaxRawSeqLen+16]) == NULL)				// 16 as a small safety margin!
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for PE2 raw reads buffering...");
		delete pPE1RawReadsBuff;
		Reset(false);
		return(eBSFerrMem);
		}
	if((pszPE2QScoresBuff = new UINT8 [cMaxQScoreLen+16]) == NULL)				// 16 as a small safety margin!
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for PE2 raw reads quality score buffering...");
		delete pPE1RawReadsBuff;
		delete pszPE1QScoresBuff;
		delete pPE2RawReadsBuff;
		Reset(false);
		return(eBSFerrMem);
		}
	}
else
	{
	pPE2RawReadsBuff = NULL;
	pszPE2QScoresBuff = NULL;
	}

int PE1QCSchema;
int PE2QCSchema;

if((Rslt=(teBSFrsltCodes)FastaPE1.FastaEstSizes(pszPE1File,NULL,NULL,NULL,NULL,NULL,&PE1QCSchema))<=1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",pszPE1File,FastaPE1.ErrText((teBSFrsltCodes)Rslt),FastaPE1.GetErrMsg());
	delete pPE1RawReadsBuff;
	delete pszPE1QScoresBuff;
	if(pPE2RawReadsBuff)
		delete pPE2RawReadsBuff;
	if(pszPE2QScoresBuff)
		delete pszPE2QScoresBuff;
	Reset(false);
	return(Rslt);
	}

if((Rslt=(teBSFrsltCodes)FastaPE1.Open(pszPE1File,true,cMaxStageBuffSize))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",pszPE1File,FastaPE1.ErrText((teBSFrsltCodes)Rslt),FastaPE1.GetErrMsg());
	delete pPE1RawReadsBuff;
	delete pszPE1QScoresBuff;
	if(pPE2RawReadsBuff)
		delete pPE2RawReadsBuff;
	if(pszPE2QScoresBuff)
		delete pszPE2QScoresBuff;	
	Reset(false);
	return(Rslt);
	}

if(FastaPE1.IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to load '%s', SOLiD colorspace not supported",pszPE1File);
	delete pPE1RawReadsBuff;
	delete pszPE1QScoresBuff;
	if(pPE2RawReadsBuff)
		delete pPE2RawReadsBuff;
	if(pszPE2QScoresBuff)
		delete pszPE2QScoresBuff;
	Reset(false);
	return(eBSFerrFileType);
	}
bIsPE1Fastq = FastaPE1.IsFastq();

if(m_Sequences.bPESeqs)
	{
	if((Rslt=(teBSFrsltCodes)FastaPE2.FastaEstSizes(pszPE2File,NULL,NULL,NULL,NULL,NULL,&PE2QCSchema))<=1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",pszPE2File,FastaPE2.ErrText((teBSFrsltCodes)Rslt),FastaPE2.GetErrMsg());
		delete pPE1RawReadsBuff;
		delete pszPE1QScoresBuff;
		delete pPE2RawReadsBuff;
		delete pszPE2QScoresBuff;
		Reset(false);
		return(Rslt);
		}

	if((Rslt=(teBSFrsltCodes)FastaPE2.Open(pszPE2File,true,cMaxStageBuffSize))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",pszPE2File,FastaPE2.ErrText((teBSFrsltCodes)Rslt),FastaPE2.GetErrMsg());
		FastaPE1.Close();
		delete pPE1RawReadsBuff;
		delete pszPE1QScoresBuff;
		delete pPE2RawReadsBuff;
		delete pszPE2QScoresBuff;
		Reset(false);
		return(Rslt);
		}

	if(FastaPE2.IsSOLiD())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to load '%s', SOLiD colorspace not supported",pszPE2File);
		FastaPE1.Close();
		FastaPE2.Close();
		delete pPE1RawReadsBuff;
		delete pszPE1QScoresBuff;
		delete pPE2RawReadsBuff;
		delete pszPE2QScoresBuff;
		Reset(false);
		return(eBSFerrFileType);
		}
	bIsPE2Fastq = FastaPE2.IsFastq();
	if(bIsPE1Fastq != 	bIsPE2Fastq)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Paired end file formats not of same type: '%s' is %s, '%s' is %s", 
										bIsPE1Fastq ? "Fastq" : "Fasta",pszPE2File,bIsPE2Fastq ? "Fastq" : "Fasta", pszPE2File);
		FastaPE1.Close();
		FastaPE2.Close();
		delete pPE1RawReadsBuff;
		delete pszPE1QScoresBuff;
		delete pPE2RawReadsBuff;
		delete pszPE2QScoresBuff;
		Reset(false);
		return(eBSFerrFileType);		
		}
	}

UINT64 CurWorkSetSize;
CurWorkSetSize = m_Sequences.AllocMemSeqs2Assemb + m_Sequences.AllocMemSeqStarts + m_Sequences.AllocMemSfx + m_Sequences.AllocMemSeqFlags;

if(CurWorkSetSize != m_CurMaxMemWorkSetBytes)
	{
	if(!SetMaxMemWorkSetSize((size_t)CurWorkSetSize))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to set maximum working set size (bytes) to %lld bytes",CurWorkSetSize);
		FastaPE1.Close();
		delete pPE1RawReadsBuff;
		delete pszPE1QScoresBuff;
		if(pPE2RawReadsBuff)
			{
			FastaPE2.Close();
			delete pPE2RawReadsBuff;
			delete pszPE2QScoresBuff;
			}
		Reset(false);
		return(eBSFerrMaxDirEls);
		}
	}

NumPE1DescrReads = 0;
NumPE1Underlen = 0;
NumPE1UnderPhredScore = 0;
NumPE1UnsupportedBases = 0;
NumPE1ExcessNs = 0;
NumPE1AcceptedReads = 0;
NumPE1ParsedReads = 0;
PE1ReadLen = 0;
PE1DescrLen = 0;


NumPE2DescrReads = 0;
NumPE2Underlen = 0;
NumPE2UnderPhredScore = 0;
NumPE2UnsupportedBases = 0;
NumPE2ExcessNs = 0;
NumPE2AcceptedReads = 0;
NumPE2ParsedReads = 0;
PE2ReadLen = 0;
PE2DescrLen = 0;
CurAcceptedMaxSeqLen = 0;
CurAcceptedMinSeqLen = 0;
CurAcceptedTotSeqLen = 0;
CurAcceptedMeanSeqLen = 0.0;
NumContamVector = 0;
NumContamFlankTrim = 0;

if(MaxNs)
	srand((unsigned)time( NULL ));
gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Parsed 0, Accepted 0 %s ...",m_Sequences.bPESeqs ? "paired sequences":"sequences");
time_t Started = time(0);
while((Rslt = (teBSFrsltCodes)(PE1ReadLen = FastaPE1.ReadSequence(pPE1RawReadsBuff,cMaxRawSeqLen,true,false))) > eBSFSuccess)
	{
	NumPE1ParsedReads += 1;
	if(!(NumPE1ParsedReads % 1000))
		{
		time_t Now = time(0);
		unsigned long ElapsedSecs = (unsigned long) (Now - Started);
		if(ElapsedSecs >= 60)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Parsed %d, Accepted %d %s",NumPE1ParsedReads,NumPE1AcceptedReads,m_Sequences.bPESeqs ? "paired sequences":"sequences");
			Started = Now;
			}
		}

	if(PE1ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		NumPE1DescrReads += 1;
		PE1DescrLen = FastaPE1.ReadDescriptor((char *)szPE1DescrBuff,sizeof(szPE1DescrBuff)-1);
		szPE1DescrBuff[sizeof(szPE1DescrBuff)-1] = '\0';

		if(bIsPE1Fastq && MinPhredScore > 0)
			{
			PE1QSLen = FastaPE1.ReadQValues((char *)pszPE1QScoresBuff, cMaxQScoreLen);
			pszPE1QScoresBuff[PE1QSLen] = '\0';
			}
		else
			{
			PE1QSLen = 0;
			pszPE1QScoresBuff[0] = '\0';
			}


		PE1ReadLen = FastaPE1.ReadSequence(pPE1RawReadsBuff,cMaxRawSeqLen);
		if(PE1ReadLen < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence after %d sequences parsed from file '%s'",NumPE1DescrReads,pszPE1File);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE1DescrBuff);
			delete pPE1RawReadsBuff;
			delete pszPE1QScoresBuff;
			FastaPE1.Close();
			if(m_Sequences.bPESeqs)
				{
				delete pPE2RawReadsBuff;
				delete pszPE2QScoresBuff;
				FastaPE2.Close();
				}
			Reset(false);
			return(eBSFerrParse);
			}

		if (PE1QSLen > 0 && PE1ReadLen != PE1QSLen)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Problem parsing sequence after %d sequences parsed from file '%s'", NumPE1DescrReads, pszPE1File);
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected sequence length to be same as quality length. Last descriptor parsed: %s", szPE1DescrBuff);
			delete pPE1RawReadsBuff;
			delete pszPE1QScoresBuff;
			FastaPE1.Close();
			if(m_Sequences.bPESeqs)
				{
				delete pPE2RawReadsBuff;
				delete pszPE2QScoresBuff;
				FastaPE2.Close();
				}
			Reset(false);
			return(eBSFerrParse);
			}

		if(m_Sequences.bPESeqs)
			{
			Rslt = (teBSFrsltCodes)(PE2ReadLen = FastaPE2.ReadSequence(pPE2RawReadsBuff,cMaxRawSeqLen,true,false));
			if(Rslt <= eBSFSuccess)
				{
				if(Rslt < eBSFSuccess)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors processing file: %s ",pszPE2File);
					while(FastaPE2.NumErrMsgs())
						gDiagnostics.DiagOut(eDLFatal,gszProcName,FastaPE2.GetErrMsg());
					}
				else
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Insufficent reads in file: %s ",pszPE2File);
				delete pPE1RawReadsBuff;
				delete pszPE1QScoresBuff;
				FastaPE1.Close();
				delete pPE2RawReadsBuff;
				delete pszPE2QScoresBuff;
				FastaPE2.Close();
				Reset(false);
				return(Rslt);
				}
			NumPE2ParsedReads += 1;
			if(PE2ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
				{
				NumPE2DescrReads += 1;
				PE2DescrLen = FastaPE2.ReadDescriptor((char *)szPE2DescrBuff,sizeof(szPE2DescrBuff)-1);
				szPE2DescrBuff[sizeof(szPE2DescrBuff)-1] = '\0'; 

				if(bIsPE2Fastq && MinPhredScore > 0)
					{
					PE2QSLen = FastaPE2.ReadQValues((char *)pszPE2QScoresBuff, cMaxQScoreLen);
					pszPE2QScoresBuff[PE2QSLen] = '\0';
					}
				else
					{
					PE2QSLen = 0;
					pszPE2QScoresBuff[0] = '\0';
					}
				PE2ReadLen = FastaPE2.ReadSequence(pPE2RawReadsBuff,cMaxRawSeqLen);
				if(PE2ReadLen < 1)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence after %d sequences parsed from file '%s'",NumPE2DescrReads,pszPE2File);
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE2DescrBuff);
					delete pPE1RawReadsBuff;
					delete pszPE1QScoresBuff;
					FastaPE1.Close();
					delete pPE2RawReadsBuff;
					delete pszPE2QScoresBuff;
					FastaPE2.Close();
					Reset(false);
					return(eBSFerrParse);
					}

				if (PE2QSLen > 0 && PE2ReadLen != PE2QSLen)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Problem parsing sequence after %d sequences parsed from file '%s'", NumPE2DescrReads, pszPE2File);
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected sequence length to be same as quality length. Last descriptor parsed: %s", szPE2DescrBuff);
					delete pPE1RawReadsBuff;
					delete pszPE1QScoresBuff;
					FastaPE1.Close();
					delete pPE2RawReadsBuff;
					delete pszPE2QScoresBuff;
					FastaPE2.Close();
					Reset(false);
					return(eBSFerrParse);
					}
				}
			else
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"raw sequence file '%s' processing error: %s ",pszPE2File, 
							Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
				delete pPE1RawReadsBuff;
				delete pszPE1QScoresBuff;
				FastaPE1.Close();
				delete pPE2RawReadsBuff;
				delete pszPE2QScoresBuff;
				FastaPE2.Close();
				Reset(false);
				return(eBSFerrParse);
				}
			}

			// if trimming any contaminant sequences then check for contaminant overlaps onto this read(s)
		if(m_pContaminants != NULL)
			{
			bool bContamVector = false;

			// currently treating any contaminant match errors as if simply there was no overlap - should really report these!!!
			if((ContamLen5PE1 = m_pContaminants->MatchContaminants(eAOF5PE1Targ,1,Trim5+1,PE1ReadLen,pPE1RawReadsBuff)) <= Trim5)
				ContamLen5PE1 = 0;
			else
				{
				if(ContamLen5PE1 == PE1ReadLen)	// will be equal if totally contained
					bContamVector = true;
				else
					if(ContamLen5PE1 > Trim5)
						ContamLen5PE1 -= Trim5;
				}
			if(!bContamVector)
				{
				if((ContamLen3PE1 = m_pContaminants->MatchContaminants(eAOF3PE1Targ,1,Trim3+1,PE1ReadLen,pPE1RawReadsBuff)) <= Trim3)
					ContamLen3PE1 = 0;
				else
					{
					if(ContamLen3PE1 == PE1ReadLen)	// will be equal if totally contained
						bContamVector = true;
					else
						if(ContamLen3PE1 > Trim3)
							ContamLen3PE1 -= Trim3;
					}
				}

			if(m_Sequences.bPESeqs && !bContamVector)
				{
				if((ContamLen5PE2 = m_pContaminants->MatchContaminants(eAOF5PE2Targ,1,Trim5+1,PE2ReadLen,pPE2RawReadsBuff)) <= Trim5)
					ContamLen5PE2 = 0;
				else
					{
					if(ContamLen5PE2 == PE2ReadLen)	// will be equal if totally contained
						bContamVector = true;
					else
						if(ContamLen5PE2 > Trim5)
							ContamLen5PE2 -= Trim5;
					}
				if(!bContamVector)
					{
					if((ContamLen3PE2 = m_pContaminants->MatchContaminants(eAOF3PE2Targ,1,Trim3+1,PE2ReadLen,pPE2RawReadsBuff)) <= Trim3)
						ContamLen3PE2 = 0;
					else
						{
						if(ContamLen3PE2 == PE2ReadLen)	// will be equal if totally contained
							bContamVector = true;
						else
						if(ContamLen3PE2 > Trim3)
							ContamLen3PE2 -= Trim3;
						}
					}
				}
			else
				{
				ContamLen5PE2 = 0;
				ContamLen3PE2 = 0;
				}
			if(bContamVector)
				{
				NumContamVector += 1;
				continue;
				}
			if(ContamLen5PE1 || ContamLen5PE2)
				NumContamFlankTrim += 1;		
			if(ContamLen3PE1 || ContamLen3PE2)
				NumContamFlankTrim += 1;
			}
		else
			{
			ContamLen5PE1 = 0;
			ContamLen3PE1 = 0;
			ContamLen5PE2 = 0;
			ContamLen3PE2 = 0;
			}

		// ensure would still have a sequence of at least cMinSeqLen after any end trims were applied
		if((Trim5 + Trim3 + ContamLen5PE1 + ContamLen3PE1 + (int)MinSeqLen) > (int)PE1ReadLen)
			{
			NumPE1Underlen += 1;
			continue;
			}

		if(m_Sequences.bPESeqs)
			{
			if((Trim5 + Trim3 + ContamLen5PE2 + ContamLen3PE2  + (int)MinSeqLen) > (int)PE2ReadLen)
				{
				NumPE2Underlen += 1;
				continue;
				}
			}


		if(bIsPE1Fastq && PE1QCSchema > 0 && MinPhredScore > 0)
			{
			pszPE1QScoresBuff[PE1QSLen - Trim3 - ContamLen3PE1] = '\0';

			if(!MeetsMinPhredScoreThres(PE1QCSchema,MinPhredScore,(char *)&pszPE1QScoresBuff[Trim5 + ContamLen5PE1]))
				{
				NumPE1UnderPhredScore += 1;
				continue;
				}
			}

		if(m_Sequences.bPESeqs && bIsPE2Fastq && PE2QCSchema > 0 && MinPhredScore > 0)
			{
			pszPE2QScoresBuff[PE2QSLen - Trim3 - ContamLen3PE2] = '\0';
			if(!MeetsMinPhredScoreThres(PE2QCSchema,MinPhredScore,(char *)&pszPE2QScoresBuff[Trim5 + ContamLen5PE2]))
				{
				NumPE2UnderPhredScore += 1;
				continue;
				}
			}


		// trim 5' and 3' as needed
		pPE1Seq = pPE1RawReadsBuff;
		if(Trim5 || ContamLen5PE1)
			{
			pPE1Seq += (Trim5 + ContamLen5PE1);
			PE1ReadLen -= (Trim5 + ContamLen5PE1);
			}
		PE1ReadLen -= (Trim3 + ContamLen3PE1);
		if(TrimSeqLen > 0 && PE1ReadLen > (UINT32)TrimSeqLen)
			PE1ReadLen = TrimSeqLen;

		if(m_Sequences.bPESeqs)
			{
			pPE2Seq = pPE2RawReadsBuff;
			if(Trim5 || ContamLen5PE2)
				{
				pPE2Seq += (Trim5 + ContamLen5PE2);
				PE2ReadLen -= (Trim5 + ContamLen5PE2);
				}
			PE2ReadLen -= (Trim3 + ContamLen3PE2);
			if(TrimSeqLen > 0 && PE2ReadLen > (UINT32)TrimSeqLen)
				PE2ReadLen = (UINT32)TrimSeqLen;
			}

		// check for excessive number of Ns 
		UINT32 Idx;
		UINT32 NumNs = 0;		// number of indeterminate bases in last 100bp window
		etSeqBase *pBase = pPE1Seq;
		for(Idx = 0; Idx < PE1ReadLen; Idx++,pBase++)
			{
			if(Idx >= 100 && (pBase[-100] & 0x07) == eBaseN)
					NumNs -= 1;
			if((*pBase & 0x07) == eBaseN)
				NumNs += 1;
			if(NumNs > (UINT32)MaxNs)
				break;
			}

		if(NumNs > (UINT32)MaxNs ||
				(PE1ReadLen <= 100 && NumNs > ((UINT32)(MaxNs * PE1ReadLen)/100)))
			{
			NumPE1ExcessNs += 1;
			continue;
			}
		if(m_Sequences.bPESeqs)
			{
			NumNs = 0;		// number of indeterminate bases in last 100bp window
			pBase = pPE2Seq;
			for(Idx = 0; Idx < (UINT32)PE2ReadLen; Idx++,pBase++)
				{
				if(Idx >= 100 && (pBase[-100] & 0x07) == eBaseN)
						NumNs -= 1;
				if((*pBase & 0x07) == eBaseN)
					NumNs += 1;
				if(NumNs > (UINT32)MaxNs)
					break;
				}

			if(NumNs > (UINT32)MaxNs ||
					(PE2ReadLen <= 100 && NumNs > (UINT32)((MaxNs * PE2ReadLen)/100)))
				{
				NumPE2ExcessNs += 1;
				continue;
				}
			}
		
		if(SampleNth > 1 && (NumPE1ParsedReads % SampleNth))
			continue;
		
		if(CurAcceptedMinSeqLen == 0 || CurAcceptedMinSeqLen > (int)PE1ReadLen)
			CurAcceptedMinSeqLen = PE1ReadLen;
		if(CurAcceptedMaxSeqLen == 0 || CurAcceptedMaxSeqLen <  (int)PE1ReadLen)
			CurAcceptedMaxSeqLen = PE1ReadLen;

		if(m_Sequences.MinSeqLen == 0 || m_Sequences.MinSeqLen > PE1ReadLen)
			m_Sequences.MinSeqLen = PE1ReadLen;
		if(m_Sequences.MaxSeqLen == 0 || m_Sequences.MaxSeqLen <  PE1ReadLen)
			m_Sequences.MaxSeqLen = PE1ReadLen;
		if((Rslt=AddSeq(PE1FileID,m_Sequences.bPESeqs ? cFlgSeqPE : 0,PE1ReadLen,pPE1Seq)) !=eBSFSuccess)
			break;
		NumPE1AcceptedReads += 1;
		CurAcceptedTotSeqLen += (INT64)PE1ReadLen;
		if(m_Sequences.bPESeqs)
			{
			if(CurAcceptedMinSeqLen == 0 || CurAcceptedMinSeqLen > (int)PE2ReadLen)
				CurAcceptedMinSeqLen = PE2ReadLen;
			if(CurAcceptedMaxSeqLen == 0 || CurAcceptedMaxSeqLen <  (int)PE2ReadLen)
				CurAcceptedMaxSeqLen = PE2ReadLen;

			if(m_Sequences.MinSeqLen == 0 || m_Sequences.MinSeqLen > PE2ReadLen)
				m_Sequences.MinSeqLen = PE2ReadLen;
			if(m_Sequences.MaxSeqLen == 0 || m_Sequences.MaxSeqLen <  PE2ReadLen)
				m_Sequences.MaxSeqLen = PE2ReadLen;
			CSeqTrans::ReverseComplement(PE2ReadLen,pPE2Seq);
			if((Rslt=AddSeq(PE2FileID,cFlgSeqPE2 | cFlgSeqPE,PE2ReadLen,pPE2Seq)) !=eBSFSuccess)
				break;
			NumPE2AcceptedReads += 1;
			CurAcceptedTotSeqLen += (INT64)PE2ReadLen;
			}

		if(Zreads > 0 && NumPE1AcceptedReads >= Zreads)
			{
			Rslt = eBSFSuccess;
			break;
			}

		if((m_Sequences.Seqs2AssembOfs * m_Sequences.SeqWrdBytes) > cMaxConcatSeqLen)			// hit limit?
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Total length of loaded concatenated and packed sequences exceeds limit of %lld bytes",cMaxConcatSeqLen);
			Rslt = eBSFerrMaxDirEls;
			break;
			}

		continue;
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"raw sequence file '%s' processing error: %s ",pszPE1File, 
				Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
		delete pPE1RawReadsBuff;
		delete pszPE1QScoresBuff;
		FastaPE1.Close();
		if(pPE2RawReadsBuff != NULL)
			{
			delete pPE2RawReadsBuff;
			delete pszPE2QScoresBuff;
			FastaPE2.Close();
			}
		Reset(false);
		return(eBSFerrParse);
		}
	}
if(pPE1RawReadsBuff != NULL)
	delete pPE1RawReadsBuff;
if(pszPE1QScoresBuff != NULL)
	delete pszPE1QScoresBuff;
if(pPE2RawReadsBuff != NULL)
	delete pPE2RawReadsBuff;
if(pszPE2QScoresBuff != NULL)
	delete pszPE2QScoresBuff;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Parsed %d, Accepted %d %s",NumPE1ParsedReads,NumPE1AcceptedReads,m_Sequences.bPESeqs ? "paired sequences":"sequences");

if(Rslt != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors processing file: %s ",pszPE1File);
	while(FastaPE1.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,FastaPE1.GetErrMsg());
	FastaPE1.Close();
	if(m_Sequences.bPESeqs)
		FastaPE2.Close();
	Reset(false);
	return(Rslt);
	}
FastaPE1.Close();
if(m_Sequences.bPESeqs)
	FastaPE2.Close();

pPE1ReadFile->SeqsParsed = NumPE1ParsedReads;
pPE1ReadFile->SeqsUnderLen = NumPE1Underlen;
pPE1ReadFile->SeqsExcessNs = NumPE1ExcessNs;
pPE1ReadFile->SeqsAccepted = NumPE1AcceptedReads;

if(m_Sequences.bPESeqs)
	{
	pPE2ReadFile->SeqsParsed = NumPE2ParsedReads;
	pPE2ReadFile->SeqsUnderLen = NumPE2Underlen;
	pPE2ReadFile->SeqsExcessNs = NumPE2ExcessNs;
	pPE2ReadFile->SeqsAccepted = NumPE2AcceptedReads;
	}

m_Sequences.TotSeqsParsed += NumPE1ParsedReads + NumPE2ParsedReads;
m_Sequences.TotSeqsUnderLen += NumPE1Underlen + NumPE2Underlen;
m_Sequences.TotSeqsExcessNs += NumPE1ExcessNs + NumPE2ExcessNs; 
m_Sequences.TotSeqs2Assemb += NumPE1AcceptedReads + NumPE2AcceptedReads;
m_Sequences.NumPE1Seqs2Assemb += NumPE1AcceptedReads;
m_Sequences.NumPE2Seqs2Assemb += NumPE2AcceptedReads;
if(m_Sequences.TotSeqs2Assemb > 0)
	m_Sequences.MeanSeqLen = (double)m_Sequences.Seqs2AssembLen / m_Sequences.TotSeqs2Assemb;
else 
	m_Sequences.MeanSeqLen = 0.0;
m_Sequences.bSeqs2AssembDirty = true;	

if(NumPE1AcceptedReads + NumPE2AcceptedReads)
	CurAcceptedMeanSeqLen = (double)CurAcceptedTotSeqLen / ((double)NumPE1AcceptedReads + NumPE2AcceptedReads);
else
	CurAcceptedMeanSeqLen = 0.0;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Completed parsing %d, accepted %d %s, accepted sequences min length: %d max length: %d mean length: %1.1f",
								NumPE1ParsedReads,NumPE1AcceptedReads,m_Sequences.bPESeqs ? "paired sequences":"sequences",CurAcceptedMinSeqLen,CurAcceptedMaxSeqLen,CurAcceptedMeanSeqLen);

if(m_pContaminants != NULL)
	{
	if(m_Sequences.bPESeqs)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Removed %u paired sequences because at least one end was totally contained in a vector contaminate sequence", NumContamVector);
	else
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Removed %u sequences because these were totally contained in a vector contaminate sequence", NumContamVector);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Trimmed %u sequences because these flank overlapped onto a contaminate sequence", NumContamFlankTrim);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: total of %d sequences not accepted as under length %d, and %d had excessive indeterminate (%d Ns)",NumPE1Underlen + NumPE2Underlen,MinSeqLen,NumPE1ExcessNs+NumPE2ExcessNs,MaxNs);

if((PE1QCSchema > 0 || PE2QCSchema > 0) && MinPhredScore > 0)
	{
	if(!m_Sequences.bPESeqs)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: total of %d sequences had Phred scores under %d",NumPE1UnderPhredScore,MinPhredScore);
	else
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: total of %d (%d PE1, %d PE2) sequences had Phred scores under %d",NumPE1UnderPhredScore + NumPE2UnderPhredScore,NumPE1UnderPhredScore,NumPE2UnderPhredScore,MinPhredScore);
	}
else
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: No Phred score filtering");

return((teBSFrsltCodes)(NumPE1AcceptedReads + NumPE2AcceptedReads));
}



// LoadSeedPEs
// Load seed PE1 and PE2 sequence fragments from fasta file
teBSFrsltCodes
CKangadna::LoadSeedPEs(char *pszPE1File,		  // optional input high confidence seed PE1 sequences file
					   char *pszPE2File,		  // optional input high confidence seed PE2 sequences file
					   int OrientatePE,			  // PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
						int TrimEnds,			  // trim input sequences, both 5' and 3' ends by this many bases
						int MinInputSeqLen)		  // only accept for assembly sequences which are, after any trimming, of at least this length
{
teBSFrsltCodes Rslt;

UINT32 MinSeqLen;
UINT32 MaxSeqLen;
INT64 TotSeqLen;


bool bIsPE1Fastq;
int NumPE1DescrReads;
int NumPE1Underlen;
int NumPE1UnsupportedBases;
int NumPE1ExcessNs;
int NumPE1AcceptedReads;
int NumPE1ParsedReads;
UINT32 PE1ReadLen;
int PE1DescrLen;
UINT8 szPE1DescrBuff[1024];

bool bIsPE2Fastq;
int NumPE2DescrReads;
int NumPE2Underlen;
int NumPE2UnsupportedBases;
int NumPE2ExcessNs;
int NumPE2AcceptedReads;
int NumPE2ParsedReads;
UINT32 PE2ReadLen;
int PE2DescrLen;
UINT8 szPE2DescrBuff[1024];

etSeqBase *pPE1Seq;
etSeqBase *pPE2Seq;
UINT8 *pPE1RawReadsBuff;
UINT8 *pPE2RawReadsBuff;
CFasta FastaPE1;
CFasta FastaPE2;

tsReadFile *pPE1ReadFile;
tsReadFile *pPE2ReadFile;

UINT32 PE1FileID;
UINT32 PE2FileID;

int MaxAllowedFiles;

m_Sequences.bPESeqs = true;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadSeedPEs: Loading PE 5' from '%s' and PE 3' from '%s'",pszPE1File,pszPE2File);

// error if memory not preallocated to hold sequences
if(m_Sequences.pSeqs2Assemb == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Expected memory to have been preallocated for holding sequences");
	return(eBSFerrInternal);
	}

PE1FileID = m_NumRawFiles + 1;

pPE1ReadFile = &m_SrcFiles[m_NumRawFiles];
pPE2ReadFile = &m_SrcFiles[m_NumRawFiles+1];
PE2FileID = m_NumRawFiles + 2;
m_NumRawFiles += 2;

MaxAllowedFiles = 2 * cKDNAMaxInFileSpecs;
if(m_NumRawFiles > MaxAllowedFiles)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Too many files to process - max files: %d files...",MaxAllowedFiles);
	Reset(false);
	return(eBSFerrParams);
	}

memset(pPE1ReadFile,0,sizeof(tsReadFile));
pPE1ReadFile->FileID = PE1FileID;
strcpy((char *)pPE1ReadFile->szFile,pszPE1File);
pPE1ReadFile->PEFileID = PE2FileID;
memset(pPE2ReadFile,0,sizeof(tsReadFile));
pPE2ReadFile->FileID = PE2FileID;
pPE2ReadFile->PEFileID = PE1FileID;
strcpy((char *)pPE2ReadFile->szFile,pszPE2File);

if((pPE1RawReadsBuff = new UINT8 [cMaxRawSeqLen+16]) == NULL)					// 16 as a small safety margin!
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Unable to allocate memory for PE1 sequence buffering...");
	Reset(false);
	return(eBSFerrMem);
	}

if((pPE2RawReadsBuff = new UINT8 [cMaxRawSeqLen+16]) == NULL)				// 16 as a small safety margin!
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Unable to allocate memory for PE2 sequence buffering...");
	delete pPE1RawReadsBuff;
	Reset(false);
	return(eBSFerrMem);
	}

if((Rslt=(teBSFrsltCodes)FastaPE1.Open(pszPE1File,true,cMaxStageBuffSize))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Unable to open '%s' [%s] %s",pszPE1File,FastaPE1.ErrText((teBSFrsltCodes)Rslt),FastaPE1.GetErrMsg());
	delete pPE1RawReadsBuff;
	if(pPE2RawReadsBuff)
		delete pPE2RawReadsBuff;
	Reset(false);
	return(Rslt);
	}

if(FastaPE1.IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Unable to load '%s', SOLiD colorspace not supported",pszPE1File);
	delete pPE1RawReadsBuff;
	if(pPE2RawReadsBuff)
		delete pPE2RawReadsBuff;
	Reset(false);
	return(eBSFerrFileType);
	}
bIsPE1Fastq = FastaPE1.IsFastq();

if((Rslt=(teBSFrsltCodes)FastaPE2.Open(pszPE2File,true,cMaxStageBuffSize))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Unable to open '%s' [%s] %s",pszPE2File,FastaPE2.ErrText((teBSFrsltCodes)Rslt),FastaPE2.GetErrMsg());
	FastaPE1.Close();
	delete pPE1RawReadsBuff;
	delete pPE2RawReadsBuff;
	Reset(false);
	return(Rslt);
	}

if(FastaPE2.IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Unable to load '%s', SOLiD colorspace not supported",pszPE2File);
	FastaPE1.Close();
	FastaPE2.Close();
	delete pPE1RawReadsBuff;
	delete pPE2RawReadsBuff;
	Reset(false);
	return(eBSFerrFileType);
	}
bIsPE2Fastq = FastaPE2.IsFastq();
if(bIsPE1Fastq != 	bIsPE2Fastq)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Paired end file formats not of same type: '%s' is %s, '%s' is %s", 
									bIsPE1Fastq ? "Fastq" : "Fasta",pszPE2File,bIsPE2Fastq ? "Fastq" : "Fasta", pszPE2File);
	FastaPE1.Close();
	FastaPE2.Close();
	delete pPE1RawReadsBuff;
	delete pPE2RawReadsBuff;
	Reset(false);
	return(eBSFerrFileType);		
	}


UINT64 CurWorkSetSize;
CurWorkSetSize = m_Sequences.AllocMemSeqs2Assemb + m_Sequences.AllocMemSeqStarts + m_Sequences.AllocMemSfx + m_Sequences.AllocMemSeqFlags;

if(CurWorkSetSize != m_CurMaxMemWorkSetBytes)
	{
	if(!SetMaxMemWorkSetSize((size_t)CurWorkSetSize))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Unable to set maximum working set size (bytes) to %lld bytes",CurWorkSetSize);
		FastaPE1.Close();
		delete pPE1RawReadsBuff;
		if(pPE2RawReadsBuff)
			{
			FastaPE2.Close();
			delete pPE2RawReadsBuff;
			}
		Reset(false);
		return(eBSFerrMaxDirEls);
		}
	}

NumPE1DescrReads = 0;
NumPE1Underlen = 0;
NumPE1UnsupportedBases = 0;
NumPE1ExcessNs = 0;
NumPE1AcceptedReads = 0;
NumPE1ParsedReads = 0;
PE1ReadLen = 0;
PE1DescrLen = 0;

NumPE2DescrReads = 0;
NumPE2Underlen = 0;
NumPE2UnsupportedBases = 0;
NumPE2ExcessNs = 0;
NumPE2AcceptedReads = 0;
NumPE2ParsedReads = 0;
PE2ReadLen = 0;
PE2DescrLen = 0;

MinSeqLen = 0;
MaxSeqLen = 0;
TotSeqLen = 0;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadSeedPEs: Parsed 0, Accepted 0 paired sequences ...");
time_t Started = time(0);
while((Rslt = (teBSFrsltCodes)(PE1ReadLen = FastaPE1.ReadSequence(pPE1RawReadsBuff,cMaxRawSeqLen,true,false))) > eBSFSuccess)
	{
	NumPE1ParsedReads += 1;
	if(!(NumPE1ParsedReads % 1000))
		{
		time_t Now = time(0);
		unsigned long ElapsedSecs = (unsigned long) (Now - Started);
		if(ElapsedSecs >= 60)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadSeedPEs: Parsed %d, Accepted %d paired sequences",NumPE1ParsedReads,NumPE1AcceptedReads);
			Started = Now;
			}
		}

	if(PE1ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		NumPE1DescrReads += 1;
		PE1DescrLen = FastaPE1.ReadDescriptor((char *)szPE1DescrBuff,sizeof(szPE1DescrBuff)-1);
		szPE1DescrBuff[sizeof(szPE1DescrBuff)-1] = '\0'; 
		PE1ReadLen = FastaPE1.ReadSequence(pPE1RawReadsBuff,cMaxRawSeqLen);
		if(PE1ReadLen < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Problem parsing sequence after %d sequences parsed from file '%s'",NumPE1DescrReads,pszPE1File);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Last descriptor parsed: %s",szPE1DescrBuff);
			delete pPE1RawReadsBuff;
			FastaPE1.Close();
			if(m_Sequences.bPESeqs)
				{
				delete pPE2RawReadsBuff;
				FastaPE1.Close();
				}
			Reset(false);
			return(eBSFerrParse);
			}

		if(m_Sequences.bPESeqs)
			{
			Rslt = (teBSFrsltCodes)(PE2ReadLen = FastaPE2.ReadSequence(pPE2RawReadsBuff,cMaxRawSeqLen,true,false));
			if(Rslt <= eBSFSuccess)
				{
				if(Rslt < eBSFSuccess)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Errors processing file: %s ",pszPE2File);
					while(FastaPE2.NumErrMsgs())
						gDiagnostics.DiagOut(eDLFatal,gszProcName,FastaPE2.GetErrMsg());
					}
				else
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Fewer PE2 sequences in file %s than in PE1 file %s",pszPE2File, pszPE1File);
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Last descriptor parsed from PE1 file %s was : %s",pszPE1File,szPE1DescrBuff);
					}
				delete pPE1RawReadsBuff;
				FastaPE1.Close();
				delete pPE2RawReadsBuff;
				FastaPE2.Close();
				Reset(false);
				return(Rslt);
				}
			NumPE2ParsedReads += 1;
			if(PE2ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
				{
				NumPE2DescrReads += 1;
				PE2DescrLen = FastaPE2.ReadDescriptor((char *)szPE2DescrBuff,sizeof(szPE2DescrBuff)-1);
				szPE2DescrBuff[sizeof(szPE2DescrBuff)-1] = '\0'; 
				PE2ReadLen = FastaPE2.ReadSequence(pPE2RawReadsBuff,cMaxRawSeqLen);
				if(PE2ReadLen < 1)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Problem parsing sequence after %d sequences parsed from file '%s'",NumPE2DescrReads,pszPE2File);
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Last descriptor parsed: %s",szPE2DescrBuff);
					delete pPE1RawReadsBuff;
					FastaPE1.Close();
					delete pPE2RawReadsBuff;
					FastaPE2.Close();
					Reset(false);
					return(eBSFerrParse);
					}
				}
			else
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Raw sequence file '%s' processing error: %s ",pszPE2File, 
							Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
				delete pPE1RawReadsBuff;
				FastaPE1.Close();
				delete pPE2RawReadsBuff;
				FastaPE2.Close();
				Reset(false);
				return(eBSFerrParse);
				}
			}

		if(PE1ReadLen < (UINT32)(2 * TrimEnds) + MinInputSeqLen)			// after any end trimming require at least MinInputSeqLen before accepting read
			{
			NumPE1Underlen += 1;
			NumPE2Underlen += 1;
			continue;
			}
		if(PE2ReadLen <  (UINT32)(2 * TrimEnds) + MinInputSeqLen)
			{
			NumPE2Underlen += 1;
			NumPE1Underlen += 1;
			continue;
			}

		if(TrimEnds)
			{

			PE1ReadLen -=  (UINT32)(2 * TrimEnds);
			PE2ReadLen -=  (UINT32)(2 * TrimEnds);
			memmove(pPE1RawReadsBuff,&pPE1RawReadsBuff[TrimEnds],PE1ReadLen);
			memmove(pPE2RawReadsBuff,&pPE2RawReadsBuff[TrimEnds],PE2ReadLen);
			}

		pPE1Seq = pPE1RawReadsBuff;
		pPE2Seq = pPE2RawReadsBuff;

		// not interested in sequences containing indeterminate Ns 
		UINT32 Idx;
		etSeqBase *pBase = pPE1Seq;
		for(Idx = 0; Idx < PE1ReadLen; Idx++,pBase++)
			{
			if((*pBase & 0x07) == eBaseN)
				break;
			}

		if(Idx <  PE1ReadLen)
			{
			NumPE1ExcessNs += 1;
			continue;
			}

		pBase = pPE2Seq;
		for(Idx = 0; Idx < (UINT32)PE2ReadLen; Idx++,pBase++)
			{
			if((*pBase & 0x07) == eBaseN)
				break;
			}

		if(Idx <  PE2ReadLen)
			{
			NumPE2ExcessNs += 1;
			continue;
			}
				
		if(m_Sequences.MinSeqLen == 0 || m_Sequences.MinSeqLen > PE1ReadLen)
			m_Sequences.MinSeqLen = PE1ReadLen;
		if(MinSeqLen == 0 || MinSeqLen > PE1ReadLen)
			MinSeqLen = PE1ReadLen;
	
		if(m_Sequences.MaxSeqLen == 0 || m_Sequences.MaxSeqLen <  PE1ReadLen)
			m_Sequences.MaxSeqLen = PE1ReadLen;
		if(MaxSeqLen < PE1ReadLen)
			MaxSeqLen = PE1ReadLen;

		if(OrientatePE == 2 || OrientatePE == 3)	// MP Illumina circularised PE1 and MP SOLiD PE1 reads are antisense so need RevCpl back to sense
			CSeqTrans::ReverseComplement(PE1ReadLen,pPE1Seq);
		if((Rslt=AddSeq(PE1FileID,cFlgSeqPE,PE1ReadLen,pPE1Seq)) !=eBSFSuccess)
			break;
		NumPE1AcceptedReads += 1;

		if(m_Sequences.MinSeqLen == 0 || m_Sequences.MinSeqLen > PE2ReadLen)
			m_Sequences.MinSeqLen = PE2ReadLen;
		if(MinSeqLen == 0 || MinSeqLen > PE2ReadLen)
			MinSeqLen = PE2ReadLen;
	
		if(m_Sequences.MaxSeqLen == 0 || m_Sequences.MaxSeqLen <  PE2ReadLen)
			m_Sequences.MaxSeqLen = PE2ReadLen;
		if(MaxSeqLen < PE2ReadLen)
			MaxSeqLen = PE2ReadLen;

		if(OrientatePE == 0 || OrientatePE == 3)   // normal short insert PE2 and SOLiD PE2 reads are antisense so need RevCpl back to sense
			CSeqTrans::ReverseComplement(PE2ReadLen,pPE2Seq);
		if((Rslt=AddSeq(PE2FileID,cFlgSeqPE2 | cFlgSeqPE,PE2ReadLen,pPE2Seq)) !=eBSFSuccess)
			break;
		NumPE2AcceptedReads += 1;

		TotSeqLen += (PE1ReadLen +  PE2ReadLen);

		if((m_Sequences.Seqs2AssembOfs * m_Sequences.SeqWrdBytes) > cMaxConcatSeqLen)			// hit limit?
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Total length of loaded concatenated and packed sequences exceeds limit of %lld bytes",cMaxConcatSeqLen);
			Rslt = eBSFerrMaxDirEls;
			break;
			}

		continue;
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Raw sequence file '%s' processing error: %s ",pszPE1File, 
				Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
		delete pPE1RawReadsBuff;
		FastaPE1.Close();
		if(pPE2RawReadsBuff != NULL)
			{
			delete pPE2RawReadsBuff;
			FastaPE2.Close();
			}
		Reset(false);
		return(eBSFerrParse);
		}
	}
if(pPE1RawReadsBuff != NULL)
	delete pPE1RawReadsBuff;
if(pPE2RawReadsBuff != NULL)
	delete pPE2RawReadsBuff;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadSeedPEs: Parsed %d, Accepted %d paired sequences",NumPE1ParsedReads,NumPE1AcceptedReads);

if(Rslt != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedPEs: Errors processing files: %s %s",pszPE1File,pszPE2File);
	while(FastaPE1.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,FastaPE1.GetErrMsg());
	FastaPE1.Close();
	if(m_Sequences.bPESeqs)
		FastaPE2.Close();
	Reset(false);
	return(Rslt);
	}
FastaPE1.Close();
if(m_Sequences.bPESeqs)
	FastaPE2.Close();

pPE1ReadFile->SeqsParsed = NumPE1ParsedReads;
pPE1ReadFile->SeqsUnderLen = NumPE1Underlen;
pPE1ReadFile->SeqsExcessNs = NumPE1ExcessNs;
pPE1ReadFile->SeqsAccepted = NumPE1AcceptedReads;

if(m_Sequences.bPESeqs)
	{
	pPE2ReadFile->SeqsParsed = NumPE2ParsedReads;
	pPE2ReadFile->SeqsUnderLen = NumPE2Underlen;
	pPE2ReadFile->SeqsExcessNs = NumPE2ExcessNs;
	pPE2ReadFile->SeqsAccepted = NumPE2AcceptedReads;
	}

m_Sequences.TotSeqsParsed += NumPE1ParsedReads + NumPE2ParsedReads;
m_Sequences.TotSeqsUnderLen += NumPE1Underlen + NumPE2Underlen;
m_Sequences.TotSeqsExcessNs += NumPE1ExcessNs + NumPE2ExcessNs; 
m_Sequences.TotSeqs2Assemb += NumPE1AcceptedReads + NumPE2AcceptedReads;
m_Sequences.NumPE1Seqs2Assemb += NumPE1AcceptedReads;
m_Sequences.NumPE2Seqs2Assemb += NumPE2AcceptedReads;
if(m_Sequences.TotSeqs2Assemb > 0)
	m_Sequences.MeanSeqLen = (double)m_Sequences.Seqs2AssembLen / m_Sequences.TotSeqs2Assemb;
else 
	m_Sequences.MeanSeqLen = 0.0;
m_Sequences.bSeqs2AssembDirty = true;	

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadSeedPEs: Completed parsing %d, accepted %d paired sequences, sequences min length: %d max length: %d mean length: %d",
								NumPE1ParsedReads,NumPE1AcceptedReads,MinSeqLen,MaxSeqLen, (NumPE1AcceptedReads + NumPE2AcceptedReads) == 0 ? 0 : (int)(TotSeqLen/((INT64)NumPE1AcceptedReads + NumPE2AcceptedReads)));

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: A total of %d sequences were under length (%dbp) and %d contained indeterminate Ns",NumPE1Underlen + NumPE2Underlen, MinInputSeqLen,NumPE1ExcessNs+NumPE2ExcessNs);

return((teBSFrsltCodes)(NumPE1AcceptedReads + NumPE2AcceptedReads));
}


// LoadSeedContigs
// Load high confidence seed contigs from fasta
// Note that there are some limits - contigs less than cMinSeedContigLen or containing any indeterminate bases will be sloughed, and
// any contigs longer than cMaxSeedContigLen will be truncated to that length
//
teBSFrsltCodes
CKangadna::LoadSeedContigs(char *pszContigsFile,	// file containing fasta seed contigs
						int TrimEnds,				// trim input sequences, both 5' and 3' ends by this many bases
						int MinInputSeqLen)			// only accept for assembly sequences which are, after any trimming, of at least this length
{
teBSFrsltCodes Rslt;

int NumCtgNBlocks;
tsBlockNsLoci *pBlockNsLoci;

int NumDescrContigs;
int NumUnderlen;
int NumTrunclen;
int NumUnsupportedBases;
int NumExcessNs;
int NumAcceptedContigs;
int NumParsedContigs;
INT64 TotContigLen;
UINT32 ContigLen;

int MinContigLen;
int MaxContigLen;

int DescrLen;
UINT8 szDescrBuff[1024];

UINT8 *pRawContigsBuff;
CFasta Fasta;

tsReadFile *pContigsFile;

UINT32 FileID;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadSeedContigs: Loading seed contigs from '%s'",pszContigsFile);


// error if memory not preallocated to hold sequences
if(m_Sequences.pSeqs2Assemb == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected memory to have been preallocated for holding compacted sequences");
	return(eBSFerrInternal);
	}

FileID = m_NumRawFiles + 1;

pContigsFile = &m_SrcFiles[m_NumRawFiles];
m_NumRawFiles += 1;

if(m_NumRawFiles > ((cKDNAMaxInFileSpecs * 2) + 1))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Too many files to process - max files allowed: %d files...",((cKDNAMaxInFileSpecs * 2) + 1));
	Reset(false);
	return(eBSFerrParams);
	}

memset(pContigsFile,0,sizeof(tsReadFile));
pContigsFile->FileID = FileID;
strcpy((char *)pContigsFile->szFile,pszContigsFile);

if((pRawContigsBuff = new UINT8 [cMaxRawSeqLen+16]) == NULL)					// 16 as a small safety margin!
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for contig buffering...");
	Reset(false);
	return(eBSFerrMem);
	}

if((Rslt=(teBSFrsltCodes)Fasta.Open(pszContigsFile,true,cMaxStageBuffSize))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedContigs: Unable to open '%s' [%s] %s",pszContigsFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	delete pRawContigsBuff;
	Reset(false);
	return(Rslt);
	}

if(Fasta.IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedContigs: Unable to load contigs from '%s', SOLiD colorspace not supported",pszContigsFile);
	delete pRawContigsBuff;
	Reset(false);
	return(eBSFerrFileType);
	}

UINT64 CurWorkSetSize;
CurWorkSetSize = m_Sequences.AllocMemSeqs2Assemb + m_Sequences.AllocMemSeqStarts + m_Sequences.AllocMemSfx + m_Sequences.AllocMemSeqFlags;

if(CurWorkSetSize != m_CurMaxMemWorkSetBytes)
	{
	if(!SetMaxMemWorkSetSize((size_t)CurWorkSetSize))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedContigs: Unable to set maximum working set size (bytes) to %lld bytes",CurWorkSetSize);
		Fasta.Close();
		delete pRawContigsBuff;
		Reset(false);
		return(eBSFerrMaxDirEls);
		}
	}

NumDescrContigs = 0;
NumUnderlen = 0;
NumTrunclen = 0;
NumUnsupportedBases = 0;
NumExcessNs = 0;
NumAcceptedContigs = 0;
NumParsedContigs = 0;
ContigLen = 0;
DescrLen = 0;
MinContigLen = 0;
MaxContigLen = 0;
TotContigLen = 0;
NumCtgNBlocks = 0;
if(MinInputSeqLen < cMinSeedContigLen)
	MinInputSeqLen = cMinSeedContigLen;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadSeedContigs: Parsed 0, Accepted 0 ...");
time_t Started = time(0);
while((Rslt = (teBSFrsltCodes)(ContigLen = Fasta.ReadSequence(pRawContigsBuff,cMaxRawSeqLen,true,false))) > eBSFSuccess)
	{
	NumParsedContigs += 1;
	if(!(NumParsedContigs % 1000))
		{
		time_t Now = time(0);
		unsigned long ElapsedSecs = (unsigned long) (Now - Started);
		if(ElapsedSecs >= 60)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadSeedContigs: Parsed %d, Accepted %d",NumParsedContigs,NumAcceptedContigs);
			Started = Now;
			}
		}

	if(ContigLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		NumDescrContigs += 1;
		DescrLen = Fasta.ReadDescriptor((char *)szDescrBuff,sizeof(szDescrBuff)-1);
		szDescrBuff[sizeof(szDescrBuff)-1] = '\0'; 
		ContigLen = Fasta.ReadSequence(pRawContigsBuff,cMaxRawSeqLen);
		if(ContigLen < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedContigs: Problem parsing sequence after %d sequences parsed from file '%s'",NumDescrContigs,pszContigsFile);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedContigs: Last descriptor parsed: %s",szDescrBuff);
			delete pRawContigsBuff;
			Fasta.Close();
			Reset(false);
			return(eBSFerrParse);
			}

			// check on contig length
		if(ContigLen < (UINT32)(MinInputSeqLen + (2 * TrimEnds)))
			{
			NumUnderlen += 1;
			continue;
			}
		if(ContigLen - (UINT32)(2 * TrimEnds) > (UINT32)(cMaxSeedContigLen))
			{
			ContigLen = cMaxSeedContigLen;
			NumTrunclen += 1;
			}

		if(TrimEnds)
			{
			ContigLen -= (UINT32)(2 * TrimEnds);
			memmove(pRawContigsBuff,&pRawContigsBuff[TrimEnds],ContigLen);
			}		

		// check for if any indeterminates; 
		// if scaffolding record positions of indeterminates and substitute with a random base, in the scaffold reporting then these will be replaced back with indeterminates   
		UINT32 Idx;
		
		pBlockNsLoci = NULL;
		etSeqBase *pBase = pRawContigsBuff;
		for(Idx = 0; Idx < ContigLen; Idx++,pBase++)
			if((*pBase & 0x07) > eBaseT)
				{
				if(m_pBlockNsLoci == NULL)
					break;
				if(pBlockNsLoci == NULL)	// new indeterminate block starting
					{
					// ensure sufficent memory has been allocated to hold this indeterminates block...
					if(m_NumBlockNsLoci >= m_AllocdBlockNsLoci)		
						{
						// try to realloc with a 20% increase
						size_t memreq;
						UINT32 Blocks;
						void *pAllocd;
						Blocks = (UINT32)(((UINT64)m_NumBlockNsLoci * 120) / (UINT64)100);
						memreq = (size_t)Blocks * sizeof(tsBlockNsLoci);
#ifdef _WIN32
						pAllocd = realloc(m_pBlockNsLoci,memreq);
#else
						pAllocd = mremap(m_pBlockNsLoci,m_AllocdBlockNsLociSize,memreq,MREMAP_MAYMOVE);
						if(pAllocd == MAP_FAILED)
							pAllocd = NULL;
#endif
						if(pAllocd == NULL)
							{
							gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedContigs: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
							return(eBSFerrMem);
							}

						m_pBlockNsLoci = (tsBlockNsLoci *)pAllocd;
						m_AllocdBlockNsLociSize = memreq;
						m_AllocdBlockNsLoci = Blocks;
						}
					pBlockNsLoci = &m_pBlockNsLoci[m_NumBlockNsLoci++];
					pBlockNsLoci->Ofs = Idx;
					pBlockNsLoci->NumNs = 0;
					pBlockNsLoci->SeqID = m_Sequences.NumSeqs2Assemb + 1;
					}
				pBlockNsLoci->NumNs += 1;
				*pBase = rand() & 0x03;
				}
			else
				pBlockNsLoci = NULL;

		if(Idx < ContigLen)
			{
			NumExcessNs += 1;
			continue;
			}
		
		if(m_Sequences.MinSeqLen == 0 || m_Sequences.MinSeqLen > ContigLen)
			m_Sequences.MinSeqLen = ContigLen;
		if(MinContigLen == 0 || MinContigLen > (int)ContigLen)
			MinContigLen = (int)ContigLen;

		if(m_Sequences.MaxSeqLen == 0 || m_Sequences.MaxSeqLen <  ContigLen)
			m_Sequences.MaxSeqLen = ContigLen;
		if(MaxContigLen == 0 || MaxContigLen < (int)ContigLen)
			MaxContigLen = (int)ContigLen;
		TotContigLen += ContigLen;
		if((Rslt=AddSeq(FileID,0,ContigLen,pRawContigsBuff)) !=eBSFSuccess)
			break;

		NumAcceptedContigs += 1;
		
		if((m_Sequences.Seqs2AssembOfs * m_Sequences.SeqWrdBytes) > cMaxConcatSeqLen)			// hit limit?
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedContigs: Total length of loaded concatenated and packed sequences exceeds limit of %lld bytes",cMaxConcatSeqLen);
			Rslt = eBSFerrMaxDirEls;
			break;
			}

		continue;
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadSeedContigs: raw high confidence contigs file '%s' processing error: %s ",pszContigsFile, 
																	Rslt > 0 ? "over length contigs" : "not a multifasta contigs or fastq file");
		Fasta.Close();
		delete pRawContigsBuff;
		Reset(false);
		return(eBSFerrParse);
		}
	}
if(pRawContigsBuff != NULL)
	delete pRawContigsBuff;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadSeedContigs: Parsed %d, Accepted %d with %d truncated as overlength (overlength limit: %d) and %d rejected due to containing indeterminate Ns",
									NumParsedContigs,NumAcceptedContigs,NumTrunclen,cMaxSeedContigLen,NumExcessNs);

if(Rslt != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors processing file: %s ",pszContigsFile);
	while(Fasta.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,Fasta.GetErrMsg());
	Fasta.Close();
	Reset(false);
	return(Rslt);
	}
Fasta.Close();

pContigsFile->SeqsParsed = NumParsedContigs;
pContigsFile->SeqsUnderLen = NumUnderlen;
pContigsFile->SeqsExcessNs = NumExcessNs;
pContigsFile->SeqsAccepted = NumAcceptedContigs;

m_Sequences.TotSeqsParsed += NumParsedContigs;
m_Sequences.TotSeqsUnderLen += NumUnderlen;
m_Sequences.TotSeqsExcessNs += NumExcessNs; 
m_Sequences.TotSeqs2Assemb += NumAcceptedContigs;
m_Sequences.NumPE1Seqs2Assemb += NumAcceptedContigs;

if(m_Sequences.TotSeqs2Assemb > 0)
	m_Sequences.MeanSeqLen = (double)m_Sequences.Seqs2AssembLen / m_Sequences.TotSeqs2Assemb;
else 
	m_Sequences.MeanSeqLen = 0.0;
m_Sequences.bSeqs2AssembDirty = true;	

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadSeedContigs: Completed parsing %d, accepted %d, min length: %d max length: %d mean length: %d",
								NumParsedContigs,NumAcceptedContigs,MinContigLen,MaxContigLen,(int)((double)TotContigLen / NumAcceptedContigs));

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadSeedContigs: total of %d sequences were under length %d and %d contained indeterminate Ns",NumUnderlen,MinInputSeqLen,NumExcessNs);

return((teBSFrsltCodes)NumAcceptedContigs);
}

	// Imortant:
    // Normally Levenshtein distance uses penalty of 0 for matches, and 1 for mismatches, inserts and deletions so that
    // the distance represents the number of edits between two sequences to make both equal
    // Because of the extreme InDel biases with PacBio long reads different penalties can be specified to the Levenshtein generating functions
bool											// returns true if QuerySeq within MaxDist from RefSeq
CKangadna::IsLevenshteinAcceptedFwd(UINT32 KMerLen,// number (1..16) of packed bases in RefSeq and QuerySeq, bits 31/30 containing 5' base
							UINT32 RefSeq,		// calculate levenshtein distance between packed sequence (2bits/base) in RefSeq, bits 1/0 containing 3' base
						   UINT32 QuerySeq,     // and packed sequence (2bits/base) in QuerySeq, bits 1/0 containing 3' base
							UINT32 MaxLev,		// will be only accepting overlapping k-mers if their Levenshteins is <= MaxDist
							int PenaltyMatch,	// override the default match penalty (0)
							int PenaltyMisMatch,// override the default mismatch penalty (1)
							int PenaltyInsert,	// override the default insert penalty (1)
							int PenaltyDelete)	// override the default deletion penalty (1)
{
UINT32 Dist;
UINT32 BaseMsk;
if(KMerLen < 1 || KMerLen > 16)
	return(0);

if(m_pAcceptLevDist != NULL && KMerLen == m_LevDistKMerLen && MaxLev == m_MaxLev &&
			PenaltyMatch == m_PenaltyMatch && PenaltyMisMatch == m_PenaltyMisMatch && PenaltyInsert == m_PenaltyInsert && PenaltyDelete == m_PenaltyDelete)
	{
	BaseMsk = 0x0ffff >> 2 * (8 - KMerLen);
	RefSeq &= BaseMsk;
	QuerySeq &= BaseMsk;
	return((m_pAcceptLevDist[(RefSeq * m_AcceptLevWidth) + (QuerySeq/8)] & (0x01 << (QuerySeq % 8))) == 0 ? false : true);
	}
Dist = GetLevenshteinDistFwd(KMerLen,RefSeq,QuerySeq,PenaltyMatch,PenaltyMisMatch,PenaltyInsert,PenaltyDelete);
return(Dist >= 0 && Dist <= MaxLev ? true : false);
}

#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))


int												// returned levenshtein distance between RefSeq and QuerySeq, -1 if errors
CKangadna::GetLevenshteinDistFwd(UINT32 KMerLen,	// number (1..16) of packed bases in both RefSeq and QuerySeq, bits 1/0 containing 3' base 
							UINT32 RefSeq,		// calculate levenshtein distance between packed sequence (2bits/base) in RefSeq, bits 1/0 containing 3' base
						   UINT32 QuerySeq,     // and packed sequence (2bits/base) in QuerySeq, bits 1/0 containing 3' base
							int PenaltyMatch,	// override the default match penalty (0)
							int PenaltyMisMatch,// override the default mismatch penalty (1)
							int PenaltyInsert,	// override the default insert penalty (1)
							int PenaltyDelete)	// override the default deletion penalty (1)
{
UINT32 x, y, lastdiag, olddiag;
UINT32 column[16+1];
UINT32 TmpSeq;
UINT32 BaseMsk;

if(KMerLen < 1 || KMerLen > 16)
	return(-1);

BaseMsk = 0x03 << 2 * (KMerLen - 1);

for (y = 1; y <= KMerLen; y++)
	column[y] = y;

for (x = 1; x <= KMerLen; x++) 
	{
    column[0] = x;
	TmpSeq =RefSeq; 
    for (y = 1, lastdiag = x-1; y <= KMerLen; y++) 
		{
        olddiag = column[y];
        column[y] = MIN3(column[y] + PenaltyInsert, column[y-1] + PenaltyDelete, lastdiag + ((TmpSeq & BaseMsk) == (QuerySeq & BaseMsk) ? PenaltyMatch : PenaltyMisMatch)); 
        lastdiag = olddiag;
		TmpSeq <<= 2;
        }
	QuerySeq <<= 2;
    }
return((int)column[KMerLen]);
}

int												// returned levenshtein distance between RefSeq and QuerySeq, -1 if errors
CKangadna::GetLevenshteinDistFwd(UINT32 KMerLen,	// number (1..cLevenshteinMaxSeqLen) of bases in both pRefSeq and pQuerySeq 
							etSeqBase *pRefSeq,	// calculate levenshtein distance between this sequence (pRefSeq[0] contains 5' base)
						   etSeqBase *pQuerySeq, // and sequence in pQuerySeq
							int PenaltyMatch,	// override the default match penalty (0)
							int PenaltyMisMatch,// override the default mismatch penalty (1)
							int PenaltyInsert,	// override the default insert penalty (1)
							int PenaltyDelete,	// override the default deletion penalty (1)
							int MaxExpDist)     // set to max expected distance to early terminate if observed distance would be greater; set 0 if no expected maximum
{
UINT32 x, y, lastdiag, olddiag;
UINT32 column[cLevenshteinMaxSeqLen+1];
etSeqBase *pQueryBase;
UINT8 RefBase;
UINT8 QueryBase;

if(KMerLen < 1 || KMerLen > cLevenshteinMaxSeqLen)
	return(-1);

for (y = 1; y <= KMerLen; y++)
	column[y] = y;

for (x = 1; x <= KMerLen; x++) 
	{
    column[0] = x;
	pQueryBase = pQuerySeq;
	if((RefBase = *pRefSeq++ & 0x07) > eBaseT)
		{
		RefBase = eBaseEOS;
		pRefSeq -= 1;
		}
    for (y = 1, lastdiag = x-1; y <= KMerLen; y++) 
		{
        olddiag = column[y];
		if((QueryBase = *pQueryBase++ & 0x07) > eBaseT)
			{
			QueryBase = eBaseEOS;
			pQuerySeq -= 1;
			}
        column[y] = MIN3(column[y] + PenaltyInsert, column[y-1] + PenaltyDelete, lastdiag + (RefBase == QueryBase ? PenaltyMatch : PenaltyMisMatch));
        lastdiag = olddiag;
        }
	if(MaxExpDist > 0 && column[KMerLen] > (UINT32)MaxExpDist && (KMerLen - x) < (column[KMerLen] - (UINT32)MaxExpDist))
		return((int)column[KMerLen]);
    }
return((int)column[KMerLen]);
}

int			// returned levenshtein distance between RefSeq and QuerySeq, -1 if errors, assumes pRefSeq, pRefSeq point to last 3' bases 
CKangadna::GetLevenshteinDistRev(UINT32 KMerLen,	// number (1..cLevenshteinMaxSeqLen) of bases in both pRefSeq and pQuerySeq 
							etSeqBase *pRefSeq,	// calculate levenshtein distance between this sequence (pRefSeq[0] contains 3' base)
						   etSeqBase *pQuerySeq, // and sequence in pQuerySeq (pQuerySeq[0] contains 3' base)
							int PenaltyMatch,	// override the default match penalty (0)
							int PenaltyMisMatch,// override the default mismatch penalty (1)
							int PenaltyInsert,	// override the default insert penalty (1)
							int PenaltyDelete,	// override the default deletion penalty (1)
							int MaxExpDist)     // set to max expected distance to early terminate if observed distance would be greater; set 0 if no expected maximum
{
UINT32 x, y, lastdiag, olddiag;
UINT32 column[cLevenshteinMaxSeqLen+1];
etSeqBase *pQueryBase;
UINT8 RefBase;
UINT8 QueryBase;

if(KMerLen < 1 || KMerLen > cLevenshteinMaxSeqLen)
	return(-1);

for (y = 1; y <= KMerLen; y++) 
	column[y] = y;

for (x = 1; x <= KMerLen; x++) 
	{
    column[0] = x;
	pQueryBase = pQuerySeq;
	if((RefBase = *pRefSeq-- & 0x07) > eBaseT)
		{
		RefBase = eBaseEOS;
		pRefSeq += 1;
		}
    for (y = 1, lastdiag = x-1; y <= KMerLen; y++) 
		{
        olddiag = column[y];
		if((QueryBase = *pQueryBase-- & 0x07) > eBaseT)
			{
			QueryBase = eBaseEOS;
			pQuerySeq += 1;
			}
        column[y] = MIN3(column[y] + PenaltyInsert, column[y-1] + PenaltyDelete, lastdiag + (RefBase == QueryBase ? PenaltyMatch : PenaltyMisMatch));
        lastdiag = olddiag;
        }
	if(MaxExpDist > 0 && column[KMerLen] > (UINT32)MaxExpDist && (KMerLen - x) < (column[KMerLen] - (UINT32)MaxExpDist))
		return((int)column[KMerLen]);
    }
return((int)column[KMerLen]);
}

// modified levenshtein to account for expected error profile with pacbio reads
// whereby substitutions are relatively rare compared to micoInDels, and microInDels are affine gap scored
// 
const UINT32 cFlgAGopened = 0x80000000;			// set if score in lower bits was the result of insert opening


int												// returned affine gap levenshtein distance between RefSeq and QuerySeq, -1 if errors
CKangadna::GetAGLevenshteinDistFwd(UINT32 KMerLen,	// number (1..cLevenshteinMaxSeqLen) of bases in both pRefSeq and pQuerySeq 
							etSeqBase *pRefSeq,	// calculate levenshtein distance between this sequence (pRefSeq[0] contains 5' base)
						   etSeqBase *pQuerySeq, // and sequence in pQuerySeq
							int PenaltyMatch,	// override the default match penalty (0)
							int PenaltyMisMatch,// override the default mismatch penalty (10)
							int PenaltyInDelOpen,	// override the default indel gap open (5 to open)
							int PenaltyInDelExtd)	// override the default indel gap extension (1 per base extension)
{
UINT32 x, y, lastdiag, olddiag;
UINT32 column[cLevenshteinMaxSeqLen+1];
etSeqBase *pQueryBase;
UINT8 RefBase;
UINT8 QueryBase;

if(KMerLen < 1 || KMerLen > cLevenshteinMaxSeqLen)
	return(-1);

for (y = 1; y <= KMerLen; y++)     // initialise with affine gap scoring
	column[y] = (((y-1) * PenaltyInDelExtd) + (UINT32)PenaltyInDelOpen) | cFlgAGopened;

for (x = 1; x <= KMerLen; x++) 
	{
    column[0] = (((x-1) * PenaltyInDelExtd) + (UINT32)PenaltyInDelOpen) | cFlgAGopened;
	pQueryBase = pQuerySeq;
	if((RefBase = *pRefSeq++ & 0x07) > eBaseT)
		{
		RefBase = eBaseEOS;
		pRefSeq -= 1;
		}

	lastdiag = x == 1 ? 0 : (UINT32)PenaltyInDelOpen + ((x-2) * PenaltyInDelExtd);
    for (y = 1; y <= KMerLen; y++) 
		{
        olddiag = column[y] & ~cFlgAGopened;
		if((QueryBase = *pQueryBase++ & 0x07) > eBaseT)
			{
			QueryBase = eBaseEOS;
			pQuerySeq -= 1;
			}

		UINT32 PutInsert = (column[y] & ~cFlgAGopened) + (column[y] & cFlgAGopened ? PenaltyInDelExtd : PenaltyInDelOpen);
		UINT32 PutDelete = (column[y-1] & ~cFlgAGopened) + (column[y-1] & cFlgAGopened ? PenaltyInDelExtd : PenaltyInDelOpen);
		UINT32 PutDiag = lastdiag + (RefBase == QueryBase ? PenaltyMatch : PenaltyMisMatch);

	    if(PutInsert < PutDelete && PutInsert < PutDiag)
			column[y] = PutInsert | cFlgAGopened;
		else
			if(PutDelete < PutDiag)
				column[y] = PutDelete | cFlgAGopened;
			else
				column[y] = PutDiag;

        lastdiag = olddiag;
        }
    }
return((int)(column[KMerLen] & ~cFlgAGopened));
}


#ifdef _WIN32
unsigned __stdcall ProcGenLevsFwdThread(void * pThreadPars)
#else
void *ProcGenLevsThread(void * pThreadPars)
#endif
{
int Rslt;
tsThreadPregenLevPars *pPars = (tsThreadPregenLevPars *)pThreadPars;			// makes it easier not having to deal with casts!
CKangadna *pKangadna = (CKangadna *)pPars->pThis;

Rslt = pKangadna->ProcGenLevsFwd(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
CKangadna::LoadPregenLevenshteinsFwd(char *pszLevenshteinsFile, // load pregenerated Levenshteins from this file
								UINT32 KMerLen,		// expecting levenshtein distances for all k-mers to be of this length (must be in range 4..8)
								UINT32 MaxLev,		// expecting that will be accepting overlapping k-mers if their Levenshteins is <= MaxLev 
							int PenaltyMatch,	// expecting match penalty (0)
							int PenaltyMisMatch,// expecting mismatch penalty (1)
							int PenaltyInsert,	// expecting insert penalty (1)
							int PenaltyDelete)	// expecting deletion penalty (1)
{
int hFile;
tsPregenLevenshteinsHdr Hdr;

if(m_pAcceptLevDist != NULL)
	{
	delete m_pAcceptLevDist;
	m_pAcceptLevDist = NULL;
	}
m_LevDistKMerLen = 0;
m_MaxLev = 0;
m_LevKMers = 0;
m_AcceptLevWidth = 0;
m_AllocAcceptLevDist = 0;

hFile = open(pszLevenshteinsFile,O_READSEQ);
if(hFile == -1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadPregenLevenshteins: Unable to open %s - %s",pszLevenshteinsFile,strerror(errno));
	return(eBSFerrOpnFile);
	}

if(_lseeki64(hFile,0,SEEK_SET)!=0)			
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadPregenLevenshteins: Seek failed to offset 0 - %s",pszLevenshteinsFile,strerror(errno));
	close(hFile);
	return(eBSFerrFileAccess);
	}

// read in header..
if(sizeof(tsPregenLevenshteinsHdr) != read(hFile,&Hdr,sizeof(tsPregenLevenshteinsHdr)))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadPregenLevenshteins: Failed reading file header - %s",pszLevenshteinsFile);
	close(hFile);
	return(eBSFerrFileAccess);	
	}

// check it is the Levenshteins file header
if(!(Hdr.Signiture[0] == 'l' && Hdr.Signiture[1] == 'e' && Hdr.Signiture[2] == 'v' && Hdr.Signiture[3] == 'd'))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadPregenLevenshteins: Failed reading file header - %s",pszLevenshteinsFile);
	close(hFile);
	return(eBSFerrFileAccess);	
	}

// and that the Levenshteins were generated with the requested KMerLen, MaxLev, PenaltyMatch, PenaltyMisMatch, PenaltyInsert, PenaltyDelete
if(!(Hdr.KMerLen == KMerLen && Hdr.MaxLev == MaxLev && Hdr.PenaltyMatch == PenaltyMatch && Hdr.PenaltyMisMatch == PenaltyMisMatch && Hdr.PenaltyInsert == PenaltyInsert && Hdr.PenaltyDelete == PenaltyDelete))
	{
	close(hFile);
	return(eBSFerrFileAccess);	
	}


if((m_pAcceptLevDist = new UINT8 [Hdr.m_AllocAcceptLevDist]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadPregenLevenshteins: Failed to allocate memory - %s",pszLevenshteinsFile);
	close(hFile);
	return(eBSFerrFileAccess);	
	}

m_LevDistKMerLen = Hdr.KMerLen;
m_LevKMers = Hdr.LevKMers;
m_MaxLev = Hdr.MaxLev;
m_PenaltyMatch = Hdr.PenaltyMatch;
m_PenaltyMisMatch = Hdr.PenaltyMisMatch;
m_PenaltyInsert = Hdr.PenaltyInsert;
m_PenaltyDelete = Hdr.PenaltyDelete;
m_AcceptLevWidth = Hdr.m_AcceptLevWidth;
m_AllocAcceptLevDist = Hdr.m_AllocAcceptLevDist;

if(m_AllocAcceptLevDist != read(hFile,m_pAcceptLevDist,m_AllocAcceptLevDist))
	{
	close(hFile);
	return(eBSFerrFileAccess);
	}
close(hFile);
return(0);
}

int
CKangadna::SavePregenLevenshteinsFwd(char *pszLevenshteinsFile) // save pregenerated Levenshteins to this file
{
teBSFrsltCodes Rslt;
int hFile;
tsPregenLevenshteinsHdr Hdr;
#ifdef _WIN32
if((hFile = open(pszLevenshteinsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hFile = open(pszLevenshteinsFile, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",pszLevenshteinsFile,strerror(errno));
	return(eBSFerrCreateFile);
	}
if(_lseeki64(hFile,0,SEEK_SET)!=0)			
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Seek failed to offset 0 - %s",pszLevenshteinsFile,strerror(errno));
	close(hFile);
	return(eBSFerrFileAccess);
	}
// initialise header..
Hdr.Signiture[0] = 'l';
Hdr.Signiture[1] = 'e';
Hdr.Signiture[2] = 'v';
Hdr.Signiture[3] = 'd';
Hdr.KMerLen	= (UINT32)m_LevDistKMerLen;
Hdr.LevKMers = m_LevKMers;
Hdr.MaxLev = m_MaxLev;
Hdr.PenaltyMatch = m_PenaltyMatch;
Hdr.PenaltyMisMatch = m_PenaltyMisMatch;
Hdr.PenaltyInsert = m_PenaltyInsert;
Hdr.PenaltyDelete = m_PenaltyDelete;
Hdr.m_AcceptLevWidth = m_AcceptLevWidth;
Hdr.m_AllocAcceptLevDist = m_AllocAcceptLevDist;

if((Rslt = ChunkedWrite(hFile,pszLevenshteinsFile,0,(UINT8 *)&Hdr,sizeof(tsPregenLevenshteinsHdr)))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"SavePregenLevenshteins: Write header to file %s failed",pszLevenshteinsFile);
	return(Rslt);
	}

if((Rslt = ChunkedWrite(hFile,pszLevenshteinsFile,sizeof(tsPregenLevenshteinsHdr),(UINT8 *)m_pAcceptLevDist,m_AllocAcceptLevDist))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"SavePregenLevenshteins: Write to file %s failed",pszLevenshteinsFile);
	close(hFile);
	return(Rslt);
	}

#ifdef _WIN32
_commit(hFile);
#else
fsync(hFile);
#endif
close(hFile);
return(eBSFSuccess);
}



UINT32											// total number of k-mer instances which have Levenshteins is <= MaxLev  
CKangadna::PregenLevenshteinsFwd(UINT32 KMerLen,	// generate levenshtein distances for all k-mers of this length (must be in range 4 to 8)
						 UINT32 MaxLev,			// will be only accepting overlapping k-mers if their Levenshteins is <= MaxLev
							int PenaltyMatch,	// override the default match penalty (0)
							int PenaltyMisMatch,// override the default mismatch penalty (1)
							int PenaltyInsert,	// override the default insert penalty (1)
							int PenaltyDelete)	// override the default deletion penalty (1)
{
tsThreadPregenLevPars Threads[cMaxWorkerThreads];
tsThreadPregenLevPars *pThread;
int KMersPerThread;
int ThreadIdx;
int NumThreads;
UINT32 StartKMer;
UINT32 NumAccepted;
UINT32 LevKMers;
if(KMerLen < 4 || KMerLen > 8)
	return(0);

if(m_pAcceptLevDist != NULL && (m_LevDistKMerLen != KMerLen || m_MaxLev != MaxLev))
	{
	delete m_pAcceptLevDist;
	m_pAcceptLevDist = NULL;
	m_LevDistKMerLen = 0;
	m_MaxLev = 0;
	m_PenaltyMatch = 0;
	m_PenaltyMisMatch = 0;
	m_PenaltyInsert = 0;
	m_PenaltyDelete = 0;
	m_AcceptLevWidth = 0;
	}
LevKMers =  0x01 << (KMerLen * 2);
m_AcceptLevWidth = (LevKMers + 7) / 8;
m_LevDistKMerLen = KMerLen; 
m_MaxLev = MaxLev;
m_LevKMers = LevKMers;
m_PenaltyMatch = PenaltyMatch;
m_PenaltyMisMatch = PenaltyMisMatch;
m_PenaltyInsert = PenaltyInsert;
m_PenaltyDelete = PenaltyDelete;
if(m_pAcceptLevDist == NULL)
	{
	m_AllocAcceptLevDist = ((size_t)LevKMers * m_AcceptLevWidth);
	if((m_pAcceptLevDist = new UINT8 [m_AllocAcceptLevDist])==NULL)
		return(eBSFerrMem);
	memset(m_pAcceptLevDist,0,m_AllocAcceptLevDist);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating Levenshteins for K-mer %d",KMerLen);

KMersPerThread = LevKMers / m_NumThreads;
if(KMersPerThread < 2)
	KMersPerThread = LevKMers;
 
StartKMer = 0;
NumThreads = 0;
pThread = Threads;
for(ThreadIdx = 0; ThreadIdx < m_NumThreads && StartKMer < LevKMers; ThreadIdx++,pThread++)
	{
	pThread->ThreadIdx = ThreadIdx;
	pThread->pThis = this;
	pThread->NumAccepted = 0;
	pThread->NumRejected = 0;
	pThread->MaxLev = MaxLev;
	pThread->KMerLen = KMerLen;
	pThread->LevKMers = LevKMers;
	pThread->PenaltyMatch = PenaltyMatch;
	pThread->PenaltyMisMatch = PenaltyMisMatch;
	pThread->PenaltyInsert = PenaltyInsert;
	pThread->PenaltyDelete = PenaltyDelete;
	pThread->StartKMer = StartKMer;
	StartKMer += KMersPerThread;
	if(StartKMer >= LevKMers || ThreadIdx + 1 == m_NumThreads)
		StartKMer = LevKMers;
	pThread->EndKMer = StartKMer - 1;
	pThread->Rslt = 0;
	NumThreads += 1;
#ifdef _WIN32
	pThread->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, ProcGenLevsFwdThread, pThread, 0, &pThread->threadID);
#else
	pThread->threadRslt = pthread_create(&pThread->threadID, NULL, ProcGenLevsThread, pThread);
#endif
	}

// allow threads a few seconds to startup
#ifdef _WIN32
Sleep(2000);
#else
sleep(2);
#endif
NumAccepted = 0;
pThread = Threads;
for (ThreadIdx = 0; ThreadIdx < NumThreads; ThreadIdx++, pThread++)
	{
#ifdef _WIN32
	while (WAIT_TIMEOUT == WaitForSingleObject(pThread->threadHandle, 60000)){};
	CloseHandle(pThread->threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60;
	while ((JoinRlt = pthread_timedjoin_np(pThread->threadID, NULL, &ts)) != 0)
		{
		ts.tv_sec += 60;
		}
#endif
	NumAccepted += pThread->NumAccepted;
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed generating Levenshteins");
return(NumAccepted);
}



int
CKangadna::ProcGenLevsFwd(tsThreadPregenLevPars *pPars) 
{
UINT32 RefSeq;
UINT32 QuerySeq;
int Dist;
for(RefSeq = pPars->StartKMer; RefSeq <= pPars->EndKMer; RefSeq++)
	{
	for(QuerySeq = 0; QuerySeq < pPars->LevKMers; QuerySeq++)
		{
		Dist = GetLevenshteinDistFwd(pPars->KMerLen,RefSeq,QuerySeq,pPars->PenaltyMatch,pPars->PenaltyMisMatch,pPars->PenaltyInsert,pPars->PenaltyDelete) > (int)pPars->MaxLev ? 0 : 1;

		m_pAcceptLevDist[(RefSeq * m_AcceptLevWidth) + (QuerySeq/8)] |= Dist << (QuerySeq % 8);
		if(Dist)
			pPars->NumAccepted += 1;
		else
			pPars->NumRejected += 1;
		}
	}
return(0);
}



teBSFrsltCodes
CKangadna::LoadLongReads(char *pszLongReadsFile) // file holding long reads (multifasta or fastq)
{
teBSFrsltCodes Rslt;
UINT32 EstNumSeqs;
UINT64 CumulativeMemory;
int SeqWrdBytes;
SeqWrdBytes = GetSeqWrdBytes();

if((Rslt = EstMemReq(pszLongReadsFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to estimate memory required for '%s",pszLongReadsFile);
	Reset(false);
	return(eBSFerrOpnFile);	// treat as though unable to open file
	}
CumulativeMemory = EstMemReq(0,0,0,0,0);
EstNumSeqs = GetEstSeqsToProc();

CumulativeMemory += EstNumSeqs * 12; // very rough estimate allowing for sparse suffix and flags requirements
if(CumulativeMemory < 1000000000)	// 100M
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Estimated total cumulative minimum required memory: %1.1f MB",(double)CumulativeMemory /1000000);
else
	{
	if(CumulativeMemory < 1000000000) // 1000M
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Estimated total cumulative minimum required memory: %1.0f MB",(double)CumulativeMemory /1000000);
	else
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Estimated total cumulative minimum required memory: %1.1f GB",(double)CumulativeMemory /1000000000);
	}

	// allocate to hold est cumulative memory upfront - may as well know now rather than later if there is insufficent memory...
if((Rslt = AllocSeqs2AssembMem((CumulativeMemory * 110)/100))!= eBSFSuccess)	// add 10% to reduce risk of later having to realloc....
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to continue");
	Reset(false);
	return(Rslt);
	}


if((Rslt = AllocBlockNsLoci(EstNumSeqs))!= eBSFSuccess)	
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to continue");
	Reset(false);
	return(Rslt);
	}


if((Rslt=LoadSeedContigs(pszLongReadsFile)) <= 0)
	{
	Reset(false);
	return(Rslt);
	}

// initialise header flags
gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadLongReads: Initialising sequence headers ...");
UpdateAllSeqHeaderFlags(0,~(cFlgNonOverlap),false);			// any PEs are now single ended
				

if((Rslt=GenRdsSfx()) < eBSFSuccess)	// generate index over all 15base words
	return((teBSFrsltCodes)Rslt);
	
if((Rslt=GenSeqStarts(true,false)) < eBSFSuccess)	// generate array of sequence starts plus array of flags from sequence headers
	return((teBSFrsltCodes)Rslt);
return(Rslt);
}



// AddSeq
// Add sequences which are to be subsequently to be suffix indexed and assembled
// Sequence bases are packed (2bits per base) into 32bits dependent on m_SeqWrdBytes
teBSFrsltCodes
CKangadna::AddSeq(UINT32 SrcFileID,				// 8bits identifies source file from which this sequence was processed
						UINT32 SeqFlags,    	// 16bits as bit flags
						UINT32 SeqLen,			// 30bit sequence length
						UINT8 *pSeq)			// ptr to sequence
{
UINT64 MemMinReq;
if(SeqLen < 1 || SeqLen > cMaxRawSeqLen || pSeq == NULL)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"AddSeq: SeqLen: %d pSeq is %s",SeqLen,pSeq == NULL ? "NULL" : "none-null");
	return(eBSFerrParams);
	}

AcquireSerialise();

MemMinReq =  4 * (1000 + m_Sequences.Seqs2AssembOfs + ((SeqLen + 14) / 15));

// ensure sufficent memory has been allocated to hold this sequence...
if(MemMinReq > m_Sequences.AllocMemSeqs2Assemb)		
	{
	// try to realloc memory with a 20% increase
	size_t memreq;
	void *pAllocd;
	memreq = (size_t)((MemMinReq * 120) / (UINT64)100);
#ifdef _WIN32
	pAllocd = realloc(m_Sequences.pSeqs2Assemb,memreq);
#else
	pAllocd = mremap(m_Sequences.pSeqs2Assemb,m_Sequences.AllocMemSeqs2Assemb,memreq,MREMAP_MAYMOVE);
	if(pAllocd == MAP_FAILED)
		pAllocd = NULL;
#endif
	if(pAllocd == NULL)
		{
		ReleaseSerialise();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeq: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}

	m_Sequences.pSeqs2Assemb = pAllocd;
	m_Sequences.AllocMemSeqs2Assemb = (UINT64)memreq;
	}

tSeqWrd4 *pSeqWrds;
tSeqWrd4 *pPackSeq;

pSeqWrds = (tSeqWrd4 *)m_Sequences.pSeqs2Assemb;
pPackSeq = &pSeqWrds[m_Sequences.Seqs2AssembOfs];
if(m_Sequences.Seqs2AssembOfs == 0)			// first sequence, so preceed with a cSeqWrd1BOS
	{
	*pPackSeq++ = cSeqWrd4BOS;
	m_Sequences.Seqs2AssembOfs += 1;
	}
m_Sequences.NumSeqs2Assemb += 1;
pPackSeq = (tSeqWrd4 *)SetSeqHeader(pPackSeq,m_Sequences.NumSeqs2Assemb,SrcFileID,SeqFlags,SeqLen,NULL);
m_Sequences.Seqs2AssembOfs += GenPackedSeqWrds(SeqLen,pSeq,0,pPackSeq);
pSeqWrds[m_Sequences.Seqs2AssembOfs] = cSeqWrd4EOS;			// currently last sequence, so finish with a cSeqWrd4EOS

m_Sequences.Seqs2AssembLen += SeqLen;
m_Sequences.bSeqs2AssembDirty;
ReleaseSerialise();
return(eBSFSuccess);
}

// AddPESeqs
// Add paired end sequences which are to be subsequently to be suffix indexed and assembled
// Sequence bases are packed (2bits per base) into 32bits dependent on m_SeqWrdBytes
teBSFrsltCodes
CKangadna::AddPESeqs(UINT32 SrcFileID,			// 8bits identifies source file from which this sequence was processed
						UINT32 PE1SeqFlags,    	// PE1 16bits as bit flags
						UINT32 PE1SeqLen,		// PE1 30bit sequence length
						UINT8 *pPE1Seq,			// PE1 ptr to sequence
						UINT32 PE2SeqFlags,    	// PE2 16bits as bit flags
						UINT32 PE2SeqLen,		// PE2 30bit sequence length
						UINT8 *pPE2Seq)			// PE2 ptr to sequence
{
UINT64 MemMinReq;
if(PE1SeqLen < 1 || PE1SeqLen > cMaxRawSeqLen || pPE1Seq == NULL)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"AddPESeqs: PE1SeqLen: %d pPE1Seq is %s",PE1SeqLen,pPE1Seq == NULL ? "NULL" : "not-null");
	return(eBSFerrParams);
	}
if(PE2SeqLen < 1 || PE2SeqLen > cMaxRawSeqLen || pPE2Seq == NULL)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"AddPESeqs: PE2SeqLen: %d pPE2Seq is %s",PE2SeqLen,pPE2Seq == NULL ? "NULL" : "not-null");
	return(eBSFerrParams);
	}


AcquireSerialise();

MemMinReq =  4 * (2000 + m_Sequences.Seqs2AssembOfs + ((PE1SeqLen + PE2SeqLen + 28) / 15)); // 2K additional for headers plus a little spare

// ensure sufficent memory has been allocated to hold this sequence...
if(MemMinReq > m_Sequences.AllocMemSeqs2Assemb)		
	{
	// try to realloc memory with a 20% increase
	size_t memreq;
	void *pAllocd;
	memreq = (size_t)((MemMinReq * 120) / (UINT64)100);
#ifdef _WIN32
	pAllocd = realloc(m_Sequences.pSeqs2Assemb,memreq);
#else
	pAllocd = mremap(m_Sequences.pSeqs2Assemb,m_Sequences.AllocMemSeqs2Assemb,memreq,MREMAP_MAYMOVE);
	if(pAllocd == MAP_FAILED)
		pAllocd = NULL;
#endif
	if(pAllocd == NULL)
		{
		ReleaseSerialise();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeq: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}

	m_Sequences.pSeqs2Assemb = pAllocd;
	m_Sequences.AllocMemSeqs2Assemb = (UINT64)memreq;
	}

tSeqWrd4 *pSeqWrds;
tSeqWrd4 *pPackSeq;


pSeqWrds = (tSeqWrd4 *)m_Sequences.pSeqs2Assemb;
pPackSeq = &pSeqWrds[m_Sequences.Seqs2AssembOfs];
if(m_Sequences.Seqs2AssembOfs == 0)			// first sequence, so preceed with a cSeqWrd1BOS
	{
	*pPackSeq++ = cSeqWrd4BOS;
	m_Sequences.Seqs2AssembOfs += 1;
	}
m_Sequences.NumSeqs2Assemb += 1;
pPackSeq = (tSeqWrd4 *)SetSeqHeader(pPackSeq,m_Sequences.NumSeqs2Assemb,SrcFileID,PE1SeqFlags,PE1SeqLen,NULL);
m_Sequences.Seqs2AssembOfs += GenPackedSeqWrds(PE1SeqLen,pPE1Seq,0,pPackSeq);

pSeqWrds = (tSeqWrd4 *)m_Sequences.pSeqs2Assemb; // reload in case was realloc'd by some other thread
pPackSeq = &pSeqWrds[m_Sequences.Seqs2AssembOfs];
m_Sequences.NumSeqs2Assemb += 1;
pPackSeq = (tSeqWrd4 *)SetSeqHeader(pPackSeq,m_Sequences.NumSeqs2Assemb,SrcFileID,PE2SeqFlags,PE1SeqLen,NULL);
m_Sequences.Seqs2AssembOfs += GenPackedSeqWrds(PE1SeqLen,pPE1Seq,0,pPackSeq);

pSeqWrds[m_Sequences.Seqs2AssembOfs] = cSeqWrd4EOS;			// currently last sequence, so finish with a cSeqWrd4EOS

m_Sequences.Seqs2AssembLen += (PE1SeqLen + PE2SeqLen);
m_Sequences.bSeqs2AssembDirty;
ReleaseSerialise();
return(eBSFSuccess);
}



int
CKangadna::GenPackedSeqWrds(int SeqLen,				// sequence length ptd at by pSeq
				etSeqBase *pSeq,		// pts to sequence requiring packing
				int MaxSeqWrds,			// truncate packed sequence if requiring more than this number of packed words, 0 if no limit
				void *pDstSeq)			// write generated packed sequence words starting at this location
{
int Rslt;
Rslt = GenPackedSeqWrds(SeqLen,pSeq,MaxSeqWrds,(tSeqWrd4 *)pDstSeq);
return(Rslt);
}


int										// sequence was packed into this number of tSeqWrd4s
CKangadna::GenPackedSeqWrds(int SeqLen,	// sequence length ptd at by pSeq
				etSeqBase *pSeq,		// pts to sequence requiring packing
				int MaxSeqWrds,			// truncate packed sequence if requiring more than this number of packed words, 0 if no limit
				tSeqWrd4 *pDstSeq)		// write generated packed sequence words starting at this location
{
int Idx;
UINT8 Base;
tSeqWrd4 PackedSeqWord;
int NumSeqWords;
int NumPackedBases;
int NumBases;
int NumSeqWrds;

NumSeqWords = (SeqLen+14)/15;			// requires this many 32bit words to contain packed sequence; each word contains max of 15 bases
NumPackedBases = NumSeqWords * 16;      // total packed bases including filler bases and word type 2MSBs
PackedSeqWord = 0;
NumSeqWrds = 0;
NumBases = 0;
for(Idx = 0; Idx < NumPackedBases; Idx++)
	{
	if((Idx & 0x0f) == 0x0f)
		{
		if(NumBases > SeqLen)
			PackedSeqWord |= cSeqWrd4PartSeq;
		*pDstSeq++ = PackedSeqWord;
		PackedSeqWord = 0;
		NumSeqWrds += 1;
		if(MaxSeqWrds && NumSeqWrds == MaxSeqWrds)
			break;
		continue;
		}
	PackedSeqWord <<= 2;
	if(NumBases < SeqLen)
		{
		Base = *pSeq++ & 0X07;
		if(Base > eBaseT)
			Base = rand() & 0x03;
		PackedSeqWord |= Base;
		}
	else
		if(NumBases > SeqLen)
			PackedSeqWord |= 0x01;
	NumBases += 1;
	}
*pDstSeq = cSeqWrd4EOS;
return(NumSeqWrds);
}

// Every packed sequence is immediately preceded by 3 x 32bit (tSeqWrd4's) header words
// Header words have their most significant bit 31 set whereas the packed sequence words have bit 31 reset so header words can easily be distingished from sequence words
// The first header word contains most significant 2 bits of a 32bit sequence identifier (in bits 24 and 25), a 8bit source file identifier (in bits 16..23), 16bit flags (in bits 0..15) 
// The second header word contains the low 30bits of the 32bit sequence identifier (in bits 0..29)
// The third header word contains the 30bit sequence length (in bits 0..29)

tSeqWrd4 *									// returned ptr to 1st word of actual packed sequence which is immediately following the header
CKangadna::GetSeqHeader(tSeqWrd4 *pSeqWrd,	// pts to a SeqWrd within the sequence
			tSeqID *pSeqID,				// returned 32 bit sequence identifier
			UINT32 *pSrcFileID,			// returned 8 bit source file identifier
			UINT32 *pFlgs,				// returned 16 bit sequence flags
			UINT32 *pSeqLen,			// returned 30 bit sequence length
			bool bSerialise)			// set true if access to headers are required to be serialised
{
UINT32 SeqID;

if(pSeqID != NULL)
	*pSeqID = 0;
if(pSrcFileID != NULL)
	*pSrcFileID = 0;
if(pFlgs != NULL)
	*pFlgs = 0;
if(pSeqLen != NULL)
	*pSeqLen = 0;

if(pSeqWrd == NULL)
	return(NULL);

tSeqWrd4 *pSeqHdrWrd;
tSeqWrd4 SeqWrd;
pSeqHdrWrd = (tSeqWrd4 *)pSeqWrd;
if(*pSeqHdrWrd == cSeqWrd4BOS || *pSeqHdrWrd == cSeqWrd4EOS) // ensure not trying to get header for BOS or EOS
	return(NULL);

while((*pSeqHdrWrd & cSeqWrd4LSWHdr) != cSeqWrd4MSWHdr)	// backup until the most significant header word
	pSeqHdrWrd -= 1;
if(bSerialise)
	AcquireSerialiseSeqHdr();

if(pFlgs != NULL || pSrcFileID != NULL || pSeqID != NULL)	// if caller requires flags or file identifier or sequence identifier
	{
	SeqWrd = *pSeqHdrWrd++;
	if(pFlgs != NULL)
		*pFlgs = (UINT32)SeqWrd & 0x000ffff;			// 16 flag bits
	if(pSrcFileID != NULL)							
		*pSrcFileID = ((UINT32)SeqWrd >> 16) & 0x00ff;	// 8bit file identifier
	if(pSeqID != NULL)						
		{
		SeqID = ((UINT32)SeqWrd << 6) & 0x0c0000000;	// bits 24,25 into bits 30,31 of identifier 
		*pSeqID = SeqID | *pSeqHdrWrd & 0x3fffffff;		// bits 0..29 into bits 0..29 of identifier
		}
	pSeqHdrWrd += 1;
	}
else
	pSeqHdrWrd += 2;

if(pSeqLen != NULL)								
	*pSeqLen = *pSeqHdrWrd &  0x3fffffff;		// 30bit sequence length
pSeqHdrWrd += 1;

if(bSerialise)
	ReleaseSerialiseSeqHdr();

return(pSeqHdrWrd);
}



tSeqWrd4 *								 // returned ptr to next word to write packed sequence 
CKangadna::SetSeqHeader(tSeqWrd4 *pSeqHdrWrd, // pts to where to write the sequence header
			tSeqID SeqID,				// sequence identifier (32 bits)
			UINT32 SrcFileID,			// identifies source file from which this sequence was processed (8bit identifier)
			UINT32 Flgs,				// sequence flags (16 bits)
			UINT32 SeqLen,				// sequence length (30 bits)
			int *pNumSeqWrds,			// used to return the number of SeqWrds used
			bool bAddingHdr,			// set true if adding header, false if updating an existing header
			bool bSerialise)			// set true if access to headers are serialised
{
UINT64 NumSeqWrds;
UINT64 AvailSeqWrds;

if(pNumSeqWrds != NULL)
	*pNumSeqWrds = 0;
if(pSeqHdrWrd == NULL || SeqID == 0 || SeqLen == 0 || SeqLen > cMaxContigLen || SrcFileID == 0 || SrcFileID > 0xff || Flgs & 0xffff0000)
	return(NULL);

if(bSerialise)
	AcquireSerialiseSeqHdr();

// ensure that sufficent memory has been allocated to hold this sequence header plus the packed sequence plus a few spare for an appended EOS word
NumSeqWrds = 1 + 3 + ((SeqLen + 14) / 15);					// header requires 3 SeqWrds plus packing 15 bases per SeqWrd plus potential EOS SeqWrd
AvailSeqWrds = (((m_Sequences.AllocMemSeqs2Assemb - (m_Sequences.Seqs2AssembOfs * 4)) + 14) / 15);

if(bAddingHdr && NumSeqWrds > AvailSeqWrds)
	{
	if(bSerialise)
		ReleaseSerialiseSeqHdr();
	return(NULL);
	}	
		
*pSeqHdrWrd++ = Flgs | (SrcFileID << 16) | ((SeqID >> 6) & 0x03000000) | cSeqWrd4MSWHdr;	    // this is the most significant word
*pSeqHdrWrd++ = (SeqID & 0x3fffffff) | cSeqWrd4LSWHdr;		// another least significant word
*pSeqHdrWrd++ = (SeqLen & 0x3fffffff)| cSeqWrd4LSWHdr;		// last least significant word, note that bits 5..0 currently unassigned
NumSeqWrds = 3;

if(pNumSeqWrds != NULL)
	*pNumSeqWrds = (int)NumSeqWrds;
if(bAddingHdr)
	m_Sequences.Seqs2AssembOfs += NumSeqWrds;
if(bSerialise)
	ReleaseSerialiseSeqHdr();
return(pSeqHdrWrd);
}

int
CKangadna::PackedRevCpl(tSeqID SeqID) // Will inplace reverse complement the specified sequence
{
tSeqWrd4 *pSeq;
if((pSeq = (tSeqWrd4 *)GetSeqHeader(SeqID,NULL,NULL,NULL,false))==NULL)
	return(eBSFerrNoEntries);	
return(PackedRevCpl(pSeq,false,0));
}


// PackedRevCplPE
// Exchanges the PE1 and PE2 sequences, RevCpl these sequences, and then updates the headers
//  
int
CKangadna::PackedRevCplPE(tSeqID PE1SeqID) // Will firstly exchange both PE1 and PE2 sequences with their headers followed by RevCpl
{
tSeqWrd4 TmpWrd;
UINT32 PE1Flags;
UINT32 PE1SeqLen;
UINT32 PE2Flags;
UINT32 PE2SeqLen;
size_t NumPE1SeqWrds;
size_t NumPE2SeqWrds;
tSeqWrd4 *pSeqPE1;
tSeqWrd4 *pSeqPE2;

if(m_Sequences.pTmpRevCplSeqs == NULL)
	{
	if((m_Sequences.pTmpRevCplSeqs = new tSeqWrd4[cMaxOvrlapSeqWrds])==NULL)
		return(eBSFerrMem);
	}

if((pSeqPE1 = (tSeqWrd4 *)GetSeqHeader(PE1SeqID,NULL,&PE1Flags,&PE1SeqLen,false))==NULL)
	return(eBSFerrNoEntries);	
	// ensure it is a PE1
if(!(PE1Flags & cFlgSeqPE) || (PE1Flags & cFlgSeqPE2))	
	return(eBSFerrNoEntries);

if((pSeqPE2 = (tSeqWrd4 *)GetSeqHeader(PE1SeqID+1,NULL,&PE2Flags,&PE2SeqLen,false))==NULL)
	return(eBSFerrNoEntries);	
	// ensure it is a PE2
if(!(PE2Flags & cFlgSeqPE2))	
	return(eBSFerrNoEntries);

pSeqPE1 -= 3;
pSeqPE2 -= 3;

	// number of tSeqWrd4's including 3 for the header in PE1 and PE2
NumPE1SeqWrds = 3 + ((PE1SeqLen + 14) / 15);
NumPE2SeqWrds = 3 + ((PE2SeqLen + 14) / 15);
m_Sequences.pSeqStarts[PE1SeqID] = m_Sequences.pSeqStarts[PE1SeqID-1] + NumPE2SeqWrds;

	// before exchanging the headers and sequences then exchange flags and sequence identifiers so
	// after headers exchanged the sequence identifiers are ascending with PE1 identifier followed by PE2 identifier
TmpWrd = *pSeqPE1;		
*pSeqPE1 = *pSeqPE2;
*pSeqPE2 = TmpWrd;
TmpWrd = pSeqPE1[1];
pSeqPE1[1]=pSeqPE2[1];
pSeqPE2[1] = TmpWrd;

	// now exchange PE1 and PE2 - because of potential overlap copy issues with memcpy, memmove is used
memcpy(m_Sequences.pTmpRevCplSeqs,pSeqPE1,sizeof(tSeqWrd4) * NumPE1SeqWrds);
memmove(pSeqPE1,pSeqPE2,sizeof(tSeqWrd4) * NumPE2SeqWrds);
memmove(pSeqPE1 + NumPE2SeqWrds,m_Sequences.pTmpRevCplSeqs,sizeof(tSeqWrd4) * NumPE1SeqWrds);

PackedRevCpl(pSeqPE1 + 3);
pSeqPE2 = pSeqPE1 + NumPE2SeqWrds + 6;
PackedRevCpl(pSeqPE2);

return(eBSFSuccess);
}


// PackedRevCplAllIncPEs
int
CKangadna::PackedRevCplAllIncPEs(void) // Will firstly exchange all PE1 and PE2 sequences with their headers followed by RevCpl all sequences
{
int Rslt;
tSeqWrd4 *pSeqToRev;
tSeqID SeqID;
tSeqWrd4 TmpWrd;
UINT32 PE1Flags;
UINT32 PE1SeqLen;
UINT32 PE2Flags;
UINT32 PE2SeqLen;
size_t NumPE1SeqWrds;
size_t NumPE2SeqWrds;
tSeqWrd4 *pSeqPE1;
tSeqWrd4 *pSeqPE2;

if(m_Sequences.NumSeqs2Assemb == 0)
	return(eBSFerrNoEntries);

if(m_Sequences.pTmpRevCplSeqs == NULL)
	{
	if((m_Sequences.pTmpRevCplSeqs = new tSeqWrd4[cMaxOvrlapSeqWrds])==NULL)
		return(eBSFerrMem);
	}

// firstly locate and exchange PE1 and PE2 sequences
pSeqToRev = NULL;
while((pSeqToRev = (tSeqWrd4 *)IterSeqHeaders(pSeqToRev,&SeqID,NULL,&PE1Flags,&PE1SeqLen,false)) != NULL)
	{
	if(!(PE1Flags & cFlgSeqPE))			// skip onto next sequence until a PE
		continue;

	pSeqPE1 = pSeqToRev - 3;			// pt to initial header word

	// have a PE1, advance onto PE2 which is expected to immediately follow the PE1
	pSeqToRev = (tSeqWrd4 *)IterSeqHeaders(pSeqToRev,&SeqID,NULL,&PE2Flags,&PE2SeqLen,false);

	pSeqPE2 = pSeqToRev - 3;			// pt to initial header word

	// before exchanging the headers and sequences then exchange flags and sequence identifiers so
	// after headers exchanged the sequence identifiers are ascending with PE1 identifier followed by PE2 identifier
	TmpWrd = *pSeqPE1;		
	*pSeqPE1 = *pSeqPE2;
	*pSeqPE2 = TmpWrd;
	TmpWrd = pSeqPE1[1];
	pSeqPE1[1]=pSeqPE2[1];
	pSeqPE2[1] = TmpWrd;

	// number of tSeqWrd4's including 3 for the header in PE1 and PE2
	NumPE1SeqWrds = 3 + ((PE1SeqLen + 14) / 15);
	NumPE2SeqWrds = 3 + ((PE2SeqLen + 14) / 15);

	// now exchange PE1 and PE2 - because of potential overlap copy issues with memcpy, memmove is used
	memcpy(m_Sequences.pTmpRevCplSeqs,pSeqPE1,sizeof(tSeqWrd4) * NumPE1SeqWrds);
	memmove(pSeqPE1,pSeqPE2,sizeof(tSeqWrd4) * NumPE2SeqWrds);
	memmove(pSeqPE1 + NumPE2SeqWrds,m_Sequences.pTmpRevCplSeqs,sizeof(tSeqWrd4) * NumPE1SeqWrds);

	pSeqToRev = pSeqPE1 + NumPE2SeqWrds + 3;
	}

// now RevCpl all sequences
Rslt = PackedRevCplAll(false,false);
return(Rslt);
}


int
CKangadna::PackedRevCplAll(bool bRevOnly,		// if true then reverse only, false to reverse complement
						bool bPE2Only)			// if true then reverse complement PE2 only
{
int Rslt;
UINT32 Flags;
tSeqWrd4 *pSeqToRev;
UINT32 HdrSeqLen;

pSeqToRev = NULL;
while((pSeqToRev = IterSeqHeaders(pSeqToRev,NULL,NULL,&Flags,&HdrSeqLen,false)) != NULL)
	{
	if(bPE2Only && !(Flags & cFlgSeqPE2))
		continue;

	Rslt = PackedRevCpl(pSeqToRev,bRevOnly,0);
	if(Rslt <= 0)
		break;
	}

return(Rslt);
}


// inplace reverse complement
// Assumes sequence has been terminated by sequence terminator EOS or following sequence header word with MSB set
// The maximum number of SeqWrds to reverse can also be specified
int
CKangadna::PackedRevCpl(tSeqWrd4 *pSeqToRev,
			bool bRevOnly,				// if true then reverse only, false to reverse complement
			UINT32 MaxSeqWrds)			// reverse for at most this many SeqWrds (0 for no limit)
{
int Idx;
int NumRevd;
int Bases2Xchg;				// number of bases to be exchanged
tSeqWrd4 *pLSeqWrd;
tSeqWrd4 *pRSeqWrd;
tSeqWrd4 SeqWrd;
UINT32 XBase;
UINT32 LBase;
UINT32 RBase;
UINT32 LSeqWrdMsk;
UINT32 RSeqWrdMsk;
int LShf;
int RShf;

if(pSeqToRev == NULL)
	return(-1);
			
if(*pSeqToRev & cSeqWrd4MSWHdr)
	return(0);

if(MaxSeqWrds == 0)
	MaxSeqWrds = 0xffffffff;

// locate the terminating right word
pLSeqWrd = pSeqToRev;
pRSeqWrd = pSeqToRev;
RSeqWrdMsk = 0x03;
RShf = 0;
NumRevd = 0;
while(MaxSeqWrds-- && !((SeqWrd = *pRSeqWrd) & cSeqWrd4MSWHdr))
	{
	if(SeqWrd & cSeqWrd4PartSeq)
		{
		RSeqWrdMsk = 0x03;
		RShf = 2;
		for(Idx = 0; Idx < 14; Idx++, RShf += 2)
			{
			if(!(SeqWrd & RSeqWrdMsk))
				break;
			RSeqWrdMsk <<= 2;
			}
		NumRevd += 14 - Idx;
		RSeqWrdMsk <<= 2;
		pRSeqWrd += 1;
		break;
		}
	NumRevd += 15;
	pRSeqWrd += 1;
	}
if(NumRevd == 1)
	return(1);

pRSeqWrd -= 1;

Bases2Xchg = NumRevd/2;
LSeqWrdMsk = 0x30000000;				// left starts with this mask
LShf = 28;
for(Idx = 0; Idx < Bases2Xchg; Idx++)
	{
	LBase = 0x03 & ((*pLSeqWrd & LSeqWrdMsk) >> LShf); // extract bases to be exchanged
	RBase = 0x03 & ((*pRSeqWrd & RSeqWrdMsk) >> RShf);
	if(!bRevOnly)							// if not reversing only then complement the bases 
		{
		LBase = ~LBase & 0x03;
		RBase = ~RBase & 0x03;
		}
	XBase = LBase;						// now exchange the bases
	LBase = RBase;
	RBase = XBase;
	*pLSeqWrd &= ~LSeqWrdMsk;			// merge the bases back
	*pLSeqWrd |= LBase << LShf;
	*pRSeqWrd &= ~RSeqWrdMsk;
	*pRSeqWrd |= RBase << RShf;

	LSeqWrdMsk >>= 2;
	LShf -= 2;
	if(LSeqWrdMsk == 0)
		{
		pLSeqWrd += 1;
		LSeqWrdMsk = 0x30000000;
		LShf = 28;
		}

	RSeqWrdMsk <<= 2;
	RShf += 2;
	if(RSeqWrdMsk == 0xc0000000)
		{
		pRSeqWrd -= 1;
		RSeqWrdMsk = 0x00000003;
		RShf = 0;
		}
	}

if((NumRevd & 0x01) && !bRevOnly)		// if odd number of bases then the central base may require complementing
	{
	LBase = (*pLSeqWrd & LSeqWrdMsk) >> LShf; // extract base to complement
	LBase = ~LBase & 0x03;
	*pLSeqWrd &= ~LSeqWrdMsk;			// merge the bases back
	*pLSeqWrd |= LBase << LShf;
	}

return(NumRevd);
}

// inplace packed sequence complement
// Assumes sequence has been terminated by sequence terminator EOS or following sequence header word with MSB set
// The maximum number of SeqWrds to complement can also be specified
int										// returned number of bases complemented
CKangadna::PackedCpl(tSeqWrd4 *pSeqToCpl,
			UINT32 MaxSeqWrds)			// complement for at most this many SeqWrds (0 for no limit)
{
int Idx;
tSeqWrd4 SeqWrd;
int NumCpld;
UINT32 Msk;
UINT32 Filler;
UINT32 FillerMsk;

if(pSeqToCpl == NULL)
	return(-1);
if(*pSeqToCpl & cSeqWrd4MSWHdr)
	return(0);
if(MaxSeqWrds == 0)
	MaxSeqWrds = 0x0ffffffff;
NumCpld = 0;
while(MaxSeqWrds-- && !((SeqWrd = *pSeqToCpl) & cSeqWrd4MSWHdr))
	{
	if(!(SeqWrd & cSeqWrd4PartSeq))
		{
		SeqWrd = ~SeqWrd;			// complements all bases in word!
		SeqWrd &= cSeqWrd4Msk;		// need to reset 2 MSBits
		NumCpld += 15;
		}
	else	// else it's a partial with 1..15 bases
		{
		SeqWrd = ~SeqWrd;			// complements all bases in word!
		Msk = 0x00000003;
		Filler = 0x00000000;
		FillerMsk = 0xfffffffc;
		for(Idx = 0; Idx < 14; Idx++)
			{
			if((SeqWrd & Msk) == Msk) // check for and reset end of partial sequence indicator
				{
				SeqWrd &= FillerMsk;
				SeqWrd |= Filler; 
				NumCpld += (14 - Idx);
				break;
				}
			Msk <<= 2;
			FillerMsk <<= 2;
			Filler <<= 2;
			Filler |= 0x02;
			}
		SeqWrd &= cSeqWrd4Msk;
    	SeqWrd |= cSeqWrd4PartSeq;
		}
	*pSeqToCpl++ = SeqWrd;
	}
return(NumCpld);
}


UINT32				// returned number of sequences in m_Sequences.pSeqs2Assemb
CKangadna::ValidateSeqs2AssembStarts(UINT32 *pNumPEs,
							  UINT32 *pNumSEs)
{
UINT32 NumSfxEls;
UINT32 NumSeqStarts;
tSeqWrd4 SeqWord;
tSeqWrd4 *pSeqWord;
UINT64 SeqWrdIdx;
UINT32 Flags;
UINT32 NumPEs;
UINT32 NumSEs;

tSeqID SeqID;
UINT32 SrcFileID;
UINT32 SeqLen;

if(pNumPEs != NULL)
	*pNumPEs = 0;
if(pNumSEs != NULL)
	*pNumSEs = 0;

if(m_Sequences.pSeqs2Assemb == NULL || m_Sequences.Seqs2AssembOfs <= 1)
	return(0);

pSeqWord = (tSeqWrd4 *)m_Sequences.pSeqs2Assemb;
SeqWord = *pSeqWord++;
if(SeqWord != cSeqWrd4BOS)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: Problem at word %llu, expected %u instead read %u ...",0,SeqWord,cSeqWrd4BOS);
	return(0);
	}

NumSfxEls = 0;
NumSeqStarts = 0;
NumPEs = 0;
NumSEs = 0;
int State = 1;		// 1 if expecting cSeqWrd4MSWHdr
int ExpectType = 0; // 0 if expecting SE or PE1, 1 if expecting PE2
bool bNoTypes = false;
UINT32 PrevFlags = 0;
UINT32 PrevSeqLen = 0;
tSeqID PrevSeqID = 0;

for(SeqWrdIdx = 1; SeqWrdIdx <= m_Sequences.Seqs2AssembOfs; SeqWrdIdx++)
	{
	SeqWord = *pSeqWord++;
	if(SeqWord == cSeqWrd4EOS)
		{
		if(SeqWrdIdx != m_Sequences.Seqs2AssembOfs)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: Problem at word %llu, unexpected %u ...",SeqWrdIdx,cSeqWrd4EOS);
			return(0);
			}
		continue;
		}

	switch(State) {
		case 1:
			if((SeqWord & 0x0e0000000) == cSeqWrd4MSWHdr)
				{
				GetSeqHeader(pSeqWord,			// pts to a SeqWrd within the sequence
						&SeqID,					// returned 32 bit sequence identifier
						&SrcFileID,				// returned 8 bit source file identifier
						&Flags,					// returned 16 bit sequence flags 
						&SeqLen,				// returned 30 bit sequence length
						false);			// set true if access to headers are required to be serialised

				if(SeqID != PrevSeqID + 1)
					gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: Problem at SeqID: %u, word %llu, expected %u for the SeqID...",SeqID,SeqWrdIdx,PrevSeqID+1);

				if(!bNoTypes)
					{
					switch(Flags & 0x03) {
						case 0:				// 0 if SE
							if(ExpectType == 0)
								break;
							gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: Problem at SeqID: %u, word %llu, expected PE2 instead processed SE ...",SeqID,SeqWrdIdx);
							bNoTypes = true;
							break;
						case 1:				// 1 if PE R1
							if(ExpectType == 0)
								{
								ExpectType = 1; // will be expecting the PE2 next
								break;
								}
							gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: Problem at SeqID: %u, word %llu, expected PE2 instead processed PE1 ...",SeqID,SeqWrdIdx);
							bNoTypes = true;
							break;

						case 3:				// 3 if PE R2
							if(ExpectType == 1)
								{
								ExpectType = 0; // will be expecting a SE or PE1 next
								break;
								}
							gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: Problem at SeqID: %u, word %llu, expected SE or PE1 instead processed PE2 ...",SeqID,SeqWrdIdx);
							bNoTypes = true;
							break;

						case 2:				// unexpected!
							gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: Problem at SeqID: %u, word %llu, expected sequence type ...",SeqID,SeqWrdIdx);
							bNoTypes = true;
							break;
						}
					}
				PrevFlags = Flags;
				PrevSeqLen = SeqLen;
				PrevSeqID = SeqID;
				State = 2;			// expecting next to be 1st of 2 cSeqWrd4LSWHdr words
				continue;
				}
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: Problem at word %llu, expected %u instead read %u ...",SeqWrdIdx,cSeqWrd4MSWHdr,SeqWord);
			return(0);
		case 2:
			if((SeqWord & 0x0c0000000) == cSeqWrd4LSWHdr)
				{
				State = 3;			// expecting 2nd of 2 cSeqWrd4LSWHdr words
				continue;
				}
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: Problem at word %llu, expected 1st %u instead read %u ...",SeqWrdIdx,cSeqWrd4LSWHdr,SeqWord);
			return(0);
		case 3:
			if((SeqWord & 0x0c0000000) == cSeqWrd4LSWHdr)
				{
				State = 4;			// expecting at least one either full or partial SeqWrd containg packed bases
				continue;
				}
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: Problem at word %llu, expected 2nd %u instead read %u ...",SeqWrdIdx,cSeqWrd4LSWHdr,SeqWord);
			return(0);

		case 4:
			if((SeqWord & 0x0c0000000) == 0 || (SeqWord & 0x0c0000000) == cSeqWrd4PartSeq)
				{
				State = 0;
				NumSeqStarts += 1;
				if(!(Flags & 0x01))
					NumSEs += 1;
				else
					if(Flags & 0x02)
						NumPEs += 1;
				continue;
				}
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: Problem at word %llu, expected 1st base containing word instead read %u ...",SeqWrdIdx,SeqWord);
			return(0);


		case 0:
			if(!(SeqWord & 0x0c0000000))
				continue;

			if((SeqWord & 0x0c0000000) == cSeqWrd4PartSeq)
				{
				State = 1;			// cSeqWrd4MSWHdr next expected
				continue;
				}

			if(SeqWord & 0x80000000)
				{
				if((SeqWord & 0x0e0000000) == cSeqWrd4MSWHdr)
					{
					pSeqWord -= 1;
					SeqWrdIdx -= 1;
					State = 1;			// cSeqWrd4MSWHdr next expected
					continue;
					}
				// any other word is unexpected - where is the cSeqWrd4MSWHdr?
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: Problem at word %llu, expected initial headr word %u instead read %u ...",SeqWrdIdx,cSeqWrd4MSWHdr,SeqWord);
				return(0);
				}
		}
	}
if(!bNoTypes && ExpectType != 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: Problem at last word %llu, was expecting PE2 ...",SeqWrdIdx);
	}

// now to validate the 
UINT32 Flags2;
UINT32 SeqLen2;
UINT32 SrcFileID2;
tSeqID SeqID2;
tSeqID PrevSeqID2;
tSeqWrd4 *pPrevSeqWrd;
tSeqWrd4 *pSeqWrd2;
pPrevSeqWrd = NULL;
for(SeqID = 1; SeqID <= m_Sequences.NumSeqs2Assemb; SeqID++)
	{
	if((pSeqWord = (tSeqWrd4 *)GetSeqHeader(SeqID,&SrcFileID,&Flags,&SeqLen,false))==NULL)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: GetSeqHeader returned NULL for SeqID %u ...",SeqID);
		return(0);
		}
	if(pPrevSeqWrd == pSeqWord)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: GetSeqHeader returned same pSeqWord for SeqID %u and %u ...",PrevSeqID2,SeqID);
		return(0);
		}
	
	if((pSeqWrd2 = GetSeqHeader(pSeqWord+2,			// pts to a SeqWrd within the sequence
						&SeqID2,					// returned 32 bit sequence identifier
						&SrcFileID2,				// returned 8 bit source file identifier
						&Flags2,					// returned 16 bit sequence flags 
						&SeqLen2,				// returned 30 bit sequence length
						false))==NULL)			// set true if access to headers are required to be serialised
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: A GetSeqHeader returned NULL for pSeqWord SeqID ...",SeqID);
		return(0);
		}
	if(SeqID != SeqID2)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidateSeqs2AssembStarts: A GetSeqHeader for SeqID %u -  PrevSeqID %u SeqID2 %u ...",SeqID, PrevSeqID2,SeqID2);
		return(0);
		}

	pPrevSeqWrd = pSeqWord;
	PrevSeqID2 = SeqID;
	
	}

if(pNumPEs != NULL)
	*pNumPEs = NumPEs;
if(pNumSEs != NULL)
	*pNumSEs = NumSEs;
return(NumSeqStarts);
}


int
CKangadna::ValidatePartialSeqsStarts(int *pPartialNumPEs,
							  int *pPartialNumSEs)
{
tSeqWrd4 *pSeqWrd;
tSeqWrd4 SeqWrd;
UINT32 SeqLen;
UINT16 SeqFlags;
UINT32 SeqID;
UINT32 NumPEs;
UINT32 NumSEs;

UINT64 WrdIdx;
UINT64 Seqs2AssembLen;
UINT32 NumSeqs2Assemb;

int NumSeqWrds;

if(pPartialNumPEs != NULL)
	*pPartialNumPEs = 0;
if(pPartialNumSEs != NULL)
	*pPartialNumSEs = 0;


if(m_NumPartialSeqs2Assemb == 0 || m_pPartialSeqs2Assemb == NULL)
	return(0);


pSeqWrd = (tSeqWrd4 *)m_pPartialSeqs2Assemb;
pSeqWrd+=1;
	
NumPEs = 0;
NumSEs = 0;
NumSeqs2Assemb = 0;
Seqs2AssembLen = 0;
NumSeqWrds = 0;
SeqID = 0;
for(WrdIdx = 1; WrdIdx <= m_PartialSeqs2AssembOfs; WrdIdx++, pSeqWrd++,NumSeqWrds--)
	{
	if((SeqWrd = *pSeqWrd) == cSeqWrd4EOS)
		break;

	// check that sequence words do not include any header words
	if(NumSeqWrds)
		{
		if(SeqWrd & cSeqWrd4MSWHdr)	// should never see any header words!
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidatePartialSeqsStarts: Header word %x at WrdIdx %u, SeqID %d ...",SeqWrd,WrdIdx,SeqID);
			return(0);
			}
		if(SeqWrd & 0x40000000 && NumSeqWrds != 1)		// if part word then expect this to be the last word
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidatePartialSeqsStarts: part word %x at WrdIdx %u, SeqID %d ...",SeqWrd,WrdIdx,SeqID);
			return(0);
			}
		}
	
	if(SeqID == 0 || NumSeqWrds == 0)
		{
		SeqID += 1;
		if(SeqWrd & cSeqWrd4MSWHdr)	// should never see any header words!
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"ValidatePartialSeqsStarts: Header word %x at WrdIdx %u instead of expected length, SeqID %d ...",SeqWrd,WrdIdx,SeqID);
			return(0);
			}

		SeqLen = SeqWrd;
		SeqFlags = SeqLen & 0x03;
		if(SeqFlags & cFlgSeqPE2)
			NumPEs += 1;
		else
			if(!(SeqFlags & cFlgSeqPE))
				NumSEs += 1;
		SeqLen >>= 2;
		Seqs2AssembLen += SeqLen;
		NumSeqWrds = 1 + ((SeqLen + 14) / 15);
		}
	}
if(SeqID != m_NumPartialSeqs2Assemb || WrdIdx != m_PartialSeqs2AssembOfs || Seqs2AssembLen != m_LenPartialSeqs2Assemb)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"CalcPartialSeqsStarts: Problem at %d ...",SeqID);
if(pPartialNumPEs != NULL)
	*pPartialNumPEs = NumPEs;
if(pPartialNumSEs != NULL)
	*pPartialNumSEs = NumSEs;

return(SeqID);
}

// returns min/mean/max lengths for PE1/PE2 and SE sequences
teBSFrsltCodes	
CKangadna::GetSeqLenDist(UINT32 *pNumPEs,	// total number of paired ends
					  int *pPE1min,	// returned PE1 min length
					  int *pPE1mean,// returned PE1 mean length
					  int *pPE1max, // returned PE1 max length
					  int *pPE2min,	// returned PE2 min length
					  int *pPE2mean,// returned PE2 mean length
					  int *pPE2max, // returned PE2 max length
					  UINT32 *pNumSEs,	// total number of single ends
					  int *pSEmin,	// returned SE min length
					  int *pSEmean,	// returned SE mean length
					  int *pSEmax)	// returned SE max length
{
tSeqID SeqID;
tSeqWrd4 *pPE1SeqWrd;
tSeqWrd4 *pPE2SeqWrd;
UINT32 PE1Flags;
UINT32 PE1SeqLen;
UINT32 PE2Flags;
UINT32 PE2SeqLen;

UINT32 TotNumSESeqs;
UINT32 TotNumPESeqs;

UINT64 SumSESeqLen;
UINT64 SumPE1SeqLen;
UINT64 SumPE2SeqLen;

int MinSELen;
int MinPE1Len;
int MinPE2Len;

int MaxSELen;
int MaxPE1Len;
int MaxPE2Len;

if(pNumPEs != NULL)
	*pNumPEs = 0;
if(pPE1min != NULL)
	*pPE1min = 0;
if(pPE1mean != NULL)
	*pPE1mean = 0;
if(pPE1max != NULL)
	*pPE1max = 0;
if(pPE2min != NULL)
	*pPE2min = 0;
if(pPE2mean != NULL)
	*pPE2mean = 0;
if(pPE2max != NULL)
	*pPE2max = 0;
if(pNumSEs != NULL)
	*pNumSEs = 0;
if(pSEmin != NULL)
	*pSEmin = 0;
if(pSEmean != NULL)
	*pSEmean = 0;
if(pSEmax != NULL)
	*pSEmax = 0;

if(m_Sequences.NumSeqs2Assemb == 0)		// can't generate length stats if nothing to assemble!
	return(eBSFSuccess);

TotNumSESeqs = 0;
TotNumPESeqs = 0;
SumSESeqLen = 0;
SumPE1SeqLen = 0;
SumPE2SeqLen = 0;
MinSELen = 0;
MinPE1Len = 0;
MinPE2Len = 0;
MaxSELen = 0;
MaxPE1Len = 0;
MaxPE2Len = 0;

for(SeqID = 1; SeqID <= m_Sequences.NumSeqs2Assemb; SeqID++)
	{
	if((pPE1SeqWrd = (tSeqWrd4 *)GetSeqHeader(SeqID,NULL,&PE1Flags,&PE1SeqLen,false))==NULL)
		return(eBSFerrObj);
	if(PE1Flags & cFlgSeqPE)
		{
		if(PE1Flags & cFlgSeqPE2)	// shouldn't be the PE2!
			return(eBSFerrObj);
		SeqID += 1;
		if((pPE2SeqWrd = (tSeqWrd4 *)GetSeqHeader(SeqID,NULL,&PE2Flags,&PE2SeqLen,false))==NULL)
			return(eBSFerrObj);
		if(!(PE2Flags & cFlgSeqPE2))	// should be the PE2!
			return(eBSFerrObj);
		TotNumPESeqs += 1;
		if(MinPE1Len == 0 || (int)PE1SeqLen < MinPE1Len)
			MinPE1Len = PE1SeqLen;
		if((int)PE1SeqLen > MaxPE1Len)
			MaxPE1Len = (int)PE1SeqLen;
		SumPE1SeqLen += PE1SeqLen;

		if(MinPE2Len == 0 || (int)PE2SeqLen < MinPE2Len)
			MinPE2Len = PE2SeqLen;
		if((int)PE2SeqLen > MaxPE2Len)
			MaxPE2Len = (int)PE2SeqLen;
		SumPE2SeqLen += PE2SeqLen;
		}
	else    // else must be a single end sequence
		{
		TotNumSESeqs += 1;
		if(MinSELen == 0 || (int)PE1SeqLen < MinSELen)
			MinSELen = (int)PE1SeqLen;
		if((int)PE1SeqLen > MaxSELen)
			MaxSELen = (int)PE1SeqLen;
		SumSESeqLen += PE1SeqLen;
		}
	}

m_Sequences.NumPE2Seqs2Assemb = TotNumPESeqs;
m_Sequences.NumPE1Seqs2Assemb = m_Sequences.NumSeqs2Assemb - m_Sequences.NumPE2Seqs2Assemb;

if(pNumPEs != NULL)
	*pNumPEs = TotNumPESeqs;
if(pPE1min != NULL)
	*pPE1min = MinPE1Len;
if(pPE1mean != NULL && TotNumPESeqs)
	*pPE1mean = (int)(SumPE1SeqLen/TotNumPESeqs);
if(pPE1max != NULL)
	*pPE1max = MaxPE1Len;
if(pPE2min != NULL)
	*pPE2min = MinPE2Len;
if(pPE2mean != NULL && TotNumPESeqs)
	*pPE2mean = (int)(SumPE2SeqLen/TotNumPESeqs);
if(pPE2max != NULL)
	*pPE2max = MaxPE2Len;
if(pNumSEs != NULL)
	*pNumSEs = TotNumSESeqs;
if(pSEmin != NULL)
	*pSEmin = MinSELen;
if(pSEmean != NULL && TotNumSESeqs)
	*pSEmean = (int)(SumSESeqLen/TotNumSESeqs);
if(pSEmax != NULL)
	*pSEmax = MaxSELen;
return(eBSFSuccess);
}

teBSFrsltCodes
CKangadna::SaveAssembSeqs(char *pszFastaFile,	// save assembled sequences - will be saved as pszFastaFile suffxed with 'SE','PE1' and 'PE2' 
						  int PassID,			// if non-zero then append to file names as the pass indicator
						  int LineBreakSeqs)	// default (0) is to line break sequences every 75bp, if this is > 10 then line break every LineBreakSeqs bases
{
int Written;
tSeqWrd4 *pCurSeq;
char *pszAsciiSeq;
tSeqID SeqID;
UINT32 SeqLen;
UINT32 SrcFileID;

char szPassID[10];
char szFastaFileSE[_MAX_PATH];
char szFastaFilePE1[_MAX_PATH];
char szFastaFilePE2[_MAX_PATH];
char *pszFastaSE;
char *pszFastaPE1;
char *pszFastaPE2;

int LineLenSE;
int FastaLenSE;
char *pszFastaSE1;

int LineLenR1;
int FastaLenR1;
char *pszFastaR1;

int LineLenR2;
int FastaLenR2;
char *pszFastaR2;

m_hOutFastaSE = -1;
m_hOutFastaR1 = -1;
m_hOutFastaR2 = - 1;

strcpy(szFastaFileSE,pszFastaFile);
if(PassID > 0)
	{
	sprintf(szPassID,".Pass%d",PassID);
	strcat(szFastaFileSE,szPassID);
	}
strcat(szFastaFileSE,".SE.fasta");

#ifdef _WIN32
if((m_hOutFastaSE = open(szFastaFileSE, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutFastaSE = open(szFastaFileSE, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",szFastaFileSE,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output SE assembly multifasta file created/truncated: '%s'",szFastaFileSE);

if((pszFastaSE = new char [10*cMaxCmpOverlapBases])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory");
	Reset();
	return(eBSFerrMem);
	}
pszFastaSE1 = pszFastaSE;

if(m_Sequences.bPESeqs && (m_Sequences.NumPE1Seqs2Assemb > 0 && m_Sequences.NumPE2Seqs2Assemb > 0))
	{
	strcpy(szFastaFilePE1,pszFastaFile);
	if(PassID)
		strcat(szFastaFilePE1,szPassID);
	strcat(szFastaFilePE1,".R1.fasta");

#ifdef _WIN32
	if((m_hOutFastaR1 = open(szFastaFilePE1, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hOutFastaR1 = open(szFastaFilePE1, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",szFastaFilePE1,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output PE1 assembly multifasta file created/truncated: '%s'",szFastaFilePE1);
	strcpy(szFastaFilePE2,pszFastaFile);
	if(PassID)
		strcat(szFastaFilePE2,szPassID);
	strcat(szFastaFilePE2,".R2.fasta");
#ifdef _WIN32
	if((m_hOutFastaR2 = open(szFastaFilePE2, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hOutFastaR2 = open(szFastaFilePE2, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",szFastaFilePE2,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output PE2 assembly multifasta file created/truncated: '%s'",szFastaFilePE2);

	if((pszFastaPE1 = new char [10*cMaxCmpOverlapBases])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory");
		delete pszFastaSE;
		Reset();
		return(eBSFerrMem);
		}

	if((pszFastaPE2 = new char [10*cMaxCmpOverlapBases])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory");
		delete pszFastaPE1;
		Reset();
		return(eBSFerrMem);
		}
	pszFastaR1 = pszFastaPE1;
	pszFastaR2 = pszFastaPE2;
	pszFastaPE1[0] = '\0';
	pszFastaPE2[0] = '\0';
	}

pCurSeq = NULL;



FastaLenSE = 0;
FastaLenR1 = 0;
FastaLenR2 = 0;

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif

UINT32 NumFastaSESeqs;
UINT32 NumFastaPE1Seqs;
UINT32 NumFastaPE2Seqs;
UINT32 SeqFlags;

NumFastaSESeqs = 0;
NumFastaPE1Seqs = 0;
NumFastaPE2Seqs = 0;

if(LineBreakSeqs == 0)			// if defaulting then line break sequences every 75bp 
	LineBreakSeqs = 75;
else
	{
	if(LineBreakSeqs < 8)		// require at least 8bp per line
		LineBreakSeqs = 8;
	else
		if(LineBreakSeqs > 32000)	// need to be reasonable!
			LineBreakSeqs = 32000;
	}

while((pCurSeq = IterSeqHeaders(pCurSeq, // iterate to next sequence following (NULL to return 1st sequence)
							&SeqID,		// returned sequence identifier
							&SrcFileID,	// returned 8 bit source file identifier
							&SeqFlags,		// returned 16 bit sequence flags
							&SeqLen,false))!=NULL)		// returned 30 bit sequence length
	{
	if(!(SeqFlags & cFlgSeqPE))	// if not a PE2 sequence then must be a SE
		{
		if((FastaLenSE + 1000) >= 10*cMaxCmpOverlapBases)
			{
			Written = write(m_hOutFastaSE,pszFastaSE,FastaLenSE);
			FastaLenSE = 0;
			pszFastaSE1 = pszFastaSE;
			}
		NumFastaSESeqs += 1;
		pszAsciiSeq = AsciifySequence(pCurSeq,SeqLen,false);
		FastaLenSE += sprintf((char *)&pszFastaSE[FastaLenSE],">Seq%u|%d %d|%d|SE\n",NumFastaSESeqs,SrcFileID,SeqLen,SrcFileID);

		pszFastaSE1 = &pszFastaSE[FastaLenSE];
		LineLenSE = 0;

		while(*pszAsciiSeq != '\0')
			{
			if(LineLenSE > LineBreakSeqs)
				{
				*pszFastaSE1++ = '\n';
				FastaLenSE += 1;
				LineLenSE = 0;
				if((FastaLenSE + LineBreakSeqs + 1000) >= 10*cMaxCmpOverlapBases)
					{
					Written = write(m_hOutFastaSE,pszFastaSE,FastaLenSE);
					FastaLenSE = 0;
					pszFastaSE1 = pszFastaSE;
					}
				}
			*pszFastaSE1++ = *pszAsciiSeq++;
			FastaLenSE += 1;
			LineLenSE += 1;
			}

		if(LineLenSE > 0)
			{
			*pszFastaSE1++ = '\n';
			FastaLenSE += 1;
			LineLenSE += 1;
			}
		continue;
		}

	// part of a PE - either PE1 or PE2
	if(!(SeqFlags & cFlgSeqPE2))	// if PE1 sequence
		{		
		if((FastaLenR1 + LineBreakSeqs + 1000) >= 10*cMaxCmpOverlapBases)
			{
			Written = write(m_hOutFastaR1,pszFastaPE1,FastaLenR1);
			FastaLenR1 = 0;
			pszFastaR1 = pszFastaPE1;
			}
		NumFastaPE1Seqs += 1;
		pszAsciiSeq = AsciifySequence(pCurSeq,SeqLen,false);
		FastaLenR1 += sprintf((char *)&pszFastaPE1[FastaLenR1],">Seq%u|%d/1 %d|%d|PE|1\n",NumFastaPE1Seqs,SrcFileID,SeqLen,SrcFileID);

		pszFastaR1 = &pszFastaPE1[FastaLenR1];
		LineLenR1 = 0;

		while(*pszAsciiSeq != '\0')
			{
			if(LineLenR1 > LineBreakSeqs)
				{
				*pszFastaR1++ = '\n';
				FastaLenR1 += 1;
				LineLenR1 = 0;
				if((FastaLenR1 + LineBreakSeqs + 1000) >= 10*cMaxCmpOverlapBases)
					{
					Written = write(m_hOutFastaR1,pszFastaPE1,FastaLenR1);
					FastaLenR1 = 0;
					pszFastaR1 = pszFastaPE1;
					}
				}
			*pszFastaR1++ = *pszAsciiSeq++;
			FastaLenR1 += 1;
			LineLenR1 += 1;
			}

		if(LineLenR1 > 0)
			{
			*pszFastaR1++ = '\n';
			FastaLenR1 += 1;
			LineLenR1 += 1;
			}
		}
	else // else PE2 so need to reverse complement the sequence
		{		
		if((FastaLenR2 + LineBreakSeqs + 1000) >= 10*cMaxCmpOverlapBases)
			{
			Written = write(m_hOutFastaR2,pszFastaPE2,FastaLenR2);
			FastaLenR2 = 0;
			pszFastaR2 = pszFastaPE2;
			}
		
		NumFastaPE2Seqs += 1;
		pszAsciiSeq = AsciifySequence(pCurSeq,SeqLen,true);
		FastaLenR2 += sprintf((char *)&pszFastaPE2[FastaLenR2],">Seq%u|%d/2 %d|%d|PE|2\n",NumFastaPE2Seqs,SrcFileID-1,SeqLen,SrcFileID);
		LineLenR2 = 0;
		pszFastaR2 = &pszFastaPE2[FastaLenR2];
		while(*pszAsciiSeq != '\0')
			{
			if(LineLenR2 > LineBreakSeqs)
				{
				*pszFastaR2++ = '\n';
				FastaLenR2 += 1;
				LineLenR2 = 0;
				if((FastaLenR2 + LineBreakSeqs + 1000) >= 10*cMaxCmpOverlapBases)
					{
					Written = write(m_hOutFastaR2,pszFastaPE2,FastaLenR2);
					FastaLenR2 = 0;
					pszFastaR2 = pszFastaPE2;
					}
				}
			*pszFastaR2++ = *pszAsciiSeq++;
			FastaLenR2 += 1;
			LineLenR2 += 1;
			}
		if(LineLenR2 > 0)
			{
			*pszFastaR2++ = '\n';
			FastaLenR2 += 1;
			LineLenR2 += 1;
			}
		}
	}

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif

if(m_hOutFastaSE != -1)
	{
	if(FastaLenSE)
		Written = write(m_hOutFastaSE,pszFastaSE,FastaLenSE);
	FastaLenSE = 0;
#ifdef _WIN32
	_commit(m_hOutFastaSE);
#else
	fsync(m_hOutFastaSE);
#endif
	close(m_hOutFastaSE);
	m_hOutFastaSE = -1;
	}
delete pszFastaSE;

if(m_Sequences.bPESeqs && (m_Sequences.NumPE1Seqs2Assemb > 0 && m_Sequences.NumPE2Seqs2Assemb > 0))
	{
	if(m_hOutFastaR1 != -1)
		{
		if(FastaLenR1)
			Written = write(m_hOutFastaR1,pszFastaPE1,FastaLenR1);
		FastaLenR1 = 0;
#ifdef _WIN32
		_commit(m_hOutFastaR1);
#else
		fsync(m_hOutFastaR1);
#endif
		close(m_hOutFastaR1);
		m_hOutFastaR1 = -1;
		}
	delete pszFastaPE1;

	if(m_hOutFastaR2 != -1)
		{
		if(FastaLenR2)
			Written = write(m_hOutFastaR2,pszFastaPE2,FastaLenR2);
		FastaLenR2 = 0;
#ifdef _WIN32
		_commit(m_hOutFastaR2);
#else
		fsync(m_hOutFastaR2);
#endif
		close(m_hOutFastaR2);
		m_hOutFastaR2 = -1;
		}
	delete pszFastaPE2;
	}

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
return(eBSFSuccess);
}

teBSFrsltCodes
CKangadna::SaveAssembAsSESeqs(char *pszFastaFile,	// save assembled sequences, PE sequences will be output as single concatenated sequence with 10 Ns separator - will be saved as pszFastaFile suffxed with '.fasta'
						  int LineBreakSeqs)	// default (0) is to line break sequences every 75bp, if this is > 10 then line break every LineBreakSeqs bases
{
int Written;
tSeqWrd4 *pCurSeq;
char *pszAsciiSeq;
tSeqID SeqID;
UINT32 SeqLen;
UINT32 SrcFileID;

char szFastaFileSE[_MAX_PATH];
char *pszFastaSE;

int LineLenSE;
int FastaLenSE;
char *pszFastaSE1;

strcpy(szFastaFileSE,pszFastaFile);
strcat(szFastaFileSE,".fasta");

#ifdef _WIN32
if((m_hOutFastaSE = open(szFastaFileSE, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutFastaSE = open(szFastaFileSE, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",szFastaFileSE,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output assembly multifasta file created/truncated: '%s'",szFastaFileSE);

if((pszFastaSE = new char [10*cMaxCmpOverlapBases])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory");
	Reset();
	return(eBSFerrMem);
	}
pszFastaSE1 = pszFastaSE;

pCurSeq = NULL;

FastaLenSE = 0;

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif

UINT32 NumFastaSESeqs;
UINT32 NumFastaPE1Seqs;
UINT32 NumFastaPE2Seqs;
UINT32 SeqFlags;

NumFastaSESeqs = 0;
NumFastaPE1Seqs = 0;
NumFastaPE2Seqs = 0;

if(LineBreakSeqs == 0)			// if defaulting then line break sequences every 75bp 
	LineBreakSeqs = 75;
else
	{
	if(LineBreakSeqs < 8)		// require at least 8bp per line
		LineBreakSeqs = 8;
	else
		if(LineBreakSeqs > 32000)	// need to be reasonable!
			LineBreakSeqs = 32000;
	}
LineLenSE = 0;
while((pCurSeq = IterSeqHeaders(pCurSeq, // iterate to next sequence following (NULL to return 1st sequence)
							&SeqID,		// returned sequence identifier
							&SrcFileID,	// returned 5 bit source file identifier in bits 4..0
							&SeqFlags,		// returned 16 bit sequence flags in bits 15..0
							&SeqLen, false))!=NULL)		// returned 30 bit sequence length
	{
	if((FastaLenSE + 1000) >= 10*cMaxCmpOverlapBases)
		{
		Written = write(m_hOutFastaSE,pszFastaSE,FastaLenSE);
		FastaLenSE = 0;
		pszFastaSE1 = pszFastaSE;
		}
	if(!(SeqFlags & cFlgSeqPE))	// if not a PE2 sequence then must be a SE
		{
		NumFastaSESeqs += 1;
		pszAsciiSeq = AsciifySequence(pCurSeq,SeqLen,false);
		FastaLenSE += sprintf((char *)&pszFastaSE[FastaLenSE],">Seq%u|SE\n",NumFastaSESeqs);
		pszFastaSE1 = &pszFastaSE[FastaLenSE];
		LineLenSE = 0;
		while(*pszAsciiSeq != '\0')
			{
			if(LineLenSE > LineBreakSeqs)
				{
				*pszFastaSE1++ = '\n';
				FastaLenSE += 1;
				LineLenSE = 0;
				if((FastaLenSE + LineBreakSeqs + 1000) >= 10*cMaxCmpOverlapBases)
					{
					Written = write(m_hOutFastaSE,pszFastaSE,FastaLenSE);
					FastaLenSE = 0;
					pszFastaSE1 = pszFastaSE;
					}
				}
			*pszFastaSE1++ = *pszAsciiSeq++;
			FastaLenSE += 1;
			LineLenSE += 1;
			}

		if(LineLenSE > 0)
			{
			*pszFastaSE1++ = '\n';
			FastaLenSE += 1;
			LineLenSE = 0;
			}
		continue;
		}

	// part of a PE - either PE1 or PE2
	if(!(SeqFlags & cFlgSeqPE2))	// if PE1 sequence
		{		
		NumFastaPE1Seqs += 1;

		pszAsciiSeq = AsciifySequence(pCurSeq,SeqLen,false);

		if(LineLenSE > 0)
			{
			*pszFastaSE1++ = '\n';
			FastaLenSE += 1;
			LineLenSE = 0;
			}

		FastaLenSE += sprintf((char *)&pszFastaSE[FastaLenSE],">Seq%u|PE\n",NumFastaSESeqs);

		pszFastaSE1 = &pszFastaSE[FastaLenSE];
		LineLenSE = 0;

		while(*pszAsciiSeq != '\0')
			{
			if(LineLenSE > LineBreakSeqs)
				{
				*pszFastaSE1++ = '\n';
				FastaLenSE += 1;
				LineLenSE = 0;
				if((FastaLenSE + LineBreakSeqs + 1000) >= 10*cMaxCmpOverlapBases)
					{
					Written = write(m_hOutFastaSE,pszFastaSE,FastaLenSE);
					FastaLenSE = 0;
					pszFastaSE1 = pszFastaSE;
					}
				}
			*pszFastaSE1++ = *pszAsciiSeq++;
			FastaLenSE += 1;
			LineLenSE += 1;
			}
		continue;
		}

		// must be PE2 of a PE
	pszAsciiSeq = AsciifySequence(pCurSeq,SeqLen,false);
	pszFastaSE1 = &pszFastaSE[FastaLenSE];
	for(int Idx = 0; Idx < 10; Idx++)
		{
		if(LineLenSE > LineBreakSeqs)
			{
			*pszFastaSE1++ = '\n';
			FastaLenSE += 1;
			LineLenSE = 0;
			if((FastaLenSE + LineBreakSeqs + 1000) >= 10*cMaxCmpOverlapBases)
				{
				Written = write(m_hOutFastaSE,pszFastaSE,FastaLenSE);
				FastaLenSE = 0;
				pszFastaSE1 = pszFastaSE;
				}
			}
		*pszFastaSE1++ = 'N';
		FastaLenSE += 1;
		LineLenSE += 1;
		}

	while(*pszAsciiSeq != '\0')
		{
		if(LineLenSE > LineBreakSeqs)
			{
			*pszFastaSE1++ = '\n';
			FastaLenSE += 1;
			LineLenSE = 0;
			if((FastaLenSE + LineBreakSeqs + 1000) >= 10*cMaxCmpOverlapBases)
				{
				Written = write(m_hOutFastaSE,pszFastaSE,FastaLenSE);
				FastaLenSE = 0;
				pszFastaSE1 = pszFastaSE;
				}
			}
		*pszFastaSE1++ = *pszAsciiSeq++;
		FastaLenSE += 1;
		LineLenSE += 1;
		}
	if(LineLenSE > 0)
		{
		*pszFastaSE1++ = '\n';
		FastaLenSE += 1;
		LineLenSE = 0;
		}
	continue;
	}

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif

if(m_hOutFastaSE != -1)
	{
	if(FastaLenSE)
		Written = write(m_hOutFastaSE,pszFastaSE,FastaLenSE);
	FastaLenSE = 0;
#ifdef _WIN32
	_commit(m_hOutFastaSE);
#else
	fsync(m_hOutFastaSE);
#endif
	close(m_hOutFastaSE);
	m_hOutFastaSE = -1;
	}
delete pszFastaSE;

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
return(eBSFSuccess);
}


teBSFrsltCodes
CKangadna::GenSeqStarts(bool bGenFlags,			// optionally also generate flags array
						bool bSerialise)		// true if serialisation required
{
size_t ReqAllocMem;

UINT32 SeqID;
UINT64 SeqWrdIdx;
UINT64 *pSeqStarts;
UINT16 *pSeqFlags;

gDiagnostics.DiagOut(eDLDiag,gszProcName,"GenSeqStarts: Generating sequence start offsets");

if(bSerialise)
	AcquireSerialise();

ReqAllocMem = (10 + m_Sequences.NumSeqs2Assemb) * sizeof(UINT64);		// 10 as a small safety factor

if(m_Sequences.pSeqStarts != NULL)
	{
	if(m_Sequences.AllocMemSeqStarts < (UINT64)ReqAllocMem || ((m_Sequences.AllocMemSeqStarts * 10)  > ((UINT64)ReqAllocMem * 12)))	// allow 20% float before reallocating
		{
#ifdef _WIN32
		free(m_Sequences.pSeqStarts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
		if(m_Sequences.pSeqStarts != MAP_FAILED)
			munmap(m_Sequences.pSeqStarts,m_Sequences.AllocMemSeqStarts);
#endif	
		m_Sequences.pSeqStarts = NULL;
		m_Sequences.AllocMemSeqStarts = 0;
		}
	}


if(m_Sequences.pSeqStarts == NULL || m_Sequences.AllocMemSeqStarts == 0)
	{
	m_Sequences.AllocMemSeqStarts = ReqAllocMem; 
#ifdef _WIN32
	m_Sequences.pSeqStarts = (UINT64 *) malloc((size_t)m_Sequences.AllocMemSeqStarts);	
	if(m_Sequences.pSeqStarts == NULL)
		{
		if(bSerialise)
			ReleaseSerialise();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenSeqStarts: Sequence start offsets array memory allocation of %llu bytes - %s",m_Sequences.AllocMemSeqStarts,strerror(errno));
		m_Sequences.AllocMemSeqStarts = 0;
		Reset(false);
		return(eBSFerrMem);
		}
#else
	if((m_Sequences.pSeqStarts = (UINT64 *)mmap(NULL,m_Sequences.AllocMemSeqStarts, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
		{
		if(bSerialise)
			ReleaseSerialise();
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenSeqStarts: Sequence start offsets array memory allocation of %llu bytes through mmap()  failed - %s",m_Sequences.AllocMemSeqStarts,strerror(errno));
		m_Sequences.pSeqStarts = NULL;
		m_Sequences.AllocMemSeqStarts = 0;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	memset(m_Sequences.pSeqStarts,0,(size_t)m_Sequences.AllocMemSeqStarts); // commits the memory!
	m_Sequences.NumSeqStarts = 0;
	UINT64 CurWorkSetSize = 0;
	CurWorkSetSize += m_Sequences.AllocMemSeqs2Assemb + m_Sequences.AllocMemSeqStarts + m_Sequences.AllocMemSfx + m_Sequences.AllocMemSeqFlags;
	if(CurWorkSetSize != m_CurMaxMemWorkSetBytes)
		{
		SetMaxMemWorkSetSize((size_t)CurWorkSetSize);
		}
	}


if(bGenFlags)
	{
	ReqAllocMem = (10 + m_Sequences.NumSeqs2Assemb) * sizeof(UINT16);		// 10 as a small safety factor

	if(m_Sequences.pSeqFlags != NULL)
		{
		if(m_Sequences.AllocMemSeqFlags < (UINT64)ReqAllocMem || ((m_Sequences.AllocMemSeqFlags * 10)  > ((UINT64)ReqAllocMem * 12)))	// allow 20% float before reallocating
			{
#ifdef _WIN32
			free(m_Sequences.pSeqFlags);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
			if(m_Sequences.pSeqFlags != MAP_FAILED)
				munmap(m_Sequences.pSeqFlags,m_Sequences.AllocMemSeqFlags);
#endif	
			m_Sequences.pSeqFlags = NULL;
			m_Sequences.AllocMemSeqFlags = 0;
			}
		}

	if(m_Sequences.pSeqFlags == NULL || m_Sequences.AllocMemSeqFlags == 0)
		{
		m_Sequences.AllocMemSeqFlags = ReqAllocMem; 
#ifdef _WIN32
		m_Sequences.pSeqFlags = (UINT16 *) malloc((size_t)m_Sequences.AllocMemSeqFlags);	
		if(m_Sequences.pSeqFlags == NULL)
			{
			if(bSerialise)
				ReleaseSerialise();
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenSeqStarts: Sequence flags array memory allocation of %llu bytes - %s",m_Sequences.AllocMemSeqFlags,strerror(errno));
			m_Sequences.AllocMemSeqFlags = 0;
			Reset(false);
			return(eBSFerrMem);
			}
#else
	if((m_Sequences.pSeqFlags = (UINT16 *)mmap(NULL,m_Sequences.AllocMemSeqFlags, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
			{
			if(bSerialise)
				ReleaseSerialise();
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenSeqStarts: Sequence flags memory allocation of %llu bytes through mmap()  failed - %s",m_Sequences.AllocMemSeqFlags,strerror(errno));
			m_Sequences.pSeqFlags = NULL;
			m_Sequences.AllocMemSeqFlags = 0;
			Reset(false);
			return(eBSFerrMem);
			}
#endif
		memset(m_Sequences.pSeqFlags,0,(size_t)m_Sequences.AllocMemSeqFlags); // commits the memory!
		m_Sequences.NumSeqFlags = 0;
		UINT64 CurWorkSetSize = 0;
		CurWorkSetSize = m_Sequences.AllocMemSeqs2Assemb + m_Sequences.AllocMemSeqStarts + m_Sequences.AllocMemSfx + m_Sequences.AllocMemSeqFlags;
		if(CurWorkSetSize != m_CurMaxMemWorkSetBytes)
			{
			SetMaxMemWorkSetSize((size_t)CurWorkSetSize);
			}
		}
	m_Sequences.NumSeqFlags = 0;
	}

pSeqStarts = m_Sequences.pSeqStarts;
pSeqFlags = m_Sequences.pSeqFlags;

UINT32 CurSeqLen;
UINT32 NumSeqWrds;
UINT32 SeqFlags;
tSeqWrd4 *pSeqWrd;
pSeqWrd = (tSeqWrd4 *)m_Sequences.pSeqs2Assemb;
pSeqWrd += 1;				
SeqWrdIdx = 1;
for(SeqID = 1; SeqID <= m_Sequences.NumSeqs2Assemb; SeqID++,pSeqStarts++)
	{
	GetSeqHeader(pSeqWrd,			// pts to a SeqWrd within the sequence
			NULL,				// returned 32 bit sequence identifier
			NULL,				// returned source file identifier
			&SeqFlags,				// returned 16 bit sequence flags
			&CurSeqLen,			// returned 30 bit sequence length
			bSerialise);				// need to serialise ?
	*pSeqStarts = SeqWrdIdx;
	if(bGenFlags)
		*pSeqFlags++ = (UINT16)SeqFlags;
	NumSeqWrds = 3 + (((CurSeqLen) + 14) / 15);
	SeqWrdIdx += NumSeqWrds;
	pSeqWrd += NumSeqWrds;
	}

m_Sequences.NumSeqStarts = m_Sequences.NumSeqs2Assemb;
if(bGenFlags)
	m_Sequences.NumSeqFlags = m_Sequences.NumSeqs2Assemb;
if(bSerialise)
	ReleaseSerialise();
gDiagnostics.DiagOut(eDLDiag,gszProcName,"GenSeqStarts: Completed sequence start offsets");
#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
return(eBSFSuccess);
}



// generate sfx over current concatenated packed sequences
// sparse suffix so generated will be at complete (any partial < 15bp tSeqWrds will not be indexed) tSeqWrd4 boundaries
// NOTE: sparse indexing (FirstNSeqWrds and ExcludeLastNSeqWrds only implemented for tSeqWrd4)
// NOTE: user can exclude any PE sequences from being indexed, useful when scaffolding
teBSFrsltCodes
CKangadna::GenRdsSfx(int FirstNSeqWrds,			// max number of SeqWrds (0 to index all), starting from 1st, to index in each read sequence
		             int ExcludeLastNSeqWrds,	// exclude last N SeqWrds in each read sequence from indexing, 0 to index FirstNSeqWrds
					 bool bExclPE)				// true to exclude sequences marked as being PE from being indexed
{
int ElSize;
UINT64 MaxSuffixEls;
UINT32 NumSfxEls;
bool bNoIndex;

gDiagnostics.DiagOut(eDLDiag,gszProcName,"GenRdsSfx: Initialising suffix array");

if(FirstNSeqWrds < 0)			// better safe than sorry
	FirstNSeqWrds = 0;
if(ExcludeLastNSeqWrds < 0)		
	ExcludeLastNSeqWrds = 0;
if(FirstNSeqWrds == 1)			// if just 1st SeqWrd to be indexed then excluding last ExcludeLastNSeqWrds becomes irrelevent
	ExcludeLastNSeqWrds = 0;

// firstly need to determine number of suffix elements required
MaxSuffixEls = 0;

tSeqWrd4 SeqWord;
tSeqWrd4 *pSeqWord;
UINT64 SeqWrdIdx;
UINT32 CurSeqLen;
pSeqWord = (tSeqWrd4 *)m_Sequences.pSeqs2Assemb;
pSeqWord += 1;   // skip initial cSeqWrd4BOS
NumSfxEls = 0;
CurSeqLen = 0;
bNoIndex = false;
UINT64 LastSeqs2AssembOfs;

LastSeqs2AssembOfs = m_Sequences.Seqs2AssembOfs;

for(SeqWrdIdx = 1; SeqWrdIdx < LastSeqs2AssembOfs; SeqWrdIdx++)
	{
	SeqWord = *pSeqWord++;
	if(SeqWord == cSeqWrd4EOS)
		break;
	if(SeqWord & cSeqWrd4LSWHdr)	// treating last partial SeqWrds as not to be indexed
		{
		if(bExclPE && ((SeqWord & cSeqWrd4MSWMsk) == cSeqWrd4MSWHdr))
			{
			if(SeqWord & 0x03)		// bit 1 set if PE sequence
				{
				CurSeqLen = 0;
				bNoIndex = true;
				continue;
				}
			bNoIndex = false;
			}

		if(!bNoIndex && ExcludeLastNSeqWrds && ((SeqWord & cSeqWrd4MSWMsk) == cSeqWrd4MSWHdr))
			{
			CurSeqLen = pSeqWord[1] & 0x3fffffff;
			CurSeqLen = (CurSeqLen + 14) / 15;					// number of sequence words including any partial final SeqWrd4
			if(CurSeqLen > (UINT32)ExcludeLastNSeqWrds)
				CurSeqLen -= ExcludeLastNSeqWrds;
			else
				CurSeqLen = 1;
			}
		else
			if(!bNoIndex && !ExcludeLastNSeqWrds)
				CurSeqLen = 0xffffffff;			// no limits on excluding last words from indexing					
		NumSfxEls = 0;
		continue;
		}
	if(CurSeqLen == 0 || (FirstNSeqWrds && NumSfxEls >= (UINT32)FirstNSeqWrds))
		continue;
	CurSeqLen -= 1;
	NumSfxEls += 1;
	MaxSuffixEls += 1;
	}

// number of suffix elements required is now known, next will be allocating memory
MaxSuffixEls += 16;							// allow for a small safety factor
ElSize = sizeof(UINT32);					// assume can use 4byte suffix elements
if(MaxSuffixEls > cMaxSfxBlkEls ||			// but if number of sfx elements > ~4G 
	LastSeqs2AssembOfs > cMaxSfxBlkEls)		// or tBaseWrds range is ~ 32bits 
	ElSize += 1;							// then need to use 5 byte elements

UINT64 ReqAllocMem;
ReqAllocMem = (MaxSuffixEls * (UINT64)ElSize);

if(m_Sequences.pSuffixArray != NULL && (m_Sequences.AllocMemSfx < ReqAllocMem || ((m_Sequences.AllocMemSfx * 10 ) > (ReqAllocMem * 12))))
	{
#ifdef _WIN32
	free(m_Sequences.pSuffixArray);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_Sequences.pSuffixArray != MAP_FAILED)
		munmap(m_Sequences.pSuffixArray,m_Sequences.AllocMemSfx);
#endif	
	m_Sequences.pSuffixArray = NULL;
	m_Sequences.AllocMemSfx = 0;
	}

if(m_Sequences.pSuffixArray == NULL || m_Sequences.AllocMemSfx == 0)
	{
	m_Sequences.AllocMemSfx = ReqAllocMem; 
#ifdef _WIN32
	m_Sequences.pSuffixArray = (void *) malloc((size_t)m_Sequences.AllocMemSfx);	
	if(m_Sequences.pSuffixArray == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenRdsSfx: Suffix array memory allocation of %llu bytes - %s",m_Sequences.AllocMemSfx,strerror(errno));
		m_Sequences.AllocMemSfx = 0;
		Reset(false);
		return(eBSFerrMem);
		}
#else
	if((m_Sequences.pSuffixArray = (void *)mmap(NULL,m_Sequences.AllocMemSfx, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenRdsSfx: Suffix array memory allocation of %llu bytes through mmap()  failed - %s",m_Sequences.AllocMemSfx,strerror(errno));
		m_Sequences.pSuffixArray = NULL;
		m_Sequences.AllocMemSfx = 0;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	memset(m_Sequences.pSuffixArray,0,(size_t)m_Sequences.AllocMemSfx); // commits the memory!

	UINT64 CurWorkSetSize = 0;
	CurWorkSetSize = m_Sequences.AllocMemSeqs2Assemb + m_Sequences.AllocMemSeqStarts + m_Sequences.AllocMemSfx + m_Sequences.AllocMemSeqFlags;
	if(CurWorkSetSize != m_CurMaxMemWorkSetBytes)
		{
		SetMaxMemWorkSetSize((size_t)CurWorkSetSize);
		}
	}

// allocation completed now into the real business of sorting
m_Sequences.SfxElSize = ElSize;

size_t SfxIdx;
pSeqWord = (tSeqWrd4 *)m_Sequences.pSeqs2Assemb;
pSeqWord += 1;
m_Sequences.NumSuffixEls = 0;
NumSfxEls = 0;
CurSeqLen = 0;	
bNoIndex = false;				
if(ElSize == 5)
	{
	UINT8 *pArr5 = (UINT8 *)m_Sequences.pSuffixArray;
	for(SeqWrdIdx = 1; SeqWrdIdx < LastSeqs2AssembOfs; SeqWrdIdx++)
		{
		SeqWord = *pSeqWord++;
		if(SeqWord == cSeqWrd4EOS)
			break;
		if(SeqWord & cSeqWrd4LSWHdr)		// could be partial SeqWrd or new header, these are not indexed
			{
			if(bExclPE && ((SeqWord & cSeqWrd4MSWMsk) == cSeqWrd4MSWHdr))
				{
				if(SeqWord & 0x03)		// bit 1 set if PE sequence
					{
					CurSeqLen = 0;
					bNoIndex = true;
					continue;
					}
				bNoIndex = false;
				}

			if(!bNoIndex && ExcludeLastNSeqWrds && ((SeqWord & cSeqWrd4MSWMsk) == cSeqWrd4MSWHdr))
				{
				CurSeqLen = pSeqWord[1] & 0x3fffffff;
				CurSeqLen = (CurSeqLen + 14) / 15;					// number of sequence words including any partial final SeqWrd4
				if(CurSeqLen > (UINT32)ExcludeLastNSeqWrds)
					CurSeqLen -= ExcludeLastNSeqWrds;
				else
					CurSeqLen = 1;
				}
			else
				if(!bNoIndex && !ExcludeLastNSeqWrds)
					CurSeqLen = 0xffffffff;			// no limits on excluding last words from indexing

			NumSfxEls = 0;
			continue;
			}
		if(CurSeqLen == 0 || (FirstNSeqWrds && NumSfxEls >= (UINT32)FirstNSeqWrds))
			continue;
		CurSeqLen -= 1;
		NumSfxEls += 1;
		pArr5 = Pack5(SeqWrdIdx,pArr5);
		m_Sequences.NumSuffixEls += 1;
		}
	Pack5(0xffffffffff,pArr5);
	m_xpConcatSeqs = (UINT8 *)m_Sequences.pSeqs2Assemb;
	gDiagnostics.DiagOut(eDLDiag,gszProcName,"GenRdsSfx: Sparse suffix array contains %lld index elements size %d bytes..",m_Sequences.NumSuffixEls,ElSize);
	gDiagnostics.DiagOut(eDLDiag,gszProcName,"GenRdsSfx: Now sorting...");
	m_MTqsort.qsort(m_Sequences.pSuffixArray,m_Sequences.NumSuffixEls,5,Sfx5SortSeqWrd4Func);
	gDiagnostics.DiagOut(eDLDiag,gszProcName,"GenRdsSfx: Sorting completed...");
	}
else
	{
	SfxIdx = 0;
	UINT32 *pArr4 = (UINT32 *)m_Sequences.pSuffixArray;
	UINT64 PrevSfxWrdIdx = 0;
	for(SeqWrdIdx = 1; SeqWrdIdx < LastSeqs2AssembOfs; SeqWrdIdx++)
		{
		SeqWord = *pSeqWord++;
		if(SeqWord  & cSeqWrd4LSWHdr)	// could be partial SeqWrd or new header, these are not indexed
			{
			if(bExclPE && ((SeqWord & cSeqWrd4MSWMsk) == cSeqWrd4MSWHdr))
				{
				if(SeqWord & 0x03)		// bit 1 set if PE sequence
					{
					CurSeqLen = 0;
					bNoIndex = true;
					continue;
					}
				bNoIndex = false;
				}
			if(!bNoIndex && ExcludeLastNSeqWrds && ((SeqWord & cSeqWrd4MSWMsk) == cSeqWrd4MSWHdr))
				{
				CurSeqLen = pSeqWord[1] & 0x3fffffff;
				CurSeqLen = (CurSeqLen + 14) / 15;					// number of sequence words including any partial final SeqWrd4
				if(CurSeqLen > (UINT32)ExcludeLastNSeqWrds)
					CurSeqLen -= ExcludeLastNSeqWrds;
				else
					CurSeqLen = 1;
				}
			else
				if(!bNoIndex && !ExcludeLastNSeqWrds)
					CurSeqLen = 0xffffffff;   // no limits on excluding last words from indexing
			NumSfxEls = 0;
			continue;
			}
		if(CurSeqLen == 0 || (FirstNSeqWrds && NumSfxEls >= (UINT32)FirstNSeqWrds))
			continue;
		CurSeqLen -= 1;
		NumSfxEls += 1;
		pArr4[SfxIdx++] = (UINT32)SeqWrdIdx;
		m_Sequences.NumSuffixEls += 1;
		}
	pArr4[SfxIdx] = 0xffffffff;
	m_xpConcatSeqs = (UINT8 *)m_Sequences.pSeqs2Assemb;
	gDiagnostics.DiagOut(eDLDiag,gszProcName,"GenRdsSfx: Sparse suffix array contains %lld index elements size %d bytes..",m_Sequences.NumSuffixEls,ElSize);
	gDiagnostics.DiagOut(eDLDiag,gszProcName,"GenRdsSfx: Now sorting...");
	m_MTqsort.qsort(m_Sequences.pSuffixArray,m_Sequences.NumSuffixEls,sizeof(UINT32),SfxSortSeqWrd4Func);
	gDiagnostics.DiagOut(eDLDiag,gszProcName,"GenRdsSfx: Sorting completed...");
	}


#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
gDiagnostics.DiagOut(eDLDiag,gszProcName,"GenRdsSfx: Suffix array generation completed");
return(eBSFSuccess);
}

// ChunkedWrite
// Seeks to specified 64bit file offset and writes to disk as chunks of no more than INT_MAX/16  
teBSFrsltCodes
CKangadna::ChunkedWrite(UINT64 WrtOfs,UINT8 *pData,UINT64 WrtLen)
{
return(ChunkedWrite(m_hOutFile,m_szOutFile,WrtOfs,pData,WrtLen));
}

teBSFrsltCodes 
CKangadna::ChunkedWrite(int hFile,char *pszFile,UINT64 WrtOfs,UINT8 *pData,UINT64 WrtLen)
{
UINT32 BlockLen;
if(_lseeki64(hFile,WrtOfs,SEEK_SET) != WrtOfs)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to seek to %ld on file %s - error %s",WrtOfs,pszFile,strerror(errno));
	return(eBSFerrFileAccess);
	}

while(WrtLen)
	{
	BlockLen = WrtLen > (UINT64)(INT_MAX/8) ? (INT_MAX/8) : (UINT32)WrtLen;
	WrtLen -= BlockLen;
	if(write(hFile,pData,(int)BlockLen)!=(int)BlockLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write to disk on file %s - error %s",pszFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	pData += BlockLen;
	}
return(eBSFSuccess);
}

// ChunkedRead
// Seeks to specified 64bit file offset and reads from disk as chunks of no more than INT_MAX/32
teBSFrsltCodes
CKangadna::ChunkedRead(UINT64 RdOfs,UINT8 *pData,UINT64 RdLen)
{
return(ChunkedRead(m_hInFile,m_szInFile,RdOfs,pData,RdLen));
}

// ChunkedRead
// Seeks to specified 64bit file offset and reads from disk as chunks of no more than INT_MAX/8
teBSFrsltCodes
CKangadna::ChunkedRead(int hFile,char *pszFile,UINT64 RdOfs,UINT8 *pData,UINT64 RdLen)
{
UINT32 BlockLen;
int BytesRead;
if(_lseeki64(hFile,(INT64)RdOfs,SEEK_SET) != RdOfs)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to seek to %llu on file - error %s",RdOfs,pszFile,strerror(errno));
	return(eBSFerrFileAccess);
	}

while(RdLen)
	{
	BlockLen = RdLen > (UINT64)(INT_MAX/4) ? (INT_MAX/4) : (UINT32)RdLen;
	RdLen -= (UINT64)BlockLen;
	if((BytesRead = read(hFile,pData,BlockLen))!=BlockLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to read from file %s - error %s",pszFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	pData += BlockLen;
	}
return(eBSFSuccess);
}

int								// returned sequence length, will be limited to MaxSeqLen
CKangadna::GetSeq(void *pStartSeqWrd,
		UINT8 *pRetSeq,			// where to copy unpacked sequence bases
		UINT32 MaxSeqLen,		// return at most MaxSeqLen bases (0 if no limit)
		bool bNoEOSTerm)		// option to not terminate unpacked string with EOS
{
UINT32 Msk;
int ShfFct;
int RetLen;
UINT32 Idx;

if(MaxSeqLen == 0)
	MaxSeqLen = 0xffffffff;
RetLen = 0;

tSeqWrd4 *pSeqWrd = (tSeqWrd4 *)pStartSeqWrd;
tSeqWrd4 SeqWrd;
while(MaxSeqLen && (((SeqWrd = *pSeqWrd++) & cSeqWrd4MSWHdr) != cSeqWrd4MSWHdr))
	{
	if(SeqWrd & cSeqWrd4PartSeq)
		{
		Msk = 0x00003;
		for(Idx = 0; Idx < 14; Idx++)
			{
			if(!(SeqWrd & Msk))
				{
				MaxSeqLen = min(MaxSeqLen,14 - Idx);
				break;
				}
			Msk <<= 2;
			}
		}
	ShfFct = 30;
	do {
		ShfFct -= 2;
		*pRetSeq++ = (SeqWrd >> ShfFct) & 0x03;
		MaxSeqLen -= 1;
		RetLen += 1;
		}
	while(MaxSeqLen && ShfFct);
	}

if(!bNoEOSTerm)
	*pRetSeq = eBaseEOS;
return(RetLen);
}


int							// returned sequence length, will be limited to MaxSeqLen
CKangadna::GetSeq(tSeqID SeqID,				// sequence identifier
			UINT8 *pRetSeq,				// where to copy unpacked sequence bases
			UINT32 MaxSeqLen,		// return at most MaxSeqLen bases (0 if no limit)
			bool bNoEOSTerm)	// option to not terminate unpacked string with EOS
{
void *pStartSeqWrd;
UINT32 SeqLen;
int ShfFct;
int RetLen;
int Idx;

if(SeqID < 1 || SeqID > m_Sequences.NumSeqs2Assemb)
	return(-1);

if((pStartSeqWrd = GetSeqHeader(SeqID,NULL,NULL,&SeqLen))==NULL)
	return(-1);

if(MaxSeqLen == 0)
	MaxSeqLen = SeqLen;

tSeqWrd4 *pSeqWrd = (tSeqWrd4 *)pStartSeqWrd;
tSeqWrd4 SeqWrd = *pSeqWrd++;
ShfFct = 28;
RetLen = min(MaxSeqLen,SeqLen);
for(Idx = 0; Idx < RetLen; Idx++)
	{
	if(Idx && !(Idx % 15))
		{
		ShfFct = 28;
		SeqWrd = *pSeqWrd++;
		}
	*pRetSeq++ = (SeqWrd >> ShfFct) & 0x03;
	ShfFct -= 2;
	}

if(!bNoEOSTerm)
	*pRetSeq = eBaseEOS;
return(RetLen);
}


// IterSeqHeaders
// Iterates to next sequence immediately following the current sequence
// If current sequence (pCurSeq) is NULL then returns the first sequence
tSeqWrd4 *								// returned ptr to 1st word of actual packed sequence (use this as pCurSeq on next invocation)
CKangadna::IterSeqHeaders(tSeqWrd4 *pCurSeq, // iterate to next sequence following this
			tSeqID *pSeqID,			// returned sequence identifier
			UINT32 *pSrcFileID,		// returned 8 bit source file identifier
			UINT32 *pFlgs,			// returned 16 bit sequence flags
			UINT32 *pSeqLen,		// returned 30 bit sequence length
    		bool bSerialise)		// true if access to headers are to be serialised
{
bool bNoHdr;
if(m_Sequences.NumSeqs2Assemb == 0)
	return(NULL);

if(pSeqID == NULL && pSrcFileID == NULL && pFlgs == NULL && pSeqLen == NULL)
	bNoHdr = true;
else
	bNoHdr = false;

tSeqWrd4 *pSeqHdrWrd;
tSeqWrd4 SeqWrd;
if(pCurSeq == NULL)
	pSeqHdrWrd = (tSeqWrd4 *)m_Sequences.pSeqs2Assemb;
else
	pSeqHdrWrd = (tSeqWrd4 *)pCurSeq;
if(*pSeqHdrWrd == cSeqWrd4EOS || pSeqHdrWrd[1] == cSeqWrd4EOS)		// ensure not trying to get header for EOS
	return(NULL);
pSeqHdrWrd += 1;
while(((SeqWrd = *pSeqHdrWrd) & cSeqWrd4LSWHdr) != cSeqWrd4MSWHdr)	// advance until the most significant header word or EOS
	{
	pSeqHdrWrd += 1;
	if(*pSeqHdrWrd == cSeqWrd4EOS)			// ensure not trying to get header for EOS
		return(NULL);
	}
if(bNoHdr)
	{
	while(*pSeqHdrWrd & cSeqWrd4MSWHdr)	// advance until sequence
		pSeqHdrWrd += 1;
	return(pSeqHdrWrd);
	}
pCurSeq = pSeqHdrWrd;

return(GetSeqHeader(pCurSeq,pSeqID,pSrcFileID,pFlgs,pSeqLen,bSerialise));
}


tSeqWrd4 *								// returned ptr to 1st word of actual packed sequence which immediately follows the header words
CKangadna::GetSeqHeader(tSeqID SeqID,	// 32 bit sequence identifier
			UINT32 *pSrcFileID,			// returned 8 bit source file identifier
			UINT32 *pFlgs,				// returned 16 bit sequence flags
			UINT32 *pSeqLen,			// returned 30 bit sequence length
			bool bSerialise)			// set true if access to headers are required to be serialise
{
tSeqWrd4 *pSeqStart;
UINT64 *pSeqStarts;

if(SeqID < 1 || SeqID > (tSeqID)m_Sequences.NumSeqs2Assemb || m_Sequences.pSeqStarts == NULL || m_Sequences.NumSeqStarts < SeqID)
	return(NULL);
pSeqStarts = &m_Sequences.pSeqStarts[SeqID-1];

tSeqWrd4 *pSeqWrd;
pSeqWrd = (tSeqWrd4 *)m_Sequences.pSeqs2Assemb;
pSeqStart = &pSeqWrd[*pSeqStarts];

pSeqStart =	GetSeqHeader(pSeqStart,		// pts to a SeqWrd within the sequence
				NULL,					
				pSrcFileID,				// returned 8bit source file identifier
				pFlgs,					// returned 16bit sequence flags 
				pSeqLen,				// returned 30bit sequence length
				bSerialise);
return(pSeqStart);
}


int // returns sequence length for sequence with specified SeqID
CKangadna::GetSeqLen(tSeqID SeqID) // sequence identifier
{
UINT32 SeqLen;
if(GetSeqHeader(SeqID,			        // 32 bit sequence identifier
			NULL,					// not interested in source file identifier
			NULL,					// not interested in sequence flags
			&SeqLen)==NULL)			// returned 30 bit sequence length
	return(-1);
return(SeqLen);
}

UINT64			// index+1 in pSfxArray of first exactly matching probe or 0 if no match				
CKangadna::LocateFirstExact(int ElSize,		// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
				  void *pProbe,				// pts to probe sequence
				  int ProbeLen,				// probe length (in bases, not tSeqWrd's) to exactly match over
				  void *pTarg,				// target sequence
				  UINT8 *pSfxArray,			// target sequence suffix array
				  UINT64 SfxLo,				// low index in pSfxArray
				  UINT64 SfxHi)				// high index in pSfxArray
{
void *pEl1;
void *pEl2;
tSeqWrd4 El1;
tSeqWrd4 El2;
int CmpRslt;
UINT64 Mark;
UINT64 TargPsn;
UINT64 TargEl;
UINT64 TargIdx;

do {
	pEl1 = pProbe;
	TargPsn = ((UINT64)SfxLo + SfxHi) / 2L;
	TargEl = TargPsn * ElSize;

	if(ElSize == 4)
		TargIdx = *(UINT32 *)&pSfxArray[TargEl];
	else
		TargIdx = Unpack5(&pSfxArray[TargEl]);

	pEl2 =  &((tSeqWrd4 *)pTarg)[TargIdx];

	if((El1 = *(tSeqWrd4 *)pEl1) < (El2 = *(tSeqWrd4 *)pEl2))  // assuming that the probe and target are both at least 16bp
		CmpRslt = -1;
	else
		if(El1 > El2)
			CmpRslt = 1;
		else
			CmpRslt = CmpPackedSeqs((tSeqWrd4 *)pEl1,(tSeqWrd4 *)pEl2,ProbeLen);

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

			TargEl = TargPsn * ElSize;

			if(ElSize == 4)
				TargIdx = *(UINT32 *)&pSfxArray[TargEl];
			else
				TargIdx = Unpack5(&pSfxArray[TargEl]);

			pEl2 =  &((tSeqWrd4 *)pTarg)[TargIdx];

			pEl1 = pProbe;
			if((El1 = *(tSeqWrd4 *)pEl1) < (El2 = *(tSeqWrd4 *)pEl2))  // assuming that the probe and target are both at least 16bp
				CmpRslt = -1;
			else
				if(El1 > El2)
					CmpRslt = 1;
				else
					CmpRslt = CmpPackedSeqs((tSeqWrd4 *)pEl1,(tSeqWrd4 *)pEl2,ProbeLen);
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
CKangadna::LocateLastExact(int ElSize,		// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
				  void *pProbe,				// pts to probe sequence
				  int ProbeLen,					// probe length (bases, not tSeqWrd4's) to exactly match over
				  void *pTarg,				// target sequence
				  UINT8 *pSfxArray,				// target sequence suffix array
				  UINT64 SfxLo,					// low index in pSfxArray
				  UINT64 SfxHi,					// high index in pSfxArray
				  UINT32 Limit)					// if non-zero then need only iterate towards last exactly matching this for this Limit iterations
{
void *pEl1;
void *pEl2;
int CmpRslt;
UINT64 TargEl;
UINT64 Mark;
UINT64 TargPsn;
UINT64 LoTargPsn;
UINT64 SfxHiMax = SfxHi;
UINT64 SfxLoMax = SfxLo;
UINT64 TargIdx;

if(Limit > 0)					// Limit is a soft limit so add a few more on 
	Limit += 10;

do {
	pEl1 = pProbe;
	TargPsn = (SfxLo + SfxHi) / 2L;
	TargEl = TargPsn * ElSize;

	if(ElSize == 4)
		TargIdx = *(UINT32 *)&pSfxArray[TargEl];
	else
		TargIdx = Unpack5(&pSfxArray[TargEl]);

	pEl2 =  &((tSeqWrd4 *)pTarg)[TargIdx];

	CmpRslt = CmpPackedSeqs((tSeqWrd4 *)pEl1,(tSeqWrd4 *)pEl2,ProbeLen);
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

			TargEl = TargPsn * ElSize;

			if(ElSize == 4)
				TargIdx = *(UINT32 *)&pSfxArray[TargEl];
			else
				TargIdx = Unpack5(&pSfxArray[TargEl]);

			pEl2 =  &((tSeqWrd4 *)pTarg)[TargIdx];

			pEl1 = pProbe;
			CmpRslt = CmpPackedSeqs((tSeqWrd4 *)pEl1,(tSeqWrd4 *)pEl2,ProbeLen);
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



// GetPackedSeqID
// Get sequence identifier for packed sequence ptd to by pSeqWrd
// Returns sequence identifier or 0 if errors
tSeqID
CKangadna::GetPackedSeqID(tSeqWrd4 *pSeqWrd)			// ptr into sequence for which sequence identifer is to be returned
{
tSeqID SeqID;

if(GetSeqHeader(pSeqWrd,&SeqID)==NULL)
	SeqID = 0;
return(SeqID);
}

// GetPackedSeq3BaseLen
// NOTE: assumes that the sequence has been terminated by EOS or a header word
int											// returns number of bases 3' right, includes initial tSeqWrd ptd at by pSeqWrd
CKangadna::GetPackedSeq3BaseLen(void *pSeqWrd)
{
return(GetPackedSeq3BaseLen((tSeqWrd4 *)pSeqWrd));
}



// GetPackedSeq3BaseLen
// Returns -1 if errors
int    
CKangadna::GetPackedSeq3BaseLen(tSeqWrd4 *pSeqWrd) // returns number of bases 3' right including tSeqWrd4 ptd at by pSeqWrd
{
int Len;
tSeqWrd4 SeqWrd;
if(pSeqWrd == NULL)
	return(-1);
if(*pSeqWrd & cSeqWrd4MSWHdr)
	return(0);
Len = 0;
while(!((SeqWrd = *pSeqWrd++) & cSeqWrd4MSWHdr))	
	{
	if(!(SeqWrd & cSeqWrd4PartSeq))
		Len += 15;
	else
		{
		Len += 14;
		while(SeqWrd & 0x03)
			{
			SeqWrd >>=2;
			Len -= 1;
			}
		break;
		}
	}
return(Len);
}


// GetPackedSeq3BaseLen
// Returns -1 if errors

int  // SfxSortFunc for UINT32 suffix elements
CKangadna::SfxSortSeqWrd4Func(const void *arg1, const void *arg2)
{
UINT64 SeqIdx1 = *(UINT32 *)arg1;
UINT64 SeqIdx2 = *(UINT32 *)arg2;
tSeqWrd4 *pSeq1 =  (tSeqWrd4 *)&m_xpConcatSeqs[SeqIdx1 * 4];
tSeqWrd4 *pSeq2 = (tSeqWrd4 *)&m_xpConcatSeqs[SeqIdx2 * 4];
return(CmpPackedSeqs(pSeq1,pSeq2,cMaxSortSfxLen));
}

int  // SfxSortFunc for 5byte suffix elements
CKangadna::Sfx5SortSeqWrd4Func(const void *arg1, const void *arg2)
{
UINT64 SeqIdx1 = Unpack5((UINT8 *)arg1);
UINT64 SeqIdx2 = Unpack5((UINT8 *)arg2);
tSeqWrd4 *pSeq1 = (tSeqWrd4 *)&m_xpConcatSeqs[SeqIdx1 * 4];
tSeqWrd4 *pSeq2 = (tSeqWrd4 *)&m_xpConcatSeqs[SeqIdx2 * 4];
return(CmpPackedSeqs(pSeq1,pSeq2,cMaxSortSfxLen));
}

int									// number of tSeqWrds returned in pDstPackedSeq (excludes optional EOS)
CKangadna::GetPackedSeq(int MaxSeqWrds,	// limit to this many tSeqWrds (0 if no limit)
			tSeqWrd4 *pSrcPackedSeq,	// get from this packed sequence
			tSeqWrd4 *pDstPackedSeq,   	// copy into this packed sequence
			bool bNoEOSTerm)			// option to not terminate packed substring with EOS
{
tSeqWrd4 SeqWrd;
UINT32 NumSeqWrds;
if(MaxSeqWrds == 0)
	MaxSeqWrds = 0xffffffff;

NumSeqWrds = 0;
while(MaxSeqWrds-- && ((SeqWrd = *pSrcPackedSeq++) & cSeqWrd4MSWHdr) != cSeqWrd4MSWHdr)
	{
	*pDstPackedSeq++ = SeqWrd;
	NumSeqWrds += 1;
	}
if(!bNoEOSTerm)
	*pDstPackedSeq = cSeqWrd4EOS;
return(NumSeqWrds);
}


// GetSeqWrdSubSeq
// NOTE: GetSeqWrdSubSeq assumes that the source sequence contains a preceding header containing the actual source length 
int												// number of packed bases returned in pDstPackedSeq (excludes optional EOS) 
CKangadna::GetSeqWrdSubSeq(int SubOfs,			// substring starting at this base offset within pSrcPackedSeq
				 int SubLen,					// substring required which contains at most this many bases (0 for maximal length)
				 tSeqWrd4 *pSrcPackedSeq,		// extract substring from this packed string
				 tSeqWrd4 *pDstPackedSeq,		// extract into this packed substring
				bool bNoEOSTerm)				// option to not terminate packed substring with EOS
{
int Rslt;
int SrcLen;
int ShfLft;
int ShfRgt;
tSeqWrd4 DstSeqWrd;
tSeqWrd4 SrcSeqWrd;
tSeqWrd4 Msk;
tSeqWrd4 Flg;

if(!bNoEOSTerm && pDstPackedSeq != NULL)		
	*pDstPackedSeq = cSeqWrd4EOS; 
if(SubOfs < 0 || SubLen < 0 || pSrcPackedSeq == NULL || pDstPackedSeq == NULL)
	return(-1);

// determine source packed sequence length
SrcLen = GetPackedSeq3BaseLen(pSrcPackedSeq);
if(SubOfs >= SrcLen)
	return(0);
if(!SubLen)
	SubLen = SrcLen;
SubLen = Rslt = min(SrcLen - SubOfs,SubLen);
pSrcPackedSeq += SubOfs / 15;
ShfLft = 2 * (SubOfs % 15);
ShfRgt = 30 - ShfLft;

while(SubLen)
	{
	DstSeqWrd = *pSrcPackedSeq++;
	if(ShfLft)
		{
		DstSeqWrd <<= ShfLft;
		if(!((SrcSeqWrd = *pSrcPackedSeq) & cSeqWrd4MSWHdr))
			DstSeqWrd |= (SrcSeqWrd & cSeqWrd4Msk) >> ShfRgt;
		}
	DstSeqWrd &= cSeqWrd4Msk;
	if(SubLen < 15)
		{
		Msk = 0x03;			
		Flg = 0x03;
		while(SubLen++ < 15)
			{
			DstSeqWrd &= ~Msk;
			if(SubLen != 15)
				DstSeqWrd |= Flg;
			Msk <<= 2;
			Flg <<= 2;
			}
		SubLen = 15;
		DstSeqWrd |= cSeqWrd4PartSeq;
		}
	SubLen -= 15;
	*pDstPackedSeq++ = DstSeqWrd;
	}
if(!bNoEOSTerm)
	*pDstPackedSeq = cSeqWrd4EOS;
return(Rslt);
}


// NOTE: GetNHSeqWrdSubSeq assumes that the source sequence does not contain a preceding header with actual source length and thus expects SubOfs + SubLen to be <= actual source length !!! 
int											// number of packed bases returned in pDstPackedSeq (excludes optional EOS)
CKangadna::GetNHSeqWrdSubSeq(int SubOfs,	// substring starting at this base offset within pSrcPackedSeq (SubOfs + SubLen expected to be <= pSrcPackedSeq length)
			int SubLen,						// substring required which contains at most this many bases
			tSeqWrd4 *pSrcPackedSeq,		// extract substring from this packed string, , may be same as pDstPackedSeq if inplace substring required
			tSeqWrd4 *pDstPackedSeq,		// extract into this packed substring,, may be same as pSrcPackedSeq if inplace substring required
			bool bNoEOSTerm)				// option to not terminate packed substring with EOS
{
int Rslt;
int ShfLft;
int ShfRgt;
tSeqWrd4 DstSeqWrd;
tSeqWrd4 SrcSeqWrd;
tSeqWrd4 Msk;
tSeqWrd4 Flg;


if(SubOfs < 0 || SubLen < 1 || pSrcPackedSeq == NULL || pDstPackedSeq == NULL)
	return(-1);
if(!bNoEOSTerm && pDstPackedSeq != pSrcPackedSeq)		
	*pDstPackedSeq = cSeqWrd4EOS; 

Rslt = SubLen;
pSrcPackedSeq += SubOfs / 15;
ShfLft = 2 * (SubOfs % 15);
ShfRgt = 30 - ShfLft;

while(SubLen)
	{
	DstSeqWrd = *pSrcPackedSeq++;
	if(ShfLft)
		{
		DstSeqWrd <<= ShfLft;
		if(!((SrcSeqWrd = *pSrcPackedSeq) & cSeqWrd4MSWHdr))
			DstSeqWrd |= (SrcSeqWrd & cSeqWrd4Msk) >> ShfRgt;
		}
	DstSeqWrd &= cSeqWrd4Msk;
	if(SubLen < 15)
		{
		Msk = 0x03;			
		Flg = 0x03;
		while(SubLen++ < 15)
			{
			DstSeqWrd &= ~Msk;
			if(SubLen != 15)
				DstSeqWrd |= Flg;
			Msk <<= 2;
			Flg <<= 2;
			}
		SubLen = 15;
		DstSeqWrd |= cSeqWrd4PartSeq;
		}
	SubLen -= 15;
	*pDstPackedSeq++ = DstSeqWrd;
	}
if(!bNoEOSTerm)
	*pDstPackedSeq = cSeqWrd4EOS;
return(Rslt);
}

// NOTE: NHSeqTrim assumes that the source sequence does not contain a preceding header with actual source length and thus expects Trim5 + Trim3 to be <= actual source length !!! 
int						// returned trimmed sequence length
CKangadna::NHSeqTrim( int Trim5,		// trim 5' by this many bases
		 int Trim3,		// trim 3' by this many bases
		 int SeqLen,	// number of bases in sequence before sequence trimmed
		 tSeqWrd4 *pSeq, // inplace trimming of this sequence
		 bool bNoEOSTerm)				// option to not terminate packed substring with EOS
{
int Rslt;
int SubLen;
int ShfLft;
int ShfRgt;
tSeqWrd4 DstSeqWrd;
tSeqWrd4 SrcSeqWrd;
tSeqWrd4 Msk;
tSeqWrd4 Flg;

tSeqWrd4 *pSrcPackedSeq;

if(Trim5 + Trim3 >= SeqLen || pSeq == NULL)
	return(-1);
SubLen = SeqLen - (Trim5 + Trim3);
Rslt = SubLen;
pSrcPackedSeq = pSeq + (Trim5 / 15);
ShfLft = 2 * (Trim5 % 15);
ShfRgt = 30 - ShfLft;

while(SubLen)
	{
	DstSeqWrd = *pSrcPackedSeq++;
	if(ShfLft)
		{
		DstSeqWrd <<= ShfLft;
		if(!((SrcSeqWrd = *pSrcPackedSeq) & cSeqWrd4MSWHdr))
			DstSeqWrd |= (SrcSeqWrd & cSeqWrd4Msk) >> ShfRgt;
		}
	DstSeqWrd &= cSeqWrd4Msk;
	if(SubLen < 15)
		{
		Msk = 0x03;			
		Flg = 0x03;
		while(SubLen++ < 15)
			{
			DstSeqWrd &= ~Msk;
			if(SubLen != 15)
				DstSeqWrd |= Flg;
			Msk <<= 2;
			Flg <<= 2;
			}
		SubLen = 15;
		DstSeqWrd |= cSeqWrd4PartSeq;
		}
	SubLen -= 15;
	*pSeq++ = DstSeqWrd;
	}
if(!bNoEOSTerm)
	*pSeq = cSeqWrd4EOS;
return(Rslt);
}


int											// 0 Seq1 == Seq2, -1 if Seq1 < Seq2, +1 if Seq1 > Seq2
CKangadna::CmpPackedSeqs(tSeqWrd4 *pProbeSeq,	// Seq1 (probe) packed sequence
			  tSeqWrd4 *pTargSeq,				// Seq2 (target) packed sequence
			  int MaxCmpLen)				// compare over at most this many bases
{
int NumSeq1Bases;
int NumSeq2Bases;
if(!MaxCmpLen || pProbeSeq == NULL || pTargSeq == NULL)
	return(0);

tSeqWrd4 BaseMsk;
tSeqWrd4 *pSeq1 = pProbeSeq;
tSeqWrd4 *pSeq2 =  pTargSeq;
tSeqWrd4 SeqWrd1;
tSeqWrd4 SeqWrd2;

SeqWrd1 = *pSeq1++;
SeqWrd2 = *pSeq2++;

while(MaxCmpLen >= 15)
	{
	if(SeqWrd1 & cSeqWrd4MSWHdr || SeqWrd2 & cSeqWrd4MSWHdr) // check if finished either sequence and starting to process next sequence
		{
		if(!(SeqWrd2 & cSeqWrd4MSWHdr))						  
			return(-1);
		if(!(SeqWrd1 & cSeqWrd4MSWHdr))						  
			return(1);
		return(0);
		}
				
	// still processing sequences
	// both words contain the max of 15 packed bases?
	if(!(SeqWrd1 & cSeqWrd4PartSeq) && !(SeqWrd2 & cSeqWrd4PartSeq)) 
		{
		if(SeqWrd1 < SeqWrd2)
			return(-1);
		if(SeqWrd1 > SeqWrd2)
			return(1);

		SeqWrd1 = *pSeq1++;
		SeqWrd2 = *pSeq2++;
		MaxCmpLen -= 15;
		}
	else	// else either or both sequence words have less than 15 bases
		break;
	}
		
// one or both words contains less than max of 15 packed bases
if(!MaxCmpLen || SeqWrd1 == SeqWrd2)								// if words still identical then sequences must have been identical
	return(0);

if(SeqWrd1 & cSeqWrd4MSWHdr || SeqWrd2 & cSeqWrd4MSWHdr) // check if finished either sequence and starting to process next sequence
	{
	if(!(SeqWrd2 & cSeqWrd4MSWHdr))						  
		return(-1);
	if(!(SeqWrd1 & cSeqWrd4MSWHdr))						  
		return(1);
	return(0);
	}

// determine number of bases in SeqWrd1
if(SeqWrd1 > cSeqWrd4Msk)
	{
	BaseMsk = 0x03;
	for(NumSeq1Bases = 14; NumSeq1Bases >= 0; NumSeq1Bases--,BaseMsk <<= 2)
		if(!(BaseMsk & SeqWrd1))
			break;
	}
else
	NumSeq1Bases = 15;

// determine number of bases in SeqWrd2
if(SeqWrd2 > cSeqWrd4Msk)
	{
	BaseMsk = 0x03;
	for(NumSeq2Bases = 14; NumSeq2Bases >= 0; NumSeq2Bases--,BaseMsk <<= 2)
		if(!(BaseMsk & SeqWrd2))
			break;
	}
else
	NumSeq2Bases = 15;

BaseMsk = 0x30000000;
do
	{
	if((SeqWrd1 & BaseMsk) < (SeqWrd2 & BaseMsk))
		return(-1);
	if((SeqWrd1 & BaseMsk) > (SeqWrd2 & BaseMsk))
		return(1);
	MaxCmpLen -= 1;
	if(!MaxCmpLen)
		return(0);
	BaseMsk >>= 2;
	NumSeq1Bases -= 1;
	NumSeq2Bases -= 1;
	}
while(NumSeq1Bases > 0 && NumSeq2Bases > 0);

if(NumSeq1Bases < NumSeq2Bases)
	return(-1);
if(NumSeq1Bases > NumSeq2Bases)
	return(1);
return(0);
}

int			// inplace shift left by one base the sequence ptd at by pSeq - note that sequence header is not modified!
CKangadna::ShfLeftPackedSeq(int CurLen,		// Seq1 currently contains this many bases, will contain CurLen-1 after shift left
							tSeqWrd4 *pSeq)		// Seq1 (probe) packed sequence to shift left
{
int Shfs;
if(pSeq == NULL || CurLen <= 1)			// need at least 2 bases to shift left and still return 1 base!
	return(0);

Shfs = CurLen;

tSeqWrd4 *pSeq1 = pSeq;
tSeqWrd4 *pSeq2;
tSeqWrd4 SeqWrd;
tSeqWrd4 SeqWrd1;

// pSeq may pt to one of the sequence header words, skip over these
while(*pSeq1 & cSeqWrd4MSWHdr)
	{
	if((*pSeq1 & 0x0e0000001) == cSeqWrd4EOS)
		return(0);
	pSeq1 += 1;
	}

//pSeq now pts into sequence words containing bases
pSeq2 = pSeq1;
while(Shfs && !((SeqWrd = *pSeq1++) & cSeqWrd4MSWHdr)) // iterate over sequence words until another sequence header word or all bases have been left shifted
	{
	SeqWrd1 = SeqWrd;
	SeqWrd <<= 2;								// left shift
	SeqWrd &= 0x3ffffffc;						// strip out top 2 bits
	if(Shfs <= 15 || SeqWrd1 & cSeqWrd4PartSeq)	// just shifted last word?
		{
		SeqWrd |= cSeqWrd4PartSeq;
		SeqWrd |= 0x03;
		Shfs = 0;								// completed shift lefts...
		}
	else										// still staying as a full sequence word
		{
		SeqWrd |= ((*pSeq1 >> 28) & 0x03);
		SeqWrd &= cSeqWrd4Msk;
		Shfs -= 15;
		}
	*pSeq2++ = SeqWrd;							// it's inplace shift left so write back...
	}

return(CurLen - 1);
}

etSeqBase						// returns a single base at Ofs 0..N relative to the SeqWrd ptd at by pSeq
CKangadna::GetBase(int Ofs,		// offset of base to return
				tSeqWrd4 *pSeq)		// Seq1 (probe) packed sequence containing base
{
int NumSeqWrdBases;
int Shf;

tSeqWrd4 *pSeq1 = pSeq;
tSeqWrd4 SeqWrd;
tSeqWrd4 BaseMsk;

// pSeq may pt to one of the sequence header words, skip over these
while(*pSeq1 & cSeqWrd4MSWHdr)
	{
	if((*pSeq1 & 0x0e0000001) == cSeqWrd4EOS)
		return(eBaseEOS);
	pSeq1 += 1;
	}
//pSeq now pts into sequence words containing bases 
while(!((SeqWrd = *pSeq1++) & cSeqWrd4MSWHdr)) // iterate over sequence words until either base located or another sequence header word
	{
	if(!(SeqWrd & cSeqWrd4PartSeq))				// if not a partial sequence word
		{
		if(Ofs >= 15)							// if requested base not in this sequence word then keep iterating words
			{
			Ofs -= 15;
			continue;
			}
		Shf = 2 * (14 - Ofs);
		SeqWrd >>= Shf;
		return(SeqWrd & 0x03);
		}

	// its a partial sequence word containing at most 14 bases
	if(Ofs >= 14)
    	return(eBaseEOS);

	// how many bases are in this partial sequence word?
	BaseMsk = 0x03;
	for(NumSeqWrdBases = 14; NumSeqWrdBases >= 0; NumSeqWrdBases--,BaseMsk <<= 2)
		if(!(BaseMsk & SeqWrd))
			break;

	if(Ofs >= NumSeqWrdBases)					// can't return a base which at an Ofs not within sequence
    	return(eBaseEOS);

	Shf = 2 * (14 - Ofs);
	SeqWrd >>= Shf;
	return(SeqWrd & 0x03);
	}

return(eBaseEOS);
}


int											// returned length of exactly matching initial (starting at ofs 0) substring between probe and target 
CKangadna::GetExactMatchLen(tSeqWrd4 *pSeqA,// Seq1 (probe) packed sequence
				tSeqWrd4 *pSeqB, 			// Seq2 (target) packed sequence
				int MaxLen)					// if non-zero then immediately return as soon as OverlapLen is at least this length - no need to exhustively find the maximal overlap
{
int NumSeq1Bases;
int NumSeq2Bases;
int MaxPartialCmpBases;
int OverlapLen;
if(pSeqA == NULL || pSeqB == NULL)
	return(0);

OverlapLen = 0;

tSeqWrd4 BaseMsk;
tSeqWrd4 *pSeq1 = pSeqA;
tSeqWrd4 *pSeq2 =  pSeqB;
tSeqWrd4 SeqWrd1;
tSeqWrd4 SeqWrd2;

while(*pSeq1 & cSeqWrd4MSWHdr)
	{
	if((*pSeq1 & 0x0e000001) == cSeqWrd4EOS)
		return(0);
	pSeq1 += 1;
	}
while(*pSeq2 & cSeqWrd4MSWHdr)
	{
	if((*pSeq1 & 0x0e0000001) == cSeqWrd4EOS)
		return(0);
	pSeq2 += 1;
	}
if(MaxLen < 0)
	MaxLen = 0;

while(!((SeqWrd1 = *pSeq1++) & cSeqWrd4MSWHdr || (SeqWrd2 = *pSeq2++) & cSeqWrd4MSWHdr))
	{
	if(!(SeqWrd1 & cSeqWrd4PartSeq) && !(SeqWrd2 & cSeqWrd4PartSeq)) // if not part words
		{
		if(SeqWrd1 == SeqWrd2)		
			{
			OverlapLen += 15;	// all bases match within word so can move onto next pair of words
			if(MaxLen && OverlapLen >= MaxLen)
				return(OverlapLen);
			continue;
			}
		}

	if(MaxLen && ((OverlapLen + 15) < MaxLen))
		return(OverlapLen);

	// either part words or there was a mismatch within words
	// determine number of bases in SeqWrd1
	if(SeqWrd1 & cSeqWrd4LSWHdr)		// SeqWrd1 will have bits 30 or 31 set if a part word
		{
		BaseMsk = 0x03;
		for(NumSeq1Bases = 14; NumSeq1Bases >= 0; NumSeq1Bases--,BaseMsk <<= 2)
			if(!(BaseMsk & SeqWrd1))
				break;
		}
	else
		NumSeq1Bases = 15;

	// determine number of bases in SeqWrd2
	if(SeqWrd2 & cSeqWrd4LSWHdr)		// SeqWrd2 will have bits 30 or 31 set if a part word
		{
		BaseMsk = 0x03;
		for(NumSeq2Bases = 14; NumSeq2Bases >= 0; NumSeq2Bases--,BaseMsk <<= 2)
			if(!(BaseMsk & SeqWrd2))
				break;
		}
	else
		NumSeq2Bases = 15;

	// do a minimal number of base compares 
	MaxPartialCmpBases = min(NumSeq1Bases,NumSeq2Bases);
	BaseMsk = 0x30000000;
	while(MaxPartialCmpBases--) {
		if((SeqWrd1 & BaseMsk) != (SeqWrd2 & BaseMsk)) 
			break;
		OverlapLen += 1;				// still matching base for base
		if(MaxLen && (OverlapLen >= MaxLen))
				return(OverlapLen);
		BaseMsk >>= 2;
		}
	break;
	}

return(OverlapLen);
}

bool										// true if Seq1 matches Seq2 for MatchLen
CKangadna::IsMatched(int ReqMatchLen,		// required match length
			int Seq1Ofs,					// base offset in Seq1 at which to to start match
			int Seq1Len,					// Seq1 is this length
			tSeqWrd4 *pSeq1,				// look for match between this sequence 1 and sequence 2
			int Seq2Ofs,					// base offset in Seq2  at which to to start match
			int Seq2Len,					// Seq 2 is this length
			tSeqWrd4 *pSeq2,				// sequence 2
			int MaxSubsPerK,				// allow upto this many substitutions per 1000 bp of overlap (expected to be in the range 0..50)
			int MaxEnd12Subs,				// allow at the initial 12bp of the 5' or 3' of match to have this many base mismatches in addition to the overall allowed MaxSubs (expected to be in range 0..6) 
			int *pNumSubs)					// total number of substitutions actually required
{
int Max3End12Subs;
int	KbpLen;
tSeqWrd4 *pStartSeq1;
tSeqWrd4 *pStartSeq2;
tSeqWrd4 CurSeq1Wrd;
tSeqWrd4 CurSeq2Wrd;
int AllowNumSubs;
int NumMMs;
int Num5MMs;
int Num3MMs;
int NumBasesChkd;

// some sanity checks
if(ReqMatchLen < 16 || (Seq1Ofs + ReqMatchLen) > Seq1Len || (Seq2Ofs + ReqMatchLen) > Seq2Len || 
   MaxSubsPerK < 0 || MaxSubsPerK > 50 ||
   MaxEnd12Subs < 0 || MaxEnd12Subs > 6 || 
   pSeq1 == NULL || pSeq2 == NULL)	
	return(0);

if(pNumSubs != NULL)
	*pNumSubs = 0;

Max3End12Subs = MaxEnd12Subs;
if(MaxEnd12Subs != 0 && ReqMatchLen < 24)		// if no end subs or minimum overlap would contain both 5'start and 3'ends then no need to adjust 3' endsubs
   Max3End12Subs = max(1,(MaxEnd12Subs * (ReqMatchLen - 11)) / 12);

KbpLen = ReqMatchLen - (2 * Max3End12Subs);
if(KbpLen > 20 && MaxSubsPerK)
	AllowNumSubs = max(1, (999 + (KbpLen * MaxSubsPerK)) / 1000);
else
	AllowNumSubs = 0;

// skip over any sequence header SeqWrds
while(*pSeq1 & cSeqWrd4MSWHdr)			
	{
	if((*pSeq1 & 0x0e000001) == cSeqWrd4EOS)
		return(false);
	pSeq1 += 1;
	}
while(*pSeq2 & cSeqWrd4MSWHdr)
	{
	if((*pSeq2 & 0x0e0000001) == cSeqWrd4EOS)
		return(false);
	pSeq2 += 1;
	}

// fetch and rotate into the MSbits the initial base
pStartSeq1 = &pSeq1[Seq1Ofs/15];
CurSeq1Wrd = *pStartSeq1++;
CurSeq1Wrd <<= 2 * (Seq1Ofs % 15); 
Seq1Ofs += 1;

pStartSeq2 = &pSeq2[Seq2Ofs/15];
CurSeq2Wrd = *pStartSeq2++;
CurSeq2Wrd <<= 2 * (Seq2Ofs % 15);
Seq2Ofs += 1;

NumBasesChkd = 0;
NumMMs = 0;
Num5MMs = 0;
Num3MMs = 0;
for(NumBasesChkd = 0; NumBasesChkd < ReqMatchLen; NumBasesChkd+=1, Seq1Ofs+=1, Seq2Ofs+=1)
	{
	if((CurSeq1Wrd & 0x030000000) != (CurSeq2Wrd  & 0x030000000))
		{
		if(MaxEnd12Subs && (NumBasesChkd < 12 || NumBasesChkd >= (ReqMatchLen - 12)))
			{
			if(NumBasesChkd < 12)
				{
				if(++Num5MMs > MaxEnd12Subs)
					break;
				}
			else
				{
				if(++Num3MMs > Max3End12Subs)
					break;
				}
			}
		else
			if(++NumMMs > AllowNumSubs)
				break;
		}

	if(!(Seq1Ofs%15))
		CurSeq1Wrd = *pStartSeq1++;
	else
		CurSeq1Wrd <<= 2;
	if(!(Seq2Ofs%15))
		CurSeq2Wrd = *pStartSeq2++;
	else
	    CurSeq2Wrd <<= 2;
	}
if(pNumSubs != NULL && NumBasesChkd == ReqMatchLen)
	*pNumSubs = NumMMs + Num5MMs + Num3MMs;
return(NumBasesChkd == ReqMatchLen);
}

// GetOverlapAB will process for overlaps in which SeqA (SeqB), completely or 5' overlaps, SeqB (SeqA)
int							// returned number of bases overlapped between SeqA and SeqB
CKangadna::GetOverlapAB(int MinOverlap,		// seq1 and seqB must overlap by at least this minimum overlap (expected to be at least 16bp)
			int *pSeqALeft,					// returned number of bases 5' before seq A is overlapped onto seq B (negative if seq B is instead overlapping onto seq A)
			int SeqALen,					// SeqA is this length
			tSeqWrd4 *pSeqA,				// look for overlap between this sequence A and sequence B
			int SeqBLen,					// Seq B is this length
			tSeqWrd4 *pSeqB,				// sequence B
			bool bNoSeqBOvrlapA,			// if true and no seq 1 onto seq B overlap then do not look for seq B overlap onto seq A
			int MaxSubsPerK,				// allow upto this many substitutions per 1000 bp of overlap (expected to be in the range 0..50)
			int MaxEnd12Subs,				// allow at the initial 12bp of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs (expected to be in range 0..6) 
			int *pNumSubs)					// number of substitutions actually required
{
int CurNumSubs;
int Max3End12Subs;
int Cur5EndSubs;
int Cur3EndSubs;
int ABNumSubs;
int SeqALeft;
int ABMaxOverlapLen;

int BANumSubs;
int SeqBLeft;
int BAMaxOverlapLen;

int FirstSeqAOfs;
int CurSeqAOfs;
int CurSeqBOfs;
tSeqWrd4 CurSeqABase;
tSeqWrd4 CurSeqBBase;
tSeqWrd4 *pStartSeqA;
tSeqWrd4 *pCurSeqA;
tSeqWrd4 *pCurSeqB;

tSeqWrd4 FirstSeqAWrd;
tSeqWrd4 FirstSeqBWrd;
tSeqWrd4 CurSeqAWrd;
tSeqWrd4 CurSeqBWrd;
tSeqWrd4 FirstSeqABase;
tSeqWrd4 FirstSeqBBase;

int KbpLen;
int CurMaxKbpNumSubs;

// some sanity checks
if(MinOverlap < 16 || MinOverlap > SeqALen || MinOverlap > SeqBLen || 
   MaxSubsPerK < 0 || MaxSubsPerK > 50 ||
   MaxEnd12Subs < 0 || MaxEnd12Subs > 6 || 
   pSeqALeft == NULL || pSeqA == NULL || pSeqB == NULL)	
	return(0);

*pSeqALeft = 0;
if(pNumSubs != NULL)
	*pNumSubs = 0;

Max3End12Subs = MaxEnd12Subs;
if(MaxEnd12Subs != 0 && MinOverlap < 24)		// if no end subs or minimum overlap would contain both 5'start and 3'ends then no need to adjust 3' endsubs
   Max3End12Subs = max(1,(MaxEnd12Subs * (MinOverlap - 11)) / 12);

// skip over any sequence header SeqWrds
while(*pSeqA & cSeqWrd4MSWHdr)			
	{
	if((*pSeqA & 0x0e000001) == cSeqWrd4EOS)
		return(0);
	pSeqA += 1;
	}
while(*pSeqB & cSeqWrd4MSWHdr)
	{
	if((*pSeqB & 0x0e0000001) == cSeqWrd4EOS)
		return(0);
	pSeqB += 1;
	}
pStartSeqA = pSeqA;

// initally look for overlaps of SeqA onto SeqB; then look for overlap of SeqB onto SeqA and choose the maximal with minimum subs as a tiebreaker
CurNumSubs = 0;
Cur5EndSubs = 0;
Cur3EndSubs = 0;
CurSeqAOfs = 0;

ABNumSubs = 0;
SeqALeft = 0;
ABMaxOverlapLen = 0;
FirstSeqBWrd = *pSeqB;
FirstSeqBBase = FirstSeqBWrd & 0x30000000; 
for(FirstSeqAOfs = 0; FirstSeqAOfs <= (SeqALen - MinOverlap); FirstSeqAOfs++)
	{
	if(!(FirstSeqAOfs % 15))
		FirstSeqAWrd = *pSeqA++ & cSeqWrd4Msk;
	else
		FirstSeqAWrd <<= 2;
	if((FirstSeqABase = (FirstSeqAWrd & 0x30000000)) != FirstSeqBBase && !(MaxEnd12Subs || MaxSubsPerK)) // if mismatches when no subs allowed then try next SeqA base
		continue;

	CurNumSubs = 0;
	Cur5EndSubs = 0;
	Cur3EndSubs = 0;
	CurSeqAWrd = FirstSeqAWrd;
	CurSeqBWrd = FirstSeqBWrd;
	CurSeqABase = FirstSeqABase;
	CurSeqBBase = FirstSeqBBase;
	CurSeqAOfs = FirstSeqAOfs;

	pCurSeqA = pSeqA;
	pCurSeqB = pSeqB;
	ABMaxOverlapLen = min(SeqALen - FirstSeqAOfs,SeqBLen);	// current max overlap possible of A onto B
	Max3End12Subs = MaxEnd12Subs;
	
	KbpLen = ABMaxOverlapLen - (2 * Max3End12Subs);
	if(KbpLen > 20 && MaxSubsPerK)
		CurMaxKbpNumSubs = max(1,(999 + (KbpLen * MaxSubsPerK))/1000); 
	else
		CurMaxKbpNumSubs = 0;
	if(MaxEnd12Subs > 0)
		{	
		if(ABMaxOverlapLen < 24)		// if no end subs or minimum overlap would contain both 5'start and 3'ends then no need to adjust 3' endsubs
			Max3End12Subs = max(1,(MaxEnd12Subs * (ABMaxOverlapLen - 11)) / 12);
		}

	for(CurSeqBOfs = 0; CurSeqBOfs < ABMaxOverlapLen; CurSeqBOfs++)
		{
		if(CurSeqBOfs > 0)
			{
			if(!(++CurSeqAOfs % 15))
				CurSeqAWrd = *pCurSeqA++ & cSeqWrd4Msk;
			else
				CurSeqAWrd <<= 2;
			CurSeqABase = CurSeqAWrd & 0x30000000;
			}

		if(!(CurSeqBOfs % 15))
			CurSeqBWrd = *pCurSeqB++ & cSeqWrd4Msk;
		else
			CurSeqBWrd <<= 2;
		CurSeqBBase = CurSeqBWrd & 0x30000000;

		if(CurSeqABase != CurSeqBBase)
			{
			if(MaxEnd12Subs)
				{
				if(CurSeqBOfs < 12)
					Cur5EndSubs += 1;
				else
					{
					if(CurSeqBOfs < (ABMaxOverlapLen - 12))
						CurNumSubs += 1;
					else
						Cur3EndSubs += 1;
					}
				}
			else
				CurNumSubs += 1;
			if(Cur5EndSubs > MaxEnd12Subs || Cur3EndSubs > Max3End12Subs || CurNumSubs > CurMaxKbpNumSubs)
				break;
			}
		}
	if(CurSeqBOfs == ABMaxOverlapLen)
		{
		ABNumSubs = Cur5EndSubs + Cur3EndSubs + CurNumSubs;
		SeqALeft = FirstSeqAOfs;
		break;
		}
	ABMaxOverlapLen = 0;
	}

if(bNoSeqBOvrlapA || ((ABMaxOverlapLen == min(SeqALen,SeqBLen)) && ABNumSubs <= 2)) // if just a few subs then accept even if may not be optimal
	{
	*pSeqALeft = SeqALeft;
	if(pNumSubs != NULL)
		*pNumSubs = ABNumSubs;

// double check!!!!!
    if(ABMaxOverlapLen > 0 &&
		ABMaxOverlapLen != min(SeqALen - SeqALeft, SeqBLen))
		printf("\nCheck this!!!!");
	return(ABMaxOverlapLen);
	}

//////////////////////////////////////////////
// now look for overlaps of SeqB onto SeqA

pSeqA = pSeqB;			// exchange the sequences and their respective lengths
CurNumSubs = SeqALen;
SeqALen = SeqBLen;
pSeqB = pStartSeqA;
SeqBLen = CurNumSubs;

CurNumSubs = 0;
Cur5EndSubs = 0;
Cur3EndSubs = 0;
CurSeqAOfs = 0;
BANumSubs = 0;
SeqBLeft = 0;
BAMaxOverlapLen = 0;
FirstSeqBWrd = *pSeqB;
FirstSeqBBase = FirstSeqBWrd & 0x30000000; 
for(FirstSeqAOfs = 0; FirstSeqAOfs <= (SeqALen - MinOverlap); FirstSeqAOfs++)
	{
	if(!(FirstSeqAOfs % 15))
		FirstSeqAWrd = *pSeqA++ & cSeqWrd4Msk;
	else
		FirstSeqAWrd <<= 2;
	if((FirstSeqABase = (FirstSeqAWrd & 0x30000000)) != FirstSeqBBase && !(MaxEnd12Subs || MaxSubsPerK)) // if mismatches when no subs allowed then try next SeqA base
		continue;

	CurNumSubs = 0;
	Cur5EndSubs = 0;
	Cur3EndSubs = 0;
	CurSeqAWrd = FirstSeqAWrd;
	CurSeqBWrd = FirstSeqBWrd;
	CurSeqABase = FirstSeqABase;
	CurSeqBBase = FirstSeqBBase;
	CurSeqAOfs = FirstSeqAOfs;

	pCurSeqA = pSeqA;
	pCurSeqB = pSeqB;
	BAMaxOverlapLen = min(SeqALen - FirstSeqAOfs,SeqBLen);	// sequences were exchanged so actually current max overlap possible of B onto A
	Max3End12Subs = MaxEnd12Subs;
	
	KbpLen = BAMaxOverlapLen - (2 * Max3End12Subs);
	if(KbpLen > 20 && MaxSubsPerK)
		CurMaxKbpNumSubs = max(1,(999 + (KbpLen * MaxSubsPerK))/1000); 
	else
		CurMaxKbpNumSubs = 0;
	if(MaxEnd12Subs > 0)
		{	
		if(BAMaxOverlapLen < 24)		// if no end subs or minimum overlap would contain both 5'start and 3'ends then no need to adjust 3' endsubs
			Max3End12Subs = max(1,(MaxEnd12Subs * (BAMaxOverlapLen - 11)) / 12);
		}

	for(CurSeqBOfs = 0; CurSeqBOfs < BAMaxOverlapLen; CurSeqBOfs++)
		{
		if(CurSeqBOfs > 0)
			{
			if(!(++CurSeqAOfs % 15))
				CurSeqAWrd = *pCurSeqA++ & cSeqWrd4Msk;
			else
				CurSeqAWrd <<= 2;
			CurSeqABase = CurSeqAWrd & 0x30000000;
			}

		if(!(CurSeqBOfs % 15))
			CurSeqBWrd = *pCurSeqB++ & cSeqWrd4Msk;
		else
			CurSeqBWrd <<= 2;
		CurSeqBBase = CurSeqBWrd & 0x30000000;

		if(CurSeqABase != CurSeqBBase)
			{
			if(MaxEnd12Subs)
				{
				if(CurSeqBOfs < 12)
					Cur5EndSubs += 1;
				else
					{
					if(CurSeqBOfs < (BAMaxOverlapLen - 12))
						CurNumSubs += 1;
					else
						Cur3EndSubs += 1;
					}
				}
			else
				CurNumSubs += 1;
			if(Cur5EndSubs > MaxEnd12Subs || Cur3EndSubs > Max3End12Subs || CurNumSubs > CurMaxKbpNumSubs)
				break;
			}
		}
	if(CurSeqBOfs == BAMaxOverlapLen)
		{
		BANumSubs = Cur5EndSubs + Cur3EndSubs + CurNumSubs;
		SeqBLeft = FirstSeqAOfs;
		break;
		}
	BAMaxOverlapLen = 0;
	}

if(!(ABMaxOverlapLen +  BAMaxOverlapLen))
	{
	*pSeqALeft = 0;
	if(pNumSubs != NULL)
		*pNumSubs = 0;
	return(0);
	}
if(ABMaxOverlapLen > BAMaxOverlapLen ||
	(ABMaxOverlapLen == BAMaxOverlapLen && 	ABNumSubs <= BANumSubs))
	{
	*pSeqALeft = SeqALeft;
	if(pNumSubs != NULL)
		*pNumSubs = ABNumSubs;
	return(ABMaxOverlapLen);
	}

*pSeqALeft = -1 * SeqBLeft;
if(pNumSubs != NULL)
	*pNumSubs = BANumSubs;
return(BAMaxOverlapLen);
}


tSeqID					// returned sequence identifer
CKangadna::SfxIdx2SeqID(UINT64 SfxWrdIdx) // index + 1
{
tSeqWrd4 *pTarg;
if((pTarg = SfxIdxToFirstSeqWrd(SfxWrdIdx))==NULL)
	return(0);
return(GetPackedSeqID(pTarg));
}


tSeqWrd4 *									// returned ptr to 1st packed sequence word
CKangadna::SfxIdxToFirstSeqWrd(UINT64 SfxWrdIdx, // index + 1
					   int *pWrdOfs)	  // returned sequence word is at this relative word offset to SfxWrdIdx; e.g if 1 then SfxWrdIdx referenced the second word in the sequence
{
int RelWrdOfs;
UINT64 TargIdx;
UINT64 TargEl;
UINT8 *pSfxEls;
if(SfxWrdIdx == 0 || SfxWrdIdx > (UINT64)m_Sequences.NumSuffixEls)
	return(NULL);

TargEl = (SfxWrdIdx - 1) * m_Sequences.SfxElSize;
pSfxEls = (UINT8 *)m_Sequences.pSuffixArray;

if(m_Sequences.SfxElSize == 4)
	TargIdx = *(UINT32 *)&pSfxEls[TargEl];
else
	TargIdx = Unpack5(&pSfxEls[TargEl]);

RelWrdOfs = 0;

tSeqWrd4 *pTarg4;
pTarg4 =  &((tSeqWrd4 *)m_Sequences.pSeqs2Assemb)[TargIdx];
while(!(pTarg4[-1] & cSeqWrd4MSWHdr))
	{
	pTarg4 -= 1;
	RelWrdOfs += 1;
	}
if(pWrdOfs != NULL)
	*pWrdOfs = RelWrdOfs;
return(pTarg4);
}


int														// number of sequence flags processed	
CKangadna::GetSeqFlags(tSeqID StartSeqID,				// starting from this sequence identifier inclusive
			tSeqID EndSeqID,							// and ending at this sequence identifier inclusive
			tsMultiSeqFlags *pMultiSeqFlags,			// holds sequence identifiers plus used to return flags
			bool bSerialise)							// set true if access to flags are required to be serialised
{
tSeqID CurSeqID;
if(StartSeqID > EndSeqID || pMultiSeqFlags == NULL)
	return(0);
if(bSerialise)
	AcquireSerialiseSeqFlags();
for(CurSeqID = StartSeqID; CurSeqID <= EndSeqID; CurSeqID++,pMultiSeqFlags++)
	{
	pMultiSeqFlags->SeqID = CurSeqID;
	pMultiSeqFlags->ResetFlags = 0;
	pMultiSeqFlags->SetFlags = 0;
	pMultiSeqFlags->CurFlags = m_Sequences.pSeqFlags[pMultiSeqFlags->SeqID-1];
	}
if(bSerialise)
	ReleaseSerialiseSeqFlags();
return((EndSeqID - StartSeqID) + 1);
}

// GetSeqFlags
// Batch fetch sequence flags for multiple sequences
int
CKangadna::GetSeqFlags(	UINT32 NumSeqs,								// number of sequences 
			tsMultiSeqFlags *pMultiSeqFlags,			// holds sequence identifiers plus used to return flags
			bool bSerialise)							// set true if access to flags are required to be serialised
{
UINT32 NumRead;
if(bSerialise)
	AcquireSerialiseSeqFlags();
for(NumRead = 0; NumRead < NumSeqs; NumRead++, pMultiSeqFlags++)
	{
	if(pMultiSeqFlags->SeqID < 1)
		{
		if(bSerialise)
			ReleaseSerialiseSeqFlags();
		return(0);
		}
	pMultiSeqFlags->CurFlags = m_Sequences.pSeqFlags[pMultiSeqFlags->SeqID-1];
	}
if(bSerialise)
	ReleaseSerialiseSeqFlags();
return(NumRead);
}


int
CKangadna::GetSeqFlags(tSeqID SeqID,			// sequence identifier (32 bits)
			bool bSerialise)	// set true if access to headers are required to be serialised
{
int SeqFlags;
if(bSerialise)
	AcquireSerialiseSeqFlags();
SeqFlags = m_Sequences.pSeqFlags[SeqID-1];
if(bSerialise)
	ReleaseSerialiseSeqFlags();
return(SeqFlags);
}



// UpdateSeqFlags
// Batch update sequence flags for multiple sequences
// Sequence flags are reset with ResetFlags before being set with SetFlags
int
CKangadna::UpdateSeqFlags(UINT32 NumSeqs,				// number of sequences requiring flag updates
			tsMultiSeqFlags *pMultiSeqFlags,			// holds sequence identifiers plus flag update operations
			bool bSerialise)							// set true if access to flags are required to be serialised
{
UINT32 NumUpdated;
UINT32 SeqFlags;
UINT32 PrevSeqFlags;
if(bSerialise)
	AcquireSerialiseSeqFlags();

for(NumUpdated = 0; NumUpdated < NumSeqs; NumUpdated++, pMultiSeqFlags++)
	{
	if(pMultiSeqFlags->SeqID < 1)
		{
		if(bSerialise)
			ReleaseSerialiseSeqFlags();
		return(0);
		}
	PrevSeqFlags = SeqFlags = m_Sequences.pSeqFlags[pMultiSeqFlags->SeqID-1];
	SeqFlags &= ~pMultiSeqFlags->ResetFlags;

	if((pMultiSeqFlags->SetFlags & cFlgOvlLenMsk) < (UINT16)(SeqFlags & cFlgOvlLenMsk))
		SeqFlags |= (pMultiSeqFlags->SetFlags & ~cFlgOvlLenMsk);			// retaining existing cFlgOvlLenMsk
	else
		{
		SeqFlags &= ~cFlgOvlLenMsk;			// replacing existing cFlgOvlLenMsk
		SeqFlags |= pMultiSeqFlags->SetFlags;								
		}
	
	if(SeqFlags != PrevSeqFlags)
		m_Sequences.pSeqFlags[pMultiSeqFlags->SeqID-1] = SeqFlags;
	pMultiSeqFlags->CurFlags = SeqFlags;
	}
if(bSerialise)
	ReleaseSerialiseSeqFlags();
return(NumUpdated);
}

// UpdateSeqFlags
// Update sequence flags for a single sequence
// Sequence flags are reset with ResetFlags before being set with SetFlags
int
CKangadna::UpdateSeqFlags(tSeqID SeqID,	// sequence identifier (32 bits)
			UINT32 SetFlags,			// flags to set, any overlap len in cFlgOvlLenMsk is only copied into flags if > existing overlap len
			UINT32 ResetFlags,			// flags to be reset
			bool bSerialise)			// set true if access to flags are required to be serialised
{
UINT32 SeqFlags;
UINT32 PrevSeqFlags;

SetFlags   &= 0x0000ffff;
ResetFlags &= 0x0000ffff;

if(bSerialise)
	AcquireSerialiseSeqFlags();
PrevSeqFlags = SeqFlags = m_Sequences.pSeqFlags[SeqID-1];

SeqFlags &= ~ResetFlags;
if((SetFlags & cFlgOvlLenMsk) < (SeqFlags & cFlgOvlLenMsk))
	SeqFlags |= (SetFlags & ~cFlgOvlLenMsk);			// retaining existing cFlgOvlLenMsk
else
	{
	SeqFlags &= ~cFlgOvlLenMsk;			// replacing existing cFlgOvlLenMsk
	SeqFlags |= SetFlags;								
	}

if(SeqFlags != PrevSeqFlags)
	m_Sequences.pSeqFlags[SeqID-1] = SeqFlags;
if(bSerialise)
	ReleaseSerialiseSeqFlags();
return(SeqFlags);
}

// Sequence flags are reset with ResetFlags before being set with SetFlags
int										// returned updated flags, < 0 if errors
CKangadna::UpdateSeqHeaderFlags(tSeqID SeqID,	// sequence identifier (32 bits)
					UINT32 SetFlags,			// flags to set, any overlap len in cFlgOvlLenMsk is only copied into flags if > existing overlap len
					UINT32 ResetFlags,			// flags to be reset
					bool bSerialise)			// set true if access to headers are required to be serialised
{
UINT64 *pSeqStart;
tSeqWrd4 *pStartSeqWrd;

if(SeqID < 1 || SeqID > (tSeqID)m_Sequences.NumSeqs2Assemb)
	return(-1);
pSeqStart = &m_Sequences.pSeqStarts[SeqID-1];

tSeqWrd4 *pSeqWrd;
pSeqWrd = (tSeqWrd4 *)m_Sequences.pSeqs2Assemb;
pStartSeqWrd = &pSeqWrd[*pSeqStart];

SetFlags   &= 0x0000ffff;
ResetFlags &= 0x0000ffff;

return(UpdateSeqHeaderFlags(pStartSeqWrd,SetFlags,ResetFlags,bSerialise));
}

// Sequence flags are reset with ResetFlags before being set with SetFlags
int					// < 0 if errors
CKangadna::UpdateAllSeqHeaderFlags(UINT32 SetFlags,	// flags to set, any overlap len in cFlgOvlLenMsk is only copied into flags if > existing overlap len
					UINT32 ResetFlags,				// flags to be reset
					bool bSerialise)				// set true if access to headers are required to be serialised
{
int Rslt;
tSeqID SeqID;
UINT32 Flags;
tSeqWrd4 *pSeqHdr;

Rslt = 0;
pSeqHdr = NULL;

SetFlags   &= 0x0000ffff;
ResetFlags &= 0x0000ffff;

while((pSeqHdr = IterSeqHeaders(pSeqHdr,&SeqID,NULL,&Flags,NULL,bSerialise))!=NULL)
	{
	if((Flags & ResetFlags) || ((Flags & SetFlags) != SetFlags))
		if((Rslt = UpdateSeqHeaderFlags(pSeqHdr,SetFlags,ResetFlags,bSerialise)) < 0)
			break;
	}

return(Rslt);
}

// Sequence flags are reset with ResetFlags before being set with SetFlags
int									// returned updated flags, or if < 0 then error 
CKangadna::UpdateSeqHeaderFlags(tSeqWrd4 *pSeqWrd,				// pts to a SeqWrd within the sequence
			UINT32 SetFlags,			// flags to set, any overlap len in cFlgOvlLenMsk is only copied into flags if > existing overlap len
			UINT32 ResetFlags,			// flags to be reset
			bool bSerialise)			// set true if access to headers are required to be serialised
{
tSeqWrd4 *pSeqHdrWrd;
tSeqWrd4 SeqWrd;

SetFlags   &= 0x0000ffff;
ResetFlags &= 0x0000ffff;

// get existing flags from sequence header
pSeqHdrWrd = (tSeqWrd4 *)pSeqWrd;
if(*pSeqHdrWrd == cSeqWrd4BOS || *pSeqHdrWrd == cSeqWrd4EOS) // ensure not trying to get header for BOS or EOS
	return(-1);
while((*pSeqHdrWrd & cSeqWrd4LSWHdr) != cSeqWrd4MSWHdr)	// backup until the most significant header word
	pSeqHdrWrd -= 1;
if(bSerialise)
	AcquireSerialiseSeqHdr();
SeqWrd = *pSeqHdrWrd;
SeqWrd &= ~ResetFlags;
if((SetFlags & cFlgOvlLenMsk) < (SeqWrd & cFlgOvlLenMsk))
	SeqWrd |= (SetFlags & ~cFlgOvlLenMsk);			// retaining existing cFlgOvlLenMsk
else
	{
	SeqWrd &= ~cFlgOvlLenMsk;			// replacing existing cFlgOvlLenMsk
	SeqWrd |= SetFlags;								
	}

*pSeqHdrWrd = SeqWrd;

if(bSerialise)
	ReleaseSerialiseSeqHdr();

return((int)(SeqWrd & 0x0000ffff));
}



// RemoveMarkedSeqs
// Removes all sequences marked for removal
//
int
CKangadna::RemoveMarkedSeqs(UINT32 RemovalFlags,		// if any of these flags set then remove this sequence
								UINT32 AllRequiredFlags,    // if containing any flags then remove this sequence unless all of these flags are set
								UINT32 AllOptionalFlags,    // if containing any flags then remove this sequence unless at least one of these flags is set
								bool bUpdateHdrFlags)		// if true, and if pSeqFlags != NULL, then replace sequence header flags with those from pSeqFlags
{
bool bRemove;
UINT32 NumRemoved;
tSeqID NewSeqID;
tSeqID CurSeqID;
UINT32 CurFlags;
UINT32 SrcFileID;
UINT32 SeqLen;
UINT32 NumSeqs2Assemb;
UINT64 Seqs2AssembLen;
UINT64 Seqs2AssembOfs;


UINT32 MinSeqLen;					// minimum sequence length
UINT32 MaxSeqLen;					// maximum sequence length

if(m_Sequences.pSeqs2Assemb == NULL || m_Sequences.AllocMemSeqs2Assemb == 0 || m_Sequences.NumSeqs2Assemb == 0)
	{
	if(m_Sequences.pSeqs2Assemb == NULL)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"RemoveMarkedSeqs, pSeqs2Assemb is NULL");
	if(m_Sequences.AllocMemSeqs2Assemb == 0)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"RemoveMarkedSeqs, AllocMemSeqs2Assemb is 0");
	if(m_Sequences.NumSeqs2Assemb == 0)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"RemoveMarkedSeqs, NumSeqs2Assemb is 0");
	return(eBSFerrParams);
	}

if(bUpdateHdrFlags && m_Sequences.pSeqFlags)
	{
	UINT16 *pSeqFlags;
	tSeqID SeqID;
	pSeqFlags = m_Sequences.pSeqFlags;
	for(SeqID = 1; SeqID <= m_Sequences.NumSeqs2Assemb; SeqID++,pSeqFlags++)
		UpdateSeqHeaderFlags(SeqID,*pSeqFlags,0,false);
	}

bRemove = false;
MinSeqLen = 0;
MaxSeqLen = 0;
NumSeqs2Assemb = 0;
Seqs2AssembLen = 0;
Seqs2AssembOfs = 1;
NewSeqID = 0;

tSeqWrd4 SeqWrd;
tSeqWrd4 *pSeqWrd = (tSeqWrd4 *)m_Sequences.pSeqs2Assemb;
tSeqWrd4 *pCopyTo = NULL;
pSeqWrd += 1;			// skip the BOS marker word
while((SeqWrd = *pSeqWrd++) != cSeqWrd4EOS)
	{
	if((SeqWrd & cSeqWrd4LSWHdr) == cSeqWrd4MSWHdr)
		{
		GetSeqHeader(pSeqWrd-1,&CurSeqID,&SrcFileID,&CurFlags,&SeqLen,false);
		bRemove = (CurFlags & RemovalFlags) ? true : false;
		if(!bRemove && AllRequiredFlags)
			{
			if((CurFlags & AllRequiredFlags) != AllRequiredFlags)
				bRemove = true;
			}

		if(!bRemove && AllOptionalFlags)
			{
			if(!(CurFlags & AllOptionalFlags))
				bRemove = true;
			}

		if(!bRemove && AllOptionalFlags)
			{
			// is trimming required?
			if((CurFlags & (cFlg5Prime | cFlg3Prime)) != (cFlg5Prime | cFlg3Prime))
				{
				int TrimTo;
				int TrimBy;
				int NumSeqWrds;

				pSeqWrd -= 1;
				if(pCopyTo == NULL)
					pCopyTo = pSeqWrd;
				TrimTo = ((SeqLen * (CurFlags & cFlgOvlLenMsk) >> 9)) / 100;
				TrimBy = SeqLen - TrimTo;
				NumSeqWrds = 3 + ((TrimTo+14)/15);

				if(MinSeqLen == 0 || MinSeqLen > (UINT32)TrimTo) // note min and max sequence accepted
					MinSeqLen = TrimTo;
				if(MaxSeqLen < (UINT32)TrimTo)
					MaxSeqLen = TrimTo;
				Seqs2AssembLen += TrimTo;
				NumSeqs2Assemb += 1;
				NewSeqID += 1;				

				SetSeqHeader(pSeqWrd,NewSeqID,SrcFileID,CurFlags | cFlg5Prime | cFlg3Prime,TrimTo,NULL,false,false);
				if(CurFlags & cFlg5Prime)		// overlap was on the 5' end of sequence so will trim back the 3' end
					NHSeqTrim(0,TrimBy,SeqLen,&pSeqWrd[3],true);
				else                            // overlap was on the 3' end of sequence so will trim back the 5' end
					NHSeqTrim(TrimBy,0,SeqLen,&pSeqWrd[3],true);

				while(NumSeqWrds--)
					{
					Seqs2AssembOfs += 1;
					*pCopyTo++ = *pSeqWrd++;
					}  
	
				while(((SeqWrd = *pSeqWrd) != cSeqWrd4EOS) && ((SeqWrd & cSeqWrd4LSWHdr) != cSeqWrd4MSWHdr)) // slough remainder of existing sequence
					pSeqWrd++;	
				continue;
				}
			}

		if(bRemove)
			{
			if(pCopyTo == NULL)
				pCopyTo = pSeqWrd - 1;
			}
		else
			{
			NumSeqs2Assemb += 1;		// accepting this sequence
			Seqs2AssembLen += SeqLen;
			if(MinSeqLen == 0 || MinSeqLen > SeqLen) // note min and max sequence accepted
					MinSeqLen = SeqLen;
			if(MaxSeqLen < SeqLen)
					MaxSeqLen = SeqLen;

			NewSeqID += 1;				// if copying over top of deleted sequences then give it a new sequence identifier
			if(pCopyTo != NULL)
				{
				SetSeqHeader(pSeqWrd-1,NewSeqID,SrcFileID,CurFlags,SeqLen,NULL,false,false);
				SeqWrd = pSeqWrd[-1];
				}	
			}
		}

	if(!bRemove)
		{
		Seqs2AssembOfs += 1;
		if(pCopyTo != NULL)
			*pCopyTo++ = SeqWrd;
		}
	}
if(pCopyTo != NULL)
	*pCopyTo = cSeqWrd4EOS;


NumRemoved = m_Sequences.NumSeqs2Assemb - NumSeqs2Assemb;
m_Sequences.NumSeqs2Assemb = NumSeqs2Assemb;
m_Sequences.Seqs2AssembLen = Seqs2AssembLen;
m_Sequences.MinSeqLen = MinSeqLen;
m_Sequences.MaxSeqLen = MaxSeqLen;
if(NumSeqs2Assemb > 0)
	m_Sequences.MeanSeqLen = (double)Seqs2AssembLen / NumSeqs2Assemb;
else
	m_Sequences.MeanSeqLen = 0.0;
m_Sequences.Seqs2AssembOfs = Seqs2AssembOfs;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"RemoveMarkedSeqs: Removed %u sequences, %u remaining",NumRemoved,m_Sequences.NumSeqs2Assemb);
if(m_Sequences.NumSeqs2Assemb == 0)
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Generated preprocessed reads file will contain 0 sequences");
return(eBSFSuccess);
}




int
CKangadna::SetAssemblyPars(bool bStrand,			// true if read strand specific assembly
				  int MinOverlap,				// minimum flank overlap length
  				  int MeanInsert,				// if paired end processing - mean insert size (defaults to 0 for single ended assembly)
				  int MinInsert,				// if paired end processing - paired end minimum insert size
				  int MaxInsert)				// if paired end processing - paired end maximum insert size
{
m_bAssembStrand = bStrand;
m_AssembMinOverlap = MinOverlap;
m_AssembMeanInsert = MeanInsert;
m_AssembMinInsert = MinInsert;
m_AssembMaxInsert = MaxInsert;
return(eBSFSuccess);
}




bool
CKangadna::SyncAllThreadsCompleted(bool bResetNextProcSeqID)
{
int ThreadsProcessing;
AcquireLock(true);
if(m_ThreadsProcessing > 0)
	m_ThreadsProcessing -= 1;
ThreadsProcessing = m_ThreadsProcessing;
if(ThreadsProcessing == 0 && bResetNextProcSeqID)
	m_NextProcSeqID = 0;
ReleaseLock(true);
if(!ThreadsProcessing)
	return(true);
do {
#ifdef _WIN32
	Sleep(2000);
#else
	sleep(2);
#endif
	AcquireLock(false);
	ThreadsProcessing = m_ThreadsProcessing;
	ReleaseLock(false);
	}
while(ThreadsProcessing);
return(true);
}


int
CKangadna::GetSeqProcRange(tSeqID *pStartingSeqID,	// where to return starting sequence identifier (inclusive)
								tSeqID *pEndingSeqID,	// where to return ending sequence identifier (inclusive)
								int MaxReq)				// at most this many identifiers limited to no more than m_NumProcSeqIDs  (0 if MaxReq is m_NumProcSeqIDs)
{
*pStartingSeqID = 0;
*pEndingSeqID = 0;

AcquireLock(true);
if(m_NextProcSeqID > m_FinalProcSeqID)
	{
	ReleaseLock(true);
	return(0);
	}
if(m_NumProcSeqIDs == 0)
	{
	if(m_ThreadsProcessing == 1)
		m_NumProcSeqIDs = cMaxMultiSeqFlags;
	else
		{
		m_NumProcSeqIDs = min(cMaxMultiSeqFlags,(m_FinalProcSeqID + m_ThreadsProcessing - 1)/(m_ThreadsProcessing * 100));
		if(m_NumProcSeqIDs < 64)
			m_NumProcSeqIDs = 64;
		}
	}
if(MaxReq == 0 || MaxReq > m_NumProcSeqIDs)
	MaxReq = m_NumProcSeqIDs;

MaxReq = (MaxReq + 1) & 0xffffe;		// ensure that if processing paired end reads then will always return a multiple of two so PE1 and PE2 will be processed by same thread
if(m_NextProcSeqID == 0)
	m_NextProcSeqID = m_StartProcSeqID;
*pStartingSeqID = m_NextProcSeqID;
m_NextProcSeqID += MaxReq;
if(m_NextProcSeqID > m_FinalProcSeqID)
	m_NextProcSeqID = m_FinalProcSeqID + 1;
*pEndingSeqID = m_NextProcSeqID - 1;
MaxReq = (int)(m_NextProcSeqID - *pStartingSeqID);
ReleaseLock(true);
return(MaxReq);
} 

// atomic merge of a source PE1 and a PE1 dest sequence with merge of a source PE2 and a PE2 dest sequence
int											// total number of bases in both pPE1DstSeqWrd and pPE2DstSeqWrd after merging or 0 if unable to merge because flags cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt already set for sequence Seq1ID or Seq2ID
CKangadna::AtomicSeqMerge(int PE1RelOverlap,	// if >= 0 then pPE1SrcSeqWrd sequence starts at this pPE1DstSeqWrd[PE1RelOverlap]; if < 0 then pPE1DstSeqWrd starts at this pPE1SrcSeqWrd[PE1RelOverlap]; PE1RelOverlap is in bases
					int PE1SrcSeqLen,			// number of bases in pPE1SrcSeqWrd sequence
				    tSeqWrd4 *pPE1SrcSeqWrd,	// merge this sequence with pPE1DstSeqWrd

					int PE2RelOverlap,			// if >= 0 then pPE2SrcSeqWrd starts at this pPE2DstSeqWrd[PE2RelOverlap]; if < 0 then pPE2DstSeqWrd starts at this pPE2SrcSeqWrd[PE2RelOverlap]; PE2RelOverlap is in bases
					int PE2SrcSeqLen,			// number of bases in pPE2SrcSeqWrd sequence
				    tSeqWrd4 *pPE2SrcSeqWrd,	// merge this sequence with pPE2DstSeqWrd				   

					int PE1DstSeqLen,			// number of bases currently in pPE1DstSeqWrd sequence
					int *pPE1DstSeqLen,			// returned length after merging
				   tSeqWrd4 *pPE1DstSeqWrd,		// merge this sequence with pPE1SrcSeqWrd; merged sequence replaces this sequence
				   tSeqWrd4 *pPE1TmpSeqWrd,		// temp sequence space for use as may be needed whilst merging
					
					int PE2DstSeqLen,			// number of bases currently in pPE2DstSeqWrd sequence
					int *pPE2DstSeqLen,			// returned length after merging
				   tSeqWrd4 *pPE2DstSeqWrd,		// merge this sequence with pPE2SrcSeqWrd; merged sequence replaces this sequence
				   tSeqWrd4 *pPE2TmpSeqWrd,		// temp sequence space for use as may be needed whilst merging

				   tSeqID Seq1ID,				// mandatory sequence identifier to test/set flags cFlgAsmbExtn | cFlgAsmbCplt
				   tSeqID Seq2ID)				// optional sequence identifier to test/set flags cFlgAsmbExtn | cFlgAsmbCplt
{
int MergePE1DstSeqLen;
int MergePE2DstSeqLen;

AcquireSerialiseSeqFlags();
if((m_Sequences.pSeqFlags[Seq1ID-1] & (cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt)) || // if already claimed or a seed then can't merge
   (Seq2ID != 0 && (m_Sequences.pSeqFlags[Seq2ID-1] & (cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt))))
	{
	ReleaseSerialiseSeqFlags();
	return(0);
	}
m_Sequences.pSeqFlags[Seq1ID-1] |= (cFlgAsmbExtn | cFlgAsmbCplt);
if(Seq2ID != 0)
	m_Sequences.pSeqFlags[Seq2ID-1] |= (cFlgAsmbExtn | cFlgAsmbCplt);
ReleaseSerialiseSeqFlags();

MergePE1DstSeqLen = SeqMerge(PE1RelOverlap,PE1SrcSeqLen,pPE1SrcSeqWrd,PE1DstSeqLen,pPE1DstSeqWrd,pPE1TmpSeqWrd,(char *)"Called from AtomicSeqMerge D");
MergePE2DstSeqLen = SeqMerge(PE2RelOverlap,PE2SrcSeqLen,pPE2SrcSeqWrd,PE2DstSeqLen,pPE2DstSeqWrd,pPE2TmpSeqWrd,(char *)"Called from AtomicSeqMerge E");

*pPE1DstSeqLen = MergePE1DstSeqLen;
*pPE2DstSeqLen = MergePE2DstSeqLen;

return(MergePE1DstSeqLen + MergePE2DstSeqLen);
}



// atomic merge of a source PE1 and PE2 with a SE dest sequence
int												// number of bases in sequence pDstSeqWrd after merging or 0 if unable to merge because flags cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt already set for sequence Seq1ID or Seq2ID
CKangadna::AtomicSeqMerge(int PE1RelOverlap,	// if >= 0 then PE1 sequence starts at this pDstSeq[PE1RelOverlap]; if < 0 then pDstSeq starts at this PE1 sequence[PE1RelOverlap]; PE1RelOverlap is in bases
					int PE1SrcSeqLen,			// number of bases in pPE1SrcSeqWrd sequence
				    tSeqWrd4 *pPE1SrcSeqWrd,	// merge this sequence with pDstSeqWrd
					int PE2RelOverlap,			// if >= 0 then pPESrcSeqWrd starts at this pDstSeq[PE2RelOverlap]; if < 0 then pDstSeq starts at this PE2 sequence[PE2RelOverlap]; PE2RelOverlap is in bases
					int PE2SrcSeqLen,			// number of bases in pPE2SrcSeqWrd sequence
				    tSeqWrd4 *pPE2SrcSeqWrd,	// merge this sequence with pDstSeqWrd				   
				    int DstSeqLen,				// number of bases currently in pDstSeqWrd sequence
				   tSeqWrd4 *pDstSeqWrd,		// merge this sequence with pPE1SrcSeqWrd and pPE2SrcSeqWrd; merged sequence replaces this sequence
				   tSeqWrd4 *pTmpSeqWrd,		// temp sequence space for use as may be needed whilst merging
				   tSeqID Seq1ID,				// mandatory sequence identifier to test/set flags cFlgAsmbExtn | cFlgAsmbCplt
					tSeqID Seq2ID)				// optional sequence identifier to test/set cFlgAsmbExtn | cFlgAsmbCplt
{
AcquireSerialiseSeqFlags();
if((m_Sequences.pSeqFlags[Seq1ID-1] & (cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt)) || // if already claimed or a seed then can't merge
   (Seq2ID != 0 && (m_Sequences.pSeqFlags[Seq2ID-1] & (cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt))))
	{
	ReleaseSerialiseSeqFlags();
	return(0);
	}
m_Sequences.pSeqFlags[Seq1ID-1] |= (cFlgAsmbExtn | cFlgAsmbCplt);
if(Seq2ID != 0)
	m_Sequences.pSeqFlags[Seq2ID-1] |= (cFlgAsmbExtn | cFlgAsmbCplt);
ReleaseSerialiseSeqFlags();

int PE1MergeLen;
int PE2MergeLen;
if(pTmpSeqWrd != NULL)
	*pTmpSeqWrd = 0;
PE1MergeLen = SeqMerge(PE1RelOverlap,PE1SrcSeqLen,pPE1SrcSeqWrd,DstSeqLen,pDstSeqWrd,pTmpSeqWrd,(char *)"Called from AtomicSeqMerge F");
if(pTmpSeqWrd != NULL)
	*pTmpSeqWrd = 1;
if(PE1RelOverlap < 0)
	PE2RelOverlap += abs(PE1RelOverlap);
PE2MergeLen = SeqMerge(PE2RelOverlap,PE2SrcSeqLen,pPE2SrcSeqWrd,PE1MergeLen,pDstSeqWrd,pTmpSeqWrd,(char *)"Called from AtomicSeqMerge G");
if(pTmpSeqWrd != NULL)
	*pTmpSeqWrd = 2;


return(PE2MergeLen);
}


int												// number of bases in sequence pDstSeqWrd after merging or 0 if unable to merge because flags cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt already set for sequence Seq1ID
CKangadna::AtomicSeqMerge(int RelOverlap,		// if >= 0 then pSrcSeq starts at this pDstSeq[RelOverlap]; if < 0 then pDstSeq starts at this pSrcSeq[RelOverlap]; RelOverlap is in bases
					int SrcSeqLen,				// number of bases in pSrcSeqWrd sequence
				   tSeqWrd4 *pSrcSeqWrd,		// merge this sequence with pDstSeqWrd
					int DstSeqLen,				// number of bases currently in pDstSeqWrd sequence
				   tSeqWrd4 *pDstSeqWrd,		// merge this sequence with pSrcSeqWrd; merged sequence replaces this sequence
				   tSeqWrd4 *pTmpSeqWrd,		// temp sequence space for use as may be needed whilst merging
				   tSeqID SeqID)				// mandatory sequence identifier to test/set flags cFlgAsmbExtn | cFlgAsmbCplt
{
int Rslt;
AcquireSerialiseSeqFlags();
if((m_Sequences.pSeqFlags[SeqID-1] & (cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt))) // if already claimed then can't merge
	{
	ReleaseSerialiseSeqFlags();
	return(0);
	}
m_Sequences.pSeqFlags[SeqID-1] |= (cFlgAsmbExtn | cFlgAsmbCplt);
ReleaseSerialiseSeqFlags();

Rslt = SeqMerge(RelOverlap,SrcSeqLen,pSrcSeqWrd,DstSeqLen,pDstSeqWrd,pTmpSeqWrd,(char *)"Called from AtomicSeqMerge H");
return(Rslt);
}

// merge of a source PE1 and PE2 with a SE dest sequence
int												// number of bases in sequence pDstSeqWrd after merging
CKangadna::SeqMergePE12ToSE(int PE1RelOverlap,			// if >= 0 then PE1 sequence starts at this pDstSeq[PE1RelOverlap]; if < 0 then pDstSeq starts at this PE1 sequence[PE1RelOverlap]; PE1RelOverlap is in bases
					int PE1SrcSeqLen,			// number of bases in pPE1SrcSeqWrd sequence
				   tSeqWrd4 *pPE1SrcSeqWrd,		// merge this sequence with pDstSeqWrd
					int PE2RelOverlap,			// if >= 0 then pPESrcSeqWrd starts at this pDstSeq[PE2RelOverlap]; if < 0 then pDstSeq starts at this PE2 sequence[PE2RelOverlap]; PE2RelOverlap is in bases
					int PE2SrcSeqLen,			// number of bases in pPE2SrcSeqWrd sequence
				   tSeqWrd4 *pPE2SrcSeqWrd,		// merge this sequence with pDstSeqWrd				   
				   int DstSeqLen,				// number of bases currently in pDstSeqWrd sequence
				   tSeqWrd4 *pDstSeqWrd,		// merge this sequence with pPE1SrcSeqWrd and pPE2SrcSeqWrd; merged sequence replaces this sequence
				   tSeqWrd4 *pTmpSeqWrd)		// temp sequence space for use as may be needed whilst merging
{
int PE1MergeLen;
int PE2MergeLen;

if(pTmpSeqWrd != NULL)
	*pTmpSeqWrd = 0;
PE1MergeLen = SeqMerge(PE1RelOverlap,PE1SrcSeqLen,pPE1SrcSeqWrd,DstSeqLen,pDstSeqWrd,pTmpSeqWrd,(char *)"SeqMergePE12ToSE A");
if(pTmpSeqWrd != NULL)
	*pTmpSeqWrd = 1;
if(PE1RelOverlap < 0)
	PE2RelOverlap += abs(PE1RelOverlap);
PE2MergeLen = SeqMerge(PE2RelOverlap,PE2SrcSeqLen,pPE2SrcSeqWrd,PE1MergeLen,pDstSeqWrd,pTmpSeqWrd,(char *)"SeqMergePE12ToSE B");
if(pTmpSeqWrd != NULL)
	*pTmpSeqWrd = 2;
return(PE2MergeLen);
}

// sequence merge; merges pSrcSeq and pDstSeq, merged sequence replaces pDstSeq 
int												// number of bases in sequence pDstSeqWrd after merging
CKangadna::SeqMerge(int RelOverlap,				// if >= 0 then pSrcSeq starts at this pDstSeq[RelOverlap]; if < 0 then pDstSeq starts at this pSrcSeq[RelOverlap]; RelOverlap is in bases
					int SrcSeqLen,				// number of bases in pSrcSeqWrd sequence
				   tSeqWrd4 *pSrcSeqWrd,		// merge this sequence with pDstSeqWrd
					int DstSeqLen,				// number of bases currently in pDstSeqWrd sequence
				   tSeqWrd4 *pDstSeqWrd,		// merge this sequence with pSrcSeqWrd; merged sequence replaces this sequence
				   tSeqWrd4 *pTmpSeqWrd,		// temp sequence space for use as may be needed whilst merging
					char *pszDiagText)			// optional diagnostics text

{
tSeqWrd4 *pSrcWrdBases;
tSeqWrd4 *pDstWrdBases;

tSeqWrd4 SrcWrd;
tSeqWrd4 DstWrd;
tSeqWrd4 SrcWrdMsk;
tSeqWrd4 DstWrdMsk;
tSeqWrd4 Base2Append;

int ExtdDstBy;

int DstWrdBaseIdx;
int SrcWrdBaseIdx;
int SrcWrdShr;
int DstWrdShl;

int SrcOfs;
int SrcLen;


int CurDstSeqLen;
int MergedDstSeqLen;

// if DstSeq already completely contains SrcSeq then no merge required
if(DstSeqLen > 0 && (RelOverlap >= 0 && (SrcSeqLen <= (DstSeqLen - RelOverlap))))
	return(DstSeqLen);

// if DstSeq completely contained within SrcSeq then just replace DstSeq with SrcSeq
if(DstSeqLen == 0 || (RelOverlap < 0 && (DstSeqLen <= (SrcSeqLen + RelOverlap))))
	{
	memmove(pDstSeqWrd,pSrcSeqWrd,sizeof(tSeqWrd4) * ((SrcSeqLen+14)/15));
	return(SrcSeqLen);
	}

// a little more processing required as neither SrcSeq or DstSeq completely contained within the other
if(RelOverlap < 0)	// if source is overlaying dest and dest will be extending then use TmpSeqWrd
	{
	memmove(pTmpSeqWrd,pDstSeqWrd,sizeof(tSeqWrd4) * ((DstSeqLen+14)/15)); // copy all words from DstSeqWrd into TmpSeqWrd
	memmove(pDstSeqWrd,pSrcSeqWrd,sizeof(tSeqWrd4) * ((SrcSeqLen+14)/15)); // and then replace DstSeqWrd with SrcSeqWrd; original DstSeqWrd bases will be used to extend
	CurDstSeqLen = SrcSeqLen;
	pSrcWrdBases = pTmpSeqWrd;	
	ExtdDstBy = abs(RelOverlap) + DstSeqLen - SrcSeqLen; // will be extending by this number of bases (1..N)
	if(ExtdDstBy < 1 || ExtdDstBy > DstSeqLen)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"SeqMerge (called from %s) Op1: ExtdDstBy %d, RelOverlap %d, DstSeqLen %d SrcSeqLen %d",pszDiagText == NULL ? "Not Specified" : pszDiagText,ExtdDstBy,RelOverlap,DstSeqLen,SrcSeqLen);
		}
	SrcOfs = DstSeqLen - ExtdDstBy;					// first base to copy comes from this base offset in SrcWrdBases
	SrcLen = DstSeqLen;
	}
else    // else dest is overlaying source and source will be extending the original DstSeqWrd
	{
	CurDstSeqLen = DstSeqLen;
	pSrcWrdBases = pSrcSeqWrd;
	ExtdDstBy = RelOverlap + SrcSeqLen - DstSeqLen;		// will be extending by this number of bases (1..N)
	if(ExtdDstBy < 1 || ExtdDstBy > SrcSeqLen)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"SeqMerge (called from %s) Op2: ExtdDstBy %d, RelOverlap %d, DstSeqLen %d SrcSeqLen %d",pszDiagText == NULL ? "Not Specified" : pszDiagText,ExtdDstBy,RelOverlap,DstSeqLen,SrcSeqLen);
		}

	SrcOfs = SrcSeqLen - ExtdDstBy;						// first base to copy comes from this base offset in SrcWrdBases
	SrcLen = SrcSeqLen;
	}

pSrcWrdBases +=  (SrcOfs / 15); // initialise to pt at src word containing 1st base to be appended onto the dst sequence
SrcWrdBaseIdx = SrcOfs % 15;  // initialise to base index (0..14) within src word
SrcWrdMsk = 0x30000000 >> (2 * SrcWrdBaseIdx); // initialise to base mask
SrcWrdShr = 2 * (14 - SrcWrdBaseIdx);

pDstWrdBases = &pDstSeqWrd[CurDstSeqLen / 15];	// initialise to pt at dst word into which 1st base is to be appended
DstWrdBaseIdx = CurDstSeqLen % 15;		        // initialise to base index (0..14) within dst word
DstWrdMsk = 0x03C000000  >> (2 * DstWrdBaseIdx);
DstWrdShl = 2 * (14 - DstWrdBaseIdx);

SrcWrd = *pSrcWrdBases++;
SrcWrd &= 0x3fffffff;
DstWrd = *pDstWrdBases;
if(DstWrdBaseIdx == 0)							// if dst starting on a word boundary then initialise as if completely empty partial
	DstWrd = 0x0fffffff;
else
	DstWrd &= 0x3fffffff;
MergedDstSeqLen = CurDstSeqLen;
while(ExtdDstBy--)
	{
	Base2Append = (SrcWrd >> SrcWrdShr) & 0x03;
	Base2Append <<= DstWrdShl;
	DstWrd &= ~DstWrdMsk;
	DstWrd |= Base2Append;
	MergedDstSeqLen += 1;
	if(!ExtdDstBy)
		{
		if(MergedDstSeqLen % 15)		// if last dst word containing < 15 bases then need to mark word as being a partial
			DstWrd |= cSeqWrd4PartSeq;
		*pDstWrdBases++ = DstWrd;
		*pDstWrdBases = cSeqWrd4EOS;
		break;
		}

	SrcWrdMsk >>= 2;
	SrcWrdShr -= 2;
	DstWrdMsk >>= 2;
	DstWrdShl -= 2;

	if(!SrcWrdMsk)
		{
		SrcWrdMsk = 0x30000000;
		SrcWrd = *pSrcWrdBases++;
		SrcWrd &= 0x3fffffff;
		SrcWrdShr = 28;
		}
	if(!DstWrdMsk)
		{
		*pDstWrdBases++ = DstWrd;
		DstWrdMsk = 0x3c000000;
		DstWrd = 0x0fffffff;
		DstWrdShl = 28;
		}
	}

return(MergedDstSeqLen);
}


INT64									// accumulated partial sequence lengths, or if <= 0 then error
CKangadna::SavePartialSeqs(int PE1Len,	// length of partially assembled sequence ptd at by pPE1
				tSeqWrd4 *pPE1,			// partially assembled PE1
				int PE2Len,				// (optional) length of partially assembled sequence ptd at by pPE2
				tSeqWrd4 *pPE2)			// (optional) partially assembled PE2
{
UINT64 PE1NumSeqWrds;
UINT64 PE2NumSeqWrds;
UINT64 LenPartialSeqs2Assemb;
int NumReqWrds;

tSeqWrd4 *pSeqWrds;
tSeqWrd4 *pPackSeq;

PE1NumSeqWrds = (PE1Len+14)/15;
if(PE2Len > 0)
	PE2NumSeqWrds = (PE2Len+14)/15;
else
	PE2NumSeqWrds = 0;

AcquireSerialiseSeqHdr();

// initial or realloc as may be required to hold these partial sequences
if(m_pPartialSeqs2Assemb == NULL || (sizeof(tSeqWrd4) * (m_PartialSeqs2AssembOfs + PE1NumSeqWrds + PE2NumSeqWrds + 10)) > m_AllocdPartialSeqs2Assemb)
	{
	if(m_AllocdPartialSeqs2Assemb == 0 || m_pPartialSeqs2Assemb == NULL)		// initial allocation required?
		{
		// assume will require 1/3 of unmerged reads - will realloc later if more memory is required
		m_AllocdPartialSeqs2Assemb = max(m_Sequences.AllocMemSeqs2Assemb/3,sizeof(tSeqWrd4) * 1000000);	

#ifdef _WIN32
		m_pPartialSeqs2Assemb = malloc((size_t)m_AllocdPartialSeqs2Assemb);	
		if(m_pPartialSeqs2Assemb == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"SavePartialSeqs: Concatenated packed sequences memory allocation of %llu bytes - %s",m_AllocdPartialSeqs2Assemb,strerror(errno));
			m_AllocdPartialSeqs2Assemb = 0;
			ReleaseSerialiseSeqHdr();
			Reset(false);
			return((INT64)eBSFerrMem);
			}
#else
		if((m_pPartialSeqs2Assemb = (void *)mmap(NULL,m_AllocdPartialSeqs2Assemb, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"SavePartialSeqs: Concatenated packed sequences memory of %llu bytes through mmap()  failed - %s",m_AllocdPartialSeqs2Assemb,strerror(errno));
			m_pPartialSeqs2Assemb = NULL;
			m_AllocdPartialSeqs2Assemb = 0;
			ReleaseSerialiseSeqHdr();
			Reset(false);
			return((INT64)eBSFerrMem);
			}
#endif

		m_NumPartialSeqs2Assemb = 0;
		m_LenPartialSeqs2Assemb = 0;
		m_PartialSeqs2AssembOfs = 1;
		*(tSeqWrd4 *)m_pPartialSeqs2Assemb = cSeqWrd4BOS; // obviously must be first sequence, so prefix with cSeqWrd4BOS
		}
	else      // else needing to extend current allocation
		{
		// realloc memory with a 25% increase over previous allocation 
		size_t memreq;
		void *pAllocd;
		memreq = (size_t)((m_AllocdPartialSeqs2Assemb * 125) / (UINT64)100);
#ifdef _WIN32
		pAllocd = realloc(m_pPartialSeqs2Assemb,memreq);
#else
		pAllocd = mremap(m_pPartialSeqs2Assemb,m_AllocdPartialSeqs2Assemb,memreq,MREMAP_MAYMOVE);
		if(pAllocd == MAP_FAILED)
			pAllocd = NULL;
#endif
		if(pAllocd == NULL)
			{
			ReleaseSerialiseSeqHdr();
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"SavePartialSeqs: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
			Reset(false);
			return(0);
			}

		m_pPartialSeqs2Assemb = pAllocd;
		m_AllocdPartialSeqs2Assemb = (UINT64)memreq;
		}
	}

// sufficent memory available to hold sequences
pSeqWrds = (tSeqWrd4 *)m_pPartialSeqs2Assemb;
pPackSeq = &pSeqWrds[m_PartialSeqs2AssembOfs];
NumReqWrds = (int)PE1NumSeqWrds + 1;
m_LenPartialSeqs2Assemb += PE1Len;
m_NumPartialSeqs2Assemb += 1;
if(PE2NumSeqWrds)
	{
	NumReqWrds += (int)PE2NumSeqWrds + 1;
	m_LenPartialSeqs2Assemb += PE2Len;
	m_NumPartialSeqs2Assemb += 1;
	}
m_PartialSeqs2AssembOfs += NumReqWrds;

pPackSeq[NumReqWrds] = cSeqWrd4EOS; // currently last sequence, so will finish with a cSeqWrd4EOS; next sequence overwrites this EOS marker

*pPackSeq++ = (tSeqWrd4)PE1Len << 2 | (tSeqWrd4)(PE2Len > 0 ? cFlgSeqPE : 0);
memmove(pPackSeq,pPE1,sizeof(tSeqWrd4) * (size_t)PE1NumSeqWrds);
pPackSeq += PE1NumSeqWrds;

if(PE2Len)
	{
	*pPackSeq++ = (tSeqWrd4)PE2Len << 2 | (tSeqWrd4)(cFlgSeqPE2 | cFlgSeqPE);
	memmove(pPackSeq,pPE2,sizeof(tSeqWrd4) * (size_t)PE2NumSeqWrds);
	}
LenPartialSeqs2Assemb = m_LenPartialSeqs2Assemb;
ReleaseSerialiseSeqHdr();
	
return((INT64)LenPartialSeqs2Assemb);
}




teBSFrsltCodes
CKangadna::SaveAsFasta(char *pszFastaFile)
{
int Written;
tSeqWrd4 *pCurSeq;
char *pszAsciiSeq;
tSeqID SeqID;
UINT32 SeqLen;
UINT32 SrcFileID;
char szFastaFilePE1[_MAX_PATH];
char szFastaFilePE2[_MAX_PATH];
static char szFastaR1[cMaxDiagAsciiSeqLen+1];
static char szFastaR2[cMaxDiagAsciiSeqLen+1];
int LineLenR1;
int FastaLenR1;
char *pszFastaR1;
int LineLenR2;
int FastaLenR2;
char *pszFastaR2;


strcpy(szFastaFilePE1,pszFastaFile);
strcat(szFastaFilePE1,".R1.fasta");
#ifdef _WIN32
if((m_hOutFastaR1 = open(szFastaFilePE1, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutFastaR1 = open(szFastaFilePE1, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",szFastaFilePE1,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output PE1 preprocessed multifasta file created/truncated: '%s'",szFastaFilePE1);
if(m_Sequences.bPESeqs)
	{
	strcpy(szFastaFilePE2,pszFastaFile);
	strcat(szFastaFilePE2,".R2.fasta");
#ifdef _WIN32
	if((m_hOutFastaR2 = open(szFastaFilePE2, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hOutFastaR2 = open(szFastaFilePE2, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",szFastaFilePE2,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output PE2 preprocessed multifasta file created/truncated: '%s'",szFastaFilePE2);
	}
pCurSeq = NULL;
pszFastaR1 = szFastaR1;
pszFastaR2 = szFastaR2;
szFastaR1[0] = '\0';
szFastaR2[0] = '\0';
FastaLenR1 = 0;
FastaLenR2 = 0;

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif

UINT32 NumFastaPE1Seqs;
UINT32 NumFastaPE2Seqs;
UINT32 FastaPE1MinLen;
UINT32 FastaPE1MaxLen;
UINT32 FastaPE2MinLen;
UINT32 FastaPE2MaxLen;
UINT64 TotFastaPE1SeqLen;
UINT64 TotFastaPE2SeqLen;
UINT32 SeqFlags;

NumFastaPE1Seqs = 0;
NumFastaPE2Seqs = 0;
FastaPE1MinLen = 0;
FastaPE1MaxLen = 0;
FastaPE2MinLen = 0;
FastaPE2MaxLen = 0;
TotFastaPE1SeqLen = 0;
TotFastaPE2SeqLen = 0;
while((pCurSeq = IterSeqHeaders(pCurSeq, // iterate to next sequence following (NULL to return 1st sequence)
							&SeqID,		// returned sequence identifier
							&SrcFileID,	// returned 8 bit source file identifier
							&SeqFlags,		// returned 16 bit sequence flags
							&SeqLen,false))!=NULL)		// returned 30 bit sequence length
	{
	if(!(SeqFlags & cFlgSeqPE2))	// if not a PE2 sequence
		{		
		if((FastaLenR1 + 200) >= cMaxDiagAsciiSeqLen)
			{
			Written = write(m_hOutFastaR1,szFastaR1,FastaLenR1);
			FastaLenR1 = 0;
			pszFastaR1 = szFastaR1;
			}
		NumFastaPE1Seqs += 1;
		if(FastaPE1MinLen == 0 || SeqLen < FastaPE1MinLen)
			FastaPE1MinLen = SeqLen;
		if(SeqLen > FastaPE1MaxLen)
			FastaPE1MaxLen = SeqLen;
		TotFastaPE1SeqLen += SeqLen;

		pszAsciiSeq = AsciifySequence(pCurSeq,SeqLen,false);
		if(!(SeqFlags & cFlgSeqPE))
			FastaLenR1 += sprintf((char *)&szFastaR1[FastaLenR1],">Seq%u|%d %d|%d|SE\n",NumFastaPE1Seqs,SrcFileID,SeqLen,SrcFileID);
		else
			FastaLenR1 += sprintf((char *)&szFastaR1[FastaLenR1],">Seq%u|%d/1 %d|%d|PE|1\n",NumFastaPE1Seqs,SrcFileID,SeqLen,SrcFileID);

		pszFastaR1 = &szFastaR1[FastaLenR1];
		LineLenR1 = 0;

		while(*pszAsciiSeq != '\0')
			{
			if(LineLenR1 > 75)
				{
				*pszFastaR1++ = '\n';
				FastaLenR1 += 1;
				LineLenR1 = 0;
				if((FastaLenR1 + 200) >= cMaxDiagAsciiSeqLen)
					{
					Written = write(m_hOutFastaR1,szFastaR1,FastaLenR1);
					FastaLenR1 = 0;
					pszFastaR1 = szFastaR1;
					}
				}
			*pszFastaR1++ = *pszAsciiSeq++;
			FastaLenR1 += 1;
			LineLenR1 += 1;
			}

		if(LineLenR1 > 0)
			{
			*pszFastaR1++ = '\n';
			FastaLenR1 += 1;
			LineLenR1 += 1;
			}
		}
	else // else PE2 so need to reverse complement the sequence
		{		
		if((FastaLenR2 + 200) >= cMaxDiagAsciiSeqLen)
			{
			Written = write(m_hOutFastaR2,szFastaR2,FastaLenR2);
			FastaLenR2 = 0;
			pszFastaR2 = szFastaR2;
			}
		
		NumFastaPE2Seqs += 1;
		if(FastaPE2MinLen == 0 || SeqLen < FastaPE2MinLen)
			FastaPE2MinLen = SeqLen;
		if(SeqLen > FastaPE2MaxLen)
			FastaPE2MaxLen = SeqLen;
		TotFastaPE2SeqLen += SeqLen;

		pszAsciiSeq = AsciifySequence(pCurSeq,SeqLen,true);
		FastaLenR2 += sprintf((char *)&szFastaR2[FastaLenR2],">Seq%u|%d/2 %d|%d|PE|2\n",NumFastaPE2Seqs,SrcFileID-1,SeqLen,SrcFileID);
		LineLenR2 = 0;
		pszFastaR2 = &szFastaR2[FastaLenR2];
		while(*pszAsciiSeq != '\0')
			{
			if(LineLenR2 > 75)
				{
				*pszFastaR2++ = '\n';
				FastaLenR2 += 1;
				LineLenR2 = 0;
				if((FastaLenR2 + 200) >= cMaxDiagAsciiSeqLen)
					{
					Written = write(m_hOutFastaR2,szFastaR2,FastaLenR2);
					FastaLenR2 = 0;
					pszFastaR2 = szFastaR2;
					}
				}
			*pszFastaR2++ = *pszAsciiSeq++;
			FastaLenR2 += 1;
			LineLenR2 += 1;
			}
		if(LineLenR2 > 0)
			{
			*pszFastaR2++ = '\n';
			FastaLenR2 += 1;
			LineLenR2 += 1;
			}
		}
	}

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif

if(m_hOutFastaR1 != -1)
	{
	if(FastaLenR1)
		Written = write(m_hOutFastaR1,szFastaR1,FastaLenR1);
	FastaLenR1 = 0;
#ifdef _WIN32
	_commit(m_hOutFastaR1);
#else
	fsync(m_hOutFastaR1);
#endif
	close(m_hOutFastaR1);
	m_hOutFastaR1 = -1;
	}
if(m_hOutFastaR2 != -1)
	{
	if(FastaLenR2)
		Written = write(m_hOutFastaR2,szFastaR2,FastaLenR2);
	FastaLenR2 = 0;
#ifdef _WIN32
	_commit(m_hOutFastaR2);
#else
	fsync(m_hOutFastaR2);
#endif
	close(m_hOutFastaR2);
	m_hOutFastaR2 = -1;
	}
#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
double MeanLen = TotFastaPE1SeqLen/(double)NumFastaPE1Seqs;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted %s Multifasta NumSeqs: %u, MaxLen: %u, MinLen: %u, TotLen: %llu, MeanLen: %1.1f",m_Sequences.bPESeqs ? "PE1" : "SE",
										NumFastaPE1Seqs,FastaPE1MaxLen,FastaPE1MinLen,TotFastaPE1SeqLen,MeanLen);
if(m_Sequences.bPESeqs)
	{
	MeanLen = TotFastaPE2SeqLen/(double)NumFastaPE2Seqs;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted PE2 Multifasta NumSeqs: %u, MaxLen: %u, MinLen: %u, TotLen: %llu, MeanLen: %1.1f",
										NumFastaPE2Seqs,FastaPE2MaxLen,FastaPE2MinLen,TotFastaPE2SeqLen,MeanLen);
	}

if(gProcessingID > 0)
	{
	char *pszEnd;
	if(!m_Sequences.bPESeqs)
		pszEnd = (char *)"SEReadsAccepted";
	else
		pszEnd = (char *)"PE1ReadsAccepted";
	MeanLen = TotFastaPE1SeqLen/(double)NumFastaPE1Seqs;
	gSQLiteSummaries.AddResult(gProcessingID,pszEnd,ePTUint32,sizeof(NumFastaPE1Seqs),"Cnt",&NumFastaPE1Seqs);
	gSQLiteSummaries.AddResult(gProcessingID,pszEnd,ePTUint32,sizeof(FastaPE1MinLen),"MinLen",&FastaPE1MinLen);
	gSQLiteSummaries.AddResult(gProcessingID,pszEnd,ePTUint32,sizeof(FastaPE1MaxLen),"MaxLen",&FastaPE1MaxLen);
	gSQLiteSummaries.AddResult(gProcessingID,pszEnd,ePTUint64,sizeof(TotFastaPE1SeqLen),"TotLen",&TotFastaPE1SeqLen);
	gSQLiteSummaries.AddResult(gProcessingID,pszEnd,ePTDouble,sizeof(MeanLen),"MeanLen",&MeanLen);

	if(m_Sequences.bPESeqs)
		{
		MeanLen = TotFastaPE2SeqLen/(double)NumFastaPE2Seqs;
		pszEnd = (char *)"PE2ReadsAccepted";
		gSQLiteSummaries.AddResult(gProcessingID,pszEnd,ePTUint32,sizeof(NumFastaPE2Seqs),"Cnt",&NumFastaPE2Seqs);
		gSQLiteSummaries.AddResult(gProcessingID,pszEnd,ePTUint32,sizeof(FastaPE2MinLen),"MinLen",&FastaPE2MinLen);
		gSQLiteSummaries.AddResult(gProcessingID,pszEnd,ePTUint32,sizeof(FastaPE2MaxLen),"MaxLen",&FastaPE2MaxLen);
		gSQLiteSummaries.AddResult(gProcessingID,pszEnd,ePTUint64,sizeof(TotFastaPE2SeqLen),"SeqLen",&TotFastaPE2SeqLen);
		gSQLiteSummaries.AddResult(gProcessingID,pszEnd,ePTDouble,sizeof(MeanLen),"MeanLen",&MeanLen);
		}
	}
return(eBSFSuccess);
}

void
CKangadna::DumpSeqs2Assemble(char *pszHeading,		// use this header text
								UINT32 Max2Dump)		// dump at most this many sequences, 0 if no limit
{
tSeqWrd4 *pCurSeq;
tSeqID SeqID;
UINT32 SrcFileID;
UINT32 Flgs;
UINT32 SeqLen;
pCurSeq = NULL;
if(Max2Dump == 0)
	Max2Dump = 0xffffffff;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sequence Dump: %s",pszHeading);

while(Max2Dump-- && (pCurSeq = IterSeqHeaders(pCurSeq, // iterate to next sequence following (NULL to return 1st sequence)
							&SeqID,		// returned sequence identifier
							&SrcFileID,		// returned 8 bit source file identifier
							&Flgs,			// returned 16 bit sequence flags
							&SeqLen,false))!=NULL)		// returned 30 bit sequence length
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"SeqID: %u FileID: %u Flags: 0x%x SeqLen: %u\n%s\n",SeqID,SrcFileID,Flgs,SeqLen,AsciifySequence(pCurSeq,200));

}


// AsciifySequence
// Use for diagnostics: will generate ascii representation of a sequence
// Sequence will be truncated if longer than MaxSeqLen (clamped to be no longer than cMaxDiagAsciiSeqLen)


char *
CKangadna::AsciifySequence(void *pSeqWrds,		// sequence words to asciify
				UINT32 MaxSeqLen,				// asciify at most MaxSeqLen bases
				bool bRevCpl)					// true if sequence to be reverse complemented before asciifying
{
static UINT8 szDump[cMaxDiagAsciiSeqLen+1];
UINT8 *pBase;
int Idx;
int SeqLen;

if(MaxSeqLen == 0)
	MaxSeqLen = cMaxDiagAsciiSeqLen;

SeqLen = GetSeq(pSeqWrds,			// sequence of interest
					(UINT8 *)szDump,	// where to copy unpacked sequence bases
					min(MaxSeqLen,sizeof(szDump)-1));	// return at most MaxSeqLen bases

if(bRevCpl)
	RevCplSeq(SeqLen,szDump);

pBase = szDump;
for(Idx = 0; Idx < SeqLen; Idx++,pBase++)	// alow for terminating '\0'
	{
	switch(*pBase & 0x07) {
		case eBaseA:
			*pBase = 'a';
			break;
		case eBaseC:
			*pBase = 'c';
			break;
		case eBaseG:
			*pBase = 'g';
			break;
		case eBaseT:
			*pBase = 't';
			break;
		default:			// treat any other as indeterminate or 'N'
			*pBase = 'n';
			break;
		}
	}
*pBase = '\0';
return((char *)szDump);
}


char *
CKangadna::AsciifySequence(tSeqID SeqID,		// sequence identifier
				UINT32 MaxSeqLen,				// asciify at most MaxSeqLen bases
				bool bRevCpl)					// true if sequence to be reverse complemented before asciifying
{
static UINT8 szDump[cMaxDiagAsciiSeqLen+1];
UINT8 *pBase;
int Idx;
int SeqLen;

if(MaxSeqLen == 0)
	MaxSeqLen = cMaxDiagAsciiSeqLen;

SeqLen = GetSeq(SeqID,				// sequence identifier
					(UINT8 *)szDump,	// where to copy unpacked sequence bases
					min(MaxSeqLen,sizeof(szDump)-1));	// return at most MaxSeqLen bases

if(m_pBlockNsLoci != NULL && m_NumBlockNsLoci > 0)
	{
	// replace any bases which were indeterminates when originally loaded with N's
	int MinBlockNsLociIdx;
	int MaxBlockNsLociIdx;
	int ProbeBlockNsLociIdx;
	tsBlockNsLoci * pBlockNsLoci;
	
	MinBlockNsLociIdx = 0;
	MaxBlockNsLociIdx = m_NumBlockNsLoci-1;

	do {
		ProbeBlockNsLociIdx = (MaxBlockNsLociIdx + MinBlockNsLociIdx) / 2;
		pBlockNsLoci = &m_pBlockNsLoci[ProbeBlockNsLociIdx];
		if(pBlockNsLoci->SeqID == SeqID)
			{
			while(pBlockNsLoci->SeqID == SeqID)
				{
				if(ProbeBlockNsLociIdx > 0 && pBlockNsLoci[-1].SeqID == SeqID)
					{
					pBlockNsLoci -= 1;
					ProbeBlockNsLociIdx -= 1;
					}
				else
					break;
				}
			while(ProbeBlockNsLociIdx < (int)m_NumBlockNsLoci && pBlockNsLoci->SeqID == SeqID)
				{
				memset(&szDump[pBlockNsLoci->Ofs],eBaseN,pBlockNsLoci->NumNs);
				pBlockNsLoci += 1;
				ProbeBlockNsLociIdx += 1;
				}
			break;
			}
		 if(pBlockNsLoci->SeqID > SeqID)
			MaxBlockNsLociIdx = ProbeBlockNsLociIdx - 1;
		else
			MinBlockNsLociIdx = ProbeBlockNsLociIdx + 1;
		}
	while(MinBlockNsLociIdx <= MaxBlockNsLociIdx);
	}


if(bRevCpl)
		RevCplSeq(SeqLen,szDump);

pBase = szDump;
for(Idx = 0; Idx < SeqLen; Idx++,pBase++)	// alow for terminating '\0'
	{
	switch(*pBase) {
		case eBaseA:
			*pBase = 'a';
			break;
		case eBaseC:
			*pBase = 'c';
			break;
		case eBaseG:
			*pBase = 'g';
			break;
		case eBaseT:
			*pBase = 't';
			break;
		default:			// treat any other as indeterminate or 'N'
			*pBase = 'n';
			break;
		}
	}
*pBase = '\0';
return((char *)szDump);
}


// AsciifySequences
// Use for diagnostics: will generate ascii representation of two sequences with non-matching bases between the sequences
// in uppercase and matching bases in lowercase
// Sequences will be truncated if longer than MaxSeqLen (clamped to be no longer than cMaxDiagAsciiSeqLen)
char *
CKangadna::AsciifySequences(tSeqID SeqID1,		// sequence identifier
				tSeqID SeqID2,		// sequence identifier
				UINT32 MaxSeqLen,	// asciify at most MaxSeqLen bases
				char **ppszSeq1,		// returned seq1
				char **ppszSeq2)		// returned seq2
{
static UINT8 szDump1[cMaxDiagAsciiSeqLen+1];		// alow for terminating '\0'
static UINT8 szDump2[cMaxDiagAsciiSeqLen+1];
UINT8 *pBase1;
UINT8 *pBase2;
int Idx;
bool bDiff;
int SeqLen1;
int SeqLen2;

SeqLen1 = GetSeq(SeqID1,		// sequence identifier
					(UINT8 *)szDump1,	// where to copy unpacked sequence bases
					min((size_t)MaxSeqLen,sizeof(szDump1)-1));	// return at most MaxSeqLen bases

SeqLen2 = GetSeq(SeqID2,		// sequence identifier
					(UINT8 *)szDump2,	// where to copy unpacked sequence bases
					min((size_t)MaxSeqLen,sizeof(szDump2)-1));	// return at most MaxSeqLen bases

MaxSeqLen = min(SeqLen1,SeqLen2);

pBase1 = szDump1;
pBase2 = szDump2;
for(Idx = 0; Idx < (int)MaxSeqLen; Idx++,pBase1++,pBase2++)
	{
	if(*pBase1 != *pBase2)
		bDiff = true;
	else
		bDiff = false;
	switch(*pBase1) {
		case eBaseA:
			*pBase1 = bDiff ? 'A' : 'a';
			break;
		case eBaseC:
			*pBase1 = bDiff ? 'C' :  'c';
			break;
		case eBaseG:
			*pBase1 = bDiff ? 'G' : 'g';
			break;
		case eBaseT:
			*pBase1 = bDiff ? 'T' : 't';
			break;
		}
	switch(*pBase2) {
		case eBaseA:
			*pBase2 = bDiff ? 'A' : 'a';
			break;
		case eBaseC:
			*pBase2 = bDiff ? 'C' :  'c';
			break;
		case eBaseG:
			*pBase2 = bDiff ? 'G' : 'g';
			break;
		case eBaseT:
			*pBase2 = bDiff ? 'T' : 't';
			break;
		}
	}
*pBase1 = '\0';
*pBase2 = '\0';

if(ppszSeq1 != NULL)
	*ppszSeq1 = (char *)szDump1;
if(ppszSeq2 != NULL)
	*ppszSeq2 = (char *)szDump2;

return((char *)szDump1);
}

// CopyPackedSeq
int							// returned number of packed sequence words copied
CKangadna::CopyPackedSeq(void *pSrcSeq,	// copy from this start tSeqWrd
			  void *pDstSeq,	// copy to this start tSeqWrd
			  int MaxLen,		// copy at most this many tSeqWrds
			  bool bTermEOS)	// if true then terminate copied sequence with EOS tSeqWrd
{
int NumSeqWrds;

NumSeqWrds = 0;

tSeqWrd4 *pSrcWrd = (tSeqWrd4 *)pSrcSeq;
tSeqWrd4 *pDstWrd = (tSeqWrd4 *)pDstSeq;
tSeqWrd4 SeqWrd;
while(MaxLen-- && (((SeqWrd = *pSrcWrd++) & cSeqWrd4MSWHdr) != cSeqWrd4MSWHdr))
	{
	NumSeqWrds += 1;
	*pDstWrd++ = SeqWrd;
	}
if(bTermEOS)
	*pDstWrd = cSeqWrd4EOS;

return(NumSeqWrds);
}


// check and return the maximal overlap of the 3' probe flank onto the target sequence
int									// length overlap of probe onto target, -1 if no overlap of at least MinOverlap
IsOverlapping(int MinOverlap,		// any overlap must be at least this length
				int ProbeLen,		// probe (putative overlapping sequence) length
				etSeqBase *pProbe,	// pts to first 5' probe base
				int TargLen,		// target (putative overlapped sequence) length
			  etSeqBase *pTarg,		// pts to first 5' target base
			  int *pProbeOfs)		// overlap starts at this offset from 5' probe
{
etSeqBase *pPBase;
etSeqBase *pTBase;
etSeqBase *pBase1Mrk;
int ProbeOfs;
int Idx1;
int Idx2;
int Idx1Max;
int Idx2Max;

*pProbeOfs = 0;
pPBase = pProbe;
pTBase = pTarg;
Idx1Max =   1 + ProbeLen - MinOverlap;
for(Idx1 = 0; Idx1 < Idx1Max; Idx1++,pPBase++)
	{
	if(*pPBase != *pTBase)			// must match to be the start of a putative overlap
		continue;

	// putative start of an overlap
	pBase1Mrk = pPBase;				// mark start of putative overlap
	ProbeOfs = Idx1;				// starting at this probe offset

	Idx2Max = min(TargLen,ProbeLen - ProbeOfs); // overlap can be at most this length
	pPBase+=1;
	pTBase+=1;
	for(Idx2 = 1; Idx2 < Idx2Max; Idx2++,pPBase++,pTBase++)
		if(*pPBase != *pTBase)	// if not matching then must have been a false overlap
			break;
	if(Idx2 != Idx2Max)			// if not matched over full Idx2Max then must have been a false overlap
		{
		pPBase = pBase1Mrk;
		continue;
		}

	// have an overlap, must be the maximal so return
	*pProbeOfs = ProbeOfs;
	return(Idx2Max);
	}
return(-1);							// no overlap of at least MinOverlap
}

// inplace reverse complement of unpacked bases
void
CKangadna::RevCplSeq(unsigned int SeqLen, // sequence to reverse complement is of this length
							etSeqBase *pSeq)	// pts to start of sequence to be reverse complemented
{
etSeqBase *pExch;
etSeqBase Tmp1;
etSeqBase Tmp2;

if(SeqLen < 1 || pSeq == NULL)
	return;
if(SeqLen == 1)
	{
	Tmp1 = *pSeq & 0x07;
	if(Tmp1 != eBaseN)
		*pSeq = (~Tmp1) & 0x03;
	return;
	}
pExch = pSeq + SeqLen-1;
SeqLen >>= 1;

while(SeqLen--)
	{
	Tmp1 = *pSeq & 0x07;
	Tmp2 = *pExch & 0x07;
	if(Tmp1 != eBaseN)
		Tmp1 = ~Tmp1 & 0x03;
	if(Tmp2 != eBaseN)
		Tmp2 = ~Tmp2 & 0x03;
	*pExch-- = Tmp1;
	*pSeq++ = Tmp2;
	}
if(pSeq == pExch && (Tmp1 = (*pSeq & 0x07)) != eBaseN)
	*pSeq = ~Tmp1 & 0x03;
}


