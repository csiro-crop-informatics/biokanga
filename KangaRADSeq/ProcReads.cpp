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

#include "StackSeqs.h"

static UINT8 *m_xpConcatSeqs;

// SetNumThreads
// Set allowed number of threads
void
CStackSeqs::SetNumThreads(int maxThreads)
{
m_NumThreads = maxThreads;
m_MTqsort.SetMaxThreads(maxThreads);
}


// CreateMutexes
// Create and initialise as appropriate all serialisation mutexes and locks
int
CStackSeqs::CreateMutexes(void)
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
#endif
	return(eBSFerrInternal);
	}

m_bMutexesCreated = true;
return(eBSFSuccess);
}

// DeleteMutexes
// Finished with serialisation mutexes and locks
void
CStackSeqs::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
#ifdef _WIN32
CloseHandle(m_hMtxIterReads);
CloseHandle(m_hMtxMHReads);

#else
pthread_mutex_destroy(&m_hMtxIterReads);
pthread_mutex_destroy(&m_hMtxMHReads);
pthread_rwlock_destroy(&m_hRwLock);

#endif
m_bMutexesCreated = false;
}

void
CStackSeqs::AcquireSerialise(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
}

void
CStackSeqs::ReleaseSerialise(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
}

void
inline CStackSeqs::AcquireLock(bool bExclusive)
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
inline CStackSeqs::ReleaseLock(bool bExclusive)
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
CStackSeqs::SetPairedEndProc(bool bIsPairedEndProc)
{
m_bIsPairedEndProc = bIsPairedEndProc;
}

void
CStackSeqs::SetCtgDescr(char *pszCtgDescr)
{
strncpy(m_szCtgDescr,pszCtgDescr,sizeof(m_szCtgDescr)-1);
m_szCtgDescr[sizeof(m_szCtgDescr)-1] = '\0';
}

bool
CStackSeqs::SetMaxMemWorkSetSize(size_t Bytes)
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



teBSFrsltCodes
CStackSeqs::AllocReadsMemory(size_t P1ReqAllocSize,size_t P2ReqAllocSize)
{
// ensure it's worth the trouble and effort!
if(P1ReqAllocSize < cMinConcatSeqLen)
	P1ReqAllocSize = cMinConcatSeqLen;
if(m_bIsPairedEndProc && (P2ReqAllocSize > 0 && P2ReqAllocSize < cMinConcatSeqLen))
	P2ReqAllocSize = cMinConcatSeqLen;

AcquireSerialise();
if(m_pP1Seqs2Assemb == NULL)		// will be NULL first time
	{
	m_AllocMemP1Seqs2Assemb = (size_t)P1ReqAllocSize;
#ifdef _WIN32
	m_pP1Seqs2Assemb = (UINT8 *) malloc((size_t)m_AllocMemP1Seqs2Assemb);
	if(m_pP1Seqs2Assemb == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Concatenated sequences P1 memory allocation of %lld bytes - %s",(INT64)m_AllocMemP1Seqs2Assemb,strerror(errno));
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pP1Seqs2Assemb = (UINT8 *)mmap(NULL,m_AllocMemP1Seqs2Assemb, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pP1Seqs2Assemb == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Concatenated sequences P1 memory of %lld bytes through mmap()  failed - %s",(INT64)m_AllocMemP1Seqs2Assemb,strerror(errno));
		m_pP1Seqs2Assemb = NULL;
		return(eBSFerrMem);
		}
#endif
	m_P1Seqs2AssembLen = 0;
	m_NumP1Seqs2Assemb = 0;
	}
else
	{
	UINT8 *pDstSeq;
	size_t memreq;
	if((m_P1Seqs2AssembLen + P1ReqAllocSize + 100) >= m_AllocMemP1Seqs2Assemb)		// 100 as a small safety margin!
		{
		memreq = (size_t)(m_P1Seqs2AssembLen + cReallocConcatSeqs + P1ReqAllocSize + 100);
#ifdef _WIN32
		pDstSeq = (UINT8 *) realloc(m_pP1Seqs2Assemb,memreq);
#else
		pDstSeq = (UINT8 *)mremap(m_pP1Seqs2Assemb,m_AllocMemP1Seqs2Assemb,memreq,MREMAP_MAYMOVE);
		if(pDstSeq == MAP_FAILED)
			pDstSeq = NULL;
#endif
		if(pDstSeq == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: P1 Memory re-allocation to %lld bytes - %s",memreq,strerror(errno));
			return(eBSFerrMem);
			}
		m_AllocMemP1Seqs2Assemb = memreq;
		m_pP1Seqs2Assemb = pDstSeq;
		}
	}

if(P2ReqAllocSize > 0)
	{
	if(m_pP2Seqs2Assemb == NULL)		// will be NULL first time
		{
		m_AllocMemP2Seqs2Assemb = (size_t)P2ReqAllocSize;
#ifdef _WIN32
		m_pP2Seqs2Assemb = (UINT8 *) malloc((size_t)m_AllocMemP2Seqs2Assemb);
		if(m_pP2Seqs2Assemb == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Concatenated P2 sequences memory allocation of %lld bytes - %s",(INT64)m_AllocMemP2Seqs2Assemb,strerror(errno));
			return(eBSFerrMem);
			}
#else
		// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
		m_pP2Seqs2Assemb = (UINT8 *)mmap(NULL,m_AllocMemP2Seqs2Assemb, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
		if(m_pP2Seqs2Assemb == MAP_FAILED)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Concatenated P2 sequences memory of %lld bytes through mmap()  failed - %s",(INT64)m_AllocMemP2Seqs2Assemb,strerror(errno));
			m_pP2Seqs2Assemb = NULL;
			return(eBSFerrMem);
			}
#endif
		m_P2Seqs2AssembLen = 0;
		m_NumP2Seqs2Assemb = 0;
		}
	else
		{
		UINT8 *pDstSeq;
		size_t memreq;
		if((m_P2Seqs2AssembLen + P2ReqAllocSize + 100) >= m_AllocMemP2Seqs2Assemb)		// 100 as a small safety margin!
			{
			memreq = (size_t)(m_P2Seqs2AssembLen + cReallocConcatSeqs + P2ReqAllocSize + 100);

#ifdef _WIN32
			pDstSeq = (UINT8 *) realloc(m_pP2Seqs2Assemb,memreq);
#else
			pDstSeq = (UINT8 *)mremap(m_pP2Seqs2Assemb,m_AllocMemP2Seqs2Assemb,memreq,MREMAP_MAYMOVE);
			if(pDstSeq == MAP_FAILED)
				pDstSeq = NULL;
#endif
			if(pDstSeq == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Memory P2 re-allocation to %lld bytes - %s",memreq,strerror(errno));
				return(eBSFerrMem);
				}
			m_AllocMemP2Seqs2Assemb = memreq;
			m_pP2Seqs2Assemb = pDstSeq;
			}
		}
	}

if(m_pSeqStarts == NULL)
	{
	m_AllocSeqStarts = cReallocNumReads;
	m_AllocMemSeqStarts = (size_t)cReallocNumReads * sizeof(tsSeqStarts);

#ifdef _WIN32
	m_pSeqStarts = (tsSeqStarts *) malloc((size_t)m_AllocMemSeqStarts);
	if(m_pSeqStarts == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: map sequence starts memory allocation of %lld bytes - %s",(INT64)m_AllocMemSeqStarts,strerror(errno));
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pSeqStarts = (tsSeqStarts *)mmap(NULL,m_AllocMemSeqStarts, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSeqStarts == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: map sequence starts memory allocation of %lld bytes through mmap()  failed - %s",(INT64)m_AllocMemSeqStarts,strerror(errno));
		m_pSeqStarts = NULL;
		return(eBSFerrMem);
		}
	#endif
	m_NumSeqStarts = 0;
	}
else
	{
	tsSeqStarts *pDstSeq;
	size_t memreq;
	if((m_NumSeqStarts + 100) >= m_AllocSeqStarts)		// 100 as a small safety margin!
		{
		m_AllocSeqStarts += cReallocNumReads;
		memreq = ((size_t)m_AllocSeqStarts * sizeof(tsSeqStarts));
#ifdef _WIN32
		pDstSeq = (tsSeqStarts *) realloc(m_pSeqStarts,memreq);
#else
		pDstSeq = (UINT8 *)mremap(m_pSeqStarts,m_AllocMemSeqStarts,memreq,MREMAP_MAYMOVE);
		if(pDstSeq == MAP_FAILED)
			pDstSeq = NULL;
#endif
		if(pDstSeq == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Memory P2 re-allocation to %lld bytes - %s",memreq,strerror(errno));
			return(eBSFerrMem);
			}
		m_AllocMemSeqStarts = memreq;
		m_pSeqStarts = pDstSeq;
		}
	}

ReleaseSerialise();
return(eBSFSuccess);
}

// LoadRawReads
// Load reads from fasta or fastq formated raw reads file
// basic assumption is that the paired ends are in the same order in their respective paired end files
// so when reads are loaded from file P1 and discarded for whatever reason, then the corresponding read can be discarded from P2
teBSFrsltCodes
CStackSeqs::LoadRawReads(int MaxNs,			// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
					int Trim5,				// trim this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int MinSeqLen,		    // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)
					int FileID,				// uniquely identifies source file P1 and P2
					char *pszP1File,		// process from this P1 file
					char *pszP2File)		// optionally process from this P2 file
{
static int FileNamesOfs = 0;

bool bIsPairedEndProc;						// true when paired end processing of read files required

teBSFrsltCodes Rslt;

CFasta P1Fasta;
int P1DescrLen;
UINT8 szP1DescrBuff[1024];
etSeqBase *pP1Seq;
UINT8 *pP1RawReadsBuff;
int P1ReadLen;
bool bP1IsFastq;
size_t P1ReqAllocSize;

CFasta P2Fasta;
int P2DescrLen;
UINT8 szP2DescrBuff[1024];
etSeqBase *pP2Seq;
UINT8 *pP2RawReadsBuff;
int P2ReadLen;
bool bP2IsFastq;
size_t P2ReqAllocSize;

int NumDescrReads;
int NumUnderlen;
int NumAcceptedReads;


int Num2Process;
int NumReads;



int NumInvalValues = 0;
int NumUnsupportedBases = 0;
int NumUnderlength = 0;
int NumExcessNs = 0;

bIsPairedEndProc = false;
if(pszP2File != NULL && pszP2File[0] != '\0')
	bIsPairedEndProc = true;
else
	bIsPairedEndProc = m_bIsPairedEndProc;

if(m_bIsPairedEndProc != bIsPairedEndProc)
	{
	if(m_bIsPairedEndProc == true)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Expected to be processing paired end reads (P1 + P2) not single ended P1s...");
	else
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Expected to be processing single ended reads (P1) not paired end (P1 + P2)...");
	Reset(false);
	return(eBSFerrInternal);
	}

if((pP1RawReadsBuff = new UINT8 [cAllocRawSeqLen]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for raw reads buffering...");
	Reset(false);
	return(eBSFerrMem);
	}

if((Rslt=(teBSFrsltCodes)P1Fasta.Open(pszP1File,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",pszP1File,P1Fasta.ErrText((teBSFrsltCodes)Rslt),P1Fasta.GetErrMsg());
	delete pP1RawReadsBuff;
	Reset(false);
	return(Rslt);
	}

if(P1Fasta.IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Only basespace stacks supported  but sequences from '%s' are in SOLiD colorspace...",pszP1File);
	delete pP1RawReadsBuff;
	P1Fasta.Close();
	Reset(false);
	return(eBSFerrOpnFile);
	}
bP1IsFastq = P1Fasta.IsFastq();

if(bIsPairedEndProc)
	{
	if((pP2RawReadsBuff = new UINT8 [cAllocRawSeqLen]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for raw reads buffering...");
		delete pP1RawReadsBuff;
		P1Fasta.Close();
		Reset(false);
		return(eBSFerrMem);
		}

	if((Rslt=(teBSFrsltCodes)P2Fasta.Open(pszP2File,true))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",pszP2File,P2Fasta.ErrText((teBSFrsltCodes)Rslt),P2Fasta.GetErrMsg());
		delete pP2RawReadsBuff;
		delete pP1RawReadsBuff;
		P1Fasta.Close();
		Reset(false);
		return(Rslt);
		}

	if(P2Fasta.IsSOLiD())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Only basespace stacks supported  but sequences from '%s' are in SOLiD colorspace...",pszP2File);
		delete pP2RawReadsBuff;
		delete pP1RawReadsBuff;
		P1Fasta.Close();
		Reset(false);
		return(eBSFerrOpnFile);
		}

	bP2IsFastq = P2Fasta.IsFastq();
	}
else
	{
	pP2RawReadsBuff = NULL;
	bP2IsFastq = false;
	}

// get file size and alloc for expected total sequences up front to try and reduce the number of reallocs which are required
// may end up allocating more memory than actually needed but on Windows HPC seems to result in a significant throughput improvement
P1ReqAllocSize = (size_t)P1Fasta.InitialFileSize();
if(bP1IsFastq)
	P1ReqAllocSize /= 2;			// quality scores effectively double the file size relative to the actual sequences
if(P1ReqAllocSize < cMinConcatSeqLen)
	P1ReqAllocSize = cMinConcatSeqLen;

if(bIsPairedEndProc)
	{
	P2ReqAllocSize = (size_t)P2Fasta.InitialFileSize();
	if(bP2IsFastq)
		P2ReqAllocSize /= 2;	// quality scores effectively double the file size relative to the actual sequences
	if(P2ReqAllocSize < cMinConcatSeqLen)
		P2ReqAllocSize = cMinConcatSeqLen;
	}
else
	P2ReqAllocSize = 0;

if((Rslt = AllocReadsMemory(P1ReqAllocSize,P2ReqAllocSize))!= eBSFSuccess)
	{
	P1Fasta.Close();
	delete pP1RawReadsBuff;
	if(bIsPairedEndProc)
		{
		delete pP2RawReadsBuff;
		P2Fasta.Close();
		}
	Reset(false);
	return(Rslt);
	}

size_t MemWorkSetBytes = ((m_AllocMemP1Seqs2Assemb + m_AllocMemP1Seqs2Assemb) * 3) / 2;

if(MemWorkSetBytes > m_CurMaxMemWorkSetBytes)
	{
	if(!SetMaxMemWorkSetSize(MemWorkSetBytes  * 2))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Length of loaded concatenated sequences exceeds limit of %lld bytes",cMaxConcatSeqLen);
		P1Fasta.Close();
		delete pP1RawReadsBuff;
		if(bIsPairedEndProc)
			{
			delete pP2RawReadsBuff;
			P2Fasta.Close();
			}
		Reset(false);
		return(eBSFerrMaxDirEls);
		}
	}

m_P1RdsSfxHdr.RdsSrcFiles[m_P1RdsSfxHdr.NumSrcFiles].SrcFileID = FileID;
strncpy((char *)m_P1RdsSfxHdr.RdsSrcFiles[m_P1RdsSfxHdr.NumSrcFiles].SrcFileName,pszP1File,sizeof(m_P1RdsSfxHdr.RdsSrcFiles[m_P1RdsSfxHdr.NumSrcFiles].SrcFileName)-1);
m_P1RdsSfxHdr.NumSrcFiles += 1;

if(bIsPairedEndProc)
	{
	m_P2RdsSfxHdr.RdsSrcFiles[m_P2RdsSfxHdr.NumSrcFiles].SrcFileID = FileID;
	strncpy((char *)m_P2RdsSfxHdr.RdsSrcFiles[m_P2RdsSfxHdr.NumSrcFiles].SrcFileName,pszP2File,sizeof(m_P2RdsSfxHdr.RdsSrcFiles[m_P2RdsSfxHdr.NumSrcFiles].SrcFileName)-1);
	m_P2RdsSfxHdr.NumSrcFiles += 1;
	}
else
	m_P2RdsSfxHdr.NumSrcFiles = 0;

NumUnsupportedBases = 0;
NumDescrReads = 0;
NumUnderlen = 0;
Num2Process = 0;
NumReads = 0;
while((Rslt = (teBSFrsltCodes)(P1ReadLen = P1Fasta.ReadSequence(pP1RawReadsBuff,cAllocRawSeqLen,true,false))) > eBSFSuccess)
	{
	if(P1ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		NumDescrReads += 1;
		P1DescrLen = P1Fasta.ReadDescriptor((char *)szP1DescrBuff,sizeof(szP1DescrBuff)-1);
		szP1DescrBuff[sizeof(szP1DescrBuff)-1] = '\0';
		P1ReadLen = P1Fasta.ReadSequence(pP1RawReadsBuff,cAllocRawSeqLen);


		// if paired end processing then try P2 - should also be a descriptor line
		if(bIsPairedEndProc)
			{
			Rslt = (teBSFrsltCodes)(P2ReadLen = P2Fasta.ReadSequence(pP2RawReadsBuff,cAllocRawSeqLen,true,false));
			if(Rslt <= 	eBSFSuccess)
				{
				while(P1Fasta.NumErrMsgs())
					gDiagnostics.DiagOut(eDLFatal,gszProcName,P2Fasta.GetErrMsg());
				delete pP1RawReadsBuff;
				delete pP2RawReadsBuff;
				P1Fasta.Close();
				P2Fasta.Close();
				}
			if(P2ReadLen != eBSFFastaDescr)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"raw sequence file '%s' processing error: %s ",pszP2File,
														Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
				delete pP1RawReadsBuff;
				delete pP2RawReadsBuff;
				P1Fasta.Close();
				P2Fasta.Close();
				return(eBSFerrParse);
				}

			P2DescrLen = P2Fasta.ReadDescriptor((char *)szP2DescrBuff,sizeof(szP2DescrBuff)-1);
			szP2DescrBuff[sizeof(szP2DescrBuff)-1] = '\0';
			}

		if(P1ReadLen < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence after %d reads parsed from '%s'",NumDescrReads-1,pszP1File);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szP1DescrBuff);
			delete pP1RawReadsBuff;
			P1Fasta.Close();
			if(bIsPairedEndProc)
				{
				delete pP2RawReadsBuff;
				P2Fasta.Close();
				}
			return(eBSFerrParse);
			}

		if(bIsPairedEndProc)
			{
			P2ReadLen = P2Fasta.ReadSequence(pP2RawReadsBuff,cAllocRawSeqLen);
			if(P2ReadLen < 0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence after %d reads parsed from '%s'",NumDescrReads-1,pszP2File);
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szP2DescrBuff);
				delete pP1RawReadsBuff;
				delete pP2RawReadsBuff;
				P1Fasta.Close();
				P2Fasta.Close();
				return(eBSFerrParse);
				}
			}

		if(P1ReadLen < (Trim5 + Trim3 + MinSeqLen))
			{
			NumUnderlen += 1;
			continue;
			}

		if(bIsPairedEndProc && (P2ReadLen < (Trim5 + Trim3 + MinSeqLen)))
			{
			NumUnderlen += 1;
			continue;
			}


		pP1Seq = pP1RawReadsBuff;

		// trim 5' and 3' as requested
		if(Trim5 > 0)
			{
			pP1Seq += Trim5;
			P1ReadLen -= Trim5;
			}

		P1ReadLen -= Trim3;

		// check for excessive number of Ns
		int Idx;
		int NumNs = 0;		// number of indeterminate bases in last 100bp window
		etSeqBase *pBase = pP1Seq;
		for(Idx = 0; Idx < P1ReadLen; Idx++,pBase++)
			{
			if(Idx >= 100)
				{
				if((pBase[-100] & 0x07) == eBaseN)
					NumNs -= 1;
				}
			if((*pBase & 0x07) == eBaseN)
				NumNs += 1;
			if(NumNs > MaxNs)
				break;
			}

		if(NumNs > MaxNs ||
				(P1ReadLen <= 100 && NumNs > ((MaxNs * P1ReadLen)/100)))
			{
			NumExcessNs += 1;
			continue;
			}

		if(bIsPairedEndProc)
			{
			pP2Seq = pP2RawReadsBuff;

			// trim 5' and 3' as requested
			if(Trim5 > 0)
				{
				pP2Seq += Trim5;
				P2ReadLen -= Trim5;
				}

			P2ReadLen -= Trim3;

			// check for excessive number of Ns
			int Idx;
			int NumNs = 0;		// number of indeterminate bases in last 100bp window
			etSeqBase *pBase = pP2Seq;
			for(Idx = 0; Idx < P2ReadLen; Idx++,pBase++)
				{
				if(Idx >= 100)
					{
					if((pBase[-100] & 0x07) == eBaseN)
						NumNs -= 1;
					}
				if((*pBase & 0x07) == eBaseN)
					NumNs += 1;
				if(NumNs > MaxNs)
					break;
				}

			if(NumNs > MaxNs ||
					(P2ReadLen <= 100 && NumNs > ((MaxNs * P2ReadLen)/100)))
				{
				NumExcessNs += 1;
				continue;
				}
			}

		NumReads += 1;

		if((Rslt=AddSeq(P1ReadLen,pP1Seq,P2ReadLen,pP2Seq)) < eBSFSuccess)
			{
			delete pP1RawReadsBuff;
			P1Fasta.Close();
			if(bIsPairedEndProc)
				{
				delete pP2RawReadsBuff;
				P2Fasta.Close();
				}
			return(Rslt);
			}

		Num2Process += 1;

		if(m_P1Seqs2AssembLen > cMaxConcatSeqLen || m_P2Seqs2AssembLen > cMaxConcatSeqLen)			// hit limit?
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Length of loaded concatenated %s sequences exceeds limit of %lld bytes",
										m_P1Seqs2AssembLen > cMaxConcatSeqLen ? "P1" : "P2", cMaxConcatSeqLen);
			delete pP1RawReadsBuff;
			P1Fasta.Close();
			if(bIsPairedEndProc)
				{
				delete pP2RawReadsBuff;
				P2Fasta.Close();
				}
			Reset(false);
			return(eBSFerrMaxDirEls);
			}
		continue;
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"raw sequence file '%s' processing error: %s ",pszP1File,
				Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
		delete pP1RawReadsBuff;
		P1Fasta.Close();
		if(bIsPairedEndProc)
			{
			delete pP2RawReadsBuff;
			P2Fasta.Close();
			}
		Reset(false);
		return(eBSFerrParse);
		}
	}
if(pP1RawReadsBuff != NULL)
	delete pP1RawReadsBuff;

if(bIsPairedEndProc)
	{
	delete pP2RawReadsBuff;
	P2Fasta.Close();
	}

if(Rslt != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors processing file: %s ",pszP1File);
	while(P1Fasta.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,P1Fasta.GetErrMsg());
	P1Fasta.Close();
	return(Rslt);
	}
P1Fasta.Close();

NumAcceptedReads = Num2Process;
m_TotSeqsParsed += NumDescrReads;
m_TotSeqsUnderLen += NumUnderlen;
m_TotSeqsExcessNs += NumExcessNs;
m_TotP1Seqs2Assemb += NumAcceptedReads;
if(bIsPairedEndProc)
	m_TotP2Seqs2Assemb += NumAcceptedReads;
m_P1RdsSfxHdr.RdsSrcFiles[m_P1RdsSfxHdr.NumSrcFiles-1].NumReads = NumAcceptedReads;
if(bIsPairedEndProc)
	m_P2RdsSfxHdr.RdsSrcFiles[m_P2RdsSfxHdr.NumSrcFiles-1].NumReads = NumAcceptedReads;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Completed parsing %d sequences from P1 '%s'",NumDescrReads,pszP1File);
if(bIsPairedEndProc)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: Completed parsing %d sequences from P2 '%s'",NumDescrReads,pszP2File);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadReads: %d were under length %d, %d excessive indeterminate (%d Ns), sequences (after sampling) accepted: %1.9d",NumUnderlen,MinSeqLen,NumExcessNs,MaxNs, NumAcceptedReads);
m_NumRawFiles += 1;
return((teBSFrsltCodes)NumAcceptedReads);
}


// AddSeq
// Add sequences which are to be subsequently suffix indexed and assembled
// basespace only sequences expected
teBSFrsltCodes								// returned SeqID or if <= 0 then error result code
CStackSeqs::AddSeq( int P1SeqLen,			// P1 sequence length
						UINT8 *pP1Seq,		// ptr to P1 sequence
						int P2SeqLen,		// P2 sequence length
						UINT8 *pP2Seq)		// ptr to P1 sequence
{
teBSFrsltCodes Rslt;
UINT64 XFormID;
int Idx;
UINT8 *pDstSeq;
tsSeqStarts *pSeqStarts;

if(P1SeqLen < 1 || pP1Seq == NULL)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"AddSeq: P1 SeqLen: %d pP1Seq is %s",P1SeqLen,pP1Seq == NULL ? "NULL" : "none-null");
	return(eBSFerrParams);
	}

if(m_bIsPairedEndProc && (P2SeqLen < 1 || pP2Seq == NULL))
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"AddSeq: P2 SeqLen: %d pP2Seq is %s",P2SeqLen,pP2Seq == NULL ? "NULL" : "none-null");
	return(eBSFerrParams);
	}

if(m_bIsPairedEndProc)
	{
	if((Rslt = AllocReadsMemory(P1SeqLen+6,P2SeqLen+6)) != eBSFSuccess)		// 6 additional to cover the 5bytes of the XForm'd PairID and the cSeqSep byte
		return(Rslt);
	}
else
	{
	if((Rslt = AllocReadsMemory(P1SeqLen+6,0)) != eBSFSuccess)		// 6 additional to cover the 5bytes of the XForm'd PairID and the cSeqSep byte
		return(Rslt);
	}

pSeqStarts = &m_pSeqStarts[m_NumSeqStarts++];
pSeqStarts->PairID = m_NumSeqStarts;
pSeqStarts->P1SeqLen = P1SeqLen;
pDstSeq = &m_pP1Seqs2Assemb[m_P1Seqs2AssembLen];
if(m_P1Seqs2AssembLen == 0)
	*pDstSeq++ = cCSeqBOS;
else
	*pDstSeq++ = cCSeqSep;
m_P1Seqs2AssembLen += 1;
pSeqStarts->P1SeqOfs = m_P1Seqs2AssembLen;
XFormID = IDtoXForm(m_NumSeqStarts);
*(UINT32 *)pDstSeq = (UINT32)XFormID;
pDstSeq += sizeof(UINT32);
*pDstSeq++ = (UINT8)(XFormID >> 32);
m_P1Seqs2AssembLen += 5;
for(Idx = 0; Idx < (int)P1SeqLen; Idx++)
	*pDstSeq++ = *pP1Seq++ & 0x07;
*pDstSeq++ = cCSeqEOS;		// will be overwritten by next AddEntry()
m_P1Seqs2AssembLen += P1SeqLen+1;
m_NumP1Seqs2Assemb = m_NumSeqStarts;
if(m_bIsPairedEndProc)
	{
	pSeqStarts->P2SeqLen = P2SeqLen;
	pDstSeq = &m_pP2Seqs2Assemb[m_P2Seqs2AssembLen];
	if(m_P2Seqs2AssembLen == 0)
		*pDstSeq++ = cCSeqBOS;
	else
		*pDstSeq++ = cCSeqSep;
	m_P2Seqs2AssembLen += 1;
	pSeqStarts->P2SeqOfs = m_P2Seqs2AssembLen;
	*(UINT32 *)pDstSeq = (UINT32)XFormID;
	pDstSeq += sizeof(UINT32);
	*pDstSeq++ = (UINT8)(XFormID >> 32);
	m_P2Seqs2AssembLen += 5;
	for(Idx = 0; Idx < (int)P2SeqLen; Idx++)
		*pDstSeq++ = *pP2Seq++ & 0x07;
	*pDstSeq++ = cCSeqEOS;		// will be overwritten by next AddEntry()
	m_P2Seqs2AssembLen += P2SeqLen+1;
	m_NumP2Seqs2Assemb = m_NumSeqStarts;
	}
else
	{
	pSeqStarts->P2SeqOfs = 0;
	pSeqStarts->P2SeqLen = 0;
	m_NumP2Seqs2Assemb = 0;
	}
return((teBSFrsltCodes)m_NumSeqStarts);
}


// generate sfx over current concatenated P1 sequences
teBSFrsltCodes
CStackSeqs::GenP1RdsSfx(void)
{
int ElSize;
UINT64 MaxSuffixEls;
UINT32 *pTmp;
MaxSuffixEls = m_P1Seqs2AssembLen;
ElSize = sizeof(UINT32);
if(MaxSuffixEls > cMaxSfxBlkEls)
	ElSize += 1;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenRdsSfx: Now generating suffix array with %lld index elements size %d bytes..",MaxSuffixEls,ElSize);

size_t ReqAllocMem;
ReqAllocMem = (size_t)((MaxSuffixEls + 10) * (INT64)ElSize);
if(m_pSuffixArray == NULL || m_AllocMemSfx == 0)
	{
	m_AllocMemSfx = ReqAllocMem;
#ifdef _WIN32
	m_pSuffixArray = (UINT32 *) malloc((size_t)m_AllocMemSfx);
	if(m_pSuffixArray == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenRdsSfx: Suffix array memory allocation of %lld bytes - %s",(INT64)m_AllocMemSfx,strerror(errno));
		m_AllocMemSfx = 0;
		Reset(false);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pSuffixArray = (UINT32 *)mmap(NULL,m_AllocMemSfx, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSuffixArray == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenRdsSfx: Suffix array memory allocation of %lld bytes through mmap()  failed - %s",(INT64)m_AllocMemSfx,strerror(errno));
		m_pSuffixArray = NULL;
		m_AllocMemSfx = 0;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	SetMaxMemWorkSetSize((size_t)(m_AllocMemP1Seqs2Assemb + m_AllocMemP1Seqs2Assemb + m_AllocMemSfx));
	}
else
	{

	// if worth the cost and suffix array can be reduced then do so - memory could be in short supply!
	if(ReqAllocMem >  (size_t)m_AllocMemSfx || ((ReqAllocMem * 10) <  ((size_t)m_AllocMemSfx * 11)))
		{
#ifdef _WIN32
		pTmp = (UINT32 *) realloc(m_pSuffixArray,ReqAllocMem);
#else
		pTmp = (UINT32 *)mremap(m_pSuffixArray,m_AllocMemSfx,ReqAllocMem,MREMAP_MAYMOVE);
		if(pTmp == MAP_FAILED)
			pTmp = NULL;
#endif
		if(pTmp == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenRdsSfx: Memory re-allocation to %ld bytes - %s",ReqAllocMem,strerror(errno));
			return(eBSFerrMem);
			}
		m_AllocMemSfx = ReqAllocMem;
		m_pSuffixArray = pTmp;
		SetMaxMemWorkSetSize((size_t)(m_AllocMemP1Seqs2Assemb + m_AllocMemP1Seqs2Assemb + m_AllocMemSfx));
		}
	}


if(ElSize == 5)
	{
	UINT8 *pArr5 = (UINT8 *)m_pSuffixArray;
	for(UINT64 El = 0; El < MaxSuffixEls; El++)
		pArr5 = Pack5(El,pArr5);
	*pArr5 = cCSeqEOS;
	m_xpConcatSeqs = &m_pP1Seqs2Assemb[1];
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenRdsSfx: Now sorting...");
	m_MTqsort.qsort(m_pSuffixArray,(INT64)MaxSuffixEls,5,Sfx5SortFunc);
	}
else
	{
	for(UINT32 El = 0; El < MaxSuffixEls; El++)
		m_pSuffixArray[El] = El;
	m_pSuffixArray[MaxSuffixEls] = cCSeqEOS;
	m_xpConcatSeqs = &m_pP1Seqs2Assemb[1];
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenRdsSfx: Now sorting...");
	m_MTqsort.qsort(m_pSuffixArray,(INT64)MaxSuffixEls,sizeof(UINT32),SfxSortFunc);
	}

m_NumSuffixEls = MaxSuffixEls;
m_SfxElSize = ElSize;


#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenRdsSfx: Suffix array generation completed");
return(eBSFSuccess);
}

// ChunkedRead
// Seeks to specified 64bit file offset and reads from disk as chunks of no more than INT_MAX/32
teBSFrsltCodes
CStackSeqs::ChunkedRead(INT64 RdOfs,UINT8 *pData,INT64 RdLen)
{
return(ChunkedRead(m_hInFile,m_szInFile,RdOfs,pData,RdLen));
}

// ChunkedRead
// Seeks to specified 64bit file offset and reads from disk as chunks of no more than INT_MAX/32
teBSFrsltCodes
CStackSeqs::ChunkedRead(int hFile,char *pszFile,INT64 RdOfs,UINT8 *pData,INT64 RdLen)
{
int BlockLen;
if(_lseeki64(hFile,RdOfs,SEEK_SET) != RdOfs)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to seek to %ld on file - error %s",RdOfs,pszFile,strerror(errno));
	return(eBSFerrFileAccess);
	}

while(RdLen)
	{
	BlockLen = RdLen > (INT64)(INT_MAX/32) ? (INT_MAX/32) : (int)RdLen;
	RdLen -= BlockLen;
	if(read(hFile,pData,BlockLen)!=BlockLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to read from disk on file %s - error %s",pszFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	pData += BlockLen;
	}
return(eBSFSuccess);
}

int  // SfxSortFunc for UINT32 suffix elements
CStackSeqs::SfxSortFunc(const void *arg1, const void *arg2)
{
etSeqBase Base1;
etSeqBase Base2;
UINT8 Byte1;
UINT8 Byte2;
int MaxCmpLen = cAbsMaxCoreLen * 500;

UINT8 *pSeq1;
UINT8 *pSeq2;
UINT32 SeqIdx1 = *(UINT32 *)arg1;
UINT32 SeqIdx2 = *(UINT32 *)arg2;
pSeq1 = &m_xpConcatSeqs[SeqIdx1];
pSeq2 = &m_xpConcatSeqs[SeqIdx2];
do {
	Base1 = (Byte1 = *pSeq1++) & 0x07;
	Base2 = (Byte2 = *pSeq2++) & 0x07;
	if(Byte1 >= cCSeqSep)
		Base1 = Byte1;
	if(Byte2 >= cCSeqSep)
		Base2 = Byte2;
	if(Base1 < Base2)
		return(-1);
	if(Base1 > Base2)
		return(1);
	}
while(MaxCmpLen-- && (Base1 < 0x07 && Base2 < 0x07));
return(0);
}


int  // SfxSortFunc for 5byte suffix elements
CStackSeqs::Sfx5SortFunc(const void *arg1, const void *arg2)
{
etSeqBase Base1;
etSeqBase Base2;
UINT8 Byte1;
UINT8 Byte2;
int MaxCmpLen = cAbsMaxCoreLen * 500;

UINT8 *pSeq1;
UINT8 *pSeq2;
UINT64 SeqIdx1 = Unpack5((UINT8 *)arg1);
UINT64 SeqIdx2 = Unpack5((UINT8 *)arg2);
pSeq1 = &m_xpConcatSeqs[SeqIdx1];
pSeq2 = &m_xpConcatSeqs[SeqIdx2];
do {
	Base1 = (Byte1 = *pSeq1++) & 0x07;
	Base2 = (Byte2 = *pSeq2++) & 0x07;
	if(Byte1 >= cCSeqSep)
		Base1 = Byte1;
	if(Byte2 >= cCSeqSep)
		Base2 = Byte2;
	if(Base1 < Base2)
		return(-1);
	if(Base1 > Base2)
		return(1);
	}
while(MaxCmpLen-- && (Base1 < 0x07 && Base2 < 0x07));
return(0);
}


// ReadsSortFunc
// Sort function for sorting reads by sequence
int
CStackSeqs::ReadsSortFunc(const void *arg1, const void *arg2)
{
tsSfxdSeq *pR1 = (tsSfxdSeq *)arg1;
tsSfxdSeq *pR2 = (tsSfxdSeq *)arg2;
etSeqBase Base1;
etSeqBase Base2;
UINT8 Byte1;
UINT8 Byte2;

UINT8 *pSeq1;
UINT8 *pSeq2;
UINT64 SeqIdx1 = pR1->ConcatSeqOfs;
UINT64 SeqIdx2 = pR2->ConcatSeqOfs;
pSeq1 = &m_xpConcatSeqs[SeqIdx1 + 5];			// pts to 1st base of sequence
pSeq2 = &m_xpConcatSeqs[SeqIdx2 + 5];			// pts to 1st base of sequence
do {
	Base1 = (Byte1 = *pSeq1++) & 0x07;
	Base2 = (Byte2 = *pSeq2++) & 0x07;
	if(Byte1 >= cCSeqSep)
		Base1 = Byte1;
	if(Byte2 >= cCSeqSep)
		Base2 = Byte2;
	if(Base1 < Base2)
		return(-1);
	if(Base1 > Base2)
		return(1);
	}
while(Base1 < 0x07 && Base2 < 0x07);
return(0);
}