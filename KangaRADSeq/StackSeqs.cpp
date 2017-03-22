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

CStackSeqs::CStackSeqs(void)
{
m_pSuffixArray = NULL;
m_pP1Seqs2Assemb = NULL;
m_pP2Seqs2Assemb = NULL;
m_pSeqStarts = NULL;
m_pPathIDs = NULL;


m_hInFile = -1;
m_hOutCtgsFile = -1;
m_hOutVCFfile = -1;

m_bMutexesCreated = false;
m_NumThreads = 0;
memset(&m_P1RdsSfxHdr,0,sizeof(m_P1RdsSfxHdr));
memset(&m_P2RdsSfxHdr,0,sizeof(m_P2RdsSfxHdr));
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


CStackSeqs::~CStackSeqs(void)
{
Reset(false);
}

void CStackSeqs::Reset(bool bSync)
{
if(m_hOutCtgsFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutCtgsFile);
#else
		fsync(m_hOutCtgsFile);
#endif
	close(m_hOutCtgsFile);
	m_hOutCtgsFile = -1;
	}

if(m_hOutVCFfile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutVCFfile);
#else
		fsync(m_hOutVCFfile);
#endif
	close(m_hOutVCFfile);
	m_hOutVCFfile = -1;
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

if(m_pP1Seqs2Assemb != NULL)
	{
#ifdef _WIN32
	free(m_pP1Seqs2Assemb);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pP1Seqs2Assemb != MAP_FAILED)
		munmap(m_pP1Seqs2Assemb,m_AllocMemP1Seqs2Assemb);
#endif	
	m_pP1Seqs2Assemb = NULL;
	m_AllocMemP1Seqs2Assemb = 0;
	}

if(m_pP2Seqs2Assemb != NULL)
	{
#ifdef _WIN32
	free(m_pP2Seqs2Assemb);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pP2Seqs2Assemb != MAP_FAILED)
		munmap(m_pP2Seqs2Assemb,m_AllocMemSeqStarts);
#endif	
	m_pP2Seqs2Assemb = NULL;
	m_AllocMemSeqStarts = 0;
	}

if(m_pSeqStarts != NULL)
	{
#ifdef _WIN32
	free(m_pSeqStarts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqStarts != MAP_FAILED)
		munmap(m_pSeqStarts,m_AllocMemP2Seqs2Assemb);
#endif	
	m_pSeqStarts = NULL;
	m_AllocMemP2Seqs2Assemb = 0;
	}

if(m_pPathIDs != NULL)
	{
	free(m_pPathIDs);
	m_pPathIDs = NULL;
	}


if(m_bMutexesCreated)
	{
	DeleteMutexes();
	m_bMutexesCreated = false;
	}

m_AllocMemSfx = 0;
m_MeanReadLen = 0;
m_CurMaxMemWorkSetBytes = 0;

m_TotSeqsParsed = 0;
m_TotSeqsUnderLen = 0;
m_TotSeqsExcessNs = 0;

m_TotP1Seqs2Assemb = 0;
m_NumP1Seqs2Assemb = 0;
m_P1Seqs2AssembLen = 0;
m_AllocMemP1Seqs2Assemb = 0; 

m_TotP2Seqs2Assemb = 0;
m_NumP2Seqs2Assemb = 0;
m_P2Seqs2AssembLen = 0;
m_AllocMemP2Seqs2Assemb = 0; 

m_AllocSeqStarts = 0;			
m_AllocMemSeqStarts = 0;		
m_NumSeqStarts = 0;				


m_NumSuffixEls = 0;
m_szInFile[0] = '\0';
m_AllocdPathIDs = 0;
m_NumPathIDs = 0;

m_NumRawFiles = 0;
m_NumRdsFiles = 0;

m_bIsPairedEndProc = false;

memset(&m_P1RdsSfxHdr,0,sizeof(m_P1RdsSfxHdr));
memset(&m_P2RdsSfxHdr,0,sizeof(m_P2RdsSfxHdr));
}

int 
CStackSeqs::Init(void)
{
Reset(true);
return(CreateMutexes());
}


teBSFrsltCodes 
CStackSeqs::Process(etPMode PMode,						// processing sensitivity mode
				int P1StackEnd,						// P1 stack maximum end float (default is 1, range 0..10)");
				int P1StackDepth,					// P1 stack minimum depth (default is 10, range 5..1000)");
				int P1StackSubRate,					// P1 stack maximum substitution rate (default is 1%, range 0..10%)");
				int P2MinOvrl,						// P2 minimum read overlap (default is 30% of read length, range 10..100)");
				int P2MaxOvrlSubRate,				// P2 maximum read overlap substitution rate (default is 5%%, range 0..10%%)");
				char *pszCtgDescr,					// generated contig descriptor prefix 
				char *pszOutCtgsFile,				// assembled contigs written to this file
				char *pszOutVCFfile,				// output Variant Call Format (VCF 4.1) to this file
				int NumThreads)					// max number of worker threads to use
{

// generate suffix array over the P1 reads
m_P1RdsSfxHdr.ConcatSeqLen = m_P1Seqs2AssembLen;
GenP1RdsSfx();

// process the P1 reads and find overlaps with other P1 reads

// with each overlap check if it's P2 reads also overlap


return(eBSFSuccess);
}