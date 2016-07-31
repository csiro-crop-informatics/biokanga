// Copyright 2013, 2014 CSIRO  ( http://www.csiro.au/ ) 
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

// PacBioUtility.cpp : utility class for PacBio processing
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

#include "../libbiokanga/bgzf.h"
#include "PacBioUtility.h"

CPacBioUtility::CPacBioUtility()
{
m_pQuerySeqs = NULL;

m_bMutexesCreated = false;
Reset();
}


CPacBioUtility::~CPacBioUtility()
{
Reset();
}

void 
CPacBioUtility::Reset(void)
{
if(m_pQuerySeqs != NULL)
	{
	delete m_pQuerySeqs;
	m_pQuerySeqs = NULL;
	}

m_TotSeqIDs = 0;
m_NumQuerySeqs = 0;
m_NxtQuerySeqIdx = 0;
m_AllocdQuerySeqs = 0;
m_bAllQuerySeqsLoaded = false;
m_LoadQuerySeqsRslt = eBSFSuccess;
m_ThreadLoadQuerySeqsRslt= 0;
m_NumInputFiles = 0;
m_ppszInputFiles = NULL;
m_szAsyncReadsFile[0]= '\0';
m_bAsyncLoading = false;

if(m_bMutexesCreated == true)
	{
	DeleteMutexes();
	m_bMutexesCreated = false;
	}
m_TermBackgoundThreads = 0;

}


int
CPacBioUtility::CreateMutexes(void)
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

m_bMutexesCreated = true;
return(eBSFSuccess);
}

void
CPacBioUtility::DeleteMutexes(void)
{
if(!m_bMutexesCreated)
	return;
#ifdef _WIN32
CloseHandle(m_hMtxIterReads);
#else
pthread_mutex_destroy(&m_hMtxIterReads);
pthread_rwlock_destroy(&m_hRwLock);
#endif
m_bMutexesCreated = false;
}

void
CPacBioUtility::AcquireSerialise(void)
{
#ifdef _WIN32
WaitForSingleObject(m_hMtxIterReads,INFINITE);
#else
pthread_mutex_lock(&m_hMtxIterReads);
#endif
}

void
CPacBioUtility::ReleaseSerialise(void)
{
#ifdef _WIN32
ReleaseMutex(m_hMtxIterReads);
#else
pthread_mutex_unlock(&m_hMtxIterReads);
#endif
}

void
CPacBioUtility::AcquireLock(bool bExclusive)
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
CPacBioUtility::ReleaseLock(bool bExclusive)
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


teBSFrsltCodes 
CPacBioUtility::StartAsyncLoadSeqs(int NumInputFiles,			// number of input files
									char **ppszInputFiles,		// names (wildcards allowed) of input files containing reads to be filtered
									int MaxReadahead)	 // read ahead for at most this number of sequences
{
if(NumInputFiles < 1 || ppszInputFiles == NULL || ppszInputFiles[0] == NULL) 
	return(eBSFerrParams);

if(MaxReadahead < 100)
	MaxReadahead = 100;
else
	if(MaxReadahead > 100000)
		MaxReadahead = 100000;

// can't start a new async load if already loading
if(m_bAsyncLoading)
	return(eBSFerrInternal);

Reset();
m_bAsyncLoading = true;
if((m_pQuerySeqs = new tsQuerySeq [MaxReadahead]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to allocate memory for query sequences");
	Reset();
	return(eBSFerrMem);
	}
m_AllocdQuerySeqs = MaxReadahead;
m_NumQuerySeqs = 0;
m_NxtQuerySeqIdx = 0;
m_TotSeqIDs = 0;

if(CreateMutexes()!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to create thread synchronisation mutexes");
	Reset();
	return(cBSFSyncObjErr);
	}

m_NumInputFiles = NumInputFiles;
m_ppszInputFiles = ppszInputFiles;

return((teBSFrsltCodes)InitLoadQuerySeqs());
}

#ifdef _WIN32
unsigned __stdcall LoadReadsFileThread(void * pThreadPars)
#else
void *LoadReadsFileThread(void * pThreadPars)
#endif
{
int Rslt;
tsLoadQuerySeqsThreadPars *pPars = (tsLoadQuerySeqsThreadPars *)pThreadPars;			// makes it easier not having to deal with casts!
CPacBioUtility *pPacBioUtility = (CPacBioUtility *)pPars->pThis;

Rslt = pPacBioUtility->ProcLoadReadsFile(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int
CPacBioUtility::InitLoadQuerySeqs(void)
{
tsLoadQuerySeqsThreadPars ThreadPars;

// initiate loading the reads
m_ThreadLoadQuerySeqsRslt = -1;

ThreadPars.pRslt = &m_ThreadLoadQuerySeqsRslt;
ThreadPars.pThis = this;
ThreadPars.Rslt = 0;

#ifdef _WIN32
m_hThreadLoadQuerySeqs = ThreadPars.threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,LoadReadsFileThread,&ThreadPars,0,&m_ThreadLoadQuerySeqsID);
#else
int ThreadRslt = ThreadPars.threadRslt = pthread_create (&m_ThreadLoadQuerySeqsID , NULL , LoadReadsFileThread , &ThreadPars );
#endif

// wait a few seconds, if major problems with loading reads then should show very quickly
#ifdef _WIN32
if(WAIT_TIMEOUT != WaitForSingleObject(m_hThreadLoadQuerySeqs, 10000))
	{
	CloseHandle(m_hThreadLoadQuerySeqs);
	m_hThreadLoadQuerySeqs = NULL;
	m_bAsyncLoading = false;
	return(m_ThreadLoadQuerySeqsRslt);
	}
#else
struct timespec ts;
int JoinRlt;
clock_gettime(CLOCK_REALTIME, &ts);
ts.tv_sec += 3;
if((JoinRlt = pthread_timedjoin_np(m_ThreadLoadQuerySeqsID, NULL, &ts)) == 0)
	{
	m_ThreadLoadQuerySeqsID = 0;
	return(m_ThreadLoadQuerySeqsRslt);
	}
#endif
return(eBSFSuccess);
}

int									// returned enqueued query identifier
CPacBioUtility::EnqueueQuerySeq(char *pszQueryIdent,    // query identifier parsed from fasta descriptor
			int QuerySeqLen,		// query sequence length
			UINT8 *pQuerySeq)       // query sequence
{
int Idx;
int SeqID;
UINT8 *pSeq;
tsQuerySeq *psQuery;
if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
	return(0);
if((pSeq = new UINT8 [QuerySeqLen]) == NULL)
	return(eBSFerrMem);
memcpy(pSeq,pQuerySeq,QuerySeqLen);

// any room left in query sequence queue?
while(1) {
	AcquireLock(true);
	if((m_NumQuerySeqs + 1) < m_AllocdQuerySeqs)
		break;
	ReleaseLock(true);
	if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
		{
		delete pSeq;
		return(0);
		}
#ifdef _WIN32
	Sleep(1000);
#else
	sleep(1);
#endif
	}
Idx = (m_NxtQuerySeqIdx + m_NumQuerySeqs) % m_AllocdQuerySeqs;
psQuery = &m_pQuerySeqs[Idx];
SeqID = ++m_TotSeqIDs;
psQuery->SeqID = SeqID;
psQuery->pQuerySeq = pSeq;
psQuery->QuerySeqLen = QuerySeqLen; 
strncpy(psQuery->szQueryIdent,pszQueryIdent,cMaxQuerySeqIdentLen);
psQuery->szQueryIdent[cMaxQuerySeqIdentLen] = '\0';
m_NumQuerySeqs += 1;
ReleaseLock(true);
return(SeqID);
}

UINT8 *												// returned dequeued sequence, caller is responsible for deleting memory allocated to hold the returned sequence (delete pRetSeq;)
CPacBioUtility::DequeueQuerySeq(int WaitSecs,		// if no sequences available to be dequeued then wait at most this many seconds for a sequence to become available
			int MaxLenQueryIdent,			// maximum length query identifier
			int *pSeqID,					// returned sequence identifier
			char *pszQueryIdent,			// where to return query identifier
			int *pQuerySeqLen)				// where to return query sequence length
{
bool bAllQuerySeqsLoaded;
UINT8 *pSeq;
tsQuerySeq *psQuery;

// any sequences available to be dequeued?
if(WaitSecs < 0)
	WaitSecs = 0;
while(1) {
	if(m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
		return(NULL);
	AcquireLock(true);
	if(m_NumQuerySeqs)
		break;
	bAllQuerySeqsLoaded = m_bAllQuerySeqsLoaded;
	ReleaseLock(true);
	if(bAllQuerySeqsLoaded || WaitSecs-- <= 0 || m_TermBackgoundThreads != 0)	// need to immediately self-terminate?
		return(NULL);
#ifdef _WIN32
	Sleep(1000);
#else
	sleep(1);
#endif
	}
psQuery = &m_pQuerySeqs[m_NxtQuerySeqIdx++];
if(m_NxtQuerySeqIdx == m_AllocdQuerySeqs)
	m_NxtQuerySeqIdx = 0;
*pSeqID = psQuery->SeqID;
pSeq = psQuery->pQuerySeq;
psQuery->pQuerySeq = NULL;
*pQuerySeqLen = psQuery->QuerySeqLen; 
strncpy(pszQueryIdent,psQuery->szQueryIdent,MaxLenQueryIdent);
pszQueryIdent[MaxLenQueryIdent-1] = '\0';
m_NumQuerySeqs -= 1;
ReleaseLock(true);
return(pSeq);
}

int
CPacBioUtility::ProcLoadReadsFile(tsLoadQuerySeqsThreadPars *pPars)
{
CFasta Fasta;
unsigned char *pSeqBuff;
unsigned char *pMskBase;
UINT32 MskIdx;
size_t BuffOfs;
size_t AllocdBuffSize;
size_t AvailBuffSize;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
UINT32 SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
bool bTruncSeq;
int Rslt;
int SeqID;
int FileIdx;
int *pRslt = pPars->pRslt;
AcquireLock(true);
m_bAllQuerySeqsLoaded = false;
m_LoadQuerySeqsRslt = eBSFSuccess;		// presumed success, changed if any processing errors
ReleaseLock(true);

char *pszInFile;
CSimpleGlob glob(SG_GLOB_FULLSORT);
for(FileIdx = 0; FileIdx < m_NumInputFiles; FileIdx++)
	{
	pszInFile = NULL;
	glob.Init();
	if(glob.Add(m_ppszInputFiles[FileIdx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"ProcLoadReadsFile: Unable to glob '%s",m_ppszInputFiles[FileIdx]);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	if(glob.FileCount() <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcLoadReadsFile:  Unable to locate any source reads file matching '%s",m_ppszInputFiles[FileIdx]);
		AcquireLock(true);
		m_bAllQuerySeqsLoaded = true;
		m_LoadQuerySeqsRslt = eBSFerrFileAccess;
		m_bAsyncLoading = false;
		ReleaseLock(true);
		return(eBSFerrOpnFile);
		}
	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
		{
		pszInFile = glob.File(FileID);
		strcpy(m_szAsyncReadsFile,pszInFile);

		gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcLoadReadsFile: Loading reads from '%s'",m_szAsyncReadsFile);

		if((Rslt=Fasta.Open(m_szAsyncReadsFile,true))!=eBSFSuccess)
			{
			if(Rslt != eBSFerrNotFasta)
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcLoadReadsFile: Unable to open '%s' [%s] %s",m_szAsyncReadsFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
			AcquireLock(true);
			m_bAllQuerySeqsLoaded = true;
			m_LoadQuerySeqsRslt = (teBSFrsltCodes)Rslt;
			m_bAsyncLoading = false;
			ReleaseLock(true);
			return(Rslt);
			}

		// note malloc is used as can then simply realloc to expand as may later be required
		AllocdBuffSize = (size_t)cAllocQuerySeqLen;
		if((pSeqBuff = (unsigned char *)malloc(AllocdBuffSize)) == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcLoadReadsFile:- Unable to allocate memory (%u bytes) for sequence buffer",(UINT32)cAllocQuerySeqLen);
			Fasta.Close();
			*pRslt = eBSFerrMem;
			AcquireLock(true);
			m_bAllQuerySeqsLoaded = true;
			m_LoadQuerySeqsRslt = (teBSFrsltCodes)Rslt;	
			m_bAsyncLoading = false;	
			ReleaseLock(true);
			return(eBSFerrMem);
			}
		AvailBuffSize = cAllocQuerySeqLen;

		bFirstEntry = true;
		bEntryCreated = false;
		bTruncSeq = false;
		SeqID = 0;
		BuffOfs = 0;
		while((Rslt = SeqLen = Fasta.ReadSequence(&pSeqBuff[BuffOfs],(int)min(AvailBuffSize,(size_t)cMaxQuerySeqLen),true,false)) > eBSFSuccess)
			{
			if(m_TermBackgoundThreads != 0)	// requested to immediately self-terminate?
				{
				Rslt = eBSErrSession;
				break;
				}

			if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
				{
				SeqID++;
				if(bEntryCreated)				// add any previous entry
					{
					if((Rslt=EnqueueQuerySeq(szName,(int)BuffOfs,pSeqBuff)) <= eBSFSuccess)
						break;
					}
				Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
				// An assumption - will one day bite real hard - is that the
				// fasta descriptor line starts with some form of unique identifier.
				// Use this identifier as the entry name.
				if(sscanf(szDescription," %s[ ,]",szName)!=1)
					sprintf(szName,"%s.%d",m_szAsyncReadsFile,++SeqID);

				bFirstEntry = false;
				bEntryCreated = true;
				bTruncSeq = false;
				BuffOfs = 0;
				continue;
				}
			else
				if(bFirstEntry)	// if there was no descriptor then dummy up one...
					{
					SeqID++;
					sprintf(szName,"%s.%d",m_szAsyncReadsFile,SeqID);
					strcpy(szDescription,"No Description provided");
					bFirstEntry = false;
					bEntryCreated = true;
					bTruncSeq = false;
					}
			if(bTruncSeq)
				continue;

			// remove any repeat masking flags and translate any noncannonical bases into a base which is sequence offset dependent - not randomised but will be reproducible 
			pMskBase = &pSeqBuff[BuffOfs];
			for(MskIdx = 0; MskIdx < SeqLen; MskIdx++,pMskBase++)
				{
				*pMskBase &= ~cRptMskFlg;
				if(*pMskBase > eBaseT)
					*pMskBase = MskIdx & 0x03;
				}

			BuffOfs += SeqLen;

			if(BuffOfs > cMaxQuerySeqLen)	// truncate at cMaxQuerySeqLen
				{
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"ProcLoadReadsFile:- Truncating overlength query sequence '%s' to %d",szName,cMaxQuerySeqLen);
				BuffOfs = cMaxQuerySeqLen;
				AvailBuffSize = AllocdBuffSize - BuffOfs;
				bTruncSeq = true;
				continue;
				}
			AvailBuffSize -= SeqLen;

			if(AvailBuffSize < (size_t)(cAllocQuerySeqLen / 2))
				{
				size_t NewSize = (size_t)cAllocQuerySeqLen + AllocdBuffSize;
				unsigned char *pTmp;
				if((pTmp = (unsigned char *)realloc(pSeqBuff,NewSize))==NULL)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcLoadReadsFile:- Unable to reallocate memory (%u bytes) for sequence buffer",(UINT32)NewSize);
					Rslt = eBSFerrMem;
					break;
					}
				pSeqBuff = pTmp;
				AllocdBuffSize = NewSize;
				AvailBuffSize = AllocdBuffSize - BuffOfs;
				}
			}
		if(Rslt < eBSFSuccess && Rslt != eBSErrSession)
			{
			while(Fasta.NumErrMsgs())
				gDiagnostics.DiagOut(eDLFatal,gszProcName,Fasta.GetErrMsg());
			}
		Fasta.Close();
		if(Rslt >= eBSFSuccess && bEntryCreated && BuffOfs > 0)			// last entry
			Rslt=EnqueueQuerySeq(szName,(int)BuffOfs,pSeqBuff);
		if(Rslt < eBSFSuccess)
			break;
		}
	}

if(pSeqBuff != NULL)
	free(pSeqBuff);
*pRslt = Rslt;
AcquireSerialise();
m_bAllQuerySeqsLoaded = true;
m_LoadQuerySeqsRslt = (teBSFrsltCodes)Rslt;
m_bAsyncLoading = false;
ReleaseSerialise();
return(Rslt);
}



// Identify PacBio ISO-seq primers at 5' and 3' ends of error corrected reads (still being developed!!!! )
int												// identified SMRTBell at this offset + 1, 0 if non-detected
CPacBioUtility::DetectIsoSeqPrimers(int StartOfs,	// search from this offset
					int *pNumTetramers,			// returned number of tetramers from SmartBell sequence detected and which were in expected order 
					int InSeqLen,				// number of bases in sequence to search for SmartBell
					etSeqBase *pInSeq,			// sequence to search for SmartBell
					int MinTetramers)			// only report SmartBell if at least this many tetramers in expected order detected
{
int Idx;
int TIdx;
int QIdx;
const etSeqBase *pQSeq;
etSeqBase *pTSeq;
int QKMers;
int PeakQKMers;
int PeakLoc;

if(StartOfs < 0 || pNumTetramers == NULL || InSeqLen < cSmartBellAdapterSeqLen || pInSeq == NULL || MinTetramers < 1 ||
	StartOfs + cSmartBellAdapterSeqLen > InSeqLen)
	return(eBSFerrParams);

PeakQKMers = 0;
PeakLoc = 0;
for(Idx = StartOfs+6; Idx < InSeqLen - cSmartBellAdapterSeqLen/3 ; Idx++)
	{
	pQSeq = &cSmartBellAdapterSeq[0];
	QKMers = 0;
	for(QIdx = 0; QIdx < cSmartBellAdapterSeqLen - 6; QIdx++,pQSeq++)
		{
		pTSeq = &pInSeq[Idx-6];
		for(TIdx = 0; TIdx < 12; TIdx++,pTSeq++)
			{
			if((*pQSeq == (*pTSeq & 0x03)) && 
				(pQSeq[1] == (pTSeq[1] & 0x03)) &&
				(pQSeq[2] == (pTSeq[2] & 0x03)) &&
				(pQSeq[3] == (pTSeq[3] & 0x03))) 
				{
				QKMers += 1;
				break;
				}
			}
		}
	if(QKMers > PeakQKMers)
		{
		PeakQKMers = QKMers;
		PeakLoc = Idx;
		}

	if(QKMers >= MinTetramers)
		{
		*pNumTetramers = QKMers;
		return(Idx);
		}
	}
*pNumTetramers = 0;
return(0);
}


int		// length after hompolymer reduction
CPacBioUtility::ReduceHomopolymers(int InSeqLen,	// number of bases in sequence
					   etSeqBase *pInSeq,			// sequence which may contain homopolymers
					   int TruncLen)				// homopolymers longer than this length (min 12) to be truncated at this length
{
UINT8 *pHomoBase;
UINT8 *pDelHomoBase;
UINT8 HomoBase;
UINT8 CurHomoBase;
UINT32 TotHomoLen;
UINT32 HomoIdx;
UINT32 StartHomoIdx;
bool bHomoRuns;
if(TruncLen < 12)
	return(InSeqLen);

// checking for near homopolymer runs and marking these bases for deletion
TotHomoLen = 0;
StartHomoIdx = 0;
bHomoRuns = false;
pHomoBase = pInSeq;
for(HomoIdx = 0; HomoIdx < (UINT32)InSeqLen; HomoIdx++,pHomoBase++)
	{
	HomoBase = *pHomoBase & 0x03;
	if(TotHomoLen == 0)
		{
		StartHomoIdx = HomoIdx;
		CurHomoBase = HomoBase;
		TotHomoLen = 1;
		}
	else
		{
		if(HomoBase == CurHomoBase)
			TotHomoLen += 1;
		else
			{
			if(TotHomoLen >= 12)				// requiring an initial seed K-mer of at least 12 bases before accepting as possibly a homopolymer
				{
				// although not matching the current homopolymer base if the next three bases would match then accept as still being part of the homopolymer
				if((HomoIdx + 4) < (UINT32)InSeqLen)
					{ 
					if(((pHomoBase[1] & 0x03) == CurHomoBase) && ((pHomoBase[2] & 0x03) == CurHomoBase) && ((pHomoBase[3] & 0x03) == CurHomoBase))
						{
						TotHomoLen += 1;
						continue;
						}
					}

				if(TotHomoLen >= (UINT32)TruncLen)				// accepting as homopolymer if at least m_FiltMinHomoLen long
					{
					pDelHomoBase = &pInSeq[StartHomoIdx+6];	// retaining the first 6 bases as these started the homopolymer run
					TotHomoLen -= 6;							// and retaining the last 6 bases as these terminated the homopolymer run
					while(TotHomoLen--) 
						*pDelHomoBase++ = eBaseInDel;			// marking base for subsequent deletion 
					bHomoRuns = true;
					}
				}
			TotHomoLen = 0;
			StartHomoIdx = 0;
			}
		}
 	}

if(TotHomoLen >= (UINT32)TruncLen)			// accepting as homopolymer if at least m_FiltMinHomoLen long
	{
	pDelHomoBase = &pInSeq[StartHomoIdx+6];	// retaining the first 5 bases as these started the homopolymer run
	TotHomoLen -= 6;					    // and retaining the last 5 bases as these terminated the homopolymer run
	while(TotHomoLen--) 
		*pDelHomoBase++ = eBaseInDel;      // marking base for subsequent deletion 
	bHomoRuns = true;
	}

if(bHomoRuns)
	{
	TotHomoLen = 0;
	pDelHomoBase = pHomoBase = pInSeq;
	for(HomoIdx = 0; HomoIdx < (UINT32)InSeqLen; HomoIdx++,pHomoBase++)
		{
		HomoBase = *pHomoBase & 0x07;
		if(HomoBase == eBaseInDel)
			continue;
		*pDelHomoBase++ = *pHomoBase;
		TotHomoLen += 1;
		}
	InSeqLen = (int)TotHomoLen;
	}

return(InSeqLen);
}

int													// number of putatve SMRTBell hairpins identified
CPacBioUtility::IdentifySMRTBells(int SMRTBellSensitivity, // sensitivity of SMRTBell detection - 5: max, 1: min
						int MaxSMRTBells,	// identify at most this many SMRTBells to return in pSMRTBells
						int SeqLen,					// length of sequence to search for SMRTBell hairpins
						etSeqBase *pSeq,			// identify all putative SMRTBell hairpins in this sequence
						tsSMRTBellHit *pSMRTBells)	// returned identified SMRTBells
{
int NumSMRTBells;
int LHE;
int RHS;
int Len;
int MinScore;
int Hits;
int NumTargsSense;
int NumTargsAnti;

int	MinSMRTBellScore;
int	MinFlankScore;

tsSMRTBellHit *pExchg;
tsSMRTBellHit Swap;
int NumSwaps;
int ExchgIdx;

tsSWQuickHit *pQuickHit;
tsSWQuickHit QuickHitsSense[50];
tsSWQuickHit QuickHitsAnti[50];

tsSWQuickHit QuickHitsExtd[50];


int NumTargsFSense;
int NumTargsFAnti;

etSeqBase AdapterAntisense[1000];
etSeqBase RHSAntisense[1000];

switch(SMRTBellSensitivity) {
	case 0: case 1:		// minimum sensitivity, maximum specificity
		MinSMRTBellScore = 17;
		MinFlankScore = 150;
		break;

	case 2:
		MinSMRTBellScore = 16;
		MinFlankScore = 125;
		break;

	case 3:						// default
		MinSMRTBellScore = 15;
		MinFlankScore = 100;
		break;

	case 4:
		MinSMRTBellScore = 14;
		MinFlankScore = 87;
		break;

	default:			// any other sensitivity treated as being maximal sensitivity but minimal specificity
		MinSMRTBellScore = 13;
		MinFlankScore = 75;
		break;
	}

NumTargsFSense = 0;
NumTargsFAnti = 0;
NumSMRTBells = 0;
NumTargsSense = SWQuick(cSmartBellAdapterSeqLen,(etSeqBase *)cSmartBellAdapterSeq,SeqLen,pSeq,15,QuickHitsSense,MinSMRTBellScore,2,-20,-4,-3,-4,-3);
if(NumTargsSense > 0)
	{
	pQuickHit = QuickHitsSense;
	for(int Idx = 0; Idx < NumTargsSense; Idx++,pQuickHit+=1)
		{
		LHE = pQuickHit->TargOfs - pQuickHit->QueryOfs;
		RHS = LHE + cSmartBellAdapterSeqLen;
		if(LHE < 100 || RHS+100 >= SeqLen)
			continue;
		Len = min(500,min(LHE,SeqLen-RHS));
		MinScore = (MinFlankScore * Len) / 500;
		memcpy(RHSAntisense,&pSeq[RHS],Len);
		CSeqTrans::ReverseComplement(Len,RHSAntisense);
		if((Hits = SWQuick(Len,&pSeq[LHE-Len],Len,RHSAntisense,1,QuickHitsExtd,MinScore,2,-20,-4,-3,-4,-3)) > 0)
			{
			pSMRTBells[NumSMRTBells].HiScore = QuickHitsExtd[0].HiScore;
			pSMRTBells[NumSMRTBells].TargOfs = (LHE + RHS)/2;
			NumSMRTBells += 1;
			MaxSMRTBells -= 1;
			// not expecting too many SMRTBells in any targeted sequence
			// ensure SMRTBells are ordered from 5' end ascending just using a simple exchange sort
			if(NumSMRTBells > 1)
				{
				do {
					pExchg = pSMRTBells;
					NumSwaps = 0;
					for(ExchgIdx = 0; ExchgIdx < NumSMRTBells - 1; ExchgIdx++, pExchg++)
						{
						if(pExchg->TargOfs > pExchg[1].TargOfs)
							{
							Swap = pExchg[1];
							pExchg[1] = *pExchg;
							*pExchg = Swap;
							NumSwaps += 1;
							}
						}
					}
				while(NumSwaps);
				}
			}
		if(MaxSMRTBells == 0)
			break;
		}
	}

if(MaxSMRTBells > 0)
	{
	memcpy(AdapterAntisense,cSmartBellAdapterSeq,cSmartBellAdapterSeqLen);
	CSeqTrans::ReverseComplement(cSmartBellAdapterSeqLen,AdapterAntisense);
	NumTargsAnti = SWQuick(cSmartBellAdapterSeqLen,(etSeqBase *)AdapterAntisense,SeqLen,pSeq,15,QuickHitsAnti,MinSMRTBellScore,2,-20,-4,-3,-4,-3);
	if(NumTargsAnti > 0)
		{
		pQuickHit = QuickHitsAnti;
		for(int Idx = 0; Idx < NumTargsAnti; Idx++,pQuickHit+=1)
			{
			LHE = pQuickHit->TargOfs - pQuickHit->QueryOfs;
			RHS = LHE + cSmartBellAdapterSeqLen;
			if(LHE < 100 || RHS+100 >= SeqLen)
				continue;
			Len = min(500,min(LHE,SeqLen-RHS));
			MinScore = (MinFlankScore * Len) / 500;
			memcpy(RHSAntisense,&pSeq[RHS],Len);
			CSeqTrans::ReverseComplement(Len,RHSAntisense);
			if((Hits = SWQuick(Len,&pSeq[LHE-Len],Len,RHSAntisense,1,&QuickHitsExtd[NumSMRTBells],MinScore,2,-20,-4,-3,-4,-3)) > 0)
				{
				pSMRTBells[NumSMRTBells].HiScore = QuickHitsExtd[0].HiScore;
				pSMRTBells[NumSMRTBells].TargOfs = (LHE + RHS)/2;
				NumSMRTBells += 1;
				MaxSMRTBells -= 1;
				// not expecting too many SMRTBells in any targeted sequence
				// ensure SMRTBells are ordered from 5' end ascending just using a simple exchange sort
				if(NumSMRTBells > 1)
					{
					do {
						pExchg = pSMRTBells;
						NumSwaps = 0;
						for(ExchgIdx = 0; ExchgIdx < NumSMRTBells - 1; ExchgIdx++, pExchg++)
							{
							if(pExchg->TargOfs > pExchg[1].TargOfs)
								{
								Swap = pExchg[1];
								pExchg[1] = *pExchg;
								*pExchg = Swap;
								NumSwaps += 1;
								}
							}
						}
					while(NumSwaps);
					}
				}
			if(MaxSMRTBells == 0)
				break;
			}
		}
	}

return(NumSMRTBells);
}


typedef struct TAG_sSWQScore {
	INT32 PathScore;		// accumulated path score
	UINT32 RunLen;			// current run length of Class type
	UINT32 RunClass;		// run length class; 0:exact match, 1: mismatches, 2: insertions into the target, 3: deletions from the target		
} tsSWQScore;

int								// number of target hits returned
CPacBioUtility::SWQuick(int QueryLen,			// relatively short query sequence to search for; typically a SMRTBell or some adapter sequence
		etSeqBase *pQuerySeq,	// query sequence
		int TargLen,			// length of sequence to search, expected to be at least the query sequence length
		etSeqBase *pTargSeq,	// identify all queries in this target sequence
		int MaxTargHits,		// return at most this many target hits
		tsSWQuickHit *pTargHits, // returned putative target hits
		int MinScore,			// only report hits having at least this score
		int MatchScore,			// score for exact base matches
		int MismatchScore,		// score for mismatches
		int InsertScore,		// score for 1st inserted base of any contiguous run of insertions into the target
		int InsertAffineScore,  // score for 2nd and subsequent inserted bases in any contiguous insertion run
		int DeleteScore,		// score for 1st deleted base of any contiguous deletions from the target
		int DeleteAffineScore)  // score for 2nd and subsequent deleted bases in any contiguous deletion run
{
static bool bFirst = true;
bool bExactMatch;
int HiScore;
int LoHiScore;
int LoHitIdx;
int HiQIdx;
int NumHiScores;
tsSWQScore CellScores[cMaxSWQuickQueryLen+1];
tsSWQScore *pCell;
tsSWQScore DiagCell;
tsSWQScore NxtDiagCell;
tsSWQScore NewMatchCell;
tsSWQScore NewInsertCell;
tsSWQScore NewDeleteCell;

int TLen, QLen;
int TIdx,QIdx;
etSeqBase *pTSeq,*pQSeq;
tsSWQuickHit *pHit;
int HitIdx;

int MaxScore;

if(QueryLen < cMinSWQuickQueryLen || QueryLen > cMaxSWQuickQueryLen || pQuerySeq == NULL ||
	TargLen < QueryLen || pTargSeq == NULL || 
	MaxTargHits < 1 || MaxTargHits > cMaxSWQuickHits || pTargHits == NULL || MinScore == 0)
	return(eBSFerrParams);

memset(pTargHits,0,sizeof(tsSWQuickHit) * MaxTargHits);
memset(CellScores,0,sizeof(tsSWQScore) * (QueryLen + 1));
NumHiScores = 0;
MaxScore = 0;
TLen = TargLen;
pTSeq = pTargSeq;
for(TIdx = 0; TIdx < TLen; TIdx++,pTSeq++)
	{
	pQSeq = pQuerySeq;
	QLen = QueryLen;
	pCell = &CellScores[1];
	DiagCell = CellScores[0];
	HiScore = 0;
	HiQIdx = -1;
	for(QIdx = 0; QIdx < QLen; QIdx++,pQSeq++,pCell++)
		{
		NxtDiagCell = *pCell;
		bExactMatch = (*pTSeq & 0x03) == (*pQSeq & 0x03) ? true : false;
		NewMatchCell = DiagCell;
		if(bExactMatch) 
			{
			NewMatchCell.PathScore += MatchScore;
			if(NewMatchCell.RunClass != 0)
				{
				NewMatchCell.RunLen = 1;
				NewMatchCell.RunClass = 0;
				}
			else 
				NewMatchCell.RunLen += 1;
			}
		else	
			{
			NewMatchCell.PathScore += MismatchScore;
			if(NewMatchCell.RunClass != 1)
				{
				NewMatchCell.RunLen = 1;
				NewMatchCell.RunClass = 1;
				}
			else
				NewMatchCell.RunLen += 1;
			}

		NewInsertCell = *pCell;
		if(NewInsertCell.RunClass != 2)
			{
			NewInsertCell.RunLen = 1;
			NewInsertCell.RunClass = 2;
			NewInsertCell.PathScore += InsertScore;
			}
		else
			{
			NewInsertCell.RunLen += 1;
			NewInsertCell.PathScore += InsertAffineScore;
			}

		NewDeleteCell = pCell[-1];
		if(NewDeleteCell.RunClass != 3)
			{
			NewDeleteCell.RunLen = 1;
			NewDeleteCell.RunClass = 3;
			NewDeleteCell.PathScore += DeleteScore;
			}
		else
			{
			NewDeleteCell.RunLen += 1;
			NewDeleteCell.PathScore += DeleteAffineScore;
			}

		if(NewMatchCell.PathScore <= 0)
			{
			NewMatchCell.PathScore = 0;
			NewMatchCell.RunClass = 0;
			NewMatchCell.RunLen = 0;
			}

		if(NewInsertCell.PathScore <= 0)
			{
			NewInsertCell.PathScore = 0;
			NewInsertCell.RunClass = 0;
			NewInsertCell.RunLen = 0;
			}

		if(NewDeleteCell.PathScore <= 0)
			{
			NewDeleteCell.PathScore = 0;
			NewDeleteCell.RunClass = 0;
			NewDeleteCell.RunLen = 0;
			}

		if(NewMatchCell.PathScore >= NewInsertCell.PathScore && NewMatchCell.PathScore >= NewDeleteCell.PathScore)
			*pCell = NewMatchCell;
		else
			{
			if(NewInsertCell.PathScore >= NewDeleteCell.PathScore)
				*pCell = NewInsertCell;
			else
				*pCell = NewDeleteCell;
			}
		DiagCell = NxtDiagCell;
		if(pCell->PathScore > HiScore)
			{
			HiScore = pCell->PathScore;
			HiQIdx = QIdx;
			}
		}
	if(HiScore > MaxScore)
		MaxScore = HiScore;

	if(HiScore >= MinScore)
		{
		// overwrite any previous highest score if that highest score was <= HiScore and within QueryLen/2 of this latest score
		// if a previous highest score is higher and within QueryLen/2 of this latest score then slough this latest HiScore
		if(NumHiScores)
			{
			pHit = &pTargHits[NumHiScores-1];
			if(pHit->TargOfs >= (TIdx - QueryLen/2))
				{
				if(HiScore >= pHit->HiScore)
					{
					pHit->TargOfs = TIdx;
					pHit->QueryOfs = HiQIdx;
					pHit->HiScore = HiScore;
					}
				continue;
				}
			if(NumHiScores == MaxTargHits)	// if score higher than lowest previously recorded then slough that lowest score to make room for this higher score
				{
				LoHiScore = HiScore;
				LoHitIdx = -1;
				pHit =  pTargHits;
				for(HitIdx = 0; HitIdx < MaxTargHits; HitIdx++,pHit++)
					{
					if(pHit->HiScore < LoHiScore)
						{
						LoHitIdx = HitIdx;
						LoHiScore = pHit->HiScore;
						}
					}
				if(LoHitIdx < 0)
					continue;
				NumHiScores -= 1;
				pHit = &pTargHits[LoHitIdx];
				if(LoHitIdx < NumHiScores)
					memmove(pHit,&pHit[1],(NumHiScores-LoHitIdx)  * sizeof(tsSWQuickHit));
				}
			pHit = &pTargHits[NumHiScores];
			}
		else
			pHit =  pTargHits;
 
		pHit->TargOfs = TIdx;
		pHit->QueryOfs = HiQIdx;
		pHit->HiScore = HiScore;
		NumHiScores += 1;
		}

	}
return(NumHiScores);
}