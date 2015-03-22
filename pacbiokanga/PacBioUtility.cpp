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
		m_LoadQuerySeqsRslt = (teBSFrsltCodes)Rslt;
		m_bAsyncLoading = false;
		ReleaseLock(true);
		return(Rslt);
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


int												// identified SMRTBell at this offset + 1, 0 if non-detected
CPacBioUtility::DetectSMRTBell(int StartOfs,	// search from this offset
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

if(StartOfs < 0 || pNumTetramers == NULL || InSeqLen < cSmartBellAdaptorSeqLen || pInSeq == NULL || MinTetramers < 1 ||
	StartOfs + cSmartBellAdaptorSeqLen > InSeqLen)
	return(eBSFerrParams);

PeakQKMers = 0;
PeakLoc = 0;
for(Idx = StartOfs+6; Idx < InSeqLen - cSmartBellAdaptorSeqLen/3 ; Idx++)
	{
	pQSeq = &cSmartBellAdaptorSeq[0];
	QKMers = 0;
	for(QIdx = 0; QIdx < cSmartBellAdaptorSeqLen - 6; QIdx++,pQSeq++)
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
					   etSeqBase *pInSeq,		// sequence which may contain homopolymers
					   etSeqBase *pOutSeq,		// homopolymer reduced sequence copied into this sequence buffer
					   int TruncLen)			// homopolymers longer than this length (min 1) to be truncated at this length
{
int InIdx;
int OutIdx;
int HomopolymerLen;

if(TruncLen < 1 || InSeqLen < TruncLen || pInSeq == NULL || *pInSeq > eBaseN || pOutSeq == NULL )
	return(-1);

OutIdx = 0;
HomopolymerLen = 0;
for(InIdx = 0; InIdx < (InSeqLen - 1); InIdx++,pInSeq++)
	{
	if(*pInSeq != pInSeq[1])
		{
		*pOutSeq++ = *pInSeq++;
		OutIdx += 1;
		HomopolymerLen = 0;
		}
	else
		{
		HomopolymerLen += 1;
		if(HomopolymerLen <= TruncLen)
			{
			*pOutSeq++ = *pInSeq++;
			OutIdx += 1;
			}
		}
	}
if(HomopolymerLen <= TruncLen || *pInSeq != pInSeq[-1])
	{
	*pOutSeq++ = *pInSeq;
	OutIdx += 1;
	}
return(OutIdx);
}

int													// number of putatve SMRTBell hairpins identified
CPacBioUtility::IdentifySMRTBells(int MaxSMRTBells,	// identify at most this many SMRTBells to return in pSMRTBells
						int SeqLen,			// length of sequence to search for SMRTBell hairpins
						etSeqBase *pSeq,	// identify all putative SMRTBell hairpins in this sequence
						tsSMRTBellHit *pSMRTBells, // returned identified SMRTBells
						int MinTetramers)	// only report SmartBell if at least this many tetramers in expected order detected
{
int NumTetramers;
int NxtLocSMRTBell;
int NumSMRTBells;

NxtLocSMRTBell = 0;

NumSMRTBells = 0;
while(NumSMRTBells < MaxSMRTBells && (NxtLocSMRTBell = DetectSMRTBell(NxtLocSMRTBell,&NumTetramers,SeqLen,pSeq,MinTetramers)) > 0)
	{
	pSMRTBells->NumTetramers = NumTetramers;
	pSMRTBells->LocOfs = NxtLocSMRTBell;
	NumSMRTBells += 1;
	NxtLocSMRTBell += cSmartBellAdaptorSeqLen;
	}
return(NumSMRTBells);
}

