/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */


#include "stdafx.h"

// Supporting at most this many concurrent TCP sessions between service requester (this server) and all service providers
// this could be increased to an internally restricted maximum of 511 but there could then be a significant performance throughput degradation
// because of the use of select() instead of ePoll or derivatives.
#ifdef WIN32
#define cMaxConcurrentSessions 100
#else
#define cMaxConcurrentSessions 100
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include <process.h>
#include "../libbiokanga/commhdrs.h"
// default on Windows is for a maximum of 64 sockets in a FD_SET select() descriptor set
// increase this limit to the maximum number of supported TCP sessions for all service provider sessions
// plus 1 for the listening socket and 1 for the control socket

#if (cMaxConcurrentSessions + 2 > 64)
#define FD_SETSIZE  (cMaxConcurrentSessions + 2)
#endif
#include <WinSock2.h>
#include <ws2tcpip.h>
#else
#include <sys/mman.h>
#include <pthread.h>
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <netdb.h>
typedef struct sockaddr_storage SOCKADDR_STORAGE;
#include "../libbiokanga/commhdrs.h"
#endif

#include "./BKScommon.h"
#include "./BKSRequester.h"
const char * cDfltListenerPort = "43123";		// default server port to listen on if not user specified

CBKSRequester::CBKSRequester()
{
m_szHostName[0] = '\0';
m_szServiceName[0] = '\0';
memset(&m_RequesterThread,0,sizeof(m_RequesterThread));
m_CreatedMutexes = 0;
m_ThreadActive = 0;
m_Terminate = 0;
m_LockEnabled = 0;
m_CASSerialise = 0;
m_bSessionTermReq = false;
m_NumInitTypes = 0;
m_NumSessions = 0;
m_NumSessEstabs = 0;	
m_CASSerialise = 0;
m_NumPendingReqs = 0;
m_NumPendingResps = 0;
m_TotRespsAvail = 0;
m_bNotifiedReqs = false;
m_bSessionTermReq = false;

m_pBKSSessEstabs = NULL;
m_pBKSTypes = NULL;

m_NumChkPtReqs = 0;
m_MaxChkPtReqs = 0;			
m_ppChkPtReqs = NULL;		

memset(m_SessionIDVect,0,sizeof(m_SessionIDVect));
memset(m_ReqIDVect, 0, sizeof(m_ReqIDVect));
size_t Size_m_Ctrl0 = sizeof(m_Ctrl[0]);
size_t Size_m_Ctrl1 = sizeof(m_Ctrl[1]);
size_t Size_m_Ctrl = sizeof(m_Ctrl);

memset(&m_Ctrl[0], 0, sizeof(m_Ctrl[0]));

memset(&m_Ctrl[1], 0, sizeof(m_Ctrl[1]));
#ifdef _WIN32
m_Ctrl[0].Socket = INVALID_SOCKET;
m_Ctrl[1].Socket = INVALID_SOCKET;
m_ListenerSock = INVALID_SOCKET;
#else
m_Ctrl[0].Socket = -1;
m_Ctrl[1].Socket = -1;
m_ListenerSock = -1;
#endif
}


CBKSRequester::~CBKSRequester()
{
Reset(false);
}

int
CBKSRequester::CreateMutexes(void)
{// check if already created...
#ifdef WIN32
if(InterlockedCompareExchange(&m_CreatedMutexes,1,1)==1)
#else
if(__sync_val_compare_and_swap(&m_CreatedMutexes,1,1)==1)
#endif
	return(eBSFSuccess);

#ifdef _WIN32
InitializeSRWLock(&m_hRWLock);
#else
if (pthread_rwlock_init(&m_hRWLock, NULL) != 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: unable to create and initialise m_hRWLock");
	return(eBSFerrInternal);
	}
#endif

m_CASSerialise = 0;
#ifdef WIN32
InterlockedExchange(&m_CreatedMutexes,1);
#else
__sync_val_compare_and_swap(&m_CreatedMutexes,0,1);
#endif
return(eBSFSuccess);
}


int	
CBKSRequester::DeleteMutexes(void)				// delete locks/mutexes used in access serialisations
{
#ifdef WIN32
if(InterlockedCompareExchange(&m_CreatedMutexes,0,0)==0)
#else
if(__sync_val_compare_and_swap(&m_CreatedMutexes,0,0)==0)
#endif
	return(eBSFSuccess);

#ifndef _WIN32
pthread_rwlock_destroy(&m_hRWLock);
#endif
#ifdef WIN32
InterlockedCompareExchange(&m_CreatedMutexes,0,1);
#else
__sync_val_compare_and_swap(&m_CreatedMutexes,1,0);
#endif
return(eBSFSuccess);
}

// acquire lock used for serialised access by multiple concurrent reader threads (bExclusive == false), or serialised access by a single thread (bExclusive == true)
bool			// returns true if lock obtained
CBKSRequester::AcquireLock(bool bExclusive)
{
#ifdef WIN32
if(InterlockedCompareExchange(&m_CreatedMutexes,1,1)==0)
#else
if(__sync_val_compare_and_swap(&m_CreatedMutexes,1,1)==0)
#endif
	return(false);

#ifdef _WIN32
if (bExclusive)
	AcquireSRWLockExclusive(&m_hRWLock);
else
	AcquireSRWLockShared(&m_hRWLock);
#else
if (bExclusive)
	pthread_rwlock_wrlock(&m_hRWLock);
else
	pthread_rwlock_rdlock(&m_hRWLock);
#endif

return(true);
}

// no longer needing acquired lock, let other threads gain access
bool		// returns true if lock released
CBKSRequester::ReleaseLock(bool bExclusive)
{
#ifdef WIN32
if(InterlockedCompareExchange(&m_CreatedMutexes,1,1)!=1)
#else
if(__sync_val_compare_and_swap(&m_CreatedMutexes,1,1)!=1)
#endif
	return(false);

#ifdef _WIN32
if (bExclusive)
	ReleaseSRWLockExclusive(&m_hRWLock);
else
	ReleaseSRWLockShared(&m_hRWLock);
#else
pthread_rwlock_unlock(&m_hRWLock);
#endif
return(true);
}

void
CBKSRequester::AcquireCASSerialise(bool bPriority)	// if bPriority true then backoff time is reduced relative to if bPriority is false, increasing the probability of acquiring the serialisation lock if there is contention
{
int SpinCnt = 10;
int BackoffMS = 1;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_CASSerialise,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 10;
	if(!bPriority)
		{
		if(BackoffMS < 50)
			BackoffMS += 1;
		else
			BackoffMS = 20 + (rand() % 31);
		}
	}
#else
while(__sync_val_compare_and_swap(&m_CASSerialise,0,1)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 10;
	if(!bPriority)
		{
		if(BackoffMS < 50)
			BackoffMS += 1;
		else
			BackoffMS = 20 + (rand() % 31);
		}
	}
#endif

}

void
CBKSRequester::ReleaseCASSerialise(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASSerialise,0,1);
#else
__sync_val_compare_and_swap(&m_CASSerialise,1,0);
#endif
}


#ifdef _WIN32
unsigned __stdcall RequesterThread(void * pThreadPars)
#else
void *RequesterThread(void * pThreadPars)
#endif
{
int Rslt;
tsRequesterThread *pPars = (tsRequesterThread *)pThreadPars;			// makes it easier not having to deal with casts!
CBKSRequester *pRequester = (CBKSRequester *)pPars->pThis;

Rslt = pRequester->StartThread(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

bool				// false if unable to begin thread execution or if executing thread took too long (default 3min) to initialise ready to accept connections
CBKSRequester::Run(int Secs)	// expecting thread to take at most this many seconds to start and initialise ready to accept connections
{
UINT32 MaxWait;
UINT32 ThreadActive;
UINT32 Started;

tsRequesterThread *pThreadPar;
m_ThreadActive = 0;
m_Terminate = 0;
ThreadActive = 1;
pThreadPar = &m_RequesterThread;
memset(pThreadPar,0,sizeof(tsRequesterThread));
pThreadPar->ThreadIdx = 1;
pThreadPar->pThis = this;
#ifdef _WIN32
pThreadPar->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, RequesterThread, pThreadPar, 0, &pThreadPar->threadID);
#else
pThreadPar->threadRslt = pthread_create(&pThreadPar->threadID, NULL, RequesterThread, pThreadPar);
#endif
if(Secs < 20)		// floor on Secs is 20 seconds
	Secs = 20;
// allow thread Secs to startup
MaxWait = Secs;		// allowing at most this many seconds for thread to startup, will check every second
do {
#ifdef WIN32
	if((Started = InterlockedCompareExchange(&m_ThreadActive,ThreadActive,ThreadActive)) != ThreadActive)
		Sleep(1000);
#else

	if((Started = __sync_val_compare_and_swap (&m_ThreadActive,ThreadActive,ThreadActive)) != ThreadActive)
		sleep(1);
#endif
	MaxWait -= 1;
	}
while(Started != ThreadActive && MaxWait > 0);
if(Started != ThreadActive)
	{
	// after Secs give up waiting and ask thread to self terminate, will allow at most 120 seconds and if thread does not self terminate then force terminate
#ifdef WIN32
	InterlockedCompareExchange(&m_Terminate,1,0);
#else
	__sync_val_compare_and_swap (&m_Terminate,0,1);
#endif

#ifdef WIN32
	if(WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, 120*1000))
		TerminateThread(pThreadPar->threadHandle,0);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 120;
	if ((JoinRlt = pthread_timedjoin_np(pThreadPar->threadID, NULL, &ts)) != 0)
		{
		pthread_cancel(pThreadPar->threadID);	
		pthread_join(pThreadPar->threadID, NULL);
		}
#endif
	Started=0;
	m_ThreadActive = 0;
	m_Terminate = 0;
	}
return(Started == 0 ? false : true);
}

void
CBKSRequester::Terminate(int Secs)			// allow this many seconds for thread to self-terminate before forcing termination
{
#ifdef WIN32
if(InterlockedCompareExchange(&m_Terminate,1,1) == 1)	// already requested to self-terminate?
	return;
#else
if(__sync_val_compare_and_swap (&m_Terminate,1,1) == 1)	// already requested to self-terminate?
	return;
#endif
	// ask thread to self terminate, will allow at most Secs seconds and if thread does not self terminate then force terminate
#ifdef WIN32
	InterlockedCompareExchange(&m_Terminate,1,0);
#else
	__sync_val_compare_and_swap (&m_Terminate,0,1);
#endif

#ifdef WIN32
	if(WAIT_TIMEOUT == WaitForSingleObject(m_RequesterThread.threadHandle, Secs * 1000))
		TerminateThread(m_RequesterThread.threadHandle,0);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += Secs;
	if ((JoinRlt = pthread_timedjoin_np(m_RequesterThread.threadID, NULL, &ts)) != 0)
		{
		pthread_cancel(m_RequesterThread.threadID);	
		pthread_join(m_RequesterThread.threadID, NULL);
		}
#endif
	m_ThreadActive = 0;
	m_Terminate = 0;

}

int
CBKSRequester::StartThread(tsRequesterThread *pPars)
{
#ifdef WIN32
InterlockedIncrement(&m_ThreadActive);
#else
__sync_fetch_and_add(&m_ThreadActive,1);
#endif
return(AcceptConnections());
}

const UINT32 cMinJobWaitSecs = 15;				// job wait times (secs) will be clamped to be within
const UINT32 cMaxJobWaitSecs = (60*60*24);		// this range - 1 day should be long enough :-)
const UINT32 cMaxParamsLen = 8192;				// job parameters can be at most this many bytes long

UINT32											// returned request identifier or 0 if all identifiers have already been allocated
CBKSRequester::AllocReqID(void)					// returns next available unused request identifier and sets that identifier as allocated
{
int ReqID;
UINT8 *pByte;
UINT8 Msk;
UINT32 *pVect;

pVect = m_ReqIDVect;
for(ReqID = 1; ReqID <= cMaxReqID; ReqID+=32,pVect+=1 )
	{
	if(*pVect == 0x0ffffffff)
		continue;
	pByte = (UINT8 *)pVect;
	while(*pByte == 0x0ff)
		{
		ReqID += 8;
		pByte += 1;
		}
	Msk = 0x01;
	while(*pByte & Msk)
		{
		ReqID += 1;
		Msk <<= 1;
		}
	*pByte |= Msk;
	return(ReqID);
	}
return(0);
}

bool				// true if ReqID was previously allocated (used), false if previously unallocated (unused)
CBKSRequester::UnallocReqID(UINT32 ReqID)   // sets the ReqID in m_ReqIDVect as unallocated and available to be reused
{
bool bAllocd;
UINT32 *pVect;
UINT32 Msk;
if(ReqID == 0 || ReqID > cMaxReqID)
	return(false);
pVect = &m_ReqIDVect[(ReqID-1)/32];
Msk = 0x01 << (ReqID-1) % 32;
bAllocd = *pVect & Msk ? true : false;
*pVect &= ~Msk;
return(bAllocd);
}

bool				// true if ReqID is currently allocated (used), false if unallocated (unused)
CBKSRequester::IsAllocReqID(UINT32 ReqID)   // checking this ReqID in m_ReqIDVect
{
UINT32 *pVect;
UINT32 Msk;
if (ReqID == 0 || ReqID > cMaxReqID)
	return(false);
pVect = &m_ReqIDVect[(ReqID - 1) / 32];
Msk = 0x01 << (ReqID - 1) % 32;
return(*pVect & Msk ? true : false);
}


UINT32											// returned session identifier or 0 if all identifiers have already been allocated
CBKSRequester::AllocSessionID(void)					// returns next available unused session identifier and sets that identifier in m_SessionIDVect as allocated
{
int SessionID;
UINT8 *pByte;
UINT8 Msk;
UINT32 *pVect;

pVect = m_SessionIDVect;
for (SessionID = 1; SessionID <= cMaxConcurrentSessions; SessionID += 32, pVect += 1)
	{
	if (*pVect == 0x0ffffffff)
		continue;
	pByte = (UINT8 *)pVect;
	while (*pByte == 0x0ff)
		{
		SessionID += 8;
		pByte += 1;
		}
	Msk = 0x01;
	while (*pByte & Msk)
		{
		SessionID += 1;
		Msk <<= 1;
		}
	*pByte |= Msk;
	return(SessionID);
	}
return(0);
}

bool				// true if SessionID was previously allocated (used), false if previously unallocated (unused)
CBKSRequester::UnallocSessionID(UINT32 SessionID)   // sets the SessionID in m_SessionIDVect as unallocated and available to be reused
{
bool bAllocd;
UINT32 *pVect;
UINT32 Msk;
if (SessionID == 0 || SessionID > cMaxConcurrentSessions)
	return(false);
pVect = &m_SessionIDVect[(SessionID - 1) / 32];
Msk = 0x01 << (SessionID - 1) % 32;
bAllocd = *pVect & Msk ? true : false;
*pVect &= ~Msk;
return(bAllocd);
}

bool				// true if SessionID is currently allocated (used), false if unallocated (unused)
CBKSRequester::IsAllocSessionID(UINT32 SessionID)   // checking this SessionID in m_SessionIDVect
{
UINT32 *pVect;
UINT32 Msk;
if (SessionID == 0 || SessionID > cMaxConcurrentSessions)
	return(false);
pVect = &m_SessionIDVect[(SessionID - 1) / 32];
Msk = 0x01 << (SessionID - 1) % 32;
return(*pVect & Msk ? true : false);
}


tJobIDEx										// packed job identifier or 0 if range errors
CBKSRequester::PackIntoJobIDEx(UINT32 ReqID,	// must be in the range 1..16777215 (24bits)
				UINT32 SessionID,				// service provider session identifier, must be in the range 1..511 (9bits)
				UINT32 InstanceID,				// service instance in the service providers session SessionID, must be in the range 1..511 (9bits)
				UINT32 TypeID,					// identifies service being provided by the service provider, must be in the range 1..15 (4bits)
				UINT32 TypeSessionID)			// index of this session in the service type, must be in the range 1..511 (9bits)
{
tJobIDEx JobIDEx;

if(TypeSessionID < 1 || TypeSessionID > 511)
	return(0);
JobIDEx = TypeSessionID;
JobIDEx <<= 4;

if (TypeID < 1 || TypeID > 15)
	return(0);
JobIDEx |= TypeID;
JobIDEx <<= 9;

if (InstanceID < 1 || InstanceID > 511)
	return(0);
JobIDEx |= InstanceID;
JobIDEx <<= 17;

if (SessionID < 1 || SessionID > 511)
	return(0);
JobIDEx |= SessionID;
JobIDEx <<= 24;

if (ReqID < 1 || ReqID > 16777215)
return(0);
JobIDEx |= ReqID;

return(JobIDEx);
}

bool				// false if any range errors whilst unpacking
CBKSRequester::UnpackFromJobIDEx(tJobIDEx JobIDEx,				// unpack from this extended job identifier
				UINT32 *pReqID,				    // returned ReqID, will be in the range 1..16777215
				UINT32 *pSessionID,				// returned SessionID, will be in the range 1..131071
				UINT32 *pInstanceID,			// returned InstanceID, will be in the range 1..511
				UINT32 *pTypeID,				// returned service TypeID, will be in the range 1..15
				UINT32 *pTypeSessionID)			// returned index of this session in the service type, will be in the range 1..511
{
tJobIDEx OrigJobIDEx;
UINT32 ReqID;
UINT32 SessionID;
UINT32 InstanceID;
UINT32 TypeID;
UINT32 TypeSessionID;
tsBKSType *pType;
tsBKSRegSessionEx *pSession;
UINT32 ReqRespInstOfs;
tsReqRespInst *pReqRespInst;

if(pReqID != NULL)
	*pReqID = 0;
if(pSessionID != NULL)
	*pSessionID = 0;
if(pInstanceID != NULL)
	*pInstanceID = 0;
if(pTypeID != NULL)
	*pTypeID = 0;
if(pTypeSessionID != NULL)
	*pTypeSessionID = 0;

if((OrigJobIDEx = JobIDEx) <= 0)
	return(false);

ReqID = JobIDEx & 0x0ffffff;
if(ReqID == 0)
	return(false);
if(pReqID != NULL)
	*pReqID = ReqID;
JobIDEx >>= 24;

SessionID = JobIDEx & 0x01ffff;
if (SessionID == 0)
	return(false);
if (pSessionID != NULL)
	*pSessionID = SessionID;
JobIDEx >>= 17;

InstanceID = JobIDEx & 0x01ff;
if (InstanceID == 0)
	return(false);
if (pInstanceID != NULL)
	*pInstanceID = InstanceID;
JobIDEx >>= 9;

TypeID = JobIDEx & 0x0f;
if(TypeID == eBKSPTUndefined || TypeID >= eBKSPTPlaceHolder)
	return(false);
if (pTypeID != NULL)
	*pTypeID = TypeID;
JobIDEx >>= 4;

TypeSessionID = JobIDEx & 0x01ff;
if (TypeSessionID == 0)
	return(false);
if (pTypeSessionID != NULL)
	*pTypeSessionID = TypeSessionID;

pType = &m_pBKSTypes[TypeID - 1];
if(pType->Detail.BKSPType != TypeID || TypeSessionID > pType->MaxSessions)
	return(false);

pSession = pType->pSessions[TypeSessionID-1];
if(pSession == NULL || pSession->Session.SessionID != SessionID || pSession->Session.BKSPState >= eBKSPSRegisteredTerm)
	return(false);

if(pSession->Session.MaxInstances < InstanceID)
	return(false);

ReqRespInstOfs = pType->ReqRespInstSize;
ReqRespInstOfs *= (InstanceID-1);
pReqRespInst = (tsReqRespInst *)&pSession->pReqResp[ReqRespInstOfs];
if(pReqRespInst->JobIDEx != OrigJobIDEx)
	return(false);
return(true);
}

int												// total number of classes 
CBKSRequester::GetNumClassInstances(teBKSPType TypeID,		// service required
						UINT32 *pCommited,			// returned number of class instances currently committed or instantiated 
						UINT32 *pUncommited)    // returned number of class instances currently not committed and available to be instantiated
{
tsBKSType *pType;
tsBKSRegSessionEx *pSession;
UINT32 NumCommited;
UINT32 NumUncommited;
if(pCommited != NULL)
	*pCommited = 0;
if(pUncommited != NULL)
	*pUncommited = 0;
NumCommited = 0;
NumUncommited = 0;
if(TypeID <= eBKSPTUndefined || TypeID >= eBKSPTPlaceHolder)
	return(-1);
pType= &m_pBKSTypes[TypeID -1];

AcquireLock(true);

if(pType->NumSessions == 0)			// need at least one service provider to continue exploring if job can be accepted
	{
	ReleaseLock(true);
	return(0);
	}
pSession = pType->pFirstSession;
do {
	if(pSession->Session.BKSPState != eBKSPSRegisteredActv || pSession->Session.MaxClassInstances == 0)
		continue;
	NumCommited += pSession->NumClassInstances;
	NumUncommited += pSession->Session.MaxClassInstances - pSession->NumClassInstances;
	}
while((pSession = pSession->pNext) != NULL);
ReleaseLock(true);

if(pCommited != NULL)
	*pCommited = NumCommited;
if(pUncommited != NULL)
	*pUncommited = NumUncommited;
return(NumCommited + NumUncommited);
}


// NOTE: processing is to attempt to distribute jobs evenly to all service providers
int				// -2: parameter errors, -1: class instance no longer exists, 0: currently no available service instance 1: if job accepted
CBKSRequester::AddJobRequest(tJobIDEx *pJobID,	// returned unique job identifier by which job can later be referenced
						teBKSPType TypeID,		// service required
						UINT64 ClassInstanceID,	// class instance on which job method is to be applied - can be 0 if no pre-existing class
						UINT32 ClassMethodID,	// class method to apply on the class instance
						UINT32 ParamsSize,		// processing parameters are this total size in bytes
						void *pParams,			// service processing parameters
						UINT32 InDataSize,		// service processing input data is this total size in bytes
						void *pInData,			// service processing input data
						UINT32 SessionID)		// if 0 then use session as specified by ClassInstanceID, otherwise use session corresponding to specific session identifier
{
tsBKSType *pType;
UINT32 ReqRespInstIdx;
tsReqRespInst *pReqRespInst;
tsBKSRegSessionEx *pSession;
tsBKSRegSessionEx *pLeastBusy;
double LeastBusy;
double CurBusy;
UINT32 Idx;
UINT32 ReqID;
UINT32 JobSessionID;
UINT32 InstanceID;
UINT64 *pClassIdentifier;
UINT32 TypeSessionID;

// validate parameters
if(pJobID == NULL || TypeID <= eBKSPTUndefined || TypeID >= eBKSPTPlaceHolder)
	return(-2);
*pJobID = 0;
pType= &m_pBKSTypes[TypeID -1];
if(pType->Detail.BKSPType != TypeID ||
   ParamsSize > pType->Detail.MaxParamLen || (ParamsSize > 0 && pParams == NULL) ||
   InDataSize >  pType->Detail.MaxQuerySeqLen || (InDataSize > 0 && pInData == NULL))
	return(-2);

#ifdef WIN32
	InterlockedIncrement(&m_NumPendingReqs);	  // letting TCP session handling thread there is at least one pending request to be processed	
#else
	__sync_fetch_and_add(&m_NumPendingReqs,1);
#endif

AcquireLock(true);

#ifdef WIN32
	InterlockedDecrement(&m_NumPendingReqs);	  // letting TCP session handling thread know this request no longer outstanding	
#else
	__sync_fetch_and_sub(&m_NumPendingReqs,1);
#endif	

if(pType->NumSessions == 0)			// need at least one service provider to continue exploring if job can be accepted
	{
	ReleaseLock(true);
	return(-1);
	}

pLeastBusy = NULL;
if(ClassInstanceID == 0)		// no pre-existing class instance if 0
	{
	// find and note provider with highest capacity to process this job request as a proportion of that providers capacity
	pSession = pType->pFirstSession;
	LeastBusy = 0.0;
	do {
		if(pSession->Session.BKSPState != eBKSPSRegisteredActv ||
			pSession->Session.NumBusy == pSession->Session.MaxInstances)
			continue;

		if(SessionID != 0 && pSession->Session.SessionID != SessionID)
			continue;
		
		if(pSession->NumClassInstances > 0 && pSession->NumClassInstances >= pSession->Session.MaxClassInstances)
			continue;

		if(pLeastBusy == NULL || ((CurBusy = (pSession->Session.NumBusy/pSession->Session.MaxInstances)) <= LeastBusy))
			{
			if(pLeastBusy == NULL || CurBusy < LeastBusy || pSession->Session.MaxInstances > pLeastBusy->Session.MaxInstances)
				{
				LeastBusy = pSession->Session.NumBusy/pSession->Session.MaxInstances;
				pLeastBusy = pSession;
				}
			}
		}
	while(pSession->Session.SessionID != SessionID && (pSession = pSession->pNext) != NULL);
	if(pLeastBusy == NULL)
		{
		ReleaseLock(true);
		return(0);		// currently no service provider capacity to accept request
		}
	pSession = pLeastBusy;
	}
else   // method is on an already instantiated class 
	{
	JobSessionID = (UINT32)((UINT64)ClassInstanceID >> 53);
	if(SessionID != 0 && SessionID != JobSessionID)
		{
		ReleaseLock(true);
		return(-1);				// class no longer exists - original session may have been terminated
		}
	pSession = pType->pFirstSession;
	do {
		if(pSession->Session.SessionID != JobSessionID)
			continue;
		pClassIdentifier = pSession->ClassInstanceIDs;
		for(Idx = 0; Idx < pSession->NumClassInstances; Idx++, pClassIdentifier++)
			{
			if(ClassInstanceID == *pClassIdentifier)
				{
				pLeastBusy = pSession;
				break;
				}
			}
		}
	while(pLeastBusy == NULL && (pSession = pSession->pNext) != NULL);
	if(pLeastBusy == NULL || pLeastBusy->Session.NumBusy >= pLeastBusy->Session.MaxInstances)
		{
		ReleaseLock(true);
		return(-1);				// class no longer exists - original session may have been terminated
		}
	pSession = pLeastBusy;
	}

ReqRespInstIdx = 0;
pReqRespInst = (tsReqRespInst *)pSession->pReqResp;
for(InstanceID = 1; InstanceID <= pSession->Session.MaxInstances; InstanceID++)
	{
	if(pReqRespInst->ReqID == 0)
		{	
		memset(pReqRespInst,0,pType->ReqRespInstSize);
		ReqID = AllocReqID();
		JobSessionID = pSession->Session.SessionID;
		TypeSessionID = pSession->Session.TypeSessionID;
		pReqRespInst->JobIDEx = PackIntoJobIDEx(ReqID, JobSessionID, InstanceID, TypeID, TypeSessionID);
		pReqRespInst->ReqID = ReqID;
		pReqRespInst->ClassInstanceID = ClassInstanceID;
		pReqRespInst->ClassMethodID = ClassMethodID;
		pReqRespInst->ParamSize = ParamsSize;
		pReqRespInst->InDataSize = InDataSize;
		if(ParamsSize > 0)
			memcpy(pReqRespInst->Data, pParams, ParamsSize);
		if(InDataSize > 0)
			memcpy(&pReqRespInst->Data[ParamsSize], pInData, InDataSize);
		pReqRespInst->FlgReq = 1;
		pReqRespInst->SubmitAt = (UINT32)time(NULL);
		pSession->Session.NumReqs += 1;
		pSession->Session.NumBusy += 1;
		pType->FlgReq = 1;
		if(!m_bNotifiedReqs)
			{
			m_bNotifiedReqs = true;
			NotifyCtrl();
			}
		*pJobID = pReqRespInst->JobIDEx;
		ReleaseLock(true);
		return(1);
		}
	pReqRespInst = (tsReqRespInst *)((UINT8 *)pReqRespInst + pType->ReqRespInstSize);
	}
ReleaseLock(true);
return(0);		// currently no service provider capacity to accept request
}

int
CBKSRequester::SendRequestFrames(void)			// iterate all sessions and if any frames ready to send and room to accept the frame in TxdBuff then initiate the sending
{
int Idx;
UINT32 TxFrameID;
UINT32 FrameLen;
tsBKSType *pType;
UINT32 ReqRespInstIdx;
tsReqRespInst *pReqRespInst;
tsBKSRegSessionEx *pSession;
tsTxdRxd *pTxdRxd;
sBKSServReq *pFrame;

pType = m_pBKSTypes;
for (Idx = 0; Idx < eBKSPTPlaceHolder-1; Idx++, pType += 1)
	{
	if (pType->Detail.BKSPType == eBKSPTUndefined || pType->NumSessions == 0)
		continue;
	pSession = pType->pFirstSession;
	do {
		if(pSession->Session.BKSPState != eBKSPSRegisteredActv || pSession->Session.NumReqs == 0)
			continue;

		pTxdRxd = &pSession->TxdRxd;
		
		pReqRespInst = (tsReqRespInst *)&pSession->pReqResp[(pSession->LastChkdReqIdx * pType->ReqRespInstSize)];
		for(ReqRespInstIdx = pSession->LastChkdReqIdx; ReqRespInstIdx < pSession->Session.MaxInstances; ReqRespInstIdx++, pReqRespInst = (tsReqRespInst *)((UINT8 *)pReqRespInst + pType->ReqRespInstSize))
			{
			if(pSession->Session.NumReqs == 0)
				{
				pSession->LastChkdReqIdx = 0;     // no outstanding requests so can start checks from the 1st instance next time this session is processed for request frames to be sent
				break;
				}
			TxFrameID = pTxdRxd->TxFrameID;						// throttle back on adding any new frames until service provider has caught up a little
			if(pTxdRxd->RxdRxFrameID > pTxdRxd->TxFrameID)
				TxFrameID += 0x07f;
			if((TxFrameID - pTxdRxd->RxdRxFrameID) > 4)
				{
				pSession->LastChkdReqIdx = ReqRespInstIdx;	// still outstanding requests but needing to throttle back, start checks from this instance next time this session is processed for request frames to be sent
				break;
				}

			if(pReqRespInst->FlgReq)		// requested to be sent?
				{
				FrameLen = sizeof(sBKSServReq) - 1 + pReqRespInst->ParamSize + pReqRespInst->InDataSize;
				if((pTxdRxd->AllocdTxdBuff - pSession->TxdRxd.TotTxd) > (FrameLen + (sizeof(tsBKSPacHdr) * 5))) // always allow spare room for some session control frames
					{
					pFrame = (sBKSServReq *)&pTxdRxd->pTxdBuff[pSession->TxdRxd.TotTxd];
					pFrame->Hdr.FrameFlags = 0;
					pFrame->Hdr.FrameLen = FrameLen;
					pFrame->Hdr.FrameType = eBKSHdrReq;
					pFrame->Hdr.RxFrameID = pTxdRxd->RxdTxFrameID;
					pFrame->Hdr.TxFrameID = pTxdRxd->TxFrameID++;
					if(pTxdRxd->TxFrameID > 0x07f)
						pTxdRxd->TxFrameID = 1;
					pFrame->Hdr.SessionID = pTxdRxd->SessionID;
					pFrame->JobIDEx = pReqRespInst->JobIDEx;
					pFrame->ClassInstanceID = pReqRespInst->ClassInstanceID;
					pFrame->ClassMethodID = pReqRespInst->ClassMethodID;
					pFrame->ParamSize = pReqRespInst->ParamSize;
					pFrame->DataSize = pReqRespInst->InDataSize;
					if(pReqRespInst->ParamSize > 0 || pReqRespInst->InDataSize > 0)
						memcpy(pFrame->ParamData,pReqRespInst->Data,pReqRespInst->ParamSize + pReqRespInst->InDataSize);
					pReqRespInst->FlgReq = 0;
					pReqRespInst->FlgProc = 1;
					pSession->Session.NumReqs -= 1;
					pSession->Session.NumProcs += 1;
					pSession->TxdRxd.TotTxd += FrameLen;
					pTxdRxd->flgSelMonWrite = 1;
					if(pTxdRxd->CurTxd == 0)
						{
						if(!TxData(pTxdRxd))
							{
							ShutdownConnection(&pSession->TxdRxd.Socket);
							pSession->Session.BKSPState = eBKSPSRegisteredTerm;
							}
						}
					}
				else                  // no room in the txd buffer so defer sending until txd buffer has been, at least partially, emptied
					{
					pSession->LastChkdReqIdx = ReqRespInstIdx;
					break;
					}
				}
			}
		if(ReqRespInstIdx == pSession->Session.MaxInstances)
			pSession->LastChkdReqIdx = 0;
		}
	while((pSession = pSession->pNext) != NULL);
	}
return(0);
}

int   // 0: accepted frame, -1: JobIDEx errors, -2 Session or type errors, -3 mismatch between instance JobIDEx's, ClassInstanceID mismatch
CBKSRequester::ProcessResponseFrame(tsBKSRegSessionEx *pSession)	// process a received response frame
{
UINT32 Diff;
UINT32 ReqID;
UINT32 SessionID;
UINT32 InstanceID;
UINT32 TypeID;
UINT32 TypeSessionID;
UINT32 InstanceOfs;
tsBKSType *pType;
tsTxdRxd *pTxdRxd;
tsReqRespInst *pInstance;
sBKSServResp *pResponse;
UINT32 Idx;
UINT64 *pClassIdentifier;

pTxdRxd = &pSession->TxdRxd;
pResponse = (sBKSServResp *)pTxdRxd->pRxdBuff;
if(!UnpackFromJobIDEx(pResponse->JobIDEx,&ReqID,&SessionID,&InstanceID,&TypeID,&TypeSessionID))
	return(-1);

if(SessionID != pTxdRxd->SessionID || TypeID != pSession->Session.BKSPType || TypeSessionID != pSession->Session.TypeSessionID)
	return(-2);

pType = &m_pBKSTypes[TypeID-1];
InstanceOfs = pType->ReqRespInstSize * (InstanceID - 1);
pInstance = (tsReqRespInst *)&pSession->pReqResp[InstanceOfs];
if(pInstance->JobIDEx != pResponse->JobIDEx || !pInstance->FlgProc)
	return(-3);

pInstance->CpltdAt = (UINT32)time(NULL);
pInstance->JobRslt = pResponse->JobRslt;
if(pResponse->ClassInstanceID != 0)
	{
	if(pInstance->ClassInstanceID != 0 && pInstance->ClassInstanceID != pResponse->ClassInstanceID)
		return(-4);

	pClassIdentifier = pSession->ClassInstanceIDs;
	for(Idx = 0; Idx < pSession->NumClassInstances; Idx++, pClassIdentifier++)
		{
		if(pResponse->ClassInstanceID == *pClassIdentifier)
			break;
		}

	if(Idx == pSession->NumClassInstances)  // new class instance?
		{
		pSession->ClassInstanceIDs[Idx] = pResponse->ClassInstanceID;
		pSession->NumClassInstances += 1;
		}
	pInstance->ClassInstanceID = pResponse->ClassInstanceID;
	}
else // class instance identifier was 0, if request specified a class identifier then that class instance is no longer valid 
	{
	if(pInstance->ClassInstanceID != 0 && pSession->NumClassInstances > 0)
		{
		pClassIdentifier = pSession->ClassInstanceIDs;
		for(Idx = 0; Idx < pSession->NumClassInstances; Idx++, pClassIdentifier++)
			{
			if(pInstance->ClassInstanceID == *pClassIdentifier)
				break;
			}
		if(Idx != pSession->NumClassInstances)
			{
			if(Idx < pSession->NumClassInstances-1)
				memmove(pClassIdentifier,&pClassIdentifier[1],(pSession->NumClassInstances-Idx - 1)*sizeof(UINT64));
			pSession->NumClassInstances -= 1;
			}
		}
	pInstance->ClassInstanceID = 0;
	}
		
pInstance->ClassMethodID = pResponse->ClassMethodID;
pInstance->OutDataSize = pResponse->DataSize;
if(pResponse->DataSize > 0)
	memcpy(pInstance->Data,pResponse->Data,pResponse->DataSize);
pInstance->FlgCpltd = 1;
pSession->Session.NumProcs -= 1;
pSession->Session.NumCpltd += 1;

#ifdef WIN32
InterlockedIncrement(&m_TotRespsAvail);	  	
#else
__sync_fetch_and_add(&m_TotRespsAvail,1);
#endif

pSession->Session.TotNumCpltd += 1;
pTxdRxd->flgKeepAliveReq = 1;

if((Diff = (pTxdRxd->TotRxd - pTxdRxd->CurPacRxd)) > 0)
	memmove(pTxdRxd->pRxdBuff,&pTxdRxd->pRxdBuff[pTxdRxd->CurPacRxd], Diff);
pTxdRxd->TotRxd = Diff;
pTxdRxd->CurPacRxd = 0;
pTxdRxd->flgRxCplt = 0;
return(0);
}

int				// < 0 if job no longer exists, 0 if job still being processed, > 0 if job completed
CBKSRequester::GetJobResponse(tJobIDEx	JobID,	// unique job identifier returned when job was originally submitted
				  UINT64 *pClassInstanceID,		// returned class instance on which job method was applied
				  UINT32 *pClassMethodID,		// returned class method applied
				   UINT32 *pJobRslt,		// job processing result as returned by service provider
				   UINT32 *pOutDataSize,	// (IN) service processing output results expected to be at most this total length, [OUT] bytes of response data copied into pOutData 
				   void *pOutData,			// service processing output results data
				   bool bRetain)		    // true if job response is to be retained and not deleted; subsequent call with bRetain==false will delete this response
{
tsBKSType *pType;
UINT32 TypeID;
UINT32 ReqID;
UINT32 SessionID;
UINT32 InstanceID;
UINT32 TypeSessionID;
UINT32 TotRespsAvail;
tsBKSRegSessionEx *pSession;
tsReqRespInst *pReqRespInst;
UINT32 ReqRespInstOfs;
UINT32 CpySize;

if(JobID < 1)					// job identifier must be supplied and it must be valid
	return(-1);

if(pJobRslt != NULL)
	*pJobRslt = 0;

AcquireLock(false);
if(!UnpackFromJobIDEx(JobID,&ReqID,&SessionID,&InstanceID,&TypeID,&TypeSessionID))
	{
	ReleaseLock(false);
	return(-1);
	}
ReleaseLock(false);

if(TypeID == eBKSPTUndefined || TypeID >= eBKSPTPlaceHolder)
	return(-1);

#ifdef WIN32
TotRespsAvail = InterlockedCompareExchange(&m_TotRespsAvail,0,0);
#else
TotRespsAvail = __sync_val_compare_and_swap(&m_TotRespsAvail,0,0);
#endif
if(TotRespsAvail == 0)
	return(0);

#ifdef WIN32
	InterlockedIncrement(&m_NumPendingResps);	  // letting TCP session handling thread there is at least one worker thread wanting to check for a job response	
#else
	__sync_fetch_and_add(&m_NumPendingResps,1);
#endif

AcquireLock(true);
#ifdef WIN32
	InterlockedDecrement(&m_NumPendingResps);	  // letting TCP session handling thread one less worker thread wanting to check for a job response	
#else
	__sync_fetch_and_sub(&m_NumPendingResps,1);
#endif
pType = &m_pBKSTypes[TypeID - 1];
if(pType->Detail.BKSPType != TypeID || TypeSessionID == 0 || TypeSessionID > pType->MaxSessions)
	{
	ReleaseLock(true);
	return(-1);
	}

pSession = pType->pSessions[TypeSessionID-1];
if(pSession == NULL || pSession->Session.SessionID != SessionID || pSession->Session.BKSPState >= eBKSPSRegisteredTerm)
	{
	ReleaseLock(true);
	return(-1);
	}

if(pSession->Session.MaxInstances < InstanceID)
	{
	ReleaseLock(true);
	return(-1);
	}

ReqRespInstOfs = pType->ReqRespInstSize;
ReqRespInstOfs *= (InstanceID-1);
pReqRespInst = (tsReqRespInst *)&pSession->pReqResp[ReqRespInstOfs];

if(pReqRespInst->FlgCpltd)		// set if job has completed
	{
	*pClassInstanceID = pReqRespInst->ClassInstanceID;
	*pClassMethodID = pReqRespInst->ClassMethodID;
	if(pOutData != NULL && pOutDataSize != NULL && *pOutDataSize > 0)
		CpySize = min(pReqRespInst->OutDataSize, *pOutDataSize);
	else
		CpySize = 0;
	if(CpySize)
		memcpy(pOutData,pReqRespInst->Data, CpySize);
	if(pOutDataSize != NULL)
		*pOutDataSize = CpySize;
	if(pJobRslt != NULL)
		*pJobRslt = pReqRespInst->JobRslt;
	if(!bRetain)
		{
		pReqRespInst->FlgCpltd = 0;
		pReqRespInst->ClassInstanceID = 0;
		pReqRespInst->ClassMethodID = 0;
		pSession->Session.NumCpltd -= 1;

#ifdef WIN32
		InterlockedDecrement(&m_TotRespsAvail);	  // letting TCP session handling thread know there is at least one worker thread wanting to check for a job response	
#else
		__sync_fetch_and_sub(&m_TotRespsAvail,1);
#endif
		pSession->Session.NumBusy -= 1;
		memset(pReqRespInst,0,pType->ReqRespInstSize);
		UnallocReqID(ReqID);
		}
	ReleaseLock(true);
	return(1);
	}

// no response available
ReleaseLock(true);
return(0);
}


int
CBKSRequester::TerminateAllSessions(void)
{
UINT32 Idx;

tsTxdRxd *pTxdRxd;
tsBKSType *pType;
tsBKSRegSessionEx *pSession;
tsBKSRegSessionEx *pNext;
tsBKSSessEstab *pSessEstab;

m_bSessionTermReq = true;		// flags that all current Session connections are being terminated

if(m_pBKSSessEstabs != NULL)
	{
	pSessEstab = m_pBKSSessEstabs;
	for (Idx = 0; Idx < cMaxSessEstab; Idx++, pSessEstab++)
		{
	#ifdef _WIN32
		if(pSessEstab->TxdRxd.Socket != INVALID_SOCKET)
			closesocket(pSessEstab->TxdRxd.Socket);
		if(pSessEstab->TxdRxd.pRxdBuff!=NULL)
			free(pSessEstab->TxdRxd.pRxdBuff);
		if (pSessEstab->TxdRxd.pTxdBuff != NULL)
			free(pSessEstab->TxdRxd.pTxdBuff);
		memset(pSessEstab,0,sizeof(tsBKSSessEstab));
		pSessEstab->TxdRxd.Socket = INVALID_SOCKET;
	#else
		if (pSessEstab->TxdRxd.Socket != -1)
			close(pSessEstab->TxdRxd.Socket);
		if (pSessEstab->TxdRxd.pRxdBuff != NULL)
			free(pSessEstab->TxdRxd.pRxdBuff);
		if (pSessEstab->TxdRxd.pTxdBuff != NULL)
			free(pSessEstab->TxdRxd.pTxdBuff);
		memset(pSessEstab, 0, sizeof(tsBKSSessEstab));
		pSessEstab->TxdRxd.Socket = -1;
	#endif
		pSessEstab->TxdRxd.TxFrameID = 1;
		pSessEstab->TxdRxd.RxdRxFrameID = 0;
		pSessEstab->TxdRxd.RxdTxFrameID = 0;
		}
	free(m_pBKSSessEstabs);
	m_pBKSSessEstabs = NULL;
	}
m_NumSessEstabs = 0;


if(m_pBKSTypes != NULL)
	{
	if(m_NumSessions && m_NumInitTypes)
		{
		// at least one Session registered, iterate over all Sessions and after closing Session socket connection then delete the Session
		pType = m_pBKSTypes;
		for (Idx = 0; Idx < eBKSPTPlaceHolder-1; Idx++, pType += 1)
			{
			if (pType->Detail.BKSPType != eBKSPTUndefined)
				{
				if ((pSession = pType->pFirstSession) != NULL)
					{
					pType->pFirstSession = NULL;
					pType->NumInstances = 0;
					pType->NumSessions = 0;
					do
						{
	#ifdef _WIN32
						if (pSession->TxdRxd.Socket != INVALID_SOCKET)
							{
							closesocket(pSession->TxdRxd.Socket);
							pSession->TxdRxd.Socket = INVALID_SOCKET;
							}
	#else
						if (pSession->TxdRxd.Socket != -1)
							{
							close(pSession->TxdRxd.Socket);
							pSession->TxdRxd.Socket = -1;
							}
	#endif
						if (pSession->TxdRxd.pRxdBuff != NULL)
							free(pSession->TxdRxd.pRxdBuff);
						if (pSession->TxdRxd.pTxdBuff != NULL)
							free(pSession->TxdRxd.pTxdBuff);
						if(pSession->pReqResp != NULL)
							free(pSession->pReqResp);
						pNext = pSession->pNext;
						free(pSession);
						}
					while ((pSession = pNext) != NULL);
					}
				}
			}
		m_NumSessions = 0;
		}
	free(m_pBKSTypes);
	m_pBKSTypes = NULL;
	}

m_NumInitTypes = 0;
	
pTxdRxd = &m_Ctrl[0];
for(Idx = 0; Idx < 2; Idx+=1, pTxdRxd+=1)
	{
#ifdef _WIN32
	if(pTxdRxd->Socket != INVALID_SOCKET)
		closesocket(pTxdRxd->Socket);
	if (pTxdRxd->pRxdBuff != NULL)
		free(pTxdRxd->pRxdBuff);
	if (pTxdRxd->pTxdBuff != NULL)
		free(pTxdRxd->pTxdBuff);
	memset(pTxdRxd, 0, sizeof(tsTxdRxd));
	pTxdRxd->Socket = INVALID_SOCKET;
#else
	if (pTxdRxd->Socket != -1)
		close(pTxdRxd->Socket);
	if (pTxdRxd->pRxdBuff != NULL)
		free(pTxdRxd->pRxdBuff);
	if (pTxdRxd->pTxdBuff != NULL)
		free(pTxdRxd->pTxdBuff);
	memset(pTxdRxd, 0, sizeof(tsTxdRxd));
	pTxdRxd->Socket = -1;
#endif
	}

#ifdef WIN32
if(m_ListenerSock != INVALID_SOCKET)
	{
	closesocket(m_ListenerSock);
	m_ListenerSock = INVALID_SOCKET;
	}
#else
if(m_ListenerSock != -1)
	{
	close(m_ListenerSock);
	m_ListenerSock = -1;
	}
#endif

if(m_ppChkPtReqs != NULL)
	{
	tsChkPtReqs *pChkPtReqs;
	for(Idx=0; Idx < m_MaxChkPtReqs; Idx++)
		{
		pChkPtReqs = m_ppChkPtReqs[Idx];
		if(pChkPtReqs != NULL)
			{
			if(pChkPtReqs->pRespData != NULL)
				{
#ifdef _WIN32
				free(pChkPtReqs->pRespData);	
#else
				if(pChkPtReqs->pRespData != MAP_FAILED)
					munmap(pChkPtReqs->pRespData,pChkPtReqs->AllocRespData);
#endif
				pChkPtReqs->pRespData = NULL;
				}
#ifdef _WIN32
			free(pChkPtReqs->pRespData);	
#else
			if(pChkPtReqs->pRespData != MAP_FAILED)
				munmap(pChkPtReqs->pRespData,pChkPtReqs->AllocRespData);
#endif
			pChkPtReqs->pRespData = NULL;
			}
		}
	free(m_ppChkPtReqs);
	m_ppChkPtReqs = NULL;
	}
m_NumChkPtReqs = 0;	
m_MaxChkPtReqs = 0;					

memset(m_SessionIDVect, 0, sizeof(m_SessionIDVect));
memset(m_ReqIDVect, 0, sizeof(m_ReqIDVect));
m_CASSerialise = 0;
m_NumPendingReqs = 0;
m_NumPendingResps = 0;
m_TotRespsAvail = 0;
m_bSessionTermReq = false;
return(eBSFSuccess);
}

int				// returns number of sessions deleted
CBKSRequester::DeleteAllSessionsInState(teBKSPEPProvState TermState)		// terminate and delete any session in specific state (normally eBKSPSRegisteredTerm) or with an invalid socket
{
int Idx;
int NumDeleted;
tsBKSType *pType;
tsBKSRegSessionEx *pSession;

tsBKSRegSessionEx *pNext;
tsBKSRegSessionEx *pPrev;
tsReqRespInst *pInstance;
UINT32 InstanceID;

NumDeleted = 0;
if (m_NumSessions)
	{
	// at least one Session registered, iterate over all sessions and if it's state is TermState then delete that session
	pType = m_pBKSTypes;
	for (Idx = 0; Idx < eBKSPTPlaceHolder - 1; Idx++, pType += 1)
		{
		if (pType->Detail.BKSPType != eBKSPTUndefined)
			{
			if ((pSession = pType->pFirstSession) != NULL)
				{
				pPrev = NULL;
				do
					{
					pNext = pSession->pNext;
					if(pSession->Session.BKSPState == TermState ||
#ifdef _WIN32
					   pSession->TxdRxd.Socket == INVALID_SOCKET)
#else
					   pSession->TxdRxd.Socket == -1)
		
#endif
						{
#ifdef _WIN32
						if (pSession->TxdRxd.Socket != INVALID_SOCKET)
							{
							closesocket(pSession->TxdRxd.Socket);
							pSession->TxdRxd.Socket = INVALID_SOCKET;
							}
#else
						if (pSession->TxdRxd.Socket != -1)
							{
							close(pSession->TxdRxd.Socket);
							pSession->TxdRxd.Socket = -1;
							}
#endif
						if (pSession->TxdRxd.pRxdBuff != NULL)
							free(pSession->TxdRxd.pRxdBuff);
						if (pSession->TxdRxd.pTxdBuff != NULL)
							free(pSession->TxdRxd.pTxdBuff);


						if(pSession->Session.NumCpltd != 0)
#ifdef WIN32
							InterlockedAdd(&m_TotRespsAvail,(int)pSession->Session.NumCpltd);	  	
#else
							__sync_fetch_and_sub(&m_TotRespsAvail,pSession->Session.NumCpltd);
#endif

						if (pSession->pReqResp != NULL)
							{
							pInstance = (tsReqRespInst *)pSession->pReqResp;
							for(InstanceID = 1; InstanceID <= pSession->Session.MaxInstances; InstanceID++)
								{
								if(pInstance->ReqID > 0)
									UnallocReqID(pInstance->ReqID);
								pInstance = (tsReqRespInst *)((UINT8 *)pInstance + pType->ReqRespInstSize);
								}
							free(pSession->pReqResp);
							}
					    if(pSession->Session.SessionID != 0)
							{
							UnallocSessionID(pSession->Session.SessionID);
							pSession->Session.SessionID = 0;
							pSession->TxdRxd.SessionID = 0;
							}
						pType->pSessions[pSession->Session.TypeSessionID-1] = NULL;
					
						if(pPrev != NULL)
							pPrev->pNext = pNext;

						if(pType->pFirstSession == pSession)
							pType->pFirstSession = pNext;

						pType->NumInstances -= pSession->Session.MaxInstances;
						pType->NumSessions -= 1;
						m_NumSessions -= 1;
				
						free(pSession);
						NumDeleted += 1;
						}
					else
						pPrev = pSession;
					}
				while ((pSession = pNext) != NULL);
				}
			}
		}
	}

return(NumDeleted);
}

int 
CBKSRequester::Reset(bool bReleaseLock)			// locked mutex requires releasing immediately prior to deleting the mutexes
{
TerminateAllSessions();

#ifdef WIN32
if(InterlockedCompareExchange(&m_CreatedMutexes,1,1)==1)
#else
if(__sync_val_compare_and_swap(&m_CreatedMutexes,1,1)==1)
#endif
	{
	if(bReleaseLock)
		ReleaseLock(true);
	DeleteMutexes();
	}
return(eBSFSuccess);
}

int					// returns total number of registered service types or teBSFrsltCodes error code if any parameterisation errors or already registered type
CBKSRequester::RegServiceType(teBKSPType BKSPType,			// registering this service type
							  UINT32 MinProviderVersion,	// service provider version must be at least this software version
							  UINT32 MaxProviderVersion,	// service provider version must be no more than this software version
							  UINT32 KeepaliveSecs,			// expecting packet activity from session peer with periodicity of no more than this number of seconds
							  UINT32 MaxSessions,			// allowing at most this many concurrent sessions of specified service type
							  UINT32 MinServiceInsts,		// any session to support at least this minimum number of service instances
							  UINT32 MaxServiceInsts,		// limit any session to support a maximum of this many service instances
							  UINT32 MaxProcSecs,           // expecting any service provider to take no more than this number of seconds to complete processing a given request
							  UINT32 MaxParamLen,			// service job parameters can be up to this length
							  UINT32 MaxQuerySeqLen,		// query sequences can be up to this length
							  UINT32 MaxTargSeqLen,			// target sequences can be up to this length
							  UINT32 MaxReqPayloadSize,		// request payloads to the service provider, including framing, can be up to this size (UINT8s)
							  UINT32 MaxRespPayloadSize)	// response payloads from the service provider, including framing, can be up to this  size (UINT8s)

{
tsBKSType *pType;
int NumInitTypes;

if(m_pBKSTypes == NULL)
	return(-1);

// ensure parameter values are acceptable
if(BKSPType <= eBKSPTUndefined || BKSPType >= eBKSPTPlaceHolder ||
   MaxSessions < 1 || MaxSessions > cMaxNumSessions ||
   KeepaliveSecs < cMinKeepaliveSecs || KeepaliveSecs > cMaxKeepaliveSecs ||
   MinServiceInsts < 1 || MinServiceInsts > MaxServiceInsts || MaxServiceInsts > cMaxServiceInsts ||
   MaxProcSecs > cAbsMaxProcSecs ||
   MaxParamLen > cAbsMaxParamLen ||
   MaxReqPayloadSize < 1 || MaxReqPayloadSize > cAbsMaxReqPayloadSize ||
   MaxRespPayloadSize < 1 || MaxRespPayloadSize > cAbsMaxRespPayloadSize ||
   MaxQuerySeqLen > cAbsMaxQuerySeqLen || MaxTargSeqLen > cAbsMaxTargSeqLen)
   return(eBSFerrParams);

if(MaxProcSecs < cAbsMinProcSecs)
	MaxProcSecs = cAbsMinProcSecs;

pType = &m_pBKSTypes[BKSPType-1];
if(pType->Detail.BKSPType != eBKSPTUndefined)				// not allowed to reinitialise if already initialised
	return(eBSFerrInternal);

memset(pType,0,sizeof(tsBKSType));
pType->Detail.BKSPType = BKSPType;
pType->Detail.MinProviderVersion = MinProviderVersion;
pType->Detail.MaxProviderVersion = MaxProviderVersion;
pType->Detail.KeepaliveSecs = KeepaliveSecs;
pType->Detail.Priority = 100;
pType->Detail.MinServiceInsts = MinServiceInsts;
pType->Detail.MaxServiceInsts = MaxServiceInsts;
pType->Detail.MaxProcSecs = MaxProcSecs;
pType->Detail.MaxParamLen = MaxParamLen;
pType->Detail.MaxQuerySeqLen = MaxQuerySeqLen;
pType->Detail.MaxTargSeqLen = MaxTargSeqLen;
pType->Detail.MaxReqPayloadSize = MaxReqPayloadSize;
pType->Detail.MaxRespPayloadSize = MaxRespPayloadSize;
pType->MaxSessions = MaxSessions;
pType->ReqRespInstSize = (UINT32)sizeof(tsReqRespInst) + max(MaxReqPayloadSize,MaxRespPayloadSize);
pType->pFirstSession = NULL;
m_NumInitTypes += 1;
NumInitTypes = m_NumInitTypes;
return(NumInitTypes);
}

// session registration process is in multiple phases
// a) tsBKSPacNegA: In the initial phase, immediately following connection acceptance, the putative provider is sent the list of acceptable services
// b) tsBKSPacNegB: the provider then responds with the service it can provide and number of service instances supported
// c) tsBKSPacNegC: the provider is then sent the service type as an acknowledgment of acceptance by the server
// NOTE: any received keepalives during the negotiation phases will be silently discarded
bool
CBKSRequester::ProgressSessEstab(tsBKSSessEstab *pSessEstab,
				 bool bCpltdWrite)	// false if frame received, true if frame sent
{
bool bRegistered;
tsBKSPacHdr *pHdr;
tsBKSOfferedService *pOfferService;
tsBKSAcceptService *pAcceptService;
teBKSPType Type;
tsBKSType *pBKSType;
UINT32 ServiceInsts;
UINT32 ClassInsts;
UINT32 Idx;
UINT32 Diff;

if(pSessEstab == NULL || pSessEstab->SEState == eSESnone ||
#ifdef WIN32
   pSessEstab->TxdRxd.Socket == INVALID_SOCKET)
#else
pSessEstab->TxdRxd.Socket == -1)
#endif
	{
	ResetSessEstab(pSessEstab);
	return(false);
	}

if(!bCpltdWrite)								// if a response frame has been received then ensure it is for the expected session and correct response type
	{
	pHdr = (tsBKSPacHdr *)pSessEstab->TxdRxd.pRxdBuff;
	if(pHdr->FrameType != eBKSHdrOfferedService ||
		pHdr->SessionID != pSessEstab->TxdRxd.SessionID ||
	   pSessEstab->TxdRxd.CurPacRxd != sizeof(tsBKSOfferedService))
		{
		ResetSessEstab(pSessEstab);
		m_NumSessEstabs -= 1;
		return(false);
		}
	pOfferService = (tsBKSOfferedService *)pHdr;
	}

switch(pSessEstab->SEState) {
	case eSESTxReqServices:				// sending list of required services to service provider, expecting frame to have been sent
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProgressSessEstab with session: %u sending list of requested services", pSessEstab->TxdRxd.SessionID);
		if(!bCpltdWrite)				// whilst frame was being sent then not expecting any received frames
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProgressSessEstab with session: %u unexpected received frame", pSessEstab->TxdRxd.SessionID);
			ResetSessEstab(pSessEstab);
			m_NumSessEstabs -= 1;
			return(false);
			}
		pSessEstab->TxdRxd.flgTxCplt = 0;
		pSessEstab->TxdRxd.CurTxd = 0;
		pSessEstab->TxdRxd.TotTxd = 0;
		pSessEstab->TxdRxd.flgSelMonExcept = 1;
		pSessEstab->TxdRxd.flgSelMonRead = 1;
		pSessEstab->TxdRxd.flgSelMonWrite = 0;
		pSessEstab->SEState = eSESRxOfferedService;			
		return(true);

	case eSESRxOfferedService:			// expecting to have received type of service offered in a tsBKSPacNegB frame
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProgressSessEstab with session: %u received offered service", pSessEstab->TxdRxd.SessionID);
		ServiceInsts = pOfferService->ServiceInsts;
		ClassInsts = pOfferService->ClassInsts;
		Type = pOfferService->BKSPType;
		if(Type != eBKSPTUndefined)
			{
			pBKSType = m_pBKSTypes;
			for(Idx = 0; Idx < eBKSPTPlaceHolder - 1; Idx+=1, pBKSType += 1)
				if(Type == pBKSType->Detail.BKSPType)
					break;
			}
		if(Idx == eBKSPTPlaceHolder-1 || Type == eBKSPTUndefined)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProgressSessEstab with session: %u received offered service but unable to match", pSessEstab->TxdRxd.SessionID);
			ResetSessEstab(pSessEstab);
			m_NumSessEstabs -= 1;
			return(false);
			}

		pSessEstab->TxdRxd.flgRxCplt = 0;
		if ((Diff = (pSessEstab->TxdRxd.TotRxd - pSessEstab->TxdRxd.CurPacRxd)) > 0)
			{
			memmove(pSessEstab->TxdRxd.pRxdBuff, &pSessEstab->TxdRxd.pRxdBuff[pSessEstab->TxdRxd.CurPacRxd], Diff);
			pSessEstab->TxdRxd.TotRxd = Diff;
			}
		else
			pSessEstab->TxdRxd.TotRxd = 0;
		pSessEstab->TxdRxd.CurPacRxd = 0;
		pSessEstab->TxdRxd.flgRxCplt = 0;

		pSessEstab->TxdRxd.CurTxd = 0;
		pSessEstab->TxdRxd.TotTxd = sizeof(tsBKSAcceptService);
		pAcceptService = (tsBKSAcceptService *)pSessEstab->TxdRxd.pTxdBuff;
		pAcceptService->Hdr.SessionID = pSessEstab->TxdRxd.SessionID;
		pAcceptService->Hdr.FrameFlags = 0;
		pAcceptService->Hdr.RxFrameID = pSessEstab->TxdRxd.RxdTxFrameID;
		pAcceptService->Hdr.TxFrameID = pSessEstab->TxdRxd.TxFrameID++;
		if(pSessEstab->TxdRxd.TxFrameID > 0x07f)
			pSessEstab->TxdRxd.TxFrameID = 1;
		pAcceptService->Hdr.FrameType = eBKSHdrAcceptService;
		pAcceptService->Hdr.FrameLen = sizeof(tsBKSAcceptService);

		if(pBKSType->NumSessions >= pBKSType->MaxSessions) // check if can support an additional Session for this type
			{
			pSessEstab->SEState = eSESTxRejectService;
			pSessEstab->TxdRxd.flgSelMonExcept = 1;
			pSessEstab->TxdRxd.flgSelMonRead = 0;
			pSessEstab->TxdRxd.flgSelMonWrite = 1;
			pAcceptService->BKSPType = eBKSPTUndefined;		// can't accept
			}
		else
			{
			pSessEstab->SEState = eSESTxAcceptService;
			pSessEstab->TxdRxd.flgSelMonExcept = 1;
			pSessEstab->TxdRxd.flgSelMonRead = 1;
			pSessEstab->TxdRxd.flgSelMonWrite = 1;
			pAcceptService->BKSPType = Type;
			pSessEstab->BKSPType = Type;
			pSessEstab->MaxInstances = ServiceInsts;
			pSessEstab->MaxClassInstances = ClassInsts;
			}

		if (!TxData(&pSessEstab->TxdRxd))
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProgressSessEstab with session: %u unable to send acceptance", pSessEstab->TxdRxd.SessionID);

			ResetSessEstab(pSessEstab);
			m_NumSessEstabs -= 1;
			return(false);
			}
		if (pSessEstab->TxdRxd.flgTxCplt != 1)		// wait for packet to be sent before registering the Session
			return(true);
		// deliberate fallthrough - acceptance or rejection was sent
	case eSESTxAcceptService:		// completed sending service offer acceptance
	case eSESTxRejectService:		// completed sending service offer rejection
		if(pSessEstab->SEState == eSESTxRejectService)	// if service was rejected then close down session
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProgressSessEstab with session: %u rejection", pSessEstab->TxdRxd.SessionID);
			ResetSessEstab(pSessEstab);
			m_NumSessEstabs -= 1;
			return(true);
			}
		// now accepting as a full session
		bRegistered = AcceptFullSession(pSessEstab);
		if(bRegistered == true)
			 gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProgressSessEstab with session: %u providing %d class instances is accepted as full session", pSessEstab->TxdRxd.SessionID,pSessEstab->MaxClassInstances);
		else
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProgressSessEstab with session: %u providing %d class instances, AcceptFullSession() returned %s", pSessEstab->TxdRxd.SessionID, pSessEstab->MaxClassInstances, bRegistered == true ? "True" : "False");
		ResetSessEstab(pSessEstab, true, bRegistered ? true : false);
		m_NumSessEstabs -= 1;
		return(bRegistered);

	default:						// if any other state then close session
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProgressSessEstab with session: %u invalid session state", pSessEstab->TxdRxd.SessionID);
		ResetSessEstab(pSessEstab, true);
		m_NumSessEstabs -= 1;
		return(false);
	}

return(true);
}

bool
CBKSRequester::AcceptFullSession(tsBKSSessEstab *pSessEstab) // accepting session being established as full session ready for normal payload request/responses
{
tsBKSType *pType;
UINT32 NumSessions;
tsBKSRegSessionEx *pSession;
UINT32 TypeSessionID;

if((pSession = (tsBKSRegSessionEx *)malloc(sizeof(tsBKSRegSessionEx)))==NULL)
	return(false);
memset(pSession,0,sizeof(tsBKSRegSessionEx));
pSession->TxdRxd = pSessEstab->TxdRxd;
#ifdef WIN32
pSessEstab->TxdRxd.Socket = INVALID_SOCKET;
#else
pSessEstab->TxdRxd.Socket = -1;
#endif

pType = &m_pBKSTypes[pSessEstab->BKSPType - 1];

if(pType->NumSessions == pType->MaxSessions)
	return(false);

// calc memory required buffering this Sessions instances
size_t MemReq;

pSession->Session.MaxInstances = pSessEstab->MaxInstances;
pSession->Session.MaxClassInstances = pSessEstab->MaxClassInstances;
MemReq = pType->Detail.MaxRespPayloadSize * 2;			// buffering up to 2 responses to try and hide some potential network latency
if ((pSession->TxdRxd.pRxdBuff = (UINT8 *)malloc(MemReq)) == NULL)
	return(false);

pSession->TxdRxd.AllocdRxdBuff = (UINT32)MemReq;

MemReq = pType->Detail.MaxReqPayloadSize * 2;			// buffering 2 requests to try and hide any potential network latency
if ((pSession->TxdRxd.pTxdBuff = (UINT8 *)malloc(MemReq)) == NULL)
	{
	free(pSession->TxdRxd.pRxdBuff);
	pSession->TxdRxd.pRxdBuff = NULL;
	return(false);
	}
pSession->TxdRxd.AllocdTxdBuff = (UINT32)MemReq;

MemReq = (pType->Detail.MaxReqPayloadSize + pType->Detail.MaxRespPayloadSize + 1000) * pSessEstab->MaxInstances;
if ((pSession->pReqResp = (UINT8 *)malloc(MemReq)) == NULL)
	{
	free(pSession->TxdRxd.pRxdBuff);
	pSession->TxdRxd.pRxdBuff = NULL;
	free(pSession->TxdRxd.pRxdBuff);
	pSession->TxdRxd.pRxdBuff = NULL;
	return(false);
	}
memset(pSession->pReqResp,0,MemReq);
pSession->AllocdReqResp = (UINT32)MemReq;
pSession->Session.BKSPState = eBKSPSRegisteredActv;
pSession->Session.BKSPType = pSessEstab->BKSPType;
pSession->Session.MaxInstances = pSessEstab->MaxInstances;
pSession->Session.MaxClassInstances = pSessEstab->MaxClassInstances;
pSession->Session.SessionID = pSession->TxdRxd.SessionID;

if(pSessEstab->TxdRxd.TotRxd > 0)
	memcpy(pSession->TxdRxd.pRxdBuff, pSessEstab->TxdRxd.pRxdBuff, pSessEstab->TxdRxd.TotRxd);
if (pSessEstab->TxdRxd.TotTxd > 0)
	memcpy(pSession->TxdRxd.pTxdBuff, pSessEstab->TxdRxd.pTxdBuff, pSessEstab->TxdRxd.CurTxd);
pType->NumSessions += 1;
pType->NumInstances += pSessEstab->MaxInstances;
for(TypeSessionID = 1; TypeSessionID <= pType->MaxSessions; TypeSessionID += 1)
	{
	if(pType->pSessions[TypeSessionID-1] != NULL)
		continue;
	pType->pSessions[TypeSessionID - 1] = pSession;
	pSession->Session.TypeSessionID = TypeSessionID;
	break;
	}

pSession->TxdRxd.flgSelMonExcept = 1;
pSession->TxdRxd.flgSelMonRead = 1;
pSession->TxdRxd.flgSelMonWrite = 0;

pSession->pNext = pType->pFirstSession;
pType->pFirstSession = pSession;
m_NumSessions += 1;
NumSessions = m_NumSessions;
return(true);
}

//
// Note: if for any reason a session establishment can't be started then the just accepted socket will be closed and false returned
// Note: if there are already cMaxNumSessions of accepted Sessions then the just accepted socket will be closed and false returned
// Note: m_NumSessEstabs will be incremented when session establishment has been started, unless a time expiring session establishment is reused

void
CBKSRequester::ResetSessEstab(tsBKSSessEstab *pSessEstab,			// reset this tsBKSSessEstab instance
				bool bKeepAllocs,									// if true then retain any existing buffer allocations
 			  bool bKeepSessionID)					// if true then retain SessionID as session has been accepted

{
UINT8 *pRxdBuff;
UINT32 AllocdRxdBuff;
UINT8 *pTxdBuff;
UINT32 AllocdTxdBuff;
if(pSessEstab == NULL)
	return;
if(pSessEstab->TxdRxd.SessionID != 0 && !bKeepSessionID)
	UnallocSessionID(pSessEstab->TxdRxd.SessionID);
pSessEstab->TxdRxd.SessionID = 0;
if(bKeepAllocs)
	{
	pRxdBuff = pSessEstab->TxdRxd.pRxdBuff;
	AllocdRxdBuff = pSessEstab->TxdRxd.AllocdRxdBuff;
	pTxdBuff = pSessEstab->TxdRxd.pTxdBuff;
	AllocdTxdBuff = pSessEstab->TxdRxd.AllocdTxdBuff;
	}
else
	{
	if(pSessEstab->TxdRxd.pRxdBuff != NULL)
		{
		free(pSessEstab->TxdRxd.pRxdBuff);
		pSessEstab->TxdRxd.pRxdBuff = NULL;
		}
	AllocdRxdBuff = 0;
	if (pSessEstab->TxdRxd.pTxdBuff != NULL)
		{
		free(pSessEstab->TxdRxd.pTxdBuff);
		pSessEstab->TxdRxd.pTxdBuff = NULL;
		}
	AllocdTxdBuff = 0;
	}
#ifdef WIN32
if (pSessEstab->TxdRxd.Socket != INVALID_SOCKET)
	closesocket(pSessEstab->TxdRxd.Socket);
#else
if (pSessEstab->TxdRxd.Socket != -1)
	close(pSessEstab->TxdRxd.Socket);
#endif
memset(pSessEstab,0,sizeof(tsBKSSessEstab));
pSessEstab->TxdRxd.pRxdBuff = pRxdBuff;
pSessEstab->TxdRxd.AllocdRxdBuff = AllocdRxdBuff;
pSessEstab->TxdRxd.pTxdBuff = pTxdBuff;
pSessEstab->TxdRxd.AllocdTxdBuff = AllocdTxdBuff;
#ifdef WIN32
pSessEstab->TxdRxd.Socket = INVALID_SOCKET;
#else
pSessEstab->TxdRxd.Socket = -1;
#endif
pSessEstab->TxdRxd.TxFrameID = 1;
pSessEstab->TxdRxd.RxdRxFrameID = 0;
pSessEstab->TxdRxd.RxdTxFrameID = 0;
}

bool
CBKSRequester::StartSessEstab(UINT32 SessionID,		// session identifier for this potential Session session
				socket_t Socket,					// communicating with Session over this connected socket
				SOCKADDR_STORAGE  *pIPaddress)		// Session is at this socket network address
{
bool bReused;
UINT32 Idx;
time_t Now;
time_t Then;
tsBKSSessEstab *pSessEstab;
tsBKSSessEstab *pOldestSessEstab;
tsServiceDetail *pDetail;
tsBKSReqServices *PacNegA;
tsBKSType *pType;


#ifdef _WIN32
if(Socket == INVALID_SOCKET)
#else
if(Socket == -1)
#endif
	{
	if(SessionID != 0)
		UnallocSessionID(SessionID);
	return(false);
	}

if(SessionID == 0 || pIPaddress == NULL ||
   m_NumSessions >= cMaxNumSessions)				// no point in negotiating if already at the limit of registered Sessions
	{
#ifdef WIN32
	closesocket(Socket);
#else
	close(Socket);
#endif
	if(SessionID != 0)
		UnallocSessionID(SessionID);
	return(false);
	}

// find unused, or if existing session establishment has taken too long then reuse that one 
Now = time(NULL);
Then = Now - 10;
pOldestSessEstab = NULL;
if((pSessEstab = m_pBKSSessEstabs)==NULL)
	{
	Reset(true);
	return(false);
	}
for(Idx = 0; Idx < cMaxSessEstab; Idx+=1, pSessEstab += 1)
	{
	if(pSessEstab->SEState == eSESnone)					// 0 if uncommitted
		break;
	if(pSessEstab->StartSecs < Then)
		{
		Then = pSessEstab->StartSecs;
		pOldestSessEstab = pSessEstab;
		}
	}
if(Idx == cMaxSessEstab)
	{
	if(pOldestSessEstab == NULL)
		{
#ifdef WIN32
		closesocket(Socket);
#else
		close(Socket);
#endif
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "StartSessEstab: Rejected, no free session establishment entries");
		UnallocSessionID(SessionID);
		return(false);
		}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "StartSessEstab: reusing previously timed out session establishment entry");
	pSessEstab = pOldestSessEstab;
	bReused = true;
	}
else
	bReused = false;

ResetSessEstab(pSessEstab,true);

if(pSessEstab->TxdRxd.pRxdBuff == NULL)
	{
	if((pSessEstab->TxdRxd.pRxdBuff = (UINT8 *)malloc(cMinTxRxBuffSize))==NULL)
		{
		if(bReused)
			m_NumSessEstabs -= 1;
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "StartSessEstab: unable to allocate memory for rxd buffering");
		UnallocSessionID(SessionID);
		return(false);
		}
	pSessEstab->TxdRxd.AllocdRxdBuff = cMinTxRxBuffSize;
	}
if (pSessEstab->TxdRxd.pTxdBuff == NULL)
	{
	if ((pSessEstab->TxdRxd.pTxdBuff = (UINT8 *)malloc(cMinTxRxBuffSize)) == NULL)
		{
		free(pSessEstab->TxdRxd.pRxdBuff);
		pSessEstab->TxdRxd.AllocdRxdBuff = 0;
		if (bReused)
			m_NumSessEstabs -= 1;
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "StartSessEstab: unable to allocate memory for txd buffering");
		UnallocSessionID(SessionID);
		return(false);
		}
	pSessEstab->TxdRxd.AllocdTxdBuff = cMinTxRxBuffSize;
	}

pSessEstab->TxdRxd.TxFrameID = 1;
pSessEstab->TxdRxd.RxdRxFrameID = 0;
pSessEstab->TxdRxd.RxdTxFrameID = 0;
pSessEstab->TxdRxd.SessionID = SessionID;
pSessEstab->TxdRxd.IPaddress = *pIPaddress;
pSessEstab->TxdRxd.Socket = Socket;
pSessEstab->StartSecs = Now;
// send list of required services to provider Session
PacNegA = (tsBKSReqServices *)pSessEstab->TxdRxd.pTxdBuff;
PacNegA->Hdr.FrameLen = sizeof(tsBKSReqServices);
PacNegA->Hdr.TxFrameID = pSessEstab->TxdRxd.TxFrameID++;
PacNegA->Hdr.RxFrameID = pSessEstab->TxdRxd.RxdTxFrameID;
PacNegA->Hdr.SessionID = SessionID;
PacNegA->Hdr.FrameFlags = 0;
PacNegA->Hdr.FrameType = eBKSHdrReqServices;
PacNegA->NumTypes = 0;
pDetail = &PacNegA->Details[0];
pType = m_pBKSTypes;
for(Idx = 0; Idx < eBKSPTPlaceHolder - 1; Idx++, pType+=1)
	{
	if(pType->NumSessions < pType->MaxSessions)
		{
		*pDetail++ = pType->Detail;
		if(PacNegA->NumTypes >= 1)
			PacNegA->Hdr.FrameLen += sizeof(tsServiceDetail);
		PacNegA->NumTypes += 1;
		}
	}
pSessEstab->TxdRxd.flgSelMonExcept = 1;
pSessEstab->TxdRxd.flgSelMonRead = 0;
pSessEstab->TxdRxd.flgSelMonWrite = 1;
pSessEstab->TxdRxd.TotTxd = PacNegA->Hdr.FrameLen;
pSessEstab->SEState = eSESTxReqServices;
if(!bReused)
	m_NumSessEstabs += 1;

if(!TxData(&pSessEstab->TxdRxd))
	{
	ResetSessEstab(pSessEstab, true);
	m_NumSessEstabs -= 1;
	UnallocSessionID(SessionID);
	return(false);
	}
if(pSessEstab->TxdRxd.flgTxCplt = 1)		// if was all sent then onto next phase whereby a response is to be expected
	{
	pSessEstab->TxdRxd.flgTxCplt = 0;
	pSessEstab->TxdRxd.CurTxd = 0;
	pSessEstab->TxdRxd.TotTxd = 0;
	pSessEstab->TxdRxd.flgSelMonExcept = 1;
	pSessEstab->TxdRxd.flgSelMonRead = 1;
	pSessEstab->TxdRxd.flgSelMonWrite = 0;
	pSessEstab->SEState = eSESRxOfferedService;
	}
return(true);
}

tsBKSRegSessionEx *
CBKSRequester::LocateSession(UINT32 SessionID) // locates an established session which is identified by SessionID
{
tsBKSRegSessionEx *pSession;
tsBKSType *pType;
int TypeIdx;

if(SessionID == 0 || !IsAllocSessionID(SessionID) ||
   m_NumSessions == 0)
	return(NULL);

// at least one Session registered, iterate over all Sessions looking for a match
pType = m_pBKSTypes;
for (TypeIdx = 0; TypeIdx < eBKSPTPlaceHolder-1; TypeIdx++, pType += 1)
	{
	if (pType->Detail.BKSPType != eBKSPTUndefined)
		{
		if ((pSession = pType->pFirstSession) != NULL)
			{
			do
				{
				if (pSession->Session.SessionID == SessionID)
					return(pSession);
				}
			while ((pSession = pSession->pNext) != NULL);
			}
		}
	}
return(NULL);
}


int						// returns number of registered Sessions which are in requested state
CBKSRequester::GetNumSessions(teBKSPType BKSPType,	// may request number of Sessions providing this specific service type, or if eBKSPTUndefined then the total
							teBKSPEPProvState BKSPState)			// Sessions must be in this specific state
{
int NumSessions;
int TypeIdx;
tsBKSType *pType;
tsBKSRegSessionEx *pSession;

if(BKSPType >= eBKSPTPlaceHolder || BKSPState == eBKSPSUndefined || BKSPState == eBKSPSRegisteredTerm)
	return(0);
NumSessions = 0;

AcquireLock(true);			// DEBUG - was false ----not needing exclusive access 
if(m_NumSessions == 0 )		// no registered Sessions?
	{
	ReleaseLock(true);      // DEBUG - was false
	return(0);
	}

pType = m_pBKSTypes;
for (TypeIdx = 0; TypeIdx < eBKSPTPlaceHolder-1; TypeIdx++, pType += 1)
	{
	if(pType->Detail.BKSPType == eBKSPTUndefined || pType->NumSessions == 0)
		continue;
	if(BKSPType == eBKSPTUndefined || pType->Detail.BKSPType == BKSPType)
		{
		if((pSession = pType->pFirstSession)!=NULL)
			{
			do {
				if(pSession->Session.BKSPState == BKSPState)
					NumSessions += 1;
				}
			while((pSession = pSession->pNext)!=NULL);
			}
		if(BKSPType != eBKSPTUndefined)
			break;
		}
	}
ReleaseLock(true);		// DEBUG - was false
return(NumSessions);
}

int						// returns total number of service instances in all registered sessions which are in requested state
CBKSRequester::GetNumInstances(teBKSPType BKSPType,					// may request number of instances providing this specific service type, or if eBKSPTUndefined then the total
				teBKSPEPProvState BKSPState)			// sessions providing service instances must be in this specific state
{
int NumInstances;
int TypeIdx;
tsBKSType *pType;
tsBKSRegSessionEx *pSession;

if (BKSPType >= eBKSPTPlaceHolder || BKSPState == eBKSPSUndefined || BKSPState >= eBKSPSPlaceHolder)
	return(0);
NumInstances = 0;

AcquireLock(true);			// // DEBUG - was false ---------not needing exclusive access 
if (m_NumSessions == 0)		// no registered Sessions?
	{
	ReleaseLock(true);      // DEBUG - was false
	return(0);
	}

pType = m_pBKSTypes;
for (TypeIdx = 0; TypeIdx < eBKSPTPlaceHolder - 1; TypeIdx++, pType += 1)
	{
	if (pType->Detail.BKSPType == eBKSPTUndefined || pType->NumSessions == 0)
		continue;
	if (BKSPType == eBKSPTUndefined || pType->Detail.BKSPType == BKSPType)
		{
		if ((pSession = pType->pFirstSession) != NULL)
			{
			do
				{
				if (pSession->Session.BKSPState == BKSPState)
					NumInstances += pSession->Session.MaxInstances;
				}
			while ((pSession = pSession->pNext) != NULL);
			}
		if (BKSPType != eBKSPTUndefined)
			break;
		}
	}
ReleaseLock(true);		// DEBUG - was false
return(NumInstances);
}

int													// actual number of service provider Sessions returned
CBKSRequester::GetSessions(teBKSPType BKSPType,	// return provider Sessions for this service type,  or if eBKSPTUndefined then all service types registered
							teBKSPEPProvState BKSPState,					// services must be in this state
			 int MaxEPs,							// return at most this many Sessions in pSessions
			 tsBKSRegSession *pSessions)		    // caller has preallocated to hold returned array snapshot of provider Sessions
{
int NumSessions;
int TypeIdx;
tsBKSType *pType;
tsBKSRegSessionEx *pSessionEx;

if(pSessions != NULL && MaxEPs > 0)
	memset(pSessions, 0, sizeof(tsBKSRegSession) * MaxEPs);

if (BKSPType >= eBKSPTPlaceHolder || BKSPState == eBKSPSUndefined || BKSPState >= eBKSPSPlaceHolder || MaxEPs < 1 || pSessions == NULL)
	return(0);

NumSessions = 0;

AcquireLock(true);				// DEBUG - was false // not needing exclusive access 
if (m_NumInitTypes == 0 || m_NumSessions == 0)
	{
	ReleaseLock(true);			// DEBUG - was false
	return(0);
	}

pType = m_pBKSTypes;
for (TypeIdx = 0; TypeIdx < eBKSPTPlaceHolder-1; TypeIdx++, pType += 1)
	{
	if (pType->Detail.BKSPType == eBKSPTUndefined || pType->NumSessions == 0)
		continue;
	if (BKSPType == eBKSPTUndefined || pType->Detail.BKSPType == BKSPType)
		{
		if ((pSessionEx = pType->pFirstSession) != NULL)
			{
			do
				{
				if (pSessionEx->Session.BKSPState == BKSPState)
					{
					*pSessions++ = pSessionEx->Session;
					NumSessions += 1;
					}
				}
			while (NumSessions < MaxEPs && (pSessionEx = pSessionEx->pNext) != NULL);
			}
		if (BKSPType != eBKSPTUndefined || NumSessions == MaxEPs)
			break;
		}
	}
ReleaseLock(true);	// DEBUG - was false
return(NumSessions);
}

#ifdef WIN32
//using namespace std;


// List of Winsock error constants mapped to an interpretation string.
// Note that this list must remain sorted by the error constants'
// values, because we do a binary search on the list when looking up
// items.
typedef struct TAG_sErrorEntry {
	int nID;
	const char* pcMessage;
} tsErrorEntry;

static tsErrorEntry gaErrorList[] = {
	{0,                  "No error"},
	{WSAEINTR,           "Interrupted system call"},
	{WSAEBADF,           "Bad file number"},
	{WSAEACCES,          "Permission denied"},
	{WSAEFAULT,          "Bad address"},
	{WSAEINVAL,          "Invalid argument"},
	{WSAEMFILE,          "Too many open sockets"},
	{ WSAEWOULDBLOCK,     "Operation would block"},
	{ WSAEINPROGRESS,     "Operation now in progress"},
	{ WSAEALREADY,        "Operation already in progress"},
	{ WSAENOTSOCK,        "Socket operation on non-socket"},
	{ WSAEDESTADDRREQ,    "Destination address required"},
	{ WSAEMSGSIZE,        "Message too long"},
	{ WSAEPROTOTYPE,      "Protocol wrong type for socket"},
	{ WSAENOPROTOOPT,     "Bad protocol option"},
	{ WSAEPROTONOSUPPORT, "Protocol not supported"},
	{ WSAESOCKTNOSUPPORT, "Socket type not supported"},
	{ WSAEOPNOTSUPP,      "Operation not supported on socket"},
	{ WSAEPFNOSUPPORT,    "Protocol family not supported"},
	{ WSAEAFNOSUPPORT,    "Address family not supported"},
	{ WSAEADDRINUSE,      "Address already in use"},
	{ WSAEADDRNOTAVAIL,   "Can't assign requested address"},
	{ WSAENETDOWN,        "Network is down"},
	{ WSAENETUNREACH,     "Network is unreachable"},
	{ WSAENETRESET,       "Net connection reset"},
	{ WSAECONNABORTED,    "Software caused connection abort"},
	{ WSAECONNRESET,      "Connection reset by peer"},
	{ WSAENOBUFS,         "No buffer space available"},
	{ WSAEISCONN,         "Socket is already connected"},
	{ WSAENOTCONN,        "Socket is not connected"},
	{ WSAESHUTDOWN,       "Can't send after socket shutdown"},
	{ WSAETOOMANYREFS,    "Too many references, can't splice"},
	{ WSAETIMEDOUT,       "Connection timed out"},
	{ WSAECONNREFUSED,    "Connection refused"},
	{ WSAELOOP,           "Too many levels of symbolic links"},
	{ WSAENAMETOOLONG,    "File name too long"},
	{ WSAEHOSTDOWN,       "Host is down"},
	{ WSAEHOSTUNREACH,    "No route to host"},
	{ WSAENOTEMPTY,       "Directory not empty"},
	{ WSAEPROCLIM,        "Too many processes"},
	{ WSAEUSERS,          "Too many users"},
	{ WSAEDQUOT,          "Disc quota exceeded"},
	{ WSAESTALE,          "Stale NFS file handle"},
	{ WSAEREMOTE,         "Too many levels of remote in path"},
	{ WSASYSNOTREADY,     "Network system is unavailable"},
	{ WSAVERNOTSUPPORTED, "Winsock version out of range"},
	{ WSANOTINITIALISED,  "WSAStartup not yet called"},
	{ WSAEDISCON,         "Graceful shutdown in progress"},
	{ WSAHOST_NOT_FOUND,  "Host not found"},
	{ WSANO_DATA,         "No host data of that type was found"},

};
const int kNumMessages = sizeof(gaErrorList) / sizeof(tsErrorEntry);


const char*
CBKSRequester::WSAGetLastErrorMessage(const char* pcMessagePrefix,
								   int nErrorID /* = 0 */)
{
int Ofs;
int Idx;
tsErrorEntry *pErrEntry;

	// Build basic error string
static char acErrorBuffer[256];
Ofs = sprintf(acErrorBuffer,"%s:", pcMessagePrefix);

	// Tack appropriate canned message onto end of supplied message prefix
pErrEntry = gaErrorList;
for(Idx = 0; Idx < kNumMessages; Idx++, pErrEntry++)
	if(pErrEntry->nID == nErrorID)
		{
		Ofs += sprintf(&acErrorBuffer[Ofs], "%s (%d)", pErrEntry->pcMessage, pErrEntry->nID);
		break;
		}
if(Idx == kNumMessages)
	Ofs += sprintf(&acErrorBuffer[Ofs], "Unkown error (%d)", pErrEntry->nID);

	// Finish error message off and return it.
return acErrorBuffer;
}
#endif

// ShutdownConnection 
// Gracefully shuts the connection Socket down. 
// Returns true if we're successful, false otherwise.
bool		// even if false is returned the socket may have been shutdown, but there was an error whilst shutting it down
CBKSRequester::ShutdownConnection(socket_t *pSocket)
{
char SloughBuffer[0x03fff];
int RxdBytes;
int TotalRxdBytes;
bool bRslt;

if(pSocket == NULL)
	return(false);

#ifdef WIN32
if(*pSocket == INVALID_SOCKET)   // not an error if socket already closed
#else
if (*pSocket == -1)
#endif
	return(true);

#ifdef WIN32
if (shutdown(*pSocket, SD_SEND) == SOCKET_ERROR)
#else
if (shutdown(*pSocket, SHUT_RD) == -1)
#endif
	{
#ifdef WIN32
	closesocket(*pSocket);
	*pSocket = INVALID_SOCKET;
#else
	close(*pSocket);
	*pSocket = -1;
#endif
	return false;
	}

	// Receive and slough any extra data still sitting on the socket
	// sloughing at most 100K bytes
TotalRxdBytes = 0;
while (TotalRxdBytes < 100000)
	{
	RxdBytes = recv(*pSocket, SloughBuffer, sizeof(SloughBuffer), 0);
#ifdef WIN32
	if (RxdBytes == SOCKET_ERROR)
#else
	if (RxdBytes == -1)
#endif
		{
#ifdef WIN32
		closesocket(*pSocket);
		*pSocket = INVALID_SOCKET;
#else
		close(*pSocket);
		*pSocket = -1;
#endif
		return false;
		}
	else 
		if (RxdBytes == 0)
			break;
	TotalRxdBytes += RxdBytes;
	}

	// Close the socket.
#ifdef WIN32
bRslt = closesocket(*pSocket) == SOCKET_ERROR ? false : true;
*pSocket = INVALID_SOCKET;
#else
bRslt = close(*pSocket) == -1 ? false : true;
*pSocket = -1;
#endif
return(bRslt);
}

int				// returns < 0 if errors, eBSFSuccess if initialisation success
CBKSRequester::Initialise(char* pszHost,				// listening on this host/IP address; NULL to use first INET IP local to this machine
						char *pszService,				// listening on this service/port; NULL to use default port 
						teBKSPType BKSPType,			// registering this service type
						UINT32 MinProviderVersion,		// service provider version must be at least this software version
						UINT32 MaxProviderVersion,		// service provider version must be no more than this software version
						UINT32 KeepaliveSecs,			// expecting packet activity from session peer with periodicity of no more than this number of seconds
						UINT32 MaxSessions,				// allowing at most this many concurrent sessions of specified service type
						UINT32 MinServiceInsts,			// any session to support at least this minimum number of service instances
						UINT32 MaxServiceInsts,			// limit any session to support a maximum of this many service instances
						UINT32 MaxProcSecs,             // expecting any service provider to take no more than this number of seconds to complete processing a given request
						UINT32 MaxParamLen,				// service job parameters can be up to this length
  						UINT32 MaxQuerySeqLen,			// query sequences can be up to this length
						UINT32 MaxTargSeqLen,			// target sequences can be up to this length
						UINT32 MaxReqPayloadSize,		// request payloads to the service provider, including framing, can be up to this size (UINT8s),
						UINT32 MaxRespPayloadSize)		// response payloads from the service provider, including framing, can be up to this  size (UINT8s)
{
int Rslt;
UINT32 Idx;
tsBKSSessEstab *pSessEstab;

Reset(false);
CreateMutexes();
AcquireLock(true);

#ifdef WIN32
alignas(8) WSADATA WsaDat;
memset(&WsaDat, 0, sizeof(WsaDat));
if ((Rslt = WSAStartup(MAKEWORD(2, 2), &WsaDat)) != 0)
	{
	Reset(true);
	return(-1);
	}
#endif

m_NumSessEstabs = 0;
if(m_pBKSSessEstabs == NULL)
	{
	m_pBKSSessEstabs = (tsBKSSessEstab *)calloc(cMaxSessEstab,sizeof(tsBKSSessEstab));
	if((pSessEstab = m_pBKSSessEstabs)==NULL)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Initialise: failed to allocate memory (%d bytes) for session negotiations", cMaxSessEstab * sizeof(tsBKSSessEstab)); 
		Reset(true);
		return(-1);
		}
	}

pSessEstab = m_pBKSSessEstabs;
for(Idx = 0; Idx < cMaxSessEstab; Idx++,pSessEstab++)
	{
#ifdef _WIN32
	pSessEstab->TxdRxd.Socket = INVALID_SOCKET;
#else
	pSessEstab->TxdRxd.Socket = -1;
#endif
	pSessEstab->TxdRxd.TxFrameID = 1;
	pSessEstab->TxdRxd.RxdRxFrameID = 0;
	pSessEstab->TxdRxd.RxdTxFrameID = 0;
	}
m_bNotifiedReqs = false;
m_NumInitTypes = 0;

if(m_pBKSTypes == NULL)
	{
	if((m_pBKSTypes = (tsBKSType *)calloc((eBKSPTPlaceHolder-1),sizeof(tsBKSType)))==NULL)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Initialise: failed to allocate memory (%d bytes) for service types", (eBKSPTPlaceHolder-1),sizeof(tsBKSType)); 
		Reset(true);
		return(-1);
		}
	}


if(m_ppChkPtReqs == NULL)
	{
	m_MaxChkPtReqs = cMaxConcurrentRequests;
	if((m_ppChkPtReqs = (tsChkPtReqs **)calloc(cMaxConcurrentRequests,sizeof(tsChkPtReqs *)))==NULL)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Initialise: failed to allocate memory (%d bytes) for concurrent requests", cMaxConcurrentRequests,sizeof(tsChkPtReqs *)); 
		Reset(true);
		return(-1);
		}
	}
m_NumChkPtReqs = 0;

if((Rslt= RegServiceType(BKSPType,MinProviderVersion,MaxProviderVersion,KeepaliveSecs,MaxSessions,MinServiceInsts,MaxServiceInsts,
								MaxProcSecs,MaxParamLen,MaxQuerySeqLen,MaxTargSeqLen,MaxReqPayloadSize,MaxRespPayloadSize)) < 0)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Registration of service type %d failed", BKSPType); 
	Reset(true);
	return(Rslt);
	}
if(pszHost != NULL && pszHost[0] != '\0')
	{
	strncpy(m_szHostName,pszHost,sizeof(m_szHostName));
	m_szHostName[sizeof(m_szHostName)-1] = '\0';
	}
else
	m_szHostName[0] = '\0';
if(pszService != NULL && pszService[0] != '\0')
	{
	strncpy(m_szServiceName,pszService,sizeof(m_szServiceName));
	m_szServiceName[sizeof(m_szServiceName)-1] = '\0';
	}
else
	m_szServiceName[0] = '\0';

if(!InitialiseListener(pszHost, pszService))
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initialise: InitialiseListener for host name '%s' and service name '%s' failed", (pszHost == NULL || pszHost[0] == '\0') ? "not specified" : pszHost,
						 (pszService == NULL || pszService[0] == '\0') ? "not specified" : pszService);
	Reset(true);
	return(-1);
	}

if(listen(m_ListenerSock,cListenBacklog))			// non-zero if errors
	{
#ifdef WIN32
	closesocket(m_ListenerSock);
	m_ListenerSock = INVALID_SOCKET;
#else
	close(m_ListenerSock);
	m_ListenerSock = -1;
#endif
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initialise: listen() with backlog %d failed",cListenBacklog);
	Reset(true);
	return(-1);
	}

if(!InitialiseCtrlSocks())
	{
#ifdef WIN32
	closesocket(m_ListenerSock);
	m_ListenerSock = INVALID_SOCKET;
#else
	close(m_ListenerSock);
	m_ListenerSock = -1;
#endif
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Initialise: InitialiseCtrlSocks() failed");
	Reset(true);
	return(-1);
	}
m_Ctrl[1].flgSelMonExcept = 1;
m_Ctrl[1].flgSelMonRead = 1;

ReleaseLock(true);
return(eBSFSuccess);
}

// control sockets are currently primarily used to ensure that select() is interruptible when there are
// sessions which require requests to be sent to peers for service provision
//
const int cCtrlbuffSize =2048;     // control socket rx/tx buffers are this size
bool 
CBKSRequester::InitialiseCtrlSocks(void) // initialise the control sockets in m_Ctrl[]
{
bool bSuccess;
socket_t ListenerSocket;
socket_t SocketPair[2];
m_bNotifiedReqs = false;

#ifdef WIN32
if(m_Ctrl[0].Socket != INVALID_SOCKET)
	closesocket(m_Ctrl[0].Socket);
if(m_Ctrl[1].Socket != INVALID_SOCKET)
	closesocket(m_Ctrl[1].Socket);
#else
if (m_Ctrl[0].Socket != -1)
	close(m_Ctrl[0].Socket);
if (m_Ctrl[1].Socket != -1)
	close(m_Ctrl[1].Socket);
#endif
if(m_Ctrl[0].pRxdBuff != NULL)
	free(m_Ctrl[0].pRxdBuff);
if(m_Ctrl[0].pTxdBuff != NULL)
	free(m_Ctrl[0].pTxdBuff);
if(m_Ctrl[1].pRxdBuff != NULL)
	free(m_Ctrl[01].pRxdBuff);
if(m_Ctrl[1].pTxdBuff != NULL)
	free(m_Ctrl[1].pTxdBuff);
memset(m_Ctrl,0,sizeof(m_Ctrl));

#ifdef WIN32
m_Ctrl[0].Socket = INVALID_SOCKET;
m_Ctrl[1].Socket = INVALID_SOCKET;

union {
    struct sockaddr_in inaddr;
    struct sockaddr addr;
} Addr;

int LastErr;
socklen_t addrlen = sizeof(Addr.inaddr);
DWORD flags =  0;
int reuse = 1;

ListenerSocket = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
if (ListenerSocket == -1)
    return(false);

memset(&Addr, 0, sizeof(Addr));
Addr.inaddr.sin_family = AF_INET;
Addr.inaddr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
Addr.inaddr.sin_port = 0;

bSuccess = false;
while(!bSuccess) {
    if (setsockopt(ListenerSocket, SOL_SOCKET, SO_REUSEADDR,
            (char*) &reuse, (socklen_t) sizeof(reuse)) == -1)
        break;
    if  (bind(ListenerSocket, &Addr.addr, sizeof(Addr.inaddr)) == SOCKET_ERROR)
        break;

    memset(&Addr, 0, sizeof(Addr));
    if  (getsockname(ListenerSocket, &Addr.addr, &addrlen) == SOCKET_ERROR)
        break;
     Addr.inaddr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    Addr.inaddr.sin_family = AF_INET;

    if (listen(ListenerSocket, 1) == SOCKET_ERROR)
        break;

    SocketPair[0] = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
    if (SocketPair[0] == -1)
        break;
    if (connect(SocketPair[0], &Addr.addr, sizeof(Addr.inaddr)) == SOCKET_ERROR)
        break;

    SocketPair[1] = accept(ListenerSocket, NULL, NULL);
    if (SocketPair[1] == -1)
        break;

    closesocket(ListenerSocket);
	ListenerSocket = INVALID_SOCKET;
	bSuccess = true;
    }


if(!bSuccess)
	{
    LastErr = WSAGetLastError();
    closesocket(ListenerSocket);
    closesocket(SocketPair[0]);
    closesocket(SocketPair[1]);
    WSASetLastError(LastErr);
    return(false);
	}
#else
m_Ctrl[0].Socket = -1;
m_Ctrl[1].Socket = -1;
if(socketpair(AF_UNIX, SOCK_STREAM | SOCK_NONBLOCK, 0, SocketPair) < 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseCtrlSocks: socketpair() error %d", errno);
	return(false);
	}
#endif
m_Ctrl[0].Socket = SocketPair[0];
m_Ctrl[1].Socket = SocketPair[1];

						// need sockets to be non-blocking
#ifdef WIN32
u_long nNoBlock = 1;
ioctlsocket(m_Ctrl[0].Socket, FIONBIO, &nNoBlock);
nNoBlock = 1;
ioctlsocket(m_Ctrl[1].Socket, FIONBIO, &nNoBlock);
#else
int Rsltz = fcntl(m_Ctrl[0].Socket, F_SETFL, fcntl(m_Ctrl[0].Socket,F_GETFL,0) | O_NONBLOCK);
if(Rsltz == -1 || fcntl(m_Ctrl[1].Socket, F_SETFL, fcntl(m_Ctrl[1].Socket,F_GETFL,0) | O_NONBLOCK) == -1)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseCtrlSocks: Unable to make sockets nonblocking");
	return(false);
	}
	
#endif
if((m_Ctrl[0].pRxdBuff = (UINT8 *)calloc(cCtrlbuffSize,1))==NULL)
	return(false);
m_Ctrl[0].AllocdRxdBuff = cCtrlbuffSize;
if((m_Ctrl[1].pRxdBuff = (UINT8 *)calloc(cCtrlbuffSize,1))==NULL)
	return(false);
m_Ctrl[1].AllocdRxdBuff = cCtrlbuffSize;
if((m_Ctrl[0].pTxdBuff = (UINT8 *)calloc(cCtrlbuffSize,1))==NULL)
	return(false);
m_Ctrl[0].AllocdTxdBuff = cCtrlbuffSize;
if((m_Ctrl[1].pTxdBuff = (UINT8 *)calloc(cCtrlbuffSize,1))==NULL)
	return(false);
m_Ctrl[1].AllocdTxdBuff = cCtrlbuffSize;
return(true);
}

bool 
CBKSRequester::ProcessCtrlMsg(int MsgLen,UINT8 *pMsg)			// process a message received by control socket m_Ctrl[1] - note: currently these messages are simply discarded
{
// currently not processing control payloads
m_Ctrl[1].flgRxCplt = 0;
m_bNotifiedReqs = false;
return(true);
}

// InitialiseListener 
bool					// true if server listening port initialised
CBKSRequester::InitialiseListener(const char* pszHost,	// listening on this host/IP address; NULL to use first INET IP local to this machine
								  char *pszService)		// listening on this service/port; NULL to use default port 
{
struct addrinfo AddrInfoHints;
struct addrinfo *pAddrInfoRes;
struct addrinfo *pNxtAddrInfoRes;

char szHost[100];
char szService[100];

socket_t ListenerSocket;
int Rslt;

#ifdef WIN32
if(m_ListenerSock != INVALID_SOCKET)
	{
	closesocket(m_ListenerSock);
	m_ListenerSock = INVALID_SOCKET;
	}
ListenerSocket = INVALID_SOCKET;
#else
if (m_ListenerSock != -1)
	{
	close(m_ListenerSock);
	m_ListenerSock = -1;
	}
m_ListenerSock = -1;
#endif

if (pszService == NULL || pszService[0] == '\0')	// use default port if not specified by caller
	pszService = (char *)cDfltListenerPort;
strncpy(m_szServiceName,pszService,sizeof(m_szServiceName));
m_szServiceName[sizeof(m_szServiceName)-1] = '\0';
memset(&AddrInfoHints, 0, sizeof(struct addrinfo));
AddrInfoHints.ai_family = AF_INET; 		// Return IPv4 and IPv6 choices
AddrInfoHints.ai_socktype = SOCK_STREAM;			 
AddrInfoHints.ai_protocol = IPPROTO_TCP;		

Rslt = getaddrinfo(pszHost == NULL || pszHost[0] == '\0' ? NULL : pszHost, m_szServiceName, &AddrInfoHints, &pAddrInfoRes);
if(Rslt != 0)		// errors if != 0
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseListener: getaddrinfo() failed to resolve host '%s' and service '%s'", pszHost, pszService);
	return(false);
	}


for (pNxtAddrInfoRes = pAddrInfoRes; pNxtAddrInfoRes != NULL; pNxtAddrInfoRes = pNxtAddrInfoRes->ai_next)
	{
	ListenerSocket = socket(pNxtAddrInfoRes->ai_family, pNxtAddrInfoRes->ai_socktype, pNxtAddrInfoRes->ai_protocol);
#ifdef WIN32
	if (ListenerSocket == INVALID_SOCKET)
#else
	if (ListenerSocket == -1)
#endif
		continue;		// try next address

	int SockOptEnable = 1;
	Rslt = setsockopt(ListenerSocket, SOL_SOCKET, SO_REUSEADDR , (char *)&SockOptEnable, sizeof(SockOptEnable));

	if ((Rslt = bind(ListenerSocket, pNxtAddrInfoRes->ai_addr, (int)pNxtAddrInfoRes->ai_addrlen)) == 0)  // 0 if bound successfully
		break;

#ifdef WIN32
	closesocket(ListenerSocket);
	ListenerSocket = INVALID_SOCKET;
#else
	close(ListenerSocket);
	ListenerSocket = -1;
#endif
	}

#ifdef WIN32
if(ListenerSocket == INVALID_SOCKET)
#else
if (ListenerSocket == -1)
#endif
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "InitialiseListener: Unable to bind socket to host '%s' and service '%s'", m_szHostName, m_szServiceName);
	freeaddrinfo(pAddrInfoRes);
	return(false);
	}


// report the listening address
getnameinfo(pNxtAddrInfoRes->ai_addr, (int)pNxtAddrInfoRes->ai_addrlen, szHost, sizeof(szHost), szService, sizeof(szService), 0);
gDiagnostics.DiagOut(eDLInfo, gszProcName, "InitialiseListener: listening for connections to host '%s' and service '%s' with protocol %d", szHost, szService, pNxtAddrInfoRes->ai_protocol);

strncpy(m_szHostName,szHost,sizeof(m_szHostName));
m_szHostName[sizeof(m_szHostName)-1] = '\0';
strncpy(m_szServiceName,szService,sizeof(m_szServiceName));
m_szServiceName[sizeof(m_szServiceName)-1] = '\0';

freeaddrinfo(pAddrInfoRes);


					// need socket to be non-blocking
#ifdef WIN32
u_long nNoBlock = 1;
ioctlsocket(ListenerSocket, FIONBIO, &nNoBlock);
#else
int Rsltz = fcntl(ListenerSocket, F_SETFL, fcntl(ListenerSocket,F_GETFL,0) | O_NONBLOCK);
if(Rsltz == -1)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseListener: Unable to make socket nonblocking");
	return(false);
	}
	
#endif
m_ListenerSock = ListenerSocket;
return(true);
}


// SetupFDSets 
// Set up the three FD sets used with select() with the sockets to be monitored for events
// returns:
// 0 if no sockets to be monitored
// otherwise on windows returns the total number of sockets to be monitored
// otherwise on linux returns the highest file descriptor plus 1 as required on linux in the select() call
int			// returns 0 if no sockets to be monitored with select(), on windows the total number of monitored sockets, on linux the highest socket file descriptor plus 1
CBKSRequester::SetupFDSets(fd_set& ReadFDs,			// select() read available socket descriptor set  
			fd_set& WriteFDs,						// select() write accepted socket descriptor set
			fd_set& ExceptFDs)						// select() exceptions descriptor set
{
int HiFD;
tsBKSType *pType;
tsTxdRxd *pTxdRxd;
tsBKSRegSessionEx *pSessionEx;
int TypeIdx;
int SessIdx;

HiFD = 0;
FD_ZERO(&ReadFDs);
FD_ZERO(&WriteFDs);
FD_ZERO(&ExceptFDs);

// Add the listener socket to the read and except FD sets
// is one.
#ifdef WIN32
if (m_ListenerSock != INVALID_SOCKET)
	{
	HiFD = 1;
#else
if (m_ListenerSock != -1) 
	{
	HiFD = m_ListenerSock + 1;
#endif
	FD_SET(m_ListenerSock, &ReadFDs);
	FD_SET(m_ListenerSock, &ExceptFDs);
	}

// Add the ControlPair 1 socket to the read and except FD sets
pTxdRxd = &m_Ctrl[1];
#ifdef WIN32
if (pTxdRxd->Socket != INVALID_SOCKET &&
#else
if (pTxdRxd->Socket != -1 &&
#endif
	(pTxdRxd->flgSelMonExcept || pTxdRxd->flgSelMonRead || pTxdRxd->flgSelMonWrite))
	{
#ifdef WIN32
	HiFD += 1;
#else
	if(HiFD <= pTxdRxd->Socket)
		HiFD = pTxdRxd->Socket + 1;
#endif
	if (pTxdRxd->flgSelMonRead)
		FD_SET(pTxdRxd->Socket, &ReadFDs);
	if (pTxdRxd->flgSelMonWrite)
		FD_SET(pTxdRxd->Socket, &WriteFDs);
	pTxdRxd->flgSelMonExcept = 1;
	FD_SET(pTxdRxd->Socket, &ExceptFDs);
	}

// Add the connections being established to the FD sets
tsBKSSessEstab *pSessEstabs;
if(m_NumSessEstabs)
	{
	if((pSessEstabs = m_pBKSSessEstabs)==NULL)
		{
		Reset(true);
		return(-1);
		}
	for (SessIdx = 0; SessIdx < cMaxSessEstab - 1; SessIdx += 1, pSessEstabs += 1)
		{
#ifdef WIN32
		if (pSessEstabs->TxdRxd.Socket != INVALID_SOCKET &&
#else
		if (pSessEstabs->TxdRxd.Socket != -1 &&
#endif
			pSessEstabs->SEState != eSESnone && (pSessEstabs->TxdRxd.flgSelMonExcept || pSessEstabs->TxdRxd.flgSelMonRead || pSessEstabs->TxdRxd.flgSelMonWrite))
			{
#ifdef WIN32
			HiFD += 1;
#else
			if (HiFD <= pSessEstabs->TxdRxd.Socket)
				HiFD = pSessEstabs->TxdRxd.Socket + 1;
#endif
			if (pSessEstabs->TxdRxd.flgSelMonRead == 1)
					FD_SET(pSessEstabs->TxdRxd.Socket, &ReadFDs);
			if (pSessEstabs->TxdRxd.flgSelMonWrite == 1)
					FD_SET(pSessEstabs->TxdRxd.Socket, &WriteFDs);
			pSessEstabs->TxdRxd.flgSelMonExcept = 1;
			FD_SET(pSessEstabs->TxdRxd.Socket, &ExceptFDs);
			}
		}
	}

pType = m_pBKSTypes;
for (TypeIdx = 0; TypeIdx < eBKSPTPlaceHolder-1; TypeIdx++, pType += 1)
	{
	if (pType->Detail.BKSPType == eBKSPTUndefined || pType->NumSessions == 0)
		continue;
	if ((pSessionEx = pType->pFirstSession) != NULL)
		{
		do
			{
#ifdef WIN32
			if(pSessionEx->TxdRxd.Socket != INVALID_SOCKET &&
#else
			if (pSessionEx->TxdRxd.Socket != -1 &&
#endif
				(pSessionEx->TxdRxd.flgSelMonExcept || pSessionEx->TxdRxd.flgSelMonRead || pSessionEx->TxdRxd.flgSelMonWrite))
				{
#ifdef WIN32
				HiFD += 1;
#else
				if (HiFD <= pSessionEx->TxdRxd.Socket)
					HiFD = pSessionEx->TxdRxd.Socket + 1;
#endif
				if(pSessionEx->TxdRxd.flgSelMonRead)
					FD_SET(pSessionEx->TxdRxd.Socket, &ReadFDs);
				if(pSessionEx->TxdRxd.flgSelMonWrite)
					FD_SET(pSessionEx->TxdRxd.Socket, &WriteFDs);
				pSessionEx->TxdRxd.flgSelMonExcept = 1;
				FD_SET(pSessionEx->TxdRxd.Socket, &ExceptFDs);
				}
			}
		while((pSessionEx = pSessionEx->pNext) != NULL);
		}
	}
return(HiFD);
}

bool
CBKSRequester::RxData(tsTxdRxd *pRxd)					// receiving session data
{
int RxdLen;
UINT32 ExpRxdLen;
tsBKSPacHdr *pRxdHdr;

// ensure socket is valid
if (pRxd == NULL)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: passed NULL pRxd");
	return(false);
	}

if (pRxd->flgErr == 1)			// no point in attempting to read from a socket already known to have error
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: SessionID: %u pRxd->flgErr previously set", pRxd->SessionID);
	return(false);
	}

#ifdef WIN32
if (pRxd->Socket == INVALID_SOCKET)
#else
if (pRxd->Socket == -1)
#endif
	{
	pRxd->flgErr = 1;
	pRxd->flgErrReason = 1;      // invalid socket 
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: SessionID: %u pRxd->Socket previously set as invalid", pRxd->SessionID);
	return(false);
	}

	// perhaps a previously received frame simply has not yet been processed ..
if (pRxd->flgRxCplt == 1)
	return(true);

	// perhaps a frame is already buffered but not yet characterised as being a complete frame
if (pRxd->TotRxd >= sizeof(tsBKSPacHdr))
	{
	pRxdHdr = (tsBKSPacHdr *)pRxd->pRxdBuff;
	// received frame identifier must never be equal to last sent frame identifier except if a keepalive
	if (pRxdHdr->FrameType != eBKSHdrKeepalive && pRxdHdr->RxFrameID == pRxd->TxFrameID)
		{
		pRxd->flgErr = 1;
		pRxd->flgErrReason = 2;      // inconsistency in frame identifier
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: SessionID: %u inconsistency in received frame identifier, peer received identifier %u same as next to send", pRxd->SessionID, pRxdHdr->RxFrameID);
		return(false);
		}
	if (pRxdHdr->FrameType <= eBKSHdrnone || pRxdHdr->FrameType >= eBKSHdrPlaceHolder)
		{
		pRxd->flgErr = 1;
		pRxd->flgErrReason = 3;      // inconsistency in frame type
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: SessionID: %u inconsistency in received frame type, FrameType = %u", pRxd->SessionID, pRxdHdr->FrameType);
		return(false);
		}
	if (pRxdHdr->FrameLen <= pRxd->TotRxd)	// accepting this frame for further processing
		{
		pRxd->RxdTxFrameID = pRxdHdr->TxFrameID;
		pRxd->RxdRxFrameID = pRxdHdr->RxFrameID;
		if(pRxdHdr->FrameType == eBKSHdrKeepalive) // if was a keep alive then note when received and slough the frame
			{
			pRxd->PacRxdAtSecs = time(NULL);
			pRxd->CurPacRxd = 0;
			if(pRxd->TotRxd > sizeof(tsBKSPacHdr))
				memmove(pRxd->pRxdBuff,&pRxd->pRxdBuff[sizeof(tsBKSPacHdr)], pRxd->TotRxd - sizeof(tsBKSPacHdr));
			pRxd->TotRxd -= sizeof(tsBKSPacHdr);
			}
		else
			{
			pRxd->CurPacRxd = pRxdHdr->FrameLen;
			pRxd->flgRxCplt = 1;
			}
		return(true);
		}
	}

	// if at least tsBKSPacHdr.FrameLen in receive buffer then the total expected header+payload size is known
	// if less then complete reading in the header only
if (pRxd->TotRxd < sizeof(UINT32))
	ExpRxdLen = sizeof(tsBKSPacHdr);
else
	{
	ExpRxdLen = *(UINT32 *)pRxd->pRxdBuff;
	if (ExpRxdLen < sizeof(tsBKSPacHdr))
		{
		pRxd->flgErr = 1;
		pRxd->flgErrReason = 4;      // inconsistency in frame sizing
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: SessionID: %u inconsistency in received frame sizing, ExpRxdLen = %u", pRxd->SessionID, ExpRxdLen);
		return(false);
		}
	}
ExpRxdLen -= pRxd->TotRxd;

	// double check that pRxdBuff has been sized to hold the ExpRxdLen
if ((ExpRxdLen + pRxd->TotRxd) > pRxd->AllocdRxdBuff)
	{
	pRxd->flgErr = 1;
	pRxd->flgErrReason = 5;      // potential for receive buffer overflow 
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: SessionID: %u inconsistency in received frame sizing, ExpRxdLen = %u, pRxd->AllocdRxdBuff = %u", pRxd->SessionID, ExpRxdLen + pRxd->TotRxd, pRxd->AllocdRxdBuff);
	return(false);
	}

	// now try receiving ...
RxdLen = recv(pRxd->Socket, (char *)&pRxd->pRxdBuff[pRxd->TotRxd], ExpRxdLen, 0);
if (RxdLen == 0)	// 0 if socket was closed by peer
	{
	pRxd->flgErr = 1;
	pRxd->flgErrReason = 6;      // socket closed by peer 
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: SessionID: %u socket closed by peer", pRxd->SessionID);
	return(false);
	}

#ifdef WIN32
if (RxdLen == SOCKET_ERROR)  // maybe the socket is simply blocking
#else
if (RxdLen == -1)
#endif
	{
	int err;
#ifdef WIN32
	err = WSAGetLastError();
	if (err == WSAEWOULDBLOCK || err == WSAEINPROGRESS)
		return(true);
#else
	err = errno;
	if (err == EWOULDBLOCK || err == EAGAIN || err == EINTR)
		return(true);
#endif
	pRxd->flgErr = 1;
	pRxd->flgErrReason = 7;      // general socket error 
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: SessionID: %u recv() returned error %d", pRxd->SessionID, err);
	return(false);
	}

	// received at least 1 byte, can now be characterised as being a complete frame?
pRxd->TotRxd += RxdLen;
if (pRxd->TotRxd >= sizeof(tsBKSPacHdr))
	{
	pRxdHdr = (tsBKSPacHdr *)pRxd->pRxdBuff;
	if (pRxdHdr->FrameType != eBKSHdrKeepalive && pRxdHdr->RxFrameID == pRxd->TxFrameID)
		{
		pRxd->flgErr = 1;
		pRxd->flgErrReason = 8;      // inconsistency in frame identifier
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: SessionID: %u inconsistency in received frame identifier, peer received identifier %u same as next to send", pRxd->SessionID, pRxdHdr->RxFrameID);
		return(false);
		}
	if (pRxdHdr->FrameType <= eBKSHdrnone || pRxdHdr->FrameType >= eBKSHdrPlaceHolder)
		{
		pRxd->flgErr = 1;
		pRxd->flgErrReason = 9;      // inconsistency in frame type
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: SessionID: %u inconsistency in received frame type, FrameType = %u", pRxd->SessionID, pRxdHdr->FrameType);
		return(false);
		}
	if (pRxdHdr->FrameLen <= pRxd->TotRxd)	// accepting this frame for further processing
		{
		pRxd->PacRxdAtSecs = time(NULL);
		pRxd->RxdTxFrameID = pRxdHdr->TxFrameID;
		pRxd->RxdRxFrameID = pRxdHdr->RxFrameID;
		if (pRxdHdr->FrameType == eBKSHdrKeepalive) // if was a keep alive then already noted when received, slough the frame
			{
			pRxd->CurPacRxd = 0;
			if (pRxd->TotRxd > sizeof(tsBKSPacHdr))
				memmove(pRxd->pRxdBuff, &pRxd->pRxdBuff[sizeof(tsBKSPacHdr)], pRxd->TotRxd - sizeof(tsBKSPacHdr));
			pRxd->TotRxd -= sizeof(tsBKSPacHdr);
			pRxd->CurPacRxd = 0;
			}
		else
			{
			pRxd->CurPacRxd = pRxdHdr->FrameLen;
			pRxd->flgRxCplt = 1;
			}
		return(true);
		}
	}
return(true);
}

bool
CBKSRequester::TxData(tsTxdRxd *pTxd)
{
	int ActTxLen;
	int ReqTxLen;

	// ensure socket is valid 
	if (pTxd == NULL)
	{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "TxData: passed NULL pTxd");
		return(false);
	}

	pTxd->flgSelMonWrite = 0;

	if (pTxd->flgErr == 1)			// no point in attempting to write a socket already known to have errors
	{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "TxData: SessionID: %u pTxd->flgErr previously set", pTxd->SessionID);
		return(false);
	}

#ifdef WIN32
	if (pTxd->Socket == INVALID_SOCKET)
#else
	if (pTxd->Socket == -1)
#endif
	{
		pTxd->flgErr = 1;
		pTxd->flgErrReason = 10;      // invalid socket 
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "TxData: SessionID: %u pTxd->Socket previously set as invalid", pTxd->SessionID);
		return(false);
	}

	if (pTxd->TotTxd == 0 || pTxd->TotTxd == pTxd->CurTxd)
	{
		pTxd->CurTxd = 0;
		pTxd->TotTxd = 0;
		pTxd->flgTxCplt = 1;
		
		return(true);
	}
	pTxd->flgTxCplt = 0;
	ReqTxLen = pTxd->TotTxd - pTxd->CurTxd;
	ActTxLen = send(pTxd->Socket, (char *)&pTxd->pTxdBuff[pTxd->CurTxd], ReqTxLen, 0);

#ifdef WIN32
	if (ActTxLen == SOCKET_ERROR)
#else
	if (ActTxLen == -1)
#endif
	{
		int err;
#ifdef WIN32
		err = WSAGetLastError();
		if (err == WSAEWOULDBLOCK || err == WSAEINPROGRESS)
			{
			pTxd->flgSelMonWrite = 1;
			return(true);
			}
#else
		err = errno;
		if (err == EWOULDBLOCK || err == EAGAIN || err == EINTR)
			{
			pTxd->flgSelMonWrite = 1;
			return(true);
			}
#endif
		pTxd->flgErr = 1;
		pTxd->flgErrReason = 11;      // general socket error 
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "TxData: SessionID: %u send() returned error %d", pTxd->SessionID, err);
		return(false);
	}

	pTxd->CurTxd += ActTxLen;
	if (pTxd->CurTxd == ReqTxLen)
		{
		pTxd->PacTxdAtSecs = time(NULL);
		pTxd->flgTxCplt = 1;
		pTxd->CurTxd = 0;
		pTxd->TotTxd = 0;
		}
	else
		pTxd->flgSelMonWrite = 1;
	return true;
}



// control messages are usually single bytes, except that if bit 7 is set then the message length is in bits 0..6
// with message immediately following the initial byte
int								// received message is this length, 0 if no outstanding messages, < 0 if errors
CBKSRequester::RcvCtrl(int BuffLen,			// available buffer into which the control message can be copied (must be at 128 bytes)
		UINT8 *pBuff)			// copy control message into this buffer
{
tsTxdRxd *pCtrl;
UINT8 *pCtrlMsg;
UINT32 ExpRxdLen;
int RxdLen;

if(BuffLen < 128 || pBuff == NULL)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "RcvCtrl: Invalid parameters");
	return(-1);
	}

pCtrl = &m_Ctrl[1];
pCtrl->flgRxCplt = 0;
m_Ctrl[1].CurPacRxd = 0;
if (pCtrl->flgErr == 1)			// no point in attempting to read from a socket already known to have errors
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "RcvCtrl: flgErr previously set");
	return(-1);
	}
#ifdef WIN32
if (pCtrl->Socket == INVALID_SOCKET)
#else
if (pCtrl->Socket == -1)
#endif
	{
	pCtrl->flgErr = 1;
	pCtrl->flgErrReason = 10;      // invalid socket 
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "RcvCtrl: Socket previously set as invalid");
	return(-1);
	}

	// perhaps a message is already buffered but not yet characterised as being a complete message
if (pCtrl->TotRxd >= sizeof(UINT8))
	{
	pCtrlMsg = (UINT8 *)pCtrl->pRxdBuff;
	if(*pCtrlMsg & 0x080)
		ExpRxdLen = *pCtrlMsg & 0x07f;
	else
		ExpRxdLen = 1;

	if (ExpRxdLen <= pCtrl->TotRxd)	// accepting this message for further processing
		{
		pCtrl->CurPacRxd = ExpRxdLen;
		pCtrl->flgRxCplt = 1;
		memcpy(pBuff,pCtrl->pRxdBuff,ExpRxdLen);
		UINT32 Diff;
		if ((Diff = (m_Ctrl[1].TotRxd - m_Ctrl[1].CurPacRxd)) > 0)
			memmove(m_Ctrl[1].pRxdBuff, &m_Ctrl[1].pRxdBuff[m_Ctrl[1].CurPacRxd], Diff);
		m_Ctrl[1].TotRxd = Diff;
		m_Ctrl[1].CurPacRxd = 0;
		return(ExpRxdLen);
		}
	}

ExpRxdLen = min(128,pCtrl->AllocdRxdBuff - pCtrl->TotRxd);

	// now try receiving ...
RxdLen = recv(pCtrl->Socket, (char *)&pCtrl->pRxdBuff[pCtrl->TotRxd], ExpRxdLen, 0);
if (RxdLen == 0)	// 0 if socket was closed by peer
	{
	pCtrl->flgErr = 1;
	pCtrl->flgErrReason = 6;      // socket closed by peer 
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "RcvCtrl: socket closed by peer");
	return(-1);
	}
#ifdef WIN32
if (RxdLen == SOCKET_ERROR)  // maybe the socket is simply blocking
#else
if (RxdLen == -1)
#endif
	{
	int err;
#ifdef WIN32
	err = WSAGetLastError();
	if (err == WSAEWOULDBLOCK || err == WSAEINPROGRESS)
		return(0);
#else
	err = errno;
	if (err == EWOULDBLOCK || err == EAGAIN || err == EINTR)
		return(0);
#endif
	pCtrl->flgErr = 1;
	pCtrl->flgErrReason = 7;      // general socket error 
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: RcvCtrl() returned error %d", err);
	return(-1);
	}

	// received at least 1 byte, can now be characterised as being a complete message?
pCtrl->TotRxd += RxdLen;
if (pCtrl->TotRxd >= sizeof(UINT8))
	{
	pCtrlMsg = (UINT8 *)pCtrl->pRxdBuff;
	if(*pCtrlMsg & 0x080)
		ExpRxdLen = *pCtrlMsg & 0x07f;
	else
		ExpRxdLen = 1;

	if (ExpRxdLen <= pCtrl->TotRxd)	// accepting this message for further processing
		{
		pCtrl->CurPacRxd = ExpRxdLen;
		pCtrl->flgRxCplt = 1;
		memcpy(pBuff,pCtrl->pRxdBuff,ExpRxdLen);
		UINT32 Diff;
		if ((Diff = (m_Ctrl[1].TotRxd - m_Ctrl[1].CurPacRxd)) > 0)
			memmove(m_Ctrl[1].pRxdBuff, &m_Ctrl[1].pRxdBuff[m_Ctrl[1].CurPacRxd], Diff);
		m_Ctrl[1].TotRxd = Diff;
		m_Ctrl[1].CurPacRxd = 0;
		return(ExpRxdLen);
		}
	}
return(0);
}



bool 
CBKSRequester::NotifyCtrl(UINT8 Msg) 
{
if(Msg & 0x080)			// reserved for multibyte messages
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "NotifyCtrl: Attempted to send a multibyte message");
	return(false);
	}
return(SendCtrlMsg(sizeof(Msg),&Msg));
}

bool
CBKSRequester::SendCtrlMsg(int Len,				// number of control message bytes to send ptd at by pCtrlMsg, can be 0
								UINT8 *pCtrlMsg)	// pCtrlMsg pts to control message bytes
{
tsTxdRxd *pCtrl;
int ActTxLen;
int ReqTxLen;

pCtrl = &m_Ctrl[0];

pCtrl->flgSelMonWrite = 0;

if (pCtrl->flgErr == 1)			// no point in attempting to write a socket already known to have errors
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "SendCtrlSeq: flgErr previously set");
	return(false);
	}

#ifdef WIN32
if (pCtrl->Socket == INVALID_SOCKET)
#else
if (pCtrl->Socket == -1)
#endif
	{
	pCtrl->flgErr = 1;
	pCtrl->flgErrReason = 10;      // invalid socket 
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "SendCtrlSeq: Socket previously set as invalid");
	return(false);
	}

if((pCtrl->TotTxd + Len) > pCtrl->AllocdTxdBuff)
	{
	pCtrl->flgErr = 1;
	pCtrl->flgErrReason = 12;      // tx buffer overflow
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "SendCtrlSeq: Requested control sequence length (%d) would overflow send buffer",Len);
	return(false);
	}

if(Len > 0 && pCtrlMsg != NULL)
	{
	memcpy(&pCtrl->pTxdBuff[pCtrl->TotTxd],pCtrlMsg,Len);
	pCtrl->TotTxd += Len;
	}
ReqTxLen = pCtrl->TotTxd - pCtrl->CurTxd;
if(ReqTxLen > 0)
	{
	pCtrl->flgTxCplt = 0;
	ActTxLen = send(pCtrl->Socket, (char *)&pCtrl->pTxdBuff[pCtrl->CurTxd], ReqTxLen, 0);
#ifdef WIN32
	if (ActTxLen == SOCKET_ERROR)
#else
	if (ActTxLen == -1)
#endif
		{
		int err;
#ifdef WIN32
		err = WSAGetLastError();
		if (err == WSAEWOULDBLOCK || err == WSAEINPROGRESS)
			{
			pCtrl->flgSelMonWrite = 1;
			return(true);
			}
#else
		err = errno;
		if (err == EWOULDBLOCK || err == EAGAIN || err == EINTR)
			{
			pCtrl->flgSelMonWrite = 1;
			return(true);
			}
#endif
		pCtrl->flgErr = 1;
		pCtrl->flgErrReason = 11;      // general socket error 
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "SendCtrlSeq: SessionID: send() returned error %d", pCtrl->SessionID, err);
		return(false);
	}

	pCtrl->CurTxd += ActTxLen;
	if (pCtrl->CurTxd == ReqTxLen)
		{
		pCtrl->flgTxCplt = 1;
		pCtrl->CurTxd = 0;
		pCtrl->TotTxd = 0;
		}
	else
		pCtrl->flgSelMonWrite = 1;
	return true;	
	}
pCtrl->CurTxd = 0;
pCtrl->TotTxd = 0;
return(true);
}

int 
CBKSRequester::AcceptConnections(void)	// start accepting connections
{
int SelectRslt;
int HiFDS;
time_t CurTimeSecs;
tsBKSType *pType;
tsBKSRegSessionEx *pSessionEx;
int TypeIdx;
struct timeval SelectTimeout;
fd_set ReadFDs, WriteFDs, ExceptFDs;
socket_t AcceptedSock;
struct sockaddr_storage SockPeerAddr;
char szPeerHost[100];
char szPeerService[100];
socklen_t SockPeerAddrSize;
int SessEstabIdx;
tsBKSSessEstab *pSessEstab;
bool bOK;

// check if requested to terminate
#ifdef WIN32
if(InterlockedCompareExchange(&m_Terminate,1,1)==1)
#else
if(__sync_val_compare_and_swap (&m_Terminate,1,1)==1)
#endif
	{
	Reset(false);
	return(0);		
	}

AcquireLock(true);

HiFDS = 0;

while (1)
	{
	// check if termination request active
#ifdef WIN32
	if(InterlockedCompareExchange(&m_Terminate,1,1)==1)
#else
	if(__sync_val_compare_and_swap (&m_Terminate,1,1)==1)
#endif
		{
		Reset(true);
		return(0);
		}

	// terminate any Session negotiations which are taking more than cMaxSessEstabTime
	CurTimeSecs = time(NULL);
	m_NumSessEstabs = 0;
	if((pSessEstab = m_pBKSSessEstabs)==NULL)
		{
		Reset(true);	
		return(false);	
		}
	for (SessEstabIdx = 0; SessEstabIdx < cMaxSessEstab; SessEstabIdx += 1, pSessEstab += 1)
		{
		if (pSessEstab->SEState != eSESnone)
			{
			if((CurTimeSecs - pSessEstab->StartSecs) > cMaxSessEstabTime)
				{
				gDiagnostics.DiagOut(eDLInfo, gszProcName, "AcceptConnections: timed out on session establishment for session: %u", pSessEstab->TxdRxd.SessionID);
				ResetSessEstab(pSessEstab);
				}
			else
				{
				m_NumSessEstabs += 1;
				
#ifndef WIN32
				if (HiFDS <= pSessEstab->TxdRxd.Socket)
					HiFDS = pSessEstab->TxdRxd.Socket + 1;
#endif
				}
			}
		}

	DeleteAllSessionsInState(eBKSPSRegisteredTerm);			// terminate and delete any sessions in state eBKSPSRegisteredTerm or with an invalid socket 

	// check if time for a keep alive on any of the established sessions
	pType = m_pBKSTypes;
	for (TypeIdx = 0; TypeIdx < eBKSPTPlaceHolder - 1; TypeIdx++, pType += 1)
		{
		if (pType->Detail.BKSPType == eBKSPTUndefined || pType->NumSessions == 0)
			continue;
		if ((pSessionEx = pType->pFirstSession) != NULL)
			{
			do {
				if(pSessionEx->Session.BKSPState == eBKSPSRegisteredActv &&  pSessionEx->TxdRxd.TotTxd == 0)
					{
					if(pSessionEx->TxdRxd.flgKeepAliveReq == 1 ||  (CurTimeSecs - pSessionEx->TxdRxd.PacTxdAtSecs) > pType->Detail.KeepaliveSecs/2)
						{
						// construct a keep alive packet and send
						pSessionEx->TxdRxd.flgKeepAliveReq = 0;
						tsBKSPacHdr *pHdr = (tsBKSPacHdr *)pSessionEx->TxdRxd.pTxdBuff;
						pHdr->FrameFlags = 0;
						pHdr->FrameLen = sizeof(tsBKSPacHdr);
						pHdr->FrameType = eBKSHdrKeepalive;
						pHdr->RxFrameID = pSessionEx->TxdRxd.RxdTxFrameID;
						if(pSessionEx->TxdRxd.TxFrameID == 1)			// note: keep alives do not have incrementing TxFrameIDs, instead using previously sent TxFrameID
							pHdr->TxFrameID = 0x07f;
						else
							pHdr->TxFrameID = pSessionEx->TxdRxd.TxFrameID-1;		
						pSessionEx->TxdRxd.TotTxd = sizeof(tsBKSPacHdr);
						if(!TxData(&pSessionEx->TxdRxd))
							{
							ShutdownConnection(&pSessionEx->TxdRxd.Socket);
							pSessionEx->Session.BKSPState = eBKSPSRegisteredTerm;
							if(pSessionEx->Session.SessionID != 0)
								{
								UnallocSessionID(pSessionEx->Session.SessionID);
								pSessionEx->Session.SessionID = 0;
								pSessionEx->TxdRxd.SessionID = 0;
								}
							}
						}
					}
				}
			while((pSessionEx = pSessionEx->pNext) != NULL);
			}
		}

	UINT32 NumPendingReqs;
	UINT32 NumPendingResps;
	UINT32 TotRespsAvail;

    // if any RMI worker threads have pending requests or wanting to check for responses then allow these access
	ReleaseLock(true);
	do
		{
#ifdef WIN32
		NumPendingReqs = InterlockedCompareExchange(&m_NumPendingReqs,0,0);
		NumPendingResps = InterlockedCompareExchange(&m_NumPendingResps,0,0);
		TotRespsAvail = InterlockedCompareExchange(&m_TotRespsAvail,0,0);
#else
		NumPendingReqs = __sync_val_compare_and_swap(&m_NumPendingReqs,0,0);
		NumPendingResps = __sync_val_compare_and_swap(&m_NumPendingResps,0,0);
		TotRespsAvail = __sync_val_compare_and_swap(&m_TotRespsAvail,0,0);
#endif
		if(NumPendingReqs == 0 && (TotRespsAvail == 0 || NumPendingResps == 0))
			break;
		CUtility::SleepMillisecs(10);
		}
	while(1);
	AcquireLock(true);

 	SendRequestFrames();
	SelectTimeout.tv_sec = 5;
	SelectTimeout.tv_usec = 0;
	HiFDS = SetupFDSets(ReadFDs, WriteFDs, ExceptFDs);
	ReleaseLock(true);
	SelectRslt = select(HiFDS, &ReadFDs, &WriteFDs, &ExceptFDs, &SelectTimeout);
	AcquireLock(true);
	if(SelectRslt == 0)			// 0 if simply timed out with no socket events occurring
		continue;
	if (SelectRslt > 0)
		{
		// at least one monitored socket event has occurred
        // event could be on the listener socket (new Session connection), session negotiation, or an established session
		if (FD_ISSET(m_ListenerSock, &ReadFDs))
			{
			SockPeerAddrSize = (int)sizeof(SockPeerAddr);
			AcceptedSock = accept(m_ListenerSock,(sockaddr*)&SockPeerAddr, &SockPeerAddrSize); // note that accept() errors are silently sloughed
#ifdef WIN32
			if (AcceptedSock != INVALID_SOCKET)
#else
			if (AcceptedSock != -1)

#endif
				{
									// report the connections peer address
				getnameinfo((sockaddr*)&SockPeerAddr, SockPeerAddrSize, szPeerHost,sizeof(szPeerHost), szPeerService,sizeof(szPeerService), 0);

				UINT32 LastSessionID;
				if((LastSessionID = AllocSessionID())==0)
					{
					// all sessions already committed, can't accept another
					gDiagnostics.DiagOut(eDLInfo, gszProcName, "All sessions already committed, terminating connection negotiation with service provider Session at : host '%s' port '%s'", szPeerHost, szPeerService);
					ShutdownConnection(&AcceptedSock);
					}
				else
					{
					gDiagnostics.DiagOut(eDLInfo, gszProcName, "Started connection (session %u) negotiation with service provider Session at : host '%s' port '%s'", LastSessionID, szPeerHost, szPeerService);

		
						// need socket to be non-blocking
#ifdef WIN32
					u_long nNoBlock = 1;
					ioctlsocket(AcceptedSock, FIONBIO, &nNoBlock);
#else
					int Rsltz = fcntl(AcceptedSock, F_SETFL, fcntl(AcceptedSock,F_GETFL,0) | O_NONBLOCK);
					if(Rsltz == -1)
						{
						gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseConnect: Unable to make socket nonblocking");
						ReleaseLock(true);
						Reset(true);
						return(-1);
						}
	
#endif
				// initiate negotiations for type of service provision by provider
					bool bStartSessEstab;
					bStartSessEstab = StartSessEstab(LastSessionID,AcceptedSock, &SockPeerAddr);
					}
				}
			}
		else if (FD_ISSET(m_ListenerSock, &ExceptFDs))   // can't ignore errors on the listening socket, these are fatal!
			{
			int err;
#ifdef WIN32
			int errlen = sizeof(err);
#else
			socklen_t errlen = sizeof(err);
#endif
			getsockopt(m_ListenerSock, SOL_SOCKET, SO_ERROR, (char*)&err, &errlen);

			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Error on listening socket: %d", err);
			Reset(true);
			return(-1);
			}

		// check for control socket event
		if(FD_ISSET(m_Ctrl[1].Socket,&ExceptFDs))
			{
			}
		if (m_Ctrl[1].flgSelMonRead && FD_ISSET(m_Ctrl[1].Socket, &ReadFDs))
			{
			FD_CLR(m_Ctrl[1].Socket, &ReadFDs);
			while (1)
				{     
				int CtrlMsgLen;
				UINT8 CtrlMsg[256];
                                                                                                                                                                                    \
				if((CtrlMsgLen = RcvCtrl(sizeof(CtrlMsg),CtrlMsg))==0)
					break;

					// handle this received control payload packet
				ProcessCtrlMsg(CtrlMsgLen,CtrlMsg);
				}
			}

		// check for any events on connections still being negotiated
		if (m_NumSessEstabs > 0)
			{
			if((pSessEstab = m_pBKSSessEstabs)==NULL)
				{
				Reset(true);
				return(false);
				}
			for (SessEstabIdx = 0; SessEstabIdx < cMaxSessEstab; SessEstabIdx += 1, pSessEstab += 1)
				{
				if (pSessEstab->SEState == eSESnone)
					continue;
				bOK = true;
				if (FD_ISSET(pSessEstab->TxdRxd.Socket, &ExceptFDs))		// if any exception then immediately terminate connection
					{
					gDiagnostics.DiagOut(eDLInfo, gszProcName, "AcceptConnections: terminating due to socket exception with session : %u", pSessEstab->TxdRxd.SessionID);
					ResetSessEstab(pSessEstab);
					m_NumSessEstabs -= 1;
					if (m_NumSessEstabs == 0)
						break;
					}
				else
					{
					if (FD_ISSET(pSessEstab->TxdRxd.Socket, &ReadFDs))
						{
						bOK = RxData(&pSessEstab->TxdRxd);
						FD_CLR(pSessEstab->TxdRxd.Socket, &ReadFDs);
						if(bOK && pSessEstab->TxdRxd.flgRxCplt == 1)
							{
							bool bProgressSessEstab;
							bProgressSessEstab = ProgressSessEstab(pSessEstab,false);
							}
						}

					if (m_NumSessEstabs == 0)
						break;

					if (FD_ISSET(pSessEstab->TxdRxd.Socket, &WriteFDs))
						{
						bOK = TxData(&pSessEstab->TxdRxd);
						FD_CLR(pSessEstab->TxdRxd.Socket, &WriteFDs);
						if (bOK && pSessEstab->TxdRxd.flgTxCplt == 1)
							{
							bool bProgressSessEstab;
							bProgressSessEstab = ProgressSessEstab(pSessEstab, true);
							}
						}
					if(!bOK)
						{
						gDiagnostics.DiagOut(eDLInfo, gszProcName, "AcceptConnections: terminating session establishment with session : %u", pSessEstab->TxdRxd.SessionID);
						ResetSessEstab(pSessEstab,true);
						m_NumSessEstabs -= 1;
						}
					if (m_NumSessEstabs == 0)
						break;
					}
				}
			}

			// any events on the session established sockets?
		pType = m_pBKSTypes;
		for (TypeIdx = 0; TypeIdx < eBKSPTPlaceHolder-1; TypeIdx++, pType += 1)
			{
			if (pType->Detail.BKSPType == eBKSPTUndefined || pType->NumSessions == 0)
				continue;
			if ((pSessionEx = pType->pFirstSession) != NULL)
				{
				do {
#ifdef WIN32
					if (pSessionEx->TxdRxd.Socket != INVALID_SOCKET &&
#else
					if (pSessionEx->TxdRxd.Socket != -1 &&
#endif
						(pSessionEx->TxdRxd.flgSelMonExcept || pSessionEx->TxdRxd.flgSelMonRead || pSessionEx->TxdRxd.flgSelMonWrite))
						{
						bOK = true;
						const char* pcErrorType = 0;

							// See if this socket's flag is set in any of the FD sets.
						if (pSessionEx->TxdRxd.flgSelMonExcept && FD_ISSET(pSessionEx->TxdRxd.Socket, &ExceptFDs))
							{
							bOK = false;
							pcErrorType = "General socket error";
							FD_CLR(pSessionEx->TxdRxd.Socket, &ExceptFDs);
							}
						else
							{
							if (pSessionEx->TxdRxd.flgSelMonRead && FD_ISSET(pSessionEx->TxdRxd.Socket, &ReadFDs))
								{
								FD_CLR(pSessionEx->TxdRxd.Socket, &ReadFDs);

								while((bOK = RxData(&pSessionEx->TxdRxd))==true)
									{                                                                                                                                                                                         \
									if(!pSessionEx->TxdRxd.flgRxCplt)
										break;
									// slough if a control packet received, the time it was received has already been noted by RxData
									if (pSessionEx->TxdRxd.CurPacRxd == sizeof(tsBKSPacHdr))
										{
										UINT32 Diff;
										if ((Diff = (pSessionEx->TxdRxd.TotRxd - pSessionEx->TxdRxd.CurPacRxd)) > 0)
											memmove(pSessionEx->TxdRxd.pRxdBuff, &pSessionEx->TxdRxd.pRxdBuff[pSessionEx->TxdRxd.CurPacRxd], Diff);
										pSessionEx->TxdRxd.TotRxd = Diff;
										pSessionEx->TxdRxd.CurPacRxd = 0;
										pSessionEx->TxdRxd.flgRxCplt = 0;
										}
									else
										{
										// handle this received payload packet
										if((ProcessResponseFrame(pSessionEx)) < 0)
											{
											bOK = false;
											gDiagnostics.DiagOut(eDLInfo, gszProcName, "AcceptConnections: Inconsistencies in received response, terminating session: %u", pSessionEx->TxdRxd.SessionID);
											break;
											}
										}
									}

								}
							if (bOK && (pSessionEx->TxdRxd.flgSelMonWrite && FD_ISSET(pSessionEx->TxdRxd.Socket, &WriteFDs)))
								{
								bOK = TxData(&pSessionEx->TxdRxd);
								FD_CLR(pSessionEx->TxdRxd.Socket, &WriteFDs);
								}
							}

						if (!bOK)
							{
							int err;
#ifdef WIN32
							err = WSAGetLastError();
#else
							err = errno;
#endif
							gDiagnostics.DiagOut(eDLInfo, gszProcName, "AcceptConnections: socket error %d, terminating established session: %u", err, pSessionEx->TxdRxd.SessionID);
							ShutdownConnection(&pSessionEx->TxdRxd.Socket);
							pSessionEx->Session.BKSPState = eBKSPSRegisteredTerm;
							UnallocSessionID(pSessionEx->TxdRxd.SessionID);
							pSessionEx->TxdRxd.SessionID = 0;
							pSessionEx->Session.SessionID = 0;
							}
						}
					}
				while((pSessionEx = pSessionEx->pNext) != NULL);
				}
			}
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: select() failed");
		ReleaseLock(true);
		Reset(true);
		return(-1);
		}
	}
}

const size_t cChkPtsReqAllocSizeIncr = 0x01fffff;	// realloc change point check points in a minimum of this sized increments
const size_t cChkPtsRespAllocSizeIncr = 0x01ffff;	// realloc change point check points in a minimum of this sized increments


UINT32			// identifies created and initialised check point requests list
CBKSRequester::CreateChkPtReqs(UINT32 SessionID,	// uniquely identifying this session between service requester and provider
				UINT32 TypeSessionID,		// used as an index within the type to reference this session
				UINT64 ClassInstanceID,		// check pointing this class instance
			    UINT32 ClassMethodID,		// identifies class method
				UINT32 ParamSize,			// instance specific parameter size
				UINT32 InDataSize,			// instance specific input data size
				UINT32 MaxOutDataSize,		// instance specific max expected output data for this method
				UINT8 *pData)				// request parameters followed by input data
{
UINT32 ChkPtID;
tsChkPtReqs *pChkPtReqs;
tsChkPtReq *pChkPtReq;
size_t MemUsed;
size_t ReallocTo;

if(m_ppChkPtReqs == NULL || m_NumChkPtReqs == m_MaxChkPtReqs)
	return(0);

MemUsed = sizeof(tsChkPtReqs) - 1 + InDataSize; 
ReallocTo = max(cChkPtsReqAllocSizeIncr,MemUsed + (cChkPtsReqAllocSizeIncr/2));

// iterate over all check point lists and locate 1st unused 
for(ChkPtID = 0; ChkPtID < m_MaxChkPtReqs; ChkPtID+=1)
	{
	pChkPtReqs = m_ppChkPtReqs[ChkPtID];

	if(pChkPtReqs == NULL || (pChkPtReqs->SessionID == 0 && pChkPtReqs->AllocMem < ReallocTo))
		{
		if(pChkPtReqs != NULL)
			{
#ifdef _WIN32
			free(pChkPtReqs);	
#else
			if(pChkPtReqs != MAP_FAILED)
				munmap(pChkPtReqs,pChkPtReqs->AllocMem);
#endif
			m_ppChkPtReqs[ChkPtID] = NULL;
			}

#ifdef _WIN32
		pChkPtReqs = (tsChkPtReqs *) malloc(ReallocTo);	// initial and perhaps the only allocation
#else
			// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
		pChkPtReqs = (tsChkPtReqs *) mmap(NULL,ReallocTo, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
		if(pChkPtReqs == MAP_FAILED)
			pChkPtReqs = NULL;
#endif
		if(pChkPtReqs == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName,"CreateChkPtReqs: Memory allocation of %lld bytes failed - %s",(INT64)ReallocTo,strerror(errno));
			return(eBSFerrMem);
			}
	    memset(pChkPtReqs,0,sizeof(tsChkPtReqs));
		pChkPtReqs->AllocMem = ReallocTo;
		m_ppChkPtReqs[ChkPtID] = pChkPtReqs;
		}
	if(pChkPtReqs->SessionID == 0)
		{
		pChkPtReqs->SessionID = SessionID;
		pChkPtReqs->TypeSessionID = TypeSessionID;
		pChkPtReqs->ClassInstanceID = ClassInstanceID;
		pChkPtReqs->NumChkPtReqs = 1;

		ReallocTo = max(cChkPtsRespAllocSizeIncr,MaxOutDataSize + cChkPtsRespAllocSizeIncr/2);

		if(pChkPtReqs->pRespData == NULL || pChkPtReqs->AllocRespData < (UINT32)ReallocTo)
			{
			if(pChkPtReqs->pRespData != NULL)
				{
#ifdef _WIN32
				free(pChkPtReqs->pRespData);	
#else
				if(pChkPtReqs->pRespData != MAP_FAILED)
					munmap(pChkPtReqs->pRespData,pChkPtReqs->AllocRespData);
#endif
				pChkPtReqs->pRespData = NULL;
				}
			pChkPtReqs->AllocRespData = 0;
			}

		if(pChkPtReqs->pRespData == NULL)
			{ 

#ifdef _WIN32
			pChkPtReqs->pRespData = (UINT8 *) malloc(ReallocTo);	// initial and perhaps the only allocation
#else
				// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
			pChkPtReqs->pRespData = (UINT8 *) mmap(NULL,ReallocTo, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
			if(pChkPtReqs->pRespData == MAP_FAILED)
				pChkPtReqs->pRespData = NULL;
#endif
			if(pChkPtReqs->pRespData == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName,"CreateChkPtReqs: Memory allocation of %lld bytes failed - %s",(INT64)ReallocTo,strerror(errno));
				return(eBSFerrMem);
				}
			memset(pChkPtReqs->pRespData,0,ReallocTo);
			pChkPtReqs->AllocRespData = (UINT32)ReallocTo;
			}

		pChkPtReqs->UsedRespData = 0;		
		pChkPtReq = &pChkPtReqs->ChkPtReqs[0];
		pChkPtReq->ClassMethodID = ClassMethodID;
		pChkPtReq->ParamSize = ParamSize;
		pChkPtReq->InDataSize = InDataSize;
		pChkPtReq->MaxOutDataSize = MaxOutDataSize;
		memcpy(pChkPtReq->Data,pData,InDataSize);
		pChkPtReqs->UsedMem = MemUsed;
		m_NumChkPtReqs += 1;
		return(ChkPtID+1);
		}
	}
return(0);
}

bool
CBKSRequester::DeleteChkPtReqs(UINT32 ChkPtID)  // identifies check point list to be deleted, memory is not dealloc'd instead is marked as reusable
{
tsChkPtReqs *pChkPtReqs;
tsChkPtReq *pChkPtReq;
if(m_ppChkPtReqs == NULL || ChkPtID == 0 || ChkPtID > m_MaxChkPtReqs)
	return(false);
pChkPtReqs = m_ppChkPtReqs[ChkPtID-1];
if(pChkPtReqs != NULL)
	{
	if(pChkPtReqs->SessionID != 0 && m_NumChkPtReqs > 0)
		m_NumChkPtReqs -= 1;
	pChkPtReqs->UsedMem = 0;
	pChkPtReqs->UsedRespData = 0;
	pChkPtReqs->SessionID = 0;
	pChkPtReqs->TypeSessionID = 0;
	pChkPtReqs->ClassInstanceID = 0;
	pChkPtReqs->NumChkPtReqs = 0;
	pChkPtReq = &pChkPtReqs->ChkPtReqs[0];
	pChkPtReq->ClassMethodID = 0;
	pChkPtReq->ParamSize = 0;
	pChkPtReq->InDataSize = 0;
	pChkPtReq->MaxOutDataSize = 0;
	}
return(true);
}

bool			// true if request was successfully check pointed		
CBKSRequester::AddChkPtReq(UINT32 ChkPtID,  // identifies check point list to extend with this request
				UINT32 SessionID,			// uniquely identifying this session between service requester and provider
			    UINT32 ClassMethodID,		// identifies class method
				UINT32 ParamSize,			// instance specific parameter size
				UINT32 InDataSize,			// instance specific input data size
				UINT32 MaxOutDataSize,		// instance specific max expected output data for this method
				UINT8 *pData)				// request parameters followed by input data
{
tsChkPtReqs *pChkPtReqs;
tsChkPtReqs *pReallocd;
tsChkPtReq *pChkPtReq;
UINT8 *pReallocdResp;
size_t AllocMem;
size_t MemUsed;
size_t ReallocTo;

if(m_ppChkPtReqs == NULL || ChkPtID == 0 || ChkPtID > m_MaxChkPtReqs)
	return(false);
pChkPtReqs = m_ppChkPtReqs[ChkPtID-1];
if(pChkPtReqs == NULL || pChkPtReqs->SessionID != SessionID)
	return(false);

// alloc more memory as may be required to hold this new request
AllocMem = pChkPtReqs->AllocMem;
MemUsed = pChkPtReqs->UsedMem + sizeof(tsChkPtReq) - 1 + InDataSize; 
if(MemUsed >= AllocMem)
	{
	ReallocTo = MemUsed + (cChkPtsReqAllocSizeIncr/2);
#ifdef _WIN32
	pReallocd = (tsChkPtReqs *)realloc(pChkPtReqs,ReallocTo);
#else
	pReallocd = (tsChkPtReqs *)mremap(pChkPtReqs,AllocMem,ReallocTo,MREMAP_MAYMOVE);
	if(pReallocd == MAP_FAILED)
		pReallocd = NULL;
#endif
	if(pReallocd == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName,"AddChkPtReq: Memory reallocation to %lld bytes failed - %s",(INT64)ReallocTo,strerror(errno));
		return(false);
		}
	pChkPtReqs = pReallocd;
	pChkPtReqs->AllocMem = ReallocTo;
	m_ppChkPtReqs[ChkPtID-1] = pChkPtReqs;
	}

if(MaxOutDataSize >= pChkPtReqs->AllocRespData)
	{
	ReallocTo = MaxOutDataSize + cChkPtsRespAllocSizeIncr/2;
	AllocMem = pChkPtReqs->AllocRespData;
#ifdef _WIN32
	pReallocdResp = (UINT8 *)realloc(pChkPtReqs->pRespData,ReallocTo);
#else
	pReallocdResp = (UINT8 *)mremap(pChkPtReqs->pRespData,AllocMem,ReallocTo,MREMAP_MAYMOVE);
	if(pReallocd == MAP_FAILED)
		pReallocd = NULL;
#endif
	if(pReallocdResp == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName,"AddChkPtReq: Memory reallocation to %lld bytes failed - %s",(INT64)ReallocTo,strerror(errno));
		return(false);
		}
	pChkPtReqs->AllocRespData = (UINT32)ReallocTo;
	pChkPtReqs->pRespData = pReallocdResp;
	}

pChkPtReq = (tsChkPtReq *)((UINT8 *)pChkPtReqs + pChkPtReqs->UsedMem);
pChkPtReq->ClassMethodID = ClassMethodID;
pChkPtReq->ParamSize = ParamSize;
pChkPtReq->InDataSize = InDataSize;
pChkPtReq->MaxOutDataSize = MaxOutDataSize;
if(InDataSize > 0)
	memcpy(pChkPtReq->Data,pData,InDataSize);
pChkPtReqs->UsedMem = MemUsed;
return(true);
}

bool			// true if response was successfully updated		
CBKSRequester::AddChkPtResp(UINT32 ChkPtID,  // identifies check point list to update with this response
				UINT32 SessionID,			// uniquely identifying this session between service requester and provider
				UINT32 JobRslt,				// job processing result as returned by service provider
				UINT32 RespSize,			// response size
				UINT8 *pData)				// response data
{
tsChkPtReqs *pChkPtReqs;
UINT8 *pReallocdResp;
size_t ReallocTo;
size_t AllocMem;
if(m_ppChkPtReqs == NULL || ChkPtID == 0 || ChkPtID > m_MaxChkPtReqs)
	return(false);
pChkPtReqs = m_ppChkPtReqs[ChkPtID-1];
if(pChkPtReqs == NULL || pChkPtReqs->SessionID != SessionID)
	return(false);
pChkPtReqs->JobRslt = JobRslt;
pChkPtReqs->UsedRespData = RespSize;
if(RespSize >= pChkPtReqs->AllocRespData)
	{
	ReallocTo = RespSize + cChkPtsRespAllocSizeIncr/2;
	AllocMem = pChkPtReqs->AllocRespData;
#ifdef _WIN32
	pReallocdResp = (UINT8 *)realloc(pChkPtReqs->pRespData,ReallocTo);
#else
	pReallocdResp = (UINT8 *)mremap(pChkPtReqs->pRespData,AllocMem,ReallocTo,MREMAP_MAYMOVE);
	if(pReallocdResp == MAP_FAILED)
		pReallocdResp = NULL;
#endif
	if(pReallocdResp == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName,"AddChkPtResp: Memory reallocation to %lld bytes failed - %s",(INT64)ReallocTo,strerror(errno));
		return(false);
		}
	pChkPtReqs->AllocRespData = (UINT32)ReallocTo;
	pChkPtReqs->pRespData = pReallocdResp;
	}

if(RespSize)
	memcpy(pChkPtReqs->pRespData,pData,RespSize);
return(true);
}

