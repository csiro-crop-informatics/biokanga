// Copyright 2016 CSIRO  ( http://www.csiro.au/ ) m_bRawDedupe
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

#include "./SSW.h"
#include "./BKScommon.h"
#include "./BKSProvider.h"
const char * cDfltServerPort = "43123";		// default server port to connect to if not user specified

CBKSProvider::CBKSProvider()
{
m_bCreatedMutexes = false;
m_bTermConnectionReq = false;
m_bNotifiedReqs = false;
m_NumInitTypes = 0;
m_MaxServInsts = 0;
m_MaxServMemGB = 0;
m_CASSerialise = 0;
m_NumPendCpltd = 0;
m_NumWorkerInsts = 0;
m_NumSWAlignReqs = 0;
m_NumClassInsts = 0;
m_MaxClassInsts = 0;
m_HiClassInstanceID = 0;							
memset(m_ClassInstances,0,sizeof(m_ClassInstances));
memset(&m_BKSConnection, 0, sizeof(m_BKSConnection));
#ifdef _WIN32
m_BKSConnection.TxdRxd.Socket = INVALID_SOCKET;
#else
m_BKSConnection.TxdRxd.Socket = -1;
#endif
memset(m_BKSTypes, 0, sizeof(m_BKSTypes));
memset(m_Ctrl, 0, sizeof(m_Ctrl));
#ifdef _WIN32
m_Ctrl[0].Socket = INVALID_SOCKET;
m_Ctrl[1].Socket = INVALID_SOCKET;
#else
m_Ctrl[0].Socket = -1;
m_Ctrl[1].Socket = -1;
#endif

Reset();
}


CBKSProvider::~CBKSProvider()
{
Reset();
}

int
CBKSProvider::CreateMutexes(void)
{
if (m_bCreatedMutexes)
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
m_bCreatedMutexes = true;
m_CASSerialise = 0;
return(eBSFSuccess);
}


int
CBKSProvider::DeleteMutexes(void)				// delete locks/mutexes used in access serialisations
{
	if (!m_bCreatedMutexes)
		return(eBSFSuccess);
#ifndef _WIN32
	pthread_rwlock_destroy(&m_hRWLock);
#endif
	m_bCreatedMutexes = false;
	return(eBSFSuccess);
}

// acquire lock used for serialised access by multiple concurrent reader threads (bExclusive == false), or serialised access by a single thread (bExclusive == true)
void
CBKSProvider::AcquireLock(bool bExclusive)
{
	if (!m_bCreatedMutexes)
		return;
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
}

// no longer needing acquired lock, let other threads gain access
void
CBKSProvider::ReleaseLock(bool bExclusive)
{
	if (!m_bCreatedMutexes)
		return;
#ifdef _WIN32
	if (bExclusive)
		ReleaseSRWLockExclusive(&m_hRWLock);
	else
		ReleaseSRWLockShared(&m_hRWLock);
#else
	pthread_rwlock_unlock(&m_hRWLock);
#endif
}

void
CBKSProvider::AcquireCASSerialise(bool bPriority)	// if bPriority true then backoff time is reduced relative to if bPriority is false, increasing the probability of acquiring the serialisation lock if there is contention
{
int SpinCnt = 1000;
int BackoffMS = 1;

#ifdef _WIN32
while(InterlockedCompareExchange(&m_CASSerialise,1,0)!=0)
	{
	if(SpinCnt -= 1)
		continue;
	CUtility::SleepMillisecs(BackoffMS);
	SpinCnt = 100;
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
	SpinCnt = 100;
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
CBKSProvider::ReleaseCASSerialise(void)
{
#ifdef _WIN32
InterlockedCompareExchange(&m_CASSerialise,0,1);
#else
__sync_val_compare_and_swap(&m_CASSerialise,1,0);
#endif
}



int				// 0 no service request job available, 1 job details are being returned, -1 job is available but MaxParamSize < required, -2 job available but MaxRequestData < required, -3 if session being terminated
CBKSProvider::GetJobToProcess(UINT32 *pInstanceID,		// MANDATORY:  service instance identifier, JobResponse() must use this identifier
	  						UINT64 *pClassInstanceID,	// MANDATORY: returned class instance to apply job to
							UINT32 *pClassMethodID,			// MANDATORY: returned class method to apply on class instance
							  UINT32 *pMaxParamsSize,	//  MANDATORY: on input, max sized parameter block accepted, on return then the actual size of the parameter block
							  UINT8 *pParams,			//  MANDATORY: returned parameter block
							  UINT32 *pMaxRequestData,	//  MANDATORY: on input the max sized request data block accepted, on return the actual size of the request data block
							  UINT8 *pRequestData)		//  MANDATORY: returned request data block
{
teBKSPProvState BKSPState;
UINT32 InstanceID;
tsReqResp *pInstance;
if(pInstanceID == NULL || pMaxParamsSize == NULL || pParams == NULL || pMaxRequestData == NULL || pRequestData == NULL)
	return(eBSFerrParams);
	// check if all threads requested to terminate
#ifdef WIN32
	if(InterlockedCompareExchange(&m_TermAllThreads,1,1)==1)
#else
	if(__sync_val_compare_and_swap (&m_TermAllThreads,1,1)==1)
#endif
		return(-1);
*pInstanceID = 0;
AcquireLock(true);

if((BKSPState = (teBKSPProvState)m_BKSConnection.BKSPState) != eBKSPSAcceptedServiceActv || m_BKSConnection.InstancesReqAvail == 0)
	{
	ReleaseLock(true);
	if(BKSPState != eBKSPSAcceptedServiceActv)
		return(-3);
	CUtility::SleepMillisecs(20);		// to allow time for a new job request to be sent to this service provider
	return(0);
	}

pInstance = (tsReqResp *)m_BKSConnection.pReqResp;
for(InstanceID = 1; InstanceID <= m_BKSConnection.NumInstances; InstanceID+=1)
	{
	if(pInstance->InstanceID == InstanceID && pInstance->flgReqAvail)
		{
		if(*pMaxParamsSize < pInstance->ParamSize)
			{
			ReleaseLock(true);
			return(-1);
			}

		if (*pMaxRequestData < pInstance->InDataSize)
			{
			ReleaseLock(true);
			return(-2);
			}

		if(pInstance->ParamSize > 0)
			memcpy(pParams,pInstance->Data,pInstance->ParamSize);
		if (pInstance->InDataSize > 0)
			memcpy(pRequestData, &pInstance->Data[pInstance->ParamSize], pInstance->InDataSize);
		*pClassInstanceID = pInstance->ClassInstanceID;
		*pClassMethodID = pInstance->ClassMethodID;
		pInstance->flgReqAvail = 0;
		pInstance->flgProc = 1;
		m_BKSConnection.InstancesProc += 1;
		m_BKSConnection.InstancesReqAvail -= 1;
		*pInstanceID = pInstance->InstanceIDEx;
		ReleaseLock(true);
		return(1);
		}
	pInstance = (tsReqResp *)((UINT8 *)pInstance + m_BKSConnection.ReqRespInstSize);
	}

ReleaseLock(true);
return(0);
}


int					// 0 if response accepted, -1 if job does not exist or parameterisation errors, -3 if session terminating
CBKSProvider::JobResponse(INT32 InstanceIDEx,	// service instance identifier returned by GetJobToProcess
				UINT64 ClassInstanceID,		// response is for this class instance
				UINT32 ProcRslt,			// service request processing result
				UINT32 ResponseSize,		// response data block size
				UINT8 *pResponseData)		// response data
{
teBKSPProvState BKSPState;
UINT32 InstanceID;
tsReqResp *pInstance;
InstanceID = InstanceIDEx & 0x0fff;

	// check if all threads requested to terminate
#ifdef WIN32
	if(InterlockedCompareExchange(&m_TermAllThreads,1,1)==1)
#else
	if(__sync_val_compare_and_swap (&m_TermAllThreads,1,1)==1)
#endif
		return(-1);

#ifdef WIN32
InterlockedIncrement(&m_NumPendCpltd);
#else
__sync_fetch_and_add(&m_NumPendCpltd,1);
#endif

AcquireLock(true);
#ifdef WIN32
	InterlockedDecrement(&m_NumPendCpltd);
#else
	__sync_fetch_and_sub(&m_NumPendCpltd,1);
#endif
if ((BKSPState = (teBKSPProvState)m_BKSConnection.BKSPState) != eBKSPSAcceptedServiceActv || 
	!m_BKSConnection.InstancesProc ||
	InstanceID == 0 || InstanceID > m_BKSConnection.NumInstances || (ResponseSize > 0 && pResponseData == NULL))
	{
	ReleaseLock(true);
	return(BKSPState == eBKSPSAcceptedServiceActv ? -3 : -1);
	}
pInstance = (tsReqResp *)((UINT8 *)m_BKSConnection.pReqResp + (m_BKSConnection.ReqRespInstSize * (InstanceID - 1)));
if(pInstance->InstanceID != InstanceID || pInstance->InstanceIDEx != InstanceIDEx || !pInstance->flgProc)
	{
	ReleaseLock(true);
	return(-1);
	}

pInstance->ClassInstanceID = ClassInstanceID;
if(ResponseSize > 0)
	memcpy(pInstance->Data, pResponseData, ResponseSize);
pInstance->OutDataSize = ResponseSize;
pInstance->JobRslt = ProcRslt;
pInstance->flgCpltd = 1;
m_BKSConnection.InstancesCpltd += 1;
m_BKSConnection.InstancesProc -= 1;
NotifyCtrl();
ReleaseLock(true);
return(0);
}

int
CBKSProvider::TerminateConnection(bool bFreeMem,		  // true if any allocated memory to be free'd
									bool bCloseCtrl,	  // close control socket pair
									bool bResetInitTypes) // reset all service types
{
int Idx;
UINT8 *pRxdBuff;
UINT8 *pTxdBuff;
UINT8 *pReqResp;
UINT32 AllocdTxdBuff;
UINT32 AllocdRxdBuff;
UINT32 AllocdReqResp;

tsTxdRxd *pTxdRxd;

m_bTermConnectionReq = true;		// flags that current Session connection is being terminated

#ifdef _WIN32
if(m_BKSConnection.TxdRxd.Socket != INVALID_SOCKET)
	{
	closesocket(m_BKSConnection.TxdRxd.Socket);
	m_BKSConnection.TxdRxd.Socket = INVALID_SOCKET;
	}
#else
	if (m_BKSConnection.TxdRxd.Socket != -1)
		{
		close(m_BKSConnection.TxdRxd.Socket);
		m_BKSConnection.TxdRxd.Socket = -1;
		}
#endif
if(bFreeMem)
	{
	if (m_BKSConnection.TxdRxd.pRxdBuff != NULL)
		free(m_BKSConnection.TxdRxd.pRxdBuff);
	if (m_BKSConnection.TxdRxd.pTxdBuff != NULL)
		free(m_BKSConnection.TxdRxd.pTxdBuff);
	if(m_BKSConnection.pReqResp != NULL)
		free(m_BKSConnection.pReqResp);
	}
else
	{
	pRxdBuff = m_BKSConnection.TxdRxd.pRxdBuff;
	pTxdBuff = m_BKSConnection.TxdRxd.pTxdBuff;
	pReqResp = m_BKSConnection.pReqResp;
	AllocdTxdBuff = m_BKSConnection.TxdRxd.AllocdTxdBuff;
	AllocdRxdBuff = m_BKSConnection.TxdRxd.AllocdRxdBuff;
	AllocdReqResp = m_BKSConnection.AllocdReqResp;
	}
memset(&m_BKSConnection,0,sizeof(m_BKSConnection));
#ifdef _WIN32
m_BKSConnection.TxdRxd.Socket = INVALID_SOCKET;
#else
m_BKSConnection.TxdRxd.Socket = -1;
#endif

m_HiClassInstanceID = 0;	
if(m_MaxClassInsts && m_NumClassInsts)
	{
	tsClassInstance *pClassInstance;
	pClassInstance = m_ClassInstances;
	for(Idx=0; Idx < (int)m_MaxClassInsts; Idx++, pClassInstance++)
		{
		if(pClassInstance->pClass != NULL)
			delete pClassInstance->pClass;
		}
	}						
memset(m_ClassInstances,0,sizeof(m_ClassInstances));
m_NumClassInsts = 0;

if(!bFreeMem)
	{
	m_BKSConnection.TxdRxd.pRxdBuff = pRxdBuff;
	m_BKSConnection.TxdRxd.pTxdBuff = pTxdBuff;
	m_BKSConnection.pReqResp = pReqResp;
	m_BKSConnection.TxdRxd.AllocdTxdBuff = AllocdTxdBuff;
	m_BKSConnection.TxdRxd.AllocdRxdBuff = AllocdRxdBuff;
	m_BKSConnection.AllocdReqResp = AllocdReqResp;
	}

if(bCloseCtrl)
	{
	pTxdRxd = &m_Ctrl[0];
	for (Idx = 0; Idx < 2; Idx += 1, pTxdRxd += 1)
		{
#ifdef _WIN32
		if (pTxdRxd->Socket != INVALID_SOCKET)
			closesocket(pTxdRxd->Socket);
		pTxdRxd->Socket = INVALID_SOCKET;
#else
		if (pTxdRxd->Socket != -1)
			close(pTxdRxd->Socket);
		pTxdRxd->Socket = -1;
#endif
		if (bFreeMem)
			{
			if (pTxdRxd->pRxdBuff != NULL)
				free(pTxdRxd->pRxdBuff);
			if (pTxdRxd->pTxdBuff != NULL)
				free(pTxdRxd->pTxdBuff);
			}
		else
			{
			pRxdBuff = pTxdRxd->pRxdBuff;
			pTxdBuff = pTxdRxd->pTxdBuff;
			AllocdTxdBuff = pTxdRxd->AllocdTxdBuff;
			AllocdRxdBuff = pTxdRxd->AllocdRxdBuff;
			}
		memset(pTxdRxd, 0, sizeof(tsTxdRxd));
#ifdef _WIN32
		pTxdRxd->Socket = INVALID_SOCKET;
#else
		pTxdRxd->Socket = -1;
#endif
		if (!bFreeMem)
			{
			pTxdRxd->pRxdBuff = pRxdBuff;
			pTxdRxd->pTxdBuff = pTxdBuff;
			pTxdRxd->AllocdTxdBuff = AllocdTxdBuff;
			pTxdRxd->AllocdRxdBuff = AllocdRxdBuff;
			}
		}
	}

if(bResetInitTypes)
	{
	m_NumInitTypes = 0;
	memset(m_BKSTypes, 0, sizeof(m_BKSTypes));
	m_MaxServInsts = 0;
	m_MaxClassInsts = 0;
	m_HiClassInstanceID = 0;
	m_MaxServMemGB = 0;
	}
m_CASSerialise = 0;
m_bNotifiedReqs = false;
m_bTermConnectionReq = false;
return(eBSFSuccess);
}

int
CBKSProvider::Reset(void)
{
TerminateConnection(true,true,true);

if (m_bCreatedMutexes)
	{
	DeleteMutexes();
	m_bCreatedMutexes = false;
	}
m_MaxServInsts = 0;
m_MaxServMemGB = 0;
m_bNotifiedReqs = false;

return(eBSFSuccess);
}

int
CBKSProvider::Initialise(int MaxConnWait,			// wait for at most this many minutes for connection
						int MaxServInsts,		// max number of service instances supported
						int MaxClassInsts,		// max number of class instances supported (may differ from MaxServInsts)
						int MaxServMemGB)		// max allocatable memory available for all service instances
{
	Reset();
	CreateMutexes();
	m_MaxConnWait = MaxConnWait;
	m_MaxServInsts = MaxServInsts;
	m_MaxClassInsts = MaxClassInsts;
	m_MaxServMemGB = MaxServMemGB;

#ifdef WIN32
	WSADATA WsaDat;
	if (WSAStartup(MAKEWORD(2, 2), &WsaDat) != 0)
		return -1;
#endif
	return(0);
}

// if at time of registration user did not specify the costing function for determining the
// number of instances to support then use this default costing function
UINT32
CostNumInstances(teBKSPType BKSPType,				// determine number of service instances to offer for this service type
				 UINT32 MaxServiceInsts,			// limited to support a maximum of this many service instances
				 UINT32 MaxQuerySeqLen,			// accepted query sequences can be up to this length
				 UINT32 MaxTargSeqLen,			// accepted target sequences can be up to this length
				 UINT32 MaxReqPayloadSize,		// request payloads from the service requester, including framing, can be up to this size (UINT8s),
				 UINT32 MaxRespPayloadSize,		// response payloads from the service provider to requester, including framing, can be up to this  size (UINT8s)
				 double ScalePayload,					// ScaledPayload = (MaxReqPayloadSize + MaxTargSeqLen) * ScalePayload
				 double ScaleTargLen,					// ScaledTargLen = MaxTargSeqLen * ScaleTargLen
				 double ScaleQueryLen,				// ScaledQueryLen = MaxTargSeqLen * ScaleQueryLen
				 double ScaleScaledTargQuery,			// ScaledScaledTargQuery = (ScaledTargLen * ScaledQueryLen) * ScaleScaledTargQuery
				 double AvailResources)				// OfferedInstances = AvailResources / (ScaledScaledTargQuery + ScalePayload)
{
	UINT32 MaxServiceInstances;

double Divisor;  // used to place a floor on divisor of 0.0001
Divisor =  ((ScalePayload * (double)(MaxReqPayloadSize + MaxRespPayloadSize)) +
														 (ScaleScaledTargQuery * (((double)MaxQuerySeqLen * ScaleQueryLen) * ((double)MaxTargSeqLen * ScaleTargLen))));
if(Divisor < 0.0001)
	Divisor = 0.0001;

	// in this implementation using same default costing function for all service types
	MaxServiceInstances = (UINT32)((AvailResources + 0.5) / Divisor);
	if (MaxServiceInstances > MaxServiceInsts)
		MaxServiceInstances = MaxServiceInsts;
	return(MaxServiceInstances);
}


int					// returns total number of registered service types or teBSFrsltCodes error code if any parameterisation errors or already registered type
CBKSProvider::RegServiceType(teBKSPType BKSPType,		// registering this service type
			   UINT32 ProviderVersion,					// service provider version 
			   UINT32 MaxServiceInsts,					// limited to support a maximum of this many service instances
			   UINT32 MaxQuerySeqLen,					// accepted query sequences can be up to this length
			   UINT32 MaxTargSeqLen,					// accepted target sequences can be up to this length
			   UINT32 MaxReqPayloadSize,				// request payloads from the service requester, including framing, can be up to this size (UINT8s),
			   UINT32 MaxRespPayloadSize,				// response payloads from the service provider to requester, including framing, can be up to this  size (UINT8s)
				double ScalePayload,					// ScaledPayload = (MaxReqPayloadSize + MaxTargSeqLen) * ScalePayload
				double ScaleTargLen,					// ScaledTargLen = MaxTargSeqLen * ScaleTargLen
				double ScaleQueryLen,					// ScaledQueryLen = MaxTargSeqLen * ScaleQueryLen
				double ScaleScaledTargQuery,			// ScaledScaledTargQuery = (ScaledTargLen * ScaledQueryLen) * ScaleScaledTargQuery;
				double AvailResources,					// OfferedInstances = AvailResources / (ScaledScaledTargQuery + ScalePayload)
				UINT32 (*  pfnCostNumInstances)(teBKSPType BKSPType,				// determine number of service instances to offer for this service type
						UINT32 MaxServiceInsts,			// limited to support a maximum of this many service instances
						UINT32 MaxQuerySeqLen,			// accepted query sequences can be up to this length
						UINT32 MaxTargSeqLen,			// accepted target sequences can be up to this length
						UINT32 MaxReqPayloadSize,		// request payloads from the service requester, including framing, can be up to this size (UINT8s),
						UINT32 MaxRespPayloadSize,		// response payloads from the service provider to requester, including framing, can be up to this  size (UINT8s)
						double ScalePayload,					// ScaledPayload = (MaxReqPayloadSize + MaxTargSeqLen) * ScalePayload
						double ScaleTargLen,					// ScaledTargLen = MaxTargSeqLen * ScaleTargLen
						double ScaleQueryLen,				// ScaledQueryLen = MaxTargSeqLen * ScaleQueryLen
						double ScaleScaledTargQuery,			// ScaledScaledTargQuery = (ScaledTargLen * ScaledQueryLen) * ScaleScaledTargQuery
						double AvailResources))				// OfferedInstances = AvailResources / (ScaledScaledTargQuery + ScalePayload)
{
tsBKSType *pType;
int NumInitTypes;

	// ensure parameter values are acceptable
if (BKSPType <= eBKSPTUndefined || BKSPType >= eBKSPTPlaceHolder ||
	MaxServiceInsts < 1 || MaxServiceInsts > cMaxServiceInsts ||
		MaxReqPayloadSize < 1 || MaxReqPayloadSize > cAbsMaxReqPayloadSize ||
		MaxRespPayloadSize < 1 || MaxRespPayloadSize > cAbsMaxRespPayloadSize ||
		MaxQuerySeqLen > cAbsMaxQuerySeqLen || MaxTargSeqLen > cAbsMaxTargSeqLen)
		return(eBSFerrParams);

if(ScalePayload < 0.001 || ScalePayload > 1000.0 ||
   ScaleTargLen < 0.001 || ScaleTargLen > 1000.0 ||
   ScaleQueryLen < 0.001 || ScaleQueryLen > 1000.0 ||
   ScaleScaledTargQuery < 0.001 || ScaleScaledTargQuery > 1000.0 ||
   AvailResources < 0.001 || AvailResources > 10000000000000.0)	
	return(eBSFerrParams);

AcquireLock(true);			// needing exclusive access
pType = &m_BKSTypes[BKSPType - 1];
if (pType->BKSPType != eBKSPTUndefined)				// not allowed to reinitialise if already initialised
	{
	ReleaseLock(true);
	return(eBSFerrInternal);
	}

if(MaxServiceInsts > (UINT32)m_MaxServInsts)
	MaxServiceInsts = (UINT32)m_MaxServInsts;
if(AvailResources > (double)m_MaxServMemGB * 0x03fffffff)
	AvailResources = (double)m_MaxServMemGB * 0x03fffffff;

memset(pType, 0, sizeof(tsBKSType));
pType->BKSPType = BKSPType;
pType->ProviderVersion = ProviderVersion;
pType->MaxServiceInsts = MaxServiceInsts;
pType->MaxQuerySeqLen = MaxQuerySeqLen;
pType->MaxTargSeqLen = MaxTargSeqLen;
pType->MaxReqPayloadSize = MaxReqPayloadSize;
pType->MaxRespPayloadSize = MaxRespPayloadSize;
pType->ScalePayload = ScalePayload;
pType->ScaleTargLen = ScaleTargLen;
pType->ScaleQueryLen = ScaleQueryLen;
pType->ScaleScaledTargQuery = ScaleScaledTargQuery;
pType->AvailResources = AvailResources;
if(pfnCostNumInstances == NULL)
	pType->pfnCostNumInstances = CostNumInstances;
else
	pType->pfnCostNumInstances = pfnCostNumInstances;
m_NumInitTypes += 1;
NumInitTypes = m_NumInitTypes;
ReleaseLock(true);
return(NumInitTypes);
}


// service provider registration process is in multiple phases
// a) eBKSPSWaitReqServices: In the initial phase, immediately following connection acceptance, then expecting a list of requested services from the server
// b) eBKSPSSendOfferedService: Next provider responds with the offered service it can provide and number of service instances supported if that service is accepted
// c) eBKSPSWaitAcceptService: Then waits for server to respond with acceptance/rejection of the offered service
// d) eBKSPSAcceptedServiceActv: If server has accepted offered service then the service becomes active
//
teBSFrsltCodes				// cBSFSuccess if no errors and registration process is continued, cBSFSocketErr if any errors and connection has been terminated, eBSFerrMem if unable to allocate memory 
CBKSProvider::ProcessSessEstab(bool bCpltdWrite)	// false if frame received from server, true if frame sent
{
UINT32 RegIdx;
UINT32 ReqIdx;
tsBKSPacHdr *pRxdHdr;
tsBKSOfferedService *pOfferedService;
tsBKSReqServices *pReqServices;
tsBKSAcceptService *pAcceptService;
tsServiceDetail *pReqService;
tsBKSType *pRegService;
UINT32 MaxServiceInstances;
UINT32 PriorityServiceInstances;
tsBKSReqServices PriorityReqService;
tsBKSType *pPriorityRegService;
UINT32 NumReqTypes;
UINT32 MemReq;
UINT32 Diff;

memset(&PriorityReqService,0,sizeof(tsBKSReqServices));
if (m_BKSConnection.BKSPState == eBKSPSUndefined ||
#ifdef WIN32
	m_BKSConnection.TxdRxd.Socket == INVALID_SOCKET)
#else
	m_BKSConnection.TxdRxd.Socket == -1)
#endif
	{
	TerminateConnection(true,false,false);
	return(cBSFSocketErr);
	}

if (!bCpltdWrite)				// if a request frame has been received then ensure it is for the expected session and of expected request type in the current state
	{
	pRxdHdr = (tsBKSPacHdr *)m_BKSConnection.TxdRxd.pRxdBuff;
	PriorityReqService.Hdr = *pRxdHdr;

	// only accepted frames whilst negotiating are either list of requested services (eBKSHdrReqServices) in state eBKSPSUndefined, or acceptance of offered service (eBKSHdrAcceptService) in state
	if (pRxdHdr->FrameType != eBKSHdrReqServices && pRxdHdr->FrameType != eBKSHdrAcceptService)
		{
		TerminateConnection(true, false, false);
		return(cBSFSocketErr);
		}
	if (pRxdHdr->FrameType == eBKSHdrReqServices)
		{
		pReqServices = (tsBKSReqServices *)pRxdHdr;
		MemReq = sizeof(tsBKSReqServices) + (sizeof(tsServiceDetail) * (pReqServices->NumTypes - 1));
		if(m_BKSConnection.BKSPState != eBKSPSWaitReqServices || m_BKSConnection.TxdRxd.CurPacRxd != MemReq)
			{
			TerminateConnection(true, false, false);
			return(cBSFSocketErr);
			}
		// payload contains list of requested services, iterate and see if any match one of the services this endpoint can provide
		if((NumReqTypes = pReqServices->NumTypes) == 0 || NumReqTypes > cMaxServiceTypes)
			{
			TerminateConnection(true, false, false);
			return(cBSFSocketErr);
			}

		PriorityServiceInstances = 0;
		pReqService = pReqServices->Details;
		for(ReqIdx = 0; ReqIdx < pReqServices->NumTypes; ReqIdx +=1, pReqService+=1)
			{
			if(pReqService->BKSPType == eBKSPTUndefined || pReqService->BKSPType >= eBKSPTPlaceHolder)
				continue;
			pRegService = &m_BKSTypes[0];
			for (RegIdx = 0; RegIdx < m_NumInitTypes; RegIdx += 1, pRegService += 1)
				{
				if(pReqService->BKSPType != pRegService->BKSPType)
					continue;
				if(pReqService->MaxProviderVersion < pRegService->ProviderVersion || pReqService->MinProviderVersion > pRegService->ProviderVersion)
					continue;
				if (pReqService->MaxQuerySeqLen > pRegService->MaxQuerySeqLen ||
					pReqService->MaxTargSeqLen > pRegService->MaxTargSeqLen ||
					pReqService->MaxReqPayloadSize > pRegService->MaxReqPayloadSize ||
					pReqService->MaxRespPayloadSize > pRegService->MaxRespPayloadSize)
					continue;
				// can supply this service, determine max number of service instances 
				MaxServiceInstances = pRegService->pfnCostNumInstances(pRegService->BKSPType, pRegService->MaxServiceInsts,
												 pReqService->MaxQuerySeqLen, pReqService->MaxTargSeqLen, pReqService->MaxReqPayloadSize, pReqService->MaxRespPayloadSize,
													   pRegService->ScalePayload, pRegService->ScaleTargLen, pRegService->ScaleQueryLen, pRegService->ScaleScaledTargQuery, pRegService->AvailResources);
				if(MaxServiceInstances == 0)
					continue;
				if(MaxServiceInstances > pReqService->MaxServiceInsts)
					MaxServiceInstances = pReqService->MaxServiceInsts;

				if(PriorityReqService.NumTypes == 0 || PriorityReqService.Details[0].Priority < pReqService->Priority)
					{
					PriorityServiceInstances = MaxServiceInstances;
					pPriorityRegService = pRegService;
					PriorityReqService.NumTypes = 1;
					PriorityReqService.Details[0] = *pReqService;
					}
				}
			}	
		UINT8 *pTmp;
		if(PriorityServiceInstances > 0)				// can offer service?
			{
			size_t AllocMem = (PriorityReqService.Details[0].MaxReqPayloadSize + PriorityReqService.Details[0].MaxRespPayloadSize) * PriorityServiceInstances * 11 / 10; // allowing 10% additional for overheads
			if (m_BKSConnection.pReqResp == NULL || m_BKSConnection.AllocdReqResp < AllocMem)
				{
				if(m_BKSConnection.pReqResp != NULL)
					pTmp = m_BKSConnection.pReqResp;
				else
					pTmp = NULL; 
				if((m_BKSConnection.pReqResp = (UINT8 *)malloc(AllocMem))==NULL)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessSessEstab: unable to allocate %lld bytes memory for pReqResp", AllocMem);
					if(pTmp != NULL)
						free(pTmp);
					TerminateConnection(true, true, true);
					return(eBSFerrMem);
					}
				if(pTmp != NULL)
					{
					memcpy(m_BKSConnection.pReqResp,pTmp, m_BKSConnection.AllocdReqResp);
					free(pTmp);
					}
				m_BKSConnection.AllocdReqResp = (UINT32)AllocMem;
				}

			AllocMem = PriorityReqService.Details[0].MaxReqPayloadSize * 15 / 10;			// providers receive requests, allowing 50% additional for overheads
			if (m_BKSConnection.TxdRxd.pRxdBuff == NULL || m_BKSConnection.TxdRxd.AllocdRxdBuff < AllocMem)
				{
				if (m_BKSConnection.TxdRxd.pRxdBuff != NULL)
					pTmp = m_BKSConnection.TxdRxd.pRxdBuff;
				else
					pTmp = NULL;
				if((m_BKSConnection.TxdRxd.pRxdBuff = (UINT8 *)malloc(AllocMem))==NULL)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessSessEstab: unable to allocate %lld bytes memory for pRxdBuff", AllocMem);
					if(pTmp != NULL)
						free(pTmp);
					TerminateConnection(true, true, true);
					return(eBSFerrMem);
					}
				if (pTmp != NULL)
					{
					memcpy(m_BKSConnection.TxdRxd.pRxdBuff, pTmp, m_BKSConnection.TxdRxd.AllocdRxdBuff);
					free(pTmp);
					}
				m_BKSConnection.TxdRxd.AllocdRxdBuff = (UINT32)AllocMem;
				}

			AllocMem = PriorityReqService.Details[0].MaxRespPayloadSize * 15 / 10;			// providers respond to requests, allowing 50% additional for overheads
			if (m_BKSConnection.TxdRxd.pTxdBuff == NULL || m_BKSConnection.TxdRxd.AllocdTxdBuff < AllocMem)
				{
				if (m_BKSConnection.TxdRxd.pTxdBuff != NULL)
					pTmp = m_BKSConnection.TxdRxd.pTxdBuff;
				else
					pTmp = NULL;
				if((m_BKSConnection.TxdRxd.pTxdBuff = (UINT8 *)malloc(AllocMem))==NULL)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "ProcessSessEstab: unable to allocate %lld bytes memory for pTxdBuff", AllocMem);
					if (pTmp != NULL)
						free(pTmp);
					TerminateConnection(true, true, true);
					return(eBSFerrMem);
					}
				if (pTmp != NULL)
					{
					memcpy(m_BKSConnection.TxdRxd.pTxdBuff, pTmp, m_BKSConnection.TxdRxd.AllocdTxdBuff);
					free(pTmp);
					}
				m_BKSConnection.TxdRxd.AllocdTxdBuff = (UINT32)AllocMem;
				}

			m_BKSConnection.InstancesBusy = 0;
			m_BKSConnection.InstancesReqAvail = 0;
			m_BKSConnection.InstancesProc = 0;
			m_BKSConnection.InstancesCpltd = 0;
			m_BKSConnection.TotNumRequests = 0;
			m_BKSConnection.BKSPType = PriorityReqService.Details[0].BKSPType;
			m_BKSConnection.BKSPState = eBKSPSSendOfferedService;
			m_BKSConnection.NumInstances = PriorityServiceInstances;
			m_MaxClassInsts = PriorityServiceInstances;
			m_BKSConnection.MaxReqPayloadSize = PriorityReqService.Details[0].MaxReqPayloadSize;
			m_BKSConnection.MaxRespPayloadSize = PriorityReqService.Details[0].MaxRespPayloadSize;
			m_BKSConnection.ReqRespInstSize = sizeof(tsReqResp) + max(m_BKSConnection.MaxReqPayloadSize, m_BKSConnection.MaxRespPayloadSize);
			m_BKSConnection.AllocdReqResp = m_BKSConnection.ReqRespInstSize * m_BKSConnection.NumInstances;
			if ((m_BKSConnection.pReqResp = (UINT8 *)malloc(m_BKSConnection.AllocdReqResp)) == NULL)
				{
				gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessSessEstab: Unable to allocate memory for %d service instances totaling %d bytes", m_BKSConnection.NumInstances, m_BKSConnection.AllocdReqResp);
				TerminateConnection(true, true, true);
				return(cBSFSocketErr);
				}
			memset(m_BKSConnection.pReqResp, 0, m_BKSConnection.AllocdReqResp);
			m_BKSConnection.KeepaliveSecs = PriorityReqService.Details[0].KeepaliveSecs;
			if(m_BKSConnection.KeepaliveSecs < cMinKeepaliveSecs)
				m_BKSConnection.KeepaliveSecs = cMinKeepaliveSecs;
			else
				if (m_BKSConnection.KeepaliveSecs > cMaxKeepaliveSecs)
					m_BKSConnection.KeepaliveSecs = cMaxKeepaliveSecs;
			m_BKSConnection.TxdRxd.SessionID = PriorityReqService.Hdr.SessionID;

			pOfferedService = (tsBKSOfferedService *)&m_BKSConnection.TxdRxd.pTxdBuff[m_BKSConnection.TxdRxd.TotTxd];
			pOfferedService->Hdr = PriorityReqService.Hdr;
			pOfferedService->Hdr.FrameFlags = 0;
			pOfferedService->Hdr.RxFrameID = m_BKSConnection.TxdRxd.RxdTxFrameID;
			pOfferedService->Hdr.TxFrameID = m_BKSConnection.TxdRxd.TxFrameID++;
			if (m_BKSConnection.TxdRxd.TxFrameID > 0x07f)
				m_BKSConnection.TxdRxd.TxFrameID = 1;
			pOfferedService->Hdr.FrameType = eBKSHdrOfferedService;
			pOfferedService->Hdr.FrameLen = sizeof(tsBKSOfferedService);
			pOfferedService->BKSPType = (teBKSPType)m_BKSConnection.BKSPType;
			pOfferedService->ProviderVersion = pPriorityRegService->ProviderVersion;
			pOfferedService->ServiceInsts = PriorityServiceInstances;
			pOfferedService->ClassInsts = PriorityServiceInstances;
			m_BKSConnection.TxdRxd.TotTxd += sizeof(tsBKSOfferedService);
			if((Diff = (m_BKSConnection.TxdRxd.TotRxd - m_BKSConnection.TxdRxd.CurPacRxd)) > 0)
				{
				memmove(m_BKSConnection.TxdRxd.pRxdBuff,&m_BKSConnection.TxdRxd.pRxdBuff[m_BKSConnection.TxdRxd.CurPacRxd], Diff);
				m_BKSConnection.TxdRxd.TotRxd = Diff;
				}
			else
				m_BKSConnection.TxdRxd.TotRxd = 0;
			m_BKSConnection.TxdRxd.CurPacRxd = 0;
			m_BKSConnection.TxdRxd.flgRxCplt = 0;

			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessSessEstab: sending response with offered service");

			// send this response 
			if (!TxData(&m_BKSConnection.TxdRxd))
				{
				TerminateConnection(true, false, false);;
				return(cBSFSocketErr);
				}
			return(eBSFSuccess);
			}

			// unable to offer any service so let server know
		pOfferedService = (tsBKSOfferedService *)&m_BKSConnection.TxdRxd.pTxdBuff[m_BKSConnection.TxdRxd.TotTxd];
		pOfferedService->Hdr = PriorityReqService.Hdr;
		pOfferedService->Hdr.FrameFlags = 0;
		pOfferedService->Hdr.RxFrameID = m_BKSConnection.TxdRxd.RxdTxFrameID;
		pOfferedService->Hdr.TxFrameID = m_BKSConnection.TxdRxd.TxFrameID++;
		if(m_BKSConnection.TxdRxd.TxFrameID > 0x07f)
			m_BKSConnection.TxdRxd.TxFrameID = 1;
		pOfferedService->Hdr.FrameType = eBKSHdrOfferedService;
		pOfferedService->Hdr.FrameLen = sizeof(tsBKSOfferedService);
		pOfferedService->BKSPType = eBKSPTUndefined;
		pOfferedService->ProviderVersion = 0;
		pOfferedService->ServiceInsts = 0;
		pOfferedService->ClassInsts = 0;
		m_BKSConnection.TxdRxd.TotTxd += sizeof(tsBKSOfferedService);
		if ((Diff = (m_BKSConnection.TxdRxd.TotRxd - m_BKSConnection.TxdRxd.CurPacRxd)) > 0)
			{
			memmove(m_BKSConnection.TxdRxd.pRxdBuff, &m_BKSConnection.TxdRxd.pRxdBuff[m_BKSConnection.TxdRxd.CurPacRxd], Diff);
			m_BKSConnection.TxdRxd.TotRxd = Diff;
			}
		else
			m_BKSConnection.TxdRxd.TotRxd = 0;
		m_BKSConnection.TxdRxd.CurPacRxd = 0;
		m_BKSConnection.TxdRxd.flgRxCplt = 0;

		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessSessEstab: sending response with no offered service");
		// send this response 
		if (!TxData(&m_BKSConnection.TxdRxd))
			{
			TerminateConnection(true, false, false);;
			return(cBSFSocketErr);
			}
		return(eBSFSuccess);
		}

	if (pRxdHdr->FrameType == eBKSHdrAcceptService || pRxdHdr->FrameType == eBKSHdrRejectService)
		{
		if((m_BKSConnection.TxdRxd.CurPacRxd != sizeof(tsBKSAcceptService)) ||
		   pRxdHdr->SessionID != m_BKSConnection.TxdRxd.SessionID)
			{
			TerminateConnection(true, false, false);
			return((teBSFrsltCodes)2);
			}
		pAcceptService = (tsBKSAcceptService *)pRxdHdr;
		if(pAcceptService->BKSPType != m_BKSConnection.BKSPType)
			{
			TerminateConnection(true, false, false);
			return(cBSFSocketErr);
			}

		if((Diff = (m_BKSConnection.TxdRxd.TotRxd - m_BKSConnection.TxdRxd.CurPacRxd)) > 0)
			{
			memmove(m_BKSConnection.TxdRxd.pRxdBuff, &m_BKSConnection.TxdRxd.pRxdBuff[m_BKSConnection.TxdRxd.CurPacRxd],Diff);
			m_BKSConnection.TxdRxd.TotRxd = Diff;
			}
		else
			m_BKSConnection.TxdRxd.TotRxd = 0;
		m_BKSConnection.TxdRxd.CurPacRxd = 0;
		m_BKSConnection.TxdRxd.flgRxCplt = 0;

		// server has accepted offered service, ready to process service requests
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ProcessSessEstab: server has accepted offered %d service instances",m_BKSConnection.NumInstances);
		m_BKSConnection.BKSPState = eBKSPSAcceptedServiceActv;
		return(eBSFSuccess);
		}
	}

// completed sending a frame
switch (m_BKSConnection.BKSPState)
	{
		case eBKSPSSendOfferedService:		// sent the offered service
			m_BKSConnection.TxdRxd.flgTxCplt = 0;
			m_BKSConnection.TxdRxd.CurTxd = 0;
			m_BKSConnection.TxdRxd.TotTxd = 0;
			m_BKSConnection.TxdRxd.flgSelMonExcept = 1;
			m_BKSConnection.TxdRxd.flgSelMonRead = 1;
			m_BKSConnection.TxdRxd.flgSelMonWrite = 0;
			m_BKSConnection.BKSPState = eBKSPSWaitAcceptService;			// server expected to respond with either eBKSHdrAcceptService or eBKSHdrRejectService in a tsBKSAcceptService frame 
			return(eBSFSuccess);

		default:						// if any other state then close session
			break;
	}

TerminateConnection(true, false, false);
return(cBSFSocketErr);
}

int							// number of received requests allocated to service instances
CBKSProvider::HandleServiceRequests(void) // handler for received requests for allocation of a service instance to process the service request
{
UINT32 InstanceID;
sBKSServReq *pServReq;
sBKSServResp *pServResp;
tsReqResp *pInstance;
UINT32 TxFrameSize;
UINT32 Diff;

pServReq = (sBKSServReq *)m_BKSConnection.TxdRxd.pRxdBuff;

if(m_BKSConnection.TxdRxd.flgRxCplt == 0 || pServReq->Hdr.FrameType != eBKSHdrReq)
	return(0);

	// can't handle if all service instances are already committed to processing previous requests
if (m_BKSConnection.InstancesBusy >= m_BKSConnection.NumInstances)
	{
		// let requester know all service instances are committed - no room at the Inn 
	TxFrameSize = sizeof(sBKSServResp);
	if((m_BKSConnection.TxdRxd.TotTxd + TxFrameSize) > m_BKSConnection.TxdRxd.AllocdTxdBuff)
		return(-1);
	pServResp = (sBKSServResp *)&m_BKSConnection.TxdRxd.pTxdBuff[m_BKSConnection.TxdRxd.TotTxd];
	pServResp->Hdr.SessionID = m_BKSConnection.TxdRxd.SessionID;
	pServResp->Hdr.FrameFlags = 0;
	pServResp->Hdr.FrameType = eBKSHdrResp;
	pServResp->Hdr.RxFrameID = m_BKSConnection.TxdRxd.RxdTxFrameID;
	pServResp->Hdr.TxFrameID = m_BKSConnection.TxdRxd.TxFrameID++;
	if(m_BKSConnection.TxdRxd.TxFrameID > 0x07f)
		m_BKSConnection.TxdRxd.TxFrameID = 1;
	pServResp->Hdr.FrameLen = TxFrameSize;
	pServResp->DataSize = 0;
	pServResp->JobIDEx = pServReq->JobIDEx;
	pServResp->ClassInstanceID = pServReq->ClassInstanceID;
	pServResp->ClassMethodID = pServReq->ClassMethodID;
	pServResp->JobRslt = 0;
	m_BKSConnection.TxdRxd.TotTxd += pServResp->Hdr.FrameLen;
	}
else  // not all service instances are committed so can accept this request
	{
	pInstance = (tsReqResp *)m_BKSConnection.pReqResp;
	for (InstanceID = 1; InstanceID <= m_BKSConnection.NumInstances; InstanceID += 1, pInstance = (tsReqResp *)((UINT8 *)pInstance + m_BKSConnection.ReqRespInstSize))
		{
		if (pInstance->InstanceID == 0)  // 0 if this instance available
			{
			memset(pInstance, 0,sizeof(tsReqResp));
			pInstance->InstanceID = InstanceID;
			m_BKSConnection.TotNumRequests += 1;
			if(m_BKSConnection.TotNumRequests == 0)
				m_BKSConnection.TotNumRequests = 1;
			pInstance->InstanceIDEx = (InstanceID & 0x0fff) | (m_BKSConnection.TotNumRequests << 12);
			pInstance->JobIDEx = pServReq->JobIDEx;
			pInstance->ClassInstanceID = pServReq->ClassInstanceID;
			pInstance->ClassMethodID = pServReq->ClassMethodID;
			pInstance->ParamSize = pServReq->ParamSize;
			pInstance->InDataSize = pServReq->DataSize;
			if ((pServReq->DataSize + pServReq->ParamSize) > 0)
				memcpy(pInstance->Data, pServReq->ParamData, pServReq->DataSize + pServReq->ParamSize);
			else
				pInstance->Data[0] = 0;
			pInstance->flgReqAvail = 1;
			m_BKSConnection.InstancesReqAvail += 1;
			m_BKSConnection.InstancesBusy += 1;
			m_BKSConnection.TxdRxd.flgKeepAliveReq = 1;
			break;
			}
		}
	}	

if ((Diff = (m_BKSConnection.TxdRxd.TotRxd - m_BKSConnection.TxdRxd.CurPacRxd)) > 0)
	{
	memmove(m_BKSConnection.TxdRxd.pRxdBuff, &m_BKSConnection.TxdRxd.pRxdBuff[m_BKSConnection.TxdRxd.CurPacRxd], Diff);
	m_BKSConnection.TxdRxd.TotRxd = Diff;
	}
else
	m_BKSConnection.TxdRxd.TotRxd = 0;
m_BKSConnection.TxdRxd.CurPacRxd = 0;
m_BKSConnection.TxdRxd.flgRxCplt = 0;
return(0);
}


int										// number of responses assembled into m_BKSConnection.TxdRxd.pTxdBuff ready to be sent back to requester 
CBKSProvider::HandleServiceResponses(void) // locate those service instances with responses ready to be sent to the requester and send these responses
{
int NumResponses;
UINT32 TxFrameID;
UINT32 InstanceID;
UINT32 TxFrameSize;
sBKSServResp *pServResp;
tsReqResp *pInstance;
tsReqResp *pNxtInstance;

NumResponses = 0;

TxFrameID = m_BKSConnection.TxdRxd.TxFrameID;						// throttle back on adding any new frames until service requester has caught up a little
if(m_BKSConnection.TxdRxd.RxdRxFrameID > m_BKSConnection.TxdRxd.TxFrameID)
	TxFrameID += 0x07f;
if((TxFrameID - m_BKSConnection.TxdRxd.RxdRxFrameID) > 4)
	return(NumResponses);

if(m_BKSConnection.InstancesCpltd > 0)
	{
	pInstance = (tsReqResp *)m_BKSConnection.pReqResp;
	for (InstanceID = 1; m_BKSConnection.InstancesCpltd > 0 && InstanceID <= m_BKSConnection.NumInstances; InstanceID += 1)
		{
		if(pInstance->InstanceID != 0 && pInstance->flgCpltd)
			{
			TxFrameID = m_BKSConnection.TxdRxd.TxFrameID;						// throttle back on adding any new frames until service requester has caught up a little
			if(m_BKSConnection.TxdRxd.RxdRxFrameID > m_BKSConnection.TxdRxd.TxFrameID)
				TxFrameID += 0x07f;
			if((TxFrameID - m_BKSConnection.TxdRxd.RxdRxFrameID) > 4)
				break;

			TxFrameSize = sizeof(sBKSServResp) - 1 + pInstance->OutDataSize;
			if((m_BKSConnection.TxdRxd.AllocdTxdBuff - m_BKSConnection.TxdRxd.TotTxd) < TxFrameSize)
				break;
			pServResp = (sBKSServResp *)&m_BKSConnection.TxdRxd.pTxdBuff[m_BKSConnection.TxdRxd.TotTxd];
			pServResp->JobIDEx = pInstance->JobIDEx;
			pServResp->ClassInstanceID = pInstance->ClassInstanceID;
			pServResp->ClassMethodID = pInstance->ClassMethodID;
			pServResp->JobRslt = pInstance->JobRslt;
			pServResp->DataSize = pInstance->OutDataSize;
			if (pInstance->OutDataSize > 0)
				memcpy(pServResp->Data, pInstance->Data, pInstance->OutDataSize);
			pServResp->Hdr.FrameLen = TxFrameSize;
			pServResp->Hdr.FrameFlags = 0;
			pServResp->Hdr.SessionID = m_BKSConnection.TxdRxd.SessionID;
			pServResp->Hdr.FrameType = eBKSHdrResp;
			pServResp->Hdr.RxFrameID = m_BKSConnection.TxdRxd.RxdTxFrameID;
			pServResp->Hdr.TxFrameID = m_BKSConnection.TxdRxd.TxFrameID++;
			if (m_BKSConnection.TxdRxd.TxFrameID > 0x07f)
				m_BKSConnection.TxdRxd.TxFrameID = 1;
			m_BKSConnection.TxdRxd.TotTxd += pServResp->Hdr.FrameLen;
			pNxtInstance = (tsReqResp *)((UINT8 *)pInstance + m_BKSConnection.ReqRespInstSize);
			memset(pInstance,0,sizeof(tsReqResp));
			pInstance = pNxtInstance;
			m_BKSConnection.InstancesCpltd -= 1;
			m_BKSConnection.InstancesBusy -= 1;
			NumResponses += 1;
			}
		else
			pInstance = (tsReqResp *)((UINT8 *)pInstance + m_BKSConnection.ReqRespInstSize);
		}
	}
if(NumResponses)
	TxData(&m_BKSConnection.TxdRxd);
return(NumResponses);
}


#ifdef WIN32
//using namespace std;


// List of Winsock error constants mapped to an interpretation string.
// Note that this list must remain sorted by the error constants'
// values, because we do a binary search on the list when looking up
// items.
typedef struct TAG_sErrorEntry
{
	int nID;
	const char* pcMessage;
} tsErrorEntry;

static tsErrorEntry gaErrorList[] = {
	{ 0,                  "No error" },
	{ WSAEINTR,           "Interrupted system call" },
	{ WSAEBADF,           "Bad file number" },
	{ WSAEACCES,          "Permission denied" },
	{ WSAEFAULT,          "Bad address" },
	{ WSAEINVAL,          "Invalid argument" },
	{ WSAEMFILE,          "Too many open sockets" },
	{ WSAEWOULDBLOCK,     "Operation would block" },
	{ WSAEINPROGRESS,     "Operation now in progress" },
	{ WSAEALREADY,        "Operation already in progress" },
	{ WSAENOTSOCK,        "Socket operation on non-socket" },
	{ WSAEDESTADDRREQ,    "Destination address required" },
	{ WSAEMSGSIZE,        "Message too long" },
	{ WSAEPROTOTYPE,      "Protocol wrong type for socket" },
	{ WSAENOPROTOOPT,     "Bad protocol option" },
	{ WSAEPROTONOSUPPORT, "Protocol not supported" },
	{ WSAESOCKTNOSUPPORT, "Socket type not supported" },
	{ WSAEOPNOTSUPP,      "Operation not supported on socket" },
	{ WSAEPFNOSUPPORT,    "Protocol family not supported" },
	{ WSAEAFNOSUPPORT,    "Address family not supported" },
	{ WSAEADDRINUSE,      "Address already in use" },
	{ WSAEADDRNOTAVAIL,   "Can't assign requested address" },
	{ WSAENETDOWN,        "Network is down" },
	{ WSAENETUNREACH,     "Network is unreachable" },
	{ WSAENETRESET,       "Net connection reset" },
	{ WSAECONNABORTED,    "Software caused connection abort" },
	{ WSAECONNRESET,      "Connection reset by peer" },
	{ WSAENOBUFS,         "No buffer space available" },
	{ WSAEISCONN,         "Socket is already connected" },
	{ WSAENOTCONN,        "Socket is not connected" },
	{ WSAESHUTDOWN,       "Can't send after socket shutdown" },
	{ WSAETOOMANYREFS,    "Too many references, can't splice" },
	{ WSAETIMEDOUT,       "Connection timed out" },
	{ WSAECONNREFUSED,    "Connection refused" },
	{ WSAELOOP,           "Too many levels of symbolic links" },
	{ WSAENAMETOOLONG,    "File name too long" },
	{ WSAEHOSTDOWN,       "Host is down" },
	{ WSAEHOSTUNREACH,    "No route to host" },
	{ WSAENOTEMPTY,       "Directory not empty" },
	{ WSAEPROCLIM,        "Too many processes" },
	{ WSAEUSERS,          "Too many users" },
	{ WSAEDQUOT,          "Disc quota exceeded" },
	{ WSAESTALE,          "Stale NFS file handle" },
	{ WSAEREMOTE,         "Too many levels of remote in path" },
	{ WSASYSNOTREADY,     "Network system is unavailable" },
	{ WSAVERNOTSUPPORTED, "Winsock version out of range" },
	{ WSANOTINITIALISED,  "WSAStartup not yet called" },
	{ WSAEDISCON,         "Graceful shutdown in progress" },
	{ WSAHOST_NOT_FOUND,  "Host not found" },
	{ WSANO_DATA,         "No host data of that type was found" },

};
const int kNumMessages = sizeof(gaErrorList) / sizeof(tsErrorEntry);


const char*
CBKSProvider::WSAGetLastErrorMessage(const char* pcMessagePrefix,
									  int nErrorID /* = 0 */)
{
	int Ofs;
	int Idx;
	tsErrorEntry *pErrEntry;

	// Build basic error string
	static char acErrorBuffer[256];
	Ofs = sprintf(acErrorBuffer, "%s:", pcMessagePrefix);

	// Tack appropriate canned message onto end of supplied message prefix
	pErrEntry = gaErrorList;
	for (Idx = 0; Idx < kNumMessages; Idx++, pErrEntry++)
		if (pErrEntry->nID == nErrorID)
		{
			Ofs += sprintf(&acErrorBuffer[Ofs], "%s (%d)", pErrEntry->pcMessage, pErrEntry->nID);
			break;
		}
	if (Idx == kNumMessages)
		Ofs += sprintf(&acErrorBuffer[Ofs], "Unkown error (%d)", pErrEntry->nID);

	// Finish error message off and return it.
	return acErrorBuffer;
}
#endif

// ShutdownConnection 
// Gracefully shuts the connection Socket down. 
// Returns true if we're successful, false otherwise.
bool		// even if false is returned the socket may have been shutdown, but there was an error whilst shutting it down
CBKSProvider::ShutdownConnection(socket_t *pSocket)
{
	char SloughBuffer[0x03fff];
	int RxdBytes;
	int TotalRxdBytes;
	bool bRslt;

	if (pSocket == NULL)
		return(false);

#ifdef WIN32
	if (*pSocket == INVALID_SOCKET)   // not an error if socket already closed
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

int		// -2 unable to initialise, -3 unable to register the service type, -1 if socket level errors, 0  if requested to terminate, 1 if connection timeout
CBKSProvider::Process(int MaxConnWait,			// wait for at most this many minutes for connection
					   int MaxServInsts,		// max number of service instances supported
						int MaxClassInsts,		// max number of class instances supported (may differ from MaxServInsts)
						int MaxServMemGB,		// max allocatable memory available for all service instances
						const char* pszHost,	// listening on this host/IP address; NULL to use first INET IP local to this machine
					   char *pszService)		// listening on this service/port; NULL to use default port 
{
int Rslt;
Reset();
if ((Rslt = Initialise(MaxConnWait,MaxServInsts,MaxClassInsts,MaxServMemGB)) < 0)
	{
	Reset();
	return(-2);
	}

if ((Rslt = RegServiceType(eBKSPTSmithWaterman)) < 0)
	{
	Reset();
	return(-3);
	}

Rslt = ConnectServer(MaxConnWait, pszHost, pszService);
Reset();
return(Rslt);
}

// control sockets are currently primarily used to ensure that select() is interruptible when there are
// sessions which require requests to be sent to peers for service provision
//
const int cCtrlbuffSize =2048;     // control socket rx/tx buffers are this size
bool 
CBKSProvider::InitialiseCtrlSocks(void) // initialise the control sockets in m_Ctrl[]
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
CBKSProvider::ProcessCtrlMsg(int MsgLen,UINT8 *pMsg)			// process a message received by control socket m_Ctrl[1] - note: currently these messages are simply discarded
{
// currrently not processing control payloads
m_Ctrl[1].flgRxCplt = 0;
m_bNotifiedReqs = false;
return(true);
}




// InitialiseConnect 
bool					// true if connecting socket initialised
CBKSProvider::InitialiseConnect(int MaxConnWait,			// wait for at most this many minutes for connection
								const char* pszHost,	// connecting to this server host/IP address; NULL to use first INET IP local to this machine
								  char *pszService)		// connecting to this server service/port; NULL to use default port 
{
UINT32 Now;
UINT32 Then;

struct addrinfo AddrInfoHints;
struct addrinfo *pAddrInfoRes;
struct addrinfo *pNxtAddrInfoRes;

char szHost[100];
char szService[100];

socket_t ConnectSocket;
int Rslt;

#ifdef WIN32
if (m_BKSConnection.TxdRxd.Socket != INVALID_SOCKET)
	{
	closesocket(m_BKSConnection.TxdRxd.Socket);
	m_BKSConnection.TxdRxd.Socket = INVALID_SOCKET;
	}
m_BKSConnection.TxdRxd.Socket = INVALID_SOCKET;
#else
if (m_BKSConnection.TxdRxd.Socket != -1)
	{
	close(m_BKSConnection.TxdRxd.Socket);
	m_BKSConnection.TxdRxd.Socket = -1;
	}
m_BKSConnection.TxdRxd.Socket = -1;
#endif


if(m_BKSConnection.TxdRxd.pRxdBuff == NULL)
	{
	if((m_BKSConnection.TxdRxd.pRxdBuff = (UINT8 *)malloc(cMinTxRxBuffSize))==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseConnect: Unable to allocate %d bytes memory", cMinTxRxBuffSize);
		return(false);
		}
	m_BKSConnection.TxdRxd.AllocdRxdBuff = cMinTxRxBuffSize;
	}
if (m_BKSConnection.TxdRxd.pTxdBuff == NULL)
	{
	if ((m_BKSConnection.TxdRxd.pTxdBuff = (UINT8 *)malloc(cMinTxRxBuffSize)) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseConnect: Unable to allocate %d bytes memory", cMinTxRxBuffSize);
		return(false);
		}
	m_BKSConnection.TxdRxd.AllocdTxdBuff = cMinTxRxBuffSize;
	}



if (pszService == NULL || pszService == '\0')	// use default port if not specified by caller
	pszService = (char *)cDfltServerPort;

Then = (UINT32)time(NULL);
MaxConnWait *= 60;			// time() returns secs

memset(&AddrInfoHints, 0, sizeof(struct addrinfo));
AddrInfoHints.ai_family = AF_INET; 		// Return IPv4 and IPv6 choices
AddrInfoHints.ai_socktype = SOCK_STREAM;
AddrInfoHints.ai_protocol = IPPROTO_TCP;

Rslt = getaddrinfo(pszHost == NULL || pszHost[0] == '\0' ? NULL : pszHost, pszService == NULL || pszService[0] == '\0' ? cDfltServerPort : pszService, &AddrInfoHints, &pAddrInfoRes);
if (Rslt != 0)		// errors if != 0
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseListener: getaddrinfo() failed to resolve host '%s' and service '%s'", pszHost, pszService);
	return(false);
	}


do {
	Rslt = -1;
	for (pNxtAddrInfoRes = pAddrInfoRes; pNxtAddrInfoRes != NULL; pNxtAddrInfoRes = pNxtAddrInfoRes->ai_next)
		{
		ConnectSocket = socket(pNxtAddrInfoRes->ai_family, pNxtAddrInfoRes->ai_socktype, pNxtAddrInfoRes->ai_protocol);
	#ifdef WIN32
		if (ConnectSocket == INVALID_SOCKET)
	#else
		if (ConnectSocket == -1)
	#endif
			continue;		// try next address

		if ((Rslt = connect(ConnectSocket, pNxtAddrInfoRes->ai_addr, (int)pNxtAddrInfoRes->ai_addrlen)) == 0)  // 0 if connected successfully
			break;

	#ifdef WIN32
		closesocket(ConnectSocket);
		ConnectSocket = INVALID_SOCKET;
	#else
		close(ConnectSocket);
		ConnectSocket = -1;
	#endif
		}
	if(Rslt == 0)
		break;

#ifdef WIN32
	Sleep(30000);
#else
	sleep(30);
#endif
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "InitialiseConnect: Retrying connection ...");
	}
while(Rslt != 0 && (((Now = (UINT32)time(NULL)) - Then) <= (UINT32)MaxConnWait));

	// report the connection address
#ifdef WIN32
if (ConnectSocket != INVALID_SOCKET)
#else
if (ConnectSocket != -1)
#endif
	{
	getnameinfo(pNxtAddrInfoRes->ai_addr, (int)pNxtAddrInfoRes->ai_addrlen, szHost, sizeof(szHost), szService, sizeof(szService), 0);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "InitialiseConnect: connected to host '%s' and service '%s' with protocol %d", szHost, szService, pNxtAddrInfoRes->ai_protocol);
	}

freeaddrinfo(pAddrInfoRes);

#ifdef WIN32
if (ConnectSocket == INVALID_SOCKET)
#else
if (ConnectSocket == -1)
#endif
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "InitialiseConnect: Unable to bind socket to host '%s' and service '%s'", pszHost, pszService);
	return(false);
	}

// need socket to be non-blocking

#ifdef WIN32
u_long nNoBlock = 1;
ioctlsocket(ConnectSocket, FIONBIO, &nNoBlock);
#else
int Rsltz = fcntl(ConnectSocket, F_SETFL, fcntl(ConnectSocket,F_GETFL,0) | O_NONBLOCK);
if(Rsltz == -1)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "InitialiseConnect: Unable to make socket nonblocking");
	return(false);
	}
	
#endif


m_BKSConnection.TxdRxd.Socket = ConnectSocket; 
m_BKSConnection.TxdRxd.flgSelMonExcept = 1;
m_BKSConnection.TxdRxd.flgSelMonRead = 1;
m_BKSConnection.TxdRxd.flgSelMonWrite = 0;
m_BKSConnection.TxdRxd.RxdRxFrameID = 0;
m_BKSConnection.TxdRxd.RxdTxFrameID = 0;
m_BKSConnection.TxdRxd.TxFrameID = 1;
m_BKSConnection.BKSPState = eBKSPSWaitReqServices;
return(true);
}


// Set up the three FD sets used with select() with the sockets to be monitored for events
// returns:
// 0 if no sockets to be monitored
// otherwise on windows returns the total number of sockets to be monitored
// otherwise on linux returns the highest file descriptor plus 1 as required on linux in the select() call
int			// returns 0 if no sockets to be monitored with select(), on windows the total number of monitored sockets, on linux the highest socket file descriptor plus 1
CBKSProvider::SetupFDSets(fd_set& ReadFDs,			// select() read available socket descriptor set  
						   fd_set& WriteFDs,		// select() write accepted socket descriptor set
						   fd_set& ExceptFDs)		// select() exceptions descriptor set
{
int HiFDs;
tsTxdRxd *pTxdRxd;
FD_ZERO(&ReadFDs);
FD_ZERO(&WriteFDs);
FD_ZERO(&ExceptFDs);
HiFDs = 0;

	// Add the connection socket to the read and except FD sets
#ifdef WIN32
if(m_BKSConnection.TxdRxd.Socket != INVALID_SOCKET)
	{
	HiFDs = 1;
#else
if (m_BKSConnection.TxdRxd.Socket != -1)
	{
	HiFDs = m_BKSConnection.TxdRxd.Socket + 1;
#endif
	FD_SET(m_BKSConnection.TxdRxd.Socket, &ReadFDs);
	if(m_BKSConnection.TxdRxd.flgSelMonWrite)
		FD_SET(m_BKSConnection.TxdRxd.Socket, &WriteFDs);
	FD_SET(m_BKSConnection.TxdRxd.Socket, &ExceptFDs);
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
	HiFDs += 1;
#else
	if(HiFDs <= pTxdRxd->Socket)
		HiFDs = pTxdRxd->Socket + 1;
#endif
	if (pTxdRxd->flgSelMonRead)
		FD_SET(pTxdRxd->Socket, &ReadFDs);
	if (pTxdRxd->flgSelMonWrite)
		FD_SET(pTxdRxd->Socket, &WriteFDs);
	pTxdRxd->flgSelMonExcept = 1;
	FD_SET(pTxdRxd->Socket, &ExceptFDs);
	}
return(HiFDs);
}

bool
CBKSProvider::RxData(tsTxdRxd *pRxd)					// receiving session data
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
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: SessionID: %u inconsistency in received frame for type %d, peer received identifier %u same as next to send", 
						pRxd->SessionID,pRxdHdr->FrameType, pRxdHdr->RxFrameID);
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
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "RxData: SessionID: %u inconsistency in received frame type %d, peer received identifier %u same as next to send", 
								pRxd->SessionID, pRxdHdr->FrameType, pRxdHdr->RxFrameID);
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
CBKSProvider::TxData(tsTxdRxd *pTxd)
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
CBKSProvider::RcvCtrl(int BuffLen,			// available buffer into which the control message can be copied (must be at 128 bytes)
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
CBKSProvider::NotifyCtrl(UINT8 Msg) 
{
if(Msg & 0x080)			// reserved for multibyte messages
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "NotifyCtrl: Attempted to send a multibyte message");
	return(false);
	}
return(SendCtrlMsg(sizeof(Msg),&Msg));
}

bool
CBKSProvider::SendCtrlMsg(int Len,				// number of control message bytes to send ptd at by pCtrlMsg, can be 0
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

int				// 1 if timed out attempting to connect to server, 0 if terminated because requested to terminate, -1 if socket level errors
CBKSProvider::ConnectServer(int MaxConnWait,			// wait for at most this many minutes for connection
							const char* pszHost,		// connect to this server host/IP address; NULL to use first INET IP local to this machine
							char *pszService)			// server expected to be listening on this service/port; NULL to use default port 
{
int SelectRslt;
UINT32 NumPendCpltd;
int HiFDs;
bool bRxDataRslt;
UINT32 Diff;
time_t CurTimeSecs;
time_t Then;
tsBKSType *pType;
struct timeval SelectTimeout;
teBKSPProvState PrevBKSPState;
fd_set ReadFDs, WriteFDs, ExceptFDs;

AcquireLock(true);
TerminateConnection(true,false,false);
m_bTermConnectionReq = false;

gDiagnostics.DiagOut(eDLInfo, gszProcName, "ConnectServer: Connecting to host name '%s' on service port '%s' allowing at most %d minutes for connection", 
									pszHost, pszService,MaxConnWait);
if (!InitialiseConnect(MaxConnWait,pszHost, pszService))
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ConnectServer: InitialiseConnect for host name '%s' on service port '%s' failed", (pszHost == NULL || pszHost[0] == '\0') ? "not specified" : pszHost,
							(pszService == NULL || pszService[0] == '\0') ? "not specified" : pszService);
	ReleaseLock(true);
	return(1);
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "ConnectServer: Creating service control socket pair"); 

if(!InitialiseCtrlSocks())
	{
#ifdef WIN32
	closesocket(m_BKSConnection.TxdRxd.Socket);
	m_BKSConnection.TxdRxd.Socket = INVALID_SOCKET;
#else
	close(m_BKSConnection.TxdRxd.Socket);
	m_BKSConnection.TxdRxd.Socket = -1;
#endif
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "ConnectServer: InitialiseCtrlSocks() failed");
	ReleaseLock(true);
	return(-1);
	}
gDiagnostics.DiagOut(eDLInfo, gszProcName, "ConnectServer: Service control socket pair created"); 

m_Ctrl[1].flgSelMonExcept = 1;
m_Ctrl[1].flgSelMonRead = 1;

PrevBKSPState = eBKSPSUndefined;
while (!m_bTermConnectionReq)		// keep processing for rxd/txd data over connection until connection is terminated
	{
	// if now providing service instances then start threads, one for each service instances
	if(PrevBKSPState != eBKSPSAcceptedServiceActv && m_BKSConnection.BKSPState == eBKSPSAcceptedServiceActv)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ConnectServer: Starting %d worker threads",m_BKSConnection.NumInstances);
		PrevBKSPState = eBKSPSAcceptedServiceActv;
		if(StartWorkerThreads(m_BKSConnection.NumInstances,m_BKSConnection.BKSPType)!=m_BKSConnection.NumInstances)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "ConnectServer: Unable to start %d worker threads",m_BKSConnection.NumInstances);
			m_BKSConnection.BKSPState = eBKSPSAcceptedServiceTerm;
			ReleaseLock(true);
			TerminateWorkerThreads();
			TerminateConnection(true, true, true);
			return(-1);
			}
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ConnectServer: Worker threads started, accepting service requests");
		Then = time(NULL);
		}

	if (m_BKSConnection.BKSPState == eBKSPSAcceptedServiceActv)   // if finished negotiating, check for any received requests or responses to be sent
		{
		HandleServiceRequests();
		HandleServiceResponses();
		}

	CurTimeSecs = time(NULL);

	if(m_BKSConnection.BKSPState == eBKSPSAcceptedServiceActv && (CurTimeSecs - Then) > 60)
		{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "ConnectServer: Processed total of %u requests (%u SW Align) with currently %u requests being progressed",m_BKSConnection.TotNumRequests,m_NumSWAlignReqs,m_BKSConnection.InstancesBusy);
		Then = CurTimeSecs;
		}

	if ((m_BKSConnection.BKSPState == eBKSPSAcceptedServiceActv) &&	m_BKSConnection.TxdRxd.TotTxd == 0)
		{
		if (m_BKSConnection.TxdRxd.flgKeepAliveReq == 1 || (CurTimeSecs - m_BKSConnection.TxdRxd.PacTxdAtSecs) > m_BKSConnection.KeepaliveSecs / 2)
			{
			m_BKSConnection.TxdRxd.flgKeepAliveReq = 0;
			// construct a keep alive packet and send
			tsBKSPacHdr *pHdr = (tsBKSPacHdr *)m_BKSConnection.TxdRxd.pTxdBuff;
			pHdr->FrameFlags = 0;
			pHdr->FrameLen = sizeof(tsBKSPacHdr);
			pHdr->FrameType = eBKSHdrKeepalive;
			pHdr->RxFrameID = m_BKSConnection.TxdRxd.RxdTxFrameID;
			if(m_BKSConnection.TxdRxd.TxFrameID == 1)		// note: keep alives do not have incrementing TxFrameIDs, instead using previously sent TxFrameID
				pHdr->TxFrameID = 0x07f;
			else
				pHdr->TxFrameID = m_BKSConnection.TxdRxd.TxFrameID - 1;
			m_BKSConnection.TxdRxd.TotTxd = sizeof(tsBKSPacHdr);
			TxData(&m_BKSConnection.TxdRxd);
			}
		}

	ReleaseLock(true);
	do {
#ifdef WIN32
		NumPendCpltd = InterlockedCompareExchange(&m_NumPendCpltd,0,0);
#else
		NumPendCpltd = __sync_val_compare_and_swap (&m_NumPendCpltd,0,0);
#endif
		if(NumPendCpltd != 0)
			CUtility::SleepMillisecs(50);
		}
	while(NumPendCpltd != 0);

	AcquireLock(true);
	SelectTimeout.tv_sec = 5;
	SelectTimeout.tv_usec = 0;
	HiFDs = SetupFDSets(ReadFDs, WriteFDs, ExceptFDs);
	ReleaseLock(true);	
	SelectRslt = select(HiFDs, &ReadFDs, &WriteFDs, &ExceptFDs, &SelectTimeout);
	AcquireLock(true);
	if (SelectRslt == 0)			// 0 if timed out with no socket events occurring
		continue;

	if(SelectRslt < 0)				// the select() call has failed? 
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: select() failed, terminating session");
		m_BKSConnection.BKSPState = eBKSPSAcceptedServiceTerm;
		ReleaseLock(true);
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating worker threads");
		TerminateWorkerThreads();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating connection with server");
		TerminateConnection(true, true, true);
		return(-1);
		}

		// select() is reporting a monitored event
		// double check for any non-recoverable errors on the connection, these are most likely errors in the packet framing
#ifdef WIN32
	if (m_BKSConnection.TxdRxd.Socket == INVALID_SOCKET ||
#else
	if (m_BKSConnection.TxdRxd.Socket == -1 ||
#endif
		m_BKSConnection.TxdRxd.flgErr)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: socket errors, terminating session");
			m_BKSConnection.BKSPState = eBKSPSAcceptedServiceTerm;
			ReleaseLock(true);
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating worker threads");
			TerminateWorkerThreads();
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating connection with server");
			TerminateConnection(true, true, true);
			return(-1);
			}

	// no connection errors

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

	pType = m_BKSTypes;
	if(m_BKSConnection.TxdRxd.flgSelMonExcept || m_BKSConnection.TxdRxd.flgSelMonRead || m_BKSConnection.TxdRxd.flgSelMonWrite)
		{
			// if an exception on the connected socket then terminate connection
		if (m_BKSConnection.TxdRxd.flgSelMonExcept && FD_ISSET(m_BKSConnection.TxdRxd.Socket, &ExceptFDs))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: socket errors, terminating session");
			m_BKSConnection.BKSPState = eBKSPSAcceptedServiceTerm;
			ReleaseLock(true);
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating worker threads");
			TerminateWorkerThreads();
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating connection with server");
			TerminateConnection(true, true, true);
			return(-1);
			}

		if (m_BKSConnection.TxdRxd.flgSelMonRead && FD_ISSET(m_BKSConnection.TxdRxd.Socket, &ReadFDs))
			{
			FD_CLR(m_BKSConnection.TxdRxd.Socket, &ReadFDs);
			bRxDataRslt = false;
			while((bRxDataRslt = RxData(&m_BKSConnection.TxdRxd))==true)    // get rxd data
				{
				if(m_BKSConnection.TxdRxd.flgRxCplt)		 // at least one complete frame received?
					{
					// slough if a keep alive received, the time it was received has already been noted by RxData
					if(m_BKSConnection.TxdRxd.CurPacRxd == sizeof(tsBKSPacHdr))
						{
						if ((Diff = (m_BKSConnection.TxdRxd.TotRxd - m_BKSConnection.TxdRxd.CurPacRxd)) > 0)
							{
							memmove(m_BKSConnection.TxdRxd.pRxdBuff, &m_BKSConnection.TxdRxd.pRxdBuff[m_BKSConnection.TxdRxd.CurPacRxd], Diff);
							m_BKSConnection.TxdRxd.TotRxd = Diff;
							}
						else
							m_BKSConnection.TxdRxd.TotRxd = 0;
						m_BKSConnection.TxdRxd.CurPacRxd = 0;
						m_BKSConnection.TxdRxd.flgRxCplt = 0;
						continue;
						}
					if(m_BKSConnection.BKSPState <= eBKSPSWaitAcceptService)	// still negotiating with server?
						ProcessSessEstab(false);
					else
						{
						tsBKSPacHdr *pHdr;
						pHdr = (tsBKSPacHdr *)m_BKSConnection.TxdRxd.pRxdBuff;
						switch(pHdr->FrameType) {
							case eBKSHdrReq:				// service request
								HandleServiceRequests();
								break;
							case eBKSHdrTermService:		// terminating this session and release all instance resources 
								gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: socket errors, terminating session");
								m_BKSConnection.BKSPState = eBKSPSAcceptedServiceTerm;
								ReleaseLock(true);
								gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating worker threads");
								TerminateWorkerThreads();
								gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating connection with server");
								TerminateConnection(true, true, true);
								return(0);
							default:						// any other frame type is treated as fatal error
								gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: socket errors, terminating session");
								m_BKSConnection.BKSPState = eBKSPSAcceptedServiceTerm;
								ReleaseLock(true);
								gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating worker threads");
								TerminateWorkerThreads();
								gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating connection with server");
								TerminateConnection(true, true, true);
								return(-1);
							}
						}
					}
				else
					break;
				}
			if(!bRxDataRslt)
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: socket errors, terminating session");
				m_BKSConnection.BKSPState = eBKSPSAcceptedServiceTerm;
				ReleaseLock(true);
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating worker threads");
				TerminateWorkerThreads();
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating connection with server");
				TerminateConnection(true, true, true);;
				return(-1);
				}

			}
		if ((m_BKSConnection.TxdRxd.flgSelMonWrite) && FD_ISSET(m_BKSConnection.TxdRxd.Socket, &WriteFDs))
			{
			if (!TxData(&m_BKSConnection.TxdRxd))
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: socket errors, terminating session");
				m_BKSConnection.BKSPState = eBKSPSAcceptedServiceTerm;
				ReleaseLock(true);
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating worker threads");
				TerminateWorkerThreads();
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating connection with server");
				TerminateConnection(true, true, true);
				return(-1);
				}

			FD_CLR(m_BKSConnection.TxdRxd.Socket, &WriteFDs);
			if (m_BKSConnection.TxdRxd.flgTxCplt)
				{
				if(m_BKSConnection.BKSPState <= eBKSPSWaitAcceptService)
					ProcessSessEstab(true);
				}
			}
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: socket errors, terminating session");
		m_BKSConnection.BKSPState = eBKSPSAcceptedServiceTerm;
		ReleaseLock(true);
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating worker threads");
		TerminateWorkerThreads();
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AcceptConnections: terminating connection with server");
		TerminateConnection(true, true, true);
		return(-1);
		}
	}
ReleaseLock(true);
return(-1);
}


#ifdef _WIN32
unsigned __stdcall WorkerInstance(void * pThreadPars)
#else
void *WorkerInstance(void * pThreadPars)
#endif
{
int Rslt;
tsWorkerInstance *pPars = (tsWorkerInstance *)pThreadPars;			// makes it easier not having to deal with casts!
CBKSProvider *pWorkerInstance = (CBKSProvider *)pPars->pThis;

Rslt = pWorkerInstance->ProcWorkerThread(pPars);
pPars->Rslt = Rslt;
#ifdef _WIN32
_endthreadex(0);
return(eBSFSuccess);
#else
pthread_exit(NULL);
#endif
}

int 
CBKSProvider::StartWorkerThreads(UINT32 NumInstances,	// this number of worker threads required
						   UINT8 BKSPType)		// workers are providing this service type
{
UINT32 MaxWait;
UINT32 Idx;
UINT32 ThreadIdx;
UINT32 StartedInstances;
tsWorkerInstance *pThreadPar;
m_TermAllThreads = 0;
m_NumPendCpltd = 0;
pThreadPar = m_WorkerInstances;
for(ThreadIdx = 0; ThreadIdx < NumInstances; ThreadIdx++,pThreadPar++)
	{
	memset(pThreadPar,0,sizeof(tsWorkerInstance));
#ifdef _WIN32
	pThreadPar->threadHandle = NULL;
#else
	pThreadPar->threadID = 0;
#endif
	if((pThreadPar->pReqData = (UINT8 *)malloc(cMaxReqDataSize))==NULL)
		break;
	if((pThreadPar->pParamData = (UINT8 *)malloc(cMaxReqParamSize))==NULL)
		break;
	if((pThreadPar->pRespData = (UINT8 *)malloc(cMaxRespDataSize))==NULL)
		break;
	if((pThreadPar->pszBuffer = (char *)malloc(cMaxMFABuffSize))==NULL)
		break;
	}
if(ThreadIdx != NumInstances)
	{
	pThreadPar = m_WorkerInstances;
	for(Idx = 0; Idx <= ThreadIdx; Idx++,pThreadPar++)
		{
		if(pThreadPar->pReqData != NULL)
			free(pThreadPar->pReqData);	
		if(pThreadPar->pParamData != NULL)
			free(pThreadPar->pParamData);	
		if(pThreadPar->pRespData != NULL)
			free(pThreadPar->pRespData);	
		if(pThreadPar->pszBuffer != NULL)
			free(pThreadPar->pszBuffer);	
		memset(pThreadPar,0,sizeof(tsWorkerInstance));
#ifdef _WIN32
		pThreadPar->threadHandle = NULL;
#else
		pThreadPar->threadID = 0;
#endif
		}
	return(0);
	}


m_ReqNumWorkerInsts = NumInstances;
m_NumWorkerInsts = 0;
m_NumSWAlignReqs = 0;
pThreadPar = m_WorkerInstances;
for (ThreadIdx = 1; ThreadIdx <= NumInstances; ThreadIdx++, pThreadPar++)
	{
	pThreadPar->ThreadIdx = ThreadIdx;
	pThreadPar->pThis = this;
#ifdef _WIN32
	pThreadPar->threadHandle = (HANDLE)_beginthreadex(NULL, 0x0fffff, WorkerInstance, pThreadPar, 0, &pThreadPar->threadID);
#else
	pThreadPar->threadRslt = pthread_create(&pThreadPar->threadID, NULL, WorkerInstance, pThreadPar);
#endif
	}

// allow threads a few seconds to startup
// check if all threads did actually startup; if they did so then m_NumWorkerInsts will have been incremented to NumInstances
MaxWait = 60;		// allowing at most 60 secs for threads to startup
do {
#ifdef WIN32
	Sleep(1000);
	StartedInstances = InterlockedCompareExchange(&m_NumWorkerInsts,NumInstances,NumInstances);
#else
	sleep(1);
	StartedInstances = __sync_val_compare_and_swap (&m_NumWorkerInsts,NumInstances,NumInstances);
#endif
	MaxWait -= 1;
	}
while(StartedInstances != NumInstances && MaxWait > 0);
if(StartedInstances != NumInstances)
	{
	TerminateWorkerThreads();
	StartedInstances = 0;
	}
return(StartedInstances);
}

// terminate all worker threads, allows a total of 2 minutes over all threads to self-terminate before force terminating
int     // number of threads requiring force termination
CBKSProvider::TerminateWorkerThreads(void)			
{
int NumForceTerminated;
UINT32 Idx;
UINT32 StartedInstances; 
tsWorkerInstance *pThreadPar;
time_t Then;
time_t Now;

#ifdef WIN32
StartedInstances = InterlockedCompareExchange(&m_NumWorkerInsts,0,0);
#else
StartedInstances = __sync_val_compare_and_swap (&m_NumWorkerInsts,0,0);
#endif
if(StartedInstances == 0)
	{
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: No worker threads to terminate");
	return(0);
	}


// request all worker threads to self terminate
#ifdef WIN32
InterlockedCompareExchange(&m_TermAllThreads,1,0);
#else
__sync_val_compare_and_swap (&m_TermAllThreads,0,1);
#endif
gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: Requesting %u worker threads to terminate",StartedInstances);
Then = time(NULL) + 120;
NumForceTerminated = 0;
pThreadPar = m_WorkerInstances;
for(Idx = 0; Idx < StartedInstances; Idx++, pThreadPar += 1)
	{
	Now = time(NULL);
	if(Now >= Then)
		Now = 1;
	else
		Now = Then - Now;

#ifdef WIN32
	if(pThreadPar->threadHandle != NULL)
		{
		if(WAIT_TIMEOUT == WaitForSingleObject(pThreadPar->threadHandle, (UINT32)Now * 1000))
			{
			NumForceTerminated += 1;
			TerminateThread(pThreadPar->threadHandle,0);
			}
		pThreadPar->threadHandle = NULL;
		}
#else
	if(pThreadPar->threadID != 0)
		{
		struct timespec ts;
		int JoinRlt;
		void *pExitRslt;
		clock_gettime(CLOCK_REALTIME, &ts);
		ts.tv_sec += Now;
		if ((JoinRlt = pthread_timedjoin_np(pThreadPar->threadID, &pExitRslt, &ts)) != 0)
			{
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "TerminateWorkerThreads: Force terminating thread %u, pthread_timedjoin_np() returned %d",pThreadPar->ThreadIdx,JoinRlt);
			NumForceTerminated += 1;
			pthread_cancel(pThreadPar->threadID);	
			pthread_join(pThreadPar->threadID, NULL);
			}
		pThreadPar->threadID = 0;
		}
#endif
	}

pThreadPar = m_WorkerInstances;
for(Idx = 0; Idx < StartedInstances; Idx++, pThreadPar += 1)
	{
	if(pThreadPar->pReqData != NULL)
		free(pThreadPar->pReqData);
	if(pThreadPar->pParamData != NULL)
		free(pThreadPar->pParamData);
	if(pThreadPar->pRespData != NULL)
		free(pThreadPar->pRespData);
	if(pThreadPar->pszBuffer != NULL)
		free(pThreadPar->pszBuffer);
	memset(pThreadPar,0,sizeof(tsWorkerInstance));
#ifdef _WIN32
	pThreadPar->threadHandle = NULL;
#else
	pThreadPar->threadID = 0;
#endif
	}

m_TermAllThreads = 0;	
return(NumForceTerminated);
}


tsClassInstance *
CBKSProvider::LocateClassInstance(UINT64 ClassInstanceID)
{
UINT32 Idx;
tsClassInstance *pInstance;
Idx = ((UINT32)(ClassInstanceID >> 40) & 0x0fff)-1;
if(Idx >= m_MaxClassInsts)
	return(NULL);
AcquireCASSerialise();
pInstance = &m_ClassInstances[Idx];
if(pInstance->ClassInstanceID == ClassInstanceID)
	{
	ReleaseCASSerialise();
	return(pInstance);
	}
ReleaseCASSerialise();
return(NULL);
}

tsClassInstance *
CBKSProvider::AllocClassInstance(void) // allocate a new tsClassInstance and initialise with class instance
{
UINT32 Idx;
tsClassInstance *pInstance;
pInstance = &m_ClassInstances[0];
AcquireCASSerialise();
for(Idx = 0; Idx < m_MaxClassInsts; Idx++, pInstance++)
	{
	if(pInstance->ClassInstanceID == 0)
		{
		m_HiClassInstanceID += 1;
		if(m_HiClassInstanceID > 0x0ffffffffff)
			m_HiClassInstanceID = 1;
		pInstance->ClassInstanceID = ((UINT64)m_BKSConnection.TxdRxd.SessionID << 53) | ((UINT64)(Idx+1) << 40) | m_HiClassInstanceID;
		m_NumClassInsts += 1;
		ReleaseCASSerialise();
		pInstance->LastAccessed = time(NULL);
		if(pInstance->pClass != NULL)		// shouldn't be required but best to be sure!
			{
			delete pInstance->pClass;
			pInstance->pClass = NULL;
			}
		pInstance->pClass = new CSSW;
		return(pInstance);
		}
	}
ReleaseCASSerialise();
return(NULL);
}

bool
CBKSProvider::FreeClassInstance(UINT64 ClassInstanceID) // free a previously allocated class instance
{
UINT32 Idx;
tsClassInstance *pInstance;
CSSW *pClass;
pClass = NULL;
Idx = ((UINT32)(ClassInstanceID >> 40) & 0x0fff)-1;
if(Idx >= m_MaxClassInsts)
	return(NULL);
AcquireCASSerialise();
pInstance = &m_ClassInstances[Idx];
if(pInstance->ClassInstanceID == ClassInstanceID)
	{
	pInstance->ClassInstanceID = 0;
	pInstance->LastAccessed = 0;
	if(pInstance->pClass != NULL)
		{
		pClass = pInstance->pClass;
		pInstance->pClass = NULL;
		}
	m_NumClassInsts -= 1;
	ReleaseCASSerialise();
	if(pClass != NULL)
		delete pClass;	
	return(true);
	}
ReleaseCASSerialise();
return(false);
}


int											// marshaled parameter required this many bytes
CBKSProvider::MarshalResp(UINT8 *pInto,					// marshal into this list
				teRMIParamType Type,			// parameter type
				void *pValue,				// parameter value
				UINT32 ValLen)				// length of parameter ptd to by pValue, only used if parameter type is pUint8

{
switch(Type) {
	case eRMIPTBool:	// boolean
		*pInto++ = (UINT8)eRMIPTBool;
		*pInto = *(bool *)pValue == true ? 1 : 0;
		return(sizeof(UINT8) + 1);

	case eRMIPTInt8:		// 8bit signed int
		*pInto++ = (UINT8)eRMIPTInt8;
		*pInto = *(INT8 *)pValue;
		return(sizeof(INT8) + 1);

	case eRMIPTUint8:       // 8bit  unsigned int
		*pInto++ = (UINT8)eRMIPTUint8;
		*pInto = *(UINT8 *)pValue;
		return(sizeof(UINT8) + 1);

	case eRMIPTInt32:		// 32bit signed int
		*pInto++ = (UINT8)eRMIPTInt32;
		*(INT32 *)pInto = *(INT32 *)pValue;
		return(sizeof(INT32) + 1);

	case eRMIPTUint32:		// 32bit unsigned int
		*pInto++ = (UINT8)eRMIPTInt32;
		*(UINT32 *)pInto = *(UINT32 *)pValue;
		return(sizeof(UINT32) + 1);

	case eRMIPTInt64:		// 64bit signed int
		*pInto++ = (UINT8)eRMIPTInt64;
		*(INT64 *)pInto = *(INT64 *)pValue;
		return(sizeof(INT64) + 1);

	case eRMIPTUint64:		// 64bit unsigned int
		*pInto++ = (UINT8)eRMIPTUint64;
		*(UINT64 *)pInto = *(UINT64 *)pValue;
		return(sizeof(UINT64) + 1);

	case eRMIPTDouble:		// floating point double
		*pInto++ = (UINT8)eRMIPTDouble;
		*(double *)pInto = *(double *)pValue;
		return(sizeof(double) + 1);

	case eRMIPTVarUint8:		// variable length
		*pInto++ = (UINT8)eRMIPTVarUint8;
		*(UINT32 *)pInto = ValLen;
		if(ValLen > 0 && pValue != NULL)
			{
			pInto += sizeof(UINT32);
			memcpy(pInto,pValue,ValLen);
			}
		return(sizeof(UINT32) + ValLen + 1);

	default:
		break;
	}
	
return(0);
}

int
CBKSProvider::UnmarshalReq(UINT32 DataLen,
				UINT8 *pFrom,		// unmarshal from this marshalled parameter list
				void *pValue)
{
UINT32 ValLen;
switch(*pFrom++) {
	case eRMIPTBool:	// boolean
		*(bool *)pValue = *pFrom == 0 ? false : true;
		return(sizeof(UINT8) + 1);

	case eRMIPTInt8:		// 8bit signed int
		*(INT8 *)pValue = *(INT8 *)pFrom;
		return(sizeof(INT8) + 1);

	case eRMIPTUint8:       // 8bit  unsigned int
		*(UINT8 *)pValue = *(UINT8 *)pFrom;
		return(sizeof(UINT8) + 1);

	case eRMIPTInt32:		// 32bit signed int
		*(INT32 *)pValue = *(INT32 *)pFrom;
		return(sizeof(INT32) + 1);

	case eRMIPTUint32:		// 32bit unsigned int
		*(UINT32 *)pValue = *(UINT32 *)pFrom;
		return(sizeof(UINT32) + 1);

	case eRMIPTInt64:		// 64bit signed int
		*(INT64 *)pValue = *(INT64 *)pFrom;
		return(sizeof(INT64) + 1);

	case eRMIPTUint64:		// 64bit unsigned int
		*(UINT64 *)pValue = *(UINT64 *)pFrom;
		return(sizeof(UINT64) + 1);

	case eRMIPTDouble:		// floating point double
		*(double *)pValue = *(double *)pFrom;
		return(sizeof(double) + 1);

	case eRMIPTVarUint8:		// variable length
		ValLen = *(UINT32 *)pFrom;
		if(ValLen > 0)
			{
			pFrom += sizeof(UINT32);
			*(void **)pValue = pFrom;
			return(sizeof(UINT32) + ValLen + 1);
			}
		*(void **)pValue = (void *)NULL;
		return(sizeof(UINT32) + 1);

	default:
		break;
	}
return(0);
}

int 
CBKSProvider::ProcWorkerThread(tsWorkerInstance *pThreadPar)
{
bool bRslt;
int iRslt;
int ReqDataOfs;
int RespDataOfs;
int JobRslt;
int NumJobsProc;
UINT32	MaxParamSize;
UINT32	MaxRequestData;
UINT32 InstanceID;

UINT64 ClassInstanceID;
UINT32 ClassMethodID;
tsClassInstance *pClassInstance;

NumJobsProc = 0;
#ifdef WIN32
InterlockedIncrement(&m_NumWorkerInsts);
#else
__sync_fetch_and_add(&m_NumWorkerInsts,1);
#endif

JobRslt = 0;
ClassInstanceID = 0;
while(JobRslt >= 0)
	{
	// check if all threads requested to terminate
#ifdef WIN32
	if(InterlockedCompareExchange(&m_TermAllThreads,1,1)==1)
#else
	if(__sync_val_compare_and_swap (&m_TermAllThreads,1,1)==1)
#endif
		{
		if(ClassInstanceID != 0)
			{
			FreeClassInstance(ClassInstanceID);
			ClassInstanceID = 0;
			}
		break;
		}
	MaxParamSize = cMaxReqParamSize;
	MaxRequestData = cMaxReqDataSize;
	if((JobRslt = GetJobToProcess(&InstanceID,&ClassInstanceID,&ClassMethodID, &MaxParamSize,pThreadPar->pParamData,&MaxRequestData,pThreadPar->pReqData)) > 0)
		{
		switch((teSWMethod)ClassMethodID) {
			case eSWMConstruct:				// instantiate new class instance
				if((pClassInstance = AllocClassInstance()) != NULL)
					ClassInstanceID = pClassInstance->ClassInstanceID;
				else
					ClassInstanceID = 0;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, 1,0, NULL);
				break;

			case eSWMDestruct:				// destroy class instance
				FreeClassInstance(ClassInstanceID);
				ClassInstanceID = 0;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, 1,0, NULL);
				break;

			case eSWMSetScores:				// SetScores()
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					int MatchScore;			// score for match
					int MismatchPenalty;	// penalty for mismatch
					int GapOpenPenalty;		// penalty for opening a gap
					int GapExtnPenalty;		// penalty if extending already opened gap
					int DlyGapExtn;			// delayed gap penalties, only apply gap extension penalty if gap at least this length
					int ProgPenaliseGapExtn;	// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles
					int AnchorLen;			// identified first and last anchors in alignment to be of at least this length
					ReqDataOfs = UnmarshalReq(sizeof(int),&pThreadPar->pReqData[0],&MatchScore);
					ReqDataOfs += UnmarshalReq(sizeof(int),&pThreadPar->pReqData[ReqDataOfs],&MismatchPenalty);
					ReqDataOfs += UnmarshalReq(sizeof(int),&pThreadPar->pReqData[ReqDataOfs],&GapOpenPenalty);
					ReqDataOfs += UnmarshalReq(sizeof(int),&pThreadPar->pReqData[ReqDataOfs],&GapExtnPenalty);
					ReqDataOfs += UnmarshalReq(sizeof(int),&pThreadPar->pReqData[ReqDataOfs],&DlyGapExtn);
					ReqDataOfs += UnmarshalReq(sizeof(int),&pThreadPar->pReqData[ReqDataOfs],&ProgPenaliseGapExtn);
					ReqDataOfs += UnmarshalReq(sizeof(int),&pThreadPar->pReqData[ReqDataOfs],&AnchorLen);
					bRslt = pClassInstance->pClass->SetScores(MatchScore,MismatchPenalty,GapOpenPenalty,GapExtnPenalty,DlyGapExtn,ProgPenaliseGapExtn,AnchorLen);
					}
				else
					bRslt = false;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, bRslt == true ? 1 : 0,0, NULL);
				break;

			case eSWMSetCPScores:			// SetCPScores()
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					int MatchScore;			// ClassifyPath() score for match
					int MismatchPenalty;	// ClassifyPath() penalty for mismatch
					int GapOpenPenalty;		// ClassifyPath() penalty for opening a gap
					int GapExtnPenalty;		// ClassifyPath() penalty if extending already opened gap
					ReqDataOfs = UnmarshalReq(sizeof(int),&pThreadPar->pReqData[0],&MatchScore);
					ReqDataOfs += UnmarshalReq(sizeof(int),&pThreadPar->pReqData[ReqDataOfs],&MismatchPenalty);
					ReqDataOfs += UnmarshalReq(sizeof(int),&pThreadPar->pReqData[ReqDataOfs],&GapOpenPenalty);
					ReqDataOfs += UnmarshalReq(sizeof(int),&pThreadPar->pReqData[ReqDataOfs],&GapExtnPenalty);
					bRslt = pClassInstance->pClass->SetCPScores(MatchScore,MismatchPenalty,GapOpenPenalty,GapExtnPenalty);
					}
				else
					bRslt = false;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, bRslt == true ? 1 : 0,0, NULL);	
			break;


			case eSWMSetMaxInitiatePathOfs:	// SetMaxInitiatePathOfs
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					int MaxInitiatePathOfs;	// require SW paths to have started within this many bp (0 to disable) on either the probe or target - effectively anchoring the SW 
					ReqDataOfs = UnmarshalReq(sizeof(int),&pThreadPar->pReqData[0],&MaxInitiatePathOfs);
					bRslt = pClassInstance->pClass->SetMaxInitiatePathOfs(MaxInitiatePathOfs);
					}
				else
					bRslt = false;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, bRslt == true ? 1 : 0,0, NULL);
				break;

			case eSWMPreAllocMaxTargLen:		// PreAllocMaxTargLen
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					UINT32 MaxTargLen;					// preallocate to process targets of this maximal length
					UINT32 MaxOverlapLen;			// allocating tracebacks for this maximal expected overlap, 0 if no tracebacks required
					ReqDataOfs = UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[0],&MaxTargLen);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&MaxOverlapLen);
					bRslt = pClassInstance->pClass->PreAllocMaxTargLen(MaxTargLen,MaxOverlapLen);
					}
				else
					bRslt = false;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, bRslt == true ? 1 : 0,0, NULL);
				break;

			case eSWMStartMultiAlignments:	// StartMultiAlignments
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					int SeqLen;						// probe sequence is this length
					etSeqBase *pProbeSeq;			// probe sequence 
					int Alignments;					// number of pairwise alignments to allocate for
					UINT8 Flags;					// flags
					ReqDataOfs = UnmarshalReq(sizeof(INT32),&pThreadPar->pReqData[0],&SeqLen);
					ReqDataOfs += UnmarshalReq(SeqLen,&pThreadPar->pReqData[ReqDataOfs],&pProbeSeq);
					ReqDataOfs += UnmarshalReq(sizeof(INT32),&pThreadPar->pReqData[ReqDataOfs],&Alignments);
					ReqDataOfs += UnmarshalReq(sizeof(UINT8),&pThreadPar->pReqData[ReqDataOfs],&Flags);

					iRslt = pClassInstance->pClass->StartMultiAlignments(SeqLen,pProbeSeq,Alignments,Flags);
					}
				else
					iRslt = -1;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, (UINT32)iRslt,0, NULL);
				break;

			case eSWMSetProbe:				// SetProbe
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					int SeqLen;						// probe sequence is this length
					etSeqBase *pProbeSeq;			// probe sequence 
					ReqDataOfs = UnmarshalReq(sizeof(INT32),&pThreadPar->pReqData[0],&SeqLen);
					ReqDataOfs += UnmarshalReq(SeqLen,&pThreadPar->pReqData[ReqDataOfs],&pProbeSeq);
					bRslt = pClassInstance->pClass->SetProbe(SeqLen,pProbeSeq);
					}
				else
					bRslt = false;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, bRslt == true ? 1 : 0,0, NULL);
				break;

			case eSWMSetTarg:				// SetTarg
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					int SeqLen;						// target sequence is this length
					etSeqBase *pTargSeq;			// target sequence 
					ReqDataOfs = UnmarshalReq(sizeof(INT32),&pThreadPar->pReqData[0],&SeqLen);
					ReqDataOfs += UnmarshalReq(SeqLen,&pThreadPar->pReqData[ReqDataOfs],&pTargSeq);
					bRslt = pClassInstance->pClass->SetTarg(SeqLen,pTargSeq);
					}
				else
					bRslt = false;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, bRslt == true ? 1 : 0,0, NULL);
				break;

			case eSWMSetAlignRange:			// SetAlignRange
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					UINT32 m_ProbeStartRelOfs;	// when aligning then start SW from this probe sequence relative offset
					UINT32 m_TargStartRelOfs; 	// and SW starting from this target sequence relative offset
					UINT32 m_ProbeRelLen;	// and SW with this probe relative length starting from m_ProbeStartRelOfs - if 0 then until end of probe sequence
					UINT32 m_TargRelLen;	// and SW with this target relative length starting from m_TargStartRelOfs - if 0 then until end of target sequence

					ReqDataOfs = UnmarshalReq(sizeof(UINT32),pThreadPar->pReqData,&m_ProbeStartRelOfs);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&m_TargStartRelOfs);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&m_ProbeRelLen);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&m_TargRelLen);
					iRslt = pClassInstance->pClass->SetAlignRange(m_ProbeStartRelOfs,m_TargStartRelOfs,m_ProbeRelLen,m_TargRelLen);
					}
				else
					iRslt = -1;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, (UINT32)iRslt,0, NULL);
				break;

			case eSWMAlign:					// Align
				tsSSWCell PeakScoreCell;
				tsSSWCell *pPeakMatchesCell;
				bool bPeakScoreCell; 
#ifdef WIN32
				InterlockedIncrement(&m_NumSWAlignReqs);		// users are likely to be interested in the number of SW alignments requested
#else
				__sync_fetch_and_add(&m_NumSWAlignReqs,1);
#endif
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					int ReqDataOfs;
					UINT32 MaxOverlapLen;			// process tracebacks for this maximal expected overlap, 0 if no tracebacks required
					ReqDataOfs = UnmarshalReq(sizeof(bool),pThreadPar->pReqData,&bPeakScoreCell);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&MaxOverlapLen);
					pPeakMatchesCell = pClassInstance->pClass->Align(bPeakScoreCell ? &PeakScoreCell : NULL,MaxOverlapLen);
					iRslt = 0;
					}
				else
					{
					bPeakScoreCell = false;
					pPeakMatchesCell = NULL;
					iRslt = -1;
					}
				if(pPeakMatchesCell != NULL)
					{
					RespDataOfs = MarshalResp(pThreadPar->pRespData,eRMIPTVarUint8,pPeakMatchesCell,sizeof(tsSSWCell));
					if(bPeakScoreCell)
						RespDataOfs += MarshalResp(&pThreadPar->pRespData[RespDataOfs],eRMIPTVarUint8,&PeakScoreCell,sizeof(tsSSWCell));
					}
				else
					RespDataOfs = 0;	

				JobRslt = JobResponse(InstanceID,ClassInstanceID, (UINT32)iRslt,RespDataOfs, pThreadPar->pRespData);
				break;

			case eSWMCombinedTargAlign:					// // method which combines the functionality of eSWMSetTarg, eSWMSetAlignRange, eSWMAlign, eSWMClassifyPath, eSWMTracebacksToAlignOps, eSWMAddMultiAlignment into a single method to reduce RMI overheads 
				tsCombinedTargAlignPars *pCombinedTargAlignPars;
				tsCombinedTargAlignRet CombinedTargAlignRet;
				etSeqBase *pTargSeq;			// target sequence 
#ifdef WIN32
				InterlockedIncrement(&m_NumSWAlignReqs);		// users are likely to be interested in the number of SW alignments requested
#else
				__sync_fetch_and_add(&m_NumSWAlignReqs,1);
#endif
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					int ReqDataOfs;
					ReqDataOfs = UnmarshalReq(sizeof(tsCombinedTargAlignPars),pThreadPar->pReqData,&pCombinedTargAlignPars);
					ReqDataOfs += UnmarshalReq(pCombinedTargAlignPars->TargSeqLen,&pThreadPar->pReqData[ReqDataOfs],&pTargSeq);
					pCombinedTargAlignPars->pTargSeq = pTargSeq;
					iRslt = pClassInstance->pClass->CombinedTargAlign(pCombinedTargAlignPars,&CombinedTargAlignRet);
					}
				else
					iRslt = -1;
				if(iRslt >= 0)
					RespDataOfs = MarshalResp(pThreadPar->pRespData,eRMIPTVarUint8,&CombinedTargAlignRet,sizeof(tsCombinedTargAlignRet));
				else
					RespDataOfs = 0;	

				JobRslt = JobResponse(InstanceID,ClassInstanceID, (UINT32)iRslt,RespDataOfs, pThreadPar->pRespData);
				break;


			case eSWMClassifyPath:			// ClassifyPath
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					int MaxArtefactDev;				// classify path as artefactual if sliding window of 500bp over any overlap deviates by more than this percentage from the overlap mean
					UINT32 ProbeStartOfs;			// alignment starts at this probe sequence offset (1..n)
					UINT32 ProbeEndOfs;				// alignment ends at this probe sequence offset
					UINT32 TargStartOfs;			// alignment starts at this target sequence offset (1..n)
					UINT32 TargEndOfs;				// alignment ends at this target sequence offset
					ReqDataOfs = UnmarshalReq(sizeof(INT32),pThreadPar->pReqData,&MaxArtefactDev);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&ProbeStartOfs);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&ProbeEndOfs);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&TargStartOfs);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&TargEndOfs);
					iRslt = pClassInstance->pClass->ClassifyPath(MaxArtefactDev,ProbeStartOfs,ProbeEndOfs,TargStartOfs,TargEndOfs);
					}
				else
					iRslt = -1;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, (UINT32)iRslt,0, NULL);
				break;

			case eSWMTracebacksToAlignOps:	// TracebacksToAlignOps
				tMAOp *pAlignOps;
				bool bAlignOps;
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					UINT32 ProbeStartOfs;			// alignment starts at this probe sequence offset (1..n)
					UINT32 ProbeEndOfs;				// alignment ends at this probe sequence offset
					UINT32 TargStartOfs;			// alignment starts at this target sequence offset (1..n)
					UINT32 TargEndOfs;				// alignment ends at this target sequence offset
					ReqDataOfs = UnmarshalReq(sizeof(UINT32),pThreadPar->pReqData,&ProbeStartOfs);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&ProbeEndOfs);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&TargStartOfs);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&TargEndOfs);
					ReqDataOfs += UnmarshalReq(sizeof(bool),&pThreadPar->pReqData[ReqDataOfs],&bAlignOps);
					iRslt = pClassInstance->pClass->TracebacksToAlignOps(ProbeStartOfs,ProbeEndOfs,TargStartOfs,TargEndOfs,bAlignOps ? &pAlignOps : NULL);
					}
				else
					{
					iRslt = -1;
					pAlignOps = NULL;
					bAlignOps = false;
					}

				if(bAlignOps && iRslt > 0)
					RespDataOfs = MarshalResp(pThreadPar->pRespData,eRMIPTVarUint8,pAlignOps,iRslt * sizeof(tMAOp));
				else
					RespDataOfs = 0;	

				JobRslt = JobResponse(InstanceID,ClassInstanceID, (UINT32)iRslt,RespDataOfs, pThreadPar->pRespData);
				break;

			case eSWMAddMultiAlignment:		// AddMultiAlignment
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					UINT32 ProbeStartOfs;			// alignment starts at this probe sequence offset (1..n)
					  UINT32 ProbeEndOfs;			// alignment ends at this probe sequence offset inclusive
					  UINT32 TargStartOfs;			// alignment starts at this target sequence offset (1..n)
					  UINT32 TargEndOfs;			// alignment ends at this target sequence offset inclusive
					  UINT32 TargSeqLen;			// target sequence length
					  etSeqBase *pTargSeq;			// alignment target sequence
					  UINT8 Flags;					// bit 7 set if target loaded as a high confidence sequence, bits 0..3 is weighting factor to apply when generating consensus bases

					ReqDataOfs = UnmarshalReq(sizeof(UINT32),pThreadPar->pReqData,&ProbeStartOfs);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&ProbeEndOfs);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&TargStartOfs);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&TargEndOfs);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&TargSeqLen);
					ReqDataOfs += UnmarshalReq(sizeof(etSeqBase **),&pThreadPar->pReqData[ReqDataOfs],&pTargSeq);
					ReqDataOfs += UnmarshalReq(sizeof(UINT8),&pThreadPar->pReqData[ReqDataOfs],&Flags);
					iRslt = pClassInstance->pClass->AddMultiAlignment(ProbeStartOfs,ProbeEndOfs,TargStartOfs,TargEndOfs,TargSeqLen,pTargSeq,Flags);
					}
				else
					iRslt = -1;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, (UINT32)iRslt,0, NULL);
				break;

			case eSWMGenMultialignConcensus:	// GenMultialignConcensus
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					iRslt = pClassInstance->pClass->GenMultialignConcensus();
				else
					iRslt = -1;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, (UINT32)iRslt,0, NULL);
				break;

			case eSWMMAlignCols2fasta:		// MAlignCols2fasta
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					UINT32 ProbeID;				// identifies sequence which was used as the probe when determining the multialignments
					INT32 MinConf;				// sequence bases averaged over 100bp must be of at least this confidence (0..9)
					INT32 MinLen;				// and sequence lengths must be of at least this length 
					UINT32 BuffSize;			// buffer allocated to hold at most this many chars
					ReqDataOfs = UnmarshalReq(sizeof(INT32),pThreadPar->pReqData,&ProbeID);
					ReqDataOfs += UnmarshalReq(sizeof(INT32),&pThreadPar->pReqData[ReqDataOfs],&MinConf);
					ReqDataOfs += UnmarshalReq(sizeof(INT32),&pThreadPar->pReqData[ReqDataOfs],&MinLen);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&BuffSize);
					if(BuffSize > cMaxMFABuffSize)
						BuffSize = cMaxMFABuffSize;				
					iRslt = pClassInstance->pClass->MAlignCols2fasta(ProbeID,MinConf,MinLen,BuffSize,pThreadPar->pszBuffer);
					}
				else
					iRslt = -1;
				if(iRslt > 0)
					RespDataOfs = MarshalResp(pThreadPar->pRespData,eRMIPTVarUint8,pThreadPar->pszBuffer,iRslt);
				else
					RespDataOfs = 0;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, (UINT32)iRslt,RespDataOfs, pThreadPar->pRespData);
				break;

			case eSWMMAlignCols2MFA:			// MAlignCols2MFA
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					{
					UINT32 ProbeID;		// identifies sequence which was used as the probe when determining the multialignments
					UINT32 BuffSize;	// buffer allocated to hold at most this many chars

					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&ProbeID);
					ReqDataOfs += UnmarshalReq(sizeof(UINT32),&pThreadPar->pReqData[ReqDataOfs],&BuffSize);
					if(BuffSize > cMaxMFABuffSize)
						BuffSize = cMaxMFABuffSize;	
					iRslt = pClassInstance->pClass->MAlignCols2MFA(ProbeID,BuffSize,pThreadPar->pszBuffer);
					}
				else
					iRslt = -1;
				if(iRslt > 0)
					RespDataOfs = MarshalResp(pThreadPar->pRespData,eRMIPTVarUint8,pThreadPar->pszBuffer,iRslt);
				else
					RespDataOfs = 0;
				JobRslt = JobResponse(InstanceID,ClassInstanceID, (UINT32)iRslt,RespDataOfs, pThreadPar->pRespData);
				break;

			default:			// currently any other method is not implemented
				if((pClassInstance = LocateClassInstance(ClassInstanceID))!=NULL)
					CUtility::SleepMillisecs(100);	// simulated some work in processing!	
				else
					ClassInstanceID = 0;
				break;
			}

		NumJobsProc += 1;
		}
	}
if(ClassInstanceID != 0)
	{
	FreeClassInstance(ClassInstanceID);
	ClassInstanceID = 0;
	}

#ifdef WIN32
InterlockedDecrement(&m_NumWorkerInsts);
#else
__sync_fetch_and_sub(&m_NumWorkerInsts,1);
#endif
return(0);
}



