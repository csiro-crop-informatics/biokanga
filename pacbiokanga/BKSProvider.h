#pragma once

#include "BKScommon.h"

const UINT32 cServiceProviderVersion = 1;			// service provider is at this version
const UINT32 cMaxServiceProviderInsts = cMaxServiceInsts;	    // limited to support a maximum of this many service instances
// when negotiating with potential service requesters then minimal buffer tx/rx buffer sizes are allocated
const int cMinTxRxBuffSize = (cMaxServiceTypes * sizeof(tsServiceDetail)) + sizeof(tsBKSReqServices) * 3;	// always allocate at least this sized TxdBuff/RxdBuffs - ensures negotiation frames fit!

const int cMaxReqDataSize =  cMaxSWReqPayloadSize;		// each worker thread allocates to process up to this much request data
const int cMaxReqParamSize = cMaxSWParamLen;			// each worker thread allocates to process up to this much parameterisation data
const int cMaxRespDataSize = cMaxSWRespPayloadSize;		// each worker thread allocates to return up to this much response data
const int cMaxMFABuffSize =  cMaxSWMAFBuffSize;			// each worker thread allocates to hold at most this sized MAlignCols2fasta/MAlignCols2MFA alignments plus row descriptor prefixes

// service providers will be in one of these exclusive states
typedef enum TAG_eBKSPProvState
{
	eBKSPSUndefined = 0,				// service provider state is undefined, yet to connect with a server
	eBKSPSWaitReqServices,				// connected with server and waiting for server to send list of required services
	eBKSPSSendOfferedService,			// sending server the offered service
	eBKSPSWaitAcceptService,			// waiting for server to accept offered service
	eBKSPSAcceptedServiceActv,			// server has accepted offered service so can now actively process service requests
	eBKSPSAcceptedServiceTerm,			// session is in the process of terminating
	eBKSPSPlaceHolder
} teBKSPProvState;

#pragma pack(1)
typedef struct TAG_sTxdRxd
{
	UINT32 SessionID;			// uniquely identifying this session between service requester and provider, requester specified
	UINT16 flgSelMonExcept : 1; // socket to be monitored for exceptions
	UINT16 flgSelMonRead : 1;	// socket to be monitored for received data available to be read
	UINT16 flgSelMonWrite : 1; // socket to be monitored for write completed
	UINT16 flgKeepAliveReq: 1; // set if a keepalive is to be sent
	UINT16 flgRxCplt : 1;	// set when complete frame received 	
	UINT16 flgTxCplt : 1;	// set when complete frame sent
	UINT16 flgErr : 1;       // set on any unrecoverable error
	UINT16 flgErrReason : 4;	// holds reason for socket level error flag set
	socket_t  Socket;		// assumed connected socket
	time_t PacRxdAtSecs;	// the time at which a frame was last received, used for determining if session still active
	time_t PacTxdAtSecs;	// the time at which a frame was last sent, used for keep alive generation

	UINT8 RxdRxFrameID;		// RxFrameID of last received frame - will be 0 if no frame yet received by session peer
	UINT8 RxdTxFrameID;		// TxFrameID of last received frame - will be 0 if no frame yet received from session peer
	UINT8 TxFrameID;		// FrameID of next frame to send, will wraparound back to 1 when incremented past 0x07f 
	UINT32 TotRxd;			// pRcvBuff contains this total number of received bytes
	UINT32 CurPacRxd;		// number of bytes in currently received packet
	UINT32 TotTxd;			// total number of bytes in pTxdBuff to be sent
	UINT32 CurTxd;			// number of bytes currently sent from pTxdBuff
	UINT32 AllocdRxdBuff;	// pRxdBuff allocated to hold at most this many received bytes
	UINT32 AllocdTxdBuff;	// pTxdBuff allocated to hold at most this many bytes
	UINT8 *pRxdBuff;		// receiving data into this buffer
	UINT8 *pTxdBuff;		// sending data from this buffer
    SOCKADDR_STORAGE  IPaddress;	// remote IP address + port of endpoint service provider (IPv4 or IPv6)
} tsTxdRxd;

typedef struct TAG_sReqResp {
	INT64 JobIDEx;						// request identifier as received from server for this request instance
	UINT64 ClassInstanceID;				// class instance referenced
	UINT32 ClassMethodID;				// identifies class method
	UINT32 InstanceID;					// service instance specific identifier
	UINT32 InstanceIDEx;				// combination of both the InstanceID (in bits 0..9) and ......
	UINT32 flgReqAvail : 1;				// this service instance is available for processing 
	UINT32 flgProc: 1;					// this service instance is currently being processed
	UINT32 flgCpltd: 1;					// service processing has completed and resultset can be sent back to service requester
	UINT32 JobRslt;						// completion result
	UINT32 ParamSize;					// instance specific parameter size
	UINT32 InDataSize;					// instance specific input data size
	UINT32 OutDataSize;					// instance specific result data size
	UINT8 Data[1];						// when service requested then parameters followed by input data, if service response then response result data
} tsReqResp;

typedef struct TAG_sBKSRegSessionEx
{
	UINT8 BKSPType;						// session is providing this teBKSType service
	UINT8 BKSPState;					// session is currently in this teBKSPProvState registration state 
	UINT32 NumInstances;				// session can process at most this many service instances
	UINT32 MaxReqPayloadSize;			// request payloads from the service requester, including framing, can be up to this size (UINT8s),
	UINT32 MaxRespPayloadSize;			// response payloads from the service provider to requester, including framing, can be up to this  size (UINT8s)
	UINT32 ReqRespInstSize;				// each service instance requires this much memory
	UINT32 KeepaliveSecs;				// expecting packet activity from endpoint with periodicity of no more than this number of seconds
	UINT32 InstancesBusy;				// currently this many requests being serviced in this session (total of instances ready to be processed, being processed, and processing has completed)
	UINT32 InstancesReqAvail;			// current number of requests ready to be processed
	UINT32 InstancesProc;				// current number of requests currently being processed
	UINT32 InstancesCpltd;				// current number of instances completed processing and ready for resultsets to be sent back to requester
	UINT32 TotNumRequests;				// total number of requests processed by this service provider in the current session
	UINT32 AllocdReqResp;				// allocation size for pReqResp
	UINT8 *pReqResp;					// allocation for NumInstances of requests and associated responses (tsReqResp's)
	tsTxdRxd TxdRxd;					// holding low level send/receive buffers + connected socket
} tsBKSRegSessionEx;


typedef struct TAG_sBKSType
{
	teBKSPType BKSPType;					// registering this service type
	UINT32 ProviderVersion;					// service provider version 
	UINT32 MaxServiceInsts;					// limited to support a maximum of this many service instances
	UINT32 MaxQuerySeqLen;					// accepting query sequences up to this length
	UINT32 MaxTargSeqLen;					// accepting target sequences up to this length
	UINT32 MaxReqPayloadSize;				// request payloads from the service requester, including framing, can be up to this size (UINT8s),
	UINT32 MaxRespPayloadSize;				// response payloads from the service provider to requester, including framing, can be up to this  size (UINT8s)

	double ScalePayload;					// ScaledPayload = (MaxReqPayloadSize + MaxTargSeqLen) * ScalePayload
	double ScaleTargLen;					// ScaledTargLen = MaxTargSeqLen * ScaleTargLen
	double ScaleQueryLen;					// ScaledQueryLen = MaxTargSeqLen * ScaleQueryLen
	double ScaleScaledTargQuery;			// ScaledScaledTargQuery = (ScaledTargLen * ScaledQueryLen) * ScaleScaledTargQuery;
	double AvailResources;					// OfferedInstances = AvailResources / (ScaledScaledTargQuery + ScalePayload);
	UINT32(*pfnCostNumInstances)(teBKSPType BKSPType,			// costing function for determining number of service instances to offer for this service type
										UINT32 MaxServiceInsts,			// limited to support a maximum of this many service instances
										UINT32 MaxQuerySeqLen,			// accepted query sequences can be up to this length
										UINT32 MaxTargSeqLen,			// accepted target sequences can be up to this length
										UINT32 MaxReqPayloadSize,		// request payloads from the service requester, including framing, can be up to this size (UINT8s),
										UINT32 MaxRespPayloadSize,		// response payloads from the service provider to requester, including framing, can be up to this  size (UINT8s)
										double ScalePayload,			// ScaledPayload = (MaxReqPayloadSize + MaxTargSeqLen) * ScalePayload
										double ScaleTargLen,			// ScaledTargLen = MaxTargSeqLen * ScaleTargLen
										double ScaleQueryLen,			// ScaledQueryLen = MaxTargSeqLen * ScaleQueryLen
										double ScaleScaledTargQuery,	// ScaledScaledTargQuery = (ScaledTargLen * ScaledQueryLen) * ScaleScaledTargQuery
										double AvailResources);			// OfferedInstances = AvailResources / (ScaledScaledTargQuery + ScalePayload)
} tsBKSType;

typedef struct TAG_sWorkerInstance {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to class instance
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	UINT32 threadID;				// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// processing result
	UINT8 *pReqData;				// malloc'd (cMaxReqDataSize) to hold request data
	UINT8 *pParamData;              // malloc'd (cMaxReqParamSize) to hold any paramertisations
	UINT8 *pRespData;               // malloc'd (cMaxRespDataSize) to hold response data
	char *pszBuffer;                // malloc'd )(cMaxMFABuffSize) to hold any textual alignment sequences
} tsWorkerInstance;

typedef struct TAG_sClassInstance {
	UINT64 ClassInstanceID;				// monotonically incremented 
	time_t LastAccessed;				// when this class instance was last accessed, used to expire unused class instances
	CSSW *pClass;						// if non-null then pts to instantiated class instance	
	} tsClassInstance;

#pragma pack()


class CBKSProvider
{

	int m_MaxConnWait;									// wait for at most this many minutes for connection
	int m_MaxServInsts;									// max number of service instances supported
	int m_MaxServMemGB;									// max allocatable memory available for all service instances
	bool m_bCreatedMutexes;								// set true after serialisation locks/mutexes initialised
	bool m_bTermConnectionReq;							// connection termination has been requested and is either being progressed or termination has completed
	bool m_bNotifiedReqs;								// set TRUE if control message has been sent to control socket 2
	UINT32 m_NumInitTypes;								// number of service types initialised in m_BKSTypes[]
	tsBKSType m_BKSTypes[eBKSPTPlaceHolder - 1];		// entry for each potentially supported service type indexed by Type-1

	tsBKSRegSessionEx m_BKSConnection;							// server connection

	UINT32 m_ReqNumWorkerInsts;							// number of worker instance threads to start
#ifdef WIN32
	alignas(4) volatile UINT32 m_NumPendCpltd;			// count of threads waiting to gain lock in JobResponse() to notify job completed response 
	alignas(4) volatile UINT32  m_NumWorkerInsts;				// number of worker instance threads actually started
	alignas(4) volatile UINT32  m_NumSWAlignReqs;				// number of SW alignments requested
	alignas(4) volatile UINT32 m_TermAllThreads;                // will be set to 1 if all worker threads are to terminate
#else
	__attribute__((aligned(4))) volatile UINT32 m_NumPendCpltd;			// count of threads waiting to gain lock in JobResponse() to notify job completed response 
	__attribute__((aligned(4))) volatile UINT32  m_NumWorkerInsts;				// number of worker instance threads actually started
	__attribute__((aligned(4)))  volatile UINT32  m_NumSWAlignReqs;				// number of SW alignments requested
	__attribute__((aligned(4))) volatile UINT32 m_TermAllThreads;                  // will be set to 1 if all worker threads are to terminate
#endif
	tsWorkerInstance m_WorkerInstances[cMaxServiceInsts];	// to hold all worker instance thread parameters

	UINT32 m_MaxClassInsts;									// at most this many class instances can be instantiated
	UINT32 m_NumClassInsts;									// currently this many class instances have been instantiated in m_ClassInstances[];
	UINT64 m_HiClassInstanceID;								// highest class instance identifier thus far allocated							
	tsClassInstance m_ClassInstances[cMaxClassInsts];		// all possible class instances

	tsTxdRxd m_Ctrl[2];									// Ctrl[0] written to by threads needing to signal select() processing thread, select() processing thread monitors m_Ctrl[1]

	teBSFrsltCodes										// cBSFSuccess if no errors and registration process is continued, cBSFSocketErr if any errors and connection has been terminated, eBSFerrMem if unable to allocate memory
		ProcessSessEstab(bool bCpltdWrite);				// false if frame received, true if frame sent

	int													// number of received requests allocated to service instances
		HandleServiceRequests(void);					// handler for received requests for allocation of a service instance to process the service request

	int                                         // number of responses assembled into m_BKSConnection.TxdRxd.pTxdBuff ready to be sent back to requester
		HandleServiceResponses(void);			// locate those service instances with responses ready to be sent to the requester and assemble responses into m_BKSConnection.TxdRxd.pTxdBuf


	int	TerminateConnection(bool bFreeMem,		  // true if any allocated memory to be free'd
							bool bCloseCtrl,	  // close control socket pair
							bool bResetInitTypes); // reset all service types

	int	CreateMutexes(void);				// create mutexes used in access serialisations
	int	DeleteMutexes(void);				// delete mutexes used in access serialisations

	void AcquireLock(bool bExclusive);	// lock for serialised access by multiple concurrent reader threads (bExclusive == false), or serialised access by single thread (bExclusive == true)
	void ReleaseLock(bool bExclusive);	// release serialised access lock

#ifdef _WIN32
	SRWLOCK m_hRWLock;					// serialising multiple reader access but single writer
#else
	pthread_rwlock_t m_hRWLock;			// serialising multiple reader access but single writer
#endif
#ifdef WIN32
	alignas(4) volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access 
#else
	__attribute__((aligned(4))) volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access
#endif
	void AcquireCASSerialise(bool bPriority = false);	// if bPriority true then backoff time is reduced relative to if bPriority is false, increasing the probability of acquiring the serialisation lock if there is contention
	void ReleaseCASSerialise(void);


	int				// 1 if timed out attempting to connect to server, 0 if terminated because requested to terminate, -1 if socket level errors
		ConnectServer(int MaxConnWait,							 // wait for at most this many minutes for connection
						const char* pszHost, char *pszService);  // connect to server at pszHost:pszService

	bool InitialiseCtrlSocks(void); // initialise the control sockets in m_Ctrl[]
																// InitialiseConnect 
	bool					// true if connecting socket initialised
		InitialiseConnect(int MaxConnWait,			// wait for at most this many minutes for connection
							const char* pszHost,	// connecting to this server host/IP address; NULL to use first INET IP local to this machine
							char *pszService);		// connecting to this server service/port; NULL to use default port 

	bool RxData(tsTxdRxd *pRxd);
	bool TxData(tsTxdRxd *pTxd);

	bool NotifyCtrl(UINT8 Msg = 0);			// notify via control sockets that, default, there is at least 1 response available to be sent to service requester
	bool SendCtrlMsg(int Len,				// number of control message bytes to send ptd at by pCtrlMsg, can be 0
								UINT8 *pCtrlMsg);	// pCtrlMsg pts to control message bytes
	int								// received message is this length, 0 if no outstanding messages, < 0 if errors
			RcvCtrl(int BuffLen,			// available buffer into which the control message can be copied 
		UINT8 *pBuff);			// copy control message into this buffer
	bool ProcessCtrlMsg(int MsgLen,UINT8 *pMsg);			// process a message received by control socket m_Ctrl[1] - note: currently these messages are simply discarded

	bool ShutdownConnection(socket_t *pSocket);


	int			// returns 0 if no sockets to be monitored with select(), on windows the total number of monitored sockets, on linux the highest socket file descriptor plus 1
		SetupFDSets(fd_set& ReadFDs,			// select() read available socket descriptor set  
					fd_set& WriteFDs,			// select() write accepted socket descriptor set
					fd_set& ExceptFDs);			// select() exceptions descriptor set
#ifdef WIN32
	const char *WSAGetLastErrorMessage(const char* pcMessagePrefix,int nErrorID = 0);
#endif
	int	Reset(void);

	int Initialise(int MaxConnWait,				// wait for at most this many minutes for connection
					int MaxServInsts,			// max number of service instances supported
						int MaxClassInsts,		// max number of class instances supported (may differ from MaxServInsts)
						int MaxServMemGB);		// max allocatable memory available for all service instances

	int StartWorkerThreads(UINT32 NumInstances,	// this number of worker threads required
						   UINT8 BKSPType);		// workers are providing this service type

	int TerminateWorkerThreads(void);			// terminate all worker threads

	int											// marshaled parameter required this many bytes
		MarshalResp(UINT8 *pInto,				// marshal into this list
				teRMIParamType Type,				// parameter type
				void *pValue,					// parameter value
				UINT32 ValLen);					// length of parameter ptd to by pValue, only used if parameter type is pUint8

	int
		UnmarshalReq(UINT32 DataLen,
					UINT8 *pFrom,		// unmarshal from this marshalled parameter list
					void *pValue);

public:
	CBKSProvider();
	~CBKSProvider();

	int ProcWorkerThread(tsWorkerInstance *pThreadPar);  // worker thread startup entry for processing requests

	int					// returns total number of registered service types or teBSFrsltCodes error code if any parameterisation errors or already registered type
		RegServiceType(teBKSPType BKSPType = eBKSPTSmithWaterman,			// registering this service type
					   UINT32 ProviderVersion = cServiceProviderVersion,	// service provider version 
					   UINT32 MaxServiceInsts = cMaxServiceProviderInsts,	// limited to support a maximum of this many service instances
					   UINT32 MaxQuerySeqLen = cMaxSWQuerySeqLen,			// accepted query sequences can be up to this length
					   UINT32 MaxTargSeqLen = cMaxSWTargSeqLen,				// accepted target sequences can be up to this length
					   UINT32 MaxReqPayloadSize = cMaxSWReqPayloadSize,		// request payloads from the service requester, including framing, can be up to this size (UINT8s),
					   UINT32 MaxRespPayloadSize = cMaxSWRespPayloadSize,	// response payloads from the service provider to requester, including framing, can be up to this  size (UINT8s)
						double ScalePayload = 1.0,					// ScaledPayload = (MaxReqPayloadSize + MaxTargSeqLen) * ScalePayload
						double ScaleTargLen = 1.0,					// ScaledTargLen = MaxTargSeqLen * ScaleTargLen
						double ScaleQueryLen = 1.0,					// ScaledQueryLen = MaxTargSeqLen * ScaleQueryLen
						double ScaleScaledTargQuery = 0.05,			// ScaledScaledTargQuery = (ScaledTargLen * ScaledQueryLen) * ScaleScaledTargQuery;
						double AvailResources = 256000000000.0,				// OfferedInstances = AvailResources / (ScaledScaledTargQuery + ScalePayload);
					   UINT32  (*pfnCostNumInstances)(teBKSPType BKSPType,				// determine number of service instances to offer for this service type
													  UINT32 MaxServiceInsts,			// limited to support a maximum of this many service instances
													  UINT32 MaxQuerySeqLen,			// accepted query sequences can be up to this length
													  UINT32 MaxTargSeqLen,			// accepted target sequences can be up to this length
													  UINT32 MaxReqPayloadSize,		// request payloads from the service requester, including framing, can be up to this size (UINT8s),
													  UINT32 MaxRespPayloadSize,		// response payloads from the service provider to requester, including framing, can be up to this  size (UINT8s)
													  double ScalePayload,					// ScaledPayload = (MaxReqPayloadSize + MaxTargSeqLen) * ScalePayload
													  double ScaleTargLen,					// ScaledTargLen = MaxTargSeqLen * ScaleTargLen
													  double ScaleQueryLen,				// ScaledQueryLen = MaxTargSeqLen * ScaleQueryLen
													  double ScaleScaledTargQuery,			// ScaledScaledTargQuery = (ScaledTargLen * ScaledQueryLen) * ScaleScaledTargQuery
													  double AvailResources) = NULL);				// OfferedInstances = AvailResources / (ScaledScaledTargQuery + ScalePayload)

	
	// after all service types offered have been registered then Process() can be invoked													
	int 
		Process(int MaxConnWait,						// wait for at most this many minutes for connection 
				int MaxServInsts,						// max number of service instances supported
				int MaxClassInsts,						// max number of class instances supported (may differ from MaxServInsts)
				int MaxServMemGB,						// max allocatable memory available for all service instances 
				const char* pszHost = NULL,				// connect to this host/IP address; NULL to use first INET IP local to this machine
				char *pszService = NULL);				// connect to this service/port; NULL to use default port 


	int			// 0 no service request job available, 1 job details are being returned, -1 job is available but MaxParamSize < required, -2 job available but MaxRequestData < required, -3 if session being terminated
			GetJobToProcess(UINT32 *pInstanceID,		// service instance identifier, the JobResponse() must use this identifier
						UINT64 *pClassInstanceID,	// returned class instance to apply job to
						UINT32 *pClassMethodID,		// returned class method to apply on class instance
					    UINT32 *pMaxParamsSize,		// on input, max sized parameter block accepted, on return then the actual size of the parameter block
						UINT8 *pParams,				// returned parameter block
						UINT32 *pMaxRequestData,	// on input the max sized request data block accepted, on return the actual size of the request data block
						UINT8 *pRequestData);		// returned request data block


	int			// 0 if response accepted, -1 if job does not exist or parameterisation errors, -3 if session terminating 
			JobResponse(INT32 InstanceID,			// service instance identifier returned by GetJobToProcess
					UINT64 ClassInstanceID,		// response is for this class instance
				    UINT32 ProcRslt,			// service request processing result
					UINT32 ResponseSize,		// response data block size
					UINT8 *pResponseData);		// response data


	tsClassInstance *LocateClassInstance(UINT64 ClassInstanceID);
	tsClassInstance *AllocClassInstance(void); // allocate a new tsClassInstance and initialise with ClassInstanceID 
	bool FreeClassInstance(UINT64 ClassInstanceID); // free a previously allocated class instance

};


