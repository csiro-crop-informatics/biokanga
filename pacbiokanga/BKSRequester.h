#pragma once

#include "BKScommon.h"

const int cListenBacklog = 5;					// maximum length of the queue of pending connections

const int cMaxConcurrentRequests = min(4095,cMaxServiceInsts * cMaxNumSessions);	// can process at most this many concurrent service requests over all session instances independent of service type
const int cMaxReqID = cMaxConcurrentRequests;	// request identifiers will range from 1..cMaxConcurrentRequests

// when negotiating with potential service providers then minimal buffer tx/rx buffer sizes are allocated
const int cMinTxRxBuffSize = (cMaxServiceTypes * sizeof(tsServiceDetail)) + sizeof(tsBKSReqServices) * 3;	// always allocate at least this sized TxdBuff/RxdBuffs - ensures negotiation frames fit!

// Session service providers will be in one of these exclusive states
typedef enum TAG_eBKSEPProvState
{
	eBKSPSUndefined = 0,				// service provider state is undefined waiting on a potential service provider to connect
	eBKSPSRegisteringNeg,				// initiated negotiation of capabilities with connected service provider
	eBKSPSRegisteredActv,				// completed capabilities negotiation and provider can now actively process service requests
	eBKSPSRegisteredTerm,				// currently terminating service provider connection
	eBKSPSPlaceHolder
} teBKSPEPProvState;

typedef enum TAG_eSessEstabState
{
	eSESnone = 0,			// not yet started
	eSESTxReqServices,		   // send required service type details to potential provider
	eSESRxOfferedService,	   // expecting to receive potential providers offered service
	eSESTxAcceptService,   // sending acceptance of offered service to provider
	eSESTxRejectService,   // sending rejection of offered service to provider
	eSESPlaceHolder
} teSessEstabState;


#pragma pack(1)

typedef INT64 tJobIDEx;			// >0 extended job identifier; 0 if none assigned, <0 if errors

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


typedef struct TAG_sBKSSessEstab
{
	UINT8 SEState;			// current session establishment state - eSESnone, eSESTxServices, eSESRxOfferService, eSESTxAcceptService 
	time_t StartSecs;		// started session establishment at this time
	teBKSPType BKSPType;	// confirmation of service type being accepted or eBKSPTUndefined if offered type not accepted 
	UINT32 MaxInstances;	// max number of instances supported
	UINT32 MaxClassInstances;// Session can instantiate at most this many class instances
	tsTxdRxd TxdRxd;		// holding low level send/receive buffers + connected socket
} tsBKSSessEstab;

// Session service providers are registered into the following
typedef struct TAG_sBKSRegSession
{
	UINT32 SessionID;					// uniquely identifying this session between service requester and provider
	UINT32 TypeSessionID;				// used as an index within the type to reference this session
	UINT8 BKSPType;						// Session is providing this teBKSType service
	UINT8 BKSPState;					// Session is currently in this teBKSPEPProvState registration state 
	UINT32 MaxInstances;				// Session can process at most this many service instances
	UINT32 MaxClassInstances;			// Session can instantiate at most this many class instances
	UINT32 NumBusy;						// total number of instances currently committed to requests, processing, or with responses
	UINT32 NumReqs;						// number of instances with service request ready to send to a service provider
	UINT32 NumProcs;					// number of instances with service requests currently being processed by service provider
	UINT32 NumCpltd;					// number of instances with service requests finished processing with response available
	UINT64 TotNumCpltd;					// total number of completed requests in this session
} tsBKSRegSession;

typedef struct TAG_sReqRespInst
{
	INT64 JobIDEx;		// unique job identifier as was generated when job was accepted for submission
	UINT32 SubmitAt;	// when job was accepted for submission
	UINT32 CpltdAt;		// when job was returned as completed
	UINT16 ReqID;		// instance specific identifier
	UINT16 FlgReq:1;	// this job instance has been initialised and is ready to be sent for processing
	UINT16 FlgProc:1;   // this job instance is currently being processed by service provider
	UINT16 FlgCpltd:1;  // service provider has completed processing and returned results are available
	UINT32 JobRslt;			// service provider completion result
	UINT64 ClassInstanceID;	// class instance referenced
	UINT32 ClassMethodID;	// identifies class method
	UINT32 ParamSize;		// instance specific parameter size
	UINT32 InDataSize;		// instance specific input data size
	UINT32 OutDataSize;		// instance specific result data size
	UINT8 Data[1];		// when service requested then parameters followed by input data, if service response then response result data
} tsReqRespInst;

typedef struct TAG_sBKSRegSessionEx
{
	struct TAG_sBKSRegSessionEx *pNext;	// Sessions are linked as will be dynamically allocated
	tsBKSRegSession Session;			// Session

	UINT32 LastChkdReqIdx;				// start checking for service requests to send from this instance index
	UINT32 NumClassInstances;			// number of class instance identifiers currently in ClassInstanceIDs
	UINT64 ClassInstanceIDs[cMaxClassInsts];	// holds all instantiated class instance identifiers for this session
	UINT32 AllocdReqResp;				// allocation size for pReqResp
	UINT8 *pReqResp;					// allocation for NumInstances of requests and associated responses
	tsTxdRxd TxdRxd;					// holding low level send/receive buffers + connected socket
} tsBKSRegSessionEx;


typedef struct TAG_sBKSType
{
	tsServiceDetail Detail;				// service type detail
	UINT32 MaxSessions;					// allowing at most this many sessions of specified service type
	UINT32  NumSessions;				// currently this number of registered sessions of this type
	UINT32 NumInstances;				// totaling this many service instances
	UINT32 ReqRespInstSize;				// each service instance requires this much memory
	UINT32 FlgReq:1;					// at least 1 session has a service request ready to send to a service provider
	tsBKSRegSessionEx *pFirstSession;	// ptr to first session of this type
	tsBKSRegSessionEx *pSessions[cMaxNumSessions]; // ptr to each session as indexed by TypeSessvIdx  
} tsBKSType;


typedef struct TAG_sChkPtReq {	// individual check pointed request
	UINT32 ClassMethodID;	// identifies class method
	UINT32 ParamSize;		// instance specific parameter size
	UINT32 InDataSize;		// instance specific input data size
	UINT32 MaxOutDataSize;	// instance specific max expected output data for this method
	UINT8 Data[1];			// parameters followed by input data
	} tsChkPtReq;

typedef struct TAG_sChkPtReqs { // list of all check pointed requests for a specific class instance
	UINT32 SessionID;			// uniquely identifying this session between service requester and provider (set to 0 if this class instance can be reused)
	UINT32 TypeSessionID;		// used as an index within the type to reference this session
	UINT64 ClassInstanceID;		// check pointing this class instance
	UINT32 JobRslt;				// last job processing result as returned by service provider
	UINT32 UsedRespData;		// pRespData currently holds this sized response for last job submitted
	UINT32 AllocRespData;		// pRespData allocated to this size
    UINT8 *pRespData;			// response data for last check pointed method successfully processed by service provider
	UINT32 NumChkPtReqs;		// number of check pointed requests in ChkPtReqs[]
	size_t AllocMem;			// total memory allocated to hold this check pointed list (includes this header)
	size_t UsedMem;				// currently using this sized memory to hold check pointed list (includes this header)
	tsChkPtReq ChkPtReqs[1];	// 1st check pointed request
	} tsChkPtReqs;

typedef struct TAG_sRequesterThread {
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
} tsRequesterThread;

#pragma pack()

class CBKSRequester
{
	char m_szHostName[cMaxHostNameLen+1];				// host on which to listen for connections
	char m_szServiceName[cMaxServiceNameLen+1];			// listening on this port

	tsRequesterThread m_RequesterThread;				// requester thread parameters

#ifdef WIN32
#ifdef _ChkLockDepth_			// if checking that serialisation locks are actually working!
	alignas(4) volatile	UINT32 m_LockDepth;
	alignas(4) volatile	UINT32 m_CASSerialiseDepth;
#endif
	alignas(4) volatile	UINT32 m_CreatedMutexes;					// set to 1 after serialisation locks/mutexes initialised
	alignas(4) volatile	UINT32 m_ThreadActive;						// set to 1 after requester thread has started and initialised
	alignas(4) volatile UINT32 m_Terminate;							// set to 1 if requester thread is required to self terminate
	alignas(4) volatile UINT32 m_LockEnabled;						// set to 1 if AcquireLock has been initialised and can be used for serialisation
	alignas(4) volatile UINT32 m_CASSerialise;						// used with synchronous compare and swap (CAS) for serialising access 
	alignas(4) volatile UINT32 m_NumPendingReqs;					// number of outstanding pending requests from RMI worker threads to be processed
	alignas(4) volatile UINT32 m_NumPendingResps;					// number of outstanding pending response checks from RMI worker threads to be processed
	alignas(4) volatile LONG m_TotRespsAvail;						// current total number of responses available to RMI worker threads over all sessions

#else
#ifdef _ChkLockDepth_		// if checking that serialisation locks are actually working!
	__attribute__((aligned(4))) volatile	UINT32 m_LockDepth;
	__attribute__((aligned(4))) volatile	UINT32 m_CASSerialiseDepth;
#endif
	__attribute__((aligned(4))) volatile UINT32 m_CreatedMutexes;	// set to 1 after serialisation locks/mutexes initialised
	__attribute__((aligned(4))) volatile UINT32 m_ThreadActive;		// set to 1 after requester thread has started and initialised
	__attribute__((aligned(4))) volatile UINT32 m_Terminate;		// set to 1 if requester thread is required to self terminate
	__attribute__((aligned(4))) volatile UINT32 m_LockEnabled;		// set to 1 if AcquireLock has been initialised and can be used for serialisation
	__attribute__((aligned(4))) volatile UINT32 m_CASSerialise;		// used with synchronous compare and swap (CAS) for serialising access 
	__attribute__((aligned(4))) volatile UINT32 m_NumPendingReqs;	// number of outstanding pending requests from RMI worker threads to be processed
	__attribute__((aligned(4))) volatile UINT32 m_NumPendingResps;	// number of outstanding pending response checks from RMI worker threads to be processed
	__attribute__((aligned(4))) volatile UINT32 m_TotRespsAvail;	// current total number of responses available to RMI worker threads over all sessions
#endif

	bool m_bNotifiedReqs;								// set TRUE if control message has been sent to control socket 2

	bool m_bSessionTermReq;								// true whilst a Session  is being deleted

	UINT32 m_NumChkPtReqs;								// m_pChkPtReqs contains this many currently active check pointed class request ptrs
	UINT32 m_MaxChkPtReqs;								// m_pChkPtReqs allocated to hold at most this many check pointed class request ptrs
	tsChkPtReqs **m_ppChkPtReqs;							// allocated to hold check pointed class request ptrs

	UINT32 m_NumInitTypes;								// number of service types initialised in m_pBKSTypes
	tsBKSType *m_pBKSTypes;								// pts to allocated array of entries for each potentially supported service type indexed by Type-1

	UINT32 m_NumSessEstabs;								// number of sessions currently in pBKSSessEstab being negotiated with
	tsBKSSessEstab *m_pBKSSessEstabs;					// pts to allocated array of session peer endpoints currently being negotiated with
 
	UINT32 m_NumSessions;								// currently there are this many registered sessions

	UINT32 m_SessionIDVect[(cMaxConcurrentSessions+31)/32];	// bit vector containing all possible session identifiers, if bit set then that session identifier has been allocated and is in use
	UINT32 m_ReqIDVect[(cMaxReqID+31) / 32];				// bit vector containing all possible request identifiers, if bit set then that request identifier has been allocated and is in use

	socket_t m_ListenerSock;							// listening on this socket for Session connections

	tsTxdRxd m_Ctrl[2];									// Ctrl[0] written to by threads needing to signal select() processing thread, select() processing thread monitors m_Ctrl[1]

	UINT32												// returned request identifier or 0 if all identifiers have already been allocated
			AllocReqID(void);							// returns next available unused request identifier and sets that identifier in m_ReqIDVect[] as now allocated

	bool												// true if ReqID was previously allocated (used), false if previously unallocated (unused)
			UnallocReqID(UINT32 ReqID);					// sets the ReqID in m_ReqIDVect as unallocated and available to be reused

	bool												// true if ReqID is currently allocated (used), false if unallocated (unused)
			IsAllocReqID(UINT32 ReqID);					// checking this ReqID in m_ReqIDVect

	UINT32												// returned session identifier or 0 if all identifiers have already been allocated
		AllocSessionID(void);							// returns next available unused session identifier and sets that identifier in m_SessionIDVect as allocated

	bool				// true if SessionID was previously allocated (used), false if previously unallocated (unused)
		UnallocSessionID(UINT32 SessionID);   // sets the SessionID in m_SessionIDVect as unallocated and available to be reused

	bool												// true if SessionID is currently allocated (used), false if unallocated (unused)
			IsAllocSessionID(UINT32 SessionID);			// checking this SessionID in m_SessionIDVect

	tsBKSRegSessionEx *LocateSession(UINT32 SessionID);	// locates an established session which is identified by SessionID

	void
		ResetSessEstab(tsBKSSessEstab *pSessEstab,			// reset this tsBKSSessEstab instance
									  bool bKeepAllocs = true,// if true then retain any existing buffer allocations
									  bool bKeepSessionID = false);	// if true then retain SessionID as session has been accepted
	bool
		StartSessEstab(UINT32 SessionID,								// Session identifier for this potential Session
									socket_t Socket,					// communicating with Session over this connected socket
									SOCKADDR_STORAGE  *pIPaddress);		// peer is at this socket network address

	bool
		ProgressSessEstab(tsBKSSessEstab *pSessEstab,
										bool bCpltdWrite);	// false if frame received, true if frame sent

	bool
		AcceptFullSession(tsBKSSessEstab *pSessEstab);		// accepting session being established as full session ready for normal payload request/responses

	int	TerminateAllSessions(void);

	int				// returns number of sessions deleted
		DeleteAllSessionsInState(teBKSPEPProvState TermState = eBKSPSRegisteredTerm);	// terminate and delete any session in specific state (normally eBKSPSRegisteredTerm) or with an invalid socket

	int	CreateMutexes(void);				// create mutexes used in access serialisations
	int	DeleteMutexes(void);				// delete mutexes used in access serialisations

	bool AcquireLock(bool bExclusive);	// lock for serialised access by multiple concurrent reader threads (bExclusive == false), or serialised access by single thread (bExclusive == true)
	bool ReleaseLock(bool bExclusive);	// release serialised access lock

	void AcquireCASSerialise(bool bPriority = false);	// if bPriority true then backoff time is reduced relative to if bPriority is false, increasing the probability of acquiring the serialisation lock if there is contention
	void ReleaseCASSerialise(void);

#ifdef _WIN32
	SRWLOCK m_hRWLock;					// serialising multiple reader access but single writer
#else
	pthread_rwlock_t m_hRWLock;			// serialising multiple reader access but single writer
#endif

	bool InitialiseListener(const char* pszHost, char *pszService);  // initialise listening socket
	bool InitialiseCtrlSocks(void); // initialise the control sockets in m_Ctrl[]
	int	AcceptConnections(void);	 // start accepting connections
	int	SendRequestFrames(void);			// iterate all sessions and if any frames ready to send and room to accept the frame in TxdBuff then initiate the sending
	int  // 0: accepted frame, -1: JobIDEx errors, -2 Session or type errors, -3 mismatch between instance JobIDEx's, ClassInstanceID mismatch
			ProcessResponseFrame(tsBKSRegSessionEx *pSession);			// process a received response frame on this session
	bool RxData(tsTxdRxd *pRxd);
	bool TxData(tsTxdRxd *pTxd);

	bool NotifyCtrl(UINT8 Msg = 0);			// notify via control sockets that, default, there is at least 1 request available to be sent to a session instance provider
	bool SendCtrlMsg(int Len,				// number of control message bytes to send ptd at by pCtrlMsg, can be 0
								UINT8 *pCtrlMsg);	// pCtrlMsg pts to control message bytes
	int								// received message is this length, 0 if no outstanding messages, < 0 if errors
			RcvCtrl(int BuffLen,			// available buffer into which the control message can be copied 
		UINT8 *pBuff);			// copy control message into this buffer
	bool ProcessCtrlMsg(int MsgLen,UINT8 *pMsg);			// process a message received by control socket m_Ctrl[1] - note: currently these messages are simply discarded

	UINT32									// identifies created and initialised check point requests list
		CreateChkPtReqs(UINT32 SessionID,	// uniquely identifying this session between service requester and provider
				UINT32 TypeSessionID,		// used as an index within the type to reference this session
				UINT64 ClassInstanceID,		// check pointing this class instance
			    UINT32 ClassMethodID,		// identifies class method
				UINT32 ParamSize,			// instance specific parameter size
				UINT32 InDataSize,			// instance specific input data size
				UINT32 MaxOutDataSize,		// instance specific max expected output data for this method
				UINT8 *pData);				// request parameters followed by input data

	bool
			DeleteChkPtReqs(UINT32 ChkPtID);  // identifies check point list to be deleted

	bool			// true if request was successfully check pointed
		AddChkPtReq(UINT32 ChkPtID,  // identifies check point list to extend with this request
				UINT32 SessionID,			// uniquely identifying this session between service requester and provider
			    UINT32 ClassMethodID,		// identifies class method
				UINT32 ParamSize,			// instance specific parameter size
				UINT32 InDataSize,			// instance specific input data size
				UINT32 MaxOutDataSize,		// instance specific max expected output data for this method
				UINT8 *pData);				// request parameters followed by input data

	bool			// true if response was successfully updated		
		AddChkPtResp(UINT32 ChkPtID,  // identifies check point list to update with this response
					UINT32 SessionID,			// uniquely identifying this session between service requester and provider
					UINT32 JobRslt,				// job processing result as returned by service provider
					UINT32 RespSize,			// response size
					UINT8 *pData);				// response data

	bool ShutdownConnection(socket_t *pSocket);

	int			// returns 0 if no sockets to be monitored with select(), on windows the total number of monitored sockets, on linux the highest socket file descriptor plus 1
		SetupFDSets(fd_set& ReadFDs,			// select() read available socket descriptor set  
					fd_set& WriteFDs,			// select() write accepted socket descriptor set
					fd_set& ExceptFDs);			// select() exceptions descriptor set


	tJobIDEx										// packed job identifier or 0 if range errors
				PackIntoJobIDEx(UINT32 ReqID,							// must be in the range 1..16777215 (24bits)
									   UINT32 SessionID,				// service provider session identifier, must be in the range 1..131071 (17bits)
									   UINT32 InstanceID,				// service instance in the service providers session SessionID, must be in the range 1..511 (9bits)
									   UINT32 TypeID,					// identifies service being provided by the service provider, must be in the range 1..15 (4bits)
									   UINT32 TypeSessionID);			// index of this session in the service type, must be in the range 1..511 (9bits)

	bool				// false if any range errors whilst unpacking
			UnpackFromJobIDEx(tJobIDEx JobIDEx,								// unpack from this extended job identifier
										 UINT32 *pReqID,				    // returned ReqID, will be in the range 1..16777215
										 UINT32 *pSessionID,				// returned SessionID, will be in the range 1..131071
										 UINT32 *pInstanceID,				// returned InstanceID, will be in the range 1..511
										 UINT32 *pTypeID,					// returned service TypeID, will be in the range 1..15
										 UINT32 *pTypeSessionID);			// returned index of this session in the service type, will be in the range 1..511

#ifdef WIN32	
	const char *WSAGetLastErrorMessage(const char* pcMessagePrefix, int nErrorID = 0);
#endif
	
	int	Reset(bool bReleaseLock);			// locked mutex requires releasing immediately prior to deleting the mutexes

	

public:
	CBKSRequester();
	~CBKSRequester();

	int StartThread(tsRequesterThread *pPars);


						// NOTE: defaults chosen are representative for Smith-Waterman service processing
	int					// returns total number of registered service types or teBSFrsltCodes error code if any parameterisation errors or already registered type
		RegServiceType(teBKSPType BKSPType = eBKSPTSmithWaterman,							// registering this service type
									  UINT32 MinProviderVersion = cMinProviderVersion,		// service provider version must be at least this software version
									  UINT32 MaxProviderVersion = cMaxProviderVersion,		// service provider version must be no more than this software version
									  UINT32 KeepaliveSecs = cDfltKeepaliveSecs,			// expecting packet activity from session peer with periodicity of no more than this number of seconds
									  UINT32 MaxSessions = cDfltNumSessions,				// allowing at most this many concurrent sessions of specified service type
									  UINT32 MinServiceInsts = cMinServiceInsts,			// any session to support at least this minimum number of service instances
									  UINT32 MaxServiceInsts = cDfltServiceInsts,			// limit any session to support a maximum of this many service instances
									  UINT32 MaxProcSecs = cMaxSWProcSecs,                   // expecting any service provider to take no more than this number of seconds to complete processing a given request
									  UINT32 MaxParamLen = cMaxSWParamLen,					// service job parameters can be up to this length
  									  UINT32 MaxQuerySeqLen = cMaxSWQuerySeqLen,			// query sequences can be up to this length
									  UINT32 MaxTargSeqLen = cMaxSWTargSeqLen,				// target sequences can be up to this length
								      UINT32 MaxReqPayloadSize = cMaxSWReqPayloadSize,		// request payloads to the service provider, including framing, can be up to this size (UINT8s),
									  UINT32 MaxRespPayloadSize = cMaxSWRespPayloadSize);	// response payloads from the service provider, including framing, can be up to this  size (UINT8s)

	int						// returns number of registered sessions which are in requested state
		GetNumSessions(teBKSPType BKSPType = eBKSPTSmithWaterman,					// may request number of sessions providing this specific service type, or if eBKSPTUndefined then the total
					 teBKSPEPProvState BKSPState = eBKSPSRegisteredActv);			// sessions must be in this specific state

	int						// returns total number of service instances in all registered sessions which are in requested state
		GetNumInstances(teBKSPType BKSPType = eBKSPTSmithWaterman,					// may request number of instances providing this specific service type, or if eBKSPTUndefined then the total
					   teBKSPEPProvState BKSPState = eBKSPSRegisteredActv);			// sessions providing service instances must be in this specific state

	int													// actual number of sessions returned
		GetSessions(teBKSPType BKSPType,				// return sessions for this service type,  or if eBKSPTUndefined then all service types registered
					 teBKSPEPProvState BKSPState,		// sessions must be in this state
				int MaxSessions,						// return at most this many sessions in pSessions[]
				tsBKSRegSession *pSessions);			// caller has preallocated to hold returned array snapshot of sessions

	int 
			Initialise(char* pszHost = NULL,			// listening on this host/IP address; NULL to use first INET IP local to this machine
						char *pszService = NULL,			// listening on this service/port; NULL to use default port 
						teBKSPType BKSPType = eBKSPTSmithWaterman,					// registering this service type
									  UINT32 MinProviderVersion = cMinProviderVersion,		// service provider version must be at least this software version
									  UINT32 MaxProviderVersion = cMaxProviderVersion,		// service provider version must be no more than this software version
									  UINT32 KeepaliveSecs = cDfltKeepaliveSecs,			// expecting packet activity from session peer with periodicity of no more than this number of seconds
									  UINT32 MaxSessions = cDfltNumSessions,				// allowing at most this many concurrent sessions of specified service type
									  UINT32 MinServiceInsts = cMinServiceInsts,			// any session to support at least this minimum number of service instances
									  UINT32 MaxServiceInsts = cDfltServiceInsts,			// limit any session to support a maximum of this many service instances
									  UINT32 MaxProcSecs = cMaxSWProcSecs,                   // expecting any service provider to take no more than this number of seconds to complete processing a given request
									  UINT32 MaxParamLen = cMaxSWParamLen,					// service job parameters can be up to this length
  									  UINT32 MaxQuerySeqLen = cMaxSWQuerySeqLen,			// query sequences can be up to this length
									  UINT32 MaxTargSeqLen = cMaxSWTargSeqLen,				// target sequences can be up to this length
								      UINT32 MaxReqPayloadSize = cMaxSWReqPayloadSize,		// request payloads to the service provider, including framing, can be up to this size (UINT8s),
									  UINT32 MaxRespPayloadSize = cMaxSWRespPayloadSize);	// response payloads from the service provider, including framing, can be up to this  size (UINT8s)

	bool Run(int Secs=180);					// false if unable to begin thread execution or if executing thread took too long (default 3min) to initialise ready to accept connections

	void Terminate(int Secs = 120);			// allow this many seconds for thread to self-terminate before forcing termination

	int												// total number of classes 
		GetNumClassInstances(teBKSPType TypeID = eBKSPTSmithWaterman,		// service required
						UINT32 *pCommited = NULL,			// returned number of class instances currently committed or instantiated 
						UINT32 *pUncommited = NULL);    // returned number of class instances currently not committed and available to be instantiated

	int				// -2: parameter errors, -1: class instance no longer exists, 0: currently no available service instance 1: if job accepted
			AddJobRequest(  tJobIDEx *pJobID,					// returned unique job identifier by which job can later be referenced
							teBKSPType TypeID,					// service type required
							  UINT64 ClassInstanceID,			// class instance on which job method is to be applied
							  UINT32 ClassMethodID,				// class method to apply on the class instance
									 UINT32 ParamsSize = 0,		// processing parameters are this total size in bytes
									 void *pParams = NULL,		// service processing parameters
									 UINT32 InDataSize = 0,		// service processing input data is this total size in bytes
									 void *pInData = NULL,		// service processing input data
									 UINT32 SessionID = 0);     // if 0 then use session as specified by ClassInstanceID, otherwise use session corresponding to specific session identifier

	int				// < 0 if job no longer exists, 0 if job still being processed, > 0 if job completed
		GetJobResponse(tJobIDEx	JobID,			// unique job identifier returned when job was submitted
				  UINT64 *pClassInstanceID,		// returned class instance on which job method was applied
				  UINT32 *pClassMethodID,		// returned class method applied
					   UINT32 *pJobRslt,		// job processing result as returned by service provider
						UINT32 *pOutDataSize,	// (IN) service processing output results expected to be at most this total length, [OUT] bytes of response data copied into pOutData 
						void *pOutData,			// service processing output results data
						bool bRetain=false);		// true if job response is to be retained and not deleted; subsequent call with bRetain==false will delete this response


};

