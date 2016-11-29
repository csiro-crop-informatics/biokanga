#pragma once

#ifndef cMaxConcurrentSessions
#define cMaxConcurrentSessions 100			// on windows, by default a max of 64 sockets are supported in FD_SET/FD_CLR and select(), 2 covers the listening and control sockets
#endif	

#ifdef _WIN32
#define socket_t SOCKET
#else
#define socket_t int
#endif

const UINT32 cMinProviderVersion = 1;			// service provider versions must be at least this software version
const UINT32 cMaxProviderVersion = 1;			// service provider version must be no more than this software version

const UINT32 cMaxHostNameLen = 80;			   // host names will be truncated to this maximal length
const UINT32 cMaxServiceNameLen = 80;		   // service names will be truncated to this maximal length


const UINT32 cMaxServiceTypes = 15;				// supporting at most this many service types (must fit within 4bits - see tJobIDEx)

const UINT32 cMinKeepaliveSecs = 15;			// session keep alive activity can be requested down to this absolute minimum (seconds)
const UINT32 cDfltKeepaliveSecs = 60;			// expecting packet activity at no more than this number of seconds since previous packet received
const UINT32 cMaxKeepaliveSecs = 600;			// session keep alive activity can be requested up to this absolute maximum (seconds)

const UINT32 cMaxNumSessions = cMaxConcurrentSessions;	// support a maximum of this many concurrent service provider sessions for a given service type (must fit within 9bits, see tJobEx)
                                                // NOTE: the total number of all sessions, independent of service type, must also total no more than cMaxConcurrentSessions  
const UINT32 cDfltNumSessions = 25;				// the default is support for at most this many concurrent sessions  a given service type (note that a single session can support multiple service instances of the same type)
const UINT32 cDfltMinSessions = 1;				// must be support for at least this many concurrent sessions for a given service type 

const UINT32 cMaxServiceInsts = 128;			// can only support at most this maximum number of service instances by single service provider (must fit within 9bits)
const UINT32 cMaxClassInsts = cMaxServiceInsts;	// any service provider can support up to this many class instances in a session

const UINT32 cDfltServiceInsts = 32;			// the default is support for at most this maximum number of service instances for any single service provider
const UINT32 cMinServiceInsts = 1;				// must be at least this minimum number of service instances for any single service provider

const UINT32 cMaxSessEstab = 20;                // allowing at most this many endpoints to be in establishment negotiation state
const UINT32 cMaxSessEstabTime = 60;            // terminate establishment negotiations if taking longer than this number of seconds

// maximal packet payload sizes depend on the actual service
// irrespective of service the following are the absolute maximum packet payload sizes supported
const UINT32 cAbsMaxReqPayloadSize = 0x03ffffff;	// request payloads to the service provider, including framing, can be up to this size (UINT8s) 
const UINT32 cAbsMaxRespPayloadSize = 0x03ffffff;	// response payloads from the service provider, including framing, can be up to this size (UINT8s) 
const UINT32 cAbsMaxQuerySeqLen = 0x01ffffff;       // maximum length of any query sequence
const UINT32 cAbsMaxTargSeqLen = 0x01ffffff;		// maximum length of any target sequence
const UINT32 cAbsMaxParamLen = 50000;				// maximum length of any parameterisation set
const UINT32 cAbsMinProcSecs = 120;					// allowing at least this many seconds for a service providers to process and return response for any individual request
const UINT32 cAbsMaxProcSecs = 36000;				// service providers expected to take less than this number of seconds to process and return response for any individual request


// following maximal sizes are for the eBKSPTSmithWaterman service providers
const UINT32 cMaxSWQuerySeqLen = 0x03ffff;        // maximum length of any query sequence
const UINT32 cMaxSWTargSeqLen  = 0x03ffff;		  // maximum length of any target sequence
const UINT32 cTypicalSWSeqLen  = 0x04fff;         // expecting the mean length of the PacBio reads to be typically around 20Kbp
const UINT32 cMaxSWParamLen = 0x03fff;            // maximum length of any SW parameterisation set
const UINT32 cMaxSWReqPayloadSize = (cMaxSWQuerySeqLen + cMaxSWTargSeqLen + cMaxSWParamLen + 0x0fff);	// request payloads to the service provider, including framing, can be up to this size (UINT8s) 
const UINT32 cMaxSWMAFBuffSize = (cTypicalSWSeqLen * 200 * 10); // allowing for MAFs to be returned with up to 200x coverage, 10x to allow for InDels, over typically sized reads 
const UINT32 cMaxSWRespPayloadSize = (cMaxSWMAFBuffSize + 0x0fff);	// maximal length response payloads from the service provider, including framing

const UINT32 cMaxSWProcSecs = 1200;				  // service providers expected to take less than this number of seconds to process and return response for any individual request


// Endpoint service providers are classified as providing one of the following exclusive service types
// NOTE: at most cMaxServiceTypes will be supported
typedef enum TAG_eBKSPType {
	eBKSPTUndefined = 0,					// service type is undefined
	eBKSPTSmithWaterman,					// provides Smith-Waterman alignments
	eBKSPTEcho,								// echo service
	eBKSPTPlaceHolder
	} teBKSPType;

typedef enum TAG_eBKSHdrType {
	eBKSHdrnone = 0,					// header type is not defined
	eBKSHdrKeepalive,					// timeout keepalive, sent in lieu of a payload carrying packet at least every service dependent time interval 
	eBKSHdrReqServices,					// negotiation phase A, notify potential service provider of service types required by server
	eBKSHdrOfferedService,				// negotiation phase B response back by potential provider with service type and number of service instances supported
	eBKSHdrAcceptService,				// negotiation phase C notify potential provider that offered service type has been accepted
	eBKSHdrRejectService,				// negotiation phase C notify potential provider that offered service type has been rejected
	eBKSHdrReq,							// service request to endpoint service provider
	eBKSHdrResp,						// service response from endpoint service provider
	eBKSHdrTermService,					// service is being terminated
	eBKSHdrPlaceHolder
	} teBKSHdrType;



#pragma pack(1)
// packet/frame structures exchanged between server as requester and endpoint service providers

// all received and sent payload data over sockets are preceded by a fixed size common header
// once the header has been received then the frame identifier (used to check for out of order and dropped frames) and following payload size is known
typedef struct TAG_sBKSPacHdr {
	UINT32 FrameLen;				// total length of this frame including this header and any following payload 
	UINT32 SessionID;				// server assigned and is uniquely identifying this session between service requester and provider, will be in the range 1..131071
	UINT8 TxFrameID;				// senders frame header identifier, monotonically incremented starting from 1 to 127 with wraparound back to 1
	UINT8 RxFrameID;				// senders last received and processed frame header identifier from current session peer, 0 if yet to receive any
	UINT8 FrameFlags;				// frame header flags, currently unused
	UINT8 FrameType;				// one of teBKSHdrType header types, specifies the payload type
} tsBKSPacHdr;
 

typedef struct TAG_sServiceDetail {
	teBKSPType BKSPType;				// service type
	UINT32 MinProviderVersion;			// service provider version must be of at least this software version
	UINT32 MaxProviderVersion;			// service provider version must be no more than this software version
	UINT32 KeepaliveSecs;				// expecting packet activity from endpoint with periodicity of no more than this number of seconds
	UINT32 MinServiceInsts;				// any endpoint to support at least this minimum number of service instances
	UINT32 MaxServiceInsts;				// limit any endpoint to support a maximum of this many service instances
	UINT32 MinClassInstances;			// any endpoint to support a minimum of this many class instances
	UINT32 MaxClassInstances;			// limit any endpoint to support a maximum of this many class instances
	UINT32 MaxProcSecs;                 // service provider expected to return a response to service request within this many seconds
	UINT32 MaxParamLen;					// service job parameters can total up to this length
	UINT32 MaxQuerySeqLen;				// query sequences can be up to this length
	UINT32 MaxTargSeqLen;				// target sequences can be up to this length
	UINT32 MaxReqPayloadSize;			// request payloads to the service provider, including framing, can be up to this size (UINT8s)
	UINT32 MaxRespPayloadSize;			// response payloads from the service provider, including framing, can be up to this  size (UINT8s)
	UINT32 Priority;				    // service provider should attempt to supply highest priority service type supported by that provider
} tsServiceDetail;

// first negotiation packet (eBKSHdrReqServices) sent by server to potential service provider with details of required services
typedef struct TAG_sBKSReqServices
	{
	tsBKSPacHdr Hdr;					// frame header (eBKSHdrReqServices)
	UINT32 NumTypes;					// server requires one of the following services; if 0 then server not requiring any services
	tsServiceDetail Details[1];			// array of service type details, one of which providers can offer to provide
	} tsBKSReqServices;

// response packet (eBKSHdrOfferedService) sent by potential service provider accepting supply of this service type and provider software of the specified version
typedef struct TAG_sBKSOfferedService
	{
	tsBKSPacHdr Hdr;					// frame header (eBKSHdrOfferedService)
	teBKSPType BKSPType;				// service type provided (eBKSPTUndefined if unable to supply any of requested services)
	UINT32 ServiceInsts;                // and is prepared to provide this many instances
	UINT32 ClassInsts;					// with this many instantiated class instances
	UINT32 Costing;						// provision of instances has an associated cost - provider may be a low resourced node so cost may be high 
	UINT32 ProviderVersion;				// service provider software is at this version
	} tsBKSOfferedService;

// in final negotiation phase server sends acceptance packet (eBKSHdrAcceptService) or rejection packet (eBKSHdrRejectService) 
typedef struct TAG_sBKSAcceptService
	{
		tsBKSPacHdr Hdr;				// frame header (eBKSHdrAcceptService or eBKSHdrRejectService)
		teBKSPType BKSPType;			// confirmation of service type being accepted, or eBKSPTUndefined if offered type rejected
		} tsBKSAcceptService;


// after service provision has been negotiated then service requests (eBKSHdrReq) are sent by server and responses (eBKSHdrResp) returned by service providers
// service request (eBKSHdrReq) sent by requesting server
typedef struct TAG_sBKSServReq {
	tsBKSPacHdr Hdr;				// frame header (eBKSHdrReq)
	INT64 JobIDEx;					// identifies this request and must be returned in the response
	UINT64 ClassInstanceID;			// class instance request applies to
	UINT32 ClassMethodID;			// identifies class method
	UINT32 ParamSize;				// parameter block is of this size in bytes
	UINT32 DataSize;				// data block is of this size in bytes
	UINT8 ParamData[1];				// parameters followed by request data 
	} sBKSServReq;

// service response (eBKSHdrResp) to a request sent by service provider
typedef struct TAG_sBKSServResp
{
	tsBKSPacHdr Hdr;				// frame header (eBKSHdrResp)
	INT64 JobIDEx;					// identifies this response and was copied from the request
	UINT64 ClassInstanceID;			// class instance response was from
	UINT32 ClassMethodID;			// identifies class method
	UINT32 JobRslt;					// service processing result
	UINT32 DataSize;				// response data block is of this size in bytes
	UINT8 Data[1];					// response data 
} sBKSServResp;

#pragma pack()
