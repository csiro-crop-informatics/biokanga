/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once

#include "./Kangadna.h"

const UINT32 cMaxNumVertices = 0x080000000;			// graphs are limited to at most 2 billion vertices
const UINT32 cMaxNumEdges    = 0x0fffffffe;			// graphs are limited to at most just under 4 billion edges

#ifdef _DEBUG
#ifdef _WIN32
const UINT32 cInitialAllocVertices =10000000;		// initially allocate for this many vertices
#else
const UINT32 cInitialAllocVertices =50000000;		// initially allocate for this many vertices
#endif
#else
const UINT32 cInitialAllocVertices =50000000;		// initially allocate for this many vertices
#endif

const double cReallocVertices   =   0.3;	    // then, as may be required, realloc in increments of this proportion of existing vertices
const UINT32 cInitialAllocEdges =  (cInitialAllocVertices * 2);	// initially allocate for this many outgoing edges
const double cReallocEdges      =   0.3;		// then, as may be required, realloc in increments of this proportion pf existing outgoing edges
const UINT32 cInitalComponentsAlloc   = 2000000;			// hopefully there will not be too many assemblies with more than 2M components - will realloc as may be required
const double cReallocComponents  = 0.3;		// realloc as may be required in this proportion of existing components

const UINT32 cTransitStackAlloc   =    5000000;		// initially allocate for transition stack to hold this many entries
const double cReallocTransitStack =    0.3;		// realloc by this proportion of existing stack entries



const UINT32 cMaxDiscRemaps = 1000;					// remap disconnected graph identifiers list limit

typedef enum TAG_eVerticesSortOrder {
	eVSOUnsorted = 0,	//  unsorted or sort order indeterminate
	eVSOVertexID,		// sorted by vertex identifier ascending
	eVSOSeqID,			// sorted by sequence identifier ascending
	eVSOVertexDiscGraphID // sorted by DiscGraphID ascending
	} eVerticesSortOrder;

#pragma pack(1)

typedef UINT32 tDiscGraphID;	// to contain disconnected subgraph identifiers, 1..D


// graph consists of vertices (representing sequences) and connecting edges (overlaying sequences) between adjacent vertices
// edges are of two types - forward edges represent overlaying sequences (verticeA is overlaying vertexB), and overlaid (verticeA is overlaid by verticeB)
// it is also important to note that the graph is likely to contain many - perhaps millions - of disconnected components 

// an instance of a Vertex which is used to represent a read sequence
// currently requires 32 bytes per vertex instance
typedef struct TAG_sGraphVertex {
	    tVertID VertexID;		// monotonically ascending (1..V) identifies this vertex
	    tVertID PEVertexID;	// if non-zero then this is the partner (paired ends) vertex identifier
		tEdgeID OutEdgeID;		// m_pGraphOutEdges[OutEdgeID-1] at which outgoing edges from this vertex start, 0 if none
		tEdgeID InEdgeID;		// m_pGraphInEdges[InEdgeID-1] identifying incoming edges to this vertex start, 0 if none
		tSeqID SeqID;			// identifies the sequence being represented by this vertex
		UINT32 SeqLen;			//  sequence is this length
		tDiscGraphID DiscGraphID; // vertex is a member of this disconnected component
		UINT32 DegreeOut:4;		// number of outgoing edges from this vertex to any adjacent vertices, clamped to max 15
		UINT32 DegreeIn:4;		// number of incoming edges from any adjacent vertices, clamped to max 15
		UINT32 flgEmitted:1;	// sequence has been emitted as part of a fragment;	
	} tsGraphVertex;


// an instance of an outgoing edge, from 'FromVertexID' to 'ToVertexID' 
// currently requires 16 bytes per edge instance
// will be sorted in FromVertexID.ToVertexID ascending order
typedef struct TAG_sGraphOutEdge {
	tVertID FromVertexID;	// edge is outgoing (overlapping) from this vertex 
	tVertID ToVertexID;	// edge is incoming (overlapped) to this Vertex
	UINT16 SeqOfs;			// overlap of FromVertexID onto ToVertexID starts at this base relative to FromVertexID 5' start
	UINT16 bSeqIDisPE2:1;	// if set then FromVertexID is PE2 of a paired end, reset if single ended or PE1 of paired end
	UINT16 bDnSeqIDisPE2:1;	// if set then FromVertexID is PE2 of a paired end, reset if single ended or PE1 of paired end
	UINT16 Contains:1;		// set if FromVertexID completely contains ToVertexID
	UINT16 OverlapSense:2;	// 0 FromVertexID sense overlaps ToVertexID sense, 1 FromVertexID antisense overlaps ToVertexID sense, 2 FromVertexID sense overlaps ToVertexID antisense
	UINT16 bRemove:1;		// set if this vertex marked for removal
	UINT16 TravFwd:1;		// set after this edge has been traversed from FromVertexID to ToVertexID
	UINT16 TravRev:1;		// set after this edge has been traversed from ToVertexID to FromVertexID
	UINT16 AuxFlgs:8;		// currently unassigned (ensure flags total to 16 so as to all fit within 16bits)
	} tsGraphOutEdge;

typedef struct TAG_sOverlappedSeq {
				tSeqID OverlappingSeqID;	// identifies the overlapping sequence
				tSeqID OverlappedSeqID;		// identifies overlapped sequence or a completely contained sequence
				int SeqOfs;					// overlap (0..n) starts at this base relative to OverlappingSeqID 5'
				bool bContains;				// true if SeqID1 completely contains DnSeqID
				int OverlapSense;			// 0 sense overlaps sense, 1 antisense overlaps sense, 2 sense overlaps antisense 
	} tsOverlappedSeq;

// disconnected graph components
typedef struct TAG_sComponent {
	tDiscGraphID ComponentID;			// identifies this component
	tVertID VertexID;					// one of the vertices which is in this component
	UINT32 NumVertices;					// there are this many vertices in this component
} tsComponent;


typedef struct TAG_sRemapDiscGraphID {
	tDiscGraphID From;		// map from	
	tDiscGraphID To;		// map to
} tsRemapDiscGraphID;

#pragma pack()

class CAssembGraph
{
	CMTqsort m_MTqsort;				// multithreaded sorting

	bool m_bTerminate;				// if set true then all threads should terminate processing
	eVerticesSortOrder m_VerticesSortOrder; // current graph vertex sort order
	bool m_bVertexEdgeSet;				// true if InEdgeID/OutEdgeID have been initialised
	UINT32 m_UsedGraphVertices;			// number of graph vertices currently used
	UINT32 m_AllocGraphVertices;		// number of graph vertices allocated
	tsGraphVertex *m_pGraphVertices;   // allocated to hold array of graph vertices

	bool m_bOutEdgeSorted;				// true if m_pGraphOutEdges has been sorted in ascending FromVertexID.ToVertexOrder
	UINT32 m_UsedGraphOutEdges;			// number of forward graph edges currently used
	UINT32 m_AllocGraphOutEdges;		// number of forward graph edges allocated
	tsGraphOutEdge *m_pGraphOutEdges;	// allocated to hold array of forward edges

	bool m_bInEdgeSorted;				// true if m_pGraphInEdges has been sorted in ascending ToVertexOrder.FromVertexID
	UINT32 m_UsedGraphInEdges;			// number of inbound graph edges currently used
	UINT32 m_AllocGraphInEdges;			// number of inbound graph edges allocated
	tEdgeID *m_pGraphInEdges;			// index onto m_pGraphOutEdges which is sorted in ToVertexID.FwdVertexID ascending order

	UINT32 m_UsedComponents;			// number of components
	UINT32 m_AllocComponents;			// number of components allocated
	tsComponent *m_pComponents;			// allocated to hold array of identified components

	UINT32 m_CurTransitDepth;			// current transition depth
	UINT32 m_MaxTransitDepth;			// deepest transition required
	UINT32 m_AllocTransitStack;			// allocation is for this many transition entries 
	tVertID *m_pTransitStack;			// allocated to hold transition entries


	UINT32 m_NumDiscRemaps;			// number of disconnected graph identifiers requiring remaps
	tsRemapDiscGraphID m_DiscRemaps[cMaxDiscRemaps];	// to hold disconnected graph identifiers requiring remaps
	UINT32 m_NumReducts;			// number of graph node reductions
	bool m_bReduceEdges;			// if true then attempt to reduce extraneous graph edges when current graph is full and more node memory needs to be allocated  

	void AcquireSerialise(void);
	void ReleaseSerialise(void);
	void AcquireLock(bool bExclusive);				// defaults as read only lock
	void ReleaseLock(bool bExclusive);

	int m_NumThreads;				// use at most this number threads for graph processing
	bool m_bMutexesCreated;			// set true if mutexes and rwlocks created/initialised
#ifdef _WIN32
	CRITICAL_SECTION m_hSCritSect;
	HANDLE m_hMtxMHReads;
	SRWLOCK m_hRwLock;
#else
	pthread_spinlock_t m_hSpinLock;
	pthread_mutex_t m_hMtxMHReads;
	pthread_rwlock_t m_hRwLock;
#endif
	static int SortOutEdgeFromVertexID(const void *arg1, const void *arg2);
	static int SortOutEdgeFromVertexIDSeqOfs(const void *arg1, const void *arg2);
	static int SortInEdgesToVertexID(const void *arg1, const void *arg2);
	static int SortOutEdgeToVertexIDSeqOfs(const void *arg1, const void *arg2);

	static int SortVerticesSeqID(const void *arg1, const void *arg2);
	static int SortVerticesVertexID(const void *arg1, const void *arg2);
	static int SortVerticesDiscGraphID(const void *arg1, const void *arg2);
	static int SortComponentNumVertices(const void *arg1, const void *arg2);

	tsGraphOutEdge *						// ptr to matching graph node, or NULL if unable to locate				
		LocateSeqIDDnSeqID(tSeqID SeqID,	// match this SeqID
				  tSeqID DnSeqID);			// with this DnSeqID

	tsGraphOutEdge *						// ptr to matching graph node, or NULL if unable to locate				
		LocateDnSeqIDSeq(tSeqID DnSeqID,	// match this DnSeqID
				  tSeqID SeqID);			// with this SeqID

	UINT32			// index+1 in m_pFwdGraphEdges of first matching FromVertexID, or 0 if non matching				
		LocateFirstVertID(tVertID FromVertexID);	// find first matching 
	
	UINT32			// index+1 in m_pFwdGraphEdges of first matching DnSeqID, or 0 if non matching				
		LocateFirstDnSeqID(tSeqID DnSeqID);		// find first matching 


	tsGraphVertex *			// ptr to vertex corresponding to SeqID , or NULL if unable to locate				
			LocateVertexSeqID(tSeqID SeqID);		// match this SeqID

	tsGraphVertex *			// ptr to vertex corresponding to VertixID , or NULL if unable to locate				
			LocateVertex(tVertID VertixID);		// match this VertixID

	UINT32					// number of replacements
		ReplaceDiscGraphID(tDiscGraphID ToReplaceID,		// replace existing disconnected graph identifiers 
							tDiscGraphID ReplaceID);		// with this identifier

	UINT32					// number of replacements
		ReplaceDiscGraphIDs(void);

	int CreateMutexes(void);
	void DeleteMutexes(void);

	int			// stack depth or < 1 if errors
		PushTransitStack(tVertID VertexID);		// push VertexID onto stack 
	tVertID PopTransitStack(void);				// popped VertexID or 0 if stack was empty
	tVertID PeekTransitStack(void);				// peeked VertexID or 0 if stack empty
	void ClearTransitStack(void);					// remove all entries from stack

	UINT32							// number of vertices with both inbound and outboud edges
		VertexConnections(void);		// identify and mark vertices which have multiple inbound edges
	UINT32 GenSeqFragment(tsGraphVertex *pVertex);		// initial seed vectex
	UINT32  TransitIdentDiscGraph(tVertID VertexID,tDiscGraphID DiscGraphID);
	UINT32  ClearEdgeTravFwdRevs(void);
	UINT32	ClearDiscCompIDs(void);

public:
	CAssembGraph(void);
	~CAssembGraph(void);

	void Reset(void);				// reset and free any allocated resources
	teBSFrsltCodes Init(void);		// initialise and allocate memory for graph nodes

	teBSFrsltCodes SetNumThreads(int maxThreads = 4);	// use 4 threads as the default 

	UINT32								// returned vertex identifier
		AddVertex(UINT32 SeqLen,		// sequence length
				  tSeqID SeqID,			// allocates and initialises a new graph vertex which references this sequence
				  tSeqID PESeqID=0);	// if paired end assembly then the paired (partner) sequence identifier, 0 if no paired end

	UINT32								// 0 if errors else number of vertices finalised
		FinaliseVertices(void);			// added vertices are sorted ready for graph edges to be added

	UINT32								// number of edges removed	 
		ReduceEdges(void);				// reduce graph by detecting and removing extraneous edges


	teBSFrsltCodes	 
		AddEdge(tSeqID OverlappingSeqID,	// identifies the overlapping sequence
				tSeqID OverlappedSeqID,		// identifies overlapped sequence or a completely contained sequence
				int SeqOfs,					// overlap (0..n) starts at this base relative to OverlappingSeqID 5'
				bool bContains,				// true if SeqID1 completely contains DnSeqID
				int OverlapSense);			// 0 sense overlaps sense, 1 antisense overlaps sense, 2 sense overlaps antisense 

	teBSFrsltCodes	 
		AddEdges(int NumSeqs,				// number of overlapped sequences
				tsOverlappedSeq *pOverlappedSeqs); // overlapped sequences

	UINT32									// 0 if errors else number of edges finalised
		FinaliseEdges(void);
	
	teBSFrsltCodes	 
		AddGraphSENode(tSeqID SeqID,		// this is the overlapping sequence
					bool bSeqIDisPE2,		// if true then SeqID is PE2 of a paired end, false if single ended or PE1 of paired end 
					tSeqID DnSeqID,			// which overlaps 5' end of this downstream sequence or contains this sequence
					bool bDnSeqIDisPE2,		// if true then DnSeqID is PE2 of a paired end, false if single ended or PE1 of paired end 
					int SeqOfs,				// overlap (0..n) starts at this base relative to SeqID 5'
					bool bContains,			// true if SeqID completely contains DnSeqID
					bool bAntisense);		// true if SeqID is antisense to DnSeqID 


	teBSFrsltCodes	 
		AddGraphPENode(tSeqID SeqID,		// this is the overlapping sequence
					bool bSeqIDisPE2,		// if true then SeqID is PE2 of a paired end, false if single ended or PE1 of paired end 
					tSeqID DnSeqID,			// which overlaps 5' end of this sequence or contains this sequence
					bool bDnSeqIDisPE2,		// if true then DnSeqID is PE2 of a paired end, false if single ended or PE1 of paired end 
					int SeqOfs,				// overlap (0..n) starts at this base relative to SeqID 5'
					bool bContains,			// true if SeqID completely contains DnSeqID
					bool bAntisense);		// true if SeqID is antisense to DnSeqID 
	
	UINT32 GetNumGraphVertices(void);		// returns current number of graph vertices

	UINT32 GetNumGraphOutEdges(void);		// returns current number of graph forward edges

	UINT32 GetNumReducts(void);				// returns current number of edge reductions

	UINT32 IdentifyDiscComponent(tVertID VertexID, // start component traversal from this vertex
				tDiscGraphID DiscGraphID);			 // mark all traversed vertices as members of this component

	UINT32	IdentifyDisconnectedSubGraphs(void);

	UINT64 WriteContigs(char *pszOutFile);  // write assembled contigs to this output file
};



