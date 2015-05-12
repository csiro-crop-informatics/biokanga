// Copyright 2013 CSIRO  ( http://www.csiro.au/ ) 
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

#pragma once

const UINT32 cMaxEdges = 250;						// any vertex can have at most this many inbound or outbound edges
const UINT32 cDfltMinAcceptScoreThres = 5000;          // default is to only accept edges for processing which have scores of at least this threshold

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
const UINT32 cInitialAllocEdges =  (cInitialAllocVertices * 3);	// initially allocate for this many outgoing edges
const double cReallocEdges      =   0.3;		// then, as may be required, realloc in increments of this proportion pf existing outgoing edges
const UINT32 cInitalComponentsAlloc   = 2000000;	// hopefully there will not be too many assemblies with more than 2M components - will realloc as may be required
const double cReallocComponents  = 0.3;		// realloc as may be required in this proportion of existing components

const UINT32 cTransitStackAlloc   =    5000000;		// initially allocate for transition stack to hold this many entries
const double cReallocTransitStack =    0.3;		// realloc by this proportion of existing stack entries



const UINT32 cMaxDiscRemaps = 1000;					// remap disconnected graph identifiers list limit

typedef enum TAG_eVerticesSortOrder {
	eVSOUnsorted = 0,	//  unsorted or sort order indeterminate
	eVSOVertexID,		// sorted by vertex identifier ascending
	eVSOSeqID,			// sorted by sequence identifier ascending
	eVSOVertexComponentID // sorted by ComponentID ascending
	} eVerticesSortOrder;

typedef enum TAG_eOverlapClass {
	eOLCOverlapping = 0,	// probe classified as overlapping target, either 5' or 3'
	eOLCcontaining,			// probe completely contains the target
	eOLCartefact			// probe contains a subsequence of target, classifying as an artefact overlap and not further processed
} eOverlapClass;

#pragma pack(1)

// graph consists of vertices (representing sequences) and connecting edges (overlaying sequences) between adjacent vertices
// edges are of two types - forward edges represent overlaying sequences (verticeA is overlaying vertexB), and overlaid (verticeA is overlaid by verticeB)
// it is also important to note that the graph is likely to contain many - perhaps millions - of disconnected components 

// an instance of a Vertex which is used to represent a read sequence
typedef struct TAG_sGraphVertex {
	    tVertID VertexID;		// monotonically ascending (1..V) identifier uniquely identifies this vertex
		tEdgeID OutEdgeID;		// m_pGraphOutEdges[OutEdgeID-1] at which outgoing edges from this vertex start, 0 if none
		tEdgeID InEdgeID;		// m_pGraphInEdges[InEdgeID-1] identifying incoming edges to this vertex start, 0 if none
		tSeqID SeqID;			// identifies the sequence being represented by this vertex
		UINT32 SeqLen;			//  sequence is this length
		tComponentID ComponentID; // vertex is a member of this disconnected component
		UINT32 DegreeOut:8;		// number of outgoing edges from this vertex to any adjacent vertices, clamped to max cMaxEdges (currently 255)
		UINT32 DegreeIn:8;		// number of incoming edges from any adjacent vertices, clamped to max cMaxEdges (currently 255)
		UINT32 flgEmitted:1;	// sequence has been emitted as part of a fragment;	
		UINT32 flgRmvEdges:1;   // sequence may contain SMRTbell hairpin as this read aligns at least twice to the same other read
	} tsGraphVertex;


// an instance of an outgoing edge, from 'FromVertexID' to 'ToVertexID' 
// normally will be sorted in FromVertexID.ToVertexID ascending order
typedef struct TAG_sGraphOutEdge {
	tVertID FromVertexID;	// edge is outgoing (overlapping) from this vertex 
	tVertID ToVertexID;		// edge is incoming (overlapped) to this Vertex
	UINT32 FromSeq5Ofs;		// overlap of FromVertexID onto ToVertexID starts at this base relative to FromVertexID 5' start
	UINT32 FromSeq3Ofs;		// overlap of FromVertexID onto ToVertexID ends at this base relative to FromVertexID 5' start
	UINT32 ToSeq5Ofs;		// overlap onto ToVertexID from FromVertexID starts at this base relative to ToVertexID 5' start
	UINT32 ToSeq3Ofs;		// overlap onto ToVertexID from FromVertexID ends at this base relative to ToVertexID 5' start
	UINT32 Score;			// score associated with this overlap, 0..0x7fff, higher scores represent higher confidence in the overlap
	UINT32 OverlapClass:2;  // original classification (eOverlapClass)
	UINT32 Contains:1;		// set if FromVertexID completely contains ToVertexID
	UINT32 Artefact:1;		// set if FromVertexID overlap onto ToVertexID classified as being artefactual
	UINT32 OverlapSense:2;	// 0 FromVertexID sense overlaps ToVertexID sense, 1 FromVertexID antisense overlaps ToVertexID sense, 2 FromVertexID sense overlaps ToVertexID antisense
	UINT32 bRemove:1;		// set if this vertex marked for removal
	UINT32 TravFwd:1;		// set after this edge has been traversed from FromVertexID to ToVertexID
	UINT32 TravRev:1;		// set after this edge has been traversed from ToVertexID to FromVertexID
	UINT32 AuxFlgs:23;		// currently unassigned (ensure flags total to <= 32 so as to all fit within UINT32)
	} tsGraphOutEdge;

typedef struct TAG_sOverlappedSeq {
				tSeqID OverlappingSeqID;	// identifies the overlapping (FromVertexID) sequence
				tSeqID OverlappedSeqID;		// identifies overlapped (ToVertexID) sequence or a completely contained sequence
				UINT32 FromSeq5Ofs;			// overlap of FromVertexID onto ToVertexID starts at this base relative to FromVertexID 5' start
				UINT32 FromSeq3Ofs;			// overlap of FromVertexID onto ToVertexID ends at this base relative to FromVertexID 5' start
				UINT32 ToSeq5Ofs;			// overlap onto ToVertexID from FromVertexID starts at this base relative to ToVertexID 5' start
				UINT32 ToSeq3Ofs;			// overlap onto ToVertexID from FromVertexID ends at this base relative to ToVertexID 5' start
				UINT16 Score;				// score associated with this overlap, 0..0x7fff, higher scores represent higher confidence in the overlap
				eOverlapClass OverlapClass;	// classification of overlap from OverlappingSeqID onto OverlappedSeqID
				bool bAntisense;			// false: sense overlaps sense, true: antisense overlaps sense 
	} tsOverlappedSeq;

// disconnected graph components
typedef struct TAG_sComponent {
	tComponentID ComponentID;			// identifies this component
	tVertID VertexID;					// one of the vertices which is in this component
	UINT32 NumVertices;					// there are this many vertices in this component
	UINT32 PathNumVertices;				// there are this many vertices along path inclusive of PathStartVertexID and PathEndVertexID
	tVertID PathStartVertexID;          // path starts from this vertex
	tVertID PathEndVertexID;            // path ends at this vertex
} tsComponent;


typedef struct TAG_sRemapComponentID {
	tComponentID From;		// map from	
	tComponentID To;		// map to
} tsRemapComponentID;

#pragma pack()

class CAssembGraph
{
	CMTqsort m_MTqsort;				// multithreaded sorting

	bool m_bTerminate;				// if set true then all threads should terminate processing

	UINT32 m_MinScaffScoreThres;      // edges (overlaps) must be at least this score to be accepted for processing

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

	UINT32 m_NumComponents;				// number of components
	UINT32 m_AllocComponents;			// number of components allocated
	tsComponent *m_pComponents;			// allocated to hold array of identified components

	UINT32 m_CurTransitDepth;			// current transition depth
	UINT32 m_MaxTransitDepth;			// deepest transition required
	UINT32 m_AllocTransitStack;			// allocation is for this many transition entries 
	tVertID *m_pTransitStack;			// allocated to hold transition entries


	UINT32 m_NumDiscRemaps;				// number of disconnected graph identifiers requiring remaps
	tsRemapComponentID m_DiscRemaps[cMaxDiscRemaps];	// to hold disconnected graph identifiers requiring remaps
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
	static int SortVerticesComponentID(const void *arg1, const void *arg2);
	static int SortComponentNumVertices(const void *arg1, const void *arg2);
	static int SortComponentID(const void *arg1, const void *arg2);

	tsGraphOutEdge *					// ptr to edge with matching FromVertexID and ToVertexID, or NULL if unable to locate				
		LocateFromToEdge(tVertID FromVertexID,	// match edge with this FromVertexID which is
				  tVertID ToVertexID);			// to this ToVertexID

	tsGraphOutEdge *							// ptr to edge with matching FromVertexID and ToVertexID, or NULL if unable to locate				
		LocateToFromEdge(tVertID ToVertexID,	// match edge with this ToVertexID which is 
				  tVertID FromVertexID);		// from this FromVertexID

	UINT32			// index+1 in m_pFwdGraphEdges of first matching FromVertexID, or 0 if non matching				
		LocateFirstVertID(tVertID FromVertexID);	// find first matching 
	
	UINT32			// index+1 in m_pFwdGraphEdges of first matching DnSeqID, or 0 if non matching				
		LocateFirstDnSeqID(tSeqID DnSeqID);		// find first matching 


	tsGraphVertex *			// ptr to vertex corresponding to SeqID , or NULL if unable to locate				
			LocateVertexSeqID(tSeqID SeqID);		// match this SeqID

	tsGraphVertex *			// ptr to vertex corresponding to VertixID , or NULL if unable to locate				
			LocateVertex(tVertID VertixID);		// match this VertixID

	UINT32					// number of replacements
		ReplaceComponentID(tComponentID ToReplaceID,		// replace existing disconnected graph identifiers 
							tComponentID ReplaceID);		// with this identifier

	UINT32					// number of replacements
		ReplaceComponentIDs(void);

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
	UINT32  TransitIdentDiscGraph(tVertID VertexID,tComponentID ComponentID);
	UINT32  ClearEdgeTravFwdRevs(void);
	UINT32	ClearDiscCompIDs(void);

public:
	CAssembGraph(void);
	~CAssembGraph(void);

	void Reset(void);								// reset and free any allocated resources
	teBSFrsltCodes		// initialise with ScaffScoreThres an  maxThreads
		Init(int ScaffScoreThres = cDfltMinAcceptScoreThres,		// accepted edges must be of at least this overlap score
						int MaxThreads = 8);						// set number of threads

	UINT32								// returned vertex identifier
		AddVertex(UINT32 SeqLen,		// sequence length
				  tSeqID SeqID);		// allocates and initialises a new graph vertex which references this sequence

	UINT32								// 0 if errors else number of vertices finalised
		FinaliseVertices(void);			// added vertices are sorted ready for graph edges to be added

	UINT32								// number of edges removed	 
		ReduceEdges(void);				// reduce graph by detecting and removing extraneous edges


	UINT32									// returns total number of edges, including this edge if accepted, thus far accepted; if more than 
		AddEdge(tSeqID OverlappingSeqID,	// identifies the overlapping sequence
				tSeqID OverlappedSeqID,		// identifies overlapped sequence or a completely contained sequence
				UINT32 Score,				// score associated with this overlap, 0..0x7fff, higher scores represent higher confidence in the overlap
				UINT32 FromSeq5Ofs,			// overlap of FromVertexID onto ToVertexID starts at this base relative to FromVertexID 5' start
				UINT32 FromSeq3Ofs,			// overlap of FromVertexID onto ToVertexID ends at this base relative to FromVertexID 5' start
				UINT32 ToSeq5Ofs,			// overlap onto ToVertexID from FromVertexID starts at this base relative to ToVertexID 5' start
				UINT32 ToSeq3Ofs,			// overlap onto ToVertexID from FromVertexID ends at this base relative to ToVertexID 5' start
				eOverlapClass OverlapClass,	// classification of overlap from OverlappingSeqID onto OverlappedSeqID
				bool bAntisense);			// false: sense overlaps sense, true: antisense overlaps sense 

	UINT32									// returns total number of edges , including any of these edges, if accepted, thus far accepted 
		AddEdges(UINT32 NumSeqs,				// number of overlapped sequences
				tsOverlappedSeq *pOverlappedSeqs); // overlapped sequences

	UINT32									// returns total number of edges accepted after finalisation
		FinaliseEdges(void);
	
	UINT32 GetNumGraphVertices(void);		// returns current number of graph vertices

	UINT32 GetNumGraphOutEdges(void);		// returns current number of graph forward edges

	UINT32 GetNumReducts(void);				// returns current number of edge reductions

	UINT32 IdentifyDiscComponent(tVertID VertexID, // start component traversal from this vertex
				tComponentID ComponentID);			 // mark all traversed vertices as members of this component

	UINT32									// returned extended length of probe sequence if probe overlap onto target was accepted, 0 if unable to extend
	OverlapAcceptable(bool bSenseDir,		// false: downstream to 3' previous as sense, true: downstream to 3' previous as antisense
					  bool *pbSenseDir,		// returned sense direction to use
					tsGraphOutEdge *pEdge);

	UINT32															// number of vertices marked as members of this component
		FindMaxShortestPaths(tComponentID ComponentID);		 // mark all traversed vertices as members of this component

	UINT32	IdentifyDiscComponents(void);

	UINT64 WriteScaffoldSeqs(char *pszOutFile);  // write scaffolded PacBio sequences to this output file
};



