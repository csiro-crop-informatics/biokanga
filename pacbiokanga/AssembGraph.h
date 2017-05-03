/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once

const UINT32 cMaxEdges = 250;						// limit any vertex to have at most this many inbound or outbound edges

													// Note: normalised overlap edge scores are independently calculated using 1 match, -1 mismatch, -3 gap open and -1 gap extension
const int cMinMin1kScore = 900;						 // when processing the multialignment overlaps for graph edges then only accept edge if normalised overlap score at least this per 1Kbp of overlap
const int cDfltMin1kScore = 980;					 // when processing the multialignment overlaps for graph edges then only accept edge if normalised overlap score at least this per 1Kbp of overlap
const int cMaxMin1kScore = 1000;					 // when processing the multialignment overlaps for graph edges then perfect match required

const UINT32 cInitialAllocVertices =500000;		// initially allocate for this many vertices, will be realloc'd if required
const double cReallocVertices   =   0.3;	    // then, as may be required, realloc in increments of this proportion of existing vertices
const UINT32 cInitialAllocEdges =  (cInitialAllocVertices * 10);	// initially allocate for this many outgoing edges
const double cReallocEdges      =   0.3;		// then, as may be required, realloc in increments of this proportion of existing outgoing edges
const UINT32 cInitalComponentsAlloc   = 50000;	// not expecting too many assemblies with more than 50K components but will realloc if required
const double cReallocComponents  = 0.3;			// realloc as may be required in this proportion of existing components

const UINT64 cMaxGraphEdges = 0x07fffffff;		// allowing at most this many edges in graph

const UINT32 cTransitStackAlloc   =    1000000;		// initially allocate for transition stack to hold this many entries
const double cReallocTransitStack =    0.3;		// realloc by this proportion of existing stack entries

const UINT32 cInitialAllocTraceBacks=100000;	// initially allocate for this many vertices, will be realloc'd if required
const double cReallocTraceBacks   =   0.3;	    // then, as may be required, realloc in increments of this proportion of existing tracebacks

const UINT32 cMaxDiscRemaps = 1000;					// remap disconnected graph identifiers list limit

typedef enum TAG_eVerticesSortOrder {
	eVSOUnsorted = 0,	//  unsorted or sort order indeterminate
	eVSOVertexID,		// sorted by vertex identifier ascending
	eVSOSeqID,			// sorted by sequence identifier ascending
	eVSOVertexComponentID // sorted by ComponentID ascending
	} eVerticesSortOrder;

typedef enum TAG_eOverlapClass {
	eOLCOverlapping = 0,	// probe classified as overlapping target, either 5' or 3'
	eOLCcontains,			// probe completely contains the target
	eOLCcontained,			// probe is completely contained within the target
	eOLCartefact			// probe contains a subsequence of target, classifying as an artefact overlap and not further processed
} eOverlapClass;

#pragma pack(1)

// graph consists of vertices (representing sequences) and connecting edges (overlaying sequences) between adjacent vertices
// edges are of two types - forward edges represent overlaying sequences (vertexA is overlaying vertexB), and overlaid (vertexA is overlaid by vertexB)
// it is also important to note that the graph is likely to contain many - perhaps millions - of disconnected components 

// an instance of a Vertex which is used to represent a read sequence
typedef struct TAG_sGraphVertex {
	    tVertID VertexID;		// monotonically ascending (1..V) identifier uniquely identifies this vertex
		tEdgeID OutEdgeID;		// m_pGraphOutEdges[OutEdgeID-1] at which outgoing edges from this vertex start, 0 if none
		tEdgeID InEdgeID;		// m_pGraphInEdges[InEdgeID-1] identifying incoming edges to this vertex start, 0 if none
		tSeqID SeqID;			// identifies the sequence being represented by this vertex
		UINT32 SeqLen;			//  sequence is this length
		UINT32 RecurseDepth;		// recurse depth at which this vertex is processed - used to determine circular reference
		tComponentID ComponentID; // vertex is a member of this disconnected component
		UINT8 DegreeOut;		// number of outgoing edges from this vertex to any adjacent vertices, clamped to max cMaxEdges (currently 250)
		UINT8 DegreeIn;			// number of incoming edges from any adjacent vertices, clamped to max cMaxEdges (currently 250)
		UINT8 flgEmitted:1;		// sequence has been emitted as part of a fragment;	
		UINT8 flgRmvEdges:1;	// sequence may contain SMRTbell hairpin as this read aligns at least twice to the same other read
		UINT8 flgPathAccepted:1;	// this vertex has been accepted as part of a highest scoring path
		UINT8 flgPathScored:1;		// set if highest scoring path already processed
		UINT8 flgPathTerm:1;		// vertex classified as path terminating
		UINT64 PathScore;			// highest score for any path originating from this vertex
		tEdgeID PathScoreEdgeID;	// highest scoring path starts with this outgoing edge
	} tsGraphVertex;


// an instance of an outgoing edge, from 'FromVertexID' to 'ToVertexID' 
// normally will be sorted in FromVertexID.ToVertexID ascending order
typedef struct TAG_sGraphOutEdge {
	tVertID FromVertexID;	// edge is outgoing from this vertex 
	tVertID ToVertexID;		// edge is incoming to this Vertex
	UINT32 FromSeqLen;		// 'From' vertex sequence length is this length - more efficent to duplicate copy of length with edges as it saves a lookup in graph vertices when length is needed
	UINT32 ToSeqLen;		// 'To' vertex sequence length is this length - more efficent to duplicate copy of length with edges as it saves a lookup in graph vertices when length is needed
	UINT32 FromSeq5Ofs;		// overlap of FromVertexID onto ToVertexID starts at this base relative to FromVertexID 5' start
	UINT32 FromSeq3Ofs;		// overlap of FromVertexID onto ToVertexID ends at this base relative to FromVertexID 5' start
	UINT32 ToSeq5Ofs;		// overlap onto ToVertexID from FromVertexID starts at this base relative to ToVertexID 5' start
	UINT32 ToSeq3Ofs;		// overlap onto ToVertexID from FromVertexID ends at this base relative to ToVertexID 5' start
	UINT32 Score;			// score associated with this overlap, higher scores represent higher confidence in the overlap
	UINT32 ScoreAlignLen;	// scored for this alignment length: (1 + ProbeAlignLen + TargAlignLen) / 2
	UINT16 flgFromAntisense:1;// set if From was revcpl as antisense probe 
	UINT16 flgToAntisense:1;// set if To was revcpl as antisense target 
	UINT16 flgRemove:1;		// set if this edge marked for removal
	UINT16 flgTravFwd:1;	// set after this edge has been traversed from FromVertexID to ToVertexID
	UINT16 flgTravRev:1;	// set after this edge has been traversed from ToVertexID to FromVertexID
	UINT16 flgInfBackEdge:1;// set if this is an inferred back edge, if 0 then edge was result of actual alignment
	} tsGraphOutEdge;

typedef struct TAG_sOverlappedSeq {
				tSeqID FromSeqID;			// identifies the 'From' overlapping (FromVertexID) sequence
				tSeqID ToSeqID;			   // identifies the 'To' overlapped (ToVertexID) sequence or a completely contained sequence
				UINT32 FromSeqLen;			// 'From' sequence length is this length
				UINT32 ToSeqLen;			// 'To' sequence length is this length
				UINT32 FromSeq5Ofs;			// overlap of FromVertexID onto ToVertexID starts at this base relative to FromVertexID 5' start
				UINT32 FromSeq3Ofs;			// overlap of FromVertexID onto ToVertexID ends at this base relative to FromVertexID 5' start
				UINT32 ToSeq5Ofs;			// overlap onto ToVertexID from FromVertexID starts at this base relative to ToVertexID 5' start
				UINT32 ToSeq3Ofs;			// overlap onto ToVertexID from FromVertexID ends at this base relative to ToVertexID 5' start
				UINT32 Score;				// score associated with this overlap, higher scores represent higher confidence in the overlap
				UINT32 ScoreAlignLen;		// scored for this alignment length: (1 + ProbeAlignLen + TargAlignLen) / 2
				UINT8 flgFromAntisense:1;	// set if From was revcpl as antisense probe 
				UINT8 flgToAntisense:1;		// set if To was revcpl as antisense target 
				eOverlapClass OverlapClass;	// classification of overlap from OverlappingSeqID onto OverlappedSeqID
	} tsOverlappedSeq;

const int cAllocPathEdges = 10000;			// initally allocate for this many path edges per component; realloc as needed

typedef struct TAG_sPathTraceBack {
				 tComponentID ComponentID;      // path is for this component
				 tVertID VertexID;				// path includes this vertex sequence
				 UINT32 SeqLen;					// vertex sequence length
				 UINT32 PathOfs;				// vertex sequence is starting at this path offset (0..CurrentPathLength)
				 UINT32 Off5;					// vertex sequence to include in path starts from this sequence 5' offset (0..SeqLen-1)
				 bool bRevCpl;					// if true then reverse complement the vertex sequence before applying Off5 and Off3
	} tsPathTraceBack;

// disconnected graph components
typedef struct TAG_sComponent {
	tComponentID ComponentID;			// identifies this component
	tVertID VertexID;					// one of the vertices which is in this component
	UINT32 NumVertices;					// there are this many vertices in this component
	tVertID PathStartVertexID;			// highest scoring path starts from this vertex
	UINT64 PathScore;					// highest scoring path has this score
	UINT32 PathLength;					// highest scoring path is this length
	UINT32 NumPathEdges;				// number of edges in path
	UINT32 NumTraceBacks;				// number of tracebacks in path
	UINT32 StartTraceBackID;			// path starts with this traceback
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

	bool m_bSenseOnlyOvlps;				// if false then both sense/sense and sense/antisense overlaps will be accepted and processed, otherwise sense/sense only overlaps accepted and processed
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

	UINT32 m_UsedTraceBacks;			// currently using this many tracebacks
	UINT32 m_AllocdTraceBacks;			// allocd to hold this many tracebacks
	tsPathTraceBack *m_pPathTraceBacks; // to hold all path tracebacks


	UINT32 m_NumDiscRemaps;				// number of disconnected graph identifiers requiring remaps
	tsRemapComponentID m_DiscRemaps[cMaxDiscRemaps];	// to hold disconnected graph identifiers requiring remaps
	UINT32 m_NumReducts;			// number of graph node reductions
	bool m_bReduceEdges;			// if true then attempt to reduce extraneous graph edges when current graph is full and more node memory needs to be allocated  

	bool m_bAcceptOrphanSeqs;		// if true then report also report sequences which have no overlap with any other sequence

	int m_NumThreads;				// use at most this number threads for graph processing
	bool m_bMutexesCreated;			// set true if mutexes and rwlocks created/initialised

#ifdef WIN32
	alignas(4)	volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access - replaces AcquireSerialise() as much more efficient
	alignas(4)	volatile unsigned int m_CASLock; // used with synchronous compare and swap (CAS) for serialising access -  - replaces AcquireLock(True) as much more efficient
#else
	__attribute__((aligned(4))) volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access - replaces AcquireSerialise() as much more efficient
	__attribute__((aligned(4))) volatile unsigned int m_CASLock; // used with synchronous compare and swap (CAS) for serialising access -  - replaces AcquireLock(True) as much more efficient
#endif
	void AcquireCASSerialise(void);
	void ReleaseCASSerialise(void);
	void AcquireCASLock(void);
	void ReleaseCASLock(void);


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
		LocateFirstFwdEdgeID(tVertID FromVertexID);	// find first matching 
	
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

	UINT32							// number of vertices with both inbound and outbound edges
		VertexConnections(void);		// identify and mark vertices which have multiple inbound edges
	UINT32 GenSeqFragment(tsGraphVertex *pVertex);		// initial seed vertex
	UINT32  TransitIdentDiscGraph(tVertID VertexID,tComponentID ComponentID);
	UINT32  ClearEdgeTravFwdRevs(void);
	UINT32	ClearDiscCompIDs(void);

public:
	CAssembGraph(void);
	~CAssembGraph(void);

	void Reset(void);								// reset and free any allocated resources
	teBSFrsltCodes		// initialise with ScaffScoreThres and  maxThreads
		Init(bool bSenseOnlyOvlps = false,				// if false then both sense/sense and sense/antisense overlaps will be accepted and processed, otherwise sense/sense only overlaps accepted and processed
			 int ScaffScoreThres = cDfltMin1kScore,		// accepted edges must be of at least this overlap score
			 bool bAcceptOrphanSeqs = false,			// also report sequences which are not overlapped or overlapping any other sequence
						int MaxThreads = 8);			// set number of threads

	UINT32								// returned vertex identifier
		AddVertex(UINT32 SeqLen,		// sequence length
				  tSeqID SeqID);		// allocates and initialises a new graph vertex which references this sequence

	UINT32								// 0 if errors else number of vertices finalised
		FinaliseVertices(void);			// added vertices are sorted ready for graph edges to be added

	UINT32								// number of edges removed	 
		ReduceEdges(void);				// reduce graph by detecting and removing extraneous edges


	UINT32									// returns total number of edges, including this edge if accepted, thus far accepted 
		AddEdge(tSeqID FromSeqID,			// identifies the 'From' or overlapping sequence
				tSeqID ToSeqID,				// identifies the 'To' overlapped sequence
				UINT32 FromSeqLen,			// 'From' sequence length is this length
				UINT32 ToSeqLen,			// 'To' sequence length is this length
				UINT32 Score,				// score associated with this overlap, higher scores represent higher confidence in the overlap
				UINT32 ScoreAlignLen,		// scored for this alignment length: (1 + ProbeAlignLen + TargAlignLen) / 2
				UINT32 FromSeq5Ofs,			// overlap of FromSeqID onto ToSeqID starts at this base relative to FromSeqID 5' start
				UINT32 FromSeq3Ofs,			// overlap of FromSeqID onto ToSeqID ends at this base relative to FromSeqID 5' start
				UINT32 ToSeq5Ofs,			// overlap onto ToSeqID from FromSeqID starts at this base relative to ToSeqID 5' start
				UINT32 ToSeq3Ofs,			// overlap onto ToSeqID from FromSeqID ends at this base relative to ToSeqID 5' start
				eOverlapClass OverlapClass,	// classification of overlap from FromSeqID onto ToSeqID, note that classification must be eOLCOverlapping
				bool bSenseOvlpAnti);	// false: 'From' sense overlaps 'To' sense, true: 'From' sense overlaps 'To' antisense 

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

	INT32											// returned From sequence extension; -1 if no sequence extension
	OverlapAcceptable(tsGraphOutEdge *pEdge,		// overlap edge
				UINT8 FromOvlpClass = 0,	// From vertex overlap classification; bit 0 set if From vertex evaluated as antisense in current path
				UINT8 *pToOvlpClass = NULL);		// returned To vertex overlap classification; bit 0 set if To vertex is evaluated as antisense in current path

	UINT64
		ScorePaths(tVertID VertexID);			// score all paths starting with outgoing edges from this vertex

	UINT64										// highest scoring of any path from pEdge
		ScorePath(UINT32 Depth,					// current recursive depth - used to detect circular paths
				  tsGraphOutEdge *pEdge,		// score paths starting with this edge
				UINT8 FromOvlpClass = 0,	// From vertex overlap classification; bit 0 set if From vertex evaluated as antisense in current path
				UINT8 *pToOvlpClass = NULL);	// returned To vertex overlap classification; bit 0 set if To vertex is evaluated as antisense in current path



	int										 // eBSFSuccess or otherwise
		FindHighestScoringPaths(void);		 // score all possible paths and record highest scoring path for each component

	int											 // eBSFSuccess or otherwise
			ReportVerticesEdgesGEXF(char *pszOutFile,		 // report as GEXF format all vertices and their edges for each component to this file
								CSeqStore *pSeqStore);	// holds sequences used to assemble contig

	int												 // eBSFSuccess or otherwise
			ReportVerticesEdgesGraphML(char *pszOutFile,		 // report as GraphML on all vertices and their edges for each component to this file
										   CSeqStore *pSeqStore);	// holds sequences used to assemble contig

	int											// returned CurrentPathLength
		AddTraceBackPath(bool bFirst,			// true if 1st vertex in path (starting a new path)
				 tComponentID ComponentID,      // path is for this component
				 tVertID VertexID,				// path includes this vertex sequence
				 UINT32 SeqLen,					// vertex sequence length
				 UINT32 PathOfs,				// vertex sequence is starting at this path offset (0..CurrentPathLength)
				 UINT32 Off5,					// vertex sequence to include in path starts from this sequence 5' offset (0..SeqLen-1)
				 bool bRevCpl);					// if true then reverse complement the vertex sequence before applying Off5

	int										// eBSFSuccess or otherwise
		GenTraceBackPath(tsComponent *pComponent); // generate traceback path for this component

	UINT32	IdentifyDiscComponents(void);

	int WriteContigSeqs(char *pszOutFile,CSeqStore *pSeqStore);  // write assmbled PacBio contig sequences to this output file
};



