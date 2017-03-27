/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once


const int cMinScaffLen = 100;			// user can specify down to this length for minimum reported scaffolded sequences
const int cDfltMinScaffLen = 300;		// default length for minimum reported scaffolded sequences
const int cMaxScaffLen = 5000;			// user can specify upto to this length for minimum reported scaffolded sequences

const int cMaxSubs100bp = 5;			// user can specify up to this many substitutions per 100bp of overlap

const size_t cSeqEdges2Alloc = 25000000; // initially allocate to hold this many overlap edges

const int cMaxOverlapsProbeTarg   = 1000;   // limit to at most this many overlaps of a probe onto a single target sequence
										     // if probe overlap limits onto target exceeded then no overlaps for that probe is accepted 							

const int cSenseHistoryDepth = 5;			// maintain sense change history to this depth

const int cDfltScaffoldGapNs = 10;			 // gaps in scaffolded sequences are represented by this many 'N's
const int cAllocScaffoldBuffSize = 25000000; // allocation size for buffering scaffold sequences to be written to file

const size_t cMainThreadStackSize = (1024*1024*10); // main processing thread may be (in a future redevelopment, but not currently) be using
										// recursive functions when graph processing and thus may be requiring plenty of stack as recursion depth can be very deep 

#pragma pack(1)

typedef struct TAG_sScaffoldPars {  // to hold thread parameters for thread doing main processing with increased stack size (cThreadStackSize) because of potential graph recursion
#ifdef _WIN32
	HANDLE threadHandle;				// handle as returned by _beginthreadex()
	unsigned int threadID;				// identifier as set by _beginthreadex()
#else
	int threadRslt;						// result as returned by pthread_create ()
	pthread_t threadID;					// identifier as set by pthread_create ()
#endif
	int Rslt;							// returned result code
	int PMode;							// processing mode: 0 - output scaffold multifasta with edge report
	int NumThreads;						// number of worker threads to use
	bool bAffinity;						// thread to core affinity
	int Subs100bp;						// allow this many induced substitutions per 100bp overlapping sequence fragments (defaults to 1, range 0..5)
	int MaxEnd12Subs;					// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
	int MinPEReadlen;					// only accept PE read lengths of at least this many bp
	int MinPEInsertSize;				// PE sequences are of this minimum insert size
	int MaxPEInsertSize;				// PE sequences are of this maximum insert size
	int MinScaffoldedSeqLen;			// reported scaffolded sequences are at least this length
	int OrientatePE;					// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
	char *pszPE1File;					// input PE1 sequences file
	char *pszPE2File;					// input PE2 sequences file
	char *pszContigsFile;				// input SE contigs file
	char *pszScaffoldedFile;			// where to write scaffolded contigs
} tsScaffoldPars;

// sequence edge overlaps FromSeqID onto ToSeqID 
typedef struct TAG_sSeqEdge {
	tEdgeID SeqEdgeID;			// uniquely identifies this overlap edge
	int ProcPhase;				// overlap was discovered during this processing phase - 0 sense/sense, 1 antisense/sense 
	tSeqID ProbeSeqID;			// probe sequence identifier
	tVertID ProbeSeqVertID;		// probes vertice
	int ProbeLen;				// probe length
	tSeqID ProbeMateSeqID;		// as probe is always a PE, then it's corresponding mate end identfier
	tSeqID TargetSeqID;			// target sequence which has been overlapped, will always be a SE
	tVertID TargetSeqVertID;	// targets vertice
	int TargetLen;				// target sequence length
	UINT16 FlgIsPE2:1;			// if edge is from PE1, 1 if edge is from PE2 of the PE
	UINT16 FlgRevCpltd:1;		// probe sequence has been reverse complemented
	UINT16 FlgSlough:1;			// slough this edge when scoring
	INT32 zFlankLen;			// sequence 5' offset at which overlap starts
	INT32 zOverlapLen;			// FlgFlankLen sequence overlaps by this many bases
	INT32 NumSubs;				// number of subs required
} tsSeqEdge;


// sequence graph vertex, references a SE contig with edge linkages to other contigs
typedef struct TAG_sSeqVertex {
	tVertID SeqVertID;			// uniquely identifies this vertex
	UINT16 FlgSenseChkd:1;		// set if this vertex has been checked for relative sense against the reference vertex
	UINT16 FlgReqRevCpl:1;		// set true if this vertex is to be RevCpl to make it's sense same as reference vertex
	UINT16 FlgCircChkd:1;		// true if this vertex has been checked for circular references
	UINT16 FlgPathed:1;			// if already commited into a path then can't be scored into a new path

	tSeqID SeqID;				// sequence identifier
	UINT32 SeqLen;				// sequence length
	tSeqID PredSeqID;			// predecessor sequence identifier
	UINT16 PredNumPEs;			// number of supporting PEs linking relative to this sequence
	INT16 PredScore;			// linkage score - negative if predecessor is antisense to this SeqID, positive if sense, higher absolute scores reflect higher confidence
	tSeqID SuccSeqID;			// successor sequence identifier
	UINT16 SuccNumPEs;			// number of supporting PEs linking relative to this sequence
	INT16 SuccScore;			// linkage score - negative if successor is antisense to this SeqID, positive if sense, higher absolute scores reflect higher confidence
} tsSeqVertex;

typedef struct TAG_sThreadProcScaffoldsPars {
	int ThreadIdx;					// index of this thread (1..m_NumThreads)
	void *pThis;					// will be initialised to pt to CKangadna instance

#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// returned result code

	int CurPhase;					// 0 if targets indexed on sense; 1 if targets have been RevCpl and then indexed
	bool bPEOnly;					// if false then identify both SE and PE (both ends) sequences which are fully contained, if true then identify contained PEs only

	int MinReqPESepDist;			// PE start sites expected to be separated by at least this many bases (minimum insert size - PE2 len)
	int MaxReqPESepDist;			// PE start sites expected to be separated by no more than this many bases (maximum insert size - PE2 len)

	int MaxSubs1K;					// allow at most this many substitutions per 1K overlapping bases
	int MaxEnd12Subs;				// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 

	bool bPESeqs;					// sequences include paired end reads (could also include SE/contigs); if reset then only SE/contig sequences
	
	bool bPESeqsRevCpl;				// true if PE1 and PE2 sequences are currently revcpl'd
	
	UINT32 AllocProbeSeqLen;		// memory allocated to pProbeSeq
	UINT32 MaxProbeSeqWrds;			// pProbeSeq allocated to hold at most this many tSeqWrds
	void *pProbeSeqWrds;			// to hold copy of probe sequence 
	void *pRevCplProbeSeqWrds;		// to hold revcpl copy of probe sequence

	UINT32 AllocMateSeqLen;			// memory allocated to pMateSeq
	UINT32 MaxMateSeqWrds;			// pMateSeq allocated to hold at most this many tSeqWrds
	void *pMateSeqWrds;				// to hold copy of mate end sequence 
	void *pRevCplMateSeqWrds;		// to hold revcpl copy of mate end sequence 
	
	UINT32 NumProcessed;			// number processed
	UINT64 NumOverlapped;			// number of sequences determined as being overlapped
	UINT32 NumPE1Overlapping;		// number of PE1 sequences which overlapped other sequences
	UINT32 NumPE2Overlapping;		// number of PE2 sequences which overlapped other sequences
} tsThreadProcScaffoldsPars;

#pragma pack()

class CScaffolder : public CdeNovoAssemb
{
	int m_NumThreads;						// use at most this number of threads
	CMTqsort m_MTqsort;						// multithreaded sorting

	tSeqID m_LowestEdgeFromSeqID;				// lowest numbered FromSeqID in m_pSeqEdges 
	tSeqID m_HighestEdgeFromSeqID;				// highest numbered FromSeqID in m_pSeqEdges
	tSeqID m_LowestEdgeToSeqID;					// lowest numbered ToSeqID in m_pSeqEdges 
	tSeqID m_HighestEdgeToSeqID;				// highest numbered ToSeqID in m_pSeqEdges
	UINT32 m_NumPEOverlapping;					// number of PE overlapping other sequences
	UINT32 m_NumSEOverlapping;					// number of SE overlapping other sequences
	UINT32 m_NumPEOverlapped;					// number of PE which are overlapped by other sequences
	UINT32 m_NumSEOverlapped;					// number of SE which are overlapped by other sequences

	UINT64 m_NumSeqEdges;					// number of edges between vertices currently in m_pSeqEdges
	size_t m_AllocMemSeqEdges;				// mem allocated for holding overlap edges
	tsSeqEdge *m_pSeqEdges;					// allocated to hold all overlap edges (normally indexed by FromSeqID ascending, FromLen descending)
	size_t m_AllocMemToSeqEdges;			// mem allocated for holding overlap edge ToSeqID ptrs
	tsSeqEdge **m_ppToSeqEdges;				// ptrs into m_pSeqEdges (m_pSeqEdges[] is ordered by FromSeqID) which are ordered by ToSeqID ascending, ToLen descending

	UINT32 m_NumSeqVertices;				// number of sequence vertices currently in m_pSeqVertices
	size_t m_AllocMemSeqVertices;			// mem allocated for holding sequence vertices
	tsSeqVertex *m_pSeqVertices;			// allocated to hold all sequence vertices

	int m_hScaffoldFasta;					// used for writing scaffolded sequence sets to multifasta
	char m_szScaffoldSetsFile[_MAX_PATH];	// scaffold sequence sets are written to this file
	int m_ScaffoldBuffLen;					// scaffold sequence buffer currently holds this many chars ready to be written to file
	int m_AllocScaffoldBuff;				// scaffold sequence buffer allocation size
	char *m_pszScaffoldBuff;				// allocated to hold buffer scaffold sequences which are to be written to file

	UINT32 m_NumScaffoldSets;				// number of scaffold sets processed
	UINT32 m_ScaffoldSetScore;				// score for current scaffold set
	UINT32 m_ScaffoldSetGaps;				// number of gaps in current scaffold set

	int m_MinScaffoldedSeqLen;				// reported scaffolded sequences are to be at least this length
	UINT32 m_AcceptedScaffolds;				// there were this number of scaffolds writen to file
	UINT32 m_UnderlenScaffolds;				// there were this number of underlength scaffolds not writen to file


	int
		ScaffoldAssemble(bool bSenseStrandOnly,			// sequences from sense strand specific
						bool bSingleEnded,				// treat all sequences as being single ended even if loaded as paired ends
						int OrientatePE,				// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
						int NumThreads,					// number of worker threads to use
						bool bAffinity,					// thread to core affinity
						int MinPEReadLen,				// PE reads must be at least this many bp long
						char *pszPE1File,		// input high confidence seed PE1 sequences file
						char *pszPE2File,		// input high confidence seed PE2 sequences file
						char *pszSeedContigsFile); // input high confidence seed SE contigs file

	int
		GenSeqEdges(int CurPhase,					// 0 if lookinging for sense PE overlaps onto sense targets; 1 if looking for antisense (RevCpl'd) PE overlaps onto sense targets
						int Subs100bp,				// allow this many induced substitutions per 100bp overlapping sequence fragments (defaults to 1, range 0..5)
						int MaxEnd12Subs,			// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
						int MinReqPESepDist,		// PE start sites expected to be separated by at least this many bases
						int MaxReqPESepDist);		// PE start sites expected be separated by no more than this many bases

	// MarkContainedSeqs
	// Mark sequences fully contained by other SE sequences
	int
		MarkContainedSeqs(int CurPhase,		// 0 if targets indexed on sense; 1 if targets have been RevCpl and then indexed
							bool bPEOnly,				// if true then also identify (cFlgContainRemove) SE sequences which are fully contained, if false then identify contained PE only 
							int SubsKbp,				// allow this many induced substitutions per 1Kbp overlapping sequence fragments (defaults to 0, range 1..10)
							int MaxEnd12Subs,			// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
							int MinReqPESepDist,		// PE start sites expected to be separated by at least this many bases
							int MaxReqPESepDist);		// PE start sites expected be separated by no more than this many bases

	teBSFrsltCodes AddOverlapEdge(tsSeqEdge *pEdge);  // pts to pre-initialised edge, SeqEdgeID is initialised within this function

	int				// number of edges (in) from other sequences onto SeqID, assumes that the SeqID sequences is sense orientated
			NumEdgesIn(tSeqID SeqID);					// edges into this sequence


	int				// number of edges (out) from SeqID, assumes that the SeqID sequences is sense orientated
			NumEdgesOut(tSeqID SeqID);					// edges out from this sequence


	int											//  0 if none, 1 if predecessor only, 2 if successor only, 3 if both predecessor and successor identified
			GetPredSuccSeq(tsSeqVertex *pSeqVertex);	// sequence vertex to update with predecessor and successor sequence identifiers and scores


	teBSFrsltCodes GenerateScaffoldGraph(void);			// generates scaffold graph (edges (overlaps) + vertices (SE or PE1/PE2s)

	int											// number of sequences marked as being fully contained
		RemoveContainedSeqs(int Subs1Kbp,				// allow this many induced substitutions per Kbp overlapping sequence fragments (defaults to 10, range 0..50)
					int MaxEnd12Subs,			// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
					int MinPEInsertSize,		// PE sequences are of this minimum insert size
					int MaxPEInsertSize,		// PE sequences are of this maximum insert size
					bool bPEOnly = false);		// if false then identify both SE and PE (both ends) sequences which are fully contained, if true then identify contained PEs only 

	int		IdentifyNonSenseSeqs(void);						// all vertices marked as requiring revcpl to make sequences sense consistent are marked as not to be processed for overlaps

	tEdgeID			// 0 if unable to locate any with matching FromSeqID, otherwise m_pSeqEdges[index-1] of lowest ordered sequence edge				 
		LocateFirstEdgeFromSeqID(tSeqID FromSeqID);	// find lowest ordered matching sequence edge with matching FromSeqID 

	tEdgeID			// 0 if unable to locate any with matching ToSeqID, otherwise m_ppToSeqEdges[index-1] of lowest ordered sequence edge				 
		LocateFirstEdgeToSeqID(tSeqID ToSeqID);	// find lowest ordered matching sequence edge with matching ToSeqID 


	tVertID							// returned unique Vertex identifier
			AddSeqVertex(tsSeqVertex *pSE);

	tsSeqVertex *											// NULL if unable to locate any with matching SeqID				 
		LocateVerticesSeqID(tSeqID SeqID);					// find vertex with this SeqID 


	int														// total number of reported scaffolds 
		ReportScaffoldSets(char *pszScaffoldFile);			// report scaffold set sequences to this multifasta file

	static int SortVerticesSeqID(const void *arg1, const void *arg2); // function for sorting vertices on their SeqID
	static int SortSeqEdgesFromSeqID(const void *arg1, const void *arg2);
	static int SortSeqEdgesFromToSeqID(const void *arg1, const void *arg2);
	static int SortSeqEdgesToSeqID(const void *arg1, const void *arg2);

public:
	CScaffolder(void);
	~CScaffolder(void);

	int ScaffolderReset(bool bSync = true);
	int ScaffolderInit(void);

	teBSFrsltCodes SetNumThreads(int maxThreads);	// max number of threads allowed

	teBSFrsltCodes
			GenScaffoldedContigs(int PMode,				// processing mode, currently unused
					int Subs100bp,						// allow this many induced substitutions per 100bp overlapping sequence fragments (defaults to 0, range 0..5)
					int MaxEnd12Subs,					// allow at the initial 12bp (2 x hexamer primer artefacts) of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs 
				    int MinPEReadlen,					// only accept PE read lengths of at least this many bp
					int MinPEInsertSize,				// PE sequences are of this minimum insert size
					int MaxPEInsertSize,				// PE sequences are of this maximum insert size
					int MinScaffoldedSeqLen,			// reported scaffolded sequences are at least this length
					int OrientatePE,					// PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
					char *pszPE1File,					// input PE1 sequences file
					char *pszPE2File,					// input PE2 sequences file
					char *pszContigsFile,				// input SE contigs file
					char *pszScaffoldedFile);			// where to write scaffolded contigs

	int ProcContained(tsThreadProcScaffoldsPars *pPars);	// threaded processing for identifying contained sequences

	int ProcScaffolds(tsThreadProcScaffoldsPars *pPars);	// threaded processing for scaffolding


	tSeqID			// 0 if no sequence identifier after FromSeqID, otherwise the next ordered edge FromSeqID 
		IterateNextFromSeqID(tSeqID FromSeqID);

	tSeqID			// 0 if no sequence identifier after ToSeqID, otherwise the next ordered edge ToSeqID 
		IterateNextToSeqID(tSeqID ToSeqID);

	tsSeqEdge *			// NULL if unable to locate any with matching FromSeqID				 
			IterateEdgeFromSeqID(tSeqID FromSeqID,	// iterate over sequence edges with matching FromSeqID
									tEdgeID *pEdgeID);	// set *pEdgeID to 0 to return 1st sequence edge, returns edge identifier to use on next iteration of IterateEdgeFromSeqID()

	tsSeqEdge *			// NULL if unable to locate any with matching ToSeqID				 
			IterateEdgeToSeqID(tSeqID ToSeqID,	// iterate over sequence edges with matching ToSeqID
									tEdgeID *pEdgeID);	// set *pEdgeID to 0 to return 1st sequence edge, returns edge identifier to use on next iteration of IterateEdgeToSeqID()

};

