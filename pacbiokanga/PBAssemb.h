#pragma once
/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "./pacbiocommon.h"

const int cMaxInFileSpecs = 100;			// user can specify upto this many input files
const int cMinSeqLen = 1000;				// minimum sequence length accepted for processing
const int cMaxDupEntries = 10;				// report 1st 10 duplicate entry names
const int cMaxAllocBuffChunk = 0x00ffffff;	// buffer for fasta sequences is realloc'd in this sized chunks

const int cDfltScaffScoreExact = 1;				// scaffolding uses a scoring system on overlaps which is independent of the scoring used when error correcting
const int cDfltScaffScoreMismatch = -1;			// expecting relatively few mismatches in error reduced scaffolding sequence overlaps
const int cDfltScaffScoreGapOpen  = -3;			// expecting relatively few gap openings in error reduced scaffolding sequence overlaps
const int cDfltScaffScoreGapExtn  = -1;			// expecting very few gap extensions in error reduced scaffolding sequence overlaps


typedef enum TAG_ePBPMode {
	ePBPMScaffold,										// scaffolding mode, uses previously generated overlap loci detail csv file and error corrected sequences
	ePBPMToGEFX,										// convert overlap loci detail into GEFX format file ready for import into a graph visualisation toolset, e.g. Gephi
	ePBPMToGraphML										// convert overlap loci detail into GraphML format file ready for import into a graph visualisation toolset, e.g. Cytoscape
	} etPBPMode;


#pragma pack(1)

// identified overlap between probe and target sequence
typedef struct TAG_sPBAOverlaps {
	UINT8 flgAntisense:1;           // probe sequence was reverse complemented
	UINT32 ProbeEntryID;            // probe sequence suffix array identifier
	UINT32 TargEntryID;				// overlap from probe was onto this target suffix array identifier
	UINT32 ProbeStartOfs;           // overlap starts at this probe offset
	UINT32 TargStartOfs;            // overlap starts at this target offset
	UINT32 ProbeOverlapLen;         // probe overlap is of this length
	UINT32 TargOverlapLen;			// target overlap is of this length
} sPBAOverlaps;

typedef struct TAG_sPBAScaffNode {
	UINT32 NodeID;					// uniquely identifies this node
	UINT32 VertexID;				// assembly graph vertex identifier
	UINT32 EntryID;					// suffix array entry identifier for indexed sequence
	UINT32 SeqLen;					// length in bp of this scaffolding node sequence
	UINT8 flgCurProc:1;				// sequence is currently being processed
	UINT8 flgContained:1;			// sequence is fully contained within another sequence
	UINT8 flgUnderlength:1;			// sequence is under length
} tsPBAScaffNode;

#pragma pack()


class CPBAssemb
{
	etPBPMode m_PMode;						// processing mode

	UINT32 m_NumOverlapProcessed;			// number of PacBio reads processed for overlapping other PacBio reads
	UINT32 m_ProvOverlapping;               // number of PacBio reads overlapping at least one other PacBio read
	UINT32 m_ProvOverlapped;				// number of PacBio reads provisionally overlapped, could be containing, another PacBio read
	UINT32 m_ProvContained;					// number of PacBio reads provisionally contained within another PacBio read
	UINT32 m_ProvArtefact;					// number of PacBio reads provisionally only partially, likely an alignment artefact, contained within another PacBio read
	UINT32 m_ProvSWchecked;					// number of times SW used to identify overlaps

	UINT32 m_OverlapFloat;					// allow up to this much float on overlaps to account for the PacBio error profile
	UINT32 m_MinScaffSeqLen;				// individual target scaffold sequences must be of at least this length (defaults to 5Kbp)
	UINT32 m_MinScaffOverlap;				// pairs of targeted scaffold sequences must overlap by at least this many bp to be considered for merging into a longer scaffold sequence (defaults to 5Kbp) 

	int m_ScaffScoreExact;					// scaffolding uses a scoring system on overlaps which is independent of the scoring used when error correcting
	int m_ScaffScoreMismatch;				// expecting relatively few mismatches in error reduced scaffolding sequence overlaps
	int m_ScaffScoreGapOpen;				// expecting relatively few gap openings in error reduced scaffolding sequence overlaps
	int m_ScaffScoreGapExtn;				// expecting very few gap extensions in error reduced scaffolding sequence overlaps
	int m_MinScaffScoreThres;				// accepted overlaps must be at least this minimum score per Kbp overlap

	UINT32 m_NumRejectedMinScaffOverlap;	// this number of overlaps rejected because overlap is less than m_MinScaffOverlap
	UINT32 m_NumRejectedScoreThres;			// this number of overlaps rejected as being below m_MinScaffScoreThres threshold
	UINT32 m_NumRejectContained;			// this number of overlaps rejected because the overlap was classified as contained
	UINT32 m_NumRejectAntisense;            // this number of overlaps rejected because the overlap was sense/antisense and only sense/sense overlaps are being accepted
	UINT32 m_NumRejectArtefact;				// this number of overlaps rejected because the overlap was classified as being an artefact
	UINT32 m_NumRejectedMinSeqLen;			// this number of overlaps rejected because either the probe or target sequence was under length
	UINT32 m_NumAcceptedOverlaps;			// this number of overlaps accepted for scaffolding

	bool m_bAcceptOrphanSeqs;				// if true then report also report sequences which have no overlap with any other sequence
	bool m_bSenseOnlyOvlps;					// if false then both sense/sense and sense/antisense overlaps will be accepted and processed, otherwise sense/sense only overlaps accepted and processed

	int m_NumErrCorrectedFiles;					// number of input error corrected file specs
	char m_szErrCorrectedFiles[cMaxInFileSpecs][_MAX_PATH];		// input error corrected files
	char m_szOutFile[_MAX_PATH];			// where to write merged scaffolded sequences


	int m_NumThreads;							// maximum number of worker threads to use

	UINT32 m_NumPBScaffNodes;					// m_pPBScaffNodes currently holds many scaffolding nodes
	UINT32 m_AllocdPBScaffNodes;				// m_pPBScaffNodes allocated to hold this many scaffolding nodes
	tsPBAScaffNode *m_pPBScaffNodes;				// allocated to hold scaffolding nodes
	UINT32 *m_pMapEntryID2NodeIDs;				// used to map from suffix array entry identifiers to the corresponding scaffolding node identifier

	CSeqStore *m_pSeqStore;						// sequence store

	CAssembGraph *m_pAssembGraph;				// used to assemble PacBio overlapping sequences into scaffolds

	void Init(void);							// initialise state to that immediately following construction
	void Reset(bool bSync);						// reset state, if bSync true then fsync before closing output file handles
	int LoadTargetSeqs(char *pszTargFile);		// load sequences in this file into in memory suffix array; file expected to contain preindexed sequences 

	int LoadTargetSeqs(int MinSeqLen,int NumTargFiles,char **pszTargFiles);		// parse, and index sequences in this file into in memory suffix array; file expected to contain either fasta or fastq sequences

	int ProcessBioseqFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
				 char *pszFile);				// file containing sequences

	int ProcessFastaFile(int MinSeqLen,			// only accept for indexing sequences of at least this length
				char *pszFile);					// file containing sequences

	UINT32										// returned tsPBScaffNode node identifier
		MapEntryID2NodeID(UINT32 EntryID);		// suffix array entry identifier

	CMTqsort m_mtqsort;				// muti-threaded qsort

static int SortLenDescending(const void *arg1, const void *arg2); // Sort scaffolding nodes by length descending

	bool m_bMutexesCreated;			// will be set true if synchronisation mutexes have been created
	int CreateMutexes(void);
	void DeleteMutexes(void);

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

public:
	CPBAssemb();
	~CPBAssemb();

	int
	Process(etPBPMode PMode,		// processing mode
		int MinScaffSeqLen,			// individual scaffold sequences must be of at least this length (defaults to 5Kbp)
		int MinScaffOverlap,		// pairs of targeted scaffold sequences must overlap by at least this many bp to be considered for merging into a longer scaffold sequence (defaults to 5Kbp) 
		int Min1kScore,             // minimum normalised 1Kbp overlap score
		bool bAcceptOrphanSeqs,		// also accepting sequences which are not overlapped or overlapping any other sequence
		bool bSenseOnlyOvlps,		// if false then both sense/sense and sense/antisense overlaps will be accepted and processed, otherwise sense/sense only overlaps accepted and processed
		char *pszMAFFile,			// pregenerated multialignment sequence overlap loci details
		int NumErrCorrectedFiles,	// number of error corrected sequence specs
		char *pszErrCorrectedFiles[],		// input error corrected sequence files
	    char *pszOutFile,			// where to write merged scaffolded sequences
		int NumThreads);				// maximum number of worker threads to use

	UINT32   //  returns number of overlaps loaded and accepted, if > cMaxValidID then cast to teBSFrsltCodes for actual error 
		LoadPacBioOvlps(char *pszPacBioOvlps,			// parse and load pregenerated PacBio sequence overlap loci CSV file 
						bool bValidateOnly = false,		// true if parse and validate only
						bool bSenseOnlyOvlps = false);				// if false then both sense/sense and sense/antisense overlaps will be accepted and processed, otherwise sense/sense only overlaps accepted and processed

};


