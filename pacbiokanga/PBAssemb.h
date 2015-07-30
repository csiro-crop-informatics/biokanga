#pragma once
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
	ePBPMScaffold										// scaffolding mode, uses previously generated overlap loci detail csv file and error corrected sequences
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
	UINT32 m_NumRejectedScoreThres;			// this number of overlaps rejected as being below m_MinScaffScoreThres threshold
	UINT32 m_NumRejectContained;			// this number of overlaps rejected because the overlap was classified as contained
	UINT32 m_NumRejectArtefact;				// this number of overlaps rejected because the overlap was classified as being an artefact
	UINT32 m_NumAcceptedOverlaps;			// this number of overlaps accepted for scaffolding


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

	volatile unsigned int m_CASSerialise; // used with synchronous compare and swap (CAS) for serialising access - replaces AcquireSerialise() as much more efficient
	void AcquireCASSerialise(void);
	void ReleaseCASSerialise(void);
	volatile unsigned int m_CASLock; // used with synchronous compare and swap (CAS) for serialising access -  - replaces AcquireLock(True) as much more efficient
	void AcquireCASLock(void);
	void ReleaseCASLock(void);

public:
	CPBAssemb();
	~CPBAssemb();

	int
	Process(etPBPMode PMode,		// processing mode
		int MinScaffSeqLen,			// individual scaffold sequences must be of at least this length (defaults to 5Kbp)
		int MinScaffOverlap,		// pairs of targeted scaffold sequences must overlap by at least this many bp to be considered for merging into a longer scaffold sequence (defaults to 5Kbp) 
		char *pszMAFFile,			// pregenerated multialignment sequence overlap loci details
		int NumErrCorrectedFiles,	// number of error corrected sequence specs
		char *pszErrCorrectedFiles[],		// input error corrected sequence files
	    char *pszOutFile,			// where to write merged scaffolded sequences
		int NumThreads);				// maximum number of worker threads to use

	UINT32   //  returns number of overlaps loaded and accepted, if > cMaxValidID then cast to teBSFrsltCodes for actual error 
		LoadPacBioOvlps(char *pszPacBioOvlps, bool bValidateOnly = false);	// load pregenerated PacBio sequence overlap loci CSV file
};

