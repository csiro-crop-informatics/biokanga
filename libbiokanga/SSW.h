#pragma once
// Striped Smith-Waterman
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

#include "./commdefs.h"

const UINT32 cSSWMinProbeOrTargLen = 5;			// require probe and target lengths to be at least this number of bp
const UINT32 cSSWMaxProbeOrTargLen = 100000000;  // require either probe or target length to be no longer than this limit
const UINT32 cSSWMaxPenaltyGap = 1000000;		 // can specify gaps of upto this length before gap extension penalties are applied

const int cSSWDfltMatchScore = 1;		// score for matching bases
const int cSSWDfltMismatchPenalty = -1;	// mismatch penalty
const int cSSWDfltGapOpenPenalty = -3;	// gap opening penalty
const int cSSWDfltGapExtnPenalty = -1;	// gap extension penalty
const int cSSWDfltDlyGapExtn = 2;       // delay applying gap extension penalty until 2nd base in gap (must be range 1..cSSWMaxPenaltyGap)
const int cSSWDfltProgPenaliseGapExtn = 0; // if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles

const int cSSWMinAnchorLen = 3;			// can specify exactly matching anchors down to this length 
const int cSSWDfltAnchorLen = 4;		// default exactly matching anchors of at least this length to be identified 
const int cSSWMaxAnchorLen = 100;      // can specify exactly matching anchors of at most this length  

const int cMaxInitiatePathOfs = 500;    // default is to require SW paths to have started within this many bp on either the probe or target - effectively anchoring the SW
const int cMinNumExactMatches = 100;    // default is only consider a path as being a peak path if that path contains at least this many exactly matching bases

const int cMaxTopNPeakMatches = 100;     // can process for at most this many peak matches in any probe vs target SW alignment

const int cMaxMAFBlockErrCorLen = 250000;		// allowing for error corrected read sequences of up to this length
const int cMaxMAFBlockLen = (cMaxMAFBlockErrCorLen * 100);	// allowing for multialignment format block buffering of up to this length

#pragma pack(1)
// each cell in stripe or column
typedef struct TAG_sSSWCell {
	UINT32 StartPOfs;			// path started at this probe sequence base offset + 1, 0 if no path
	UINT32 StartTOfs;			// path started at this target sequence base offset + 1, 0 if no path
	UINT32 EndPOfs;				// path ended at this probe sequence base offset + 1, 0 if no path
	UINT32 EndTOfs;				// path ended at this target sequence base offset + 1, 0 if no path

	UINT32 PFirstAnchorStartOfs; // first anchor of at least m_AnchorLen exactly matching bases started at this probe sequence base offset + 1, 0 if no anchor
	UINT32 TFirstAnchorStartOfs; // first anchor of at least m_AnchorLen exactly matching bases started at this target sequence base offset + 1, 0 if no anchor
	UINT32 PLastAnchorEndOfs;	 // last anchor of at least m_AnchorLen exactly matching bases ended at this probe sequence base offset + 1, 0 if no anchor
	UINT32 TLastAnchorEndOfs;	 // last anchor of at least m_AnchorLen exactly matching bases ended at this probe sequence base offset + 1, 0 if no anchor

	UINT32 NumMatches;			// total number of bases accepted along the path as matching, either exactly or as mismatches
	UINT32 NumExacts;			// total number of bases accepted along the path as exactly matching
	UINT32 NumGapsIns;			// total number of insertion gaps opened along the path
	UINT32 NumGapsDel;		    // total number of deletion gaps along the path
	UINT32 NumBasesIns;		    // total number of bases, over all gaps, inserted into probe
	UINT32 NumBasesDel;		    // total number of bases, over all gaps, deleted from probe

	INT32 PeakScore;			// peak score along path
	INT32 CurScore;				// current score
	UINT32 LeftInDelLen;		// current left (insert relative to target) gap length
	UINT32 DownInDelLen;		// current down (deletion relative to target) gap length
	UINT32 CurExactLen;			// current exactly matching sequence length
	
} tsSSWCell;


typedef struct TAG_sSSWTraceback {
	UINT32 IdxP;							// current index into m_Probe[]
	UINT32 IdxT;							// current index into m_Targ[]
} tsSSWTraceback;
const UINT32 cTrBkIdxMsk = 0x3fffffff;		// IdxT is in lower 30bits of tsSSWTraceback.IdxT, IdxP is in lower 30bits of tsSSWTraceback.IdxP  
const UINT32 cTrBkFlgsMsk = 0xc0000000;		// traceback direction is in bits 30..31 of IdxP, supplementary flags for traceback reduction and substution when matching are in bits 30..31 of IdxT
const UINT32 cTrBkFlgStart = 0x00000000;	// traceback direction: 5' start of alignment, no further traceback
const UINT32 cTrBkFlgMatch = 0x40000000;	// traceback direction: matching base, could be a sub so check on cTrBkFlgSub in IdxT, trace back to IdxT-1, IdxP-1
const UINT32 cTrBkFlgIns = 0x80000000;	    // traceback direction: base inserted into probe relative to target, trace back to IdxT, IdxP-1
const UINT32 cTrBkFlgDel = 0xc0000000;	    // traceback direction: base deleted from probe relative to target, trace back to IdxT-1, IdxP
const UINT32 cTrBkFlgRetain = 0x80000000;	  // if set then retain this traceback when reducing tracebacks
const UINT32 cTrBkFlgSub = 0x40000000;	     // if set then substution was required to match

typedef UINT8 tMAOp;				// alignment operations can be one of the following
const UINT8 cMAMatch = 0x00;		// base match between probe and target; note that may not be an exact match
const UINT8 cMAInsert = 0x01;		// base inserted into probe relative to target - or could be base deleted from target relative to probe
const UINT8 cMADelete = 0x02;		// base deleted from probe relative to target- or could be base inserted into target relative to probe

const UINT8 cMACompleted = 0x03;	// alignment completed


const int cMaxMultiAlignSeqs = 200;			// can process at most this number of sequences in a multialignment
const int cMaxParsimoniousAlignLen = 20;	// only processing for parsimonious multiple alignments for sequences which are no longer than this length
const UINT32 cMaxAcceptExpCombs = 10000000; // only accepting for parsimonious multiple alignments if the expected number of combinations of InDels over all unique sequences is no more than this threshold

typedef struct TAG_sPermInDels {
	int SeqID;					// identifies this sequence instance within the multialignment
	int SeqLen;					// sequence length including InDels
	int NxtSeqID;				// SeqID of next sequence with identical sequence, 0 if this sequence is unique
	UINT8 flgCharacterised:1;   // 0 until sequence has been characterised, if 1 then flgNoPermutate, flgAllInDels, flgNoInDels, flgRepInstance, and NumCopies have been set appropriately
	UINT8 flgNoPermutate:1;		// if 1 then do not permutate this sequence
	UINT8 flgUndef:1;			// sequence contains at least one undefined base, this sequence will not be further processed
	UINT8 flgAllInDels:1;		// sequence contains all InDels, no bases
	UINT8 flgNoInDels:1;		// sequence contains no InDels, is all bases
	UINT8 flgRepInstance:1;     // this is a representative instance, NumCopies only valid if a representative instance
	UINT32 NumCopies;           // number of exact copies of this sequence, includes this instance, only valid if a representative instance
	UINT64 SeqInDelMsk;			// mask for detecting InDel position wrap around required
	int NumInDels;				// number of InDels in sequence
	UINT64 InitialInDelPsns;	// bitmap of initial InDel positions in sequence
	UINT64 CurInDelPsns;		// bitmap of current InDel positions in permutated sequence
	etSeqBase InPermSeq[cMaxParsimoniousAlignLen];		// sequence instance curently being permutated
	etSeqBase ParsePermSeq[cMaxParsimoniousAlignLen];	// most parsimonious permutated instance
} tsPermInDels;

// multiple alignment columns
typedef struct TAG_sMAlignCol {
	UINT32 ColIdx;			// index + 1 in m_pMACols[] at which this alignment column starts, 1 if first
    UINT32 NxtColIdx;		// index + 1 in m_pMACols[] at which next alignment column starts, 0 if this is the last alignment column
    UINT32 PrevColIdx;		// index + 1 in m_pMACols[] at which previous alignment column starts, 0 if this is the first alignment column
	UINT32 Ofs;				// alignment col is for this probe relative offset (1..n)
	UINT16  Extn;			// if > 0 then relative offset to RefLoci of an inserted column 
	UINT16  Depth;			// currently with this coverage 
	UINT8 ConsConf;			// nearly arbitary confidence in the ConsBase; scale is from 0 for no confidence with increasing values representing more confidence
	UINT8  ConsBase;        // consensus base
	UINT8  Bases[1];		// probe base in  Bases[0] with target bases following 	
} tsMAlignCol;

typedef struct TAG_sConfBase {   // associating consensus base with corresponding concensus confidence
	UINT8 Base;		// consensus base 
	UINT8 Conf;		// confidence in base
} tsConfBase;

#pragma pack()

class CSSW
{
	bool m_bStartedMultiAlignments;     // set true after 1st call to StartMultiAlignments() and reset 
	int m_ErrCorSeqID;					// number of error corrected sequences reported by this class instance
	etSeqBase *m_pProbe;				// alloc'd probe sequence memory
	UINT32 m_ProbeAllocd;				// actual m_pProbe alloc'd size (in teSWSeqBases's)
	UINT32 m_ProbeLen;					// current probe length

	etSeqBase *m_pTarg;					// alloc'd target sequence memory
	UINT32 m_TargAllocd;				// actual m_pTarg alloc'd size (in teSWSeqBases's)
	UINT32 m_TargLen;					// current targ length

	int m_MatchScore;					// score for matching bases (0..100)
	int m_MismatchPenalty;				// mismatch penalty (-100..0)
	int m_GapOpenPenalty;				// gap opening penalty (-100..0)
	int m_GapExtnPenalty;				// gap extension penalty (-100..0)
	int m_DlyGapExtn;					// delayed gap penalties, only apply gap extension penalty if gap at least this length
    int m_ProgPenaliseGapExtn;			// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles

	int m_MaxInitiatePathOfs;			// if non-zero then only allow new paths to start if within that offset (0 to disable) on either probe or target - effectively an anchored SW

	UINT32 m_AnchorLen;				// identified anchors between aligned probe and target must be at least this length

	UINT32 m_UsedCells;				// number of currently allocated cells used
	UINT32 m_AllocdCells;			// number of currently allocated cells
	size_t m_AllocdCellSize;        // total current allocation size for m_pAllocdCells 
	tsSSWCell *m_pAllocdCells;		// allocated to hold cells	

	UINT32 m_UsedTracebacks;		// number of currently allocated Tracebacks used
	UINT32 m_AllocdTracebacks;		// number of currently allocated Tracebacks
	size_t m_AllocdTracebacksSize;  // total current allocation size for m_pTracebacks 
	tsSSWTraceback *m_pAllocdTracebacks;	// allocated to hold Tracebacks	

	UINT32 m_MAAlignOps;			// number of currently allocated alignment operators used
	size_t m_AllocdMAAlignOpsSize;		// total current allocation size for alignment operators  
	tMAOp *m_pMAAlignOps;			// remapped from tracebacks the alignment operators used when merging multiple alignments into a consensus sequence

	int m_MinNumExactMatches;		// peak matches must contain at least this many exact matches to qualify as a peak match
	tsSSWCell m_PeakMatchesCell;	// cell identified as path containing highest number of matches 
	tsSSWCell m_PeakScoreCell;		// cell identified as usual SW peak scoring as conventional in most SW scoring schemes

	int m_MaxTopNPeakMatches;		// can identify at most this many top peak matches for current probe vs target
	int m_NumTopNPeakMatches;		// currently identified this many top peak matches for current probe vs target
	tsSSWCell m_TopPeakMatches[cMaxTopNPeakMatches];  // used to hold the top peak matches

	tsPermInDels m_PermIndels[cMaxMultiAlignSeqs]; // used when determining the maximally scoring multialigned deletions relative to the probe sequence 

	UINT32 m_MAProbeSeqLen;			// probe sequence length used when initialising multialignment columns
	UINT32 m_MACols;					// actual number of multialignment columns used
	UINT32 m_MADepth;					// each column supports multialignments of this maximal depth
	UINT32 m_MACoverage;               // actual coverage depth over all columns
	UINT32 m_MAColSize;				// each column is this size

	UINT32 m_AllocMACols;				// this number of columns have been allocated
	size_t m_AllocMAColsSize;		// allocated memory size in bytes to hold multialignment cols
	UINT8 *m_pMACols;				// pts to allocated memory for multialignment cols 


	int m_AllocParsimoniousSize;	// total current allocation size for holding bases whilst processing for most parsimonious alignments within probe InDel gaps 
	UINT8 *m_pParsimoniousBuff;		// allocated to hold sequences when determining most most parsimonious alignments

	bool m_bIsGZ;					// true if processing a gz compressed multialignment file
	gzFile m_gzFile;				// opened for reading (if compressed) a multialignment file as may have been generated by MAlignCols2MFA()
	int m_hMAFFile;					// opened for reading (if not compressed) a multialignment file as may have been generated by MAlignCols2MFA()
	int m_hConsSeqFile;				// opened for writing the consensus sequences to
	INT64 m_MAFFileOfs;				// will contain file offset corresponding to last read from file
	int m_MAFAlignBuffIdx;			// offset in m_pszMAFAlignBlock at which next char is being parsed
	int m_MAFAlignBuffered;			// number of chars currently buffered in m_pszMAFAlignBuff
	INT64 m_MAFFileLineNum;			// line number in MAF file currently being parsed
	int m_AllocMAFAlignBuffSize;	// m_pszMAFAlignBlock allocated to hold this many chars
	char *m_pszMAFAlignBuff;		// allocated to buffer the MAF alignment blocks whilst parsing
	tsConfBase *m_pConsConfSeq;		// allocated to hold the parsed consensus bases and consensus confidence scores
	char *m_pszConsensusBuff;		// allocated to buffer the generated consensus sequence

	int    // parse out the next multialignment block consensus bases and consensus confidence scores into m_pConsConfSeq
			ParseConsConfSeq(bool bCpltdReadMAF,	// true if m_pConsConfSeq contains all remaining multialignment blocks loaded from file 
						 int MinErrCorrectLen,		// error corrected sequences must be at least this minimum length	
						 UINT32 *pProbeID);			// returned probe identifier

	int    // total number of returned chars in pszBuffer for the textual representation of error corrected consensus sequence (could be multiple consensus sequences)
		ConsConfSeq2Seqs(UINT32 ProbeID,	// identifies sequence which was used as the probe when determining the multialignments
				  int MinConf,				// sequence bases averaged over a 100bp window must be of at least this confidence (0..9) with the initial and final bases having at least this confidence
				  int MinLen);				// and sequence lengths must be of at least this length 

	UINT64	MSBBitMsk(UINT64 BitsSet);			// get bit mask of most significant bit set in BitsSet
	UINT64	NSBBitMsk(int Nth,UINT64 BitsSet);	// get bit mask of Nth (1 if MSB) significant bit set in BitsSet
	int	NumBitsSet(UINT64 BitsSet);				// get total number of bits set in BitsSet

	tsMAlignCol *								// inserted column or NULL if errors inserting
		InsertCol(UINT32 PrevCol);				// allocate and insert new column after PrevCol

	int								// higest score for path anchored at the pProbe[0] and pTarg[0]  
		ChkStartPath(int Scorethres, // looking for a path of at least this SW score
			  int SeqLen,			// both pProbe and pTarg are at least this length, max of 100 will be processed
			  etSeqBase *pProbe,	// SW this sequence against
			  etSeqBase *pTarg,		// this sequence
			  int SeedScore=10,		// seed the path with this initial score
			  int MatchScore=cSSWDfltMatchScore,			// exact match score
			  int MismatchPenalty=cSSWDfltMismatchPenalty,	// mismatch penalty
			  int GapOpenPenalty=cSSWDfltGapOpenPenalty,		// gap open penalty
			  int GapExtnPenalty=cSSWDfltGapExtnPenalty);// gap extension penalty

public:
	CSSW();
	~CSSW();
	void Reset(void);			    // reset state back to that immediately following instantiation

	bool SetScores(int MatchScore= cSSWDfltMatchScore,			// score for match
				int MismatchPenalty  = cSSWDfltMismatchPenalty,	// penalty for mismatch
				int GapOpenPenalty  = cSSWDfltGapOpenPenalty,	// penalty for opening a gap
				int GapExtnPenalty  = cSSWDfltGapOpenPenalty,	// penalty if extending already opened gap
				int DlyGapExtn = cSSWDfltDlyGapExtn,			// delayed gap penalties, only apply gap extension penalty if gap at least this length
				int ProgPenaliseGapExtn = cSSWDfltProgPenaliseGapExtn,	// if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles
				int AnchorLen = cSSWDfltAnchorLen);				// identified first and last anchors in alignment to be of at least this length

	bool SetMaxInitiatePathOfs(int MaxInitiatePathOfs = cMaxInitiatePathOfs);	// require SW paths to have started within this many bp (0 to disable) on either the probe or target - effectively anchoring the SW 

	bool SetMinNumExactMatches(int MinNumExactMatches = cMinNumExactMatches);		// require at least this many exactly matching in path to further process that path
	bool SetTopNPeakMatches(int MaxTopNPeakMatches = cMaxTopNPeakMatches/2);		// can process for at most this many peak matches in any probe vs target SW alignment

	int GetTopNPeakMatches(tsSSWCell **pPeakMatches = NULL);		// returns number of peak matches in any probe vs target SW alignment


	bool SetProbe( UINT32 Len,etSeqBase *pSeq);					// set probe sequence to use in subsequent alignments
	bool SetTarg( UINT32 Len,etSeqBase *pSeq);					// set target sequence to use in subsequent alignments

	bool PreAllocMaxTargLen( UINT32 MaxTargLen, bool bNoTracebacks = false);				// preallocate to process targets of this maximal length, with option to note allocate for tracebacks


	tsSSWTraceback *InitiateTraceback(UINT32 IdxP, UINT32 IdxT);
	tsSSWTraceback *								// next traceback along path; NULL if at 5' start of path
		NxtTraceback(tsSSWTraceback *pCur);			// use to iterate over tracebacks

	UINT32											// number of tracebacks marked, can be less than actual path length if some tracebacks already marked
		 MarkTracebackPath(UINT32 MarkFlag,			// mark the traceback path which 3' ends at IdxP and IdxT with this marker flag(s) into the tsSSWTraceback.IdxT
				UINT32 IdxP, UINT32 IdxT);

	UINT32											// after reduction there are this many tracebacks retained
		ReduceTracebacks(UINT32 RetainFlag,			// reduce tracebacks by removing tracebacks which have NOT been marked with this flag in IdxT
					 UINT32 ResetFlags);		    // and reset these flags in the retained traceback in IdxT

	int ValidateTracebacks(tsSSWCell *pCell);		// valdate that the path can be traced back for this cell

	tsSSWCell *										// smith-waterman style local alignment, returns highest accumulated exact matches scoring cell
				Align(tsSSWCell *pPeakScoreCell = NULL,	// optionally also return conventional peak scoring cell
						bool bNoTracebacks = false);	// if false then not using tracebacks

	double											// parsimony of multiple alignment, 0 (min) to 1.0 (max) 
		ParsimoniousMultialign( int Depth,			// parsimony for this number of bases in each column
								int NumCols,		// alignment parsimony over this many columns
								tsMAlignCol *pInDelStartCol); // starting at this column

	double										// parsimony of returned (pSequences) multiple alignment, 0 (min) to 1.0 (max) 
		ParsimoniousBasesMultialign(int NumSeqs,		// number of sequences
		   int SeqLen,							// each sequence is this length 
		   etSeqBase *pSequences);				// input: ptr to each sequence concatenated together, sequences may contain eBaseA, eBaseC,eBaseG,eBaseT,eBaseN and eBaseInDel
												// output:pSequences updated with most parsimonious   				

	bool									// false: sequence has been completely permutated, true: new permutation 
		PermuteInDels(tsPermInDels *pPerm);

	UINT64										// returned bit mask of next combination 
		Nxt_nCr(int nPossibilities,				// from n possibilities
				int rOutComes,					// choose r outcomes 
			   UINT64 PrevCombination);			// bit mask of previous combination generated


	int												// number of alignment ops generated
		TracebacksToAlignOps(UINT32 ProbeStartOfs,	// alignment starts at this probe sequence offset (1..n)
					UINT32 ProbeEndOfs,				// alignment ends at this probe sequence offset
					UINT32 TargStartOfs,			// alignment starts at this target sequence offset (1..n)
					UINT32 TargEndOfs);				// alignment ends at this target sequence offset

	int
		StartMultiAlignments(int SeqLen,			// probe sequence is this length
					etSeqBase *pProbeSeq,			// probe sequence 
					int Alignments);				// number of pairwise alignments to allocate for

	int
		AddMultiAlignment(UINT32 ProbeStartOfs,		// alignment starts at this probe sequence offset (1..n)
					  UINT32 ProbeEndOfs,			// alignment ends at this probe sequence offset inclusive
					  UINT32 TargStartOfs,			// alignment starts at this target sequence offset (1..n)
					  UINT32 TargEndOfs,			// alignment ends at this target sequence offset inclusive
					  etSeqBase *pTargSeq);			// alignment target sequence

	int
	GenMultialignConcensus(void);

	int      // total number of returned chars in pszBuffer for the textual representation of the error corrected consensus sequence (could be multiple consensus sequences)
		MAlignCols2fasta(UINT32 ProbeID,	// identifies sequence which was used as the probe when determining the multialignments
				  int MinConf,	// sequence bases averaged over 100bp must be of at least this confidence (0..9)
				  int MinLen,			// and sequence lengths must be of at least this length 
				  UINT32 BuffSize,		// buffer allocated to hold at most this many chars
				  char *pszBuffer);		// output error corrected sequences to this buffer

	int
		GenConsensusFromMAF(int MinErrCorrectLen,	// error corrected sequences must be at least this minimum length
				 int MinConcScore,			// error corrected sequences trimmed until mean 100bp concensus score is at least this threshold
				char *pszErrCorFile,		// name of file into which write error corrected sequences
				char *pszMultiAlignFile);	// name of file containing multiple alignments to process

	int      // total number of returned chars in pszBuffer for the textual representation of the multialignment 
		MAlignCols2MFA( UINT32 ProbeID,		// identifies sequence which was used as the probe when determining the multialignments
					    UINT32 BuffSize,	// buffer allocated to hold at most this many chars
					    char *pszBuffer);	// output multialignment textual representation to this buffer


};

