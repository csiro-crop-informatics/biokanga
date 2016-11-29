#pragma once
#include "./commdefs.h"

// new release
const int cSFXVersion = 5;				// current file structure version
const int cSFXVersionBack = 3;			// can handle previous file structures back to this version

const int cSigWaitSecs = 5;				// background readahead thread wakes every cSigWaitSecs sec just in case a signalling event missed
const int cSigMaxWaitSecs = 600;		// master thread times out after cSigMaxWaitSecs
const int cSigTermWaitSecs = 2;			// background readahead thread is allowed at most cSigTermWaitSecs to self terminate if requested to terminate

const int cDfltMaxIter = 50000;			// default max iterations per subsegmented sequence when matching that subsegment
const int cMaxKmerLen = 17;			// limit on length of KMers which can be frequency counted when checking for over-occurances

const int cMaxNumIdentNodes = 1024000;	// allow at most this many TargSeqIDs to be hash linked per thread
const int cHashEntries = 0x03fff;		// TargSeqID start loci are hashed into this many entries

const int cMaxMMExploreInDel = 7;		// only if more than this number of mismatches then worth the effort for exploring microInDels
const int cMinInDelSeqLen = 7;			// only explore for InDel if there would be at least this length probe which may be InDel'd
const int cMaxMicroInDelLen  = 20;		// allow microInDels to be at most this length
const int cMaxMicroInDelMM  = 2;		// micro indels may only have at most this total number of mismatches
const int cMinJunctAlignSep  = 25;		// putative splice junctions if read segments separated by at least this distance
const int cDfltJunctAlignSep = 50000;	// default maximum separation
const int cMaxJunctAlignSep = 100000;	// absolute maximum separation
const int cMaxJunctAlignMM  = 2;		// splice junction alignments may only have at most this total number of mismatches
const int cMinJunctSegLen  = 10;		// splice junction segments must be at least this length

const int cMaxPutInDelOfss = 80;	    // absolute maximum number of putative InDels that can be explored in any probe

// scoring for alignments, InDels and splice junctions
const int cBaseScore = 500;			// base score before add/sub matches etc
const int cMaxScore = 1000;			// clamp scores to be no larger than this value
const int cScoreMatch = 3;			// score (add) matches as this
const int cScoreMismatch = 5;		// score (subtract) mismatches as this
const int cSpliceDonorAccept = 50; // score (add) for presence of splice donor+acceptor sites
const int cSpliceLen = 10;			// score (subtract) for every 1K of separation between splice junctions
const int cScoreInDelOpn = 20;		// score (subtract) opening an InDel as this
const int cScoreInDelExtn = 1;		// score (subtract) extending an opened InDel as this
const int cMaxAlignOpcodesLen = 256;  // limit alignment opcodes + length bytes to be at most this length

// PacBio processing
const int cMinPacBioSeedCoreLen = 8;							// user can specify seed cores down to this minimum length
const int cMaxPacBioSeedCoreLen = 100;							// user can specify seed cores up to to this maximum length
const int cMaxPacBioSeedExtn = 25;                              // if seed core < this limit then extend seed matches looking for near matches

const int cPacBioSeedCoreExtn = 100;							// looking for matches over this length seed core extension
const int cPacBiokExtnKMerLen = 4;								// matches in seed core extension must be of at least this length
const int cPacBioMinKmersExtn = 15;							    // require at least this many cPacBiokExtnKMerLen-mer matches over cPacBioSeedCoreExtn core extension

const int cMaxQualTargs = 10000;				// can process at most this many qualified targets

typedef enum etALStrand {
	eALSboth,							// align to both the watson and crick strand
	eALSWatson,							// align to the watson strand only
	eALSCrick,							// align to the crick strand only
	eALSnone							// align to neither strand
} eALStrand;

typedef enum etHRslt {
	eHRnone = 0,						// no change to that of previous search or no hits
	eHRhits,							// hits, MMDelta criteria met and within the max allowed number of hits
	eHRMMDelta,							// same or a new LowMMCnt unique hit but MMDelta criteria was not met
	eHRHitInsts,						// same or a new LowMMCnt but now too many multiple hit instances
	eHRRMMDelta						    // same but with a reduced MMDelta differential
} tHRslt;

#pragma pack(1)

// each entry for sequences is described by the following fixed size structure
typedef struct TAG_sSfxEntry {
	UINT32 EntryID;						// identifies each entry (1..n), unique within this suffix file
	UINT32 fBlockID;					// which SfxBlock contains sequence for this entry
	UINT8 szSeqName[cMaxDatasetSpeciesChrom];	// entry name - typically a chromosome name
	UINT16 NameHash;					// hash over szSeqName
	UINT32 SeqLen;						// sequence length - excludes sequence terminator eSeqEOS
	UINT64 StartOfs;					// offset into concatenated sequences of where this sequence starts
	UINT64 EndOfs;						// offset into concatenated sequences of where this sequence ends
} tsSfxEntry;


// file contains 0 or 1 entry blocks containing all entry descriptors
typedef struct TAG_sSfxEntriesBlock {
	UINT32 NumEntries;					// current number of entries
	UINT32 MaxEntries;					// can contain at most this number of entries
	tsSfxEntry Entries[1];				// entry descriptors
	} tsSfxEntriesBlock;

// file contains 0 or 1 SfxBlocks
typedef struct TAG_sSfxBlock {
	UINT32 BlockID;						// identifies (1..n) this suffix block within this file
	UINT32 NumEntries;					// number of entries contained in this block
	UINT64 ConcatSeqLen;				// total length (including terminators) of all concatenated sequences (is also the number of elements in suffix array)
	UINT32 SfxElSize;					// number of bytes per suffix array element = will be either 4 or 5
	UINT8 SeqSuffix[1];					// concatenated sequences followed by suffix array
	} tsSfxBlock;

// each entry for sequences is described by the following fixed size structure
typedef struct TAG_sSfxEntryV3 {
	UINT32 EntryID;						// identifies each entry (1..n), unique within this suffix file
	UINT32 fBlockID;					// which SfxBlock contains sequence for this entry
	UINT8 szSeqName[cMaxDatasetSpeciesChromV3];	// entry name - typically a chromosome name
	UINT16 NameHash;					// hash over szSeqName
	UINT32 SeqLen;						// sequence length - excludes sequence terminator eSeqEOS
	UINT64 StartOfs;					// offset into concatenated sequences of where this sequence starts
	UINT64 EndOfs;						// offset into concatenated sequences of where this sequence ends
} tsSfxEntryV3;

// file contains 0 or 1 entry blocks containing all entry descriptors
typedef struct TAG_sSfxEntriesBlockV3 {
	UINT32 NumEntries;					// current number of entries
	UINT32 MaxEntries;					// can contain at most this number of entries
	tsSfxEntryV3 Entries[1];				// entry descriptors
	} tsSfxEntriesBlockV3;


typedef struct TAG_sIdentNode {							// TargSeqIDs are hash linked into tsIdentNodes
	UINT32 TargSeqID;
	struct TAG_sIdentNode *pNxt;						// will be NULL if last linked in current hash chain
} tsIdentNode;

typedef struct TAG_sQueryAlignNodes {					
	UINT32 AlignNodeID;									// alignment node identifier
	UINT8 FlgStrand:1;									// 0 if aligning query sense, 1 if aligning query antisense
	UINT8 FlgFirst2tRpt:1;								// set 1 if this node determined to be the 1st node in path which meets minimum scoring critera
	UINT8 Flg2Rpt:1;									// set 1 if this node part of path which meets minimum scoring critera
	UINT8 FlgScored:1;									// set 0 if unscored, 1 if scored and HiScore, HiScorePathNextIdx are valid
	INT32 QueryID;										// query sequence identifer
	UINT32 QueryStartOfs;								// starting at this loci
	UINT32 TargSeqID;									// aligning to thid target sequence identifer
	UINT32 TargStartOfs;								// starting at this loci
	UINT32 AlignLen;									// alignment length
	UINT32 NumMismatches;								// aligns with this number of mismatching nucleotides
	UINT32 CurCoreSegOfs;								// query core used to seed this alignment node was from this sequence offset
	INT32 HiScore;										// highest cumulative score for nodes along highest scoring path starting with this node (includes this node's score)
	UINT32 HiScorePathNextIdx;							// if > 0 then idex-1 of next node on currently highest scoring path; 0 if no other nodes on path
} tsQueryAlignNodes;

typedef struct TAG_sQualTarg {
	UINT32 TargEntryID;		// identifies sequence
	UINT8  Hits;			// against which there are this many hits (clamped to be at most 255)
	} tsQualTarg;

#pragma pack()

// file can contain at most one suffix block, this block contains one suffix array with this array indexing multiple sequences
// suffix blocks need to be of a size such that the complete block can be loaded into memory on the targeted machine

const UINT32 cAllocSfxEntries = 500000;		// alloc/realloc for sequence entries in this many increments
const UINT32 cMaxSfxEntries = 100000000;		// can handle upto this number of sequence entries
const UINT32 cMaxAllowSeqLen = 0xfff00000;  // maximum length (nearly 4G) of any individual sequence allowed 
const UINT64 cMaxAllowConcatSeqLen = 1000000000000; // max supported concatenation length of all sequences (must fit within 40bits)
const UINT32 cReallocBlockEls = cMaxAllowSeqLen/10;  // minimum realloc for sfxblock elements
const UINT64 cThres8ByteSfxEls = 4000000000;  // if concatenated sequence length >= this threshold then use 5bytes per suffix element instead of 4 when creating suffix index

const int cMaxCultivars = 100;	// can handle at most this many different cultivars
const int cMinCultivarPreSufLen = 5;	// prefix or suffix length must be at least this many bases
const int cMaxCultivarPreSufLen = 200;  // prefix or suffix length must be no longer than this many bases
const int cTotCultivarKMerLen   = 300;	// prefix + suffix length combined no longer than this many bases

#pragma pack(4)

// fixed size V3 file header with cMaxDatasetSpeciesChrom increased to 81
typedef struct TAG_sSfxHeaderV3 {
	unsigned char Magic[4];			 		// magic chars to identify this file as a SfxArrayV2 file
	INT32 Version;							// file structure version
	UINT32 Attributes;						// file attributes - currently bit 0 bisulfite processing, bit 1 for colorspace
	UINT64 FileLen;							// current file length (write psn for newly created blocks)
	UINT64 EntriesOfs;						// offset in file to tsSfxEntryBlock, 0 if no entries
	UINT32 EntriesSize;						// allocation size for loading tsSfxEntriesBlock into memory from file
	UINT32 NumSfxBlocks;					// current number of tsSfxBlock blocks in file, either 0 or 1
	UINT64 SfxBlockSize;						// size of SfxBlock in file - used if preallocing memory to hold SfxBlock
	UINT64 SfxBlockOfs;						// file offset at which tsSfxBlock starts
	UINT8 szDatasetName[cMaxDatasetSpeciesChrom]; // dataset name - usually the genome species name
	UINT8 szDescription[cMBSFFileDescrLen];	// describes contents of file
	UINT8 szTitle[cMBSFShortFileDescrLen];	// short title by which this file can be distingished from other files in dropdown lists etc
} tsSfxHeaderV3;

// original fixed size V3 file header with cMaxDatasetSpeciesChrom set to 36
typedef struct TAG_sSfxHeaderVv {
	unsigned char Magic[4];			 		// magic chars to identify this file as a SfxArrayV2 file
	INT32 Version;							// file structure version
	UINT32 Attributes;						// file attributes - currently bit 0 bisulfite processing, bit 1 for colorspace
	UINT64 FileLen;							// current file length (write psn for newly created blocks)
	UINT64 EntriesOfs;						// offset in file to tsSfxEntryBlock, 0 if no entries
	UINT32 EntriesSize;						// allocation size for loading tsSfxEntriesBlock into memory from file
	UINT32 NumSfxBlocks;					// current number of tsSfxBlock blocks in file, either 0 or 1
	UINT64 SfxBlockSize;						// size of SfxBlock in file - used if preallocing memory to hold SfxBlock
	UINT64 SfxBlockOfs;						// file offset at which tsSfxBlock starts
	UINT8 szDatasetName[cMaxDatasetSpeciesChromV3]; // dataset name - usually the genome species name
	UINT8 szDescription[cMBSFFileDescrLen];	// describes contents of file
	UINT8 szTitle[cMBSFShortFileDescrLen];	// short title by which this file can be distingished from other files in dropdown lists etc
} tsSfxHeaderVv;

#pragma pack()

#pragma pack(1)

// returned multiple hit loci, note that can have 2 segments if microInDel of read spans splice junction
typedef struct TAG_sSegLoci {
	    UINT16 ReadOfs;				// Original match started at this read sequence offset
		UINT8 Strand;				// matched on this strand
		UINT32 ChromID;				// suffix entry (chromosome) matched on
		UINT64 MatchLoci;			// original match loci
		UINT16 MatchLen;			// original match length
		UINT8 Mismatches;			// original number of mismatches
		UINT16 TrimLeft;			// left flank trimming removes this many bases
        UINT16 TrimRight;			// right flank trimming removes this many bases
		UINT8 TrimMismatches;		// after trimming there are this many mismatches
}	tsSegLoci;

typedef struct TAG_tsHitLoci {
	 etSeqBase BisBase;				// bisulfite mapping base if appropriate
	 UINT8 FlgChimeric:1;			// set if chimeric alignment with TrimLeft and/or TrimRight flank trimming required
	 UINT8 FlgInDel:1;				// set if microInDel alignment with 2 segments
	 UINT8 FlgInsert:1;				// set if alignment with insertion into read relative to reference, false if deletion from read
	 UINT8 FlgSplice:1;				// set if splice junction with 2 segments
 	 UINT8 FlgNonOrphan:1;			// set if determined that at least two splice or InDels sharing intersegment start/ends
	 UINT16 Score;					// alignment match score
	 tsSegLoci Seg[2];				// both used if InDel or splice junction	
 } tsHitLoci;

typedef struct TAG_sKMerCultCnts {
	UINT32 EntryID;						// counts are for this sfx entry
	UINT32 SenseCnts;					// total number of K-Mers on the sense strand
	UINT32 AntisenseCnts;				// total number of K-Mers on the antisense strand
} tsKMerCultDist;

typedef struct TAG_sKMerCultsCnts {
	UINT64 SenseCnts;					// total number of K-Mers on the sense strand over all cultivars
	UINT64 AntisenseCnts;				// total number of K-Mers on the antisense strand over all cultivars
	UINT32 NumCultivars;				// number of cultivars with counts
	UINT32 TotCultivars;				// total number of cultivars with associated EntryIDs in CultCnts
	etSeqBase KMerSeq[cTotCultivarKMerLen + 1];	// the K-Mer sequence which potentially could contain both prefix + suffix of maximal length
	tsKMerCultDist CultCnts[cMaxCultivars];	// individual cultivar counts
} tsKMerCultsCnts;


// cigar sequence style opcodes
// The following opcodes are in bits 0..3 of the opcode byte with bit 7 set to indicate that this is an opcode
// Values associated with the opcodes are in subsequent bytes with bit 7 always reset
//
typedef enum TAG_eAlignOPcodes {
	eAOPempty = 0x080,					// used to mark that there are no more opcodes to be popped or returned 
	eAOPMatch = 0x081,					// probe bases exactly matches target sequence bases for length specified in following byte(s)
	eAOPMismatch = 0x082,				// probe base does not match target sequence base, probe base is in next byte bits 0..3 and targ sequence base is in bits 4..7
	eAOPInsert = 0x083,					// probe has inserted bases relative to target, length of insertion is specified in following byte(s)
	eAOPDelete = 0x084,					// probe has deleted bases relative to target, length of insertion is specified in following byte(s)
	eAOPSkip = 0x085,					// skip this number of bases in target, number of bases to skip is in specified in following byte(s)
} eAlignOPcodes;

class CAlignOpCode : public CErrorCodes, protected CEndian {
	UINT8 m_AlignOpCodes[cMaxAlignOpcodesLen];	// to hold alignment opcodes + associated values for lowest scored alignment
	int m_MaxAlignOpCodes;					// how many bytes of opcodes + associated values can be stored 
	int m_CurNumAlignOpCodes;				// current number of opcodes + associated values in m_pAlignOpCodes
public:
	CAlignOpCode(void);
	~CAlignOpCode(void);
	void Reset(void);							// clears m_CurNumAlignOpCodes
	void Clone(CAlignOpCode &CloneFrom);		// clone from
	bool Init(int Len,UINT8 *pOpCodes);			// initialise using these pre-existing opcodes
	bool TruncNthMM(int NthMM);						// truncate at the Nth mismatch
	bool Append(int Len,UINT8 *pOpCodes);
	bool Append(CAlignOpCode &ApendFrom);
	bool AppendNthMM(int NthMM,CAlignOpCode &AppendFrom); // first truncate at the Nth mismatch then append AppendFrom 
	eAlignOPcodes								// opcode popped, or -1 if errors
		PopAlignOpCode(UINT32 *pCounts);		// opcode's associated count or value
	bool										// success or otherwise
		PushAlignOpCode(eAlignOPcodes OpCode,	// opcode to push
				UINT32 Counts);					// opcode associated count or value
	bool										// success or otherwise
		PushAlignMismatch(etSeqBase ProbeBase,	// probe base
				  etSeqBase RefBase);			// target or reference base
};

#pragma pack()


class CSfxArrayV3 : public CErrorCodes, protected CEndian
{
	friend class CLocKMers;
	friend class CMarkerKMers;
	bool m_bInMemSfx;							// true if in-memory suffix processing only - no file I/O
	int m_hFile;							    // opened/created file handle
	bool m_bV3File;								// suffix file opened was a V3 file
	char m_szFile[_MAX_PATH+1];				    // file name as opened/created
	bool m_bHdrDirty;							// TRUE if header has been updated and should be written to disk
	bool m_bCreate;								// TRUE if file opened in create mode
	bool m_bBisulfite;							// TRUE if bisulfite processing
	bool m_bColorspace;							// TRUE if colorspace (SOLiD) processing

	int m_MaxIter;								// max allowed iterations (depth) per subsegmented sequence when matching that subsegment

	int m_MaxMMExploreInDel;					// if more than this number of mismatches then explore microInDels
	int m_MaxInDelLen;							// any microInDel InDel must be <= this length
	int m_MinInDelSeqLen;						// only explore for InDel if there is at least this length sequence which may be InDel'd

	static UINT8 CvrtSOLiDmap[5][7];
	
	int m_ReqBlockID;							// block requested to be loaded by background thread
	teBSFrsltCodes m_ReqBlockRslt;				// set by background processing thread with block loading result
	bool m_bTermThread;							// set true if background processing threads are to terminate
	bool m_bThreadActive;						// set true if any background processing threads have been started

	int m_MaxQSortThreads;						// max number of threads to use when sorting
	CMTqsort m_MTqsort;							// multithreaded qsort

	UINT32 m_MaxKMerOccs;						// if there are more than MaxKMerOccs instances of a Kmer then these will be classified as an over-occurance
	size_t m_AllocOccKMerClasMem;				// allocation memory size for m_pOccKMerClas 
	UINT8 *m_pOccKMerClas;						// to hold Kmer instance classifications packed 4 per byte
	int m_OccKMerLen;							// processing KMers of this length for over occurance against m_MaxKMerOccs


#ifdef _WIN32
static	unsigned __stdcall ThreadedPrereadBlocks(void * pThreadPars);
	HANDLE m_threadHandle;						// handle as returned by beginthreadex
	UINT32 m_threadID;							// identifier as set by beginthreadex
	HANDLE m_JobReqEvent;
	HANDLE m_JobAckEvent;						// used by background thread to ack the job request when job processed
	HANDLE m_JobMutex;							// used to serialise access by main and background thread to shared resources
#else
	int m_threadRslt;								// result as returned by pthread_create ()
	pthread_t m_threadID;							// identifier as set by pthread_create ()
	pthread_mutex_t m_JobMutex;					// used to serialise access by main and background thread to shared resources
	pthread_cond_t m_JobReqEvent;
	pthread_cond_t m_JobAckEvent;
static	void *ThreadedPrereadBlocks(void * pThreadPars);
#endif

	volatile unsigned int m_CASSeqFlags; // used with synchronous compare and swap (CAS) for serialising access to base flags


	void SerialiseBaseFlags(void);
	void ReleaseBaseFlags(void);

	UINT32	m_MaxAllowSeqLen;					// any individual sequence can be at most this length
	UINT64 m_MaxSfxBlockEls;					// limit suffix blocks to hold at most this many index elements

	tsSfxHeaderV3 m_SfxHeader;					// loaded suffix file header
	tsSfxEntriesBlock *m_pEntriesBlock;			// loaded entries block
	tsSfxBlock *m_pSfxBlock;					// loaded suffix block
	UINT64 m_AllocSfxBlockMem;					// memory allocation size for loaded suffix blocks
	UINT64 m_AllocEntriesBlockMem;				// memory allocation size for loaded entry block
	UINT64 m_AllocBisulfiteMem;					// memory allocation size for loaded bisulfite
    UINT8 *m_pBisulfateBases;					// used whilst constructing sfx array if bisulfite processing
	UINT64 m_EstSfxEls;							// estimate of how mant sfx array elements may be required when creating sfx array 

	teBSFrsltCodes ChunkedWrite(INT64 WrtOfs,UINT8 *pData,INT64 WrtLen); // handles WrtLen > INT_MAX
	teBSFrsltCodes ChunkedRead(INT64 RdOfs,UINT8 *pData,INT64 RdLen);  // handles RdLen > INT_MAX
	void InitHdr(void);
	teBSFrsltCodes Hdr2Disk(void);				// writes header to file
	teBSFrsltCodes Disk2Hdr(char *pszFile);		// loads header from file
	teBSFrsltCodes Disk2Entries(void);			// loads entries from file
	teBSFrsltCodes Entries2Disk(void);			// writes entries to file
	teBSFrsltCodes SfxBlock2Disk(void);			// writes sfx block to file
	teBSFrsltCodes Disk2SfxBlock(int BlockID);	// loads specified sfx block from file

	teBSFrsltCodes Flush2Disk(void);			// flush and commit to disk

	tsSfxEntry *MapChunkHit2Entry(UINT64 ChunkOfs); // Maps the chunk hit loci to the relevant sequence entry

	int TransformToColorspace(UINT8 *pSrcBases,	// basespace sequence to transform into colorspace (SOLiD)
					INT64 SeqLen,				// sequence length
					UINT8 *pDstColorspace,		// where to write colorspace
					bool bHiNibble=false,		// if false then write into bits 0..3, otherwise bits 4..7, of destination
					etSeqBase Seed=eBaseN);		// use this base as the seed base immediately preceding the first base

	int TransformToBasespace(UINT8 *pSrcColorspace,	// colorspace (SOLiD) sequence to transform
					INT64 SeqLen,				// sequence length
					UINT8 *pDstBases,			// where to write basespace
					bool bHiNibble=false,		// if false then write into bits 0..3, otherwise bits 4..7, of destination
					etSeqBase Seed=eBaseN);		// use this base as the seed base immediately preceding the first base

	bool AcceptMasked(etSeqBase *pTarg,etSeqBase *pProbe,int MatchLen); // determine if repeat masked within limits

	int										// number (upto Lim) of non-canonical bases in pSeq 
			NumNonCanonicals(int Lim, // process for at most this many non-canonical bases
					UINT32 SeqLen,		// pSeq is of this max length
					etSeqBase *pSeq);	// pSeq to process for count of non-canonical

	int CmpProbeTarg(etSeqBase *pEl1,etSeqBase *pEl2,int Len);
	int BSCmpProbeTarg(etSeqBase *pEl1,etSeqBase *pEl2,int Len);

	int											// length of exact match
		ExactMatchLen(etSeqBase *pProbe,		// determine exactly matching length between probe
							etSeqBase *pTarg,	// and target sequences
							int MaxMatchLen);	// for up to this length

	etSeqBase GetBisBase(int TargMatchLen,etSeqBase *pTargBase,etSeqBase *pProbeBase);

	etSeqBase *GetPtrSeq(int EntryID,UINT32 Loci);

	INT64			// index+1 in pSfxArray of first exactly matching probe or 0 if no match
		LocateFirstExact(etSeqBase *pProbe,  // pts to probe sequence
				  int ProbeLen,					// probe length to exactly match over
				  etSeqBase *pTarg,				// target sequence
				  int SfxElSize,				// size in bytes of suffix element - expected to be either 4 or 5
				  void *pSfxArray,				// target sequence suffix array
				  INT64 TargStart,				// position in pTarg (0..n) corresponding to start of suffix array
				  INT64 SfxLo,					// low index in pSfxArray
				  INT64 SfxHi);					// high index in pSfxArray

	INT64			// index+1 in pSfxArray of last exactly matching probe or 0 if no match
		LocateLastExact(etSeqBase *pProbe, // pts to probe sequence
				  int ProbeLen,					// probe length to exactly match over
				  etSeqBase *pTarg,				// target sequence
				  int SfxElSize,				// size in bytes of suffix element - expected to be either 4 or 5
				  void *pSfxArray,				// target sequence suffix array
				  INT64 TargStart,				// position in pTarg (0..n) corresponding to start of suffix array
				  INT64 SfxLo,					// low index in pSfxArray
				  INT64 SfxHi);					// high index in pSfxArray


	int									// -1 if errors, 0 if no match, 1 if mismatches only, 2 InDel and insertion into probe, 2 InDel and deletion from probe
		ExploreInDelMatchRight(UINT32 ExtdProcFlags,	// flags indicating if lower levels need to do any form of extended processing with this specific read...
			UINT32 ReadID,				// identifies this read
			char CurStrand,				// aligning on this strand	
		    int microInDelLen,		   	// microInDel can be upto (inclusive) this length
			int MaxTotMM,				// can be at most this many mismatches in total
 		    int ProbeLen,				// length of probe excluding any eBaseEOS
		    etSeqBase *pProbe,			// pts to 5' start of probe sequence
		    tsSfxEntry *pEntry,			// target at TargIdx is in this sfx entry
		    INT64 TargIdx,				// pTarg corresponds to this suffix index
		    etSeqBase *pTarg,			// pts to 5' start of target sequence
		    tsHitLoci *pHit);			// where to return hit loci	

		int									// -1 if errors, 0 if no match, 1 if mismatches only, 2 InDel and insertion into probe, 2 InDel and deletion from probe
		ExploreInDelMatchLeft(UINT32 ExtdProcFlags,			// flags indicating if lower levels need to do any form of extended processing with this specific read...
			UINT32 ReadID,				// identifies this read
			char CurStrand,				// aligning on this strand	
		    int microInDelLen,		   	// microInDel can be upto (inclusive) this length
			int MaxTotMM,				// can be at most this many mismatches in total
 		    int ProbeLen,				// length of probe excluding any eBaseEOS
		    etSeqBase *pProbe,			// pts to 5' start of probe sequence
			tsSfxEntry *pEntry,		    // target at TargIdx is in this sfx entry
		    INT64 TargIdx,				// pTarg corresponds to this suffix index
		    etSeqBase *pTarg,			// pts to 5' start of target sequence
		    tsHitLoci *pHit);			// where to return hit loci	

		int									// -1 if errors, 0 if no match, 1 if mismatches only, 2 InDel and insertion into probe, 3 InDel and deletion from probe
		ExploreSpliceRight(char CurStrand,					// aligning on this strand
			int MaxSpliceJunctLen,		// junction has to be no longer than this length
			int MaxTotMM,			// can be at most this many mismatches in total
			int CoreLen,			// core length used
 		   int ProbeLen,			// length of probe excluding any eBaseEOS
		   etSeqBase *pProbe,		// pts to 5' start of probe sequence
		   INT64 TargOfs,			// pTarg corresponds to this suffix offset
		   INT64 TargLen,			// max length of target to be explored 
		   etSeqBase *pTarg,		// pts to target sequence, 5' 1st base if bRight, 3' last base if !bRight
		   tsHitLoci *pHit);		// where to return hit loci	

		int									// -1 if errors, 0 if no match, 1 if mismatches only, 2 InDel and insertion into probe, 3 InDel and deletion from probe
			ExploreSpliceLeft(char CurStrand,					// aligning on this strand
			int MaxSpliceJunctLen,		// junction has to be no longer than this length
			int MaxTotMM,			// can be at most this many mismatches in total
			int CoreLen,			// core length used
 		   int ProbeLen,			// length of probe excluding any eBaseEOS
		   etSeqBase *pProbe,		// pts to 5' start of probe sequence
		   INT64 TargOfs,			// pTarg corresponds to this suffix offset
		   INT64 TargLen,			// max length of target to be explored 
		   etSeqBase *pTarg,		// pts to target sequence, 5' 1st base if bRight, 3' last base if !bRight
		   tsHitLoci *pHit);		// where to return hit loci

		INT64												// returned number of antisense K-Mers identified
			AntisenseKMerCultsCnts(int PrefixKMerLen,			// report on K-Mers having this prefix sequence length
			   int SuffixKMerLen,		// and allow for the K-mers containing suffixes of this length (can be 0) 
			  tsKMerCultsCnts *pKMerCultsCnts); // update these cultivar counts

public:
	CSfxArrayV3(void);
	~CSfxArrayV3(void);
	int Reset(bool bFlush = true);			// reset state back to that immediately following instantiation
	int Open(char *pszSeqFile,				// specifies file to open or create
			   bool bCreate = false,		// create file if it doesn't already exist, truncate if it does
			   bool bBisulfite = false,		// if true then bisulfite processing
               bool bColorspace = false);	// if true then colorspace (SOLiD) processing

	// Memory resident suffix array
	// User creates the suffix array in memory and can immediately access for alignments etc without needing to write to file
	int Open(bool bBisulfite = false,				// if true then bisulfite processing
					bool bColorspace = false);		// if true then colorspace (SOLiD) processing


	int Close(bool bFlush = true);			// closes opened file

    // obtain a copy of the header for external diagnostics
	tsSfxHeaderV3 *GetSfxHeader(tsSfxHeaderV3 *pCopyTo);

	int SetMaxIter(int MaxIter);			// set maximum iterations on identically matching subsequences
	int GetMaxIter(void);					// get maximum iterations on identically matching subsequences

	teBSFrsltCodes
		InitOverOccKMers(int KMerLen,			// will be processing for over occuring KMers of this length (max cMaxKmerLen)
					int MaxKMerOccs);			// which if there are more than MaxKMerOccs instances will be classified as an over-occurance

	int											// 0: no instances, 1: number occurances <= m_MaxKMerOccs, 2: number occurances > m_MaxKMerOccs
		OverOccKMerClas(int KMerLen,			// KMer length, checked and must be equal to m_OccKMerLen
					etSeqBase *pSeq);			// return over occurance classification for this sequence

	void SetInitalSfxAllocEls(INT64 NumEls);	// estimated number of elements to allocate when creating suffix array - hint only!

	bool IsSOLiD(void);						// returns true if index created in colorspace (SOLiD)
	teBSFrsltCodes SetDatasetName(char *pszDataset);	// sets file dataset name
	char *GetDatasetName(void);						// returns file dataset name
	teBSFrsltCodes SetDescription(char *pszDescription);
	teBSFrsltCodes GetDescription(int MaxLen,char *pszDescription);
	teBSFrsltCodes SetTitle(char *pszTitle);
	teBSFrsltCodes GetTitle(int MaxLen,char *pszTitle);

	teBSFrsltCodes	AddEntry(char *pszSeqIdent,	// sequence identifier, typically chromosome name
				etSeqBase   *pSeq,				// sequence to add to suffix array
				UINT32 SeqLen,					// sequence length excluding any eBaseEOS
				UINT16 Flags = 0);				// user specified flags


	int Finalise(void);						// finalise and commit to disk the suffix array after all entries added with AddEntry()

	int										// if non-zero then returned number of identifiers 
		ChkDupEntries(int MaxIdents,		// maximum number of identifers to return in pIdents (caller allocates to hold returned identifiers)
					  UINT32 *pIdents);		// checks if there are duplicate entry names and reports identifier

	int										// returns number of matching sequences
		Locate(etSeqBase *pProbe,			// sequence containing probe
				  UINT32 ProbePsn,			// start position in probe sequence at which to search
				  etSeqBase *pTarg,			// target sequence
				  UINT32 *pSfxArray,		// target sequence suffix array
				  UINT32 TargStart,			// position in pTarg (0..n) corresponding to start of suffix array
				 UINT32 SfxLen);			// number of suffixs in pSfxArray

	INT64						// returned number of cultivar sequences
		LocateMultiCultivarMarkers(int PrefixKMerLen,	// prefix K-mer length
						   int SuffixMerLen,			// prefix K-mer length
						   int NumCultivars,			// number of cultivars in pCultivars
						   UINT32 *pEntryIDs,			// sequences of interest will be on these
						   void *pThis,					// class instance for callback
						   int Callback(void *pThis,UINT32 EntryID,etSeqBase *pSeq)); // callback on cultivar marker sequences

	int												// < 0 if errors, 0 if no processing required by thread, 1 if *pSfxIdxRange initilised to number of elements to process 
		GenKMerCultThreadRange(int KMerLen,	// processing K-Mers of this length
					  int ThreadInst,	// starting and range required for this thread instance
					  int NumThreads,				// total number of threads which will be used for generating cultivar K-Mers
					  INT64 StartSfxIdx,			// given starting suffix index (0..N) for the thread instance
					  INT64 *pEndSfxIdx);			// returned end suffix (inclusive) index to be processed by this thread

	INT64						// < 0 if errors, otherwise the total number of identified K-Mers over both strands (could be more than reported if MinCultivars > 1
		GenKMerCultsCnts(bool bSenseOnly,				// true if sense strand only processing, default is to process both sense and antisense
						INT64 StartSfxIdx,				// starting suffix index
							INT64 EndSfxIdx,			// finish processing at this suffix index, inclusive - if 0 then process all remaining
							int PrefixKMerLen,			// report on K-Mers having this prefix sequence length
						    int SuffixKMerLen,			// and allow for the K-mers containing suffixes of this length (can be 0) 
							int MinCultivars,			// only report if K-Mers present in at least this many different cultivars (0 if must be present in all cultivars)
							int MaxHomozygotic,			// only report prefixes if K-Mer suffixes are homozygotic between a maximum of this many cultivars, if 1 then no other cultivars
							void *pThis,				// callers class instance
						    int Callback(void *pThis,tsKMerCultsCnts *pCultsCnts)); // callback on K-Mers


			int											// < 0 if errors, 0 if no matches or change, 1 if mumber of matches accepted, 2 MMDelta criteria not met, 3 too many match instances  
				AlignReads(UINT32 ExtdProcFlags,		// flags indicating if lower levels need to do any form of extended processing with this specific read...
						 UINT32 ReadID,					// identifies this read
						 int MinChimericLen,			// if checking for chimerics then minimim length required, set to 0 if not checking for chimerics
						 int TotMM,						// max number of mismatches allowed including any core
						 int CoreLen,					// core window length 
						 int CoreDelta,					// core window offset increment (1..n)
						 int MaxNumCoreSlides,			// limit number of core slides to this many per strand over the read length
						 int MinCoreLen,				// do not reduce core below this length
						 int MMDelta,					// minimum (1..n) mismatch difference between the best and next best core alignment
						 eALStrand Align2Strand,		// align to this strand
  						 int microInDelLen,				// microInDel length maximum
						 int MaxSpliceJunctLen,			// junction has to be no longer than this length 
  						 int *pLowHitInstances,			// In/Out number of match instances for lowest number of mismatches thus far for this read
						 int *pLowMMCnt,				// In/Out lowest number of mismatches thus far for this read
						 int *pNxtLowMMCnt,				// In/Out next to lowest number of mismatches thus far for this read
						 etSeqBase *pProbeSeq,			// probe
						 int ProbeLen,					// probe length
						 int MaxHits,					// (IN) process for at most this number of hits
						 tsHitLoci *pHits,				// where to return hits (at most MaxHits)
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes); // memory allocated by caller for holding tsIdentNodes

			int						// < 0 if errors, 0 if no alignmentse, 1 if at least one alignment
				IfAnyAlignments(UINT32 ReadID,					// identifies this read
				int MaxTotMM,			        // max number of mismatches allowed
				int CoreLen,					// core window length
				int CoreDelta,					// core window offset increment (1..n)
				int MaxNumCoreSlides,			// max number of times to slide core on each strand
				eALStrand Align2Strand,			// align to this strand
				etSeqBase *pProbeSeq,			// probe sequence
				int ProbeLen,					// probe length
				int CurMaxIter,					// max allowed iterations per subsegmented sequence when matching that subsegment
				int NumAllocdIdentNodes,		// memory has been allocated by caller for holding up to this many tsIdentNodes
				tsIdentNode *pAllocsIdentNodes); // memory allocated by caller for holding tsIdentNodes



			// Returns true if sequence is to be filtered out either because it contains a non-canonical base or the percentage of repeat masked bases is too high
			bool FiltOutRepeats(etRPTMasking RPTMasking,	// how to interpret cRptMskFlg'd bases
						 UINT32 MaskPercent,		// filter out if percentage of repeats is above this percentage (0-100)
						 etSeqBase *pSeq,			// pts to sequence
						 UINT32 SeqLen);				// sequence length


			INT64	BSIterateExacts(etSeqBase *pProbeSeq,	// bisulfite probe
 						 UINT32 ProbeLen,			// probe length
						 INT64 PrevHitIdx,			// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
 						 UINT32 *pTargEntryID,		// if match then where to return suffix entry (chromosome) matched on
						 UINT32 *pHitLoci);			// if match then where to return loci

			 INT64 IterateExacts(etSeqBase *pProbeSeq,// probe
 						 UINT32 ProbeLen,			// probe length
						 INT64 PrevHitIdx,			// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
 						 UINT32 *pTargEntryID,		// if match then where to return suffix entry (chromosome) matched on
						 UINT32 *pHitLoci,			// if match then where to return loci
						 UINT32 *pTargSeqLen = NULL, // optionally update with the matched target sequence length				
						 etSeqBase **ppTargSeq = NULL);	 // optionally update with ptr to start of the target sequence, exact match will have started at &ppTargSeq[*pHitLoci]  

			 INT64 IterateExacts(etSeqBase *pProbeSeq,// probe
 						 UINT32 ProbeLen,			// probe length
						 INT64 PrevHitIdx,			// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
 						 UINT32 *pTargEntryID,		// if match then where to return suffix entry (chromosome) matched on
						 UINT32 *pHitLoci,			// if match then where to return loci
						 UINT32 MaxExtdLen,			// attempt to extend exact match of ProbeLen out to MaxExtdLen and report extended match length in *pHitExtdLen
						int *pHitExtdLen);			// was able to extend core out to this length

			int											// number of target sequences which were prequalified
				PreQualTargs(UINT32 ProbeEntryID,		// probe sequence entry identifier used to detect self hit which are sloughed
						int ProbeSeqLen,				// probe sequence length
						etSeqBase *pProbeSeq,			// prequalify with cores from this probe sequence
						int QualCoreLen,				// prequalifying core length
						int DeltaCoreOfs,				// offset core windows of QualCoreLen along the probe sequence when checking for overlaps
						int TargLenDiffBp,				// 0 if disabled, accepted prequalified targets must have length within this differential (bp) of probe sequence length (used in transcriptome processing)
						int MaxPreQuals,				// max number of target sequences to pre-qualify
						tsQualTarg *pQualTargs);		// into this prequalified list

			INT64						// returned hit idex (1..n) or 0 if no hits
				IteratePacBio(etSeqBase *pProbeSeq,				// probe sequence
 										UINT32 ProbeLen,		// probe sequence length
 									 UINT32 SeedCoreLen,		// using this seed core length
									 UINT32 SloughEntryID,		// if > 0 then hits to this entry (normally would be the probe sequence entry identifier) are to be sloughed
									 UINT32 MinTargLen,			// accepted hit target sequences must be at least this length
									 UINT32 MaxTargLen,			// if > 0 then accepted hit target sequences must be no longer than this length
									 INT64 PrevHitIdx,			// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
 									 UINT32 *pTargEntryID,		// if match then where to return suffix entry (chromosome) matched on
									 UINT32 *pHitLoci,			// if match then where to return loci
									  int NumPreQuals = 0,			// number of pre-qualified sequences in pQualTargs
									tsQualTarg *pQualTargs = NULL,	// holds prequalified target sequence identifiers
									 int PacBioMinKmersExtn = cPacBioMinKmersExtn);		// accepting as putative overlap if extension matches at least this many cPacBiokExtnKMerLen (currently 4bp)	

			int                  // estimated identity (0..100) of overlap between pProbe and pTarg derived by dividing final score by the MatchScore and then normalising by length
				QuickScoreOverlap(int SeedCoreLen,          // initial seed core exact matching length between probe and target was this many bp 
									int SeqLen,				// quick score over this many remaining bases in pProbe and pTarg following the exactly matching SeedCoreLen
									etSeqBase *pProbe,		// probe sequence scoring onto
									etSeqBase *pTarg,		// this target sequence
								    int MatchScore = 3,     // exact match score  ((in any Pacbio sequence then expecting ~85% of all base calls to be correct), note that between any 2 sequences then the relative error rate doubles!
									int InDelPenalty = 2,  // insertions and deletions much more likely than substitutions (in any Pacbio sequence then expecting ~14% of all error events to be insertions)
							        int SubstitutePenalty = 7); // base call subs are relatively rare (in any Pacbio sequence then expecting ~1% of all error events to be substitutions)

			int												// much quicker, but far less accurate, than the QuickScoreOverlap SW banded function (above) as scoring only based on the number of 4-mers matching within a 16bp window along the two overlaping sequences
				QuickScoreOverlap(int SeqLen,				// both probe and target are of this minimum length (must be at least 16bp, will be clamped to be no more than 100bp)
							  etSeqBase *pProbe,			// scoring overlap of probe sequence onto
							  etSeqBase *pTarg);		    // this target sequence


			 int									// < 0 if errors, 0 if no matches to any other chroms other than ProbeChromID allowing up to MaxTotMM, 1 if matches to any other chrom
				MatchesOtherChroms( int ProbeChromID,	// probe is from this chrom
				    int MaxTotMM,		// allow for at most this number of missmatches
				    int ProbeLen,		// probe sequence is this length
					etSeqBase *pProbe);	// check for matches from this probe sequence

			 int									// < 0 if errors, 0 if no matches to any other chroms not in pProbeChromIDs[] allowing up to MaxTotMM, 1 if matches to any other chrom
				MatchesOtherChroms(int NumProbeIDs, // number of chroms in pProbeChromIDs
					int *pProbeChromIDs,	// probe is from one of these chroms
				    int MaxTotMM,		// allow for at most this number of missmatches
				    int ProbeLen,		// probe sequence is this length
					etSeqBase *pProbe);	// check for matches from this probe sequence


			 int						// < 0 if errors, 0 if no matches, 1 if a unique match, 2 if multiple matches
				LocateApproxUniques(int AccumReadHits,	// how many reads have already been matched, must be 0 or 1 
   						etSeqBase *pProbeSeq,			// probe
 						 etSeqBase *pPatternSeq,		// contains pattern to match with, etBaseN represents wildcard bases to match against any in the target
														// will be updated on return with etBaseN's changed to actual subsitution bases
						 UINT32 ProbeLen,				// probe, and also pattern, length
						 int ExpMismatches,				// expected pattern wildcard matches
 						 UINT32 *pTargEntryID,			// if match then where to return suffix entry (chromosome) matched on
						 UINT32 *pHitLoci,				// if unique match then where to return loci
						int CurMaxIter);				// max allowed iterations per subsegmented sequence when matching that subsegment

	int						// < 0 if errors, 0 if no matches or change, 1 if mumber of matches accepted, 2 MMDelta criteria not met, 3 too many match instances  
		LocateCoreMultiples(UINT32 ExtdProcFlags,		// flags indicating if lower levels need to do any form of extended processing with this specific read...
						 UINT32 ReadID,					// identifies this read
						 int MinChimericLen,			// if checking for chimerics then minimim length required, set to 0 if not checking for chimerics
						 int TotMM,						// max number of mismatches allowed including any core
						 int CoreLen,					// core window length 
						 int CoreDelta,					// core window offset increment (1..n)
						 int MaxNumCoreSlides,			// max number of times to slide core on each strand
						 int MMDelta,					// minimum (1..n) mismatch difference between the best and next best core alignment
 					 	 eALStrand Align2Strand,		// align to this strand
						 int *pLowHitInstances,			// In/Out number of match instances for lowest number of mismatches thus far for this read
						 int *pLowMMCnt,				// In/Out lowest number of mismatches thus far for this read
						 int *pNxtLowMMCnt,				// In/Out next to lowest number of mismatches thus far for this read

						 etSeqBase *pProbeSeq,			// probe
						 int ProbeLen,					// probe length
						 int MaxHits,					// (IN) process for at most this number of hits
						 tsHitLoci *pHits,				// where to return hits
						 int CurMaxIter,				// max allowed iterations per subsegmented sequence when matching that subsegment
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes); // memory allocated by caller for holding tsIdentNodes

		int						// < 0 if errors, 0 if no matches or change, 1 if mumber of matches accepted, 2 MMDelta criteria not met, 3 too many match instances
			LocateQuerySeqs(UINT32 QuerySeqID,			// identifies this query sequence
  						 etSeqBase *pProbeSeq,			// probe
						 UINT32 ProbeLen,				// probe length
						 UINT32 m_MaxMMThres,			// seed extensions terminate if the mismatches score is increased above this threshold; when extending seed matches are scored with 1 and mismatches scored with 2
						 UINT32 CoreLen,				// core window length
						 UINT32 CoreDelta,				// core window offset increment (1..n)
						 eALStrand Align2Strand,		// align to this strand
						 UINT32 MinMatchLen,			// putative alignments must be at least this length
						 UINT32 MaxHits,				// (IN) process for at most this number of hits
						 tsQueryAlignNodes *pHits,		// where to return hits (at most MaxHits)
						 UINT32 CurMaxIter);			// max allowed iterations per subsegmented sequence when matching that subsegment

			
	int						// < 0 if errors, 0 if no matches, 1..MaxHits, or MaxHits+1 if additional matches have been sloughed
		LocateBestMatches(UINT32 ReadID,			// identifies this read
						 int MaxTotMM,			        // return matches having at most this number of mismatches
						 int CoreLen,					// core window length
						 int CoreDelta,					// core window offset increment (1..n)
						 int MaxNumCoreSlides,			// max number of times to slide core on each strand
						 eALStrand Align2Strand,		// align to this strand
 						 etSeqBase *pProbeSeq,			// probe
						 int ProbeLen,					// probe length
						 int MaxHits,					// return at most this number of best hits
						 int *pHitInstances,			// returned number of match instances in pHits
						 tsHitLoci *pHits,				// where to return best hits (at most MaxHits)
						 int CurMaxIter,				// max allowed iterations per subsegmented sequence when matching that subsegment
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding up to this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes); // memory allocated by caller for holding tsIdentNodes


	int						// < 0 if errors or if probe sequence non-unique local to ProbeEntry, 0 if no matches, otherwise the total match instances  (see genzygosity.cpp)
		LocateAllNearMatches(int MaxTotMM,				// max number of mismatches allowed
						 int CoreLen,					// core window length 
						 int CoreDelta,					// core window offset increment (1..n)
						 int MaxNumCoreSlides,			// max number of times to slide core on each strand
						 int MaxMatches,				// if more than 0 and encountering more matches than this limit then return match count 0
 					 	 int NumEntries,				// number of entries in pEntryMatch
						 UINT32 *pEntryMatches,			// return number of matches to each entry in this array
						 UINT32 ProbeEntry,				// probe was from this entry
						 UINT32 ProbeOfs,				// starting at this offset
						 int ProbeLen,					// probe length
						 etSeqBase *pProbeSeq,			// probe
						 int CurMaxIter,				// max allowed iterations per subsegmented sequence when matching that subsegment
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes); // memory allocated by caller for holding tsIdentNodes

	int						// < 0 if errors, eHRnone if none or too many matches or change, eHRhits if mumber of matches accepted  
		LocateInDels(  UINT32 ExtdProcFlags,			// flags indicating if lower levels need to do any form of extended processing with this specific read...
						 UINT32 ReadID,					// identifies this read
						 int microInDelLen,				// microInDel length maximum
						 int MaxTotMM,					// max number of mismatches allowed
						 int CoreLen,					// core window length 
						 eALStrand Align2Strand,		// align to this strand
  						 int *pLowHitInstances,			// In/Out number of match instances for lowest number of mismatches thus far for this read
						 int *pLowMMCnt,				// In/Out lowest number of mismatches thus far for this read
						 int *pNxtLowMMCnt,				// In/Out next to lowest number of mismatches thus far for this read
						 etSeqBase *pProbeSeq,			// probe
						 int ProbeLen,					// probe length
						 int MaxHits,					// (IN) process for at most this number of hits
						 tsHitLoci *pHits,				// where to return any hit
 						 int *pScore,					// where to return InDel score
						 int CurMaxIter,				// max allowed iterations per subsegmented sequence when matching that subsegment
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes); // memory allocated by caller for holding tsIdentNodes

		int						// < 0 if errors, eHRnone if no matches or too many, eHRhits if mumber of matches accepted  
		LocateSpliceJuncts(UINT32 ExtdProcFlags,			// flags indicating if lower levels need to do any form of extended processing with this specific read...
						 UINT32 ReadID,					// identifies this read
						 int SpliceJunctLen,		// junction has to be no longer than this length 
						 int MaxTotMM,					// max number of mismatches allowed
						 int CoreLen,					// core window length 
						 eALStrand Align2Strand,		// align to this strand
  						 int *pLowHitInstances,			// In/Out number of match instances for lowest number of mismatches thus far for this read
						 int *pLowMMCnt,				// In/Out lowest number of mismatches thus far for this read
						 int *pNxtLowMMCnt,				// In/Out next to lowest number of mismatches thus far for this read
						 etSeqBase *pProbeSeq,			// probe
						 int ProbeLen,					// probe length
						 int MaxHits,					// (IN) process for at most this number of hits
						 tsHitLoci *pHits,				// where to return any hit
 						 int *pScore,					// where to return score
						 int CurMaxIter,				// max allowed iterations per subsegmented sequence when matching that subsegment
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes); // memory allocated by caller for holding tsIdentNodes


		
		int									// -1 if errors, 0 if no match, 1 if mismatches only, 2 InDel and insertion into probe, 3 InDel and deletion from probe
		AlignPairedRead(bool b3primeExtend,	// if false extend towards 5' of targeted chrom, if true extend towards 3' of targeted chrom
					bool bAntisense,		// if false extend with read sense, if true extend with read RevCpl'd (antisense)
					UINT32 ChromID,		  // accepted aligned read was on this chromosome
					UINT32 StartLoci,	  // accepted aligned read started at this loci
					UINT32 EndLoci,		  // and ending at this loci
					int MinDistance,	  // expecting partner to align at least this distance away from accepted aligned read (inclusive of read lengths)
					int MaxDistance,	  // but no more than this distance away (inclusive of read lengths)
					int MaxAllowedMM,	  // any accepted alignment can have at most this many mismatches
					int MinHamming,		  // and must be at least this Hamming away from the next best putative alignment 
					int ReadLen,		  // length of read excluding any eBaseEOS
					etSeqBase *pRead,	  // pts to 5' start of read sequence
					tsHitLoci *pAlign);	  // where to return any paired read alignment loci	

		int						// < 0 if errors, otherwise minimum hamming of probe to target
		LocateSfxHammings(int RHamm,					// restricted hammings limit
                         bool bSAHammings,				// if true then hammings on both sense and antisense required
						 int KMerLen,					// hammings for k-mers of this length
						 int SampleNth,					// sample every Nth k-mer
						 int SrcSeqLen,					// generate hammings over sequence of this length
						 UINT32 MinCoreDepth,			// initially explore cores to at least this depth and only continue exploring if depth would be no more than MaxCoreDepth
						 UINT32 MaxCoreDepth,			// explore cores to at most this depth
						 UINT32 EntryID,  				// KMer sequences are to be from this suffix entry
						 UINT32 EntryLoci,				// and starting at this loci
				 		int IntraInterBoth,	    // 0: hammings over both intra (same sequence as probe K-mer drawn from) and inter (different sequences to that from which probe K-mer drawn), 1: Intra only, 2: Inter only
						 UINT8 *pHammings,				// returned hammings
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes); // memory allocated by caller for holding tsIdentNodes

		int						// < 0 if errors, otherwise minimum hamming of probe to target
		LocateHammings(int RHamm,					// restricted hammings limit
                         bool bSAHammings,				// if true then hammings on both sense and antisense required
						 int KMerLen,					// hammings for k-mers of this length
						 int SampleNth,					// sample every Nth k-mer
						 int SrcSeqLen,					// sequence containing k-mers is of this length
						 UINT32 MinCoreDepth,		// initially explore cores to at least this depth and only continue exploring if depth would be no more than MaxCoreDepth
						 UINT32 MaxCoreDepth,		// explore cores to at most this depth
						 etSeqBase *pSrcSeq,			// hamming k-mers from this sequence
						 UINT32 SeqEntry,  				// sequence was from this suffix entry (0 if probe not from same assembly)
						 UINT32 SeqLoci,				// and starting at this loci
						 int IntraInterBoth,	    // 0: hammings over both intra (same sequence as probe K-mer drawn from) and inter (different sequences to that from which probe K-mer drawn), 1: Intra only, 2: Inter only
						 UINT8 *pHammings,				// returned hammings
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes); // memory allocated by caller for holding tsIdentNodes

		int						// < 0 if errors, otherwise minimum hamming of probe to target
		LocateHamming(int RHammMin,						// process for Hammings of at least this minimum
						 int RHammMax,					// process for Hammings upto this limit
						 UINT32 MinCoreDepth,			// initially explore cores to at least this depth and only continue exploring if depth would be no more than MaxCoreDepth
						 UINT32 MaxCoreDepth,			// explore cores to at most this depth
						 bool bSAHammings,				// if true then hammings relative to either/both sense and antisense required
 						 etSeqBase *pProbeSeq,			// probe sequence
						 UINT32 ProbeLen,				// probe length
 						 UINT32 ProbeEntry,  			// probe was from this suffix entry (0 if probe not from same assembly)
						 UINT32 ProbeLoci,				// and starting at this loci
						 int IntraInterBoth,	    // 0: hammings over both intra (same sequence as probe K-mer drawn from) and inter (different sequences to that from which probe K-mer drawn), 1: Intra only, 2: Inter only
 						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes); // memory allocated by caller for holding tsIdentNodes

		//
		int						// < 0 if errors, 0 if no matches or change, 1 if mumber of matches accepted, 2 MMDelta criteria not met, 3 too many match instances  
		LocateHammingX(int RHamm,					// restricted hammings limit
						 UINT32 MinCoreDepth,		// initially explore cores to at least this depth and only continue exploring if depth would be no more than MaxCoreDepth
						 UINT32 MaxCoreDepth,		// explore cores to at most this depth
                         bool bSAHammings,				// if true then hammings on both sense and antisense required
						 etSeqBase *pProbeSeq,			// probe sequence
						 UINT32 ProbeLen,				// probe length
 						 UINT32 ProbeEntry,  			// probe was from this suffix entry (0 if probe not from same assembly)
						 UINT32 ProbeLoci,				// and starting at this loci
						 int NumAllocdIdentNodes,		// memory has been allocated by caller for holding upto this many tsIdentNodes
						 tsIdentNode *pAllocsIdentNodes); // memory allocated by caller for holding tsIdentNodes

	int Next(int PrevBlockID = 0);						// iterates over block identifiers 
	int SetTargBlock(int BlockID);						// specifies which block is to be loaded for subsequent sequence matches
	int GetNumEntries(void);							// returns number of entries

	void InitAllIdentFlags(UINT16 Flags = 0);			// initialise all entries to have Flags value
	UINT16 GetIdentFlags(UINT32 EntryID);				// returns flags for EntryID
	UINT16												// returned flags prior to this call
		SetResetIdentFlags(UINT32 EntryID,				// atomic set/reset flags for this entry
					UINT16 SetFlags = 0,				// set these flags then
					UINT16 ResetFlags = 0);				// reset these flags

	int GetIdentName(UINT32 EntryID,int MaxLen,char *pszSeqIdent); // get sequence name for specified entry identifier
    int GetIdent(char *pszSeqIdent);					// returns identifier for specified sequence name
	UINT32 GetSeqLen(UINT32 EntryID);						// get length of sequence for requested entry identifier
	UINT64 GetTotSeqsLen(void);							// returns total length of all sequences in currently loaded entries block
	UINT32 GetMinSeqLen(void);							// returns minimum length of any sequence in currently loaded entries block
	UINT32 GetMaxSeqLen(void);							// returns maximum length of any sequence in currently loaded entries block

	
	int GetBaseFlags(UINT32 EntryID,						// identifies sequence containing loci flags to be returned - flags returned are in bits 0..3
				UINT32 Loci);							// offset within sequence of base flags to return

	int									// NOTE: returns original flags (EntryID flags in bits 3..0) immediately prior to the set set/reset
		SetBaseFlags(UINT32 EntryID,		// identifies sequence containing loci flags to be set
				UINT32 Loci,			// offset within sequence of flags to set
	            int SetFlgs=0,			// set flags in bits 0..3
				int ResetFlgs=0);		// reset flags in bits 0..3

	int									// NOTE: returns original flags (EntryID1 flags in bits 3..0 and EntryID2 in bits 7..4) immediately prior to the set set/reset
		SetBaseFlags(int EntryID1,		// identifies 1st sequence containing loci flags to be set
				UINT32 Loci1,			// offset within sequence of flags to set
				int EntryID2,			// identifies 2nd sequence containing loci flags to be set
				UINT32 Loci2,			// offset within sequence of flags to set
	            int SetFlgs=0,			// set flags in bits 0..3
				int ResetFlgs=0);			// reset flags in bits 0..3

	int GetBase(int EntryID,			// identifies sequence containing base to be returned
				UINT32 Loci);			// offset within sequence of base to return

	UINT32					// returned sequence length (may be shorter than that requested) or 0 if errors
		GetSeq(int EntryID,UINT32 Loci,etSeqBase *pRetSeq,UINT32 Len);	// get sequence for entry starting at Loci and of length Len
	
	UINT32					// returned sequence length (may be shorter than that requested) or 0 if errors
		GetColorspaceSeq(int EntryID,UINT32 Loci,etSeqBase *pRetSeq,UINT32 Len); // if indexed in colorspace then returns sequence as colorspace

	int	QSortSeq(INT64 SeqLen,		// total concatenated sequence length
						etSeqBase *pSeq,	// pts to start of concatenated sequences
						int SfxElSize,		// suffix element size (will be either 4 or 8)
						void *pArray);		// allocated to hold suffix elements
	void SetMaxQSortThreads(int MaxThreads);			// sets maximum number of threads to use in multithreaded qsorts

	int						// returns the previously utilised MaxBaseCmpLen
		SetMaxBaseCmpLen(int MaxBaseCmpLen);		// sets maximum number of bases which need to be compared for equality in multithreaded qsorts, will be clamped to be in range 10..(5*cMaxReadLen)

	UINT32			// number of exactly matching KMers up to specified MaxToMatch inclusive, or 0 if no match, MaxToMatch+1 if more than MaxToMatch Kmers
		NumExactKMers(UINT32 MaxToMatch,		// if 0 then no limit on matches, otherwise match up to this number of instances and if more then return MaxToMatch+1 
					int KMerLen,				// KMer length to exactly match over
					etSeqBase *pKMerSeq);		// pts to KMer sequence


	static int CompareExtdMatches( const void *arg1, const void *arg2);


};
