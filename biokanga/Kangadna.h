// Copyright 2013 CSIRO  ( http://www.csiro.au/ ) m_bRawDedupe
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

const UINT32 cPPCRdsFileVersion = 2;			// current file header version
const UINT32 cPPCRdsFileVersionBack= 2;			// backward compatible to this version

typedef UINT32 tSeqWrd4;					// used to hold upto 15 packed 2bit canonical bases plus 2 flg bits in bits 30 and 31 (32bits or 4 bytes total)
typedef UINT32 tSeqWrd1;					// used to hold single canonical base in bits 0..1 plus plus 2 flg bits in bits 6 and 7
typedef UINT32 tSeqID;						// sequence identifiers are 32bits
typedef UINT32 tVertID;						// vertice identifiers are 32bits
typedef UINT64 tEdgeID;						// edges identifiers are 64bits
typedef UINT32 tSubGraphID;					// disconnected subgraph identifiers are 32bits 

const UINT32 cMinSeqs2Assemb = 100;			// require at least this many sequences for preprocessing

const UINT32 cMaxRawSeqLen   = 0x03ffffff;	 // allow for raw seed SE sequences loaded from file of no more than 64Mbp (26bits required)
const UINT32 cMaxQScoreLen   = 0x03ffff;	 // if filtering on the basis of minimum Phred scores then allow for fastq sequence lengths up to this many bases and quality scores
const UINT32 cMaxContigLen   = 0x3fffffff;	 // allow for assembled contigs of no more than 1Gbp (30bits required)
const UINT32 cMaxSeqID       = 0xffffffff;	 // allow for sequence identifiers of upto ~4Billion (32bits required)

// Every packed sequence is immediately preceded by 3 x 32bit (tSeqWrd4's) header words
// Header words have their most significant bit 31 set whereas the packed sequence words have bit 31 reset so header words can easily be distingished from sequence words
// The first header word contains most significant 2 bits of a 32bit sequence identifier (in bits 24 and 25), a 8bit source file identifier (in bits 16..23), 16bit flags (in bits 0..15) 
// The second header word contains the low 30bits of the 32bit sequence identifier (in bits 0..29)
// The third header word contains the 30bit sequence length (in bits 0..29)

// Top 2 bits (bits 31,30) as used to flag the payload type in the following 30 bits of the current word
const tSeqWrd4 cSeqWrd4BOS = 0xA0000000;	// bit 31 set, bit 30 reset, bit 29 set, bit 0 reset, flags flags this as the all concatenated sequences BOS marker word 
const tSeqWrd4 cSeqWrd4EOS = 0xA0000001;	// bit 31 set, bit 30 reset, bit 29 set, bit 0 set, flags flags this as the all concatenated sequences EOS marker word 
const tSeqWrd4 cSeqWrd4MSWHdr = 0x80000000;	// bit 31 set, bit 30 reset, bit 29 reset flags this as the most significant header word
const tSeqWrd4 cSeqWrd4MSWMsk = 0xe0000000;	// mask suitable for determining if a cSeqWrd4MSWHdr
const tSeqWrd4 cSeqWrd4LSWHdr = 0xC0000000;	// bit 31 and bit 30 both set flags this as the next significant header word
const tSeqWrd4 cSeqWrd4PartSeq=0x40000000;	// bit 31 reset and bit 30 set flags low order 29..0 bits as containing 1..14 bases
const tSeqWrd4 cSeqWrd4Msk  = 0x3fffffff;	// payload mask for bits 29..0 - sequence identifier, sequence length, or 15 bases 

const tSeqWrd1 cSeqWrd1PartSeq=0x40000000;	// bit 7 reset and bit 6 set flags low order 1..0 bits as containing 1 bases
const tSeqWrd1 cSeqWrd1Msk  = 0x3f;			// payload mask for bits 5..0  

// sequences being processed for duplicates/filtered are marked with these flags (can have a maximum of 16 defined as must fit within UINT16
const UINT16 cFlgSeqPE        = 0x001;		// sequence is paired ended; if reset then sequence is single ended
const UINT16 cFlgSeqPE2       = 0x002;		// sequence is PE2; if reset then sequence is single ended or PE1 with partner flagged as being cFlgSeqPE2
const UINT16 cFlg5Prime		  = 0x004;		// the 5' end of this sequence is overlapped, or overlaps at least 1 other sequence
const UINT16 cFlg3Prime		  = 0x008;		// the 3' end of this sequence is overlapped, or overlaps at least 1 other sequence
const UINT16 cFlgNoProc		  = 0x010;		// sequence is marked as not requiring further processing
const UINT16 cFlgSeqUnique	  = 0x020;		// sequence is marked as being unique with no other exactly matching sequences; or if duplicates not filtered then cFlgSeqUnique, cFlgSeq1stDup and cFlgSeqNthDup will be reset
const UINT16 cFlgSeq1stDup	  = 0x040;		// sequence is marked as being the representative duplicate sequence with other indentically matching sequences marked as being cFlgSeqNthDup
const UINT16 cFlgSeqNthDup	  = 0x080;		// sequence is marked as being a duplicate of one or more another sequences with one representative duplicate sequence marked as cFlgSeq1stDup
const UINT16 cFlgSeqRemove	 = 0x0100;		// sequence is marked as being for removal, either because it is a duplicate and duplicates being removed, or because it is not fully overlapped
const UINT16 cFlgOvlLenMsk	 = 0xFE00;		// sequence overlays, or is overlaid by, at least one other sequence by this percentage of the sequence length 

// sequences being processed for de Novo assembly are marked with these flags (can have a maximum of 16 defined as must fit within UINT16
//const UINT16 cFlgSeqPE        = 0x001;		// sequence is paired ended; if reset then sequence is single ended
//const UINT16 cFlgSeqPE2       = 0x002;		// sequence is PE2; if reset then sequence is single ended or PE1 with partner flagged as being cFlgSeqPE2
const UINT16 cFlgAsmbSeed       = 0x004;		// sequence being used as seed to be extended by other sequences
const UINT16 cFlgAsmbExtn       = 0x008;		// sequence is being used to extend some other cFlgAsmbSeed sequence
const UINT16 cFlgAsmbCplt       = 0x010;		// sequence accepted and not available for further assembly
const UINT16 cFlgNonOverlap     = 0x020;		// sequence determined to have inconsistent sense relative to other sequences and not to be used when attempting to overlap with other sequences
const UINT16 cFlgContainRemove     = 0x040;		// sequence has been marked for removal, sequence is completely contained by some other SE sequence

// Suffix sparsity applies to the word boundaries at which suffixes will be indexed
// Word boundaries depend on the base packing which will be one of the user specified following; 
// 1 bases per byte (8bits)
// 15 bases per unsigned int (32bits)
typedef enum TAG_eSfxSparsity {
	eSSparsity1 = 1,		    // every base is indexed (1 base per byte) tSeqWrd1
	eSSparsity15 = 15		    // every 15th base is indexed (15 bases packed into 32bits) tSeqWrd4
	} etSfxSparsity;
const int cDfltSuffixSparsity = eSSparsity15;	// default suffix sparsity supported is at 15bp (32bit word) boundaries

const UINT32 cMaxSfxBlkEls = 4000000000;	// construct suffix block arrays with 4byte array elements if no more than this, otherwise use 5 byte array elements

const UINT64 cMaxConcatSeqLen = (UINT64)0x0ffffffffff; // arbitary limit to all concatenated read sequences lengths - 1Tbp should be enough!   

const int cMaxWorkerThreads = 128;			// limiting max number of threads to this many
const int cMaxDupInstances = 2500;			// maintain instance counts up this max number of instances

const int cMaxMultiSeqFlags = 4000;			// limit local copy of sequence flags to at most this many

const int cKDNAMaxInFileSpecs = 100;		// currently allowing up to this many input files each for PE1 (or single ended) and PE2, 150 total max plus potentially one seed contig file
										// could be increased but a major limitation to extending max input files is that packed sequence headers are constructed with a 8bit source file identifer limit 
const int cMinSEOvlp  = 20;				// user can specify down to this minimal SE overlap required to merge SEs
const int cMaxSEOvlp  = 300;			// user can specify up to this maximum minimal SE overlap required to merge SEs
const int cDfltInitSEOvlp = 150;		// default initial minimal SE overlap required to merge SEs
const int cDfltFinSEOvlp  = 25;			// default final minimal SE overlap required to merge SEs

const int cMinPEOvlp = 35;				// user can specify down to this minimal sum of end overlaps required to merge PEs
const int cMaxPEOvlp = 200;				// user can specify up to this maximum sum of end overlaps required to merge PEs
const int cDfltInitPEOvlp = 150;		// initial minimal PE sum of end overlaps required to merge PEs
const int cDfltFinPEOvlp  = 35;			// default final minimal PE sum of end overlaps required to merge PEs
const int cMinPE1PE2ToSEOvlp  = 16;		// user can specify down to this minimal overlap of PE1 onto PE2 required to merge as SE
const int cMaxPE1PE2ToSEOvlp  = 100;	// user can specify up to this maximum overlap of PE1 onto PE2 required to merge as SE
const int cDfltMinPE1PE2ToSEOvlp  = 20;	// minimal overlap of PE1 onto PE2 required to merge as SE

// if PEs are unable to overlap then this could be because of allelic issues
const int cMinPESeqLen2SE = 125;			// individual ends must be of at least this length and
const int cMinPETotSeqLen2SE = 600;			// if any individual or sum PE1 + PE2 lengths sum to >= cMinPETotSeqLen2SE then reclassify as being two separate SE sequences in final assembly phases

const int cMaxOvrlapSeqWrds = (1024 * 1024);		// process at most this many overlapping tSeqWrd4s
const int cMaxCmpOverlapBases = (cMaxOvrlapSeqWrds*15);	// handle sequences of this maximal unpacked length when comparing for overlaps
const int cMaxSortSfxLen = cMaxCmpOverlapBases;		// sort suffixes out to this maximal suffix (bases)
const int cMinSeedContigLen = 100;					// will accept seed contigs which are >= this minimum length
const int cMaxSeedContigLen = cMaxCmpOverlapBases;	// truncate seed contigs which are > this length to this length
const  int cMinDfltSeqLenToAssemb = 80;			// default is to require any sequence, after any trimming, loaded for assembly to be at least this length in bp 

const int cOVLnone = 0x01;					// Overlay classification flags - flags no overlay between Seq1 and Seq2


											
const int cMinReqPESepDist = 5;				// require PEs when mapped onto SE sequences to have a minimum separation (insert - read length) distance
const int cMaxReqPESepDist = 1500;			// require PEs when mapped onto SE sequences to have a maximum separation (insert - read length) distance

const int cMaxDiagAsciiSeqLen = 25000000;	// allocate to hold upto this many bases in ascii when reporting on contig sequences

const int cReadsThreadStackSize = 0x01fffff; // increase stack size to 2MB for working threads

// Suffix array indexes are sparse, every Nth base is suffix indexed  
// processing modes
typedef enum TAG_ePMode {
	ePMdefault,					// default processing mode is to load and assemble raw reads
	ePMfileLoadTypeSeqs,		// load sparse suffix reads from file and assemble
	ePMfileSaveTypeSeqs,		// save preprocessed sparse suffix reads to file
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

typedef enum TAG_eOvlFlankPhase {
	eOvlpSenseToSense = 0,				// it's sense overlap sense flank processing (probe cFlg5Prime and target cFlg3Prime set if overlap) 
	eOvlpAntiSenseToSense,				// it's antisense overlap sense flank processing (probe cFlg5Prime and target cFlg5Prime set if overlap)
	eOvlpSenseToAntiSense				// it's sense overlap antisense flank processing (probe cFlg3Prime and target cFlg3Prime set if overlap)
	} etOvlFlankPhase;

#pragma pack(1)

typedef struct TAG_sPregenLevenshteinsHdr {
	UINT8 Signiture[4];  // will be initialised to contain 'levd'
	UINT32 LevKMers;	// number of Kmers, matrix size will be LevKMers * LevKMers
	UINT32 KMerLen;		// generate levenshtein distances for all k-mers of this length (currently limited to 8bp)
	UINT32 MaxLev;		// will be only accepting overlapping k-mers if their Levenshteins is <= MaxLev
	INT32 PenaltyMatch;	// match penalty
	INT32 PenaltyMisMatch;	// override the default mismatch penalty
	INT32 PenaltyInsert;	// override the default insert penalty
	INT32 PenaltyDelete;	// override the default deletion penalty
	UINT32 m_AcceptLevWidth;	// each matrix row is this width in bytes, each byte contains bit flags packed 8 per byte 
	UINT32 m_AllocAcceptLevDist;  // allocation size for m_pAcceptLevDist
} tsPregenLevenshteinsHdr;

typedef struct TAG_sThreadPregenLevPars {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to class instance
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;				// levelshtein processing result
	UINT32 NumAccepted;		// number of KMers accepted as being within MaxLev
	UINT32 NumRejected;     // number of KMers not accepted as being within MaxLev
	UINT32 NumRows;		// number of KMers to be processd by this thread
	UINT32 StartKMer;	// start KMer inclusive
	UINT32 EndKMer;		// end KMer inclusive
	UINT32 LevKMers;	// number of Kmers, matrix size will be LevKMers * LevKMers
	UINT32 KMerLen;		// generate levenshtein distances for all k-mers of this length (currently limited to 8bp)
	UINT32 MaxLev;		// will be only accepting overlapping k-mers if their Levenshteins is <= MaxLev
	INT32 PenaltyMatch;	// match penalty
	INT32 PenaltyMisMatch;	// override the default mismatch penalty
	INT32 PenaltyInsert;	// override the default insert penalty
	INT32 PenaltyDelete;	// override the default deletion penalty
} tsThreadPregenLevPars;


typedef struct TAG_sProcReadsCtrl {
	teBSFrsltCodes Rslt;			// if not eBSFSuccess then errors whilst reading from file
	bool bPE1PE2AllReadsProc;		// set true when all reads for current SE or PE1 and PE2 files have been completely processed
	bool bProcPE;					// true if processing PE 
	UINT32 NumPE1ParsedReads;		// number of reads parsed from the current SE or PE1 fasta file
	UINT32 NumPE1DescrReads;		// number of reads with descriptors parsed from the current SE or PE1 fasta file
	UINT32 NumPE2ParsedReads;		// number of reads parsed from the current PE2 fasta file
	UINT32 NumPE2DescrReads;		// number of reads with descriptors parsed from the current PE2 fasta file

	UINT32 CurTotPE1ReadsParsed;	// current number of reads, either SE or PE, parsed by all threads
	UINT32 CurTotPE1ReadsAccepted;	// current number of reads, either SE or PE, accepted by all threads

	int MinPhredScore;				// minimum acceptable Phred score
	bool bIsPE1Fastq;				// true if filtering SE or PE1 on fastq quality scores
	bool bIsPE2Fastq;				// true if filtering PE2 on fastq quality scores
	int PE1FileID;					// PE1 reads are being processed fron this file
	int PE2FileID;					// PE2 reads are being processed fron this file
	char szPE1File[_MAX_PATH];		// name of SE or PE1 reads file
	CFasta *pFastaPE1;				// SE or PE1 reads are being processed from this fasta file
	char szPE2File[_MAX_PATH];		// name of PE2 reads file
	CFasta *pFastaPE2;				// if PE processing then PE2 reads are being processed from this fasta file
} tsProcReadsCtrl;

typedef struct TAG_sThreadFiltReadsPars {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to class instance
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	tsProcReadsCtrl *pProcReadsCtrl;    // processing reads from these files and using these thresholds

	int Rslt;				// processing result
	int MaxNs;				// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
	int MinPhredScore;		// filter out input sequences with mean Phred score lower than this threshold
	int Trim5;				// trim this number of 5' bases from input sequences (default is 0, range 0..20)
	int Trim3;				// trim this number of 3' bases from input sequences (default is 0, range 0..20)
	int MinSeqLen;		    // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 20..10000)
	int TrimSeqLen;			// trim sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...10000)
	int PE1QCSchema;		// SE/PE1 using this quality scoring schema - guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger
	int PE2QCSchema;		// PE2 using this quality scoring schema - guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger

	int SampleNth;			// process every Nth reads
	int Zreads;				// maximum number of reads to accept for processing from any file

	int PE1FileID;			// PE1 reads are being processed fron this file
	int PE2FileID;			// PE2 reads are being processed fron this file

	UINT32 PE1ReadLen;
	UINT8 PE1RawReadsBuff[cMaxRawSeqLen+1];
	UINT32 PE1QSLen;
	UINT8 szPE1QScoresBuff[cMaxQScoreLen+1];
	UINT32 PE2ReadLen;
	UINT8 *PE2RawReadsBuff[cMaxRawSeqLen+1];
	UINT32 PE2QSLen;
	UINT8 szPE2QScoresBuff[cMaxQScoreLen+1];

	int NumPE1ParsedReads;   // thread parsed this many PE1 reads
	int NumPE1AcceptedReads; // and this many were accepted after filtering
	int NumPE2ParsedReads;   // thread parsed this many PE2 reads
	int NumPE2AcceptedReads; // and this many were accepted after filtering

	INT64 AcceptedTotSeqLen; // sum of all sequence lengths, both PE1 and PE2 accepted by this thread
	int MaxAcceptedReadLen;	 // max length of any read, either PE1 or PE2, accepted by this thread
	int MinAcceptedReadLen;	 // min length of any read, either PE1 or PE2, accepted by this thread

	int NumPE1Underlen;		// number of PE1's not accepted as after any trimming they were under length
	int NumPE2Underlen;		// number of PE2's not accepted as after any trimming they were under length

	int NumPE1UnderPhredScore; // number of PE1's not accepted as under MinPhredScore threshold
	int NumPE2UnderPhredScore; // number of PE2's not accepted as under MinPhredScore threshold

	int NumPE1ExcessNs;		   // number of PE1's not accepte as more than this number of MaxNs indeterminates
	int NumPE2ExcessNs;		   // number of PE2's not accepte as more than this number of MaxNs indeterminates

	UINT32 NumContamVector;		// number of reads contained within contaminate vector
	UINT32 NumContamFlankTrim;  // number of reads which were flank trimmed because of overlap onto contaminate adaptor
} tsThreadFiltReadsPars;
 
typedef struct TAG_sSeqOverlay {
	tSeqID SeqID1;		// overlay is between this sequence
	int Seq1WrdBytes;	//    bases are packed into either 1 or 4 bytes corresponding respectively to eSSparsity1 or eSSparsity15  
	int Seq1Len;		//    sequence 1 unpacked length
	tSeqID SeqID2;		// and this sequence
	int Seq2Len;		//    sequence 2 unpacked length
	int Seq2WrdBytes;	//    bases are packed into either 1 or 4 bytes corresponding respectively to eSSparsity1 or eSSparsity15  
	UINT16 OvlFlgs;		// overlay classification flags
	int OverlayLen;		// overlay is of this length
	int NumSubs;		// with this many substitutions required
	int Seq1StartOfs;	// overlay starts at this offset on Seq1
	int Seq2StartOfs;	// overlay starts at this offset on Seq2
	int Seq1EndOfs;		// overlay finishes at this offset on Seq1
	int Seq2EndOfs;		// overlay finishes at this offset on Seq2
	bool bSeq1Packed;	// set true if Seq1 still in packed format
	bool bSeq2Packed;	// set true if Seq2 still in packed format
	UINT8 Seq1[cMaxCmpOverlapBases+1];	// holds unpacked (bSeq1Packed == false) or packed (bSeq1Packed == true) sequence 1	
	UINT8 Seq2[cMaxCmpOverlapBases+1];  // holds unpacked (bSeq2Packed == false) or packed (bSeq2Packed == true) sequence 1
} tsSeqOverlay;



typedef struct TAG_sSequences {		// contains sequences and sparse suffix for sequences
	UINT8 bSeqs2AssembDirty;		// set true whenever any changes to NumSeqs2Assemb, TotSeqs2Assemb, Seqs2AssembLen or pSeqs2Assemb have been modified
									// and new suffix and descriptors/starts should be generated  
	UINT8  bPESeqs;					// non-zero if sequences are from paired end reads
	UINT32 TotSeqsParsed;			// total number of sequences parsed
    UINT32 TotSeqsUnderLen;			// total number of sequences filtered out because underlength
	UINT32 TotSeqsExcessNs;			// total number of sequences filtered out because too many Ns
	UINT32 NumProcessed;			// number processed
	UINT32 NumDuplicates;			// of which this number are duplicates
	UINT32 NumOverlapping;			// number of sequences which overlapped other sequences
	UINT64 NumOverlapped;			// number of sequences determined as being overlapped
	UINT32 TotSeqs2Assemb;			// original total number of sequences to be assembled after length filtering
	UINT32 NumSeqs2Assemb;			// total number of sequences accepted into m_pSeqs2Assemb
	UINT32 NumPE1Seqs2Assemb;		// total number of sequences accepted into m_pSeqs2Assemb from PE1 (or single ended reads)
	UINT32 NumPE2Seqs2Assemb;		// total number of sequences accepted into m_pSeqs2Assemb from PE2
	UINT64 Seqs2AssembLen;			// the sequences in m_pSeqs2Assemb total to this number of bases
	UINT8  SfxSparsity;				// suffix sparsity - expected to be one of eSSparsity1 or eSSparsity15
	UINT32 SeqWrdBytes;				// bases are packed into either 4 bytes corresponding to eSSparsity15, or a single byte eSSparsity1 if sequences are unpacked  
	UINT32 SeqHdrLen;				// length (in tSeqWrds) of sequence headers; currently 3 (tSeqWrd4) for both eSSparsity1 and eSSparsity4
	UINT64 Seqs2AssembOfs;			// offset in m_pSeqs2Assemb to write next sequence
	UINT64 AllocMemSeqs2Assemb;		// memory (bytes) currently allocated for m_pSeqs2Assemb 
	UINT32 NumSeqStarts;			// number of sequence start offsets in pSeqStarts
	UINT64 AllocMemSeqStarts;		// memory (bytes) currently allocated to m_pSeqStarts 
	UINT32 NumSeqFlags;				// number of sequence flags used in pSeqFlags
	UINT64 AllocMemSeqFlags;		// memory (bytes) currently allocated to pSeqFlags
	double MeanSeqLen;				// mean length of all sequences
	UINT32 MinSeqLen;				// minimum sequence length
	UINT32 MaxSeqLen;				// maximum sequence length
	UINT32	SfxElSize;				// suffix element size - either 4, <= cMaxSfxBlkEls, or 5 if > cMaxSfxBlkEls 
	UINT64 NumSuffixEls;			// number of elements in suffix array
	UINT64 AllocMemSfx;				// allocated memory size for suffix array
	UINT64 OfsSeqs2Assemb;			// when loading from/to disk then holds the file offset at which the concatenated packed sequences start
	union {							// as a union to ensure that size of this structure remains constant independent of 32/64 bit compilation 
		struct {
			void *pTmpRevCplSeqs;       // used to hold sequences whilst revcpl PEs, either as tSeqWrd4s if packed or tSeqWrd1s if unpacked sequences
			void *pSeqs2Assemb;				// holds sequences concatenated with UINT32 lengths prepended to each sequence which are to be subsequently suffix indexed and assembled
			void *pSuffixArray;				// to hold suffix array for concatenated read/contig sequences, element sizes may be either 4 or 5 bytes as specified by m_SfxElSize
			UINT64 *pSeqStarts;				// holds array of sequence start tSeqWrd1 or tSeqWrd4 offsets (indexed by sequence identifiers)
			UINT16  *pSeqFlags;				// holds array of sequence flags (indexed by sequence identifiers)
			};
		UINT64 Padd4Ptrs[5];
		};
} tsSequences;

typedef struct TAG_sEstSeqs {
	bool bCalcMeanSeqLen;			// set true when ever the mean sequence length needs to be re-calculated
	UINT32 NumFiles;				// number of files from which these estimates were derived 
	UINT64 TotSeqLens;				// estimates of total sequence length  
	UINT32 TotNumSeqs;				// estimates of total number of sequences to be parsed
	UINT32 NumSeqsToProc;			// estimate of total number of sequences, after any sampling to be further processed
	UINT32 MeanSeqLen;				// estimates of mean sequence length
} tsEstSeqs;


typedef struct TAG_sMultiSeqFlags {
	tSeqID SeqID;				// sequence identifier (32 bits)
	UINT16 SetFlags;			// flags to set
	UINT16 ResetFlags;			// flags to be reset
	UINT16 CurFlags;			// after flag processing then contains the updated flags		
} tsMultiSeqFlags;

typedef struct TAG_sReadFile {
	UINT8 FileID;					// uniquely identifies this sequence source file
	UINT8 PEFileID;					// if paired end processing then it's partner source file identifier
	UINT32 SeqsParsed;				// number of sequences parsed
    UINT32 SeqsUnderLen;			// number of sequences filtered out because underlength
	UINT32 SeqsExcessNs;			// number of sequences filtered out because too many Ns
	UINT32 SeqsAccepted;			// number of sequences accepted
	UINT8 szFile[_MAX_PATH];		// source file name
	} tsReadFile;

typedef struct TAG_sSparseSfxFileHdr {
	UINT8 Magic[4];					// magic chars to identify this file as a preprocessed and packed read sequences file
	UINT32 Version;					// file structure version (cPPCRdsFileVersion)
	UINT64 FileSize;				// file should be of this size
	UINT32 NumRawFiles;				// number of raw reads files processed
	tsReadFile SrcFiles[1 + (cKDNAMaxInFileSpecs * 2)]; // original source files from which sequences were parsed, allows for PE1 and PE2 source file names plus an optional seed contig file
	tsSequences Sequences;			// to hold sequences for PE1 and PE2
} tsPPCRdsFileHdr;

const int cSizeofPPCRdsFileHdr = sizeof(tsPPCRdsFileHdr);


// because the base packing is assuming that only canonical bases ACGT are packed into 2bits then indeterminates 'N' are substituted with a random canonical base
// and the starting offsets of contiguous blocks 'N's of these substitions recorded so that after processing the 'N's can be restored
// currently this indeterminate processing is only utilised in scaffolding as it is enables progressive scaffolding with different mate pair libraries 
typedef struct TAG_sBlockNsLoci {
	tSeqID SeqID;				// sequence identifier (32 bits) containing this block of indeterminates
	UINT32 Ofs;					// offset (0..n) in the sequence at which the indeterminates start
	UINT32 NumNs;				// number of indeterminates in this block
	} tsBlockNsLoci;

#pragma pack()

class CScaffolder;
class CdeNovoAssemb;
class CArtefactReduce;
class CPacBio;

class CKangadna
{
friend class CScaffolder;				// CScaffolder requires access to all members of this class including private members - overhead of setter/getters would be too much of a performance hit
friend class CdeNovoAssemb;				// CdeNovoAssemb requires access to all members of this class including private members - overhead of setter/getters would be too much of a performance hit
friend class CArtefactReduce;			// CArtefactReduce requires access to all members of this class including private members - overhead of setter/getters would be too much of a performance hit
friend class CPacBio;					// CPacBio requires access to all members of this class including private members - overhead of setter/getters would be too much of a performance hit

	int m_PMode;				// processing mode (actually it's the requested output mode in preprocessing phase)

	CMTqsort m_MTqsort;			// multithreaded sorting

	CContaminants *m_pContaminants;	// instantiated if filtering read sequences for contaminant adaptors or primers

	int m_NumRawFiles;			// number of raw reads files processed
	int m_NumRdsFiles;			// number of preprocessed (kangar) reads files processed
	tsReadFile m_SrcFiles[1 + (cKDNAMaxInFileSpecs * 2)];		// to hold names of all source sequence files, allows for PE1 and PE2 to each have max allowed plus an optional seed contig file

	int m_hInFile;				// input file handle
	char m_szInFile[_MAX_PATH];	// input file

	int m_hOutFile;				// output file handle
	char m_szOutFile[_MAX_PATH];// output file

	int m_hOrphansFile;			// output orphans file handle
	char m_szOrphansFile[_MAX_PATH];// output orphans file

	int m_hInSeqTypesFile;	// input processed seq types file handle
	char m_szProcSeqTypesFile[_MAX_PATH];// input sparse suffix file

	int m_hOutSeqTypesFile;				// output sparse suffix file handle
	char m_szOutSeqTypesFile[_MAX_PATH];// output sparse suffix file	
		
	int m_hOutFastaSE;			// used for writing deduped and error reduced reads to multifasta SE
	int m_hOutFastaR1;			// used for writing deduped and error reduced reads to multifasta PE1
	int m_hOutFastaR2;			// used for writing deduped and error reduced reads to multifasta PE2
	char m_szCtgDescr[80];		// contig descriptor prefix


	bool m_bRawDedupe;				// removing all duplicate reads
	bool m_bDedupeIndependent;		// if paired end preprocessing then treat as if single ended when deduping

	bool m_bAssembStrand;			// true if read strand specific assembly
	int m_AssembMinOverlap;			// minimum flank overlap length
  	int m_AssembMeanInsert;			// if paired end processing - mean insert size (defaults to 0 for single ended assembly)
	int m_AssembMinInsert;			// if paired end processing - paired end minimum insert size
	int m_AssembMaxInsert;			// if paired end processing - paired end maximum insert size

	int m_MinReqMergeLen;				// only merge PE ends if would result in SE of at least this length

	int m_MinOverlapReductFact;		// factor/100 (50..99) by which to reduce required overlaps on each assembly pass

	int m_NReduceThresSteps;		// how many threshold lowering steps over which to reduce the initial overlap thresholds down to their minimums
	int m_InitialReqPEPrimOverlap;	// 1st pass, if primary probe is overlapping onto a PE then the initial primary overlap (onto PE1 or PE2) must be of at least this length
	int m_InitialReqPESecOverlap;	// 1st pass, if primary probe was overlapping onto a PE then the secondary probe overlap (onto PE1 or PE2) must be of at least this length
	int m_InitialReqPESumOverlap;	// 1st pass, if primary probe was overlapping onto a PE then the sum of the PE1 and PE2 overlap must be of at least this length
	int m_InitialReqSEPrimOverlap;	// 1st pass, if primary probe is overlapping onto a SE then the overlap must be of at least this length
	// initial minimum overlaps are reduced by m_MinOverlapReductFact in subsequent assembly passes with lower limits according to the following
	int m_MinReqPEPrimOverlap;		// if primary probe is overlapping onto a PE then the initial primary overlap (onto PE1 or PE2) must be of at least this length
	int m_MinReqPESecOverlap;		// if primary probe was overlapping onto a PE then the secondary probe overlap (onto PE1 or PE2) must be of at least this length
	int m_MinReqPESumOverlap;		// if primary probe was overlapping onto a PE then the sum of the PE1 and PE2 overlap must be of at least this length
	int m_MinReqSEPrimOverlap;		// if primary probe is overlapping onto a SE then the overlap must be of at least this length
	int m_MinPE2SEOvlp;				// if probe PE1 and probe PE2 being considered for merging then there must be an overlap of at least this many bases
	int m_MinReqPESepDist;			// PE start sites expected to be separated by at least this many bases
	int m_MaxReqPESepDist;			// PE start sites expected be separated by no more than this many bases



	int m_NumThreads;			// max number of processing threads to use 
	int m_ThreadsProcessing;	// number of threads dispatched to handle current processing task, decremented as each thread completes current task phase
	
	bool m_bAffinity;			// thread to core affinity

	char *m_pszLineBuff;			// allocd for buffering of output assembled contigs
	int m_LineBuffLen;				// current number of chars buffered in m_pszLineBuff
	int m_AllocLineBuff;			// m_pszLineBuff allocated to hold at most this number of chars

	UINT64 m_CurMaxMemWorkSetBytes;  // currently set max working set in bytes
	UINT32 m_WinPageSize;			// windows memory page size in bytes (0 if process not on Windows)
	UINT64 m_BaseWinMinMem;			// windows base min working set memory in bytes when process initially started
	UINT64 m_BaseWinMaxMem;			// windows base max working set memory in bytes when process initially started

	tsSequences m_Sequences;			//  5' (eSTypePE1) and 3' (eSTtypePE2) sequences

	tsEstSeqs m_SeqEsts;			// estimates of sequence lengths + total number of sequences for each sequence type

	UINT32 m_NumPartialSeqs2Assemb;	 // total number of partial sequences to assemble
	UINT64 m_LenPartialSeqs2Assemb;  // total length of partial sequences
	UINT64 m_PartialSeqs2AssembOfs;	 // offset into m_pPartialSeqs2Assemb at which next sequence can be copied
	UINT64 m_AllocdPartialSeqs2Assemb; // bytes memory allocated to m_pPartialSeqs2Assemb;
	void *m_pPartialSeqs2Assemb;		// to hold partially assembled sequences during each assembly phase


	UINT32 m_NumBlockNsLoci;		// number of blocks currently used
	UINT32 m_AllocdBlockNsLoci;		// currently this number of blocks have been allocated - may be on demand reallocated
	size_t m_AllocdBlockNsLociSize;  // memory allocation size
	tsBlockNsLoci *m_pBlockNsLoci;   // pts to blocks of indeterminate start loci + len

	int								// returned sequence length, will be limited to MaxSeqLen
		GetSeq(tSeqID SeqID,		// sequence identifier
			UINT8 *pRetSeq,			// where to copy unpacked sequence bases
			UINT32 MaxSeqLen = 0,	// return at most MaxSeqLen bases (0 if no limit)
			bool bNoEOSTerm = false);	// option to not terminate unpacked string with EOS

	int								// returned sequence length, will be limited to MaxSeqLen
		GetSeq(void *pStartSeqWrd,
				UINT8 *pRetSeq,				// where to copy unpacked sequence bases
				UINT32 MaxSeqLen = 0,		// return at most MaxSeqLen bases (0 if no limit)
				bool bNoEOSTerm = false);	// option to not terminate unpacked string with EOS
	
	int // returns sequence length for sequence with specified SeqID
			GetSeqLen(tSeqID SeqID); // sequence identifier

	teBSFrsltCodes	// returns min/mean/max lengths for PE1/PE2 and SE sequences
		GetSeqLenDist(UINT32 *pNumPEs = NULL,	// total number of paired ends (number of pairs!)
					  int *pPE1min = NULL,	// returned PE1 min length
					  int *pPE1mean = NULL,// returned PE1 mean length
					  int *pPE1max = NULL, // returned PE1 max length
					  int *pPE2min = NULL,	// returned PE2 min length
					  int *pPE2mean = NULL,// returned PE2 mean length
					  int *pPE2max = NULL, // returned PE2 max length
					  UINT32 *pNumSEs = NULL,	// total number of single ends
					  int *pSEmin = NULL,	// returned SE min length
					  int *pSEmean = NULL,	// returned SE mean length
					  int *pSEmax = NULL);	// returned SE max length


	teBSFrsltCodes ChunkedWrite(UINT64 WrtOfs,UINT8 *pData,UINT64 WrtLen);
	teBSFrsltCodes ChunkedWrite(int hFile,char *pszFile,UINT64 WrtOfs,UINT8 *pData,UINT64 WrtLen);
	teBSFrsltCodes ChunkedRead(int hFile,char *pszFile,UINT64 RdOfs,UINT8 *pData,UINT64 RdLen);
	teBSFrsltCodes ChunkedRead(UINT64 RdOfs,UINT8 *pData,UINT64 RdLen);

				


	void ResetThreadedIterReads(void); // must be called by master thread prior to worker threads calling ThreadedIterReads()
	void AcquireSerialise(void);
	void ReleaseSerialise(void);
	void AcquireSerialiseNxtProcRead(void);
	void ReleaseSerialiseNxtProcRead(void);
	void AcquireSerialiseSeqHdr(void);
	void AcquireSerialiseReadsCtrl(void);
	void ReleaseSerialiseReadsCtrl(void);
	void ReleaseSerialiseSeqHdr(void);
	void AcquireSerialiseSeqFlags(void);
	void ReleaseSerialiseSeqFlags(void);
	void AcquireLock(bool bExclusive);				// defaults as read only lock
	void ReleaseLock(bool bExclusive);

	CStopWatch m_StopWatch;

	bool m_bMutexesCreated;						// set true if mutexes and rwlocks created/initialised

#ifdef _WIN32
	HANDLE m_hMtxIterReads;
	HANDLE m_hMtxIterNxtProcRead;
	HANDLE m_hMtxMHReads;
	CRITICAL_SECTION m_hSCritSectSeqHdrs;
	CRITICAL_SECTION m_hSCritSectSeqFlags;
	SRWLOCK m_hRwLock;
	
#else
	pthread_mutex_t m_hMtxIterReads;
	pthread_mutex_t m_hMtxIterNxtProcRead;
	pthread_mutex_t m_hMtxMHReads;
	pthread_spinlock_t m_hSpinLockSeqHdrs;
	pthread_spinlock_t m_hSpinLockSeqFlags;
	pthread_rwlock_t m_hRwLock;
#endif

#ifdef WIN32
	alignas(4)	volatile unsigned int m_CASSeqFlags; // used with synchronous compare and swap (CAS) for serialising access to sequence header flags
	alignas(4)	volatile unsigned int m_CASNxtProcRead; // used with synchronous compare and swap (CAS) for serialising access to next read for loading from file
	alignas(4)	volatile unsigned int m_CASReadsCtrl; // used with synchronous compare and swap (CAS) for serialising access to shared read control when loading from file
#else
	__attribute__((aligned(4))) volatile unsigned int m_CASSeqFlags; // used with synchronous compare and swap (CAS) for serialising access to sequence header flags
	__attribute__((aligned(4))) volatile unsigned int m_CASNxtProcRead; // used with synchronous compare and swap (CAS) for serialising access to next read for loading from file
	__attribute__((aligned(4))) volatile unsigned int m_CASReadsCtrl; // used with synchronous compare and swap (CAS) for serialising access to shared read control when loading from file
#endif

	int CreateMutexes(void);
	void DeleteMutexes(void);

	static inline UINT64		// unpacks the 5 bytes ptd at by p5Bytes and returns as UINT64 
		Unpack5(UINT8 *p5Bytes)
		{
		// ensures that only 5 bytes are actually accessed, can't just access as a UINT64 and mask retain bottom 40 bits...
		return((UINT64)((UINT64)p5Bytes[4] << 32 | *(UINT32 *)p5Bytes));
		}

	static inline UINT8 *Pack5(UINT64 Value, UINT8 *p5Bytes)
		{
		*(UINT32 *)p5Bytes = (UINT32)Value;
		p5Bytes[4] = (UINT8)(Value >> 32);
		return(p5Bytes + 5);
		}

	teBSFrsltCodes AddSeq(UINT32 SrcFileID,		// 8bits identifies source file from which this sequence was processed
						UINT32 SeqFlags,    	// 16bits as bit flags
						UINT32 SeqLen,			// 30bit sequence length
						UINT8 *pSeq);			// ptr to sequence

	teBSFrsltCodes AddPESeqs(UINT32 SrcFileID,			// 8bits identifies source file from which this sequence was processed
						UINT32 PE1SeqFlags,    	// PE1 16bits as bit flags
						UINT32 PE1SeqLen,		// PE1 30bit sequence length
						UINT8 *pPE1Seq,			// PE1 ptr to sequence
						UINT32 PE2SeqFlags,    	// PE2 16bits as bit flags
						UINT32 PE2SeqLen,		// PE2 30bit sequence length
						UINT8 *pPE2Seq);			// PE2 ptr to sequence

	int
		GenPackedSeqWrds(int SeqLen,	// sequence length ptd at by pSeq
				etSeqBase *pSeq,		// pts to sequence requiring packing
				int MaxSeqWrds,			// truncate packed sequence if requiring more than this number of packed words, 0 if no limit
				void *pDstSeq);			// write generated packed sequence words starting at this location

	int									// returned number of tSeqWrd4s used to hold packed sequence
		GenPackedSeqWrds(int SeqLen,	// sequence length ptd at by pSeq
				etSeqBase *pSeq,		// pts to sequence requiring packing
				int MaxSeqWrds,			// truncate packed sequence if requiring more than this number of packed words, 0 if no limit
				tSeqWrd4 *pDstSeq);		// write generated packed sequence words starting at this location

	tSeqID	GetPackedSeqID(tSeqWrd4 *pSeqWrd);	// ptr into sequence for which sequence identifer is to be returned
	int GetPackedSeqBaseLen(tSeqWrd4 *pSeqWrd);
	int GetPackedSeqFlags(tSeqWrd4 *pSeqWrd);
	int GetPackedSeq5BaseLen(tSeqWrd4 *pSeqWrd);	   // returns number of bases 5' left of tSeqWrd32 ptd at by pSeqWrd

	int										// returns number of bases 3' right including tSeqWrd4 ptd at by pSeqWrd
		GetPackedSeq3BaseLen(void *pSeqWrd);    
	int GetPackedSeq3BaseLen(tSeqWrd4 *pSeqWrd);    // returns number of bases 3' right including tSeqWrd4 ptd at by pSeqWrd

	int									// number of tSeqWrds returned in pDstPackedSeq (excludes optional EOS)
		GetPackedSeq(int MaxSeqWrds,			// limit to this many tSeqWrds (0 if no limit)
			 tSeqWrd4 *pSrcPackedSeq,		// get from this packed sequence
			 tSeqWrd4 *pDstPackedSeq,   	// copy into this packed sequence
			 bool bNoEOSTerm = false);		// option to not terminate packed substring with EOS


	// NOTE: GetSeqWrdSubSeq assumes that the source sequence contains a preceding header containing the actual source length 
	int										// number of packed bases returned in pDstPackedSeq (excludes optional EOS)
		GetSeqWrdSubSeq(int SubOfs,			// substring starting at this base offset within pSrcPackedSeq
			 int SubLen,					// substring required which contains at most this many bases (0 for maximal length)
			 tSeqWrd4 *pSrcPackedSeq,		// extract substring from this packed string
			 tSeqWrd4 *pDstPackedSeq,		// extract into this packed substring
			 bool bNoEOSTerm = false);		// option to not terminate packed substring with EOS

	// NOTE: GetNHSeqWrdSubSeq assumes that the source sequence does not contain a preceding header with actual source length and thus expects SubOfs + SubLen to be <= actual source length !!! 
	int										// number of packed bases returned in pDstPackedSeq (excludes optional EOS)
		GetNHSeqWrdSubSeq(int SubOfs,		// substring starting at this base offset within pSrcPackedSeq (SubOfs + SubLen expected to be <= pSrcPackedSeq length)
			 int SubLen,					// substring required which contains at most this many bases
			 tSeqWrd4 *pSrcPackedSeq,		// extract substring from this packed string
			 tSeqWrd4 *pDstPackedSeq,		// extract into this packed substring, may be same as pSrcPackedSeq if inplace substring required
			 bool bNoEOSTerm = false);		// option to not terminate packed substring with EOS
	// NOTE: NHSeqTrim assumes that the source sequence does not contain a preceding header with actual source length and thus expects Trim5 + Trim3 to be <= actual source length !!! 
	int						// returned trimmed sequence length
		NHSeqTrim( int Trim5,		// trim 5' by this many bases
			 int Trim3,		// trim 3' by this many bases
			 int SeqLen,	// number of bases in sequence before sequence trimmed
			 tSeqWrd4 *pSeq, // inplace trimming of this sequence
			 bool bNoEOSTerm = false);				// option to not terminate packed substring with EOS


	bool SetMaxMemWorkSetSize(size_t Bytes);

	
	void	// inplace reverse complement of unpacked bases
		RevCplSeq(unsigned int SeqLen, // sequence to reverse complement is of this length
							etSeqBase *pSeq);	// pts to start of sequence to be reverse complemented

	int										// returned number of bases overlapped between SeqA and SeqB
		GetOverlapAB(int MinOverlap,		// seqA and seqB must overlap by at least this minimum overlap (expected to be at least 16bp)
			int *pSeqALeft,					// returned number of bases 5' before seq A is overlapped onto seq B (negative if seq B is instead overlapping onto seq A)
			int SeqALen,					// SeqA is this length
			tSeqWrd4 *pSeqA,				// look for overlap between this sequence A and sequence B
			int SeqBLen,					// Seq B is this length
			tSeqWrd4 *pSeqB,				// sequence B
			bool bNoSeqBOvrlapA = false,	// if true and no seq 1 onto seq B overlap then do not look for seq B overlap onto seq A
			int MaxSubsPerK = 0,			// allow upto this many substitutions per 1000 bp of overlap (expected to be in the range 0..20)
			int MaxEnd12Subs = 0,			// allow at the initial 12bp of the 5' or 3' of overlap to have this many base mismatches in addition to the overall allowed MaxSubs (expected to be in range 0..6) 
			int *pNumSubs = NULL);			// total number of substitutions actually required

	static int								// 0 Seq1 == Seq2, -1 if Seq1 < Seq2, +1 if Seq1 > Seq2
		CmpPackedSeqs(tSeqWrd4 *pSeq1,			// Seq1 (probe) packed sequence
				tSeqWrd4 *pSeq2,				// Seq2 (target) packed sequence
				int MaxCmpLen);				// compare over at most this many bases

	int										// returned length of exactly matching initial (starting at ofs 0) substring between probe and target 
		GetExactMatchLen(tSeqWrd4 *pSeq1,	// Seq1 (probe) packed sequence
				tSeqWrd4 *pSeq2,			// Seq2 (target) packed sequence
				int MaxLen=0);				// if non-zero then immediately return as soon as OverlapLen is at least this length - no need to exhustively find the maximal overlap

	bool									// true if Seq1 matches Seq2 for MatchLen
	IsMatched(int ReqMatchLen,				// required match length
			int Seq1Ofs,					// base offset in Seq1 at which to to start match
			int Seq1Len,					// Seq1 is this length
			tSeqWrd4 *pSeq1,				// look for match between this sequence 1 and sequence 2
			int Seq2Ofs,					// base offset in Seq2  at which to to start match
			int Seq2Len,					// Seq 2 is this length
			tSeqWrd4 *pSeq2,				// sequence 2
			int MaxSubsPerK,				// allow upto this many substitutions per 1000 bp of overlap (expected to be in the range 0..50)
			int MaxEnd12Subs,				// allow at the initial 12bp of the 5' or 3' of match to have this many base mismatches in addition to the overall allowed MaxSubs (expected to be in range 0..6) 
			int *pNumSubs);					// total number of substitutions actually required

	etSeqBase								// returns a single base at Ofs 0..N relative to the SeqWrd ptd at by pSeq
			GetBase(int Ofs,				// offset of base to return
				tSeqWrd4 *pSeq);				// Seq1 (probe) packed sequence containing bases

	int			// inplace shift left by one base the sequence ptd at by pSeq - note that any sequence header is not modified!
		ShfLeftPackedSeq(int CurLen,		// Seq1 currently contains this many bases, will contain CurLen-1 after shift left
						tSeqWrd4 *pSeq);		// Seq1 (probe) packed sequence to shift left


	static int SfxSortSeqWrd4Func(const void *arg1, const void *arg2);
	static int Sfx5SortSeqWrd4Func(const void *arg1, const void *arg2);


	UINT32 m_LevDistKMerLen;                    // Levenshtein distances have been precalculated for sequences of this K-mer length
	UINT32 m_MaxLev;							// and if <= this distance have been bit set as accepted in m_pAcceptLevDist
	UINT32 m_LevKMers;							// number of K-mers of len
	int m_PenaltyMatch;							// override the default match penalty (0)
	int m_PenaltyMisMatch;						// override the default mismatch penalty (1)
	int m_PenaltyInsert;						// override the default insert penalty (1)
	int m_PenaltyDelete;						// override the default deletion penalty (1)
	UINT32 m_AcceptLevWidth;					// each matrix row is this width in bytes, each byte contains bit flags packed 8 per byte 
	UINT32 m_AllocAcceptLevDist;            // allocation size for m_pAcceptLevDist
	UINT8 *m_pAcceptLevDist;				// allocated to hold Levenshtein accepted as bit flags packed 8 per byte 


	tSeqWrd4 *									// returned ptr to 1st packed sequence word
			SfxIdxToFirstSeqWrd(UINT64 SfxWrdIdx, // index + 1
							   int *pWrdOfs = NULL);	  // returned sequence word is at this relative word offset to SfxWrdIdx; e.g if 1 then SfxWrdIdx referenced the second word in the sequence 
							
	tSeqID									// returned sequence identifer
		SfxIdx2SeqID(UINT64 SfxWrdIdx);		// (offset + 1) into pTypeSeqs  
			
	UINT64			// index+1 pSfxArray of first exactly matching probe or 0 if no match				
		LocateFirstExact(int ElSize,					// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
				  void *pProbe,					// pts to probe sequence
				  int ProbeLen,					// probe length (in bases, not tSeqWrds) to exactly match over
				  void *pTarg,					// target sequence
				  UINT8 *pSfxArray,				// target sequence suffix array
				  UINT64 SfxLo,					// low index in pSfxArray
				  UINT64 SfxHi);				// high index in pSfxArray

	UINT64			// index+1 in pSfxArray of last exactly matching probe or 0 if no match					
		LocateLastExact(int ElSize,					// sizeof elements in pSfxArray - currently will be either 4 or 5 bytes
				  void *pProbe,					// pts to probe sequence
				  int ProbeLen,					// probe length (in bases, not tSeqWrds) to exactly match over
				  void *pTarg,					// target sequence
				  UINT8 *pSfxArray,				// target sequence suffix array
				  UINT64 SfxLo,					// low index in pSfxArray
				  UINT64 SfxHi,					// high index in pSfxArray
				  UINT32 Limit);				// if non-zero then need only iterate towards last exactly matching this for this Limit iterations

	tSeqWrd4 *								// returned ptr to 1st word of actual packed sequence (use this as pCurSeq on next invocation)
		IterSeqHeaders(tSeqWrd4 *pCurSeq,	// iterate to next sequence following (NULL to return 1st sequence)
			tSeqID *pSeqID = NULL,			// returned 32bit sequence identifier
			UINT32 *pSrcFileID = NULL,		// returned 8 bit source file identifier
			UINT32 *pFlgs = NULL,			// returned 16 bit sequence flags
			UINT32 *pSeqLen = NULL,		    // returned 30 bit sequence length
			bool bSerialise = false);		// true if access to headers are to be serialised, false if access is not serialised

	tSeqWrd4 *							// returned ptr to 1st word of actual packed sequence
		GetSeqHeader(tSeqWrd4 *pSeqWrd,	// pts to a SeqWrd within the sequence
			tSeqID *pSeqID = NULL,		// returned 32 bit sequence identifier
			UINT32 *pSrcFileID = NULL,	// returned 8 bit source file identifier
			UINT32 *pFlgs = NULL,		// returned 16 bit sequence flags
			UINT32 *pSeqLen = NULL, 	// returned 30 bit sequence length
			bool bSerialise = false);	// set true if access to headers are required to be serialised

	tSeqWrd4 *								// returned ptr to 1st word of actual packed sequence
		GetSeqHeader(tSeqID SeqID,		// 32 bit sequence identifier
			UINT32 *pSrcFileID = NULL,	// returned 8 bit source file identifier
			UINT32 *pFlgs = NULL,		// returned 16 bit sequence flags
			UINT32 *pSeqLen = NULL,		// returned 30 bit sequence length
			bool bSerialise = false);	// set true if access to headers are required to be serialised

	tSeqWrd4 *									// returned ptr to next word to write packed sequence 
		SetSeqHeader(tSeqWrd4 *pSeqHdrWrd,		// pts to where to write the sequence header
			tSeqID SeqID,				// sequence identifier (32 bits)
			UINT32 SrcFileID,			// identifies source file from which this sequence was processed (8 bits)
			UINT32 Flgs,				// sequence flags (16 bits)
			UINT32 SeqLen,				// sequence length (30 bits)
			int *pNumSeqWrds = NULL,	// used to return the number of SeqWrds used
			bool bAddingHdr = true,     // set true if adding header, false if updating an existing header
			bool bSerialise = false);	// set true if access to headers are required to be serialised

	int									// returned updated flags, or if < 0 then error 
		UpdateSeqHeaderFlags(tSeqID SeqID,				// sequence identifier (32 bits)
			UINT32 SetFlags,			// flags to set
			UINT32 ResetFlags,			// flags to be reset
			bool bSerialise = false);	// set true if access to headers are required to be serialised

	int									// returned updated flags, or if < 0 then error 
		UpdateSeqHeaderFlags(tSeqWrd4 *pSeqWrd,				// pts to a SeqWrd within the sequence
			UINT32 SetFlags,			// flags to set
			UINT32 ResetFlags,			// flags to be reset
			bool bSerialise = false);	// set true if access to headers are required to be serialised

		int					// < 0 if errors
		UpdateAllSeqHeaderFlags(UINT32 SetFlags,	// flags to set
					UINT32 ResetFlags,				// flags to be reset
					bool bSerialise);				// set true if access to headers are required to be serialised

	int
		GetSeqFlags(tSeqID StartSeqID,			// starting from this sequence identifier inclusive
			tSeqID EndSeqID,					// and ending at this sequence identifier inclusive
			tsMultiSeqFlags *pMultiSeqFlags,	// holds sequence identifiers plus used to return flags
			bool bSerialise = false);			// set true if access to flags are required to be serialised

	int
		GetSeqFlags(UINT32 NumSeqs,						// number of sequences 
			tsMultiSeqFlags *pMultiSeqFlags,	// holds sequence identifiers plus used to return flags
			bool bSerialise = false);			// set true if access to flags are required to be serialised

	int
		GetSeqFlags(tSeqID SeqID,			// sequence identifier (32 bits)
			bool bSerialise = false);	// set true if access to headers are required to be serialised

	int
		UpdateSeqFlags(UINT32 NumSeqs,								// number of sequences requiring flag updates
			tsMultiSeqFlags *pMultiSeqFlags,			// holds sequence identifiers plus flag update operations
			bool bSerialise =false);							// set true if access to flags are required to be serialised

	int									// returned updated flags, or if < 0 then error 
		UpdateSeqFlags(tSeqID SeqID,				// sequence identifier (32 bits)
			UINT32 SetFlags,			// flags to set
			UINT32 ResetFlags,			// flags to be reset
			bool bSerialise = false);	// set true if access to headers are required to be serialised


	tSeqID m_StartProcSeqID;			// start sequence processing from this sequence identifer
	tSeqID m_FinalProcSeqID;			// finish sequence processing at this sequence identifer
	tSeqID m_NextProcSeqID;				// next sequence identifier for processing
	int m_NumProcSeqIDs;				// each thread to process at most this many sequences as a block

	INT64									// accumulated partial sequence lengths	
		SavePartialSeqs(int PE1Len,		// length of partially assembled sequence ptd at by pPE1
				tSeqWrd4 *pPE1,			// partially assembled PE1
				int PE2Len,				// (optional) length of partially assembled sequence ptd at by pPE2
				tSeqWrd4 *pPE2);		// (optional) partially assembled PE2

	int	GetSeqProcRange(tSeqID *pStartingSeqID,	// where to return starting sequence identifier (inclusive)
								tSeqID *pEndingSeqID,	// where to return ending sequence identifier (inclusive)
								int MaxReq = 0);		// at most this many identifiers limited to no more than m_NumProcSeqIDs  (0 if MaxReq is m_NumProcSeqIDs)

	int
		PackedRevCplAllIncPEs(void);		 // Will firstly exchange all PE1 and PE2 sequences with their headers followed by RevCpl all sequences

	int
		PackedRevCplAll(bool bRevOnly = false,			// if true then reverse only, false to reverse complement
						bool bPE2 = false);				// if true then reverse complement PE2 only
	int
		PackedRevCpl(tSeqWrd4 *pSeqToRev,
			bool bRevOnly = false,			// if true then reverse only, false to reverse complement
			UINT32 MaxSeqWrds = 0);			// reverse for at most this many SeqWrds (0 for no limit)

	int										// returned number of bases complemented
		PackedCpl(tSeqWrd4 *pSeqToCpl,
			UINT32 MaxSeqWrds = 0);			// complement for at most this many SeqWrds (0 for no limit)

	int										// returned number of packed sequence words copied
		CopyPackedSeq(void *pSrcSeq,		// copy from this start tSeqWrd
			  void *pDstSeq,				// copy to this start tSeqWrd
			  int MaxLen,					// copy at most this many tSeqWrds
			  bool bTermEOS);				// if true then terminate copied sequence with EOS tSeqWrd

	int											// number of bases in sequence pDstSeqWrd after merging or 0 if unable to merge because flags cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt already set for sequence SrcSeqID
		AtomicSeqMerge(int RelOverlap,			// if >= 0 then pSrcSeq starts at this pDstSeq[RelOverlap]; if < 0 then pDstSeq starts at this pSrcSeq[RelOverlap]; RelOverlap is in bases
					int SrcSeqLen,				// number of bases in pSrcSeqWrd sequence
				   tSeqWrd4 *pSrcSeqWrd,		// merge this sequence with pDstSeqWrd
					int DstSeqLen,				// number of bases currently in pDstSeqWrd sequence
				   tSeqWrd4 *pDstSeqWrd,		// merge this sequence with pSrcSeqWrd; merged sequence replaces this sequence
				   tSeqWrd4 *pTmpSeqWrd,		// temp sequence space for use as may be needed whilst merging
				   tSeqID  SeqID);			    // mandatory sequence identifier to test/set flags cFlgAsmbExtn | cFlgAsmbCplt

	int											// total number of bases in both pPE1DstSeqWrd and pPE2DstSeqWrd after merging or 0 if unable to merge because flags cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt already set for sequence Seq1ID or Seq2ID
			AtomicSeqMerge(int PE1RelOverlap,	// if >= 0 then pPE1SrcSeqWrd sequence starts at this pPE1DstSeqWrd[PE1RelOverlap]; if < 0 then pPE1DstSeqWrd starts at this pPE1SrcSeqWrd[PE1RelOverlap]; PE1RelOverlap is in bases
					int PE1SrcSeqLen,			// number of bases in pPE1SrcSeqWrd sequence
				    tSeqWrd4 *pPE1SrcSeqWrd,	// merge this sequence with pPE1DstSeqWrd

					int PE2RelOverlap,			// if >= 0 then pPE2SrcSeqWrd starts at this pPE2DstSeqWrd[PE2RelOverlap]; if < 0 then pPE2DstSeqWrd starts at this pPE2SrcSeqWrd[PE2RelOverlap]; PE2RelOverlap is in bases
					int PE2SrcSeqLen,			// number of bases in pPE2SrcSeqWrd sequence
				    tSeqWrd4 *pPE2SrcSeqWrd,	// merge this sequence with pPE2DstSeqWrd				   

					int PE1DstSeqLen,			// number of bases currently in pPE1DstSeqWrd sequence
					int *pPE1DstSeqLen,			// returned length after merging
				   tSeqWrd4 *pPE1DstSeqWrd,		// merge this sequence with pPE1SrcSeqWrd; merged sequence replaces this sequence
				   tSeqWrd4 *pPE1TmpSeqWrd,		// temp sequence space for use as may be needed whilst merging
					
					int PE2DstSeqLen,			// number of bases currently in pPE2DstSeqWrd sequence
					int *pPE2DstSeqLen,			// returned length after merging
				   tSeqWrd4 *pPE2DstSeqWrd,		// merge this sequence with pPE2SrcSeqWrd; merged sequence replaces this sequence
				   tSeqWrd4 *pPE2TmpSeqWrd,		// temp sequence space for use as may be needed whilst merging

				   tSeqID Seq1ID,				// mandatory sequence identifier to test/set flags cFlgAsmbExtn | cFlgAsmbCplt
				   tSeqID Seq2ID = 0);			// optional sequence identifier to test/set flags cFlgAsmbExtn | cFlgAsmbCplt

	int											// number of bases in sequence pDstSeqWrd after merging or 0 if unable to merge because flags cFlgAsmbSeed | cFlgAsmbExtn | cFlgAsmbCplt already set for sequence Seq1ID or Seq2ID
		AtomicSeqMerge(int PE1RelOverlap,		// if >= 0 then PE1 sequence starts at this pDstSeq[PE1RelOverlap]; if < 0 then pDstSeq starts at this PE1 sequence[PE1RelOverlap]; PE1RelOverlap is in bases
					int PE1SrcSeqLen,			// number of bases in pPE1SrcSeqWrd sequence
				    tSeqWrd4 *pPE1SrcSeqWrd,	// merge this sequence with pDstSeqWrd
					int PE2RelOverlap,			// if >= 0 then pPESrcSeqWrd starts at this pDstSeq[PE2RelOverlap]; if < 0 then pDstSeq starts at this PE2 sequence[PE2RelOverlap]; PE2RelOverlap is in bases
					int PE2SrcSeqLen,			// number of bases in pPE2SrcSeqWrd sequence
				    tSeqWrd4 *pPE2SrcSeqWrd,	// merge this sequence with pDstSeqWrd				   
				    int DstSeqLen,				// number of bases currently in pDstSeqWrd sequence
				   tSeqWrd4 *pDstSeqWrd,		// merge this sequence with pPE1SrcSeqWrd and pPE2SrcSeqWrd; merged sequence replaces this sequence
				   tSeqWrd4 *pTmpSeqWrd,		// temp sequence space for use as may be needed whilst merging
				   tSeqID Seq1ID,				// mandatory sequence identifier to test/set flags cFlgAsmbExtn | cFlgAsmbCplt
					tSeqID Seq2ID = 0);			// optional sequence identifier to test/set cFlgAsmbExtn | cFlgAsmbCplt


				// sequence merge; merges pSrcSeq and pDstSeq, merged sequence replaces pDstSeq 
	int											// number of bases in sequence pDstSeqWrd after merging
		SeqMerge(int RelOverlap,				// if >= 0 then pSrcSeq starts at this pDstSeq[RelOverlap]; if < 0 then pDstSeq starts at this pSrcSeq[RelOverlap]; RelOverlap is in bases
					int SrcSeqLen,				// number of bases in pSrcSeqWrd sequence
				   tSeqWrd4 *pSrcSeqWrd,		// merge this sequence with pDstSeqWrd
					int DstSeqLen,				// number of bases currently in pDstSeqWrd sequence
				   tSeqWrd4 *pDstSeqWrd,		// merge this sequence with pSrcSeqWrd; merged sequence replaces this sequence
				   tSeqWrd4 *pTmpSeqWrd,		// temp sequence space for use as may be needed whilst merging
					char *pszDiagText = NULL);			// optional diagnostics text


				// merge of a source PE1 and PE2 with a SE dest sequence
	int												// number of bases in sequence pDstSeqWrd after merging
		SeqMergePE12ToSE(int PE1RelOverlap,				// if >= 0 then PE1 sequence starts at this pDstSeq[PE1RelOverlap]; if < 0 then pDstSeq starts at this PE1 sequence[PE1RelOverlap]; PE1RelOverlap is in bases
					int PE1SrcSeqLen,			// number of bases in pPE1SrcSeqWrd sequence
				   tSeqWrd4 *pPE1SrcSeqWrd,		// merge this sequence with pDstSeqWrd
				
					int PE2RelOverlap,			// if >= 0 then pPESrcSeqWrd starts at this pDstSeq[PE2RelOverlap]; if < 0 then pDstSeq starts at this PE2 sequence[PE2RelOverlap]; PE2RelOverlap is in bases
					int PE2SrcSeqLen,			// number of bases in pPE2SrcSeqWrd sequence
				   tSeqWrd4 *pPE2SrcSeqWrd,		// merge this sequence with pDstSeqWrd				   
				   
				   int DstSeqLen,				// number of bases currently in pDstSeqWrd sequence
				   tSeqWrd4 *pDstSeqWrd,		// merge this sequence with pPE1SrcSeqWrd and pPE2SrcSeqWrd; merged sequence replaces this sequence
				   tSeqWrd4 *pTmpSeqWrd);		// temp sequence space for use as may be needed whilst merging

	char *
		AsciifySequence(void *pSeqWrds,			// sequence words to asciify
				UINT32 MaxSeqLen = 0,			// asciify at most cMaxDiagAsciiSeqLen bases
				bool bRevCpl = false);			// true if sequence to be reverse complemented before asciifying
	
	char *
		AsciifySequence(tSeqID SeqID,			// sequence identifier
				UINT32 MaxSeqLen = 0,			// asciify at most MaxSeqLen bases (limited to at most cMaxDiagAsciiSeqLen)
				bool bRevCpl = false);			// true if sequence to be reverse complemented before asciifying

	char *
		AsciifySequences(tSeqID SeqID1,			// sequence identifier
				tSeqID SeqID2,					// sequence identifier
				UINT32 MaxSeqLen,				// asciify at most MaxSeqLen bases
				char **ppszSeq1,				// returned seq1
				char **ppszSeq2);				// returned seq2

	void	DumpSeqs2Assemble(char *pszHeading,		// use this header text
								UINT32 Max2Dump = 100);		// dump at most this many sequences, 0 if no limit

	teBSFrsltCodes DumpHeader(char *pszTypeSeqFile);	// opens, dumps file header, and then closes pszTypeSeqFile

public:
	CKangadna(void);
	~CKangadna(void);

	bool SyncAllThreadsCompleted(bool bResetNextProcSeqID = false);	// returns true when all threads have completed current processing phase, option is to also reset m_NextProcSeqID

	int Reset(bool bSync = true);
	void ResetTypeSeqs(void);
	void SetCtgDescr(char *pszCtgDescr);	// set contig descriptor prefix
	teBSFrsltCodes SetNumThreads(int maxThreads,bool bAffinity=false);

	void SetPMode(int PMode = 0);							// set processing mode

	void SetHugePages(bool bHugePages=false);				// use huge pages if kernel supported 
	void SetDedupePE(bool bDedupeIndependent=false);		// dedupe policy on paired ends

	teBSFrsltCodes SetSfxSparsity(etSfxSparsity SfxSparsity);		// set suffix sparsity


	// Imortant:
    // Normally Levenshtein distance uses penalty of 0 for matches, and 1 for mismatches, inserts and deletions so that
    // the distance represents the number of edits between two sequences to make both equal
    // Because of the extreme InDel biases with PacBio long reads different penalties can be specified to the Levenshtein generating functions
static const int cLevenshteinMaxSeqLen = 2000;   // can calc levenshtein distances between sequences of upto this length
static const int cLevenshteinDefault = 1;         // levenshtein distance default penalty for mismatches, insertions and deletions, matches are 0 
static const int cLevenshteinPacBioMismatch = 2;  // PacBio has very low rates of base substitutions, most errors are InDels, so mismatches have penalty applied

	int
       LoadPregenLevenshteinsFwd(char *pszLevenshteinsFile, // load pregenerated Levenshteins from this file
				UINT32 KMerLen,		// expecting levenshtein distances for all k-mers to be of this length (must be in range 4..8)
				UINT32 MaxLev,		// expecting that will be accepting overlapping k-mers if their Levenshteins is <= MaxLev 
				int PenaltyMatch = 0,	// override the default match penalty
				int PenaltyMisMatch = cLevenshteinDefault,// override the default mismatch penalty
				int PenaltyInsert = cLevenshteinDefault,	// override the default insert penalty
				int PenaltyDelete = cLevenshteinDefault);	// override the default deletion penalty

	int
       SavePregenLevenshteinsFwd(char *pszLevenshteinsFile); // save pregenerated Levenshteins to this file

	UINT32										// total number of k-mer instances which have Levenshteins is <= MaxLev  
		PregenLevenshteinsFwd(UINT32 KMerLen,		// generate levenshtein distances for all k-mers of this length (must be in range 4..8)
						 UINT32 MaxLev,			// will be only accepting overlapping k-mers if their Levenshteins is <= MaxLev 
  					    int PenaltyMatch = 0,	// override the default match penalty
						int PenaltyMisMatch = cLevenshteinDefault,	// override the default mismatch penalty
						int PenaltyInsert = cLevenshteinDefault,	// override the default insert penalty
						int PenaltyDelete = cLevenshteinDefault);	// override the default deletion penalty

	static int										// returned levenshtein distance between RefSeq and QuerySeq, -1 if errors
		GetLevenshteinDistFwd(UINT32 KMerLen,		// number (1..16) of packed bases in both RefSeq and QuerySeq, bits 31/30 containing 5' base 
							UINT32 RefSeq,		// calculate levenshtein distance between packed sequence (2bits/base) in RefSeq, bits 31/30 containing 5' base 
						   UINT32 QuerySeq,     // and packed sequence (2bits/base) in QuerySeq, bits 31/30 containing 5' base
						   int PenaltyMatch = 0,	// override the default match penalty
						   int PenaltyMisMatch = cLevenshteinDefault,	// override the default mismatch penalty
						   int PenaltyInsert = cLevenshteinDefault,	// override the default insert penalty
						   int PenaltyDelete = cLevenshteinDefault);	// override the default deletion penalty


	static int										// returned levenshtein distance between RefSeq and QuerySeq, -1 if errors
		GetLevenshteinDistFwd(UINT32 KMerLen,		// number (1..1000) of bases in both pRefSeq and pQuerySeq 
							etSeqBase *pRefSeq,	// calculate levenshtein distance between this sequence (pRefSeq[0] contains 5' base)
						   etSeqBase *pQuerySeq, // and sequence in pQuerySeq
						   int PenaltyMatch = 0,	// override the default match penalty
						   int PenaltyMisMatch = cLevenshteinDefault,	// override the default mismatch penalty
						   int PenaltyInsert = cLevenshteinDefault,	// override the default insert penalty
						   int PenaltyDelete = cLevenshteinDefault,	// override the default deletion penalty
						   int MaxExpDist = 0);     // set to max expected distance to early terminate if observed distance would be greater; set 0 if no expected maximum

	static int			// returned levenshtein distance between RefSeq and QuerySeq, -1 if errors, assumes pRefSeq, pRefSeq point to last 3' bases 
		GetLevenshteinDistRev(UINT32 KMerLen,	// number (1..cLevenshteinMaxSeqLen) of bases in both pRefSeq and pQuerySeq 
							etSeqBase *pRefSeq,	// calculate levenshtein distance between this sequence (pRefSeq[0] contains 3' base)
						   etSeqBase *pQuerySeq, // and sequence in pQuerySeq (pQuerySeq[0] contains 3' base)
							int PenaltyMatch = 0,	// override the default match penalty (0)
						   int PenaltyMisMatch = cLevenshteinDefault,	// override the default mismatch penalty
						   int PenaltyInsert = cLevenshteinDefault,	// override the default insert penalty
						   int PenaltyDelete = cLevenshteinDefault,	// override the default deletion penalty
						   int MaxExpDist = 0);     // set to max expected distance to early terminate if observed distance would be greater; set 0 if no expected maximum

	static int												// returned affine gap levenshtein distance between RefSeq and QuerySeq, -1 if errors
		GetAGLevenshteinDistFwd(UINT32 KMerLen,	// number (1..cLevenshteinMaxSeqLen) of bases in both pRefSeq and pQuerySeq 
							etSeqBase *pRefSeq,	// calculate levenshtein distance between this sequence (pRefSeq[0] contains 5' base)
						   etSeqBase *pQuerySeq, // and sequence in pQuerySeq
							int PenaltyMatch = 0,	// override the default match penalty (0)
							int PenaltyMisMatch = 10,// override the default mismatch penalty (10)
							int PenaltyInDelOpen = 5,	// override the default indel gap open (5 to open)
							int PenaltyInDelExtd = 1);	// override the default indel gap extension (1 per base extension)

	bool										// returns true if QuerySeq within MaxDist from RefSeq
		IsLevenshteinAcceptedFwd(UINT32 KMerLen,	// number (1..16) of packed bases in RefSeq and QuerySeq, bits 31/30 containing 5' base
							UINT32 RefSeq,		// calculate levenshtein distance between packed sequence (2bits/base) in RefSeq, bits 31/30 containing 5' base 
						   UINT32 QuerySeq,     // and packed sequence (2bits/base) in QuerySeq, bits 31/30 containing 5' base
							UINT32 MaxLev,		// will be only accepting overlapping k-mers if their Levenshteins is <= MaxDist
						   int PenaltyMatch = 0,	// override the default match penalty
						   int PenaltyMisMatch = cLevenshteinDefault,	// override the default mismatch penalty
						   int PenaltyInsert = cLevenshteinDefault,	// override the default insert penalty
						   int PenaltyDelete = cLevenshteinDefault);	// override the default deletion penalty

	int ProcGenLevsFwd(tsThreadPregenLevPars *pPars); // threaded processing of PregenLevenshteins

	
	int GetSeqWrdBytes(void);			// returns size of tSeqWd 

		UINT32									// total number of reads accepted for processing into next phase
			GetNumReads(UINT32 *pNumPE1Reads=NULL,	// returned total number of single ended or 5' paired end read sequences accepted parsed
			UINT32 *pNumPE2Reads=NULL,				// returned total number of 3' paired end read sequences accepted parsed
			UINT32 *pTotPE1Seqs=NULL,				// number of PE1 sequences remaining at end of each phase completion
			UINT32 *pTotPE2Seqs=NULL,				// number of PE2 sequences remaining at end of each phase completion
			UINT32 *pTotSeqsParsed=NULL,			// returned total number of sequences parsed for 3' and 5' combined
			UINT32 *pTotSeqsUnderLen=NULL,			// returned total number of under length sequences parsed for 3' and 5' combined
			UINT32 *pTotSeqsExcessNs=NULL,			// returned total number of under length sequences parsed for 3' and 5' combined
			UINT32 *pMeanSeqLen = NULL,				// returned mean (rounded down) length of all sequences
			UINT32 *pMinSeqLen = NULL,				// returned minimum sequence length
			UINT32 *pMaxSeqLen = NULL);				// returned maximum sequence length			

	teBSFrsltCodes			// load previously saved concatenated and packed sequences from file 
		LoadPackedSeqsFromFile(char *pszTypeSeqFile);	// loading is from this file

	teBSFrsltCodes			// save concatenated and packed sequences to file 
		SavePackedSeqsToFile(char *pszTypeSeqFile);	// save to this file

	UINT32				// returned number of sequences in m_Sequences.pSeqs2Assemb
			ValidateSeqs2AssembStarts(UINT32 *pNumPEs = NULL,
							  UINT32 *pNumSEs = NULL);
	int ValidatePartialSeqsStarts(int *pPartialNumPEs = NULL,
							  int *pPartialNumSEs = NULL);


	int PackedRevCplPE(tSeqID PE1SeqID);	 // Will firstly exchange all PE1 and PE2 sequences with their headers followed by RevCpl all sequences

	int PackedRevCpl(tSeqID SeqID);			 // Will inplace reverse complement the specified sequence


	// accumulate estimates of read lengths, total number of sequences for each file type allowing for efficent subsequent memory allocations
	// memory estimates must be generated prior to calling LoadReads()
	teBSFrsltCodes 
			EstMemReq(char *pszInFile);      // accumulate estimates from this file

	UINT64			// returned estimate of memory size (bytes) required to hold packed sequences plus headers 
		EstMemReq(	int Trim5,				// will be trimming this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// will be trimming trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int TrimSeqLen,			// will be trimming sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...10000)
					int SampleNth,			// will only be processing every Nth reads
					int Zreads);				// will only accept this number of reads for processing from any file

	UINT32 GetEstSeqsToProc(void);			// returned estimate of number of sequences to process 
	
	tsEstSeqs * 
			GetEstMemReq(int Trim5,				// will be trimming this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// will be trimming trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int TrimSeqLen,			// will be trimming sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...10000)
					int SampleNth,			// will only be processing every Nth reads
					int Zreads);			// will only accept this number of reads for processing from any file
							
	static bool				// false if mean Phred score is below minimum threshold or any individual base Phred is < 10, true if Phred accepted
		MeetsMinPhredScoreThres(int QSSchema,	// quality scoring schema - guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger,
					int MinMeanPhredScore,		// minimum allowed mean (over all read bases) Phred score
					char *pszPhredScores);		// read ascii Phred scores

	// LoadReads
	// Load reads from fasta, fastq or csfasta formated raw reads file
	teBSFrsltCodes LoadReads(int MaxNs,		// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
					int MinPhredScore,		// filter out input sequences with mean Phred score lower than this threshold
					int Trim5,				// trim this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int MinSeqLen,		    // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)
					int TrimSeqLen,			// trim sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...10000)
					int SampleNth,			// process every Nth reads
					int Zreads,				// maximum number of reads to accept for processing from any file
					char *pszPE1File,		// file containing reads (kangar or raw fasta/fastq)
					char *pszPE2File);		// if paired end processing then PE2 3' file containing reads

	teBSFrsltCodes
	LoadReadsThreaded(int MaxNs,				// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
					int MinPhredScore,		// filter out input sequences with mean Phred score lower than this threshold
					int Trim5,				// trim this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int MinSeqLen,		    // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 20..10000)
					int TrimSeqLen,			// trim sequences to be no longer than this length (default is 0 for no length trimming, MinSeqLen...10000)
					int SampleNth,			// process every Nth reads
					int Zreads,				// maximum number of reads to accept for processing from any file
					char *pszPE1File,		// file containing reads (kangar or raw fasta/fastq)
					char *pszPE2File,		// if paired end processing then PE2 3' file containing reads
					int MaxNumThreads = 0);		// use at most this any threads, 0 if no limit

	teBSFrsltCodes										// returns < eBSFSuccess if any errors, eBSFSuccess if all reads already processed, 1 if SE or 2 if PE reads being returned for processing
					GetNxtProcRead(tsThreadFiltReadsPars *pPars);		// uniquely identifies calling thread, used to manage buffers for quality scores and sequences

	teBSFrsltCodes	LoadLongReads(char *pszLongReadsFile); // file holding long reads (multifasta or fastq)


	teBSFrsltCodes	// load high confidence seed contigs used when de Novo assembling 
					LoadSeedContigs(char *pszContigsFile, // file containing fasta seed contigs
					int TrimEnds = 0,				// trim input sequences, both 5' and 3' ends by this many bases
					int MinInputSeqLen = cMinDfltSeqLenToAssemb);	// seed contigs, after any trimming, must be of at least this length to be accepted


	teBSFrsltCodes	// load high confidence seed PE1 and PE2 used when de Novo assembling
			LoadSeedPEs(char *pszPE1File,		  // high confidence seed PE1 sequences file
					   char *pszPE2File,		  // high confidence seed PE2 sequences file
					   int OrientatePE = 0, 	  // PE end orientations 0: sense/antisense, 1: sense/sense, 2: antisense/sense, 3: antisense/antisense 
						int TrimEnds = 0,		  // trim input sequences, both 5' and 3' ends by this many bases
					  int MinInputSeqLen = cMinDfltSeqLenToAssemb);	// seed contigs, after any trimming, must be of at least this length to be accepted


	teBSFrsltCodes 
			AllocLoadBlock(char *pszInFile, // loading from this file
					 UINT64 FileOfs,				// load from this file offset
					 UINT64  AllocBlockSize,		// load, and allocate for, this block size from disk
					 void **ppLoadedBlock,		// returned ptr to allocated memory
					 UINT64 *pAllocBlockSize);	// size of allocated memory


	teBSFrsltCodes GenSeqStarts(bool bGenFlags=false,	// optionally generate flags array
								bool bSerialise = false);	// set true if access to headers are required to be serialised

	teBSFrsltCodes GenRdsSfx(int FirstNSeqWrds = 0, 	// max number of SeqWrds (0 to index all), starting from 1st, to index in each read sequence
		                     int ExcludeLastNSeqWrds = 0,	// exclude last N SeqWrds in each read sequence from indexing, 0 to index FirstNSeqWrds
							 bool bExclPE = false);	    // true to exclude sequences marked as being PE from being indexed

	teBSFrsltCodes AllocSeqs2AssembMem(UINT64 ReqAllocSize);	// alloc/realloc to at least ReqAllocSize (bytes)

	teBSFrsltCodes AllocBlockNsLoci(UINT32 ReqAllocBlocks);		// alloc/realloc to at least ReqAllocSize (tsBlockNsLoci)

	
	int GenContigSeqs(void);	// process overlapped reads and generate contig sequences


	int
		RemoveMarkedSeqs(UINT32 RemovalFlags,		// if any of these flags set then remove this sequence
								UINT32 AllRequiredFlags = 0,        // if containing any flags then remove this sequence unless all of these flags are set
								UINT32 AllOptionalFlags = 0,        // if containing any flags then remove this sequence unless at least one of these flags is set
								bool bUpdateHdrFlags = false);		// if true, and if pSeqFlags != NULL, then replace sequence header flags with those from pSeqFlags


						 

	int
		SetAssemblyPars(bool bStrand,			// true if read strand specific assembly
				  int MinOverlap,				// minimum flank overlap length
  				  int MeanInsert,				// if paired end processing - mean insert size (defaults to 0 for single ended assembly)
				  int MinInsert,				// if paired end processing - paired end minimum insert size
				  int MaxInsert);				// if paired end processing - paired end maximum insert size

	
	int	FreeSfx(void);
	int FreeSeqStarts(bool bFreeFlags = true);	// optionally also free flags array

	teBSFrsltCodes
		SaveAssembSeqs(char *pszFastaFile,		// save assembled sequences - will be saved as pszFastaFile suffxed with 'SE','PE1' and 'PE2'
						int PassID = 0,			// if non-zero then append to file names as the pass indicator
						int LineBreakSeqs = 0);	// default is to line break sequences every 75bp, if this is > 10 then line break every LineBreakSeqs bases

	teBSFrsltCodes
	SaveAssembAsSESeqs(char *pszFastaFile,	// save assembled sequences, PE sequences will be output as single concatenated sequence with 10 Ns separator - will be saved as pszFastaFile suffxed with '.fasta'
						int LineBreakSeqs = 0);	// default is to line break sequences every 75bp, if this is > 10 then line break every LineBreakSeqs bases


	// LoadReads() passes over each read (SE) or pair of reads (PE) to ProcReadsThread for filtering; ProcReadsThread will call AddSeq with reads which pass the filter
	int ProcReadsThread(tsThreadFiltReadsPars *pPars);

	teBSFrsltCodes SaveAsFasta(char *pszFastaFile);		// save as multifasta the deduped and error reduced input sequences
	int Test(int LocateTest = 0,		// set to max number of reads to be located (0 if no test required, -1 for all reads)
					 int SubSeqTest = 0,		// set to max number of subsequence iterations (0 if no test required, -1 for all reads)
					 int SampleNth = 1);		// locate test every Nth read
				
};



