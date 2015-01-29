#pragma once

const int cBSFRdsVersion = 6;			// latest file header version which can be handled
const int cBSFRdsVersionBack= 5;		// backward compatible to this version

const unsigned int cMaxInFileSpecs = 20;	// allow user to specify upto this many input file specs

const int cMaxWorkerThreads = 64;		// allow for at most 64 worker threads (will be clamped to the max available)
const int cMaxReadsPerBlock = 256;		// max number of reads allocated for processing per thread as a block (could increase but may end up with 1 thread doing more than fair share of workload)

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

const int cDfltAllowedSubs = 5;			// default max number of aligner induced substitutions
const int cMaxAllowedSubs = 15;			// allow at most this many aligner induced substitutions
const int cMinCoreLen = 5;				// absolute minimum core length supported
const int cMaxCoreLen = 100;			// absolute maximum core length supported

const int cSNPBkgndRateWindow = 51;		// window size within which the substitution background rate is determined for SNP detection processing
										// should be odd size so that the centric base has same number of flanking bases over which background is calc'd

const int cSNPCentfFlankLen = 3;			// SNP centroid flank length
const int cSNPCentroidLen = ((cSNPCentfFlankLen * 2) + 1); // SNP centroid k-mer length
const int cSNPCentroidEls = ( 4 * 4 * 4 * 4 * 4 * 4 * 4);  // number of SNP centroid elements - could develop a macro but not worth the effort???

const int cDfltSensCoreIters  = 5000;	// default sensitivity core explore depth 
const int cMoreSensCoreIters  = 25000;	// more sensitivity core explore depth
const int cUltraSensCoreIters = 50000;	// ultra sensitivity core explore depth
const int cMinSensCoreIters   = 1000;	// min sensitivity core explore depth

const int cPriorityExacts = 10;			// when attempting to prioritorise exactly matching reads to higher confidence sequence then allow for this many more multiloci exact (no subs) hits

const int cDfltMaxNs = 1;				// default max number of indeterminate 'N's in the read or aligned to subsequence
const int cMaxNs = 5;					// allow user to specify at most this number of 'N's in the read or aligned to subsequence


const double cMinSNPproportion = 0.25;	// only check for SNP if at least this proportion of bases are not ref bases
const double cDfltQValueSNP = 0.05;     // default - used in Benjamini–Hochberg 
const int cDfltMinSNPreads = 5;         // default - must be at least this number of reads covering the loci before exploring as a SNP
const double cMaxBkgdNoiseThres = 0.20;	// background noise must be <= this value before exploring as being a SNP
const int cMinSNPreads = 2;             // user can specify down to this many reads covering the loci covering the loci before exploring as a SNP (assumes must be processing error free reads!)
const int cMaxSNPreads = 100;           // user can specify a max of at least this many reads covering the loci covering the loci before exploring as a SNP
const double cMinSeqErrRate = 0.01;		// sets a floor on minimum sequencing error rate per base - is used in binominal calculation
const int cAllocLociPValues = 100000;   // allocate for putative SNP loci in this many increments

const int cAllocLineBuffSize = 8196000; // when writing to results file then can buffer up to this many chars

const int cDfltMaxMultiHits = 5;		// default is to process at most this number of per read multihits
const int cMaxMultiHits = 500;			// user can specify at most this many multihits
const int cMaxAllHits = 100000;			// but if reporting all multihit loci then limit is increased to this value

const int cAllocMultihits = 250000;	 // alloc/realloc for multihit loci in this many instance increments
const int cDfltReadLen = 128;			 // assume reads of of this length with descriptors - not critical as actual read lengths are processed
const size_t cReadsHitReAlloc = 10000000; // realloc allocation to hold this many read instances

const int cPairMinLen =	25;				// apparent paired reads sequences must be of at least this length
const int cPairMaxLen = 5000;			// apparent paired reads sequences must be no longer than this length
const int cDfltPairMinLen = 100;		// default apparent paired reads sequences must be of at least this length
const int cDfltPairMaxLen = 300;		// default apparent paired reads sequences must be no longer than this length

// following constants are very empirical, relate to determination of multimatch read loci association, and need to be fine tuned..
const UINT16 cUniqueClustFlg = 0x08000;	// used to flag read as being clustered to reads which are uniquely aligned
const int cClustMultiOverLap = 10;		// to be clustered, reads must overlap by at least this number of bp 
const int cClustUniqueScore = 5;		// unique reads given this score if in current cluster window
const int cClustMultiScore = 1;			// multimatch read loci given this score if in current cluster window
const int cClustScaleFact = 10;			// scores are scaled down by this factor
const UINT32 cMHminScore = 50;			// any putative multimatch alignment score must be at least this to be accepted as the alignment for that read

const int cReadHitBuffLen = 0xfffff;		// sets per thread buffer size for holding string output read hit records
const int cDataBuffAlloc = 0x07ffffff;		// alloc to hold reads in this byte sized increments
const int cRdsBuffAlloc =   0x07fffff;		// alloc to hold preprocessed reads (for stats) in this byte sized allocation

const int cMaxDescrLen = 128;				// allow for descriptors of upto this length

const unsigned int cMinSeqLen = 15;			// sequences must be at least this length otherwise user is warned and sequence sloughed

const int cDfltRelSiteStartOfs = -4;		// default relative octamer start site offset, default gives 4bp 5' of read start
const int cMaxSitePrefOfs = 100;			// allow octamer starts to be at most this away from read start sites 

#pragma pack(1)

typedef enum eHLalign {
	eHLunique = 0,				// read aligned to a unique hit loci
	eHLrandom = 1,				// read multialigned and loci was randomly choosen
	eHLclustunique = 2,			// read multialigned and loci near other uniquely aligned reads was choosen
	eHLclustany = 3			    // read multialigned and loci near other multialigned reads was choosen
} teHLalign;


// a specific read may align to multiple loci and each of these aligned to loci may be fragmented due to InDels
// 
typedef struct TAG_sReadAlignLoci {
	UINT32 ReadHitIdx;			// read aligning
	UINT32 NxtReadAlign;		// next alignment for same read
	UINT16  FlagHL:2;			// how read hit loci was choosen (see enum eHLalign) 
	UINT16  FlagMH:1;			// set if this hit is part of a multihit read, reset if a unique read hit
	UINT16  FlagMHA:1;			// set if this hit is provisionally accepted as the alignment to report
	UINT16  FlagSegs:1;			// set if this hit has been segmented into 2 hit loci - e.g. contains an InDel or splice junctions
	UINT16  FlagCA:1;			// set if this hit has been determined as being correctly aligned (original loci known as read was generated by kangas)
	UINT16  FlagIA:1;			// set if this hit has been determined as being incorrectly aligned (original loci known as read was generated by kangas)
	UINT16  FlagTR:1;			// set if this hit has been end trimmed
	UINT16 ReadOfs;				// alignment starts at this read sequence offset
	UINT8 Strand;				// alignment to this strand
	UINT32 ChromID;				// suffix entry (chromosome) matched on
	UINT32 MatchLoci;			// original match loci
	UINT16 MatchLen;			// original match length
	UINT8 Mismatches;			// original number of mismatches
	UINT8 TrimLeft;				// left flank trimming removes this many bases
    UINT8 TrimRight;			// right flank trimming removes this many bases
	UINT8 TrimMismatches;		// after trimming there are this many mismatches
} tsReadAlignLoci;

// reads may have a loci to which they align, if multiple loci then multiple instances of tsHitLoci for each alignment
typedef struct TAG_sReadHitLoci {
	tsHitLoci Hit;				// hits can contain 2 aligned segments so as to accommodate InDel or splice junctions
	UINT8  FlagHL:2;			// how read hit loci was choosen (see enum eHLalign) 
	UINT8  FlagMH:1;			// set if this hit is part of a multihit read, reset if a unique read hit
	UINT8  FlagMHA:1;			// set if this hit is provisionally accepted as the alignment to report
	UINT8  FlagSegs:1;			// set if this hit has been segmented into 2 hit loci - e.g. contains an InDel or splice junctions
	UINT8  FlagCA:1;			// set if this hit has been determined as being correctly aligned (original loci known as read was generated by kangas)
	UINT8  FlagIA:1;			// set if this hit has been determined as being incorrectly aligned (original loci known as read was generated by kangas)
	UINT8  FlagTR:1;			// set if this hit has been end trimmed
} tsReadHitLoci;

// read hit loci, descriptor, sequence, quality
typedef struct TAG_sReadHit {
	UINT32 ReadHitIdx;			// current read hit index + 1 for this read
	UINT32 ReadID;				// read identifier from the preprocessed read (tsRawReadV5)
	UINT32 PairReadID;			// this read's paired raw read identifier (0 if not a paired read) bit32 reset if 5' fwd read, bit32 set if 3' read of pair
	UINT32 NumReads;			// number of source reads merged into this read
	UINT8 DescrLen;				// descriptor length
	UINT16 ReadLen;				// read length of sequence packed into Read following the descriptor
	INT16  LowHitInstances;		// number of hit target loci instances at LowMMCnt mismatches
	INT8   LowMMCnt;			// lowest number of mismatches for this read thus far
	INT8   NxtLowMMCnt;			// next lowest number of mismatches for this read thus far
	INT16  NumHits;				// number of target hits for this read, currently limited to be at most 1
	tsReadHitLoci HitLoci;		// this is currently the best hit for this read
	UINT8  SiteIdx;				// read has been characterised into this index into m_OctSitePrefs[]
	UINT8  Read[1];				// packed read descriptor, sequence and quality values for this read
	} tsReadHit;

// Hamming specific structures
const int cMaxHammingChroms = 200;	// can handle at most this many chromosomes with precalculated hammings
typedef struct TAG_sHamChrom {
	UINT32 ChromID;					// uniquely identifies this chromosome
	UINT8  szChrom[cMaxDatasetSpeciesChrom];	// chrom name
	UINT32 NumEls;					// number of subsequences with hammings on this chrom
	UINT8 Dists[1];					// array, in ascending loci order, of hamming distances
} tsHamChrom;

typedef struct TAG_sHamHdr {
	UINT8 Magic[4];		        // magic chars 'bham' to identify this file as a biosequence file containing hamming edit distances
	UINT32 Version;				// structure version 
	INT32 Len;					// file length, also allocation size required when loading hammings into memory
	UINT16 NumChroms;		    // number of chromosomes with Hammings
	UINT32 ChromOfs[cMaxHammingChroms];	// offsets to each chromosomes respective tsHamChrom
} tsHamHdr;

#pragma pack()
