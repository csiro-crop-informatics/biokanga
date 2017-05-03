/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once

const int cBSFRdsVersion = 6;			// latest file header version which can be handled
const int cBSFRdsVersionBack= 5;		// backward compatible to this version

const unsigned int cMaxInFileSpecs = 100;	// allow user to specify upto this many input file specs

const int cMaxWorkerThreads = 128;			// limiting max number of threads to this many
const int cMaxReadsPerBlock = 4096;		// max number of reads allocated for processing per thread as a block (could increase but may end up with 1 thread doing more than fair share of workload)

const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

const int cDfltAllowedSubs = 10;		// default max number of aligner induced substitutions per 100bp
const int cMaxAllowedSubs = 15;			// allow at most this many aligner induced substitutions per 100bp
const int cMaxTotAllowedSubs = 63;		// irrespective of length, limit total number of subs to 63
const int cMinCoreLen = 8;				// absolute minimum core length supported
const int cMaxCoreLen = 100;			// absolute maximum core length supported

const int cSNPBkgndRateWindow = 51;		// window size within which the substitution background rate is determined for SNP detection processing
										// should be odd size so that the centric base has same number of flanking bases over which background is calc'd

const int cSNPCentfFlankLen = 3;			// SNP centroid flank length
const int cSNPCentroidLen = ((cSNPCentfFlankLen * 2) + 1); // SNP centroid k-mer length
const int cSNPCentroidEls = ( 4 * 4 * 4 * 4 * 4 * 4 * 4);  // number of SNP centroid elements - could develop a macro but not worth the effort???
const int cDfltMaxDiSNPSep = 300;       // DiSNPs (TriSNPs) are two (three) individual SNPs separated by at most this many bp 

const int cMinMarkerLen = 25;			// minimum allowed marker length
const int cMaxMarkerLen = 500;			// maximum allowed marker length

const int cDfltSensCoreIters  = 5000;	// default sensitivity core explore depth 
const int cMoreSensCoreIters  = 10000;	// more sensitivity core explore depth
const int cUltraSensCoreIters = 20000;	// ultra sensitivity core explore depth
const int cMinSensCoreIters   = 2500;	// min sensitivity core explore depth

const int cPriorityExacts = 10;			// when attempting to prioritorise exactly matching reads to higher confidence sequence then allow for this many more multiloci exact (no subs) hits

const int cDfltMaxNs = 1;				// default max number of indeterminate 'N's in the read or aligned to subsequence
const int cMaxNs = 5;					// allow user to specify at most this number of 'N's in the read or aligned to subsequence


const double cDfltQValueSNP = 0.05;     // default - used in Benjamini–Hochberg 
const int cDfltMinSNPreads = 5;         // default - must be at least this number of reads covering the loci before exploring as a SNP
const double cMaxBkgdNoiseThres = 0.20;	// background noise must be <= this value before exploring as being a SNP
const int cMinSNPreads = 1;             // user can specify down to this many reads covering the loci covering the loci before exploring as a SNP (assumes must be processing error free reads!)
const int cMaxSNPreads = 100;           // user can specify a max of at least this many reads covering the loci covering the loci before exploring as a SNP
const double cMinSeqErrRate = 0.01;		// sets a floor on minimum sequencing error rate per base - is used in binominal calculation
const double cDfltMinMarkerSNPProp = (1.0/3.0);	// polymorphic bases within marker sequences must be at no more than this proportion of total bases covering the marker loci to be accepted

const int cAllocLociPValues = 100000;   // allocate for putative SNP loci in this many increments

const int cAllocLineBuffSize = 0x01fffffff; // 512MB buffer - when writing to results file then allow for buffering up to this many chars so as to reduce write frequency

const int cDfltMaxMultiHits = 5;		// default is to process at most this number of per read multihits
const int cMaxMultiHits = 500;			// user can specify at most this many multihits
const int cMaxAllHits = 100000;			// but if reporting all multihit loci then limit is increased to this value

const int cAllocMultihits = 25000000;		// alloc/realloc for multihit loci in this many instance increments
const int cDfltReadLen = 200;			 // assume reads plus descriptors combined of of this length - not critical as actual read lengths are processed
const size_t cReadsHitReAlloc = 50000000; // realloc allocation to hold this many read instances

const int cPairMinLen =	25;				// apparent paired reads sequences must be of at least this minimum length
const int cPairMaxLen = 100000;			// apparent paired reads sequences are restricted be of this maximum length
const int cDfltPairMinLen = 100;		// default apparent paired reads sequences insert size must be of at least this length
const int cDfltPairMaxLen = 1000;		// default apparent paired reads sequences insert size must be no longer than this length

const int cMaxConstrainedChroms = 64;   // at most this many chroms can have loci base constraints
const int cMaxConstrainedLoci = (cMaxConstrainedChroms * 100); // allow an average of 100 constraints per constrained chrom

// following constants are very empirical, relate to determination of multimatch read loci association, and need to be fine tuned..
const UINT16 cUniqueClustFlg = 0x08000;	// used to flag read as being clustered to reads which are uniquely aligned
const int cClustMultiOverLap = 10;		// to be clustered, reads must overlap by at least this number of bp 
const int cClustUniqueScore = 5;		// unique reads given this score if in current cluster window
const int cClustMultiScore = 1;			// multimatch read loci given this score if in current cluster window
const int cClustScaleFact = 10;			// scores are scaled down by this factor
const UINT32 cMHminScore = 50;			// any putative multimatch alignment score must be at least this to be accepted as the alignment for that read

const int cReadHitBuffLen = 0xfffff;		// sets per thread buffer size for holding string output read hit records
const int cDataBuffAlloc = 0x0fffffff;		// alloc to hold reads in this byte sized increments
const int cRdsBuffAlloc =   0x07fffff;		// alloc to hold preprocessed reads (for stats) in this byte sized allocation

const int cMaxDescrLen = 128;				// allow for descriptors of upto this length

const unsigned int cMinSeqLen = 15;			// sequences must be at least this length otherwise user is warned and sequence sloughed
const unsigned int cDfltMinAcceptReadLen = 50; // by default reads after any end trimming must be at least this length to be accepted for alignment processing
const unsigned int cDfltMaxAcceptReadLen = 500; // by default reads after any end trimming must be no longer than this length to be accepted for alignment processing
const unsigned int cMaxSeqLen = 2000;	// sequences must be no longer than this length otherwise user is warned and sequence sloughed

const int cPCRPrimerSubs = 5;				// user can specify for upto this many PCR hexamer primer subs in 5' flank over 12bp

const int cDfltRelSiteStartOfs = -4;		// default relative octamer start site offset, default gives 4bp 5' of read start
const int cMaxSitePrefOfs = 100;			// allow octamer starts to be at most this away from read start sites 

const int cNumOctamers = 0x010000;			// number of different octamers (4^8) in m_OctSitePrefs[]

#pragma pack(1)

// each read will be marked with the reason as to why that read was not accepted as being aligned
typedef enum eNAR {
	eNARUnaligned = 0,			// read has yet to be aligned
	eNARAccepted,				// read has been accepted as being aligned
	eNARNs,						// not accepted because contains too many indeterminate bases
	eNARNoHit,					// not accepted as aligned because unable to find any potential hits
	eNARMMDelta,				// not accepted as aligned because MMDelta criteria not met
	eNARMultiAlign,				// not accepted as aligned because was multiloci aligned
	eNARTrim,					// not accepted as aligned because aligned read excessively trimmed
	eNARSpliceJctn,				// not accepted as aligned because aligned read orphan splice junction
	eNARmicroInDel,				// not accepted as aligned because aligned read orphan microInDel
	eNARPCRdup,					// not accepted as aligned because aligned read was a PCR duplicate
	eNARNonUnique,				// not accepted as aligned because aligned read was a duplicate sequence
	eNARChromFilt,				// not accepted as aligned because aligned to a filtered chrom
	eNARRegionFilt,				// not accepted as aligned because not aligned to a priority region
	eNARPEInsertMin,			// not accepted as aligned because both PE ends align but less than minimum insert size
	eNARPEInsertMax,			// not accepted as aligned because both PE ends align but more than maximum insert size
	eNARPENoHit,				// not accepted as aligned because PE aligning and although able to align this read was unable to align partner end
	eNARPEStrand,				// not accepted as aligned because PE aligning and although able to align this read other read was aligned to inconsistent strand
	eNARPEChrom,				// not accepted as aligned because PE aligning and although able to align this read other read was aligned to different chromosome
	eNARPEUnalign,				// not accepted as aligned because PE aligning and unable to accept this alignment
	eNARLociConstrained,		// not accepted as aligned because alignment violated loci base constraints
	eNARundefined				// used as place holder setting 
} teNAR;

// read alignment non-acceptance status descriptive text
typedef struct TAG_sNAR {
	teNAR NAR;		// NAR enumeration
	char *pszNAR;   // NAR abbreviation
	char *pszNARdescr; // NAR descriptive text 
} sNAR;

typedef enum eHLalign {
	eHLunique = 0,				// read aligned to a unique hit loci
	eHLrandom = 1,				// read multialigned and loci was randomly choosen
	eHLclustunique = 2,			// read multialigned and loci near other uniquely aligned reads was choosen
	eHLclustany = 3			    // read multialigned and loci near other multialigned reads was choosen
} teHLalign;

typedef struct TAG_sQScoreDist {
	UINT32 Subs;				// count of aligner induced substitutions
	UINT32 QInsts;				// count of instances of this Phred quality score
} tsQScoreDist;


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
	UINT16 TrimLeft;			// left flank trimming removes this many bases
    UINT16 TrimRight;			// right flank trimming removes this many bases
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
	UINT32 PrevSizeOf;			// size of the previously loaded tsReadHit - allows easy referencing of partner pairs
	UINT32 ReadHitIdx;			// current read hit index + 1 for this read
	UINT32 ReadID;				// read identifier from the preprocessed read (tsRawReadV5)
	UINT32 PairReadID;			// this read's paired raw read identifier (0 if not a paired read) bit31 reset if 5' fwd read, bit31 set if 3' read of pair
	UINT32 NumReads;			// number of source reads merged into this read
	UINT8 DescrLen;				// descriptor length
	UINT16 ReadLen;				// read length of sequence packed into Read following the descriptor
	INT16  LowHitInstances;		// number of hit target loci instances at LowMMCnt mismatches
	INT8   LowMMCnt;			// lowest number of mismatches for this read thus far
	INT8   NxtLowMMCnt;			// next lowest number of mismatches for this read thus far
	INT16  NumHits;				// number of target hits for this read, currently limited to be at most 1
	UINT8  NAR:6;				// if not eNARAccepted then reason (etNAR) for this read not being accepted as being aligned
	UINT8  FlgPEAligned:1;      // PE read accepted as PE alignment
	tsReadHitLoci HitLoci;		// this is currently the best hit for this read
	UINT8  SiteIdx;				// read has been characterised into this index into m_OctSitePrefs[]
	UINT8  Read[1];				// packed read descriptor, sequence and quality values for this read
	} tsReadHit;

// alignments can be filtered if they are to a target sequence containing loci which are defined by the user as being base constrained
// reads which align to a constrained target loci are checked against the constraint and treated as unaligned if not meeting the constaint 
typedef struct TAG_sConstraintLoci {
	UINT32 ChromID;				// constraint applies to this target chrom/transcript
	UINT32 StartLoci;			// starting from this loci inclusive
	UINT32 EndLoci;				// to this loci inclusive
	UINT8  Constraint;			// constraint - must be a combination of one of the following bit mapped bases A: 0x01, C: 0x02:, G: 0x04, T: 0x08, target base: 0x10
} tsConstraintLoci;

// processing modes
typedef enum TAG_ePMode {
	ePMdefault,					// default processing mode
	ePMMoreSens,				// more sensitive - slower
	ePMUltraSens,				// ultra sensitive - much slower
	ePMLessSens,				// less sensitive - quicker
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// multiloci read handling modes
typedef enum TAG_eMLMode {
	eMLdefault,					// default is to simply slough those reads which match to multiple loci
	eMLdist,					// acumulate distribution stats only
	eMLrand,					// randomly select one of the aligned loci
	eMLuniq,					// cluster multiply aligned reads with those reads which are uniquely aligned
	eMLcluster,					// cluster multiply aligned reads with both uniquely (high priority) aligned and then with other multiloci alignments
	eMLall,						// report all multihit loci
	eMLplaceholder				// used to set the enumeration range
} etMLMode;

// output format modes
typedef enum TAG_eFMode {
	eFMdefault,					// default is for CSV match loci only
	eFMmatch,					// CSV loci + match sequence
	eFMread,					// CSV loci + read sequence
	eFMreadmatch,				// CSV loci + read + match sequence
	eFMbed,						// UCSC BED format
	eFMsam,						// SAM toolset format
	eFMsamAll,					// SAM toolset format, includes all reads even if not accepted as aligned
	eFMplaceholder				// used to set the enumeration range
	} etFMode;

// SAM can be generated as SAM, BAM, or BAM bgzf compressed
typedef enum eSAMFormat {
	etSAMFformat = 0,			// output SAM as SAM
	etSAMFBAM					// output as BAM compressed with bgzf
} teSAMFormat;


// paired reads alignment processing mode
typedef enum TAG_ePEproc {
	ePEdefault = 0,				// default is not to process for paired ends, single end alignment only
	ePEorphan,					// process for paired ends and if one end is orphaned because multiple aligned then try to locate unique alignment downstream
	ePEunique,					// process for paired ends but only accept putative if both ends uniquely aligned
	ePEorphanSE,				// ePEorphan processing but if unable to accept as PE then process ends independently as if SE
	ePEuniqueSE,				// ePEunique processing but if unable to accept as PE then process ends independently as if SE
	ePEplaceholder				// used to set the enumeration range
} etPEproc;

typedef enum TAG_eQFiltType {
	eQFnone,					// no filtering by quality scores
	eQFabs,						// filter by absolute quality scores
	eQFrel,						// filter by relative quality scores
	eQFplaceholder				// used to set the enumeration range
	} etQFiltType;

typedef enum TAG_eReadsSortMode {
	    eRSMunsorted = 0,		// reads are initially unsorted
		eRSMReadID,				// index by ascending ReadID
		eRSMPairReadID,			// index by ascending PairReadID
		eRSMPEHitMatch,			// index by ascending hit count, chrom, PairReadID
		eRSMHitMatch,			// index by ascending hit count,chrom, loci, strand, level
		eRSMSeq,				// index by ascending sequence then ReadID
		eRSMplaceholder			// used to limit the enumeration range
} etReadsSortMode;

typedef struct TAG_sSNPcnts {
	etSeqBase RefBase;	// reference base
	UINT32 NumRefBases;		// counts of reference base in reads covering this loci
	UINT32 NumNonRefBases;	// total count of non-reference bases in reads covering this loci
	UINT32 NonRefBaseCnts[5]; // counts of non-reference bases a,c,g,t,n covering this loci
	} tsSNPcnts;

typedef struct TAG_sAdjacentSNPs {
	int StartLoci;		// currently iterating reads which overlap between StartLoci and	
	int EndLoci;		// EndLoci inclusive
	tsReadHit *pFirstIterReadHit; // this read was the first iterated returned read overlapping StartLoci and EndLoci
	tsReadHit *pPrevIterReadHit; // this read was the previously iterated returned read overlapping StartLoci and EndLoci
	} tsAdjacentSNPs;

typedef struct TAG_sChromSNPs {
	UINT32 ChromID;		// uniquely identifies this chromosome
	UINT32 ChromLen;	// this chromosome length
	tsAdjacentSNPs AdjacentSNPs[2]; // allowing for both DiSNPs and TriSNPs
 	tsReadHit *pFirstReadHit; // 1st read on chromosome which was accepted for SNP processing
	tsReadHit *pLastReadHit; // last read on chromosome which was accepted for SNP processing
	UINT32 AllocChromLen; // cnts can be for a chromosome of at most this length
	INT64 TotMatch;	// total number of aligned read bases which exactly matched corresponding chrom sequence bases
	INT64 TotMismatch;	// total number of aligned read bases which mismatched corresponding chrom sequence base
	UINT32 MeanReadLen;  // mean length of all reads used for identifying putative SNPs, determines max separation used for Di/TriSNP counts
	UINT32 NumReads;     // number of reads used for identifying putative SNPs from which MeanReadLen was calculated 
	UINT64 TotReadLen;	 // total length, in bp, of all reads used for identifying putative SNPs from which MeanReadLen was calculated
	tsSNPcnts Cnts[1]; // will be allocated to hold base cnts at each loci in this chromosome
	} tsChromSNPs;

typedef struct TAG_sSegJuncts {
	tsReadHit *pRead;	// read containing this RNA-seq splice or microInDel junction
	UINT32 Cnt;			// number of reads sharing this splice junction or microInDel junction
	UINT32 ChromID;		// junction is on this chrom
	UINT32 Starts;		// junction starts
	UINT32 Ends;		// junction ends
	} tsSegJuncts;


typedef struct TAG_sLociPValues {
	UINT32 Loci;		// putative SNP at this loci
	double PValue;      // having this PValue
	UINT32 Rank;		// and this ordered rank
	double LocalBkGndSubRate; // local background substitution rate
	UINT32 LocalReads;  // total number of aligned bases within the local background
	UINT32 LocalSubs;	// total number local aligner induced substitutions within the local background
	tsSNPcnts SNPcnts;	// counts of each base a,c,g,t,n at this loci
	UINT32 NumReads;	// number of reads aligned at this loci
	UINT32 NumSubs;		// number of aligner induced substitutions at this loci
	UINT32 MarkerID;	// generated marker sequence will have this identifier as '>Marker<MarkerID>'
	UINT32 NumPolymorphicSites; // number of polymorthic sites within the marker sequence
} tsLociPValues;

typedef struct TAG_sSNPCentroid {
	UINT32 CentroidID;		// uniquely identifies this centroid sequence
	UINT32 NumInsts;		// number of instances of this centroid sequence (qualified by having at least min number of reads covering)
	int NumSNPs;			// number of SNPs identified centroid to this sequence
	UINT32 RefBaseCnt;		// counts of reference bases covering SNP
	UINT32 NonRefBaseCnts[5]; // counts of each non-reference bases covering SNP
} tsSNPCentroid;

#pragma pack()

#pragma pack(4)
typedef struct TAG_sThreadMatchPars {
	int ThreadIdx;					// index of this thread (1..m_NumThreads)
	void *pThis;					// will be initialised to pt to CAligner instance
	int NumIdentNodes;				// number of ident nodes allocd for use by this thread
	tsIdentNode *pIdentNodes;		// thread to use these nodes

#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int CurBlockID;					// current suffix block identifier
	int ChromID;					// hit chrom identifier
	int MinChimericLen;				// if checking for chimerics then minimim percentage of read length required, set to 0 if not checking for chimerics
    int NumAllowedSubs;				// number of allowed substitutions per 100bp of read length
	int MaxNumSlides;				// limit on number of times core window can be moved or slide to right over read per 100bp of read length
	int	MinCoreLen;					// minimum core length allowed
    etQFiltType FiltQuality;		// filtering type
	bool bQAbove;					// if true then filter out those below or equal to FiltThres
	int FiltThres;					// filter threshold
	eALStrand AlignStrand;			// align on to watson, crick or both strands of target
	int microInDelLen;				// microInDel length maximum
	int SpliceJunctLen;				// maximum splice junction length when aligning RNAseq reads
	int Rslt;						// returned result code
	int MinEditDist;				// any matches must have at least this Hamming edit distance to the next best match
	int MaxSubs;					// maximum number of substitutions allowed per 100bp of actual read length
	int PlusHits;					// returned number of hits on to plus strand
	int MinusHits;					// returned number of hits on to minus strand
	int ChimericHits;				// returned number of hits which were chimeric
	int NumReadsProc;				// returned number of reads processed by this thread instance
	int OutBuffIdx;					// index at which to write next formated hit into szOutBuff
	UINT8 *pszOutBuff;				// used to buffer multiple hit formated output records prior to writing to disk
	tsHitLoci *pMultiHits;			// allocated to hold read multihit loci
} tsThreadMatchPars;

typedef struct TAG_sClusterThreadPars {
	int ThreadIdx;						// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to CAligner instance

#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// returned result code
} tsClusterThreadPars;

typedef struct TAG_sLoadReadsThreadPars {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to CAligner instance

#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	UINT32 SampleNthRawRead;		// sample every Nth raw read (or read pair) for processing (1..10000)
	int *pRslt;						// write intermediate result codes to this location
	int Rslt;						// returned result code
} tsLoadReadsThreadPars;


typedef struct TAG_sReadsHitBlock {
	int NumReads;			// number of reads for processing in this block
	int MaxReads;			// block can hold at most this number of reads
	tsReadHit *pReadHits[cMaxReadsPerBlock]; // reads for processing
} tsReadsHitBlock;
#pragma pack()

// SOLiDmap
// Used for mapping from base to colorspace and the reverse
static UINT8 SOLiDmap[5][5] = {
	{0,1,2,3,4},	// a
	{1,0,3,2,4},	// c
	{2,3,0,1,4},    // g
	{3,2,1,0,4},    // t
	{0,1,2,3,4}};	// n

typedef struct TAG_sOctSitePrefs {
	int Octamer;	// identifies this octamer (0 == bases:aaaaaaaa, 1 == bases:aaaaaaac, 0xFFFF == tttttttt)
	int NumSites;	// number of unique loci at which at least 1 read starting with this octamer aligned
	int NumOccs;   // total number of reads aligning to assembly with this start octamer
	double RelScale; // relative normalisation scaling factor to use
	} tsOctSitePrefs;



class CAligner
{

	CMTqsort m_mtqsort;				// muti-threaded qsort

	CContaminants *m_pContaminants; // for use when trimming reads containing contaminants

	tsBSFRdsHdr m_FileHdr;			// processed reads file header

	bool m_bPEInsertLenDist;		// true if stats file to include PE insert length distributions for each transcript
	bool m_bPEcircularised;			// CAUTION: experimental - true if processing for PE spaning circularised fragments
	bool m_bAllReadsLoaded;			// set true when all reads have been parsed and loaded
	teBSFrsltCodes m_LoadReadsRslt;	// set with exit code from background reads load thread, read after checking if m_bAllReadsLoaded has been set
	size_t m_DataBuffOfs;			// offset at which to read in next read
	UINT32 m_NumDescrReads;			// number of reads thus far parsed

	tsReadHit *m_pReadHits;			// memory allocated to hold reads, reads are written contiguously into this memory
									// caution: m_pReadHits is allocated/freed with malloc/realloc/free or mmap/mremap/munmap
	size_t m_AllocdReadHitsMem;		// how many bytes of memory  for reads have been allocated
	size_t m_UsedReadHitsMem;		// how many bytes of allocated reads memory is currently used
	UINT32 m_NumReadsLoaded;		// m_pReadHits contains this many reads
	UINT32 m_OrigNumReadsLoaded;	// if multiloci read alignments being treated as if each was a separate read then this is a copy of m_NumReadsLoaded prior to overwriting with the multiloci read count 
	UINT32 m_FinalReadID;			// final read identifier loaded as a preprocessed read (tsProcRead)
	UINT32 m_PrevSizeOf;			// size (UINT8's) of the previously loaded tsReadHit - allows easy referencing of partner pairs

	tsReadHit **m_ppReadHitsIdx;	// memory allocated to hold array of ptrs to read hits in m_pReadHits - usually sorted by some critera
	UINT32 m_AllocdReadHitsIdx;		// how many elements for m_pReadHitsIdx have been allocated
	etReadsSortMode	m_CurReadsSortMode;	// sort mode last used on m_ppReadHitsIdx

	size_t m_AllocLociPValuesMem;   // total memory currently allocated to m_pLociPValues
	UINT32 m_NumLociPValues;		// current number of LociPValues
	tsLociPValues *m_pLociPValues; // allocated to hold putative SNP loci and their associated PValues
	double m_QValue;				// QValue used in

	etMLMode m_MLMode;				// how to process multiloci matching reads
	bool m_bClampMaxMLmatches;		// normally if too many multiloci hits for a read that read is not accepted as aligned; if m_bAcceptAllMultihits true then these reads are accepted as aligned
	bool m_bLocateBestMatches;		// align for best matches, not just those 1H better than any other, upto at most MaxMLmatches

	int m_ClustWinSize;				// multimatch clustering window size
	int m_SMscore;					// score for uniquely matched reads in clustering window
	int m_MMscore;					// score for multimatched matched reads in clustering window
	tsReadHit *m_pMultiHits;		// holds multihit loci instances
	UINT32 m_AllocdMultiHits;		// how many elements for m_pMultiHits have been allocated
	size_t m_AllocdMultiHitsMem;	// how much memory for m_pMultiHits has been allocated
	UINT32 m_NumMultiHits;			// number of multihit loci elements currently used
	UINT32 m_NumUniqueMultiHits;	// number of multihit reads which were initially mapped to a single loci
	UINT32 m_NumProvMultiAligned;	// number of reads provisionally with more than one aligned to loci
	UINT32 m_NARAccepted;			// number of reads being reported as accepted aligned

	UINT32 m_NumMultiAll;				// number of alignments in m_pMultiAll
	size_t m_NxtMultiAllOfs;			// offset in m_pMultiSAM at which to append next alignment
	size_t m_AllocMultiAllMem;			// current allocation size 			
	tsReadHit *m_pMultiAll;				// allocated memory for holding multiloci alignments which are all to be reported

	UINT32 m_NumSloughedNs;			// number of reads which were not aligned because they contained excessive number of indeterminate bases
	UINT32 m_TotNonAligned;			// number of reads for which no alignment was discovered
	UINT32 m_TotAcceptedAsAligned;	// number of reads accepted as being aligned
	UINT32 m_TotAcceptedAsUniqueAligned;  // number of reads accepted as aligned which were uniquely aligned
	UINT32 m_TotAcceptedAsMultiAligned;   // number of reads accepted as aligned which aligned to multiple loci


	UINT32 m_TotAcceptedHitInsts;	// number of reads aligned and accepted (included in m_TotAcceptedAsAligned) even though there were more than m_MaxMLmatches instances
	UINT32 m_TotLociAligned;		// number of loci aligned which have been reported
	UINT32 m_TotNotAcceptedDelta;	// number of reads aligned but not accepted because of insufficient hamming

	UINT32 m_SampleNthRawRead;		// sample every Nth raw read (or read pair) for processing (1..10000)
	UINT32 m_MaxReadsLen;			// longest read processed
	UINT32 m_MinReadsLen;			// shortest read processed
	UINT32 m_AvReadsLen;			// average length read processed

	int m_ThreadCoredApproxRslt;    // normally >= 0, but will be set < 0 if any errors by threads within the threaded cored approximate alignment processing

	int m_MinSNPreads;				// before SNP can be called there must be at least this number of reads covering the loci
	double m_SNPNonRefPcnt;			// only process for SNP if more/equal than this percentage number of reads are non-ref at putative SNP loci (defaults to 25) 
	int m_MaxDiSNPSep;				// putative DiSNPs only processed if separation between the SNPs <= this many bp

	tBSFEntryID m_PrevSAMTargEntry; // used to determine when generating SAM output if the target chrom has already been loaded
	char m_szSAMTargChromName[128];	// holds previously loaded SAM chrom

	UINT64 m_BlockTotSeqLen;		// total sequence length in currently loaded suffix block
	int	m_MinCoreLen;				// minimum core length allowed

	eALStrand m_AlignStrand;		// which strand to align to - both, watson '+' or crick '-'


	int m_NumConstrainedChroms;		// number of chroms having at least one constraint
	UINT32 m_ConstrainedChromIDs[cMaxConstrainedChroms];    // chromosomes having at least one loci base constraint 

	int m_NumConstraintLoci;		// number of constrained alignment loci in m_pConstraintLoci[]
	tsConstraintLoci *m_pConstraintLoci; // allocated to hold any read alignment loci constraints

	int m_MinChimericLen;			// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence
	bool m_bReportChimerics;        // true if individual chimeric sequences to also be reported for diagnostics
	int m_microInDelLen;			// microInDel length maximum
	int m_SpliceJunctLen;			// minimum splice junction length when aligning RNAseq reads

	int m_PerThreadAllocdIdentNodes;    // each thread can use this many tsIdentNodes
	int m_TotAllocdIdentNodes;			// total number of tsIdentNodes allocated
	tsIdentNode *m_pAllocsIdentNodes;	// memory allocated to hold tsIdentNodes required by all threads

	teSAMFormat m_SAMFormat;		// output SAM as SAM, BAM or BAM compressed with bgzf

	int m_hInFile;			// input file handle

	int m_hOutFile;			// results output file handle
	int m_hBAIFile;			// used when outputing alignments as BAM, this handle is for the BAM index file

	int m_hIndOutFile;		// optional output microInDel file handle
	int m_hJctOutFile;		// optional output splice junction file handle
	int m_hSitePrefsFile;	// optional output site octamer preferencing file
	int m_hStatsFile;		// optional output stats file handle
	int m_hInsertLensFile;  // optional output PE inserts dist file handle
	int m_hNoneAlignFile;	// optional output none-alignable file handle
	int m_hMultiAlignFile;	// optional output multialigned reads file handle
	int m_hSNPfile;			// file handle used if SNPs are being processed
	int m_hDiSNPfile;       // file handle used if DiSNPs are being processed
	int m_hTriSNPfile;       // file handle used if TriSNPs are being processed
	int m_hSNPCentsfile;	// file handle used if SNP centroids are being processed

	gzFile m_gzOutFile;			// results output when compressing the output as gzip
	gzFile m_gzIndOutFile;		// optional output microInDel when compressing the output as gzip
	gzFile m_gzJctOutFile;		// optional output splice junction when compressing the output as gzip
	gzFile m_gzNoneAlignFile;	// optional output none-alignable when compressing the output as gzip
	gzFile m_gzMultiAlignFile;	// optional output multialigned reads when compressing the output as gzip
	gzFile m_gzSNPfile;			// when compressing the output as gzip used if SNPs are being processed
	gzFile m_gzSNPCentsfile;	// fwhen compressing the output as gzip if SNP centroids are being processed

	char m_szDiSNPFile[_MAX_PATH];		// DiSNP results to this file
	char m_szTriSNPFile[_MAX_PATH];		// TriSNP results to this file
	char m_szIndRsltsFile[_MAX_PATH];		// microIndel results to this file
	char m_szJctRsltsFile[_MAX_PATH];		// splice junction results to this file

	tsSNPCentroid *m_pSNPCentroids;		// allocated to hold centroid distributions if SNP processing

	char *m_pszSitePrefsFile;		// site octamer preferencing to this file
	int m_SitePrefsOfs;				// offset read start sites when processing  octamer preferencing, range -100..100

	bool m_bgzOutFile;				// true if alignments to be generated compressed with gzopen/gzwrite/gzclose
	bool m_bgzNoneAlignFile;		// true if none aligned to be generated compressed with gzopen/gzwrite/gzclose
	bool m_bgzMultiAlignFile;		// true if multialigned to be generated compressed with gzopen/gzwrite/gzclose

	etPEproc m_PEproc;				// paired reads alignment processing mode
	int m_NumPE2InputFiles;			// number of input PE2 file specs
	char **m_ppszPE2InputFiles;		// optional raw PE2 paired reads are in these files
	char *m_pszOutFile;				// write alignments to this file
	char *m_pszOutBAIFile;			// write BAM indexes to this file
	char *m_pszSNPRsltsFile;		// write SNPs to this file
	char *m_pszSNPCentroidFile;		// write SNP centroids to this file

	char *m_pszSfxFile;				// target as suffix array
	char *m_pszStatsFile;			// aligner induced substitutions stats file or paired end length distributions
	char *m_pszPEInsertDistFile;	// PE insert size distributions for each transcript or assembly sequence
	char *m_pszMultiAlignFile;		// file to contain reads which are aligned to multiple locations
	char *m_pszNoneAlignFile;		// file to contain reads which were non-alignable

	int m_NumPE1InputFiles;			// number of input PE1 file specs
	char **m_ppszPE1InputFiles;		// names of inputs file (wildcards allowed unless in dump mode) containing raw reads

	char *m_pszMarkerFile;			// marker file
	int m_hMarkerFile;				// if outputing marker sequences then the file handle to use
	int m_MarkerID;					// current marker identifier
	int m_Marker5Len;				// marker sequence 5' flank to be of this length
	int m_Marker3Len;				// marker sequence 3' flank to be of this length
	double m_MarkerPolyThres;		// maximum allowed marker sequence base polymorphism independent of centroid SNP (default 0.333, range 0.0 to 0.5)



	tsChromSNPs *m_pChromSNPs;		// allocated for use in SNP processing
	int m_TotNumSNPs;				// total number of SNPs discovered
	
	INT64 m_LociBasesCovered;		// total number of targeted loci (bases) covered by aligned reads when SNP processing - could be used for to determine fold coverage
	INT64 m_LociBasesCoverage;		// total number of read bases covering aligned to loci bases 

	int m_szLineBuffIdx;			// offset into m_pszLineBuff at which to next write
	char *m_pszLineBuff;			// allocated to hold output line buffering

	CSfxArrayV3 *m_pSfxArray;		// suffix array holds genome of interest
	char m_szTargSpecies[cMaxDatasetSpeciesChrom+1]; // suffix array was generated over this targeted species

	CBEDfile *m_pPriorityRegionBED;	// to hold exact match priority regions
	bool m_bFiltPriorityRegions;	// if priority regions requested then can also request that only alignments into these regions be reported

	tsOctSitePrefs m_OctSitePrefs[2][cNumOctamers];	// to contain all 64k possible octamers for each aligned to strand

	int *m_pLenDist;				// allocated to hold paired read sequence length counts

	int m_Trim5;					// trim this number of bases from 5' end of reads when loading the reads
	int m_Trim3;					// trim this number of bases from 3' end of reads when loading the reads

	int m_MinAcceptReadLen;					// only accepting reads for alignment if at least this length after any end trimming
	int m_MaxAcceptReadLen;					// only accepting reads for alignment if no longer than this length after any end trimming

	int m_MaxSubs;					// accepted aligned reads must have at most this many subs per 100bp aligned sequence
	int m_InitalAlignSubs;			// putatively accepted as aligned reads can have at most this many subs per 100bp aligned sequence

	int m_MaxRptSAMSeqsThres;		// report all SAM chroms or sequences if number of reference chroms <= this limit (defaults to 10000)

		// note that sub distributions are characterised into 4 bands: Phred 0..8, 9..19, 20..29, and 30+
	tsQScoreDist m_AlignQSubDist[4][cMaxFastQSeqLen+1];	// to hold aligner induced substitution count distribution
	int m_AlignMSubDist[cMaxFastQSeqLen+1];  // to hold multiple aligner induced substitutions counts
	int m_MaxMSubDist;		// actual maximum number of multiple aligner induced substitutions in any read
	int m_MaxAlignLen;		// actual maximum length aligned read length

	etPMode m_PMode;		// processing mode
	etFMode m_FMode;		// output format mode
	bool m_bSNPsVCF;		// true if SNPs reporting as VCF instead of the default CSV
	int m_NumThreads;		// number of worker threads to use
	int m_MaxNs;			// max number of indeterminate bases 'N' to acept in read before deeming read as unalignable
	int m_MaxMLmatches;			// allow at most this many multihits by any single read before accepting read as being aligned
	int m_MultiHitDist[cMaxMultiHits];	// used to record the accepted as aligned multihit distribution

	tsHitLoci *m_pAllocsMultiHitLoci; // allocated to hold all multihit loci for all threads
	char *m_pszTrackTitle;			// track title if output format is UCSC BED
	UINT8 *m_pAllocsMultiHitBuff;	  // allocated to hold per thread buffered output when processing all multihit loci

	int m_NumIncludeChroms;				// number of RE include chroms
	int m_NumExcludeChroms;				// number of RE exclude chroms

	#ifdef _WIN32
	Regexp *m_IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	Regexp *m_ExcludeChromsRE[cMaxExcludeChroms];
	#else
	regex_t m_IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	regex_t m_ExcludeChromsRE[cMaxExcludeChroms];
	#endif

	UINT32 m_NumReadsProc;		// number of reads thus far processed - note this is total reads handed out to processing threads
								// and should be treated as a guide only
	size_t m_NxtReadProcOfs;	// byte offset into m_pReadHits of next read to be processed

	bool m_bBisulfite;			// true if bisulfite methylation patterning processing
	bool m_bIsSOLiD;			// true if SOLiD or colorspace processing
	etFQMethod m_QMethod;		// fastq quality value method

	int m_ElimPlusTrimed;		// number of aligned reads flank trimmed on the '+' strand
	int m_ElimMinusTrimed;		// number of aligned reads flank trimmed on the '-' strand

	bool m_bMutexesCreated;		// will be set true if synchronisation mutexes have been created

	unsigned long m_ProcessingStartSecs;	
	UINT8 m_TermBackgoundThreads; // if non-zero then all background threads are to immediately terminate processing

	int m_ThreadLoadReadsRslt;

	UINT32 m_CurClusterFrom;

	static sNAR m_NARdesc[eNARundefined];	// NAR deescriptive text

#ifdef _WIN32
	unsigned int m_ThreadLoadReadsID;
#else
	pthread_t m_ThreadLoadReadsID;
#endif


	void Init(void);			// initialise state to that immediately following construction
	void Reset(bool bSync);		// reset state, if bSync true then fsync before closing output file handles

	char *ReplaceTabs(char *pszTabTxt); // Inplace replacement of any tabs with a single space char

	int AppendStr(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if != '\0' then prefix with this separator (usually ',' or '\t')
		  char Quote,		// if != '\0' then quote string with this char (usually single or double quote char)
		  char *pStr,		// '\0' terminated string
		  char TrailSep);	// if != '\0' then suffix with this separator (usually ',' or '\t' or '\n')

	int AppendChrs(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if != '\0' then prefix with this separator (usually ',' or '\t')
		  char Quote,		// if != '\0' then quote chars with this char (usually single or double quote char)
		  int NumChrs,		// number of chars to append
		  char *Chrs,		// pts to chars to append
		  char TrailSep);	// if != '\0' then suffix with this separator (usually ',' or '\t' or '\n')

	int							// length written
		AppendUInt(char *pszBuff,	// write to this buffer
		  char LeadSep,		// if > '\0' then prefix with this separator (usually ',' or '\t')
		  UINT32 Value,
		  char TrailSep);	// if > '\0' then suffix with this separator (usually ',' or '\t' or '\n')


	teBSFrsltCodes LoadLociConstraints(char *pszLociConstraints);	// load loci constraints from file

	int LoadReads(char *pszRdsFile);	// file containing preprocessed reads (genreads output)

	teBSFrsltCodes Disk2Hdr(char *pszRdsFile);	// read from disk and validate header

	// Following DiagSequences and LogReadHit are purely for use during low level debugging
#ifdef USEWHILEDBUG1
	char *DiagCmpSequences(int SeqLen,etSeqBase *pProbeSeq,etSeqBase *pTargSeq);
	char *DiagSequences(int ProbeLen,etSeqBase *pProbeSeq,char Strand,int ChromID,int MatchLoci,int MatchLen);
	void LogReadHit(tsReadHit *pReadHit);
#endif


	static UINT32 AdjStartLoci(tsSegLoci *pSeg);
	static UINT32 AdjEndLoci(tsSegLoci *pSeg);
	static UINT32 AdjHitLen(tsSegLoci *pSeg);
	static UINT32 AdjAlignStartLoci(tsHitLoci *pHit);
	static UINT32 AdjAlignEndLoci(tsHitLoci *pHit);
	static UINT32 AdjAlignHitLen(tsHitLoci *pHit);
	static int AdjAlignMismatches(tsHitLoci *pHit);

// NOTE: will only return SNP bases, e.g. aligned read sequence must not contain microInDel or span splice junctions
	static eSeqBase AdjAlignSNPBase(tsReadHit *pReadHit,	// aligned read
		   UINT32 ChromID,			// read expected to have aligned to this chromosome
			UINT32 Loci);            // base to be returned is at this alignment loci, base will be complemented if antisense alignment

	int ReportChimerics(char *pszChimericSeqFile = NULL);			// report chimerically trimmed read sequences to file pszChimericSeqFile

	int AutoTrimFlanks(int MinFlankExacts); // Autotrim back aligned read flanks until there are at least MinFlankExacts exactly matching bases in the flanks

	int PCR5PrimerCorrect(int MaxAllowedSubRate, int KLen = 12); // correct substitutions in first Klen 5' read bases - assumed to be PCR random primer artefacts - until overall read meets MaxAllowedSubRate

	int NumUniqueAlignedReads(void);		// return count of uniquely aligned reads

	int DedupeNonuniqueReads(void);	// dedupe reads

	int ReducePCRduplicates(int WinLen);		// remove potential PCR artefacts

	int RemoveOrphanSpliceJuncts(int SpliceJunctLen);	// remove unsupported orphan splice junctions

	int	RemoveOrphanMicroInDels(int microInDelLen);	// remove any unsupported orphan microInDels

	bool					// true if chrom is accepted, false if chrom not accepted
		AcceptThisChromID(UINT32 ChromID);

	bool				  // true if read alignment meets any loci base constraints
		AcceptLociConstraints(tsReadHit *pReadHit);   // read alignment to check

	int								// -1: base not meeting constraints, 0: chromID has no constraints, 1: ChromID constrained but base accepted
		AcceptBaseConstraint(UINT32 ChromID,	  // base aligned to this chrom/sequence
                                  UINT32 Loci,		  // aligned to this loci
								  etSeqBase Base);	  // base in read
	int						// number of alignments violating a loci base constraint
		IdentifyConstraintViolations(bool bPEread); // true if processing for PE's, false if processing for SE


	int				// > 0: accepted with return value the fragment (insert) size
		AcceptProvPE(int PairMinLen,  // only accept paired reads with a combined sequence length of at least this
			 int PairMaxLen,  // only accept paired reads with a combined sequence length of no more than this
			 bool bPairStrand,				// accept paired ends if on same strand
			 tsReadHit *pFwdReadHit,
			 tsReadHit *pRevReadHit);

	int				// returned PE insert size, <= 0 if errors
		PEInsertSize(int PairMinLen,	// only accept paired reads with a combined sequence length of at least this
			 int PairMaxLen,			// only accept paired reads with a combined sequence length of no more than this
			 bool bPairStrand,			// accept paired ends if on same strand		
			 UINT8 PE1Strand,			// PE1 aligned on to this strand
		     UINT32 PE1StartLoci,		// PE read starts at this loci
		     UINT32 PE1EndLoci,			// PE1 read ends at this loci
             UINT8 PE2Strand,			// PE2 aligned on to this strand
		     UINT32 PE2StartLoci,		// PE2 read starts at this loci
		     UINT32 PE2EndLoci);		// PE2 read ends at this loci

	int ProcessPairedEnds(etPEproc PEproc, // paired reads alignment processing mode
				  int MinEditDist, // accepted alignments must be at least this Hamming away from other putative alignments
				  int PairMinLen,  // only accept paired reads with a combined sequence length of at least this
				  int PairMaxLen,  // only accept paired reads with a combined sequence length of no more than this
				  bool bPairStrand,	// accept paired ends if on same strand
				  int MaxSubs);	   // aligned reads can have at most this many subs per 100bp

	int ReportAlignStats(void);		// report basic alignment statistics

	
	UINT32							// Median insert length
		MedianInsertLen(UINT32 NumInserts,				// number of insert lengths in pInsertLens
				UINT32 *pInsertLens);		// insert lengths
	int ReportPEInsertLenDist(void);  // Report PE insert length distributions for each transcript or assembly contig/sequence

	int WriteSubDist(tsReadHit *pReadHit);

	int WriteReadHits(bool bPEProc);		   // true if processing paired ends

	// Write alignments as SAM or BAM format
	int WriteBAMReadHits(etFMode ProcMode,	   // eFMsam or eFMsamAll
							teSAMFormat SAMFormat, // if SAM output format then could be SAM or BAM compressed dependent on the file extension used
							bool bPEProc,		   // true if processing paired ends
							int ComprLev);		   // BAM to be BGZF compressed at the requested level (0..9)

	int ReportBAMread(tsReadHit *pReadHit,		// read to report
			int RefID,							// read aligns to this BAM refID (-1) if unaligned
			int  ReadIs,						// 0 if SE, 1 if PE1 of a PE, 2 if PE2 of a PE
			tsBAMalign *pBAMalign);				// BAM alignment to return

	// calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open)
	int BAMreg2bin(int beg, int end);

	// calculate the list of bins that may overlap with region [beg,end) (zero-based)
	int BAMreg2bins(int beg, int end, UINT16 *plist);

	int				// normally NumHits, but will be actual number of hits if unable to accept any loci hit because of chromosome filtering
		WriteHitLoci(tsThreadMatchPars *pThreadPars,tsReadHit *pReadHit,int NumHits,tsHitLoci *pHits);

	int ReportMultiAlign(void);
	int ReportNoneAligned(void);

	char *Octamer2Txt(int Octamer);		 // Report on site octamer site preferencing distribution

	int OutputSNPs(void);
	int ProcessSNPs(void);

	int ProcessSiteProbabilites(int RelSiteStartOfs); // offset the site octamer by this relative start offset (read start base == 0)
	int WriteSitePrefs(void);

	int FiltByChroms(void);

	int FiltByPriorityRegions(void); // remove hits not aligned into a priority region

	
	int ReportTargHitCnts(void);	 // if PEs then each end is separately counted
	
	int WriteBasicCountStats(void);


	bool								// true if file to be generated compressed with gzopen/gzwrite/gzclose
		FileReqWriteCompr(char *pszFile); // If last 3 chars of file name is ".gz" then this file is assumed to require compression
	int CreateOrTruncResultFiles(void);	// Create, or if file already exists then truncate, the multitude of results files which user may have requested

	int CompileChromRegExprs(int	NumIncludeChroms,	// number of chromosome regular expressions to include
		char **ppszIncludeChroms,		// array of include chromosome regular expressions
		int	NumExcludeChroms,			// number of chromosome expressions to exclude
		char **ppszExcludeChroms);		// array of exclude chromosome regular expressions

	int InitiateLoadingReads(void);		// reads are loaded asynchronously to the alignments

	int AssignMultiMatches(void); // false to cluster with uniques, true to cluster with multimatches

	int AddMultiHit(tsReadHit *pReadHit);

	int AddEntry(bool bIsPairRead,		// true if this is the paired read PE2
		 UINT32 PairReadID,		// identifies partner of this read if paired read processing
		 UINT8 FileID,			// identifies file from which this read was parsed
		 int DescrLen,			// length of following descriptor
		 char *pszReadDescr,	// copy of descriptor, used to pair reads with matching descriptors
		 int ReadLen,			// length of following read
		 UINT8 *pszReadBuff);	// packed read + phred score


	double										// returned prob of read being error free
		GenProbErrFreeRead(int QSSchema,		// guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger  
				   int ReadLen,					// read length
				   UINT8 *pQScores);			// Phred quality scores over the read


	teBSFrsltCodes LoadRawReads(bool bIsPairReads,				// true if paired end processing - PE1 reads in pszPE1File and PE2 reads in pszPE2File
		  int FileID,						// uniquely identifies source file for PE1, FileID + 1 uniquely identifies PE2 file
		  char *pszPE1File,					// process PE1 reads from this file
		  char *pszPE2File);					// optionally process PE2 reads from this file

	int CreateMutexes(void);
	void DeleteMutexes(void);
	void AcquireSerialise(void);
    void ReleaseSerialise(void);
    void AcquireSerialiseMH(void);
	void ReleaseSerialiseMH(void);
	void AcquireLock(bool bExclusive);
	void ReleaseLock(bool bExclusive);

	int LocateCoredApprox(int MinEditDist,	// any matches must have at least this edit distance to the next best match
		int MaxSubs);						// maximum number of substitutions allowed per 100bp of actual read length

	tsReadHit *LocateRead(UINT32 ReadID);	 // Locate read with requested ReadID

	int AddMHitReads(UINT32 NumHits,		// number of multimatches loci in pHits
		tsReadHit *pHits);					// pts to array of hit loci

	int SortReadHits(etReadsSortMode SortMode,		// sort mode required
				bool bSeqSorted = false,			// used to optimise eRSMSeq processing, if it is known that reads are already sorted in sequence order (loaded from pre-processed .rds file)
				bool bForce = false);				// if true then force sort

	void ResetThreadedIterReads(void);		 // must be called by master thread prior to worker threads calling ThreadedIterReads()

	UINT32		// Returns the number of reads thus far loaded and processed for alignment
		ApproxNumReadsProcessed(UINT32 *pNumProcessed,UINT32 *pNumLoaded);

	tsReadHit *IterReads(tsReadHit *pCurReadHit);			// to start from first read then pass in NULL as pCurReadHit
	tsReadHit *IterSortedReads(tsReadHit *pCurReadHit);		// iterate over sorted reads, to start from first read then pass in NULL as pCurReadHit
	tsReadHit *		// returned read which overlaps StartLoci and EndLoci, NULL if no read located
			IterateReadsOverlapping(bool bTriSNPs, // false if iterating DiSNPs, true if iterating TriSNPs
						tsChromSNPs *pChromSNPs, // processing SNPs on this chromosome
						int StartLoci,				// returned reads are required to overlap both this starting and
						int EndLoci);				// this ending loci

	bool	// returns false if no more reads availing for processing by calling thread
		ThreadedIterReads(tsReadsHitBlock *pRetBlock);

	int											// returned number of unique reads
		NumDnUniques(tsReadHit *pCurReadHit,		// current read
				int WinLen,					// only interested in unique reads starting within this window
				bool bStrandDep);			// if true then unique loci reads must be on current read stand

	int											// returned number of unique reads
		NumUpUniques(tsReadHit *pCurReadHit,		// current read
				int WinLen,					// only interested in unique reads starting within this window
				bool bStrandDep);			// if true then unique loci reads must be on current read stand

	int RunClusteringThreads(int NumThreads);
	int										// returns 0 if finished clustering or cnt of multihits to be processed by this thread
		GetClusterStartEnd(UINT32 *pMatchFrom,			// cluster from this inclusive index
					UINT32 *pMatchUntil);		// until this inclusive index


	static int SortReadIDs(const void *arg1, const void *arg2);
	static int SortPairReadIDs(const void *arg1, const void *arg2);
	static int SortPEHitMatch(const void *arg1, const void *arg2);
	static int SortHitMatch(const void *arg1, const void *arg2);
	static int SortMultiHits(const void *arg1, const void *arg2);
	static int SortMultiHitReadIDs(const void *arg1, const void *arg2);
	static int SortSegJuncts(const void *arg1, const void *arg2);
	static int SortReadSeqs(const void *arg1, const void *arg2);
	static int SortLociPValues(const void *arg1, const void *arg2);
	static int SortPValuesLoci(const void *arg1, const void *arg2);
	static int SortSiteRelScale(const void *arg1, const void *arg2);
	static int SortSiteRelOctamer(const void *arg1, const void *arg2);
	static int SortConstraintLoci(const void *arg1, const void *arg2);

#ifdef _WIN32
	HANDLE m_hMtxIterReads;
	HANDLE m_hMtxMHReads;
	HANDLE m_hMtxMultiMatches;
	SRWLOCK m_hRwLock;
	HANDLE m_hThreadLoadReads;
#else
	pthread_mutex_t m_hMtxIterReads;
	pthread_mutex_t m_hMtxMHReads;
	pthread_mutex_t m_hMtxMultiMatches;
	pthread_rwlock_t m_hRwLock;
#endif

	

public:
	CAligner(void);
	~CAligner(void);

	int
		Align(etPMode PMode,					// processing mode
				UINT32 SampleNthRawRead,		// sample every Nth raw read for processing (1..N)
				etFQMethod Quality,				// quality scoring for fastq sequence files
				bool bSOLiD,					// if true then processing in colorspace
				bool bBisulfite,				// if true then process for bisulfite methylation patterning
				etPEproc PEproc,				// paired reads alignment processing mode
				int PairMinLen,					// accept paired end alignments with observed insert size of at least this (default = 100)
				int PairMaxLen,					// accept paired end alignments with observed insert size of at most this (default = 1000)
				bool bPairStrand,				// accept paired ends if on same strand
				bool bPEcircularised,			// experimental - true if processing for PE spanning circularised fragments
				bool bPEInsertLenDist,			// experimental - true if stats file to include PE insert length distributions for each transcript
				eALStrand AlignStrand,			// align on to watson, crick or both strands of target
				int MinChimericLen,				// minimum chimeric length as a percentage (0 to disable, otherwise 50..99) of probe sequence length: negative if chimeric diagnostics to be reported
				bool bChimericRpt,				// report chimeric trimming detail for individual reads (default is not to report)
				int microInDelLen,				// microInDel length maximum
				int SpliceJunctLen,				// maximum splice junction length when aligning RNAseq reads
				int MinSNPreads,				// must be at least this number of reads covering any loci before processing for SNPs at this loci
				double QValue,					// QValue controlling FDR (Benjamini–Hochberg) SNP prediction
				double SNPNonRefPcnt,			// only process for SNP if more/equal than this percentage number of reads are non-ref at putative SNP loci (defaults to 25.0) 
				int MarkerLen,					// marker sequences of this length
				double MarkerPolyThres,			// maximum allowed marker sequence base polymorphism independent of centroid SNP (default 0.333, range 0.0 to 0.5)
				int PCRartefactWinLen,			// if >= 0 then window size to use when attempting to reduce the number of  PCR differential amplification artefacts (reads stacking to same loci)
				etMLMode MLMode,				// multimatch loci reads processing mode
				int MaxMLmatches,				// accept at most this number of multimatched alignments for any read
				bool bClampMaxMLmatches,		// accept as if MaxMLmatches even if more than this number of MaxMLmatches multimached alignments
				bool bLocateBestMatches,		// align for best matches, not just those 1H better than any other, upto at most MaxMLmatches
				int MaxNs,					    // allow at most this number of indeterminate eBaseNs in read before deeming as nonalignable
				int MinEditDist,				// any matches must have at least this edit distance to the next best match
				int MaxSubs,					// maximum number of substitutions allowed per 100bp of actual read length
				int Trim5,						// trim this number of bases from 5' end of reads when loading the reads
				int Trim3,						// trim this number of bases from 3' end of reads when loading the reads
				int MinAcceptReadLen,					// only accepting reads for alignment if at least this length after any end trimming
				int MaxAcceptReadLen,					// only accepting reads for alignment if no longer than this length after any end trimming
				int MinFlankExacts,				// trim matched reads on 5' and 3' flanks until at least this number of exactly matching bases in flanks
				int PCRPrimerCorrect,			// initially align with MaxSubs+PCRPrimerCorrect subs allowed but then correct substitutions in 5' 12bp until overall sub rate within MaxSubs
				int MaxRptSAMSeqsThres,			// report all SAM chroms or sequences if number of reference chroms <= this limit (defaults to 10000)
				etFMode FMode,					// output format mode
				teSAMFormat SAMFormat,			// if SAM output format then could be SAM, BAM or BAM compressed dependent on the file extension used
				int SitePrefsOfs,				// offset read start sites when processing  octamer preferencing, range -100..100
				int NumThreads,					// number of worker threads to use
				char *pszTrackTitle,			// track title if output format is UCSC BED
				int NumPE1InputFiles,			// number of input PE1 or single ended file specs
				char *pszPE1InputFiles[],		// names of input files (wildcards allowed unless processing paired ends) containing raw reads
				int NumPE2InputFiles,			// number of input PE2 file specs
				char *pszPE2InputFiles[],		// optional raw paired reads are in these files
				char *pszPriorityRegionFile,	// optional high priority BED file contains prioritised region loci
				bool bFiltPriorityRegions,		// true if non-priority alignments to be filtered out 
				char *pszOutFile,				// where to write alignments
				char *pszSNPFile,				// Output SNPs (CSV format) to this file (default is to output file name with '.snp' appended)
				char *pszMarkerFile,			// Output markers to this file
				char *pszSNPCentroidFile,		// Output SNP centorids (CSV format) to this file (default is for no centroid processing)
				char *pszSfxFile,				// target as suffix array
				char *pszStatsFile,				// aligner induced substitutions stats file
				char *pszMultiAlignFile,		// file to contain reads which are aligned to multiple locations
				char *pszNoneAlignFile,			// file to contain reads which were non-alignable
				char *pszSitePrefsFile,			// file to contain aligned reads octamer preferencing
				char *pszLociConstraintsFile,	// loci base constraints file
				char *pszContamFile,			// file containing contaminant sequences
				int	NumIncludeChroms,			// number of chromosome regular expressions to include
				char **ppszIncludeChroms,		// array of include chromosome regular expressions
				int	NumExcludeChroms,			// number of chromosome expressions to exclude
				char **ppszExcludeChroms);		// array of exclude chromosome regular expressions

		int ProcAssignMultiMatches(tsClusterThreadPars *pPars);
		int ProcCoredApprox(tsThreadMatchPars *pPars);
		int ProcLoadReadFiles(tsLoadReadsThreadPars *pPars);

};

