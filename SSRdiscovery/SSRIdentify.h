#pragma once

const size_t cMaxAllocBuffChunk = 0x0ffffff;		// allocate for input sequences buffering in these sized chunks
const int cMaxAllocRptSSRs = 0x07fffff;				// SSR reporting buffer size

const int cMinRepElLen = 1;				// minimum element K-mer length
const int cDfltMinRepElLen = 2;			// default minimum element K-mer length
const int cDfltMaxRepElLen = 5;			// default maximum element k-mer length
const int cMaxRepElLen = 20;			// maximum element k-mer length

const int cMinTandemRpts = 2;			// min number of tandem repeats
const int cDfltMinTandemRpts = 5;		// default min tandem repeats
const int cDfltMaxTandemRpts = 10;		// defualt max tandem repeats
const int cMaxTandemRpts = 50;			// max number of tandem repeats

// reporting of SSR loci is in one of the following formats
typedef enum eRptSSRsFromat {
	eRFCsv = 0,			// report as CSV
	eRFBed,				// BED file format
	eRFSam				// SAM file format
	} teRptSSRsFromat;

#pragma pack(1)
typedef struct TAG_sKMerDist {
	UINT32 Cnt;							// number of occurances of SSR with this repeating K-mer element
	UINT32 TandemRpts[1];				// number of tandem repeats (extended to user specified max repeats)
	} tsKMerDist;
#pragma pack()

class CSSRIdentify
{
	CStopWatch m_CurTime;			// used for progress messaging
	teRptSSRsFromat m_RptSSRsFormat;	// format in which to report SSRs
	int m_hOutFile;				// output processing results to this file
	int m_hOutKMerFreqFile;		// output repeating element KMer freqs to this file 

	int m_IdxRptSSRs;			// current index into m_pszRptSSRsBuff for buffered SSRs reporting
	char *m_pszRptSSRsBuff;		// allocated to buffer SSRs reporting

	size_t m_SeqBuffLen;		// number of bases currently buffered in m_pSeqBuff
	size_t m_AllocdSeqBuffMem;  // size of memory currently allocated to m_pSeqBuff
	UINT8 *m_pSeqBuff;			// buffers sequences as read from file

	int m_MinRepElLen;			// identify repeating elements of this minimum length
	int m_MaxRepElLen;			// ranging upto this maximum length
	int m_MinTandemRpts;		// minimum number of tandem repeats
	int m_MaxTandemRpts;		// maximum number of repeats

	UINT32 m_TotNumAcceptedSSRs;	// total number of SSRs accepted
	UINT32 m_TotNumExcessiveTandemSSRs;    // total number of putative SSRs not accepted because m_MaxTandemRpts

	UINT32 m_TotNumAcceptedKmerSSRs[cMaxRepElLen+1];

	int m_KMerFreqLen;				// if > 0 then counting this length KMer
	int m_NumKMers;					// number of KMers
	int m_SizeOfKMerDist;			// size of each tsKMerDist allowing for max number of tandem repeats requested by user
	size_t m_AllocdKMerFreqMem;		// size of memory currently allocated to m_pKMerFreq
	tsKMerDist *m_pKMerDist;		// allocated to hold KMer sequence instance counts; only used if single KMer element length being processed

	int ProcessBioseqFile(char *pszFile);	// load and process a bioseq file for SSRs
	int ProcessFastaFile(char *pszFile);	// load and process a multifasta file for SSRs

	int Report(int RepElLen,			// identified SSR contains repeat elements of this length
			   int NumTandemEls,		// each repeat element is repeated this many times
			   INT64 SSRStartOfs,		// repeat element starts at this offset within the targeted sequence
				char *pszDescr,			// descriptor for the targeted sequence
				char *pszInFile,		// sequence parsed from this file
			   INT64 TargSeqLen,		// targeted sequence is this length
			   etSeqBase *pTargSeq);	// targeted sequence within which the SSR has been located

	int
		ReportCSV(int RepElLen,			// identified SSR contains repeat elements of this length
			int NumTandemEls,			// each repeat element is repeated this many times
			INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
			char *pszDescr,				// descriptor for the targeted sequence
			char *pszInFile,			// sequence parsed from this file
			INT64 TargSeqLen,			// targeted sequence is this length
			etSeqBase *pTargSeq);		// targeted sequence within which the SSR has been located

	int
		ReportBED(int RepElLen,			// identified SSR contains repeat elements of this length
			int NumTandemEls,			// each repeat element is repeated this many times
			INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
			char *pszDescr,				// descriptor for the targeted sequence
			char *pszInFile,			// sequence parsed from this file
			INT64 TargSeqLen,			// targeted sequence is this length
			etSeqBase *pTargSeq);		// targeted sequence within which the SSR has been located

	int
		ReportSAM(int RepElLen,			// identified SSR contains repeat elements of this length
			int NumTandemEls,			// each repeat element is repeated this many times
			INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
			char *pszDescr,				// descriptor for the targeted sequence
			char *pszInFile,			// sequence parsed from this file
			INT64 TargSeqLen,			// targeted sequence is this length
			etSeqBase *pTargSeq);		// targeted sequence within which the SSR has been located

	int	ReportKMers(char *pszKMerFreqFile);	// report SSR repeating element K-mer frequencies to this file

	int	ReportProgress(bool bForce = false);	// let user know that there is processing activity, normally progress reportde evry 60 sec unless bForce set true

	int
		IdentifySSRs(
			 char *pszDescr,			// descriptor for the targeted sequence
			 char *pszInFile,			// sequence parsed from this file
			 INT64 TargSeqLen,			// sequence length of targeted sequence within which to search for SSRs
			 etSeqBase *pTargSeq);		// identify SSRs in this targeted sequence

	etSeqBase *AllocSeqBuff(size_t SeqLen);	// allocate for at least this sequence length

	int CntKMer(int KMerLen,			// count this KMer
				int Rpts,				// tandem repeat counts
				etSeqBase *pSeq);		// sequence 

public:
	CSSRIdentify(void);
	~CSSRIdentify(void);

	void Init(void);						// initialisation
	void Reset(void);					// resets state back to that imediately following initialisation
	int
		Process(int PMode,			// procesisng mode - currently unused..
			teRptSSRsFromat RptSSRsFormat,	// report SSRs in this file format
			int MinRepElLen,	// identify repeating elements of this minimum length
			int MaxRepElLen,		// ranging upto this maximum length
			int MinTandemRpts,		// minimum number of tandem repeats
			int MaxTandemRpts,		// maximum number of repeats
			int NumInFileSpecs,		// number of input, could be wildcarded, file specs
			char *pszInFile[],		// files to be processed
			char *pszKMerFreqFile,	// optional, output element KMer freq to this file
			char *pszOutFile);		// SSRs to this file

};

