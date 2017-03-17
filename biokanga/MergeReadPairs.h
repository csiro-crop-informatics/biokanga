
#pragma once

const int cMaxBarCodeLen = 30;	// if barcode processing then barcodes can be at most this length
const int cMaxLenMIDsName = 40;  // individual well identifiers or MIDs names can be at most this length
const int cMaxNumBarcodes = 1000;	// able to process at most this many different barcodes
const int cMinBCReadSeqLen = 50;   // only read sequences of at least this length will be processed for barcodes

// Merge read pairs where there is an expected 3' read overlap onto the 5' read

// output merged file format
typedef enum TAG_ePMode {
	ePMdefault = 0,				// Standard processing
	ePMcombined,				// combine overlaps and orphan reads into same file
	ePMseparate,				// overlaps and orphans into separate files
	ePMAmplicon,				// amplicon processing with 5' and 3' barcodes to identify originating plate well, overlap merging
	ePMAmpliconNoMerge,			// amplicon processing with 5' and 3' barcodes to identify originating plate well, no overlap merging
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

// output merged file format
typedef enum TAG_eOFormat {
	eOFauto = 0,				// default is to merge and output in same format as inputs
	eOFfasta,					// merge and output as fasta
	eOFfastq,					// merge and output as fastq
	eOFplaceholder				// used to set the enumeration range
	} etOFormat;


const int cMinOverlapLen = 1;		// minimium allowed used specified number of overlap bases when merging PE sequences
const int cMaxOverlapPercSubs = 10;	// maximum allowed substitutions as a percentage of the overlap between merged PE sequences

const int cAllocOutBuffLen = 0x0fffff; // 1M output buffering when writing to file

const char cPhredLowScore = 'A';	// generated Phred score used when there was a merge induced substitution in overlay and input sequences were non-fastq
const char cPhredOrphanScore = 'H';	// generated Phred score used if orphan reads and input sequences were non-fastq
const char cPhredHiScore = 'J';		// generated Phred score used when there was no merge substitution in overlay and input sequences were non-fastq

typedef enum TAG_eProcPhase {
	ePPUninit,					// uninitialised
	ePPReset,					// reset ready for file processing
	ePPRdyOpen,					// ready for input files to be opened 
	ePRdyProc					// files opened, ready for processing
} etProcPhase;

#pragma pack(1)

typedef struct TAG_sMIDsBarcode  // MIDs (Multiplex Identifiers) barcode
	{
	int MIDsID;				// monotonically increasing MIDs identifier
	char szMIDsName[cMaxLenMIDsName + 1];	// MIDs name
	int PE1BarcodeLen;				// SE or PE1 MIDs barcode sequence length
	etSeqBase PE1Barcode[cMaxBarCodeLen];	// SE or PE1 MIDs barcode sequence
	int PE2BarcodeLen;				// if PE MIDs then PE2 MIDs barcode sequence length
	etSeqBase PE2Barcode[cMaxBarCodeLen];	// if PE MIDs then PE2 MIDs barcode sequence
	} tsMIDsBarcode;

typedef struct TAG_sBarcode {
	UINT8 ColRow;		// 0 if barcode at 5' end (column) and 1 if 3' end (row) of amplicon
	UINT8 Psn;		    // row or column position 1..N
	UINT8 Len;			// barcode length
	UINT32 BarCode;		// packed barcode bases	
	UINT32 RevCplBarCode;  // revcpl packed barcode bases	
	} tsBarcode;

typedef struct TAG_sAmpliconWellFile
	{
	char szOutFile[_MAX_PATH];	// writing to this file
	int hOutFile;			    // file handle for opened szOutFile
	int CurBuffLen;				// currently pOutBuffer holds this many chars ready to be written to file
	int AllocdOutBuff;			// pOutBuffer can hold at most this many chars
	char *pOutBuffer;			// allocated to buffer output
	} tsAmpliconWellFile;

typedef struct TAG_sAmpliconWell {
	int WellID;					// identifies well 1..96, or if 454 reads then the mapped to BarcodeID
	int NumASequences;			// number of sequences attributed to this well
	tsAmpliconWellFile WellFile[2]; // if SE then output to 1 file, if PE then each end written to separate files
} tsAmpliconWell;

#pragma pack()

class CMergeReadPairs
{
	etProcPhase m_ProcPhase;	// current processing phase
	etPMode m_PMode;			// processing mode 
	bool m_bAmpliconNoMerge;	// true if amplicon overlaps not to be merged and output as SE (will be output as PE end reads)
	etOFormat m_OFormat;		// output file format
	bool m_bIsFastq;			// if basespace then input reads are fastq
	bool m_bAppendOut;			// if true then append if output files exist, otherwise trunctate existing output files 
	int m_MinOverlap;			// reads must overlap by at least this many bases
    int m_MaxOverlapPropSubs;	// and overlap can have at most this proportion of substitutions

	CFasta m_PE5Fasta;			// 5' (PE1) reads being processed
	CFasta m_PE3Fasta;          // 3' (PE2) reads being processed

	int m_hOutMerged;			// file handle for output merged reads if not processing barcoded reads
	int m_hOut5Unmerged;		// file handle for output 5' unmerged reads
	int m_hOut3Unmerged;		// file handle for output 3' unmerged reads
	
	char m_szMergeOutFile[_MAX_PATH];	// use this file name as the merged output file name, or prefix filename if processing barcoded reads
	char m_szIn5ReadsFile[_MAX_PATH];	// 5' reads currently being processed are in this file
	char m_szIn3ReadsFile[_MAX_PATH];	// 3' reads currently being processed are in this file
	char m_sz5UnmergedOutFile[_MAX_PATH];	// write 5' unmerged reads to this file
	char m_sz3UnmergedOutFile[_MAX_PATH];   // write 3' unmerged reads to this file

	int m_CurMSeqLen;			// m_pszMSeqs currently holds this many chars
	int m_AllocdMSeqs;			// m_pszMSeqs can hold at most this many chars
	char *m_pszMSeqs;			// will be allocated for holding merged sequences ready for writing to file if not processing barcoded reads
	int m_CurUnmergedP1Seqs;	// m_pszUnmergedP1Seqs currently holds this many chars
	int m_AllocdUnmergedP1Seqs;	// m_pszUnmergedP1Seqs can hold at most this many chars
	char *m_pszUnmergedP1Seqs;  // will be allocated for holding unmerged P1 sequences ready for writing to file
	int m_CurUnmergedP2Seqs;	// m_pszUnmergedP2Seqs currently holds this many chars
	int m_AllocdUnmergedP2Seqs;	// m_pszUnmergedP2Seqs can hold at most this many chars
	char *m_pszUnmergedP2Seqs;  // will be allocated for holding unmerged P2 sequences ready for writing to file

	int m_MaxBarcode5Len;		// maximun length of any 5' barcode - max 16
	int m_MaxBarcode3Len;		// maximun length of any 3' barcode - max 16
	int m_NumBarcodes;			// number of barcodes in m_pBarcodes[]
	tsBarcode *m_pBarcodes;		// well barcoding
	int m_NumWells;				// processing is for this number of wells
	tsAmpliconWell m_WellFiles[cMaxNumBarcodes];    // to hold merged well sequence file output buffers and respective file handles

	int m_454Leadin;				// start checking for barcode after this many leadin bases (default 4)
	int m_Max454BarcodeSubs;		// allowing at most this many subs (default 1) when matching MIDs barcodes
	int m_NumMIDsBarcodes;			// number of MIDs barcodes actually loaded
	tsMIDsBarcode m_MIDsBarcodes[cMaxNumBarcodes];   // holds MIDs names and associated barcodes

	int											// returned best uniquely matching MIDsID or 0 if unable to match any MIDs barcode allowing at most MaxAllowSubs substitutions
		MapSEMIDsBarcode(int SeqLen,			// number of bases in sequence putatively containing MIPs barcode 
						 int Ofs,				// look for matching MIDs barcode starting at this pSeq 5' offset
						 int MaxAllowSubs,		// max allowed substitutions
						 int *pReqSubs,			// best match required this number of substitutions
						 etSeqBase *pSeq);		// sequence

	int											// returned best uniquely matching MIDsID or 0 if unable to match any MIDs barcode allowing at most MaxAllowSubs substitutions
		MapPEMIDsBarcode(int SeqLen,			// minimum length of either pPE1Seq or pPE2Seq  
						 int Ofs,				// look for matching MIDs barcode starting at this pPE1Seq or pPE2Seq 5' offset
						 int MaxAllowSubs,		// max allowed substitutions for either PE1 or PE2 (not total over both PE1 and PE2)
						 int *pReqSubs,			// best match required this total number of substitutions over both PE1 and PE2 barcodes
						 etSeqBase *pPE1Seq,	// map barcodes at 5' of this sequence
						 etSeqBase *pPE2Seq);	// and barcodes at 5'of this sequence to the MIDs

	int							// returned well number (1..96) or 0 if unable to identify well from the barcodes
		MapSEBarcodesToWell(int SeqLen,		// num bases in SE pSeq
				etSeqBase *pSeq);			// map barcodes at 5' and 3' end of this SE sequence to the well

	int							// returned well number (1..96) or 0 if unable to identify well from the barcodes
		MapPEBarcodesToWell(int SeqLen,		//minimum length of either pPE1Seq or pPE2Seq
											 etSeqBase *pPE1Seq,// map barcodes at 5' of this sequence
											 etSeqBase *pPE2Seq); // and barcodes at 5'of this sequence to the well

	char *RemoveQuotes(char *pszRawText);

public:
	CMergeReadPairs(void);
	~CMergeReadPairs(void);
	void Reset(bool bSync = true);   // if bDync true then opened output files will be sync'd (commited) to disk before closing

	int				// return number of MIDs initialised
		InitMIDs(bool bPE,		// true if processing for PE MIDs else if false then processing SE MIDs
				 char *pszMIDsBarcodeFile); // initialise with SE MIDs barcodes and associated name identifiers from a MIDs CSV file  

	int				// return number of wells initialised
		InitDfltWells(bool bNoMerge = false); // initialise with default well barcodes and well identifiers, if bNoMerge then do not merge PE reads and report PE1/PE2 instead of merged SE 

	int									// returns number of sequences merged and written to output files
			MergeOverlaps(etPMode PMode,		// processing mode 
			  etOFormat OFormat,				// output file format
			  int MinOverlap,					// reads must overlap by at least this many bases
              int MaxOverlapPropSubs,			// and overlap can have at most this proportion of substitutions, if > 0 then floor of 1 sub allowed
			  int NumInPE5Files,				// number of input single ended or 5' end if paired end reads file
			  char **pszInPE5Files,				// input single ended or 5' end if paired end reads files
			  int NumInPE3Files,				// number of input input 3' end if paired end reads files
			  char **pszInPE3Files,				// input 3' end if paired end reads files
			  char *pszMergeOutFile,			// write merged overlaping reads to this file
  			  int StartNum = 1,					// use this initial starting sequence identifier
  			  bool bAppendOut = false);			// if true then append if output files exist, otherwise trunctate existing output files

	int													// returns number of sequences merged and written to output file
		ProcOverlapPairs(int StartNum);		// initial starting sequence identifier

	int 
		OpenFiles(void);					// open files (as set by MergeOverlaps) for processing
};



