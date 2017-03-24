#pragma once
#include "./commdefs.h"

/*
Fastq scoring schema (from http://en.wikipedia.org/wiki/FASTQ_format )

SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
.................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
..LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
|                         |    |        |                              |                     |
33                        59   64       73                            104                   126
0........................26...31.......40
-5....0........9.............................40
0........9.............................40
3.....9.............................40
0.2......................26...31........41

S - Sanger        Phred+33,  raw reads typically (0, 40)
X - Solexa        Solexa+64, raw reads typically (-5, 40)
I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
(Note: See discussion above).
L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
*/

const unsigned long cDfltStageBuffSize = 0x03ffffff;  // 64M buffer as default
const unsigned long cMaxStageBuffSize = 0x07ffffff;	 // 128M buffer as maximum
const unsigned long cMinStageBuffSize = 0x0fffff;	 // 1M buffer as minimum

const int cNumFastaBlocks = 1;						 // when readahead disk buffering implemented then set to be 2 or more

const unsigned int cMaxGenFastaLineLen = 79;		// limit generated Fasta lines to this length
const unsigned int cMaxFastaDescrLen   = 8192;	    // Fasta descriptor lines can be concatenated..
const unsigned int cMaxFastQSeqLen = cMaxReadLen;	// Allow FastQ sequences to be of this max length

const int cgzAllocInBuffer = 0x1ffffff;				// gz processing input buffer size
const int cgzAllocOutBuffer = 0x1ffffff;			// gz processing output buffer size

#pragma pack(1)
typedef struct TAG_sFastaBlock
	{
	INT64  FileOfs;				// file offset from where fasta block (m_pCurBuffer) last read from file
	INT32 BuffIdx;				// process next char from pBlock[BuffIdx]
	INT32 BuffCnt;				// block currently loaded with this many chars
	INT32 AllocSize;			// block was allocated to buffer at most this many chars in Fasta[]
	UINT8 *pBlock;				// allocated to hold a block of fasta file content
	} tsFastaBlock;
#pragma pack()

class CFasta : public CErrorCodes
{
	int m_hFile;				// opened for write fasta
	gzFile m_gzFile;			// opened for read (could be compressed) fasta or fastq
	char m_szFile[_MAX_PATH];	// to hold fasta file path+name
	UINT64 m_StatFileSize;		// file size as returned by stat() when file initially opened
	bool m_bIsGZ;				// true if processing a gz compressed file
	bool m_bIsFastQ;			// true if processing a fastq file
	bool m_bIscsfasta;			// sequences are in SOLiD csfasta format
	static UINT8 m_SOLiDmap[5][5]; // used for transforming from SOLiD colorspace into basespace
	bool m_bRead;				// TRUE if reading fasta file, FALSE if write to fasta file

	tsFastaBlock *m_pCurFastaBlock;    // buffered fasta block currently being processed
	tsFastaBlock m_FastaBlocks[cNumFastaBlocks];    // allow for at most cNumFastaBlocks buffered fasta file blocks; currently not implemented but in future will allow for readahead of blocks

	bool m_DescrAvail;			// true if NEW descriptor available, reset by ReadDescriptor()
	char m_szDescriptor[cMaxFastaDescrLen+1];	// to hold last descriptor parsed
	char m_szFastqSeq[cMaxFastQSeqLen+1];	    // to hold last FastQ sequence line parsed
	char m_szFastqSeqQ[cMaxFastQSeqLen+1];	    // to hold last FastQ quality line parsed
	int m_FastqSeqLen;
	int m_FastqSeqIdx;
	int	m_FastqSeqQLen;
	int m_CurFastQParseLine;	// current FastQ line being processed

	INT64 m_FileDescrOfs;		// file offset from where last descriptor was parsed
	INT64 m_FileReadDescrOfs;   // offset of last descriptor returned by ReadDescriptor()
	unsigned int m_DescriptorLen;
	unsigned int m_CurLineLen;

	int CheckIsFasta(void);		// checks if file contents are likely to be fasta or fastq format
	int	ParseFastQblockQ(void); // Parses a fastq block (seq identifier + sequence + quality scores)

public:
	CFasta(void);
	~CFasta(void);
	void Cleanup(void);
	int Reset(INT64 FileOfs = 0l);				// reset context to that following an Open() with option to start processing at FileOfs
	int Open(char *pszFile,bool Read = true,unsigned long BufferSize = cDfltStageBuffSize);
	bool IsFastq(void);							// true if opened file is in fastq format
	bool IsSOLiD(void);							// true if opened file is in SOLiD or colorspace format
	UINT64 InitialFileSize(void);				// file size when initially opened for reading

	UINT32										// returns estimated number of sequences in fasta/fastq file
		FastaEstSizes(char *pszFile,			// fasta or fastq file path+name to estimate sizes
			  INT64 *pFileSize = NULL,			// file is this size on disk
			  INT32 *pEstMaxDescrLen = NULL,	// with estimated maximum descriptor length
			  INT32 *pEstMeanDescrLen = NULL,	// estimated mean descriptor length
			  INT32 *pEstMaxSeqLen = NULL,		// and estimated maximum sequence length
			  INT32 *pEstMeanSeqLen = NULL,		// estimated mean sequence length
			  INT32 *pEstScoreSchema = NULL);	// guestimated scoring schema - 0: no scoring, 1: Solexa, 2: Illumina 1.3+, 3: Illumina 1.5+, 4: Illumina 1.8+ or could be Sanger 
	
	int Close(void);							 // close opened fasta file
	int	LocateDescriptor(char *pszPrefixToMatch=NULL,// prefix to match or NULL if match any
						 bool bFromStart=false);	 // true: from start, false: from next descriptor

	int											// number actually read
		ReadSubsequence(etSeqBase *pSeq,		// where to return subsequence (NO TERMINATING eBaseEOS)
					 unsigned int MaxLen,		// reads upto MaxLen sequence bases from file
					 unsigned int SeqOfs=0,		// relative offset in sequence at which to read
					 bool bFromStart=false);	// true: from fasta file start, false: from next sequence

	int ReadSequence(void *pRetSeq = NULL,		// where to return sequence, can be NULL if only interested in the sequence length
					 int Max2Ret = 0,			// max to return, ignored if pRetSeq is NULL
					 bool bSeqBase = true,		// if false then return ascii, if true then return as etSeqBase
					 bool RptMskUpperCase = false);	// default is false, UCSC softmasked use lowercase when repeat masked


	int ReadQValues(char *pRetQValues,	// where to return quality values
					 int Max2Ret);			// max to return

	int ReadDescriptor(char *pszDescriptor,int MaxLen); // copies last descriptor processed into pszDescriptor and returns copied length
	INT64 GetDescrFileOfs(void);				// returns file offset at which descriptor returned by ReadDescriptor() was parsed from

	static int Ascii2Sense(char *pAscii,		  // expected to be '\0' terminated, or SeqLen long
					int MaxSeqLen,				  // maximal sized sequence that pSeq can hold
					etSeqBase *pSeq,			  // where to return sequence
					bool RptMskUpperCase=false);  // true if bases are softmasked as uppercase instead of default lowercase

	int Write(char Symbol);							// writes sequence char to fasta file
	int Write(char *pSymbols,unsigned int Cnt);		// writes Cnt sequence chars to fasta file
	int WriteDescriptor(char *pszDescriptor);		// terminates any current sequence and writes fasta file descriptor starting new sequence
};

