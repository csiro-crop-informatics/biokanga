#pragma once
#include "commdefs.h"

const int cFileClassifyBuffLen = 8196;	
const int cMinFileClassifyLen = 100;			// aribarily - don't bother if less than this length

typedef enum eClassifyFileType {
	eCFTCSV = 0,		// file has been classified as being CSV
	eCFTBED,			// file has been classified as being BED
	eCFTSAM,			// file has been classified as being SAM
	eCFTopenerr,		// unable to open file for reading
	eCFTlenerr,			// file length is insufficent to classify type
	eCFTunknown,		// unable to reliably classify
	} etClassifyFileType;

class CUtility
{

public:
	CUtility(void){};
	~CUtility(void){};


	// ClassifyFileType
// Attempt to classify the alignment file as one of CSV, BED or SAM from it's initial 8k char contents
// Currently processes CSV, BED and SAM format file types
// Assumes must be SAM if initial lines have at least one prefixed by a '@' followed by a 2 letter record type code 
//

	static etClassifyFileType ClassifyFileType(char *pszFileName);	// attempt to classify file as being either CSV, BED or SAM

	static void splitpath(char *pszFullPath, char *pszDir,
                 char *pFname);

	static int TrimQuotes(char *pszTxt);
	static char *TrimWhitespc(char *pszText);
	static char *TrimWhitespcExtd(char *pszText);
	static char *TrimQuotedWhitespcExtd(char *pszText);
	static char *ReduceWhitespace(char *pszText);

	static UINT16			// generated 16bit hash over the lowercased chromosome name; hashes as 0 if pszName == null or is empty
		GenHash16(char *pszName);	// name will be lowercased whilst hashing

	static int						// generated 24bit hash over the lowercased chromosome name; hashes as 0 if pszName == null or is empty
		GenHash24(char *pszName);	// name will be lowercased whilst hashing


	bool						// true if requested number of bytes was written to pOutFile
		SafeWrite(void *pOutFile,		// either gzFile or file handle
					void *pBuff,		// pts to buffer containing data to write
					size_t Max2Write,	// write this many bytes from pBuff
					bool bgzOut=false);		// true if pOutFile pts to a gzFile, false if pOutFile is file handle

	static bool				// true if requested number of bytes was written to hOutFile 
		SafeWrite(int hOutFile,void *pBuff,size_t Max2Write);

	static bool	// true if requested number of bytes have been compressed to pgzFile
		SafeWrite_gz(gzFile pgzFile,void *pBuff,size_t Max2Write);

	static int arg_parsefromfile(int argc,		// cnt of args ptd to by *argv[]
							char *argv[],		// pts to array of args
							char **retargv[]);	// returned array of args

	static char *			// returns current and maximum resource limits ( getrlimit ), returns NULL in WIndows
		ReportResourceLimits(void);

	static int
		GetNumSubseqs(int AlignLen,		// alignment length incl InDels
			   int NumSeqs,				// number of sequences
			   etSeqBase *pSeqs[1]);		// pts to array of ptrs to the sequences to process

	static int						// returned subsequence length
	GetFilteredSubSeq(int AlignLen, // remaining alignment length incl InDels
		  int StartIdx,			// initial starting index into pRef[]/pRel[]
		  int *pFirstIdx,		// returned index into pRef[]/pRel[] of first base in returned subsequence
		  int *pLastIdx,		// returned index into pRef[]/pRel[] of last base in returned subsequence
		  etSeqBase *pRef,		// reference sequence
		  etSeqBase *pRel,		// relative sequence
  		  int ReqMinLen=1,		// subsequence must be of at least this length
		  int ReqMinIdent=0,		// and have at least this overall identity (0..100%)
		  bool bStartEndIdent=false);  // and start/end on base that is identical

	static int					// returned number of subsequences
		   GetNumFilteredSubseqs(int AlignLen,		// alignment length incl InDels
		   etSeqBase *pRef,		// reference sequence
		   etSeqBase *pRel,		// relative sequence
   		   int ReqMinLen=1,		// subsequence must be of at least this length
		   int ReqMinIdent=0,		// and have at least this overall identity (0..100%)
		   bool bStartEndIdent=false);  // and start/end on base that is identical

	// ChkTargDepend
	// Check list of source files against a target file
	// If target file does not exist, or any source file is newer than target then returns 1
	// If target file is newer or same age as all the source files then returns 0
	// If any source file is specified but does not exist then returns -1 and copies file name into pszInaccessible
	static int ChkTargDepend(char *pszInaccessible,	// where to copy name of the first detected inaccessible source file (if any)
						int MaxLen,				// max length to copy
						char *pszTarg, ...);    // var number of source files, use NULL or empty (pszSource[0] == '\0') to terminate

	// Chk2TargDepend
	// Check list of source files against two target files
	// If target file does not exist, or any source file is newer than target then returns 1
	// If target file is newer or same age as all the source files then returns 0
	// If any source file is specified but does not exist then returns -1 and copies file name into pszInaccessible
	static int Chk2TargDepend(char *pszInaccessible,	// where to copy name of the first detected inaccessible source file (if any)
						int MaxLen,				// max length to copy
						char *pszTarg1,char *pszTarg2, ...); // var number of source files, use NULL or empty (pszSource[0] == '\0') to terminate

	static void SleepMillisecs(UINT32 milliseconds); // cross-platform sleep function


};
