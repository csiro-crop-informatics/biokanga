#pragma once
#include "./commdefs.h"

const int cBSFMaxErrMsgs	= 10;		// number of error messages retained
const int cBSFMaxErrMsgLen	= 512;		// max length of any error message
const int cBSFMaxErrMsgSrc	= 25;		// message source max length
const int cBSFMaxErrMsgTime = 50;		// timestamp + message source max length

// global error codes
// Success is represented by positive values >= 0
// Errors are always < 0



typedef enum eBSFrsltCodes {
	eBSFSuccess = 0,		// success!
	eBSFFastaDescr = INT_MAX,// fasta descriptor line avail
	eBSFerrParams=-100,		// parameter error
	eBSFerrFileClosed,		// must open file before calling this method
	eBSFerrOfs,				// offset specified past end of data
	eBSFerrCvrtType,		// unable to convert to specified data type
	eBSFerrEntryCreate,		// Can't seal entry if none started
	eBSFerrMem,				// unable to alloc memory
	eBSFerrNotBioseq,		// file exists but not a bioseq file
	eBSFerrNotFasta,		// file exists but not a fasta
	eBSFerrFileType,		// file is a bioseq file, but does not contain specified type
	eBSFerrFileDPT,			// file is a bioseq file, but does not contain Data Points
	eBSFerrOpnFile,			// unable to open file
	eBSFerrCreateFile,		// unable to create file
	eBSFerrExecFile,		// unable to execute external process
	eBSFerrClosed,			// file is closed
	eBSFerrFileVer,			// file version error
	eBSFerrFileAccess,		// file access (seek/read/write) failed
	eBSFerrRead,			// file not opened for read
	eBSFerrWrite,			// file not opened for write
	eBSFerrTypeChange,		// data type changed
	eBSFerrUserType,		// can't currently handle user defined dataset types
	eBSFerrConfParam,		// illegal structure parameter specified
	eBSFerrFastaDescr,		// can't locate fasta descriptor
	eBSFerrFastaSymb,		// illegal fasta sequence symbol
	eBSFerrFastqChr,		// illegal fastq character
	eBSFerrFastqSeqID,		// problem in fastq sequence identifier
	eBSFerrFastqSeq,		// problem in fastq sequence
	eBSFerrFastqDescr,		// problem in fastq descriptor
	eBSFerrFastqQScores,	// problem in fastq quality scores

	eBSFerrLocField,		// Expected field not present
	eBSFerrStructParm,		// Error processing structural parameters file
	eBSFerrStructStep,		// Unable to determine conformation for step
	eBSFerrCentroidParam,	// Error processing centroids file
	eBSFerrNumSrcFiles,		// too many source files
	eBSFerrNoSrcFiles,		// no source file
	eBSFerrNumAlignBlks,	// can't add any more alignment blocks
	eBSFerrAlignBlk,		// can't locate alignment block or block not started
	eBSFerrObj,				// unable to instantiate internal object
	eBSFerrRefDataset,		// reference dataset error
	eBSFerrDataPtType,		// unsupported data point type
	eBSFerrNoFeatures,		// feature file is empty - no features
	eBSFerrNoEntries,		// bioseq file is empty - no entries
	eBSFerrMaxDirEls,		// hit max datasegs
	eBSFerrMaxChroms,		// hit max chromosomes
	eBSFerrMaxDatasets,		// hit max datasets limit
	eBSFerrMaxFeatures,		// hit max features
	eBSFerrMaxEntries,		// hit max entries
	eBSFerrFeature,			// can't locate feature
	eBSFerrExon,			// can't locate exon
	eBSFerrIntron,			// can't locate intron
	eBSFerrEntry,			// can't locate entry
	eBSFerrDataset,			// can't locate dataset
	eBSFerrChrom,			// can't locate chromosome
	eBSFerrSpecies,			// can't locate species
	eBSFerrParse,			// error whilst parsing file - unexpected format
	eBSFerrGene,			// BEDfile file does not contain gene details - utr,exons,introns etc
	eBSFerrProbMatrices,	// Error processing transitional probabilities matrices file
	eBSFErrBase,			// unrecognised base
	eBSFerrFileOpened,		// File must be closed before calling this method
	eBSFerrFileName,		// Unable to parse filename for ref+rel species names
	eBSFerrFieldCnt,		// CSV number of fields parsed not number expected
	eBSFerrRowCnt,			// CSV number of rows parsed not number expected
	eBSFerrFieldID,			// CSV field identifier outside of range
	eBSFerrNumRange,		// numeric conversion range error
	eBSFerrQuoteErrs,		// additional chars followed end of quoted text before field terminator
	eBSFerrQuoteIncomplete,	// quoted CSV field was not terminated by another quote
	eBSFerrRegion,			// requested region not supported
	eBSFerrDupGOTerm,		// duplicate GO:Term
	eBSFerrGOID,			// unable to locate requested GO:Ident identifier
	eBSFerrGOTagVal,		// unable to locate requested GO: tag value
	eBSErrSession,			// session requested not active
	cBSFSyncObjErr,			// unable to synchronise access to object

	cBSFSocketErr,			// TCP socket level error
	cBSFNWSProtErr,			// Services protocol error


	eBSFerrInternal = -1	// internal processing error, inconsistency detected
} teBSFrsltCodes; 

class CErrorCodes
{
	int m_NumMsgs;					// number of error messages
	int m_OldestMsg;				// index of oldest error message (0..cBSFMaxErrMsgs-1)
	char m_szErrMsgs[cBSFMaxErrMsgs][cBSFMaxErrMsgLen + cBSFMaxErrMsgTime];
public:
	CErrorCodes(void);
	~CErrorCodes(void);
	static const char *ErrText(teBSFrsltCodes ErrCode);		// returns error text for specified rslt code
	int AddErrMsg(const char *pszSource,	// used to identify message source
				  const char *pszFormat,...); // message text (sprintf format) plus any parameters
	char *GetErrMsg(void);			// pops and returns oldest error message
	int NumErrMsgs(void);			// returns number of error messages
	void ClearErrs(void);			// clears all error messages

};
