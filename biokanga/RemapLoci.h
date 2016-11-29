#pragma once

const int cAllocOutBuff = 0x100000;		// output buffering allocation size if BED or CSV format output

class CRemapLoci
{
	CSAMfile *m_pInBAMfile;				// if SAM/BAM input
	CSAMfile *m_pOutBAMfile;			// if SAM/BAM output

	CBEDfile *m_pMappingBED;			// BED containing remapping locii


	int m_hOutFile;						// output file handle for BED 
	int m_OutBuffIdx;					// current index into m_pszOutBuff at which to next write output formated for BED 
	int m_AllocOutBuff;				    // output buffer allocated to hold this many chars
	char *m_pszOutBuff;					// allocated output buffer

	char *TrimWhitespace(char *pTxt);	// trim whitespace
	bool								// true if file to be generated compressed with gzopen/gzwrite/gzclose
		FileReqWriteCompr(char *pszFile); // If last 3 chars of file name is ".gz" then this file is assumed to require compression

public:
	CRemapLoci();
	~CRemapLoci();

	void Reset(void);

	int
	RemapLocii(int PMode,					// processing mode
					 int FType,				// alignment file type
					char *pszInAlignFile,	// alignment file with loci to be remapped
					char *pszInBEDFile,     // BED file containing loci remapping
					char *pszRemappedFile);	// write remapped alignments to this file

	int
	RemapBEDLocii(char *pszInAlignFile,		// BED alignment file with loci to be remapped
				char *pszRemappedFile);		// write remapped alignments to this file

	int
	RemapSAMLocii(char *pszInAlignFile,		// SAM or BAM alignment file with loci to be remapped
				char *pszRemappedFile);		// write remapped alignments to this file

};

