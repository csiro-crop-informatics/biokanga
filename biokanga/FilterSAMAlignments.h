#pragma once

const int cMaxExcludeChroms = 20;		// allow upto this many regexpr for specifying chroms to exclude
const int cMaxIncludeChroms = 20;		// allow upto this many regexpr for specifying chroms to include

class CFilterSAMAlignments
{
	CSAMfile *m_pInBAMfile;				// SAM/BAM input file
	CSAMfile *m_pOutBAMfile;			// SAM/BAM output file

	int m_NumIncludeChroms;			// number of chromosomes explicitly defined to be included
	char **m_ppszIncludeChroms;		// ptr to array of reg expressions defining chroms to include
	int m_NumExcludeChroms;			// number of chromosomes explicitly defined to be excluded
	char **m_ppszExcludeChroms;		// ptr to array of reg expressions defining chroms to include
	#ifdef _WIN32
	Regexp *m_IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	Regexp *m_ExcludeChromsRE[cMaxExcludeChroms];
	#else
	regex_t m_IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	regex_t m_ExcludeChromsRE[cMaxExcludeChroms];
	#endif

	char m_szFiltChrom[_MAX_PATH];	// used to cache last chrom processed	
	bool m_bFiltChrom;				// and it's filtered status


	int
	SetChromFilters(int NumIncludeChroms,		 // number of chromosome regular expressions to include
						char **ppszIncludeChroms,	 // array of include chromosome regular expressions
						int NumExcludeChroms,		 // number of chromosome expressions to exclude
						char **ppszExcludeChroms);	 // array of exclude chromosome regular expressions

// ExcludeThisChrom
// Returns true if pszChrom is to be excluded from processing
	bool	ExcludeThisChrom(char *pszChrom);

	char *TrimWhitespace(char *pTxt);	// trim whitespace
	bool								// true if file to be generated compressed with gzopen/gzwrite/gzclose
		FileReqWriteCompr(char *pszFile); // If last 3 chars of file name is ".gz" then this file is assumed to require compression

public:
	CFilterSAMAlignments();
	~CFilterSAMAlignments();

	void Reset(void);
	void Init(void);

	int								// number of alignments which were retained and written to output file after filtering was applied
		FilterSAMbyChrom(int NumIncludeChroms,		// number of retained chromosomes regular expressions
		char **ppszIncludeChroms,	// array of include chromosome regular expressions
		int NumExcludeChroms,		// number of chromosome expressions to explicitly exclude
		char **ppszExcludeChroms,	// array of exclude chromosome regular expressions
		char *pszInFile,			// input file containing alignments to be filtered
		char *pszOutFile);			// write filtered alignments to this output file
};

