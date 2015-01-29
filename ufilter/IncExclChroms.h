#pragma once

const int cAllocIncExclChroms = 100;		// allocate for this many chromosomes in an allocation block
const int cMaxIncludeChroms = 20;		// max number of include chromosomes regular expressions
const int cMaxExcludeChroms = 20;		// max number of exclude chromosomes regular expressions

typedef struct TAG_sIncExclChrom {
	bool bInc;								// true if chromsome to be included, false if chromosome to be excluded
	int ChromID;							// used to uniquely identify this chrom
	char szChrom[cMaxDatasetSpeciesChrom];	// chromosome name
} tsIncExclChrom;


class CIncExclChroms {
	tsIncExclChrom *m_pLastmatch;				// pts to last added or matched chromosome
#ifdef _WIN32
	Regexp *m_IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	Regexp *m_ExcludeChromsRE[cMaxExcludeChroms];
#else
	regex_t m_IncludeChromsRE[cMaxIncludeChroms];	// compiled regular expressions
	regex_t m_ExcludeChromsRE[cMaxExcludeChroms];
#endif

	int m_NumIncludeChroms;			// number of chromosome regular expressions to include
	int	m_NumExcludeChroms;			// number of chromosome expressions to exclude

	int m_AllocdIncExclChroms;					// number allocated
	int m_IncExclChroms;						// number in use
	tsIncExclChrom *m_pIncExclChroms;			// allocated array
	
	tsIncExclChrom *LocateChrom(char *pszChrom);  // returns previously cached chrom processing state
	int AddChrom(char *pszChrom,bool bProcChrom); // caches chrom processing state

public:
	CIncExclChroms();
	CIncExclChroms(int	NumIncludeChroms,	// number of chromosome regular expressions to include
		char **ppszIncludeChroms,	// array of include chromosome regular expressions
		int	NumExcludeChroms,		// number of chromosome expressions to exclude
		char **ppszExcludeChroms);

	~CIncExclChroms();

	void Reset(void);
	int InitChromExclusion(int	NumIncludeChroms,	// number of chromosome regular expressions to include
		char **ppszIncludeChroms,	// array of include chromosome regular expressions
		int	NumExcludeChroms,		// number of chromosome expressions to exclude
		char **ppszExcludeChroms);	// array of exclude chromosome regular expressions
	int IncludeThisChrom(char *pszChrom);
	char *LocateChrom(int ChromID);	// returns chrom corresponding to specified ChromID

};

