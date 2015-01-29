#pragma once

const int cFiltLociAllocChunk = 10000;	// allocate for loci in this many instances

typedef struct TAG_sLociChrom {
	int ChromID;							// unique identifier for this chromosome
	UINT16 Hash;							// hash over chromosome name
	char szChrom[cMaxDatasetSpeciesChrom];	// chromosome
} tsLociChrom;

typedef struct TAG_sFilterLoci {
	int ChromID;				
	int StartLoci;
	int EndLoci;
} tsFilterLoci;

class CFilterLoci  : public CCSVFile
{
	int m_NumFilterLocii;			// number of loci in m_pFilterLocii
	int m_AllocdFilterLocii;		// allocd FilterLocii
	tsFilterLoci *m_pFilterLocii;	// array of FilterLocii

	int m_CachLociChromID;			// last located chromosome
	int m_NumLociChroms;			// number of chromosomes
	int m_AllocdLociChroms;			// number allocated
	tsLociChrom *m_pLociChroms;		// loci chromosomes

	int m_MaxLociLen;				// maximum length of any single loci
	int LocateChrom(char *pszChrom);
	int AddChrom(char *pszChrom);

	UINT16 GenNameHash(char *pszName);
	static int SortLocii( const void *arg1, const void *arg2);
public:
	CFilterLoci(void);
	~CFilterLoci(void);
	void Reset(void);
	int Load(char *pszFile);
	bool Locate(char *pszChrom, int StartLoci, int EndLoci);

};
