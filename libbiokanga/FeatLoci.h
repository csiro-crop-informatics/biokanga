#pragma once

typedef struct TAG_sLoci {
	int RefID;			// reference identifier
	int ChromID;		// which chromosome
	int StartLoci;		// start
	int EndLoci;		// end
} tsLoci;

typedef struct TAG_sChrom {
	int ChromID;		 // uniquely identifies chromosome
	int LociID;			 // index + 1 of Loci on this chromosome
	unsigned short Hash; // hash over chromosome name
	char szName[80];	 // chromosome name
} tsChrom;

class CFeatLoci
{
	int m_NumLoci;
	int m_AllocLoci;
	tsLoci *m_pLocii;
	int m_NumChroms;
	int m_AllocChroms;
	tsChrom *m_pChroms;
	int AddChrom(char *pszChrom);
	unsigned short GenNameHash(char *pszName);
	static int SortLoci( const void *arg1, const void *arg2);

public:
	CFeatLoci(void);
	~CFeatLoci(void);
	int GetNumLoci(void);
	int LocateChrom(char *pszChrom);
	int AddLoci(int RefID,char *pszChrom, int StartLoci, int EndLoci);
	int ClusterLoci(int MaxDistance);
	int GetLoci(int LociID,int *pRefID,char **ppszChrom,int *pStartLoci,int *pEndLoci);
	int Locate1stLoci(int ChromID);			// returns 1st loci (1..n) on specified chromosome
	int Locate1stLoci(char *pszChrom);		// returns 1st loci (1..n) on specified chromosome
	int LocateNxtLociChrom(int LociID);		// returns next loci on same chromosome as that for LociID (1..n)

};

