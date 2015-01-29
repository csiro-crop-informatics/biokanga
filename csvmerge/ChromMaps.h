#pragma once

#pragma pack(1)
typedef struct TAG_sChromMap {
	struct TAG_sChromMap *pNext;	// next allocated chrom map
	char szChrom[cMaxDatasetSpeciesChrom];
	int StartLoci;
	int EndLoci;
	int MapLen;						// number of unsigned int's allocated to pMap
	unsigned int *pMap;			
} tsChromMap;
#pragma pack()

typedef enum {
	eElIntersect = 0,				// Ref and Rel intersect (ref,rel must both be present)
	eElRefExclusive,				// Ref exclusive	(ref only)
	eElRelExclusive,				// Rel exclusive	(rel only)
	eElRefRelUnion,					// Ref and Rel union  (ref or rel)
	eElRefNotRefRel					// Ref or Rel not present
	} etElMapType;

class CChromMaps
{
	int m_NumChromMaps;				// number of chrom maps in use
	tsChromMap *m_pChromMaps;		// pts to linked list of chrom maps
	tsChromMap **m_ppSortedChromMaps; // pts to array of ptrs to chrom maps sorted by chrom name
	tsChromMap *m_pLastLocChrom;	// pts to last located chrom
	int InitMaps(void);				// initialise chrom maps
	tsChromMap *LocateChromMap(char *pszChrom);	// locates chrom map for specified chromosome
	int SetChromMap(bool bIsRel,char *pszChrom,int Start,int End);
	int AddChrom(char *pszChrom,int StartLoci,int EndLoci);
	int LoadChroms(bool bSetMaps,bool bIsRel,bool bSkipFirst,int MinLength,int MaxLength,int FlankExtend,char *pszCSVFile);
	tsChromMap *LocateChrom(char *pszChrom);
	static int SortChromNames( const void *arg1, const void *arg2);

public:
	CChromMaps(void);
	~CChromMaps(void);
	int Load(bool bSkipFirst,int MinLength,int MaxLength,int RefFlankExtend,char *pCSVRef,int RelFlankExtend=0,char *pCSVRel=NULL);
	int Process(etElMapType Mode,int (*Handler)(char *pszChrom,int StartLoci,int EndLoci,void *pParm),void *pParm);

};
