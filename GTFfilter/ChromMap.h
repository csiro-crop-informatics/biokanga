#pragma once


const int cAllocChromMaps = 1000;	// allocate for chrom maps in this many increments

typedef struct TAG_sChromMap {
	char szChrom[cMaxGeneNameLen];
	char szContig[cMaxGeneNameLen];
	UINT32 Start;
	UINT32 End;
} tsChromMap;

class CChromMap
{
	CCSVFile *m_pCSV;		// chrom map file opened as CSV
	char m_szMapFile[_MAX_PATH+1];	// chrom map file name
	int m_CurLineNum;		// line currenly being parsed
	int m_CurLineLen;		// length of current line

	int m_AllocdChromMaps;		// m_pChromMaps can hold at most this many chrom maps
	size_t  m_AllocdChromMapsMem;	// mem allocated for m_pChromMaps
	int m_CurNumChromMaps;		// m_pChromMaps currently holds this many chrom maps
	tsChromMap *m_pChromMaps;	// will be allocated to hold chrom maps

	char *TrimWhitespace(char *pTxt); // trim leading and trailing whitespace

	static int SortMapEntries(const void *pEl1,const void *pEl2);

public:
	CChromMap(void);
	~CChromMap(void);
	void Reset(void);
	int LoadMap(char *pszMapFile);		// chrom mapping file containing contig to chrom mappings
	tsChromMap *LocateMapEntry(char *pszContig); // locates map entry for specified contig
	int Map(char *pszContig,		// input: map this contig - return: chrom contig is on
			int *pStart,			// input: starting at this contig offset (1..n) - return: starts at this chrom loci (1..n)
			int *pEnd);				// input: and ending at this contig offset (start..n) - return: ends at this chrom loci (start..n)

};
