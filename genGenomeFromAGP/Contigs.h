#pragma once

const int cMaxContigDescrLen = 128;		// accept contig descriptors of this length max incl of terminating '\0'

const int cContigAllocEntries = 10000;	// alloc in increments of this many contig entries
const int cContigSeqChunk = 0x0fffff;	// alloc for holding contig sequence chunks
const int cContigSeqAlloc = cContigSeqChunk * 50;	// alloc to hold contig sequences in this byte sized increments

typedef struct TAG_sContigEntry {
	UINT32 ContigID;					// uniquely identifes this contig
	char szDescr[cMaxContigDescrLen];	// full contig descriptor line
	char szContig[cMaxDatasetSpeciesChrom];	 // parsed out contig identifier
	UINT32 ContigLen;		// length of this contig
	UINT32 ContigSeqOfs;	// offset in m_pContigSeqs at which this contig sequence starts
	} tsContigEntry;

class CContigs
{
	UINT8 *m_pContigChunk;				// used to buffer contig sequence chunks
	size_t m_ContigSeqsOfs;				// offset in m_pContigSeqs at which to next write
	size_t m_AllocdContigSeqMem;		// total memory allocated to m_pContigSeqs 
	UINT8 *m_pContigSeqs;				// allocated to hold concatenated contig sequences

	UINT32 m_NumContigEntries;			// number of contig entries
	UINT32 m_AllocContigEntries;		// this many contig entries have been allocated
	size_t m_AllocContigEntriesMem;	    // total memory allocated for contig entries
	tsContigEntry *m_pContigEntries;	// allocated to hold contig entries

	teBSFrsltCodes LoadContigFile(char *pszFile);

	static int SortContigs(const void *pEl1,const void *pEl2);

public:
	CContigs(void);
	~CContigs(void);
	void Reset();

	teBSFrsltCodes LoadContigs(char *pszFile);		// process for contigs from this (can be wildcarded) file(s)
	tsContigEntry *LocateContig(char *pszContig);   // returns specified contig 
	UINT8 *LocateSeq(char *pszCompID,UINT32 Start); // returns ptr to sequence


};

