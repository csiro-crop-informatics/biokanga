#pragma once
#include "./commdefs.h"

const unsigned int cBlockExtdMatches = 0x0ffff;				 // number of matches per sSfxHeaderBlock
const unsigned int cMaxMatchSeqSeq = 20000;					 // all sequence match lengths will be clamped to no longer than this length :-)

#pragma pack(1)
// matches are saved into the following nodes which are chained in order of addition
typedef struct TAG_sSfxMatch {
	struct TAG_sSfxMatch *pPrev;			  // previous node
	struct TAG_sSfxMatch *pNext;			  // next node
	bool	IsPartner;						  // true if identified as a partner duplicate
	unsigned int NumInCluster;				  // number of other matches which overlap this instance
	unsigned int ProbePsn;					  // probe locus position (0..ProbeLen-1) within probe entry
	unsigned int TargPsn;					  // target locus position (0..n) within target entry
	unsigned int MatchLen;					  // length of match (0..L)
	unsigned int IdentCnt;					  // count of bases which exactly match between probe and target ( <= MatchLen)
	unsigned int Mode;						  // bit 0 fwd (sense), bit 1 fwd rev, bit 2 cpl, bit 3 cpl rev (antisense)
} tsSfxMatch;

typedef struct sSfxHeaderBlock {
	sSfxHeaderBlock *pNext;					// next alloc'd block
	sSfxHeaderBlock *pPrev;					// previously alloc'd block
	tsSfxMatch Matches[cBlockExtdMatches];	// all matches in this block
} sSfxHeaderBlock;

#pragma pack()

class CRsltsFile : public CErrorCodes
{
	char m_szRsltsFile[_MAX_PATH];			// results file name
	int m_hRsltsFile;						// opened results file handle
	bool m_bRsltsXML;						// results file format - false: CSV, true: XML
	bool m_bClusterMatches;					// true if matches are to be clustered (matches overlap on either probe or target)

	unsigned int m_NumExtdMatches;			// number of matches in m_pAllocdExtdMatches
	sSfxHeaderBlock *m_pFirstAllocdBlock;	// pts to first allocated matches block 
	sSfxHeaderBlock *m_pLastAllocdBlock;	// pts to last allocated matches block 
	tsSfxMatch *m_pFreeMatches;				// pts to first free match avail for allocation

	tsSfxMatch *m_pMRAMatch;				// pts to most recently added match
		
	unsigned int m_UniqueResultID;			// used to uniquely results within all result sets

	char m_szProbeDescr[cBSFSourceSize+cBSFDescriptionSize+1];
	char m_szTargDescr[cBSFSourceSize+cBSFDescriptionSize+1];
	unsigned int m_ProbeEntryID;
	unsigned int m_TargEntryID;

public:
	CRsltsFile(void);
	~CRsltsFile(void);

	int Reset(void);			// reset state back to that immediately following instantiation
	int Open(char *pszRsltsFile, // specifies file to truncate/create
				 bool bRsltsXML = false,	// false: CSV, true: XML
				 bool bClusterMatches = false); // true if matches are to be clustered (matches overlap on either probe or target)

	int Close(void);			// closes opened file

	int StartRsltSet(unsigned int TargEntryID, char *pszTarget,	// descriptive text about the target
						 unsigned int ProbeEntryID, char *pszProbe);// descriptive text about the probe
	
	int EndRsltSet(unsigned int ProbeLen,unsigned char *pProbeSeq);

	int	AddMatch(unsigned int ProbePsn,
					 unsigned int TargPsn,
					 unsigned int MatchLen,
					 unsigned int IdentCnt,
					 unsigned int Mode);

	bool IsAlreadyMatched(unsigned int ProbePsn,
						unsigned int TargPsn,
						unsigned int MatchLen,
						unsigned int Mode);

	void MarkDuplicateMatches(void);

	bool NearlyExtend(unsigned int ProbePsn,
                        unsigned int TargPsn,
					    unsigned int MatchLen,
					    unsigned int IdentCnt,		// how many bases in match were exact matches
					    unsigned int Mode);

	void ClusterMatches(void);
	tsSfxMatch *LocateFirstPPsnMatch(unsigned int psnstart,unsigned int psnend);
	tsSfxMatch *LocateNextPPsnMatch(tsSfxMatch *pExtdMatch,unsigned int psnstart,unsigned int psnend);
	tsSfxMatch *LocateFirstTPsnMatch(unsigned int psnstart,unsigned int psnend);
	tsSfxMatch *LocateNextTPsnMatch(tsSfxMatch *pExtdMatch,unsigned int psnstart,unsigned int psnend);


};
