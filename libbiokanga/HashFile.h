#pragma once
#include "./commdefs.h"

const unsigned int cDfltHFBuffSize = 0x03ffff;  // default buffer allocation size
const unsigned int cMaxHFBuffSize = 0x07ffff;	// max buffer allocation size

const unsigned int cMaxHashWinSize = 200;			// maximum hash window size
const unsigned int cMinHashWinSize = 10;			// minimum hash window size
const unsigned int cMaxHashRange = 26;				// hashes are generated (min(Windowsize,cMaxHashRange)) over this many bits - 1

const unsigned int cHashesPerBlock = 0x0ffff;	// number of hash buckets per hash allocation block
const unsigned int cMaxExtdMatches = 0x0fffff;  // max number of extended matches per chromosome v chromosome

typedef unsigned int tHashID;					// uniquely identifies a hash during construction (1..n)

#pragma pack(1)
// hash file entries contain the following header by which file type and construction parameters can be identified
typedef struct TAG_sHashHeader {
	unsigned int EntryID;
	unsigned int WindowSize;					// window size over which hashes were calculated
	unsigned int HashSlots;						// number of hash slots
	unsigned int NumHashes;						// total number of hashes
	unsigned int SeqLen;						// sequence length
} tsHashHeader;


typedef struct TAG_sHashEntry {
	unsigned int Psn;						// locii psn in chromosome of sequence with same hash
	tHashID NextHashID;						// reference to next sequence with same hash
} tsHashEntry;

typedef struct TAG_sHashArraySlot {
	unsigned int NumEntries;				// number of entries having same hash
	tHashID FirstHashID;					// reference to first hash entry 
}tsHashArraySlot;

typedef struct TAG_sExtdMatch {
	unsigned int ProbePsn;
	unsigned int TargPsn;
	unsigned int MatchLen;
	int Mode;
} tsExtdMatch;

#pragma pack()

class CHashFile : public CBioSeqFile
{
	tsHashHeader m_CurHashHeader;				// file header

	unsigned int m_MaxSeqLen;			// max number of bases that could be held in m_pSeq
	etSeqBase *m_pSeq;					// allocd memory to hold genomic sequence
	unsigned int m_MaxNumSlots;			// max number of slots in m_pSlots
	tsHashArraySlot *m_pSlots;			// allocd memory to hold hash array slots
	unsigned int m_MaxNumHashes;		// max number of hashes that could be held in pHashes
	tsHashEntry *m_pHashes;				// allocd memory to hold sequence psns for each hash
	unsigned int m_HashMask;

	unsigned int m_SeqEntryID;			// current sequence
	unsigned int m_HashesEntryID;
	unsigned int m_SlotsEntryID;

	tsExtdMatch *m_pExtdMatches;
	unsigned int m_NumExtdMatches;

	int m_hStatsXML;

	int m_Mode0Total;
	int m_Mode1Total;
	int m_Mode2Total;
	int m_Mode3Total;

	static etSeqBase m_Sense2AntiMap[7];

	tsHashHeader *LoadHashHeader(unsigned int EntryID);
	tsHashArraySlot *LoadHashSlots(unsigned int EntryID);
	tsHashEntry *LoadHashes(unsigned int EntryID);
	etSeqBase *LoadSeq(unsigned int EntryID);

public:
	CHashFile(void);
	~CHashFile(void);
	int Reset(bool bFlush = true);			// reset state back to that immediately following instantiation
	int Open(char *pszSeqFile,				// specifies file to open or create
			   bool bCreate = false);		// create file if it doesn't already exist, truncate if it does
	int Close(bool bFlush = true);			// closes opened file


	void GenSenseHash(unsigned int *pHash,etSeqBase *pSeqBase,unsigned int Len);
	void GenAntisenseHash(unsigned int *pHash,etSeqBase *pSeqBase,unsigned int Len);

	bool SeqSenseMatches(int ProbeLen,etSeqBase *pProbe,etSeqBase *pTarg);
	bool SeqSenseRevMatches(int ProbeLen,etSeqBase *pProbe,etSeqBase *pTarg);

	bool SeqAntisenseMatches(int ProbeLen,etSeqBase *pProbe,etSeqBase *pTarg);
	bool SeqAntisenseRevMatches(int ProbeLen,etSeqBase *pProbe,etSeqBase *pTarg);

	static unsigned int SeqSenseMatchesLeft(unsigned int ProbeLen,etSeqBase *pProbe,unsigned int ProbePsn,unsigned int TargLen,
							   etSeqBase *pTarg,unsigned int TargPsn);
	static unsigned int SeqSenseMatchesRight(unsigned int ProbeLen,etSeqBase *pProbe,unsigned int ProbePsn,unsigned int TargLen,
							   etSeqBase *pTarg,unsigned int TargPsn);
	static unsigned int SeqAntisenseMatchesLeft(unsigned int ProbeLen,etSeqBase *pProbe,unsigned int ProbePsn,unsigned int TargLen,
							   etSeqBase *pTarg,unsigned int TargPsn);
	static unsigned int SeqAntisenseMatchesRight(unsigned int ProbeLen,etSeqBase *pProbe,unsigned int ProbePsn,unsigned int TargLen,
							   etSeqBase *pTarg,unsigned int TargPsn);

	static unsigned int SeqSenseRevMatchesLeft(unsigned int ProbeLen,etSeqBase *pProbe,unsigned int ProbePsn,unsigned int TargLen,
							   etSeqBase *pTarg,unsigned int TargPsn);
	static unsigned int SeqSenseRevMatchesRight(unsigned int ProbeLen,etSeqBase *pProbe,unsigned int ProbePsn,unsigned int TargLen,
							   etSeqBase *pTarg,unsigned int TargPsn);
	static unsigned int SeqAntisenseRevMatchesLeft(unsigned int ProbeLen,etSeqBase *pProbe,unsigned int ProbePsn,unsigned int TargLen,
							   etSeqBase *pTarg,unsigned int TargPsn);
	static unsigned int SeqAntisenseRevMatchesRight(unsigned int ProbeLen,etSeqBase *pProbe,unsigned int ProbePsn,unsigned int TargLen,
							   etSeqBase *pTarg,unsigned int TargPsn);



	unsigned int								// (EntryID) uniquely identifies this hashed sequence
		AddEntry(bool bChkComplexity,		 // true if low complexity sequences are to be filtered out
				char *pszSource,		// from where was this sequence originated
				char *pszDescription,		// descriptive text about this sequence
				etSeqBase   *pSeq,			// sequence to generate hashes over
				unsigned int SeqLen,		// sequence length
				unsigned int WindowSize);	// hashes are generated over this size non-overlapping window 

	bool Save(void);						// save to disk
	bool LocateMatches(bool bChkComplexity,		 // true if low complexity sequences are to be filtered out
						unsigned int EntryID,
					   char *pszTargDescr,
					   unsigned int ProbeID,  // used to uniquely identify this probe when ProcessThisMatch is called
					   char *pszProbeDescr,
					   unsigned int MinMatchLen, // matches located must of at least this length
					   etSeqBase *pProbe,
					   unsigned int ProbeLen,
					   int hXMLdumpFile = -1);			// where to write matches as XML

	bool ExtendMatch(unsigned int EntryID,
					   unsigned int ProbeID,
   					   unsigned int ThresLen,		// threshold extended length
					   etSeqBase *pProbe,
							 unsigned int ProbeLen,
							 unsigned int ProbeMatchPsn,
							 etSeqBase *pTarget,
							 unsigned int TargLen,
							 unsigned int TargMatchPsn,
							 unsigned int MatchLen,
							 int MatchMode);

	bool IsLowComplexity(etSeqBase *pSeq,unsigned int Len);
					   
	static int CompareExtdMatches( const void *arg1, const void *arg2 );
};
