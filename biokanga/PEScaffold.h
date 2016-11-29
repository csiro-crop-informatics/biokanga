#pragma once

const int cMaxContigNameLen = 50;		// handle Contigs with names up to this length
const int cMaxPENameLen = 50;			// handle paired end names up to this length

const int cMinSeqLength = 10;			// minimal accepted sequence length
const int cMaxSeqLength = 10000;		// maximal accepted sequence length

const int cAllocContigNames	= 250000;	// alloc/realloc Contigosome names in increments of this many
const int cAllocPENames		= 25000000;	// alloc/realloc paired end identifier names in increments of this many
const int cAllocScafolds	= 25000000;	// alloc/realloc scaffolds in increments of this many

const int cHashSize = 0x0ffffff;			// use this sized hash (24bits) 

#pragma pack(1)
typedef struct TAG_sPEScaffoldContig {
	INT32 ContigID;						// uniquely identifies this contig
	char szContig[cMaxContigNameLen+1];	// which has this name
	INT32 NumClustered;					// Total number of contigs in this cluster 
	INT32 ClusterID;					// identifies the cluster containing this contig
	INT32 HashNext;						// if non-zero then identifier of next contig with same hash
	} tsPEScaffoldContig;

typedef struct TAG_sPEIdent {
	INT32 IdentID;						// uniquely identifies this PE identifier
	char szIdent[cMaxPENameLen+1];		// which has this name
	INT32 PEScafoldID;					// associated scaffold (invalid after scaffolds sorted)
	INT32 HashNext;						// if non-zero then identifier of next PEIdent with same hash
	} tsPEIdent;

typedef struct TAG_sPEScaffold {
	INT32 PEScafoldID;					// uniquely identifies this scaffold
	INT32 PE12SeqID;					// PE1/PE2 sequence identifier
	INT32 PE1ContigID;					// PE1 is aligned onto this contig (0 if alignment unknown)
	INT32 PE2ContigID;					// PE2 is aligned onto this contig (0 if alignment unknown)
	UINT8 PE1Sense:1;					// 1 if PE1 aligned sense onto PE1ContigID; 0 if aligned antisense
	UINT8 PE2Sense:1;					// 1 if PE2 aligned sense onto PE2ContigID; 0 if aligned antisense
} tsPEScaffold; 
#pragma pack()

class CPEScaffold
{
	CMTqsort m_qsort;						// muti-threaded qsort

	char m_szSeqIDTermChrs[50];				// to hold set of chars used if identifying sequence identifiers to be right trimmed

	int m_AllocdNumScaffoldContigs;			// current allocation can hold at most this many Contigs
	int m_NumScaffoldContigs;				// this many scaffold Contigs are currently in m_pScaffoldContigs
	size_t m_AllocdScaffoldContigsMem;		// m_pScaffoldContigs current memory allocation size
	tsPEScaffoldContig *m_pScaffoldContigs;	// allocated to hold scaffold Contigosome names
	int *m_pHashContigs;					// holds hashes mapping to Contig identifers 

	int m_AllocdNumPEIdents;				// current allocation can hold at most this many PE identifiers
	int m_NumPEIdents;						// this many PE identifiers are currently in m_pNumPEIdents
	size_t m_AllocdPEIdentsMem;				// m_pNumPEIdents current memory allocation size
	tsPEIdent *m_pPEIdents;					// allocated to hold PE identifier names
	int *m_pHashPEIdents;					// holds hashes mapping to PEIdent identifers 

	int m_AllocdNumScaffolds;				// current allocation can hold at most this many scaffolds
	int m_NumScaffolds;						// this many scaffolds are currently in m_pScaffoldContigs
	size_t m_AllocdScaffoldsMem;			// m_pScaffolds current memory allocation size
	tsPEScaffold *m_pScaffolds;				// allocated to hold scaffolds

	size_t m_AllocdPE2ScaffoldsMem;			// m_ppPE2Scaffolds memory allocation size - allocated to hold m_NumScaffolds ptrs
	tsPEScaffold **m_ppPE2Scaffolds;		// scaffolds sorted by PE2Contig.PE1ContigID ascending

	int m_NumClusters;						// number of scaffolded clusters
	int m_MaxNumClustered;					// largest cluster contains this many contigs

	int m_hOutFile;							// corelations to this file

	void Init(void);						// initialise state during class instantiation
	void Reset(void);						// initialise state to that immediately following class instantiation

	char *TrimWhitespace(char *pTxt);		// inplace whitespace trimming

	int AddContigName(char *pszContigName);	// register this contig name

	char *GetContigName(int ContigID);		// returns ptr to contig name

	int AddPEIdent(char *pszIdentName);		// register this sequence name

	char *GetSeqName(int SeqID);			// returns ptr to sequence name

	int
		AddScaffold(bool bPE2,				// if false then PE1, if true then PE2
			char *pszPEIdent,				// paired end indentifier used to corelate paired ends
			char *pszContig,					// PE aligns onto this Contig
			char Strand);					// '+' or '-'

	int LoadSAM(bool bPE2,					// false if loading PE1, true if loading PE2
		char *pszSAMFile);					// load alignments from this SAM file
	tsPEScaffold * LocateMateScaffold(int PE1ContigID,int PE2ContigID); // locate first instance of this scaffold
	int LocatePE1Scaffold(int ContigID); // locate 1st instance of scaffold having this Contig identifier as PE1, returns 0 if unable to locate
	int LocatePE2Scaffold(int ContigID); // locate 1st instance of scaffold having this Contig identifier as PE2, returns 0 if unable to locate
	int										// returned next CurIdx to use on a subsequent call to IterPE1s (0 if all scaffords have been iterated)
		IterPE1s(int CurIdx,				// iterates, in PE1ContigID.PE2ContigID ascending order, all tsPEScaffolds, set to 0 for 1st scaffold
				int ContigID,				// iterate for this PE1 contig
			tsPEScaffold **ppPEScaffold);	// returned scafford or NULL if all scaffords have been iterated
	int										// returned next CurIdx to use on a subsequent call to IterPE2s (0 if all scaffords have been iterated)
		IterPE2s(int CurIdx,				// iterates, in PE2ContigID.PE1ContigID ascending order, all tsPEScaffolds
				int ContigID,				// iterate for this PE2 contig
				tsPEScaffold **ppPEScaffold); // returned scafford or NULL if all scaffords have been iterated


	int	ReportCorelationships(char *pszOutFile); // report corelationships to file

	int IdentifyClusters(void);				// indentify cluster sizes
	int RecurseCluster(int ClusterID,		// identifies the cluster to be associated with all scaffolded contigs
			   int ContigID);				// starting from this contig
	
	static int SortScaffolds(const void *arg1, const void *arg2); // Sort scaffolds by PE1ContigID.PE2ContigID ascending
	static int SortPE2Scaffolds(const void *arg1, const void *arg2); // Sort scaffolds by PE2ContigID.PE1ContigID ascending

public:
	CPEScaffold(void);
	~CPEScaffold(void);

	int Process(int PMode,					// processing mode
		char *pszSeqIDTerm,					// pair sequence identifiers until this terminating character(s) - defaults to none terminating
		char *pszInPE1File,					// input PE1 file
		char *pszInPE2File,					// input PE2 file
		char *pszOutFile);					// output corelations file
	

};

