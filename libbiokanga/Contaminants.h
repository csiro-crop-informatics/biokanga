#pragma once
// Contaminants are processed as one of two contaminate classes, class of a sequence is specified with naming convention in the contaminate file
// Flank: Can form either a prefix or suffix of the query sequence - an example of these could be adaptor sequence contamination
// Vector: Can completely contain the query sequence - an example of these could be cloning vector contamination

#include "commdefs.h"

const int cMinContamQuerySeqLen = 20;				// only process query (usually a NGS read sequence) sequences which are at least this length
const int cMaxContamQuerySeqLen = 2000;			// only process query (usually a NGS read sequence) sequences which are at most this length

const int cMinContaminantLen = 4;			// minimum length adaptor contaminant accepted
const int cMaxContaminantLen = 200;			// maximum length adaptor contaminant accepted
const int cMaxNumContaminants = (200 * 8);	// allow at most this many adaptor contaminant sequences to be loaded, need to allow for the multiple orientations a user may request
const int cAllocNumContaminants = ((cMaxNumContaminants+7)/8);	// initially alloc, then realloc as may be required, for this many Contaminants at a time		
const int cAllocNumContamNodes = (cAllocNumContaminants*cMaxContaminantLen); // initially alloc, then realloc as may be required, for this many contaminant nodes	

const int cMaxNumVectors = 10;				// allow at most this many vector contaminant sequences to be loaded
const int cMinVectorSeqLen = 100;			// vector sequences must be of at least this length
const int cMaxVectorSeqLen = 0x0ffffff;		// vector sequences can be up to this length - allows for yeast and complete bacteria

// contaminants are classed as one of the following
typedef enum TAG_eContamClass {
	eCCAllContam = 0,	// covers both flank and vector contaminates
	eCCVectContam,		// vector contaminate, contaminate (ex. BAC) expected to completely contain the read
	eCCFlankContam		// flank sequence contaminate, contaminate (ex. adaptor sequence) partially overlaps onto read flank
} teContamClass;


// Contaminant sequences are used for overlay processing according to one of the following 
typedef enum TAG_eContamType {
	eAOF5PE1Targ = 0,	// match for Contaminant sequences valid for overlaying onto 5' of a SE/PE1 target sequence
	eAOF5PE2Targ,		// match for Contaminant sequences valid for overlaying onto 5' of a PE2 target sequence
	eAOF3PE1Targ,		// match for Contaminant sequences valid for overlaying onto 3' of a SE/PE1 target sequence
	eAOF3PE2Targ,		// match for Contaminant sequences valid for overlaying onto 3' of a PE2 target sequence
	eAOFVector,			// match for contaminant sequence completely containing the target sequence
	eAOFPlaceholder
} teContamType;

#pragma pack(1)

// Vector sequence characterisation
typedef struct TAG_sVectContam {
	UINT16 ContamID;					// uniquely identifies this vector sequence
	char szName[cMaxGeneNameLen];		// vector name
	UINT8 FlgPE1Sense:1;				// check for sense overlaps of PE1 reads
	UINT8 FlgPE1Antisense:1;			// check for antisense overlaps	of PE1 reads			
	UINT8 FlgPE2Sense:1;				// check for sense overlaps of PE2 reads
	UINT8 FlgPE2Antisense:1;			// check for antisense overlaps	of PE2 reads				
	UINT32 HitTot;						// number of times this vector sequence contained a query read sequence
	INT32 ContamLen;					// length of vector sequence
	etSeqBase *pBases;					// allocated to hold the vector sequence bases
	INT32 *pSfxIdx;						// allocated to hold suffix index over pBases
} tsVectContam;

// Flank contaminant characterisation
typedef struct TAG_sFlankContam {
	UINT16 ContamID;					// uniquely identifies this Contaminant sequence
	teContamType Type;					// flank type of this contaminant
	UINT8 FlgRevCpl:1;					// 1 if contaminant sequence has been ReCpl'd relative to when loaded from contaminants file
	char szName[cMaxGeneNameLen];		// Contaminant name
	UINT32 HitTot;						// number of times this Contaminant sequence was overlapping onto a target sequence
	UINT8 ContamLen;					// length of Contaminant sequence
	UINT32 HitDist[cMaxContaminantLen+1];	// overlap length hit count distribution
	etSeqBase Bases[cMaxContaminantLen+1];// holds the Contaminant sequence bases
} tsFlankContam;

typedef struct TAG_sContaminantType {
	teContamType Type;						// identifies the overlay type of contaminants
	UINT32 NumChecks;						// number of times this type was checked for an overlap onto a target sequence
	UINT32 HitTot;							// number of times this type was overlapping onto a target sequence
	UINT32 HitDist[cMaxContaminantLen+1];	// overlap length hit count distribution for all contaminants of this type
	int NumContaminants;					// number of contaminants of this type
	int MaxContamSeqLen;					// longest contaminant sequence length of this type
	int MinContamSeqLen;					// shortest Contaminant sequence length	of this type
	tsFlankContam *pFirstContam;			// pts to first contaminant of this type
	tsFlankContam *pLastContam;				// pts to last contaminant of this type
	UINT32 RootContamSeqNodeIdx;			// index + 1 (0 if none) into m_pContamSeqNodes[] of root node for contaminant sequences of this type
} tsContaminantType;

typedef struct TAG_sContamSeqNodeBase {
	UINT8 Base;								// contaminant base
	UINT8 MaxContamSfxLen;					// maximum suffix length of any contaminate prefix sequence ending with current base
	UINT16 ContamID;						// if overlay onto target prefix starts with this contaminant suffix then attribute to this contaminate
	UINT32 ChildNodeIdx;					// index of child contaminate sequence node
} tsContamSeqNodeBase;

typedef struct TAG_sContamSeqNode {
	UINT8 NumNodeBases;						// number of node bases in this contaminant sequence node
	tsContamSeqNodeBase Bases[5];			// upto 5 node bases (a,c,g,t,n)
} tsContamSeqNode;

#pragma pack()


class CContaminants
{
	char szContaminantFile[_MAX_PATH];		// file from which Contaminant sequences were loaded
	UINT8 *m_pSeqBuff;						// allocated to buffer sequences when parsing from file
	tsContaminantType m_ContaminantTypes[eAOFPlaceholder];	// indexed by teContamType
	int m_TotNumContaminants;				// total number of Contaminants loaded into either m_pContaminants or m_ContaminantVectors
	int m_NumFlankContaminates;				// number of flanking contaminants
	int m_NumVectContaminates;				// number of contamiant sequences in m_ContaminantVectors
	int m_MaxFlankContamSeqLen;				// longest flank contaminant sequence length of any type
	int m_MinFlankContamSeqLen;				// shortest flank contaminant sequence length of any type
	int m_MaxVectContamSeqLen;				// longest Vect contaminant sequence length of any type
	int m_MinVectContamSeqLen;				// shortest Vect contaminant sequence length of any type


	int m_AllocContaminants;				// m_pContaminants is allocated to hold this many Contaminants
	size_t m_AllocContaminantsMem;			// current memory allocation size for holding Contaminants
	tsFlankContam *m_pContaminants;			 // allocated to hold loaded Contaminants, note that after loading these are sorted by length decending
	CMTqsort m_mtqsort;						// muti-threaded qsort
	
	UINT32 m_NumContamSeqNodes;				// currently using this many nodes in contamiant sequences B-Tree
	size_t m_AllocdContamSeqNodes;			// allocated to hold this many contamiant nodes
	size_t m_AllocdContamSeqNodesMem;		// memory size allocated to hold contamiant nodes
	tsContamSeqNode *m_pContamSeqNodes;		// ptr to allocated nodes

	tsVectContam m_ContaminantVectors[cMaxNumVectors]; // to hold vector contamiant sequences

	int m_CacheContamID;				    // if non-zero then last accessed contaminate identifier
	teContamClass m_CacheContamIDClass;     // if m_CacheContamID non-zero then relevant contaminate class
	int m_CacheContamIDIdx;					// if m_CacheContamID non-zero then relevant index for contaminate class

	
	static int SortContamTypeLen(const void *arg1, const void *arg2); // sort Contaminants by Type ascending then ContamLen descending
	static int IndexVectorSeq(const void *arg1, const void *arg2);	// generate index over vector contaiminant sequence in m_pCurVector2Index
	static tsVectContam *m_pCurVector2Index; // current vector sequence being indexed with IndexVectorSeq()

	int										// returned contaminate identifier
		AddVectContam(bool bPE1Sense,		// check for PE1 contained sense
					bool bPE1Antisense,		// check for PE1 contained antisense
					bool bPE2Sense,			// check for PE2 contained sense
					bool bPE2Antisense,		// check for PE2 contained antisense
					char *pszName,			// vector name
					int ContamLen,			// vector sequence is of this length
					etSeqBase *pContamSeq);	// the vector sequence


	int										// returned contaminate identifier
		AddFlankContam(teContamType Type,	// contaminant is of this overlay type
					bool bRevCpl,			// contaminant sequence has been RevCpl'd
					char *pszName,			// Contaminant name
					int ContamLen,			// Contaminant sequence is of this length
					etSeqBase *pContamSeq);	// the Contaminant sequence
	
	int										// When last Contaminant has been added with AddContaminant then Finalise() must be called to
		FinaliseContaminants(void);			// sort Contaminants by length descending and to initialise of m_ContaminantTypes[]
		

	int IndexContamSeq(teContamType Type,	// contaminant sequence is of this overlay type
					 UINT16 ContamID,		// identifies contaminant sequence
					 UINT8 SeqLen,			// contaminant sequence is of this length
					 etSeqBase *pSeq);		// contaminant sequence

	int
		MatchVectContam(int AllowSubsRate,			// if non-zero then allow substitutions in the overlapping Contaminants at this rate per 25bp of overlap length if overlap >= 10bp
					int QueryLen,			// query sequence length
					etSeqBase *pQuerySeq,	// attempt to locate a maximal containment onto this query read sequence
					tsVectContam *pVectContam); // checking for containment from this vector sequence


	int			// index+1 in pSfxArray of first exactly matching probe or 0 if no match
		LocateFirstExact(etSeqBase *pProbe,  // pts to probe sequence
				  int ProbeLen,					// probe length to exactly match over
				  etSeqBase *pTarg,				// target sequence
				  int *pSfxArray,				// target sequence suffix array
				  int SfxLo,					// low index in pSfxArray
				  int SfxHi);					// high index in pSfxArray

	int			// index+1 in pSfxArray of last exactly matching probe or 0 if no match
		LocateLastExact(etSeqBase *pProbe, // pts to probe sequence
				  int ProbeLen,					// probe length to exactly match over
				  etSeqBase *pTarg,				// target sequence
				  int *pSfxArray,				// target sequence suffix array
				  int SfxLo,					// low index in pSfxArray
				  int SfxHi);					// high index in pSfxArray

	int			// index+1 in pSfxArray of last exactly matching probe or 0 if no match
		LocateNextExact(int CurSfxIdx,			// previously returned index either from LocateFirstExact() or LocateNextExact()
				  etSeqBase *pProbe,			// pts to probe sequence
				  int ProbeLen,					// probe length to exactly match over
				  etSeqBase *pTarg,				// target sequence
				  int *pSfxArray,				// target sequence suffix array
				  int SfxLo,					// low index in pSfxArray
				  int SfxHi);					// high index in pSfxArray


	int			// < 0 if no matches, otherwise the number of substitutions required for the match 
		RecursiveMatch(bool bSuffixOverlaps,  // true if processing for contaminant suffix overlaps onto target prefix, false if processing for contaminant prefix overlaps onto target suffix
				int MaxSubs,				// maximum allowed substitutions
				int TargLen,				// target sequence length 
				etSeqBase *pTargSeq,		// target sequence	
				tsContamSeqNode *pCurNode,  // current node
				int *pContamID);				// match was to this contamination sequence

	bool m_bSerialCreated;					// set true if serialisation rwlocks created/initialised
#ifdef _WIN32
	CRITICAL_SECTION m_hSCritSect;
#else
	pthread_spinlock_t m_hSpinLock;
#endif
	void AcquireSerialise(void);					// serialise access
	void ReleaseSerialise(void);

public:
	CContaminants(void);
	~CContaminants(void);

	void Reset(void);					// reset back to state immediately following class instantiation
	int Init(void);						// initialise class state

	char *									// returned text for Type
		ContaminateType2Txt(teContamType Type);	// contaminant is of this overlay type
		
	int									// < 0 if errors, 0..N the number of contaminant sequences loaded
		LoadContaminantsFile(char *pszContaminantsFile);	// load contaminant sequences from this multifasta file

	int			// 0 if no Contaminant contains this target sequence
		MatchVectContams(bool bIsPE2,		// true if target sequence is a PE2, else if false then target is either a SE or PE1
					int AllowSubsRate,		// if non-zero then allow substitutions in the overlapping Contaminants at this rate per 25bp of overlap length if overlap >= 10bp
					int QueryLen,			// query sequence length
					etSeqBase *pQuerySeq);	// attempt to locate a contaminate vector sequence containing this query read sequence

	int			// 0 if no Contaminant overlap, 1..N number of Contaminant suffix bases overlaping onto target pTargSeq
		MatchContaminants(teContamType Type,		// process for this overlay type
					int AllowSubsRate,		// if non-zero then allow substitutions in the overlapping Contaminants at this rate per 25bp of overlap length with a minimum of 1 subs allowed
					int MinOverlap,			// minimum required overlap
					int QueryLen,			// query sequence length
					etSeqBase *pQuerySeq);	// attempt to locate a contaminate flanking sequence which overlays onto this query read sequence



	int NumOfContaminants(teContamClass ComtamClass= eCCAllContam);				// returns total number of both flank and vector contaminants
	int	MaxContaminantLen(teContamClass ComtamClass= eCCAllContam);				// returns longest length of any contaminant sequence
	int	MinContaminantLen(teContamClass ComtamClass= eCCAllContam);				// returns shortest length of any contaminant sequence

	UINT32 NumChecks(teContamType Type);	// returns number of times this contaminant type was checked for an overlap onto a target sequence

	teContamType							// returned contaminant type
		ContaminantType(int ContamID);		// contaminant identifier
	char *									// returned contaminant name
		ContaminantName(int ContamID);		// contaminant identifier

	int										// returned contaminant identifier
		ContaminantID(char *pszName);		// returns contaminant identifer for name

	teContamClass							// returned contaminant class
		ContaminantClass(int ContamID);		// contaminant identifier

	int										// returned contaminant length
		ContaminantLen(int ContamID);		// contaminant identifier
		
	int										// total number of overlaps 	 					
		ContaminantDist(int ContamID,		// contaminant distributions for this contaminant
						int LenCnts = 0,	// return for up to this length overlaps, ignored for vector contaminates
						int *pCnts = NULL);	// returned overlap length counts, ignored for vector contaminates

};

