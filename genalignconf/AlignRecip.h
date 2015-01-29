#pragma once
#include "../conservlib/commhdrs.h"

const int cARMaxNameLen = 25;				// max length of any species or chromosome name
const int cARAllocAlignSeqIncr = 40000000;	// alloc memory for m_pAlignSeqs in this sized (tsARalignSeq) block
const int cARAllocChromNamesIncr = 10000;		// alloc memory for m_pChroms in this sized (tsARchromname) increments
											// NOTE: species may not be fully asembled and contain scaffolds instead of chromosomes - there can be 100K+ scaffolds!
const int cARAllocSpeciesNamesIncr = 50;		// alloc memory for m_pSpecies in this sized (tsARspeciesname) increments

#pragma pack(1)

typedef struct TAG_sARspeciesname {
	char SpeciesID;							// uniquely identifies this species
	__int64	AssembLen;						// total length of all species specific chromosomes 
	unsigned short int Hash;				// hash on szName
	char szName[cARMaxNameLen+1];			// place holder for species name + '\0'
	} tsARspeciesname;

typedef struct TAG_sARchromname {
	int ChromID;							// uniquely identifies this chromosome name
	char SpeciesID;							// chromosome in this species
	int ChromLen;							// chromosome length
	unsigned short int Hash;				// hash on szName
	char szName[cARMaxNameLen+1];			// place holder for chromosome name + '\0'
	} tsARchromname;

typedef struct TAG_sARalignSeq {
	short int SrcFileID;	// identifies original alignment file (1..2)
	int BlockID;			// identifies original alignment block (1..n)
	int SubSeqID;			// identifies subsequence alignment within the original alignment block (1..n)
	int Len;				// length of subsequence
	int Matches;			// number of exact matches in this subsequence
	int RefChromID;			// reference species chromosome identifier
	int RefOfs;				// offset in reference chromosome at which subsequence starts
	int RelChromID;			// relative species chromosome identifier
	int RelOfs;				// offset in relative chromosome at which subsequence starts adjusted for UCSC strand offset
	char Strand;			// on which strand relative to the reference is the relative alignment
} tsARalignSeq;

const int cARMaxSeqLenStats = 1000000;	// maximum aligned sequence length for distribution stats processing

typedef struct TAG_sARSeqStat {
	int	Cnt;						// number of sequences in this length bin
	int Ident25;					// number of sequences with >= 25% identity
	int Ident50;					// number of sequences with >= 50% identity
	int Ident75;					// number of sequences with >= 75% identity
	int Idents;						// total number of bases claimed identical
	int CntReciprocal;				// number of subsequences in this length bin with at least 1 base reciprocal
	int NumPutative;				// number of bases which map to an alignment in BA - may be reciprocal
	int NumPutativeMiss[9];			// number of putative bases which missed by less than 2,5,10,25,50,100,1000,10000,100000
	int Reciprocal;					// number of bases which are reciprocal 
	int ErrChrom;					// number of sequences which mapped back to a different chromosome
	int ErrStrand;					// number of sequences which mapped back to a different strand
	int ErrChromStrand;				// number of sequences which mapped back to either a different strand or chromosome
} tsARSeqStat;

typedef struct TAG_sARReciprocalStats {
	int NumClaimedAligned;	// total claimed as aligned from species A onto species B
	int NumClaimedMatched;  // total claimed as exact matches
	int NumReciprocal;		// total which are reciprocally aligned
	tsARSeqStat SeqsPlusStrand[cARMaxSeqLenStats];
	tsARSeqStat SeqsMinusStrand[cARMaxSeqLenStats];

} tsARReciprocalStats;

#pragma pack()

class CAlignRecip
{
	int m_NumAligns;			// number of alignment instances in m_pAlignSeqs
	int m_AllocdAlignSeqsInsts;	// how many instances have been allocd to m_pAlignSeqs
	tsARalignSeq *m_pAlignSeqs;	// pts to array holding instances of tsAlignSeq's
	tsARalignSeq **m_ppAlignSeqsIdx; // pts to index on m_pAlignSeqs sorted by relative chromosomes.ofs
	
	int m_NumSpecies;				// number of species
	int m_AllocdSpeciesInsts;		// how many instances have been allocd to m_pSpecies
	tsARspeciesname *m_pSpecies;

	int m_NumChroms;				// number of chromosomes
	int m_AllocdChromsInsts;		// how many instances have been allocd to m_pChroms
	tsARchromname *m_pChroms;

	tsARReciprocalStats *m_pStats;	// generated stats

	unsigned short GenNameHash(char *pszName);
	
	static int SortAlignRef( const void *arg1, const void *arg2);
	static int SortAlignRel( const void *arg1, const void *arg2);
	static int SortSpeciesChroms( const void *arg1, const void *arg2);

public:
	CAlignRecip(void);
	~CAlignRecip(void);
	int Reset(void);
	int ProcessChroms(char *pszSpeciesChromFile,bool bChkDups=false);
	int Open(char *pszAlignFile);
	int OutputStats(char *pszRsltsFile);

	int Close(void);
	int AddSpeciesChrom(char *pszSpecies,char *pszChrom,int ChromLen,bool bChkDups=true);
	int GetChromID(char *pszSpecies,char *pszChrom);
	int GetChromID(int SpeciesID,char *pszChrom);
	int GetSpeciesID(char *pszSpecies);
	int GetChromLen(int ChromID);
	__int64 GetSpeciesLen(int SpeciesID);	// returns sum of all chromosome lengths for species
	int AXTStrand2Other(char *pszSpecies,char *pszChrom,int Psn);
	int AXTStrand2Other(int ChromID,int Psn);
	int MAFStrand2Other(char *pszSpecies,char *pszChrom,int Psn);
	int MAFStrand2Other(int ChromID,int Psn);

	int	// returned number of bases which are reciprocal
		CalcSeqStats(FILE *pMissesFile,
					 int SrcFileID,				// original source file from which alignment was processed
					 int BlockID,				// alignment block identifier in source file
					 int SubSeqID,				// which subsequence in alignment block
						int Len,				// alignment length
						int Matches,			// claimed number of exact identities in this aligned sequence
						char *pszRefSpecies,	// A species
						char *pszRefChrom,      // source chromosome
						int RefOfs,				// source offset
						char *pszRelSpecies,	// B species
						char *pszRelChrom,		// B chromosome
						int RelOfs,				// B offset
						char Strand);			// B strand

	void OutputMiss(FILE *pMissesFile,
						   char Strand,
						   int Distance,
						   int SrcFileID,
						   int BlockID,
						   int SubSeqID,
						   int InitLen,
						   int RefOfs,
						   int RelOfs,
						   int Ofs,
						   tsARalignSeq *pProbe);

	int LocateReciprocals(char *pszAlignFile,char *pszErrAlignFile=NULL,char *pszMissesFile=NULL);
};


