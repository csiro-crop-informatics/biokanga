#pragma once
#include "./commdefs.h"

const int cAllocAlignSeqIncr = 40000000;	// alloc memory for m_pAlignSeqs in this sized (tsASalignSeq) block
const int cAllocChromNamesIncr = 10000;		// alloc memory for m_pChroms in this sized (tsASchromname) increments
											// NOTE: species may not be fully asembled and contain scaffolds instead of chromosomes - there can be 100K+ scaffolds!
const int cAllocSpeciesNamesIncr = 50;		// alloc memory for m_pSpecies in this sized (tsASspeciesname) increments

#pragma pack(1)

typedef struct TAG_sASspeciesname {
	char SpeciesID;							// uniquely identifies this species
	INT64	AssembLen;						// total length of all species specific chromosomes 
	unsigned short int Hash;				// hash on szName
	char szName[cMaxDatasetSpeciesChrom];			// place holder for species name
	} tsASspeciesname;

typedef struct TAG_sASchromname {
	int ChromID;							// uniquely identifies this chromosome name
	char SpeciesID;							// chromosome in this species
	int ChromLen;							// chromosome length
	unsigned short int Hash;				// hash on szName
	char szName[cMaxDatasetSpeciesChrom];			// place holder for chromosome name
	} tsASchromname;

typedef struct TAG_sASalignSeq {
	short int SrcFileID;	// identifies original AXT or MAF file (1..n)
	int BlockID;			// identifies original alignment block (1..n)
	int SubSeqID;			// identifies subsequence alignment within the original alignment block (1..n)
	int Len;				// length of subsequence
	int Matches;			// number of exact matches in this subsequence
	int RefChromID;			// reference chromosome identifier
	int RefOfs;				// offset in reference chromosome at which subsequence starts
	int RelChromID;			// relative chromosome identifier
	int RelOfs;				// offset in relative chromosome at which subsequence starts adjusted for UCSC strand offset
	char Strand;			// on which strand relative to the reference is the alignment
} tsASalignSeq;

const int cMaxSeqLenStats = 20000;	// maximum aligned sequence length for distribution stats processing

typedef struct sSeqStat {
	int	Cnt;						// number of sequences in this length bin
	int Ident25;					// number of sequences with >= 25% identity
	int Ident50;					// number of sequences with >= 50% identity
	int Ident75;					// number of sequences with >= 75% identity
	int Idents;						// total number of bases claimed identical
	int CntReciprocal;				// number of sequences in this length bin with at least 1 base reciprocal
	int NumPutative;				// number of bases which map to an alignment in BA - may be reciprocal
	int NumPutativeMiss[9];			// number of putative bases which missed by less than 2,5,10,25,50,100,1000,10000,100000
	int Reciprocal;					// number of bases which are reciprocal
	int ErrChrom;					// number of sequences which mapped back to a different chromosome
	int ErrStrand;					// number of sequences which mapped back to a different strand
	int ErrChromStrand;				// number of sequences which mapped back to either a different strand or chromosome
} tsSeqStat;

typedef struct sReciprocalStats {
	int NumClaimedAligned;	// total claimed as aligned from species A onto species B
	int NumClaimedMatched;  // total claimed as exact matches
	int NumReciprocal;		// total which are reciprocally aligned
	tsSeqStat SeqsPlusStrand[cMaxSeqLenStats];
	tsSeqStat SeqsMinusStrand[cMaxSeqLenStats];

} tsReciprocalStats;

#pragma pack()

class CAlignValidate   : public CErrorCodes
{
	int m_NumAligns;			// number of alignment instances in m_pAlignSeqs
	int m_AllocdAlignSeqsInsts;	// how many instances have been allocd to m_pAlignSeqs
	tsASalignSeq *m_pAlignSeqs;	// pts to array holding instances of tsAlignSeq's
	tsASalignSeq **m_ppAlignSeqsIdx; // pts to index on m_pAlignSeqs sorted by relative chromosomes.ofs
	
	int m_NumSpecies;				// number of species
	int m_AllocdSpeciesInsts;		// how many instances have been allocd to m_pSpecies
	tsASspeciesname *m_pSpecies;

	int m_NumChroms;				// number of chromosomes
	int m_AllocdChromsInsts;		// how many instances have been allocd to m_pChroms
	tsASchromname *m_pChroms;

	tsReciprocalStats *m_pStats;	// generated stats

	unsigned short int GenNameHash(char *pszName);
	
	static int SortAlignRef( const void *arg1, const void *arg2);
	static int SortAlignRel( const void *arg1, const void *arg2);
	static int SortSpeciesChroms( const void *arg1, const void *arg2);

public:
	CAlignValidate(void);
	~CAlignValidate(void);
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
	INT64 GetSpeciesLen(int SpeciesID);	// returns sum of all chromosome lengths for species
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
						   tsASalignSeq *pProbe);

	int LocateReciprocals(char *pszAlignFile,char *pszErrAlignFile=NULL,char *pszMissesFile=NULL);
};


