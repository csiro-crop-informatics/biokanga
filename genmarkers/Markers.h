#pragma once
#include "../libbiokanga/commdefs.h"

const int cMaxLenName = 100;			// accept species or chrom/contig names of at most this length
const int cMaxMarkerSpecies = 25;		// allow at most 25 different species
const int cMaxSeqID = 100000000;		// allow at most 10^8 different target species sequences
const int cDfltThreads = 8;				// by default use 8 threads when sorting
const int cMaxWorkerThreads = 128;			// limiting max number of threads to this many

// following two constants are used when determing species specific SNPs
const int cAltMaxBases = 1;			// must be no more than this number of bases in alternative species
const double cMinBaseThres = 0.50;  // and species base must be at least this proportion of total bases

#pragma pack(1)
typedef struct TAG_sAlignLoci {
	UINT32 AlignID;		// uniquely identifies this alignment instance (1..N)
	UINT8 TargSpeciesID;// identifies aligned to cultivar or species
	UINT32 TargSeqID;	// identifies aligned to sequence - could be a chrom/contig/transcript
	UINT32 TargLoci;	// loci within SeqID at which SNPs observed
	UINT8 TargRefBase;	// loci is this reference base
	UINT8 ProbeSpeciesID;	// identifies cultivar or species with sequences aligning to TargSpecies
	UINT8 FiltLowTotBases:1;	// 1 if this alignment has fewer TotBases than reporting threshold
	UINT8 NumSpeciesWithCnts;	// number of species at this loci which have TotBases >= reporting threshold
	UINT32 TotBases;			// sum of base counts in ProbeBaseCnts
	UINT8 CultSpecBase;			// cultivar specific base allowing identification of this cultivar
	UINT8 CultSpecBaseConf;		// confidence to ascribe to CultSpecBase (0..100)
	UINT32 ProbeBaseCnts[5];	// indexed by A,C,G,T,N : number instances probe base aligned to TargRefBase 
} tsAlignLoci;

const int cAllocAlignLoci = 2000000;	// initially allocate to hold this many alignments

typedef struct TAG_sSNPSSpecies {
	UINT32 SpeciesID;	// uniquely identifies this species instance
	UINT8 szSpecies[cMaxLenName+1];	// species name (allow for trailing '\0' terminator)
	UINT8 IsRefSpecies:1;	// set if this is a reference or target species, reset if a probe species
} tsSNPSSpecies;


typedef struct TAG_sSeqName {
	UINT32 SeqID;			// uniquely identifies this sequence name (1..N)
	UINT8 Len;				// length of this tsSeqName instance 
	UINT64 NxtSeqOfs;		// offset into m_pAllocSeqNames at which next sequence with same name hash starts
	UINT8 szSeqName[1];		// sequence name including terminating '\0'
} tsSeqName;
const size_t cAllocSeqNames = 50000;	// allocate to incrementally hold this many sequence names
const size_t cAllocMemSeqNames = (sizeof(tsSeqName) + cMaxLenName) * cAllocSeqNames; // allocate in this sized increments memory (m_pAllocSeqNames) for holding sequence names
const size_t cAllocMinDiffSeqNames = (sizeof(tsSeqName) + cMaxLenName) * 5; // reallocate if less than this many bytes remaining 
#pragma pack()

class CMarkers
{

	int m_MaxQSortThreads;						// max number of threads to use when sorting
	CMTqsort m_MTqsort;							// multithreaded qsort
	bool m_bSorted;								// set true if alignments sorted

	CHyperEls *m_pHypers;					// to hold alignments from csv, bed or sam file

	tsSNPSSpecies *m_pCurSpecies;			// currently processed species
	UINT8 m_NumSpecies;						// current number of species in m_Species (also includes the reference species)
	UINT8 m_RefSpeciesID;						// identifer for species identified as being the reference species
	tsSNPSSpecies m_Species[cMaxMarkerSpecies];	// array of currently known species

	UINT32 m_NumSeqNames;			// currently there are this many sequence names
	UINT64 m_UsedMemSeqNames;		// memory currently used for sequence names
	UINT64 m_AllocMemSeqNames;		// memory allocated for sequence names
	tsSeqName *m_pAllocSeqNames;	// allocated to hold sequence names

	UINT64  m_AllocMemSeqNameIDsOfs;	// memory allocated for sequence m_pSeqNameIDsOfs
	UINT64 *m_pAllocSeqNameIDsOfs;		// allocated to hold sequence identifiers to sequence offsets in m_pAllocSeqNames

	UINT32 m_UsedNameHashArray;			// currently using this entries in the SeqNameHashArray
	UINT64 *m_pSeqNameHashArray;		// allocated to hold offsets into m_pAllocSeqNames for sequence name hashes

	UINT8 m_szCurSeqName[cMaxLenName+1];	// holds last processed sequence name
	UINT32 m_CurSeqNameID;				// and it's corresponding identifer

	UINT32 m_UsedAlignLoci;				// currently using this many alignment loci
	UINT32 m_AllocAlignLoci;			// allocated to hold this many alignment loci
	size_t m_AllocMemAlignLoci;			// allocation memory size
	tsAlignLoci *m_pAllocAlignLoci;		// allocated to hold alignment loci 

	int AddLoci(UINT8 TargSpeciesID,		// reads were aligned to this cultivar or species
				UINT32 TargSeqID,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				UINT32 TargLoci,			// loci within target sequence at which SNPs observed
				etSeqBase TargRefBase,	// loci has this reference base
				UINT8 ProbeSpeciesID,	// reads were aligned from this cultivar or species
				UINT32 ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				UINT32 ProbeCntC,		// number instances probe base C aligned to TargRefBase
				UINT32 ProbeCntG,		// number instances probe base G aligned to TargRefBase
				UINT32 ProbeCntT,		// number instances probe base T aligned to TargRefBase
				UINT32 ProbeCntN);		// number instances probe base N aligned to TargRefBase

public:
	CMarkers(void);
	~CMarkers(void);
	void Reset(void);	// clears all allocated resources

	int SetMaxThreads(int MaxThreads);		// maximum number of threads to be utilised
	int		// qsorts alignment loci by TargSeqID,TargLoci,ProbeSpeciesID ascending
		SortTargSeqLociSpecies(void);

	int 
		LoadSNPFile(int MinBases,			// accept SNPs with at least this number covering bases
					  double MaxPValue,		// accept SNPs with at most this P-value
					  char *pszRefSpecies,	// this is the reference species 
					  char *pszProbeSpecies, // this species reads were aligned to the reference species from which SNPs were called 
					  char *pszSNPFile);	// SNP file to parse and load

	int 
		AddImputedAlignments(int MinBases,			// must be at least this number of reads covering the SNP loci
					  char *pszRefSpecies,				// this is the reference species 
					char *pszProbeSpecies,				// this species reads were aligned to the reference species from which SNPs were called 
					char *pszAlignFile,					// file containing alignments
					int FType = 0);							// input alignment file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM)

	UINT8									// returned species identifer (1..cMaxSpecies)
		AddSpecies(char *pszSpecies,bool IsRefSpecies = false);	// cultivar or species

	char *								// returned species name corresponding to SpeciesID
		SpeciesIDtoName(UINT8 SpeciesID); // species identifier

	UINT8 NameToSpeciesID(char *pszSpecies); // returned species identifier for specified name, returns 0 if name not previously added with AddSpecies)

	bool								// true if a reference or target species
		IsRefSpecies(UINT8 SpeciesID);	// species to test

	UINT32								// returned sequence identifier (1..cMaxSeqID)
		AddTargSeq(char *pszSeqName);	// sequence name - could be chrom, contig, transcript name

	UINT32 NameToSeqID(char *pszSeqName); // returned sequence identifier for specified name, returns 0 if name not previously added with AddTargSeq)

	char *								// returned sequence name
		SeqIDtoName(UINT32 SeqID);		// sequence identifier for which name is to be returned

	int AddLoci(char *pszTargSpecies,	// reads were aligned to this cultivar or species
				char *pszTargSeq,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				UINT32 TargLoci,		// loci within target sequence at which SNPs observed
				etSeqBase TargRefBase,	// loci has this reference base
				char *pszProbeSpecies,	// reads were aligned from this cultivar or species
				UINT32 ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				UINT32 ProbeCntC,		// number instances probe base C aligned to TargRefBase
				UINT32 ProbeCntG,		// number instances probe base G aligned to TargRefBase
				UINT32 ProbeCntT,		// number instances probe base T aligned to TargRefBase
				UINT32 ProbeCntN);		// number instances probe base N aligned to TargRefBase

	int IdentSpeciesSpec(int AltMaxCnt,	// max count allowed for base being processed in any other species
						int MinCnt,		// min count required for base being processed in species
						double Thres,	// to be processed a base must be at least this proportion of total
						int MinSpeciesWithCnts = 0,			// must be at least this number of species with base counts more than MinSpeciesTotCntThres - 0 if no limit 
						int MinSpeciesTotCntThres = 0);		// individual species must have at least this number of total bases at SNP loci - 0 if no threshold


	int Report(char *pszRefGenome,			    // reference genome assembly against which other species were aligned
			int NumRelGenomes,					// number of relative genome names
			char *pszRelGenomes[],				// relative genome names
			char *pszReportFile,				// report to this file
			int MinSpeciesWithCnts = 0,			// must be at least this number of species with base counts more than MinSpeciesTotCntThres - 0 if no limit 
			int MinSpeciesTotCntThres = 0);		// individual species must have at least this number of total bases at SNP loci - 0 if no threshold
};


