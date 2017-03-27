/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once
#include "../libbiokanga/commdefs.h"

const int cMaxLenName = 100;			// accept species or chrom/contig names of at most this length
const int cMaxMarkerSpecies = 1000;		// allow at most 1000 different species or cultivars
const int cMaxSeqID = 100000000;		// allow at most 10^8 different target species sequences

// following two constants are used when determining species specific SNPs
const int cAltMaxBases = 1;			// must be no more than this number of bases in alternative species
const double cMinBaseThres = 0.50;  // and species base must be at least this proportion of total bases

const UINT16 cFlgSNPcnts = 0x01;		// alignment base counts parsed from SNP file
const UINT16 cFlgImputCnts = 0x02;      // alignment base counts imputed from coverage and assumed to be the target reference base
const UINT16 cFlgAlignCnts = 0x04;		// alignment base counts imputed from alignment sequences


#pragma pack(1)
typedef struct TAG_sAlignLoci {
	INT64 AlignID;				// uniquely identifies this alignment instance (1..N)
	UINT16 TargSpeciesID;		// identifies aligned to cultivar or species
	UINT32 TargSeqID;			// identifies aligned to sequence - could be a chrom/contig/transcript
	UINT32 TargLoci;			// loci within SeqID at which SNPs observed
	UINT8 TargRefBase;			// loci is this reference base
	UINT16 ProbeSpeciesID;		// identifies cultivar or species with sequences aligning to TargSpecies
	UINT8 FiltLowTotBases:1;	// 1 if this alignment has fewer TotBases than reporting threshold
	UINT16 NumSpeciesWithCnts;	// number of species at this loci which have TotBases >= reporting threshold
	UINT32 TotBases;			// sum of base counts in ProbeBaseCnts
	UINT8 CultSpecBase;			// cultivar specific base allowing identification of this cultivar
	UINT8 CultSpecBaseConf;		// confidence to ascribe to CultSpecBase (0..100)
	UINT16 Flags;				// any loci associated flags
	UINT32 ProbeBaseCnts[5];	// indexed by A,C,G,T,N : number instances probe base aligned to TargRefBase 
} tsAlignLoci;

const int cAllocAlignLoci = 40000000;	// initially allocate to hold this many alignments
const int cReAllocAlignPerc = 125;	    // if required to realloc then realloc by this percentage of the current allocation

typedef struct TAG_sSNPSSpecies {
	UINT16 SpeciesID;	// uniquely identifies this species instance
	UINT8 szSpecies[cMaxLenName+1];	// species name (allow for trailing '\0' terminator)
	UINT8 IsRefSpecies:1;	// set if this is a reference or target species, reset if a probe species
} tsSNPSSpecies;


typedef struct TAG_sSeqName {
	UINT32 SeqID;			// uniquely identifies this sequence name (1..N)
	UINT8 Len;				// length of this tsSeqName instance 
	UINT64 NxtSeqOfs;		// offset into m_pAllocSeqNames at which next sequence with same name hash starts
	UINT8 szSeqName[1];		// sequence name including terminating '\0'
} tsSeqName;
const size_t cAllocSeqNames = 5000000;	// allocate to incrementally hold this many sequence names
const size_t cAllocMemSeqNames = (sizeof(tsSeqName) + cMaxLenName) * cAllocSeqNames; // allocate in this sized increments memory (m_pAllocSeqNames) for holding sequence names
const size_t cAllocMinDiffSeqNames = (sizeof(tsSeqName) + cMaxLenName) * 100; // reallocate if less than this many bytes remaining 
#pragma pack()

class CMarkers
{

	CHyperEls *m_pHypers;					// to hold alignments from csv, bed or sam file

	tsSNPSSpecies *m_pCurSpecies;			// currently processed species
	UINT16 m_NumSpecies;						// current number of species in m_Species (also includes the reference species)
	UINT16 m_RefSpeciesID;						// identifer for species identified as being the reference species
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

	INT64 m_UsedAlignLoci;				// currently using this many alignment loci
	INT64 m_AllocAlignLoci;				// allocated to hold this many alignment loci
	size_t m_AllocMemAlignLoci;			// allocation memory size
	tsAlignLoci *m_pAllocAlignLoci;		// allocated to hold alignment loci 


	INT64 AddLoci(UINT16 TargSpeciesID,		// reads were aligned to this cultivar or species
				UINT32 TargSeqID,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				UINT32 TargLoci,			// loci within target sequence at which SNPs observed
				etSeqBase TargRefBase,	// loci has this reference base
				UINT16 ProbeSpeciesID,	// reads were aligned from this cultivar or species
				UINT32 ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				UINT32 ProbeCntC,		// number instances probe base C aligned to TargRefBase
				UINT32 ProbeCntG,		// number instances probe base G aligned to TargRefBase
				UINT32 ProbeCntT,		// number instances probe base T aligned to TargRefBase
				UINT32 ProbeCntN,		// number instances probe base N aligned to TargRefBase
				UINT16 Flags);			// any loci associated flags

	bool m_bSorted;								// set true if alignments sorted
	static int QSortAlignSeqLociSpecies(const void *arg1, const void *arg2); // qsorts alignment loci by TargSeqID,TargLoci,ProbeSpeciesID ascending

public:
	CMarkers(void);
	~CMarkers(void);
	void Reset(void);	// clears all allocated resources

	INT64		// qsorts alignment loci by TargSeqID,TargLoci,ProbeSpeciesID ascending
		SortTargSeqLociSpecies(void);

	int 
		LoadSNPFile(int MinBases,			// accept SNPs with at least this number covering bases
					  double MaxPValue,		// accept SNPs with at most this P-value
					  char *pszRefSpecies,	// this is the reference species 
					  char *pszProbeSpecies, // this species reads were aligned to the reference species from which SNPs were called 
					  char *pszSNPFile);	// SNP file to parse and load

	INT64 
		AddImputedAlignments(int MinBases,		// must be at least this number of reads covering the SNP loci
					  char *pszRefSpecies,		// this is the reference species 
					char *pszProbeSpecies,		// this species reads were aligned to the reference species from which SNPs were called 
					char *pszAlignFile,			// file containing alignments
					int FType = 0,				// input alignment file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM)
					bool bSeqs = false,			// if alignment file contains the read sequence then impute bases from the actual sequences	
					int EstNumSeqs = 0,			// estimated number of sequences (0 if no estimate)
					int EstSeqLen = 0);			// estimated mean sequence length (0 if no estimate)			

	UINT16									// returned species identifer (1..cMaxSpecies)
		AddSpecies(char *pszSpecies,bool IsRefSpecies = false);	// cultivar or species

	char *								// returned species name corresponding to SpeciesID
		SpeciesIDtoName(UINT16 SpeciesID); // species identifier

	UINT16 NameToSpeciesID(char *pszSpecies); // returned species identifier for specified name, returns 0 if name not previously added with AddSpecies)

	bool								// true if a reference or target species
		IsRefSpecies(UINT16 SpeciesID);	// species to test

	UINT32								// returned sequence identifier (1..cMaxSeqID)
		AddTargSeq(char *pszSeqName);	// sequence name - could be chrom, contig, transcript name

	UINT32 NameToSeqID(char *pszSeqName); // returned sequence identifier for specified name, returns 0 if name not previously added with AddTargSeq)

	char *								// returned sequence name
		SeqIDtoName(UINT32 SeqID);		// sequence identifier for which name is to be returned

	int PreAllocSNPs(INT64 EstNumSNPS);	// preallocate memory for this many estimated SNP loci
	int PreAllocSeqs(int EstNumSeqs, int MeanSeqLen); // preallocate memory for this estimate of number sequences having this mean sequence length


	INT64 AddLoci(char *pszTargSpecies,	// reads were aligned to this cultivar or species
				char *pszTargSeq,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				UINT32 TargLoci,		// loci within target sequence at which SNPs observed
				etSeqBase TargRefBase,	// loci has this reference base
				char *pszProbeSpecies,	// reads were aligned from this cultivar or species
				UINT32 ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				UINT32 ProbeCntC,		// number instances probe base C aligned to TargRefBase
				UINT32 ProbeCntG,		// number instances probe base G aligned to TargRefBase
				UINT32 ProbeCntT,		// number instances probe base T aligned to TargRefBase
				UINT32 ProbeCntN);		// number instances probe base N aligned to TargRefBase

	int IdentSpeciesSpec(int AltMaxCnt,	// max count allowed for base being processed in any other species, 0 if no limit
						int MinCnt,		// min count required for base being processed in species
						double SNPMmajorPC,		// to be processed major putative SNP base must be at least this percentage of total
						int MinSpeciesWithCnts = 0,			// must be at least this number of species with base counts more than MinSpeciesTotCntThres - 0 if no limit 
						int MinSpeciesTotCntThres = 0);		// individual species must have at least this number of total bases at SNP loci - 0 if no threshold

	INT64 NumAlignLoci(void);					// returns current number of alignment/SNP loci

	INT64											// number of markers reported
		Report(char *pszRefGenome,			    // reference genome assembly against which other species were aligned
			int NumRelGenomes,					// number of relative genome names
			char *pszRelGenomes[],				// relative genome names
			char *pszReportFile,				// report to this file
			int MinSpeciesWithCnts = 0,			// must be at least this number of species with base counts more than MinSpeciesTotCntThres - 0 if no limit 
			int MinSpeciesTotCntThres = 0,  	// individual species must have at least this number of total bases at SNP loci - 0 if no threshold
			bool bSloughRefOnly = false);		// do not report if no inter-cultivar SNP marker, i.e if cultivars all same with the polymorthic site relative to reference only 
};


