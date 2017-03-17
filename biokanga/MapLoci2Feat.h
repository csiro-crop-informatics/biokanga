#pragma once

const int cDfltJoinOverlap = 0;			// join cores which overlap and have starts which at most differ by cMaxJoinOverlap
const int cMaxJoinOverlap = 100;		// max allowed join overlap

const int cMaxLengthRange = 1000000;	// maximal element length

										// processing modes
typedef enum TAG_eMLFPMode
{
	ePMdefault,					// default is to associate features along full length of cores
	ePMstarts,					// associate features with starts only
	ePMdyad,					// associate features with assumed dyad positioned 73nt downstream of 5' core start
	ePMplaceholder				// used to set the enumeration range
} etMLFPMode;

// strand processing modes
typedef enum TAG_eStrandProc
{
	eStrandDflt,			// default is to ignore the strand
	eStrandSense,			// process for sense
	eStrandAnti,			// process for antisense
	eStrandPlaceholder
} etStrandProc;

// feature isoform processing
typedef enum TAG_eISOFProc
{
	eISOFPRPKM,					// default is to report the feature isoform with maximal RPKM
	eISOFReads,					// report the feature isoform with maximal total reads
	eISOFPall,				    // report all feature isoforms
	eISOFPlaceholder
} etISOFProc;

#pragma pack(1)
typedef struct TAG_sFeatCntDist
{
	int FeatID;							// identifies this feature
	int GeneLen;						// gene length
	int TranscribedLen;					// transcribed length
	char szName[80];					// gene name
	int RegionCnts[8];					// region counts
	double RPKM;						// RPKM for this feature
	double RelAbundance;				// relative abundance (sum of all the reciprocals for each read's RelScale)
	double LenNormRelAbundance;			// RelAbundance normalised to transcript length
	double TransRelAbundance;			// RelAbundance as a proportion of the sum of all LenNormRelAbundance
	int UniqueReadLociHits;				// number of unique loci in this transcript to which at least one read mapped
	double UniqueHitsRelAbundance;		// RelAbundance normalised to number of UniqueReadLociHits in this transcript
	double TransUniqueHitsRelAbundance; // RelAbundance as a proportion of the sum of all UniqueHitsRelAbundance
	int NumExonReads;					// number of reads in UTR's plus CDS
	bool bIsIsoform;					// set true if feature is an isoform, assumed if name at least 9 long and suffixed by ".[0-99]"
	bool bMaxRPKM;						// set true if feature is an isoform with the maximal RPKM
	bool bMaxExonReads;				    // set true if feature is an isoform with the maximal reads in UTR's+CDS

} tsFeatCntDist;

typedef struct TAG_sChromRegionCnts
{
	int CDS;							// cnts in CDS
	int UTR5;							// cnts in 5'UTR
	int UTR3;							// cnts in 3'UTR
	int Intronic;						// cnts in intronic
	int upstream5;						// cnts in 5'upstream
	int downstream3;					// cnts in 3'downstream
	int intron3exon;					// cnts in intron3'/5'exon
	int exon3intron;					// cnts in exon3'/5'intron
	int Intergenic;						// cnts in Intergenic
} tsChromRegionCnts;

#pragma pack()

class CMapLoci2Feat
{
	int m_AllocdChromRegionCnts;				// m_pChromRegionCnts allocated for this many chroms
	tsChromRegionCnts	*m_pChromRegionCnts;	// to hold per chrom region counts

	CHyperEls *m_pHypers;
	tsFeatCntDist *m_pFeatCntDists;		// to hold feature count distributions
	CBEDfile *m_pBiobed;

	etMLFPMode m_MLFPMode;					// processing mode
	etStrandProc m_StrandProc;			// strand processing
	etISOFProc m_IsoformRprt;    		// feature isoform reporting mode
	int m_RegRegionLen;					// regulatory region length
	UINT32 m_NumEls;					// total number of elements to be processed

	int m_NumSplitEls;					// number of elements identified as being split elements
	bool m_bFeatinsts;					// true if counts are to be accumulated on a feature unique basis instead of any features covering the genomic locus
	bool m_bOneCntRead;					// true if only single counts to be accumulated per read
	int m_hRsltFile;

	int m_hFeatRsltFile;

public:
	CMapLoci2Feat();
	~CMapLoci2Feat();

	void MLFReset(void);
	int MapLoci2Features(char *pszRsltsFile);
	bool IsSameFeature(char *pszFeatA, char *pszFeatB);
	bool			// true if suffix was trimmed
		TrimNameIso(char *pszName); // Inplace remove any name isoform suffix of the form '.[0-99]'
									
	bool IsIsoform(char *pszName);	// assumes that if the feature name is at least cMinNameRootLen long and suffixed by '.[0-99]' then thats an isoform
	int MapFeatures2Loci(char *pszFeatRsltsFile);

	int MLFProcess(etMLFPMode PMode,				// processing mode
								  bool bDedupe,				// true if input elements are to be deduped
								  bool bFeatinsts,			// true if input elements are to be associated to individual features, false if to all features at that locus
								  bool bOneCntRead,			// true if one count per read rule to be applied (functional regions are prioritised with CDS as the highest) 
								  etISOFProc IsoformRprt,		// feature isoform reporting mode
								  etStrandProc StrandProc,	// how to process read + element strand
								  int FType,					// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM
								  teCSVFormat CSVFormat,		// if CSV input then expected file format
								  char *pszInLociFile,		// input CSV, BED or SAM loci file
								  char *pszInBEDFile,			// input BED file
								  char *pszRsltsFile,			// output loci to feature mapping file
								  char *pszFeatRsltsFile,		// optional feature mapping results file
								  char *pszSummRsltsFile,		// optional output chrom summary results file
								  int RegRegionLen,			// regulatory region length
								  int MinLength,				// minimum element length
								  int MaxLength,				// maximum element length
								  int JoinOverlap);			// deduping join overlap

};

