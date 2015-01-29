#pragma once

const int cMaxNumAssocFiles = 20;			// alllow for at most this many feature association files 

const int cChromSeqReAlloc = 5000000;	// realloc chrom sequence size
const UINT32 cMaxAccumCnt= 0x07ffffff;  // truncate accumulated counts to be at most this value
const UINT32 cRetainCnt  = 0x08000000;  // accumulated counts >= this value are counts marked as being to be retained

const size_t cROIAlloc = 10000;		// incrementally allocate for this many ROIs

const int cMinMedianCov  = 1;			// minimum user specified region median coverage
const int cDfltMedianCov = 2;			// default minimum region median coverage
const int cMaxMedianCov = 1000;			// maximum user specified minimum median coverage
const int cMinRegionLen  = 50;			// minimum user specified region length
const int cDfltRegionLen = 100;			// default minimum region length
const int cMaxRegionLen = 100000;		// maximum user specified minimum region length

const int cMinGapLen = 0;				// minimum user specified gap without any aligned reads in reported regions
const int cDfltGapLen = 10;				// default user specified gap without any aligned reads in reported regions
const int cMaxGapLen = 100;				// maximum user specified gap without any aligned reads in reported regions

const int cDfltMaxDist = 1000000000;	// used to flag that any feature is at least this distance away
const char *cpszNoFeat = "FxNxL";		// unlocated features are assigned this name

const int cMaxChromCov = 1000;			// can handle at most this many chromosomes
const int cAllocCovCnts = 0x07fffff;	// allocate for chrom coverage cnts in this sized increments

// processing modes
typedef enum TAG_ePROIMode {		
	ePROIMdefault,					// Auto determine
	ePROIMCsv,						// CSV loci
	ePROIMBed,						// UCSC BED format
	ePROIMSam,						// SAM format
	ePROIMplaceholder				// used to set the enumeration range
	} etPROIMode;

// strand processing modes
typedef enum TAG_eStrandProc {
		eStrandDflt,			// default is to ignore the strand
		eStrandWatson,			// process for Watson
		eStrandCrick,			// process for Crick
		eStrandPlaceholder
} etStrandProc;

// output format modes
typedef enum TAG_eFROIMode {
	eFROIsumCSV,					// default is for summary CSV in which the minimum distance from all BEDs is reported 
	eFROIallCSV,					// CSV in which the minimum distance for each BED is reported 
	eFROIbed,						// UCSC BED format
	eFROIplaceholder				// used to set the enumeration range
	} etFROIMode;

typedef enum eBEDRegion {
	eMEGRAny = 0,				// process any region
	eMEGRIntergenic,			// only process intergenic
	eMEGRExons,					// only process exons
	eMEGRIntrons,				// only process introns
	eMEGRCDS,					// only process CDSs
	eMEGUTR,					// only process UTRs
	eMEG5UTR,					// only process 5'UTRs
	eMEG3UTR					// only process 3'UTRs
} etBEDRegion;

#pragma pack(1)

typedef struct TAG_sChromCnts {
	char szChrom[cMaxDatasetSpeciesChrom+1];	// coverage is on this chromosome
	int AllocCovCnts;							// allocated to hold cnts for this sized chromosome
	int StartOfs;								// pCovCnts[offset] of first coverage cnt
	int EndOfs;									// pCovCnts[offset] of last coverage cnt
	UINT32 *pCovCnts;							// coverage counts
} tsChromCnts;


typedef struct TAG_sROI {
	int RegionID;					// unique region identifier
	int ChromIdx;					// ROI is on this chromosome (m_ChromCnts[])
	int StartOfRegion;				// region starts at this loci
	int EndOfRegion;				// region ends at this loci
	char Strand;					// and is on this strand
	double BPKM;					// bases per thousand per Million aligned bases

	char USFeatStrand;				// nearest upstream feature is on this strand
	int USFeatDist;					// nearest upstream feature distance
	char szUSFeatID[cMaxGeneNameLen];	// nearest upstream feature name
	int USFeatFileID;				// nearest upstream feature file (in m_pszFeatFiles[USFeatFileID-1])

	char DSFeatStrand;				// nearest downstream feature is on this strand
	int DSFeatDist;					// nearest downstream feature distance
	char szDSFeatID[cMaxGeneNameLen];	// nearest downstream feature name
	int DSFeatFileID;				// nearest downstream feature file (in m_pszFeatFiles[DSFeatFileID-1])
} tsROI;
#pragma pack()




class CLocateROI
{
	int m_NumOfROIs;					// current total number of ROIs in m_pROIs
	int m_AllocNumROIs;					// m_pROIs has been allocated to hold at most this number of ROIs
	size_t	m_AllocNumROIsMem;			// actual memory size allocated
	tsROI *m_pROIs;						// allocated to hold all ROIs
	UINT32 m_NumAcceptedReads;			// total number of accepted reads or element loci

	tsChromCnts m_ChromCnts[cMaxChromCov];
	int m_NumChromsCov;		

	int m_hRsltsFile;						// handle for opened results file

	char m_ReadStrand;						// accept reads on this strand
	char m_RestainStrand;					// filter ROI by filter elements on this strand 
	char m_FeatDistStrand;					// distances to features on this strand 

	int m_NumFeatFiles;						// number of feature file names in m_pszFeatFiles[] following
	char *m_pszFeatFiles[cMaxNumAssocFiles];	// feature file names

	CCSVFile *m_pCSVFile;					// used if processing input CSV file for coverage
	CBEDfile *m_pBEDFile;					// used if processing input BED files for coverage

	CBEDfile *m_pDistBEDFile;				// current distance annotation bed file

	char *TrimWhitespace(char *pTxt);
	void Init(void);
	void Reset(void);
	int
		BuildReadCoverage(char *pszChrom,		// coverage is onto this chrom
				  int StartOfs,				// coverage start at this offset 
				  int EndOfs,				// and ends at this offset inclusive
				  int Cnt);					// increment coverage by this
	int
		RetainReadCoverage(char *pszChrom,		// filtering is onto this chrom
			  int StartOfs,				// filtering is to start at this offset 
			  int EndOfs);				// and ends at this offset inclusive

	int
		IdentROI(etFROIMode FMode,					// output in this format to m_hRsltsFile
		  int MinMedianCov,				// minimum median coverage required in reported regions
			int MinRegionLen,			// report regions which are of at least this length
			int MaxGapLen,				// report regions containing no read gaps of <= this length 
			bool bNoFilter);				// true if all counts to be processed

	int
		WriteRegions(bool bFinal,					// true if last call to this function
			 char *pszAnnoFile,				// file containing annotated features for distance calculations
			 int AnnoFileID,				// annotated file identifer
			 etFROIMode FMode,					// output in this format to m_hRsltsFile
			char *pszTitle,					// CSV species or title used if BED output format
			 char FeatDistStrand);			// distance to features on this strand 

	int											// returns number of elements accepted or error result if < 0
		GenCSVWiggle(int Limit,					// limit (0 if no limit) processing to this many elements total
				 char Strand,					// ROI strand
				 char *pszInFile);				// CSV loci input file

	int 
	GenBEDWiggle(int Limit,						// limit (0 if no limit) processing to this many bases total
					 char ROIStrand,			// ROI strand
					 char *pszInFile);			// UCSC BED input file

	int
		GenSAMWiggle(int Limit,					// limit (0 if no limit) processing to this many bases total
				 char ROIStrand,				// ROI strand
				 char *pszInFile);				// SAM input file

	bool								// returns false if errors
		DistanceToFeatures(int AnnoFileID,		// annotated file identifer
				   char FiltStrand,				// only interested in distances to nearest feature which is on this strand
				   int ROIChromID,				// ROI is on this chromosome
				   char *pszROIChrom,			// ROI chrom name
				   tsROI *pROI);

	int 
		FilterRegions(etBEDRegion Region,		// which regions are of interest
				 char ROIStrand,			// region on this strand
				 char *pszInFile);			// UCSC BED containing regions

public:
	CLocateROI();
	~CLocateROI();

	int
		Process(etPROIMode PMode,					// processing mode
			etFROIMode FMode,					// output format - CSV or BED
			char *pszTitle,					// CSV species or title used if BED output format
			char ReadStrand,				// only accept reads from this strand ('*' if either)
			char RestainStrand,				// filter ROI by filter elements on this strand ('*' if either)
			char FeatDistStrand,			// distances to features on this strand ('*' if either)
			etBEDRegion Region,				// filter by retaining this functional region
			int Limit,						// limit (0 if no limit) processing to this many bases total
			int MinMedianCov,				// minimin median coverage required in reported regions
			int MinRegionLen,				// report regions which are of at least this length
			int MaxGapLen,					// report regions containing no read gaps of <= this length 
			char *pszRetainFile,			// file containing annotated regions to retain
			int NumAssocFiles,				// number of association files in pszAssocFiles[]
			char *pszAssocFiles[],			// distance association files
			char *pszInFile,				// input CSV or BED file containing read alignment loci
			char *pszRsltsFile);				// output CSV region loci file

};

