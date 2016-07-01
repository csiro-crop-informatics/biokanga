#pragma once
#include "./commdefs.h"

const UINT32 cBSFeatVersion = 12;		// increment each time the file structure is changed
const UINT32 cBSFFeatVersionBack = 11;	// can handle all versions starting from this minimum supported version

const int cMaxNumChroms   = 50000000;		 // maximum number of chromosomes or contigs supported
const int cMaxFeatNameLen = cMaxGeneNameLen;	// max feature name length
const int cMaxNumFeats    = 100000000;	 // max number of features supported
const int cMaxNumExons	  = 8000;		 // max number of exons supported for gene features
const int cLineBuffLen    = 128 + (cMaxNumExons * 8); // max length BED line expected
const int cAllocFeatIncr = 100000000;	 // alloc/realloc memory for features in this increment (bytes) 
const int cAllocBEDChromNamesIncr = 250000000; // alloc/realloc memory for chromosome names in this increment (bytes)

const int cMinRegLen  = 0;			// min regulatory length
const int cDfltRegLen = 2000;		// default regulatory region length
const int cMaxRegLen  = 1000000;	// max regulatory region length
const int cMinSpliceOverlap = 2;	 // overlaps of features between introns and exons must be at least this to count as a splice site overlap

const int cBFSessionSyncTimeout = 10000L;// 10 second timeout on synchronising access to session instance data
const int cMinNumBFSessions = 50;		 // minimum number of sessions supported
const int cMaxNumBFSessions = 5000;		 // maximum number of sessions supported

// feature bits - note that features can concurrently overlap multiple gene components
const int cFeatBitCDS = 0x01;			// part of feature overlaps CDS
const int cFeatBit5UTR  = 0x02;			// part of feature overlaps 5'UTR
const int cFeatBit3UTR  = 0x04;			// part of feature overlaps 3'UTR
const int cFeatBitIntrons = 0x08;		// part of feature overlaps Intron
const int cFeatBitUpstream = 0x10;		// part of feature overlaps 5'upstream regulatory 	
const int cFeatBitDnstream = 0x020;		// part of feature overlaps 3'downstream regulatory
const int cIntronExonSpliceSite = 0x040; // part of feature overlaps intron3'/5'exon splice site
const int cExonIntronSpliceSite = 0x080; // part of feature overlaps exon3'/5'intron splice site
// usually none of the above feature bits set if intergenic so region == 0, but some applications which use region filtering may use the following
const int cFeatBitIG = 0x100;			// feature is intergenic


// Number of characterised genomic regions: IG,5'UP,5'UTR,CDS,INT,3'UTR,3'DN,3'SPLICE,5'SPLICE
const int cNumCharRegs = 9; 

// following are not priority ordered and are handled separately
const int cFeatBitFeature = 0x0100;		// note: feature could be a gene
const int cFeatBitExons   = 0x0200;		// 5'UTR or CDS or 3'UTR
									
const int cAnyFeatBits = cFeatBitUpstream | cFeatBit5UTR | cFeatBitCDS | cFeatBitIntrons | cFeatBit3UTR | cFeatBitDnstream | cFeatBitFeature;
const int cRegionFeatBits = cFeatBitUpstream | cFeatBit5UTR | cFeatBitCDS | cFeatBitIntrons | cFeatBit3UTR | cFeatBitDnstream;
const int cOverlaysSpliceSites = cIntronExonSpliceSite | cExonIntronSpliceSite;

const UINT32 cFeatFiltIn = 0x01;		// this feature is returned when locating features from chromosome offsets
const UINT32 cFeatFiltOut = 0x02;		// this feature is excluded when locating feature bits

typedef enum eBEDFeatureType {
	eBTGeneExons = 1,					// BED file contains features which are genes + exons. e.g knowngenes or refseq
	eBTGeneral,							// BED file contains features only with any supplementary information unparsed. e.g repeats
	eBTAnyBed,							// accept either eBTGeneExons or eBTGeneral when opening for filtering
	eBTLast								// marker setting the range of teBEDFeatureTypes
} teBEDFeatureType;

#pragma pack(1)

typedef struct TAG_sUnsortToSortID {    // used whilst intialising feature to chrom identifier mapping
	UINT32 OldID;		// previous identifier
	UINT32 NewID;		// which is this new identifier
} tsUnsortToSortID;

typedef struct TAG_sBEDchromname {
	INT32 ChromID;							// uniquely identifies this chromosome name
	INT32 NumFeatures;						// number of features on this chromosome
	INT32 FirstStartID;						// identifier of of first feature (by start) on this chromosome (also used as a temp ref to next chrom same hash)
	INT32 LastStartID;						// identifier of of last feature (by start) on this chromosome
	INT32 MaxFeatLen;						// maximum size of any feature on this chromosome
	UINT32 Hash;							// 24bit hash on szName
	char szName[cMaxDatasetSpeciesChrom];	// place holder for chromosome name
	} tsBEDchromname;

// supplemental information known to be genes with exon info (eBTGeneExons) is parsed into the following structure
// Note: if gene is on '-' strand then starts will be > ends as offsets are normalised to the '+' strand
typedef struct TAG_sGeneStruct {
	INT32 Size;							// size of this instance
	INT32	NumExons;						// number of exons in this gene
	INT32 thickStart;						// where coding starts 
	INT32 thickEnd;						// where coding ends
	INT32 ExonStartEnds[cMaxNumExons * 2]; // array of exon starts+exon ends
	} tsGeneStructure;

typedef struct TAG_sBEDfeature {
	INT32 FeatureID;					// uniquely identifies this feature
	INT32 Size;							// size in bytes of this instance when concatenated
	INT32 UserClass;					// user application feature classification (0 if unclassified)
	UINT32 FiltFlags;					// used when filtering features in/out from feature bit processing
	INT32 NameInst;						// name instance (1..n)
	tChromID ChromID;					// feature is on this chromosome
	INT32 Start;						// feature starts at this chromosome offset (0..n)				
	INT32 End;							// feature ends at this chromosome offset (Start + Len -1)
	INT32 Score;						// feature score
	UINT32 Hash;						// hash on szName
	UINT8 FeatNameLen;					// strlen(szName)
	char Strand;						// which strand feature is located on
	char szName[1];						// place holder for feature name + '\0' which is followed by any supplemental information e.g tsGeneStructure  
} tsBEDfeature;
#pragma pack()


#pragma pack(8)
typedef struct TAG_sBEDFileHdr {
	UINT8 Magic[4];						// magic chars to identify this file as a biosequence file
	UINT64 FileLen;						// current file length
	INT64 ChromNamesOfs;				// file offset to chromosome names
	INT64 FeaturesOfs;					// file offset to features
	INT64 FeatureNamesOfs;				// file offset to sorted (by name->chrom->start->end)
	INT64 FeatureChromStartsOfs;		// file offset to sorted (by chrom->start->end)
	UINT32 Type;						// biosequence file type 
	UINT32 Version;						// header version, incremented if structure changes with later releases
	UINT32 SizeOfHdr;					// total size of this header
	INT32 MaxChroms;					// maximum number of chromosomes supported
	INT32 NumChroms;					// actual number of chromosomes
	INT32 MaxFeatures;					// maximum number of features supported
	INT32 NumFeatures;					// actual number of features
	INT32 FeatType;						// what type of features are in this file
	UINT32 FeaturesSize;					// disk/memory space required to hold concatenated Features
	INT32 ChromNamesSize;				// disk/memory space required to hold concatenated chromosome names
	char szDescription[cMBSFFileDescrLen];// describes contents of file
	char szTitle[cMBSFShortFileDescrLen];	// short title by which this file can be distingished from other files in dropdown lists etc
}tsBEDFileHdr;

#pragma pack()


class CBEDfile  : protected CEndian, public CErrorCodes
{	
	CMTqsort m_mtqsort;					// muti-threaded qsort
	int m_hFile;						// opened/created file handle
	char m_szFile[_MAX_PATH+1];			// file name as opened/created
	bool m_bCreate;					    // TRUE if file opened for create 
	bool m_bHdrDirty;					// TRUE if header needs to be written to file
	bool m_bFeaturesAvail;				// TRUE if features are loaded and can be accessed through indexes

	tsBEDFileHdr m_FileHdr;

	int m_MaxParseFeats;				// if > 0 then sets the limit on number of features parsed

	char m_szLineBuff[cLineBuffLen];		// used to buffer lines of text whilst parsing input BED or GFF file text
	char m_szAttributes[cLineBuffLen];		// used to buffer BED SuppInfo or GFF3 attributes

	size_t  m_AllocChromNamesSize;		// disk/memory space required to hold concatenated chromosome names
	
	tsBEDchromname *m_pChromNames;		// pts to array of tsBEDchromname's
	UINT32 *m_pChromHashes;				// used to hold chrom name hash indexes into m_pChromNames
	
	size_t m_AllocFeaturesSize;			// disk/memory space required to hold concatenated Features
	tsBEDfeature *m_pFeatures;			// pts to  array of tsBEDfeature's - note that features vary in size so can't use m_pFeatures[x]!!!

	tsBEDfeature **m_ppFeatureNames;        // sorted (by name->chrom->start->end) array of ptrs into m_pFeatures
	tsBEDfeature **m_ppFeatureChromStarts;  // sorted (by chrom->start->end) array of ptrs into m_pFeatures

	 int m_MinScore;					// current score cutoffs
	 int m_MaxScore;					// scores on any feature must be between these scores
	 tsBEDchromname *m_pCacheChromName;  // caches last chromosome name located through LocateChromName()
	 char m_OnStrand;					// which strand features must be on '+'/'-' or '*' for either

	void InitHdr(void);
	void Reset(bool Flush = false);

	teBSFrsltCodes Disk2Hdr(char *pszBioBed, teBEDFeatureType FeatType);
	teBSFrsltCodes Hdr2Disk(void);			// write header to disk

	tsUnsortToSortID *					// matching or NULL if unable to locate OldID
		LocateU2S(UINT32 OldID,			// original feature chromid to locate
				  int NumOldIDs,		// number of old chromids to binary search over
					tsUnsortToSortID *pSortedU2s); // sorted old to new chrom id mappings

	tsBEDchromname *LocateChromName(char *pszChromName);
	tsBEDchromname *LocateChromName(int ChromID);


	tsBEDfeature *LocateFeature(char *pszFeatName);

	teBSFrsltCodes LoadFeatures(bool bCloseFile = true);
	teBSFrsltCodes ReadDisk(INT64 DiskOfs,int Len,void *pTo);
	teBSFrsltCodes Flush2Disk(void);
	teBSFrsltCodes SortFeatures(void);
	int	LocateStart(int ChromID, // feature is on which chromosome
					int OverLapsOfs, // a point on the chromosome which returned feature starts above chrom.ofs
					 int FiltInFlags, // filter out any features which do not have at least one of the specified filter flags set
					 int FiltOutFlags); // filter out any features which have at least one of the specified filter flags set

	

	 bool InInternFeat(int FeatBits,int ChromID,int StartOfs,int EndOfs);

	 char *TrimWhitespace(char *pTxt);		// trim leading and trailing whitespace, returns ptr to 1st non-whitspace char

    static int SortFeatureNames(const void *arg1, const void *arg2); // used to sort by feature name->chromid->start->end
	static int SortChromStarts(const void *arg1, const void *arg2);  // used to sort by feature chromid->start->end
	static int SortChromNames( const void *arg1, const void *arg2);  // used to sort chromosome names
	static int SortU2S( const void *arg1, const void *arg2);	// used to sort chrom ids when mapping features on to chrom initialisation

public:
	CBEDfile(void);
	~CBEDfile(void);
	
	teBSFrsltCodes Open(char *pszFileName,teBEDFeatureType FeatType=eBTAnyBed,bool bCreate=false,int MaxParseFeats=0);

	teBSFrsltCodes SetMinScore(int Score);	// sets minimum score threshold
	teBSFrsltCodes SetMaxScore(int Score);	// sets maximum score threshold
	int GetMinScore(void);					// gets minimum score threshold
	int GetMaxScore(void);					// gets maximum score threshold
	
	teBSFrsltCodes SetDescription(char *pszDescription);
	teBSFrsltCodes GetDescription(int MaxLen,char *pszDescription);
	teBSFrsltCodes SetTitle(char *pszTitle);
	teBSFrsltCodes GetTitle(int MaxLen,char *pszTitle);

	teBSFrsltCodes Close(bool bFlush2Disk = false);

	teBSFrsltCodes ProcessBedFile(char *pszFileName);		// process/parse BED format file
	teBSFrsltCodes ProcessGFF3File(char *pszFileName);		// process/parse as a GFF3 format file

	teBSFrsltCodes ProcessUltraCoreFile(char *pszFileName); // process/parse UltraCore format file
	teBSFrsltCodes ProcessGroupCoreFile(char *pszFileName,int MinLen = INT_MIN, int MaxLen = INT_MAX); // process/parse group - multispecies - UltraCores  format file

    
	teBSFrsltCodes AddFeature(char *pszFeatName,	// feature name
					 char *pszChromName,		// chromosome name
					 int Start,					// start offset (0..n) on chromosome
					 int End,					// end offset on chromosome
					 int Score,					// feature score
					 char Strand,				// on which strand feature is located
					 char *pszSuppInfo);		// supplementary information


	teBSFrsltCodes SetFilter(char *pszFeatName,UINT32 FiltFlags = cFeatFiltIn);	// sets specified feature to have FiltFlags
	teBSFrsltCodes SetFilter(int FeatID,UINT32 FiltFlags = cFeatFiltIn);	// sets specified feature to have FiltFlags
	teBSFrsltCodes SetFilters(UINT32 FiltFlags = cFeatFiltIn);			// sets all features to have specified FiltFlags
	bool Filter(int FeatID,UINT32 FiltFlags);								// returns true if feature has any of specified FiltFlags set
	bool SetStrand(char Strand);				// sets which strand features must be on in subsequent processing '+'/'-' or '*' for either

	bool ContainsGeneDetail(void);			// returns true if file contains gene detail (utr/cds/intron etc)
	int GetNextFeatureID(int CurFeatureID,int MinScore = INT_MIN,int MaxScore = INT_MAX);
	int GetPrevFeatureID(int CurFeatureID,int MinScore = INT_MIN,int MaxScore = INT_MAX);

	int GetNameInst(int FeatureID);		// returns name instance for specified feature (1..n) will be 1 if feature name is first instance

	int LocateFeatureAfter(int ChromID,				// feature is on this chromosome
					 int ChromOfs,					// feature starts on or immediately after this offset
 					 int FiltInFlags=cFeatFiltIn,		// filter out any features which do not have at least one of the specified filter flags set
					 int FiltOutFlags=cFeatFiltOut);	// filter out any features which have at least one of the specified filter flags set

	int LocateFeatureBefore(int ChromID,	// feature is on this chromosome
					 int ChromOfs,			// feature ends on or immediately before this offset
 					 int FiltInFlags=cFeatFiltIn,		// filter out any features which do not have at least one of the specified filter flags set
					 int FiltOutFlags=cFeatFiltOut); 	// filter out any features which have at least one of the specified filter flags set

	int											// returned feature identifer
		LocateFeatureIDbyName(char *pszFeatName,// feature to locate
							int Ith = 1);		// Ith instance to locate (1..n)


	teBSFrsltCodes GetFeature(int FeatureID,// feature instance identifier
					 char *pszName=NULL,	// where to return feature name
					 char *pszChrom=NULL,	// where to return chromosome name
					 int *pStart=NULL,		// where to return feature start on chromosome (0..n) 
					 int *pEnd=NULL,		// where to return feature end on chromosome
 					 int *pScore=NULL,		// where to return score
 					 char *pStrand=NULL,	// where to return strand
					 int MaxSuppInfoLen = 0,	// how much memory has been allocated to pSuppInfo
					 void *pSuppInfo = NULL);	// where to return any supplementary information


	teBSFrsltCodes InitUserClass(int DefltClass);   // initialise all features to have the specified user class
	teBSFrsltCodes SetUserClass(char *pszFeatName, int UserClass);	// set user classification for feature instance identifier
	int GetUserClass(char *pszFeatName);		// get user classification for feature instance identifier

	teBSFrsltCodes SetUserClass(int FeatureID, int UserClass);	// set user classification for feature instance identifier
	int GetUserClass(int FeatureID);		// get user classification for feature instance identifier

	int GetBEDFormatedFeature(int FeatureID,// feature instance identifier
					int BuffLen,			// allocated buffer length (expected to be at least 128 chars) 
					char *pszBuff);			// where to return the feature formated as TAB delimited BED element

	int					// error or strlen of pszBuff returned if no errors
		GetRemappedBEDFormatedFeature(int FeatureID,// feature instance identifier
						int BuffLen,			  // allocated buffer length (expected to be at least 128 chars) 
						char *pszBuff,			  // where to return the feature formated as TAB delimited BED element
						char *pszRemappedChrom,	  // remapping onto this chrom (NULL if retaining existing chrom)
						char *pszRemappedFeat,	  // remapping as this feature (NULL if retaining existing feature name)
						int  RemappedRelOfs);	  // relative loci


	teBSFrsltCodes GetChromosome(int ChromID,			// chromosome identifier
					 char *pszChrom = NULL,		// where to return chromosome name
					 int *pNumFeatures = NULL,		// where to return number of features on this chromosome
					 int *pFirstStartID = NULL,		// where to return identifier of first feature on this chromosome
					 int *pLastStartID = NULL);		// where to return identifier of last feature on this chromosome

	int											    // returned chromosome identifer
		LocateChromIDbyName(char *pszChromName);	// chromosome name to locate

	int GetFeatureChromID(int FeatureID);			// returns chromosome identifer on which feature lies

	int										 // returned feature identifier
		LocateFeatureIDonChrom(int ChromID,	 // feature is on which chromosome
							 int OverLapsOfs, // a point on the chromosome which returned feature is to overlap by at least one base
							 int Ith,		 // Ith instance to overlap (1..n)
 							 int FiltInFlags=cFeatFiltIn, // filter out any features which do not have at least one of the specified filter flags set
  							 int FiltOutFlags=cFeatFiltOut); // filter out any features which have at least one of the specified filter flags set

	int										  // returned feature identifier
		LocateFeatureIDinRangeOnChrom(int ChromID,	  // feature is on which chromosome
							 int Start,       // feature must end on or after Start
							 int End,		  // and start on or before End 
							 int Ith,		  // Ith instance to return (1..n)
 							 int FiltInFlags=cFeatFiltIn, // filter out any features which do not have at least one of the specified filter flags set
  							 int FiltOutFlags=cFeatFiltOut); // filter out any features which have at least one of the specified filter flags set

	int										  // returned number of features
		GetNumFeatures(int ChromID,			  // features are on which chromosome
					   int Start,			  // features must end on or after Start
					   int End, 		      // and start on or before End
    				 int FiltInFlags=cFeatFiltIn,    // filter out any features which do not have at least one of the specified filter flags set
    				 int FiltOutFlags=cFeatFiltOut); // filter out any features which have at least one of the specified filter flags set

	int GetFirstFeatureID(int ChromID,				// first feature is on this chromosome
    				 int FiltInFlags=cFeatFiltIn, // filter out any features which do not have at least one of the specified filter flags set
			 		int FiltOutFlags=cFeatFiltOut); // filter out any features which have at least one of the specified filter flags set

	int GetLastFeatureID(int ChromID,				// last feature is on this chromosome
    				 int FiltInFlags=cFeatFiltIn, // filter out any features which do not have at least one of the specified filter flags set
			 		int FiltOutFlags=cFeatFiltOut); // filter out any features which have at least one of the specified filter flags set

	
	int GetNumChromosomes(void);			  // returns the number of chromosomes
	int GetNumFeatures(int ChromID = 0);	  // return number of features on specified chromosome, if ChromID == 0 then returns total over all chromosomes

	int GetProvNumChromosomes(void);			  // returns provisional number of chromosomes
	int GetProvNumFeatures(void);	  // return provisional number of features over all chromosomes

	int GetFeatStart(int ChromID,int Ith);	  // returns start of the Ith feature on chromosome Ith (1..n)
	int GetFeatEnd(int ChromID,int Ith);	  // returns end of the Ith feature on chromosome Ith (1..n)

	int										  // Maps feature or transcript relative ofs to chrom loci, assumes always that the RelOfs is on the '+' strand and returns '+' strand loci 
		MapTransOfs2Loci(int FeatureID,		  // identifies which feature/transcript
		            int RelOfs,				  // relative offset from start of feature
					char *pStrand,			  // feature is on this strand
					tChromID *pChromID,			  // feature is on this chrom
					int *pLoci);			  // feature relative offset is at this chrom loci
					
	int GetFeatureOverlaps(int ChkFeatBits,	// which internal feature bits to check for
				   int FeatID,				// feature identifier
				   int StartOfs,			// start (0..n)
				   int EndOfs,				// end along chromosome start + len - 1
				   int Distance = -1);		// used if cFeatBitUpstream or cFeatBitDnstream to specify max distance away from gene start/end

	// following functions assume co-ordinates are specified along the '+' strand and are not feature specific
	int GetFeatureBits(int ChromID,int StartOfs,int EndOfs,int FeatBits = cAnyFeatBits,int Updnstream = cDfltRegLen); // returns any feature bits between start/end
	bool InAnyFeature(int ChromID, int StartOfs,int EndOfs);		  // returns true if any point between StartOfs and EndOfs overlaps any gene
	bool InAnyExon(int ChromID,int StartOfs,int EndOfs);		  // returns true if any point between StartOfs and EndOfs overlaps any exon
	bool InAnyCDS(int ChromID,int StartOfs,int EndOfs);			  // returns true if any point between StartOfs and EndOfs overlaps any CDS
	bool InAnyIntron(int ChromID,int StartOfs,int EndOfs);		  // returns true if any point between StartOfs and EndOfs overlaps any intron	
	bool InAny5UTR(int ChromID,int StartOfs,int EndOfs);		  // returns true if any point between StartOfs and EndOfs overlaps 5' UTR of any gene	
	bool InAny3UTR(int ChromID,int StartOfs,int EndOfs);		  // returns true if any point between StartOfs and EndOfs 3' UTR of any gene
	bool InAny5Upstream(int ChromID,int StartOfs,int EndOfs,int Distance); // returns true if any point between StartOfs and EndOfs overlaps 5' upstream of gene within Distance of gene start
	bool InAny3Dnstream(int ChromID,int StartOfs,int EndOfs,int Distance); // returns true if any point between StartOfs and EndOfs overlaps 5' upstream of gene within Distance of gene end

	// following functions are feature specific
	bool InFeature(int FeatureID, int StartOfs,int EndOfs);		  // returns true if any point between StartOfs and EndOfs overlaps specified feature
	bool InExon(int FeatureID, int StartOfs,int EndOfs);		  // returns true if any point between StartOfs and EndOfs overlaps any exon
	bool InCDS(int FeatureID, int StartOfs,int EndOfs);			  // returns true if any point between StartOfs and EndOfs overlaps any CDS
	bool InIntron(int FeatureID, int StartOfs,int EndOfs);		  // returns true if any point between StartOfs and EndOfs overlaps any intron	
	bool In5UTR(int FeatureID, int StartOfs,int EndOfs);		  // returns true if any point between StartOfs and EndOfs overlaps 5' UTR of gene	
	bool In3UTR(int FeatureID,int StartOfs,int EndOfs);		  // returns true if any point between StartOfs and EndOfs 3' UTR of gene
	bool In5Upstream(int FeatureID, int StartOfs,int EndOfs,int Distance); // returns true if any point between StartOfs and EndOfs overlaps 5' upstream of gene within Distance of gene start 
	bool In3Dnstream(int FeatureID, int StartOfs,int EndOfs,int Distance); // returns true if any point between StartOfs and EndOfs overlaps 5' upstream of gene within Distance of gene end 

	int GetSpliceSiteBits(int ChromID,int StartOfs,int EndOfs,int OverlapDistance);
	int GetFeatureBitsSpliceOverlaps(int FeatID, // feature identifier
				   int StartOfs,			// start along cromosome 0..n
				   int EndOfs,				// end along chromosome start + len - 1
				   int OverlapDistance = cMinSpliceOverlap);	// must overlap by at least this many bases (1..n)

	int GetNumExons(int FeatureID);	// returns number of exons - includes UTRs + CDS
	int GetNumIntrons(int FeatureID);  // returns number of introns
	int GetFeatScore(int FeatureID);		// returns score for this feature

	int GetTranscribedLen(int FeatureID);		// returns total transcribed length for specified feature
	int GetFeatLen(int FeatureID);				// returns total length for specified feature
	int GetCDSLen(int FeatureID);				// returns CDS transcribed length (0 if no CDS) for specified feature
	int Get3UTRLen(int FeatureID);				// returns 3' UTR transcribed length (0 if no 3' UTR) for specified feature
	int Get5UTRLen(int FeatureID);				// returns 5' UTR transcribed length (0 if no 5' UTR) for specified feature

	int GetCDSStart(int FeatureID);				// returns relative start offset of CDS - NOTE add to '+' strand gene start, subtract on '-' strand gene start
	int GetCDSEnd(int FeatureID);				// returns relative end offset of CDS - NOTE add to '+' strand gene start, subtract on '-' strand gene start
	int GetExonStart(int FeatureID,int ExonID); // returns relative start offset of the ExonID'th - NOTE add to '+' strand gene start, subtract on '-' strand gene start
	int GetExonEnd(int FeatureID,int ExonID);   // returns relative end offset of the ExonID'th - NOTE add to '+' strand gene start, subtract on '-' strand gene start
	int GetIntronStart(int FeatureID,int IntronID); // returns start offset of the IntronID'th - NOTE add to '+' strand gene start, subtract on '-' strand gene start
	int GetIntronEnd(int FeatureID,int IntronID); // returns end offset of the IntronID'th - NOTE add to '+' strand gene start, subtract on '-' strand gene start

	teFuncRegion MapFeatureBits2Idx(int FeatureBits);	// returns highest priority functional region from feature bits
};
