#pragma once
#include "./commdefs.h"

const int cDSFVersion = 4;				// file header + structure version
const unsigned int cDSFVersionBack = 4;	// minimum version that is this version is still compatible with

const int cDSFAvgSegs = 512000;			// can be this many datasegs per dataset (an average over all datasets)
const int cDSFMaxDatasets = 50;			// can be this many datasets (species)
const int cDSFMaxChromsDataset = 512000;	// at most this many chromosomes in any single dataset
const int cDSFMaxTotalChroms   = 2048000;	// and this many chromosomes over all datasets
const int cDSFMaxDirEls  = cDSFAvgSegs * cDSFMaxDatasets;
const int cDSFMaxDataPts = 512000;	// max number of data points in any dataset segment


#pragma pack(1)

typedef enum TAG_etDataPtType {
	eDPTPackedBases = 0,					// packed bases - 2 per byte
	eDPTbyte,								// 8 bit
	eDPTword,								// 16bit
	eDPTdword,								// 32 bit
	eDPTdouble,								// 8 byte
	eDPTundef								// undefined
} etDataPtType;

typedef struct TAG_sFDSChromDataset {
	INT32 NumDatasegs;				// number of datasegs in this dataset on this chromosome
	INT32 ChromOfs;					// lowest offset of any dataseg in this dataset on this chromosome
	INT32 NumPts;					// total number of datapoints covered by all datasegs including between datasegs
} sFDSChromDataset;

typedef struct TAG_sFDSRelChrom {
		INT32 ChromLen;				// chromosome length (1..n) or 0 if length unknown
		INT32 ChromID;				// globally unique (1..n) chromosome identifier
		INT16 DatasetID;			// to which dataset this chromosome belongs
		UINT16 Hash;				// hash used to quickly eliminate chromosome names that can't match in searches
	    INT32 NxtChromID;			// identifier of next chromosome having same hashed name or 0 if none
		char szName[cMaxDatasetSpeciesChrom]; // chromosome name
} tsFDSRelChrom;

typedef struct TAG_sFDSRefChrom {
		INT32 NumDatasets;			// number of datasets (ref + rel) on this chromosome	
		INT32 NumDatasegs;			// number of datasegs in all datasets on this chromosome
		INT32 ChromOfs;				// lowest offset of any dataseg in any dataset on this chromosome
		INT32 NumPts;				// total number of datapoints covered by all datasegs including between datasegs
		INT32 ChromID;				// matches chromosome identifier in tsFDSRelChrom
		sFDSChromDataset Datasets[cDSFMaxDatasets]; // all dataset start/ends on this chromosome
} tsFDSRefChrom;

typedef struct TAG_sFDSDataset {
		INT32 NumSegs;						// total number of segments in this dataset
		INT32 NumChroms;					// total number of chromosomes in this dataset
		INT16 DatasetID;					// globally unique (1..n) dataset identifier
		UINT16 Hash;						// hash used to quickly eliminate names that can't match in searches
	    char szName[cMaxDatasetSpeciesChrom];	// dataset name
} tsFDSDataset;

typedef struct TAG_sDSFileHdr {
	unsigned char Magic[4];				// magic chars to identify this file as a biosequence file
	unsigned int Type;					// biosequence file type 
	unsigned int Version;				// header version, incremented if structure changes with later releases
	UINT64 FileLen;					    // current file length
	unsigned int SizeOfHdr;				// total size of this header - alignment blocks (sAlignBlocks) immediately follow
	int MaxChroms;						// max number of chromosomes supported
	int MaxDatasets;					// max number of datasets supported
	int MaxDirEls;						// max number of dataseg directory elements supported
	int MaxDataPts;						// max number of data points supported per dataseg
	int NumChroms;						// actual number of chromosomes
	int NumDatasets;					// actual number of datasets
	int NumDirEls;						// actual number of dataseg directory elements
	int DataSegSize;					// actual maximal sized dataseg
	int RefDatsetID;					// identifer for reference dataset
	int NumRefChroms;					// number of reference chromosomes in RefChroms
	etDataPtType DataPtType;			// data point type for all datapoints in this file
	int DataPtSize;						// number of bytes per DataPtType - note that eDPTPackedBases == 1 
	INT64 DirElOfs;					// file offset to dataseg directory
	tsFDSDataset Datasets[cDSFMaxDatasets];	 // directory of all datasets - reference + relative
	tsFDSRelChrom Chroms[cDSFMaxTotalChroms];// directory of all chromosomes
	tsFDSRefChrom RefChroms[cDSFMaxChromsDataset];	 // directory of all reference chromosome datasets
// added in version 2
	char szDescription[cMBSFFileDescrLen];// describes contents of file
	char szTitle[cMBSFShortFileDescrLen];	// short title by which this file can be distingished from other files in dropdown lists etc

}tsDSFileHdr;

typedef struct TAG_sDSFDirEl {
	UINT32	LoFileOfs;				// low 32bits of where on disk the associated data points are
	UINT8	HiFileOfs;				// high 8 bits of where on disk the associated data points are
	INT32 RefChromOfs;				// (sort order 3) starting offset on reference dataset chromosome
	INT32 NumPts;					// (sort order 4) number of data point values
	INT8  RelDatasetID;				// (sort order 1) relative dataset identifer
	INT8  RefChromID;				// (sort order 2) maps to this reference dataset chromosome
	INT32 RelChromOfs;				// starting offset on relative dataset chromosome
	UINT32 RelChromID:30;			// is from this relative dataset chromosome
	UINT32 RelStrand:1;				// relative strand '+' (0) or '-' (1)
} sDSFDirEl;

typedef struct TAG_sDSFCtx {
	void *m_pCachedDataset;			// used when returning structure conformation datapoints
	int m_LoadedDataChromID;
	int m_LoadedDataDatasetID;
	int m_LoadedDataChromOfs; 
	int m_LoadedDataNumPts;
	union {							// union to contain current default loaded data marker value
		  etSeqBase Base;
		  unsigned char Byte;
		  unsigned short Word;
		  unsigned int DWord;
		  double Double;
		} m_LoadedDataMarker;
} sDSFCtx;


#pragma pack()

class CDataPoints : public CConformation
{
	char m_szFile[_MAX_PATH];		// file containing this instance
	int m_hFile;					// file handle for opened alignment file
	tsDSFileHdr m_FileHdr;			// file header
	sDSFDirEl *m_pDirEls;			// array of data segment directory elements

	INT32 m_ChromHashes[0x0ffff];	// hashed indexes into chromosomes

	int m_AllocdCacheSize;
	void *m_pCachedDataset;			// used when returning structure conformation datapoints
	int m_LoadedDataChromID;
	int m_LoadedDataDatasetID;
	int m_LoadedDataChromOfs; 
	int m_LoadedDataNumPts;
	union {							// union to contain current default loaded data marker value
		  etSeqBase Base;
		  unsigned char Byte;
		  unsigned short Word;
		  unsigned int DWord;
		  double Double;
		} m_LoadedDataMarker;

	bool m_bOpenRead;				// opened in read mode
	bool m_bOpenWrite;				// opened in write mode


	void InitHdr(void);

protected:
	unsigned short GenNameHash(char *pszName);

	int UnpackBases(unsigned int SeqLen,	// unpacked sequence length 
		  etSeqBase *pUnpacked,				// pts to sequence to unpack into
		  unsigned int NibbleOfs,			// nibble to start unpack from (0..1)
		  unsigned char *pPacked);			// where to unpack from

	int PackBases(unsigned int SeqLen,	// unpacked sequence length 
		  etSeqBase *pUnpacked,				// pts to sequence to pack from
		  unsigned int NibbleOfs,			// nibble to start pack into (0..SeqLen-1)
		  unsigned char *pPacked);			// where to pack

	tsFDSRelChrom *GetRelChrom(int RelChromID);
	tsFDSRefChrom *GetRefChrom(int RefChromID);
	sFDSChromDataset *GetRefChromDataset(int RelDatasetID,int RefChromID);

	sDSFDirEl *LocateRelDirEl(int RelDatasetID,int RefChromID,int RefChromOfs,bool bClosest = false);
	int FillDataPts(void *pDataPts,void *pFiller,int DataPtSize,int NumDataPts);
	static int CompareDirEls( const void *arg1, const void *arg2 );
	int PreloadDataset(int RelDatasetID,		// relative dataset identifier
						int RefChromID,			// reference chromosome identifier
						int RefChromOfs,		// start on reference chromosome 
						int NumPts,				// number of datapoints to get
						void *pMarker = NULL);			// (optional) pts to value to use as missing data marker

public:
	CDataPoints(void);
	~CDataPoints(void);

	int Open(char *pszFile,					// specifies file to open or create
			   bool bCreate = false);		// create file or truncate if it already exists

	teBSFrsltCodes SetDescription(char *pszDescription);
	teBSFrsltCodes GetDescription(int MaxLen,char *pszDescription);
	teBSFrsltCodes SetTitle(char *pszTitle);
	teBSFrsltCodes GetTitle(int MaxLen,char *pszTitle);

	etDataPtType GetDataType(void);

	int Close(bool bWrtDirHdr=false);

	int GetNumDatasets(void);				// returns number of datasets
	int GetNumChroms(void);					// returns total number of chromosomes over all datasets
	int GetNumRefChroms(void);				// returns number of chromosomes in reference dataset
	char *GetChromName(int ChromID);		// returns chromosome name
	char *GetDatasetName(int DatasetID);	// returns dataset name
	int GetDatasetID(char *pszDataset);
	int GetChromID(int DatasetID,char *pszChromName);
	int GetNxtRefChromID(int CurID);
	
	int GetRefNumDatasets(int RefChromID);
	int GetRefNumSegs(int RelDataset,int RefChromID);
	int GetRefChromOfs(int RelDataset,int RefChromID);
	int GetRefNumPts(int RelDataset,int RefChromID);
	int GetFirstOfs(int RelDatasetID,int RefChromID);
	int GetLastOfs(int RelDatasetID,int RefChromID);
	char GetStrand(int RelDatasetID,int RefChromID,int RefChromOfs);

	int SetRefDataset(char *pszRefDataset,etDataPtType DataPtType);
	int AddRelDataset(char *pszRelDataset);
    
	int AddRefChrom(char *pszRefChrom,int ChromLen);				  // returns reference dataset chromosome identifier, adds chromosome if it does not already exist
	int AddRelChrom(char *pszRelChrom,int RelDatasetID,int ChromLen); // returns relative dataset chromosome identifier, adds chromosome if it does not already exist

	int AddDataseg(char *pszRefChrom,			// on which reference dataset chromosome
					 int RefChromOfs,			// start on reference chromosome (0..RefChromLen-1)
 					 char *pszRelDataset,		// which relative dataset
					 char *pszRelChrom,			// from which relative dataset chromosome
					 char RelStrand,			// on which relative dataset strand '+' or '-'
					 int RelChromOfs,			// start on relative chromosome ('+' == 0..RelChromLen-1,'-' == RelChromLen-1..0)
					 int NumPts,				// number of datapoints
					 void *pDataPts);			// array of data points

	int	AddDataseg(int RefChromID,		// which reference dataset chromosome
					 int RefChromOfs,		// start on reference chromosome (0..n)
					 int RelDatasetID,		// which relative dataset
					 int RelChromID,		// which relative dataset chromosome
					 char RelStrand,		// on which relative strand '+' or '-'
					 int RelChromOfs,		// start on relative chromosome 
					 int NumPts,			// number of datapoints
					 void *pDataPts);		// array of data points

	
	int GetDataPoints(int RelDatasetID,			// relative dataset identifier
						int RefChromID,			// reference chromosome identifier
						int RefChromOfs,		// start on reference chromosome (0..n) 
						int NumPts,				// number of datapoints in dataset
						void *pRetDataPts,		// where to return datapoint values
						void *pMissingMarker);	// value to use if there are missing datapoints


	int GetValuesBase(int RelDatasetID,		// relative dataset identifier
					  int RefChromID,		// reference chromosome identifier
					  int RefChromOfs,		// start on reference chromosome (0..n) 
					  int NumPts,			// number of points to return
					  unsigned char *pToRet,	// where to return etSeqBases
					  etSeqBase MissingMarker = eBaseUndef);	// value to use if missing datapoints
	
	int GetShuffledValuesBase(int RelDatasetID,// relative dataset identifier
					  int RefChromID,		// reference chromosome identifier
					  int RefChromOfs,		// start on reference chromosome (0..n) 
					  int NumPts,			// number of points to return
					  unsigned char *pToRet,	// where to return etSeqBases
					  etSeqBase MissingMarker = eBaseUndef);	// value to use if missing datapoints

	void ShuffleDinucs(etSeqBase *pSeq,int SeqLen);

	int GetValuesChar(int RelDatasetID,		// relative dataset identifier
					  int RefChromID,		// reference chromosome identifier
					  int RefChromOfs,		// start on reference chromosome (0..n) 
					  int NumPts,			// number of points to return
					  char *pToRet,			// where to return chars
					  char MissingMarker = SCHAR_MIN);	// value to use if missing datapoints

	int GetValuesShort(int RelDatasetID,	// relative dataset identifier
					  int RefChromID,		// reference chromosome identifier
					  int RefChromOfs,		// start on reference chromosome (0..n) 
					  int NumPts,			// number of points to return
					  short int *pToRet,	// where to return shorts
					  short int MissingMarker = SHRT_MIN);	// value to use if missing datapoints

	int GetValuesInt(int RelDatasetID,		// relative dataset identifier
					  int RefChromID,		// reference chromosome identifier
					  int RefChromOfs,		// start on reference chromosome (0..n) 
					  int NumPts,			// number of points to return
					  int *pToRet,	// where to return ints
					  int MissingMarker = INT_MIN);	// value to use if missing datapoints

	int GetValuesDouble(int RelDatasetID,		// relative dataset identifier
						int RefChromID,			// reference chromosome identifier
						int RefChromOfs,		// start on reference chromosome (0..n) 
						int NumPts,				// number of points to return
 					  double *pToRet,	// where to return doubles
						double MissingMarker = DBL_MIN); // value to use if missing datapoints

};
