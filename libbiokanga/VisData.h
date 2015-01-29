#pragma once


const int cVisDataPtVersion = 1;			// file header + structure version
const int cMaxVisDatasegs = 1000000;		// can be this many datasegs per dataset per chromosome
const int cMaxVisDatasets = 1000;			// can be this many datasets per chromosome
const int cMaxVisChroms   = 100;			// and will handle this many chromosomes

const int cMaxVisDataPts = 512000;		// max number of data points in any dataset segment
const int cMaxVizDatasetNameLen = 25;		// max length of any dataset name
const int cMaxVizChromNameLen = 25;			// max length of any chromosome name

#pragma pack(1)

// Datasets may be sparse and hence are represented as segments - a dataset contains at least one segment
typedef struct TAG_sDataseg {
	INT64 Nxt;							// file offset of next related tsDataseg 
	INT64 Prv;							// file offset of previous related tsDataseg
	unsigned int Size;						// total size of this instance
	unsigned __int16 DatasetID;				// which dataset
	unsigned __int16 ChromID;				// chromosome
	unsigned int ChromOfs;					// starting offset on chromosome
	unsigned int NumPts;					// number of data point values
	double DataPts[1];						// data point value
} tsDataseg;


typedef struct TAG_sDataset {
	unsigned int DatasetID;					// identifies the dataset
	unsigned int NumSegs;					// total number of segments
	INT64	 FirstSeg;						// file offset of first dataseg in this dataset
	INT64  LastSeg;						// file offset of last dataseg in this dataset
	unsigned int FirstChromOfs;				// first chromosome offset at which any segment starts
	unsigned int LastChromOfs;				// last chromosome offset at which any segment starts
} tsDataset;

typedef struct TAG_sVisChrom {
		unsigned __int16 ChromID;				// globally unique (1..n) chromosome identifier
		unsigned __int16 Hash;					// hash used to quickly eliminate chromosome names that can't match in searches
		unsigned __int16 NameLen;				// strlen(szName)
	    char szName[cMaxVizChromNameLen+1];		// chromosome name
		unsigned __int16 NumDatasets;			// number of datasets represented in this chromosome
		tsDataset Datasets[cMaxVisDatasets];	
} tsVisChrom;

typedef struct TAG_sVisDataset {
		unsigned __int16 DatasetID;				// globally unique (1..n) dataset identifier
		unsigned __int16 Hash;					// hash used to quickly eliminate names that can't match in searches
		unsigned __int16 NameLen;				// strlen(szName)
	    char szName[cMaxVizDatasetNameLen+1];	// dataset name
} tsVisDatasetName;

typedef struct TAG_sVisHdr {
	unsigned char Magic[4];				// magic chars to identify this file as a biosequence file
	unsigned __int32 Type;				// biosequence file type 
	unsigned __int32 Version;			// header version, incremented if structure changes with later releases
	UINT64 FileLen;						// current file length
	unsigned __int32 SizeOfHdr;			// total size of this header - alignment blocks (sAlignBlocks) immediately follow
	teDataType DataType;				// datatype of values
	unsigned __int16 MaxChroms;			// max number of chromosomes supported
	unsigned __int16 MaxDatasets;		// max number of datasets supported
	unsigned __int32 MaxDatasegs;		// max number of datasegs supported
	unsigned __int32 MaxDataPts;
	unsigned __int16 NumChroms;			// actual number of chromosomes
	unsigned __int16 NumDatasets;		// actual number of datasets
	unsigned __int32 NumDatasegs;		// actual number of datasegs
	unsigned __int32 DatasegSize;		// actual maximal sized dataseg
	tsVisDatasetName Datasets[cMaxVisDatasets];	// directory of all datasets names
	tsVisChrom Chroms[cMaxVisChroms];	// directory of all chromosomes
}tsVisHdr;


#pragma pack()

//const int cAlloc4VisBlock =   ((cMaxVisDataPtsBlock * sizeof(tsVisDataPt)) * cMaxVisSpecies) + sizeof(tsVisBlock);


class CVisData
{
	TCHAR m_szFile[_MAX_PATH];		// file containing this instance
	int m_hFile;					// file handle for opened alignment file
	bool m_HdrDirty;				// true if header (m_AlignHdr) needs to be written to disk 
	bool m_SegDirty;				// true if m_pSeg needs to be written to disk
	tsVisHdr m_VisHdr;				// file header
	int m_AllocdSegSize;			// how much memory was allocd for m_pDataseg
	int  m_SegSize;					// how much of the allocd memory is currently being used
	tsDataseg *m_pSeg;				// pts to memory holding a dataseg

public:
	CVisData(void);
	~CVisData(void);

	bool Open(char *pszVisFile,	// specifies file to open or create
			   bool bCreate);		// create file or truncate if it already exists
	bool Close(void);
	bool InitHdr(void);
	bool Reset(bool bFlush);
	bool Flush2Disk(bool BlockOnly);		// flush and commit to disk

	unsigned __int32 StartSeg(void);
	bool WriteSeg(tsDataseg *pSeg);
	bool EndSeg(void);	

	unsigned __int16 GenNameHash(TCHAR *pszName);
	unsigned __int16 GetDatasetID(TCHAR *pszDataset);
	unsigned __int16 GetChromID(TCHAR *pszChromName);

	bool LoadDataseg(INT64 FilePsn,tsDataseg *pSeg,bool bNoDataPts = false);

	bool AddDataseg(TCHAR *pszDataset,			// dataset name
					 TCHAR *pszChromName,
					 unsigned __int32 ChromOfs,	// start on chromosome 
					 unsigned int NumPts,		// number of datapoints in dataset
					 double *pDataPts);			// array of data points

	bool GetDataPoints(CHAR *pszDataset,			// dataset name
						TCHAR *pszChromName,
						unsigned __int32 ChromOfs,	// start on chromosome 
						unsigned int NumPts,		// number of datapoints in dataset
						double MissingMarker,		// value to use if there are missing datapoints
						double *pRetValues);		// where to return datapoint values


};
