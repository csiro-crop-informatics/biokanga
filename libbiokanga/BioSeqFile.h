#pragma once
#include "./commdefs.h"

const UINT32 cBSFVersion = 10;					// increment each time the file strucure for bioseq or suffix is changed
const UINT32 cBSFVersionBack = 10;				// minimum version that is this version is still compatible with

const UINT32 cBSFMaxDirEntries = 20000000;	// entry directory can contain upto this many entries
const UINT32 cBSFDataBuffSize     = 0x07ffff; // buffer size used when creating entries
const UINT32 cAllocEntriesIncr = 0x0fffff;	// alloc this size when reallocing directory
const UINT32 cBSFSourceSize = 255;			 // max string size for source of each entry 					
const UINT32 cBSFDescriptionSize = 4095;		 // max string size for describing each entry

const UINT32 cMaxFileReadChunk = 2000000000;	 // read in these chunk sized blocks to overcome 2G read limit

typedef INT32 tBSFEntryID;							// entry identifiers

#pragma pack(1)

typedef enum eDataType {
	eUndefDataType = 0,
	eSeqBaseType,				// data is etSeqBase packed to 4bits/base
	eAsciiType,					// data is 7bit ascii
	eBinDataType,				// data is raw 8bit binary
	eDataPts					// data is visualisation data points
} teDataType;
typedef unsigned char etDataType;

typedef struct TAG_sBSFDirEntry {
	INT64 DataPsn;						    // absolute file psn at which biosequence data starts
	UINT32 Size;							// size in bytes of this instance when concatenated
	INT32 EntryID;							// 1..n uniquely identifies this entry instance
	INT32 NameInst;							// name instance (1..n) -- can have multiple entries with duplicate names
	UINT32 DataLen;							// data length when unpacked
	UINT16 Hash;							// hash on szName
	UINT8 DType:4;							// data type - teDataType
	UINT8 MskSense:1;						// 1 if uppercase represents repeat masked, 0 if lowercase represents repeat masked
	char szName[1];							// place holder for entry name + '\0' + optional description + '\0'
} tsBSFDirEntry;
#pragma pack()

#pragma pack(8)
typedef struct TAG_sBSFHeader {
	unsigned char Magic[4];			 		// magic chars to identify this file as a biosequence file
	INT64 FileLen;							// current file length (write psn for newly created entries)
	INT64 SeqOfs;							// where concatenated sequences start on disk
	INT64 SeqSize;							// disk/memory space required to hold concatenated sequences
	INT64 DirEntriesOfs;					// where directory entries start on disk
	INT64 DirEntryIDIdxOfs;					// where index by EntryID into directory entries start on disk
	INT64 DirEntryNameIdxOfs;				// where index by name into directory entries start on disk
	INT32 Type;								// biosequence file type
	INT32 Version;							// file structure version
	INT32 MaxEntries;						// max number of entries supported 
	INT32 NumEntries;						// actual number of directory entries
	INT32 DirEntriesSize;					// actual disk/memory space required to hold directory entries
	char szDatasetName[cMaxDatasetSpeciesChrom]; // dataset name
	char szDescription[cMBSFFileDescrLen];	// describes contents of file
	char szTitle[cMBSFShortFileDescrLen];	// short title by which this file can be distingished from other files in dropdown lists etc
} tsBSFHeader;

#pragma pack()


class CBioSeqFile : protected CEndian, public CConformation
{
	int m_hFile;							   // opened/created file handle
	char m_szFile[_MAX_PATH+1];				   // file name	as opened/created
	bool m_bHdrDirty;						   // TRUE if current header needs to be written to disk 
	bool m_bCreate;								// TRUE if file opened in create mode
	tsBSFHeader m_FileHdr;					   // file header

	int m_AllocDirEntriesSize;					// disk/memory space required to hold concatenated directory entries
	tsBSFDirEntry *m_pDirEntries;				// pts to  array of tsBSFDirectory's
	tsBSFDirEntry **m_ppDirEntryNames;			// sorted (by name) array of ptrs into m_pDirEntries
	tsBSFDirEntry **m_ppDirEntryIDs;			// sorted (by EntryID) array of ptrs into m_pDirEntries

	tsBSFDirEntry *m_pCreatingEntry;		   // entry currently being created, NULL if none
	UINT32 m_DataBuffLen;				   // number of bytes of data currently in m_pDataBuff
	UINT32 m_DataBuffNibs;			   // number of nibbles currently in m_pDataBuff
	unsigned char *m_pDataBuff;				   // data buffer for use whilst entries are created
	UINT64 m_DataWrtPsn;						// file offset at which to next write contents of m_pDataBuff
	teBSFrsltCodes LoadEntries(void);			// load entries directory into memory
	teBSFrsltCodes ReadDisk(INT64 DiskOfs,int Len,void *pTo); // reads disk block into memory
	teBSFrsltCodes Hdr2Disk(void);			// write header to disk
	teBSFrsltCodes Disk2Hdr(char *pszSeqFile,int FileType); // read header from disk

protected:
	void InitHdr(UINT32 Type);			// initialise file header with biosequence file type
	teBSFrsltCodes Flush2Disk(void);					   // flush and commit to disk
	tsBSFDirEntry *LocateEntry(tBSFEntryID EntryID); // locate entry by EntryID
	teBSFrsltCodes SortEntries(void);
	unsigned short GenNameHash(char *pszName);
	static int SortEntryNames(const void *arg1, const void *arg2); // used to sort by entry name->EntryID


public:
	CBioSeqFile(void);
	virtual ~CBioSeqFile(void);
	teBSFrsltCodes Reset(bool bFlush = true);		// reset state back to that immediately following instantiation before file opened/created
	virtual teBSFrsltCodes Open(char *pszSeqFile,	// specifies file to open or create
    			UINT32 Type = eSeqBaseType,	// biosequence file type 
				bool bCreate = false);				// create file if it doesn't already exist, truncate if it does

	teBSFrsltCodes SetDescription(char *pszDescription);
	teBSFrsltCodes GetDescription(int MaxLen,char *pszDescription);
	teBSFrsltCodes SetTitle(char *pszTitle);
	teBSFrsltCodes GetTitle(int MaxLen,char *pszTitle);

	int GetType(void);
	teBSFrsltCodes Close(void);
	teBSFrsltCodes SetDatasetName(char *pszDataset);	// sets file dataset name
	char *GetDatasetName(void);						// returns file dataset name
	int LocateEntryIDbyName(char *pszName,			// entry name to locate
							int Ith = 1);			// Ith instance to locate (1..n)
	teBSFrsltCodes Exists(tBSFEntryID EntryID);	// returns eBSFSuccess if specified entry exists
	int NumEntries(void);						// returns number of entries
	tBSFEntryID Next(tBSFEntryID Current);		// returns next entry identifer after current
	tBSFEntryID Prev(tBSFEntryID Current);		// returns previous entry identifer before current

	int GetName(tBSFEntryID EntryID,int MaxLen,char *pszName); // get entry name
	int GetDescription(tBSFEntryID EntryID,int MaxLen,char *pDescr); // get entry description
	teBSFrsltCodes GetNameDescription(tBSFEntryID EntryID,int MaxNameLen,char *pszName,int MaxDescrLen,char *pszDescr);
	
	etDataType GetDataType(tBSFEntryID EntryID);
	UINT32 GetDataLen(tBSFEntryID EntryID);		// get entry data length (returns 0 if no associated data or data unavailable)
	int GetData(tBSFEntryID EntryID,etDataType ReqDataType,UINT32 Ofs,unsigned char *pBuffer,UINT32 MaxLen);

	tBSFEntryID CreateEntry(char *pszName, char *pszDescription,etDataType DataType, bool RptMskUpperCase=false);

	teBSFrsltCodes AddData(UINT32 DataLen, unsigned char *pData);		// add data to the currently created entry
	teBSFrsltCodes SealEntry(void);							// seal currently created entry and commits to file - no more data can be added to that entry						

	int GetShuffledBases(int ChromID,	// reference chromosome identifier
					  int ChromOfs,			// start on reference chromosome (0..n) 
					  int NumPts,			// number of points to return
					  unsigned char *pToRet,	// where to return etSeqBases
					  etSeqBase MissingMarker);	// value to use if missing bases

	void ShuffleDinucs(etSeqBase *pSeq,int SeqLen);
	bool Dump2XML(char *pszXMLfile,  UINT32 MaxDumpSeqLen);

	static int  LocateBase(etSeqBase Probe, int TargLen,etSeqBase *pTarget);
	static int LocateSequence(int ProbeLen, etSeqBase *pProbe,int TargLen,etSeqBase *pTarget);

	static int GetNumMasked(int ProbeLen, etSeqBase *pProbe);

	static int
		PackBases(UINT32 SeqLen,	// unpacked sequence length 
		  etSeqBase *pUnpacked,				// pts to sequence to pack from
		  UINT32 NibbleOfs,			// nibble to start pack into (0..SeqLen-1)
		  unsigned char *pPacked);			// where to pack into
	static int
		UnpackBases(UINT32 SeqLen,	// unpacked sequence length 
		  etSeqBase *pUnpacked,				// pts to sequence to unpack into
		  UINT32 NibbleOfs,			// nibble to start unpack from (0..SeqLen-1)
		  unsigned char *pPacked);				// where to unpack from


	static UINT32 ClassComplexity(etSeqBase *pSeq,
							 UINT32 Len,
							 UINT32 ThresPercent);

	static UINT32 ScoreComplexity(etSeqBase *pSeq,
							 UINT32 Len);

	static UINT32								// number of aligned bases
		QuickAlignRight(UINT32 AllowedMismatches,	// total allowed mismatches
		   UINT32 MaxMismatchSeqLen,	// max number of mismatches allowed in any run 
		   UINT32 AllowedProbeInDels,	// total allowed indel events on probe
		   UINT32 AllowedTargInDels,	// total allowed indel events on target
		   UINT32 MaxInDelExtension,	// max length of any InDel extension 
		   UINT32 ProbeLen,			// remaining probe length
		   unsigned char *pProbe,
   		   UINT32 TargLen,			// remaining target length
		   unsigned char *pTarg);

	static UINT32						// number of aligned bases
		QuickAlignLeft(UINT32 AllowedMismatches,	// total allowed mismatches
		   UINT32 MaxMismatchSeqLen,	// max number of mismatches allowed in any run 
		   UINT32 AllowedProbeInDels,	// total allowed indel events on probe
		   UINT32 AllowedTargInDels,	// total allowed indel events on target
		   UINT32 MaxInDelExtension,	// max length of any InDel extension 
		   UINT32 ProbeLen,			// remaining probe length
		   unsigned char *pProbe,
   		   UINT32 TargLen,			// remaining target length
		   unsigned char *pTarg);

	static char * MakeXMLsafe(char *pszStr);
};

