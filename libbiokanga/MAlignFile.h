#pragma once
#include "./commdefs.h"

const int cMALGNVersion = 12;			// file header + structure version
const int cMALGNVersionBack = 10;		// backwards compatiable to this version
const int cMaxAlignedSpecies = 100;		// can handle upto this many species in an alignment
const int cMaxAlignedChroms = 0x0ffffff;  	// can handle a total of this many chromosomes or contigs
const int cMaxSrcFiles = 100;		    	// can handle upto this many alignment source files
const int cMaxAlignBlocks =  0x07ffffff;   	// max number of alignment blocks
const int cAllocAlignBlockDirs =  0x0fffff;   // allocate in increments of this many alignment block directory entries
const int cMaxAlignSeqLen =   0x03fffff;	// max length of any aligned sequence per block
const int cMaxMetaDataLen = 1024;		// max total length of any header metadata per source file
const int cMaxParamLen = 4096;			// max total length of any parameter data per source file

const int cChromIDHashSize = 0x08000;		// chromosome hash table size - must be power of 2


// typedefs to make it easy to change common identifier types
typedef INT32 tChromID;			// chromosome identifier
typedef INT32 tAlignBlockID;	// alignment block identifier
typedef UINT8  tSpeciesID;		// species identifier

typedef enum TAG_eMAOpen {
	eMAPOReadOnly = 0,			// open for read only access, updates not allowed
	eMAPOUpdate,				// allow limited number of updates to existing file (currently only confidence scores)
	eMAPOCreate					// create new file or truncate existing
} teMAOpen;

#pragma pack(1)

typedef struct TAG_sAlignSpecies {
	INT32 AlignSpeciesLen;		// total size of this instance
	INT32 ChromID;				// chromosome identifier (species.chrom unique)
	INT32 ChromOfs;				// start offset (0..ChromLen-1) on ChromID (if '-' strand then ChromLen-1..0) 
	INT32 AlignXInDelLen;		// alignment length (1..n) in relative chromosome excluding InDel'-' markers
	INT8  Strand;				// "+" or "-". If the value is "-", the sequence aligns to the reverse-complemented source
	INT8 ConfScore;				// confidence score (0..127) or flags (0..6) for this species in this alignment block
	// note that the actual aligned sequence is concatenated on to each tsAlignSpecies instance
} tsAlignSpecies;

const INT8 cMASeqOverlayFlg = 0x01;	// flag in tsAlignSpecies used to indicate that this sequence is multiple aligned 

typedef struct TAG_sAlignBlock {
	INT32 BlockLenWithSpecies;		// actual total size of this alignment block including all concatenated tsAlignSpecies
	INT32 AlignIncInDelLen;			// alignment sequence length (1..n) , includes '-' insertion markers					
	INT32 AlgnScore;				// alignment block score
	INT32  NumSpecies;				// number of species/chromosomes in represented in this block
	// note that an array of tsAlignSpecies is concatenated on to each tsAlignBlock instance
} tsAlignBlock;

typedef struct TAG_sSpeciesName {
		UINT16 Hash;				// hash used to quickly eliminate species names that can't match in searches
		INT8 SpeciesID;				// globally unique (1..n) species identifier
		INT8 NameLen;				// strlen(szSpeciesName)
	    INT8 szSpeciesName[cMaxDatasetSpeciesChrom]; // species name
} tsSpeciesName;

typedef struct TAG_sChromName {
		INT32 ChromID;				// globally unique (1..n) chromosome identifier
		INT32 NxtHashChromID;		// next chromosome which shares same hash (0 if this is last)
		INT32 ChromLen;				// chromosome length if known - 0 if unknown
		UINT16 Hash;				// hash used to quickly eliminate chromosome names that can't match in searches
		INT8 SpeciesID;				// for this species
		INT8 NameLen;				// strlen(szChromName)
	    INT8 szChromName[cMaxDatasetSpeciesChrom]; // chromosome name
} tsChromName;

typedef struct TAG_sSrcFile {
	INT16 SrcFileID;					// globally unique (1..n), identifies file from which alignments are sourced
	INT16 MetaDataLen;					// header ("##") metadata length
	INT16 AlignParamsLen;				// alignment ("#") parameter length
	INT8 szSrcFile[_MAX_PATH];			// source file path
	INT8 szHeaderMetaData[cMaxMetaDataLen];	// any header ("##") metadata
	INT8 szAlignParams[cMaxParamLen];			// any alignment ("#") parameters
} tsSrcFile;



typedef struct TAG_sBlockDirEl {
	INT64	FileOfs;				// where on disk the associated alignment block starts
	INT32 BlockID;					// block identifer
	INT32 ChromID;					// (sort order 1)reference chromosome
	INT32   ChromOfs;				// (sort order 2) starting offset (0..ChromLen-1)
	INT32  AlignXInDelLen;			// (sort order 3 - longest) reference chromosome block alignment length excluding any InDel '-' (1..n)
} tsBlockDirEl;


typedef struct TAG_sSegCache {
	 char szRelDataset[cMaxDatasetSpeciesChrom];	// which relative dataset
	 unsigned short RelDatasetHash;				// hash on szRelDataset
	 char szRefChrom[cMaxDatasetSpeciesChrom];		// which reference chromosome
	 unsigned short RefChromHash;				// hash on szRefChrom
	 int RefChromOfs;							// start on reference chromosome
	 char szRelChrom[cMaxDatasetSpeciesChrom];		// relative chromosome
	 int RelChromOfs;							// relative chromosome
	 char RelStrand;							// relative strand '+' or '-'
	 int NumPts;								// number of datapoints
	 int AllocdSize;							// max number of data points that can be currently held by pDataPts 
	 unsigned char *pDataPts;					// pts to array of data points
} tsSegCache;

// session specific context
typedef struct TAG_sMAFCtx {
	tChromID m_CurRefChromID;		// aligned block (m_pAlignBlock) reference chromosome identifier
	int m_CurRefChromOfs;			// aligned block (m_pAlignBlock) reference chromosome starting offset 
	int m_CurRefChromXInDelLen;	    // aligned block (m_pAlignBlock) reference chromosome length excluding InDels
	INT64 m_AllocdBlockSize;			// how much memory was allocd for m_pAlignBlock
	INT64 m_CurBlocksSize;			// how much of the allocd memory is currently being used
	tAlignBlockID m_AlignBlockID;	// identifier for current block in m_pAlignBlock
	tsAlignBlock *m_pAlignBlock;	// pts to memory allocated to hold an alignment block
	int m_NumCachedSegs;			// number of cached segments in m_pSegCache
	tsSegCache *m_pSegCache;		// pts to memory allocated to hold segments when joining adjacent segments
} tsMAFCtx;

#pragma pack()

#pragma pack(8)
typedef struct TAG_sAlignHdr {
	UINT8 Magic[4];					// magic chars to identify this file as a biosequence file
	UINT64 FileLen;					// current file length
	INT64 DirElOfs;					// file offset to block directory
	INT64 ChromNamesOfs;			// file offset to ChromNames directory
	UINT32 Type;					// biosequence file type 
	UINT32 Version;					// header version, incremented if structure changes with later releases
	UINT32 SizeOfHdr;				// total size of this header - alignment blocks (sAlignBlocks) immediately follow
	INT32 MaxChroms;				// max number of species.chromosomes supported
	INT32 MaxAlignBlocks;			// max number of alignment blocks supported
	INT32 MaxAlignSeqLen;			// max length of any aligned sequence (including InDels) supported
	INT32 NumChroms;				// actual number of aligned chromosomes
	INT32 NumAlignBlocks;			// actual number of alignment blocks
	INT32 AlignIncInDelLen;			// actual longest aligned sequence (incuding InDels) in any alignment
	INT32 AlignBlockLen;			// actual longest alignment block
	tsSpeciesName SpeciesNames[cMaxAlignedSpecies];	// directory of all aligned species
	INT16 NumSrcFiles;				// actual number of files from which alignments were sourced
	INT16 MaxSrcFiles;				// maximum number of source alignment files supported
	INT8 DataType;					// datatype of sequences
	INT8 NumSpecies;				// actual number of aligned species
	INT8 RefSpeciesID;				// identifer for reference species - all other are relative
	INT8 MaxSpecies;				// max number of species supported
	tsSrcFile SrcFiles[cMaxSrcFiles];			// directory of source files containing sequence alignments
	INT8 szDescription[cMBSFFileDescrLen];		// describes contents of file
	INT8 szTitle[cMBSFShortFileDescrLen];		// short title by which this file can be distingished from other files in dropdown lists etc
}tsAlignHdr;
#pragma pack()
// best initial guestimate of an upper limit on a maximal sized block
// if an actual block is larger then error reported
const INT64 cAlloc4Block = (10 * (sizeof(tsAlignSpecies) + cMaxAlignSeqLen)) + sizeof(tsAlignBlock);

class CMAlignFile : protected CEndian,public CErrorCodes
{
	char m_szFile[_MAX_PATH];		// file containing this instance
	int m_hFile;					// file handle for opened alignment file
	bool m_bHdrDirty;				// true if header (m_FileHdr) needs to be written to disk 
	teMAOpen m_AccessMode;			// eMAPOReadOnly, eMAPOUpdate or eMAPOCreate 
	bool m_BlockStarted;			// true whiltst an alignment block is started until it is closed 
	int m_AlignBlockSeqs;			// cnt of number of sequences associated with alignment block

	tsAlignHdr m_FileHdr;			// alignment header
	INT64 m_AllocdDirElsMem;		// size of memory allocation for m_pDirEls
	int m_NumAllocdDirEls;			// actual many dir elements allocated 
	tsBlockDirEl *m_pDirEls;		// pts to array of block directory elements

	tsChromName *m_pChromNames;		// directory of all aligned species chromosomes


	tChromID m_CurRefChromID;		// aligned block (m_pAlignBlock) reference chromosome identifier
	int m_CurRefChromOfs;			// aligned block (m_pAlignBlock) reference chromosome starting offset 
	int m_CurRefChromXInDelLen;	    // aligned block (m_pAlignBlock) reference chromosome length excluding InDels
	INT64 m_AllocdBlockSize;			// how much memory was allocd for m_pAlignBlock
	INT64 m_CurBlocksSize;			// how much of the allocd memory is currently being used

	tAlignBlockID m_AlignBlockID;	// identifier for current block in m_pAlignBlock
	tsAlignBlock *m_pAlignBlock;	// pts to memory allocated to hold an alignment block
	int m_NumCachedSegs;			// number of cached segments in m_pSegCache
	tsSegCache *m_pSegCache;		// pts to memory allocated to hold segments when joining adjacent segments

	int m_NumFiltSpecies;			// number of species marked to be filtered in m_FiltSpecies
	INT8	m_FiltSpecies[cMaxAlignedSpecies]; // if cMASeqOverlayFlg set then blocks with cMASeqOverlayFlg set for species are skiped by NxtBlock()

	tChromID m_ChromHshTbl[cChromIDHashSize]; // hash table
	
	// ChunkedWrite
	// Seeks to specified 64bit file offset and writes to disk as chunks of no more than INT_MAX/16
	teBSFrsltCodes
		ChunkedWrite(INT64 WrtOfs,UINT8 *pData,INT64 WrtLen);

	// ChunkedRead
	// Seeks to specified 64bit file offset and reads from disk as chunks of no more than INT_MAX/32
	// This allows checks for thread termination requests at reasonable time intervals
	// returns (teBSFrsltCodes)1 if m_bTermThread was set by another thread
	teBSFrsltCodes
		ChunkedRead(INT64 RdOfs,UINT8 *pData,INT64 RdLen);

	teBSFrsltCodes Hdr2Disk(void);	
	teBSFrsltCodes Disk2Hdr(char *pszSeqFile,int FileType);

	tsBlockDirEl *LocateDirEl(tChromID ChromID,int ChromOfs,bool bClosest = false);
	static int CompareBlockDirEls( const void *arg1, const void *arg2 );

public:
	CMAlignFile(void);
	~CMAlignFile(void);
	int Reset(bool bFlush = true);
	int Flush2Disk(bool BlockOnly = true);
	int InitHdr(void);

	int Open(char *pszMAlignFile,
			teMAOpen AccessMode = eMAPOReadOnly,			// eMAPOReadOnly, eMAPOUpdate or eMAPOCreate 
		   char *pszRefSpecies = NULL);	// reference species 

	teBSFrsltCodes SetDescription(char *pszDescription);
	teBSFrsltCodes GetDescription(int MaxLen,char *pszDescription);
	teBSFrsltCodes SetTitle(char *pszTitle);
	teBSFrsltCodes GetTitle(int MaxLen,char *pszTitle);

	
	char *GetRefSpeciesName(void);		// get reference name
	int GetRefSpeciesID(void);			// get reference species identifer

	int Close(bool bWrtDirHdr = true);
	int LocateSpeciesID(char *pszSpeciesName);			// returns species identifer for specified species name
	int LocateChromID(char *pszSpeciesName,char *pszChromName); //returns chromosome identifier
	int LocateChromID(tSpeciesID SpeciesID,char *pszChromName); //returns chromosome identifier

	int StartFile(char *pszSrcFile);		// start processing of specified source file
	int EndFile(void);						// end processing of current source file
	int AddHeader(char*pszHeader);			// adds header line as parsed
	int AddParameters(char *pszParams);		// adds parameter line as parsed
	int StartAlignBlock(long Score);		// starts new alignment block with autogenerated block identier
	int EndAlignBlock(void);				// ends current block
	int WriteBlock(int BlockID,tsAlignBlock *pBlock);
	int	AddAlignSeq(char *pszSpeciesName,	    // species being added to alignment
						 char *pszChromName,	// alignment is on this species chromosome
  						 int ChromLen,			// chromosome length or 0 if unknown
 						 int ChromOfs,			// aligns from this relative offset ('+' == 0..ChromLen-1, '-' == ChromLen-1..0)
						 int IncInDelLen,		// alignment length (1..n) including any InDels
						 char Strand,			// alignment strand
						 char *pszSeq,			// species.chromosome specific alignment sequence
						 bool RptMskUpperCase=false, // true if uppercase represents softmasked repeats
						 INT8 ConfScore=0);		// confidence score or flags associated with this alignment sequence

	 int AddAlignSeqBases(char *pszSpeciesName,	// species being added to alignment
						 char *pszChromName,	// alignment is on species chromosome
 						 int ChromLen,			// chromosome length or 0 if unknown
						 int ChromOfs,			// aligns from this relative offset ('+' == 0..ChromLen-1, '-' == ChromLen-1..0)
					 	 int IncInDelLen,		// alignment length (1..n) in relative chromosome incl InDels
						 char Strand,			// aligns from this strand
						 etSeqBase *pBases,	// species.chromosome specific alignment sequence
						 INT8 ConfScore=0);		// confidence score or flags associated with this alignment sequence

 	int GetNumBlocks(void);				// returns number of blocks in alignment
 	int GetNumChroms(void);				// returns number of chromosomes in alignment

	int LoadBlock(tAlignBlockID BlockID); // load specified block into memory
	int LoadBlock(tsBlockDirEl *pDirEl);  // directory element with block file offset

	
	int NxtBlock(tAlignBlockID CurBlockID);	// returns next alignment block identifier - 0 returns 1st block identifier

	int NxtBlockChrom(tAlignBlockID CurBlockID,tChromID ChromID); // returns next alignment block identifier - 0 returns 1st block identifier
	int NxtChromID(tSpeciesID SpeciesID,tChromID CurChromID);
	tsAlignSpecies *LoadAlignSpecies(tAlignBlockID BlockID,tSpeciesID SpeciesID); //locate and return ptr to alignment for specified block and species
    tsAlignSpecies *LoadAlignChrom(tAlignBlockID BlockID,tChromID ChromID);		// which chromosome to locate alignment for
	
	int GetNumSpecies(tAlignBlockID BlockID = 0);	// returns number of species aligned in specified block ( if 0 then globally)
	int GetNxtBlockSpeciesID(tAlignBlockID BlockID,tSpeciesID SpeciesID);	// returns the next species identifier in specified block
	int GetMaxAlignLen(void);				// returns maximal length alignment (includes '-' indel markers)
	int GetScore(tAlignBlockID BlockID);	// returns score currently associated with the specified block (as processed from orginal MAF file)

	int GetConfScore(tAlignBlockID BlockID,tSpeciesID SpeciesID);		// returns confidence score (0..127) or flags (0..6) for specified block + species
	int SetConfScore(tAlignBlockID BlockID,tSpeciesID SpeciesID,INT8 Score=0,bool bFinal = true);	// sets block.species confidence score (0..127) or flags (0..6), if Final true then commits to disk

	int SetConfFilt(tSpeciesID SpeciesID,bool bFilt); // sets filtering state - cMASeqOverlayFlg - for NxtBlock() filtering for this species
	int SetAllConfFilt(bool bFilt);					  // sets filtering state - cMASeqOverlayFlg - for NxtBlock() filtering for all species
	bool GetConfFilt(tSpeciesID SpeciesID);			  // returns filtering state - cMASeqOverlayFlg - for NxtBlock() filtering for this species 

	char *GetSpeciesName(tSpeciesID SpeciesID);		// returns species name corresponding to species identifier (1..n)
	int GetAlignLen(tAlignBlockID BlockID,tSpeciesID SpeciesID);	// get alignment length (includes '-' indel markers)
	int GetRelChromEndOfs(tAlignBlockID BlockID,tSpeciesID RelSpeciesID);		// get alignment end psn
	char GetStrand(tAlignBlockID BlockID,tSpeciesID SpeciesID);	// get strand - '+' or '-'
	char *GetRelChrom(tAlignBlockID BlockID,tSpeciesID RelSpeciesID);	// get chromosome
	char *GetChromName(tChromID ChromID);  
	int  GetChromLen(tChromID ChromID);			// returns chromosome length  
	int  MapChromOfsToOtherStrand(tChromID ChromID,int ChromOfs);
	int GetRelChromID(tAlignBlockID BlockID,tSpeciesID SpeciesID);	// get chromosome for species
	int GetRelChromOfs(tAlignBlockID CurBlockID,tSpeciesID RelSpeciesID);

	int GetSpeciesID(tChromID ChromID);	// returns species identifier for specified chromosome	

	int GetBlockID(tChromID RefChromID,int ChromOfs);		// returns blockid which contains chromosome.ChromOfs
	etSeqBase *GetSeq(tAlignBlockID BlockID,tSpeciesID SpeciesID);   // returns alignment sequence for specified block and species
	unsigned short GenNameHash(char *pszName);	// generate 16bit hash on specified name
	void RehashSpeciesNames(void);				// rehash species names
	void RehashChromNames(void);				// rehash chromosome names
	int MultiAlignAll2DataPts(char *pszMAF,		// source bioseq mutialignment file
				char *pszDataPointsFile,		// DataPoints file
				char *pszDescr,					// file description
				char *pszTitle);					// file title

	int FlushSegs2Dataset(CDataPoints *pDataPoints);

	int NewSegEntry(tsSegCache *pEntry,
					unsigned short RefChromHash,
					unsigned short RelDatasetHash,
				     char *pszRefChrom,			// reference species chromosome
					 int RefChromOfs,			// reference species chromosome offset
					 char *pszRelDataset,		// into which species relative dataset
					 char *pszRelChrom,			// on which chromosome
					 int RelChromOfs,			// relative species chromosome offset
					 char RelStrand,				// relative strand '+' or '-'
					 int NumPts,				// number of datapoints
					 unsigned char *pDataPts);	// array of data points


	int JoinSegments(CDataPoints *pDataPoints,
				     char *pszRefChrom,			// reference species chromosome
					 int RefChromOfs,			// reference species chromosome offset
					 char *pszRelDataset,		// into which species relative dataset
					 char *pszRelChrom,			// on which chromosome
					 int RelChromOfs,			// relative species chromosome offset 
					 char RelStrand,			// relative strand '+' or '-'
					 int NumPts,				// number of datapoints
					 unsigned char *pDataPts);	// array of data points


	int											// returns number of bases (incl InDels and NotAlignedBase)
		LoadAlignment(char *pszRefChrom,	// alignment is to this ref species chromosome
					  int StartLoci,			// starting from this ref chromosome loci
					  int EndLoci,				// and ending at this ref chromosome loci
					  char *pszRelSpecies,		// aligned bases are to be returned for this species
					  int MaxBases,				// max number of bases to be returned
					  etSeqBase *pRefBases,		// where to return ref alignment bases
					  etSeqBase *pRelBases,		// where to return rel alignment bases
					  etSeqBase RefUnalignedBase = eBaseUndef, // base to return if no ref species alignment
					  etSeqBase RelUnalignedBase = eBaseInDel); // base to return if ref species aligned but no alignment onto rel species

	int												// returns number of bases (incl InDels and NotAlignedBase)
		LoadAlignment(INT32 RefChromID,	// alignment is to this ref species chromosome
					  int StartLoci,			// starting from this ref chromosome loci
					  int EndLoci,				// and ending at this ref chromosome loci
					  int RelSpeciesID,			// aligned bases are to be returned for this species
					  int MaxBases,				// max number of bases to be returned
					  etSeqBase *pRefBases,		// where to return ref alignment bases
					  etSeqBase *pRelBases,		// where to return rel alignment bases
					  etSeqBase RefUnalignedBase = eBaseUndef, // base to return if no ref species alignment
					  etSeqBase RelUnalignedBase = eBaseInDel); // base to return if ref species aligned but no alignment onto rel species

	int											// returns number of bases (incl InDels and NotAlignedBase)
		LoadAlignments(INT32 RefChromID,		// alignment is to this ref species chromosome
					  int StartLoci,			// starting from this ref chromosome loci
					  int EndLoci,				// and ending at this ref chromosome loci
					  int MaxBases,				// max number of bases for each species to be returned
					  int NumSpecies,			// number of species to return alignments for
					  int *pSpeciesIDs,			// pts to array of species identifiers (ref species must be first)
					  int RelOfs,	   		    // offset in pSeqBases[species][] at which to return alignment bases
					  etSeqBase *pSeqBases[],	// pts to array of ptrs of where to return alignment bases for each species			
					  etSeqBase RefUnalignedBase = eBaseUndef, // base to return if no ref species alignment
					  etSeqBase RelUnalignedBase = eBaseInDel); // base to return if ref species aligned but no alignment onto rel species

};
