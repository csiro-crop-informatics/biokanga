/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */
#include "stdafx.h"
#if _WIN32
#include <process.h>
#include "../libbiokanga/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libbiokanga/commhdrs.h"
#endif

CMAlignFile::CMAlignFile(void)
{
m_hFile = -1;
m_pAlignBlock = NULL;
m_pChromNames = NULL;
m_AllocdBlockSize = 0;
m_BlockStarted = false;
m_pSegCache = NULL;
m_pDirEls = NULL;
m_AllocdDirElsMem = 0;
m_NumAllocdDirEls = 0;
m_NumCachedSegs = 0;
Reset(false);
}

CMAlignFile::~CMAlignFile(void)
{
if(m_hFile != -1)
	close(m_hFile);

if(m_pDirEls)
	{
#ifdef _WIN32
	free(m_pDirEls);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pDirEls != MAP_FAILED)
		munmap(m_pDirEls,m_AllocdDirElsMem);
#endif
	}

if(m_pAlignBlock)
	delete m_pAlignBlock;
if(m_pChromNames)
	delete m_pChromNames;

if(m_pSegCache != NULL)
	{
	for(int Idx = 0; Idx < cDSFMaxDatasets; Idx++)
		if(m_pSegCache[Idx].pDataPts != NULL)
			delete m_pSegCache[Idx].pDataPts;
	delete m_pSegCache;
	}
}

// Hdr2Disk
// Writes file header out to disk
teBSFrsltCodes 
CMAlignFile::Hdr2Disk(void)		
{
tsAlignHdr FileHdr;
tsAlignHdr *pHdr;
int WrtLen;
int Idx;

WrtLen = sizeof(tsAlignHdr);

if(m_bIsBigEndian)	// if on a big-endian machine then need to make little endian as that is our native file format
	{
	memmove(&FileHdr,&m_FileHdr,WrtLen);
	FileHdr.FileLen = SwapUI64Endians(m_FileHdr.FileLen);					// current file length
	FileHdr.DirElOfs = SwapUI64Endians(m_FileHdr.DirElOfs);					// file offset to block directory
	FileHdr.ChromNamesOfs = SwapUI64Endians(m_FileHdr.ChromNamesOfs);			// file offset to ChromNames directory
	FileHdr.Type = SwapUI32Endians(m_FileHdr.Type);					// biosequence file type 
	FileHdr.Version = SwapUI32Endians(m_FileHdr.Version);					// header version, incremented if structure changes with later releases
	FileHdr.SizeOfHdr = SwapUI32Endians(m_FileHdr.SizeOfHdr);				// total size of this header - alignment blocks (sAlignBlocks) immediately follow
	FileHdr.MaxChroms = SwapUI32Endians(m_FileHdr.MaxChroms);				// max number of species.chromosomes supported
	FileHdr.MaxAlignBlocks = SwapUI32Endians(m_FileHdr.MaxAlignBlocks);			// max number of alignment blocks supported
	FileHdr.MaxAlignSeqLen = SwapUI32Endians(m_FileHdr.MaxAlignSeqLen);			// max length of any aligned sequence (including InDels) supported
	FileHdr.NumChroms = SwapUI32Endians(m_FileHdr.NumChroms);				// actual number of aligned chromosomes
	FileHdr.NumAlignBlocks = SwapUI32Endians(m_FileHdr.NumAlignBlocks);	// actual number of alignment blocks
	FileHdr.AlignIncInDelLen = SwapUI32Endians(m_FileHdr.AlignIncInDelLen);			// actual longest aligned sequence (incuding InDels) in any alignment
	FileHdr.AlignBlockLen = SwapUI32Endians(m_FileHdr.AlignBlockLen);			// actual longest alignment block
	FileHdr.NumSrcFiles = SwapUI16Endians(m_FileHdr.NumSrcFiles);				// actual number of files from which alignments were sourced
	FileHdr.MaxSrcFiles = SwapUI16Endians(m_FileHdr.MaxSrcFiles);				// maximum number of source alignment files supported

	tsSpeciesName *pSpecies;
	pSpecies = FileHdr.SpeciesNames;
	for(Idx = 0; Idx < m_FileHdr.NumSpecies; Idx++, pSpecies++)
		pSpecies->Hash = SwapUI16Endians(pSpecies->Hash);

	tsSrcFile *pSrc;
	pSrc = FileHdr.SrcFiles;
	for(Idx=0; Idx < m_FileHdr.NumSrcFiles; Idx++,pSrc++)
		{
		pSrc->SrcFileID = SwapUI16Endians(pSrc->SrcFileID);					// globally unique (1..n), identifies file from which alignments are sourced
		pSrc->MetaDataLen = SwapUI16Endians(pSrc->MetaDataLen);					// header ("##") metadata length
		pSrc->AlignParamsLen = SwapUI16Endians(pSrc->AlignParamsLen);				// alignment ("#") parameter length
		}

	pHdr = &FileHdr;
	}
else
	pHdr = &m_FileHdr;
if(_lseeki64(m_hFile,0,SEEK_SET) ||
		write(m_hFile,pHdr,WrtLen)!=WrtLen)
	{
	AddErrMsg("CDataPoints::Flush2Disk","Unable to write file header to disk on file %s - error %s",m_szFile,strerror(errno));
	Reset(false);
	return(eBSFerrFileAccess);
	}
m_bHdrDirty = false;
return(eBSFSuccess);
}

// Disk2Hdr
// Reads file header from disk
teBSFrsltCodes 
CMAlignFile::Disk2Hdr(char *pszSeqFile,int FileType)		
{
int Idx;

if(_lseeki64(m_hFile,0,SEEK_SET)!=0)			// read in header..
	{
	AddErrMsg("CMAlignFile::Disk2Hdr","Seek failed to offset 0 on %s - %s",pszSeqFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

if(sizeof(tsAlignHdr) != read(m_hFile,&m_FileHdr,sizeof(tsAlignHdr)))
	{
	AddErrMsg("CMAlignFile::Disk2Hdr","Read of file header failed on %s - %s",pszSeqFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

// header read, validate it as being a sequence file header
if(tolower(m_FileHdr.Magic[0]) != 'b' || 
	tolower(m_FileHdr.Magic[1]) != 'i' || 
	tolower(m_FileHdr.Magic[2]) != 'o' || 
	tolower(m_FileHdr.Magic[3]) != 's')
	{
	AddErrMsg("CMAlignFile::Disk2Hdr","%s opened but no magic signature - not a bioseq file",pszSeqFile);
	Reset(false);			// closes opened file..
	return(eBSFerrNotBioseq);
	}


if(m_bIsBigEndian)	// file was written with little-endian ordering
	{
	m_FileHdr.FileLen = SwapUI64Endians(m_FileHdr.FileLen);					// current file length
	m_FileHdr.DirElOfs = SwapUI64Endians(m_FileHdr.DirElOfs);					// file offset to block directory
	m_FileHdr.ChromNamesOfs = SwapUI64Endians(m_FileHdr.ChromNamesOfs);			// file offset to ChromNames directory
	m_FileHdr.Type = SwapUI32Endians(m_FileHdr.Type);					// biosequence file type 
	m_FileHdr.Version = SwapUI32Endians(m_FileHdr.Version);					// header version, incremented if structure changes with later releases
	m_FileHdr.SizeOfHdr = SwapUI32Endians(m_FileHdr.SizeOfHdr);				// total size of this header - alignment blocks (sAlignBlocks) immediately follow
	m_FileHdr.MaxChroms = SwapUI32Endians(m_FileHdr.MaxChroms);				// max number of species.chromosomes supported
	m_FileHdr.MaxAlignBlocks = SwapUI32Endians(m_FileHdr.MaxAlignBlocks);			// max number of alignment blocks supported
	m_FileHdr.MaxAlignSeqLen = SwapUI32Endians(m_FileHdr.MaxAlignSeqLen);			// max length of any aligned sequence (including InDels) supported
	m_FileHdr.NumChroms = SwapUI32Endians(m_FileHdr.NumChroms);				// actual number of aligned chromosomes
	m_FileHdr.NumAlignBlocks = SwapUI32Endians(m_FileHdr.NumAlignBlocks);	// actual number of alignment blocks
	m_FileHdr.AlignIncInDelLen = SwapUI32Endians(m_FileHdr.AlignIncInDelLen);			// actual longest aligned sequence (incuding InDels) in any alignment
	m_FileHdr.AlignBlockLen = SwapUI32Endians(m_FileHdr.AlignBlockLen);			// actual longest alignment block
	m_FileHdr.NumSrcFiles = SwapUI16Endians(m_FileHdr.NumSrcFiles);				// actual number of files from which alignments were sourced
	m_FileHdr.MaxSrcFiles = SwapUI16Endians(m_FileHdr.MaxSrcFiles);				// maximum number of source alignment files supported

	tsSpeciesName *pSpecies;
	pSpecies = m_FileHdr.SpeciesNames;
	for(Idx = 0; Idx < m_FileHdr.NumSpecies; Idx++, pSpecies++)
		pSpecies->Hash = SwapUI16Endians(pSpecies->Hash);

	tsSrcFile *pSrc;
	pSrc = m_FileHdr.SrcFiles;
	for(Idx=0; Idx < m_FileHdr.NumSrcFiles; Idx++,pSrc++)
		{
		pSrc->SrcFileID = SwapUI16Endians(pSrc->SrcFileID);					// globally unique (1..n), identifies file from which alignments are sourced
		pSrc->MetaDataLen = SwapUI16Endians(pSrc->MetaDataLen);					// header ("##") metadata length
		pSrc->AlignParamsLen = SwapUI16Endians(pSrc->AlignParamsLen);				// alignment ("#") parameter length
		}
	}

	// check alignment file is the type we are expecting
if(FileType != cBSFTypeAny && m_FileHdr.Type != FileType)
	{
	AddErrMsg("CMAlignFile::Disk2Hdr","%s opened as a multialignment file - expected type %d, file type is %d",pszSeqFile,FileType,m_FileHdr.Type);
	Reset(false);			// closes opened file..
	return(eBSFerrFileType);
	}

	// can we handle this version?
if(m_FileHdr.Version < cMALGNVersionBack || m_FileHdr.Version > cMALGNVersion)
	{
	AddErrMsg("CMAlignFile::Disk2Hdr","%s opened as a multialignment file - expected between version %d and %d, file version is %d",pszSeqFile,
			cBSFVersionBack,cBSFVersion,m_FileHdr.Version);
	Reset(false);			// closes opened file..
	return(eBSFerrFileVer);
	}


m_bHdrDirty = false;
return(eBSFSuccess);
}


// Reset()
// Resets alignment file class context back to that immediately following instantiation
int
CMAlignFile::Reset(bool bFlush)		// true (default) is to write any pending header writes to disk before closing opened file
{
if(m_hFile != -1)
	{
	if(bFlush)
		Flush2Disk(false);
	close(m_hFile);
	m_hFile = -1;
	}
m_BlockStarted = false;
m_szFile[0] = '\0';
InitHdr();
m_bHdrDirty = false;
m_AlignBlockSeqs = 0;

if(m_pAlignBlock)
	{
	delete m_pAlignBlock;
	m_pAlignBlock = NULL;
	}
m_AllocdBlockSize = 0;
m_AlignBlockID = 0;

if(m_pChromNames != NULL)
	{
	delete m_pChromNames;
	m_pChromNames = NULL;
	}

if(m_pSegCache != NULL)
	{
	for(int Idx = 0; Idx < cDSFMaxDatasets; Idx++)
		if(m_pSegCache[Idx].pDataPts != NULL)
			delete m_pSegCache[Idx].pDataPts;
	delete m_pSegCache;
	m_pSegCache = NULL;
	}

m_NumFiltSpecies = 0;
memset(m_FiltSpecies,0,sizeof(m_FiltSpecies));

memset(m_ChromHshTbl,0,sizeof(m_ChromHshTbl));

return(eBSFSuccess);
}

int
CMAlignFile::InitHdr(void)
{
memset(&m_FileHdr,0,sizeof(tsAlignHdr));
m_FileHdr.Magic[0] = 'b';
m_FileHdr.Magic[1] = 'i';
m_FileHdr.Magic[2] = 'o';
m_FileHdr.Magic[3] = 's';
m_FileHdr.Type = cBSFTypeMultiAlign;		 // it's a multialignment file
m_FileHdr.Version = cMALGNVersion;		     // file structure version
m_FileHdr.FileLen = sizeof(tsAlignHdr);	     // current file length
m_FileHdr.SizeOfHdr = sizeof(tsAlignHdr);    // size of this header
m_FileHdr.DataType = eSeqBaseType;			 // store sequences as enumerated bases
m_FileHdr.MaxSrcFiles = cMaxSrcFiles;	     // maximum number of source alignment files supported
m_FileHdr.MaxSpecies=cMaxAlignedSpecies;	 // max number of species supported
m_FileHdr.MaxChroms=cMaxAlignedChroms;		 // max number of species.chromosomes supported
m_FileHdr.MaxAlignBlocks = cMaxAlignBlocks;  // max number of alignment blocks supported
m_FileHdr.MaxAlignSeqLen = cMaxAlignSeqLen;  // max aligned seq length supported
m_FileHdr.szDescription[0] = '\0';
m_FileHdr.szTitle[0] = '\0';
m_bHdrDirty = true;
return(eBSFSuccess);
}

// Close
// Close file, and if bWrtDirHdr==true then sort and write out block directory header
// NOTE: After sorting a check is made for blocks which overlap other blocks as AXT/MAF files
//		 can contain overlapping blocks which are separated in the source AXT/MAF files by other
//		 blocks which do not overlap. These blocks can only be identified once the complete AXT/MAF
//	     file has been processed.
//
//		 
int
CMAlignFile::Close(bool bWrtDirHdr) // if true (default) then write out block directory 
{
int Rslt;
int WrtLen;
int Idx;
int NumAlignBlocks;
tsBlockDirEl *pDirEl;
tsBlockDirEl *pPrvDirEl;


if(m_hFile != -1)
	{
	if(bWrtDirHdr && m_AccessMode == eMAPOCreate)
		{

		if(m_pDirEls != NULL && m_FileHdr.NumAlignBlocks)
			{
			pPrvDirEl = m_pDirEls;
			if(m_FileHdr.NumAlignBlocks == 1)
				pPrvDirEl->BlockID = 1;
			else
				{
				// sort directory by ref species chrom->ofs
				qsort(m_pDirEls,m_FileHdr.NumAlignBlocks,sizeof(tsBlockDirEl),CompareBlockDirEls);
			
				// now try and identify any overlap blocks
				pDirEl = &m_pDirEls[1];
				pPrvDirEl->BlockID = NumAlignBlocks = 1;
				for(Idx = 2; Idx <= (int)m_FileHdr.NumAlignBlocks; Idx++,pDirEl++)
					{	
					if(pDirEl->ChromID == pPrvDirEl->ChromID &&				
						pDirEl->ChromOfs < (pPrvDirEl->ChromOfs + pPrvDirEl->AlignXInDelLen))	// have overlap?
						continue;							// slough this block
					// this block can be kept as it doesn't overlap
					pPrvDirEl += 1;
					if(pPrvDirEl != pDirEl)
						*pPrvDirEl = *pDirEl;
					pPrvDirEl->BlockID = ++NumAlignBlocks;		
					}
				m_FileHdr.NumAlignBlocks = NumAlignBlocks;
				}

			if(m_bIsBigEndian)
				{
				pDirEl = m_pDirEls;
				for(Idx = 0; Idx < m_FileHdr.NumAlignBlocks; Idx++, pDirEl++)
					{
					pDirEl->FileOfs=SwapUI64Endians(pDirEl->FileOfs); // where on disk the associated alignment block starts
					pDirEl->BlockID=SwapUI32Endians(pDirEl->BlockID);			// block identifer
					pDirEl->ChromID=SwapUI32Endians(pDirEl->ChromID);				// (sort order 1)reference chromosome
					pDirEl->ChromOfs=SwapUI32Endians(pDirEl->ChromOfs);				// (sort order 2) starting offset (0..ChromLen-1)
					pDirEl->AlignXInDelLen=SwapUI32Endians(pDirEl->AlignXInDelLen);			// (sort order 3 - longest) reference chromosome block alignment length excluding any InDel '-' (1..n)
					}
				}

			WrtLen = m_FileHdr.NumAlignBlocks * sizeof(tsBlockDirEl);	
			if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
				write(m_hFile,m_pDirEls,WrtLen)!=WrtLen)
				{
				close(m_hFile);
				m_hFile = -1;
				return(eBSFerrFileAccess);
				}
			m_FileHdr.DirElOfs = m_FileHdr.FileLen;
			m_FileHdr.FileLen += WrtLen;
			}
		else
			m_FileHdr.NumAlignBlocks = 0;


		if(m_pChromNames != NULL && m_FileHdr.NumChroms)
			{
			if(m_bIsBigEndian)
				{
				tsChromName *pChrom = m_pChromNames;
				for(Idx = 0; Idx < m_FileHdr.NumChroms; Idx++, pChrom++)
					{
					pChrom->ChromID = SwapUI32Endians(pChrom->ChromID); // globally unique (1..n) chromosome identifier
					pChrom->NxtHashChromID = SwapUI32Endians(pChrom->NxtHashChromID);		// next chromosome which shares same hash (0 if this is last)
					pChrom->ChromLen = SwapUI32Endians(pChrom->ChromLen);					// chromosome length if known - 0 if unknown
					pChrom->Hash = SwapUI16Endians(pChrom->Hash);					// hash used to quickly eliminate chromosome names that can't match in searches
					}
				}

			WrtLen = m_FileHdr.NumChroms * sizeof(tsChromName);	
			if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
				write(m_hFile,m_pChromNames,WrtLen)!=WrtLen)
				{
				close(m_hFile);
				m_hFile = -1;
				return(eBSFerrFileAccess);
				}
			m_FileHdr.ChromNamesOfs = m_FileHdr.FileLen;
			m_FileHdr.FileLen += WrtLen;
			}
		else
			{
			m_FileHdr.ChromNamesOfs = 0;
			m_FileHdr.NumChroms = 0;
			}

		if((Rslt=Hdr2Disk())!= eBSFSuccess)
			return(Rslt);
		}
	close(m_hFile);
	m_hFile = -1;

	}
m_szFile[0] = '\0';
if(m_pDirEls)
	{
#ifdef _WIN32
	free(m_pDirEls);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pDirEls != MAP_FAILED)
		munmap(m_pDirEls,m_AllocdDirElsMem);
#endif
	m_AllocdDirElsMem = 0;
	m_NumAllocdDirEls = 0;
	m_pDirEls = NULL;
	}

InitHdr();
m_AccessMode = eMAPOReadOnly;
return(eBSFSuccess);
}

// ChunkedWrite
// Seeks to specified 64bit file offset and writes to disk as chunks of no more than INT_MAX/16
teBSFrsltCodes
CMAlignFile::ChunkedWrite(INT64 WrtOfs,UINT8 *pData,INT64 WrtLen)
{
int BlockLen;
if(_lseeki64(m_hFile,WrtOfs,SEEK_SET) != WrtOfs)
	{
	AddErrMsg("CMAlignFile::ChunkedWrite","Unable to seek to %ld on file %s - error %s",WrtOfs,m_szFile,strerror(errno));
	Reset(false);
	return(eBSFerrFileAccess);
	}

while(WrtLen)
	{
	BlockLen = WrtLen > (INT64)(INT_MAX/16) ? (INT_MAX/16) : (int)WrtLen;
	WrtLen -= BlockLen;
	if(write(m_hFile,pData,BlockLen)!=BlockLen)
		{
		AddErrMsg("CMAlignFile::ChunkedWrite","Unable to write to disk on file %s - error %s",m_szFile,strerror(errno));
		Reset(false);
		return(eBSFerrFileAccess);
		}
	pData += BlockLen;
	}
return(eBSFSuccess);
}

// ChunkedRead
// Seeks to specified 64bit file offset and reads from disk as chunks of no more than INT_MAX/32
teBSFrsltCodes
CMAlignFile::ChunkedRead(INT64 RdOfs,UINT8 *pData,INT64 RdLen)
{
UINT32 BlockLen;
UINT32 ActualRdLen;
if(_lseeki64(m_hFile,RdOfs,SEEK_SET) != RdOfs)
	{
	AddErrMsg("CSfxArrayV3::ChunkedRead","Unable to seek to %ld on file %s - error %s",RdOfs,m_szFile,strerror(errno));
	return(eBSFerrFileAccess);
	}

while(RdLen)
	{
	BlockLen = RdLen > (INT64)(INT_MAX/2) ? (INT_MAX/2) : (UINT32)RdLen;
	RdLen -= BlockLen;
	if((ActualRdLen = read(m_hFile,pData,BlockLen))!=BlockLen)
		{
		AddErrMsg("CMAlignFile::ChunkedRead","Unable to read from disk file %s - error %s",m_szFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	pData += BlockLen;
	}
return(eBSFSuccess);
}


// Open()
// Open specified multiple alignment file pszMAlignFile
// Option to create or truncate pszMAlignFile
int 
CMAlignFile::Open(char *pszMAlignFile,	// specifies file to open or create
			  teMAOpen AccessMode,	// eMAPOReadOnly, eMAPOUpdate or eMAPOCreate 
			   char *pszRefSpeciesName)	// reference species 
{
int Rslt;
INT64 AllocLen;
int ReadLen;
INT64 FileOfs;

Close(false);						// reset context in case file still opened

if(pszMAlignFile == NULL || *pszMAlignFile == '\0') // validate parameters
	return(eBSFerrParams);
if(m_AccessMode == eMAPOCreate && (pszRefSpeciesName == NULL || pszRefSpeciesName[0] == '\0'))
	return(eBSFerrParams);

switch(AccessMode) {
#ifdef _WIN32
	case eMAPOReadOnly:				// read only, no updates
		m_hFile = open(pszMAlignFile,O_READSEQ ); // file access is normally sequential..
		break;
	case eMAPOUpdate:				// mostly reads but could also be updates
		m_hFile = open(pszMAlignFile,O_READORWRITERAND); 
		break;
	case eMAPOCreate:				// create/truncate
		m_hFile = open(pszMAlignFile, O_CREATETRUNC );
		break;
#else
	case eMAPOReadOnly:				// read only, no updates
		m_hFile = open64(pszMAlignFile,O_READSEQ); // file access is normally sequential..
		break;
	case eMAPOUpdate:				// mostly reads but could also be updates
		m_hFile = open(pszMAlignFile,O_READORWRITERAND); 
		break;
	case eMAPOCreate:				// create/truncate
		if((m_hFile = open64(pszMAlignFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	        if(ftruncate(m_hFile,0)!=0)
				{
				AddErrMsg("CMAlignFile::Open","Unable to truncate %s - %s",pszMAlignFile,strerror(errno));
				Rslt = eBSFerrCreateFile;
				return(Rslt);
				}
		break;
#endif
	}

if(m_hFile == -1)					// check if file open succeeded
	{
	AddErrMsg("CMAlignFile::Open","Unable to open %s - %s",pszMAlignFile,strerror(errno));
	Rslt = AccessMode == eMAPOCreate ? eBSFerrCreateFile : eBSFerrOpnFile;
	Reset(false);
	return(Rslt);
	}

strncpy(m_szFile,pszMAlignFile,_MAX_PATH);
m_szFile[_MAX_PATH-1] = '\0';

if(AccessMode == eMAPOCreate)
	{
	InitHdr();
	// initialise with specified reference species
	tsSpeciesName *pSpecies = &m_FileHdr.SpeciesNames[0];
	m_FileHdr.RefSpeciesID = pSpecies->SpeciesID = m_FileHdr.NumSpecies = 1;
	strncpy((char *)pSpecies->szSpeciesName,pszRefSpeciesName,cMaxDatasetSpeciesChrom);
	pSpecies->szSpeciesName[cMaxDatasetSpeciesChrom-1] = '\0';
	pSpecies->Hash = GenNameHash(pszRefSpeciesName);
	m_AlignBlockSeqs = 0;
	Flush2Disk(false);
		// allocate memory to hold an initially sized block directory, will be extended as may be required
	m_NumAllocdDirEls = cAllocAlignBlockDirs;

	AllocLen = m_NumAllocdDirEls * sizeof(tsBlockDirEl);
	m_AllocdDirElsMem = AllocLen;

#ifdef _WIN32
	m_pDirEls = (tsBlockDirEl *) malloc((size_t)m_AllocdDirElsMem);	// initial and will be realloc'd as may be required for additional blocks
	if(m_pDirEls == NULL)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes through malloc()  failed - %s",(INT64)m_AllocdDirElsMem,strerror(errno));
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pDirEls = (tsBlockDirEl *)mmap(NULL,(size_t)m_AllocdDirElsMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pDirEls == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)m_AllocdDirElsMem,strerror(errno));
		m_pDirEls = NULL;
		}
#endif
	if(m_pDirEls == NULL)
		{
		m_AllocdDirElsMem = 0;
		Close();
		return(eBSFerrMem);
		}
	m_pDirEls[0].BlockID = 0;
	m_pDirEls[0].ChromID = eBSFerrChrom;

	// allocate memory to hold maximal number of chrom names
	AllocLen = m_FileHdr.MaxChroms * sizeof(tsChromName);	
	if((m_pChromNames = (tsChromName *)new unsigned char [(int)AllocLen])==NULL)
		{
		Close();
		return(eBSFerrMem);
		}
	m_AccessMode = AccessMode;
	}
else // else opening existing file
	{
	if((Rslt = Disk2Hdr(pszMAlignFile,cBSFTypeMultiAlign))!=eBSFSuccess)
		{
		Reset(false);			// closes opened file..
		return(Rslt);
		}

	// ensure alignment file has a reference species
	if(m_FileHdr.RefSpeciesID == 0)
		{
		AddErrMsg("CMAlignFile::Open","expected reference species to be defined, %s has none...",pszMAlignFile);
		Reset(false);
		return(eBSFerrDataset);
		}

	// ensure contains at least one alignment block
	if(m_FileHdr.NumAlignBlocks < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Alignment file opened but contains no alignment blocks");
		Close();
		return(eBSFerrNumAlignBlks);
		}

		// allocate memory to hold block directory + guard element
	AllocLen = (m_FileHdr.NumAlignBlocks + 1) * sizeof(tsBlockDirEl);
	m_AllocdDirElsMem = AllocLen;
	m_NumAllocdDirEls = m_FileHdr.NumAlignBlocks + 1;

#ifdef _WIN32
	m_pDirEls = (tsBlockDirEl *) malloc((size_t)m_AllocdDirElsMem);	// initial and expected to be the only allocation
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pDirEls = (tsBlockDirEl *)mmap(NULL,(size_t)m_AllocdDirElsMem, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pDirEls == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)m_AllocdDirElsMem,strerror(errno));
		m_pDirEls = NULL;
		}
#endif
	if(m_pDirEls == NULL)
		{
		m_AllocdDirElsMem = 0;
		Close();
		return(eBSFerrMem);
		}

	AllocLen -= sizeof(tsBlockDirEl);
	if((Rslt=ChunkedRead(m_FileHdr.DirElOfs,(UINT8 *)m_pDirEls,AllocLen))!=eBSFSuccess)
		{
		Close();
		return(Rslt);
		}

    m_pDirEls[m_FileHdr.NumAlignBlocks].BlockID = 0;			// initial guard
	m_pDirEls[m_FileHdr.NumAlignBlocks].ChromID = eBSFerrChrom;
	if(m_bIsBigEndian)
		{
		int Idx;
		tsBlockDirEl *pDirEl;
		pDirEl = m_pDirEls;
		for(Idx = 0; Idx < m_FileHdr.NumAlignBlocks; Idx++, pDirEl++)
			{
			pDirEl->FileOfs=SwapUI64Endians(pDirEl->FileOfs); // where on disk the associated alignment block starts
			pDirEl->BlockID=SwapUI32Endians(pDirEl->BlockID);			// block identifer
			pDirEl->ChromID=SwapUI32Endians(pDirEl->ChromID);				// (sort order 1)reference chromosome
			pDirEl->ChromOfs=SwapUI32Endians(pDirEl->ChromOfs);				// (sort order 2) starting offset (0..ChromLen-1)
			pDirEl->AlignXInDelLen=SwapUI32Endians(pDirEl->AlignXInDelLen);			// (sort order 3 - longest) reference chromosome block alignment length excluding any InDel '-' (1..n)
			}
		}

			// allocate memory to hold chrom names
	AllocLen = m_FileHdr.NumChroms * sizeof(tsChromName);	
	if((m_pChromNames = (tsChromName *)new unsigned char [(int)AllocLen])==NULL)
		{
		Close();
		return(eBSFerrMem);
		}

	if((FileOfs=_lseeki64(m_hFile,m_FileHdr.ChromNamesOfs,SEEK_SET)) != m_FileHdr.ChromNamesOfs)
		{
		Close();
		return(eBSFerrFileAccess);
		}
	if((ReadLen = read(m_hFile,m_pChromNames,(unsigned int)AllocLen)) != AllocLen)
		{
		Close();
		return(eBSFerrFileAccess);
		}
	if(m_bIsBigEndian)
		{
		int Idx;
		tsChromName *pChrom = m_pChromNames;
		for(Idx = 0; Idx < m_FileHdr.NumChroms; Idx++, pChrom++)
			{
			pChrom->ChromID = SwapUI32Endians(pChrom->ChromID); // globally unique (1..n) chromosome identifier
			pChrom->NxtHashChromID = SwapUI32Endians(pChrom->NxtHashChromID);		// next chromosome which shares same hash (0 if this is last)
			pChrom->ChromLen = SwapUI32Endians(pChrom->ChromLen);					// chromosome length if known - 0 if unknown
			pChrom->Hash = SwapUI16Endians(pChrom->Hash);					// hash used to quickly eliminate chromosome names that can't match in searches
			}
		}


	// rehash in case hash name function has changed
	RehashSpeciesNames();
	RehashChromNames();
	m_AccessMode = AccessMode;
	}
return(eBSFSuccess);
}


teBSFrsltCodes 
CMAlignFile::SetDescription(char *pszDescription)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(m_AccessMode != eMAPOCreate)
	return(eBSFerrRead);
if(m_FileHdr.Version >= 2)
	{
	strncpy((char *)m_FileHdr.szDescription,pszDescription,sizeof(m_FileHdr.szDescription));
	m_FileHdr.szDescription[sizeof(m_FileHdr.szDescription)-1] = '\0';
	return(eBSFSuccess);
	}
return(eBSFerrFileVer);
}

teBSFrsltCodes 
CMAlignFile::GetDescription(int MaxLen,char *pszDescription)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(m_AccessMode == eMAPOCreate)
	return(eBSFerrWrite);
if(m_FileHdr.Version >= 2)
	strncpy(pszDescription,(char *)m_FileHdr.szDescription,MaxLen);
else
	strncpy(pszDescription,m_szFile,MaxLen);
pszDescription[MaxLen-1] = '\0';
return(eBSFSuccess);
}

teBSFrsltCodes 
CMAlignFile::SetTitle(char *pszTitle)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(m_AccessMode != eMAPOCreate)
	return(eBSFerrRead);
if(m_FileHdr.Version >= 2)
	{
	strncpy((char *)m_FileHdr.szTitle,pszTitle,sizeof(m_FileHdr.szTitle));
	m_FileHdr.szTitle[sizeof(m_FileHdr.szTitle)-1] = '\0';
	return(eBSFSuccess);
	}
return(eBSFerrFileVer);
}

teBSFrsltCodes 
CMAlignFile::GetTitle(int MaxLen,char *pszTitle)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(m_AccessMode == eMAPOCreate)
	return(eBSFerrWrite);
if(m_FileHdr.Version >= 2)
	strncpy(pszTitle,(char *)m_FileHdr.szTitle,MaxLen);
else
	{
	char szFname[_MAX_FNAME];
#ifdef _WIN32
	_splitpath(m_szFile,NULL,NULL,szFname,NULL);
#else
	CUtility::splitpath(m_szFile,NULL,szFname);
#endif
	}
pszTitle[MaxLen-1] = '\0';
return(eBSFSuccess);
}


char *
CMAlignFile::GetRefSpeciesName(void)		// get reference name
{
if(m_hFile == -1)
	{
	AddErrMsg("CMAlignFile::GetRefSpeciesName","file not opened");
	return(NULL);
	}
if(m_FileHdr.RefSpeciesID < 1 || m_FileHdr.RefSpeciesID > m_FileHdr.NumSpecies)
	return(NULL);
return((char *)m_FileHdr.SpeciesNames[m_FileHdr.RefSpeciesID-1].szSpeciesName);
}

int 
CMAlignFile::GetRefSpeciesID(void)			// get reference species identifer
{
if(m_hFile == -1)
	{
	AddErrMsg("CMAlignFile::GetRefSpeciesID","file not opened");
	return(eBSFerrClosed);
	}
if(m_FileHdr.RefSpeciesID < 1 || m_FileHdr.RefSpeciesID > m_FileHdr.NumSpecies)
	return(eBSFerrDataset);
return(m_FileHdr.RefSpeciesID);
}

int
CMAlignFile::Flush2Disk(bool BlockOnly)						// flush and commit to disk
{
int Rslt;
if(m_hFile == -1)
	{
	AddErrMsg("CMAlignFile::Flush2Disk","file not opened");
	return(eBSFerrClosed);
	}

if(m_BlockStarted && m_AlignBlockSeqs)									    // does block contain at least one sequence?
	{
	Rslt = WriteBlock(m_AlignBlockID,m_pAlignBlock);
	if(Rslt != eBSFSuccess)
		{
		AddErrMsg("CMAlignFile::Flush2Disk","WriteBlock(BlockID=%d) failed - %s",m_AlignBlockID,CErrorCodes::ErrText((teBSFrsltCodes)Rslt));
		return(Rslt);
		}
	m_BlockStarted = false;
	m_pAlignBlock->BlockLenWithSpecies = sizeof(tsAlignBlock);
	m_pAlignBlock->NumSpecies = 0;
	m_pAlignBlock->AlignIncInDelLen = 0;
	m_pAlignBlock->AlgnScore = 0;
	m_bHdrDirty = true;
	m_AlignBlockSeqs = 0;
	m_CurRefChromID = 0;		
	m_CurRefChromOfs = 0;		
	m_CurRefChromXInDelLen = 0;
	m_AlignBlockID = 0;
	}

if(m_bHdrDirty && !BlockOnly)
	{
	if((Rslt=Hdr2Disk())!=eBSFSuccess)
		return(Rslt);
	m_bHdrDirty = false;
	}

return(eBSFSuccess);
}

// StartFile
// Start parsing of new source alignment file
int	 
CMAlignFile::StartFile(char *pszSrcFile)
{
tsSrcFile *pSrcFile;
if(m_FileHdr.NumSrcFiles == cMaxSrcFiles)
	return(eBSFerrNumSrcFiles);
pSrcFile = &m_FileHdr.SrcFiles[m_FileHdr.NumSrcFiles++];
pSrcFile->SrcFileID=m_FileHdr.NumSrcFiles;
strncpy((char *)pSrcFile->szSrcFile,pszSrcFile,_MAX_PATH);
pSrcFile->szSrcFile[_MAX_PATH-1] = '\0';
pSrcFile->AlignParamsLen = 0;
pSrcFile->MetaDataLen = 0;
pSrcFile->szAlignParams[0] = '\0';
pSrcFile->szHeaderMetaData[0]= '\0';
m_bHdrDirty = true;
return(eBSFSuccess);
}

// EndFile
// End processing of current source alignment file
int 
CMAlignFile::EndFile(void) // no real need for this function except that it balances StartFile -:)
{
return(eBSFSuccess);
}

int
CMAlignFile::AddHeader(char *pszHeaderLine)		// adds header line as parsed
{
int Len;
int MaxLen;
tsSrcFile *pSrcFile;
if(!m_FileHdr.NumSrcFiles)
	return(eBSFerrNoSrcFiles);
Len = (int)strlen(pszHeaderLine);
pSrcFile = &m_FileHdr.SrcFiles[m_FileHdr.NumSrcFiles-1];

MaxLen = cMaxMetaDataLen - pSrcFile->MetaDataLen;
if(MaxLen > 10)
	{
	if(pSrcFile->MetaDataLen)
		{
		pSrcFile->szHeaderMetaData[pSrcFile->MetaDataLen++] = ' ';
		MaxLen -= 1;
		}
	strncpy((char *)&pSrcFile->szHeaderMetaData[pSrcFile->MetaDataLen],pszHeaderLine,MaxLen);
	pSrcFile->szHeaderMetaData[cMaxMetaDataLen-1] = '\0';
	}
pSrcFile->MetaDataLen = (INT16)strlen((char *)pSrcFile->szHeaderMetaData);
m_bHdrDirty = true;
return(eBSFSuccess);
}

int 
CMAlignFile::AddParameters(char *pszParamsLine)	// adds parameter line as parsed
{
int Len;
int MaxLen;
tsSrcFile *pSrcFile;
if(!m_FileHdr.NumSrcFiles)
	return(eBSFerrNoSrcFiles);
Len = (int)strlen(pszParamsLine);
pSrcFile = &m_FileHdr.SrcFiles[m_FileHdr.NumSrcFiles-1];

MaxLen = cMaxParamLen - pSrcFile->AlignParamsLen;
if(MaxLen > 10)
	{
	if(pSrcFile->AlignParamsLen)
		{
		pSrcFile->szAlignParams[pSrcFile->AlignParamsLen++] = ' ';
		MaxLen -= 1;
		}
	strncpy((char *)&pSrcFile->szAlignParams[pSrcFile->AlignParamsLen],pszParamsLine,MaxLen);
	pSrcFile->szAlignParams[cMaxParamLen-1] = '\0';
	}
pSrcFile->AlignParamsLen = (INT16)strlen((char *)pSrcFile->szAlignParams);
m_bHdrDirty = true;
return(eBSFSuccess);
}

// Some heuristic rules are applied due to some artifacts in the UCSC MAF and AXF files -
// a) Can have alignments of a single sequence
// b) Can have multiple alignments to same locii
// c) Can have alignment blocks which overlap other blocks
// Rules are -
// a) If block closed with a single sequence then block is simply reused - single sequences are sloughed
// b) If block closed has a ChromOfs == previous block ChromOfs then block simply reused - new block is sloughed
// c) Overlaping blocks are detected when directory header is sorted just prior to file being closed
// Because these heuristics are rather arbitary they are documented here as well as elsewhere in this file
// Note also that if an alignment block is closed = EndAlignBlock() = without any intevening call to AddAlignSeq() then
// the block is not written to disk and the alignment block identifier is reused
int											// returns alignment block identifier or error code
CMAlignFile::StartAlignBlock(long Score)	// starts new alignment block with autogenerated block identifier
{
int Rslt;
if(m_FileHdr.NumAlignBlocks >= cMaxAlignBlocks)
	return(eBSFerrNumAlignBlks);

if(m_pAlignBlock == NULL)
	{
	if(NULL == (m_pAlignBlock = (tsAlignBlock *)new unsigned char[cAlloc4Block]))
		return(eBSFerrMem);
	m_AllocdBlockSize = cAlloc4Block;
	m_AlignBlockID = 0;
	}
else
	if(m_BlockStarted)			// write out any previously existing block
		{
		if(!m_AlignBlockSeqs)		// if block started but no associated sequences then reuse this block
			{								// including the block identifier
			m_AlignBlockID = 0;
			m_BlockStarted = false;
			m_FileHdr.NumAlignBlocks--;		// block is being sloughed!
			m_AlignBlockID = 0;
			}
		else
			{
			Rslt = WriteBlock(m_AlignBlockID,m_pAlignBlock);
			if(Rslt != eBSFSuccess)
				{
				AddErrMsg("CMAlignFile::StartAlignBlock","WriteBlock(BlockID=%d) failed - %s",m_AlignBlockID,CErrorCodes::ErrText((teBSFrsltCodes)Rslt));
				return(Rslt);
				}
			}
		}

m_pAlignBlock->AlgnScore = Score;
m_pAlignBlock->BlockLenWithSpecies = sizeof(tsAlignBlock);
m_pAlignBlock->NumSpecies = 0;
m_CurRefChromID = (tChromID)eBSFerrChrom;
m_BlockStarted = true;
m_AlignBlockSeqs = 0;
m_AlignBlockID = ++m_FileHdr.NumAlignBlocks;
return(m_AlignBlockID);
}

int
CMAlignFile::EndAlignBlock(void)			// ends current block
{
int Rslt;
if(!m_BlockStarted)							// was a prior call to StartAlignBlock() made?
	return(eBSFSuccess);					// treat as success as there was no alignment block processing

if(!m_AlignBlockSeqs)				// if block started but no associated sequences then reuse this block
	{								// including the block identifier
	m_AlignBlockID = 0;
	m_BlockStarted = false;
	m_FileHdr.NumAlignBlocks--;		// block is being sloughed!
	return(eBSFSuccess);			
	}

Rslt = WriteBlock(m_AlignBlockID,m_pAlignBlock);
m_AlignBlockID = 0;
m_BlockStarted = false;
return(Rslt);
}

// WriteBlock
// adds specified block to block directory then writes block to file
// then resets block
// Some heuristic rules are applied due to some artifacts in the UCSC MAF files -
// a) Have observed block alignments of a single sequence
// b) Have observed multiple block alignments to same locii
// c) Have observed blocks which are contained within other blocks
// Rules are -
// a) If block closed with a single alignment then block is simply reused - single alignments are sloughed
// b) If block closed is contained within previous block then block simply reused - new block is sloughed
// Note that these heuristics are rather arbitary!
int
CMAlignFile::WriteBlock(int BlockID,tsAlignBlock *pBlock) // block containing all species alignments to reference
{
tsBlockDirEl *pDirEl;
tsAlignBlock Block;
int SpeciesIdx;
tsAlignSpecies *pSpecies;
UINT8 *pByte;

if(pBlock == NULL || !pBlock->BlockLenWithSpecies || !m_BlockStarted)
	return(eBSFerrAlignBlk);

if(m_CurRefChromID == (tChromID)eBSFerrChrom)
	return(eBSFerrAlignBlk);

if((BlockID + 10) > m_NumAllocdDirEls)		// need to realloc the block directory?
	{
	size_t memreq = ((m_NumAllocdDirEls + cAllocAlignBlockDirs) * sizeof(tsBlockDirEl));
#ifdef _WIN32
	pDirEl = (tsBlockDirEl *) realloc(m_pDirEls,memreq);
	if(pDirEl == NULL)
		{
#else
	pDirEl = (tsBlockDirEl *)mremap(m_pDirEls,m_AllocdDirElsMem,memreq,MREMAP_MAYMOVE);
	if(pDirEl == MAP_FAILED)
		{
#endif
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"WriteBlock: Memory for block directory reallocation to %lld bytes failed - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pDirEls = pDirEl;
	m_AllocdDirElsMem = memreq;
	m_NumAllocdDirEls += cAllocAlignBlockDirs;
	}

// check to see if this block overlaps as previous block or if it is a single species block
// if so then block will be sloughed
if(BlockID > 1)
	pDirEl = &m_pDirEls[BlockID-2];
else
	pDirEl = NULL;

if((pDirEl != NULL  && 
   pDirEl->ChromOfs <= m_CurRefChromOfs && 
   (pDirEl->AlignXInDelLen + pDirEl->ChromOfs) >=  (m_CurRefChromOfs + m_CurRefChromXInDelLen)) ||
		pBlock->NumSpecies == 1)
	{
	pBlock->BlockLenWithSpecies = sizeof(tsAlignBlock);
	pBlock->NumSpecies = 0;
	pBlock->AlignIncInDelLen = 0;
	pBlock->AlgnScore = 0;
	m_BlockStarted = false;
	m_FileHdr.NumAlignBlocks--;		// block is being sloughed!
	m_AlignBlockID = 0;
	m_AlignBlockSeqs = 0;
	return(eBSFSuccess);			// too many artifacts in MAF/AXT files to treat them as fatal or even warning errors!
	}

// initialise block directory element
pDirEl = &m_pDirEls[BlockID-1];
pDirEl->ChromID = m_CurRefChromID;
pDirEl->ChromOfs = m_CurRefChromOfs;
pDirEl->AlignXInDelLen = m_CurRefChromXInDelLen;

pDirEl->BlockID = BlockID;
pDirEl->FileOfs= m_FileHdr.FileLen;

if(m_bIsBigEndian)
	{
	memmove(&Block,pBlock,sizeof(tsAlignBlock));
	pBlock->BlockLenWithSpecies=SwapUI32Endians(pBlock->BlockLenWithSpecies);// actual total size of this alignment block including all concatenated tsAlignSpecies
	pBlock->AlignIncInDelLen=SwapUI32Endians(pBlock->AlignIncInDelLen);		// alignment sequence length (1..n) , includes '-' insertion markers					
	pBlock->AlgnScore=SwapUI32Endians(pBlock->AlgnScore);					// alignment block score
	pBlock->NumSpecies=SwapUI32Endians(pBlock->NumSpecies);					// number of species/chromosomes in represented in this block
	pSpecies = (tsAlignSpecies *)(((UINT8 *)pBlock) + sizeof(tsAlignBlock));
	pByte = (UINT8 *)pSpecies;
	for(SpeciesIdx = 0; SpeciesIdx < Block.NumSpecies; SpeciesIdx++)
		{
		pByte += pSpecies->AlignSpeciesLen;
		pSpecies->AlignSpeciesLen=SwapUI32Endians(pSpecies->AlignSpeciesLen);	// total size of this instance
		pSpecies->ChromID=SwapUI32Endians(pSpecies->ChromID);					// chromosome identifier (species.chrom unique)
		pSpecies->ChromOfs=SwapUI32Endians(pSpecies->ChromOfs);					// start offset (0..ChromLen-1) on ChromID (if '-' strand then ChromLen-1..0) 
		pSpecies->AlignXInDelLen=SwapUI32Endians(pSpecies->AlignXInDelLen);		// alignment length (1..n) in relative chromosome excluding InDel'-' markers
		pSpecies = (tsAlignSpecies *)pByte;
		}

	if(_lseeki64(m_hFile,pDirEl->FileOfs,SEEK_SET) != pDirEl->FileOfs ||			
		write(m_hFile,pBlock,Block.BlockLenWithSpecies)!=Block.BlockLenWithSpecies)
		{
		Close();			// closes opened file..
		return(eBSFerrFileAccess);
		}

	memmove(pBlock,&Block,sizeof(tsAlignBlock));
	pSpecies = (tsAlignSpecies *)(((UINT8 *)pBlock) + sizeof(tsAlignBlock));
	pByte = (UINT8 *)pSpecies;
	for(SpeciesIdx = 0; SpeciesIdx < pBlock->NumSpecies; SpeciesIdx++)
		{
		pSpecies = (tsAlignSpecies *)pByte;
		pSpecies->AlignSpeciesLen=SwapUI32Endians(pSpecies->AlignSpeciesLen);	// total size of this instance
		pSpecies->ChromID=SwapUI32Endians(pSpecies->ChromID);					// chromosome identifier (species.chrom unique)
		pSpecies->ChromOfs=SwapUI32Endians(pSpecies->ChromOfs);					// start offset (0..ChromLen-1) on ChromID (if '-' strand then ChromLen-1..0) 
		pSpecies->AlignXInDelLen=SwapUI32Endians(pSpecies->AlignXInDelLen);		// alignment length (1..n) in relative chromosome excluding InDel'-' markers
		pByte += pSpecies->AlignSpeciesLen;
		}
	}
else
	{
	if(_lseeki64(m_hFile,pDirEl->FileOfs,SEEK_SET) != pDirEl->FileOfs ||			
		write(m_hFile,pBlock,pBlock->BlockLenWithSpecies)!=pBlock->BlockLenWithSpecies)
		{
		Close();			// closes opened file..
		return(eBSFerrFileAccess);
		}
	}

m_FileHdr.FileLen += pBlock->BlockLenWithSpecies;
pDirEl+=1;					// initialise new end of array marker
pDirEl->BlockID = 0;
pDirEl->ChromOfs = 0;
pDirEl->AlignXInDelLen = 0;
pDirEl->ChromID = eBSFerrChrom;

if(m_FileHdr.AlignBlockLen < pBlock->BlockLenWithSpecies)
	m_FileHdr.AlignBlockLen = pBlock->BlockLenWithSpecies;
if(m_FileHdr.AlignIncInDelLen < pBlock->AlignIncInDelLen)
	m_FileHdr.AlignIncInDelLen = pBlock->AlignIncInDelLen;
pBlock->BlockLenWithSpecies = sizeof(tsAlignBlock);
pBlock->NumSpecies = 0;
pBlock->AlignIncInDelLen = 0;
pBlock->AlgnScore = 0;
m_bHdrDirty = true;
m_BlockStarted = false;
m_AlignBlockSeqs = 0;
return(eBSFSuccess);
}

// GenNameHash
// Generates a 16bit hash on specified lowercased name
// This hash can then be used to quickly eliminate probe names which can't match a target name by comparing hashes
UINT16 
CMAlignFile::GenNameHash(char *pszName)
{
unsigned long hash = 5381;
char Chr;
while (Chr = *pszName++)
	hash = ((hash << 5) + hash) + tolower(Chr);
return ((UINT16)hash);
}


// RehashSpeciesNames
// Rehashes species names - can be useful if name hash function has been updated!
void
CMAlignFile::RehashSpeciesNames(void)
{
tsSpeciesName *pSpecies;
int Idx;
pSpecies = &m_FileHdr.SpeciesNames[0];
for(Idx = 0; Idx < m_FileHdr.NumSpecies; Idx++, pSpecies++)
	pSpecies->Hash = GenNameHash((char *)pSpecies->szSpeciesName);
}

// RehashChromNames
void 
CMAlignFile::RehashChromNames(void)
{
tsChromName *pChrom;
tChromID HashID;
int Idx;
if(m_pChromNames == NULL)
	return;
pChrom = m_pChromNames;
for(Idx = 0; Idx < (int)m_FileHdr.NumChroms; Idx++,pChrom++)
	{
	pChrom->Hash = GenNameHash((char *)pChrom->szChromName);
	HashID = m_ChromHshTbl[pChrom->Hash & (cChromIDHashSize-1)];
	pChrom->NxtHashChromID = HashID;
	m_ChromHshTbl[pChrom->Hash & (cChromIDHashSize-1)] = pChrom->ChromID;
	}
}

// LocateSpeciesID
// Returns SpeciesID which is used to identify the species from it's name
// eBSFerrDataset is returned if unable to locate species
// Uses linear search but shouldn't impact on performance as limited number of species
int
CMAlignFile::LocateSpeciesID(char *pszSpeciesName)
{
tsSpeciesName *pSpecies;
int Idx;
unsigned short Hash;
Hash = GenNameHash(pszSpeciesName);
pSpecies = &m_FileHdr.SpeciesNames[0];
for(Idx = 0; Idx < m_FileHdr.NumSpecies; Idx++, pSpecies++)
	{
	if(pSpecies->Hash == Hash &&
		!stricmp((char *)pSpecies->szSpeciesName,pszSpeciesName))
		return(pSpecies->SpeciesID);
	}
return(eBSFerrDataset);
}

//returns chromosome identifier
int
CMAlignFile::LocateChromID(char *pszSpeciesName,char *pszChromName)
{
int SpeciesID;
int Rslt;
if((Rslt = SpeciesID=LocateSpeciesID(pszSpeciesName)) <= eBSFSuccess)
	return(Rslt);
return(LocateChromID(SpeciesID,pszChromName));
}


int
CMAlignFile::SetConfScore(tAlignBlockID BlockID,tSpeciesID SpeciesID,INT8 ConfScore,bool bFinal)
{
tsAlignSpecies *pSpecies;

if(!SpeciesID || SpeciesID > m_FileHdr.NumSpecies || !BlockID || BlockID > m_FileHdr.NumAlignBlocks)
	return(eBSFerrParams);

if((pSpecies = LoadAlignSpecies(BlockID,SpeciesID))==NULL)
	return(eBSFerrSpecies);		// assume it's due to the Species not in specified alignment block

pSpecies->ConfScore = ConfScore;
if(!bFinal)
	return(eBSFSuccess);

tsBlockDirEl *pDirEl = &m_pDirEls[BlockID-1];
if(_lseeki64(m_hFile,pDirEl->FileOfs,SEEK_SET) != pDirEl->FileOfs ||			
	write(m_hFile,m_pAlignBlock,m_pAlignBlock->BlockLenWithSpecies)!=m_pAlignBlock->BlockLenWithSpecies)
		{
		Close();			// closes opened file..
		return(eBSFerrFileAccess);
		}
return(eBSFSuccess);
}

// returns chromosome identifier
int 
CMAlignFile::LocateChromID(tSpeciesID SpeciesID,char *pszChromName)
{
tChromID HashID;
unsigned short Hash;
tsChromName *pChrom;

if(SpeciesID < 1 || SpeciesID > m_FileHdr.NumSpecies || m_pChromNames == NULL)
	return(eBSFerrDataset);

Hash = GenNameHash(pszChromName);
HashID = m_ChromHshTbl[Hash  & (cChromIDHashSize-1)];
while(HashID != 0)
	{
	pChrom = &m_pChromNames[HashID-1]; 
	if(pChrom->Hash == Hash &&
		pChrom->SpeciesID == SpeciesID &&
		!stricmp((char *)pChrom->szChromName,pszChromName))
			return(pChrom->ChromID);
	HashID = pChrom->NxtHashChromID;
	}
return(eBSFerrChrom);
}



// NxtChromID
// Returns the next allocated ChromID for the specified SpeciesID
// If CurChromID == 0 then the returned ChromID is the lowest allocated ChromID
// Will return eBSFerrChrom when all ChromIDs have been returned
int
CMAlignFile::NxtChromID(tSpeciesID SpeciesID,
						tChromID CurChromID)
{
int Idx;
int NxtID = INT_MAX;
tsChromName *pChrom = m_pChromNames;
for(Idx = 0; Idx < m_FileHdr.NumChroms; Idx++,pChrom++)
	{
	if(pChrom->SpeciesID == SpeciesID)
		{
		if(CurChromID < pChrom->ChromID &&
			pChrom->ChromID < NxtID)
				NxtID = pChrom->ChromID;
		}
	}
return(NxtID == INT_MAX ? eBSFerrChrom : NxtID);
}

// AddAlignSeq
// Add an aligned sequence into current block
// Some heuristic rules are applied due to some artifacts in the UCSC MAF files -
// a) Can have alignments of a single sequence
// b) Can have multiple alignments to same locii
// Rules are -
// a) If block closed with a single aligment then block is simply reused - single alignment is sloughed
// b) If block closed has a ChromOfs == previous block ChromOfs then block simply reused - new block is sloughed
// Because these heuristics are rather arbitary they are documented here 
int	
CMAlignFile::AddAlignSeq(char *pszSpeciesName,	    // species being added to alignment
						 char *pszChromName,	// alignment is on this species chromosome
  						 int ChromLen,			// chromosome length or 0 if unknown
 						 int ChromOfs,			// aligns from this relative offset ('+' == 0..ChromLen-1, '-' == ChromLen-1..0)
						 int IncInDelLen,		// alignment length (1..n) including any InDels
						 char Strand,			// alignment strand
						 char *pszSeq,			// species.chromosome specific alignment sequence
						 bool RptMskUpperCase,  // true if uppercase represents softmasked repeats
 						 INT8 ConfScore)		// confidence associated with this alignment sequence
{
static etSeqBase Ascii2SeqBuff[0x07fff];
etSeqBase *pSeq;
int Rslt;
if(IncInDelLen >= 0x07fff)	// need to use dynamically allocated buffer?
	{
	pSeq = new unsigned char [IncInDelLen+1];
	CSeqTrans::MapAscii2Sense(pszSeq,IncInDelLen,pSeq,RptMskUpperCase);
	Rslt = AddAlignSeqBases(pszSpeciesName,pszChromName,ChromLen,ChromOfs,IncInDelLen,Strand,pSeq,ConfScore);
	delete pSeq;
	}
else
	{
	pSeq = Ascii2SeqBuff;
	CSeqTrans::MapAscii2Sense(pszSeq,IncInDelLen,pSeq,RptMskUpperCase);
	Rslt = AddAlignSeqBases(pszSpeciesName,pszChromName,ChromLen,ChromOfs,IncInDelLen,Strand,pSeq,ConfScore);
	}
return(Rslt);
}

int 
CMAlignFile::AddAlignSeqBases(char *pszSpeciesName,	// species being added to alignment
						 char *pszChromName,	// alignment is on species chromosome
 						 int ChromLen,			// chromosome length or 0 if unknown
						 int ChromOfs,			// aligns from this relative offset ('+' == 0..ChromLen-1, '-' == ChromLen-1..0)
					 	 int IncInDelLen,		// alignment length (1..n) in relative chromosome incl InDels
						 char Strand,			// aligns from this strand
						 etSeqBase *pBases,		// species.chromosome specific alignment sequence
						 INT8 ConfScore)		// confidence associated with this alignment sequence
{
tsSpeciesName *pSpecies;
tsChromName *pChrom;
tsAlignSpecies *pAlignSpecies;
tSpeciesID SpeciesID;
tChromID ChromID;
int XInDelLen;
etSeqBase *pSeq;
int Cnt;

if(!m_BlockStarted || m_pAlignBlock == NULL ||
   pszSpeciesName == NULL || pszSpeciesName[0] == '\0' ||
   pszChromName == NULL || pszChromName[0] == '\0' ||
   ChromOfs < 0 ||
   IncInDelLen <= 0 ||
   !(Strand == '+' || Strand == '-') ||
   pBases == NULL)
	return(eBSFerrParams);

// determine length without InDels
XInDelLen = 0;
pSeq = pBases;
for(Cnt=0;Cnt < IncInDelLen; Cnt++)
	if(*pSeq++ != eBaseInDel)
		XInDelLen++;

SpeciesID = (tSpeciesID)LocateSpeciesID(pszSpeciesName);
if(SpeciesID == (tSpeciesID)eBSFerrDataset)		// if a new species
	{
	if(m_FileHdr.NumSpecies >= cMaxAlignedSpecies)
		return(eBSFerrMaxDatasets);
	pSpecies = &m_FileHdr.SpeciesNames[m_FileHdr.NumSpecies++];
	memset(pSpecies,0,sizeof(tsSpeciesName));
	SpeciesID = pSpecies->SpeciesID = m_FileHdr.NumSpecies;
	strncpy((char *)pSpecies->szSpeciesName,pszSpeciesName,cMaxDatasetSpeciesChrom);
	pSpecies->szSpeciesName[cMaxDatasetSpeciesChrom-1] = '\0';
	pSpecies->Hash = GenNameHash(pszSpeciesName);
	}

ChromID = LocateChromID(SpeciesID,pszChromName);
if(ChromID == (tChromID)eBSFerrChrom)			// if a new chromosome
	{
	if(m_FileHdr.NumChroms >= cMaxAlignedChroms)
		{
		AddErrMsg("CMAlignFile::AddAlignSeqBases","Too many chromosomes, can handle at most %d",cMaxAlignedChroms);
		return(eBSFerrMaxChroms);
		}
	pChrom = &m_pChromNames[m_FileHdr.NumChroms++];
	memset(pChrom,0,sizeof(tsChromName));
	ChromID = pChrom->ChromID = m_FileHdr.NumChroms;
	strncpy((char *)pChrom->szChromName,pszChromName,cMaxDatasetSpeciesChrom);
	pChrom->szChromName[cMaxDatasetSpeciesChrom-1] = '\0';
	pChrom->Hash = GenNameHash(pszChromName);
	pChrom->SpeciesID = SpeciesID;
	pChrom->ChromLen = ChromLen;

	pChrom->NxtHashChromID = m_ChromHshTbl[pChrom->Hash & (cChromIDHashSize-1)];
	m_ChromHshTbl[pChrom->Hash & (cChromIDHashSize-1)] = pChrom->ChromID;
	}
else
	pChrom = &m_pChromNames[ChromID-1];

if(SpeciesID == m_FileHdr.RefSpeciesID)
	{
	m_CurRefChromID = ChromID;
	m_CurRefChromOfs = ChromOfs;
	m_CurRefChromXInDelLen = XInDelLen;
	}

pAlignSpecies = (tsAlignSpecies *)((unsigned char *)m_pAlignBlock + m_pAlignBlock->BlockLenWithSpecies);
pAlignSpecies->ChromID  = ChromID;
pAlignSpecies->AlignXInDelLen = XInDelLen;
pAlignSpecies->ChromOfs = ChromOfs;
pAlignSpecies->Strand = Strand;
pAlignSpecies->ConfScore = ConfScore;
pSeq = (etSeqBase *)pAlignSpecies + sizeof(tsAlignSpecies);
memmove(pSeq,pBases,IncInDelLen);

pAlignSpecies->AlignSpeciesLen = sizeof(tsAlignSpecies) + IncInDelLen;
m_pAlignBlock->BlockLenWithSpecies += pAlignSpecies->AlignSpeciesLen;
m_pAlignBlock->AlignIncInDelLen = IncInDelLen;
m_pAlignBlock->NumSpecies+=1;

m_AlignBlockSeqs++;			// at least one seq in this alignment block
m_bHdrDirty = true;
return(eBSFSuccess);
}


// load specified alignment block into memory
int
CMAlignFile::LoadBlock(tAlignBlockID BlockID)	// if directory element not known
{
tsBlockDirEl *pDirEl;

// check if block exists...
if(BlockID < 1 || !m_FileHdr.NumAlignBlocks || BlockID > m_FileHdr.NumAlignBlocks)
	return(eBSFerrAlignBlk);

// check if block already loaded...
if(m_pAlignBlock != NULL && m_AlignBlockID == BlockID)
	return(eBSFSuccess);

pDirEl = &m_pDirEls[BlockID-1];
return(LoadBlock(pDirEl));
}

int
CMAlignFile::LoadBlock(tsBlockDirEl *pDirEl)// directory element with block file offset
{
UINT32 BlockRemaining;
int SpeciesIdx;
tsAlignSpecies *pSpecies;
UINT8 *pByte;

if(pDirEl->FileOfs < m_FileHdr.SizeOfHdr)	// should never have a block which starts before the header finishes!
	return(eBSFerrInternal);				// but bugs or incomplete block writes may result in incorrect offset

	// if memory not already allocated to hold block then allocate sufficent memory
	// to hold largest block in file
if(m_pAlignBlock == NULL)
	{
	if(NULL == (m_pAlignBlock = (tsAlignBlock *)new unsigned char[m_FileHdr.AlignBlockLen]))
		return(eBSFerrMem);
	memset(m_pAlignBlock,0,sizeof(tsAlignBlock));
	m_AllocdBlockSize = m_FileHdr.AlignBlockLen;
	m_AlignBlockID = 0;
	}
else
	// check if block already loaded...
	if(m_AlignBlockID == pDirEl->BlockID)
		return(eBSFSuccess);
m_AlignBlockID = 0;
if(pDirEl->FileOfs != _lseeki64(m_hFile,pDirEl->FileOfs,SEEK_SET))
	return(eBSFerrFileAccess);
if(sizeof(tsAlignBlock) != read(m_hFile,m_pAlignBlock,sizeof(tsAlignBlock)))
	return(eBSFerrFileAccess);	

if(m_bIsBigEndian)
	{
	m_pAlignBlock->BlockLenWithSpecies=SwapUI32Endians(m_pAlignBlock->BlockLenWithSpecies);	// actual total size of this alignment block including all concatenated tsAlignSpecies
	m_pAlignBlock->AlignIncInDelLen=SwapUI32Endians(m_pAlignBlock->AlignIncInDelLen);		// alignment sequence length (1..n) , includes '-' insertion markers					
	m_pAlignBlock->AlgnScore=SwapUI32Endians(m_pAlignBlock->AlgnScore);						// alignment block score
	m_pAlignBlock->NumSpecies=SwapUI32Endians(m_pAlignBlock->NumSpecies);					// number of species/chromosomes in represented in this block
	}

m_AlignBlockID = pDirEl->BlockID;
// fixed size header has been read, read in balance of block
BlockRemaining = m_pAlignBlock->BlockLenWithSpecies - sizeof(tsAlignBlock);
if(BlockRemaining)
	{
	if(BlockRemaining != read(m_hFile,(unsigned char *)m_pAlignBlock + sizeof(tsAlignBlock),BlockRemaining))
		return(eBSFerrFileAccess);
	}

if(m_bIsBigEndian)
	{
	pSpecies = (tsAlignSpecies *)(((UINT8 *)m_pAlignBlock) + sizeof(tsAlignBlock));
	pByte = (UINT8 *)pSpecies;
	for(SpeciesIdx = 0; SpeciesIdx < m_pAlignBlock->NumSpecies; SpeciesIdx++)
		{
		pSpecies = (tsAlignSpecies *)pByte;
		pSpecies->AlignSpeciesLen=SwapUI32Endians(pSpecies->AlignSpeciesLen);	// total size of this instance
		pSpecies->ChromID=SwapUI32Endians(pSpecies->ChromID);					// chromosome identifier (species.chrom unique)
		pSpecies->ChromOfs=SwapUI32Endians(pSpecies->ChromOfs);					// start offset (0..ChromLen-1) on ChromID (if '-' strand then ChromLen-1..0) 
		pSpecies->AlignXInDelLen=SwapUI32Endians(pSpecies->AlignXInDelLen);		// alignment length (1..n) in relative chromosome excluding InDel'-' markers
		pByte += pSpecies->AlignSpeciesLen;
		}
	}
return(eBSFSuccess);
}

tsAlignSpecies *										// NULL if species not in block alignment
CMAlignFile::LoadAlignSpecies(tAlignBlockID BlockID,	// which alignment block
							  tSpeciesID SpeciesID)		// which species to locate alignment for
{
tsAlignSpecies *pSpecies;
int NumSpecies;
tChromID ChromID;

if(!SpeciesID || SpeciesID > m_FileHdr.NumSpecies || !BlockID || BlockID > m_FileHdr.NumAlignBlocks)
	return(NULL);

if(LoadBlock(BlockID)!=eBSFSuccess)
	return(NULL);

// check if species is represented in the alignment block
NumSpecies = m_pAlignBlock->NumSpecies;
if(!NumSpecies)
	return(NULL);

pSpecies = (tsAlignSpecies *)((unsigned char *)m_pAlignBlock + sizeof(tsAlignBlock));
while(ChromID = pSpecies->ChromID) 
	{
	if(SpeciesID == m_pChromNames[ChromID-1].SpeciesID)
		break;
	if(!--NumSpecies)
		return(NULL);
	pSpecies = (tsAlignSpecies *)((unsigned char *)pSpecies + pSpecies->AlignSpeciesLen);
	}
return(pSpecies);
}

tsAlignSpecies *										// NULL if chromosome not in block alignment
CMAlignFile::LoadAlignChrom(tAlignBlockID BlockID,		// which alignment block
							tChromID ChromID)			// which chromosome to locate alignment for
{
tsAlignSpecies *pSpecies;
UINT32 NumSpecies;

if(!ChromID || ChromID > m_FileHdr.NumChroms || !BlockID || BlockID > m_FileHdr.NumAlignBlocks)
	return(NULL);

if(LoadBlock(BlockID)!=eBSFSuccess)
	return(NULL);

// check if chromosome is represented in the alignment block
NumSpecies = m_pAlignBlock->NumSpecies;
if(!NumSpecies)
	return(NULL);

pSpecies = (tsAlignSpecies *)((unsigned char *)m_pAlignBlock + sizeof(tsAlignBlock));
while(ChromID != pSpecies->ChromID) 
	{
	if(!--NumSpecies)
		return(NULL);
	pSpecies = (tsAlignSpecies *)((unsigned char *)pSpecies + pSpecies->AlignSpeciesLen);
	}
return(pSpecies);
}

// returns score currently associated with the specified block (as processed from orginal MAF file)
int 
CMAlignFile::GetScore(tAlignBlockID BlockID)
{
int Rslt;
if((Rslt=LoadBlock(BlockID))!=eBSFSuccess)
	return(Rslt);
return(m_pAlignBlock->AlgnScore);
}


// returns confidence score currently associated with the specified block and species
int 
CMAlignFile::GetConfScore(tAlignBlockID BlockID,tSpeciesID SpeciesID)
{
tsAlignSpecies *pSpecies;

if(!SpeciesID || SpeciesID > m_FileHdr.NumSpecies || !BlockID || BlockID > m_FileHdr.NumAlignBlocks)
	return(eBSFerrParams);

if((pSpecies = LoadAlignSpecies(BlockID,SpeciesID))==NULL)
	return(eBSFerrSpecies);		// assume it's due to the Species not in specified alignment block
return(pSpecies->ConfScore);
}

// sets filtering state - cMASeqOverlayFlg - for NxtBlock() filtering for this species
int 
CMAlignFile::SetConfFilt(tSpeciesID SpeciesID,bool bFilt)
{
if(!SpeciesID || SpeciesID > m_FileHdr.NumSpecies)
	return(eBSFerrParams);

if(bFilt)
	{
	if(m_FiltSpecies[SpeciesID-1] & cMASeqOverlayFlg)	// if already set then no further action required
		return(eBSFSuccess);
	m_FiltSpecies[SpeciesID-1] |= cMASeqOverlayFlg;
	m_NumFiltSpecies += 1;
	}
else
	{
	if(!(m_FiltSpecies[SpeciesID-1] & cMASeqOverlayFlg))	// if already reset then no further action required
		return(eBSFSuccess);

	m_FiltSpecies[SpeciesID-1] &= ~cMASeqOverlayFlg;
	m_NumFiltSpecies -= 1;
	}
return(eBSFSuccess);
}

// sets filtering state - cMASeqOverlayFlg - for NxtBlock() filtering for all species
int 
CMAlignFile::SetAllConfFilt(bool bFilt)
{
int SpeciesIdx;
for(SpeciesIdx = 0; SpeciesIdx < m_FileHdr.NumSpecies; SpeciesIdx++)
	{
	if(bFilt)
		m_FiltSpecies[SpeciesIdx] |= cMASeqOverlayFlg;
	else
		m_FiltSpecies[SpeciesIdx] &= ~cMASeqOverlayFlg;
	}
if(bFilt)
	m_NumFiltSpecies = m_FileHdr.NumSpecies;
else
	m_NumFiltSpecies = 0;
return(eBSFSuccess);
}

// returns filtering state - cMASeqOverlayFlg - for NxtBlock() filtering for this species 
bool 
CMAlignFile::GetConfFilt(tSpeciesID SpeciesID)
{
if(!SpeciesID || SpeciesID > m_FileHdr.NumSpecies)
	return(false);
return(m_FiltSpecies[SpeciesID-1] & cMASeqOverlayFlg ? true : false);
}

// GetNumSpecies
// Get total number of species which are represented as being aligned
int 
CMAlignFile::GetNumSpecies(tAlignBlockID BlockID)			// returns number of species aligned
{
int Rslt;
if(!BlockID)
	return(m_FileHdr.NumSpecies);
if((Rslt=LoadBlock(BlockID))!=eBSFSuccess)
	return(Rslt);
return(m_pAlignBlock->NumSpecies);
}

// returns the next species identifier in specified block, 0 if no more to return
int 
CMAlignFile::GetNxtBlockSpeciesID(tAlignBlockID CurBlockID,tSpeciesID CurSpeciesID)	
{
int Rslt;
int NumSpecies;
int ChromID;
int NxtSpeciesID;
int LoSpeciesID;
tsAlignSpecies *pSpecies;
if(!m_FileHdr.NumAlignBlocks || CurBlockID < 1 || CurBlockID > m_FileHdr.NumAlignBlocks)
	return(eBSFerrAlignBlk);
if((Rslt=LoadBlock(CurBlockID))!=eBSFSuccess)
	return(Rslt);
NumSpecies = m_pAlignBlock->NumSpecies;
if(!NumSpecies)
	return(0);

pSpecies = (tsAlignSpecies *)((unsigned char *)m_pAlignBlock + sizeof(tsAlignBlock));
NxtSpeciesID = cMaxAlignedSpecies+1;		// used to indicate that NxtSpeciesID does not exist, if less than (cMaxAlignedSpecies+1) then it holds the next species identifier
while(ChromID = pSpecies->ChromID) 
	{
	LoSpeciesID = m_pChromNames[ChromID-1].SpeciesID;
	if(LoSpeciesID > CurSpeciesID && LoSpeciesID < NxtSpeciesID)
		NxtSpeciesID = LoSpeciesID;
	if(!--NumSpecies)
		break;
	pSpecies = (tsAlignSpecies *)((unsigned char *)pSpecies + pSpecies->AlignSpeciesLen);
	}
return(NxtSpeciesID > cMaxAlignedSpecies ? 0 : NxtSpeciesID);
}


// GetSpeciesName
// Get species name for specified species identifier
char *
CMAlignFile::GetSpeciesName(tSpeciesID SpeciesID)	// returns species name corresponding to species identifier (1..n)
{
if(!SpeciesID || SpeciesID > m_FileHdr.NumSpecies)
	return(NULL);
return((char *)m_FileHdr.SpeciesNames[SpeciesID-1].szSpeciesName);
}


// GetSpeciesID
// Returns species identifer for specified chromid
int
CMAlignFile::GetSpeciesID(tChromID ChromID)	
{
tsChromName *pChrom;
int Idx;
if(!ChromID || ChromID > m_FileHdr.NumChroms)
	return(-1);
pChrom = m_pChromNames;
for(Idx = 0; Idx < m_FileHdr.NumChroms; Idx++,pChrom++)
	{
	if(pChrom->ChromID == ChromID)
		return(pChrom->SpeciesID);
	}
return(eBSFerrChrom);
}



// NxtBlock
// Returns next block identifier relative to specified block identifier
int 
CMAlignFile::NxtBlock(tAlignBlockID CurBlockID)				// returns next alignment block identifier - 0 returns 1st block identifier
{
int SpeciesIdx;
tsAlignSpecies *pSpecies;

if(!m_FileHdr.NumAlignBlocks || CurBlockID >= m_FileHdr.NumAlignBlocks)
	return(eBSFerrAlignBlk);

if(!m_NumFiltSpecies)
	return(CurBlockID + 1);

// need to check if block needs to be skipped
while(CurBlockID < m_FileHdr.NumAlignBlocks)	
	{
	CurBlockID += 1;
	for(SpeciesIdx = 0; SpeciesIdx < m_FileHdr.NumSpecies; SpeciesIdx++)
		{
			// if this is a species to be filtered...
		if(m_FiltSpecies[SpeciesIdx] & cMASeqOverlayFlg)
			{
			if((pSpecies = LoadAlignSpecies(CurBlockID,SpeciesIdx+1))==NULL)
				continue; // can't filter on a species which is not aligned!

			if(pSpecies->ConfScore & cMASeqOverlayFlg)
				break;		// can't return this block
			}
		}
	if(SpeciesIdx == m_FileHdr.NumSpecies)
		return(CurBlockID);
	}

return(eBSFerrAlignBlk);
}

// NxtBlockChrom
// Returns next block identifier on specified reference chromosome
int 
CMAlignFile::NxtBlockChrom(tAlignBlockID CurBlockID, // current block identifier (0 returns 1st)
						   tChromID ChromID)  // reference chromosome
{
tsBlockDirEl *pDirEl;
int NxtBlockID;

if(!m_FileHdr.NumAlignBlocks || CurBlockID > m_FileHdr.NumAlignBlocks)
	return(eBSFerrAlignBlk);
if(ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);

// locate block with lowest ofs in chromosome
if((pDirEl = LocateDirEl(ChromID,0,true)) == NULL)
	return(eBSFerrChrom);

// iterate through chromosome until blockid following current is located
NxtBlockID = INT_MAX;
do {
	if(pDirEl->ChromID != ChromID)	// will break also on last array element
		break;
	if(pDirEl->BlockID <= CurBlockID)
		continue;
	if(pDirEl->BlockID < NxtBlockID)
		NxtBlockID = pDirEl->BlockID;
	}
while(pDirEl++); // safe not to check for end of array as last element is marked by ChromID == eBSFerrChrom
return(NxtBlockID == INT_MAX ? eBSFerrChrom : NxtBlockID);
}

int 
CMAlignFile::GetMaxAlignLen(void)					// returns maximal length alignment (includes '-' indel markers)
{
return(m_FileHdr.AlignIncInDelLen);
}

int 
CMAlignFile::GetNumBlocks(void)				// returns number of blocks in alignment
{
return(m_FileHdr.NumAlignBlocks);
}

int 
CMAlignFile::GetNumChroms(void)				// returns number of chromosomes in alignment
{
return(m_FileHdr.NumChroms);
}

int 
CMAlignFile::GetAlignLen(tAlignBlockID BlockID,tSpeciesID SpeciesID)	// get alignment length (includes '-' indel markers)
{
tsAlignSpecies *pAlignSpecies;
if((pAlignSpecies = LoadAlignSpecies(BlockID,SpeciesID))==NULL) // ensures species is represented plus loads block
	return(eBSFerrAlignBlk);
return(m_pAlignBlock->AlignIncInDelLen);
}



// GetRelChromOfs
// Gets genomic alignment start position for species in specified block
// 
int 
CMAlignFile::GetRelChromOfs(tAlignBlockID BlockID,tSpeciesID SpeciesID)	// get alignment start psn
{
tsAlignSpecies *pAlignSpecies;
if((pAlignSpecies = LoadAlignSpecies(BlockID,SpeciesID))==NULL)
	return(eBSFerrAlignBlk);
return(pAlignSpecies->ChromOfs);
}

// GetRelChromEndOfs
// Gets genomic alignment end position for species in specified block
int 
CMAlignFile::GetRelChromEndOfs(tAlignBlockID BlockID,tSpeciesID SpeciesID)		// get alignment end psn
{
tsAlignSpecies *pAlignSpecies;
if((pAlignSpecies = LoadAlignSpecies(BlockID,SpeciesID))==NULL)
	return(eBSFerrAlignBlk);
return(pAlignSpecies->AlignXInDelLen + pAlignSpecies->ChromOfs - 1);
}

// GetStrand
// Gets strand for species in specified block
char 
CMAlignFile::GetStrand(tAlignBlockID BlockID,tSpeciesID SpeciesID) // set strand - '+' or '-'
{
tsAlignSpecies *pAlignSpecies;
if((pAlignSpecies = LoadAlignSpecies(BlockID,SpeciesID))==NULL)
	return('+');	// default to '+' strand - no error returned!
return(pAlignSpecies->Strand);
}

// GetRelChrom
// Gets chromosome for species in specified block
char *
CMAlignFile::GetRelChrom(tAlignBlockID BlockID,tSpeciesID RelSpeciesID)  // get alignment chromosome
{
tsAlignSpecies *pAlignSpecies;
if((pAlignSpecies = LoadAlignSpecies(BlockID,RelSpeciesID))==NULL)
	return(NULL);
return((char *)m_pChromNames[pAlignSpecies->ChromID-1].szChromName);
}

// GetRelChromID
// get alignment chromosome for species
int 
CMAlignFile::GetRelChromID(tAlignBlockID BlockID,tSpeciesID RelSpeciesID)
{
tsAlignSpecies *pAlignSpecies;
if((pAlignSpecies = LoadAlignSpecies(BlockID,RelSpeciesID))==NULL)
	return(eBSFerrChrom);
return(pAlignSpecies->ChromID);
}

// GetChromName
// returns chromosome name
char *
CMAlignFile::GetChromName(tChromID ChromID)  
{
if(ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(NULL);
return((char *)m_pChromNames[ChromID-1].szChromName);
}

// GetChromLen
// returns chromosome length
int  
CMAlignFile::GetChromLen(tChromID ChromID)			  
{
if(ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
return(m_pChromNames[ChromID-1].ChromLen);
}

// MapChromOfsToOtherStrand
// Maps chromosome offset (0..ChromLen-1) on specified chromosome to other strand
// Usually used if offset is returned for '-' strand and this offset needs to be mapped to '+' strand
int  
CMAlignFile::MapChromOfsToOtherStrand(tChromID ChromID,int ChromOfs)			  
{
int ChromLen;
if(ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
ChromLen = m_pChromNames[ChromID-1].ChromLen;
if(ChromOfs < 0 || ChromOfs >= ChromLen)
	return(eBSFerrParams);
return((ChromLen - 1) - ChromOfs);
}


int												// returns number of bases (incl InDels and NotAlignedBase)
CMAlignFile::LoadAlignment(char *pszRefChrom,	// alignment is to this ref species chromosome
					  int StartLoci,			// starting from this ref chromosome loci
					  int EndLoci,				// and ending at this ref chromosome loci
					  char *pszRelSpecies,		// aligned bases are to be returned for this species
					  int MaxBases,				// max number of bases to be returned
					  etSeqBase *pRefBases,		// where to return ref alignment bases
					  etSeqBase *pRelBases,		// where to return rel alignment bases
					  etSeqBase RefUnalignedBase, // base to return if no ref species alignment
					  etSeqBase RelUnalignedBase) // base to return if ref species aligned but no alignment onto rel species
{
int Rslt;
INT32 RefChromID;
tSpeciesID RelSpeciesID;
if((Rslt = RefChromID = LocateChromID(GetRefSpeciesID(),pszRefChrom))<0)
	return(Rslt);
if((Rslt = RelSpeciesID = LocateSpeciesID(pszRelSpecies))<0)
	return(Rslt);
return(LoadAlignment(RefChromID,StartLoci,EndLoci,RelSpeciesID,MaxBases,pRefBases,pRelBases,RefUnalignedBase,RelUnalignedBase));
}

int												// returns number of bases (incl InDels and NotAlignedBase)
CMAlignFile::LoadAlignment(INT32 RefChromID,	// alignment is to this ref species chromosome
					  int StartLoci,			// starting from this ref chromosome loci
					  int EndLoci,				// and ending at this ref chromosome loci
					  int RelSpeciesID,			// aligned bases are to be returned for this species
					  int MaxBases,				// max number of bases to be returned
					  etSeqBase *pRefBases,		// where to return ref alignment bases
					  etSeqBase *pRelBases,		// where to return rel alignment bases
					  etSeqBase RefUnalignedBase, // base to return if no ref species alignment
					  etSeqBase RelUnalignedBase) // base to return if ref species aligned but no alignment onto rel species
{
int RefSpeciesID;
int BasesReturned;
tsBlockDirEl *pDirEl;
etSeqBase *pSeq;
int	SeqOfs;
int	SeqLoci;
int	SeqLen;
int AlignLen;
etSeqBase Base;

// check parameters
if(RefChromID < 0 || RefChromID > m_FileHdr.NumChroms ||
   MaxBases < 1 || StartLoci < 0 || StartLoci > EndLoci ||
   RelSpeciesID < 1 || RelSpeciesID > m_FileHdr.NumSpecies ||
   pRefBases == NULL || pRelBases == NULL)
	return(eBSFerrParams);

RefSpeciesID = GetRefSpeciesID();
BasesReturned = 0;

do {
	if((pDirEl = LocateDirEl(RefChromID,StartLoci,true))==NULL)
		{
		int NumUndef = EndLoci - StartLoci + 1;
		if(NumUndef > MaxBases - BasesReturned)
			NumUndef = MaxBases - BasesReturned;

		memset(pRefBases,RefUnalignedBase,NumUndef);
		memset(pRelBases,RefUnalignedBase,NumUndef);
		BasesReturned += NumUndef;
		break;
		}
	if(StartLoci < pDirEl->ChromOfs)		// if block starts after requested start loci then no alignment yet...
		{
		while(StartLoci <= EndLoci && StartLoci < pDirEl->ChromOfs && BasesReturned < MaxBases)
			{
			*pRefBases++ = RefUnalignedBase;			// no alignment
			*pRelBases++ = RefUnalignedBase;
			StartLoci++;
			BasesReturned += 1;
			}
		}

	if(BasesReturned == MaxBases || pDirEl->ChromOfs > EndLoci)
		return(BasesReturned);

	// load block sequence for reference
	AlignLen = GetAlignLen(pDirEl->BlockID,RefSpeciesID);
	pSeq = GetSeq(pDirEl->BlockID,RefSpeciesID);
	SeqOfs = 0;
	SeqLoci = pDirEl->ChromOfs;
	SeqLen = 0;
	while(AlignLen-- && StartLoci <= EndLoci && BasesReturned < MaxBases)
		{
		Base = *pSeq++;
		if(!SeqLen && SeqLoci < StartLoci)
			{
			if((Base  & ~cRptMskFlg) <= eBaseN)
				SeqLoci++;
			SeqOfs++;
			continue;
			}
		*pRefBases++ = Base;
		SeqLen += 1;
		BasesReturned += 1;
		if((Base  & ~cRptMskFlg) <= eBaseN)
			StartLoci++;
		}

	// load block sequence for relative species
	if((pSeq = GetSeq(pDirEl->BlockID,RelSpeciesID))!=NULL)
		memmove(pRelBases,&pSeq[SeqOfs],SeqLen);
	else
		memset(pRelBases,RelUnalignedBase,SeqLen);
	pRelBases += SeqLen;
	}
while(BasesReturned < MaxBases && StartLoci <= EndLoci);

return(BasesReturned);
}

// LoadAlignments
// Loads multiple species alignments
int												// returns number of bases (incl InDels and NotAlignedBase)
CMAlignFile::LoadAlignments(INT32 RefChromID,	// alignment is to this ref species chromosome
					  int StartLoci,			// starting from this ref chromosome loci
					  int EndLoci,				// and ending at this ref chromosome loci
					  int MaxBases,				// max number of bases for each species to be returned
					  int NumSpecies,			// number of species to return alignments for
					  int *pSpeciesIDs,			// pts to array of species identifiers (ref species must be first)
					  int RelOfs,	   		    // offset in pSeqBases[species][] at which to return alignment bases
					  etSeqBase *pSeqBases[],	// pts to array of ptrs of where to return alignment bases for each species			
					  etSeqBase RefUnalignedBase, // base to return if no ref species alignment
					  etSeqBase RelUnalignedBase) // base to return if ref species aligned but no alignment onto rel species
{
int RefSpeciesID;
int BasesReturned;
tsBlockDirEl *pDirEl;
etSeqBase *pLocSeqBases[cMaxAlignedSpecies];
etSeqBase *pSeq;
int	SeqOfs;
int	SeqLoci;
int	SeqLen;
int AlignLen;
int CurSpecies;
etSeqBase Base;

// check parameters
if(RefChromID < 0 || RefChromID > m_FileHdr.NumChroms ||
   NumSpecies < 1 || NumSpecies > m_FileHdr.NumSpecies ||
   MaxBases < 1 || StartLoci < 0 || StartLoci > EndLoci ||
   pSpeciesIDs == NULL || pSeqBases == NULL)
	return(eBSFerrParams);
for(CurSpecies = 0; CurSpecies < NumSpecies; CurSpecies++)
	{
	if(pSpeciesIDs[CurSpecies] < 1 || pSpeciesIDs[CurSpecies] > m_FileHdr.NumSpecies)
		return(eBSFerrSpecies);
	if((pLocSeqBases[CurSpecies] = pSeqBases[CurSpecies]) == NULL)
		return(eBSFerrParams);
	pLocSeqBases[CurSpecies] += RelOfs;
	 }

RefSpeciesID = GetRefSpeciesID();
if(RefSpeciesID != pSpeciesIDs[0])					// double check!
	return(eBSFerrSpecies);

BasesReturned = 0;
do {
	if((pDirEl = LocateDirEl(RefChromID,StartLoci,true))==NULL)
		{
		int NumUndef = EndLoci - StartLoci + 1;
		if(NumUndef > MaxBases - BasesReturned)
			NumUndef = MaxBases - BasesReturned;
		for(CurSpecies = 0; CurSpecies < NumSpecies; CurSpecies++)
			memset(pLocSeqBases[CurSpecies],RefUnalignedBase,NumUndef);
		BasesReturned += NumUndef;
		break;
		}

	if(StartLoci < pDirEl->ChromOfs)		// if block starts after requested start loci then no alignment yet...
		{
		while(StartLoci <= EndLoci && StartLoci < pDirEl->ChromOfs && BasesReturned < MaxBases)
			{
			for(CurSpecies = 0; CurSpecies < NumSpecies; CurSpecies++)
				*pLocSeqBases[CurSpecies]++ = RefUnalignedBase;
			StartLoci++;
			BasesReturned += 1;
			}
		}

	if(BasesReturned == MaxBases || pDirEl->ChromOfs > EndLoci)
		return(BasesReturned);

	// load block sequence for reference
	AlignLen = GetAlignLen(pDirEl->BlockID,RefSpeciesID);
	pSeq = GetSeq(pDirEl->BlockID,RefSpeciesID);
	SeqOfs = 0;
	SeqLoci = pDirEl->ChromOfs;
	SeqLen = 0;
	while(AlignLen-- && StartLoci <= EndLoci && BasesReturned < MaxBases)
		{
		Base = *pSeq++;
		if(!SeqLen && SeqLoci < StartLoci)
			{
			if((Base  & ~cRptMskFlg) <= eBaseN)
				SeqLoci++;
			SeqOfs++;
			continue;
			}
		*pLocSeqBases[0]++ = Base;
		SeqLen += 1;
		BasesReturned += 1;
		if((Base  & ~cRptMskFlg) <= eBaseN)
			StartLoci++;
		}

	// now load block sequences for each of the relative species
	for(CurSpecies = 1; CurSpecies < NumSpecies; CurSpecies++)
		{
		if((pSeq = GetSeq(pDirEl->BlockID,pSpeciesIDs[CurSpecies]))!=NULL)
			memmove(pLocSeqBases[CurSpecies],&pSeq[SeqOfs],SeqLen);
		else
			memset(pLocSeqBases[CurSpecies],RelUnalignedBase,SeqLen);
		pLocSeqBases[CurSpecies] += SeqLen;
		}
	}
while(BasesReturned < MaxBases && StartLoci <= EndLoci);

return(BasesReturned);
}



// GetSeq
// Gets sequence for species in specified block 
etSeqBase *
CMAlignFile::GetSeq(tAlignBlockID BlockID,tSpeciesID SpeciesID)   // returns alignment sequence for specified block and species
{
tsAlignSpecies *pAlignSpecies;
if((pAlignSpecies = LoadAlignSpecies(BlockID,SpeciesID))==NULL)
	return(NULL);
return((etSeqBase *)((char *)pAlignSpecies + sizeof(tsAlignSpecies)));
}


// GetBlockID
// returns BlockID for block in which RefChromID has a datapoint at ChromOfs
int 
CMAlignFile::GetBlockID(tChromID RefChromID,int ChromOfs)
{
tsBlockDirEl *pDirEl;
pDirEl = LocateDirEl(RefChromID,ChromOfs);
return(pDirEl == NULL ? eBSFerrChrom : pDirEl->BlockID);
}



// MultiAlignAll2DataPts
// Process bioseq multialignment file into a bioseq data points file
// Major difference between these file types is that the multialignment file retains the
// concept of blocks of alignments from the orginal source MAF files where as the data points file
// contains datasegments which are species sequences derived from blocks and hence which may extend over
// multiple blocks. These datasegments are normalised to the reference sequence such that the reference 
// alignment has indels '-' removed so it is contiguous, and the relative sequences realigned to the
// contiguous reference. 
int 
CMAlignFile::MultiAlignAll2DataPts(char *pszMAF,// source bioseq mutialignment file
				char *pszDataPointsFile,		// DataPoints file
				char *pszDescr,					// file description
				char *pszTitle)					// file title
{
int RefSpeciesID,RelSpeciesID,RefChromID, RelChromID;
char *pszRelSpecies,*pszRefSpecies;
char *pszRefChrom,*pszRelChrom;
int RefChromOfs,RelChromOfs;

etSeqBase *pRefSeq;
etSeqBase *pRelSeq;
etSeqBase *pRefBase;
etSeqBase *pRelBase;
etSeqBase *pRef2Struct;
etSeqBase *pRel2Struct;
etSeqBase *pSeq2proc;
etSeqBase RefBase;
etSeqBase RelBase;
char RefStrand;

int NumSpecies;									// number of aligned species
int Rslt;

UINT32 CurBlockID;
int MaxAlignLen;
UINT32 AlignPsn;
UINT32 RefAlignStart;
UINT32 RelAlignStart;
UINT32 AlignLen;
UINT32 RefAlignLen;
UINT32 RefStructLen;
UINT32 RelStructLen;
char RelStrand;

int PrevRefChromID;
int RefChromLen;
int RelDatasetID;
int RelChromLen;


tsBlockDirEl *pDirEl;
int DirElIdx;


UINT32 AllocLen = 0;
unsigned char *pSeqBuff = NULL;

if((Rslt=Open(pszMAF))!=eBSFSuccess)
	{
	AddErrMsg("CMAlignFile::MultiAlignAll2DataPts","Unable to open bioalign MAF file %s\n",pszMAF);
	return(Rslt);
	}

// get the reference species
pszRefSpecies = GetRefSpeciesName();

// ensure ref species is represented in alignment file
if((Rslt = RefSpeciesID = LocateSpeciesID(pszRefSpecies))<= eBSFSuccess)
	{
	AddErrMsg("CMAlignFile::MultiAlignAll2DataPts","Unable to locate ref species %s in multialignment file %s\n",pszRefSpecies,pszMAF);
	return(Rslt);
	}

CDataPoints *pDataPoints = new CDataPoints;
if(pDataPoints->Open(pszDataPointsFile,true) != eBSFSuccess)
	{
	delete pDataPoints;
	return(false);
	}
pDataPoints->SetDescription(pszDescr);
pDataPoints->SetTitle(pszTitle);

// grab enough memory to hold maximal length aligned sequences + guard
pRef2Struct = pRel2Struct = NULL;
MaxAlignLen = GetMaxAlignLen();
pRef2Struct = new unsigned char [MaxAlignLen + 1];
pRel2Struct = new unsigned char [MaxAlignLen + 1];
if(pRef2Struct == NULL || pRel2Struct == NULL)
	{
	AddErrMsg("CMAlignFile::MultiAlignAll2DataPts","Unable to allocate memory to hold aligned sequences");
	pDataPoints->Close(false);
	delete pDataPoints;
	pDataPoints = NULL;
	
	if(pRef2Struct != NULL)
		delete pRef2Struct;
	if(pRel2Struct != NULL)
		delete pRel2Struct;
	return(false);
	}

// datapoints are always going to be packed bases - 2 per byte
if((Rslt=pDataPoints->SetRefDataset(pszRefSpecies,eDPTPackedBases)) <= eBSFSuccess)
	{
	AddErrMsg("CMAlignFile::MultiAlignAll2DataPts","Unable to set reference sequence %s into pDataPoints",pszRefSpecies);
	pDataPoints->Close(false);
	delete pDataPoints;
	pDataPoints = NULL;
	
	if(pRef2Struct != NULL)
		delete pRef2Struct;
	if(pRel2Struct != NULL)
		delete pRel2Struct;
	return(false);
	}

// determine the number of species so we know the range of species identifiers ( 1..NumSpecies)
NumSpecies = GetNumSpecies();

// iterate over reference block directory array which contains elements sorted by chrom then offset
PrevRefChromID = -1;
pDirEl = m_pDirEls;
for(DirElIdx = 0; DirElIdx < (int)m_FileHdr.NumAlignBlocks; DirElIdx++,pDirEl++)
	{
	RefChromID = pDirEl->ChromID;
	CurBlockID = pDirEl->BlockID;
	pszRefChrom = GetChromName(RefChromID);
	
	LoadBlock(pDirEl);		// pull block into memory from disk 

	// not interested in alignments in which there is only a single sequence - Yep, maf's sometimes do contain single sequences!
	if(GetNumSpecies(CurBlockID)==1)
		continue;

	RefAlignStart= GetRelChromOfs(CurBlockID,RefSpeciesID);
	RefAlignLen =  GetAlignLen(CurBlockID,RefSpeciesID);
	RefStrand = GetStrand(CurBlockID,RefSpeciesID);

		// get the reference sequence
	if((pRefSeq = GetSeq(CurBlockID,RefSpeciesID))==NULL)
		continue;

	pRefBase = pRefSeq;
	AlignPsn = RefAlignStart;
	AlignLen = RefAlignLen;
	pSeq2proc = pRef2Struct;
	RefStructLen = 0;
	while(AlignLen--)
		{
		RefBase = *pRefBase++ & ~cRptMskFlg;
		if(RefBase != eBaseInDel)
			{
			*pSeq2proc++ = RefBase;
			if(!RefStructLen)
				RefChromOfs = AlignPsn;
			RefStructLen++;
			AlignPsn++;
			}
		}

// debug!!
static int PrvOfs = -1;
static int PrvLen = -1;
static int PrvChromID = -1;
static int PrvBlockID = -1;
if(PrvBlockID != -1 && PrvChromID == RefChromID)
	{
	if(PrvOfs > RefChromOfs)
		printf("\nPrevBlock:%d %s:%d is after CurBlock:%d Ofs:%d",PrvBlockID,pszRefChrom,PrvOfs,pDirEl->BlockID,RefChromOfs);
	else
		if((PrvOfs + PrvLen) > RefChromOfs)
			printf("\nPrevBlock:%d %s:%d len:%d overlaps CurBlock:%d Ofs:%d",PrvBlockID,pszRefChrom,PrvOfs,PrvLen,pDirEl->BlockID,RefChromOfs);
	}
else
	if(PrvChromID != RefChromID)
		printf("\nProcessing reference chromosome %s",pszRefChrom);

PrvOfs = RefChromOfs;
PrvLen = RefStructLen;
PrvChromID = RefChromID;
PrvBlockID = pDirEl->BlockID;
// end debug


	if(PrevRefChromID != RefChromID)
		{
		PrevRefChromID = RefChromID;
		RefChromLen = GetChromLen(RefChromID);
		pDataPoints->AddRefChrom(pszRefChrom,RefChromLen);
		}
	Rslt=JoinSegments(pDataPoints,
				    pszRefChrom,			// reference species chromosome
					RefChromOfs,			// reference species chromosome offset
					pszRefSpecies,			// into which species relative dataset
					pszRefChrom,			// on which chromosome
					RefChromOfs,			// relative species chromosome offset
					RefStrand,				// relative species strand
					RefStructLen,			// number of datapoints
					pRef2Struct);			// array of data points
	if(Rslt != eBSFSuccess)
			{
			if(pRef2Struct != NULL)
				delete pRef2Struct;
			if(pRel2Struct != NULL)
				delete pRel2Struct;
			pDataPoints->Close(false);
			delete pDataPoints;
			return(Rslt);
			}

	// iterate over relative species aligned in current block
	for(RelSpeciesID = 1;RelSpeciesID <= NumSpecies;RelSpeciesID++) 
		{
		if(RelSpeciesID == RefSpeciesID)
			continue;

		// get each relative sequence
		if((pRelSeq = GetSeq(CurBlockID,RelSpeciesID))==NULL) // some species not present in all blocks
			continue;

		pszRelSpecies = GetSpeciesName(RelSpeciesID);
		RelChromID = GetRelChromID(CurBlockID,RelSpeciesID);	
		pszRelChrom = GetChromName(RelChromID);

		RelChromLen = GetChromLen(RelChromID);
		RelDatasetID = pDataPoints->AddRelDataset(pszRelSpecies);
		Rslt = pDataPoints->AddRelChrom(pszRelChrom,RelDatasetID,RelChromLen);
		if(Rslt < eBSFSuccess)
			{
			AddErrMsg("CMAlignFile::MultiAlignAll2DataPts","Unable add relative chromosome %s",pszRelChrom);
			return(Rslt);
			}
		pRelBase = pRelSeq;
		RelAlignStart= GetRelChromOfs(CurBlockID,RelSpeciesID);

		RelChromOfs = GetRelChromOfs(CurBlockID,RelSpeciesID);
		RelStrand = GetStrand(CurBlockID,RelSpeciesID);

		pRefBase = pRefSeq;
		AlignPsn = RefAlignStart;
		AlignLen = RefAlignLen;

		pSeq2proc = pRef2Struct;
		RelStructLen = 0;
		while(AlignLen--)
			{
			RelBase = *pRelBase++ & ~cRptMskFlg;
			RefBase = *pRefBase++ & ~cRptMskFlg;
			if(RefBase != eBaseInDel)
				{
				*pSeq2proc++ = RelBase;
				if(!RelStructLen)
					RefChromOfs = AlignPsn;
				RelStructLen++;
				}
			}
		Rslt=JoinSegments(pDataPoints,
				    pszRefChrom,			// reference species chromosome
					RefChromOfs,			// reference species chromosome offset
					pszRelSpecies,			// into which species relative dataset
					pszRelChrom,			// on which chromosome
					RelChromOfs,			// relative species chromosome offset
					RelStrand,				// relative strand
					RefStructLen,			// number of datapoints
					pRef2Struct);			// array of data points
		if(Rslt != eBSFSuccess)
			{
			if(pRef2Struct != NULL)
				delete pRef2Struct;
			if(pRel2Struct != NULL)
				delete pRel2Struct;
			pDataPoints->Close(false);
			delete pDataPoints;
			return(Rslt);
			}
		}
	}

// commit segments to disk
FlushSegs2Dataset(pDataPoints);

if(pRef2Struct != NULL)
	delete pRef2Struct;
if(pRel2Struct != NULL)
	delete pRel2Struct;
if(pDataPoints != NULL)
	{
	pDataPoints->Close(true);
	delete pDataPoints;
	}

return(eBSFSuccess);
}


int
CMAlignFile::NewSegEntry(tsSegCache *pEntry,
					unsigned short RefChromHash,
					unsigned short RelDatasetHash,
				     char *pszRefChrom,			// reference species chromosome
					 int RefChromOfs,			// reference species chromosome offset
					 char *pszRelDataset,		// into which species relative dataset
					 char *pszRelChrom,			// on which chromosome
					 int RelChromOfs,			// relative species chromosome offset
					 char RelStrand,				// relative strand '+' or '-'
					 int NumPts,				// number of datapoints
					 unsigned char *pDataPts)	// array of data points
{
pEntry->RelDatasetHash = RelDatasetHash;
strncpy(pEntry->szRelDataset,pszRelDataset,cMaxDatasetSpeciesChrom);
pEntry->szRelDataset[cMaxDatasetSpeciesChrom-1] = '\0';
pEntry->RefChromHash = RefChromHash;
strncpy(pEntry->szRefChrom,pszRefChrom,cMaxDatasetSpeciesChrom);
pEntry->szRefChrom[cMaxDatasetSpeciesChrom-1] = '\0';
pEntry->RefChromOfs = RefChromOfs;
strncpy(pEntry->szRelChrom,pszRelChrom,cMaxDatasetSpeciesChrom);
pEntry->szRelChrom[cMaxDatasetSpeciesChrom-1] = '\0';
pEntry->RelChromOfs = RefChromOfs;
pEntry->RelStrand = RelStrand;

if(pEntry->pDataPts == NULL ||pEntry->AllocdSize < NumPts)
	{
	int ReqSize;
	if(pEntry->pDataPts != NULL)
		delete pEntry->pDataPts;
    pEntry->AllocdSize = 0;
	ReqSize = 10000;
	if(NumPts > 10000)
		ReqSize = (NumPts * 4) / 3;
	if((pEntry->pDataPts = new unsigned char [ReqSize])==NULL)
		return(eBSFerrMem);
	pEntry->AllocdSize = ReqSize;
	}
memmove(pEntry->pDataPts,pDataPts,NumPts);
pEntry->NumPts = NumPts;
return(eBSFSuccess);
}

int
CMAlignFile::FlushSegs2Dataset(CDataPoints *pDataPoints)
{
tsSegCache *pEntry;
int Cnt;
int Rslt;
if(!m_NumCachedSegs)
	return(eBSFSuccess);
pEntry = m_pSegCache;
for(Cnt = 0; Cnt < m_NumCachedSegs; Cnt++,pEntry++)
	{
	if((Rslt=pDataPoints->AddDataseg(pEntry->szRefChrom,pEntry->RefChromOfs,
							pEntry->szRelDataset,pEntry->szRelChrom,pEntry->RelChromOfs,pEntry->RelStrand,
						pEntry->NumPts,pEntry->pDataPts))!=eBSFSuccess)
		return(Rslt);
	pEntry->RefChromHash = 0;
	pEntry->RefChromOfs = 0;
	pEntry->RelStrand = '+';
	pEntry->RelChromOfs = 0;
	pEntry->RelDatasetHash = 0;
	pEntry->NumPts = 0;
	pEntry->szRefChrom[0] = '\0';
	pEntry->szRelChrom[0] = '\0';
	pEntry->szRelDataset[0] = '\0';
	}
m_NumCachedSegs = 0;
return(eBSFSuccess);
}

int
CMAlignFile::JoinSegments(CDataPoints *pDataPoints,
				     char *pszRefChrom,			// reference species chromosome
					 int RefChromOfs,			// reference species chromosome offset
					 char *pszRelDataset,		// into which species relative dataset
					 char *pszRelChrom,			// on which chromosome
					 int RelChromOfs,			// relative species chromosome offset
					 char RelStrand,			// relative strand '+' or '-'
					 int NumPts,				// number of datapoints
					 unsigned char *pDataPts)	// array of data points
{
unsigned short RelDatasetHash = GenNameHash(pszRelDataset);
unsigned short RefChromHash   = GenNameHash(pszRefChrom);
tsSegCache *pEntry;
int Cnt;
int Rslt;
if(m_pSegCache == NULL)
	{
	if((m_pSegCache = new tsSegCache[cDSFMaxDatasets])==NULL)
		return(eBSFerrMem);
	memset(m_pSegCache,0,sizeof(tsSegCache) * cDSFMaxDatasets);
	m_NumCachedSegs = 0;
	}

pEntry = m_pSegCache;	
if(m_NumCachedSegs > 0)
	{
	for(Cnt = 0; Cnt < m_NumCachedSegs; Cnt++, pEntry++)
		{
		if(pEntry->RelDatasetHash == RelDatasetHash &&
			!stricmp(pEntry->szRelDataset,pszRelDataset))
			break;
		}

	if(Cnt > cDSFMaxDatasets)
		return(eBSFerrMaxDatasets);
	}
else
	Cnt = 0;

if(!m_NumCachedSegs || Cnt >= m_NumCachedSegs)	// new dataset starting...
	{
	m_NumCachedSegs += 1;
	return(NewSegEntry(pEntry,RefChromHash,RelDatasetHash,pszRefChrom,RefChromOfs,
								pszRelDataset,pszRelChrom,RelChromOfs,RelStrand,NumPts,pDataPts));
	}

// located existing dataset, see if can join new segment onto existing
if(pEntry->RefChromHash != RefChromHash ||		// make sure it's the same chromosome!!
	stricmp(pEntry->szRefChrom,pszRefChrom) ||
	pEntry->RelStrand != RelStrand ||
	pEntry->RefChromOfs + pEntry->NumPts != RefChromOfs) // and segments are contiguous
	{
	// different chromosome or segments not contiguous...
	// write out existing and create new
	
	if((Rslt=pDataPoints->AddDataseg(pEntry->szRefChrom,pEntry->RefChromOfs,
					pEntry->szRelDataset,pEntry->szRelChrom,pEntry->RelChromOfs,pEntry->RelStrand,
						pEntry->NumPts,pEntry->pDataPts))!=eBSFSuccess)
		return(Rslt);
	return(NewSegEntry(pEntry,RefChromHash,RelDatasetHash,pszRefChrom,RefChromOfs,
								pszRelDataset,pszRelChrom,RelChromOfs,RelStrand,NumPts,pDataPts));
	}

// same chromosome, no gap so extend...
if(pEntry->AllocdSize < pEntry->NumPts + NumPts) // need to alloc more memory?
	{
	unsigned char *pNewAlloc;
	int ReqSize = 10000;
	if(pEntry->NumPts + NumPts > 10000)
		ReqSize = ((pEntry->NumPts + NumPts) * 4) / 3;
	if((pNewAlloc = new unsigned char [ReqSize])== NULL)
		return(eBSFerrMem);
	memmove(pNewAlloc,pEntry->pDataPts,pEntry->NumPts);
	delete pEntry->pDataPts;
	pEntry->pDataPts = pNewAlloc;
	pEntry->AllocdSize = ReqSize;
	}
memmove(&pEntry->pDataPts[pEntry->NumPts],pDataPts,NumPts);
pEntry->NumPts += NumPts;

return(eBSFSuccess);
}


// Locates and returns a ptr to matching tsBlockDirEl where the RefChromID exactly match, and
// the specified RefChromOfs is contained within the matching block
// If bClosest is true then if no exact match on ChromID.ChromOfs will reurn
// a ptr to the DirEl which has a ChromOfs closest but after the ChromOfs specified.
tsBlockDirEl *
CMAlignFile::LocateDirEl(tChromID ChromID,int ChromOfs,bool bClosest)
{
tsBlockDirEl *pProbe;
tsBlockDirEl *pClosest = NULL;
int Left = 0;
int Right = m_FileHdr.NumAlignBlocks - 1;
int MidPt;

while(Right >= Left) {
	MidPt = (Right + Left)/2;
	pProbe = &m_pDirEls[MidPt];

	if(pProbe->ChromID > ChromID)
		{
		Right = MidPt - 1;
		continue;
		}
	else
		if(pProbe->ChromID < ChromID)
			{
			Left = MidPt + 1;
			continue;
			}

	if(pProbe->ChromOfs > ChromOfs)
		{
		pClosest = pProbe;
		Right = MidPt - 1;
		continue;
		}
	else
		if((pProbe->ChromOfs + pProbe->AlignXInDelLen - 1) < ChromOfs)
			{
			Left = MidPt + 1;
			continue;
			}

	return(pProbe);
	}
return(bClosest ? pClosest : NULL);
}

int 
CMAlignFile::CompareBlockDirEls( const void *arg1, const void *arg2)
{
tsBlockDirEl *pEl1 = (tsBlockDirEl *)arg1;
tsBlockDirEl *pEl2 = (tsBlockDirEl *)arg2;
if(pEl1->ChromID == pEl2->ChromID)
	{
	if(pEl1->ChromOfs < pEl2->ChromOfs)
		return(-1);
	if(pEl1->ChromOfs > pEl2->ChromOfs)
		return(1);
	if(pEl1->AlignXInDelLen > pEl2->AlignXInDelLen) 
		return(-1);									
	if(pEl1->AlignXInDelLen < pEl2->AlignXInDelLen)
		return(1);
	return(0);
	}
return(pEl1->ChromID < pEl2->ChromID ? -1 : 1);
}

