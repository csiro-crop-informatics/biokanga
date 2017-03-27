/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

CBioSeqFile::CBioSeqFile(void)
{
m_hFile = -1;
m_pDataBuff = NULL;
m_pCreatingEntry = NULL;
m_pDirEntries = NULL;
m_ppDirEntryNames = NULL;
m_ppDirEntryIDs = NULL;
Reset(false);
}

CBioSeqFile::~CBioSeqFile(void)
{
if(m_hFile != -1)
	close(m_hFile);
if(m_pDataBuff != NULL)
	delete m_pDataBuff;
if(m_pDirEntries != NULL)
	delete m_pDirEntries;
if(m_ppDirEntryNames != NULL)
	delete m_ppDirEntryNames;
if(m_ppDirEntryIDs != NULL)
	delete m_ppDirEntryIDs;
}

// Hdr2Disk
// Writes file header out to disk
teBSFrsltCodes 
CBioSeqFile::Hdr2Disk(void)		
{
tsBSFHeader FileHdr;
tsBSFHeader *pHdr;
int WrtLen;

WrtLen = sizeof(tsBSFHeader);

if(m_bIsBigEndian)	// if on a big-endian machine then need to make little endian as that is our native file format
	{
	memmove(&FileHdr,&m_FileHdr,WrtLen);
	FileHdr.FileLen = SwapUI64Endians(m_FileHdr.FileLen);		// current file length
	FileHdr.SeqOfs = SwapUI64Endians(m_FileHdr.SeqOfs);			// where concatenated sequences start on disk
	FileHdr.SeqSize = SwapUI64Endians(m_FileHdr.SeqSize);		// disk/memory space required to hold concatenated sequences
	FileHdr.DirEntriesOfs = SwapUI64Endians(m_FileHdr.DirEntriesOfs);		// where directory entries start on disk
	FileHdr.DirEntryIDIdxOfs = SwapUI64Endians(m_FileHdr.DirEntryIDIdxOfs);	// where index by EntryID into directory entries start on disk
	FileHdr.DirEntryNameIdxOfs = SwapUI64Endians(m_FileHdr.DirEntryNameIdxOfs);	// where index by name into directory entries start on disk
	FileHdr.Type = SwapUI32Endians(m_FileHdr.Type);								// biosequence file type
	FileHdr.Version = SwapUI32Endians(m_FileHdr.Version);						// file structure version
	FileHdr.MaxEntries = SwapUI32Endians(m_FileHdr.MaxEntries);					// max number of entries supported 
	FileHdr.NumEntries = SwapUI32Endians(m_FileHdr.NumEntries);					// actual number of directory entries
	FileHdr.DirEntriesSize = SwapUI32Endians(m_FileHdr.DirEntriesSize);			// actual disk/memory space required to hold directory entries	}
	pHdr = &FileHdr;
	}
else
	pHdr = &m_FileHdr;
if(_lseeki64(m_hFile,0,SEEK_SET) ||
		write(m_hFile,pHdr,WrtLen)!=WrtLen)
	{
	AddErrMsg("CBioSeqFile::Hdr2Disk","Unable to write file header to disk on file %s - error %s",m_szFile,strerror(errno));
	Reset(false);
	return(eBSFerrFileAccess);
	}
m_bHdrDirty = false;
return(eBSFSuccess);
}

teBSFrsltCodes
CBioSeqFile::Disk2Hdr(char *pszSeqFile,int FileType)
{
if(_lseeki64(m_hFile,0,SEEK_SET)!=0)			// read in header..
	{
	AddErrMsg("CBioSeqFile::Disk2Hdr","Seek failed to offset 0 on %s - %s",pszSeqFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

if(sizeof(tsBSFHeader) != read(m_hFile,&m_FileHdr,sizeof(tsBSFHeader)))
	{
	AddErrMsg("CBioSeqFile::Disk2Hdr","Read of file header failed on %s - %s",pszSeqFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

// header read, validate it as being a sequence file header
if(tolower(m_FileHdr.Magic[0]) != 'b' || 
	tolower(m_FileHdr.Magic[1]) != 'i' || 
	tolower(m_FileHdr.Magic[2]) != 'o' || 
	tolower(m_FileHdr.Magic[3]) != 's')
	{
	AddErrMsg("CBioSeqFile::Disk2Hdr","%s opened but no magic signature - not a bioseq file",pszSeqFile);
	Reset(false);			// closes opened file..
	return(eBSFerrNotBioseq);
	}

if(m_bIsBigEndian)	// file was written with little-endian ordering
	{
	m_FileHdr.FileLen = SwapUI64Endians(m_FileHdr.FileLen);		// current file length
	m_FileHdr.SeqOfs = SwapUI64Endians(m_FileHdr.SeqOfs);			// where concatenated sequences start on disk
	m_FileHdr.SeqSize = SwapUI64Endians(m_FileHdr.SeqSize);		// disk/memory space required to hold concatenated sequences
	m_FileHdr.DirEntriesOfs = SwapUI64Endians(m_FileHdr.DirEntriesOfs);		// where directory entries start on disk
	m_FileHdr.DirEntryIDIdxOfs = SwapUI64Endians(m_FileHdr.DirEntryIDIdxOfs);	// where index by EntryID into directory entries start on disk
	m_FileHdr.DirEntryNameIdxOfs = SwapUI64Endians(m_FileHdr.DirEntryNameIdxOfs);	// where index by name into directory entries start on disk
	m_FileHdr.Type = SwapUI32Endians(m_FileHdr.Type);								// biosequence file type
	m_FileHdr.Version = SwapUI32Endians(m_FileHdr.Version);						// file structure version
	m_FileHdr.MaxEntries = SwapUI32Endians(m_FileHdr.MaxEntries);					// max number of entries supported 
	m_FileHdr.NumEntries = SwapUI32Endians(m_FileHdr.NumEntries);					// actual number of directory entries
	m_FileHdr.DirEntriesSize = SwapUI32Endians(m_FileHdr.DirEntriesSize);			// actual disk/memory space required to hold directory entries	}
	}

	// check bioseq file is the type we are expecting
if(FileType != cBSFTypeAny && m_FileHdr.Type != FileType)
	{
	AddErrMsg("CBioSeqFile::Disk2Hdr","%s opened as a bioseq file - expected type %d, file type is %d",pszSeqFile,FileType,m_FileHdr.Type);
	Reset(false);			// closes opened file..
	return(eBSFerrFileType);
	}

	// can we handle this version?
if(m_FileHdr.Version < cBSFVersionBack || m_FileHdr.Version > cBSFVersion)
	{
	AddErrMsg("CBioSeqFile::Disk2Hdr","%s opened as a bioseq file - expected between version %d and %d, file version is %d",pszSeqFile,
			cBSFVersionBack,cBSFVersion,m_FileHdr.Version);
	Reset(false);			// closes opened file..
	return(eBSFerrFileVer);
	}

return(eBSFSuccess);
}



teBSFrsltCodes
CBioSeqFile::Close(void)
{
return(Reset(true));
}

// Reset()
// Resets sequence file class context back to that immediately following instantiation
teBSFrsltCodes
CBioSeqFile::Reset(bool bFlush)		// true (default) is to write any pending header/directory/sequence to disk before closing opened file
{
teBSFrsltCodes Rslt = eBSFSuccess;
if(m_hFile != -1)
	{
	if(bFlush)
		Rslt = Flush2Disk();
	close(m_hFile);
	m_hFile = -1;
	}

if(m_pDirEntries != NULL)
	{
	delete m_pDirEntries;
	m_pDirEntries = NULL;
	}
if(m_ppDirEntryNames != NULL)
	{
	delete m_ppDirEntryNames;
	m_ppDirEntryNames = NULL;
	}
if(m_ppDirEntryIDs != NULL)
	{
	delete m_ppDirEntryIDs;
	m_ppDirEntryIDs = NULL;
	}
if(m_pDataBuff != NULL)
	{
	delete m_pDataBuff;
	m_pDataBuff = NULL;
	}
m_szFile[0] = '\0';	    			  	
m_bHdrDirty = false;					   
m_DataBuffLen = 0;
m_DataBuffNibs = 0;
m_AllocDirEntriesSize = 0;
m_bCreate = false;
m_DataWrtPsn = 0;
InitHdr(0);
return(Rslt);
}

void
CBioSeqFile::InitHdr(unsigned int Type) // biosequence file type
{
memset(&m_FileHdr,0,sizeof(tsBSFHeader));
m_FileHdr.Magic[0] = 'b';
m_FileHdr.Magic[1] = 'i';
m_FileHdr.Magic[2] = 'o';
m_FileHdr.Magic[3] = 's';
m_FileHdr.Type = Type;
m_FileHdr.Version = cBSFVersion;	        // file structure version
m_FileHdr.FileLen = sizeof(tsBSFHeader);	// current file length (nxt write psn)
strcpy(m_FileHdr.szDatasetName,"dataset1"); // default dataset name if none supplied
m_FileHdr.MaxEntries = cBSFMaxDirEntries;   // max number of entries
m_FileHdr.szDescription[0] = '\0';
m_FileHdr.szTitle[0] = '\0';
m_bHdrDirty = true; 
}

// Open()
// Open specified sequence file pszSeqFile with staging buffer of BuffSize
// Option to create or truncate pszSeqFile
teBSFrsltCodes 
CBioSeqFile::Open(char *pszSeqFile,	// specifies file to open or create
    			unsigned int Type,	// biosequence file type 
			   bool bCreate)		// create file or truncate if it already exists
{
teBSFrsltCodes Rslt;
if(pszSeqFile == NULL || *pszSeqFile == '\0') // validate parameters
	return(eBSFerrParams);

Reset(false);						// reset context in case file still opened

#ifdef _WIN32
if(!bCreate)
	m_hFile = open(pszSeqFile, O_READSEQ ); // file access is normally sequential..
else
	m_hFile = open(pszSeqFile, O_CREATETRUNC );
#else
if(!bCreate)
	m_hFile = open64(pszSeqFile, O_READSEQ ); // file access is normally sequential..
else
	if((m_hFile = open64(pszSeqFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=1)
		{
		if(ftruncate(m_hFile,0)!=0)
			{
			AddErrMsg("CBioSeqFile::Open","Unable to truncate %s - %s",pszSeqFile,strerror(errno));
			Rslt = bCreate ? eBSFerrCreateFile : eBSFerrOpnFile;
			return(Rslt);
			}
		}
#endif

if(m_hFile == -1)					// check if file open succeeded
	{
	AddErrMsg("CBioSeqFile::Open","Unable to open %s - %s",pszSeqFile,strerror(errno));
	Rslt = bCreate ? eBSFerrCreateFile : eBSFerrOpnFile;
	Reset(false);
	return(Rslt);
	}

strncpy(m_szFile,pszSeqFile,_MAX_PATH);
m_szFile[_MAX_PATH-1] = '\0';

if(bCreate)
	{
	InitHdr(Type);
	m_bHdrDirty = true;
	m_bCreate = true;
	if((Rslt = Flush2Disk()) != eBSFSuccess)
		{
		Reset(false);			// closes opened file..
		return(Rslt);
		}
	}
else // else opening existing file
	{
	if((Rslt=Disk2Hdr(pszSeqFile,Type))!=eBSFSuccess)
		{
		Reset(false);			// closes opened file..
		return(eBSFerrFileAccess);
		}

	// if not empty then load directory entries
	if(m_FileHdr.NumEntries)
		LoadEntries();
	}

return(eBSFSuccess);
}




teBSFrsltCodes 
CBioSeqFile::SetDescription(char *pszDescription)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(!m_bCreate)
	return(eBSFerrRead);
if(m_FileHdr.Version >= 5)
	{
	strncpy(m_FileHdr.szDescription,pszDescription,sizeof(m_FileHdr.szDescription));
	m_FileHdr.szDescription[sizeof(m_FileHdr.szDescription)-1] = '\0';
	return(eBSFSuccess);
	}
return(eBSFerrFileVer);
}

teBSFrsltCodes 
CBioSeqFile::GetDescription(int MaxLen,char *pszDescription)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(m_bCreate)
	return(eBSFerrWrite);
if(m_FileHdr.Version >= 5)
	strncpy(pszDescription,m_FileHdr.szDescription,MaxLen);
else
	strncpy(pszDescription,m_szFile,MaxLen);
pszDescription[MaxLen-1] = '\0';
return(eBSFSuccess);
}

teBSFrsltCodes 
CBioSeqFile::SetTitle(char *pszTitle)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(!m_bCreate)
	return(eBSFerrRead);
if(m_FileHdr.Version >= 5)
	{
	strncpy(m_FileHdr.szTitle,pszTitle,sizeof(m_FileHdr.szTitle));
	m_FileHdr.szTitle[sizeof(m_FileHdr.szTitle)-1] = '\0';
	return(eBSFSuccess);
	}
return(eBSFerrFileVer);
}

teBSFrsltCodes 
CBioSeqFile::GetTitle(int MaxLen,char *pszTitle)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(m_bCreate)
	return(eBSFerrWrite);
if(m_FileHdr.Version >= 5)
	strncpy(pszTitle,m_FileHdr.szTitle,MaxLen);
else
	{
	char szFname[_MAX_FNAME];
#ifdef _WIN32
	_splitpath(m_szFile,NULL,NULL,szFname,NULL);
#else
	CUtility::splitpath(m_szFile,NULL,szFname);
#endif
	strncpy(pszTitle,szFname,MaxLen);
	}
pszTitle[MaxLen-1] = '\0';
return(eBSFSuccess);
}

// LoadEntries
// Allocates memory for entry dictionary and indexes, then reads them in from disk
teBSFrsltCodes
CBioSeqFile::LoadEntries(void)
{
teBSFrsltCodes Rslt;
int Idx;
if(m_hFile == -1)
	return(eBSFerrClosed);
if(!m_FileHdr.NumEntries)		// any features to load?
	return(eBSFerrNoEntries);

if(m_pDirEntries != NULL)
	{
	delete m_pDirEntries;
	m_pDirEntries = NULL;
	m_AllocDirEntriesSize = 0;
	}
if(m_ppDirEntryNames != NULL)
	{
	delete m_ppDirEntryNames;
	m_ppDirEntryNames = NULL;
	}
if(m_ppDirEntryIDs != NULL)
	{
	delete m_ppDirEntryIDs;
	m_ppDirEntryIDs = NULL;
	}
// allocate memory to hold directory entries
if((m_pDirEntries = (tsBSFDirEntry *)new unsigned char [m_FileHdr.DirEntriesSize])==NULL)
	{
	m_AllocDirEntriesSize = 0;
	return(eBSFerrMem);
	}
m_AllocDirEntriesSize = m_FileHdr.DirEntriesSize;

if((m_ppDirEntryNames = (tsBSFDirEntry **)new tsBSFDirEntry *[m_FileHdr.NumEntries])==NULL)
	return(eBSFerrMem);

if((m_ppDirEntryIDs = (tsBSFDirEntry **)new tsBSFDirEntry *[m_FileHdr.NumEntries])==NULL)
	return(eBSFerrMem);

// all required memory has been allocated, read from disk
Rslt = ReadDisk(m_FileHdr.DirEntriesOfs,m_FileHdr.DirEntriesSize,m_pDirEntries);
if(Rslt == eBSFSuccess)
	{
	if(m_bIsBigEndian)	// if on a big-endian machine then need to make little endian as that is our native file format
		{
		tsBSFDirEntry *pDir = m_pDirEntries;
		UINT8 *pByte = (UINT8 *)pDir;
		for(Idx = 0; Idx < m_FileHdr.NumEntries; Idx++)
			{
			pDir = (tsBSFDirEntry *)pByte;
			pDir->DataPsn = SwapUI64Endians(pDir->DataPsn);	   // absolute file psn at which biosequence data starts
			pDir->Size= SwapUI32Endians(pDir->Size);			// size in bytes of this instance when concatenated
			pDir->EntryID= SwapUI32Endians(pDir->EntryID);		// 1..n uniquely identifies this entry instance
			pDir->NameInst= SwapUI32Endians(pDir->NameInst);	// name instance (1..n) -- can have multiple entries with duplicate names
			pDir->DataLen= SwapUI32Endians(pDir->DataLen);		// data length when unpacked
			pDir->Hash= SwapUI16Endians(pDir->Hash);			// hash on szName
			pByte += pDir->Size;
			}
		}
	}

if(Rslt == eBSFSuccess)
	Rslt = ReadDisk(m_FileHdr.DirEntryNameIdxOfs,m_FileHdr.NumEntries * sizeof(INT32),m_ppDirEntryNames);
if(Rslt == eBSFSuccess)
	Rslt = ReadDisk(m_FileHdr.DirEntryIDIdxOfs,m_FileHdr.NumEntries * sizeof(INT32),m_ppDirEntryIDs);
#pragma warning(push)
#pragma warning(disable: 4311)

INT32 *pInt32;
tsBSFDirEntry **ppEl;
pInt32 = (INT32 *)m_ppDirEntryNames;
pInt32 += m_FileHdr.NumEntries-1;
ppEl = &m_ppDirEntryNames[m_FileHdr.NumEntries-1];
for(Idx = m_FileHdr.NumEntries-1; Idx >= 0; Idx--,pInt32--,ppEl--)
	{
	if(m_bIsBigEndian)
			*pInt32 = SwapUI32Endians(*pInt32);
	 *ppEl = (tsBSFDirEntry *)((char *)m_pDirEntries + *pInt32);
	}

pInt32 = (INT32 *)m_ppDirEntryIDs;
pInt32 += m_FileHdr.NumEntries-1;
ppEl = &m_ppDirEntryIDs[m_FileHdr.NumEntries-1];
for(Idx = m_FileHdr.NumEntries-1; Idx >= 0; Idx--,pInt32--,ppEl--)
	{
	if(m_bIsBigEndian)
		*pInt32 = SwapUI32Endians(*pInt32);
	*ppEl = (tsBSFDirEntry *)((char *)m_pDirEntries + *pInt32);
	}

#pragma warning(pop)
return(Rslt);
}

// ReadDisk
// Reads block of size 'Len' from disk starting at 'DiskOfs' into preallocated memory at 'pTo'
teBSFrsltCodes
CBioSeqFile::ReadDisk(INT64 DiskOfs,int Len,void *pTo)
{
if(_lseeki64(m_hFile,DiskOfs,SEEK_SET)!=DiskOfs)
	{
	AddErrMsg("CBioSeqFile::ReadDisk","Seek failed on %s - %s",m_szFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}
if(read(m_hFile,pTo,Len)!=Len)
	{
	AddErrMsg("CBioSeqFile::ReadDisk","Read failed on %s - %s",m_szFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}
return(eBSFSuccess);
}



// GetType
// Returns file header type
int 
CBioSeqFile::GetType(void)
{
return(m_hFile == -1 ? eBSFerrClosed : m_FileHdr.Type);
}

// GenNameHash
// Generates a 16bit hash on specified lowercased name
// This hash can then be used to quickly eliminate probe names which can't match a target name by comparing hashes
unsigned short
CBioSeqFile::GenNameHash(char *pszName)
{
unsigned short Hash = 0;
unsigned short MSB;
char chr;
while(chr = *pszName++)
	{
	MSB = (Hash & 0x0c000) >> 13;
	Hash <<= 2;
	Hash ^= tolower(chr);
	Hash ^= MSB;
	}
return(Hash);
}




teBSFrsltCodes 
CBioSeqFile::Flush2Disk(void)						// flush and commit to disk
{
int WrtLen;
teBSFrsltCodes Rslt;
int *pOfs;
int Idx;

if(m_hFile != -1 && m_bCreate)		// if file opened for write
	{
	if(m_FileHdr.NumEntries && m_FileHdr.DirEntriesSize)
		{
		if((Rslt=SortEntries())!=eBSFSuccess)
			{
			AddErrMsg("CBioSeqFile::Flush2Disk","Unable to sort entries %s",m_szFile);
			Reset(false);
			return(Rslt);
			}

		if(m_bIsBigEndian)	// if on a big-endian machine then need to make little endian as that is our native file format
			{
			tsBSFDirEntry *pDir = m_pDirEntries;
			UINT8 *pByte = (UINT8 *)pDir;
			for(Idx = 0; Idx < m_FileHdr.NumEntries; Idx++)
				{
				pByte += pDir->Size;
				pDir->DataPsn = SwapUI64Endians(pDir->DataPsn);						    // absolute file psn at which biosequence data starts
				pDir->Size= SwapUI32Endians(pDir->Size);;								// size in bytes of this instance when concatenated
				pDir->EntryID= SwapUI32Endians(pDir->EntryID);							// 1..n uniquely identifies this entry instance
				pDir->NameInst= SwapUI32Endians(pDir->NameInst);							// name instance (1..n) -- can have multiple entries with duplicate names
				pDir->DataLen= SwapUI32Endians(pDir->DataLen);							// data length when unpacked
				pDir->Hash= SwapUI16Endians(pDir->Hash);							// hash on szName
				pDir = (tsBSFDirEntry *)pByte;
				}
			}

		WrtLen = m_FileHdr.DirEntriesSize;	
		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
			write(m_hFile,m_pDirEntries,WrtLen)!=WrtLen)
			{
			AddErrMsg("CBioSeqFile::Flush2Disk","Unable to write directory entries to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}

		m_FileHdr.DirEntriesOfs = m_FileHdr.FileLen;
		m_FileHdr.FileLen += WrtLen;

		// change index ptrs into offsets
		pOfs = (int *)m_ppDirEntryNames;
		for(Idx = 0; Idx < (int)m_FileHdr.NumEntries; Idx++,pOfs++)
			{
			*pOfs = (int)((char *)m_ppDirEntryNames[Idx] - (char *)m_pDirEntries);
			if(m_bIsBigEndian)
				*pOfs = SwapUI32Endians(*pOfs);
			}
		WrtLen = m_FileHdr.NumEntries * sizeof(tsBSFDirEntry *);	
		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
			write(m_hFile,m_ppDirEntryNames,WrtLen)!=WrtLen)
			{
			AddErrMsg("CBioSeqFile::Flush2Disk","Unable to write entry names index to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		m_FileHdr.DirEntryNameIdxOfs = m_FileHdr.FileLen;
		m_FileHdr.FileLen += WrtLen;

		pOfs = (int *)m_ppDirEntryIDs;
		for(Idx = 0; Idx < (int)m_FileHdr.NumEntries; Idx++,pOfs++)
			{
			*pOfs = (int)((char *)m_ppDirEntryIDs[Idx] - (char *)m_pDirEntries);
			if(m_bIsBigEndian)
				*pOfs = SwapUI32Endians(*pOfs);
			}

		WrtLen = m_FileHdr.NumEntries * sizeof(tsBSFDirEntry *);	
		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
			write(m_hFile,m_ppDirEntryIDs,WrtLen)!=WrtLen)
			{
			AddErrMsg("CBioSeqFile::Flush2Disk","Unable to write entry names index to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		m_FileHdr.DirEntryIDIdxOfs = m_FileHdr.FileLen;
		m_FileHdr.FileLen += WrtLen;
		}

		// now write the header to disk
	if((Rslt=Hdr2Disk())!= eBSFSuccess)
		return(Rslt);
	m_bHdrDirty = false;
	}
return(eBSFSuccess);
}

teBSFrsltCodes 
CBioSeqFile::SetDatasetName(char *pszDataset)			// sets file dataset name
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(pszDataset != NULL && pszDataset[0] != '\0')
	{
	strncpy(m_FileHdr.szDatasetName,pszDataset,cMaxDatasetSpeciesChrom);
	m_FileHdr.szDatasetName[cMaxDatasetSpeciesChrom-1] = '\0';
	return(eBSFSuccess);
	}
return(eBSFerrParams);
}

char *
CBioSeqFile::GetDatasetName(void)						// returns file dataset name
{
if(m_hFile == -1)
	return(NULL);
return(m_FileHdr.szDatasetName);
}

teBSFrsltCodes
CBioSeqFile::Exists(tBSFEntryID EntryID)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
return(LocateEntry(EntryID)!=NULL ? eBSFSuccess : eBSFerrEntry);
}

int
CBioSeqFile::NumEntries(void)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
return(m_FileHdr.NumEntries);
}

tsBSFDirEntry *
CBioSeqFile::LocateEntry(tBSFEntryID EntryID)
{
tsBSFDirEntry *pEntry;
if(m_hFile == -1)
	return(NULL);
if(m_bCreate)				// can't locate entries if file opened for create
	return(NULL);
if(EntryID < 1 || EntryID > (int)m_FileHdr.NumEntries)
	return(NULL);
pEntry = m_ppDirEntryIDs[EntryID-1];
return(pEntry->EntryID == EntryID ? pEntry : NULL);
}

// LocateEntryIDbyName
// Returns feature identifier for Ith feature having specified feature name as there could be multiple features with same name
int
CBioSeqFile::LocateEntryIDbyName(char *pszName,	// entry to locate
							int Ith)			// Ith instance to locate (1..n)
{
char *pszName2;
char c1;
char c2;
int CmpRslt;
tsBSFDirEntry *pProbe;
if(Ith < 1 || pszName == NULL || pszName[0] == '\0')
	return(eBSFerrParams);
if(m_hFile == -1)
	return(eBSFerrClosed);
if(m_bCreate)
	return(eBSFerrRead);

if( m_FileHdr.NumEntries == 0)
   return(eBSFerrEntry);

int Left = 0;
int Right = m_FileHdr.NumEntries - 1;
int MidPt;
c1 = tolower(*pszName);

while(Right >= Left) {
	MidPt = (Right + Left)/2;
	pProbe = m_ppDirEntryNames[MidPt];
	pszName2 = &pProbe->szName[0];
	c2 = tolower(*pszName2);
	if(c1 < c2)
		{
		Right = MidPt - 1;
		continue;
		}

	if(c1 > c2)
		{
		Left = MidPt + 1;
		continue;
		}

	CmpRslt = stricmp(pszName,pProbe->szName);
	if(CmpRslt < 0)
		{
		Right = MidPt - 1;
		continue;
		}
	else
		if(CmpRslt > 0)
			{
			Left = MidPt + 1;
			continue;
			}


	// have a match, but this match may not be match of required name instance (Ith)
	if(pProbe->NameInst == Ith)
		return(pProbe->EntryID);

	// do linear search forward/backwards until Ith instance of name located
	// performance shouldn't normally be an issue as many name duplications are not expected
	if(pProbe->NameInst > Ith)
		{	
		while(--MidPt >= Left) // linear search backwards
			{
			pProbe = m_ppDirEntryNames[MidPt];
			if(stricmp(pszName,pProbe->szName))
				break;
			
			if(pProbe->NameInst == Ith)
				return(pProbe->EntryID);

			if(pProbe->NameInst < Ith)
				break;
			}
		}
	else
		{	
		while(++MidPt <= Right) // // linear search forwards
			{
			pProbe = m_ppDirEntryNames[MidPt];
			if(stricmp(pszName,pProbe->szName))
				break;
			if(pProbe->NameInst == Ith)
				return(pProbe->EntryID);
			if(pProbe->NameInst > Ith)
				break;
			}
		}
	return((int)eBSFerrEntry);	// name exists but couldn't locate Ith instance
	}
return((int)eBSFerrEntry); // couldn't locate any instances of name
}


tBSFEntryID 
CBioSeqFile::Next(tBSFEntryID Current)		// returns next entry identifer after current
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(m_bCreate)				// can't locate entries if file opened for create
	return(eBSFerrRead);
if(Current >= m_FileHdr.NumEntries )
	return(eBSFerrEntry);

return(++Current);
}

tBSFEntryID 
CBioSeqFile::Prev(tBSFEntryID Current)		// returns previous entry identifer before current
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(m_bCreate)				// can't locate entries if file opened for create
	return(eBSFerrRead);
if(Current < 2)
	return(eBSFerrEntry);

if(Current > m_FileHdr.NumEntries)
	return(m_FileHdr.NumEntries);
return(--Current);
}


// returns entry name	
int 
CBioSeqFile::GetName(tBSFEntryID EntryID,int MaxLen,char *pszName)
{
int NameLen;
tsBSFDirEntry *pEntry;
if(m_hFile == -1)
	return(eBSFerrClosed);
if(m_bCreate)				// can't locate entries if file opened for create
	return(eBSFerrRead);
if(pszName == NULL || MaxLen < 1 )
	return(eBSFerrParams);
pEntry = LocateEntry(EntryID);
if(pEntry == NULL)
	return(eBSFerrEntry);

strncpy(pszName,pEntry->szName,MaxLen);
pszName[MaxLen-1] = '\0';
NameLen = (int)strlen(pszName);
return(NameLen);
}


// returns description	
int 
CBioSeqFile::GetDescription(tBSFEntryID EntryID,int MaxLen,char *pszDescr)
{
int NameLen;
int DescrLen;
tsBSFDirEntry *pEntry;
if(m_hFile == -1)
	return(eBSFerrClosed);
if(pszDescr == NULL || MaxLen < 1 )
	return(eBSFerrParams);
pEntry = LocateEntry(EntryID);
if(pEntry == NULL)
	return(eBSFerrEntry);

NameLen = (int)strlen(pEntry->szName);
if((UINT32)(NameLen + (int)sizeof(tsBSFDirEntry)) < pEntry->Size)
	{
	strncpy(pszDescr,&pEntry->szName[NameLen+1],MaxLen);
	pszDescr[MaxLen-1] = '\0';
	DescrLen = (int)strlen(pszDescr);
	}
else
	{
	*pszDescr = '\0';
	DescrLen = 0;
	}
return(DescrLen);
}


teBSFrsltCodes 
CBioSeqFile::GetNameDescription(tBSFEntryID EntryID,
								int MaxNameLen,
								char *pszName,
								int MaxDescrLen,
								char *pszDescr)
{
int NameLen = 0;
tsBSFDirEntry *pEntry = NULL;
if(m_hFile == -1)
	return(eBSFerrClosed);

if(pszName != NULL)
	{
	if(MaxNameLen < 1 )
		return(eBSFerrParams);
	pEntry = LocateEntry(EntryID);
	if(pEntry == NULL)
		return(eBSFerrEntry);
	strncpy(pszName,pEntry->szName,MaxNameLen);
	pszName[MaxNameLen-1] = '\0';
	NameLen = (int)strlen(pszName);
	}

if(pszDescr != NULL)
	{
	if(MaxDescrLen < 1)
		return(eBSFerrParams);

	if(pEntry == NULL)
		{
		pEntry = LocateEntry(EntryID);
		if(pEntry == NULL)
			return(eBSFerrEntry);
		}
	if(!NameLen)
		NameLen = (int)strlen(pEntry->szName);
	if((UINT32)(NameLen + (int)sizeof(tsBSFDirEntry)) < pEntry->Size)
		{
		strncpy(pszDescr,&pEntry->szName[NameLen+1],MaxDescrLen);
		pszDescr[MaxDescrLen-1] = '\0';
		}
	else
		*pszDescr = '\0';
	}
return(eBSFSuccess);
}


// returns unpacked data length	
UINT32 
CBioSeqFile::GetDataLen(tBSFEntryID EntryID)
{
tsBSFDirEntry *pEntry;
if(m_hFile == -1)
	return(0);

pEntry = LocateEntry(EntryID);
return(pEntry == NULL ? 0 : pEntry->DataLen);
}


etDataType
CBioSeqFile::GetDataType(tBSFEntryID EntryID)
{
tsBSFDirEntry *pEntry;
if(m_hFile == -1)			// if file closed then treat as
	return(eUndefDataType); // undefined datatype

pEntry = LocateEntry(EntryID);
if(pEntry == NULL)
	return(eUndefDataType);
return(pEntry->DType);
}


int 
CBioSeqFile::GetData(tBSFEntryID EntryID,
				etDataType ReqDataType,
				 UINT32 Ofs,				// start base offset (0..MaxLen-1)
				 unsigned char *pBuffer,	// where to return bases
				 UINT32 MaxLen)				// max length to return (1..n)
{
unsigned char DataBuff[0x03fff];
int BuffLen;
unsigned int FileSeqLen;
int CurReadFileSeqLen;
int BuffOfs;

INT64 DataPsn;
unsigned int ReadLen;
tsBSFDirEntry *pEntry;
if(m_hFile == -1)
	return(eBSFerrClosed);

if(pBuffer == NULL || MaxLen < 1)
	return(eBSFerrParams);
pEntry = LocateEntry(EntryID);
if(pEntry == NULL)
	return(eBSFerrEntry);

if(!pEntry->DataLen || Ofs >= pEntry->DataLen)
	return(eBSFerrOfs);

ReadLen = pEntry->DataLen - Ofs;
if(ReadLen > (int)MaxLen)
	ReadLen = MaxLen;
MaxLen = ReadLen;

if(pEntry->DType == eSeqBaseType) // eSeqBaseType are packed 2 per byte
	{
	DataPsn = pEntry->DataPsn + (Ofs / 2); // where to start disk read from
	FileSeqLen = (ReadLen+1)/2;		 // how many bytes to read from disk
	if(Ofs & 0x01 && !(ReadLen & 0x01))	 // if to start with hi nibble and end on low then
		FileSeqLen += 1;			 // need to read additional byte
	}
else
	{
	DataPsn = pEntry->DataPsn + Ofs;
	FileSeqLen = ReadLen;
	}

if(DataPsn != _lseeki64(m_hFile,DataPsn,SEEK_SET))
	{
	AddErrMsg("CBioSeqFile::GetData","Seek failed to offset %I64d on %s- %s",DataPsn,m_szFile,strerror(errno));
	return(eBSFerrFileAccess);
	}
if(ReqDataType == pEntry->DType)
	{
	BuffOfs = 0;
	while(FileSeqLen)
		{
		CurReadFileSeqLen = min(FileSeqLen,cMaxFileReadChunk);
		if(CurReadFileSeqLen != (BuffLen = read(m_hFile,&pBuffer[BuffOfs],CurReadFileSeqLen)))	
			{
			AddErrMsg("CBioSeqFile::GetData","Read(1) failed for %d bytes on %s- %s",CurReadFileSeqLen,m_szFile,strerror(errno));
			return(eBSFerrFileAccess);
			}
		BuffOfs += CurReadFileSeqLen;
		FileSeqLen -= CurReadFileSeqLen;
		}
	if(ReqDataType == eSeqBaseType) // if eSeqBaseType then an inplace unpack...
		{
		UnpackBases(ReadLen,			// unpacked sequence length 
					pBuffer,			// pts to sequence to unpack into
					Ofs & 0x01,			// nibble to start unpack from
					pBuffer);	
		}

	return(ReadLen);
	}

if(ReqDataType == eAsciiType && pEntry->DType == eSeqBaseType)
	{
	UINT32 NumCvrted = 0;
	while(NumCvrted < MaxLen && FileSeqLen)
		{
		CurReadFileSeqLen = min(FileSeqLen,sizeof(DataBuff));

		if(CurReadFileSeqLen != (BuffLen = read(m_hFile,DataBuff,CurReadFileSeqLen)))	
			{
			AddErrMsg("CBioSeqFile::GetData","Read(1) failed for %d bytes on %s- %s",CurReadFileSeqLen,m_szFile,strerror(errno));
			return(eBSFerrFileAccess);
			}

		BuffLen = min(MaxLen,(UINT32)CurReadFileSeqLen*2);
		CSeqTrans::MapPackedSeq2Ascii(Ofs,DataBuff,BuffLen,(char *)pBuffer);
		
		Ofs = 0;
		pBuffer += BuffLen;
		FileSeqLen -= CurReadFileSeqLen;
		NumCvrted += BuffLen;
		}
	return(MaxLen);
	}

if(ReqDataType == eSeqBaseType && pEntry->DType == eAsciiType)
	{
	while(ReadLen)
		{
		BuffLen = ReadLen > sizeof(DataBuff) ? sizeof(DataBuff) : ReadLen;
		if(BuffLen != read(m_hFile,&DataBuff,BuffLen))
			{
			AddErrMsg("CBioSeqFile::GetData","Read(3) failed for %d bytes on %s- %s",BuffLen,m_szFile,strerror(errno));
			return(eBSFerrFileAccess);
			}
		CSeqTrans::MapAscii2Sense((char *)DataBuff,BuffLen,(etSeqBase *)pBuffer);
		pBuffer += BuffLen;
		ReadLen -= BuffLen;
		}
	return(MaxLen);
	}
// can only be executing here if unable to convert to requested data type
return(eBSFerrCvrtType);
}


tBSFEntryID 
CBioSeqFile::CreateEntry(char *pszSource,
						 char *pszDescription,
						 etDataType DataType,
						 bool RptMskUpperCase)	// default is false, UCSC softmasked use lowercase when repeat masked
{
int Rslt;
int SourceLen;
int DescrLen;

tsBSFDirEntry *pEntry;
if(m_hFile == -1)
	return(eBSFerrClosed);
if(!m_bCreate)
	return(eBSFerrWrite);

if(pszSource == NULL || *pszSource == '\0')
	return(eBSFerrParams);

if((SourceLen = (int)strlen(pszSource)) > cBSFSourceSize)
	return(eBSFerrParams);

if(pszDescription == NULL || *pszDescription == '\0')
	DescrLen = 0;
else
	{
	if((DescrLen = (int)strlen(pszDescription)) > cBSFDescriptionSize)
		return(eBSFerrParams);
	DescrLen += 1;
	}
if(m_pCreatingEntry != NULL)			// seal any entry currently being created
	if((Rslt=SealEntry())!=eBSFSuccess)
		return(Rslt);

// Need to extend directory?
if(m_pDirEntries == NULL || 
   ((m_FileHdr.DirEntriesSize + (int)sizeof(tsBSFDirEntry) + SourceLen + DescrLen) > m_AllocDirEntriesSize))
	{
	pEntry = (tsBSFDirEntry *)new char [m_FileHdr.DirEntriesSize + cAllocEntriesIncr];
	if(pEntry == NULL)
		return(eBSFerrMem);
	if(m_pDirEntries)
		{
		if(m_FileHdr.DirEntriesSize)
			memmove(pEntry,m_pDirEntries,m_FileHdr.DirEntriesSize);
		delete m_pDirEntries;
		}
	m_pDirEntries = pEntry;
	m_AllocDirEntriesSize = m_FileHdr.DirEntriesSize + cAllocEntriesIncr;
	}
pEntry = (tsBSFDirEntry *)((char *)m_pDirEntries + m_FileHdr.DirEntriesSize);
pEntry->EntryID = ++m_FileHdr.NumEntries;
pEntry->Size = sizeof(tsBSFDirEntry) + SourceLen + DescrLen;
m_FileHdr.DirEntriesSize += pEntry->Size;
pEntry->Hash = GenNameHash(pszSource);
strcpy(pEntry->szName,pszSource);
if(DescrLen)
	strcpy(&pEntry->szName[SourceLen+1],pszDescription);
pEntry->DataLen = 0;			// data length
pEntry->DType = DataType;		// data type
pEntry->MskSense = RptMskUpperCase ? 1 : 0;
pEntry->DataPsn = 0;			// absolute file psn at which biosequence data starts
m_pCreatingEntry = pEntry;
m_DataBuffLen = 0;
m_DataBuffNibs = 0;
m_DataWrtPsn = 0;
m_bHdrDirty = true;
return(pEntry->EntryID);
}


// PackBases
// Pack bases (pUnpacked) into nibbles (pPacked) with bits 0..3 holding pUnpacked[n] and bits 4..7 holding pUnpacked[n+1]
// if pUnpacked pts to pPacked then inplace packing will be performed
int
CBioSeqFile::PackBases(unsigned int SeqLen,	// unpacked sequence length 
		  etSeqBase *pUnpacked,				// pts to sequence to pack from
		  unsigned int NibbleOfs,			// nibble to start pack into (0..SeqLen-1)
		  unsigned char *pPacked)			// where to pack
{
unsigned char byte;
bool bStartHi=false;	// true if first base is to be packed into hi nibble
bool bEndLo=false;		// true id last base is to be packed into lo nibble

if(pUnpacked == NULL || pPacked == NULL || !SeqLen)
	return(eBSFerrParams);	// error

if((SeqLen + NibbleOfs) & 0x01) // if last nibble to write is odd then only low nibble to be written
	{
	bEndLo = true;
	SeqLen -= 1;
	}

if(SeqLen && NibbleOfs & 0x01)// if start is odd then first base is to go into hi nibble
	{
	bStartHi = true;		
	SeqLen -= 1;
	}

pPacked += NibbleOfs/2;   // pt to first byte containing packed bases
SeqLen /= 2;
if(bStartHi)		// if first base is to go into hi nibble 
	{
	byte = *pPacked & 0x0f;	// preserve any existing lo nibble
	*pPacked++ = byte | ((*pUnpacked++ << 4) & 0x0f0); // combine
	}

while(SeqLen--)
	{
	*pPacked = *pUnpacked++ & 0x0f;
	*pPacked++ |= (*pUnpacked++ << 4) & 0x0f0;
	}

if(bEndLo)		// finish with lo nibble?
	{
	byte = *pPacked & 0x0f0;				// get any existing hi nibble and
	*pPacked = byte | (*pUnpacked & 0x0f); // combine with new lo nibble
	}
return(eBSFSuccess);
}


// UnpackBases
// unpack (packed into nibbles) from pPacked into pUnpack
// if pUnpacked pts to pPacked then inplace unpacking will be performed
int
CBioSeqFile::UnpackBases(unsigned int SeqLen,	// unpacked sequence length 
		  etSeqBase *pUnpacked,				// pts to sequence to unpack into
		  unsigned int NibbleOfs,			// nibble to start unpack from (0..1)
		  unsigned char *pPacked)			// where to unpack from
{
bool bEndHi=false;			// true if last base is to be unpacked from hi nibble
bool bStartLo=false;		// true if first base is to be unpacked from lo nibble

if(pUnpacked == NULL || pPacked == NULL || !SeqLen || NibbleOfs > 1)
	return(eBSFerrParams);	// error

pUnpacked += SeqLen - 1;
pPacked += (SeqLen + NibbleOfs - 1) / 2;

if(NibbleOfs & 0x01)	// if starting on high nibble
	{
	bEndHi = true;		// remember - reverse processing so will handle at end			
	if(!(SeqLen & 0x01))// if ending on low nibble
	   bStartLo = true; // remember - reverse processing so will handle at start
	else
	   bStartLo = false;
	}
else					// else starting with low nibble
	{	
	bEndHi = false;		// remember - reverse processing so will handle at end			
	if(SeqLen & 0x01)	// if ending on low nibble
	   bStartLo = true; // remember - reverse processing so will handle at start
	else
	  bStartLo = false;
	}

SeqLen = (SeqLen - NibbleOfs)/2;

if(bStartLo)				// if first base is from lo nibble 
	*pUnpacked-- = *pPacked-- & 0x0f; 

while(SeqLen--)
	{
	*pUnpacked-- = (*pPacked >> 4) & 0x0f;
	*pUnpacked-- = *pPacked-- & 0x0f;
	}

if(bEndHi)		// finish with hi nibble?
	*pUnpacked = (*pPacked >> 4) & 0x0f; 

return(eBSFSuccess);
}


// AddData
// Add data to specified entry
// NOTE - can only add data to currently created entry
teBSFrsltCodes 
CBioSeqFile::AddData(unsigned int DataLen, unsigned char *pData)
{
unsigned int BlockLen;
unsigned int PackedLen;
if(m_pCreatingEntry == NULL || pData == NULL || !DataLen)
	return(eBSFerrParams);
	
if(m_pDataBuff == NULL)
	{
	m_pDataBuff = new unsigned char [cBSFDataBuffSize + 1];	// allow for additional low nibble write
	if(m_pDataBuff == NULL)
		{
		AddErrMsg("CBioSeqFile::AddData","Memory allocation of %d bytes for %s- %s",cBSFDataBuffSize+1,m_szFile,strerror(errno));
		return(eBSFerrMem);
		}
	m_DataBuffLen = 0;
	m_DataBuffNibs = 0;
	}

m_bHdrDirty = true;
if(m_DataWrtPsn == 0)
	{
	m_DataWrtPsn = m_FileHdr.FileLen;
	m_pCreatingEntry->DataPsn = m_FileHdr.FileLen;
	}

// can pack 2x if bases added...
if(m_pCreatingEntry->DType == eSeqBaseType)
	{
	PackedLen = (DataLen+1)/2;
	BlockLen = cBSFDataBuffSize * 2;
	}
else
	{
	PackedLen = DataLen;
	BlockLen = cBSFDataBuffSize;
	}

// if too much data to buffer then process as blocks
while(PackedLen > cBSFDataBuffSize)
	{
	AddData(BlockLen,pData);
	DataLen -= BlockLen;
	if(!DataLen)
		return(eBSFSuccess);
	if(m_pCreatingEntry->DType == eSeqBaseType)
		PackedLen = (DataLen+1)/2;
	else
		PackedLen = DataLen;
	pData += BlockLen;
	}

// if can't fit into remaining buffer then write current buffer to disk
if(m_DataBuffLen >= cBSFDataBuffSize || (m_DataBuffLen >= 1 && (PackedLen > (cBSFDataBuffSize - m_DataBuffLen))))
	{
	if(m_DataBuffNibs & 0x01)		// write out complete bytes with both nibbles containing bases
		m_DataBuffLen -= 1;
	if(m_DataWrtPsn != _lseeki64(m_hFile,m_DataWrtPsn,SEEK_SET))
		{
		AddErrMsg("CBioSeqFile::AddData","Seek to %I64d on %s - %s",m_DataWrtPsn,m_szFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	if((BlockLen = write(m_hFile,m_pDataBuff,m_DataBuffLen))!= m_DataBuffLen)
		{
		AddErrMsg("CBioSeqFile::AddData","Write of %d bytes on %s - %s",m_DataBuffLen,m_szFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	m_DataWrtPsn += m_DataBuffLen;
	if(m_DataBuffNibs & 0x01)		// if last low nibble not written then make it the first
		{
		m_pDataBuff[0] = m_pDataBuff[m_DataBuffLen];
		m_DataBuffLen = 1;
		m_DataBuffNibs = 1;
		}
	else
		{
		m_DataBuffNibs = 0;
		m_DataBuffLen = 0;
		}
	}

if(m_pCreatingEntry->DType == eSeqBaseType)
	{
	PackBases(DataLen,pData,m_DataBuffNibs,m_pDataBuff);
	m_DataBuffNibs += DataLen;
	m_DataBuffLen = (m_DataBuffNibs+1)/2;
	}
else
	{
	memmove(&m_pDataBuff[m_DataBuffLen],pData,DataLen);
	m_DataBuffLen += DataLen;
	m_DataBuffNibs = m_DataBuffLen * 2;
	}

m_pCreatingEntry->DataLen += DataLen;
return(eBSFSuccess);
}



teBSFrsltCodes
CBioSeqFile::SealEntry(void)
{
int Len;
if(m_pCreatingEntry == NULL)			// error if entry not currently being created..
	return(eBSFerrEntryCreate);

if(m_DataBuffLen != 0)
	{
	// write out data sitting in buffer
	if(m_DataWrtPsn != _lseeki64(m_hFile,m_DataWrtPsn,SEEK_SET))
		{
		AddErrMsg("CBioSeqFile::SealEntry","Seek to %I64d on %s - %s",m_DataWrtPsn,m_szFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	
	Len = write(m_hFile,m_pDataBuff,m_DataBuffLen);
	if(Len != m_DataBuffLen)
		{
		AddErrMsg("CBioSeqFile::SealEntry","Write of %d bytes on %s - %s",m_DataBuffLen,m_szFile,strerror(errno));
		return(eBSFerrFileAccess);
		}

	m_DataWrtPsn += m_DataBuffLen;
	}

if(m_DataWrtPsn != 0)
	{
	if(m_FileHdr.SeqOfs == 0)
		m_FileHdr.SeqOfs = m_FileHdr.FileLen;
	m_FileHdr.FileLen = m_DataWrtPsn;
	m_FileHdr.SeqSize = m_DataWrtPsn - m_FileHdr.SeqOfs;
	m_DataWrtPsn = 0;
	}

m_DataBuffLen = 0;
m_DataBuffNibs = 0;
m_pCreatingEntry = NULL;
return(eBSFSuccess);
}


teBSFrsltCodes
CBioSeqFile::SortEntries(void)
{
int EntryID;
int NameInst;
tsBSFDirEntry *pEntry;
tsBSFDirEntry *pPrevEntry;

if(!m_FileHdr.NumEntries || m_pDirEntries == NULL)
	return(eBSFerrEntry);

if(m_ppDirEntryNames != NULL)
	{
	delete m_ppDirEntryNames;
	m_ppDirEntryNames = NULL;
	}
if(m_ppDirEntryIDs != NULL)
	{
	delete m_ppDirEntryIDs;
	m_ppDirEntryIDs = NULL;
	}

// create array of ptrs into concatenated entries which will be ordered by EntryID
m_ppDirEntryIDs = new tsBSFDirEntry * [m_FileHdr.NumEntries];
if(m_ppDirEntryIDs == NULL)
	{
	delete m_ppDirEntryNames;
	m_ppDirEntryNames = NULL;
	return(eBSFerrMem);
	}

pEntry = m_pDirEntries;
for(EntryID = 0; EntryID < m_FileHdr.NumEntries; EntryID++)
	{
	m_ppDirEntryIDs[EntryID] = pEntry;
	pEntry = (tsBSFDirEntry *)(((char *)pEntry) + pEntry->Size);
	}

// finally make the entry ids into ascending order
for(EntryID = 1; EntryID <= m_FileHdr.NumEntries; EntryID++)
	{
	pEntry = m_ppDirEntryIDs[EntryID-1];
	pEntry->EntryID = EntryID;
	}

// create array of ptrs into concatenated features which will be used to sort
// by name then EntryID
m_ppDirEntryNames = new tsBSFDirEntry * [m_FileHdr.NumEntries];
if(m_ppDirEntryNames == NULL)
	return(eBSFerrMem);

pEntry = m_pDirEntries;
for(EntryID = 0; EntryID < m_FileHdr.NumEntries; EntryID++)
	{
	m_ppDirEntryNames[EntryID] = pEntry;
	pEntry = (tsBSFDirEntry *)(((char *)pEntry) + pEntry->Size);
	}
qsort(m_ppDirEntryNames,m_FileHdr.NumEntries,sizeof(tsBSFDirEntry *),SortEntryNames);

// now that they are sorted by name then determine the number of name instances
NameInst = 1;
pPrevEntry = m_ppDirEntryNames[0];
pPrevEntry->NameInst = NameInst++;
for(EntryID = 1; EntryID < m_FileHdr.NumEntries; EntryID++)
	{
	pEntry = m_ppDirEntryNames[EntryID];
	if(pPrevEntry->Hash != pEntry->Hash ||
		stricmp(pEntry->szName,pPrevEntry->szName))
		NameInst = 1;
	pEntry->NameInst = NameInst++;
	pPrevEntry = pEntry;
	}

return(eBSFSuccess);
}

//LocateBase
// Locates the first instance of the specified base and returns the position (0..TargLen-1) of that base
// Returns -1 if no instance of the specified base located
// Note that bases are compared independent of any repeat masking flag - cRptMskFlg
int
CBioSeqFile::LocateBase(etSeqBase Probe, int TargLen,etSeqBase *pTarget)
{
int Psn;
etSeqBase Unmasked = Probe & ~cRptMskFlg;
for(Psn = 0 ; Psn < TargLen; Psn++)
	if((*pTarget++ &  ~cRptMskFlg) == Unmasked)
		return(Psn);
return(-1);
}


int
CBioSeqFile::LocateSequence(int ProbeLen, etSeqBase *pProbe, int TargLen,etSeqBase *pTarget)
{
int TargPsn;
int ProbePsn;
etSeqBase *pP;
etSeqBase *pT;
if(!ProbeLen || !TargLen || ProbeLen > TargLen ||
   pTarget == NULL || pProbe == NULL)
   return(-3);

etSeqBase Unmasked = *pProbe & ~cRptMskFlg;
for(TargPsn = 0; TargPsn < (TargLen - ProbeLen); TargPsn++,pTarget++)
	{
	if(Unmasked == (*pTarget & ~cRptMskFlg))
		{
		if(ProbeLen == 1)
			return(TargPsn);
	
		// have an anchor..
		pT = pTarget+1;
		pP = pProbe+1;
		for(ProbePsn = 1; ProbePsn < ProbeLen; ProbePsn++,pT++,pP++)
			{
			if((*pT & ~cRptMskFlg) != (*pP & ~cRptMskFlg))
				break;
			}

		if(ProbePsn == ProbeLen)
			return(TargPsn);
		}
	}
return(-1);
} 

// SortEntryNames
// Used to sort by entry names -> EntryID
int 
CBioSeqFile::SortEntryNames( const void *arg1, const void *arg2)
{
int Rslt;
tsBSFDirEntry *pEl1 = *(tsBSFDirEntry **)arg1;
tsBSFDirEntry *pEl2 = *(tsBSFDirEntry **)arg2;

char *pszName1 = &pEl1->szName[0];
char *pszName2 = &pEl2->szName[0];
char c1 = tolower(*pszName1);
char c2 = tolower(*pszName2);
if(c1 < c2)
	return(-1);
if(c1 > c2)
	return(1);
if((Rslt = stricmp(pszName1,pszName2))!=0)
	return(Rslt);
if(pEl1->EntryID < pEl2->EntryID)
	return(-1);
if(pEl1->EntryID > pEl2->EntryID)
	return(1);
return(0);
}



// GetShuffledBases
// Gets the aligned sequence and then shuffles it whilst preserving dinucleotide frequencies
// Returns shuffled dinucleotide values as etSeqBases
// Note that dinucleotide frequencies are preserved within each aligned segment
int 
CBioSeqFile::GetShuffledBases(int ChromID,	// reference chromosome identifier
					  int ChromOfs,			// start on reference chromosome (0..n) 
					  int NumPts,			// number of points to return
					  unsigned char *pToRet,	// where to return etSeqBases
					  etSeqBase MissingMarker)	// value to use if missing bases
{
int SegStart;
int SegLen;
etSeqBase *pSeg;
int Psn;
int Len;

Len = GetData(ChromID,eSeqBaseType,ChromOfs,pToRet,NumPts);
if(Len == eBSFerrOfs)
	Len = 0;
else
	if(Len < 0)
		return(Len);

if(Len > 0)
	{
	srand((unsigned int)ChromOfs);
	SegStart = SegLen = 0;
	pSeg = pToRet;
	for(Psn =  0; Psn < Len; Psn++,pSeg++)
		{
		if(*pSeg == MissingMarker || (*pSeg & ~cRptMskFlg) > eBaseT)
			{
			if(SegLen)			// segment started?
				ShuffleDinucs(&pToRet[SegStart-1],SegLen);
			SegLen = 0;
			SegStart = Psn+1;
			}
		else
			{
			*pSeg = *pSeg & ~cRptMskFlg; // can only shuffle unmasked 
			SegLen++;
			}
		}
	ShuffleDinucs(&pToRet[SegStart],SegLen);
	}

if(Len < NumPts)
	memset(&pToRet[Len],eBaseUndef,NumPts - Len);
return(eBSFSuccess);
}

// ShuffleDinucs
// a) If SeqLen < 30 then simple random swaps are used which preserves base frequencies only
// b) If SeqLen >= 30 then 1st order markov used which preserves dinucleotide frequencies
// 
void
CBioSeqFile::ShuffleDinucs(etSeqBase *pSeq,int SeqLen)
{
int Psn;
etSeqBase Base;
int TotDinucs = 0;
int Accum;

// if number of dinucleotides < 500 then just do simple random swaps
// as otherwise dinucleotide shuffling has too constraints and generated sequence
// can end up being too close to the orginal
if(SeqLen < 500)	
	{
	for(Accum=0; Accum < SeqLen; Accum++)
		{
		Base = pSeq[Accum];
		Psn = rand() % SeqLen;
		pSeq[Accum] = pSeq[Psn];
		pSeq[Psn] = Base;
		}
	return;
	}
CShuffle Shuffle;
Shuffle.SeqDPShuffle(SeqLen,pSeq,pSeq);
}


// dump entries as XML - actual sequence dump is optional
bool
CBioSeqFile::Dump2XML(char *pszXMLfile, 
					  unsigned int MaxDumpSeqLen)
{
int hDumpXMLFile;
char szBuff[1000];
char szSource[cBSFSourceSize+1];
char szDescr[cBSFDescriptionSize+1];
unsigned int EntryID;
unsigned int DataType;
unsigned int DataLen;
int Len;

unsigned char *pSeqBuff = NULL;
unsigned int MaxAllocdSeqLen = 0;
unsigned int DumpSeqLen = 0;

#ifdef _WIN32
hDumpXMLFile = open(pszXMLfile,O_CREATETRUNC );
#else
if((hDumpXMLFile = open(pszXMLfile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	{
	if(ftruncate(hDumpXMLFile,0)!=0)
		{
		AddErrMsg("CBioSeqFile::Dump2XML","Unable to truncate %s - %s",pszXMLfile,strerror(errno));
		return(false);
		}
	}
#endif

if(hDumpXMLFile < 0)
	{
	AddErrMsg("CBioSeqFile::Dump2XML","Unable to create or truncate dump file %s error: %s",pszXMLfile,strerror(errno));
	return(false);
	}
Len = sprintf(szBuff,"<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n");
Len += sprintf(&szBuff[Len],"<dataroot xmlns:od=\"urn:schemas-microsoft-com:officedata\" ");
Len += sprintf(&szBuff[Len],"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" ");
Len += sprintf(&szBuff[Len],"xsi:noNamespaceSchemaLocation=\"bioseqfile.xsd\" generated=\"2004-04-25T15:30:02\" >");
if(write(hDumpXMLFile,szBuff,Len)!=Len)
	{
	AddErrMsg("CBioSeqFile::Dump2XML","Write to file %s failed with error: %s",pszXMLfile,strerror(errno));
	return(false);
	}

// iterate over each entry
EntryID = 0;
while((EntryID = Next(EntryID))!=0)
	{
	GetNameDescription(EntryID,sizeof(szSource)-1,szSource,sizeof(szDescr)-1,szDescr);
	DataLen = GetDataLen(EntryID);
	DataType = GetDataType(EntryID);
	Len = sprintf(szBuff,"\n<entry>\n<entryid>%d</entryid>\n<source>%s</source>\n<descriptor>%s</descriptor>\n<length>%u</length>\n<type>%d</type>\n",
			EntryID,szSource,szDescr,DataLen,DataType);
	if(DataLen > 0 && MaxDumpSeqLen && DataType == cBSFTypeSeq)
		{
		DumpSeqLen = DataLen > MaxDumpSeqLen ? MaxDumpSeqLen : DataLen;
		if(pSeqBuff == NULL || MaxAllocdSeqLen < (DumpSeqLen * 2) + 1)
			{
			if(pSeqBuff != NULL)
				delete pSeqBuff;
			MaxAllocdSeqLen = (DumpSeqLen * 2) + 10000; // note that additional is allocated to reduce prob of a later alloc being required
			pSeqBuff = new unsigned char [MaxAllocdSeqLen];
			if(pSeqBuff == NULL)
				{
				close(hDumpXMLFile);
				return(false);
				}
			}
		GetData(EntryID,DataType,0,pSeqBuff,DumpSeqLen);
		CSeqTrans::MapSeq2Ascii(pSeqBuff,DumpSeqLen,(char *)&pSeqBuff[DumpSeqLen]);
		Len = sprintf(szBuff,"\n<seq>%s</seq>",&pSeqBuff[DumpSeqLen]);
		if(write(hDumpXMLFile,szBuff,Len)!=Len)
			{
			AddErrMsg("CBioSeqFile::Dump2XML","Write to file %s failed with error: %s",pszXMLfile,strerror(errno));
			return(false);
			}
		Len = 0;
		}
	Len += sprintf(&szBuff[Len],"</entry>");
	if(write(hDumpXMLFile,szBuff,Len)!=Len)
		{
		AddErrMsg("CBioSeqFile::Dump2XML","Write to file %s failed with error: %s",pszXMLfile,strerror(errno));
		return(false);
		}
	}
Len = sprintf(szBuff,"\n</dataroot>\n");
if(write(hDumpXMLFile,szBuff,Len)!=Len)
	AddErrMsg("CBioSeqFile::Dump2XML","Write to file %s failed with error: %s",pszXMLfile,strerror(errno));
close(hDumpXMLFile);
if(pSeqBuff != NULL)
	delete pSeqBuff;
return(true);
}

int 
CBioSeqFile::GetNumMasked(int ProbeLen, etSeqBase *pProbe)
{
int NumMasked = 0;
etSeqBase Base;
while(ProbeLen-- && (Base = *pProbe++) != eBaseEOS)
	if(Base & cRptMskFlg)
		NumMasked++;
return(NumMasked);
}


// ClassComplexity
// Attempt to classify the complexity of a sequence on the basis of base frequency
// 3 == all bases are represented
// 2 == 1 base is under representated as less than ThresPercent of Len
// 1 == 2 bases are grossly underrepresentated
// 0 == 3 bases are grossly underrepresentated
unsigned int
CBioSeqFile::ClassComplexity(etSeqBase *pSeq,
							 unsigned int Len,
							 unsigned int ThresPercent) // bases below this are under representated
{
unsigned int BaseCnts[4];				// used to hold counts of individual nucleotides
unsigned int UnderRep = (Len * ThresPercent) / 100;
unsigned int Cnt;
unsigned int Ofs;

if(pSeq == NULL || ThresPercent > 99 || Len < 10)
	return(true);

memset(BaseCnts,0,sizeof(BaseCnts));
for(Cnt = 0; Cnt < Len; Cnt++)
	{
	Ofs = *pSeq++ & ~cRptMskFlg;
	if(Ofs >= eBaseN)					// any eBaseN's? treat as though it's a low complexity sequence
		return(0);
	BaseCnts[Ofs & 0x03]++;				// update counts of individual nucleotides
	}

// now analyse the number of bases under represented
int UnderRepCnt = 0;
for(Cnt = 0; Cnt < 4; Cnt++)
	{
	if(BaseCnts[Cnt] < UnderRep)
		UnderRepCnt++;
	}
return(3-UnderRepCnt);
}



unsigned int
CBioSeqFile::ScoreComplexity(etSeqBase *pSeq,
							 unsigned int Len)
							 
{
int TriFreqCnts[64+1];		// used to hold counts of trimers plus a end marker
unsigned char Base;
int CntA;
int *pCntA;

#ifdef USETHISCODE
int CntB;
int *pCntB;
#endif

unsigned int BaseIdx;

int HighCnts[10];
unsigned int Score;

if(pSeq == NULL || Len < 6)
	return(0);

memset(TriFreqCnts,0,sizeof(TriFreqCnts));
BaseIdx = 0;
Base = *pSeq++ & ~cRptMskFlg;
if(Base < eBaseN)
	{
	BaseIdx = Base;
	BaseIdx <<= 2;
	}
Base = *pSeq++ & ~cRptMskFlg;
if(Base < eBaseN)
	BaseIdx |= Base;
for(CntA = 2; CntA < (int)Len; CntA++)
	{
	BaseIdx <<= 2;
	BaseIdx &= 0x03c;
	Base = *pSeq++ & ~cRptMskFlg;
	if(Base >= eBaseN)					// slough eBaseN's
		continue;
	BaseIdx |= Base;
	TriFreqCnts[BaseIdx]++;				// update counts of individual nucleotides
	}


// group all non-zero counts together
TriFreqCnts[64] = -1;					// mark the end of the freq cnts
#ifdef USETHISCODE
pCntA = TriFreqCnts;
while((CntA = *pCntA++) != -1)
	{
	if(CntA == 0)
		{
		pCntB = pCntA;
		while((CntB = *pCntB++)== 0);
		if(CntB == -1)
			break;
		pCntA[-1] = CntB;
		pCntA = pCntB;
		}
	}
#endif

// find highest counts 
pCntA = TriFreqCnts;
HighCnts[0] = HighCnts[1] = HighCnts[2] = 0;
while((CntA = *pCntA++) != -1)
	{
	if(!CntA)
		continue;
	if(CntA >= HighCnts[0])
		{
		HighCnts[2] = HighCnts[1];
		HighCnts[1] = HighCnts[0];
		HighCnts[0] = CntA;
		continue;
		}
	if(CntA >= HighCnts[1])
		{
		HighCnts[2] = HighCnts[1];
		HighCnts[1] = CntA;
		continue;
		}
	if(CntA > HighCnts[2])
		{
		HighCnts[2] = CntA;
		continue;
		}
	}
Score = HighCnts[0] + HighCnts[1] + HighCnts[2];
return((Score * 100) / (Len - 2));
}

// QuickAlignRight
// An anchored quick alignment which can handle a limited number of mismatches and InDels
// Returns the maximal number of exactly alignable bases
// Treats indeterminate bases 'N' as always matching but does not include them in the returned score
unsigned int								// number of aligned bases
CBioSeqFile::QuickAlignRight(unsigned int AllowedMismatches,	// total allowed mismatches
		   unsigned int MaxMismatchSeqLen,	// max number of mismatches allowed in any run 
		   unsigned int AllowedProbeInDels,	// total allowed indel events on probe
		   unsigned int AllowedTargInDels,	// total allowed indel events on target
		   unsigned int MaxInDelExtension,	// max length of any InDel extension 
		   unsigned int ProbeLen,			// remaining probe length
		   unsigned char *pProbe,
   		   unsigned int TargLen,			// remaining target length
		   unsigned char *pTarg)
{
unsigned int Score = 0;
unsigned int MismatchSeqLen = 0;
unsigned int Psn;
unsigned int NewScore = 0;
unsigned int IndelScore = 0;
unsigned int InDelPsn;
unsigned int LastInDelPsn;
etSeqBase UnmaskedProbe;
etSeqBase UnmaskedTarg;

for(Psn = 0; Psn < ProbeLen &&	Psn < TargLen; Psn++)
	{
	UnmaskedProbe = pProbe[Psn] & ~cRptMskFlg;
	UnmaskedTarg = pTarg[Psn] & ~cRptMskFlg; 
	if(UnmaskedProbe < eBaseN  && UnmaskedProbe == UnmaskedTarg)
		{
		Score++;
		MismatchSeqLen = 0;
		}
	else
		{
		if(UnmaskedProbe >= eBaseN || UnmaskedTarg >= eBaseN)
			continue;

		// too many mismatches or current sequence of mismatches too long?
		if(!AllowedMismatches-- || MismatchSeqLen++ >= MaxMismatchSeqLen)	
			break;

		// if first mismatch in sequence then it maybe worth exploring indel events if next few bases also mismatch..
		if(MismatchSeqLen == 1)			
			{
			if(ProbeLen < 10 || TargLen < 10 || ((pProbe[Psn+1] & ~cRptMskFlg) == (pTarg[Psn+1] & ~cRptMskFlg) &&
				(pProbe[Psn+2] &  ~cRptMskFlg) == (pTarg[Psn+2] &  ~cRptMskFlg)))
				continue;

			if(AllowedProbeInDels)
				{
				AllowedProbeInDels--;
				
				// lookahead and see if following an InDel event an alignment of at least 3nt can be started...
				InDelPsn = Psn + 1;
				LastInDelPsn = InDelPsn + MaxInDelExtension;
				if(LastInDelPsn > ProbeLen - 7)
					LastInDelPsn = ProbeLen - 7;

				for(;InDelPsn <= LastInDelPsn; InDelPsn++)
					if((pProbe[InDelPsn] &  ~cRptMskFlg) == (pTarg[Psn] &  ~cRptMskFlg) &&
						(pProbe[InDelPsn+1] &  ~cRptMskFlg) == (pTarg[Psn+1] &  ~cRptMskFlg) &&
						(pProbe[InDelPsn+2] &  ~cRptMskFlg) == (pTarg[Psn+2] &  ~cRptMskFlg))
						break;

				if(InDelPsn <= LastInDelPsn)
					{
					NewScore = 3 + Score + QuickAlignRight(AllowedMismatches, 
							MaxMismatchSeqLen,
							AllowedProbeInDels,
							AllowedTargInDels,
							MaxInDelExtension,
							ProbeLen - (InDelPsn + 3),
							&pProbe[InDelPsn + 3],
							TargLen - (Psn + 3), 
							&pTarg[Psn + 3]);

					if(NewScore > IndelScore)
						IndelScore = NewScore;
					}
				}

			if(AllowedTargInDels)
				{
				AllowedTargInDels--;
				
				// lookahead and see if following an InDel event an alignment of at least 3nt can be started...
				InDelPsn = Psn + 1;
				LastInDelPsn = InDelPsn + MaxInDelExtension;
				if(LastInDelPsn > TargLen - 7)
					LastInDelPsn = TargLen - 7;

				for(;InDelPsn <= LastInDelPsn; InDelPsn++)
					if((pProbe[Psn]  &  ~cRptMskFlg) == (pTarg[InDelPsn] &  ~cRptMskFlg) &&
						(pProbe[Psn+1] &  ~cRptMskFlg) == (pTarg[InDelPsn+1] &  ~cRptMskFlg) &&
						(pProbe[Psn+2] &  ~cRptMskFlg) == (pTarg[InDelPsn+2] &  ~cRptMskFlg))
						break;

				if(InDelPsn <= LastInDelPsn)
					{
					NewScore = 3 + Score + QuickAlignRight(AllowedMismatches, 
							MaxMismatchSeqLen,
							AllowedProbeInDels,
							AllowedTargInDels,
							MaxInDelExtension,
							ProbeLen - (Psn + 3),
							&pProbe[Psn + 3],
							TargLen - (InDelPsn + 3), 
							&pTarg[InDelPsn + 3]);

					if(NewScore > IndelScore)
						IndelScore = NewScore;
					}
				}
			}
		}
	}
return(Score > IndelScore ? Score : IndelScore);
}


unsigned int								// number of aligned bases
CBioSeqFile::QuickAlignLeft(unsigned int AllowedMismatches,	// total allowed mismatches
		   unsigned int MaxMismatchSeqLen,	// max number of mismatches allowed in any run 
		   unsigned int AllowedProbeInDels,	// total allowed indel events on probe
		   unsigned int AllowedTargInDels,	// total allowed indel events on target
		   unsigned int MaxInDelExtension,	// max length of any InDel extension 
		   unsigned int ProbeLen,			// remaining probe length
		   unsigned char *pProbe,
   		   unsigned int TargLen,			// remaining target length
		   unsigned char *pTarg)
{
unsigned int Score = 0;
unsigned int MismatchSeqLen = 0;
unsigned int Psn;
unsigned int NewScore = 0;
unsigned int IndelScore = 0;
unsigned int InDelPsn;
unsigned int LastInDelPsn;
etSeqBase UnmaskedProbe;
etSeqBase UnmaskedTarg;

for(Psn = 0; Psn < ProbeLen &&	Psn < TargLen; Psn++)
	{
	UnmaskedProbe = pProbe[-(int)Psn] & ~cRptMskFlg;
	UnmaskedTarg = pTarg[-(int)Psn] & ~cRptMskFlg; 
	if(UnmaskedProbe < eBaseN  && UnmaskedProbe == UnmaskedTarg)
		{
		Score++;
		MismatchSeqLen = 0;
		}
	else
		{
		if(UnmaskedProbe >= eBaseN || UnmaskedTarg >= eBaseN)
			continue;

		// too many mismatches or current sequence of mismatches too long?
		if(!AllowedMismatches-- || MismatchSeqLen++ >= MaxMismatchSeqLen)	
			break;

		// if first mismatch in sequence then it maybe worth exploring indel events if next few bases also mismatch..
		if(MismatchSeqLen == 1)			
			{
			if(ProbeLen < 10 || TargLen < 10 || ((pProbe[-(int)(Psn+1)] &  ~cRptMskFlg) == (pTarg[-(int)(Psn+1)] &  ~cRptMskFlg) &&
				(pProbe[-(int)(Psn+2)]  &  ~cRptMskFlg) == (pTarg[-(int)(Psn+2)] &  ~cRptMskFlg)))
				continue;

			if(AllowedProbeInDels)
				{
				AllowedProbeInDels--;
				
				// lookahead and see if following an InDel event an alignment of at least 3nt can be started...
				InDelPsn = Psn + 1;
				LastInDelPsn = InDelPsn + MaxInDelExtension;
				if(LastInDelPsn > ProbeLen - 7)
					LastInDelPsn = ProbeLen - 7;

				for(;InDelPsn <= LastInDelPsn; InDelPsn++)
					if((pProbe[-(int)InDelPsn]  &  ~cRptMskFlg) == (pTarg[-(int)Psn]  &  ~cRptMskFlg) &&
						(pProbe[-(int)(InDelPsn+1)] &  ~cRptMskFlg) == (pTarg[-(int)(Psn+1)] &  ~cRptMskFlg) &&
						(pProbe[-(int)(InDelPsn+2)] &  ~cRptMskFlg) == (pTarg[-(int)(Psn+2)] &  ~cRptMskFlg))
						break;

				if(InDelPsn <= LastInDelPsn)
					{
					NewScore = 3 + Score + QuickAlignRight(AllowedMismatches, 
							MaxMismatchSeqLen,
							AllowedProbeInDels,
							AllowedTargInDels,
							MaxInDelExtension,
							ProbeLen - (InDelPsn + 3),
							&pProbe[-(int)(InDelPsn + 3)],
							TargLen - (Psn + 3), 
							&pTarg[-(int)(Psn + 3)]);

					if(NewScore > IndelScore)
						IndelScore = NewScore;
					}
				}

			if(AllowedTargInDels)
				{
				AllowedTargInDels--;
				
				// lookahead and see if following an InDel event an alignment of at least 3nt can be started...
				InDelPsn = Psn + 1;
				LastInDelPsn = InDelPsn + MaxInDelExtension;
				if(LastInDelPsn > TargLen - 7)
					LastInDelPsn = TargLen - 7;

				for(;InDelPsn <= LastInDelPsn; InDelPsn++)
					if((pProbe[-(int)Psn] &  ~cRptMskFlg) == (pTarg[-(int)InDelPsn] &  ~cRptMskFlg) &&
						(pProbe[-(int)(Psn+1)] &  ~cRptMskFlg) == (pTarg[-(int)(InDelPsn+1)] &  ~cRptMskFlg) &&
						(pProbe[-(int)(Psn+2)] &  ~cRptMskFlg) == (pTarg[-(int)(InDelPsn+2)] &  ~cRptMskFlg))
						break;

				if(InDelPsn <= LastInDelPsn)
					{
					NewScore = 3 + Score + QuickAlignRight(AllowedMismatches, 
							MaxMismatchSeqLen,
							AllowedProbeInDels,
							AllowedTargInDels,
							MaxInDelExtension,
							ProbeLen - (Psn + 3),
							&pProbe[-(int)(Psn + 3)],
							TargLen - (InDelPsn + 3), 
							&pTarg[-(int)(InDelPsn + 3)]);

					if(NewScore > IndelScore)
						IndelScore = NewScore;
					}
				}
			}
		}
	}
return(Score > IndelScore ? Score : IndelScore);
}


// MakeXMLsafe
// Replaces any chars which may cause XML to baff with a safe character '.'
char *
CBioSeqFile::MakeXMLsafe(char *pszStr)
{
char *pszTrans = pszStr;
char Chr;

while(Chr=*pszStr++)
	{
	switch(Chr) {
		case '&':
		case '<':
		case '>':
			pszStr[-1] = '.';
		default:
			continue;
		}
	}
return(pszTrans);
}

