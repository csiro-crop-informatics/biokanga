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

CDataPoints::CDataPoints(void)
{
m_hFile = -1;
m_pDirEls = NULL;
m_pCachedDataset = NULL;
m_AllocdCacheSize = 0;
m_bOpenRead = false;
m_bOpenWrite = false;
memset(m_ChromHashes,0,sizeof(m_ChromHashes));
}

CDataPoints::~CDataPoints(void)
{
if(m_hFile != -1)
	close(m_hFile);
if(m_pDirEls)
	delete m_pDirEls;
if(m_pCachedDataset)
	delete (char *)m_pCachedDataset;
}


void
CDataPoints::InitHdr(void)
{
memset(&m_FileHdr,0,sizeof(m_FileHdr));
m_FileHdr.Magic[0] = 'b';
m_FileHdr.Magic[1] = 'i';
m_FileHdr.Magic[2] = 'o';
m_FileHdr.Magic[3] = 's';
m_FileHdr.Type =	cBSFTypeDataPts;			// it's a dataset file
m_FileHdr.Version =	cDSFVersion;				// file structure version
m_FileHdr.SizeOfHdr =sizeof(tsDSFileHdr);		// size of this header
m_FileHdr.FileLen =	sizeof(tsDSFileHdr);		// current file length
m_FileHdr.MaxDatasets=cDSFMaxDatasets;			// max number of datasets supported
m_FileHdr.MaxChroms = cDSFMaxTotalChroms;		// max number of chromosomes supported
m_FileHdr.MaxDirEls = cDSFMaxDirEls;			// max number of datasegs
m_FileHdr.MaxDataPts = cDSFMaxDataPts;
m_FileHdr.RefDatsetID = 0;
m_FileHdr.DirElOfs = 0;
m_bOpenRead = false;
m_bOpenWrite = false;
m_LoadedDataChromID = 0;
m_LoadedDataDatasetID = 0;
m_LoadedDataChromOfs = 0; 
m_LoadedDataNumPts = 0;
m_FileHdr.szDescription[0] = '\0';
m_FileHdr.szTitle[0] = '\0';
}

int
CDataPoints::Close(bool bWrtDirHdr)
{
int WrtLen;
if(m_hFile != -1)
	{
	if(bWrtDirHdr && m_bOpenWrite)
		{
		if(m_pDirEls != NULL && m_FileHdr.NumDirEls)
			{
			qsort(m_pDirEls,m_FileHdr.NumDirEls,sizeof(sDSFDirEl),CompareDirEls);

			WrtLen = m_FileHdr.NumDirEls * sizeof(sDSFDirEl);	
			if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
				write(m_hFile,m_pDirEls,WrtLen)!=WrtLen)
				{
				AddErrMsg("CDataPoints::Close","Unable to write data segments directory to disk on file %s - error %s",m_szFile,strerror(errno));
				close(m_hFile);
				m_hFile = -1;
				return(eBSFerrFileAccess);
				}
			m_FileHdr.DirElOfs = m_FileHdr.FileLen;
			m_FileHdr.FileLen += WrtLen;
			}
		else
			m_FileHdr.NumDirEls = 0;

		if(_lseeki64(m_hFile,0,SEEK_SET) ||
				write(m_hFile,&m_FileHdr,sizeof(tsDSFileHdr))!=sizeof(tsDSFileHdr))
			{
			AddErrMsg("CDataPoints::Close","Unable to write file header to disk on file %s - error %s",m_szFile,strerror(errno));
			close(m_hFile);
			m_hFile = -1;
			return(eBSFerrFileAccess);
			}
		}
		close(m_hFile);
	m_hFile = -1;

	}
m_szFile[0] = '\0';
if(m_pDirEls)
	{
	delete m_pDirEls;
	m_pDirEls = NULL;
	}
if(m_pCachedDataset)
	{
	delete (char *)m_pCachedDataset;
	m_pCachedDataset = NULL;
	}
m_AllocdCacheSize = 0;
InitHdr();
m_bOpenRead = false;
m_bOpenWrite = false;
return(eBSFSuccess);
}



int
CDataPoints::Open(char *pszFile,	// specifies file to open or create
			   bool bCreate)	// create file or truncate if it already exists
{
int AllocLen;
INT64 FileOfs;
int Len;

Close();						// reset context in case file still opened

// clear any errors that may have been previously generated
// otherwise the user may get confused as which errors are significant!
ClearErrs();

if(pszFile == NULL || *pszFile == '\0') // validate parameters
	{
	AddErrMsg("CDataPoints::Open","No file to open specified");
	return(eBSFerrParams);
	}

if(!bCreate)
#ifdef _WIN32
	m_hFile = open(pszFile, _O_RDWR | _O_BINARY ); // file access is normally random..
else
	m_hFile = open(pszFile,O_CREATETRUNC);
#else
	m_hFile = open64(pszFile,O_RDWR); // file access is normally random..
else
	{
	if((m_hFile = open64(pszFile,O_RDWR | O_CREAT, S_IREAD | S_IWRITE))!=-1)
	   if(ftruncate(m_hFile,0)!=0)
			{
			AddErrMsg("Open","Unable to truncate %s - %s",pszFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
	}
#endif

if(m_hFile == -1)					// check if file open succeeded
	{
	AddErrMsg("CDataPoints::Open","Unable to %s file %s - error %s",bCreate ? "create" : "open",pszFile,strerror(errno));
	return(eBSFerrOpnFile);
	}

strncpy(m_szFile,pszFile,_MAX_PATH);
m_szFile[_MAX_PATH-1] = '\0';

if(bCreate)
	{
	// allocate memory to hold segment directory
	AllocLen = (1 + m_FileHdr.MaxDirEls) * sizeof(sDSFDirEl);	// additional is for safety
	if((m_pDirEls = (sDSFDirEl *)new unsigned char [AllocLen])==NULL)
		{
		AddErrMsg("CDataPoints::Open","Unable to alloc memory (%d byte) for data segments directory from file %s ",AllocLen,pszFile);
		Close();
		return(eBSFerrMem);
		}
	m_bOpenWrite = true;
	}
else // else opening existing file
	{
	_lseeki64(m_hFile,0,SEEK_SET);			// read in header..
	if(sizeof(tsDSFileHdr) != read(m_hFile,&m_FileHdr,sizeof(tsDSFileHdr)))
		{
		AddErrMsg("CDataPoints::Open","Unable to read header from opened file  %s - error %s",pszFile,strerror(errno));
		Close();			// closes opened file..
		return(eBSFerrFileAccess);
		}

	// header read, validate it as being a bioseq data point file header
	if(tolower(m_FileHdr.Magic[0]) != 'b' || 
		tolower(m_FileHdr.Magic[1]) != 'i' || 
		tolower(m_FileHdr.Magic[2]) != 'o' || 
		tolower(m_FileHdr.Magic[3]) != 's')
		{
		AddErrMsg("CDataPoints::Open","Opened file %s is not a bioseq file",pszFile);
		Close();			// closes opened file..
		return(eBSFerrNotBioseq);
		}

	// is it a dataset file?
	if(m_FileHdr.Type != cBSFTypeDataPts)
		{
		AddErrMsg("CDataPoints::Open","Opened file %s is a bioseq file but does not contain data points",pszFile);
		Close();			// closes opened file..
		return(eBSFerrFileDPT);
		}

	// can we handle this version?
	if(m_FileHdr.Version < cDSFVersionBack || m_FileHdr.Version > cDSFVersion)
		{
		AddErrMsg("CDataPoints::Open","Incompatible file version %d in opened file %s, expected between version %d and %d",
			m_FileHdr.Version,pszFile,cDSFVersionBack,cDSFVersion);
		Close();			// closes opened file..
		return(eBSFerrFileVer);
		}

		// if version was version 2 or later, and no description/title supplied then generate the description/title from the file name...
	if(m_FileHdr.Version >= 2)
		{
		if(m_FileHdr.szDescription[0] == '\0')
			{
			strncpy(m_FileHdr.szDescription,pszFile,sizeof(m_FileHdr.szDescription));
			m_FileHdr.szDescription[sizeof(m_FileHdr.szDescription)-1] = '\0';
			}

		if(m_FileHdr.szTitle[0] == '\0')
			{
			char szFname[_MAX_FNAME];
#ifdef _WIN32
			_splitpath(pszFile,NULL,NULL,szFname,NULL);
#else
			CUtility::splitpath(pszFile,NULL,szFname);
#endif
			strncpy(m_FileHdr.szTitle,szFname,sizeof(m_FileHdr.szTitle));
			m_FileHdr.szTitle[sizeof(m_FileHdr.szTitle)-1] = '\0';
			}
		}

	// allocate memory to hold segment directory
	AllocLen = m_FileHdr.NumDirEls * sizeof(sDSFDirEl);	
	if((m_pDirEls = (sDSFDirEl *)new unsigned char [AllocLen])==NULL)
		{
		AddErrMsg("CDataPoints::Open","Unable to alloc memory (%d byte) for data segments directory from opened file %s ",AllocLen,pszFile);
		Close();
		return(eBSFerrMem);
		}
	if((FileOfs=_lseeki64(m_hFile,m_FileHdr.DirElOfs,SEEK_SET)) != m_FileHdr.DirElOfs)
		{
		AddErrMsg("CDataPoints::Open","Unable to seek to %I64d prior to reading data segments directory on opened file %s - error %s",m_FileHdr.DirElOfs,pszFile,strerror(errno));
		Close();
		return(eBSFerrFileAccess);
		}

	if((Len = read(m_hFile,m_pDirEls,AllocLen)) != AllocLen)
		{
		AddErrMsg("CDataPoints::Open","Unable to read data segments directory from opened file %s - error %s",pszFile,strerror(errno));
		Close();
		return(eBSFerrFileAccess);
		}
	m_bOpenRead = true;
	}

return(eBSFSuccess);
}



teBSFrsltCodes 
CDataPoints::SetDescription(char *pszDescription)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(!m_bOpenWrite)
	return(eBSFerrRead);
if(m_FileHdr.Version >= 2)
	{
	strncpy(m_FileHdr.szDescription,pszDescription,sizeof(m_FileHdr.szDescription));
	m_FileHdr.szDescription[sizeof(m_FileHdr.szDescription)-1] = '\0';
	return(eBSFSuccess);
	}
return(eBSFerrFileVer);
}

teBSFrsltCodes 
CDataPoints::GetDescription(int MaxLen,char *pszDescription)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(m_bOpenWrite)
	return(eBSFerrWrite);
if(m_FileHdr.Version >= 2)
	strncpy(pszDescription,m_FileHdr.szDescription,MaxLen);
else
	strncpy(pszDescription,m_szFile,MaxLen);
pszDescription[MaxLen-1] = '\0';
return(eBSFSuccess);
}

teBSFrsltCodes 
CDataPoints::SetTitle(char *pszTitle)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(!m_bOpenWrite)
	return(eBSFerrRead);
if(m_FileHdr.Version >= 2)
	{
	strncpy(m_FileHdr.szTitle,pszTitle,sizeof(m_FileHdr.szTitle));
	m_FileHdr.szTitle[sizeof(m_FileHdr.szTitle)-1] = '\0';
	return(eBSFSuccess);
	}
return(eBSFerrFileVer);
}

teBSFrsltCodes 
CDataPoints::GetTitle(int MaxLen,char *pszTitle)
{
if(m_hFile == -1)
	return(eBSFerrClosed);
if(m_bOpenWrite)
	return(eBSFerrWrite);
if(m_FileHdr.Version >= 2)
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

etDataPtType 
CDataPoints::GetDataType(void)
{
if(m_hFile == -1)
	return((etDataPtType)eBSFerrClosed);
return(m_FileHdr.DataPtType);
}

// GenNameHash
// Generates a 16bit hash on specified lowercased name
// This hash can then be used to quickly eliminate probe names which can't match a target name by comparing hashes
unsigned short 
CDataPoints::GenNameHash(char *pszName)
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

tsFDSRelChrom *
CDataPoints::GetRelChrom(int RelChromID)
{
if(m_hFile == -1)
	return(NULL);
if(RelChromID < 1 || RelChromID > m_FileHdr.NumChroms)
	return(NULL);
return(&m_FileHdr.Chroms[RelChromID-1]);
}



tsFDSRefChrom *
CDataPoints::GetRefChrom(int RefChromID)
{
tsFDSRefChrom *pRefChrom;
int Idx;
if(m_hFile == -1)
	return(NULL);
if(RefChromID < 1 || RefChromID > m_FileHdr.NumChroms)
	return(NULL);
if(m_FileHdr.Chroms[RefChromID-1].DatasetID != m_FileHdr.RefDatsetID)
	return(NULL);
pRefChrom = &m_FileHdr.RefChroms[0];
for(Idx = 0; Idx < m_FileHdr.NumRefChroms; Idx++, pRefChrom++)
	if(pRefChrom->ChromID == RefChromID)
		return(pRefChrom);
return(NULL);
}

// SetRefDataset
// Sets the dataset designated as being the reference dataset
// Must be the only reference dataset and the first dataset created
int 
CDataPoints::SetRefDataset(char *pszRefDataset,
						    etDataPtType DataPtType) // data point type of all data points in this file
{
tsFDSDataset *pDataset;

if(m_FileHdr.RefDatsetID || m_FileHdr.NumDatasets)	
	return(eBSFerrRefDataset);
switch(DataPtType) {
	case eDPTPackedBases:					// packed bases - 2 per byte
	case eDPTbyte:							// 8 bit
		m_FileHdr.DataPtSize = 1;
		break;
	case eDPTword:							// 16bit
		m_FileHdr.DataPtSize = 2;
		break;
	case eDPTdword:							// 32 bit
		m_FileHdr.DataPtSize = 4;
		break;
	case eDPTdouble:						// 8 byte
		m_FileHdr.DataPtSize = 8;
		break;
	default:
		return(eBSFerrDataPtType);
	}

pDataset = &m_FileHdr.Datasets[0];
memset(pDataset,0,sizeof(tsFDSDataset));
m_FileHdr.RefDatsetID = pDataset->DatasetID = ++m_FileHdr.NumDatasets;
m_FileHdr.DataPtType  = DataPtType;
strncpy(pDataset->szName,pszRefDataset,cMaxDatasetSpeciesChrom);
pDataset->szName[cMaxDatasetSpeciesChrom-1] = '\0';
pDataset->Hash = GenNameHash(pszRefDataset);
pDataset->NumSegs = 0;
return(m_FileHdr.RefDatsetID);
}

int
CDataPoints::AddRefChrom(char *pszRefChrom,int ChromLen)
{
tsFDSRefChrom *pRefChrom;
tsFDSRelChrom *pNewRefChrom;
int RefChromID;

RefChromID = GetChromID(m_FileHdr.RefDatsetID,pszRefChrom);
if(RefChromID > eBSFSuccess)	// already added?
	return(RefChromID);
if(m_FileHdr.NumChroms >= m_FileHdr.MaxChroms)
	{
	AddErrMsg("CDataPoints::AddRefChrom","Too many chromosomes, can handle at most %d for file %s ",m_FileHdr.MaxChroms,m_szFile);
	return(eBSFerrMaxChroms);
	}
if(m_FileHdr.NumRefChroms >= cDSFMaxChromsDataset)
	{
	AddErrMsg("CDataPoints::AddRefChrom","Too many reference chromosomes, can handle at most %d for file %s ",cDSFMaxChromsDataset,m_szFile);
	return(eBSFerrMaxChroms);
	}

pNewRefChrom = &m_FileHdr.Chroms[m_FileHdr.NumChroms++];
memset(pNewRefChrom,0,sizeof(tsFDSRelChrom));
RefChromID = pNewRefChrom->ChromID = m_FileHdr.NumChroms;
strncpy(pNewRefChrom->szName,pszRefChrom,cMaxDatasetSpeciesChrom);
pNewRefChrom->szName[cMaxDatasetSpeciesChrom-1] = '\0';
pNewRefChrom->Hash = GenNameHash(pszRefChrom);
pNewRefChrom->DatasetID = m_FileHdr.RefDatsetID;
pNewRefChrom->ChromLen = ChromLen;

if(m_ChromHashes[pNewRefChrom->Hash] != 0)
	pNewRefChrom->NxtChromID = m_ChromHashes[pNewRefChrom->Hash];
m_ChromHashes[pNewRefChrom->Hash] = RefChromID;

pRefChrom = &m_FileHdr.RefChroms[m_FileHdr.NumRefChroms++];
memset(pRefChrom,0,sizeof(tsFDSRefChrom));
pRefChrom->ChromID = RefChromID;
pRefChrom->ChromOfs = -1;
pRefChrom->NumPts = -1;
m_FileHdr.Datasets[0].NumChroms = m_FileHdr.NumRefChroms;
return(RefChromID);
}

int
CDataPoints::AddRelChrom(char *pszRelChrom,int RelDatasetID,int ChromLen)
{
int RelChromID;
tsFDSRelChrom *pRelChrom;
tsFDSDataset *pDataset;

if(RelDatasetID == m_FileHdr.RefDatsetID)
	return(AddRefChrom(pszRelChrom,ChromLen));

RelChromID = GetChromID(RelDatasetID,pszRelChrom);
if(RelChromID > eBSFSuccess)	// already known?
	return(RelChromID);

if(m_FileHdr.NumChroms >= m_FileHdr.MaxChroms)
	{
	AddErrMsg("CDataPoints::AddRelChrom","Too many chromosomes, can handle at most %d for file %s ",m_FileHdr.MaxChroms,m_szFile);
	return(eBSFerrMaxChroms);
	}
pDataset = &m_FileHdr.Datasets[RelDatasetID-1];
if(pDataset->NumChroms >= cDSFMaxChromsDataset)
	{
	AddErrMsg("CDataPoints::AddRelChrom","Too many chromosomes for dataset '%s', can handle at most %d",pDataset->szName,cDSFMaxChromsDataset);
	return(eBSFerrMaxChroms);
	}
pDataset->NumChroms += 1;
pRelChrom = &m_FileHdr.Chroms[m_FileHdr.NumChroms++];
memset(pRelChrom,0,sizeof(tsFDSRelChrom));
RelChromID = pRelChrom->ChromID = m_FileHdr.NumChroms;
strncpy(pRelChrom->szName,pszRelChrom,cMaxDatasetSpeciesChrom);
pRelChrom->szName[cMaxDatasetSpeciesChrom-1] = '\0';
pRelChrom->Hash = GenNameHash(pszRelChrom);
pRelChrom->ChromLen = ChromLen;
pRelChrom->DatasetID = RelDatasetID;
if(m_ChromHashes[pRelChrom->Hash] != 0)
	pRelChrom->NxtChromID = m_ChromHashes[pRelChrom->Hash];
m_ChromHashes[pRelChrom->Hash] = RelChromID;
return(RelChromID);
}


int
CDataPoints::AddRelDataset(char *pszRelDataset)
{
tsFDSDataset *pDataset;
int RelDatasetID = GetDatasetID(pszRelDataset);
if(RelDatasetID < eBSFSuccess)		// if a new dataset being started - can't be reference as that was set..
	{
	if(m_FileHdr.NumDatasets >= cDSFMaxDatasets)
		{
		AddErrMsg("CDataPoints::AddRelDataset","Too many datasets (species), can handle at most %d for file %s ",cDSFMaxDatasets,m_szFile);
		return(eBSFerrMaxDatasets);
		}
	pDataset = &m_FileHdr.Datasets[m_FileHdr.NumDatasets++];
	memset(pDataset,0,sizeof(tsFDSDataset));
	RelDatasetID = pDataset->DatasetID = m_FileHdr.NumDatasets;
	strncpy(pDataset->szName,pszRelDataset,cMaxDatasetSpeciesChrom);
	pDataset->szName[cMaxDatasetSpeciesChrom-1] = '\0';
	pDataset->Hash = GenNameHash(pszRelDataset);
	pDataset->NumSegs = 0;
	}
return(RelDatasetID);
}

int
CDataPoints::AddDataseg(char *pszRefChrom,		// which reference dataset chromosome
					 int RefChromOfs,			// start on reference chromosome (0..n)
					 char *pszRelDataset,		// which relative dataset
					 char *pszRelChrom,			// which relative dataset chromosome
					 char RelStrand,			// on which relative strand '+' or '-'
					 int RelChromOfs,			// start on relative chromosome 
					 int NumPts,				// number of datapoints
					 void *pDataPts)
{
int RelDatasetID;
int RelChromID;
int RefChromID;
if(pszRelDataset == NULL || pszRelDataset[0] == '\0' ||
   pszRefChrom == NULL || pszRefChrom[0] == '\0' ||
   pszRelChrom == NULL || pszRelChrom[0] == '\0' ||
   RelChromOfs < 0 || RefChromOfs < 0) 
	return(eBSFerrParams);

if(NumPts <= 0 || pDataPts == NULL)
	return(eBSFerrParams);

if(!m_FileHdr.RefDatsetID)
	return(eBSFerrRefDataset);

RelDatasetID = GetDatasetID(pszRelDataset);
if(RelDatasetID < eBSFSuccess)
	return(RelDatasetID);
RelChromID   = GetChromID(RelDatasetID,pszRelChrom);
if(RelChromID < eBSFSuccess)
	return(RelChromID);
RefChromID   = GetChromID(m_FileHdr.RefDatsetID,pszRefChrom);
if(RefChromID < eBSFSuccess)
	return(RefChromID);

return(AddDataseg(RefChromID,RefChromOfs,RelDatasetID,RelChromID,RelStrand,RelChromOfs,NumPts,pDataPts));
}


// AddDataseg
// Add new data segment
// Note that reference dataset chromosome and relative dataset chromome must have been previously added
int
CDataPoints::AddDataseg(int RefChromID,		// which reference dataset chromosome
					 int RefChromOfs,		// start on reference chromosome (0..n)
					 int RelDatasetID,		// which relative dataset
					 int RelChromID,		// which relative dataset chromosome
					 char RelStrand,		// on which relative strand '+' or '-'
					 int RelChromOfs,		// start on relative chromosome 
					 int NumPts,			// number of datapoints
					 void *pDataPts)		// array of data points
{
tsFDSRelChrom *pRelChrom;
tsFDSRefChrom *pRefChrom;
tsFDSDataset *pDataset;

sFDSChromDataset *pChromDataset;
sDSFDirEl *pDirEl;
tsOctStructParam *pStructConf = (tsOctStructParam *)pDataPts;

int DataSegLen;
int PackedBytes;
int Num2Alloc;
void *pData2Save;

if(!m_bOpenWrite)
	return(eBSFerrWrite);

if(RelChromOfs < 0 || RefChromOfs < 0) 
	return(eBSFerrParams);

if(NumPts <= 0 || pDataPts == NULL)
	return(eBSFerrParams);

if(!m_FileHdr.RefDatsetID)
	return(eBSFerrRefDataset);

if(m_FileHdr.NumDirEls == cDSFMaxDirEls)
	{
	AddErrMsg("CDataPoints::AddDataseg","Too many data segments, can handle at most %d for file %s ",cDSFMaxDirEls,m_szFile);
	return(eBSFerrMaxDirEls);
	}
if(m_FileHdr.DataPtType == eDPTPackedBases)
	{
	PackedBytes = (NumPts + 1)/2;
	if(m_pCachedDataset == NULL || PackedBytes > m_AllocdCacheSize)
		{
		if(m_pCachedDataset != NULL)
			delete (char *)m_pCachedDataset;
		Num2Alloc = (PackedBytes * 3)/2;
		if(Num2Alloc < 100000)
			Num2Alloc = 100000;
		if(Num2Alloc > 1000000)
			Num2Alloc = PackedBytes;
		if((m_pCachedDataset = new unsigned char [Num2Alloc])==NULL)
			{
			AddErrMsg("CDataPoints::AddDataseg","Unable to alloc memory (%d byte) for packed base segment cache in file %s ",Num2Alloc,m_szFile);
			return(eBSFerrMem);
			}
		m_AllocdCacheSize = Num2Alloc;
		}

	PackBases(NumPts,				// unpacked sequence length 
		(etSeqBase *)pDataPts,		// pts to sequence to pack from
		0,							// nibble to start pack into (0..SeqLen-1)
		(unsigned char *)m_pCachedDataset);	// where to pack into

	DataSegLen = PackedBytes;
	pData2Save = m_pCachedDataset;
	}
else
	{
	DataSegLen = NumPts * m_FileHdr.DataPtSize;
	pData2Save = pDataPts;
	}

pDataset = &m_FileHdr.Datasets[RelDatasetID-1];
pRelChrom = &m_FileHdr.Chroms[RelChromID-1];

pRefChrom = GetRefChrom(RefChromID);
if(pRefChrom->NumPts == -1)		// first time for this chromosome?
	{
	pRefChrom->ChromOfs = RefChromOfs;
	pRefChrom->NumPts = NumPts;
	}
pChromDataset = &pRefChrom->Datasets[RelDatasetID-1];
if(!pChromDataset->NumDatasegs)	// 1st time for dataset on this chromosome?
	{
	pChromDataset->ChromOfs = INT_MAX;
	pRefChrom->NumDatasets++;
	}
pChromDataset->NumDatasegs++;
pRefChrom->NumDatasegs++;
if(RelChromOfs < pChromDataset->ChromOfs)
	pChromDataset->ChromOfs = RelChromOfs;
if(RelChromOfs + NumPts > pChromDataset->ChromOfs + pChromDataset->NumPts)
	pChromDataset->NumPts = RelChromOfs + NumPts - pChromDataset->ChromOfs;

if(RefChromOfs < pRefChrom->ChromOfs)
	pRefChrom->ChromOfs = RelChromOfs;
if(RefChromOfs + NumPts > pRefChrom->ChromOfs + pRefChrom->NumPts)
	pRefChrom->NumPts = RefChromOfs + NumPts - pRefChrom->ChromOfs;

pDirEl = &m_pDirEls[m_FileHdr.NumDirEls++];
pDirEl->RefChromID = RefChromID;
pDirEl->RefChromOfs = RefChromOfs;
pDirEl->RelChromID = RelChromID;
pDirEl->RelChromOfs = RelChromOfs;
pDirEl->RelStrand = RelStrand == '+' ? 0 : 1;
pDirEl->RelDatasetID = RelDatasetID;
pDirEl->NumPts = NumPts;
pDirEl->LoFileOfs= (UINT32)(m_FileHdr.FileLen & 0x0ffffffff);
pDirEl->HiFileOfs= (UINT8)((m_FileHdr.FileLen >> 32L) & 0x0ff);

if(DataSegLen > (int)m_FileHdr.DataSegSize)
	m_FileHdr.DataSegSize = DataSegLen;

if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||			
	write(m_hFile,pData2Save,DataSegLen)!=DataSegLen)
		{
		AddErrMsg("CDataPoints::AddDataseg","Unable to write data segment to disk on file %s - error %s",m_szFile,strerror(errno));
		Close();			// closes opened file..
		return(eBSFerrFileAccess);
		}
m_FileHdr.FileLen += DataSegLen;
pDataset->NumSegs++;	
return(eBSFSuccess);
}

// LocateRelDirEl
// Locates and returns a ptr to matching DirEl where the RelDatasetID and RefChromID exactly match, and
// the specified RefChromOfs is on the matching data point segment.
// If bClosest is true then if no exact match on ChromID.DatasetID.ChromOfs will reurn
// a ptr to the DirEl which has a ChromOfs closest but after than the ChromOfs specified.
sDSFDirEl *
CDataPoints::LocateRelDirEl(int RelDatasetID,int RefChromID,int RefChromOfs,bool bClosest)
{
sDSFDirEl *pProbe;
sDSFDirEl *pClosest = NULL;
int Left = 0;
int Right = m_FileHdr.NumDirEls - 1;
int MidPt;

while(Right >= Left) {
	MidPt = (Right + Left)/2;
	pProbe = &m_pDirEls[MidPt];

	if(pProbe->RelDatasetID > RelDatasetID)
		{
		Right = MidPt - 1;
		continue;
		}
	else
		if(pProbe->RelDatasetID < RelDatasetID)
			{
			Left = MidPt + 1;
			continue;
			}

	if(pProbe->RefChromID > RefChromID)
		{
		Right = MidPt - 1;
		continue;
		}
	else
		if(pProbe->RefChromID < RefChromID)
			{
			Left = MidPt + 1;
			continue;
			}

	if(pProbe->RefChromOfs > RefChromOfs)
		{
		pClosest = pProbe;
		Right = MidPt - 1;
		continue;
		}
	else
		if((pProbe->RefChromOfs + pProbe->NumPts) <= RefChromOfs)
			{
			Left = MidPt + 1;
			continue;
			}

	return(pProbe);
	}
return(bClosest ? pClosest : NULL);
}

int
CDataPoints::FillDataPts(void *pDataPts,void *pFiller,int DataPtSize,int NumDataPts)
{
double *pDouble;
UINT32 *pInt;
unsigned short *pShort;
unsigned char *pByte;
if(pDataPts == NULL || !DataPtSize || !NumDataPts)
	return(eBSFerrParams);
if(pFiller == NULL) // nothing to do - treat as success
	return(eBSFSuccess);

switch(DataPtSize) {
	case 1:				// do a memset
		memset((unsigned char *)pDataPts,*(unsigned char *)pFiller,NumDataPts);
		break;

	case 2:				 
		pShort = (unsigned short *)pDataPts;
		while(NumDataPts--)
			*pShort++ = *(unsigned short *)pFiller;
		break;

	case 4:
		pInt = (UINT32 *)pDataPts;
		while(NumDataPts--)
			*pInt++ = *(UINT32 *)pFiller;
		break;

	case 8:
		pDouble = (double *)pDataPts;
		while(NumDataPts--)
			*pDouble++ = *(double *)pFiller;
		break;

	default:
		pByte = (unsigned char *)pDataPts;
		while(NumDataPts--)
			{
			memmove(pByte,pFiller,DataPtSize);
			pByte += DataPtSize;
			}
		break;
	}
return(eBSFSuccess);
}

int 
CDataPoints::GetNumDatasets(void)
{
if(m_hFile != -1) 
	return(m_FileHdr.NumDatasets);
return(eBSFerrClosed);
}

// GetNumChroms
// Returns the total number of chromosomes over all datasets
int 
CDataPoints::GetNumChroms(void)
{
if(m_hFile != -1) 
	return(m_FileHdr.NumChroms);
return(eBSFerrClosed);
}

// GetNumChroms
// Returns the number of chromosomes in the reference dataset
int 
CDataPoints::GetNumRefChroms(void)
{
if(m_hFile != -1) 
	return(m_FileHdr.NumRefChroms);
return(eBSFerrClosed);
}

// GetNxtRefChromID
// Returns the next reference dataset chromosome identifier which is immediately higher in value than that specified
// Returns 0 when all reference chromIDs have been returned
int
CDataPoints::GetNxtRefChromID(int CurID)	// intialise with 0 to get first chromID
{
int Idx;
int MaxCurID = INT_MAX;
tsFDSRefChrom *pRefChrom;
if(m_hFile == -1) 
	return(eBSFerrClosed);

pRefChrom = m_FileHdr.RefChroms;
for(Idx = 0; Idx < m_FileHdr.NumRefChroms; Idx++,pRefChrom++)
	{
	if(pRefChrom->ChromID > CurID &&
		pRefChrom->ChromID < MaxCurID)
		MaxCurID = pRefChrom->ChromID;
	}
return(MaxCurID == INT_MAX ? 0 : MaxCurID);
}

char *
CDataPoints::GetChromName(int ChromID)
{
if(m_hFile != -1) 
	if(ChromID > 0 && ChromID <= m_FileHdr.NumChroms)
		return(m_FileHdr.Chroms[ChromID-1].szName);
return(NULL);
}

char *
CDataPoints::GetDatasetName(int DatasetID)
{
if(m_hFile != -1)
	if(DatasetID > 0 && DatasetID <= m_FileHdr.NumDatasets)
		return(m_FileHdr.Datasets[DatasetID-1].szName);
return(NULL);
}

// GetDatasetID
// Returns DatasetID (1..n) which is used to identify the dataset from it's name
// eBSFerrDataset is returned if unable to locate dataset
int
CDataPoints::GetDatasetID(char *pszName)
{
tsFDSDataset *pDataset;
int Idx;
unsigned short Hash;
if(m_hFile == -1)
	return(eBSFerrClosed);

Hash = GenNameHash(pszName);
pDataset = &m_FileHdr.Datasets[0];
for(Idx = 0; Idx < m_FileHdr.NumDatasets; Idx++, pDataset++)
	{
	if(pDataset->Hash == Hash &&
		!stricmp(pDataset->szName,pszName))
		return(pDataset->DatasetID);
	}
return(eBSFerrDataset);
}

// GetChromID
// Returns ChromID (1..n) which is used to identify the chromosome from it's name
// eBSFerrChrom is returned if unable to locate chromosome
int
CDataPoints::GetChromID(int DatasetID,char *pszName)
{
tsFDSRelChrom *pChrom;
int Idx;
unsigned short Hash;
if(m_hFile == -1)
	return(eBSFerrClosed);

if(!m_FileHdr.NumChroms)
	return(eBSFerrChrom);

Hash = GenNameHash(pszName);
if(!(Idx = m_ChromHashes[Hash]))
	return(eBSFerrChrom);
do {
	pChrom = &m_FileHdr.Chroms[Idx-1];
	if(pChrom->DatasetID == DatasetID &&
		pChrom->Hash == Hash &&
		!stricmp(pChrom->szName,pszName))
			return(pChrom->ChromID);
	}
while((Idx = pChrom->NxtChromID) > 0);
return(eBSFerrChrom);
}

sFDSChromDataset *
CDataPoints::GetRefChromDataset(int RelDatasetID,int RefChromID)
{
tsFDSRefChrom *pRefChrom;
int Idx;
if(RelDatasetID < 1 || RelDatasetID > m_FileHdr.NumDatasets)
	return(NULL);
pRefChrom = m_FileHdr.RefChroms;
for(Idx = 0 ; Idx < m_FileHdr.NumRefChroms; Idx++, pRefChrom++)
	if(pRefChrom->ChromID == RefChromID)
		{

		return(&pRefChrom->Datasets[RelDatasetID-1]);
		}
return(NULL);
}

// GetRefNumDatasets
// Returns number of relative datasets which have at least 1 dataseg on reference chromosome
int CDataPoints::GetRefNumDatasets(int RefChromID)
{
tsFDSRefChrom *pRefChrom;
if((pRefChrom =  GetRefChrom(RefChromID))!=NULL)
	return(pRefChrom->NumDatasets);
return(eBSFerrChrom);
}

// GetRefNumSegs
// Returns number of segments on RefChromID for specified relative dataset 
int CDataPoints::GetRefNumSegs(int RelDatasetID,int RefChromID)
{
sFDSChromDataset *pChromDataset;
if((pChromDataset = GetRefChromDataset(RelDatasetID,RefChromID))!=NULL)
	return(pChromDataset->NumDatasegs);
return(eBSFerrDataset);
}

// GetRefChromOfs
// Returns offset of first datapoint on ChromID for specified relative dataset 
int CDataPoints::GetRefChromOfs(int RelDatasetID,int RefChromID)
{
sFDSChromDataset *pChromDataset;
if((pChromDataset = GetRefChromDataset(RelDatasetID,RefChromID))!=NULL)
	return(pChromDataset->ChromOfs);
return(eBSFerrDataset);
}

// GetRefNumPts
// Returns number of data points (including missing) from first datapoint on RefChromID to last
// for specified relative dataset 
int CDataPoints::GetRefNumPts(int RelDatasetID,int RefChromID)
{
sFDSChromDataset *pChromDataset;
if((pChromDataset = GetRefChromDataset(RelDatasetID,RefChromID))!=NULL)
	return(pChromDataset->NumPts);
return(eBSFerrDataset);
}

// GetFirstOfs
// Returns offset of first datapoint on specified chromosome
int CDataPoints::GetFirstOfs(int RelDatasetID,int RefChromID)
{
sFDSChromDataset *pChromDataset;
if((pChromDataset = GetRefChromDataset(RelDatasetID,RefChromID))!=NULL)
	return(pChromDataset->ChromOfs);
return(eBSFerrDataset);
}

// GetLastOfs
// returns offset of last datapoint on specified chromosome
int CDataPoints::GetLastOfs(int RelDatasetID,int RefChromID)
{
sFDSChromDataset *pChromDataset;
if((pChromDataset = GetRefChromDataset(RelDatasetID,RefChromID))!=NULL)
	return(pChromDataset->ChromOfs + pChromDataset->NumPts - 1);
return(eBSFerrDataset);
}

// GetStrand
// Returns strand for datapoint - '+' or '-'
// If datapoint does not exist then returns '?'
char
CDataPoints::GetStrand(int RelDatasetID,int RefChromID,int RefChromOfs)
{
sDSFDirEl *pDirEl;
pDirEl = LocateRelDirEl(RelDatasetID,RefChromID,RefChromOfs,false);
if(pDirEl == NULL)
	return('?');
return(pDirEl->RelStrand == 0 ? '+' : '-');
}

// GetShuffledValuesBase
// Gets the aligned sequence and then shuffles it whilst preserving dinucleotide frequencies
// Returns shuffled dinucleotide values as etSeqBases
// Note that dinucleotide frequencies are preserved within each aligned segment
int 
CDataPoints::GetShuffledValuesBase(int RelDatasetID,// relative dataset identifier
					  int RefChromID,		// reference chromosome identifier
					  int RefChromOfs,		// start on reference chromosome (0..n) 
					  int NumPts,			// number of points to return
					  unsigned char *pToRet,	// where to return etSeqBases
					  etSeqBase MissingMarker)	// value to use if missing datapoints
{
int Rslt;
int SegStart;
int SegLen;
etSeqBase *pSeg;
int Psn;

Rslt = GetValuesBase(RelDatasetID,RefChromID,RefChromOfs,NumPts,pToRet,MissingMarker);
if(Rslt < eBSFSuccess)
	return(Rslt);

srand((UINT32)RefChromOfs);
SegStart = SegLen = 0;
pSeg = pToRet;
for(Psn =  0; Psn < NumPts; Psn++,pSeg++)
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
return(eBSFSuccess);
}

// ShuffleDinucs
// a) If SeqLen < 30 then simple random swaps are used which preserves base frequencies only
// b) If SeqLen >= 30 then 1st order markov used which preserves dinucleotide frequencies
// 
void
CDataPoints::ShuffleDinucs(etSeqBase *pSeq,int SeqLen)
{
int Psn;
etSeqBase Base;
int TotDinucs = 0;
int Accum;

// if number of dinucleotides < 500 then just do simple random swaps
// as otherwise dinucleotide shuffling has too constraints and generated sequence
// can end up being too close to the orginal
if(SeqLen < 500)	// caution - if upper limit changed to be more than RAND_MAX (32767) then use CShuffle::SeqShuffle()
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

// GetValuesBase
// Gets values as etSeqBases
int 
CDataPoints::GetValuesBase(int RelDatasetID,// relative dataset identifier
					  int RefChromID,		// reference chromosome identifier
					  int RefChromOfs,		// start on reference chromosome (0..n) 
					  int NumPts,			// number of points to return
					  unsigned char *pToRet,	// where to return etSeqBases
					  etSeqBase MissingMarker)	// value to use if missing datapoints
{
int Idx;
char *pChar;
short *pWord;
int *pDWord;
double *pDouble;
int Ofs;
int iRslt;

if(m_hFile == -1)
	return(eBSFerrClosed);
if(RefChromID <= eBSFSuccess || RefChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
if(RelDatasetID <= eBSFSuccess || RelDatasetID > m_FileHdr.NumDatasets)
	return(eBSFerrDataset);

// ensure datapoints are sitting in cache
if((iRslt=PreloadDataset(RelDatasetID,RefChromID,RefChromOfs,NumPts,&MissingMarker)) < eBSFSuccess)
	return(iRslt);

Ofs = RefChromOfs - m_LoadedDataChromOfs;
switch(m_FileHdr.DataPtType) {
	case eDPTPackedBases:		// bases packed 2 per byte
	case eDPTbyte:				// 8 bit byte or char
		pChar = (char *)m_pCachedDataset;
		pChar += Ofs;
		memmove(pToRet,pChar,NumPts);
		break;

	case eDPTword:				// 16 bit
		pWord = (short *)m_pCachedDataset;
		pWord += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pWord++)
			pToRet[Idx] = *(unsigned char *)pWord;
		break;

	case eDPTdword:				// 32 bit
		pDWord = (int *)m_pCachedDataset;
		pDWord += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pDWord++)
			pToRet[Idx] = *(unsigned char *)pDWord;
		break;

	case eDPTdouble:			// double
		pDouble = (double *)m_pCachedDataset;
		pDouble += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pDouble++)
			pToRet[Idx] = *(unsigned char *)pDouble;
		break;

	default:
		return(eBSFerrUserType);
	}
return(eBSFSuccess);
}

int 
CDataPoints::GetValuesChar(int RelDatasetID,		// relative dataset identifier
					  int RefChromID,		// reference chromosome identifier
					  int RefChromOfs,		// start on reference chromosome (0..n) 
					  int NumPts,			// number of points to return
					  char *pToRet,			// where to return chars
					  char MissingMarker)	// value to use if missing datapoints
{
int Idx;
char *pChar;
short *pWord;
int *pDWord;
double *pDouble;
int Ofs;
int iRslt;

if(m_hFile == -1)
	return(eBSFerrClosed);
if(RefChromID <= eBSFSuccess || RefChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
if(RelDatasetID <= eBSFSuccess || RelDatasetID > m_FileHdr.NumDatasets)
	return(eBSFerrDataset);

if((iRslt=PreloadDataset(RelDatasetID,RefChromID,RefChromOfs,NumPts,&MissingMarker)) < eBSFSuccess)
	return(iRslt);

Ofs = RefChromOfs - m_LoadedDataChromOfs;
switch(m_FileHdr.DataPtType) {
	case eDPTPackedBases:		// bases packed 2 per byte
	case eDPTbyte:				// 8 bit byte or char
		pChar = (char *)m_pCachedDataset;
		pChar += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pChar++)
			pToRet[Idx] = *(unsigned char *)pChar;
		break;

	case eDPTword:				// 16 bit
		pWord = (short *)m_pCachedDataset;
		pWord += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pWord++)
			pToRet[Idx] = *(unsigned char *)pWord;
		break;

	case eDPTdword:				// 32 bit
		pDWord = (int *)m_pCachedDataset;
		pDWord += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pDWord++)
			pToRet[Idx] = *(unsigned char *)pDWord;
		break;

	case eDPTdouble:			// double
		pDouble = (double *)m_pCachedDataset;
		pDouble += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pDouble++)
			pToRet[Idx] = (unsigned char)(int)*pDouble;
		break;

	default:
		return(eBSFerrUserType);
	}
return(eBSFSuccess);
}

int 
CDataPoints::GetValuesShort(int RelDatasetID,	// relative dataset identifier
					  int RefChromID,		// reference chromosome identifier
					  int RefChromOfs,		// start on reference chromosome (0..n) 
					  int NumPts,			// number of points to return
					  short int *pToRet,	// where to return shorts
					  short int MissingMarker)	// value to use if missing datapoints
{
int Idx;
char *pChar;
short *pWord;
int *pDWord;
double *pDouble;
int Ofs;
int iRslt;

if(m_hFile == -1)
	return(eBSFerrClosed);
if(RefChromID <= eBSFSuccess || RefChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
if(RelDatasetID <= eBSFSuccess || RelDatasetID > m_FileHdr.NumDatasets)
	return(eBSFerrDataset);

if((iRslt=PreloadDataset(RelDatasetID,RefChromID,RefChromOfs,NumPts,&MissingMarker)) < eBSFSuccess)
	return(iRslt);

Ofs = RefChromOfs - m_LoadedDataChromOfs;
switch(m_FileHdr.DataPtType) {
	case eDPTPackedBases:		// bases packed 2 per byte
	case eDPTbyte:				// 8 bit byte or char
		pChar = (char *)m_pCachedDataset;
		pChar += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pChar++)
			pToRet[Idx] = (short)*pChar;
		break;

	case eDPTword:				// 16 bit
		pWord = (short *)m_pCachedDataset;
		pWord += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pWord++)
			pToRet[Idx] = (short)*pWord;
		break;

	case eDPTdword:				// 32 bit
		pDWord = (int *)m_pCachedDataset;
		pDWord += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pDWord++)
			pToRet[Idx] = (short)*pDWord;
		break;

	case eDPTdouble:			// double
		pDouble = (double *)m_pCachedDataset;
		pDouble += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pDouble++)
			pToRet[Idx] = (short)*pDouble;
		break;

	default:
		return(eBSFerrUserType);
	}
return(eBSFSuccess);
}

int
CDataPoints::GetValuesInt(int RelDatasetID,		// relative dataset identifier
						  int RefChromID,		// reference chromosome identifier
						  int RefChromOfs,		// start on reference chromosome (0..n) 
						  int NumPts,			// number of points to return
						  int *pToRet,			// where to return ints
						  int MissingMarker)	// value to use if missing datapoints
{
int Idx;
char *pChar;
short *pWord;
int *pDWord;
double *pDouble;
int Ofs;
int iRslt;

if(m_hFile == -1)
	return(eBSFerrClosed);
if(RefChromID <= eBSFSuccess || RefChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
if(RelDatasetID <= eBSFSuccess || RelDatasetID > m_FileHdr.NumDatasets)
	return(eBSFerrDataset);

if((iRslt=PreloadDataset(RelDatasetID,RefChromID,RefChromOfs,NumPts,&MissingMarker)) < eBSFSuccess)
	return(iRslt);

Ofs = RefChromOfs - m_LoadedDataChromOfs;
switch(m_FileHdr.DataPtType) {
	case eDPTPackedBases:		// bases 
	case eDPTbyte:				// 8 bit byte or char
		pChar = (char *)m_pCachedDataset;
		pChar += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pChar++)
			pToRet[Idx] = (int)*pChar;
		break;

	case eDPTword:				// 16 bit
		pWord = (short *)m_pCachedDataset;
		pWord += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pWord++)
			pToRet[Idx] = (int)*pWord;
		break;

	case eDPTdword:				// 32 bit
		pDWord = (int *)m_pCachedDataset;
		pDWord += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pDWord++)
			pToRet[Idx] = (int)*pDWord;
		break;

	case eDPTdouble:			// double
		pDouble = (double *)m_pCachedDataset;
		pDouble += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pDouble++)
			pToRet[Idx] = (int)*pDouble;
		break;

	default:
		return(eBSFerrUserType);
	}
return(eBSFSuccess);
}

// returns values as doubles
int
CDataPoints::GetValuesDouble(int RelDatasetID,		// relative dataset identifier
						int RefChromID,			// reference chromosome identifier
						int RefChromOfs,		// start on reference chromosome (0..n) 
						int NumPts,				// number of points to return
					  double *pToRet,			// where to return doubles
					double MissingMarker)   // value to use if missing datapoints
{
int Idx;
char *pChar;
short *pWord;
int *pDWord;
double *pDouble;
int Ofs;
int iRslt;

if(m_hFile == -1)
	return(eBSFerrClosed);
if(RefChromID <= eBSFSuccess || RefChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
if(RelDatasetID <= eBSFSuccess || RelDatasetID > m_FileHdr.NumDatasets)
	return(eBSFerrDataset);

if((iRslt=PreloadDataset(RelDatasetID,RefChromID,RefChromOfs,NumPts,&MissingMarker)) < eBSFSuccess)
	return(iRslt);

Ofs = RefChromOfs - m_LoadedDataChromOfs;
switch(m_FileHdr.DataPtType) {
	case eDPTPackedBases:		// why would anyone want bases returned as doubles?		
	case eDPTbyte:				// 8 bit byte or char
		pChar = (char *)m_pCachedDataset;
		pChar += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pChar++)
			pToRet[Idx] = (double)*pChar;
		break;

	case eDPTword:				// 16 bit
		pWord = (short *)m_pCachedDataset;
		pWord += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pWord++)
			pToRet[Idx] = (double)*pWord;
		break;

	case eDPTdword:				// 32 bit
		pDWord = (int *)m_pCachedDataset;
		pDWord += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pDWord++)
			pToRet[Idx] = (double)*pDWord;
		break;

	case eDPTdouble:			// double
		pDouble = (double *)m_pCachedDataset;
		pDouble += Ofs;
		for(Idx = 0; Idx < NumPts; Idx++,pDouble++)
			pToRet[Idx] = *pDouble;
		break;

	default:
		return(eBSFerrUserType);
	}
return(eBSFSuccess);
}

//
// Preload and cache dataset
int
CDataPoints::PreloadDataset(int RelDatasetID,	// relative dataset identifier
						int RefChromID,			// reference chromosome identifier
						int RefChromOfs,		// start on reference chromosome 
						int NumPts,				// number of datapoints to get
						void *pMarker)			// (optional) pts to value to use as missing data marker
{
int Num2Alloc;
etSeqBase BaseMarker;
char ByteMarker;
short WordMarker;
int DWordMarker;
double DoubleMarker;
int Mem2Hold;
int iRslt;
bool bNewMarker;

if(m_hFile == -1)
	return(eBSFerrClosed);

if(RefChromID <= eBSFSuccess || RefChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
if(RelDatasetID <= eBSFSuccess || RelDatasetID > m_FileHdr.NumDatasets)
	return(eBSFerrDataset);

if(pMarker == NULL)
	{
	switch(m_FileHdr.DataPtType) {
		case eDPTPackedBases:		// bases packed 2 per byte
			BaseMarker = eBaseUndef;
			pMarker = &BaseMarker;
			break;

		case eDPTbyte:				// 8 bit byte or char
			ByteMarker = SCHAR_MIN;
			pMarker = &ByteMarker;
			break;

		case eDPTword:				// 16 bit
			WordMarker = SHRT_MIN;
			pMarker = &WordMarker;
			break;

		case eDPTdword:				// 32 bit
			DWordMarker = INT_MIN;
			pMarker = &DWordMarker;
			break;

		case eDPTdouble:			// double
			DoubleMarker = DBL_MIN;
			pMarker = &DoubleMarker;
			break;

		default:				// oops..
			return(eBSFerrTypeChange);
		}
	}

// check if marker has changed
bNewMarker = false;			// default is that marker has not changed
switch(m_FileHdr.DataPtType) {
	case eDPTPackedBases:		// bases packed 2 per byte
		bNewMarker = *(etSeqBase *)pMarker != m_LoadedDataMarker.Base ? true : false;
		m_LoadedDataMarker.Base = *(etSeqBase *)pMarker;
		break;

	case eDPTbyte:				// 8 bit byte or char
		bNewMarker = *(unsigned char *)pMarker != m_LoadedDataMarker.Byte ? true : false;
		m_LoadedDataMarker.Byte = *(unsigned char *)pMarker;
		break;

	case eDPTword:				// 16 bit
		bNewMarker = *(unsigned short int *)pMarker != m_LoadedDataMarker.Word ? true : false;
		m_LoadedDataMarker.Word = *(unsigned short int *)pMarker;
		break;

	case eDPTdword:				// 32 bit
		bNewMarker = *(UINT32 *)pMarker != m_LoadedDataMarker.DWord ? true : false;
		m_LoadedDataMarker.DWord = *(UINT32 *)pMarker;
		break;

	case eDPTdouble:			// double
		bNewMarker = *(double *)pMarker != m_LoadedDataMarker.Double ? true : false;
		m_LoadedDataMarker.Double = *(double *)pMarker;
		break;

	default:				// oops..
		return(eBSFerrTypeChange);
	}

// if dataset covering requested chrom + ofs + numdps is already cached then return
if(m_pCachedDataset != NULL &&
   !bNewMarker &&
   m_LoadedDataChromID == RefChromID && 
   m_LoadedDataDatasetID == RelDatasetID &&
   RefChromOfs >= m_LoadedDataChromOfs && 
	(RefChromOfs + NumPts) <= (m_LoadedDataChromOfs + m_LoadedDataNumPts))
		return(eBSFSuccess);


Mem2Hold = NumPts * m_FileHdr.DataPtSize;
if(m_pCachedDataset == NULL || Mem2Hold > m_AllocdCacheSize)
	{
	if(m_pCachedDataset != NULL)
		delete (char *)m_pCachedDataset;
	Num2Alloc = (NumPts * 3)/2;
	if(Num2Alloc < 100000)
		Num2Alloc = 100000;
	if(Num2Alloc > 1000000)
		Num2Alloc = NumPts;
	Mem2Hold = Num2Alloc * m_FileHdr.DataPtSize;
	if((m_pCachedDataset = new unsigned char [Mem2Hold])==NULL)
		{
		AddErrMsg("CDataPoints::PreloadDataset","Unable to alloc memory (%d byte) for packed base segment cache in file %s ",Mem2Hold,m_szFile);
		return(eBSFerrMem);
		}
	m_AllocdCacheSize = Num2Alloc;
	}


if((iRslt=GetDataPoints(RelDatasetID,RefChromID,RefChromOfs,NumPts,m_pCachedDataset,pMarker)) == eBSFSuccess)
	{
	m_LoadedDataChromID = RefChromID;
	m_LoadedDataDatasetID = RelDatasetID;
	m_LoadedDataChromOfs = RefChromOfs; 
	m_LoadedDataNumPts = NumPts;
	return(eBSFSuccess);
	}
m_LoadedDataChromID = 0;
m_LoadedDataDatasetID = 0;
m_LoadedDataChromOfs = 0; 
m_LoadedDataNumPts = 0;
return(iRslt);
}

int												// number of datapoints returned
CDataPoints::GetDataPoints(int RelDatasetID,	// relative dataset identifier
						int RefChromID,			// reference chromosome identifier
						int RefChromOfs,		// start on reference chromosome (0..n) 
						int NumPts,				// number of datapoints in dataset
						void *pRetDataPts,		// where to return datapoint values
						void *pMissingMarker)	// value to use if there are missing datapoints
{
sDSFDirEl *pDirEl;
unsigned char *pByte = (unsigned char *)pRetDataPts;			// just to make it a little easier instead of constant type casts..

INT64 FileOfs;
int Len;
int Num2Read;
int Num2Skip;
int NumMarkers;

if(!m_bOpenRead)
	return(eBSFerrRead);

if(RefChromID <= eBSFSuccess || RefChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
if(RelDatasetID <= eBSFSuccess || RelDatasetID > m_FileHdr.NumDatasets)
	return(eBSFerrDataset);
	
while(NumPts) {
	// get ptr to segment which either contains the starting ChromPsn or immediately follows 
	pDirEl = LocateRelDirEl(RelDatasetID,RefChromID,RefChromOfs,true);
	if(pDirEl == NULL)
		{
		NumMarkers = NumPts;
		break;
		}
	if(pDirEl->RefChromOfs > RefChromOfs)
		{
		NumMarkers = pDirEl->RefChromOfs - RefChromOfs;
		if(NumMarkers > (int)NumPts)
			NumMarkers = NumPts;
		FillDataPts(pByte,pMissingMarker,m_FileHdr.DataPtSize,NumMarkers);
		NumPts -= NumMarkers;
		if(!NumPts)
			break;
		RefChromOfs += NumMarkers;
		pByte += (NumMarkers * m_FileHdr.DataPtSize);
		}
	
	Num2Skip = RefChromOfs - pDirEl->RefChromOfs;
	Num2Read = pDirEl->NumPts - Num2Skip;
	if(Num2Read > NumPts)
		Num2Read = NumPts;

	FileOfs = (INT64)(UINT64)pDirEl->LoFileOfs + ((UINT64)pDirEl->HiFileOfs << 32L);

	if(m_FileHdr.DataPtType != eDPTPackedBases)
		{
		FileOfs += (m_FileHdr.DataPtSize * Num2Skip);
		Len = m_FileHdr.DataPtSize * Num2Read;
		}
	else
		{
		FileOfs += Num2Skip/2;
		Len = (Num2Read+1)/2;
		if(Num2Skip & 0x01 && !(Num2Read & 0x01))
			Len+=1;
		}
	if(_lseeki64(m_hFile,FileOfs,SEEK_SET)!=FileOfs ||			
		Len != read(m_hFile,pByte,Len))
		{
		AddErrMsg("CDataPoints::GetDataPoints","Unable to read data segment from opened file  %s - error %s",m_szFile,strerror(errno));
		Close();			// closes opened file..
		return(eBSFerrFileAccess);
		}
	
	NumPts -= Num2Read;
	RefChromOfs += Num2Read;
	if(m_FileHdr.DataPtType != eDPTPackedBases)
		pByte += Len;
	else
		{
		UnpackBases(Num2Read,			// unpacked sequence length 
		  pByte,						// pts to sequence to unpack into
		  Num2Skip & 0x01,				// nibble to start unpack from (0..1)
		  pByte);						// where to unpack from
		pByte += Num2Read;
		}
	}
if(NumPts)
	FillDataPts(pByte,pMissingMarker,m_FileHdr.DataPtSize,NumMarkers);
return(eBSFSuccess);
}


// UnpackBases
// unpack (packed into nibbles) from pPacked into pUnpack
// if pUnpacked pts to pPacked then inplace unpacking will be performed
int
CDataPoints::UnpackBases(UINT32 SeqLen,	// unpacked sequence length 
		  etSeqBase *pUnpacked,				// pts to sequence to unpack into
		  UINT32 NibbleOfs,			// nibble to start unpack from (0..1)
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


// PackBases
// Pack bases (pUnpacked) into nibbles (pPacked) with bits 0..3 holding pUnpacked[n] and bits 4..7 holding pUnpacked[n+1]
// if pUnpacked pts to pPacked then inplace packing will be performed
int
CDataPoints::PackBases(UINT32 SeqLen,	// unpacked sequence length 
		  etSeqBase *pUnpacked,				// pts to sequence to pack from
		  UINT32 NibbleOfs,			// nibble to start pack into (0..SeqLen-1)
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



int 
CDataPoints::CompareDirEls( const void *arg1, const void *arg2 )
{
sDSFDirEl *pEl1 = (sDSFDirEl *)arg1;
sDSFDirEl *pEl2 = (sDSFDirEl *)arg2;
if(pEl1->RelDatasetID == pEl2->RelDatasetID )
	{
	if(pEl1->RefChromID < pEl2->RefChromID)
		return(-1);
	if(pEl1->RefChromID > pEl2->RefChromID)
		return(1);
	if(pEl1->RefChromOfs < pEl2->RefChromOfs)
		return(-1);
	if(pEl1->RefChromOfs > pEl2->RefChromOfs)
		return(1);
	if(pEl1->NumPts < pEl2->NumPts)
		return(-1);
	if(pEl1->NumPts > pEl2->NumPts)
		return(1);
	return(0);
	}

return(pEl1->RelDatasetID < pEl2->RelDatasetID ? -1 : 1);
}






