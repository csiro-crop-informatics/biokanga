#include "StdAfx.h"
#include <fcntl.h>      /* Needed only for _O_RDWR definition */
#include <io.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include "./conservlib.h"
#include "./visdata.h"

CVisData::CVisData(void)
{
m_hFile = -1;
m_pSeg = NULL;
}

CVisData::~CVisData(void)
{
if(m_hFile != -1)
	close(m_hFile);

if(m_pSeg)
	delete m_pSeg;
}

// Reset()
// Resets data points file class context back to that immediately following instantiation
bool
CVisData::Reset(bool bFlush)		// true (default) is to write any pending header writes to disk before closing opened file
{
if(m_hFile != -1)
	{
	if(bFlush)
		Flush2Disk(false);
	close(m_hFile);
	m_hFile = -1;
	}
m_szFile[0] = '\0';
InitHdr();
m_HdrDirty = false;
m_SegDirty = false;

if(m_pSeg)
	{
	delete m_pSeg;
	m_pSeg = NULL;
	}
m_AllocdSegSize = 0;
return(true);
}

bool
CVisData::InitHdr(void)
{
memset(&m_VisHdr,0,sizeof(m_VisHdr));
m_VisHdr.Magic[0] = 'b';
m_VisHdr.Magic[1] = 'i';
m_VisHdr.Magic[2] = 'o';
m_VisHdr.Magic[3] = 's';
m_VisHdr.Type =		cBSFTypeDataPts;			// it's a visualisation dataset file
m_VisHdr.Version =	cVisDataPtVersion;		    // file structure version
m_VisHdr.SizeOfHdr =sizeof(tsVisHdr);			// size of this header
m_VisHdr.FileLen =	sizeof(tsVisHdr);			// current file length
m_VisHdr.DataType = eDataPts;					// store points as doubles
m_VisHdr.MaxDatasets=cMaxVisDatasets;			// max number of datasets supported
m_VisHdr.MaxChroms = cMaxVisChroms;				// max number of chromosomes supported
m_VisHdr.MaxDatasegs = cMaxVisDatasegs;			// max number of datasegs
m_VisHdr.MaxDataPts = cMaxVisDataPts;
m_HdrDirty = true;
return(true);
}

bool
CVisData::Close(void)
{
return(Reset(true));
}



bool 
CVisData::Open(char *pszVisFile,	// specifies file to open or create
			   bool bCreate)		// create file or truncate if it already exists
{
Reset(false);						// reset context in case file still opened

if(pszVisFile == NULL || *pszVisFile == '\0') // validate parameters
	{
	CDiagnostics::DiagOut(eDLFatal,"CVisFile::Open","No data point file specified!");
	return(false);
	}

if(!bCreate)
	m_hFile = open(pszVisFile, _O_RDWR | _O_BINARY ); // file access is normally random..
else
	m_hFile = open(pszVisFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE );

if(m_hFile == -1)					// check if file open succeeded
	{
	if(!bCreate)
		CDiagnostics::DiagOut(eDLFatal,"CVisFile::Open","Unable to open file - '%s' - %s\n",	pszVisFile,strerror(errno));
	else
		CDiagnostics::DiagOut(eDLFatal,"CVisFile::Open","Unable to create or truncate file - '%s' - %s\n",	pszVisFile,strerror(errno));
	Reset(false);
	return(false);
	}

strncpy(m_szFile,pszVisFile,_MAX_PATH);
m_szFile[_MAX_PATH-1] = '\0';

if(bCreate)
	{
	InitHdr();
	m_SegDirty = false;
	Flush2Disk(false);
	}
else // else opening existing file
	{
	_lseeki64(m_hFile,0,SEEK_SET);			// read in header..
	if(sizeof(tsVisHdr) != read(m_hFile,&m_VisHdr,sizeof(tsVisHdr)))
		{
		CDiagnostics::DiagOut(eDLFatal,"CVisData::Open","unable to read complete file header from opened file %s",pszVisFile);
		Reset(false);			// closes opened file..
		return(false);
		}

	// header read, validate it as being a bioseq visualisation data point file header
	if(tolower(m_VisHdr.Magic[0]) != 'b' || 
		tolower(m_VisHdr.Magic[1]) != 'i' || 
		tolower(m_VisHdr.Magic[2]) != 'o' || 
		tolower(m_VisHdr.Magic[3]) != 's')
		{
		CDiagnostics::DiagOut(eDLFatal,"CVizData::Open","opened file %s not a visualisation data point file (invalid magic signiture)",pszVisFile);
		Reset(false);			// closes opened file..
		return(false);
		}

	// is it a data point file?
	if(m_VisHdr.Type != cBSFTypeDataPts)
		{
		CDiagnostics::DiagOut(eDLFatal,"CVizData::Open","file %s does not contain visualisation data point",pszVisFile);
		Reset(false);			// closes opened file..
		return(false);
		}

	// can we handle this version?
	if(m_VisHdr.Version != cVisDataPtVersion)
		{
		CDiagnostics::DiagOut(eDLFatal,"CVizData::Open","file %s structure created as version %d which can't be processed by this application (version %d) - upgrade required",
			pszVisFile,m_VisHdr.Version,cVisDataPtVersion);
		Reset(false);			// closes opened file..
		return(false);
		}
	}

return(true);
}

bool 
CVisData::Flush2Disk(bool BlockOnly)						// flush and commit to disk
{
if(m_hFile == -1)
	return(true);
if(m_HdrDirty) {
	_lseeki64(m_hFile,0,SEEK_SET);
	if(write(m_hFile,&m_VisHdr,sizeof(tsVisHdr))!=sizeof(tsVisHdr))
		{
		CDiagnostics::DiagOut(eDLFatal,"CVisData::Flush2Disk","unable to write file header");
		return(false);
		}
	}
_commit(m_hFile);
return(true);
}


// GenNameHash
// Generates a 16bit hash on specified lowercased name
// This hash can then be used to quickly eliminate probe names which can't match a target name by comparing hashes
unsigned __int16 
CVisData::GenNameHash(TCHAR *pszName)
{
unsigned __int16 Hash = 0;
unsigned __int16 MSB;
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

// GetDatasetID
// Returns GetDatasetID which is used to identify the dataset from it's name
// 0 is returned if unable to locate dataset
unsigned __int16
CVisData::GetDatasetID(TCHAR *pszName)
{
tsVisDatasetName *pDataset;
unsigned int Idx;
unsigned __int16 Hash;
Hash = GenNameHash(pszName);
pDataset = &m_VisHdr.Datasets[0];
for(Idx = 0; Idx < m_VisHdr.NumDatasets; Idx++, pDataset++)
	{
	if(pDataset->Hash == Hash &&
		!stricmp(pDataset->szName,pszName))
		return(pDataset->DatasetID);
	}
return(0);
}

unsigned __int16
CVisData::GetChromID(TCHAR *pszName)
{
tsVisChrom *pChrom;
unsigned int Idx;
unsigned __int16 Hash;
Hash = GenNameHash(pszName);
pChrom = &m_VisHdr.Chroms[0];
for(Idx = 0; Idx < m_VisHdr.NumChroms; Idx++,pChrom++)
	{
	if(pChrom->Hash == Hash &&
		!stricmp(pChrom->szName,pszName))
			return(pChrom->ChromID);
	}
return(0);
}


bool
CVisData::AddDataseg(TCHAR *pszDataset,
					 TCHAR *pszChrom,
					 unsigned __int32 ChromOfs,	// start on chromsome 
					 unsigned int NumPts,		// number of datapoints
					 double *pDataPts)			// array of data points
{
tsVisChrom *pChrom;
tsVisDatasetName *pDatasetName;
tsDataset *pDataset;
tsDataseg LastSeg;
unsigned __int32 CpyLen;
unsigned __int16 DatasetID;
unsigned __int16 ChromID;
INT64 FileOfs;
int Len;

if(m_pSeg == NULL)
	return(false);


ChromID = GetChromID(pszChrom);
if(!ChromID)
	{
	if(m_VisHdr.NumChroms >= cMaxVisChroms)
		return(0);
	pChrom = &m_VisHdr.Chroms[m_VisHdr.NumChroms++];
	memset(pChrom,0,sizeof(tsVisChrom));
	ChromID = pChrom->ChromID = m_VisHdr.NumChroms;
	strncpy(pChrom->szName,pszChrom,cMaxVizChromNameLen);
	pChrom->szName[cMaxVizChromNameLen-1] = _T('\0');
	pChrom->Hash = GenNameHash(pszChrom);
	}
else
	pChrom = &m_VisHdr.Chroms[ChromID-1];

DatasetID = GetDatasetID(pszDataset);
if(!DatasetID)		// if a new dataset
	{
	if(m_VisHdr.NumDatasets >= cMaxVisDatasets)
		return(0);
	pDatasetName = &m_VisHdr.Datasets[m_VisHdr.NumDatasets++];
	memset(pDatasetName,0,sizeof(tsVisDatasetName));
	DatasetID = pDatasetName->DatasetID = m_VisHdr.NumDatasets;
	strncpy(pDatasetName->szName,pszDataset,cMaxVizDatasetNameLen);
	pDatasetName->szName[cMaxVizDatasetNameLen-1] = _T('\0');
	pDatasetName->Hash = GenNameHash(pszDataset);
	}
pDataset = &pChrom->Datasets[DatasetID-1];
if(pDataset->NumSegs)
	{
	FileOfs = pDataset->LastSeg;
	Len = offsetof(tsDataseg,DataPts);
	_lseeki64(m_hFile,FileOfs,SEEK_SET);			
	if(Len != read(m_hFile,&LastSeg,Len))
		{
		CDiagnostics::DiagOut(eDLFatal,"CVisData::AddDataseg","unable to read dataseg header");
		Reset(false);			// closes opened file..
		return(false);
		}
	LastSeg.Nxt = m_VisHdr.FileLen;
	_lseeki64(m_hFile,FileOfs,SEEK_SET);
	write(m_hFile,&LastSeg,Len);
	m_pSeg->Prv = FileOfs;
	}
else
	{
	pDataset->FirstChromOfs = INT_MAX;
	pDataset->LastChromOfs = INT_MIN;
	pDataset->NumSegs = 0;
	pDataset->FirstSeg = m_VisHdr.FileLen; 
	m_pSeg->Prv = 0;
	}
pDataset->LastSeg = m_VisHdr.FileLen;
if(pDataset->FirstChromOfs > ChromOfs)
	pDataset->FirstChromOfs = ChromOfs;
if(pDataset->LastChromOfs < ChromOfs)
	pDataset->LastChromOfs = ChromOfs;
pDataset->NumSegs++;	
m_pSeg->Nxt = 0;
m_pSeg->ChromID = ChromID;
m_pSeg->DatasetID = DatasetID;
m_pSeg->ChromOfs = ChromOfs;
m_pSeg->NumPts = NumPts;
CpyLen = NumPts * sizeof(double);
memcpy(m_pSeg->DataPts,pDataPts,CpyLen);
m_pSeg->Size = offsetof(tsDataseg,NumPts) + CpyLen;

_lseeki64(m_hFile,m_VisHdr.FileLen,SEEK_SET);
write(m_hFile,m_pSeg,m_pSeg->Size);

m_VisHdr.FileLen += m_pSeg->Size;
m_HdrDirty = true;
return(true);
}

bool 
CVisData::GetDataPoints(CHAR *pszDataset,			// dataset name
						TCHAR *pszChrom,
						unsigned __int32 ChromOfs,	// start on chromsome 
						unsigned int NumPts,		// number of datapoints in dataset
						double MissingMarker,		// value to use if there are missing datapoints
						double *pRetValues)			// where to return datapoint values
{
unsigned __int16 DatasetID;
unsigned __int16 ChromID;
INT64 FileOfs;
int Len;
unsigned int Num2Read;
unsigned int Num2Skip;
int NumMarkers;
tsDataseg LastSeg;
tsVisChrom *pChrom;
tsDataset *pDataset;
if(!(ChromID = GetChromID(pszChrom)))
	{
	return(false);
	}

if(!(DatasetID = GetDatasetID(pszDataset)))
	{
	return(false);
	}

pChrom = &m_VisHdr.Chroms[ChromID-1];
pDataset = &pChrom->Datasets[DatasetID-1];
LastSeg.Nxt = pDataset->FirstSeg;

while(NumPts) {
	FileOfs = LastSeg.Nxt;
	if(FileOfs == 0)
		break;
	Len = offsetof(tsDataseg,DataPts);
	_lseeki64(m_hFile,FileOfs,SEEK_SET);			
	if(Len != read(m_hFile,&LastSeg,Len))
		{
		CDiagnostics::DiagOut(eDLFatal,"CVisData::AddDataseg","unable to read dataseg header");
		Reset(false);			// closes opened file..
		return(false);
		}

	if((LastSeg.ChromOfs + LastSeg.NumPts - 1) < ChromOfs)
		continue;
	if(LastSeg.ChromOfs >= ChromOfs + NumPts)
		break;

	if(LastSeg.ChromOfs > ChromOfs)		// need to do some filling?
		{
		NumMarkers = LastSeg.ChromOfs - ChromOfs;
		NumPts -= NumMarkers;
		ChromOfs += NumMarkers;
		while(NumMarkers--)
			*pRetValues++ = MissingMarker;
		}

	Num2Skip = ChromOfs - LastSeg.ChromOfs;
	Num2Read = LastSeg.NumPts - Num2Skip;
	if(Num2Read > NumPts)
		Num2Read = NumPts;

	FileOfs += Len + (sizeof(double) * Num2Skip);
	Len = sizeof(double) * Num2Read;
	_lseeki64(m_hFile,FileOfs,SEEK_SET);			
	if(Len != read(m_hFile,pRetValues,Len))
		{
		CDiagnostics::DiagOut(eDLFatal,"CVisData::AddDataseg","unable to read dataseg header");
		Reset(false);			// closes opened file..
		return(false);
		}
	
	NumPts -= Num2Read;
	pRetValues += Num2Read;
	}

while(NumPts--)
	*pRetValues++ = MissingMarker;
return(true);
}

