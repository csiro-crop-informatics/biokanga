// Copyright 2013 CSIRO  ( http://www.csiro.au/ ) 
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License
//   Please contact stuart.stephen@csiro.au for support or 
//   to submit modifications to this source

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include <sys/mman.h>
#include "./commhdrs.h"
#endif

CBEDfile::CBEDfile(void)
{
short word = 0x4321;
m_pChromNames = NULL;		// pts to array of tsBEDchromname's sorted by chromid
m_pFeatures = NULL;			// pts to array of tsBEDfeature's sorted by name-->chromid-->start-->end
m_ppFeatureNames = NULL; // sorted (by name->chrom->start->end) array of ptrs into m_pFeatures
m_ppFeatureChromStarts = NULL; // sorted (by chrom->start->end) array of ptrs into m_pFeatures
m_pChromHashes = NULL;
m_hFile = -1;
Reset(false);
}

CBEDfile::~CBEDfile(void)
{
Reset(false);
}

teBSFrsltCodes
CBEDfile::Disk2Hdr(char *pszBioBed,teBEDFeatureType FeatType)
{
if(_lseeki64(m_hFile,0,SEEK_SET)!=0)			// read in header..
	{
	AddErrMsg("CBEDfile::Disk2Hdr","Seek failed to offset 0 on %s - %s",pszBioBed,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

if(sizeof(tsBEDFileHdr) != read(m_hFile,&m_FileHdr,sizeof(tsBEDFileHdr)))
	{
	AddErrMsg("CBEDfile::Disk2Hdr","Read of file header failed on %s - %s",pszBioBed,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

	// header read, validate it as being a BED file header
if(tolower(m_FileHdr.Magic[0]) != 'b' || 
	tolower(m_FileHdr.Magic[1]) != 'i' || 
	tolower(m_FileHdr.Magic[2]) != 'o' || 
	tolower(m_FileHdr.Magic[3]) != 's')
	{
	AddErrMsg("CBEDfile::Disk2Hdr","%s opened but no magic signature - not a biobed file",pszBioBed);
	Reset(false);			// closes opened file..
	return(eBSFerrNotBioseq);
	}



if(m_bIsBigEndian)	// file was written with little-endian ordering
	{
	m_FileHdr.FileLen = SwapUI64Endians(m_FileHdr.FileLen);				// current file length
	m_FileHdr.ChromNamesOfs = SwapUI64Endians(m_FileHdr.ChromNamesOfs);	// file offset to chromosome names
	m_FileHdr.FeaturesOfs  = SwapUI64Endians(m_FileHdr.FeaturesOfs);		// file offset to features
	m_FileHdr.FeatureNamesOfs  = SwapUI64Endians(m_FileHdr.FeatureNamesOfs);	// file offset to sorted (by name->chrom->start->end)
	m_FileHdr.FeatureChromStartsOfs  = SwapUI64Endians(m_FileHdr.FeatureChromStartsOfs);		// file offset to sorted (by chrom->start->end)
	m_FileHdr.Type  = SwapUI32Endians(m_FileHdr.Type);						// biosequence file type 
	m_FileHdr.Version  = SwapUI32Endians(m_FileHdr.Version);				// header version, incremented if structure changes with later releases
	m_FileHdr.SizeOfHdr  = SwapUI32Endians(m_FileHdr.SizeOfHdr);			// total size of this header
	m_FileHdr.MaxChroms  = SwapUI32Endians(m_FileHdr.MaxChroms);			// maximum number of chromosomes supported
	m_FileHdr.NumChroms  = SwapUI32Endians(m_FileHdr.NumChroms);			// actual number of chromosomes
	m_FileHdr.MaxFeatures  = SwapUI32Endians(m_FileHdr.MaxFeatures);		// maximum number of features supported
	m_FileHdr.NumFeatures  = SwapUI32Endians(m_FileHdr.NumFeatures);		// actual number of features
	m_FileHdr.FeatType  = SwapUI32Endians(m_FileHdr.FeatType);				// what type of features are in this file
	m_FileHdr.FeaturesSize  = SwapUI32Endians(m_FileHdr.FeaturesSize);		// disk/memory space required to hold concatenated Features
	m_FileHdr.ChromNamesSize  = SwapUI32Endians(m_FileHdr.ChromNamesSize);	// disk/memory space required to hold concatenated chromosome names
	}

	// check bioseq file is the type we are expecting
if(m_FileHdr.Type != cBSFTypeFeat)
	{
	AddErrMsg("CBEDfile::Disk2Hdr","%s opened as a bioseq file - expected type %d, file type is %d",pszBioBed,cBSFTypeFeat,m_FileHdr.Type);
	Reset(false);			// closes opened file..
	return(eBSFerrFileType);
	}

if(FeatType != eBTAnyBed && m_FileHdr.FeatType != FeatType)
	{
	AddErrMsg("CBEDfile::Disk2Hdr","%s opened as a bioseq file - expected type %d, file type is %d",pszBioBed,FeatType,m_FileHdr.FeatType);
	Reset(false);			// closes opened file..
	return(eBSFerrFileType);
	}

	// can we handle this version?
if(m_FileHdr.Version > cBSFeatVersion || m_FileHdr.Version < cBSFFeatVersionBack)
	{
	AddErrMsg("CBEDfile::Disk2Hdr","%s opened as a bioseq file - can only handle versions %d to %d, file version is %d",pszBioBed,
			cBSFFeatVersionBack,cBSFeatVersion,m_FileHdr.Version);
	Reset(false);			// closes opened file..
	return(eBSFerrFileVer);
	}
return(eBSFSuccess);
}

teBSFrsltCodes
CBEDfile::Hdr2Disk(void)
{
tsBEDFileHdr FileHdr;
tsBEDFileHdr *pHdr;
int WrtLen;

WrtLen = sizeof(tsBEDFileHdr);

if(m_bIsBigEndian)	// if on a big-endian machine then need to make little endian as that is our native file format
	{
	memmove(&FileHdr,&m_FileHdr,WrtLen);
	FileHdr.FileLen = SwapUI64Endians(m_FileHdr.FileLen);				// current file length
	FileHdr.ChromNamesOfs = SwapUI64Endians(m_FileHdr.ChromNamesOfs);	// file offset to chromosome names
	FileHdr.FeaturesOfs  = SwapUI64Endians(m_FileHdr.FeaturesOfs);		// file offset to features
	FileHdr.FeatureNamesOfs  = SwapUI64Endians(m_FileHdr.FeatureNamesOfs);	// file offset to sorted (by name->chrom->start->end)
	FileHdr.FeatureChromStartsOfs  = SwapUI64Endians(m_FileHdr.FeatureChromStartsOfs);		// file offset to sorted (by chrom->start->end)
	FileHdr.Type  = SwapUI32Endians(m_FileHdr.Type);						// biosequence file type 
	FileHdr.Version  = SwapUI32Endians(m_FileHdr.Version);				// header version, incremented if structure changes with later releases
	FileHdr.SizeOfHdr  = SwapUI32Endians(m_FileHdr.SizeOfHdr);			// total size of this header
	FileHdr.MaxChroms  = SwapUI32Endians(m_FileHdr.MaxChroms);			// maximum number of chromosomes supported
	FileHdr.NumChroms  = SwapUI32Endians(m_FileHdr.NumChroms);			// actual number of chromosomes
	FileHdr.MaxFeatures  = SwapUI32Endians(m_FileHdr.MaxFeatures);		// maximum number of features supported
	FileHdr.NumFeatures  = SwapUI32Endians(m_FileHdr.NumFeatures);		// actual number of features
	FileHdr.FeatType  = SwapUI32Endians(m_FileHdr.FeatType);				// what type of features are in this file
	FileHdr.FeaturesSize  = SwapUI32Endians(m_FileHdr.FeaturesSize);		// disk/memory space required to hold concatenated Features
	FileHdr.ChromNamesSize  = SwapUI32Endians(m_FileHdr.ChromNamesSize);	// disk/memory space required to hold concatenated chromosome names
	pHdr = &FileHdr;
	}
else
	pHdr = &m_FileHdr;

if(_lseeki64(m_hFile,0,SEEK_SET) ||
			write(m_hFile,pHdr,WrtLen)!=WrtLen)
	{
	AddErrMsg("CBEDfile::Flush2Disk","Unable to write file header to disk on file %s - error %s",m_szFile,strerror(errno));
	Reset(false);
	return(eBSFerrFileAccess);
	}
m_bHdrDirty = false;
return(eBSFSuccess);
}


void
CBEDfile::Reset(bool bFlush)
{
if(bFlush)
	Flush2Disk();

if(m_hFile != -1)
	{
	close(m_hFile);
	m_hFile = -1;
	}

if(m_pChromNames != NULL)
	{
#ifdef _WIN32
	free(m_pChromNames);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pChromNames != MAP_FAILED)
		munmap(m_pChromNames,m_AllocChromNamesSize);
#endif
	m_pChromNames = NULL;
	m_AllocChromNamesSize = 0;
	}

if(m_pFeatures != NULL)
	{
#ifdef _WIN32
	free(m_pFeatures);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pFeatures != MAP_FAILED)
		munmap(m_pFeatures,m_AllocFeaturesSize);
#endif
	m_pFeatures = NULL;
	m_AllocFeaturesSize = 0;
	}

if(m_ppFeatureNames != NULL)
	{
	delete m_ppFeatureNames;
	m_ppFeatureNames = NULL;
	}

if(m_ppFeatureChromStarts != NULL)
	{
	delete m_ppFeatureChromStarts;
	m_ppFeatureChromStarts = NULL;
	}

if(m_pChromHashes != NULL)
	{
	delete m_pChromHashes;
	m_pChromHashes = NULL;
	}
m_pCacheChromName = NULL;
m_AllocChromNamesSize = 0;			// disk/memory space required to hold concatenated chromosome names
m_AllocFeaturesSize = 0;			// disk/memory space required to hold concatenated Features
m_MinScore = INT_MIN;
m_MaxScore = INT_MAX;
m_OnStrand = '*';

m_szFile[0] = '\0';				    // file name as opened/created
m_bCreate = false;					// TRUE if file opened for create or update 
m_bHdrDirty = false;				// TRUE if header needs to be written to file
m_bFeaturesAvail = false;			// features not available until read from disk

m_MaxParseFeats = 0;				// initially no restriction on number of features to parse
InitHdr();
}

void 
CBEDfile::InitHdr(void)
{
memset(&m_FileHdr,0,sizeof(m_FileHdr));
m_FileHdr.Magic[0] = 'b';
m_FileHdr.Magic[1] = 'i';
m_FileHdr.Magic[2] = 'o';
m_FileHdr.Magic[3] = 's';
m_FileHdr.Type = cBSFTypeFeat;			// biosequence file type 
m_FileHdr.Version = cBSFeatVersion;		// header version, incremented if structure changes with later releases
m_FileHdr.FileLen = sizeof(m_FileHdr);	// current file length
m_FileHdr.SizeOfHdr = sizeof(m_FileHdr);// total size of this header
m_FileHdr.MaxChroms = cMaxNumChroms;	// max number of chromosomes supported
m_FileHdr.MaxFeatures = cMaxNumFeats;	// maximum number of features supported
m_FileHdr.szDescription[0] = '\0';
m_FileHdr.szTitle[0] = '\0';
m_bHdrDirty = true;						// TRUE if header needs to be written to file
}


// Open
// Opens and loads -
// 1) Attempts to open as a preprocessed biobed formated file, if that fails then
// 2) Attempts to open as a UCSC BED format file, if that fails then
// 3) Attempts to open as a GFF3 format file  
teBSFrsltCodes
CBEDfile::Open(char *pszBioBed,					// BioBed, BED or GFF3 input file
				teBEDFeatureType FeatType,		// expected to contain 
				bool bCreate,					// if true then will be creating BioBed file with pszBioBed name
				int MaxParseFeats)				// parse at most this many features
{
teBSFrsltCodes Rslt;
if(pszBioBed == NULL || *pszBioBed == '\0') // validate parameters
	{
	AddErrMsg("CBEDfile::Open","Parameter errors - %s",pszBioBed);
	return(eBSFerrParams);
	}


Reset(false);						// reset context in case file was previously opened
m_MaxParseFeats = MaxParseFeats;

#ifdef _WIN32
if(!bCreate)
	m_hFile = open(pszBioBed, O_READSEQ ); // file access is normally sequential..
else
	m_hFile = open(pszBioBed, O_CREATETRUNC);
#else
if(!bCreate)
	m_hFile = open64(pszBioBed, O_READSEQ ); // file access is normally sequential..
else
	{
     if((m_hFile = open64(pszBioBed,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
          if(ftruncate(m_hFile,0)!=0)
			{
			AddErrMsg("CBEDfile::Open","Unable to truncate %s - %s",pszBioBed,strerror(errno));
			Rslt = bCreate ? eBSFerrCreateFile : eBSFerrOpnFile;
			return(Rslt);
			}
	}
#endif

if(m_hFile == -1)					// check if file open succeeded
	{
	AddErrMsg("CBEDfile::Open","Unable to open %s - %s",pszBioBed,strerror(errno));
	Rslt = bCreate ? eBSFerrCreateFile : eBSFerrOpnFile;

	return(Rslt);
	}

strncpy(m_szFile,pszBioBed,_MAX_PATH);
m_szFile[_MAX_PATH-1] = '\0';

if(m_pChromHashes == NULL)
	{
	if((m_pChromHashes = (UINT32 *) new UINT32 [ 0x1000000 ]) == NULL)
		{
		AddErrMsg("CBEDfile::Open","Unable to allocate memory for chrom hashes");
		Reset(false);
		return(eBSFerrMem);
		}
	memset(m_pChromHashes,0,sizeof(UINT32) * 0x1000000);
	}

if(bCreate)
	{
	m_bCreate = true;
	m_bFeaturesAvail = false;
	InitHdr();
	m_FileHdr.FeatType = FeatType;
	if((Rslt = Flush2Disk()) != eBSFSuccess)
		{
		Reset(false);			// closes opened file..
		return(Rslt);
		}
	}
else // else opening existing file which could be a bioseq bed, raw bed file or GFF3
	{
	// check first few lines of file for likely file format
	char szFileHdr[5];
	int FileHdrLen;
	if(_lseeki64(m_hFile,0,SEEK_SET)!=0)			// read in potential signiture
		{
		AddErrMsg("CBEDfile::Open","Seek failed to offset 0 on %s - %s",pszBioBed,strerror(errno));
		Reset(false);			// closes opened file..
		return(eBSFerrFileAccess);
		}
	if((FileHdrLen=read(m_hFile,szFileHdr,sizeof(szFileHdr))) < 4)
		{
		AddErrMsg("CBEDfile::Open","Read of file header failed on %s - %s",pszBioBed,strerror(errno));
		Reset(false);			// closes opened file..
		return(eBSFerrFileAccess);
		}
	if(!strnicmp(szFileHdr,"bios",4))		// signiture for BioBed format file
		{
		if((Rslt=Disk2Hdr(pszBioBed,FeatType))!=eBSFSuccess)
			{
			Reset(false);			// closes opened file..
			return(eBSFerrFileAccess);
			}

		// if not empty then load chromosome names and features + indexes..
		if(m_FileHdr.NumFeatures)
			{
			if((Rslt=LoadFeatures(true))!=eBSFSuccess)	// close file after features loaded
				{
				AddErrMsg("CBEDfile::Open","Error loading features from %s",pszBioBed);
				Reset(false);
				return(Rslt);
				}
			}
		}
	else		// else it's expected to be raw ascii BED or GFF3
		{
		// first attempt to parse as BED, if that fails then try to parse as GFF3
		close(m_hFile);
		m_hFile = -1;
		InitHdr();
		m_FileHdr.FeatType = FeatType;
		m_bHdrDirty = false;

		// initially assume BED
		if((Rslt=ProcessBedFile(pszBioBed))!=eBSFSuccess)	// file to process
			{
#ifdef PROCESSGFF3s
			if(Rslt != eBSFerrFileType)
#endif
				{
				AddErrMsg("CBEDfile::Open","Unable to process raw BED file %s",pszBioBed);
				Reset(false);
				return(Rslt);
				}
#ifdef PROCESSGFF3s
			// failed to parse as BED, try as GFF3
			if((Rslt=ProcessGFF3File(pszBioBed))!=eBSFSuccess)	// file to process
				{
				AddErrMsg("CBEDfile::Open","Unable to process raw BED or GFF3 file %s",pszBioBed);
				Reset(false);
				return(Rslt);
				}
#endif
			}

		if((Rslt=SortFeatures())!=eBSFSuccess)
			{
			AddErrMsg("CBEDfile::Flush2Disk","Unable to sort features %s",m_szFile);
			Reset(false);
			return(Rslt);
			}

		char szFname[_MAX_FNAME];
#ifdef _WIN32
		_splitpath(pszBioBed,NULL,NULL,szFname,NULL);
#else
		CUtility::splitpath(pszBioBed,NULL,szFname);
#endif
		strncpy(m_FileHdr.szTitle,szFname,sizeof(m_FileHdr.szTitle)-1);
		}
	m_bFeaturesAvail = true;		// features are in memory and can now be accessed
	}
return(eBSFSuccess);
}



// sets minimum score threshold
// forces minimum score to 0 if < 0
teBSFrsltCodes 
CBEDfile::SetMinScore(int Score)
{
m_MinScore = Score;
return(eBSFSuccess);
}

// sets maximum score threshold
teBSFrsltCodes 
CBEDfile::SetMaxScore(int Score)
{
m_MaxScore = Score;
return(eBSFSuccess);
}

// gets minimum score threshold
int 
CBEDfile::GetMinScore(void)
{
return(m_MinScore);
}

// gets maximum score threshold
int 
CBEDfile::GetMaxScore(void)
{
return(m_MaxScore);
}


teBSFrsltCodes 
CBEDfile::SetDescription(char *pszDescription)
{
if(m_FileHdr.Version >= 2)
	{
	strncpy(m_FileHdr.szDescription,pszDescription,sizeof(m_FileHdr.szDescription));
	m_FileHdr.szDescription[sizeof(m_FileHdr.szDescription)-1] = '\0';
	return(eBSFSuccess);
	}
return(eBSFerrFileVer);
}

teBSFrsltCodes 
CBEDfile::GetDescription(int MaxLen,char *pszDescription)
{
if(m_FileHdr.Version >= 2 && m_FileHdr.szDescription[0] != '\0')
	strncpy(pszDescription,m_FileHdr.szDescription,MaxLen);
else
	strncpy(pszDescription,m_szFile,MaxLen);
pszDescription[MaxLen-1] = '\0';
return(eBSFSuccess);
}

teBSFrsltCodes 
CBEDfile::SetTitle(char *pszTitle)
{
strncpy(m_FileHdr.szTitle,pszTitle,sizeof(m_FileHdr.szTitle));
m_FileHdr.szTitle[sizeof(m_FileHdr.szTitle)-1] = '\0';
return(eBSFSuccess);
}

teBSFrsltCodes 
CBEDfile::GetTitle(int MaxLen,char *pszTitle)
{
if(m_FileHdr.Version >= 2 && m_FileHdr.szTitle[0] != '\0')
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

char *
CBEDfile::TrimWhitespace(char *pTxt)
{
char *pStart;
char Chr;
	// strip leading whitespace
while(Chr = *pTxt++)
	if(!isspace(Chr))
			break;
if(Chr == '\0')					// empty line?
	return(pTxt-1);
pStart = pTxt-1;
while(Chr = *pTxt)			// fast forward to line terminator
	pTxt++;
pTxt-=1;
while(Chr = *pTxt--)
	if(!isspace(Chr))
		break;
pTxt[2] = '\0';
return(pStart);
}


typedef struct TAG_sGFF3FeatureType {
	int FeatTypeID;		// uniquely identifies this feature type
	const char *pszFeatType;	// feature type
	} tsGFF3FeatureType;

tsGFF3FeatureType GFF3FeatureTypes[] = {
	{1,"gene"},
	{2,"mRNA"},
	{3,"exon"},
	{4,"intron"},
	{5,"CDS"},
	{6,"three_prime_UTR"},
	{7,"five_prime_UTR"}
	};
const int cNumGFF3FeatTypes = (sizeof(GFF3FeatureTypes)/sizeof(tsGFF3FeatureType));

// defines GFF3 attributes of interest
typedef struct TAG_sGFF3AttribType {
	int AttribTypeID;		// uniquely identifies this attribute type
	const char *pszAttribType;	// attribute type
	} tsGFF3AttribType;

tsGFF3AttribType GFF3AttribTypes[] = {
	{1,"ID"},
	{2,"Parent"},
	{3,"Name"},
	};
const int cNumGFF3AttribTypes = (sizeof(GFF3AttribTypes)/sizeof(tsGFF3AttribType));

const int cMaxGFF3AttribValueLen = 200;	// truncate GFF3 attribute values to this length
typedef struct TAG_sGFF3AttribValue {
	int AttribTypeID;		// attribute type identifier
	char szAttribValue[cMaxGFF3AttribValueLen+1];	// value associated with this attribute
	} tsGFF3AttribValue;


// ProcessGFF3File
// Opens and parses GFF3 format file
// <seqid>\t<source>\t<type>\t<start>\t<end>\t<score>\t<strand>\t<phase>\t<attributes>
teBSFrsltCodes 
CBEDfile::ProcessGFF3File(char *pszFileName)		// process/parse as a GFF2 or GFF3 format file
{
teBSFrsltCodes Rslt;
FILE *pGFFStream;
int LineNum;
int NumFields;
int FeatNum;
int FeatTypeID;
tsGFF3FeatureType *pFeatType;

int AttribTypeID;
tsGFF3AttribValue AttribValues[cNumGFF3AttribTypes];
tsGFF3AttribType *pAttribType;
tsGFF3AttribValue *pAttribValue;

char *pTxt;
char szSeqName[cMaxFeatNameLen+1];
char szSource[cMaxFeatNameLen+1];
char szFeature[cMaxFeatNameLen+1];
char szScore[cMaxFeatNameLen+1];
char cStrand;
char cFrame;

int iScore;
int iFrame;
int ChromStart;
int ChromEnd;
int AttribStart;
int AttribLen;
int Idx;
char Chr;

if(pszFileName == NULL || *pszFileName == '\0')
	return(eBSFerrParams);
if((pGFFStream = fopen(pszFileName,"r"))==NULL)
	{
	AddErrMsg("CBEDfile::ProcessGFF23File","Unable to fopen GFF2/3 format file %s error: %s",pszFileName,strerror(errno));
	return(eBSFerrOpnFile);
	}
LineNum = 0;
FeatNum = 0;
while(fgets(m_szLineBuff,sizeof(m_szLineBuff)-1,pGFFStream)!= NULL)
	{
	LineNum += 1;
	if(!FeatNum && LineNum >= 20)	// if can't load at least one feature from 1st 20 lines then can't be GFF format
		{
		fclose(pGFFStream);
		return(eBSFerrFileType);
		}
	pTxt = TrimWhitespace(m_szLineBuff);
	if(*pTxt=='\0' || *pTxt=='#')	// simply slough lines which were just whitespace or start with '#'
		continue;
	NumFields = sscanf(pTxt," %50s %50s %50s %d %d %50s %c %c %n",
			szSeqName,szSource,szFeature,&ChromStart,&ChromEnd,szScore,&cStrand,&cFrame,&AttribStart);
	if(NumFields < 8)
		{
		AddErrMsg("CBEDfile::ProcessGFF23File","Errors whilst parsing - %s",m_szFile);
		Rslt = eBSFerrParse; // if at least one feature processed then treat rubbish as error
		break;
		}

	// only interested in a subset of feature types - those relating to genes and transcripts
	pFeatType = GFF3FeatureTypes;
	for(Idx=0; Idx < cNumGFF3FeatTypes; Idx++,pFeatType++)
		{
		if(!stricmp(szFeature,pFeatType->pszFeatType))
			break;
		}
	if(Idx == cNumGFF3FeatTypes)
		continue;

	// accepting this feature type
	FeatTypeID = pFeatType->FeatTypeID;

	// some fields may have '.' to indicate that that field value is unknown or not applicable
	// replace these instances with defaults

	// Source could be '.'; default these to be 'unknown'
	if(szSource[0] == '.')
		strcpy(szSource,"Unknown");

	// Score allowed to be '.' instead of a floating point if no score available; default 0
	if(szScore[0] == '.')
		iScore = 0;
	else
		{
		iScore = (int)atof(szScore);
		if(iScore < 0)
			iScore = 0;
		else
			if(iScore > 999)
				iScore = 999;
		}

	// Strand allowed to be '.' or '?'; if not '-' then default '+'
	if(cStrand != '-')
		cStrand = '+';
	// Frame allowed to be '.'; if '.' then default to be '0'
	switch(cFrame) {
		case '1':
			iFrame = 1;
			break;
		case '2':
			iFrame = 2;
			break;
		default:
			iFrame = 0;
			break;
		}

	// parse out attributes of interest and their associated values
	pTxt = &pTxt[AttribStart];
	AttribLen = (int)strlen(pTxt);

	int ParseState;		// 0 whilst looking for attribute; 1 looking for '='; 2 looking for associated value(s)
	int IdentLen;
	pAttribValue = AttribValues;
	for(Idx = 0; Idx < cNumGFF3AttribTypes; Idx++,pAttribValue++)
		{
		pAttribValue->AttribTypeID = Idx+1;
		pAttribValue->szAttribValue[0] = '\0';
		}

	ParseState = 0;
	for(Idx = 0; Idx < AttribLen; Idx++,pTxt++)
		{
		Chr = *pTxt;
		switch(Chr) {
			case '=':					
				if(ParseState != 1)				// if not in this parse state then reset state back to expecting to parse attribute
					ParseState = 0;
				else
					{
					ParseState = 2;				// now expecting to parse out the attribute value
					pAttribValue = &AttribValues[AttribTypeID-1];
					IdentLen = 0;
					}
				continue;

			case ' ':							// in any parse state just slough whitespace
			case '\t':
				continue;

			default:
				switch(ParseState) {
					case 0:							// trying to parse out attribute identifier
						pAttribType = GFF3AttribTypes;
						for(Idx=0; Idx < cNumGFF3AttribTypes; Idx++,pAttribType++)
							{
							if(!strcmp(pTxt,pAttribType->pszAttribType))
								break;
							}
						if(Idx == cNumGFF3FeatTypes)
							continue;
						IdentLen = (int)strlen(pAttribType->pszAttribType);
						AttribTypeID = pAttribType->AttribTypeID;
						Idx += IdentLen - 1; 
						pTxt += IdentLen - 1;
						ParseState = 1;
						continue;
					case 1:							// expecting '='
						ParseState = 0;				// so back to expecting attribute identifier
						continue;
					default:						// parsing out the attribute value
						if(Chr == ';' || Chr == '=')	// not allowed within attribute values
							{
							ParseState = 0;				// back to expecting attribute identifier
							continue;
							}
						break;
					}
				break;
			}

		if(IdentLen < cMaxGFF3AttribValueLen)
			{
			pAttribValue->szAttribValue[IdentLen++] = Chr;
			pAttribValue->szAttribValue[IdentLen] = '\0';
			}
		}
	// this is where this entry needs to be added to other entries
	}
// sort the entries and then iterate building BED entries

fclose(pGFFStream);
return(FeatNum ? Rslt : eBSFerrParse);
}

// ProcessBedFile
// Opens, parses (calling AddFeature()), closes specified file assumed to be in BED style format
// BED style format
// BED may start with some header lines, basic assumption is that these lines are not formated as <text>\t<int>\t<int>, other lines may start with '#' as comment lines
// Header and comment lines are simply sloughed, but it is expected that feature lines will be present within the first 20 lines parsed
//
// <chrom><tab><chromStart><tab><chromEnd><tab><name><tab><score><tab><strand>...
// chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or contig (e.g. ctgY1). 
// chromStart - The starting position of the feature in the chromosome or contig. The first base in a chromosome is numbered 0. 
// chromEnd - The ending position of the feature in the chromosome or contig. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 
// The 9 additional optional BED fields are:
// name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode. 
// score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). 
// strand - Defines the strand - either '+' or '-'. 
// thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). 
// thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays). 
// reserved - This should always be set to zero. 
// blockCount - The number of blocks (exons) in the BED line. 
// blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount. 
// blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount. 
teBSFrsltCodes 
CBEDfile::ProcessBedFile(char *pszFileName)	// file to process
{
FILE *pBEDStream;
int LineNum;
int NumOctParams;
char szChrom[(cMaxDatasetSpeciesChrom * 2) + 1]; // much longer allocation than usually required, parsed out chrom names will be truncated to be no longer than cMaxDatasetSpeciesChrom-1 
int chromStart;
int chromEnd;
char szName[(cMaxFeatNameLen * 2) + 1];  // much longer allocation than usually required, parsed out feature names will be truncated to be no longer than cMaxFeatNameLen-1
int SupInfoStart;
int Score;
char Strand;
int Cnt;
int FeatNum;
teBSFrsltCodes Rslt;
bool bCSV;
char *pTxt;

if(pszFileName == NULL || *pszFileName == '\0')
	return(eBSFerrParams);
if((pBEDStream = fopen(pszFileName,"r"))==NULL)
	{
	AddErrMsg("CBEDfile::ProcessBedFile","Unable to fopen BED format file %s error: %s",pszFileName,strerror(errno));
	return(eBSFerrOpnFile);
	}
LineNum = 0;
NumOctParams = 0;
FeatNum = 0;
bCSV = false;	// BED files by default used TAB as the field separator so initially assume TABs and not commas
Rslt = eBSFerrFileType; // assume the worst..., no features parsed
while(fgets(m_szLineBuff,sizeof(m_szLineBuff)-1,pBEDStream)!= NULL)
	{
	LineNum += 1;
	if(!FeatNum && LineNum >= 20)	// if can't load at least one feature from 1st 20 lines then can't be BED format
		{
		fclose(pBEDStream);
		return(eBSFerrFileType);
		}

	pTxt = TrimWhitespace(m_szLineBuff);
	if(*pTxt=='\0' || *pTxt=='#')	// simply slough lines which were just whitespace or start with '#'
		continue;

	m_szAttributes[0] = '\0';
	if(!bCSV)
		{
		// not CSV so assume it is a feature line to be processed
		Cnt = sscanf(pTxt," %140s %d %d %140s %d %c %n",
			szChrom,&chromStart,&chromEnd,szName,&Score,&Strand,&SupInfoStart);
		if(!FeatNum && Cnt < 3)	// if parse fails with TABs on 1st expected feature then explore ',' being used as the separator
			bCSV = true;
		}
	if(bCSV)	// failed on tabs as field separator, now expecting ',' as field separator
	    Cnt = sscanf(pTxt," %140s , %d , %d , %140s , %d , %c , %n",
			szChrom,&chromStart,&chromEnd,szName,&Score,&Strand,&SupInfoStart);
	if(Cnt < 3)		// must be rubbish or perhaps a header on this line
		{
		if(!FeatNum)
			{
			bCSV = false;	// could be a header line in BED file
			continue;		// no features processed, hope for the best and keep trying!
			}
		AddErrMsg("CBEDfile::ProcessBedFile","Errors whilst parsing - %s",m_szFile);
		Rslt = eBSFerrParse; // if at least one feature processed then treat rubbish as error
		break;
		}
		
	szChrom[cMaxDatasetSpeciesChrom-1] = '\0'; // have seen some chrom and feature names with exactly cMaxDatasetSpeciesChrom or cMaxFeatNameLen chars excluding the terminator
	szName[cMaxFeatNameLen-1] ='\0';	
	
	switch(m_FileHdr.FeatType) {
		case eBTGeneExons:		// check if line parsed does have supplementary info covering exon starts etc, if not then dummy up one with single exon
			switch(Cnt) {
				case 3:		// name, score,strand,suppinfo missing 
					sprintf(szName,"feat%d",FeatNum+1);
					Score = 0; Strand = '+'; 
					break;
				case 4:		// score,strand,suppinfo missing 
					Score = 0; Strand = '+'; 
					break;
				case 5:		// strand,suppinfo missing
					Strand = '+'; 
					break;
				}
			if(Cnt < 6 || pTxt[SupInfoStart] == '\0' ||	// if no exon starts etc then dummy up...
				pTxt[SupInfoStart] == '\n' || pTxt[SupInfoStart] == '\r')
				{
				sprintf(m_szAttributes,"%d\t%d\t0\t1\t%d,\t0,",chromStart,chromEnd,chromEnd-chromStart);
				}
			else
				{
				strncpy(m_szAttributes,&pTxt[SupInfoStart],sizeof(m_szAttributes));
				m_szAttributes[sizeof(m_szAttributes)-1] = '\0';
				}
			break;

		case eBTGeneral:		// may, or may not, have supplementary info
		default:
			switch(Cnt) {
				case 3:		// name, score,strand,suppinfo missing 
					sprintf(szName,"feat%d",FeatNum+1);
					Score = 0; Strand = '+'; m_szAttributes[0] = '\0';
					break;
				case 4:		// score,strand,suppinfo missing 
					Score = 0; Strand = '+'; m_szAttributes[0] = '\0';
					break;
				case 5:		// strand,suppinfo missing
					Strand = '+'; m_szAttributes[0] = '\0';
					break;
				case 6:		// suppinfo may be present, 
					if(pTxt[SupInfoStart] == '\0' || 
						pTxt[SupInfoStart] == '\n' || pTxt[SupInfoStart] == '\r')
						m_szAttributes[0] = '\0';
					else
						{
						strncpy(m_szAttributes,&pTxt[SupInfoStart],sizeof(m_szAttributes));
						m_szAttributes[sizeof(m_szAttributes)-1] = '\0';
						m_FileHdr.FeatType = eBTGeneExons;
						}
					break;
				default:
					continue;
				}
			break;
		}
	
	if(Strand == '.')		// seen cases where strand in BED was '.' so assume these were referencing the sense strand
		Strand = '+';
	else
		{
		if(!(Strand == '+' || Strand == '-' || Strand == '?'))	// treat any other non-'+' or '-' as being that strand is unknown
			Strand = '?';
		}

	if((Rslt=AddFeature(szName,szChrom,chromStart,chromEnd,Score,Strand,m_szAttributes))!=eBSFSuccess)
		{
		fclose(pBEDStream);
		return(Rslt);
		}
	FeatNum++;
	if(m_MaxParseFeats > 0 && FeatNum >= m_MaxParseFeats)
		break;
	}
fclose(pBEDStream);
return(FeatNum ? Rslt : eBSFerrParse);
}


// ProcessUltraCoreFile
// Opens, parses (calling AddFeature()), closes specified file assumed to be in UltraCoreFile style format
// UltraCoreFile style format
// <CoreID>,<RefSpecies>,<RelSpecies>,<TotLen>,<LeftFlankLen>,<CoreLen>,<RightFlankLen>,<RefChrom>,<RefChromOfs>,<RelChrom>,<RelChromOfs>,<Strand>,<Sequence>
// CoreID - uniquely identifies core in this file
// RefSpecies - name of reference species e.g "hg17"
// RelSpecies - name of relative species e.g "tetNig1" 
// TotLen - total length Left+Core+right e.g 150
// LeftFlankLen - size of left flanking region containing missmatches 
// CoreLen - size of core with no missmatches
// RightFlankLen - size of right flanking region containing missmatches
// RefChrom - reference species chromosome e.g "chr10"
// RefChromOfs - offset on RefChrom at which core starts e.g 123456789
// RelChrom - relative chromosome e.g "chr13"
// RelChromOfs - offset on RelChrom at which core starts e.g 123456789
// Strand - relative species strand "+" or "-"
// Sequence - the actual core sequence
teBSFrsltCodes 
CBEDfile::ProcessUltraCoreFile(char *pszFileName)	// file to process
{
FILE *pBEDStream;
char *pChr;
char *pTxt;
char Chr;
int LineNum;
int TotLen;
int LeftFlankLen;
int CoreLen;
int RightFlankLen;
char szLineBuff[cLineBuffLen];
char szRefSpecies[51];
char szRefChrom[51];
int RefChromOfs;
char szRelSpecies[51];
char szRelChrom[51];
int RelChromOfs;
char szName[cMaxFeatNameLen+1];
char Strand;
int Cnt;
int FeatNum;
teBSFrsltCodes Rslt;

if(pszFileName == NULL || *pszFileName == '\0')
	return(eBSFerrParams);
if((pBEDStream = fopen(pszFileName,"r"))==NULL)
	{
	AddErrMsg("CBEDfile::ProcessBedFile","Unable to open UltraCoreFile format file %s error: %s",pszFileName,strerror(errno));
	return(eBSFerrOpnFile);
	}
LineNum = 0;
FeatNum = 1;
while(fgets(szLineBuff,sizeof(szLineBuff),pBEDStream)!= NULL)
	{
	pTxt = TrimWhitespace(szLineBuff);
	if(*pTxt=='\0' || *pTxt=='#')	// simply slough lines which were just whitespace or comment lines
		continue;

	// replace all instances of \" with ' ', makes it easier to parse
	pChr = pTxt;
	while(Chr = *pChr++)
		{
		if(Chr == '\"')
			pChr[-1] = ' ';
		}
    Cnt = sscanf(pTxt," %d , %50[^ ,] , %50[^ ,] , %d , %d , %d , %d, %50[^ ,] , %d , %50[^ ,] , %d , %c ,",
			&FeatNum,szRefSpecies,szRelSpecies,&TotLen,&LeftFlankLen,&CoreLen,&RightFlankLen,szRefChrom,&RefChromOfs,szRelChrom,&RelChromOfs,&Strand);
	if(Cnt != 12)	// slough lines which do not contain expected format
		continue;
	sprintf(szName,"core%d",FeatNum);
	if((Rslt=AddFeature(szName,szRefChrom,RefChromOfs,RefChromOfs+CoreLen,TotLen,Strand,NULL))!=eBSFSuccess)
		{
		fclose(pBEDStream);
		return(Rslt);
		}
	}
fclose(pBEDStream);
return(eBSFSuccess);
}

// ProcessGroupCoreFile
// Opens, parses (calling AddFeature()), closes specified file assumed to be in group Core style format
// group core style format
// <CoreID>,<Group>,<RefSpecies>,<RefChrom>,<RefChromOfs>,<CoreLen>
// CoreID - uniquely identifies core in this file (1..n)
// Group - core is common to all species within this group e.g "mammals"
// RefSpecies - name of reference species e.g "hg17"
// RefChrom - reference species chromosome e.g "chr10"
// RefStartOfs - offset on RefChrom at which core starts e.g 123456789
// CoreLen - size of core with no missmatches
// Sequence - the actual core sequence
teBSFrsltCodes 
CBEDfile::ProcessGroupCoreFile(char *pszFileName,	// file to process
							   int MinLen,			// length range accepted 
							   int MaxLen)	
{
FILE *pBEDStream;
char *pChr;
char *pTxt;
char Chr;
int LineNum;
int CoreLen;
int RefChromEndOfs;
char szGroup[51];
char szLineBuff[cLineBuffLen];
char szRefSpecies[51];
char szRefChrom[51];
int RefChromOfs;
char szName[cMaxFeatNameLen+1];
int Cnt;
int FeatNum;
teBSFrsltCodes Rslt;
bool bFormatKnown = false;
bool bFormatExtd = false;

if(pszFileName == NULL || *pszFileName == '\0')
	return(eBSFerrParams);
if((pBEDStream = fopen(pszFileName,"r"))==NULL)
	{
	AddErrMsg("CBEDfile::ProcessBedFile","Unable to open Group UltraCoreFile format file %s error: %s",pszFileName,strerror(errno));
	return(eBSFerrOpnFile);
	}
LineNum = 0;
FeatNum = 1;
while(fgets(szLineBuff,sizeof(szLineBuff),pBEDStream)!= NULL)
	{
	pTxt = TrimWhitespace(szLineBuff);
	if(*pTxt=='\0')	// simply slough lines which were just whitespace
		continue;
	// replace all instances of \" with ' ', makes it easier to parse
	pChr = pTxt;
	while(Chr = *pChr++)
		{
		if(Chr == '\"')
			pChr[-1] = ' ';
		}
	// first try the extended format in which there is an additional field containing the end offset
    if(!bFormatKnown || bFormatExtd)
		{
		Cnt = sscanf(pTxt," %d , %50[^ ,] , %50[^ ,] ,  %50[^ ,] , %d , %d , %d",
			&FeatNum,szGroup,szRefSpecies,szRefChrom,&RefChromOfs,&RefChromEndOfs,&CoreLen);
		if(Cnt == 7)
			bFormatKnown = bFormatExtd = true;
		else
			Cnt = 0;
		}
	if(!bFormatKnown || !bFormatExtd)
		{
		Cnt = sscanf(pTxt," %d , %50[^ ,] , %50[^ ,] ,  %50[^ ,] ,%d , %d",
			&FeatNum,szGroup,szRefSpecies,szRefChrom,&RefChromOfs,&CoreLen);
		if(Cnt == 6)
			{
			bFormatKnown = true;
			bFormatExtd = false;
			}
		else
			Cnt = 0;
		}
	if(Cnt == 0)
		continue;

	// slough cores whose length is outside of specified range
	if(CoreLen < MinLen || CoreLen > MaxLen)
		continue;
	sprintf(szName,"gcore%d",FeatNum);
	if((Rslt=AddFeature(szName,szRefChrom,RefChromOfs,RefChromOfs+CoreLen,CoreLen,'+',NULL))!=eBSFSuccess)
		{
		fclose(pBEDStream);
		return(Rslt);
		}
	}
fclose(pBEDStream);
return(eBSFSuccess);
}


// SortFeatures
// Sorts feature names and start offsets so that searches on these can be optimised with
// binary rather than linear searches
teBSFrsltCodes
CBEDfile::SortFeatures(void)
{
int Idx;
int FeatLen;
int FeatureID;
int ChromID;
int NameInst;
tsBEDfeature *pFeature;
tsBEDfeature *pPrevFeature;
int NumThreads;

tsBEDchromname *pChrom;

if(!m_FileHdr.NumFeatures || m_pFeatures == NULL)
	return(eBSFerrNoFeatures);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sort optimising loaded chromosomes (%d) and features (%d) ...",m_FileHdr.NumChroms,m_FileHdr.NumFeatures);

if(m_ppFeatureNames != NULL)
	{
	delete m_ppFeatureNames;
	m_ppFeatureNames = NULL;
	}
if(m_ppFeatureChromStarts != NULL)
	{
	delete m_ppFeatureChromStarts;
	m_ppFeatureChromStarts = NULL;
	}


#ifdef _WIN32
SYSTEM_INFO SystemInfo;
GetSystemInfo(&SystemInfo);
NumThreads = SystemInfo.dwNumberOfProcessors;
#else
NumThreads = sysconf(_SC_NPROCESSORS_CONF);
#endif
NumThreads = min(8,NumThreads);
m_mtqsort.SetMaxThreads(NumThreads);

// sort chromosome names into ascending order
m_mtqsort.qsort(m_pChromNames,m_FileHdr.NumChroms,sizeof(tsBEDchromname),SortChromNames);


tsUnsortToSortID *pAllocdU2S;
tsUnsortToSortID *pU2S;

if((pAllocdU2S = new tsUnsortToSortID [m_FileHdr.NumChroms])==NULL)
	return(eBSFerrMem);

// remap ChromIDs in chromnames to be Idx+1
// at same time initialise fields used later
pChrom = m_pChromNames;
pU2S = pAllocdU2S;
for(ChromID=1;ChromID <= m_FileHdr.NumChroms;ChromID++,pChrom++,pU2S++)
	{
	pU2S->OldID = pChrom->ChromID;
	pU2S->NewID = ChromID;
	pChrom->ChromID = ChromID;
	pChrom->FirstStartID = 0;
	pChrom->LastStartID = 0;
	pChrom->NumFeatures = 0;
	pChrom->MaxFeatLen = 0;
	}
if(m_FileHdr.NumChroms > 1)
	m_mtqsort.qsort(pAllocdU2S,m_FileHdr.NumChroms,sizeof(tsUnsortToSortID),SortU2S);

pFeature = m_pFeatures;
for(Idx = 0; Idx < m_FileHdr.NumFeatures; Idx++)
	{
	pU2S = LocateU2S(pFeature->ChromID,m_FileHdr.NumChroms,pAllocdU2S);
	pFeature->ChromID = pU2S->NewID;
	pFeature = (tsBEDfeature *)(((char *)pFeature) + pFeature->Size);
	}

delete pAllocdU2S;

// create array of ptrs (m_ppFeatureNames) into concatenated features which will be used to sort
// by name->chrom->start->end
m_ppFeatureNames = new tsBEDfeature * [m_FileHdr.NumFeatures];
if(m_ppFeatureNames == NULL)
	return(eBSFerrMem);

pFeature = m_pFeatures;
for(FeatureID = 0; FeatureID < m_FileHdr.NumFeatures; FeatureID++)
	{
	m_ppFeatureNames[FeatureID] = pFeature;
	pFeature = (tsBEDfeature *)(((char *)pFeature) + pFeature->Size);
	}
m_mtqsort.qsort(m_ppFeatureNames,m_FileHdr.NumFeatures,sizeof(tsBEDfeature *),SortFeatureNames);

// now that they are sorted then determine the number of name instances
NameInst = 1;
pPrevFeature = m_ppFeatureNames[0];
pPrevFeature->NameInst = NameInst++;
for(FeatureID = 1; FeatureID < m_FileHdr.NumFeatures; FeatureID++)
	{
	pFeature = m_ppFeatureNames[FeatureID];
	if(pPrevFeature->Hash != pFeature->Hash ||
#ifdef _WIN32
		stricmp(pFeature->szName,pPrevFeature->szName))
#else
		strcasecmp(pFeature->szName,pPrevFeature->szName))
#endif

		NameInst = 1;
	pFeature->NameInst = NameInst++;
	pPrevFeature = pFeature;
	}

// create array of ptrs (m_ppFeatureChromStarts) into concatenated features which will be used to sort
// by chrom->start->end->nameinst
m_ppFeatureChromStarts = new tsBEDfeature * [m_FileHdr.NumFeatures];
if(m_ppFeatureChromStarts == NULL)
	{
	delete m_ppFeatureNames;
	m_ppFeatureNames = NULL;
	return(eBSFerrMem);
	}

pFeature = m_pFeatures;
for(FeatureID = 0; FeatureID < m_FileHdr.NumFeatures; FeatureID++)
	{
	m_ppFeatureChromStarts[FeatureID] = pFeature;
	pFeature = (tsBEDfeature *)(((char *)pFeature) + pFeature->Size);
	}
m_mtqsort.qsort(m_ppFeatureChromStarts,m_FileHdr.NumFeatures,sizeof(tsBEDfeature *),SortChromStarts);


// make the feature ids map to the chrom start indexes, determine per chromosome maximum feature lengths/number of features
for(Idx = 0; Idx < m_FileHdr.NumFeatures; Idx++)
	{
	pFeature = m_ppFeatureChromStarts[Idx];
	pFeature->FeatureID = Idx+1;
	pChrom = &m_pChromNames[pFeature->ChromID-1];
	if(pChrom->FirstStartID == 0)
		{
		pChrom->FirstStartID = Idx+1;
		pChrom->NumFeatures = 1;
		}
	else
		pChrom->NumFeatures += 1;
	pChrom->LastStartID = Idx+1;
	FeatLen = pFeature->End - pFeature->Start + 1;
	if(FeatLen > pChrom->MaxFeatLen)
		pChrom->MaxFeatLen = FeatLen;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sort optimisation completed");
return(eBSFSuccess);
}

bool 
CBEDfile::SetStrand(char Strand)	// sets globally which strand features must be on in subsequent processing '+'/'-' or '*' for either
{
if(Strand != '*' && Strand != '+' && Strand != '-')
	return(false);
m_OnStrand = Strand;
return(true);
}


tsUnsortToSortID *
CBEDfile::LocateU2S(UINT32 OldID,	// original feature chromid to locate
		  int NumOldIDs,			// number of old chromids to binary search over
		  tsUnsortToSortID *pSortedU2s) // sorted old to new chrom id mappings
{
tsUnsortToSortID *pEl;
INT32 StartIdx;
INT32 MidIdx;
INT32 EndIdx;
if(NumOldIDs < 25)				// if just a few then do a simple linear search
	{
	while(NumOldIDs--)
		{
		if(pSortedU2s->OldID == OldID)
			return(pSortedU2s);
		pSortedU2s += 1;
		}
	return(NULL);
	}

StartIdx = 0;
EndIdx = NumOldIDs-1;
while(EndIdx >= StartIdx) {
	MidIdx = (EndIdx + StartIdx)/2;
	pEl = &pSortedU2s[MidIdx];
	if(pEl->OldID == OldID)
		return(pEl);
	if(pEl->OldID > OldID)	
		EndIdx = MidIdx - 1;
	else
		StartIdx = MidIdx + 1;
	}
return(NULL);
}

// LocateChromName
// checks if chromosome was previously cached as result of successful LocateChromName() and
// if not a cache match then a binary lookup
tsBEDchromname *
CBEDfile::LocateChromName(char *pszChromName)
{
int Rslt;
int StartIdx;
int EndIdx;
int MidIdx;
int Cnt;
tsBEDchromname *pChromName;
UINT32 Hash;
if(pszChromName == NULL || *pszChromName == '\0' ||
   m_FileHdr.NumChroms == 0 || m_pChromNames == NULL)
	return(NULL);

Hash = (UINT32)CUtility::GenHash24(pszChromName);

// check if name matches that previously returned
if(m_pCacheChromName != NULL)
	{
	if(Hash == m_pCacheChromName->Hash &&
			!stricmp(m_pCacheChromName->szName,pszChromName))
		return(m_pCacheChromName);
	}

// name not in cache
// If less than 25
if(m_FileHdr.NumChroms <= 25)
	{
	pChromName = m_pChromNames;
	for(Cnt = 0; Cnt < m_FileHdr.NumChroms; Cnt++,pChromName++)
		{
		if(Hash == pChromName->Hash &&
			!stricmp(pChromName->szName,pszChromName))
			return(m_pCacheChromName = pChromName);	// cache the successful locate
		}
	return(NULL);	// no chromosome with specified name located
	}

// more than 25, do a binary search
StartIdx = 0;
EndIdx = m_FileHdr.NumChroms-1;

while(EndIdx >= StartIdx) {
	MidIdx = (EndIdx + StartIdx)/2;
	pChromName = &m_pChromNames[MidIdx];
	Rslt = stricmp(pszChromName,pChromName->szName);
	if(Rslt < 0)	
		{
		EndIdx = MidIdx - 1;
		continue;
		}

	if(Rslt > 0)	
		{
		StartIdx = MidIdx + 1;
		continue;
		}

	return(pChromName);
	}
return(NULL);
}

// LocateChromName
// Locates chromosome name by simply using ChromID as the index into m_pChromNames
tsBEDchromname *
CBEDfile::LocateChromName(int ChromID)
{
if(ChromID < 1 || ChromID > m_FileHdr.NumChroms ||
   m_FileHdr.NumChroms == 0 || m_pChromNames == NULL)
	return(NULL);
return(&m_pChromNames[ChromID - 1]);
}

// AddFeature
// Adds new feature
// If first time chromosome name referenced then adds new chromosome
teBSFrsltCodes
CBEDfile::AddFeature(char *pszFeatName,			// feature name
					 char *pszChromName,		// chromosome name
					 int Start,					// start offset (0..n) on chromosome
					 int End,					// end offset on chromosome - note: will subtract 1 to make it an inclusive end
					 int Score,					// feature score
					 char Strand,				// on which strand feature is located
					 char *pszSuppInfo)			// any supplementary information
{
int ChromNameLen;
int FeatNameLen;
int SuppInfoLen;
int thickStart;			// The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). 
int thickEnd;			// The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays). 
int	Reserved;			// This should always be set to zero. 
int blockCount;			// The number of blocks (exons) in the BED line.
int ExonStarts[cMaxNumExons];	// to hold exon starts as parsed from BED file in pszSuppInfo
int ExonSizes[cMaxNumExons]; //to hold exon sizes as parsed from BED file in pszSuppInfo
int Psn,RelPsn;
int Cnt,Idx;
UINT32 Hash;
static 	int MaxLinkedLen = 0;

static int MaxNumFeatures = 0;

tsGeneStructure GeneStructure;

tsBEDchromname *pChromName;
tsBEDfeature *pFeature;
if(pszFeatName == NULL || pszFeatName[0] == '\0' ||
	pszChromName == NULL || pszChromName[0] == '\0' ||
	Start < 0 || End < Start || (Strand != '+' && Strand != '-' && Strand != '?'))
	return(eBSFerrParams);

if((ChromNameLen = (int)strlen(pszChromName)) > cMaxDatasetSpeciesChrom - 1)
	return(eBSFerrParams);
if((FeatNameLen = (int)strlen(pszFeatName)) > cMaxFeatNameLen - 1)
	return(eBSFerrParams);
if(pszSuppInfo == NULL || *pszSuppInfo == '\0')
	SuppInfoLen = 0;
else
	SuppInfoLen = (int)strlen(pszSuppInfo) + 1;

if(m_FileHdr.NumFeatures == m_FileHdr.MaxFeatures)
	return(eBSFerrMaxFeatures);

// realloc memory as may be required to hold at least one new  tsBEDchromname instance in m_pChromNames
if(m_pChromNames == NULL || 
   (((m_FileHdr.ChromNamesSize + (UINT32)sizeof(tsBEDchromname)) >= (UINT32)m_AllocChromNamesSize)))
	{
	size_t ReallocTo;
	if(m_pChromNames == NULL)
		{
		ReallocTo = cAllocBEDChromNamesIncr;
#ifdef _WIN32
		m_pChromNames = (tsBEDchromname *) malloc(ReallocTo);	// initial and perhaps the only allocation
#else
			// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
		m_pChromNames = (tsBEDchromname *)mmap(NULL,ReallocTo, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
		if(m_pChromNames == MAP_FAILED)
			m_pChromNames = NULL;
#endif
		if(m_pChromNames == NULL)
			{
			AddErrMsg("CBEDfile::AddFeature","Memory allocation of %lld bytes failed - %s",(INT64)ReallocTo,strerror(errno));
			m_AllocChromNamesSize = 0;
			return(eBSFerrMem);
			}
		m_AllocChromNamesSize = ReallocTo;
		m_FileHdr.ChromNamesSize = 0;
	    memset(m_pChromNames,0,sizeof(tsBEDchromname));
		}
	else
		{
		ReallocTo = m_AllocChromNamesSize + cAllocBEDChromNamesIncr;
#ifdef _WIN32
		pChromName = (tsBEDchromname *) realloc(m_pChromNames,ReallocTo);
#else
		pChromName = (tsBEDchromname *)mremap(m_pChromNames,m_AllocChromNamesSize,ReallocTo,MREMAP_MAYMOVE);
		if(pChromName == MAP_FAILED)
			pChromName = NULL;
#endif
		if(pChromName == NULL)
			return(eBSFerrMem);
		m_pChromNames = pChromName;
		m_AllocChromNamesSize = ReallocTo;
		}
	m_pCacheChromName = NULL;
	}


switch(m_FileHdr.FeatType) {
	case eBTGeneExons:				// need to parse supplementary information into a tsGeneStructure
		if(pszSuppInfo != NULL && pszSuppInfo[0] != '\0')
			Cnt = sscanf(pszSuppInfo," %d %d %d %d %n",&thickStart,&thickEnd,&Reserved,&blockCount,&Psn);
		else
			Cnt = 0;	
		if(Cnt == 0)	// if unable to get supplementary info because there was none then use the gene start and end as representing a single exon
			{
			blockCount = 1;
			GeneStructure.NumExons = 1;
			thickStart = Start;
			thickEnd = End;
			GeneStructure.ExonStartEnds[0] = 0;
			GeneStructure.ExonStartEnds[1] = End - Start - 1;
			}
		else
			if(Cnt != 4 || blockCount > cMaxNumExons)
				return(eBSFerrFeature);

		if(blockCount > MaxNumFeatures)
			MaxNumFeatures = blockCount;

		if(Cnt == 4)
			{
			for(Idx = 0; Idx < blockCount; Idx++)
				{
				sscanf(&pszSuppInfo[Psn]," %d , %n",&ExonSizes[Idx],&RelPsn);
				Psn += RelPsn;
				}
			for(Idx = 0; Idx < blockCount-1; Idx++)
				{
				sscanf(&pszSuppInfo[Psn]," %d , %n",&ExonStarts[Idx],&RelPsn);
				Psn += RelPsn;
				}
			sscanf(&pszSuppInfo[Psn]," %d",&ExonStarts[Idx]);
			for(Idx = 0; Idx < blockCount; Idx++)
				{
				GeneStructure.ExonStartEnds[Idx*2] = ExonStarts[Idx];
				GeneStructure.ExonStartEnds[(Idx*2)+1] = ExonStarts[Idx] + ExonSizes[Idx] - 1;
				}
			}

		GeneStructure.NumExons = blockCount;
		GeneStructure.thickStart = thickStart - Start;
		GeneStructure.thickEnd   = thickEnd - Start;
		GeneStructure.Size = offsetof(tsGeneStructure,ExonStartEnds) + (blockCount * sizeof(int) * 2);
        SuppInfoLen = GeneStructure.Size;
		break;
	default:
		break;
	}

Hash = (UINT32)CUtility::GenHash24(pszChromName);
pChromName = NULL;

// check if name matches prev cached chrom name
if(m_pCacheChromName != NULL)
	{
	if(Hash == m_pCacheChromName->Hash &&
			!stricmp(m_pCacheChromName->szName,pszChromName))
			pChromName = m_pCacheChromName;
	}

if(pChromName == NULL) // NULL if not prev cached
	{
	int LinkedLen;

	bool bNewChrom = false;
	pChromName = m_pChromNames;
	if(m_FileHdr.NumChroms == 0)		// first so can't be any others with same hash..
		{
		memset(m_pChromHashes,0,sizeof(UINT32) * 0x1000000); // ensure hashes are reset - 0 if no chrom has that hash
		bNewChrom = true;
		}
	else
		{
		LinkedLen = 0;
		if(m_pChromHashes[Hash-1] != 0)		// any chrom same hash?
			{
			pChromName = &m_pChromNames[m_pChromHashes[Hash-1]-1];
			while(1) 
				{
				if(Hash == pChromName->Hash && !stricmp(pszChromName,pChromName->szName))
					break;			
				if(pChromName->FirstStartID == 0)
					{
					bNewChrom = true;
					break;
					}
				pChromName = &m_pChromNames[pChromName->FirstStartID - 1];
				LinkedLen += 1;
				}
			if(LinkedLen > MaxLinkedLen)
				MaxLinkedLen = LinkedLen;
			}
		else
			bNewChrom = true;
		}

	if(bNewChrom)
		{
		if(m_pChromHashes[Hash-1] == 0)	
			m_pChromHashes[Hash-1] = m_FileHdr.NumChroms + 1;	
		else
			pChromName->FirstStartID = m_FileHdr.NumChroms + 1;
		pChromName = &m_pChromNames[m_FileHdr.NumChroms++];
		strcpy(pChromName->szName,pszChromName);
		pChromName->ChromID = m_FileHdr.NumChroms;
		pChromName->NumFeatures = 0;
		pChromName->FirstStartID = 0;
		pChromName->LastStartID = 0;
		pChromName->MaxFeatLen = 0;
		pChromName->Hash = Hash;
		m_FileHdr.ChromNamesSize += sizeof(tsBEDchromname);
		}
	m_pCacheChromName = pChromName;
	}

pChromName->NumFeatures++;

// realloc memory as may be required to hold a new tsBEDfeature instance
if(m_pFeatures == NULL || 
   ((m_FileHdr.FeaturesSize + (int)sizeof(tsBEDfeature) + FeatNameLen + SuppInfoLen+100) > m_AllocFeaturesSize))
	{
	size_t ReallocTo;
	if(m_pFeatures == NULL)
		{
		ReallocTo = cAllocFeatIncr*10;
#ifdef _WIN32
		m_pFeatures = (tsBEDfeature *) malloc(ReallocTo);	// initial and perhaps the only allocation
#else
			// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
		m_pFeatures = (tsBEDfeature *)mmap(NULL,ReallocTo, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
		if(m_pFeatures == MAP_FAILED)
			m_pFeatures = NULL;
#endif

		if(m_pFeatures == NULL)
			{
			AddErrMsg("CBEDfile::LoadFeatures","Memory allocation of %lld bytes failed - %s",(INT64)ReallocTo,strerror(errno));
			m_AllocFeaturesSize = 0;
			return(eBSFerrMem);
			}
		m_AllocFeaturesSize = ReallocTo;
		}
	else
		{
		ReallocTo = m_FileHdr.FeaturesSize + cAllocFeatIncr;
#ifdef _WIN32
		pFeature = (tsBEDfeature *) realloc(m_pFeatures,ReallocTo);
#else
		pFeature = (tsBEDfeature *)mremap(m_pFeatures,m_AllocFeaturesSize,ReallocTo,MREMAP_MAYMOVE);
		if(pFeature == MAP_FAILED)
			pFeature = NULL;
#endif
		if(pFeature == NULL)
			return(eBSFerrMem);
		m_pFeatures = pFeature;
		m_AllocFeaturesSize = ReallocTo;
		}
	}
pFeature = (tsBEDfeature *)((char *)m_pFeatures + m_FileHdr.FeaturesSize);
pFeature->FeatureID = ++m_FileHdr.NumFeatures;
pFeature->Size = sizeof(tsBEDfeature) + FeatNameLen + SuppInfoLen;
m_FileHdr.FeaturesSize += pFeature->Size;
pFeature->FiltFlags = cFeatFiltIn;		// by default all features are locatable when searching by chromosome offset
pFeature->UserClass = 0;
pFeature->Strand = Strand;
pFeature->Score = Score;
pFeature->ChromID = pChromName->ChromID;
pFeature->Start = Start;
pFeature->End = End - 1;

pFeature->FeatNameLen = FeatNameLen;
strcpy(pFeature->szName,pszFeatName);	// supplementary info is concatenated onto the feature name
if(SuppInfoLen)
	{
	switch(m_FileHdr.FeatType) {
		case eBTGeneExons:				// need to parse supplementary information into a tsGeneStructure
			memmove(&pFeature->szName[FeatNameLen+1],&GeneStructure,SuppInfoLen);
			break;
		default:
			strcpy(&pFeature->szName[FeatNameLen+1],pszSuppInfo);
			break;
		}
	}
pFeature->Hash = (UINT32)CUtility::GenHash24(pszFeatName);
return(eBSFSuccess);
}

teBSFrsltCodes
CBEDfile::Flush2Disk(void)
{
int WrtLen;
teBSFrsltCodes Rslt;
int *pOfs;
int Idx;
int Idy;

if(m_hFile != -1 && m_bCreate)		// if file opened for write
	{
	if(m_FileHdr.NumChroms && m_FileHdr.NumFeatures && m_FileHdr.FeaturesSize && m_FileHdr.ChromNamesSize)
		{
		if((Rslt=SortFeatures())!=eBSFSuccess)
			{
			AddErrMsg("CBEDfile::Flush2Disk","Unable to sort features %s",m_szFile);
			Reset(false);
			return(Rslt);
			}

		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen)
			{
			AddErrMsg("CBEDfile::Flush2Disk","Unable to seek chromosome names to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}

		if(m_bIsBigEndian)
			{
			tsBEDchromname *pChromName = m_pChromNames;
			for(Idx = 0; Idx < m_FileHdr.NumChroms; Idx++, pChromName++)
				{
				pChromName->ChromID = SwapUI32Endians(pChromName->ChromID);
				pChromName->FirstStartID = SwapUI32Endians(pChromName->FirstStartID);
				pChromName->Hash = SwapUI32Endians(pChromName->Hash);
				pChromName->LastStartID = SwapUI32Endians(pChromName->LastStartID);
				pChromName->MaxFeatLen = SwapUI32Endians(pChromName->MaxFeatLen);
				pChromName->NumFeatures=SwapUI32Endians(pChromName->NumFeatures);
				}
			}
		WrtLen = m_FileHdr.ChromNamesSize;	
		if(write(m_hFile,m_pChromNames,WrtLen)!=WrtLen)
			{
			AddErrMsg("CBEDfile::Flush2Disk","Unable to write chromosome names to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}

		m_FileHdr.ChromNamesOfs = m_FileHdr.FileLen;
		m_FileHdr.FileLen += m_FileHdr.ChromNamesSize;

		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen)
			{
			AddErrMsg("CBEDfile::Flush2Disk","Unable to seek features to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}

		if(m_bIsBigEndian)
			{
			tsBEDfeature *pFeature;
			tsGeneStructure *pGene;
			UINT8 *pBytes;

			pFeature = m_pFeatures;
			pBytes = (UINT8 *)pFeature;
			for(Idx = 0; Idx < m_FileHdr.NumFeatures; Idx++)
				{
				if(m_FileHdr.FeatType == eBTGeneExons)
					pGene = (tsGeneStructure *)&pFeature->szName[pFeature->FeatNameLen+1];
				else
					pGene = NULL;
				pBytes += pFeature->Size;
				pFeature->FeatureID = SwapUI32Endians(pFeature->FeatureID);					// uniquely identifies this feature
				pFeature->Size = SwapUI32Endians(pFeature->Size);						// size in bytes of this instance when concatenated
				pFeature->FiltFlags = SwapUI32Endians(pFeature->FiltFlags);					// used when filtering features in/out from feature bit processing
				pFeature->NameInst = SwapUI32Endians(pFeature->NameInst);					// name instance (1..n)
				pFeature->ChromID = SwapUI32Endians(pFeature->ChromID);					// feature is on this chromosome
				pFeature->Start = SwapUI32Endians(pFeature->Start);						// feature starts at this chromosome offset (0..n)				
				pFeature->End = SwapUI32Endians(pFeature->End);						// feature ends at this chromosome offset (Start + Len -1)
				pFeature->Score = SwapUI32Endians(pFeature->Score);						// feature score
				pFeature->Hash = SwapUI32Endians(pFeature->Hash);						// hash on szName
				if(pGene != NULL)
					{
					INT32 *pStartEnd;
					pStartEnd = pGene->ExonStartEnds;
					for(Idy = 0; Idy < pGene->NumExons*2; Idy++,pStartEnd++)
						*pStartEnd = SwapUI32Endians(*pStartEnd);
					pGene->Size  = SwapUI32Endians(pGene->Size);							// size of this instance
					pGene->NumExons = SwapUI32Endians(pGene->NumExons);						// number of exons in this gene
					pGene->thickStart = SwapUI32Endians(pGene->thickStart);						// where coding starts 
					pGene->thickEnd = SwapUI32Endians(pGene->thickEnd);						// where coding ends
					}
				pFeature = (tsBEDfeature *)pBytes;
				}
			}

		WrtLen = m_FileHdr.FeaturesSize;	
		if(write(m_hFile,m_pFeatures,WrtLen)!=WrtLen)
			{
			AddErrMsg("CBEDfile::Flush2Disk","Unable to write features to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}

		m_FileHdr.FeaturesOfs = m_FileHdr.FileLen;
		m_FileHdr.FileLen += WrtLen;

		// change index ptrs into 32 bit offsets relative to m_pFeatures
		pOfs = (int *)m_ppFeatureNames;
		for(Idx = 0; Idx < m_FileHdr.NumFeatures; Idx++, pOfs++)
			{
			*pOfs = (int)((char *)m_ppFeatureNames[Idx] - (char *)m_pFeatures);
			if(m_bIsBigEndian)
				*pOfs = SwapUI32Endians(*pOfs);
			}

		WrtLen = m_FileHdr.NumFeatures * sizeof(int);	
		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
			write(m_hFile,m_ppFeatureNames,WrtLen)!=WrtLen)
			{
			AddErrMsg("CBEDfile::Flush2Disk","Unable to write feature name index to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		m_FileHdr.FeatureNamesOfs = m_FileHdr.FileLen;
		m_FileHdr.FileLen += WrtLen;



		pOfs = (int *)m_ppFeatureChromStarts;
		for(Idx = 0; Idx < m_FileHdr.NumFeatures; Idx++,pOfs++)
			{
			*pOfs = (int)((char *)m_ppFeatureChromStarts[Idx] - (char *)m_pFeatures);
			if(m_bIsBigEndian)
				*pOfs = SwapUI32Endians(*pOfs);
			}
			
		WrtLen = m_FileHdr.NumFeatures * sizeof(int);	
		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
			write(m_hFile,m_ppFeatureChromStarts,WrtLen)!=WrtLen)
			{
			AddErrMsg("CBEDfile::Flush2Disk","Unable to write feature chrom->start index to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		m_FileHdr.FeatureChromStartsOfs = m_FileHdr.FileLen;
		m_FileHdr.FileLen += WrtLen;
		}

		// now write the header to disk
	Hdr2Disk();

	}
m_bHdrDirty = false;
return(eBSFSuccess);
}



teBSFrsltCodes
CBEDfile::Close(bool bFlush2Disk)
{
teBSFrsltCodes Rslt = eBSFSuccess;
if(bFlush2Disk)
	Rslt = Flush2Disk();
if(m_hFile != -1)		// if file still open
	{
	close(m_hFile);
	m_hFile = -1;
	}
m_szFile[0] = '\0';
Reset(false);
m_bCreate = false;
return(Rslt);
}

// LoadFeatures
// Allocates memory for features and chromosome names, then reads them in from disk
// If specified then file is closed after any errors or when features loaded
teBSFrsltCodes
CBEDfile::LoadFeatures(bool bCloseFile)
{
teBSFrsltCodes Rslt;
int Idx;
int Idy;
int *pOfs;
if(m_hFile == -1)
	return(eBSFerrClosed);
if(!m_FileHdr.NumFeatures)		// any features to load?
	{
	AddErrMsg("CBEDfile::LoadFeatures","No features to load - %s",m_szFile);
	if(bCloseFile)
		{
		close(m_hFile);
		m_hFile = -1;
		}
	return(eBSFerrNoFeatures);
	}

if(m_pFeatures != NULL)
	{
#ifdef _WIN32
	free(m_pFeatures);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pFeatures != MAP_FAILED)
		munmap(m_pFeatures,m_AllocFeaturesSize);
#endif
	m_pFeatures = NULL;
	m_AllocFeaturesSize = 0;
	}

#ifdef _WIN32
m_pFeatures = (tsBEDfeature *) malloc(m_FileHdr.FeaturesSize);	// initial and perhaps the only allocation
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pFeatures = (tsBEDfeature *)mmap(NULL,m_FileHdr.FeaturesSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pFeatures == MAP_FAILED)
	m_pFeatures = NULL;
#endif

if(m_pFeatures == NULL)
	{
	AddErrMsg("CBEDfile::LoadFeatures","Memory allocation of %lld bytes failed - %s",(INT64)m_FileHdr.FeaturesSize,strerror(errno));
	if(bCloseFile)
		{
		close(m_hFile);
		m_hFile = -1;
		}
	m_AllocFeaturesSize = 0;
	return(eBSFerrMem);
	}
m_AllocFeaturesSize = m_FileHdr.FeaturesSize;	

// allocate memory to hold chromosome names
if(m_pChromNames != NULL)
	{
	delete m_pChromNames;
	m_pChromNames = NULL;
	m_AllocChromNamesSize = 0;
	}

if((m_pChromNames = (tsBEDchromname *)new unsigned char [m_FileHdr.ChromNamesSize])==NULL)
	{
	AddErrMsg("CBEDfile::LoadFeatures","Unable to alloc %d bytes of memory to hold feature chrom names - %s",m_FileHdr.ChromNamesSize,m_szFile);
	m_AllocChromNamesSize = 0;
	if(bCloseFile)
		{
		close(m_hFile);
		m_hFile = -1;
		}
	return(eBSFerrMem);
	}
m_AllocChromNamesSize = m_FileHdr.ChromNamesSize;

if(m_ppFeatureNames != NULL)
	{
	delete m_ppFeatureNames;
	m_ppFeatureNames = NULL;
	}

if((m_ppFeatureNames = (tsBEDfeature **)new tsBEDfeature *[m_FileHdr.NumFeatures])==NULL)
	{
	AddErrMsg("CBEDfile::LoadFeatures","Unable to alloc memory to hold %d feature ptrs - %s",m_FileHdr.NumFeatures,m_szFile);
	if(bCloseFile)
		{
		close(m_hFile);
		m_hFile = -1;
		}
	return(eBSFerrMem);
	}
if(m_ppFeatureChromStarts != NULL)
	{
	delete m_ppFeatureChromStarts;
	m_ppFeatureChromStarts = NULL;
	}

if((m_ppFeatureChromStarts = (tsBEDfeature **)new tsBEDfeature *[m_FileHdr.NumFeatures])==NULL)
	{
	AddErrMsg("CBEDfile::LoadFeatures","Unable to alloc memory to hold %d feature chrom ptrs - %s",m_FileHdr.NumFeatures,m_szFile);
	if(bCloseFile)
		{
		close(m_hFile);
		m_hFile = -1;
		}
	return(eBSFerrMem);
	}

// all required memory has been allocated, read from disk
Rslt = ReadDisk(m_FileHdr.ChromNamesOfs,m_FileHdr.ChromNamesSize,m_pChromNames);
if(Rslt == eBSFSuccess)
	{
	if(m_bIsBigEndian)
		{
		tsBEDchromname *pChromName = m_pChromNames;
		for(Idx = 0; Idx < m_FileHdr.NumChroms; Idx++, pChromName++)
			{
			pChromName->ChromID = SwapUI32Endians(pChromName->ChromID);
			pChromName->FirstStartID = SwapUI32Endians(pChromName->FirstStartID);
			pChromName->Hash = SwapUI32Endians(pChromName->Hash);
			pChromName->LastStartID = SwapUI32Endians(pChromName->LastStartID);
			pChromName->MaxFeatLen = SwapUI32Endians(pChromName->MaxFeatLen);
			pChromName->NumFeatures=SwapUI32Endians(pChromName->NumFeatures);
			}
		}
	}

if(Rslt == eBSFSuccess)
	{
	Rslt = ReadDisk(m_FileHdr.FeaturesOfs,m_FileHdr.FeaturesSize,m_pFeatures);
	
	if(m_bIsBigEndian)
		{
		tsBEDfeature *pFeature;
		tsGeneStructure *pGene;
		UINT8 *pBytes;

		pFeature = m_pFeatures;
		pBytes = (UINT8 *)pFeature;
		for(Idx = 0; Idx < m_FileHdr.NumFeatures; Idx++)
			{
			pFeature = (tsBEDfeature *)pBytes;
			pFeature->FeatureID = SwapUI32Endians(pFeature->FeatureID);				// uniquely identifies this feature
			pFeature->Size = SwapUI32Endians(pFeature->Size);						// size in bytes of this instance when concatenated
			pFeature->FiltFlags = SwapUI32Endians(pFeature->FiltFlags);				// used when filtering features in/out from feature bit processing
			pFeature->NameInst = SwapUI32Endians(pFeature->NameInst);				// name instance (1..n)
			pFeature->ChromID = SwapUI32Endians(pFeature->ChromID);					// feature is on this chromosome
			pFeature->Start = SwapUI32Endians(pFeature->Start);						// feature starts at this chromosome offset (0..n)				
			pFeature->End = SwapUI32Endians(pFeature->End);							// feature ends at this chromosome offset (Start + Len -1)
			pFeature->Score = SwapUI32Endians(pFeature->Score);						// feature score
			pFeature->Hash = SwapUI32Endians(pFeature->Hash);						// hash on szName
			if(m_FileHdr.FeatType == eBTGeneExons)
				{
				pGene = (tsGeneStructure *)&pFeature->szName[pFeature->FeatNameLen+1];
				INT32 *pStartEnd;
				pStartEnd = pGene->ExonStartEnds;
				pGene->NumExons = SwapUI32Endians(pGene->NumExons);						// number of exons in this gene
				for(Idy = 0; Idy < pGene->NumExons*2; Idy++,pStartEnd++)
					*pStartEnd = SwapUI32Endians(*pStartEnd);
				pGene->Size  = SwapUI32Endians(pGene->Size);							// size of this instance
				pGene->thickStart = SwapUI32Endians(pGene->thickStart);						// where coding starts 
				pGene->thickEnd = SwapUI32Endians(pGene->thickEnd);						// where coding ends
				}
			pBytes += pFeature->Size;
			}
		}
	}
	

if(Rslt == eBSFSuccess)
	Rslt = ReadDisk(m_FileHdr.FeatureNamesOfs,m_FileHdr.NumFeatures * sizeof(int),m_ppFeatureNames);

#pragma warning(push)
#pragma warning(disable: 4311)
if(Rslt == eBSFSuccess)
        {
        pOfs = (int *)m_ppFeatureNames;
        pOfs += m_FileHdr.NumFeatures-1;
        for(Idx = m_FileHdr.NumFeatures-1;Idx >= 0; Idx--,pOfs--)
			{
			if(m_bIsBigEndian)
				*pOfs = SwapUI32Endians(*pOfs);
            m_ppFeatureNames[Idx] = (tsBEDfeature *)((char *)m_pFeatures + *pOfs);
			}
        }

if(Rslt == eBSFSuccess)
        Rslt = ReadDisk(m_FileHdr.FeatureChromStartsOfs,m_FileHdr.NumFeatures * sizeof(int),m_ppFeatureChromStarts);

if(Rslt == eBSFSuccess)
        {
        pOfs = (int *)m_ppFeatureChromStarts;
        pOfs += m_FileHdr.NumFeatures-1;
        for(Idx = m_FileHdr.NumFeatures-1;Idx >= 0; Idx--, pOfs--)
			{
			if(m_bIsBigEndian)
				*pOfs = SwapUI32Endians(*pOfs);
            m_ppFeatureChromStarts[Idx] = (tsBEDfeature *)((char *)m_pFeatures + *pOfs);
	        }
		}
#pragma warning(pop)
if(bCloseFile)
	{
	close(m_hFile);
	m_hFile = -1;
	}
return(Rslt);
}

// ReadDisk
// Reads block of size 'Len' from disk starting at 'DiskOfs' into preallocated memory at 'pTo'
teBSFrsltCodes
CBEDfile::ReadDisk(INT64 DiskOfs,int Len,void *pTo)
{
if(_lseeki64(m_hFile,DiskOfs,SEEK_SET)!=DiskOfs)
	{
	AddErrMsg("CBEDfile::ReadDisk","Seek failed on %s - %s",m_szFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}
if(read(m_hFile,pTo,Len)!=Len)
	{
	AddErrMsg("CBEDfile::ReadDisk","Read failed on %s - %s",m_szFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}
return(eBSFSuccess);
}

// ContainsGeneDetail
// Returns true if file contains gene structure detail
bool
CBEDfile::ContainsGeneDetail(void)
{
if(!m_bFeaturesAvail)
	return(false);
return(m_FileHdr.FeatType == eBTGeneExons ? true : false);
}

// GetNextFeatureID
// Returns next feature identifier in ChromID->Start-End order
// If CurFeatureID == 0 then returns 1st feature identifier
int
CBEDfile::GetNextFeatureID(int CurFeatureID,int MinScore,int MaxScore)
{
int Idx;
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || CurFeatureID >= m_FileHdr.NumFeatures)
	return(eBSFerrFeature);

if(CurFeatureID < 0)
	CurFeatureID = 0;

for(Idx = CurFeatureID; Idx < m_FileHdr.NumFeatures; Idx++)
	{
	pProbe = m_ppFeatureChromStarts[Idx];
	if(pProbe->Score >= MinScore && pProbe->Score <= MaxScore)
		
		return(pProbe->FeatureID);
	}
return(eBSFerrFeature);
}

// GetPrevFeatureID
// Returns previous feature identifier
// If CurFeatureID >  m_FileHdr.NumFeatures then returns the last feature identifer
int
CBEDfile::GetPrevFeatureID(int CurFeatureID,int MinScore,int MaxScore)
{
int Idx;
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || CurFeatureID < 2)
	return(eBSFerrFeature);

if(CurFeatureID > m_FileHdr.NumFeatures)
	CurFeatureID = m_FileHdr.NumFeatures+1;
for(Idx = CurFeatureID-2; Idx >= 0; Idx--)
	{
	pProbe = m_ppFeatureChromStarts[Idx];
	if(pProbe->Score >= MinScore && pProbe->Score <= MaxScore)
		return(pProbe->FeatureID);
	}
return(eBSFerrFeature);
}


int 
CBEDfile::LocateFeatureBefore(int ChromID,	// feature is on this chromosome
					 int ChromOfs,			// feature ends on or immediately before this offset
 					 int FiltInFlags,		// filter out any features which do not have at least one of the specified filter flags set
					 int FiltOutFlags)		// filter out any features which have at least one of the specified filter flags set
{
tsBEDchromname *pChrom;
int Idx;
tsBEDfeature *pProbe;

int PrevFeatID = eBSFerrFeature;
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(eBSFerrFeature);

if(ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
pChrom = LocateChromName(ChromID);
if(pChrom == NULL)
	return(eBSFerrChrom);
if(pChrom->NumFeatures < 1 || ChromOfs == 0)
	return(eBSFerrFeature);

for(Idx = pChrom->FirstStartID; Idx <= pChrom->LastStartID; Idx++)
	{
	pProbe = m_ppFeatureChromStarts[Idx-1];
	if(!(pProbe->FiltFlags & FiltInFlags) || (pProbe->FiltFlags & FiltOutFlags))
		continue;
	if(pProbe->Start > ChromOfs)
		break;

	if(pProbe->Score >= m_MinScore && pProbe->Score <= m_MaxScore)
		PrevFeatID = pProbe->FeatureID;
	}
return(PrevFeatID);
}

// GetNameInst
// Returns name instance (1..n) for requested feature identifier
int 
CBEDfile::GetNameInst(int FeatureID)
{
tsBEDfeature *pFeature;
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);

pFeature = m_ppFeatureChromStarts[FeatureID-1];
return(pFeature->NameInst);
}

int 
CBEDfile::GetFirstFeatureID(int ChromID,			// first feature is on this chromosome
    				 int FiltInFlags,	// filter out any features which do not have at least one of the specified filter flags set
			 		int FiltOutFlags)  // filter out any features which have at least one of the specified filter flags set
{
tsBEDchromname *pChrom;
int Idx;
tsBEDfeature *pProbe;

int FirstFeatID = eBSFerrFeature;
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(eBSFerrFeature);

if(ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
pChrom = LocateChromName(ChromID);
if(pChrom == NULL)
	return(eBSFerrChrom);

if(pChrom->NumFeatures < 1)
	return(eBSFerrFeature);

for(Idx = pChrom->FirstStartID; Idx <= pChrom->LastStartID; Idx++)
	{
	pProbe = m_ppFeatureChromStarts[Idx-1];
	if(!(pProbe->FiltFlags & FiltInFlags) || (pProbe->FiltFlags & FiltOutFlags))
		continue;
	if(pProbe->Score >= m_MinScore && pProbe->Score <= m_MaxScore)
		{
		FirstFeatID = pProbe->FeatureID;
		break;
		}
	}
return(FirstFeatID);
}

int 
CBEDfile::GetLastFeatureID(int ChromID,				// last feature is on this chromosome
    				 int FiltInFlags,				// filter out any features which do not have at least one of the specified filter flags set
			 		int FiltOutFlags)				// filter out any features which have at least one of the specified filter flags set
{
tsBEDchromname *pChrom;
tsBEDfeature *pProbe;

int Idx;
int LastFeatID = eBSFerrFeature;
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(eBSFerrFeature);

if(ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
pChrom = LocateChromName(ChromID);
if(pChrom == NULL)
	return(eBSFerrChrom);
if(pChrom->NumFeatures < 1)
	return(eBSFerrFeature);

for(Idx = pChrom->LastStartID; Idx >= pChrom->FirstStartID; Idx--)
	{
	pProbe = m_ppFeatureChromStarts[Idx-1];
	if(!(pProbe->FiltFlags & FiltInFlags) || (pProbe->FiltFlags & FiltOutFlags))
		continue;
	if(pProbe->Score >= m_MinScore && pProbe->Score <= m_MaxScore)
		{
		LastFeatID = pProbe->FeatureID;
		break;
		}
	}
return(LastFeatID);
}


int 
CBEDfile::LocateFeatureAfter(int ChromID,	// locate feature is on this chromosome
					 int ChromOfs,			// feature starts on or immediately after this offset
 					 int FiltInFlags,		// filter out any features which do not have at least one of the specified filter flags set
					 int FiltOutFlags)	// filter out any features which have at least one of the specified filter flags set
{
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
if(ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
return(LocateStart(ChromID,ChromOfs,FiltInFlags,FiltOutFlags));
}


// LocateFeatureIDbyName
// Returns feature identifier for Ith feature having specified feature name 
// as there could be multiple features with same name
// The initial feature returned will be that feature having the lowest chrom.start
int
CBEDfile::LocateFeatureIDbyName(char *pszFeatName,	// feature to locate
							int Ith)			    // Ith instance to locate (1..n)
{
tsBEDfeature *pProbe;
char *pszName2;
char c1;
char c2;
int Rslt;

if(Ith < 1 || pszFeatName == NULL || pszFeatName[0] == '\0')
	return(eBSFerrParams);
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(eBSFerrFeature);

int Left = 0;
int Right = m_FileHdr.NumFeatures - 1;
int MidPt;
c1 = tolower(*pszFeatName);

while(Right >= Left) {
	MidPt = (Right + Left)/2;
	pProbe = m_ppFeatureNames[MidPt];
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

	Rslt = stricmp(pszFeatName,pProbe->szName);
	if(Rslt < 0)
		{
		Right = MidPt - 1;
		continue;
		}
	else
		if(Rslt > 0)
			{
			Left = MidPt + 1;
			continue;
			}


	// have a match, but this match may not be match of required name instance (Ith)
	if(pProbe->NameInst == Ith)
		return(pProbe->FeatureID);

	// do linear search forward/backwards until Ith instance of name located
	// performance shouldn't normally be an issue as many name duplications are not expected
	if(pProbe->NameInst > Ith)
		{	
		while(--MidPt >= Left) // linear search backwards
			{
			pProbe = m_ppFeatureNames[MidPt];
			if(stricmp(pszFeatName,pProbe->szName))
				break;
			
			if(pProbe->NameInst == Ith)
				return(pProbe->FeatureID);

			if(pProbe->NameInst < Ith)
				break;
			}
		}
	else
		{	
		while(++MidPt <= Right) // // linear search forwards
			{
			pProbe = m_ppFeatureNames[MidPt];
			if(stricmp(pszFeatName,pProbe->szName))
				break;
			if(pProbe->NameInst == Ith)
				return(pProbe->FeatureID);
			if(pProbe->NameInst > Ith)
				break;
			}
		}
	return(eBSFerrFeature);	// name exists but couldn't locate Ith instance
	}
return(eBSFerrFeature); // couldn't locate any instances of name
}


int					// error or strlen of pszBuff returned if no errors
CBEDfile::GetRemappedBEDFormatedFeature(int FeatureID,// feature instance identifier
					int BuffLen,			  // allocated buffer length (expected to be at least 128 chars) 
					char *pszBuff,			  // where to return the feature formated as TAB delimited BED element
					char *pszRemappedChrom,	  // remapping onto this chrom (NULL if retaining existing chrom)
					char *pszRemappedFeat,	  // remapping as this feature (NULL if retaining existing feature name)
					int  RemappedRelOfs)	  // relative loci (0 if no change)
{
char szBuff[0x03fff];						// BED feature formated into this buffer and then strncpy'd into user supplied pszBuff
int Len;
int Idx;
int SuppInfoLen;
tsGeneStructure *pGene;
tsBEDfeature *pFeature;
tsBEDchromname *pChromName;
if(BuffLen < 128 || pszBuff == NULL)
	return(eBSFerrParams);
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
pFeature = m_ppFeatureChromStarts[FeatureID-1];
pChromName = &m_pChromNames[pFeature->ChromID-1];

if(!(pszRemappedChrom != NULL && pszRemappedChrom[0] != '\0'))
	pszRemappedChrom = pChromName->szName;
if(!(pszRemappedFeat != NULL && pszRemappedFeat[0] != '\0'))
	pszRemappedFeat = pChromName->szName;


Len = sprintf(szBuff,"%s\t%d\t%d\t%s\t%d\t%c",pszRemappedChrom,RemappedRelOfs + pFeature->Start,RemappedRelOfs + pFeature->End + 1,pszRemappedFeat,pFeature->Score,pFeature->Strand);
if(m_FileHdr.FeatType == eBTGeneExons)
	{
	SuppInfoLen = pFeature->Size - (pFeature->FeatNameLen + (int)sizeof(tsBEDfeature));
	pGene = (tsGeneStructure *)&pFeature->szName[pFeature->FeatNameLen+1];

	Len += sprintf(&szBuff[Len],"\t%d\t%d\t%d\t%d\t",pGene->thickStart + RemappedRelOfs + pFeature->Start,pGene->thickEnd + RemappedRelOfs + pFeature->Start,0,pGene->NumExons);
	for(Idx = 0; Idx < pGene->NumExons; Idx++)
		{
		Len += sprintf(&szBuff[Len],"%d,",1 + pGene->ExonStartEnds[(Idx*2)+1] - pGene->ExonStartEnds[Idx*2]); // exon sizes
		}
	Len += sprintf(&szBuff[Len],"\t");
	for(Idx = 0; Idx < pGene->NumExons-1; Idx++)
		{
		Len += sprintf(&szBuff[Len],"%d,",pGene->ExonStartEnds[Idx*2]); // exon starts
		}
	Len += sprintf(&szBuff[Len],"%d",pGene->ExonStartEnds[Idx*2]); // exon starts
	}
strncpy(pszBuff,szBuff,min(BuffLen,Len+1));
pszBuff[min(BuffLen,Len+1)-1] = '\0';
return((int)strlen(pszBuff)); 
}

int					// error or strlen of pszBuff returned if no errors
CBEDfile::GetBEDFormatedFeature(int FeatureID,// feature instance identifier
					int BuffLen,			// allocated buffer length (expected to be at least 128 chars) 
					char *pszBuff)			// where to return the feature formated as TAB delimited BED element
{
char szBuff[0x03fff];						// BED feature formated into this buffer and then strncpy'd into user supplied pszBuff
int Len;
int Idx;
int SuppInfoLen;
tsGeneStructure *pGene;
tsBEDfeature *pFeature;
tsBEDchromname *pChromName;
if(BuffLen < 128 || pszBuff == NULL)
	return(eBSFerrParams);
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
pFeature = m_ppFeatureChromStarts[FeatureID-1];
pChromName = &m_pChromNames[pFeature->ChromID-1];
Len = sprintf(szBuff,"%s\t%d\t%d\t%s\t%d\t%c",pChromName->szName,pFeature->Start,pFeature->End + 1,pFeature->szName,pFeature->Score,pFeature->Strand);
if(m_FileHdr.FeatType == eBTGeneExons)
	{
	SuppInfoLen = pFeature->Size - (pFeature->FeatNameLen + (int)sizeof(tsBEDfeature));
	pGene = (tsGeneStructure *)&pFeature->szName[pFeature->FeatNameLen+1];

	Len += sprintf(&szBuff[Len],"\t%d\t%d\t%d\t%d\t",pGene->thickStart + pFeature->Start,pGene->thickEnd + pFeature->Start,0,pGene->NumExons);
	for(Idx = 0; Idx < pGene->NumExons; Idx++)
		{
		Len += sprintf(&szBuff[Len],"%d,",1 + pGene->ExonStartEnds[(Idx*2)+1] - pGene->ExonStartEnds[Idx*2]); // exon sizes
		}
	Len += sprintf(&szBuff[Len],"\t");
	for(Idx = 0; Idx < pGene->NumExons-1; Idx++)
		{
		Len += sprintf(&szBuff[Len],"%d,",pGene->ExonStartEnds[Idx*2]); // exon starts
		}
	Len += sprintf(&szBuff[Len],"%d",pGene->ExonStartEnds[Idx*2]); // exon starts
	}
strncpy(pszBuff,szBuff,min(BuffLen,Len+1));
pszBuff[min(BuffLen,Len+1)-1] = '\0';
return((int)strlen(pszBuff)); 
}

// GetFeature
// Returns feature details for specified FeatureID
// If NULL is passed in for any specific parameter then that specific feature detail is not returned
teBSFrsltCodes
CBEDfile::GetFeature(int FeatureID,			// feature to return detail for
					 char *pszName,			// where to return feature name
					 char *pszChrom,		// where to return chromosome
					 int *pStart,			// where to return feature start on chromosome (0..n) 
					 int *pEnd,				// where to return feature end on chromosome
 					 int *pScore,			// where to return score
					 char *pStrand,			// where to return strand
					 int MaxSuppInfoLen,	// how much memory has been allocated to pSuppInfo (must be sufficent to contain supp info else pSuppInfo treated as if NULL)
					 void *pSuppInfo)		// where to return any supplementary information
{
tsBEDfeature *pFeature;
tsBEDchromname *pChromName;
int SuppInfoLen;
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);

pFeature = m_ppFeatureChromStarts[FeatureID-1];
pChromName = &m_pChromNames[pFeature->ChromID-1];
if(pszName != NULL)
	strcpy(pszName,pFeature->szName);

if(pSuppInfo != NULL && MaxSuppInfoLen > 1)
	{
	SuppInfoLen = pFeature->Size - (pFeature->FeatNameLen + (int)sizeof(tsBEDfeature));
	if(SuppInfoLen <= MaxSuppInfoLen)
		memmove((char *)pSuppInfo,&pFeature->szName[pFeature->FeatNameLen+1],SuppInfoLen);
	else
		*(char *)pSuppInfo = '\0';
	}
if(pszChrom != NULL)
	strcpy(pszChrom,pChromName->szName);
if(pStart != NULL)
	*pStart = pFeature->Start;
if(pEnd != NULL)
	*pEnd = pFeature->End;
if(pStrand != NULL)
	*pStrand = pFeature->Strand;
if(pScore != NULL)
	*pScore = pFeature->Score;
return(eBSFSuccess);
}

int 
CBEDfile::GetFeatureChromID(int FeatureID)			// returns chromosome identifer on which feature lies
{
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
return(m_ppFeatureChromStarts[FeatureID-1]->ChromID);
}


// MapTransOfs2Loci
// Maps feature or transcript relative ofs to chrom loci, assumes always that the RelOfs is on the '+' strand and returns '+' strand loci 
int										  
CBEDfile::MapTransOfs2Loci(int FeatureID,	 // identifies which feature/transcript
		            int RelOfs,				 // relative offset from start of feature
					char *pStrand,			 // feature is on this strand
					tChromID *pChromID,		 // returned: feature is on this chrom
					int *pLoci)			     // returned: feature relative offset is at this chrom loci
{
tsBEDfeature *pFeature;
tsGeneStructure *pGene;
INT32 *pExonStartEnd;
int ChromLoci;
int Exon;
int ExonLen;

// default to '+' strand and loci 0 with ChromID outside of usual range
if(pStrand != NULL)
	*pStrand = '+';
if(pChromID != NULL)
	*pChromID = 0;
if(pLoci != NULL)
	*pLoci = 0;

// double check that features have been loaded and requested feature is within range
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);

pFeature = m_ppFeatureChromStarts[FeatureID-1];
if(pChromID != NULL)
	*pChromID = pFeature->ChromID;
if(pStrand != NULL)
	*pStrand = pFeature->Strand;
if(pLoci == NULL)
	return(eBSFSuccess);



// if RelOfs before start of feature (could be in upstream promoter region), or features not gene exons then can still easily return the chrom loci 
if(RelOfs < 0 || m_FileHdr.FeatType != eBTGeneExons)	
	{
	ChromLoci = pFeature->Start + RelOfs;
	if(ChromLoci < 0)
		return(eBSFerrFeature);
	*pLoci = ChromLoci;
	return(eBSFSuccess);
	}

pGene = (tsGeneStructure *)&pFeature->szName[pFeature->FeatNameLen+1];
pExonStartEnd = pGene->ExonStartEnds;
for(Exon = 0; Exon < pGene->NumExons; Exon++, pExonStartEnd += 2)
	{
	ChromLoci = pExonStartEnd[0];
	ExonLen = 1 + pExonStartEnd[1] - pExonStartEnd[0];
	if(ExonLen <= RelOfs)
		{
		RelOfs -= ExonLen;
		ChromLoci = pExonStartEnd[1] + 1;
		continue;
		}
	break;
	}
*pLoci = pFeature->Start + ChromLoci + RelOfs;
return(eBSFSuccess);
}

// GetChromosome
// returns chromosome detail
// If NULL is passed in for any specific parameter then that specific chromosome detail is not returned
teBSFrsltCodes 
CBEDfile::GetChromosome(int ChromID,		// chromosome identifier
					 char *pszChrom,		// where to return chromosome name
					 int *pNumFeatures,		// where to return number of features on this chromosome
					 int *pFirstStartID,		// where to return identifier of first feature on this chromosome
					 int *pLastStartID)		// where to return identifier of last feature on this chromosome
{
if(!m_bFeaturesAvail)
	return(eBSFerrFeature);
if(ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);

tsBEDchromname *pChrom = LocateChromName(ChromID);
if(pChrom == NULL)
	return(eBSFerrChrom);
if(pszChrom != NULL)
	strcpy(pszChrom,pChrom->szName);
if(pNumFeatures != NULL)
	*pNumFeatures = pChrom->NumFeatures;
if(pFirstStartID != NULL)
	*pFirstStartID = pChrom->FirstStartID; // start indexes are same as feature identifiers 
if(pLastStartID != NULL)
	*pLastStartID = pChrom->LastStartID;
return(eBSFSuccess);
}

// LocateChromIDbyName
// Returns chromosome identifier for specified name
// NOTE: treats 'chrm' and 'mitichondria' as being synomyous
// NOTE: treats 'chrc' and 'chloroplast' as being synomyous
int													// returned chromosome identifer
CBEDfile::LocateChromIDbyName(char *pszChromName)	// chromosome name to locate
{
tsBEDchromname *pChrom;
if(!m_bFeaturesAvail)
	return(eBSFerrFeature);
if(pszChromName == NULL || *pszChromName == '\0')
	return(eBSFerrParams);
pChrom = LocateChromName(pszChromName);
if(pChrom == NULL)
	{
	if(!stricmp(pszChromName,"chloroplast"))
		pChrom = LocateChromName((char *)"ChrC");
	else
		if(!stricmp(pszChromName,"mitochondria"))
			pChrom = LocateChromName((char *)"ChrM");
	if(pChrom == NULL)
		{
		if(!stricmp(pszChromName,"ChrC"))
			pChrom = LocateChromName((char *)"chloroplast");
		else
			if(!stricmp(pszChromName,"ChrM"))
				pChrom = LocateChromName((char *)"mitochondria");
		}
	}
return(pChrom == NULL ? eBSFerrChrom : pChrom->ChromID);
}


// LocateStart
// Returns feature identifier which has a Start which is immediately above OverLapsOfs
// Returns:
// <0       Error
// 0		No feature has an Start which is immediately above OverLapsOfs
// >0		feature identifier
// Note that there could be multiple features which start at same locus immediately above OverLapsOfs
int										// returned feature identifier
CBEDfile::LocateStart(int ChromID,		// feature is on which chromosome
					int OverLapsOfs,	// a point on the chromosome which returned feature starts above chrom.ofs
					 int FiltInFlags, // filter out any features which do not have at least one of the specified filter flags set
					 int FiltOutFlags) // filter out any features which have at least one of the specified filter flags set
{
int StartIdx;
int EndIdx;
int MidIdx;
tsBEDchromname *pChrom;
tsBEDfeature *pProbe;

//sanity checks whilst debugging
#ifdef _DEBUG
if(OverLapsOfs < 0 || ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(eBSFerrParams);
if(!m_bFeaturesAvail)
	return(eBSFerrFeature);
if(ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
#endif

pChrom = LocateChromName(ChromID);
if(pChrom == NULL)
	return(eBSFerrChrom);


StartIdx = pChrom->FirstStartID - 1;
EndIdx = pChrom->LastStartID - 1;

while(EndIdx >= StartIdx) {
	MidIdx = (EndIdx + StartIdx)/2;
	pProbe = (tsBEDfeature *)m_ppFeatureChromStarts[MidIdx];
	if(pProbe->Start > OverLapsOfs)	
		{
		EndIdx = MidIdx - 1;
		continue;
		}

	if(pProbe->Start <= OverLapsOfs)	
		{
		StartIdx = MidIdx + 1;
		continue;
		}
	}

if(pProbe->Start > OverLapsOfs && (pProbe->Score >= m_MinScore && pProbe->Score <= m_MaxScore))
	{
	return(pProbe->FeatureID);
	}
// ended up on element which has Start <= OverLapsOfs so advance to 1st element filtering as appropriate
// which has Start >  OverLapsOfs
while(++MidIdx <= (pChrom->LastStartID - 1))
	{
	pProbe = (tsBEDfeature *)m_ppFeatureChromStarts[MidIdx];
	if(pProbe->Score < m_MinScore || pProbe->Score > m_MaxScore)
		continue;
	if((m_OnStrand != '*' && m_OnStrand != pProbe->Strand) || !(pProbe->FiltFlags & FiltInFlags) || pProbe->FiltFlags & FiltOutFlags)
		continue;
	if(pProbe->Start > OverLapsOfs)
		return(pProbe->FeatureID);
	}
return(0);
}



// LocateFeatureIDonChrom
// Returns feature identifier which overlaps a specific chrom.ofs
// As there can be multiple features overlaping then set Ith to 1..n to return the nth overlaping feature
// Returns:
// <0	Error
// 0    No feature which overlaps
// >0   Feature identifier
int										 // returned feature identifier
CBEDfile::LocateFeatureIDonChrom(int ChromID, // feature is on which chromosome
							 int OverLapsOfs, // a point on the chromosome on which returned feature is to overlap by at least one base
							 int Ith,		  // Ith instance to overlap (1..n)
 							 int FiltInFlags, // filter out any features which do not have at least one of the specified filter flags set
  							 int FiltOutFlags) // filter out any features which have at least one of the specified filter flags set
{
int StartIdx;
int EndIdx;
int StartID;
tsBEDchromname *pChrom;
tsBEDfeature *pProbe;

//sanity checks whilst debugging
#ifdef _DEBUG
if(OverLapsOfs < 0 || Ith < 1 || ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(eBSFerrParams);
if(!m_bFeaturesAvail)
	return(eBSFerrFeature);
if(ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
#endif


pChrom = LocateChromName(ChromID);
if(pChrom == NULL)
	return(eBSFerrChrom);
if(!pChrom->NumFeatures)
	return(0);

// determine initial lower limit to start linear search from which is a feature ending before OverLapsOfs
// search terminates when start > OverLapsOfs
if(pChrom->NumFeatures > 50 && OverLapsOfs > pChrom->MaxFeatLen) 
	{
	StartID = LocateStart(ChromID,OverLapsOfs - pChrom->MaxFeatLen,FiltInFlags,FiltOutFlags);
	if(StartID == 0)
		StartIdx = pChrom->FirstStartID - 1;
	else
		StartIdx = StartID - 1;
	}
else
	StartIdx = pChrom->FirstStartID - 1;
EndIdx = pChrom->LastStartID - 1;
while(StartIdx <= EndIdx)
	{
	pProbe = m_ppFeatureChromStarts[StartIdx];
	if(pProbe->Start <= OverLapsOfs &&
	   pProbe->End >= OverLapsOfs)
		{
		if(pProbe->Score < m_MinScore || pProbe->Score > m_MaxScore)
			continue;
		if((m_OnStrand != '*' && m_OnStrand != pProbe->Strand) || !(pProbe->FiltFlags & FiltInFlags) || pProbe->FiltFlags & FiltOutFlags)
			continue;
		if(--Ith)
			{
			return(pProbe->FeatureID);
			}
		continue;
		}
	if(pProbe->Start > OverLapsOfs)
		break;
	}

return(0);
}


int											  // returned feature identifier
CBEDfile::LocateFeatureIDinRangeOnChrom(int ChromID, // feature is on which chromosome
							 int StartOfs,       // feature must end on or after Start
							 int EndOfs,		  // and start on or before End 
							 int Ith,		  // Ith instance to return (1..n)
 							 int FiltInFlags, // filter out any features which do not have at least one of the specified filter flags set
  							 int FiltOutFlags) // filter out any features which have at least one of the specified filter flags set
{
int StartIdx;
int EndIdx;
int StartID;
tsBEDfeature *pProbe;

// sanity checks whilst debug
#ifdef _DEBUG
if(StartOfs < 0 || StartOfs > EndOfs || Ith < 1)
	return(eBSFerrParams);
if(!m_bFeaturesAvail)
	return(eBSFerrFeature);
if(ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
#endif

tsBEDchromname *pChrom = LocateChromName(ChromID);
if(pChrom == NULL)
	return(eBSFerrChrom);

// determine initial lower limit to start linear search from which is a feature ending before StartOfs
// search terminates when start > StartOfs
if(pChrom->NumFeatures > 20 && StartOfs > pChrom->MaxFeatLen) 
	{
	StartID = LocateStart(ChromID,StartOfs - pChrom->MaxFeatLen,FiltInFlags,FiltOutFlags);
	if(StartID == 0)
		StartIdx = pChrom->FirstStartID - 1;
	else
		StartIdx = StartID - 1;
	}
else
	StartIdx = pChrom->FirstStartID - 1;
EndIdx = pChrom->LastStartID - 1;
while(StartIdx <= EndIdx)
	{
	pProbe = m_ppFeatureChromStarts[StartIdx++];
	if(pProbe->Start <= EndOfs &&
	   pProbe->End >= StartOfs)
		{
		if(pProbe->Score < m_MinScore || pProbe->Score > m_MaxScore)
			continue;
		if((m_OnStrand != '*' && m_OnStrand != pProbe->Strand) || !(pProbe->FiltFlags & FiltInFlags) || pProbe->FiltFlags & FiltOutFlags)
			continue;
		if(--Ith == 0)
			{
			return(pProbe->FeatureID);
			}
		continue;
		}
	if(pProbe->Start > EndOfs)
		break;
	}
return(0);
}



// returns the number of chromosomes
int 
CBEDfile::GetNumChromosomes(void)
{
if(!m_bFeaturesAvail)
	return(eBSFerrFeature);
return(m_FileHdr.NumChroms);
}



// return number of features on specified chromosome
// if ChromID == 0 then returns total over all chromosomes
int 
CBEDfile::GetNumFeatures(int ChromID)
{
tsBEDchromname *pChrom;
#ifdef _DEBUG
// sanity checks!
if(!m_bFeaturesAvail)
	return(eBSFerrFeature);
if(ChromID < 0 || ChromID > m_FileHdr.NumChroms)
	return(eBSFerrParams);
#endif

if(ChromID == 0)
	return(m_FileHdr.NumFeatures);
pChrom = LocateChromName(ChromID);
if(pChrom == NULL)
	return(eBSFerrChrom);
return(pChrom->NumFeatures);
}

// returns the provisional number of chromosomes
int 
CBEDfile::GetProvNumChromosomes(void)
{
return(m_FileHdr.NumChroms);
}

// return provisional number of features on all chromosome
int 
CBEDfile::GetProvNumFeatures(void)
{
return(m_FileHdr.NumFeatures);
}

int										  // returned number of features
CBEDfile::GetNumFeatures(int ChromID,	  // features are on which chromosome
					   int StartOfs,	  // features must end on or after Start
					   int EndOfs,	      // and start on or before End 
					 int FiltInFlags, // filter out any features which do not have at least one of the specified filter flags set
					 int FiltOutFlags) // filter out any features which have at least one of the specified filter flags set
{
int NumFeatures;
int StartIdx;
int EndIdx;
int StartID;
tsBEDfeature *pProbe;

// sanity checks whilst debug
#ifdef _DEBUG
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
if(StartOfs < 0 || StartOfs > EndOfs)
	return(eBSFerrParams);
if(ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
#endif

tsBEDchromname *pChrom = LocateChromName(ChromID);
if(pChrom == NULL)
	return(eBSFerrChrom);

// determine initial lower limit to start linear search from which is a feature ending before StartOfs
// search terminates when start > StartOfs
if(pChrom->NumFeatures > 20 && StartOfs > pChrom->MaxFeatLen) 
	{
	StartID = LocateStart(ChromID,StartOfs - pChrom->MaxFeatLen,FiltInFlags,FiltOutFlags);
	if(StartID == 0)
		StartIdx = pChrom->FirstStartID - 1;
	else
		StartIdx = StartID - 1;
	}
else
	StartIdx = pChrom->FirstStartID - 1;
EndIdx = pChrom->LastStartID - 1;
NumFeatures = 0;
while(StartIdx <= EndIdx)
	{
	pProbe = m_ppFeatureChromStarts[StartIdx];
	if((m_OnStrand == '*' || m_OnStrand == pProbe->Strand) &&
		pProbe->FiltFlags & FiltInFlags && !(pProbe->FiltFlags & FiltOutFlags) &&
			pProbe->Start <= EndOfs &&
			pProbe->End >= StartOfs)
		{
		if(pProbe->Score >= m_MinScore && pProbe->Score <= m_MaxScore)
			NumFeatures++;
		continue;
		}
	if(pProbe->Start > EndOfs)
		break;
	}
return(NumFeatures);
}

// returns start of the Ith (1..n) feature on chromosome 
int 
CBEDfile::GetFeatStart(int ChromID,int Ith)
{
tsBEDchromname *pChrom;
tsBEDfeature *pProbe;
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
if(Ith < 1 || ChromID < 1)
	return(eBSFerrParams);
if(ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
pChrom = LocateChromName(ChromID);
if(pChrom == NULL)
	return(eBSFerrChrom);
if(Ith > pChrom->NumFeatures)
	return(eBSFerrFeature);
pProbe = m_ppFeatureChromStarts[pChrom->FirstStartID + Ith - 2];
return(pProbe->Start);
}

// returns end of the Ith  (1..n) feature on chromosome ChromID
int 
CBEDfile::GetFeatEnd(int ChromID,int Ith)
{
tsBEDchromname *pChrom;
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
if(Ith < 1 || ChromID < 1)
	return(eBSFerrParams);
if(ChromID > m_FileHdr.NumChroms)
	return(eBSFerrChrom);
pChrom = LocateChromName(ChromID);
if(pChrom == NULL)
	return(eBSFerrChrom);
if(Ith > pChrom->NumFeatures)
	return(eBSFerrFeature);
pProbe = m_ppFeatureChromStarts[pChrom->FirstStartID + Ith - 2];
return(pProbe->End);
}


// SortChromStarts
// Used to sort by ChromID ---> Start --->End ---> NameInst
int 
CBEDfile::SortChromStarts( const void *arg1, const void *arg2)
{
tsBEDfeature *pEl1 = *(tsBEDfeature **)arg1;
tsBEDfeature *pEl2 = *(tsBEDfeature **)arg2;

if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->Start < pEl2->Start)
	return(-1);
if(pEl1->Start > pEl2->Start)
	return(1);
if(pEl1->End < pEl2->End)
	return(-1);
if(pEl1->End > pEl2->End)
	return(1);
if(pEl1->NameInst < pEl2->NameInst)
	return(-1);
if(pEl1->NameInst > pEl2->NameInst)
	return(1);
return(0);
}


// SortFeatureNames
// Used to sort by feature names --> ChromID ---> Start ---> End
int 
CBEDfile::SortFeatureNames( const void *arg1, const void *arg2)
{
int Rslt;
tsBEDfeature *pEl1 = *(tsBEDfeature **)arg1;
tsBEDfeature *pEl2 = *(tsBEDfeature **)arg2;

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
if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->Start < pEl2->Start)
	return(-1);
if(pEl1->Start > pEl2->Start)
	return(1);
if(pEl1->End < pEl2->End)
	return(-1);
if(pEl1->End > pEl2->End)
	return(1);
return(0);
}


// SortChromNames
// Used to sort chromosome names
int 
CBEDfile::SortChromNames( const void *arg1, const void *arg2)
{
tsBEDchromname *pEl1 = (tsBEDchromname *)arg1;
tsBEDchromname *pEl2 = (tsBEDchromname *)arg2;

char *pszName1 = &pEl1->szName[0];
char *pszName2 = &pEl2->szName[0];
char c1 = tolower(*pszName1);
char c2 = tolower(*pszName2);
if(c1 < c2)
	return(-1);
if(c1 > c2)
	return(1);
return(stricmp(pszName1,pszName2));
}

int
CBEDfile::SortU2S( const void *arg1, const void *arg2)
{
tsUnsortToSortID *pEl1 = (tsUnsortToSortID *)arg1;
tsUnsortToSortID *pEl2 = (tsUnsortToSortID *)arg2;
if(pEl1->OldID < pEl2->OldID)
	return(-1);
if(pEl1->OldID > pEl2->OldID)
	return(1);
return(0);
}

// returns number of exons - includes UTRs + CDS
int 
CBEDfile::GetNumExons(int FeatureID)
{
tsGeneStructure *pGene;
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
if(m_FileHdr.FeatType != eBTGeneExons)
	return(eBSFerrGene);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
pGene = (tsGeneStructure *)&pProbe->szName[pProbe->FeatNameLen+1];
return(pGene->NumExons);
}

int 
CBEDfile::GetNumIntrons(int FeatureID)  // returns number of introns
{
int Rslt = GetNumExons(FeatureID);
return(Rslt <= 0 ? Rslt : Rslt -1);
}

int 
CBEDfile::GetTranscribedLen(int FeatureID)		// returns total transcribed length for specified feature
{
int Exon;
int Len;
tsGeneStructure *pGene;
tsBEDfeature *pProbe;
INT32 *pExonStartEnd;
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
if(m_FileHdr.FeatType != eBTGeneExons)
	return(eBSFerrGene);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
pGene = (tsGeneStructure *)&pProbe->szName[pProbe->FeatNameLen+1];
Len = 0;
pExonStartEnd = pGene->ExonStartEnds;
for(Exon = 0; Exon < pGene->NumExons; Exon++, pExonStartEnd += 2)
	Len += 1 + pExonStartEnd[1] - pExonStartEnd[0];
return(Len);
}

int 
CBEDfile::GetCDSLen(int FeatureID)				// returns CDS transcribed length for specified feature
{
int Exon;
int Len;
tsGeneStructure *pGene;
tsBEDfeature *pProbe;
INT32 *pExonStartEnd;
UINT32 StartLoci;
UINT32 EndLoci;
UINT32 CDSstart;
UINT32 CDSend;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
if(m_FileHdr.FeatType != eBTGeneExons)
	return(eBSFerrGene);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
pGene = (tsGeneStructure *)&pProbe->szName[pProbe->FeatNameLen+1];

CDSstart = pGene->thickStart;
CDSend = pGene->thickEnd;
if(CDSstart == CDSend)
	return(0);

pExonStartEnd = pGene->ExonStartEnds;
Len = 0;
for(Exon = 0; Exon < pGene->NumExons; Exon++, pExonStartEnd += 2)
	{
	StartLoci = pExonStartEnd[0];
	EndLoci   = pExonStartEnd[1];
	if(EndLoci < CDSstart || StartLoci > CDSend)
		continue;
	if(StartLoci < CDSstart)
		StartLoci = CDSstart;
	if(EndLoci > CDSend)
		EndLoci = CDSend;
	if(StartLoci <= EndLoci)
		Len += 1 + EndLoci - StartLoci;
	}

return(Len);
}

int 
CBEDfile::Get3UTRLen(int FeatureID)				// returns 3' UTR transcribed length for specified feature
{
int Exon;
int Len;
tsGeneStructure *pGene;
tsBEDfeature *pProbe;
INT32 *pExonStartEnd;
UINT32 StartLoci;
UINT32 EndLoci;
UINT32 CDSstart;
UINT32 CDSend;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
if(m_FileHdr.FeatType != eBTGeneExons)
	return(eBSFerrGene);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
pGene = (tsGeneStructure *)&pProbe->szName[pProbe->FeatNameLen+1];

CDSstart = pGene->thickStart;
CDSend = pGene->thickEnd;
if(CDSstart == CDSend)
	return(0);

pExonStartEnd = pGene->ExonStartEnds;
Len = 0;
for(Exon = 0; Exon < pGene->NumExons; Exon++, pExonStartEnd += 2)
	{
	StartLoci = pExonStartEnd[0];
	EndLoci   = pExonStartEnd[1];

	if(EndLoci <= CDSend && StartLoci >= CDSstart) // is exon CDS only?
		continue;

	if(pProbe->Strand == '-')
		{
		// check if 3' UTR on '-' strand 
		if(StartLoci < CDSstart)
			{
			if(EndLoci >= CDSstart)
				EndLoci = CDSstart - 1;
			}
		}
	else 
		{	
		// check if 3'UTR on '+' strand
		if(EndLoci > CDSend)
			{
			if(StartLoci <= CDSend)
				StartLoci = CDSend+1;
			}
		}
	Len += 1 + EndLoci - StartLoci;
	}
return(Len);
}

int 
CBEDfile::Get5UTRLen(int FeatureID)				// returns 5' transcribed length for specified feature
{
int Exon;
int Len;
tsGeneStructure *pGene;
tsBEDfeature *pProbe;
INT32 *pExonStartEnd;
UINT32 StartLoci;
UINT32 EndLoci;
UINT32 CDSstart;
UINT32 CDSend;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
if(m_FileHdr.FeatType != eBTGeneExons)
	return(eBSFerrGene);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
pGene = (tsGeneStructure *)&pProbe->szName[pProbe->FeatNameLen+1];

CDSstart = pGene->thickStart;
CDSend = pGene->thickEnd;
if(CDSstart == CDSend)
	return(0);

pExonStartEnd = pGene->ExonStartEnds;
Len = 0;
for(Exon = 0; Exon < pGene->NumExons; Exon++, pExonStartEnd += 2)
	{
	StartLoci = pExonStartEnd[0];
	EndLoci   = pExonStartEnd[1];

	if(EndLoci <= CDSend && StartLoci >= CDSstart) // is exon CDS only?
		continue;

	if(pProbe->Strand != '-')
		{
		// check if 5' UTR on '+' strand 
		if(StartLoci < CDSstart)
			{
			if(EndLoci >= CDSstart)
				EndLoci = CDSstart - 1;
			}
		}
	else 
		{	
		// check if 5'UTR on '-' strand
		if(EndLoci > CDSend)
			{
			if(StartLoci <= CDSend)
				StartLoci = CDSend+1;
			}
		}
	Len += 1 + EndLoci - StartLoci;
	}
return(Len);
}

int 
CBEDfile::GetFeatLen(int FeatureID)				// returns total length for specified feature
{
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
return(1 + pProbe->End - pProbe->Start);
}

teBSFrsltCodes
CBEDfile::InitUserClass(int DefltClass)
{
tsBEDfeature *pProbe;
int Idx;
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(eBSFerrFeature);

for(Idx = 0; Idx < m_FileHdr.NumFeatures; Idx++)
	{
	pProbe = m_ppFeatureChromStarts[Idx];
	pProbe->UserClass = (DefltClass & 0x07fffffff);
	}
return(eBSFSuccess);
}

teBSFrsltCodes 
CBEDfile::SetUserClass(char *pszFeatName,int UserClass)	
{
tsBEDfeature *pProbe;

bool bSet;
int FeatureID;
int Nth = 1;
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
bSet = false;
while((FeatureID = LocateFeatureIDbyName(pszFeatName,Nth++)) > 0)
	{
	pProbe = m_ppFeatureChromStarts[FeatureID-1];
	pProbe->UserClass = UserClass & 0x7fffffff;
	bSet = true;
	}
return(bSet ? eBSFSuccess : eBSFerrFeature);
}

int 
CBEDfile::GetUserClass(char *pszFeatName)	
{
tsBEDfeature *pProbe;
int FeatureID;
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(eBSFerrFeature);

if((FeatureID = LocateFeatureIDbyName(pszFeatName,1)) > 0)
	{
	pProbe = m_ppFeatureChromStarts[FeatureID-1];
	return((teBSFrsltCodes)(pProbe->UserClass & 0x7fffffff));
	}
return(eBSFerrFeature);
}

teBSFrsltCodes 
CBEDfile::SetUserClass(int FeatureID, int UserClass)		// set user classification for feature instance identifier
{
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
pProbe->UserClass = UserClass & 0x7fffffff;
return(eBSFSuccess);
}

int 
CBEDfile::GetUserClass(int FeatureID)		// get user classification for feature instance identifier
{
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
return((teBSFrsltCodes)(pProbe->UserClass & 0x7fffffff));
}

int 
CBEDfile::GetFeatScore(int FeatureID)				// returns score for this feature
{
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
return(pProbe->Score);
}

// returns relative offset at which CDS starts in gene
// NOTE: offset is strand specific, add to gene start for genes on '+' strand, subtract from gene start on '-' strand
int 
CBEDfile::GetCDSStart(int FeatureID)
{
tsGeneStructure *pGene;
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
if(m_FileHdr.FeatType != eBTGeneExons)
	return(eBSFerrGene);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
pGene = (tsGeneStructure *)&pProbe->szName[pProbe->FeatNameLen+1];
return(pGene->thickStart);
}

// returns relative offset at which CDS ends in gene
// NOTE: offset is strand specific, add to gene start for genes on '+' strand, subtract from gene start on '-' strand
int 
CBEDfile::GetCDSEnd(int FeatureID)
{
tsGeneStructure *pGene;
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
if(m_FileHdr.FeatType != eBTGeneExons)
	return(eBSFerrGene);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
pGene = (tsGeneStructure *)&pProbe->szName[pProbe->FeatNameLen+1];
return(pGene->thickEnd);
}


// returns relative start offset of the ExonID'th (1..n)
// NOTE: offset is strand specific, add to gene start for genes on '+' strand, subtract from gene start on '-' strand
int 
CBEDfile::GetExonStart(int FeatureID,int ExonID)
{
tsGeneStructure *pGene;
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);

if(ExonID < 1)
	return(eBSFerrExon);

pProbe = m_ppFeatureChromStarts[FeatureID-1];
if(m_FileHdr.FeatType != eBTGeneExons)
	{
	if(ExonID != 1)
		return(eBSFerrExon);
	return(pProbe->Start);
	}

pGene = (tsGeneStructure *)&pProbe->szName[pProbe->FeatNameLen+1];
if(ExonID > pGene->NumExons)
	return(eBSFerrExon);
return(pGene->ExonStartEnds[(ExonID-1) * 2] + pProbe->Start);
}

// returns relative end offset of the ExonID'th (1..n)
// NOTE: offset is strand specific, add to gene start for genes on '+' strand, subtract from gene start on '-' strand
int 
CBEDfile::GetExonEnd(int FeatureID,int ExonID)
{
tsGeneStructure *pGene;
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);

if(ExonID < 1)
	return(eBSFerrExon);

pProbe = m_ppFeatureChromStarts[FeatureID-1];
if(m_FileHdr.FeatType != eBTGeneExons)
	{
	if(ExonID != 1)
		return(eBSFerrExon);
	return(pProbe->End);
	}

pGene = (tsGeneStructure *)&pProbe->szName[pProbe->FeatNameLen+1];
if(ExonID > pGene->NumExons)
	return(eBSFerrExon);
return(pGene->ExonStartEnds[((ExonID-1) * 2)+1] + pProbe->Start);
}

// returns relative start offset of the IntronID'th (1..n)
// NOTE: offset is strand specific, add to gene start for genes on '+' strand, subtract from gene start on '-' strand
int 
CBEDfile::GetIntronStart(int FeatureID,int IntronID)
{
tsGeneStructure *pGene;
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
if(m_FileHdr.FeatType != eBTGeneExons)
	return(eBSFerrGene);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
pGene = (tsGeneStructure *)&pProbe->szName[pProbe->FeatNameLen+1];
if(IntronID < 1 || IntronID >= pGene->NumExons)
	return(eBSFerrIntron);
return(pGene->ExonStartEnds[((IntronID-1) * 2)+1] + pProbe->Start + 1);
}


// returns relative end offset of the IntronID'th (1..n)
// NOTE: offset is strand specific, add to gene start for genes on '+' strand, subtract from gene start on '-' strand
int 
CBEDfile::GetIntronEnd(int FeatureID,int IntronID)
{
tsGeneStructure *pGene;
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);

if(m_FileHdr.FeatType != eBTGeneExons)
	return(eBSFerrGene);

pProbe = m_ppFeatureChromStarts[FeatureID-1];
pGene = (tsGeneStructure *)&pProbe->szName[pProbe->FeatNameLen+1];
if(IntronID < 1 || IntronID >= pGene->NumExons)
	return(eBSFerrIntron);
return(pGene->ExonStartEnds[(IntronID) * 2] + pProbe->Start - 1);
} 

// following functions assume co-ordinates are specified along the '+' strand

// returns true if any point between StartOfs and EndOfs overlaps any feature
bool 
CBEDfile::InAnyFeature(int ChromID,int StartOfs,int EndOfs)
{
int FeatID;
// sanity checks!
#ifdef _DEBUG
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(false);
if(!ChromID || ChromID > m_FileHdr.NumChroms)
	return(false);
if(StartOfs < 0 || StartOfs > EndOfs)
	return(false);
#endif

if((FeatID = LocateFeatureIDinRangeOnChrom(ChromID,StartOfs,EndOfs,1,cFeatFiltIn,cFeatFiltOut))>0)
	return(true);
return(false);
}


int
CBEDfile::GetFeatureOverlaps(int ChkFeatBits,// which feature bits to check for
				   int FeatID,				// feature identifier
				   int StartOfs,			// start along cromosome 0..n
				   int EndOfs,				// end along chromosome start + len - 1
				   int Distance)			// used if cFeatBitUpstream or cFeatBitDnstream to specify max distance away from gene start/end
{
int OverlapFlgs;
int RelStartOfs;
int RelEndOfs;
int Idx;
tsGeneStructure *pGene;

// no point in determining if in up/dnstream when length of up/dnstream region in <= 0
if(Distance <= 0)
	ChkFeatBits &= ~(cFeatBitUpstream | cFeatBitDnstream);
if(!(ChkFeatBits & cAnyFeatBits))
	return(0);

OverlapFlgs = 0;
tsBEDfeature *pProbe = m_ppFeatureChromStarts[FeatID-1];

if(pProbe->Score < m_MinScore || pProbe->Score > m_MaxScore)
	return(OverlapFlgs);
if((m_OnStrand != '*' && m_OnStrand != pProbe->Strand))
	return(OverlapFlgs);

// check for upstream overlap
if(ChkFeatBits & cFeatBitUpstream)
	{
	if(pProbe->Strand == '+')
		{
		if(StartOfs < pProbe->Start &&
			EndOfs >= pProbe->Start - Distance)
			OverlapFlgs |= cFeatBitUpstream;
		}
	else	// else on '-' strand
		{
		if(StartOfs <= (pProbe->End + Distance) &&
			EndOfs > pProbe->End)
			OverlapFlgs |= cFeatBitUpstream;
		}
	}

// check for downstream overlap
if(ChkFeatBits & cFeatBitDnstream)
	{
	if(pProbe->Strand == '+')
		{
		if(StartOfs <= pProbe->End + Distance &&
			EndOfs > pProbe->End)
			OverlapFlgs |= cFeatBitDnstream;
		}
	else	// else on '-' strand
		{	
		if(StartOfs < pProbe->Start &&
				EndOfs > pProbe->Start - Distance)
				OverlapFlgs |= cFeatBitDnstream;
		}
	}


// none of the remaining feature bit tests will succeed if no overlaps onto the gene/feature so may as well check
if(StartOfs > pProbe->End || EndOfs < pProbe->Start)	
	return(OverlapFlgs);

if(ChkFeatBits & cFeatBitFeature)
	OverlapFlgs |= cFeatBitFeature;

if(m_FileHdr.FeatType != eBTGeneExons)
	return(OverlapFlgs);

// make start/end relative to gene start
if(StartOfs <= pProbe->Start)
	RelStartOfs = 0;
else
	RelStartOfs = StartOfs - pProbe->Start;
if(EndOfs > pProbe->End)
    RelEndOfs = pProbe->End - pProbe->Start;
else
	RelEndOfs = EndOfs - pProbe->Start;

pGene = (tsGeneStructure *)&pProbe->szName[pProbe->FeatNameLen+1];

// check for overlaps with exons  
if(ChkFeatBits & (cFeatBit5UTR | cFeatBitCDS | cFeatBit3UTR | cFeatBitExons))
	{
	for(Idx = 0; Idx < pGene->NumExons; Idx++)
		{
		if(RelEndOfs < pGene->ExonStartEnds[Idx*2])
			break;
		if(RelStartOfs <= pGene->ExonStartEnds[(Idx*2)+1] &&	// does at least one base overlap this exon?
			RelEndOfs >= pGene->ExonStartEnds[Idx*2])
			{
			// at least one base is inside an exon
			if(ChkFeatBits & cFeatBitExons)
				OverlapFlgs |= (ChkFeatBits & cFeatBitExons);
			if(ChkFeatBits & cFeatBit5UTR && !(OverlapFlgs & cFeatBit5UTR))
				{
				if(pProbe->Strand == '+')
					{
					if(RelStartOfs < pGene->thickStart &&	// does at least one base UTR extend into this exon?
						pGene->thickStart > pGene->ExonStartEnds[Idx*2])
						OverlapFlgs |= cFeatBit5UTR;
					}
				else									// feature is on '-' strand
					{
					if(RelEndOfs > pGene->thickEnd &&	// does at least one base overlap 5' UTR?
						pGene->thickEnd < pGene->ExonStartEnds[(Idx*2)+1])
						OverlapFlgs |= cFeatBit5UTR;
					}
				}
			if(ChkFeatBits & cFeatBit3UTR && !(OverlapFlgs & cFeatBit3UTR))
				{
				if(pProbe->Strand == '+')
					{
					if(RelEndOfs > pGene->thickEnd &&	// does at least one base overlap 3' UTR
						pGene->thickEnd < pGene->ExonStartEnds[(Idx*2)+1])
						OverlapFlgs |= cFeatBit3UTR;
					}
				else									// feature is on '-' strand
					{
					if(RelStartOfs < pGene->thickStart &&	// does at least one base overlap 3' UTR
						pGene->thickStart > pGene->ExonStartEnds[Idx*2])
						OverlapFlgs |= cFeatBit3UTR;
					}
				}

			// check if any part in a CDS
			if(RelStartOfs <= pGene->thickEnd &&
				RelEndOfs >= pGene->thickStart)
				OverlapFlgs |= cFeatBitCDS;

			// no point in iterating once all requested bits have been set
			if((OverlapFlgs & (cFeatBit5UTR | cFeatBitCDS | cFeatBit3UTR | cFeatBitExons)) == 
					(ChkFeatBits & (cFeatBit5UTR | cFeatBitCDS | cFeatBit3UTR | cFeatBitExons)))
				break;
			}

		}
	}

	// check for intron overlaps
if(ChkFeatBits & cFeatBitIntrons && pGene->NumExons > 1)
	{
	for(Idx = 0; Idx < pGene->NumExons-1; Idx++)
		{
		if(RelStartOfs < pGene->ExonStartEnds[(Idx+1)*2] &&
			RelEndOfs > pGene->ExonStartEnds[(Idx*2) + 1])
			{
			OverlapFlgs |= cFeatBitIntrons;
			break;						// only need to find one intron overlap
			}
		else
			if(RelEndOfs < pGene->ExonStartEnds[(Idx*2) + 1])
				break;
		}
	}

return(OverlapFlgs);
}

// returns splice site overlaps
// Overlaps must be of at least Overlap bases from splice site into intron and exon to be characterised as an overlap
// 0:	no splice site overlaps
// cIntronExonSpliceSite	overlaps intron/exon on intron3'-5'exon splice site
// cExonIntronSpliceSite	overlaps exon/intron on exon3'-5'intron splice site
int
CBEDfile::GetFeatureBitsSpliceOverlaps(int FeatID,	// feature identifier
				   int StartOfs,			// start along cromosome 0..n
				   int EndOfs,				// end along chromosome start + len - 1
				   int Overlap)			    // must overlap by at least this many bases (1..n)
{
int OverlapFlgs;
int RelStartOfs;
int RelEndOfs;
int Idx;
tsGeneStructure *pGene;
// sanity checks!
#ifdef _DEBUG
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(0);

if(StartOfs < 0 || StartOfs > EndOfs)
	return(0);
if(FeatID < 1 || FeatID > m_FileHdr.NumFeatures || Overlap < 1)
	return(0);
#endif
// not interested in the trivial case where can't have an overlap
if(Overlap < 1 || ((EndOfs - StartOfs + 1) < (2 * Overlap)))
	return(0);
OverlapFlgs = 0;
tsBEDfeature *pProbe = m_ppFeatureChromStarts[FeatID-1];
if(pProbe->Score < m_MinScore || pProbe->Score > m_MaxScore)
	return(OverlapFlgs);
if((m_OnStrand != '*' && m_OnStrand != pProbe->Strand))
	return(OverlapFlgs);

// none of the remaining feature bit tests will succeed if no overlaps onto the gene/feature so may as well check
if(StartOfs > pProbe->End || EndOfs < pProbe->Start)
	return(0);

if(m_FileHdr.FeatType != eBTGeneExons)
	return(0);

// make start/end relative to gene start
if(StartOfs <= pProbe->Start)
	RelStartOfs = 0;
else
	RelStartOfs = StartOfs - pProbe->Start;
if(EndOfs > pProbe->End)
    RelEndOfs = pProbe->End - pProbe->Start;
else
	RelEndOfs = EndOfs - pProbe->Start;

pGene = (tsGeneStructure *)&pProbe->szName[pProbe->FeatNameLen+1];

if(pGene->NumExons > 1)				// must be at least 2 exons before an intron can exist
	{
	int ExonStart;
	int ExonEnd;
	for(Idx = 0; Idx < pGene->NumExons-1; Idx++)
		{
		ExonEnd = pGene->ExonStartEnds[(Idx * 2) + 1];
		ExonStart = pGene->ExonStartEnds[(Idx + 1) * 2];
		if(RelStartOfs <= (ExonEnd - (Overlap - 1)) &&
			RelEndOfs >= (ExonEnd + Overlap))
			OverlapFlgs |= pProbe->Strand == '+' ? cExonIntronSpliceSite : cIntronExonSpliceSite;
		if(RelStartOfs <= (ExonStart - Overlap) &&
			RelEndOfs >= (ExonStart + (Overlap - 1)))
			OverlapFlgs |= pProbe->Strand == '+' ? cIntronExonSpliceSite : cExonIntronSpliceSite;
		if(RelEndOfs < ExonStart)
			break;
		}
	}
return(OverlapFlgs);
}


bool 
CBEDfile::InInternFeat(int FeatBits,int ChromID,int StartOfs,int EndOfs)
{
int FeatID;
int FeatNxt;
int OverLaps;

// sanity checks only whilst debugging
#ifdef _DEBUG
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(false);
if(!ChromID || ChromID > m_FileHdr.NumChroms)
	return(false);
if(StartOfs < 0 || StartOfs > EndOfs)
	return(false);
if(!(FeatBits & cAnyFeatBits))
	return(false);
#endif

FeatNxt = 1;
while((FeatID = LocateFeatureIDinRangeOnChrom(ChromID,StartOfs,EndOfs,FeatNxt++,cFeatFiltIn,cFeatFiltOut))>0)
	{
	OverLaps = GetFeatureOverlaps(FeatBits,FeatID,StartOfs,EndOfs);
	if(OverLaps & FeatBits)
		return(true);
	}
return(false);
}

int
CBEDfile::GetFeatureBits(int ChromID,int StartOfs,int EndOfs,int FeatBits,int Updnstream)
{
int FeatID;
int FeatNxt;
int Overlaps;
int OverlapStartOfs;
int OverlapEndOfs;

// no point in determining if in up/dnstream when length of up/dnstream region in <= 0
if(Updnstream <= 0)
	FeatBits &= ~(cFeatBitUpstream | cFeatBitDnstream);
if(!(FeatBits & cAnyFeatBits))
	return(0);

if(FeatBits & (cFeatBitUpstream | cFeatBitDnstream))
	{
	OverlapStartOfs = StartOfs - Updnstream;
	if(OverlapStartOfs < 0)
		OverlapStartOfs = 0;
	OverlapEndOfs = EndOfs + Updnstream;
	}
else
	{
	OverlapStartOfs = StartOfs;
	OverlapEndOfs = EndOfs;
	}

FeatNxt = 1;
Overlaps = 0;
while((FeatID = LocateFeatureIDinRangeOnChrom(ChromID,OverlapStartOfs,OverlapEndOfs,FeatNxt++,cFeatFiltIn,cFeatFiltOut))>0)
	Overlaps |= GetFeatureOverlaps(FeatBits,FeatID,StartOfs,EndOfs,Updnstream);
return(Overlaps);
}

int
CBEDfile::GetSpliceSiteBits(int ChromID,int StartOfs,int EndOfs,int OverlapDistance)
{
int FeatID;
int FeatNxt;
int Overlaps;
int OverlapStartOfs;
int OverlapEndOfs;

// sanity checks only whilst debugging
#ifdef _DEBUG
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(0);

if(!ChromID || ChromID > m_FileHdr.NumChroms || OverlapDistance < 1 ||
	StartOfs < 0 || StartOfs > EndOfs)
	return(0);
#endif

OverlapStartOfs = StartOfs;
OverlapEndOfs = EndOfs;

FeatNxt = 1;
Overlaps = 0;
while((FeatID = LocateFeatureIDinRangeOnChrom(ChromID,OverlapStartOfs,OverlapEndOfs,FeatNxt++,cFeatFiltIn,cFeatFiltOut))>0)
	Overlaps |= GetFeatureBitsSpliceOverlaps(FeatID,StartOfs,EndOfs,OverlapDistance);
return(Overlaps);
}



// returns true if any point between StartOfs and EndOfs overlaps any CDS
bool 
CBEDfile::InAnyCDS(int ChromID,int StartOfs,int EndOfs)
{
return(InInternFeat(cFeatBitCDS,ChromID,StartOfs,EndOfs));
}

// returns true if any point between StartOfs and EndOfs overlaps any exon
bool 
CBEDfile::InAnyExon(int ChromID,int StartOfs,int EndOfs)
{
return(InInternFeat(cFeatBitExons,ChromID,StartOfs,EndOfs));
}

// returns true if any point between StartOfs and EndOfs overlaps any intron
bool 
CBEDfile::InAnyIntron(int ChromID,int StartOfs,int EndOfs)
{
return(InInternFeat(cFeatBitIntrons,ChromID,StartOfs,EndOfs));
}


// returns true if any point between StartOfs and EndOfs overlaps 5' UTR of any gene
bool CBEDfile::InAny5UTR(int ChromID,int StartOfs,int EndOfs)
{
return(InInternFeat(cFeatBit5UTR,ChromID,StartOfs,EndOfs));
}

// returns true if any point between StartOfs and EndOfs overlaps 3' UTR of any gene
bool CBEDfile::InAny3UTR(int ChromID,int StartOfs,int EndOfs)
{
return(InInternFeat(cFeatBit3UTR,ChromID,StartOfs,EndOfs));
}

// returns true if any point between StartOfs and EndOfs overlaps 5' upstream of gene within Distance of any gene start
bool CBEDfile::InAny5Upstream(int ChromID,int StartOfs,int EndOfs,int Distance)
{
int FeatID;
int CurIth;
int RelStartOfs;
int RelEndOfs;
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(false);
// sanity checks!
#ifdef _DEBUG
if(Distance < 1 || StartOfs < 0 || StartOfs > EndOfs)
	return(false);
if(ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(false);
#endif

if(StartOfs >= Distance)
	RelStartOfs = StartOfs - Distance;
else
	RelStartOfs = 0;
RelEndOfs = EndOfs + Distance;
CurIth = 1;
while((FeatID = LocateFeatureIDinRangeOnChrom(ChromID,RelStartOfs,RelEndOfs,CurIth++,cFeatFiltIn,cFeatFiltOut))>0)
	{
	if(In5Upstream(FeatID,StartOfs,EndOfs,Distance)==true)
		return(true);
	if(In3Dnstream(FeatID,StartOfs,EndOfs,Distance)==true)
		return(true);
	}
return(false);
}

// returns true if any point between StartOfs and EndOfs overlaps 3' downstream of any gene within Distance of gene end
bool CBEDfile::InAny3Dnstream(int ChromID,int StartOfs,int EndOfs,int Distance)
{
int FeatID;
int CurIth;
int RelStartOfs;
int RelEndOfs;
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(false);

// sanity checks!
#ifdef _DEBUG
if(Distance < 1 || StartOfs < 0 || StartOfs > EndOfs)
	return(false);

if(ChromID < 1 || ChromID > m_FileHdr.NumChroms)
	return(false);
#endif

if(StartOfs >= Distance)
	RelStartOfs = StartOfs - Distance;
else
	RelStartOfs = 0;
RelEndOfs = EndOfs + Distance;
CurIth = 1;

// try from gene starts
while((FeatID = LocateFeatureIDinRangeOnChrom(ChromID,RelStartOfs,RelEndOfs,CurIth++,cFeatFiltIn,cFeatFiltOut))>0)
	if(In3Dnstream(FeatID,StartOfs,EndOfs,Distance)==true)
		return(true);

// try from gene ends
CurIth = 1;
while((FeatID = LocateFeatureIDinRangeOnChrom(ChromID,RelStartOfs,RelEndOfs,CurIth++,cFeatFiltIn,cFeatFiltOut))>0)
	if(In3Dnstream(FeatID,StartOfs,EndOfs,Distance)==true)
		return(true);
return(false);
}

// returns true if any point between StartOfs and EndOfs overlaps specified gene
bool 
CBEDfile::InFeature(int FeatureID, int StartOfs,int EndOfs)		  
{
int OverLaps;
OverLaps = GetFeatureOverlaps(cFeatBitFeature,FeatureID,StartOfs,EndOfs);
return(OverLaps & cFeatBitFeature ? true : false);
}

// returns true if any point between StartOfs and EndOfs overlaps any CDS
bool 
CBEDfile::InCDS(int FeatureID, int StartOfs,int EndOfs)		  
{
int OverLaps;
OverLaps = GetFeatureOverlaps(cFeatBitCDS,FeatureID,StartOfs,EndOfs);
return(OverLaps & cFeatBitCDS ? true : false);
}

// returns true if any point between StartOfs and EndOfs overlaps any exon
bool 
CBEDfile::InExon(int FeatureID, int StartOfs,int EndOfs)		  
{
int OverLaps;
OverLaps = GetFeatureOverlaps(cFeatBitExons,FeatureID,StartOfs,EndOfs);
return(OverLaps & cFeatBitExons ? true : false);
}

// returns true if any point between StartOfs and EndOfs overlaps any intron	
bool 
CBEDfile::InIntron(int FeatureID, int StartOfs,int EndOfs)		  
{
int OverLaps;
OverLaps = GetFeatureOverlaps(cFeatBitIntrons,FeatureID,StartOfs,EndOfs);
return(OverLaps & cFeatBitIntrons ? true : false);
}

// returns true if any point between StartOfs and EndOfs overlaps 5' UTR of specified gene	
bool 
CBEDfile::In5UTR(int FeatureID, int StartOfs,int EndOfs)		  
{
int OverLaps;
OverLaps = GetFeatureOverlaps(cFeatBit5UTR,FeatureID,StartOfs,EndOfs);
return(OverLaps & cFeatBit5UTR ? true : false);
}

// returns true if any point between StartOfs and EndOfs 3' UTR of specified gene
bool CBEDfile::In3UTR(int FeatureID,int StartOfs,int EndOfs)		  
{
int OverLaps;
OverLaps = GetFeatureOverlaps(cFeatBit3UTR,FeatureID,StartOfs,EndOfs);
return(OverLaps & cFeatBit3UTR ? true : false);
}

// returns true if any point between StartOfs and EndOfs overlaps 5' upstream of gene within Distance
bool 
CBEDfile::In5Upstream(int FeatureID, int StartOfs,int EndOfs,int Distance) 
{
int OverLaps;
OverLaps = GetFeatureOverlaps(cFeatBitUpstream,FeatureID,StartOfs,EndOfs,Distance);
return(OverLaps & cFeatBitUpstream ? true : false);
}

// returns true if any point between StartOfs and EndOfs overlaps 3' downstream of gene within Distance 
bool 
CBEDfile::In3Dnstream(int FeatureID, int StartOfs,int EndOfs,int Distance) 
{
int OverLaps;
OverLaps = GetFeatureOverlaps(cFeatBitDnstream,FeatureID,StartOfs,EndOfs,Distance);
return(OverLaps & cFeatBitDnstream ? true : false);
}

// sets specified feature to have FiltFlags
// 
teBSFrsltCodes 
CBEDfile::SetFilter(char *pszFeatName,unsigned int FiltFlags)	
{
tsBEDfeature *pProbe;

bool bSet;
int FeatureID;
int Nth = 1;
if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
bSet = false;
while((FeatureID = LocateFeatureIDbyName(pszFeatName,Nth++)) > 0)
	{
	pProbe = m_ppFeatureChromStarts[FeatureID-1];
	pProbe->FiltFlags = FiltFlags;
	bSet = true;
	}
return(bSet ? eBSFSuccess : eBSFerrFeature);
}

// sets specified feature to have FiltFlags
teBSFrsltCodes 
CBEDfile::SetFilter(int FeatureID,unsigned int FiltFlags)	
{
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(eBSFerrFeature);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
pProbe->FiltFlags = FiltFlags;
return(eBSFSuccess);
}

teBSFrsltCodes 
CBEDfile::SetFilters(unsigned int FiltFlags)				// sets all features to have specified FiltFlags
{
int FeatureIdx;
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || m_FileHdr.NumFeatures < 1)
	return(eBSFerrFeature);
for(FeatureIdx = 0;FeatureIdx < m_FileHdr.NumFeatures; FeatureIdx++)
	{
	pProbe = m_ppFeatureChromStarts[FeatureIdx];
	pProbe->FiltFlags = FiltFlags;
	}
return(eBSFSuccess);
}

bool 
CBEDfile::Filter(int FeatureID,unsigned int FiltFlags)			// returns true if feature has any of specified FiltFlags set
{
tsBEDfeature *pProbe;

if(!m_bFeaturesAvail || !m_FileHdr.NumFeatures || FeatureID < 1 || FeatureID > m_FileHdr.NumFeatures)
	return(false);
pProbe = m_ppFeatureChromStarts[FeatureID-1];
return(pProbe->FiltFlags & FiltFlags ? true : false);
}

// returns highest priority functional region from feature bits
teFuncRegion 
CBEDfile::MapFeatureBits2Idx(int FeatureBits)
{
if(FeatureBits & cFeatBitCDS)	// CDS has the highest priority
	return(eFRCDS);
if(FeatureBits & cFeatBit5UTR)
	return(eFR5UTR);
if(FeatureBits & cFeatBit3UTR)
	return(eFR3UTR);
if(FeatureBits & cFeatBitIntrons)
	return(eFRIntronic);
if(FeatureBits & cFeatBitUpstream)
	return(eFRUpstream);
if(FeatureBits & cFeatBitDnstream)
	return(eFRDnstream);
return(eFRIntergenic);			// and intergenic the lowest
}


