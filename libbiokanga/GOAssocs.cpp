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


CGOAssocs::CGOAssocs(void)
{
m_pszGOIdents = NULL;		// memory allocated to hold conactenated '\0' separated GO:term names
m_pszGeneIdents = NULL;		// memory allocated to hold conactenated '\0' separated gene names
m_pGenes=NULL;				// memory allocated to hold tsGOGene structures
m_pAssocs=NULL;				// memory allocated to hold association structures
m_pGeneMappings = NULL;		// memory allocated for gene mappings
m_pGeneFilters = NULL;		// memory allocated for gene filters
m_pGeneGOitems = NULL;		// memory allocated for gene GO items
m_pszGeneGoItemsTxt = NULL;	// memory allocated for gene GO items text
m_NumGeneMaps = 0;			// number of gene name mappings
m_AllocdGeneMaps = 0;		// number of allocated gene name mappings
m_NumGeneFilters = 0;		// number of gene name filterings
m_AllocdGeneFilters = 0;	// number of allocated gene name filterings

m_hFile = -1;
Reset(false);
}

CGOAssocs::~CGOAssocs(void)
{
Reset(false);
}


void
CGOAssocs::Reset(bool bFlush)
{
if(bFlush)
	Flush2Disk();

if(m_hFile != -1)
	{
	close(m_hFile);
	m_hFile = -1;
	}

ClearAssociations();

m_NumGeneMaps = 0;			// number of gene name mappings
m_AllocdGeneMaps = 0;		// number of allocated gene name mappings
if(m_pGeneMappings != NULL)
	{
	delete m_pGeneMappings;
	m_pGeneMappings = NULL;
	}

m_NumGeneFilters = 0;			// number of gene Filters
m_AllocdGeneFilters = 0;		// number of allocated gene Filters
if(m_pGeneFilters != NULL)
	{
	delete m_pGeneFilters;
	m_pGeneFilters = NULL;
	}

m_NumGeneGOitems = 0;		// number of gene GO items
m_AllocdGeneGOitems = 0;	// number of allocated gene GO items
if(m_pGeneGOitems != NULL) // pts to gene GO items
	{
	delete m_pGeneGOitems;
	m_pGeneGOitems = NULL;
	}

m_NumGeneGOitemsTxt = 0;	// chrs used in m_pszGeneGoItemsTxt
m_AllocdGeneGoItemsTxt = 0;	// current allocation
if(m_pszGeneGoItemsTxt != NULL)	// allocated to hold all gene GO items text		
	{
	delete m_pszGeneGoItemsTxt;
	m_pszGeneGoItemsTxt = NULL;
	}

m_szFile[0] = '\0';				    // file name as opened/created
m_bCreate = false;					// TRUE if file opened for create or update 
m_bHdrDirty = false;				// TRUE if header needs to be written to file
m_DfltEvidence = eGOEnone;
m_bGeneNamesAsPtrs = false;
m_bDoGeneMaps = false;
m_bGOFlybase = false;
InitHdr();
}

teBSFrsltCodes
CGOAssocs::Disk2Hdr(char *pszGoFile,int FileType)
{
if(_lseeki64(m_hFile,0,SEEK_SET)!=0)			// read in header..
	{
	AddErrMsg("CGOAssocs::Open","Seek failed to offset 0 on %s - %s",pszGoFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

if(sizeof(tsGOAssocFileHdr) != read(m_hFile,&m_FileHdr,sizeof(tsGOAssocFileHdr)))
	{
	AddErrMsg("CGOAssocs::Open","Read of file header failed on %s - %s",pszGoFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

	// header read, validate it as being a GOAssocs file header
if(tolower(m_FileHdr.Magic[0]) != 'b' || 
	tolower(m_FileHdr.Magic[1]) != 'i' || 
	tolower(m_FileHdr.Magic[2]) != 'o' || 
	tolower(m_FileHdr.Magic[3]) != 's')
	{
	AddErrMsg("CGOAssocs::Open","%s opened but no magic signature - not a GOAssocs file",pszGoFile);
	Reset(false);			// closes opened file..
	return(eBSFerrNotBioseq);
	}

if(m_bIsBigEndian)	// file was written with little-endian ordering
	{
	m_FileHdr.FileLen = SwapUI64Endians(m_FileHdr.FileLen);					// current file length
	m_FileHdr.GOGenesOfs = SwapUI64Endians(m_FileHdr.GOGenesOfs);			// file offset to genes
	m_FileHdr.GOAssocsOfs = SwapUI64Endians(m_FileHdr.GOAssocsOfs);			// file offset to associations
	m_FileHdr.GeneIdentsOfs = SwapUI64Endians(m_FileHdr.GeneIdentsOfs);		// file offset to gene identifiers
	m_FileHdr.GOIdentsOfs = SwapUI64Endians(m_FileHdr.GOIdentsOfs);			// file offset to GO identifiers

	m_FileHdr.Type = SwapUI32Endians(m_FileHdr.Type);						// GOAssoc file type 
	m_FileHdr.Version = SwapUI32Endians(m_FileHdr.Version);					// header version, incremented if structure changes with later releases
	m_FileHdr.SizeOfHdr = SwapUI32Endians(m_FileHdr.SizeOfHdr);				// total size of this header

	m_FileHdr.GOGenesCnt = SwapUI32Endians(m_FileHdr.GOGenesCnt);			// number of genes
	m_FileHdr.GOGenesSize = SwapUI32Endians(m_FileHdr.GOGenesSize);			// size (bytes) genes

	m_FileHdr.GOAssocsCnt = SwapUI32Endians(m_FileHdr.GOAssocsCnt);			// number of associations
	m_FileHdr.GOAssocsSize = SwapUI32Endians(m_FileHdr.GOAssocsSize);		// size (bytes) associations

	m_FileHdr.GeneIdentsCnt = SwapUI32Endians(m_FileHdr.GeneIdentsCnt);		// number of gene identifiers
	m_FileHdr.GeneIdentsSize = SwapUI32Endians(m_FileHdr.GeneIdentsSize);	// size (bytes) gene identifiers
	
	m_FileHdr.GOIdentsCnt = SwapUI32Endians(m_FileHdr.GOIdentsCnt);			// number of GO identifiers
	m_FileHdr.GOIdentsSize = SwapUI32Endians(m_FileHdr.GOIdentsSize);		// size (bytes) GO identifiers
	}

	// check GOAssocs file is the type we are expecting
if(m_FileHdr.Type != FileType)
	{
	AddErrMsg("CGOAssocs::Open","%s opened as a GOAssocs file - expected type %d, file type is %d",pszGoFile,FileType,m_FileHdr.Type);
	Reset(false);			// closes opened file..
	return(eBSFerrFileType);
	}

	// can we handle this version?
if(m_FileHdr.Version > cBSGOAssocVersion || m_FileHdr.Version < cBSGOAssocVersionBack)
	{
	AddErrMsg("CGOAssocs::Open","%s opened as a GOAssocs file - can only handle versions %d to %d, file version is %d",pszGoFile,
			cBSGOAssocVersionBack,cBSGOAssocVersion,m_FileHdr.Version);
	Reset(false);			// closes opened file..
	return(eBSFerrFileVer);
	}
return(eBSFSuccess);
}

teBSFrsltCodes
CGOAssocs::Hdr2Disk(void)
{
tsGOAssocFileHdr FileHdr;
tsGOAssocFileHdr *pHdr;
int WrtLen;

WrtLen = sizeof(tsGOAssocFileHdr);

if(m_bIsBigEndian)	// if on a big-endian machine then need to make little endian as that is our native file format
	{
	memmove(&FileHdr,&m_FileHdr,WrtLen);
	FileHdr.FileLen = SwapUI64Endians(m_FileHdr.FileLen);					// current file length
	FileHdr.GOGenesOfs = SwapUI64Endians(m_FileHdr.GOGenesOfs);			// file offset to genes
	FileHdr.GOAssocsOfs = SwapUI64Endians(m_FileHdr.GOAssocsOfs);			// file offset to associations
	FileHdr.GeneIdentsOfs = SwapUI64Endians(m_FileHdr.GeneIdentsOfs);		// file offset to gene identifiers
	FileHdr.GOIdentsOfs = SwapUI64Endians(m_FileHdr.GOIdentsOfs);			// file offset to GO identifiers

	FileHdr.Type = SwapUI32Endians(m_FileHdr.Type);						// GOAssoc file type 
	FileHdr.Version = SwapUI32Endians(m_FileHdr.Version);					// header version, incremented if structure changes with later releases
	FileHdr.SizeOfHdr = SwapUI32Endians(m_FileHdr.SizeOfHdr);				// total size of this header

	FileHdr.GOGenesCnt = SwapUI32Endians(m_FileHdr.GOGenesCnt);			// number of genes
	FileHdr.GOGenesSize = SwapUI32Endians(m_FileHdr.GOGenesSize);			// size (bytes) genes

	FileHdr.GOAssocsCnt = SwapUI32Endians(m_FileHdr.GOAssocsCnt);			// number of associations
	FileHdr.GOAssocsSize = SwapUI32Endians(m_FileHdr.GOAssocsSize);		// size (bytes) associations

	FileHdr.GeneIdentsCnt = SwapUI32Endians(m_FileHdr.GeneIdentsCnt);		// number of gene identifiers
	FileHdr.GeneIdentsSize = SwapUI32Endians(m_FileHdr.GeneIdentsSize);	// size (bytes) gene identifiers
	
	FileHdr.GOIdentsCnt = SwapUI32Endians(m_FileHdr.GOIdentsCnt);			// number of GO identifiers
	FileHdr.GOIdentsSize = SwapUI32Endians(m_FileHdr.GOIdentsSize);		// size (bytes) GO identifiers
	pHdr = &FileHdr;
	}
else
	pHdr = &m_FileHdr;

if(_lseeki64(m_hFile,0,SEEK_SET) ||
			write(m_hFile,pHdr,WrtLen)!=WrtLen)
	{
	AddErrMsg("CGOAssocs::Hdr2Disk","Unable to write file header to disk on file %s - error %s",m_szFile,strerror(errno));
	Reset(false);
	return(eBSFerrFileAccess);
	}
m_bHdrDirty = false;
return(eBSFSuccess);
}



void
CGOAssocs::ClearAssociations(void)
{
if(m_pszGOIdents != NULL)
	{
	delete m_pszGOIdents;
	m_pszGOIdents = NULL;
	}

if(m_pszGeneIdents != NULL)
	{
	delete m_pszGeneIdents;
	m_pszGeneIdents = NULL;
	}

if(m_pGenes != NULL)
	{
	delete m_pGenes;
	m_pGenes = NULL;
	}

if(m_pAssocs != NULL)
	{
	delete m_pAssocs;
	m_pAssocs = NULL;
	}

m_GOIdentsCnt=0;		// number of GO identifiers
m_NxtGOIdentOfs=0;		// offset into m_pszGOIdents where to next write GO:term identifiers
m_MaxAllocGOIdents=0;	// size of memory currently allocated to m_pszGOIdents

m_GeneIdentsCnt = 0;	// current number of gene name identifiers
m_NxtGeneIdentOfs=0;	// offset into m_pszGeneIdents where to next write gene names
m_MaxAllocGeneIdents=0;	// size of memory currently allocated to m_pszGeneIdents

m_GOGenesCnt=0;			// current number of genes
m_MaxAllocGOGene=0;		// max number of tsGOGene that m_pGenes can hold

m_GOAssocCnt=0;			// current number of associations
m_MaxAllocGOAssoc=0;	// max number of tsGOAssoc that m_pAssocs can hold
}


void 
CGOAssocs::InitHdr(void)
{
memset(&m_FileHdr,0,sizeof(m_FileHdr));
m_FileHdr.Magic[0] = 'b';
m_FileHdr.Magic[1] = 'i';
m_FileHdr.Magic[2] = 'o';
m_FileHdr.Magic[3] = 's';
m_FileHdr.Type = cBSFTypeGOAssoc;		// biosequence file type 
m_FileHdr.Version = cBSGOAssocVersion;	// header version, incremented if structure changes with later releases
m_FileHdr.FileLen = sizeof(m_FileHdr);	// current file length
m_FileHdr.SizeOfHdr = sizeof(m_FileHdr);// total size of this header
m_FileHdr.szDescription[0] = '\0';
m_FileHdr.szTitle[0] = '\0';
m_bHdrDirty = true;						// TRUE if header needs to be written to file
}

teBSFrsltCodes 
CGOAssocs::Open(char *pszFileName,bool bCreate)
{
teBSFrsltCodes Rslt;
if(pszFileName == NULL || *pszFileName == '\0') // validate parameters
	{
	AddErrMsg("CGOAssocs::Open","Parameter errors");
	return(eBSFerrParams);
	}
Reset(false);						// reset context in case file was previously opened

#ifdef _WIN32
if(!bCreate)
	m_hFile = open(pszFileName, O_READSEQ ); // file access is normally sequential..
else
	m_hFile = open(pszFileName, O_CREATETRUNC);
#else
if(!bCreate)
m_hFile = open64(pszFileName, O_READSEQ); // file access is normally sequential..
else
	{
     if((m_hFile = open64(pszFileName,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	   if(ftruncate(m_hFile,0)!=0)
			{
			AddErrMsg("Open","Unable to truncate %s - %s",pszFileName,strerror(errno));
			return(eBSFerrCreateFile);
			}
	}
#endif

if(m_hFile == -1)					// check if file open succeeded
	{
	AddErrMsg("CGOAssocs::Open","Unable to open %s - %s",pszFileName,strerror(errno));
	Rslt = bCreate ? eBSFerrCreateFile : eBSFerrOpnFile;
	return(Rslt);
	}

strncpy(m_szFile,pszFileName,_MAX_PATH);
m_szFile[_MAX_PATH-1] = '\0';
if(bCreate)
	{
	m_bCreate = true;
	m_bGOAssocAvail = false;
	InitHdr();
	if((Rslt = Flush2Disk()) != eBSFSuccess)
		{
		Reset(false);			// closes opened file..
		return(Rslt);
		}
	}
else // else opening existing file
	{
	if((Rslt=Disk2Hdr(pszFileName,cBSFTypeGOAssoc))!=eBSFSuccess)
		{
		Reset(false);			// closes opened file..
		return((teBSFrsltCodes)Rslt);
		}

	// if not empty then load genes and associations...
	if(m_FileHdr.GeneIdentsCnt)
		{
		if((Rslt=LoadAssociations())!=eBSFSuccess)	// close file after associations loaded
			{
			AddErrMsg("CGOAssocs::Open","Error loading GO associations from %s",pszFileName);
			Reset(false);
			return(Rslt);
			}
		}
	m_bGOAssocAvail = true;		// associations are in memory and can now be accessed
	}
return(eBSFSuccess);
}

teBSFrsltCodes 
CGOAssocs::Close(bool bFlush2Disk)
{
if(bFlush2Disk) {
	if(m_GOGenesCnt > 1 && m_pGenes != NULL) {
		SwitchGeneNameIDX2PTR();
		qsort(m_pGenes,m_GOGenesCnt,sizeof(tsGOGene),SortGenes);
		SwitchGeneNamePTR2IDX();
		}
	Flush2Disk();
	}
Reset();
return(eBSFSuccess);
}

teBSFrsltCodes 
CGOAssocs::Flush2Disk(void)
{
int Rslt;
int WrtLen;
int Idx;

if(m_hFile != -1 && m_bCreate)	// if file opened for write because creating..
	{
	if(m_GOGenesCnt)			// are there genes and associations to write?
		{
		// ensure genes are sorted by identifier
		SwitchGeneNameIDX2PTR();
		qsort(m_pGenes,m_GOGenesCnt,sizeof(tsGOGene),SortGenes);
		SwitchGeneNamePTR2IDX();

		// write out m_GOIdentsCnt GO identifiers
		m_FileHdr.GOIdentsCnt = m_GOIdentsCnt;
		if(m_GOIdentsCnt)
			{
			m_FileHdr.GOIdentsSize = m_NxtGOIdentOfs;
			m_FileHdr.GOIdentsOfs = m_FileHdr.FileLen;
			WrtLen = m_FileHdr.GOIdentsSize;	
			if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
				write(m_hFile,(char *)m_pszGOIdents,WrtLen)!=WrtLen)
				{
				AddErrMsg("CGOTerms::Flush2Disk","Unable to write GO identifiers to disk on file %s - error %s",m_szFile,strerror(errno));
				Reset(false);
				return(eBSFerrFileAccess);
				}
			m_FileHdr.FileLen += WrtLen;
			}
		else
			{
			m_FileHdr.GOIdentsSize = 0;
			m_FileHdr.GOIdentsOfs = 0;
			}

		// write out m_GeneIdentsCnt identifiers
		m_FileHdr.GeneIdentsCnt = m_GeneIdentsCnt;
		if(m_GeneIdentsCnt)
			{
			m_FileHdr.GeneIdentsSize = m_NxtGeneIdentOfs;
			m_FileHdr.GeneIdentsOfs = m_FileHdr.FileLen;
			WrtLen = m_FileHdr.GeneIdentsSize;	
			if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
				write(m_hFile,(char *)m_pszGeneIdents,WrtLen)!=WrtLen)
				{
				AddErrMsg("CGOTerms::Flush2Disk","Unable to write GO identifiers to disk on file %s - error %s",m_szFile,strerror(errno));
				Reset(false);
				return(eBSFerrFileAccess);
				}
			m_FileHdr.FileLen += WrtLen;
			}
		else
			{
			m_FileHdr.GeneIdentsSize = 0;
			m_FileHdr.GeneIdentsOfs = 0;
			}



		// write out gene identifiers
		m_FileHdr.GOGenesCnt = m_GOGenesCnt;
		m_FileHdr.GOGenesSize = sizeof(tsGOGene) * m_GOGenesCnt;
		m_FileHdr.GOGenesOfs = m_FileHdr.FileLen;
		WrtLen = m_FileHdr.GOGenesSize;	

		if(m_bIsBigEndian)
			{
			tsGOGene *pGOGene = m_pGenes;
			for(Idx = 0; Idx < m_GOGenesCnt; Idx++,pGOGene++)
				{
				pGOGene->Name.Pad64 = SwapUI64Endians(pGOGene->Name.Pad64);			// gene name
				pGOGene->GeneID = SwapUI32Endians(pGOGene->GeneID);					// unique identifier for this gene or transcriptional loci
				pGOGene->NumGOAssoc = SwapUI32Endians(pGOGene->NumGOAssoc);			// number of associated GO terms
				pGOGene->GOAssocIdx = SwapUI32Endians(pGOGene->GOAssocIdx);			// index at which first GOAssoc (tsGOAssoc) for this gene starts
				}
			}


		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
			write(m_hFile,(char *)m_pGenes,WrtLen)!=WrtLen)
			{
			AddErrMsg("CGOTerms::Flush2Disk","Unable to write GO identifiers to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		m_FileHdr.FileLen += WrtLen;


		// write out m_GOAssocCnt associations
		m_FileHdr.GOAssocsCnt = m_GOAssocCnt;
		if(m_GOAssocCnt)
			{
			m_FileHdr.GOAssocsSize = sizeof(tsGOAssoc) * m_GOAssocCnt;
			m_FileHdr.GOAssocsOfs = m_FileHdr.FileLen;
			WrtLen = m_FileHdr.GOAssocsSize;	

			if(m_bIsBigEndian)
				{
				tsGOAssoc *pGOAssoc = m_pAssocs;
				for(Idx = 0; Idx < m_GOAssocCnt; Idx++,pGOAssoc++)
					{
					pGOAssoc->GeneID = SwapUI32Endians(pGOAssoc->GeneID);				// associated with which gene
					pGOAssoc->IdentIdx = SwapUI32Endians(pGOAssoc->IdentIdx);			// index at which associated GO:term identifier starts
					pGOAssoc->Evidence = SwapUI16Endians(pGOAssoc->Evidence);			// eGOEvidence
					}
				}

			if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
				write(m_hFile,(char *)m_pAssocs,WrtLen)!=WrtLen)
				{
				AddErrMsg("CGOTerms::Flush2Disk","Unable to write GO associations to disk on file %s - error %s",m_szFile,strerror(errno));
				Reset(false);
				return(eBSFerrFileAccess);
				}
			m_FileHdr.FileLen += WrtLen;
			}
		else
			{
			m_FileHdr.GOAssocsSize = 0;
			m_FileHdr.GOAssocsOfs = 0;
			}
		}

		// now write the header to disk
	if((Rslt = Hdr2Disk())!=eBSFSuccess)
		{
		Reset(false);
		return((teBSFrsltCodes)Rslt);
		}
	m_bHdrDirty = false;
	}
return(eBSFSuccess);
}


// LoadAssociations
// Load GO associations from previously opened file
// Any already loaded associations are deleted
teBSFrsltCodes
CGOAssocs::LoadAssociations(void)
{
teBSFrsltCodes Rslt;
int Idx;

if(m_hFile == -1)
	return(eBSFerrClosed);


ClearAssociations();

m_bGeneNamesAsPtrs=false;	// TRUE if TermGOID currently set as ptrs instead of idx into m_pGOIDs (required when sorting) 

if(!m_FileHdr.GOGenesCnt)	// any genes and associations in file to load?
	return(eBSFerrNoFeatures);

// allocate all memory up front before loading from file
if(m_FileHdr.GOIdentsCnt)
	{
	if((m_pszGOIdents = (char *)new unsigned char [m_FileHdr.GOIdentsSize])==NULL)
		return(eBSFerrMem);
	}
m_GOIdentsCnt =	m_FileHdr.GOIdentsCnt;
m_NxtGOIdentOfs = m_FileHdr.GOIdentsSize;
m_MaxAllocGOIdents = m_FileHdr.GOIdentsSize;

if(m_FileHdr.GeneIdentsCnt)
	{
	if((m_pszGeneIdents = (char *)new unsigned char [m_FileHdr.GeneIdentsSize])==NULL)
		return(eBSFerrMem);
	}
m_GeneIdentsCnt =	m_FileHdr.GeneIdentsCnt;
m_NxtGeneIdentOfs = m_FileHdr.GeneIdentsSize;
m_MaxAllocGeneIdents = m_FileHdr.GeneIdentsSize;

if(m_FileHdr.GOGenesCnt)
	{
	if((m_pGenes = (tsGOGene *)new unsigned char [m_FileHdr.GOGenesSize])==NULL)
		return(eBSFerrMem);
	}
m_GOGenesCnt =	m_FileHdr.GOGenesCnt;
m_MaxAllocGOGene = m_FileHdr.GOGenesCnt;

if(m_FileHdr.GOAssocsCnt)
	{
	if((m_pAssocs = (tsGOAssoc *)new unsigned char [m_FileHdr.GOAssocsSize])==NULL)
		return(eBSFerrMem);
	}
m_GOAssocCnt =	m_FileHdr.GOAssocsCnt;
m_MaxAllocGOAssoc = m_FileHdr.GOAssocsCnt;

// all required memory has been allocated, now read from disk
if(m_pszGOIdents != NULL)
	Rslt = ReadDisk(m_FileHdr.GOIdentsOfs,m_FileHdr.GOIdentsSize,m_pszGOIdents);
	
if(m_pszGeneIdents != NULL && Rslt >= eBSFSuccess)
	Rslt = ReadDisk(m_FileHdr.GeneIdentsOfs,m_FileHdr.GeneIdentsSize,m_pszGeneIdents);

	
if(m_pGenes != NULL && Rslt >= eBSFSuccess)	
	{
	Rslt = ReadDisk(m_FileHdr.GOGenesOfs,m_FileHdr.GOGenesSize,m_pGenes);
	if(Rslt == eBSFSuccess && m_bIsBigEndian)
		{
		tsGOGene *pGOGene = m_pGenes;
		for(Idx = 0; Idx < m_GOGenesCnt; Idx++,pGOGene++)
			{
			pGOGene->Name.Pad64 = SwapUI64Endians(pGOGene->Name.Pad64);			// gene name
			pGOGene->GeneID = SwapUI32Endians(pGOGene->GeneID);					// unique identifier for this gene or transcriptional loci
			pGOGene->NumGOAssoc = SwapUI32Endians(pGOGene->NumGOAssoc);			// number of associated GO terms
			pGOGene->GOAssocIdx = SwapUI32Endians(pGOGene->GOAssocIdx);			// index at which first GOAssoc (tsGOAssoc) for this gene starts
			}
		}
	}

	
if(m_pAssocs != NULL && Rslt >= eBSFSuccess)
	{
	Rslt = ReadDisk(m_FileHdr.GOAssocsOfs,m_FileHdr.GOAssocsSize,m_pAssocs);
	if(Rslt == eBSFSuccess && m_bIsBigEndian)
		{
		tsGOAssoc *pGOAssoc = m_pAssocs;
		for(Idx = 0; Idx < m_GOAssocCnt; Idx++,pGOAssoc++)
			{
			pGOAssoc->GeneID = SwapUI32Endians(pGOAssoc->GeneID);				// associated with which gene
			pGOAssoc->IdentIdx = SwapUI32Endians(pGOAssoc->IdentIdx);			// index at which associated GO:term identifier starts
			pGOAssoc->Evidence = SwapUI16Endians(pGOAssoc->Evidence);			// eGOEvidence
			}
		}
	}

if(Rslt != eBSFSuccess)
	ClearAssociations();
return(Rslt);
}

// ReadDisk
// Reads block of size 'Len' from disk starting at 'DiskOfs' into preallocated memory at 'pTo'
teBSFrsltCodes
CGOAssocs::ReadDisk(INT64 DiskOfs,int Len,void *pTo)
{
if(_lseeki64(m_hFile,DiskOfs,SEEK_SET)!=DiskOfs)
	{
	AddErrMsg("CGOTerms::ReadDisk","Seek failed on %s - %s",m_szFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}
if(read(m_hFile,pTo,Len)!=Len)
	{
	AddErrMsg("CGOTerms::ReadDisk","Read failed on %s - %s",m_szFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}
return(eBSFSuccess);
}

int 
CGOAssocs::Parse(etGOAssocParseType Type,  // expected file type
				 char *pszAssocFile,    	// GO association file to parse
				 char *pszGeneMapFile,		// optional gene map file
				 char *pszGeneFiltFile)		// optional gene filter file
{
int Rslt;
switch(Type) {
	case eGOAPMUCSU:
		Rslt = ParseUCSC(pszAssocFile,pszGeneFiltFile);
		break;
	case eGOAPGO: case eGOAPGOTAIR: case eGOAPGOFB:
		Rslt = ParseGOAnnotation(Type,pszAssocFile,pszGeneMapFile,pszGeneFiltFile);
		break;
	default:
		return(eBSFerrParams);
	}
GenGONumGenes();
return(Rslt);
}

// ParseGOAnnotation
// Expected geneontology.com annotation file format (tab delimited)
//<DB><tab><DB_Object_ID><tab><DB_Object_Symbol><tab><Qualifier><tab><GOID><tab><DBReference><tab>
// <Evidence><tab><WithOrFrom><tab><Aspect><tab><DB_Object_Name><tab><DB_Object_Synonym><tab>
// <DB_Object_Type><tab><taxon><tab><Date><tab><Assigned_by>
// NOTE: All fields are mandatory with the exception of <Qualifier>,<WithOrFrom>,<DB_Object_Name> and <DB_Object_Synonym>
int
CGOAssocs::ParseGOAnnotation(etGOAssocParseType Type,  // expected file type
							 char *pszGOAnnotation,	// GO annotation file to parse
								char *pszGeneMap,	// optional gene mapping file
								char *pszGeneFilters) // optional gene filters file
{
FILE *pAssocStream;
int Len;
int LineNum;
int Idx;
char szLineBuff[cLineBuffLen];
char *pTxt;
char szDB[cMaxLenGO_DB];
char szDB_Object_ID[cMaxGOFieldLen];
char szDB_Object_Symbol[cMaxGOFieldLen];
char szQualifier[cMaxGOFieldLen];
char szGOID[cMaxGOFieldLen];
char szDBReference[cMaxGOFieldLen];
char szEvidence[cMaxGOFieldLen];
char szWithOrFrom[cMaxGOFieldLen];
char szAspect[cMaxGOFieldLen];
char szDB_Object_Name[cMaxGOFieldLen];
char szDB_Object_Synonym[cMaxGOFieldLen];
char szDB_Object_Type[cMaxGOFieldLen];
char szTaxon[cMaxGOFieldLen];
char szDate[cMaxGOFieldLen];
char szAssigned_by[cMaxGOFieldLen];



int UnmappedGenes = 0;
int MappedGenes = 0;

int Rslt = eBSFSuccess;

m_bGOFlybase = Type == eGOAPGOFB ? true : false;
m_FileParseType = Type;

if(pszGeneFilters != NULL && pszGeneFilters[0] != '\0')
	{
	if((Idx = ParseGeneFilters(pszGeneFilters)) < 0)
		return(Idx);
	if(!Idx)
		{
		AddErrMsg("CGOAssocs::Open","Gene filters file parsed but no genes to filter! %s",pszGeneFilters);
		return(eBSFerrOpnFile);
		}
	m_bDoGeneFilters = true;
	}
else
	m_bDoGeneFilters = false;


if((m_bGOFlybase && m_bDoGeneFilters) || pszGeneMap != NULL && pszGeneMap[0] != '\0')
	{
	if((Idx = ParseGeneMapping(pszGeneMap)) < 0)
		return(Idx);
	if(!Idx && (pszGeneMap != NULL && pszGeneMap[0] != '\0'))
		{
		AddErrMsg("CGOAssocs::Open","Gene mappings file parsed but no mappings! %s",pszGeneMap);
		return(eBSFerrOpnFile);
		}
	m_bDoGeneMaps = true;
	}
else
	m_bDoGeneMaps = false;

if((pAssocStream = fopen(pszGOAnnotation,"r"))==NULL)
	{
	AddErrMsg("CGOAssocs::Open","Unable to open GO annotation file %s error: %s",pszGOAnnotation,strerror(errno));
	return(eBSFerrOpnFile);
	}

LineNum = 0;

while(fgets(szLineBuff,sizeof(szLineBuff),pAssocStream)!= NULL)
	{
	LineNum++;
	// simply slough lines which were just whitespace or start with '!' as comment line
	Len = TrimWhitespace(szLineBuff);
	if(!Len || szLineBuff[0] == '!')	
		continue;

	szDB[0] = '\0';
	szDB_Object_ID[0] = '\0';
	szDB_Object_Symbol[0] = '\0';
	szQualifier[0] = '\0';
	szGOID[0] = '\0';
	szDBReference[0] = '\0';
	szEvidence[0] = '\0';
	szWithOrFrom[0] = '\0';
	szAspect[0] = '\0';
	szDB_Object_Name[0] = '\0';
	szDB_Object_Synonym[0] = '\0';
	szDB_Object_Type[0] = '\0';
	szTaxon[0] = '\0';
	szDate[0] = '\0';
	szAssigned_by[0] = '\0';

	// parse out each field individually as some are optional
	pTxt = szLineBuff;
	Idx = ParseNxtField(pTxt,true,sizeof(szDB)-1,szDB);
	if(Idx > 0)
		{
		pTxt += Idx;
		Idx = ParseNxtField(pTxt,true,sizeof(szDB_Object_ID)-1,szDB_Object_ID);
		if(Idx > 0)
			{
			pTxt += Idx;
			Idx = ParseNxtField(pTxt,true,sizeof(szDB_Object_Symbol)-1,szDB_Object_Symbol);
			if(Idx > 0)
				{
				pTxt += Idx;
				Idx = ParseNxtField(pTxt,false,sizeof(szQualifier)-1,szQualifier);
				if(Idx > 0)
					{
					pTxt += Idx;
					Idx = ParseNxtField(pTxt,true,sizeof(szGOID)-1,szGOID);
					if(Idx > 0)
						{
						pTxt += Idx;
						Idx = ParseNxtField(pTxt,true,sizeof(szDBReference)-1,szDBReference);
						if(Idx > 0)
							{
							pTxt += Idx;
							Idx = ParseNxtField(pTxt,true,sizeof(szEvidence)-1,szEvidence);
							if(Idx > 0)
								{
								pTxt += Idx;
								Idx = ParseNxtField(pTxt,false,sizeof(szWithOrFrom)-1,szWithOrFrom);
								if(Idx > 0)
									{
									pTxt += Idx;
									Idx = ParseNxtField(pTxt,true,sizeof(szAspect)-1,szAspect);
									if(Idx > 0)
										{
										pTxt += Idx;	
										Idx = ParseNxtField(pTxt,false,sizeof(szDB_Object_Name)-1,szDB_Object_Name);
										if(Idx > 0)
											{
											pTxt += Idx;	
											Idx = ParseNxtField(pTxt,false,sizeof(szDB_Object_Synonym)-1,szDB_Object_Synonym);
											if(Idx > 0)
												{
												pTxt += Idx;	
												Idx = ParseNxtField(pTxt,true,sizeof(szDB_Object_Type)-1,szDB_Object_Type);
												if(Idx > 0)
													{
													pTxt += Idx;
													Idx = ParseNxtField(pTxt,true,sizeof(szTaxon)-1,szTaxon);
													if(Idx > 0)
														{
														pTxt += Idx;
														Idx = ParseNxtField(pTxt,true,sizeof(szDate)-1,szDate);
														if(Idx > 0)
															{
															pTxt += Idx;
															Idx = ParseNxtField(pTxt,true,sizeof(szAssigned_by)-1,szAssigned_by);
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	
	// Idx will contain last parse error if any were detected 
	if(Idx < 0)
		{
		AddErrMsg("CGOAssocs::ParseGOAnnotation","GO association file line %d parse error",LineNum);
		fclose(pAssocStream);
		return(eBSFerrParse);
		}

	if(m_FileParseType == eGOAPGOTAIR)	// For TAIR GO associations, need to trim back the suffix if suffix is of the form '.[-9]'
		{
		TrimNameIso(szDB_Object_ID);
		TrimNameIso(szDB_Object_Symbol);
		TrimNameIso(szDB_Object_Name);
		}
	// no parse errors, can provisionally add to gene GO list 
	if((Rslt = ProvAddGeneGOlist(m_bDoGeneMaps,szGOID,szEvidence,szQualifier,szAspect,szDB_Object_Type,szDB_Object_ID,szDB_Object_Symbol,szDB_Object_Name,szDB_Object_Synonym)) < 0)
		{
		AddErrMsg("CGOAssocs::ParseGOAnnotation","GO association file line %d ProvAddGeneGOlist() error",LineNum);
		fclose(pAssocStream);
		return(Rslt);
		}
	}

fclose(pAssocStream);
ProcessGeneGOlist();
return(Rslt);
}

// TrimNameIso
// Inplace remove any name isoform suffix of the form '.[0-99]'
bool			// true if isoform suffix was trimmed
CGOAssocs::TrimNameIso(char *pszName)
{
char *pSfx;
int NameLen;
if(pszName == NULL || pszName[0] == '\0')
	return(false);
NameLen = (int)strlen(pszName);
if(NameLen < cMinIsonameLen)
	return(false);
pSfx = &pszName[NameLen-1];
if(*pSfx >= '0' && *pSfx <= '9')
	{
	pSfx -= 1;
	if(*pSfx >= '0' && *pSfx <= '9')
		pSfx -= 1;
	if(*pSfx == '.')
		{
		*pSfx = '\0';
		return(true);
		}
	}
return(false);
}

int
CGOAssocs::ProvAddGeneGOlist(bool bDoGeneMaps,
							 char *pszGOID,
							 char *pszEvidence,
							 char *pszQualifier,
							 char *pszAspect,
							 char *pszObjectType,
							 char *pszObjectID,
							 char *pszObjectSymbol,
							 char *pszObjectName,
							 char *pszSynonyms)
{
int nThInstance;
int Rslt;
tsGeneMap *pMappedGene;
int Idx;
char szSynonym[cMaxGOFieldLen];
int Cnt = 0;

if(!bDoGeneMaps)	// if not requested to map GO gene names
	{
	if((Rslt = AddGeneGOlist(MapGOobjType(pszObjectType),pszObjectID,pszQualifier,pszAspect,pszEvidence,pszGOID)) < 0)
		return(Rslt);
	if((Rslt = AddGeneGOlist(MapGOobjType(pszObjectType),pszObjectSymbol,pszQualifier,pszAspect,pszEvidence,pszGOID)) < 0)
		return(Rslt);
	if(pszObjectName != NULL && pszObjectName[0])
		if((Rslt = AddGeneGOlist(MapGOobjType(pszObjectType),pszObjectName,pszQualifier,pszAspect,pszEvidence,pszGOID)) < 0)
			return(Rslt);
	if(pszSynonyms != NULL && pszSynonyms[0])
		{
		while((Idx = ParseNxtSynonym(pszSynonyms,szSynonym)) > 0)
			{
			if(m_FileParseType == eGOAPGOTAIR)	// For TAIR GO associations, need to trim back the suffix if suffix is of the form '.[-9]'
				TrimNameIso(szSynonym);
			if((Rslt = AddGeneGOlist(MapGOobjType(pszObjectType),szSynonym,pszQualifier,pszAspect,pszEvidence,pszGOID)) < 0)
					return(Rslt);
			pszSynonyms += Idx;
			}
		}
	return(1);
	}

// need to map the GO gene names onto user specified gene names
nThInstance = 1;
while((pMappedGene = LocateGeneMap(pszObjectID,nThInstance++))!=NULL)
	{
	if((Rslt = AddGeneGOlist(MapGOobjType(pszObjectType),pMappedGene->szToGene,pszQualifier,pszAspect,pszEvidence,pszGOID)) < 0)
		return(Rslt);
	Cnt += 1;
	}

nThInstance = 1;
while((pMappedGene = LocateGeneMap(pszObjectSymbol,nThInstance++))!=NULL)
	{
	if((Rslt = AddGeneGOlist(MapGOobjType(pszObjectType),pMappedGene->szToGene,pszQualifier,pszAspect,pszEvidence,pszGOID)) < 0)
		return(Rslt);
	Cnt += 1;
	}

if(pszObjectName != NULL && pszObjectName[0])
	{
	nThInstance = 1;
	while((pMappedGene = LocateGeneMap(pszObjectName,nThInstance++))!=NULL)
		{
		if((Rslt = AddGeneGOlist(MapGOobjType(pszObjectType),pMappedGene->szToGene,pszQualifier,pszAspect,pszEvidence,pszGOID)) < 0)
			return(Rslt);
		Cnt += 1;
		}
	}

if(pszSynonyms != NULL && pszSynonyms[0])
	{
	while((Idx = ParseNxtSynonym(pszSynonyms,szSynonym)) > 0)
		{
		if(m_FileParseType == eGOAPGOTAIR)	// For TAIR GO associations, need to trim back the suffix if suffix is of the form '.[-9]'
			TrimNameIso(szSynonym);
		nThInstance = 1;
		while((pMappedGene = LocateGeneMap(szSynonym,nThInstance++))!=NULL)
			{
			if((Rslt = AddGeneGOlist(MapGOobjType(pszObjectType),pMappedGene->szToGene,pszQualifier,pszAspect,pszEvidence,pszGOID)) < 0)
				return(Rslt);
			Cnt += 1;
			}
		pszSynonyms += Idx;
		}
	}

return(Cnt);
}

int
CGOAssocs::ParseNxtSynonym(char *pszToParse,char *pszRetField)
{
char *pszField = pszRetField;
char Chr;
int Idx;
int SymLen;

if(pszRetField == NULL || pszToParse == NULL)
	return(-1);

pszField[0] = '\0';
if(*pszToParse == '\0')
	return(0);

Idx = 0;
SymLen = 0;
pszField = pszRetField;
while((Chr = (*pszField++ = *pszToParse++)) != '\0') 
	{
	Idx++;
	if(Chr == ' ' && !SymLen)		// need to trim any prefix spaces
		{
		pszField -= 1;
		continue;
		}

	if(Chr == '|')
		{
		if(!SymLen)					// slough empty symbols
			{
			pszField -= 1;
			continue;
			}
		pszField[-1] = '\0';
		break;
		}

	SymLen += 1;
	}

if(SymLen)
	{
	pszField -= 2;
	while(*pszField == ' ')		// trim any suffix spaces
		pszField -= 1;
	pszField[1] = '\0';
	return(Idx);
	}
return(0);	
}

int 
CGOAssocs::AddGeneGOlist(etGOobjType ObjType,
   			   const char *pszGene,			// gene name
			   const char *pszQuals,			// list of GO qualifiers
			   const char *pszAspects,		// list of GO aspects
			   const char *pszEvidences,		// list of GO evidence
			   const char *pszGOIDs)			// list of GO:id
{
int NameLen;
int QualLen;
int AspectsLen;
int EvidencesLen;
int GOIDsLen;
int ReqLen;

tsGeneGOitem *pTmp;
char *pChr;

if(m_bDoGeneFilters && LocateGeneFilter(pszGene) == NULL)
	return(0);

if(m_pGeneGOitems == NULL || m_NumGeneGOitems == m_AllocdGeneGOitems)
	{
	if((pTmp = new tsGeneGOitem [m_AllocdGeneGOitems+10000])==NULL)
		{
		AddErrMsg("CGOAssocs::AddGeneGOlist","Unable to allocate memory");
		return(eBSFerrMem);
		}

	if(m_pGeneGOitems != NULL)
		{
		if(m_NumGeneGOitems)
			memmove(pTmp,m_pGeneGOitems,sizeof(tsGeneGOitem) * m_NumGeneGOitems);
		delete m_pGeneGOitems;
		}
	else
		{
		m_AllocdGeneGOitems = 0;
		m_NumGeneGOitems = 0;
		}
	m_pGeneGOitems = pTmp;
	m_AllocdGeneGOitems += 10000;
	}

NameLen = (int)strlen(pszGene) + 1;			// gene name
QualLen = (int)strlen(pszQuals) + 1;		// list of GO qualifiers
AspectsLen=(int)strlen(pszAspects) + 1;		// list of GO aspects
EvidencesLen=(int)strlen(pszEvidences) + 1;	// list of GO evidence
GOIDsLen=(int)strlen(pszGOIDs) + 1;			// list of GO:id
ReqLen = NameLen + QualLen + AspectsLen + EvidencesLen + GOIDsLen;

if(m_pszGeneGoItemsTxt == NULL || m_NumGeneGOitemsTxt + ReqLen >= m_AllocdGeneGoItemsTxt)
	{
	if((pChr = new char [m_AllocdGeneGoItemsTxt+100000000])==NULL)  // this will one day bite hard, reallocating but not updating all ptrs into m_pszGeneGoItemsTxt
		{
		AddErrMsg("CGOAssocs::AddGeneGOlist","Unable to allocate memory");
		return(eBSFerrMem);
		}

	if(m_pszGeneGoItemsTxt != NULL)
		{
		if(m_NumGeneGOitemsTxt)
			memmove(pChr,m_pszGeneGoItemsTxt,m_NumGeneGOitemsTxt);
		delete m_pszGeneGoItemsTxt;
		}
	else
		{
		m_AllocdGeneGoItemsTxt = 0;
		m_NumGeneGOitemsTxt = 0;
		}
	m_pszGeneGoItemsTxt = pChr;
	m_AllocdGeneGoItemsTxt += 100000000;
	}

pTmp = &m_pGeneGOitems[m_NumGeneGOitems++];
pChr = &m_pszGeneGoItemsTxt[m_NumGeneGOitemsTxt];
pTmp->ObjType = ObjType;
pTmp->pszGeneName = pChr;
strcpy(pChr,pszGene);
pChr += NameLen;

pTmp->pszQuals = pChr;
strcpy(pChr,pszQuals);
pChr += QualLen;

pTmp->pszAspects = pChr;
strcpy(pChr,pszAspects);
pChr += AspectsLen;

pTmp->pszEvidences = pChr;
strcpy(pChr,pszEvidences);
pChr += EvidencesLen;

pTmp->pszGOIDs = pChr;
strcpy(pChr,pszGOIDs);
pChr += GOIDsLen;

m_NumGeneGOitemsTxt += ReqLen;
return(eBSFSuccess);
}

// process gene GO list, dedupe and combine
int
CGOAssocs::ProcessGeneGOlist(void)
{
int Cnt;
int NumUnique;
tsGeneGOitem *pTmp;
tsGeneGOitem *pTmp1;
char *pChr;
char szQuals[cMaxGOIDsPerGene * 10];
char szAspects[cMaxGOIDsPerGene * 10];
char szEvidences[cMaxGOIDsPerGene * 10];
char szGOIDs[cMaxGOIDsPerGene * 20];

if(m_pGeneGOitems == NULL || !m_NumGeneGOitems)
	return(eBSFSuccess);

NumUnique = 1;
if(m_NumGeneGOitems > 1)
	{
	qsort(m_pGeneGOitems,m_NumGeneGOitems,sizeof(tsGeneGOitem),SortGeneGoItems);
	pTmp1 = m_pGeneGOitems;
	pTmp = &m_pGeneGOitems[1];
	strcpy(szQuals,pTmp1->pszQuals);
	strcpy(szAspects,pTmp1->pszAspects);
	strcpy(szEvidences,pTmp1->pszEvidences);
	strcpy(szGOIDs,pTmp1->pszGOIDs);
	for(Cnt = 1; Cnt < m_NumGeneGOitems; Cnt++, pTmp++)
		{
		if(!stricmp(pTmp->pszGeneName,pTmp1->pszGeneName) &&
			pTmp->ObjType == pTmp1->ObjType)
				{
				pChr = &szQuals[strlen(szQuals)];
				*pChr++ = '|';
				strcpy(pChr,pTmp->pszQuals);

				pChr = &szAspects[strlen(szAspects)];
				*pChr++ = '|';
				strcpy(pChr,pTmp->pszAspects);

				pChr = &szEvidences[strlen(szEvidences)];
				*pChr++ = '|';
				strcpy(pChr,pTmp->pszEvidences);

				pChr = &szGOIDs[strlen(szGOIDs)];
				*pChr++ = '|';
				strcpy(pChr,pTmp->pszGOIDs);
				continue;
				}

		// different gene name and object type
		DeDupeGOIDs(szQuals,szAspects,szEvidences,szGOIDs);	// remove duplcate GOIDs

		Add(pTmp1->ObjType,pTmp1->pszGeneName,szQuals,szAspects,szEvidences,szGOIDs);
		pTmp1 = pTmp;
		strcpy(szQuals,pTmp1->pszQuals);
		strcpy(szAspects,pTmp1->pszAspects);
		strcpy(szEvidences,pTmp1->pszEvidences);
		strcpy(szGOIDs,pTmp1->pszGOIDs);
		NumUnique++;
		}
	Add(pTmp1->ObjType,pTmp1->pszGeneName,szQuals,szAspects,szEvidences,szGOIDs);
	}
else
	Add(m_pGeneGOitems->ObjType,m_pGeneGOitems->pszGeneName,m_pGeneGOitems->pszQuals,m_pGeneGOitems->pszAspects,m_pGeneGOitems->pszEvidences,m_pGeneGOitems->pszGOIDs);
return(NumUnique);
}

typedef struct TAG_sDedupedGOID {
	char szGOID[16];
	char szGOQual[16];
	char szGOAspect[16];
	char szGOEvidence[16];
} tsDedupedGOID;

// remove duplcate GOIDs
int
CGOAssocs::DeDupeGOIDs(char *pszQuals,char *pszAspects,char *pszEvidences,char *pszGOIDs)
{
tsDedupedGOID Deduped[cMaxGOIDsPerGene];
tsDedupedGOID *pDeduped;
int nthField;
int ParseState;
int NumGOIDs;

tsDedupedGOID *pTmp;
tsDedupedGOID *pTmp1;
int Cnt;

nthField = 1;
NumGOIDs = 0;
pDeduped = Deduped;
do {
	ParseState = ParseNthField(pszGOIDs,pDeduped->szGOID,nthField);
	if(ParseState == 2 && stricmp(pDeduped->szGOID,"n/a"))
		{
		ParseNthField(pszQuals,pDeduped->szGOQual,nthField);
		ParseNthField(pszAspects,pDeduped->szGOAspect,nthField);
		ParseNthField(pszEvidences,pDeduped->szGOEvidence,nthField);
		pDeduped += 1;
		NumGOIDs += 1;
		}
	}
while(nthField++ && ParseState > 0);

// now can sort and dedupe
if(NumGOIDs > 1)
	{
	qsort(Deduped,NumGOIDs,sizeof(tsDedupedGOID),SortDedupeGOIDs);
	int NumUnique = 1;
	pTmp1 = Deduped;
	pTmp = &Deduped[1];
	for(Cnt = 1; Cnt < NumGOIDs; Cnt++, pTmp++)
		{
		if(!stricmp(pTmp->szGOID,pTmp1->szGOID))
			continue;
		if(NumUnique != Cnt)
			pTmp1[1] = *pTmp;
		pTmp1 += 1;
		NumUnique++;
		}
	NumGOIDs = NumUnique;
	}

// all deduped, back to list formats
*pszQuals = '\0';
*pszAspects = '\0';
*pszEvidences = '\0';
*pszGOIDs = '\0';
pTmp = Deduped;
for(Cnt = 0; Cnt < NumGOIDs; Cnt++,pTmp++)
	{
	if(Cnt >= 1)
		{
		strcat(pszQuals,"|");
		strcat(pszAspects,"|");
		strcat(pszEvidences,"|");
		strcat(pszGOIDs,"|");
		}

	strcat(pszQuals,pTmp->szGOQual);
	strcat(pszAspects,pTmp->szGOAspect);
	strcat(pszEvidences,pTmp->szGOEvidence);
	strcat(pszGOIDs,pTmp->szGOID);
	}

return(NumGOIDs);
}

int
CGOAssocs::ParseGeneMapping(char *pszGeneMapping)	// gene mapping file to parse
{
FILE *pGeneStream;
int Len;
int LineNum;
int Idx;
char szLineBuff[cLineBuffLen];
char *pTxt;
char *pChr;
char Chr;
char szMapTo[cMaxGeneNameLen];
char szMapFrom[cMaxGeneNameLen * cMaxMapFroms];	// allow cMaxMapFroms maximally sized map from genes
char szMapFroms[cMaxMapFroms][cMaxGeneNameLen];
int NumFroms;
int Cnt;

tsGeneMap *pTmp;
tsGeneMap *pTmp1;
int Rslt = eBSFSuccess;

if(m_bGOFlybase && m_NumGeneFilters && m_pGeneFilters != NULL)
	{
	tsGeneFilter *pFBGeneTrns;
	pTmp = new tsGeneMap [m_NumGeneFilters];
	if(pTmp == NULL)
		{
		AddErrMsg("CGOAssocs::ParseGeneMapping","Unable to allocate memory");
		return(eBSFerrMem);
		}
	m_AllocdGeneMaps = m_NumGeneFilters;
	m_pGeneMappings = pTmp;
	pFBGeneTrns = m_pGeneFilters;
	for(m_NumGeneMaps = 0; m_NumGeneMaps < m_NumGeneFilters; pTmp++,pFBGeneTrns++,m_NumGeneMaps++)
		{
		strcpy(pTmp->szToGene,pFBGeneTrns->szGene);
		pChr = pTmp->szFromGene;
		pTxt = pTmp->szToGene;
		while(Chr=(*pChr++ = *pTxt++))
			if(Chr == '-')	// used to separate cannonical gene name from transcript id e.g CG10000-RC
				{
				pChr[-1] = '\0';
				break;
				}
		}
	}

if(pszGeneMapping != NULL && *pszGeneMapping != '\0')
	{
	if((pGeneStream = fopen(pszGeneMapping,"r"))==NULL)
		{
		AddErrMsg("CGOAssocs::ParseGeneMapping","Unable to open gene mapping file %s error: %s",pszGeneMapping,strerror(errno));
		return(eBSFerrOpnFile);
		}

	LineNum = 0;
	szMapFrom[0] = '\0';
	szMapTo[0] = '\0';

	while(fgets(szLineBuff,sizeof(szLineBuff),pGeneStream)!= NULL)
		{
		LineNum++;
		// simply slough lines which were just whitespace or start with '#' as comment line
		Len = TrimWhitespace(szLineBuff);
		if(!Len || szLineBuff[0] == '#')	
			continue;
		
		pTxt = szLineBuff;
		
		// initial field expected to be the 'to' name
		if((Idx = ParseNxtField(pTxt,true,sizeof(szMapTo)-1,szMapTo,true)) < 0)
			{
			AddErrMsg("CGOAssocs::ParseGeneMapping","gene mapping file '%s' line %d parse error",pszGeneMapping,LineNum);
			fclose(pGeneStream);
			return(eBSFerrParse);
			}

		// next and subsequent expected to be the 'from' names
		// note that there may not be any 'from' names in which case the 'to' is simply sloughed
		pTxt += Idx;
		if((Idx = ParseNxtField(pTxt,false,sizeof(szMapFrom)-1,szMapFrom,false)) < 1)
			continue;
		
		if(!stricmp(szMapFrom,"n/a"))
			continue;

		pTxt = szMapFrom;
		pChr = szMapFroms[0];
		*pChr = '\0';
		NumFroms = 0;
		Len = 0;
		while(Chr = *pTxt++)
			{
			switch(Chr) {
				case ',':			// must be multiple mapped from's
					if(!Len)
						continue;
					*pChr = '\0';	// terminate current
					pChr = szMapFroms[++NumFroms];
					Len = 0;
					continue;

				case ' ':			// just slough any whitespace
				case '\t':
					continue;

				default:			// accumulate
					*pChr++ = Chr;
					Len += 1;
					continue;
				}
			}
		if(Len)
			{
			*pChr = '\0';
			NumFroms += 1;
			}
		if(m_pGeneMappings == NULL || (m_NumGeneMaps + NumFroms) >= m_AllocdGeneMaps)
			{
			pTmp = new tsGeneMap [m_AllocdGeneMaps + 10000];
			if(pTmp == NULL)
				{
				AddErrMsg("CGOAssocs::ParseGeneMapping","Unable to allocate memory");
				fclose(pGeneStream);
				return(eBSFerrMem);
				}

			if(m_pGeneMappings != NULL)
				{
				if(m_NumGeneMaps)
					memmove(pTmp,m_pGeneMappings,sizeof(tsGeneMap) * m_NumGeneMaps);
				delete m_pGeneMappings;
				}
			else
				{
				m_NumGeneMaps = 0;
				m_AllocdGeneMaps = 0;
				}
			m_AllocdGeneMaps += 10000;
			m_pGeneMappings = pTmp;
			}
		while(NumFroms--)
			{
			pTmp = &m_pGeneMappings[m_NumGeneMaps++];
			strcpy(pTmp->szFromGene,szMapFroms[NumFroms]);
			strcpy(pTmp->szToGene,szMapTo);
			}
		}
	fclose(pGeneStream);
	}

if(m_NumGeneMaps > 1) // sort and dedupe
	{
	qsort(m_pGeneMappings,m_NumGeneMaps,sizeof(tsGeneMap),SortMapGenes);
	int NumUnique = 1;
	pTmp1 = m_pGeneMappings;
	pTmp = &m_pGeneMappings[1];
	for(Cnt = 1; Cnt < m_NumGeneMaps; Cnt++, pTmp++)
		{
		if(!stricmp(pTmp->szFromGene,pTmp1->szFromGene) &&
			!stricmp(pTmp->szToGene,pTmp1->szToGene))
			continue;
		if(NumUnique != Cnt)
			pTmp1[1] = *pTmp;
		pTmp1 += 1;
		NumUnique++;
		}
	m_NumGeneMaps = NumUnique;
	}
m_pGeneMappings[0].nThInstance = 1;
if(m_NumGeneMaps > 1)
	{
	pTmp1 = m_pGeneMappings;
	pTmp = &m_pGeneMappings[1];
	for(Cnt = 1; Cnt < m_NumGeneMaps; Cnt++, pTmp++,pTmp1++)
		{
		if(!stricmp(pTmp->szFromGene,pTmp1->szFromGene))
			pTmp->nThInstance = pTmp1->nThInstance + 1;
		else
			pTmp->nThInstance = 1;
		}
	}
return(m_NumGeneMaps);
}

int
CGOAssocs::ParseGeneFilters(char *pszGeneFilters)	// gene filter file to parse
{
FILE *pFilterStream;
int Len;
int LineNum;
int Idx;
char szLineBuff[cLineBuffLen];
char *pTxt;
char szChrom[cMaxDatasetSpeciesChrom];
int chromStart;
int chromEnd;
char szGene[cMaxGeneNameLen];
bool bIsBED;

int Cnt;

tsGeneFilter *pTmp;
tsGeneFilter *pTmp1;
int Rslt = eBSFSuccess;

if((pFilterStream = fopen(pszGeneFilters,"r"))==NULL)
	{
	AddErrMsg("CGOAssocs::ParseGeneFilters","Unable to open gene filter file %s error: %s",pszGeneFilters,strerror(errno));
	return(eBSFerrOpnFile);
	}

LineNum = 0;
szGene[0] = '\0';

bIsBED = true;	// assume BED file, otherwise use 1st field as the gene name

while(fgets(szLineBuff,sizeof(szLineBuff),pFilterStream)!= NULL)
	{
	LineNum++;
	// simply slough lines which were just whitespace or start with '#' as comment line
	Len = TrimWhitespace(szLineBuff);
	if(!Len || szLineBuff[0] == '#')	
		continue;
	
	pTxt = szLineBuff;
	
	if(bIsBED)
		{
		Cnt = sscanf(pTxt," %36s %d %d %50s",
				szChrom,&chromStart,&chromEnd,szGene);
		if(Cnt != 4)
			{
			if(!m_NumGeneFilters)	// if parse fails for BED on 1st gene then assume 1st field is gene
				bIsBED = false;
			else
				{
				AddErrMsg("CGOAssocs::ParseGeneFilters","gene filter file (BED) '%s' line %d parse error",pszGeneFilters,LineNum);
				fclose(pFilterStream);
				return(eBSFerrParse);
				}
			}
		}
	
	if(!bIsBED)
		{
		// initial field expected to be the 'to' name, subsequent fields are simply sloughed
		if((Idx = ParseNxtField(pTxt,true,sizeof(szGene)-1,szGene,true)) < 0)
			{
			AddErrMsg("CGOAssocs::ParseGeneFilters","gene filter file (non-BED) '%s' line %d parse error",pszGeneFilters,LineNum);
			fclose(pFilterStream);
			return(eBSFerrParse);
			}
		}

	if(!stricmp(szGene,"n/a"))
		continue;

	if(m_pGeneFilters == NULL || (m_NumGeneFilters + 1) >= m_AllocdGeneFilters)
		{
		pTmp = new tsGeneFilter [m_AllocdGeneFilters + 10000];
		if(pTmp == NULL)
			{
			AddErrMsg("CGOAssocs::ParseGeneFilters","Unable to allocate memory");
			fclose(pFilterStream);
			return(eBSFerrMem);
			}

		if(m_pGeneFilters != NULL)
			{
			if(m_NumGeneFilters)
				memmove(pTmp,m_pGeneFilters,sizeof(tsGeneFilter) * m_NumGeneFilters);
			delete m_pGeneFilters;
			}
		else
			{
			m_NumGeneFilters = 0;
			m_AllocdGeneFilters = 0;
			}
		m_AllocdGeneFilters += 10000;
		m_pGeneFilters = pTmp;
		}
	pTmp = &m_pGeneFilters[m_NumGeneFilters++];
	strcpy(pTmp->szGene,szGene);
	}
fclose(pFilterStream);
if(m_NumGeneFilters > 1) // sort and dedupe
	{
	qsort(m_pGeneFilters,m_NumGeneFilters,sizeof(tsGeneFilter),SortGeneFilters);
	int NumUnique = 1;
	pTmp1 = m_pGeneFilters;
	pTmp = &m_pGeneFilters[1];
	for(Cnt = 1; Cnt < m_NumGeneFilters; Cnt++, pTmp++)
		{
		if(!stricmp(pTmp->szGene,pTmp1->szGene))
			continue;
		if(NumUnique != Cnt)
			pTmp1[1] = *pTmp;
		pTmp1 += 1;
		NumUnique++;
		}
	m_NumGeneFilters = NumUnique;
	}
return(m_NumGeneFilters);
}

bool								// true if mapped
CGOAssocs::MapGene(char *pszGene,int nThInstance)	// gene to be mapped
{
tsGeneMap *pMatch;
if((pMatch = LocateGeneMap(pszGene,nThInstance))==NULL)
	return(false);
strcpy(pszGene,pMatch->szToGene);
return(true);
}

// ParseNxtField
// Parses next field value from pszLine
// Field values are expected to be delimited by tab or (if bAcceptComma true) comma
int											// total number of chrs parsed from pszLine or -1 if errors
CGOAssocs::ParseNxtField(char *pszLine,
						 bool bReqValue,	// value is mandatory
						 int MaxLen,		// max allowed ret value
						 char *pszRetValue,	// where to return parsed value
						 bool bAcceptComma) // accept comma delimited
{
int ValLen = 0;
int ChrCnt = 0;
while(*pszLine != '\0' && *pszLine != '\n' && *pszLine != '\r')
	{
	ChrCnt++;
	if(*pszLine == '\t' || (bAcceptComma && *pszLine == ','))
		break;
	// note following strips leading spaces...
	if((ValLen > 0 || *pszLine != ' ') && ValLen < MaxLen)
		{
		*pszRetValue++ = *pszLine;
		ValLen++;
		}
	pszLine++;
	}
*pszRetValue = '\0';
if(bReqValue && !ValLen)
	return(-1);
return(ChrCnt);
}
 

// ParseUCSC
// Expected file structure:
// <genename><tab><SwisProtIDs><tab><GOquals><tab><GOIDs><tab><GOAspects>
// Note that there may be one or more SwissProtIDs,GOquals,GOIDs and GOAspects which will be comma delimited
// Note that if there are multiple comma delimited terms then the last term has a trailing comma,
// e.g. "NM_145243\tQ5T3G7,Q5T3G6,\t,,,\tGO:0004222,GO:0006508,GO:0016020,\tF,P,C,\n"
// Note in example above how there are no GOquals specified but commas are used because there are multiple GOIDs
// e.g. "NM_007324\tQ5T0F6\t\tGO:0008270\tF\n"
// Note in example above that the missing GOqual is marked by 2 tabs because there is a single GOIDs
// 
// If gene has no SwisProtIDs then 'n/a' is expected
// If gene has no associated GO terms then 'n/a\tn/a\tn/a' is expected
// 
int
CGOAssocs::ParseUCSC(char *pszUCSCGoAssoc,char *pszGeneFilters)	// UCSC file to parse
{
FILE *pAssocStream;
int Len;
int Cnt;
int LineNum;
int Idx;
char szLineBuff[cLineBuffLen];
char *pTxt;
char szGeneID[50];
char szSwissProtIDs[2000];
char szQuals[2000];
char szAspects[2000];
char szgoIDs[2000];
int Rslt = eBSFSuccess;

if(pszGeneFilters != NULL && pszGeneFilters[0] != '\0')
	{
	if((Idx = ParseGeneFilters(pszGeneFilters)) < 0)
		return(Idx);
	if(!Idx)
		{
		AddErrMsg("CGOAssocs::ParseUCSC","Gene filters file parsed but no genes to filter! %s",pszGeneFilters);
		return(eBSFerrOpnFile);
		}
	m_bDoGeneFilters = true;
	}
else
	m_bDoGeneFilters = false;

if((pAssocStream = fopen(pszUCSCGoAssoc,"r"))==NULL)
	{
	AddErrMsg("CGOAssocs::Open","Unable to open GO association file %s error: %s",pszUCSCGoAssoc,strerror(errno));
	return(eBSFerrOpnFile);
	}

LineNum = 0;
while(fgets(szLineBuff,sizeof(szLineBuff),pAssocStream)!= NULL)
	{
	LineNum++;
	// simply slough lines which were just whitespace or which start with '#' pr ';'
	Len = TrimWhitespace(szLineBuff);
	if(!Len || szLineBuff[0] == '#'  || szLineBuff[0] == ';')	
		continue;
	
	pTxt = szLineBuff;
	// parse out gene and SwissProtIDs
	Cnt = sscanf(pTxt,"%50s\t%2000s%n",
		szGeneID,szSwissProtIDs,&Idx);
	if(Cnt != 2)
		{
		AddErrMsg("CGOAssocs::Open","GO association file line %d parse error",LineNum);
		fclose(pAssocStream);
		return(eBSFerrParse);
		}
	
	// check if missing GOQuals (double tab)
	if(pTxt[Idx] == '\t' && pTxt[Idx+1] == '\t')
		{
		szQuals[0] = '\0';
		Cnt += sscanf(&pTxt[Idx],"\t%2000s\t%2000s",
			szgoIDs,szAspects);
		}
	else
		{
		Cnt += sscanf(&pTxt[Idx],"\t%2000s\t%2000s\t%2000s",
			szQuals,szgoIDs,szAspects);
		}
	if(Cnt < 4)
		{
		AddErrMsg("CGOAssocs::Open","GO association file line %d parse error",LineNum);
		fclose(pAssocStream);
		return(eBSFerrParse);
		}

	if((Rslt = AddGeneGOlist(eGOTgene,szGeneID,szQuals,szAspects,"IEA",szgoIDs)) < 0)
		{
		fclose(pAssocStream);
		return(Rslt);
		}
	}
fclose(pAssocStream);
ProcessGeneGOlist();
return(Rslt);
}


int 
CGOAssocs::Add(etGOobjType ObjType,
   			   const char *pszGene,			// gene name
			   const char *pszGOQuals,		// list of GO qualifiers
			   const char *pszGOAspects,		// list of GO aspects
			   const char *pszGOEvidences,	// list of GO evidence
			   const char *pszGOIDs)			// list of GO:id
{
int Genelen;
int GOIDslen;

int nthField;
char szGOQual[100];		// currently parsed out from list of GO qualifiers
char szGOAspect[100];	// currently parsed out from list of GO aspects
char szGOEvidence[100];	// currently parsed out from list of GO evidence
char szGOID[100];		// currently parsed out from list of GO:id
int ParseState;

int Mem2Alloc;
int Assoc2Alloc;
int Gene2Alloc;
char *pTmp;
tsGOAssoc *pAssoc;
tsGOGene *pGene;

SwitchGeneNamePTR2IDX();	// if gene names are represented as being ptrs then need to switch to relative offsets
							// this is required in case reallocations are required

Genelen = (int)strlen(pszGene);

if(pszGOIDs != NULL && *pszGOIDs != '\0')
	GOIDslen = (int)strlen(pszGOIDs);
else
	GOIDslen = 0;

// time to allocate or realloc?
if(m_pszGOIdents == NULL || (m_NxtGOIdentOfs + GOIDslen + 1) >= m_MaxAllocGOIdents)
	{
	if(m_pszGOIdents == NULL)
		{
		m_NxtGOIdentOfs = 0;
		Mem2Alloc = 2 * cAllocIdents;
		}
	else
		Mem2Alloc = m_MaxAllocGOIdents + cAllocIdents;
	pTmp = new char[Mem2Alloc];
	if(pTmp == NULL)
		return(eBSFerrMem);
	if(m_pszGOIdents != NULL && m_NxtGOIdentOfs)
		memmove(pTmp,m_pszGOIdents,m_NxtGOIdentOfs);
	if(m_pszGOIdents != NULL)
		delete m_pszGOIdents;
	m_pszGOIdents = pTmp;
	m_MaxAllocGOIdents = Mem2Alloc;
	}

// time to allocate or realloc?
if(m_pszGeneIdents == NULL || (m_NxtGeneIdentOfs + Genelen + 1) >= m_MaxAllocGeneIdents)
	{
	if(m_pszGeneIdents == NULL)
		{
		m_NxtGeneIdentOfs = 0;
		Mem2Alloc = 2 * cAllocIdents;
		}
	else
		Mem2Alloc = m_MaxAllocGeneIdents + cAllocIdents;
	pTmp = new char[Mem2Alloc];
	if(pTmp == NULL)
		return(eBSFerrMem);
	if(m_pszGeneIdents != NULL && m_NxtGeneIdentOfs)
		memmove(pTmp,m_pszGeneIdents,m_NxtGeneIdentOfs);
	if(m_pszGeneIdents != NULL)
		delete m_pszGeneIdents;
	m_pszGeneIdents = pTmp;
	m_MaxAllocGeneIdents = Mem2Alloc;
	}


// time to allocate?
if(m_pAssocs == NULL || (m_GOAssocCnt + cMaxGOIDsPerGene) >= m_MaxAllocGOAssoc)
	{
	if(m_pAssocs == NULL || !m_GOAssocCnt)
		Assoc2Alloc = 2 * cAllocNumAssoc;
	else
		Assoc2Alloc = (m_MaxAllocGOAssoc + cAllocNumAssoc);
	Mem2Alloc = Assoc2Alloc * sizeof(tsGOAssoc);
	pTmp = new char[Mem2Alloc];
	if(pTmp == NULL)
		return(eBSFerrMem);
	if(m_pAssocs != NULL && m_GOAssocCnt)
		memmove(pTmp,(char *)m_pAssocs,m_GOAssocCnt * sizeof(tsGOAssoc));
	if(m_pAssocs != NULL)
		delete (char *)m_pAssocs;
	m_pAssocs = (tsGOAssoc *)pTmp;
	m_MaxAllocGOAssoc = Assoc2Alloc;
	}

// time to allocate?
if(m_pGenes == NULL || m_GOGenesCnt >= m_MaxAllocGOGene)
	{
	if(m_pGenes == NULL || !m_GOGenesCnt)
		Gene2Alloc = 2 * cAllocNumGenes;
	else
		Gene2Alloc = (m_MaxAllocGOGene + cAllocNumGenes);
	Mem2Alloc = Gene2Alloc * sizeof(tsGOGene);
	pTmp = new char[Mem2Alloc];
	if(pTmp == NULL)
		return(eBSFerrMem);
	if(m_pGenes != NULL && m_GOGenesCnt)
		memmove(pTmp,(char *)m_pGenes,m_GOGenesCnt * sizeof(tsGOGene));
	if(m_pGenes != NULL)
		delete (char *)m_pGenes;
	m_pGenes = (tsGOGene *)pTmp;
	m_MaxAllocGOGene = Gene2Alloc;
	}

pGene = &m_pGenes[m_GOGenesCnt++];
pGene->GeneID = m_GOGenesCnt;
pGene->Name.Idx = m_NxtGeneIdentOfs;
strcpy(&m_pszGeneIdents[m_NxtGeneIdentOfs],pszGene);
m_NxtGeneIdentOfs += ((int)strlen(pszGene)+1);
m_GeneIdentsCnt++;
pGene->GOAssocIdx = -1;
pGene->NumGOAssoc = 0;

if(GOIDslen)
	{
	nthField = 1;
	do {
		ParseState = ParseNthField(pszGOIDs,szGOID,nthField);
		if(ParseState == 2 && stricmp(szGOID,"n/a"))
			{
			ParseNthField(pszGOQuals,szGOQual,nthField);
			ParseNthField(pszGOAspects,szGOAspect,nthField);
			ParseNthField(pszGOEvidences,szGOEvidence,nthField);
			pAssoc = &m_pAssocs[m_GOAssocCnt];
			pAssoc->Qual = MapGOQual(szGOQual);
			pAssoc->Aspect = MapGOaspect(szGOAspect);
			if((pAssoc->Evidence = MapGOEvidence(szGOEvidence)) == eGOEnone)
				pAssoc->Evidence = m_DfltEvidence;
			pTmp = &m_pszGOIdents[m_NxtGOIdentOfs];
			strcpy(pTmp,szGOID);
			if(!pGene->NumGOAssoc)
				pGene->GOAssocIdx = m_GOAssocCnt;
			pAssoc->IdentIdx = m_NxtGOIdentOfs;
			m_NxtGOIdentOfs += ((int)strlen(szGOID) + 1);
			m_GOIdentsCnt++;
			pAssoc->GeneID = pGene->GeneID;
			pGene->NumGOAssoc++;
			m_GOAssocCnt++;
			}
		}
	while(nthField++ && ParseState > 0);
	}
return(eBSFSuccess);
}

void
CGOAssocs::SwitchGeneNameIDX2PTR(void)
{
tsGOGene *pGene;

/// if already ptrs, or no genes then easy return
if(m_bGeneNamesAsPtrs)
	return;
if((pGene = m_pGenes)==NULL || !m_GOGenesCnt)
	{
	m_bGeneNamesAsPtrs = true;
	return;
	}

for(int Cnt=0;Cnt < m_GOGenesCnt; Cnt++,pGene++)
	pGene->Name.ptr = &m_pszGeneIdents[pGene->Name.Idx];
m_bGeneNamesAsPtrs = true;
}

void
CGOAssocs::SwitchGeneNamePTR2IDX(void)
{
tsGOGene *pGene;

// if already idx, or no genes then easy return
if(!m_bGeneNamesAsPtrs)
	return;
if((pGene = m_pGenes)==NULL || !m_GOGenesCnt)
	{
	m_bGeneNamesAsPtrs = false;
	return;
	}

for(int Cnt=0;Cnt < m_GOGenesCnt; Cnt++,pGene++)
	pGene->Name.Idx = (int)((char *)pGene->Name.ptr - m_pszGeneIdents);
m_bGeneNamesAsPtrs = false;
}


// TrimWhitespace
// Inplace trim any leading/trailing whitespace
// Returns number of chrs in trimmed string excluding terminating '\0'
int
CGOAssocs::TrimWhitespace(char *pTxt)
{
int Len = 0;
char *pStart = pTxt;
char Chr;
if(pTxt == NULL || *pTxt == '\0')
	return(0);

	// strip any leading whitespace
while((Chr = *pTxt) && isspace(Chr))
	pTxt++;
if(Chr == '\0')					// empty line?
	{
	*pStart = '\0';
	return(0);
	}

while(*pTxt)			// fast forward to string terminator '\0' 
	{
	*pStart++ = *pTxt++;
	Len++;
	}
pStart--;				// backup to last chr 
while(isspace(*pStart))
	{
	pStart--;
	Len--;
	}
pStart[1] = '\0';
return(Len);
}

// StripQuotesWS
// Inplace strips any bracketing double quotes plus leading/trailing whitespace
int							// returns length of quote and whitespace stripped *pszTxt
CGOAssocs::StripQuotesWS(char *pszTxt)
{
char QuoteChr = '\0';
int Len = 0;
char *pDst = pszTxt;

if(pszTxt == NULL || *pszTxt == '\0')
	return(0);

// first slough any leading whitespace
while(*pszTxt && isspace(*pszTxt))
	pszTxt++;
// check if left with empty string
if(*pszTxt == '\0')
	{
	*pDst = '\0';
	return(0);
	}

// check if single or double quote
if(*pszTxt == '"' || *pszTxt == '\'')
	{
	QuoteChr = *pszTxt;		// note which quote was used
	pszTxt++;
	// could be more whitespace, slough 
	while(*pszTxt && isspace(*pszTxt))
		pszTxt++;
	// check if left with empty string
	if(*pszTxt == '\0')
		{
		*pDst = '\0';
		return(0);
		}
	}

// any leading whitespace plus quote char have been trimed
while(*pszTxt)
	{
	Len++;
	*pDst++ = *pszTxt++;
	}

pDst--;	// backup to pt at last chr, not at terminating '\0'
while(isspace(*pDst))	// slough any trailing whitespace
	{
	pDst--;
	Len--;
	}

// check for quote chr of same type which started string
if(Len && *pDst == QuoteChr)	
	{
	pDst--;				// slough trailing quote chr
	Len--;
	// could be more whitespace, slough 
	while(Len && isspace(*pDst))
		{
		Len--;
		pDst--;
		}
	}
if(!Len)
	*pDst = '\0';
else
	pDst[1] = '\0';
return(Len);
}


// ParseNthField
// Parses out and returns the nth field from pszToParse
// If pszToParse contains multiple fields then fields are expected to be either comma or '|' separated
// with last field optionally containing trailing comma e.g. "value,value," 
// Returns:
// -1   parameter error
// 0	no more fields to parse - nthField > number of fields in pszToParse
// 1	if field present but empty e.g. "value,,value"
// 2    "fieldcontents"
//
int
CGOAssocs::ParseNthField(const char *pszToParse,char *pszRetField,int nthField)
{
char *pszField;
char Chr;
if(pszRetField == NULL || pszToParse == NULL || nthField < 1)
	return(-1);
pszRetField[0] = '\0';
if(*pszToParse == '\0')
	return(0);
pszField = pszRetField;
while((Chr = *pszToParse++) != '\0') 
	{
	if(nthField == 1)
		{
		if(Chr == ',' || Chr == '|')
			{
			*pszField = '\0';
			break;
			}
		*pszField++ = Chr;
		continue;
		}

	if(Chr == ',' || Chr == '|')
		nthField--;
	}

if(Chr == '\0' && pszRetField[0] == '\0')
	return(0);
*pszField = '\0';
return(pszRetField[0] == '\0' ? 1 : 2);	
}


etGOaspect
CGOAssocs::MapGOaspect(char *pszTxt)
{
int Len = TrimWhitespace(pszTxt);
if(!Len || !stricmp(pszTxt,"n/a"))
	return(eGOaspectUndefined);
switch(tolower(*pszTxt)) {
	case 'p': return(eGOaspectP);
	case 'f': return(eGOaspectF);
	case 'c': return(eGOaspectC);
	default:
		break;
	}
return(eGOaspectUndefined);
}

const char *
CGOAssocs::MapGOaspect2Txt(etGOaspect GOaspect)
{
switch(GOaspect) {
	case eGOaspectP: return("P biological process");
	case eGOaspectF: return("F molecular function");
	case eGOaspectC: return("C cellular process");
	default: break;
	}
return("none");
}


etGOQual
CGOAssocs::MapGOQual(char *pszTxt)
{
int Qual;
char *pszQual;
char szQual[1024];
int QualLen;

if(pszTxt == NULL || *pszTxt == '\0')
	return(eGOQnone);

strncpy(szQual,pszTxt,sizeof(szQual));
szQual[sizeof(szQual)-1] = '\0';

QualLen = TrimWhitespace(szQual);
if(!QualLen || !stricmp(szQual,"n/a"))
	return(eGOQnone);

Qual = (int)eGOQnone;
pszQual = szQual;
while(*pszQual != '\0')
	{
	if(*pszQual == '\t' || *pszQual == ' ' || *pszQual == '|' || *pszQual == ',')
		{
		pszQual++;
		continue;
		}
	if(!strnicmp(pszQual,"not",3))
		{
		Qual |= (int)eGOQNOT;
		pszQual += 3;
		continue;
		}
	if(!strnicmp(pszQual,"contributes_to",14))
		{
		Qual |= (int)eGOQcontributes_to;
		pszQual += 14;
		continue;
		}
	if(!strnicmp(pszQual,"colocalizes_with",16))
		{
		Qual |= (int)eGOQcolocalizes_with;
		pszQual += 16;
		continue;
		}
	break;
	}

return((etGOQual)Qual);
}

const char *
CGOAssocs::MapGOQual2Txt(etGOQual GOqual)
{
static char szQual[512];
if(GOqual == eGOQnone)
	return("none");
szQual[0] = '\0';
if(GOqual & (int)eGOQNOT)
	strcpy(szQual,"NOT");
if(GOqual & (int)eGOQcontributes_to)
	{
	if(GOqual & (int)eGOQNOT)
		strcat(szQual,"|Contributes_to");
	else
		strcpy(szQual,"Contributes_to");
	}
if(GOqual & (int)eGOQcolocalizes_with)
	{
	if(GOqual & ((int)eGOQNOT | (int)eGOQcontributes_to))
		strcat(szQual,"|Colocalizes_with");
	else
		strcpy(szQual,"Colocalizes_with");
	}
return(szQual);
}

etGOobjType 
CGOAssocs::MapGOobjType(char *pszTxt)
{
int Len = TrimWhitespace(pszTxt);
if(!Len || !stricmp(pszTxt,"n/a"))
	return(eGOTundefined);

if(!stricmp(pszTxt,"Gene"))
	return(eGOTgene);
if(!stricmp(pszTxt,"Transcript"))
	return(eGOTtranscript);
if(!stricmp(pszTxt,"Protein"))
	return(eGOTprotein);
if(!stricmp(pszTxt,"Protein_structure"))
	return(eGOTprotein_structure);
if(!stricmp(pszTxt,"Complex"))
	return(eGOTcomplex);

return(eGOTundefined);
}

const char *
CGOAssocs::MapGOobjType2Txt(etGOobjType GOobjType)
{
switch(GOobjType) {
	case eGOTgene: return("Gene");
	case eGOTtranscript: return("Transcript");
	case eGOTprotein: return("Protein");
	case eGOTprotein_structure: return("Protein_structure");
	case eGOTcomplex: return("Complex");
	default: break;
	}
return("?none");
}

etGOEvidence 
CGOAssocs::MapGOEvidence(char *pszTxt)
{
int Evidence;
char *pszEvidence;
char szEvidence[512];
int Len;
if(pszTxt == NULL || *pszTxt == '\0')
	return(eGOEnone);
strncpy(szEvidence,pszTxt,sizeof(szEvidence));
szEvidence[sizeof(szEvidence)-1] = '\0';
Len = TrimWhitespace(szEvidence);
if(!Len || !stricmp(szEvidence,"n/a"))
	return(eGOEnone);

pszEvidence = szEvidence;
Evidence = (int)eGOEnone;
while(*pszEvidence != '\0')
	{
	if(*pszEvidence == '\t' || *pszEvidence == ' ' || *pszEvidence == '|' || *pszEvidence == ',')
		{
		pszEvidence++;
		continue;
		}
	if(!strnicmp(pszEvidence,"IMP",3))
		{
		Evidence |= (int)eGOEIMP;
		pszEvidence += 3;
		continue;
		}
	if(!strnicmp(pszEvidence,"IGI",3))
		{
		Evidence |= (int)eGOEIGI;
		pszEvidence += 3;
		continue;
		}
	if(!strnicmp(pszEvidence,"IPI",3))
		{
		Evidence |= (int)eGOEIPI;
		pszEvidence += 3;
		continue;
		}

	if(!strnicmp(pszEvidence,"ISS",3))
		{
		Evidence |= (int)eGOEISS;
		pszEvidence += 3;
		continue;
		}
	
	if(!strnicmp(pszEvidence,"IDA",3))
		{
		Evidence |= (int)eGOEIDA;
		pszEvidence += 3;
		continue;
		}

	if(!strnicmp(pszEvidence,"IEP",3))
		{
		Evidence |= (int)eGOEIEP;
		pszEvidence += 3;
		continue;
		}
	
	
	if(!strnicmp(pszEvidence,"IEA",3))
		{
		Evidence |= (int)eGOEIEA;
		pszEvidence += 3;
		continue;
		}
	
	if(!strnicmp(pszEvidence,"TAS",3))
		{
		Evidence |= (int)eGOETAS;
		pszEvidence += 3;
		continue;
		}
	
	if(!strnicmp(pszEvidence,"NAS",3))
		{
		Evidence |= (int)eGOENAS;
		pszEvidence += 3;
		continue;
		}
	
	
	if(!strnicmp(pszEvidence,"END",3))
		{
		Evidence |= (int)eGOEND;
		pszEvidence += 3;
		continue;
		}

	if(!strnicmp(pszEvidence,"EIC",3))
		{
		Evidence |= (int)eGOEIC;
		pszEvidence += 3;
		continue;
		}

	if(!strnicmp(pszEvidence,"ERC",3))
		{
		Evidence |= (int)eGOERC;
		pszEvidence += 3;
		continue;
		}	
	break;
	}
return((etGOEvidence)Evidence);
}

char *
CGOAssocs::MapGOEvidence2Txt(etGOEvidence GOEvidence)
{
int Idx;
static char szEvidence[512];
const char *pszCurEvidence;
int EvidenceMsk = 0x01;
bool bPrev = false;
if(GOEvidence == eGOEnone)
	return((char *)"none");
szEvidence[0] = '\0';

for(Idx=0;Idx < 12 && EvidenceMsk <= (int)GOEvidence;Idx++,EvidenceMsk <<= 1)
	{
	switch((etGOEvidence)((int)GOEvidence & EvidenceMsk)) {
		case eGOEIMP:				//inferred from mutant phenotype
			pszCurEvidence = "IMP";
			break;
		case eGOEIGI:				//inferred from genetic interaction [with <database:gene_symbol[allele_symbol]>]
			pszCurEvidence = "IGI";
			break;
		case eGOEIPI:				//inferred from physical interaction [with <database:protein_name>]
			pszCurEvidence = "IPI";
			break;
		case eGOEISS:				//inferred from sequence similarity [with <database:sequence_id>]
			pszCurEvidence = "ISS";
			break;
		case eGOEIDA:				//inferred from direct assay
			pszCurEvidence = "IDA";
			break;
		case eGOEIEP:				//inferred from expression pattern
			pszCurEvidence = "IEP";
			break;
		case eGOEIEA:				//inferred from electronic annotation [with <database:id>]
			pszCurEvidence = "IEA";
			break;
		case eGOETAS:				//traceable author statement
			pszCurEvidence = "TAS";
			break;
		case eGOENAS:				//non-traceable author statement
			pszCurEvidence = "NAS";
			break;
		case eGOEND:				//no biological data available
			pszCurEvidence = "END";
			break;
		case eGOEIC:				//inferred from reviewed computational analysis
			pszCurEvidence = "EIC";
			break;
		case eGOERC:				//inferred by curator [from <GO:id>]
			pszCurEvidence = "ERC";
			break;
		default:
			continue;
		}
	if(bPrev)
		strcat(szEvidence,"|");
	strcat(szEvidence,pszCurEvidence);
	bPrev = true;
	}
return((char *)szEvidence);
}

// LocateGene performs a binary search for specified gene
// Returns a ptr to term identified by specified gene or NULL if unable to locate term
// Uses binary search
tsGOGene *
CGOAssocs::LocateGene(char *pszGene)
{
tsGOGene *pProbe;
char *pszID;
int Rslt;
if(pszGene == NULL || *pszGene == '\0')
	return(NULL);
int Lo,Mid,Hi;	// search limits
Lo = 0; Hi = m_GOGenesCnt-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pProbe = &m_pGenes[Mid];
	pszID = ((char *)m_pszGeneIdents + pProbe->Name.Idx);
	Rslt = stricmp(pszGene,pszID);
	if(Rslt < 0)	
		{
		Hi = Mid - 1;
		continue;
		}
	if(Rslt > 0)	
		{
		Lo = Mid + 1;
		continue;
		}
	return(pProbe);
	}
return(NULL);
}

// LocateGeneFilter performs a binary search for specified gene
// Returns ptr to located gene filter if gene filter present or NULL if unable to locate filter gene
// Uses binary search
tsGeneFilter *
CGOAssocs::LocateGeneFilter(const char *pszGeneFilter)
{
tsGeneFilter *pProbe;
int Rslt;
if(pszGeneFilter == NULL || *pszGeneFilter == '\0' || m_pGeneFilters == NULL || m_NumGeneFilters <= 0)
	return(NULL);
int Lo,Mid,Hi;	// search limits
Lo = 0; Hi = m_NumGeneFilters-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pProbe = &m_pGeneFilters[Mid];
	Rslt = stricmp(pszGeneFilter,pProbe->szGene);
	if(Rslt < 0)	
		{
		Hi = Mid - 1;
		continue;
		}
	if(Rslt > 0)	
		{
		Lo = Mid + 1;
		continue;
		}
	return(pProbe);
	}
return(NULL);
}


tsGeneMap *
CGOAssocs::LocateGeneMap(char *pszGene, int nThInstance)
{
tsGeneMap *pProbe;
char *pszFromGene;
int Rslt;
if(pszGene == NULL || *pszGene == '\0' || nThInstance == 0)
	return(NULL);

int Lo,Mid,Hi;	// search limits
Lo = 0; Hi = m_NumGeneMaps-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pProbe = &m_pGeneMappings[Mid];
	pszFromGene = pProbe->szFromGene;
	Rslt = stricmp(pszGene,pszFromGene);
	if(Rslt < 0)	
		{
		Hi = Mid - 1;
		continue;
		}
	if(Rslt > 0)	
		{
		Lo = Mid + 1;
		continue;
		}
	if(nThInstance == pProbe->nThInstance)
		return(pProbe);

	if(nThInstance > pProbe->nThInstance)
		Lo = Mid + 1;
	else
		Hi = Mid - 1;
	}
return(NULL);
}


int
CGOAssocs::GetNumGOGenes(etOntologies OntologyClass)
{
if(m_pGenes == NULL || !m_GOGenesCnt || m_pAssocs == NULL || !m_GOAssocCnt)
	return(0);

switch(OntologyClass) {
	case eONTBiological:	// biological process P
		return(m_GOBIOGenesCnt);
	case eONTMolecular:		// molecular function F
		return(m_GOMOLGenesCnt);
	case eONTCellular:		// cellular process C
		return(m_GOCELGenesCnt);
	default:
		break;
	}
return(0);
}

int
CGOAssocs::GenGONumGenes(void)
{
int Idx;
int AssocIdx;
tsGOAssoc *pAssoc;
tsGOGene *pGene;
int BIOcnt;
int MOLcnt;
int CELcnt;

if(m_pGenes == NULL || !m_GOGenesCnt || m_pAssocs == NULL || !m_GOAssocCnt)
	return(-1);

m_GOBIOGenesCnt = m_GOMOLGenesCnt = m_GOCELGenesCnt = 0;
pGene = m_pGenes;
for(Idx = 0; Idx < m_GOGenesCnt; Idx++,pGene++)
	{
	if(pGene->NumGOAssoc)
		{
		BIOcnt = MOLcnt = CELcnt = 0;
		pAssoc = &m_pAssocs[pGene->GOAssocIdx];
		for(AssocIdx = 0; AssocIdx < pGene->NumGOAssoc; AssocIdx++,pAssoc++)
			{
			switch(pAssoc->Aspect & 0x07) {
				case eGOaspectP:	// biological process P
					BIOcnt++;
					break;
				case eGOaspectF:	// molecular function F
					MOLcnt++;
					break;
				case eGOaspectC:	// cellular process C
					CELcnt++;
					break;
				default:
					break;
				}
			}
		if(BIOcnt)
			m_GOBIOGenesCnt += 1;
		if(MOLcnt)
			m_GOMOLGenesCnt += 1;
		if(CELcnt)
			m_GOCELGenesCnt += 1;
		}
	}
return(0);
}

int 
CGOAssocs::GetNumGenes(void)					// returns total number of genes
{
return(m_GOGenesCnt);
}

int 
CGOAssocs::GetNumGeneGOitems(void)					// returns total number of gene GO items
{
return(m_NumGeneGOitems);
}

int 
CGOAssocs::GetNumGeneMaps(void)					// returns total number of gene mappings
{
return(m_NumGeneMaps);
}

char *
CGOAssocs::GetGene(int NthGene)				// returns Nth gene name
{
tsGOGene *pGene;
if(NthGene < 1 || NthGene > m_GOGenesCnt)
	return(NULL);
pGene = &m_pGenes[NthGene-1];

return((char *)m_pszGeneIdents + pGene->Name.Idx);
}

int 
CGOAssocs::GetNumGOIDs(char *pszGeneName)			// returns number of GO Identifiers associated with gene
{
tsGOGene *pGene;

// firstly locate the gene
if((pGene = LocateGene(pszGeneName)) == NULL)
	return(eBSFerrGOID);
return(pGene->NumGOAssoc);
}

char *
CGOAssocs::GetGOID(char *pszGeneName,int NthTerm)	// returns Nth GO Identifier associated with gene or NULL if none
{
tsGOGene *pGene;
tsGOAssoc *pAssoc;

// firstly locate the gene
if((pGene = LocateGene(pszGeneName)) == NULL)
	return(NULL);
if(NthTerm < 1 || NthTerm > pGene->NumGOAssoc)
	return(NULL);

pAssoc = &m_pAssocs[pGene->GOAssocIdx + NthTerm - 1];
return(&m_pszGOIdents[pAssoc->IdentIdx]);
}

int 
CGOAssocs::GetGOAttribs(char *pszGeneName,int NthTerm, etGOEvidence *pEvidence,etGOobjType *pType,etGOQual *pQual,etGOaspect *pAspect)
{
tsGOGene *pGene;
tsGOAssoc *pAssoc;

// firstly locate the gene
if((pGene = LocateGene(pszGeneName)) == NULL)
	return(eBSFerrGOID);
if(NthTerm < 1 || NthTerm > pGene->NumGOAssoc)
	return(eBSFerrGOID);
pAssoc = &m_pAssocs[pGene->GOAssocIdx + NthTerm - 1];
if(pEvidence != NULL)
	*pEvidence = (etGOEvidence)pAssoc->Evidence;
if(pType != NULL)
	*pType = (etGOobjType)pGene->Type;
if(pQual != NULL)
	*pQual = (etGOQual)pAssoc->Qual;
if(pAspect != NULL)
	*pAspect = (etGOaspect)pAssoc->Aspect;

return(eBSFSuccess);
}

// SortDedupeGOIDs
// Used to sort by GOID
int 
CGOAssocs::SortDedupeGOIDs( const void *arg1, const void *arg2)
{
tsDedupedGOID *pEl1 = (tsDedupedGOID *)arg1;
tsDedupedGOID *pEl2 = (tsDedupedGOID *)arg2;
char *pszGOID1 = (char *)pEl1->szGOID;
char *pszGOID2 = (char *)pEl2->szGOID;
char c1 = tolower(*pszGOID1);
char c2 = tolower(*pszGOID2);
if(c1 < c2)
	return(-1);
if(c1 > c2)
	return(1);
return(stricmp(pszGOID1,pszGOID2));
}


// SortGeneGoItems
// Used to sort by gene names and ObjType
int 
CGOAssocs::SortGeneGoItems( const void *arg1, const void *arg2)
{
int Cmp;
tsGeneGOitem *pEl1 = (tsGeneGOitem *)arg1;
tsGeneGOitem *pEl2 = (tsGeneGOitem *)arg2;
char *pszName1 = (char *)pEl1->pszGeneName;
char *pszName2 = (char *)pEl2->pszGeneName;
char c1 = tolower(*pszName1);
char c2 = tolower(*pszName2);
if(c1 < c2)
	return(-1);
if(c1 > c2)
	return(1);
if((Cmp = stricmp(pszName1,pszName2))!=0)
	return(Cmp);
if(pEl1->ObjType < pEl2->ObjType)
	return(-1);
if(pEl1->ObjType > pEl2->ObjType)
	return(1);
return(0);
}


// SortMapGenes
// Used to sort by gene names
int 
CGOAssocs::SortMapGenes( const void *arg1, const void *arg2)
{
int Cmp;
tsGeneMap *pEl1 = (tsGeneMap *)arg1;
tsGeneMap *pEl2 = (tsGeneMap *)arg2;
char *pszName1 = (char *)pEl1->szFromGene;
char *pszName2 = (char *)pEl2->szFromGene;
char c1 = tolower(*pszName1);
char c2 = tolower(*pszName2);
if(c1 < c2)
	return(-1);
if(c1 > c2)
	return(1);
if((Cmp = stricmp(pszName1,pszName2))!=0)
	return(Cmp);

pszName1 = (char *)pEl1->szToGene;
pszName2 = (char *)pEl2->szToGene;
c1 = tolower(*pszName1);
c2 = tolower(*pszName2);
if(c1 < c2)
	return(-1);
if(c1 > c2)
	return(1);
return(stricmp(pszName1,pszName2));
}

// SortGeneFilters
// Used to sort by gene names
int 
CGOAssocs::SortGeneFilters( const void *arg1, const void *arg2)
{
tsGeneFilter *pEl1 = (tsGeneFilter *)arg1;
tsGeneFilter *pEl2 = (tsGeneFilter *)arg2;
char *pszName1 = (char *)pEl1->szGene;
char *pszName2 = (char *)pEl2->szGene;
char c1 = tolower(*pszName1);
char c2 = tolower(*pszName2);
if(c1 < c2)
	return(-1);
if(c1 > c2)
	return(1);
return(stricmp(pszName1,pszName2));
}

// SortGenes
// Used to sort by gene names
int 
CGOAssocs::SortGenes( const void *arg1, const void *arg2)
{
tsGOGene *pEl1 = (tsGOGene *)arg1;
tsGOGene *pEl2 = (tsGOGene *)arg2;
char *pszName1 = (char *)pEl1->Name.ptr;
char *pszName2 = (char *)pEl2->Name.ptr;
char c1 = tolower(*pszName1);
char c2 = tolower(*pszName2);
if(c1 < c2)
	return(-1);
if(c1 > c2)
	return(1);
return(stricmp(pszName1,pszName2));
}




