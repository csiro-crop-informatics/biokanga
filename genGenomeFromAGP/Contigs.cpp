#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include <process.h>
#include "../libbiokanga/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libbiokanga/commhdrs.h"
#endif

#include "Contigs.h"

CContigs::CContigs(void)
{
m_pContigSeqs = NULL;
m_pContigEntries = NULL;
m_pContigChunk = NULL;
Reset();
}

CContigs::~CContigs(void)
{
if(m_pContigSeqs != NULL)
	{
#ifdef _WIN32
	free(m_pContigSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pContigSeqs != MAP_FAILED)
		munmap(m_pContigSeqs,m_AllocdContigSeqMem);
#endif
	}
if(m_pContigEntries != NULL)
	{
#ifdef _WIN32
	free(m_pContigEntries);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pContigEntries != MAP_FAILED)
		munmap(m_pContigSeqs,m_AllocContigEntriesMem);
#endif
	}
if(m_pContigChunk != NULL)
	delete m_pContigChunk;
}

void
CContigs::Reset(void)
{
if(m_pContigSeqs != NULL)
	{
#ifdef _WIN32
	free(m_pContigSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pContigSeqs != MAP_FAILED)
		munmap(m_pContigSeqs,m_AllocdContigSeqMem);
#endif
	m_pContigSeqs = NULL;
	}
if(m_pContigEntries != NULL)
	{
#ifdef _WIN32
	free(m_pContigEntries);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pContigEntries != MAP_FAILED)
		munmap(m_pContigSeqs,m_AllocContigEntriesMem);
#endif
	m_pContigEntries = NULL;
	}
if(m_pContigChunk != NULL)
	{
	delete m_pContigChunk;
	m_pContigChunk = NULL;
	}

m_AllocdContigSeqMem = 0;			// total memory allocated to m_pContigSeqs 
m_ContigSeqsOfs = 0;				// offset in m_pContigSeqs at which to next write
m_NumContigEntries = 0;				// number of contig entries
m_AllocContigEntriesMem = 0;	    // total memory allocated for contig entries
}



// LoadContigFiles
// Parse input fasta format file(s) containing contig sequences
teBSFrsltCodes  
CContigs::LoadContigs(char *pszContigFiles)
{
teBSFrsltCodes Rslt;
char *pszInfile;
int NumInputFilesProcessed;
CSimpleGlob glob(SG_GLOB_FULLSORT);

glob.Init();
if(glob.Add((const char *)pszContigFiles) < SG_SUCCESS)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszContigFiles);
	return(eBSFerrOpnFile);	// treat as though unable to open file
	}

if(glob.FileCount() <= 0)
	{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigs: Unable to locate any source multifasta contig file matching '%s",pszContigFiles);
	Reset();
	return(eBSFerrOpnFile);
	}

Rslt = eBSFSuccess;
NumInputFilesProcessed = 0;
for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
	{
	pszInfile = glob.File(FileID);
	NumInputFilesProcessed += 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing multifasta file for contigs '%s'",pszInfile);
	Rslt = LoadContigFile(pszInfile);
	if(Rslt != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigs: Load failed for file '%s'\n",pszInfile);
		Reset();
		return(Rslt);
		}
	}
if(m_NumContigEntries > 1)
	qsort(m_pContigEntries,m_NumContigEntries,sizeof(tsContigEntry),SortContigs);
return(eBSFSuccess);
}

// LoadContigFile
// Parse input fasta format file containing contig sequences
teBSFrsltCodes  
CContigs::LoadContigFile(char *pszFile)
{
CFasta Fasta;
tsContigEntry *pCurContig;
char szName[cMaxDatasetSpeciesChrom+1];			// to hold contig name
char szDescription[cBSFDescriptionSize+1];			// to hold contig descriptor
int SeqLen;
int Descrlen;
bool bEntryStarted;
int Rslt;
size_t memreq;
UINT8 *pAllocSeqs;

// alloc for contig sequence chunk buffering
if(m_pContigChunk == NULL)
	{
	if((m_pContigChunk = new unsigned char [cContigSeqChunk+1]) == NULL)	// why the additional 1 is required????
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigFile: Unable to allocate memory (%d bytes) for contig sequence chunks",cContigSeqChunk+1);
		Reset();
		return(eBSFerrMem);
		}
	}

if(m_pContigEntries == NULL)
	{
	// initial allocation for contig entries
	memreq = (size_t)sizeof(tsContigEntry) * cContigAllocEntries;

#ifdef _WIN32
	m_pContigEntries = (tsContigEntry *) malloc(memreq);	// initial and perhaps the only allocation

	if(m_pContigEntries == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigFile: Memory allocation of %lld bytes for contig entries - %s",(INT64)memreq,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
		// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pContigEntries = (tsContigEntry *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pContigEntries == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigFile: Memory allocation of %lld bytes through mmap() for contig entries failed - %s",(INT64)memreq,strerror(errno));
		m_pContigEntries = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_NumContigEntries = 0;
	m_AllocContigEntries = cContigAllocEntries;
	m_AllocContigEntriesMem = memreq;
	}

if(m_pContigSeqs == NULL)
	{
	// initial allocation for contig sequences
	memreq = (size_t)sizeof(UINT8) * cContigSeqAlloc;

#ifdef _WIN32
	m_pContigSeqs = (UINT8 *) malloc(memreq);	// initial and perhaps the only allocation

	if(m_pContigSeqs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigFile: Memory allocation of %lld bytes for contig sequences - %s",(INT64)memreq,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
		// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pContigSeqs = (UINT8 *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pContigSeqs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigFile: Memory allocation of %lld bytes through mmap() for contig sequences failed - %s",(INT64)memreq,strerror(errno));
		m_pContigSeqs = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_ContigSeqsOfs = 0;
	m_AllocdContigSeqMem = memreq;
	}


gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadContigFile: Processing contigs from %s..",pszFile);
if((Rslt=(teBSFrsltCodes)Fasta.Open(pszFile,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigFile: Unable to open '%s' [%s] %s",pszFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	Reset();
	return(eBSFerrOpnFile);
	}
bEntryStarted = false;
while((Rslt = SeqLen = Fasta.ReadSequence(m_pContigChunk,cContigSeqChunk,true,false)) > eBSFSuccess)
	{

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		if(m_NumContigEntries == m_AllocContigEntries)
			{
			memreq = sizeof(tsContigEntry) * (m_AllocContigEntries + cContigAllocEntries); 
#ifdef _WIN32
			pCurContig = (tsContigEntry *) realloc(m_pContigEntries,memreq);
#else
			pCurContig = (tsContigEntry *)mremap(m_pContigEntries,m_AllocContigEntriesMem,memreq,MREMAP_MAYMOVE);
			if(pCurContig == MAP_FAILED)
				pCurContig = NULL;
#endif
			if(pCurContig == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessContigFile: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
				return(eBSFerrMem);
				}
			m_AllocContigEntriesMem = memreq;
			m_AllocContigEntries += cContigAllocEntries;
			m_pContigEntries = pCurContig;
			}

		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		// An assumption - will one day bite real hard - is that the
		// fasta descriptor line starts with some form of unique identifier.
		// Use this identifier as the contig entry name.
		if(sscanf(szDescription," %s[ ,]",szName) != 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessContigFile: Unable to parse contig identifer from descriptor '%s'",szDescription);
			Reset();
			return(eBSFerrFastaDescr);
			}

		pCurContig = &m_pContigEntries[m_NumContigEntries];
	    memset(pCurContig,0,sizeof(tsContigEntry));
		strncpy(pCurContig->szDescr,szDescription,sizeof(pCurContig->szDescr));
		pCurContig->szDescr[sizeof(pCurContig->szDescr)-1] = '\0';
		strncpy(pCurContig->szContig,szName,sizeof(pCurContig->szContig));
		pCurContig->szContig[sizeof(pCurContig->szContig)-1] = '\0';
		pCurContig->ContigID = ++m_NumContigEntries;
		pCurContig->ContigSeqOfs = (UINT32)m_ContigSeqsOfs;
		bEntryStarted = true;

		continue;
		}
	else
		if(!bEntryStarted)	// if there was no descriptor then thats a fatal error...
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigFile: Expected sequence to be preceded by descriptor containing contig identifier in contig file %s",pszFile);
			Reset();
			return(eBSFerrFastaDescr);
			}

	// have sequence chunk, copy to m_pContigSeqs
	if((m_AllocdContigSeqMem - m_ContigSeqsOfs) < (size_t)SeqLen)		// realloc as may be required
		{
		memreq = m_AllocdContigSeqMem + (sizeof(UINT8) * cContigSeqAlloc); 
#ifdef _WIN32
		pAllocSeqs = (UINT8 *) realloc(m_pContigSeqs,memreq);
#else
		pAllocSeqs = (UINT8 *)mremap(m_pContigSeqs,m_AllocdContigSeqMem,memreq,MREMAP_MAYMOVE);
		if(pAllocSeqs == MAP_FAILED)
			pAllocSeqs = NULL;
#endif
		if(pAllocSeqs == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessContigFile: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
			return(eBSFerrMem);
			}
		m_AllocdContigSeqMem = memreq;
		m_pContigSeqs = pAllocSeqs;
		}
	pCurContig = &m_pContigEntries[m_NumContigEntries-1];
	memcpy(&m_pContigSeqs[m_ContigSeqsOfs],m_pContigChunk,SeqLen);
	m_ContigSeqsOfs += SeqLen;
	pCurContig->ContigLen += SeqLen;
	}
if(Rslt < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessContigFile [%s] %s",Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	Reset();
	return((teBSFrsltCodes)Rslt);
	}
return(eBSFSuccess);
}

tsContigEntry *
CContigs::LocateContig(char *pszContig)
{
int Cmp;
int Lo,Mid,Hi;	// search limits
tsContigEntry *pContig;

if(!m_NumContigEntries)
	return(NULL);

Lo = 0; Hi = m_NumContigEntries-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pContig = &m_pContigEntries[Mid];
	if((Cmp = stricmp(pszContig,pContig->szContig))==0)
		return(pContig);
	if(Cmp < 0)	
		Hi = Mid - 1;
	else	
		Lo = Mid + 1;
	}
return(NULL);
}

// returns ptr to sequence for pszContig starting at Start (1..ContigLen)
UINT8 *
CContigs::LocateSeq(char *pszContig,UINT32 Start) 
{
tsContigEntry *pContig;
if((pContig = LocateContig(pszContig))==NULL)
	return(NULL);
if(pContig->ContigLen < Start)
	return(NULL);
return(&m_pContigSeqs[pContig->ContigSeqOfs + Start - 1]);
}

int
CContigs::SortContigs(const void *pEl1,const void *pEl2)
{
tsContigEntry *pC1 = (tsContigEntry *)pEl1;
tsContigEntry *pC2 = (tsContigEntry *)pEl2;
return(stricmp(pC1->szContig,pC2->szContig));
}

