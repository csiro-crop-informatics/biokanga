#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include "../libbiokanga/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libbiokanga/commhdrs.h"
#endif

#include "ChromMap.h"

CChromMap::CChromMap(void)
{
m_pCSV = NULL;
m_pChromMaps = NULL;
Reset();
}

CChromMap::~CChromMap(void)
{
if(m_pCSV != NULL)
	delete m_pCSV;
if(m_pChromMaps != NULL)
	{
#ifdef _WIN32
	free(m_pChromMaps);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pChromMaps != MAP_FAILED)
		munmap(m_pChromMaps,m_AllocdChromMapsMem);
#endif
	}
}

char *
CChromMap::TrimWhitespace(char *pTxt)
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

void 
CChromMap::Reset(void)
{
if(m_pCSV != NULL)
	{
	delete m_pCSV;
	m_pCSV = NULL;
	}
if(m_pChromMaps != NULL)
	{
#ifdef _WIN32
	free(m_pChromMaps);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pChromMaps != MAP_FAILED)
		munmap(m_pChromMaps,m_AllocdChromMapsMem);
#endif
	m_pChromMaps = NULL;
	}
m_AllocdChromMapsMem = 0;
m_AllocdChromMaps = 0;	
m_CurNumChromMaps = 0;	
m_CurLineNum = 0;
}

int
CChromMap::LoadMap(char *pszChromMap)
{
int Rslt;
int NumFields;
size_t memreq;
tsChromMap ChromMap;
tsChromMap *pTmpAlloc;

char *pszContig;
char *pszChrom;
int Start;
int End;

if(pszChromMap == NULL || *pszChromMap == '\0')
	return(eBSFerrParams);

Reset();
m_pCSV = new CCSVFile;
if((m_pCSV->Open(pszChromMap))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open chrom map file %s",pszChromMap);
	return(eBSFerrOpnFile);
	}
strcpy(m_szMapFile,pszChromMap);

memreq = (size_t)sizeof(tsChromMap) * cAllocChromMaps;

#ifdef _WIN32
m_pChromMaps = (tsChromMap *) malloc(memreq);	// initial and perhaps the only allocation

if(m_pChromMaps == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadChromMap: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pChromMaps = (tsChromMap *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pChromMaps == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadChromMap::Load: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pChromMaps = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif

m_CurNumChromMaps = 0;
m_AllocdChromMaps = cAllocChromMaps;
m_AllocdChromMapsMem = memreq;
m_CurLineNum = 0;

while((Rslt=m_pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pCSV->GetCurFields();
	if(NumFields < 4)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least 4 fields in '%s', GetCurFields() returned '%d'",pszChromMap,NumFields);
		Reset();
		return(eBSFerrFieldCnt);
		}
	if(!m_CurNumChromMaps && m_pCSV->IsLikelyHeaderLine())
		continue;

	if(m_CurNumChromMaps == m_AllocdChromMaps)
		{
		memreq = (size_t)sizeof(tsChromMap) * (m_AllocdChromMaps + cAllocChromMaps);
#ifdef _WIN32
		pTmpAlloc = (tsChromMap *) realloc(m_pChromMaps,memreq);
#else
		pTmpAlloc = (tsChromMap *)mremap(m_pChromMaps,m_AllocdChromMapsMem,memreq,MREMAP_MAYMOVE);
		if(pTmpAlloc == MAP_FAILED)
			pTmpAlloc = NULL;
#endif
		if(pTmpAlloc == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
			Reset();	
			return(eBSFerrMem);
			}
		m_pChromMaps = pTmpAlloc;
		m_AllocdChromMapsMem = memreq;
		m_AllocdChromMaps += cAllocChromMaps;
		}


	// parse out the expected fields
	m_pCSV->GetText(1,&pszContig);
	m_pCSV->GetText(2,&pszChrom);
	m_pCSV->GetInt(3,&Start);
	m_pCSV->GetInt(4,&End);
	strcpy(ChromMap.szContig,pszContig);
	strcpy(ChromMap.szChrom,pszChrom);
	ChromMap.Start = Start;
	ChromMap.End = End;
	m_pChromMaps[m_CurNumChromMaps++]=ChromMap;
	}
delete m_pCSV;
m_pCSV = NULL;
if(Rslt == 0 && m_CurNumChromMaps > 1)
	qsort(m_pChromMaps,m_CurNumChromMaps,sizeof(tsChromMap),SortMapEntries);
return(Rslt == 0 ? m_CurNumChromMaps : Rslt);
}

int 
CChromMap::Map(char *pszContig,		// input: map this contig - return: chrom contig is on
			int *pStart,			// input: starting at this contig offset (1..n) - return: starts at this chrom loci (1..n)
			int *pEnd)				// input: and ending at this contig offset (start..n) - return: ends at this chrom loci (start..n)
{
tsChromMap *pMap;
if(pszContig == NULL || pszContig[0] == '\0' || (pStart != NULL && *pStart < 0) || (pEnd != NULL && *pEnd < 0))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Map: errors in function parameters");
	return(eBSFerrParams);
	}
if((pMap=LocateMapEntry(pszContig)) == NULL)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Map: unable to locate entry for '%s'",pszContig);
	return(eBSFerrEntry);
	}
strcpy(pszContig,pMap->szChrom);
if(pStart != NULL)
	*pStart = (*pStart + pMap->Start) - 1;
if(pEnd != NULL)
	*pEnd = (*pEnd + pMap->Start) - 1;
return(eBSFSuccess);
}

tsChromMap *
CChromMap::LocateMapEntry(char *pszContig)
{
int Cmp;
int Lo,Mid,Hi;	// search limits
tsChromMap *pContig;

if(!m_CurNumChromMaps)
	return(NULL);

Lo = 0; Hi = m_CurNumChromMaps-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pContig = &m_pChromMaps[Mid];
	if((Cmp = stricmp(pszContig,pContig->szContig))==0)
		return(pContig);
	if(Cmp < 0)	
		Hi = Mid - 1;
	else	
		Lo = Mid + 1;
	}
return(NULL);
}

// SortMapEntries
// used to sort entries by object name, chrom, start, ascending
int
CChromMap::SortMapEntries(const void *pEl1,const void *pEl2)
{
int Cmp;
tsChromMap *pC1 = (tsChromMap *)pEl1;
tsChromMap *pC2 = (tsChromMap *)pEl2;
if((Cmp = stricmp(pC1->szContig,pC2->szContig))!=0)
	return(Cmp);
if((Cmp = stricmp(pC1->szChrom,pC2->szChrom))!=0)
	return(Cmp);
if(pC1->Start > pC2->Start)
	return(1);
if(pC1->Start < pC2->Start)
	return(-1);
return(0);
}