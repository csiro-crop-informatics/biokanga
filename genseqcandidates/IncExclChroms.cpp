#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

#include "IncExclChroms.h"


CIncExclChroms::CIncExclChroms()
{
int Idx;
#ifdef _WIN32
for(Idx = 0; Idx < cMaxIncludeChroms; Idx++)
	m_IncludeChromsRE[Idx] = NULL;
for(Idx = 0; Idx < cMaxExcludeChroms; Idx++)
	m_ExcludeChromsRE[Idx] = NULL;
#endif
m_AllocdIncExclChroms = 0;
m_IncExclChroms = 0;
m_pIncExclChroms = NULL;
m_pLastmatch = NULL;
m_NumIncludeChroms = 0;
m_NumExcludeChroms = 0;
}

CIncExclChroms::~CIncExclChroms()
{
int Idx;
#ifdef _WIN32
for(Idx = 0; Idx < cMaxIncludeChroms; Idx++)
	{
	if(m_IncludeChromsRE[Idx] != NULL)
		delete m_IncludeChromsRE[Idx];
	}
#endif
if(m_pIncExclChroms != NULL)
	delete m_pIncExclChroms;
}

void CIncExclChroms::Reset(void)
{
int Idx;
#ifdef _WIN32
for(Idx = 0; Idx < cMaxIncludeChroms; Idx++)
	{
	if(m_IncludeChromsRE[Idx] != NULL)
		{
		delete m_IncludeChromsRE[Idx];
		m_IncludeChromsRE[Idx] = NULL;
		}
	}
for(Idx = 0; Idx < cMaxExcludeChroms; Idx++)
	{
	if(m_ExcludeChromsRE[Idx] != NULL)
		{
		delete m_ExcludeChromsRE[Idx];
		m_ExcludeChromsRE[Idx] = NULL;
		}
	}
#endif
if(m_pIncExclChroms != NULL)
	{
	delete m_pIncExclChroms;
	m_pIncExclChroms = NULL;
	}
m_AllocdIncExclChroms = 0;
m_IncExclChroms = 0;
m_pLastmatch = NULL;
m_NumIncludeChroms = 0;
m_NumExcludeChroms = 0;
}

int
CIncExclChroms::InitChromExclusion(int	NumIncludeChroms,	// number of chromosome regular expressions to include
		char **ppszIncludeChroms,	// array of include chromosome regular expressions
		int	NumExcludeChroms,		// number of chromosome expressions to exclude
		char **ppszExcludeChroms)	// array of exclude chromosome regular expressions
{
int Idx;
Reset();

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];		// to hold RegErr as textual representation ==> regerror();
#endif

m_NumIncludeChroms = NumIncludeChroms;
m_NumExcludeChroms = NumExcludeChroms;

if(NumIncludeChroms == 0 && NumExcludeChroms == 0)
	return(eBSFSuccess);

#ifdef _WIN32
try {
	for(Idx=0;Idx < NumIncludeChroms;Idx++)
		{
		m_IncludeChromsRE[Idx] = new Regexp();
		m_IncludeChromsRE[Idx]->Parse((const STRCHAR *)ppszIncludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include regexpr chrom '%s'",ppszIncludeChroms[Idx]);
	Reset();
	return(eBSFerrMem);
	}
try {
	for(Idx=0;Idx < NumExcludeChroms;Idx++)
		{
		m_ExcludeChromsRE[Idx] = new Regexp();
		m_ExcludeChromsRE[Idx]->Parse((const STRCHAR *)ppszExcludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude regexpr chrom '%s'",ppszExcludeChroms[Idx]);
	Reset();
	return(eBSFerrMem);
	}

#else
for(Idx=0;Idx < NumIncludeChroms;Idx++)
	{
	RegErr=regcomp(&m_IncludeChromsRE[Idx],ppszIncludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&m_IncludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include chrom '%s' error: %s",ppszIncludeChroms[Idx],szRegErr);
		Reset();
		return(eBSFerrMem);
		}
	}
for(Idx=0;Idx < NumExcludeChroms;Idx++)
	{
	RegErr = regcomp(&m_ExcludeChromsRE[Idx],ppszExcludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&m_ExcludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude chrom '%s' error: %s",ppszExcludeChroms[Idx],szRegErr);
		Reset();
		return(eBSFerrMem);
		}
	}
#endif
return(eBSFSuccess);
}

//IncludeThisChrom
//Returns 0 if to be excluded, > 1 as the chromid when to be included or < 0 if error
//
int
CIncExclChroms::IncludeThisChrom(char *pszChrom)
{
int ChromIdx;
bool bProcChrom;
tsIncExclChrom *pIncExcl;

// if no include/exclude then chromosome is to be included
if(m_NumIncludeChroms == 0 && m_NumExcludeChroms == 0)
	{
	if((pIncExcl = LocateChrom(pszChrom))!=NULL)
		return(pIncExcl->ChromID);
	return(AddChrom(pszChrom,true));
	}

// optimisation is to check if chromosome previously known to be included/excluded
if((pIncExcl = LocateChrom(pszChrom))!=NULL)
	return((int)pIncExcl->bInc ? pIncExcl->ChromID : 0);

// not previously processed... 
#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
int RegErr;					// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

// check if to be excluded
bProcChrom = true;
for(ChromIdx = 0; ChromIdx < m_NumExcludeChroms; ChromIdx++)
	{
#ifdef _WIN32	
	if(m_ExcludeChromsRE[ChromIdx]->Match((const STRCHAR *)pszChrom,&mc))
#else
	if(!regexec(&m_ExcludeChromsRE[ChromIdx],pszChrom,1,&mc,0))
#endif
		{
		bProcChrom = false;
		break;
		}
	}

			// to be included?
if(bProcChrom && m_NumIncludeChroms)
	{
	bProcChrom = false;
	for(ChromIdx = 0; ChromIdx < m_NumIncludeChroms; ChromIdx++)
		{
#ifdef _WIN32
		if(m_IncludeChromsRE[ChromIdx]->Match((const STRCHAR *)pszChrom,&mc))
#else
		if(!regexec(&m_IncludeChromsRE[ChromIdx],pszChrom,1,&mc,0))
#endif
			{									
			bProcChrom = true;
			break;
			}
		}
	}

ChromIdx = AddChrom(pszChrom,bProcChrom);
return(bProcChrom ? ChromIdx : 0);
}

tsIncExclChrom *
CIncExclChroms::LocateChrom(char *pszChrom)
{
tsIncExclChrom *pTmp;
int Idx;
if(m_pIncExclChroms == NULL || m_IncExclChroms == 0)
	return(NULL);
if(m_pLastmatch != NULL && !stricmp(m_pLastmatch->szChrom,pszChrom))
	return(m_pLastmatch);

pTmp = m_pIncExclChroms;
for(Idx = 0; Idx < m_IncExclChroms; Idx++, pTmp++)
	if(!stricmp(pTmp->szChrom,pszChrom))
		{
		m_pLastmatch = pTmp;
		return(pTmp);
		}
return(NULL);
}

char *
CIncExclChroms::LocateChrom(int ChromID)	// returns chrom corresponding to specified ChromID
{
if(m_pIncExclChroms == NULL || ChromID < 1 || ChromID > m_IncExclChroms)
	return(NULL);
return(m_pIncExclChroms[ChromID-1].szChrom);
}


int					// returns total number of chromosomes added, becomes the chrom identifier! ( < 0 if error)
CIncExclChroms::AddChrom(char *pszChrom,bool bProcChrom) // caches chrom processing state
{
tsIncExclChrom *pTmp;
if(m_pIncExclChroms == NULL || m_IncExclChroms >= m_AllocdIncExclChroms)
	{
	if(m_pIncExclChroms == NULL)
		{
		m_IncExclChroms = 0;
		m_AllocdIncExclChroms = 0;
		}
	if((pTmp = new tsIncExclChrom [m_AllocdIncExclChroms + cAllocIncExclChroms])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding IncExcl chromosomes",sizeof(tsIncExclChrom) * (m_AllocdIncExclChroms + cAllocIncExclChroms));
		return(eBSFerrMem);
		}
	if(m_pIncExclChroms != NULL)
		{
		memcpy(pTmp,m_pIncExclChroms,sizeof(tsIncExclChrom) * m_IncExclChroms);
		delete m_pIncExclChroms;
		}
	m_AllocdIncExclChroms += cAllocIncExclChroms;
	m_pIncExclChroms = pTmp;
	}
pTmp = &m_pIncExclChroms[m_IncExclChroms++];
pTmp->ChromID = m_IncExclChroms;
pTmp->bInc = bProcChrom;
strncpy(pTmp->szChrom,pszChrom,cMaxDatasetSpeciesChrom);
pTmp->szChrom[cMaxDatasetSpeciesChrom-1] = '\0';
m_pLastmatch = pTmp;
return(m_IncExclChroms);
}


