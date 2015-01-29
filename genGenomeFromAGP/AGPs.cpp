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

#include "AGPs.h"

CAGPs::CAGPs(void)
{
m_pAGPentries = NULL;
Reset();
}

CAGPs::~CAGPs(void)
{
if(m_pAGPentries != NULL)
	{
#ifdef _WIN32
	free(m_pAGPentries);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAGPentries != MAP_FAILED)
		munmap(m_pAGPentries,m_AllocdAGPMem);
#endif
	}
}

char *
CAGPs::TrimWhitespace(char *pTxt)
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
CAGPs::Reset(void)
{
if(m_pAGPentries != NULL)
	{
#ifdef _WIN32
	free(m_pAGPentries);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAGPentries != MAP_FAILED)
		munmap(m_pAGPentries,m_AllocdAGPMem);
#endif
	m_pAGPentries = NULL;
	}
m_AllocdAGPMem = 0;
m_NumAGPentries = 0;
}


int 
CAGPs::LoadAGPs(char *pszAGPFile)
{
int Rslt;
FILE *pAGPStream;
int Cnt;
int Idx;
int scannxt;
char *pTxt;
char szLineBuff[cLineBuffLen];
char cEntryType;
char szGapType[30];
char szLinkage[10];
char szOrientation[10];
tsAGPentry Entry;
tsAGPentry *pEntry;
size_t memreq;

if(pszAGPFile == NULL || *pszAGPFile == '\0')
	return(eBSFerrParams);

memreq = (size_t)sizeof(tsAGPentry) * cAllocAGPentries;

#ifdef _WIN32
m_pAGPentries = (tsAGPentry *) malloc(memreq);	// initial and perhaps the only allocation

if(m_pAGPentries == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessAGP: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pAGPentries = (tsAGPentry *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pAGPentries == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessAGP: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pAGPentries = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_NumAGPentries = 0;
m_AllocAGPentries = cAllocAGPentries;
m_AllocdAGPMem = memreq;

if((pAGPStream = fopen(pszAGPFile,"r"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to fopen AGP format file %s error: %s",pszAGPFile,strerror(errno));
	return(eBSFerrOpnFile);
	}

Rslt = eBSFSuccess;
while(fgets(szLineBuff,sizeof(szLineBuff),pAGPStream)!= NULL)
	{
	pTxt = TrimWhitespace(szLineBuff);
	if(*pTxt=='\0' || *pTxt=='#')	// simply slough lines which were just whitespace or start with '#'
		continue;
	Cnt = sscanf(pTxt,"%50s %d %d %d %c %n",Entry.szObjIdent,&Entry.Start,&Entry.End,&Entry.PartNum,&cEntryType,&scannxt);
	if(Cnt != 5)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst parsing - %s - %s",pszAGPFile, pTxt);
		Rslt = eBSFerrParse; // if at least one feature processed then treat rubbish as error
		break;
		}
	switch(tolower(cEntryType)) {
		case 'a':
			Entry.Type = eAGPSSA;			// Active Finishing
			break;
		case 'd':
			Entry.Type = eAGPSSD;		// Draft HTG (often phase1 and phase2 are called Draft, whether or not they have the draft keyword).
			break;
		case 'f':
			Entry.Type = eAGPSSF;		// Finished HTG (phase 3)
			break;
		case 'g':
			Entry.Type = eAGPSSG;		// Whole Genome Finishing
			break;
		case 'n':
			Entry.Type = eAGPSSN;		// gap with specified size
			break;
		case 'o':
			Entry.Type = eAGPSSO;		// Other sequence (typically means no HTG keyword)
			break;
		case 'p':
			Entry.Type = eAGPSSP;		// Pre Draft
			break;
		case 'u':
			Entry.Type = eAGPSSU;		// gap of unknown size, typically defaulting to predefined values.
			break;
		case 'w':
			Entry.Type = eAGPSSW;		// WGS contig
			break;
		default:
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst parsing entry type '%c' - %s",cEntryType,pszAGPFile);
			Rslt = eBSFerrParse; // if at least one feature processed then treat rubbish as error
			break;
		}
	if(Rslt != eBSFSuccess)
		break;

	if(Entry.Type == eAGPSSN)		// gap?
		{
		Cnt = sscanf(&pTxt[scannxt],"%d %30s %10s", &Entry.Gap.GapLen, szGapType, szLinkage);
		if(Cnt != 3)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst parsing line '%s' - %s",pTxt, pszAGPFile);
			Rslt = eBSFerrParse; // if at least one feature processed then treat rubbish as error
			break;
			}
		if(!stricmp(szLinkage,"yes"))
			Entry.Gap.bLinkage = true;
		else
			if(!stricmp(szLinkage,"no"))
				Entry.Gap.bLinkage = false;
			else
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst parsing line gap linkage '%s' - %s", pTxt, pszAGPFile);
				Rslt = eBSFerrParse; // if at least one feature processed then treat rubbish as error
				break;
				}


		if(!stricmp(szGapType,"fragment"))
			Entry.Gap.GapType = AGPfragment;
		else
			if(!stricmp(szGapType,"clone"))
				Entry.Gap.GapType = AGPclone;
			else
				if(!stricmp(szGapType,"contig"))
					Entry.Gap.GapType = AGPcontig;
				else
					if(!stricmp(szGapType,"centromere"))
						Entry.Gap.GapType = AGPcentromere;
					else
						if(!stricmp(szGapType,"hort_arm"))
							Entry.Gap.GapType = AGPshort_arm;
						else
							if(!stricmp(szGapType,"heterochromatin"))
								Entry.Gap.GapType = AGPheterochromatin;
							else
								if(!stricmp(szGapType,"telomere"))
									Entry.Gap.GapType = AGPtelomere;
								else
									if(!stricmp(szGapType,"repeat"))
										Entry.Gap.GapType = AGPrepeat;
									else
										{
										gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst parsing gap type in line '%s' - %s",pTxt,pszAGPFile);
										Rslt = eBSFerrParse; // if at least one feature processed then treat rubbish as error
										break;
										}

		}
	else
		{
		Cnt = sscanf(&pTxt[scannxt],"%50s %d %d %10s", Entry.Comp.szCompID, &Entry.Comp.Start, &Entry.Comp.End, szOrientation);
		if(Cnt != 4)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst parsing line '%s' - %s",pTxt, pszAGPFile);
			Rslt = eBSFerrParse; // if at least one feature processed then treat rubbish as error
			break;
			}
		if(!stricmp(szOrientation,"+"))
			Entry.Comp.Orientation = eAGPOplus;
		else
			if(!stricmp(szOrientation,"-"))
				Entry.Comp.Orientation = eAGPOminus;
		else
			if(!stricmp(szOrientation,"0"))
				Entry.Comp.Orientation = eAGPOunknown;
		else
			if(!stricmp(szOrientation,"na"))
				Entry.Comp.Orientation = eAGPna;
			else
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst parsing orientation in line '%s' - %s",pTxt,pszAGPFile);
				Rslt = eBSFerrParse; // if at least one feature processed then treat rubbish as error
				break;
				}

		}
	if(m_NumAGPentries == m_AllocAGPentries)
		{
		tsAGPentry *pNewEntry;
		memreq = (size_t)sizeof(tsAGPentry) * (m_AllocAGPentries + cAllocAGPentries);
#ifdef _WIN32
		pNewEntry = (tsAGPentry *) realloc(m_pAGPentries,memreq);
#else
		pNewEntry = (tsAGPentry *)mremap(m_pAGPentries,m_AllocdAGPMem,memreq,MREMAP_MAYMOVE);
		if(pNewEntry == MAP_FAILED)
			pNewEntry = NULL;
#endif
		if(pNewEntry == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddMHitReads: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
			return(eBSFerrMem);
			}
		m_pAGPentries = pNewEntry;
		m_AllocAGPentries += cAllocAGPentries;
		m_AllocdAGPMem = memreq;
		}
	m_pAGPentries[m_NumAGPentries++] = Entry;
	}

fclose(pAGPStream);
if(Rslt == eBSFSuccess)
	{
	if(m_NumAGPentries > 1)
		qsort(m_pAGPentries,m_NumAGPentries,sizeof(tsAGPentry),SortAGPentries);
	pEntry = &m_pAGPentries[0];
	for(Idx = 1; Idx <= (int)m_NumAGPentries; Idx++, pEntry++)
		pEntry->EntryID = Idx;
	}
return(Rslt);
}

// locate CompID
// currently a linear search as number of entries expected to be relatively small, maybe a few hundred at most
tsAGPentry *
CAGPs::LocateCompID(char *pszCompID)
{
UINT32 Idx;
tsAGPentry *pEntry = m_pAGPentries;
if(!m_NumAGPentries || pszCompID == NULL || pszCompID[0] == '\0')
	return(NULL);
for(Idx = 0; Idx < m_NumAGPentries; Idx++, pEntry++)
	{
	if(pEntry->Type != eAGPSSN && !stricmp(pszCompID,pEntry->Comp.szCompID))
		return(pEntry);
	}
return(NULL);
}

// returns next entry (NULL if none) after pPrev. If pPrev == NULL then returns first
tsAGPentry *
CAGPs::Next(tsAGPentry *pPrev)
{
if(!m_NumAGPentries)
	return(NULL);
if(pPrev == NULL)
	return(&m_pAGPentries[0]);
if(pPrev->EntryID < 1 || pPrev->EntryID >= m_NumAGPentries)
	return(NULL);
return(&m_pAGPentries[pPrev->EntryID]);
}

tsAGPentry *
CAGPs::Entry(UINT32 EntryID)
{
if(!m_NumAGPentries)
	return(NULL);
if(EntryID < 1 || EntryID > m_NumAGPentries)
	return(NULL);
return(&m_pAGPentries[EntryID-1]);
}

// SortAGPentries
// used to sort AGP entries by object name, start, ascending
int
CAGPs::SortAGPentries(const void *pEl1,const void *pEl2)
{
int Cmp;
tsAGPentry *pC1 = (tsAGPentry *)pEl1;
tsAGPentry *pC2 = (tsAGPentry *)pEl2;
if((Cmp = stricmp(pC1->szObjIdent,pC2->szObjIdent))!=0)
	return(Cmp);
return(pC1->Start - pC2->Start);
}