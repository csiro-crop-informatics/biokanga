#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stdafx.h"

#if _WIN32
#include "../conservlib/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

#include "FilterRefIDs.h"

CFilterRefIDs::CFilterRefIDs(void)
{
m_pFilterRefIDs = NULL;
m_NumFilterRefIDs = 0;
m_AllocdFilterRefIDs = 0;
}

CFilterRefIDs::~CFilterRefIDs(void)
{
if(m_pFilterRefIDs != NULL)
	delete m_pFilterRefIDs;
}

void
CFilterRefIDs::Reset(void)
{
if(m_pFilterRefIDs != NULL)
	{
	delete m_pFilterRefIDs;
	m_pFilterRefIDs= NULL;
	}
m_NumFilterRefIDs = 0;
m_AllocdFilterRefIDs = 0;
}


int
CFilterRefIDs::Open(char *pszFile)
{
int Rslt;
int *pRefIDs;
int RefID;
int NumFields;

if((Rslt=CCSVFile::Open(pszFile))!=eBSFSuccess)
	{
	AddErrMsg("CFilterRefIDs::Open","Unable to open '%s'",pszFile);
	return(Rslt);
	}

// now parse each RefID which is to be filtered out
while((Rslt=NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = GetCurFields();
	if(NumFields < 1)
		{
		AddErrMsg("CFilterRefIDs::Open","Expected at least 1 fields in '%s', GetCurFields() returned '%d'",pszFile,NumFields);
		return(eBSFerrFieldCnt);
		}
	GetInt(1,&RefID);

	if(m_pFilterRefIDs == NULL || m_NumFilterRefIDs == m_AllocdFilterRefIDs)
		{
		if((pRefIDs = new int [m_AllocdFilterRefIDs + cFiltRefIDsAllocChunk])==NULL)
			{
			AddErrMsg("CFilterRefIDs::Open","Unable to alloc memory for filtered RefIDs");
			return(eBSFerrMem);
			}

		if(m_pFilterRefIDs != NULL)
			{
			memcpy(pRefIDs,m_pFilterRefIDs,sizeof(int) * m_NumFilterRefIDs);
			delete m_pFilterRefIDs;
			}
		else
			{
			m_AllocdFilterRefIDs = 0;
			m_NumFilterRefIDs = 0;
			}
		m_AllocdFilterRefIDs += cFiltRefIDsAllocChunk;
		m_pFilterRefIDs = pRefIDs;
		}
	m_pFilterRefIDs[m_NumFilterRefIDs++] = RefID;
	}

CCSVFile::Close();
if(m_NumFilterRefIDs > 1)
	qsort(m_pFilterRefIDs,m_NumFilterRefIDs,sizeof(int),SortRefIDs);

return(m_NumFilterRefIDs);
}

// returns true if RefID is in RefIDs loaded
bool 
CFilterRefIDs::Locate(int RefID)
{
int *pProbe;
int Lo,Mid,Hi;	// search limits
Lo = 0; Hi = m_NumFilterRefIDs-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pProbe = &m_pFilterRefIDs[Mid];
	if(RefID == *pProbe)
		return(true);
	if(RefID < *pProbe)	
		{
		Hi = Mid - 1;
		continue;
		}
	Lo = Mid + 1;
	}
return(NULL);
}


int 
CFilterRefIDs::SortRefIDs( const void *arg1, const void *arg2)
{
return(*(int *)arg2 - *(int *)arg1);
}



