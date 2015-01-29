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
#include "./commhdrs.h"
#endif


CFilterLoci::CFilterLoci(void)
{
m_pFilterLocii = NULL;
m_pLociChroms = NULL;
Reset();
}

CFilterLoci::~CFilterLoci(void)
{
if(m_pFilterLocii != NULL)
	delete m_pFilterLocii;
if(m_pLociChroms != NULL)
	delete m_pLociChroms;
}

void
CFilterLoci::Reset(void)
{
if(m_pFilterLocii != NULL)
	{
	delete m_pFilterLocii;
	m_pFilterLocii= NULL;
	}
if(m_pLociChroms != NULL)
	{
	delete m_pLociChroms;
	m_pLociChroms= NULL;
	}
m_NumFilterLocii = 0;
m_AllocdFilterLocii = 0;
m_MaxLociLen = 0;
m_CachLociChromID = 0;			// last located chromosome
m_NumLociChroms = 0;			// number of chromosomes
m_AllocdLociChroms = 0;			// number allocated
}


int
CFilterLoci::Load(char *pszFile)
{
int Rslt;
tsFilterLoci *pLocii;
int ChromID;
int StartLoci;
int EndLoci;
int LociLen;
int NumFields;
char *pszChrom;

if((Rslt=CCSVFile::Open(pszFile))!=eBSFSuccess)
	{
	AddErrMsg("CFilterLoci::Open","Unable to open '%s'",pszFile);
	return(Rslt);
	}

// now parse each loci which is to be filtered out
while((Rslt=NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = GetCurFields();
	if(NumFields < 1)
		{
		AddErrMsg("CFilterLoci::Open","Expected at least 6 fields in '%s', GetCurFields() returned '%d'",pszFile,NumFields);
		CCSVFile::Close();
		return(eBSFerrFieldCnt);
		}
	if(!m_NumFilterLocii && IsLikelyHeaderLine())
		continue;

	GetText(4,&pszChrom);
	GetInt(5,&StartLoci);
	GetInt(6,&EndLoci);
	
	if(m_pFilterLocii == NULL || m_NumFilterLocii == m_AllocdFilterLocii)
		{
		if((pLocii = new tsFilterLoci [m_AllocdFilterLocii + cFiltLociAllocChunk])==NULL)
			{
			AddErrMsg("CFilterLoci::Open","Unable to alloc memory for filtered Locii");
			CCSVFile::Close();		
			return(eBSFerrMem);
			}

		if(m_pFilterLocii != NULL)
			{
			if(m_NumFilterLocii)
				memmove(pLocii,m_pFilterLocii,sizeof(tsFilterLoci) * m_NumFilterLocii);
			delete m_pFilterLocii;
			}
		else
			{
			m_AllocdFilterLocii = 0;
			m_NumFilterLocii = 0;
			}
		m_AllocdFilterLocii += cFiltLociAllocChunk;
		m_pFilterLocii = pLocii;
		}

	if((ChromID = AddChrom(pszChrom)) < 1)
		{
		CCSVFile::Close();
		return(ChromID);
		}

	m_pFilterLocii[m_NumFilterLocii].ChromID = ChromID;
	m_pFilterLocii[m_NumFilterLocii].StartLoci = StartLoci;
	m_pFilterLocii[m_NumFilterLocii++].EndLoci = EndLoci;
	LociLen=1 + EndLoci - StartLoci;
	if(LociLen > m_MaxLociLen)
		m_MaxLociLen = LociLen;
	}

CCSVFile::Close();
if(m_NumFilterLocii > 1)
	{
	qsort(m_pFilterLocii,m_NumFilterLocii,sizeof(tsFilterLoci),SortLocii);
	tsFilterLoci *pTmp1;
	tsFilterLoci *pTmp2;
	int NumDeduped = 1;
	bool bDeduped = false;
	pTmp1 = m_pFilterLocii;
	pTmp2 = &m_pFilterLocii[1];
	for(int Idx = 1; Idx < m_NumFilterLocii; Idx++, pTmp2++)
		{
		if(pTmp1->ChromID == pTmp2->ChromID)
			{
			if(pTmp1->EndLoci >= pTmp2->StartLoci)
				{
				if(pTmp2->EndLoci > pTmp1->EndLoci)
					pTmp1->EndLoci = pTmp2->EndLoci;
				bDeduped = true;
				continue;
				}
			}
		if(bDeduped)
			*pTmp1 = *pTmp2;
		pTmp1 += 1;
		NumDeduped += 1;
		}
	m_NumFilterLocii = NumDeduped;
	}
return(m_NumFilterLocii);
}

UINT16 
CFilterLoci::GenNameHash(char *pszName)
{
unsigned long hash = 5381;
char Chr;
while (Chr = *pszName++)
	hash = ((hash << 5) + hash) + tolower(Chr);
return ((UINT16)hash);
}

int
CFilterLoci::LocateChrom(char *pszChrom)
{
tsLociChrom *pChrom;
int Idx;
UINT16 Hash;
Hash = GenNameHash(pszChrom);

if(m_pLociChroms != NULL && m_NumLociChroms)
	{
	if(m_CachLociChromID > 0)
		{
		pChrom = &m_pLociChroms[m_CachLociChromID-1];
		if((pChrom->Hash == Hash) && !stricmp(pszChrom,pChrom->szChrom))
			return(pChrom->ChromID);
		}

	pChrom = m_pLociChroms;
	for(Idx = 0; Idx < m_NumLociChroms; Idx++,pChrom++)
		{
		if((pChrom->Hash == Hash) && !stricmp(pszChrom,pChrom->szChrom))
			return(m_CachLociChromID = pChrom->ChromID);
		}
	}
return(0);
}

int
CFilterLoci::AddChrom(char *pszChrom)
{
tsLociChrom *pChrom;
int Idx;
UINT16 Hash;
Hash = GenNameHash(pszChrom);

if(m_pLociChroms != NULL && m_NumLociChroms)
	{
	if(m_CachLociChromID > 0)
		{
		pChrom = &m_pLociChroms[m_CachLociChromID-1];
		if((pChrom->Hash == Hash) && !stricmp(pszChrom,pChrom->szChrom))
			return(pChrom->ChromID);
		}

	pChrom = m_pLociChroms;
	for(Idx = 0; Idx < m_NumLociChroms; Idx++,pChrom++)
		{
		if((pChrom->Hash == Hash) && !stricmp(pszChrom,pChrom->szChrom))
			return(m_CachLociChromID = pChrom->ChromID);
		}
	}

if(m_pLociChroms == NULL || m_NumLociChroms == m_AllocdLociChroms)
	{
	if((pChrom = new tsLociChrom[m_AllocdLociChroms + cFiltLociAllocChunk])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to alloc memory for chrom names");
		return(eBSFerrMem);
		}
	if(m_pLociChroms != NULL)
		{
		memmove(pChrom,m_pLociChroms,sizeof(tsLociChrom) * m_NumLociChroms);
		delete m_pLociChroms;
		}
	else
		{
		m_AllocdLociChroms = 0;
		m_NumLociChroms = 0;
		}
	m_pLociChroms = pChrom;
	m_AllocdLociChroms += cFiltLociAllocChunk;
	}
pChrom = &m_pLociChroms[m_NumLociChroms++];
pChrom->ChromID = m_NumLociChroms;
pChrom->Hash = Hash;
strcpy(pChrom->szChrom,pszChrom);
return(m_CachLociChromID = pChrom->ChromID);
}

// returns true if chrom.start-chrom.end overlaps any Loci
bool 
CFilterLoci::Locate(char *pszChrom, int StartLoci, int EndLoci)
{
int ChromID;
tsFilterLoci *pProbe;
int Lo,Mid,Hi;	// search limits
if((ChromID = LocateChrom(pszChrom)) < 1)
	return(false);

Lo = 0; Hi = m_NumFilterLocii-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pProbe = &m_pFilterLocii[Mid];
	if(ChromID < pProbe->ChromID)
		{
		Hi = Mid - 1;
		continue;
		}

	if(ChromID > pProbe->ChromID)
		{
		Lo = Mid + 1;
		continue;
		}

	// on matching chromosome
	// check if have an overlap
	if(EndLoci >= pProbe->StartLoci && StartLoci <= pProbe->EndLoci)
		return(true);

	if(StartLoci < pProbe->StartLoci)
		{
		Hi = Mid - 1;
		continue;
		}

	Lo = Mid + 1;
	continue;
	}
return(false);
}


int 
CFilterLoci::SortLocii( const void *arg1, const void *arg2)
{
tsFilterLoci *pEl1 = (tsFilterLoci *)arg1;
tsFilterLoci *pEl2 = (tsFilterLoci *)arg2;
if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->StartLoci < pEl2->StartLoci)
	return(-1);
if(pEl1->StartLoci > pEl2->StartLoci)
	return(1);
if(pEl1->EndLoci < pEl2->EndLoci)
	return(-1);
if(pEl1->EndLoci > pEl2->EndLoci)
	return(1);
return(0);
}



