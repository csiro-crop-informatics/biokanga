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
#include "Markers.h"

static int QSortAlignSeqLociSpecies(const void *p1,const void *p2);

CMarkers::CMarkers(void)
{
m_pAllocSeqNames = NULL;
m_pAllocSeqNameIDsOfs = NULL;
m_pSeqNameHashArray = NULL;
m_pAllocAlignLoci = NULL;
m_pHypers = NULL;
Reset();
}


CMarkers::~CMarkers(void)
{
Reset();
}


// Reset
// Release all allocated resources
void 
CMarkers::Reset(void)	// clears all allocated resources
{
if(m_pAllocSeqNames != NULL)
	{
#ifdef _WIN32
	free(m_pAllocSeqNames);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocSeqNames != MAP_FAILED)
		munmap(m_pAllocSeqNames,m_AllocMemSeqNames);
#endif
	m_pAllocSeqNames = NULL;
	}

if(m_pAllocSeqNameIDsOfs != NULL)
	{
#ifdef _WIN32
	free(m_pAllocSeqNameIDsOfs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocSeqNameIDsOfs != MAP_FAILED)
		munmap(m_pAllocSeqNameIDsOfs,m_AllocMemSeqNameIDsOfs);
#endif
	m_pAllocSeqNameIDsOfs = NULL;
	}

if(m_pAllocAlignLoci != NULL)
	{
#ifdef _WIN32
	free(m_pAllocAlignLoci);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pAllocAlignLoci != MAP_FAILED)
		munmap(m_pAllocAlignLoci,m_AllocMemAlignLoci);
#endif
	m_pAllocAlignLoci = NULL;
	}

if(m_pSeqNameHashArray != NULL)
	{
	delete m_pSeqNameHashArray;
	m_pSeqNameHashArray = NULL;
	}

if(m_pHypers != NULL)
	{
	delete m_pHypers;
	m_pHypers = NULL;
	}
m_NumSpecies = 0; 
m_RefSpeciesID = 0;
m_NumSeqNames = 0;		
m_UsedMemSeqNames = 0;		
m_AllocMemSeqNames = 0;	
m_UsedNameHashArray = 0;	
m_AllocMemSeqNameIDsOfs = 0;
m_UsedAlignLoci = 0;	
m_AllocAlignLoci = 0;	
m_AllocMemAlignLoci = 0;	
memset(m_Species,0,sizeof(m_Species));
m_pCurSpecies = NULL;
m_szCurSeqName[0] = '\0';
m_CurSeqNameID = 0;
m_bSorted = false; 
}

UINT16											// returned species identifier (1..cMaxSpecies)
CMarkers::AddSpecies(char *pszSpecies,bool IsRefSpecies)	// cultivar or species name
{
int SpeciesID;
tsSNPSSpecies *pSpecies;

// linear search should be ok as normally only expect a few species involved in marker processing
if(m_NumSpecies > 0)
	{
	if(m_pCurSpecies != NULL && !stricmp((char *)m_pCurSpecies->szSpecies,pszSpecies))
		{
		if(IsRefSpecies && m_RefSpeciesID != m_pCurSpecies->SpeciesID)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSpecies: Attempted to add multiple reference species (%s)",pszSpecies);
			return(0);
			}
		return(m_pCurSpecies->SpeciesID);
		}
	pSpecies = &m_Species[0];
	for(SpeciesID = 0; SpeciesID < m_NumSpecies; SpeciesID++,pSpecies++)
		{
		if(!stricmp((char *)pSpecies->szSpecies,pszSpecies))
			{
			if(IsRefSpecies && m_RefSpeciesID != pSpecies->SpeciesID)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSpecies: Attempted to add multiple reference species (%s)",pszSpecies);
				return(0);
				}
			m_pCurSpecies = pSpecies;
			return(pSpecies->SpeciesID);
			}
		}
	}

if(m_NumSpecies == cMaxMarkerSpecies)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSpecies: Attempted to add more (%s) than the max (%d) supported species names",pszSpecies,cMaxMarkerSpecies);
	return(0);
	}
pSpecies = &m_Species[m_NumSpecies++];
pSpecies->SpeciesID = m_NumSpecies;
pSpecies->IsRefSpecies = IsRefSpecies ? 1 : 0;
if(IsRefSpecies)
	m_RefSpeciesID = pSpecies->SpeciesID;
strncpy((char *)pSpecies->szSpecies,pszSpecies,sizeof(pSpecies->szSpecies)-1);
pSpecies->szSpecies[sizeof(pSpecies->szSpecies)-1] = '\0';
m_pCurSpecies = pSpecies;
return(pSpecies->SpeciesID);
}

char *
CMarkers::SpeciesIDtoName(UINT16 SpeciesID)
{
if(SpeciesID < 1 || SpeciesID > m_NumSpecies)
	return(NULL);
return((char *)&m_Species[SpeciesID-1].szSpecies[0]);
}

bool
CMarkers::IsRefSpecies(UINT16 SpeciesID)
{
if(SpeciesID < 1 || SpeciesID > m_NumSpecies || m_RefSpeciesID < 1)
	return(false);
return(m_RefSpeciesID == SpeciesID);
}

UINT16					// returned species identifier for specified name, returns 0 if name not previously added with AddSpecies)
CMarkers::NameToSpeciesID(char *pszSpecies)
{
int SpeciesID;
tsSNPSSpecies *pSpecies;

if(m_NumSpecies == 0)
	return(0);

// linear search should be ok as normally only expect a few species involved in marker processing
if(m_pCurSpecies != NULL && !stricmp((char *)m_pCurSpecies->szSpecies,pszSpecies))
	return(m_pCurSpecies->SpeciesID);

pSpecies = &m_Species[0];
for(SpeciesID = 0; SpeciesID < m_NumSpecies; SpeciesID++,pSpecies++)
	{
	if(!stricmp((char *)pSpecies->szSpecies,pszSpecies))
		{
		m_pCurSpecies = pSpecies;
		return(pSpecies->SpeciesID);
		}
	}
return(0);
}

char *								// returned sequence name
CMarkers::SeqIDtoName(UINT32 SeqID)	// sequence identifier (1..m_NumSeqNames) for which name is to be returned
{
tsSeqName *pSeqName;
UINT8 *pSeqNames;
UINT64 Ofs;
if(SeqID < 1 || SeqID > m_NumSeqNames)
	return(NULL);
Ofs = m_pAllocSeqNameIDsOfs[SeqID-1];
pSeqNames = (UINT8 *)m_pAllocSeqNames;
pSeqName = (tsSeqName *)&pSeqNames[Ofs];
return((char *)pSeqName->szSeqName);
}

UINT32 
CMarkers::NameToSeqID(char *pszSeqName) // returned sequence identifier for specified name, returns 0 if name not previously added with AddTargSeq)
{
UINT16 Hash;
UINT64 SeqNameOfs;
UINT64 NxtSeqOfs;
tsSeqName *pSeqName;

if(m_pSeqNameHashArray == NULL || m_pAllocSeqNames == NULL || m_NumSeqNames == 0)
	return(0);
// may have struck it lucky and sequence name same as previously added ...
if(m_szCurSeqName[0] != '\0' && !stricmp(pszSeqName,(char *)m_szCurSeqName))
	return(m_CurSeqNameID);
// hash sequence name and use as index into hash array
Hash = CUtility::GenHash16(pszSeqName);
SeqNameOfs = m_pSeqNameHashArray[Hash];	// SeqNameOfs will be 0 if this is a new hash not previously processed

UINT8 *pSeqNames = (UINT8 *)m_pAllocSeqNames;	// used as a convenience instead of requiring many casts when subsequently deriving pSeqName

if(SeqNameOfs == 0)	// no sequence name with this hash previously added?
	return(0);

// seen same hash previously, iterate along names with same hash and check to see if name already known
NxtSeqOfs = SeqNameOfs;
while(SeqNameOfs != 0)
	{
	pSeqName = (tsSeqName *)&pSeqNames[SeqNameOfs - 1];
	if(!stricmp((char *)pSeqName->szSeqName,pszSeqName))
		{
		strcpy((char *)m_szCurSeqName,pszSeqName);
		m_CurSeqNameID = pSeqName->SeqID;
		return(m_CurSeqNameID);
		}
	SeqNameOfs = pSeqName->NxtSeqOfs;
	}
return(0);
}

UINT32										// returned sequence identifier (1..cMaxSeqID)
CMarkers::AddTargSeq(char *pszSeqName)	// sequence name - could be a chrom, contig, transcript name
{
int SeqNameLen;
UINT16 Hash;
UINT64 SeqNameOfs;
UINT64 NxtSeqOfs;
tsSeqName *pSeqName;

// allocate for 16bit name hashes
if(m_pSeqNameHashArray == NULL)
	{
	if((m_pSeqNameHashArray = new UINT64 [0x010001])==NULL)	// allowing for case whereby hash was generated exceeding 16bits!!!!!
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Memory allocation of %lld bytes - %s",(INT64)0x010001 * sizeof(UINT64),strerror(errno));
		return(eBSFerrMem);
		}
	m_UsedNameHashArray = 0;
	memset((size_t *)m_pSeqNameHashArray,0,(size_t)(sizeof(UINT64) * 0x010001));
	}

if(m_pAllocSeqNames == NULL)		// initial allocation?
	{
	size_t memreq = cAllocMemSeqNames;

#ifdef _WIN32
	m_pAllocSeqNames = (tsSeqName *) malloc(memreq);	// initial and perhaps the only allocation
	if(m_pAllocSeqNames == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pAllocSeqNames = (tsSeqName *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pAllocSeqNames == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
		m_pAllocSeqNames = NULL;
		return(eBSFerrMem);
		}
#endif
	m_AllocMemSeqNames = (UINT64)cAllocMemSeqNames;
	m_UsedMemSeqNames = 0;
	m_NumSeqNames = 0;
	}
else
	{
	if(m_AllocMemSeqNames <= (cAllocMinDiffSeqNames + m_UsedMemSeqNames))	// play safe and increase allocation?
		{
		size_t memreq;
		memreq = cAllocMemSeqNames + (size_t)m_AllocMemSeqNames;
#ifdef _WIN32
		pSeqName = (tsSeqName *) realloc(m_pAllocSeqNames,memreq);
#else
		pSeqName = (tsSeqName *)mremap(m_pAllocSeqNames,m_AllocMemSeqNames,memreq,MREMAP_MAYMOVE);
		if(pSeqName == MAP_FAILED)
			pSeqName = NULL;
#endif
		if(pSeqName == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Sequence names memory re-allocation to %lld from %lld bytes - %s",(INT64)memreq,m_AllocMemSeqNames,strerror(errno));
			return(eBSFerrMem);
			}
		m_pAllocSeqNames = pSeqName;
		m_AllocMemSeqNames = memreq;
		}
	}


if(m_pAllocSeqNameIDsOfs == NULL)		// initial allocation?
	{
	size_t memreq = cAllocSeqNames * sizeof(UINT64);

#ifdef _WIN32
	m_pAllocSeqNameIDsOfs = (UINT64 *) malloc(memreq);	// initial and perhaps the only allocation
	if(m_pAllocSeqNameIDsOfs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Memory allocation to %lld bytes - %s",(INT64)memreq,strerror(errno));
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pAllocSeqNameIDsOfs = (UINT64 *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pAllocSeqNameIDsOfs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
		m_pAllocSeqNameIDsOfs = NULL;
		return(eBSFerrMem);
		}
#endif
	m_AllocMemSeqNameIDsOfs = (UINT64)memreq;
	}
else
	{
	if(m_AllocMemSeqNameIDsOfs <= ((m_NumSeqNames + 10) * sizeof(UINT64)))	// play safe and increase allocation?
		{
		UINT64 *pRealloc;
		size_t memreq;
		memreq = (size_t)(m_AllocMemSeqNames + (cAllocSeqNames * sizeof(UINT64)));
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AddTargSeq: memory re-allocation to %lld from %lld bytes",(INT64)memreq,m_AllocMemSeqNameIDsOfs);

#ifdef _WIN32
		pRealloc = (UINT64 *) realloc(m_pAllocSeqNameIDsOfs,memreq);
#else
		pRealloc = (UINT64 *)mremap(m_pAllocSeqNameIDsOfs,m_AllocMemSeqNameIDsOfs,memreq,MREMAP_MAYMOVE);
		if(pRealloc == MAP_FAILED)
			pRealloc = NULL;
#endif
		if(pRealloc == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Sequence name idex memory re-allocation to %lld from %lld bytes - %s",(INT64)memreq,m_AllocMemSeqNameIDsOfs,strerror(errno));
			return(eBSFerrMem);
			}
		m_pAllocSeqNameIDsOfs = (UINT64 *)pRealloc;
		m_AllocMemSeqNameIDsOfs = memreq;
		}
	}

// may have struck it lucky and sequence name same as previously added ...
if(m_szCurSeqName[0] != '\0' && !stricmp(pszSeqName,(char *)m_szCurSeqName))
	return(m_CurSeqNameID);

// hash sequence name and use as index into hash array
Hash = CUtility::GenHash16(pszSeqName);
SeqNameOfs = m_pSeqNameHashArray[Hash];	// SeqNameOfs will be 0 if this is a new hash not previously processed

UINT8 *pSeqNames = (UINT8 *)m_pAllocSeqNames;	// used as a convenience instead of requiring many casts when subsequently deriving pSeqName

if(SeqNameOfs == 0)	// no sequence name with this hash previously added?
	{
	m_UsedNameHashArray += 1;
	NxtSeqOfs = 0;
	}
else
	{
	// seen same hash previously, iterate along names with same hash and check to see if name already known
	NxtSeqOfs = SeqNameOfs;
	while(SeqNameOfs != 0)
		{
		pSeqName = (tsSeqName *)&pSeqNames[SeqNameOfs - 1];
		if(!stricmp((char *)pSeqName->szSeqName,pszSeqName))
			{
			strcpy((char *)m_szCurSeqName,pszSeqName);
			m_CurSeqNameID = pSeqName->SeqID;
			return(m_CurSeqNameID);
			}
		SeqNameOfs = pSeqName->NxtSeqOfs;
		}
	}
pSeqName = (tsSeqName *)&pSeqNames[m_UsedMemSeqNames];
SeqNameLen = (int)min(cMaxLenName-1,strlen(pszSeqName));
pSeqName->Len = (UINT8)(sizeof(tsSeqName) + SeqNameLen);
strncpy((char *)pSeqName->szSeqName,pszSeqName,SeqNameLen);
pSeqName->szSeqName[SeqNameLen] = '\0';
pSeqName->SeqID = m_NumSeqNames + 1;
pSeqName->NxtSeqOfs = NxtSeqOfs;
m_pSeqNameHashArray[Hash] = m_UsedMemSeqNames + 1;
m_pAllocSeqNameIDsOfs[m_NumSeqNames] = m_UsedMemSeqNames;
m_UsedMemSeqNames += pSeqName->Len;
m_NumSeqNames += 1;
return(pSeqName->SeqID);
}

int CMarkers::PreAllocSeqs(int EstNumSeqs, int MeanSeqLen) // preallocate memory for this estimate of number sequences having this mean sequence length
{
int Rslt;
Rslt = m_pHypers->PreAllocMem(EstNumSeqs,MeanSeqLen);
return(Rslt);
}

int
CMarkers::PreAllocSNPs(INT64 EstNumSNPS)	// estimating will be required to process this many SNP loci
{
size_t memreq; 
INT64 AllocdAlignLoci;

AllocdAlignLoci = ((99 + EstNumSNPS) * (INT64)105) / (INT64)100;  // allowing for an extra 5% to reduce probability of a realloc being subsequently required if estimate was a little low
memreq = (size_t)AllocdAlignLoci * sizeof(tsAlignLoci); 
#ifdef _WIN32
m_pAllocAlignLoci = (tsAlignLoci *) malloc(memreq);	// initial and perhaps the only allocation
if(m_pAllocAlignLoci == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"PreAllocSNPs: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pAllocAlignLoci = (tsAlignLoci *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pAllocAlignLoci == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"PreAllocSNPs: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pAllocAlignLoci = NULL;
	return(eBSFerrMem);
	}
#endif
m_AllocMemAlignLoci = memreq;
m_UsedAlignLoci = 0;
m_AllocAlignLoci = AllocdAlignLoci;
memset(m_pAllocAlignLoci,0,memreq);
return(eBSFSuccess);
}

INT64 
CMarkers::NumAlignLoci(void)					// returns current number of alignment/SNP loci
{
return(m_UsedAlignLoci);
}

INT64 
CMarkers::AddLoci(UINT16 TargSpeciesID,	// reads were aligned to this cultivar or species
				UINT32 TargSeqID,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				UINT32 TargLoci,		// loci within target sequence at which SNPs observed
				etSeqBase TargRefBase,	// loci has this reference base
				UINT16 ProbeSpeciesID,	// reads were aligned from this cultivar or species
				UINT32 ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				UINT32 ProbeCntC,		// number instances probe base C aligned to TargRefBase
				UINT32 ProbeCntG,		// number instances probe base G aligned to TargRefBase
				UINT32 ProbeCntT,		// number instances probe base T aligned to TargRefBase
				UINT32 ProbeCntN,		// number instances probe base U aligned to TargRefBase
				UINT16 Flags)			// any loci associated flags
{
tsAlignLoci *pLoci;

if(m_pAllocAlignLoci == NULL)		// initial allocation?
	{
	size_t memreq = cAllocAlignLoci * sizeof(tsAlignLoci);

#ifdef _WIN32
	m_pAllocAlignLoci = (tsAlignLoci *) malloc(memreq);	// initial and perhaps the only allocation
	if(m_pAllocAlignLoci == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddLoci: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pAllocAlignLoci = (tsAlignLoci *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pAllocAlignLoci == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddLoci: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
		m_pAllocAlignLoci = NULL;
		return(eBSFerrMem);
		}
#endif
	m_AllocMemAlignLoci = memreq;
	m_UsedAlignLoci = 0;
	m_AllocAlignLoci = cAllocAlignLoci;
	}
else
	{
	if(m_AllocAlignLoci <= (m_UsedAlignLoci + 100))	// play safe and increase allocation?
		{
		size_t memreq;
		INT64 AllocTo;
		AllocTo = (INT64)cReAllocAlignPerc * m_AllocAlignLoci; 
		memreq = (AllocTo * sizeof(tsAlignLoci))/100;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"AddLoci: memory re-allocation to %lld from %lld bytes",(INT64)memreq,(INT64)m_AllocMemAlignLoci);

#ifdef _WIN32
		pLoci = (tsAlignLoci *) realloc(m_pAllocAlignLoci,memreq);
#else
		pLoci = (tsAlignLoci *)mremap(m_pAllocAlignLoci,m_AllocMemAlignLoci,memreq,MREMAP_MAYMOVE);
		if(pLoci == MAP_FAILED)
			pLoci = NULL;
#endif
		if(pLoci == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddLoci: Memory re-allocation to %lld bytes - %s",(INT64)memreq,strerror(errno));
			return(eBSFerrMem);
			}
		m_pAllocAlignLoci = pLoci;
		m_AllocMemAlignLoci = memreq;
		m_AllocAlignLoci = AllocTo; 
		}
	}

pLoci = &m_pAllocAlignLoci[m_UsedAlignLoci++];
pLoci->AlignID = m_UsedAlignLoci;
pLoci->TargSpeciesID = TargSpeciesID;
pLoci->TargSeqID = TargSeqID;
pLoci->TargLoci = TargLoci;
pLoci->TargRefBase = TargRefBase;
pLoci->Flags = Flags;
pLoci->ProbeSpeciesID = ProbeSpeciesID;
pLoci->ProbeBaseCnts[0] = ProbeCntA;
pLoci->ProbeBaseCnts[1] = ProbeCntC;
pLoci->ProbeBaseCnts[2] = ProbeCntG;
pLoci->ProbeBaseCnts[3] = ProbeCntT;
pLoci->ProbeBaseCnts[4] = ProbeCntN;
m_bSorted = false;
return(m_UsedAlignLoci);
}

// AddLoci
// Add loci on target which has identified SNP when aligned to from probe sequences
INT64  
CMarkers::AddLoci(char *pszTargSpecies,	// reads were aligned to this cultivar or species
				char *pszTargSeq,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				UINT32 TargLoci,			// loci within target sequence at which SNPs observed
				etSeqBase TargRefBase,	// loci has this reference base
				char *pszProbeSpecies,	// reads were aligned from this cultivar or species
				UINT32 ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				UINT32 ProbeCntC,		// number instances probe base C aligned to TargRefBase
				UINT32 ProbeCntG,		// number instances probe base G aligned to TargRefBase
				UINT32 ProbeCntT,		// number instances probe base T aligned to TargRefBase
				UINT32 ProbeCntN)		// number instances probe base U aligned to TargRefBase
{
UINT16 TargID;
UINT16 ProbeID;
UINT32 TargSeqID;
INT64 Rslt;

TargID = AddSpecies(pszTargSpecies,true);
if(TargID == 0 || TargID != m_RefSpeciesID)
	return(-1);

ProbeID = AddSpecies(pszProbeSpecies,false);
if(ProbeID == 0 || ProbeID == m_RefSpeciesID)
	return(-1);

TargSeqID = AddTargSeq(pszTargSeq);
if(TargSeqID == 0)
	return(-1);

Rslt = AddLoci(TargID,TargSeqID,TargLoci,TargRefBase,ProbeID,ProbeCntA,ProbeCntC,ProbeCntG,ProbeCntT,ProbeCntN,cFlgSNPcnts);
return(Rslt);
}

int 
CMarkers::LoadSNPFile(int MinBases,			// accept SNPs with at least this number covering bases
					  double MaxPValue,		// accept SNPs with at most this P-value
					  char *pszRefSpecies,	// this is the reference species 
					  char *pszProbeSpecies, // this species reads were aligned to the reference species from which SNPs were called 
					  char *pszSNPFile)		// SNP file to parse and load
{
int Rslt;
INT64 Rslt64;
int NumFields;
int NumElsParsed;

char *pszRefSeq;
etSeqBase RefBase;
int StartLoci;
int Mismatches;
char *pszRefBase;
int MMBaseA;
int MMBaseC;
int MMBaseG;
int MMBaseT;
int MMBaseN;

int CoveringBases;
int NumFilteredOut;
int FilteredCovBases;
int FilteredPValue;

CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszSNPFile))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open SNP file: %s",pszSNPFile);
	delete pCSV;
	return(Rslt);
	}

double PValue;

NumElsParsed = 0;
NumFilteredOut = 0;
FilteredCovBases = 0;
FilteredPValue = 0;

while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(!NumElsParsed && (NumFields >= 23 && pCSV->IsLikelyHeaderLine())) // expecting at least 23 fields so only then check for header line
		continue;

	NumElsParsed += 1;

	// apply any filtering
	pCSV->GetInt(11,&CoveringBases);
	if(CoveringBases < MinBases)
		{
		NumFilteredOut += 1;
		FilteredCovBases += 1;
		continue;
		}

	pCSV->GetDouble(10,&PValue);
	if(PValue > MaxPValue)
		{
		NumFilteredOut += 1;
		FilteredPValue += 1;
		continue;
		}

	pCSV->GetText(4,&pszRefSeq);
	pCSV->GetInt(5,&StartLoci);
	pCSV->GetInt(12,&Mismatches);
	pCSV->GetText(13,&pszRefBase);
	pCSV->GetInt(14,&MMBaseA);
	pCSV->GetInt(15,&MMBaseC);
	pCSV->GetInt(16,&MMBaseG);
	pCSV->GetInt(17,&MMBaseT);
	pCSV->GetInt(18,&MMBaseN);

	switch(pszRefBase[0]) {
		case 'a': case 'A':
			RefBase = eBaseA;
			MMBaseA = CoveringBases - Mismatches;
			break;
		case 'c': case 'C':
			RefBase = eBaseC;
			MMBaseC = CoveringBases - Mismatches;
			break;
		case 'g': case 'G':
			RefBase = eBaseG;
			MMBaseG = CoveringBases - Mismatches;
			break;
		case 't': case 'T':
			RefBase = eBaseT;
			MMBaseT = CoveringBases - Mismatches;
			break;
		default:
			RefBase = eBaseN;
			MMBaseN = CoveringBases - Mismatches;
			break;
		}

	Rslt64 = AddLoci(pszRefSpecies,		// reads were aligned to this cultivar or species
				pszRefSeq,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				StartLoci,		// loci within target sequence at which SNPs observed
				RefBase,		// loci has this reference base
				pszProbeSpecies,// reads were aligned from this cultivar or species
				MMBaseA,		// number instances probe base A aligned to TargRefBase 
				MMBaseC,		// number instances probe base C aligned to TargRefBase
				MMBaseG,		// number instances probe base G aligned to TargRefBase
				MMBaseT,		// number instances probe base T aligned to TargRefBase
				MMBaseN);		// number instances probe base U aligned to TargRefBase

	if(Rslt64 < 1)
		{
		if(pCSV != NULL)
			delete pCSV;
		return(-1);
		}
	}
if(pCSV != NULL)
	delete pCSV;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d SNPs from file: %s",NumElsParsed,pszSNPFile);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted %d SNPs, filtered out %d high P-Value (> %.3f), filtered out %d low coverage ( < %d bases)",
					NumElsParsed - NumFilteredOut, FilteredPValue, MaxPValue, FilteredCovBases, MinBases);

return(NumElsParsed - NumFilteredOut);
}

// AddImputedAlignments
// Add alignments for species where no snp was called but other species do have snp called
// The no call could be because there were none or insufficent reads covering the loci, or there was coverage but no snp!
INT64 
CMarkers::AddImputedAlignments(int MinBases,			// must be at least this number of reads covering the SNP loci
					  char *pszRefSpecies,				// this is the reference species 
					char *pszProbeSpecies,				// this species reads were aligned to the reference species from which SNPs were called 
					char *pszAlignFile,					// file containing alignments
					int FType,							// input alignment file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM)
					bool bSeqs,							// if alignment file contains the read sequence then impute bases from the actual sequences	
					int EstNumSeqs,						// estimated number of sequences (0 if no estimate)
					int EstSeqLen)						// estimated mean sequence length (0 if no estimate)           			

{
int Rslt;
INT64 Rslt64;
int NumEls;
UINT16 RefSpeciesID;
UINT16 ProbeSpeciesID;
int MinLength = 50;
int MaxLength = 1000;
UINT16 ImputFlags;

if((ProbeSpeciesID = NameToSpeciesID(pszProbeSpecies)) < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate identifier for probe species '%s'",pszProbeSpecies);
	return(eBSFerrInternal);
	}
if((RefSpeciesID = NameToSpeciesID(pszRefSpecies)) < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate identifier for probe species '%s'",pszRefSpecies);
	return(eBSFerrInternal);
	}


SortTargSeqLociSpecies();	// must be sorted ....

if(m_pHypers != NULL)
	{
	m_pHypers->Reset();
	delete m_pHypers;
	m_pHypers = NULL;
	}

if((m_pHypers = new CHyperEls)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CHyperEls");
	return(eBSFerrObj);
	}

if(EstNumSeqs != 0)
	if((Rslt = m_pHypers->PreAllocMem(EstNumSeqs,EstSeqLen)) != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to prealloc allocate memory for sequences");
		return(eBSFerrObj);
		}

etClassifyFileType FileType;

if(FType == 0)
	FileType = CUtility::ClassifyFileType(pszAlignFile);
else
	FileType = (etClassifyFileType)(FType - 1);

ImputFlags = cFlgImputCnts;
switch(FileType) {
	case eCFTopenerr:		// unable to open file for reading
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: '%s'",pszAlignFile);
		return(eBSFerrOpnFile);

	case eCFTlenerr:		// file length is insufficent to classify type
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to classify file type (insuffient data points): '%s'",pszAlignFile);
		return(eBSFerrFileAccess);

	case eCFTunknown:		// unable to reliably classify
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to reliably classify file type: '%s'",pszAlignFile);
		return(eBSFerrFileType);

	case eCFTCSV:			// file has been classified as being CSV
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing CSV file (%d..%d): '%s'",MinLength,MaxLength,pszAlignFile);
		if((Rslt = m_pHypers->ParseCSVFileElements(pszAlignFile,MinLength,MaxLength,eCSVFdefault)) < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in CSV file (%d..%d): '%s'",MinLength,MaxLength,pszAlignFile);
			Reset();
			return(Rslt);
			}
		break;

	case eCFTBED:			// file has been classified as being BED
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing BED file (%d..%d): '%s'",MinLength,MaxLength,pszAlignFile);
		if((Rslt = m_pHypers->ParseBEDFileElements(pszAlignFile,MinLength,MaxLength)) < 0)
			{	
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in BED file (%d..%d): '%s'",MinLength,MaxLength,pszAlignFile);
			Reset();
			return(Rslt);
			}
		break;

	case eCFTSAM:			// file has been classified as being SAM
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing SAM file (%d..%d): '%s'",MinLength,MaxLength,pszAlignFile);
		if((Rslt = m_pHypers->ParseSAMFileElements(pszAlignFile,MinLength,MaxLength,bSeqs)) < 0)
			{	
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in SAM file (%d..%d): '%s'",MinLength,MaxLength,pszAlignFile);
			Reset();
			return(Rslt);
			}
		ImputFlags = cFlgAlignCnts;
		break;

	default:
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to classify file type: '%s'",pszAlignFile);
		Reset();
		return(eBSFerrFileType);
	}


NumEls = m_pHypers->NumEls();
if(NumEls == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No elements with length range %d..%d in file: '%s'",MinLength,MaxLength,pszAlignFile);
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loaded and parsed %d elements",NumEls);

//gDiagnostics.DiagOut(eDLInfo,gszProcName,"Now identifying split (microInDels or splice junction spanning?) elements...",NumEls);
//NumSplitEls = m_pHypers->IdentifySplitElements();					// identify any split elements which may be present
//gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identified %d split (microInDels or splice junction spanning?) elements",NumSplitEls);

// try to find overlaps on to SNPs called on other species where this probe species is not represented!
bool bProbeAligned;
INT64 AlignIdx;
INT64 UsedAlignLoci;
tsAlignLoci *pAlign;
UINT32 CurTargSeqID;
UINT32 CurTargLociLoci;
UINT32 PrevTargSeqID;
int NumOverlapping;
int HyperChromID;
char *pszTargSeq;
etSeqBase TargRefBase;

INT64 TotNonOverlapping;
INT64 TotOverlapping;

UINT32 ProbeCntA;		// number instances probe base A aligned to TargRefBase 
UINT32 ProbeCntC;		// number instances probe base C aligned to TargRefBase
UINT32 ProbeCntG;		// number instances probe base G aligned to TargRefBase
UINT32 ProbeCntT;		// number instances probe base T aligned to TargRefBase
UINT32 ProbeCntN;		// number instances probe base U aligned to TargRefBase

TotOverlapping = 0;
TotNonOverlapping = 0;
AlignIdx = 0;
bProbeAligned = false;
PrevTargSeqID = 0;
pszTargSeq = NULL;
UsedAlignLoci = m_UsedAlignLoci;
for(AlignIdx = 0; AlignIdx < UsedAlignLoci; AlignIdx++)
	{
	pAlign = &m_pAllocAlignLoci[AlignIdx];			// m_pAllocAlignLoci could be realloc'd so best to take the address each iteration....
	if(AlignIdx == 0)
		{
		CurTargSeqID = pAlign->TargSeqID;
		CurTargLociLoci = pAlign->TargLoci;
		TargRefBase = pAlign->TargRefBase;
		if(pAlign->ProbeSpeciesID == ProbeSpeciesID)
			bProbeAligned = true;
		else
			bProbeAligned = false;
		continue;
		}


	if(bProbeAligned == false && (CurTargSeqID != pAlign->TargSeqID || CurTargLociLoci != pAlign->TargLoci))
		{
				// no snp called for probe species, check if there were reads aligned to reference
		if(CurTargSeqID != PrevTargSeqID || pszTargSeq == NULL)
			{
			pszTargSeq = SeqIDtoName(CurTargSeqID);
			HyperChromID = m_pHypers->GetChromID(pszTargSeq);
			PrevTargSeqID = CurTargSeqID;
			}

		ProbeCntA = 0;
		ProbeCntC = 0;
		ProbeCntG = 0;
		ProbeCntT = 0;
		ProbeCntN = 0;
		if(HyperChromID < 1)		// will be < 1 if no target sequence alignments in alignment file
			NumOverlapping = 0;
		else
			{
			NumOverlapping = m_pHypers->LocateLociBaseCnts(HyperChromID,CurTargLociLoci,&ProbeCntA,&ProbeCntC,&ProbeCntG,&ProbeCntT,&ProbeCntN);

			if(NumOverlapping == 0 || (ProbeCntA == 0 && ProbeCntC == 0 && ProbeCntG == 0 && ProbeCntT == 0 && ProbeCntN == 0))
				{
				ProbeCntA = 0;
				ProbeCntC = 0;
				ProbeCntG = 0;
				ProbeCntT = 0;
				ProbeCntN = 0;
				switch(TargRefBase) {
					case eBaseA:
						ProbeCntA = NumOverlapping;
						break;
					case eBaseC:
						ProbeCntC = NumOverlapping;
						break;
					case eBaseG:
						ProbeCntG = NumOverlapping;
						break;
					case eBaseT:
						ProbeCntT = NumOverlapping;
						break;
					case eBaseN:
						ProbeCntN = NumOverlapping;
						break;
					}
				}
			}


		if(NumOverlapping == 0)
			TotNonOverlapping += 1;
		else
			TotOverlapping += 1;

		// add number overlapping to a new loci record with refbase count set to NumOverlapping
		Rslt64 = AddLoci(RefSpeciesID,	// reads were aligned to this cultivar or species
				CurTargSeqID,			// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				CurTargLociLoci,		// loci within target sequence at which SNPs observed
				TargRefBase,			// loci has this reference base
				ProbeSpeciesID,			// reads were aligned from this cultivar or species
				ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				ProbeCntC,		// number instances probe base C aligned to TargRefBase
				ProbeCntG,		// number instances probe base G aligned to TargRefBase
				ProbeCntT,		// number instances probe base T aligned to TargRefBase
				ProbeCntN,		// number instances probe base U aligned to TargRefBase
				ImputFlags);			// user flag to indicate these are imputed counts, not from the SNP file
		if(Rslt64 < 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddImputedAlignments: AddLoci returned error %d",(int)Rslt64);
			Reset();
			return(Rslt64);			
			}

		pAlign = &m_pAllocAlignLoci[AlignIdx];
		}

	if(CurTargSeqID != pAlign->TargSeqID || CurTargLociLoci != pAlign->TargLoci)
		{
		CurTargSeqID = pAlign->TargSeqID;
		CurTargLociLoci = pAlign->TargLoci;
		TargRefBase = pAlign->TargRefBase;
		bProbeAligned = false;
		}
	if(pAlign->ProbeSpeciesID == ProbeSpeciesID)
		bProbeAligned = true;
	}

if(!bProbeAligned)
	{
	if(CurTargSeqID != PrevTargSeqID || pszTargSeq == NULL)
		{
		pszTargSeq = SeqIDtoName(CurTargSeqID);
		HyperChromID = m_pHypers->GetChromID(pszTargSeq);
		PrevTargSeqID = CurTargSeqID;
		}
	ProbeCntA = 0;
	ProbeCntC = 0;
	ProbeCntG = 0;
	ProbeCntT = 0;
	ProbeCntN = 0;
	if(HyperChromID < 1)
		NumOverlapping = 0;
	else
		{
		NumOverlapping = m_pHypers->LocateLociBaseCnts(HyperChromID,CurTargLociLoci,&ProbeCntA,&ProbeCntC,&ProbeCntG,&ProbeCntT,&ProbeCntN);

		if(NumOverlapping == 0 || (ProbeCntA == 0 && ProbeCntC == 0 && ProbeCntG == 0 && ProbeCntT == 0 && ProbeCntN == 0))
			{
			ProbeCntA = 0;
			ProbeCntC = 0;
			ProbeCntG = 0;
			ProbeCntT = 0;
			ProbeCntN = 0;
			switch(TargRefBase) {
				case eBaseA:
					ProbeCntA = NumOverlapping;
					break;
				case eBaseC:
					ProbeCntC = NumOverlapping;
					break;
				case eBaseG:
					ProbeCntG = NumOverlapping;
					break;
				case eBaseT:
					ProbeCntT = NumOverlapping;
					break;
				case eBaseN:
					ProbeCntN = NumOverlapping;
					break;
				}
			}
		}
	if(NumOverlapping < MinBases)
		NumOverlapping = 0;
	if(NumOverlapping == 0)
		TotNonOverlapping += 1;
	else
		TotOverlapping += 1;

		// add number overlapping to a new loci record with refbase count set to NumOverlapping
	Rslt64 = AddLoci(RefSpeciesID,	// reads were aligned to this cultivar or species
				CurTargSeqID,			// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				CurTargLociLoci,		// loci within target sequence at which SNPs observed
				TargRefBase,			// loci has this reference base
				ProbeSpeciesID,			// reads were aligned from this cultivar or species
				ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				ProbeCntC,		// number instances probe base C aligned to TargRefBase
				ProbeCntG,		// number instances probe base G aligned to TargRefBase
				ProbeCntT,		// number instances probe base T aligned to TargRefBase
				ProbeCntN,		// number instances probe base U aligned to TargRefBase
				ImputFlags);			// user flag to indicate these are imputed counts, not from the SNP file

	if(Rslt64 < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddImputedAlignments: AddLoci returned error %d",(int)Rslt64);
		Reset();
		return(Rslt64);
		}	
	pAlign = &m_pAllocAlignLoci[AlignIdx];
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Added %d loci with alignments, %d loci with no alignments",TotOverlapping,TotNonOverlapping);
return(TotOverlapping);
}


int
CMarkers::IdentSpeciesSpec(int AltMaxCnt,	// max count allowed for base being processed in any other species, 0 if no limit
						int MinCnt,		// min count required for base being processed in species
						double SNPMmajorPC,		// to be processed major putative SNP base must be at least this proportion of total
						int MinSpeciesWithCnts,			// must be at least this number of species with base counts more than MinSpeciesTotCntThres - 0 if no limit 
						int MinSpeciesTotCntThres)		// individual species must have at least this number of total bases at SNP loci - 0 if no threshold

{
UINT32 AlignIdx;
tsAlignLoci *pAlign;
tsAlignLoci *pAlignSpecies;
tsAlignLoci *pAlignSpeciesA;
int NumSpecies = m_NumSpecies-1;
int SpeciesIdx;
int SpeciesIdxA;
int BaseIdx;
bool bAcceptSpec;
etSeqBase BestAcceptBase;
double BestAcceptConf;
double CurConf;
int NumSpeciesWithCnts;

SortTargSeqLociSpecies();

pAlign = &m_pAllocAlignLoci[0];
for(AlignIdx = 0; AlignIdx < m_UsedAlignLoci; AlignIdx += NumSpecies, pAlign += NumSpecies)
	{
	NumSpeciesWithCnts = 0;
	pAlignSpecies = pAlign;
	for(SpeciesIdx = 0; SpeciesIdx < NumSpecies; SpeciesIdx++,pAlignSpecies += 1)
		{
		pAlignSpecies->TotBases = pAlignSpecies->ProbeBaseCnts[0]+pAlignSpecies->ProbeBaseCnts[1]+pAlignSpecies->ProbeBaseCnts[2]+pAlignSpecies->ProbeBaseCnts[3]+pAlignSpecies->ProbeBaseCnts[4];
		if(pAlignSpecies->TotBases == 0 || (UINT32)MinSpeciesTotCntThres > pAlignSpecies->TotBases)
			{
			pAlignSpecies->CultSpecBase = eBaseN;
			pAlignSpecies->CultSpecBaseConf = 0;
			pAlignSpecies->FiltLowTotBases = 1;
			continue;
			}
		pAlignSpecies->FiltLowTotBases = 0;
		NumSpeciesWithCnts += 1;

		// if proportion of major SNP base above a threshold then check if any of the other species have any bases
		BestAcceptBase = eBaseN;
		BestAcceptConf = 0.0;
		for(BaseIdx = 0; BaseIdx < 4; BaseIdx++)
			{
			CurConf = (pAlignSpecies->ProbeBaseCnts[BaseIdx] / (double)pAlignSpecies->TotBases);
			if(pAlignSpecies->ProbeBaseCnts[BaseIdx] >= (UINT32)MinCnt && (CurConf >= SNPMmajorPC))
				{
				bAcceptSpec = true;
				pAlignSpeciesA = pAlign;
				for(SpeciesIdxA = 0; SpeciesIdxA < NumSpecies; SpeciesIdxA++,pAlignSpeciesA += 1)
					{
					if(SpeciesIdxA == SpeciesIdx)
						continue;

					if(AltMaxCnt > 0 && pAlignSpeciesA->ProbeBaseCnts[BaseIdx] >= (UINT32)AltMaxCnt)
						{
						bAcceptSpec = false;
						break;
						}
					}

				if(bAcceptSpec)
					{
					if(CurConf > BestAcceptConf)
						{
						BestAcceptBase = BaseIdx;
						BestAcceptConf = (pAlignSpecies->ProbeBaseCnts[BaseIdx] / (double)pAlignSpecies->TotBases);
						}
					}
				}

			pAlignSpecies->CultSpecBase = BestAcceptBase;
			pAlignSpecies->CultSpecBaseConf = (UINT8)(100 * BestAcceptConf);
			}
		}
	pAlignSpecies = pAlign;
	for(SpeciesIdx = 0; SpeciesIdx < NumSpecies; SpeciesIdx++,pAlignSpecies += 1)
		pAlignSpecies->NumSpeciesWithCnts = NumSpeciesWithCnts; 
	}
return(0);
}

const int cRptBuffSize = 0x03ffff;          // will allocate this sized reporting buffer
INT64 
CMarkers::Report(char *pszRefGenome,		// reference genome assembly against which other species were aligned
			int NumRelGenomes,				// number of relative genome names
			char *pszRelGenomes[],			// relative genome names
			char *pszReportFile,			// report to this file
			int MinSpeciesWithCnts,			// must be at least this number of species with base counts more than MinSpeciesTotCntThres - 0 if no limit 
			int MinSpeciesTotCntThres,		// individual species must have at least this number of total bases at SNP loci - 0 if no limit
			bool bSloughRefOnly)			// do not report if no inter-cultivar SNP marker, i.e if cultivars all same with the polymorthic site relative to reference only 
{
UINT32 Idx;
int m_hOutFile;
char *pszBuff;
int BuffIdx;
char cBase;
tsAlignLoci *pAlign;
UINT32 CurTargLoci;
UINT32 CurTargSeqID;
char *pszRefSeq;
UINT32 PrevTargSeqID;
bool bFirstReport;
INT64 NumSloughed;
INT64 NumReported;

SortTargSeqLociSpecies();

#ifdef _WIN32
m_hOutFile = open(pszReportFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(pszReportFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszReportFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszReportFile);
	Reset();
	return(eBSFerrCreateFile);
	}

if((pszBuff = new char [cRptBuffSize]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Report: unable to allocate memory for buffering reports");
	Reset();
	return(eBSFerrMem);
	}

BuffIdx = 0;
BuffIdx += sprintf(&pszBuff[BuffIdx],"\"%s:TargSeq\",\"Loci\",\"TargBase\",\"NumSpeciesWithCnts\"",pszRefGenome);
for(Idx = 0; Idx < (UINT32)m_NumSpecies-1; Idx++)
	BuffIdx += sprintf(&pszBuff[BuffIdx],",\"%s:CntsSrc\",\"%s:Base\",\"%s:Score\",\"%s:BaseCntTot\",\"%s:BaseCntA\",\"%s:BaseCntC\",\"%s:BaseCntG\",\"%s:BaseCntT\",\"%s:BaseCntN\"",
				pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx]);

CUtility::SafeWrite(m_hOutFile,pszBuff,BuffIdx);
BuffIdx = 0;

pAlign = &m_pAllocAlignLoci[0];
CurTargLoci = pAlign->TargLoci;
CurTargSeqID = pAlign->TargSeqID;
PrevTargSeqID = 0;
pszRefSeq = NULL;
bFirstReport = true;
NumSloughed = 0;
NumReported = 0;
for(Idx = 0; Idx < m_UsedAlignLoci; Idx++, pAlign++)
	{
	if(MinSpeciesWithCnts > pAlign->NumSpeciesWithCnts)
		continue;

	// user may have requesed that only variants between the cultivars are of interest; if all cultivars have same variant, even if different to reference, then slough
	if(bSloughRefOnly && (bFirstReport || pAlign->TargLoci != CurTargLoci || pAlign->TargSeqID != CurTargSeqID))	// starting new row
		{
		tsAlignLoci *pTmpAlign = pAlign;
		UINT8 CultBase = pAlign->CultSpecBase;
		int NumSpeciesWithCnts = (int)pAlign->NumSpeciesWithCnts;
		int ChkIdx;
		for(ChkIdx = 0; ChkIdx < NumSpeciesWithCnts; ChkIdx++,pTmpAlign++)
			{
			if(CultBase != pTmpAlign->CultSpecBase)
				break;
			}
		if(ChkIdx == NumSpeciesWithCnts)	// all cultivar bases same?
			{
			pAlign = pTmpAlign - 1;
			Idx += ChkIdx - 1;
			NumSloughed += 1;
			continue;
			}
		}

	if(bFirstReport || pAlign->TargLoci != CurTargLoci || pAlign->TargSeqID != CurTargSeqID)		// start new row...
		{
		bFirstReport = false;
		CurTargLoci = pAlign->TargLoci;
		CurTargSeqID = pAlign->TargSeqID;
		// new row start here 
		switch(pAlign->TargRefBase) {
			case eBaseA:
				cBase = 'A';
				break;
			case eBaseC:
				cBase = 'C';
				break;
			case eBaseG:
				cBase = 'G';
				break;
			case eBaseT:
				cBase = 'T';
				break;
			case eBaseN:
				cBase = 'N';
				break;
			}

		if(pszRefSeq == NULL || pAlign->TargSeqID != PrevTargSeqID)
			{
			pszRefSeq = SeqIDtoName(pAlign->TargSeqID);
			PrevTargSeqID = pAlign->TargSeqID;
			}
		BuffIdx += sprintf(&pszBuff[BuffIdx],"\n\"%s\",%d,\"%c\",%d",pszRefSeq,pAlign->TargLoci,cBase,pAlign->NumSpeciesWithCnts);
		if((BuffIdx + 500) > cRptBuffSize)
			{
			CUtility::SafeWrite(m_hOutFile,pszBuff,BuffIdx);
			BuffIdx = 0;
			}
		NumReported += 1;
		}
	switch(pAlign->CultSpecBase) {
			case eBaseA:
				cBase = 'A';
				break;
			case eBaseC:
				cBase = 'C';
				break;
			case eBaseG:
				cBase = 'G';
				break;
			case eBaseT:
				cBase = 'T';
				break;
			case eBaseN:
				cBase = 'N';
				break;
			}
	BuffIdx += sprintf(&pszBuff[BuffIdx],",\"%c\",\"%c\",%d,%d,%d,%d,%d,%d,%d",pAlign->Flags & cFlgSNPcnts ? 'S' : 'I',
						cBase,pAlign->CultSpecBaseConf,pAlign->TotBases,pAlign->ProbeBaseCnts[0],pAlign->ProbeBaseCnts[1],pAlign->ProbeBaseCnts[2],pAlign->ProbeBaseCnts[3],pAlign->ProbeBaseCnts[4]);
	if((BuffIdx + 500) > cRptBuffSize)
		{
		CUtility::SafeWrite(m_hOutFile,pszBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
if(BuffIdx && pszBuff != NULL)
	{
	CUtility::SafeWrite(m_hOutFile,pszBuff,BuffIdx);
	BuffIdx = 0;
	}
#ifdef _WIN32
_commit(m_hOutFile);
#else
fsync(m_hOutFile);
#endif
close(m_hOutFile);
m_hOutFile = -1;
if(pszBuff != NULL)
	delete pszBuff;
return(NumReported);
}

static tsAlignLoci *gpAllocAlignLoci;		// as allocated to hold alignment loci;

INT64
CMarkers::SortTargSeqLociSpecies(void)	
{
// check if anything to sort ....
if(m_pAllocAlignLoci == NULL || m_UsedAlignLoci == 0)
	return(0);

if(m_UsedAlignLoci == 1)
	return(1);
	
if(!m_bSorted)
	{
	gpAllocAlignLoci = m_pAllocAlignLoci;
	qsort(gpAllocAlignLoci,m_UsedAlignLoci,sizeof(tsAlignLoci),QSortAlignSeqLociSpecies);
	}
return(m_UsedAlignLoci);
}

// QSortAlignSeqLociSpecies
// qsorts alignment loci by TargSeqID,TargLoci,ProbeSpeciesID ascending
int CMarkers::QSortAlignSeqLociSpecies(const void *p1,const void *p2)
{
tsAlignLoci *pAlign1;
tsAlignLoci *pAlign2;

pAlign1 = (tsAlignLoci *)p1;
pAlign2 = (tsAlignLoci *)p2;
if(pAlign1->TargSeqID > pAlign2->TargSeqID)
	return(1);
if(pAlign1->TargSeqID < pAlign2->TargSeqID)
	return(-1);
if(pAlign1->TargLoci > pAlign2->TargLoci)
	return(1);
if(pAlign1->TargLoci < pAlign2->TargLoci)
	return(-1);
if(pAlign1->ProbeSpeciesID > pAlign2->ProbeSpeciesID)
	return(1);
if(pAlign1->ProbeSpeciesID < pAlign2->ProbeSpeciesID)
	return(-1);
return(0);
}