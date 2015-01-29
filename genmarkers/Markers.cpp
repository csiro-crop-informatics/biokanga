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
m_MaxQSortThreads = cDfltThreads;
}

UINT8											// returned species identifier (1..cMaxSpecies)
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
CMarkers::SpeciesIDtoName(UINT8 SpeciesID)
{
if(SpeciesID < 1 || SpeciesID > m_NumSpecies)
	return(NULL);
return((char *)&m_Species[SpeciesID-1].szSpecies[0]);
}

bool
CMarkers::IsRefSpecies(UINT8 SpeciesID)
{
if(SpeciesID < 1 || SpeciesID > m_NumSpecies || m_RefSpeciesID < 1)
	return(false);
return(m_RefSpeciesID == SpeciesID);
}

UINT8					// returned species identifier for specified name, returns 0 if name not previously added with AddSpecies)
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
	if((m_pSeqNameHashArray = new UINT64 [0x010000])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Memory allocation of %lld bytes - %s",(INT64)0x010000 * sizeof(UINT64),strerror(errno));
		return(eBSFerrMem);
		}
	m_UsedNameHashArray = 0;
	memset((size_t *)m_pSeqNameHashArray,0,(size_t)(sizeof(UINT64) * 0x010000));
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
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTargSeq: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
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
	if(m_AllocMemSeqNameIDsOfs <= ((m_UsedMemSeqNames + 10) * sizeof(UINT64)))	// play safe and increase allocation?
		{
		UINT64 *pRealloc;
		size_t memreq;
		memreq = (size_t)(m_AllocMemSeqNames + (cAllocSeqNames * sizeof(UINT64)));
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

int 
CMarkers::AddLoci(UINT8 TargSpeciesID,	// reads were aligned to this cultivar or species
				UINT32 TargSeqID,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				UINT32 TargLoci,		// loci within target sequence at which SNPs observed
				etSeqBase TargRefBase,	// loci has this reference base
				UINT8 ProbeSpeciesID,	// reads were aligned from this cultivar or species
				UINT32 ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				UINT32 ProbeCntC,		// number instances probe base C aligned to TargRefBase
				UINT32 ProbeCntG,		// number instances probe base G aligned to TargRefBase
				UINT32 ProbeCntT,		// number instances probe base T aligned to TargRefBase
				UINT32 ProbeCntN)		// number instances probe base U aligned to TargRefBase
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
	if(m_AllocAlignLoci <= (m_UsedAlignLoci + 10))	// play safe and increase allocation?
		{
		size_t memreq;
		memreq = (cAllocAlignLoci + m_AllocAlignLoci) * sizeof(tsAlignLoci);
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
		m_AllocAlignLoci += cAllocAlignLoci; 
		}
	}

pLoci = &m_pAllocAlignLoci[m_UsedAlignLoci++];
pLoci->AlignID = m_UsedAlignLoci;
pLoci->TargSpeciesID = TargSpeciesID;
pLoci->TargSeqID = TargSeqID;
pLoci->TargLoci = TargLoci;
pLoci->TargRefBase = TargRefBase;
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
int 
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
UINT8 TargID;
UINT8 ProbeID;
UINT32 TargSeqID;
int Rslt;

TargID = AddSpecies(pszTargSpecies,true);
if(TargID == 0 || TargID != m_RefSpeciesID)
	return(-1);

ProbeID = AddSpecies(pszProbeSpecies,false);
if(ProbeID == 0 || ProbeID == m_RefSpeciesID)
	return(-1);

TargSeqID = AddTargSeq(pszTargSeq);
if(TargSeqID == 0)
	return(-1);

if((Rslt = AddLoci(TargID,TargSeqID,TargLoci,TargRefBase,ProbeID,ProbeCntA,ProbeCntC,ProbeCntG,ProbeCntT,ProbeCntN)) < 0)
	return(Rslt);

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
	if(!NumElsParsed && (NumFields >= 21 && pCSV->IsLikelyHeaderLine())) // expecting at least 21 fields so only then check for header line
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

	Rslt = AddLoci(pszRefSpecies,		// reads were aligned to this cultivar or species
				pszRefSeq,		// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				StartLoci,		// loci within target sequence at which SNPs observed
				RefBase,		// loci has this reference base
				pszProbeSpecies,// reads were aligned from this cultivar or species
				MMBaseA,		// number instances probe base A aligned to TargRefBase 
				MMBaseC,		// number instances probe base C aligned to TargRefBase
				MMBaseG,		// number instances probe base G aligned to TargRefBase
				MMBaseT,		// number instances probe base T aligned to TargRefBase
				MMBaseN);		// number instances probe base U aligned to TargRefBase

	if(Rslt < 1)
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
int 
CMarkers::AddImputedAlignments(int MinBases,			// must be at least this number of reads covering the SNP loci
					  char *pszRefSpecies,				// this is the reference species 
					char *pszProbeSpecies,				// this species reads were aligned to the reference species from which SNPs were called 
					char *pszAlignFile,					// file containing alignments
					int FType)							// input alignment file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM)
{
int Rslt;
int NumEls;
UINT8 RefSpeciesID;
UINT8 ProbeSpeciesID;
int MinLength = 50;
int MaxLength = 1000;

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
	delete m_pHypers;
	m_pHypers = NULL;
	}

if((m_pHypers = new CHyperEls)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CHyperEls");
	return(eBSFerrObj);
	}

etClassifyFileType FileType;

if(FType == 0)
	FileType = CUtility::ClassifyFileType(pszAlignFile);
else
	FileType = (etClassifyFileType)(FType - 1);

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
		if((Rslt = m_pHypers->ParseSAMFileElements(pszAlignFile,MinLength,MaxLength)) < 0)
			{	
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in SAM file (%d..%d): '%s'",MinLength,MaxLength,pszAlignFile);
			Reset();
			return(Rslt);
			}
		break;
	}


NumEls = m_pHypers->NumEls();
if(NumEls == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No elements with length range %d..%d in CSV file: '%s'",MinLength,MaxLength,pszAlignFile);
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loaded and parsed %d elements",NumEls);

//gDiagnostics.DiagOut(eDLInfo,gszProcName,"Now identifying split (microInDels or splice junction spanning?) elements...",NumEls);
//NumSplitEls = m_pHypers->IdentifySplitElements();					// identify any split elements which may be present
//gDiagnostics.DiagOut(eDLInfo,gszProcName,"Identified %d split (microInDels or splice junction spanning?) elements",NumSplitEls);

// try to find overlaps on to SNPs called on other species where this probe species is not represented!
bool bProbeAligned;
UINT32 AlignIdx;
UINT32 UsedAlignLoci;
tsAlignLoci *pAlign;
UINT32 CurTargSeqID;
UINT32 CurTargLociLoci;
UINT32 PrevTargSeqID;
int NumOverlapping;
int HyperChromID;
char *pszTargSeq;
etSeqBase TargRefBase;

UINT32 TotNonOverlapping;
UINT32 TotOverlapping;

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
		if(HyperChromID < 1)		// will be < 1 if no target sequence alignments in alignment file
			NumOverlapping = 0;
		else
			NumOverlapping = m_pHypers->LocateNumOverlapping(HyperChromID,CurTargLociLoci);

		if(NumOverlapping < MinBases)
			NumOverlapping = 0;
		if(NumOverlapping == 0)
			TotNonOverlapping += 1;
		else
			TotOverlapping += 1;

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

		// add number overlapping to a new loci record with refbase count set to NumOverlapping
		Rslt = AddLoci(RefSpeciesID,	// reads were aligned to this cultivar or species
				CurTargSeqID,			// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				CurTargLociLoci,		// loci within target sequence at which SNPs observed
				TargRefBase,			// loci has this reference base
				ProbeSpeciesID,			// reads were aligned from this cultivar or species
				ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				ProbeCntC,		// number instances probe base C aligned to TargRefBase
				ProbeCntG,		// number instances probe base G aligned to TargRefBase
				ProbeCntT,		// number instances probe base T aligned to TargRefBase
				ProbeCntN);		// number instances probe base U aligned to TargRefBase
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
	if(HyperChromID < 1)
		NumOverlapping = 0;
	else
		NumOverlapping = m_pHypers->LocateNumOverlapping(HyperChromID,CurTargLociLoci);
	if(NumOverlapping < MinBases)
		NumOverlapping = 0;
	if(NumOverlapping == 0)
		TotNonOverlapping += 1;
	else
		TotOverlapping += 1;

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
		// add number overlapping to a new loci record with refbase count set to NumOverlapping
	Rslt = AddLoci(RefSpeciesID,	// reads were aligned to this cultivar or species
				CurTargSeqID,			// alignments to this sequence - could be a chrom/contig/transcript - from pszSpecies
				CurTargLociLoci,		// loci within target sequence at which SNPs observed
				TargRefBase,			// loci has this reference base
				ProbeSpeciesID,			// reads were aligned from this cultivar or species
				ProbeCntA,		// number instances probe base A aligned to TargRefBase 
				ProbeCntC,		// number instances probe base C aligned to TargRefBase
				ProbeCntG,		// number instances probe base G aligned to TargRefBase
				ProbeCntT,		// number instances probe base T aligned to TargRefBase
				ProbeCntN);		// number instances probe base U aligned to TargRefBase
	pAlign = &m_pAllocAlignLoci[AlignIdx];
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Added %d loci with alignments, %d loci with no alignments",TotOverlapping,TotNonOverlapping);
return(0);
}


int
CMarkers::IdentSpeciesSpec(int AltMaxCnt,	// max count allowed for base being processed in any other species
						int MinCnt,		// min count required for base being processed in species
						double Thres,	// to be processed a base must be at least this proportion of total
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

		// if proportion of bases above a threshold then check if any of the other species have any bases
		BestAcceptBase = eBaseN;
		BestAcceptConf = 0.0;
		for(BaseIdx = 0; BaseIdx < 4; BaseIdx++)
			{
			CurConf = (pAlignSpecies->ProbeBaseCnts[BaseIdx] / (double)pAlignSpecies->TotBases);
			if(pAlignSpecies->ProbeBaseCnts[BaseIdx] >= (UINT32)MinCnt && (CurConf >= Thres))
				{
				bAcceptSpec = true;
				pAlignSpeciesA = pAlign;
				for(SpeciesIdxA = 0; SpeciesIdxA < NumSpecies; SpeciesIdxA++,pAlignSpeciesA += 1)
					{
					if(SpeciesIdxA == SpeciesIdx)
						continue;

					if(pAlignSpeciesA->ProbeBaseCnts[BaseIdx] >= (UINT32)AltMaxCnt)
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

int 
CMarkers::Report(char *pszRefGenome,		// reference genome assembly against which other species were aligned
			int NumRelGenomes,				// number of relative genome names
			char *pszRelGenomes[],			// relative genome names
			char *pszReportFile,			// report to this file
			int MinSpeciesWithCnts,			// must be at least this number of species with base counts more than MinSpeciesTotCntThres - 0 if no limit 
			int MinSpeciesTotCntThres)		// individual species must have at least this number of total bases at SNP loci - 0 if no limit
{
UINT32 Idx;
int m_hOutFile;
char szBuff[4096];
int BuffIdx;
char cBase;
tsAlignLoci *pAlign;
UINT32 CurTargLoci;
UINT32 CurTargSeqID;
char *pszRefSeq;
UINT32 PrevTargSeqID;
bool bFirstReport;

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

BuffIdx = 0;
BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s:TargSeq\",\"Loci\",\"TargBase\",\"NumSpeciesWithCnts\"",pszRefGenome);
for(Idx = 0; Idx < (UINT32)m_NumSpecies-1; Idx++)
	BuffIdx += sprintf(&szBuff[BuffIdx],",\"%s:Base\",\"%s:Score\",\"%s:BaseCntTot\",\"%s:BaseCntA\",\"%s:BaseCntC\",\"%s:BaseCntG\",\"%s:BaseCntT\",\"%s:BaseCntN\"",
				pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx],pszRelGenomes[Idx]);

pAlign = &m_pAllocAlignLoci[0];
CurTargLoci = pAlign->TargLoci;
CurTargSeqID = pAlign->TargSeqID;
PrevTargSeqID = 0;
pszRefSeq = NULL;
bFirstReport = true;
for(Idx = 0; Idx < m_UsedAlignLoci; Idx++, pAlign++)
	{
	if(MinSpeciesWithCnts > pAlign->NumSpeciesWithCnts)
		continue;

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
		BuffIdx += sprintf(&szBuff[BuffIdx],"\n\"%s\",%d,\"%c\",%d",pszRefSeq,pAlign->TargLoci,cBase,pAlign->NumSpeciesWithCnts);
		if((BuffIdx + 250) > sizeof(szBuff))
			{
			CUtility::SafeWrite(m_hOutFile,szBuff,BuffIdx);
			BuffIdx = 0;
			}
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
	BuffIdx += sprintf(&szBuff[BuffIdx],",\"%c\",%d,%d,%d,%d,%d,%d,%d",cBase,pAlign->CultSpecBaseConf,pAlign->TotBases,pAlign->ProbeBaseCnts[0],pAlign->ProbeBaseCnts[1],pAlign->ProbeBaseCnts[2],pAlign->ProbeBaseCnts[3],pAlign->ProbeBaseCnts[4]);
	if((BuffIdx + 250) > sizeof(szBuff))
		{
		CUtility::SafeWrite(m_hOutFile,szBuff,BuffIdx);
		BuffIdx = 0;
		}
	}
if(BuffIdx)
	{
	CUtility::SafeWrite(m_hOutFile,szBuff,BuffIdx);
	BuffIdx = 0;
	}
close(m_hOutFile);
return(0);
}

// SetMaxThreads
// Set maximum number of threads to be utilised
int
CMarkers::SetMaxThreads(int MaxThreads)		// maximum number of threads to be utilised
{
if(MaxThreads > 64)
	MaxThreads = 64;
m_MaxQSortThreads = MaxThreads;
m_MTqsort.SetMaxThreads(MaxThreads);
return(MaxThreads);
}

int
CMarkers::SortTargSeqLociSpecies(void)	
{
// check if anything to sort ....
if(m_pAllocAlignLoci == NULL || m_UsedAlignLoci == 0)
	return(0);
if(m_UsedAlignLoci == 1)
	return(1);
if(!m_bSorted)
	{
	m_MTqsort.SetMaxThreads(m_MaxQSortThreads);
	m_MTqsort.qsort(m_pAllocAlignLoci,m_UsedAlignLoci,sizeof(tsAlignLoci),QSortAlignSeqLociSpecies);
	m_bSorted = true;
	}
return(m_UsedAlignLoci);
}

// QSortAlignSeqLociSpecies
// qsorts alignment loci by TargSeqID,TargLoci,ProbeSpeciesID ascending
static int QSortAlignSeqLociSpecies(const void *p1,const void *p2)
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