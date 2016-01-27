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

CHyperEls::CHyperEls(void)
{
m_pElements = NULL;
m_pChroms = NULL;
m_pChromHashes = NULL;
m_pSeqBases = NULL;
Reset();
}

CHyperEls::~CHyperEls(void)
{
if(m_pSeqBases != NULL)
	{
#ifdef _WIN32
	free(m_pSeqBases);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqBases != MAP_FAILED)
		munmap(m_pSeqBases,m_MemAllocSeqBases);
#endif
	}

if(m_pElements != NULL)
	{
#ifdef _WIN32
	free(m_pElements);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pElements != MAP_FAILED)
		munmap(m_pElements,m_MemAllocEls);
#endif
	}
if(m_pChroms != NULL)
	{
#ifdef _WIN32
	free(m_pChroms);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pChroms != MAP_FAILED)
		munmap(m_pChroms,m_MemAllocChroms);
#endif
	}

if(m_pChromHashes != NULL)
	delete m_pChromHashes;
}


void
CHyperEls::Reset(void)
{
if(m_pSeqBases != NULL)
	{
#ifdef _WIN32
	free(m_pSeqBases);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqBases != MAP_FAILED)
		munmap(m_pSeqBases,m_MemAllocSeqBases);
#endif
	m_pSeqBases = NULL;
	}

if(m_pElements != NULL)
	{
#ifdef _WIN32
	free(m_pElements);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pElements != MAP_FAILED)
		munmap(m_pElements,m_MemAllocEls);
#endif
	m_pElements = NULL;
	}
m_MemAllocSeqBases = 0;
m_MemUsedSeqBases = 0;
m_NumSeqs = 0;
m_MemAllocEls = 0;
m_NumElsAllocd = 0;
m_NumEls = 0;
m_NumElsPlus = 0;
m_NumElsMinus = 0;
m_MaxLenEl = 0;
m_MinLenEl = -1;


if(m_pChroms != NULL)
	{
#ifdef _WIN32
	free(m_pChroms);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pChroms != MAP_FAILED)
		munmap(m_pChroms,m_MemAllocChroms);
#endif
	m_pChroms = 0;
	}

if(m_pChromHashes != NULL)
	{
	delete m_pChromHashes;
	m_pChromHashes = NULL;
	}
m_NumChromsAllocd = 0;
m_MemAllocChroms = 0;
m_NumChroms = 0;
m_MRAChromID = 0;
m_MRAElTypeID = 0;		// most recently accessed element type identifier
m_NumElTypes = 0;
m_MRARefSpeciesID=0;		// most recently accessed ref species identifier
m_NumRefSpecies=0;
m_MRARelSpeciesID=0;		// most recently accessed rel species identifier
m_NumRelSpecies=0;

m_NumElsParsed=0;		// number of elements parsed before filtering
m_NumFiltLen=0;			// number elements filtered out because of length
m_FiltDeduped=0;		// number of elements filter out because of deduping

m_EstNumEls = 0;
}

int
CHyperEls::AddChrom(char *pszChrom)	// get unique chromosome identifier
{
tsHyperChrom *pNewChrom;
UINT16 aHash = CUtility::GenHash16(pszChrom);

if(m_pChroms == NULL)
	{
	m_MemAllocChroms = cChromsInitalAllocNum * sizeof(tsHyperChrom);
#ifdef _WIN32
	m_pChroms = (tsHyperChrom *) malloc(m_MemAllocChroms);	// initial and perhaps the only allocation

	if(m_pChroms == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddChrom: Memory allocation of %lld bytes - %s",(INT64)m_MemAllocChroms,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pChroms = (tsHyperChrom *)mmap(NULL,m_MemAllocChroms, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pChroms == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddChrom: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)m_MemAllocChroms,strerror(errno));
		m_pChroms = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_NumChromsAllocd = cChromsInitalAllocNum;
	m_NumChroms = 0;
	}

if(m_pChromHashes == NULL)
	{
	m_pChromHashes = (UINT64 *) new UINT64 [0x010001];
	memset(m_pChromHashes,0,sizeof(UINT64)*0x010001);
	}

if((m_NumChroms+5) >= m_NumChromsAllocd)
	{
	size_t memreq = m_MemAllocChroms + (cChromsGrowAllocNum * sizeof(tsHyperChrom));
#ifdef _WIN32
	pNewChrom = (tsHyperChrom *) realloc(m_pChroms,memreq);
	if(pNewChrom == NULL)
		{
#else
	pNewChrom = (tsHyperChrom *)mremap(m_pChroms,m_MemAllocChroms,memreq,MREMAP_MAYMOVE);
	if(pNewChrom == MAP_FAILED)
		{
#endif
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddElCore: Memory reallocation to %lld bytes failed - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pChroms = pNewChrom;
	m_MemAllocChroms = memreq;
	m_NumChromsAllocd += cChromsGrowAllocNum;
	}

if(m_NumChroms != 0)
	{
	// high probability that chromosome will be one which was last accessed
	if(m_MRAChromID > 0)
		{
		pNewChrom = &m_pChroms[m_MRAChromID-1];
		if(aHash == pNewChrom->Hash && !stricmp(pszChrom,pNewChrom->szChromName))
			return(pNewChrom->ChromID);
		}

	if(m_pChromHashes[aHash-1] != 0)		// at least one other chrom hashed to same hash?
		{
		pNewChrom = &m_pChroms[m_pChromHashes[aHash-1]-1];
		do
			{
			if(aHash == pNewChrom->Hash && !stricmp(pszChrom,pNewChrom->szChromName))
				return(pNewChrom->ChromID);
			if(pNewChrom->NxtHashChromIdx == 0)
				break;
			pNewChrom = &m_pChroms[pNewChrom->NxtHashChromIdx - 1];
			}
		while(1);
		pNewChrom->NxtHashChromIdx = m_NumChroms + 1;
		}
	else
		m_pChromHashes[aHash-1] = m_NumChroms + 1;
	}
else
	m_pChromHashes[aHash-1] = 1;

// new chromosome entry
pNewChrom = &m_pChroms[m_NumChroms++];
pNewChrom->Hash = aHash;
pNewChrom->NxtHashChromIdx = 0;
strcpy(pNewChrom->szChromName,pszChrom);
pNewChrom->ChromID = m_NumChroms;
m_MRAChromID = m_NumChroms;
return(m_NumChroms);
}

// AddElType
// returns unique element type identifier
int
CHyperEls::AddElType(char *pszElType)
{
tsItemName *pNewType;
unsigned short aHash = CUtility::GenHash16(pszElType);

if(m_NumElTypes == cMaxElTypes)
	return(-1);

if(m_NumElTypes != 0)
	{
	// high probability that type will be one which was last accessed
	if(m_MRAElTypeID > 0)
		{
		pNewType = &m_ElTypes[m_MRAElTypeID-1];
		if(aHash == pNewType->Hash && !stricmp(pszElType,(char *)pNewType->szName))
			return(pNewType->NameID);
		}
	// not most recently accessed, need to do a linear search
	pNewType = &m_ElTypes[0];
	for(m_MRAElTypeID = 1; m_MRAElTypeID <= m_NumElTypes; m_MRAElTypeID++,pNewType++)
		if(aHash == pNewType->Hash && !stricmp(pszElType,(char *)pNewType->szName))
			return(pNewType->NameID);
	}
// new element type entry
pNewType = &m_ElTypes[m_NumElTypes++];
pNewType->Hash = aHash;
strcpy((char *)pNewType->szName,pszElType);
pNewType->NameID = m_NumElTypes;
m_MRAElTypeID = m_NumElTypes;
return(m_NumElTypes);
}

// AddRefSpecies
// returns unique ref species identifier
int
CHyperEls::AddRefSpecies(char *pszSpecies)
{
tsItemName *pNewSpecies;
if(pszSpecies == NULL || pszSpecies[0] == '\0')
	return(-1);

unsigned short aHash = CUtility::GenHash16(pszSpecies);

if(m_NumRefSpecies == cMaxSpecies)
	return(-1);

if(m_NumRefSpecies != 0)
	{
	// high probability that species will be one which was last accessed
	if(m_MRARefSpeciesID > 0)
		{
		pNewSpecies = &m_RefSpecies[m_MRARefSpeciesID-1];
		if(aHash == pNewSpecies->Hash && !stricmp(pszSpecies,(char *)pNewSpecies->szName))
			return(pNewSpecies->NameID);
		}
	// not most recently accessed, need to do a linear search
	pNewSpecies = &m_RefSpecies[0];
	for(m_MRARefSpeciesID = 1; m_MRARefSpeciesID <= m_NumElTypes; m_MRARefSpeciesID++,pNewSpecies++)
		if(aHash == pNewSpecies->Hash && !stricmp(pszSpecies,(char *)pNewSpecies->szName))
			return(pNewSpecies->NameID);
	}
// new element type entry
pNewSpecies = &m_RefSpecies[m_NumRefSpecies++];
pNewSpecies->Hash = aHash;
strcpy((char *)pNewSpecies->szName,pszSpecies);
pNewSpecies->NameID = m_NumRefSpecies;
m_MRARefSpeciesID = m_NumRefSpecies;
return(m_NumRefSpecies);
}

// AddRelSpecies
// returns unique rel species identifier
int
CHyperEls::AddRelSpecies(char *pszSpecies)
{
tsItemName *pNewSpecies;
unsigned short aHash = CUtility::GenHash16(pszSpecies);

if(m_NumRelSpecies == cMaxSpecies)
	return(-1);

if(m_NumRelSpecies != 0)
	{
	// high probability that species will be one which was last accessed
	if(m_MRARelSpeciesID > 0)
		{
		pNewSpecies = &m_RelSpecies[m_MRARelSpeciesID-1];
		if(aHash == pNewSpecies->Hash && !stricmp(pszSpecies,(char *)pNewSpecies->szName))
			return(pNewSpecies->NameID);
		}
	// not most recently accessed, need to do a linear search
	pNewSpecies = &m_RelSpecies[0];
	for(m_MRARelSpeciesID = 1; m_MRARelSpeciesID <= m_NumElTypes; m_MRARelSpeciesID++,pNewSpecies++)
		if(aHash == pNewSpecies->Hash && !stricmp(pszSpecies,(char *)pNewSpecies->szName))
			return(pNewSpecies->NameID);
	}
// new element type entry
pNewSpecies = &m_RelSpecies[m_NumRelSpecies++];
pNewSpecies->Hash = aHash;
strcpy((char *)pNewSpecies->szName,pszSpecies);
pNewSpecies->NameID = m_NumRelSpecies;
m_MRARelSpeciesID = m_NumRelSpecies;
return(m_NumRelSpecies);
}

int 
CHyperEls::PreAllocMem(int EstNumEls, int MeanSeqLen) // preallocate memory for this estimate of number elements having this mean sequence length
{
size_t memreqseqs;
size_t memreqels;
m_EstNumEls = ((size_t)EstNumEls * 105)/100;						// allowing additional to reduce chances of requiring a subsequent realloc
memreqseqs = (size_t)m_EstNumEls * sizeof(etSeqBase) * MeanSeqLen; 
memreqels = (size_t)m_EstNumEls * sizeof(tsHyperElement);

if(m_pSeqBases != NULL && m_MemAllocSeqBases >= memreqseqs && 		// why alloc if already allocated?
   m_pElements != NULL && m_MemAllocEls >= memreqels)
	{
	m_MemUsedSeqBases = 0;
	m_NumSeqs = 0;
	m_NumEls = 0;
	m_NumElsPlus = 0;
	m_NumElsMinus = 0;
	m_MaxLenEl = 0;
	m_MinLenEl = -1;
	return(eBSFSuccess);
	}

if(m_pSeqBases != NULL && m_MemAllocSeqBases < memreqseqs)
	{
#ifdef _WIN32
	free(m_pSeqBases);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqBases != MAP_FAILED)
		munmap(m_pSeqBases,m_MemAllocSeqBases);
#endif
	m_pSeqBases = NULL;	
	m_MemAllocSeqBases = 0;
	}

if(m_pElements != NULL &&  m_MemAllocEls < memreqels)
	{
#ifdef _WIN32
	free(m_pElements);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pElements != MAP_FAILED)
		munmap(m_pElements,m_MemAllocEls);
#endif
	m_pElements = NULL;
	m_MemAllocEls = 0;
	}


if(m_pSeqBases == NULL)
	{
	m_MemAllocSeqBases = memreqseqs; 

#ifdef _WIN32
	m_pSeqBases = (etSeqBase *) malloc(m_MemAllocSeqBases);	// initial and perhaps the only allocation
	if(m_pSeqBases == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"PreAllocMem: Memory allocation of %lld bytes - %s",(INT64)m_MemAllocSeqBases,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
		// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pSeqBases = (etSeqBase *)mmap(NULL,m_MemAllocSeqBases, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSeqBases == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"PreAllocMem: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)m_MemAllocSeqBases,strerror(errno));
		m_pSeqBases = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	}
memset(m_pSeqBases,0,m_MemAllocSeqBases);
m_MemUsedSeqBases = 0;
m_NumSeqs = 0;

if(m_pElements == NULL)
	{
	m_NumElsAllocd = m_EstNumEls;
	m_MemAllocEls = memreqels;
#ifdef _WIN32
	m_pElements = (tsHyperElement *) malloc(m_MemAllocEls);	// initial and perhaps the only allocation
	if(m_pElements == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddElCore: Memory allocation of %lld bytes - %s",(INT64)m_MemAllocEls,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
		// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pElements = (tsHyperElement *)mmap(NULL,m_MemAllocEls, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pElements == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddElCore: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)m_MemAllocEls,strerror(errno));
		m_pElements = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	}
memset(m_pElements,0,m_MemAllocEls);	
m_NumEls = 0;
m_NumElsPlus = 0;
m_NumElsMinus = 0;
m_MaxLenEl = 0;
m_MinLenEl = -1;
return(eBSFSuccess);
}

int
CHyperEls::AddElCore(int SrcID,	// element identifier as parsed from source file
		  char *pszElType,	// element type
		  char *pszRefSpecies, // reference species
		  char *pszRelSpecies, // relative species list
		  char *pszChrom,	// element is on this chromosome
		  int StartLoci,	// element starts at this loci
		  int Len,			// and is of this length
		  int Features,		// element overlays these features
		  char Strand,		// on this strand
		  int RelScale,		// relative scaling factor for this read
		  int NumSNPBases,  // up to a max of 4 variant bases can be specified when SNP calling
		  etSeqBase *pSeq) // any associated sequence, NULL if none,  sequence is expected to be of length NumSNPBases (if not 0) otherwise Len
{
tsHyperElement *pEl;
UINT8 *pBase;
UINT8 Base;
UINT8 SNPBases;
size_t memreq;
int Idx;

if(pszRelSpecies == NULL || pszRelSpecies[0] == '\0' ||
	pszRefSpecies == NULL || pszRefSpecies[0] == '\0' ||
	pszChrom == NULL || pszChrom[0] == '\0')
	return(eBSFerrInternal);

if(NumSNPBases < 0 || NumSNPBases > 4 || (NumSNPBases > 0 && pSeq == NULL))
	return(eBSFerrInternal);

SNPBases = 0;
if(NumSNPBases)
	{
	for(Idx = 0; Idx < NumSNPBases; Idx++,pSeq += 1)
		{
		Base = *pSeq & 0x03;
		if(Base > eBaseT)
			return(eBSFErrBase);
		SNPBases |= Base << (Idx * 2);
		}
	}


if(!NumSNPBases && pSeq != NULL && m_pSeqBases == NULL)
	{
	if(m_EstNumEls > 0)
		memreq = (size_t)(((size_t)m_EstNumEls * 110)/(size_t)100) * sizeof(etSeqBase) * 100;
	else
		memreq = (size_t)cElInitalAllocNum  * sizeof(etSeqBase) * (size_t)100;
	m_MemAllocSeqBases = memreq;
#ifdef _WIN32
	m_pSeqBases = (etSeqBase *) malloc(m_MemAllocSeqBases);	// initial and perhaps the only allocation

	if(m_pSeqBases == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddElCore: Memory allocation of %lld bytes - %s",(INT64)m_MemAllocSeqBases,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pSeqBases = (etSeqBase *)mmap(NULL,m_MemAllocSeqBases, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSeqBases == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddElCore: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)m_MemAllocSeqBases,strerror(errno));
		m_pSeqBases = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_MemUsedSeqBases = 0;
	m_NumSeqs = 0;
	}

if(!NumSNPBases && pSeq != NULL && (m_MemUsedSeqBases + Len + 10000)  >= m_MemAllocSeqBases)
	{
	memreq = m_MemAllocSeqBases + (cElGrowAllocNum * sizeof(etSeqBase) * 200);
#ifdef _WIN32
	pBase = (etSeqBase *) realloc(m_pSeqBases,memreq);
	if(pBase == NULL)
		{
#else
	pBase = (etSeqBase *)mremap(m_pSeqBases,m_MemAllocSeqBases,memreq,MREMAP_MAYMOVE);
	if(pBase == MAP_FAILED)
		{
#endif
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddElCore: Memory reallocation to %lld bytes failed - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pSeqBases = pBase;
	m_MemAllocSeqBases = memreq;
	}

if(m_pElements == NULL)
	{
	if(m_EstNumEls > 0)
		{
		m_NumElsAllocd = (int)(((size_t)m_EstNumEls * 110)/(size_t)100);
		memreq = (size_t)m_NumElsAllocd * sizeof(tsHyperElement);
		}
	else
		{
		memreq = (size_t)cElInitalAllocNum * sizeof(tsHyperElement);
		m_NumElsAllocd = cElInitalAllocNum;
		}
	m_MemAllocEls = memreq;
#ifdef _WIN32
	m_pElements = (tsHyperElement *) malloc(m_MemAllocEls);	// initial and perhaps the only allocation

	if(m_pElements == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddElCore: Memory allocation of %lld bytes - %s",(INT64)m_MemAllocEls,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pElements = (tsHyperElement *)mmap(NULL,m_MemAllocEls, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pElements == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddElCore: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)m_MemAllocEls,strerror(errno));
		m_pElements = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	
	m_NumEls = 0;
	m_NumElsPlus = 0;
	m_NumElsMinus = 0;
	m_MaxLenEl = 0;
	m_MinLenEl = -1;
	}

if((m_NumEls + 1000)  >= m_NumElsAllocd)
	{
	memreq = m_MemAllocEls + (cElGrowAllocNum * sizeof(tsHyperElement));
#ifdef _WIN32
	pEl = (tsHyperElement *) realloc(m_pElements,memreq);
	if(pEl == NULL)
		{
#else
	pEl = (tsHyperElement *)mremap(m_pElements,m_MemAllocEls,memreq,MREMAP_MAYMOVE);
	if(pEl == MAP_FAILED)
		{
#endif
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddElCore: Memory reallocation to %lld bytes failed - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pElements = pEl;
	m_MemAllocEls = memreq;
	m_NumElsAllocd += cElGrowAllocNum;
	}

pEl = &m_pElements[m_NumEls];
m_NumEls += 1;
pEl->SrcID = SrcID;
pEl->ElTypeID = AddElType(pszElType);
pEl->RefSpeciesID = AddRefSpecies(pszRefSpecies);
pEl->RelSpeciesID = AddRelSpecies(pszRelSpecies);
pEl->ElID = m_NumEls;
pEl->ChromID = AddChrom(pszChrom);	// get chromosome identifier
pEl->StartLoci = StartLoci;
pEl->Len = Len;
pEl->NumSNPBases = NumSNPBases;
pEl->SNPBases = SNPBases;
pEl->Features = Features;
pEl->RelScale = RelScale;
pEl->PlusStrand = Strand == '-' ? 0 :  1;
if(pEl->PlusStrand)
	m_NumElsPlus += 1;
else
	m_NumElsMinus += 1;
if(Len > m_MaxLenEl)
	m_MaxLenEl = Len;
if(m_MinLenEl == -1 || m_MinLenEl > Len)
	m_MinLenEl = Len;

if(!NumSNPBases && pSeq != NULL && (*pSeq & 0x07) < eBaseEOS)
	{
	pEl->SeqBasesOfs = m_MemUsedSeqBases + 1;
	pBase = &m_pSeqBases[m_MemUsedSeqBases++];
	
	while(Len--)
		{
		Base = *pSeq++ & 0x07;
		if(Len > 0)
			Base |= (*pSeq++ & 0x07) << 4;
		*pBase++ = Base;
		m_MemUsedSeqBases += 1;
		}
	m_NumSeqs += 1;
	}
else
	pEl->SeqBasesOfs = 0;

return(m_NumEls);
}

int				// returns number of elements parsed
CHyperEls::ParseCSVFileElements(char *pszFile,	  // file (in CSV format) containing element loci
							 int MinLength,  // slough any elements which are less MinLength length
							 int MaxLength,  // slough any elements which are more than MaxLength length
							 teCSVFormat CSVFormat) // expected CSV format
{
int NumFields;
int NumEls = 0;
int Rslt;
int SrcID;
char *pszChrom;
char *pszElType;
char *pszRefSpecies;
char *pszRelSpecies;
char *pszStrand;
int StartLoci;
int EndLoci;
int Len;
int Features;
char Strand;
int RelScale;

int AutoFormat;

char *pSrc;
char *pDst;
char Chr;

char *pszChkFormat;
char *pszSNPbases;
etSeqBase SNPbases[100];
etSeqBase *pSNPbases;
int NumSNPbases;
int NumCultivars;
int NumSNPCultivars;

char szTargAssemblyName[cMaxDatasetSpeciesChrom+1]; // targeted assembly against which SNPs in a 'biokanga snpmarkers' generated file were called
int NumCultivarNames;
char szCultivarNames[100][cMaxDatasetSpeciesChrom+1];		// name of this cultivar, can be at most 100 cultivars in a 'biokanga snpmarkers' generated file
int CultEls;


CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszFile))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszFile);
	delete pCSV;
	return(Rslt);
	}

m_EstNumEls = pCSV->EstNumRows();

m_NumFiltLen = 0;
m_NumElsParsed = 0;
SrcID = 0;
AutoFormat = 0;
szTargAssemblyName[0] = '\0';
NumCultivarNames = 0;
CultEls = 0;
memset(szCultivarNames,0,sizeof(szCultivarNames));

while((Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = pCSV->GetCurFields();
	if(NumFields < 8)				// all formats are 8+ fields so error if less
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 8+ fields in '%s', GetCurFields() returned '%d at row %d'",pszFile,NumFields, m_NumElsParsed);
		return(eBSFerrFieldCnt);
		}

	if(!m_NumElsParsed && (NumFields > 7 && pCSV->IsLikelyHeaderLine())) // expecting at least 8 fields so only then check for header line
		{
		if(NumFields < 12)
			continue;
				// at least 12 fields then parse header assuming it is a 'biokanga snpmarkers' format
		if((((NumFields - 4) % 8) != 0 && (((NumFields - 4) % 9) != 0)))
			continue;
		pCSV->GetText(2,&pszChkFormat);
		if(stricmp(pszChkFormat,"Loci"))
			continue;
		pCSV->GetText(3,&pszChkFormat);
		if(stricmp(pszChkFormat,"TargBase"))
			continue;
		pCSV->GetText(4,&pszChkFormat);
		if(stricmp(pszChkFormat,"NumSpeciesWithCnts"))
			continue;
		if(((NumFields - 4) % 9) == 0)
			{
			CultEls = 9;
			AutoFormat = 3;
			NumCultivars = (NumFields - 4) / 9;
			}
		else
			{
			CultEls = 8;
			AutoFormat = 4;
			NumCultivars = (NumFields - 4) / 8;
			}
		pCSV->GetText(1,&pSrc);
		pDst = szTargAssemblyName;
		Len = 0;
		while(Len < sizeof(szTargAssemblyName)-1 && (Chr = *pSrc++) && Chr != ':')
			{
			*pDst++ = Chr;
			Len++;
			*pDst='\0';
			}

		for(NumCultivarNames = 0; NumCultivarNames < NumCultivars; NumCultivarNames++)
			{
			memset(&szCultivarNames[NumCultivarNames],0,sizeof(szCultivarNames[NumCultivarNames]));
			pCSV->GetText(5+(NumCultivarNames*CultEls),&pSrc);
			pDst = szCultivarNames[NumCultivarNames];
			Len = 0;
			while(Len < sizeof(&szCultivarNames[NumCultivarNames])-1 && (Chr = *pSrc++) && Chr != ':')
				{
				*pDst++ = Chr;
				Len++;
				*pDst='\0';
				}
			}
		continue;
		}

	if(AutoFormat == 0 && CSVFormat == eCSVFdefault)							// default layout and not checked for SNP call formats?
		{
		AutoFormat = 1;			// 1 flags that there has been a check for the CSV format, then will be set to be 2 if simple SNp loci, 3 if snpmarkers post-3.0.6, and 4 if snpmarkers pre-3.0.6
		NumSNPCultivars = 0;
		NumCultivars = 0;

		// how to determine which CSV format to be parsed?
		// if 'biokanga snpmarkers' generated then there will be (4 + ({8 or 9}) * numcultivars) fields, {8} if pre-3.0.6 or {9} if post-3.0.6
        // if a simple SNP loci then there will be 8 if no variant bases, and 9 if variant bases are specified
        // if any other CSV format then there will be 8 or more fields
        // to identify as a 'biokanga snpmarkers' then check if num fields as expected for this format and that field 4 is not text chrom but is integer 1..numcultivars
        // to identify as simple SNP loci then check if 8 or 9 fields, and if 9 then the this field only contains cannonical bases A, C, G, or T.
		// if not identifed as either snpmarkers or simple SNP loci, then process as other - eCSVFdefault, eCSVFprobe or eCSVFAlignM3  

		// try to determine if the CSV file format is either snpmarkers or simple SNPs, if not then will be processing as other
		if(NumFields == 8 || NumFields == 9)	// could be simple SNP loci calls
			{
			pCSV->GetText(2,&pszChkFormat);		// expecting the element type to be "SNP" if simple SNP loci
			if(!stricmp(pszChkFormat,"SNP"))
				AutoFormat = 2;					// accepting format as being simple SNP loci with optional variant bases
			}
		else      // else check if snpmarkers
			{
			if(NumCultivarNames > 0 &&  NumFields >= 12)	
				{
				pCSV->GetInt(4,&NumSNPCultivars);			// expecting the number of cultivars with called SNPs, if was text (chrom name) then not a snpmarker
				if(NumSNPCultivars > 0 && NumSNPCultivars <= 100) // needs to be at least 1 and less than max number of supported marker cultivars
					{
					pCSV->GetText(5,&pszChkFormat);		// if snpmarkers then expect non-integer text, if not snpmarkers then will be an integer
					if(!isdigit(*pszChkFormat))			// if snpmarkers then will not start with a digit
						{
						pCSV->GetText(6,&pszChkFormat);	// if snpmarkers post-3.0.8 then expect non-integer text, if pre-3.0.8 then will be an integer
						if(!isdigit(*pszChkFormat))
							{
							if(((NumFields - 4) % 9) == 0)
								{
								NumCultivars = (NumFields - 4) / 9;
								if(NumSNPCultivars <= NumCultivars)
									AutoFormat = 3;
								}
							}
						}
					else
						{
						if(((NumFields - 4) % 8) == 0)
							{
							NumCultivars = (NumFields - 4) / 8;
							if(NumSNPCultivars <= NumCultivars)
								AutoFormat = 4;
							}
						}
					}
				}
			}
		if(AutoFormat == 1)
			{
			NumCultivars = 0;
			NumSNPCultivars = 0;
			}
		}

	pSNPbases = NULL;
	NumSNPbases = 0;
	m_NumElsParsed += 1;

	switch(CSVFormat) {
		case eCSVFdefault:							// default layout
			if(AutoFormat == 2)			// simple SNP loci with optional variant bases
				{
								// example row: 1,"SNP","AFNW1","Chrom123",1011,1011,1,"+","AG"  : note that the final field containing variant bases is optional
				pCSV->GetText(2,&pszElType);		// expecting the element type to be "SNP" if simple SNP loci
				if(stricmp(pszElType,"SNP"))
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expecting 'SNP' as the element type, parsed '%s' near line %d",pszElType,m_NumElsParsed);
					delete pCSV;
					return(eBSFerrParse);
					}
				pCSV->GetInt(7,&Len);
				if(Len != 1)			// SNPs expected to be single bases
					{
					m_NumFiltLen += 1;
					continue;
					}
				pCSV->GetInt(1,&SrcID);
				pCSV->GetText(3,&pszRefSpecies);
				pCSV->GetText(4,&pszChrom);
				pCSV->GetInt(5,&StartLoci);
				pCSV->GetInt(6,&EndLoci);
				if(StartLoci != EndLoci)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expecting start and end loci to be equal, parsed %d and %d near line %d",StartLoci,EndLoci,m_NumElsParsed);
					delete pCSV;
					return(eBSFerrParse);
					}
				pCSV->GetText(8,&pszStrand);
				if(*pszStrand != '+')
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expecting strand to be '+', parsed '+' near line %d",pszStrand,m_NumElsParsed);
					delete pCSV;
					return(eBSFerrParse);
					}
				Strand = '+';
				pszRelSpecies = pszRefSpecies;
				RelScale = 0;
				Features = 0;
				// SNP marker simple CSV will have 8 or 9 fields, if 9 fields then that is expected to contain list of variant SNP bases
				NumSNPbases = 0;
				SNPbases[0] = eBaseEOS;
				if(NumFields == 9)									// field 9 if present contains list of variant SNP bases
					{
					pCSV->GetText(9,&pszSNPbases);					// assume variant bases
					if(pszSNPbases != NULL)
						{
						// only accept the cannonical variant bases and limit to 4
						int Idx;
						int NumBasesParsed = (int)strlen(pszSNPbases);
						char *pszVarBases = pszSNPbases;
						for(Idx = 0; Idx < NumBasesParsed && NumSNPbases < 4; Idx++,pszVarBases++) 
							{
							switch(*pszVarBases) {
								case 'a': case 'A':
									SNPbases[NumSNPbases++] = eBaseA;
									continue;
								case 'c': case 'C':
									SNPbases[NumSNPbases++] = eBaseC;
									continue;
								case 'g': case 'G':
									SNPbases[NumSNPbases++] = eBaseG;
									continue;
								case 't': case 'T': case 'u': case 'U':
									SNPbases[NumSNPbases++] = eBaseT;
									continue;
								case ' ': case '\t': case ',': case '|':
									continue;
								default:	
									gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expecting variant bases to be one or more ACGT, parsed '%s' near line %d",pszSNPbases,m_NumElsParsed);
									delete pCSV;
									return(eBSFerrParse);
								}
							}
						}
					}
				if(NumSNPbases > 0)
					pSNPbases = SNPbases;
				else
					pSNPbases = NULL;
				break;
				}
			else
				{
				if(AutoFormat == 3 ||				// 'biokanga snpmarkers' generated with 9 fields per cultivar
					AutoFormat == 4)				// 'biokanga snpmarkers' generated with 8 fields per cultivar
					{
					SrcID += 1;
					if(MinLength > 1)
						{
						m_NumFiltLen += 1;
						continue;
						}
					pszElType = (char *)"SNP";
					pszRefSpecies = (char *)"RefSpecies";
					pszRelSpecies = (char *)"RelAllCultivars";
					pCSV->GetText(1,&pszChrom);
					pCSV->GetInt(2,&StartLoci);
					EndLoci = StartLoci;
					RelScale = 0;
					Features = 0;
					Strand = '+';
					Len = 1;
					pSNPbases = NULL; // using ref base as the SNP base
					break;
					}
				}
			// note fallthrough if AutoFormat not for simple SNP loci or for snpmarkers

		case eCSVFprobe:							// use probe loci
		case eCSVFAlignM3:							// use GenMAlignScore initial fields
			pCSV->GetInt(7,&Len);
			if(Len < MinLength || Len > MaxLength)
				{
				m_NumFiltLen += 1;
				continue;
				}

			pCSV->GetInt(1,&SrcID);
			pCSV->GetText(2,&pszElType);
			pCSV->GetText(3,&pszRefSpecies);
			pCSV->GetText(4,&pszChrom);
			pCSV->GetInt(5,&StartLoci);
			pCSV->GetInt(6,&EndLoci);
			pCSV->GetText(8,&pszRelSpecies);


			if(!(m_NumElsParsed % 1000000))
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processed %d elements",m_NumElsParsed);

			if(pszRelSpecies[0] == '+' || pszRelSpecies[0] == '-')	// instead of the relative species, field 8 may contain the strand
				{
				Strand = *pszRelSpecies;
				pszRelSpecies = pszRefSpecies;
				Features = 0;
				RelScale = 0;
				if(NumFields >= 9)					
					pCSV->GetInt(9,&RelScale);
				else
					RelScale = 0;
				}
			else
				{
				Strand = '+';
				if(NumFields > 8 && CSVFormat == eCSVFdefault)
					pCSV->GetInt(9,&Features);
				else
					Features = 0;
				}
			break;
	
		case eCSVFtarget:				// use target loci
			if(NumFields < 12)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 12+ fields in '%s', GetCurFields() returned '%d'",pszFile,NumFields);
				return(eBSFerrFieldCnt);
				}
			m_NumElsParsed += 1;
			pCSV->GetInt(12,&Len);
			if(Len < MinLength || Len > MaxLength)
				{
				m_NumFiltLen += 1;
				continue;
				}
			pCSV->GetInt(1,&SrcID);
			pCSV->GetText(2,&pszElType);
			pCSV->GetText(3,&pszRelSpecies);
			pCSV->GetText(8,&pszRefSpecies);
			pCSV->GetText(9,&pszChrom);
			pCSV->GetInt(10,&StartLoci);
			pCSV->GetInt(11,&EndLoci);
			RelScale = 0;
			Features = 0;
			Strand = '+';
			break;
		}

	Rslt = AddElCore(SrcID,pszElType,pszRefSpecies,pszRelSpecies,pszChrom,StartLoci,Len,Features,Strand,RelScale,NumSNPbases,pSNPbases);
	if(Rslt < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddElCore failed");
		break;
		}
	NumEls++;
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d elements",NumEls);
return(Rslt >= 0 ? NumEls : Rslt);
}

// ParseBEDFileElements
// If de Novo then expecting user to have aligned against a psuedo genome in which contigs were concatenated with spacer N's such that there is a single chromosome
int				// returns number of elements parsed
CHyperEls::ParseBEDFileElements(char *pszFile,	  // file (in CSV format) containing element loci
							 int MinLength,  // slough any elements which are less MinLength length
							 int MaxLength)  // slough any elements which are more than MaxLength length
{
teBSFrsltCodes Rslt;
CBEDfile *pBedFile;

int StartLoci;
int EndLoci;
int Score;
char szTitle[128];
char szChrom[128];
char szFeatName[128];
char FeatStrand;
int CurFeatureID;
int NumEls;
int NumFeatures;

if((pBedFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load from %s",pszFile);
if((Rslt=pBedFile->Open(pszFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(pBedFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pBedFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszFile);
	delete pBedFile;
	Reset();
	return(eBSFerrOpnFile);
	}
NumFeatures = pBedFile->GetNumFeatures();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed load from '%s' of %d features",pszFile,NumFeatures);
if(NumFeatures < 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Nothing to do, no features to map read loci against!");
	delete pBedFile;
	return(eBSFerrNoFeatures);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed load from %s",pszFile);

NumEls = 0;
CurFeatureID = 0;
pBedFile->GetTitle(sizeof(szTitle),szTitle);
while(Rslt >= eBSFSuccess && (CurFeatureID = pBedFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	if(!(NumEls % 5000000))
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d elements",NumEls);
	pBedFile->GetFeature(CurFeatureID,	// feature instance identifier
				szFeatName,				// where to return feature name
				szChrom,				// where to return chromosome name
				&StartLoci,				// where to return feature start on chromosome (0..n)
				&EndLoci,				// where to return feature end on chromosome
 				&Score,					// where to return score
 				&FeatStrand);			// where to return strand

	Rslt = (teBSFrsltCodes)AddElCore(CurFeatureID,(char *)"feat",szTitle,szTitle,szChrom,StartLoci,1 + EndLoci - StartLoci,0,FeatStrand,Score);
	if(Rslt < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddElCore failed");
		break;
		}
	NumEls++;
	}
delete pBedFile;
m_NumElsParsed = NumEls;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d elements",NumEls);
return(Rslt >= 0 ? NumEls : Rslt);
}

char *
CHyperEls::TrimWhitespace(char *pTxt)
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

// SAM format processing
// If de Novo then expecting user to have aligned against a psuedo genome in which contigs were concatenated with spacer N's such that there is a single chromosome
int				// returns number of elements parsed
CHyperEls::ParseSAMFileElements(char *pszFile,	  // file (in SAM format) containing element loci
							 int MinLength,  // slough any elements which are less MinLength length
							 int MaxLength,  // slough any elements which are more than MaxLength length
							 bool bRetainSeq) // true if sequences are to be retained
{
teBSFrsltCodes Rslt;
int NumParsedElLines;
int NumAcceptedEls;
int LineLen;
char szLine[cMaxReadLen  * 3];				// buffer input lines
char szSeqBases[cMaxReadLen+10];	// and may need to retain the alignment sequences
etSeqBase SeqBases[cMaxReadLen+10];

int Tmp = 0;
int SeqLenz;							// sequence length as parsed from SAM/BAM file
int SeqOfs;							// bases starting at this offset are being associated

char *pTxt;
char szDescriptor[128];			// parsed out descriptor
int Flags;						// parsed out flags
char szChrom[128];				// parsed out chrom
int StartLoci;					// start loci
char szCigar[128];
char szRNext[128];
int MAPQ;
int PNext;
int TLen;
int NumUnmappedEls;

CSAMfile BAMfile;

// open SAM for reading
if(pszFile == NULL || *pszFile == '\0')
	return(eBSFerrParams);

if((Rslt = (teBSFrsltCodes)BAMfile.Open(pszFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"ParseSAMFileElements: Unable to open SAM format file %s",pszFile);
	return((teBSFrsltCodes)Rslt);
	}

NumParsedElLines = 0;
NumAcceptedEls = 0;
NumUnmappedEls = 0;

while((LineLen = BAMfile.GetNxtSAMline(szLine)) > 0)
	{
	NumParsedElLines += 1;
	if(!(NumParsedElLines % 5000000) || NumParsedElLines == 1)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d element lines",NumParsedElLines);

	szLine[sizeof(szLine)-1] = '\0';
	pTxt = TrimWhitespace(szLine);
	if(*pTxt=='\0' || *pTxt=='@')	// simply slough lines which were just whitespace or start with '@'
		continue;
	
	// expecting to parse as "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t", szDescriptor, Flags, m_szSAMTargChromName, StartLoci+1,MAPQ,szCigar,pszRNext,PNext,TLen);
	// interest is in the chromname, startloci, length, and optionally the sequence

	if(bRetainSeq)
		{
		sscanf(szLine,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%16383s\t",szDescriptor, &Flags, szChrom, &StartLoci,&MAPQ,szCigar,szRNext,&PNext,&TLen,szSeqBases);
		SeqLenz = 0;
		SeqOfs = 0;
		SeqLenz = (int)strlen(szSeqBases);
		CSeqTrans::MapAscii2Sense(szSeqBases,0,SeqBases);
		SeqBases[SeqLenz] = eBaseEOS;
		}
	else
		{
		sscanf(szLine,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t",szDescriptor, &Flags, szChrom, &StartLoci,&MAPQ,szCigar,szRNext,&PNext,&TLen);
		SeqLenz = 0;
		szSeqBases[0] = '\0';
		SeqBases[0] = eBaseEOS;
		}

		// check if element has been mapped, if not then slough ...
	if(StartLoci == 0 || Flags & 0x04 || szCigar[0] == '*')	// set if unmapped or Cigar is unknown
		{
		NumUnmappedEls += 1;
	    continue;
		}

	// need to estimate  alignment length  from the Cigar string ...
	StartLoci -= 1;

	char *pCgr;
	char Cgr;
	int SegLen;
	int SeqLen;

	SeqOfs = 0;
	SegLen = 0;
	SeqLen = 0;
	pCgr = szCigar;
	Rslt = (teBSFrsltCodes)1;
	while(Rslt >= 1 && (Cgr = *pCgr++) != '\0')
		{
		if(Cgr >= '0' && Cgr <= '9')
			{
			SegLen *= 10;
			SegLen += Cgr - '0';
			continue;
			}
		switch(Cgr) {
			case 'm': case 'M': case '=': case 'x': case 'X':	// match or mismatches
				SeqLen += SegLen;
				continue;

			case 'i': case 'I':				// insertion to the reference
			case 'n': case 'N':				// mRNA to DNA alignment intron, or skipped target region?
				if(SeqLen >=  MinLength && SeqLen <=  MaxLength)
					{
					NumAcceptedEls += 1;
					Rslt = (teBSFrsltCodes)AddElCore(NumAcceptedEls,(char *)"feat",(char *)"NoSpecies",(char *)"NoSpecies",szChrom,StartLoci,SeqLen,0,Flags & 0x010 ? '-' : '+',0,0,SeqBases[0] == eBaseEOS ? NULL : &SeqBases[SeqOfs]);
					if(Rslt < 1)
						{
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"ParseSAMFileElements:AddElCore failed");
						BAMfile.Close();
						return(Rslt);
						}
					}
				else
					NumUnmappedEls += 1;
				SeqOfs += (SeqLen + SegLen);
				StartLoci += (SeqLen + SegLen);
				SegLen = 0;
				SeqLen = 0;
				continue;

			default:
				SegLen = 0;
				SeqLen = 0;
				continue;
			}
		}
	if(Rslt >= 1)
		{
		if(SeqLen >=  MinLength && SeqLen <=  MaxLength)
			{
			NumAcceptedEls += 1;
			Rslt = (teBSFrsltCodes)AddElCore(NumAcceptedEls,(char *)"feat",(char *)"NoSpecies",(char *)"NoSpecies",szChrom,StartLoci,SeqLen,0,Flags & 0x010 ? '-' : '+',0,0,SeqBases[0] == eBaseEOS ? NULL : &SeqBases[SeqOfs]);
			}
		else
			{
			NumUnmappedEls += 1;
			continue;
			}
		}

	if(Rslt < 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ParseSAMFileElements:AddElCore failed");
		break;
		}
	}

m_NumElsParsed = NumAcceptedEls;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d element lines, unmapped %d, accepted %d, ",NumParsedElLines,NumUnmappedEls,NumAcceptedEls);
BAMfile.Close();
return(Rslt >= 0 ? NumAcceptedEls : Rslt);
}

// returns number of elements parsed and accepted
int
CHyperEls::NumEls(void)
{
return(m_NumEls);
}

// returns number of elements filtered out because of length
int
CHyperEls::NumElsFiltLen(void)
{
return(m_NumFiltLen);
}

// returns number of elements parsed before filtering
int
CHyperEls::NumElsParsed(void)
{
return(m_NumElsParsed);
}


// returns number of elements filtered out because of deduping
int
CHyperEls::NumElsDeduped(void)
{
return(m_FiltDeduped);
}

// returns ptr to element identified by ElID
tsHyperElement *
CHyperEls::GetElement(int ElID)
{
if(ElID < 1 || ElID > m_NumEls)
	return(NULL);
return(&m_pElements[ElID-1]);
}

char *
CHyperEls::GetChrom(int ChromID)	// gets ptr to chromosome identified by ChromID
{
if(ChromID < 1 || ChromID > m_NumChroms)
	return(NULL);
return(&m_pChroms[ChromID-1].szChromName[0]);
}


char *CHyperEls::GetType(unsigned char ElTypeID)			// get element type
{
if(ElTypeID < 1 || ElTypeID > m_NumElTypes)
	return(NULL);
return((char *)m_ElTypes[ElTypeID-1].szName);
}

char *
CHyperEls::GetRelSpecies(unsigned char SpeciesID)		// get relative species
{
if(SpeciesID < 1 || SpeciesID > m_NumRelSpecies)
	return(NULL);
return((char *)m_RelSpecies[SpeciesID-1].szName);
}


char *CHyperEls::GetRefSpecies(unsigned char SpeciesID)		// get reference species
{
if(SpeciesID < 1 || SpeciesID > m_NumRefSpecies)
	return(NULL);
return((char *)m_RefSpecies[SpeciesID-1].szName);
}

int
CHyperEls::GetChromID(char *pszChrom) // get chromosome identifier for pszChrom
{
tsHyperChrom *pNewChrom;
if(!m_NumChroms)
	return(eBSFerrChrom);
unsigned short aHash = CUtility::GenHash16(pszChrom);
	// high probability that chromosome will be one which was last accessed
if(m_MRAChromID > 0)
	{
	pNewChrom = &m_pChroms[m_MRAChromID-1];
	if(aHash == pNewChrom->Hash && !stricmp(pszChrom,pNewChrom->szChromName))
		return(pNewChrom->ChromID);
	}
	// not most recently accessed, need to do a linear search
pNewChrom = m_pChroms;
for(m_MRAChromID = 1; m_MRAChromID <= m_NumChroms; m_MRAChromID++,pNewChrom++)
	if(aHash == pNewChrom->Hash && !stricmp(pszChrom,pNewChrom->szChromName))
		return(pNewChrom->ChromID);
m_MRAChromID = 0;
return(eBSFerrChrom);
}

// Use to iterate through elements on specified chrom, returning ptr to element instance
tsHyperElement *
CHyperEls::LocateChromNthElement(int ChromID,		// iterate elements on this chrom
					int nTh)						// return this element instance (1..N)
{
tsHyperElement *pTarg;
tsHyperElement *pMarkTarg;
int MarkIdx;
int Hi;
int Lo;
int Mid;

if(nTh < 1 || nTh > m_NumEls)
	return(NULL);
if(ChromID < 1 || ChromID > m_NumChroms)
	return(NULL);

pMarkTarg = NULL;
Hi = m_NumEls - 1;
Lo = 0;
while(Hi >= Lo)
	{
	Mid = (Hi + Lo)/2;
	pTarg = &m_pElements[Mid];
	if(pTarg->ChromID < ChromID)
		{
		Lo = Mid + 1;
		continue;
		}
	else
		if(pTarg->ChromID > ChromID)
			{
			Hi = Mid - 1;
			continue;
			}
	// same chromosome, getting close...
	if(pMarkTarg == NULL)
		{
		pMarkTarg = pTarg;
		MarkIdx = Mid;
		}
	else
		if(pTarg->StartLoci < pMarkTarg->StartLoci)
			{
			pMarkTarg = pTarg;
			MarkIdx = Mid;
			}
	Hi = Mid - 1;
	}
if(pMarkTarg != NULL)
	{
	while(nTh > 1 && pMarkTarg->ChromID == ChromID)
		{
		nTh -= 1;
		pMarkTarg += 1;
		MarkIdx += 1;
		if(MarkIdx == m_NumEls)
			return(NULL);
		}
	if(pMarkTarg->ChromID != ChromID)
		return(NULL);
	}
return(pMarkTarg);				
}

// Use to iterate through elements on specified chrom, returning the element loci and length
int
CHyperEls::LocateChromNthElement(int ChromID,		// iterate elements on this chrom
					int nTh,						// return this element instance (1..N)
					int *pElLoci,					// element is at this loci on the chrom
					int *pElLen)					// and is this length
{
tsHyperElement *pTarg;

if(pElLoci != NULL)
	*pElLoci = 0;
if(pElLen != NULL)
	*pElLen = 0;
if((pTarg = LocateChromNthElement(ChromID,nTh))==NULL)
	return(0);
if(pElLoci != NULL)
	*pElLoci = pTarg->StartLoci;
if(pElLen != NULL)
	*pElLen = pTarg->Len;
return(pTarg->ElID);    
}


// Locate
// Returns the nTh instance of element which overlaps by at least
// MinOverlap bases on ChromID.StartLoci to ChromID.StartLoci+Len-1
// nTh must be 1..n
int
CHyperEls::Locate(int ChromID,int StartLoci,int Len, int MinOverlap,int nTh,int *pOverlapStartLoci,int *pOverlapLen)
{
tsHyperElement *pTarg;
int Hi;
int Lo;
int Mid;
int ProbeEndLoci;
int	RelEndLoci;
int	OverlapStartLoci;
int	OverlapEndLoci;
int	OverlapLen;

if(nTh < 1 || StartLoci < 0 || MinOverlap < 1 || MinOverlap > m_MaxLenEl)
	return(eBSFerrParams);
if(ChromID < 1 || ChromID > m_NumChroms)
	return(eBSFerrChrom);

// first find an element which has a StartLoci which is immediately
// after the probe ChromID.StartLoci+Len, this ensures that there is no overlap
ProbeEndLoci = StartLoci + Len - 1;
Hi = m_NumEls - 1;
Lo = 0;
while(Hi >= Lo)
	{
	Mid = (Hi + Lo)/2;
	pTarg = &m_pElements[Mid];
	if(pTarg->ChromID < ChromID)
		{
		Lo = Mid + 1;
		continue;
		}
	else
		if(pTarg->ChromID > ChromID)
			{
			Hi = Mid - 1;
			continue;
			}
	// same chromosome, getting close...
	// now looking for targ which starts immediately after probe
	if(pTarg->StartLoci <= ProbeEndLoci)
		{
		Lo = Mid + 1;
		continue;
		}
	else
		if(pTarg->StartLoci > ProbeEndLoci)
			{
			Hi = Mid - 1;
			continue;
			}
	}

// iterate up through m_pElements until certain that pTarg->ChromID > ChromID or pTarg->StartLoci > ProbeEndLoci
if(Lo >= m_NumEls)
	Mid = m_NumEls-1;
else
	if(Hi < 0)
		Mid = 0;
pTarg = &m_pElements[Mid];
while((Mid < m_NumEls-1) && pTarg->ChromID <= ChromID && pTarg->StartLoci < ProbeEndLoci)
	{
	pTarg++;
	Mid++;
	}
pTarg = &m_pElements[Mid];
// now assured that pTarg->StartLoci > ProbeEndLoci
while(Mid-- >= 0)
	{
		// check if backed up to lower ordered chromosome or loci that could never overlap
	if(pTarg->ChromID < ChromID || (pTarg->ChromID == ChromID && (pTarg->StartLoci + m_MaxLenEl) < StartLoci))
		break;

	RelEndLoci = pTarg->StartLoci +  pTarg->Len - 1;
	if(pTarg->ChromID == ChromID &&
		pTarg->StartLoci <= ProbeEndLoci &&
		RelEndLoci >= StartLoci &&
		pTarg->Len >= MinOverlap)
		{
		OverlapStartLoci = pTarg->StartLoci < StartLoci ? StartLoci : pTarg->StartLoci;
		OverlapEndLoci   = RelEndLoci > ProbeEndLoci ? ProbeEndLoci : RelEndLoci;
		OverlapLen = 1 + OverlapEndLoci - OverlapStartLoci;
		if(OverlapLen >= MinOverlap)
			{
			if(! --nTh)
				{
				if(pOverlapStartLoci != NULL)
					*pOverlapStartLoci = OverlapStartLoci;
				if(pOverlapLen != NULL)
					*pOverlapLen = OverlapLen;
				return(pTarg->ElID);
				}
			}
		}
	pTarg--;
	}
if(pOverlapStartLoci != NULL)
	*pOverlapStartLoci = 0;
if(pOverlapLen != NULL)
	*pOverlapLen = 0;
return(0);
}

// LocateNumOverlapping
// Returns the number of elements which overlaps by at least
// MinOverlap bases on ChromID.StartLoci to ChromID.StartLoci+Len-1
int
CHyperEls::LocateNumOverlapping(int ChromID,int StartLoci,int Len, int MinOverlap)
{
tsHyperElement *pTarg;
int NumOverlapping;
int Hi;
int Lo;
int Mid;
int ProbeEndLoci;
int	RelEndLoci;
int	OverlapStartLoci;
int	OverlapEndLoci;
int	OverlapLen;

if(StartLoci < 0 || MinOverlap < 1 || MinOverlap > m_MaxLenEl)
	return(eBSFerrParams);
if(ChromID < 1 || ChromID > m_NumChroms)
	return(eBSFerrChrom);

// first find an element which has a StartLoci which is immediately
// after the probe ChromID.StartLoci+Len, this ensures that there is no overlap
ProbeEndLoci = StartLoci + Len - 1;
Hi = m_NumEls - 1;
Lo = 0;
while(Hi >= Lo)
	{
	Mid = (Hi + Lo)/2;
	pTarg = &m_pElements[Mid];
	if(pTarg->ChromID < ChromID)
		{
		Lo = Mid + 1;
		continue;
		}
	else
		if(pTarg->ChromID > ChromID)
			{
			Hi = Mid - 1;
			continue;
			}
	// same chromosome, getting close...
	// now looking for targ which starts immediately after probe
	if(pTarg->StartLoci <= ProbeEndLoci)
		{
		Lo = Mid + 1;
		continue;
		}
	else
		if(pTarg->StartLoci > ProbeEndLoci)
			{
			Hi = Mid - 1;
			continue;
			}
	}

// iterate up through m_pElements until certain that pTarg->ChromID > ChromID or pTarg->StartLoci > ProbeEndLoci
if(Lo >= m_NumEls)
	Mid = m_NumEls-1;
else
	if(Hi < 0)
		Mid = 0;
pTarg = &m_pElements[Mid];
while((Mid < m_NumEls-1) && pTarg->ChromID <= ChromID && pTarg->StartLoci < ProbeEndLoci)
	{
	pTarg++;
	Mid++;
	}
pTarg = &m_pElements[Mid];
// now assured that pTarg->StartLoci > ProbeEndLoci
NumOverlapping = 0;
while(Mid-- >= 0)
	{
		// check if backed up to lower ordered chromosome or loci that could never overlap
	if(pTarg->ChromID < ChromID || (pTarg->ChromID == ChromID && (pTarg->StartLoci + m_MaxLenEl) < StartLoci))
		break;

	RelEndLoci = pTarg->StartLoci +  pTarg->Len - 1;
	if(pTarg->ChromID == ChromID &&
		pTarg->StartLoci <= ProbeEndLoci &&
		RelEndLoci >= StartLoci &&
		pTarg->Len >= MinOverlap)
		{
		OverlapStartLoci = pTarg->StartLoci < StartLoci ? StartLoci : pTarg->StartLoci;
		OverlapEndLoci   = RelEndLoci > ProbeEndLoci ? ProbeEndLoci : RelEndLoci;
		OverlapLen = 1 + OverlapEndLoci - OverlapStartLoci;
		if(OverlapLen >= MinOverlap)
			NumOverlapping += 1;
		}
	pTarg--;
	}
return(NumOverlapping);
}

int
CHyperEls::LocateLociBaseCnts(int ChromID,int StartLoci,UINT32 *pCntA, UINT32 *pCntC,UINT32 *pCntG,UINT32 *pCntT,UINT32 *pCntN) // returns counts of all bases aligning to specified chrom+loci
{
tsHyperElement *pTarg;
etSeqBase SeqBase;
int NumOverlapping;
int Hi;
int Lo;
int Mid;
int ProbeEndLoci;
int	RelEndLoci;
int	OverlapStartLoci;
int	OverlapEndLoci;
int	OverlapLen;

if(StartLoci < 0 || pCntA == NULL ||pCntC == NULL || pCntG == NULL || pCntT == NULL || pCntN == NULL)
	return(eBSFerrParams);
*pCntA = 0;
*pCntC = 0;
*pCntG = 0;
*pCntT = 0;
*pCntN = 0;

if(ChromID < 1 || ChromID > m_NumChroms)
	return(eBSFerrChrom);

// first find an element which has a StartLoci which is immediately
// after the probe ChromID.StartLoci+Len, this ensures that there is no overlap
ProbeEndLoci = StartLoci;
Hi = m_NumEls - 1;
Lo = 0;
while(Hi >= Lo)
	{
	Mid = (Hi + Lo)/2;
	pTarg = &m_pElements[Mid];
	if(pTarg->ChromID < ChromID)
		{
		Lo = Mid + 1;
		continue;
		}
	else
		if(pTarg->ChromID > ChromID)
			{
			Hi = Mid - 1;
			continue;
			}
	// same chromosome, getting close...
	// now looking for targ which starts immediately after probe
	if(pTarg->StartLoci <= ProbeEndLoci)
		{
		Lo = Mid + 1;
		continue;
		}
	else
		if(pTarg->StartLoci > ProbeEndLoci)
			{
			Hi = Mid - 1;
			continue;
			}
	}

// iterate up through m_pElements until certain that pTarg->ChromID > ChromID or pTarg->StartLoci > ProbeEndLoci
if(Lo >= m_NumEls)
	Mid = m_NumEls-1;
else
	if(Hi < 0)
		Mid = 0;
pTarg = &m_pElements[Mid];
while((Mid < m_NumEls-1) && pTarg->ChromID <= ChromID && pTarg->StartLoci < ProbeEndLoci)
	{
	pTarg++;
	Mid++;
	}
pTarg = &m_pElements[Mid];
// now assured that pTarg->StartLoci > ProbeEndLoci
NumOverlapping = 0;
while(Mid-- >= 0)
	{
		// check if backed up to lower ordered chromosome or loci that could never overlap
	if(pTarg->ChromID < ChromID || (pTarg->ChromID == ChromID && (pTarg->StartLoci + m_MaxLenEl) < StartLoci))
		break;

	RelEndLoci = pTarg->StartLoci +  pTarg->Len - 1;
	if(pTarg->ChromID == ChromID &&
		pTarg->StartLoci <= ProbeEndLoci &&
		RelEndLoci >= StartLoci &&
		pTarg->Len >= 1)
		{
		OverlapStartLoci = pTarg->StartLoci < StartLoci ? StartLoci : pTarg->StartLoci;
		OverlapEndLoci   = RelEndLoci > ProbeEndLoci ? ProbeEndLoci : RelEndLoci;
		OverlapLen = 1 + OverlapEndLoci - OverlapStartLoci;
		if(OverlapLen >= 1)
			{
			NumOverlapping += 1;
			if(pTarg->SeqBasesOfs == 0)
				continue;
			SeqBase  = m_pSeqBases[pTarg->SeqBasesOfs-1 + (OverlapStartLoci - pTarg->StartLoci)/2];
			if((OverlapStartLoci - pTarg->StartLoci) & 0x01)
				SeqBase >>= 4;
			
			switch(SeqBase & 0x07) {
				case eBaseA:
					*pCntA += 1;
					break;
				case eBaseC:
					*pCntC += 1;
					break;
				case eBaseG:
					*pCntG += 1;
					break;
				case eBaseT:
					*pCntT += 1;
					break;
				default:
					*pCntN += 1;
					break;
				}
			}
		}
	pTarg--;
	}
return(NumOverlapping);
}


int
CHyperEls::FiltPutCoresByStrand(int MaxJoinOverlap,		// join or merge two or more core elements if end of one is within this window of start of another
								bool bPlus,				// false if minus, true if plus strand
								bool bDedupe)			// true if cores are to be deduped
{
tsHyperElement *pCoreA;
tsHyperElement *pCoreB;
int AEnd;
int BEnd;
int ABStartDif;

int ElIdx;
int Num2Delete;

int NumEls;
tsHyperElement *pStartEl;
int NumElsDeleted;
int StartMinusIdx;

if((bPlus && m_NumElsPlus < 2) ||
   (!bPlus && m_NumElsMinus < 2))
   return(0);

if(bPlus)
	{
	pStartEl = m_pElements;
	NumEls = m_NumElsPlus;
	StartMinusIdx = m_NumElsPlus;
	}
else
	{
	pStartEl = &m_pElements[m_NumElsPlus];
	NumEls = m_NumElsMinus;
	}

NumElsDeleted = 0;
do {
		// sort cores on chrom.start.len ascending
	qsort(pStartEl,NumEls,sizeof(tsHyperElement),SortElements);
	Num2Delete = 0;
	pCoreA = pStartEl;
	pCoreB = pCoreA + 1;
	for(ElIdx = 0; ElIdx < NumEls-1; ElIdx++,pCoreB++)
		{
		if(pCoreB->Len == 0)		// skip if already deleted
				continue;
		if(pCoreA->ChromID != pCoreB->ChromID) // onto a different chromosome?
			{
			pCoreA = pCoreB;
			continue;
			}
		ABStartDif = pCoreB->StartLoci - pCoreA->StartLoci;
		AEnd = pCoreA->StartLoci + pCoreA->Len - 1;
		BEnd = pCoreB->StartLoci + pCoreB->Len - 1;

		if(bDedupe)
			{
			if(!ABStartDif)	// if same start loci then AEnd will always be <= BEnd as cores sorted length ascending
				{
				pCoreA->Len = pCoreB->Len;
				pCoreA->Features = pCoreB->Features;
				pCoreB->Len = 0;
				Num2Delete+=1;
				continue;
				}
			else					// pCoreB starts after pCoreA
				if(BEnd <= AEnd)	// is pCoreB completely contained in pCoreA?
					{
					pCoreB->Len = 0;
					Num2Delete+=1;
					continue;
					}
			}

			// should join if core start loci only differ by a few bases
		if(MaxJoinOverlap &&
			ABStartDif <= MaxJoinOverlap)
			{
			pCoreA->Len = 1 + BEnd - pCoreA->StartLoci;
			pCoreB->Len = 0;
			Num2Delete+=1;
			continue;
			}

		do {
			pCoreA += 1;
			}
		while(pCoreA->Len == 0);
		}

	if(Num2Delete)
		{
		pCoreA = pStartEl;
		pCoreB = pCoreA;
		for(ElIdx = 0; ElIdx < NumEls; ElIdx++,pCoreB++)
			{
			if(pCoreB->Len != 0)
				{
				if(pCoreA != pCoreB)
					memmove(pCoreA,pCoreB,sizeof(tsHyperElement));
				pCoreA += 1;
				}
			}
		NumEls -= Num2Delete;
		NumElsDeleted += Num2Delete;
		}
	}
while(Num2Delete && NumEls > 1);

if(NumElsDeleted)
	{
	m_NumEls -= NumElsDeleted;
	if(bPlus)
		{
		m_NumElsPlus -= NumElsDeleted;
		if(m_NumElsMinus)
			memmove(&m_pElements[m_NumElsPlus],&m_pElements[StartMinusIdx],m_NumElsMinus * sizeof(tsHyperElement));
		}
	else
		{
		m_NumElsMinus -= NumElsDeleted;
		}
	}
return(NumElsDeleted);
}

int
CHyperEls::FiltPutCores(int MaxJoinOverlap,bool bDedupe)
{
tsHyperElement *pCoreA;
int NumDeleted;
int MinusEndLoci;
int ElIdx;
if(m_NumEls <= 1)
	return(m_NumEls);

// sort cores by strand.chrom.start.len ascending so can process elements strand specific
qsort(m_pElements,m_NumEls,sizeof(tsHyperElement),SortElementsStrand);
if(m_NumElsPlus)
	NumDeleted = FiltPutCoresByStrand(MaxJoinOverlap,true,bDedupe);

if(m_NumElsMinus)
	{
	MinusEndLoci = 0;
	pCoreA = &m_pElements[m_NumElsPlus];
	for(ElIdx = 0; ElIdx < m_NumElsMinus; ElIdx++,pCoreA++)
		{
		pCoreA->StartLoci += pCoreA->Len;
		if(pCoreA->StartLoci > MinusEndLoci)
			MinusEndLoci = pCoreA->StartLoci;
		}
	pCoreA = &m_pElements[m_NumElsPlus];
	for(ElIdx = 0; ElIdx < m_NumElsMinus; ElIdx++,pCoreA++)
		pCoreA->StartLoci = MinusEndLoci - pCoreA->StartLoci;

	NumDeleted += FiltPutCoresByStrand(MaxJoinOverlap,false);
	pCoreA = &m_pElements[m_NumElsPlus];
	for(ElIdx = 0; ElIdx < m_NumElsMinus; ElIdx++,pCoreA++)
		{
		pCoreA->StartLoci = MinusEndLoci - pCoreA->StartLoci;
		pCoreA->StartLoci -= pCoreA->Len;
		}
	}

qsort(m_pElements,m_NumEls,sizeof(tsHyperElement),SortElements);
return(m_NumEls);
}

// Identify split elements by virtue of having shared SrcID's
// If reads have been split (microInDels or splice junctions for example) then
// these reads need to be identified
int
CHyperEls::IdentifySplitElements(void)
{
int ElIdx;
tsHyperElement *pEl1;
tsHyperElement *pEl2;
int SplitOrder;
int NumSplits;

if(m_NumEls == 1)
	{
	pEl1 = m_pElements;
	pEl1->SplitElement = 0;
	pEl1->SplitLast = 0;
	pEl1->SplitOrder = 0;
	}
if(m_NumEls <= 1)
	return(0);

// sort cores by SrcID.strand.chrom.start.len ascending so can process elements strand specific
qsort(m_pElements,m_NumEls,sizeof(tsHyperElement),SortSplitElements);
pEl1 = m_pElements;
pEl2 = pEl1 + 1;
NumSplits = 0;
SplitOrder = 0;
for(ElIdx = 1; ElIdx < m_NumEls; ElIdx++,pEl2++,pEl1++)
	{
	pEl1->SplitElement = 1;			// assume not a split element
	pEl1->SplitFirst = 0;
	pEl1->SplitLast = 0;
	pEl1->SplitElement = 0;

	// now check if a split
	if(pEl1->SrcID != pEl2->SrcID)	// if not sharing same SrcID then no split
		{
		if(SplitOrder != 0)
			{
			pEl1->SplitElement = 1;
			pEl1->SplitLast = 1;
			pEl1->SplitOrder = SplitOrder;
			}
		SplitOrder = 0;
		continue;
		}

	// split
	pEl1->SplitElement = 1;
	if(SplitOrder == 0)
		{
		NumSplits += 1;
		pEl1->SplitFirst = 1;
		}
	pEl1->SplitOrder = SplitOrder++;
	}

if(SplitOrder != 0)
	{
	pEl1->SplitElement = 1;
	pEl1->SplitLast = 1;
	pEl1->SplitOrder = SplitOrder;
	}
else
	{
	pEl1->SplitElement = 0;
	pEl1->SplitLast = 0;
	pEl1->SplitOrder = 0;
	}
qsort(m_pElements,m_NumEls,sizeof(tsHyperElement),SortElements);
return(NumSplits);
}


// DedupeSort
// Removes any duplicate elements and sorts
// Because the core ElID is used as an index into m_pElements then after dedupe the core ElIDs are
// updated to reflect their changed offsets in m_pElements
int			// returns number of elements
CHyperEls::DedupeSort(int MaxJoinOverlap,bool bDedup)
{
tsHyperElement *pCore;
int ElIdx;
int OrigNumELs = m_NumEls;
if(m_pElements == NULL)
	return(0);
if(m_NumEls >= 2)
	{
	FiltPutCores(MaxJoinOverlap,bDedup);
	pCore = m_pElements;
	for(ElIdx = 1; ElIdx <= m_NumEls; ElIdx++,pCore++)
		pCore->ElID = ElIdx;
	}
m_FiltDeduped = OrigNumELs - m_NumEls;
return(m_NumEls);
}

// SortElementsStrand
// Used to sort elements by feature Strand --> ChromID ---> StartLoci --> Len ascending
// Plus strand before minus strand
int
CHyperEls::SortElementsStrand( const void *arg1, const void *arg2)
{
tsHyperElement *pEl1 = (tsHyperElement *)arg1;
tsHyperElement *pEl2 = (tsHyperElement *)arg2;
if(pEl1->PlusStrand > pEl2->PlusStrand)
	return(-1);
if(pEl1->PlusStrand < pEl2->PlusStrand)
	return(1);
if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->StartLoci < pEl2->StartLoci)
	return(-1);
if(pEl1->StartLoci > pEl2->StartLoci)
	return(1);
if(pEl1->Len < pEl2->Len)
	return(-1);
if(pEl1->Len > pEl2->Len)
	return(1);
return(0);
}

// SortElements
// Used to sort elements by feature ChromID ---> StartLoci --> Len --> Strand ascending
// Plus strand before minus strand
int
CHyperEls::SortElements( const void *arg1, const void *arg2)
{
tsHyperElement *pEl1 = (tsHyperElement *)arg1;
tsHyperElement *pEl2 = (tsHyperElement *)arg2;
if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->StartLoci < pEl2->StartLoci)
	return(-1);
if(pEl1->StartLoci > pEl2->StartLoci)
	return(1);
if(pEl1->Len < pEl2->Len)
	return(-1);
if(pEl1->Len > pEl2->Len)
	return(1);
if(pEl1->PlusStrand > pEl2->PlusStrand)
	return(-1);
if(pEl1->PlusStrand < pEl2->PlusStrand)
	return(1);
return(0);
}

// SortSplitElements
// Used to sort elements by feature SrcID --> Strand   --> ChromID ---> StartLoci --> Len ascending
// Plus strand before minus strand
int
CHyperEls::SortSplitElements( const void *arg1, const void *arg2)
{
tsHyperElement *pEl1 = (tsHyperElement *)arg1;
tsHyperElement *pEl2 = (tsHyperElement *)arg2;
if(pEl1->SrcID < pEl2->SrcID)
	return(-1);
if(pEl1->SrcID > pEl2->SrcID)
	return(1);
if(pEl1->PlusStrand > pEl2->PlusStrand)
	return(-1);
if(pEl1->PlusStrand < pEl2->PlusStrand)
	return(1);
if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->StartLoci < pEl2->StartLoci)
	return(-1);
if(pEl1->StartLoci > pEl2->StartLoci)
	return(1);
if(pEl1->Len < pEl2->Len)
	return(-1);
if(pEl1->Len > pEl2->Len)
	return(1);
return(0);
}
