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

#include "pacbiocommon.h"
#include "SSW.h"
#include "SWAlign.h"


CSWAlign::CSWAlign()
{
m_pSfxArray = NULL;
m_pSWInstances = NULL;
m_bLocksCreated = false;
Reset();
}


CSWAlign::~CSWAlign()
{
Reset();
}

void
CSWAlign::Reset(void)
{
if(m_pSfxArray != NULL)
	{
	delete m_pSfxArray;
	m_pSfxArray = NULL;
	}

if(m_pSWInstances != NULL)
	{
	UINT32 InstIdx;
	tsPBSSWInstance *pSWInstance;
	pSWInstance = m_pSWInstances;
	for(InstIdx = 0; InstIdx < m_NumSWInstances; InstIdx++,pSWInstance++)
		{
		if(pSWInstance->pCoreHits != NULL)
			{
#ifdef _WIN32
			free(pSWInstance->pCoreHits);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
			if(pSWInstance->pCoreHits != MAP_FAILED)
				munmap(pSWInstance->pCoreHits,pSWInstance->AllocdCoreHitsSize);
#endif	
			}

		if(pSWInstance->pProbeSeq != NULL)
			delete pSWInstance->pProbeSeq;

		if(pSWInstance->pTargSeq != NULL)
			delete pSWInstance->pTargSeq;

		if(pSWInstance->pSW != NULL)
			delete pSWInstance->pSW;
		
		if(pSWInstance->pmtqsort != NULL)
			delete pSWInstance->pmtqsort;
		}
	delete m_pSWInstances;
	m_pSWInstances = NULL;
	}

m_NumSWInstances = 0;
m_AllocdSWInstances = 0;
if(m_bLocksCreated)
	{
	DeleteLocks();
	m_bLocksCreated = false;
	}

m_MinOverlapLen=0;		
m_MaxSeedCoreDepth=0;		
m_DeltaCoreOfs=0;			 
m_MinNumCores=0;			
m_MaxAcceptHitsPerSeedCore=0; 

m_DfltMaxProbeSeqLen=0;	

m_NumSWInstances=0;			
m_AllocdSWInstances=0;			
m_MaxTargSeqLen=0;		


m_SWMatchScore = cSSWDfltMatchScore;
m_SWMismatchPenalty = cSSWDfltMismatchPenalty;
m_SWGapOpenPenalty = cSSWDfltGapOpenPenalty;
m_SWGapExtnPenalty = cSSWDfltGapExtnPenalty;
m_SWDlyGapExtn = cSSWDfltDlyGapExtn;
m_SWProgPenaliseGapExtn = cSSWDfltProgPenaliseGapExtn;
m_SWProgExtnPenaltyLen = 0;
m_SWAnchorLen = cSSWDfltAnchorLen;

m_CPMatchScore = cSSWDfltMatchScore;
m_CPMismatchPenalty = cSSWDfltMismatchPenalty;
m_CPGapOpenPenalty = cSSWDfltGapOpenPenalty;
m_CPGapExtnPenalty = cSSWDfltGapOpenPenalty;
m_MaxInitiatePathOfs = cDfltSSWDInitiatePathOfs;

m_OverlapFloat = cDfltSSWDOvlpFloat;
m_MinPBSeqLen = cDfltMinPBSeqLen;	
m_MinPBSeqOverlap = cDfltMinErrCorrectLen;
m_MaxArtefactDev = cDfltMaxArtefactDev;

m_DeltaCoreOfs = cDfltDeltaCoreOfs;
m_MaxSeedCoreDepth = cDfltMaxSeedCoreDepth;
m_SeedCoreLen = cDfltSeedCoreLen;
m_MinNumSeedCores = cDfltNumSeedCores;
}

int
CSWAlign::CreateLocks(void)
{
if(m_bLocksCreated)
	return(eBSFSuccess);

#ifdef _WIN32
InitializeSRWLock(&m_hRwLock);
#else
if(pthread_rwlock_init (&m_hRwLock,NULL)!=0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to create rwlock");
	return(eBSFerrInternal);
	}
#endif
m_bLocksCreated = true;
return(eBSFSuccess);
}

void
CSWAlign::DeleteLocks(void)
{
if(!m_bLocksCreated)
	return;
#ifndef _WIN32
pthread_rwlock_destroy(&m_hRwLock);
#endif
m_bLocksCreated = false;
}

void
CSWAlign::AcquireLock(bool bExclusive)
{
#ifdef _WIN32
if(bExclusive)
	AcquireSRWLockExclusive(&m_hRwLock);
else
	AcquireSRWLockShared(&m_hRwLock);
#else
if(bExclusive)
	pthread_rwlock_wrlock(&m_hRwLock);
else
	pthread_rwlock_rdlock(&m_hRwLock);
#endif
}

void
CSWAlign::ReleaseLock(bool bExclusive)
{
#ifdef _WIN32
if(bExclusive)
	ReleaseSRWLockExclusive(&m_hRwLock);
else
	ReleaseSRWLockShared(&m_hRwLock);
#else
pthread_rwlock_unlock(&m_hRwLock);
#endif
}


int 
CSWAlign::Initialise(UINT32 MaxSWAInstances,	// initialise for this many instances
			UINT32 MinOverlapLen,				// the putative overlap would be of at least this length
			UINT32 MaxSeedCoreDepth,			// only further extend a seed core if there are no more than this number of matching cores in all targeted sequences
			UINT32 DeltaCoreOfs,				// offset core windows of coreSeqLen along the probe sequence when checking for overlaps 
			UINT32 CoreSeqLen,					// putative overlaps are explored if there are cores of at least this length in any putative overlap
			UINT32 MinNumCores,					// and if the putative overlap contains at least this many cores
			UINT32 MaxAcceptHitsPerSeedCore,	// limit accepted hits per seed core to no more than this many
			UINT32 DfltMaxProbeSeqLen)			// initially allocate for this length probe sequence to be aligned, will be realloc'd as may be required
{
Reset();
CreateLocks();
m_MinOverlapLen = MinOverlapLen;
m_MaxSeedCoreDepth = MaxSeedCoreDepth;
m_DeltaCoreOfs = DeltaCoreOfs;
m_SeedCoreLen = CoreSeqLen;
m_MinNumCores = MinNumCores;
m_MaxAcceptHitsPerSeedCore = MaxAcceptHitsPerSeedCore;
m_DfltMaxProbeSeqLen = DfltMaxProbeSeqLen;

if((m_pSWInstances = new tsPBSSWInstance[MaxSWAInstances])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Init: Failed memory allocation for %u instances",MaxSWAInstances);
	Reset();
	return(eBSFerrMem);
	}
memset(m_pSWInstances,0,sizeof(tsPBSSWInstance) * MaxSWAInstances);
m_AllocdSWInstances = MaxSWAInstances;
m_NumSWInstances = 0;
return(eBSFSuccess);
}

int										// returns number of target sequences loaded and indexed, or error result if < 0
CSWAlign::LoadTargetSeqs(int MinSeqLen,	// only accept target sequences of at least this length
			char *pszTargSeqsFile,		// load target sequences from this file
			int NumThreads)				// max number of threads when indexing
{
int Rslt;
int NumEntries;
int NumDupEntries;
UINT32 DupEntries[20];
char szDupEntry[100];

INT64 SumFileSizes;
SumFileSizes = 0;
#ifdef _WIN32
struct _stat64 st;
if(!_stat64(pszTargSeqsFile,&st))
#else
struct stat64 st;
if(!stat64(pszTargSeqsFile,&st))
#endif
	SumFileSizes = (INT64)st.st_size;
if(SumFileSizes < (INT64)MinSeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Loading target sequences, unable to access '%s' or contaminate sequence length smaller than %d",pszTargSeqsFile,MinSeqLen);
	Reset();
	return(eBSFerrFileAccess);
	}

if(m_pSfxArray != NULL)
	delete m_pSfxArray;
if((m_pSfxArray = new CSfxArrayV3) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Loading target sequences, unable to instantiate instance of CSfxArrayV3");
	return(eBSFerrObj);
	}
m_pSfxArray->Reset(false);
m_pSfxArray->SetMaxQSortThreads(NumThreads);

Rslt=m_pSfxArray->Open(false,false);
if(Rslt !=eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Loading target sequences, unable to create in-memory suffix array - error %d %s",Rslt,m_pSfxArray->GetErrMsg());
	Reset();
	return(Rslt);
	}

if((Rslt=m_pSfxArray->SetDescription((char *)"inmem")) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Loading target sequences, Unable to set description 'inmem'");
	Reset();
	return(Rslt);
	}
if((Rslt=m_pSfxArray->SetTitle((char *)"inmem")) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Loading target sequences, Unable to set title 'inmem'");
	Reset();
	return(Rslt);
	}

if((Rslt = m_pSfxArray->SetDatasetName((char *)"inmem")) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Loading target sequences, Unable to set dataset name 'inmem'");
	Reset();
	return(Rslt);
	}
m_pSfxArray->SetInitalSfxAllocEls(SumFileSizes);	// just a hint which is used for initial allocations by suffix processing
Rslt=LoadTargFastaFile(MinSeqLen,pszTargSeqsFile);
if(Rslt < eBSFSuccess)
	{
	m_pSfxArray->Close(false);
	Reset();
	return(Rslt);
	}
NumEntries = m_pSfxArray->GetNumEntries();
if(NumEntries < 1)
	return(NumEntries);
m_MaxTargSeqLen = m_pSfxArray->GetMaxSeqLen();

	// check for duplicate entry names
if(NumEntries > 1 && (NumDupEntries = m_pSfxArray->ChkDupEntries(20,&DupEntries[0])) > 0)
	{
	while(NumDupEntries--)
		{
		m_pSfxArray->GetIdentName(DupEntries[NumDupEntries],sizeof(szDupEntry)-1,szDupEntry); // get sequence name for specified entry identifier
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Loading contaminate sequences, duplicate sequence entry name '%s' in file '%s'",szDupEntry,pszTargSeqsFile);
		}
	m_pSfxArray->Close(false);
	Reset();
	return(eBSFerrFastqSeqID);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"CreateBioseqSuffixFile: sorting suffix array...");
m_pSfxArray->Finalise();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"CreateBioseqSuffixFile: sorting completed");

if(m_SeedCoreLen <= cMaxKmerLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Initialising for over occurring ( > %d) K-mers of length %d",m_MaxSeedCoreDepth,m_SeedCoreLen);
	if((Rslt = m_pSfxArray->InitOverOccKMers((int)m_SeedCoreLen,m_MaxSeedCoreDepth))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Failed to initialise for over occurring K-mers");
		Reset();
		return(Rslt);
		}
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Initialised for over occurring K-mers");
	}

return(NumEntries);
}

// SetScores
// Set match, mismatch, gap opening and gap extension scores
bool
CSWAlign::SetScores(int MatchScore,			// score for match
				int MismatchPenalty,	// penalty for mismatch
				int GapOpenPenalty,		// penalty for opening a gap
				int GapExtnPenalty,		// penalty if extending already opened gap
				int DlyGapExtn,			// delayed gap penalties, only apply gap extension penalty if gap at least this length
				int ProgPenaliseGapExtn, // if non-zero then progressively increment gap extension penalty for gaps of at least this length, 0 to disable, used for PacBio style error profiles
				int AnchorLen)			 // identified first and last anchors in alignment to be of at least this length
{
if (MatchScore <= 0 || MatchScore > 100 ||
	DlyGapExtn < 1 || DlyGapExtn > cSSWMaxPenaltyGap ||
	ProgPenaliseGapExtn < 0 || ProgPenaliseGapExtn > cSSWMaxPenaltyGap ||
	MismatchPenalty < -100 || MismatchPenalty > 0 ||
	GapOpenPenalty < -100 || GapOpenPenalty > 0 ||
	GapExtnPenalty < -100 || GapExtnPenalty > 0 ||
	AnchorLen < cSSWMinAnchorLen || AnchorLen > cSSWMaxAnchorLen)
	return(false);
m_SWMatchScore = MatchScore;
m_SWMismatchPenalty = MismatchPenalty;
m_SWGapOpenPenalty = GapOpenPenalty;
m_SWGapExtnPenalty = GapExtnPenalty;
m_SWDlyGapExtn = DlyGapExtn;
if (ProgPenaliseGapExtn > 0 && ProgPenaliseGapExtn < DlyGapExtn)
	ProgPenaliseGapExtn = DlyGapExtn;
m_SWProgPenaliseGapExtn = ProgPenaliseGapExtn;
m_SWAnchorLen = AnchorLen;
return(true);
}


// SetCPScores
// Set match, mismatch, gap opening and gap extension scores
bool
CSWAlign::SetCPScores(int MatchScore,			// score for match
				  int MismatchPenalty,	// penalty for mismatch
				  int GapOpenPenalty,		// penalty for opening a gap
				  int GapExtnPenalty)		// penalty if extending already opened gap
{
if (MatchScore <= 0 || MatchScore > 100 ||
	MismatchPenalty < -100 || MismatchPenalty > 0 ||
	GapOpenPenalty < -100 || GapOpenPenalty > 0 ||
	GapExtnPenalty < -100 || GapExtnPenalty > 0)
	return(false);
m_CPMatchScore = MatchScore;
m_CPMismatchPenalty = MismatchPenalty;
m_CPGapOpenPenalty = GapOpenPenalty;
m_CPGapExtnPenalty = GapExtnPenalty;
return(true);
}

bool
CSWAlign::SetMaxInitiatePathOfs(int MaxInitiatePathOfs,	// require SW paths to have started within this many bp (0 to disable) on either the probe or target - effectively anchoring the SW 
					int OverlapFloat)		// with this overlap float 
{
if (MaxInitiatePathOfs < 0 || MaxInitiatePathOfs > 10000)
	return(false);
m_MaxInitiatePathOfs = MaxInitiatePathOfs;
m_OverlapFloat = OverlapFloat;
return(true);
}

UINT32			// returned instance identifier
CSWAlign::InitInstance(void)
{
tsPBSSWInstance *pSWInstance;
AcquireLock(true);
if(m_NumSWInstances == m_AllocdSWInstances)
	{
	ReleaseLock(true);
	return(0);
	}
pSWInstance = &m_pSWInstances[m_NumSWInstances++];
memset(pSWInstance,0,sizeof(tsPBSSWInstance));
pSWInstance->InstanceID = m_NumSWInstances;
ReleaseLock(true);

pSWInstance->AllocdTargSeqSize = m_MaxTargSeqLen + 10;
if((pSWInstance->pTargSeq = new etSeqBase [pSWInstance->AllocdTargSeqSize])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"InitInstance: target sequence memory allocation of %d bytes - %s",pSWInstance->AllocdTargSeqSize,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}

pSWInstance->AllocdProbeSeqSize = m_DfltMaxProbeSeqLen + 10;
if((pSWInstance->pProbeSeq = new etSeqBase [pSWInstance->AllocdProbeSeqSize])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"InitInstance: probe sequence memory allocation of %d bytes - %s",pSWInstance->AllocdProbeSeqSize,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}

if((pSWInstance->pmtqsort = new CMTqsort) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"InitInstance: Core hits instantiation of CMTqsort failed");
	Reset();
	return(eBSFerrObj);
	}
pSWInstance->pmtqsort->SetMaxThreads(4);

pSWInstance->pSW = new CSSW;
pSWInstance->pSW->SetScores(m_SWMatchScore,m_SWMismatchPenalty,m_SWGapOpenPenalty,m_SWGapExtnPenalty,m_SWProgExtnPenaltyLen,min(63,m_SWProgExtnPenaltyLen+3), m_SWAnchorLen);
pSWInstance->pSW->SetCPScores(m_CPMatchScore, m_CPMismatchPenalty, m_CPGapOpenPenalty, m_CPGapExtnPenalty);
pSWInstance->pSW->SetMaxInitiatePathOfs(m_MaxInitiatePathOfs); 
pSWInstance->pSW->PreAllocMaxTargLen(pSWInstance->AllocdTargSeqSize);

pSWInstance->AllocdCoreHits = cAllocdNumCoreHits;
pSWInstance->AllocdCoreHitsSize = cAllocdNumCoreHits * sizeof(tsPBSSWACoreHit);
#ifdef _WIN32
pSWInstance->pCoreHits = (tsPBSSWACoreHit *)malloc(pSWInstance->AllocdCoreHitsSize);	
if(pSWInstance->pCoreHits == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"InitInstance: probe sequence memory allocation of %d bytes - %s",pSWInstance->AllocdCoreHitsSize,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
if((pSWInstance->pCoreHits = (tsPBSSWACoreHit *)mmap(NULL,pSWInstance->AllocdCoreHitsSize, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0)) == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"InitInstance: probe sequence memory allocation of %d bytes - %s",pSWInstance->AllocdCoreHitsSize,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#endif

return(pSWInstance->InstanceID);
}

// LoadTargFastaFile
// Parse input fasta format target sequence file into a biosequence suffix array file
int
CSWAlign::LoadTargFastaFile(int MinSeqLen,		// only accept for indexing sequences of at least this length
				char *pszFile,						// file containing sequences
				int Flags)							// default is for flags = cFlgLCSeq used with PacBio read sequences
{
CFasta Fasta;
unsigned char *pSeqBuff;
unsigned char *pMskBase;
UINT32 MskIdx;
size_t BuffOfs;
size_t AllocdBuffSize;
size_t AvailBuffSize;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
UINT32 SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
int Rslt;
int SeqID;
int NumSeqsAccepted;
size_t TotAcceptedLen;
UINT32 NumSeqsUnderlength;

if((Rslt=Fasta.Open(pszFile,true))!=eBSFSuccess)
	{
	if(Rslt != eBSFerrNotFasta)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile: Unable to open '%s' [%s] %s",pszFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	return(Rslt);
	}

AllocdBuffSize = (size_t)cAllocTargetSeqSize * 16;
// note malloc is used as can then simply realloc to expand as may later be required
if((pSeqBuff = (unsigned char *)malloc(AllocdBuffSize)) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to allocate memory (%u bytes) for sequence buffer",(UINT32)AllocdBuffSize);
	Fasta.Close();
	return(eBSFerrMem);
	}
AvailBuffSize = AllocdBuffSize;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile:- Adding %s..",pszFile);

bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
BuffOfs = 0;
NumSeqsUnderlength = 0;
NumSeqsAccepted = 0;
TotAcceptedLen = 0;
while((Rslt = SeqLen = Fasta.ReadSequence(&pSeqBuff[BuffOfs],(int)min(AvailBuffSize,(size_t)cAllocTargetSeqSize),true,false)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if(bEntryCreated)				// add any previous entry
			{
			if(BuffOfs < (size_t)MinSeqLen)
				NumSeqsUnderlength += 1;
			else
				{
				if((Rslt=m_pSfxArray->AddEntry(szName,pSeqBuff,(UINT32)BuffOfs,Flags)) < eBSFSuccess)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d %s",Rslt,m_pSfxArray->GetErrMsg());
					break;
					}
				else
					{
					NumSeqsAccepted += 1;
					TotAcceptedLen += BuffOfs;
					}
				}
			Rslt = eBSFSuccess;
			}
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		// An assumption - will one day bite real hard - is that the
		// fasta descriptor line starts with some form of unique identifier.
		// Use this identifier as the entry name.
		if(sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",pszFile,++SeqID);

		bFirstEntry = false;
		bEntryCreated = true;
		BuffOfs = 0;
		continue;
		}
	else
		if(bFirstEntry)	// if there was no descriptor then dummy up one...
			{
			SeqID++;
			sprintf(szName,"%s.%d",pszFile,SeqID);
			strcpy(szDescription,"No Description provided");
			bFirstEntry = false;
			bEntryCreated = true;
			}

	// remove any repeat masking flags so that sorts can actually sort
	// if run of more than 25 Ns and at least 5 Ns to end of buffer then randomly mutate
	// every 13th N
	//	e.g <25Ns>r<12Ns>r<12Ns> where r is a pseudorandom base
	pMskBase = &pSeqBuff[BuffOfs];
	int SeqNs = 0;
	for(MskIdx = 0; MskIdx < SeqLen; MskIdx++,pMskBase++)
		{
		*pMskBase &= ~cRptMskFlg;
		if(*pMskBase == eBaseN && (MskIdx+5) < SeqLen)
			{
			if(++SeqNs > 25 &&
				pMskBase[1] == eBaseN &&
				pMskBase[2] == eBaseN &&
				pMskBase[3] == eBaseN &&
				pMskBase[4] == eBaseN)
				{
				if(!(SeqNs % 13))	// mutate every 13th
					*pMskBase = rand() % 4;
				}
			}
		else
			SeqNs = 0;
		}

	BuffOfs += SeqLen;
	AvailBuffSize -= SeqLen;
	if(AvailBuffSize < (size_t)(cAllocTargetSeqSize / 8))
		{
		size_t NewSize = (size_t)cAllocTargetSeqSize + AllocdBuffSize;
		unsigned char *pTmp;
		if((pTmp = (unsigned char *)realloc(pSeqBuff,NewSize))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile:- Unable to reallocate memory (%u bytes) for sequence buffer",(UINT32)NewSize);
			return(eBSFerrMem);
			}
		pSeqBuff = pTmp;
		AllocdBuffSize = NewSize;
		AvailBuffSize = AllocdBuffSize - BuffOfs;
		}
	}
if(Rslt < eBSFSuccess && Rslt != eBSErrSession)
	{
	while(Fasta.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,Fasta.GetErrMsg());
	}

if(Rslt >= eBSFSuccess && bEntryCreated && BuffOfs > 0)			// close entry
	{
	if(BuffOfs < (size_t)MinSeqLen)
		{
		NumSeqsUnderlength += 1;
		Rslt = eBSFSuccess;
		}
	else
		{
		if((Rslt=m_pSfxArray->AddEntry(szName,pSeqBuff,(UINT32)BuffOfs,Flags)) < eBSFSuccess)
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d %s",Rslt,m_pSfxArray->GetErrMsg());
		else
			{
			Rslt = eBSFSuccess;
			NumSeqsAccepted += 1;
			TotAcceptedLen += BuffOfs;
			}
		}
	}
if(pSeqBuff != NULL)
	free(pSeqBuff);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile - %d parsed, %d accepted, %dbp mean length, %d sequences not accepted for indexing as length under %dbp ",
					SeqID,NumSeqsAccepted,NumSeqsAccepted == 0 ? 0 :(int)(TotAcceptedLen/NumSeqsAccepted),NumSeqsUnderlength,MinSeqLen);
return(Rslt);
}

int					// returns index 1..N of just added core hit or -1 if errors
CSWAlign::AddCoreHit(tsPBSSWInstance *pInstance,		// using this instance
				bool bRevCpl,					// true if core sequence was revcpl'd before matching
			   UINT32 ProbeOfs,                 // hit started at this probe offset
			   UINT32 TargSeqID,                // probe core matched onto this target sequence
			   UINT32 TargOfs,                  // probe core matched starting at this target loci
			   UINT32 HitLen)					// hit was of this length
{
tsPBSSWACoreHit *pCoreHit;

if((pInstance->NumCoreHits + 5) > pInstance->AllocdCoreHits)	// need to realloc memory to hold additional cores?
	{
		// realloc memory with a 25% increase over previous allocation 
	int coresreq;
	size_t memreq;
	void *pAllocd;
	coresreq = (int)(((INT64)pInstance->AllocdCoreHits * 125) / (INT64)100);
	memreq = coresreq * sizeof(tsPBSSWACoreHit);

#ifdef _WIN32
		pAllocd = realloc(pInstance->pCoreHits,memreq);
#else
		pAllocd = mremap(pInstance->pCoreHits,pInstance->AllocdCoreHitsSize,memreq,MREMAP_MAYMOVE);
		if(pAllocd == MAP_FAILED)
			pAllocd = NULL;
#endif
		if(pAllocd == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"SavePartialSeqs: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
			Reset();
			return(eBSFerrMem);
			}
		pInstance->pCoreHits = (tsPBSSWACoreHit *)pAllocd;
		pInstance->AllocdCoreHitsSize = memreq;
		pInstance->AllocdCoreHits = coresreq; 
		}
		
pCoreHit = &pInstance->pCoreHits[pInstance->NumCoreHits++];
pCoreHit->flgRevCpl = bRevCpl ? 1 : 0;
pCoreHit->flgMulti = 0;
pCoreHit->ProbeOfs = ProbeOfs;
pCoreHit->TargSeqID = TargSeqID;
pCoreHit->HitLen = HitLen;
pCoreHit->TargOfs = TargOfs;
memset(&pCoreHit[1],0,sizeof(tsPBSSWACoreHit));	// ensuring that used cores are always terminated with a marker end of cores initialised to 0
return(pInstance->NumCoreHits);
}

int									// returns index 1..N of core hits remaining or -1 if errors
CSWAlign::RemoveAddedCoreHits(tsPBSSWInstance *pInstance,		// using this instance
				int NumToRemove)                   // removing the last NumToRemove AddCoreHit() added
{

if(pInstance->NumCoreHits > (UINT32)NumToRemove)
	{
	pInstance->NumCoreHits -= NumToRemove;
	memset(&pInstance->pCoreHits[pInstance->NumCoreHits],0,sizeof(tsPBSSWACoreHit));	// ensuring that used cores are always terminated with a marker end of cores initialised to 0
	return(pInstance->NumCoreHits);
	}
pInstance->NumCoreHits = 0;
memset(pInstance->pCoreHits,0,sizeof(tsPBSSWACoreHit));	// ensuring that used cores are always terminated with a marker end of cores initialised to 0
return(-1);
}

int
CSWAlign::IdentifyCoreHits(tsPBSSWInstance *pInstance,		// using this instance
				bool bRevCpl)		// true if probe sequence to be reverse complemented
{
INT64 PrevHitIdx;
INT64 NextHitIdx;
UINT32 HitEntryID;
UINT32 HitLoci;
UINT32 HitsThisCore;
UINT32 HighHitsThisCore;
UINT32 TotHitsAllCores;
UINT32 ProbeOfs;
UINT32 LastProbeOfs;

UINT32 ChkOvrLapCoreProbeOfs;
UINT32 LastCoreProbeOfs;
int ChkOvrLapCoreStartIdx;

int TooManyHits = 0;

etSeqBase *pCoreSeq;

if(bRevCpl)
	CSeqTrans::ReverseComplement(pInstance->ProbeSeqLen,pInstance->pProbeSeq);
pCoreSeq = pInstance->pProbeSeq;
ChkOvrLapCoreProbeOfs = 0;
ChkOvrLapCoreStartIdx = 0;
LastCoreProbeOfs = 0;
PrevHitIdx = 0;
HitsThisCore = 0;
HighHitsThisCore = 0;
TotHitsAllCores = 0;
LastProbeOfs = 1 + pInstance->ProbeSeqLen - m_SeedCoreLen;
if(m_SeedCoreLen < cMaxPacBioSeedExtn)
	LastProbeOfs -= 120;

for(ProbeOfs = 0; ProbeOfs < LastProbeOfs; ProbeOfs+=m_DeltaCoreOfs,pCoreSeq+=m_DeltaCoreOfs)
	{
	PrevHitIdx = 0;
	HitsThisCore = 0;

    while((NextHitIdx = m_pSfxArray->IteratePacBio(pCoreSeq,pInstance->ProbeSeqLen - ProbeOfs,m_SeedCoreLen,0,1,PrevHitIdx,&HitEntryID,&HitLoci)) > 0)
		{
		PrevHitIdx = NextHitIdx;
  		AddCoreHit(pInstance,bRevCpl,ProbeOfs,HitEntryID,HitLoci,m_SeedCoreLen);
		HitsThisCore += 1;
		if(HitsThisCore > m_MaxAcceptHitsPerSeedCore)
			{
			TooManyHits += 1;
			RemoveAddedCoreHits(pInstance,HitsThisCore);
			HitsThisCore = 0;
			break;
			}

		}
	if(HitsThisCore)	// if at least one hit from this core
		{
		if(HitsThisCore > HighHitsThisCore)
			HighHitsThisCore = HitsThisCore;
		TotHitsAllCores += HitsThisCore;
		}
	}
if(bRevCpl)
	CSeqTrans::ReverseComplement(pInstance->ProbeSeqLen,pInstance->pProbeSeq);

return(pInstance->NumCoreHits);
}

int												    // 0 if probe aligns to no target sequence, otherwise the target sequence ID
CSWAlign::AlignProbeSeq(UINT32 SWAInstance,         // alignment instance
						UINT32 ProbeSeqLen,			// sequence to align is this length
						etSeqBase *pProbeSeq,       // probe sequence to align
						bool bSenseOnly,		   // true if to align probe sense only, false to align both sense and antisense	
    					tsSSWCell *pRetMatched)    // optional (if not NULL) returned match detail

{
tsPBSSWInstance *pSWAInstance;

UINT32 TargLen;
UINT32 CurTargCoreHitCnts;

UINT32 ProbeAlignLength;
UINT32 TargAlignLength;
UINT32 TargSeqLen;
UINT32 LongSAligns;
UINT32 LongAAligns;
tsSSWCell *pPeakMatchesCell;
tsSSWCell PeakMatchesCell;
#ifdef _PEAKSCOREACCEPT
tsSSWCell PeakScoreCell;
#endif
bool bTargSense;
UINT32 Idx;
UINT32 CurSummaryHitCnts;
UINT32 LowestSummaryHitCnts;
sPBSSWCoreHitCnts *pLowestSummaryHitCnts;
UINT32 HitIdx;

tsPBSSWACoreHit *pCoreHit;
UINT32 CurTargSeqID;

UINT32 CurTargHitOfs;
UINT32 CurProbeHitOfs;
UINT32 CurSEntryIDHits;
UINT32 CurAEntryIDHits;
UINT32	CurSTargStartOfs;
UINT32	CurSTargEndOfs;
UINT32	CurATargStartOfs;
UINT32	CurATargEndOfs;
UINT32	CurSProbeStartOfs;
UINT32	CurSProbeEndOfs;
UINT32	CurAProbeStartOfs;
UINT32	CurAProbeEndOfs;
UINT32 ProvOverlapping;
UINT32 ProvOverlapped;
UINT32 ProvContained;
UINT32 ProvArtefact;
UINT32 ProvSWchecked;
UINT32 MinOverlapLen;
UINT32 AdjOverlapFloat;

tsPBSSWACoreHit *pFirstCoreHit;
tsPBSSWACoreHit *pNxtCoreHit;
tsPBSSWACoreHit *pMaxCoreHit;
UINT32 MaxWinSize;
UINT32 RelWinSize;
bool bFirstHitNewTargSeq;

sPBSSWCoreHitCnts *pSummaryCnts;
int NumInMultiAlignment;
int Class;
NumInMultiAlignment = 0;

if(pRetMatched != NULL)
	memset(pRetMatched,0,sizeof(tsSSWCell));
if(SWAInstance < 1 || SWAInstance > m_NumSWInstances)
	return(eBSFerrEntry);
pSWAInstance = &m_pSWInstances[SWAInstance-1];
AcquireLock(true);
if(pSWAInstance->FlgActive)
	{
	ReleaseLock(true);
	return(eBSFerrInternal);
	}
pSWAInstance->FlgActive = 1;
ReleaseLock(true);

pSWAInstance->NumCoreHits = 0;
pSWAInstance->NumTargCoreHitCnts = 0;
pSWAInstance->TargSeqLen = 0;
memset(pSWAInstance->TargCoreHitCnts,0,sizeof(sPBSSWCoreHitCnts));

if(pSWAInstance->pTargSeq == NULL || pSWAInstance->AllocdTargSeqSize < m_MaxTargSeqLen + 1)
	{
	if(pSWAInstance->pTargSeq != NULL)
		{
		delete pSWAInstance->pTargSeq;
		pSWAInstance->AllocdTargSeqSize = m_MaxTargSeqLen + 1;
		if((pSWAInstance->pTargSeq = new etSeqBase [pSWAInstance->AllocdTargSeqSize]) == NULL)
			return(eBSFerrMem);
		}
	}

if(pSWAInstance->pProbeSeq == NULL || pSWAInstance->AllocdProbeSeqSize < ProbeSeqLen + 1)
	{
	if(pSWAInstance->pProbeSeq != NULL)
		{
		delete pSWAInstance->pProbeSeq;
		pSWAInstance->AllocdProbeSeqSize = (ProbeSeqLen * 110) / 100;		// allow 10% more so as to reduce chances of requiring future reallocs
		if((pSWAInstance->pProbeSeq = new etSeqBase [pSWAInstance->AllocdProbeSeqSize]) == NULL)
			return(eBSFerrMem);
		}
	}
memcpy(pSWAInstance->pProbeSeq,pProbeSeq,ProbeSeqLen);
pSWAInstance->pProbeSeq[ProbeSeqLen] = eBaseEOS;
pSWAInstance->ProbeSeqLen = ProbeSeqLen;

IdentifyCoreHits(pSWAInstance,false);
if(!bSenseOnly)
	IdentifyCoreHits(pSWAInstance,true);

ProvOverlapping = 0;
ProvOverlapped = 0;
ProvContained = 0;
ProvArtefact = 0;
ProvSWchecked = 0;
AdjOverlapFloat = m_OverlapFloat + m_SeedCoreLen;
if(m_SeedCoreLen < cMaxPacBioSeedExtn)
	AdjOverlapFloat += 120;

pSWAInstance->NumTargCoreHitCnts = 0;
memset(pSWAInstance->TargCoreHitCnts,0,sizeof(pSWAInstance->TargCoreHitCnts));
if(pSWAInstance->NumCoreHits >= m_MinNumCores)
	{
			// resort core hits by TargNodeID.TargOfs.ProbeNodeID.ProbeOfs ascending
	pSWAInstance->pmtqsort->qsort(pSWAInstance->pCoreHits,pSWAInstance->NumCoreHits,sizeof(tsPBSSWACoreHit),SortCoreHitsByTargProbeOfs);
	pCoreHit = pSWAInstance->pCoreHits;
	for(HitIdx = 0; HitIdx < pSWAInstance->NumCoreHits; HitIdx++, pCoreHit++)
		pCoreHit->flgMulti = 0;
		
	CurSEntryIDHits = 0;
	CurAEntryIDHits = 0;
	CurTargHitOfs = 0;
	CurProbeHitOfs = 0;
	CurSTargStartOfs = 0;
	CurSTargEndOfs = 0;
	CurATargStartOfs = 0;
	CurATargEndOfs = 0;
	CurSProbeStartOfs = 0;
	CurSProbeEndOfs = 0;
	CurAProbeStartOfs = 0;
	CurAProbeEndOfs = 0;

	// with large target sequences then can have many artifactual core hits
	// process and mark these probable artifact hits by identifying the most spatially related cluster of hits; hits outside of
	// the most spatially related cluster are marked as being artefactual 

	UINT32 PrevAcceptedProbeOfs;
	UINT32 PrevAcceptedTargOfs;
	UINT32 MaxNoHitGapLen;
	UINT32 NumNoHitGaps;
	UINT32 SumNoHitGapLens;

	MaxNoHitGapLen = 1000;                 // expecting cores to be distributed such that there shouldn't be many separated by more than this intercore bp gap
	CurTargSeqID = 0;
	pFirstCoreHit = NULL;
	MaxWinSize = (ProbeSeqLen * 115) / 100;
	pCoreHit = pSWAInstance->pCoreHits;
	pMaxCoreHit = NULL;
	for (HitIdx = 0; HitIdx < pSWAInstance->NumCoreHits; HitIdx++, pCoreHit++)
		{
		if (CurTargSeqID == 0)    // 0 if 1st hit about to be processed for a new target sequence
			{
			pMaxCoreHit = NULL;
			pFirstCoreHit = pCoreHit;
			CurTargSeqID = pCoreHit->TargSeqID;
			TargSeqLen = m_pSfxArray->GetSeqLen(CurTargSeqID);
			}

		// if just checked last core hit for the current target ...
		if (HitIdx + 1 == pSWAInstance->NumCoreHits || pCoreHit[1].TargSeqID != CurTargSeqID)
			{
			while (pFirstCoreHit <= pCoreHit)
				{
				pFirstCoreHit->WinHits = 0;
				pFirstCoreHit->flgClustered = 0;
				pNxtCoreHit = pFirstCoreHit;
				RelWinSize = 0;
				PrevAcceptedProbeOfs = 0;
				PrevAcceptedTargOfs = 0;
				SumNoHitGapLens = 0;
				NumNoHitGaps = 0;

				if(pFirstCoreHit->ProbeOfs > MaxNoHitGapLen && pFirstCoreHit->TargOfs > MaxNoHitGapLen)
					{
					pFirstCoreHit += 1;
					continue;
					}
					
				while (pNxtCoreHit <= pCoreHit)
					{
					RelWinSize = pNxtCoreHit->TargOfs - pFirstCoreHit->TargOfs;
					if (RelWinSize > MaxWinSize)
						break;
					if (pNxtCoreHit->flgRevCpl == pFirstCoreHit->flgRevCpl)	// only interested in hits which are same sense as the first hit in window
						{
						if(PrevAcceptedProbeOfs == 0 || pNxtCoreHit->ProbeOfs >= PrevAcceptedProbeOfs + m_SeedCoreLen) // a single matched seed core extension may have resulted in multiple hits, reduce counts by requiring a differential of at least m_SeedCoreLen
							{
							if((pNxtCoreHit->ProbeOfs - PrevAcceptedProbeOfs) >= MaxNoHitGapLen)
								{
								NumNoHitGaps += 1;
								SumNoHitGapLens += pNxtCoreHit->ProbeOfs - PrevAcceptedProbeOfs;
								}
							PrevAcceptedProbeOfs = pNxtCoreHit->ProbeOfs;
							PrevAcceptedTargOfs = pNxtCoreHit->TargOfs;
							pFirstCoreHit->WinHits += 1;
							}
						}
					pNxtCoreHit += 1;
					}

				if(pFirstCoreHit->WinHits >= m_MinNumCores && ((PrevAcceptedProbeOfs + MaxNoHitGapLen) > ProbeSeqLen || (PrevAcceptedTargOfs + MaxNoHitGapLen) > TargSeqLen))
					{
					UINT32 PutativeOverlapLen = PrevAcceptedProbeOfs - pFirstCoreHit->ProbeOfs;
					if(((SumNoHitGapLens * 100) / PutativeOverlapLen) <= 20)		// only accepting if no more than 20% of alignment sums to gaps > MaxNoHitGapLen
						{
						if (pMaxCoreHit == NULL || pFirstCoreHit->WinHits > pMaxCoreHit->WinHits)
							pMaxCoreHit = pFirstCoreHit;
						}
					}
				pFirstCoreHit += 1;
				}

			if (pMaxCoreHit != NULL)
			{
				RelWinSize = 0;
				pNxtCoreHit = pMaxCoreHit;
				while (pNxtCoreHit <= pCoreHit)
				{
					RelWinSize = pNxtCoreHit->TargOfs - pMaxCoreHit->TargOfs;
					if (RelWinSize > MaxWinSize)
						break;
					if (pNxtCoreHit->flgRevCpl == pMaxCoreHit->flgRevCpl)	// only interested in hits which are same sense as the first hit in window
						pNxtCoreHit->flgClustered = 1;
					pNxtCoreHit->flgMulti = 0;
					pNxtCoreHit += 1;
				}
			}
			CurTargSeqID = 0;  // looking for cores on a new target
		}
	}

	// iterate and count hits for each TargNodeID whilst recording the loci of the first and last hit so can determine if overlap is a likely artefact
	CurSEntryIDHits = 0;
	CurAEntryIDHits = 0;
	CurTargSeqID = 0;
	CurTargHitOfs = 0;
	CurProbeHitOfs = 0;
	CurSTargStartOfs = 0;
	CurSTargEndOfs = 0;
	CurATargStartOfs = 0;
	CurATargEndOfs = 0;
	CurSProbeStartOfs = 0;
	CurSProbeEndOfs = 0;
	CurAProbeStartOfs = 0;
	CurAProbeEndOfs = 0;

	bFirstHitNewTargSeq = false;
	pCoreHit = pSWAInstance->pCoreHits;
	for(HitIdx = 0; HitIdx < pSWAInstance->NumCoreHits; HitIdx++,pCoreHit++)
		{
		if(CurTargSeqID == 0)    // 0 if 1st hit about to be processed for a new target sequence
			{
			bFirstHitNewTargSeq = true;
			CurTargSeqID = pCoreHit->TargSeqID;
			CurSEntryIDHits = 0;
			CurAEntryIDHits = 0;
			CurSTargStartOfs = 0;
			CurSTargEndOfs = 0;
			CurATargStartOfs = 0;
			CurATargEndOfs = 0;
			CurSProbeStartOfs = 0;
			CurSProbeEndOfs = 0;
			CurAProbeStartOfs = 0;
			CurAProbeEndOfs = 0;
			}

		if(pCoreHit->flgClustered && pCoreHit->TargSeqID == CurTargSeqID) // same target sequence so check for starting/ending offsets and accumulate hit counts 
			{
			CurTargHitOfs = pCoreHit->TargOfs;
			CurProbeHitOfs = pCoreHit->ProbeOfs;
			if(pCoreHit->flgRevCpl == 0)
				{
				if(bFirstHitNewTargSeq == true || CurTargHitOfs < CurSTargStartOfs)
					CurSTargStartOfs = CurTargHitOfs;
				if(CurTargHitOfs > CurSTargEndOfs)
					CurSTargEndOfs = CurTargHitOfs;	
				if(bFirstHitNewTargSeq == true || CurProbeHitOfs < CurSProbeStartOfs)
					CurSProbeStartOfs = CurProbeHitOfs;
				if(CurProbeHitOfs > CurSProbeEndOfs)
					CurSProbeEndOfs = CurProbeHitOfs;
				}
			else
				{
				if(bFirstHitNewTargSeq == true || CurTargHitOfs < CurATargStartOfs)
					CurATargStartOfs = CurTargHitOfs;
				if(CurTargHitOfs > CurATargEndOfs)
					CurATargEndOfs = CurTargHitOfs;	
				if(bFirstHitNewTargSeq == true || CurProbeHitOfs < CurAProbeStartOfs)
					CurAProbeStartOfs = CurProbeHitOfs;
				if(CurProbeHitOfs > CurAProbeEndOfs)
					CurAProbeEndOfs = CurProbeHitOfs;
				}
			bFirstHitNewTargSeq = false;
			if(pCoreHit->flgMulti != 1)
				{
				if(pCoreHit->flgRevCpl == 0)
					CurSEntryIDHits += 1;
				else
					CurAEntryIDHits += 1;
				}
			}

		// if just processed last core hit for the current target ...
		if(HitIdx + 1 == pSWAInstance->NumCoreHits || pCoreHit[1].TargSeqID != CurTargSeqID)
			{
			// checking here that the first and last hit are consistent with either a completely contained or overlapped
			TargLen = m_pSfxArray->GetSeqLen(pCoreHit->TargSeqID);
			if(CurSEntryIDHits >= m_MinNumCores) 
				{
				if((CurSProbeStartOfs >= m_OverlapFloat &&  CurSTargStartOfs >= m_OverlapFloat) ||
						((TargLen - CurSTargEndOfs) >= AdjOverlapFloat && (ProbeSeqLen - CurSProbeEndOfs) >= AdjOverlapFloat))
					CurSEntryIDHits = 0;
				}
			else
				CurSEntryIDHits = 0;

			if(CurAEntryIDHits >= m_MinNumCores)
				{
				if((CurAProbeStartOfs >= m_OverlapFloat && CurATargStartOfs >= m_OverlapFloat) ||
					((TargLen - CurATargEndOfs) >= AdjOverlapFloat && (ProbeSeqLen - CurAProbeEndOfs) >= AdjOverlapFloat))
					CurAEntryIDHits = 0;
				}
			else
				CurAEntryIDHits = 0;

			if(CurSEntryIDHits >= m_MinNumCores || CurAEntryIDHits >= m_MinNumCores)
				{
				if(pSWAInstance->NumTargCoreHitCnts == cPBSSWSummaryTargCoreHitCnts)
					{
					LowestSummaryHitCnts = 0;
					pLowestSummaryHitCnts = NULL;
					pSummaryCnts = pSWAInstance->TargCoreHitCnts; 
					for(Idx = 0; Idx < cPBSSWSummaryTargCoreHitCnts; Idx++, pSummaryCnts++)
						{
						CurSummaryHitCnts = pSummaryCnts->NumSHits + pSummaryCnts->NumAHits;
						if(LowestSummaryHitCnts == 0 || CurSummaryHitCnts < LowestSummaryHitCnts)
							{
							LowestSummaryHitCnts = CurSummaryHitCnts;
							pLowestSummaryHitCnts = pSummaryCnts;
							}
						}

					if((CurSEntryIDHits + CurAEntryIDHits) <= LowestSummaryHitCnts)
						{
						CurTargSeqID = 0; 
						continue;
						}
					pSummaryCnts = pLowestSummaryHitCnts;
					}
				else
					pSummaryCnts = &pSWAInstance->TargCoreHitCnts[pSWAInstance->NumTargCoreHitCnts++];
				pSummaryCnts->TargSeqID = CurTargSeqID;
				pSummaryCnts->STargStartOfs = CurSTargStartOfs;
				pSummaryCnts->STargEndOfs = CurSTargEndOfs;
				pSummaryCnts->ATargStartOfs = CurATargStartOfs;
				pSummaryCnts->ATargEndOfs = CurATargEndOfs;
				pSummaryCnts->SProbeStartOfs = CurSProbeStartOfs;
				pSummaryCnts->SProbeEndOfs = CurSProbeEndOfs;
				pSummaryCnts->AProbeStartOfs = CurAProbeStartOfs;
				pSummaryCnts->AProbeEndOfs = CurAProbeEndOfs;
				pSummaryCnts->NumSHits = CurSEntryIDHits;
				pSummaryCnts->NumAHits = CurAEntryIDHits;
				}
			CurTargSeqID = 0; 
			}
		}
	}
 
	// can't process, SW over all would be too resource intensive, all targets which meet the minimum number of core hits requested so choose the top cMaxProbePBSSWs as ranked by the number of core hits
	if(pSWAInstance->NumTargCoreHitCnts > 1)
		{
		pSWAInstance->pmtqsort->qsort(pSWAInstance->TargCoreHitCnts,pSWAInstance->NumTargCoreHitCnts,sizeof(sPBSSWCoreHitCnts),SortCoreHitsDescending);
		if(pSWAInstance->NumTargCoreHitCnts > cMaxProbePBSSWs)		// clamp to no more than this many SW alignments
			pSWAInstance->NumTargCoreHitCnts = cMaxProbePBSSWs;
		}

	NumInMultiAlignment = 0;
	if(pSWAInstance->NumTargCoreHitCnts > 0)
		{
		LongSAligns = 0;
		LongAAligns = 0;

		pSWAInstance->pSW->SetProbe(pSWAInstance->ProbeSeqLen,pSWAInstance->pProbeSeq);
		pSummaryCnts = &pSWAInstance->TargCoreHitCnts[0];
		NumInMultiAlignment = 0;
		for(CurTargCoreHitCnts = 0; CurTargCoreHitCnts < pSWAInstance->NumTargCoreHitCnts; CurTargCoreHitCnts++,pSummaryCnts++)
			{
			if(pSummaryCnts->NumSHits < m_MinNumCores && pSummaryCnts->NumAHits < m_MinNumCores)
				continue;
			bTargSense = pSummaryCnts->NumSHits >= pSummaryCnts->NumAHits ? true :  false;
			MinOverlapLen = m_MinOverlapLen;
			TargSeqLen = m_pSfxArray->GetSeqLen(pSummaryCnts->TargSeqID); 
			m_pSfxArray->GetSeq(pSummaryCnts->TargSeqID,0,pSWAInstance->pTargSeq,TargSeqLen);
			if(!bTargSense)
				CSeqTrans::ReverseComplement(TargSeqLen,pSWAInstance->pTargSeq);
			pSWAInstance->pTargSeq[TargSeqLen] = eBaseEOS;
			pSWAInstance->pSW->SetTarg(TargSeqLen,pSWAInstance->pTargSeq);
			
			// restrict the range over which the SW will be processed to that of the overlap +/- m_OverlapFloat
			int Rslt;
			if(bTargSense)
				{
				if(pSummaryCnts->SProbeStartOfs < m_OverlapFloat)
					pSummaryCnts->SProbeStartOfs = 0;
				else
					pSummaryCnts->SProbeStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->SProbeEndOfs + AdjOverlapFloat >= ProbeSeqLen)
					pSummaryCnts->SProbeEndOfs = ProbeSeqLen - 1;
				else
					pSummaryCnts->SProbeEndOfs += AdjOverlapFloat;
				if(pSummaryCnts->STargStartOfs < m_OverlapFloat)
					pSummaryCnts->STargStartOfs = 0;
				else
					pSummaryCnts->STargStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->STargEndOfs + AdjOverlapFloat >= TargSeqLen)
					pSummaryCnts->STargEndOfs = TargSeqLen - 1;
				else
					pSummaryCnts->STargEndOfs += AdjOverlapFloat;
				Rslt = pSWAInstance->pSW->SetAlignRange(pSummaryCnts->SProbeStartOfs,pSummaryCnts->STargStartOfs,
											pSummaryCnts->SProbeEndOfs + 1 - pSummaryCnts->SProbeStartOfs,pSummaryCnts->STargEndOfs + 1 - pSummaryCnts->STargStartOfs);
				}
			else
				{
				UINT32 Xchg;
				Xchg = pSummaryCnts->AProbeStartOfs;
				pSummaryCnts->AProbeStartOfs = ProbeSeqLen - (pSummaryCnts->AProbeEndOfs + 1);
				pSummaryCnts->AProbeEndOfs = ProbeSeqLen - (Xchg + 1);
				Xchg = pSummaryCnts->ATargStartOfs;
				pSummaryCnts->ATargStartOfs = TargSeqLen - (pSummaryCnts->ATargEndOfs + 1);
				pSummaryCnts->ATargEndOfs = TargSeqLen - (Xchg + 1);

				if(pSummaryCnts->AProbeStartOfs < m_OverlapFloat)
					pSummaryCnts->AProbeStartOfs = 0;
				else
					pSummaryCnts->AProbeStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->AProbeEndOfs + AdjOverlapFloat >= ProbeSeqLen)
					pSummaryCnts->AProbeEndOfs = ProbeSeqLen - 1;
				else
					pSummaryCnts->AProbeEndOfs += AdjOverlapFloat;
				if(pSummaryCnts->ATargStartOfs < m_OverlapFloat)
					pSummaryCnts->ATargStartOfs = 0;
				else
					pSummaryCnts->ATargStartOfs -= m_OverlapFloat;
				if(pSummaryCnts->ATargEndOfs + AdjOverlapFloat >= TargSeqLen)
					pSummaryCnts->ATargEndOfs = TargSeqLen - 1;
				else
					pSummaryCnts->ATargEndOfs += AdjOverlapFloat;

				Rslt = pSWAInstance->pSW->SetAlignRange(pSummaryCnts->AProbeStartOfs,pSummaryCnts->ATargStartOfs,
											pSummaryCnts->AProbeEndOfs + 1 - pSummaryCnts->AProbeStartOfs,pSummaryCnts->ATargEndOfs + 1 - pSummaryCnts->ATargStartOfs);
				}

			pPeakMatchesCell = pSWAInstance->pSW->Align(NULL, min(ProbeSeqLen, m_MaxTargSeqLen));
			ProvSWchecked += 1;
			if(pPeakMatchesCell != NULL && pPeakMatchesCell->NumMatches >= (MinOverlapLen/2))
				{
				PeakMatchesCell = *pPeakMatchesCell;
				ProbeAlignLength = PeakMatchesCell.EndPOfs - PeakMatchesCell.StartPOfs + 1;
				TargAlignLength = PeakMatchesCell.EndTOfs - PeakMatchesCell.StartTOfs + 1;
				}
			else
				{
				memset(&PeakMatchesCell,0,sizeof(PeakMatchesCell));
				ProbeAlignLength = 0;
				TargAlignLength = 0;
				}

			if(((1+ ProbeAlignLength + TargAlignLength) / 2) >= MinOverlapLen)
				{
				ProvOverlapping += 1;

				// characterise the overlapped target
				// eOLCOverlapping if probe accepted as overlapping, either 5' or 3'
				// eOLCcontaining if both ends of target completely contained within probe
                // eOLCartefact if target is only partially contained
				int PathClass;
				Class = (int)ePBSSWOLCOverlapping;		
				if((PeakMatchesCell.StartTOfs >= m_OverlapFloat &&  PeakMatchesCell.StartPOfs >= m_OverlapFloat) ||
					 ((TargSeqLen - PeakMatchesCell.EndTOfs) >= m_OverlapFloat && (ProbeSeqLen - PeakMatchesCell.EndPOfs) >= m_OverlapFloat))
					{
					Class = (int)ePBSSWOLCartefact;
					ProvArtefact += 1;
					}

				if(Class == ePBSSWOLCOverlapping && (PathClass = pSWAInstance->pSW->ClassifyPath(m_MaxArtefactDev,
																					PeakMatchesCell.PFirstAnchorStartOfs,PeakMatchesCell.PLastAnchorEndOfs,
																					PeakMatchesCell.TFirstAnchorStartOfs,PeakMatchesCell.TLastAnchorEndOfs)) > 0)
					{
					Class = (int)ePBSSWOLCartefact;
					ProvArtefact += 1;
					}

				if(Class == ePBSSWOLCOverlapping)
					{
					if(PeakMatchesCell.StartTOfs < m_OverlapFloat && (TargSeqLen - PeakMatchesCell.EndTOfs) < m_OverlapFloat) // is target completely contained by probe?
						Class = (int)ePBSSWOLCcontains;
					else
						if(PeakMatchesCell.StartPOfs < m_OverlapFloat && (ProbeSeqLen - PeakMatchesCell.EndPOfs) < m_OverlapFloat) // or is probe completely contained within target?
							Class = (int)ePBSSWOLCcontained;
					
					AcquireLock(true);
					pSWAInstance->FlgActive = 0;
					ReleaseLock(true);
					if(pRetMatched != NULL)
						*pRetMatched = PeakMatchesCell;
					return(pSummaryCnts->TargSeqID);
					}
				ProvOverlapped += 1;
				}
			}
		}

AcquireLock(true);
if(pSWAInstance->FlgActive == 0)
	{
	ReleaseLock(true);
	return(eBSFerrInternal);
	}
pSWAInstance->FlgActive = 0;
ReleaseLock(true);
return(0);
}


// SortCoreHitsByProbeNodeID
// Sort core hits by ProbeNodeID.TargNodeID.ProbeOfs.TargOfs.flgRevCpl ascending
int
CSWAlign::SortCoreHitsByProbeTargOfs(const void *arg1, const void *arg2)
{
tsPBSSWACoreHit *pEl1 = (tsPBSSWACoreHit *)arg1;
tsPBSSWACoreHit *pEl2 = (tsPBSSWACoreHit *)arg2;

if(pEl1->TargSeqID < pEl2->TargSeqID)
	return(-1);
if(pEl1->TargSeqID > pEl2->TargSeqID)
	return(1);
if(pEl1->ProbeOfs < pEl2->ProbeOfs)	
	return(-1);
if(pEl1->ProbeOfs > pEl2->ProbeOfs)
	return(1);
if(pEl1->TargOfs < pEl2->TargOfs)	
	return(-1);
if(pEl1->TargOfs > pEl2->TargOfs)
	return(1);
if(pEl1->flgRevCpl != pEl2->flgRevCpl)
	return(1);
return(0);
}

// SortCoreHitsByTargNodeID
// Sort core hits by ProbeNodeID.TargNodeID.TargOfs.ProbeOfs.flgRevCpl ascending
int
CSWAlign::SortCoreHitsByTargProbeOfs(const void *arg1, const void *arg2)
{
tsPBSSWACoreHit *pEl1 = (tsPBSSWACoreHit *)arg1;
tsPBSSWACoreHit *pEl2 = (tsPBSSWACoreHit *)arg2;

if(pEl1->TargSeqID < pEl2->TargSeqID)
	return(-1);
if(pEl1->TargSeqID > pEl2->TargSeqID)
	return(1);
if(pEl1->TargOfs < pEl2->TargOfs)	
	return(-1);
if(pEl1->TargOfs > pEl2->TargOfs)
	return(1);
if(pEl1->ProbeOfs < pEl2->ProbeOfs)	
	return(-1);
if(pEl1->ProbeOfs > pEl2->ProbeOfs)
	return(1);
if(pEl1->flgRevCpl != pEl2->flgRevCpl)
	return(1);
return(0);
}

// SortCoreHitsDescending
// Sort target core hits by number of hits descending
int
CSWAlign::SortCoreHitsDescending(const void *arg1, const void *arg2)
{
sPBSSWCoreHitCnts *pEl1 = (sPBSSWCoreHitCnts *)arg1;
sPBSSWCoreHitCnts *pEl2 = (sPBSSWCoreHitCnts *)arg2;

UINT32 El1NumHits;
UINT32 El2NumHits;
El1NumHits = max(pEl1->NumSHits,pEl1->NumAHits);
El2NumHits = max(pEl2->NumSHits,pEl2->NumAHits);
if(El1NumHits > El2NumHits)
	return(-1);
if(El1NumHits < El2NumHits)
	return(1);
return(0);
}



