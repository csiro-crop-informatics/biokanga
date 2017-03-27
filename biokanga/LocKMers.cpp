/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

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

#include "./biokanga.h"
#include "LocKMers.h"

CLocKMers::CLocKMers(void)
{
m_pSfxArray = NULL;
m_pBlockSeqBuff = NULL;
m_pMarkerBuff = NULL;
m_pMarkerLenDist = NULL;
m_hOutFile = -1;
m_hOutReadsFile = -1;
#ifdef _WIN32
InitializeSRWLock(&m_hRwLock);
#else
pthread_rwlock_init( &m_hRwLock,NULL);
#endif
Reset();
}


CLocKMers::~CLocKMers(void)
{
if(m_hOutFile != -1)
	close(m_hOutFile);
if(m_hOutReadsFile != -1)
	close(m_hOutReadsFile);
if(m_pSfxArray != NULL)
	delete m_pSfxArray;
if(m_pBlockSeqBuff != NULL)
	delete m_pBlockSeqBuff;
if(m_pMarkerBuff != NULL)
	delete m_pMarkerBuff;
if(m_pMarkerLenDist != NULL)
	delete m_pMarkerLenDist;

#ifndef _WIN32
pthread_rwlock_destroy(&m_hRwLock);
#endif
}

void
CLocKMers::Reset(bool bSync)
{
if(m_pSfxArray != NULL)
	{
	delete m_pSfxArray;
	m_pSfxArray = NULL;
	}
if(m_pBlockSeqBuff != NULL)
	{
	delete m_pBlockSeqBuff;
	m_pBlockSeqBuff = NULL;
	}

if(m_hOutFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutFile);
#else
		fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if(m_hOutReadsFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutReadsFile);
#else
		fsync(m_hOutReadsFile);
#endif
	close(m_hOutReadsFile);
	m_hOutReadsFile = -1;
	}

if(m_pMarkerBuff != NULL)
	{
	delete m_pMarkerBuff;
	m_pMarkerBuff = NULL;
	}

if(m_pMarkerLenDist != NULL)
	{
	delete m_pMarkerLenDist;
	m_pMarkerLenDist = NULL;
	}

m_szDataset[0] = '\0';
m_szMarkerFile[0] = '\0';
m_szMarkerReadsFile[0] = '\0';
m_bKMerReads = false;
m_NumSfxEntries = 0;
m_NxtBlockSeqBuffOfs = 0;
m_BlockSeqBuffLen = 0;
m_AllocBlockSeqBuffSize = 0;
m_bInitSeqBuffering = true;
m_CurPartialCultivarID = 0;
m_CurPartialCultivarOfs = 0;
m_bAllEntrySeqsLoaded = false;
m_bAllBlocksReturned = false;
m_NumBlocks = 0;				
m_NumKMers = 0;					
m_NumPutativeKMers = 0;			
m_NumAcceptedKMers = 0;	
m_NumAcceptedExtdKMers = 0;
m_MarkerID = 0;		
m_MarkerBuffOfs = 0;			
m_AllocMarkerBuffSize = 0;	
m_MaxMarkerLen = 0;
m_MinMarkerLen = 0;
m_szTargetCultivar[0] = '\0';
m_NumPartialCultivars = 0;
memset(m_PartialCultivars,0,sizeof(m_PartialCultivars));
}

// AcquireLock
// Aquire lock, either exclusive or read only, on class shared instance vars
void 
CLocKMers::AcquireLock(bool bExclusive)
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


// ReleaseLock
// Release lock, must be exclusive or read only as used when AcquireLock was used
void
CLocKMers::ReleaseLock(bool bExclusive)
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
CLocKMers::ReportMarker(int MarkerLen,			// number of bases in marker sequence
			 etSeqBase *pMarkerSeq)				// marker sequence
{
int TruncMarkerLen;
char szMarkerSeq[cMaxExtKMerLen+1];				// additional for '\0' terminator

if(m_hOutFile == -1 || MarkerLen == 0 || pMarkerSeq == NULL || m_pMarkerBuff == NULL)	// better safe than sorry...
	return(0);

TruncMarkerLen = MarkerLen > cMaxExtKMerLen ? cMaxExtKMerLen : MarkerLen;
CSeqTrans::MapSeq2UCAscii(pMarkerSeq,TruncMarkerLen,szMarkerSeq);
AcquireLock(true);
if((m_MarkerBuffOfs + (2 * cMaxExtKMerLen)) > m_AllocMarkerBuffSize)
	{
	CUtility::SafeWrite(m_hOutFile,m_pMarkerBuff,m_MarkerBuffOfs);
	m_MarkerBuffOfs = 0;
	}
m_MarkerID += 1;
if(MarkerLen != TruncMarkerLen)					// if needed to clamp
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"ReportMarker: Truncated overlength (%d) Marker%d down to length %d",MarkerLen,m_MarkerID,cMaxExtKMerLen);
	MarkerLen = cMaxExtKMerLen;
	}
m_MarkerBuffOfs += sprintf((char *)&m_pMarkerBuff[m_MarkerBuffOfs],">%sM%d %d|%d|%d\n%s\n",m_szTargetCultivar,m_MarkerID,TruncMarkerLen,m_KMerLen,m_MinHamming,szMarkerSeq);
m_pMarkerLenDist[TruncMarkerLen-1] += 1;
if(m_MaxMarkerLen < TruncMarkerLen)
	m_MaxMarkerLen = TruncMarkerLen;
if(m_MinMarkerLen == 0 || m_MinMarkerLen > TruncMarkerLen)
	m_MinMarkerLen = TruncMarkerLen;
ReleaseLock(true);
return(eBSFSuccess);
}

int
CLocKMers::ReportMarkerRead(int ReadLen,		// number of bases in read sequence
			 etSeqBase *pReadSeq)				// read sequence
{
int TruncMarkerLen;
char szMarkerSeq[cMaxExtKMerLen+1];				// additional for '\0' terminator

if(m_hOutReadsFile == -1 || ReadLen == 0 || pReadSeq == NULL || m_pMarkerBuff == NULL)	// better safe than sorry...
	return(0);

TruncMarkerLen = ReadLen > cMaxExtKMerLen ? cMaxExtKMerLen : ReadLen;
CSeqTrans::MapSeq2UCAscii(pReadSeq,TruncMarkerLen,szMarkerSeq);
AcquireLock(true);
if((m_MarkerBuffOfs + (2 * cMaxExtKMerLen)) > m_AllocMarkerBuffSize)
	{
	CUtility::SafeWrite(m_hOutReadsFile,m_pMarkerBuff,m_MarkerBuffOfs);
	m_MarkerBuffOfs = 0;
	}
m_MarkerID += 1;
if(ReadLen != TruncMarkerLen)					// if needed to clamp
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"ReportMarkerRead: Truncated overlength (%d) Marker%d down to length %d",ReadLen,m_MarkerID,cMaxExtKMerLen);
	ReadLen = cMaxExtKMerLen;
	}
m_MarkerBuffOfs += sprintf((char *)&m_pMarkerBuff[m_MarkerBuffOfs],">%sS%d %d|%d|%d\n%s\n",m_szTargetCultivar,m_MarkerID,TruncMarkerLen,m_KMerLen,m_MinHamming,szMarkerSeq);
m_pMarkerLenDist[TruncMarkerLen-1] += 1;
if(m_MaxMarkerLen < TruncMarkerLen)
	m_MaxMarkerLen = TruncMarkerLen;
if(m_MinMarkerLen == 0 || m_MinMarkerLen > TruncMarkerLen)
	m_MinMarkerLen = TruncMarkerLen;
ReleaseLock(true);
return(eBSFSuccess);
}

int												// number of bases still unprocessed for K-mers in sequence buffer										
CLocKMers::FillBlockSeqBuff(void)				// maximally fill sequence buffer		
{
int CntNs;
int Idx;
int SeqLen;
int BlockSeqLen;
int NumUnprocessed;
tsPartialCultivar *pPartialCultivar;
UINT8 *pDstBase;
UINT8 *pSrcBase;

// only bother filling if sequence buffer getting lower than 3 * cBlockReqSize bases
if(!m_bAllEntrySeqsLoaded && (m_CurPartialCultivarID == 0 || m_CurPartialCultivarID > m_NumPartialCultivars))
	return(eBSFerrInternal);	

NumUnprocessed = (int)m_BlockSeqBuffLen - (int)m_NxtBlockSeqBuffOfs;
if(m_bAllEntrySeqsLoaded || NumUnprocessed >= (3 * cBlockReqSize))
	return(NumUnprocessed);

if(m_NxtBlockSeqBuffOfs > 0)
	{
	memmove(m_pBlockSeqBuff,&m_pBlockSeqBuff[m_NxtBlockSeqBuffOfs],NumUnprocessed);
	m_NxtBlockSeqBuffOfs = 0;
	}
m_BlockSeqBuffLen = NumUnprocessed;

	// try to maximally fill buffer
pPartialCultivar = &m_PartialCultivars[m_CurPartialCultivarID - 1];
if(m_CurPartialCultivarOfs >= pPartialCultivar->EntryLen)		// ensure haven't already completed buffering sequences from this partial cultivar...
	{
	m_CurPartialCultivarOfs = 0;
	if(m_CurPartialCultivarID == m_NumPartialCultivars)			// no more partial cultivars to load sequences from?
		{
		m_CurPartialCultivarID = 0;
		m_bAllEntrySeqsLoaded = true;
		return(m_BlockSeqBuffLen);
		}
	m_CurPartialCultivarID += 1;								// try next partial cultivar
	m_pBlockSeqBuff[m_NxtBlockSeqBuffOfs++] = eBaseEOS;			// ensure current buffered sequence is terminated
	}

do {
	pPartialCultivar = &m_PartialCultivars[m_CurPartialCultivarID - 1];
	SeqLen = m_pSfxArray->GetSeq(pPartialCultivar->EntryID,m_CurPartialCultivarOfs,&m_pBlockSeqBuff[m_BlockSeqBuffLen],(m_AllocBlockSeqBuffSize - 10) - m_BlockSeqBuffLen);
	if(SeqLen > 0)			// loaded any sequence?
		{
		m_CurPartialCultivarOfs += (UINT32)SeqLen;
		m_BlockSeqBuffLen += (UINT32)SeqLen;
		if(m_CurPartialCultivarOfs >= pPartialCultivar->EntryLen)		// need to next load from another partial cultivar?
			{
			m_pBlockSeqBuff[m_BlockSeqBuffLen++] = eBaseEOS;			// ensure current buffered sequence is terminated
			m_CurPartialCultivarOfs = 0;
			if(m_CurPartialCultivarID == m_NumPartialCultivars)			// no more partial cultivars to load sequences from?
				{
				m_CurPartialCultivarID = 0;
				m_bAllEntrySeqsLoaded = true;
				}
			else
				m_CurPartialCultivarID += 1;								// try next partial cultivar
			}
		}
	else
		return(eBSFerrInternal);
	}	
while(!m_bAllEntrySeqsLoaded && m_BlockSeqBuffLen < (m_AllocBlockSeqBuffSize - 10));

// any run of Ns or eBaseEOS's are assumed to be inter-sequence markers and are to be replaced by a single eBaseEOS
pDstBase = pSrcBase = m_pBlockSeqBuff;
BlockSeqLen = 0;
for(CntNs = 0, Idx = 0; Idx < (int)m_BlockSeqBuffLen; Idx++,pSrcBase++)
	{
	if(*pSrcBase == eBaseN || *pSrcBase == eBaseEOS)				
		{
		CntNs += 1;
		continue;
		}
	if(CntNs > 0)	
		{
		*pDstBase++ = eBaseEOS;
		BlockSeqLen += 1;
		CntNs = 0;
		}
	*pDstBase++ = *pSrcBase;
	BlockSeqLen += 1;
	}
if(CntNs)
	{
	*pDstBase++ = eBaseEOS;
	BlockSeqLen += 1;
	}

m_BlockSeqBuffLen = BlockSeqLen;
m_pBlockSeqBuff[m_BlockSeqBuffLen] = eBaseEOG;
m_NxtBlockSeqBuffOfs = 0;
return(m_BlockSeqBuffLen);
}

static UINT32 m_NumRetBlocks = 0;

int										    // returned block of concatenated sequences total length
CLocKMers::GetBlockSeqs(int MaxLength,		//  maximum total block length to return
						UINT8 *pBlockSeqs)	// copy block of sequences into this buffer
{
int MaxAvailBlockSize;
UINT8 *pSrcBase;
UINT8 *pDstBase;
UINT8 CurBase;
int CurRetLen;
int EOSRetLen;	

if(m_pBlockSeqBuff == NULL || m_NumPartialCultivars == 0 || m_NumSfxEntries == 0)
	return(eBSFerrInternal);

AcquireLock(true);
if(m_bInitSeqBuffering)						// just starting and buffering requires initialisation?
	{
	m_NumRetBlocks = 0;
	m_NxtBlockSeqBuffOfs = 0;
	m_BlockSeqBuffLen = 0;
	m_bAllEntrySeqsLoaded = false;
	m_bAllBlocksReturned = false;
	m_CurPartialCultivarID = 1;		// start with this partial cultivar
	m_CurPartialCultivarOfs = 0;
	FillBlockSeqBuff();
	m_bInitSeqBuffering = false;
	if(!m_BlockSeqBuffLen)						// should alway be able to load initial sequences...
		{
		m_CurPartialCultivarOfs = 0;
		m_NxtBlockSeqBuffOfs = 0;
		m_BlockSeqBuffLen = 0;
		m_bAllEntrySeqsLoaded = true;
		m_bAllBlocksReturned = true;
		ReleaseLock(true);
		return(0);
		}
	}
else
	{
	if(m_bAllBlocksReturned || (m_bAllEntrySeqsLoaded && m_NxtBlockSeqBuffOfs >= m_BlockSeqBuffLen))	// have all partial cultivar chrom sequence blocks been returned for processing?
		{
		m_NxtBlockSeqBuffOfs = 0;
		m_BlockSeqBuffLen = 0;
		m_bInitSeqBuffering = false;
		m_bAllEntrySeqsLoaded = true;
		m_bAllBlocksReturned = true;
		ReleaseLock(true);
		return(0);
		}
	}

// determine the maximum potential size of block which could be returned 
if(m_NxtBlockSeqBuffOfs < m_BlockSeqBuffLen)
	MaxAvailBlockSize = m_BlockSeqBuffLen - m_NxtBlockSeqBuffOfs;
else
	MaxAvailBlockSize = 0;

// is buffer getting low, if so then try to top up
if(!m_bAllEntrySeqsLoaded && MaxAvailBlockSize < (2 * MaxLength))
	{
	FillBlockSeqBuff();
	if(m_BlockSeqBuffLen == 0)					// sequence buffer still empty even after trying to top it up?
		{
		m_CurPartialCultivarOfs = 0;
		m_NxtBlockSeqBuffOfs = 0;
		m_BlockSeqBuffLen = 0;
		m_bAllEntrySeqsLoaded = true;
		m_bAllBlocksReturned = true;
		ReleaseLock(true);
		return(0);
		}
	}

// ensure returned blocks must contain complete rather than partial sequences
pSrcBase = &m_pBlockSeqBuff[m_NxtBlockSeqBuffOfs];
// slough any initial eBaseEOS's
while(m_NxtBlockSeqBuffOfs < m_BlockSeqBuffLen && (CurBase = *pSrcBase) == eBaseEOS)		
	{
	pSrcBase += 1;
	m_NxtBlockSeqBuffOfs += 1;			
	}

// check if removing initial eBaseEOS's has emptied buffer
if(*pSrcBase == eBaseEOG || m_NxtBlockSeqBuffOfs == m_BlockSeqBuffLen)
	{
	m_CurPartialCultivarOfs = 0;
	m_NxtBlockSeqBuffOfs = 0;
	m_BlockSeqBuffLen = 0;
	m_bAllEntrySeqsLoaded = true;
	m_bAllBlocksReturned = true;
	ReleaseLock(true);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"GetBlockSeqs: At block %d, removing eBaseEOSs emptied buffer",m_NumRetBlocks);
	return(0);
	}

// buffer not empty return upto MaxLength bases ensuring last returned sequence is not a partial sequence
pDstBase = pBlockSeqs;
CurRetLen = 0;
EOSRetLen = 0;
while((CurBase = *pSrcBase++) != eBaseEOG)			
	{
	CurRetLen += 1;
	*pDstBase++ = CurBase;
	if(CurBase == eBaseEOS)		
		EOSRetLen = CurRetLen;
	if(CurRetLen == MaxLength)
		break;
	}

m_NxtBlockSeqBuffOfs += EOSRetLen;
if(EOSRetLen == 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"GetBlockSeqs: At block %d, EOSRetLen:0 CurRetLen:%d MaxLength: %d",m_NumRetBlocks,CurRetLen,MaxLength);
else
	m_NumRetBlocks += 1;
ReleaseLock(true);
return(EOSRetLen);	
}
	

bool	// returns true if suffix array entry identifier is for a pseudo-chromosome which is part of the target cultivar  
CLocKMers::IsTargetCultivar(int EntryID)		// suffix array entry identifier 
{
return(LocTargetCultivar(EntryID) == NULL ? false : true);
}

tsPartialCultivar *				// NULL if EntryID not a targeted cultivar
CLocKMers::LocTargetCultivar(int EntryID)	// suffix array entry identifier
{
int Idx;
tsPartialCultivar *pPartialCultivar;
pPartialCultivar = m_PartialCultivars;
for(Idx = 0; Idx < m_NumPartialCultivars; Idx++, pPartialCultivar++)
	if(pPartialCultivar->EntryID == EntryID)
		return(pPartialCultivar);
return(NULL);
}

UINT32												// returns number of accepted cultivar specific extended K-mers
CLocKMers::GetKMerProcProgress(UINT32 *pNumBlocks,	// number of blocks processed
					INT64 *pNumKMers,				// these blocks contain this number of K-mers
					UINT32 *pNumPutativeKMers,		// putatively - before check for Hamming - there are this number of unique K-mers mapping to target cultivar
					UINT32 *pNumAcceptedKMers)		// after Hamming check this number of K-mers have been accepted
{
UINT32 NumAcceptedExtdKMers;
AcquireLock(true);
if(pNumBlocks != NULL)
	*pNumBlocks = m_NumBlocks;
if(pNumKMers != NULL)
	*pNumKMers = m_NumKMers;
if(pNumPutativeKMers != NULL)
	*pNumPutativeKMers = m_NumPutativeKMers;
if(pNumAcceptedKMers != NULL)
	*pNumAcceptedKMers = m_NumAcceptedKMers;
NumAcceptedExtdKMers = m_NumAcceptedExtdKMers;
ReleaseLock(true);
return(NumAcceptedExtdKMers);
}

int 
CLocKMers::MarkersCallback(void *pThis,UINT32 EntryID,etSeqBase *pSeq)
{

return(0);
}

int
CLocKMers::LocKMers(etPMode PMode,				// processing mode - defaults to 0
		  int KMerLen,					// this length K-mers
	  	  int PrefixLen,				// inter-cultivar shared prefix length
		  int SuffixLen,				// cultivar specific suffix length
		  int MinWithPrefix,			// minimum number of cultivars required to have the shared prefix
		  int MinHamming,				// must be at least this Hamming away from any other K-mer in other cultivars
		  char *pszCultivarName,		// targeted cultivar name for which K-mer markers are required
		  int NumPartialCultivars,		// there are this many pseudo chromosomes for targeted cultivar for which K-mer markers are required
		  char *ppszPartialCultivars[],	// pseudo chromosome names which identify targeted cultivar
		  char *pszSfxPseudoGenome,		// contains pregenerated suffix over psuedochromosomes for each cultivar
		  char *pszMarkerFile,			// output potential markers to this file
		  char *pszMarkerReadsFile,		// optionally output reads containing potential markers to this file
		  int NumThreads)				// max number of threads allowed
{
int Rslt;
int Idx;
UINT32 NumBlocks;
INT64 NumKMers;
UINT32 NumPutativeKMers;
UINT32 NumAcceptedKMers;
UINT32 NumAcceptedExtdKMers;
tsPartialCultivar *pPartialCultivar;

Reset();

m_PMode = PMode;
m_KMerLen = KMerLen;
m_PrefixLen = PrefixLen;
m_SuffixLen = SuffixLen;
m_MinWithPrefix = MinWithPrefix;

m_MinHamming = MinHamming;
m_NumThreads = NumThreads;
strncpy(m_szMarkerFile,pszMarkerFile,sizeof(m_szMarkerFile));
m_szMarkerFile[sizeof(m_szMarkerFile)-1] = '\0';
if(pszMarkerReadsFile != NULL && pszMarkerReadsFile[0] != '\0')
	{
	m_bKMerReads = true;
	strncpy(m_szMarkerReadsFile,pszMarkerReadsFile,sizeof(m_szMarkerReadsFile));
	m_szMarkerReadsFile[sizeof(m_szMarkerReadsFile)-1] = '\0';
	}
else
	{
	m_bKMerReads = false;
	m_szMarkerReadsFile[0] = '\0';
	}
strncpy(m_szTargetCultivar,pszCultivarName,sizeof(m_szTargetCultivar));
m_szTargetCultivar[sizeof(m_szTargetCultivar)-1] = '\0';
m_NumPartialCultivars = NumPartialCultivars;
pPartialCultivar = m_PartialCultivars;
for(Idx = 0; Idx < NumPartialCultivars; Idx++,pPartialCultivar++)
	{
	pPartialCultivar->PartialCultivarID = 0;
	pPartialCultivar->EntryID = 0;
	pPartialCultivar->EntryLen = 0;
	strncpy(pPartialCultivar->szEntryName,ppszPartialCultivars[Idx],sizeof(pPartialCultivar->szEntryName));
	pPartialCultivar->szEntryName[sizeof(pPartialCultivar->szEntryName)-1] = '\0';
	}

// sort so can check for duplicate cultivar pseudo-chromosome
if(NumPartialCultivars > 1)
	{
	qsort(m_PartialCultivars,NumPartialCultivars,sizeof(tsPartialCultivar),sortpartcultnames);
	pPartialCultivar = m_PartialCultivars;
	for(Idx = 0; Idx < NumPartialCultivars-1; Idx++,pPartialCultivar++)
		if(!stricmp(pPartialCultivar->szEntryName,pPartialCultivar[1].szEntryName))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Duplicate partial cultivar pseudo-chromosome name : '%s'",pPartialCultivar->szEntryName);
			return(eBSFerrParams);
			}
	}

// allocate buffers
if((m_pMarkerBuff = new UINT8 [cMarkerSeqBuffSize])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Unable to allocate (%d bytes) for marker buffering",cMarkerSeqBuffSize);
	Reset();
	return(eBSFerrMem);
	}
m_AllocMarkerBuffSize = cMarkerSeqBuffSize;

if((m_pMarkerLenDist = new UINT32 [cMaxExtKMerLen])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Unable to allocate (%d bytes) for marker length distributions",sizeof(UINT32) * cMaxExtKMerLen);
	Reset();
	return(eBSFerrMem);
	}


if((m_pBlockSeqBuff = new UINT8 [cConcatSeqBuffSize])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Unable to allocate (%d bytes) for sequence buffering",cConcatSeqBuffSize);
	Reset();
	return(eBSFerrMem);
	}
m_AllocBlockSeqBuffSize = cConcatSeqBuffSize;

if((m_pSfxArray = new CSfxArrayV3)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Error: Unable to instantiate instance of CSfxArrayV3");
	Reset();
	return(eBSFerrObj);
	}

if((Rslt = m_pSfxArray->Open(pszSfxPseudoGenome))!=eBSFSuccess)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open input bioseq suffix array file '%s'",pszSfxPseudoGenome);
	Reset();
	return(Rslt);
	}

// report to user some sfx array metadata for user conformation the targeted assembly is correct
strcpy(m_szDataset,m_pSfxArray->GetDatasetName());
if(m_pSfxArray->IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process SOLiD colorspace suffix array file '%s'",pszSfxPseudoGenome);
	Reset();
	return(eBSFerrRefDataset);
	}
tsSfxHeaderV3 SfxHeader;
m_pSfxArray->GetSfxHeader(&SfxHeader);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Psuedo-assembly Name: '%s' Descr: '%s' Title: '%s' Version: %d",
					 m_szDataset,SfxHeader.szDescription,SfxHeader.szTitle,SfxHeader.Version);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assembly has blocks: %d, max block size: %llu",SfxHeader.NumSfxBlocks,SfxHeader.SfxBlockSize);

// ensure assembly and suffix array has been fully loaded into memory 
int CurBlockID = 1;		// currently only single block supported
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading genome assembly suffix array...");
if((Rslt=m_pSfxArray->SetTargBlock(CurBlockID))<0)
	{
	while(m_pSfxArray->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pSfxArray->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fatal: unable to load assembly suffix array");
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Assembly suffix array loaded");


// ensure that all the targeted pseudo-chromosomes are in the suffix array, and that there are other chromosomes (other species) in addition to those targeted
// at same time the pseudo-chromosome ID's and lengths are retreived
m_NumSfxEntries = m_pSfxArray->GetNumEntries();
if(m_NumSfxEntries < 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No entries in '%s'",pszSfxPseudoGenome);
	Reset();
	return(eBSFerrEntry);
	}
if(m_NumSfxEntries <= m_NumPartialCultivars)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Only %d entries in '%s', expected additional ( > %d) to targeted cultivar partial pseudo-chromosomes ",m_NumSfxEntries,pszSfxPseudoGenome,m_NumPartialCultivars);
	Reset();
	return(eBSFerrEntry);
	}
pPartialCultivar = m_PartialCultivars;
for(Idx = 0; Idx < NumPartialCultivars; Idx++,pPartialCultivar++)
	{
	if((pPartialCultivar->EntryID = m_pSfxArray->GetIdent(pPartialCultivar->szEntryName)) < 1)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to locate '%s' in '%s'",pPartialCultivar->szEntryName,pszSfxPseudoGenome);
		Reset();
		return(eBSFerrEntry);
		}
	pPartialCultivar->PartialCultivarID = Idx+1;
	m_CultChromIDs[Idx] = pPartialCultivar->EntryID;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Located target cultivar pseudo-chrom '%s' in '%s'",pPartialCultivar->szEntryName,pszSfxPseudoGenome);
	pPartialCultivar->EntryLen = m_pSfxArray->GetSeqLen(pPartialCultivar->EntryID);
	}

// now report the non-targeted chrom names whilst populating m_AllCultivars[]
int EntryID;
char szSfxEntryName[100];
tsCultivar *pCultivar;
tsPartialCultivar *pPartial;
pCultivar = m_AllCultivars;

for(EntryID = 1; EntryID <= m_NumSfxEntries; EntryID++, pCultivar++)
	{
	pCultivar->Status = 0;
	m_pSfxArray->GetIdentName(EntryID,sizeof(szSfxEntryName),szSfxEntryName);
	pCultivar->Cultivar.EntryID = EntryID;
	strncpy(pCultivar->Cultivar.szEntryName,szSfxEntryName,sizeof(pCultivar->Cultivar.szEntryName));
	pCultivar->Cultivar.szEntryName[sizeof(pCultivar->Cultivar.szEntryName)-1] = '\0';
	pCultivar->Cultivar.EntryLen = m_pSfxArray->GetSeqLen(EntryID);
	if((pPartial = LocTargetCultivar(EntryID))!=NULL)
		{
		pCultivar->Cultivar.PartialCultivarID = pPartial->PartialCultivarID;
		continue;
		}
	pCultivar->Cultivar.PartialCultivarID = 0;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"   Processing K-mers against non-target pseudo-chromosome '%s'",szSfxEntryName);
	}

if(PMode == ePMPrefixKMers && (m_MinWithPrefix == 0 || m_MinWithPrefix > (m_NumSfxEntries - NumPartialCultivars)))
	m_MinWithPrefix = m_NumSfxEntries - NumPartialCultivars;

// looks good to go so create/truncate output marker sequence file
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Truncating/Creating output multi-fasta marker file '%s'",m_szMarkerFile);
#ifdef _WIN32
m_hOutFile = open(m_szMarkerFile,O_CREATETRUNC );
#else
if((m_hOutFile = open(m_szMarkerFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_szMarkerFile,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",m_szMarkerFile);
	m_hOutFile = -1;
	Reset();
	return(eBSFerrCreateFile);
	}


// user also wanting K-mer marker containing reads to be reported?
if(m_bKMerReads)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Truncating/Creating output marker containing multi-fasta reads file '%s'",m_szMarkerReadsFile);
	#ifdef _WIN32
	m_hOutReadsFile = open(m_szMarkerReadsFile,O_CREATETRUNC );
	#else
	if((m_hOutReadsFile = open(m_szMarkerReadsFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hOutReadsFile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_szMarkerReadsFile,strerror(errno));
				Reset();
				return(eBSFerrCreateFile);
				}
	#endif

	if(m_hOutReadsFile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output marker containing multi-fasta reads file '%s'",m_szMarkerReadsFile);
		m_hOutReadsFile = -1;
		Reset();
		return(eBSFerrCreateFile);
		}
	}

memset(m_pMarkerLenDist,0,sizeof(UINT32) * cMaxExtKMerLen);
m_bInitSeqBuffering = true;

// TIME TO PLAY!!!!
#ifdef TRYCULTIVARMARKERS
UINT32 EntryIDs[50];
int EntryIDIdx;
pCultivar = m_AllCultivars;
for(EntryIDIdx = 0; EntryIDIdx < m_NumSfxEntries; EntryIDIdx++, pCultivar++)
	{
	pCultivar->Status = 0;
	EntryIDs[EntryIDIdx] = pCultivar->Cultivar.EntryID;
	}
m_pSfxArray->LocateMultiCultivarMarkers(30,30,m_NumSfxEntries,EntryIDs,this,MarkersCallback);

// END OF TIME TO PLAY
#endif

tsKMerThreadPars WorkerThreads[cMaxWorkerThreads];			// allow for max possible user configured number of threads
int ThreadIdx;

// initialise and startup K-mer processing worker threads
memset(WorkerThreads,0,sizeof(WorkerThreads));
for(ThreadIdx = 0; ThreadIdx < m_NumThreads; ThreadIdx++)
	{
	WorkerThreads[ThreadIdx].ThreadIdx = ThreadIdx + 1;
	WorkerThreads[ThreadIdx].pThis = this;
	WorkerThreads[ThreadIdx].pBlockSeqs = new UINT8 [cBlockReqSize];
	WorkerThreads[ThreadIdx].AllocBlockSeqsSize = cBlockReqSize;
#ifdef _WIN32
	WorkerThreads[ThreadIdx].threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,KMerThreadStart,&WorkerThreads[ThreadIdx],0,&WorkerThreads[ThreadIdx].threadID);
#else
	WorkerThreads[ThreadIdx].threadRslt =	pthread_create (&WorkerThreads[ThreadIdx].threadID , NULL , KMerThreadStart , &WorkerThreads[ThreadIdx] );
#endif
	}

// allow threads a few seconds to startup
#ifdef _WIN32
	Sleep(10000);
#else
	sleep(10);
#endif

// let user know that this K-mer processing process is working hard...
NumAcceptedExtdKMers = GetKMerProcProgress(&NumBlocks,&NumKMers,&NumPutativeKMers, &NumAcceptedKMers);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress - K-Mers processed: %lld, Cultivar specific : %u, Hamming retained: %u, Accepted: %u",NumKMers,NumPutativeKMers, NumAcceptedKMers,NumAcceptedExtdKMers);

// wait for all threads to have completed
for(ThreadIdx = 0; ThreadIdx < m_NumThreads; ThreadIdx++)
	{
#ifdef _WIN32
	while(WAIT_TIMEOUT == WaitForSingleObject( WorkerThreads[ThreadIdx].threadHandle, 60000 * 10))
		{
		NumAcceptedExtdKMers = GetKMerProcProgress(&NumBlocks,&NumKMers,&NumPutativeKMers, &NumAcceptedKMers);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress - K-Mers processed: %lld, Cultivar specific : %u, Hamming retained: %u, Accepted: %u",NumKMers,NumPutativeKMers, NumAcceptedKMers,NumAcceptedExtdKMers);
		}
	CloseHandle( WorkerThreads[ThreadIdx].threadHandle);
#else
	struct timespec ts;
	int JoinRlt;
	clock_gettime(CLOCK_REALTIME, &ts);
	ts.tv_sec += 60 * 10;
	while((JoinRlt = pthread_timedjoin_np(WorkerThreads[ThreadIdx].threadID, NULL, &ts)) != 0)
		{
		NumAcceptedExtdKMers = GetKMerProcProgress(&NumBlocks,&NumKMers,&NumPutativeKMers, &NumAcceptedKMers);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress - K-Mers processed: %lld, Cultivar specific : %u, Hamming retained: %u, Accepted: %u",NumKMers,NumPutativeKMers, NumAcceptedKMers,NumAcceptedExtdKMers);
		ts.tv_sec += 60;
		}
#endif
	if(WorkerThreads[ThreadIdx].pBlockSeqs!=NULL)
		{
		delete WorkerThreads[ThreadIdx].pBlockSeqs;
		WorkerThreads[ThreadIdx].pBlockSeqs = NULL;
		}
	}

NumAcceptedExtdKMers = GetKMerProcProgress(&NumBlocks,&NumKMers,&NumPutativeKMers, &NumAcceptedKMers);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed - K-Mers processed: %lld, Cultivar specific : %u, Hamming retained: %u, Accepted: %u",NumKMers,NumPutativeKMers, NumAcceptedKMers,NumAcceptedExtdKMers);

if(gProcessingID > 0)
	{
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Kmers",ePTInt64,sizeof(NumKMers),"Processed",&NumKMers);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Kmers",ePTUint32,sizeof(NumPutativeKMers),"Cultivar",&NumPutativeKMers);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Kmers",ePTUint32,sizeof(NumAcceptedKMers),"HammingRetained",&NumAcceptedKMers);
	gSQLiteSummaries.AddResult(gProcessingID,(char *)"Kmers",ePTUint32,sizeof(NumAcceptedExtdKMers),"Accepted",&NumAcceptedExtdKMers);
	}

if(m_hOutFile != -1)
	{
	if(m_MarkerBuffOfs)
		CUtility::SafeWrite(m_hOutFile,m_pMarkerBuff,m_MarkerBuffOfs);
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Closed output multi-fasta marker file '%s'",m_szMarkerFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Number of reported markers: %u",m_MarkerID);
if(m_MarkerID)
	{
	// log the marker length distributions
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Marker length distribution");
	if(m_MinMarkerLen)
		for(Idx = m_MinMarkerLen; Idx <= m_MaxMarkerLen; Idx++)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Length: %d Count: %u",Idx,m_pMarkerLenDist[Idx-1]);
			if(gProcessingID)
				gSQLiteSummaries.AddResultXY(gProcessingID,(char *)"Markers",ePTInt32,sizeof(Idx),"Length",&Idx,ePTInt32,sizeof(m_pMarkerLenDist[Idx-1]),"Cnt",&m_pMarkerLenDist[Idx-1]);
			}
	else
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Length: %d Count: 0",m_KMerLen);
	}

// user requested marker containing reads also be reported?
if(m_bKMerReads)
	{
	bool bMarkerRead;
	int CultID;
	tsPartialCultivar *pCult;
	UINT8 *pSeq;
	UINT8 *pMarkStart;
	UINT8 SeqBase;
	UINT32 Loci;
	UINT32 ReadLen;

	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Writing reads containing markers to '%s'",m_szMarkerReadsFile);
	m_MarkerID = 0;
	m_MarkerBuffOfs = 0;
	m_MinMarkerLen = 0;
	m_MaxMarkerLen = 0;

	memset(m_pMarkerLenDist,0,sizeof(UINT32) * cMaxExtKMerLen);
	pCult = m_PartialCultivars;
	for(CultID = 1; CultID <= m_NumPartialCultivars; CultID++,pCult++)
		{
		bMarkerRead = false;
		ReadLen = 0;
		pSeq = (UINT8 *)m_pSfxArray->GetPtrSeq(pCult->EntryID,0);
		pMarkStart = pSeq;
		for(Loci = 0; Loci < pCult->EntryLen; Loci++)
			{
			SeqBase = *pSeq++;
			if((SeqBase & 0x0f) > eBaseT || Loci == pCult->EntryLen-1)
				{
				if(bMarkerRead)	
					{
					// write to file this read which starts at pMarkStart and is of length ReadLen
					ReportMarkerRead(ReadLen,pMarkStart);
					bMarkerRead = false;
					}
				pMarkStart = pSeq;
				ReadLen = 0;
				continue;
				}
			ReadLen += 1;
			if((SeqBase >> 4) & 0x02)	// marked as containing K-mer marker?
				bMarkerRead = true;
			}
		}

	if(m_hOutReadsFile != -1)
		{
		if(m_MarkerBuffOfs)
			CUtility::SafeWrite(m_hOutReadsFile,m_pMarkerBuff,m_MarkerBuffOfs);
#ifdef _WIN32
		_commit(m_hOutReadsFile);
#else
		fsync(m_hOutReadsFile);
#endif
		close(m_hOutReadsFile);
		m_hOutReadsFile = -1;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Closed output reads containing marker multi-fasta reads file '%s'",m_szMarkerReadsFile);
		}

		// log the marker reads length distributions
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Number of reported reads containing markers: %u",m_MarkerID);
	if(m_MarkerID)
		{
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Marker containing reads length distribution");
		if(m_MinMarkerLen)
			for(Idx = m_MinMarkerLen; Idx <= m_MaxMarkerLen; Idx++)
				{
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Length: %d Count: %u",Idx,m_pMarkerLenDist[Idx-1]);
				if(gProcessingID)
					gSQLiteSummaries.AddResultXY(gProcessingID,(char *)"MarkerReads",ePTInt32,sizeof(Idx),"Length",&Idx,ePTInt32,sizeof(m_pMarkerLenDist[Idx-1]),"Cnt",&m_pMarkerLenDist[Idx-1]);
				}
		else
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Length: %d Count: 0",m_KMerLen);
		}
	}


Reset();
return(eBSFSuccess);
}

// Thread startup
#ifdef WIN32
unsigned int __stdcall CLocKMers::KMerThreadStart(void *args)
{
#else
void * CLocKMers::KMerThreadStart(void *args)
{
#endif
tsKMerThreadPars *pArgs = (tsKMerThreadPars *)args;
pArgs->pThis->LocateSpeciesUniqueKMers(pArgs);
#ifdef WIN32
ExitThread(1);
#else
return NULL;
#endif
}

int							// marking reads containing the identified marker, both sense and antisense reads marked
CLocKMers::MarkContainingReads(int MarkKMerLen,						// marker is of this length
							   etSeqBase *pMarkStartKMerBase)		// marker sequence	
{
UINT8 RevCplKMerSeq[cMaxKMerLen+1];
INT64 PrevHitIdx;
INT64 NxtHitIdx;
UINT32 TargHitID;
UINT32 TargHitLoci;

memmove(RevCplKMerSeq,pMarkStartKMerBase,MarkKMerLen);
CSeqTrans::ReverseComplement(MarkKMerLen,RevCplKMerSeq);

PrevHitIdx = 0;
while((NxtHitIdx = m_pSfxArray->IterateExacts(pMarkStartKMerBase,MarkKMerLen,PrevHitIdx,&TargHitID,&TargHitLoci)) > 0)
	{
	m_pSfxArray->SetBaseFlags(TargHitID,TargHitLoci,0x02);
	PrevHitIdx = NxtHitIdx;
	}

PrevHitIdx = 0;
pMarkStartKMerBase = RevCplKMerSeq;
while((NxtHitIdx = m_pSfxArray->IterateExacts(pMarkStartKMerBase,MarkKMerLen,PrevHitIdx,&TargHitID,&TargHitLoci)) > 0)
	{
	m_pSfxArray->SetBaseFlags(TargHitID,TargHitLoci,0x02);
	PrevHitIdx = NxtHitIdx;
	}

return(eBSFSuccess);	
}


int
CLocKMers::LocateSpeciesUniqueKMers(tsKMerThreadPars *pPars)
{
UINT8 RevCplKMer[cMaxKMerLen+1];			// to hold reverse compliment of current K-mer being processed
UINT8 PrefixedCultivars[cMaxTargCultivarChroms];
int SubSeqsLen;
int StartIdx;
int	MarkStartIdx;
etSeqBase *pMarkStartKMerBase;
int	MarkKMerLen;
int CurSeqLen;
INT64 PrevHitIdx;
INT64 NxtHitIdx;
INT64 RcplHitIdx;
UINT32 TargHitID;
UINT32 TargHitLoci;
UINT32 RcplTargHitID;
UINT32 RcplTargHitLoci;
UINT32 NumBlocks;
INT64 NumKMers;
UINT32 NumPutativeKMers;
UINT32 NumAcceptedKMers;
UINT32 NumAcceptedExtdKMers;
int MatchesOthers;

int NumCultsPrefixed;

etSeqBase *pStartKMerBase;
etSeqBase *pEndKMerBase;

bool bTargHit;
bool bNonTargHit;
bool bIsDup;
int LociFlags;

int Cntr = 0;
NumPutativeKMers = 0;
NumAcceptedKMers = 0;
NumAcceptedExtdKMers = 0;
NumBlocks = 0;
NumKMers = 0;
while((SubSeqsLen=GetBlockSeqs(pPars->AllocBlockSeqsSize - 1,pPars->pBlockSeqs)) > 0)  // -1 because a eBaseEOG will be later appended
	{
	if(NumBlocks++ > 3)			// could report more frequently but as long as reporting a few times a minute then that's enough
		{
		AcquireLock(true);
		m_NumBlocks += NumBlocks;
		m_NumKMers += NumKMers;
		m_NumPutativeKMers += NumPutativeKMers;
		m_NumAcceptedKMers += NumAcceptedKMers;
		m_NumAcceptedExtdKMers += NumAcceptedExtdKMers;
		ReleaseLock(true);
		NumBlocks = 0;
		NumKMers = 0;
		NumPutativeKMers = 0;
		NumAcceptedKMers = 0;
		NumAcceptedExtdKMers = 0;
		}

	if(SubSeqsLen <= m_KMerLen)    // blocks returned will always be terminated by eBaseEOS hence <= compare
		continue;

	pPars->pBlockSeqs[SubSeqsLen] = eBaseEOG;

		// processing here on subsequences
	pEndKMerBase = pPars->pBlockSeqs;
	do {
		pStartKMerBase = pEndKMerBase;
		CurSeqLen = 0;
		while(*pEndKMerBase++ <= eBaseT)		// ensure K-mers do not contain any eBaseN's (should have already been removed) or eBaseEOS's
			CurSeqLen += 1;
		if(CurSeqLen < m_KMerLen)
			continue;

		MarkKMerLen = 0;
		for(StartIdx = 0; StartIdx <= CurSeqLen - m_KMerLen; StartIdx++, pStartKMerBase++)
			{
			NumKMers += 1;
			bTargHit = false;
			bNonTargHit = false;
			bIsDup = false;
			PrevHitIdx = 0;
			NxtHitIdx = 0;
			RcplHitIdx = 0;
			memmove(RevCplKMer,pStartKMerBase,m_KMerLen);
			CSeqTrans::ReverseComplement(m_KMerLen,RevCplKMer);
			while((NxtHitIdx = m_pSfxArray->IterateExacts(pStartKMerBase,m_KMerLen,PrevHitIdx,&TargHitID,&TargHitLoci)) > 0)
				{
				if(!PrevHitIdx && IsTargetCultivar(TargHitID))
					{
					RcplHitIdx = m_pSfxArray->IterateExacts(RevCplKMer,m_KMerLen,0,&RcplTargHitID,&RcplTargHitLoci);
					if(RcplHitIdx != 0)
						LociFlags = m_pSfxArray->SetBaseFlags(TargHitID,TargHitLoci,RcplTargHitID,RcplTargHitLoci,0x01);
					else
						LociFlags = m_pSfxArray->SetBaseFlags(TargHitID,TargHitLoci,0x01);
					// check if already known that any previously processed K-mer hits this loci; if so this K-mer is a duplicate or revcpl of already processed K-mer
					if(LociFlags & 0x011)
						{
						bIsDup = true;			// treat as dup even though may have been revcpl
						break;
						}
					}

				if(IsTargetCultivar(TargHitID))			// was hit onto any of the target cultivar pseudo-chroms?
					bTargHit = true;
				else
					{
					bNonTargHit = true;
					break;
					}
				PrevHitIdx = NxtHitIdx;
				}
			if(bNonTargHit || bIsDup)
				continue;

			PrevHitIdx = RcplHitIdx;
			while((NxtHitIdx = m_pSfxArray->IterateExacts(RevCplKMer,m_KMerLen,PrevHitIdx,&TargHitID,&TargHitLoci)) > 0)
				{
				if(IsTargetCultivar(TargHitID))
					bTargHit = true;
				else
					{
					bNonTargHit = true;
					break;
					}
				PrevHitIdx = NxtHitIdx;
				}
			if(bNonTargHit || !bTargHit)
				continue;

			// this is a unique K-mer which only exists in the target species
			NumPutativeKMers += 1;

			// check if this K-mer is at least the minimum requested Hammings away from any K-mer in any of the other species
			MatchesOthers = 0;
			if(m_MinHamming > 1)
				{
				MatchesOthers = m_pSfxArray->MatchesOtherChroms(m_NumPartialCultivars, m_CultChromIDs, m_MinHamming - 1, m_KMerLen,pStartKMerBase);
				if(!MatchesOthers)
					MatchesOthers = m_pSfxArray->MatchesOtherChroms(m_NumPartialCultivars, m_CultChromIDs, m_MinHamming - 1, m_KMerLen,  RevCplKMer);
				}

			// if unable to match allowing subs in other cultivars then this K-mer is putatively accepted
			// if must share prefix sequence with other cultivars then check this condition can be met
			if(!MatchesOthers)
				{
				NumAcceptedKMers += 1;											// accepted on Hammings...

				if(m_PMode == ePMPrefixKMers && m_PrefixLen)					// if must share prefix sequence with other cultivars then check this condition can be met
					{
					PrevHitIdx = 0;
					NumCultsPrefixed = 0;
					memset(PrefixedCultivars,0,sizeof(PrefixedCultivars[0]) * m_NumSfxEntries); 
					while(NumCultsPrefixed < m_MinWithPrefix && ((NxtHitIdx = m_pSfxArray->IterateExacts(pStartKMerBase,m_PrefixLen,PrevHitIdx,&TargHitID,&TargHitLoci)) > 0))
						{
						PrevHitIdx = NxtHitIdx;
						if(IsTargetCultivar(TargHitID))				// not interested in self hits onto targeted cultivars
							continue;
						// is this a new cultivar not previously seen..
						if(!PrefixedCultivars[TargHitID-1])
							{
							PrefixedCultivars[TargHitID-1] = 1;
							NumCultsPrefixed += 1;
							}
						}
					if(NumCultsPrefixed < m_MinWithPrefix)
						{
						// retry with reverse complemented prefix sequence
						PrevHitIdx = 0;
						while(NumCultsPrefixed < m_MinWithPrefix && ((NxtHitIdx = m_pSfxArray->IterateExacts(RevCplKMer,m_PrefixLen,PrevHitIdx,&TargHitID,&TargHitLoci)) > 0))
							{
							PrevHitIdx = NxtHitIdx;
							if(IsTargetCultivar(TargHitID))				// not interested in self hits onto targeted cultivars
								continue;
							// is this a new cultivar not previously seen..
							if(!PrefixedCultivars[TargHitID-1])
								{
								PrefixedCultivars[TargHitID-1] = 1;
								NumCultsPrefixed += 1;
								}
							}
						if(NumCultsPrefixed < m_MinWithPrefix)
							continue;
						}
					}

				if(m_bKMerReads)
					MarkContainingReads(m_KMerLen,pStartKMerBase);

				// report K-mers of m_KMerLen only?
				if(m_PMode != ePMExtdKMers)
					{
					ReportMarker(m_KMerLen,pStartKMerBase);
					NumAcceptedExtdKMers+=1;
					continue;
					}

				// if would be a an extension of previous accepted then extend
				if((MarkKMerLen && (StartIdx > (MarkStartIdx + 1))))
					{
					ReportMarker(MarkKMerLen,pMarkStartKMerBase);
					NumAcceptedExtdKMers += 1;
					MarkKMerLen = 0;
					}
				if(MarkKMerLen == 0)
					{
					MarkKMerLen = m_KMerLen;
					pMarkStartKMerBase = pStartKMerBase;
					}
				else
					MarkKMerLen += 1;
				MarkStartIdx = StartIdx;
				}
			}
		if(m_PMode == ePMExtdKMers && MarkKMerLen)
			{
			ReportMarker(MarkKMerLen,pMarkStartKMerBase);
			NumAcceptedExtdKMers += 1;
			MarkKMerLen = 0;
			}
		}
	while(*pEndKMerBase != eBaseEOG);
	}
		
AcquireLock(true);
m_NumBlocks += NumBlocks;
m_NumKMers += NumKMers;
m_NumPutativeKMers += NumPutativeKMers;
m_NumAcceptedKMers += NumAcceptedKMers;
m_NumAcceptedExtdKMers += NumAcceptedExtdKMers;
ReleaseLock(true);
return(eBSFSuccess);
}


// LocateSharedUniqueKMers
// locate all unique K-mers of specified length which are common to all cultivars
int
CLocKMers::LocateSharedUniqueKMers(tsKMerThreadPars *pPars)
{
UINT8 RevCplKMer[cMaxKMerLen+1];			// to hold reverse compliment of current K-mer being processed

int SubSeqsLen;
int StartIdx;
int	MarkStartIdx;
etSeqBase *pMarkStartKMerBase;
int	MarkKMerLen;
int CurSeqLen;
INT64 PrevHitIdx;
INT64 NxtHitIdx;
INT64 RcplHitIdx;
UINT32 TargHitID;
UINT32 TargHitLoci;
UINT32 RcplTargHitID;
UINT32 RcplTargHitLoci;
UINT32 NumBlocks;
INT64 NumKMers;
UINT32 NumPutativeKMers;
UINT32 NumAcceptedKMers;
UINT32 NumAcceptedExtdKMers;
int MatchesOthers;

etSeqBase *pStartKMerBase;
etSeqBase *pEndKMerBase;

bool bTargHit;
bool bNonTargHit;
bool bIsDup;
int LociFlags;

int Cntr = 0;
NumPutativeKMers = 0;
NumAcceptedKMers = 0;
NumAcceptedExtdKMers = 0;
NumBlocks = 0;
NumKMers = 0;
while((SubSeqsLen=GetBlockSeqs(pPars->AllocBlockSeqsSize - 1,pPars->pBlockSeqs)) > 0)  // -1 because a eBaseEOG will be later appended
	{
	if(NumBlocks++ > 3)			// could report more frequently but as long as reporting a few times a minute then that's enough
		{
		AcquireLock(true);
		m_NumBlocks += NumBlocks;
		m_NumKMers += NumKMers;
		m_NumPutativeKMers += NumPutativeKMers;
		m_NumAcceptedKMers += NumAcceptedKMers;
		m_NumAcceptedExtdKMers += NumAcceptedExtdKMers;
		ReleaseLock(true);
		NumBlocks = 0;
		NumKMers = 0;
		NumPutativeKMers = 0;
		NumAcceptedKMers = 0;
		NumAcceptedExtdKMers = 0;
		}

	if(SubSeqsLen <= m_KMerLen)    // blocks returned will always be terminated by eBaseEOS hence <= compare
		continue;

	pPars->pBlockSeqs[SubSeqsLen] = eBaseEOG;

		// processing here on subsequences
	pEndKMerBase = pPars->pBlockSeqs;
	do {
		pStartKMerBase = pEndKMerBase;
		CurSeqLen = 0;
		while(*pEndKMerBase++ <= eBaseT)		// ensure K-mers do not contain any eBaseN's (should have already been removed) or eBaseEOS's
			CurSeqLen += 1;
		if(CurSeqLen < m_KMerLen)
			continue;

		MarkKMerLen = 0;
		for(StartIdx = 0; StartIdx <= CurSeqLen - m_KMerLen; StartIdx++, pStartKMerBase++)
			{
			NumKMers += 1;
			bTargHit = false;
			bNonTargHit = false;
			bIsDup = false;
			PrevHitIdx = 0;
			NxtHitIdx = 0;
			RcplHitIdx = 0;
			memmove(RevCplKMer,pStartKMerBase,m_KMerLen);
			CSeqTrans::ReverseComplement(m_KMerLen,RevCplKMer);
			while((NxtHitIdx = m_pSfxArray->IterateExacts(pStartKMerBase,m_KMerLen,PrevHitIdx,&TargHitID,&TargHitLoci)) > 0)
				{
				if(!PrevHitIdx && IsTargetCultivar(TargHitID))
					{
					RcplHitIdx = m_pSfxArray->IterateExacts(RevCplKMer,m_KMerLen,0,&RcplTargHitID,&RcplTargHitLoci);
					if(RcplHitIdx != 0)
						LociFlags = m_pSfxArray->SetBaseFlags(TargHitID,TargHitLoci,RcplTargHitID,RcplTargHitLoci,0x01);
					else
						LociFlags = m_pSfxArray->SetBaseFlags(TargHitID,TargHitLoci,0x01);
					// check if already known that any previously processed K-mer hits this loci; if so this K-mer is a duplicate or revcpl of already processed K-mer
					if(LociFlags & 0x011)
						{
						bIsDup = true;			// treat as dup even though may have been revcpl
						break;
						}
					}

				if(IsTargetCultivar(TargHitID))			// was hit onto any of the target cultivar pseudo-chroms?
					bTargHit = true;
				else
					{
					bNonTargHit = true;
					break;
					}
				PrevHitIdx = NxtHitIdx;
				}
			if(bNonTargHit || bIsDup)
				continue;

			PrevHitIdx = RcplHitIdx;
			while((NxtHitIdx = m_pSfxArray->IterateExacts(RevCplKMer,m_KMerLen,PrevHitIdx,&TargHitID,&TargHitLoci)) > 0)
				{
				if(IsTargetCultivar(TargHitID))
					bTargHit = true;
				else
					{
					bNonTargHit = true;
					break;
					}
				PrevHitIdx = NxtHitIdx;
				}
			if(bNonTargHit || !bTargHit)
				continue;

			// this is a unique K-mer which only exists in the target species
			NumPutativeKMers += 1;

			// check if this K-mer is at least the minimum requested Hammings away from any K-mer in any of the other species
			MatchesOthers = 0;
			if(m_MinHamming > 1)
				{
				MatchesOthers = m_pSfxArray->MatchesOtherChroms(m_NumPartialCultivars, m_CultChromIDs, m_MinHamming - 1, m_KMerLen,pStartKMerBase);
				if(!MatchesOthers)
					MatchesOthers = m_pSfxArray->MatchesOtherChroms(m_NumPartialCultivars, m_CultChromIDs, m_MinHamming - 1, m_KMerLen,  RevCplKMer);
				}

			// if unable to match allowing subs in other cultivars then this K-mer is accepted
			if(!MatchesOthers)
				{
				NumAcceptedKMers += 1;

				if(m_bKMerReads)
					MarkContainingReads(m_KMerLen,pStartKMerBase);

				// report K-mers of m_KMerLen only?
				if(m_PMode == ePMNoExtdKMers)
					{
					ReportMarker(m_KMerLen,pStartKMerBase);
					NumAcceptedExtdKMers+=1;
					continue;
					}

				// if would be a an extension of previous accepted then extend
				if((MarkKMerLen && (StartIdx > (MarkStartIdx + 1))))
					{
					ReportMarker(MarkKMerLen,pMarkStartKMerBase);
					NumAcceptedExtdKMers += 1;
					MarkKMerLen = 0;
					}
				if(MarkKMerLen == 0)
					{
					MarkKMerLen = m_KMerLen;
					pMarkStartKMerBase = pStartKMerBase;
					}
				else
					MarkKMerLen += 1;
				MarkStartIdx = StartIdx;
				}
			}
		if(m_PMode == ePMExtdKMers && MarkKMerLen)
			{
			ReportMarker(MarkKMerLen,pMarkStartKMerBase);
			NumAcceptedExtdKMers += 1;
			MarkKMerLen = 0;
			}
		}
	while(*pEndKMerBase != eBaseEOG);
	}
		
AcquireLock(true);
m_NumBlocks += NumBlocks;
m_NumKMers += NumKMers;
m_NumPutativeKMers += NumPutativeKMers;
m_NumAcceptedKMers += NumAcceptedKMers;
m_NumAcceptedExtdKMers += NumAcceptedExtdKMers;
ReleaseLock(true);
return(eBSFSuccess);
}

// used when sorting cultivar partial pseudo-chrom names
int 
CLocKMers::sortpartcultnames(const void *pEl1,const void *pEl2)
{
tsPartialCultivar *p1 = (tsPartialCultivar *)pEl1;
tsPartialCultivar *p2 = (tsPartialCultivar *)pEl2;
return(stricmp(p1->szEntryName,p2->szEntryName));
}
