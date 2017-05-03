/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif



CSfxNMerFreq::CSfxNMerFreq(void)
{
m_pInstances = NULL; // pts to array of NMer instances
m_pSeq = NULL;		// holds all NMer sequences concatenated
Reset();
}

CSfxNMerFreq::CSfxNMerFreq(bool bOverRep,int MinFreqCutoff,int MaxFreqCutoff, int MaxInstances, int NMerLen)
{
m_pInstances = NULL; // pts to array of NMer instances
m_pSeq = NULL;		// holds all NMer sequences concatenated
Init(bOverRep,MinFreqCutoff,MaxFreqCutoff,MaxInstances,NMerLen);
}

CSfxNMerFreq::~CSfxNMerFreq(void)
{
if(m_pInstances != NULL)
	delete m_pInstances;
if(m_pSeq != NULL)
	delete m_pSeq;
}

void 
CSfxNMerFreq::Reset(void)
{
if(m_pInstances != NULL)	// pts to array of NMer instances
	{
	delete m_pInstances;
	m_pInstances = NULL; 
	}

if(m_pSeq != NULL)			// holds all NMer sequences concatenated
	{
	delete m_pSeq;
	m_pSeq = NULL;		
	}
m_NMerLen = 0;			// freqs for NMers of this length
m_MinFreqCutoff = 0;	// minimum freq cutoff
m_MaxFreqCutoff = 0;	// maximum freq cutoff
m_MaxInstances = 0;		// maximum number of NMer instances
m_CurNumInstances= 0;	// current number of NMer instances
m_bOverRep = false;		// true if processing for over represented instances
}

int 
CSfxNMerFreq::Init(bool bOverRep,int MinFreqCutoff,int MaxFreqCutoff, int MaxInstances, int NMerLen)
{
int Idx;
etSeqBase *pSeq;
tsNMerFreq *pNMerFreq;

Reset();
if((m_pInstances = new tsNMerFreq[MaxInstances])==NULL)
	return(eBSFerrMem);
if((m_pSeq = (etSeqBase *) new unsigned char [MaxInstances * NMerLen])==NULL)
	{
	delete m_pInstances;
	return(eBSFerrMem);
	}
pSeq = m_pSeq;
pNMerFreq = m_pInstances;
for(Idx=0; Idx < MaxInstances; Idx++,pNMerFreq++,pSeq+=NMerLen)
	{
	pNMerFreq->Freq = 0;
	pNMerFreq->pSeq = pSeq;
	}
m_NMerLen = NMerLen;				// freqs for NMers of this length
m_MinFreqCutoff = MinFreqCutoff;	// minimum freq cutoff
m_MaxFreqCutoff = MaxFreqCutoff;	// maximum freq cutoff
m_MaxInstances = MaxInstances;		// maximum number of NMer instances
m_CurNumInstances= 0;				// current number of NMer instances
m_bOverRep = bOverRep;				// true if processing for over representated instances
return(eBSFSuccess);
}

int 
CSfxNMerFreq::Add(int Freq,etSeqBase *pSeq)
{
etSeqBase *pFreeSeq;
tsNMerFreq *pEl;
tsNMerFreq Tmp,TmpA;
int Lo,Mid,Hi;	// search limits

// freq must be in range between min and max freq cutoffs
if(Freq > m_MaxFreqCutoff || Freq < m_MinFreqCutoff)
	return(0);

if(!m_CurNumInstances)
	{
	pEl = m_pInstances;
	pEl->Freq = Freq;
	pEl->pSeq = m_pSeq;
	memmove(pEl->pSeq,pSeq,m_NMerLen);
	m_CurNumInstances = 1;
	return(m_CurNumInstances);
	}

if(m_bOverRep) 
	{
	if(m_CurNumInstances == m_MaxInstances && Freq <= m_pInstances[0].Freq)
		return(0);

	if(m_CurNumInstances > 0)
		{
		// can do a binary search because instances are in decreasing Freq order
		Lo = 0; Hi = m_CurNumInstances-1;
		while(Hi >= Lo) 
			{
			Mid = (Hi + Lo)/2;
			pEl = &m_pInstances[Mid];
			if(Freq == pEl->Freq)
				break;
			if(Freq > pEl->Freq)
				{
				Hi = Mid - 1;
				continue;
				}
			Lo = Mid + 1;
			}

		// binary search completed, check if appending, replacing or inserting 
		if(Mid == m_CurNumInstances-1 && // could be appending or replacing
			Freq < pEl->Freq)
			{
			if(m_CurNumInstances != m_MaxInstances) // appending if space, otherwise replacing
				pEl += 1;
			pEl->Freq = Freq;
			memmove(pEl->pSeq,pSeq,m_NMerLen);
			}
		else						// else inserting
			{
			if(m_CurNumInstances == m_MaxInstances)
				pFreeSeq = m_pInstances[m_MaxInstances-1].pSeq;
			else
				pFreeSeq = m_pInstances[m_CurNumInstances].pSeq;
			TmpA.Freq = Freq;
			TmpA.pSeq = pFreeSeq;
			memmove(TmpA.pSeq,pSeq,m_NMerLen);
			for(;Mid < m_CurNumInstances;Mid++)
				{
				Tmp = *pEl;
				*pEl++ = TmpA;
				TmpA = Tmp;
				}
			}
		}

	if(m_CurNumInstances < m_MaxInstances)
		m_CurNumInstances += 1;
	}
else
	{
	if(m_CurNumInstances == m_MaxInstances && Freq >= m_pInstances[m_CurNumInstances-1].Freq)
		return(0);

	if(m_CurNumInstances > 0)
		{
		// can do a binary search because instances are in increasing Freq order
		Lo = 0; Hi = m_CurNumInstances-1;
		while(Hi >= Lo) 
			{
			Mid = (Hi + Lo)/2;
			pEl = &m_pInstances[Mid];
			if(Freq == pEl->Freq)
				break;
			if(Freq < pEl->Freq)
				{
				Hi = Mid - 1;
				continue;
				}
			Lo = Mid + 1;
			}
		// binary search completed, check if appending, replacing or inserting 
		if(Mid == m_CurNumInstances-1 && // could be appending or replacing
			Freq > pEl->Freq)
			{
			if(m_CurNumInstances != m_MaxInstances) // appending if space, otherwise replacing
				pEl += 1;
			pEl->Freq = Freq;
			memmove(pEl->pSeq,pSeq,m_NMerLen);
			}
		else						// else inserting
			{
			if(m_CurNumInstances == m_MaxInstances)
				pFreeSeq = m_pInstances[m_MaxInstances-1].pSeq;
			else
				pFreeSeq = m_pInstances[m_CurNumInstances].pSeq;
			TmpA.Freq = Freq;
			TmpA.pSeq = pFreeSeq;
			memmove(TmpA.pSeq,pSeq,m_NMerLen);
			for(;Mid < m_CurNumInstances;Mid++)
				{
				Tmp = *pEl;
				*pEl++ = TmpA;
				TmpA = Tmp;
				}
			}
		}

	if(m_CurNumInstances < m_MaxInstances)
		m_CurNumInstances += 1;
	}
return(m_CurNumInstances);
}

int
CSfxNMerFreq::GetNumInstances(void)
{
return(m_CurNumInstances);
}

int
CSfxNMerFreq::GetNMerLen(void)
{
return(m_NMerLen);
}

// constructor
CSfxArray::CSfxArray(void)
{
m_pSfxArray = NULL;
m_pSfxSeq = NULL;
m_pPutCores = NULL;
Reset(false);
}

//destructor
CSfxArray::~CSfxArray(void)
{
if(m_pSfxArray != NULL) 							// currently loaded suffix array
	delete m_pSfxArray;
if(m_pSfxSeq != NULL)
	delete m_pSfxSeq;
if(m_pPutCores != NULL)
	delete m_pPutCores;
}

int
CSfxArray::Reset(bool bFlush)			// reset state back to that immediately following instantiation
{
int Rslt;
if((Rslt = Close(bFlush)) != eBSFSuccess)
	return(Rslt);

memset(&m_CurSfxHeader,0,sizeof(m_CurSfxHeader));	// currently loaded suffix header
if(m_pSfxArray != NULL) 							// currently loaded suffix array
	{
	delete m_pSfxArray;
	m_pSfxArray = NULL;
	}

if(m_pSfxSeq != NULL) 								// currently loaded suffix sequence
	{
	delete m_pSfxSeq;
	m_pSfxSeq = NULL;
	}

if(m_pPutCores != NULL) 							// pts to memory allocated to hold putative matching cores
	{
	delete m_pPutCores;
	m_pPutCores = NULL;
	}

m_NumPutCores = 0;
m_MaxPutCores = 0;
m_CurMaxSfxSize = 0;								
m_CurMaxSeqSize = 0;
m_CurSfxSeqChromID = 0;
m_CurSfxArrayChromID = 0;
m_ProbeRepeatPC = 100;
m_TargRepeatPC = 100;
return(eBSFSuccess);
}

// Open()
// Opens specified suffix file pszFile
// Option to create or truncate pszFile
int
CSfxArray::Open(char *pszFile,
					  bool bCreate)
{
Reset(false);
return(CBioSeqFile::Open(pszFile,cBSFTypeSfx,bCreate));
}



int 
CSfxArray::Close(bool bFlush)
{
return(CBioSeqFile::Flush2Disk());
}


int											
CSfxArray::AddEntry(char *pszSource,		// from where was this sequence originated
				char *pszDescription,		// descriptive text about this sequence
				etSeqBase   *pSeq,			// sequence to generate suffix indexes over (MUST be terminated by eBaseEOS)
				unsigned int SeqLen)		// sequence length excluding eBaseEOS
{
unsigned int BlockID;
unsigned int SeqOfs;
unsigned int *pSfxArray;

int Rslt;

CSAIS sais;

pSfxArray = new unsigned int [SeqLen];	// array to hold suffix ptrs
if(pSfxArray == NULL)
	return(eBSFerrMem);

BlockID = 1;

m_CurSfxHeader.ChromID  = CBioSeqFile::CreateEntry(pszSource,pszDescription,eBinDataType);
m_CurSfxHeader.SeqLen   = SeqLen;
m_CurSfxHeader.ElSize = sizeof(UINT32);

if(m_bIsBigEndian)
	{
	m_CurSfxHeader.ChromID=SwapUI32Endians(m_CurSfxHeader.ChromID);		// identifies chromosome (1..n) within this suffix file
	m_CurSfxHeader.ElSize=SwapUI32Endians(m_CurSfxHeader.ElSize);		// size of each suffix array element (in bytes)
	m_CurSfxHeader.SeqLen=SwapUI32Endians(m_CurSfxHeader.SeqLen);		// sequence length (also suffix array size in # of elements)
	}
CBioSeqFile::AddData(sizeof(tsSfxHeader),(unsigned char *)&m_CurSfxHeader);
if(m_bIsBigEndian)
	{
	m_CurSfxHeader.ChromID=SwapUI32Endians(m_CurSfxHeader.ChromID);		// identifies chromosome (1..n) within this suffix file
	m_CurSfxHeader.ElSize=SwapUI32Endians(m_CurSfxHeader.ElSize);		// size of each suffix array element (in bytes)
	m_CurSfxHeader.SeqLen=SwapUI32Endians(m_CurSfxHeader.SeqLen);		// sequence length (also suffix array size in # of elements)
	}

SeqOfs = 0;

sais.sais(pSeq,(int *)pSfxArray,		// sort into this prealloc'd suffix array
				SeqLen);	

if(m_bIsBigEndian)
	{
	unsigned int Idx;
	UINT32 *p32 = pSfxArray;
	for(Idx = 0; Idx < SeqLen; Idx++, p32++)
		*p32 = SwapUI32Endians(*p32);
	}
if((Rslt = CBioSeqFile::AddData(SeqLen * sizeof(unsigned int),(unsigned char *)pSfxArray)) >= eBSFSuccess)
	{
	if((Rslt = CBioSeqFile::AddData(SeqLen,pSeq)) >= eBSFSuccess)
		Rslt = CBioSeqFile::SealEntry();
	}

if(pSfxArray)
	delete pSfxArray;
	
return(Rslt);
}


tsSfxHeader *
CSfxArray::LoadSfxHeader(unsigned int EntryID)
{
unsigned int DataReq;
unsigned int DataLen;
if(m_CurSfxHeader.ChromID == EntryID)
	return(&m_CurSfxHeader);
m_CurSfxSeqChromID = 0;
DataReq = sizeof(tsSfxHeader);
if((DataLen = GetDataLen(EntryID)) < DataReq)
	return(NULL);
if(DataReq != GetData(EntryID,eBinDataType,0,(unsigned char *)&m_CurSfxHeader,DataReq))
	return(NULL);

if(m_bIsBigEndian)
	{
	m_CurSfxHeader.ChromID=SwapUI32Endians(m_CurSfxHeader.ChromID);		// identifies chromosome (1..n) within this suffix file
	m_CurSfxHeader.ElSize=SwapUI32Endians(m_CurSfxHeader.ElSize);		// size of each suffix array element (in bytes)
	m_CurSfxHeader.SeqLen=SwapUI32Endians(m_CurSfxHeader.SeqLen);		// sequence length (also suffix array size in # of elements)
	}

return(&m_CurSfxHeader);
}


unsigned int *
CSfxArray::LoadSfxArray(unsigned int EntryID)
{
unsigned int SfxPsn;
unsigned int DataReq;

if(LoadSfxHeader(EntryID)==NULL)
	return(NULL);

if(m_CurSfxArrayChromID == EntryID)
	return(m_pSfxArray);

m_CurSfxArrayChromID = 0;
SfxPsn = sizeof(tsSfxHeader);
DataReq = m_CurSfxHeader.SeqLen * m_CurSfxHeader.ElSize;
if(m_pSfxArray == NULL || DataReq > m_CurMaxSfxSize)
	{
	if(m_pSfxArray != NULL)
		delete m_pSfxArray;
	m_pSfxArray = new unsigned int [m_CurSfxHeader.SeqLen];
	if(m_pSfxArray == NULL)
		return(NULL);
	m_CurMaxSfxSize = DataReq;
	}

if(DataReq != GetData(EntryID,eBinDataType,SfxPsn,(unsigned char *)m_pSfxArray,DataReq))
	return(NULL);

if(m_bIsBigEndian)
	{
	unsigned int Idx;
	UINT32 *p32 = m_pSfxArray;
	for(Idx = 0; Idx < m_CurSfxHeader.SeqLen; Idx++, p32++)
		*p32 = SwapUI32Endians(*p32);
	}
m_CurSfxArrayChromID = EntryID;
return(m_pSfxArray);
}

int
CSfxArray::GetSfxSeqLen(unsigned int EntryID)
{
if(LoadSfxHeader(EntryID) == NULL)
	return(eBSFerrEntry);
return(m_CurSfxHeader.SeqLen);
}

int
CSfxArray::GetSfxSeq(unsigned int EntryID,UINT32 Offset,etSeqBase *pRetSeq,UINT32 Len,bool bClearSeqMsk)
{
etSeqBase *pSeq;
if((pSeq = LoadSfxSeq(EntryID,bClearSeqMsk)) == NULL)
	return(eBSFerrInternal);
if(Offset >= m_CurSfxHeader.SeqLen)
	return(eBSFerrInternal);
Len = min(Len,m_CurSfxHeader.SeqLen - Offset);
memmove(pRetSeq,&pSeq[Offset],Len);
return(Len);
}

etSeqBase *
CSfxArray::LoadSfxSeq(unsigned int EntryID,bool bClearSeqMsk)
{
unsigned int SeqPsn;
unsigned int DataReq;
etSeqBase *pSeq;

if(LoadSfxHeader(EntryID) == NULL)
	return(NULL);

if(m_CurSfxSeqChromID == EntryID)
	return(m_pSfxSeq);
m_CurSfxSeqChromID = 0;

SeqPsn = sizeof(tsSfxHeader);
SeqPsn += m_CurSfxHeader.SeqLen * m_CurSfxHeader.ElSize; 
DataReq = m_CurSfxHeader.SeqLen;
if(m_pSfxSeq == NULL || DataReq >= m_CurMaxSeqSize)
	{
	if(m_pSfxSeq != NULL)
		delete m_pSfxSeq;
	m_pSfxSeq = new etSeqBase[DataReq + 1];
	if(m_pSfxSeq == NULL)
		return(NULL);
	m_CurMaxSeqSize = DataReq + 1;
	}

if(DataReq != GetData(EntryID,eBinDataType,SeqPsn,(unsigned char *)m_pSfxSeq,DataReq))
	return(NULL);

m_pSfxSeq[DataReq] = eBaseEOS;
m_CurSfxSeqChromID = EntryID;

if(bClearSeqMsk)	// if need to clear any sequence repeat flags
	{
	pSeq = m_pSfxSeq;
	for(SeqPsn = 0; SeqPsn < DataReq; SeqPsn++,pSeq++)
		*pSeq &= ~cRptMskFlg;
	}
return(m_pSfxSeq);
}


bool
CSfxArray::AcceptMasked(etSeqBase *pTarg,etSeqBase *pProbe,int MatchLen)
{
etSeqBase CurBase;
int NumMasked;
int MaxAllowed;
int Cnt;
if(m_ProbeRepeatPC == 100 && m_TargRepeatPC == 100)
	return(true);
if(m_ProbeRepeatPC < 100)
	{
	if(m_ProbeRepeatPC > 0)
		MaxAllowed = 1 + ((m_ProbeRepeatPC * MatchLen)/100);
	else
		MaxAllowed = 0;
	NumMasked = 0;
	for(Cnt = 0; Cnt < MatchLen; Cnt++,pProbe++)
		{
		CurBase = *pProbe;
		if((!m_bProbeRepeatSense && CurBase & cRptMskFlg) ||
			(m_bProbeRepeatSense && !(CurBase & cRptMskFlg)))
			NumMasked++;
		if(NumMasked > MaxAllowed)
			return(false);
		}
	}

if(m_TargRepeatPC < 100)
	{
	if(m_TargRepeatPC > 0)
		MaxAllowed = 1 + ((m_TargRepeatPC * MatchLen)/100);
	else
		MaxAllowed = 0;
	NumMasked = 0;
	for(Cnt = 0; Cnt < MatchLen; Cnt++,pTarg++)
		{
		CurBase = *pTarg;
		if((!m_bTargRepeatSense && CurBase & cRptMskFlg) ||
			(m_bTargRepeatSense && !(CurBase & cRptMskFlg)))
			NumMasked++;
		if(NumMasked > MaxAllowed)
			return(false);
		}
	}
return(true);
}


// ProcessMatch
// Process a putative match
// Extends match right and left as appropriate
// Then adds match to results file (m_pRsltsFile->AddMatch or m_pRsltsFile->NearlyExtend)
// Note: If putative match does not meet repeat mask (m_ProbeRepeatPC and m_TargRepeatPC) then match is sloughed and false returned
//
bool													// true if match added to results file
CSfxArray::ProcessMatch(etSeqBase *pProbe,				// probe sequence
						unsigned int ProbePsn,			// offset in probe sequence at which seed match starts
						etSeqBase *pTarg,				// target sequence
						unsigned int TargPsn,			// offset in target at which seed match starts
						unsigned int MatchLen,			// match length
						unsigned int MaximalMatchLen)	// extend matches left/right to this maximum
{
unsigned int LeftExtension;
unsigned int RightExtension;
unsigned int IdentCnt;
etSeqBase *pP;
etSeqBase *pT;

// check for self hit (same chromosome against itself)
// assumption is that if match is at the limit of what is feasable and if at same
// position on chromosomes with same identifier then it must be a self match
if(m_ProbeEntryID == m_TargEntryID && 
   TargPsn == ProbePsn && 
   MatchLen == m_CurMaxMatchLen)
	{
	m_ChkSelfMatches = true;	// assume self against self for remainder of this chromosome
	return(false);
	}

// right extend the match only if match allows for mutations
if(m_bDontExtend || !(m_MaxMutsCnt | m_MaxNsCnt))
	RightExtension = 0;
else
	RightExtension = ExtendRight(&pProbe[ProbePsn+MatchLen], &pTarg[TargPsn+MatchLen], MaximalMatchLen - MatchLen);

// left extend the match
if(m_bDontExtend || MatchLen >= MaximalMatchLen || ProbePsn == 0 || TargPsn == 0)
	LeftExtension = 0;
else
	LeftExtension = ExtendLeft(pProbe, ProbePsn, pTarg, TargPsn, MaximalMatchLen - MatchLen) - 1;

// if match length less than minimum required then time to packup and go home..
MatchLen += (RightExtension + LeftExtension);
if(MatchLen < m_CurMinMatchLen)
	return(false);

TargPsn -= LeftExtension;
ProbePsn -= LeftExtension;

		// adjust probe psns according to the mode
		// in antisense and reverse binding the probe sequence was reversed
unsigned int AdjProbePsn =	m_CurMatchMode >= 2 ? m_ProbeLen - ProbePsn - MatchLen : ProbePsn;

// slough if contained within any previous match
if(!m_bAllowContained && (*m_pRsltsFileHandler)(eEMRIsAlreadyMatched,AdjProbePsn,TargPsn,MatchLen,0,m_CurMatchMode))
//if(!m_bAllowContained && m_pRsltsFile->IsAlreadyMatched(AdjProbePsn,TargPsn,MatchLen,m_CurMatchMode))
	return(false);

// ensure meets repeat masked constraints
if(!AcceptMasked(&pTarg[TargPsn],&pProbe[ProbePsn],MatchLen))
	return(false);

// calculate the identity (number of bases exactly matching between probe and target, will be <= MatchLen)
if(!(m_MaxMutsCnt | m_MaxNsCnt))
	IdentCnt = MatchLen;
else
	{
	pT = &pTarg[TargPsn];
	pP = &pProbe[ProbePsn];
	IdentCnt = 0;
	for(unsigned int Cnt = 0; Cnt < MatchLen; Cnt++)

		if((*pP++ & ~cRptMskFlg) == (*pT++ & ~cRptMskFlg))
			IdentCnt++;

	// ensure exact identies >= min
	if(IdentCnt < m_CurMinMatchLen)
		return(false);

if((*m_pRsltsFileHandler)(eEMRNearlyExtend,AdjProbePsn,TargPsn,MatchLen,IdentCnt,m_CurMatchMode))
//	if(m_pRsltsFile->NearlyExtend(AdjProbePsn,TargPsn,MatchLen,IdentCnt,m_CurMatchMode))
		return(true);
	}
return((*m_pRsltsFileHandler)(eEMRAddMatch,AdjProbePsn,TargPsn,MatchLen,IdentCnt,m_CurMatchMode)==eBSFSuccess?true:false);

//return(m_pRsltsFile->AddMatch(AdjProbePsn,TargPsn,MatchLen,IdentCnt,m_CurMatchMode)==eBSFSuccess?true:false);
}

unsigned int								// #bases extended 
CSfxArray::ExtendLeft(etSeqBase *pProbe,	// sequence containing probe
		  unsigned int ProbePsn,			// start position in probe sequence at which to extend left 
			  etSeqBase *pTarg,				// target sequence
			  unsigned int TargPsn,			// start position in target sequence at which to extend left
			  unsigned int Max2Extend)		// extend left for at most this many bases
{
unsigned int CurMatchLen;
unsigned int Max2Compare;
etSeqBase b1;
etSeqBase b2;

if(!Max2Extend)
	return(0);

Max2Compare = 1 + (ProbePsn <= TargPsn ? ProbePsn : TargPsn);
if(Max2Compare > Max2Extend)
	Max2Compare = Max2Extend;

CurMatchLen = 0;
pProbe += ProbePsn;
pTarg += TargPsn;

// a) extension is on exact matches - no mutations allowed
if((m_MaxMutsCnt | m_MaxNsCnt) == 0)
	{
#ifdef thiscode
	while(Max2Compare-- && 
		((b1 = (*pProbe-- & ~cRptMskFlg)) == (b2 = (*pTarg-- & ~cRptMskFlg)) || 
		  bRNArules && ((b1 == eBaseG && b2 == eBaseT) || (b1 == eBaseT && b2 == eBaseG))) && 
		b1 < eBaseN && 
		b2 < eBaseN)
		CurMatchLen++;
#else
	while(Max2Compare-- && 
		(b1 = (*pProbe-- & ~cRptMskFlg)) == (b2 = (*pTarg-- & ~cRptMskFlg)) && 
		b1 < eBaseN && b2 < eBaseN)
		CurMatchLen++;
#endif
	}
else
	{
	// b) extension allows for limited number of mutations m_MaxMutsLen and indeterminant bases 'N'
	// If a run of consequative mutations is longer then m_MaxMutBlockLen then match is terminated
	unsigned int TotMuts = 0;
	unsigned int TotNs = 0;
	unsigned int MutBlockLen = 0;
	unsigned int ProvMatchLen = 0;


	while(Max2Compare-- && 
		((b1 = (*pProbe-- & ~cRptMskFlg)) <= eBaseN) && 
		((b2 = (*pTarg-- & ~cRptMskFlg)) <= eBaseN))
		{
		if(b1 == eBaseN || b2 == eBaseN)
			{
			if(++TotNs > m_MaxNsCnt)
				break;
			}
		else
			{

			if(b1 != b2)
				{
				if(++TotMuts > m_MaxMutsCnt ||
					++MutBlockLen > m_MaxMutBlockLen)
					break;
				}
			else			// else is an exact match so record matchlen as being provisional match length
				{
				MutBlockLen = 0;
				CurMatchLen = ProvMatchLen + 1;
				}
			}
		ProvMatchLen++;
		}
	}

return(CurMatchLen);
}


unsigned int								// #bases extended 
CSfxArray::ExtendRight(etSeqBase *pProbe,	// sequence containing probe
			  etSeqBase *pTarg,				// target sequence
			  unsigned int Max2Extend)		// extend right for at most this many bases
{
unsigned int CurMatchLen;
unsigned int Max2Compare;
etSeqBase b1;
etSeqBase b2;

if(!Max2Extend)
	return(0);

Max2Compare = Max2Extend;

CurMatchLen = 0;

// a) extension is on exact matches - no mutations allowed
if((m_MaxMutsCnt | m_MaxNsCnt) == 0)
	{

	while(Max2Compare-- && 
		(b1 = (*pProbe++ & ~cRptMskFlg)) == (b2 = (*pTarg++ & ~cRptMskFlg)) && 
		b1 < eBaseN && 
		b2 < eBaseN)
		CurMatchLen++;
	}
else
	{
	// b) extension allows for limited number of mutations m_MaxMutsLen and indeterminant bases 'N'
	// If a run of consequative mutations is longer then m_MaxMutBlockLen then match is terminated
	unsigned int TotMuts = 0;
	unsigned int TotNs = 0;
	unsigned int MutBlockLen = 0;
	unsigned int ProvMatchLen = 0;
	while(Max2Compare-- && ((b1 = (*pProbe++ & ~cRptMskFlg)) <= eBaseN) && ((b2 = (*pTarg++ & ~cRptMskFlg)) <= eBaseN))
	{
		if(b1 == eBaseN || b2 == eBaseN)
			{
			if(++TotNs > m_MaxNsCnt)
				break;
			}
		else
			{
			if(b1 != b2)
				{
				if(++TotMuts > m_MaxMutsCnt ||
					++MutBlockLen > m_MaxMutBlockLen)
					break;
				}
			else			// else is an exact match so record matchlen as being provisional match length
				{
				MutBlockLen = 0;
				CurMatchLen = ProvMatchLen + 1;
				}
			}
		ProvMatchLen++;
		}
	}

return(CurMatchLen);
}


int											// returns number of matches					
CSfxArray::Locate( etSeqBase *pProbe,		// sequence containing probe
				  unsigned int ProbePsn,	// start position in probe sequence at which to search 
				  etSeqBase *pTarg,			// target sequence
				  unsigned int *pSfxArray,	// target sequence suffix array
				  unsigned int TargStart,   // position in pTarg (0..n) corresponding to start of suffix array
				  unsigned int SfxLen)		// number of suffixs in pSfxArray
{
etSeqBase *pEl1;
etSeqBase *pEl2;
etSeqBase *pSeq2Loc;

etSeqBase b1;
etSeqBase b2;
int CurMatchLen;
int Max2Compare;
int NumMatches;
int SfxLo;
int SfxHi;
int TargPsn;
int MatchPsn;

pSeq2Loc = &pProbe[ProbePsn];

if((*pSeq2Loc & ~cRptMskFlg) >= eBaseN)
	return(0);

NumMatches = 0;
SfxLo = 0;
SfxHi = SfxLen - 1;
do {
	TargPsn = (SfxLo + SfxHi) / 2;
	pEl2 = &pTarg[pSfxArray[TargPsn + TargStart]];
	if(m_bDontExtend || (m_ChkSelfMatches && 
			ProbePsn == pSfxArray[TargPsn + TargStart] &&
			m_ProbeEntryID == m_TargEntryID))
		Max2Compare = m_CurSampleSize;
	else
		Max2Compare = m_CurMaxMatchLen;
	CurMatchLen = 0;
	pEl1 = pSeq2Loc;

	while(Max2Compare-- && 
			(b1 = (*pEl1++ & ~cRptMskFlg)) == (b2 = (*pEl2++ & ~cRptMskFlg)) && 
			b1 < eBaseN && 
			b2 < eBaseN)
		CurMatchLen++;

	if(CurMatchLen >= (int)m_CurSampleSize)
		{
		if(m_ChkSelfMatches && 
			ProbePsn == pSfxArray[TargPsn + TargStart] &&
			m_ProbeEntryID == m_TargEntryID)
			{
			NumMatches = 0;
			}
		else
			{
			if(ProcessMatch(pProbe,
							ProbePsn,
							pTarg,
							pSfxArray[TargPsn + TargStart],
							CurMatchLen)) // process the match
				NumMatches = 1;
			}

		MatchPsn = TargPsn;
		while(TargPsn > 0 && NumMatches < (int)m_CurMaxNumMatches && --TargPsn >= SfxLo)
			{
			if(m_ChkSelfMatches && 
				ProbePsn == pSfxArray[TargPsn + TargStart] &&
				m_ProbeEntryID == m_TargEntryID)
				continue;

			Max2Compare = m_bDontExtend ? m_CurSampleSize : m_CurMaxMatchLen;
			CurMatchLen = 0;
			pEl1 = pSeq2Loc;
			pEl2 = &pTarg[pSfxArray[TargPsn + TargStart]];
			while(Max2Compare-- && 
				(b1 = (*pEl1++ & ~cRptMskFlg)) == (b2 = (*pEl2++ & ~cRptMskFlg)) && 
				b1 < eBaseN && 
				b2 < eBaseN)
				CurMatchLen++;
			if(CurMatchLen >= (int)m_CurSampleSize)
				{
				if(ProcessMatch(pProbe,ProbePsn,pTarg,pSfxArray[TargPsn + TargStart],CurMatchLen)) // process the match
					NumMatches++;
				}
			else
				break;
			}

		TargPsn = MatchPsn;	
		while(NumMatches < (int)m_CurMaxNumMatches && ++TargPsn <= SfxHi)
			{
			if(m_ChkSelfMatches && 
				ProbePsn == pSfxArray[TargPsn + TargStart] &&
				m_ProbeEntryID == m_TargEntryID)
				continue;

			Max2Compare = m_bDontExtend ? m_CurSampleSize : m_CurMaxMatchLen;
			CurMatchLen = 0;
			pEl1 = pSeq2Loc;
			pEl2 = &pTarg[pSfxArray[TargPsn + TargStart]];
			while(Max2Compare-- && 
				(b1 = (*pEl1++ & ~cRptMskFlg)) == (b2 = (*pEl2++ & ~cRptMskFlg)) && 
				b1 < eBaseN && 
				b2 < eBaseN)
				CurMatchLen++;
			if(CurMatchLen >= (int)m_CurSampleSize)
				{
				if(ProcessMatch(pProbe,ProbePsn,pTarg,pSfxArray[TargPsn + TargStart],CurMatchLen)) // process the match
					NumMatches++;
				}
			else
				break;
			}

		return(NumMatches);				// return number of matches..
		}

	switch(b1) {
		case eBaseN:
			if(b2 < eBaseN)
				{
				SfxHi = TargPsn - 1;
				continue;
				}
			SfxLo = TargPsn+1;
			continue;

		case eBaseEOS:
			SfxHi = TargPsn - 1;
			continue;

		default:
			if(b2 >= eBaseN)
				SfxLo = TargPsn+1;
			else
				{
				if(b1 < b2)
					SfxHi = TargPsn - 1;
				else
					SfxLo = TargPsn+1;
				}
			continue;
		}
	}
while(SfxHi >= SfxLo);
return(0);
}


int														// returns number of matches
CSfxArray::LocateMatches(unsigned int EntryID,			// target (suffix) entry identifier
					     unsigned int ProbeID,			// probe entry identifier, uniquely identify this probe when ProcessThisMatch is called
					     unsigned int MinMatchLen,		// matches located must be at least this length
 						 unsigned int MinExactLen,		// and contain a subsequence of at least this 100% identity length
						 unsigned int MaxNs,			// maximum number of indeterminate bases allowed
						 unsigned int MaxMutations,		// maximum number of mutations allowed
						 unsigned int MaxMutBlockLen,	// max length of any mutated sequence block
   						 unsigned int MaxNumMatchesPPS,	// max number of matches per MinMatchLen probe subsequence
						 etSeqBase *pProbe,				// probe sequence (ProbeLen bases)
						unsigned int ProbeLen,			// number of bases ptd to by pProbe
						unsigned int Modes,				// bit 0 fwd (sense), bit 1 fwd rev, bit 2 cpl, bit 3 cpl rev (antisense)
						int RsltsFileHandler(int Mode,unsigned int ProbePsn,unsigned int TargPsn,unsigned int MatchLen,unsigned int IdentCnt,unsigned int FHMode),
						bool bChkSelfMatches,			// true if self matches to be filtered out
						bool bFiltComplexity,			// true if low complexity sequences to be filtered
					    bool bNotSampled,				// true if sampling not to be used (use sliding window)
						bool bAllowContained,			// true if contained sequences are to be reported
						bool bDontExtend,				// true if matches are not to be maximally extended
					    bool bSeqOut,					// true if matching sequences to be output in results
						bool bProbeRepeatSense,			// false if bases marked as masked to be treated as repeats, otherwise unmarked bases treated as repeats
						bool bTargRepeatSense,			// false if bases marked as masked to be treated as repeats, otherwise unmarked bases treated as repeats
						unsigned int ProbeRepeatPC, 	// percentage of bases in any probe matches allowed to have been repeat masked (0..100)
    					unsigned int TargRepeatPC)		// percentage of bases in any (Suffix) target matches allowed to have been repeat masked (0..100)
{
unsigned int ProbeOfs;
unsigned int SfxBlockID;
unsigned int TargStart;
unsigned int ProbeLenRemaining;
unsigned int Idx;

unsigned int OverRepresented;
unsigned int UnderRepresented;
int NumMasked;											// current number of bases marked as repeats in window
int MaxMasked;											// maximum allowed repeat bases in current window

etSeqBase CurBase;
unsigned int NumFiltered;
int NumMatches;
unsigned int BaseCnt[6];
int CurMode;

bool bComplemented;
bool bReversed;
bool bWindowSampled;

if(MinMatchLen < 5)
	MinMatchLen = 5;
 if(MinExactLen < 5)
	 MinExactLen = 5;
 else
	 if(MinExactLen > MinMatchLen)
		MinExactLen = MinMatchLen;

if(MaxMutations == 0 && MaxNs == 0)		// if no mismatches allowed then the minimum exact
	MinExactLen = MinMatchLen;			// length is forced to be same as MinMatchLen

if(ProbeLen < MinMatchLen)
	return(0);

if(m_MaxMutsCnt && !MaxMutBlockLen)
	MaxMutBlockLen = m_MaxMutsCnt;

if(MaxMutBlockLen > m_MaxMutsCnt)
	MaxMutBlockLen = m_MaxMutsCnt;

if(ProbeRepeatPC > 100)
	ProbeRepeatPC = 100;
if(TargRepeatPC > 100)
	TargRepeatPC = 100;
m_ProbeRepeatPC = ProbeRepeatPC;
m_TargRepeatPC = TargRepeatPC;
m_bProbeRepeatSense = bProbeRepeatSense;
m_bTargRepeatSense = bTargRepeatSense;

m_ProbeLen = ProbeLen;

m_ProbeEntryID = ProbeID;
m_TargEntryID = EntryID;

m_ChkSelfMatches = bChkSelfMatches;

m_Modes  = Modes;
m_CurMinMatchLen = MinMatchLen;			// minimal required match length
m_CurMaxNumMatches = MaxNumMatchesPPS;	// only look for this many matches

m_MaxNsCnt = MaxNs;						// maximum number of indeterminate nts allowed
m_MaxMutsCnt = MaxMutations;			// maximum number of mutations allowed when extending
m_MaxMutBlockLen = MaxMutBlockLen;		// maximum number of mutations alllowed in any sequence block
m_MinExactSeedLen = MinExactLen;		// minimum seed length for identity extension 
m_bAllowContained = bAllowContained;	// true if contained sequences are to be reported
m_bDontExtend = bDontExtend;			// true if matches are not to be maximally extended
m_bSeqOut = bSeqOut;					// true if matching sequences to be output in results

m_pRsltsFileHandler = RsltsFileHandler;

if(LoadSfxHeader(EntryID) == NULL)
	return(-1);

SfxBlockID = 1;

if(LoadSfxSeq(EntryID) == NULL)
	return(-1);

if(bDontExtend)							// if not allowing extension of matches then force sliding fixed window 
	bNotSampled = true;

if(bNotSampled || MaxNs || MaxMutations || 
   m_MinExactSeedLen < 12)				// if min match length is too small then it is
	{									// more efficent to use a sliding window rather than a sampling window
	m_CurSampleSize = m_MinExactSeedLen;// as a sampling window looks for seed matches of 0.5 the window size
	bWindowSampled = false;				// then with small windows there are too many false positives on the seed matches
	}									// when extended out to MinMatchLen
else
	{
	bWindowSampled = true;
	m_CurSampleSize = m_MinExactSeedLen / 2;
	}

if(bFiltComplexity)
	{
	OverRepresented = 1 + ((m_CurSampleSize * 7) / 10);   // approx 70%
	UnderRepresented = 1 + ((m_CurSampleSize * 4) / 100); // approx 4%
	}

if(ProbeRepeatPC < 100)
	{
	if(ProbeRepeatPC > 0)
		MaxMasked = 1+((ProbeRepeatPC * m_MinExactSeedLen)/100);
	else
		MaxMasked = 0;
	}
else
	MaxMasked = -1;

NumMatches = 0;
bComplemented = false;
bReversed = false;

if(LoadSfxArray(EntryID) == NULL)
	return(-1);

for(CurMode = 0; CurMode < 4; CurMode++)
	{
	switch(CurMode) {
		case 0:				// (sense) forward p:agctg against t:agctg
			if(!(Modes & 0x01))
				continue;

			if(bComplemented)
				{
				CSeqTrans::ComplementStrand(ProbeLen,pProbe);
				bComplemented = false;
				}
			if(bReversed)
				{
				CSeqTrans::ReverseSeq(ProbeLen,pProbe);
				bReversed = false;
				}
			break;

		case 1:				// (ps 5'->3',5'->3') forward complement p:tcgac against t:agctg 
			if(!(Modes & 0x02))
				continue;
			if(!bComplemented)
				{
				CSeqTrans::ComplementStrand(ProbeLen,pProbe);
				bComplemented = true;
				}
			if(bReversed)
				{
				CSeqTrans::ReverseSeq(ProbeLen,pProbe);
				bReversed = false;
				}
			break;

		case 2:				// (antisense) reverse complement p:cagct against t:agctg
			if(!(Modes & 0x04))
				continue;

			if(!bComplemented)
				{
				CSeqTrans::ComplementStrand(ProbeLen,pProbe);
				bComplemented = true;
				}
			if(!bReversed)
				{
				CSeqTrans::ReverseSeq(ProbeLen,pProbe);
				bReversed = true;
				}
			break;
				
		case 3:				// (inversion) reverse p:gtcga against t:agctg
			if(!(Modes & 0x08))
				continue;
			if(!bReversed)
				{
				CSeqTrans::ReverseSeq(ProbeLen,pProbe);
				bReversed = true;
				}
			if(bComplemented)
				{
				CSeqTrans::ComplementStrand(ProbeLen,pProbe);
				bComplemented = false;
				}
			break;
		}

	m_CurMatchMode = CurMode;
	TargStart = 0;
	ProbeOfs = 0;		
	ProbeLenRemaining = ProbeLen;
	NumMasked = 0;
	memset(BaseCnt,0,sizeof(BaseCnt));
	NumFiltered = 0;
	for(Idx = 0; Idx < (m_CurSampleSize-1); Idx++)
		{
		CurBase = pProbe[Idx];

		if((!bProbeRepeatSense && CurBase & cRptMskFlg) ||
			(bProbeRepeatSense && !(CurBase & cRptMskFlg)))
				NumMasked++;
		CurBase &= ~cRptMskFlg;
		if(CurBase > eBaseN)
			CurBase = eBaseN;
		BaseCnt[CurBase]++;
		}
				
	while(ProbeLenRemaining >= m_CurSampleSize) 
		{
		m_CurMaxMatchLen = ProbeLenRemaining;
		if(ProbeLenRemaining > m_CurSampleSize)
			{
			CurBase = pProbe[ProbeOfs+m_CurSampleSize-1];
			if((!bProbeRepeatSense && CurBase & cRptMskFlg) ||
				(bProbeRepeatSense && !(CurBase & cRptMskFlg)))
				NumMasked++;
			CurBase &= ~cRptMskFlg;
			if(CurBase > eBaseN)
				CurBase = eBaseN;
			BaseCnt[CurBase]++;
			}

		if(!bWindowSampled || !(ProbeOfs % m_CurSampleSize))
			{
			bool DoIt = true;					// assume will call Locate unless some condition fails
				// if left edge of probe base is indeterminate, or  sampling and any base in window is indeterminate
				// then we already know we won't get an exact match so don't wast time on the call to Locate
			if(BaseCnt[eBaseN])	
				DoIt = false;

			if(DoIt && bFiltComplexity)
				{
				for(Idx = 0; DoIt && Idx < eBaseN; Idx++)
					{
					if(BaseCnt[Idx] >= OverRepresented ||
						BaseCnt[Idx] <= UnderRepresented)
						DoIt = false;
					}
				}
			
				// if window has more repeats than allowed then no point in trying to locate any matches
			if(DoIt && ProbeRepeatPC < 100)
				{
				if(NumMasked > MaxMasked)
					DoIt = false;
				}

			if(DoIt)
				{
				NumMatches += Locate(pProbe,					// sequence to locate
						ProbeOfs,				// offset in pProbe at which to match from
						m_pSfxSeq,				// target sequence
						m_pSfxArray,			// target sequence suffix array
						TargStart,				// position in pTarg (0..n) corresponding to start of suffix array
						m_CurSfxHeader.SeqLen);// number of suffixs in pSfxArray
				}
			else
				NumFiltered++;
			}

		ProbeOfs+=1;
		ProbeLenRemaining -= 1;
		if(ProbeLenRemaining >= m_CurSampleSize)
			{
			CurBase = pProbe[ProbeOfs-1];
			if((!bProbeRepeatSense && CurBase & cRptMskFlg) ||
				(bProbeRepeatSense && !(CurBase & cRptMskFlg)))
				NumMasked--;
			
			CurBase &= ~cRptMskFlg;
			if(CurBase > eBaseN)
				CurBase = eBaseN;
			if(BaseCnt[CurBase])
				BaseCnt[CurBase]--;
			}
		}
	}

if(bReversed)
	CSeqTrans::ReverseSeq(ProbeLen,pProbe);
if(bComplemented)
	CSeqTrans::ComplementStrand(ProbeLen,pProbe);

return(NumMatches);
}

//LocateNoMatches
//Locate longest non-overlapping sequences in probe which are not in target
bool 
CSfxArray::LocateNoMatches(unsigned int EntryID,		// target entry identifier
					     unsigned int ProbeID,			// probe entry identifier, uniquely identify this probe when ProcessThisMatch is called
					     unsigned int MinMatchLen,		// non-matches located must be at least this length
   						 unsigned int MaxNumMatchesPPS,	// max number of matches per MinMatchLen probe subsequence
						 etSeqBase *pProbe,				// probe sequence (ProbeLen bases)
						unsigned int ProbeLen,			// number of bases ptd to by pProbe
						unsigned int Modes,				// bit 0 fwd (sense), bit 1 fwd rev, bit 2 cpl, bit 3 cpl rev (antisense)
						int RsltsFileHandler(int Mode,unsigned int ProbePsn,unsigned int TargPsn,unsigned int MatchLen,unsigned int IdentCnt,unsigned int FHMode),
						bool bChkSelfMatches)			// true if self matches to be filtered out
{
unsigned int ProbeOfs;
unsigned int SfxBlockID;
unsigned int TargStart;
unsigned int ProbeLenRemaining;
bool bComplemented;
bool bReversed;
int CurMode;

if(MinMatchLen < 8)
	MinMatchLen = 8;
 
if(ProbeLen < MinMatchLen)
	return(false);

m_ProbeLen = ProbeLen;

m_ProbeEntryID = ProbeID;
m_TargEntryID = EntryID;

m_ChkSelfMatches = bChkSelfMatches;

m_Modes  = Modes;
m_CurMinMatchLen = MinMatchLen;			// minimal required match length
m_CurMaxNumMatches = MaxNumMatchesPPS;	// only look for this many matches


m_pRsltsFileHandler = RsltsFileHandler;

if(LoadSfxHeader(EntryID) == NULL)
	return(false);

SfxBlockID = 1;

if(LoadSfxSeq(EntryID) == NULL)
	return(false);

bComplemented = false;
bReversed = false;

if(LoadSfxArray(EntryID) == NULL)
	return(false);

for(CurMode = 0; CurMode < 4; CurMode++)
	{
		switch(CurMode) {
		case 0:				// (sense) forward p:agctg against t:agctg
			if(!(Modes & 0x01))
				continue;

			if(bComplemented)
				{
				CSeqTrans::ComplementStrand(ProbeLen,pProbe);
				bComplemented = false;
				}
			if(bReversed)
				{
				CSeqTrans::ReverseSeq(ProbeLen,pProbe);
				bReversed = false;
				}
			break;

		case 1:				// forward complement p:tcgac against t:agctg 
			if(!(Modes & 0x02))
				continue;
			if(!bComplemented)
				{
				CSeqTrans::ComplementStrand(ProbeLen,pProbe);
				bComplemented = true;
				}
			if(bReversed)
				{
				CSeqTrans::ReverseSeq(ProbeLen,pProbe);
				bReversed = false;
				}
			break;

		case 2:				// (antisense) reverse complement p:cagct against t:agctg
			if(!(Modes & 0x04))
				continue;

			if(!bComplemented)
				{
				CSeqTrans::ComplementStrand(ProbeLen,pProbe);
				bComplemented = true;
				}
			if(!bReversed)
				{
				CSeqTrans::ReverseSeq(ProbeLen,pProbe);
				bReversed = true;
				}
			break;
			
		case 3:				// reverse p:gtcga against t:agctg
			if(!(Modes & 0x08))
				continue;
			if(!bReversed)
				{
				CSeqTrans::ReverseSeq(ProbeLen,pProbe);
				bReversed = true;
				}
			if(bComplemented)
				{
				CSeqTrans::ComplementStrand(ProbeLen,pProbe);
				bComplemented = false;
				}
			break;
		}

	m_CurMatchMode = CurMode;
	TargStart = 0;
	ProbeOfs = 0;		
	ProbeLenRemaining = ProbeLen;

		
	while(ProbeLenRemaining >= MinMatchLen) 
		{

		bool bNotLocated= !IsLocated(pProbe,				// sequence to locate
						ProbeOfs,				// offset in pProbe at which to match from
						MinMatchLen,			// minimum match length
						m_pSfxSeq,				// target sequence
						m_pSfxArray,			// target sequence suffix array
						TargStart,				// position in pTarg (0..n) corresponding to start of suffix array
						m_CurSfxHeader.SeqLen);

		if(bNotLocated)
			(*m_pRsltsFileHandler)(eEMRAddMatch,ProbeOfs,0,MinMatchLen,MinMatchLen,CurMode);
			
		ProbeOfs+=1;
		ProbeLenRemaining -= 1;
		}
	}


if(bReversed)
	CSeqTrans::ReverseSeq(ProbeLen,pProbe);
if(bComplemented)
	CSeqTrans::ComplementStrand(ProbeLen,pProbe);
return(true);
}

// IsLocated
// Returns true if specified probe sequence of minimal ProbeLen is located in target sequence 
bool						
CSfxArray::IsLocated( etSeqBase *pProbe,	// sequence containing probe
				  unsigned int ProbePsn,	// start position in probe sequence at which to search
				  unsigned int MinNoMatchLen,	// minimal probe length 
				  etSeqBase *pTarg,			// target sequence
				  unsigned int *pSfxArray,	// target sequence suffix array
				  unsigned int TargStart,   // position in pTarg (0..n) corresponding to start of suffix array
				  unsigned int SfxLen)		// number of suffixs in pSfxArray
{
etSeqBase *pEl1;
etSeqBase *pEl2;
etSeqBase *pSeq2Loc;

etSeqBase b1;
etSeqBase b2;
int CurMatchLen;
int Max2Compare;
int NumMatches;
int SfxLo;
int SfxHi;
int TargPsn;

pSeq2Loc = &pProbe[ProbePsn];

if((*pSeq2Loc & ~cRptMskFlg) >= eBaseN)
	return(true);

NumMatches = 0;
SfxLo = 0;
SfxHi = SfxLen - 1;
do {
	TargPsn = (SfxLo + SfxHi) / 2;
	pEl2 = &pTarg[pSfxArray[TargPsn + TargStart]];
	Max2Compare = MinNoMatchLen;
	CurMatchLen = 0;
	pEl1 = pSeq2Loc;
	while(Max2Compare-- && (b1 = (*pEl1++ & ~cRptMskFlg)) == (b2 = (*pEl2++ & ~cRptMskFlg)) && b1 <= eBaseT)
		CurMatchLen++;
	if(CurMatchLen == (int)MinNoMatchLen)
		return(true);
	if(b1 > eBaseT)
		return(false);

	if(b2 > eBaseT)
		SfxHi = TargPsn-1;
	else
		{
		if(b1 < b2)
			SfxHi = TargPsn - 1;
		else
			SfxLo = TargPsn+1;
		}
	}
while(SfxHi >= SfxLo);
return(false);
}



// LocateNearExacts
// Returns count of all nearly exact maximally extended matching sequences
// These matching sequences must contain an exactly matching core of MinCoreLen and
// left/right flanking sequences which are maximally extended with at most LeftMaxMismatches/RightMaxMismatches
// The flanking sequences are extended for at most LeftMaxExtend and RightMaxExtend bases.
int						
CSfxArray::LocateNearExacts(sAMProbe *pProbe,	// probe
 				  unsigned int ProbeLen,	    // probe length
				  etSeqBase *pProbeSeq,			// probe sequence	
				  unsigned int MinCoreLen,	    // minimal core length from which to extend left + right
				  unsigned int LeftMaxExtend,   // maximum length to extend left
				  unsigned int RightMaxExtend,  // maximum length to extend right
  				  unsigned int LeftMaxMismatches,  // total mismatches allowed in left flank
   				  unsigned int RightMaxMismatches, // total mismatches allowed in right flank
				  unsigned int MinHitLen,	       // putative hits must be of at least this length
				  
				  // following is call back on all putative hits of probe into suffix entry, allows caller to filter putative hits
				  int (* FilterCore)(sAMProbe *pProbe,				// probe
										unsigned int ProbeOfs,		// probe offset at which probe hit starts 
										unsigned int HitLen,		// number of bases hit
										unsigned int ProbeCoreOfs,	// probe offset at which exact core starts
										unsigned int ProbeCoreLen,	// length of exact core
										unsigned int NumLeftMismatches, // number of mismatches in left flank
										unsigned int NumRightMismatches, // number of mismatches in right flank
										unsigned int TargetEntryID,	// target suffix array entry 
										unsigned int TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]
				  
				 // following is call back on all accepted hits of probe into suffix entry
				  int (* ProcessCore)(sAMProbe *pProbe,			// probe
										unsigned int ProbeOfs,		// probe offset at which hit starts 
										unsigned int HitLen,		// number of bases hit
										unsigned int TargEntryID,	// target suffix array entry 
										unsigned int TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]
				unsigned int TargEntryID,							// which suffix entry to match probe on
				etRPTMasking RPTMasking,							// how to interpret probe bases with cRptMskFlg set
				unsigned int MaskPercent)							// if matches have more than this percent repeats then slough
{
int Rslt;
int TargPsn;	// index in pSfxArray of current suffix sequence being processed
int StartProbeOfs;
unsigned int TargOfs;
unsigned int CurCoreOfs;
unsigned int CurProbeOfs;
unsigned int CurMaxMatchLen;
unsigned int CurMatchLen;
unsigned int LeftCoreStart;
unsigned int ProbeCoreLen;
unsigned int MinSeedLen;
unsigned int DeltaProbeOfs;

int NumProcessed;	// total number of cores after filtering
bool bTerm;
int Idx;

etSeqBase *pEl1;
etSeqBase *pEl2;
etSeqBase *pEl1r;
etSeqBase *pEl2r;

etSeqBase *pSeq2Loc;

etSeqBase *pTarg;			// target sequence
unsigned int *pSfxArray;	// target sequence suffix array
unsigned int TargStart = 0;   // position in pTarg (0..n) corresponding to start of suffix array
unsigned int SfxLen;		// number of suffixs in pSfxArray

etSeqBase b1;
etSeqBase b2;

int Max2Compare;

unsigned int LeftMismatches;
unsigned int RightMismatches;


if(MinCoreLen > ProbeLen)
	return(-1);
if(ProbeLen > (unsigned int)pProbe->ProbeLen)
	return(-1);

if(LoadSfxHeader(TargEntryID) == NULL)
	return(-1);

if(LoadSfxSeq(TargEntryID) == NULL)
	return(-1);

if(LoadSfxArray(TargEntryID) == NULL)
	return(-1);

if(MinCoreLen <= 12)
	{
	MinSeedLen = MinCoreLen;
	DeltaProbeOfs = 1;
	}
else		// else minhitlen must be at >= 13
	{
	if(MinCoreLen <= 49)
		{
		DeltaProbeOfs = MinCoreLen / 4;
		MinSeedLen = MinCoreLen - DeltaProbeOfs;
		if(MinSeedLen < 12)	// can only be if minhitlen < 16
			{
			MinSeedLen = 12;
			DeltaProbeOfs = MinCoreLen - MinSeedLen;
			}
		}
	else	// else minhitlen must be at >= 50 so MinSeedLen will be at least 25
		{
		DeltaProbeOfs = MinCoreLen / 2;
		MinSeedLen = MinCoreLen - DeltaProbeOfs;
		}
	}

m_NumPutCores = 0;			// clears out any previous putative cores
pTarg = m_pSfxSeq;
pSfxArray = m_pSfxArray;
SfxLen = m_CurSfxHeader.SeqLen;
NumProcessed = 0;
StartProbeOfs = 0;
CurProbeOfs = 0;
bTerm = false;

do {
	pSeq2Loc = &pProbeSeq[CurProbeOfs];			// initial start

	// putative cores which contain any bases other than canonical (a,c,g or t) are simply sloughed
	if(FiltOutRepeats(eRPTHignor,0,pSeq2Loc,MinSeedLen))
		continue;

	// locate first core of at least MinSeedLen which exactly matches
	if((TargPsn = LocateFirstExact(pSeq2Loc,MinSeedLen,pTarg,pSfxArray,TargStart,0,SfxLen-1)) < 0)
		continue;

	do  {
		TargOfs = (unsigned int)((pEl2 = &pTarg[pSfxArray[TargPsn + TargStart]]) - pTarg);
		pEl1r = &pProbeSeq[CurProbeOfs + MinSeedLen];
		pEl2r = pEl2 + MinCoreLen;

		// firstly try to extend left
		ProbeCoreLen = CurMatchLen = MinSeedLen;
		CurCoreOfs = CurProbeOfs;
		LeftCoreStart = CurProbeOfs;
		LeftMismatches = 0;
		if(CurCoreOfs > 0 && TargOfs > 0)
			{
			pEl1 = &pProbeSeq[CurCoreOfs - 1];
			pEl2 -= 1;
			do	{
				b1 = (*pEl1-- & ~cRptMskFlg);
				b2 = (*pEl2-- & ~cRptMskFlg);
				if(b1 > eBaseT || b2 > eBaseT)			// can't handle non-canonical
					break;
				if((b1 != b2) && ++LeftMismatches > LeftMaxMismatches)
					{
					LeftMismatches -= 1;
					break;
					}
				if(b1 == b2 && !LeftMismatches)
					{
					ProbeCoreLen++;
					LeftCoreStart--;
					}
				CurMatchLen++;
				CurCoreOfs -= 1;
				TargOfs -= 1;
				}
			while(CurCoreOfs > 0 && TargOfs > 0);
			}

				// after extending left then try to extend right
		CurMaxMatchLen = ProbeLen - CurCoreOfs;
		RightMismatches = 0;
		while(CurMatchLen < (int)CurMaxMatchLen)
			{
			b1 = (*pEl1r++ & ~cRptMskFlg);
			b2 = (*pEl2r++ & ~cRptMskFlg);
			if(b1 > eBaseT || b2 > eBaseT)			// can't handle non-canonical
				break;
			if((b1 != b2) && ++RightMismatches > RightMaxMismatches)
				{
				RightMismatches -= 1;
				break;
				}
			if(b1 == b2 && !RightMismatches)
				ProbeCoreLen += 1;
			CurMatchLen++;
			}

		if(ProbeCoreLen >= MinCoreLen && CurMatchLen >= MinHitLen)
			{
				// apply filtering for repeats and non-canonical bases
			if(!FiltOutRepeats(RPTMasking,MaskPercent,&pProbeSeq[StartProbeOfs],CurMatchLen))
				{			
					// check if caller wants to filter this putative approximate match out
				if(FilterCore != NULL)
					{
					if((Rslt = (*FilterCore)(pProbe,StartProbeOfs,
									CurMatchLen,
									LeftCoreStart,
									ProbeCoreLen,
									LeftMismatches,RightMismatches,
									TargEntryID,
									TargOfs,
									m_pSfxSeq)) < 0)
							return(Rslt);
					if(!Rslt)
						{
						bTerm = true;
						continue;
						}
					}
				else
					Rslt = 1;
				if(Rslt > 0 && (Rslt =AddPutCore(StartProbeOfs,	// probe offset at which hit starts
							CurMatchLen,						// number of bases hit
							TargOfs)) < 1)						// target loci hit
					return(Rslt);
				}
			}

		// any more suffixes to check
		if(++TargPsn >= (int)SfxLen)
			break;
		// check if probe matches next suffix
		Max2Compare = MinSeedLen;
		CurMatchLen = 0;
		pEl2 = &pTarg[pSfxArray[TargPsn + TargStart]];
		pEl1 = pSeq2Loc;
		while(Max2Compare-- && (b1 = (*pEl1++ & ~cRptMskFlg)) == (b2 = (*pEl2++ & ~cRptMskFlg)) && b1 < eBaseN && b2 < eBaseN)
			CurMatchLen++;
		}
	while(!bTerm && CurMatchLen >= (int)MinSeedLen);
	}
while(!bTerm && (CurProbeOfs += DeltaProbeOfs) <= (ProbeLen - MinSeedLen));

NumProcessed = FiltPutCores();
if(NumProcessed <= 0)
	return(NumProcessed);
tsPutCore *pCore = m_pPutCores;
int NumAccepted = 0;
for(Idx = 0; Idx < NumProcessed; Idx++,pCore++)
	{
	if(pCore->HitLen)
		{
		if((Rslt = (*ProcessCore)(pProbe,pCore->ProbeOfs,pCore->HitLen,TargEntryID,pCore->TargOfs,m_pSfxSeq))<=0)
			return(Rslt < 0 ? Rslt : NumAccepted);
		NumAccepted += 1;
		}
	}
return(NumAccepted);
}

int
CSfxArray::AddPutCore(unsigned int ProbeOfs,	// probe offset at which hit starts 
	unsigned int HitLen,					// number of bases hit
	unsigned int TargOfs)					// target loci hit	
{
tsPutCore *pCore;
int CoreIdx;
int ProbeEndOfs = ProbeOfs + HitLen - 1;
int TargEndOfs = TargOfs + HitLen - 1;

// iterate cores to determine if core is contained within existing core, or contains an existing core
// if contained within an existing then don't add, if containing then replace existing 
if(m_pPutCores != NULL && m_NumPutCores)
	{
	pCore = m_pPutCores;
	for(CoreIdx = 0; CoreIdx < m_NumPutCores; CoreIdx++,pCore++)
		{
		if(pCore->HitLen == 0)
			continue;
		if(pCore->ProbeOfs <= (UINT32)ProbeOfs && (pCore->ProbeOfs + pCore->HitLen - 1) >= (UINT32)ProbeEndOfs)
			{
			if(pCore->TargOfs <= TargOfs && (pCore->TargOfs + pCore->HitLen - 1) >= (UINT32)TargEndOfs)
				return(m_NumPutCores);
			}

		if(pCore->ProbeOfs >= (UINT32)ProbeOfs && (pCore->ProbeOfs + pCore->HitLen - 1) <= (UINT32)ProbeEndOfs)
			{
			if(pCore->TargOfs >= TargOfs && (pCore->TargOfs + pCore->HitLen - 1) <= (UINT32)TargEndOfs)
				{
				pCore->ProbeOfs = ProbeOfs;
				pCore->HitLen = HitLen;
				pCore->TargOfs = TargOfs;
				return(m_NumPutCores);
				}
			}
		}
	}

// core not contained within any other core so add
if(m_pPutCores == NULL || m_NumPutCores >= m_MaxPutCores)
	{
	if((pCore = new tsPutCore [cAllocPutCores + m_MaxPutCores])==NULL)
		return(eBSFerrMem);

	if(m_pPutCores != NULL)
		{
		memmove(pCore,m_pPutCores,m_NumPutCores * sizeof(tsPutCore));
		delete m_pPutCores;
		m_MaxPutCores += cAllocPutCores;
		}
	else
		{
		m_MaxPutCores = cAllocPutCores;
		m_NumPutCores = 0;
		}
	m_pPutCores = pCore;
	}
pCore = &m_pPutCores[m_NumPutCores++];
pCore->ProbeOfs = ProbeOfs;
pCore->HitLen = HitLen;
pCore->TargOfs = TargOfs;
return(m_NumPutCores);
}

// FiltContainedPutCores
// Marks for deletion (sets core length to 0) any cores starting from NxtProbeIdx which are 
// internally contained within cores starting at CurProbeIdx
int								// returns number of cores to be deleted
CSfxArray::FiltContainedPutCores(int CurProbeIdx, 
								 int NxtProbeIdx)
{
int Num2Mark;

unsigned int CurStartOfs;
unsigned int CurEndOfs;
unsigned int HitLen;
unsigned int CurTargOfs;
unsigned int CurTargEndOfs;

tsPutCore *pCurCore;
tsPutCore *pNxtCore;
pCurCore = &m_pPutCores[CurProbeIdx];

// nothing to do if probe already marked for removal?
if(!(HitLen = pCurCore->HitLen))
	return(0);

CurStartOfs = pCurCore->ProbeOfs;
CurEndOfs = CurStartOfs + pCurCore->HitLen - 1;
CurTargOfs = pCurCore->TargOfs;
CurTargEndOfs = CurTargOfs + pCurCore->HitLen - 1;
pNxtCore = &m_pPutCores[NxtProbeIdx];
Num2Mark = 0;
do	{
	if(pNxtCore->HitLen == 0)				// if target already marked then skip to next
		{
		pNxtCore+=1;
		continue;
		}

	if(pNxtCore->ProbeOfs > CurEndOfs)		// if next starts after current ends then subsequent probes can't be contained
		break;


	if(CurStartOfs <= pNxtCore->ProbeOfs &&
		CurEndOfs >= (pNxtCore->ProbeOfs + pNxtCore->HitLen - 1) && // if target core is contained within probe core
		CurTargOfs <= pNxtCore->TargOfs &&							// and target loci is contained
		CurTargEndOfs >= (pNxtCore->TargOfs + pNxtCore->HitLen - 1))
		{
		pNxtCore->HitLen = 0;	// mark target core to be removed
		Num2Mark++;
		}

	pNxtCore+=1;
	}
while(++NxtProbeIdx < m_NumPutCores);

return(Num2Mark);
}

// FiltOverlapPutCores
// Marks for deletion the shorter of any cores which overlap both on the CurProbeIdx and NxtProbeIdx..m_NumPutCores
int							// returns number of cores marked for deletion
CSfxArray::FiltOverlapPutCores(int CurProbeIdx, int NxtProbeIdx)
{
int Num2Mark;

unsigned int CurStartOfs;
unsigned int CurEndOfs;
unsigned int HitLen;
unsigned int CurTargOfs;
tsPutCore *pCurCore;
tsPutCore *pNxtCore;
pCurCore = &m_pPutCores[CurProbeIdx];
if(!(HitLen = pCurCore->HitLen))
	return(0);
CurStartOfs = pCurCore->ProbeOfs;
CurEndOfs = CurStartOfs + HitLen - 1;
CurTargOfs = pCurCore->TargOfs;
pNxtCore = &m_pPutCores[NxtProbeIdx];
Num2Mark = 0;

do	{
	if(pNxtCore->HitLen == 0)				// on to next if already marked for deletion
		{
		pNxtCore+=1;
		continue;
		}

	if(pNxtCore->ProbeOfs > CurEndOfs)		// if next core start is after current core ends then can't be any overlap
		break;


	if((CurStartOfs <= pNxtCore->ProbeOfs) && // already know that pNxtCore->ProbeOfs <= CurEndOfs so here we are checking for any overlap on the probe cores
		((pNxtCore->TargOfs + pNxtCore->HitLen - 1) >= CurTargOfs) && // assuming probe core overlaps then here we are checking for overlap on the target
		(pNxtCore->TargOfs <= (CurTargOfs + HitLen - 1)))
		{									  // have an overlap on both the probe core and target sequence loci 
		// retain longest core, shortest to be deleted
		if(pNxtCore->HitLen > HitLen)
			{
			pCurCore->HitLen = 0;
			Num2Mark++;
			break;
			}
		else
			pNxtCore->HitLen = 0;	// mark target core to be removed
		Num2Mark++;
		}
	pNxtCore+=1;
	}
while(++NxtProbeIdx < m_NumPutCores);

return(Num2Mark);
}

// FiltPutCores
// Filter cores for contained hits from probes onto target sequence, optionally also remove overlaps
int
CSfxArray::FiltPutCores(bool bRemoveOverlaps)
{
int MarkIdx;
int NumMarked;
bool bResort;			// will be set true if any cores marked and resort is required
tsPutCore *pCore;

if(m_NumPutCores <= 1)
	return(m_NumPutCores);

qsort(m_pPutCores,m_NumPutCores,sizeof(tsPutCore),SortPutCores);

// mark any cores contained within other cores
bResort = false;
do {
	NumMarked = 0;
	for(MarkIdx = 0; MarkIdx < m_NumPutCores-1; MarkIdx++)
		NumMarked += FiltContainedPutCores(MarkIdx,MarkIdx+1);
	if(NumMarked)
		bResort = true;
	}
while(NumMarked > 0);

if(bRemoveOverlaps)	// are overlaps to be removed?
	{
	if(bResort)
		{
		// sort to put marked cores last and then update m_NumPutCores with number of cores not marked
		qsort(m_pPutCores,m_NumPutCores,sizeof(tsPutCore),SortPutCores);
		pCore = m_pPutCores;
		for(MarkIdx = 0; MarkIdx < m_NumPutCores; MarkIdx++,pCore++)
			if(pCore->HitLen == 0)
				break;
		m_NumPutCores = MarkIdx;
		}

	bResort = false;
	// mark any overlap cores
	do {
		NumMarked = 0;
		for(MarkIdx = 0; MarkIdx < m_NumPutCores-1; MarkIdx++)
			NumMarked += FiltOverlapPutCores(MarkIdx,MarkIdx+1);
		if(NumMarked)
			bResort = true;
		}
	while(NumMarked > 0);
	}

if(bResort)
	{
	// sort to put marked cores last and then update m_NumPutCores with number of cores not marked
	qsort(m_pPutCores,m_NumPutCores,sizeof(tsPutCore),SortPutCores);
	pCore = m_pPutCores;
	for(MarkIdx = 0; MarkIdx < m_NumPutCores; MarkIdx++,pCore++)
		if(pCore->HitLen == 0)
			break;
	m_NumPutCores = MarkIdx;
	}
return(m_NumPutCores);
}

// FiltOutRepeats
// Returns true if sequence is to be filtered out either because it contains a non-canonical base or the percentage of repeat masked bases is too high
bool
CSfxArray::FiltOutRepeats(etRPTMasking RPTMasking,	// how to interpret cRptMskFlg'd bases
						 unsigned int MaskPercent,	// filter out if percentage of repeats is above this percentage (0-100)
						 etSeqBase *pSeq,			// pts to sequence
						 int SeqLen)				// sequence length
{
etSeqBase Base;
unsigned int Rpts = 0;
unsigned int Percent;
int Len = SeqLen;
switch(RPTMasking) {
	case eRPTHignor:	// treat all bases as not a repeat but still check for non-canonical bases
		while(Len--)
			if((*pSeq++ & ~cRptMskFlg)  > eBaseT)
				return(true);
		return(false);

	case eRPTHmasked:	// treat any base cRptMskFlg as a repeat masked base
		while(Len--)
			{
			if(((Base = *pSeq++) & ~cRptMskFlg)  > eBaseT)
				return(true);
			if(Base & cRptMskFlg)
				Rpts++;
			}
		break;

	case eRPTHunmasked:		// treat any base without cRptMskFlg as a repeat masked base 
		while(Len--)
			{
			if(((Base = *pSeq++) & ~cRptMskFlg)  > eBaseT)
				return(true);
			if(!(Base & cRptMskFlg))
				Rpts++;
			}
		break;
	}

Percent = (100*Rpts) / SeqLen;
return(Percent >= MaskPercent ? true : false);
}

// LocateExacts
// Returns count of all exact matching sequences
// These matching sequences are exact (100% identity) and contain no mismatches
int	
CSfxArray::LocateExacts(sAMProbe *pProbe,		// probe
				  unsigned int ProbeLen,	    // probe sequence length
				  etSeqBase *pProbeSeq,			// probe sequence
				  unsigned int MinHitLen,	    // minimum target hit length required

  				  // following is call back on all putative hits of probe into suffix entry, allows caller to filter putative hits
				  int (* FilterCore)(sAMProbe *pProbe,				// probe
										unsigned int ProbeOfs,		// probe offset at which exact core starts 
										unsigned int ProbeCoreLen,	// length of exact core
										unsigned int TargetEntryID,	// target suffix array entry 
										unsigned int TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]

 				  // following is call back on all hits of probe into suffix entry
				  int (* ProcessCore)(sAMProbe *pProbe,			// probe
										unsigned int ProbeOfs,		// probe offset at which hit starts 
										unsigned int HitLen,		// number of bases hit
										unsigned int TargEntryID,	// target suffix array entry 
										unsigned int TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]
				unsigned int TargEntryID,							// which suffix entry to match probe on
				etRPTMasking RPTMasking,							// how to interpret probe bases with cRptMskFlg set
				unsigned int MaskPercent)							// if matches have more than this percent repeats then slough
{
int Rslt;
int StartProbeOfs;
int TargPsn;	// index in pSfxArray of current suffix sequence being processed
unsigned int TargOfs;
unsigned int CurCoreOfs;
unsigned int CurProbeOfs;
unsigned int CurMaxMatchLen;
unsigned int MaxMatchLen;
unsigned int CurMatchLen;
unsigned int MinSeedLen;
unsigned int DeltaProbeOfs;
int NumProcessed;	// total number of cores after filtering
tsPutCore *pCore;
bool bTerm;
int Idx;

etSeqBase *pEl1;
etSeqBase *pEl2;
etSeqBase *pEl1r;
etSeqBase *pEl2r;

etSeqBase *pSeq2Loc;

etSeqBase *pTarg;			// target sequence
unsigned int *pSfxArray;	// target sequence suffix array
unsigned int TargStart = 0;   // position in pTarg (0..n) corresponding to start of suffix array
unsigned int SfxLen;		// number of suffixs in pSfxArray

etSeqBase b1;
etSeqBase b2;

int Max2Compare;

if(MinHitLen > ProbeLen)
   return(0);

if(LoadSfxHeader(TargEntryID) == NULL)
	return(-1);

if(LoadSfxSeq(TargEntryID) == NULL)
	return(-1);

if(LoadSfxArray(TargEntryID) == NULL)
	return(-1);

m_NumPutCores = 0;			// clears out any previous putative cores
MaxMatchLen = 0;
pTarg = m_pSfxSeq;
pSfxArray = m_pSfxArray;
SfxLen = m_CurSfxHeader.SeqLen;
NumProcessed = 0;

if(MinHitLen <= 12)
	{
	MinSeedLen = MinHitLen;
	DeltaProbeOfs = 1;
	}
else		// else minhitlen must be at >= 13
	{
	if(MinHitLen <= 49)
		{
		DeltaProbeOfs = MinHitLen / 4;
		MinSeedLen = MinHitLen - DeltaProbeOfs;
		if(MinSeedLen < 12)	// can only be if minhitlen < 16
			{
			MinSeedLen = 12;
			DeltaProbeOfs = MinHitLen - MinSeedLen;
			}
		}
	else	// else minhitlen must be at >= 50 so MinSeedLen will be at least 25
		{
		DeltaProbeOfs = MinHitLen / 2;
		MinSeedLen = MinHitLen - DeltaProbeOfs;
		}
	}

StartProbeOfs = 0;
CurProbeOfs = 0;
bTerm = false;

do {
	pSeq2Loc = &pProbeSeq[CurProbeOfs];			// initial start
	// putative cores which contain any bases other than canonical (a,c,g or t) are simply sloughed
	if(FiltOutRepeats(eRPTHignor,0,pSeq2Loc,MinSeedLen))
		continue;

	// locate first core of at least MinSeedLen which exactly matches
	if((TargPsn = LocateFirstExact(pSeq2Loc,MinSeedLen,pTarg,pSfxArray,TargStart,0,SfxLen-1)) < 0)
		continue;

	// have at least one core exactly matching MinSeedLen which could be extended left+right to MinHitLen
	do  {
		
		TargOfs = (unsigned int)((pEl2 = &pTarg[pSfxArray[TargPsn + TargStart]]) - pTarg);
		pEl1r = &pProbeSeq[CurProbeOfs + MinSeedLen];
		pEl2r = pEl2 + MinSeedLen;

		// firstly try to extend left
		CurMatchLen = MinSeedLen;
		CurCoreOfs = CurProbeOfs;
		if(CurCoreOfs > 0 && TargOfs > 0)
			{
			pEl1 = &pProbeSeq[CurCoreOfs - 1];
			pEl2 -= 1;
			do	{
				b1 = (*pEl1-- & ~cRptMskFlg);
				b2 = (*pEl2-- & ~cRptMskFlg);
				if((b1 != b2) || (b1 > eBaseT || b2 > eBaseT))
					break;
				CurMatchLen++;
				CurCoreOfs -= 1;
				TargOfs -= 1;
				}
			while(CurCoreOfs > 0 && TargOfs > 0);
			}

		// after extending left then try to extend right
		CurMaxMatchLen = pProbe->ProbeLen - CurCoreOfs;
		while(CurMatchLen < (int)CurMaxMatchLen)
			{
			b1 = (*pEl1r++ & ~cRptMskFlg);
			b2 = (*pEl2r++ & ~cRptMskFlg);
			if((b1 != b2) || (b1 > eBaseT || b2 > eBaseT))
				break;
			CurMatchLen++;
			}

		if(CurMatchLen > MaxMatchLen)
			MaxMatchLen = CurMatchLen;

		if(CurMatchLen >= MinHitLen)
			{
			// apply filtering for repeats and non-canonical bases before calling AddPutCore
			if(!FiltOutRepeats(RPTMasking,MaskPercent,&pProbeSeq[CurCoreOfs],CurMatchLen))
				{	
				// check if caller wants to filter this putative approximate match out
				if(FilterCore != NULL)
					{
					if((Rslt = (*FilterCore)(pProbe,CurCoreOfs,
								CurMatchLen,
								TargEntryID,
								TargOfs,
								m_pSfxSeq)) < 0)
						return(Rslt);
					if(!Rslt)
						{
						bTerm = true;
						continue;
						}
					}
				else
					Rslt = 1;
				if(Rslt > 0 && (Rslt =AddPutCore(CurCoreOfs,	// probe offset at which hit starts
							CurMatchLen,						// number of bases hit
							TargOfs)) < 1)						// target loci hit
					return(Rslt);
				}
			}

		// any more suffixes to check
		if(++TargPsn >= (int)SfxLen)
			break;

		// check if probe matches next suffix
		Max2Compare = MinSeedLen;
		CurMatchLen = 0;
		pEl2 = &pTarg[pSfxArray[TargPsn + TargStart]];
		pEl1 = pSeq2Loc;
		while(Max2Compare-- && (b1 = (*pEl1++ & ~cRptMskFlg)) == (b2 = (*pEl2++ & ~cRptMskFlg)) && b1 <= eBaseT)
			CurMatchLen++;
		}
	while(!bTerm && CurMatchLen >= (int)MinSeedLen);
	}
while(!bTerm && (CurProbeOfs += DeltaProbeOfs) <= (ProbeLen - MinSeedLen));

NumProcessed = FiltPutCores();
if(NumProcessed <= 0)
	return(NumProcessed);
pCore = m_pPutCores;
int NumAccepted = 0;
for(Idx = 0; Idx < NumProcessed; Idx++,pCore++)
	{
	if(pCore->HitLen)
		{
		if((Rslt = (*ProcessCore)(pProbe,pCore->ProbeOfs,pCore->HitLen,TargEntryID,pCore->TargOfs,m_pSfxSeq))<=0)
			return(Rslt < 0 ? Rslt : NumAccepted);
		NumAccepted += 1;
		}
	}
return(NumAccepted);
}

int
CSfxArray::SetTargEntry(unsigned int TargEntryID)		// which suffix entry (chromosome) to process with IterateExacts
{
if(LoadSfxHeader(TargEntryID) == NULL)
	return(-1);

if(LoadSfxSeq(TargEntryID,true) == NULL)
	return(-1);

if(LoadSfxArray(TargEntryID) == NULL)
	return(-1);
return(TargEntryID);
}

int						// returned loci (0..n) or -1 if errors
CSfxArray::IterateExacts(etSeqBase *pProbeSeq,			// probe
 						 unsigned int ProbeLen,			// probe length
						 unsigned int PrevHitIdx,		// 0 if starting new sequence, otherwise set to return value of previous successful iteration return
						 unsigned int *pHitLoci,		// if match then where to return loci
	 				     int MaxNumExacts)				// only interested in at most this number of exact matches (0 if all exacts are of interest)
{
int TargPsn;

etSeqBase *pEl1;
etSeqBase *pEl2;

etSeqBase *pTarg;			// target sequence
unsigned int *pSfxArray;	// target sequence suffix array
unsigned int TargStart = 0;   // position in pTarg (0..n) corresponding to start of suffix array
unsigned int SfxLen;		// number of suffixs in pSfxArray

*pHitLoci = 0;
	   // any more suffixes to check
if(PrevHitIdx >= (int)m_CurSfxHeader.SeqLen)
	return(0);

pTarg = m_pSfxSeq;
pSfxArray = m_pSfxArray;
SfxLen = m_CurSfxHeader.SeqLen;

if(!PrevHitIdx)
	{
	// locate first core of at least MinSeedLen which exactly matches
	if((TargPsn = LocateFirstExact(pProbeSeq,ProbeLen,pTarg,pSfxArray,0,0,SfxLen-1,MaxNumExacts)) < 0)
		return(0);	// no match
	*pHitLoci = pSfxArray[TargPsn];
	return(TargPsn + 1);
	}

  // check if probe matches next suffix
pEl2 = &pTarg[pSfxArray[PrevHitIdx]];
pEl1 = pProbeSeq;

if(!memcmp(pEl1,pEl2,ProbeLen))
	{
	*pHitLoci = pSfxArray[PrevHitIdx];
	return(PrevHitIdx+1);
	}

return(0);
}

int						// < 0 if errors, 0 if no matches, 1 if a unique match, 2 if multiple matches
CSfxArray::LocateApproxUniques(int AccumReadHits,				// how many reads have already been matched 
						 etSeqBase *pProbeSeq,			// probe
 						 etSeqBase *pPatternSeq,		// contains pattern to match with, etBaseN represents wildcard bases to match against any in the target
														// will be updated on return with etBaseN's changed to actual subsitution bases  - NULL if all bases can be wildcarded
						 unsigned int ProbeLen,			// probe, and also pattern, length
						 int MinSubCnt,					// minimum allowed pattern wildcard matches
						 int MaxSubCnt,					// maximum allowed pattern wildcard matches
						 unsigned int *pHitLoci,		// if unique match then where to return loci
						 int *pSubCnt)					// used to return number of pattern wildcard substitutions actually required
{
int CurProbeSegOfs;
int CurProbeSegLen;
int SegIdx;
int IterCnt;

etSeqBase *pProbeBase;
etSeqBase *pPatternBase;
etSeqBase *pTargBase;

int TargIdx;

int CurSubCnt;
int PatIdx;
int TargSeqLeftIdx;
int ProbeSeqLeftIdx;
int TargMatchLen;
int PutativeTargLoci;
int CurInstances;
int CurMaxSubCnt;

etSeqBase *pTarg;			// target sequence
unsigned int *pSfxArray;	// target sequence suffix array
unsigned int SfxLen;		// number of suffixs in pSfxArray

*pHitLoci = 0;

pTarg = m_pSfxSeq;
pSfxArray = m_pSfxArray;
SfxLen = m_CurSfxHeader.SeqLen;

CurProbeSegOfs = 0;
CurInstances = AccumReadHits;
CurMaxSubCnt = MaxSubCnt;
for(SegIdx = 0; SegIdx <= MaxSubCnt; SegIdx++,CurProbeSegOfs += CurProbeSegLen)
	{
	CurProbeSegLen = (ProbeLen - CurProbeSegOfs) / (1 + (MaxSubCnt - SegIdx));
	if((TargIdx = LocateFirstExact(&pProbeSeq[CurProbeSegOfs],CurProbeSegLen,pTarg,pSfxArray,0,0,SfxLen-1,0)) < 0)
		{
		// no segment match, onto next segment
		continue;
		}

	IterCnt = 0;
	while(1)
		{
		if(IterCnt++) {				// if iterating subsequent (to LocateFirstExact) targets
			// ensure not about to iterate past end of suffix array!
			if((TargIdx + 1) >= (int)m_CurSfxHeader.SeqLen || (pSfxArray[TargIdx+1] + CurProbeSegLen) > m_CurSfxHeader.SeqLen) 
				break;
			pTargBase = &pTarg[pSfxArray[TargIdx+1]]; // then check that target is still matching
			pProbeBase = &pProbeSeq[CurProbeSegOfs];
			if(memcmp(pProbeBase,pTargBase,CurProbeSegLen))
				break;					// break if target no longer matches
			TargIdx += 1;
			}

		if(!SegIdx)			// if leftmost segment then check pattern whilst extending right
			{
			TargSeqLeftIdx = pSfxArray[TargIdx] + CurProbeSegLen;
			ProbeSeqLeftIdx = CurProbeSegLen;
			TargMatchLen = ProbeLen - CurProbeSegLen;
			PutativeTargLoci = pSfxArray[TargIdx];
			}
		else
			if(SegIdx == MaxSubCnt)	// if rightmost segment then check pattern upto start of segment
				{
				TargMatchLen = ProbeLen - CurProbeSegLen;
				TargSeqLeftIdx = pSfxArray[TargIdx] - TargMatchLen;
				ProbeSeqLeftIdx = 0;
				PutativeTargLoci = TargSeqLeftIdx;
 				}
			else
				{
				TargMatchLen = ProbeLen; 
				TargSeqLeftIdx = pSfxArray[TargIdx] - CurProbeSegOfs;
				ProbeSeqLeftIdx = 0;
				PutativeTargLoci = TargSeqLeftIdx;
				}

		// ensure comparisons are within start/end range of target sequence/assembly
		if(!TargMatchLen || TargSeqLeftIdx < 0 || (TargSeqLeftIdx + TargMatchLen) > (int)m_CurSfxHeader.SeqLen)
			break;

		// now do the matching
		pTargBase = &pTarg[TargSeqLeftIdx]; // then check that target is still matching
		pProbeBase = &pProbeSeq[ProbeSeqLeftIdx];
		pPatternBase = pPatternSeq == NULL ? NULL : &pPatternSeq[ProbeSeqLeftIdx];
		CurSubCnt = 0;
		for(PatIdx = 0; PatIdx < TargMatchLen; PatIdx++,pTargBase++,pProbeBase++)
			{
			if(*pProbeBase != *pTargBase)
				{
				if(pPatternBase != NULL && *pPatternBase != eBaseN)
					break;
				if(++CurSubCnt > MaxSubCnt)
					break;
				}
			if(pPatternBase != NULL)
					pPatternBase += 1;
			}
		if(PatIdx != TargMatchLen || MinSubCnt > CurSubCnt)
			continue;

		if(CurSubCnt < CurMaxSubCnt)
			{
			CurMaxSubCnt = CurSubCnt;
			CurInstances = AccumReadHits;
			}
		CurInstances += 1;
		if(CurInstances == 1)
			*pHitLoci = PutativeTargLoci;
		else
			if(MinSubCnt == CurMaxSubCnt)
				break;
		}
	if(CurInstances > 1 && MinSubCnt == CurMaxSubCnt)
		break;
	}
*pSubCnt = CurMaxSubCnt;
return(CurInstances);
}


int
CSfxArray::CmpProbeTarg(etSeqBase *pEl1,etSeqBase *pEl2,int Len)
{
int Psn;
for(Psn=0; Psn < Len; Psn++,pEl1++,pEl2++)
	{
	if(*pEl1 > *pEl2 || *pEl2 == eBaseEOS)
		return(1);
	if(*pEl1 < *pEl2)
		return(-1);
	}
return(0);
}

int			// index in pSfxArray of exactly matching probe or -1 if no match					
CSfxArray::LocateFirstExact(etSeqBase *pProbe,  // pts to probe sequence
				  unsigned int ProbeLen,		// probe length to exactly match over
				  etSeqBase *pTarg,				// target sequence
				  unsigned int *pSfxArray,		// target sequence suffix array
				  unsigned int TargStart,		// position in pTarg (0..n) corresponding to start of suffix array
				  int SfxLo,					// low index in pSfxArray
				  int SfxHi,					// high index in pSfxArray
			     int MaxNumExacts)				// only interested in at most this number of exact matches (0 if all exacts are of interest)
{
etSeqBase *pEl1;
etSeqBase *pEl2;

int CmpRslt;

int TargPsn;
int LowPsn;

pEl1 = pProbe;
do {
	TargPsn = (SfxLo + SfxHi) / 2;
	pEl2 = &pTarg[pSfxArray[TargPsn + TargStart]];

	if(!(CmpRslt = CmpProbeTarg(pEl1,pEl2,ProbeLen))) // note - can't use simple memcmp() as target contains eBaseEOS
		{
		if(TargPsn == 0 || SfxLo == TargPsn)
			return(TargPsn);
		if(MaxNumExacts && !--MaxNumExacts)
			return(TargPsn);
		LowPsn = LocateFirstExact(pProbe,ProbeLen,pTarg,pSfxArray,TargStart,SfxLo,TargPsn - 1,MaxNumExacts);
		return(LowPsn < 0 ? TargPsn : LowPsn);
		}

	if(CmpRslt < 0)
		SfxHi = TargPsn - 1;
	else
		SfxLo = TargPsn+1;
	}
while(SfxHi >= SfxLo);

return(-1);	// unable to locate any instance of pProbe
}

// GetNMerFreq
// Counts under or over represented NMers of specified length
int	
CSfxArray::GetNMerFreq(unsigned int TargEntryID,		// which suffix entry
					  CSfxNMerFreq *pNMerFreq,			// pre-initialised NMerFreq with min/max freq cutoff etc
						  			// following is call back on all putative hits into suffix entry
					  bool		// false if instance hit to be sloughed
								(* Filter)(unsigned int NMerLen,		// number of bases hit
								unsigned int TargEntryID,	// target suffix array entry 
								unsigned int TargOfs,		// target loci hit	
								etSeqBase *pSeq))			// target sequence, hit starts at &pSeq[TargOfs]
{
unsigned int SfxIdx;
int NumMatches;
etSeqBase b1;
etSeqBase b2;
int CurMatchLen;
int Max2Compare;
int NMerLen;
etSeqBase *pEl;
etSeqBase *pEl1;
etSeqBase *pEl2;

etSeqBase *pTarg;			// target sequence
unsigned int *pSfxArray;	// target sequence suffix array
unsigned int SfxLen;		// number of suffixs in pSfxArray
unsigned int InitialNumInstances;

if(LoadSfxHeader(TargEntryID) == NULL)
	return(-1);

if(LoadSfxSeq(TargEntryID) == NULL)
	return(-1);

if(LoadSfxArray(TargEntryID) == NULL)
	return(-1);

if((NMerLen = pNMerFreq->m_NMerLen) < 1)
	return(-1);

InitialNumInstances = pNMerFreq->m_CurNumInstances;

pTarg = m_pSfxSeq;
pSfxArray = m_pSfxArray;
SfxLen = m_CurSfxHeader.SeqLen;
if(NMerLen > (int)SfxLen)
	return(0);

NumMatches = 0;
pEl = NULL;

for(SfxIdx = 0; SfxIdx <= (SfxLen - NMerLen); SfxIdx++)
	{
	if(!Filter(NMerLen,TargEntryID,pSfxArray[SfxIdx],pTarg))
		continue;
	if(!NumMatches)
		{
		pEl = &pTarg[pSfxArray[SfxIdx]];
		NumMatches = 1;
		continue;
		}

	pEl1 = pEl;
	pEl2 = &pTarg[pSfxArray[SfxIdx]];
	Max2Compare = NMerLen;
	CurMatchLen = 0;
	while(Max2Compare-- && (b1 = (*pEl1++ & ~cRptMskFlg)) == (b2 = (*pEl2++ & ~cRptMskFlg)) && b1 <= eBaseT)
		CurMatchLen++;

	if(CurMatchLen == NMerLen) // still an exact match?
		{
		NumMatches += 1;
		continue;
		}

	// mismatch
	pNMerFreq->Add(NumMatches,pEl);
	NumMatches = 0;
	SfxIdx -= 1;
	}

if(NumMatches > 0)
	pNMerFreq->Add(NumMatches,pEl);
return(pNMerFreq->m_CurNumInstances - InitialNumInstances);
}

// SortPutCores
// Sorts putative cores by:
// ascending probe offset
// decending hit length
// ascending target offset
// BUT if hit length is 0 (core marked for deletion) then that core will be sorted as being last irrespective of probe and target offsets
int CSfxArray::SortPutCores(const void *arg1, const void *arg2)
{
tsPutCore *pEl1 = (tsPutCore *)arg1;
tsPutCore *pEl2 = (tsPutCore *)arg2;

if(pEl1->HitLen == 0 || pEl2->HitLen == 0)
	{
	if(pEl1->HitLen)
		return(-1);
	if(pEl2->HitLen)
		return(1);
	return(0);
	}

// sort by ascending probe start loci
if(pEl1->ProbeOfs != pEl2->ProbeOfs)
	return(pEl1->ProbeOfs < pEl2->ProbeOfs ? -1 : 1);

// sort by descending hit length
if(pEl1->HitLen != pEl2->HitLen)
	return(pEl1->HitLen < pEl2->HitLen ? 1 : -1);

// sort by ascending target start loci
if(pEl1->TargOfs != pEl2->TargOfs)
	return(pEl1->TargOfs < pEl2->TargOfs ? -1 : 1);

return(0);
}


int CSfxArray::CompareExtdMatches( const void *arg1, const void *arg2 )
{
tsSfxMatch *pEl1 = (tsSfxMatch *)arg1;
tsSfxMatch *pEl2 = (tsSfxMatch *)arg2;
if(pEl1->ProbePsn == pEl2->ProbePsn )
	{
	if(pEl1->TargPsn < pEl2->TargPsn)
		return(-1);
	if(pEl1->TargPsn > pEl2->TargPsn)
		return(1);
	if(pEl1->MatchLen < pEl2->MatchLen)
		return(-1);
	if(pEl1->MatchLen > pEl2->MatchLen)
		return(1);
	return(0);
	}

return(pEl1->ProbePsn < pEl2->ProbePsn ? -1 : 1);
}
