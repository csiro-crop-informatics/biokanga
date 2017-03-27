/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "stdafx.h"

#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

etSeqBase CHashFile::m_Sense2AntiMap[7] = {eBaseT, eBaseG, eBaseC, eBaseA,eBaseN, eBaseUndef, eBaseInDel };

CHashFile::CHashFile(void)
{
m_pSeq = NULL;				// allocd memory to hold genomic sequence
m_pSlots = NULL;			// allocd memory to hold hash array slots
m_pHashes = NULL;			// allocd memory to hold sequence psns for each hash
m_pExtdMatches = NULL;
Reset(false);
}

CHashFile::~CHashFile(void)
{
Reset(false);
}

int
CHashFile::Reset(bool bFlush)
{
Close(bFlush);
if(m_pSeq != NULL)
	{
	delete m_pSeq;
	m_pSeq = NULL;
	}
if(m_pSlots != NULL)
	{
	delete m_pSlots;
	m_pSlots = NULL;
	}
if(m_pHashes != NULL)
	{
	delete m_pHashes;
	m_pHashes = NULL;
	}
if(m_pExtdMatches != NULL)
	{
	delete m_pExtdMatches;
	m_pExtdMatches = NULL;
	}

memset(&m_CurHashHeader,0,sizeof(m_CurHashHeader));

m_MaxSeqLen = 0;			// max number of bases that could be currenty held in m_pSeq
m_MaxNumSlots=0;			// max number of slots be currently in m_pSlots
m_MaxNumHashes=0;			// max number of hashes that could be currently held in pHashes

m_SeqEntryID = 0;
m_HashesEntryID = 0;
m_SlotsEntryID = 0;

m_NumExtdMatches = 0;

return(eBSFSuccess);
}

// Open()
// Opens specified hash file pszHashFile
// Option to create or truncate pszSeqFile
int CHashFile::Open(char *pszFile,
					  bool bCreate)
{
Reset(false);
return(CBioSeqFile::Open(pszFile,cBSFTypeHash,bCreate));
}



int 
CHashFile::Close(bool bFlush)
{
return(CBioSeqFile::Flush2Disk());
}


void
CHashFile::GenSenseHash(unsigned int *pHash,
					 etSeqBase *pSeqBase,
					 unsigned int Len)
{
unsigned int Hash = 0;
unsigned int Cnt;

if(Len > cMaxHashRange)
	Len = cMaxHashRange;

for(Cnt = 0; Cnt < Len; Cnt++)
	{
	Hash <<= 1;
	Hash |= *pSeqBase++ & 0x001;
	}
*pHash = Hash & m_HashMask;
return;
}

void
CHashFile::GenAntisenseHash(unsigned int *pHash,
					 etSeqBase *pSeqBase,
					 unsigned int Len)
{
unsigned int Hash = 0;
unsigned int Cnt;

if(Len > cMaxHashRange)
	Len = cMaxHashRange;

for(Cnt = 0; Cnt < Len; Cnt++)
	{
	Hash <<= 1;
	Hash |= m_Sense2AntiMap[*pSeqBase++] & 0x001;
	}
*pHash = Hash & m_HashMask;
return;
}

// SeqSenseMatches matches pProbe forward against pTarg forward
bool
CHashFile::SeqSenseMatches(int ProbeLen,etSeqBase *pProbe,etSeqBase *pTarg)
{
while(ProbeLen--)
	{
	if(*pProbe >= eBaseN || *pTarg >= eBaseN)
		return(false);
	if(*pProbe++ != *pTarg++)
		return(false);
	}
return(true);
}

// SeqSenseRevMatches matches pProbe forward against pTarg reverse
bool
CHashFile::SeqSenseRevMatches(int ProbeLen,etSeqBase *pProbe,etSeqBase *pTarg)
{
pTarg += ProbeLen - 1;
while(ProbeLen--)
	{
	if(*pProbe >= eBaseN || *pTarg >= eBaseN)
		return(false);
	if(*pProbe++ != *pTarg--)
		return(false);
	}
return(true);
}


// SeqAntisenseMatches matches  antisense pProbe forward against pTarg forward
bool
CHashFile::SeqAntisenseMatches(int ProbeLen,etSeqBase *pProbe,etSeqBase *pTarg)
{
while(ProbeLen--)
	{
	if(*pProbe >= eBaseN || *pTarg >= eBaseN)
		return(false);
	if(m_Sense2AntiMap[*pProbe++] != *pTarg++)
		return(false);
	}
return(true);
}

// SeqAntisenseRevMatches matches antisense pProbe forward against pTarg reverse
bool
CHashFile::SeqAntisenseRevMatches(int ProbeLen,etSeqBase *pProbe,etSeqBase *pTarg)
{
pTarg += ProbeLen - 1;
while(ProbeLen--)
	{
	if(*pProbe >= eBaseN || *pTarg >= eBaseN)
		return(false);
	if(m_Sense2AntiMap[*pProbe++] != *pTarg--)
		return(false);
	}
return(true);
}




unsigned int								// (EntryID) uniquely identifies this hashed sequence
CHashFile::AddEntry(bool bChkComplexity,	// true if low complexity sequences are to be filtered out
				char *pszSource,			// from where was this sequence originated
				char *pszDescription,		// descriptive text about this sequence
				etSeqBase   *pSeq,			// sequence to generate hashes over
				unsigned int SeqLen,		// sequence length
				unsigned int WindowSize)	// hashes are generated over this size non-overlapping window 
{
unsigned int SenseHash;
tHashID NxtHashID;
tHashID HashEntryID;
unsigned int SeqOfs;
unsigned int HashRange;
unsigned int MaxNumHashes;
unsigned int MaxEntries;
unsigned int ExactMatchCnt;

AddErrMsg("CHashFile::AddEntry","Hashing sequence of length %9.9d, window size %d from %s|%s...         ",SeqLen,WindowSize,pszSource,pszDescription); 
etSeqBase *pTargWindow;
etSeqBase *pProbeWindow;

if(WindowSize < cMinHashWinSize || WindowSize > cMaxHashWinSize)
	return(false);
if(SeqLen < WindowSize)
	return(false);

if(WindowSize < cMaxHashRange)
	HashRange = 1 << WindowSize;
else
	HashRange = 1 << cMaxHashRange;

m_HashMask = HashRange - 1;

tsHashArraySlot *pSlot;				// pts to array slot for hash
tsHashEntry *pCurHash;				// pts to just generated entry for hash  

tsHashEntry *pEntry;				// pts to entry being checked for repeats  


Reset(false);

if(m_MaxNumSlots < HashRange || m_pSlots == NULL)
	{
	if(m_pSlots != NULL)
		delete m_pSlots;
	m_pSlots = new tsHashArraySlot[HashRange];
	if(m_pSlots == NULL)
		{
		m_MaxNumSlots = 0;
		return(false);
		}
	m_MaxNumSlots = HashRange;
	}

// calculate the maximum number of hash buckets required
// could be less because masked (containing N) sequence hashes are not saved
MaxNumHashes = (SeqLen + WindowSize - 1) / WindowSize;

if(m_MaxNumHashes < MaxNumHashes || m_pHashes == NULL)
	{
	if(m_pHashes != NULL)
		delete m_pHashes;
	m_pHashes = new tsHashEntry[MaxNumHashes];
	if(m_pHashes == NULL)
		{
		delete m_pSlots;
		m_pSlots = NULL;
		m_MaxNumSlots = 0;
		return(false);
		}
	m_MaxNumHashes = MaxNumHashes;
	}

m_CurHashHeader.HashSlots = HashRange;
m_CurHashHeader.SeqLen	  = SeqLen;
m_CurHashHeader.WindowSize = WindowSize;

memset(m_pSlots,0,HashRange * sizeof(tsHashArraySlot));
MaxEntries = 0;
SeqOfs = 0;
NxtHashID = 1;
pCurHash = &m_pHashes[0];

while(SeqLen >= WindowSize)
	{
	pTargWindow = &pSeq[SeqOfs];

	if(LocateBase(eBaseN,WindowSize,pTargWindow)< 0) // slough if contains any masked or indeterminate bases
		{
		if(!bChkComplexity || !IsLowComplexity(pTargWindow,WindowSize))
			{
			GenSenseHash(&SenseHash,pTargWindow,WindowSize);
			pSlot = &m_pSlots[SenseHash];
			if(pSlot->NumEntries)						// check if same sequence already hashed - not interested in repeats..
				{
				HashEntryID = pSlot->FirstHashID;
				ExactMatchCnt = 0;
				while(HashEntryID != 0)
					{
					pEntry = &m_pHashes[HashEntryID -1];
					pProbeWindow = &pSeq[pEntry->Psn];
					if(*pProbeWindow == *pTargWindow && 
						SeqSenseMatches(m_CurHashHeader.WindowSize,pProbeWindow,pTargWindow))
						ExactMatchCnt++;
					HashEntryID = pEntry->NextHashID;
					}
				if(ExactMatchCnt > 10)
					{
					SeqLen -= WindowSize;
					SeqOfs += WindowSize;
					continue;
					}
				}
		
			pSlot->NumEntries++;
			if(pSlot->NumEntries > MaxEntries)
				MaxEntries = pSlot->NumEntries;
			pCurHash->Psn = SeqOfs;
			pCurHash->NextHashID = pSlot->FirstHashID;
			pSlot->FirstHashID = NxtHashID++;
			pCurHash++;
			}
		}
	SeqLen -= WindowSize;
	SeqOfs += WindowSize;
	}

m_CurHashHeader.NumHashes = NxtHashID;

m_CurHashHeader.EntryID = CBioSeqFile::CreateEntry(pszSource,pszDescription,eBinDataType);
CBioSeqFile::AddData(sizeof(tsHashHeader),(unsigned char *)&m_CurHashHeader);
CBioSeqFile::AddData(m_CurHashHeader.HashSlots * sizeof(tsHashArraySlot),(unsigned char *)m_pSlots);
CBioSeqFile::AddData(m_CurHashHeader.NumHashes * sizeof(tsHashEntry),(unsigned char *)m_pHashes);
CBioSeqFile::AddData(m_CurHashHeader.SeqLen,pSeq);
CBioSeqFile::SealEntry();
AddErrMsg("CHashFile::AddEntry","Max number entries for any hash slot was: %d",MaxEntries);
return(true);
}


tsHashHeader *
CHashFile::LoadHashHeader(unsigned int EntryID)
{
unsigned int DataReq;
unsigned int DataLen;
if(m_CurHashHeader.EntryID == EntryID)
	return(&m_CurHashHeader);
DataReq = sizeof(tsHashHeader);
if((DataLen = GetDataLen(EntryID)) < DataReq)
	return(NULL);
if(DataReq != GetData(EntryID,eBinDataType,0,(unsigned char *)&m_CurHashHeader,DataReq))
	return(NULL);
return(&m_CurHashHeader);
}

tsHashArraySlot *
CHashFile::LoadHashSlots(unsigned int EntryID)
{
unsigned int DataOfs;
unsigned int DataReq;
unsigned int DataLen;

if(m_CurHashHeader.EntryID != EntryID)
	LoadHashHeader(EntryID);
if(m_SlotsEntryID == EntryID)
	return(m_pSlots);
m_SlotsEntryID = 0;

DataOfs = sizeof(tsHashHeader);
DataReq = m_CurHashHeader.HashSlots * sizeof(tsHashArraySlot);
if((DataLen = GetDataLen(EntryID)) < DataOfs + DataReq)
	return(NULL);

if(m_pSlots != NULL && m_MaxNumSlots < m_CurHashHeader.HashSlots)
	{
	delete m_pSlots;
	m_pSlots = NULL;
	m_MaxNumSlots = 0;
	}

if(m_pSlots == NULL)
	{
	m_pSlots = new tsHashArraySlot[m_CurHashHeader.HashSlots];
	if(m_pSlots == NULL)
		return(NULL);
	m_MaxNumSlots = m_CurHashHeader.HashSlots;
	}

if(GetData(EntryID,eBinDataType,DataOfs,(unsigned char *)m_pSlots,DataReq) != DataReq)
	return(NULL);
m_SlotsEntryID = EntryID;
return(m_pSlots);
}

tsHashEntry *
CHashFile::LoadHashes(unsigned int EntryID)
{
unsigned int DataOfs;
unsigned int DataReq;
unsigned int DataLen;

if(m_CurHashHeader.EntryID != EntryID)
	LoadHashHeader(EntryID);
if(m_HashesEntryID == EntryID)
	return(m_pHashes);
m_HashesEntryID = 0;

DataOfs = sizeof(tsHashHeader) + (m_CurHashHeader.HashSlots * sizeof(tsHashArraySlot));
DataReq = m_CurHashHeader.NumHashes * sizeof(tsHashEntry);
if((DataLen = GetDataLen(EntryID)) < DataOfs + DataReq)
	return(NULL);

if(m_pHashes != NULL && m_MaxNumHashes < m_CurHashHeader.NumHashes)
	{
	delete m_pHashes;
	m_pHashes = NULL;
	m_MaxNumHashes = 0;
	}

if(m_pHashes == NULL)
	{
	m_pHashes = new tsHashEntry[ m_CurHashHeader.NumHashes];
	if(m_pHashes == NULL)
		return(NULL);
	m_MaxNumHashes = m_CurHashHeader.NumHashes;
	}

if(GetData(EntryID,eBinDataType,DataOfs,(unsigned char *)m_pHashes,DataReq) != DataReq)
	return(NULL);
m_HashesEntryID = EntryID;
return(m_pHashes);
}


unsigned char *
CHashFile::LoadSeq(unsigned int EntryID)
{
unsigned int DataOfs;
unsigned int DataReq;
unsigned int DataLen;

if(m_CurHashHeader.EntryID != EntryID)
	LoadHashHeader(EntryID);
if(m_SeqEntryID == EntryID)
	return(m_pSeq);
m_SeqEntryID = 0;

DataOfs = sizeof(tsHashHeader) + (m_CurHashHeader.HashSlots * sizeof(tsHashArraySlot)) + (m_CurHashHeader.NumHashes * sizeof(tsHashEntry));
DataReq = m_CurHashHeader.SeqLen;
if((DataLen = GetDataLen(EntryID)) < DataOfs + DataReq)
	return(NULL);

if(m_pSeq != NULL && m_MaxSeqLen < m_CurHashHeader.SeqLen)
	{
	delete m_pSeq;
	m_pSeq = NULL;
	m_MaxSeqLen = 0;
	}

if(m_pSeq == NULL)
	{
	m_pSeq = new unsigned char[m_CurHashHeader.SeqLen];
	if(m_pSeq == NULL)
		return(NULL);
	m_MaxSeqLen = m_CurHashHeader.SeqLen;
	}

if(GetData(EntryID,eBinDataType,DataOfs,m_pSeq,DataReq) != DataReq)
	return(NULL);

m_SeqEntryID = EntryID;
return(m_pSeq);
}





bool 
CHashFile::LocateMatches(bool bChkComplexity,		 // true if low complexity sequences are to be filtered out
						unsigned int EntryID,
 					   char *pszTargDescr,		 // describes target
					   unsigned int ProbeID,	 // used to uniquely identify this probe when ProcessThisMatch is called
					   char *pszProbeDescr,		 // describes probe
					   unsigned int MinMatchLen, // matches must be of at least this length
					   etSeqBase *pProbe,
					   unsigned int ProbeLen,
					   int hXMLdumpFile)
{
tsHashArraySlot *pSlot;				// pts to array slot for hash
tsHashEntry *pEntry;				// pts to entry for hash 
unsigned int ProbeOfs;
tHashID SenseHash;
tHashID AntisenseHash;
unsigned int HashEntryID;
etSeqBase *pProbeWindow;
etSeqBase *pProbeWindowRight;
etSeqBase *pTargWindow;


char szBuff[20000];
int DumpSize;
int NumExtdMatches;

if(m_pExtdMatches == NULL)
	{
	m_pExtdMatches = new tsExtdMatch[cMaxExtdMatches];
	if(m_pExtdMatches == NULL)
		return(false);
	}
m_NumExtdMatches = 0;
m_Mode0Total=0;
m_Mode1Total=0;
m_Mode2Total=0;
m_Mode3Total=0;


if(LoadHashHeader(EntryID) == NULL)
	return(false);

if(LoadHashSlots(EntryID) == NULL)
	return(false);

if(LoadHashes(EntryID) == NULL)
	return(false);

if(LoadSeq(EntryID) == NULL)
	return(false);

if(ProbeLen < m_CurHashHeader.WindowSize)
	return(false);

AddErrMsg("CHashFile::LocateMatches"," (%d)%s vs (%d)%s processing..            ",EntryID,pszTargDescr,ProbeID,pszProbeDescr);
m_HashMask = m_CurHashHeader.HashSlots - 1;
ProbeOfs = 0;		// start probe window at left edge of probe
pProbeWindow = pProbe;
pProbeWindowRight = pProbe + m_CurHashHeader.WindowSize - 1;
GenSenseHash(&SenseHash,pProbeWindow,m_CurHashHeader.WindowSize);
GenAntisenseHash(&AntisenseHash,pProbeWindow,m_CurHashHeader.WindowSize);

while((ProbeLen - ProbeOfs) >= m_CurHashHeader.WindowSize) 
	{
	if(*pProbeWindow < eBaseN)
		{
		if(!bChkComplexity || !IsLowComplexity(pProbeWindow,m_CurHashHeader.WindowSize)) 
			{
			pSlot = &m_pSlots[SenseHash];
			if(pSlot->NumEntries != 0 && pSlot->NumEntries < 10000)
				{
				HashEntryID = pSlot->FirstHashID;
				while(HashEntryID != 0)
					{
					pEntry = &m_pHashes[HashEntryID -1];
					pTargWindow = &m_pSeq[pEntry->Psn];
					if(*pProbeWindow == *pTargWindow && 
						SeqSenseMatches(m_CurHashHeader.WindowSize,pProbeWindow,pTargWindow))
							ExtendMatch(EntryID,ProbeID,MinMatchLen,pProbe,ProbeLen,ProbeOfs,m_pSeq,m_CurHashHeader.SeqLen,pEntry->Psn,m_CurHashHeader.WindowSize,0);
					if(*pProbeWindowRight == *pTargWindow && SeqSenseRevMatches(m_CurHashHeader.WindowSize,pProbeWindow,pTargWindow))
							ExtendMatch(EntryID,ProbeID,MinMatchLen,pProbe,ProbeLen,ProbeOfs,m_pSeq,m_CurHashHeader.SeqLen,pEntry->Psn,m_CurHashHeader.WindowSize,1);
					HashEntryID = pEntry->NextHashID;
					}
				}

			pSlot = &m_pSlots[AntisenseHash];
			if(pSlot->NumEntries != 0 && pSlot->NumEntries < 20000)
				{
				HashEntryID = pSlot->FirstHashID;
				while(HashEntryID != 0)
					{
					pEntry = &m_pHashes[HashEntryID -1];
					pTargWindow = &m_pSeq[pEntry->Psn];
					if(*pProbeWindow == m_Sense2AntiMap[*pTargWindow] && SeqAntisenseMatches(m_CurHashHeader.WindowSize,pProbeWindow,pTargWindow))
							ExtendMatch(EntryID,ProbeID,MinMatchLen,pProbe,ProbeLen,ProbeOfs,m_pSeq,m_CurHashHeader.SeqLen,pEntry->Psn,m_CurHashHeader.WindowSize,2);

					if(*pProbeWindowRight == m_Sense2AntiMap[*pTargWindow] && SeqAntisenseRevMatches(m_CurHashHeader.WindowSize,pProbeWindow,pTargWindow))
							ExtendMatch(EntryID,ProbeID,MinMatchLen,pProbe,ProbeLen,ProbeOfs,m_pSeq,m_CurHashHeader.SeqLen,pEntry->Psn,m_CurHashHeader.WindowSize,3);
					HashEntryID = pEntry->NextHashID;
					}
				}
			}
		}
	ProbeOfs+=1;
	pProbeWindow+=1;
	pProbeWindowRight+=1;

	SenseHash <<= 1;
	AntisenseHash <<= 1;
	SenseHash |= pProbeWindow[m_CurHashHeader.WindowSize-2] & 0x001;
	AntisenseHash |= m_Sense2AntiMap[pProbeWindow[m_CurHashHeader.WindowSize-2]] & 0x001;
	SenseHash &= m_HashMask;
	AntisenseHash &= m_HashMask;

	
	}
// Order extended matches by ProbePsn, TargPsn, MatchLen, then Mode so that duplicates can be filtered out..
	AddErrMsg("CHashFile::LocateMatches"," %d vs %d processed, mode 0: %d mode 1: %d mode 2: %d mode 3: %d..",
					EntryID,ProbeID,m_Mode0Total,m_Mode1Total,m_Mode2Total,m_Mode3Total);

if(!m_NumExtdMatches)
	return(true);
AddErrMsg("CHashFile::LocateMatches"," sorting %d extended matches ready for deduping..",m_NumExtdMatches);
if(m_NumExtdMatches > 1)
	qsort(m_pExtdMatches,m_NumExtdMatches,sizeof(tsExtdMatch),CompareExtdMatches);
tsExtdMatch *pEl1 = m_pExtdMatches;
tsExtdMatch *pEl2 = pEl1;
NumExtdMatches = 0;
char szProbeMatched[16000];
char szTargetMatched[16000];

for(unsigned int Idx = 0; Idx < m_NumExtdMatches; Idx++,pEl2++)
	{
	if(Idx > 0 && (pEl1->ProbePsn == pEl2->ProbePsn &&
		pEl1->TargPsn == pEl2->TargPsn &&
		pEl1->MatchLen == pEl2->MatchLen))
		continue;
	NumExtdMatches++;
	CSeqTrans::MapSeq2Ascii(&m_pSeq[pEl2->TargPsn],pEl2->MatchLen,szTargetMatched);
	CSeqTrans::MapSeq2Ascii(&pProbe[pEl2->ProbePsn],pEl2->MatchLen,szProbeMatched);
	szProbeMatched[pEl2->MatchLen] = '\0';
	szTargetMatched[pEl2->MatchLen] = '\0';
	DumpSize = sprintf(szBuff,"\n<match>\n<mode>%d</mode>\n<tid>%d</tid>\n<tdescr>%s</tdescr>\n<pid>%d</pid>\n<pdescr>%s</pdescr>\n<tpsn>%d</tpsn>\n<ppsn>%d</ppsn>\n<len>%d</len>\n<tseq>%s</tseq>\n<pseq>%s</pseq>\n</match>",
			pEl2->Mode,EntryID,pszTargDescr,ProbeID,pszProbeDescr,pEl2->TargPsn,pEl2->ProbePsn,pEl2->MatchLen,szTargetMatched,szProbeMatched);
	if(write(hXMLdumpFile,szBuff,DumpSize)!=DumpSize)
		{
		AddErrMsg("CHashFile::LocateMatches"," write failed..%s",strerror(errno));
		return(false);
		}
	pEl1 = pEl2;
	}
#ifdef _WIN32
_commit(hXMLdumpFile);
#else
fsync(hXMLdumpFile);
#endif

AddErrMsg("CHashFile::LocateMatches"," total of %d matches written to XML..",NumExtdMatches);


return(true);
}

int CHashFile::CompareExtdMatches( const void *arg1, const void *arg2 )
{
tsExtdMatch *pEl1 = (tsExtdMatch *)arg1;
tsExtdMatch *pEl2 = (tsExtdMatch *)arg2;
if(pEl1->ProbePsn == pEl2->ProbePsn)
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



// SeqSenseMatchesLeft attempts to extend the match left
unsigned int
CHashFile::SeqSenseMatchesLeft(unsigned int ProbeLen,
							   etSeqBase *pProbe,
							   unsigned int ProbePsn,	// current probe match left start
   							   unsigned int TargLen,
							   etSeqBase *pTarg,
							   unsigned int TargPsn)	// current target match left start 
{
unsigned int LeftDelta = 0;
etSeqBase *pPLeft = &pProbe[ProbePsn-1];
etSeqBase *pTLeft = &pTarg[TargPsn-1];

while(LeftDelta < ProbePsn && LeftDelta < TargPsn)
	{
	if(*pPLeft == eBaseN || *pTLeft == eBaseN)
		return(LeftDelta);

	if(*pPLeft-- != *pTLeft--)
		return(LeftDelta);
	LeftDelta++;
	}
return(LeftDelta);
}

// SeqSenseMatchesRight attempts to extend the match left
unsigned int
CHashFile::SeqSenseMatchesRight(unsigned int ProbeLen,
							   etSeqBase *pProbe,
							   unsigned int ProbePsn,	// current probe match right
							   unsigned int TargLen,
							   etSeqBase *pTarg,
							   unsigned int TargPsn)	// current target match right 
{
unsigned int RightDelta = 0;
etSeqBase *pPRight = &pProbe[ProbePsn];
etSeqBase *pTRight = &pTarg[TargPsn];

ProbeLen -= ProbePsn + 1;
TargLen -= TargPsn + 1;

while(RightDelta < ProbeLen && RightDelta < TargLen)
	{
	if(*pPRight == eBaseN || *pTRight == eBaseN)
		return(RightDelta);

	if(*pPRight++ != *pTRight++)
		return(RightDelta);
	RightDelta++;
	}
return(RightDelta);
}


// SeqAntisenseMatchesLeft attempts to extend the match left
unsigned int
CHashFile::SeqAntisenseMatchesLeft(unsigned int ProbeLen,
							   etSeqBase *pProbe,
							   unsigned int ProbePsn,	// current probe match left start
   							   unsigned int TargLen,
							   etSeqBase *pTarg,
							   unsigned int TargPsn)	// current target match left start 
{
unsigned int LeftDelta = 0;
etSeqBase *pPLeft = &pProbe[ProbePsn-1];
etSeqBase *pTLeft = &pTarg[TargPsn-1];

while(LeftDelta < ProbePsn && LeftDelta < TargPsn)
	{
	if(*pPLeft == eBaseN || *pTLeft == eBaseN)
		return(LeftDelta);


	if(m_Sense2AntiMap[*pPLeft--] != *pTLeft--)
		return(LeftDelta);
	LeftDelta++;
	}
return(LeftDelta);
}

// SeqAntisenseMatchesRight attempts to extend the match right
unsigned int
CHashFile::SeqAntisenseMatchesRight(unsigned int ProbeLen,
							   etSeqBase *pProbe,
							   unsigned int ProbePsn,	// current probe match right start
							   unsigned int TargLen,
							   etSeqBase *pTarg,
							   unsigned int TargPsn)	// current target match right start 
{
unsigned int RightDelta = 0;
etSeqBase *pPRight = &pProbe[ProbePsn];
etSeqBase *pTRight = &pTarg[TargPsn];

ProbeLen -= ProbePsn + 1;
TargLen -=  TargPsn + 1;

while(RightDelta < ProbeLen && RightDelta < TargLen)
	{
	if(*pPRight == eBaseN || *pTRight == eBaseN)
		return(RightDelta);

	if(m_Sense2AntiMap[*pPRight++] != *pTRight++)
		return(RightDelta);
	RightDelta++;
	}
return(RightDelta);
}


// SeqSenseRevMatchesLeft attempts to extend the match left
unsigned int
CHashFile::SeqSenseRevMatchesLeft(unsigned int ProbeLen,
							   etSeqBase *pProbe,
							   unsigned int ProbePsn,	// current probe match left
   							   unsigned int TargLen,
							   etSeqBase *pTarg,
							   unsigned int TargPsn)	// current target match right 
{
unsigned int LeftDelta = 0;
etSeqBase *pPLeft = &pProbe[ProbePsn-1];
etSeqBase *pTRight = &pTarg[TargPsn];
TargLen -= TargPsn + 1;

while(LeftDelta < ProbePsn && LeftDelta < TargLen)
	{
	if(*pTRight == eBaseN || *pPLeft == eBaseN)
		return(LeftDelta);


	if(*pTRight++ != *pPLeft--)
		return(LeftDelta);
	LeftDelta++;
	}
return(LeftDelta);
}

// SeqSenseRevMatchesRight attempts to extend the match right
unsigned int
CHashFile::SeqSenseRevMatchesRight(unsigned int ProbeLen,
							   etSeqBase *pProbe,
							   unsigned int ProbePsn,	// current probe match right start
							   unsigned int TargLen,
							   etSeqBase *pTarg,
							   unsigned int TargPsn)	// current target match right start 
{
unsigned int RightDelta = 0;
etSeqBase *pPRight = &pProbe[ProbePsn];
etSeqBase *pTLeft =  &pTarg[TargPsn-1];

ProbeLen -= ProbePsn + 1;

while(RightDelta < ProbeLen && RightDelta < TargPsn)
	{
	if(*pPRight == eBaseN || *pTLeft == eBaseN)
		return(RightDelta);

	if(*pPRight++ != *pTLeft--)
		return(RightDelta);
	RightDelta++;
	}
return(RightDelta);
}


// SeqAntisenseRevMatchesLeft attempts to extend the match left
unsigned int
CHashFile::SeqAntisenseRevMatchesLeft(unsigned int ProbeLen,
							  etSeqBase *pProbe,
							   unsigned int ProbePsn,	// current probe match left
   							   unsigned int TargLen,
							   etSeqBase *pTarg,
							   unsigned int TargPsn)	// current target match left 
{
unsigned int LeftDelta = 0;
etSeqBase *pPLeft = &pProbe[ProbePsn-1];
etSeqBase *pTRight = &pTarg[TargPsn];
TargLen -= TargPsn + 1;

while(LeftDelta < ProbePsn && LeftDelta < TargLen)
	{
	if(*pPLeft == eBaseN || *pTRight == eBaseN)
		return(LeftDelta);

	if(m_Sense2AntiMap[*pPLeft--] != *pTRight++)
		return(LeftDelta);
	LeftDelta++;
	}
return(LeftDelta);
}

// SeqAntisenseMatchesRight attempts to extend the match right
unsigned int
CHashFile::SeqAntisenseRevMatchesRight(unsigned int ProbeLen,
							   etSeqBase *pProbe,
							   unsigned int ProbePsn,	// current probe match right start
							   unsigned int TargLen,
							   etSeqBase *pTarg,
							   unsigned int TargPsn)	// current target match right start 
{
unsigned int RightDelta = 0;
etSeqBase *pPRight = &pProbe[ProbePsn];
etSeqBase *pTLeft =  &pTarg[TargPsn-1];

ProbeLen -= ProbePsn + 1;

while(RightDelta < ProbeLen && RightDelta < TargPsn)
	{
	if(*pPRight == eBaseN || *pTLeft == eBaseN)
		return(RightDelta);

	if(m_Sense2AntiMap[*pPRight++] != *pTLeft--)
		return(RightDelta);
	RightDelta++;
	}
return(RightDelta);
}



// ProcessMatch
// An match of length MatchLen has been determined to exist between &pProbe[ProbeMatchPsn] and &pTarget[TargMatchPsn]
// Processing attempts to maximise the match length by extending left and right
// Returns true if extended matchlength >= ThresLen
bool 
CHashFile::ExtendMatch(unsigned int EntryID,
					   unsigned int ProbeID,
					   unsigned int ThresLen,		// threshold extended length
					    etSeqBase *pProbe,
						unsigned int ProbeLen,
				    	 unsigned int ProbeMatchPsn,
						etSeqBase *pTarget,
						 unsigned int TargLen,
						 unsigned int TargMatchPsn,
					     unsigned int MatchLen,
						 int MatchMode)
{
unsigned int LeftPsnDelta;
unsigned int RightPsnDelta;
tsExtdMatch *pExtdMatch;
//char szSeqBuffer[4096];

if(EntryID == ProbeID && ProbeMatchPsn == TargMatchPsn) // not interested in selfmatches!!!
	return(false);

switch(MatchMode) {
	case 0:				// matches forward on probe strand against forward on target sense strand (p:acgtcg matches t:acgtcg)
		m_Mode0Total++;
		LeftPsnDelta = SeqSenseMatchesLeft(ProbeLen,pProbe,ProbeMatchPsn,TargLen,pTarget,TargMatchPsn);
		RightPsnDelta = SeqSenseMatchesRight(ProbeLen,pProbe,ProbeMatchPsn+MatchLen,TargLen,pTarget,TargMatchPsn+MatchLen);
		ProbeMatchPsn -= LeftPsnDelta;
		TargMatchPsn -= LeftPsnDelta;
		MatchLen += LeftPsnDelta + RightPsnDelta;
		break;
		
	case 1:				// matches forward on probe strand against reverse on target sense strand (p:acgtcg matches t:gctgca)
		m_Mode1Total++;
		LeftPsnDelta  = SeqSenseRevMatchesLeft(ProbeLen,pProbe,ProbeMatchPsn,TargLen,pTarget,TargMatchPsn+MatchLen);
		RightPsnDelta = SeqSenseRevMatchesRight(ProbeLen,pProbe,ProbeMatchPsn+MatchLen,TargLen,pTarget,TargMatchPsn);
		ProbeMatchPsn -= LeftPsnDelta;
		TargMatchPsn -= RightPsnDelta;
		MatchLen += LeftPsnDelta + RightPsnDelta;
		break;

	case 2:				// matches forward on probe strand against forward on target antisense strand (p:acgtcg matches t:tgcagc)
		m_Mode2Total++;
		LeftPsnDelta = SeqAntisenseMatchesLeft(ProbeLen,pProbe,ProbeMatchPsn,TargLen,pTarget,TargMatchPsn);
		RightPsnDelta = SeqAntisenseMatchesRight(ProbeLen,pProbe,ProbeMatchPsn+MatchLen,TargLen,pTarget,TargMatchPsn+MatchLen);
		ProbeMatchPsn -= LeftPsnDelta;
		TargMatchPsn -= LeftPsnDelta;
		MatchLen += LeftPsnDelta + RightPsnDelta;
		break;

	case 3:				// matches forward on probe strand against reverse on target antisense strand (p:acgtcg matches t:cgacgt)
		m_Mode3Total++;
		LeftPsnDelta = SeqAntisenseRevMatchesLeft(ProbeLen,pProbe,ProbeMatchPsn,TargLen,pTarget,TargMatchPsn+MatchLen);
		RightPsnDelta = SeqAntisenseRevMatchesRight(ProbeLen,pProbe,ProbeMatchPsn+MatchLen,TargLen,pTarget,TargMatchPsn);
		ProbeMatchPsn -= LeftPsnDelta;
		TargMatchPsn -= RightPsnDelta;
		MatchLen += LeftPsnDelta + RightPsnDelta;
		break;
	
	default:
		return(false);
	}

if(MatchLen < ThresLen)
	return(false);

if(MatchLen > 3000)				// clamp matchlen to a reasonable limit, otherwise will start overwriting buffers!
	MatchLen = 3000;			// extra long frozen sequences are assembly artifacts or repeats..

if(m_NumExtdMatches == cMaxExtdMatches)
	return(false);
pExtdMatch = &m_pExtdMatches[m_NumExtdMatches++];
pExtdMatch->ProbePsn = ProbeMatchPsn;
pExtdMatch->TargPsn = TargMatchPsn;
pExtdMatch->MatchLen = MatchLen;
pExtdMatch->Mode = MatchMode;
return(true);
}

bool
CHashFile::IsLowComplexity(etSeqBase *pSeq,unsigned int Len)
{
unsigned int BaseCnts[4];				// used to hold counts of individual nucleotides
unsigned int DimerCnts[16];				// used to hold counts of nucleotide pairs
unsigned int DistinctBases;
unsigned int DistinctDimers;
unsigned int Idx;
unsigned int Cnt;
unsigned int Ofs;

if(pSeq == NULL || Len < 10)
	return(true);

memset(BaseCnts,0,sizeof(BaseCnts));
memset(DimerCnts,0,sizeof(DimerCnts));
DistinctBases = DistinctDimers = 0;
for(Idx = 0, Cnt = 0; Cnt < Len; Cnt++)
	{
	Ofs = *pSeq++;
	if(Ofs >= 4)					    // any eBaseN's? treat as though it's a low complexity sequence
		return(true);
	Idx |= Ofs;							// update counts of dimers
	if(!BaseCnts[Idx & 0x03])
		DistinctBases++;
	BaseCnts[Idx & 0x03]++;				// update counts of individual nucleotides
	if(Cnt)								// after 1st then update dimer counts
		{
		if(!DimerCnts[Idx & 0x0f])
			DistinctDimers++;
		DimerCnts[Idx & 0x0f]++;
		}
	Idx <<= 2;
	}

// now analyse the counts..
// if only 2 bases are counted then treat as low complexity
if(DistinctBases <= 2)
	return(true);
// if 5 or less distinct dimers then treat as low complexity
if(DistinctDimers <= 5)
	return(true);

return(false);
}
