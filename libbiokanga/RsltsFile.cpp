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

CRsltsFile::CRsltsFile(void)
{
m_pFirstAllocdBlock = NULL;		// pts to first allocated matches block 
m_hRsltsFile = -1;
Reset();
}

CRsltsFile::~CRsltsFile(void)
{
Reset();
}

int 
CRsltsFile::Open(char *pszRsltsFile, // specifies file to truncate/create
				 bool bRsltsXML,		// false: CSV, true: XML
				 bool bClusterMatches) // true if matches are to be clustered (matches overlap on either probe or target)
{
char szBuff[1000];
int Len;
Reset();
#ifdef _WIN32
m_hRsltsFile = open(pszRsltsFile,O_CREATETRUNC);
#else
if((m_hRsltsFile = open(pszRsltsFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
      if(ftruncate(m_hRsltsFile,0)!=0)
			{
			AddErrMsg("CRsltsFile::Open","Unable to truncate %s - %s",pszRsltsFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif
if(m_hRsltsFile < 0)
	{
	AddErrMsg("CRsltsFile::Open","Unable to open %s - %s",pszRsltsFile,strerror(errno));
	return(eBSFerrOpnFile);
	}
strcpy(m_szRsltsFile,pszRsltsFile);
m_bRsltsXML = bRsltsXML;
m_bClusterMatches = bClusterMatches;
if(bRsltsXML)
	{
	Len = sprintf(szBuff,"<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n");
	Len += sprintf(&szBuff[Len],"<dataroot xmlns:od=\"urn:schemas-microsoft-com:officedata\" ");
	Len += sprintf(&szBuff[Len],"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" ");
	Len += sprintf(&szBuff[Len],"xsi:noNamespaceSchemaLocation=\"match.xsd\" generated=\"2004-04-25T15:30:02\" >");
	}
else
	{
	Len = sprintf(szBuff,"\"id\",\"mode\",\"tid\",\"tdescr\",\"pid\",\"pdescr\",");
	Len+= sprintf(&szBuff[Len],"\"plen\",\"tpsn\",\"ppsn\",\"len\",\"ident\",\"numclustered\",\"ispartner\",\"mskd\",\"cplx\",\"pseq\"");
	}
if(write(m_hRsltsFile,szBuff,Len)!=Len)
	{
	AddErrMsg("CRsltsFile::Open","Write to %s - %s",pszRsltsFile,strerror(errno));
	close(m_hRsltsFile);
	m_hRsltsFile = -1;
	return(eBSFerrFileAccess);
	}
return(eBSFSuccess);
}

int 
CRsltsFile::Reset(void)			// reset state back to that immediately following instantiation
{
if(m_hRsltsFile != -1)
	{
	close(m_hRsltsFile);
	m_hRsltsFile = -1;
	}

if(m_pFirstAllocdBlock != NULL)
	{
	do {
		m_pLastAllocdBlock = m_pFirstAllocdBlock->pNext;
		delete m_pFirstAllocdBlock;
		}
	while((m_pFirstAllocdBlock = m_pLastAllocdBlock)!= NULL);
	m_pFirstAllocdBlock = NULL;
	m_pLastAllocdBlock = NULL;
	}

m_pMRAMatch = NULL;
m_pFreeMatches = NULL;
m_NumExtdMatches = 0;
m_UniqueResultID = 0;
return(eBSFSuccess);
}

int 
CRsltsFile::Close(void)				// closes opened results file
{
int Len;
char szBuff[100];
if(m_hRsltsFile == -1)
	return(eBSFerrClosed);
if(m_bRsltsXML)
	{
	Len = sprintf(szBuff,"\n</dataroot>\n");
	if(write(m_hRsltsFile,szBuff,Len)!=Len)
		{
		AddErrMsg("CRsltsFile::Close","Write failed - %s",strerror(errno));
		close(m_hRsltsFile);
		m_hRsltsFile = -1;
		Reset();
		return(eBSFerrFileAccess);
		}
	}
close(m_hRsltsFile);
m_hRsltsFile = -1;
Reset();
return(eBSFSuccess);
}

int
CRsltsFile::StartRsltSet(unsigned int TargEntryID, char *pszTarget,	// descriptive text about the target
						 unsigned int ProbeEntryID, char *pszProbe)	// descriptive text about the probe
{
sSfxHeaderBlock *pBlock;
tsSfxMatch *pMatch;

m_NumExtdMatches = 0;
m_ProbeEntryID = ProbeEntryID;
m_TargEntryID = TargEntryID;

strncpy(m_szTargDescr,pszTarget,sizeof(m_szTargDescr));
m_szTargDescr[sizeof(m_szTargDescr)-1] = '\0';
strncpy(m_szProbeDescr,pszProbe,sizeof(m_szProbeDescr));
m_szProbeDescr[sizeof(m_szProbeDescr)-1] = '\0';

CBioSeqFile::MakeXMLsafe(m_szTargDescr);
CBioSeqFile::MakeXMLsafe(m_szProbeDescr);

// ensure that there is one  block allocation avail to hold extended matches
if(m_pFirstAllocdBlock == NULL)
	{
	m_pFirstAllocdBlock = new sSfxHeaderBlock;		
	if(m_pFirstAllocdBlock == NULL)
		return(eBSFerrMem);
	m_pFirstAllocdBlock->pNext = NULL;
	m_pFirstAllocdBlock->pPrev = NULL;
	}

	// if more than one block then delete all but the first to free up memory
if((pBlock = m_pFirstAllocdBlock->pNext) != NULL)
	{
	m_pFirstAllocdBlock->pNext = NULL;
	do {
		m_pLastAllocdBlock = pBlock->pNext;
		delete pBlock;
		}
	while((pBlock = m_pLastAllocdBlock)!= NULL);
	}
m_pLastAllocdBlock = m_pFirstAllocdBlock;

	// link tsSfxMatch'es ready for subsequent allocations from m_pFreeMatches
pMatch = &m_pFirstAllocdBlock->Matches[0];
m_pFreeMatches = pMatch;
m_pMRAMatch = NULL;
for(int Idx = 0; Idx < (cBlockExtdMatches - 1); Idx++,pMatch+=1)
	{
	pMatch->pPrev = NULL;
	pMatch->pNext = pMatch+1;
	}
pMatch->pNext = NULL;
return(eBSFerrClosed);
}


int
CRsltsFile::EndRsltSet(unsigned int ProbeLen,
					   unsigned char *pProbeSeq)	// NULL if no sequence to be output
{
char szBuff[cMaxMatchSeqSeq + 5001];
char szTargSeq[cMaxMatchSeqSeq+1];
int Len;
int TextSeqLen;
int NumMasked;
int PctMasked;
int CplxScore;
tsSfxMatch *pEl1;

// iterate over each result set
if(pEl1 = m_pMRAMatch)
	{
	szTargSeq[0] = '\0';
	do {
		TextSeqLen = pEl1->MatchLen > cMaxMatchSeqSeq ? cMaxMatchSeqSeq : pEl1->MatchLen;
		PctMasked = 0;
		CplxScore = 0;
		szTargSeq[0] = '\0';
		if(pProbeSeq != NULL)
			{
			CSeqTrans::MapSeq2Ascii(&pProbeSeq[pEl1->ProbePsn],TextSeqLen,szTargSeq);
			NumMasked = CBioSeqFile::GetNumMasked(TextSeqLen,&pProbeSeq[pEl1->ProbePsn]);
			if(NumMasked > 0)
				PctMasked = (NumMasked * 100) / TextSeqLen;
			CplxScore = CBioSeqFile::ScoreComplexity(&pProbeSeq[pEl1->ProbePsn],TextSeqLen);
			}
		szTargSeq[TextSeqLen] = '\0';
	
		if(m_bRsltsXML)
			{
			Len = sprintf(szBuff,"\n<match>\n<id>%d</id>\n<mode>%d</mode>\n<tid>%d</tid>\n<tdescr>%s</tdescr>\n<pid>%d</pid>\n<pdescr>%s</pdescr>\n<plen>%d</plen>\n<tpsn>%d</tpsn>\n<ppsn>%d</ppsn>\n<len>%d</len>\n<ident>%d</ident>\n<numclustered>%d</numclustered>\n<ispartner>%s</ispartner>\n<mskd>%d</mskd>\n<cplx>%d</cplx>\n<pseq>%s</pseq>\n</match>",
				++m_UniqueResultID,
				pEl1->Mode,
				m_TargEntryID,
				m_szTargDescr,
				m_ProbeEntryID,
				m_szProbeDescr,
				ProbeLen,
				pEl1->TargPsn,
				pEl1->ProbePsn,
				pEl1->MatchLen,
				pEl1->IdentCnt,
				pEl1->NumInCluster,
				pEl1->IsPartner ? "yes" : "no",
				PctMasked,
				CplxScore,
				szTargSeq);
			}
		else
			{
			Len = sprintf(szBuff,"\n%d,%d,%d,\"%s\",%d,\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,\"%s\"",
				++m_UniqueResultID,
				pEl1->Mode,
				m_TargEntryID,
				m_szTargDescr,
				m_ProbeEntryID,
				m_szProbeDescr,
				ProbeLen,
				pEl1->TargPsn,
				pEl1->ProbePsn,
				pEl1->MatchLen,
				pEl1->IdentCnt,
				pEl1->NumInCluster,
				pEl1->IsPartner ? 1 : 0,
				PctMasked,
				CplxScore,
				szTargSeq);
			}
		if(write(m_hRsltsFile,szBuff,Len)!=Len)
			{
			AddErrMsg("CRsltsFile::EndRsltSet","Write to %s - %s",m_szRsltsFile,strerror(errno));
			close(m_hRsltsFile);
			m_hRsltsFile = -1;
			return(eBSFerrFileAccess);
			}
		}
	while(pEl1 = pEl1->pPrev);
	}
#ifdef _WIN32
_commit(m_hRsltsFile);
#else
fsync(m_hRsltsFile);
#endif
return(eBSFSuccess);
}

int
CRsltsFile::AddMatch(unsigned int ProbePsn,			// psn (0..n) in probe at which match starts
					 unsigned int TargPsn,			// psn (0..n) in target at which match starts
					 unsigned int MatchLen,			// match length
					 unsigned int IdentCnt,			// how many bases in match were exact matches
					 unsigned int Mode)
{
tsSfxMatch *pExtdMatch;
sSfxHeaderBlock *pAllocdBlock;
if(MatchLen > cMaxMatchSeqSeq)
	MatchLen = cMaxMatchSeqSeq;

// need more tsSfxMatch nodes to hold this match?
if(m_pFreeMatches == NULL)
	{
	pAllocdBlock = new sSfxHeaderBlock;		
	if(pAllocdBlock == NULL)
		return(eBSFerrMem);

	pAllocdBlock->pNext = NULL;
	pAllocdBlock->pPrev = m_pLastAllocdBlock;
	m_pLastAllocdBlock->pNext = pAllocdBlock;
	m_pLastAllocdBlock = pAllocdBlock;

	pExtdMatch = &m_pLastAllocdBlock->Matches[0];
	m_pFreeMatches = pExtdMatch;
	for(int Idx = 0; Idx < (cBlockExtdMatches-1); Idx++, pExtdMatch += 1)
		{
		pExtdMatch->pPrev = NULL;
		pExtdMatch->pNext = pExtdMatch+1;
		}
	pExtdMatch->pNext = NULL;
	}

pExtdMatch = m_pFreeMatches;
m_pFreeMatches = m_pFreeMatches->pNext;
pExtdMatch->pNext = NULL;
pExtdMatch->pPrev = m_pMRAMatch;
if(m_pMRAMatch != NULL)
	m_pMRAMatch->pNext = pExtdMatch;

m_pMRAMatch = pExtdMatch;
m_pMRAMatch->ProbePsn = ProbePsn;
m_pMRAMatch->TargPsn = TargPsn;
m_pMRAMatch->MatchLen = MatchLen;
m_pMRAMatch->IdentCnt = IdentCnt;
m_pMRAMatch->IsPartner = false;
m_pMRAMatch->NumInCluster = 0;
m_pMRAMatch->Mode = Mode;
return(eBSFSuccess);
}



// check if putative match is contained within another already matched sequence
// if putative match is contained within another then it is treated as if already matched
// if putative match contains another match then it replaces the existing match
bool 
CRsltsFile::IsAlreadyMatched(unsigned int ProbePsn,
								  unsigned int TargPsn,
								  unsigned int MatchLen,
								  unsigned int Mode)
{
bool bContained;
tsSfxMatch *pExtdMatch;
unsigned int ProbeEnd;
unsigned int ExtdProbeEnd;
unsigned int TargEnd;
unsigned int ExtdTargEnd;

if((pExtdMatch = m_pMRAMatch)==NULL)
	return(false);
if(MatchLen > cMaxMatchSeqSeq)
	MatchLen = cMaxMatchSeqSeq;
bContained = false;
ProbeEnd = ProbePsn + MatchLen - 1;
TargEnd = TargPsn + MatchLen - 1;
while(pExtdMatch != NULL)
	{
	if(Mode == pExtdMatch->Mode)
	  {
	   ExtdProbeEnd = pExtdMatch->ProbePsn + pExtdMatch->MatchLen - 1;
	   ExtdTargEnd = pExtdMatch->TargPsn + pExtdMatch->MatchLen - 1;

	   	  // if match completely contained within an existing then slough..
	   if((ProbePsn >= pExtdMatch->ProbePsn &&
			ProbeEnd <= ExtdProbeEnd) &&
			(TargPsn >= pExtdMatch->TargPsn && 
			TargEnd <= ExtdTargEnd))
			return(true);

	      // if match completely contains existing then simply replace existing
	      // if first, or if subsequent then remove...
	   if((ProbePsn <= pExtdMatch->ProbePsn &&
			ProbeEnd >= ExtdProbeEnd) &&
			(TargPsn <= pExtdMatch->TargPsn && 
			TargEnd >= ExtdTargEnd))
			{
			if(bContained)
				{
				if(pExtdMatch->pNext != NULL)
					pExtdMatch->pNext->pPrev = pExtdMatch->pPrev;
				if(pExtdMatch->pPrev != NULL)
					pExtdMatch->pPrev->pNext = pExtdMatch->pNext;
	
				pExtdMatch->pNext = m_pFreeMatches;
				m_pFreeMatches = pExtdMatch;
				}
			else
				{
				bContained = true;
				pExtdMatch->MatchLen = MatchLen;
				pExtdMatch->ProbePsn = ProbePsn;
				pExtdMatch->TargPsn = TargPsn;
				}
			}
		}
	pExtdMatch = pExtdMatch->pPrev;
	}
return(bContained);
}


// cluster
// cluster matches by overlap on either ppsn or tpsn
void 
CRsltsFile::ClusterMatches(void)
{
bool bRemoved = false;
tsSfxMatch *pExtdMatch;
tsSfxMatch *pScanMatch;

unsigned int Mode;
unsigned int ProbePsn;
unsigned int ProbeEnd;
unsigned int TargPsn; 
unsigned int MatchLen;
unsigned int NumInCluster;

if(!m_bClusterMatches)	// made optional as this can be very resource intensive
	return;

if((pExtdMatch = m_pMRAMatch)==NULL)
	return;

while(pExtdMatch != NULL)
	{
	Mode = pExtdMatch->Mode;
	ProbePsn = pExtdMatch->ProbePsn;
	MatchLen = pExtdMatch->MatchLen;
	ProbeEnd = ProbePsn + MatchLen - 1;

	TargPsn = pExtdMatch->TargPsn;
	pScanMatch = m_pMRAMatch;
	NumInCluster = 0;
	while(pScanMatch != NULL)
		{
		if(((pScanMatch->ProbePsn >= ProbePsn && pScanMatch->ProbePsn <= ProbeEnd) ||
			((pScanMatch->ProbePsn + MatchLen) >= ProbePsn && (pScanMatch->ProbePsn + MatchLen) <= ProbeEnd)))
			NumInCluster++;
		pScanMatch = pScanMatch->pPrev;
		}
	pExtdMatch->NumInCluster = NumInCluster;
	pExtdMatch = pExtdMatch->pPrev;
	}
	
}


tsSfxMatch *
CRsltsFile::LocateFirstPPsnMatch(unsigned int psnstart,unsigned int psnend)
{
tsSfxMatch *pExtdMatch;
if((pExtdMatch = m_pMRAMatch)==NULL)
	return(NULL);
while(pExtdMatch != NULL)
	if(pExtdMatch->ProbePsn >= psnstart &&
		pExtdMatch->ProbePsn >= psnend)
		return(pExtdMatch);
return(NULL);
}

tsSfxMatch *
CRsltsFile::LocateNextPPsnMatch(tsSfxMatch *pExtdMatch,unsigned int psnstart,unsigned int psnend)
{
if(pExtdMatch == NULL)
	return(NULL);
while((pExtdMatch = pExtdMatch->pPrev) != NULL)
	{
	if(pExtdMatch->ProbePsn >= psnstart &&
		pExtdMatch->ProbePsn >= psnend)
		return(pExtdMatch);
	}
return(NULL);
}


tsSfxMatch *
CRsltsFile::LocateFirstTPsnMatch(unsigned int psnstart,unsigned int psnend)
{
tsSfxMatch *pExtdMatch;
if((pExtdMatch = m_pMRAMatch)==NULL)
	return(NULL);
while(pExtdMatch != NULL)
	if(pExtdMatch->TargPsn >= psnstart &&
		pExtdMatch->TargPsn >= psnend)
		return(pExtdMatch);
return(NULL);
}

tsSfxMatch *
CRsltsFile::LocateNextTPsnMatch(tsSfxMatch *pExtdMatch,unsigned int psnstart,unsigned int psnend)
{
if(pExtdMatch == NULL)
	return(NULL);
while((pExtdMatch = pExtdMatch->pPrev) != NULL)
	{
	if(pExtdMatch->TargPsn >= psnstart &&
		pExtdMatch->TargPsn >= psnend)
		return(pExtdMatch);
	}
return(NULL);
}



// when looking for self matches (matches within the same sequence) then
// there will be duplicate matches in which a 2nd match will be made from a probe position back to an earlier match
// that had been made to the current probe postion.
void 
CRsltsFile::MarkDuplicateMatches(void)
{
bool bRemoved = false;
tsSfxMatch *pExtdMatch;
tsSfxMatch *pScanMatch;

unsigned int Mode;
unsigned int ProbePsn; 
unsigned int TargPsn; 
unsigned int MatchLen;

if((pExtdMatch = m_pMRAMatch)==NULL)
	return;

while(pExtdMatch != NULL)
	{
	if(!pExtdMatch->IsPartner)
		{
		Mode = pExtdMatch->Mode;
		ProbePsn = pExtdMatch->ProbePsn;
		TargPsn = pExtdMatch->TargPsn;
		MatchLen = pExtdMatch->MatchLen;
		pScanMatch = m_pMRAMatch;
		while(pScanMatch != NULL)
			{
			if(!pScanMatch->IsPartner && 
				pScanMatch != pExtdMatch &&
				Mode == pScanMatch->Mode &&
				ProbePsn == pScanMatch->TargPsn &&
				TargPsn == pScanMatch->ProbePsn &&
				MatchLen == pScanMatch->MatchLen)
				pScanMatch->IsPartner = true;
			pScanMatch = pScanMatch->pPrev;
			}
		}
	pExtdMatch = pExtdMatch->pPrev;
	}
return;
}





// NearlyExtend
// Checks and extends an existing match which can be extended left or right
// Returns true if match was extended or if match already contained within a prev match
bool CRsltsFile::NearlyExtend(unsigned int ProbePsn,
                              unsigned int TargPsn,
							  unsigned int MatchLen,
							  unsigned int IdentCnt,		// how many bases in match were exact matches
							  unsigned int Mode)
{
tsSfxMatch *pExtdMatch;
unsigned int ProbeEnd;
unsigned int ExtdProbeEnd;
unsigned int TargEnd;
unsigned int ExtdTargEnd;

if((pExtdMatch = m_pMRAMatch)==NULL)
	return(false);
if(MatchLen > cMaxMatchSeqSeq)
	MatchLen = cMaxMatchSeqSeq;
while(pExtdMatch != NULL) 
	{	
	if(Mode == pExtdMatch->Mode)
	   {
	   ProbeEnd = ProbePsn + MatchLen - 1;
	   ExtdProbeEnd = pExtdMatch->ProbePsn + pExtdMatch->MatchLen - 1;
	   TargEnd = TargPsn + MatchLen - 1;
	   ExtdTargEnd = pExtdMatch->TargPsn + pExtdMatch->MatchLen - 1;

	   	  // if match completely contained within existing then slough..
	   if((ProbePsn >= pExtdMatch->ProbePsn &&
			ProbeEnd <= ExtdProbeEnd) &&
			(TargPsn >= pExtdMatch->TargPsn && 
			TargEnd <= ExtdTargEnd))
			return(true);

	   if((ProbePsn == pExtdMatch->ProbePsn &&
		  TargPsn == pExtdMatch->TargPsn) ||
			(ProbeEnd == ExtdProbeEnd &&
			TargEnd == ExtdTargEnd))
			{
			if(MatchLen <= pExtdMatch->MatchLen)
				return(true);

			pExtdMatch->MatchLen = MatchLen;
			pExtdMatch->ProbePsn = ProbePsn;
			pExtdMatch->TargPsn = TargPsn;
			pExtdMatch->IdentCnt = IdentCnt;
			return(true);
			}
		}
	pExtdMatch = pExtdMatch->pPrev;
	}
return(false);
}



