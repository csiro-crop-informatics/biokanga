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

#include "SSRIdentify.h"


CSSRIdentify::CSSRIdentify(void)
{
m_pSeqBuff = NULL;
m_pKMerDist = NULL;
m_pszRptSSRsBuff = NULL;
Init();
}

CSSRIdentify::~CSSRIdentify(void)
{
if(m_pSeqBuff != NULL)
	{
#ifdef _WIN32
	free(m_pSeqBuff);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqBuff != MAP_FAILED)
		munmap(m_pSeqBuff,m_AllocdSeqBuffMem);
#endif
	m_pSeqBuff = NULL;
	}
if(m_pKMerDist != NULL)
	{
#ifdef _WIN32
	free(m_pKMerDist);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pKMerDist != MAP_FAILED)
		munmap(m_pKMerDist,m_AllocdKMerFreqMem);
#endif
	m_pKMerDist = NULL;
	}
if(m_pszRptSSRsBuff != NULL)
	delete m_pszRptSSRsBuff;
}


void
CSSRIdentify::Init(void)
{
m_pSeqBuff = NULL;
m_pKMerDist = NULL;
m_pszRptSSRsBuff = NULL;
m_hOutFile = -1;
m_hOutKMerFreqFile = -1;
Reset();
m_CurTime.Start();
}

void
CSSRIdentify::Reset(void)
{
if(m_hOutFile != -1)
	{
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if(m_hOutKMerFreqFile != -1)
	{
#ifdef _WIN32
	_commit(m_hOutKMerFreqFile);
#else
	fsync(m_hOutKMerFreqFile);
#endif
	close(m_hOutKMerFreqFile);
	m_hOutKMerFreqFile = -1;
	}

if(m_pszRptSSRsBuff != NULL)
	{
	delete m_pszRptSSRsBuff;
	m_pszRptSSRsBuff = NULL;
	}

if(m_pSeqBuff != NULL)
	{
#ifdef _WIN32
	free(m_pSeqBuff);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeqBuff != MAP_FAILED)
		munmap(m_pSeqBuff,m_AllocdSeqBuffMem);
#endif
	m_pSeqBuff = NULL;
	}
if(m_pKMerDist != NULL)
	{
#ifdef _WIN32
	free(m_pKMerDist);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pKMerDist != MAP_FAILED)
		munmap(m_pKMerDist,m_AllocdKMerFreqMem);
#endif
	m_pKMerDist = NULL;
	}

m_AllocdSeqBuffMem = 0;
m_SeqBuffLen = 0;
m_AllocdKMerFreqMem = 0;
m_KMerFreqLen = 0;
m_IdxRptSSRs = 0;
m_CurTime.Stop();
}


etSeqBase *
CSSRIdentify::AllocSeqBuff(size_t SeqLen)				// allocate for at least this sequence length
{
size_t memreq;
etSeqBase *pTmp;

if(m_pSeqBuff != NULL && m_AllocdSeqBuffMem >= SeqLen)
	return(m_pSeqBuff);

if(m_pSeqBuff == NULL)
	{
	memreq = max(SeqLen,(size_t)cMaxAllocBuffChunk);
#ifdef _WIN32
	m_pSeqBuff = (etSeqBase *) malloc(SeqLen);	// initial and perhaps the only allocation
	if(m_pSeqBuff == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqBuff: Memory allocation of %lld bytes - %s",(INT64)SeqLen,strerror(errno));
		return(NULL);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pSeqBuff = (etSeqBase *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSeqBuff == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqBuff: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
		m_pSeqBuff = NULL;
		return(NULL);
		}
#endif
	}
else
	{
	memreq = SeqLen + cMaxAllocBuffChunk;
#ifdef _WIN32
	pTmp = (etSeqBase *) realloc(m_pSeqBuff,memreq);
#else
	pTmp = (etSeqBase *)mremap(m_pSeqBuff,m_AllocdSeqBuffMem,memreq,MREMAP_MAYMOVE);
	if(pTmp == MAP_FAILED)
		pTmp = NULL;
#endif
	if(pTmp == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AllocSeqBuff: Memory re-allocation to %lld bytes - %s",(INT64)memreq,strerror(errno));
		return(NULL);
		}
	m_pSeqBuff = pTmp;
	}
m_AllocdSeqBuffMem = memreq;
return(m_pSeqBuff);
}

// ProcessBioseqFile
// Process input biosequence file
int
CSSRIdentify::ProcessBioseqFile(char *pszFile)
{
CBioSeqFile BioSeqFile;
char szSource[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
UINT32 SeqLen;
int Rslt;
tBSFEntryID CurEntryID;

if((Rslt=BioSeqFile.Open(pszFile,cBSFTypeSeq,false))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d %s",Rslt,BioSeqFile.GetErrMsg());
	return(Rslt);
	}

CurEntryID = 0;
while((Rslt = CurEntryID = BioSeqFile.Next(CurEntryID)) > eBSFSuccess)
	{
	BioSeqFile.GetNameDescription(CurEntryID,cBSFSourceSize-1,(char *)&szSource,
											cBSFDescriptionSize-1,(char *)&szDescription);
	SeqLen = BioSeqFile.GetDataLen(CurEntryID);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing %s|%s",szSource,szDescription);

	if(!SeqLen)
		continue;

	if(AllocSeqBuff(SeqLen) == NULL)
		{
		Rslt = eBSFerrMem;
		break;
		}

	if((Rslt = BioSeqFile.GetData(CurEntryID,eSeqBaseType,0,m_pSeqBuff,SeqLen)) != SeqLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d %s",Rslt,BioSeqFile.GetErrMsg());
		break;
		}

	if((Rslt=IdentifySSRs(szSource,pszFile,SeqLen,m_pSeqBuff)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessBioseqFile - error %d",Rslt);
		break;
		}

	ReportProgress();
	}
if(Rslt == eBSFerrEntry)
	Rslt = eBSFSuccess;

BioSeqFile.Close();
return(Rslt);
}

// ProcessFastaFile
// Parse input fasta format file
int
CSSRIdentify::ProcessFastaFile(char *pszFile)
{
CFasta Fasta;

size_t AvailBuffSize;
char szName[cBSFSourceSize];
char szDescription[cBSFDescriptionSize];
UINT32 SeqLen;
int Descrlen;
bool bFirstEntry;
bool bEntryCreated;
int Rslt;
int SeqID;

if((Rslt=Fasta.Open(pszFile,true))!=eBSFSuccess)
	{
	if(Rslt != eBSFerrNotFasta)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile: Unable to open '%s' [%s] %s",pszFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	return(Rslt);
	}

if(m_pSeqBuff == NULL)				// if not already allocated then allocate to hold cMaxAllocBuffChunk bases 
	{
	SeqLen = cMaxAllocBuffChunk;
	if(AllocSeqBuff(SeqLen) == NULL)
		{
		Rslt = eBSFerrMem;
		Fasta.Close();
		return(Rslt);
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessFastaFile:- Adding %s..",pszFile);

bFirstEntry = true;
bEntryCreated = false;
SeqID = 0;
m_SeqBuffLen = 0;
AvailBuffSize = m_AllocdSeqBuffMem;
while((Rslt = SeqLen = Fasta.ReadSequence(&m_pSeqBuff[m_SeqBuffLen],(int)min(AvailBuffSize,(size_t)cMaxAllocBuffChunk),true,false)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		SeqID++;
		if(bEntryCreated)				// add any previous entry
			{
			if((Rslt=IdentifySSRs(szName,pszFile,m_SeqBuffLen,m_pSeqBuff)) < eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessFastaFile - error %d",Rslt);
				break;
				}
			}
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		// An assumption - will one day bite real hard - is that the
		// fasta descriptor line starts with some form of unique identifier.
		// Use this identifier as the entry name.
		if(sscanf(szDescription," %s[ ,]",szName)!=1)
			sprintf(szName,"%s.%d",pszFile,++SeqID);

		bFirstEntry = false;
		bEntryCreated = true;
		m_SeqBuffLen = 0;
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

	m_SeqBuffLen += SeqLen;
	AvailBuffSize -= SeqLen;
	if(AvailBuffSize < (size_t)(cMaxAllocBuffChunk / 8))
		{
		if(AllocSeqBuff(m_AllocdSeqBuffMem + SeqLen) == NULL)
			{
			Rslt = eBSFerrMem;
			Fasta.Close();
			return(Rslt);
			}
		AvailBuffSize = m_AllocdSeqBuffMem - m_SeqBuffLen;
		}
	}

if(Rslt >= eBSFSuccess && bEntryCreated && m_SeqBuffLen > 0)			// close entry
	if((Rslt=IdentifySSRs(szName,pszFile,m_SeqBuffLen,m_pSeqBuff)) < eBSFSuccess)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"IdentifySSRs - error %d",Rslt);
	else
		Rslt = eBSFSuccess;

return(Rslt);
}


int
CSSRIdentify::ReportCSV(int RepElLen,	// identified SSR contains repeat elements of this length
	   int NumTandemEls,			// each repeat element is repeated this many times
	   INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
	    char *pszDescr,				// descriptor for the targeted sequence
		char *pszInFile,			// sequence parsed from this file
	   INT64 TargSeqLen,			// targeted sequence is this length
	   etSeqBase *pTargSeq)			// targeted sequence within which the SSR has been located
{
static bool bFirst = true;

if(bFirst)
	{
	m_IdxRptSSRs += sprintf(m_pszRptSSRsBuff,"\"ID\",\"Proc\",\"Species\",\"Chrom\",\"Start\",\"End\",\"Len\",\"Strand\",\"KMer\",\"Rpts\",\"Seq\"\n");
	bFirst = false;
	}
if(m_IdxRptSSRs > (cMaxAllocRptSSRs - 2000))
	{
	CUtility::SafeWrite(m_hOutFile,m_pszRptSSRsBuff,m_IdxRptSSRs);
	m_IdxRptSSRs = 0;
	}

m_IdxRptSSRs += sprintf(&m_pszRptSSRsBuff[m_IdxRptSSRs],"%d,\"SSRs\",\"N/A\",\"%s\",%lld,%lld,%d,\"+\",%d,%d,\"",m_TotNumAcceptedSSRs,
						pszDescr,SSRStartOfs,SSRStartOfs+(RepElLen*NumTandemEls)-1,RepElLen*NumTandemEls,RepElLen,NumTandemEls);
CSeqTrans::MapSeq2UCAscii(&pTargSeq[SSRStartOfs],RepElLen * NumTandemEls,&m_pszRptSSRsBuff[m_IdxRptSSRs]);
m_IdxRptSSRs += RepElLen  * NumTandemEls;
m_IdxRptSSRs += sprintf(&m_pszRptSSRsBuff[m_IdxRptSSRs],"\"\n");
return(m_IdxRptSSRs);
}

int
CSSRIdentify::ReportBED(int RepElLen,	// identified SSR contains repeat elements of this length
	   int NumTandemEls,			// each repeat element is repeated this many times
	   INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
	    char *pszDescr,				// descriptor for the targeted sequence
		char *pszInFile,			// sequence parsed from this file
	   INT64 TargSeqLen,			// targeted sequence is this length
	   etSeqBase *pTargSeq)			// targeted sequence within which the SSR has been located
{
return(0);
}

int
CSSRIdentify::ReportSAM(int RepElLen,	// identified SSR contains repeat elements of this length
	   int NumTandemEls,			// each repeat element is repeated this many times
	   INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
	    char *pszDescr,				// descriptor for the targeted sequence
		char *pszInFile,			// sequence parsed from this file
	   INT64 TargSeqLen,			// targeted sequence is this length
	   etSeqBase *pTargSeq)			// targeted sequence within which the SSR has been located
{
return(0);
}

// Reporting as CSV, BED or SAM
int
CSSRIdentify::Report(int RepElLen,	// identified SSR contains repeat elements of this length
	   int NumTandemEls,			// each repeat element is repeated this many times
	   INT64 SSRStartOfs,			// repeat element starts at this offset within the targeted sequence
	    char *pszDescr,				// descriptor for the targeted sequence
		char *pszInFile,			// sequence parsed from this file
	   INT64 TargSeqLen,			// targeted sequence is this length
	   etSeqBase *pTargSeq)			// targeted sequence within which the SSR has been located
{
int Rslt;
switch(m_RptSSRsFormat) {
	case eRFCsv:			// CSV
		Rslt = ReportCSV(RepElLen,NumTandemEls,SSRStartOfs,pszDescr,pszInFile,TargSeqLen,pTargSeq);
		break;

	case eRFBed:			// BED
		Rslt = ReportBED(RepElLen,NumTandemEls,SSRStartOfs,pszDescr,pszInFile,TargSeqLen,pTargSeq);
		break;

	case eRFSam:			// SAM
		Rslt = ReportSAM(RepElLen,NumTandemEls,SSRStartOfs,pszDescr,pszInFile,TargSeqLen,pTargSeq);
		break;
		break;
	}
return(Rslt);
}


UINT16
GenSeqHash16(int SeqLen,	// hash this length sequence
			etSeqBase *pSeq) // sequence to hash
{
UINT32 Hash;
int Idx;

Hash = 19937;			// circular prime as hash seed

// hash over the element
for(Idx = 0; Idx < SeqLen; Idx++)
	{
	Hash = (Hash ^ (*pSeq & 0x07)) * 3119;	// a circular prime
	Hash ^= (Hash >> 13);
	Hash &= 0x0ffff;
	}
if(Hash == 0)			// 0 reserved as an error indicator so force hash to be 1
	Hash += 1;
return(Hash);
}


// if kmer is less than 10bp then count number of instances of KMer

int
CSSRIdentify::CntKMer(int KMerLen,	 // count this KMer
					int Rpts,		 // tandem repeat counts
					etSeqBase *pSeq) // sequence 
{
int Idx;
UINT32 FreqIdx;
tsKMerDist *pKMerDist;

if(m_pKMerDist == NULL || KMerLen != m_KMerFreqLen)
	return(0);

FreqIdx = 0;
for(Idx = 0; Idx < KMerLen; Idx++)
	{
	FreqIdx <<= 2;
	FreqIdx |= 0x03 & *pSeq++;
	}

pKMerDist = (tsKMerDist *)((UINT8 *)m_pKMerDist + (m_SizeOfKMerDist * FreqIdx));

pKMerDist->Cnt += 1;
pKMerDist->TandemRpts[Rpts - m_MinTandemRpts] += 1;
return(pKMerDist->Cnt);
}


int
CSSRIdentify::IdentifySSRs(char *pszDescr,	// descriptor for the targeted sequence
			 char *pszInFile,			// sequence parsed from this file
			 INT64 TargSeqLen,			// sequence length of targeted sequence within which to search for SSRs
			 etSeqBase *pTargSeq)		// identify SSRs in this targeted sequence
{
bool bSlough;
etSeqBase BaseA;
etSeqBase BaseB;
INT64 Ofs;
INT64 SSRStartOfs;
int RepElLen;
int RepElOfs;
int NumTandemEls;
int NumSSRs;
etSeqBase *pBase;
etSeqBase *pRepElBase;

NumSSRs = 0;
for(RepElLen = m_MinRepElLen; RepElLen <= m_MaxRepElLen; RepElLen++)
	{
	pBase = pTargSeq;
	SSRStartOfs = 0;
	NumTandemEls = 0;
	bSlough = false;
	for(Ofs = 0; Ofs < (TargSeqLen - RepElLen); Ofs++,pBase++)
		{
		pRepElBase = pBase;
		for(RepElOfs = 0; RepElOfs < RepElLen; RepElOfs++,pRepElBase++)
			{
			BaseA = *pRepElBase & 0x07;
			BaseB = pRepElBase[RepElLen] & 0x07;
			if(BaseA > eBaseT || BaseB > eBaseT)	// not interested in potential SSR if it contains any non-cannonical base
				{
				bSlough = true;
				continue;
				}
			if(BaseA != BaseB)
				break;
			}
		if(RepElOfs == RepElLen)
			{
			if(NumTandemEls == 0)
				{
				SSRStartOfs = Ofs;
				NumTandemEls = 1;
				}
			NumTandemEls += 1;
			Ofs += RepElLen - 1;
			pBase = &pTargSeq[Ofs];
			}
		else
			{
			if(!bSlough && NumTandemEls >= m_MinTandemRpts)
				{
				if(NumTandemEls > m_MaxTandemRpts)
					m_TotNumExcessiveTandemSSRs += 1;
				else
					{
					m_TotNumAcceptedSSRs += 1;
					Report(RepElLen,NumTandemEls,SSRStartOfs,pszDescr,pszInFile,TargSeqLen,pTargSeq);
					NumSSRs += 1;
					m_TotNumAcceptedKmerSSRs[RepElLen] += 1;
					CntKMer(RepElLen,NumTandemEls,&pTargSeq[SSRStartOfs]);
					}
				}
			NumTandemEls = 0;
			bSlough = false;
			}
		if(!(Ofs % 100))
			ReportProgress();
		}
	}
return(NumSSRs);
}	

int
CSSRIdentify::ReportKMers(char *pszKMerFreqFile)	// report SSR repeating element K-mer frequencies to this file
{
int Idx;
tsKMerDist *pFreqDist;
UINT32 SeqIdx;
char Base;
int Idy;
int Rpts;
char szKmer[100];
int BuffIdx;
char szBuff[0x03fff];

if(m_hOutKMerFreqFile != -1 && m_KMerFreqLen != 0)
	{
	BuffIdx = 0;
	pFreqDist = m_pKMerDist;
	for(Idx = 0; Idx < m_NumKMers; Idx+=1)
		{
		SeqIdx = Idx;
		for(Idy = 0; Idy < m_KMerFreqLen; Idy++)
			{
			switch(SeqIdx & 0x03) {
				case 0:
					Base = 'A';
					break;
				case 1:
					Base = 'C';
					break;
				case 2:
					Base = 'G';
					break;
				case 3:
					Base = 'T';
					break;
				}
			szKmer[m_KMerFreqLen - (1+Idy)] = Base;
			SeqIdx >>= 2;
			}
		szKmer[m_KMerFreqLen] = '\0';
		BuffIdx += sprintf(&szBuff[BuffIdx],"\"%s\",%d",szKmer,pFreqDist->Cnt);
		for(Rpts = cMinTandemRpts; Rpts <= m_MaxTandemRpts; Rpts++)
			if(Rpts < m_MinTandemRpts)
				BuffIdx += sprintf(&szBuff[BuffIdx],",0");
			else
				BuffIdx += sprintf(&szBuff[BuffIdx],",%d",pFreqDist->TandemRpts[Rpts - m_MinTandemRpts]);
		BuffIdx += sprintf(&szBuff[BuffIdx],"\n");
		if((BuffIdx + 1000) > (int)sizeof(szBuff))
			{
			CUtility::SafeWrite(m_hOutKMerFreqFile,szBuff,BuffIdx);
			BuffIdx = 0;
			}
		pFreqDist = (tsKMerDist *)((UINT8 *)pFreqDist + m_SizeOfKMerDist);
		}
	if(BuffIdx > 0)
		CUtility::SafeWrite(m_hOutKMerFreqFile,szBuff,BuffIdx);
#ifdef _WIN32
	_commit(m_hOutKMerFreqFile);
#else
	fsync(m_hOutKMerFreqFile);
#endif
	close(m_hOutKMerFreqFile);
	m_hOutKMerFreqFile = -1;
	}
return(eBSFSuccess);
}


int
CSSRIdentify::ReportProgress(bool bForce)	// let user know that there is processing activity, normally progress reported every 30 sec unless bForce set true
{
static long PrevSecs = 0;
long CurSecs = m_CurTime.ReadUSecs();
if(bForce || CurSecs < PrevSecs || ((CurSecs - PrevSecs) > 60))
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Progress: Identified %u SSRs",m_TotNumAcceptedSSRs);
	PrevSecs = CurSecs;
	}
return(0);
}

int
CSSRIdentify::Process(int PMode,			// procesisng mode - currently unused..
		teRptSSRsFromat RptSSRsFormat,		// report SSRs in this file format
		int MinRepElLen,					// identify repeating elements of this minimum length
		int MaxRepElLen,					// ranging upto this maximum length
		int MinTandemRpts,					// minimum number of tandem repeats
		int MaxTandemRpts,					// maximum number of repeats
		int NumInFileSpecs,					// number of input, could be wildcarded, file specs
		char *pszInFiles[],					// files to be processed
		char *pszKMerFreqFile,				// optional, output element KMer freq to this file
		char *pszOutFile)					// SSRs to this file
{
int Rslt;
int Idx;
UINT32 TotNumAcceptedSSRs;
UINT32 TotNumExcessiveTandemSSRs;

char *pszInFile;
CSimpleGlob glob(SG_GLOB_FULLSORT);

m_RptSSRsFormat = RptSSRsFormat;
m_MinRepElLen = MinRepElLen;
m_MaxRepElLen = MaxRepElLen;
m_MinTandemRpts = MinTandemRpts;
m_MaxTandemRpts = MaxTandemRpts;

m_TotNumExcessiveTandemSSRs = 0;
m_TotNumAcceptedSSRs = 0;
TotNumAcceptedSSRs = 0;
TotNumExcessiveTandemSSRs = 0;

memset(m_TotNumAcceptedKmerSSRs,0,sizeof(m_TotNumAcceptedKmerSSRs));

// if single KMer element length being processed and that KMer is <= 10 then allocate to hold instance counts
if(pszKMerFreqFile != NULL && pszKMerFreqFile[0] != '\0' &&
   MinRepElLen == MaxRepElLen && MinRepElLen <= 10)
	{
	size_t memreq;
	int Power2 = MinRepElLen;
	m_SizeOfKMerDist = (int)(sizeof(tsKMerDist) + (sizeof(UINT32) * (MaxTandemRpts - MinTandemRpts)));
	m_NumKMers = 1;
	while(Power2--)
		m_NumKMers <<= 2;
	memreq = (size_t)m_NumKMers * (size_t)m_SizeOfKMerDist;
#ifdef _WIN32
	m_pKMerDist = (tsKMerDist *) malloc(memreq);	// initial and only allocation
	if(m_pKMerDist == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pKMerDist = (tsKMerDist *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pKMerDist == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
		m_pKMerDist = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocdKMerFreqMem = memreq;
	memset(m_pKMerDist,0,memreq);
	m_KMerFreqLen = MaxRepElLen;
#ifdef _WIN32
	if((m_hOutKMerFreqFile = open(pszKMerFreqFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hOutKMerFreqFile = open(pszKMerFreqFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",pszKMerFreqFile,strerror(errno));
		Reset();
		return(eBSFerrCreateFile);
		}
	char szBuff[2000];
	int BuffIdx;
	BuffIdx = sprintf(szBuff,"\"K-Mer\",\"Cnt\"");
	for(Idx = cMinTandemRpts; Idx <= m_MaxTandemRpts; Idx++)
		BuffIdx += sprintf(&szBuff[BuffIdx],",\"Rpts:%d\"",Idx);
	BuffIdx += sprintf(&szBuff[BuffIdx],"\n");
	CUtility::SafeWrite(m_hOutKMerFreqFile,szBuff,BuffIdx);
	}
else
	{
	m_AllocdKMerFreqMem = 0;
	m_SizeOfKMerDist = 0;
	m_KMerFreqLen = 0;
	m_NumKMers = 0;
	}

if((m_pszRptSSRsBuff = new char [cMaxAllocRptSSRs]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d bytes buffering for SSRs output",cMaxAllocRptSSRs);
	Reset();
	return(eBSFerrMem);
	}
m_IdxRptSSRs = 0;

#ifdef _WIN32
if((m_hOutFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutFile = open(pszOutFile,O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Output results file created/truncated: '%s'",pszOutFile);

for(Idx = 0; Idx < NumInFileSpecs; Idx++)
	{
	glob.Init();
	if(glob.Add(pszInFiles[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"Unable to glob '%s",pszInFiles[Idx]);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	if(glob.FileCount() <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any input file matching '%s",pszInFiles[Idx]);
		continue;
		}
	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
		{
		pszInFile = glob.File(FileID);
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing file '%s'",pszInFile);
		if((Rslt = ProcessFastaFile(pszInFile)) < eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLWarn,gszProcName,"Failed procesisng file '%s",pszInFile);
			return(Rslt);
			}
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted %u SSRs, rejected %u excessively long",
							m_TotNumAcceptedSSRs - TotNumAcceptedSSRs,
							m_TotNumExcessiveTandemSSRs - TotNumExcessiveTandemSSRs);
		TotNumAcceptedSSRs = m_TotNumAcceptedSSRs;
		TotNumExcessiveTandemSSRs = m_TotNumExcessiveTandemSSRs;
		}
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total SSRs accepted %u, rejected %u tandem repeats as excessively long",m_TotNumAcceptedSSRs,m_TotNumExcessiveTandemSSRs);
for(Idx=MinRepElLen; Idx <= MaxRepElLen; Idx++)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"K-mer %d, Number %u",Idx,m_TotNumAcceptedKmerSSRs[Idx]);



if(m_hOutFile != -1)
	{
	if(m_IdxRptSSRs)
		{
		CUtility::SafeWrite(m_hOutFile,m_pszRptSSRsBuff,m_IdxRptSSRs);
		m_IdxRptSSRs = 0;
		}
#ifdef _WIN32
	_commit(m_hOutFile);
#else
	fsync(m_hOutFile);
#endif
	close(m_hOutFile);
	m_hOutFile = -1;
	}


Rslt = ReportKMers(pszKMerFreqFile);
Reset();
return(Rslt);
}