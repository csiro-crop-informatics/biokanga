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

#include "Transcriptome.h"


CTranscriptome::CTranscriptome(void)
{
m_pCSV = NULL;
m_pTransRegions = NULL;
m_pTranscripts = NULL;
m_hFragsFile = -1;
m_hTransFile = -1;
m_hDistFile = -1;
Reset(false);
}


CTranscriptome::~CTranscriptome(void)
{
Reset(false);
}

int
CTranscriptome::Reset(bool bSync)
{
if(m_hFragsFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hFragsFile);
#else
		fsync(m_hFragsFile);
#endif
	close(m_hFragsFile);
	m_hFragsFile = -1;
	}

if(m_hDistFile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hDistFile);
#else
		fsync(m_hDistFile);
#endif
	close(m_hDistFile);
	m_hDistFile = -1;
	}

if(m_hTransFile != -1)
	{
	close(m_hTransFile);
	m_hTransFile = -1;
	}

if(m_pTransRegions != NULL)
	{
#ifdef _WIN32
	free(m_pTransRegions);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pTransRegions != MAP_FAILED)
		munmap(m_pTransRegions,m_MemAllocdTransRegions);
#endif
	m_pTransRegions = NULL;
	}
if(m_pTranscripts != NULL)
	{
#ifdef _WIN32
	free(m_pTranscripts);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pTranscripts != MAP_FAILED)
		munmap(m_pTranscripts,m_MemAllocdTrnscripts);
#endif
	m_pTranscripts = NULL;
	}

if(m_pCSV != NULL)
	{
	delete m_pCSV;
	m_pCSV = NULL;
	}

m_NumRegionsParsed = 0;
m_NumTransRegions = 0;			// number of transcribed regions
m_NumAllocdTransRegions = 0;	// currently allocated to hold this many transcribed regions
m_MemAllocdTransRegions = 0;	// memory size currently allocated to hold m_NumAllocdTransRegions transcribed regions

m_NumTranscripts = 0;			// current total number of transcripts
m_CurTotTransLen = 0;			// current total concatenated transcript length
m_MemAllocdTrnscripts = 0;		// memory size currently allocated to hold 

m_szTransFile[0] = '\0';
m_szFragsFile[0] = '\0';
m_szDistFile[0] = '\0';
m_TargFragLen = 0;
m_NumFragLenBins = 0;
m_MinFragLen = 0;
m_MaxFragLen = 0;
m_StericOcc5Len = 0;
m_StericOcc3Len = 0;
return(eBSFSuccess);
}


TRandomCombined<TRanrotWGenerator,TRandomMersenne> RGseeds((int)time(0));

int
CTranscriptome::Fragmentate(char *pszTransFile,	// transcripts and their abundance are from this CSV file
                     char *pszFragsFile,	// fragmentation counts to this file
					 char *pszDistFile,		// output fragment length distributions (CSV) file
					 int TargFragLen,		// fragment until mean fragmentation length is <= this	
					 int NumFragLenBins,	// bin fragment lengths into this many bins so as to obtain a distribution
					 int MinFragLen,		// count fragments with >= this length and
					 int MaxFragLen,		// <= this length
			 		 int StericOcc5Len,		// 5' steric occlusion region length
					 int StericOcc3Len)		// 3' steric occlusion region length
{
int Rslt;

int LenDists[1001];		// length distibutions from 0 to 9999 and over in 10bp increments
int NumFields;
int ElID;
char *pszElName;
int ElLen;
int ElInsts;
double RandSiteVal;
size_t RandSiteOfs;
size_t TotalNumFrags;
UINT8 *pBase;
UINT8 Flgs;
int MeanFragLen;
int TotFragsInRange;
int OccIdx;

int NumUnderlen;
int NumUnderInsts;
int NumOverlen;
int NumOverInsts;
int NumRegionOverlen;
int	CurFragsUnderLen;
int	CurFragsOverLen;

char szRsltsBuff[8196];
int RsltBuffIdx;
int CurFragLen;
int CurFragsInRange;
size_t ExpTotSeqLen;
size_t CurTotSeqLen;
tsTransRegion *pTransRegion;
UINT8 *pSeq;

memset(LenDists,0,sizeof(LenDists));

Reset(false);
strncpy(m_szTransFile,pszTransFile,sizeof(m_szTransFile)-1);
m_szTransFile[sizeof(m_szTransFile)-1] = '\0';
strncpy(m_szFragsFile,pszFragsFile,sizeof(m_szFragsFile)-1);
m_szTransFile[sizeof(m_szFragsFile)-1] = '\0';
strncpy(m_szDistFile,pszDistFile,sizeof(m_szFragsFile)-1);
m_szDistFile[sizeof(m_szDistFile)-1] = '\0';

m_TargFragLen = TargFragLen;
m_NumFragLenBins = NumFragLenBins;
m_MinFragLen = MinFragLen;
m_MaxFragLen = MaxFragLen;
m_StericOcc5Len = StericOcc5Len;
m_StericOcc3Len = StericOcc3Len;

// open and process transcripts from CSV file pszTransFile
if((m_pCSV = new CCSVFile) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrMem);
	}

if((Rslt=m_pCSV->Open(m_szTransFile))!=eBSFSuccess)
	{
	while(m_pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",m_szTransFile);
	Reset(false);
	return(Rslt);
	}

#ifdef _WIN32
m_hFragsFile = open(pszFragsFile,O_CREATETRUNC );
#else
if((m_hFragsFile = open(pszFragsFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
if(ftruncate(m_hFragsFile,0)!=0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszFragsFile,strerror(errno));
	Reset(false);
	return(eBSFerrCreateFile);
	}
#endif

if(m_hFragsFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszFragsFile);
	Reset(false);
	return(eBSFerrCreateFile);
	}
RsltBuffIdx = sprintf(szRsltsBuff,"\"ID\",\"Feat\",\"Len\",\"Insts\",\"NumFragsInLen\",\"NumFragsUnderlen\",\"NumFragsOverLen\"\n");
if(write(m_hFragsFile,szRsltsBuff,RsltBuffIdx)!=RsltBuffIdx)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fragmentate: Write %d bytes to file '%s' failed - %s",RsltBuffIdx,pszFragsFile,strerror(errno));
		Reset(false);
		return(eBSFerrFileAccess);
		}
RsltBuffIdx = 0;


#ifdef _WIN32
m_hDistFile = open(pszDistFile,O_CREATETRUNC );
#else
if((m_hDistFile = open(pszDistFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
if(ftruncate(m_hDistFile,0)!=0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszDistFile,strerror(errno));
	Reset(false);
	return(eBSFerrCreateFile);
	}
#endif

if(m_hDistFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszDistFile);
	Reset(false);
	return(eBSFerrCreateFile);
	}
RsltBuffIdx = sprintf(szRsltsBuff,"\"MinLen\",\"MaxLen\",\"Insts\"\n");
if(write(m_hDistFile,szRsltsBuff,RsltBuffIdx)!=RsltBuffIdx)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fragmentate: Write %d bytes to file '%s' failed - %s",RsltBuffIdx,pszDistFile,strerror(errno));
		Reset(false);
		return(eBSFerrFileAccess);
		}
RsltBuffIdx = 0;

m_NumRegionsParsed = 0;
NumUnderlen = 0;
NumUnderInsts = 0;
NumOverlen = 0;
NumOverInsts = 0;
NumRegionOverlen = 0;
while((Rslt=m_pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pCSV->GetCurFields();
	if(!m_NumRegionsParsed && (NumFields >= cExpNumInputFields && m_pCSV->IsLikelyHeaderLine())) // expecting at least cExpNumInputFields fields so only then check for header line
		continue;

	if(NumFields < cExpNumInputFields)	// must have as a minimum the  <FeatID>,<FeatName>,<FeatLen>,<FeatInsts>
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least %d fields in '%s', GetCurFields() returned '%d'",cExpNumInputFields,m_szTransFile,NumFields);
		Reset(false);
		return(eBSFerrFieldCnt);
		}


	m_pCSV->GetInt(1,&ElID);
	m_pCSV->GetText(2,&pszElName);
	m_pCSV->GetInt(3,&ElLen);
	m_pCSV->GetInt(4,&ElInsts);
	if(ElLen < 20)
		{
		if(NumUnderlen++ <= 10)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sloughed a underlength ( < 20 ) transcribed region '%d,%s,%d,%d' in '%s'",ElID,pszElName, ElLen,ElInsts,m_szTransFile);
		continue;
		}

	if(ElLen > 200000)
		{
		if(NumOverlen++ <= 10)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sloughed overlength ( > 200Kb ) region '%d,%s,%d,%d' in '%s'",ElID,pszElName, ElLen,ElInsts,m_szTransFile);
		continue;
		}

	if(ElInsts < 1)
		{
		if(NumUnderInsts++ <= 10)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sloughed a 0 occurances transcribed region '%d,%s,%d,%d' in '%s'",ElID,pszElName, ElLen,ElInsts,m_szTransFile);
		continue;
		}

	if(ElInsts > 1000000)
		{
		if(NumOverInsts++ <= 10)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sloughed 1M+ occurances transcribed region '%d,%s,%d,%d' in '%s'",ElID,pszElName, ElLen,ElInsts,m_szTransFile);
		continue;
		}

	// have to ensure that no one transcribed region totally dominates all memory resources
	if(((size_t)ElLen * ElInsts) > (size_t)4000000000)
		{
		if(NumRegionOverlen++ <= 10)
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sloughed a region for which the sum sequence length is > 4Gbp '%d,%s,%d,%d' in '%s'",ElID,pszElName, ElLen,ElInsts,m_szTransFile);
		continue;
		}

	m_NumRegionsParsed += 1;
	if((Rslt = AddTranscripts(m_NumRegionsParsed,ElLen,ElInsts))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors adding transcribed region '%s'",pszElName);
		Reset(false);
		return(Rslt);
		}
	}
m_pCSV->Close();

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total sloughed %d under length, and %d over length",NumUnderlen,NumOverlen);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total sloughed %d 0 instances, and %d >= 1M instances",NumUnderInsts,NumOverInsts);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Total sloughed %d in which the total region sequence length sums to >= 4Gbp",NumRegionOverlen);


if(!m_NumRegionsParsed)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Nothing to do, no transcribed rgions parsed and accepted from file '%s'",m_szTransFile);
	Reset();
	return(eBSFSuccess);
	}

// randomly select a fragmentation site and if not sterically occluded then fragment at that site
// when fragmenting then need to set occlusion up/down stream regions
// when fragmenting need to update mean fragment lengths

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Fragmentation started...");


TotalNumFrags = m_NumTranscripts;
MeanFragLen = (int)(m_CurTotTransLen / TotalNumFrags);
while(MeanFragLen > m_TargFragLen)
	{
	RandSiteVal = RGseeds.Random();;
	RandSiteOfs = (size_t)(m_CurTotTransLen * RandSiteVal);			// random selection of potential frgamentation site
	if(RandSiteOfs >= (m_CurTotTransLen - StericOcc3Len) ||
		RandSiteOfs < (size_t)StericOcc5Len)
		continue;
	pBase = &m_pTranscripts[RandSiteOfs];
	Flgs = *pBase;
	if(Flgs & (cFragStericOcc5Flag | cFragStericOcc3Flag))		// if already sterically occluded then need to try different fragmentation site
		continue;

	*pBase++ = (cFragStericOcc5Flag | cFragStartFlag);
	for(OccIdx = 1; OccIdx < StericOcc5Len && !(*pBase & cFragStericOcc3Flag); OccIdx++)
		*pBase++ |= cFragStericOcc5Flag;
	pBase = &m_pTranscripts[RandSiteOfs-1];
	*pBase-- = (cFragStericOcc3Flag | cFragEndFlag);
	for(OccIdx = 1; OccIdx < StericOcc3Len && !(*pBase & cFragStericOcc5Flag); OccIdx++)
		*pBase-- = cFragStericOcc3Flag;
	
	TotalNumFrags += 1;
	MeanFragLen = (int)(m_CurTotTransLen / TotalNumFrags);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Simulated %d fragments with mean fragment length of %d",TotalNumFrags,MeanFragLen);
// determine total number of fragments within user specified range
m_NumRegionsParsed = 0;
TotFragsInRange = 0;
do {
	pTransRegion = &m_pTransRegions[m_NumRegionsParsed++];
	pSeq = &m_pTranscripts[pTransRegion->TransSeqOfs];
	ExpTotSeqLen = pTransRegion->NumInsts * pTransRegion->TransLen;
	CurTotSeqLen = 0;

	do {
		CurFragLen = 0;
		if(!(*pSeq++ & cFragStartFlag))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Internal processing error, expected fragment start flag");
			Reset(false);
			return(eBSFerrInternal);
			}
		CurFragLen = 1;
		while(!(*pSeq++ & cFragEndFlag))
			CurFragLen += 1;
		CurFragLen += 1;
		if(CurFragLen >= MinFragLen && CurFragLen <= m_MaxFragLen)
			TotFragsInRange += 1;
		CurTotSeqLen += CurFragLen;
		if(CurFragLen > 10000)
			CurFragLen = 10000;
		LenDists[CurFragLen / 10] += 1;
		}
	while(CurTotSeqLen < ExpTotSeqLen);
	}
while(m_NumRegionsParsed < m_NumTransRegions);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Fragmentation completed, filtering on length range and reporting...");

// iterate over generated fragments and report on those within the user requested length bounds
if((Rslt=m_pCSV->Open(m_szTransFile))!=eBSFSuccess)
	{
	while(m_pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",m_szTransFile);
	Reset(false);
	return(Rslt);
	}


m_NumRegionsParsed = 0;
while((Rslt=m_pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	NumFields = m_pCSV->GetCurFields();
	if(!m_NumRegionsParsed && (NumFields >= cExpNumInputFields && m_pCSV->IsLikelyHeaderLine())) // expecting at least cExpNumInputFields fields so only then check for header line
		continue;

	if(NumFields < cExpNumInputFields)	// must have as a minimum the  <FeatID>,<FeatName>,<FeatLen>,<FeatInsts>
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected at least %d fields in '%s', GetCurFields() returned '%d'",cExpNumInputFields,m_szTransFile,NumFields);
		delete m_pCSV;
		Reset(false);
		return(eBSFerrFieldCnt);
		}

	m_pCSV->GetInt(1,&ElID);
	m_pCSV->GetText(2,&pszElName);
	m_pCSV->GetInt(3,&ElLen);
	m_pCSV->GetInt(4,&ElInsts);

	if(ElLen < 20 || ElLen > 200000 || ElInsts < 1 || ElInsts > 1000000)
		continue;

	if(((size_t)ElLen * ElInsts) > (size_t)4000000000)
		continue;

	pTransRegion = &m_pTransRegions[m_NumRegionsParsed++];

	pSeq = &m_pTranscripts[pTransRegion->TransSeqOfs];
	ExpTotSeqLen = pTransRegion->NumInsts * pTransRegion->TransLen;
	CurTotSeqLen = 0;
	CurFragsInRange = 0;
	CurFragsUnderLen = 0;
	CurFragsOverLen = 0;
	do {
		CurFragLen = 0;
		if(!(*pSeq++ & cFragStartFlag))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Internal processing error, expected fragment start flag");
			Reset(false);
			return(eBSFerrInternal);
			}
		CurFragLen = 1;
		while(!(*pSeq++ & cFragEndFlag))
			CurFragLen += 1;
		CurFragLen += 1;

		if(CurFragLen < MinFragLen)
			CurFragsUnderLen += 1;
		else
			{
			if(CurFragLen > MaxFragLen)
				CurFragsOverLen += 1;
			else
				CurFragsInRange += 1;
			}
		CurTotSeqLen += CurFragLen;
		}
	while(CurTotSeqLen < ExpTotSeqLen);
	RsltBuffIdx += sprintf(&szRsltsBuff[RsltBuffIdx],"%d,\"%s\",%d,%d,%d,%d,%d\n",ElID,pszElName,ElLen,ElInsts,CurFragsInRange,CurFragsUnderLen,CurFragsOverLen);
	if(RsltBuffIdx > (sizeof(szRsltsBuff) - 500))
		{
		if(write(m_hFragsFile,szRsltsBuff,RsltBuffIdx)!=RsltBuffIdx)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fragmentate: Write %d bytes to file '%s' failed - %s",RsltBuffIdx,pszFragsFile,strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		RsltBuffIdx = 0;
		}
	}
if(RsltBuffIdx > 0)
	{
	if(write(m_hFragsFile,szRsltsBuff,RsltBuffIdx)!=RsltBuffIdx)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fragmentate: Write %d bytes to file '%s' failed - %s",RsltBuffIdx,pszFragsFile,strerror(errno));
		Reset();
		return(eBSFerrFileAccess);
		}
	RsltBuffIdx = 0;
	}
m_pCSV->Close();

#ifdef _WIN32
_commit(m_hFragsFile);
#else
fsync(m_hFragsFile);
#endif
close(m_hFragsFile);
m_hFragsFile = -1;


// now write out the overall length histogram
for(int Idx = 0; Idx < 1001; Idx++)
	{
	if(Idx < 1001)
		RsltBuffIdx += sprintf(&szRsltsBuff[RsltBuffIdx],"%d,%d,%d\n",Idx*10,(Idx*10)+9,LenDists[Idx]);
	else
		RsltBuffIdx += sprintf(&szRsltsBuff[RsltBuffIdx],"10000,100000,%d\n",LenDists[Idx]);

	if(RsltBuffIdx > (sizeof(szRsltsBuff) - 500))
		{
		if(write(m_hDistFile,szRsltsBuff,RsltBuffIdx)!=RsltBuffIdx)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fragmentate: Write %d bytes to file '%s' failed - %s",RsltBuffIdx,pszDistFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		RsltBuffIdx = 0;
		}
	}

if(RsltBuffIdx > 0)
	{
	if(write(m_hDistFile,szRsltsBuff,RsltBuffIdx)!=RsltBuffIdx)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fragmentate: Write %d bytes to file '%s' failed - %s",RsltBuffIdx,pszDistFile,strerror(errno));
		Reset();
		return(eBSFerrFileAccess);
		}
	RsltBuffIdx = 0;
	}

#ifdef _WIN32
_commit(m_hDistFile);
#else
fsync(m_hDistFile);
#endif
close(m_hDistFile);
m_hDistFile = -1;

Reset(true);
return(eBSFSuccess);
}


int
CTranscriptome::AddTranscripts(int RegionID, // uniquely identifies transcribed region 
			int TransLen,				// transcripts are of this length 
			int NumInsts)				// and there are this many instances
{
int Inst;
int SeqOfs;
size_t MinConcatSeqLen;
tsTransRegion *pRegion;
UINT8 *pTransSeqs;
UINT8 Flgs;

if(m_pTransRegions == NULL)
	{
	m_MemAllocdTransRegions = cTransRegionAllocNum * sizeof(tsTransRegion);
#ifdef _WIN32
	m_pTransRegions = (tsTransRegion *) malloc(m_MemAllocdTransRegions);	// initial and perhaps the only allocation

	if(m_pTransRegions == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTranscripts: Region memory allocation of %Iu bytes - %s",(INT64)m_MemAllocdTransRegions,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pTransRegions = (tsTransRegion *)mmap(NULL,m_MemAllocdTransRegions, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pTransRegions == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTranscripts: Region memory allocation of %lld bytes through mmap()  failed - %s",(INT64)m_MemAllocdTransRegions,strerror(errno));
		m_pTransRegions = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_NumAllocdTransRegions = cTransRegionAllocNum;
	m_NumTransRegions = 0;
	}

MinConcatSeqLen = (((size_t)TransLen + 10) * NumInsts);	// scaleup slightly so as to reduce reallocs which may be required
if(m_pTranscripts == NULL)
	{
	m_MemAllocdTrnscripts = max((size_t)cTransSeqAllocSize,MinConcatSeqLen);

#ifdef _WIN32
	m_pTranscripts = (UINT8 *) malloc(m_MemAllocdTrnscripts);	// initial and perhaps the only allocation

	if(m_pTranscripts == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTranscripts: Sequence memory allocation of %Iu bytes - %s",(INT64)m_MemAllocdTrnscripts,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pTranscripts = (UINT8 *)mmap(NULL,m_MemAllocdTrnscripts, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pTranscripts == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTranscripts: Sequence allocation of %lld bytes through mmap()  failed - %s",(INT64)m_MemAllocdTrnscripts,strerror(errno));
		m_pTranscripts = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_MemAllocdTrnscripts = m_MemAllocdTrnscripts;
	m_CurTotTransLen = 0;
	m_NumTranscripts = 0;
	}


// need to reallocate to hold more transcribed regions?
if((m_NumTransRegions + 10)  >= m_NumAllocdTransRegions)
	{
	size_t memreq = m_MemAllocdTransRegions + (cTransRegionAllocNum * sizeof(tsTransRegion));
#ifdef _WIN32
	pRegion = (tsTransRegion *) realloc(m_pTransRegions,memreq);
	if(pRegion == NULL)
		{
#else
	pRegion = (tsTransRegion *)mremap(m_pTransRegions,m_MemAllocdTransRegions,memreq,MREMAP_MAYMOVE);
	if(pRegion == MAP_FAILED)
		{
		pRegion = NULL;
#endif
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTranscripts: Region memory reallocation to %Iu bytes failed - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pTransRegions = pRegion;
	m_MemAllocdTransRegions = memreq;
	m_NumAllocdTransRegions += cTransRegionAllocNum;
	}

// need to reallocate to hold more transcript sequences?
if((m_CurTotTransLen + MinConcatSeqLen)  >= m_MemAllocdTrnscripts)
	{
	size_t memreq = m_MemAllocdTrnscripts + max((size_t)cTransSeqAllocSize,MinConcatSeqLen);
#ifdef _WIN32
	pTransSeqs = (UINT8 *) realloc(m_pTranscripts,memreq);
	if(pTransSeqs == NULL)
		{
#else
	pTransSeqs = (UINT8 *)mremap(m_pTranscripts,m_MemAllocdTrnscripts,memreq,MREMAP_MAYMOVE);
	if(pTransSeqs == MAP_FAILED)
		{
		pTransSeqs = NULL;
#endif
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddTranscripts: Sequence memory reallocation to %Iu bytes failed - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_pTranscripts = pTransSeqs;
	m_MemAllocdTrnscripts = memreq;
	}

pRegion = &m_pTransRegions[m_NumTransRegions++];
pRegion->ElIdx = RegionID;
pRegion->TransLen = TransLen;
pRegion->NumInsts = NumInsts;
pRegion->TransSeqOfs = m_CurTotTransLen;
pTransSeqs = &m_pTranscripts[m_CurTotTransLen];

for(Inst = 0; Inst < NumInsts; Inst++)
	{
	m_NumTranscripts += 1;
	for(SeqOfs = 0; SeqOfs < TransLen; SeqOfs++)
		{
		Flgs = 0;						    // this is an unfragmented sequence base
		if(SeqOfs == 0)
			Flgs |= cFragStartFlag;			// this base starts the transcript
		else
			if(SeqOfs == (TransLen - 1))
				Flgs |= cFragEndFlag;		// this base ends the transcript
		if(SeqOfs < m_StericOcc5Len)
			Flgs |= cFragStericOcc5Flag;			// this base is part of the 5' steric occlusion	
		if((SeqOfs + m_StericOcc3Len) >= TransLen)
			Flgs |= cFragStericOcc3Flag;			// this base is part of the 3' steric occlusion			
		*pTransSeqs++ = Flgs;
		m_CurTotTransLen += 1;
		}
	}

return(eBSFSuccess);
}