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

#include "./ChromFeatures.h"

CChromFeatures::CChromFeatures(void)
{
m_pBedFile = NULL;
m_NumChroms = 0;
m_hOutFile = -1;						// file handle for generated BED file to contain merged features
Reset();
}

CChromFeatures::~CChromFeatures(void)
{
tsChromFeatures *pChrom;
int Idx;
if(m_hOutFile != -1)
	close(m_hOutFile);
if(m_pBedFile != NULL)
	delete m_pBedFile;
if(m_NumChroms > 0)
	{
	pChrom = &m_Chroms[0];
	for(Idx = 0; Idx < m_NumChroms; Idx++,pChrom++)
		{
		if(pChrom->pFeatures != NULL)
			{
#ifdef _WIN32
			free(pChrom->pFeatures);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
			if(pChrom->pFeatures != MAP_FAILED)
				munmap(pChrom->pFeatures,pChrom->NumFeaturesAllocMem);
#endif
			}	
		}
	}
}

void
CChromFeatures::Reset(void)
{
int Idx;
tsChromFeatures *pChrom;
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}
if(m_pBedFile != NULL)
	{
	delete m_pBedFile;
	m_pBedFile = NULL;
	}
if(m_NumChroms > 0)
	{
	pChrom = &m_Chroms[0];
	for(Idx = 0; Idx < m_NumChroms; Idx++,pChrom++)
		{
		if(pChrom->pFeatures != NULL)
			{
#ifdef _WIN32
			free(pChrom->pFeatures);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
			if(pChrom->pFeatures != MAP_FAILED)
				munmap(pChrom->pFeatures,pChrom->NumFeaturesAllocMem);
#endif
			pChrom->pFeatures = NULL;
			}	
		}
	}
memset(m_Chroms,0,sizeof(m_Chroms));
m_NumChroms = 0;
m_pCurBedFile = NULL;
m_pCurLocChrom = NULL;
m_NumBedFiles = 0;
szOutFile[0] = '\0';
m_LineBuffIdx = 0;
m_NumIncludeChroms = 0;
m_NumExcludeChroms = 0;
m_ppszIncludeChroms = NULL;
m_ppszExcludeChroms = NULL;
m_szFiltChrom[0]='\0';
m_bFiltChrom = false;
}

int
CChromFeatures::SetChromFilters(int NumIncludeChroms,		 // number of chromosome regular expressions to include
						char **ppszIncludeChroms,	 // array of include chromosome regular expressions
						int NumExcludeChroms,		 // number of chromosome expressions to exclude
						char **ppszExcludeChroms)	 // array of exclude chromosome regular expressions
{
int Idx;

#ifndef _WIN32
int RegErr;				// regular expression parsing error
char szRegErr[128];			// to hold RegErr as textual representation ==> regerror();
#endif

m_NumIncludeChroms = NumIncludeChroms;		// number of chromosomes explicitly defined to be included
m_ppszIncludeChroms = ppszIncludeChroms;	// ptr to array of reg expressions defining chroms to include - overides exclude
m_NumExcludeChroms = NumExcludeChroms;		// number of chromosomes explicitly defined to be excluded
m_ppszExcludeChroms = ppszExcludeChroms;	// ptr to array of reg expressions defining chroms to include

if(NumIncludeChroms == 0 && NumExcludeChroms == 0)
	return(eBSFSuccess);

#ifdef _WIN32
try {
	for(Idx=0;Idx < NumIncludeChroms;Idx++)
		{
		m_IncludeChromsRE[Idx] = new Regexp();
		m_IncludeChromsRE[Idx]->Parse(ppszIncludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include regexpr chrom '%s'",ppszIncludeChroms[Idx]);
	Reset();
	return(eBSFerrMem);
	}
try {
	for(Idx=0;Idx < NumExcludeChroms;Idx++)
		{
		m_ExcludeChromsRE[Idx] = new Regexp();
		m_ExcludeChromsRE[Idx]->Parse(ppszExcludeChroms[Idx],false);	// note case insensitive
		}
	}
catch(...)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude regexpr chrom '%s'",ppszExcludeChroms[Idx]);
	Reset();
	return(eBSFerrMem);
	}

#else
for(Idx=0;Idx < NumIncludeChroms;Idx++)
	{

	RegErr=regcomp(&m_IncludeChromsRE[Idx],ppszIncludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&m_IncludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process include chrom '%s' error: %s",ppszIncludeChroms[Idx],szRegErr);
		Reset();
		return(eBSFerrMem);
		}
	}
for(Idx=0;Idx < NumExcludeChroms;Idx++)
	{
	RegErr = regcomp(&m_ExcludeChromsRE[Idx],ppszExcludeChroms[Idx],REG_EXTENDED | REG_ICASE);	// note case insensitive
	if(RegErr)
		{
		regerror(RegErr,&m_ExcludeChromsRE[Idx],szRegErr,sizeof(szRegErr));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to process exclude chrom '%s' error: %s",ppszExcludeChroms[Idx],szRegErr);
		Reset();
		return(eBSFerrMem);
		}
	}
#endif
return(eBSFSuccess);
}

// ExcludeThisChrom
// Returns true if pszChrom is to be excluded from processing
bool
CChromFeatures::ExcludeThisChrom(char *pszChrom)
{
#ifdef _WIN32
RegexpMatch mc;
#else
regmatch_t mc;
#endif

int Idx;
if(!m_NumExcludeChroms && !m_NumIncludeChroms)
	return(false);

if(m_szFiltChrom[0] != 0)
	{
	if(!stricmp(m_szFiltChrom,pszChrom))
		return(m_bFiltChrom);
	}
strcpy(m_szFiltChrom,pszChrom);
m_bFiltChrom = false;

// explicitly to be excluded?
for(Idx = 0; Idx < m_NumExcludeChroms; Idx++)
#ifdef _WIN32	
	if(m_ExcludeChromsRE[Idx]->Match(pszChrom,&mc))
#else
	if(!regexec(&m_ExcludeChromsRE[Idx],pszChrom,1,&mc,0))
#endif
		return(m_bFiltChrom = true);

// explicitly to be included?
for(Idx = 0; Idx < m_NumIncludeChroms; Idx++)
	{
#ifdef _WIN32
	if(m_IncludeChromsRE[Idx]->Match(pszChrom,&mc))
#else
	if(!regexec(&m_IncludeChromsRE[Idx],pszChrom,1,&mc,0))
#endif
		return(m_bFiltChrom = false);
	}


// if chromosomes were defined as to explicitly include then this chrom is to be filtered out
m_bFiltChrom = m_NumIncludeChroms > 0 ? true : false;
return(m_bFiltChrom);
}


int
CChromFeatures::AddFeature(char *pszBedFile,	// feature is from this BED file
						   char *pszChrom,		// and is on this chromosome
						   int StartLoci,		// starts at this loci
						   int EndLoci,			// and finishes at this loci
						   char Strand)			// and is on this strand
{
int Idx;
size_t memreq;
tsChromFeatures *pChrom;
tsFeatureEl *pFeature;
tsBED *pBed;

if(m_pCurBedFile == NULL ||
   stricmp(m_pCurBedFile->szBedFile,pszBedFile))
	{
	m_pCurBedFile = NULL;
	if(m_NumBedFiles > 0)
		{
		pBed = &m_BedFiles[0];
		for(Idx = 0; Idx < m_NumBedFiles; Idx++)
			{
			if(!stricmp(pBed->szBedFile,pszBedFile))
				{
				m_pCurBedFile = pBed;
				break;
				}
			}
		}
	if(m_pCurBedFile == NULL)
		{
		if(m_NumBedFiles == cMaxNumBedFiles)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddFeature: too many BED files");
			Reset();
			return(eBSFerrMaxEntries);
			}
		m_pCurBedFile = &m_BedFiles[m_NumBedFiles++];
		strcpy(m_pCurBedFile->szBedFile,pszBedFile);
		m_pCurBedFile->BedID = m_NumBedFiles;
		m_pCurBedFile->NumEls = 0;
		}
	}

if(m_pCurLocChrom != NULL &&
   !stricmp(m_pCurLocChrom->szChrom,pszChrom))
	{
	pChrom = m_pCurLocChrom;
	if(pChrom->StartLoci == -1 || pChrom->StartLoci > StartLoci)
		pChrom->StartLoci = StartLoci;
	if(pChrom->EndLoci < EndLoci)
		pChrom->EndLoci = EndLoci;
	}
else
	{
	// iterate previously loaded chroms and see if chrom previously added
	if(m_NumChroms > 0)
		{
		pChrom = &m_Chroms[0];
		for(Idx= 0; Idx < m_NumChroms; Idx++,pChrom++)
			if(!stricmp(pChrom->szChrom,pszChrom))
				break;
		}
	if(m_NumChroms == 0 || Idx ==  m_NumChroms)
		{
		pChrom = &m_Chroms[m_NumChroms++];
		memset(pChrom,0,sizeof(tsChromFeatures));
		pChrom->StartLoci = -1;
		pChrom->ChromID = m_NumChroms;
		strncpy(pChrom->szChrom,pszChrom,sizeof(pChrom->szChrom));
		pChrom->szChrom[sizeof(pChrom->szChrom)-1] = '\0';
		pChrom->NumFeatures = 0;
		pChrom->NumFeaturesAlloc = 0;
		pChrom->NumFeaturesAllocMem = 0;
		pChrom->pFeatures = NULL;
		memreq = 2 * cNumAllocEls * sizeof(tsFeatureEl);

#ifdef _WIN32
		pFeature = (tsFeatureEl *) malloc(memreq);	// initial and perhaps the only allocation
		if(pFeature == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddFeature: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
			Reset();
			return(eBSFerrMem);
			}
#else
		// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
		pFeature = (tsFeatureEl *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
		if(pFeature == MAP_FAILED)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddFeature: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
			Reset();
			return(eBSFerrMem);
			}
#endif
		pChrom->pFeatures = pFeature;
		pChrom->NumFeaturesAllocMem = memreq;;
		pChrom->NumFeaturesAlloc = cNumAllocEls;
		}
	if(pChrom->StartLoci == -1 || pChrom->StartLoci > StartLoci)
		pChrom->StartLoci = StartLoci;
	if(pChrom->EndLoci < EndLoci)
		pChrom->EndLoci = EndLoci;
	m_pCurLocChrom = pChrom;
	}

// add feature to this chrom
if((pChrom->NumFeatures + 10) >=  pChrom->NumFeaturesAlloc)	// 10 is to give a little float
	{
	memreq = 2 * (pChrom->NumFeaturesAlloc + cNumAllocEls) * sizeof(tsFeatureEl);
#ifdef _WIN32
	pFeature = (tsFeatureEl *) realloc(pChrom->pFeatures,memreq);
#else
	pFeature = (tsFeatureEl *)mremap(pChrom->pFeatures,pChrom->NumFeaturesAllocMem,memreq,MREMAP_MAYMOVE);
	if(pFeature == MAP_FAILED)
		pFeature = NULL;
#endif
	if(pFeature == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadReads: Memory re-allocation to %d bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	pChrom->pFeatures = pFeature;
	pChrom->NumFeaturesAlloc += cNumAllocEls;
	pChrom->NumFeaturesAllocMem = memreq;
	}

pFeature = &pChrom->pFeatures[pChrom->NumFeatures*2];
pChrom->NumFeatures += 1;
pFeature->BedID = m_pCurBedFile->BedID;
pFeature->FeatureID = pChrom->NumFeatures;
pFeature->ChromID = pChrom->ChromID;
pFeature->Loci = StartLoci;
pFeature->Strand = Strand;
pFeature->FlgStart = 1;
pFeature->FlgEnd = 0;
pFeature += 1;
pFeature->BedID = m_pCurBedFile->BedID;
pFeature->FeatureID = pChrom->NumFeatures;
pFeature->ChromID = pChrom->ChromID;
pFeature->Loci = EndLoci;
pFeature->Strand = Strand;
pFeature->FlgStart = 0;
pFeature->FlgEnd = 1;
m_pCurBedFile->NumEls += 1;

return(eBSFSuccess);
}

tsChromFeatures *
CChromFeatures::LocateChrom(char *pszChrom)
{
int Idx;
tsChromFeatures *pChrom;

if(m_pCurLocChrom != NULL &&
   !stricmp(m_pCurLocChrom->szChrom,pszChrom))
	return(m_pCurLocChrom);
pChrom = &m_Chroms[0];
for(Idx = 0; Idx < m_NumChroms; Idx++,pChrom++)
	if(!stricmp(pChrom->szChrom,pszChrom))
		{
		m_pCurLocChrom = pChrom;
		return(pChrom);
		}
return(NULL);
}

teBSFrsltCodes 
CChromFeatures::LoadBED(etBEDRegion Region,		// functional region of interest
				char ROIStrand,					// region on this strand
				 char *pszInFile)				// UCSC BED containing regions
{
int Rslt;
int CurFeatureID;
int NumExons;
int NumIntrons;
int CDSstart;
int CDSend;
int StartLoci;
int EndLoci;
int Score;
char szChrom[128];
char szPrevChrom[128];
char szFeatName[128];
char Strand;
int RefID;
int IntergenicStart;
int Idx;
int NumEls;
tsChromFeatures *pChrom;

if((m_pBedFile = new CBEDfile)==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CBEDfile");
	Reset();
	return(eBSFerrObj);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to load %s",pszInFile);
if((Rslt=m_pBedFile->Open(pszInFile,eBTAnyBed)) !=eBSFSuccess)
	{
	while(m_pBedFile->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBedFile->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open '%s' for processing",pszInFile);
	Reset();
	return(eBSFerrOpnFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Load completed, iterating features...");
if(!m_pBedFile->ContainsGeneDetail() && Region > eMEGRIntergenic)			// returns true if file contains gene detail (utr/cds/intron etc)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"filter bed file %s does not contain gene features",pszInFile);
	Reset();
	return(eBSFerrFeature);
	}

// now iterate over the features, filtering as may be appropriate
szPrevChrom[0] = '\0';
CurFeatureID = 0;
IntergenicStart = 0;
RefID = 0;
Rslt = eBSFSuccess;
NumEls = 0;
while(Rslt == eBSFSuccess && (CurFeatureID = m_pBedFile->GetNextFeatureID(CurFeatureID)) > 0)
	{
	NumEls += 1;
	m_pBedFile->GetFeature(CurFeatureID,	// feature instance identifier
				szFeatName,				// where to return feature name
				szChrom,				// where to return chromosome name
				&StartLoci,				// where to return feature start on chromosome (0..n) 
				&EndLoci,				// where to return feature end on chromosome
 				&Score,					// where to return score
 				&Strand);				// where to return strand

	if(CurFeatureID == 1 || stricmp(szChrom,szPrevChrom))	// if new chromosome then reset IntergenicStart
		{
		strcpy(szPrevChrom,szChrom);
		IntergenicStart = 0;
		}

	// exclude this chrom?

	if(ExcludeThisChrom(szChrom))
		continue;

	if(ROIStrand != '*')
		{
		if(ROIStrand != Strand)
			continue;
		}
	
	if(Region != eMEGRAny)
		{
		NumExons = m_pBedFile->GetNumExons(CurFeatureID);					// returns number of exons - includes UTRs + CDS
		NumIntrons = m_pBedFile->GetNumIntrons(CurFeatureID);
		CDSstart = StartLoci + m_pBedFile->GetCDSStart(CurFeatureID);		// returns relative start offset of CDS - NOTE add to '+' strand gene start, subtract on '-' strand gene start
		CDSend = StartLoci + m_pBedFile->GetCDSEnd(CurFeatureID);			// returns relative end offset of CDS - NOTE add to '+' strand gene start, subtract on '-' strand gene start
		}
	Rslt = eBSFSuccess;

	switch(Region) {
		case eMEGRAny:			// retain any region
			Rslt=AddFeature(pszInFile,szChrom,StartLoci,EndLoci,Strand);
			continue;

		case eMEGRIntergenic:	// only retain intergenic
			if(IntergenicStart < StartLoci)
				Rslt=AddFeature(pszInFile,szChrom,IntergenicStart,StartLoci-1,Strand);
			if(IntergenicStart <= EndLoci)
				IntergenicStart = EndLoci+1;
			continue;

		case eMEGRExons:		// only retain exons
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBedFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBedFile->GetExonEnd(CurFeatureID,Idx);
				if(StartLoci <= EndLoci)
					Rslt=AddFeature(pszInFile,szChrom,StartLoci,EndLoci,Strand);
				}
			continue;

		case eMEGRIntrons:		// only retain introns
			for(Idx = 1; Idx <= NumIntrons; Idx++)
				{
				StartLoci = m_pBedFile->GetIntronStart(CurFeatureID,Idx);
				EndLoci = m_pBedFile->GetIntronEnd(CurFeatureID,Idx);
				if(StartLoci <= EndLoci)
					Rslt=AddFeature(pszInFile,szChrom,StartLoci,EndLoci,Strand);
				}
			continue;

		case eMEGRCDS:			// only retain CDSs
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBedFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBedFile->GetExonEnd(CurFeatureID,Idx);
				if(EndLoci < CDSstart || StartLoci > CDSend)
					continue;
				if(StartLoci < CDSstart)
					StartLoci = CDSstart;
				if(EndLoci > CDSend)
					EndLoci = CDSend;
				if(StartLoci <= EndLoci)
					AddFeature(pszInFile,szChrom,StartLoci,EndLoci,Strand);
				}
			continue;

		case eMEGUTR:			// only process UTRs - single exon may have both 5' and 3' UTRs
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBedFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBedFile->GetExonEnd(CurFeatureID,Idx);
				if(EndLoci <= CDSend && StartLoci >= CDSstart) // is exon CDS only?
					continue;

				// check if 5' UTR 
				if(StartLoci < CDSstart)
					{
					if(EndLoci >= CDSstart)
						Rslt=AddFeature(pszInFile,szChrom,StartLoci,CDSstart-1,Strand);
					else
						Rslt=AddFeature(pszInFile,szChrom,StartLoci,EndLoci,Strand);
					}
					
				// check if 3'UTR
				if(EndLoci > CDSend)
					{
					if(StartLoci <= CDSend)
						StartLoci = CDSend+1;
					Rslt=AddFeature(pszInFile,szChrom,StartLoci,EndLoci,Strand);
					}
				}
			continue;

		case eMEG5UTR:			// only process 5'UTRs - strand sensitive
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBedFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBedFile->GetExonEnd(CurFeatureID,Idx);
				if(EndLoci <= CDSend && StartLoci >= CDSstart) // is exon CDS only?
					continue;

				if(Strand != '-')
					{
					// check if 5' UTR on '+' strand 
					if(StartLoci < CDSstart)
						{
						if(EndLoci >= CDSstart)
							EndLoci = CDSstart - 1;
						}
					}
				else 
					{	
					// check if 5'UTR on '-' strand
					if(EndLoci > CDSend)
						{
						if(StartLoci <= CDSend)
							StartLoci = CDSend+1;
						}
					}
				Rslt=AddFeature(pszInFile,szChrom,StartLoci,EndLoci,Strand);
				}
			continue;

		case eMEG3UTR:			// only process 3'UTRs  - strand sensitive
			for(Idx = 1; Idx <= NumExons; Idx++)
				{
				StartLoci = m_pBedFile->GetExonStart(CurFeatureID,Idx);
				EndLoci   = m_pBedFile->GetExonEnd(CurFeatureID,Idx);
				if(EndLoci <= CDSend && StartLoci >= CDSstart) // is exon CDS only?
					continue;

				if(Strand == '-')
					{
					// check if 3' UTR on '-' strand 
					if(StartLoci < CDSstart)
						{
						if(EndLoci >= CDSstart)
							EndLoci = CDSstart - 1;
						}
					}
				else 
					{	
					// check if 3'UTR on '+' strand
					if(EndLoci > CDSend)
						{
						if(StartLoci <= CDSend)
							StartLoci = CDSend+1;
						}
					}
				Rslt=AddFeature(pszInFile,szChrom,StartLoci,EndLoci,Strand);
				}
			continue;
		}
	}
pChrom = &m_Chroms[0];
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processed %d features",NumEls);
return((teBSFrsltCodes)Rslt);
}


void 
CChromFeatures::DumpFeat(char *pszMsg,tsFeatureEl *pFeat)
{
char szBuff[200];
sprintf(szBuff,"%s BedID: %d FeatID: %d ChrID: %d Strand: %c Loci: %d FlgStart: %d FlgEnd: %d\n",pszMsg,
														pFeat->BedID,				// BED file from which this feature was loaded
														pFeat->FeatureID,			// feature start/ends share same feature identifier
														pFeat->ChromID,				// chrom on which this feature is located
														pFeat->Strand,				// feature is on this strand
														pFeat->Loci,				// feature starts or ends at this loci (0..N)
														pFeat->FlgStart,			// will be 1 if feature start
														pFeat->FlgEnd);				// will be 1 if feature end
gDiagnostics.DiagOut(eDLInfo,gszProcName,szBuff);
}

int
CChromFeatures::ReportMergedFeat(UINT32 FeatID,				// uniquely identifies this merged feature
								 char *pszChrom,			// feature is on this chrom
								 UINT32 MergeStartLoci,		// starting at this loci
								 UINT32 MergeEndLoci,		// ending at this loci
								 char Strand)				// and is on this strand
{
if((m_LineBuffIdx + 100) > sizeof(m_LineBuff))
	{
	CUtility::SafeWrite(m_hOutFile,m_LineBuff,m_LineBuffIdx);
	m_LineBuffIdx = 0;
	}

m_LineBuffIdx += sprintf(&m_LineBuff[m_LineBuffIdx],"%s\t%d\t%d\tMRG%d\t0\t%c\n",pszChrom,MergeStartLoci,MergeEndLoci+1,FeatID,Strand);

return(0);
}

// merge features which overlap, or at most have a maximum separation
int 
CChromFeatures::MergeFeatures(char *pszOutFile,			// file to which merged features are to be written
							  bool bStrandSpecific,		// if true then merging is strand specific
							  int MaxSep,				// merge features which are separated by up to this separation, if 0 then don't merge
							  int MinLen)				// only report merged features if they span at least this distance
{
int StrandIdx;					// current strand index into ElDepth, MergeStartLoci and MergeEndLoci
int ElDepth[2];					// current overlap depth for each strand
UINT32 MergeStartLoci[2];		// current merge start for each strand
UINT32 MergeEndLoci[2];			// current merge end for each strand
bool bJoinFeats[2];				// true if features are separated by at most MaxSep and thus are being joined
char Strand;					// report merge as being on this strand

UINT32 NumMerged;
int ChromIdx;
int FeatIdx;
tsChromFeatures *pChrom;
tsFeatureEl *pFeat;


int IdxSep;
tsFeatureEl *pFeatSep;


#ifdef _WIN32
if((m_hOutFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutFile = open(pszOutFile, O_RDWR | O_CREAT,S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}

// next sort the features to be merged by ascending chrom, loci
gDiagnostics.DiagOut(eDLInfo,gszProcName,"now sorting features...");
pChrom = &m_Chroms[0];
for(ChromIdx = 0; ChromIdx < m_NumChroms; ChromIdx++,pChrom++)
	qsort(pChrom->pFeatures,2*pChrom->NumFeatures,sizeof(* pChrom->pFeatures),SortFeatures);

// features sorted, now start the merging..
gDiagnostics.DiagOut(eDLInfo,gszProcName,"now merging features...");
NumMerged = 0;
pChrom = &m_Chroms[0];
for(ChromIdx = 0; ChromIdx < m_NumChroms; ChromIdx++,pChrom++)
	{
	pFeat = pChrom->pFeatures;
	ElDepth[0] = 0;
	ElDepth[1] = 0;
	bJoinFeats[0] = false;
	bJoinFeats[1] = false;
	for(FeatIdx = 0; FeatIdx < (2 * pChrom->NumFeatures); FeatIdx++,pFeat++)
		{
		if(bStrandSpecific)			// set if merging is strand specific
			{
			Strand = pFeat->Strand ;
			StrandIdx = pFeat->Strand == '+' ? 0 : 1;
			}
		else
			{
			Strand = '+';
			StrandIdx = 0;
			}
		if(pFeat->FlgStart == 1)	// if starting a feature
			{
			if(ElDepth[StrandIdx] == 0)
				{
				MergeStartLoci[StrandIdx] = pFeat->Loci;
				ElDepth[StrandIdx] = 1;
				}
			else
				if(bJoinFeats[StrandIdx] == false) 
					ElDepth[StrandIdx] += 1;
			bJoinFeats[StrandIdx] = false;
			}
		else 
			if(pFeat->FlgEnd == 1)	// if finish of a feature
				{
				if(ElDepth[StrandIdx] == 1)
					{
					if(MaxSep > 0	// need to look ahead and check if another feature starting on same strand within MaxSep?
						&& FeatIdx < ((2 * pChrom->NumFeatures)-1))
						{
						pFeatSep = pFeat + 1;
						bJoinFeats[StrandIdx] = false;
						for(IdxSep = FeatIdx + 1; IdxSep < (2 * pChrom->NumFeatures); IdxSep++, pFeatSep++)
							{
							if(bStrandSpecific && pFeatSep->Strand != Strand)
								continue;
							if((pFeatSep->Loci - pFeat->Loci) <= (UINT32)MaxSep)
								bJoinFeats[StrandIdx] = true;
							break;
							}
						if(bJoinFeats[StrandIdx])
							continue;
						}
					MergeEndLoci[StrandIdx] = pFeat->Loci;
					if((1 + MergeEndLoci[StrandIdx] - MergeStartLoci[StrandIdx]) >= (UINT32)MinLen)
						ReportMergedFeat(++NumMerged,pChrom->szChrom,MergeStartLoci[StrandIdx],MergeEndLoci[StrandIdx],Strand);
					}
				ElDepth[StrandIdx] -= 1;
				}
		}
	}
if(m_LineBuffIdx > 0)
	{
	CUtility::SafeWrite(m_hOutFile,m_LineBuff,m_LineBuffIdx);
	m_LineBuffIdx = 0;
	}
close(m_hOutFile);
m_hOutFile = -1;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"merging features completed");
Reset();
return(NumMerged);
}

// SortFeatures
// sort by ascending chrom, loci, start, FlgStart (feature starts)
int
CChromFeatures::SortFeatures(const void *arg1, const void *arg2)
{
tsFeatureEl *pEl1 = (tsFeatureEl *)arg1;
tsFeatureEl *pEl2 = (tsFeatureEl *)arg2;

if(pEl1->ChromID < pEl2->ChromID)
		return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);

if(pEl1->Loci < pEl2->Loci)
	return(-1);
if(pEl1->Loci > pEl2->Loci)
	return(1);

if(pEl1->FlgStart > pEl2->FlgStart)
	return(-1);
if(pEl1->FlgStart < pEl2->FlgStart)
	return(1);
return(0);
}
