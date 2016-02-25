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

#include "biokanga.h"
#include "MarkerSeq.h"

CMarkerSeq::CMarkerSeq(void)
{
m_pHypers = NULL;
m_pFasta = NULL;
m_pSeq = NULL;
m_pMarkerBuff = NULL;
m_hOutMarkerFile = -1;
Reset();
}


CMarkerSeq::~CMarkerSeq(void)
{
Reset();
}

int 
CMarkerSeq::Reset(void)
{
if(m_hOutMarkerFile != -1)
	{
	close(m_hOutMarkerFile);
	m_hOutMarkerFile = -1;
	}
if(m_pHypers != NULL)
	{
	delete m_pHypers;
	m_pHypers = NULL;
	}
if(m_pFasta != NULL)
	{
	delete m_pFasta;
	m_pFasta = NULL;
	}
if(m_pSeq != NULL)
	{
#ifdef _WIN32
	free(m_pSeq);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pSeq != MAP_FAILED)
		munmap(m_pSeq,m_AllocdSeqMem);
#endif
	m_pSeq = NULL;
	}
if(m_pMarkerBuff != NULL)
	{
	delete m_pMarkerBuff;
	m_pMarkerBuff = NULL;
	}
m_PMode = 0;	
m_Ftype = 0;
m_CSVFormat = (teCSVFormat)0;
m_MaxNs = 0;				
m_Extend5= 0;				
m_Extend3= 0;				
m_MinSeqLen=0;		
m_MaxSeqLen=0;
m_bNoOverlaps = false;
m_NumEls = 0;				
m_AllocdSeqMem = 0;		
m_AllocdMarkerBuff = 0;
MarkerBuffOfs = 0;

m_szInLociFile[0] = '\0';		
m_szInFastaFile[0] = '\0';		
m_szOutFile[0] = '\0';			
return(eBSFSuccess);
}

int
CMarkerSeq::ProcessMarkerSeqs(int PMode,		// currently default processing only is supported
					int Ftype,					// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM
					teCSVFormat CSVFormat,		// expected input CSV loci file format
					int MaxNs,					// filter out marker sequences having higher than this number of indeterminate bases per 100bp (default is 0, range 0..5)
					int Extend5,				// extend markers 5' bases from centroid of loci (default is 0, range 0..250)
					int Extend3,				// extend markers 3' bases from centroid of loci (default is 0, range 5..250)
					int MinSeqLen,				// filter out marker sequences which are less than this length (default is 10bp, range 10..1000)
					int MaxSeqLen,				// filter out marker sequences which are longer than this length (default is 1000bp, minseqlen..1000)
					bool bNoOverlaps,			// filter out marker sequences which overlap with other marker sequences
					char *pszInLociFile,		// Loci file specifying the SNP or region loci in assembly from which marker sequences are to be generated (CSV, BED or SAM)
					char *pszInFastaFile,		// multifasta assembly file from which marker sequences are to be generated
					char *pszOutFile)			// marker sequences written to this file
{
int Rslt;
CFasta Fasta;
INT64 FastaSize;
etClassifyFileType FileType;

Reset();

m_PMode = PMode;
m_Ftype = Ftype;
m_CSVFormat = CSVFormat;
m_MaxNs = MaxNs;
m_Extend5= Extend5;
m_Extend3= Extend3;
m_MinSeqLen= MinSeqLen;
m_MaxSeqLen = MaxSeqLen;
m_bNoOverlaps = bNoOverlaps;
strncpy(m_szInLociFile,pszInLociFile,sizeof(m_szInLociFile));
m_szInLociFile[sizeof(m_szInLociFile)-1] = '\0';
strncpy(m_szInFastaFile,pszInFastaFile,sizeof(m_szInFastaFile));
m_szInFastaFile[sizeof(m_szInFastaFile)-1] = '\0';
strncpy(m_szOutFile,pszOutFile,sizeof(m_szOutFile));
m_szOutFile[sizeof(m_szOutFile)-1] = '\0';

// allocate to buffer generated marker sequences
if((m_pMarkerBuff = new UINT8 [cMarkerBuffAlloc])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessMarkerSeqs: Unable to allocate %d bytes memory for marker buffer",cMarkerBuffAlloc);
	Reset();
	return(eBSFerrMem);
	}
m_AllocdMarkerBuff = cMarkerBuffAlloc;
MarkerBuffOfs = 0;

// allocate to hold max sized sequences loaded from input fasta file
UINT32 EstNumSeqs;
if((EstNumSeqs = Fasta.FastaEstSizes(pszInFastaFile,&FastaSize)) == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessMarkerSeqs: Nothing to do, fasta file '%s' seems empty or inaccessable",pszInFastaFile);
	Reset();
	return(1);
	}

size_t memreq;
memreq = (size_t)min(FastaSize,0x05fffffff);			// limited to processing contigous fasta sequences which are less than 1.6G

#ifdef _WIN32
m_pSeq = (UINT8 *) malloc(memreq);	// initial and perhaps the only allocation

if(m_pSeq == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessMarkerSeqs: Memory allocation of %lld bytes - %s",(INT64)memreq,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}
#else
// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
m_pSeq = (UINT8 *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
if(m_pSeq == MAP_FAILED)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessMarkerSeqs: Memory allocation of %lld bytes through mmap()  failed - %s",(INT64)memreq,strerror(errno));
	m_pSeq = NULL;
	Reset();
	return(eBSFerrMem);
	}
#endif
m_AllocdSeqMem = memreq;

// firstly load loci and sort ascending
if((m_pHypers = new CHyperEls) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CHyperEls");
	return(eBSFerrObj);
	}

if(Ftype == 0)
	FileType = CUtility::ClassifyFileType(pszInLociFile);
else
	FileType = (etClassifyFileType)(Ftype - 1);

switch(FileType) {
	case eCFTopenerr:		// unable to open file for reading
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: '%s'",pszInLociFile);
		return(eBSFerrOpnFile);

	case eCFTlenerr:		// file length is insufficent to classify type
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to classify file type (insuffient data points): '%s'",pszInLociFile);
		return(eBSFerrFileAccess);

	case eCFTunknown:		// unable to reliably classify
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Unable to reliably classify file type: '%s'",pszInLociFile);
		return(eBSFerrFileType);

	case eCFTCSV:			// file has been classified as being CSV
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing CSV file (%d..%d): '%s'",1,m_MaxSeqLen,pszInLociFile);
		if((Rslt = m_pHypers->ParseCSVFileElements(pszInLociFile,1,m_MaxSeqLen,CSVFormat)) < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in CSV file (%d..%d): '%s'",1,m_MaxSeqLen,pszInLociFile);
			Reset();
			return(Rslt);
			}
		break;

	case eCFTBED:			// file has been classified as being BED
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing BED file (%d..%d): '%s'",1,m_MaxSeqLen,pszInLociFile);
		if((Rslt = m_pHypers->ParseBEDFileElements(pszInLociFile,1,m_MaxSeqLen)) < 0)
			{	
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in BED file (%d..%d): '%s'",1,m_MaxSeqLen,pszInLociFile);
			Reset();
			return(Rslt);
			}
		break;

	case eCFTSAM:			// file has been classified as being SAM
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsing SAM file (%d..%d): '%s'",1,m_MaxSeqLen,pszInLociFile);
		if((Rslt = m_pHypers->ParseSAMFileElements(pszInLociFile,1,m_MaxSeqLen)) < 0)
			{	
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Parse errors in SAM file (%d..%d): '%s'",1,m_MaxSeqLen,pszInLociFile);
			Reset();
			return(Rslt);
			}
		break;
	}


m_NumEls = m_pHypers->NumEls();
if(m_NumEls == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No elements with length range %d..%d in CSV file: '%s'",1,m_MaxSeqLen,pszInLociFile);
	Reset();
	return(Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loaded and parsed %d elements, removing duplicates",m_NumEls);

m_NumEls = m_pHypers->DedupeSort(0,true);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"After duplicates removed there are %d elements",m_NumEls);

// iterate the assembly sequences
if((m_pFasta = new CFasta())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create CFasta object");
	return(eBSFerrObj);
	}

if((Rslt=m_pFasta->Open(m_szInFastaFile,true)) < eBSFSuccess)
	{
	while(m_pFasta->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pFasta->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open fasta file '%s'",m_szInFastaFile);
	Reset();
	return(Rslt);
	}

// create the output marker file
#ifdef _WIN32
if((m_hOutMarkerFile = open(pszOutFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((m_hOutMarkerFile = open(pszOutFile, O_RDWR | O_CREAT |O_TRUNC, S_IREAD | S_IWRITE))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create or truncate %s - %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrCreateFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Marker output multifasta file created/truncated: '%s'",pszOutFile);

// check if descriptor identifies a sequence which contains marker
int CurSeqLen;
UINT32 CurEntryID;
int ChromID;
UINT8 *pSeqStart;
UINT8 *pSeqEnd;
char *pSrc;
char *pBase;
char *pDst;
char Chr;
char szDescription[cBSFDescriptionSize];
char szChrom[cMaxDatasetSpeciesChrom+1];
char szMarkerSeq[cMaxMarkerSeqLen+1];
bool bDescriptor;
int NumMarkers;
int Num5Filtered;
int Num3Filtered;
int NumUnderlen;
int NumOverlen;
int NumNs;
int NumNfiltered;
int NumMatchingPutativeChroms;
int NumNonMatchingPutativeChroms;
int NumAcceptedMarkerChroms;
tsHyperElement *pHypEl;

int NumOvrlapFiltered;
int PrevMarkerChromID;
int PrevMarkerEnd;

bDescriptor = false;
CurEntryID = 0;
NumMarkers = 0;

Num5Filtered = 0;
Num3Filtered = 0;
NumUnderlen = 0;
NumOverlen = 0;
NumNs = 0;
NumNfiltered = 0;
NumMatchingPutativeChroms = 0;
NumNonMatchingPutativeChroms = 0;
NumAcceptedMarkerChroms = 0;
NumOvrlapFiltered = 0;
PrevMarkerChromID = 0;
PrevMarkerEnd = 0;

while((Rslt = CurSeqLen = m_pFasta->ReadSequence(m_pSeq,(int)(m_AllocdSeqMem-1))) > eBSFSuccess)
	{
	if(Rslt == eBSFFastaDescr)		// just read a descriptor line
		{
		m_pFasta->ReadDescriptor(szDescription,cBSFDescriptionSize);
		// chrom is identified as being all chars in descriptor line upto the first non-whitespace char
		pSrc = szDescription;
		pDst = szChrom;
		while((Chr = *pSrc++) != '\0' && !isspace(Chr))
			*pDst++ = Chr;
		*pDst = '\0';
		bDescriptor = true;
		continue;
		}
	else									// just read sequence
		if(!bDescriptor)					// if there was no descriptor then dummy up one...
			{
			sprintf(szDescription,"Probe%d",CurEntryID+1);
			bDescriptor = true;
			}

	// have a sequence
	// check if the descriptor matches any loci chromosome
	if((ChromID = m_pHypers->GetChromID(szChrom)) < 1)
		{
		NumNonMatchingPutativeChroms += 1;
		continue;
		}
	NumMatchingPutativeChroms += 1;

	// iterate the loci on this chromosome
	int Nth = 1;
	int ElStart;
	int ElLen;
	int MarkerSeqLen;
	int NumMarkersThisChrom = 0;
	int NumSNPBases;
	UINT8 SNPBases;
	int VariantID;
	NumMarkersThisChrom = 0;
	PrevMarkerEnd = 0;
	while((pHypEl = m_pHypers->LocateChromNthElement(ChromID,Nth)) != NULL)
		{
		Nth += 1;
		ElStart = pHypEl->StartLoci;
		ElLen = pHypEl->Len;
		if(ElStart < Extend5)
			{
			Num5Filtered += 1;
			continue;
			}
		if(ElStart+ElLen+Extend3 > CurSeqLen)
			{
			Num3Filtered += 1;
			continue;
			}
		pSeqStart = &m_pSeq[ElStart - Extend5];
		pSeqEnd = &m_pSeq[ElStart + ElLen + Extend3 - 1];
		MarkerSeqLen = 1 + (int)(pSeqEnd - pSeqStart);
		if(MarkerSeqLen < MinSeqLen)
			{
			NumUnderlen += 1;
			continue;
			}
		if(MarkerSeqLen > MaxSeqLen)
			{
			NumOverlen += 1;
			continue;
			}
		CSeqTrans::MapSeq2LCAscii(pSeqStart,MarkerSeqLen,szMarkerSeq);
		szMarkerSeq[MarkerSeqLen] = '\0';

		// check if marker contains too many indeterminate bases
		NumNs = 0;
		pBase = szMarkerSeq;
		while(Chr = *pBase++)
			{
			if(Chr == 'N')
				{
				NumNs += 1;
				if(NumNs > MaxNs)
					break;
				}
			}
		if(NumNs > MaxNs)
			{
			NumNfiltered += 1;
			continue;
			}

		// check if meeting non-overlap requirement WRT the previously accepted marker
		if(m_bNoOverlaps && PrevMarkerChromID == ChromID && PrevMarkerEnd > 0 && (ElStart - Extend5) <= PrevMarkerEnd)
			{
			NumOvrlapFiltered += 1;
			continue;
			}
		PrevMarkerChromID = ChromID; 
		PrevMarkerEnd = ElStart + ElLen + Extend3 - 1;

		NumSNPBases = pHypEl->NumSNPBases;
		SNPBases = pHypEl->SNPBases;
		VariantID = 0;
		do {
				// any replacements of the central loci base(s) required?
			if(NumSNPBases > 0)
				{
				szMarkerSeq[Extend5] = CSeqTrans::MapBase2Ascii((etSeqBase)(SNPBases & 0x03));
				SNPBases >>= 2;
				NumSNPBases -= 1;
				VariantID += 1;
				}
			else
				VariantID = 0;

			// write out the sequence here
			NumMarkers += 1;
			NumMarkersThisChrom += 1;
			if(NumMarkers > 1)
				MarkerBuffOfs += sprintf((char *)&m_pMarkerBuff[MarkerBuffOfs],"\n");
			if(VariantID == 0)
				MarkerBuffOfs += sprintf((char *)&m_pMarkerBuff[MarkerBuffOfs],">Marker%d %s|%d|%d|%d|%d\n%s",NumMarkers,szChrom,ElStart - Extend5,MarkerSeqLen,Extend5,Extend3,szMarkerSeq);
			else
				MarkerBuffOfs += sprintf((char *)&m_pMarkerBuff[MarkerBuffOfs],">Marker%d|%d %s|%d|%d|%d|%d\n%s",NumMarkers,VariantID,szChrom,ElStart - Extend5,MarkerSeqLen,Extend5,Extend3,szMarkerSeq);

			if(MarkerBuffOfs > (m_AllocdMarkerBuff - 5000))
				{
				CUtility::SafeWrite(m_hOutMarkerFile,m_pMarkerBuff,MarkerBuffOfs);
				MarkerBuffOfs = 0;
				}
			}
		while(NumSNPBases > 0);


		}
	if(NumMarkersThisChrom)
		NumAcceptedMarkerChroms += 1;
	}

if(m_hOutMarkerFile != -1)
	{
	if(MarkerBuffOfs > 0)
		CUtility::SafeWrite(m_hOutMarkerFile,m_pMarkerBuff,MarkerBuffOfs);
#ifdef _WIN32
	_commit(m_hOutMarkerFile);
#else
	fsync(m_hOutMarkerFile);
#endif
	close(m_hOutMarkerFile);
	m_hOutMarkerFile = -1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Number of markers generated: %d",NumMarkers);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Putative markers on %d chromosomes, no putative markers on %d chromosomes",NumMatchingPutativeChroms,NumNonMatchingPutativeChroms);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Accepted markers on %d chromosomes",NumAcceptedMarkerChroms);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Filtered out %d Ns, %d underlength, %d overlength, %d 5' flank, %d 3' flank, %d overlapping",NumNfiltered,NumUnderlen,NumOverlen,Num5Filtered,Num3Filtered,NumOvrlapFiltered);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Marker output multifasta file closed: '%s'",pszOutFile);
	}

Reset();
return(0);
}