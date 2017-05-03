/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

// MapLoci2Feat.cpp : contains the CMapLoci2Feat class implementation

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

#include "MapLoci2Feat.h"


CMapLoci2Feat::CMapLoci2Feat()
{
m_hRsltFile = -1;
m_hFeatRsltFile = -1;
m_pChromRegionCnts = NULL;
m_pBiobed = NULL;
m_pHypers = NULL;
m_pFeatCntDists = NULL;
MLFReset();
}


CMapLoci2Feat::~CMapLoci2Feat()
{
MLFReset();
}

void
CMapLoci2Feat::MLFReset(void)
{
	if (m_hRsltFile != -1)
	{
		close(m_hRsltFile);
		m_hRsltFile = -1;
	}
	if (m_hFeatRsltFile != -1)
	{
		close(m_hFeatRsltFile);
		m_hFeatRsltFile = -1;
	}
	if (m_pBiobed != NULL)
	{
		delete m_pBiobed;
		m_pBiobed = NULL;
	}
	if (m_pHypers != NULL)
	{
		delete m_pHypers;
		m_pHypers = NULL;
	}
	if (m_pFeatCntDists != NULL)
	{
		delete m_pFeatCntDists;
		m_pFeatCntDists = NULL;
	}

	if (m_pChromRegionCnts != NULL)
	{
		delete m_pChromRegionCnts;
		m_pChromRegionCnts = NULL;
	}

	m_MLFPMode = ePMdefault;
	m_StrandProc = eStrandDflt;
	m_IsoformRprt = eISOFPRPKM;
	m_RegRegionLen = cDfltRegLen;
	m_AllocdChromRegionCnts = 0;
	m_NumEls = 0;
	m_NumSplitEls = 0;
	m_bFeatinsts = false;
	m_bOneCntRead = false;
}

// MapLoci2Features
// Maps element loci to features and writes mapping for each element loci to pszRsltsFile
// Assumes that the bed file containing features (m_pBiobed) has been opened and that all element loci have been parsed into m_pHypers
int
CMapLoci2Feat::MapLoci2Features(char *pszRsltsFile)
{
	int Rslt;
	char szLineBuff[0x03fff];
	int BuffIdx;
	int SrcID;
	char Strand;
	char *pszChrom;
	int ChromID;
	char *pszElType;
	char *pszRefSpecies;
	char *pszRelSpecies;
	int PrevChromID;
	int PrevElTypeID;
	int PrevRefSpeciesID;
	int PrevRelSpeciesID;
	int StartLoci;
	int EndLoci;
	int Len;
	int Features;
	int AccumFeatures;
	int NumFeatsOverlap;
	int FeatMsk;
	int FeatIdx;
	UINT32 ElID;
	int FeatID;
	int NxtFeatID;
	int NxtFeatStart;
	int NxtFeatEnd;
	int PrvFeatID;
	int PrvFeatStart;
	int PrvFeatEnd;
	int RelScale;
	int NumMissingChroms;

	int PrevStartChromID = -1;
	int	PrevStartLoci = -1;
	char PrevStrand = '*';
	bool bStartUniqLoci = true;

	int MaxChromID;
	int CoreStartLoci;
	int CoreEndLoci;
	int TotNumFeatures;
	tsFeatCntDist *pCurFeatCntDist;		// to hold currently being processed feature count distribution
	tsHyperElement *pEl;
	int *pChromRegionCnts;

	if (m_hRsltFile != -1)				// ensure closed
	{
		close(m_hRsltFile);
		m_hRsltFile = -1;
	}
	if (pszRsltsFile != NULL && pszRsltsFile[0] != '\0')
	{
#ifdef _WIN32
		if ((m_hRsltFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE)) == -1)
#else
		if ((m_hRsltFile = open(pszRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE)) == -1)
#endif
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to create %s - %s", pszRsltsFile, strerror(errno));
			MLFReset();
			return(eBSFerrCreateFile);
		}
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Output file created/truncated: '%s'", pszRsltsFile);
	}

	pszChrom = NULL;
	PrevChromID = -1;
	pszElType = NULL;
	PrevElTypeID = -1;
	pszRefSpecies = NULL;
	PrevRefSpeciesID = -1;
	pszRelSpecies = NULL;
	PrevRelSpeciesID = -1;
	BuffIdx = 0;

	MaxChromID = 0;
	AccumFeatures = 0;


	m_AllocdChromRegionCnts = m_pBiobed->GetNumChromosomes();
	if(m_pChromRegionCnts != NULL)
		{
		delete m_pChromRegionCnts;
		m_pChromRegionCnts = NULL;
		}
	m_AllocdChromRegionCnts += 1;      // for safety!
	if((m_pChromRegionCnts = new tsChromRegionCnts[m_AllocdChromRegionCnts]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to allocate for %d tsChromRegionCnts instances", m_AllocdChromRegionCnts);
		MLFReset();
		return(eBSFerrMem);
		}

	memset(m_pChromRegionCnts, 0, sizeof(tsChromRegionCnts) * m_AllocdChromRegionCnts);

	// determine how many features
	TotNumFeatures = m_pBiobed->GetNumFeatures();
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Feature file contains %1.9d features", TotNumFeatures);
	if (m_pFeatCntDists == NULL)
	{
		if ((m_pFeatCntDists = new tsFeatCntDist[TotNumFeatures]) == NULL)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to allocate for %d tsFeatCntDist instances", TotNumFeatures);
			MLFReset();
			return(eBSFerrMem);
		}
		memset(m_pFeatCntDists, 0, sizeof(tsFeatCntDist) * TotNumFeatures);
		pCurFeatCntDist = m_pFeatCntDists;
		for (FeatID = 1; FeatID <= TotNumFeatures; FeatID++, pCurFeatCntDist++)
		{
			pCurFeatCntDist->FeatID = FeatID;
			pCurFeatCntDist->TranscribedLen = m_pBiobed->GetTranscribedLen(FeatID);
			pCurFeatCntDist->GeneLen = m_pBiobed->GetFeatLen(FeatID);
			m_pBiobed->GetFeature(FeatID, pCurFeatCntDist->szName);
			if (pCurFeatCntDist->szName[0] == '\0' || pCurFeatCntDist->GeneLen < 1)
			{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Inconsistency in feature: %d", FeatID);
				MLFReset();
				return(eBSFerrInternal);
			}


		}
	}

	NumMissingChroms = 0;
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Associating (from %d) element %9.9d", m_NumEls, 1);
	// Note: alignments will have been sorted ascending by chrom, start loci
	for (ElID = 1; ElID <= m_NumEls; ElID++)
	{
		AccumFeatures = 0;
		if (!(ElID % 100000))
			printf("\b\b\b\b\b\b\b\b\b%9.9d", ElID);
		pEl = m_pHypers->GetElement(ElID);
		if (pEl == NULL)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to get details for element: %d", ElID);
			MLFReset();
			return(eBSFerrInternal);
		}

		// if only associating the actual start loci, and if a split element then need to process first if element on '+' and last if element on '-' strand
		if (m_MLFPMode == ePMstarts && pEl->SplitElement == 1)
		{
			if (pEl->PlusStrand == 1 && pEl->SplitFirst != 1)
				continue;
			if (pEl->PlusStrand == 0 && pEl->SplitLast != 1)
				continue;
		}

		if (pszChrom == NULL || PrevChromID != pEl->ChromID)
		{
			pszChrom = m_pHypers->GetChrom(pEl->ChromID);
			if (pszChrom == NULL || pszChrom[0] == '\0')
			{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to get chrom text for element: %d ChromID: %d", ElID, pEl->ChromID);
				MLFReset();
				return(eBSFerrInternal);
			}

			// some old datasets may be referencing ChrM as mitochondria, or ChrC as chloroplast
			// so need to check for these
			if (!stricmp(pszChrom, "chloroplast"))
				pszChrom = (char *)"ChrC";
			else
				if (!stricmp(pszChrom, "mitochondria"))
					pszChrom = (char *)"ChrM";
			PrevChromID = pEl->ChromID;
			PrevStartChromID = -1;
			PrevStartLoci = -1;
			PrevStrand = '*';
			bStartUniqLoci = true;
		}

		if (pszElType == NULL || PrevElTypeID != pEl->ElTypeID)
		{
			pszElType = m_pHypers->GetType(pEl->ElTypeID);
			PrevElTypeID = pEl->ElTypeID;
		}

		if (pszRefSpecies == NULL || PrevRefSpeciesID != pEl->RefSpeciesID)
		{
			pszRefSpecies = m_pHypers->GetRefSpecies(pEl->RefSpeciesID);
			PrevRefSpeciesID = pEl->RefSpeciesID;
		}

		if (pszRelSpecies == NULL || PrevRelSpeciesID != pEl->RelSpeciesID)
		{
			pszRelSpecies = m_pHypers->GetRelSpecies(pEl->RelSpeciesID);
			PrevRelSpeciesID = pEl->RelSpeciesID;
		}


		StartLoci = pEl->StartLoci;
		EndLoci = pEl->StartLoci + pEl->Len - 1;

		Len = pEl->Len;
		Strand = pEl->PlusStrand ? '+' : '-';

		if (PrevStartChromID != pEl->ChromID || StartLoci != PrevStartLoci || Strand != PrevStrand) // not processed this loci previously?
		{
			PrevStartChromID = pEl->ChromID;
			PrevStartLoci = StartLoci;
			PrevStrand = Strand;
			bStartUniqLoci = true;
		}
		else
			bStartUniqLoci = false;

		SrcID = pEl->SrcID;
		Features = pEl->Features;
		RelScale = pEl->RelScale;
		switch (m_StrandProc)
		{
			case eStrandDflt:
				break;
			case eStrandSense:
				m_pBiobed->SetStrand(Strand);
				break;
			case eStrandAnti:
				m_pBiobed->SetStrand(Strand == '+' ? '-' : '+');
				break;
		}
		FeatID = 0;

		if (m_pBiobed != NULL)
		{
			Rslt = ChromID = m_pBiobed->LocateChromIDbyName(pszChrom);
			if (Rslt == eBSFerrChrom)
			{
				// some old datasets may be referencing ChrM as mitochondria, or ChrC as chloroplast
				// so need to check for these
				if (!stricmp(pszChrom, "ChrM"))
					Rslt = ChromID = m_pBiobed->LocateChromIDbyName((char *)"mitochondria");
				else
					if (!stricmp(pszChrom, "ChrC"))
						Rslt = ChromID = m_pBiobed->LocateChromIDbyName((char *)"chloroplast");
			}

			if (Rslt == eBSFerrChrom)
			{
				if (NumMissingChroms++ < 10)
					gDiagnostics.DiagOut(eDLInfo, gszProcName, "Unable to locate chromosome %s in BED file", pszChrom);
				if (NumMissingChroms == 10)
					gDiagnostics.DiagOut(eDLInfo, gszProcName, "Not reporting any additional missing chromosomes");
				Rslt = 0;
				continue;
			}

			if (Rslt > 0)
			{
				if (MaxChromID < ChromID)
					MaxChromID = ChromID;

				CoreStartLoci = StartLoci;
				CoreEndLoci = EndLoci;

				switch (m_MLFPMode)
				{
					case ePMdefault:
						break;

					case ePMstarts:
						if (Strand == '+')
							CoreEndLoci = CoreStartLoci;
						else
							CoreStartLoci = CoreEndLoci;
						break;

					case ePMdyad:
						if (Strand == '+')
						{
							CoreStartLoci += 73;
							CoreEndLoci = CoreStartLoci;
						}
						else
						{
							CoreEndLoci -= 73;
							if (CoreEndLoci < 0)
								CoreEndLoci = 0;
							CoreStartLoci = CoreEndLoci;
						}
						break;
				}

				// see if overlapping any features
				AccumFeatures = 0;
				NumFeatsOverlap = 0;
				do
				{
					FeatID = m_pBiobed->LocateFeatureIDinRangeOnChrom(ChromID,	// feature is on which chromosome
																	  CoreStartLoci,						// feature must end on or after Start
																	  CoreEndLoci,						// and start on or before End 
																	  NumFeatsOverlap + 1);				// Ith instance to return (1..n)

					if (FeatID > 0 && m_pFeatCntDists != NULL)
					{
						if (m_bFeatinsts)
						{
							Features = m_pBiobed->GetFeatureOverlaps(cRegionFeatBits, FeatID, CoreStartLoci, CoreEndLoci, m_RegRegionLen);
							Features |= m_pBiobed->GetFeatureBitsSpliceOverlaps(FeatID, CoreStartLoci, CoreEndLoci, cMinSpliceOverlap);
						}
						else
						{
							Features = m_pBiobed->GetFeatureBits(ChromID,			// feature is on which chromosome
																 CoreStartLoci,							// feature must end on or after Start
																 CoreEndLoci,							// and start on or before End
																 cRegionFeatBits,
																 m_RegRegionLen);

							Features |= m_pBiobed->GetSpliceSiteBits(ChromID, CoreStartLoci, CoreEndLoci, cMinSpliceOverlap);
						}

						if (m_bOneCntRead)
						{
							for (FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1)
								if (Features & FeatMsk)
								{
									Features = FeatMsk;
									break;
								}
						}
						AccumFeatures |= Features;
						pCurFeatCntDist = &m_pFeatCntDists[FeatID - 1];
						for (FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1)
							if (Features & FeatMsk)
								pCurFeatCntDist->RegionCnts[FeatIdx] += 1;

						// only accumulate relative abundance for reads in exons
						if (Features & (cFeatBitCDS | cFeatBit5UTR | cFeatBit3UTR))
						{
							pCurFeatCntDist->RelAbundance += 999.0 / (double)max(1, RelScale);
							if (bStartUniqLoci)
								pCurFeatCntDist->UniqueReadLociHits += 1;
						}
					}
					if (FeatID > 0)
						NumFeatsOverlap += 1;
				}
				while (FeatID > 0);

				if (!NumFeatsOverlap) // if not overlapping or not contained in any feature then locate nearest feature up/dnstream
				{
					// find feature starting after core end loci
					NxtFeatID = m_pBiobed->LocateFeatureAfter(ChromID,	// feature is on this chromosome
															  CoreEndLoci);					         // feature starts on or immediately after this offset
					if (NxtFeatID > 0)
						m_pBiobed->GetFeature(NxtFeatID,		// feature instance identifier
											  NULL,							// where to return feature name
											  NULL,							// where to return chromosome name
											  &NxtFeatStart,					// where to return feature start on chromosome (0..n) 
											  &NxtFeatEnd);					// where to return feature end on chromosome

																			// find feature ending before or at core start loci
					PrvFeatID = m_pBiobed->LocateFeatureBefore(ChromID,	// feature is on this chromosome
															   CoreStartLoci);			// feature ends on or immediately before this offset
					if (PrvFeatID > 0)
						m_pBiobed->GetFeature(PrvFeatID,		// feature instance identifier
											  NULL,	// where to return feature name
											  NULL,	// where to return chromosome name
											  &PrvFeatStart,		// where to return feature start on chromosome (0..n) 
											  &PrvFeatEnd);		// where to return feature end on chromosome

					if (NxtFeatID < 1)
						FeatID = PrvFeatID;
					else
					{
						if (PrvFeatID < 1)
							FeatID = NxtFeatID;
						else
						{
							if ((NxtFeatStart - CoreEndLoci) < (CoreStartLoci - PrvFeatEnd))
								FeatID = NxtFeatID;
							else
								FeatID = PrvFeatID;
						}
					}


					if (m_bFeatinsts)
						AccumFeatures = m_pBiobed->GetFeatureOverlaps(cRegionFeatBits, FeatID, CoreStartLoci, CoreEndLoci, m_RegRegionLen);
					else
						AccumFeatures = m_pBiobed->GetFeatureBits(ChromID,			// feature is on which chromosome
																  CoreStartLoci,							// feature must end on or after Start
																  CoreEndLoci,							// and start on or before End
																  cRegionFeatBits,
																  m_RegRegionLen);

					if (m_bOneCntRead)
					{
						for (FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1)
							if (AccumFeatures & FeatMsk)
							{
								AccumFeatures = FeatMsk;
								break;
							}
					}
					if (AccumFeatures && FeatID > 0 && m_pFeatCntDists != NULL)
					{
						pCurFeatCntDist = &m_pFeatCntDists[FeatID - 1];

						for (FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1)
							if (AccumFeatures & FeatMsk)
								pCurFeatCntDist->RegionCnts[FeatIdx] += 1;
					}
				}
			}
		}

		if (m_hRsltFile != -1)
			BuffIdx += sprintf(&szLineBuff[BuffIdx], "%d,\"%s\",\"%s\",\"%s\",%d,%d,%d,\"%c\",%d,%d\n",
							   SrcID, pszElType, pszRefSpecies, pszChrom, StartLoci, EndLoci, Len, Strand, AccumFeatures, RelScale);

		if (m_bOneCntRead)
		{
			for (FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1)
				if (AccumFeatures & FeatMsk)
				{
					AccumFeatures = FeatMsk;
					break;
				}
		}

		if (!AccumFeatures)
			m_pChromRegionCnts[ChromID - 1].Intergenic += 1;
		else
			{
			pChromRegionCnts = (int *)&m_pChromRegionCnts[ChromID - 1];
			for (FeatMsk = 0x01, FeatIdx = 0; FeatIdx < 8; FeatIdx++, FeatMsk <<= 1, pChromRegionCnts+=1)
				if (AccumFeatures & FeatMsk)
					*pChromRegionCnts += 1;
			}

		if (m_hRsltFile != -1 && ((BuffIdx + 1000) > sizeof(szLineBuff)))
		{
			CUtility::SafeWrite(m_hRsltFile, szLineBuff, BuffIdx);
			BuffIdx = 0;
		}
	}
	printf("\b\b\b\b\b\b\b\b\b%9.9d", ElID - 1);
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "All elements (%d) now associated", ElID - 1);
	if (NumMissingChroms > 0)
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "There were %d elements not associated because of missing chromosomes", NumMissingChroms);

	if (m_hRsltFile != -1)
	{
		if (BuffIdx > 0)
			CUtility::SafeWrite(m_hRsltFile, szLineBuff, BuffIdx);
#ifdef _WIN32
		_commit(m_hRsltFile);
#else
		fsync(m_hRsltFile);
#endif
		close(m_hRsltFile);
		m_hRsltFile = -1;
	}
	return(eBSFSuccess);
}

// CompareFeatName
// Used to sort feature names ascending
static int
CompareFeatName(const void *arg1, const void *arg2)
{
	tsFeatCntDist *pEl1 = (tsFeatCntDist *)arg1;
	tsFeatCntDist *pEl2 = (tsFeatCntDist *)arg2;
	return(stricmp(pEl1->szName, pEl2->szName));
}


//
// IsSameFeature
// Feature names must be at least cMinNameRootLen chars long
// To be an isoform they must have suffixes of ".[0-99]"
// feature names first have any suffix trimmed off and then the remaining root names are compared for equality
// 
bool
CMapLoci2Feat::IsSameFeature(char *pszFeatA, char *pszFeatB)
{
	char ChrA;
	char ChrB;
	char *pSfxA;
	char *pSfxB;
	int NameLen;
	bool bIsIsoform;
	if (pszFeatA == NULL || *pszFeatA == '\0')
		return(false);

	if (pszFeatB == NULL || *pszFeatB == '\0')
		return(false);

	NameLen = (int)strlen(pszFeatA);
	if (NameLen < cMinIsonameLen)
		return(false);
	pSfxA = &pszFeatA[NameLen - 1];
	if (*pSfxA >= '0' && *pSfxA <= '9')
	{
		pSfxA -= 1;
		if (*pSfxA >= '0' && *pSfxA <= '9')
			pSfxA -= 1;
		if (*pSfxA == '.')
		{
			*pSfxA = '\0';
			ChrA = '.';
		}
		else
		{
			ChrA = '\0';
			pSfxA = &pszFeatA[NameLen];
		}
	}
	else
	{
		pSfxA = &pszFeatA[NameLen];
		ChrA = '\0';
	}

	NameLen = (int)strlen(pszFeatB);
	if (NameLen < cMinIsonameLen)
		return(false);
	pSfxB = &pszFeatB[NameLen - 1];
	if (*pSfxB >= '0' && *pSfxB <= '9')
	{
		pSfxB -= 1;
		if (*pSfxB >= '0' && *pSfxB <= '9')
			pSfxB -= 1;
		if (*pSfxB == '.')
		{
			*pSfxB = '\0';
			ChrB = '.';
		}
		else
		{
			ChrB = '\0';
			pSfxB = &pszFeatB[NameLen];
		}
	}
	else
	{
		ChrB = '\0';
		pSfxB = &pszFeatB[NameLen];
	}

	bIsIsoform = stricmp(pszFeatA, pszFeatB) == 0 ? true : false;
	*pSfxA = ChrA;
	*pSfxB = ChrB;
	return(bIsIsoform);
}


// TrimNameIso
// Inplace remove any name isoform suffix of the form '.[0-99]'
bool			// true if suffix was trimmed
CMapLoci2Feat::TrimNameIso(char *pszName)
{
	char *pSfx;
	int NameLen;
	if (pszName == NULL || pszName[0] == '\0')
		return(false);
	NameLen = (int)strlen(pszName);
	if (NameLen < cMinIsonameLen)
		return(false);
	pSfx = &pszName[NameLen - 1];
	if (*pSfx >= '0' && *pSfx <= '9')
	{
		pSfx -= 1;
		if (*pSfx >= '0' && *pSfx <= '9')
			pSfx -= 1;
		if (*pSfx == '.')
		{
			*pSfx = '\0';
			return(true);
		}
	}
	return(false);
}

// IsIsoform
// assumes that if the feature name is at least cMinNameRootLen long and suffixed by '.[0-99]' then thats an isoform
bool
CMapLoci2Feat::IsIsoform(char *pszName)
{
	int NameLen;
	char *pSfx;
	NameLen = (int)strlen(pszName);
	if (NameLen < cMinIsonameLen)
		return(false);
	pSfx = &pszName[NameLen - 1];
	if (*pSfx >= '0' && *pSfx <= '9')
	{
		pSfx -= 1;
		if (*pSfx >= '0' && *pSfx <= '9')
			pSfx -= 1;
		if (*pSfx == '.')
			return(true);
	}
	return(false);

}

int
CMapLoci2Feat::MapFeatures2Loci(char *pszFeatRsltsFile)
{
	char szLineBuff[0x03fff];
	int BuffIdx;
	int TotNumFeatures;
	int FeatID;
	int FeatIdx;
	double RPKM;
	double LenRelAbundance;
	double SumTransLenRelAbundance;
	double UniqueHitsRelAbundance;
	double SumTransUniqueHitsRelAbundance;

	tsFeatCntDist *pCurFeatCntDist;		// to hold currently being processed feature count distribution
	tsFeatCntDist *pMaxRPKM;		// best RPKM isoform instance
	tsFeatCntDist *pMaxExonReads;		// best reads isoform instance
	bool bSameMaxRPKMFeature;		// true if processing same feature, different isoforms, for maximal RPKMs
	bool bSameMaxExonReadsFeature;  // true if processing same feature,  different isoforms, for maximal exonic reads

	if (m_hFeatRsltFile != -1)		// shouldn't be open but let's make sure...
	{
		close(m_hFeatRsltFile);
		m_hFeatRsltFile = -1;
	}
	if (m_pBiobed != NULL && pszFeatRsltsFile != NULL && pszFeatRsltsFile[0] != '\0')
	{
		// determine how many features
		TotNumFeatures = m_pBiobed->GetNumFeatures();
		if (m_pFeatCntDists == NULL)
		{
			if ((m_pFeatCntDists = new tsFeatCntDist[TotNumFeatures]) == NULL)
			{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to allocate for %d tsFeatCntDist instances", TotNumFeatures);
				MLFReset();
				return(eBSFerrMem);
			}
			memset(m_pFeatCntDists, 0, sizeof(tsFeatCntDist) * TotNumFeatures);
			pCurFeatCntDist = m_pFeatCntDists;
			for (FeatID = 1; FeatID <= TotNumFeatures; FeatID++, pCurFeatCntDist++)
			{
				pCurFeatCntDist->FeatID = FeatID;
				pCurFeatCntDist->TranscribedLen = m_pBiobed->GetTranscribedLen(FeatID);
				pCurFeatCntDist->GeneLen = m_pBiobed->GetFeatLen(FeatID);
				m_pBiobed->GetFeature(FeatID, pCurFeatCntDist->szName);
			}
		}

#ifdef _WIN32
		if ((m_hFeatRsltFile = open(pszFeatRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE)) == -1)
#else
		if ((m_hFeatRsltFile = open(pszFeatRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE)) == -1)
#endif
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to create %s - %s", pszFeatRsltsFile, strerror(errno));
			MLFReset();
			return(eBSFerrCreateFile);
		}
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Features output file created/truncated: '%s'", pszFeatRsltsFile);
	}
	else
		TotNumFeatures = 0;



	if (m_hFeatRsltFile != -1 && m_pFeatCntDists != NULL)
	{
		// sort by feature name 
		qsort(m_pFeatCntDists, TotNumFeatures, sizeof(tsFeatCntDist), CompareFeatName);

		// generate RPKMs and ExonReads over all features and flag if feature is an isoform
		SumTransLenRelAbundance = 0.0;
		SumTransUniqueHitsRelAbundance = 0.0;
		pCurFeatCntDist = m_pFeatCntDists;
		pMaxRPKM = NULL;
		pMaxExonReads = NULL;
		for (FeatID = 1; FeatID <= TotNumFeatures; FeatID++, pCurFeatCntDist++)
		{
			pCurFeatCntDist = &m_pFeatCntDists[FeatID - 1];
			pCurFeatCntDist->NumExonReads = pCurFeatCntDist->RegionCnts[0] + pCurFeatCntDist->RegionCnts[1] + pCurFeatCntDist->RegionCnts[2];
			if (pCurFeatCntDist->NumExonReads)	// if at least 1 read mapping to feature...
			{
				RPKM = ((double)pCurFeatCntDist->NumExonReads * 1000.0f);
				RPKM /= (double)pCurFeatCntDist->TranscribedLen;
				RPKM *= 1000000.0f / (double)m_NumEls;
				LenRelAbundance = pCurFeatCntDist->RelAbundance / (double)pCurFeatCntDist->TranscribedLen;
				UniqueHitsRelAbundance = pCurFeatCntDist->RelAbundance / (double)pCurFeatCntDist->UniqueReadLociHits;
			}
			else
			{
				RPKM = 0.0f;
				LenRelAbundance = 0;
				UniqueHitsRelAbundance = 0;
			}
			pCurFeatCntDist->RPKM = RPKM;
			pCurFeatCntDist->LenNormRelAbundance = LenRelAbundance;
			pCurFeatCntDist->UniqueHitsRelAbundance = UniqueHitsRelAbundance;
			SumTransLenRelAbundance += LenRelAbundance;
			SumTransUniqueHitsRelAbundance += UniqueHitsRelAbundance;
			pCurFeatCntDist->bIsIsoform = IsIsoform(pCurFeatCntDist->szName);
			if (m_IsoformRprt == eISOFPall)		// if reporting all isoforms then that's easy...
			{
				pCurFeatCntDist->bMaxRPKM = true;
				pCurFeatCntDist->bMaxExonReads = true;
				continue;
			}


			bSameMaxRPKMFeature = pMaxRPKM == NULL ? false : IsSameFeature(pCurFeatCntDist->szName, pMaxRPKM->szName);
			bSameMaxExonReadsFeature = pMaxExonReads == NULL ? false : IsSameFeature(pCurFeatCntDist->szName, pMaxExonReads->szName);

			// will report on 'maximal' isoforms
			// if first iteration (pMaxRPKM will be NULL) or not an isoform of current maximal RPKM ...
			if (pMaxRPKM == NULL || !bSameMaxRPKMFeature)
			{
				pCurFeatCntDist->bMaxRPKM = true;
				pMaxRPKM = pCurFeatCntDist;
			}

			// if first iteration (pMaxExonReads will be NULL) or not an isoform of current maximal MaxExonReads ...
			if (pMaxExonReads == NULL || !bSameMaxExonReadsFeature)
			{
				pCurFeatCntDist->bMaxExonReads = true;
				pMaxExonReads = pCurFeatCntDist;
			}

			// same feature - but different isoforms
			if (bSameMaxRPKMFeature && pCurFeatCntDist->RPKM > pMaxRPKM->RPKM)
			{
				pMaxRPKM->bMaxRPKM = false;
				pCurFeatCntDist->bMaxRPKM = true;
				pMaxRPKM = pCurFeatCntDist;
			}

			if (bSameMaxExonReadsFeature && pCurFeatCntDist->NumExonReads > pMaxExonReads->NumExonReads)
			{
				pMaxExonReads->bMaxExonReads = false;
				pCurFeatCntDist->bMaxExonReads = true;
				pMaxExonReads = pCurFeatCntDist;
			}
		}

		// now total is known, can determine abundance proportions for each transcript 
		if (SumTransLenRelAbundance > 0.0)
		{
			for (FeatID = 1; FeatID <= TotNumFeatures; FeatID++, pCurFeatCntDist++)
			{
				pCurFeatCntDist = &m_pFeatCntDists[FeatID - 1];
				pCurFeatCntDist->TransRelAbundance = pCurFeatCntDist->LenNormRelAbundance / SumTransLenRelAbundance;
				pCurFeatCntDist->TransUniqueHitsRelAbundance = pCurFeatCntDist->UniqueHitsRelAbundance / SumTransUniqueHitsRelAbundance;
			}
		}

		BuffIdx = sprintf(szLineBuff, "\"FeatID\",\"Feature\",\"GeneLen\",\"TransLen\",\"CDS\",\"5'UTR\",\"3'UTR\",\"Introns\",\"5'upstream\",\"3'downstream\",\"intron3'/5'exon\",\"exon3'/5'intron\",\"RPKM\",\"ExonReads\",\"UniqueLociHits\",\"RelAbundance\",\"UniqueRelAbundance\",\"TransUniqueRelAbundance\",\"LenRelAbundance\",\"TransLenRelAbundance\"");
		pCurFeatCntDist = m_pFeatCntDists;

		for (FeatID = 1; FeatID <= TotNumFeatures; FeatID++, pCurFeatCntDist++)
		{
			switch (m_IsoformRprt)
			{
				case eISOFPRPKM:
					if (!pCurFeatCntDist->bMaxRPKM)
						continue;
					break;
				case eISOFReads:
					if (!pCurFeatCntDist->bMaxExonReads)
						continue;
					break;
				default:
					break;
			}
			BuffIdx += sprintf(&szLineBuff[BuffIdx], "\n%d,\"%s\",%d,%d",
							   FeatID, pCurFeatCntDist->szName, pCurFeatCntDist->GeneLen, pCurFeatCntDist->TranscribedLen);
			pCurFeatCntDist = &m_pFeatCntDists[FeatID - 1];
			for (FeatIdx = 0; FeatIdx < 8; FeatIdx++)
				BuffIdx += sprintf(&szLineBuff[BuffIdx], ",%d", pCurFeatCntDist->RegionCnts[FeatIdx]);
			BuffIdx += sprintf(&szLineBuff[BuffIdx], ",%1.6f,%d,%d,%1.9f,%1.9f,%1.9f,%1.9f,%1.9f",
							   pCurFeatCntDist->RPKM, pCurFeatCntDist->NumExonReads,
							   pCurFeatCntDist->UniqueReadLociHits,
							   pCurFeatCntDist->RelAbundance,
							   pCurFeatCntDist->UniqueHitsRelAbundance,
							   pCurFeatCntDist->TransUniqueHitsRelAbundance,
							   pCurFeatCntDist->LenNormRelAbundance,
							   pCurFeatCntDist->TransRelAbundance);

			// something really strange on Ubuntu - sometimes a number of features are replicated at the end of the generated file
			// almost as though the final write was occuring 2x but I can't determine why so am playing safe now and manually writing out the
			// buffer on the final iteration and also doing a fsync immediately before closing the handle
			if (FeatID == TotNumFeatures || ((BuffIdx + 1000) > sizeof(szLineBuff)))
			{
				if (!CUtility::SafeWrite(m_hFeatRsltFile, szLineBuff, BuffIdx))
				{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to write %s", pszFeatRsltsFile);
					MLFReset();
					return(eBSFerrCreateFile);
				}
				BuffIdx = 0;
			}
		}
		if (BuffIdx)
			CUtility::SafeWrite(m_hFeatRsltFile, szLineBuff, BuffIdx);
#ifdef _WIN32
		_commit(m_hFeatRsltFile);
#else
		fsync(m_hFeatRsltFile);
#endif
		BuffIdx = 0;
		close(m_hFeatRsltFile);
		m_hFeatRsltFile = -1;
	}
	if (m_pFeatCntDists != NULL)
	{
		delete m_pFeatCntDists;
		m_pFeatCntDists = NULL;
	}
	return(eBSFSuccess);
}


int
CMapLoci2Feat::MLFProcess(etMLFPMode PMode,				// processing mode
						  bool bDedupe,				// true if input elements are to be deduped
						  bool bFeatinsts,			// true if input elements are to be associated to individual features, false if to all features at that locus
						  bool bOneCntRead,			// true if one count per read rule to be applied (functional regions are prioritised with CDS as the highest) 
						  etISOFProc IsoformRprt,		// feature isoform reporting mode
						  etStrandProc StrandProc,	// how to process read + element strand
						  int FType,					// input element file format: 0 - auto, 1 - CSV, 2 - BED, 3 - SAM
						  teCSVFormat CSVFormat,		// if CSV input then expected file format
						  char *pszInLociFile,		// input CSV, BED or SAM loci file
						  char *pszInBEDFile,			// input BED file
						  char *pszRsltsFile,			// output loci to feature mapping file
						  char *pszFeatRsltsFile,		// optional feature mapping results file
						  char *pszSummRsltsFile,		// optional output chrom summary results file
						  int RegRegionLen,			// regulatory region length
						  int MinLength,				// minimum element length
						  int MaxLength,				// maximum element length
						  int JoinOverlap)			// deduping join overlap
{
	int Rslt;
	char *pszChrom;
	int ChromID;
	char szChrom[cMaxDatasetSpeciesChrom];
	int TmpNumEls;

	MLFReset();
	m_MLFPMode = PMode;
	m_StrandProc = StrandProc;
	m_IsoformRprt = IsoformRprt;
	m_RegRegionLen = RegRegionLen;
	m_bFeatinsts = bFeatinsts;
	m_bOneCntRead = bOneCntRead;

	if ((m_pHypers = new CHyperEls) == NULL)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CHyperEls");
		return(eBSFerrObj);
	}

	etClassifyFileType FileType;

	if (pszInBEDFile != NULL && pszInBEDFile[0] != '\0')
	{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Loading features file '%s'", pszInBEDFile);
		if ((m_pBiobed = new CBEDfile) == NULL)
		{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CBEDfile gene/exon file '%s'", pszInBEDFile);
			MLFReset();
			return(eBSFerrObj);
		}

		if ((Rslt = m_pBiobed->Open(pszInBEDFile, eBTGeneExons)) != eBSFSuccess)
		{
			while (m_pBiobed->NumErrMsgs())
				gDiagnostics.DiagOut(eDLFatal, gszProcName, m_pBiobed->GetErrMsg());
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open features file '%s'", pszInBEDFile);
			MLFReset();
			return(eBSFerrObj);
		}
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Opened features file '%s'", pszInBEDFile);
	}
	else
	{
		m_pBiobed = NULL;
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No gene/exon BED file specified");
		MLFReset();
		return(eBSFerrObj);
	}

	if (FType == 0)
		FileType = CUtility::ClassifyFileType(pszInLociFile);
	else
		FileType = (etClassifyFileType)(FType - 1);

	switch (FileType)
	{
		case eCFTopenerr:		// unable to open file for reading
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open file: '%s'", pszInLociFile);
			return(eBSFerrOpnFile);

		case eCFTlenerr:		// file length is insufficent to classify type
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Unable to classify file type (insufficent data points): '%s'", pszInLociFile);
			return(eBSFerrFileAccess);

		case eCFTunknown:		// unable to reliably classify
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Unable to reliably classify file type: '%s'", pszInLociFile);
			return(eBSFerrFileType);

		case eCFTCSV:			// file has been classified as being CSV
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Parsing CSV file (%d..%d): '%s'", MinLength, MaxLength, pszInLociFile);
			if ((Rslt = m_pHypers->ParseCSVFileElements(pszInLociFile, MinLength, MaxLength, CSVFormat)) < 0)
			{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Parse errors in CSV file (%d..%d): '%s'", MinLength, MaxLength, pszInLociFile);
				MLFReset();
				return(Rslt);
			}
			break;

		case eCFTBED:			// file has been classified as being BED
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Parsing BED file (%d..%d): '%s'", MinLength, MaxLength, pszInLociFile);
			if ((Rslt = m_pHypers->ParseBEDFileElements(pszInLociFile, MinLength, MaxLength)) < 0)
			{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Parse errors in BED file (%d..%d): '%s'", MinLength, MaxLength, pszInLociFile);
				MLFReset();
				return(Rslt);
			}
			break;

		case eCFTSAM:			// file has been classified as being SAM
			gDiagnostics.DiagOut(eDLInfo, gszProcName, "Parsing SAM file (%d..%d): '%s'", MinLength, MaxLength, pszInLociFile);
			if ((Rslt = m_pHypers->ParseSAMFileElements(pszInLociFile, MinLength, MaxLength)) < 0)
			{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Parse errors in SAM file (%d..%d): '%s'", MinLength, MaxLength, pszInLociFile);
				MLFReset();
				return(Rslt);
			}
			break;
	}


	m_NumEls = m_pHypers->NumEls();
	if (m_NumEls == 0)
	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "No elements with length range %d..%d in file: '%s'", MinLength, MaxLength, pszInLociFile);
		MLFReset();
		return(Rslt);
	}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Loaded and parsed %d elements", m_NumEls);

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Now identifying split (microInDels or splice junction spanning?) elements...", m_NumEls);
	m_NumSplitEls = m_pHypers->IdentifySplitElements();					// identify any split elements which may be present
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Identified %d split (microInDels or splice junction spanning?) elements", m_NumSplitEls);

	if (JoinOverlap > 0 || bDedupe)
	{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Now joining/deduping elements...", m_NumEls);
		TmpNumEls = m_NumEls;
		m_NumEls = m_pHypers->DedupeSort(JoinOverlap, bDedupe);				// dedupe, join and sort elements ascending
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "After join/dedupe, %d elements removed, %d elements remaining",
							 TmpNumEls - m_NumEls, m_NumEls);
	}

	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Mapping elements to features...");
	if ((Rslt = MapLoci2Features(pszRsltsFile)) < eBSFSuccess)
	{
		MLFReset();
		return(Rslt);
	}
	gDiagnostics.DiagOut(eDLInfo, gszProcName, "Mapping elements to features completed");

	if (pszFeatRsltsFile != NULL && pszFeatRsltsFile[0] != '\0')
	{
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Mapping features to elements...");
		if ((Rslt = MapFeatures2Loci(pszFeatRsltsFile)) < eBSFSuccess)
		{
			MLFReset();
			return(Rslt);
		}
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Mapping features to elements completed");
	}

	if (pszSummRsltsFile != NULL && pszSummRsltsFile[0] != '\0')
	{
		char szOutBuff[8000];
		int BufIdx;
		int hOutFile;
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Generating summary chromosome results file...");

#ifdef _WIN32
		if ((hOutFile = open(pszSummRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE)) == -1)
#else
		if ((hOutFile = open(pszSummRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE)) == -1)
#endif
		{
			gDiagnostics.DiagOutMsgOnly(eDLFatal, "Unable to create or truncate  output file %s error: %s", pszSummRsltsFile, strerror(errno));
			return(-1);
		}
		BufIdx = sprintf(szOutBuff, "Chrom, CDS, 5'UTR, 3'UTR, Introns, 5'upstream, 3'downstream, intron3'/5'exon, exon3'/5'intron, Intergenic\n");
		for (ChromID = 0; ChromID < m_pBiobed->GetNumChromosomes(); ChromID++)
		{
			m_pBiobed->GetChromosome(ChromID + 1, szChrom);
			pszChrom = szChrom;
			if (!stricmp(pszChrom, "chloroplast"))
				pszChrom = (char *)"ChrC";
			else
				if (!stricmp(pszChrom, "mitochondria"))
					pszChrom = (char *)"ChrM";
			BufIdx += sprintf(&szOutBuff[BufIdx], "\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
							  pszChrom, m_pChromRegionCnts[ChromID].CDS, m_pChromRegionCnts[ChromID].UTR5, m_pChromRegionCnts[ChromID].UTR3,
							 m_pChromRegionCnts[ChromID].Intronic, m_pChromRegionCnts[ChromID].upstream5, m_pChromRegionCnts[ChromID].downstream3,
							 m_pChromRegionCnts[ChromID].intron3exon, m_pChromRegionCnts[ChromID].exon3intron, m_pChromRegionCnts[ChromID].Intergenic);
			if ((BufIdx + 1000) > sizeof(szOutBuff))
			{
				CUtility::SafeWrite(hOutFile, szOutBuff, BufIdx);
				BufIdx = 0;
			}
		}
		if (BufIdx)
			CUtility::SafeWrite(hOutFile, szOutBuff, BufIdx);
#ifdef _WIN32
		_commit(hOutFile);
#else
		fsync(hOutFile);
#endif
		close(hOutFile);
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Generating summary chromosome results file completed");
	}

	MLFReset();
	return(Rslt < 0 ? m_NumEls : Rslt);
}


