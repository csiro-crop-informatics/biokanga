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

#include "./biokanga.h"
#include "SNPLoci.h"


CSNPLoci::CSNPLoci(void)
{
m_pCSV = NULL;
Reset();
}


CSNPLoci::~CSNPLoci(void)
{
}

void
CSNPLoci::Reset(void)
{
if(m_pCSV != NULL)
	{
	delete m_pCSV;
	m_pCSV = NULL;
	}
}


// Input CSV format expected for SNP calls as generated by 'biokanga align'
// "SNP_ID","ElType","Species","Chrom","StartLoci","EndLoci","Len","Strand","Rank","PValue","Bases","Mismatches","RefBase",
// "MMBaseA","MMBaseC","MMBaseG","MMBaseT","MMBaseN","BackgroundSubRate","TotWinBases","TotWinMismatches"
//
// Input CSV format expected for marker calls as generated by 'biokanga snpmarkers'
// "UniGene57K:TargSeq","Loci","TargBase","NumSpeciesWithCnts",
//     "Baxter:Base","Baxter:Score","Baxter:BaseCntTot","Baxter:BaseCntA","Baxter:BaseCntC","Baxter:BaseCntG","Baxter:BaseCntT","Baxter:BaseCntN",
//     "Baxter:SrcCnts","Chara:Base","Chara:Score","Chara:BaseCntTot","Chara:BaseCntA","Chara:BaseCntC","Chara:BaseCntG","Chara:BaseCntT","Chara:BaseCntN",
//     ... repeated for as many cultivars as processed ....
//     "Yitpi:SrcCnts","Yitpi:Base","Yitpi:Score","Yitpi:BaseCntTot","Yitpi:BaseCntA","Yitpi:BaseCntC","Yitpi:BaseCntG","Yitpi:BaseCntT","Yitpi:BaseCntN"

typedef struct TAG_sBaseCnts {
	etSeqBase Base;     // counts are for this base
	double Proportion;	// this base counts as proportion of all counts
	int Cnts;			// number of counts for this base
	int AllelicRank;	// this bases allelic ranking
} tsBaseCnts;

int 
CSNPLoci::LoadSNPs(char *pszSNPFile)				// load SNP calls from this CSV file
{
int Rslt;
int NumFields;
int CultivarIdx;
UINT32 EstNumSNPs;
UINT32 NumSNPsParsed;
UINT32 RowNumber;
bool bMarkerFormat;
int ExpNumFields;
int BaseCntsIdx;
tsBaseCnts *pBaseCnts1;
tsBaseCnts *pBaseCnts2;
tsBaseCnts *pMajorAllele;
tsBaseCnts *pMinorAllele;
tsBaseCnts BaseCnts[5];
char szChrom[cMaxDatasetSpeciesChrom+1];
eSeqBase  RefBase;
int SNPLoci;
int CntBases;
int CntRef;
int CntMM;
int NumSpeciesWithCnts;

if(m_pCSV != NULL)
	delete m_pCSV;
if((m_pCSV = new CCSVFile)== NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}

if((Rslt=m_pCSV->Open(pszSNPFile))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszSNPFile);
	Reset();
	return(Rslt);
	}

EstNumSNPs = m_pCSV->EstNumRows();
NumSNPsParsed = 0;
RowNumber = 0;
ExpNumFields = 0;
m_NumCultivars = 0;
bMarkerFormat = false;
while((Rslt=m_pCSV->NextLine()) > 0)				// onto next line containing fields
	{
	RowNumber += 1;
	NumFields = m_pCSV->GetCurFields();	
	if(ExpNumFields && ExpNumFields != NumFields)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Inconsistency in number of fields, previously %d but %d fields parsed from '%s' near line %d",ExpNumFields,NumFields,pszSNPFile,RowNumber);
		Reset();
		return(eBSFerrFieldCnt);
		}
	if(!ExpNumFields)
		{
		if(NumFields < cSNPMarkerNumFields)								// must be at least this many if 'biokanga snpmarkers' format with single cultivar
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Too few (%d) fields parsed from '%s' near line %d, is this file generated by 'biokanga align/snpmarkers'?",NumFields,pszSNPFile,RowNumber);
			Reset();
			return(eBSFerrFieldCnt);
			}
		if(NumFields == cAlignNumSNPfields)								// assume 'biokanga align'?
			{
			m_NumCultivars = 1;
			bMarkerFormat = false;
			}
		else
			{
			m_NumCultivars = (NumFields - 4)/9;
			if(((m_NumCultivars * 9) + 4) != NumFields)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Invalid number (%d) of fields parsed from '%s' near line %d, is this file generated by 'biokanga align/snpmarkers'?",NumFields,pszSNPFile,RowNumber);
				Reset();
				return(eBSFerrFieldCnt);
				}
			if(m_NumCultivars > cMaxCultivars)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Too many cultivars (%d) in '%s', max allowed is %d",m_NumCultivars,pszSNPFile,cMaxCultivars);
				Reset();
				return(eBSFerrFieldCnt);
				}
			bMarkerFormat = true;
			}
		ExpNumFields = NumFields;
		}

	if(RowNumber == 1)
		{
		if(m_pCSV->IsLikelyHeaderLine()) 
			{
			if(!NumSNPsParsed && bMarkerFormat)
				{
				// parse out target assembly name against which alignments were made and the cultivar names
				char *pSrc;
				char *pDst;
				char Chr;
				int Len;

				m_pCSV->GetText(1,&pSrc);
				pDst = m_TargAssemblyName;
				Len = 0;
				while(Len < sizeof(m_TargAssemblyName)-1 && (Chr = *pSrc++) && Chr != ':')
					{
					*pDst++ = Chr;
					Len++;
					*pDst='\0';
					}

				for(CultivarIdx = 0; CultivarIdx < m_NumCultivars; CultivarIdx++)
					{
					memset(&m_Cultivars[CultivarIdx],0,sizeof(m_Cultivars[CultivarIdx]));
					m_pCSV->GetText(4+(CultivarIdx*9),&pSrc);
					pDst = m_Cultivars[CultivarIdx].szName;
					Len = 0;
					while(Len < sizeof(m_Cultivars[CultivarIdx].szName)-1 && (Chr = *pSrc++) && Chr != ':')
						{
						*pDst++ = Chr;
						Len++;
						*pDst='\0';
						}
					}
				}
			continue;
			}
		else
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected CSV file '%s' first line to be a header line with  fully quoted field names",pszSNPFile);
			Reset();
			return(eBSFerrFieldCnt);
			}
		}

	if(!NumSNPsParsed && !bMarkerFormat)
		{
		// parse out target assembly name against which alignments were made and the cultivar names
		char *pSrc;
		char *pDst;
		char Chr;
		int Len;

		m_pCSV->GetText(3,&pSrc);
		pDst = m_TargAssemblyName;
		Len = 0;
		while(Len < sizeof(m_TargAssemblyName)-1 && (Chr = *pSrc++))
			{
			*pDst++ = Chr;
			Len++;
			*pDst='\0';
			}
		memset(&m_Cultivars[0],0,sizeof(m_Cultivars[0]));
		strcpy(m_Cultivars[0].szName,(char *)"Unknown");
		}

	char *pszTxt;


	memset(BaseCnts,0,sizeof(BaseCnts));
	pBaseCnts1 = &BaseCnts[0];
	for(BaseCntsIdx = 0; BaseCntsIdx < 5; BaseCntsIdx++,pBaseCnts1++)
		pBaseCnts1->Base = (etSeqBase)BaseCntsIdx;

	if(!bMarkerFormat)
		{
		int SNPlen;
		m_pCSV->GetInt(7,&SNPlen);			// check that we are processing SNPs!
		if(SNPlen != 1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected SNP CSV file '%s' to only contain 1 base SNPs, 'Len' = %d near line %d",pszSNPFile,SNPlen,RowNumber);
			Reset();
			return(eBSFerrFieldCnt);
			}
		m_pCSV->GetText(4,&pszTxt);			// get "Chrom"
		strncpy(szChrom,pszTxt,sizeof(szChrom));
		szChrom[sizeof(szChrom)-1] = '\0';
		m_pCSV->GetInt(5,&SNPLoci);			// get "StartLoci"
		m_pCSV->GetInt(11,&CntBases);		// get "Bases"
		m_pCSV->GetInt(12,&CntMM);			// get "Mismatches"
		CntRef = CntBases - CntMM;
		m_pCSV->GetText(13,&pszTxt);		// get "RefBase"
		m_pCSV->GetInt(14,&BaseCnts[0].Cnts);	// get "MMBaseA"
		m_pCSV->GetInt(15,&BaseCnts[1].Cnts);	// get "MMBaseC"
		m_pCSV->GetInt(16,&BaseCnts[2].Cnts);	// get "MMBaseG"
		m_pCSV->GetInt(17,&BaseCnts[3].Cnts);	// get "MMBaseT"
		m_pCSV->GetInt(18,&BaseCnts[4].Cnts);	// get "MMBaseN"
		switch(*pszTxt) {
			case 'a': case 'A':
				RefBase = eBaseA;
				BaseCnts[0].Cnts = CntRef;
				break;
			case 'c': case 'C':
				RefBase = eBaseC;
				BaseCnts[1].Cnts = CntRef;
				break;
			case 'g': case 'G':
				RefBase = eBaseG;
				BaseCnts[2].Cnts = CntRef;
				break;
			case 't': case 'T': case 'u': case 'U':	// U in case RNA alignments..
				RefBase = eBaseT;
				BaseCnts[3].Cnts = eBaseT;
				break;
			case 'n': case 'N':				// unlikely to have a SNP against an indeterminate base but you never know...
				RefBase = eBaseN;
				BaseCnts[4].Cnts = CntRef;
				break;
			default:
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected SNP RefBase ('%s') to be one of 'ACGTN' in CSV file '%s' near line %d",pszTxt,pszSNPFile,RowNumber);
				Reset();
				return(eBSFerrFieldCnt);
			}

		// Which is the major allele, and possibly also a minor allele
		// sort by count ascending
		tsBaseCnts TmpCnts;
		pBaseCnts1 = BaseCnts;
		pBaseCnts2 = pBaseCnts1+1;
		for(BaseCntsIdx = 0; BaseCntsIdx < 4; BaseCntsIdx++,pBaseCnts1++,pBaseCnts2++)
			if(pBaseCnts1->Cnts > pBaseCnts2->Cnts)
				{
				TmpCnts = *pBaseCnts1;
				*pBaseCnts1 = *pBaseCnts2;
				*pBaseCnts2 = TmpCnts;
				}
		pBaseCnts1 = BaseCnts;
		pBaseCnts2 = pBaseCnts1+1;
		for(BaseCntsIdx = 0; BaseCntsIdx < 3; BaseCntsIdx++,pBaseCnts1++,pBaseCnts2++)
			if(pBaseCnts1->Cnts > pBaseCnts2->Cnts)
				{
				TmpCnts = *pBaseCnts1;
				*pBaseCnts1 = *pBaseCnts2;
				*pBaseCnts2 = TmpCnts;
				}
		pMajorAllele = &BaseCnts[5];
		pMinorAllele = &BaseCnts[4];
		}
	else     // else parsing for snpmarker format
		{
		m_pCSV->GetText(1,&pszTxt);			// get ":TargSeq"
		strncpy(szChrom,pszTxt,sizeof(szChrom));
		szChrom[sizeof(szChrom)-1] = '\0';
		m_pCSV->GetInt(2,&SNPLoci);			// get "Loci"
		m_pCSV->GetText(3,&pszTxt);			// get "TargBase"
		switch(*pszTxt) {
			case 'a': case 'A':
				RefBase = eBaseA;
				break;
			case 'c': case 'C':
				RefBase = eBaseC;
				break;
			case 'g': case 'G':
				RefBase = eBaseG;
				break;
			case 't': case 'T': case 'u': case 'U':	// U in case RNA alignments..
				RefBase = eBaseT;
				break;
			case 'n': case 'N':				// unlikely to have a SNP against an indeterminate base but you never know...
				RefBase = eBaseN;
				break;
			default:
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected SNP TargBase ('%s') to be one of 'ACGTN' in CSV file '%s' near line %d",pszTxt,pszSNPFile,RowNumber);
				Reset();
				return(eBSFerrFieldCnt);
			}
		m_pCSV->GetInt(4,&NumSpeciesWithCnts);			// get "NumSpeciesWithCnts"
		if(NumSpeciesWithCnts < 1 || NumSpeciesWithCnts > m_NumCultivars)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected NumSpeciesWithCnts (%d) to be between 1 and %d in CSV file '%s' near line %d",NumSpeciesWithCnts,m_NumCultivars,pszSNPFile,RowNumber);
			Reset();
			return(eBSFerrNumRange);
			}

		int FieldIdx = 5;
		for(CultivarIdx = 0; CultivarIdx < m_NumCultivars; CultivarIdx++)
			{
			m_pCSV->GetText(FieldIdx++,&pszTxt);		// get cultivar ":Base"
			switch(*pszTxt) {
				case 'a': case 'A':
					RefBase = eBaseA;
					break;
				case 'c': case 'C':
					RefBase = eBaseC;
					break;
				case 'g': case 'G':
					RefBase = eBaseG;
					break;
				case 't': case 'T': case 'u': case 'U':	// U in case RNA alignments..
					RefBase = eBaseT;
					break;
				case 'n': case 'N':				// unlikely to have a SNP against an indeterminate base but you never know...
					RefBase = eBaseN;
					break;
				default:
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected SNP TargBase ('%s') to be one of 'ACGTN' in CSV file '%s' near line %d",pszTxt,pszSNPFile,RowNumber);
					Reset();
					return(eBSFerrFieldCnt);
				}
			m_pCSV->GetInt(FieldIdx++,&CntBases);		// get ":Score"
			m_pCSV->GetInt(FieldIdx++,&CntBases);		// get ":BaseCntTot"
			m_pCSV->GetInt(FieldIdx++,&CntBases);		// get ":BaseCntA"
			m_pCSV->GetInt(FieldIdx++,&CntBases);		// get ":BaseCntC"
			m_pCSV->GetInt(FieldIdx++,&CntBases);		// get ":BaseCntG"
			m_pCSV->GetInt(FieldIdx++,&CntBases);		// get ":BaseCntT"
			m_pCSV->GetInt(FieldIdx++,&CntBases);		// get ":BaseCntN"
			}

		// Which is the major allele, and possibly also a minor allele
		// sort by count ascending
		tsBaseCnts TmpCnts;
		pBaseCnts1 = BaseCnts;
		pBaseCnts2 = pBaseCnts1+1;
		for(BaseCntsIdx = 0; BaseCntsIdx < 4; BaseCntsIdx++,pBaseCnts1++,pBaseCnts2++)
			if(pBaseCnts1->Cnts > pBaseCnts2->Cnts)
				{
				TmpCnts = *pBaseCnts1;
				*pBaseCnts1 = *pBaseCnts2;
				*pBaseCnts2 = TmpCnts;
				}
		pBaseCnts1 = BaseCnts;
		pBaseCnts2 = pBaseCnts1+1;
		for(BaseCntsIdx = 0; BaseCntsIdx < 3; BaseCntsIdx++,pBaseCnts1++,pBaseCnts2++)
			if(pBaseCnts1->Cnts > pBaseCnts2->Cnts)
				{
				TmpCnts = *pBaseCnts1;
				*pBaseCnts1 = *pBaseCnts2;
				*pBaseCnts2 = TmpCnts;
				}
		pMajorAllele = &BaseCnts[5];
		pMinorAllele = &BaseCnts[4];
		}

	}

return(0);
}

int 
CSNPLoci::LoadSeqs(char *pszSeqFile)				// load sequences from this multifasta file, SNPs were called relatative to these sequences
{
return(0);
}

int 
CSNPLoci::Filter(int MinSep)						// filter out SNPs which are not separated by at least this many bases from any other SNP loci 
{
return(0);
}

int 
CSNPLoci::Dedupe(bool bSenseOnly)			// remove any SNP sequences which are duplicates of other SNP sequences
{
return(0);
}

int 
CSNPLoci::Report(char *pszOutFile)				// report SNP sequences to this file
{
return(0);
}

int 
CSNPLoci::Process(char *pszSNPFile,					// load SNP calls from this CSV file
				char *pszSeqFile,			// load sequences from this multifasta file, SNPs were called relatative to these sequences
				char *pszOutFile,			// report SNP sequences to this file
				int  Extd5,					// extend SNP 5' this many bases
				int  Extd3,					// extend SNP 3' this many bases
				int MinSep,					// filter out SNPs which are not separated by at least this many bases from any other SNP loci
				bool bSenseOnly)			// remove any SNP sequences which are duplicates of other SNP sequences
{
return(0);
}