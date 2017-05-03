/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */
// Extracts fasta sequences from a multifasta file
// Sequences to be extracted are identified by their descriptors

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

#include "MergeReadPairs.h"

// a couple of macros for packing both sense and antisense fixed length barcodes
#define PACKEDBARCODE(b1,b2,b3,b4,b5,b6) (UINT32)((b1 << 10) | (b2 << 8) | (b3 << 6) | (b4 << 4) | (b5 << 2) | b6)
#define REVCPLPACKEDBARCODE(b1,b2,b3,b4,b5,b6) (UINT32)(((0x03 ^ b6) << 10) | ((0x03 ^ b5) << 8) | ((0x03 ^ b4) << 6) | ((0x03 ^ b3) << 4) | ((0x03 ^ b2) << 2) | (0x03 ^ b1))

static tsBarcode Barcodes[] = { // default hardcoded barcodes 
		// column barcodes
		{0,1,6,PACKEDBARCODE(eBaseA,eBaseA,eBaseC,eBaseC,eBaseA,eBaseA),REVCPLPACKEDBARCODE(eBaseA,eBaseA,eBaseC,eBaseC,eBaseA,eBaseA)},
		{0,2,6,PACKEDBARCODE(eBaseA,eBaseC,eBaseC,eBaseC,eBaseC,eBaseC),REVCPLPACKEDBARCODE(eBaseA,eBaseC,eBaseC,eBaseC,eBaseC,eBaseC)},
		{0,3,6,PACKEDBARCODE(eBaseA,eBaseG,eBaseC,eBaseC,eBaseG,eBaseG),REVCPLPACKEDBARCODE(eBaseA,eBaseG,eBaseC,eBaseC,eBaseG,eBaseG)},
		{0,4,6,PACKEDBARCODE(eBaseA,eBaseT,eBaseC,eBaseC,eBaseT,eBaseT),REVCPLPACKEDBARCODE(eBaseA,eBaseT,eBaseC,eBaseC,eBaseT,eBaseT)},
		{0,5,6,PACKEDBARCODE(eBaseA,eBaseA,eBaseG,eBaseG,eBaseA,eBaseA),REVCPLPACKEDBARCODE(eBaseA,eBaseA,eBaseG,eBaseG,eBaseA,eBaseA)},
		{0,6,6,PACKEDBARCODE(eBaseA,eBaseC,eBaseG,eBaseG,eBaseC,eBaseC),REVCPLPACKEDBARCODE(eBaseA,eBaseC,eBaseG,eBaseG,eBaseC,eBaseC)},
		{0,7,6,PACKEDBARCODE(eBaseA,eBaseG,eBaseG,eBaseG,eBaseG,eBaseG),REVCPLPACKEDBARCODE(eBaseA,eBaseG,eBaseG,eBaseG,eBaseG,eBaseG)},
		{0,8,6,PACKEDBARCODE(eBaseA,eBaseT,eBaseG,eBaseG,eBaseT,eBaseT),REVCPLPACKEDBARCODE(eBaseA,eBaseT,eBaseG,eBaseG,eBaseT,eBaseT)},
		{0,9,6,PACKEDBARCODE(eBaseA,eBaseA,eBaseT,eBaseT,eBaseA,eBaseA),REVCPLPACKEDBARCODE(eBaseA,eBaseA,eBaseT,eBaseT,eBaseA,eBaseA)},
		{0,10,6,PACKEDBARCODE(eBaseA,eBaseC,eBaseT,eBaseT,eBaseC,eBaseC),REVCPLPACKEDBARCODE(eBaseA,eBaseC,eBaseT,eBaseT,eBaseC,eBaseC)},
		{0,11,6,PACKEDBARCODE(eBaseA,eBaseG,eBaseT,eBaseT,eBaseG,eBaseG),REVCPLPACKEDBARCODE(eBaseA,eBaseG,eBaseT,eBaseT,eBaseG,eBaseG)},
		{0,12,6,PACKEDBARCODE(eBaseA,eBaseT,eBaseT,eBaseT,eBaseT,eBaseT),REVCPLPACKEDBARCODE(eBaseA,eBaseT,eBaseT,eBaseT,eBaseT,eBaseT)},
		// row barcodes
		{1,1,6,PACKEDBARCODE(eBaseT,eBaseA,eBaseA,eBaseT,eBaseA,eBaseA),REVCPLPACKEDBARCODE(eBaseT,eBaseA,eBaseA,eBaseT,eBaseA,eBaseA)},
		{1,2,6,PACKEDBARCODE(eBaseT,eBaseC,eBaseA,eBaseT,eBaseC,eBaseC),REVCPLPACKEDBARCODE(eBaseT,eBaseC,eBaseA,eBaseT,eBaseC,eBaseC)},
		{1,3,6,PACKEDBARCODE(eBaseT,eBaseG,eBaseA,eBaseT,eBaseG,eBaseG),REVCPLPACKEDBARCODE(eBaseT,eBaseG,eBaseA,eBaseT,eBaseG,eBaseG)},
		{1,4,6,PACKEDBARCODE(eBaseT,eBaseT,eBaseA,eBaseT,eBaseT,eBaseT),REVCPLPACKEDBARCODE(eBaseT,eBaseT,eBaseA,eBaseT,eBaseT,eBaseT)},
		{1,5,6,PACKEDBARCODE(eBaseT,eBaseA,eBaseT,eBaseA,eBaseA,eBaseA),REVCPLPACKEDBARCODE(eBaseT,eBaseA,eBaseT,eBaseA,eBaseA,eBaseA)},
		{1,6,6,PACKEDBARCODE(eBaseT,eBaseC,eBaseT,eBaseA,eBaseC,eBaseC),REVCPLPACKEDBARCODE(eBaseT,eBaseC,eBaseT,eBaseA,eBaseC,eBaseC)},
		{1,7,6,PACKEDBARCODE(eBaseT,eBaseG,eBaseT,eBaseA,eBaseG,eBaseG),REVCPLPACKEDBARCODE(eBaseT,eBaseG,eBaseT,eBaseA,eBaseG,eBaseG)},
		{1,8,6,PACKEDBARCODE(eBaseT,eBaseT,eBaseT,eBaseA,eBaseT,eBaseT),REVCPLPACKEDBARCODE(eBaseT,eBaseT,eBaseT,eBaseA,eBaseT,eBaseT)}
	};
static int NumBarcodes = (sizeof(Barcodes)/sizeof(tsBarcode));  // number of barcodes for which row and column barcodes are declared
static int NumPlateWells = 96;    // number of plate wells identified by barcodes (rows * columns)



char *
CMergeReadPairs::RemoveQuotes(char *pszRawText)
{
	char *pSrcChr;
	char *pDstChr;
	bool bInSpace;
	char Chr;
	CUtility::TrimQuotedWhitespcExtd(pszRawText);
	pSrcChr = pszRawText;
	pDstChr = pSrcChr;
	bInSpace = false;
	while ((Chr = *pSrcChr++) != '\0')
	{
		if (Chr == '\'' || Chr == '"')
			continue;
		if (Chr == ' ' || Chr == '\t')
		{
			if (bInSpace)
				continue;
			bInSpace = true;
		}
		else
			bInSpace = false;
		*pDstChr++ = Chr;
	}
	*pDstChr = '\0';
	return(pszRawText);
}

int										// returns number of MIDs initialised or eBSFErrxxxx
CMergeReadPairs::InitMIDs(bool bPE,		// true if processing for PE MIDs else if false then processing SE MIDs
				 char *pszMIDsBarcodeFile) // initialise with SE MIDs barcodes and associated name identifiers from a MIDs CSV file 
{
int Rslt;
int Idx;

int NumFields;
int NumElsParsed;
int m_NumMIDs;
int MIDsNameLen;
int PE1BarcodeLen;
etSeqBase PE1Barcode[cMaxBarCodeLen + 1];
int PE2BarcodeLen;
etSeqBase PE2Barcode[cMaxBarCodeLen + 1];
char *pszMIDsName;
char *pszMIDsBarcode;
char *pBase;
tsMIDsBarcode *pMIDsBarcode;


CCSVFile *pCSV = new CCSVFile;
if (pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}

if ((Rslt = pCSV->Open(pszMIDsBarcodeFile)) != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to open MIDs barcode file: %s", pszMIDsBarcodeFile);
	delete pCSV;
	return(Rslt);
}

NumElsParsed = 0;
NumFields = 0;
m_NumMIDs = 0;

// expected MIDs format is CSV with:
// Col1: 'Name',
// Col2: "PE1 Barcode" SE or 5'PE1 if PE processing
// Col3: "PE2 Barcode" PE2 if PE processing
while((Rslt = pCSV->NextLine()) > 0)
	{
	NumFields = pCSV->GetCurFields();
	if (!NumElsParsed && (((!bPE && NumFields >= 2) || (bPE && NumFields >= 3)) && pCSV->IsLikelyHeaderLine())) // check for header line
		continue;
	if(!bPE && NumFields < 2)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "SE processing, expected MIDs name and SE 5' barcode at line %d in file '%s'", pCSV->GetLineNumber(), pszMIDsBarcodeFile);
		delete pCSV;
		return(eBSFerrParse);
		}
	if(bPE && NumFields < 3)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "PE processing, expected MIDs name, PE1 barcode and PE2 barcode at line %d in file '%s'", pCSV->GetLineNumber(), pszMIDsBarcodeFile);
		delete pCSV;
		return(eBSFerrParse);
		}
	
	NumElsParsed += 1;
	pCSV->GetText(1, &pszMIDsName);
	pCSV->GetText(2, &pszMIDsBarcode);
	RemoveQuotes(pszMIDsName);
	RemoveQuotes(pszMIDsBarcode);

	if(((MIDsNameLen = (int)strlen(pszMIDsName)) > cMaxLenMIDsName) || MIDsNameLen < 2)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected MIDs name '%s' at line %d in file '%s' to be in length range 3..%d", pszMIDsName, pCSV->GetLineNumber(), pszMIDsBarcodeFile, cMaxLenMIDsName);
		delete pCSV;
		return(eBSFerrParse);
		}
	if (((PE1BarcodeLen = (int)strlen(pszMIDsBarcode)) > cMaxBarCodeLen) || PE1BarcodeLen < 4)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected MIDs barcode '%s' at line %d in file '%s' to be in length range 4..%d", pszMIDsBarcode, pCSV->GetLineNumber(), pszMIDsBarcodeFile, cMaxBarCodeLen);
		delete pCSV;
		return(eBSFerrParse);
		}

	pBase = pszMIDsBarcode;
	for(Idx = 0; Idx < PE1BarcodeLen; pBase++)
		{
		switch(*pBase) {
			case 'a': case 'A': case 'c': case 'C': case 'g': case 'G': case 't': case 'T':
				continue;
			default:
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected MIDs barcode '%s' at line %d in file '%s' to only contain canonical bases ACGT", pszMIDsBarcode, pCSV->GetLineNumber(), pszMIDsBarcodeFile);
				delete pCSV;
				return(eBSFerrParse);
			}
		}

	CSeqTrans::MapAscii2Sense(pszMIDsBarcode, PE1BarcodeLen, PE1Barcode);
	PE1Barcode[PE1BarcodeLen] = eBaseEOS;

	PE2Barcode[0] = eBaseEOS;
	PE2BarcodeLen = 0;
 
	if(bPE)
		{
		pCSV->GetText(3, &pszMIDsBarcode);
		RemoveQuotes(pszMIDsBarcode);
		if (((PE2BarcodeLen = (int)strlen(pszMIDsBarcode)) > cMaxBarCodeLen) || PE2BarcodeLen < 4)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected MIDs PE2 barcode '%s' at line %d in file '%s' to be in length range 4..%d", pszMIDsBarcode, pCSV->GetLineNumber(), pszMIDsBarcodeFile, cMaxBarCodeLen);
			delete pCSV;
			return(eBSFerrParse);
			}
		pBase = pszMIDsBarcode;
		for (Idx = 0; Idx < PE2BarcodeLen; pBase++)
			{
			switch (*pBase)
				{
				case 'a': case 'A': case 'c': case 'C': case 'g': case 'G': case 't': case 'T':
					continue;
				default:
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Expected MIDs PE2 barcode '%s' at line %d in file '%s' to only contain canonical bases ACGT", pszMIDsBarcode, pCSV->GetLineNumber(), pszMIDsBarcodeFile);
					delete pCSV;
					return(eBSFerrParse);
				}
			}

		CSeqTrans::MapAscii2Sense(pszMIDsBarcode, PE2BarcodeLen, PE2Barcode);
		PE2Barcode[PE2BarcodeLen] = eBaseEOS;
		if (!memcmp(PE1Barcode, PE2Barcode, min(PE1BarcodeLen, PE2BarcodeLen)))
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "Duplicate MIDs PE1 and PE2 barcodes at line %d in file '%s'", pCSV->GetLineNumber(), pszMIDsBarcodeFile);
			delete pCSV;
			return(eBSFerrParse);
			}
		}


	// check that this is the only instance of either name or barcodes - these must be globally unique
	if(m_NumMIDs > 0)
		{
		pMIDsBarcode = m_MIDsBarcodes;
		for(Idx = 0; Idx < m_NumMIDs; Idx+=1, pMIDsBarcode += 1)
			{
			if(!memcmp(pMIDsBarcode->PE1Barcode,PE1Barcode,min(pMIDsBarcode->PE1BarcodeLen,PE1BarcodeLen)) || 
			   !memcmp(pMIDsBarcode->PE2Barcode, PE2Barcode, min(pMIDsBarcode->PE2BarcodeLen, PE2BarcodeLen)) ||
			   !memcmp(pMIDsBarcode->PE1Barcode, PE2Barcode, min(pMIDsBarcode->PE1BarcodeLen, PE2BarcodeLen)) ||
			   !memcmp(pMIDsBarcode->PE2Barcode, PE1Barcode, min(pMIDsBarcode->PE2BarcodeLen, PE1BarcodeLen)) ||
				!stricmp(pMIDsBarcode->szMIDsName, pszMIDsName))
				{
				gDiagnostics.DiagOut(eDLFatal, gszProcName, "Duplicate name or MIDs barcode at line %d in file '%s'", pCSV->GetLineNumber(), pszMIDsBarcodeFile);
				delete pCSV;
				return(eBSFerrParse);
				}
			}
		}

	pMIDsBarcode = &m_MIDsBarcodes[m_NumMIDs++];
	pMIDsBarcode->MIDsID = m_NumMIDs;
	pMIDsBarcode->PE1BarcodeLen = PE1BarcodeLen;
	memcpy(pMIDsBarcode->PE1Barcode, PE1Barcode, PE1BarcodeLen + 1);
	strcpy(pMIDsBarcode->szMIDsName, pszMIDsName);

	if(bPE)
		{
		pMIDsBarcode->PE2BarcodeLen = PE2BarcodeLen;
		memcpy(pMIDsBarcode->PE2Barcode, PE2Barcode, PE2BarcodeLen + 1);
		}
	else
		{
		pMIDsBarcode->PE2BarcodeLen = 0;
		pMIDsBarcode->PE2Barcode[0] = eBaseEOS;
		}
	}

delete pCSV;

return(m_NumMIDs);
}



int				// returned best uniquely matching MIDsID or 0 if unable to match any MIDs barcode allowing at most m_MaxMIPsBarcodeSubs substitutions
CMergeReadPairs::MapSEMIDsBarcode(int SeqLen,			// number of bases in sequence putatively containing MIDs barcode 
				 int Ofs,				// look for matching MIDs barcode starting at this pSeq 5' offset
				 int MaxAllowSubs,		// max allowed substitutions
				 int *pReqSubs,			// returned best match required this number of substitutions
				 etSeqBase *pSeq)		// sequence
{
int CurMIDsID;
int BestMIDsID;
int CurBarcodeLen;
int CurNumSubs;
int LowestNumSubs;
bool bUniqueLowestSubs;
tsMIDsBarcode *pMIDs;
etSeqBase *pBase1;
etSeqBase *pBase;
etSeqBase *pBarcodeBase;

if(pReqSubs != NULL)
	*pReqSubs = 0;
if(pSeq == NULL || SeqLen < (Ofs + cMinBCReadSeqLen) || MaxAllowSubs < 0 || MaxAllowSubs > 5 || m_NumMIDsBarcodes < 1)
	return(0);
	
BestMIDsID = 0;
LowestNumSubs = 0;
bUniqueLowestSubs = true;
pMIDs = m_MIDsBarcodes;
pBase1 = &pSeq[Ofs];
for(CurMIDsID = 1; CurMIDsID <= m_NumMIDsBarcodes; CurMIDsID +=1, pMIDs += 1)
	{
	pBase = pBase1;
	pBarcodeBase = pMIDs->PE1Barcode;
	CurNumSubs = 0;
	for(CurBarcodeLen = 0; CurBarcodeLen < pMIDs->PE1BarcodeLen; CurBarcodeLen+=1, pBase+=1, pBarcodeBase+=1)
		{
		if((*pBase & 0x07) != (*pBarcodeBase & 0x07))
			{
			if(++CurNumSubs > MaxAllowSubs)
				break;
			}
		}
	if (CurNumSubs == 0)  // can't do any better than an exact match with no subs required
		return(CurMIDsID);
	if(CurNumSubs > MaxAllowSubs)	// too many subs required?
		continue;

	if(BestMIDsID == 0 || CurNumSubs <= LowestNumSubs) // will be selecting match with lowest number of subs to be returned
		{
		bUniqueLowestSubs = CurNumSubs < LowestNumSubs ? true : false; // ensure that only unique matches are returned
		LowestNumSubs = CurNumSubs;
		BestMIDsID = CurMIDsID;
		}	
	}

if(bUniqueLowestSubs)
	{
	if(pReqSubs != NULL)
		*pReqSubs = LowestNumSubs;
	return(BestMIDsID);
	}
return(0);
}

int											// returned best uniquely matching MIDsID or 0 if unable to match any MIDs barcode allowing at most MaxAllowSubs substitutions
CMergeReadPairs::MapPEMIDsBarcode(int SeqLen,			// minimum length of either pPE1Seq or pPE2Seq  
				 int Ofs,				// look for matching MIDs barcode starting at this pPE1Seq or pPE2Seq 5' offset
				 int MaxAllowSubs,		// max allowed substitutions for either PE1 or PE2 (not total over both PE1 and PE2)
				 int *pReqSubs,			// best match required this total number of substitutions over both PE1 and PE2 barcodes
				 etSeqBase *pPE1Seq,	// map barcodes at 5' of this sequence
				 etSeqBase *pPE2Seq)	// and barcodes at 5'of this sequence to the MIDs
{
	int CurMIDsID;
	int BestMIDsID;
	int CurBarcodeLen;
	int CurNumSubs;
	int CurTotNumSubs;
	int LowestNumSubs;
	bool bUniqueLowestSubs;
	tsMIDsBarcode *pMIDs;
	etSeqBase *pBase1;
	etSeqBase *pBase2;
	etSeqBase *pBase;
	etSeqBase *pBarcodeBase;

	if (pReqSubs != NULL)
		*pReqSubs = 0;
	if (pPE1Seq == NULL || pPE2Seq == NULL || SeqLen < (Ofs + cMinBCReadSeqLen) || MaxAllowSubs < 0 || MaxAllowSubs > 5 || m_NumMIDsBarcodes < 1)
		return(0);

	BestMIDsID = 0;
	LowestNumSubs = 0;
	CurTotNumSubs = 0;
	bUniqueLowestSubs = true;
	pMIDs = m_MIDsBarcodes;
	pBase1 = &pPE1Seq[Ofs];
	pBase2 = &pPE2Seq[Ofs];
	for (CurMIDsID = 1; CurMIDsID <= m_NumMIDsBarcodes; CurMIDsID += 1, pMIDs += 1)
		{
		CurTotNumSubs = 0;
		CurNumSubs = 0;
		pBase = pBase1;
		pBarcodeBase = pMIDs->PE1Barcode;
		for (CurBarcodeLen = 0; CurBarcodeLen < pMIDs->PE1BarcodeLen; CurBarcodeLen += 1, pBase += 1, pBarcodeBase += 1)
			{
			if ((*pBase & 0x07) != (*pBarcodeBase & 0x07))
				{
				if (++CurNumSubs > MaxAllowSubs)
					break;
				}
			}

		if (CurNumSubs > MaxAllowSubs)
			continue;
		CurTotNumSubs = CurNumSubs;
		CurNumSubs = 0;
		pBase = pBase2;
		pBarcodeBase = pMIDs->PE2Barcode;
		for (CurBarcodeLen = 0; CurBarcodeLen < pMIDs->PE2BarcodeLen; CurBarcodeLen += 1, pBase += 1, pBarcodeBase += 1)
			{
			if ((*pBase & 0x07) != (*pBarcodeBase & 0x07))
				{
				if (++CurNumSubs > MaxAllowSubs)
					break;
				}
			}
		if (CurNumSubs > MaxAllowSubs)
			continue;
		CurTotNumSubs += CurNumSubs;

		if (CurTotNumSubs == 0)  // can't do any better than an exact match with no subs required
			return(CurMIDsID);

		if (BestMIDsID == 0 || CurTotNumSubs <= LowestNumSubs) // will be selecting match with lowest number of subs to be returned
			{
			bUniqueLowestSubs = CurTotNumSubs < LowestNumSubs ? true : false; // ensure that only unique matches are returned
			LowestNumSubs = CurTotNumSubs;
			BestMIDsID = CurMIDsID;
			}
	}

	if (bUniqueLowestSubs)
	{
		if (pReqSubs != NULL)
			*pReqSubs = LowestNumSubs;
		return(BestMIDsID);
	}
	return(0);
}



int							// returned well number (1..96) or 0 if unable to identify well from the barcodes
CMergeReadPairs::MapSEBarcodesToWell(int SeqLen,		// num bases in SE pSeq
				etSeqBase *pSeq)	// map barcodes at 5' and 3' end of this SE sequence to the well
{
UINT32 Pack5;		// to hold the 5' extracted assumed barcode from pSeq
UINT32 Pack3;		// to hold the 3' extracted assumed barcode from pSeq
int Idx;
etSeqBase *pBase;
tsBarcode *pBarcode;
tsBarcode *pBarcode5;
tsBarcode *pBarcode3;

if(m_pBarcodes == NULL || (m_MaxBarcode5Len + m_MaxBarcode3Len) == 0 || SeqLen < (m_MaxBarcode5Len + m_MaxBarcode3Len))
	return(0);

Pack5 = 0;
Pack3 = 0;
// extract the 5' barcode 
pBase= pSeq;
for(Idx = 0; Idx < m_MaxBarcode5Len; Idx++, pBase++)
	{
	if(*pBase > eBaseT)			// only accepting canonical bases
		return(0);
	Pack5 <<= 2;
	Pack5 |= *pBase & 0x03;
	}		

// extract the 3' barcode
pBase= &pSeq[SeqLen-m_MaxBarcode3Len];
for(Idx = 0; Idx < m_MaxBarcode3Len; Idx++, pBase++)
	{
	if(*pBase > eBaseT)			// only accepting canonical bases
		return(0);
	Pack3 <<= 2;
	Pack3 |= *pBase & 0x03;
	}

// try matching the barcodes
pBarcode5 = NULL;
pBarcode3 = NULL;
pBarcode = m_pBarcodes;
for(Idx = 0; Idx < m_NumBarcodes; Idx++,pBarcode++)
	{
	if(pBarcode->BarCode == Pack5 && !pBarcode->ColRow)	
		pBarcode5 = pBarcode;

	if(pBarcode->RevCplBarCode == Pack3 && pBarcode->ColRow)
		pBarcode3 = pBarcode;

	if(pBarcode5 != NULL && pBarcode3 != NULL)
		return(pBarcode5->Psn + ((pBarcode3->Psn-1) * 12));
	}

pBarcode5 = NULL;
pBarcode3 = NULL;
pBarcode = m_pBarcodes;
for(Idx = 0; Idx < m_NumBarcodes; Idx++,pBarcode++)
	{
	if(pBarcode->BarCode == Pack5 && pBarcode->ColRow)
		pBarcode5 = pBarcode;

	if(pBarcode->RevCplBarCode == Pack3 && !pBarcode->ColRow)
		pBarcode3 = pBarcode;

	if(pBarcode5 != NULL && pBarcode3 != NULL)
		return(pBarcode3->Psn + ((pBarcode5->Psn-1) * 12));
	}

return(0);
}


int							// returned well number (1..96) or 0 if unable to identify well from the barcodes
CMergeReadPairs::MapPEBarcodesToWell(int SeqLen,		//minimum length of either p5Seq or P3Seq
								   etSeqBase *pPE1Seq,// map barcodes at 5' of this sequence
								   etSeqBase *pPE2Seq)// and barcodes at 5'of this sequence to the well
{
	UINT32 Pack5;		// to hold the 5' extracted assumed barcode from pSeq
	UINT32 Pack3;		// to hold the 3' extracted assumed barcode from pSeq
	int Idx;
	etSeqBase *pBase;
	tsBarcode *pBarcode;
	tsBarcode *pBarcode5;
	tsBarcode *pBarcode3;

	if (m_pBarcodes == NULL || (m_MaxBarcode5Len + m_MaxBarcode3Len) == 0 || SeqLen < m_MaxBarcode5Len || SeqLen < m_MaxBarcode3Len)
		return(0);

	Pack5 = 0;
	Pack3 = 0;
	// extract the PE1 5' barcode 
	pBase = pPE1Seq;
	for (Idx = 0; Idx < m_MaxBarcode5Len; Idx++, pBase++)
		{
		if (*pBase > eBaseT)			// only accepting cannonical bases
			return(0);
		Pack5 <<= 2;
		Pack5 |= *pBase & 0x03;
		}

	// extract the PE2 5' barcode
	pBase = pPE2Seq;
	for (Idx = 0; Idx < m_MaxBarcode3Len; Idx++, pBase++)
		{
		if (*pBase > eBaseT)			// only accepting cannonical bases
			return(0);
		Pack3 <<= 2;
		Pack3 |= *pBase & 0x03;
		}

	// try matching the barcodes
	// firstly with amplicon 5' sense and amplicon 3' antisense; if no match then try 5' antisense and 3' sense
pBarcode5 = NULL;
pBarcode3 = NULL;
pBarcode = m_pBarcodes;

for (Idx = 0; Idx < m_NumBarcodes; Idx++, pBarcode++)
	{
	if (pBarcode->BarCode == Pack5 && !pBarcode->ColRow)
		pBarcode5 = pBarcode;

	if (pBarcode->BarCode == Pack3 && pBarcode->ColRow)
		pBarcode3 = pBarcode;

	if (pBarcode5 != NULL && pBarcode3 != NULL)
		return(pBarcode5->Psn + ((pBarcode3->Psn - 1) * 12));
	}

pBarcode5 = NULL;
pBarcode3 = NULL;
pBarcode = m_pBarcodes;
for (Idx = 0; Idx < m_NumBarcodes; Idx++, pBarcode++)
	{
	if (pBarcode->BarCode == Pack5 && pBarcode->ColRow)
		pBarcode5 = pBarcode;

	if (pBarcode->BarCode == Pack3 && !pBarcode->ColRow)
		pBarcode3 = pBarcode;
	
	if (pBarcode5 != NULL && pBarcode3 != NULL)
		return(pBarcode3->Psn + ((pBarcode5->Psn - 1) * 12));
	}

return(0);
}

int				// return number of wells initialised
CMergeReadPairs::InitDfltWells(bool bNoMerge) // initialise with default well barcodes and well identifiers, if bNoMerge then do not merge PE reads and report PE1/PE2 instead of merged SE 
{
tsAmpliconWell *pWell;
int WellIdx;
m_MaxBarcode5Len = 6;		// maximun length of any 5' barcode - max 16
m_MaxBarcode3Len = 6;		// maximun length of any 3' barcode - max 16
m_NumBarcodes = NumBarcodes;	// number of barcodes in m_pBarcodes[]
m_pBarcodes = Barcodes;
m_NumWells = 0;
pWell = m_WellFiles;
memset(pWell,0,sizeof(m_WellFiles));
for(WellIdx = 0; WellIdx < NumPlateWells; WellIdx++,pWell++)
	{
	if(bNoMerge)
		{
		sprintf(pWell->WellFile[1].szOutFile, "%s.Well%d.PE2.%s", m_szMergeOutFile, WellIdx + 1, m_OFormat == eOFfasta ? "fasta" : "fastq");
		if ((pWell->WellFile[1].pOutBuffer = new char[cAllocOutBuffLen]) == NULL)
			return(eBSFerrMem);
		pWell->WellFile[1].AllocdOutBuff = cAllocOutBuffLen;
		}

	sprintf(pWell->WellFile[0].szOutFile,"%s.Well%d.%s.%s",m_szMergeOutFile,WellIdx+1, bNoMerge == true ? "PE1" : "SE",m_OFormat == eOFfasta ? "fasta" : "fastq");
	pWell->WellFile[0].hOutFile = -1;
	pWell->WellFile[1].hOutFile = -1;
	if((pWell->WellFile[0].pOutBuffer = new char [cAllocOutBuffLen]) == NULL)
		return(eBSFerrMem);
	pWell->WellFile[0].AllocdOutBuff = cAllocOutBuffLen;

	pWell->WellID = WellIdx + 1;
	m_NumWells += 1;
	}

return(m_NumWells);
}

CMergeReadPairs::CMergeReadPairs(void)
{
m_pszMSeqs = NULL;
m_pszUnmergedP1Seqs = NULL;
m_pszUnmergedP2Seqs = NULL;
m_hOutMerged = -1;
m_hOut5Unmerged = -1;
m_hOut3Unmerged = -1;
m_ProcPhase = ePPUninit;
m_NumWells = 0;
memset(m_WellFiles,0,sizeof(m_WellFiles));
for(int WellIdx = 0; WellIdx < m_NumWells; WellIdx++)
	{
	m_WellFiles[WellIdx].WellFile[0].hOutFile = -1;
	m_WellFiles[WellIdx].WellFile[1].hOutFile = -1;
	}
Reset(false);
}


CMergeReadPairs::~CMergeReadPairs(void)
{
Reset(false);
}

void
CMergeReadPairs::Reset(bool bSync)
{
int WellIdx;
int WellFileIdx;
tsAmpliconWellFile *pWellFile;

tsAmpliconWell *pWell;
if(m_ProcPhase != ePPUninit)
	{
	if(m_NumWells > 0)
		{
		pWell = m_WellFiles;
		for(WellIdx = 0; WellIdx < m_NumWells; WellIdx++,pWell++)
			{
			pWellFile = &pWell->WellFile[0];
			for(WellFileIdx = 0; WellFileIdx < 2; WellFileIdx += 1, pWellFile++)
				{
				if(pWellFile->pOutBuffer != NULL)
					{
					if(bSync && pWellFile->CurBuffLen && pWellFile->hOutFile != -1)
						CUtility::SafeWrite(pWellFile->hOutFile, pWellFile->pOutBuffer, pWellFile->CurBuffLen);
					delete pWellFile->pOutBuffer;
					pWellFile->pOutBuffer = NULL;
					}
				pWellFile->AllocdOutBuff = 0;
				if(pWellFile->hOutFile != -1)
					{
					if(bSync)
#ifdef _WIN32
					_commit(pWellFile->hOutFile);
#else
					fsync(pWellFile->hOutFile);
#endif
					close(pWellFile->hOutFile);
					pWellFile->hOutFile = -1;
					}
				}
			}
		}

	if(m_pszMSeqs != NULL)
		{
		if(bSync && m_CurMSeqLen && m_hOutMerged != -1)
			CUtility::SafeWrite(m_hOutMerged,m_pszMSeqs,m_CurMSeqLen);
		delete m_pszMSeqs;
		m_pszMSeqs = NULL;
		}
	m_AllocdMSeqs = 0;

	if(m_pszUnmergedP1Seqs != NULL)
		{
		if(bSync && m_CurUnmergedP1Seqs && m_hOut5Unmerged != -1)
			CUtility::SafeWrite(m_hOut5Unmerged,m_pszUnmergedP1Seqs,m_CurUnmergedP1Seqs);
		delete m_pszUnmergedP1Seqs;
		m_pszUnmergedP1Seqs = NULL;
		}
	m_AllocdUnmergedP1Seqs = 0;

	if(bSync && m_pszUnmergedP2Seqs != NULL)
		{
		if(m_CurUnmergedP2Seqs && m_hOut3Unmerged != -1)
			CUtility::SafeWrite(m_hOut3Unmerged,m_pszUnmergedP2Seqs,m_CurUnmergedP2Seqs);
		delete m_pszUnmergedP2Seqs;
		m_pszUnmergedP2Seqs = NULL;
		}
	m_AllocdUnmergedP2Seqs = 0;

	if(m_hOutMerged != -1)
		{
		if(bSync)
#ifdef _WIN32
			_commit(m_hOutMerged);
#else
			fsync(m_hOutMerged);
#endif
		close(m_hOutMerged);
		m_hOutMerged = -1;
		}

	if(m_hOut5Unmerged != -1)
		{
		if(bSync)
#ifdef _WIN32
			_commit(m_hOut5Unmerged);
#else
			fsync(m_hOut5Unmerged);
#endif
		close(m_hOut5Unmerged);
		m_hOut5Unmerged = -1;
		}

	if(m_hOut3Unmerged != -1)
		{
		if(bSync)
#ifdef _WIN32
			_commit(m_hOut3Unmerged);
#else
			fsync(m_hOut3Unmerged);
#endif
		close(m_hOut3Unmerged);
		m_hOut3Unmerged = -1;
		}
	}

memset(m_WellFiles,0,sizeof(m_WellFiles));
pWell = m_WellFiles;
for(WellIdx = 0; WellIdx < cMaxNumBarcodes; WellIdx++,pWell++)
	{
	pWell->WellFile[0].hOutFile = -1;
	pWell->WellFile[1].hOutFile = -1;
	}
m_bAmpliconNoMerge = false;
m_hOutMerged = -1;			// file handle for output merged read microcontigs
m_hOut5Unmerged = -1;		// file handle for output 5' unmerged reads
m_hOut3Unmerged = -1;		// file handle for output 3' unmerged reads
m_pszMSeqs = NULL;			// will be allocated for holding merged sequences ready for writing to file
m_pszUnmergedP1Seqs = NULL;  // will be allocated for holding unmerged P1 sequences ready for writing to file
m_pszUnmergedP2Seqs = NULL;  // will be allocated for holding unmerged P2 sequences ready for writing to file
m_AllocdMSeqs = 0;
m_AllocdUnmergedP1Seqs = 0;
m_AllocdUnmergedP2Seqs = 0;

m_szIn5ReadsFile[0] = 0;	// 5' reads are in this file
m_szIn3ReadsFile[0] = 0;	// 3' reads are in this file
m_szMergeOutFile[0] = 0;	// write merged overlaping reads to this file
m_sz5UnmergedOutFile[0] = 0;	// write 5' unmerged reads to this file
m_sz3UnmergedOutFile[0] = 0;   // write 3' unmerged reads to this file

m_PMode = ePMdefault;
m_OFormat = eOFauto;		
m_bIsFastq = false;			// if true then input reads are fastq
m_bAppendOut = false;		// if true then append if output files exist, otherwise trunctate existing output files 
m_MinOverlap = 0;			// reads must overlap by at least this many bases
m_MaxOverlapPropSubs = 0;	// and overlap can have at most this proportion of substitutions; if non-zero then a floor of 1 sub is allowed

m_CurMSeqLen = 0;			// m_pszMSeqs currently holds this many chars
m_CurUnmergedP1Seqs = 0;	// m_pszUnmergedP1Seqs currently holds this many chars
m_CurUnmergedP2Seqs = 0;	// m_pszUnmergedP2Seqs currently holds this many chars
m_ProcPhase = ePPReset;
}

// OpenFiles
// Initialise and open files for processing
int
CMergeReadPairs::OpenFiles(void)
{
int WellIdx;
tsAmpliconWell *pWellFileOut;

// check if files specified for open/create/append...
if(m_ProcPhase != ePPRdyOpen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"OpenFiles: unexpected state error");
	return(eBSFerrInternal);
	}

if(m_bAppendOut)
	{
	if(m_PMode < ePMAmplicon)
		{
#ifdef _WIN32
		m_hOutMerged = open(m_szMergeOutFile,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT),(_S_IREAD | _S_IWRITE));
#else
		m_hOutMerged = open(m_szMergeOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE);
#endif
		if(m_hOutMerged == -1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_szMergeOutFile);
			Reset(false);
			return(eBSFerrCreateFile);
			}
		// seek to end of file ready for appending
		lseek(m_hOutMerged,0,SEEK_END);
		}
	else
		{
		pWellFileOut = m_WellFiles;
		for(WellIdx = 0; WellIdx < m_NumWells; WellIdx++,pWellFileOut++)
			{
#ifdef _WIN32
			pWellFileOut->WellFile[0].hOutFile = open(pWellFileOut->WellFile[0].szOutFile,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT),(_S_IREAD | _S_IWRITE));
#else
				pWellFileOut->WellFile[0].hOutFile = open(pWellFileOut->WellFile[0].szOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE);
#endif
			if(pWellFileOut->WellFile[0].hOutFile == -1)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output sequences file '%s'", pWellFileOut->WellFile[0].szOutFile);
				Reset(false);
				return(eBSFerrCreateFile);
				}
			// seek to end of file ready for appending
			lseek(pWellFileOut->WellFile[0].hOutFile,0,SEEK_END);

			if(m_bAmpliconNoMerge == true)
				{
#ifdef _WIN32
				pWellFileOut->WellFile[1].hOutFile = open(pWellFileOut->WellFile[1].szOutFile, (_O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT), (_S_IREAD | _S_IWRITE));
#else
				pWellFileOut->WellFile[1].hOutFile = open(pWellFileOut->WellFile[1].szOutFile, O_RDWR | O_CREAT, S_IREAD | S_IWRITE);
#endif
				if (pWellFileOut->WellFile[1].hOutFile == -1)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: unable to create/truncate output sequences file '%s'", pWellFileOut->WellFile[1].szOutFile);
					Reset(false);
					return(eBSFerrCreateFile);
					}
				// seek to end of file ready for appending
				lseek(pWellFileOut->WellFile[1].hOutFile, 0, SEEK_END);
				}
			else
				pWellFileOut->WellFile[1].hOutFile = -1;
			}
		}

	if(m_sz5UnmergedOutFile[0] != '\0')
		{
#ifdef _WIN32
		m_hOut5Unmerged = open(m_sz5UnmergedOutFile,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT),(_S_IREAD | _S_IWRITE));
#else
		m_hOut5Unmerged = open(m_sz5UnmergedOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE);
#endif
		if(m_hOut5Unmerged == -1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_sz5UnmergedOutFile);
			Reset(false);
			return(eBSFerrCreateFile);
			}
		// seek to end of file ready for appending
		lseek(m_hOut5Unmerged,0,SEEK_END);
		}

	if(m_sz3UnmergedOutFile[0] != '\0')
		{
#ifdef _WIN32
		m_hOut3Unmerged = open(m_sz3UnmergedOutFile,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT),(_S_IREAD | _S_IWRITE));
#else
		m_hOut3Unmerged = open(m_sz3UnmergedOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE);
#endif
		if(m_hOut3Unmerged == -1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output 3' unmerged file '%s'",m_sz3UnmergedOutFile);
			Reset(false);
			return(eBSFerrCreateFile);
			}
		// seek to end of file ready for appending
		lseek(m_hOut3Unmerged,0,SEEK_END);
		}
	}
else
	{
	if(m_PMode < ePMAmplicon)
		{
#ifdef _WIN32
		m_hOutMerged = open(m_szMergeOutFile,O_CREATETRUNC );
#else
		if((m_hOutMerged = open(m_szMergeOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hOutMerged,0)!=0)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s merged PE microcontigs - %s",m_szMergeOutFile,strerror(errno));
					Reset(false);
					return(eBSFerrCreateFile);
					}
#endif
		if(m_hOutMerged == -1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_szMergeOutFile);
			Reset(false);
			return(eBSFerrCreateFile);
			}
		lseek(m_hOutMerged,0,SEEK_SET);
		}
	else
		{
		pWellFileOut = m_WellFiles;
		for(WellIdx = 0; WellIdx < m_NumWells; WellIdx++,pWellFileOut++)
			{
#ifdef _WIN32
			pWellFileOut->WellFile[0].hOutFile = open(pWellFileOut->WellFile[0].szOutFile,O_CREATETRUNC);
#else
			if((pWellFileOut->WellFile[0].hOutFile = open(pWellFileOut->WellFile[0].szOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
				if(ftruncate(pWellFileOut->WellFile[0].hOutFile,0)!=0)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s output sequences file - %s", pWellFileOut->WellFile[0].szOutFile,strerror(errno));
					Reset(false);
					return(eBSFerrCreateFile);
					}
#endif
			if(pWellFileOut->WellFile[0].hOutFile == -1)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output sequences file '%s'", pWellFileOut->WellFile[0].szOutFile);
				Reset(false);
				return(eBSFerrCreateFile);
				}
			// seek to end of file ready for appending
			lseek(pWellFileOut->WellFile[0].hOutFile,0,SEEK_END);


			if (m_bAmpliconNoMerge == true)
				{
#ifdef _WIN32
				pWellFileOut->WellFile[1].hOutFile = open(pWellFileOut->WellFile[1].szOutFile, O_CREATETRUNC);
#else
				if ((pWellFileOut->WellFile[1].hOutFile = open(pWellFileOut->WellFile[1].szOutFile, O_RDWR | O_CREAT, S_IREAD | S_IWRITE)) != -1)
					if (ftruncate(pWellFileOut->WellFile[1].hOutFile, 0) != 0)
						{
						gDiagnostics.DiagOut(eDLFatal, gszProcName, "Unable to truncate %s output sequences file - %s", pWellFileOut->WellFile[1].szOutFile, strerror(errno));
						Reset(false);
						return(eBSFerrCreateFile);
						}
#endif
				if (pWellFileOut->WellFile[1].hOutFile == -1)
					{
					gDiagnostics.DiagOut(eDLFatal, gszProcName, "Process: unable to create/truncate output sequences file '%s'", pWellFileOut->WellFile[1].szOutFile);
					Reset(false);
					return(eBSFerrCreateFile);
					}
					// seek to end of file ready for appending
				lseek(pWellFileOut->WellFile[1].hOutFile, 0, SEEK_END);
				}
			else
				pWellFileOut->WellFile[1].hOutFile = -1;
			}
		}

	if(m_sz5UnmergedOutFile[0] != '\0')
		{
#ifdef _WIN32
		m_hOut5Unmerged = open(m_sz5UnmergedOutFile,O_CREATETRUNC );
#else
		if((m_hOut5Unmerged = open(m_sz5UnmergedOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hOut5Unmerged,0)!=0)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s merged PE microcontigs - %s",m_sz5UnmergedOutFile,strerror(errno));
					Reset(false);
					return(eBSFerrCreateFile);
					}
#endif
		if(m_hOut5Unmerged == -1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_sz5UnmergedOutFile);
			Reset(false);
			return(eBSFerrCreateFile);
			}
		lseek(m_hOut5Unmerged,0,SEEK_SET);
		}

	if(m_sz3UnmergedOutFile[0] != '\0')
		{
#ifdef _WIN32
		m_hOut3Unmerged = open(m_sz3UnmergedOutFile,O_CREATETRUNC );
#else
		if((m_hOut3Unmerged = open(m_sz3UnmergedOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hOut3Unmerged,0)!=0)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s merged PE microcontigs - %s",m_sz3UnmergedOutFile,strerror(errno));
					Reset(false);
					return(eBSFerrCreateFile);
					}
#endif
		if(m_hOut3Unmerged == -1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_sz3UnmergedOutFile);
			Reset(false);
			return(eBSFerrCreateFile);
			}
		lseek(m_hOut3Unmerged,0,SEEK_SET);
		}
	}
m_ProcPhase = ePRdyProc;
return(eBSFSuccess);
}


int												// returns number of sequences merged and written to output files
CMergeReadPairs::MergeOverlaps(etPMode PMode,	// processing mode 
			  etOFormat OFormat,				// output file format
			  int MinOverlap,					// reads must overlap by at least this many bases
              int MaxOverlapPropSubs,			// and overlap can have at most this proportion of substitutions, if > 0 then floor of 1 sub allowed
			  int NumInPE5Files,				// number of input single ended or 5' end if paired end reads file
			  char **pszInPE5Files,				// input single ended or 5' end if paired end reads files
			  int NumInPE3Files,				// number of input input 3' end if paired end reads files
			  char **pszInPE3Files,				// input 3' end if paired end reads files
			  char *pszMergeOutFile,			// write merged overlaping reads to this file
  			  int StartNum,						// use this initial starting sequence identifier
  			  bool bAppendOut)				// if true then append if output files exist, otherwise trunctate existing output files

{
int Rslt;
int Idx;

if(MinOverlap < 1 ||					// must be at least a 1 base overlap required!
   MaxOverlapPropSubs < 0 ||			// must be reasonable number of allowed substitutions in overlap
   NumInPE5Files < 1 || NumInPE3Files < 0 ||
   NumInPE5Files != NumInPE3Files)
   	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"MergeOverlaps: inconsistency in function parameters");
	return(eBSFerrParams);
	}
Reset(false);	
m_PMode = PMode;
m_bAmpliconNoMerge = PMode == ePMAmpliconNoMerge ? true : false;

m_OFormat = OFormat;

m_bAppendOut = bAppendOut;
m_MinOverlap = MinOverlap;
m_MaxOverlapPropSubs = MaxOverlapPropSubs;
strncpy(m_szMergeOutFile,pszMergeOutFile,_MAX_PATH);
m_szMergeOutFile[_MAX_PATH-1] = '\0';

if(m_PMode >= ePMAmplicon)
	{
	if((Rslt = InitDfltWells(m_bAmpliconNoMerge)) < eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"MergeOverlaps: unable to initialise well barcode sequences");
		return(eBSFerrParams);
		}	
	}

if(PMode == ePMseparate)
	{
	strncpy(m_sz5UnmergedOutFile,m_szMergeOutFile,_MAX_PATH-20);
	strcat(m_sz5UnmergedOutFile,".UPE1");
	strncpy(m_sz3UnmergedOutFile,m_szMergeOutFile,_MAX_PATH-20);
	strcat(m_sz3UnmergedOutFile,".UPE2");
	}
else
	{
	m_sz5UnmergedOutFile[0] = '\0';
	m_sz3UnmergedOutFile[0] = '\0';
	}

m_ProcPhase = ePPRdyOpen;
for(Idx = 0; Idx < NumInPE5Files; Idx++)
	{
	strncpy(m_szIn5ReadsFile,pszInPE5Files[Idx],_MAX_PATH);
	m_szIn5ReadsFile[_MAX_PATH-1] = '\0';
	strncpy(m_szIn3ReadsFile,pszInPE3Files[Idx],_MAX_PATH);
	m_szIn3ReadsFile[_MAX_PATH-1] = '\0';

	if((Rslt=(teBSFrsltCodes)m_PE5Fasta.Open(m_szIn5ReadsFile,true))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"MergeOverlaps: Unable to open '%s' [%s] %s",m_szIn5ReadsFile,m_PE5Fasta.ErrText((teBSFrsltCodes)Rslt),m_PE5Fasta.GetErrMsg());
		Reset(false);
		return(Rslt);
		}

	if(m_PE5Fasta.IsSOLiD())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"MergeOverlaps: SOLiD colorspace processing is not supported, '%s' is colorspace...",m_szIn5ReadsFile);
		m_PE5Fasta.Close();
		Reset(false);
		return(eBSFerrOpnFile);
		}

	if((Rslt=(teBSFrsltCodes)m_PE3Fasta.Open(m_szIn3ReadsFile,true))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"MergeOverlaps: Unable to open '%s' [%s] %s",m_szIn3ReadsFile,m_PE3Fasta.ErrText((teBSFrsltCodes)Rslt),m_PE3Fasta.GetErrMsg());
		m_PE5Fasta.Close();
		Reset(false);
		return(Rslt);
		}

	if(m_PE3Fasta.IsSOLiD())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"MergeOverlaps: SOLiD colorspace processing is not supported, '%s' is colorspace...",m_szIn3ReadsFile);
		m_PE5Fasta.Close();
		m_PE3Fasta.Close();
		Reset(false);
		return(eBSFerrOpnFile);
		}

	m_bIsFastq = m_PE5Fasta.IsFastq();
	if(m_bIsFastq != m_PE3Fasta.IsFastq())
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"MergeOverlaps: Paired end files are not both either fastq or fasta type\n '%s' is %s, '%s' is %s",
									m_szIn5ReadsFile,m_bIsFastq ? "fastq" : "fasta",  m_szIn3ReadsFile,m_PE3Fasta.IsFastq() ? "fastq" : "fasta");
		m_PE5Fasta.Close();
		m_PE3Fasta.Close();
		Reset(false);
		return(eBSFerrOpnFile);
		}

	if(m_OFormat == eOFauto)
		m_OFormat = m_bIsFastq ? eOFfastq : eOFfasta;

	if(m_ProcPhase == ePPRdyOpen)
		{
		if((Rslt = OpenFiles()) < eBSFSuccess)
			{
			Reset(false);
			return(Rslt);
			}
		}

	if(PMode <= ePMAmplicon)
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"MergeOverlaps: Processing '%s' for overlaps with '%s'",m_szIn5ReadsFile,m_szIn3ReadsFile);
	else
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "MergeOverlaps: Processing '%s' for barcodes with '%s'", m_szIn5ReadsFile, m_szIn3ReadsFile);

	if((Rslt = ProcOverlapPairs(StartNum)) < eBSFSuccess)
		{
		m_PE5Fasta.Close();
		m_PE3Fasta.Close();
		Reset(false);
		return(Rslt);
		}

	if (PMode <= ePMAmplicon)
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "MergeOverlaps: Completed processing '%s' for overlaps with '%s'", m_szIn5ReadsFile, m_szIn3ReadsFile);
	else
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "MergeOverlaps: Completed processing '%s' for barcodes with '%s'", m_szIn5ReadsFile, m_szIn3ReadsFile);


	m_PE5Fasta.Close();
	m_PE3Fasta.Close();

	StartNum += Rslt;
	}

Reset(true);
return((int)Rslt);
}

const int cMinReadSeqLen = 20;			// raw read sequences must be at least this length
const int cMaxReadDesrLen = 200;		// allow for read descriptors of no longer than this length
const int ccMaxFastQSeqLen = 4000;		// allow for amplicon PE sequences of this length

int													// returns number of sequences merged and written to output file
CMergeReadPairs::ProcOverlapPairs(int StartNum)		// initial starting sequence identifier
{
int Rslt;
int Idx;

int WellIdx;

UINT8 PE5Seq[ccMaxFastQSeqLen+1];
UINT8 PE5Qual[ccMaxFastQSeqLen+1];
UINT8 PE3Seq[ccMaxFastQSeqLen+1];
UINT8 PE3Qual[ccMaxFastQSeqLen+1];
char szPE5DescrBuff[cMaxReadDesrLen+1];
char szPE3DescrBuff[cMaxReadDesrLen+1];
UINT8 MSeq[(2*ccMaxFastQSeqLen)+1];
UINT8 szMSeq[(2*ccMaxFastQSeqLen)+1];
UINT8 szMQual[(2*ccMaxFastQSeqLen)+1];
int PlateCell;

UINT64 Spacer1 = 0;
int PE5DescrLen;
int PE3DescrLen;

int NumPE5Reads;
int NumPE3Reads;
int MinObsOverlap;
int MaxObsOverlap;

int OverlapDistCnts[ccMaxFastQSeqLen];
int OverlapDistSubs[cMaxOverlapPercSubs+1];

int MergeIdx;
UINT8 *pMSeq;
UINT8 *pMSeqQual;
UINT8 *pSeq5Qual;
UINT8 *pSeq3Qual;
etSeqBase *pSeq5;
etSeqBase *pSeq3;
etSeqBase Base5;
etSeqBase Base3;
int OL5Idx;
int OL3Idx;
int MaxOverlap;
int OLStartsIdx;
int ReqOverlap3;
int MaxSubs;
int CurSubs;
int CurOvlpScore;
int BestOvlpScore;
int AllowedSubs;

int NumOverlapping;
int NumUnmapped;
int NumPENoBarcode;
int NumPEWithBarcode;
int PE5SeqLen;
int PE3SeqLen;
int PE5QualLen;
int PE3QualLen;


if(m_hOutMerged != -1 && m_pszMSeqs == NULL)
	{
	if((m_pszMSeqs = new char [cAllocOutBuffLen])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d) for buffering file writes",cAllocOutBuffLen);
		Reset(false);
		return(eBSFerrMem);
		}
	m_AllocdMSeqs = cAllocOutBuffLen;
	m_CurMSeqLen = 0;
	}

if(m_PMode == ePMseparate)
	{
	if(m_hOut5Unmerged != -1 && m_pszUnmergedP1Seqs == NULL)
		{
		if((m_pszUnmergedP1Seqs = new char [cAllocOutBuffLen])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d) for buffering file writes",cAllocOutBuffLen);
			Reset(false);
			return(eBSFerrMem);
			}
		m_AllocdUnmergedP1Seqs = cAllocOutBuffLen;
		m_CurUnmergedP1Seqs = 0;
		}

	if(m_hOut3Unmerged != -1 && m_pszUnmergedP2Seqs == NULL)
		{
		if((m_pszUnmergedP2Seqs = new char [cAllocOutBuffLen])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d) for buffering file writes",cAllocOutBuffLen);
			Reset(false);
			return(eBSFerrMem);
			}
		m_AllocdUnmergedP2Seqs = cAllocOutBuffLen;
		m_CurUnmergedP2Seqs = 0;
		}
	}

m_ProcPhase = ePRdyProc;

memset(OverlapDistCnts,0,sizeof(OverlapDistCnts));
memset(OverlapDistSubs,0,sizeof(OverlapDistSubs));
NumPE5Reads = 0;
NumPE3Reads = 0;
NumOverlapping = 0;
MinObsOverlap = -1;
MaxObsOverlap = -1;
NumUnmapped = 0;
NumPENoBarcode = 0;
NumPEWithBarcode = 0;
time_t Started = time(0);
while((Rslt = (teBSFrsltCodes)(PE5SeqLen = m_PE5Fasta.ReadSequence(PE5Seq,ccMaxFastQSeqLen,true,false))) > eBSFSuccess)
	{
	if(!(NumPE5Reads % 10000) && NumPE5Reads > 0)
		{
		time_t Now = time(0);
		unsigned long ElapsedSecs = (unsigned long) (Now - Started);
		if(ElapsedSecs >= 60)
			{
			switch(m_PMode) {
				case ePMAmplicon:
					gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processed %d paired reads, %d overlapped, %d reported with associated barcodes", NumPE5Reads, NumOverlapping, NumOverlapping - NumUnmapped);
					break;

				case ePMAmpliconNoMerge:
					gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processed %d paired reads, %d reported with associated barcodes", NumPE5Reads, NumPEWithBarcode);
					break;

				default:
					gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processed %d paired reads of which %d overlapped", NumPE5Reads, NumOverlapping);
					break;
				}

			Started = Now;

			if(m_PMode == ePMseparate)
				{
				if(m_CurUnmergedP1Seqs > 0)
					{
					CUtility::SafeWrite(m_hOut5Unmerged,m_pszUnmergedP1Seqs,(size_t)m_CurUnmergedP1Seqs);
					m_CurUnmergedP1Seqs = 0;
					}
				if(m_CurUnmergedP2Seqs > 0)
					{
					CUtility::SafeWrite(m_hOut3Unmerged,m_pszUnmergedP2Seqs,(size_t)m_CurUnmergedP2Seqs);
					m_CurUnmergedP2Seqs = 0;
					}
				}


			if(m_PMode >= ePMAmplicon && m_NumWells)
				{
				tsAmpliconWell *pWell;
				pWell = m_WellFiles;
				for(WellIdx = 0; WellIdx < m_NumWells; WellIdx++,pWell++)
					{
					if(pWell->WellFile[0].CurBuffLen > 0)
						{
						CUtility::SafeWrite(pWell->WellFile[0].hOutFile, pWell->WellFile[0].pOutBuffer,(size_t)pWell->WellFile[0].CurBuffLen);
						pWell->WellFile[0].CurBuffLen = 0;
						}
					if(pWell->WellFile[1].hOutFile != -1 && pWell->WellFile[1].CurBuffLen > 0)
						{
						CUtility::SafeWrite(pWell->WellFile[1].hOutFile, pWell->WellFile[1].pOutBuffer, (size_t)pWell->WellFile[1].CurBuffLen);
						pWell->WellFile[1].CurBuffLen = 0;
						}
					}
				}
			}
		}
 
	NumPE5Reads += 1;
	if(PE5SeqLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		PE5DescrLen = m_PE5Fasta.ReadDescriptor((char *)szPE5DescrBuff,sizeof(szPE5DescrBuff)-1);
		szPE5DescrBuff[sizeof(szPE5DescrBuff)-1] = '\0';

		PE5SeqLen = m_PE5Fasta.ReadSequence(PE5Seq,ccMaxFastQSeqLen);
		if(PE5SeqLen < cMinReadSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence from '%s' after %d reads parsed",m_szIn5ReadsFile,NumPE5Reads);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE5DescrBuff);
			m_PE5Fasta.Close();
			m_PE3Fasta.Close();
			return(eBSFerrParse);
			}
		if(m_bIsFastq)
			{
			PE5QualLen = m_PE5Fasta.ReadQValues((char *)PE5Qual,ccMaxFastQSeqLen);
			PE5Qual[PE5QualLen] = '\0';
			
			if(PE5QualLen != PE5SeqLen)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing quality from '%s' after %d reads parsed",m_szIn5ReadsFile,NumPE5Reads);
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Sequence length: %d Quality length: %d",PE5SeqLen,PE5QualLen);
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE5DescrBuff);
				m_PE5Fasta.Close();
				m_PE3Fasta.Close();
				return(eBSFerrParse);
				}
			}
		}

	Rslt = (teBSFrsltCodes)(PE3SeqLen = m_PE3Fasta.ReadSequence(PE3Seq,ccMaxFastQSeqLen,true,false));
	if(Rslt <= eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fewer reads in '%s' than '%s' after %d reads parsed",m_szIn3ReadsFile,m_szIn5ReadsFile,NumPE3Reads);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE3DescrBuff);
		m_PE5Fasta.Close();
		m_PE3Fasta.Close();
		return(eBSFerrParse);
		}	
	NumPE3Reads += 1;
	if(PE3SeqLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		PE3DescrLen = m_PE3Fasta.ReadDescriptor((char *)szPE3DescrBuff,sizeof(szPE3DescrBuff)-1);
		szPE3DescrBuff[sizeof(szPE3DescrBuff)-1] = '\0'; 
		PE3SeqLen = m_PE3Fasta.ReadSequence(PE3Seq,ccMaxFastQSeqLen);
		if(PE3SeqLen < cMinReadSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence from '%s' after %d reads parsed",m_szIn3ReadsFile,NumPE3Reads);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE3DescrBuff);
			m_PE5Fasta.Close();
			m_PE3Fasta.Close();
			return(eBSFerrParse);
			}
		if(m_bIsFastq)
			{
			PE3QualLen = m_PE3Fasta.ReadQValues((char *)PE3Qual,ccMaxFastQSeqLen);
			PE3Qual[PE3QualLen] = '\0';
			if(PE3QualLen != PE3SeqLen)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing quality from '%s' after %d reads parsed",m_szIn3ReadsFile,NumPE3Reads);
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Sequence length: %d Quality length: %d",PE3SeqLen,PE3QualLen);
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE5DescrBuff);
				m_PE5Fasta.Close();
				m_PE3Fasta.Close();
				return(eBSFerrParse);
				}
			}
		}

	  int PEPlateCell;
  	if(m_PMode >= ePMAmplicon && m_bAmpliconNoMerge)
		{
		if ((PEPlateCell = MapPEBarcodesToWell(min(PE5SeqLen,PE3SeqLen), PE5Seq, PE3Seq)) < 1)
			NumPENoBarcode += 1;
		else
			{
			tsAmpliconWell *pWell;
			MergeIdx = PE5SeqLen - 6;		// trim barcode
			CSeqTrans::MapSeq2UCAscii(&PE5Seq[6], MergeIdx, (char *)szMSeq);

			if (!m_bIsFastq && m_OFormat != eOFfasta)
				memset(PE5Qual, cPhredOrphanScore, PE5SeqLen);
			PE5Qual[PE5SeqLen] = '\0';
			szMSeq[MergeIdx] = '\0';
			pWell = &m_WellFiles[PEPlateCell - 1];
			if (m_OFormat == eOFfasta)
				pWell->WellFile[0].CurBuffLen += sprintf((char *)&pWell->WellFile[0].pOutBuffer[pWell->WellFile[0].CurBuffLen], ">Seq%d PE1 %s\n%s\n", StartNum + NumPEWithBarcode, szPE5DescrBuff, (char *)szMSeq);
			else
				pWell->WellFile[0].CurBuffLen += sprintf((char *)&pWell->WellFile[0].pOutBuffer[pWell->WellFile[0].CurBuffLen], "@Seq%d PE1 %s\n%s\n+\n%s\n", StartNum + NumPEWithBarcode, szPE5DescrBuff, (char *)szMSeq, &PE5Qual[6]);
				
			if (pWell->WellFile[0].CurBuffLen > pWell->WellFile[0].AllocdOutBuff - (4 * ccMaxFastQSeqLen))
				{
				CUtility::SafeWrite(pWell->WellFile[0].hOutFile, pWell->WellFile[0].pOutBuffer, (size_t)pWell->WellFile[0].CurBuffLen);
				pWell->WellFile[0].CurBuffLen = 0;
				}

			MergeIdx = PE3SeqLen - 6;		// trim 3' barcode
			CSeqTrans::MapSeq2UCAscii(&PE3Seq[6], MergeIdx, (char *)szMSeq);

			if (!m_bIsFastq && m_OFormat != eOFfasta)
				memset(PE3Qual, cPhredOrphanScore, PE3SeqLen);
			PE3Qual[PE3SeqLen] = '\0';
			szMSeq[MergeIdx] = '\0';
			if (m_OFormat == eOFfasta)
				pWell->WellFile[1].CurBuffLen += sprintf((char *)&pWell->WellFile[1].pOutBuffer[pWell->WellFile[1].CurBuffLen], ">Seq%d PE2 %s\n%s\n", StartNum + NumPEWithBarcode, szPE3DescrBuff, (char *)szMSeq);
			else
				pWell->WellFile[1].CurBuffLen += sprintf((char *)&pWell->WellFile[1].pOutBuffer[pWell->WellFile[1].CurBuffLen], "@Seq%d PE2 %s\n%s\n+\n%s\n", StartNum + NumPEWithBarcode, szPE3DescrBuff, (char *)szMSeq, &PE3Qual[6]);

			if (pWell->WellFile[1].CurBuffLen > pWell->WellFile[1].AllocdOutBuff - (4 * ccMaxFastQSeqLen))
				{
				CUtility::SafeWrite(pWell->WellFile[1].hOutFile, pWell->WellFile[1].pOutBuffer, (size_t)pWell->WellFile[1].CurBuffLen);
				pWell->WellFile[1].CurBuffLen = 0;
				}
			NumPEWithBarcode += 1;
			}

		continue;
		}

	// PE3 needs to be revcpl'd
	CSeqTrans::ReverseComplement(PE3SeqLen,PE3Seq);
	if(m_bIsFastq)
		CSeqTrans::ReverseSeq(PE3QualLen,PE3Qual);

	// now try for maximal overlap of at least m_MinOverlap allowing (if user specified) for sequencer base call errors
	// matches scored += 1, and mismatches scored -= 2 
	MaxSubs = m_MaxOverlapPropSubs; 
	MaxOverlap = 0;
	OLStartsIdx = 0;
	BestOvlpScore = -1;
	for(OL5Idx = m_PMode >= ePMAmplicon ? m_MaxBarcode5Len : 0; OL5Idx <= (PE5SeqLen - m_MinOverlap); OL5Idx++)
		{
		if(m_PMode >= ePMAmplicon && ((PE3SeqLen - m_MaxBarcode3Len) < (PE5SeqLen - OL5Idx)))
			continue;

		pSeq5 = &PE5Seq[OL5Idx];
		pSeq3 = PE3Seq;
		CurSubs = 0;
		CurOvlpScore = 0;
		ReqOverlap3 = min(PE5SeqLen - OL5Idx,PE3SeqLen);
		if(MaxSubs != 0)
			{
			if(ReqOverlap3 >= 20)
				AllowedSubs = min(MaxSubs,1 + ((ReqOverlap3 * m_MaxOverlapPropSubs) / 100));
			else
				{
				if(ReqOverlap3 >= 10)
					AllowedSubs = 2;
				else
					{
					if(ReqOverlap3 >= 5)
						AllowedSubs = 1;
					else
						AllowedSubs = 0;
					}
				AllowedSubs = min(MaxSubs,AllowedSubs);
				}
			}
		else
			AllowedSubs = 0;
		for(OL3Idx=0; OL3Idx < ReqOverlap3 && CurSubs <= AllowedSubs; OL3Idx++,pSeq5++,pSeq3++)
			{
			Base5 = *pSeq5 & 0x07;
			Base3 = *pSeq3 & 0x07;
			if(Base5 > eBaseT || Base3 > eBaseT || Base5 != Base3)
				{
				CurOvlpScore -= 2;
				CurSubs += 1;
				if(CurSubs > AllowedSubs)
					break;
				}
			else
				CurOvlpScore += 1;
			}

		if(OL3Idx == ReqOverlap3 && CurSubs <= MaxSubs)
			{
			if(CurOvlpScore > BestOvlpScore)
				{
				BestOvlpScore = CurOvlpScore;
				if(CurSubs == 0 || CurSubs < MaxSubs)
					{
					MaxOverlap = OL3Idx;
					OLStartsIdx = OL5Idx;
					}
				MaxSubs = CurSubs;
				if(MaxSubs == 0)
					break;
				}
			}
		}

	if(MaxOverlap < 1 && (m_PMode == ePMcombined || m_PMode == ePMseparate))
		{
		CSeqTrans::ReverseComplement(PE3SeqLen,PE3Seq);
		CSeqTrans::MapSeq2UCAscii(PE3Seq,PE3SeqLen,(char *)szMSeq);
		szMSeq[PE3SeqLen] = '\0';
		if(m_OFormat == eOFfastq)
			{
			if(m_bIsFastq)
				CSeqTrans::ReverseSeq(PE3QualLen,PE3Qual);
			else
				memset(PE3Qual,cPhredOrphanScore,PE3SeqLen);
			PE3Qual[PE3SeqLen] = '\0';
			}

		if(m_PMode == ePMseparate)
			{
			if(m_OFormat == eOFfasta)
				m_CurUnmergedP2Seqs += sprintf((char *)&m_pszUnmergedP2Seqs[m_CurUnmergedP2Seqs],">%s\n%s\n",szPE3DescrBuff,(char *)szMSeq);
			else
				m_CurUnmergedP2Seqs += sprintf((char *)&m_pszUnmergedP2Seqs[m_CurUnmergedP2Seqs],"@%s\n%s\n+\n%s\n",szPE3DescrBuff,(char *)szMSeq,(char *)PE3Qual);
				
			if(m_CurUnmergedP2Seqs > (cAllocOutBuffLen - (8 * ccMaxFastQSeqLen)))
				{
				CUtility::SafeWrite(m_hOut3Unmerged,m_pszUnmergedP2Seqs,(size_t)m_CurUnmergedP2Seqs);
				m_CurUnmergedP2Seqs = 0;
				}
			}
		else
			{
			if(m_OFormat == eOFfasta)
				m_CurMSeqLen += sprintf((char *)&m_pszMSeqs[m_CurMSeqLen],">%s\n%s\n",szPE3DescrBuff,(char *)szMSeq);
			else
				m_CurMSeqLen += sprintf((char *)&m_pszMSeqs[m_CurMSeqLen],"@%s\n%s\n+\n%s\n",szPE3DescrBuff,(char *)szMSeq,(char *)PE3Qual);

			if(m_CurMSeqLen > (cAllocOutBuffLen - (8 * ccMaxFastQSeqLen)))
				{
				CUtility::SafeWrite(m_hOutMerged,m_pszMSeqs,(size_t)m_CurMSeqLen);
				m_CurMSeqLen = 0;
				}	
			}

		CSeqTrans::MapSeq2UCAscii(PE5Seq,PE5SeqLen,(char *)szMSeq);
		szMSeq[PE5SeqLen] = '\0';
		if(m_OFormat == eOFfastq)
			{
			if(!m_bIsFastq)
				memset(PE5Qual,cPhredOrphanScore,PE5SeqLen);
			PE5Qual[PE5SeqLen] = '\0';
			}

			
		if(m_PMode == ePMseparate)
			{
			if(m_OFormat == eOFfasta)
				m_CurUnmergedP1Seqs += sprintf((char *)&m_pszUnmergedP1Seqs[m_CurUnmergedP1Seqs],">%s\n%s\n",szPE5DescrBuff,(char *)szMSeq);
			else
				m_CurUnmergedP1Seqs += sprintf((char *)&m_pszUnmergedP1Seqs[m_CurUnmergedP1Seqs],"@%s\n%s\n+\n%s\n",szPE5DescrBuff,(char *)szMSeq,(char *)PE5Qual);

			if(m_CurUnmergedP1Seqs > (cAllocOutBuffLen - (8 * ccMaxFastQSeqLen)))
				{
				CUtility::SafeWrite(m_hOut5Unmerged,m_pszUnmergedP1Seqs,(size_t)m_CurUnmergedP1Seqs);
				m_CurUnmergedP1Seqs = 0;
				}
			}
		else
			{
			if(m_OFormat == eOFfasta)
				m_CurMSeqLen += sprintf((char *)&m_pszMSeqs[m_CurMSeqLen],">%s\n%s\n",szPE5DescrBuff,(char *)szMSeq);
			else
				m_CurMSeqLen += sprintf((char *)&m_pszMSeqs[m_CurMSeqLen],"@%s\n%s\n+\n%s\n",szPE5DescrBuff,(char *)szMSeq,(char *)PE5Qual);

			if(m_CurMSeqLen > (cAllocOutBuffLen - (8 * ccMaxFastQSeqLen)))
				{
				CUtility::SafeWrite(m_hOutMerged,m_pszMSeqs,(size_t)m_CurMSeqLen);
				m_CurMSeqLen = 0;
				}	
			}

		continue;
		}

	if(MaxOverlap > 0)
		{
		if(MinObsOverlap == -1 || MaxOverlap < MinObsOverlap)
			MinObsOverlap = MaxOverlap;
		if(MaxObsOverlap == -1 || MaxOverlap > MaxObsOverlap)
			MaxObsOverlap = MaxOverlap;
		OverlapDistSubs[MaxSubs] += 1;
		OverlapDistCnts[MaxOverlap] += 1;
		NumOverlapping += 1;
		// output to file...
		pMSeq = MSeq;
		pSeq5 = PE5Seq;
		pMSeqQual = szMQual;
		pSeq5Qual = PE5Qual;
		for(MergeIdx = 0; MergeIdx < OLStartsIdx; MergeIdx++,pSeq5++,pMSeq++,pMSeqQual++,pSeq5Qual++)
			{
			*pMSeq = *pSeq5;
			if(m_OFormat == eOFfastq)
				{
				if(m_bIsFastq)
					*pMSeqQual = *pSeq5Qual;
				else
					*pMSeqQual = cPhredHiScore;
				}
			}
		
		pSeq3 = PE3Seq;
		pSeq3Qual = PE3Qual;

		if(MergeIdx < PE5SeqLen)
			{
			for(; MergeIdx < min(PE5SeqLen,PE3SeqLen + OLStartsIdx); MergeIdx++,pSeq5++,pMSeq++,pSeq3++,pMSeqQual++,pSeq5Qual++,pSeq3Qual++)
				{
				if(*pSeq5 == *pSeq3)
					{
					*pMSeq = *pSeq3;
					if(m_OFormat == eOFfastq)
						{
						if(m_bIsFastq)
							{
							if(*pSeq3Qual >= *pSeq5Qual)
								*pMSeqQual = *pSeq3Qual;
							else
								*pMSeqQual = *pSeq5Qual;
							}
						else
							*pMSeqQual = cPhredHiScore;
						}
					}
				else	// base difference: which one to choose... 
					{
					if(m_bIsFastq)	// choose base with highest Phred score
						{
						if(*pSeq3Qual >= *pSeq5Qual)
							{
							*pMSeq = *pSeq3;
							if(m_OFormat == eOFfastq)
								*pMSeqQual = *pSeq3Qual;
							}
						else
							{
							*pMSeq = *pSeq5;
							if(m_OFormat == eOFfastq)
								*pMSeqQual = *pSeq5Qual;
							}
						}
					else
						{
						*pMSeq = *pSeq3;
						if(m_OFormat == eOFfastq)
							*pMSeqQual = cPhredLowScore;
						}
					}
				}
			}
		


		if(MergeIdx < PE5SeqLen)
			{
			for( ;MergeIdx < PE5SeqLen;MergeIdx++,pMSeq++,pSeq5++,pMSeqQual++,pSeq5Qual++)
				{
				*pMSeq = *pSeq5;
				if(m_OFormat == eOFfastq)
					{
					if(m_bIsFastq)
						*pMSeqQual = *pSeq5Qual;
					else
						*pMSeqQual = cPhredHiScore;
					}
				}

			}
		else
			{
			if(MaxOverlap < PE3SeqLen)
				{
				for( ;MaxOverlap < PE3SeqLen;MaxOverlap++,MergeIdx++,pMSeq++,pSeq3++,pMSeqQual++,pSeq3Qual++)
					{
					*pMSeq = *pSeq3;
					if(m_OFormat == eOFfastq)
						{
						if(m_bIsFastq)
							*pMSeqQual = *pSeq3Qual;
						else
							*pMSeqQual = cPhredHiScore;
						}
					}
				}
			}

		if(m_PMode == ePMAmplicon)
			{
			if((PlateCell = MapSEBarcodesToWell(MergeIdx,MSeq)) < 1)
				{
				NumUnmapped += 1;
				continue;
				}
			}

		CSeqTrans::MapSeq2UCAscii(MSeq,MergeIdx,(char *)szMSeq);
		szMSeq[MergeIdx] = '\0';

		if(m_PMode == ePMAmplicon)
			{
			tsAmpliconWell *pWell;
			MergeIdx -= 6;				// trim 3' barcode
			szMQual[MergeIdx] = '\0';
			szMSeq[MergeIdx] = '\0';	
			pWell = &m_WellFiles[PlateCell-1];
			if(m_OFormat == eOFfasta)
				pWell->WellFile[0].CurBuffLen += sprintf((char *)&pWell->WellFile[0].pOutBuffer[pWell->WellFile[0].CurBuffLen],">MSeq%d %s\n%s\n",StartNum+NumOverlapping,szPE5DescrBuff,(char *)&szMSeq[6]);
			else
				{
				szMQual[MergeIdx] = '\0';
				pWell->WellFile[0].CurBuffLen += sprintf((char *)&pWell->WellFile[0].pOutBuffer[pWell->WellFile[0].CurBuffLen],"@MSeq%d %s\n%s\n+\n%s\n",StartNum+NumOverlapping,szPE5DescrBuff,(char *)&szMSeq[6],&szMQual[6]);
				}
			if(pWell->WellFile[0].CurBuffLen > pWell->WellFile[0].AllocdOutBuff - (4 * ccMaxFastQSeqLen))
				{
				CUtility::SafeWrite(pWell->WellFile[0].hOutFile, pWell->WellFile[0].pOutBuffer,(size_t)pWell->WellFile[0].CurBuffLen);
				pWell->WellFile[0].CurBuffLen = 0;
				}
			}
		else
			{
			if(m_OFormat == eOFfasta)
				m_CurMSeqLen += sprintf((char *)&m_pszMSeqs[m_CurMSeqLen],">MSeq%d %s\n%s\n",StartNum+NumOverlapping,szPE5DescrBuff,(char *)szMSeq);
			else
				{
				szMQual[MergeIdx] = '\0';
				m_CurMSeqLen += sprintf((char *)&m_pszMSeqs[m_CurMSeqLen],"@MSeq%d %s\n%s\n+\n%s\n",StartNum+NumOverlapping,szPE5DescrBuff,(char *)szMSeq,szMQual);
				}
			if(m_CurMSeqLen > (cAllocOutBuffLen - (8 * ccMaxFastQSeqLen)))
				{
				CUtility::SafeWrite(m_hOutMerged,m_pszMSeqs,(size_t)m_CurMSeqLen);
				m_CurMSeqLen = 0;
				}
			}
		}
	}

if(m_PMode == ePMAmplicon && m_NumWells)
	{
	tsAmpliconWell *pWell;
	pWell = m_WellFiles;
	for(WellIdx = 0; WellIdx < m_NumWells; WellIdx++,pWell++)
		{
		if(pWell->WellFile[0].CurBuffLen > 0)
			{
			CUtility::SafeWrite(pWell->WellFile[0].hOutFile, pWell->WellFile[0].pOutBuffer,(size_t)pWell->WellFile[0].CurBuffLen);
			pWell->WellFile[0].CurBuffLen = 0;
			}

		if (pWell->WellFile[1].hOutFile != -1 && pWell->WellFile[1].CurBuffLen > 0)
			{
			CUtility::SafeWrite(pWell->WellFile[1].hOutFile, pWell->WellFile[1].pOutBuffer, (size_t)pWell->WellFile[1].CurBuffLen);
			pWell->WellFile[1].CurBuffLen = 0;
			}
		}
	}

if(m_CurMSeqLen > 0)
	{
	CUtility::SafeWrite(m_hOutMerged,m_pszMSeqs,(size_t)m_CurMSeqLen);
	m_CurMSeqLen = 0;
	}

if(m_CurUnmergedP1Seqs > 0)
	{
	CUtility::SafeWrite(m_hOut5Unmerged,m_pszUnmergedP1Seqs,(size_t)m_CurUnmergedP1Seqs);
	m_CurUnmergedP1Seqs = 0;
	}

if(m_CurUnmergedP2Seqs > 0)
	{
	CUtility::SafeWrite(m_hOut3Unmerged,m_pszUnmergedP2Seqs,(size_t)m_CurUnmergedP2Seqs);
	m_CurUnmergedP2Seqs = 0;
	}

switch (m_PMode)
	{
	case ePMAmplicon:
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processed %d paired reads, %d overlapped, %d reported with associated barcodes", NumPE5Reads, NumOverlapping, NumOverlapping - NumUnmapped);
		break;

	case ePMAmpliconNoMerge:
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processed %d paired reads, %d reported with associated barcodes", NumPE5Reads, NumPEWithBarcode);
		return(NumPEWithBarcode);

	default:
		gDiagnostics.DiagOut(eDLInfo, gszProcName, "Processed %d paired reads of which %d overlapped", NumPE5Reads, NumOverlapping);
		break;
	}

if(NumOverlapping)
	{
	for(Idx = MinObsOverlap; Idx <= MaxObsOverlap; Idx++)
		gDiagnostics.DiagOut(eDLDebug,gszProcName,"OverlapLen: %d Counts: %d",Idx,OverlapDistCnts[Idx]);
	for(Idx = 0; Idx <= cMaxOverlapPercSubs; Idx++)
		gDiagnostics.DiagOut(eDLDebug,gszProcName,"Substitutions: %d Counts: %d",Idx,OverlapDistSubs[Idx]);
	}

return(NumOverlapping);
}

