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


CContaminants::CContaminants(void)
{
m_pContaminants = NULL;
m_pSeqBuff = NULL;
m_TotNumContaminants = 0;
m_NumFlankContaminates = 0;
m_NumVectContaminates = 0;
memset(m_ContaminantVectors,0,sizeof(m_ContaminantVectors));
m_AllocContaminantsMem = 0;
m_pContamSeqNodes = 0;
m_AllocdContamSeqNodesMem = 0;	
m_bSerialCreated = false;
Reset();
}


CContaminants::~CContaminants(void)
{
if(m_pSeqBuff != NULL)
	delete m_pSeqBuff;

if (m_pContaminants != NULL)
	{
#ifdef _WIN32
	free(m_pContaminants);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pContaminants != MAP_FAILED)
		munmap(m_pContaminants, m_AllocContaminantsMem);
#endif
	}

if (m_pContamSeqNodes != NULL)
	{
#ifdef _WIN32
	free(m_pContamSeqNodes);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pContamSeqNodes != MAP_FAILED)
		munmap(m_pContamSeqNodes, m_AllocdContamSeqNodesMem);
#endif
	}
if(m_NumVectContaminates > 0)
	{
	for(int Idx = 0; Idx < m_NumVectContaminates; Idx++)
		{
		if(m_ContaminantVectors[Idx].pBases != NULL)
			delete m_ContaminantVectors[Idx].pBases;
		if(m_ContaminantVectors[Idx].pSfxIdx != NULL)
			delete m_ContaminantVectors[Idx].pSfxIdx;
		}
	}
}

void
CContaminants::AcquireSerialise(void)
{
int SpinCnt = 1000;
#ifdef _WIN32
while(!TryEnterCriticalSection(&m_hSCritSect))
	{
	if (SpinCnt -= 1)
		continue;
	SwitchToThread();
	SpinCnt = 100;
	}
#else
while (pthread_spin_trylock(&m_hSpinLock) == EBUSY)
	{
	if (SpinCnt -= 1)
		continue;
	pthread_yield();
	SpinCnt = 100;
	}
#endif
}

void
CContaminants::ReleaseSerialise(void)
{
#ifdef _WIN32
LeaveCriticalSection(&m_hSCritSect);
#else
pthread_spin_unlock(&m_hSpinLock);
#endif
}

int
CContaminants::Init(void)
{
m_NumFlankContaminates = 0;
m_TotNumContaminants = 0;
m_AllocContaminants = 0;
m_NumContamSeqNodes = 0;			
szContaminantFile[0] = '\0';

m_MaxFlankContamSeqLen = 0;
m_MinFlankContamSeqLen = 0;
m_MaxVectContamSeqLen = 0;
m_MinVectContamSeqLen = 0;

m_CacheContamID = 0;
m_CacheContamIDClass = eCCAllContam;
m_CacheContamIDIdx = 0;

memset(m_ContaminantTypes,0,sizeof(m_ContaminantTypes));

m_NumVectContaminates = 0;
memset(m_ContaminantVectors,0,sizeof(m_ContaminantVectors));

if(!m_bSerialCreated)
	{
#ifdef _WIN32
	if (!InitializeCriticalSectionAndSpinCount(&m_hSCritSect, 1000))
		{
#else
	if(pthread_spin_init(&m_hSpinLock,PTHREAD_PROCESS_PRIVATE)!=0)
		{
#endif
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "Fatal: unable to create serialisation locks");
		return(eBSFerrInternal);
		}
	m_bSerialCreated = true;
	}
return(eBSFSuccess);
}

void
CContaminants::Reset(void)		// reset back to state immediately following class instantiation
{
int Idx;
if (m_pContaminants != NULL)
	{
#ifdef _WIN32
	free(m_pContaminants);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pContaminants != MAP_FAILED)
		munmap(m_pContaminants, m_AllocContaminantsMem);
#endif
	m_pContaminants = NULL;
	m_AllocContaminantsMem = 0;
	}

if (m_pContamSeqNodes != NULL)
	{
#ifdef _WIN32
	free(m_pContamSeqNodes);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if (m_pContamSeqNodes != MAP_FAILED)
		munmap(m_pContamSeqNodes, m_AllocdContamSeqNodesMem);
#endif
	m_pContamSeqNodes = NULL;
	m_AllocdContamSeqNodesMem = 0;
	}

if(m_NumVectContaminates > 0)
	{
	for(Idx = 0; Idx < m_NumVectContaminates; Idx++)
		{
		if(m_ContaminantVectors[Idx].pBases != NULL)
			delete m_ContaminantVectors[Idx].pBases;
		if(m_ContaminantVectors[Idx].pSfxIdx != NULL)
			delete m_ContaminantVectors[Idx].pSfxIdx;
		}
	}
m_NumVectContaminates = 0;
memset(m_ContaminantVectors,0,sizeof(m_ContaminantVectors));
m_NumFlankContaminates = 0;
m_TotNumContaminants = 0;
Init();
}


// Sequence naming convention
// It is expected that the sequence names are suffixed with a contaminant overlay code and if not present then defaults to that sequence applying to 5' PE1 and PE2 putative overlays
// For adaptor (or ends of vectors) type processing the suffix starts with the '@' separator char followed by any combination of '1' (5PE1), '2' (5PE2), '3' (3PE1), '4' (3PE2) '5' (5PE1 RevCpl), '6' (5PE2 RevCpl), '7' (3PE1 RevCpl) and '8' (3PE2 RevCpl)
// For vector type processing then the suffix starts with the '&' separator char following by any combination of '1' (PE1), '2' (PE2), '5' (PE1 RevCpl), '6' (PE2 RevCpl) 
// Thus if the sequence name is:
// >contamABC@12
// Then the forgoing specifies that the contaminant sequence contamABC is used when checking for overlays onto the 5' end of both PE1 and PE2 target sequences
// Sequence ordering in the contaminant file determines the priority of that sequence, later occuring sequences have lower priority than earlier occuring sequences
// 
int									// < 0 if errors, 0..N the number of contaminant sequences loaded
CContaminants::LoadContaminantsFile(char *pszContaminantFile)		// load contaminant sequences from this multifasta file
{
int Rslt;
CFasta Fasta;
UINT8 *pSeqBase;
int SeqIdx;
char szName[cMaxGeneNameLen + 50];				// additional is to allow for user name concatened contaminant codes
int Idx;
char *pChar;
char Chr;
char szDescription[cMaxGeneNameLen + 25];
int NameLen;
int SeqLen;
int Descrlen;
int SeqID;

bool b5PE1;				// (1) true if Contaminant valid for overlaps onto 5' SE/PE1 sequences
bool b5PE2;				// (2) true if Contaminant valid for overlaps onto 5' PE2 sequences
bool b3PE1;				// (3) true if Contaminant valid for overlaps onto 3' SE/PE1 sequences
bool b3PE2;				// (4) true if Contaminant valid for overlaps onto 3' PE2 sequences
bool b5PE1RC;			// (5) true if Contaminant after RevCpl is valid for overlaps onto 5' SE/PE1 sequences
bool b5PE2RC;			// (6) true if Contaminant after RevCpl is valid for overlaps onto 5' PE2 sequences
bool b3PE1RC;			// (7) true if Contaminant after RevCpl is valid for overlaps onto 3' SE/PE1 sequences
bool b3PE2RC;			// (8) true if Contaminant after RevCpl is valid for overlaps onto 3' PE2 sequences
bool bIsVector;			

if((m_pSeqBuff = new UINT8 [cMaxVectorSeqLen+1]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadContaminantsFile: Unable to allocate %d bytes memory",cMaxVectorSeqLen);
	Reset();
	return(eBSFerrMem);
	}

if ((Rslt = Fasta.Open(pszContaminantFile, true)) != eBSFSuccess)
	{
	if (Rslt != eBSFerrNotFasta)
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadContaminantsFile: Unable to open '%s' [%s] %s", pszContaminantFile, Fasta.ErrText((teBSFrsltCodes)Rslt), Fasta.GetErrMsg());
	while(Fasta.NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,Fasta.GetErrMsg());
	Reset();
	return(Rslt);
	}

gDiagnostics.DiagOut(eDLInfo, gszProcName, "LoadContaminantsFile:- Processing %s..", pszContaminantFile);

SeqID = 0;
b5PE1 = false;
b5PE2 = false;
b3PE1 = false;
b3PE2 = false;
bIsVector = false;
while ((Rslt = SeqLen = Fasta.ReadSequence(m_pSeqBuff, cMaxVectorSeqLen)) > eBSFSuccess)
	{
	if(SeqLen == eBSFFastaDescr)		// just read a descriptor line so assume about to start a new contaminant sequence
		{
		b5PE1 = false;
		b5PE2 = false;
		b3PE1 = false;
		b3PE2 = false;
		b5PE1RC = false;
		b5PE2RC = false;
		b3PE1RC = false;
		b3PE2RC = false;
		SeqID++;
		Descrlen = Fasta.ReadDescriptor(szDescription, sizeof(szDescription)-1);
		if (sscanf(szDescription, " %s[ ,]", szName) != 1)
			{
			sprintf(szName, "ContamSeq.%d", SeqID);
			b5PE1 = true;
			b5PE2 = true;
			b5PE1RC = true;
			b5PE2RC = true;
			}
		else
			{
			NameLen = (int)strlen(szName);		
			pChar = &szName[NameLen-1];
			for(Idx = NameLen; Idx > 1; Idx--,pChar--)	
				{
				if(*pChar == '@' || *pChar == '&' || !(*pChar >= '1' && *pChar <= '8'))
					break;
				}
			bIsVector = *pChar == '&' ? true : false;
			if((*pChar == '@' || *pChar == '&') && pChar[1] != '\0')
				{
				*pChar = '\0';
				pChar += 1;
				while(Chr = *pChar++)
					{
					switch(Chr) {
						case '1':
							b5PE1 = true;
							break;
						case '2':
							b5PE2 = true;
							break;
						case '3':
							b3PE1 = true;
							break;
						case '4':
							b3PE2 = true;
							break;
						case '5':
							b5PE1RC = true;
							break;
						case '6':
							b5PE2RC = true;
							break;
						case '7':
							b3PE1RC = true;
							break;
						case '8':
							b3PE2RC = true;
							break;
						}
					}
				}
			else
				{
				b5PE1 = true;
				b5PE2 = true;
				b5PE1RC = true;
				b5PE2RC = true;
				}
			}
		szName[cMaxGeneNameLen-1] = '\0';
		continue;
		}

	if(!SeqID)	// if there was no descriptor then dummy up one...
		{
		SeqID++;
		sprintf(szName, "ContamSeq.%d", SeqID);
		strcpy(szDescription, "No Description provided");
		b5PE1 = true;
		b5PE2 = true;
		b5PE1RC = true;
		b5PE2RC = true;
		b3PE1 = false;
		b3PE2 = false;
		b3PE1RC = false;
		b3PE2RC = false;
		}

	// check contaminate sequence is within accepted length range for type - adaptor or vector
	if(!bIsVector && (SeqLen < cMinContaminantLen || SeqLen > cMaxContaminantLen))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadContamiantsFile: Sequence for '%s' outside of accepted length range %d..%d", szName, cMinContaminantLen, cMaxContaminantLen);
		Rslt = Rslt = eBSFerrFastqSeq;
		break;
		}
	if(bIsVector && (SeqLen < cMinVectorSeqLen || SeqLen > cMaxVectorSeqLen))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadContamiantsFile: Vector sequence for '%s' outside of accepted length range %d..%d", szName, cMinContaminantLen, cMaxContaminantLen);
		Rslt = Rslt = eBSFerrFastqSeq;
		break;
		}

	// only accepting a,c,g,t,n bases in contaminate sequences so check for these
	pSeqBase = m_pSeqBuff;
	for (SeqIdx = 0; SeqIdx < SeqLen; SeqIdx++, pSeqBase++)
		{
		*pSeqBase &= 0x07;		
		if(*pSeqBase > eBaseN)
			{
			gDiagnostics.DiagOut(eDLFatal, gszProcName, "LoadContaminantsFile: Illegal base in % sequence, only bases A,C,G,T,N accepted",szName);
			Rslt = eBSFErrBase;
			break;
			}
		}

	if(bIsVector)
		{
		if ((Rslt = AddVectContam(b5PE1 || b3PE1,b5PE1RC || b3PE1RC,b5PE2 || b3PE2, b5PE2RC || b3PE2RC,szName, SeqLen,m_pSeqBuff)) < eBSFSuccess)
			break;
		}
	else
		{
		if(b5PE1)
			if ((Rslt = AddFlankContam(eAOF5PE1Targ,false,szName, SeqLen,m_pSeqBuff)) < eBSFSuccess)
				break;

		if(b5PE2)
			if ((Rslt = AddFlankContam(eAOF5PE2Targ,false,szName, SeqLen,m_pSeqBuff)) < eBSFSuccess)
				break;

		if(b3PE1)
			if ((Rslt = AddFlankContam(eAOF3PE1Targ,false,szName, SeqLen,m_pSeqBuff)) < eBSFSuccess)
				break;

		if(b3PE2)
			if ((Rslt = AddFlankContam(eAOF3PE2Targ,false,szName, SeqLen,m_pSeqBuff)) < eBSFSuccess)
				break;

		if(b5PE1RC || b5PE2RC || b3PE1RC || b3PE2RC)		// allowing RevCpl contaminants overlaps onto targets??
			{
			CSeqTrans::ReverseComplement(SeqLen,m_pSeqBuff);
			strcat(szName,"xRC");
			if(b5PE1RC)
				if ((Rslt = AddFlankContam(eAOF5PE1Targ,true,szName, SeqLen,m_pSeqBuff)) < eBSFSuccess)
					break;

			if(b5PE2RC)
				if ((Rslt = AddFlankContam(eAOF5PE2Targ,true,szName, SeqLen,m_pSeqBuff)) < eBSFSuccess)
					break;

			if(b3PE1RC)
				if ((Rslt = AddFlankContam(eAOF3PE1Targ,true,szName, SeqLen,m_pSeqBuff)) < eBSFSuccess)
					break;

			if(b3PE2RC)
				if ((Rslt = AddFlankContam(eAOF3PE2Targ,true,szName, SeqLen,m_pSeqBuff)) < eBSFSuccess)
					break;
			}
		}
	}

if(Rslt >= eBSFSuccess)
	Rslt = FinaliseContaminants();
if(Rslt < eBSFSuccess)
	Reset();
return(Rslt);
}


// When last Contaminant has been added with AddContaminant then Finalise() must be called to
// sort Contaminants by length descending and for the initialisation of m_ContaminantTypes[]
int
CContaminants::FinaliseContaminants(void)
{
int ContamIdx;

tsFlankContam *pPrevContaminants[4];

tsFlankContam *pContaminant;
tsContaminantType *pContaminantType;

if(!m_TotNumContaminants)					// can't check for Contaminants if none loaded!
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "FinaliseContaminant: No Contaminants loaded");
	Reset();
	return(eBSFerrParams);
	}

if(m_NumFlankContaminates)
	{
	if(m_NumFlankContaminates > 1)
		{
		m_mtqsort.SetMaxThreads(2);
		m_mtqsort.qsort(m_pContaminants,m_NumFlankContaminates,sizeof(tsFlankContam),SortContamTypeLen);
		}

	memset(m_ContaminantTypes,0,sizeof(m_ContaminantTypes));
	memset(pPrevContaminants,0,sizeof(pPrevContaminants));

	pContaminant = m_pContaminants;
	for(ContamIdx=0; ContamIdx < m_NumFlankContaminates; ContamIdx++,pContaminant++)
		{
		pContaminant->ContamID = ContamIdx+1; 
		pContaminantType = &m_ContaminantTypes[pContaminant->Type];	
		if(pContaminantType->pFirstContam == NULL)
			pContaminantType->pFirstContam = pContaminant;
		pContaminantType->pLastContam = pContaminant;
		pContaminantType->NumContaminants += 1;
		if(pContaminantType->MaxContamSeqLen < pContaminant->ContamLen)
			pContaminantType->MaxContamSeqLen = pContaminant->ContamLen;
		if(m_MaxFlankContamSeqLen < pContaminant->ContamLen)
			m_MaxFlankContamSeqLen = pContaminant->ContamLen;
		if(pContaminantType->MinContamSeqLen == 0 || pContaminantType->MinContamSeqLen >= pContaminant->ContamLen)
			pContaminantType->MinContamSeqLen = pContaminant->ContamLen;
		if(m_MinFlankContamSeqLen == 0 || m_MinFlankContamSeqLen > pContaminant->ContamLen)
			m_MinFlankContamSeqLen = pContaminant->ContamLen;
		}

	// now generate the b-tree like index over the contaminate sequences
	pContaminant = m_pContaminants;
	for(ContamIdx=0; ContamIdx < m_NumFlankContaminates; ContamIdx++,pContaminant++)
		{
		IndexContamSeq(pContaminant->Type,			// contaminant sequence is of this overlay type
						 pContaminant->ContamID,	// identifies contaminant sequence
						 pContaminant->ContamLen,	// contaminant sequence is of this length
						 pContaminant->Bases);		// contaminant sequence
		}
	}
return(m_TotNumContaminants);
}

char *									// returned text for Type
CContaminants::ContaminateType2Txt(teContamType Type)	// contaminant is of this overlay type)
{
char *pTxt;
switch(Type) {
	case eAOF5PE1Targ:		// match for Contaminant sequences valid for overlaying onto 5' of a SE/PE1 target sequence
		pTxt = (char *)"5'SE/PE1";
		break;
	case eAOF5PE2Targ:		// match for Contaminant sequences valid for overlaying onto 5' of a PE2 target sequence
		pTxt = (char *)"5'PE2";
		break;
	case eAOF3PE1Targ:		// match for Contaminant sequences valid for overlaying onto 3' of a SE/PE1 target sequence
		pTxt = (char *)"3'SE/PE1";
		break;
	case eAOF3PE2Targ:		// match for Contaminant sequences valid for overlaying onto 3' of a PE2 target sequence
		pTxt = (char *)"3'PE2";
		break;
	case eAOFVector:		// match for contaminant sequence completely containing the target sequence
		pTxt = (char *)"Vector";
		break;
	default:
		pTxt = NULL;
		break;	
	}
return(pTxt);
}

int												// returned contaminate identifier
CContaminants::AddVectContam(bool bPE1Sense,// check for PE1 contained sense
					bool bPE1Antisense,		// check for PE1 contained antisense
					bool bPE2Sense,			// check for PE2 contained sense
					bool bPE2Antisense,		// check for PE2 contained antisense
					char *pszName,			// vector name
					int ContamLen,			// vector sequence is of this length
					etSeqBase *pContamSeq)	// the vector sequence
{
etSeqBase ContamBase;
etSeqBase  *pContamBase;
tsVectContam *pVect;
INT32 *pSfxIdx;
int VectIdx;
int Idx;

if(!(bPE1Sense || bPE1Antisense || bPE2Sense || bPE2Antisense) || pszName == NULL || pszName[0] == '\0' || ContamLen < cMinVectorSeqLen || ContamLen > cMaxVectorSeqLen)
	return(-1);

if(m_NumVectContaminates >= cMaxNumVectors)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddVectContam: Too many vector contaminants (max allowed %d) current contaminant is '%s'",cMaxNumVectors,pszName);
	Reset();
	return(eBSFerrParams);
	}

if(strlen(pszName) > cMaxGeneNameLen-1)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddVectContam: Contaminant name '%s' length must be <= %d ",pszName,cMaxGeneNameLen-1);
	Reset();
	return(eBSFerrFastqSeqID);
	}

if(ContamLen < cMinVectorSeqLen || ContamLen > cMaxVectorSeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddVectContam: Contaminant '%s' sequence length expected to be in range %d..%d ",pszName,cMinVectorSeqLen,cMaxVectorSeqLen);
	Reset();
	return(eBSFerrFastqSeq);
	}

if(ContaminantID(pszName) > 0)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddVectContam: contaminate vector name '%s' already in use",pszName);
	Reset();
	return(eBSFerrFastaDescr);
	}

pContamBase = pContamSeq;
for(Idx = 0; Idx < ContamLen; Idx++,pContamBase++)
	{
	ContamBase = *pContamBase;
	if(ContamBase > eBaseN)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddVectContam: Contaminant name '%s' sequence contains illegal Contaminant base",pszName);
		Reset();
		return(eBSFerrFastqSeq);
		}
	}

pVect = &m_ContaminantVectors[0];
for(VectIdx = 0; VectIdx < m_NumVectContaminates; VectIdx++,pVect++)
	{
	if(!stricmp(pVect->szName,pszName))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddVectContam: contaminate vector name '%s' already in use",pszName);
		Reset();
		return(eBSFerrMem);
		}
	if((ContamLen == pVect->ContamLen) && !memcmp(pVect->pBases,pContamSeq,ContamLen))
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddVectContam: duplicate vector sequence for '%s' and '%s'",pszName,pVect->szName);
		Reset();
		return(eBSFerrMem);
		}
	}

// accepting this vector sequence
pVect = &m_ContaminantVectors[m_NumVectContaminates++];
memset(pVect,0,sizeof(*pVect));
pVect->ContamID = ++m_TotNumContaminants;
strcpy(pVect->szName,pszName);
pVect->FlgPE1Sense = bPE1Sense;
pVect->FlgPE1Antisense = bPE1Antisense;
pVect->FlgPE2Sense = bPE2Sense;
pVect->FlgPE2Antisense = bPE2Antisense;
pVect->HitTot = 0;
pVect->ContamLen = ContamLen;
if((pVect->pBases = new UINT8 [ContamLen+1])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddVectContam: unable to allocate for vector sequence memory %d", ContamLen+1);
	Reset();
	return(eBSFerrMem);
	}

memcpy(pVect->pBases,pContamSeq,ContamLen);
pVect->pBases[ContamLen] = eBaseEOS;
if((pVect->pSfxIdx = new int [ContamLen+1])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddVectContam: unable to allocate for vector index memory %d", sizeof(int) * (ContamLen+1));
	Reset();
	return(eBSFerrMem);
	}
pSfxIdx = pVect->pSfxIdx;
for(VectIdx=0;VectIdx <= ContamLen; VectIdx++)
	*pSfxIdx++ = VectIdx;
m_pCurVector2Index = pVect;
m_mtqsort.SetMaxThreads(8);
m_mtqsort.qsort(pVect->pSfxIdx,ContamLen,sizeof(UINT32),IndexVectorSeq);

return(pVect->ContamID);
}


int												// returned contaminate identifier
CContaminants::AddFlankContam(teContamType Type,	// contaminant is of this overlay type
					bool bRevCpl,			// contaminant sequence has been RevCpl'd
					char *pszName,			// Contaminant name
					int ContamLen,			// Contaminant sequence is of this length
					etSeqBase *pContamSeq)	// the Contaminant sequence
{
size_t AllocMem;
char *pszType;
int Idx;
etSeqBase ContamBase;
etSeqBase  *pContamBase;
tsFlankContam *pFlankContam;

// validate the parameters
if (pszName == NULL || pszName[0] == '\0' || pContamSeq == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFlankContam: Parameter errors");
	Reset();
	return(eBSFerrParams);
	}

if(m_NumFlankContaminates >= cMaxNumContaminants)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFlankContam: Too many flank contaminants (max allowed %d) current contaminant is '%s'",cMaxNumContaminants,pszName);
	Reset();
	return(eBSFerrParams);
	}

if((pszType = ContaminateType2Txt(Type)) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFlankContam: Parameter error, unrecognised type %d",Type);
	Reset();
	return(eBSFerrParams);
	}

if(strlen(pszName) > cMaxGeneNameLen-1)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFlankContam: Contaminant name '%s' of overlap type %s length must be <= %d ",pszName,pszType,cMaxGeneNameLen-1);
	Reset();
	return(eBSFerrFastqSeqID);
	}
if(ContamLen < cMinContaminantLen || ContamLen > cMaxContaminantLen)
	{
	gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFlankContam: Contaminant '%s' of overlap type %s sequence length expected to be in range %d..%d ",pszName,pszType,cMinContaminantLen,cMaxContaminantLen);
	Reset();
	return(eBSFerrFastqSeq);
	}

pContamBase = pContamSeq;
for(Idx = 0; Idx < ContamLen; Idx++,pContamBase++)
	{
	ContamBase = *pContamBase;
	if(ContamBase > eBaseN)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFlankContam: Contaminant name '%s' of overlap type %s contains illegal Contaminant base",pszName,pszType);
		Reset();
		return(eBSFerrFastqSeq);
		}
	}

// if first Contaminant for any type added then allocate to hold initial cAllocNumContaminants Contaminants
if(m_pContaminants == NULL)
	{
	// initial allocation for the tsContaminant's; will be realloc'd later as may be required
	AllocMem = (size_t)cAllocNumContaminants * sizeof(tsFlankContam);
#ifdef _WIN32
	m_pContaminants = (tsFlankContam *)malloc(AllocMem);
	if (m_pContaminants == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFlankContam: unable to allocate for Contaminant memory %llu", AllocMem);
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pContaminants = (tsFlankContam *)mmap(NULL, (size_t)AllocMem, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pContaminants == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFlankContam: unable to allocate for Contaminant memory %llu", AllocMem);
		m_pContaminants = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocContaminantsMem = AllocMem;
	m_AllocContaminants = cAllocNumContaminants;
	m_NumFlankContaminates = 0;
	memset(m_pContaminants,0,AllocMem);
	}

if((m_NumFlankContaminates + 1) > m_AllocContaminants)		// time to realloc some more Contaminants? 
	{
	size_t ReallocSize;
	tsFlankContam *pRealloc;
	ReallocSize = (m_AllocContaminants + cAllocNumContaminants) * sizeof(tsFlankContam);
#ifdef _WIN32
	pRealloc = (tsFlankContam *)realloc(m_pContaminants, ReallocSize);
#else
	pRealloc = (tsFlankContam *)mremap(m_pContaminants, m_AllocContaminantsMem, ReallocSize, MREMAP_MAYMOVE);
	if (pRealloc == MAP_FAILED)
		pRealloc = NULL;
#endif
	if (pRealloc == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFlankContam: Memory re-allocation to %lld bytes - %s", ReallocSize, strerror(errno));
		return(eBSFerrMem);
		}

	m_pContaminants = pRealloc;
	memset(&m_pContaminants[m_AllocContaminants],0, ReallocSize - m_AllocContaminantsMem);
	m_AllocContaminantsMem = ReallocSize;
	m_AllocContaminants += cAllocNumContaminants;
	}

// do quick check to see if this sequence is a duplicate in name or sequence of a previously loaded Contaminant of same type
// if Contaminant with same name already added, but different sequences then error
// if another Contaminant exists, different name, with same sequence then error
// otherwise create a new Contaminant
pFlankContam = m_pContaminants;
for(Idx=0; Idx < m_NumFlankContaminates; Idx++,pFlankContam++)
	{
	if(pFlankContam->Type != Type)
		continue;
	int bSameSeq;
	int bSameName;
	if(pFlankContam->ContamLen == ContamLen)
		bSameSeq = memcmp(pFlankContam->Bases,pContamSeq,ContamLen) == 0 ? true : false;
	else
		bSameSeq = false;
	 bSameName = stricmp(pFlankContam->szName,pszName)  == 0 ? true : false;

	 if(!bSameName && !bSameSeq)
		continue;

	 if(bSameName && bSameSeq)				// if contaminant of same name and sequence then can't accept duplicates
	 	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFlankContam: Contaminant name '%s' of overlap type %s duplicated with same sequences",pszName,pszType);
		Reset();
		return(eBSFerrFastqSeqID);
		}

	 if(bSameName && !bSameSeq)
	 	{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFlankContam: Contaminant name '%s' of overlap type %s duplicated with different sequences",pszName,pszType);
		Reset();
		return(eBSFerrFastqSeqID);
		}
	
	if(!bSameName && bSameSeq)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "AddFlankContam: Contaminant names '%s' and '%s' of overlap type %s with duplicated sequence",pszName,pFlankContam->szName,pszType);
		Reset();
		return(eBSFerrFastqSeqID);
		}
	}

// accepting this Contaminant
pFlankContam = &m_pContaminants[m_NumFlankContaminates++];
memset(pFlankContam,0,sizeof(tsFlankContam));
pFlankContam->ContamID = ++m_TotNumContaminants;
pFlankContam->Type = Type;
pFlankContam->FlgRevCpl = bRevCpl ? 1 : 0;
pFlankContam->ContamLen = ContamLen;
strcpy(pFlankContam->szName,pszName);
memcpy(pFlankContam->Bases,pContamSeq,ContamLen);
pFlankContam->Bases[ContamLen] = eBaseEOS;
return(pFlankContam->ContamID);
}


int 
CContaminants::IndexContamSeq(teContamType Type,	// contaminant sequence is of this overlay type
					 UINT16 ContamID,		// identifies contaminant sequence
					 UINT8 SeqLen,			// contaminant sequence is of this length
					 etSeqBase *pSeq)		// contaminant sequence
{
size_t AllocMem;
bool bSuffixOverlaps;
etSeqBase *pOverlapBase;
tsContaminantType *pContaminateType;
tsContamSeqNode *pCurNode;
tsContamSeqNodeBase *pCurNodeBase;
tsContamSeqNodeBase *pParentNodeBase;
int NodeBaseIdx;
int RemainSeqLen;

// validate parameter ranges
if(Type < eAOF5PE1Targ || Type > eAOF3PE2Targ || pSeq == NULL || SeqLen > cMaxContaminantLen || SeqLen < cMinContaminantLen)
	{
	Reset();
	return(eBSFerrParams);
	}


// if first time then allocate to hold initial contaminant sequence nodes
if(m_pContamSeqNodes == NULL)
	{
	// initial allocation for the tsContaminant's; will be realloc'd later as may be required
	AllocMem = (size_t)cAllocNumContamNodes * sizeof(tsContamSeqNode);
#ifdef _WIN32
	m_pContamSeqNodes = (tsContamSeqNode *)malloc(AllocMem);
	if (m_pContamSeqNodes == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "IndexContamSeq: unable to allocate for ContamSeqNodes memory %llu", AllocMem);
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and seems to have issues if more than 2GB allocation
	m_pContamSeqNodes = (tsContamSeqNode *)mmap(NULL, (size_t)AllocMem, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (m_pContamSeqNodes == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "IndexContamSeq: unable to allocate for ContamSeqNodes memory %llu", AllocMem);
		m_pContamSeqNodes = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_AllocdContamSeqNodesMem = AllocMem;
	m_AllocdContamSeqNodes = cAllocNumContamNodes;
	m_NumContamSeqNodes = 0;
	memset(m_pContamSeqNodes,0, m_AllocdContamSeqNodesMem);
	}

if(m_NumContamSeqNodes + cMaxContaminantLen + 10 >= m_AllocdContamSeqNodes)
	{
	tsContamSeqNode *pRealloc;
	AllocMem = (m_AllocdContamSeqNodes + cAllocNumContamNodes) * sizeof(tsContamSeqNode);
#ifdef _WIN32
	pRealloc = (tsContamSeqNode *)realloc(m_pContamSeqNodes, AllocMem);
#else
	pRealloc = (tsContamSeqNode *)mremap(m_pContamSeqNodes, m_AllocdContamSeqNodesMem, AllocMem, MREMAP_MAYMOVE);
	if (pRealloc == MAP_FAILED)
		pRealloc = NULL;
#endif
	if (pRealloc == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal, gszProcName, "IndexContamSeq: Memory re-allocation to %lld bytes - %s", AllocMem, strerror(errno));
		return(eBSFerrMem);
		}

	m_pContamSeqNodes = pRealloc;
	memset(&m_pContamSeqNodes[m_AllocdContamSeqNodes],0, AllocMem - m_AllocdContamSeqNodesMem);
	m_AllocdContamSeqNodesMem = AllocMem;
	m_AllocdContamSeqNodes += cAllocNumContamNodes;
	}

RemainSeqLen = SeqLen;
bSuffixOverlaps = (Type == eAOF5PE1Targ || Type == eAOF5PE2Targ) ? true : false;
if(bSuffixOverlaps)							// expecting to overlap the suffixes of the contaminants onto the prefixes of the the targeted sequences
	pOverlapBase = &pSeq[SeqLen-1];
else
	pOverlapBase = pSeq;					// expecting to overlap the prefixes of the contaminants onto the suffixes of the the targeted sequences

// if first of type then start a root node
pContaminateType = &m_ContaminantTypes[Type];
if(pContaminateType->RootContamSeqNodeIdx == 0)
	{
	pContaminateType->RootContamSeqNodeIdx = m_NumContamSeqNodes + 1;
	pCurNode = &m_pContamSeqNodes[m_NumContamSeqNodes++];
	memset(pCurNode,0,sizeof(tsContamSeqNode));
	pCurNode->NumNodeBases = 1;
	pCurNodeBase = pCurNode->Bases;
	pCurNodeBase->Base = *pOverlapBase;
	pCurNodeBase->ContamID = ContamID;
	pCurNodeBase->ChildNodeIdx = 0;
	pCurNodeBase->MaxContamSfxLen = RemainSeqLen;
	pParentNodeBase = pCurNodeBase;
	if(bSuffixOverlaps)
		pOverlapBase -= 1;
	else
		pOverlapBase += 1;
	RemainSeqLen -= 1;
	}
else
	{
	pCurNode = &m_pContamSeqNodes[pContaminateType->RootContamSeqNodeIdx-1];
	pParentNodeBase = NULL;
	}

do {
	if(pParentNodeBase != NULL)		// not null if starting a child node which is to be linked from the previous node base
		{
		pCurNode = &m_pContamSeqNodes[m_NumContamSeqNodes++];
		memset(pCurNode,0,sizeof(tsContamSeqNode));
		pCurNode->NumNodeBases = 1;
		pCurNodeBase = pCurNode->Bases;
		pCurNodeBase->Base = *pOverlapBase;
		pCurNodeBase->ContamID = ContamID;
		pCurNodeBase->MaxContamSfxLen = RemainSeqLen;
		pCurNodeBase->ChildNodeIdx = 0;
		pParentNodeBase->ChildNodeIdx = m_NumContamSeqNodes;
		pParentNodeBase = pCurNodeBase;
		}
	else
		{
		pCurNodeBase = pCurNode->Bases;
		for(NodeBaseIdx = 0; NodeBaseIdx < pCurNode->NumNodeBases; NodeBaseIdx++,pCurNodeBase++)
			{
			if(pCurNodeBase->Base == *pOverlapBase)
				{
				if(pCurNodeBase->MaxContamSfxLen < RemainSeqLen)
					pCurNodeBase->MaxContamSfxLen = RemainSeqLen;
				break;
				}
			}
		if(NodeBaseIdx == pCurNode->NumNodeBases)	// will be equal if base not already in current node
			{
			pCurNode->NumNodeBases += 1;
			pCurNodeBase->Base = *pOverlapBase;
			pCurNodeBase->ContamID = ContamID;
			pCurNodeBase->MaxContamSfxLen = RemainSeqLen;
			pCurNodeBase->ChildNodeIdx = 0;
			pParentNodeBase = pCurNodeBase;
			}
		else
			pCurNode = &m_pContamSeqNodes[pCurNodeBase->ChildNodeIdx-1];
		}
	if(bSuffixOverlaps)
		pOverlapBase -= 1;
	else
		pOverlapBase += 1;
	}
while(--RemainSeqLen);
return(0);
}

int			// < 0 if no matches, otherwise the number of substitutions required for the match 
CContaminants::RecursiveMatch(bool bSuffixOverlaps,  // true if processing for contaminant suffix overlaps onto target prefix, false if processing for contaminant prefix overlaps onto target suffix
				 int MaxSubs,				// maximum allowed substitutions
				int TargLen,				// target sequence prefix/suffix length 
				etSeqBase *pTargSeq,		// target sequence 	
				tsContamSeqNode *pCurNode,  // current node
				int *pContamID)				// match was to this contamination sequence

{
int Rslt;
int CurMaxSubs;
int CurContamID;
int BestReqSubs;
int BestContamID;
int NodeBaseIdx;
tsContamSeqNodeBase *pCurContamBase;
etSeqBase ContamBase;
etSeqBase TargBase;

CurContamID = 0;
BestContamID = 0;
BestReqSubs = 0;
TargBase = *pTargSeq & 0x07;

pCurContamBase = pCurNode->Bases;
for(NodeBaseIdx = 0; NodeBaseIdx < pCurNode->NumNodeBases; NodeBaseIdx++,pCurContamBase++)
	{
	if(pCurContamBase->MaxContamSfxLen < TargLen)	// check next base in node if the longest suffix of any current contaminant is shorter then the current target prefix
		continue;				
	CurMaxSubs = MaxSubs;
	ContamBase = pCurContamBase->Base & 0x07;
	if(TargBase == eBaseN || (ContamBase != eBaseN && TargBase != ContamBase))
		{
		if(CurMaxSubs == 0)			// if no more subs allowed then try for exact match on another nodebase
			continue;
		CurMaxSubs -= 1;			// one less sub allowed
		}
	
	// have accepted the target sequence base as matching the contaminant base
	if(TargLen == 1)	// if meeting overlap requirements then if required a substitution perhaps another contaminate base will exactly match
		{
		if(CurMaxSubs == MaxSubs)	// no sub reuired so immediately acceptable
			{
			*pContamID = pCurContamBase->ContamID;
			return(0);
			}
		if(CurContamID == 0)		// no other base matched so provisionally accepting this match but continue checking other bases in this node
			{
			BestReqSubs = 1;
			BestContamID = pCurContamBase->ContamID;
			}
		continue;
		}

	// still more target prefix bases to match on contaminant suffix bases, recurse down ...
	Rslt=RecursiveMatch(bSuffixOverlaps,CurMaxSubs,TargLen-1,bSuffixOverlaps ? pTargSeq-1 : pTargSeq + 1,&m_pContamSeqNodes[pCurContamBase->ChildNodeIdx-1],&CurContamID);
	if(Rslt < 0)		// if no matches then try next base in current node
		continue;

	if(Rslt == 0 && CurMaxSubs == MaxSubs)		// match which required no substitutions; can't better so accept ...
		{
		*pContamID = CurContamID;
		return(0);
		}

	if(BestContamID == 0 || Rslt < BestReqSubs)
		{
		BestContamID = CurContamID;
		BestReqSubs = Rslt + MaxSubs - CurMaxSubs;
		}
	}
*pContamID = BestContamID;
return(BestContamID == 0 ? -1 : BestReqSubs);
}


const int cMaxSeedIters = 25000;		// max number of iterations (exploration depth) with current QueryWinLen bases looking for marching seed window to extend
int					// returns number of mismatches or -1 if no hits within AllowSubsRate
CContaminants::MatchVectContam(int AllowSubsRate,// if non-zero then allow substitutions in the containg contaminants at this rate per 25bp of overlap length if overlap >= 10bp
					int QueryLen,			// query sequence length
					etSeqBase *pQuerySeq,	// attempt to locate a maximal containment onto this query read sequence
					tsVectContam *pVectContam) // checking for containment from this vector sequence
{
int SeedIters;
int AllowedMMs;
int CmpIdx;
int NumMMs;
int LowestMMs;
int HitIdx;
int BaseIdx;
etSeqBase QBase;
etSeqBase CBase;
int QueryWinLen;
int Ofs;
etSeqBase *pQueryBase;
etSeqBase *pContamBase;

if(QueryLen > 25 && AllowSubsRate > 0)
	{
	AllowedMMs = (QueryLen / 25) * AllowSubsRate;
	QueryWinLen = QueryLen / (AllowedMMs+1);
	LowestMMs = AllowedMMs + 1;
	for(Ofs = 0; Ofs < QueryLen; Ofs += QueryWinLen)
		{
		if((Ofs + QueryWinLen) > QueryLen)
			Ofs = QueryLen - QueryWinLen;
		SeedIters = 0;
		while(SeedIters++ <= cMaxSeedIters)
			{
			if(SeedIters == 1)
				{
				HitIdx = LocateFirstExact(&pQuerySeq[Ofs],		// pts to probe sequence
								 QueryWinLen,			// probe length to exactly match over
						  pVectContam->pBases,			// target sequence
						  pVectContam->pSfxIdx,			// target sequence suffix array
						  0,							// low index in pSfxArray
						  pVectContam->ContamLen-1);	// high index in pSfxArray
				if(HitIdx == 0)
					break;
				}
			else
				{
				HitIdx = LocateNextExact(HitIdx,						// previously hit this 
											&pQuerySeq[Ofs],			// pts to probe sequence
												 QueryWinLen,			// probe length to exactly match over
										  pVectContam->pBases,			// target sequence
										  pVectContam->pSfxIdx,			// target sequence suffix array
										  0,							// low index in pSfxArray
										  pVectContam->ContamLen-1);	// high index in pSfxArray
				if(HitIdx == 0)
					break;

				}
		
			BaseIdx = pVectContam->pSfxIdx[HitIdx-1];
			if(BaseIdx < Ofs || (BaseIdx + QueryLen - Ofs) > pVectContam->ContamLen)
				continue;
			BaseIdx -= Ofs;
			pQueryBase = pQuerySeq;
			pContamBase = &pVectContam->pBases[BaseIdx];
			NumMMs = 0;

			for(CmpIdx = 0; CmpIdx < QueryLen; CmpIdx++,pQueryBase++,pContamBase++)
				{
				if(((QBase = *pQueryBase & 0x0f) == eBaseEOS) || ((CBase = *pContamBase & 0x0f) == eBaseEOS))
					break;
				if(QBase != CBase)
					{
					NumMMs += 1;
					if(NumMMs >= LowestMMs)
						break;
					}
				}
			if(CmpIdx == QueryLen && NumMMs < LowestMMs)	
				LowestMMs = NumMMs;
			if(LowestMMs == 0)
				return(0);
			}
		}
	return(LowestMMs > AllowedMMs ? -1 : LowestMMs);
	}
else
	{
	HitIdx =		// index+1 in pSfxArray of first exactly matching probe or 0 if no match
		LocateFirstExact(pQuerySeq,				// pts to probe sequence
						 QueryLen,				// probe length to exactly match over
				  pVectContam->pBases,			// target sequence
				  pVectContam->pSfxIdx,			// target sequence suffix array
				  0,							// low index in pSfxArray
				  pVectContam->ContamLen-1);	// high index in pSfxArray
	}

return(HitIdx > 0 ? 0 : -1);
}

int			// 0 if no Contaminant contains this target sequence
CContaminants::MatchVectContams(bool bIsPE2,		// true if target sequence is a PE2, else if false then target is either a SE or PE1
					int AllowSubsRate,			// if non-zero then allow substitutions in the overlapping Contaminants at this rate per 25bp of overlap length if overlap >= 10bp
					int QueryLen,			// query sequence length
					etSeqBase *pQuerySeq)	// attempt to locate a maximal overlap onto this query read sequence
{
etSeqBase QuerySeq[cMaxContamQuerySeqLen+1];
etSeqBase RevCplQuerySeq[cMaxContamQuerySeqLen+1];
int NumSubs;
int LowestNumSubs;
tsVectContam *pBestVectContam;
int VectIdx;
tsVectContam *pVectContam;

if(m_NumVectContaminates == 0 || AllowSubsRate < 0 || QueryLen < cMinContamQuerySeqLen || QueryLen > cMaxContamQuerySeqLen || pQuerySeq == NULL)
	return(0);
memcpy(QuerySeq,pQuerySeq,QueryLen);
QuerySeq[QueryLen] = eBaseEOS;
memcpy(RevCplQuerySeq,pQuerySeq,QueryLen);
CSeqTrans::ReverseComplement(QueryLen,RevCplQuerySeq);
RevCplQuerySeq[QueryLen] = eBaseEOS;

LowestNumSubs = -1;
pBestVectContam = NULL;
pVectContam = m_ContaminantVectors;
for(VectIdx = 0; VectIdx < m_NumVectContaminates; VectIdx++, pVectContam++)
	{
	if(!bIsPE2 && !(pVectContam->FlgPE1Antisense || pVectContam->FlgPE1Sense))
		continue;
	if(bIsPE2 && !(pVectContam->FlgPE2Antisense || pVectContam->FlgPE2Sense))
		continue;
	if(pVectContam->ContamLen < QueryLen)
		continue;

	// attempt to find a match
	AcquireSerialise();
	m_ContaminantTypes[eAOFVector].NumChecks += 1;
	ReleaseSerialise();
	if(!bIsPE2 && pVectContam->FlgPE1Sense || bIsPE2 && pVectContam->FlgPE2Sense)
		{
		if((NumSubs = MatchVectContam(AllowSubsRate,QueryLen,QuerySeq,pVectContam)) >= 0)
			{
			if(NumSubs == 0)
				{
				AcquireSerialise();
				m_ContaminantTypes[eAOFVector].HitTot += 1;
				pVectContam->HitTot += 1;
				ReleaseSerialise();
				return(QueryLen);
				}
			if(pBestVectContam == NULL || NumSubs < LowestNumSubs)
				{
				LowestNumSubs = NumSubs;
				pBestVectContam = pVectContam;
				}
			}
		}

	// if antisense allowed then try for a match
	if(!bIsPE2 && pVectContam->FlgPE1Antisense || bIsPE2 && pVectContam->FlgPE2Antisense)
		{
		if((NumSubs = MatchVectContam(AllowSubsRate,QueryLen,RevCplQuerySeq,pVectContam)) >= 0)
			{
			if(NumSubs == 0)
				{
				AcquireSerialise();
				m_ContaminantTypes[eAOFVector].HitTot += 1;
				pVectContam->HitTot += 1;
				ReleaseSerialise();
				return(QueryLen);
				}
			if(pBestVectContam == NULL || NumSubs < LowestNumSubs)
				{
				LowestNumSubs = NumSubs;
				pBestVectContam = pVectContam;
				}
			}
		}
	}
if(pBestVectContam != NULL)
	{
	AcquireSerialise();
	m_ContaminantTypes[eAOFVector].HitTot += 1;
	pBestVectContam->HitTot += 1;
	ReleaseSerialise();
	return(QueryLen);
	}
return(0);
}


int			// < 0 if errors, 0 if no Contaminant overlap, 1..N number of Contaminant suffix bases overlapping
CContaminants::MatchContaminants(teContamType Type,	// process for this overlay type
					int AllowSubsRate,		// if non-zero then allow substitutions in the overlapping Contaminants at this rate per 25bp of overlap length if overlap >= 10bp
					int MinOverlap,			// minimum required overlap
					int QueryLen,			// query sequence length
					etSeqBase *pQuerySeq)	// attempt to locate a contaminate flanking sequence which overlays onto this query read sequence
{
int Rslt;
bool bSuffixOverlaps;
int CurOverlapLen;
int OverlapIdx;
int ContamID;
int MaxAcceptedSubs;

etSeqBase *pTargBase;
tsContamSeqNode *pTypeRootNode;
tsContaminantType *pContaminantType;
tsFlankContam *pContaminant;

if(QueryLen < cMinContamQuerySeqLen || QueryLen > cMaxContamQuerySeqLen)
	return(0);

if(pQuerySeq == NULL || MinOverlap < 0 || AllowSubsRate < 0 || AllowSubsRate > 3)
	return(eBSFerrInternal);

// check if anything to do, any contaminants of requested type and if then are they long enough to overlap by at least the requested minimum required overlap?
if (m_TotNumContaminants == 0)	// nothing to do?
	return(0);

// any vector contaminates are processed for containment of this target sequence before processing for contaminate overlays
if(m_NumVectContaminates)
	{
	if((Rslt = MatchVectContams((Type == eAOF5PE1Targ || Type == eAOF5PE1Targ) ? false : true,AllowSubsRate,QueryLen,pQuerySeq)) != 0)
		return(Rslt);
	}

if(MinOverlap < 1)
	MinOverlap = 1;

pContaminantType = &m_ContaminantTypes[Type];
if (pContaminantType->NumContaminants == 0 || pContaminantType->RootContamSeqNodeIdx == 0)
	return(0);
CurOverlapLen = min(QueryLen, pContaminantType->MaxContamSeqLen);
if (CurOverlapLen < MinOverlap)
	return(0);
AcquireSerialise();
pContaminantType->NumChecks += 1;		
ReleaseSerialise();
bSuffixOverlaps = (Type == eAOF5PE1Targ || Type == eAOF5PE2Targ) ? true : false;
if(bSuffixOverlaps)
	pTargBase = &pQuerySeq[CurOverlapLen-1];
else
	pTargBase = &pQuerySeq[QueryLen - CurOverlapLen];

pTypeRootNode = &m_pContamSeqNodes[pContaminantType->RootContamSeqNodeIdx-1];

// have at least one contaminant of requested type which is at least the minimum required length to be an otherlap
for (OverlapIdx = CurOverlapLen; OverlapIdx >= MinOverlap; OverlapIdx--,bSuffixOverlaps ? pTargBase -= 1 : pTargBase += 1 )
	{
	if(AllowSubsRate > 0 && OverlapIdx >= 10)
		MaxAcceptedSubs = (AllowSubsRate * (OverlapIdx + 15)) / 25;
	else
		MaxAcceptedSubs = 0;
		
	if((Rslt=RecursiveMatch(bSuffixOverlaps,MaxAcceptedSubs,OverlapIdx,pTargBase,pTypeRootNode,&ContamID)) >= 0)
		{
		pContaminant = &m_pContaminants[ContamID - 1];
		AcquireSerialise();
		pContaminantType->HitTot += 1;
		pContaminantType->HitDist[OverlapIdx] += 1;
		pContaminant->HitTot += 1;
		pContaminant->HitDist[OverlapIdx] += 1;
		ReleaseSerialise();
		return(OverlapIdx);
		}
	}

return(0);
}


int 
CContaminants::NumOfContaminants(teContamClass ComtamClass)			// returns number of contaminants loaded 
{
if(m_TotNumContaminants == 0)
	return(0);
switch(ComtamClass) {
	case eCCFlankContam:			// flank contaminant only
		return(m_NumFlankContaminates);
	case eCCVectContam:				// vector contaminant only
		return(m_NumVectContaminates);
	case eCCAllContam:
	default:						// total number of both flank and vector contaminants
		break;
	}
return(m_TotNumContaminants);
}

int 
CContaminants::MaxContaminantLen(teContamClass ComtamClass)			// returns longest length of any contaminant sequence 
{
if(m_TotNumContaminants == 0)
	return(0);
switch(ComtamClass) {
	case eCCFlankContam:			// flank contaminant only
		return(m_MaxFlankContamSeqLen);
	case eCCVectContam:				// vector contaminant only
		return(m_MaxVectContamSeqLen);
	case eCCAllContam:
	default:						// longest of either  flank and vector contaminants
		break;
	}
return(m_MaxFlankContamSeqLen >= m_MaxVectContamSeqLen ? m_MaxFlankContamSeqLen : m_MaxVectContamSeqLen);
}

int 
CContaminants::MinContaminantLen(teContamClass ComtamClass)			// returns shortest length of any contaminant sequence 
{
if(m_TotNumContaminants == 0)
	return(0);
switch(ComtamClass) {
	case eCCFlankContam:			// flank contaminant only
		return(m_MinFlankContamSeqLen);
	case eCCVectContam:				// vector contaminant only
		return(m_MinVectContamSeqLen);
	case eCCAllContam:
	default:						// shortest of either  flank and vector contaminants > 0
		break;
	}
if(m_NumFlankContaminates == 0)
	return(m_MinVectContamSeqLen);
if(m_NumVectContaminates == 0)
	return(m_MinFlankContamSeqLen);
return(m_MinFlankContamSeqLen <= m_MaxVectContamSeqLen ? m_MinFlankContamSeqLen : m_MaxVectContamSeqLen);
}

UINT32 
CContaminants::NumChecks(teContamType Type)				// returns number of times this contaminant type was checked for an overlap onto a target sequence
{
if(m_TotNumContaminants == 0 || Type < eAOF5PE1Targ || Type > eAOFVector)
	return(0);

return(m_ContaminantTypes[(int)Type].NumChecks);
}

teContamType											// returned contaminant type ( -1 if unable to locate ContamID)
CContaminants::ContaminantType(int ContamID)			// contaminant identifier
{
int Idx;
if(m_TotNumContaminants == 0 || ContamID < 1 || ContamID > m_TotNumContaminants)
	return((teContamType)-1);

if(m_CacheContamID == ContamID)
	{
	if(m_CacheContamIDClass == eCCFlankContam)
		return(m_pContaminants[m_CacheContamIDIdx].Type);
	else
		return(eAOFVector);
	}
m_CacheContamID = 0;
if(m_NumFlankContaminates)
	{
	tsFlankContam *pFlankContam; 
	pFlankContam = m_pContaminants;
	for(Idx = 0; Idx < m_NumFlankContaminates; Idx++,pFlankContam++)
		if(pFlankContam->ContamID == ContamID)
			{
			m_CacheContamID = ContamID;
			m_CacheContamIDClass = eCCFlankContam;
			m_CacheContamIDIdx = Idx;
			return(pFlankContam->Type);
			}
	}
if(m_NumVectContaminates)
	{
	tsVectContam *pVectContam;
	pVectContam = m_ContaminantVectors;
	for(Idx = 0; Idx < m_NumVectContaminates; Idx++,pVectContam++)
		if(pVectContam->ContamID == ContamID)
			{
			m_CacheContamID = ContamID;
			m_CacheContamIDClass = eCCVectContam;
			m_CacheContamIDIdx = Idx;
			return(eAOFVector);
			}
	}
return((teContamType)-1);
}

char *												// returned contaminant name (NULL if unable to locate ContamID)
CContaminants::ContaminantName(int ContamID)		// contaminant identifier
{
int Idx;
if(m_TotNumContaminants == 0 || ContamID < 1 || ContamID > m_TotNumContaminants)
	return(NULL);
if(m_CacheContamID == ContamID)
	{
	if(m_CacheContamIDClass == eCCFlankContam)
		return(m_pContaminants[m_CacheContamIDIdx].szName);
	else
		return(m_ContaminantVectors[m_CacheContamIDIdx].szName);
	}
m_CacheContamID = 0;
if(m_NumFlankContaminates)
	{
	tsFlankContam *pFlankContam; 
	pFlankContam = m_pContaminants;
	for(Idx = 0; Idx < m_NumFlankContaminates; Idx++,pFlankContam++)
		if(pFlankContam->ContamID == ContamID)
			{
			m_CacheContamID = ContamID;
			m_CacheContamIDClass = eCCFlankContam;
			m_CacheContamIDIdx = Idx;
			return(pFlankContam->szName);
			}
	}
if(m_NumVectContaminates)
	{
	tsVectContam *pVectContam;
	pVectContam = m_ContaminantVectors;
	for(Idx = 0; Idx < m_NumVectContaminates; Idx++,pVectContam++)
		if(pVectContam->ContamID == ContamID)
			{
			m_CacheContamID = ContamID;
			m_CacheContamIDClass = eCCVectContam;
			m_CacheContamIDIdx = Idx;
			return(pVectContam->szName);
			}
	}
return(NULL);
}

teContamClass										// returned contaminant class (-1 if unable to locate ContamID)
CContaminants::ContaminantClass(int ContamID)		// contaminant identifier
{
int Idx;
if(m_TotNumContaminants == 0 || ContamID < 1 || ContamID > m_TotNumContaminants)
	return((teContamClass)-1);
if(m_CacheContamID == ContamID)
	return(m_CacheContamIDClass);
m_CacheContamID = 0;
if(m_NumFlankContaminates)
	{
	tsFlankContam *pFlankContam; 
	pFlankContam = m_pContaminants;
	for(Idx = 0; Idx < m_NumFlankContaminates; Idx++,pFlankContam++)
		if(pFlankContam->ContamID == ContamID)
			{
			m_CacheContamID = ContamID;
			m_CacheContamIDClass = eCCFlankContam;
			m_CacheContamIDIdx = Idx;
			return(eCCFlankContam);
			}
	}
if(m_NumVectContaminates)
	{
	tsVectContam *pVectContam;
	pVectContam = m_ContaminantVectors;
	for(Idx = 0; Idx < m_NumVectContaminates; Idx++,pVectContam++)
		if(pVectContam->ContamID == ContamID)
			{
			m_CacheContamID = ContamID;
			m_CacheContamIDClass = eCCVectContam;
			m_CacheContamIDIdx = Idx;
			return(eCCVectContam);
			}
	}
return((teContamClass)-1);
}


int												// returned contaminant identifier  (-1 if unable to locate name)
CContaminants::ContaminantID(char *pszName)		// returns contaminant identifer for name
{
int Idx;
if(!m_TotNumContaminants || pszName == NULL || pszName[0] == '\0')
	return(-1);

if(m_CacheContamID > 0)
	{
	char *pCacheName;
	if(m_CacheContamIDClass == eCCFlankContam)
		pCacheName = m_pContaminants[m_CacheContamIDIdx].szName;
	else
		pCacheName = m_ContaminantVectors[m_CacheContamIDIdx].szName;
	if(!stricmp(pCacheName,pszName))
		return(m_CacheContamID);
	}
m_CacheContamID = 0;
if(m_NumFlankContaminates)
	{
	tsFlankContam *pFlankContam; 
	pFlankContam = m_pContaminants;
	for(Idx = 0; Idx < m_NumFlankContaminates; Idx++,pFlankContam++)
		if(!stricmp(pFlankContam->szName,pszName))
			{
			m_CacheContamID = pFlankContam->ContamID;
			m_CacheContamIDClass = eCCFlankContam;
			m_CacheContamIDIdx = Idx;
			return(m_CacheContamID);
			}
	}
if(m_NumVectContaminates)
	{
	tsVectContam *pVectContam;
	pVectContam = m_ContaminantVectors;
	for(Idx = 0; Idx < m_NumVectContaminates; Idx++,pVectContam++)
		if(!stricmp(pVectContam->szName,pszName))
			{
			m_CacheContamID = pVectContam->ContamID;
			m_CacheContamIDClass = eCCVectContam;
			m_CacheContamIDIdx = Idx;
			return(m_CacheContamID);
			}
	}
return(-1);
}

int												// returned contaminant length (-1 if unable to locate ContamID)
CContaminants::ContaminantLen(int ContamID)		// contaminant identifier
{
int Idx;
if(m_TotNumContaminants == 0 || ContamID < 1 || ContamID > m_TotNumContaminants)
	return(-1);
if(m_CacheContamID == ContamID)
	{
	if(m_CacheContamIDClass == eCCFlankContam)
		return(m_pContaminants[m_CacheContamIDIdx].ContamLen);
	else
		return(m_ContaminantVectors[m_CacheContamIDIdx].ContamLen);
	}
m_CacheContamID = 0;
if(m_NumFlankContaminates)
	{
	tsFlankContam *pFlankContam; 
	pFlankContam = m_pContaminants;
	for(Idx = 0; Idx < m_NumFlankContaminates; Idx++,pFlankContam++)
		if(pFlankContam->ContamID == ContamID)
			{
			m_CacheContamID = pFlankContam->ContamID;
			m_CacheContamIDClass = eCCFlankContam;
			m_CacheContamIDIdx = Idx;
			return(pFlankContam->ContamLen);
			}
	}
if(m_NumVectContaminates)
	{
	tsVectContam *pVectContam;
	pVectContam = m_ContaminantVectors;
	for(Idx = 0; Idx < m_NumVectContaminates; Idx++,pVectContam++)
		if(pVectContam->ContamID == ContamID)
			{
			m_CacheContamID = pVectContam->ContamID;
			m_CacheContamIDClass = eCCVectContam;
			m_CacheContamIDIdx = Idx;
			return(pVectContam->ContamLen);
			}
	}
return(-1);
}
		
int										// total number of overlaps 	 					
CContaminants::ContaminantDist(int ContamID,		// contaminant distributions for this contaminant
						int LenCnts,		// return for up to this length overlaps, ignored for vector contaminates
						int *pCnts)		// returned overlap length counts, ignored for vector contaminates 
{
int Idx;
if(LenCnts && pCnts != NULL)
	memset(pCnts,0,sizeof(int)*LenCnts);
if(m_TotNumContaminants == 0 || ContamID < 1 || ContamID > m_TotNumContaminants)
	return(0);

if(m_CacheContamID == ContamID)
	{
	if(m_CacheContamIDClass == eCCFlankContam)
		{
		if(LenCnts > 0 && pCnts != NULL)
			{
			if(LenCnts > m_pContaminants[m_CacheContamIDIdx].ContamLen)
				 LenCnts = m_pContaminants[m_CacheContamIDIdx].ContamLen;
			memcpy(pCnts,&m_pContaminants[m_CacheContamIDIdx].HitDist[1],LenCnts * sizeof(UINT32));
			}
		return(m_pContaminants[m_CacheContamIDIdx].HitTot);
		}
	else
		return(m_ContaminantVectors[m_CacheContamIDIdx].HitTot);
	}

m_CacheContamID = 0;
if(m_NumFlankContaminates)
	{
	tsFlankContam *pFlankContam; 
	pFlankContam = m_pContaminants;
	for(Idx = 0; Idx < m_NumFlankContaminates; Idx++,pFlankContam++)
		if(pFlankContam->ContamID == ContamID)
			{
			if(LenCnts > 0 && pCnts != NULL)
				{
				if(LenCnts > pFlankContam->ContamLen)
					 LenCnts = pFlankContam->ContamLen;
				memcpy(pCnts,&pFlankContam->HitDist[1],LenCnts * sizeof(UINT32));
				}
			m_CacheContamID = pFlankContam->ContamID;
			m_CacheContamIDClass = eCCFlankContam;
			m_CacheContamIDIdx = Idx;
			return(pFlankContam->HitTot);
			}
	}
if(m_NumVectContaminates)
	{
	tsVectContam *pVectContam; 
	pVectContam = m_ContaminantVectors;
	for(Idx = 0; Idx < m_NumVectContaminates; Idx++,pVectContam++)
		if(pVectContam->ContamID == ContamID)
			{
			m_CacheContamID = pVectContam->ContamID;
			m_CacheContamIDClass = eCCFlankContam;
			m_CacheContamIDIdx = Idx;
			return(pVectContam->HitTot);
			}
	}

return(0);
}

int												// index+1 in pSfxArray of first exactly matching probe or 0 if no match
CContaminants::LocateFirstExact(etSeqBase *pProbe,  // pts to probe sequence
				  int ProbeLen,					// probe length to exactly match over
				  etSeqBase *pTarg,				// target sequence
				  int *pSfxArray,				// target sequence suffix array
				  int SfxLo,					// low index in pSfxArray
				  int SfxHi)					// high index in pSfxArray
{
etSeqBase *pEl1;
etSeqBase *pEl2;
UINT8 El1;
UINT8 El2;

int CmpRslt;
int Ofs;
int Mark;
int TargPsn;
do {
	pEl1 = pProbe;
	TargPsn = (int)(((INT64)SfxLo + SfxHi) / 2L);
	pEl2 = &pTarg[pSfxArray[TargPsn]];
	CmpRslt = 0;
	for(Ofs=0; Ofs < ProbeLen; Ofs++)
		{
		El2 = *pEl2++ & 0x0f;
		if(El2 == eBaseEOS)
			{
			CmpRslt = -1;
			break;
			}
		El1 = *pEl1++ & 0x0f;
		if(El1 > El2)
			{
			CmpRslt = 1;
			break;
			}
		if(El1 < El2)
			{
			CmpRslt = -1;
			break;
			}
		}

	if(!CmpRslt)	// if a match then may not be the lowest indexed match
		{
		if(TargPsn == 0 || SfxLo == TargPsn) // check if already lowest
			return(TargPsn + 1);
		// iterate until lowest located
		while(1) {
			if(CmpRslt == 0)
				{
				Mark = TargPsn;
				if(Mark == 0)
					return(1);
				SfxHi = TargPsn - 1;
				}
			TargPsn = (int)(((INT64)SfxLo + SfxHi) / 2L);
			pEl2 = &pTarg[pSfxArray[TargPsn]];

			pEl1 = pProbe;
			CmpRslt = 0;
			for(Ofs=0; Ofs < ProbeLen; Ofs++)
				{
				El2 = *pEl2++ & 0x0f;
				if(El2 == eBaseEOS)
					{
					CmpRslt = -1;
					break;
					}
				El1 = *pEl1++ & 0x0f;
				if(El1 > El2)
					{
					CmpRslt = 1;
					break;
					}
				if(El1 < El2)
					{
					CmpRslt = -1;
					break;
					}
				}

			if(CmpRslt == 0)				// 0 if still matching
				continue;
			SfxLo = TargPsn + 1;
			if(SfxLo == Mark)
				return(Mark+1);
			}
		}

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		SfxHi = TargPsn - 1;
		}
	else
		SfxLo = TargPsn+1;
	}
while(SfxHi >= SfxLo);

return(0);	// unable to locate any instance of pProbe
}

int			// index+1 in pSfxArray of last exactly matching probe or 0 if no match
CContaminants::LocateLastExact(etSeqBase *pProbe, // pts to probe sequence
				  int ProbeLen,					// probe length to exactly match over
				  etSeqBase *pTarg,				// target sequence
				  int *pSfxArray,			// target sequence suffix array
				  int SfxLo,					// low index in pSfxArray
				  int SfxHi)					// high index in pSfxArray
{
etSeqBase *pEl1;
etSeqBase *pEl2;
UINT8 El1;
UINT8 El2;

int CmpRslt;
int Ofs;
int Mark;
int TargPsn;
int SfxHiMax = SfxHi;
int SfxLoMax = SfxLo;
do {
	pEl1 = pProbe;
	TargPsn = (int)(((INT64)SfxLo + SfxHi) / 2L);
	pEl2 = &pTarg[pSfxArray[TargPsn]];
	CmpRslt = 0;
	for(Ofs=0; Ofs < ProbeLen; Ofs++)
		{
		El2 = *pEl2++ & 0x0f;
		if(El2 == eBaseEOS)
			{
			CmpRslt = -1;
			break;
			}
		El1 = *pEl1++ & 0x0f;
		if(El1 > El2)
			{
			CmpRslt = 1;
			break;
			}
		if(El1 < El2)
			{
			CmpRslt = -1;
			break;
			}
		}

	if(!CmpRslt)	// if a match then may not be the highest indexed match
		{
		if(TargPsn == SfxHiMax || SfxHi == TargPsn) // check if already highest
			return(TargPsn + 1);
		// iterate until highest located
		while(1) {
			if(CmpRslt == 0)
				{
				Mark = TargPsn;
				if(Mark == SfxHi)
					return(Mark+1);
				SfxLo = TargPsn + 1;
				}
			TargPsn = (int)(((INT64)SfxLo + SfxHi) / 2L);
			pEl2 = &pTarg[pSfxArray[TargPsn]];

			pEl1 = pProbe;
			CmpRslt = 0;
			for(Ofs=0; Ofs < ProbeLen; Ofs++)
				{
				El2 = *pEl2++ & 0x0f;
				if(El2 == eBaseEOS)
					{
					CmpRslt = -1;
					break;
					}
				El1 = *pEl1++ & 0x0f;
				if(El1 > El2)
					{
					CmpRslt = 1;
					break;
					}
				if(El1 < El2)
					{
					CmpRslt = -1;
					break;
					}
				}

			if(CmpRslt == 0)				// 0 if still matching
				continue;
			SfxHi = TargPsn - 1;
			if(SfxHi == Mark)
				return(Mark+1);
			}
		}

	if(CmpRslt < 0)
		{
		if(TargPsn == 0)
			break;
		SfxHi = TargPsn - 1;
		}
	else
		SfxLo = TargPsn+1;
	}
while(SfxHi >= SfxLo);

return(0);	// unable to locate any instance of pProbe
}

int			// index+1 in pSfxArray of last exactly matching probe or 0 if no match
CContaminants::LocateNextExact(int CurSfxIdx,	// previously returned index either from LocateFirstExact() or LocateNextExact()
				  etSeqBase *pProbe,			// pts to probe sequence
				  int ProbeLen,					// probe length to exactly match over
				  etSeqBase *pTarg,				// target sequence
				  int *pSfxArray,				// target sequence suffix array
				  int SfxLo,					// low index in pSfxArray
				  int SfxHi)					// high index in pSfxArray
{
etSeqBase *pTargBase;
etSeqBase TargBase;
etSeqBase ProbeBase;
if(CurSfxIdx == 0)
	return(LocateFirstExact(pProbe,ProbeLen,pTarg,pSfxArray,SfxLo,SfxHi));
if((CurSfxIdx - 1) > SfxHi)
	return(0);
pTargBase = &pTarg[pSfxArray[CurSfxIdx]];
while(ProbeLen && (((TargBase = (*pTargBase & 0x0f)) <= eBaseN) && ((ProbeBase = (*pProbe & 0x0f)) <= eBaseN)))
	{
	ProbeLen -= 1;
	if(ProbeBase != TargBase)
		return(0);
	}
return(ProbeLen == 0 ? CurSfxIdx+1 : 0);
}

// SortContamTypeLen
// Sort contaminants by type then ContamLen descending then by ContamID ascending
int
CContaminants::SortContamTypeLen(const void *arg1, const void *arg2)
{
tsFlankContam *pEl1 = (tsFlankContam *)arg1;
tsFlankContam *pEl2 = (tsFlankContam *)arg2;

if(pEl1->Type > pEl2->Type)
	return(1);
if(pEl1->Type < pEl2->Type)
	return(-1);

if(pEl1->ContamLen < pEl2->ContamLen)
	return(1);
if(pEl1->ContamLen > pEl2->ContamLen )
		return(-1);

if(pEl1->ContamID > pEl2->ContamID)
	return(1);
if(pEl1->ContamID < pEl2->ContamID)
	return(-1);
return(0);
}

tsVectContam *CContaminants::m_pCurVector2Index; // current vector sequence to index with IndexVectorSeq() 

int
CContaminants::IndexVectorSeq(const void *arg1, const void *arg2)
{
UINT8 *pSeq1;
UINT8 *pSeq2;
etSeqBase b1;
etSeqBase b2;
if(m_pCurVector2Index == NULL)
	return(0);
pSeq1 = &m_pCurVector2Index->pBases[*(INT32 *)arg1];
pSeq2 = &m_pCurVector2Index->pBases[*(INT32 *)arg2];

while(1)
	{
	b1 = (*pSeq1++ & 0x0f);
	b2 = (*pSeq2++ & 0x0f);
	if(b1 > eBaseN && b2 > eBaseN)
		break;
	if(b1 != b2)
		return(b1 < b2 ? -1 : 1);
	}
return(0);
}