// CSVMerge.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32
#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif
#include "ChromMaps.h"

// Ref + Rel feature bits are stored packed 2 bits
const int cRefRelBitsPerInt = 16;							// asumes sizeof(unsigned int) is 32, and 2bits for ref+rel packed
const unsigned int cMSBitIntMsk = 0xc0000000;				// MSBit mask

CChromMaps::CChromMaps(void)
{
m_pLastLocChrom = NULL;
m_NumChromMaps=0;			// number of chrom maps in use
m_pChromMaps=NULL;			// linked list of chrom maps
m_ppSortedChromMaps = NULL;	// pts to array of ptrs to chrom maps sorted by chrom name
}

CChromMaps::~CChromMaps(void)
{
tsChromMap *pNxtMap;
while(m_pChromMaps != NULL)
	{
	if(m_pChromMaps->pMap != NULL)
		delete m_pChromMaps->pMap;
	pNxtMap = m_pChromMaps->pNext;
	delete m_pChromMaps;
	m_pChromMaps = pNxtMap;
	}

if(m_ppSortedChromMaps != NULL)
	delete m_ppSortedChromMaps;
}


int
CChromMaps::AddChrom(char *pszChrom,int StartLoci,int EndLoci)
{
tsChromMap *pChrom;
tsChromMap *pPrvChrom;

pChrom = m_pChromMaps;
pPrvChrom = NULL;

if(m_pLastLocChrom != NULL &&
   !stricmp(m_pLastLocChrom->szChrom,pszChrom))
	{
	pChrom = m_pLastLocChrom;
	if(pChrom->StartLoci == -1 || pChrom->StartLoci > StartLoci)
		pChrom->StartLoci = StartLoci;
	if(pChrom->EndLoci < EndLoci)
		pChrom->EndLoci = EndLoci;
	return(eBSFSuccess);
	}

// see if chrom previously added
while(pChrom != NULL)
	{
	if(!stricmp(pChrom->szChrom,pszChrom))
		break;
	pPrvChrom = pChrom;
	pChrom = pChrom->pNext;
	}

// need to allocate more maps?
if(pChrom == NULL)
	{
	if((pChrom = new tsChromMap)==NULL)
		return(eBSFerrMem);
	memset(pChrom,0,sizeof(tsChromMap));
	pChrom->StartLoci = -1;
	if(pPrvChrom != NULL)
		pPrvChrom->pNext = pChrom;
	else
		m_pChromMaps = pChrom;
	m_NumChromMaps += 1;
	strncpy(pChrom->szChrom,pszChrom,sizeof(pChrom->szChrom));
	pChrom->szChrom[sizeof(pChrom->szChrom)-1] = '\0';
	}
if(pChrom->StartLoci == -1 || pChrom->StartLoci > StartLoci)
	pChrom->StartLoci = StartLoci;
if(pChrom->EndLoci < EndLoci)
	pChrom->EndLoci = EndLoci;
m_pLastLocChrom = pChrom;
return(eBSFSuccess);
}

int
CChromMaps::InitMaps(void)
{
tsChromMap *pChrom;
tsChromMap **pChromSort;

if(m_NumChromMaps == 0 || m_pChromMaps == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"No elements to process");
	return(eBSFerrNoEntries);
	}

if((m_ppSortedChromMaps = new tsChromMap *[m_NumChromMaps])==NULL)
	return(eBSFerrMem);
pChromSort = m_ppSortedChromMaps;
m_pLastLocChrom = NULL;
pChrom = m_pChromMaps;
while(pChrom != NULL)
	{
	int MapLen = (cRefRelBitsPerInt + pChrom->EndLoci) / cRefRelBitsPerInt;
	if((pChrom->pMap = new unsigned int [MapLen])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory %d bytes for holding element merge map",MapLen * sizeof(unsigned int));
		return(eBSFerrMem);
		}
	memset(pChrom->pMap,0,sizeof(unsigned int) * MapLen);
	pChrom->MapLen = MapLen;
	*pChromSort++ = pChrom;
	pChrom = pChrom->pNext;
	}
qsort(m_ppSortedChromMaps,m_NumChromMaps,sizeof(tsChromMap *),SortChromNames);
return(eBSFSuccess);
}

tsChromMap *
CChromMaps::LocateChromMap(char *pszChrom)
{
tsChromMap *pProbe;
char *pszTarg;
int Rslt;
if(pszChrom == NULL || *pszChrom == '\0')
	return(NULL);
if(m_pLastLocChrom != NULL && !stricmp(m_pLastLocChrom->szChrom,pszChrom))
	return(m_pLastLocChrom);
int Lo,Mid,Hi;	// search limits
Lo = 0; Hi = m_NumChromMaps-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pProbe = m_ppSortedChromMaps[Mid];
	pszTarg = pProbe->szChrom;
	Rslt = stricmp(pszChrom,pszTarg);
	if(Rslt < 0)	
		{
		Hi = Mid - 1;
		continue;
		}
	if(Rslt > 0)	
		{
		Lo = Mid + 1;
		continue;
		}
	m_pLastLocChrom = pProbe;
	return(pProbe);
	}
return(NULL);
}


// SetChromMap
// pChrom->pMap |= bits between Start and End inclusive
// Returns number of bits set or < 0 if parameter errors
int
CChromMaps::SetChromMap(bool bIsRel,char *pszChrom,int Start,int End)
{
int StartBit;
unsigned int Bit2Set;
unsigned int Bits;
unsigned int Int2Set;
unsigned int LSBit;
int NumBits;
int NumSet;
int StartOfs;				
unsigned int *pBits;
tsChromMap *pChrom;

if(pszChrom == NULL || pszChrom[0] == '\0' || Start < 0 || Start > End)
	return(eBSFerrParams);

if((pChrom = LocateChromMap(pszChrom))==NULL)
	return(eBSFerrChrom);
if(End > pChrom->EndLoci || Start < pChrom->StartLoci)
	return(eBSFerrInternal);

if(bIsRel)
	{
	Bit2Set = 0x02;
	Int2Set = 0xaaaaaaaa;
	}
else
	{
	Bit2Set = 0x01;
	Int2Set = 0x55555555;
	}

NumSet = NumBits = 1 + End - Start;
StartOfs = Start/cRefRelBitsPerInt;				
pBits = &pChrom->pMap[StartOfs];

if(StartBit = (Start % cRefRelBitsPerInt))		// initial bit to set, 0 if on an int boundary
	{
	Bits = LSBit = Bit2Set << (2*StartBit);	// StartBit will be 1..BitsPerInt-1
	while(--NumBits && !(Bits & cMSBitIntMsk)) 
		{
		Bits <<= 2;
		Bits |= LSBit;
		}
	*pBits++ |= Bits;
	}

while(NumBits >= cRefRelBitsPerInt)
	{
	*pBits++ |= Int2Set;
	NumBits -= cRefRelBitsPerInt;
	}

if(NumBits > 0)
	{
	Bits = Bit2Set;
	while(--NumBits)
		{
		Bits <<= 2;
		Bits |= Bit2Set;
		}
	*pBits |= Bits;
	}
return(NumSet);
}

tsChromMap *
CChromMaps::LocateChrom(char *pszChrom)
{
tsChromMap *pChrom;
pChrom = m_pChromMaps;
while(pChrom != NULL)
	{
	if(!stricmp(pChrom->szChrom,pszChrom))
		return(pChrom);
	pChrom = pChrom->pNext;
	}
return(NULL);
}

int 
CChromMaps::Load(bool bSkipFirst,int MinLength,int MaxLength,int RefExtend,char *pszCSVRef,int RelExtend,char *pszCSVRel)
{
int Rslt;

// iterate over ref + rel CSV creating chrom maps and determining element lengths
if((Rslt = LoadChroms(false,false,bSkipFirst,MinLength,MaxLength,RefExtend,pszCSVRef)) < 0)
	return(Rslt);

if(pszCSVRel != NULL && pszCSVRel[0] != '\0')
	if((Rslt = LoadChroms(false,true,bSkipFirst,MinLength,MaxLength,RelExtend,pszCSVRel)) < 0)
		return(Rslt);

// sort chrom maps and allocate memory for map vectors
if((Rslt = InitMaps())<0)
	return(Rslt);

// iterate over ref + rel CSV setting maps
if((Rslt = LoadChroms(true,false,bSkipFirst,MinLength,MaxLength,RefExtend,pszCSVRef))<0)
	return(Rslt);
if(pszCSVRel != NULL && pszCSVRel[0] != '\0')
	Rslt = LoadChroms(true,true,bSkipFirst,MinLength,MaxLength,RelExtend,pszCSVRel);
return(Rslt);
}


int
CChromMaps::LoadChroms(bool bSetMaps,bool bIsRel,bool bSkipFirst,int MinLength,int MaxLength,int FlankExtend,char *pszCSVFile)
{
int Rslt;
int NumFields;
int Len;
int NumAdded;
int NumProcessed;
char *pszChrom;
int StartLoci;
int EndLoci;
int NumUnderLen;
int NumOverLen;
int NumFlankUnderLen;
int NumFlankOverLen;

CCSVFile *pCSV;

if((pCSV = new CCSVFile)== NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	return(eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszCSVFile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszCSVFile);
	delete pCSV;
	return(Rslt);
	}


NumUnderLen = 0;
NumOverLen = 0;
NumFlankUnderLen = 0;
NumFlankOverLen = 0;

Rslt = eBSFSuccess;
NumAdded = 0;
NumProcessed = 0;
while(Rslt >= 0 && (Rslt=pCSV->NextLine()) > 0)	// onto next line containing fields
	{
	if(bSkipFirst)					// skip header line?
		{
		bSkipFirst = false;
		continue;
		}
	
	NumFields = pCSV->GetCurFields();
	if(NumFields < 7)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected 7+ fields in '%s', GetCurFields() returned '%d'",pszCSVFile,NumFields);
		Rslt = eBSFerrFieldCnt;
		break;
		}

	NumProcessed += 1;
	pCSV->GetInt(7,&Len);
	if(FlankExtend <= 0 && Len < MinLength)
		{
		NumUnderLen += 1;
		continue;
		}
		
	if(FlankExtend >= 0 && Len > MaxLength)
		{
		NumOverLen += 1;
		continue;
		}

	pCSV->GetText(4,&pszChrom);
	pCSV->GetInt(5,&StartLoci);
	pCSV->GetInt(6,&EndLoci);
	if(FlankExtend > 0)
		{
		StartLoci -= FlankExtend;
		if(StartLoci < 0)
			StartLoci = 0;
		EndLoci += FlankExtend;
		Len = 1 + EndLoci - StartLoci;
		}
	else
		if(FlankExtend < 0)
			{
			StartLoci += abs(FlankExtend);
			EndLoci -= abs(FlankExtend);
			if(EndLoci < 0)
				EndLoci = 0;
			if(EndLoci < StartLoci)
				Len = 0;
			else
				Len = 1 + EndLoci - StartLoci;
			}
	if(Len > MaxLength)
		{
		NumFlankOverLen += 1;
		continue;
		}
	if(!Len || Len < MinLength)
		{
		NumFlankUnderLen += 1;
		continue;
		}

	if(bSetMaps)
		{
		Rslt = SetChromMap(bIsRel,pszChrom,StartLoci,EndLoci);
		}
	else
		{
		Rslt = AddChrom(pszChrom,StartLoci,EndLoci);
		}
	if(Rslt < 0)
		printf("\nProblem!");
	else
		NumAdded += 1;
	if(!(NumAdded % 1000))
		printf("\rChr: %s Count: %d        ",pszChrom,NumAdded);
	}
printf("\rChr: %s Count: %d        ",pszChrom,NumAdded);
if(pCSV != NULL)
	delete pCSV;

if(!bSetMaps)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"%s - Accepted: %d Processed: %d Underlength: %d Overlength: %d UnderFlankLength: %d OverFlankLength: %d",
			pszCSVFile,NumAdded,NumProcessed,NumUnderLen,NumOverLen,NumFlankUnderLen,NumFlankOverLen);
	}
return(Rslt);
}

int
CChromMaps::Process(etElMapType Mode,int (*Handler)(char *pszChrom,int StartLoci,int EndLoci,void *pParm),void *pParm)
{
int Rslt;
int Idx;
int Loci;
int ElLen;
bool bInEl;
int StartLoci;
int EndLoci;
unsigned int NormBits;
unsigned int *pMap;
tsChromMap *pChrom;

Rslt = eBSFSuccess;
pChrom = *m_ppSortedChromMaps;
for(Idx = 0; Rslt >= eBSFSuccess && Idx < m_NumChromMaps; Idx++)
	{
	pChrom = m_ppSortedChromMaps[Idx];
	pMap = pChrom->pMap;
	ElLen = 0;
	bInEl = false;
	for(Loci = 0; Loci <= pChrom->EndLoci; Loci++)
		{
		if(Loci && !(Loci % cRefRelBitsPerInt))
			pMap++;
		bInEl = false;
		NormBits = 0x03 & (*pMap >> (2*(Loci % cRefRelBitsPerInt)));
		switch(Mode) {
			case eElIntersect:				// intersect, both Ref and Rel must be present
				if(NormBits == 0x03)
					bInEl = true;
				break;
			case eElRefExclusive:			// ref exclusive, Ref must be present and Rel absent
				if(NormBits == 0x01)
					bInEl = true;
				break;
			case eElRelExclusive:			// rel exclusive, Ref must be absent and Rel present
				if(NormBits == 0x02)
					bInEl = true;
				break;
			case eElRefRelUnion:			// union, Ref or Rel must be present
				if(NormBits != 0x00)
					bInEl = true;
				break;
			case eElRefNotRefRel:			// neither present
				if(NormBits == 0x00)
					bInEl = true;
				break;
			}
		if(bInEl)
			{
			if(!ElLen++)
				StartLoci = Loci;
			EndLoci = Loci;
			}
		else
			if(ElLen)
				{
				Rslt = Handler(pChrom->szChrom,StartLoci,EndLoci,pParm);
				ElLen = 0;
				}
		}
	if(ElLen)
		Rslt = Handler(pChrom->szChrom,StartLoci,EndLoci,pParm);
	}
return(Rslt);
}



// SortChromNames
// Used to sort by chromosome names
int 
CChromMaps::SortChromNames( const void *arg1, const void *arg2)
{
tsChromMap *pEl1 = *(tsChromMap **)arg1;
tsChromMap *pEl2 = *(tsChromMap **)arg2;

char *pszName1 = &pEl1->szChrom[0];
char *pszName2 = &pEl2->szChrom[0];
char c1 = tolower(*pszName1);
char c2 = tolower(*pszName2);
if(c1 < c2)
	return(-1);
if(c1 > c2)
	return(1);
return(stricmp(pszName1,pszName2));
}
