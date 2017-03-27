/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

CAlignValidate::CAlignValidate(void)
{
m_pAlignSeqs = NULL;
m_ppAlignSeqsIdx = NULL;
m_pSpecies = NULL;
m_pChroms = NULL;
m_pStats = NULL;
Reset();
}

CAlignValidate::~CAlignValidate(void)
{
if(m_pAlignSeqs != NULL)
	delete m_pAlignSeqs;

if(m_ppAlignSeqsIdx != NULL)
	delete m_ppAlignSeqsIdx;

if(m_pSpecies != NULL)
	delete m_pSpecies;

if(m_pChroms != NULL)
	delete m_pChroms;

if(m_pStats != NULL)
	delete m_pStats;
}

int
CAlignValidate::Reset(void)
{
if(m_pAlignSeqs != NULL)
	{
	delete m_pAlignSeqs;
	m_pAlignSeqs = NULL;
	}
m_NumAligns = 0;
m_AllocdAlignSeqsInsts = 0;

if(m_ppAlignSeqsIdx != NULL)
	{
	delete m_ppAlignSeqsIdx;
	m_ppAlignSeqsIdx = NULL;
	}
	
if(m_pSpecies != NULL)
	{
	delete m_pSpecies;
	m_pSpecies = NULL;
	}
m_NumSpecies = 0;
m_AllocdSpeciesInsts = 0;

if(m_pChroms != NULL)
	{
	delete m_pChroms;
	m_pChroms = NULL;
	}
m_NumChroms = 0;
m_AllocdChromsInsts = 0;

if(m_pStats != NULL)
	{
	delete m_pStats;
	m_pStats = NULL;
	}
return(eBSFSuccess);
}

int
CAlignValidate::ProcessChroms(char *pszSpeciesChromFile,bool bChkDups)
{
FILE *pStream;
int LineNum;
char szLineBuff[512];
int Cnt;
char chr;
char *pSrc;
int Rslt;
int Idx;
char szSpecies[cMaxDatasetSpeciesChrom];
char szChrom[cMaxDatasetSpeciesChrom];
int ChromLen;

if((pStream = fopen(pszSpeciesChromFile,"r"))==NULL)
	{
	printf("Unable to open species chromosome file %s error: %s",pszSpeciesChromFile,strerror(errno));
	return(eBSFerrOpnFile);
	}
LineNum = 0;
while(fgets(szLineBuff,sizeof(szLineBuff),pStream)!= NULL)
	{
	LineNum++;
	if(strlen(szLineBuff) < 5)	// simply slough lines which are too short to contain minimal chrom spec
		continue;
	// strip any leading whitespace
	pSrc = szLineBuff;
	while(chr = *pSrc)
		{
		if(!isspace(chr))
			break;
		pSrc++;
		}
	if(chr == '\0')
		continue;

	 Cnt = sscanf(pSrc," %25[^ \t,] , %25[^ \t,] , %d",
			szSpecies,szChrom,&ChromLen);

	 if(Cnt != 3)
		{
		printf("Error parsing species chromosome file %s at line %d, expected species then chrom then len %d parameters\n%s\n",pszSpeciesChromFile,LineNum,Cnt,szLineBuff);
		fclose(pStream);
		return(eBSFerrStructParm); 
		}
	if((Rslt=AddSpeciesChrom(szSpecies,szChrom,ChromLen,bChkDups)) < 1)
		{
		fclose(pStream);
		return(Rslt);
		}
	}
fclose(pStream);

if(m_pChroms != NULL && m_NumChroms > 1)
	qsort(m_pChroms,m_NumChroms,sizeof(tsASchromname),SortSpeciesChroms);
for(Idx = 0; Idx < m_NumChroms; Idx++)
	m_pChroms[Idx].ChromID = Idx+1;

return(eBSFSuccess);
}

int
CAlignValidate::Open(char *pszAlignFile)
{
FILE *pStream;
int LineNum;
char szLineBuff[512];
int Cnt;
char chr;
char *pSrc;
char *pDst;
int SrcFileID;


char szRefSpecies[cMaxDatasetSpeciesChrom];
char szRefChrom[cMaxDatasetSpeciesChrom];
char szRelSpecies[cMaxDatasetSpeciesChrom];
char szRelChrom[cMaxDatasetSpeciesChrom];
tsASalignSeq *pAlign;

if(m_pAlignSeqs == NULL)
	{
	m_pAlignSeqs = (tsASalignSeq *)new tsASalignSeq [cAllocAlignSeqIncr];
	if(m_pAlignSeqs == NULL)
		return(eBSFerrMem);
	m_AllocdAlignSeqsInsts = cAllocAlignSeqIncr;
	}
m_NumAligns = 0;

if((pStream = fopen(pszAlignFile,"r"))==NULL)
	{
	printf("Unable to open alignment file %s error: %s",pszAlignFile,strerror(errno));
	delete m_pAlignSeqs;
	m_pAlignSeqs = NULL;
	return(eBSFerrOpnFile);
	}
LineNum = 0;
m_NumAligns = 0;
pAlign = m_pAlignSeqs;

while(fgets(szLineBuff,sizeof(szLineBuff),pStream)!= NULL)
	{
	LineNum++;
	if(strlen(szLineBuff) < 5)	// simply slough lines which are too short to contain anything worth parsing
		continue;

	// strip any whitespace and quotes
	pDst = pSrc = szLineBuff;
	while(chr = *pSrc++)
		if(!isspace(chr) && chr != '\'' && chr != '"')
			*pDst++ = chr;
	*pDst = '\0';
	if(szLineBuff[0] == '\0')
		continue;

	// grow m_pAlignSeqs when needed...
	if(m_NumAligns >= m_AllocdAlignSeqsInsts)
		{
		pAlign = (tsASalignSeq *)new tsASalignSeq [m_AllocdAlignSeqsInsts + cAllocAlignSeqIncr];
		if(pAlign == NULL)
			{
			fclose(pStream);
			printf("Failed reallocing memory for extending m_pAlignBA");
			delete m_pAlignSeqs;
			m_pAlignSeqs = NULL;
			return(eBSFerrMem);
			}
		if(m_NumAligns)
			memmove(pAlign,m_pAlignSeqs,m_NumAligns * sizeof(tsASalignSeq));
		delete m_pAlignSeqs;
		m_pAlignSeqs = pAlign;
		m_AllocdAlignSeqsInsts += cAllocAlignSeqIncr;
		}

	 Cnt = sscanf(szLineBuff,"%d,%d,%d,%d,%d,%25[^,],%25[^,],%d,%25[^,],%25[^,],%d,%c",
			&SrcFileID,&pAlign->BlockID,&pAlign->SubSeqID,&pAlign->Len,&pAlign->Matches,szRefSpecies,szRefChrom,&pAlign->RefOfs,szRelSpecies,szRelChrom,&pAlign->RelOfs,&pAlign->Strand);
	 if(Cnt != 12)
		{
		printf("Error parsing alignment file %s at line %d, expected 12 but only parsed %d parameters\n%s\n",pszAlignFile,LineNum,Cnt,szLineBuff);
		fclose(pStream);
		delete m_pAlignSeqs;
		m_pAlignSeqs = NULL;
		return(eBSFerrStructParm); 
		}
	pAlign->SrcFileID = (short int)SrcFileID;
	if((pAlign->RefChromID = GetChromID(szRefSpecies,szRefChrom))<1)
		return(pAlign->RefChromID);

	if((pAlign->RelChromID = GetChromID(szRelSpecies,szRelChrom)) < 1)
		return(pAlign->RefChromID);
	pAlign++;
	m_NumAligns++;
	}
fclose(pStream);

// now sort m_pAlignSeqs by RefChromID.RefOfs
qsort(m_pAlignSeqs,m_NumAligns,sizeof(tsASalignSeq),SortAlignRef);

return(eBSFSuccess);
}

int
CAlignValidate::Close(void)
{
return(Reset());
}

int
CAlignValidate::OutputStats(char *pszRsltsFile)
{
FILE *pStream;
int Len;
tsSeqStat *pSeqStat;

if((pStream = fopen(pszRsltsFile,"w+"))==NULL)
	{
	printf("Unable to open/create results file %s error: %s",pszRsltsFile,strerror(errno));
	return(eBSFerrOpnFile);
	}
fprintf(pStream,"Strand,Length,Seqs,Bases,Matches,Ident75,Ident50,Ident25,BasesPutative,ErrChromOrStrands,ErrChroms,ErrStrands,Miss1,Miss5,Miss10,Miss25,Miss50,Miss100,Miss1K,Miss10K,Miss100K,SeqsReciprocal,BasesReciprocal");
pSeqStat = &m_pStats->SeqsPlusStrand[0];
for(Len = 1; Len <= cMaxSeqLenStats; Len++,pSeqStat++)
	if(pSeqStat->Cnt)
		{
		fprintf(pStream,"\n'+',%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
			    Len,pSeqStat->Cnt,Len * pSeqStat->Cnt,pSeqStat->Idents,
				pSeqStat->Ident75,pSeqStat->Ident50,pSeqStat->Ident25,
				pSeqStat->NumPutative,pSeqStat->ErrChromStrand,pSeqStat->ErrChrom,pSeqStat->ErrStrand,
				pSeqStat->NumPutativeMiss[0],pSeqStat->NumPutativeMiss[1],pSeqStat->NumPutativeMiss[2],pSeqStat->NumPutativeMiss[3],
				pSeqStat->NumPutativeMiss[4],pSeqStat->NumPutativeMiss[5],pSeqStat->NumPutativeMiss[6],pSeqStat->NumPutativeMiss[7],pSeqStat->NumPutativeMiss[8],
				pSeqStat->CntReciprocal,pSeqStat->Reciprocal);
		}
pSeqStat = &m_pStats->SeqsMinusStrand[0];
for(Len = 1; Len <= cMaxSeqLenStats; Len++,pSeqStat++)
	if(pSeqStat->Cnt)
		{
		fprintf(pStream,"\n'-',%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
			    Len,pSeqStat->Cnt,Len * pSeqStat->Cnt,pSeqStat->Idents,
				pSeqStat->Ident75,pSeqStat->Ident50,pSeqStat->Ident25,
				pSeqStat->NumPutative,pSeqStat->ErrChromStrand,pSeqStat->ErrChrom,pSeqStat->ErrStrand,
				pSeqStat->NumPutativeMiss[0],pSeqStat->NumPutativeMiss[1],pSeqStat->NumPutativeMiss[2],pSeqStat->NumPutativeMiss[3],
				pSeqStat->NumPutativeMiss[4],pSeqStat->NumPutativeMiss[5],pSeqStat->NumPutativeMiss[6],pSeqStat->NumPutativeMiss[7],pSeqStat->NumPutativeMiss[8],
				pSeqStat->CntReciprocal,pSeqStat->Reciprocal);
		}
fprintf(pStream,"\n");
fclose(pStream);
return(eBSFSuccess);
}

// AddSpeciesChrom
// If Species.Chrom previously added then returns the existing identifer
// Otherwise adds Species.Chrom and returns a count of all Species.Chrom added
// Option to check if species/chrom already added
int
CAlignValidate::AddSpeciesChrom(char *pszSpecies,char *pszChrom,int ChromLen,bool bChkDups)
{
int SpeciesID;
int ChromID;
int Idx;
tsASchromname *pChrom;
tsASspeciesname *pSpecies;
unsigned short Hash;

if(m_pSpecies == NULL)
	{
	m_pSpecies = (tsASspeciesname *)new tsASspeciesname [cAllocSpeciesNamesIncr];
	if(m_pSpecies == NULL)
		return(eBSFerrMem);
	m_AllocdSpeciesInsts = cAllocSpeciesNamesIncr;
	m_NumSpecies = 0;
	}
if(m_pChroms == NULL)
	{
	m_pChroms = (tsASchromname *)new tsASchromname [cAllocChromNamesIncr];
	if(m_pChroms == NULL)
		return(eBSFerrMem);
	m_AllocdChromsInsts = cAllocChromNamesIncr;
	m_NumChroms = 0;
	}

// Species already known?
if((SpeciesID = GetSpeciesID(pszSpecies)) < 1)
	{
	// new species...
	if(m_NumSpecies >= m_AllocdSpeciesInsts)
		{
		pSpecies = (tsASspeciesname *)new tsASspeciesname [m_AllocdSpeciesInsts + cAllocSpeciesNamesIncr];
		if(pSpecies == NULL)
			return(eBSFerrMem);
		if(m_NumSpecies)
			memmove(pSpecies,m_pSpecies,m_NumSpecies * sizeof(tsASspeciesname));
		delete m_pSpecies;
		m_pSpecies = pSpecies;
		m_AllocdSpeciesInsts += cAllocSpeciesNamesIncr;
		}

	pSpecies = &m_pSpecies[m_NumSpecies++];
	SpeciesID = pSpecies->SpeciesID = m_NumSpecies;
	pSpecies->Hash = GenNameHash(pszSpecies);
	strncpy(pSpecies->szName,pszSpecies,cMaxDatasetSpeciesChrom);
	pSpecies->szName[cMaxDatasetSpeciesChrom-1] = '\0';
	pSpecies->AssembLen = 0;
	}
else
	pSpecies = &m_pSpecies[SpeciesID-1];

// check if chromosome already known
Hash = GenNameHash(pszChrom);
ChromID = 0;
if(bChkDups)
	{
	if(m_pChroms != NULL && m_NumChroms)
		{
		pChrom = m_pChroms;
		for(Idx = 0; Idx < m_NumChroms; Idx++,pChrom++)
			if(pChrom->SpeciesID == SpeciesID && pChrom->Hash == Hash && !stricmp(pChrom->szName,pszChrom))
				{
				ChromID = pChrom->ChromID;
				break;
				}
		}
	}

if(!ChromID)
	{
	if(m_NumChroms >= m_AllocdChromsInsts)
		{
		pChrom = (tsASchromname *)new tsASchromname [m_AllocdChromsInsts + cAllocChromNamesIncr];
		if(pChrom == NULL)
			return(eBSFerrMem);
		if(m_NumChroms)
			memmove(pChrom,m_pChroms,m_NumChroms * sizeof(tsASchromname));
		delete m_pChroms;
		m_pChroms = pChrom;
		m_AllocdChromsInsts += cAllocChromNamesIncr;
		}

	pChrom = &m_pChroms[m_NumChroms++];
	ChromID = pChrom->ChromID = m_NumChroms;
	pChrom->SpeciesID = SpeciesID;
	pChrom->ChromLen = ChromLen;
	pSpecies->AssembLen += ChromLen;
	pChrom->Hash = Hash;
	strncpy(pChrom->szName,pszChrom,cMaxDatasetSpeciesChrom);
	pChrom->szName[cMaxDatasetSpeciesChrom-1] = '\0';
	}

return(m_NumChroms);
}


int
CAlignValidate::GetChromID(char *pszSpecies,char *pszChrom)
{
int SpeciesID;
if((SpeciesID = GetSpeciesID(pszSpecies)) >= 1)
	return(GetChromID(SpeciesID,pszChrom));
return(0);
}

// GetChromID
// Binary search
int
CAlignValidate::GetChromID(int SpeciesID,char *pszChrom)
{
tsASchromname *pChrom;
int Left;
int Right;
int MidPt;
int Rslt;
if(!m_NumChroms || m_pChroms == NULL || pszChrom == NULL || *pszChrom == '\0')
	return(0);
Left = 0;
Right = m_NumChroms - 1;
MidPt;
while(Right >= Left) {
	MidPt = (Right + Left)/2;
	pChrom = &m_pChroms[MidPt];
	if(pChrom->SpeciesID > SpeciesID)
		{
		Right = MidPt - 1;
		continue;
		}
	if(pChrom->SpeciesID < SpeciesID)
		{
		Left = MidPt + 1;
		continue;
		}
	// matching on species identifier, check on chromosome name

	Rslt = stricmp(pszChrom,pChrom->szName);

	if(Rslt > 0)
		{
		Left = MidPt + 1;
		continue;
		}
	if(Rslt < 0)
		{
		Right = MidPt - 1;
		continue;
		}
	// located species.chromosome
	return(pChrom->ChromID);
	}
return(0);	// could not locate species.chromosome
}

// GetSpeciesID
// Linear search as number of species will generally be low
int
CAlignValidate::GetSpeciesID(char *pszSpecies)
{
int Idx;
unsigned short int Hash;
tsASspeciesname *pSpecies;
if(!m_NumSpecies || m_pSpecies == NULL || pszSpecies == NULL || *pszSpecies == '\0')
	return(0);
Hash = GenNameHash(pszSpecies);
pSpecies = m_pSpecies;
for(Idx = 0; Idx < m_NumSpecies; Idx++,pSpecies++)
	if(pSpecies->Hash == Hash &&
		!strncmp(pSpecies->szName,pszSpecies,cMaxDatasetSpeciesChrom))
		return(pSpecies->SpeciesID);
return(0);
}

int
CAlignValidate::GetChromLen(int ChromID)
{
if(ChromID < 1 || ChromID > m_NumChroms)
	return(0);
return(m_pChroms[ChromID-1].ChromLen);
}

INT64 
CAlignValidate::GetSpeciesLen(int SpeciesID)	// returns sum of all chromosome lengths for species
{
if(SpeciesID < 1 || SpeciesID > m_NumSpecies)
	return(0);
return(m_pSpecies[SpeciesID-1].AssembLen);
}


// Strand2Other
// Returns '-' strand psn for Psn supplied on '+' strand,
// or '+' strand psn for Psn supplied on '-' strand
//
int
CAlignValidate::MAFStrand2Other(char *pszSpecies,char *pszChrom,int Psn)
{
int SpeciesID;
if((SpeciesID = GetSpeciesID(pszSpecies)) >= 1)
	return(MAFStrand2Other(GetChromID(SpeciesID,pszChrom),Psn));
return(0);
}

int
CAlignValidate::AXTStrand2Other(char *pszSpecies,char *pszChrom,int Psn)
{
int SpeciesID;
if((SpeciesID = GetSpeciesID(pszSpecies)) >= 1)
	return(AXTStrand2Other(GetChromID(SpeciesID,pszChrom),Psn));
return(0);
}

int
CAlignValidate::MAFStrand2Other(int ChromID,int Psn)
{
int ChromLen;
if(ChromID < 1 || ChromID > m_NumChroms)
	return(0);
ChromLen = m_pChroms[ChromID-1].ChromLen;
if(Psn < 0)			// be nice, assume < 0 is 0
	Psn = ChromLen-1;
else
	{
	if(Psn >= ChromLen)
		Psn = 0;
	else
		Psn = (ChromLen - 1) - Psn;
	}
return(Psn);
}

int
CAlignValidate::AXTStrand2Other(int ChromID,int Psn)
{
int ChromLen;
if(ChromID < 1 || ChromID > m_NumChroms)
	return(0);
ChromLen = m_pChroms[ChromID-1].ChromLen;
if(Psn <= 1)			// be nice and assume that caller thinks it's a 0..n coordinate system
	Psn = ChromLen;
else
	{
	if(Psn >= ChromLen)
		Psn = 1;
	else
		Psn = (ChromLen - Psn) + 1;
	}
return(Psn);
}

// GenNameHash
// Generates a 16bit hash on specified lowercased name
// This hash can then be used to quickly eliminate probe names which can't match a target name by comparing hashes
unsigned short
CAlignValidate::GenNameHash(char *pszName)
{
unsigned short Hash = 0;
unsigned short MSB;
char chr;
while(chr = *pszName++)
	{
	MSB = (Hash & 0x0c000) >> 13;
	Hash <<= 2;
	Hash ^= tolower(chr);
	Hash ^= MSB;
	}
return(Hash);
}

int 
CAlignValidate::LocateReciprocals(char *pszAlignFile,char *pszErrAlignFile,char *pszMissesFile)
{
FILE *pStream;
FILE *pErrStream;
FILE *pMissesStream;
int LineNum;
char szLineBuff[512];
int Cnt;
char chr;
char *pSrc;
char *pDst;

int SrcFileID;
int BlockID;
int SubSeqID;
int Len;
int Matches;
int RefOfs;
int RelOfs;
char Strand;
int NumReciprocal;

char szRefSpecies[cMaxDatasetSpeciesChrom];
char szRefChrom[cMaxDatasetSpeciesChrom];
char szRelSpecies[cMaxDatasetSpeciesChrom];
char szRelChrom[cMaxDatasetSpeciesChrom];



if(m_pStats == NULL)
	{
	m_pStats = new tsReciprocalStats;
	if(m_pStats == NULL)
		return(eBSFerrMem);
	}
memset(m_pStats,0,sizeof(tsReciprocalStats));

if((pStream = fopen(pszAlignFile,"r"))==NULL)
	{
	printf("Unable to open alignment file %s error: %s",pszAlignFile,strerror(errno));
	delete m_pAlignSeqs;
	m_pAlignSeqs = NULL;
	return(eBSFerrOpnFile);
	}
if(pszErrAlignFile != NULL && pszErrAlignFile[0] != '\0')
	{
	if((pErrStream = fopen(pszErrAlignFile,"w"))==NULL)
		{
		printf("Unable to create alignment errors file %s error: %s",pszErrAlignFile,strerror(errno));
		fclose(pStream);
		delete m_pAlignSeqs;
		m_pAlignSeqs = NULL;
		return(eBSFerrOpnFile);
		}
	}
else
	pErrStream = NULL;
	

if(pszMissesFile != NULL && pszMissesFile[0] != '\0')
	{
	if((pMissesStream = fopen(pszMissesFile,"w"))==NULL)
		{
		printf("Unable to create alignment misses file %s error: %s",pszMissesFile,strerror(errno));
		fclose(pStream);
		if(pErrStream != NULL)
			fclose(pErrStream);
		delete m_pAlignSeqs;
		m_pAlignSeqs = NULL;
		return(eBSFerrOpnFile);
		}
	}
else
	pMissesStream = NULL;

LineNum = 0;
while(fgets(szLineBuff,sizeof(szLineBuff),pStream)!= NULL)
	{
	LineNum++;
	if(strlen(szLineBuff) < 5)	// simply slough lines which are too short to contain anything worth parsing
		continue;

	// strip any whitespace and quotes
	pDst = pSrc = szLineBuff;
	while(chr = *pSrc++)
		if(!isspace(chr) && chr != '\'' && chr != '"')
			*pDst++ = chr;
	*pDst = '\0';
	if(szLineBuff[0] == '\0')
		continue;

	 Cnt = sscanf(szLineBuff," %d , %d , %d , %d , %d , %25[^ \t,] , %25[^ \t,] , %d , %25[^ \t,] , %25[^ \t,] , %d , %c",
			&SrcFileID,&BlockID,&SubSeqID,&Len,&Matches,szRefSpecies,szRefChrom,&RefOfs,szRelSpecies,szRelChrom,&RelOfs,&Strand);
	 if(Cnt != 12)
		{
		printf("Error parsing alignment file %s at line %d, expected 12 but only parsed %d parameters\n%s\n",pszAlignFile,LineNum,Cnt,szLineBuff);
		fclose(pStream);
		fclose(pErrStream);
		return(eBSFerrStructParm); 
		}

	NumReciprocal = CalcSeqStats(pMissesStream,
								SrcFileID,BlockID,SubSeqID,
								Len,		// alignment length
								Matches,		// claimed number of exact identities in this aligned sequence
								szRefSpecies,	// A species
								szRefChrom,      // source chromosome
								RefOfs,				// source offset
								szRelSpecies,	// B species
								szRelChrom,		// B chromosome
								RelOfs,				// B offset
								Strand);			// B strand

	// this is where we could write out those blocks which did/did'nt reciprocally align for
	// further analysis by some other program
	if(pErrStream != NULL && NumReciprocal != Len)
		{
		// this block has at least one base which is not reciprocally aligned...
		fprintf(pErrStream,"%d,%d,%d,%d,%d,%s,%s,%d,%s,%s,%d,%c,%d\n",
			SrcFileID,BlockID,SubSeqID,Len,Matches,szRefSpecies,szRefChrom,RefOfs,szRelSpecies,szRelChrom,RelOfs,Strand,NumReciprocal);
		}
	}
fclose(pStream);
if(pErrStream != NULL)
	fclose(pErrStream);
if(pMissesStream != NULL)
	fclose(pMissesStream);
printf("\nNumClaimedAligned: %d, NumClaimedMatched: %d, NumReciprocal: %d",m_pStats->NumClaimedAligned,m_pStats->NumClaimedMatched,m_pStats->NumReciprocal);
return(eBSFSuccess);
}


// CalcSeqStats
// Returns number of bases which are reciprocal between alignment AB and BA
int						// returns number of bases reciprocally aligned
CAlignValidate::CalcSeqStats(FILE *pMisses,				// where to write detail on misses
							 int SrcFileID,				// original source file from which alignment was processed
							 int BlockID,				// alignment block identifier in source file
							 int SubSeqID,				// which subsequence in alignment block
							 int Len,					// alignment length
								int Matches,			// claimed number of exact identities in this aligned sequence
								char *pszRefSpecies,	// A species
								char *pszRefChrom,      // source chromosome
								int RefOfs,				// source offset
								char *pszRelSpecies,	// B species
								char *pszRelChrom,		// B chromosome
								int RelOfs,				// B offset
								char Strand)			// B strand
{
tsSeqStat *pStat;
bool bHitDiffChrom;
bool bHitDiffStrand;
bool bHitDiffChromStrand;
int MissDistance;
int NearMissDistance[11];
int Idents;
int NumReciprocal;
int NumPutative;
int Ofs;
tsASalignSeq *pProbe;
int Diff;
int RelSpeciesID;
int RelChromID;
int RefSpeciesID;
int RefChromID;
int Left;
int Right;
int MidPt;
int InitLen = Len;
int InitRefOfs = RefOfs;
int InitRelOfs = RelOfs;
//int ChromLen;

if(Len > cMaxSeqLenStats)
	return(0);
if(!(RelSpeciesID = GetSpeciesID(pszRelSpecies)))
   return(0);
if(!(RelChromID = GetChromID(RelSpeciesID,pszRelChrom)))
	return(0);
if(!(RefSpeciesID = GetSpeciesID(pszRefSpecies)))
   return(0);
if(!(RefChromID = GetChromID(RefSpeciesID,pszRefChrom)))
	return(0);


Idents = (Matches*100)/Len;				// calc sequence identity
m_pStats->NumClaimedAligned += Len;		// update global stats
m_pStats->NumClaimedMatched += Matches;

if(Strand == '+')						// update stats for sequences of this length and strand
	pStat = &m_pStats->SeqsPlusStrand[Len-1];
else
	pStat = &m_pStats->SeqsMinusStrand[Len-1];
pStat->Cnt++;
pStat->Idents += Matches;
if(Idents >= 75)						// >= 75% identity
	pStat->Ident75++;
if(Idents >= 50)						// >= 50% identity
	pStat->Ident50++;
if(Idents >= 25)						// >= 25% identity
	pStat->Ident25++;

// locate on BA alignment all bases which have same reference chrom.offset as AB relative chrom.ofs
NumReciprocal = 0;
NumPutative = 0;
pProbe = NULL;
bHitDiffChrom = false;
bHitDiffStrand = false;
bHitDiffChromStrand = false;
memset(NearMissDistance,0,sizeof(NearMissDistance));
for(Ofs = 0; Ofs < Len; Ofs++,RelOfs++,RefOfs++)
	{
	if(pProbe != NULL)		// optimisation: if had success last probe then chances are that pProbe is still valid!
		{
		if(pProbe->RefOfs <= RelOfs &&	
			((pProbe->RefOfs + pProbe->Len) - 1) >= RelOfs)
			{
			NumPutative++;
			NumReciprocal++;
			continue;
			}
		}
	
	Left = 0;
	Right = m_NumAligns - 1;
	MidPt;
	while(Right >= Left) {
		MidPt = (Right + Left)/2;
		pProbe = &m_pAlignSeqs[MidPt];
		if(pProbe->RefChromID > RelChromID)
			{
			Right = MidPt - 1;
			pProbe = NULL;
			continue;
			}
		if(pProbe->RefChromID < RelChromID)
			{
			Left = MidPt + 1;
			pProbe = NULL;
			continue;
			}
		// matching on chromosomes, now need to locate on the offset
		if(pProbe->RefOfs > RelOfs)	
			{
			Right = MidPt - 1;
			pProbe = NULL;
			continue;
			}

		if(((pProbe->RefOfs + pProbe->Len) - 1) < RelOfs)
			{
			Left = MidPt + 1;
			pProbe = NULL;
			continue;
			}

		// have a BA alignment which has it's reference in same range as the AB relative
		// check if BA alignment pts back to AB alignment - e.g is reciprocal
		NumPutative++;
		Diff = RelOfs - pProbe->RefOfs;
		if(pProbe->RelChromID == RefChromID &&
			pProbe->RelOfs + Diff == RefOfs &&
			pProbe->Strand == Strand)
			NumReciprocal++;
		else
			{
			// this base is not reciprocal - why not?
			if(pProbe->RelChromID != RefChromID ||  // note if wrong on either strand or chromosome 
				pProbe->Strand != Strand)
				{
				bHitDiffChromStrand = true;
				if(pProbe->RelChromID != RefChromID)// note if at least one hit is on a different chromosome
					bHitDiffChrom = true;
				if(pProbe->Strand != Strand)		// note if at least one hit is on a different strand
					bHitDiffStrand = true;
				}
			else		// else at least on same chromosome/strand!
				{
				MissDistance = abs(RefOfs - (pProbe->RelOfs + Diff)); // if only missed by N then count as a near miss
				if(pMisses != NULL)
					OutputMiss(pMisses,Strand,MissDistance,SrcFileID,BlockID,SubSeqID,InitLen,InitRefOfs,InitRelOfs,Ofs,pProbe);

				if(MissDistance <= 1)
					pStat->NumPutativeMiss[0]++;
				else
					{
					if(MissDistance <= 5)
						pStat->NumPutativeMiss[1]++;
					else
						{
						if(MissDistance <= 10)
							pStat->NumPutativeMiss[2]++;
						else
							{
							if(MissDistance <= 25)
								pStat->NumPutativeMiss[3]++;
							else
								{
								if(MissDistance <= 50)
									pStat->NumPutativeMiss[4]++;
								else
									{
									if(MissDistance <= 100)
										pStat->NumPutativeMiss[5]++;
									else
										{
										if(MissDistance <= 1000)
											pStat->NumPutativeMiss[6]++;
										else
											{
											if(MissDistance <= 10000)
												pStat->NumPutativeMiss[7]++;
											else
												if(MissDistance <= 100000)
													pStat->NumPutativeMiss[8]++;
											}
										}
									}
								}
							}
						}
					}
				}

			pProbe = NULL;
			}
		break;
		}
	}

m_pStats->NumReciprocal += NumReciprocal;
if(NumReciprocal) 
	pStat->CntReciprocal++;
pStat->Reciprocal += NumReciprocal;
pStat->NumPutative += NumPutative;

if(bHitDiffChrom)
	pStat->ErrChrom++;
if(bHitDiffStrand)
	pStat->ErrStrand++;
if(bHitDiffChromStrand)
	pStat->ErrChromStrand++;
return(NumReciprocal);
}

void
CAlignValidate::OutputMiss(FILE *pMisses,
						   char Strand,
						   int MissDistance,
						   int SrcFileID,
						   int BlockID,
						   int SubSeqID,
						   int InitLen,
						   int RefOfs,
						   int RelOfs,
						   int Ofs,
						   tsASalignSeq *pProbe)
{
static bool bFirst = true;
static char PrevStrand = ' ';
static int Num2Output = 1000000;
static int PrevDistance = 0;
static int PrevSrcFileID = 0;
static int PrevBlockID = 0;
static int PrevSubSeqID = 0;
static tsASalignSeq *pPrevProbe = NULL;

if(bFirst)
	{
	fprintf(pMisses,"\"Strand\",\"MissDistance\",\"hgsrcfile\",\"hgblock\",\"hgsseq\",\"hglen\",\"hgrefofs\",\"hgrelofs\",\"ofs\",\"mmsrcfile\",\"mmblock\",\"mmsseq\",\"mmlen\",\"mmrefofs\",\"mmrelofs\"");
	bFirst = false;
	}
if(Num2Output == 0)
	return;
if(PrevStrand != Strand ||
   PrevDistance != MissDistance ||
	PrevSrcFileID != SrcFileID ||
	PrevBlockID != BlockID ||
	PrevSubSeqID != SubSeqID ||
	pPrevProbe != pProbe)
	{
	fprintf(pMisses,"\n%c,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",Strand,MissDistance,SrcFileID,BlockID,SubSeqID,InitLen,RefOfs,RelOfs,Ofs,pProbe->SrcFileID,pProbe->BlockID,pProbe->SubSeqID,pProbe->Len,pProbe->RefOfs,pProbe->RelOfs);
	PrevDistance = MissDistance;
	PrevSrcFileID = SrcFileID;
	PrevBlockID = BlockID;
	PrevSubSeqID = SubSeqID;
	pPrevProbe = pProbe;
	PrevStrand = Strand;
	Num2Output-=1;
	}
}


// SortAlignRel
// Used to sort by RelChromID -> RelOfs
int 
CAlignValidate::SortAlignRel( const void *arg1, const void *arg2)
{
tsASalignSeq *pEl1 = *(tsASalignSeq **)arg1;
tsASalignSeq *pEl2 = *(tsASalignSeq **)arg2;

if(pEl1->RelChromID < pEl2->RelChromID)
	return(-1);
if(pEl1->RelChromID > pEl2->RelChromID)
	return(1);
if(pEl1->RelOfs < pEl2->RelOfs)
	return(-1);
if(pEl1->RelOfs > pEl2->RelOfs)
	return(1);
return(0);
}

// SortAlignRef
// Used to sort by RefChromID -> RefOfs
int 
CAlignValidate::SortAlignRef( const void *arg1, const void *arg2)
{
tsASalignSeq *pEl1 = (tsASalignSeq *)arg1;
tsASalignSeq *pEl2 = (tsASalignSeq *)arg2;

if(pEl1->RefChromID < pEl2->RefChromID)
	return(-1);
if(pEl1->RefChromID > pEl2->RefChromID)
	return(1);
if(pEl1->RefOfs < pEl2->RefOfs)
	return(-1);
if(pEl1->RefOfs > pEl2->RefOfs)
	return(1);
if(pEl1->RelChromID < pEl2->RelChromID)
	return(-1);
if(pEl1->RelChromID > pEl2->RelChromID)
	return(1);
if(pEl1->RelOfs < pEl2->RelOfs)
	return(-1);
if(pEl1->RelOfs > pEl2->RelOfs)
	return(1);
return(0);
}

// SortSpeciesChroms
// Used to sort species+chroms
int 
CAlignValidate::SortSpeciesChroms( const void *arg1, const void *arg2)
{
tsASchromname *pEl1 = (tsASchromname *)arg1;
tsASchromname *pEl2 = (tsASchromname *)arg2;
if(pEl1->SpeciesID < pEl2->SpeciesID)
	return(-1);
if(pEl1->SpeciesID > pEl2->SpeciesID)
	return(1);
return(stricmp(pEl1->szName,pEl2->szName));
}


