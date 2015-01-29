#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif
#include "./SNPsFile.h"

tsSNPEnumValue 
CSNPsFile::m_SNPMoleTypeEnums[] = {
	{"unknown",eSNPMTUnknown},
	{"genomic",eSNPNTGenomic},
	{"cDNA",eSNPMTcDNA},
	{"mito",eSNPMTmito},
	{"chloro",eSNPMTchloro},
	{"",-1}
};

tsSNPEnumValue 
CSNPsFile::m_SNPClassEnums[] = {
	{"unknown",eSNPCunknown},
	{"snp",eSNPCsnp},
	{"in-del",eSNPCindel},
	{"het",eSNPChet},
	{"microsat",eSNPCmicrosat},
	{"named",eSNPCnamed},
	{"no-variation",eSNPCnovariation},
	{"mixed",eSNPCmixed},
	{"mnp",eSNPCmnp},
	{"",-1}
};

tsSNPEnumValue 
CSNPsFile::m_SNPValidEnums[] = {
	{"unknown",eSNPVunknown},
	{"'other-pop",eSNPVotherpop},
	{"by-frequency",eSNPVbyfrequency},
	{"by-cluster",eSNPVbycluster},
	{"by-2hit-2allele",eSNPVby2hit2allele},
	{"by-hapmap",eSNPVbyhapmap},
	{"genotype",eSNPVgenotype},
	{"",-1}
};

tsSNPEnumValue 
CSNPsFile::m_SNPFuncEnums[] = {
	{"unknown",eSNPFunknown},
	{"locus-region",eSNPFlocusregion},
	{"coding",eSNPFcoding},
	{"coding-synon",eSNPFcodingsynon},
	{"coding-nonsynon",eSNPFcodingnonsynon},
	{"mrna-utr",eSNPFmrnautr},
	{"intron",eSNPFintron},
	{"splice-site",eSNPFsplicesite},
	{"reference",eSNPFreference},
	{"exception",eSNPFexception},
	{"",-1}
	};

tsSNPEnumValue 
CSNPsFile::m_SNPSourceEnums[] = {
	{"unknown",eSNPSunknown},
	{"dbSnp",eSNPSdbSnp},
	{"Affy10K",eSNPSAffy10K},
	{"Affy10Kv2",eSNPSAffy10Kv2},
	{"Affy50K_HindIII",eSNPSAffy50K_HindIII},
	{"Affy50K_XbaI",eSNPSAffy50K_XbaI},
	{"",-1}
};

tsSNPEnumValue 
CSNPsFile::m_SNPLocTypeEnums[] = {
	{"unknown",eSNPLunknown},
	{"range",eSNPLrange},
	{"exact",eSNPLexact},
	{"between",eSNPLbetween},
	{"",-1}
};


CSNPsFile::CSNPsFile(void)
{
m_pPMs = NULL;				// pts to array of tsSNPs
m_hFile = -1;
Reset();
}

CSNPsFile::~CSNPsFile(void)
{
Reset();
}


void
CSNPsFile::Reset(void)
{
if(m_hFile != -1)
	{
	close(m_hFile);
	m_hFile = -1;
	}

if(m_pPMs != NULL)
	{
	delete m_pPMs;
	m_pPMs = NULL;				// pts to array of tsSNPs
	}

m_NumPMs = 0;				// current number of polymorphisms in m_pPMs
m_AllocdPMs = 0;			// how many instances of tsSNPs m_pPMs can hold before realloc required
m_NumChroms = 0;			// no chromosomes
}

// ProcessSNPsFile
// Opens, parses (calling AddSNP()), closes specified file assumed to be in UCSC downloaded table style format
// Expected format
// <bin><sep><chrom><sep><chromStart><sep><chromEnd><sep><name><sep><score><sep><strand><sep><observed><sep><molType><sep><class><sep><valid><sep><avHet><sep><avHetSE><sep><func><sep><locType><sep><source><sep><exception><NL>
// where <sep> is the tab character
// bin 585 smallint(5) unsigned 
// chrom chr1 enum('unknown','chr','chr1','chr1_random','chr2','chr2_random','chr3','chr3_random','chr4','chr4_random','chr5','chr5_random','chr6','chr6_random','chr7','chr7_random','chr8','chr8_random','chr9','chr9_random','chr10','chr10_random','chr11','chr11_random','chr12','chr12_random','chr13','chr13_random','chr14','chr14_random','chr15','chr15_random','chr16','chr16_random','chr17','chr17_random','chr18','chr18_random','chr19','chr19_random','chr20','chr20_random','chr21','chr21_random','chr22','chr22_random','chr23','chr23_random','chr24','chr24_random','chr25','chr25_random','chr26','chr26_random','chr27','chr27_random','chr28','chr28_random','chr29','chr29_random','chr30','chr30_random','chr31','chr31_random','chr32','chr32_random','chr33','chr33_random','chr34','chr34_random','chr35','chr35_random','chr36','chr36_random','chr37','chr37_random','chr38','chr38_random','chr6_hla_hap1','chr6_hla_hap2','chrFinished','chr2L','chr2L_random','chr2R','chr2R_random','chr2h','chr2h_random','chr3L','chr3L_random','chr3R','chr3R_random','chr3h','chr3h_random','chrI','chrII','chrIII','chrIV','chrM','chrM_random','chrNA','chrNA_random','chrU','chrU_random','chrUn','chrUn_random','chrV','chrV_random','chrW','chrW_random','chrX','chrX_random','chrXh','chrXh_random','chrY','chrY_random','chrYh','chrYh_random','chrZ','chrZ_random') 
// chromStart 690 int(10) unsigned 
// chromEnd 691 int(10) unsigned 
// name rs13441248 varchar(255) 
// score 0 int(10) unsigned 
// strand + enum('?','+','-') 
// observed A/G varchar(255) 
// molType genomic enum('unknown','genomic','cDNA','mito','chloro') 
// class snp enum('unknown','snp','in-del','het','microsat','named','no-variation','mixed','mnp') 
// valid unknown set('unknown','other-pop','by-frequency','by-cluster','by-2hit-2allele','by-hapmap','genotype') 
// avHet 0 float 
// avHetSE 0 float 
// func unknown set('unknown','locus-region','coding','coding-synon','coding-nonsynon','mrna-utr','intron','splice-site','reference','exception') 
// locType exact enum('unknown','range','exact','between') 
// source dbSnp enum('unknown','dbSnp','Affy10K','Affy10Kv2','Affy50K_HindIII','Affy50K_XbaI') 
// exception   set('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24') 
teBSFrsltCodes 
CSNPsFile::ProcessSNPsFile(char *pszFileName)	// file to process
{
FILE *pSNPsStream;
int LineNum;
int ParsePsn;
char szLineBuff[cLineBuffLen];
char szChrom[51];
char szClass[51];
char szLocType[51];
char szMolType[51];
char szObserved[1024];

int chromStart;
char szName[cMaxFeatNameLen+1];
char Strand;
int Cnt;
int Idx;
int SNPsNum;
teBSFrsltCodes Rslt;
tsSNPs *pPMs;
int CurChromID;

if(pszFileName == NULL || *pszFileName == '\0')
	return(eBSFerrParams);
if((pSNPsStream = fopen(pszFileName,"r"))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open SNPs raw ascii table format file %s error: %s",pszFileName,strerror(errno));
	return(eBSFerrOpnFile);
	}
LineNum = 0;
SNPsNum = 1;
while(fgets(szLineBuff,sizeof(szLineBuff),pSNPsStream)!= NULL)
	{
	if(strlen(szLineBuff) < 5)	// simply slough lines which are too short to contain anything worth parsing
		continue;

	Cnt = sscanf(szLineBuff," %*d %50s %d %*d %50s 0 %c %n",
			szChrom,&chromStart,szName,&Strand,&ParsePsn);
	if(Cnt != 4)
		{
		if(!Cnt && !LineNum)	// assume its the header line
			continue;
		fclose(pSNPsStream);
		return(eBSFerrParse);
		}
	Cnt = sscanf(&szLineBuff[ParsePsn]," %s %s %s %*s 0 0 %*s %s ",
			szObserved,szMolType,szClass,szLocType);
			
	SNPsNum++;

	if((Rslt = (teBSFrsltCodes)AddSNP(szChrom,		    // on which chromosome
			chromStart,		// start offset (0..n)
			Strand,		// which strand ('?','+','-')
			szObserved,	// observed mutations separated by '/'
			szMolType,	// enum: Sample used to find this variant
			szClass,		// enum: Variant classification
			szLocType)) < 0)	// enum: Descr
		{
		fclose(pSNPsStream);
		return(Rslt);
		}
	LineNum++;
	}

// sort them by ChromoID->StartOfs...
gDiagnostics.DiagOut(eDLInfo,gszProcName,"All parsed, now sorting...");
qsort(m_pPMs,m_NumPMs,sizeof(tsSNPs),SortPMs);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"All sorted...");

// determinine first SNP for each chromosome
pPMs = m_pPMs;
CurChromID = 0;
for(Idx = 1; Idx <= m_NumPMs; Idx++,pPMs++)
	{
	if(CurChromID != pPMs->ChromID)
		{
		CurChromID = pPMs->ChromID;
		m_Chroms[CurChromID-1].FirstSNPID = Idx;
		}
	}
fclose(pSNPsStream);
return(eBSFSuccess);
}


// If chromosome not already in table then adds chromosome to table
int		// returned chromosome identifier (1..n) or -1 if error
CSNPsFile::AddChrom(char *pszChrom)
{
int Idx;
char ChrChrom;
tsPMChrom *pChrom;
if(pszChrom == NULL || pszChrom[0] == '\0')
	return(-1);
ChrChrom = tolower(*pszChrom);
pChrom = m_Chroms;
if(m_NumChroms > 0)
	{
	for(Idx = 0; Idx < m_NumChroms; Idx++,pChrom++)
		{
		if(ChrChrom == tolower(pChrom->szName[0]) &&
			!stricmp(pszChrom,pChrom->szName))
			return(pChrom->ChromID);
		}
	if(Idx >= cMaxPMNumChroms)
		return(-1);
	}
else
	Idx = 0;
Idx++;
pChrom->ChromID = Idx;
m_NumChroms = Idx;
strcpy(pChrom->szName,pszChrom);
return(Idx);
}


// If chromosome located in table then returns a chromosome identifier (1..n) or -1 if unable to locate
// A simple linear search, but shouldn't be too many entries in table so not a major performance hit
int		// returned chromosome identifier (1..n) or -1 if error
CSNPsFile::LocateChrom(char *pszChrom)
{
int Idx;
char ChrChrom;
tsPMChrom *pChrom;
if(pszChrom == NULL || pszChrom[0] == '\0' || !m_NumChroms)
	return(-1);
ChrChrom = tolower(*pszChrom);
pChrom = m_Chroms;
for(Idx = 0; Idx < m_NumChroms; Idx++,pChrom++)
	{
	if(ChrChrom == tolower(pChrom->szName[0]) &&
		!stricmp(pszChrom,pChrom->szName))
		return(pChrom->ChromID);
	}
return(-1);
}


// LocateChrom
// returns chromsome name or NULL if unable to locate
// A simple linear search, but shouldn't be too many entries in table so not a major performance hit
char *
CSNPsFile::LocateChrom(int ChromID)			// expected to be 1..m_NumChroms
{
tsPMChrom *pChrom;
int Idx;
if(ChromID < 1 || ChromID > m_NumChroms)
	return(NULL);
pChrom = m_Chroms;
for(Idx = 0; Idx < m_NumChroms; Idx++,pChrom++)
	if(pChrom->ChromID == ChromID)
		return(pChrom->szName);
return(NULL);
}

// AddSNPs
// Parse and add SNP
// currently only interested in a subset of SNPs - must be:
// molType == genomic
// class == snp
// locType == exact
//
int 
CSNPsFile::AddSNP(char *pszChrom,		// on which chromosome
				int iStartOfs,		// start offset (0..n)
				char Strand,		// which strand ('?','+','-')
				char *pszObserved,	// observed mutations separated by '/'
				char *pszMoleType,	// enum: Sample used to find this variant
				char *pszClass,		// enum: Variant classification
				char *pszLocType)	// enum: Describes how a segment of the reference assembly must be altered to represent the variant SNP allele
{
char *pChr;
char Chr;
tsSNPs *pTmp;
tsSNPs *pPM;
int iChromID;
int iSNPs;

// apply filtering...
if(MolType2Enum(pszMoleType)!=eSNPNTGenomic)
	return(0);
if(Class2Enum(pszClass) != eSNPCsnp)
	return(0);
if(LocType2Enum(pszLocType)!=eSNPLexact)
	return(0);

if(Strand != '+' && Strand != '-')
	return(0);

// ensure m_NumPMs has room to accept this instance
if(m_pPMs == NULL || m_NumPMs == m_AllocdPMs)
	{
	if(m_pPMs != NULL)
		{
		pTmp = (tsSNPs *)new tsSNPs [m_NumPMs + cAllocPMChunk];
		if(pTmp == NULL)
			return(-1);
		memcpy(pTmp,m_pPMs,sizeof(tsSNPs) * m_NumPMs);
		delete m_pPMs;
		m_pPMs = pTmp;
		m_AllocdPMs = m_NumPMs + cAllocPMChunk;
		}
	else
		{
		m_pPMs = new tsSNPs [cAllocPMChunk];
		if(m_pPMs == NULL)
			return(-1);
		m_AllocdPMs = cAllocPMChunk;
		m_NumPMs = 0;
		}
	}

// get chromosome identifier (adds if a new chromosome)
iChromID = AddChrom(pszChrom);
if(iChromID < 1)
	return(-1);

// parse the observed polymorphisms
iSNPs = 0;
pChr = pszObserved;
while((Chr = *pChr++) != '\0')
	{
	switch(Chr) {
		case 'a': case 'A':
			if(Strand == '+')
				iSNPs |= 0x01;
			else
				iSNPs |= 0x08;
			continue;
		case 'c': case 'C':
			if(Strand == '+')
				iSNPs |= 0x02;
			else
				iSNPs |= 0x04;
			continue;
		case 'g': case 'G':
			if(Strand == '+')
				iSNPs |= 0x04;
			else
				iSNPs |= 0x02;
			continue;
		case 't': case 'T':
			if(Strand == '+')
				iSNPs |= 0x08;
			else
				iSNPs |= 0x01;
			continue;
		default:				// supposed to be '/' but accept anything as a separator
			continue;
		}
	}

pPM = &m_pPMs[m_NumPMs++];
pPM->ChromID = iChromID;
pPM->StartOfs = iStartOfs;
pPM->Strand = Strand;
pPM->PMs = iSNPs;
return(m_NumPMs);
}

// returns next SNP identifier, identifiers are ordered by chromosome then offset
int 
CSNPsFile::GetNextSNP(int CurSNPID,  // 0 returns the 1st
					  int ChromID)	 // 0 if any chromosome
{
if(CurSNPID < 0 || CurSNPID >= m_NumPMs)
	return(-1);
if(ChromID < 0 || ChromID > m_NumChroms)
	return(eBSFerrChrom);

if(CurSNPID == 0 && ChromID != 0)	 // return 1st SNP on specified chromosome
	return(m_Chroms[ChromID-1].FirstSNPID);
CurSNPID+=1;
if(ChromID != 0 && m_pPMs[CurSNPID-1].ChromID != ChromID)
	return(0);
return(CurSNPID);
}

// returns chromosome identifier on which SNPID is located
int 
CSNPsFile::GetChromID(int SNPID)
{
if(SNPID < 1 || SNPID > m_NumPMs)
	return(-1);
return(m_pPMs[SNPID-1].ChromID);
}

// returns offset on chromosome at which SNPID is located
int 
CSNPsFile::GetOffset(int SNPID)
{
if(SNPID < 1 || SNPID > m_NumPMs)
	return(-1);
return(m_pPMs[SNPID-1].StartOfs);
}

char 
CSNPsFile::GetStrand(int SNPID)		// returns strand for this SNP
{
if(SNPID < 1 || SNPID > m_NumPMs)
	return(-1);
return(m_pPMs[SNPID-1].Strand);
}

// returns Base polymorphism (true/false)
bool 
CSNPsFile::IsPM(int SNPID,teSeqBases Base)	
{
if(SNPID < 1 || SNPID > m_NumPMs)
	return(false);
switch(Base) {
	case eBaseA:
		return(m_pPMs[SNPID-1].PMs & 0x01 ? true : false);
	case eBaseC:
		return(m_pPMs[SNPID-1].PMs & 0x02 ? true : false);
	case eBaseG:
		return(m_pPMs[SNPID-1].PMs & 0x04 ? true : false);
	case eBaseT:
		return(m_pPMs[SNPID-1].PMs & 0x08 ? true : false);
	default:
		break;	
	}
return(false);
}


teSNPMoleType
CSNPsFile::MolType2Enum(char *pszTxt)
{
tsSNPEnumValue *pEnum = m_SNPMoleTypeEnums;
while(pEnum->EnumValue != -1)
	{
	if(!stricmp(pszTxt,pEnum->pszName))
		break;
	pEnum++;
	}
return((teSNPMoleType)pEnum->EnumValue);
}

teSNPClass
CSNPsFile::Class2Enum(char *pszTxt)
{
tsSNPEnumValue *pEnum = m_SNPClassEnums;
while(pEnum->EnumValue != -1)
	{
	if(!stricmp(pszTxt,pEnum->pszName))
		break;
	pEnum++;
	}
return((teSNPClass)pEnum->EnumValue);
}

teSNPValid
CSNPsFile::Valid2Enum(char *pszTxt)
{
tsSNPEnumValue *pEnum = m_SNPValidEnums;
while(pEnum->EnumValue != -1)
	{
	if(!stricmp(pszTxt,pEnum->pszName))
		break;
	pEnum++;
	}
return((teSNPValid)pEnum->EnumValue);
}

teSNPFunc
CSNPsFile::Funct2Enum(char *pszTxt)
{
tsSNPEnumValue *pEnum = m_SNPFuncEnums;
while(pEnum->EnumValue != -1)
	{
	if(!stricmp(pszTxt,pEnum->pszName))
		break;
	pEnum++;
	}
return((teSNPFunc)pEnum->EnumValue);
}

teSNPSource
CSNPsFile::Source2Enum(char *pszTxt)
{
tsSNPEnumValue *pEnum = m_SNPSourceEnums;
while(pEnum->EnumValue != -1)
	{
	if(!stricmp(pszTxt,pEnum->pszName))
		break;
	pEnum++;
	}
return((teSNPSource)pEnum->EnumValue);
}

eSNPLocType
CSNPsFile::LocType2Enum(char *pszTxt)
{
tsSNPEnumValue *pEnum = m_SNPLocTypeEnums;
while(pEnum->EnumValue != -1)
	{
	if(!stricmp(pszTxt,pEnum->pszName))
		break;
	pEnum++;
	}
return((eSNPLocType)pEnum->EnumValue);
}

char *
CSNPsFile::MolTypeEnum2Txt(teSNPMoleType aEnum)
{
tsSNPEnumValue *pEnum = m_SNPMoleTypeEnums;
if(aEnum < 0)
	return(NULL);
while(pEnum->EnumValue != -1)
	{
	if(pEnum->EnumValue == aEnum)
		break;
	pEnum++;
	}
return((char *)pEnum->pszName);
}

char *
CSNPsFile::ClassEnum2Txt(teSNPClass aEnum)
{
tsSNPEnumValue *pEnum = m_SNPClassEnums;
if(aEnum < 0)
	return(NULL);
while(pEnum->EnumValue != -1)
	{
	if(pEnum->EnumValue == aEnum)
		break;
	pEnum++;
	}
return((char *)pEnum->pszName);
}

char *
CSNPsFile::ValidEnum2Txt(teSNPValid aEnum)
{
tsSNPEnumValue *pEnum = m_SNPValidEnums;
if(aEnum < 0)
	return(NULL);
while(pEnum->EnumValue != -1)
	{
	if(pEnum->EnumValue == aEnum)
		break;
	pEnum++;
	}
return((char *)pEnum->pszName);
}

char *
CSNPsFile::FunctEnum2Txt(teSNPFunc aEnum)
{
tsSNPEnumValue *pEnum = m_SNPFuncEnums;
if(aEnum < 0)
	return(NULL);
while(pEnum->EnumValue != -1)
	{
	if(pEnum->EnumValue == aEnum)
		break;
	pEnum++;
	}
return((char *)pEnum->pszName);
}

char *
CSNPsFile::SourceEnum2Txt(teSNPSource aEnum)
{
tsSNPEnumValue *pEnum = m_SNPSourceEnums;
if(aEnum < 0)
	return(NULL);
while(pEnum->EnumValue != -1)
	{
	if(pEnum->EnumValue == aEnum)
		break;
	pEnum++;
	}
return((char *)pEnum->pszName);
}

char *
CSNPsFile::LocTypeEnum2Txt(eSNPLocType aEnum)
{
tsSNPEnumValue *pEnum = m_SNPLocTypeEnums;
if(aEnum < 0)
	return(NULL);
while(pEnum->EnumValue != -1)
	{
	if(pEnum->EnumValue == aEnum)
		break;
	pEnum++;
	}
return((char *)pEnum->pszName);
}

// used to sort by ChromID->Start
int 
CSNPsFile::SortPMs(const void *arg1, const void *arg2)
{
tsSNPs *pEl1 = (tsSNPs *)arg1;
tsSNPs *pEl2 = (tsSNPs *)arg2;

if(pEl1->Strand != '+' || pEl2->Strand != '+')
	{
	int aaaa = 111;
	}

if(pEl1->ChromID < pEl2->ChromID)
	return(-1);
if(pEl1->ChromID > pEl2->ChromID)
	return(1);
if(pEl1->StartOfs < pEl2->StartOfs)
	return(-1);
if(pEl1->StartOfs > pEl2->StartOfs)
	return(1);
if(pEl1->Strand != '+' || pEl2->Strand != '+')
	return(0);
return(0);
}
