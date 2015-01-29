
#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32

#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

#include "AlignSubSeqs.h"

CAlignSubSeqs::CAlignSubSeqs(void)
{
m_pSubSeqs=NULL;		// pts to array holding instances of tsASSubSeq's
m_pSpecies=NULL;
m_pChroms=NULL;
m_pSubSeqStats=NULL;
m_pBiobed = NULL;
m_hRsltsFile = -1;
m_hUCSCMAF = -1;
m_hMuscleMAF = -1;
m_hClustalWMAF = -1;
Reset();
}

CAlignSubSeqs::~CAlignSubSeqs(void)
{
CloseMAFs();
if(m_pSubSeqs != NULL)
	delete m_pSubSeqs;

if(m_pSpecies != NULL)
	delete m_pSpecies;

if(m_pChroms != NULL)
	delete m_pChroms;

if(m_pSubSeqStats != NULL)
	delete m_pSubSeqStats;

if(m_pBiobed != NULL)
	delete m_pBiobed;
if(m_hRsltsFile != -1)
	close(m_hRsltsFile);
}

int
CAlignSubSeqs::Reset(void)
{
CloseMAFs();

if(m_pSubSeqs != NULL)
	delete m_pSubSeqs;

if(m_pSpecies != NULL)
	delete m_pSpecies;

if(m_pChroms != NULL)
	delete m_pChroms;

if(m_pSubSeqStats != NULL)
	delete m_pSubSeqStats;

if(m_pBiobed != NULL)
	delete m_pBiobed;
m_pBiobed = NULL;
if(m_hRsltsFile != -1)
	close(m_hRsltsFile);
m_hRsltsFile = -1;
m_NumSubSeqs = 0;				// number of alignment subsequence instances in m_pAlignSeqs
m_AllocdSubSeqs = 0;			// memory has been allocated to hold this many instances
m_pSubSeqs = NULL;				// pts to array holding instances of tsASSubSeq's
m_NumSpecies = 0;				// number of species
m_AllocdSpeciesInsts = 0;		// how many instances have been allocd to m_pSpecies
m_pSpecies = NULL;
m_NumChroms = 0;				// number of chromosomes
m_AllocdChromsInsts = 0;		// how many instances have been allocd to m_pChroms
m_pChroms = NULL;
m_NumSubSeqStats = 0;			// number of subsequence stats instances in m_pSubSeqStats
m_AllocdSubSeqStats = 0;		// memory has been allocated to hold this many instances
m_pSubSeqStats = NULL;			// pts to array holding instances of tsSSStats's
m_TotalNumReciprocal = 0;
m_TotalRefHits = 0;
m_TotalRefMisses = 0;
m_TotalRelMisses = 0;
return(eBSFSuccess);
}


int
CAlignSubSeqs::ProcessAlignment(bool bSecondary,	// true if secondary alignment being processed
				  char *pszRefSpecies,			// reference species
				  char *pszRelSpecies,			// relative species
				  char *pszAlignmentFile)		// alignment file
{
int Rslt;
int SeqIdx;
int BlockID;
int RefSpeciesID;
int RefChromID;
char RefStrand;
int RefAlignLen;
int RefChromOfs;
int RelSpeciesID;
int RelChromID;
char RelStrand;
int RelChromOfs;
etSeqBase *pRefSeq;
etSeqBase *pRelSeq;
etSeqBase *pRefBase;
etSeqBase *pRelBase;
etSeqBase RefBase;
etSeqBase RelBase;

int	SubSeqLen;
int	Matches;
int	NxtRefChromOfs;
int	NxtRelChromOfs;

int PrvRefChromID;
int PrvRelChromID;

char *pszChrom;
int SpeciesRefChromID;
int SpeciesRelChromID;

tsSSStats SubSeqStats;
tsSSStats TotalSeqStats;


CMAlignFile *pAlignments;

if( pszRefSpecies == NULL || pszRefSpecies[0]=='\0' ||
	pszRelSpecies == NULL || pszRelSpecies[0]=='\0' ||
	pszAlignmentFile == NULL || pszAlignmentFile[0]=='\0')
	return(eBSFerrParams);

if((pAlignments = new CMAlignFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create new instance of CMAlignFile");
	return(eBSFerrObj);
	}

if((Rslt=pAlignments->Open(pszAlignmentFile))!=eBSFSuccess)
	{
	while(pAlignments->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pAlignments->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open alignment file %s",pszAlignmentFile);
	return(Rslt);
	}

// ensure ref species is present in alignment
if((RefSpeciesID = pAlignments->LocateSpeciesID(pszRefSpecies))<1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate ref species '%s' in alignment file '%s'",pszRefSpecies,pszAlignmentFile);
	delete pAlignments;
	return(eBSFerrSpecies);
	}

// ensure ref species is the reference species
if(RefSpeciesID != 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Species %s is not the reference species in alignment file '%s'",pszRefSpecies,pszAlignmentFile);
	delete pAlignments;
	return(eBSFerrSpecies);
	}

// ensure rel species is present in alignment
if((RelSpeciesID = pAlignments->LocateSpeciesID(pszRelSpecies))<1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate rel species '%s' in alignment file '%s'",pszRelSpecies,pszAlignmentFile);
	delete pAlignments;
	return(eBSFerrSpecies);
	}

if((pRefSeq = new unsigned char [cASSSMaxBlockSize])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding reference species sequences",cASSSMaxBlockSize);
	delete pAlignments;
	return(eBSFerrMem);
	}

if((pRelSeq = new unsigned char [cASSSMaxBlockSize])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding relative species sequences",cASSSMaxBlockSize);
	delete pAlignments;
	delete pRefSeq;
	return(eBSFerrMem);
	}

// iterate over all blocks in alignment
Rslt = eBSFSuccess;
BlockID = 0;
PrvRefChromID = 0;
PrvRelChromID = 0;
TotalSeqStats.NumReciprocal = 0;
TotalSeqStats.RefHits = 0;
TotalSeqStats.RefMisses = 0;
TotalSeqStats.RelMisses = 0;

while(Rslt >= eBSFSuccess && (BlockID = pAlignments->NxtBlock(BlockID)) > 0)
	{

	RefChromID  = pAlignments->GetRelChromID(BlockID,RefSpeciesID);
	if(RefChromID < 1)	// if no chrom for ref species then species is missing...
		continue;
	RelChromID  = pAlignments->GetRelChromID(BlockID,RelSpeciesID);
	if(RelChromID < 1)	// if no chrom for rel species then species is missing...
		continue;
	RefStrand   = pAlignments->GetStrand(BlockID,RefSpeciesID);
	if(RefStrand != '+')	// reference alignment always is expected to be from '+' strand
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Reference block not on '+' strand for species '%s' in block %d, file '%s'",pszRefSpecies,BlockID,pszAlignmentFile);
		Rslt = eBSFerrInternal;
		continue;
		}

	RefChromOfs = pAlignments->GetRelChromOfs(BlockID,RefSpeciesID);
	RefAlignLen = pAlignments->GetAlignLen(BlockID,RefSpeciesID);
	RelStrand   = pAlignments->GetStrand(BlockID,RelSpeciesID);
    RelChromOfs = pAlignments->GetRelChromOfs(BlockID,RelSpeciesID);

	// get the actual block alignment sequences
	if((pRefBase = pAlignments->GetSeq(BlockID,RefSpeciesID))==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Missing sequence for species '%s' in block %d, file '%s'",pszRefSpecies,BlockID,pszAlignmentFile);
		Rslt = eBSFerrInternal;
		continue;
		}
	memcpy(pRefSeq,pRefBase,RefAlignLen);
	if((pRelBase = pAlignments->GetSeq(BlockID,RelSpeciesID))==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Missing sequence for species '%s' in block %d, file '%s'",pszRefSpecies,BlockID,pszAlignmentFile);
		Rslt = eBSFerrInternal;
		continue;
		}
	memcpy(pRelSeq,pRelBase,RefAlignLen);


	if(PrvRefChromID != RefChromID)
		{
		if(PrvRefChromID && !bSecondary)
			OutputResults();
		pszChrom = pAlignments->GetChromName(RefChromID);
		SpeciesRefChromID = AddSpeciesChrom(pszRefSpecies,pszChrom,true);
		PrvRefChromID = RefChromID;
		}

	if(PrvRelChromID != RelChromID)
		{
		pszChrom = pAlignments->GetChromName(RelChromID);
		SpeciesRelChromID = AddSpeciesChrom(pszRelSpecies,pszChrom,true);
		PrvRelChromID = RelChromID;
		}

	// process the alignment sequences into subsequences (as bounded by InDels)
	SubSeqLen = 0;
	Matches = 0;
	pRefBase = pRefSeq;
	pRelBase = pRelSeq;
	NxtRefChromOfs = RefChromOfs;
	NxtRelChromOfs = RelChromOfs;

	for(SeqIdx = 0; SeqIdx < RefAlignLen; SeqIdx++)
		{
		RefBase = *pRefBase++ & ~cRptMskFlg;
		RelBase = *pRelBase++ & ~cRptMskFlg;
		if(RefBase == eBaseInDel && RelBase == eBaseInDel) // slough InDels if present in both species
			continue;

		// subsequence if neither species is InDel
		if(RefBase != eBaseInDel && RelBase != eBaseInDel)
			{
			if(!SubSeqLen)				// is this subsequence just starting?
				{
				RefChromOfs = NxtRefChromOfs;
				RelChromOfs = NxtRelChromOfs;
				Matches = 0;
				}

			if(RefBase == RelBase)	// count number of matches in subsequence
				Matches++;
			NxtRefChromOfs++;
			NxtRelChromOfs++;
			SubSeqLen++;
			continue;
			}

		// one species has InDel
		if(RelBase == eBaseInDel)
			NxtRefChromOfs++;
		else
			NxtRelChromOfs++;

		// did InDel terminate previously started subsequence?
		if(SubSeqLen)	
			{
			if(bSecondary)
				Rslt = AddSubSeq(SpeciesRefChromID,RefStrand,RefChromOfs,SpeciesRelChromID,RelStrand,RelChromOfs,SubSeqLen);
			else
				{
				SubSeqStats.Matches = Matches;
				SubSeqStats.SpeciesRefChromID = SpeciesRefChromID;
				SubSeqStats.RefChromOfs = RefChromOfs;
				SubSeqStats.RefStrand = RefStrand;
				SubSeqStats.SpeciesRelChromID = SpeciesRelChromID;
                SubSeqStats.RelChromOfs = RelChromOfs;
				SubSeqStats.RelStrand = RelStrand;
				SubSeqStats.SubSeqLen = SubSeqLen;
				Rslt = CalcSeqStats(&SubSeqStats);
				m_TotalNumReciprocal += SubSeqStats.NumReciprocal;
				m_TotalRefHits += SubSeqStats.RefHits;
				m_TotalRefMisses += SubSeqStats.RefMisses;
				m_TotalRelMisses += SubSeqStats.RelMisses;
				AddSeqStats(&SubSeqStats);
				}
			SubSeqLen = 0;
			}
		}
	if(SubSeqLen)	// had a subsequence started?
		{
		if(bSecondary)
			Rslt = AddSubSeq(SpeciesRefChromID,RefStrand,RefChromOfs,SpeciesRelChromID,RelStrand,RelChromOfs,SubSeqLen);
		else
			{
			SubSeqStats.Matches = Matches;
			SubSeqStats.SpeciesRefChromID = SpeciesRefChromID;
			SubSeqStats.RefChromOfs = RefChromOfs;
			SubSeqStats.RefStrand = RefStrand;
			SubSeqStats.SpeciesRelChromID = SpeciesRelChromID;
            SubSeqStats.RelChromOfs = RelChromOfs;
			SubSeqStats.RelStrand = RelStrand;
			SubSeqStats.SubSeqLen = SubSeqLen;
			Rslt = CalcSeqStats(&SubSeqStats);
			AddSeqStats(&SubSeqStats);
			TotalSeqStats.NumReciprocal += SubSeqStats.NumReciprocal;
			TotalSeqStats.RefHits += SubSeqStats.RefHits;
			TotalSeqStats.RefMisses += SubSeqStats.RefMisses;
			TotalSeqStats.RelMisses += SubSeqStats.RelMisses;	
			}
		}
	}

if(bSecondary)
	SortSubSeqsByRef();
else
	OutputResults();

delete pAlignments;
delete pRefSeq;
delete pRelSeq;

return(Rslt);
}


// ProcessCompareAlignment
// Processes alignment blocks between MinLen and MaxLen inclusive by
// passing them to the MUSCLE and ClustalW alignment programs
int
CAlignSubSeqs::ProcessCompareAlignment(char *pszChrom,			// only process this chromosome - "*" for all
									   int Region, // blocks must be completely contained within this region
									   int MinBlockLen,	  // blocks must of at least this length
									   int MaxBlockLen,	  // and no longer than this
									   int MinNumSpecies, // blocks must have at least this number of species
									   int MaxNumSpecies, // and no more than this
									   char *pszTmpFilePrfx,	// intermediate file prefix e.g tmp 
									   char *pszTmpFilePath,	// intermediate file path e.g d:\\results\expr1 
									   char *pszExePath,	//where to load external aligners - qscore, clustalw and muscle from 
									   char *pszRefSpecies,// reference species
					   				   char *pszAlignmentFile) // alignment file
{
int Rslt;
int BlockID;
int RefSpeciesID;
int RefChromID;
char RefStrand;
int RefAlignLen;
etSeqBase *pSeq;

char *pszRefChrom;
int RefChromOfs;
int	RefChromEndOfs;
int	RegionIn;

int PrvRefChromID;
int PrvRelChromID;

int MaxNumAlignedSpecies;
int NumBlockSpecies;
int CurSpeciesID;
int NumSpeciesSeqs;
char *pszSpeciesName;
int FiltOnChromID;

CMAlignFile *pAlignments;
tsASSubSeqAligns *pSubSeqAligns;
tsASSpeciesSubSeqAlign *pSubSeqAlign;

m_bMultipleFeatBits = false;		// can't accept blocks in overlapping regions

if( pszRefSpecies == NULL || pszRefSpecies[0]=='\0' ||
	pszAlignmentFile == NULL || pszAlignmentFile[0]=='\0')
	return(eBSFerrParams);

if((pAlignments = new CMAlignFile())==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create new instance of CMAlignFile");
	return(eBSFerrObj);
	}

if((Rslt=pAlignments->Open(pszAlignmentFile))!=eBSFSuccess)
	{
	while(pAlignments->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pAlignments->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open alignment file %s",pszAlignmentFile);
	return(Rslt);
	}

MaxNumAlignedSpecies = pAlignments->GetNumSpecies();
if(MaxNumAlignedSpecies < MinNumSpecies)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Number of species - %d - in '%s' is less - %d - than minimum required",MaxNumAlignedSpecies,pszAlignmentFile,MinNumSpecies);
	delete pAlignments;
	return(Rslt);
	}

// ensure ref species is present in alignment
if((RefSpeciesID = pAlignments->LocateSpeciesID(pszRefSpecies))<1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate ref species '%s' in alignment file '%s'",pszRefSpecies,pszAlignmentFile);
	delete pAlignments;
	return(eBSFerrSpecies);
	}

// ensure ref species is the reference species
if(RefSpeciesID != 1)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Species %s is not the reference species in alignment file '%s'",pszRefSpecies,pszAlignmentFile);
	delete pAlignments;
	return(eBSFerrSpecies);
	}

// if chromosome is not all - "*" - then ensure chromosome is present
if(pszChrom[0] != '*')
	{
	FiltOnChromID = pAlignments->LocateChromID(pszRefSpecies,pszChrom);
	if(FiltOnChromID <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate species '%s' chromosome '%s' in alignment file '%s'",pszRefSpecies,pszChrom,pszAlignmentFile);
		delete pAlignments;
		return(eBSFerrSpecies);
		}
	}
else
	FiltOnChromID = 0;

// allocate memory to be able to hold all species alignment sequences
int MemReq = sizeof(tsASSubSeqAligns) + (sizeof(tsASSpeciesSubSeqAlign) * MaxNumAlignedSpecies-1);
pSubSeqAligns = (tsASSubSeqAligns *) new char [MemReq];
if(pSubSeqAligns == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory (%d bytes) for holding tsMUSCLEparams - %d species",MemReq,MaxNumAlignedSpecies);
	delete pAlignments;
	return(eBSFerrMem);
	}

pSubSeqAligns->MaxNumSpecies = MaxNumAlignedSpecies;
pSubSeqAligns->CurNumSpecies = 0;
sprintf(pSubSeqAligns->szClustalwExe,"%s\\clustalw.exe",pszExePath);
sprintf(pSubSeqAligns->szMuscleExe,"%s\\muscle.exe",pszExePath);
sprintf(pSubSeqAligns->szQScoreExe,"%s\\qscore.exe",pszExePath);
sprintf(pSubSeqAligns->szQScoreCSV,"%s\\%sqscore.csv",pszTmpFilePath,pszTmpFilePrfx);
sprintf(pSubSeqAligns->szUCSCfa,"%s\\%sUCSCBlock.fa",pszTmpFilePath,pszTmpFilePrfx);

sprintf(pSubSeqAligns->szUCSCMAFfile,"%s\\%sUSCSMAF.txt",pszTmpFilePath,pszTmpFilePrfx);
sprintf(pSubSeqAligns->szClustalWMAFfile,"%s\\%sClustalWMAF.txt",pszTmpFilePath,pszTmpFilePrfx);
sprintf(pSubSeqAligns->szMuscleMAFfile,"%s\\%sMuscleMAF.txt",pszTmpFilePath,pszTmpFilePrfx);

sprintf(pSubSeqAligns->Blocks[0].szMFA,"%s\\%sUCSC.mfa",pszTmpFilePath,pszTmpFilePrfx);
sprintf(pSubSeqAligns->Blocks[1].szMFA,"%s\\%sMuscle.mfa",pszTmpFilePath,pszTmpFilePrfx);
sprintf(pSubSeqAligns->Blocks[2].szMFA,"%s\\%sClustalw.mfa",pszTmpFilePath,pszTmpFilePrfx);

if((Rslt = OpenMAFs(pSubSeqAligns)) != eBSFSuccess)
	{
	delete pSubSeqAligns;
	delete pAlignments;
	return(Rslt);
	}

// iterate over all blocks in alignment
Rslt = eBSFSuccess;
BlockID = 0;
PrvRefChromID = 0;
PrvRelChromID = 0;

while(Rslt >= eBSFSuccess && (BlockID = pAlignments->NxtBlock(BlockID)) > 0)
	{
	NumBlockSpecies = pAlignments->GetNumSpecies(BlockID);
	if(NumBlockSpecies < MinNumSpecies || NumBlockSpecies > MaxNumSpecies)
		continue;
	RefAlignLen = pAlignments->GetAlignLen(BlockID,RefSpeciesID);
	if(RefAlignLen < MinBlockLen || RefAlignLen > MaxBlockLen)
		continue;

	// block has acceptable number of species and is within length range
	RefChromID  = pAlignments->GetRelChromID(BlockID,RefSpeciesID);
	// check if chromosome is to be processed
	if(RefChromID < 1 || (FiltOnChromID > 0 && (FiltOnChromID != RefChromID)))
		continue;

	RefStrand   = pAlignments->GetStrand(BlockID,RefSpeciesID);
	if(RefStrand != '+')	// reference alignment always is expected to be from '+' strand
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Reference block not on '+' strand for species '%s' in block %d, file '%s'",pszRefSpecies,BlockID,pszAlignmentFile);
		Rslt = eBSFerrInternal;
		continue;
		}

	// check if block completely contained within region
	pszRefChrom = pAlignments->GetChromName(RefChromID);
	RefChromOfs = pAlignments->GetRelChromOfs(BlockID,RefSpeciesID);
	RefChromEndOfs = pAlignments->GetRelChromEndOfs(BlockID,RefSpeciesID);
	RegionIn = CharacteriseRegion(pszRefChrom,RefChromOfs,RefChromEndOfs);
	if(RegionIn < 0 || (Region != -1 && RegionIn != Region))
		continue;

	pSubSeqAligns->Region = RegionIn;

	// get all sequences in this block
	CurSpeciesID = 0;
	NumSpeciesSeqs = 0;
	while((CurSpeciesID = pAlignments->GetNxtBlockSpeciesID(BlockID,CurSpeciesID)) > 0)
		{
		pSubSeqAlign = &pSubSeqAligns->Species[NumSpeciesSeqs];
		pSubSeqAlign->SpeciesID = CurSpeciesID;
		strcpy(pSubSeqAlign->szSpeciesName,pAlignments->GetSpeciesName(CurSpeciesID));
		pSubSeqAlign->ChromID = pAlignments->GetRelChromID(BlockID,CurSpeciesID);
		strcpy(pSubSeqAlign->szChromName,pAlignments->GetChromName(pSubSeqAlign->ChromID));
		pSubSeqAlign->Strand   = pAlignments->GetStrand(BlockID,CurSpeciesID);
		pSubSeqAlign->ChromOfs = pAlignments->GetRelChromOfs(BlockID,CurSpeciesID);
		pSubSeqAlign->srcSize = pAlignments->GetChromLen(pSubSeqAlign->ChromID);
		if((pSeq = pAlignments->GetSeq(BlockID,CurSpeciesID))==NULL)
			{
			pszSpeciesName = pAlignments->GetSpeciesName(CurSpeciesID);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Missing sequence for species '%s' in block %d, file '%s'",pszSpeciesName,BlockID,pszAlignmentFile);
			Rslt = eBSFerrInternal;
			continue;
			}
		memcpy(pSubSeqAlign->UCSCSeq,pSeq,RefAlignLen);
		NumSpeciesSeqs++;
		}
	pSubSeqAligns->CurNumSpecies = NumSpeciesSeqs;
	pSubSeqAligns->Blocks[0].SeqInDelsLen = RefAlignLen;
	pSubSeqAligns->BlockID = BlockID;

	// now see what clustalw and muscle think of this alignment!
	Rslt = ProcessExtern(pAlignments,
				pSubSeqAligns);
	if(Rslt == eBSFSuccess)
		OutputExternAlignResults(pSubSeqAligns,pAlignments);
	else
		if(Rslt == eBSFErrBase)	// alignments may contain 'N's so don't want to stop
			Rslt = eBSFSuccess; // processing of subsequent blocks!
	}
CloseMAFs();
delete pSubSeqAligns;
delete pAlignments;
return(Rslt);
}


int 
CAlignSubSeqs::ProcessExtern(CMAlignFile *pAlignments,
			tsASSubSeqAligns *pSubSeqAligns)

{
int SpeciesIdx;
int BuffIdx;
int Idx;
int LinePsn;
int hUCSCfile;

int Rslt;

char *pszInBuff;
etSeqBase *pSeq;
etSeqBase Base;
char BaseChr;
tsASSpeciesSubSeqAlign *pAlignSubSeq;
CFasta Fasta;

BuffIdx = 0;
LinePsn = 0;
pszInBuff = pSubSeqAligns->FileBuffer;
for(SpeciesIdx = 0; SpeciesIdx < pSubSeqAligns->CurNumSpecies; SpeciesIdx++)
	{
	pAlignSubSeq = &pSubSeqAligns->Species[SpeciesIdx];
	if(SpeciesIdx)
		BuffIdx += sprintf(&pszInBuff[BuffIdx],"\n\n");
	BuffIdx += sprintf(&pszInBuff[BuffIdx],">%s\n",pAlignSubSeq->szSpeciesName);
	LinePsn = 0;
	pAlignSubSeq->SeqLen = 0;
	pSeq = pAlignSubSeq->UCSCSeq;
	for(Idx=0; Idx < pSubSeqAligns->Blocks[0].SeqInDelsLen; Idx++,pSeq++)
		{
		if((Base = *pSeq) == eBaseInDel)
			continue;
		switch(Base & ~cRptMskFlg) {
			case eBaseA:
				BaseChr = 'A';
				break;
			case eBaseC:
				BaseChr = 'C';
				break;
			case eBaseG:
				BaseChr = 'G';
				break;
			case eBaseT:
				BaseChr = 'T';
				break;
			default:		// don't process blocks in which there are unrecognised bases - most likely 'N's
				return(eBSFErrBase);
			}
		pAlignSubSeq->SeqLen++;
		pszInBuff[BuffIdx++] = BaseChr;
		if(++LinePsn > 60)
			{
			BuffIdx += sprintf(&pszInBuff[BuffIdx],"\n");
			LinePsn = 0;
			}
		}
	}
pszInBuff[BuffIdx] = '\0';

#ifdef _WIN32
if((hUCSCfile = open(pSubSeqAligns->szUCSCfa, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hUCSCfile = open(pSubSeqAligns->szUCSCfa, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CAlignSubSeqs::ProcessExtern Unable to open %s - %s",pSubSeqAligns->szUCSCfa,strerror(errno));
	delete pszInBuff;
	return(eBSFerrCreateFile);
	}
CUtility::SafeWrite(hUCSCfile,pszInBuff,BuffIdx);
close(hUCSCfile);

BuffIdx = 0;
LinePsn = 0;
for(SpeciesIdx = 0; SpeciesIdx < pSubSeqAligns->CurNumSpecies; SpeciesIdx++)
	{
	pAlignSubSeq = &pSubSeqAligns->Species[SpeciesIdx];

	if(SpeciesIdx)
		BuffIdx += sprintf(&pszInBuff[BuffIdx],"\n\n");
	BuffIdx += sprintf(&pszInBuff[BuffIdx],">%s\n",pAlignSubSeq->szSpeciesName);
	LinePsn = 0;
	pSeq = pAlignSubSeq->UCSCSeq;

	for(Idx=0; Idx < pSubSeqAligns->Blocks[0].SeqInDelsLen; Idx++,pSeq++)
		{
		switch(*pSeq & ~cRptMskFlg) {
			case eBaseA:
				BaseChr = 'A';
				break;
			case eBaseC:
				BaseChr = 'C';
				break;
			case eBaseG:
				BaseChr = 'G';
				break;
			case eBaseT:
				BaseChr = 'T';
				break;
			case eBaseInDel:
				BaseChr = '-';
				break;
			default:		// don't process blocks in which there are unrecognised bases - most likely 'N's
				return(eBSFErrBase);
			}
		pszInBuff[BuffIdx++] = BaseChr;
		if(++LinePsn > 60)
			{
			BuffIdx += sprintf(&pszInBuff[BuffIdx],"\n");
			LinePsn = 0;
			}
		}
	}
pszInBuff[BuffIdx] = '\0';

#ifdef _WIN32
if((hUCSCfile = open(pSubSeqAligns->Blocks[0].szMFA, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
if((hUCSCfile = open(pSubSeqAligns->Blocks[0].szMFA, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CAlignSubSeqs::Muscle Unable to open %s - %s",pSubSeqAligns->Blocks[0].szMFA,strerror(errno));
	delete pszInBuff;
	return(eBSFerrCreateFile);
	}
CUtility::SafeWrite(hUCSCfile,pszInBuff,BuffIdx);
close(hUCSCfile);

// generate muscle alignment
Rslt = ExternAlign(true,			// false: clustalw, true: muscle
			pAlignments,
			pSubSeqAligns);

// generate clustal alignment
int ErrState;
if(Rslt >= eBSFSuccess)
	Rslt = ExternAlign(false,// false: clustalw, true: muscle
     		pAlignments,
			pSubSeqAligns);
else
	ErrState = Rslt;		// just a debug entry!

	// output aligned blocks as UCSC mfa
if(Rslt >= eBSFSuccess)
	Rslt = GenerateMAFBlocks(pAlignments,pSubSeqAligns);
else
	ErrState = Rslt;		// just a debug entry!
// score UCSC vs Clustalw, UCSC vs Muscle, and Clustalw vs Muscle
if(Rslt >= eBSFSuccess)
	Rslt = ExternScore(pSubSeqAligns);
else
	ErrState = Rslt;		// just a debug entry!
return(Rslt);
}

int 
CAlignSubSeqs::ExternAlign(bool bMuscle,	// false: clustalw, true: muscle
			CMAlignFile *pAlignments,
			tsASSubSeqAligns *pSubSeqAligns)
{
char szExecFileCmd[2000];
int SpeciesIdx;
int BuffIdx;
int LinePsn;
char *pszOutBuff;
int RIdx;
int SysRslt;

tsASSpeciesSubSeqAlign *pAlignSubSeq;
CFasta Fasta;
int FastaSeqLen;
int Descrlen;
char szDescription[cBSFDescriptionSize+1];

if(bMuscle)
	RIdx = 1;
else
    RIdx = 2;
if(bMuscle)
	sprintf(szExecFileCmd,"%s -in %s -out %s -stable", pSubSeqAligns->szMuscleExe,pSubSeqAligns->szUCSCfa,pSubSeqAligns->Blocks[RIdx].szMFA);
else	
	sprintf(szExecFileCmd,"%s /infile=%s /output=fasta /outfile=%s /outorder=input", pSubSeqAligns->szClustalwExe,pSubSeqAligns->szUCSCfa,pSubSeqAligns->Blocks[RIdx].szMFA);
if((SysRslt = system(szExecFileCmd)) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"CAlignSubSeqs::ExternAlign Unable to system('%s'), returned %d",pSubSeqAligns->szMuscleExe,SysRslt);
		return(eBSFerrExecFile);
		}

// grab the external alignment
pszOutBuff = pSubSeqAligns->FileBuffer;
if(Fasta.Open(pSubSeqAligns->Blocks[RIdx].szMFA,true)!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CAlignSubSeqs::ExternAlign Unable to open %s",pSubSeqAligns->Blocks[RIdx].szMFA);
	return(eBSFerrCreateFile);
	}

while((FastaSeqLen = Fasta.ReadSequence(pszOutBuff,sizeof(pSubSeqAligns->FileBuffer)-1)) > eBSFSuccess)
	{
	if(FastaSeqLen == eBSFFastaDescr)		// just read a descriptor line
		{
		Descrlen = Fasta.ReadDescriptor(szDescription,cBSFDescriptionSize);
		BuffIdx = 0;
		LinePsn = 0;
		for(SpeciesIdx = 0; SpeciesIdx < pSubSeqAligns->CurNumSpecies; SpeciesIdx++)
			{
			pAlignSubSeq = &pSubSeqAligns->Species[SpeciesIdx];
			if(!stricmp(szDescription,pAlignments->GetSpeciesName(pAlignSubSeq->SpeciesID)))
				break;
			}
		if(SpeciesIdx == pSubSeqAligns->CurNumSpecies)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"CAlignSubSeqs::ExternAlign Unable to match returned species %s",szDescription);
			delete pszOutBuff;
			Fasta.Close();
			return(eBSFerrSpecies);
			}
		continue;
		}
	if(bMuscle)
		{
		pSubSeqAligns->Blocks[1].SeqInDelsLen = FastaSeqLen;
		memcpy(pAlignSubSeq->MuscleSeq,pszOutBuff,FastaSeqLen);
		}
	else
		{
		pSubSeqAligns->Blocks[2].SeqInDelsLen = FastaSeqLen;
		memcpy(pAlignSubSeq->ClustalwSeq,pszOutBuff,FastaSeqLen);
		}
	}
Fasta.Close();
return(eBSFSuccess);
}

int
CAlignSubSeqs::ExternScore(tsASSubSeqAligns *pSubSeqAligns)
{
int Rslt;
int RefBlock;
int TestBlock;
int Idx;

char szExecFileCmd[_MAX_PATH];
CCSVFile QScore;

// score external aligners vs UCSC and against themselves
for(Idx = 0; Idx < 3; Idx++)
	{
	switch(Idx) {
		case 0:					// muscle vs UCSC
			RefBlock = 1;
			TestBlock = 0;
			break;
		case 1:					// clustalw vs UCSC
			RefBlock = 2;
			TestBlock = 0;
			break;
		case 2:					// muscle vs clustalw
			RefBlock = 1;
			TestBlock = 2;
			break;
		}

	sprintf(szExecFileCmd,"%s -ref %s -test %s -out %s", 
				pSubSeqAligns->szQScoreExe,
				pSubSeqAligns->Blocks[RefBlock].szMFA,
				pSubSeqAligns->Blocks[TestBlock].szMFA,
				pSubSeqAligns->szQScoreCSV);
	if(system(szExecFileCmd) < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"CAlignSubSeqs::ExternScore Unable to system('%s')",szExecFileCmd);
		return(eBSFerrExecFile);
		}

	// grab the QScore results
	Rslt = QScore.Open(pSubSeqAligns->szQScoreCSV);
	if(Rslt < eBSFSuccess)
		{
		printf("\nUnable to score alignments!");
		Rslt = eBSFSuccess;
		continue;
		}
	if(Rslt >= eBSFSuccess)
		Rslt = QScore.NextLine();
	if(Rslt >= eBSFSuccess)
		Rslt = QScore.GetDouble(3,&pSubSeqAligns->Scores[Idx].PREFABQScore);
	if(Rslt >= eBSFSuccess)
		Rslt = QScore.GetDouble(4,&pSubSeqAligns->Scores[Idx].ModelerScore);
	if(Rslt >= eBSFSuccess)
		Rslt = QScore.GetDouble(5,&pSubSeqAligns->Scores[Idx].ClineScore);
	if(Rslt >= eBSFSuccess)	
		Rslt = QScore.GetDouble(6,&pSubSeqAligns->Scores[Idx].TCScore);
	QScore.Close();

	// remove the QScore output file
	remove(pSubSeqAligns->szQScoreCSV);
	}
return(Rslt);
}



int
CAlignSubSeqs::SortSubSeqsByRef(void)
{
// now sort m_pSubSeqs by RefChromID.RefOfs
qsort(m_pSubSeqs,m_NumSubSeqs,sizeof(tsASSubSeq),SortAlignRef);

return(eBSFSuccess);
}

// AddSubSeq
// Adds the subsequence processed from secondary alignment to array ptd at by m_pSubSeqs
// When all subsequences processed then array is sorted
int					// returned m_NumSubSeqs
CAlignSubSeqs::AddSubSeq(int SpeciesRefChromID,		// reference chromosome identifier
		  char RefStrand,		// reference strand
		  int RefChromOfs,		// reference offset
		  int SpeciesRelChromID,		// relative chromosome identifier
		  char RelStrand,		// relative strand
		  int RelChromOfs,		// relative offset
		  int SubSeqLen)		// subsequence length
{
tsASSubSeq *pAlign;
if(m_pSubSeqs == NULL)
	{
	m_pSubSeqs = (tsASSubSeq *)new tsASSubSeq [cASSAllocSubSeqIncr];
	if(m_pSubSeqs == NULL)
		return(eBSFerrMem);
	m_AllocdSubSeqs = cASSAllocSubSeqIncr;
	m_NumSubSeqs = 0;
	}
else
	{
	// grow m_pAlignSeqs when needed...
	if(m_NumSubSeqs >= m_AllocdSubSeqs)
		{
		pAlign = (tsASSubSeq *)new tsASSubSeq [m_AllocdSubSeqs + cASSAllocSubSeqIncr];
		if(pAlign == NULL)
			{
			delete m_pSubSeqs;
			m_pSubSeqs = NULL;
			m_AllocdSubSeqs = 0;
			m_NumSubSeqs = 0;
			return(eBSFerrMem);
			}
		memcpy(pAlign,m_pSubSeqs,m_NumSubSeqs * sizeof(tsASSubSeq));
		delete m_pSubSeqs;
		m_pSubSeqs = pAlign;
		m_AllocdSubSeqs += cASSAllocSubSeqIncr;
		}
	}
pAlign = &m_pSubSeqs[m_NumSubSeqs++];
pAlign->SpeciesRefChromID = SpeciesRefChromID;
pAlign->RefStrand = RefStrand;
pAlign->RefChromOfs = RefChromOfs;
pAlign->SpeciesRelChromID = SpeciesRelChromID;
pAlign->RelStrand = RelStrand;
pAlign->RelChromOfs = RelChromOfs;
pAlign->SubSeqLen = SubSeqLen;
return(m_NumSubSeqs);
}

// AddSeqStats
int					// returned m_NumSubSeqStats
CAlignSubSeqs::AddSeqStats(tsSSStats *pSeqStats)		// subsequence alignment stats
{
tsSSStats *pStats;
if(m_pSubSeqStats == NULL)
	{
	m_pSubSeqStats = (tsSSStats *)new tsSSStats [cASSAllocSeqStatsIncr];
	if(m_pSubSeqStats == NULL)
		return(eBSFerrMem);
	m_AllocdSubSeqStats = cASSAllocSeqStatsIncr;
	m_NumSubSeqStats = 0;
	}
else
	{
	// grow m_pAlignSeqs when needed...
	if(m_NumSubSeqStats >= m_AllocdSubSeqStats)
		{
		pStats = (tsSSStats *)new tsSSStats [m_AllocdSubSeqStats + cASSAllocSeqStatsIncr];
		if(pStats == NULL)
			{
			delete m_pSubSeqStats;
			m_pSubSeqStats = NULL;
			m_AllocdSubSeqStats = 0;
			m_NumSubSeqStats = 0;
			return(eBSFerrMem);
			}
		memcpy(pStats,m_pSubSeqStats,m_NumSubSeqStats * sizeof(tsSSStats));
		delete m_pSubSeqStats;
		m_pSubSeqStats = pStats;
		m_AllocdSubSeqStats += cASSAllocSeqStatsIncr;
		}
	}
pStats = &m_pSubSeqStats[m_NumSubSeqStats++];
memcpy(pStats,pSeqStats,sizeof(tsSSStats));
return(m_NumSubSeqStats);
}



// AddSpeciesChrom
// If Species.Chrom previously added then returns the existing identifer
// Otherwise adds Species.Chrom and returns a count of all Species.Chrom added
// Option to check if species/chrom already added
int
CAlignSubSeqs::AddSpeciesChrom(char *pszSpecies,char *pszChrom,bool bChkDups)
{
int SpeciesID;
int ChromID;
int Idx;
tsASSChromName *pChrom;
tsASSSpeciesName *pSpecies;
unsigned short Hash;

if(m_pSpecies == NULL)
	{
	m_pSpecies = (tsASSSpeciesName *)new tsASSSpeciesName [cASSAllocSpeciesNamesIncr];
	if(m_pSpecies == NULL)
		return(eBSFerrMem);
	m_AllocdSpeciesInsts = cASSAllocSpeciesNamesIncr;
	m_NumSpecies = 0;
	}
if(m_pChroms == NULL)
	{
	m_pChroms = (tsASSChromName *)new tsASSChromName [cASSAllocChromNamesIncr];
	if(m_pChroms == NULL)
		return(eBSFerrMem);
	m_AllocdChromsInsts = cASSAllocChromNamesIncr;
	m_NumChroms = 0;
	}

// Species already known?
if((SpeciesID = GetSpeciesID(pszSpecies)) < 1)
	{
	// new species...
	if(m_NumSpecies >= m_AllocdSpeciesInsts)
		{
		pSpecies = (tsASSSpeciesName *)new tsASSSpeciesName [m_AllocdSpeciesInsts + cASSAllocSpeciesNamesIncr];
		if(pSpecies == NULL)
			return(eBSFerrMem);
		memcpy(pSpecies,m_pSpecies,m_NumSpecies * sizeof(tsASSSpeciesName));
		delete m_pSpecies;
		m_pSpecies = pSpecies;
		m_AllocdSpeciesInsts += cASSAllocSpeciesNamesIncr;
		}

	pSpecies = &m_pSpecies[m_NumSpecies++];
	SpeciesID = pSpecies->SpeciesID = m_NumSpecies;
	pSpecies->Hash = GenNameHash(pszSpecies);
	strncpy(pSpecies->szName,pszSpecies,cASSMaxNameLen);
	pSpecies->szName[cASSMaxNameLen] = '\0';
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
				return(pChrom->ChromID);
			
		}
	}

if(m_NumChroms >= m_AllocdChromsInsts)
	{
	pChrom = (tsASSChromName *)new tsASSChromName [m_AllocdChromsInsts + cASSAllocChromNamesIncr];
	if(pChrom == NULL)
		return(eBSFerrMem);
	memcpy(pChrom,m_pChroms,m_NumChroms * sizeof(tsASSChromName));
	delete m_pChroms;
	m_pChroms = pChrom;
	m_AllocdChromsInsts += cASSAllocChromNamesIncr;
	}

pChrom = &m_pChroms[m_NumChroms++];
ChromID = pChrom->ChromID = m_NumChroms;
pChrom->SpeciesID = SpeciesID;
pChrom->Hash = Hash;
strncpy(pChrom->szName,pszChrom,cASSMaxNameLen);
pChrom->szName[cASSMaxNameLen] = '\0';

return(ChromID);
}


int
CAlignSubSeqs::GetChromID(char *pszSpecies,char *pszChrom)
{
int SpeciesID;
if((SpeciesID = GetSpeciesID(pszSpecies)) >= 1)
	return(GetChromID(SpeciesID,pszChrom));
return(0);
}

// GetChromID
// Binary search
int
CAlignSubSeqs::GetChromID(int SpeciesID,char *pszChrom)
{
tsASSChromName *pChrom;
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
CAlignSubSeqs::GetSpeciesID(char *pszSpecies)
{
int Idx;
unsigned short int Hash;
tsASSSpeciesName *pSpecies;
if(!m_NumSpecies || m_pSpecies == NULL || pszSpecies == NULL || *pszSpecies == '\0')
	return(0);
Hash = GenNameHash(pszSpecies);
pSpecies = m_pSpecies;
for(Idx = 0; Idx < m_NumSpecies; Idx++,pSpecies++)
	if(pSpecies->Hash == Hash &&
		!strncmp(pSpecies->szName,pszSpecies,cASSMaxNameLen))
		return(pSpecies->SpeciesID);
return(0);
}


char *
CAlignSubSeqs::GetChrom(int SpeciesRelChromID)
{
tsASSChromName *pChrom;
if(SpeciesRelChromID < 1 || SpeciesRelChromID > m_NumChroms)
	return(NULL);
pChrom = &m_pChroms[SpeciesRelChromID-1];
return(pChrom->szName);
}

char *
CAlignSubSeqs::GetSpecies(int SpeciesRelChromID)
{
tsASSChromName *pChrom;
tsASSSpeciesName *pSpecies;
if(SpeciesRelChromID < 1 || SpeciesRelChromID > m_NumChroms)
	return(NULL);
pChrom = &m_pChroms[SpeciesRelChromID-1];
pSpecies = &m_pSpecies[pChrom->SpeciesID-1];
return(pSpecies->szName);
}


// CalcSeqStats
// 
// Returns number of bases which have identical alignements in subsequence
int						// returns number of bases reciprocally aligned
CAlignSubSeqs::CalcSeqStats(tsSSStats *pSubSeqStats)
{
int Ofs;
static tsASSubSeq *pProbe = NULL;
int Diff;
int Left;
int Right;
int MidPt;
bool bMatchOnRef;
int RefChromOfs;			// primary alignment reference offset
int RelChromOfs;		    // primary alignment relative chromosome offset

pSubSeqStats->NumReciprocal = 0;
pSubSeqStats->RefHits = 0;
pSubSeqStats->RefMisses = 0;
pSubSeqStats->RelMisses = 0;
RefChromOfs = pSubSeqStats->RefChromOfs;
RelChromOfs = pSubSeqStats->RelChromOfs;

pProbe = NULL;
for(Ofs = 0; Ofs < pSubSeqStats->SubSeqLen; Ofs++, RefChromOfs++, RelChromOfs++)
	{
	bMatchOnRef = false;
	if(pProbe != NULL)		// optimisation: if had success last probe then chances are that pProbe is still valid!
		{
		if(pProbe->RefStrand == pSubSeqStats->RefStrand &&
			pProbe->SpeciesRefChromID == pSubSeqStats->SpeciesRefChromID &&
			pProbe->RefChromOfs <= RefChromOfs && 
			(pProbe->RefChromOfs + pProbe->SubSeqLen) > RefChromOfs)
			bMatchOnRef = true;
		}

	if(!bMatchOnRef)
		{
		// need to a binary search...
		Left = 0;
		Right = m_NumSubSeqs - 1;
		MidPt;
		while(!bMatchOnRef && Right >= Left) {
			MidPt = (Right + Left)/2;
			pProbe = &m_pSubSeqs[MidPt];
			if(pProbe->SpeciesRefChromID > pSubSeqStats->SpeciesRefChromID)
				{
				Right = MidPt - 1;
				continue;
				}
			if(pProbe->SpeciesRefChromID < pSubSeqStats->SpeciesRefChromID)
				{
				Left = MidPt + 1;
				continue;
				}
			// matching on chromsomes, now need to locate on the offset
			if(pProbe->RefChromOfs > RefChromOfs)	
				{
				Right = MidPt - 1;
				continue;
				}

			if(((pProbe->RefChromOfs + pProbe->SubSeqLen) - 1) < RefChromOfs)
				{
				Left = MidPt + 1;
				continue;
				}

			// now check the strand
			if(pProbe->RefStrand > pSubSeqStats->RefStrand)	
				{
				Right = MidPt - 1;
				continue;
				}

			if(pProbe->RefStrand < pSubSeqStats->RefStrand)
				{
				Left = MidPt + 1;
				continue;
				}

			bMatchOnRef = true;
			}
		}

	if(!bMatchOnRef)	// couldn't even match on the reference!
		{
		pSubSeqStats->RefMisses++;
		pProbe = NULL;
		continue;
		}

	pSubSeqStats->RefHits++;
	Diff = RefChromOfs - pProbe->RefChromOfs;
	if(pProbe->SpeciesRelChromID == pSubSeqStats->SpeciesRelChromID &&
			pProbe->RelChromOfs + Diff == RelChromOfs &&
			pProbe->RelStrand == pSubSeqStats->RelStrand)
		pSubSeqStats->NumReciprocal++;
	else
		{
		// this base is not reciprocal - why not?

		pSubSeqStats->RelMisses++;
		}
	}

return(pSubSeqStats->NumReciprocal);
}


// GenNameHash
// Generates a 16bit hash on specified lowercased name
// This hash can then be used to quickly eliminate probe names which can't match a target name by comparing hashes
unsigned short
CAlignSubSeqs::GenNameHash(char *pszName)
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

// SortAlignRel
// Used to sort by RelChromID -> RelChromOfs --> RefChromID -> RefChromOfs
int 
CAlignSubSeqs::SortAlignRel( const void *arg1, const void *arg2)
{
tsASSubSeq *pEl1 = *(tsASSubSeq **)arg1;
tsASSubSeq *pEl2 = *(tsASSubSeq **)arg2;
if(pEl1->SpeciesRelChromID < pEl2->SpeciesRelChromID)
	return(-1);
if(pEl1->SpeciesRelChromID > pEl2->SpeciesRelChromID)
	return(1);
if(pEl1->RelChromOfs < pEl2->RelChromOfs)
	return(-1);
if(pEl1->RelChromOfs > pEl2->RelChromOfs)
	return(1);
if(pEl1->SpeciesRefChromID < pEl2->SpeciesRefChromID)
	return(-1);
if(pEl1->SpeciesRefChromID > pEl2->SpeciesRefChromID)
	return(1);
if(pEl1->RefChromOfs < pEl2->RefChromOfs)
	return(-1);
if(pEl1->RefChromOfs > pEl2->RefChromOfs)
	return(1);
return(0);
}

// SortAlignRef
// Used to sort by AlignID --> RefChromID -> RefChromOfs --> RefStrand --> RelChromID -> RelChromOfs -> RelStrand
int 
CAlignSubSeqs::SortAlignRef( const void *arg1, const void *arg2)
{
tsASSubSeq *pEl1 = (tsASSubSeq *)arg1;
tsASSubSeq *pEl2 = (tsASSubSeq *)arg2;

if(pEl1->SpeciesRefChromID < pEl2->SpeciesRefChromID)
	return(-1);
if(pEl1->SpeciesRefChromID > pEl2->SpeciesRefChromID)
	return(1);
if(pEl1->RefChromOfs < pEl2->RefChromOfs)
	return(-1);
if(pEl1->RefChromOfs > pEl2->RefChromOfs)
	return(1);
if(pEl1->RefStrand < pEl2->RefStrand)
	return(-1);
if(pEl1->RefStrand > pEl2->RefStrand)
	return(1);
if(pEl1->SpeciesRelChromID < pEl2->SpeciesRelChromID)
	return(-1);
if(pEl1->SpeciesRelChromID > pEl2->SpeciesRelChromID)
	return(1);
if(pEl1->RelChromOfs < pEl2->RelChromOfs)
	return(-1);
if(pEl1->RelChromOfs > pEl2->RelChromOfs)
	return(1);
if(pEl1->RelStrand < pEl2->RelStrand)
	return(-1);
if(pEl1->RelStrand > pEl2->RelStrand)
	return(1);
return(0);
}

// SortSpeciesChroms
// Used to sort species+chroms
int 
CAlignSubSeqs::SortSpeciesChroms( const void *arg1, const void *arg2)
{
tsASSChromName *pEl1 = (tsASSChromName *)arg1;
tsASSChromName *pEl2 = (tsASSChromName *)arg2;
if(pEl1->SpeciesID < pEl2->SpeciesID)
	return(-1);
if(pEl1->SpeciesID > pEl2->SpeciesID)
	return(1);
return(stricmp(pEl1->szName,pEl2->szName));
}


int 
CAlignSubSeqs::OpenMAFs(tsASSubSeqAligns *pSubSeqAligns)
{
char szBuffer[512];
int BuffLen;

BuffLen = sprintf(szBuffer,"##maf version=1 scoring=maf_project.v10\n#Generated by genalignconf.exe\n");

if(pSubSeqAligns->szUCSCMAFfile[0] != '\0')
	{
#ifdef _WIN32
	if((m_hUCSCMAF = open(pSubSeqAligns->szUCSCMAFfile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hUCSCMAF = open(pSubSeqAligns->szUCSCMAFfile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"CAlignSubSeqs::OpenMAFs Unable to create or truncate UCSC MAF %s - %s",pSubSeqAligns->szUCSCMAFfile,strerror(errno));
		return(eBSFerrCreateFile);
		}
	CUtility::SafeWrite(m_hUCSCMAF,szBuffer,BuffLen);
	}
else
	m_hUCSCMAF = -1;
if(pSubSeqAligns->szMuscleMAFfile[0] != '\0')
	{
#ifdef _WIN32
	if((m_hMuscleMAF = open(pSubSeqAligns->szMuscleMAFfile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hMuscleMAF = open(pSubSeqAligns->szMuscleMAFfile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"CAlignSubSeqs::OpenMAFs Unable to create or truncate Muscle MAF %s - %s",pSubSeqAligns->szMuscleMAFfile,strerror(errno));
		if(m_hUCSCMAF != -1)
			{
			close(m_hUCSCMAF);
			m_hUCSCMAF = -1;
			}
		return(eBSFerrCreateFile);
		}
	CUtility::SafeWrite(m_hMuscleMAF,szBuffer,BuffLen);
	}
else
	m_hMuscleMAF = -1;

if(pSubSeqAligns->szClustalWMAFfile[0] != '\0')
	{
#ifdef _WIN32
	if((m_hClustalWMAF = open(pSubSeqAligns->szClustalWMAFfile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hClustalWMAF = open(pSubSeqAligns->szClustalWMAFfile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"CAlignSubSeqs::OpenMAFs Unable to create or truncate ClustalW MAF %s - %s",pSubSeqAligns->szClustalWMAFfile,strerror(errno));
		if(m_hUCSCMAF != -1)
			{
			close(m_hUCSCMAF);
			m_hUCSCMAF = -1;
			}

		if(m_hMuscleMAF != -1)
			{
			close(m_hMuscleMAF);
			m_hMuscleMAF = -1;
			}
		return(eBSFerrCreateFile);
		}
	CUtility::SafeWrite(m_hClustalWMAF,szBuffer,BuffLen);
	}
else
	m_hClustalWMAF = -1;
return(eBSFSuccess);
}

int
CAlignSubSeqs::CloseMAFs(void)
{
char szBuffer[512];
int BuffLen;
BuffLen = sprintf(szBuffer,"\n\n");
if(m_hUCSCMAF != -1)
	{
	CUtility::SafeWrite(m_hUCSCMAF,szBuffer,BuffLen);
	close(m_hUCSCMAF);
	m_hUCSCMAF = -1;
	}

if(m_hMuscleMAF != -1)
	{
	CUtility::SafeWrite(m_hMuscleMAF,szBuffer,BuffLen);
	close(m_hMuscleMAF);
	m_hMuscleMAF = -1;
	}
if(m_hClustalWMAF != -1)
	{
	CUtility::SafeWrite(m_hClustalWMAF,szBuffer,BuffLen);
	close(m_hClustalWMAF);
	m_hClustalWMAF = -1;
	}
return(eBSFSuccess);
}

// GenerateMAFBlocks
// Generates MAF formated file for specified sequence
//
//Each sequence line begins with "s" and contains the following required fields: 
//src -- Name of one of the source sequences included in the alignment. For sequences that are resident in a browser assembly, the form database.chromosome allows automatic creation of links to other assemblies. Non-browser sequences are typically referenced by the species name alone. Species names must not contain spaces: concatenate multi-word names or replace the space with an underscore. 
//start -- Start of the aligning region in the source sequence, using zero-based position coordinates. If the strand value is "-", this field defines the start relative to the reverse-complemented source sequence. 
//size -- Size of the aligning region in the source sequence. This is equal to the number of non-dash characters in the text field (see below). 
//strand -- "+" or "-". If the value is "-", the sequence aligns to the reverse-complemented source. 
//srcSize -- Size of the entire source sequence, not just the portions involved in the alignment. 
//text -- Nucleotides (or amino acids) in the alignment are represented by upper-case letters; repeats are shown in lower case. Insertions are indicated by "-". 
//
// Notes:
// MAF start coordinates are zero based (0..n) but when displayed in the UCSC browser they are one based (1..n+1)
int 
CAlignSubSeqs::GenerateMAFBlocks(CMAlignFile *pAlignments,
			tsASSubSeqAligns *pSubSeqAligns)
{
char szBuffer[512];
char szSpeciesChrom[100];
int BuffLen;
int Idx;
tsASSpeciesSubSeqAlign *pAligns;
pAligns = pSubSeqAligns->Species;

if(m_hUCSCMAF == -1 && m_hClustalWMAF == -1 && m_hMuscleMAF == -1)
	return(eBSFSuccess);
BuffLen = sprintf(szBuffer,"\n\na score=%d.0",pSubSeqAligns->BlockID);

if(m_hUCSCMAF != -1)
	CUtility::SafeWrite(m_hUCSCMAF,szBuffer,BuffLen);

if(m_hMuscleMAF != -1)
	CUtility::SafeWrite(m_hMuscleMAF,szBuffer,BuffLen);

if(m_hClustalWMAF != -1)
	CUtility::SafeWrite(m_hClustalWMAF,szBuffer,BuffLen);

for(Idx = 0; Idx < pSubSeqAligns->CurNumSpecies; Idx++,pAligns++)
	{
	sprintf(szSpeciesChrom,"%s.%s",pAligns->szSpeciesName,pAligns->szChromName);
	BuffLen = sprintf(szBuffer,"\ns %-20.20s %9.9d %9.9d %c %9.9d  \t",
		szSpeciesChrom,pAligns->ChromOfs,pAligns->SeqLen,pAligns->Strand,pAligns->srcSize);

	if(m_hUCSCMAF != -1)
		{
		CUtility::SafeWrite(m_hUCSCMAF,szBuffer,BuffLen);
		CSeqTrans::MapSeq2Ascii(pAligns->UCSCSeq,pSubSeqAligns->Blocks[0].SeqInDelsLen,pSubSeqAligns->MAFseqBuff);
		CUtility::SafeWrite(m_hUCSCMAF,pSubSeqAligns->MAFseqBuff,pSubSeqAligns->Blocks[0].SeqInDelsLen);
		}

	if(m_hClustalWMAF != -1)
		{
		CUtility::SafeWrite(m_hClustalWMAF,szBuffer,BuffLen);
		CSeqTrans::MapSeq2Ascii(pAligns->ClustalwSeq,pSubSeqAligns->Blocks[2].SeqInDelsLen,pSubSeqAligns->MAFseqBuff);
		CUtility::SafeWrite(m_hClustalWMAF,pSubSeqAligns->MAFseqBuff,pSubSeqAligns->Blocks[2].SeqInDelsLen);
		}
	if(m_hMuscleMAF != -1)
		{
		CUtility::SafeWrite(m_hMuscleMAF,szBuffer,BuffLen);
		CSeqTrans::MapSeq2Ascii(pAligns->MuscleSeq,pSubSeqAligns->Blocks[1].SeqInDelsLen,pSubSeqAligns->MAFseqBuff);
		CUtility::SafeWrite(m_hMuscleMAF,pSubSeqAligns->MAFseqBuff,pSubSeqAligns->Blocks[1].SeqInDelsLen);
		}
	}
return(eBSFSuccess);
}

int 
CAlignSubSeqs::InitResults(int Mode,char *pszRsltsFile,char *pszBiobedFile, bool bMultipleFeatBits,int RegLen)
{
int Rslt;
char szLineBuff[4096];
int NumChrs;
if(pszBiobedFile != NULL && pszBiobedFile[0] != '\0')
	{
	if((m_pBiobed = (CBEDfile *)new CBEDfile())==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"\nUnable to instantiate CBEDfile");
		return(eBSFerrObj);
		}

	if((Rslt = m_pBiobed->Open(pszBiobedFile,eBTAnyBed))!=eBSFSuccess)
		{
		while(m_pBiobed->NumErrMsgs())
			gDiagnostics.DiagOut(eDLFatal,gszProcName,m_pBiobed->GetErrMsg());
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open biobed file %s",pszBiobedFile);
		delete m_pBiobed;
		m_pBiobed = NULL;
		return(Rslt);
		}
	m_RegLen = RegLen;
	m_bMultipleFeatBits = bMultipleFeatBits;
	}
else
	{
	m_pBiobed = NULL;
	m_RegLen = 0;
	m_bMultipleFeatBits = false;
	}

#ifdef _WIN32
	if((m_hRsltsFile = open(pszRsltsFile, _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC, _S_IREAD | _S_IWRITE ))==-1)
#else
	if((m_hRsltsFile = open(pszRsltsFile, O_RDWR | O_CREAT | O_TRUNC, S_IREAD | S_IWRITE ))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to create %s - %s",pszRsltsFile,strerror(errno));
	delete m_pBiobed;
	m_pBiobed = NULL;
	return(eBSFerrCreateFile);
	}
if(!Mode)
	NumChrs = sprintf(szLineBuff,"\"Chrom\",\"Region\",\"RelStrand\",\"TotSubSeqs\",\"TotBases\",\"Matches\",\"Reciprocal\",\"RefMisses\",\"RelMisses\",\"RefHits\"\n");
else
	{
	NumChrs = sprintf(szLineBuff,"\"BlockID\",\"NumSpecies\",\"Region\",\"Chrom\",\"Ofs\",\"RefSeqLen\",");
	NumChrs += sprintf(&szLineBuff[NumChrs],"\"UCSCBlockLen\",\"MuscleBlockLen\",\"ClustalwBlockLen\",");
	NumChrs += sprintf(&szLineBuff[NumChrs],"\"MvUCline\",\"MvUModeler\",\"MvUPREFABQ\",\"MvUTC\",");
	NumChrs += sprintf(&szLineBuff[NumChrs],"\"CvUCline\",\"CvUModeler\",\"CvUPREFABQ\",\"CvUTC\",");
	NumChrs += sprintf(&szLineBuff[NumChrs],"\"MvCCline\",\"MvCModeler\",\"MvCPREFABQ\",\"MvCTC\"");
	}
NumChrs += sprintf(&szLineBuff[NumChrs],"\n");
CUtility::SafeWrite(m_hRsltsFile,szLineBuff,NumChrs);
return(eBSFSuccess);
}

int
CAlignSubSeqs::CharacteriseRegion(char *pszChrom,int ChromOfs,int ChromOfsEnd)
{
int RegionIdx;
int FeatIdx;
int BEDChromID = 0;
int FeatureBits;
int BitMsk;

if(m_pBiobed != NULL)
	{
	BEDChromID = m_pBiobed->LocateChromIDbyName(pszChrom);

	if(BEDChromID > 0)
		FeatureBits = m_pBiobed->GetFeatureBits(BEDChromID,ChromOfs,ChromOfsEnd,cRegionFeatBits,m_RegLen);
	else
		FeatureBits = 0;
	RegionIdx = 0;		// default to intergenic if no feature bits set
	if(FeatureBits)		
		{
		BitMsk = cFeatBitCDS;
		for(FeatIdx = 1; FeatIdx < 7; FeatIdx++,BitMsk <<= 1)
			{
			if(BitMsk & FeatureBits)
				{
				if(RegionIdx)			// if already have feature
					{
					RegionIdx = -1;	// although was sequence of interest, more than one feature bit so can't contribute to stats
					break;
					}
				RegionIdx = FeatIdx;
				if(m_bMultipleFeatBits)	// if multiple features allowed then don't check for any additional
					break;
				}
			}
		if(RegionIdx == -1)
			return(-1);
		switch(RegionIdx) {
			case 0: case 2: case 4: case 6:		// IG,5'UTR, Introns and 3'DS
				break;
			case 1:								// CDS
				RegionIdx = 3;
				break;
			case 3:								// 3'UTR
				RegionIdx = 5;
				break;
			case 5:								// 5'US
				RegionIdx = 1;
				break;
				}
		}
	else
		RegionIdx = 0;
	}
return(RegionIdx);
}

char *
CAlignSubSeqs::Region2Text(int RegionIdx)
{
char *pszRegion;
switch(RegionIdx) {
	case 0: pszRegion = (char *)"IG"; break;
	case 1: pszRegion = (char *)"5'US"; break;
	case 2: pszRegion = (char *)"5'UTR"; break;
	case 3: pszRegion = (char *)"CDS"; break;
	case 4: pszRegion = (char *)"INTRON"; break;
	case 5: pszRegion = (char *)"3'UTR"; break;
	case 6: pszRegion = (char *)"3'DS"; break;
	default: pszRegion = (char *)"???"; break;
	}
return(pszRegion);
}

int
CAlignSubSeqs::OutputResults(void)
{
char szLineBuff[4096];
int NumChrs;
tsSSStats *pSubSeqStat;
int RegionIdx;
int StatIdx;

int BEDChromID = 0;
int SpeciesRefChromID = 0;
char *pszChrom;


int NumRegions;
int Counts[7][4][7];
int StrandIdx;

char *pszStrand;
char *pszRegion;

if(!m_NumSubSeqStats)
	return(eBSFSuccess);

memset(Counts,0,sizeof(Counts));

NumChrs = 0;
pSubSeqStat = m_pSubSeqStats;
for(StatIdx = 0; StatIdx < m_NumSubSeqStats; StatIdx++, pSubSeqStat++)
	{
	if(pSubSeqStat->SpeciesRefChromID != SpeciesRefChromID)
		{
		SpeciesRefChromID = pSubSeqStat->SpeciesRefChromID;
		pszChrom = GetChrom(SpeciesRefChromID);
		}
	if(m_pBiobed != NULL)
		{
		RegionIdx = CharacteriseRegion(pszChrom,pSubSeqStat->RefChromOfs,pSubSeqStat->RefChromOfs + pSubSeqStat->SubSeqLen - 1);
		if(RegionIdx == -1)
			continue;	
		NumRegions = 7;
		}
	else
		{
		RegionIdx = 0;
		NumRegions = 1;
		}

	if(pSubSeqStat->RelStrand == '+')
		StrandIdx = 0;
	else
		StrandIdx = 1;

	Counts[RegionIdx][StrandIdx][0] += 1;
	Counts[RegionIdx][StrandIdx][1] += pSubSeqStat->SubSeqLen;
	Counts[RegionIdx][StrandIdx][2] += pSubSeqStat->Matches;
	Counts[RegionIdx][StrandIdx][3] += pSubSeqStat->NumReciprocal;
	Counts[RegionIdx][StrandIdx][4] += pSubSeqStat->RefMisses;	
	Counts[RegionIdx][StrandIdx][5] += pSubSeqStat->RelMisses;
	Counts[RegionIdx][StrandIdx][6] += pSubSeqStat->RefHits;
	}
NumChrs = 0;
for(RegionIdx =0; RegionIdx < NumRegions; RegionIdx++)
	{
	pszRegion = Region2Text(RegionIdx);


	for(StrandIdx = 0; StrandIdx < 2; StrandIdx++)
		{
		switch(StrandIdx) {
			case 0: pszStrand = (char *)"+"; break;
			case 1: pszStrand = (char *)"-"; break;
			}
		NumChrs += sprintf(&szLineBuff[NumChrs],"\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d\n",
			pszChrom,
			pszRegion,
			pszStrand,
			Counts[RegionIdx][StrandIdx][0],	//NumSeqs;
			Counts[RegionIdx][StrandIdx][1],	//NumBases;
			Counts[RegionIdx][StrandIdx][2],	//NumMatches;
			Counts[RegionIdx][StrandIdx][3],	//NumReciprocal,
			Counts[RegionIdx][StrandIdx][4],	//RefMisses,	
			Counts[RegionIdx][StrandIdx][5],	//RelMisses,
			Counts[RegionIdx][StrandIdx][6]);	//RefHits);
		}
						
	if(NumChrs > 2000)
		{
		CUtility::SafeWrite(m_hRsltsFile,szLineBuff,NumChrs);
		NumChrs = 0;
		}
	}
if(NumChrs)
	CUtility::SafeWrite(m_hRsltsFile,szLineBuff,NumChrs);

m_NumSubSeqStats = 0; // reset seqstats ready for next primary chromosome
return(eBSFSuccess);
}


int
CAlignSubSeqs::OutputExternAlignResults(tsASSubSeqAligns *pSubSeqAligns,
										CMAlignFile *pAlignments)
{
tsASSpeciesSubSeqAlign *pAlignSubSeq;
char szLineBuff[4096];
int Idx;
int NumChrs = 0;
pAlignSubSeq = &pSubSeqAligns->Species[0];
NumChrs = sprintf(szLineBuff,"%d,%d,%d,\"%s\",%d,%d,%d,%d,%d",
			pSubSeqAligns->BlockID,
			pSubSeqAligns->CurNumSpecies,
			pSubSeqAligns->Region,
			pAlignSubSeq->szChromName,
			pAlignSubSeq->ChromOfs,
			pAlignSubSeq->SeqLen,
			pSubSeqAligns->Blocks[0].SeqInDelsLen,
			pSubSeqAligns->Blocks[1].SeqInDelsLen,
			pSubSeqAligns->Blocks[2].SeqInDelsLen);
for(Idx=0;Idx< 3; Idx++)
	{
	NumChrs += sprintf(&szLineBuff[NumChrs],",%.3g,%.3g,%.3g,%.3g",
			pSubSeqAligns->Scores[Idx].ClineScore,
			pSubSeqAligns->Scores[Idx].ModelerScore,
			pSubSeqAligns->Scores[Idx].PREFABQScore,
			pSubSeqAligns->Scores[Idx].TCScore);
	}
NumChrs += sprintf(&szLineBuff[NumChrs],"\n");
CUtility::SafeWrite(m_hRsltsFile,szLineBuff,NumChrs);
return(eBSFSuccess);
}


int 
CAlignSubSeqs::EndResults(void)
{
if(m_hRsltsFile != -1)
	close(m_hRsltsFile);
m_hRsltsFile = -1;
if(m_pBiobed != NULL)
	delete m_pBiobed;
m_pBiobed = NULL;
return(eBSFSuccess);
}
